#!/usr/bin/env python3
"""Connected canonical continuation of the 600-cell dust slab.

Prior-art commit: 52a6d50.
Protocol commits: 393e528, c2e942d, cf00b38.
No alternate root seed or branch search is performed.
"""

from collections import Counter
import contextlib
import importlib.util
import io
from itertools import combinations
import json
import multiprocessing as mp_pool
from pathlib import Path
import sys

import mpmath as arb


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_600cell_dust_canonical_continuation.json"
PRIOR_ART_COMMIT = "52a6d50"
PROTOCOL_COMMIT = "393e528"
OUTCOME_CLARIFICATION_COMMIT = "c2e942d"
PARITY_CLARIFICATION_COMMIT = "cf00b38"
RANK_RESULT_COMMIT = "715b6ad"
RANK_INPUT = HERE / "gravity_600cell_dust_canonical_legendre_rank.json"
GLUING_INPUT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
DPS = 100
arb.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational": (arb.mpf("1e-20"), arb.mpf("1e-15")),
    "validation": (arb.mpf("3e-20"), arb.mpf("3e-15")),
}
ARITHMETIC_FLOOR = arb.mpf("1e-70")
RESIDUAL_TOLERANCE = arb.mpf("1e-50")
FULL_RESIDUAL_TOLERANCE = arb.mpf("1e-40")
TYPE_SPREAD_TOLERANCE = arb.mpf("1e-50")
ENTRY_FACTOR = arb.mpf(10)
NONZERO_FACTOR = arb.mpf(100)
MAX_NEWTON_ITERATIONS = 30
MAX_BACKTRACKING = 20
MAX_BISECTIONS = 20
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def mean(values):
    return sum(values, arb.mpf(0))/len(values)


def infinity_norm(values):
    return max(abs(value) for value in values)


def orbit_sort_key(orbit, phase):
    representative = min(orbit)
    logical = tuple(vertex % 120 for vertex in representative)
    phase_pair = tuple(sorted(phase[vertex] for vertex in logical))
    return phase_pair, tuple(sorted(orbit))


spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_canonical_continuation",
    HERE / "verify_gravity_global_regge_orbits.py",
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise

check(
    "the imported one-slab geometric/action core retains all 43 certificates",
    gro.tests == gro.passed == 43,
)


def augment_boundary_orbits(base):
    old_orbits = tuple(sorted(
        gro.orbit_partition(base["old_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    final_orbits = tuple(sorted(
        gro.orbit_partition(base["new_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    return {
        **base,
        "old_orbits": old_orbits,
        "final_orbits": final_orbits,
        "old_lookup": {
            edge: index for index, orbit in enumerate(old_orbits) for edge in orbit
        },
        "final_lookup": {
            edge: index
            for index, orbit in enumerate(final_orbits) for edge in orbit
        },
    }


models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}
for parity, model in models.items():
    check(
        f"{parity}: the reduced carrier has 30 old, 35 internal and 30 final orbits",
        len(model["old_orbits"]) == 30
        and len(model["edge_orbits"]) == 35
        and len(model["pole_orbits"]) == 5
        and len(model["final_orbits"]) == 30
        and Counter(map(len, model["old_orbits"])) == Counter({24: 30})
        and Counter(map(len, model["edge_orbits"])) == Counter({24: 35})
        and Counter(map(len, model["final_orbits"])) == Counter({24: 30}),
    )


rank_input = json.loads(RANK_INPUT.read_text())
gluing_input = json.loads(GLUING_INPUT.read_text())
check(
    "the committed rank and gluing controls authorize continuation",
    rank_input["outcomes"] == ["CANONICAL_LEGENDRE_REGULAR"]
    and rank_input["nonlinear_continuation_accepted"] is True
    and rank_input["passed"] == rank_input["tests"] == 17
    and gluing_input["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing_input["passed"] == gluing_input["tests"]
    and set(rank_input["parities"])
        == set(gluing_input["parities"])
        == {"even", "odd"},
)


ARB_I = arb.mpc(0, 1)
ARB_M_STAR = arb.mpf(10)
ARB_ZETA = (arb.pi**2*arb.sqrt(2)/50)**(arb.mpf(1)/3)
ARB_R0 = 4*ARB_M_STAR/(3*arb.pi)
ARB_L0 = ARB_ZETA*ARB_R0
ARB_L0_SQUARE = ARB_L0**2
ARB_EPSILON_3 = 2*arb.pi-5*arb.acos(arb.mpf(1)/3)
ARB_MASS = (90/arb.pi)*ARB_EPSILON_3*ARB_L0
ARB_TAU = arb.mpf("0.0102")
ARB_RHO = ARB_TAU**2
ARB_SLANT_SQUARE = ARB_L0_SQUARE-ARB_RHO
ARB_BASE_OLD = tuple(ARB_L0_SQUARE for _ in range(30))
ARB_BASE_STATE = arb.matrix([
    arb.log(ARB_SLANT_SQUARE),
    arb.log(ARB_RHO),
    arb.log(ARB_L0_SQUARE),
])


def arb_log_minus(value):
    scale = max(arb.mpf(1), abs(value))
    if abs(arb.im(value)) < arb.mpf("1e-80")*scale:
        real = arb.re(value)
        if real < 0:
            return arb.log(-real)-ARB_I*arb.pi
        return arb.log(real)
    return arb.log(value)


def arb_signed_volume_square(squared, local_vertices):
    vertices = list(local_vertices)
    dimension = len(vertices)-1
    if dimension == 0:
        return arb.mpf(1)
    base = vertices[0]
    others = vertices[1:]
    gram = arb.matrix([
        [
            (squared[base][left]+squared[base][right]-squared[left][right])/2
            for right in others
        ]
        for left in others
    ])
    return arb.det(gram)/(arb.factorial(dimension)**2)


def arb_angle_record(squared):
    gram = arb.matrix([
        [
            (squared[0][left]+squared[0][right]-squared[left][right])/2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    inverse = gram**-1
    simplex_volume_square = arb_signed_volume_square(squared, range(5))
    facet_volume_squares = {
        omitted: arb_signed_volume_square(
            squared, [vertex for vertex in range(5) if vertex != omitted]
        )
        for omitted in range(5)
    }
    leading_minors = []
    for size in range(1, 5):
        principal = arb.matrix([
            [gram[left, right] for right in range(size)]
            for left in range(size)
        ])
        leading_minors.append(arb.det(principal))
    signs = [1]
    for value in leading_minors:
        signs.append(1 if value > 0 else -1 if value < 0 else 0)
    negative_directions = None if 0 in signs else sum(
        left != right for left, right in zip(signs, signs[1:])
    )

    angles = {}
    minimum_argument = arb.inf
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge_vertices = tuple(
            vertex for vertex in range(5)
            if vertex not in (omitted_a, omitted_b)
        )
        hinge_volume_square = arb_signed_volume_square(squared, hinge_vertices)
        gram_derivative = arb.matrix(4, 4)
        opposite_edge = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                gram_derivative[left-1, right-1] = (
                    int({0, left} == opposite_edge)
                    + int({0, right} == opposite_edge)
                    - int(left != right and {left, right} == opposite_edge)
                )/2
        product = inverse*gram_derivative
        volume_derivative = simplex_volume_square*sum(
            product[index, index] for index in range(4)
        )
        denominator = (
            arb.sqrt(arb.mpc(facet_volume_squares[omitted_a]))
            * arb.sqrt(arb.mpc(facet_volume_squares[omitted_b]))
        )
        cosine = 16*volume_derivative/denominator
        sine = -arb.mpf(4)/3*(
            arb.sqrt(arb.mpc(hinge_volume_square))
            * arb.sqrt(arb.mpc(simplex_volume_square))
        )/denominator
        argument = cosine+ARB_I*sine
        minimum_argument = min(minimum_argument, abs(argument))
        angles[hinge_vertices] = -ARB_I*arb_log_minus(argument)
    return {
        "angles": angles,
        "negative_directions": negative_directions,
        "minimum_leading_minor": min(abs(value) for value in leading_minors),
        "minimum_argument": minimum_argument,
    }


def triangle_area_square(values):
    x, y, z = values
    return (2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16


def triangle_area_square_partials(values):
    x, y, z = values
    return ((y+z-x)/8, (x+z-y)/8, (x+y-z)/8)


def edge_data(model, edge, q_old, x, q_new):
    edge = tuple(sorted(edge))
    if edge in model["old_lookup"]:
        index = model["old_lookup"][edge]
        return q_old[index], index
    if edge in model["edge_to_variable"]:
        internal = model["edge_to_variable"][edge]
        value = int(model["edge_jacobian"][edge])*x[internal]
        return value, 30+internal
    if edge in model["final_lookup"]:
        index = model["final_lookup"][edge]
        return q_new[index], 65+index
    raise ValueError(f"edge absent from slab carrier: {edge}")


def simplex_squared(model, simplex, q_old, x, q_new):
    squared = [[arb.mpf(0) for _ in range(5)] for _ in range(5)]
    for left, right in combinations(range(5), 2):
        value, _ = edge_data(
            model, (simplex[left], simplex[right]), q_old, x, q_new
        )
        squared[left][right] = squared[right][left] = value
    return squared


def action_and_gradient(model, q_old, x, q_new):
    curvature = [
        arb.pi if min(orbit) in model["boundary_triangles"] else 2*arb.pi
        for orbit in model["triangle_orbits"]
    ]
    negative_counts = Counter()
    minimum_leading_minor = arb.inf
    minimum_argument = arb.inf
    for simplex_orbit in model["simplex_orbits"]:
        simplex = min(simplex_orbit)
        record = arb_angle_record(simplex_squared(model, simplex, q_old, x, q_new))
        negative_counts[record["negative_directions"]] += len(simplex_orbit)
        minimum_leading_minor = min(
            minimum_leading_minor, record["minimum_leading_minor"]
        )
        minimum_argument = min(minimum_argument, record["minimum_argument"])
        for local_hinge, angle in record["angles"].items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            curvature[model["triangle_to_orbit"][triangle]] += angle

    action_sum = arb.mpc(0)
    gradient = [arb.mpc(0) for _ in range(95)]
    for orbit_index, orbit in enumerate(model["triangle_orbits"]):
        triangle = min(orbit)
        edges = tuple(tuple(sorted(edge)) for edge in combinations(triangle, 2))
        data = tuple(edge_data(model, edge, q_old, x, q_new) for edge in edges)
        values = tuple(item[0] for item in data)
        area = arb.sqrt(arb.mpc(triangle_area_square(values)))
        triangle_curvature = curvature[orbit_index]
        multiplicity = arb.mpf(len(orbit))
        action_sum += multiplicity*area*triangle_curvature
        for (squared_value, variable), partial in zip(
            data, triangle_area_square_partials(values)
        ):
            gradient[variable] += (
                -ARB_I*(multiplicity/24)*triangle_curvature
                * partial*squared_value/(2*area)
            )
    gravitational = -ARB_I*action_sum
    dust = -(8*arb.pi*ARB_MASS/5)*sum(arb.sqrt(value) for value in x[30:35])
    for internal in range(30, 35):
        gradient[30+internal] += (
            -(4*arb.pi*ARB_MASS/5)*arb.sqrt(x[internal])/24
        )
    return gravitational+dust, tuple(gradient), {
        "negative_counts": negative_counts,
        "minimum_leading_minor": minimum_leading_minor,
        "minimum_argument": minimum_argument,
        "representative_simplices": len(model["simplex_orbits"]),
    }


def state_geometry(state):
    slant, pole, final = (arb.exp(state[index]) for index in range(3))
    return (
        tuple([slant]*30+[pole]*5),
        tuple([final]*30),
    )


def branch_pass(record, maximum_imaginary):
    return bool(
        record["negative_counts"] == Counter({1: 2400})
        and record["minimum_leading_minor"] > 0
        and record["minimum_argument"] > arb.mpf("1e-6")
        and maximum_imaginary < arb.mpf("1e-70")
    )


def evaluate_state(model, state, target):
    x, q_new = state_geometry(state)
    action, gradient, branch = action_and_gradient(
        model, ARB_BASE_OLD, x, q_new
    )
    pre = tuple(-gradient[index] for index in range(30))
    full_residual = (
        tuple(gradient[30:65])
        + tuple(pre[index]-target[index] for index in range(30))
    )
    reduced = (
        mean(gradient[30:60]),
        mean(gradient[60:65]),
        mean(tuple(pre[index]-target[index] for index in range(30))),
    )
    type_groups = (
        full_residual[:30], full_residual[30:35], full_residual[35:65]
    )
    type_spread = max(
        max(abs(value-mean(group)) for value in group) for group in type_groups
    )
    maximum_imaginary = max(
        abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient)
    )
    return {
        "action": action,
        "gradient": gradient,
        "pre": pre,
        "reduced": reduced,
        "full_residual": full_residual,
        "full_norm": infinity_norm(full_residual),
        "type_spread": type_spread,
        "branch": branch,
        "maximum_imaginary": maximum_imaginary,
        "branch_pass": branch_pass(branch, maximum_imaginary),
    }


def core_three(model, state):
    x, q_new = state_geometry(state)
    action, gradient, branch = action_and_gradient(
        model, ARB_BASE_OLD, x, q_new
    )
    values = (
        mean(gradient[30:60]),
        mean(gradient[60:65]),
        mean(tuple(-gradient[index] for index in range(30))),
    )
    maximum_imaginary = max(
        abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient)
    )
    return values, branch, maximum_imaginary


def pack_core(result):
    values, branch, maximum_imaginary = result
    return {
        "values": [
            (arb.nstr(arb.re(value), DPS+5), arb.nstr(arb.im(value), DPS+5))
            for value in values
        ],
        "negative_counts": {
            "NONE" if key is None else str(key): value
            for key, value in branch["negative_counts"].items()
        },
        "minimum_leading_minor": arb.nstr(
            branch["minimum_leading_minor"], DPS+5
        ),
        "minimum_argument": arb.nstr(branch["minimum_argument"], DPS+5),
        "maximum_imaginary": arb.nstr(maximum_imaginary, DPS+5),
    }


def unpack_core(raw):
    values = tuple(arb.mpc(arb.mpf(real), arb.mpf(imag))
                   for real, imag in raw["values"])
    branch = {
        "negative_counts": Counter({
            None if key == "NONE" else int(key): value
            for key, value in raw["negative_counts"].items()
        }),
        "minimum_leading_minor": arb.mpf(raw["minimum_leading_minor"]),
        "minimum_argument": arb.mpf(raw["minimum_argument"]),
    }
    return values, branch, arb.mpf(raw["maximum_imaginary"])


_WORKER_MODEL = None


def initialize_worker(model):
    global _WORKER_MODEL
    arb.mp.dps = DPS
    _WORKER_MODEL = model


def core_worker(task):
    state = arb.matrix([arb.mpf(value) for value in task])
    return pack_core(core_three(_WORKER_MODEL, state))


def spectral_norm(matrix):
    singular = arb.svd_r(matrix, compute_uv=False)
    return singular[0] if singular.rows else arb.mpf(0)


def maximum_entry(matrix):
    return max(abs(matrix[row, column])
               for row in range(matrix.rows) for column in range(matrix.cols))


def perturb_state(state, coordinate, delta):
    result = arb.matrix(state)
    result[coordinate] += delta
    return result


def jacobian_record(model, state, pool):
    task_keys = []
    tasks = []
    for coordinate in range(3):
        for pair in DERIVATIVE_STEPS.values():
            for h in pair:
                for sign in (1, -1):
                    key = (coordinate, arb.nstr(sign*h, 20))
                    task_keys.append(key)
                    displaced = perturb_state(state, coordinate, sign*h)
                    tasks.append(tuple(arb.nstr(displaced[index], DPS+5)
                                       for index in range(3)))
    raw_results = pool.map(core_worker, tasks, chunksize=1)
    values = {}
    all_branch_pass = True
    minimum_minor = arb.inf
    minimum_argument = arb.inf
    maximum_imaginary = arb.mpf(0)
    for key, raw in zip(task_keys, raw_results):
        vector, branch, imaginary = unpack_core(raw)
        values[key] = vector
        maximum_imaginary = max(maximum_imaginary, imaginary)
        minimum_minor = min(minimum_minor, branch["minimum_leading_minor"])
        minimum_argument = min(minimum_argument, branch["minimum_argument"])
        all_branch_pass &= branch_pass(branch, imaginary)

    matrices = {}
    for name, pair in DERIVATIVE_STEPS.items():
        matrices[name] = {}
        for label, h in zip(("primary", "shadow"), pair):
            matrix = arb.matrix(3, 3)
            positive = arb.nstr(h, 20)
            negative = arb.nstr(-h, 20)
            for column in range(3):
                plus = values[(column, positive)]
                minus = values[(column, negative)]
                for row in range(3):
                    matrix[row, column] = arb.re((plus[row]-minus[row])/(2*h))
            matrices[name][label] = matrix

    op = matrices["operational"]["primary"]
    op_shadow = matrices["operational"]["shadow"]
    val = matrices["validation"]["primary"]
    val_shadow = matrices["validation"]["shadow"]
    d_op = op-op_shadow
    d_val = val-val_shadow
    d_cross = op-val
    entry_pass = bool(all(
        abs(d_cross[row, column]) <= ENTRY_FACTOR*(
            abs(d_op[row, column])+abs(d_val[row, column])+ARITHMETIC_FLOOR
        )
        for row in range(3) for column in range(3)
    ))
    epsilon = (
        spectral_norm(d_op)+spectral_norm(d_val)+spectral_norm(d_cross)
        + ARITHMETIC_FLOOR
    )
    singular = arb.svd_r(op, compute_uv=False)
    singular_values = [singular[index] for index in range(3)]
    resolved = singular_values[-1] > NONZERO_FACTOR*epsilon
    return {
        "matrix": op,
        "validation_matrix": val,
        "singular_values": singular_values,
        "epsilon": epsilon,
        "entry_pass": entry_pass,
        "branch_pass": bool(all_branch_pass),
        "minimum_minor": minimum_minor,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
        "maximum_cross": maximum_entry(d_cross),
        "maximum_operational_proxy": maximum_entry(d_op),
        "maximum_validation_proxy": maximum_entry(d_val),
        "resolved": bool(resolved and entry_pass and all_branch_pass),
    }


def target_at(pre, forward, lam):
    return tuple(pre[index]+lam*(forward[index]-pre[index]) for index in range(30))


def state_record(state, evaluation, jacobian, lam, iterations, correction_proxy):
    x, q_new = state_geometry(state)
    rho_ratio = x[30]/ARB_RHO
    lapse_relation = (x[0]+x[30]-ARB_L0_SQUARE)/ARB_L0_SQUARE
    boundary_log = state[2]-arb.log(ARB_L0_SQUARE)
    return {
        "lambda": arb.nstr(lam, 40),
        "state_logs": [arb.nstr(state[index], 70) for index in range(3)],
        "state_log_displacement": [
            arb.nstr(state[index]-ARB_BASE_STATE[index], 50)
            for index in range(3)
        ],
        "rho_over_rho0": arb.nstr(rho_ratio, 50),
        "lapse_relation_relative": arb.nstr(lapse_relation, 50),
        "boundary_square_log_change": arb.nstr(boundary_log, 50),
        "spatial_length_log_change": arb.nstr(boundary_log/2, 50),
        "reduced_residual_norm": arb.nstr(
            infinity_norm(evaluation["reduced"]), 40
        ),
        "full_residual_norm": arb.nstr(evaluation["full_norm"], 40),
        "within_type_spread": arb.nstr(evaluation["type_spread"], 40),
        "newton_iterations": iterations,
        "final_correction_proxy": arb.nstr(correction_proxy, 40),
        "jacobian_singular_values": [
            arb.nstr(value, 50) for value in jacobian["singular_values"]
        ],
        "jacobian_epsilon": arb.nstr(jacobian["epsilon"], 40),
        "minimum_leading_minor": arb.nstr(
            min(evaluation["branch"]["minimum_leading_minor"],
                jacobian["minimum_minor"]), 40
        ),
        "minimum_angle_argument": arb.nstr(
            min(evaluation["branch"]["minimum_argument"],
                jacobian["minimum_argument"]), 40
        ),
        "maximum_imaginary": arb.nstr(
            max(evaluation["maximum_imaginary"],
                jacobian["maximum_imaginary"]), 40
        ),
    }


def solve_target(model, seed, target, lam, pool):
    state = arb.matrix(seed)
    trace = []
    correction_proxy = arb.inf
    for iteration in range(MAX_NEWTON_ITERATIONS+1):
        evaluation = evaluate_state(model, state, target)
        residual_norm = infinity_norm(evaluation["reduced"])
        if not evaluation["branch_pass"]:
            return {
                "success": False,
                "reason": "BRANCH_FAILURE",
                "state": state,
                "trace": trace,
            }
        jacobian = jacobian_record(model, state, pool)
        if not jacobian["resolved"]:
            return {
                "success": False,
                "reason": "JACOBIAN_UNRESOLVED",
                "state": state,
                "trace": trace,
            }
        if residual_norm < RESIDUAL_TOLERANCE:
            final_delta = arb.lu_solve(
                jacobian["matrix"],
                arb.matrix([-arb.re(value) for value in evaluation["reduced"]]),
            )
            correction_proxy = max(abs(final_delta[index]) for index in range(3))
            full_pass = bool(
                evaluation["full_norm"] < FULL_RESIDUAL_TOLERANCE
                and evaluation["type_spread"] < TYPE_SPREAD_TOLERANCE
            )
            return {
                "success": full_pass,
                "reason": "CONVERGED" if full_pass else "FULL_RESIDUAL_FAILED",
                "state": state,
                "evaluation": evaluation,
                "jacobian": jacobian,
                "iterations": iteration,
                "correction_proxy": correction_proxy,
                "trace": trace,
                "record": state_record(
                    state, evaluation, jacobian, lam, iteration,
                    correction_proxy,
                ),
            }
        if iteration >= MAX_NEWTON_ITERATIONS:
            return {
                "success": False,
                "reason": "ITERATION_LIMIT",
                "state": state,
                "trace": trace,
            }
        try:
            delta = arb.lu_solve(
                jacobian["matrix"],
                arb.matrix([-arb.re(value) for value in evaluation["reduced"]]),
            )
        except (ZeroDivisionError, ValueError):
            return {
                "success": False,
                "reason": "JACOBIAN_SOLVE_FAILED",
                "state": state,
                "trace": trace,
            }
        correction_proxy = max(abs(delta[index]) for index in range(3))
        accepted = False
        accepted_alpha = None
        valid_trials = 0
        for power in range(MAX_BACKTRACKING+1):
            alpha = arb.mpf(2)**(-power)
            trial = state+alpha*delta
            trial_evaluation = evaluate_state(model, trial, target)
            if not trial_evaluation["branch_pass"]:
                continue
            valid_trials += 1
            trial_norm = infinity_norm(trial_evaluation["reduced"])
            if trial_norm <= (1-alpha/4)*residual_norm:
                state = trial
                accepted = True
                accepted_alpha = alpha
                break
        trace.append({
            "iteration": iteration,
            "residual_norm": arb.nstr(residual_norm, 30),
            "correction_norm": arb.nstr(correction_proxy, 30),
            "accepted_alpha": (
                arb.nstr(accepted_alpha, 20) if accepted_alpha is not None else None
            ),
            "smallest_singular": arb.nstr(jacobian["singular_values"][-1], 30),
            "jacobian_epsilon": arb.nstr(jacobian["epsilon"], 30),
        })
        if not accepted:
            return {
                "success": False,
                "reason": (
                    "BRANCH_FAILURE" if valid_trials == 0
                    else "NO_ARMIJO_STEP"
                ),
                "state": state,
                "trace": trace,
            }
    raise AssertionError("unreachable Newton loop exit")


def predictor(state, jacobian, delta_lambda, momentum_delta_mean):
    tangent = arb.lu_solve(
        jacobian["matrix"],
        arb.matrix([0, 0, momentum_delta_mean]),
    )
    return state+delta_lambda*tangent


def target_vectors(parity):
    momentum = gluing_input["parities"][parity]["momenta"]
    pre = tuple(arb.mpf(item["real"]) for item in momentum["pre"])
    post = tuple(arb.mpf(item["real"]) for item in momentum["post"])
    old_to_final = gluing_input["parities"][parity]["geometry"][
        "old_to_final_orbit_map"
    ]
    forward = tuple(post[old_to_final[old]] for old in range(30))
    return pre, forward, old_to_final


def classify_endpoint(record):
    length_change = arb.mpf(record["spatial_length_log_change"])
    if length_change > arb.mpf("1e-12"):
        return "CANONICAL_FORWARD_ROOT_EXPANDING"
    if length_change < -arb.mpf("1e-12"):
        return "CANONICAL_FORWARD_ROOT_CONTRACTING"
    return "CANONICAL_FORWARD_ROOT_STATIC"


print("="*78)
print("CONNECTED CANONICAL CONTINUATION")
print("="*78)

fork_context = mp_pool.get_context("fork")
records = {}

for parity in ("even", "odd"):
    print(f"[{parity}] starting reproduction control", flush=True)
    model = models[parity]
    pre, forward, orbit_map = target_vectors(parity)
    momentum_delta_mean = mean(tuple(
        forward[index]-pre[index] for index in range(30)
    ))
    raw_w = arb.matrix([-ARB_RHO/ARB_SLANT_SQUARE, 1, 0])
    norm_w = arb.sqrt(30*raw_w[0]**2+5*raw_w[1]**2+30*raw_w[2]**2)
    seed = ARB_BASE_STATE+arb.mpf("1e-3")*raw_w/norm_w

    with fork_context.Pool(
        processes=8,
        initializer=initialize_worker,
        initargs=(model,),
    ) as pool:
        reproduction = solve_target(
            model, seed, target_at(pre, forward, arb.mpf(0)), arb.mpf(0), pool
        )
        recovery_error = (
            max(abs(reproduction["state"][index]-ARB_BASE_STATE[index])
                for index in range(3))
            if reproduction["success"] else arb.inf
        )
        reproduction_pass = bool(
            reproduction["success"] and recovery_error < arb.mpf("1e-10")
        )
        check(
            f"{parity}: the perturbed weakest-direction seed reproduces the published slab",
            reproduction_pass,
            f"reason={reproduction['reason']}, recovery={arb.nstr(recovery_error, 8)}",
        )

        accepted_points = []
        attempts = []
        terminal_failure = None
        if reproduction_pass:
            current_state = reproduction["state"]
            current_jacobian = reproduction["jacobian"]
            current_lambda = arb.mpf(0)
            accepted_points.append(reproduction["record"])
            for coarse_index in range(1, 65):
                coarse_lambda = arb.mpf(coarse_index)/64
                coarse_seed = predictor(
                    current_state, current_jacobian,
                    coarse_lambda-current_lambda, momentum_delta_mean,
                )
                result = solve_target(
                    model, coarse_seed,
                    target_at(pre, forward, coarse_lambda), coarse_lambda, pool,
                )
                attempts.append({
                    "kind": "coarse",
                    "lambda": arb.nstr(coarse_lambda, 30),
                    "success": result["success"],
                    "reason": result["reason"],
                    "trace": result["trace"],
                })
                if result["success"]:
                    current_state = result["state"]
                    current_jacobian = result["jacobian"]
                    current_lambda = coarse_lambda
                    accepted_points.append(result["record"])
                    print(
                        f"[{parity}] lambda={float(current_lambda):.6f} "
                        f"rho/rho0={float(arb.mpf(result['record']['rho_over_rho0'])):.3e} "
                        f"dlogL={float(arb.mpf(result['record']['spatial_length_log_change'])):.3e}",
                        flush=True,
                    )
                    continue

                upper_lambda = coarse_lambda
                last_failure = result
                for bisection in range(1, MAX_BISECTIONS+1):
                    midpoint = (current_lambda+upper_lambda)/2
                    midpoint_seed = predictor(
                        current_state, current_jacobian,
                        midpoint-current_lambda, momentum_delta_mean,
                    )
                    midpoint_result = solve_target(
                        model, midpoint_seed,
                        target_at(pre, forward, midpoint), midpoint, pool,
                    )
                    attempts.append({
                        "kind": "bisection",
                        "index": bisection,
                        "lambda": arb.nstr(midpoint, 30),
                        "success": midpoint_result["success"],
                        "reason": midpoint_result["reason"],
                        "trace": midpoint_result["trace"],
                    })
                    if midpoint_result["success"]:
                        current_state = midpoint_result["state"]
                        current_jacobian = midpoint_result["jacobian"]
                        current_lambda = midpoint
                        accepted_points.append(midpoint_result["record"])
                    else:
                        upper_lambda = midpoint
                        last_failure = midpoint_result
                terminal_failure = {
                    "failed_coarse_lambda": arb.nstr(coarse_lambda, 30),
                    "last_accepted_lambda": arb.nstr(current_lambda, 30),
                    "upper_failed_lambda": arb.nstr(upper_lambda, 30),
                    "reason": last_failure["reason"],
                }
                break

    if not reproduction_pass:
        outcome = "CANONICAL_CONTINUATION_CONTROL_FAILED"
    elif current_lambda == 1:
        outcome = classify_endpoint(accepted_points[-1])
    else:
        last = accepted_points[-1]
        collapse = bool(
            arb.mpf("0.49") <= current_lambda < arb.mpf("0.5")
            and arb.mpf(last["rho_over_rho0"]) < arb.mpf("1e-12")
            and abs(arb.mpf(last["boundary_square_log_change"])) < arb.mpf("1e-8")
            and abs(arb.mpf(last["lapse_relation_relative"])) < arb.mpf("1e-10")
            and sum(item["kind"] == "bisection" for item in attempts)
                == MAX_BISECTIONS
        )
        if collapse:
            outcome = "CONNECTED_BRANCH_APPROACHES_ZERO_LAPSE"
        elif terminal_failure and terminal_failure["reason"] == "BRANCH_FAILURE":
            outcome = "CANONICAL_CONTINUATION_BRANCH_TERMINATED"
        else:
            outcome = "CANONICAL_CONTINUATION_NUMERICALLY_OPEN"

    check(
        f"{parity}: the connected-continuation outcome follows the frozen hierarchy",
        outcome in {
            "CANONICAL_CONTINUATION_CONTROL_FAILED",
            "CANONICAL_CONTINUATION_BRANCH_TERMINATED",
            "CANONICAL_CONTINUATION_NUMERICALLY_OPEN",
            "CONNECTED_BRANCH_APPROACHES_ZERO_LAPSE",
            "CANONICAL_FORWARD_ROOT_EXPANDING",
            "CANONICAL_FORWARD_ROOT_CONTRACTING",
            "CANONICAL_FORWARD_ROOT_STATIC",
        },
        f"outcome={outcome}, accepted_points={len(accepted_points)}",
    )
    records[parity] = {
        "outcome": outcome,
        "orbit_map": orbit_map,
        "momentum": {
            "pre_mean": arb.nstr(mean(pre), 60),
            "forward_mean": arb.nstr(mean(forward), 60),
            "delta_mean": arb.nstr(momentum_delta_mean, 60),
            "pre_spread": arb.nstr(max(abs(value-mean(pre)) for value in pre), 30),
            "forward_spread": arb.nstr(
                max(abs(value-mean(forward)) for value in forward), 30
            ),
        },
        "reproduction": {
            "success": reproduction_pass,
            "reason": reproduction["reason"],
            "recovery_error": arb.nstr(recovery_error, 40),
            "trace": reproduction["trace"],
            "record": reproduction.get("record"),
        },
        "accepted_points": accepted_points,
        "attempts": attempts,
        "terminal_failure": terminal_failure,
    }


outcomes = {record["outcome"] for record in records.values()}
forward_outcomes = {
    "CANONICAL_FORWARD_ROOT_EXPANDING",
    "CANONICAL_FORWARD_ROOT_CONTRACTING",
    "CANONICAL_FORWARD_ROOT_STATIC",
}
forward_accepted = bool(outcomes and outcomes <= forward_outcomes)
parity_agreement = None
if forward_accepted:
    even_record = records["even"]["accepted_points"][-1]
    odd_record = records["odd"]["accepted_points"][-1]
    scale_difference = abs(
        arb.mpf(even_record["spatial_length_log_change"])
        - arb.mpf(odd_record["spatial_length_log_change"])
    )
    correction_band = (
        arb.mpf(even_record["final_correction_proxy"])
        + arb.mpf(odd_record["final_correction_proxy"])
    )*10
    agreement_band = max(arb.mpf("1e-12"), correction_band)
    parity_agreement = {
        "scale_difference": arb.nstr(scale_difference, 40),
        "agreement_band": arb.nstr(agreement_band, 40),
        "pass": bool(scale_difference <= agreement_band),
    }
    forward_accepted &= parity_agreement["pass"]

if forward_accepted:
    verdict = (
        "DERIVED COMPUTATIONAL LOCAL/NONLINEAR: both connected canonical "
        "branches reach the forward momentum target and agree within the "
        "frozen kinematic/error band.  This is one symmetric reduced tick, "
        "not a continuum-gravity or physical-clock result."
    )
elif outcomes == {"CONNECTED_BRANCH_APPROACHES_ZERO_LAPSE"}:
    verdict = (
        "DERIVED COMPUTATIONAL NEGATIVE ON THE FROZEN CONNECTED PATH: both "
        "branches approach the zero-lapse boundary before the forward target. "
        "No canonically selected next frame is obtained; unrelated roots are "
        "not searched."
    )
else:
    verdict = (
        "OPEN/NEGATIVE CANONICAL CONTINUATION: parity outcomes are "
        f"{sorted(outcomes)}; the first forward frame is not established."
    )

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "outcome_clarification_commit": OUTCOME_CLARIFICATION_COMMIT,
    "parity_clarification_commit": PARITY_CLARIFICATION_COMMIT,
    "rank_result_commit": RANK_RESULT_COMMIT,
    "precision_digits": DPS,
    "derivative_steps": {
        name: [arb.nstr(value, 10) for value in pair]
        for name, pair in DERIVATIVE_STEPS.items()
    },
    "coarse_grid": "lambda=k/64, k=1,...,64",
    "maximum_bisections": MAX_BISECTIONS,
    "parities": records,
    "outcomes": sorted(outcomes),
    "forward_frame_accepted": forward_accepted,
    "parity_agreement": parity_agreement,
    "labels": {
        "connected_continuation": (
            "DERIVED COMPUTATIONAL" if all(
                record["reproduction"]["success"] for record in records.values()
            ) else "CONTROL FAILURE"
        ),
        "physical_time": "NOT ESTABLISHED",
        "continuum_gravity": "NOT ESTABLISHED",
        "perturbative_shape_dynamics": "NOT TESTED",
    },
    "verdict": verdict,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2)+"\n")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOMES: {sorted(outcomes)}")
print(verdict)
raise SystemExit(0 if passed == tests else 1)

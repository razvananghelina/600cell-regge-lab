#!/usr/bin/env python3
"""Second canonical homothetic tick for the fixed-mass 600-cell dust slab.

Prior-art commit: 1865e13.
Protocol commit: fb119d1.
First-tick result commit: 46a7361.

The sole predictor repeats the first logarithmic scale increment and carries
forward its lapse.  No alternate seed, target, mass or coefficient is fitted.
"""

import ast
import contextlib
import io
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
RESPONSE_SOURCE = (
    HERE / "verify_gravity_600cell_dust_homothetic_mass_conservation.py"
)
FIRST_TICK_ARTIFACT = (
    HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
)
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_second_homothetic_tick.json"
PRIOR_ART_COMMIT = "1865e13"
PROTOCOL_COMMIT = "fb119d1"
FIRST_TICK_RESULT_COMMIT = "46a7361"
DPS = 100
arb.mp.dps = DPS
STEP_SETS = {
    "operational_primary": arb.mpf("1e-20"),
    "operational_shadow": arb.mpf("1e-15"),
    "validation_primary": arb.mpf("3e-20"),
    "validation_shadow": arb.mpf("3e-15"),
}
ARITHMETIC_FLOOR = arb.mpf("1e-60")
ENTRY_FACTOR = arb.mpf(10)
NONZERO_FACTOR = arb.mpf(100)
RESIDUAL_TOLERANCE = arb.mpf("1e-25")
MAX_ITERATIONS = 8
MAX_DAMPING = 10
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def load_response_prefix():
    tree = ast.parse(RESPONSE_SOURCE.read_text(), filename=str(RESPONSE_SOURCE))
    cut = None
    for index, node in enumerate(tree.body):
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "scale_points"
            for target in node.targets
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError("response evaluator cutoff was not found")
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    namespace = {
        "__file__": str(RESPONSE_SOURCE),
        "__name__": "second_tick_imported_response_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(RESPONSE_SOURCE), "exec"), namespace)
    return namespace


def text(value, digits=50):
    return arb.nstr(value, digits)


def mean(values):
    return sum(values, arb.mpf(0))/len(values)


def vector_norm(values):
    return arb.sqrt(sum(abs(value)**2 for value in values))


def infinity_norm(values):
    return max(abs(value) for value in values)


def spread(values):
    average = mean(values)
    return max(abs(value-average) for value in values)


def matrix_spectral_norm(matrix):
    singular = arb.svd_r(matrix, compute_uv=False)
    return singular[0]


def maximum_matrix_entry(matrix):
    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows) for column in range(matrix.cols)
    )


def real_vector(raw):
    return tuple(arb.mpf(item) for item in raw)


def maximum_imaginary(action, gradient):
    return max(abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient))


print("Second canonical homothetic 600-cell dust tick", flush=True)
response_core = load_response_prefix()
models = response_core["models"]
core = response_core["core"]
action_and_gradient = core["action_and_gradient"]
branch_pass = core["branch_pass"]

check(
    "the imported action/evaluator core retains its six controls",
    response_core["tests"] == response_core["passed"] == 6,
)
check(
    "the second-tick provenance and parity carrier are frozen",
    PRIOR_ART_COMMIT == "1865e13"
    and PROTOCOL_COMMIT == "fb119d1"
    and FIRST_TICK_RESULT_COMMIT == "46a7361"
    and set(models) == {"even", "odd"},
)

first_tick = json.loads(FIRST_TICK_ARTIFACT.read_text())
gluing = json.loads(GLUING_ARTIFACT.read_text())
upstream_ok = bool(
    first_tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and first_tick["passed"] == first_tick["tests"] == 7
    and first_tick["precision_protocol_commit"] == "3c34a59"
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing["passed"] == gluing["tests"]
    and set(first_tick["solutions"])
        == set(gluing["parities"])
        == {"even", "odd"}
)
check("the accepted first tick and gluing artifacts authorize iteration", upstream_ok)


ARB_L0_SQUARE = core["ARB_L0_SQUARE"]
ARB_RHO = core["ARB_RHO"]


def evaluate_general(parity, lower_log, increment_log, rho_relative_log):
    """Evaluate one homothetic slab with an arbitrary lower spatial scale."""
    model = models[parity]
    q_old_value = arb.exp(2*lower_log)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*(lower_log+increment_log))*ARB_L0_SQUARE
    rho_value = ARB_RHO*arb.exp(rho_relative_log)
    diagonal = arb.exp(2*lower_log+increment_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("homothetic magnitude left the positive domain")
    q_old = tuple([q_old_value]*30)
    internal = tuple([diagonal]*30+[rho_value]*5)
    q_new = tuple([q_new_value]*30)
    action, gradient, branch = action_and_gradient(
        model, q_old, internal, q_new
    )
    local = tuple(arb.re(value) for value in gradient[30:65])
    pre = tuple(-arb.re(value) for value in gradient[:30])
    post = tuple(arb.re(value) for value in gradient[65:95])
    imaginary = maximum_imaginary(action, gradient)
    return {
        "lower_log": lower_log,
        "increment_log": increment_log,
        "rho_relative_log": rho_relative_log,
        "q_old": q_old_value,
        "q_new": q_new_value,
        "rho": rho_value,
        "diagonal": diagonal,
        "action": action,
        "gradient": gradient,
        "local": local,
        "pre": pre,
        "post": post,
        "branch": branch,
        "maximum_imaginary": imaginary,
        "branch_pass": branch_pass(branch, imaginary),
    }


first_states = {}
targets = {}
target_bounds = {}
maps = {}
map_target_ok = True
for parity in models:
    solution = first_tick["solutions"][parity]
    a1, r1 = (arb.mpf(value) for value in solution["state"])
    first_states[parity] = (a1, r1)
    mapping = tuple(
        gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"]
    )
    maps[parity] = mapping
    stored_post = real_vector(solution["post_momentum"])
    targets[parity] = tuple(stored_post[index] for index in mapping)
    target_bounds[parity] = arb.mpf(solution["junction_bound"])
    map_target_ok &= bool(
        sorted(mapping) == list(range(30))
        and len(stored_post) == len(targets[parity]) == 30
        and all(arb.isfinite(value) for value in targets[parity])
        and target_bounds[parity] > 0
    )
check(
    "the vertex-derived maps give complete finite second-tick targets",
    map_target_ok,
)


reconstructed_first = {}
reconstruction_ok = True
for parity in models:
    a1, r1 = first_states[parity]
    record = evaluate_general(parity, arb.mpf(0), a1, r1)
    reconstructed_first[parity] = record
    stored = first_tick["solutions"][parity]
    stored_pre = real_vector(stored["pre_momentum"])
    stored_post = real_vector(stored["post_momentum"])
    reconstruction_ok &= bool(
        record["branch_pass"]
        and infinity_norm(record["local"]) < arb.mpf("1e-24")
        and infinity_norm(tuple(
            left-right for left, right in zip(record["pre"], stored_pre)
        )) < arb.mpf("1e-45")
        and infinity_norm(tuple(
            left-right for left, right in zip(record["post"], stored_post)
        )) < arb.mpf("1e-45")
    )
reconstruction_ok &= bool(
    infinity_norm(tuple(
        left-right for left, right in zip(
            reconstructed_first["even"]["local"],
            reconstructed_first["odd"]["local"],
        )
    )) < arb.mpf("1e-24")
    and infinity_norm(tuple(
        left-right for left, right in zip(
            reconstructed_first["even"]["post"],
            reconstructed_first["odd"]["post"],
        )
    )) < arb.mpf("1e-24")
)
check(
    "the generalized evaluator reconstructs the committed first tick",
    reconstruction_ok,
    "max local even={}, odd={}".format(
        text(infinity_norm(reconstructed_first["even"]["local"]), 8),
        text(infinity_norm(reconstructed_first["odd"]["local"]), 8),
    ),
)


def evaluate_reduced(parity, state):
    increment_log, rho_step_log = state
    a1, r1 = first_states[parity]
    record = evaluate_general(
        parity, a1, increment_log, r1+rho_step_log
    )
    momentum_residual = tuple(
        value-target for value, target in zip(record["pre"], targets[parity])
    )
    reduced = (mean(record["local"][30:]), mean(momentum_residual))
    return {
        **record,
        "state": tuple(state),
        "momentum_residual": momentum_residual,
        "reduced": reduced,
        "reduced_norm": infinity_norm(reduced),
    }


seeds = {
    parity: (first_states[parity][0], arb.mpf(0)) for parity in models
}
precontrols_ok = bool(upstream_ok and map_target_ok and reconstruction_ok)
seed_records = {}
if precontrols_ok:
    seed_records = {
        parity: evaluate_reduced(parity, seeds[parity]) for parity in models
    }
seed_branch_ok = bool(
    precontrols_ok
    and all(record["branch_pass"] for record in seed_records.values())
)
check(
    "the sole constant-increment predictor remains on the Lorentzian branch",
    seed_branch_ok,
    "even F=({},{}); odd F=({},{})".format(
        *(text(value, 8) for value in seed_records.get("even", {"reduced": (arb.nan, arb.nan)})["reduced"]),
        *(text(value, 8) for value in seed_records.get("odd", {"reduced": (arb.nan, arb.nan)})["reduced"]),
    ),
)


def perturb(state, coordinate, amount):
    result = list(state)
    result[coordinate] += amount
    return tuple(result)


def calibrated_jacobian(parity, state):
    values = {}
    all_branch_pass = True
    minimum_minor = arb.inf
    minimum_argument = arb.inf
    maximum_imaginary_value = arb.mpf(0)
    for name, step in STEP_SETS.items():
        for coordinate in range(2):
            for direction in (-1, 1):
                evaluation = evaluate_reduced(
                    parity, perturb(state, coordinate, direction*step)
                )
                values[(name, coordinate, direction)] = evaluation["reduced"]
                all_branch_pass &= evaluation["branch_pass"]
                minimum_minor = min(
                    minimum_minor,
                    evaluation["branch"]["minimum_leading_minor"],
                )
                minimum_argument = min(
                    minimum_argument,
                    evaluation["branch"]["minimum_argument"],
                )
                maximum_imaginary_value = max(
                    maximum_imaginary_value, evaluation["maximum_imaginary"]
                )

    matrices = {}
    for name, step in STEP_SETS.items():
        matrix = arb.matrix(2, 2)
        for coordinate in range(2):
            plus = values[(name, coordinate, 1)]
            minus = values[(name, coordinate, -1)]
            for row in range(2):
                matrix[row, coordinate] = (plus[row]-minus[row])/(2*step)
        matrices[name] = matrix

    operational = matrices["operational_primary"]
    d_op = operational-matrices["operational_shadow"]
    d_val = matrices["validation_primary"]-matrices["validation_shadow"]
    d_cross = operational-matrices["validation_primary"]
    entry_pass = bool(all(
        abs(d_cross[row, column]) <= ENTRY_FACTOR*(
            abs(d_op[row, column])+abs(d_val[row, column])+ARITHMETIC_FLOOR
        )
        for row in range(2) for column in range(2)
    ))
    epsilon = (
        matrix_spectral_norm(d_op)+matrix_spectral_norm(d_val)
        + matrix_spectral_norm(d_cross)+ARITHMETIC_FLOOR
    )
    singular = arb.svd_r(operational, compute_uv=False)
    singular_values = (singular[0], singular[1])
    resolved = bool(
        entry_pass and all_branch_pass
        and singular_values[1] > NONZERO_FACTOR*epsilon
    )
    return {
        "matrices": matrices,
        "operational": operational,
        "d_op": d_op,
        "d_val": d_val,
        "d_cross": d_cross,
        "entry_pass": entry_pass,
        "branch_pass": bool(all_branch_pass),
        "epsilon": epsilon,
        "singular_values": singular_values,
        "determinant": arb.det(operational),
        "condition": singular_values[0]/singular_values[1],
        "minimum_leading_minor": minimum_minor,
        "minimum_angle_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary_value,
        "resolved": resolved,
    }


solve_records = {}
all_solver_branch_ok = seed_branch_ok
all_jacobians_resolved = True

if seed_branch_ok:
    for parity in models:
        print(f"  {parity}: deterministic second-tick Newton solve", flush=True)
        state = seeds[parity]
        evaluation = seed_records[parity]
        history = []
        jacobian_attempts = []
        rejected_branch_trials = 0
        converged = False
        failure = None
        accepted_iterations = 0

        while accepted_iterations < MAX_ITERATIONS:
            if evaluation["reduced_norm"] < RESIDUAL_TOLERANCE:
                converged = True
                break
            jacobian = calibrated_jacobian(parity, state)
            jacobian_attempts.append({"state": state, "jacobian": jacobian})
            all_solver_branch_ok &= jacobian["branch_pass"]
            all_jacobians_resolved &= jacobian["resolved"]
            print(
                "    J attempt {}: singular=({}, {}) epsilon={} resolved={}".format(
                    len(jacobian_attempts),
                    text(jacobian["singular_values"][0], 8),
                    text(jacobian["singular_values"][1], 8),
                    text(jacobian["epsilon"], 8),
                    jacobian["resolved"],
                ),
                flush=True,
            )
            if not jacobian["branch_pass"]:
                failure = "JACOBIAN_BRANCH_FAILURE"
                break
            if not jacobian["resolved"]:
                failure = "JACOBIAN_ERROR_BAND_FAILURE"
                break

            rhs = arb.matrix((-evaluation["reduced"][0], -evaluation["reduced"][1]))
            correction = arb.lu_solve(jacobian["operational"], rhs)
            accepted = None
            for exponent in range(MAX_DAMPING+1):
                alpha = arb.mpf(2)**(-exponent)
                trial_state = (
                    state[0]+alpha*correction[0],
                    state[1]+alpha*correction[1],
                )
                trial = evaluate_reduced(parity, trial_state)
                if not trial["branch_pass"]:
                    rejected_branch_trials += 1
                    continue
                if trial["reduced_norm"] <= (
                    1-alpha/4
                )*evaluation["reduced_norm"]:
                    accepted = (alpha, trial_state, trial)
                    break
            if accepted is None:
                failure = "NO_ARMIJO_DAMPING"
                break

            alpha, state, evaluation = accepted
            accepted_iterations += 1
            history.append({
                "iteration": accepted_iterations,
                "state": state,
                "residual": evaluation["reduced"],
                "residual_norm": evaluation["reduced_norm"],
                "alpha": alpha,
                "correction": (correction[0], correction[1]),
                "jacobian": jacobian,
            })
            print(
                "    iter={} alpha={} u={} v={} ||F||inf={}".format(
                    accepted_iterations, text(alpha, 5), text(state[0], 14),
                    text(state[1], 14), text(evaluation["reduced_norm"], 8),
                ),
                flush=True,
            )

        if evaluation["reduced_norm"] < RESIDUAL_TOLERANCE:
            converged = True
        if not converged and failure is None:
            failure = "ITERATION_LIMIT"

        endpoint_jacobian = None
        if converged:
            endpoint_jacobian = calibrated_jacobian(parity, state)
            all_solver_branch_ok &= endpoint_jacobian["branch_pass"]
            all_jacobians_resolved &= endpoint_jacobian["resolved"]
            if not endpoint_jacobian["branch_pass"]:
                failure = "ENDPOINT_BRANCH_FAILURE"
                converged = False
            elif not endpoint_jacobian["resolved"]:
                failure = "ENDPOINT_JACOBIAN_ERROR_BAND_FAILURE"
                converged = False

        diagonal = evaluation["local"][:30]
        poles = evaluation["local"][30:]
        momentum_residual = evaluation["momentum_residual"]
        full_gate = bool(
            converged
            and infinity_norm(diagonal) < arb.mpf("1e-60")
            and infinity_norm(poles) < arb.mpf("1e-25")
            and infinity_norm(evaluation["local"]) < arb.mpf("1e-25")
            and spread(diagonal) < arb.mpf("1e-60")
            and spread(poles) < arb.mpf("1e-60")
            and vector_norm(momentum_residual) <= target_bounds[parity]
            and spread(momentum_residual) <= target_bounds[parity]
        )
        solve_records[parity] = {
            "state": state,
            "evaluation": evaluation,
            "history": history,
            "jacobian_attempts": jacobian_attempts,
            "endpoint_jacobian": endpoint_jacobian,
            "converged": converged,
            "failure": failure,
            "iterations": accepted_iterations,
            "rejected_branch_trials": rejected_branch_trials,
            "full_gate": full_gate,
            "momentum_bound": target_bounds[parity],
        }


newton_ok = bool(
    seed_branch_ok
    and len(solve_records) == 2
    and all(record["converged"] for record in solve_records.values())
)
if newton_ok:
    even = solve_records["even"]
    odd = solve_records["odd"]
    u_difference = abs(even["state"][0]-odd["state"][0])
    v_difference = abs(even["state"][1]-odd["state"][1])
    pre_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even["evaluation"]["pre"], odd["evaluation"]["pre"]
        )
    ))
    post_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even["evaluation"]["post"], odd["evaluation"]["post"]
        )
    ))
    parity_gate = bool(
        u_difference < arb.mpf("1e-25")
        and v_difference < arb.mpf("1e-25")
        and pre_difference < arb.mpf("1e-22")
        and post_difference < arb.mpf("1e-22")
    )
else:
    parity_gate = False
    u_difference = v_difference = None
    pre_difference = post_difference = None


controls_ok = bool(precontrols_ok and seed_branch_ok and all_solver_branch_ok)
full_ok = bool(
    newton_ok and all(record["full_gate"] for record in solve_records.values())
)

if not controls_ok:
    outcome = "SECOND_TICK_CONTROL_FAILED"
elif not all_jacobians_resolved:
    outcome = "SECOND_TICK_JACOBIAN_OPEN"
elif not newton_ok:
    outcome = "SECOND_TICK_NEWTON_OPEN"
elif not full_ok:
    outcome = "SECOND_TICK_FULL_SUBSTITUTION_FAILED"
elif not parity_gate:
    outcome = "SECOND_TICK_SCHEDULE_DEPENDENT"
else:
    accepted_u = solve_records["even"]["state"][0]
    if abs(accepted_u) <= arb.mpf("1e-20"):
        outcome = "SECOND_TICK_STATIONARY"
    elif accepted_u < 0:
        outcome = "SECOND_TICK_CONTINUED_CONTRACTION"
    else:
        outcome = "SECOND_TICK_TURNED_TO_EXPANSION"

check(
    "all accepted seed, Jacobian and Newton states retain the Lorentzian branch",
    all_solver_branch_ok,
)
check(
    "the frozen hierarchy assigns exactly one second-tick outcome",
    outcome in {
        "SECOND_TICK_CONTROL_FAILED",
        "SECOND_TICK_JACOBIAN_OPEN",
        "SECOND_TICK_NEWTON_OPEN",
        "SECOND_TICK_FULL_SUBSTITUTION_FAILED",
        "SECOND_TICK_SCHEDULE_DEPENDENT",
        "SECOND_TICK_STATIONARY",
        "SECOND_TICK_CONTINUED_CONTRACTION",
        "SECOND_TICK_TURNED_TO_EXPANSION",
    },
    outcome,
)


def serialize_matrix(matrix):
    return [
        [text(matrix[row, column], 50) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def serialize_jacobian(record):
    if record is None:
        return None
    return {
        "matrices": {
            name: serialize_matrix(matrix)
            for name, matrix in record["matrices"].items()
        },
        "epsilon": text(record["epsilon"], 40),
        "singular_values": [
            text(value, 50) for value in record["singular_values"]
        ],
        "determinant": text(record["determinant"], 50),
        "condition": text(record["condition"], 40),
        "entry_pass": record["entry_pass"],
        "branch_pass": record["branch_pass"],
        "resolved": record["resolved"],
        "minimum_leading_minor": text(record["minimum_leading_minor"], 40),
        "minimum_angle_argument": text(record["minimum_angle_argument"], 40),
        "maximum_imaginary": text(record["maximum_imaginary"], 40),
        "maximum_d_op": text(maximum_matrix_entry(record["d_op"]), 30),
        "maximum_d_val": text(maximum_matrix_entry(record["d_val"]), 30),
        "maximum_d_cross": text(maximum_matrix_entry(record["d_cross"]), 30),
    }


serialized_solutions = {}
for parity, result in solve_records.items():
    u, v = result["state"]
    a1, r1 = first_states[parity]
    evaluation = result["evaluation"]
    serialized_solutions[parity] = {
        "first_state": [text(a1, 60), text(r1, 60)],
        "seed": [text(value, 60) for value in seeds[parity]],
        "state": [text(u, 60), text(v, 60)],
        "scale_ratio_step": text(arb.exp(u), 60),
        "scale_ratio_from_initial": text(arb.exp(a1+u), 60),
        "rho_ratio_step": text(arb.exp(v), 60),
        "rho_ratio_from_initial": text(arb.exp(r1+v), 60),
        "tau_ratio_step": text(arb.exp(v/2), 60),
        "tau_ratio_from_initial": text(arb.exp((r1+v)/2), 60),
        "iterations": result["iterations"],
        "rejected_branch_trials": result["rejected_branch_trials"],
        "converged": result["converged"],
        "failure": result["failure"],
        "full_gate": result["full_gate"],
        "reduced_residual": [text(value, 50) for value in evaluation["reduced"]],
        "reduced_residual_norm": text(evaluation["reduced_norm"], 40),
        "diagonal_residuals": [
            text(value, 40) for value in evaluation["local"][:30]
        ],
        "pole_residuals": [
            text(value, 40) for value in evaluation["local"][30:]
        ],
        "pre_momentum": [text(value, 50) for value in evaluation["pre"]],
        "post_momentum": [text(value, 50) for value in evaluation["post"]],
        "target_momentum": [text(value, 50) for value in targets[parity]],
        "junction_residual": [
            text(value, 40) for value in evaluation["momentum_residual"]
        ],
        "junction_norm": text(vector_norm(evaluation["momentum_residual"]), 40),
        "junction_bound": text(result["momentum_bound"], 40),
        "endpoint_jacobian": serialize_jacobian(result["endpoint_jacobian"]),
        "jacobian_attempts": [
            {
                "state": [text(value, 50) for value in item["state"]],
                "jacobian": serialize_jacobian(item["jacobian"]),
            }
            for item in result["jacobian_attempts"]
        ],
        "history": [
            {
                "iteration": item["iteration"],
                "state": [text(value, 50) for value in item["state"]],
                "residual": [text(value, 40) for value in item["residual"]],
                "residual_norm": text(item["residual_norm"], 30),
                "alpha": text(item["alpha"], 20),
                "correction": [text(value, 40) for value in item["correction"]],
                "jacobian": serialize_jacobian(item["jacobian"]),
            }
            for item in result["history"]
        ],
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "first_tick_result_commit": FIRST_TICK_RESULT_COMMIT,
    "steps": {name: text(value, 20) for name, value in STEP_SETS.items()},
    "maps": {parity: list(mapping) for parity, mapping in maps.items()},
    "solutions": serialized_solutions,
    "parity_gate": {
        "u_difference": None if u_difference is None else text(u_difference, 40),
        "v_difference": None if v_difference is None else text(v_difference, 40),
        "pre_momentum_infinity_difference": (
            None if pre_difference is None else text(pre_difference, 40)
        ),
        "post_momentum_infinity_difference": (
            None if post_difference is None else text(post_difference, 40)
        ),
        "passed": parity_gate,
    },
    "classification": {
        "second_local_homothetic_move": bool(full_ok and parity_gate),
        "continued_contraction": outcome == "SECOND_TICK_CONTINUED_CONTRACTION",
        "absolute_time_unit": False,
        "full_65_variable_global_uniqueness": False,
        "anisotropic_stability": "OPEN",
        "refinement_stability": "OPEN",
        "label": "STRUCTURAL / candidate pseudo-constraint",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

for parity, result in solve_records.items():
    u, v = result["state"]
    evaluation = result["evaluation"]
    print(
        "  {}: u={} v={} L2/L1={} tau2/tau1={} ||F||={} ||junction||={}".format(
            parity, text(u, 16), text(v, 16), text(arb.exp(u), 16),
            text(arb.exp(v/2), 16), text(evaluation["reduced_norm"], 8),
            text(vector_norm(evaluation["momentum_residual"]), 8),
        ),
        flush=True,
    )
print(f"\nOUTCOME: {outcome}", flush=True)
print(f"RESULT: {passed}/{tests}", flush=True)
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Fixed 32-step canonical-target homotopy to the second 600-cell dust tick.

Prior-art commit: e760462.
Protocol commit: a2564d5.
Accepted start commit: 46a7361.
Direct-solver boundary commit: 6346ad0.
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
DIRECT_ARTIFACT = HERE / "gravity_600cell_dust_second_homothetic_tick.json"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_second_tick_homotopy.json"
PRIOR_ART_COMMIT = "e760462"
PROTOCOL_COMMIT = "a2564d5"
START_RESULT_COMMIT = "46a7361"
DIRECT_BOUNDARY_COMMIT = "6346ad0"
DPS = 100
arb.mp.dps = DPS
NODE_COUNT = 32
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
MAX_ITERATIONS = 6
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
        "__name__": "second_tick_homotopy_imported_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(RESPONSE_SOURCE), "exec"), namespace)
    return namespace


def text(value, digits=50):
    return arb.nstr(value, digits)


def mean(values):
    return sum(values, arb.mpf(0))/len(values)


def infinity_norm(values):
    return max(abs(value) for value in values)


def vector_norm(values):
    return arb.sqrt(sum(abs(value)**2 for value in values))


def spread(values):
    average = mean(values)
    return max(abs(value-average) for value in values)


def matrix_spectral_norm(matrix):
    return arb.svd_r(matrix, compute_uv=False)[0]


def maximum_matrix_entry(matrix):
    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows) for column in range(matrix.cols)
    )


def real_vector(raw):
    return tuple(arb.mpf(value) for value in raw)


def maximum_imaginary(action, gradient):
    return max(abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient))


print("Fixed canonical-target homotopy to the second dust tick", flush=True)
response_core = load_response_prefix()
models = response_core["models"]
core = response_core["core"]
action_and_gradient = core["action_and_gradient"]
branch_pass = core["branch_pass"]
ARB_L0_SQUARE = core["ARB_L0_SQUARE"]
ARB_RHO = core["ARB_RHO"]

check(
    "the imported action/evaluator core retains its six controls",
    response_core["tests"] == response_core["passed"] == 6,
)
check(
    "the fixed homotopy provenance, node count and parity pair are frozen",
    PRIOR_ART_COMMIT == "e760462"
    and PROTOCOL_COMMIT == "a2564d5"
    and START_RESULT_COMMIT == "46a7361"
    and DIRECT_BOUNDARY_COMMIT == "6346ad0"
    and NODE_COUNT == 32
    and set(models) == {"even", "odd"},
)

first_tick = json.loads(FIRST_TICK_ARTIFACT.read_text())
direct = json.loads(DIRECT_ARTIFACT.read_text())
gluing = json.loads(GLUING_ARTIFACT.read_text())
upstream_ok = bool(
    first_tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and first_tick["passed"] == first_tick["tests"] == 7
    and direct["outcome"] == "SECOND_TICK_NEWTON_OPEN"
    and direct["passed"] == direct["tests"] == 8
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing["passed"] == gluing["tests"]
)
check("the accepted start and direct-solver boundary authorize homotopy", upstream_ok)


first_states = {}
target_zero = {}
target_one = {}
target_bounds = {}
maps = {}
target_control_ok = True
for parity in models:
    solution = first_tick["solutions"][parity]
    first_states[parity] = tuple(arb.mpf(value) for value in solution["state"])
    target_zero[parity] = real_vector(solution["target_momentum"])
    stored_post = real_vector(solution["post_momentum"])
    mapping = tuple(
        gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"]
    )
    maps[parity] = mapping
    target_one[parity] = tuple(stored_post[index] for index in mapping)
    target_bounds[parity] = arb.mpf(solution["junction_bound"])
    target_control_ok &= bool(
        sorted(mapping) == list(range(30))
        and len(target_zero[parity]) == len(target_one[parity]) == 30
        and all(
            arb.isfinite(value)
            for value in target_zero[parity]+target_one[parity]
        )
        and target_bounds[parity] > 0
    )
check("the endpoint targets and vertex-derived maps are complete", target_control_ok)


def interpolated_target(parity, parameter):
    return tuple(
        (1-parameter)*left+parameter*right
        for left, right in zip(target_zero[parity], target_one[parity])
    )


def evaluate_state(parity, parameter, state):
    upper_log, rho_log = state
    first_scale_log, _ = first_states[parity]
    lower_log = parameter*first_scale_log
    q_old_value = arb.exp(2*lower_log)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*upper_log)*ARB_L0_SQUARE
    rho_value = ARB_RHO*arb.exp(rho_log)
    diagonal = arb.exp(lower_log+upper_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("homotopy magnitude left the positive domain")
    q_old = tuple([q_old_value]*30)
    internal = tuple([diagonal]*30+[rho_value]*5)
    q_new = tuple([q_new_value]*30)
    action, gradient, branch = action_and_gradient(
        models[parity], q_old, internal, q_new
    )
    local = tuple(arb.re(value) for value in gradient[30:65])
    pre = tuple(-arb.re(value) for value in gradient[:30])
    post = tuple(arb.re(value) for value in gradient[65:95])
    target = interpolated_target(parity, parameter)
    momentum_residual = tuple(
        value-wanted for value, wanted in zip(pre, target)
    )
    reduced = (mean(local[30:]), mean(momentum_residual))
    imaginary = maximum_imaginary(action, gradient)
    return {
        "parameter": parameter,
        "state": tuple(state),
        "lower_log": lower_log,
        "q_old": q_old_value,
        "q_new": q_new_value,
        "rho": rho_value,
        "diagonal": diagonal,
        "action": action,
        "gradient": gradient,
        "local": local,
        "pre": pre,
        "post": post,
        "target": target,
        "momentum_residual": momentum_residual,
        "reduced": reduced,
        "reduced_norm": infinity_norm(reduced),
        "branch": branch,
        "maximum_imaginary": imaginary,
        "branch_pass": branch_pass(branch, imaginary),
    }


def perturb(state, coordinate, amount):
    result = list(state)
    result[coordinate] += amount
    return tuple(result)


def calibrated_jacobian(parity, parameter, state):
    values = {}
    all_branch_pass = True
    minimum_minor = arb.inf
    minimum_argument = arb.inf
    maximum_imaginary_value = arb.mpf(0)
    for name, step in STEP_SETS.items():
        for coordinate in range(2):
            for direction in (-1, 1):
                evaluation = evaluate_state(
                    parity, parameter,
                    perturb(state, coordinate, direction*step),
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


def full_node_gate(parity, evaluation, jacobian):
    diagonal = evaluation["local"][:30]
    poles = evaluation["local"][30:]
    momentum = evaluation["momentum_residual"]
    return bool(
        evaluation["branch_pass"]
        and jacobian is not None
        and jacobian["resolved"]
        and infinity_norm(diagonal) < arb.mpf("1e-60")
        and infinity_norm(poles) < arb.mpf("1e-25")
        and infinity_norm(evaluation["local"]) < arb.mpf("1e-25")
        and spread(diagonal) < arb.mpf("1e-60")
        and spread(poles) < arb.mpf("1e-60")
        and vector_norm(momentum) <= target_bounds[parity]
        and spread(momentum) <= target_bounds[parity]
    )


start_records = {}
start_jacobians = {}
start_ok = bool(upstream_ok and target_control_ok)
if start_ok:
    for parity in models:
        start_records[parity] = evaluate_state(
            parity, arb.mpf(0), first_states[parity]
        )
        start_jacobians[parity] = calibrated_jacobian(
            parity, arb.mpf(0), first_states[parity]
        )
        start_ok &= full_node_gate(
            parity, start_records[parity], start_jacobians[parity]
        )
check(
    "lambda=0 reconstructs the complete accepted rank-two first-tick root",
    start_ok,
    "even ||F||={}, odd ||F||={}".format(
        text(start_records.get("even", {"reduced_norm": arb.inf})["reduced_norm"], 8),
        text(start_records.get("odd", {"reduced_norm": arb.inf})["reduced_norm"], 8),
    ),
)


paths = {}
all_branch_ok = start_ok
all_jacobians_resolved = start_ok
any_newton_failure = False
any_full_failure = False

if start_ok:
    for parity in models:
        print(f"  {parity}: advancing {NODE_COUNT} fixed homotopy nodes", flush=True)
        path = [{
            "node": 0,
            "parameter": arb.mpf(0),
            "state": first_states[parity],
            "evaluation": start_records[parity],
            "endpoint_jacobian": start_jacobians[parity],
            "attempts": [],
            "iterations": 0,
            "rejected_branch_trials": 0,
            "full_gate": True,
            "failure": None,
        }]
        state = first_states[parity]

        for node in range(1, NODE_COUNT+1):
            parameter = arb.mpf(node)/NODE_COUNT
            evaluation = evaluate_state(parity, parameter, state)
            all_branch_ok &= evaluation["branch_pass"]
            attempts = []
            iterations = 0
            rejected_branch_trials = 0
            failure = None
            converged = evaluation["reduced_norm"] < RESIDUAL_TOLERANCE

            while not converged and iterations < MAX_ITERATIONS:
                jacobian = calibrated_jacobian(parity, parameter, state)
                attempts.append({"state": state, "jacobian": jacobian})
                all_branch_ok &= jacobian["branch_pass"]
                all_jacobians_resolved &= jacobian["resolved"]
                if not jacobian["branch_pass"]:
                    failure = "JACOBIAN_BRANCH_FAILURE"
                    break
                if not jacobian["resolved"]:
                    failure = "JACOBIAN_ERROR_BAND_FAILURE"
                    break

                rhs = arb.matrix((
                    -evaluation["reduced"][0], -evaluation["reduced"][1]
                ))
                correction = arb.lu_solve(jacobian["operational"], rhs)
                accepted = None
                for exponent in range(MAX_DAMPING+1):
                    alpha = arb.mpf(2)**(-exponent)
                    trial_state = (
                        state[0]+alpha*correction[0],
                        state[1]+alpha*correction[1],
                    )
                    trial = evaluate_state(parity, parameter, trial_state)
                    if not trial["branch_pass"]:
                        rejected_branch_trials += 1
                        continue
                    if trial["reduced_norm"] <= (
                        1-alpha/4
                    )*evaluation["reduced_norm"]:
                        accepted = (trial_state, trial)
                        break
                if accepted is None:
                    failure = "NO_ARMIJO_DAMPING"
                    break
                state, evaluation = accepted
                iterations += 1
                converged = evaluation["reduced_norm"] < RESIDUAL_TOLERANCE

            if not converged and failure is None:
                failure = "ITERATION_LIMIT"
            endpoint_jacobian = None
            if converged:
                endpoint_jacobian = calibrated_jacobian(
                    parity, parameter, state
                )
                all_branch_ok &= endpoint_jacobian["branch_pass"]
                all_jacobians_resolved &= endpoint_jacobian["resolved"]
                if not endpoint_jacobian["branch_pass"]:
                    failure = "ENDPOINT_BRANCH_FAILURE"
                    converged = False
                elif not endpoint_jacobian["resolved"]:
                    failure = "ENDPOINT_JACOBIAN_ERROR_BAND_FAILURE"
                    converged = False

            node_full = bool(
                converged
                and full_node_gate(parity, evaluation, endpoint_jacobian)
            )
            any_newton_failure |= not converged
            any_full_failure |= bool(converged and not node_full)
            path.append({
                "node": node,
                "parameter": parameter,
                "state": state,
                "evaluation": evaluation,
                "endpoint_jacobian": endpoint_jacobian,
                "attempts": attempts,
                "iterations": iterations,
                "rejected_branch_trials": rejected_branch_trials,
                "full_gate": node_full,
                "failure": failure,
            })
            smin = (
                arb.nan if endpoint_jacobian is None
                else endpoint_jacobian["singular_values"][1]
            )
            print(
                "    node={:02d}/{} lambda={} iter={} b={} r={} ||F||={} smin={} full={}".format(
                    node, NODE_COUNT, text(parameter, 5), iterations,
                    text(state[0], 12), text(state[1], 12),
                    text(evaluation["reduced_norm"], 7), text(smin, 7),
                    node_full,
                ),
                flush=True,
            )
            if not node_full:
                break
        paths[parity] = path


paths_complete = bool(
    start_ok
    and all(len(path) == NODE_COUNT+1 for path in paths.values())
    and all(all(node["full_gate"] for node in path) for path in paths.values())
)
if paths_complete:
    maximum_state_difference = arb.mpf(0)
    for even_node, odd_node in zip(paths["even"], paths["odd"]):
        maximum_state_difference = max(
            maximum_state_difference,
            abs(even_node["state"][0]-odd_node["state"][0]),
            abs(even_node["state"][1]-odd_node["state"][1]),
        )
    even_endpoint = paths["even"][-1]["evaluation"]
    odd_endpoint = paths["odd"][-1]["evaluation"]
    endpoint_pre_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even_endpoint["pre"], odd_endpoint["pre"]
        )
    ))
    endpoint_post_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even_endpoint["post"], odd_endpoint["post"]
        )
    ))
    schedule_gate = bool(
        maximum_state_difference < arb.mpf("1e-25")
        and endpoint_pre_difference < arb.mpf("1e-22")
        and endpoint_post_difference < arb.mpf("1e-22")
    )
else:
    schedule_gate = False
    maximum_state_difference = None
    endpoint_pre_difference = endpoint_post_difference = None


if not start_ok or not all_branch_ok:
    outcome = "SECOND_TICK_HOMOTOPY_CONTROL_FAILED"
elif not all_jacobians_resolved:
    outcome = "SECOND_TICK_HOMOTOPY_JACOBIAN_OPEN"
elif any_newton_failure:
    outcome = "SECOND_TICK_HOMOTOPY_NEWTON_OPEN"
elif any_full_failure or not paths_complete:
    outcome = "SECOND_TICK_HOMOTOPY_FULL_SUBSTITUTION_FAILED"
elif not schedule_gate:
    outcome = "SECOND_TICK_HOMOTOPY_SCHEDULE_DEPENDENT"
else:
    first_scale_log = first_states["even"][0]
    endpoint_increment = paths["even"][-1]["state"][0]-first_scale_log
    if abs(endpoint_increment) <= arb.mpf("1e-20"):
        outcome = "SECOND_TICK_HOMOTOPY_STATIONARY"
    elif endpoint_increment < 0:
        outcome = "SECOND_TICK_HOMOTOPY_CONTINUED_CONTRACTION"
    else:
        outcome = "SECOND_TICK_HOMOTOPY_TURNED_TO_EXPANSION"

check(
    "all accepted homotopy evaluations retain the Lorentzian branch",
    all_branch_ok,
)
check(
    "the frozen hierarchy assigns exactly one homotopy outcome",
    outcome in {
        "SECOND_TICK_HOMOTOPY_CONTROL_FAILED",
        "SECOND_TICK_HOMOTOPY_JACOBIAN_OPEN",
        "SECOND_TICK_HOMOTOPY_NEWTON_OPEN",
        "SECOND_TICK_HOMOTOPY_FULL_SUBSTITUTION_FAILED",
        "SECOND_TICK_HOMOTOPY_SCHEDULE_DEPENDENT",
        "SECOND_TICK_HOMOTOPY_STATIONARY",
        "SECOND_TICK_HOMOTOPY_CONTINUED_CONTRACTION",
        "SECOND_TICK_HOMOTOPY_TURNED_TO_EXPANSION",
    },
    outcome,
)


def serialize_matrix(matrix):
    return [
        [text(matrix[row, column], 45) for column in range(matrix.cols)]
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
        "epsilon": text(record["epsilon"], 35),
        "singular_values": [
            text(value, 45) for value in record["singular_values"]
        ],
        "determinant": text(record["determinant"], 45),
        "condition": text(record["condition"], 35),
        "entry_pass": record["entry_pass"],
        "branch_pass": record["branch_pass"],
        "resolved": record["resolved"],
        "minimum_leading_minor": text(record["minimum_leading_minor"], 35),
        "minimum_angle_argument": text(record["minimum_angle_argument"], 35),
        "maximum_imaginary": text(record["maximum_imaginary"], 35),
        "maximum_d_op": text(maximum_matrix_entry(record["d_op"]), 25),
        "maximum_d_val": text(maximum_matrix_entry(record["d_val"]), 25),
        "maximum_d_cross": text(maximum_matrix_entry(record["d_cross"]), 25),
    }


def serialize_node(parity, node):
    evaluation = node["evaluation"]
    state = node["state"]
    lower_log = evaluation["lower_log"]
    first_rho_log = first_states[parity][1]
    return {
        "node": node["node"],
        "parameter": text(node["parameter"], 20),
        "state": [text(value, 55) for value in state],
        "scale_ratio_upper_lower": text(arb.exp(state[0]-lower_log), 55),
        "rho_ratio_from_first": text(arb.exp(state[1]-first_rho_log), 55),
        "tau_ratio_from_first": text(arb.exp((state[1]-first_rho_log)/2), 55),
        "iterations": node["iterations"],
        "rejected_branch_trials": node["rejected_branch_trials"],
        "failure": node["failure"],
        "full_gate": node["full_gate"],
        "reduced_residual": [text(value, 45) for value in evaluation["reduced"]],
        "reduced_residual_norm": text(evaluation["reduced_norm"], 35),
        "diagonal_residuals": [
            text(value, 35) for value in evaluation["local"][:30]
        ],
        "pole_residuals": [
            text(value, 35) for value in evaluation["local"][30:]
        ],
        "momentum_residuals": [
            text(value, 35) for value in evaluation["momentum_residual"]
        ],
        "junction_norm": text(vector_norm(evaluation["momentum_residual"]), 35),
        "pre_momentum": [text(value, 45) for value in evaluation["pre"]],
        "post_momentum": [text(value, 45) for value in evaluation["post"]],
        "target_momentum": [text(value, 45) for value in evaluation["target"]],
        "endpoint_jacobian": serialize_jacobian(node["endpoint_jacobian"]),
        "attempts": [
            {
                "state": [text(value, 45) for value in attempt["state"]],
                "jacobian": serialize_jacobian(attempt["jacobian"]),
            }
            for attempt in node["attempts"]
        ],
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "start_result_commit": START_RESULT_COMMIT,
    "direct_boundary_commit": DIRECT_BOUNDARY_COMMIT,
    "node_count": NODE_COUNT,
    "steps": {name: text(value, 20) for name, value in STEP_SETS.items()},
    "maps": {parity: list(mapping) for parity, mapping in maps.items()},
    "paths": {
        parity: [serialize_node(parity, node) for node in path]
        for parity, path in paths.items()
    },
    "schedule_gate": {
        "maximum_state_difference": (
            None if maximum_state_difference is None
            else text(maximum_state_difference, 35)
        ),
        "endpoint_pre_difference": (
            None if endpoint_pre_difference is None
            else text(endpoint_pre_difference, 35)
        ),
        "endpoint_post_difference": (
            None if endpoint_post_difference is None
            else text(endpoint_post_difference, 35)
        ),
        "passed": schedule_gate,
    },
    "classification": {
        "connected_second_homothetic_tick": bool(paths_complete and schedule_gate),
        "intermediate_nodes_are_physical_ticks": False,
        "absolute_time_unit": False,
        "anisotropic_stability": "OPEN",
        "refinement_stability": "OPEN",
        "label": "STRUCTURAL / candidate pseudo-constraint",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

if paths_complete:
    for parity in models:
        endpoint = paths[parity][-1]
        b, r = endpoint["state"]
        a1, r1 = first_states[parity]
        print(
            "  {} endpoint: u={} v={} L2/L1={} tau2/tau1={}".format(
                parity, text(b-a1, 16), text(r-r1, 16),
                text(arb.exp(b-a1), 16), text(arb.exp((r-r1)/2), 16),
            ),
            flush=True,
        )
print(f"\nOUTCOME: {outcome}", flush=True)
print(f"RESULT: {passed}/{tests}", flush=True)
if passed != tests:
    raise SystemExit(1)

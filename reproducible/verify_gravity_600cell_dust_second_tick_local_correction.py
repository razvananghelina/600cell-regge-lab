#!/usr/bin/env python3
"""Local canonical correction of the contracting second-tick root.

Prior-art commit: fcc4d7c.
Protocol commit: 54dd336.
"""

import ast
import contextlib
import hashlib
import io
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
RESPONSE_SOURCE = HERE / "verify_gravity_600cell_dust_homothetic_mass_conservation.py"
ROOT_ARTIFACT = HERE / "gravity_600cell_dust_stationary_root_enumeration.json"
TICK_ARTIFACT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
COMPARISON_ARTIFACT = HERE / "gravity_600cell_dust_second_tick_stationary_target.json"
OUTPUT = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
PRIOR_ART_COMMIT = "fcc4d7c"
PROTOCOL_COMMIT = "54dd336"
ROOT_SHA256 = "0ec5ba520ea25b39dd6cfd3c349d49fe480df2abee359854e1316b5af4d9fa2f"
TICK_SHA256 = "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9"
GLUING_SHA256 = "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77"
COMPARISON_SHA256 = "1769f9677e57d83a7a87ba200acdf6519f808d77d5bfddd7bc595fd6a23b27e7"
STEP_SETS = {
    "operational_primary": arb.mpf("1e-20"),
    "operational_shadow": arb.mpf("1e-15"),
    "validation_primary": arb.mpf("3e-20"),
    "validation_shadow": arb.mpf("3e-15"),
}
DPS = 100
arb.mp.dps = DPS
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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    return arb.svd_r(matrix, compute_uv=False)[0]


def maximum_matrix_entry(matrix):
    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows) for column in range(matrix.cols)
    )


def maximum_imaginary(action, gradient):
    return max(abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient))


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
        "__name__": "second_tick_local_imported_response_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(RESPONSE_SOURCE), "exec"), namespace)
    return namespace


print("Second homothetic tick: local canonical correction", flush=True)
response_core = load_response_prefix()
models = response_core["models"]
core = response_core["core"]
action_and_gradient = core["action_and_gradient"]
branch_pass = core["branch_pass"]
ARB_L0_SQUARE = core["ARB_L0_SQUARE"]
ARB_RHO = core["ARB_RHO"]

check(
    "the imported complete-action evaluator retains its six controls",
    response_core["tests"] == response_core["passed"] == 6,
)

roots = json.loads(ROOT_ARTIFACT.read_text())
tick = json.loads(TICK_ARTIFACT.read_text())
gluing = json.loads(GLUING_ARTIFACT.read_text())
comparison = json.loads(COMPARISON_ARTIFACT.read_text())
hashes_ok = bool(
    digest(ROOT_ARTIFACT) == ROOT_SHA256
    and digest(TICK_ARTIFACT) == TICK_SHA256
    and digest(GLUING_ARTIFACT) == GLUING_SHA256
    and digest(COMPARISON_ARTIFACT) == COMPARISON_SHA256
)
provenance_ok = bool(
    PRIOR_ART_COMMIT == "fcc4d7c"
    and PROTOCOL_COMMIT == "54dd336"
    and roots["outcome"] == "STATIONARY_ROOTS_ENUMERATED"
    and roots["passed"] == roots["tests"] == 5
    and tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and tick["passed"] == tick["tests"] == 7
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing["passed"] == gluing["tests"] == 25
    and comparison["outcome"] == "STATIONARY_SECOND_TICK_NO_HIT"
    and comparison["passed"] == comparison["tests"] == 5
    and comparison["hit_count"] == 0
)
check(
    "the four frozen inputs and target-before-correction provenance pass",
    hashes_ok and provenance_ok,
    "root={} tick={} glue={} comparison={}".format(
        digest(ROOT_ARTIFACT), digest(TICK_ARTIFACT),
        digest(GLUING_ARTIFACT), digest(COMPARISON_ARTIFACT),
    ),
)

A1 = arb.mpf(roots["geometry"]["a1"])
R1 = arb.mpf(roots["geometry"]["r1"])
contracting_root = roots["roots"][0]
SEED = (arb.mpf(contracting_root["upper_log"]), R1)
seed_selection_ok = bool(
    contracting_root["index"] == 0
    and contracting_root["kind"] == "sign_bracket"
    and SEED[0] < A1
    and roots["roots"][1]["kind"] == "node"
    and arb.mpf(roots["roots"][1]["upper_log"]) == 0
)

targets = {}
bounds = {}
maps = {}
target_ok = True
for parity in ("even", "odd"):
    mapping = tuple(gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"])
    post = tuple(arb.mpf(value) for value in tick["solutions"][parity]["post_momentum"])
    maps[parity] = mapping
    targets[parity] = tuple(post[index] for index in mapping)
    bounds[parity] = arb.mpf(tick["solutions"][parity]["junction_bound"])
    target_ok &= bool(
        sorted(mapping) == list(range(30))
        and len(post) == len(targets[parity]) == 30
        and bounds[parity] > 0
    )
check(
    "the sole contracting seed and both complete canonical targets are fixed",
    seed_selection_ok and target_ok,
    f"seed=({text(SEED[0], 30)},{text(SEED[1], 30)})",
)


def evaluate_reduced(parity, state):
    upper_log, rho_log = state
    q_old_value = arb.exp(2*A1)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*upper_log)*ARB_L0_SQUARE
    rho_value = ARB_RHO*arb.exp(rho_log)
    diagonal = arb.exp(A1+upper_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("second-tick homothetic magnitude left positive domain")
    action, gradient, branch = action_and_gradient(
        models[parity], tuple([q_old_value]*30),
        tuple([diagonal]*30+[rho_value]*5), tuple([q_new_value]*30),
    )
    local = tuple(arb.re(value) for value in gradient[30:65])
    pre = tuple(-arb.re(value) for value in gradient[:30])
    post = tuple(arb.re(value) for value in gradient[65:95])
    momentum_residual = tuple(
        value-target for value, target in zip(pre, targets[parity])
    )
    reduced = (mean(local[30:]), mean(momentum_residual))
    imaginary = maximum_imaginary(action, gradient)
    return {
        "state": tuple(state),
        "action": action,
        "gradient": gradient,
        "local": local,
        "pre": pre,
        "post": post,
        "momentum_residual": momentum_residual,
        "reduced": reduced,
        "reduced_norm": infinity_norm(reduced),
        "branch": branch,
        "maximum_imaginary": imaginary,
        "branch_pass": branch_pass(branch, imaginary),
    }


seed_records = {parity: evaluate_reduced(parity, SEED) for parity in models}
stored_p = arb.mpf(contracting_root["P"])
seed_control_ok = all(
    record["branch_pass"]
    and infinity_norm(record["local"]) < arb.mpf("1e-25")
    and abs(mean(record["pre"])-stored_p) < arb.mpf("1e-25")
    for record in seed_records.values()
)
seed_control_ok &= bool(
    infinity_norm(tuple(
        left-right for left, right in zip(
            seed_records["even"]["local"], seed_records["odd"]["local"]
        )
    )) < arb.mpf("1e-24")
)
check(
    "a fresh complete-action evaluation reproduces the committed sole seed",
    seed_control_ok,
    "F_even=({}, {}), F_odd=({}, {})".format(
        text(seed_records["even"]["reduced"][0], 10),
        text(seed_records["even"]["reduced"][1], 10),
        text(seed_records["odd"]["reduced"][0], 10),
        text(seed_records["odd"]["reduced"][1], 10),
    ),
)


def perturb(state, coordinate, amount):
    changed = list(state)
    changed[coordinate] += amount
    return tuple(changed)


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
                    minimum_minor, evaluation["branch"]["minimum_leading_minor"]
                )
                minimum_argument = min(
                    minimum_argument, evaluation["branch"]["minimum_argument"]
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
all_solver_branch_ok = True
all_jacobians_resolved = True
for parity in models:
    print(f"  {parity}: deterministic local Newton correction", flush=True)
    state = SEED
    evaluation = seed_records[parity]
    history = []
    jacobian_attempts = []
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
            "    J{} singular=({}, {}) epsilon={} resolved={}".format(
                len(jacobian_attempts), text(jacobian["singular_values"][0], 9),
                text(jacobian["singular_values"][1], 9),
                text(jacobian["epsilon"], 9), jacobian["resolved"],
            ), flush=True,
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
            all_solver_branch_ok &= trial["branch_pass"]
            if not trial["branch_pass"]:
                failure = "TRIAL_BRANCH_FAILURE"
                break
            if trial["reduced_norm"] <= (
                1-alpha/4
            )*evaluation["reduced_norm"]:
                accepted = (alpha, trial_state, trial)
                break
        if failure == "TRIAL_BRANCH_FAILURE":
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
            "    iter={} alpha={} b={} r={} ||F||={}".format(
                accepted_iterations, text(alpha, 4), text(state[0], 16),
                text(state[1], 16), text(evaluation["reduced_norm"], 9),
            ), flush=True,
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
    momentum = evaluation["momentum_residual"]
    full_gate = bool(
        converged
        and infinity_norm(diagonal) < arb.mpf("1e-60")
        and infinity_norm(poles) < arb.mpf("1e-25")
        and infinity_norm(evaluation["local"]) < arb.mpf("1e-25")
        and spread(diagonal) < arb.mpf("1e-60")
        and spread(poles) < arb.mpf("1e-60")
        and vector_norm(momentum) <= bounds[parity]
        and spread(momentum) <= bounds[parity]
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
        "full_gate": full_gate,
    }

if all(record["converged"] for record in solve_records.values()):
    even = solve_records["even"]
    odd = solve_records["odd"]
    b_difference = abs(even["state"][0]-odd["state"][0])
    r_difference = abs(even["state"][1]-odd["state"][1])
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
        b_difference < arb.mpf("1e-25")
        and r_difference < arb.mpf("1e-25")
        and pre_difference < arb.mpf("1e-22")
        and post_difference < arb.mpf("1e-22")
    )
else:
    parity_gate = False
    b_difference = r_difference = pre_difference = post_difference = None

controls_ok = bool(
    hashes_ok and provenance_ok and seed_selection_ok and target_ok
    and seed_control_ok and all_solver_branch_ok
)
newton_ok = all(record["converged"] for record in solve_records.values())
full_ok = all(record["full_gate"] for record in solve_records.values())
if not controls_ok:
    outcome = "SECOND_TICK_LOCAL_CONTROL_FAILED"
elif not all_jacobians_resolved:
    outcome = "SECOND_TICK_LOCAL_JACOBIAN_OPEN"
elif not newton_ok:
    outcome = "SECOND_TICK_LOCAL_NEWTON_OPEN"
elif not full_ok:
    outcome = "SECOND_TICK_LOCAL_FULL_SUBSTITUTION_FAILED"
elif not parity_gate:
    outcome = "SECOND_TICK_LOCAL_SCHEDULE_DEPENDENT"
else:
    outcome = "SECOND_HOMOTHETIC_TICK_ACCEPTED"

check("every evaluated state retains the Lorentzian branch", all_solver_branch_ok)
check(
    "the frozen hierarchy assigns one local second-tick outcome",
    outcome in {
        "SECOND_TICK_LOCAL_CONTROL_FAILED",
        "SECOND_TICK_LOCAL_JACOBIAN_OPEN",
        "SECOND_TICK_LOCAL_NEWTON_OPEN",
        "SECOND_TICK_LOCAL_FULL_SUBSTITUTION_FAILED",
        "SECOND_TICK_LOCAL_SCHEDULE_DEPENDENT",
        "SECOND_HOMOTHETIC_TICK_ACCEPTED",
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
        "singular_values": [text(value, 50) for value in record["singular_values"]],
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
    state = result["state"]
    evaluation = result["evaluation"]
    u2 = state[0]-A1
    v2 = state[1]-R1
    serialized_solutions[parity] = {
        "seed": [text(value, 60) for value in SEED],
        "state_absolute": [text(value, 60) for value in state],
        "relative_state": [text(u2, 60), text(v2, 60)],
        "scale_ratio_L2_over_L1": text(arb.exp(u2), 60),
        "rho_ratio_rho2_over_rho1": text(arb.exp(v2), 60),
        "tau_ratio_tau2_over_tau1": text(arb.exp(v2/2), 60),
        "scale_increment_ratio_u2_over_a1": text(u2/A1, 60),
        "iterations": result["iterations"],
        "converged": result["converged"],
        "failure": result["failure"],
        "full_gate": result["full_gate"],
        "reduced_residual": [text(value, 50) for value in evaluation["reduced"]],
        "reduced_residual_norm": text(evaluation["reduced_norm"], 40),
        "diagonal_residuals": [text(value, 40) for value in evaluation["local"][:30]],
        "pole_residuals": [text(value, 40) for value in evaluation["local"][30:]],
        "pre_momentum": [text(value, 50) for value in evaluation["pre"]],
        "post_momentum": [text(value, 50) for value in evaluation["post"]],
        "target_momentum": [text(value, 50) for value in targets[parity]],
        "junction_residual": [text(value, 40) for value in evaluation["momentum_residual"]],
        "junction_norm": text(vector_norm(evaluation["momentum_residual"]), 40),
        "junction_bound": text(bounds[parity], 40),
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
    "input_sha256": {
        "stationary_roots": digest(ROOT_ARTIFACT),
        "first_tick": digest(TICK_ARTIFACT),
        "gluing": digest(GLUING_ARTIFACT),
        "stationary_target_comparison": digest(COMPARISON_ARTIFACT),
    },
    "fixed_mass": True,
    "mass_recomputed_from_later_scale": False,
    "sole_seed_root_index": 0,
    "steps": {name: text(value, 20) for name, value in STEP_SETS.items()},
    "solutions": serialized_solutions,
    "parity_gate": {
        "b_difference": None if b_difference is None else text(b_difference, 40),
        "r_difference": None if r_difference is None else text(r_difference, 40),
        "pre_momentum_infinity_difference": (
            None if pre_difference is None else text(pre_difference, 40)
        ),
        "post_momentum_infinity_difference": (
            None if post_difference is None else text(post_difference, 40)
        ),
        "passed": parity_gate,
    },
    "classification": {
        "local_second_tick": outcome == "SECOND_HOMOTHETIC_TICK_ACCEPTED",
        "arbitrary_tick_recurrence": "OPEN",
        "absolute_clock": "OPEN",
        "refinement_stability": "OPEN",
        "label": "DERIVED COMPUTATIONAL if accepted",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

for parity, result in solve_records.items():
    state = result["state"]
    u2 = state[0]-A1
    v2 = state[1]-R1
    print(
        "  {}: b={} r={} u2={} v2={} L2/L1={} tau2/tau1={} ||F||={} ||J||={}".format(
            parity, text(state[0], 17), text(state[1], 17), text(u2, 17),
            text(v2, 17), text(arb.exp(u2), 17), text(arb.exp(v2/2), 17),
            text(result["evaluation"]["reduced_norm"], 9),
            text(vector_norm(result["evaluation"]["momentum_residual"]), 9),
        ), flush=True,
    )
print(f"OUTCOME: {outcome}", flush=True)
print(f"Tests passed: {passed}/{tests}", flush=True)
raise SystemExit(0 if passed == tests else 1)

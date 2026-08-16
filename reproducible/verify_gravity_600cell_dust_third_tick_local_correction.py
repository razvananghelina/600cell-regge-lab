#!/usr/bin/env python3
"""Cleanly preregistered local canonical correction of the third dust tick.

Prior-art commit: 7b9a676.
Protocol commit: 1782b29.
"""

import ast
import contextlib
import hashlib
import io
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
BASE_SOLVER_SOURCE = HERE / "verify_gravity_600cell_dust_second_tick_local_correction.py"
RESPONSE_SOURCE = HERE / "verify_gravity_600cell_dust_homothetic_mass_conservation.py"
ROOT_ARTIFACT = HERE / "gravity_600cell_dust_third_tick_stationary_roots.json"
SECOND_TICK_ARTIFACT = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
COMPARISON_ARTIFACT = HERE / "gravity_600cell_dust_third_tick_stationary_target.json"
OUTPUT = HERE / "gravity_600cell_dust_third_tick_local_correction.json"
PRIOR_ART_COMMIT = "7b9a676"
PROTOCOL_COMMIT = "1782b29"
BASE_SOLVER_SHA256 = "cef59fa0bc3a1c8fa3be0193234371b7dda303a0ec72683ddcdd88bcb40f3725"
ROOT_SHA256 = "02d4589a7df0851c67a31fc0a41c5ef8851a82c758214c1c5e8729afddfe479f"
SECOND_TICK_SHA256 = "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70"
GLUING_SHA256 = "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77"
COMPARISON_SHA256 = "4d1f81dafcab9d3aa40ff08fdaaad90b80235809dd32becd790bdee1704ab6cf"
R1_TEXT = "-0.00000355925313517063343725030533963917396571974345422547402551491"
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


def file_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_solver_functions():
    """Load only audited function definitions; never execute the old solve."""
    tree = ast.parse(BASE_SOLVER_SOURCE.read_text(), filename=str(BASE_SOLVER_SOURCE))
    body = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    module = ast.Module(body=body, type_ignores=[])
    namespace = {
        "__file__": str(BASE_SOLVER_SOURCE),
        "__name__": "third_tick_imported_solver_functions",
        "arb": arb,
        "ast": ast,
        "contextlib": contextlib,
        "io": io,
        "json": json,
        "Path": Path,
        "RESPONSE_SOURCE": RESPONSE_SOURCE,
        "STEP_SETS": STEP_SETS,
        "ARITHMETIC_FLOOR": ARITHMETIC_FLOOR,
        "ENTRY_FACTOR": ENTRY_FACTOR,
        "NONZERO_FACTOR": NONZERO_FACTOR,
        "tests": 0,
        "passed": 0,
    }
    exec(compile(module, str(BASE_SOLVER_SOURCE), "exec"), namespace)
    return namespace


solver = load_solver_functions()
check = solver["check"]
text = solver["text"]
mean = solver["mean"]
vector_norm = solver["vector_norm"]
infinity_norm = solver["infinity_norm"]
spread = solver["spread"]
maximum_imaginary = solver["maximum_imaginary"]
load_response_prefix = solver["load_response_prefix"]
calibrated_jacobian = solver["calibrated_jacobian"]
serialize_jacobian = solver["serialize_jacobian"]

print("Third homothetic tick: local canonical correction", flush=True)
response_core = load_response_prefix()
models = response_core["models"]
core = response_core["core"]
action_and_gradient = core["action_and_gradient"]
branch_pass = core["branch_pass"]
ARB_L0_SQUARE = core["ARB_L0_SQUARE"]
ARB_RHO = core["ARB_RHO"]
solver["models"] = models

check(
    "the imported complete-action evaluator retains its six controls",
    response_core["tests"] == response_core["passed"] == 6,
)

roots = json.loads(ROOT_ARTIFACT.read_text())
second_tick = json.loads(SECOND_TICK_ARTIFACT.read_text())
gluing = json.loads(GLUING_ARTIFACT.read_text())
comparison = json.loads(COMPARISON_ARTIFACT.read_text())
hashes_ok = bool(
    file_digest(BASE_SOLVER_SOURCE) == BASE_SOLVER_SHA256
    and file_digest(ROOT_ARTIFACT) == ROOT_SHA256
    and file_digest(SECOND_TICK_ARTIFACT) == SECOND_TICK_SHA256
    and file_digest(GLUING_ARTIFACT) == GLUING_SHA256
    and file_digest(COMPARISON_ARTIFACT) == COMPARISON_SHA256
)
provenance_ok = bool(
    PRIOR_ART_COMMIT == "7b9a676"
    and PROTOCOL_COMMIT == "1782b29"
    and roots["outcome"] == "THIRD_TICK_STATIONARY_ROOTS_ENUMERATED"
    and roots["passed"] == roots["tests"] == 5
    and roots["target_parsed"] is False
    and roots["candidate_count"] == 2
    and second_tick["outcome"] == "SECOND_HOMOTHETIC_TICK_ACCEPTED"
    and second_tick["passed"] == second_tick["tests"] == 6
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing["passed"] == gluing["tests"] == 25
    and comparison["outcome"] == "STATIONARY_THIRD_TICK_NO_HIT"
    and comparison["passed"] == comparison["tests"] == 5
    and comparison["hit_count"] == 0
)
check(
    "all frozen artifacts and the audited solver-function source pass",
    hashes_ok and provenance_ok,
    "roots={} second={} comparison={} base={}".format(
        file_digest(ROOT_ARTIFACT), file_digest(SECOND_TICK_ARTIFACT),
        file_digest(COMPARISON_ARTIFACT), file_digest(BASE_SOLVER_SOURCE),
    ),
)

A1 = arb.mpf(roots["geometry"]["A1"])
B2 = arb.mpf(roots["geometry"]["B2"])
R2 = arb.mpf(roots["geometry"]["R2"])
R1 = arb.mpf(R1_TEXT)
contracting_root = roots["roots"][0]
C_SEED = arb.mpf(contracting_root["parities"]["even"]["upper_log"])
SEED = (C_SEED, R2)
seed_selection_ok = bool(
    contracting_root["index"] == 0
    and contracting_root["kind"] == "sign_bracket"
    and C_SEED < B2
    and contracting_root["parities"]["even"]["full_gate"]
    and contracting_root["parities"]["even"]["derivative_gate"]
    and roots["roots"][1]["kind"] == "node_cluster"
)

targets = {}
bounds = {}
maps = {}
target_ok = True
for parity in ("even", "odd"):
    mapping = tuple(gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"])
    post = tuple(
        arb.mpf(value)
        for value in second_tick["solutions"][parity]["post_momentum"]
    )
    maps[parity] = mapping
    targets[parity] = tuple(post[index] for index in mapping)
    bounds[parity] = arb.mpf(second_tick["solutions"][parity]["junction_bound"])
    target_ok &= bool(
        sorted(mapping) == list(range(30)) and len(post) == 30
        and len(targets[parity]) == 30 and bounds[parity] > 0
    )
check(
    "the sole forward seed and both complete canonical targets are fixed",
    seed_selection_ok and target_ok,
    f"seed=({text(SEED[0], 30)},{text(SEED[1], 30)})",
)


def evaluate_reduced(parity, state):
    upper_log, rho_log = state
    q_old_value = arb.exp(2*B2)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*upper_log)*ARB_L0_SQUARE
    rho_value = ARB_RHO*arb.exp(rho_log)
    diagonal = arb.exp(B2+upper_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("third-tick corrected state left positive domain")
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


solver["evaluate_reduced"] = evaluate_reduced
seed_records = {parity: evaluate_reduced(parity, SEED) for parity in models}
seed_control_ok = True
for parity, record in seed_records.items():
    stored = contracting_root["parities"][parity]["evaluation"]
    seed_control_ok &= bool(
        record["branch_pass"]
        and infinity_norm(record["local"]) < arb.mpf("1e-25")
        and abs(mean(record["pre"])-arb.mpf(stored["P"])) < arb.mpf("1e-25")
        and abs(mean(record["local"][30:])-arb.mpf(stored["G"])) < arb.mpf("1e-25")
    )
seed_control_ok &= bool(
    infinity_norm(tuple(
        left-right for left, right in zip(
            seed_records["even"]["local"], seed_records["odd"]["local"]
        )
    )) < arb.mpf("1e-24")
)
check(
    "fresh complete-action evaluations reproduce the committed sole seed",
    seed_control_ok,
    "F_even=({}, {}), F_odd=({}, {})".format(
        text(seed_records["even"]["reduced"][0], 10),
        text(seed_records["even"]["reduced"][1], 10),
        text(seed_records["odd"]["reduced"][0], 10),
        text(seed_records["odd"]["reduced"][1], 10),
    ),
)

solve_records = {}
all_solver_branch_ok = True
all_jacobians_resolved = True
for parity in models:
    print(f"  {parity}: deterministic third-tick Newton correction", flush=True)
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
            "    iter={} alpha={} C={} R={} ||F||={}".format(
                accepted_iterations, text(alpha, 4), text(state[0], 17),
                text(state[1], 17), text(evaluation["reduced_norm"], 9),
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
    c_difference = abs(even["state"][0]-odd["state"][0])
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
        c_difference < arb.mpf("1e-25")
        and r_difference < arb.mpf("1e-25")
        and pre_difference < arb.mpf("1e-22")
        and post_difference < arb.mpf("1e-22")
    )
else:
    parity_gate = False
    c_difference = r_difference = pre_difference = post_difference = None

controls_ok = bool(
    hashes_ok and provenance_ok and seed_selection_ok and target_ok
    and seed_control_ok and all_solver_branch_ok
)
newton_ok = all(record["converged"] for record in solve_records.values())
full_ok = all(record["full_gate"] for record in solve_records.values())
if not controls_ok:
    outcome = "THIRD_TICK_LOCAL_CONTROL_FAILED"
elif not all_jacobians_resolved:
    outcome = "THIRD_TICK_LOCAL_JACOBIAN_OPEN"
elif not newton_ok:
    outcome = "THIRD_TICK_LOCAL_NEWTON_OPEN"
elif not full_ok:
    outcome = "THIRD_TICK_LOCAL_FULL_SUBSTITUTION_FAILED"
elif not parity_gate:
    outcome = "THIRD_TICK_LOCAL_SCHEDULE_DEPENDENT"
else:
    outcome = "THIRD_HOMOTHETIC_TICK_ACCEPTED"

check("every evaluated state retains the Lorentzian branch", all_solver_branch_ok)
check(
    "the frozen hierarchy assigns one third-tick outcome",
    outcome in {
        "THIRD_TICK_LOCAL_CONTROL_FAILED",
        "THIRD_TICK_LOCAL_JACOBIAN_OPEN",
        "THIRD_TICK_LOCAL_NEWTON_OPEN",
        "THIRD_TICK_LOCAL_FULL_SUBSTITUTION_FAILED",
        "THIRD_TICK_LOCAL_SCHEDULE_DEPENDENT",
        "THIRD_HOMOTHETIC_TICK_ACCEPTED",
    },
    outcome,
)

serialized_solutions = {}
for parity, result in solve_records.items():
    state = result["state"]
    evaluation = result["evaluation"]
    u3 = state[0]-B2
    v3 = state[1]-R2
    serialized_solutions[parity] = {
        "seed": [text(value, 60) for value in SEED],
        "state_absolute": [text(value, 60) for value in state],
        "relative_state": [text(u3, 60), text(v3, 60)],
        "scale_ratio_L3_over_L2": text(arb.exp(u3), 60),
        "rho_ratio_rho3_over_rho2": text(arb.exp(v3), 60),
        "tau_ratio_tau3_over_tau2": text(arb.exp(v3/2), 60),
        "u3_over_A1": text(u3/A1, 60),
        "v3_over_R1": text(v3/R1, 60),
        "C_over_A1": text(state[0]/A1, 60),
        "R_over_R1": text(state[1]/R1, 60),
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
    "clean_new_state_preregistration": True,
    "input_sha256": {
        "base_solver_functions": file_digest(BASE_SOLVER_SOURCE),
        "stationary_roots": file_digest(ROOT_ARTIFACT),
        "second_tick": file_digest(SECOND_TICK_ARTIFACT),
        "gluing": file_digest(GLUING_ARTIFACT),
        "target_comparison": file_digest(COMPARISON_ARTIFACT),
    },
    "fixed_mass": True,
    "mass_recomputed_from_later_scale": False,
    "sole_seed_root_index": 0,
    "steps": {name: text(value, 20) for name, value in STEP_SETS.items()},
    "solutions": serialized_solutions,
    "parity_gate": {
        "C_difference": None if c_difference is None else text(c_difference, 40),
        "R_difference": None if r_difference is None else text(r_difference, 40),
        "pre_momentum_infinity_difference": (
            None if pre_difference is None else text(pre_difference, 40)
        ),
        "post_momentum_infinity_difference": (
            None if post_difference is None else text(post_difference, 40)
        ),
        "passed": parity_gate,
    },
    "classification": {
        "local_third_tick": outcome == "THIRD_HOMOTHETIC_TICK_ACCEPTED",
        "integer_sequence": "PATTERN",
        "arbitrary_tick_recurrence": "OPEN",
        "absolute_clock": "OPEN",
        "refinement_stability": "OPEN",
    },
    "outcome": outcome,
    "tests": solver["tests"],
    "passed": solver["passed"],
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

for parity, result in solve_records.items():
    state = result["state"]
    u3 = state[0]-B2
    v3 = state[1]-R2
    print(
        "  {}: C={} R={} u3/A1={} v3/R1={} L3/L2={} tau3/tau2={} ||F||={} ||J||={}".format(
            parity, text(state[0], 18), text(state[1], 18), text(u3/A1, 16),
            text(v3/R1, 16), text(arb.exp(u3), 17), text(arb.exp(v3/2), 17),
            text(result["evaluation"]["reduced_norm"], 9),
            text(vector_norm(result["evaluation"]["momentum_residual"]), 9),
        ), flush=True,
    )
print(f"OUTCOME: {outcome}", flush=True)
print(f"Tests passed: {solver['passed']}/{solver['tests']}", flush=True)
raise SystemExit(0 if solver["passed"] == solver["tests"] else 1)

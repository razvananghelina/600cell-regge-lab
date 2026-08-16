#!/usr/bin/env python3
"""Weak-lapse jet of the three-step canonical 600-cell dust map.

Prior-art commit: 90757a6.
Protocol commit: 1c26fdb.
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
STATIC_ARTIFACT = HERE / "gravity_600cell_dust_regular_lapse_identity.json"
TICK1_ARTIFACT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
TICK2_ARTIFACT = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
TICK3_ARTIFACT = HERE / "gravity_600cell_dust_third_tick_local_correction.json"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_weak_lapse_recurrence.json"
PRIOR_ART_COMMIT = "90757a6"
PROTOCOL_COMMIT = "1c26fdb"
BASE_SOLVER_SHA256 = "cef59fa0bc3a1c8fa3be0193234371b7dda303a0ec72683ddcdd88bcb40f3725"
INPUT_HASHES = {
    "static": "5079428fade247f730ebc07e5e2eae388b48045cd5201e84afb3186bfc248a51",
    "tick1": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "tick2": "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70",
    "tick3": "ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0",
    "gluing": "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77",
}
STEP_SETS = {
    "operational_primary": arb.mpf("1e-20"),
    "operational_shadow": arb.mpf("1e-15"),
    "validation_primary": arb.mpf("3e-20"),
    "validation_shadow": arb.mpf("3e-15"),
}
LAMBDA_TEXTS = ("0.5", "0.25", "0.125")
DPS = 100
arb.mp.dps = DPS
ARITHMETIC_FLOOR = arb.mpf("1e-60")
ENTRY_FACTOR = arb.mpf(10)
NONZERO_FACTOR = arb.mpf(100)
RESIDUAL_TOLERANCE = arb.mpf("1e-25")
JUNCTION_TOLERANCE = arb.mpf("1e-24")
ASYMPTOTIC_FLOOR = arb.mpf("1e-40")
MAX_ITERATIONS = 8
MAX_DAMPING = 10


def file_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_solver_functions():
    tree = ast.parse(BASE_SOLVER_SOURCE.read_text(), filename=str(BASE_SOLVER_SOURCE))
    body = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "__file__": str(BASE_SOLVER_SOURCE),
        "__name__": "weak_lapse_imported_solver_functions",
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
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(BASE_SOLVER_SOURCE), "exec"),
        namespace,
    )
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

print("Weak-lapse recurrence audit of the canonical dust map", flush=True)
response_core = load_response_prefix()
models = response_core["models"]
core = response_core["core"]
action_and_gradient = core["action_and_gradient"]
branch_pass = core["branch_pass"]
ARB_L0 = core["ARB_L0"]
ARB_L0_SQUARE = core["ARB_L0_SQUARE"]
ARB_RHO = core["ARB_RHO"]
ARB_TAU = core["ARB_TAU"]
ARB_EPSILON_3 = core["ARB_EPSILON_3"]
K0 = ARB_EPSILON_3*ARB_L0*ARB_TAU/4
solver["models"] = models

check(
    "the imported complete-action evaluator retains its six controls",
    response_core["tests"] == response_core["passed"] == 6,
)

paths = {
    "static": STATIC_ARTIFACT,
    "tick1": TICK1_ARTIFACT,
    "tick2": TICK2_ARTIFACT,
    "tick3": TICK3_ARTIFACT,
    "gluing": GLUING_ARTIFACT,
}
hashes = {name: file_digest(path) for name, path in paths.items()}
static_artifact = json.loads(STATIC_ARTIFACT.read_text())
tick1 = json.loads(TICK1_ARTIFACT.read_text())
tick2 = json.loads(TICK2_ARTIFACT.read_text())
tick3 = json.loads(TICK3_ARTIFACT.read_text())
gluing = json.loads(GLUING_ARTIFACT.read_text())
hashes_ok = bool(
    file_digest(BASE_SOLVER_SOURCE) == BASE_SOLVER_SHA256
    and hashes == INPUT_HASHES
)
provenance_ok = bool(
    PRIOR_ART_COMMIT == "90757a6"
    and PROTOCOL_COMMIT == "1c26fdb"
    and static_artifact["outcome"] == "REGULAR_LAPSE_IDENTITY_PROVED"
    and static_artifact["passed"] == static_artifact["tests"] == 13
    and tick1["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and tick1["passed"] == tick1["tests"] == 7
    and tick2["outcome"] == "SECOND_HOMOTHETIC_TICK_ACCEPTED"
    and tick2["passed"] == tick2["tests"] == 6
    and tick3["outcome"] == "THIRD_HOMOTHETIC_TICK_ACCEPTED"
    and tick3["passed"] == tick3["tests"] == 6
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing["passed"] == gluing["tests"] == 25
    and LAMBDA_TEXTS == ("0.5", "0.25", "0.125")
)
check(
    "all frozen artifacts, lambdas and the audited solver source pass",
    hashes_ok and provenance_ok,
    f"hashes={hashes}; base={file_digest(BASE_SOLVER_SOURCE)}",
)

accepted_states = {
    1: tuple(arb.mpf(value) for value in tick1["solutions"]["even"]["state"]),
    2: tuple(
        arb.mpf(value) for value in tick2["solutions"]["even"]["state_absolute"]
    ),
    3: tuple(
        arb.mpf(value) for value in tick3["solutions"]["even"]["state_absolute"]
    ),
}
maps = {
    parity: tuple(
        gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"]
    )
    for parity in ("even", "odd")
}
maps_ok = all(sorted(mapping) == list(range(30)) for mapping in maps.values())


def evaluate_raw(parity, lam, lower_log, upper_log, relative_rho_log):
    q_old_value = arb.exp(2*lower_log)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*upper_log)*ARB_L0_SQUARE
    rho_value = lam**2*ARB_RHO*arb.exp(relative_rho_log)
    diagonal = arb.exp(lower_log+upper_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("weak-lapse homothetic state left positive domain")
    action, gradient, branch = action_and_gradient(
        models[parity], tuple([q_old_value]*30),
        tuple([diagonal]*30+[rho_value]*5), tuple([q_new_value]*30),
    )
    local = tuple(arb.re(value) for value in gradient[30:65])
    pre = tuple(-arb.re(value) for value in gradient[:30])
    post = tuple(arb.re(value) for value in gradient[65:95])
    imaginary = maximum_imaginary(action, gradient)
    return {
        "action": action,
        "gradient": gradient,
        "local": local,
        "pre": pre,
        "post": post,
        "G": mean(local[30:]),
        "P": mean(pre),
        "branch": branch,
        "maximum_imaginary": imaginary,
        "branch_pass": branch_pass(branch, imaginary),
    }


static_controls = {}
static_ok = maps_ok
all_branch_ok = True
for lam_text in LAMBDA_TEXTS:
    lam = arb.mpf(lam_text)
    expected = lam*K0
    static_controls[lam_text] = {}
    for parity in ("even", "odd"):
        record = evaluate_raw(parity, lam, arb.mpf(0), arb.mpf(0), arb.mpf(0))
        local_error = infinity_norm(record["local"])
        pre_error = infinity_norm(tuple(value+expected for value in record["pre"]))
        post_error = infinity_norm(tuple(value-expected for value in record["post"]))
        passed_static = bool(
            record["branch_pass"] and local_error < arb.mpf("1e-60")
            and pre_error < arb.mpf("1e-60")
            and post_error < arb.mpf("1e-60")
        )
        static_ok &= passed_static
        all_branch_ok &= record["branch_pass"]
        static_controls[lam_text][parity] = {
            "local_error": local_error,
            "pre_error": pre_error,
            "post_error": post_error,
            "branch_pass": record["branch_pass"],
            "passed": passed_static,
        }
check(
    "all six scaled static states reproduce the exact all-lapse identity",
    static_ok,
)


def make_reduced_evaluator(parity, lam, lower_log, target):
    def evaluate_reduced(_parity, state):
        if _parity != parity:
            raise ValueError("parity context mismatch")
        raw = evaluate_raw(parity, lam, lower_log, state[0], state[1])
        momentum_residual = tuple(
            value-wanted for value, wanted in zip(raw["pre"], target)
        )
        reduced = (raw["G"], mean(momentum_residual))
        return {
            **raw,
            "state": tuple(state),
            "momentum_residual": momentum_residual,
            "reduced": reduced,
            "reduced_norm": infinity_norm(reduced),
        }
    return evaluate_reduced


def solve_one(parity, lam, lower_log, target, seed, label):
    global all_branch_ok
    evaluator = make_reduced_evaluator(parity, lam, lower_log, target)
    solver["evaluate_reduced"] = evaluator
    state = seed
    evaluation = evaluator(parity, state)
    all_branch_ok &= evaluation["branch_pass"]
    history = []
    jacobian_attempts = []
    all_jacobians_resolved = True
    converged = False
    failure = None
    accepted_iterations = 0
    while accepted_iterations < MAX_ITERATIONS:
        if evaluation["reduced_norm"] < RESIDUAL_TOLERANCE:
            converged = True
            break
        jacobian = calibrated_jacobian(parity, state)
        jacobian_attempts.append({"state": state, "jacobian": jacobian})
        all_branch_ok &= jacobian["branch_pass"]
        all_jacobians_resolved &= jacobian["resolved"]
        print(
            "    {} J{} smin={} eps={} resolved={}".format(
                label, len(jacobian_attempts),
                text(jacobian["singular_values"][1], 8),
                text(jacobian["epsilon"], 8), jacobian["resolved"],
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
            trial = evaluator(parity, trial_state)
            all_branch_ok &= trial["branch_pass"]
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
            "    {} iter={} alpha={} a={} r={} ||F||={}".format(
                label, accepted_iterations, text(alpha, 4), text(state[0], 14),
                text(state[1], 14), text(evaluation["reduced_norm"], 8),
            ), flush=True,
        )
    if evaluation["reduced_norm"] < RESIDUAL_TOLERANCE:
        converged = True
    if not converged and failure is None:
        failure = "ITERATION_LIMIT"
    endpoint_jacobian = None
    if converged:
        endpoint_jacobian = calibrated_jacobian(parity, state)
        all_branch_ok &= endpoint_jacobian["branch_pass"]
        all_jacobians_resolved &= endpoint_jacobian["resolved"]
        if not endpoint_jacobian["branch_pass"]:
            converged = False
            failure = "ENDPOINT_BRANCH_FAILURE"
        elif not endpoint_jacobian["resolved"]:
            converged = False
            failure = "ENDPOINT_JACOBIAN_ERROR_BAND_FAILURE"
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
        and vector_norm(momentum) < JUNCTION_TOLERANCE
        and spread(momentum) < JUNCTION_TOLERANCE
    )
    return {
        "state": state,
        "evaluation": evaluation,
        "target": target,
        "seed": seed,
        "history": history,
        "jacobian_attempts": jacobian_attempts,
        "endpoint_jacobian": endpoint_jacobian,
        "jacobians_resolved": all_jacobians_resolved,
        "converged": converged,
        "failure": failure,
        "iterations": accepted_iterations,
        "full_gate": full_gate,
    }


results = {}
all_jacobians_resolved = True
all_newton_ok = True
all_full_ok = True
for lam_text in LAMBDA_TEXTS:
    lam = arb.mpf(lam_text)
    print(f"  lambda={lam_text}", flush=True)
    results[lam_text] = {}
    for parity in ("even", "odd"):
        lower_log = arb.mpf(0)
        target = tuple([lam*K0]*30)
        step_records = []
        for step_index in (1, 2, 3):
            accepted = accepted_states[step_index]
            seed = (lam**2*accepted[0], lam**2*accepted[1])
            label = f"lambda={lam_text} {parity} n={step_index}"
            record = solve_one(parity, lam, lower_log, target, seed, label)
            step_records.append(record)
            all_jacobians_resolved &= record["jacobians_resolved"]
            all_newton_ok &= record["converged"]
            all_full_ok &= record["full_gate"]
            lower_log = record["state"][0]
            target = tuple(
                record["evaluation"]["post"][index] for index in maps[parity]
            )
        results[lam_text][parity] = step_records

parity_ok = True
parity_records = {}
for lam_text in LAMBDA_TEXTS:
    parity_records[lam_text] = []
    for step_index in range(3):
        even = results[lam_text]["even"][step_index]
        odd = results[lam_text]["odd"][step_index]
        a_difference = abs(even["state"][0]-odd["state"][0])
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
        passed_parity = bool(
            a_difference < arb.mpf("1e-25")
            and r_difference < arb.mpf("1e-25")
            and pre_difference < arb.mpf("1e-22")
            and post_difference < arb.mpf("1e-22")
        )
        parity_ok &= passed_parity
        parity_records[lam_text].append({
            "a_difference": a_difference,
            "r_difference": r_difference,
            "pre_difference": pre_difference,
            "post_difference": post_difference,
            "passed": passed_parity,
        })

targets = {
    "u2_over_u1": arb.mpf(2),
    "u3_over_u1": arb.mpf(3),
    "a2_over_u1": arb.mpf(3),
    "a3_over_u1": arb.mpf(6),
    "v2_over_v1": arb.mpf(3),
    "v3_over_v1": arb.mpf(5),
    "r2_over_v1": arb.mpf(4),
    "r3_over_v1": arb.mpf(9),
    "post1_over_k": arb.mpf(3),
    "post2_over_k": arb.mpf(5),
    "post3_over_k": arb.mpf(7),
}
observables = {}
leading_coefficients = {}
for lam_text in LAMBDA_TEXTS:
    lam = arb.mpf(lam_text)
    steps = results[lam_text]["even"]
    a = [arb.mpf(0)]+[record["state"][0] for record in steps]
    r = [arb.mpf(0)]+[record["state"][1] for record in steps]
    u = [a[index]-a[index-1] for index in (1, 2, 3)]
    v = [r[index]-r[index-1] for index in (1, 2, 3)]
    k_lam = lam*K0
    post = [mean(record["evaluation"]["post"]) for record in steps]
    observables[lam_text] = {
        "u2_over_u1": u[1]/u[0],
        "u3_over_u1": u[2]/u[0],
        "a2_over_u1": a[2]/u[0],
        "a3_over_u1": a[3]/u[0],
        "v2_over_v1": v[1]/v[0],
        "v3_over_v1": v[2]/v[0],
        "r2_over_v1": r[2]/v[0],
        "r3_over_v1": r[3]/v[0],
        "post1_over_k": post[0]/k_lam,
        "post2_over_k": post[1]/k_lam,
        "post3_over_k": post[2]/k_lam,
    }
    leading_coefficients[lam_text] = {
        "u1_over_lambda2": u[0]/lam**2,
        "v1_over_lambda2": v[0]/lam**2,
    }

asymptotic = {}
all_trend = True
all_quadratic = True
all_richardson = True
for name, target_value in targets.items():
    values = [observables[lam_text][name] for lam_text in LAMBDA_TEXTS]
    errors = [abs(value-target_value) for value in values]
    resolved = all(error > ASYMPTOTIC_FLOOR for error in errors)
    trend = bool(resolved and errors[2] < errors[1] < errors[0])
    if resolved:
        orders = (
            arb.log(errors[0]/errors[1])/arb.log(2),
            arb.log(errors[1]/errors[2])/arb.log(2),
        )
    else:
        orders = (arb.nan, arb.nan)
    quadratic = bool(
        trend and all(arb.mpf("1.8") <= order <= arb.mpf("2.2") for order in orders)
    )
    coarse = (4*values[1]-values[0])/3
    fine = (4*values[2]-values[1])/3
    richardson_epsilon = abs(fine-coarse)+ASYMPTOTIC_FLOOR
    richardson = abs(fine-target_value) <= 10*richardson_epsilon
    all_trend &= trend
    all_quadratic &= quadratic
    all_richardson &= richardson
    asymptotic[name] = {
        "target": target_value,
        "values": values,
        "errors": errors,
        "resolved": resolved,
        "trend": trend,
        "orders": orders,
        "quadratic": quadratic,
        "richardson_coarse": coarse,
        "richardson_fine": fine,
        "richardson_epsilon": richardson_epsilon,
        "richardson_integer_consistent": richardson,
    }

controls_ok = bool(hashes_ok and provenance_ok and maps_ok and static_ok)
if not controls_ok or not all_branch_ok:
    outcome = "WEAK_LAPSE_RECURRENCE_CONTROL_FAILED"
elif not all_jacobians_resolved:
    outcome = "WEAK_LAPSE_RECURRENCE_JACOBIAN_OPEN"
elif not all_newton_ok:
    outcome = "WEAK_LAPSE_RECURRENCE_NEWTON_OPEN"
elif not all_full_ok:
    outcome = "WEAK_LAPSE_RECURRENCE_FULL_GATE_FAILED"
elif not parity_ok:
    outcome = "WEAK_LAPSE_RECURRENCE_SCHEDULE_DEPENDENT"
elif all_trend and all_quadratic and all_richardson:
    outcome = "WEAK_LAPSE_QUADRATIC_INTEGER_LAW"
else:
    outcome = "WEAK_LAPSE_INTEGER_TREND_ONLY"

check("all scaled action, derivative and trial states retain the branch", all_branch_ok)
check(
    "the frozen hierarchy assigns one weak-lapse outcome",
    outcome in {
        "WEAK_LAPSE_RECURRENCE_CONTROL_FAILED",
        "WEAK_LAPSE_RECURRENCE_JACOBIAN_OPEN",
        "WEAK_LAPSE_RECURRENCE_NEWTON_OPEN",
        "WEAK_LAPSE_RECURRENCE_FULL_GATE_FAILED",
        "WEAK_LAPSE_RECURRENCE_SCHEDULE_DEPENDENT",
        "WEAK_LAPSE_INTEGER_TREND_ONLY",
        "WEAK_LAPSE_QUADRATIC_INTEGER_LAW",
    },
    outcome,
)


def serialize_solve(record):
    evaluation = record["evaluation"]
    return {
        "seed": [text(value, 60) for value in record["seed"]],
        "state": [text(value, 60) for value in record["state"]],
        "target": [text(value, 50) for value in record["target"]],
        "iterations": record["iterations"],
        "converged": record["converged"],
        "failure": record["failure"],
        "jacobians_resolved": record["jacobians_resolved"],
        "full_gate": record["full_gate"],
        "reduced_residual": [text(value, 50) for value in evaluation["reduced"]],
        "reduced_residual_norm": text(evaluation["reduced_norm"], 40),
        "diagonal_residuals": [text(value, 40) for value in evaluation["local"][:30]],
        "pole_residuals": [text(value, 40) for value in evaluation["local"][30:]],
        "pre_momentum": [text(value, 50) for value in evaluation["pre"]],
        "post_momentum": [text(value, 50) for value in evaluation["post"]],
        "junction_residual": [text(value, 40) for value in evaluation["momentum_residual"]],
        "junction_norm": text(vector_norm(evaluation["momentum_residual"]), 40),
        "endpoint_jacobian": serialize_jacobian(record["endpoint_jacobian"]),
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
            for item in record["history"]
        ],
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": {
        **hashes,
        "base_solver_functions": file_digest(BASE_SOLVER_SOURCE),
    },
    "lambdas": list(LAMBDA_TEXTS),
    "fixed_mass": True,
    "tick4_target_parsed": False,
    "static_controls": {
        lam_text: {
            parity: {
                key: (value if isinstance(value, bool) else text(value, 40))
                for key, value in record.items()
            }
            for parity, record in parity_records_for_lambda.items()
        }
        for lam_text, parity_records_for_lambda in static_controls.items()
    },
    "solves": {
        lam_text: {
            parity: [serialize_solve(record) for record in records]
            for parity, records in parity_records_for_lambda.items()
        }
        for lam_text, parity_records_for_lambda in results.items()
    },
    "parity": {
        lam_text: [
            {
                key: (value if isinstance(value, bool) else text(value, 40))
                for key, value in record.items()
            }
            for record in records
        ]
        for lam_text, records in parity_records.items()
    },
    "observables": {
        lam_text: {name: text(value, 60) for name, value in values.items()}
        for lam_text, values in observables.items()
    },
    "leading_coefficients": {
        lam_text: {name: text(value, 60) for name, value in values.items()}
        for lam_text, values in leading_coefficients.items()
    },
    "asymptotic": {
        name: {
            "target": text(record["target"], 20),
            "values": [text(value, 60) for value in record["values"]],
            "errors": [text(value, 40) for value in record["errors"]],
            "resolved": record["resolved"],
            "trend": record["trend"],
            "orders": [text(value, 30) for value in record["orders"]],
            "quadratic": record["quadratic"],
            "richardson_coarse": text(record["richardson_coarse"], 50),
            "richardson_fine": text(record["richardson_fine"], 50),
            "richardson_epsilon": text(record["richardson_epsilon"], 35),
            "richardson_integer_consistent": record["richardson_integer_consistent"],
        }
        for name, record in asymptotic.items()
    },
    "classification": {
        "all_trends": all_trend,
        "all_quadratic_orders": all_quadratic,
        "all_richardson_integer_consistent": all_richardson,
        "leading_integer_law": outcome == "WEAK_LAPSE_QUADRATIC_INTEGER_LAW",
        "exact_all_order_recurrence": False,
        "spatial_refinement": "NOT TESTED",
        "emergent_time": "OPEN",
    },
    "outcome": outcome,
    "tests": solver["tests"],
    "passed": solver["passed"],
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

for name, record in asymptotic.items():
    print(
        "  {} target={} values=({}, {}, {}) orders=({}, {}) rich={} quad={}".format(
            name, text(record["target"], 8),
            *(text(value, 13) for value in record["values"]),
            *(text(value, 8) for value in record["orders"]),
            record["richardson_integer_consistent"], record["quadratic"],
        ), flush=True,
    )
print(f"OUTCOME: {outcome}", flush=True)
print(f"Tests passed: {solver['passed']}/{solver['tests']}", flush=True)
raise SystemExit(0 if solver["passed"] == solver["tests"] else 1)

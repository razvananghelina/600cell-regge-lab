#!/usr/bin/env python3
"""Out-of-sample fourth weak-lapse tick of the canonical 600-cell dust map.

Prior-art commit: 40d77c7.
Protocol commit: 9eae63b.
"""

import ast
import contextlib
import hashlib
import io
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
WEAK_SOURCE = HERE / "verify_gravity_600cell_dust_weak_lapse_recurrence.py"
BASE_SOLVER_SOURCE = HERE / "verify_gravity_600cell_dust_second_tick_local_correction.py"
RESPONSE_SOURCE = HERE / "verify_gravity_600cell_dust_homothetic_mass_conservation.py"
WEAK_ARTIFACT = HERE / "gravity_600cell_dust_weak_lapse_recurrence.json"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_fourth_tick.json"

PRIOR_ART_COMMIT = "40d77c7"
PROTOCOL_COMMIT = "9eae63b"
WEAK_SOURCE_SHA256 = "fe2298b6adcc934f04a995ce9d37a30b99b84f6064f1baeb6e544d002d2e97e5"
BASE_SOLVER_SHA256 = "cef59fa0bc3a1c8fa3be0193234371b7dda303a0ec72683ddcdd88bcb40f3725"
WEAK_ARTIFACT_SHA256 = "500be1c4e2d7ec4104b9773bc1cfc71065c9d930607119eb616367d18fa5d8f9"
GLUING_ARTIFACT_SHA256 = "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77"

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
REPRODUCTION_TOLERANCE = arb.mpf("1e-45")
TRAINING_BAND = arb.mpf("4.6222921056804246599831556548181231e-10")
MAX_ITERATIONS = 8
MAX_DAMPING = 10


def file_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_weak_functions():
    """Load only audited function definitions; do not run the 3-tick script."""
    wanted = {
        "load_solver_functions",
        "evaluate_raw",
        "make_reduced_evaluator",
        "solve_one",
        "serialize_solve",
    }
    tree = ast.parse(WEAK_SOURCE.read_text(), filename=str(WEAK_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    if {node.name for node in body} != wanted:
        raise RuntimeError("audited weak-lapse function set is incomplete")
    namespace = {
        "__file__": str(WEAK_SOURCE),
        "__name__": "fourth_tick_imported_weak_functions",
        "arb": arb,
        "ast": ast,
        "contextlib": contextlib,
        "io": io,
        "json": json,
        "Path": Path,
        "BASE_SOLVER_SOURCE": BASE_SOLVER_SOURCE,
        "RESPONSE_SOURCE": RESPONSE_SOURCE,
        "STEP_SETS": STEP_SETS,
        "ARITHMETIC_FLOOR": ARITHMETIC_FLOOR,
        "ENTRY_FACTOR": ENTRY_FACTOR,
        "NONZERO_FACTOR": NONZERO_FACTOR,
        "RESIDUAL_TOLERANCE": RESIDUAL_TOLERANCE,
        "JUNCTION_TOLERANCE": JUNCTION_TOLERANCE,
        "MAX_ITERATIONS": MAX_ITERATIONS,
        "MAX_DAMPING": MAX_DAMPING,
    }
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(WEAK_SOURCE), "exec"),
        namespace,
    )
    return namespace


weak = load_weak_functions()
solver = weak["load_solver_functions"]()
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

print("Out-of-sample fourth weak-lapse dust tick", flush=True)
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

weak.update({
    "solver": solver,
    "check": check,
    "text": text,
    "mean": mean,
    "vector_norm": vector_norm,
    "infinity_norm": infinity_norm,
    "spread": spread,
    "maximum_imaginary": maximum_imaginary,
    "calibrated_jacobian": calibrated_jacobian,
    "serialize_jacobian": serialize_jacobian,
    "models": models,
    "action_and_gradient": action_and_gradient,
    "branch_pass": branch_pass,
    "ARB_L0": ARB_L0,
    "ARB_L0_SQUARE": ARB_L0_SQUARE,
    "ARB_RHO": ARB_RHO,
    "ARB_TAU": ARB_TAU,
    "ARB_EPSILON_3": ARB_EPSILON_3,
    "all_branch_ok": True,
})
evaluate_raw = weak["evaluate_raw"]
solve_one = weak["solve_one"]
serialize_solve = weak["serialize_solve"]

check(
    "the imported complete-action evaluator retains its six controls",
    response_core["tests"] == response_core["passed"] == 6,
)

weak_artifact = json.loads(WEAK_ARTIFACT.read_text())
gluing = json.loads(GLUING_ARTIFACT.read_text())
hashes = {
    "weak_source": file_digest(WEAK_SOURCE),
    "base_solver": file_digest(BASE_SOLVER_SOURCE),
    "weak_artifact": file_digest(WEAK_ARTIFACT),
    "gluing": file_digest(GLUING_ARTIFACT),
}
expected_hashes = {
    "weak_source": WEAK_SOURCE_SHA256,
    "base_solver": BASE_SOLVER_SHA256,
    "weak_artifact": WEAK_ARTIFACT_SHA256,
    "gluing": GLUING_ARTIFACT_SHA256,
}
training_band_from_artifact = 10*max(
    arb.mpf(record["richardson_epsilon"])
    for record in weak_artifact["asymptotic"].values()
)
maps = {
    parity: tuple(
        gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"]
    )
    for parity in ("even", "odd")
}
provenance_ok = bool(
    hashes == expected_hashes
    and PRIOR_ART_COMMIT == "40d77c7"
    and PROTOCOL_COMMIT == "9eae63b"
    and weak_artifact["outcome"] == "WEAK_LAPSE_QUADRATIC_INTEGER_LAW"
    and weak_artifact["tests"] == weak_artifact["passed"] == 5
    and weak_artifact["tick4_target_parsed"] is False
    and tuple(weak_artifact["lambdas"]) == LAMBDA_TEXTS
    and training_band_from_artifact == TRAINING_BAND
    and all(sorted(mapping) == list(range(30)) for mapping in maps.values())
)
check(
    "all frozen inputs and the precommitted training band pass",
    provenance_ok,
    f"hashes={hashes}; B_train={text(training_band_from_artifact, 40)}",
)

reproduction = {}
reproduction_ok = True
all_branch_ok = True
for lam_text in LAMBDA_TEXTS:
    lam = arb.mpf(lam_text)
    reproduction[lam_text] = {}
    for parity in ("even", "odd"):
        records = weak_artifact["solves"][lam_text][parity]
        second_state = tuple(arb.mpf(value) for value in records[1]["state"])
        third_state = tuple(arb.mpf(value) for value in records[2]["state"])
        raw = evaluate_raw(
            parity, lam, second_state[0], third_state[0], third_state[1]
        )
        stored_pre = tuple(
            arb.mpf(value) for value in records[2]["pre_momentum"]
        )
        stored_post = tuple(
            arb.mpf(value) for value in records[2]["post_momentum"]
        )
        pre_error = infinity_norm(tuple(
            value-stored for value, stored in zip(raw["pre"], stored_pre)
        ))
        post_error = infinity_norm(tuple(
            value-stored for value, stored in zip(raw["post"], stored_post)
        ))
        passed_reproduction = bool(
            raw["branch_pass"]
            and pre_error < REPRODUCTION_TOLERANCE
            and post_error < REPRODUCTION_TOLERANCE
        )
        reproduction_ok &= passed_reproduction
        all_branch_ok &= raw["branch_pass"]
        reproduction[lam_text][parity] = {
            "pre_error": pre_error,
            "post_error": post_error,
            "branch_pass": raw["branch_pass"],
            "passed": passed_reproduction,
        }
check(
    "all committed third endpoints reproduce before extrapolation",
    reproduction_ok,
)

results = {}
all_jacobians_resolved = True
all_newton_ok = True
all_full_ok = True
for lam_text in LAMBDA_TEXTS:
    lam = arb.mpf(lam_text)
    print(f"  lambda={lam_text}", flush=True)
    results[lam_text] = {}
    for parity in ("even", "odd"):
        records = weak_artifact["solves"][lam_text][parity]
        first_state = tuple(arb.mpf(value) for value in records[0]["state"])
        third_state = tuple(arb.mpf(value) for value in records[2]["state"])
        stored_post = tuple(
            arb.mpf(value) for value in records[2]["post_momentum"]
        )
        target = tuple(stored_post[index] for index in maps[parity])
        seed = (10*first_state[0], 16*first_state[1])
        label = f"lambda={lam_text} {parity} n=4"
        record = solve_one(
            parity, lam, third_state[0], target, seed, label
        )
        results[lam_text][parity] = record
        all_jacobians_resolved &= record["jacobians_resolved"]
        all_newton_ok &= record["converged"]
        all_full_ok &= record["full_gate"]

all_branch_ok &= weak["all_branch_ok"]
parity_ok = True
parity_records = {}
for lam_text in LAMBDA_TEXTS:
    even = results[lam_text]["even"]
    odd = results[lam_text]["odd"]
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
    parity_records[lam_text] = {
        "a_difference": a_difference,
        "r_difference": r_difference,
        "pre_difference": pre_difference,
        "post_difference": post_difference,
        "passed": passed_parity,
    }

targets = {
    "u4_over_u1": arb.mpf(4),
    "a4_over_u1": arb.mpf(10),
    "v4_over_v1": arb.mpf(7),
    "r4_over_v1": arb.mpf(16),
    "post4_over_k": arb.mpf(9),
}
observables = {parity: {} for parity in ("even", "odd")}
for parity in ("even", "odd"):
    for lam_text in LAMBDA_TEXTS:
        lam = arb.mpf(lam_text)
        prior = weak_artifact["solves"][lam_text][parity]
        first_state = tuple(arb.mpf(value) for value in prior[0]["state"])
        third_state = tuple(arb.mpf(value) for value in prior[2]["state"])
        fourth = results[lam_text][parity]
        u1 = first_state[0]
        v1 = first_state[1]
        u4 = fourth["state"][0]-third_state[0]
        v4 = fourth["state"][1]-third_state[1]
        post4 = mean(fourth["evaluation"]["post"])
        observables[parity][lam_text] = {
            "u4_over_u1": u4/u1,
            "a4_over_u1": fourth["state"][0]/u1,
            "v4_over_v1": v4/v1,
            "r4_over_v1": fourth["state"][1]/v1,
            "post4_over_k": post4/(lam*K0),
        }

asymptotic = {parity: {} for parity in ("even", "odd")}
all_trend = True
all_quadratic = True
all_internal_richardson = True
all_external_band = True
for parity in ("even", "odd"):
    for name, target_value in targets.items():
        values = [
            observables[parity][lam_text][name] for lam_text in LAMBDA_TEXTS
        ]
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
            trend
            and all(
                arb.mpf("1.8") <= order <= arb.mpf("2.2")
                for order in orders
            )
        )
        coarse = (4*values[1]-values[0])/3
        fine = (4*values[2]-values[1])/3
        epsilon4 = abs(fine-coarse)+ASYMPTOTIC_FLOOR
        internal = abs(fine-target_value) <= 10*epsilon4
        external = abs(fine-target_value) <= TRAINING_BAND
        all_trend &= trend
        all_quadratic &= quadratic
        all_internal_richardson &= internal
        all_external_band &= external
        asymptotic[parity][name] = {
            "target": target_value,
            "values": values,
            "errors": errors,
            "resolved": resolved,
            "trend": trend,
            "orders": orders,
            "quadratic": quadratic,
            "richardson_coarse": coarse,
            "richardson_fine": fine,
            "epsilon4": epsilon4,
            "internal_richardson": internal,
            "external_training_band": external,
        }

controls_ok = bool(provenance_ok and reproduction_ok)
if not controls_ok or not all_branch_ok:
    outcome = "FOURTH_TICK_CONTROL_FAILED"
elif not all_jacobians_resolved:
    outcome = "FOURTH_TICK_JACOBIAN_OPEN"
elif not all_newton_ok:
    outcome = "FOURTH_TICK_NEWTON_OPEN"
elif not all_full_ok:
    outcome = "FOURTH_TICK_FULL_GATE_FAILED"
elif not parity_ok:
    outcome = "FOURTH_TICK_SCHEDULE_DEPENDENT"
elif all_trend and not (
    all_quadratic and all_internal_richardson and all_external_band
):
    outcome = "FOURTH_TICK_INTEGER_TREND_ONLY"
elif not (
    all_trend and all_quadratic
    and all_internal_richardson and all_external_band
):
    outcome = "FOURTH_TICK_WEAK_LAPSE_PREDICTION_REFUTED"
else:
    outcome = "FOURTH_TICK_WEAK_LAPSE_PREDICTION_CONFIRMED"

check(
    "all fourth-step action, derivative and trial states retain the branch",
    all_branch_ok,
)
allowed_outcomes = {
    "FOURTH_TICK_CONTROL_FAILED",
    "FOURTH_TICK_JACOBIAN_OPEN",
    "FOURTH_TICK_NEWTON_OPEN",
    "FOURTH_TICK_FULL_GATE_FAILED",
    "FOURTH_TICK_SCHEDULE_DEPENDENT",
    "FOURTH_TICK_INTEGER_TREND_ONLY",
    "FOURTH_TICK_WEAK_LAPSE_PREDICTION_REFUTED",
    "FOURTH_TICK_WEAK_LAPSE_PREDICTION_CONFIRMED",
}
check(
    "the frozen hierarchy assigns one fourth-tick outcome",
    outcome in allowed_outcomes,
    outcome,
)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "lambdas": list(LAMBDA_TEXTS),
    "fixed_mass": True,
    "lambda1_fourth_state_evaluated": False,
    "training_band": text(TRAINING_BAND, 45),
    "training_band_from_artifact": text(training_band_from_artifact, 45),
    "reproduction": {
        lam_text: {
            parity: {
                key: (value if isinstance(value, bool) else text(value, 40))
                for key, value in record.items()
            }
            for parity, record in parity_records_for_lambda.items()
        }
        for lam_text, parity_records_for_lambda in reproduction.items()
    },
    "solves": {
        lam_text: {
            parity: serialize_solve(record)
            for parity, record in parity_records_for_lambda.items()
        }
        for lam_text, parity_records_for_lambda in results.items()
    },
    "parity": {
        lam_text: {
            key: (value if isinstance(value, bool) else text(value, 40))
            for key, value in record.items()
        }
        for lam_text, record in parity_records.items()
    },
    "observables": {
        parity: {
            lam_text: {name: text(value, 60) for name, value in values.items()}
            for lam_text, values in records.items()
        }
        for parity, records in observables.items()
    },
    "asymptotic": {
        parity: {
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
                "epsilon4": text(record["epsilon4"], 40),
                "internal_richardson": record["internal_richardson"],
                "external_training_band": record["external_training_band"],
            }
            for name, record in records.items()
        }
        for parity, records in asymptotic.items()
    },
    "classification": {
        "all_trends": all_trend,
        "all_quadratic_orders": all_quadratic,
        "all_internal_richardson": all_internal_richardson,
        "all_external_training_band": all_external_band,
        "out_of_sample_iteration_index": True,
        "exact_all_order_recurrence": False,
        "spatial_refinement": "NOT TESTED",
        "emergent_time": "OPEN",
    },
    "outcome": outcome,
    "tests": solver["tests"],
    "passed": solver["passed"],
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

for parity in ("even", "odd"):
    for name, record in asymptotic[parity].items():
        print(
            "  {} {} target={} values=({}, {}, {}) orders=({}, {}) "
            "internal={} external={} quadratic={}".format(
                parity, name, text(record["target"], 8),
                *(text(value, 13) for value in record["values"]),
                *(text(value, 8) for value in record["orders"]),
                record["internal_richardson"],
                record["external_training_band"],
                record["quadratic"],
            ),
            flush=True,
        )
print(f"OUTCOME: {outcome}", flush=True)
print(f"Tests passed: {solver['passed']}/{solver['tests']}", flush=True)
raise SystemExit(0 if solver["passed"] == solver["tests"] else 1)

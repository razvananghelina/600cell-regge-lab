#!/usr/bin/env python3
"""Solve and compare the 32 preregistered nonlinear canonical input pairs.

Prior-art commit: 526a202.
Protocol commit: 05f76c3.
Frozen-seed commit: b6370bd.
"""

import ast
from collections import Counter
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
import json
from pathlib import Path
import sys

import mpmath as arb


HERE = Path(__file__).resolve().parent
SEED_INPUT = HERE / "gravity_600cell_dust_nonlinear_boundary_covariance_seeds.json"
CANONICAL_SOURCE = HERE / "verify_gravity_600cell_dust_canonical_legendre_rank.py"
ACTION_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
OUTPUT = HERE / "gravity_600cell_dust_nonlinear_boundary_covariance.json"

PRIOR_ART_COMMIT = "526a202"
PROTOCOL_COMMIT = "05f76c3"
SEED_COMMIT = "b6370bd"
SEED_SHA256 = "2104c69ba6b21d3a3d92c7071d7f2702cb7d33f7f0e3ff17954f64c469f0c01d"
CANONICAL_SOURCE_SHA256 = "396c491fe51a9f5e04fa8402e2e5b16884fe23fc5057d8ded325e6064fbd3b9e"
ACTION_SOURCE_SHA256 = "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf"
DPS = 100
arb.mp.dps = DPS
RESIDUAL_TOLERANCE = arb.mpf("1e-55")
CORRECTION_TOLERANCE = arb.mpf("1e-45")
IMAGINARY_TOLERANCE = arb.mpf("1e-70")
ANGLE_TOLERANCE = arb.mpf("1e-6")
MAX_ITERATIONS = 20
MAX_BACKTRACKING = 12
CONSISTENT_FACTOR = arb.mpf(10)
BROKEN_FACTOR = arb.mpf(100)
UNCERTAINTY_FLOOR = arb.mpf("1e-70")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_audited_functions():
    wanted = {
        "orbit_sort_key",
        "augment_boundary_orbits",
        "arb_log_minus",
        "arb_signed_volume_square",
        "arb_angle_record",
        "triangle_area_square",
        "triangle_area_square_partials",
        "edge_data",
        "simplex_squared",
        "action_and_gradient",
    }
    tree = ast.parse(CANONICAL_SOURCE.read_text(), filename=str(CANONICAL_SOURCE))
    body = [node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in wanted]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch: missing={wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]),
                str(CANONICAL_SOURCE), "exec"),
        globals(),
    )


def vector_norm(vector):
    return arb.sqrt(sum(abs(value)**2 for value in vector))


def infinity_norm(vector):
    return max(abs(value) for value in vector)


def matrix_from_rows(rows):
    return arb.matrix([[arb.mpf(value) for value in row] for row in rows])


def permute_output(permutation, vector):
    result = [arb.mpf(0) for _ in range(60)]
    for source, target in enumerate(permutation):
        result[target] = vector[source]
        result[30+target] = vector[30+source]
    return tuple(result)


def number(value, digits=70):
    return arb.nstr(value, digits, strip_zeros=False)


seed = json.loads(SEED_INPUT.read_text())
hashes = {
    "seeds": digest(SEED_INPUT),
    "canonical_source": digest(CANONICAL_SOURCE),
    "action_source": digest(ACTION_SOURCE),
}
provenance_ok = bool(
    hashes == {
        "seeds": SEED_SHA256,
        "canonical_source": CANONICAL_SOURCE_SHA256,
        "action_source": ACTION_SOURCE_SHA256,
    }
    and seed.get("prior_art_commit") == PRIOR_ART_COMMIT
    and seed.get("protocol_commit") == PROTOCOL_COMMIT
    and seed.get("precision_correction_commit") == "8981046"
    and seed.get("outcome") == "NONLINEAR_BOUNDARY_COVARIANCE_CASES_FROZEN"
    and seed.get("passed") == seed.get("tests") == 10
    and seed.get("number_of_paired_cases") == 32
    and seed.get("nonlinear_perturbed_action_evaluations") == 0
    and seed.get("nonlinear_outputs_compared") is False
    and seed.get("continuum_target_parsed") is False
    and seed.get("speed_target_parsed") is False
    and seed.get("full_720_edge_carrier") is False
)

load_audited_functions()
spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_nonlinear_covariance", ACTION_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
models = {
    parity: augment_boundary_orbits(model)
    for parity, model in gro.models.items()
}
carrier_ok = bool(
    gro.tests == gro.passed == 43
    and all(
        len(model["old_orbits"]) == 30
        and len(model["edge_orbits"]) == 35
        and len(model["final_orbits"]) == 30
        for model in models.values()
    )
)

ARB_I = arb.mpc(0, 1)
ARB_M_STAR = arb.mpf(10)
ARB_ZETA = (arb.pi**2*arb.sqrt(2)/50)**(arb.mpf(1)/3)
ARB_R0 = 4*ARB_M_STAR/(3*arb.pi)
ARB_L0 = ARB_ZETA*ARB_R0
ARB_L0_SQUARE = ARB_L0**2
ARB_EPSILON_3 = 2*arb.pi-5*arb.acos(arb.mpf(1)/3)
ARB_MASS = (90/arb.pi)*ARB_EPSILON_3*ARB_L0

p_star = arb.mpf(seed["p_star"])
physical_map = tuple(seed["physical_edge_permutation"])
base = {}
jacobians = {}
for parity in ("even", "odd"):
    source = seed["parities"][parity]
    base[parity] = {
        "old": tuple(arb.mpf(value) for value in source["base_old"]),
        "x": tuple(arb.mpf(value) for value in source["base_x"]),
        "new": tuple(arb.mpf(value) for value in source["base_new"]),
        "pre": tuple(arb.mpf(value) for value in source["base_pre"]),
    }
    jacobians[parity] = {
        "operational": matrix_from_rows(source["canonical_operational"]),
        "validation": matrix_from_rows(source["canonical_validation"]),
    }


def evaluate(parity, old_logs, target, state):
    q_old = tuple(arb.exp(value) for value in old_logs)
    x = tuple(arb.exp(state[index]) for index in range(35))
    q_new = tuple(arb.exp(state[35+index]) for index in range(30))
    action, gradient, branch = action_and_gradient(
        models[parity], q_old, x, q_new
    )
    residual = tuple(arb.re(gradient[30+index]) for index in range(35)) + tuple(
        -arb.re(gradient[index])-target[index] for index in range(30)
    )
    post = tuple(arb.re(gradient[65+index]) for index in range(30))
    output = tuple(state[35+index] for index in range(30)) + tuple(
        value/p_star for value in post
    )
    maximum_imaginary = max(
        abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient)
    )
    branch_ok = bool(
        branch["negative_counts"] == Counter({1: 2400})
        and branch["minimum_leading_minor"] > 0
        and branch["minimum_argument"] > ANGLE_TOLERANCE
        and maximum_imaginary < IMAGINARY_TOLERANCE
    )
    return {
        "residual": residual,
        "residual_norm": infinity_norm(residual),
        "output": output,
        "branch": branch,
        "maximum_imaginary": maximum_imaginary,
        "branch_ok": branch_ok,
    }


def solve(parity, case, variant):
    factor = arb.mpf(case["factor"])
    record = case["parities"][parity]
    ray = tuple(arb.mpf(value) for value in record["input_ray"])
    old_logs = tuple(
        arb.log(base[parity]["old"][index])+factor*ray[index]
        for index in range(30)
    )
    target = tuple(
        base[parity]["pre"][index]+factor*ray[30+index]
        for index in range(30)
    )
    seed_delta = tuple(
        arb.mpf(value) for value in record[f"unknown_seed_delta_{variant}"]
    )
    state = arb.matrix(
        [arb.log(value) for value in base[parity]["x"]]
        + [arb.log(value) for value in base[parity]["new"]]
    )
    state += arb.matrix(seed_delta)
    jacobian = jacobians[parity][variant]
    trace = []
    branch_all = True
    evaluation_error = None
    for iteration in range(MAX_ITERATIONS+1):
        try:
            current = evaluate(parity, old_logs, target, state)
        except Exception as error:
            evaluation_error = repr(error)
            branch_all = False
            break
        branch_all &= current["branch_ok"]
        if not current["branch_ok"]:
            return {
                "success": False,
                "reason": "BRANCH_FAILURE",
                "branch_all": False,
                "trace": trace,
                "evaluation_error": evaluation_error,
            }
        residual_norm = current["residual_norm"]
        try:
            delta = arb.lu_solve(
                jacobian,
                arb.matrix([-value for value in current["residual"]]),
            )
        except (ZeroDivisionError, ValueError) as error:
            return {
                "success": False,
                "reason": "JACOBIAN_SOLVE_FAILED",
                "branch_all": branch_all,
                "trace": trace,
                "evaluation_error": repr(error),
            }
        correction_norm = max(abs(delta[index]) for index in range(65))
        if (residual_norm < RESIDUAL_TOLERANCE
                and correction_norm < CORRECTION_TOLERANCE):
            corrected_state = state+delta
            try:
                corrected = evaluate(parity, old_logs, target, corrected_state)
            except Exception as error:
                return {
                    "success": False,
                    "reason": "CORRECTION_EVALUATION_FAILED",
                    "branch_all": False,
                    "trace": trace,
                    "evaluation_error": repr(error),
                }
            branch_all &= corrected["branch_ok"]
            correction_output = vector_norm(tuple(
                corrected["output"][index]-current["output"][index]
                for index in range(60)
            ))
            success = bool(
                corrected["branch_ok"]
                and current["maximum_imaginary"] < IMAGINARY_TOLERANCE
            )
            return {
                "success": success,
                "reason": "CONVERGED" if success else "FINAL_CONTROL_FAILED",
                "branch_all": branch_all,
                "iterations": iteration,
                "residual_norm": residual_norm,
                "correction_norm": correction_norm,
                "correction_output": correction_output,
                "output": current["output"],
                "minimum_leading_minor": min(
                    current["branch"]["minimum_leading_minor"],
                    corrected["branch"]["minimum_leading_minor"],
                ),
                "minimum_argument": min(
                    current["branch"]["minimum_argument"],
                    corrected["branch"]["minimum_argument"],
                ),
                "maximum_imaginary": max(
                    current["maximum_imaginary"], corrected["maximum_imaginary"]
                ),
                "trace": trace,
                "evaluation_error": None,
            }
        if iteration >= MAX_ITERATIONS:
            return {
                "success": False,
                "reason": "ITERATION_LIMIT",
                "branch_all": branch_all,
                "residual_norm": residual_norm,
                "correction_norm": correction_norm,
                "trace": trace,
                "evaluation_error": None,
            }
        accepted = False
        accepted_alpha = None
        valid_trials = 0
        for power in range(MAX_BACKTRACKING+1):
            alpha = arb.mpf(2)**(-power)
            trial_state = state+alpha*delta
            try:
                trial = evaluate(parity, old_logs, target, trial_state)
            except Exception:
                branch_all = False
                continue
            branch_all &= trial["branch_ok"]
            if not trial["branch_ok"]:
                continue
            valid_trials += 1
            if trial["residual_norm"] <= (1-alpha/4)*residual_norm:
                state = trial_state
                accepted = True
                accepted_alpha = alpha
                break
        trace.append({
            "iteration": iteration,
            "residual_norm": number(residual_norm, 35),
            "correction_norm": number(correction_norm, 35),
            "accepted_alpha": (
                number(accepted_alpha, 20) if accepted_alpha is not None else None
            ),
        })
        if not accepted:
            return {
                "success": False,
                "reason": "BRANCH_FAILURE" if valid_trials == 0 else "NO_ARMIJO_STEP",
                "branch_all": branch_all,
                "residual_norm": residual_norm,
                "correction_norm": correction_norm,
                "trace": trace,
                "evaluation_error": None,
            }
    return {
        "success": False,
        "reason": "EVALUATION_EXCEPTION",
        "branch_all": branch_all,
        "trace": trace,
        "evaluation_error": evaluation_error,
    }


def compact_solve(result):
    payload = {
        key: value for key, value in result.items()
        if key not in {
            "output", "residual_norm", "correction_norm", "correction_output",
            "minimum_leading_minor", "minimum_argument", "maximum_imaginary",
        }
    }
    for key in (
        "residual_norm", "correction_norm", "correction_output",
        "minimum_leading_minor", "minimum_argument", "maximum_imaginary",
    ):
        if key in result:
            payload[key] = number(result[key])
    if "output" in result:
        payload["output"] = [number(value) for value in result["output"]]
    return payload


case_results = []
all_branch_controls = True
for case_index, case in enumerate(seed["cases"], 1):
    print(f"[{case_index:02d}/32] {case['id']}", flush=True)
    solves = {
        parity: {
            variant: solve(parity, case, variant)
            for variant in ("operational", "validation")
        }
        for parity in ("even", "odd")
    }
    all_branch_controls &= all(
        result["branch_all"]
        for parity in solves.values() for result in parity.values()
    )
    all_success = all(
        result["success"]
        for parity in solves.values() for result in parity.values()
    )
    if not all_success:
        classification = "OPEN_SOLVE"
        defect = uncertainty = ratio = None
    else:
        operational_even = solves["even"]["operational"]["output"]
        operational_odd = solves["odd"]["operational"]["output"]
        mapped_even = permute_output(physical_map, operational_even)
        defect = vector_norm(tuple(
            operational_odd[index]-mapped_even[index] for index in range(60)
        ))
        uncertainty = UNCERTAINTY_FLOOR
        for parity in ("even", "odd"):
            op = solves[parity]["operational"]
            val = solves[parity]["validation"]
            uncertainty += vector_norm(tuple(
                op["output"][index]-val["output"][index] for index in range(60)
            ))
            uncertainty += max(op["correction_output"], val["correction_output"])
        ratio = defect/uncertainty
        if defect <= CONSISTENT_FACTOR*uncertainty:
            classification = "COVARIANT"
        elif defect > BROKEN_FACTOR*uncertainty:
            classification = "BROKEN"
        else:
            classification = "OPEN"
    case_results.append({
        "id": case["id"],
        "direction_index": case["direction_index"],
        "sector": case["sector"],
        "sign": case["sign"],
        "level": case["level"],
        "amplitude": case["amplitude"],
        "classification": classification,
        "defect": number(defect) if defect is not None else None,
        "uncertainty": number(uncertainty) if uncertainty is not None else None,
        "ratio": number(ratio) if ratio is not None else None,
        "solves": {
            parity: {variant: compact_solve(result)
                     for variant, result in variants.items()}
            for parity, variants in solves.items()
        },
    })

classification_counts = Counter(record["classification"] for record in case_results)
order_diagnostics = []
for direction in range(1, 5):
    for sector in ("POSITION", "MOMENTUM"):
        for sign in (-1, 1):
            pair = [record for record in case_results
                    if record["direction_index"] == direction
                    and record["sector"] == sector
                    and record["sign"] == sign]
            half = next(record for record in pair if record["level"] == "0.5")
            full = next(record for record in pair if record["level"] == "1.0")
            diagnostic = {
                "direction_index": direction,
                "sector": sector,
                "sign": sign,
                "available": False,
                "label": "NOT_AVAILABLE",
            }
            if half["classification"] == full["classification"] == "BROKEN":
                d_half = arb.mpf(half["defect"])
                d_full = arb.mpf(full["defect"])
                u_half = arb.mpf(half["uncertainty"])
                u_full = arb.mpf(full["uncertainty"])
                observed = arb.log(d_full/d_half, 2)
                interval_ok = d_half > u_half and d_full > u_full
                if interval_ok:
                    lower = arb.log((d_full-u_full)/(d_half+u_half), 2)
                    upper = arb.log((d_full+u_full)/(d_half-u_half), 2)
                    label = (
                        "QUADRATIC_COMPATIBLE"
                        if lower >= arb.mpf("1.5") and upper <= arb.mpf("2.5")
                        else "NONQUADRATIC_RESOLVED"
                    )
                else:
                    lower = upper = None
                    label = "ORDER_OPEN"
                diagnostic.update({
                    "available": True,
                    "observed_order": number(observed),
                    "lower_order": number(lower) if lower is not None else None,
                    "upper_order": number(upper) if upper is not None else None,
                    "label": label,
                })
            order_diagnostics.append(diagnostic)

cases_ok = bool(
    len(case_results) == 32
    and len({record["id"] for record in case_results}) == 32
    and sum(classification_counts.values()) == 32
)
implementation_ok = all(
    not result.get("evaluation_error")
    for record in case_results
    for parity in record["solves"].values()
    for result in parity.values()
)
classification_ok = set(classification_counts) <= {
    "COVARIANT", "BROKEN", "OPEN", "OPEN_SOLVE"
}

control_ok = bool(
    provenance_ok and carrier_ok and cases_ok and implementation_ok
    and all_branch_controls and classification_ok
)
if not control_ok:
    outcome = "NONLINEAR_BOUNDARY_COVARIANCE_CONTROL_FAILED"
elif classification_counts["BROKEN"]:
    outcome = "NONLINEAR_BOUNDARY_COVARIANCE_BROKEN_ON_FROZEN_CASES"
elif classification_counts["OPEN"] or classification_counts["OPEN_SOLVE"]:
    outcome = "NONLINEAR_BOUNDARY_COVARIANCE_OPEN"
else:
    outcome = "NONLINEAR_BOUNDARY_COVARIANCE_CONSISTENT_ON_FROZEN_CASES"

tests = [
    ("frozen seed and audited source hashes", provenance_ok),
    ("imported carrier retains 43/43 and quotient dimensions", carrier_ok),
    ("all 32 frozen case identifiers evaluated exactly once", cases_ok),
    ("no action evaluation raised an implementation exception", implementation_ok),
    ("all evaluated geometries retain the Lorentzian branch", all_branch_controls),
    ("every paired case receives exactly one frozen label", classification_ok),
    ("no continuum, speed, or full-carrier target parsed", True),
    ("mechanical outcome assigned", outcome in {
        "NONLINEAR_BOUNDARY_COVARIANCE_CONTROL_FAILED",
        "NONLINEAR_BOUNDARY_COVARIANCE_BROKEN_ON_FROZEN_CASES",
        "NONLINEAR_BOUNDARY_COVARIANCE_OPEN",
        "NONLINEAR_BOUNDARY_COVARIANCE_CONSISTENT_ON_FROZEN_CASES",
    }),
]
passed = sum(bool(ok) for _, ok in tests)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "seed_commit": SEED_COMMIT,
    "input_sha256": hashes,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "full_720_edge_carrier": False,
    "number_of_paired_cases": 32,
    "classification_counts": dict(classification_counts),
    "order_diagnostics": order_diagnostics,
    "cases": case_results,
    "passed": passed,
    "tests": len(tests),
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"classifications={dict(classification_counts)}")
if any(record["defect"] is not None for record in case_results):
    available = [record for record in case_results if record["defect"] is not None]
    print(
        "defect range {} ... {}".format(
            number(min(arb.mpf(record["defect"]) for record in available), 10),
            number(max(arb.mpf(record["defect"]) for record in available), 10),
        )
    )
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) and control_ok else 1)

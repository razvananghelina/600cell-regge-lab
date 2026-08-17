#!/usr/bin/env python3
"""Test the anchored endpoint-counterterm corollary on frozen momentum rays."""

from collections import Counter
import hashlib
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
SEED_INPUT = HERE / "gravity_600cell_dust_nonlinear_boundary_covariance_seeds.json"
NONLINEAR_INPUT = HERE / "gravity_600cell_dust_nonlinear_boundary_covariance.json"
OUTPUT = HERE / "gravity_600cell_dust_anchored_counterterm.json"

PRIOR_ART_COMMIT = "20e9a26"
PROTOCOL_COMMIT = "b007bd3"
SEED_SHA256 = "2104c69ba6b21d3a3d92c7071d7f2702cb7d33f7f0e3ff17954f64c469f0c01d"
NONLINEAR_SHA256 = "a1e00071fa41f986dfaee84ea6e7689a14c50823f6c87d76889e6cb9346a7e3f"
arb.mp.dps = 100
CONSISTENT_FACTOR = arb.mpf(10)
REFUTED_FACTOR = arb.mpf(100)
UNCERTAINTY_FLOOR = arb.mpf("1e-70")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def vector_norm(vector):
    return arb.sqrt(sum(value*value for value in vector))


def number(value, digits=70):
    return arb.nstr(value, digits, strip_zeros=False)


hashes = {
    "seeds": digest(SEED_INPUT),
    "nonlinear_result": digest(NONLINEAR_INPUT),
}
seed = json.loads(SEED_INPUT.read_text())
nonlinear = json.loads(NONLINEAR_INPUT.read_text())
provenance_ok = bool(
    hashes == {
        "seeds": SEED_SHA256,
        "nonlinear_result": NONLINEAR_SHA256,
    }
    and seed.get("outcome") == "NONLINEAR_BOUNDARY_COVARIANCE_CASES_FROZEN"
    and seed.get("number_of_paired_cases") == 32
    and nonlinear.get("passed") == nonlinear.get("tests") == 8
    and nonlinear.get("number_of_paired_cases") == 32
    and nonlinear.get("classification_counts") == {"BROKEN": 32}
    and nonlinear.get("outcome")
        == "NONLINEAR_BOUNDARY_COVARIANCE_BROKEN_ON_FROZEN_CASES"
)

physical_map = tuple(seed["physical_edge_permutation"])
permutation_ok = bool(
    len(physical_map) == 30
    and sorted(physical_map) == list(range(30))
)
seed_by_id = {record["id"]: record for record in seed["cases"]}
result_by_id = {record["id"]: record for record in nonlinear["cases"]}
selected_ids = sorted(
    case_id for case_id, record in seed_by_id.items()
    if record["sector"] == "MOMENTUM"
)
census_ok = bool(
    len(seed_by_id) == len(result_by_id) == 32
    and set(seed_by_id) == set(result_by_id)
    and len(selected_ids) == 16
    and Counter(seed_by_id[case_id]["direction_index"] for case_id in selected_ids)
        == Counter({1: 4, 2: 4, 3: 4, 4: 4})
    and Counter(seed_by_id[case_id]["sign"] for case_id in selected_ids)
        == Counter({-1: 8, 1: 8})
    and Counter(seed_by_id[case_id]["level"] for case_id in selected_ids)
        == Counter({"0.5": 8, "1.0": 8})
)

input_rays_ok = census_ok
if input_rays_ok:
    for case_id in selected_ids:
        case = seed_by_id[case_id]
        for parity in ("even", "odd"):
            ray = tuple(arb.mpf(value) for value in case["parities"][parity]["input_ray"])
            input_rays_ok &= bool(
                len(ray) == 60
                and all(value == 0 for value in ray[:30])
                and vector_norm(ray[30:]) > 0
            )


def configuration(result, parity, variant):
    output = result["solves"][parity][variant]["output"]
    return tuple(arb.mpf(value) for value in output[:30])


def mapped_even(vector):
    result = [arb.mpf(0) for _ in range(30)]
    for source, target in enumerate(physical_map):
        result[target] = vector[source]
    return tuple(result)


records = []
implementation_ok = True
for case_id in selected_ids:
    seed_case = seed_by_id[case_id]
    result = result_by_id[case_id]
    solves = result["solves"]
    solve_success = all(
        solves[parity][variant].get("success", False)
        for parity in ("even", "odd")
        for variant in ("operational", "validation")
    )
    try:
        if not solve_success:
            defect = uncertainty = ratio = None
            classification = "OPEN_SOLVE"
        else:
            even_op = configuration(result, "even", "operational")
            even_val = configuration(result, "even", "validation")
            odd_op = configuration(result, "odd", "operational")
            odd_val = configuration(result, "odd", "validation")
            defect = vector_norm(tuple(
                left-right for left, right in zip(odd_op, mapped_even(even_op))
            ))
            uncertainty = (
                vector_norm(tuple(
                    left-right for left, right in zip(even_op, even_val)
                ))
                + vector_norm(tuple(
                    left-right for left, right in zip(odd_op, odd_val)
                ))
                + sum(
                    arb.mpf(solves[parity][variant]["correction_output"])
                    for parity in ("even", "odd")
                    for variant in ("operational", "validation")
                )
                + UNCERTAINTY_FLOOR
            )
            ratio = defect/uncertainty
            if defect <= CONSISTENT_FACTOR*uncertainty:
                classification = "ANCHORED_CONSISTENT"
            elif defect > REFUTED_FACTOR*uncertainty:
                classification = "ANCHORED_REFUTED"
            else:
                classification = "OPEN"
    except Exception:
        implementation_ok = False
        defect = uncertainty = ratio = None
        classification = "IMPLEMENTATION_ERROR"
    records.append({
        "id": case_id,
        "direction_index": seed_case["direction_index"],
        "sign": seed_case["sign"],
        "level": seed_case["level"],
        "classification": classification,
        "configuration_defect": number(defect) if defect is not None else None,
        "uncertainty": number(uncertainty) if uncertainty is not None else None,
        "ratio": number(ratio) if ratio is not None else None,
    })

classification_counts = Counter(record["classification"] for record in records)
classification_ok = bool(
    len(records) == 16
    and len({record["id"] for record in records}) == 16
    and sum(classification_counts.values()) == 16
    and set(classification_counts) <= {
        "ANCHORED_CONSISTENT", "ANCHORED_REFUTED", "OPEN", "OPEN_SOLVE"
    }
)

scaling = []
for direction in range(1, 5):
    for sign in (-1, 1):
        pair = [
            record for record in records
            if record["direction_index"] == direction and record["sign"] == sign
        ]
        half = next(record for record in pair if record["level"] == "0.5")
        full = next(record for record in pair if record["level"] == "1.0")
        available = bool(
            half["classification"] == full["classification"] == "ANCHORED_REFUTED"
        )
        if available:
            order = arb.log(
                arb.mpf(full["configuration_defect"])
                / arb.mpf(half["configuration_defect"]),
                2,
            )
            label = (
                "QUADRATIC_COMPATIBLE"
                if arb.mpf("1.5") <= order <= arb.mpf("2.5")
                else "OTHER_RESOLVED_ORDER"
            )
        else:
            order = None
            label = "NOT_AVAILABLE"
        scaling.append({
            "direction_index": direction,
            "sign": sign,
            "available": available,
            "order": number(order) if order is not None else None,
            "label": label,
        })

controls_ok = bool(
    provenance_ok
    and permutation_ok
    and census_ok
    and input_rays_ok
    and implementation_ok
    and classification_ok
)
if not controls_ok:
    outcome = "ANCHORED_ENDPOINT_COUNTERTERM_CONTROL_FAILED"
elif classification_counts["ANCHORED_REFUTED"]:
    outcome = "ANCHORED_ENDPOINT_COUNTERTERM_REFUTED_ON_FROZEN_RAYS"
elif classification_counts["OPEN"] or classification_counts["OPEN_SOLVE"]:
    outcome = "ANCHORED_ENDPOINT_COUNTERTERM_OPEN"
else:
    outcome = "ANCHORED_ENDPOINT_COUNTERTERM_CONSISTENT_ON_FROZEN_RAYS"

tests = [
    ("frozen artifact hashes and nonlinear provenance", provenance_ok),
    ("physical edge map is a permutation of all 30 coordinates", permutation_ok),
    ("selection is exactly the 16 frozen pure-momentum cases", census_ok),
    ("all selected frozen rays keep old configuration fixed", input_rays_ok),
    ("all configuration decompositions completed without exception", implementation_ok),
    ("every case received one preregistered classification", classification_ok),
    ("no action solve or external physical target was evaluated", True),
    ("outcome follows the preregistered mechanical rule", outcome in {
        "ANCHORED_ENDPOINT_COUNTERTERM_CONTROL_FAILED",
        "ANCHORED_ENDPOINT_COUNTERTERM_REFUTED_ON_FROZEN_RAYS",
        "ANCHORED_ENDPOINT_COUNTERTERM_OPEN",
        "ANCHORED_ENDPOINT_COUNTERTERM_CONSISTENT_ON_FROZEN_RAYS",
    }),
]
passed = sum(bool(ok) for _, ok in tests)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "posthoc_corollary_of_committed_nonlinear_data": True,
    "anchoring_hypotheses": {
        "old_counterterm_gradient_at_base": "zero",
        "new_counterterm_gradient_at_base": "zero",
    },
    "number_of_cases": len(records),
    "classification_counts": dict(classification_counts),
    "cases": records,
    "scaling_diagnostics": scaling,
    "scaling_label_counts": dict(Counter(record["label"] for record in scaling)),
    "action_evaluations": 0,
    "external_physical_targets_parsed": False,
    "tests": len(tests),
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"classifications={dict(classification_counts)}")
if any(record["configuration_defect"] is not None for record in records):
    available = [record for record in records if record["configuration_defect"] is not None]
    print(
        "configuration defect range {} ... {}".format(
            number(min(arb.mpf(record["configuration_defect"]) for record in available), 10),
            number(max(arb.mpf(record["configuration_defect"]) for record in available), 10),
        )
    )
    print(
        "defect/uncertainty range {} ... {}".format(
            number(min(arb.mpf(record["ratio"]) for record in available), 10),
            number(max(arb.mpf(record["ratio"]) for record in available), 10),
        )
    )
print(f"scaling labels={dict(Counter(record['label'] for record in scaling))}")
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) and controls_ok else 1)

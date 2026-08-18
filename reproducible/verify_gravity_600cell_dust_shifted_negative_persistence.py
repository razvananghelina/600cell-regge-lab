#!/usr/bin/env python3
"""Target-disclosed rank gate for shifted negative-shape persistence.

Prior-art commit: 33da8dd.
Protocol commit: 62810b4.
Open-vs-changed correction commit: c9ccf0f.

The blind shifted census was committed at 5b474c2 before this comparison.
No spectral projector or reduced product is constructed unless the complete
shifted 15+10 sign split is resolved.
"""

from collections import Counter
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
OLD_JSON = HERE / "gravity_600cell_dust_shape_stiffness.json"
OLD_SOURCE = HERE / "verify_gravity_600cell_dust_shape_stiffness.py"
SHIFTED_JSON = HERE / "gravity_600cell_dust_shifted_shape_stiffness.json"
SHIFTED_SOURCE = HERE / "verify_gravity_600cell_dust_shifted_shape_stiffness.py"
OUTPUT = HERE / "gravity_600cell_dust_shifted_negative_persistence.json"

PRIOR_ART_COMMIT = "33da8dd"
PROTOCOL_COMMIT = "62810b4"
BLIND_RESULT_COMMIT = "5b474c2"
OUTCOME_CORRECTION_COMMIT = "c9ccf0f"
EXPECTED_HASHES = {
    "old_json": "03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868",
    "old_source": "d4f0a9a805910de37011ba70f407907daa2d11c650aeea22e571ab867282a44c",
    "shifted_json": "14fe5bc91e3ae4712c6ea19b8120785e2facd364e1ceb194009123fa353a4315",
    "shifted_source": "031d0dd1cab45d0093015fcab7ce7b56e098a5742895eed71f5a531aee31c2a6",
}
PARITIES = ("even", "odd")
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
TARGET_SECTORS = (4, 5)
TARGET_NEGATIVE_RANK = 15
TARGET_POSITIVE_RANK = 10
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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sf(value):
    return f"{float(value):.17e}"


def cells(payload):
    for parity in PARITIES:
        sectors = payload["parities"][parity]
        for sector in sectors:
            sector_index = int(sector["sector_index"])
            dimension = int(sector["irrep_dimension"])
            for variant_record in sector["variants"]:
                yield parity, sector_index, dimension, variant_record


print("=" * 78)
print("SHIFTED NEGATIVE-SHAPE PERSISTENCE: TARGET-DISCLOSED RANK GATE")
print("=" * 78)

paths = {
    "old_json": OLD_JSON,
    "old_source": OLD_SOURCE,
    "shifted_json": SHIFTED_JSON,
    "shifted_source": SHIFTED_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
old = json.loads(OLD_JSON.read_text())
shifted = json.loads(SHIFTED_JSON.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and old["outcome"] == "SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED"
    and old["passed"] == old["tests"] == 12
    and shifted["outcome"] == "SHIFTED_SHAPE_STIFFNESS_SIGN_OPEN"
    and shifted["passed"] == shifted["tests"] == 12
    and old["enumeration"]["shape_pencils"] == 56
    and shifted["enumeration"]["shape_pencils"] == 56
)
check("blind old and shifted inputs have exact frozen provenance", provenance_ok, str(hashes))

old_cells = list(cells(old))
shifted_cells = list(cells(shifted))
enumeration_ok = bool(
    len(old_cells) == len(shifted_cells) == 56
    and {(p, s, d, v["variant"]) for p, s, d, v in old_cells}
    == {(p, s, d, v["variant"]) for p, s, d, v in shifted_cells}
)
check("both blind artifacts enumerate the same complete 56 cells", enumeration_ok)


def aggregate_labels(entries):
    counter = Counter()
    physical = Counter()
    for _, _, dimension, record in entries:
        labels = record["hermitian_pencil"]["A_sign_labels"]
        counter.update(labels)
        for label in labels:
            physical[label] += dimension
    return counter, physical


old_direct, old_physical = aggregate_labels(old_cells)
shifted_direct, shifted_physical = aggregate_labels(shifted_cells)
ledger_ok = bool(
    dict(old_direct) == old["pencil_sign_counts"]
    and dict(old_physical) == old["full_multiplicity_pencil_sign_counts"]
    and dict(shifted_direct) == shifted["pencil_sign_counts"]
    and dict(shifted_physical)
    == shifted["full_multiplicity_pencil_sign_counts"]
)
check("all sign labels independently reproduce both aggregate ledgers", ledger_ok)

old_by_key = {
    (p, s, v["variant"]): (d, v)
    for p, s, d, v in old_cells
}
shifted_by_key = {
    (p, s, v["variant"]): (d, v)
    for p, s, d, v in shifted_cells
}

records = []
status_counts = Counter()
old_target_ok = True
midpoint_negative_counts = Counter()
maximum_ordered_midpoint_distance = 0.0
minimum_error_ratio = math.inf
maximum_error_ratio = 0.0

for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for variant in VARIANTS:
            key = (parity, sector_index, variant)
            old_dimension, old_record = old_by_key[key]
            new_dimension, new_record = shifted_by_key[key]
            old_h = old_record["hermitian_pencil"]
            new_h = new_record["hermitian_pencil"]
            old_counts = Counter(old_h["A_sign_labels"])
            new_counts = Counter(new_h["A_sign_labels"])
            cell_dimension = len(new_h["A_sign_labels"])
            old_target_this = bool(
                old_dimension == new_dimension == 1
                and cell_dimension == 25
                and old_counts["NEGATIVE_RESOLVED"] == TARGET_NEGATIVE_RANK
                and old_counts["POSITIVE_RESOLVED"] == TARGET_POSITIVE_RANK
                and old_counts["OPEN"] == 0
                and old_counts["ZERO_CONSISTENT"] == 0
            )
            old_target_ok &= old_target_this

            negative_lower = new_counts["NEGATIVE_RESOLVED"]
            negative_upper = cell_dimension - new_counts["POSITIVE_RESOLVED"]
            exact_target = bool(
                negative_lower == negative_upper == TARGET_NEGATIVE_RANK
                and new_counts["POSITIVE_RESOLVED"] == TARGET_POSITIVE_RANK
                and new_counts["OPEN"] == 0
                and new_counts["ZERO_CONSISTENT"] == 0
            )
            if not (negative_lower <= TARGET_NEGATIVE_RANK <= negative_upper):
                status = "RANK_CHANGED_RESOLVED"
            elif exact_target:
                status = "RANK_PERSISTENCE_RESOLVED"
            else:
                status = "RANK_PERSISTENCE_OPEN"
            status_counts[status] += 1

            old_values = [float(x) for x in old_h["A_eigenvalues"]]
            new_values = [float(x) for x in new_h["A_eigenvalues"]]
            midpoint_negative = sum(x < 0 for x in new_values)
            midpoint_negative_counts[midpoint_negative] += 1
            ordered_distance = max(
                abs(left - right) for left, right in zip(old_values, new_values)
            )
            maximum_ordered_midpoint_distance = max(
                maximum_ordered_midpoint_distance, ordered_distance
            )
            old_error = float(old_h["restricted_A_error"])
            new_error = float(new_h["restricted_A_error"])
            error_ratio = new_error / old_error
            minimum_error_ratio = min(minimum_error_ratio, error_ratio)
            maximum_error_ratio = max(maximum_error_ratio, error_ratio)
            gap = new_values[TARGET_NEGATIVE_RANK] - new_values[TARGET_NEGATIVE_RANK - 1]

            records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "old_sign_counts": dict(old_counts),
                "shifted_sign_counts": dict(new_counts),
                "negative_rank_lower": negative_lower,
                "negative_rank_upper": negative_upper,
                "status": status,
                "shifted_midpoint_negative_count": midpoint_negative,
                "shifted_midpoint_minimum": sf(min(new_values)),
                "shifted_midpoint_cluster_gap": sf(gap),
                "ordered_midpoint_spectrum_distance": sf(ordered_distance),
                "old_restricted_error": sf(old_error),
                "shifted_restricted_error": sf(new_error),
                "shifted_to_old_error_ratio": sf(error_ratio),
            })

check("all 16 disclosed old cells have the exact certified 15+10 split", old_target_ok)
rank_gate_complete = bool(
    len(records) == 16
    and sum(status_counts.values()) == 16
    and all(
        record["negative_rank_lower"]
        <= TARGET_NEGATIVE_RANK
        <= record["negative_rank_upper"]
        for record in records
    )
)
check("all 16 shifted target rank intervals are classified", rank_gate_complete, str(dict(status_counts)))

off_target_resolved_negative = sum(
    Counter(record["hermitian_pencil"]["A_sign_labels"])["NEGATIVE_RESOLVED"]
    for _, sector, _, record in shifted_cells
    if sector not in TARGET_SECTORS
)
check(
    "the blind shifted census has no unreported resolved-negative sector",
    off_target_resolved_negative == 0,
    f"off_target_resolved_negative={off_target_resolved_negative}",
)

if not (provenance_ok and enumeration_ok and ledger_ok and old_target_ok):
    outcome = "SHIFTED_NEGATIVE_PERSISTENCE_CONTROL_FAILED"
elif status_counts["RANK_CHANGED_RESOLVED"]:
    outcome = "SHIFTED_NEGATIVE_RANK_OR_SECTOR_CHANGED"
elif status_counts["RANK_PERSISTENCE_OPEN"]:
    outcome = "SHIFTED_NEGATIVE_RANK_OR_SECTOR_OPEN"
else:
    # The frozen hashes do not enter this branch.  A future exact-rank result
    # requires the projector/dynamics stages specified in the protocol rather
    # than silently declaring persistence here.
    outcome = "SHIFTED_NEGATIVE_PROJECTOR_STAGE_REQUIRED"

projector_comparisons_performed = 0
reduced_products_constructed = 0
gate_obeyed = bool(
    outcome == "SHIFTED_NEGATIVE_RANK_OR_SECTOR_OPEN"
    and status_counts == Counter({"RANK_PERSISTENCE_OPEN": 16})
    and projector_comparisons_performed == 0
    and reduced_products_constructed == 0
)
check("the preregistered open branch stops before projector fitting", gate_obeyed, outcome)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "blind_shifted_result_commit": BLIND_RESULT_COMMIT,
    "open_vs_changed_correction_commit": OUTCOME_CORRECTION_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "target": {
        "sectors": list(TARGET_SECTORS),
        "negative_rank_per_sector": TARGET_NEGATIVE_RANK,
        "positive_rank_per_sector": TARGET_POSITIVE_RANK,
        "cells": 16,
    },
    "old_direct_sign_counts": dict(old_direct),
    "shifted_direct_sign_counts": dict(shifted_direct),
    "old_full_multiplicity_sign_counts": dict(old_physical),
    "shifted_full_multiplicity_sign_counts": dict(shifted_physical),
    "target_status_counts": dict(status_counts),
    "target_records": records,
    "midpoint_diagnostic": {
        "classification": "POST-BLIND STRUCTURAL PATTERN ONLY",
        "negative_count_counts": dict(midpoint_negative_counts),
        "maximum_ordered_spectrum_distance": sf(maximum_ordered_midpoint_distance),
        "shifted_to_old_error_ratio_minimum": sf(minimum_error_ratio),
        "shifted_to_old_error_ratio_maximum": sf(maximum_error_ratio),
        "used_to_select_outcome": False,
    },
    "projector_comparisons_performed": projector_comparisons_performed,
    "reduced_products_constructed": reduced_products_constructed,
    "interpretation_limits": {
        "negative_rank_15_refuted": False,
        "negative_rank_15_certified_at_shifted_tick": False,
        "common_carrier_certified": False,
        "bundle_rotation_measured": False,
        "physical_wave_mode_claimed": False,
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("SCIENTIFIC OUTCOME:", outcome)
print("target rank statuses:", dict(status_counts))
print("shifted midpoint negative counts:", dict(midpoint_negative_counts))
print("shifted/old error ratio:", minimum_error_ratio, "...", maximum_error_ratio)
print(f"{passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

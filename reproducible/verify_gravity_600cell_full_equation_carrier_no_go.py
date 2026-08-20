#!/usr/bin/env python3
"""Consolidate the complete full-equation scale+strut carrier no-go."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIOR = ROOT / "docs/gravity/gravity_600cell_full_equation_carrier_no_go_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_equation_carrier_no_go_protocol.md"
NONHOM = HERE / "gravity_600cell_full_scale_strut_canonical_precision_adversarial.json"
HOM_PRIMARY = HERE / "gravity_600cell_full_scale_strut_homogeneous_resolution.json"
HOM_ADVERSARIAL = HERE / "gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial_p200g.json"
POLE = HERE / "gravity_600cell_homogeneous_pole_transversality.json"
OUTPUT = HERE / "gravity_600cell_full_equation_carrier_no_go.json"

PROTOCOL_COMMIT = "5514c0a"
EXPECTED_HASHES = {
    "prior": "8c6df3380e5feef40e18b90b4a451fbfb3597f5e53cdab0d7e97be48b354892b",
    "protocol": "bf1d080f0b465cac54ae952cef2d408c80bf0f4b910e3197561e7dcfee149740",
    "nonhom": "ecf02fd76b0c1d4d95cd206c639a027400c2053bdb1850018d57ff2721861db3",
    "hom_primary": "70d7583756acdbee77893f98d57054ab074d9353a86247840cc1eb2c7b6be931",
    "hom_adversarial": "fab74a26ae940cf0e65f26a4f6f167285cc269e282c40d7a630f37d65ba7ab07",
    "pole": "d8fd2b0cd71d428d6cef5874b0cd6cf0496f174db13471bdb818a0803d182e0a",
}
INPUTS = {
    "prior": PRIOR,
    "protocol": PROTOCOL,
    "nonhom": NONHOM,
    "hom_primary": HOM_PRIMARY,
    "hom_adversarial": HOM_ADVERSARIAL,
    "pole": POLE,
}

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


print("=" * 78)
print("COMPLETE FULL-EQUATION SCALE+STRUT CARRIER NO-GO")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
nonhom = json.loads(NONHOM.read_text())
hom_primary = json.loads(HOM_PRIMARY.read_text())
hom_adversarial = json.loads(HOM_ADVERSARIAL.read_text())
pole = json.loads(POLE.read_text())

provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and nonhom["outcome"] == "NONHOMOGENEOUS_DIRECT_MINOR_REPLICATED"
    and nonhom["passed"] == nonhom["tests"] == 7
    and hom_primary["outcome"] == "HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE"
    and hom_primary["passed"] == hom_primary["tests"] == 10
    and hom_adversarial["outcome"]
    == "HOMOGENEOUS_WEAK_POLE_LINE_REPLICATED_AFTER_CONTROL_REPAIR"
    and hom_adversarial["passed"] == hom_adversarial["tests"] == 7
    and pole["outcome"] == "HOMOGENEOUS_WEAK_LINE_TRANSVERSE_TO_POLE_EQUATION"
    and pole["passed"] == pole["tests"] == 6
)
check("all no-go inputs retain frozen provenance", provenance_ok)

nonhom_count = int(nonhom["nonhomogeneous_sector_count"])
minor_count = int(nonhom["direct_minor_certificate_count"])
nonhom_coverage_ok = bool(
    nonhom_count == 12
    and minor_count == 48
    and len(nonhom["records"]) == 48
    and nonhom["classification"]["nonhomogeneous_canonical_intersection"]
    == "ZERO; ADVERSARIALLY REPLICATED"
)
check("all twelve nonhomogeneous cells have 48/48 direct certificates",
      nonhom_coverage_ok)

homogeneous_weak_dimensions = {}
homogeneous_full_dimensions = {}
homogeneous_transverse = {}
homogeneous_coverage_ok = True
for parity in ("even", "odd"):
    rank = hom_primary["rank_closure"][parity]
    weak_dimension = int(rank["D_exact_nullity"])
    transverse = bool(
        pole["parities"][parity]["determinant_certified"]
        and pole["parities"][parity]["convention_agrees"]
    )
    homogeneous_weak_dimensions[parity] = weak_dimension
    homogeneous_transverse[parity] = transverse
    homogeneous_full_dimensions[parity] = 0 if transverse else weak_dimension
    homogeneous_coverage_ok &= bool(
        weak_dimension == 1
        and rank["K_exact_nullity"] == 1
        and hom_adversarial["parities"][parity]["rank_ok"]
        and hom_adversarial["parities"][parity]["line_ok"]
        and transverse
    )
check("both unique homogeneous weak lines are pole-transverse",
      homogeneous_coverage_ok)

full_dimensions = [0] * nonhom_count + list(homogeneous_full_dimensions.values())
complete_zero_ok = bool(len(full_dimensions) == 14 and all(
    dimension == 0 for dimension in full_dimensions
))
check("all fourteen full-equation parity/sector intersections are zero",
      complete_zero_ok)

# The negative control must retain the two weak homogeneous lines if the pole
# result is replaced by a pole-null hypothesis.
pole_null_control = [0] * nonhom_count + list(homogeneous_weak_dimensions.values())
negative_control_ok = bool(
    len(pole_null_control) == 14
    and pole_null_control.count(1) == 2
    and pole_null_control.count(0) == 12
)
check("the pole-null corruption leaves exactly two homogeneous lines",
      negative_control_ok)

if not provenance_ok or not negative_control_ok:
    outcome = "FULL_EQUATION_CARRIER_NO_GO_CONTROL_FAILED"
elif not nonhom_coverage_ok or not homogeneous_coverage_ok:
    outcome = "FULL_EQUATION_CARRIER_COVERAGE_OPEN"
elif not complete_zero_ok:
    outcome = "FULL_EQUATION_CARRIER_INTERSECTION_OPEN"
else:
    outcome = "FULL_SCALE_STRUT_FULL_EQUATION_INTERSECTION_ZERO"

allowed = {
    "FULL_EQUATION_CARRIER_NO_GO_CONTROL_FAILED",
    "FULL_EQUATION_CARRIER_COVERAGE_OPEN",
    "FULL_EQUATION_CARRIER_INTERSECTION_OPEN",
    "FULL_SCALE_STRUT_FULL_EQUATION_INTERSECTION_ZERO",
}
check("the preregistered carrier no-go hierarchy assigns one verdict",
      outcome in allowed, outcome)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "coverage": {
        "parities": 2,
        "sectors_per_parity": 7,
        "nonhomogeneous_cells": nonhom_count,
        "nonhomogeneous_direct_certificates": minor_count,
        "homogeneous_cells": 2,
        "full_equation_cells": len(full_dimensions),
    },
    "homogeneous_weak_dimensions": homogeneous_weak_dimensions,
    "homogeneous_pole_transverse": homogeneous_transverse,
    "homogeneous_full_dimensions": homogeneous_full_dimensions,
    "full_equation_intersection_dimensions": full_dimensions,
    "negative_control_pole_null_dimensions": pole_null_control,
    "classification": {
        "complete_scale_strut_full_equation_intersection": (
            "ZERO IN ALL PARITY/SECTOR CELLS"
            if outcome == "FULL_SCALE_STRUT_FULL_EQUATION_INTERSECTION_ZERO"
            else "OPEN"
        ),
        "carrier_intersection_selection_route": (
            "CLOSED ON THE FIXED CARRIER"
            if outcome == "FULL_SCALE_STRUT_FULL_EQUATION_INTERSECTION_ZERO"
            else "OPEN"
        ),
        "unrestricted_canonical_map": "NOT REFUTED",
        "accepted_homogeneous_four_tick_map": "NOT REFUTED",
        "gravity_or_perturbations_in_general": "NOT TESTED",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"TOTAL: {passed}/{tests} tests PASSED")
print(f"Artifact: {OUTPUT.name}")
if passed != tests:
    raise SystemExit(1)


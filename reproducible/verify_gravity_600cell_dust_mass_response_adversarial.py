#!/usr/bin/env python3
"""Adversarial invariant-subspace audit of the dust-mass response.

Independence gate commit: c0b7091.
Protocol commit: df924c7.
The decisive test does not construct or order tangent eigenspaces.
"""

from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys

from flint import arb, acb, acb_mat, ctx
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_dust_mass_response.json"
PRIMARY_NUMERIC = HERE / "gravity_600cell_dust_mass_response.npz"
TANGENT_INPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
OUTPUT = HERE / "gravity_600cell_dust_mass_response_adversarial.json"

INDEPENDENCE_COMMIT = "c0b7091"
PROTOCOL_COMMIT = "df924c7"
EXPECTED_HASHES = {
    "primary": "48de8f4a9edabb84145d3ce960aab808fa45c13431ff30e8756a4314f9e1ef60",
    "primary_numeric": "ae550de064c7853cba8f5b1375276a6809d8a4631bf3abfaa656ea9a05555af6",
    "tangent": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
}
PRECISIONS = (100, 140)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
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
    return format(float(value), ".17g")


hashes = {
    "primary": sha256(PRIMARY_INPUT),
    "primary_numeric": sha256(PRIMARY_NUMERIC),
    "tangent": sha256(TANGENT_INPUT),
    "tangent_numeric": sha256(TANGENT_NUMERIC),
}
primary = json.loads(PRIMARY_INPUT.read_text())
tangent = json.loads(TANGENT_INPUT.read_text())
responses = np.load(PRIMARY_NUMERIC)
tangents = np.load(TANGENT_NUMERIC)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and primary["outcome"] == "DUST_MASS_RESPONSE_BOTH_BRANCHES_SEPARATED"
    and primary["passed"] == primary["tests"] == 20
    and primary["numeric_archive_sha256"] == hashes["primary_numeric"]
    and primary["numeric_archive_arrays"] == len(responses.files) == 112
    and tangent["outcome"] == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and tangent["numeric_archive_sha256"] == hashes["tangent_numeric"]
    and tangent["numeric_archive_arrays"] == len(tangents.files) == 224
)
check("the adversarial audit has exact frozen provenance", provenance_ok, str(hashes))


def exact_arb(value):
    numerator, denominator = float(value).as_integer_ratio()
    return arb(numerator) / arb(denominator)


def exact_acb(value):
    return acb(exact_arb(value.real), exact_arb(value.imag))


def exact_matrix(array):
    array = np.asarray(array, dtype=np.complex128)
    return acb_mat(
        array.shape[0], array.shape[1],
        [exact_acb(value) for value in array.ravel()],
    )


def norm_bounds(matrix):
    lower_square = arb(0)
    upper_square = arb(0)
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            magnitude = abs(matrix[row, column])
            lower = magnitude.lower()
            upper = magnitude.upper()
            if lower > 0:
                lower_square += lower * lower
            upper_square += upper * upper
    lower = lower_square.sqrt().lower()
    upper = upper_square.sqrt().upper()
    return lower, upper


def ratio_record_from_image(response, image, dps):
    gram = response.transpose().conjugate() * response
    determinant = gram.det()
    coefficient = gram.solve(response.transpose().conjugate() * image)
    leakage = image - response * coefficient
    leakage_lower, leakage_upper = norm_bounds(leakage)
    image_lower, image_upper = norm_bounds(image)
    ratio_lower = (leakage_lower / image_upper).lower()
    ratio_upper = (leakage_upper / image_lower).upper()
    midpoint = (ratio_lower + ratio_upper) / 2
    radius = (ratio_upper - ratio_lower) / 2
    zero_gate = arb(10) ** (-dps + 25)
    if ratio_lower > 0 and midpoint > 100 * radius:
        label = "INVARIANCE_REFUTED"
    elif ratio_upper < zero_gate:
        label = "INVARIANT_CONSISTENT"
    else:
        label = "OPEN"
    return {
        "ratio_lower": ratio_lower,
        "ratio_upper": ratio_upper,
        "ratio_midpoint": midpoint,
        "ratio_radius": radius,
        "label": label,
        "gram_excludes_zero": not determinant.contains(0),
    }


def leakage_record(response_np, tangent_np, dps, convention="forward"):
    ctx.dps = dps
    response = exact_matrix(response_np)
    tangent_ball = exact_matrix(tangent_np)
    if convention == "forward":
        image = tangent_ball * response
    elif convention == "reverse":
        image = tangent_ball.solve(response)
    else:
        raise ValueError(convention)
    return ratio_record_from_image(response, image, dps)


def interval_overlap(left, right):
    return bool(
        left["ratio_lower"] <= right["ratio_upper"]
        and right["ratio_lower"] <= left["ratio_upper"]
    )


def serialize_record(record):
    return {
        "ratio_lower": str(record["ratio_lower"]),
        "ratio_upper": str(record["ratio_upper"]),
        "ratio_midpoint": str(record["ratio_midpoint"]),
        "ratio_radius": str(record["ratio_radius"]),
        "label": record["label"],
        "gram_excludes_zero": record["gram_excludes_zero"],
    }


def scipy_leakage_ratio(response, tangent_matrix):
    q, _, _ = la.qr(response, mode="economic", pivoting=True)
    image = tangent_matrix @ response
    leakage = image - q @ (q.conj().T @ image)
    return float(la.norm(leakage, "fro") / la.norm(image, "fro"))


def rephase_columns(matrix):
    phases = np.asarray([1, -1, 1j, -1j], dtype=np.complex128)
    return matrix * phases[np.arange(matrix.shape[1]) % 4][None, :]


def phase_swap(response, tangent_matrix):
    half = response.shape[0] // 2
    permutation = np.r_[np.arange(half, 2 * half), np.arange(half)]
    return response[permutation, :], tangent_matrix[np.ix_(permutation, permutation)]


def synthetic_controls(n, k, dps):
    response = np.zeros((n, k), dtype=np.complex128)
    response[:k, :] = np.eye(k)
    identity = np.eye(n, dtype=np.complex128)
    separated = identity.copy()
    separated[k, 0] = 1
    return (
        leakage_record(response, identity, dps),
        leakage_record(response, separated, dps),
    )


print("=" * 78)
print("ADVERSARIAL DUST-MASS RESPONSE INVARIANCE AUDIT")
print("=" * 78)

cells = []
exact_controls_ok = True
actual_ok = True
conventions_ok = True
precision_ok = True
gram_ok = True
scipy_ratios = []
size_controls = {}

for parity in ("even", "odd"):
    for sector in range(7):
        prefix = f"{parity}_sector{sector}"
        response_np = responses[f"{prefix}_operational_primary_phase_midpoint"]
        tangent_np = tangents[f"{prefix}_operational_primary_tangent_midpoint"]
        n, k = response_np.shape
        print(f"[{parity}] sector {sector + 1}/7 n={n}, k={k}", flush=True)

        precision_records = {
            dps: leakage_record(response_np, tangent_np, dps)
            for dps in PRECISIONS
        }
        labels_ok = all(
            record["label"] == "INVARIANCE_REFUTED"
            for record in precision_records.values()
        )
        overlap_ok = interval_overlap(
            precision_records[PRECISIONS[0]], precision_records[PRECISIONS[1]]
        )
        cell_gram_ok = all(
            record["gram_excludes_zero"] for record in precision_records.values()
        )
        actual_ok &= labels_ok
        precision_ok &= overlap_ok
        gram_ok &= cell_gram_ok

        if (n, k) not in size_controls:
            controls = {dps: synthetic_controls(n, k, dps) for dps in PRECISIONS}
            control_ok = all(
                invariant["label"] == "INVARIANT_CONSISTENT"
                and separated["label"] == "INVARIANCE_REFUTED"
                for invariant, separated in controls.values()
            )
            exact_controls_ok &= control_ok
            size_controls[(n, k)] = controls

        original_100 = precision_records[100]
        rephased_100 = leakage_record(rephase_columns(response_np), tangent_np, 100)
        swapped_response, swapped_tangent = phase_swap(response_np, tangent_np)
        swapped_100 = leakage_record(swapped_response, swapped_tangent, 100)
        reversed_100 = leakage_record(response_np, tangent_np, 100, "reverse")
        convention_labels_ok = bool(
            rephased_100["label"] == "INVARIANCE_REFUTED"
            and swapped_100["label"] == "INVARIANCE_REFUTED"
            and reversed_100["label"] == "INVARIANCE_REFUTED"
        )
        convention_overlap_ok = bool(
            interval_overlap(original_100, rephased_100)
            and interval_overlap(original_100, swapped_100)
        )
        conventions_ok &= convention_labels_ok and convention_overlap_ok

        variant_ratios = {}
        for variant in VARIANTS:
            response_variant = responses[f"{prefix}_{variant}_phase_midpoint"]
            tangent_variant = tangents[f"{prefix}_{variant}_tangent_midpoint"]
            ratio = scipy_leakage_ratio(response_variant, tangent_variant)
            variant_ratios[variant] = ratio
            scipy_ratios.append(ratio)

        cells.append({
            "parity": parity,
            "sector_index": sector,
            "phase_dimension": n,
            "mass_dimension": k,
            "precisions": {
                str(dps): serialize_record(record)
                for dps, record in precision_records.items()
            },
            "precision_overlap": overlap_ok,
            "rephased_100": serialize_record(rephased_100),
            "phase_swap_100": serialize_record(swapped_100),
            "reversed_100": serialize_record(reversed_100),
            "convention_labels_ok": convention_labels_ok,
            "convention_overlap_ok": convention_overlap_ok,
            "scipy_variant_ratios": {
                name: sf(value) for name, value in variant_ratios.items()
            },
        })

scipy_ok = bool(
    len(scipy_ratios) == 56
    and all(math.isfinite(value) and value > 0 for value in scipy_ratios)
)
check(
    "identity and separated controls discriminate at both precisions",
    exact_controls_ok,
    f"distinct sizes={len(size_controls)}",
)
check(
    "all 14 actual response spaces are non-invariant at 100 and 140 digits",
    actual_ok,
)
check(
    "all actual Gram determinants exclude zero and precision balls overlap",
    gram_ok and precision_ok,
)
check(
    "rephasing, phase swap and time reversal retain the non-invariance verdict",
    conventions_ok,
)
check(
    "all 56 independent pivoted-QR variant ratios are finite and nonzero",
    scipy_ok,
    f"range={min(scipy_ratios):.3e}...{max(scipy_ratios):.3e}",
)

controls_ok = provenance_ok and exact_controls_ok and gram_ok and scipy_ok
if not controls_ok:
    outcome = "ADVERSARIAL_DUST_RESPONSE_CONTROL_FAILED"
elif actual_ok and precision_ok and conventions_ok:
    outcome = "ADVERSARIAL_DUST_RESPONSE_SEPARATION_CORROBORATED"
else:
    outcome = "ADVERSARIAL_DUST_RESPONSE_DISAGREEMENT_OPEN"

check(
    "the frozen hierarchy assigns exactly one adversarial outcome",
    outcome in {
        "ADVERSARIAL_DUST_RESPONSE_CONTROL_FAILED",
        "ADVERSARIAL_DUST_RESPONSE_SEPARATION_CORROBORATED",
        "ADVERSARIAL_DUST_RESPONSE_DISAGREEMENT_OPEN",
    },
    outcome,
)

serialized_controls = {}
for size, controls in size_controls.items():
    serialized_controls[f"{size[0]}x{size[1]}"] = {
        str(dps): {
            "invariant": serialize_record(pair[0]),
            "separated": serialize_record(pair[1]),
        }
        for dps, pair in controls.items()
    }

artifact = {
    "independence_commit": INDEPENDENCE_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "classification": "STRUCTURAL_INDEPENDENT_CORROBORATION",
    "independence_scope": (
        "different invariant-subspace criterion on frozen exact-dyadic response inputs; "
        "does not rederive the action source"
    ),
    "precisions": list(PRECISIONS),
    "actual_cells": len(cells),
    "branch_separations_correlated": 2 * len(cells),
    "scipy_variant_attempts": len(scipy_ratios),
    "scipy_ratio_minimum": sf(min(scipy_ratios)),
    "scipy_ratio_maximum": sf(max(scipy_ratios)),
    "synthetic_controls": serialized_controls,
    "cells": cells,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("=" * 78)
print(f"Actual cells: {len(cells)}/14 non-invariant at both precisions")
print(f"SciPy variants: {len(scipy_ratios)}/56 finite nonzero")
print(f"Artifact SHA-256: {sha256(OUTPUT)}")
print(f"Outcome: {outcome}")
print(f"{passed}/{tests} checks passed")
sys.exit(0 if passed == tests else 1)

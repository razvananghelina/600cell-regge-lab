#!/usr/bin/env python3
"""High-precision resolution of the full scale--strut Gram disagreement."""

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import linalg as sla


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_full_scale_strut_precision.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_precision_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_precision_protocol.md"
CORRECTION = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_precision_correction_protocol.md"
FIRST_FAILURE_NOTE = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_precision_first_failure.md"
SERIALIZATION_FAILURE = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_precision_serialization_failure.md"
FIRST_RESULT = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_carrier_first_result.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_carrier.py"
PRIMARY_JSON = HERE / "gravity_600cell_full_scale_strut_carrier.json"
FIRST_FAILURE_JSON = HERE / "gravity_600cell_full_scale_strut_precision_first_failure.json"

PROTOCOL_COMMIT = "624fc96"
EXPECTED_HASHES = {
    "prior_art": "4eb4556e2c38671554db8eece1c4701fa6099dd56624c2803110a4ff9c09d015",
    "protocol": "603aa2bd2c54de143df3598b7d5d03cac07338de51d5b646d068dcef2498d7e2",
    "correction": "7fb1fe2a4a5a2785ba485283e8f5958b40e53b43b6f1f763d7b205cad6cb8394",
    "first_failure_note": "400d6c8565c5deda57cf638f5af7802d1f789e257d1be6ea192e1e3e9a491faa",
    "serialization_failure": "d90f0e375c43893d806ba1bf6bcace067998756790bafddeee8789c61b215118",
    "first_result": "5753375ca2a6c4f5152f134474176501b580a1c55b7a871b3a39fa6321d82f61",
    "primary_source": "e68105df4058f7d2ed39a6913f29e88cd9fe88e123ff52260acf698a2bd7da49",
    "primary_json": "6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d",
    "first_failure_json": "23199cf8da5ed4b41d3022174e75e3035e85ddb1af8b2b9ba5aadf03132d2c68",
}
INPUTS = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "correction": CORRECTION,
    "first_failure_note": FIRST_FAILURE_NOTE,
    "serialization_failure": SERIALIZATION_FAILURE,
    "first_result": FIRST_RESULT,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "first_failure_json": FIRST_FAILURE_JSON,
}

VERTEX_COUNT = 120
DATA_COUNT = 240
ROW_COUNT = 1560
mp.mp.dps = 80
tests = 0
passed = 0


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
    return sha256(path.read_bytes()).hexdigest()


def json_default(value):
    """Repair only NumPy's non-standard boolean scalar serialization."""
    if isinstance(value, np.bool_):
        return bool(value)
    raise TypeError(f"unsupported JSON scalar: {type(value).__module__}.{type(value).__name__}")


def relative(left, right, floor=0.0):
    return abs(left - right) / max(floor, abs(left), abs(right))


def exact_float(value):
    numerator, denominator = float(value).as_integer_ratio()
    return mp.mpf(numerator) / denominator


def build_carrier(record, background):
    lam = mp.mpf(background["lambda"])
    rho = mp.mpf(background["rho"])
    L0_square = mp.mpf(background["L0_square"])
    q_diag = lam * L0_square - rho
    A = -16 * rho / (L0_square * (lam - 1) ** 2)
    B = 8 + 16 * rho / (L0_square * (lam - 1) ** 2)
    scale_factor = L0_square / (8 * q_diag)
    kappa = rho / ((lam - 1) * q_diag)
    matrix = np.zeros((ROW_COUNT, DATA_COUNT), dtype=np.float64)
    sparse_rows = []
    internal_edges = tuple(tuple(map(int, edge)) for edge in record["internal_edge_order"])
    final_edges = tuple(tuple(map(int, edge)) for edge in record["final_edge_order"])
    for row, edge in enumerate(internal_edges):
        lower, upper = edge
        sparse = {}
        if upper == lower + VERTEX_COUNT:
            sparse[VERTEX_COUNT + lower] = 1.0
        else:
            target = upper - VERTEX_COUNT
            sparse = {
                lower: float(scale_factor * A),
                target: float(scale_factor * B),
                VERTEX_COUNT + lower: float(kappa),
                VERTEX_COUNT + target: float(-lam * kappa),
            }
        for column, value in sparse.items():
            matrix[row, column] = value
        sparse_rows.append(sparse)
    for offset, edge in enumerate(final_edges):
        row = len(internal_edges) + offset
        left = edge[0] - VERTEX_COUNT
        right = edge[1] - VERTEX_COUNT
        sparse = {left: float(1 / lam), right: float(1 / lam)}
        for column, value in sparse.items():
            matrix[row, column] = value
        sparse_rows.append(sparse)
    return matrix, sparse_rows, internal_edges


def binary_spectra(matrix):
    gesdd = sla.svd(
        matrix, compute_uv=False, overwrite_a=False,
        check_finite=True, lapack_driver="gesdd",
    )
    gesvd = sla.svd(
        matrix, compute_uv=False, overwrite_a=False,
        check_finite=True, lapack_driver="gesvd",
    )
    gram_eigenvalues = np.linalg.eigvalsh(matrix.T @ matrix)[::-1]
    gram = np.sqrt(np.maximum(0.0, gram_eigenvalues))
    return gesdd, gesvd, gram, gram_eigenvalues


def high_precision_gram_spectrum(sparse_rows, size=DATA_COUNT):
    gram = mp.matrix(size, size)
    for sparse in sparse_rows:
        exact = [(column, exact_float(value)) for column, value in sparse.items()]
        for left, left_value in exact:
            for right, right_value in exact:
                gram[left, right] += left_value * right_value
    eigenvalues = list(mp.eigsy(gram, eigvals_only=True))
    eigenvalues.sort(reverse=True)
    singular = [mp.sqrt(value) if value >= 0 else -mp.sqrt(-value) for value in eigenvalues]
    return singular, eigenvalues


def maximum_relative(left, right):
    return max(
        relative(float(a), float(b), np.finfo(float).tiny)
        for a, b in zip(left, right)
    )


def synthetic_control():
    delta = 1e-7
    c = 1 / np.sqrt(2.0)
    rotation = np.array([[c, -c], [c, c]], dtype=np.float64)
    matrix = np.diag([1.0, delta]) @ rotation
    gesdd, gesvd, gram, _ = binary_spectra(matrix)
    sparse = [
        {column: float(matrix[row, column]) for column in range(2) if matrix[row, column]}
        for row in range(2)
    ]
    high, eigenvalues = high_precision_gram_spectrum(sparse, size=2)
    truth = np.array([1.0, delta])
    direct_truth_error = max(
        maximum_relative(gesdd, truth), maximum_relative(gesvd, truth)
    )
    binary_gram_error = maximum_relative(gram, truth)
    high_direct_error = maximum_relative(high, gesvd)
    passed_control = bool(
        direct_truth_error < 1e-8
        and binary_gram_error > 1e-4
        and high_direct_error < 1e-8
        and all(value > 0 for value in eigenvalues)
    )
    return passed_control, {
        "delta": delta,
        "known_singular_values": [f"{value:.17e}" for value in truth],
        "gesdd": [f"{value:.17e}" for value in gesdd],
        "gesvd": [f"{value:.17e}" for value in gesvd],
        "binary_gram": [f"{value:.17e}" for value in gram],
        "high_precision_gram": [mp.nstr(value, 50) for value in high],
        "direct_truth_max_relative_error": f"{direct_truth_error:.17e}",
        "binary_gram_max_relative_error": f"{binary_gram_error:.17e}",
        "high_direct_max_relative_error": f"{high_direct_error:.17e}",
    }


print("=" * 78)
print("HIGH-PRECISION FULL SCALE--STRUT CARRIER AUDIT")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all precision-audit inputs retain frozen provenance", provenance_ok, str(hashes))

primary = json.loads(PRIMARY_JSON.read_text())
first_failure = json.loads(FIRST_FAILURE_JSON.read_text())
frozen_ok = bool(
    primary.get("outcome") == "FULL_SCALE_STRUT_NUMERICALLY_OPEN"
    and primary.get("passed") == primary.get("tests") == 18
    and primary.get("finite_formula_agrees") is True
    and all(
        record["candidate_mismatch_count"] == 0
        and record["full_residual_nonzero_count"] == 0
        for record in primary["new_exact_global_controls"]
    )
    and all(
        primary["parities"][parity]["complete_exact_rank"] == DATA_COUNT
        and 0.02 < float(primary["parities"][parity]["spectral_diagnostic"]["maximum_relative_discrepancy"]) < 0.04
        for parity in ("even", "odd")
    )
)
check("the frozen first artifact retains its exact positive and numeric-open content", frozen_ok)

first_failure_ok = bool(
    first_failure.get("outcome") == "FULL_SCALE_STRUT_PRECISION_CONTROL_FAILED"
    and first_failure.get("passed") == 6
    and first_failure.get("tests") == 8
    and all(
        first_failure["parities"][parity]["decisive_precision_criteria"] is True
        and first_failure["parities"][parity]["gesdd_gesvd_max_relative_discrepancy"] == "0.00000000000000000e+00"
        and first_failure["parities"][parity]["high_precision_all_positive"] is True
        for parity in ("even", "odd")
    )
)
check("the first 6/8 control failure is preserved literally", first_failure_ok)

synthetic_ok, synthetic_record = synthetic_control()
check("the calibrated 1e7-condition control separates direct SVD from binary Gram", synthetic_ok, str(synthetic_record))

records = {}
all_reconstruction_ok = True
all_row_order_ok = True
all_corruption_ok = True
all_decisive_ok = True
for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing frozen carrier", flush=True)
    frozen = primary["parities"][parity]
    matrix, sparse_rows, internal_edges = build_carrier(frozen, primary["background"])
    support_histogram = Counter(len(row) for row in sparse_rows)
    reconstruction_shape = bool(
        matrix.shape == (ROW_COUNT, DATA_COUNT)
        and len(sparse_rows) == ROW_COUNT
        and support_histogram == Counter({4: 720, 1: 120, 2: 720})
        and np.all(np.isfinite(matrix))
    )
    gesdd, gesvd, binary_gram, binary_gram_eigenvalues = binary_spectra(matrix)
    condition = float(gesdd[0] / gesdd[-1])
    binary_error = maximum_relative(gesdd, binary_gram)
    frozen_condition = float(frozen["spectral_diagnostic"]["direct_condition_number"])
    frozen_error = float(frozen["spectral_diagnostic"]["maximum_relative_discrepancy"])
    reconstruction_match = bool(
        reconstruction_shape
        and relative(condition, frozen_condition) < 1e-10
        and abs(binary_error - frozen_error)
        < np.finfo(float).eps * condition**2
    )
    all_reconstruction_ok &= reconstruction_match

    reversed_singular = sla.svd(
        matrix[::-1, :], compute_uv=False, check_finite=True,
        lapack_driver="gesvd",
    )
    row_order_error = maximum_relative(reversed_singular, gesvd)
    row_order_bound = DATA_COUNT * np.finfo(float).eps * float(gesvd[0] / gesvd[-1])
    row_order_ok = row_order_error < row_order_bound
    all_row_order_ok &= row_order_ok

    corrupted = matrix.copy()
    first_pole_row = next(
        row for row, edge in enumerate(internal_edges)
        if edge[1] == edge[0] + VERTEX_COUNT
    )
    pole_vertex = internal_edges[first_pole_row][0]
    corrupted[first_pole_row, VERTEX_COUNT + pole_vertex] = 0.0
    corrupted_singular = sla.svd(
        corrupted, compute_uv=False, check_finite=True,
        lapack_driver="gesvd",
    )
    corruption_ok = bool(corrupted_singular[-1] < gesvd[-1])
    all_corruption_ok &= corruption_ok

    print(f"[{parity}] accumulating exact-binary 80-decimal Gram", flush=True)
    high, high_eigenvalues = high_precision_gram_spectrum(sparse_rows)
    driver_error = maximum_relative(gesdd, gesvd)
    high_error = maximum_relative(high, gesvd)
    positive = all(value > 0 for value in high_eigenvalues)
    positivity_margin = (
        high_eigenvalues[-1]
        / (mp.mpf("1e-80") * high_eigenvalues[0])
        if positive else mp.mpf(0)
    )
    eps_kappa_square = np.finfo(float).eps * float(gesvd[0] / gesvd[-1]) ** 2
    explanation_ratio = binary_error / eps_kappa_square
    decisive = bool(
        driver_error < 1e-8
        and high_error < 1e-8
        and positive
        and positivity_margin > mp.mpf("1e40")
        and 0.05 < explanation_ratio < 20
    )
    all_decisive_ok &= decisive
    records[parity] = {
        "matrix_shape": list(matrix.shape),
        "support_histogram": {str(key): value for key, value in sorted(support_histogram.items())},
        "reconstruction": {
            "direct_condition_number": f"{condition:.17e}",
            "frozen_condition_number": f"{frozen_condition:.17e}",
            "binary_gram_error": f"{binary_error:.17e}",
            "frozen_binary_gram_error": f"{frozen_error:.17e}",
            "matches": reconstruction_match,
        },
        "gesdd_singular_values": [f"{value:.17e}" for value in gesdd],
        "gesvd_singular_values": [f"{value:.17e}" for value in gesvd],
        "high_precision_gram_singular_values": [mp.nstr(value, 60) for value in high],
        "binary_gram_singular_values": [f"{value:.17e}" for value in binary_gram],
        "binary_gram_minimum_raw_eigenvalue": f"{binary_gram_eigenvalues[-1]:.17e}",
        "gesdd_gesvd_max_relative_discrepancy": f"{driver_error:.17e}",
        "gesvd_high_precision_max_relative_discrepancy": f"{high_error:.17e}",
        "high_precision_all_positive": positive,
        "high_precision_positivity_margin_over_1e_minus_80": mp.nstr(positivity_margin, 60),
        "epsilon_kappa_squared": f"{eps_kappa_square:.17e}",
        "old_error_over_epsilon_kappa_squared": f"{explanation_ratio:.17e}",
        "row_reversal_max_relative_discrepancy": f"{row_order_error:.17e}",
        "row_reversal_conditioned_bound": f"{row_order_bound:.17e}",
        "row_reversal_control": row_order_ok,
        "corruption": {
            "deleted_pole_row": first_pole_row,
            "baseline_smallest_singular": f"{gesvd[-1]:.17e}",
            "corrupted_smallest_singular": f"{corrupted_singular[-1]:.17e}",
            "strictly_smaller": corruption_ok,
        },
        "decisive_precision_criteria": decisive,
    }

check("both carrier matrices reproduce the frozen binary64 diagnostics", all_reconstruction_ok)
check("row-order reversal preserves both direct spectra", all_row_order_ok)
check("deleting the first pole coefficient strictly weakens both carriers", all_corruption_ok)
check("both high-precision comparisons are evaluated and classified", all(record["high_precision_gram_singular_values"] for record in records.values()), f"decisive={all_decisive_ok}")

controls_ok = bool(
    provenance_ok and frozen_ok and first_failure_ok and synthetic_ok and all_reconstruction_ok
    and all_row_order_ok and all_corruption_ok
)
if not controls_ok:
    outcome = "FULL_SCALE_STRUT_PRECISION_CONTROL_FAILED"
elif all_decisive_ok:
    outcome = "FULL_SCALE_STRUT_PRECISION_RESOLVED"
elif not all_decisive_ok:
    outcome = "FULL_SCALE_STRUT_PRECISION_DISAGREEMENT"
else:
    outcome = "FULL_SCALE_STRUT_PRECISION_OPEN"

allowed = {
    "FULL_SCALE_STRUT_PRECISION_CONTROL_FAILED",
    "FULL_SCALE_STRUT_PRECISION_RESOLVED",
    "FULL_SCALE_STRUT_PRECISION_DISAGREEMENT",
    "FULL_SCALE_STRUT_PRECISION_OPEN",
}
check("the preregistered precision hierarchy assigns one outcome", outcome in allowed, outcome)

payload = {
    "prior_art_commit": "16e4380",
    "protocol_commit": PROTOCOL_COMMIT,
    "correction_commit": "3c96ce9",
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "synthetic_control": synthetic_record,
    "parities": records,
    "classification": {
        "old_binary64_gram_disagreement": (
            "NORMAL_EQUATION_PRECISION LOSS" if outcome == "FULL_SCALE_STRUT_PRECISION_RESOLVED" else "OPEN"
        ),
        "frozen_binary64_matrix_rank": "EXACTLY 240 FROM PRIMARY COMBINATORIAL PROOF",
        "generic_geometric_formula": "OPEN PENDING SYMBOLIC ADVERSARIAL REPLICATION",
        "physical_interpretation": "NOT EVALUATED",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(
    json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n"
)

print("-" * 78)
print(outcome)
print(f"TOTAL: {passed}/{tests} tests PASSED")
if passed != tests or outcome == "FULL_SCALE_STRUT_PRECISION_CONTROL_FAILED":
    raise SystemExit(1)

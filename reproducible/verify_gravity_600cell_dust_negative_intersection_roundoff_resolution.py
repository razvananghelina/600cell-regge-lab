#!/usr/bin/env python3
"""Fixed-binary-input resolution of the negative-intersection audit.

Framing gates: cf64f34, ba42d40.
Preregistered protocol: b54f80a.
"""

from collections import Counter
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
NEGATIVE_SOURCE = HERE / "verify_gravity_600cell_dust_negative_fiber_transport.py"
NEGATIVE_ARTIFACT = HERE / "gravity_600cell_dust_negative_fiber_transport.json"
TANGENT_ARCHIVE = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"
PRIMARY_ARTIFACT = (
    HERE / "gravity_600cell_dust_negative_transported_intersection.json"
)
AUDIT_SOURCE = (
    HERE
    / "verify_gravity_600cell_dust_negative_transported_intersection_adversarial.py"
)
AUDIT_ARTIFACT = (
    HERE
    / "gravity_600cell_dust_negative_transported_intersection_adversarial.json"
)
OUTPUT = (
    HERE / "gravity_600cell_dust_negative_intersection_roundoff_resolution.json"
)

PRIOR_ART_COMMITS = ("cf64f34", "ba42d40")
PROTOCOL_COMMIT = "b54f80a"
EXPECTED_HASHES = {
    "negative_source": (
        "f462e507500d7f02ecf799f0d4b320e05795216a36a0d10eb908d6dc67b48181"
    ),
    "negative_artifact": (
        "d630bf07066f88c35eee5a62a80ec1f43399a95ea882a43528289220c67f4599"
    ),
    "tangent_archive": (
        "ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d"
    ),
    "primary_artifact": (
        "c490431bdaeae3026692cd358f60d0b47ef5d63aa59217e400daac807ed21be0"
    ),
    "audit_source": (
        "6aa7e841d31bdc87568a6e4370ed334b6f5c09884669ee13fcd68d46ea4b3162"
    ),
    "audit_artifact": (
        "d5074507326bb981ad7573bd562c1aa9f0af4e1eb6b6924e3ac959a5fa1d3340"
    ),
}

PRECISIONS = (100, 140)
PARITIES = ("even", "odd")
TIMES = ("old", "shifted")
TARGET_SECTORS = (4, 5)
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


def smp(value, digits=40):
    value = mp.mpf(value)
    if mp.isinf(value):
        return "inf" if value > 0 else "-inf"
    if mp.isnan(value):
        return "nan"
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def mp_diag(values):
    result = mp.matrix(len(values), len(values))
    for index, value in enumerate(values):
        result[index, index] = value
    return result


def mp_operator_norm(matrix):
    if matrix.rows == 0 or matrix.cols == 0:
        return mp.mpf(0)
    gram = (matrix.H * matrix + (matrix.H * matrix).H) / 2
    values = mp.eighe(gram, eigvals_only=True)
    largest = max(mp.re(values[index]) for index in range(values.rows))
    return mp.sqrt(max(mp.mpf(0), largest))


def exact_float(value):
    numerator, denominator = float(value).as_integer_ratio()
    return mp.mpf(numerator) / mp.mpf(denominator)


def exact_numpy_matrix(matrix):
    array = np.asarray(matrix)
    rows, columns = array.shape
    result = mp.matrix(rows, columns)
    for row in range(rows):
        for column in range(columns):
            value = complex(array[row, column])
            result[row, column] = mp.mpc(
                exact_float(value.real), exact_float(value.imag)
            )
    return result


def phase_basis(configuration):
    result = mp.matrix(60, 30)
    conjugate = configuration.conjugate()
    for row in range(30):
        for column in range(15):
            result[row, column] = configuration[row, column]
            result[30 + row, 15 + column] = conjugate[row, column]
    return result


def gauge_basis(basis, denominator, sign):
    result = mp.matrix(basis.rows, basis.cols)
    for column in range(basis.cols):
        source_column = basis.cols - 1 - column
        phase = mp.exp(
            sign * mp.j * mp.mpf(column) / mp.mpf(denominator)
        )
        for row in range(basis.rows):
            result[row, column] = basis[row, source_column] * phase
    return result


def direct_singulars(matrix):
    values = mp.svd(matrix, compute_uv=False)
    return [mp.mpf(mp.re(values[index])) for index in range(values.rows)]


def gram_singulars(matrix):
    gram = (matrix.H * matrix + (matrix.H * matrix).H) / 2
    values = mp.eighe(gram, eigvals_only=True)
    singulars = [
        mp.sqrt(max(mp.mpf(0), mp.re(values[index])))
        for index in range(values.rows)
    ]
    return sorted(singulars, reverse=True)


def fixed_label(value, error):
    if not mp.isfinite(value) or not mp.isfinite(error):
        return "FIXED_SINGULAR_OPEN"
    if value <= 10 * error:
        return "FIXED_SINGULAR_ZERO_CONSISTENT"
    if value > 100 * error:
        return "FIXED_SINGULAR_NONZERO_RESOLVED"
    return "FIXED_SINGULAR_OPEN"


hashes = {
    "negative_source": sha256(NEGATIVE_SOURCE),
    "negative_artifact": sha256(NEGATIVE_ARTIFACT),
    "tangent_archive": sha256(TANGENT_ARCHIVE),
    "primary_artifact": sha256(PRIMARY_ARTIFACT),
    "audit_source": sha256(AUDIT_SOURCE),
    "audit_artifact": sha256(AUDIT_ARTIFACT),
}
negative_artifact = json.loads(NEGATIVE_ARTIFACT.read_text())
primary_artifact = json.loads(PRIMARY_ARTIFACT.read_text())
audit_artifact = json.loads(AUDIT_ARTIFACT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and negative_artifact["passed"] == negative_artifact["tests"] == 8
    and negative_artifact["outcome"]
    == "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED"
    and primary_artifact["passed"] == primary_artifact["tests"] == 10
    and primary_artifact["outcome"]
    == "NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL"
    and audit_artifact["passed"] == audit_artifact["tests"] == 11
    and audit_artifact["outcome"]
    == "ADVERSARIAL_NEGATIVE_INTERSECTION_DISAGREEMENT_OPEN"
)
check("all frozen roundoff-resolution inputs retain exact provenance",
      provenance_ok, str(hashes))

print("[setup] replaying the earlier binary negative-fiber source", flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    negative = runpy.run_path(str(NEGATIVE_SOURCE))
replay_ok = bool(
    negative["passed"] == negative["tests"] == 8
    and negative["outcome"] == "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED"
    and sha256(NEGATIVE_ARTIFACT) == EXPECTED_HASHES["negative_artifact"]
    and sum(len(negative["projectors"][time]) for time in TIMES) == 32
)
check("the binary source replays byte-identically", replay_ok)

# Freeze NumPy source arrays before changing mpmath precision.
binary_projectors = {}
for time_name in TIMES:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for variant in VARIANTS:
                key = (time_name, parity, sector_index, variant)
                binary_projectors[key] = np.array(
                    negative["projectors"][time_name][
                        (parity, sector_index, variant)
                    ]["projector"],
                    dtype=np.complex128,
                    copy=True,
                )

binary_tangents = {}
with np.load(TANGENT_ARCHIVE, allow_pickle=False) as tangent_archive:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for variant in VARIANTS:
                prefix = f"{parity}_sector{sector_index}_t2_{variant}"
                binary_tangents[(parity, sector_index, variant)] = np.array(
                    tangent_archive[f"{prefix}_midpoint"],
                    dtype=np.complex128,
                    copy=True,
                )

projector_records = []
precision_cells = {}
singular_records = []
spectra = {}
errors = {}
projector_ok = True
spectra_ok = True
crosscheck_ok = True
inverse_ok = True
synthetic_ok = True
gauge_ok = True
global_counts = {precision: Counter() for precision in PRECISIONS}

for precision in PRECISIONS:
    mp.mp.dps = precision
    arithmetic_floor = mp.mpf("1e-75")
    fiber_data = {}
    for time_name in TIMES:
        for parity in PARITIES:
            for sector_index in TARGET_SECTORS:
                for variant in VARIANTS:
                    key = (time_name, parity, sector_index, variant)
                    projector = exact_numpy_matrix(binary_projectors[key])
                    projector = (projector + projector.H) / 2
                    values, vectors = mp.eighe(projector)
                    bottom = vectors[:, :15]
                    top = vectors[:, 15:]
                    gap = mp.re(values[15] - values[14])
                    top_values = [mp.re(values[index]) for index in range(15, 30)]
                    eigen_residual = mp_operator_norm(
                        projector * top - top * mp_diag(top_values)
                    )
                    eta = (
                        2 * eigen_residual / (gap - 2 * eigen_residual)
                        + arithmetic_floor
                        if gap > 2 * eigen_residual else mp.inf
                    )
                    full_orthogonality = mp_operator_norm(
                        vectors.H * vectors - mp.eye(30)
                    )
                    fiber = phase_basis(top)
                    complement = phase_basis(bottom)
                    phase_orthogonality = max(
                        mp_operator_norm(fiber.H * fiber - mp.eye(30)),
                        mp_operator_norm(
                            complement.H * complement - mp.eye(30)
                        ),
                        mp_operator_norm(complement.H * fiber),
                    )
                    ordered = all(
                        mp.re(values[index]) <= mp.re(values[index + 1])
                        for index in range(29)
                    )
                    split_ok = bool(
                        mp.re(values[14]) < mp.mpf("1e-10")
                        and mp.re(values[15]) > 1 - mp.mpf("1e-10")
                    )
                    complete = bool(
                        ordered and split_ok and gap > 2 * eigen_residual
                        and full_orthogonality < mp.mpf("1e-70")
                        and phase_orthogonality < mp.mpf("1e-70")
                        and mp.isfinite(eta)
                    )
                    projector_ok &= complete
                    fiber_data[key] = {
                        "fiber": fiber,
                        "complement": complement,
                        "eta": eta,
                    }
                    projector_records.append({
                        "precision": precision,
                        "time": time_name,
                        "parity": parity,
                        "sector_index": sector_index,
                        "variant": variant,
                        "lower_projector_edge": smp(values[14]),
                        "upper_projector_edge": smp(values[15]),
                        "spectral_gap": smp(gap),
                        "eigen_residual": smp(eigen_residual),
                        "projector_error": smp(eta),
                        "full_orthogonality": smp(full_orthogonality),
                        "phase_orthogonality": smp(phase_orthogonality),
                        "split_ok": split_ok,
                        "complete": complete,
                    })

    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for variant in VARIANTS:
                cell_key = (parity, sector_index, variant)
                old = fiber_data[("old",) + cell_key]
                shifted = fiber_data[("shifted",) + cell_key]
                tangent = exact_numpy_matrix(binary_tangents[cell_key])
                tangent_norm = mp_operator_norm(tangent)
                source = old["fiber"]
                target = shifted["fiber"]
                target_complement = shifted["complement"]
                leakage_matrix = target_complement.H * tangent * source
                error = (
                    (old["eta"] + shifted["eta"]
                     + old["eta"] * shifted["eta"])
                    * tangent_norm
                    + arithmetic_floor * max(mp.mpf(1), tangent_norm)
                )
                direct = direct_singulars(leakage_matrix)
                gram = gram_singulars(leakage_matrix)
                spectrum_difference = max(
                    abs(left - right) for left, right in zip(direct, gram)
                )
                invertible = True
                try:
                    inverse_smallest = 1 / mp_operator_norm(
                        mp.inverse(leakage_matrix)
                    )
                except (ZeroDivisionError, ValueError, ZeroDivisionError):
                    inverse_smallest = mp.mpf(0)
                    invertible = False
                inverse_difference = abs(direct[-1] - inverse_smallest)
                labels = [fixed_label(value, error) for value in direct]
                counts = Counter(labels)
                global_counts[precision].update(counts)

                full_map = target * source.H
                full_leakage = target_complement.H * full_map * source
                full_values = direct_singulars(full_leakage)
                full_labels = [fixed_label(value, error) for value in full_values]
                zero_map = target_complement * source.H
                zero_leakage = target_complement.H * zero_map * source
                zero_values = direct_singulars(zero_leakage)
                zero_labels = [fixed_label(value, error) for value in zero_values]
                cell_synthetic = bool(
                    Counter(full_labels)[
                        "FIXED_SINGULAR_ZERO_CONSISTENT"
                    ] == 30
                    and Counter(zero_labels)[
                        "FIXED_SINGULAR_NONZERO_RESOLVED"
                    ] == 30
                )
                synthetic_ok &= cell_synthetic

                gauged_source = gauge_basis(source, 7, 1)
                gauged_complement = gauge_basis(target_complement, 11, -1)
                gauged_matrix = (
                    gauged_complement.H * tangent * gauged_source
                )
                gauged = direct_singulars(gauged_matrix)
                gauge_difference = max(
                    abs(left - right) for left, right in zip(direct, gauged)
                )
                cell_gauge = gauge_difference <= 10 * error
                gauge_ok &= cell_gauge

                finite = all(mp.isfinite(value) for value in direct + gram)
                ordered = all(
                    direct[index] >= direct[index + 1]
                    for index in range(29)
                )
                dimension_ok = bool(
                    source.rows == target.rows == target_complement.rows == 60
                    and source.cols == target.cols
                    == target_complement.cols == 30
                    and leakage_matrix.rows == leakage_matrix.cols == 30
                )
                cell_crosscheck = spectrum_difference <= 10 * error
                cell_inverse = bool(
                    invertible and inverse_difference <= 10 * error
                )
                crosscheck_ok &= cell_crosscheck
                inverse_ok &= cell_inverse
                spectra_ok &= bool(
                    finite and ordered and dimension_ok
                    and len(direct) == len(gram) == 30
                )
                resolved_rank = counts[
                    "FIXED_SINGULAR_NONZERO_RESOLVED"
                ]
                spectra[(precision,) + cell_key] = direct
                errors[(precision,) + cell_key] = error
                precision_cells[(precision,) + cell_key] = {
                    "precision": precision,
                    "tangent_norm": smp(tangent_norm),
                    "leakage_norm": smp(direct[0]),
                    "smallest_singular": smp(direct[-1]),
                    "error": smp(error),
                    "smallest_error_units": smp(direct[-1] / error),
                    "resolved_rank": resolved_rank,
                    "zero_consistent": counts[
                        "FIXED_SINGULAR_ZERO_CONSISTENT"
                    ],
                    "open": counts["FIXED_SINGULAR_OPEN"],
                    "svd_gram_difference": smp(spectrum_difference),
                    "inverse_smallest": smp(inverse_smallest),
                    "inverse_difference": smp(inverse_difference),
                    "gauge_difference": smp(gauge_difference),
                    "full_intersection_control_rank": Counter(full_labels)[
                        "FIXED_SINGULAR_NONZERO_RESOLVED"
                    ],
                    "zero_intersection_control_rank": Counter(zero_labels)[
                        "FIXED_SINGULAR_NONZERO_RESOLVED"
                    ],
                    "complete": bool(
                        finite and ordered and dimension_ok
                        and cell_crosscheck and cell_inverse
                        and cell_synthetic and cell_gauge
                    ),
                }
                for index, (value, label) in enumerate(zip(direct, labels)):
                    singular_records.append({
                        "precision": precision,
                        "parity": parity,
                        "sector_index": sector_index,
                        "variant": variant,
                        "index": index,
                        "value": smp(value),
                        "error": smp(error),
                        "error_units": smp(value / error),
                        "label": label,
                    })

check("all 64 fixed-input projector reconstructions certify",
      projector_ok and len(projector_records) == 64)
check("all 32 per-precision leakage spectra are finite and complete",
      spectra_ok and len(precision_cells) == 32)
check("all direct SVD and Gram spectra agree within residual bounds",
      crosscheck_ok)
check("all inverse norms independently recover the smallest singular value",
      inverse_ok)
check("all fixed-input full/zero-intersection controls discriminate",
      synthetic_ok)
check("all high-precision spectra survive basis-gauge stress", gauge_ok)

cell_records = []
precision_stability_ok = True
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for variant in VARIANTS:
            cell_key = (parity, sector_index, variant)
            low = spectra[(PRECISIONS[0],) + cell_key]
            high = spectra[(PRECISIONS[1],) + cell_key]
            stability_error = (
                errors[(PRECISIONS[0],) + cell_key]
                + errors[(PRECISIONS[1],) + cell_key]
            )
            differences = [abs(left - right) for left, right in zip(low, high)]
            maximum_difference = max(differences)
            stable = maximum_difference <= 10 * stability_error
            precision_stability_ok &= stable
            cell_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "precision_100": precision_cells[(100,) + cell_key],
                "precision_140": precision_cells[(140,) + cell_key],
                "maximum_precision_difference": smp(maximum_difference),
                "precision_stability_error": smp(stability_error),
                "precision_stable": stable,
                "complete": bool(
                    stable
                    and precision_cells[(100,) + cell_key]["complete"]
                    and precision_cells[(140,) + cell_key]["complete"]
                ),
            })

check("all 16 complete spectra are stable from 100 to 140 digits",
      precision_stability_ok and len(cell_records) == 16)
census_ok = bool(
    len(singular_records) == 960
    and all(sum(global_counts[p].values()) == 480 for p in PRECISIONS)
)
check("all 960 fixed-input singular values receive frozen labels",
      census_ok, str({p: dict(global_counts[p]) for p in PRECISIONS}))

all_rank30 = all(
    precision_cells[(precision, parity, sector_index, variant)][
        "resolved_rank"
    ] == 30
    for precision in PRECISIONS
    for parity in PARITIES
    for sector_index in TARGET_SECTORS
    for variant in VARIANTS
)
controls_ok = bool(
    provenance_ok and replay_ok and projector_ok and spectra_ok
    and crosscheck_ok and inverse_ok and synthetic_ok and gauge_ok
    and precision_stability_ok and census_ok
    and len(projector_records) == 64 and len(cell_records) == 16
)
if not controls_ok:
    outcome = "NEGATIVE_INTERSECTION_ROUNDOFF_RESOLUTION_CONTROL_FAILED"
elif not all_rank30:
    outcome = "NEGATIVE_INTERSECTION_ROUNDOFF_DISAGREEMENT_REMAINS_OPEN"
else:
    outcome = "NEGATIVE_INTERSECTION_ROUNDOFF_DISAGREEMENT_RESOLVED"

allowed = {
    "NEGATIVE_INTERSECTION_ROUNDOFF_RESOLUTION_CONTROL_FAILED",
    "NEGATIVE_INTERSECTION_ROUNDOFF_DISAGREEMENT_REMAINS_OPEN",
    "NEGATIVE_INTERSECTION_ROUNDOFF_DISAGREEMENT_RESOLVED",
}
check("the preregistered roundoff-resolution hierarchy assigns one outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commits": list(PRIOR_ART_COMMITS),
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "primary_numeric_singular_values_read": False,
    "first_audit_threshold_changed": False,
    "fitted_threshold_or_alignment_used": False,
    "precisions": list(PRECISIONS),
    "arithmetic_floor": "1e-75",
    "projector_count": len(projector_records),
    "projector_records": projector_records,
    "cell_count": len(cell_records),
    "cell_records": cell_records,
    "global_label_counts": {
        str(precision): dict(global_counts[precision])
        for precision in PRECISIONS
    },
    "singular_records": singular_records,
    "classification": {
        "first_float64_audit": "INCONCLUSIVE UNDER ITS FROZEN ENVELOPE",
        "fixed_binary_input_rank": (
            "DERIVED COMPUTATIONAL RANK 30"
            if outcome
            == "NEGATIVE_INTERSECTION_ROUNDOFF_DISAGREEMENT_RESOLVED"
            else "OPEN"
        ),
        "source_certified_primary_result": (
            "ELIGIBLE FOR CONSOLIDATION"
            if outcome
            == "NEGATIVE_INTERSECTION_ROUNDOFF_DISAGREEMENT_RESOLVED"
            else "OPEN"
        ),
        "physical_pre_post_constraint_surface": "NOT ESTABLISHED",
        "graph_lagrangian_propagator_dispersion_mass_inertia_or_speed": (
            "NOT COMPUTED"
        ),
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"labels: {artifact['global_label_counts']}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

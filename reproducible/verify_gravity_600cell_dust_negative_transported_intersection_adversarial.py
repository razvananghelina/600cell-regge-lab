#!/usr/bin/env python3
"""Independent binary64 audit of the transported negative-fiber intersection.

Independence gate: 8d532e8.
Preregistered protocol: 411df3f.
"""

from collections import Counter
import contextlib
import hashlib
import io
import json
import math
from pathlib import Path
import runpy

import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
NEGATIVE_SOURCE = HERE / "verify_gravity_600cell_dust_negative_fiber_transport.py"
NEGATIVE_ARTIFACT = HERE / "gravity_600cell_dust_negative_fiber_transport.json"
TANGENT_ARCHIVE = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"
PRIMARY_ARTIFACT = (
    HERE / "gravity_600cell_dust_negative_transported_intersection.json"
)
OUTPUT = (
    HERE
    / "gravity_600cell_dust_negative_transported_intersection_adversarial.json"
)

PRIOR_ART_COMMIT = "8d532e8"
PROTOCOL_COMMIT = "411df3f"
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
}

PARITIES = ("even", "odd")
TIMES = ("old", "shifted")
TARGET_SECTORS = (4, 5)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
MACHINE_EPSILON = np.finfo(float).eps
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


def operator_norm(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    return float(singular[0]) if len(singular) else 0.0


def audit_floor(tangent, leakage):
    return float(
        1000 * MACHINE_EPSILON * 60
        * max(1.0, operator_norm(tangent), operator_norm(leakage))
    )


def singular_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error):
        return "AUDIT_SINGULAR_OPEN"
    if value <= 10 * error:
        return "AUDIT_SINGULAR_ZERO_CONSISTENT"
    if value > 100 * error:
        return "AUDIT_SINGULAR_NONZERO_RESOLVED"
    return "AUDIT_SINGULAR_OPEN"


def phase_basis(configuration_basis):
    result = np.zeros((60, 30), dtype=np.complex128)
    result[:30, :15] = configuration_basis
    result[30:, 15:] = configuration_basis.conj()
    return result


def qr_complement(basis):
    complete, _ = np.linalg.qr(basis, mode="complete")
    return complete[:, 30:]


def leakage(tangent, source, target_complement, driver="numpy"):
    matrix = target_complement.conj().T @ tangent @ source
    if driver == "numpy":
        singular = np.linalg.svd(matrix, compute_uv=False)
    elif driver == "gesvd":
        singular = la.svd(
            matrix, compute_uv=False, check_finite=True,
            lapack_driver="gesvd",
        )
    else:
        raise ValueError(driver)
    return matrix, np.asarray(singular, dtype=float)


hashes = {
    "negative_source": sha256(NEGATIVE_SOURCE),
    "negative_artifact": sha256(NEGATIVE_ARTIFACT),
    "tangent_archive": sha256(TANGENT_ARCHIVE),
    "primary_artifact": sha256(PRIMARY_ARTIFACT),
}
negative_artifact = json.loads(NEGATIVE_ARTIFACT.read_text())
primary_artifact = json.loads(PRIMARY_ARTIFACT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and negative_artifact["passed"] == negative_artifact["tests"] == 8
    and negative_artifact["outcome"]
    == "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED"
    and primary_artifact["passed"] == primary_artifact["tests"] == 10
    and primary_artifact["outcome"]
    == "NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL"
)
check("the frozen audit inputs retain exact provenance",
      provenance_ok, str(hashes))

print("[setup] replaying the earlier binary negative-fiber construction",
      flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    negative = runpy.run_path(str(NEGATIVE_SOURCE))
replay_ok = bool(
    negative["passed"] == negative["tests"] == 8
    and negative["outcome"] == "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED"
    and sha256(NEGATIVE_ARTIFACT) == EXPECTED_HASHES["negative_artifact"]
    and sum(len(negative["projectors"][time]) for time in TIMES) == 32
)
check("the binary negative-fiber construction replays byte-identically",
      replay_ok)

# Extract fresh orthonormal bases from the already committed binary projector
# midpoints.  The exact primary project's mpmath bases are never loaded.
bases = {}
basis_records = []
bases_ok = True
for time_name in TIMES:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for variant in VARIANTS:
                key = (parity, sector_index, variant)
                projector = np.asarray(
                    negative["projectors"][time_name][key]["projector"],
                    dtype=np.complex128,
                )
                hermitian = (projector + projector.conj().T) / 2
                values, vectors = la.eigh(
                    hermitian, check_finite=True, driver="evr"
                )
                configuration = vectors[:, -15:]
                phase = phase_basis(configuration)
                configuration_defect = operator_norm(
                    configuration.conj().T @ configuration - np.eye(15)
                )
                phase_defect = operator_norm(
                    phase.conj().T @ phase - np.eye(30)
                )
                split_ok = bool(
                    values[14] < 1e-10 and values[15] > 1 - 1e-10
                )
                complete = bool(
                    split_ok
                    and configuration_defect < 1e-12
                    and phase_defect < 1e-12
                    and np.all(np.isfinite(values))
                    and np.all(np.isfinite(phase))
                )
                bases_ok &= complete
                bases[(time_name, parity, sector_index, variant)] = phase
                basis_records.append({
                    "time": time_name,
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "lower_projector_edge": sf(values[14]),
                    "upper_projector_edge": sf(values[15]),
                    "configuration_orthogonality_defect": sf(
                        configuration_defect
                    ),
                    "phase_orthogonality_defect": sf(phase_defect),
                    "split_ok": split_ok,
                    "complete": complete,
                })

check("all 32 binary projectors yield clean rank-15 phase bases",
      bases_ok and len(bases) == len(basis_records) == 32)

cell_records = []
singular_records = []
global_counts = Counter()
spectra_ok = True
synthetic_ok = True
gauge_ok = True
complement_ok = True
driver_ok = True
ordering_ok = True

swap = np.zeros((60, 60), dtype=np.complex128)
swap[:30, 30:] = np.eye(30)
swap[30:, :30] = np.eye(30)

with np.load(TANGENT_ARCHIVE, allow_pickle=False) as tangent_archive:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for variant in VARIANTS:
                key = (parity, sector_index, variant)
                prefix = f"{parity}_sector{sector_index}_t2_{variant}"
                tangent = np.asarray(
                    tangent_archive[f"{prefix}_midpoint"],
                    dtype=np.complex128,
                )
                source = bases[("old",) + key]
                target = bases[("shifted",) + key]
                target_complement = qr_complement(target)
                actual_matrix, actual_singular = leakage(
                    tangent, source, target_complement
                )
                error = audit_floor(tangent, actual_matrix)
                labels = [
                    singular_label(value, error) for value in actual_singular
                ]
                counts = Counter(labels)
                global_counts.update(counts)

                # Known full-intersection and zero-intersection controls.
                full_map = target @ source.conj().T
                _, full_singular = leakage(
                    full_map, source, target_complement
                )
                zero_map = target_complement @ source.conj().T
                _, zero_singular = leakage(
                    zero_map, source, target_complement
                )
                full_labels = [
                    singular_label(value, error) for value in full_singular
                ]
                zero_labels = [
                    singular_label(value, error) for value in zero_singular
                ]
                cell_synthetic = bool(
                    Counter(full_labels)[
                        "AUDIT_SINGULAR_ZERO_CONSISTENT"
                    ] == 30
                    and Counter(zero_labels)[
                        "AUDIT_SINGULAR_NONZERO_RESOLVED"
                    ] == 30
                )
                synthetic_ok &= cell_synthetic

                # Basis-gauge stress.
                phases0 = np.exp(1j * np.arange(30) / 7.0)
                phases1 = np.exp(-1j * np.arange(30) / 11.0)
                gauged_source = source[:, ::-1] * phases0
                gauged_target = target[:, ::-1] * phases1
                gauged_complement = qr_complement(gauged_target)
                _, gauged_singular = leakage(
                    tangent, gauged_source, gauged_complement
                )
                gauge_difference = float(np.max(np.abs(
                    actual_singular - gauged_singular
                )))
                cell_gauge = gauge_difference <= 10 * error
                gauge_ok &= cell_gauge

                # A separately constructed target complement.
                null_complement = la.null_space(
                    target.conj().T, rcond=None, check_finite=True
                )
                _, null_singular = leakage(
                    tangent, source, null_complement
                )
                complement_difference = float(np.max(np.abs(
                    actual_singular - null_singular
                )))
                cell_complement = bool(
                    null_complement.shape == (60, 30)
                    and complement_difference <= 10 * error
                )
                complement_ok &= cell_complement

                # Independent LAPACK SVD driver.
                _, gesvd_singular = leakage(
                    tangent, source, target_complement, driver="gesvd"
                )
                driver_difference = float(np.max(np.abs(
                    actual_singular - gesvd_singular
                )))
                cell_driver = driver_difference <= 10 * error
                driver_ok &= cell_driver

                # Simultaneous q/p ordering convention change.
                swapped_source = swap @ source
                swapped_target = swap @ target
                swapped_tangent = swap @ tangent @ swap
                swapped_complement = qr_complement(swapped_target)
                _, swapped_singular = leakage(
                    swapped_tangent, swapped_source, swapped_complement
                )
                ordering_difference = float(np.max(np.abs(
                    actual_singular - swapped_singular
                )))
                cell_ordering = ordering_difference <= 10 * error
                ordering_ok &= cell_ordering

                finite = bool(
                    np.all(np.isfinite(actual_singular))
                    and math.isfinite(error)
                    and math.isfinite(gauge_difference)
                    and math.isfinite(complement_difference)
                    and math.isfinite(driver_difference)
                    and math.isfinite(ordering_difference)
                )
                ordered = all(
                    actual_singular[index] >= actual_singular[index + 1]
                    for index in range(29)
                )
                dimension_ok = bool(
                    source.shape == target.shape == (60, 30)
                    and target_complement.shape == (60, 30)
                    and actual_matrix.shape == (30, 30)
                    and len(actual_singular) == 30
                )
                spectra_ok &= bool(finite and ordered and dimension_ok)
                resolved_rank = counts[
                    "AUDIT_SINGULAR_NONZERO_RESOLVED"
                ]
                cell_records.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "leakage_norm": sf(actual_singular[0]),
                    "smallest_singular": sf(actual_singular[-1]),
                    "error": sf(error),
                    "smallest_error_units": sf(
                        actual_singular[-1] / error
                    ),
                    "numerical_rank": resolved_rank,
                    "zero_consistent": counts[
                        "AUDIT_SINGULAR_ZERO_CONSISTENT"
                    ],
                    "open": counts["AUDIT_SINGULAR_OPEN"],
                    "full_intersection_control_rank": Counter(full_labels)[
                        "AUDIT_SINGULAR_NONZERO_RESOLVED"
                    ],
                    "zero_intersection_control_rank": Counter(zero_labels)[
                        "AUDIT_SINGULAR_NONZERO_RESOLVED"
                    ],
                    "gauge_spectrum_difference": sf(gauge_difference),
                    "complement_spectrum_difference": sf(
                        complement_difference
                    ),
                    "svd_driver_spectrum_difference": sf(driver_difference),
                    "ordering_spectrum_difference": sf(ordering_difference),
                    "complete": bool(
                        finite and ordered and dimension_ok
                        and cell_synthetic and cell_gauge and cell_complement
                        and cell_driver and cell_ordering
                    ),
                })
                for index, (value, label) in enumerate(zip(
                    actual_singular, labels
                )):
                    singular_records.append({
                        "parity": parity,
                        "sector_index": sector_index,
                        "variant": variant,
                        "index": index,
                        "value": sf(value),
                        "error": sf(error),
                        "error_units": sf(value / error),
                        "label": label,
                    })

check("all 16 independent leakage spectra are finite and complete",
      spectra_ok and len(cell_records) == 16)
check("all full- and zero-intersection controls discriminate",
      synthetic_ok)
check("all spectra survive deterministic basis gauge stress", gauge_ok)
check("complete-QR and null-space complements agree", complement_ok)
check("NumPy and GESVD singular spectra agree", driver_ok)
check("simultaneous q/p ordering changes preserve every spectrum", ordering_ok)
census_ok = bool(
    len(singular_records) == 480
    and sum(global_counts.values()) == 480
)
check("all 480 independent singular values receive frozen labels",
      census_ok, str(dict(global_counts)))

controls_ok = bool(
    provenance_ok and replay_ok and bases_ok and spectra_ok
    and synthetic_ok and gauge_ok and complement_ok and driver_ok
    and ordering_ok and census_ok and len(cell_records) == 16
)
full_rank_cells = sum(
    item["numerical_rank"] == 30 for item in cell_records
)
if not controls_ok:
    outcome = "ADVERSARIAL_NEGATIVE_INTERSECTION_CONTROL_FAILED"
elif full_rank_cells < 16:
    outcome = "ADVERSARIAL_NEGATIVE_INTERSECTION_DISAGREEMENT_OPEN"
else:
    outcome = "ADVERSARIAL_NEGATIVE_INTERSECTION_ZERO_CORROBORATED"

allowed = {
    "ADVERSARIAL_NEGATIVE_INTERSECTION_CONTROL_FAILED",
    "ADVERSARIAL_NEGATIVE_INTERSECTION_DISAGREEMENT_OPEN",
    "ADVERSARIAL_NEGATIVE_INTERSECTION_ZERO_CORROBORATED",
}
check("the preregistered adversarial hierarchy assigns one outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "primary_numeric_singular_values_read": False,
    "fitted_alignment_or_threshold_used": False,
    "basis_count": len(basis_records),
    "basis_records": basis_records,
    "cell_count": len(cell_records),
    "full_rank_cells": full_rank_cells,
    "global_label_counts": dict(global_counts),
    "cell_records": cell_records,
    "singular_records": singular_records,
    "classification": {
        "primary_exact_certificate": "UNDER ADVERSARIAL AUDIT",
        "independent_binary_rank": (
            "STRUCTURAL INDEPENDENT CORROBORATION"
            if outcome
            == "ADVERSARIAL_NEGATIVE_INTERSECTION_ZERO_CORROBORATED"
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
print(f"global singular labels: {dict(global_counts)}")
print(f"full-rank cells: {full_rank_cells}/16")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

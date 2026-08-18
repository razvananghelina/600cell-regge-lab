#!/usr/bin/env python3
"""Mechanically independent audit of generalized phase transport.

Prior-art/independence gate: f071366.
Preregistered audit protocol: 8f0e7f6.
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
DIRECT_SOURCE = (
    HERE / "verify_gravity_600cell_dust_generalized_bundle_direct_precision.py"
)
DIRECT_ARTIFACT = (
    HERE / "gravity_600cell_dust_generalized_bundle_direct_precision.json"
)
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_two_step_full_tangent.py"
TANGENT_METADATA = HERE / "gravity_600cell_dust_two_step_full_tangent.json"
TANGENT_ARCHIVE = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"
PHASE_ARTIFACT = HERE / "gravity_600cell_dust_generalized_phase_transport.json"
OUTPUT = (
    HERE / "gravity_600cell_dust_generalized_phase_transport_adversarial.json"
)

PRIOR_ART_COMMIT = "f071366"
PROTOCOL_COMMIT = "8f0e7f6"
EXPECTED_HASHES = {
    "direct_source": (
        "01479fcaa7e5354ea3bb72306ac8cd433a87b11a539f912075d69273a014b510"
    ),
    "direct_artifact": (
        "8ded406366dbf291da02dfbf995c4e37036cc6ce745d9240d14905664ba6042a"
    ),
    "tangent_source": (
        "c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717"
    ),
    "tangent_metadata": (
        "f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc"
    ),
    "tangent_archive": (
        "ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d"
    ),
    "phase_artifact": (
        "45eb9a3e80ead758d9b3c2f8e1eccff44b06e2759251ab00c447aa53e6705743"
    ),
}

PARITIES = ("even", "odd")
TARGET_SECTORS = (4, 5)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
TIMES = ("old", "shifted")
ZERO_THRESHOLD = 1.0e-10
NONZERO_THRESHOLD = 1.0e-6
IMAGE_RANK_RCOND = 1.0e-12
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
    return float(la.svdvals(matrix)[0]) if matrix.size else 0.0


def mp_to_numpy(matrix):
    return np.array([
        [complex(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


def projector(basis):
    return basis @ basis.conj().T


def phase_basis(basis, convention):
    if convention == "canonical_dual":
        lower = basis.conj()
        upper = basis
    elif convention == "same_representation":
        upper = lower = basis
    elif convention == "swapped_dual":
        upper = basis.conj()
        lower = basis
    else:
        raise ValueError(convention)
    return la.block_diag(upper, lower)


def leakage_label(value):
    if not math.isfinite(value):
        return "INDEPENDENT_OPEN"
    if value <= ZERO_THRESHOLD:
        return "INDEPENDENT_CLOSED"
    if value > NONZERO_THRESHOLD:
        return "INDEPENDENT_NONCLOSING"
    return "INDEPENDENT_OPEN"


def image_metrics(tangent, source, target):
    image = tangent @ source
    singular = la.svdvals(image)
    threshold = IMAGE_RANK_RCOND * singular[0]
    rank = int(np.sum(singular > threshold))
    image_basis, _ = la.qr(image, mode="economic")
    complement = la.null_space(target.conj().T)
    sine_max = operator_norm(complement.conj().T @ image_basis)
    residual = image - target @ (target.conj().T @ image)
    relative_least_squares = operator_norm(residual) / operator_norm(image)
    projector_crosscheck = operator_norm(
        projector(image_basis) - projector(target)
    )
    return {
        "rank": rank,
        "rank_ratio": float(singular[-1] / singular[0]),
        "sine_max": sine_max,
        "relative_least_squares": relative_least_squares,
        "projector_crosscheck": projector_crosscheck,
        "crosscheck_difference": abs(sine_max - projector_crosscheck),
        "target_complement": complement,
    }


hashes = {
    "direct_source": sha256(DIRECT_SOURCE),
    "direct_artifact": sha256(DIRECT_ARTIFACT),
    "tangent_source": sha256(TANGENT_SOURCE),
    "tangent_metadata": sha256(TANGENT_METADATA),
    "tangent_archive": sha256(TANGENT_ARCHIVE),
    "phase_artifact": sha256(PHASE_ARTIFACT),
}
tangent_metadata = json.loads(TANGENT_METADATA.read_text())
phase_result = json.loads(PHASE_ARTIFACT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and tangent_metadata["passed"] == tangent_metadata["tests"] == 16
    and tangent_metadata["outcome"]
    == "TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED"
    and tangent_metadata["numeric_archive_sha256"]
    == EXPECTED_HASHES["tangent_archive"]
    and phase_result["passed"] == phase_result["tests"] == 7
    and phase_result["outcome"] == "GENERALIZED_PHASE_TRANSPORT_REFUTED"
)
check("all frozen independent inputs retain exact provenance",
      provenance_ok, str(hashes))

print("[setup] replaying the earlier direct-precision pencil source", flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    direct = runpy.run_path(str(DIRECT_SOURCE))
direct_replay_ok = bool(
    direct["passed"] == direct["tests"] == 13
    and direct["outcome"] == "DIRECT_GENERALIZED_COMMON_BUNDLE_RESOLVED"
    and sha256(DIRECT_ARTIFACT) == EXPECTED_HASHES["direct_artifact"]
    and len(direct["matrix_cells"]) == 32
    and len(direct["projectors"]) == 32
)
check("the earlier direct-precision source replays byte-identically",
      direct_replay_ok)

# Rebuild the conformal carriers rather than taking the direct verifier's
# stored carrier bases.
conformal = {}
conformal_ok = True
for parity in PARITIES:
    incidence = direct["incidences"][parity]["incidence"].astype(
        np.complex128
    )
    for sector_index in TARGET_SECTORS:
        sector_basis = mp_to_numpy(
            direct["sector_data"][sector_index]["basis"]
        )
        edge_basis = np.kron(
            np.eye(30, dtype=np.complex128), sector_basis
        )
        compressed = edge_basis.conj().T @ incidence
        left, singular, _ = np.linalg.svd(compressed, full_matrices=False)
        basis = left[:, :5]
        conformal[(parity, sector_index)] = basis
        conformal_ok &= bool(
            basis.shape == (30, 5)
            and singular[4] > 1.0e-12 * singular[0]
            and operator_norm(basis.conj().T @ basis - np.eye(5))
            <= ZERO_THRESHOLD
        )

bases = {}
basis_records = []
bases_ok = conformal_ok
for time_name in TIMES:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            carrier = conformal[(parity, sector_index)]
            for variant in VARIANTS:
                key = (time_name, parity, sector_index, variant)
                cell = direct["matrix_cells"][key]
                h_m = (
                    cell["M"]["midpoint"]
                    + cell["M"]["midpoint"].conj().T
                ) / 2
                h_v = (
                    cell["V"]["midpoint"]
                    + cell["V"]["midpoint"].conj().T
                ) / 2
                shape = la.null_space(carrier.conj().T @ h_m)
                m_shape = shape.conj().T @ h_m @ shape
                v_shape = shape.conj().T @ h_v @ shape
                b_matrix = -(m_shape + m_shape.conj().T) / 2
                a_matrix = -(v_shape + v_shape.conj().T) / 2
                b_values = np.linalg.eigvalsh(b_matrix)
                positive_metric = bool(b_values[0] > 0)
                if positive_metric:
                    lower = np.linalg.cholesky(b_matrix)
                    left_reduced = la.solve_triangular(
                        lower, a_matrix, lower=True
                    )
                    whitened = la.solve_triangular(
                        lower, left_reduced.conj().T, lower=True
                    ).conj().T
                    whitened = (whitened + whitened.conj().T) / 2
                    values, whitened_vectors = np.linalg.eigh(whitened)
                    vectors = la.solve_triangular(
                        lower.conj().T, whitened_vectors, lower=False
                    )
                else:
                    values = np.full(25, np.nan)
                    vectors = np.full((25, 25), np.nan, dtype=np.complex128)
                negative = int(np.sum(values < 0))
                positive = int(np.sum(values > 0))
                full_residual = (
                    a_matrix @ vectors
                    - b_matrix @ vectors @ np.diag(values)
                )
                residual_scale = max(
                    1.0,
                    operator_norm(a_matrix @ vectors),
                    operator_norm(b_matrix @ vectors @ np.diag(values)),
                )
                generalized_residual = operator_norm(full_residual) / residual_scale
                lifted = shape @ vectors[:, :15]
                basis, _ = np.linalg.qr(lifted, mode="reduced")
                orthonormality = operator_norm(
                    basis.conj().T @ basis - np.eye(15)
                )
                independent_projector = projector(basis)
                upstream_projector = direct["projectors"][key]["projector"]
                projector_agreement = operator_norm(
                    independent_projector - upstream_projector
                )
                complete = bool(
                    shape.shape == (30, 25)
                    and positive_metric
                    and negative == 15 and positive == 10
                    and orthonormality <= ZERO_THRESHOLD
                    and generalized_residual <= 1.0e-9
                    and projector_agreement <= 1.0e-8
                    and np.all(np.isfinite(basis))
                )
                bases_ok &= complete
                bases[key] = basis
                basis_records.append({
                    "time": time_name,
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "shape_dimension": int(shape.shape[1]),
                    "kinetic_minimum": sf(b_values[0]),
                    "negative": negative,
                    "positive": positive,
                    "generalized_residual": sf(generalized_residual),
                    "orthonormality_defect": sf(orthonormality),
                    "upstream_projector_distance": sf(projector_agreement),
                    "complete": complete,
                })

check("all 32 generalized bases reconstruct by Cholesky whitening",
      bases_ok and len(bases) == len(basis_records) == 32)

canonical_records = []
stress_records = []
label_counts = Counter()
stress_counts = Counter()
tangent_ok = True
image_controls_ok = True
synthetic_controls_ok = True
all_finite = True
with np.load(TANGENT_ARCHIVE) as tangent_archive:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for variant in VARIANTS:
                prefix = f"{parity}_sector{sector_index}_t2_{variant}"
                tangent = tangent_archive[f"{prefix}_midpoint"]
                radii = tangent_archive[f"{prefix}_radii"]
                tangent_ok &= bool(
                    tangent.shape == radii.shape == (60, 60)
                    and np.all(np.isfinite(tangent))
                    and np.all(np.isfinite(radii))
                )
                old_basis = bases[("old", parity, sector_index, variant)]
                shifted_basis = bases[
                    ("shifted", parity, sector_index, variant)
                ]
                source = phase_basis(old_basis, "canonical_dual")
                target = phase_basis(shifted_basis, "canonical_dual")
                metrics = image_metrics(tangent, source, target)
                label = leakage_label(metrics["sine_max"])
                label_counts[label] += 1

                phases0 = np.exp(1j * np.arange(30) / 7.0)
                phases1 = np.exp(-1j * np.arange(30) / 11.0)
                gauged_source = source[:, ::-1] * phases0
                gauged_target = target[:, ::-1] * phases1
                gauged = image_metrics(tangent, gauged_source, gauged_target)
                gauge_difference = abs(
                    metrics["sine_max"] - gauged["sine_max"]
                )

                positive_map = target @ source.conj().T
                positive = image_metrics(positive_map, source, target)
                negative_target = target.copy()
                negative_target[:, 0] = metrics["target_complement"][:, 0]
                negative_map = negative_target @ source.conj().T
                negative = image_metrics(negative_map, source, target)

                image_controls_ok &= bool(
                    metrics["rank"] == 30
                    and metrics["crosscheck_difference"] <= ZERO_THRESHOLD
                    and gauge_difference <= ZERO_THRESHOLD
                )
                synthetic_controls_ok &= bool(
                    positive["rank"] == 30
                    and positive["sine_max"] <= ZERO_THRESHOLD
                    and negative["rank"] == 30
                    and negative["sine_max"] >= 0.5
                )
                finite_values = (
                    metrics["rank_ratio"], metrics["sine_max"],
                    metrics["relative_least_squares"],
                    metrics["projector_crosscheck"], gauge_difference,
                    positive["sine_max"], negative["sine_max"],
                )
                all_finite &= all(math.isfinite(value) for value in finite_values)
                canonical_records.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "archive_radius_frobenius": sf(la.norm(radii, "fro")),
                    "image_rank": metrics["rank"],
                    "image_rank_ratio": sf(metrics["rank_ratio"]),
                    "sine_max": sf(metrics["sine_max"]),
                    "relative_least_squares": sf(
                        metrics["relative_least_squares"]
                    ),
                    "projector_crosscheck": sf(
                        metrics["projector_crosscheck"]
                    ),
                    "crosscheck_difference": sf(
                        metrics["crosscheck_difference"]
                    ),
                    "basis_gauge_difference": sf(gauge_difference),
                    "positive_control_sine": sf(positive["sine_max"]),
                    "negative_control_sine": sf(negative["sine_max"]),
                    "label": label,
                })

                for convention in ("same_representation", "swapped_dual"):
                    stress_source = phase_basis(old_basis, convention)
                    stress_target = phase_basis(shifted_basis, convention)
                    stress = image_metrics(
                        tangent, stress_source, stress_target
                    )
                    stress_label = leakage_label(stress["sine_max"])
                    stress_counts[(convention, stress_label)] += 1
                    all_finite &= bool(
                        math.isfinite(stress["sine_max"])
                        and math.isfinite(stress["relative_least_squares"])
                    )
                    stress_records.append({
                        "parity": parity,
                        "sector_index": sector_index,
                        "variant": variant,
                        "convention": convention,
                        "image_rank": stress["rank"],
                        "sine_max": sf(stress["sine_max"]),
                        "relative_least_squares": sf(
                            stress["relative_least_squares"]
                        ),
                        "label": stress_label,
                    })

check("all 16 archived second-slab tangents are finite and complete",
      tangent_ok and len(canonical_records) == 16)
check("all canonical image, crosscheck and basis-gauge controls pass",
      image_controls_ok and all_finite)
check("all 16 positive and negative synthetic controls discriminate",
      synthetic_controls_ok)
check("all 16 canonical and 32 convention-stress records exist",
      len(canonical_records) == 16 and len(stress_records) == 32,
      f"canonical={dict(label_counts)}, stress={dict(stress_counts)}")

controls_ok = bool(
    provenance_ok and direct_replay_ok and bases_ok and tangent_ok
    and image_controls_ok and synthetic_controls_ok and all_finite
    and len(canonical_records) == 16 and len(stress_records) == 32
)
if not controls_ok:
    outcome = "ADVERSARIAL_PHASE_TRANSPORT_CONTROL_FAILED"
elif label_counts["INDEPENDENT_CLOSED"]:
    outcome = "ADVERSARIAL_PHASE_TRANSPORT_CONTRADICTED"
elif label_counts["INDEPENDENT_OPEN"]:
    outcome = "ADVERSARIAL_PHASE_TRANSPORT_DISAGREEMENT_OPEN"
elif label_counts["INDEPENDENT_NONCLOSING"] == 16:
    outcome = "ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED"
else:
    outcome = "ADVERSARIAL_PHASE_TRANSPORT_CONTROL_FAILED"

allowed = {
    "ADVERSARIAL_PHASE_TRANSPORT_CONTROL_FAILED",
    "ADVERSARIAL_PHASE_TRANSPORT_CONTRADICTED",
    "ADVERSARIAL_PHASE_TRANSPORT_DISAGREEMENT_OPEN",
    "ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED",
}
check("the preregistered adversarial hierarchy assigns one outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "accepted_result_numeric_values_read": False,
    "fitted_alignment_used": False,
    "thresholds": {
        "zero": sf(ZERO_THRESHOLD),
        "nonzero": sf(NONZERO_THRESHOLD),
        "image_rank_rcond": sf(IMAGE_RANK_RCOND),
    },
    "basis_records": basis_records,
    "canonical_label_counts": dict(label_counts),
    "canonical_records": canonical_records,
    "convention_stress_counts": {
        f"{convention}:{label}": count
        for (convention, label), count in sorted(stress_counts.items())
    },
    "convention_stress_records": stress_records,
    "classification": {
        "high_precision_refutation": (
            "DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED"
            if outcome
            == "ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED"
            else "OPEN PENDING DISAGREEMENT RESOLUTION"
        ),
        "this_replication": "STRUCTURAL INDEPENDENT CORROBORATION",
        "transported_intersection_or_graph": "NOT COMPUTED",
        "propagator_dispersion_mass_or_speed": "NOT COMPUTED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"canonical labels: {dict(label_counts)}")
print(f"convention stress: {dict(stress_counts)}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)


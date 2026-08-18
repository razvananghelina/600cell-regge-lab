#!/usr/bin/env python3
"""Certified rank census for the transported negative-stiffness phase fiber.

Prior-art/framing gate: 65419f6.
Preregistered protocol: 2927fb9.
"""

from collections import Counter
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy

import mpmath as mp
from flint import ctx


HERE = Path(__file__).resolve().parent
PHASE_SOURCE = (
    HERE / "verify_gravity_600cell_dust_generalized_phase_transport.py"
)
PHASE_ARTIFACT = HERE / "gravity_600cell_dust_generalized_phase_transport.json"
PHASE_ADVERSARIAL_ARTIFACT = (
    HERE / "gravity_600cell_dust_generalized_phase_transport_adversarial.json"
)
NEGATIVE_SOURCE = HERE / "verify_gravity_600cell_dust_negative_fiber_transport.py"
NEGATIVE_ARTIFACT = HERE / "gravity_600cell_dust_negative_fiber_transport.json"
OUTPUT = HERE / "gravity_600cell_dust_negative_transported_intersection.json"

PRIOR_ART_COMMIT = "65419f6"
PROTOCOL_COMMIT = "2927fb9"
EXPECTED_HASHES = {
    "phase_source": (
        "9c4c36b463a8faaa8d40b7db1b6b1852e3c04155c1b6ada4d02fbda747f6fcf3"
    ),
    "phase_artifact": (
        "45eb9a3e80ead758d9b3c2f8e1eccff44b06e2759251ab00c447aa53e6705743"
    ),
    "phase_adversarial_artifact": (
        "c33615ac6d0f3133e53077f46c5ee766b9c633d4d64c32124c24839c9c84c880"
    ),
    "negative_source": (
        "f462e507500d7f02ecf799f0d4b320e05795216a36a0d10eb908d6dc67b48181"
    ),
    "negative_artifact": (
        "d630bf07066f88c35eee5a62a80ec1f43399a95ea882a43528289220c67f4599"
    ),
}

DPS = 100
BALL_DPS = 80
MP_FLOOR = mp.mpf("1e-70")
CONTROL_FACTOR = mp.mpf("1e-60")
PARITIES = ("even", "odd")
TIMES = ("old", "shifted")
TARGET_SECTORS = (4, 5)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)

mp.mp.dps = DPS
ctx.dps = BALL_DPS
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


def smp(value, digits=30):
    value = mp.mpf(value)
    if mp.isinf(value):
        return "inf" if value > 0 else "-inf"
    if mp.isnan(value):
        return "nan"
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def numpy_to_mp(matrix):
    rows, columns = matrix.shape
    result = mp.matrix(rows, columns)
    for row in range(rows):
        for column in range(columns):
            value = complex(matrix[row, column])
            result[row, column] = mp.mpc(repr(value.real), repr(value.imag))
    return result


def singular_label(value, error):
    if not mp.isfinite(value) or not mp.isfinite(error):
        return "SINGULAR_OPEN"
    if value <= 10 * error:
        return "SINGULAR_ZERO_CONSISTENT"
    if value > 100 * error:
        return "SINGULAR_NONZERO_RESOLVED"
    return "SINGULAR_OPEN"


def singular_values(matrix):
    values = mp.svd(matrix, compute_uv=False)
    return sorted(
        (mp.mpf(mp.re(values[index])) for index in range(values.rows)),
        reverse=True,
    )


hashes = {
    "phase_source": sha256(PHASE_SOURCE),
    "phase_artifact": sha256(PHASE_ARTIFACT),
    "phase_adversarial_artifact": sha256(PHASE_ADVERSARIAL_ARTIFACT),
    "negative_source": sha256(NEGATIVE_SOURCE),
    "negative_artifact": sha256(NEGATIVE_ARTIFACT),
}
phase_artifact = json.loads(PHASE_ARTIFACT.read_text())
phase_adversarial_artifact = json.loads(PHASE_ADVERSARIAL_ARTIFACT.read_text())
negative_artifact = json.loads(NEGATIVE_ARTIFACT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and phase_artifact["passed"] == phase_artifact["tests"] == 7
    and phase_artifact["outcome"] == "GENERALIZED_PHASE_TRANSPORT_REFUTED"
    and phase_adversarial_artifact["passed"]
    == phase_adversarial_artifact["tests"] == 8
    and phase_adversarial_artifact["outcome"]
    == "ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED"
    and negative_artifact["passed"] == negative_artifact["tests"] == 8
    and negative_artifact["outcome"]
    == "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED"
)
check("all frozen source and control artifacts retain exact provenance",
      provenance_ok, str(hashes))

print("[setup] replaying the exact source-ball and tangent chain", flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    phase = runpy.run_path(str(PHASE_SOURCE))
residual = phase["residual"]
phase_replay_ok = bool(
    phase["passed"] == phase["tests"] == 7
    and phase["outcome"] == "GENERALIZED_PHASE_TRANSPORT_REFUTED"
    and sha256(PHASE_ARTIFACT) == EXPECTED_HASHES["phase_artifact"]
    and len(residual["exact_cells"]) == 32
    and len(residual["conformal"]) == 4
    and len(phase["tangent_cells"]) == 16
)
check("the accepted high-precision source and tangent chain replays",
      phase_replay_ok)

print("[setup] replaying binary negative projectors as overlap controls",
      flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    negative = runpy.run_path(str(NEGATIVE_SOURCE))
mp.mp.dps = DPS
ctx.dps = BALL_DPS
negative_replay_ok = bool(
    negative["passed"] == negative["tests"] == 8
    and negative["outcome"] == "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED"
    and sha256(NEGATIVE_ARTIFACT) == EXPECTED_HASHES["negative_artifact"]
    and sum(len(negative["projectors"][time]) for time in TIMES) == 32
)
check("the binary negative-projector controls replay byte-identically",
      negative_replay_ok)

mp_operator_norm = residual["mp_operator_norm"]
mp_hermitian_norm = residual["mp_hermitian_norm"]
acb_midpoint_to_mp = residual["acb_midpoint_to_mp"]
acb_radius_frobenius = residual["acb_radius_frobenius"]
shape_basis = residual["shape_basis"]
block_diagonal_phase = phase["block_diagonal_phase"]

# Exact synthetic controls for the rank/nullity interpretation.
q_first = mp.matrix(60, 60)
q_last = mp.matrix(60, 60)
for index in range(30):
    q_first[index, index] = 1
    q_last[30 + index, 30 + index] = 1
identity = mp.eye(60)
control_error = mp.mpf("1e-50")
full_control_values = singular_values((identity - q_first) * identity * q_first)
zero_control_values = singular_values((identity - q_last) * identity * q_first)
full_control_labels = [
    singular_label(value, control_error) for value in full_control_values
]
zero_control_labels = [
    singular_label(value, control_error) for value in zero_control_values
]
synthetic_controls_ok = bool(
    Counter(full_control_labels)["SINGULAR_ZERO_CONSISTENT"] == 60
    and Counter(zero_control_labels)["SINGULAR_NONZERO_RESOLVED"] == 30
    and Counter(zero_control_labels)["SINGULAR_ZERO_CONSISTENT"] == 30
)
check("known full- and zero-intersection synthetic controls classify exactly",
      synthetic_controls_ok)

# Reconstruct the ordinary Euclidean negative-stiffness projectors directly
# from source M,V balls.  No generalized-pencil projector enters this loop.
projectors = {}
projector_records = []
projector_ok = True
binary_overlap_ok = True
for time_name in TIMES:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            u_basis = residual["conformal"][(parity, sector_index)]
            for variant in VARIANTS:
                key = (time_name, parity, sector_index, variant)
                balls = residual["exact_cells"][key]
                m_matrix = acb_midpoint_to_mp(balls["M"])
                v_matrix = acb_midpoint_to_mp(balls["V"])
                m_matrix = (m_matrix + m_matrix.H) / 2
                v_matrix = (v_matrix + v_matrix.H) / 2
                m_norm = mp_hermitian_norm(m_matrix)
                v_norm = mp_hermitian_norm(v_matrix)
                epsilon_m = (
                    acb_radius_frobenius(balls["M"])
                    + MP_FLOOR * max(mp.mpf(1), m_norm)
                )
                epsilon_v = (
                    acb_radius_frobenius(balls["V"])
                    + MP_FLOOR * max(mp.mpf(1), v_norm)
                )
                shape = shape_basis(u_basis, m_matrix)
                if shape is None:
                    projector_ok = False
                    continue
                w_basis = shape["basis"]
                row_gap = shape["row_gap"]
                eta_s = (
                    2 * epsilon_m / (row_gap - 2 * epsilon_m) + MP_FLOOR
                    if row_gap > 2 * epsilon_m else mp.inf
                )
                epsilon_a = (
                    epsilon_v
                    + 2 * eta_s * (v_norm + epsilon_v)
                    + MP_FLOOR * max(mp.mpf(1), v_norm)
                )
                a_matrix = -(w_basis.H * v_matrix * w_basis)
                a_matrix = (a_matrix + a_matrix.H) / 2
                values, vectors = mp.eighe(a_matrix)
                midpoint_negative = sum(
                    mp.re(values[index]) < 0 for index in range(25)
                )
                midpoint_positive = sum(
                    mp.re(values[index]) > 0 for index in range(25)
                )
                negative_edge = mp.re(values[14])
                positive_edge = mp.re(values[15])
                gap = positive_edge - negative_edge
                sign_resolved = bool(
                    negative_edge < -100 * epsilon_a
                    and positive_edge > 100 * epsilon_a
                )
                q_minus = vectors[:, :15]
                q_plus = vectors[:, 15:]
                cross_residual = mp_operator_norm(
                    q_plus.H * a_matrix * q_minus
                )
                eta_residual = (
                    mp.sin(mp.mpf("0.5") * mp.atan(2 * cross_residual / gap))
                    if gap > 0 else mp.inf
                )
                eta_source = (
                    2 * epsilon_a / (gap - 2 * epsilon_a) + MP_FLOOR
                    if gap > 2 * epsilon_a else mp.inf
                )
                lifted = w_basis * q_minus
                projector = lifted * lifted.H
                projector = (projector + projector.H) / 2
                eta_p = 2 * eta_s + eta_residual + eta_source + MP_FLOOR
                hermiticity = mp_operator_norm(projector - projector.H)
                idempotence = mp_operator_norm(projector * projector - projector)
                relative_null = shape["null_residual"] / max(
                    mp.mpf(1), mp_operator_norm(u_basis.H * m_matrix)
                )
                binary = negative["projectors"][time_name][
                    (parity, sector_index, variant)
                ]
                binary_distance = mp_operator_norm(
                    projector - numpy_to_mp(binary["projector"])
                )
                binary_overlap = bool(binary_distance <= mp.mpf(binary["eta"]))
                complete = bool(
                    len(shape["pivots"]) == 5
                    and shape["orthogonality"] < CONTROL_FACTOR
                    and relative_null < CONTROL_FACTOR
                    and midpoint_negative == 15
                    and midpoint_positive == 10
                    and sign_resolved
                    and gap > 2 * epsilon_a
                    and hermiticity < CONTROL_FACTOR
                    and idempotence < CONTROL_FACTOR
                    and mp.isfinite(eta_p)
                )
                projector_ok &= complete
                binary_overlap_ok &= binary_overlap
                projectors[key] = {"projector": projector, "eta": eta_p}
                projector_records.append({
                    "time": time_name,
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "shape_pivots": list(shape["pivots"]),
                    "shape_row_gap": smp(row_gap),
                    "shape_relative_null_residual": smp(relative_null),
                    "shape_orthogonality": smp(shape["orthogonality"]),
                    "source_M_error": smp(epsilon_m),
                    "source_V_error": smp(epsilon_v),
                    "restricted_stiffness_error": smp(epsilon_a),
                    "negative_edge": smp(negative_edge),
                    "positive_edge": smp(positive_edge),
                    "sign_error_units": smp(
                        min(abs(negative_edge), positive_edge) / epsilon_a
                    ),
                    "cluster_gap": smp(gap),
                    "eigen_cross_residual": smp(cross_residual),
                    "projector_shape_term": smp(2 * eta_s),
                    "projector_residual_term": smp(eta_residual),
                    "projector_source_term": smp(eta_source),
                    "projector_error": smp(eta_p),
                    "projector_hermiticity": smp(hermiticity),
                    "projector_idempotence": smp(idempotence),
                    "midpoint_negative": midpoint_negative,
                    "midpoint_positive": midpoint_positive,
                    "sign_resolved": sign_resolved,
                    "binary_projector_distance": smp(binary_distance),
                    "binary_projector_error": smp(binary["eta"]),
                    "binary_overlap": binary_overlap,
                    "complete": complete,
                })

check("all 32 high-precision negative projectors certify",
      projector_ok and len(projectors) == len(projector_records) == 32)
check("all high-precision projectors overlap their binary controls",
      binary_overlap_ok)

cell_records = []
singular_records = []
global_counts = Counter()
spectra_ok = True
structural_ok = True
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for variant in VARIANTS:
            old_key = ("old", parity, sector_index, variant)
            shifted_key = ("shifted", parity, sector_index, variant)
            tangent_key = (parity, sector_index, variant)
            if old_key not in projectors or shifted_key not in projectors:
                spectra_ok = False
                continue
            old = projectors[old_key]
            shifted = projectors[shifted_key]
            q0 = block_diagonal_phase(old["projector"])
            q1 = block_diagonal_phase(shifted["projector"])
            tangent_ball = phase["tangent_cells"][tangent_key]
            tangent = acb_midpoint_to_mp(tangent_ball)
            tangent_norm = mp_operator_norm(tangent)
            epsilon_t = (
                acb_radius_frobenius(tangent_ball)
                + MP_FLOOR * max(mp.mpf(1), tangent_norm)
            )
            residual_matrix = (mp.eye(60) - q1) * tangent * q0
            residual_error = (
                epsilon_t
                + (old["eta"] + shifted["eta"]
                   + old["eta"] * shifted["eta"])
                * (tangent_norm + epsilon_t)
                + MP_FLOOR * max(mp.mpf(1), tangent_norm)
            )
            singulars = singular_values(residual_matrix)
            labels = [
                singular_label(value, residual_error) for value in singulars
            ]
            counts = Counter(labels)
            global_counts.update(counts)
            finite = all(mp.isfinite(value) for value in singulars)
            ordered = all(
                singulars[index] >= singulars[index + 1]
                for index in range(59)
            )
            resolved_rank = counts["SINGULAR_NONZERO_RESOLVED"]
            lower_half_zero = all(
                label == "SINGULAR_ZERO_CONSISTENT"
                for label in labels[30:]
            )
            cell_structural = bool(
                resolved_rank <= 30 and lower_half_zero
            )
            structural_ok &= cell_structural
            spectra_ok &= bool(
                len(singulars) == 60 and finite and ordered
            )
            certified_dimension = 0 if resolved_rank == 30 else None
            cell_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "tangent_norm": smp(tangent_norm),
                "tangent_error": smp(epsilon_t),
                "residual_norm": smp(singulars[0]),
                "residual_error": smp(residual_error),
                "resolved_rank_lower_bound": resolved_rank,
                "structural_rank_upper_bound": 30,
                "certified_intersection_dimension": certified_dimension,
                "zero_consistent": counts["SINGULAR_ZERO_CONSISTENT"],
                "open": counts["SINGULAR_OPEN"],
                "nonzero_resolved": resolved_rank,
                "lower_30_zero_consistent": lower_half_zero,
                "complete": bool(finite and ordered and cell_structural),
            })
            for index, (value, label) in enumerate(zip(singulars, labels)):
                singular_records.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "index": index,
                    "value": smp(value),
                    "error": smp(residual_error),
                    "error_units": smp(value / residual_error),
                    "label": label,
                })

check("all 16 high-precision singular spectra are finite and ordered",
      spectra_ok and len(cell_records) == 16)
census_ok = bool(
    len(singular_records) == 960
    and sum(global_counts.values()) == 960
)
check("all 960 singular values receive frozen labels",
      census_ok, str(dict(global_counts)))
check("all spectra obey the exact rank-30 structural upper bound",
      structural_ok)

controls_ok = bool(
    provenance_ok and phase_replay_ok and negative_replay_ok
    and synthetic_controls_ok and projector_ok and binary_overlap_ok
    and spectra_ok and census_ok and structural_ok
    and len(projectors) == 32 and len(cell_records) == 16
)
zero_cells = sum(
    item["certified_intersection_dimension"] == 0
    for item in cell_records
)
if not controls_ok:
    outcome = "NEGATIVE_TRANSPORTED_INTERSECTION_CONTROL_FAILED"
elif zero_cells == 16:
    outcome = "NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL"
else:
    outcome = "NEGATIVE_TRANSPORTED_INTERSECTION_DIMENSION_OPEN"

allowed = {
    "NEGATIVE_TRANSPORTED_INTERSECTION_CONTROL_FAILED",
    "NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL",
    "NEGATIVE_TRANSPORTED_INTERSECTION_DIMENSION_OPEN",
}
check("the preregistered negative-intersection hierarchy assigns one outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "acceptance_status": "AWAITING_ADVERSARIAL_REPLICATION",
    "desired_rank_inspected_before_protocol": False,
    "fitted_metric_alignment_or_graph_used": False,
    "synthetic_controls": {
        "full_intersection": dict(Counter(full_control_labels)),
        "zero_intersection": dict(Counter(zero_control_labels)),
        "passed": synthetic_controls_ok,
    },
    "projector_count": len(projector_records),
    "projector_records": projector_records,
    "cell_count": len(cell_records),
    "zero_intersection_cells": zero_cells,
    "global_label_counts": dict(global_counts),
    "cell_records": cell_records,
    "singular_records": singular_records,
    "classification": {
        "negative_stiffness_inertia": "DERIVED COMPUTATIONAL",
        "negative_spectral_fiber": (
            "STRUCTURAL RELATIVE TO FROZEN CARRIER HILBERT METRIC"
        ),
        "transported_intersection_primary_result": (
            "PROVISIONAL DERIVED COMPUTATIONAL ZERO"
            if outcome
            == "NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL"
            else "OPEN"
        ),
        "adversarial_acceptance": "NOT YET PERFORMED",
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
print(f"zero-intersection cells: {zero_cells}/16")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

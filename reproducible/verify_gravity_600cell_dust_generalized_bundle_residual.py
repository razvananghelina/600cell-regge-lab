#!/usr/bin/env python3
"""Residual-certified finite-family generalized bundle comparison.

Prior-art/framing commit: 7e815de.
Protocol commits: 1b48e52, ad2eadc.
"""

from collections import Counter
from itertools import combinations
import contextlib
import hashlib
import io
import json
import math
from pathlib import Path
import runpy

import mpmath as mp
import numpy as np
import scipy.linalg as la
from flint import ctx


HERE = Path(__file__).resolve().parent
DIRECT_SOURCE = (
    HERE / "verify_gravity_600cell_dust_generalized_bundle_direct_precision.py"
)
DIRECT_ARTIFACT = (
    HERE / "gravity_600cell_dust_generalized_bundle_direct_precision.json"
)
OUTPUT = HERE / "gravity_600cell_dust_generalized_bundle_residual.json"

PRIOR_ART_COMMIT = "7e815de"
PROTOCOL_COMMITS = ("1b48e52", "ad2eadc")
EXPECTED_HASHES = {
    "direct_source": (
        "01479fcaa7e5354ea3bb72306ac8cd433a87b11a539f912075d69273a014b510"
    ),
    "direct_artifact": (
        "8ded406366dbf291da02dfbf995c4e37036cc6ce745d9240d14905664ba6042a"
    ),
}

DPS = 100
BALL_DPS = 80
MP_FLOOR = mp.mpf("1e-70")
RANK_FACTOR = mp.mpf("1e-50")
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


def mp_diag(values):
    result = mp.matrix(len(values), len(values))
    for index, value in enumerate(values):
        result[index, index] = value
    return result


def mp_hermitian_norm(matrix):
    if matrix.rows == 0 or matrix.cols == 0:
        return mp.mpf(0)
    hermitian = (matrix + matrix.H) / 2
    values = mp.eighe(hermitian, eigvals_only=True)
    return max(abs(values[index]) for index in range(values.rows))


def mp_operator_norm(matrix):
    if matrix.rows == 0 or matrix.cols == 0:
        return mp.mpf(0)
    gram = (matrix.H * matrix + (matrix.H * matrix).H) / 2
    values = mp.eighe(gram, eigvals_only=True)
    largest = max(mp.re(values[index]) for index in range(values.rows))
    return mp.sqrt(max(mp.mpf(0), largest))


def mp_to_numpy(matrix):
    return np.array([
        [complex(float(mp.re(matrix[row, column])),
                 float(mp.im(matrix[row, column])))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


def acb_midpoint_to_mp(matrix):
    result = mp.matrix(matrix.nrows(), matrix.ncols())
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            value = matrix[row, column]
            real = value.real.mid().str(BALL_DPS, radius=False, more=True)
            imag = value.imag.mid().str(BALL_DPS, radius=False, more=True)
            result[row, column] = mp.mpc(real, imag)
    return result


def acb_radius_frobenius(matrix):
    # A factor two safely dominates the final binary conversion of each
    # nonnegative radius upper bound.  MP_FLOOR dominates in the present cells.
    squares = []
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            value = matrix[row, column]
            real = mp.mpf(repr(float(value.real.rad().upper())))
            imag = mp.mpf(repr(float(value.imag.rad().upper())))
            radius = 2 * mp.sqrt(real**2 + imag**2)
            squares.append(radius**2)
    return mp.sqrt(mp.fsum(squares))


def lexicographic_orthonormal_columns(matrix, target_rank):
    scale = max(mp.mpf(1), mp_operator_norm(matrix))
    threshold = RANK_FACTOR * scale
    vectors = []
    pivots = []
    for column in range(matrix.cols):
        vector = mp.matrix(matrix.rows, 1)
        for row in range(matrix.rows):
            vector[row] = matrix[row, column]
        for _ in range(2):
            for basis in vectors:
                vector -= basis * (basis.H * vector)[0]
        norm = mp.sqrt(max(mp.mpf(0), mp.re((vector.H * vector)[0])))
        if norm > threshold:
            vectors.append(vector / norm)
            pivots.append(column)
            if len(vectors) == target_rank:
                break
    if len(vectors) != target_rank:
        return mp.matrix(matrix.rows, 0), tuple(pivots), mp.inf, mp.inf
    basis = mp.matrix(matrix.rows, target_rank)
    for column, vector in enumerate(vectors):
        for row in range(matrix.rows):
            basis[row, column] = vector[row]
    orthogonality = mp_operator_norm(basis.H * basis - mp.eye(target_rank))
    residual = mp_operator_norm((mp.eye(matrix.rows) - basis * basis.H) * matrix)
    return basis, tuple(pivots), orthogonality, residual / scale


def compressed_incidence(prior, parity, sector_index):
    incidence = prior["incidences"][parity]["incidence"]
    sector = prior["sector_data"][sector_index]
    basis = sector["basis"]
    if sector["dimension"] != 1 or basis.rows != 24 or basis.cols != 1:
        raise RuntimeError("target sector is not a 24x1 minimal carrier")
    compressed = mp.matrix(30, 120)
    for orbit in range(30):
        for group_element in range(24):
            coefficient = mp.conj(basis[group_element, 0])
            edge_row = 24 * orbit + group_element
            for vertex in np.flatnonzero(incidence[edge_row]):
                compressed[orbit, int(vertex)] += coefficient
    return compressed


def shape_basis(u_basis, m_matrix):
    row = u_basis.H * m_matrix
    _, pivots, _, _ = lexicographic_orthonormal_columns(row, 5)
    if len(pivots) != 5:
        return None
    free = tuple(index for index in range(30) if index not in pivots)
    pivot_block = mp.matrix([[row[r, c] for c in pivots] for r in range(5)])
    free_block = mp.matrix([[row[r, c] for c in free] for r in range(5)])
    coefficients = -mp.inverse(pivot_block) * free_block
    raw = mp.matrix(30, 25)
    for local_row, global_row in enumerate(pivots):
        for column in range(25):
            raw[global_row, column] = coefficients[local_row, column]
    for column, global_row in enumerate(free):
        raw[global_row, column] = 1
    gram = (raw.H * raw + (raw.H * raw).H) / 2
    gram_values, gram_vectors = mp.eighe(gram)
    if min(gram_values[index] for index in range(gram_values.rows)) <= 0:
        return None
    inverse_sqrt = gram_vectors * mp_diag([
        1 / mp.sqrt(gram_values[index]) for index in range(gram_values.rows)
    ]) * gram_vectors.H
    basis = raw * inverse_sqrt
    row_gram = (row * row.H + (row * row.H).H) / 2
    row_values = mp.eighe(row_gram, eigvals_only=True)
    row_gap = mp.sqrt(min(mp.re(row_values[index])
                          for index in range(row_values.rows)))
    return {
        "basis": basis,
        "pivots": pivots,
        "row_gap": row_gap,
        "null_residual": mp_operator_norm(row * basis),
        "orthogonality": mp_operator_norm(basis.H * basis - mp.eye(25)),
    }


def projector_label(distance, error):
    if not mp.isfinite(distance) or not mp.isfinite(error):
        return "OPEN"
    if distance <= 10 * error:
        return "ZERO_CONSISTENT"
    if distance > 100 * error:
        return "ROTATION_RESOLVED"
    return "OPEN"


def projector_distance(left, right):
    return mp_hermitian_norm(left - right)


hashes = {
    "direct_source": sha256(DIRECT_SOURCE),
    "direct_artifact": sha256(DIRECT_ARTIFACT),
}
provenance_ok = hashes == EXPECTED_HASHES
check("the committed direct-bundle inputs retain exact provenance",
      provenance_ok, str(hashes))

print("[setup] replaying the committed direct verifier", flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    prior = runpy.run_path(str(DIRECT_SOURCE))
prior_replay_ok = bool(
    prior["passed"] == prior["tests"] == 13
    and prior["outcome"] == "DIRECT_GENERALIZED_COMMON_BUNDLE_RESOLVED"
    and sha256(DIRECT_ARTIFACT) == EXPECTED_HASHES["direct_artifact"]
)
check("the committed direct verifier replays byte-identically",
      prior_replay_ok)

# Replay only the exact source-ball construction, retaining the acb matrices
# that the previous verifier intentionally reduced to binary64 for its audit.
print("[setup] reconstructing schedule-local Flint source balls", flush=True)
exact_cells = {}
reconstruction_ok = True
for parity in PARITIES:
    a1, r1 = [mp.mpf(value) for value in
              prior["first_tick"]["solutions"][parity]["state"]]
    a2, r2 = [mp.mpf(value) for value in
              prior["second_tick"]["solutions"][parity]["state_absolute"]]
    a3, r3 = [mp.mpf(value) for value in
              prior["third_tick"]["solutions"][parity]["state_absolute"]]
    slab_specs = (
        ("slab1", mp.mpf(0), a1, r1),
        ("slab2", a1, a2, r2),
        ("slab3", a2, a3, r3),
    )
    slab_data = {}
    for slab_name, a_old, a_new, r_value in slab_specs:
        index_data, kind_values = prior["slab_index_data"](
            prior["models"][parity], a_old, a_new, r_value
        )
        geometry = prior["prepare_geometry"](
            prior["models"][parity], index_data
        )
        mapping = prior["boundary_mapping"](index_data)
        sectors, _ = prior["high_precision_sector_bases"](index_data)
        pattern_cache, branch = prior["high_precision_pattern_cache"](
            geometry["patterns"], kind_values
        )
        kernels, kernel_control = prior["assemble_full_representative_kernels"](
            index_data, geometry, pattern_cache
        )
        reconstruction_ok &= bool(
            branch["entry_pass"]
            and kernel_control["maximum_imaginary"] < prior["ARITHMETIC_FLOOR"]
        )
        slab_data[slab_name] = {
            "mapping": mapping,
            "sectors": sectors,
            "kernels": kernels,
        }

    for sector_index in TARGET_SECTORS:
        projected = {name: {} for name in ("slab1", "slab2", "slab3")}
        for slab_name in projected:
            sector = slab_data[slab_name]["sectors"][sector_index]
            reconstruction_ok &= sector["dimension"] == 1
            for variant in VARIANTS:
                raw = prior["project_full_kernel"](
                    slab_data[slab_name]["kernels"][variant], sector
                )
                block = (raw + raw.H) / 2
                _, determinant, tangent, _ = prior["build_tangent_ball"](
                    block, 1, slab_data[slab_name]["mapping"]
                )
                reconstruction_ok &= not determinant.contains(0)
                principal, residuals = prior["reconstruct_principal"](
                    prior["split_tangent"](tangent, 30)
                )
                reconstruction_ok &= all(
                    residual[row, column].contains(0)
                    for residual in residuals.values()
                    for row in range(residual.nrows())
                    for column in range(residual.ncols())
                )
                projected[slab_name][variant] = principal

        for time_name, lower, upper in (
            ("old", "slab1", "slab2"),
            ("shifted", "slab2", "slab3"),
        ):
            for variant in VARIANTS:
                lower_p = projected[lower][variant]
                upper_p = projected[upper][variant]
                kminus = lower_p["10"]
                kzero = lower_p["11"] + upper_p["00"]
                kplus = upper_p["01"]
                m_ball = prior["scale"](kminus + kplus, "0.5")
                v_ball = kminus + kzero + kplus
                exact_cells[(time_name, parity, sector_index, variant)] = {
                    "M": m_ball,
                    "V": v_ball,
                }

check("all 32 schedule-local Flint source pairs reconstruct",
      reconstruction_ok and len(exact_cells) == 32)

# Geometry-only conformal bases.
conformal = {}
conformal_records = []
conformal_ok = True
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        compressed = compressed_incidence(prior, parity, sector_index)
        basis, pivots, orthogonality, residual = (
            lexicographic_orthonormal_columns(compressed, 5)
        )
        binary_basis = prior["conformal"][(parity, sector_index)]["basis"]
        binary_projector = binary_basis @ binary_basis.conj().T
        high_projector = mp_to_numpy(basis * basis.H)
        binary_distance = la.svdvals(high_projector - binary_projector)[0]
        complete = bool(
            basis.cols == 5 and len(pivots) == 5
            and orthogonality < CONTROL_FACTOR
            and residual < CONTROL_FACTOR
            and binary_distance <= 10 * prior["conformal"][
                (parity, sector_index)
            ]["eta"]
        )
        conformal_ok &= complete
        conformal[(parity, sector_index)] = basis
        conformal_records.append({
            "parity": parity,
            "sector_index": sector_index,
            "pivots": list(pivots),
            "orthogonality": smp(orthogonality),
            "relative_image_residual": smp(residual),
            "binary_projector_distance": f"{binary_distance:.17e}",
            "complete": complete,
        })

check("all four lexicographic conformal carriers certify rank five",
      conformal_ok and len(conformal_records) == 4)

# High-precision schedule-local projectors.
projectors = {}
projector_records = []
projector_ok = True
binary_overlap_ok = True
for time_name in TIMES:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            u_basis = conformal[(parity, sector_index)]
            for variant in VARIANTS:
                key = (time_name, parity, sector_index, variant)
                balls = exact_cells[key]
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
                epsilon_b = (
                    epsilon_m + 2 * eta_s * (m_norm + epsilon_m)
                    + MP_FLOOR * max(mp.mpf(1), m_norm)
                )
                epsilon_a = (
                    epsilon_v + 2 * eta_s * (v_norm + epsilon_v)
                    + MP_FLOOR * max(mp.mpf(1), v_norm)
                )
                b_matrix = -(w_basis.H * m_matrix * w_basis)
                a_matrix = -(w_basis.H * v_matrix * w_basis)
                b_matrix = (b_matrix + b_matrix.H) / 2
                a_matrix = (a_matrix + a_matrix.H) / 2
                b_values = mp.eighe(b_matrix, eigvals_only=True)
                b_min = min(mp.re(b_values[i]) for i in range(25))
                b_max = max(mp.re(b_values[i]) for i in range(25))
                b_lower = b_min - epsilon_b
                positive = b_lower > 0 and b_min > 100 * epsilon_b
                if not positive:
                    projector_ok = False
                    continue
                lower = mp.cholesky(b_matrix)
                lower_inverse = mp.inverse(lower)
                h_matrix = lower_inverse * a_matrix * lower_inverse.H
                h_matrix = (h_matrix + h_matrix.H) / 2
                values, vectors = mp.eighe(h_matrix)
                negative = sum(mp.re(values[i]) < 0 for i in range(25))
                positive_count = sum(mp.re(values[i]) > 0 for i in range(25))
                gap = mp.re(values[15] - values[14])
                q_minus = vectors[:, :15]
                q_plus = vectors[:, 15:]
                residual_norm = mp_operator_norm(q_plus.H * h_matrix * q_minus)
                eta_residual = (
                    mp.sin(mp.mpf("0.5") * mp.atan(2 * residual_norm / gap))
                    if gap > 0 else mp.inf
                )
                a_norm = mp_hermitian_norm(a_matrix)
                epsilon_pencil = (
                    epsilon_a / b_lower
                    + a_norm * epsilon_b / (b_min * b_lower)
                    + MP_FLOOR * max(mp.mpf(1), mp_hermitian_norm(h_matrix))
                )
                eta_source = (
                    2 * epsilon_pencil / (gap - 2 * epsilon_pencil)
                    + MP_FLOOR
                    if gap > 2 * epsilon_pencil else mp.inf
                )
                lifted = w_basis * lower_inverse.H * q_minus
                gram = (lifted.H * lifted + (lifted.H * lifted).H) / 2
                projector = lifted * mp.inverse(gram) * lifted.H
                projector = (projector + projector.H) / 2
                eta_shape_term = 2 * eta_s
                eta_residual_term = mp.sqrt(b_max / b_lower) * eta_residual
                eta_source_term = mp.sqrt(b_max / b_lower) * eta_source
                eta_metric_term = epsilon_b / b_lower
                eta_p = (
                    eta_shape_term + eta_residual_term + eta_source_term
                    + eta_metric_term + MP_FLOOR
                )
                hermiticity = mp_operator_norm(projector - projector.H)
                idempotence = mp_operator_norm(projector * projector - projector)
                relative_null = shape["null_residual"] / max(
                    mp.mpf(1), mp_operator_norm(u_basis.H * m_matrix)
                )
                complete = bool(
                    len(shape["pivots"]) == 5
                    and shape["orthogonality"] < CONTROL_FACTOR
                    and relative_null < CONTROL_FACTOR
                    and positive and negative == 15 and positive_count == 10
                    and gap > 2 * epsilon_pencil
                    and hermiticity < CONTROL_FACTOR
                    and idempotence < CONTROL_FACTOR
                    and mp.isfinite(eta_p)
                )
                projector_ok &= complete
                prior_projector = prior["projectors"][key]
                binary_distance = la.svdvals(
                    mp_to_numpy(projector) - prior_projector["projector"]
                )[0]
                binary_overlap = binary_distance <= prior_projector["eta"]
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
                    "source_M_radius": smp(epsilon_m),
                    "source_V_radius": smp(epsilon_v),
                    "kinetic_minimum": smp(b_min),
                    "kinetic_lower_bound": smp(b_lower),
                    "generalized_gap": smp(gap),
                    "whitened_residual": smp(residual_norm),
                    "pencil_source_error": smp(epsilon_pencil),
                    "projector_shape_term": smp(eta_shape_term),
                    "projector_residual_term": smp(eta_residual_term),
                    "projector_source_term": smp(eta_source_term),
                    "projector_metric_term": smp(eta_metric_term),
                    "projector_error": smp(eta_p),
                    "projector_hermiticity": smp(hermiticity),
                    "projector_idempotence": smp(idempotence),
                    "negative": negative,
                    "positive": positive_count,
                    "binary_projector_distance": f"{binary_distance:.17e}",
                    "binary_overlap": binary_overlap,
                    "complete": complete,
                })

check("all 32 high-precision generalized projectors certify",
      projector_ok and len(projectors) == len(projector_records) == 32)
check("all high-precision projectors overlap the committed binary controls",
      binary_overlap_ok)

# Control the 16 previously disclosed matched distances.
prior_identity = {
    (item["parity"], item["sector_index"], item["variant"]): item
    for item in prior["identity_records"]
}
matched_control_records = []
matched_distance_ok = True
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for variant in VARIANTS:
            old = projectors[("old", parity, sector_index, variant)]
            shifted = projectors[("shifted", parity, sector_index, variant)]
            distance = projector_distance(old["projector"], shifted["projector"])
            committed = prior_identity[(parity, sector_index, variant)]
            difference = abs(distance - mp.mpf(committed["distance"]))
            overlap = difference <= mp.mpf(committed["error"])
            matched_distance_ok &= overlap
            matched_control_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "distance": smp(distance),
                "committed_distance": committed["distance"],
                "difference": smp(difference),
                "committed_error": committed["error"],
                "overlap": overlap,
            })

check("all 16 matched distances overlap the committed comparisons",
      matched_distance_ok and len(matched_control_records) == 16)

# Complete within-time and all-pairs old/shifted comparisons.
within_records = []
within_counts = Counter()
for time_name in TIMES:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for left_variant, right_variant in combinations(VARIANTS, 2):
                left = projectors[(time_name, parity, sector_index, left_variant)]
                right = projectors[(time_name, parity, sector_index, right_variant)]
                distance = projector_distance(left["projector"], right["projector"])
                error = left["eta"] + right["eta"] + MP_FLOOR
                label = projector_label(distance, error)
                within_counts[label] += 1
                within_records.append({
                    "time": time_name,
                    "parity": parity,
                    "sector_index": sector_index,
                    "left_variant": left_variant,
                    "right_variant": right_variant,
                    "distance": smp(distance),
                    "error": smp(error),
                    "error_units": smp(distance / error),
                    "label": label,
                })

cross_records = []
cross_counts = Counter()
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for old_variant in VARIANTS:
            for shifted_variant in VARIANTS:
                old = projectors[("old", parity, sector_index, old_variant)]
                shifted = projectors[
                    ("shifted", parity, sector_index, shifted_variant)
                ]
                distance = projector_distance(
                    old["projector"], shifted["projector"]
                )
                error = old["eta"] + shifted["eta"] + MP_FLOOR
                label = projector_label(distance, error)
                cross_counts[label] += 1
                cross_records.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "old_variant": old_variant,
                    "shifted_variant": shifted_variant,
                    "distance": smp(distance),
                    "error": smp(error),
                    "error_units": smp(distance / error),
                    "label": label,
                })

comparison_counts_ok = (
    len(within_records) == sum(within_counts.values()) == 48
    and len(cross_records) == sum(cross_counts.values()) == 64
)
check("the complete 48 within-time and 64 cross-time comparisons exist",
      comparison_counts_ok,
      f"within={dict(within_counts)}, cross={dict(cross_counts)}")

family_records = []
family_robust = True
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        within_subset = [
            item for item in within_records
            if item["parity"] == parity and item["sector_index"] == sector_index
        ]
        cross_subset = [
            item for item in cross_records
            if item["parity"] == parity and item["sector_index"] == sector_index
        ]
        diameter_upper = max(
            mp.mpf(item["distance"]) + mp.mpf(item["error"])
            for item in within_subset
        )
        separation_lower = min(
            max(mp.mpf(0), mp.mpf(item["distance"]) - mp.mpf(item["error"]))
            for item in cross_subset
        )
        robustness = (
            separation_lower / diameter_upper
            if diameter_upper > 0 else mp.inf
        )
        robust = robustness > 100
        family_robust &= robust
        family_records.append({
            "parity": parity,
            "sector_index": sector_index,
            "family_diameter_upper": smp(diameter_upper),
            "cross_separation_lower": smp(separation_lower),
            "robustness_ratio": smp(robustness),
            "robust_above_100": robust,
        })

check("all four finite-family robustness ratios are assigned",
      len(family_records) == 4 and all(
          mp.isfinite(mp.mpf(item["family_diameter_upper"]))
          and mp.isfinite(mp.mpf(item["cross_separation_lower"]))
          for item in family_records
      ))

controls_ok = bool(
    provenance_ok and prior_replay_ok and reconstruction_ok
    and len(exact_cells) == 32 and conformal_ok and projector_ok
    and binary_overlap_ok and matched_distance_ok and comparison_counts_ok
    and len(projectors) == 32 and len(family_records) == 4
)
if not controls_ok:
    outcome = "RESIDUAL_BUNDLE_CONTROL_FAILED"
elif cross_counts["ROTATION_RESOLVED"] == 64 and family_robust:
    outcome = "RESIDUAL_FINITE_FAMILY_ROTATION_RESOLVED"
elif cross_counts["ROTATION_RESOLVED"] == 64:
    outcome = "RESIDUAL_FINITE_FAMILY_SCHEME_COMPETITION_OPEN"
elif cross_counts["ZERO_CONSISTENT"] == 64:
    outcome = "RESIDUAL_FINITE_FAMILY_ZERO_CONSISTENT"
else:
    outcome = "RESIDUAL_FINITE_FAMILY_ROTATION_OPEN"

allowed = {
    "RESIDUAL_BUNDLE_CONTROL_FAILED",
    "RESIDUAL_FINITE_FAMILY_ROTATION_RESOLVED",
    "RESIDUAL_FINITE_FAMILY_SCHEME_COMPETITION_OPEN",
    "RESIDUAL_FINITE_FAMILY_ZERO_CONSISTENT",
    "RESIDUAL_FINITE_FAMILY_ROTATION_OPEN",
}
check("the preregistered residual hierarchy assigns one outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commits": list(PROTOCOL_COMMITS),
    "input_sha256": hashes,
    "target_disclosed_upstream": True,
    "finite_derivative_family": list(VARIANTS),
    "outcome": outcome,
    "controls_ok": controls_ok,
    "conformal_records": conformal_records,
    "projector_records": projector_records,
    "matched_control_records": matched_control_records,
    "within_counts": dict(within_counts),
    "within_records": within_records,
    "cross_counts": dict(cross_counts),
    "cross_records": cross_records,
    "family_records": family_records,
    "classification": {
        "finite_family_result": "DERIVED COMPUTATIONAL",
        "analytic_hessian_rotation": "OPEN",
        "action_selected_connection": "OPEN",
        "reduced_propagator": "NOT COMPUTED",
        "dispersion_mass_inertia_or_speed": "NOT COMPUTED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"within labels: {dict(within_counts)}")
print(f"cross labels: {dict(cross_counts)}")
print("robustness ratios: " + str([
    item["robustness_ratio"] for item in family_records
]))
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

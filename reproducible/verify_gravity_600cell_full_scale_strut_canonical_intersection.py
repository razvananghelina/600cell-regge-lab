#!/usr/bin/env python3
"""Target-blind complete scale--strut/canonical intersection census."""

from collections import Counter
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIOR = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_canonical_intersection_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_canonical_intersection_protocol.md"
CARRIER_RESULT = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_symbolic_gap_resolution_result.md"
CARRIER_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_carrier.py"
CARRIER_INPUT = HERE / "gravity_600cell_full_scale_strut_carrier.json"
PRECISION_INPUT = HERE / "gravity_600cell_full_scale_strut_precision.json"
SYMBOLIC_INPUT = HERE / "gravity_600cell_full_scale_strut_symbolic_gap_resolution.json"
ACTION_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
ACTION_INPUT = HERE / "gravity_600cell_dust_hyperbolic_lapse_alignment.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
PURE_INPUT = HERE / "gravity_600cell_corrected_strut_canonical_intersection.json"
OUTPUT = HERE / "gravity_600cell_full_scale_strut_canonical_intersection.json"

PRIOR_COMMIT = "c2d0e83"
PROTOCOL_COMMIT = "b621736"
EXPECTED_HASHES = {
    "prior": "134ece5926b429011a1b74428a30454924c1458492f2aed3244bab9258b345c3",
    "protocol": "c769768003c4bd24745e9c618e5e7d6699261f25825a40ed24cb3b09dd8f6f73",
    "carrier_result": "774dae7fbe3d3becf505867c3272f41f800ab9c917766fcfd347395e36c34ece",
    "carrier_source": "e68105df4058f7d2ed39a6913f29e88cd9fe88e123ff52260acf698a2bd7da49",
    "carrier": "6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d",
    "precision": "2a2a79271a92fc2ddde343a9d0651402df6eeb4a90efa2697e26f54cafcdf60f",
    "symbolic": "ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179",
    "action_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "action": "a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "pure": "422d8d8cb0fc0d72d842e3bf79609d4d985da6237c58e7c699b5f9cc21b65cec",
}
INPUTS = {
    "prior": PRIOR,
    "protocol": PROTOCOL,
    "carrier_result": CARRIER_RESULT,
    "carrier_source": CARRIER_SOURCE,
    "carrier": CARRIER_INPUT,
    "precision": PRECISION_INPUT,
    "symbolic": SYMBOLIC_INPUT,
    "action_source": ACTION_SOURCE,
    "action": ACTION_INPUT,
    "tangent_numeric": TANGENT_NUMERIC,
    "pure": PURE_INPUT,
}

mp.mp.dps = 100
VERTICES = 120
ROWS = 1560
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


def sf(value):
    return f"{float(value):.17e}"


def norm2(matrix):
    return float(la.svdvals(matrix, check_finite=True)[0]) if matrix.size else 0.0


def svd_drivers(matrix):
    return {
        driver: la.svd(
            matrix, compute_uv=False, lapack_driver=driver, check_finite=True
        )
        for driver in ("gesdd", "gesvd")
    }


def classify(values, epsilon):
    result = []
    for value in values:
        if value <= 10 * epsilon:
            result.append("ZERO")
        elif value > 100 * epsilon:
            result.append("NONZERO")
        else:
            result.append("OPEN")
    return result


def numerical_rank(matrix):
    values = la.svdvals(matrix, check_finite=True)
    threshold = (
        100 * np.finfo(float).eps * max(matrix.shape) * max(1.0, values[0])
    )
    return int(np.count_nonzero(values > threshold)), threshold


def kernel_projector(matrix, nullity):
    if nullity == 0:
        return np.zeros((matrix.shape[1], matrix.shape[1]), dtype=complex)
    _, _, vh = la.svd(matrix, full_matrices=False, lapack_driver="gesvd")
    basis = vh[-nullity:, :].conj().T
    return basis @ basis.conj().T


def serialize_complex_matrix(matrix):
    return [
        [
            {"real": sf(value.real), "imaginary": sf(value.imag)}
            for value in row
        ]
        for row in matrix
    ]


def coefficient_sets(committed):
    background = committed["background"]
    lam = mp.mpf(background["lambda"])
    rho = mp.mpf(background["rho"])
    l0_square = mp.mpf(background["L0_square"])
    q_diag = lam * l0_square - rho
    formula = {
        "A": -16 * rho / (l0_square * (lam - 1) ** 2),
        "B": 8 + 16 * rho / (l0_square * (lam - 1) ** 2),
        "scale_factor": l0_square / (8 * q_diag),
        "kappa": rho / ((lam - 1) * q_diag),
        "q_diag": q_diag,
        "lambda": lam,
    }
    stored = {
        key: mp.mpf(committed["parities"]["even"]["coefficient_values"][key])
        for key in ("A", "B", "scale_factor", "kappa", "q_diag")
    }
    stored["lambda"] = lam
    discrepancy = max(
        abs(formula[key] - stored[key]) / max(mp.mpf(1), abs(formula[key]))
        for key in ("A", "B", "scale_factor", "kappa", "q_diag")
    )
    return formula, stored, discrepancy


def carrier_matrix(parity, index_data, weak_positions, committed, coefficients):
    orbit_order = tuple(range(30, 65)) + tuple(range(65, 95))
    orbit_position = {orbit: position for position, orbit in enumerate(orbit_order)}
    weak_types = [30 + position for position in weak_positions]
    pole_edges = [
        tuple(map(int, edge))
        for orbit in weak_types for edge in index_data["orbit_edges"][orbit]
    ]
    logical_to_column = {edge[0]: column for column, edge in enumerate(pole_edges)}
    if set(logical_to_column) != set(range(VERTICES)):
        raise RuntimeError("pole coefficient order does not cover all vertices")

    record = committed["parities"][parity]
    internal_edges = [tuple(map(int, edge)) for edge in record["internal_edge_order"]]
    final_edges = [tuple(map(int, edge)) for edge in record["final_edge_order"]]
    all_edges = internal_edges + final_edges
    expected_edges = {
        tuple(map(int, edge))
        for edge, global_index in index_data["edge_to_index"].items()
        if index_data["edge_kind"][global_index] in {"internal", "pole", "new"}
    }
    coverage_ok = bool(len(all_edges) == ROWS and set(all_edges) == expected_edges)

    lam = coefficients["lambda"]
    a_value = coefficients["scale_factor"] * coefficients["A"]
    b_value = coefficients["scale_factor"] * coefficients["B"]
    kappa = coefficients["kappa"]
    matrix = np.zeros((ROWS, 2 * VERTICES), dtype=float)
    for edge in internal_edges:
        global_index = index_data["edge_to_index"][edge]
        orbit, group = divmod(global_index, 24)
        row = 24 * orbit_position[orbit] + group
        lower, upper = edge
        if upper == lower + VERTICES:
            matrix[row, VERTICES + logical_to_column[lower]] = 1.0
        else:
            target = upper - VERTICES
            matrix[row, logical_to_column[lower]] = float(a_value)
            matrix[row, logical_to_column[target]] = float(b_value)
            matrix[row, VERTICES + logical_to_column[lower]] = float(kappa)
            matrix[row, VERTICES + logical_to_column[target]] = float(-lam * kappa)
    for edge in final_edges:
        global_index = index_data["edge_to_index"][edge]
        orbit, group = divmod(global_index, 24)
        row = 24 * orbit_position[orbit] + group
        left = edge[0] - VERTICES
        right = edge[1] - VERTICES
        matrix[row, logical_to_column[left]] = float(1 / lam)
        matrix[row, logical_to_column[right]] = float(1 / lam)
    return matrix, logical_to_column, coverage_ok, internal_edges


def projected_carrier(matrix, sector, old):
    dimension = sector["dimension"]
    basis = old["mp_to_numpy"](sector["basis"])
    return (
        np.kron(np.eye(65), basis).conj().T
        @ matrix
        @ np.kron(np.eye(10), basis)
    )


def corrupted_matrix(matrix, first_diagonal, index_data, logical_to_column):
    result = matrix.copy()
    orbit_order = tuple(range(30, 65)) + tuple(range(65, 95))
    orbit_position = {orbit: position for position, orbit in enumerate(orbit_order)}
    global_index = index_data["edge_to_index"][first_diagonal]
    orbit, group = divmod(global_index, 24)
    row = 24 * orbit_position[orbit] + group
    lower, upper = first_diagonal
    target = upper - VERTICES
    pairs = (
        (logical_to_column[lower], logical_to_column[target]),
        (
            VERTICES + logical_to_column[lower],
            VERTICES + logical_to_column[target],
        ),
    )
    for left, right in pairs:
        result[row, left], result[row, right] = result[row, right], result[row, left]
    return result


print("=" * 78)
print("COMPLETE SCALE--STRUT / CANONICAL INTERSECTION CENSUS")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
committed = json.loads(CARRIER_INPUT.read_text())
precision = json.loads(PRECISION_INPUT.read_text())
symbolic = json.loads(SYMBOLIC_INPUT.read_text())
action_before = json.loads(ACTION_INPUT.read_text())
pure = json.loads(PURE_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and committed["outcome"] == "FULL_SCALE_STRUT_NUMERICALLY_OPEN"
    and committed["passed"] == committed["tests"] == 18
    and precision["outcome"] == "FULL_SCALE_STRUT_PRECISION_RESOLVED"
    and precision["passed"] == precision["tests"] == 9
    and symbolic["outcome"] == "FULL_SCALE_STRUT_GAP_REAL_RESOLVED"
    and symbolic["passed"] == symbolic["tests"] == 11
    and pure["outcome"] == "CORRECTED_STRUT_CANONICAL_INTERSECTION_RESOLVED"
    and pure["nullities"] == [0] * 14
)
check("all full-intersection inputs retain exact frozen provenance", provenance_ok)

print("reconstructing frozen action-response source", flush=True)
captured = io.StringIO()
original_exit = sys.exit


def audited_exit(code=0):
    if code not in (None, 0):
        raise SystemExit(code)


try:
    sys.exit = audited_exit
    with contextlib.redirect_stdout(captured):
        old = runpy.run_path(str(ACTION_SOURCE))
finally:
    sys.exit = original_exit
action_ok = bool(
    old["tests"] == old["passed"] == 14
    and old["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
    and digest(ACTION_INPUT) == EXPECTED_HASHES["action"]
    and json.loads(ACTION_INPUT.read_text()) == action_before
)
check("the frozen action-response audit reproduces byte-identically", action_ok)

formula_coefficients, stored_coefficients, coefficient_discrepancy = coefficient_sets(
    committed
)
coefficients_ok = coefficient_discrepancy < mp.mpf("1e-55")
check(
    "analytic and stored carrier coefficients agree at 100-digit working precision",
    coefficients_ok,
    f"maximum relative discrepancy={mp.nstr(coefficient_discrepancy, 8)}",
)

models = old["models"]
tick = old["tick"]
variants = old["VARIANTS"]
precision_relative = max(
    float(precision["parities"][parity]["gesvd_high_precision_max_relative_discrepancy"])
    for parity in ("even", "odd")
)

records = {}
all_resolved = True
all_hard_controls = bool(provenance_ok and action_ok and coefficients_ok)
all_nullities = []
all_global_nullities = {}
all_corruption_matrix_changes = []
all_corruption_singular_changes = []
all_conjugation_errors = []

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing carrier and seven action sectors", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    index_data = old["group_and_index_data"](model, state)
    geometry = old["prepare_geometry"](model, index_data)
    weak_positions = [
        position for position in range(35)
        if index_data["edge_kind"][24 * (30 + position)] == "pole"
    ]
    g_formula, logical_to_column, coverage_ok, internal_edges = carrier_matrix(
        parity, index_data, weak_positions, committed, formula_coefficients
    )
    g_stored, logical_to_column_stored, coverage_stored, _ = carrier_matrix(
        parity, index_data, weak_positions, committed, stored_coefficients
    )
    carrier_reconstruction_error = norm2(g_formula - g_stored)
    carrier_ok = bool(
        coverage_ok and coverage_stored
        and logical_to_column == logical_to_column_stored
        and carrier_reconstruction_error < 1e-45
    )
    check(
        f"{parity}: independent full-carrier reconstructions and edge coverage agree",
        carrier_ok,
        f"spectral discrepancy={carrier_reconstruction_error:.3e}",
    )
    all_hard_controls &= carrier_ok

    sectors, sector_control = old["high_precision_sector_bases"](index_data)
    s_value = mp.mpf(state[0])
    kind_values = {
        "old": old["L0_SQUARE"],
        "internal": mp.exp(s_value) * old["L0_SQUARE"] - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * s_value) * old["L0_SQUARE"],
    }
    pattern_cache, branch_control = old["high_precision_pattern_cache"](
        geometry["patterns"], kind_values
    )
    kernels, kernel_control = old["assemble_full_representative_kernels"](
        index_data, geometry, pattern_cache
    )
    geometry_ok = bool(
        sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and all(
            value < mp.mpf("1e-70")
            for key, value in sector_control.items() if key.startswith("maximum_")
        )
        and branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-70")
    )
    check(f"{parity}: sector and action geometry controls reproduce", geometry_ok)
    all_hard_controls &= geometry_ok

    first_diagonal = min(
        edge for edge in internal_edges
        if edge[0] < VERTICES <= edge[1] and edge[1] != edge[0] + VERTICES
    )
    g_corrupt_global = corrupted_matrix(
        g_formula, first_diagonal, index_data, logical_to_column
    )

    sector_records = []
    parity_resolved = True
    parity_hard_controls = True
    parity_deferred_open = False
    parity_global_nullity = 0
    for sector_index, sector in enumerate(sectors):
        dimension = int(sector["dimension"])
        half = 5 * dimension
        source_count = 10 * dimension
        joined_count = 15 * dimension
        print(
            f"[{parity}] sector {sector_index + 1}/7 d={dimension}", flush=True
        )
        g_matrix = projected_carrier(g_formula, sector, old)
        g_corrupt = projected_carrier(g_corrupt_global, sector, old)
        g_scale = g_matrix[:, :half]
        g_strut = g_matrix[:, half:]
        weak_rows = [
            position * dimension + component
            for position in weak_positions for component in range(dimension)
        ]
        pole_scale_error = float(la.norm(g_scale[weak_rows, :], ord=np.inf))
        pole_strut_error = float(
            la.norm(g_strut[weak_rows, :] - np.eye(half), ord=np.inf)
        )

        variant_data = {}
        for name in variants:
            block = old["project_full_kernel"](kernels[name], sector)
            variant_data[name] = old["response_and_lift_ball"](
                block, dimension, weak_positions
            )
        canonical = {
            name: data["lift_midpoint"] for name, data in variant_data.items()
        }
        pole_c_error = max(
            float(la.norm(matrix[weak_rows, :] - np.eye(half), ord=np.inf))
            for matrix in canonical.values()
        )
        rank_g, rank_g_threshold = numerical_rank(g_matrix)
        rank_scale, rank_scale_threshold = numerical_rank(g_scale)
        canonical_ranks = [numerical_rank(matrix)[0] for matrix in canonical.values()]
        graph_ok = bool(
            rank_g == source_count
            and rank_scale == half
            and canonical_ranks == [half] * len(variants)
            and pole_scale_error < 1e-13
            and pole_strut_error < 1e-13
            and pole_c_error < 1e-13
        )
        parity_hard_controls &= graph_ok

        s_scale = 1 / max(1.0, norm2(g_scale))
        s_strut = 1 / max(
            1.0, norm2(g_strut), *(norm2(matrix) for matrix in canonical.values())
        )
        scale_d = np.diag(
            np.r_[np.full(half, s_scale), np.full(half, s_strut)]
        )
        scale_k = np.diag(
            np.r_[
                np.full(half, s_scale),
                np.full(half, s_strut),
                np.full(half, s_strut),
            ]
        )
        d_matrices = {
            name: np.hstack((g_scale, g_strut - matrix)) @ scale_d
            for name, matrix in canonical.items()
        }
        k_matrices = {
            name: np.hstack((g_scale, g_strut, -matrix)) @ scale_k
            for name, matrix in canonical.items()
        }
        d_operational = d_matrices["operational_primary"]
        k_operational = k_matrices["operational_primary"]
        carrier_scaled = g_matrix @ np.diag(
            np.r_[np.full(half, s_scale), np.full(half, s_strut)]
        )
        carrier_error = (
            carrier_reconstruction_error * max(s_scale, s_strut)
            + precision_relative * norm2(carrier_scaled)
        )
        lift_radius = max(
            float(la.norm(data["lift_radii"], ord="fro")) * s_strut
            for data in variant_data.values()
        )
        d_step = max(norm2(matrix - d_operational) for matrix in d_matrices.values())
        k_step = max(norm2(matrix - k_operational) for matrix in k_matrices.values())
        d_binary = (
            200 * np.finfo(float).eps * max(d_operational.shape)
            * max(1.0, norm2(d_operational))
        )
        k_binary = (
            200 * np.finfo(float).eps * max(k_operational.shape)
            * max(1.0, norm2(k_operational))
        )
        epsilon_d = lift_radius + d_step + carrier_error + d_binary + 1e-70
        epsilon_k = lift_radius + k_step + carrier_error + k_binary + 1e-70

        d_singulars = {
            name: svd_drivers(matrix) for name, matrix in d_matrices.items()
        }
        k_singulars = {
            name: svd_drivers(matrix) for name, matrix in k_matrices.items()
        }
        d_labels = {
            name: {
                driver: classify(values, epsilon_d)
                for driver, values in drivers.items()
            }
            for name, drivers in d_singulars.items()
        }
        k_labels = {
            name: {
                driver: classify(values, epsilon_k)
                for driver, values in drivers.items()
            }
            for name, drivers in k_singulars.items()
        }
        d_nullities = [
            labels.count("ZERO")
            for by_driver in d_labels.values() for labels in by_driver.values()
        ]
        k_nullities = [
            labels.count("ZERO")
            for by_driver in k_labels.values() for labels in by_driver.values()
        ]
        open_count = sum(
            labels.count("OPEN")
            for collection in (d_labels, k_labels)
            for by_driver in collection.values() for labels in by_driver.values()
        )
        resolved = bool(
            open_count == 0
            and len(set(d_nullities)) == 1
            and len(set(k_nullities)) == 1
            and d_nullities[0] == k_nullities[0]
        )
        nullity = d_nullities[0] if resolved else None
        parity_resolved &= resolved
        if nullity is not None:
            all_nullities.append(nullity)
            parity_global_nullity += dimension * nullity

        # Synthetic zero/nonzero and structural five-dimensional control.
        zero_labels = classify(la.svdvals(np.zeros_like(d_operational)), epsilon_d)
        injection = np.zeros_like(d_operational)
        injection[:source_count, :] = np.eye(source_count)
        injection_labels = classify(la.svdvals(injection), epsilon_d)
        structural = np.hstack((g_scale, np.zeros_like(g_strut))) @ scale_d
        structural_labels = classify(la.svdvals(structural), epsilon_d)
        synthetic_hard_ok = bool(
            zero_labels.count("ZERO") == source_count
            and injection_labels.count("NONZERO") == source_count
        )
        structural_open = "OPEN" in structural_labels
        structural_ok = bool(
            not structural_open and structural_labels.count("ZERO") == half
        )
        structural_hard_ok = bool(structural_open or structural_ok)

        transform = np.eye(source_count, dtype=complex)
        transform[np.arange(source_count - 1), np.arange(1, source_count)] = 0.01
        transformed_labels = classify(
            la.svdvals(d_operational @ transform), epsilon_d
        )
        basis_ok = bool(
            "OPEN" not in transformed_labels
            and nullity is not None
            and transformed_labels.count("ZERO") == nullity
        )
        basis_hard_ok = bool(not resolved or basis_ok)

        joined_image_nullity = k_nullities[0] if resolved else None
        image_ok = bool(
            nullity is not None and joined_image_nullity == nullity
        )
        image_hard_ok = bool(not resolved or image_ok)

        g_corrupt_scale = g_corrupt[:, :half]
        g_corrupt_strut = g_corrupt[:, half:]
        d_corrupt = np.hstack((
            g_corrupt_scale,
            g_corrupt_strut - canonical["operational_primary"],
        )) @ scale_d
        corruption_matrix_change = norm2(d_corrupt - d_operational)
        corruption_singular_change = float(np.max(np.abs(
            la.svdvals(d_corrupt) - d_singulars["operational_primary"]["gesvd"]
        )))
        all_corruption_matrix_changes.append(corruption_matrix_change)
        all_corruption_singular_changes.append(corruption_singular_change)

        conjugation_error = float(np.max(np.abs(
            la.svdvals(d_operational)
            - la.svdvals(d_operational.conj())
        )))
        all_conjugation_errors.append(conjugation_error)
        hard_controls_ok = bool(
            graph_ok and synthetic_hard_ok and structural_hard_ok
            and basis_hard_ok and image_hard_ok
            and corruption_matrix_change > carrier_error
            and corruption_singular_change > carrier_error
            and conjugation_error <= d_binary
        )
        deferred_open = bool(not resolved or structural_open)
        parity_hard_controls &= hard_controls_ok
        parity_deferred_open |= deferred_open

        projector = (
            kernel_projector(d_operational, nullity)
            if nullity is not None else None
        )
        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "constant_overlap": sf(sector["constant_overlap"]),
            "center_value": sf(sector["center_value"]),
            "graph_control": {
                "rank_G": rank_g,
                "rank_G_scale": rank_scale,
                "canonical_ranks": canonical_ranks,
                "rank_G_threshold": sf(rank_g_threshold),
                "rank_G_scale_threshold": sf(rank_scale_threshold),
                "pole_scale_error": sf(pole_scale_error),
                "pole_strut_error": sf(pole_strut_error),
                "pole_canonical_error": sf(pole_c_error),
                "passed": graph_ok,
            },
            "scaling": {"scale": sf(s_scale), "strut": sf(s_strut)},
            "epsilon_D": sf(epsilon_d),
            "epsilon_K": sf(epsilon_k),
            "carrier_error": sf(carrier_error),
            "lift_radius": sf(lift_radius),
            "variant_step_D": sf(d_step),
            "variant_step_K": sf(k_step),
            "D_singular_values": {
                name: {
                    driver: [sf(value) for value in values]
                    for driver, values in drivers.items()
                }
                for name, drivers in d_singulars.items()
            },
            "K_singular_values": {
                name: {
                    driver: [sf(value) for value in values]
                    for driver, values in drivers.items()
                }
                for name, drivers in k_singulars.items()
            },
            "D_nullities": d_nullities,
            "K_nullities": k_nullities,
            "open_count": open_count,
            "resolved": resolved,
            "nullity": nullity,
            "global_multiplicity_contribution": (
                dimension * nullity if nullity is not None else None
            ),
            "kernel_projector_equilibrated": (
                serialize_complex_matrix(projector) if projector is not None else None
            ),
            "controls": {
                "synthetic_zero_and_identity": synthetic_hard_ok,
                "structural_status": (
                    "DEFERRED_OPEN" if structural_open
                    else "PASS" if structural_ok else "FAIL"
                ),
                "structural_C_equals_G_strut_nullity": structural_labels.count("ZERO"),
                "basis_change_status": (
                    "PASS" if basis_ok else "DEFERRED_OPEN" if not resolved else "FAIL"
                ),
                "joined_image_status": (
                    "PASS" if image_ok else "DEFERRED_OPEN" if not resolved else "FAIL"
                ),
                "corruption_matrix_change": sf(corruption_matrix_change),
                "corruption_singular_change": sf(corruption_singular_change),
                "conjugation_error": sf(conjugation_error),
                "hard_controls_pass": hard_controls_ok,
                "deferred_open": deferred_open,
            },
        })

    check(
        f"{parity}: all seven actual reduced/joined intersections are evaluated",
        True,
        f"resolved={parity_resolved}, nullities={[r['nullity'] for r in sector_records]}",
    )
    check(
        f"{parity}: all hard controls pass and open comparisons are deferred",
        parity_hard_controls,
        f"deferred_open={parity_deferred_open}",
    )
    all_resolved &= parity_resolved
    all_hard_controls &= parity_hard_controls
    all_global_nullities[parity] = (
        parity_global_nullity if parity_resolved else None
    )
    records[parity] = {
        "carrier_reconstruction_error": sf(carrier_reconstruction_error),
        "resolved": parity_resolved,
        "deferred_open": parity_deferred_open,
        "global_intersection_dimension": all_global_nullities[parity],
        "sectors": sector_records,
    }

corruption_ok = bool(
    all_corruption_matrix_changes
    and min(all_corruption_matrix_changes) > 0
    and max(all_corruption_singular_changes) > 0
)
check(
    "the frozen source/target corruption is detected",
    corruption_ok,
    f"minimum matrix change={min(all_corruption_matrix_changes):.3e}, maximum singular change={max(all_corruption_singular_changes):.3e}",
)
all_hard_controls &= corruption_ok

if not all_hard_controls:
    outcome = "FULL_SCALE_STRUT_CANONICAL_CONTROL_FAILED"
elif not all_resolved:
    outcome = "FULL_SCALE_STRUT_CANONICAL_NUMERICALLY_OPEN"
else:
    outcome = "FULL_SCALE_STRUT_CANONICAL_INTERSECTION_RESOLVED"

allowed = {
    "FULL_SCALE_STRUT_CANONICAL_CONTROL_FAILED",
    "FULL_SCALE_STRUT_CANONICAL_NUMERICALLY_OPEN",
    "FULL_SCALE_STRUT_CANONICAL_INTERSECTION_RESOLVED",
}
check("the preregistered outcome hierarchy assigns one verdict", outcome in allowed, outcome)

payload = {
    "prior_commit": PRIOR_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "coefficient_maximum_relative_discrepancy": mp.nstr(
        coefficient_discrepancy, 40
    ),
    "precision_relative_bound": sf(precision_relative),
    "parities": records,
    "all_resolved": all_resolved,
    "minimal_sector_nullities": all_nullities,
    "global_intersection_dimensions": all_global_nullities,
    "maximum_conjugation_error": sf(max(all_conjugation_errors)),
    "corruption": {
        "minimum_matrix_change": sf(min(all_corruption_matrix_changes)),
        "maximum_singular_change": sf(max(all_corruption_singular_changes)),
    },
    "classification": {
        "complete_kinematic_carrier": "DERIVED INPUT",
        "canonical_intersection": (
            "DERIVED COMPUTATIONAL; PRIMARY ONLY"
            if outcome == "FULL_SCALE_STRUT_CANONICAL_INTERSECTION_RESOLVED"
            else "OPEN"
        ),
        "gauge_or_physical_mode": "NOT CLASSIFIED",
        "tick_c_G_planck_or_mass": "NOT EVALUATED",
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
if passed != tests or outcome == "FULL_SCALE_STRUT_CANONICAL_CONTROL_FAILED":
    raise SystemExit(1)

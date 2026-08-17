#!/usr/bin/env python3
"""Blind 2T-resolved 1440-dimensional boundary tangent census.

Prior-art commit: 5dfc2a7.
Protocol commit: bc114bf.
No continuum harmonic, speed, Planck or particle target is loaded.
"""

import ast
from collections import Counter, defaultdict
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
import json
import math
from pathlib import Path
import sys
import zipfile

from flint import acb, acb_mat, ctx
import mpmath as mp
import numpy as np
import scipy.linalg as la
from scipy.optimize import linear_sum_assignment


HERE = Path(__file__).resolve().parent
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
GLUING_INPUT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
REDUCED_INPUT = HERE / "gravity_600cell_dust_dynamic_tangent.json"
RANK_INPUT = HERE / "gravity_600cell_dust_full_anisotropic_legendre_rank.json"
SCHUR_INPUT = HERE / "gravity_600cell_dust_full_lapse_schur.json"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
OUTPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
NUMERIC_OUTPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"

PRIOR_ART_COMMIT = "5dfc2a7"
PROTOCOL_COMMIT = "bc114bf"
EXPECTED_HASHES = {
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "gluing": "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77",
    "reduced_tangent": "1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5",
    "full_rank": "7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226",
    "full_schur": "4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}

DPS = 100
BALL_DPS = 80
mp.mp.dps = DPS
ctx.dps = BALL_DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-20"),
    "operational_shadow": mp.mpf("1e-15"),
    "validation_primary": mp.mpf("3e-20"),
    "validation_shadow": mp.mpf("3e-15"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {hinge: index for index, hinge in enumerate(LOCAL_HINGES)}
I = mp.mpc(0, 1)
ARITHMETIC_FLOOR = mp.mpf("1e-70")
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


def norm2(matrix):
    values = la.svdvals(np.asarray(matrix))
    return float(values[0]) if len(values) else 0.0


hashes = {
    "tick": sha256(TICK_INPUT),
    "gluing": sha256(GLUING_INPUT),
    "reduced_tangent": sha256(REDUCED_INPUT),
    "full_rank": sha256(RANK_INPUT),
    "full_schur": sha256(SCHUR_INPUT),
    "rank_source": sha256(RANK_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
}
tick = json.loads(TICK_INPUT.read_text())
gluing = json.loads(GLUING_INPUT.read_text())
reduced_input = json.loads(REDUCED_INPUT.read_text())
rank_input = json.loads(RANK_INPUT.read_text())
schur_input = json.loads(SCHUR_INPUT.read_text())

provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and tick["passed"] == tick["tests"] == 7
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and reduced_input["outcome"] in {
        "DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT",
        "DYNAMIC_SHAPE_TANGENT_SCHEDULE_ROBUST",
        "DYNAMIC_SHAPE_TANGENT_SCHEDULE_OPEN",
    }
    and rank_input["outcome"] == "FULL_CANONICAL_LEGENDRE_REGULAR"
    and all(
        item["full_resolved_rank"] == 1560
        and item["full_error_consistent_nullity"] == 0
        and item["full_numerically_open_count"] == 0
        for item in rank_input["parities"].values()
    )
    and schur_input["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
    and all(
        item["resolved_schur_rank"] == 120
        and item["schur_zero_count"] == 0
        and item["schur_open_count"] == 0
        for item in schur_input["parities"].values()
    )
)
check("all preregistered inputs have exact frozen provenance", provenance_ok, str(hashes))


spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_full_boundary_tangent", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check(
    "the direct one-slab geometry import retains all 43 certificates",
    gro.tests == gro.passed == 43,
)


M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON3 * L0
RHO0 = mp.mpf("0.0102") ** 2


def load_audited_functions():
    wanted = {
        "orbit_sort_key",
        "augment_boundary_orbits",
        "log_minus",
        "signed_volume_square",
        "angle_record",
        "area_data",
        "extended_edge_image",
        "group_and_index_data",
        "prepare_geometry",
    }
    tree = ast.parse(RANK_SOURCE.read_text(), filename=str(RANK_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch: missing={wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(RANK_SOURCE), "exec"),
        globals(),
    )


load_audited_functions()
models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}


def mp_frobenius(matrix):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in matrix))


def mp_submatrix(matrix, rows, columns):
    return mp.matrix([[matrix[row, column] for column in columns] for row in rows])


def cluster_sorted(values, tolerance=mp.mpf("1e-70")):
    clusters = []
    for index, value in enumerate(values):
        if not clusters or abs(value - values[clusters[-1][0]]) > tolerance:
            clusters.append([index])
        else:
            clusters[-1].append(index)
    return clusters


def high_precision_sector_bases(index_data):
    table = index_data["table"]
    regular = []
    for element in range(24):
        matrix = mp.matrix(24, 24)
        for column in range(24):
            matrix[int(table[element, column]), column] = 1
        regular.append(matrix)

    central = mp.matrix(24, 24)
    for class_index, conjugacy_class in enumerate(index_data["classes"]):
        for element in conjugacy_class:
            central += (2**class_index) * regular[element]
    hermitian_center = (
        (central + central.H) / 2
        + mp.sqrt(2) * (central - central.H) / (2 * I)
    )
    center_values_matrix, center_vectors = mp.eighe(hermitian_center)
    center_values = [center_values_matrix[index] for index in range(24)]
    center_clusters = cluster_sorted(center_values)
    sectors = []
    maxima = {
        "orthonormal": mp.mpf(0),
        "center": mp.mpf(0),
        "splitter": mp.mpf(0),
        "right_leakage": mp.mpf(0),
        "conjugate_pair": mp.mpf(0),
    }

    for cluster in center_clusters:
        isotypic = mp_submatrix(center_vectors, range(24), cluster)
        dimension = math.isqrt(len(cluster))
        if dimension**2 != len(cluster):
            raise RuntimeError("central eigenspace is not a square isotypic block")
        splitter = None
        selected = None
        splitter_value = mp.mpf(0)
        splitter_matrix = mp.matrix(24, 24)
        if dimension == 1:
            selected = isotypic
        else:
            for element in range(24):
                if element == index_data["identity"]:
                    continue
                y_matrix = I * (regular[element] - regular[element].H)
                restricted = isotypic.H * y_matrix * isotypic
                split_values_matrix, split_vectors = mp.eighe(restricted)
                split_values = [
                    split_values_matrix[index] for index in range(dimension**2)
                ]
                split_clusters = cluster_sorted(split_values)
                if len(split_clusters) == dimension and all(
                    len(item) == dimension for item in split_clusters
                ):
                    splitter = element
                    selected = isotypic * mp_submatrix(
                        split_vectors, range(dimension**2), split_clusters[0]
                    )
                    splitter_value = split_values[split_clusters[0][0]]
                    splitter_matrix = y_matrix
                    break
        if selected is None:
            raise RuntimeError("deterministic minimal-sector splitter failed")

        center_value = mp.fsum(center_values[index] for index in cluster) / len(cluster)
        orthonormal = selected.H * selected - mp.eye(dimension)
        center_residual = hermitian_center * selected - selected * center_value
        splitter_residual = (
            splitter_matrix * selected - selected * splitter_value
            if dimension > 1 else mp.matrix(24, dimension)
        )
        compressed_center = selected.H * central * selected
        old_center = mp.fsum(
            compressed_center[index, index] for index in range(dimension)
        ) / dimension

        right_representations = []
        projector = selected * selected.H
        right_leakage = mp.mpf(0)
        for element in range(24):
            right = mp.matrix(24, 24)
            for row in range(24):
                column = int(table[row, element])
                right[row, column] = 1
            compressed = selected.H * right * selected
            right_representations.append(compressed)
            right_leakage = max(
                right_leakage,
                mp_frobenius(right * selected - selected * compressed),
            )

        constant = mp.matrix([1 / mp.sqrt(24) for _ in range(24)])
        constant_overlap = mp_frobenius(selected.H * constant) ** 2
        conjugate_pair = selected.T * selected.apply(mp.conj) - mp.eye(dimension)
        maxima["orthonormal"] = max(maxima["orthonormal"], mp_frobenius(orthonormal))
        maxima["center"] = max(maxima["center"], mp_frobenius(center_residual))
        maxima["splitter"] = max(maxima["splitter"], mp_frobenius(splitter_residual))
        maxima["right_leakage"] = max(maxima["right_leakage"], right_leakage)
        maxima["conjugate_pair"] = max(
            maxima["conjugate_pair"], mp_frobenius(conjugate_pair)
        )
        sectors.append({
            "dimension": dimension,
            "basis": selected,
            "right_representations": right_representations,
            "splitter": splitter,
            "center_value": center_value,
            "old_central_eigenvalue": old_center,
            "constant_overlap": constant_overlap,
        })
    return sectors, {
        "isotypic_dimensions": sorted(len(item) for item in center_clusters),
        "irrep_dimensions": sorted(item["dimension"] for item in sectors),
        **{f"maximum_{key}": value for key, value in maxima.items()},
    }


def high_precision_pattern_cache(patterns, kind_values):
    cache = {}
    minimum_minor = mp.inf
    minimum_argument = mp.inf
    maximum_raw_imaginary = mp.mpf(0)
    maximum_cross = mp.mpf(0)
    maximum_proxy = mp.mpf(0)
    entry_pass = True
    base_negative_counts = Counter()
    displaced_negative_counts = Counter()
    for pattern in sorted(patterns):
        values = [kind_values[kind] for kind in pattern]
        base_angles, base_branch = angle_record(values)
        base_negative_counts[base_branch["negative_directions"]] += patterns[pattern]
        minimum_minor = min(minimum_minor, base_branch["minimum_leading_minor"])
        minimum_argument = min(minimum_argument, base_branch["minimum_argument"])
        maximum_raw_imaginary = max(
            maximum_raw_imaginary, *(abs(mp.im(value)) for value in base_angles)
        )
        derivatives = {}
        for name, step in DERIVATIVE_STEPS.items():
            matrix = mp.matrix(10, 10)
            for column in range(10):
                plus = list(values)
                minus = list(values)
                plus[column] *= mp.exp(step)
                minus[column] *= mp.exp(-step)
                plus_angles, plus_branch = angle_record(plus)
                minus_angles, minus_branch = angle_record(minus)
                for branch in (plus_branch, minus_branch):
                    displaced_negative_counts[branch["negative_directions"]] += 1
                    minimum_minor = min(minimum_minor, branch["minimum_leading_minor"])
                    minimum_argument = min(minimum_argument, branch["minimum_argument"])
                for row in range(10):
                    matrix[row, column] = (
                        plus_angles[row] - minus_angles[row]
                    ) / (2 * step)
                    maximum_raw_imaginary = max(
                        maximum_raw_imaginary, abs(mp.im(matrix[row, column]))
                    )
            derivatives[name] = matrix
        for row in range(10):
            for column in range(10):
                op = derivatives["operational_primary"][row, column]
                op_shadow = derivatives["operational_shadow"][row, column]
                val = derivatives["validation_primary"][row, column]
                val_shadow = derivatives["validation_shadow"][row, column]
                cross = abs(op - val)
                proxy = abs(op - op_shadow) + abs(val - val_shadow) + ARITHMETIC_FLOOR
                maximum_cross = max(maximum_cross, cross)
                maximum_proxy = max(maximum_proxy, proxy)
                entry_pass &= bool(cross <= 10 * proxy)
        cache[pattern] = {
            "base_angles": tuple(base_angles),
            "derivatives": derivatives,
        }
    return cache, {
        "minimum_leading_minor": minimum_minor,
        "minimum_argument": minimum_argument,
        "maximum_raw_angle_or_derivative_imaginary": maximum_raw_imaginary,
        "maximum_cross": maximum_cross,
        "maximum_proxy": maximum_proxy,
        "entry_pass": bool(entry_pass),
        "base_negative_counts": base_negative_counts,
        "displaced_negative_counts": displaced_negative_counts,
    }


def assemble_full_representative_kernels(index_data, geometry, pattern_cache):
    signed_base = index_data["signed_base"]
    identity = index_data["identity"]
    curvature = [
        mp.pi if record["boundary"] else 2 * mp.pi
        for record in geometry["triangle_records"]
    ]
    for simplex in geometry["simplex_records"]:
        angles = pattern_cache[simplex["pattern"]]["base_angles"]
        for hinge_index, triangle in enumerate(simplex["hinge_triangles"]):
            curvature[triangle] += angles[hinge_index]

    kernels = {name: defaultdict(lambda: mp.mpc(0)) for name in VARIANTS}

    def location(global_index):
        return divmod(global_index, 24)

    def add_all(row_index, column_index, value):
        row_type, row_group = location(row_index)
        if row_group != identity:
            return
        column_type, column_group = location(column_index)
        key = (row_type, column_type, column_group)
        for kernel in kernels.values():
            kernel[key] += value

    def add_angle(name, row_index, column_index, value):
        row_type, row_group = location(row_index)
        if row_group != identity:
            return
        column_type, column_group = location(column_index)
        kernels[name][(row_type, column_type, column_group)] += value

    triangle_local = []
    for triangle_index, record in enumerate(geometry["triangle_records"]):
        indices = record["indices"]
        values = [signed_base[index] for index in indices]
        _, area_gradient, area_hessian = area_data(values)
        coefficient = -I * curvature[triangle_index]
        for row, row_index in enumerate(indices):
            for column, column_index in enumerate(indices):
                add_all(
                    row_index, column_index,
                    coefficient * area_hessian[row][column],
                )
        triangle_local.append(tuple(area_gradient))

    dust_hessian = -(2 * mp.pi * MASS / 120) * mp.sqrt(index_data["rho"])
    for pole_index in geometry["pole_indices"]:
        add_all(pole_index, pole_index, dust_hessian)

    for simplex in geometry["simplex_records"]:
        simplex_indices = simplex["indices"]
        derivatives = pattern_cache[simplex["pattern"]]["derivatives"]
        for hinge_index, triangle in enumerate(simplex["hinge_triangles"]):
            triangle_indices = geometry["triangle_records"][triangle]["indices"]
            area_gradient = triangle_local[triangle]
            for row, row_index in enumerate(triangle_indices):
                if location(row_index)[1] != identity:
                    continue
                for column, column_index in enumerate(simplex_indices):
                    for name in VARIANTS:
                        add_angle(
                            name, row_index, column_index,
                            -I * area_gradient[row]
                            * derivatives[name][hinge_index, column],
                        )

    maximum_imaginary = max(
        (abs(mp.im(value)) for kernel in kernels.values() for value in kernel.values()),
        default=mp.mpf(0),
    )
    return kernels, {
        "nonzero_entries": {name: len(kernel) for name, kernel in kernels.items()},
        "maximum_imaginary": maximum_imaginary,
    }


def project_full_kernel(kernel, sector):
    dimension = sector["dimension"]
    block = mp.matrix(95 * dimension, 95 * dimension)
    right = sector["right_representations"]
    for (row_type, column_type, group), value in kernel.items():
        representation = right[group]
        for row in range(dimension):
            for column in range(dimension):
                block[
                    row_type * dimension + row,
                    column_type * dimension + column,
                ] += value * representation[row, column]
    return block


def mp_to_numpy(matrix):
    return np.array([
        [complex(float(mp.re(matrix[row, column])), float(mp.im(matrix[row, column])))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


def mp_to_acb(value):
    return acb(mp.nstr(mp.re(value), 90), mp.nstr(mp.im(value), 90))


def mp_matrix_to_acb(matrix):
    return acb_mat(matrix.rows, matrix.cols, [mp_to_acb(value) for value in matrix])


def acb_midpoint_and_radii(matrix):
    midpoint = np.empty((matrix.nrows(), matrix.ncols()), dtype=np.complex128)
    radii = np.empty((matrix.nrows(), matrix.ncols()), dtype=float)
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            value = matrix[row, column]
            midpoint[row, column] = complex(
                float(value.real.mid()), float(value.imag.mid())
            )
            radii[row, column] = math.hypot(
                float(value.real.rad().upper()), float(value.imag.rad().upper())
            )
    return midpoint, radii


def expanded_types(start, stop, dimension):
    return [
        kind * dimension + component
        for kind in range(start, stop) for component in range(dimension)
    ]


def canonical_from_full_block(block, dimension):
    old = expanded_types(0, 30, dimension)
    internal = expanded_types(30, 65, dimension)
    new = expanded_types(65, 95, dimension)
    top = mp_submatrix(block, internal, internal + new)
    bottom = -mp_submatrix(block, old, internal + new)
    result = mp.matrix(len(internal) + len(old), len(internal) + len(new))
    for row in range(len(internal)):
        for column in range(len(internal) + len(new)):
            result[row, column] = top[row, column]
    for row in range(len(old)):
        for column in range(len(internal) + len(new)):
            result[len(internal) + row, column] = bottom[row, column]
    return result


def build_tangent_ball(block, dimension, mapping):
    old = expanded_types(0, 30, dimension)
    internal = expanded_types(30, 65, dimension)
    new = expanded_types(65, 95, dimension)
    nd = 30 * dimension
    xd = 35 * dimension

    k_xx = mp_submatrix(block, internal, internal)
    k_xn = mp_submatrix(block, internal, new)
    k_ox = mp_submatrix(block, old, internal)
    k_on = mp_submatrix(block, old, new)
    j_matrix = mp.matrix(xd + nd, xd + nd)
    for row in range(xd):
        for column in range(xd):
            j_matrix[row, column] = k_xx[row, column]
        for column in range(nd):
            j_matrix[row, xd + column] = k_xn[row, column]
    for row in range(nd):
        for column in range(xd):
            j_matrix[xd + row, column] = -k_ox[row, column]
        for column in range(nd):
            j_matrix[xd + row, xd + column] = -k_on[row, column]

    rhs = mp.matrix(xd + nd, 2 * nd)
    k_xo = mp_submatrix(block, internal, old)
    k_oo = mp_submatrix(block, old, old)
    for row in range(xd):
        for column in range(nd):
            rhs[row, column] = -k_xo[row, column]
    for row in range(nd):
        for column in range(nd):
            rhs[xd + row, column] = k_oo[row, column]
        rhs[xd + row, nd + row] = 1

    j_ball = mp_matrix_to_acb(j_matrix)
    rhs_ball = mp_matrix_to_acb(rhs)
    det_j = j_ball.det()
    solved = j_ball.solve(rhs_ball)
    y_x = acb_mat(xd, 2 * nd, [
        solved[row, column] for row in range(xd) for column in range(2 * nd)
    ])
    y_n = acb_mat(nd, 2 * nd, [
        solved[xd + row, column]
        for row in range(nd) for column in range(2 * nd)
    ])

    k_no = mp_matrix_to_acb(mp_submatrix(block, new, old))
    k_nx = mp_matrix_to_acb(mp_submatrix(block, new, internal))
    k_nn = mp_matrix_to_acb(mp_submatrix(block, new, new))
    direct = acb_mat(nd, 2 * nd)
    for row in range(nd):
        for column in range(nd):
            direct[row, column] = k_no[row, column]
    p_post = direct + k_nx * y_x + k_nn * y_n
    raw = acb_mat(2 * nd, 2 * nd, [
        (y_n[row, column] if row < nd else p_post[row - nd, column])
        for row in range(2 * nd) for column in range(2 * nd)
    ])

    tangent = acb_mat(2 * nd, 2 * nd)
    for next_type, final_type in enumerate(mapping):
        for component in range(dimension):
            target = next_type * dimension + component
            source = final_type * dimension + component
            for column in range(2 * nd):
                tangent[target, column] = raw[source, column]
                tangent[nd + target, column] = raw[nd + source, column]

    omega = acb_mat(2 * nd, 2 * nd)
    for index in range(nd):
        omega[index, nd + index] = 1
        omega[nd + index, index] = -1
    defect = tangent.transpose().conjugate() * omega * tangent - omega
    return j_matrix, det_j, tangent, defect


def central_match(sector, stored_sectors, used):
    target = complex(
        float(mp.re(sector["old_central_eigenvalue"])),
        float(mp.im(sector["old_central_eigenvalue"])),
    )
    choices = []
    for index, stored in enumerate(stored_sectors):
        if index in used:
            continue
        center = complex(
            float(stored["central_eigenvalue"]["real"]),
            float(stored["central_eigenvalue"]["imaginary"]),
        )
        choices.append((abs(center - target), index))
    _, index = min(choices)
    used.add(index)
    return index, stored_sectors[index]


def optimal_spectral_distance(left, right):
    left = np.asarray(left, dtype=np.complex128)
    right = np.asarray(right, dtype=np.complex128)
    costs = np.abs(left[:, None] - right[None, :])
    rows, columns = linear_sum_assignment(costs)
    return float(np.max(costs[rows, columns])) if len(rows) else 0.0


def tangent_analysis(ball_records):
    matrices = {name: item["midpoint"] for name, item in ball_records.items()}
    radii = {name: item["radii"] for name, item in ball_records.items()}
    defects = {name: item["defect_midpoint"] for name, item in ball_records.items()}
    defect_radii = {name: item["defect_radii"] for name, item in ball_records.items()}
    op = matrices["operational_primary"]
    op_shadow = matrices["operational_shadow"]
    val = matrices["validation_primary"]
    val_shadow = matrices["validation_shadow"]

    epsilon_t = (
        norm2(op - op_shadow) + norm2(val - val_shadow) + norm2(op - val)
        + max(la.norm(value, "fro") for value in radii.values()) + 1e-70
    )
    epsilon_sym = (
        norm2(defects["operational_primary"] - defects["operational_shadow"])
        + norm2(defects["validation_primary"] - defects["validation_shadow"])
        + norm2(defects["operational_primary"] - defects["validation_primary"])
        + max(la.norm(value, "fro") for value in defect_radii.values()) + 1e-70
    )
    symplectic_norm = norm2(defects["operational_primary"])
    symplectic_ok = bool(symplectic_norm <= 10 * epsilon_sym)

    singular = {
        name: la.svd(value, compute_uv=False, lapack_driver="gesvd")
        for name, value in matrices.items()
    }
    reciprocal = {
        name: values * values[::-1] - 1 for name, values in singular.items()
    }
    tangent_radius = max(la.norm(value, "fro") for value in radii.values())
    epsilon_reciprocal = (
        float(np.max(np.abs(reciprocal["operational_primary"]
                            - reciprocal["operational_shadow"])))
        + float(np.max(np.abs(reciprocal["validation_primary"]
                              - reciprocal["validation_shadow"])))
        + float(np.max(np.abs(reciprocal["operational_primary"]
                              - reciprocal["validation_primary"])))
        + 2 * max(1.0, singular["operational_primary"][0]) * tangent_radius
        + 10 * np.finfo(float).eps
          * max(1.0, float(np.max(np.abs(reciprocal["operational_primary"]))))
    )
    reciprocal_norm = float(np.max(np.abs(reciprocal["operational_primary"])))
    reciprocal_ok = bool(reciprocal_norm <= 10 * epsilon_reciprocal)

    determinant_moduli = {}
    determinant_log_moduli = {}
    for name, matrix in matrices.items():
        _, log_abs = np.linalg.slogdet(matrix)
        determinant_log_moduli[name] = float(log_abs)
        determinant_moduli[name] = float(np.exp(log_abs)) if log_abs < 700 else math.inf
    condition_t = float(np.linalg.cond(op))
    epsilon_det = (
        abs(determinant_log_moduli["operational_primary"]
            - determinant_log_moduli["operational_shadow"])
        + abs(determinant_log_moduli["validation_primary"]
              - determinant_log_moduli["validation_shadow"])
        + abs(determinant_log_moduli["operational_primary"]
              - determinant_log_moduli["validation_primary"])
        + 10 * np.finfo(float).eps * op.shape[0] * max(1.0, condition_t)
    )
    determinant_ok = bool(
        abs(determinant_log_moduli["operational_primary"]) <= 10 * epsilon_det
    )

    eigen = {}
    eigenvectors = None
    eigen_residual = 0.0
    eigenvector_condition = math.inf
    for name, matrix in matrices.items():
        if name == "operational_primary":
            values, vectors = la.eig(matrix)
            eigenvectors = vectors
            eigen[name] = values
            eigenvector_condition = float(np.linalg.cond(vectors))
            scale = max(1.0, norm2(matrix))
            eigen_residual = max(
                float(np.linalg.norm(matrix @ vectors[:, index]
                                     - values[index] * vectors[:, index]) / scale)
                for index in range(len(values))
            )
        else:
            eigen[name] = la.eigvals(matrix)
    variation = max(
        optimal_spectral_distance(eigen["operational_primary"], eigen[name])
        for name in ("operational_shadow", "validation_primary", "validation_shadow")
    )
    eig_floor = (
        10 * np.finfo(float).eps * max(1.0, norm2(op))
        * max(1.0, eigenvector_condition)
    )
    epsilon_eigenvalue = variation + eig_floor
    distances = np.abs(np.abs(eigen["operational_primary"]) - 1)
    if math.isfinite(epsilon_eigenvalue):
        unit = distances < 10 * epsilon_eigenvalue
        resolved = distances > 100 * epsilon_eigenvalue
        open_flags = ~(unit | resolved)
    else:
        unit = np.zeros(len(distances), dtype=bool)
        resolved = np.zeros(len(distances), dtype=bool)
        open_flags = np.ones(len(distances), dtype=bool)
    reciprocal_eigen_distance = optimal_spectral_distance(
        eigen["operational_primary"],
        1 / np.conjugate(eigen["operational_primary"]),
    )
    epsilon_singular = (
        max(float(np.max(np.abs(singular["operational_primary"] - singular[name])))
            for name in ("operational_shadow", "validation_primary", "validation_shadow"))
        + tangent_radius
        + 10 * np.finfo(float).eps * max(1.0, singular["operational_primary"][0])
    )
    canonicality_ok = bool(symplectic_ok and reciprocal_ok and determinant_ok)
    return {
        "matrices": matrices,
        "singular_arrays": singular,
        "eigen_arrays": eigen,
        "epsilon_t": epsilon_t,
        "epsilon_sym": epsilon_sym,
        "symplectic_norm": symplectic_norm,
        "symplectic_ok": symplectic_ok,
        "epsilon_reciprocal": epsilon_reciprocal,
        "reciprocal_norm": reciprocal_norm,
        "reciprocal_ok": reciprocal_ok,
        "determinant_log_moduli": determinant_log_moduli,
        "determinant_moduli": determinant_moduli,
        "epsilon_determinant_log_modulus": epsilon_det,
        "determinant_ok": determinant_ok,
        "condition_tangent": condition_t,
        "eigenvector_condition": eigenvector_condition,
        "maximum_eigenpair_residual": eigen_residual,
        "epsilon_eigenvalue": epsilon_eigenvalue,
        "epsilon_singular": epsilon_singular,
        "reciprocal_eigenvalue_distance": reciprocal_eigen_distance,
        "unit_count": int(np.sum(unit)),
        "resolved_off_unit_count": int(np.sum(resolved)),
        "open_unit_count": int(np.sum(open_flags)),
        "canonicality_ok": canonicality_ok,
    }


def serialize_float(value):
    return f"{float(value):.17e}"


def serialize_mp(value, digits=70):
    return mp.nstr(value, digits)


def serialize_complex(value):
    return {
        "real": serialize_float(np.real(value)),
        "imaginary": serialize_float(np.imag(value)),
    }


def deterministic_npz(path, arrays):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, np.asarray(arrays[name]), allow_pickle=False)
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            archive.writestr(info, buffer.getvalue())
    temporary.replace(path)


print("=" * 78)
print("BLIND FULL 1440-DIMENSIONAL 600-CELL BOUNDARY TANGENT CENSUS")
print("=" * 78)

action_sets = {
    parity: tuple(sorted(tuple(int(value) for value in action)
                         for action in models[parity]["stabilizer"]))
    for parity in ("even", "odd")
}
same_actions = action_sets["even"] == action_sets["odd"]
check(
    "both schedule parities have literally the same frozen 2T action set",
    same_actions and len(action_sets["even"]) == 24,
)

records = {}
numeric_arrays = {}
global_controls = provenance_ok and same_actions and gro.tests == gro.passed == 43
common_sector_signature = None

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing full carrier and boundary identification", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    index_data = group_and_index_data(model, state)
    geometry = prepare_geometry(model, index_data)
    mapping = tuple(gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"])

    edge_map_ok = bool(
        sorted(mapping) == list(range(30))
        and all(
            tuple(vertex + 120 for vertex in index_data["orbit_edges"][old_type][group])
            == index_data["orbit_edges"][65 + mapping[old_type]][group]
            for old_type in range(30) for group in range(24)
        )
    )
    carrier_ok = bool(
        len(model["old_edges"]) == 720
        and len(model["internal_edges"]) == 840
        and len(model["new_edges"]) == 720
        and len(model["slab"]) == 2400
        and len(geometry["triangle_records"]) == 6240
        and len(geometry["patterns"]) == 20
        and len(index_data["edge_to_index"]) == 2280
    )
    check(
        f"{parity}: literal layer shift fixes all 720 boundary edge labels",
        edge_map_ok and carrier_ok,
        f"mapping={mapping}, patterns={len(geometry['patterns'])}",
    )

    sectors, sector_control = high_precision_sector_bases(index_data)
    sector_signature = tuple(
        (sector["dimension"], mp.nstr(sector["old_central_eigenvalue"], 70),
         sector["splitter"])
        for sector in sectors
    )
    if common_sector_signature is None:
        common_sector_signature = sector_signature
    basis_ok = bool(
        sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and all(
            sector_control[key] < mp.mpf("1e-70")
            for key in (
                "maximum_orthonormal", "maximum_center", "maximum_splitter",
                "maximum_right_leakage", "maximum_conjugate_pair",
            )
        )
        and sum(60 * sector["dimension"] ** 2 for sector in sectors) == 1440
        and sector_signature == common_sector_signature
    )
    check(
        f"{parity}: seven high-precision minimal sectors exhaust 1440 real phase dimensions",
        basis_ok,
        "dims=" + str([sector["dimension"] for sector in sectors])
        + ", max residual=" + mp.nstr(max(
            value for key, value in sector_control.items() if key.startswith("maximum_")
        ), 5),
    )

    s = mp.mpf(state[0])
    kind_values = {
        "old": L0_SQUARE,
        "internal": mp.exp(s) * L0_SQUARE - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * s) * L0_SQUARE,
    }
    print(f"[{parity}] differentiating 20 local Lorentzian patterns", flush=True)
    pattern_cache, branch_control = high_precision_pattern_cache(
        geometry["patterns"], kind_values
    )
    kernels, kernel_control = assemble_full_representative_kernels(
        index_data, geometry, pattern_cache
    )
    kernel_ok = bool(
        branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and branch_control["minimum_leading_minor"] > 0
        and branch_control["minimum_argument"] > mp.mpf("1e-6")
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-70")
        and len(set(kernel_control["nonzero_entries"].values())) == 1
    )
    check(
        f"{parity}: full representative kernels pass branch, step and reality gates",
        kernel_ok,
        f"entries={kernel_control['nonzero_entries']}, "
        f"imag={mp.nstr(kernel_control['maximum_imaginary'], 5)}",
    )

    stored_used = set()
    sector_records = []
    all_reproduced = True
    all_reciprocity = True
    all_determinants = True
    all_canonicality = True
    trivial_error = math.inf

    for sector_index, sector in enumerate(sectors):
        dimension = sector["dimension"]
        print(
            f"[{parity}] sector {sector_index + 1}/7, d={dimension}: "
            f"projecting K and solving {(65*dimension)}x{(65*dimension)} J",
            flush=True,
        )
        blocks = {
            name: project_full_kernel(kernel, sector)
            for name, kernel in kernels.items()
        }
        numpy_blocks = {name: mp_to_numpy(block) for name, block in blocks.items()}
        op_block = numpy_blocks["operational_primary"]
        hermitian_norm = norm2(op_block - op_block.conj().T)
        hermitian_epsilon = (
            norm2(op_block - numpy_blocks["operational_shadow"])
            + norm2(numpy_blocks["validation_primary"]
                    - numpy_blocks["validation_shadow"])
            + norm2(op_block - numpy_blocks["validation_primary"])
            + 1e-70
        )
        reciprocity_ok = bool(hermitian_norm <= 10 * hermitian_epsilon)
        all_reciprocity &= reciprocity_ok

        stored_index, stored_sector = central_match(
            sector, rank_input["parities"][parity]["sectors"], stored_used
        )
        canonical_op = canonical_from_full_block(
            blocks["operational_primary"], dimension
        )
        current_singular = la.svd(
            mp_to_numpy(canonical_op), compute_uv=False, lapack_driver="gesvd"
        )
        stored_singular = np.asarray([
            float(value) for value in stored_sector["singular_values"]
        ])
        spectrum_error = float(np.max(
            np.abs(current_singular - stored_singular)
            / np.maximum(1.0, np.abs(stored_singular))
        ))
        reproduced = bool(spectrum_error < 2e-10)
        all_reproduced &= reproduced

        ball_records = {}
        determinant_flags = {}
        for name, block in blocks.items():
            _, det_j, tangent_ball, defect_ball = build_tangent_ball(
                block, dimension, mapping
            )
            midpoint, radii = acb_midpoint_and_radii(tangent_ball)
            defect_midpoint, defect_radii = acb_midpoint_and_radii(defect_ball)
            determinant_flags[name] = not det_j.contains(0)
            ball_records[name] = {
                "midpoint": midpoint,
                "radii": radii,
                "defect_midpoint": defect_midpoint,
                "defect_radii": defect_radii,
                "det_j": det_j,
            }
            prefix = f"{parity}_sector{sector_index}_{name}"
            numeric_arrays[f"{prefix}_tangent_midpoint"] = midpoint
            numeric_arrays[f"{prefix}_tangent_radii"] = radii
            numeric_arrays[f"{prefix}_defect_midpoint"] = defect_midpoint
            numeric_arrays[f"{prefix}_defect_radii"] = defect_radii
        determinant_ok = all(determinant_flags.values())
        all_determinants &= determinant_ok

        analysis = tangent_analysis(ball_records)
        all_canonicality &= analysis["canonicality_ok"]
        if sector["constant_overlap"] > 1 - mp.mpf("1e-70"):
            stored_tangent = np.asarray([
                [float(value) for value in row]
                for row in reduced_input["parities"][parity]["tangent_matrix"]
            ])
            current_tangent = analysis["matrices"]["operational_primary"].real
            trivial_error = float(
                np.linalg.norm(current_tangent - stored_tangent)
                / max(1.0, np.linalg.norm(current_tangent), np.linalg.norm(stored_tangent))
            )

        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "center": sector["old_central_eigenvalue"],
            "splitter": sector["splitter"],
            "constant_overlap": sector["constant_overlap"],
            "stored_sector_index": stored_index,
            "canonical_spectrum_error": spectrum_error,
            "hermitian_norm": hermitian_norm,
            "hermitian_epsilon": hermitian_epsilon,
            "reciprocity_ok": reciprocity_ok,
            "determinant_flags": determinant_flags,
            "determinant_balls": {name: str(item["det_j"])
                                  for name, item in ball_records.items()},
            "analysis": analysis,
            "maximum_tangent_radius": max(
                float(np.max(item["radii"])) for item in ball_records.values()
            ),
            "maximum_defect_radius": max(
                float(np.max(item["defect_radii"])) for item in ball_records.values()
            ),
        })

    reproduction_ok = bool(
        len(stored_used) == 7 and all_reproduced and trivial_error < 2e-8
    )
    check(
        f"{parity}: all full J spectra and the trivial tangent reproduce frozen audits",
        reproduction_ok,
        f"max J error={max(item['canonical_spectrum_error'] for item in sector_records):.3e}, "
        f"trivial tangent={trivial_error:.3e}",
    )
    check(
        f"{parity}: all projected complete Hessians are reciprocal within calibration",
        all_reciprocity,
        "max margin ratio=" + f"{max(
            item['hermitian_norm'] / max(1e-300, 10 * item['hermitian_epsilon'])
            for item in sector_records
        ):.3e}",
    )
    check(
        f"{parity}: all 28 Flint pre-Legendre determinants exclude zero",
        all_determinants,
    )
    check(
        f"{parity}: all seven tangent blocks pass calibrated canonicality gates",
        all_canonicality,
        "symplectic max ratio=" + f"{max(
            item['analysis']['symplectic_norm']
            / max(1e-300, 10 * item['analysis']['epsilon_sym'])
            for item in sector_records
        ):.3e}",
    )

    controls_ok = bool(
        edge_map_ok and carrier_ok and basis_ok and kernel_ok
        and reproduction_ok and all_reciprocity
    )
    global_controls &= controls_ok
    records[parity] = {
        "controls_ok": controls_ok,
        "all_determinants_exclude_zero": all_determinants,
        "all_canonicality_gates_pass": all_canonicality,
        "mapping": mapping,
        "basis_control": sector_control,
        "branch_control": branch_control,
        "kernel_control": kernel_control,
        "trivial_tangent_error": trivial_error,
        "sectors": sector_records,
    }


schedule_records = []
for sector_index in range(7):
    even = records["even"]["sectors"][sector_index]
    odd = records["odd"]["sectors"][sector_index]
    even_analysis = even["analysis"]
    odd_analysis = odd["analysis"]
    eigen_distance = optimal_spectral_distance(
        even_analysis["eigen_arrays"]["operational_primary"],
        odd_analysis["eigen_arrays"]["operational_primary"],
    )
    singular_distance = float(np.max(np.abs(
        even_analysis["singular_arrays"]["operational_primary"]
        - odd_analysis["singular_arrays"]["operational_primary"]
    )))
    eigen_uncertainty = (
        even_analysis["epsilon_eigenvalue"] + odd_analysis["epsilon_eigenvalue"]
    )
    singular_uncertainty = (
        even_analysis["epsilon_singular"] + odd_analysis["epsilon_singular"]
    )
    if (
        math.isfinite(eigen_uncertainty)
        and eigen_distance <= 10 * eigen_uncertainty
        and singular_distance <= 10 * singular_uncertainty
    ):
        label = "SCHEDULE_ROBUST"
    elif (
        (math.isfinite(eigen_uncertainty)
         and eigen_distance > 100 * eigen_uncertainty)
        or singular_distance > 100 * singular_uncertainty
    ):
        label = "SCHEDULE_DEPENDENT"
    else:
        label = "SCHEDULE_OPEN"
    schedule_records.append({
        "sector_index": sector_index,
        "dimension": even["dimension"],
        "eigenvalue_distance": eigen_distance,
        "eigenvalue_uncertainty": eigen_uncertainty,
        "singular_distance": singular_distance,
        "singular_uncertainty": singular_uncertainty,
        "label": label,
    })

schedule_labels = Counter(item["label"] for item in schedule_records)
if schedule_labels["SCHEDULE_DEPENDENT"]:
    schedule_outcome = "SCHEDULE_DEPENDENT"
elif schedule_labels["SCHEDULE_OPEN"]:
    schedule_outcome = "SCHEDULE_OPEN"
else:
    schedule_outcome = "SCHEDULE_ROBUST"
check(
    "the preregistered sectorwise schedule comparison assigns all seven labels",
    sum(schedule_labels.values()) == 7,
    f"global={schedule_outcome}, labels={dict(schedule_labels)}",
)

all_determinants = all(
    record["all_determinants_exclude_zero"] for record in records.values()
)
all_canonicality = all(
    record["all_canonicality_gates_pass"] for record in records.values()
)
if not global_controls:
    outcome = "FULL_BOUNDARY_TANGENT_CONTROL_FAILED"
elif not all_determinants:
    outcome = "FULL_BOUNDARY_TANGENT_RANK_OPEN"
elif not all_canonicality:
    outcome = "FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED"
else:
    outcome = "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
check(
    "the frozen hierarchy assigns exactly one full-boundary tangent outcome",
    outcome in {
        "FULL_BOUNDARY_TANGENT_CONTROL_FAILED",
        "FULL_BOUNDARY_TANGENT_RANK_OPEN",
        "FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED",
        "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED",
    },
    f"outcome={outcome}",
)


def public_sector(record):
    analysis = record["analysis"]
    eigenvalues = analysis["eigen_arrays"]["operational_primary"]
    singular = analysis["singular_arrays"]["operational_primary"]
    return {
        "sector_index": record["sector_index"],
        "dimension": record["dimension"],
        "real_phase_dimension_with_multiplicity": 60 * record["dimension"] ** 2,
        "central_eigenvalue": {
            "real": serialize_mp(mp.re(record["center"])),
            "imaginary": serialize_mp(mp.im(record["center"])),
        },
        "splitter_group_index": record["splitter"],
        "constant_overlap": serialize_mp(record["constant_overlap"]),
        "stored_sector_index": record["stored_sector_index"],
        "canonical_spectrum_error": serialize_float(record["canonical_spectrum_error"]),
        "hermitian_norm": serialize_float(record["hermitian_norm"]),
        "hermitian_epsilon": serialize_float(record["hermitian_epsilon"]),
        "reciprocity_ok": record["reciprocity_ok"],
        "pre_legendre_determinants_exclude_zero": record["determinant_flags"],
        "pre_legendre_determinant_balls": record["determinant_balls"],
        "maximum_tangent_ball_radius": serialize_float(record["maximum_tangent_radius"]),
        "maximum_defect_ball_radius": serialize_float(record["maximum_defect_radius"]),
        "epsilon_t": serialize_float(analysis["epsilon_t"]),
        "symplectic_norm": serialize_float(analysis["symplectic_norm"]),
        "epsilon_sym": serialize_float(analysis["epsilon_sym"]),
        "symplectic_ok": analysis["symplectic_ok"],
        "reciprocal_singular_norm": serialize_float(analysis["reciprocal_norm"]),
        "epsilon_reciprocal_singular": serialize_float(
            analysis["epsilon_reciprocal"]
        ),
        "reciprocal_singular_ok": analysis["reciprocal_ok"],
        "determinant_modulus": serialize_float(
            analysis["determinant_moduli"]["operational_primary"]
        ),
        "determinant_log_modulus": serialize_float(
            analysis["determinant_log_moduli"]["operational_primary"]
        ),
        "epsilon_determinant_log_modulus": serialize_float(
            analysis["epsilon_determinant_log_modulus"]
        ),
        "determinant_ok": analysis["determinant_ok"],
        "tangent_condition": serialize_float(analysis["condition_tangent"]),
        "eigenvector_condition": serialize_float(analysis["eigenvector_condition"]),
        "maximum_eigenpair_residual": serialize_float(
            analysis["maximum_eigenpair_residual"]
        ),
        "epsilon_eigenvalue": serialize_float(analysis["epsilon_eigenvalue"]),
        "epsilon_singular": serialize_float(analysis["epsilon_singular"]),
        "reciprocal_conjugate_eigenvalue_distance": serialize_float(
            analysis["reciprocal_eigenvalue_distance"]
        ),
        "spectral_radius": serialize_float(np.max(np.abs(eigenvalues))),
        "minimum_eigenvalue_modulus": serialize_float(np.min(np.abs(eigenvalues))),
        "unit_consistent_count_minimal": analysis["unit_count"],
        "resolved_off_unit_count_minimal": analysis["resolved_off_unit_count"],
        "open_unit_count_minimal": analysis["open_unit_count"],
        "canonicality_ok": analysis["canonicality_ok"],
        "singular_values": [serialize_float(value) for value in singular],
        "eigenvalues": [serialize_complex(value) for value in eigenvalues],
    }


deterministic_npz(NUMERIC_OUTPUT, numeric_arrays)
numeric_hash = sha256(NUMERIC_OUTPUT)
artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "numeric_archive": NUMERIC_OUTPUT.name,
    "numeric_archive_sha256": numeric_hash,
    "numeric_archive_arrays": len(numeric_arrays),
    "full_boundary_phase_dimension": 1440,
    "full_boundary_configuration_dimension": 720,
    "minimal_sector_dimensions": [sector["dimension"] for sector in records["even"]["sectors"]],
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "planck_target_parsed": False,
    "particle_target_parsed": False,
    "outcome": outcome,
    "schedule_outcome": schedule_outcome,
    "schedule_sectors": [
        {
            key: (serialize_float(value) if isinstance(value, float) else value)
            for key, value in item.items()
        }
        for item in schedule_records
    ],
    "parities": {
        parity: {
            "controls_ok": record["controls_ok"],
            "all_determinants_exclude_zero": record["all_determinants_exclude_zero"],
            "all_canonicality_gates_pass": record["all_canonicality_gates_pass"],
            "old_to_final_orbit_map": list(record["mapping"]),
            "trivial_tangent_error": serialize_float(record["trivial_tangent_error"]),
            "basis_control": {
                key: ([int(value) for value in item] if isinstance(item, list)
                      else serialize_mp(item))
                for key, item in record["basis_control"].items()
            },
            "branch_control": {
                "minimum_leading_minor": serialize_mp(
                    record["branch_control"]["minimum_leading_minor"]
                ),
                "minimum_argument": serialize_mp(
                    record["branch_control"]["minimum_argument"]
                ),
                "maximum_raw_angle_or_derivative_imaginary": serialize_mp(
                    record["branch_control"]["maximum_raw_angle_or_derivative_imaginary"]
                ),
                "maximum_cross": serialize_mp(record["branch_control"]["maximum_cross"]),
                "maximum_proxy": serialize_mp(record["branch_control"]["maximum_proxy"]),
                "entry_pass": record["branch_control"]["entry_pass"],
                "base_negative_counts": dict(record["branch_control"]["base_negative_counts"]),
                "displaced_negative_counts": dict(
                    record["branch_control"]["displaced_negative_counts"]
                ),
            },
            "kernel_control": {
                "nonzero_entries": record["kernel_control"]["nonzero_entries"],
                "maximum_imaginary": serialize_mp(
                    record["kernel_control"]["maximum_imaginary"]
                ),
            },
            "sectors": [public_sector(item) for item in record["sectors"]],
        }
        for parity, record in records.items()
    },
    "classification": {
        "canonical_map": (
            "DERIVED COMPUTATIONAL"
            if outcome == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
            else "OPEN_OR_FAILED"
        ),
        "one_step_spectrum": "DERIVED COMPUTATIONAL",
        "physical_stability": "OPEN",
        "gauge_invariant_modes": "OPEN",
        "multi_tick_dynamics": "OPEN",
        "dispersion_relation": "OPEN",
        "limiting_speed": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Schedule: {schedule_outcome} {dict(schedule_labels)}")
for parity, record in records.items():
    replicated = Counter()
    radii = []
    for sector in record["sectors"]:
        dimension = sector["dimension"]
        analysis = sector["analysis"]
        replicated["unit"] += dimension * analysis["unit_count"]
        replicated["off"] += dimension * analysis["resolved_off_unit_count"]
        replicated["open"] += dimension * analysis["open_unit_count"]
        radii.append(float(np.max(np.abs(
            analysis["eigen_arrays"]["operational_primary"]
        ))))
    print(
        f"{parity}: replicated unit/off/open="
        f"{replicated['unit']}/{replicated['off']}/{replicated['open']}, "
        f"sector spectral radii={[f'{value:.6g}' for value in radii]}"
    )
print(f"Numeric archive: {NUMERIC_OUTPUT.name} sha256={numeric_hash}")
print(f"Results: {passed}/{tests} tests passed.")
sys.exit(0 if passed == tests else 1)

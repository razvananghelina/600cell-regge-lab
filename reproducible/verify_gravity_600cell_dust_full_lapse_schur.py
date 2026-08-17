#!/usr/bin/env python3
"""High-precision 2T Schur audit of all 120 vertex-lapse directions.

Prior-art commit: 58f14e1.
Protocol commit: 15a0699.
No continuum spectrum, wave speed, or experimental target is evaluated.
"""

from collections import Counter, defaultdict
import contextlib
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import sys

from flint import acb, acb_mat, ctx
import mpmath as mp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
RANK_INPUT = HERE / "gravity_600cell_dust_full_anisotropic_legendre_rank.json"
RANK_VERIFIER = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
OUTPUT = HERE / "gravity_600cell_dust_full_lapse_schur.json"
PRIOR_ART_COMMIT = "58f14e1"
PROTOCOL_COMMIT = "15a0699"
EXPECTED_HASHES = {
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "rank_artifact": "7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226",
    "rank_verifier": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
}
DPS = 100
BALL_DPS = 80
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-20"),
    "operational_shadow": mp.mpf("1e-15"),
    "validation_primary": mp.mpf("3e-20"),
    "validation_shadow": mp.mpf("3e-15"),
}
I = mp.mpc(0, 1)
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


hashes = {
    "tick": sha256(TICK_INPUT),
    "rank_artifact": sha256(RANK_INPUT),
    "rank_verifier": sha256(RANK_VERIFIER),
}
tick = json.loads(TICK_INPUT.read_text())
rank_input = json.loads(RANK_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and rank_input["outcome"] == "FULL_CANONICAL_LEGENDRE_REGULAR"
    and rank_input["passed"] == rank_input["tests"] == 18
    and all(
        record["outcome"] == "FULL_CANONICAL_LEGENDRE_REGULAR"
        and record["full_resolved_rank"] == 1560
        and record["full_error_consistent_nullity"] == 0
        and record["full_numerically_open_count"] == 0
        for record in rank_input["parities"].values()
    )
)
check("the frozen full-rank inputs have exact committed provenance", provenance_ok, str(hashes))


spec = importlib.util.spec_from_file_location(
    "full_anisotropic_rank_for_lapse_schur", RANK_VERIFIER
)
full = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = full
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(full)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check(
    "the imported complete 2280-edge verifier retains all 18 controls",
    full.tests == full.passed == 18,
)


def mp_max_abs(matrix):
    return max((abs(value) for value in matrix), default=mp.mpf(0))


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
    maximum_orthonormal = mp.mpf(0)
    maximum_center_residual = mp.mpf(0)
    maximum_splitter_residual = mp.mpf(0)
    maximum_right_leakage = mp.mpf(0)

    for cluster in center_clusters:
        isotypic = mp_submatrix(center_vectors, range(24), cluster)
        dimension = math.isqrt(len(cluster))
        if dimension**2 != len(cluster):
            raise RuntimeError("high-precision central eigenspace is not square")
        splitter = None
        selected = None
        splitter_value = None
        splitter_matrix = None
        if dimension == 1:
            selected = isotypic
            splitter_value = mp.mpf(0)
            splitter_matrix = mp.matrix(24, 24)
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
            raise RuntimeError("high-precision deterministic splitter failed")

        orthonormal = selected.H * selected - mp.eye(dimension)
        center_value = mp.fsum(center_values[index] for index in cluster) / len(cluster)
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
        right_leakage = mp.mpf(0)
        projector = selected * selected.H
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
        maximum_orthonormal = max(maximum_orthonormal, mp_frobenius(orthonormal))
        maximum_center_residual = max(
            maximum_center_residual, mp_frobenius(center_residual)
        )
        maximum_splitter_residual = max(
            maximum_splitter_residual, mp_frobenius(splitter_residual)
        )
        maximum_right_leakage = max(maximum_right_leakage, right_leakage)
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
        "maximum_orthonormal_error": maximum_orthonormal,
        "maximum_center_residual": maximum_center_residual,
        "maximum_splitter_residual": maximum_splitter_residual,
        "maximum_right_invariance_leakage": maximum_right_leakage,
    }


def high_precision_pattern_cache(patterns, kind_values):
    cache = {}
    minimum_minor = mp.inf
    minimum_argument = mp.inf
    maximum_imaginary = mp.mpf(0)
    base_negative_counts = Counter()
    displaced_negative_counts = Counter()
    for pattern in sorted(patterns):
        values = [kind_values[kind] for kind in pattern]
        base_angles, base_branch = full.angle_record(values)
        base_negative_counts[base_branch["negative_directions"]] += patterns[pattern]
        minimum_minor = min(minimum_minor, base_branch["minimum_leading_minor"])
        minimum_argument = min(minimum_argument, base_branch["minimum_argument"])
        maximum_imaginary = max(
            maximum_imaginary,
            *(abs(mp.im(value)) for value in base_angles),
        )
        derivatives = {}
        for name, step in DERIVATIVE_STEPS.items():
            matrix = mp.matrix(10, 10)
            for column in range(10):
                plus = list(values)
                minus = list(values)
                plus[column] *= mp.exp(step)
                minus[column] *= mp.exp(-step)
                plus_angles, plus_branch = full.angle_record(plus)
                minus_angles, minus_branch = full.angle_record(minus)
                for branch in (plus_branch, minus_branch):
                    displaced_negative_counts[branch["negative_directions"]] += 1
                    minimum_minor = min(
                        minimum_minor, branch["minimum_leading_minor"]
                    )
                    minimum_argument = min(
                        minimum_argument, branch["minimum_argument"]
                    )
                for row in range(10):
                    matrix[row, column] = (
                        plus_angles[row] - minus_angles[row]
                    ) / (2 * step)
                    maximum_imaginary = max(
                        maximum_imaginary, abs(mp.im(matrix[row, column]))
                    )
            derivatives[name] = matrix
        cache[pattern] = {
            "base_angles": tuple(base_angles),
            "derivatives": derivatives,
        }
    return cache, {
        "minimum_leading_minor": minimum_minor,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
        "base_negative_counts": base_negative_counts,
        "displaced_negative_counts": displaced_negative_counts,
    }


ROW_TYPES = tuple(range(30, 65)) + tuple(range(0, 30))
COLUMN_TYPES = tuple(range(30, 65)) + tuple(range(65, 95))
ROW_POSITION = {kind: index for index, kind in enumerate(ROW_TYPES)}
COLUMN_POSITION = {kind: index for index, kind in enumerate(COLUMN_TYPES)}


def assemble_representative_kernels(index_data, geometry, pattern_cache):
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

    kernels = {name: defaultdict(lambda: mp.mpc(0)) for name in DERIVATIVE_STEPS}

    def canonical_location(global_index, positions):
        orbit_type, group = divmod(global_index, 24)
        if orbit_type not in positions:
            return None
        return positions[orbit_type], group

    def add_all(row_index, column_index, value):
        row_location = canonical_location(row_index, ROW_POSITION)
        column_location = canonical_location(column_index, COLUMN_POSITION)
        if (
            row_location is None or column_location is None
            or row_location[1] != identity
        ):
            return
        row_position = row_location[0]
        column_position, column_group = column_location
        sign = -1 if row_position >= 35 else 1
        key = (row_position, column_position, column_group)
        for kernel in kernels.values():
            kernel[key] += sign * value

    def add_angle(name, row_index, column_index, value):
        row_location = canonical_location(row_index, ROW_POSITION)
        column_location = canonical_location(column_index, COLUMN_POSITION)
        if (
            row_location is None or column_location is None
            or row_location[1] != identity
        ):
            return
        row_position = row_location[0]
        column_position, column_group = column_location
        sign = -1 if row_position >= 35 else 1
        kernels[name][(row_position, column_position, column_group)] += sign * value

    triangle_local = []
    for triangle_index, record in enumerate(geometry["triangle_records"]):
        indices = record["indices"]
        values = [signed_base[index] for index in indices]
        _, area_gradient, area_hessian = full.area_data(values)
        coefficient = -I * curvature[triangle_index]
        for row, row_index in enumerate(indices):
            for column, column_index in enumerate(indices):
                add_all(
                    row_index, column_index,
                    coefficient * area_hessian[row][column],
                )
        triangle_local.append(tuple(area_gradient))

    dust_hessian = -(2 * mp.pi * full.MASS / 120) * mp.sqrt(index_data["rho"])
    for pole_index in geometry["pole_indices"]:
        add_all(pole_index, pole_index, dust_hessian)

    for simplex in geometry["simplex_records"]:
        simplex_indices = simplex["indices"]
        derivatives = pattern_cache[simplex["pattern"]]["derivatives"]
        for hinge_index, triangle in enumerate(simplex["hinge_triangles"]):
            triangle_indices = geometry["triangle_records"][triangle]["indices"]
            area_gradient = triangle_local[triangle]
            for row, row_index in enumerate(triangle_indices):
                row_location = canonical_location(row_index, ROW_POSITION)
                if row_location is None or row_location[1] != identity:
                    continue
                for column, column_index in enumerate(simplex_indices):
                    for name in DERIVATIVE_STEPS:
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
        "nonzero_kernel_entries": {
            name: len(kernel) for name, kernel in kernels.items()
        },
        "maximum_imaginary": maximum_imaginary,
    }


def project_kernel(kernel, sector):
    dimension = sector["dimension"]
    block = mp.matrix(65 * dimension, 65 * dimension)
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


def expanded_indices(positions, dimension):
    return [
        position * dimension + component
        for position in positions for component in range(dimension)
    ]


def ball_schur(block, dimension, weak_positions):
    strong_positions = [index for index in range(65) if index not in weak_positions]
    strong = expanded_indices(strong_positions, dimension)
    weak = expanded_indices(weak_positions, dimension)
    a_mp = mp_submatrix(block, strong, strong)
    b_mp = mp_submatrix(block, strong, weak)
    c_mp = mp_submatrix(block, weak, strong)
    d_mp = mp_submatrix(block, weak, weak)
    a = mp_matrix_to_acb(a_mp)
    b = mp_matrix_to_acb(b_mp)
    c = mp_matrix_to_acb(c_mp)
    d = mp_matrix_to_acb(d_mp)
    det_a = a.det()
    solve = a.solve(b)
    schur = d - c * solve
    det_schur = schur.det()
    a_midpoint, a_radii = acb_midpoint_and_radii(a)
    schur_midpoint, schur_radii = acb_midpoint_and_radii(schur)
    solve_midpoint, solve_radii = acb_midpoint_and_radii(solve)
    return {
        "a_midpoint": a_midpoint,
        "a_radii": a_radii,
        "schur_midpoint": schur_midpoint,
        "schur_radii": schur_radii,
        "solve_midpoint": solve_midpoint,
        "solve_radii": solve_radii,
        "det_a": det_a,
        "det_schur": det_schur,
        "strong_positions": strong_positions,
        "weak_positions": weak_positions,
    }


def calibrated_singular_classification(records, field):
    matrices = {name: record[f"{field}_midpoint"] for name, record in records.items()}
    radii = {name: record[f"{field}_radii"] for name, record in records.items()}
    operational = matrices["operational_primary"]
    epsilon_step = (
        la.norm(operational - matrices["operational_shadow"], 2)
        + la.norm(matrices["validation_primary"] - matrices["validation_shadow"], 2)
        + la.norm(operational - matrices["validation_primary"], 2)
    )
    epsilon_ball = max(la.norm(value, "fro") for value in radii.values())
    left, singular, right_h = la.svd(
        operational, full_matrices=False, lapack_driver="gesdd"
    )
    singular_gesvd = la.svd(
        operational, compute_uv=False, lapack_driver="gesvd"
    )
    residuals = []
    for index, value in enumerate(singular):
        u = left[:, index]
        v = right_h[index, :].conj()
        residuals.append(float(
            np.linalg.norm(operational @ v - value * u)
            + np.linalg.norm(operational.conj().T @ u - value * v)
        ))
    epsilon_svd = max(residuals, default=0.0) + float(
        np.max(np.abs(singular - singular_gesvd))
    )
    epsilon_global = epsilon_step + epsilon_ball + epsilon_svd
    other_singular = {
        name: la.svd(matrix, compute_uv=False, lapack_driver="gesvd")
        for name, matrix in matrices.items() if name != "operational_primary"
    }
    resolved = singular > 100 * epsilon_global
    zero = np.array([
        max([singular[index]] + [
            values[index] for values in other_singular.values()
        ]) < 10 * epsilon_global
        for index in range(len(singular))
    ], dtype=bool)
    open_flags = ~(resolved | zero)
    return {
        "singular_values": singular,
        "left": left,
        "right_h": right_h,
        "epsilon_step": epsilon_step,
        "epsilon_ball": epsilon_ball,
        "epsilon_svd": epsilon_svd,
        "epsilon_global": epsilon_global,
        "resolved": resolved,
        "zero": zero,
        "open": open_flags,
    }


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


def geometric_lapse_matrix(index_data, weak_positions):
    diagonal = mp.exp(mp.mpf(_CURRENT_STATE[0])) * full.L0_SQUARE - index_data["rho"]
    coefficient = float(-index_data["rho"] / diagonal)
    weak_orbit_types = [30 + position for position in weak_positions]
    pole_edges = []
    for orbit_type in weak_orbit_types:
        pole_edges.extend(index_data["orbit_edges"][orbit_type])
    pole_to_column = {edge: index for index, edge in enumerate(pole_edges)}
    rows = []
    columns = []
    values = []
    for pole, column in pole_to_column.items():
        pole_index = index_data["edge_to_index"][pole]
        pole_type, pole_group = divmod(pole_index, 24)
        canonical_pole_row = COLUMN_POSITION[pole_type] * 24 + pole_group
        rows.append(canonical_pole_row)
        columns.append(column)
        values.append(1.0)
        new_vertex = pole[1]
        for edge, global_index in index_data["edge_to_index"].items():
            if (
                index_data["edge_kind"][global_index] == "internal"
                and edge[1] == new_vertex
            ):
                orbit_type, group = divmod(global_index, 24)
                rows.append(COLUMN_POSITION[orbit_type] * 24 + group)
                columns.append(column)
                values.append(coefficient)
    matrix = sp.csr_matrix((values, (rows, columns)), shape=(1560, 120))
    expected = np.zeros(1560)
    for position, orbit_type in enumerate(COLUMN_TYPES):
        kind = index_data["edge_kind"][24 * orbit_type]
        if kind == "internal":
            expected[24 * position:24 * (position + 1)] = coefficient
        elif kind == "pole":
            expected[24 * position:24 * (position + 1)] = 1
    summed = np.asarray(matrix.sum(axis=1)).ravel()
    collective_error = float(
        np.linalg.norm(summed / np.linalg.norm(summed) - expected / np.linalg.norm(expected))
    )
    return matrix, {
        "rank": int(np.linalg.matrix_rank(matrix.toarray())),
        "collective_normalized_error": collective_error,
        "coefficient": coefficient,
    }


def projected_geometric_lapse(geometry_matrix, sector, weak_positions):
    dimension = sector["dimension"]
    basis = np.array([
        [complex(float(mp.re(sector["basis"][row, column])),
                 float(mp.im(sector["basis"][row, column])))
         for column in range(dimension)]
        for row in range(24)
    ])
    full_basis = np.kron(np.eye(65), basis)
    pole_basis = np.kron(np.eye(5), basis)
    projected = full_basis.conj().T @ (geometry_matrix @ pole_basis)
    strong_positions = [index for index in range(65) if index not in weak_positions]
    order = expanded_indices(strong_positions + weak_positions, dimension)
    return projected[order, :]


def subspace_distance(left, right):
    q_left = la.qr(left, mode="economic")[0]
    q_right = la.qr(right, mode="economic")[0]
    overlap = la.svd(q_left.conj().T @ q_right, compute_uv=False)
    minimum = min(1.0, max(0.0, float(np.min(overlap))))
    return math.sqrt(max(0.0, 1 - minimum**2)), minimum


def subspace_label(distance):
    if distance < 1e-8:
        return "IDENTIFIED"
    if distance > 1e-4:
        return "SEPARATED"
    return "NUMERICALLY_OPEN"


def serialize_float(value):
    return f"{float(value):.17e}"


def serialize_mp(value, digits=70):
    return mp.nstr(value, digits)


def serialize_complex(value):
    return {
        "real": serialize_float(np.real(value)),
        "imaginary": serialize_float(np.imag(value)),
    }


def serialize_ball(value):
    return str(value)


print("=" * 78)
print("FULL 120-DIMENSIONAL VERTEX-LAPSE SCHUR AUDIT")
print("=" * 78)

ctx.dps = BALL_DPS
records = {}
global_controls = provenance_ok and full.tests == full.passed == 18

for parity in ("even", "odd"):
    print(f"[{parity}] preparing high-precision representative kernel", flush=True)
    model = full.models[parity]
    _CURRENT_STATE = tick["solutions"][parity]["state"]
    index_data = full.group_and_index_data(model, _CURRENT_STATE)
    geometry = full.prepare_geometry(model, index_data)
    weak_positions = [
        position for position in range(35)
        if index_data["edge_kind"][24 * (30 + position)] == "pole"
    ]
    pole_orbit_control = bool(
        len(weak_positions) == 5
        and len(geometry["pole_indices"]) == 120
        and all(
            len(index_data["global_orbits"][30 + position]) == 24
            for position in weak_positions
        )
        and all(
            index_data["edge_kind"][24 * (30 + position)] == "internal"
            for position in range(35) if position not in weak_positions
        )
    )
    check(
        f"{parity}: five geometry-selected free pole orbits define the 120-dimensional partition",
        pole_orbit_control,
        f"weak positions={weak_positions}",
    )

    sectors, sector_control = high_precision_sector_bases(index_data)
    basis_ok = bool(
        sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and sector_control["maximum_orthonormal_error"] < mp.mpf("1e-70")
        and sector_control["maximum_center_residual"] < mp.mpf("1e-70")
        and sector_control["maximum_splitter_residual"] < mp.mpf("1e-70")
        and sector_control["maximum_right_invariance_leakage"] < mp.mpf("1e-70")
        and sum(item["constant_overlap"] > 1 - mp.mpf("1e-70") for item in sectors) == 1
    )
    check(
        f"{parity}: the target-independent high-precision 2T basis passes all residual gates",
        basis_ok,
        "max residuals=" + ", ".join(
            f"{key}:{mp.nstr(value, 4)}" for key, value in sector_control.items()
            if key.startswith("maximum_")
        ),
    )

    s = mp.mpf(_CURRENT_STATE[0])
    kind_values = {
        "old": full.L0_SQUARE,
        "internal": mp.exp(s) * full.L0_SQUARE - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * s) * full.L0_SQUARE,
    }
    pattern_cache, branch = high_precision_pattern_cache(
        geometry["patterns"], kind_values
    )
    kernels, kernel_control = assemble_representative_kernels(
        index_data, geometry, pattern_cache
    )
    branch_kernel_ok = bool(
        branch["base_negative_counts"] == Counter({1: 2400})
        and branch["displaced_negative_counts"] == Counter({1: 1600})
        and branch["minimum_leading_minor"] > 0
        and branch["minimum_argument"] > mp.mpf("1e-6")
        and branch["maximum_imaginary"] < mp.mpf("1e-70")
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-70")
        and len(set(kernel_control["nonzero_kernel_entries"].values())) == 1
    )
    check(
        f"{parity}: all representative-kernel branch and reality controls pass",
        branch_kernel_ok,
        f"entries={kernel_control['nonzero_kernel_entries']}, "
        f"imag={mp.nstr(kernel_control['maximum_imaginary'], 5)}",
    )

    blocks_by_sector = []
    stored_used = set()
    spectrum_errors = []
    for sector in sectors:
        blocks = {
            name: project_kernel(kernel, sector) for name, kernel in kernels.items()
        }
        stored_index, stored_sector = central_match(
            sector, rank_input["parities"][parity]["sectors"], stored_used
        )
        singular = la.svd(
            mp_to_numpy(blocks["operational_primary"]),
            compute_uv=False, lapack_driver="gesvd",
        )
        stored_singular = np.asarray([
            float(value) for value in stored_sector["singular_values"]
        ])
        spectrum_error = float(np.max(
            np.abs(singular - stored_singular)
            / np.maximum(1.0, np.abs(stored_singular))
        ))
        spectrum_errors.append(spectrum_error)
        blocks_by_sector.append({
            "sector": sector,
            "blocks": blocks,
            "stored_index": stored_index,
            "spectrum_error": spectrum_error,
        })
    spectrum_control_ok = bool(
        len(stored_used) == 7 and max(spectrum_errors) < 2e-10
    )
    check(
        f"{parity}: high-precision representative blocks reproduce all stored spectra",
        spectrum_control_ok,
        f"maximum normalized discrepancy={max(spectrum_errors):.3e}",
    )

    geometry_lapse, geometry_control = geometric_lapse_matrix(
        index_data, weak_positions
    )
    geometry_ok = bool(
        geometry_control["rank"] == 120
        and geometry_control["collective_normalized_error"] < 1e-14
    )
    check(
        f"{parity}: the 120 frozen vertex-lapse columns are independent and sum to the collective lapse",
        geometry_ok,
        f"rank={geometry_control['rank']}, "
        f"sum error={geometry_control['collective_normalized_error']:.3e}",
    )

    sector_records = []
    all_strong_resolved = True
    any_schur_open = False
    any_schur_zero = False
    all_schur_resolved = True
    all_determinants_exclude_zero = True
    for item in blocks_by_sector:
        sector = item["sector"]
        dimension = sector["dimension"]
        print(
            f"[{parity}] Flint solves for d={dimension}, "
            f"A={60*dimension}, S={5*dimension}",
            flush=True,
        )
        ball_records = {
            name: ball_schur(block, dimension, weak_positions)
            for name, block in item["blocks"].items()
        }
        determinant_a = {
            name: not record["det_a"].contains(0)
            for name, record in ball_records.items()
        }
        determinant_s = {
            name: not record["det_schur"].contains(0)
            for name, record in ball_records.items()
        }
        strong_class = calibrated_singular_classification(ball_records, "a")
        schur_class = calibrated_singular_classification(ball_records, "schur")
        strong_resolved = bool(
            all(determinant_a.values())
            and np.all(strong_class["resolved"])
            and not np.any(strong_class["open"] | strong_class["zero"])
        )
        schur_resolved = bool(
            all(determinant_s.values())
            and np.all(schur_class["resolved"])
            and not np.any(schur_class["open"] | schur_class["zero"])
        )
        all_strong_resolved &= strong_resolved
        all_schur_resolved &= schur_resolved
        any_schur_open |= bool(np.any(schur_class["open"]))
        any_schur_zero |= bool(np.any(schur_class["zero"]))
        all_determinants_exclude_zero &= all(determinant_s.values())

        operational = ball_records["operational_primary"]
        strong_positions = operational["strong_positions"]
        order = expanded_indices(strong_positions + weak_positions, dimension)
        full_block = mp_to_numpy(item["blocks"]["operational_primary"])
        ordered_block = full_block[np.ix_(order, order)]
        solve_midpoint = operational["solve_midpoint"]
        lift = np.vstack((-solve_midpoint, np.eye(5 * dimension)))
        geometric = projected_geometric_lapse(
            geometry_lapse, sector, weak_positions
        )
        _, _, right_h = la.svd(
            ordered_block, full_matrices=False, lapack_driver="gesdd"
        )
        weak_singular = right_h[-5 * dimension:, :].conj().T
        lift_geometry_distance, lift_geometry_overlap = subspace_distance(
            lift, geometric
        )
        lift_singular_distance, lift_singular_overlap = subspace_distance(
            lift, weak_singular
        )
        geometry_singular_distance, geometry_singular_overlap = subspace_distance(
            geometric, weak_singular
        )
        subspaces = {
            "canonical_vs_geometric": {
                "distance": lift_geometry_distance,
                "minimum_overlap": lift_geometry_overlap,
                "label": subspace_label(lift_geometry_distance),
            },
            "canonical_vs_weak_singular": {
                "distance": lift_singular_distance,
                "minimum_overlap": lift_singular_overlap,
                "label": subspace_label(lift_singular_distance),
            },
            "geometric_vs_weak_singular": {
                "distance": geometry_singular_distance,
                "minimum_overlap": geometry_singular_overlap,
                "label": subspace_label(geometry_singular_distance),
            },
        }

        sector_records.append({
            "dimension": dimension,
            "old_central_eigenvalue": sector["old_central_eigenvalue"],
            "splitter": sector["splitter"],
            "constant_overlap": sector["constant_overlap"],
            "stored_sector_index": item["stored_index"],
            "stored_spectrum_error": item["spectrum_error"],
            "determinant_a_excludes_zero": determinant_a,
            "determinant_schur_excludes_zero": determinant_s,
            "determinant_a_balls": {
                name: record["det_a"] for name, record in ball_records.items()
            },
            "determinant_schur_balls": {
                name: record["det_schur"] for name, record in ball_records.items()
            },
            "strong": strong_class,
            "schur": schur_class,
            "schur_midpoint": operational["schur_midpoint"],
            "schur_radii": operational["schur_radii"],
            "solve_radii": operational["solve_radii"],
            "strong_resolved": strong_resolved,
            "schur_resolved": schur_resolved,
            "subspaces": subspaces,
        })

    strong_ok = bool(all_strong_resolved)
    check(
        f"{parity}: all seven 60d strong blocks are precision-resolved and invertible",
        strong_ok,
        "minimum margins=" + ", ".join(
            f"{record['strong']['singular_values'][-1]/(100*record['strong']['epsilon_global']):.3e}"
            for record in sector_records
        ),
    )

    if not (
        pole_orbit_control and basis_ok and branch_kernel_ok
        and spectrum_control_ok and geometry_ok
    ):
        outcome = "FULL_LAPSE_SCHUR_ASSEMBLY_CONTROL_FAILED"
    elif not strong_ok:
        outcome = "FULL_LAPSE_SCHUR_STRONG_BLOCK_OPEN"
    elif any_schur_open:
        outcome = "FULL_LAPSE_SCHUR_NUMERICALLY_OPEN"
    elif all_schur_resolved and all_determinants_exclude_zero:
        outcome = "FULL_LAPSE_SCHUR_REGULAR"
    elif any_schur_zero:
        outcome = "FULL_LAPSE_SCHUR_DEGENERATE"
    else:
        outcome = "FULL_LAPSE_SCHUR_NUMERICALLY_OPEN"
    check(
        f"{parity}: the preregistered hierarchy assigns the full pole-Schur outcome",
        outcome in {
            "FULL_LAPSE_SCHUR_ASSEMBLY_CONTROL_FAILED",
            "FULL_LAPSE_SCHUR_STRONG_BLOCK_OPEN",
            "FULL_LAPSE_SCHUR_NUMERICALLY_OPEN",
            "FULL_LAPSE_SCHUR_REGULAR",
            "FULL_LAPSE_SCHUR_DEGENERATE",
        },
        f"outcome={outcome}, resolved="
        f"{sum(np.sum(item['schur']['resolved'])*item['dimension'] for item in sector_records)}/120, "
        f"zero={sum(np.sum(item['schur']['zero'])*item['dimension'] for item in sector_records)}, "
        f"open={sum(np.sum(item['schur']['open'])*item['dimension'] for item in sector_records)}",
    )

    subspace_labels = Counter(
        comparison["label"]
        for record in sector_records
        for comparison in record["subspaces"].values()
    )
    check(
        f"{parity}: all 21 preregistered subspace comparisons receive a mechanical label",
        sum(subspace_labels.values()) == 21
        and set(subspace_labels) <= {"IDENTIFIED", "SEPARATED", "NUMERICALLY_OPEN"},
        str(dict(subspace_labels)),
    )

    controls_ok = bool(
        pole_orbit_control and basis_ok and branch_kernel_ok
        and spectrum_control_ok and geometry_ok
    )
    global_controls &= controls_ok
    records[parity] = {
        "outcome": outcome,
        "controls_pass": controls_ok,
        "weak_orbit_positions": weak_positions,
        "geometry_lapse": geometry_control,
        "basis_control": sector_control,
        "branch_control": branch,
        "kernel_control": kernel_control,
        "maximum_stored_spectrum_error": max(spectrum_errors),
        "resolved_schur_rank": int(sum(
            np.sum(item["schur"]["resolved"]) * item["dimension"]
            for item in sector_records
        )),
        "schur_zero_count": int(sum(
            np.sum(item["schur"]["zero"]) * item["dimension"]
            for item in sector_records
        )),
        "schur_open_count": int(sum(
            np.sum(item["schur"]["open"]) * item["dimension"]
            for item in sector_records
        )),
        "subspace_labels": dict(subspace_labels),
        "sectors": sector_records,
    }


if not global_controls:
    combined_outcome = "FULL_LAPSE_SCHUR_ASSEMBLY_CONTROL_FAILED"
elif any(
    record["outcome"] == "FULL_LAPSE_SCHUR_STRONG_BLOCK_OPEN"
    for record in records.values()
):
    combined_outcome = "FULL_LAPSE_SCHUR_STRONG_BLOCK_OPEN"
elif any(
    record["outcome"] == "FULL_LAPSE_SCHUR_NUMERICALLY_OPEN"
    for record in records.values()
):
    combined_outcome = "FULL_LAPSE_SCHUR_NUMERICALLY_OPEN"
elif all(
    record["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
    for record in records.values()
):
    combined_outcome = "FULL_LAPSE_SCHUR_REGULAR"
elif any(
    record["outcome"] == "FULL_LAPSE_SCHUR_DEGENERATE"
    for record in records.values()
):
    combined_outcome = "FULL_LAPSE_SCHUR_DEGENERATE"
else:
    combined_outcome = "FULL_LAPSE_SCHUR_NUMERICALLY_OPEN"


def serialize_sector(record):
    return {
        "irrep_dimension": record["dimension"],
        "old_central_eigenvalue": {
            "real": serialize_mp(mp.re(record["old_central_eigenvalue"])),
            "imaginary": serialize_mp(mp.im(record["old_central_eigenvalue"])),
        },
        "splitter_group_index": record["splitter"],
        "constant_overlap": serialize_mp(record["constant_overlap"]),
        "stored_sector_index": record["stored_sector_index"],
        "stored_spectrum_error": serialize_float(record["stored_spectrum_error"]),
        "determinant_a_excludes_zero": record["determinant_a_excludes_zero"],
        "determinant_schur_excludes_zero": record["determinant_schur_excludes_zero"],
        "determinant_a_balls": {
            name: serialize_ball(value)
            for name, value in record["determinant_a_balls"].items()
        },
        "determinant_schur_balls": {
            name: serialize_ball(value)
            for name, value in record["determinant_schur_balls"].items()
        },
        "strong": {
            "resolved_count": int(np.sum(record["strong"]["resolved"])),
            "zero_count": int(np.sum(record["strong"]["zero"])),
            "open_count": int(np.sum(record["strong"]["open"])),
            "singular_values": [
                serialize_float(value) for value in record["strong"]["singular_values"]
            ],
            "epsilon_step": serialize_float(record["strong"]["epsilon_step"]),
            "epsilon_ball": serialize_float(record["strong"]["epsilon_ball"]),
            "epsilon_svd": serialize_float(record["strong"]["epsilon_svd"]),
            "epsilon_global": serialize_float(record["strong"]["epsilon_global"]),
        },
        "schur": {
            "resolved_count": int(np.sum(record["schur"]["resolved"])),
            "zero_count": int(np.sum(record["schur"]["zero"])),
            "open_count": int(np.sum(record["schur"]["open"])),
            "singular_values": [
                serialize_float(value) for value in record["schur"]["singular_values"]
            ],
            "epsilon_step": serialize_float(record["schur"]["epsilon_step"]),
            "epsilon_ball": serialize_float(record["schur"]["epsilon_ball"]),
            "epsilon_svd": serialize_float(record["schur"]["epsilon_svd"]),
            "epsilon_global": serialize_float(record["schur"]["epsilon_global"]),
            "midpoint_matrix": [
                [serialize_complex(value) for value in row]
                for row in record["schur_midpoint"]
            ],
            "maximum_entry_radius": serialize_float(
                np.max(record["schur_radii"])
            ),
            "maximum_solve_entry_radius": serialize_float(
                np.max(record["solve_radii"])
            ),
        },
        "strong_resolved": record["strong_resolved"],
        "schur_resolved": record["schur_resolved"],
        "subspaces": {
            name: {
                "projector_distance": serialize_float(value["distance"]),
                "minimum_principal_overlap": serialize_float(
                    value["minimum_overlap"]
                ),
                "label": value["label"],
            }
            for name, value in record["subspaces"].items()
        },
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "precision": {"mpmath_decimal_digits": DPS, "flint_decimal_digits": BALL_DPS},
    "derivative_steps": {
        name: mp.nstr(step, 20) for name, step in DERIVATIVE_STEPS.items()
    },
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "outcome": combined_outcome,
    "parities": {
        parity: {
            **{key: value for key, value in record.items()
               if key not in {"basis_control", "branch_control", "kernel_control", "sectors"}},
            "basis_control": {
                key: value if isinstance(value, list) else serialize_mp(value)
                for key, value in record["basis_control"].items()
            },
            "branch_control": {
                "minimum_leading_minor": serialize_mp(
                    record["branch_control"]["minimum_leading_minor"]
                ),
                "minimum_argument": serialize_mp(
                    record["branch_control"]["minimum_argument"]
                ),
                "maximum_imaginary": serialize_mp(
                    record["branch_control"]["maximum_imaginary"]
                ),
                "base_negative_counts": dict(
                    record["branch_control"]["base_negative_counts"]
                ),
                "displaced_negative_counts": dict(
                    record["branch_control"]["displaced_negative_counts"]
                ),
            },
            "kernel_control": {
                "nonzero_kernel_entries": record["kernel_control"]["nonzero_kernel_entries"],
                "maximum_imaginary": serialize_mp(
                    record["kernel_control"]["maximum_imaginary"]
                ),
            },
            "sectors": [serialize_sector(item) for item in record["sectors"]],
        }
        for parity, record in records.items()
    },
    "classification": {
        "schur_rank": (
            "DERIVED COMPUTATIONAL"
            if combined_outcome == "FULL_LAPSE_SCHUR_REGULAR"
            else "OPEN_OR_DEGENERATE"
        ),
        "geometric_lapse_identification": "STRUCTURAL COMPARISON",
        "pseudo_constraint_interpretation": "OPEN",
        "gauge_identification": "OPEN",
        "graviton_interpretation": "OPEN",
        "limiting_speed": "OPEN",
    },
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Summary: {passed}/{tests} checks passed")
print(f"Outcome: {combined_outcome}")
for parity, record in records.items():
    print(
        f"  {parity}: Schur rank={record['resolved_schur_rank']}/120, "
        f"zero={record['schur_zero_count']}, open={record['schur_open_count']}, "
        f"subspaces={record['subspace_labels']}"
    )
print(f"Artifact: {OUTPUT}")

if passed != tests:
    raise SystemExit(1)

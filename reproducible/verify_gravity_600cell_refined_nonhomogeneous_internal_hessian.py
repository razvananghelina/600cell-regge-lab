#!/usr/bin/env python3
"""Complete local internal-Hessian census on the refined 600-cell slab.

Prior-art commit: d4dc6c7.
Protocol commit: fdf6f89.
No spatial-mode, continuum, speed, or particle target is loaded.
"""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_refined_nonhomogeneous_internal_hessian.json"
ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
CURVATURE = HERE / "gravity_600cell_refined_local_curvature_mass.json"
CURVATURE_SOURCE = HERE / "verify_gravity_600cell_refined_local_curvature_mass.py"
NULL_ARTIFACT = HERE / "gravity_600cell_refined_h4_null_coupling_adversarial.json"
NULL_RESULT = ROOT / "docs/gravity/gravity_600cell_refined_h4_null_coupling_result.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_internal_hessian_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_internal_hessian_protocol.md"

PRIOR_ART_COMMIT = "d4dc6c7"
PROTOCOL_COMMIT = "fdf6f89"
EXPECTED_HASHES = {
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py":
        "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "reproducible/gravity_600cell_refined_local_curvature_mass.json":
        "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "reproducible/verify_gravity_600cell_refined_local_curvature_mass.py":
        "c54f17708a2678b925cfce96fcfc7d6baaeeaf0577bedbf22b5d0435c069fae6",
    "reproducible/gravity_600cell_refined_h4_null_coupling_adversarial.json":
        "5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9",
    "docs/gravity/gravity_600cell_refined_h4_null_coupling_result.md":
        "660a3707f24f44d0393e6a1804e407fa45aa4782a98960438a296da50c35825a",
    "docs/gravity/gravity_600cell_refined_nonhomogeneous_internal_hessian_prior_art.md":
        "6e766435caa20404fcc9d30403cc27a969b827ab4d3fc232ce3ee58b1f90cd38",
}

PAIR4 = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR4)}
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_TRIANGLES = tuple(combinations(range(5), 3))
TRIANGLE_EDGES = tuple(combinations(range(3), 2))
TAU_TEXT = "0.0102"
DECIMAL_PRECISIONS = (100, 140)
DIFFERENCE_STEPS = ("1e-10", "5e-11", "2.5e-11")
EXPECTED_F = (2640, 17040, 28800, 14400)
EXPECTED_COUNTS = {
    "pentachora": 57600,
    "triangles": 149280,
    "boundary_edges": 34080,
    "internal_edges": 19680,
    "cross_edges": 17040,
    "vertical_edges": 2640,
    "all_edges": 53760,
}

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
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=50):
    return mp.nstr(value, digits, strip_zeros=False)


def float_text(value):
    return format(float(value), ".17e")


def load_action_definitions():
    tree = ast.parse(ACTION_SOURCE.read_text(), filename=str(ACTION_SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "mp": mp,
        "np": np,
        "json": json,
        "Path": Path,
        "sha256": sha256,
        "combinations": combinations,
        "permutations": permutations,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "HERE": HERE,
        "ROOT": ROOT,
        "PAIR4": PAIR4,
        "LOCAL_TRIANGLES": np.asarray(LOCAL_TRIANGLES, dtype=np.int8),
        "TAU_TEXT": TAU_TEXT,
        "VARIABLES": (
            tuple(("old",) + pair for pair in PAIR4)
            + tuple(("new",) + pair for pair in PAIR4)
            + tuple(("cross",) + pair for pair in PAIR4)
            + tuple(("rho", rank) for rank in range(4))
        ),
        "INTERNAL_VARIABLES": (
            tuple(("cross",) + pair for pair in PAIR4)
            + tuple(("rho", rank) for rank in range(4))
        ),
        "FD_STEP_TEXTS": ("1e-15", "5e-16"),
        "FD_GATE_TEXT": "1e-24",
        "EXPECTED_F": EXPECTED_F,
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(ACTION_SOURCE), "exec"), namespace)
    return namespace


def kind_from_states(left_state, right_state):
    left_rank, left_layer = left_state
    right_rank, right_layer = right_state
    if left_layer == right_layer:
        return ("spatial",) + tuple(sorted((left_rank, right_rank)))
    if left_rank == right_rank:
        return ("vertical", left_rank)
    return ("cross",) + tuple(sorted((left_rank, right_rank)))


def kind_value(kind, geometry):
    family = kind[0]
    if family == "vertical":
        return -(mp.mpf(TAU_TEXT) ** 2)
    pair = (kind[1], kind[2])
    spatial = geometry["s0"] ** 2 * geometry["unit_squares"][pair]
    if family == "cross":
        return spatial - mp.mpf(TAU_TEXT) ** 2
    return spatial


def squared_from_values(values):
    squared = [[mp.mpf(0) for _ in range(5)] for _ in range(5)]
    for value, (left, right) in zip(values, LOCAL_EDGES):
        squared[left][right] = squared[right][left] = value
    return squared


def gram_negative_count(values):
    squared = squared_from_values(values)
    gram = mp.matrix([
        [
            (squared[0][left] + squared[0][right] - squared[left][right]) / 2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    leading = []
    for size in range(1, 5):
        leading.append(mp.det(mp.matrix([
            [gram[row, column] for column in range(size)]
            for row in range(size)
        ])))
    signs = [1] + [1 if value > 0 else -1 if value < 0 else 0 for value in leading]
    if 0 in signs:
        return None, min(abs(value) for value in leading)
    return sum(left != right for left, right in zip(signs, signs[1:])), min(
        abs(value) for value in leading
    )


def angle_vector(actions, values):
    angles, identity, argument = actions["angle_record"](
        squared_from_values(values)
    )
    return mp.matrix([angles[hinge] for hinge in LOCAL_TRIANGLES]), identity, argument


def angle_stencil_at_precision(actions, kinds, geometry, dps):
    with mp.workdps(dps):
        values = [kind_value(kind, geometry) for kind in kinds]
        base, identity, argument = angle_vector(actions, values)
        negative, minor = gram_negative_count(values)
        central = []
        diagnostics = [(negative, minor, argument, identity)]
        for step_text in DIFFERENCE_STEPS:
            step = mp.mpf(step_text)
            derivative = mp.matrix(10, 10)
            for column in range(10):
                plus = list(values)
                minus = list(values)
                plus[column] *= mp.exp(step)
                minus[column] *= mp.exp(-step)
                a_plus, i_plus, z_plus = angle_vector(actions, plus)
                a_minus, i_minus, z_minus = angle_vector(actions, minus)
                n_plus, m_plus = gram_negative_count(plus)
                n_minus, m_minus = gram_negative_count(minus)
                diagnostics.extend((
                    (n_plus, m_plus, z_plus, i_plus),
                    (n_minus, m_minus, z_minus, i_minus),
                ))
                for row in range(10):
                    derivative[row, column] = (
                        a_plus[row] - a_minus[row]
                    ) / (2 * step)
            central.append(derivative)
        coarse = (4 * central[1] - central[0]) / 3
        fine = (4 * central[2] - central[1]) / 3
        return {
            "base": [mp.mpc(value) for value in base],
            "coarse": coarse,
            "fine": fine,
            "diagnostics": diagnostics,
        }


def build_angle_cache(actions, pattern_kinds, geometries):
    cache = []
    maximum_envelope = mp.mpf(0)
    maximum_identity = mp.mpf(0)
    minimum_argument = mp.inf
    minimum_minor = mp.inf
    all_lorentzian = True
    for pattern_index, kinds in enumerate(pattern_kinds):
        print(
            f"  differentiating local pentachoron pattern "
            f"{pattern_index + 1}/{len(pattern_kinds)}",
            flush=True,
        )
        low = angle_stencil_at_precision(
            actions, kinds, geometries[DECIMAL_PRECISIONS[0]],
            DECIMAL_PRECISIONS[0],
        )
        high = angle_stencil_at_precision(
            actions, kinds, geometries[DECIMAL_PRECISIONS[1]],
            DECIMAL_PRECISIONS[1],
        )
        selected = np.empty((10, 10), dtype=np.complex128)
        envelope = np.empty((10, 10), dtype=np.float64)
        for row in range(10):
            for column in range(10):
                value = high["fine"][row, column]
                error = 100 * (
                    abs(high["fine"][row, column] - high["coarse"][row, column])
                    + abs(high["fine"][row, column] - low["fine"][row, column])
                    + mp.mpf("1e-70")
                )
                selected[row, column] = complex(value)
                envelope[row, column] = float(error)
                maximum_envelope = max(maximum_envelope, error)
        for record in low["diagnostics"] + high["diagnostics"]:
            negative, minor, argument, identity = record
            all_lorentzian &= negative == 1
            minimum_minor = min(minimum_minor, minor)
            minimum_argument = min(minimum_argument, argument)
            maximum_identity = max(maximum_identity, identity)
        cache.append({
            "base": np.asarray([complex(value) for value in high["base"]]),
            "derivative": selected,
            "envelope": envelope,
        })
    return tuple(cache), {
        "all_displaced_lorentzian": bool(all_lorentzian),
        "maximum_angle_identity_residual": mp_text(maximum_identity, 30),
        "minimum_logarithm_argument": mp_text(minimum_argument, 30),
        "minimum_gram_leading_minor": mp_text(minimum_minor, 30),
        "maximum_local_derivative_envelope": mp_text(maximum_envelope, 30),
    }


def area_data(values):
    x = list(values)
    q = (
        2 * (x[0] * x[1] + x[0] * x[2] + x[1] * x[2])
        - x[0] ** 2 - x[1] ** 2 - x[2] ** 2
    ) / 16
    area = mp.sqrt(mp.mpc(q))
    q_plain = (
        (x[1] + x[2] - x[0]) / 8,
        (x[0] + x[2] - x[1]) / 8,
        (x[0] + x[1] - x[2]) / 8,
    )
    q_log = [q_plain[index] * x[index] for index in range(3)]
    q_hessian = [[mp.mpf(0) for _ in range(3)] for _ in range(3)]
    for row in range(3):
        for column in range(3):
            plain_second = mp.mpf(-1) / 8 if row == column else mp.mpf(1) / 8
            q_hessian[row][column] = plain_second * x[row] * x[column]
            if row == column:
                q_hessian[row][column] += q_plain[row] * x[row]
    gradient = [value / (2 * area) for value in q_log]
    hessian = [[
        q_hessian[row][column] / (2 * area)
        - q_log[row] * q_log[column] / (4 * area ** 3)
        for column in range(3)
    ] for row in range(3)]
    return area, gradient, hessian


def build_area_cache(triangle_kinds, geometry):
    result = []
    with mp.workdps(DECIMAL_PRECISIONS[1]):
        for kinds in triangle_kinds:
            _, gradient, hessian = area_data([
                kind_value(kind, geometry) for kind in kinds
            ])
            result.append({
                "gradient": np.asarray([complex(value) for value in gradient]),
                "hessian": np.asarray([
                    [complex(value) for value in row] for row in hessian
                ]),
            })
    return tuple(result)


def encoded_edges(rows, vertex_total):
    rows = np.asarray(rows, dtype=np.int64)
    return rows[:, 0] * vertex_total + rows[:, 1]


def schedule_geometry(actions, top, colours, order):
    slab = actions["staircase_slab"](top, colours, order)
    vertex_count = len(colours)
    slab_vertex_count = 2 * vertex_count

    simplex_edge_rows = np.sort(slab[:, np.asarray(LOCAL_EDGES)], axis=2)
    all_edges, edge_inverse = np.unique(
        simplex_edge_rows.reshape(-1, 2), axis=0, return_inverse=True
    )
    simplex_edges = edge_inverse.reshape(len(slab), 10).astype(np.int32)
    lower = all_edges[:, 1] < vertex_count
    upper = all_edges[:, 0] >= vertex_count
    boundary_mask = lower | upper
    internal_mask = ~boundary_mask
    internal_edges = all_edges[internal_mask]
    global_to_internal = np.full(len(all_edges), -1, dtype=np.int32)
    global_to_internal[internal_mask] = np.arange(
        len(internal_edges), dtype=np.int32
    )
    simplex_internal = global_to_internal[simplex_edges]

    triangle_rows = np.sort(
        slab[:, np.asarray(LOCAL_TRIANGLES)].reshape(-1, 3), axis=1
    )
    triangles, triangle_inverse = np.unique(
        triangle_rows, axis=0, return_inverse=True
    )
    simplex_triangles = triangle_inverse.reshape(len(slab), 10).astype(np.int32)
    triangle_edge_rows = np.sort(
        triangles[:, np.asarray(TRIANGLE_EDGES)], axis=2
    )
    edge_codes = encoded_edges(all_edges, slab_vertex_count)
    triangle_edge_codes = encoded_edges(
        triangle_edge_rows.reshape(-1, 2), slab_vertex_count
    )
    triangle_edges = np.searchsorted(edge_codes, triangle_edge_codes).reshape(-1, 3)
    if not np.array_equal(
        edge_codes[triangle_edges.ravel()], triangle_edge_codes
    ):
        raise RuntimeError("triangle edge lookup failed")
    triangle_internal = global_to_internal[triangle_edges]
    triangle_boundary = (
        (triangles[:, 2] < vertex_count)
        | (triangles[:, 0] >= vertex_count)
    )

    extended_ranks = np.concatenate((colours, colours)).astype(np.int8)
    extended_layers = np.concatenate((
        np.zeros(vertex_count, dtype=np.int8),
        np.ones(vertex_count, dtype=np.int8),
    ))
    state_codes = extended_ranks + 4 * extended_layers
    simplex_state_rows = state_codes[slab]
    unique_simplex_states, simplex_pattern = np.unique(
        simplex_state_rows, axis=0, return_inverse=True
    )
    simplex_kinds = []
    for states in unique_simplex_states:
        decoded = tuple((int(value % 4), int(value // 4)) for value in states)
        simplex_kinds.append(tuple(
            kind_from_states(decoded[left], decoded[right])
            for left, right in LOCAL_EDGES
        ))

    triangle_state_rows = state_codes[triangles]
    unique_triangle_states, triangle_pattern = np.unique(
        triangle_state_rows, axis=0, return_inverse=True
    )
    triangle_kinds = []
    for states in unique_triangle_states:
        decoded = tuple((int(value % 4), int(value // 4)) for value in states)
        triangle_kinds.append(tuple(
            kind_from_states(decoded[left], decoded[right])
            for left, right in TRIANGLE_EDGES
        ))

    internal_orbits = np.empty(len(internal_edges), dtype=np.int8)
    internal_is_vertical = np.zeros(len(internal_edges), dtype=bool)
    for index, (left, right) in enumerate(internal_edges):
        base_left = int(left % vertex_count)
        base_right = int(right % vertex_count)
        rank_left = int(colours[base_left])
        rank_right = int(colours[base_right])
        if base_left == base_right:
            internal_orbits[index] = 6 + rank_left
            internal_is_vertical[index] = True
        else:
            internal_orbits[index] = PAIR_INDEX[tuple(sorted((rank_left, rank_right)))]

    counts = {
        "pentachora": len(slab),
        "triangles": len(triangles),
        "boundary_edges": int(boundary_mask.sum()),
        "internal_edges": len(internal_edges),
        "cross_edges": int((~internal_is_vertical).sum()),
        "vertical_edges": int(internal_is_vertical.sum()),
        "all_edges": len(all_edges),
    }
    return {
        "order": tuple(map(int, order)),
        "slab": slab,
        "all_edges": all_edges,
        "internal_edges": internal_edges,
        "simplex_edges": simplex_edges,
        "simplex_internal": simplex_internal,
        "triangles": triangles,
        "triangle_edges": triangle_edges,
        "triangle_internal": triangle_internal,
        "triangle_boundary": triangle_boundary,
        "simplex_triangles": simplex_triangles,
        "simplex_pattern": simplex_pattern.astype(np.int16),
        "simplex_kinds": tuple(simplex_kinds),
        "triangle_pattern": triangle_pattern.astype(np.int16),
        "triangle_kinds": tuple(triangle_kinds),
        "internal_orbits": internal_orbits,
        "internal_is_vertical": internal_is_vertical,
        "counts": counts,
    }


def global_pattern_catalogue(schedule_geometries):
    simplex_kinds = sorted({
        kinds for geometry in schedule_geometries for kinds in geometry["simplex_kinds"]
    })
    triangle_kinds = sorted({
        kinds for geometry in schedule_geometries for kinds in geometry["triangle_kinds"]
    })
    simplex_index = {kinds: index for index, kinds in enumerate(simplex_kinds)}
    triangle_index = {kinds: index for index, kinds in enumerate(triangle_kinds)}
    for geometry in schedule_geometries:
        local_to_global = np.asarray([
            simplex_index[kinds] for kinds in geometry["simplex_kinds"]
        ], dtype=np.int16)
        geometry["simplex_pattern_global"] = local_to_global[
            geometry["simplex_pattern"]
        ]
        tri_local_to_global = np.asarray([
            triangle_index[kinds] for kinds in geometry["triangle_kinds"]
        ], dtype=np.int16)
        geometry["triangle_pattern_global"] = tri_local_to_global[
            geometry["triangle_pattern"]
        ]
    return tuple(simplex_kinds), tuple(triangle_kinds)


def add_sparse_component(matrix, count_matrix, rows, columns, data, size):
    if len(data) == 0:
        return matrix, count_matrix
    component = sp.coo_matrix(
        (np.asarray(data, dtype=np.complex128), (rows, columns)),
        shape=(size, size),
    ).tocsr()
    component.sum_duplicates()
    counts = sp.coo_matrix(
        (np.ones(len(data), dtype=np.int16), (rows, columns)),
        shape=(size, size),
    ).tocsr()
    counts.sum_duplicates()
    return matrix + component, count_matrix + counts


def assemble_internal(
    schedule, angle_cache, area_cache, per_vertex_masses, geometry
):
    size = schedule["counts"]["internal_edges"]
    simplex_pattern = schedule["simplex_pattern_global"]
    triangle_pattern = schedule["triangle_pattern_global"]
    simplex_triangles = schedule["simplex_triangles"]
    triangle_internal = schedule["triangle_internal"]
    simplex_internal = schedule["simplex_internal"]

    angle_base = np.asarray([record["base"] for record in angle_cache])
    angle_derivative = np.asarray([record["derivative"] for record in angle_cache])
    angle_envelope = np.asarray([record["envelope"] for record in angle_cache])
    area_gradient = np.asarray([record["gradient"] for record in area_cache])
    area_hessian = np.asarray([record["hessian"] for record in area_cache])

    curvature = np.where(
        schedule["triangle_boundary"], np.pi, 2 * np.pi
    ).astype(np.complex128)
    simplex_angles = angle_base[simplex_pattern]
    np.add.at(curvature, simplex_triangles.ravel(), simplex_angles.ravel())

    gradient = np.zeros(size, dtype=np.complex128)
    gradient_abs_rows = np.zeros(size, dtype=np.float64)
    matrix = sp.csr_matrix((size, size), dtype=np.complex128)
    count_matrix = sp.csr_matrix((size, size), dtype=np.int16)
    absolute_row_terms = np.zeros(size, dtype=np.float64)
    derivative_error_rows = np.zeros(size, dtype=np.float64)

    coefficient = -1j * curvature
    tri_grad = area_gradient[triangle_pattern]
    tri_hess = area_hessian[triangle_pattern]

    for row_position in range(3):
        rows = triangle_internal[:, row_position]
        valid = rows >= 0
        data = coefficient[valid] * tri_grad[valid, row_position]
        np.add.at(gradient, rows[valid], data)
        np.add.at(gradient_abs_rows, rows[valid], np.abs(data))

    area_rows = []
    area_columns = []
    area_values = []
    for row_position in range(3):
        rows = triangle_internal[:, row_position]
        for column_position in range(3):
            columns = triangle_internal[:, column_position]
            valid = (rows >= 0) & (columns >= 0)
            values = coefficient[valid] * tri_hess[
                valid, row_position, column_position
            ]
            area_rows.append(rows[valid])
            area_columns.append(columns[valid])
            area_values.append(values)
            np.add.at(absolute_row_terms, rows[valid], np.abs(values))
    matrix, count_matrix = add_sparse_component(
        matrix,
        count_matrix,
        np.concatenate(area_rows),
        np.concatenate(area_columns),
        np.concatenate(area_values),
        size,
    )

    for hinge in range(10):
        triangle_indices = simplex_triangles[:, hinge]
        hinge_pattern = triangle_pattern[triangle_indices]
        rows_all = triangle_internal[triangle_indices]
        derivative = angle_derivative[simplex_pattern, hinge]
        derivative_envelope = angle_envelope[simplex_pattern, hinge]
        chunk_rows = []
        chunk_columns = []
        chunk_values = []
        for row_position in range(3):
            rows = rows_all[:, row_position]
            area_values = tri_grad[triangle_indices, row_position]
            for column_position in range(10):
                columns = simplex_internal[:, column_position]
                valid = (rows >= 0) & (columns >= 0)
                values = -1j * area_values[valid] * derivative[
                    valid, column_position
                ]
                errors = np.abs(area_values[valid]) * derivative_envelope[
                    valid, column_position
                ]
                chunk_rows.append(rows[valid])
                chunk_columns.append(columns[valid])
                chunk_values.append(values)
                np.add.at(absolute_row_terms, rows[valid], np.abs(values))
                np.add.at(derivative_error_rows, rows[valid], errors)
        matrix, count_matrix = add_sparse_component(
            matrix,
            count_matrix,
            np.concatenate(chunk_rows),
            np.concatenate(chunk_columns),
            np.concatenate(chunk_values),
            size,
        )

    tau = float(mp.mpf(TAU_TEXT))
    vertical_indices = np.flatnonzero(schedule["internal_is_vertical"])
    rank_by_vertical = schedule["internal_orbits"][vertical_indices] - 6
    dust_gradient = np.asarray([
        -4 * np.pi * float(per_vertex_masses[rank]) * tau
        for rank in rank_by_vertical
    ])
    dust_hessian = dust_gradient / 2
    gradient[vertical_indices] += dust_gradient
    gradient_abs_rows[vertical_indices] += np.abs(dust_gradient)
    dust = sp.coo_matrix(
        (dust_hessian, (vertical_indices, vertical_indices)),
        shape=(size, size),
    ).tocsr()
    matrix = matrix + dust
    count_matrix = count_matrix + sp.coo_matrix(
        (np.ones(len(vertical_indices), dtype=np.int16),
         (vertical_indices, vertical_indices)),
        shape=(size, size),
    ).tocsr()
    np.add.at(absolute_row_terms, vertical_indices, np.abs(dust_hessian))

    matrix.sum_duplicates()
    matrix.sort_indices()
    count_matrix.sum_duplicates()
    maximum_term_count = int(count_matrix.data.max())
    unit_roundoff = np.finfo(np.float64).eps / 2
    operations = maximum_term_count + 32
    gamma = operations * unit_roundoff / (1 - operations * unit_roundoff)
    rounding_rows = (gamma + 4 * unit_roundoff) * absolute_row_terms
    error_rows = derivative_error_rows + rounding_rows
    operator_error = float(np.max(error_rows))
    gradient_error_rows = (gamma + 4 * unit_roundoff) * gradient_abs_rows
    gradient_error = float(np.max(gradient_error_rows))

    antisymmetric = matrix - matrix.T
    antisymmetry = float(np.max(np.asarray(
        np.abs(antisymmetric).sum(axis=1)
    ).ravel()))
    imaginary = float(np.max(np.asarray(
        np.abs(matrix.imag).sum(axis=1)
    ).ravel()))
    gradient_imaginary = float(np.max(np.abs(gradient.imag)))
    gradient_residual = float(np.max(np.abs(gradient.real)))

    real = matrix.real.tocsr()
    real = ((real + real.T) * 0.5).tocsr()
    real.sum_duplicates()
    real.sort_indices()
    return real, {
        "curvature_maximum_imaginary": float(np.max(np.abs(curvature.imag))),
        "gradient_maximum_absolute": gradient_residual,
        "gradient_maximum_imaginary": gradient_imaginary,
        "gradient_forward_error": gradient_error,
        "raw_hessian_imaginary_row_norm": imaginary,
        "raw_hessian_antisymmetric_row_norm": antisymmetry,
        "operator_error_row_bound": operator_error,
        "maximum_term_count": maximum_term_count,
        "nnz": int(real.nnz),
    }


def aggregate_gradient(actions, combinatorics, geometry, masses, coordinates):
    evaluation = actions["evaluate_schedule"](
        combinatorics, geometry, coordinates
    )
    result = [
        evaluation["gradient"][("cross",) + pair] for pair in PAIR4
    ]
    for rank in range(4):
        rho = coordinates["rho", rank]
        result.append(
            evaluation["gradient"]["rho", rank]
            - 4 * mp.pi * masses[rank] * mp.sqrt(rho)
        )
    return mp.matrix(result)


def aggregate_hessian_at_precision(
    actions, combinatorics, geometry, masses, dps
):
    with mp.workdps(dps):
        base = actions["base_coordinates"](geometry)
        keys = tuple(("cross",) + pair for pair in PAIR4) + tuple(
            ("rho", rank) for rank in range(4)
        )
        central = []
        for step_text in DIFFERENCE_STEPS:
            step = mp.mpf(step_text)
            derivative = mp.matrix(10, 10)
            for column, key in enumerate(keys):
                plus = dict(base)
                minus = dict(base)
                plus[key] *= mp.exp(step)
                minus[key] *= mp.exp(-step)
                g_plus = aggregate_gradient(
                    actions, combinatorics, geometry, masses, plus
                )
                g_minus = aggregate_gradient(
                    actions, combinatorics, geometry, masses, minus
                )
                for row in range(10):
                    derivative[row, column] = (
                        g_plus[row] - g_minus[row]
                    ) / (2 * step)
            central.append(derivative)
        return (4 * central[1] - central[0]) / 3, (
            4 * central[2] - central[1]
        ) / 3


def aggregate_hessian(actions, combinatorics, geometries, masses):
    low_coarse, low_fine = aggregate_hessian_at_precision(
        actions, combinatorics, geometries[100], masses, 100
    )
    high_coarse, high_fine = aggregate_hessian_at_precision(
        actions, combinatorics, geometries[140], masses, 140
    )
    selected = np.empty((10, 10), dtype=np.float64)
    envelope = np.empty((10, 10), dtype=np.float64)
    for row in range(10):
        for column in range(10):
            selected[row, column] = float(mp.re(
                high_fine[row, column] + high_fine[column, row]
            ) / 2)
            envelope[row, column] = float(100 * (
                abs(high_fine[row, column] - high_coarse[row, column])
                + abs(high_fine[row, column] - low_fine[row, column])
                + mp.mpf("1e-70")
            ))
    return selected, float(np.max(envelope))


def orbit_pullback(matrix, orbit_indices):
    size = len(orbit_indices)
    indicators = sp.coo_matrix(
        (
            np.ones(size),
            (np.arange(size, dtype=np.int32), orbit_indices),
        ),
        shape=(size, 10),
    ).tocsr()
    return np.asarray((indicators.T @ matrix @ indicators).todense())


def product_tangent(schedule, geometry):
    rho = float(mp.mpf(TAU_TEXT) ** 2)
    values = np.empty(len(schedule["internal_edges"]), dtype=np.float64)
    for index, orbit in enumerate(schedule["internal_orbits"]):
        if orbit >= 6:
            values[index] = 1.0
        else:
            pair = PAIR4[int(orbit)]
            cross = float(
                geometry["s0"] ** 2 * geometry["unit_squares"][pair]
                - mp.mpf(TAU_TEXT) ** 2
            )
            values[index] = -rho / cross
    return values / np.linalg.norm(values)


def sparse_row_norm(matrix):
    if matrix.nnz == 0:
        return 0.0
    return float(np.max(np.asarray(np.abs(matrix).sum(axis=1)).ravel()))


def csr_digest(matrix):
    matrix = matrix.copy().tocsr()
    matrix.sum_duplicates()
    matrix.sort_indices()
    payload = b"".join((
        np.asarray(matrix.indptr, dtype="<i8").tobytes(),
        np.asarray(matrix.indices, dtype="<i4").tobytes(),
        np.asarray(matrix.data, dtype="<f8").tobytes(),
    ))
    return sha256(payload).hexdigest()


def reversal_permutation(forward, reverse, vertex_count):
    reverse_codes = encoded_edges(reverse["internal_edges"], 2 * vertex_count)
    swapped = forward["internal_edges"].copy()
    swapped = np.where(
        swapped < vertex_count, swapped + vertex_count, swapped - vertex_count
    )
    swapped.sort(axis=1)
    codes = encoded_edges(swapped, 2 * vertex_count)
    positions = np.searchsorted(reverse_codes, codes)
    if not np.array_equal(reverse_codes[positions], codes):
        raise RuntimeError("time-reversal internal-edge mapping failed")
    return positions.astype(np.int32)


def fixed_probe(size, probe_index):
    indices = np.arange(size, dtype=np.int64)
    if probe_index == 0:
        return np.ones(size)
    modulus = (3, 5, 7, 11, 13, 17, 19)[probe_index - 1]
    return np.where(((indices + 1) * (probe_index + 2)) % modulus < modulus // 2,
                    1.0, -1.0)


def bordered_spectrum(matrix, tangent, operator_error):
    size = matrix.shape[0]
    column = sp.csr_matrix(tangent.reshape(-1, 1))
    bordered = sp.bmat(
        [[matrix, column], [column.T, sp.csr_matrix((1, 1))]],
        format="csc",
    )
    norm_inf = sparse_row_norm(bordered)
    result = {
        "factorized": False,
        "resolved_nonsingular": False,
        "matrix_infinity_norm": float_text(norm_inf),
        "operator_error": float_text(operator_error),
    }
    try:
        factor = spla.splu(
            bordered, permc_spec="COLAMD", diag_pivot_thresh=1.0
        )
    except RuntimeError as error:
        result["factorization_error"] = str(error)
        return result
    result["factorized"] = True
    inverse = spla.LinearOperator(
        bordered.shape,
        matvec=factor.solve,
        rmatvec=lambda value: factor.solve(value, trans="T"),
        dtype=np.float64,
    )
    eigen_runs = []
    maximum_ritz = 0.0
    for run_index, tolerance in enumerate((1e-10, 1e-12)):
        v0 = fixed_probe(size + 1, run_index + 1)
        try:
            eigenvalues, eigenvectors = spla.eigsh(
                bordered,
                k=8,
                sigma=0.0,
                which="LM",
                OPinv=inverse,
                v0=v0,
                tol=tolerance,
                maxiter=20000,
            )
            order = np.argsort(np.abs(eigenvalues))
            eigenvalues = eigenvalues[order]
            eigenvectors = eigenvectors[:, order]
            residuals = np.asarray([
                np.linalg.norm(
                    bordered @ eigenvectors[:, index]
                    - eigenvalues[index] * eigenvectors[:, index]
                )
                for index in range(8)
            ])
            maximum_ritz = max(maximum_ritz, float(residuals.max()))
            eigen_runs.append({
                "converged": True,
                "tolerance": float_text(tolerance),
                "eigenvalues_nearest_zero": [float_text(value) for value in eigenvalues],
                "ritz_residuals": [float_text(value) for value in residuals],
            })
        except Exception as error:  # ARPACK reports convergence in exception state.
            eigen_runs.append({
                "converged": False,
                "tolerance": float_text(tolerance),
                "error": repr(error),
            })

    solve_residuals = []
    for probe_index in range(8):
        rhs = fixed_probe(size + 1, probe_index)
        solution = factor.solve(rhs)
        residual = np.linalg.norm(bordered @ solution - rhs, ord=np.inf) / max(
            1.0, np.linalg.norm(rhs, ord=np.inf)
        )
        solve_residuals.append(float(residual))
    maximum_solve = max(solve_residuals)
    uncertainty = operator_error + maximum_ritz + maximum_solve * norm_inf

    paired = False
    minimum_absolute = 0.0
    if len(eigen_runs) == 2 and all(run["converged"] for run in eigen_runs):
        first = np.asarray([
            float(value) for value in eigen_runs[0]["eigenvalues_nearest_zero"]
        ])
        second = np.asarray([
            float(value) for value in eigen_runs[1]["eigenvalues_nearest_zero"]
        ])
        paired = bool(np.max(np.abs(first - second)) <= 100 * uncertainty)
        minimum_absolute = float(min(np.min(np.abs(first)), np.min(np.abs(second))))
    solved = maximum_solve <= 100 * uncertainty / max(1.0, norm_inf)
    separated = minimum_absolute > 100 * uncertainty
    result.update({
        "eigen_runs": eigen_runs,
        "solve_residuals": [float_text(value) for value in solve_residuals],
        "maximum_solve_residual": float_text(maximum_solve),
        "maximum_ritz_residual": float_text(maximum_ritz),
        "uncertainty": float_text(uncertainty),
        "paired_eigensolves": paired,
        "solves_within_gate": bool(solved),
        "minimum_absolute_eigenvalue": float_text(minimum_absolute),
        "spectrally_separated": bool(separated),
        "resolved_nonsingular": bool(paired and solved and separated),
        "lu_permutation_row_sha256": sha256(
            np.asarray(factor.perm_r, dtype="<i4").tobytes()
        ).hexdigest(),
        "lu_permutation_column_sha256": sha256(
            np.asarray(factor.perm_c, dtype="<i4").tobytes()
        ).hexdigest(),
    })
    return result


def corruption_control(schedule, angle_cache, area_cache, operator_error):
    # Deterministically use the first actual incidence with both an internal
    # hinge edge and an internal simplex edge.  Corrupt exactly that one local
    # derivative coefficient, never a complete symmetry orbit.
    for simplex in range(len(schedule["slab"])):
        pattern = int(schedule["simplex_pattern_global"][simplex])
        for hinge in range(10):
            triangle = int(schedule["simplex_triangles"][simplex, hinge])
            tri_pattern = int(schedule["triangle_pattern_global"][triangle])
            for row_position in range(3):
                row = int(schedule["triangle_internal"][triangle, row_position])
                if row < 0:
                    continue
                for column_position in range(10):
                    column = int(schedule["simplex_internal"][simplex, column_position])
                    if column < 0:
                        continue
                    delta = float(abs(
                        -1j
                        * area_cache[tri_pattern]["gradient"][row_position]
                        * 1e-4
                    ))
                    return {
                        "simplex": simplex,
                        "pattern": pattern,
                        "hinge": hinge,
                        "row_position": row_position,
                        "column_position": column_position,
                        "matrix_entry_change": float_text(delta),
                        "detection_gate": float_text(100 * operator_error),
                        "detected": bool(delta > 100 * operator_error),
                    }
    raise RuntimeError("no internal corruption incidence found")


print("=" * 78)
print("COMPLETE REFINED NONHOMOGENEOUS INTERNAL HESSIAN")
print("=" * 78)

actual_hashes = {name: digest(ROOT / name) for name in EXPECTED_HASHES}
provenance_ok = check(
    "all frozen geometry, action, matter and null-line inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIOR_ART_COMMIT == "d4dc6c7"
    and PROTOCOL_COMMIT == "fdf6f89",
    str(actual_hashes),
)

curvature = json.loads(CURVATURE.read_text())
null_artifact = json.loads(NULL_ARTIFACT.read_text())
upstream_ok = check(
    "the frozen inputs carry the required accepted but scoped outcomes",
    curvature["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and curvature["tests"] == {"passed": 15, "total": 15}
    and null_artifact["outcome"]
        == "ADVERSARIAL_REFINED_H4_NULL_COUPLING_CORROBORATED"
    and null_artifact["tests"] == {"passed": 11, "total": 11},
)

actions = load_action_definitions()
definitions_ok = check(
    "only frozen action definitions are loaded and its top level is not executed",
    {
        "tetrahedra_from_adjacency",
        "barycentric_chambers",
        "all_simplices",
        "staircase_slab",
        "schedule_combinatorics",
        "exact_geometry",
        "base_coordinates",
        "angle_record",
        "evaluate_schedule",
    } <= set(actions)
    and "OUTPUT" not in actions,
)

_, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
_, top, colours = actions["barycentric_chambers"](coarse_top)
spatial_cells = actions["all_simplices"](tuple(map(tuple, top)))
orders = tuple(permutations(range(4)))
print("Rebuilding all 24 labelled local slab geometries...", flush=True)
schedules = tuple(
    schedule_geometry(actions, top, colours, order) for order in orders
)
simplex_kinds, triangle_kinds = global_pattern_catalogue(schedules)
topology_ok = check(
    "the carrier and every complete local slab have the preregistered census",
    tuple(len(layer) for layer in spatial_cells) == EXPECTED_F
    and len(schedules) == 24
    and all(schedule["counts"] == EXPECTED_COUNTS for schedule in schedules),
    f"simplex patterns={len(simplex_kinds)}, triangle patterns={len(triangle_kinds)}",
)

geometries = {
    dps: actions["exact_geometry"](dps) for dps in DECIMAL_PRECISIONS
}
for geometry in geometries.values():
    geometry["mass"] = mp.mpf(0)
masses_total = tuple(
    mp.mpf(value)
    for value in curvature["selected_rank_matter"]["total_masses"]
)
masses_per_vertex = tuple(
    mp.mpf(value)
    for value in curvature["selected_rank_matter"]["per_vertex_masses"]
)

print("Building target-free high-precision local stencil catalogue...", flush=True)
angle_cache, stencil_diagnostics = build_angle_cache(
    actions, simplex_kinds, geometries
)
area_cache = build_area_cache(triangle_kinds, geometries[140])
stencil_ok = check(
    "all local derivative stencils retain the preregistered Lorentzian branch",
    stencil_diagnostics["all_displaced_lorentzian"]
    and mp.mpf(stencil_diagnostics["minimum_logarithm_argument"]) > 0
    and mp.mpf(stencil_diagnostics["minimum_gram_leading_minor"]) > 0,
    str(stencil_diagnostics),
)

aggregate_combinatorics = tuple(
    actions["schedule_combinatorics"](top, colours, order) for order in orders
)

schedule_records = [None] * 24
pair_records = []
visited = set()
all_stationary = True
all_real_symmetric = True
all_aggregate = True
all_null = True
all_reversal = True
all_corruption = True
all_resolved = True
maximum_gradient_fraction = 0.0
maximum_reality_fraction = 0.0
maximum_aggregate_fraction = 0.0
maximum_null_fraction = 0.0
maximum_reversal_fraction = 0.0

for forward_index, order in enumerate(orders):
    if forward_index in visited:
        continue
    reverse_order = tuple(reversed(order))
    reverse_index = orders.index(reverse_order)
    visited.update((forward_index, reverse_index))
    pair_data = []
    for schedule_index in (forward_index, reverse_index):
        schedule = schedules[schedule_index]
        print(
            f"Assembling schedule {schedule_index + 1}/24 order="
            f"{schedule['order']}...",
            flush=True,
        )
        matrix, diagnostics = assemble_internal(
            schedule,
            angle_cache,
            area_cache,
            masses_per_vertex,
            geometries[140],
        )
        aggregate, aggregate_error = aggregate_hessian(
            actions,
            aggregate_combinatorics[schedule_index],
            geometries,
            masses_total,
        )
        pullback = orbit_pullback(matrix, schedule["internal_orbits"])
        pullback_difference = float(np.max(np.abs(pullback - aggregate)))
        pullback_gate = 100 * (
            diagnostics["operator_error_row_bound"] + aggregate_error
        )
        tangent = product_tangent(schedule, geometries[140])
        null_residual = float(np.linalg.norm(matrix @ tangent, ord=np.inf))
        multiplication_error = (
            np.finfo(np.float64).eps
            * sparse_row_norm(matrix)
            * np.linalg.norm(tangent, ord=1)
        )
        null_gate = 100 * (
            diagnostics["operator_error_row_bound"] + multiplication_error
        )
        gradient_gate = 100 * max(
            diagnostics["gradient_forward_error"], np.finfo(float).tiny
        )
        reality_gate = 100 * max(
            diagnostics["operator_error_row_bound"], np.finfo(float).tiny
        )
        stationary = (
            diagnostics["gradient_maximum_absolute"] <= gradient_gate
            and diagnostics["gradient_maximum_imaginary"] <= gradient_gate
        )
        real_symmetric = (
            diagnostics["raw_hessian_imaginary_row_norm"] <= reality_gate
            and diagnostics["raw_hessian_antisymmetric_row_norm"] <= reality_gate
        )
        aggregate_ok = pullback_difference <= pullback_gate
        null_ok = null_residual <= null_gate
        maximum_gradient_fraction = max(
            maximum_gradient_fraction,
            diagnostics["gradient_maximum_absolute"] / gradient_gate,
            diagnostics["gradient_maximum_imaginary"] / gradient_gate,
        )
        maximum_reality_fraction = max(
            maximum_reality_fraction,
            diagnostics["raw_hessian_imaginary_row_norm"] / reality_gate,
            diagnostics["raw_hessian_antisymmetric_row_norm"] / reality_gate,
        )
        maximum_aggregate_fraction = max(
            maximum_aggregate_fraction, pullback_difference / pullback_gate
        )
        maximum_null_fraction = max(
            maximum_null_fraction, null_residual / null_gate
        )
        all_stationary &= stationary
        all_real_symmetric &= real_symmetric
        all_aggregate &= aggregate_ok
        all_null &= null_ok
        record = {
            "order": list(schedule["order"]),
            "counts": schedule["counts"],
            "csr_sha256": csr_digest(matrix),
            "diagnostics": {
                key: float_text(value) if isinstance(value, float) else value
                for key, value in diagnostics.items()
            },
            "gradient_gate": float_text(gradient_gate),
            "reality_gate": float_text(reality_gate),
            "aggregate_pullback_maximum_difference": float_text(
                pullback_difference
            ),
            "aggregate_pullback_gate": float_text(pullback_gate),
            "product_tangent_residual": float_text(null_residual),
            "product_tangent_gate": float_text(null_gate),
            "stationary": bool(stationary),
            "real_symmetric": bool(real_symmetric),
            "aggregate_pullback_matches": bool(aggregate_ok),
            "product_tangent_is_null": bool(null_ok),
        }
        schedule_records[schedule_index] = record
        pair_data.append((schedule, matrix, tangent, diagnostics))

    forward, reverse = pair_data
    permutation = reversal_permutation(
        forward[0], reverse[0], len(colours)
    )
    reverse_pulled = reverse[1][permutation, :][:, permutation]
    reversal_difference = sparse_row_norm(forward[1] - reverse_pulled)
    reversal_gate = 100 * (
        forward[3]["operator_error_row_bound"]
        + reverse[3]["operator_error_row_bound"]
    )
    reversal_ok = reversal_difference <= reversal_gate
    all_reversal &= reversal_ok
    maximum_reversal_fraction = max(
        maximum_reversal_fraction, reversal_difference / reversal_gate
    )

    corruption = corruption_control(
        forward[0], angle_cache, area_cache,
        forward[3]["operator_error_row_bound"],
    )
    all_corruption &= corruption["detected"]
    spectrum = bordered_spectrum(
        forward[1], forward[2], forward[3]["operator_error_row_bound"]
    )
    all_resolved &= spectrum["resolved_nonsingular"]
    if spectrum.get("eigen_runs"):
        forward_values = spectrum["eigen_runs"][0].get(
            "eigenvalues_nearest_zero", []
        )
    else:
        forward_values = []
    reverse_spectrum_match = bool(reversal_ok)
    pair_records.append({
        "forward_schedule": forward_index,
        "reverse_schedule": reverse_index,
        "forward_order": list(orders[forward_index]),
        "reverse_order": list(orders[reverse_index]),
        "reversal_maximum_row_difference": float_text(reversal_difference),
        "reversal_gate": float_text(reversal_gate),
        "reversal_covariant": bool(reversal_ok),
        "reverse_spectrum_inherited_by_exact_congruence": reverse_spectrum_match,
        "corruption_control": corruption,
        "bordered_spectrum": spectrum,
        "nearest_zero_values_for_pair": forward_values,
    })

construction_ok = all(
    (provenance_ok, upstream_ok, definitions_ok, topology_ok, stencil_ok,
     all_stationary, all_real_symmetric, all_aggregate, all_null,
     all_reversal, all_corruption)
)
if not construction_ok:
    outcome = "LOCAL_EXTENSION_INVALID"
elif all_resolved:
    outcome = "COMPLETE_INTERNAL_KERNEL_IS_PRODUCT_DURATION_LINE"
else:
    outcome = "COMPLETE_INTERNAL_KERNEL_NUMERICALLY_OPEN"

check(
    "all 24 individual local internal gradients are stationary inside error",
    all_stationary,
    f"maximum used fraction={maximum_gradient_fraction:.6e}",
)
check(
    "all 24 raw local Hessians are real and reciprocal inside error",
    all_real_symmetric,
    f"maximum used fraction={maximum_reality_fraction:.6e}",
)
check(
    "all 24 complete local pullbacks reproduce the aggregate ten-orbit blocks",
    all_aggregate,
    f"maximum used fraction={maximum_aggregate_fraction:.6e}",
)
check(
    "the analytic product-duration tangent is null in all 24 complete blocks",
    all_null,
    f"maximum used fraction={maximum_null_fraction:.6e}",
)
check(
    "all 12 explicit schedule/reverse pairs are congruent under layer reversal",
    all_reversal,
    f"maximum used fraction={maximum_reversal_fraction:.6e}",
)
check(
    "the preregistered one-incidence stencil corruption is detected in every pair",
    all_corruption,
)
classification_ok = check(
    "the complete bordered census follows the frozen verdict hierarchy",
    outcome in {
        "LOCAL_EXTENSION_INVALID",
        "COMPLETE_INTERNAL_KERNEL_NUMERICALLY_OPEN",
        "COMPLETE_INTERNAL_KERNEL_IS_PRODUCT_DURATION_LINE",
    }
    and len(pair_records) == 12
    and len(schedule_records) == 24
    and all(record is not None for record in schedule_records),
    outcome,
)

artifact = {
    "title": "Complete refined nonhomogeneous internal-Hessian census",
    "date": "2026-08-21",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "carrier": "K0=P(sd K_600)",
        "f_vector": list(EXPECTED_F),
        "tau0": TAU_TEXT,
        "coordinates": "logarithmic absolute signed squared local edge lengths",
        "boundary_dimension": 34080,
        "internal_dimension": 19680,
        "cross_dimension": 17040,
        "vertical_dimension": 2640,
        "schedule_count": 24,
        "time_reversal_pair_count": 12,
        "decimal_precisions": list(DECIMAL_PRECISIONS),
        "difference_steps": list(DIFFERENCE_STEPS),
        "simplex_pattern_count": len(simplex_kinds),
        "triangle_pattern_count": len(triangle_kinds),
        "matter_perturbations": "none; conserved curvature-selected masses frozen",
    },
    "stencil_diagnostics": stencil_diagnostics,
    "controls": {
        "all_stationary": bool(all_stationary),
        "all_real_symmetric": bool(all_real_symmetric),
        "all_aggregate_pullbacks_match": bool(all_aggregate),
        "all_product_tangents_null": bool(all_null),
        "all_time_reversal_pairs_covariant": bool(all_reversal),
        "all_corruptions_detected": bool(all_corruption),
        "maximum_gradient_gate_fraction": float_text(maximum_gradient_fraction),
        "maximum_reality_gate_fraction": float_text(maximum_reality_fraction),
        "maximum_aggregate_gate_fraction": float_text(maximum_aggregate_fraction),
        "maximum_null_gate_fraction": float_text(maximum_null_fraction),
        "maximum_reversal_gate_fraction": float_text(maximum_reversal_fraction),
    },
    "census": {
        "schedules": schedule_records,
        "time_reversal_pairs": pair_records,
        "resolved_pair_count": sum(
            record["bordered_spectrum"]["resolved_nonsingular"]
            for record in pair_records
        ),
        "open_pair_count": sum(
            not record["bordered_spectrum"]["resolved_nonsingular"]
            for record in pair_records
        ),
    },
    "scope": {
        "all_internal_edge_directions_included": True,
        "boundary_response_computed": False,
        "spatial_laplacian_or_continuum_target_loaded": False,
        "wave_speed_or_physical_constant_extracted": False,
        "schedule_selected_or_averaged": False,
        "full_suite_run": False,
    },
    "status_labels": {
        "kernel_census": (
            "DERIVED COMPUTATIONAL PRIMARY" if all_resolved
            else "OPEN" if construction_ok else "INVALID"
        ),
        "graviton_or_propagation": "OPEN / NOT TESTED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"OUTCOME: {outcome}")
print(f"Resolved bordered pairs: {artifact['census']['resolved_pair_count']}/12")
print(f"Artifact: {OUTPUT}")
print(f"Tests: {passed}/{tests}")
if passed != tests or not classification_ok:
    raise SystemExit(1)

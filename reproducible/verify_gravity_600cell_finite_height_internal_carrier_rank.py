#!/usr/bin/env python3
"""Target-free internal-equation rank on the finite-height 240 carrier.

Prior-art/framing commit: 3da6ec6.
Protocol commit: 941e36f.
Registry commit: 4ac19f2.

The complete ambient Hessian is deliberately never materialized.
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

import mpmath as mp
import networkx as nx
import numpy as np
from scipy.linalg import qr


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_finite_height_internal_carrier_rank.json"
MATRIX_OUTPUT = (
    HERE / "gravity_600cell_finite_height_internal_carrier_rank_matrices.npz"
)
PRIOR_ART = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_internal_carrier_rank_prior_art.md"
)
PRIMARY_JSON = HERE / "gravity_600cell_finite_height_carrier_quadratic.json"
ADVERSARIAL_JSON = (
    HERE / "gravity_600cell_finite_height_carrier_quadratic_adversarial.json"
)
ADVERSARIAL_MATRICES = (
    HERE / "gravity_600cell_finite_height_carrier_quadratic_adversarial_matrices.npy"
)
ADVERSARIAL_SOURCE = (
    HERE / "verify_gravity_600cell_finite_height_carrier_quadratic_adversarial.py"
)
HIGH_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RANK_SOURCE = (
    HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
)
SYMBOLIC_INPUT = (
    HERE / "gravity_600cell_full_scale_strut_symbolic_gap_resolution.json"
)
RUN_ALL = HERE / "run_all.py"

INPUTS = {
    "prior_art": PRIOR_ART,
    "primary_json": PRIMARY_JSON,
    "adversarial_json": ADVERSARIAL_JSON,
    "adversarial_matrices": ADVERSARIAL_MATRICES,
    "adversarial_source": ADVERSARIAL_SOURCE,
    "high_source": HIGH_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "symbolic_input": SYMBOLIC_INPUT,
}
EXPECTED_HASHES = {
    "prior_art": (
        "f3d04db084a63944e6963687747dcbe510d910d12f8b65b612ae50f6a1d89696"
    ),
    "primary_json": (
        "0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9"
    ),
    "adversarial_json": (
        "54915cf364c36af6bbc8e1dbd36433079269d293453478bfdf589e547d462ad6"
    ),
    "adversarial_matrices": (
        "8a3ea0c3b8ee720d8ffdf07e7486aefdd0247ca1cfdbeb99f443091376f31729"
    ),
    "adversarial_source": (
        "8d37012f556ce5be0bb863ad12d4572d197c90a5b96974912e81a98c1956a8f8"
    ),
    "high_source": (
        "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571"
    ),
    "geometry_source": (
        "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf"
    ),
    "symbolic_input": (
        "ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179"
    ),
}

PRIOR_ART_COMMIT = "3da6ec6"
PROTOCOL_COMMIT = "941e36f"
REGISTRY_COMMIT = "4ac19f2"
VERIFIER_NAME = Path(__file__).name

DPS = 180
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-25"),
    "operational_shadow": mp.mpf("5e-26"),
    "validation_primary": mp.mpf("2.5e-26"),
    "validation_shadow": mp.mpf("1.25e-26"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
ARITHMETIC_FLOOR = mp.mpf("1e-150")
RANK_FLOOR = mp.mpf("1e-135")
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {
    hinge: index for index, hinge in enumerate(LOCAL_HINGES)
}
I = mp.mpc(0, 1)
VERTICES = 120
OLD = 720
INTERNAL = 840
NEW = 720
FULL = OLD + INTERNAL + NEW
ACTIVE = INTERNAL + NEW
DATA = 2 * VERTICES

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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=80):
    return mp.nstr(value, digits)


def load_named_functions(path, names, namespace):
    tree = ast.parse(path.read_text(), filename=str(path))
    selected = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    found = {node.name for node in selected}
    missing = set(names) - found
    if missing:
        raise RuntimeError(f"missing audited functions: {sorted(missing)}")
    module = ast.Module(body=selected, type_ignores=[])
    ast.fix_missing_locations(module)
    exec(compile(module, str(path), "exec"), namespace)


def registry_inventory(path):
    tree = ast.parse(path.read_text(), filename=str(path))
    scripts = None
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == "scripts"
            for target in node.targets
        ):
            scripts = ast.literal_eval(node.value)
            break
    if scripts is None:
        raise RuntimeError("run_all.py has no literal scripts registry")
    counts = Counter(scripts)
    duplicates = sorted(name for name, count in counts.items() if count != 1)
    return scripts, duplicates


def mp_frobenius(matrix):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in matrix))


def matrix_difference(left, right):
    return mp.matrix([
        [
            left[row, column] - right[row, column]
            for column in range(left.cols)
        ]
        for row in range(left.rows)
    ])


def matrix_scale_columns(matrix, scale_count, scale_value, strut_value):
    result = mp.matrix(matrix.rows, matrix.cols)
    for row in range(matrix.rows):
        for column in range(matrix.cols):
            factor = scale_value if column < scale_count else strut_value
            result[row, column] = matrix[row, column] * factor
    return result


def richardson(coarse, fine):
    return mp.matrix([
        [
            (4 * fine[row, column] - coarse[row, column]) / 3
            for column in range(coarse.cols)
        ]
        for row in range(coarse.rows)
    ])


def mp_to_numpy(matrix):
    return np.asarray([
        [
            complex(
                float(mp.re(matrix[row, column])),
                float(mp.im(matrix[row, column])),
            )
            for column in range(matrix.cols)
        ]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


def finite_height_formula_control(primary):
    q = mp.mpf(primary["background"]["q"])
    h = mp.mpf(primary["background"]["h"])
    lam = mp.mpf(primary["background"]["lambda"])
    rho = mp.mpf(primary["background"]["rho"])
    mass = mp.mpf(primary["background"]["mass"])
    v = mp.mpf(3) / 2

    def epsilon(value):
        square = value**2
        return 2 * mp.pi - 5 * mp.acos(
            (square + 2) / (2 * (square + 3))
        )

    def mu(value):
        return 180 * epsilon(value) / (mp.pi * mp.sqrt(value**2 + 4))

    def momentum(value):
        square = value**2
        return (
            180 * value * epsilon(value) / mp.sqrt(square + 4)
            - 600 * mp.sqrt(3)
            * mp.asinh(value / mp.sqrt(8 * (square + 3)))
        )

    elimination = (
        4 * mp.pi * (mu(q) - mu(v))
        + q * (momentum(q) - momentum(v))
    )
    reconstructed_h = (
        momentum(q) - momentum(v)
    ) / (2 * mp.pi * mu(q))
    return {
        "q": q,
        "h": h,
        "lambda": lam,
        "rho": rho,
        "mass": mass,
        "elimination": elimination,
        "h_error": abs(h - reconstructed_h),
        "lambda_error": abs(lam - (1 + h * q)),
        "rho_error": abs(rho - h**2),
        "mass_error": abs(mass - mu(v)),
    }


def build_sparse_carrier(model, index_data, lam, rho):
    q_diag = lam - rho
    coefficient_a = -16 * rho / (lam - 1) ** 2
    coefficient_b = 8 + 16 * rho / (lam - 1) ** 2
    a_value = coefficient_a / (8 * q_diag)
    b_value = coefficient_b / (8 * q_diag)
    kappa = rho / ((lam - 1) * q_diag)

    rows = [dict() for _ in range(ACTIVE)]
    first_diagonal = None
    for edge, global_index in index_data["edge_to_index"].items():
        if global_index < OLD:
            continue
        edge = tuple(int(value) for value in edge)
        row = global_index - OLD
        kind = index_data["edge_kind"][global_index]
        if kind == "pole":
            lower, upper = edge
            if upper != lower + VERTICES:
                raise RuntimeError("invalid pole endpoint convention")
            rows[row][VERTICES + lower] = mp.mpf(1)
        elif kind == "internal":
            lower, upper = edge
            target = upper - VERTICES
            rows[row] = {
                lower: a_value,
                target: b_value,
                VERTICES + lower: kappa,
                VERTICES + target: -lam * kappa,
            }
            orbit_type = global_index // 24
            candidate = (min(index_data["orbit_edges"][orbit_type]), orbit_type)
            if first_diagonal is None or candidate[0] < first_diagonal[0]:
                first_diagonal = candidate
        elif kind == "new":
            left, right = edge[0] - VERTICES, edge[1] - VERTICES
            rows[row] = {left: 1 / lam, right: 1 / lam}
        else:
            raise RuntimeError(f"unexpected active edge kind {kind}")

    pole_types = [
        orbit_type for orbit_type in range(30, 65)
        if index_data["edge_kind"][24 * orbit_type] == "pole"
    ]
    vertex_orbit_group = {}
    for orbit_position, orbit_type in enumerate(pole_types):
        for group, edge in enumerate(index_data["orbit_edges"][orbit_type]):
            lower, upper = tuple(map(int, edge))
            if upper != lower + VERTICES or lower in vertex_orbit_group:
                raise RuntimeError("pole orbit does not define vertex order")
            vertex_orbit_group[lower] = (orbit_position, group)
    data_orbit_group = {}
    for vertex in range(VERTICES):
        orbit, group = vertex_orbit_group[vertex]
        data_orbit_group[vertex] = (orbit, group)
        data_orbit_group[VERTICES + vertex] = (5 + orbit, group)

    support = Counter(len(row) for row in rows)
    scale_support = Counter()
    strut_support = Counter()
    for row in rows:
        for column in row:
            if column < VERTICES:
                scale_support[column] += 1
            else:
                strut_support[column - VERTICES] += 1

    graph = nx.Graph()
    graph.add_nodes_from(range(VERTICES))
    for edge in model["new_edges"]:
        graph.add_edge(edge[0] - VERTICES, edge[1] - VERTICES)
    graph_ok = bool(
        graph.number_of_edges() == NEW
        and nx.is_connected(graph)
        and not nx.is_bipartite(graph)
    )
    pole_identity = bool(
        len(pole_types) == 5
        and len(vertex_orbit_group) == VERTICES
        and all(
            rows[index_data["edge_to_index"][(vertex, vertex + VERTICES)] - OLD]
            == {VERTICES + vertex: mp.mpf(1)}
            for vertex in range(VERTICES)
        )
    )
    algebra_error = max(
        abs(a_value + b_value - 1 / q_diag),
        abs(kappa - lam * kappa + rho / q_diag),
    )
    controls_ok = bool(
        support == Counter({1: 120, 2: 720, 4: 720})
        and all(scale_support[index] == 24 for index in range(VERTICES))
        and all(strut_support[index] == 13 for index in range(VERTICES))
        and pole_identity
        and graph_ok
        and algebra_error < mp.mpf("1e-160")
        and len(data_orbit_group) == DATA
    )
    return {
        "rows": rows,
        "pole_types": tuple(pole_types),
        "data_orbit_group": data_orbit_group,
        "first_diagonal_orbit": first_diagonal[1],
        "controls_ok": controls_ok,
        "controls": {
            "support": dict(sorted(support.items())),
            "scale_support": [min(scale_support.values()), max(scale_support.values())],
            "strut_support": [min(strut_support.values()), max(strut_support.values())],
            "pole_types": pole_types,
            "pole_identity": pole_identity,
            "connected_nonbipartite_upper_graph": graph_ok,
            "exact_rank_240": pole_identity and graph_ok,
            "algebra_error": mp_text(algebra_error),
            "coefficients": {
                "a": mp_text(a_value),
                "b": mp_text(b_value),
                "kappa": mp_text(kappa),
            },
        },
    }


def combine_kernels(coarse, fine):
    keys = set(coarse) | set(fine)
    return {
        key: (4 * fine.get(key, 0) - coarse.get(key, 0)) / 3
        for key in keys
    }


def compose_response(kernel, index_data, carrier_rows, reverse=False):
    table = index_data["table"]
    result = [defaultdict(lambda: mp.mpc(0)) for _ in range(ACTIVE)]
    maximum_kernel_imaginary = mp.mpf(0)
    expanded_entries = 0
    for (row_type, column_type, relative), value in kernel.items():
        if row_type < 30 or column_type < 30:
            continue
        maximum_kernel_imaginary = max(
            maximum_kernel_imaginary, abs(mp.im(value))
        )
        for row_group in range(24):
            column_group = int(
                table[relative, row_group]
                if reverse else table[row_group, relative]
            )
            row = 24 * (row_type - 30) + row_group
            column = 24 * (column_type - 30) + column_group
            for data_column, carrier_value in carrier_rows[column].items():
                result[row][data_column] += value * carrier_value
                expanded_entries += 1
    cleaned = [
        {column: value for column, value in row.items() if value != 0}
        for row in result
    ]
    return cleaned, {
        "expanded_products": expanded_entries,
        "nonzero_response_entries": sum(len(row) for row in cleaned),
        "maximum_kernel_imaginary": maximum_kernel_imaginary,
    }


def sparse_linear_combination(left_factor, left, right_factor, right):
    result = []
    for left_row, right_row in zip(left, right):
        keys = set(left_row) | set(right_row)
        row = {
            key: left_factor * left_row.get(key, 0)
            + right_factor * right_row.get(key, 0)
            for key in keys
        }
        result.append({key: value for key, value in row.items() if value != 0})
    return result


def sparse_richardson(coarse, fine):
    return sparse_linear_combination(mp.mpf(-1) / 3, coarse, mp.mpf(4) / 3, fine)


def sparse_frobenius(rows, selected_rows=None, imaginary=False):
    indices = range(len(rows)) if selected_rows is None else selected_rows
    if imaginary:
        return mp.sqrt(mp.fsum(
            abs(mp.im(value)) ** 2
            for row in indices for value in rows[row].values()
        ))
    return mp.sqrt(mp.fsum(
        abs(value) ** 2
        for row in indices for value in rows[row].values()
    ))


def sparse_difference_frobenius(left, right, selected_rows=None):
    indices = range(len(left)) if selected_rows is None else selected_rows
    return mp.sqrt(mp.fsum(
        abs(left[row].get(key, 0) - right[row].get(key, 0)) ** 2
        for row in indices
        for key in set(left[row]) | set(right[row])
    ))


def sparse_to_numpy(rows, selected_rows):
    result = np.zeros((len(selected_rows), DATA), dtype=np.complex128)
    for target, source in enumerate(selected_rows):
        for column, value in rows[source].items():
            result[target, column] = complex(
                float(mp.re(value)), float(mp.im(value))
            )
    return result


def project_sparse_response(
    rows, orbit_positions, data_orbit_group, sector
):
    dimension = int(sector["dimension"])
    basis = sector["basis"]
    result = mp.matrix(len(orbit_positions) * dimension, 10 * dimension)
    orbit_target = {
        orbit_position: target
        for target, orbit_position in enumerate(orbit_positions)
    }
    for global_row in range(INTERNAL):
        orbit_position, row_group = divmod(global_row, 24)
        if orbit_position not in orbit_target:
            continue
        target_orbit = orbit_target[orbit_position]
        for data_column, value in rows[global_row].items():
            data_orbit, data_group = data_orbit_group[data_column]
            for left_component in range(dimension):
                left = mp.conj(basis[row_group, left_component])
                if left == 0:
                    continue
                for right_component in range(dimension):
                    right = basis[data_group, right_component]
                    if right != 0:
                        result[
                            target_orbit * dimension + left_component,
                            data_orbit * dimension + right_component,
                        ] += left * value * right
    return result


def quadratic_from_response(carrier_rows, response_rows):
    result = mp.matrix(DATA, DATA)
    for row, carrier in enumerate(carrier_rows):
        response = response_rows[row]
        for left_column, left_value in carrier.items():
            for right_column, right_value in response.items():
                result[left_column, right_column] += (
                    mp.conj(left_value) * right_value
                )
    return mp.matrix([
        [
            (result[row, column] + mp.conj(result[column, row])) / 2
            for column in range(DATA)
        ]
        for row in range(DATA)
    ])


def column_block_frobenius(matrix, start, stop):
    return mp.sqrt(mp.fsum(
        abs(matrix[row, column]) ** 2
        for row in range(matrix.rows)
        for column in range(start, stop)
    ))


def singular_values(matrix):
    values = mp.svd(matrix, compute_uv=False)
    return [mp.mpf(value) for value in values]


def classify_values(values, normalization, error_total):
    labels = []
    normalized = []
    for value in values:
        ratio = abs(value) / normalization
        normalized.append(ratio)
        if ratio <= 10 * error_total:
            labels.append("ZERO")
        elif ratio > 100 * error_total:
            labels.append("NONZERO")
        else:
            labels.append("OPEN")
    return labels, normalized


def triplet_nullity(matrices):
    normalization = max(mp.mpf(1), mp_frobenius(matrices["M12"]))
    step_error = max(
        mp_frobenius(matrix_difference(matrices["M01"], matrices["M12"])),
        mp_frobenius(matrix_difference(matrices["M12"], matrices["M23"])),
    ) / normalization
    total_error = step_error + RANK_FLOOR
    level_records = {}
    nullities = []
    resolved = True
    for level in ("M01", "M12", "M23"):
        values = singular_values(matrices[level])
        labels, normalized = classify_values(
            values, normalization, total_error
        )
        resolved &= "OPEN" not in labels
        nullity = labels.count("ZERO")
        nullities.append(nullity)
        level_records[level] = {
            "singular_values": [mp_text(value) for value in values],
            "normalized_singular_values": [
                mp_text(value) for value in normalized
            ],
            "labels": labels,
            "nullity": nullity,
        }
    resolved &= len(set(nullities)) == 1
    return {
        "normalization": normalization,
        "step_error": step_error,
        "total_error": total_error,
        "resolved": bool(resolved),
        "nullity": nullities[0] if resolved else None,
        "rank": matrices["M12"].cols - nullities[0] if resolved else None,
        "levels": level_records,
    }


def deterministic_minor(matrices, rank):
    if rank == 0:
        return {
            "rank": 0,
            "rows": [],
            "columns": [],
            "determinants": {level: "1.0" for level in matrices},
            "error": "0.0",
            "gate": "0.0",
            "passed": True,
        }
    binary = mp_to_numpy(matrices["M12"])
    _, _, column_pivots = qr(binary, pivoting=True, mode="economic")
    columns = [int(value) for value in column_pivots[:rank]]
    _, _, row_pivots = qr(
        binary[:, columns].T, pivoting=True, mode="economic"
    )
    rows = [int(value) for value in row_pivots[:rank]]
    determinants = {}
    for level, matrix in matrices.items():
        minor = mp.matrix([
            [matrix[row, column] for column in columns]
            for row in rows
        ])
        determinants[level] = mp.det(minor)
    determinant_error = max(
        abs(determinants[left] - determinants[right])
        for left, right in (("M01", "M12"), ("M12", "M23"), ("M01", "M23"))
    ) + mp.mpf("1e-120")
    gate = 100 * determinant_error
    passed_minor = min(abs(value) for value in determinants.values()) > gate
    return {
        "rank": rank,
        "rows": rows,
        "columns": columns,
        "determinants": {
            level: mp_text(value) for level, value in determinants.items()
        },
        "error": mp_text(determinant_error),
        "gate": mp_text(gate),
        "passed": bool(passed_minor),
    }


def upper_triangular_source_change(size):
    transform = mp.eye(size)
    for index in range(size):
        transform[index, index] = 1 + mp.mpf(index + 1) / 1000
        if index + 1 < size:
            transform[index, index + 1] = mp.mpf(1) / 100
        if index + 3 < size:
            transform[index, index + 3] = mp.mpf(1) / 200
    return transform


def analyze_block(raw_matrices, dimension):
    scale_count = 5 * dimension
    scale_norm = column_block_frobenius(
        raw_matrices["M12"], 0, scale_count
    )
    strut_norm = column_block_frobenius(
        raw_matrices["M12"], scale_count, 10 * dimension
    )
    scale_factor = 1 / max(mp.mpf(1), scale_norm)
    strut_factor = 1 / max(mp.mpf(1), strut_norm)
    matrices = {
        level: matrix_scale_columns(
            matrix, scale_count, scale_factor, strut_factor
        )
        for level, matrix in raw_matrices.items()
    }
    census = triplet_nullity(matrices)
    minor = (
        deterministic_minor(matrices, census["rank"])
        if census["resolved"] else {
            "rank": None,
            "rows": [],
            "columns": [],
            "determinants": {},
            "error": None,
            "gate": None,
            "passed": False,
        }
    )
    transform = upper_triangular_source_change(matrices["M12"].cols)
    changed = {
        level: matrix * transform for level, matrix in matrices.items()
    }
    changed_census = triplet_nullity(changed)
    upper_ok = bool(
        census["resolved"]
        and changed_census["resolved"]
        and census["nullity"] == changed_census["nullity"]
    )
    zero_labels, _ = classify_values(
        [mp.mpf(0)] * matrices["M12"].cols,
        census["normalization"], census["total_error"],
    )
    injection_labels, _ = classify_values(
        [mp.mpf(1)] * matrices["M12"].cols,
        census["normalization"], census["total_error"],
    )
    synthetic_zero = zero_labels == ["ZERO"] * matrices["M12"].cols
    synthetic_injection = (
        injection_labels == ["NONZERO"] * matrices["M12"].cols
    )
    public = {
        "shape": [matrices["M12"].rows, matrices["M12"].cols],
        "source_scaling": {
            "scale_norm": mp_text(scale_norm),
            "strut_norm": mp_text(strut_norm),
            "scale_factor": mp_text(scale_factor),
            "strut_factor": mp_text(strut_factor),
        },
        "normalization": mp_text(census["normalization"]),
        "step_error": mp_text(census["step_error"]),
        "total_error": mp_text(census["total_error"]),
        "resolved": census["resolved"],
        "rank": census["rank"],
        "nullity": census["nullity"],
        "levels": census["levels"],
        "minor": minor,
        "hostile_controls": {
            "synthetic_zero_full_nullity": synthetic_zero,
            "synthetic_injection_zero_nullity": synthetic_injection,
            "upper_triangular_resolved": changed_census["resolved"],
            "upper_triangular_nullity": changed_census["nullity"],
            "upper_triangular_rank_invariant": upper_ok,
        },
    }
    numerical_ok = bool(census["resolved"] and minor["passed"])
    hostile_ok = bool(synthetic_zero and synthetic_injection and upper_ok)
    return public, matrices, numerical_ok, hostile_ok


def copy_sparse_rows(rows):
    return [dict(row) for row in rows]


def corrupt_carrier(carrier):
    rows = copy_sparse_rows(carrier["rows"])
    orbit_position = carrier["first_diagonal_orbit"] - 30
    touched = []
    for group in range(24):
        row = 24 * orbit_position + group
        scale_columns = sorted(column for column in rows[row] if column < VERTICES)
        if not scale_columns:
            raise RuntimeError("selected diagonal carrier row has no scale coefficient")
        column = scale_columns[0]
        rows[row][column] = rows[row].get(column, 0) + mp.mpf(1) / 10
        touched.append((row, column))
    return rows, touched


def corrupt_hessian_response(rows, carrier):
    result = copy_sparse_rows(rows)
    orbit_position = carrier["first_diagonal_orbit"] - 30
    touched = []
    for group in range(24):
        row = 24 * orbit_position + group
        for column, value in carrier["rows"][row].items():
            result[row][column] = (
                result[row].get(column, 0) + mp.mpf(1) * value / 10
            )
        touched.append(row)
    return result, touched


def right_kernel_projector(matrix, nullity):
    if nullity == 0:
        return np.zeros((matrix.shape[1], matrix.shape[1]), dtype=np.complex128)
    _, _, right = np.linalg.svd(matrix, full_matrices=False)
    vectors = right[-nullity:, :].conj().T
    return vectors @ vectors.conj().T


def projector_comparison(global_levels, nullities):
    if any(value is None for value in nullities.values()):
        return {
            "classification": "OPEN",
            "nullities": nullities,
            "reason": "at least one sector census is unresolved",
            "difference_two_norm": None,
            "uncertainty": None,
        }, None
    if nullities["even"] == nullities["odd"] == 0:
        return {
            "classification": "ZERO_KERNEL",
            "nullity": 0,
            "difference_two_norm": "0.0",
            "uncertainty": "0.0",
        }, None
    if nullities["even"] != nullities["odd"]:
        return {
            "classification": "DEPENDENT",
            "nullities": nullities,
            "difference_two_norm": None,
            "uncertainty": None,
        }, None
    nullity = nullities["even"]
    projectors = {}
    uncertainty = 0.0
    for parity in ("even", "odd"):
        projectors[parity] = {
            level: right_kernel_projector(
                global_levels[parity][level], nullity
            )
            for level in ("M01", "M12", "M23")
        }
        uncertainty = max(
            uncertainty,
            *(
                float(np.linalg.norm(
                    projectors[parity][left] - projectors[parity][right],
                    ord=2,
                ))
                for left, right in (
                    ("M01", "M12"), ("M12", "M23"), ("M01", "M23")
                )
            ),
        )
    uncertainty += 500 * np.finfo(float).eps * max(INTERNAL, DATA)
    difference = float(np.linalg.norm(
        projectors["even"]["M12"] - projectors["odd"]["M12"], ord=2
    ))
    if difference <= 10 * uncertainty:
        classification = "AGREE"
    elif difference > 100 * uncertainty:
        classification = "DEPENDENT"
    else:
        classification = "OPEN"
    return {
        "classification": classification,
        "nullity": nullity,
        "difference_two_norm": f"{difference:.17e}",
        "uncertainty": f"{uncertainty:.17e}",
        "agreement_gate": f"{10 * uncertainty:.17e}",
        "dependence_gate": f"{100 * uncertainty:.17e}",
    }, projectors


def relative_sparse_change(reference, changed, selected_rows):
    return sparse_difference_frobenius(
        reference, changed, selected_rows
    ) / max(mp.mpf(1), sparse_frobenius(reference, selected_rows))


def json_counter(counter):
    return {str(key): int(value) for key, value in sorted(counter.items())}


print("[setup] checking frozen provenance and registry", flush=True)
primary = json.loads(PRIMARY_JSON.read_text())
adversarial = json.loads(ADVERSARIAL_JSON.read_text())
input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
provenance_ok = bool(input_hashes == EXPECTED_HASHES)
check("all frozen inputs retain their preregistered hashes", provenance_ok)

scripts, duplicate_registry_entries = registry_inventory(RUN_ALL)
registry_ok = bool(
    scripts.count(VERIFIER_NAME) == 1 and not duplicate_registry_entries
)
check(
    "the verifier is registered exactly once and the registry has no duplicates",
    registry_ok,
    f"entries={len(scripts)}, duplicates={duplicate_registry_entries}",
)

accepted_inputs_ok = bool(
    primary.get("outcome") == "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_PRIMARY"
    and primary.get("tests") == primary.get("passed") == 22
    and adversarial.get("outcome")
    == "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_ADVERSARIALLY_REPLICATED"
    and adversarial.get("tests") == adversarial.get("passed") == 18
    and all(
        mp.mpf(primary["parities"][parity]["gradient"]["internal_maximum"])
        < mp.mpf("1e-25")
        for parity in ("even", "odd")
    )
)
check(
    "the accepted quadratic artifacts and stationary internal gradients remain intact",
    accepted_inputs_ok,
)

state = finite_height_formula_control(primary)
background_ok = bool(
    abs(state["elimination"]) < mp.mpf("1e-70")
    and state["h_error"] < mp.mpf("1e-70")
    and state["lambda_error"] < mp.mpf("1e-70")
    and state["rho_error"] < mp.mpf("1e-70")
    and state["mass_error"] < mp.mpf("1e-70")
    and state["h"] > 0
    and state["lambda"] > 0
    and state["rho"] > 0
)
check(
    "the finite-height background independently satisfies the exact homogeneous formulas",
    background_ok,
    (
        f"E={mp_text(state['elimination'], 6)}, "
        f"h_error={mp_text(state['h_error'], 6)}"
    ),
)

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_internal_rank", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
geometry_import_ok = bool(gro.tests == gro.passed == 43)
check(
    "the complete one-slab geometry retains all 43 upstream certificates",
    geometry_import_ok,
)

library = {
    "ARITHMETIC_FLOOR": ARITHMETIC_FLOOR,
    "Counter": Counter,
    "DPS": DPS,
    "DERIVATIVE_STEPS": DERIVATIVE_STEPS,
    "I": I,
    "LOCAL_EDGES": LOCAL_EDGES,
    "LOCAL_HINGES": LOCAL_HINGES,
    "LOCAL_HINGE_INDEX": LOCAL_HINGE_INDEX,
    "L0_SQUARE": mp.mpf(1),
    "MASS": state["mass"],
    "RHO0": state["rho"],
    "VARIANTS": VARIANTS,
    "cluster_sorted": None,
    "combinations": combinations,
    "defaultdict": defaultdict,
    "gro": gro,
    "math": math,
    "mp": mp,
    "mp_frobenius": mp_frobenius,
    "mp_submatrix": None,
    "np": np,
}
load_named_functions(
    RANK_SOURCE,
    {
        "log_minus",
        "signed_volume_square",
        "angle_record",
        "area_data",
        "extended_edge_image",
        "orbit_sort_key",
        "augment_boundary_orbits",
        "group_and_index_data",
        "prepare_geometry",
    },
    library,
)
load_named_functions(
    HIGH_SOURCE,
    {
        "mp_submatrix",
        "cluster_sorted",
        "high_precision_sector_bases",
        "high_precision_pattern_cache",
        "assemble_full_representative_kernels",
    },
    library,
)
models = {
    parity: library["augment_boundary_orbits"](model)
    for parity, model in gro.models.items()
}

records = {}
runtime = {}
for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing orbit kernels and responses", flush=True)
    model = models[parity]
    index_data = library["group_and_index_data"](
        model, (mp.log(state["lambda"]), mp.mpf(0))
    )
    geometry = library["prepare_geometry"](model, index_data)
    carrier = build_sparse_carrier(
        model, index_data, state["lambda"], state["rho"]
    )
    carrier_geometry_ok = bool(
        len(model["slab"]) == 2400
        and len(index_data["edge_to_index"]) == FULL
        and len(geometry["patterns"]) == 20
        and carrier["controls_ok"]
    )
    check(
        f"{parity}: geometry and exact sparse rank-240 carrier pass",
        carrier_geometry_ok,
        f"support={carrier['controls']['support']}",
    )

    sectors, sector_control = library["high_precision_sector_bases"](
        index_data
    )
    sector_ok = bool(
        sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        and sum(item["dimension"] ** 2 for item in sectors) == 24
        and all(
            sector_control[f"maximum_{key}"] < mp.mpf("1e-140")
            for key in (
                "orthonormal", "center", "splitter", "right_leakage",
                "conjugate_pair",
            )
        )
    )
    check(
        f"{parity}: seven deterministic minimal 2T sectors pass exact controls",
        sector_ok,
        f"dimensions={sector_control['irrep_dimensions']}",
    )

    kind_values = {
        "old": mp.mpf(1),
        "internal": state["lambda"] - state["rho"],
        "pole": -state["rho"],
        "new": state["lambda"] ** 2,
    }
    pattern_cache, pattern_control = library["high_precision_pattern_cache"](
        geometry["patterns"], kind_values
    )
    pattern_ok = bool(
        pattern_control["entry_pass"]
        and pattern_control["base_negative_counts"] == Counter({1: 2400})
        and set(pattern_control["displaced_negative_counts"]) == {1}
        and pattern_control["minimum_leading_minor"] > 0
        and pattern_control["minimum_argument"] > mp.mpf("1e-6")
    )
    check(
        f"{parity}: all derivative levels retain the Lorentzian branch and step hierarchy",
        pattern_ok,
        (
            f"cross={mp_text(pattern_control['maximum_cross'], 6)}, "
            f"proxy={mp_text(pattern_control['maximum_proxy'], 6)}"
        ),
    )

    kernels, kernel_control = library["assemble_full_representative_kernels"](
        index_data, geometry, pattern_cache
    )
    kernel_ok = bool(
        set(kernels) == set(VARIANTS)
        and all(len(kernel) > 0 for kernel in kernels.values())
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-140")
    )
    check(
        f"{parity}: physical identity-row Hessian kernels are real to the frozen gate",
        kernel_ok,
        f"imag={mp_text(kernel_control['maximum_imaginary'], 6)}",
    )

    raw_responses = {}
    compose_controls = {}
    for name in VARIANTS:
        raw_responses[name], compose_controls[name] = compose_response(
            kernels[name], index_data, carrier["rows"]
        )
    responses = {
        "M01": sparse_richardson(
            raw_responses["operational_primary"],
            raw_responses["operational_shadow"],
        ),
        "M12": sparse_richardson(
            raw_responses["operational_shadow"],
            raw_responses["validation_primary"],
        ),
        "M23": sparse_richardson(
            raw_responses["validation_primary"],
            raw_responses["validation_shadow"],
        ),
    }
    response_norm = max(mp.mpf(1), sparse_frobenius(responses["M12"]))
    response_imaginary = max(
        sparse_frobenius(rows, imaginary=True) / response_norm
        for rows in responses.values()
    )
    response_ok = bool(response_imaginary < mp.mpf("1e-140"))
    check(
        f"{parity}: complete active responses have negligible physical imaginary residue",
        response_ok,
        f"relative_imag={mp_text(response_imaginary, 6)}",
    )

    pole_positions = tuple(sorted(
        orbit_type - 30 for orbit_type in carrier["pole_types"]
    ))
    diagonal_positions = tuple(
        position for position in range(35) if position not in pole_positions
    )
    split_ok = bool(
        len(pole_positions) == 5
        and len(diagonal_positions) == 30
        and set(pole_positions).isdisjoint(diagonal_positions)
    )
    check(
        f"{parity}: physical labels give the frozen 720 diagonal plus 120 pole row split",
        split_ok,
        f"pole_positions={pole_positions}",
    )

    sector_records = []
    block_runtime = []
    parity_blocks_resolved = True
    parity_block_controls_ok = True
    for sector_index, sector in enumerate(sectors):
        dimension = int(sector["dimension"])
        item = {
            "sector_index": sector_index,
            "dimension": dimension,
            "center_value": mp_text(sector["center_value"]),
            "constant_overlap": mp_text(sector["constant_overlap"]),
            "scopes": {},
        }
        runtime_item = {}
        for scope, positions in (
            ("diagonal", diagonal_positions), ("full", tuple(range(35)))
        ):
            projected = {
                level: project_sparse_response(
                    rows, positions, carrier["data_orbit_group"], sector
                )
                for level, rows in responses.items()
            }
            block_record, scaled, numerical_ok, hostile_ok = analyze_block(
                projected, dimension
            )
            item["scopes"][scope] = block_record
            runtime_item[scope] = {
                "raw": projected,
                "scaled": scaled,
            }
            parity_blocks_resolved &= numerical_ok
            parity_block_controls_ok &= hostile_ok
        sector_records.append(item)
        block_runtime.append(runtime_item)
        print(
            f"[{parity}] sector {sector_index} d={dimension}: "
            f"diag nullity={item['scopes']['diagonal']['nullity']}, "
            f"full nullity={item['scopes']['full']['nullity']}",
            flush=True,
        )
    check(
        f"{parity}: every sector classifier passes its synthetic and source-change controls",
        parity_block_controls_ok,
    )

    global_levels = {
        level: sparse_to_numpy(rows, list(range(INTERNAL)))
        for level, rows in responses.items()
    }
    diagonal_rows = [
        24 * position + group
        for position in diagonal_positions for group in range(24)
    ]
    global_diagonal = {
        level: sparse_to_numpy(rows, diagonal_rows)
        for level, rows in responses.items()
    }

    kernel_r12 = combine_kernels(
        kernels["operational_shadow"], kernels["validation_primary"]
    )
    corrupted_carrier_rows, carrier_touched = corrupt_carrier(carrier)
    corrupted_response, _ = compose_response(
        kernel_r12, index_data, corrupted_carrier_rows
    )
    carrier_corruption_effect = relative_sparse_change(
        responses["M12"], corrupted_response, range(INTERNAL)
    )
    hessian_corrupted, hessian_touched = corrupt_hessian_response(
        responses["M12"], carrier
    )
    hessian_corruption_effect = relative_sparse_change(
        responses["M12"], hessian_corrupted, range(INTERNAL)
    )
    wrong_response, _ = compose_response(
        kernel_r12, index_data, carrier["rows"], reverse=True
    )
    wrong_convention_effect = relative_sparse_change(
        responses["M12"], wrong_response, range(INTERNAL)
    )
    corruption_ok = bool(
        carrier_corruption_effect > mp.mpf("1e-12")
        and hessian_corruption_effect > mp.mpf("1e-12")
    )
    check(
        f"{parity}: frozen carrier and Hessian corruptions are detected",
        corruption_ok,
        (
            f"carrier={mp_text(carrier_corruption_effect, 6)}, "
            f"Hdiag={mp_text(hessian_corruption_effect, 6)}"
        ),
    )

    full_nullity = None
    diagonal_nullity = None
    if parity_blocks_resolved:
        full_nullity = sum(
            record["dimension"] * record["scopes"]["full"]["nullity"]
            for record in sector_records
        )
        diagonal_nullity = sum(
            record["dimension"]
            * record["scopes"]["diagonal"]["nullity"]
            for record in sector_records
        )

    records[parity] = {
        "carrier": carrier["controls"],
        "sector_controls": {
            "isotypic_dimensions": sector_control["isotypic_dimensions"],
            "irrep_dimensions": sector_control["irrep_dimensions"],
            **{
                key: mp_text(value)
                for key, value in sector_control.items()
                if key.startswith("maximum_")
            },
        },
        "pattern": {
            "minimum_leading_minor": mp_text(
                pattern_control["minimum_leading_minor"]
            ),
            "minimum_argument": mp_text(pattern_control["minimum_argument"]),
            "maximum_cross": mp_text(pattern_control["maximum_cross"]),
            "maximum_proxy": mp_text(pattern_control["maximum_proxy"]),
            "base_negative_counts": json_counter(
                pattern_control["base_negative_counts"]
            ),
            "displaced_negative_counts": json_counter(
                pattern_control["displaced_negative_counts"]
            ),
            "raw_boost_imaginary_diagnostic": mp_text(
                pattern_control["maximum_raw_angle_or_derivative_imaginary"]
            ),
        },
        "kernel": {
            "nonzero_entries": kernel_control["nonzero_entries"],
            "maximum_imaginary": mp_text(kernel_control["maximum_imaginary"]),
            "response_relative_imaginary": mp_text(response_imaginary),
        },
        "row_split": {
            "pole_orbit_positions": list(pole_positions),
            "diagonal_orbit_positions": list(diagonal_positions),
        },
        "sectors": sector_records,
        "global_nullities": {
            "diagonal": diagonal_nullity,
            "full": full_nullity,
        },
        "hostile_controls": {
            "carrier_corruption_relative_effect": mp_text(
                carrier_corruption_effect
            ),
            "carrier_corruption_touched": [list(item) for item in carrier_touched],
            "hessian_diagonal_corruption_relative_effect": mp_text(
                hessian_corruption_effect
            ),
            "hessian_diagonal_touched_rows": hessian_touched,
            "wrong_group_product_diagnostic": mp_text(wrong_convention_effect),
        },
    }
    runtime[parity] = {
        "carrier": carrier,
        "responses": responses,
        "global_levels": global_levels,
        "global_diagonal": global_diagonal,
        "blocks_resolved": parity_blocks_resolved,
        "block_controls_ok": parity_block_controls_ok,
        "control_ok": bool(
            carrier_geometry_ok and sector_ok and pattern_ok and kernel_ok
            and response_ok and split_ok and parity_block_controls_ok
            and corruption_ok
        ),
    }

even_signature = records["even"]["sectors"]
odd_signature = records["odd"]["sectors"]
sector_signatures_equal = bool(
    len(even_signature) == len(odd_signature)
    and all(
        left["dimension"] == right["dimension"]
        and abs(mp.mpc(left["center_value"]) - mp.mpc(right["center_value"]))
        < mp.mpf("1e-140")
        and abs(
            mp.mpc(left["constant_overlap"])
            - mp.mpc(right["constant_overlap"])
        ) < mp.mpf("1e-140")
        for left, right in zip(even_signature, odd_signature)
    )
)
check(
    "even and odd reconstructions retain the same deterministic sector signature",
    sector_signatures_equal,
)

global_nullities = {
    parity: records[parity]["global_nullities"]["full"]
    for parity in ("even", "odd")
}
projector_record, projectors = projector_comparison(
    {
        parity: runtime[parity]["global_levels"]
        for parity in ("even", "odd")
    },
    global_nullities,
)
check(
    "the frozen physical-data kernel comparison follows its declared classifier",
    projector_record["classification"]
    in {"ZERO_KERNEL", "AGREE", "DEPENDENT", "OPEN"},
    str(projector_record),
)

print("[closure] reading accepted quadratic matrices only after rank census", flush=True)
accepted_quadratics = np.load(ADVERSARIAL_MATRICES, allow_pickle=False)
quadratic_errors = {}
for parity_index, parity in enumerate(("even", "odd")):
    reconstructed = quadratic_from_response(
        runtime[parity]["carrier"]["rows"],
        runtime[parity]["responses"]["M12"],
    )
    reconstructed_numpy = mp_to_numpy(reconstructed)
    quadratic_errors[parity] = float(
        np.linalg.norm(
            reconstructed_numpy - accepted_quadratics[parity_index], ord="fro"
        )
        / max(1.0, np.linalg.norm(accepted_quadratics[parity_index], ord="fro"))
    )
quadratic_closure_ok = max(quadratic_errors.values()) < 1e-10
check(
    "the response reconstruction closes on both accepted quadratic matrices",
    quadratic_closure_ok,
    str(quadratic_errors),
)

all_controls = bool(
    provenance_ok and registry_ok and accepted_inputs_ok and background_ok
    and geometry_import_ok and sector_signatures_equal and quadratic_closure_ok
    and all(runtime[parity]["control_ok"] for parity in ("even", "odd"))
)
actual_blocks_resolved = bool(
    all(runtime[parity]["blocks_resolved"] for parity in ("even", "odd"))
)
if not all_controls:
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_CONTROL_FAILED"
elif not actual_blocks_resolved or projector_record["classification"] == "OPEN":
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_NUMERICALLY_OPEN"
elif global_nullities == {"even": 0, "odd": 0}:
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_FULL_COLUMN_RANK_PRIMARY"
elif (
    global_nullities["even"] == global_nullities["odd"]
    and global_nullities["even"] > 0
    and projector_record["classification"] == "AGREE"
):
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_SELECTED_PRIMARY"
else:
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_SCHEDULE_DEPENDENT_PRIMARY"

check(
    "the primary outcome follows the frozen hierarchy",
    outcome in {
        "FINITE_HEIGHT_INTERNAL_CARRIER_CONTROL_FAILED",
        "FINITE_HEIGHT_INTERNAL_CARRIER_NUMERICALLY_OPEN",
        "FINITE_HEIGHT_INTERNAL_CARRIER_FULL_COLUMN_RANK_PRIMARY",
        "FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_SELECTED_PRIMARY",
        "FINITE_HEIGHT_INTERNAL_CARRIER_SCHEDULE_DEPENDENT_PRIMARY",
    },
    outcome,
)

matrix_payload = {
    "even_R12_full": runtime["even"]["global_levels"]["M12"],
    "odd_R12_full": runtime["odd"]["global_levels"]["M12"],
    "even_R12_diagonal": runtime["even"]["global_diagonal"]["M12"],
    "odd_R12_diagonal": runtime["odd"]["global_diagonal"]["M12"],
}
if projectors is not None:
    matrix_payload["even_kernel_projector"] = projectors["even"]["M12"]
    matrix_payload["odd_kernel_projector"] = projectors["odd"]["M12"]
np.savez_compressed(MATRIX_OUTPUT, **matrix_payload)
matrix_hash = sha256(MATRIX_OUTPUT)

artifact = {
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "status": "PRIMARY_ONLY_REPLICATION_REQUIRED_FOR_ANY_MATERIAL_OUTCOME",
    "method": (
        "180-digit identity-row 2T Hessian kernels composed with the exact "
        "rank-240 finite-height carrier, seven minimal sectors, singular "
        "intervals, and direct determinant minors"
    ),
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "input_sha256": input_hashes,
    },
    "background": {
        key: mp_text(value) for key, value in state.items()
    },
    "global_nullities": global_nullities,
    "kernel_comparison": projector_record,
    "quadratic_closure_relative_errors": {
        key: f"{value:.17e}" for key, value in quadratic_errors.items()
    },
    "parities": records,
    "matrices": {
        "path": MATRIX_OUTPUT.name,
        "sha256": matrix_hash,
        "arrays": {
            key: list(value.shape) for key, value in matrix_payload.items()
        },
    },
    "interpretation": {
        "scope": (
            "one exact rank-240 scale-plus-strut carrier at one finite-height "
            "600-cell slab background"
        ),
        "full_column_rank_meaning": (
            "bounded negative for this carrier and background only"
        ),
        "finite_step_invariant_continuation": "SEPARATE_DERIVED_RESULT",
        "infinite_proper_time_evolution": "OPEN",
        "gravitons": "NOT_DERIVED",
        "wave_equation": "NOT_DERIVED",
        "tick_c_G_planck_particle_masses": "NOT_DERIVED",
    },
    "firewall": {
        "accepted_quadratic_matrices_read_after_rank_census": True,
        "continuum_target_parsed": False,
        "desired_rank_or_nullity_parsed": False,
        "full_suite_run": False,
        "binary64_singular_values_used_for_rank_classification": False,
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"OUTCOME: {outcome}")
print(f"RESULT: {passed}/{tests} PASS")
print(f"MATRIX SHA: {matrix_hash}")
if passed != tests:
    raise SystemExit(1)

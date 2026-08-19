#!/usr/bin/env python3
"""Exact tangent-admissibility audit for schedule-free 600-cell frusta."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import networkx as nx
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_canonical_data_admissibility.json"
PRIOR_ART = (
    ROOT / "docs/gravity/gravity_600cell_canonical_data_admissibility_prior_art.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_canonical_data_admissibility_protocol.md"
)
GLOBAL_PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_global_variable_face_closure_protocol.md"
)
GLOBAL_SOURCE = HERE / "verify_gravity_600cell_global_variable_face_closure.py"
GLOBAL_JSON = HERE / "gravity_600cell_global_variable_face_closure.json"
LOCAL_RESULT = ROOT / "docs/gravity/gravity_600cell_variable_face_connection_result.md"
LOCAL_SOURCE = HERE / "verify_gravity_600cell_variable_face_connection.py"
LOCAL_JSON = HERE / "gravity_600cell_variable_face_connection.json"

PROTOCOL_COMMIT = "4b007c2"
EXPECTED_HASHES = {
    "prior_art": "8acf8e29a809f065033d49dd03c1e858f1ed7a219d1be6acfb27851f78f1ce56",
    "protocol": "8db29cb9af699da660b969988eeb76c5e605e67c5ec65716795ada2e34674185",
    "global_protocol": "ed79c6a15ade377ae09854b3cad3028eb1c0f43cf8e85789d46993fe25ed1b49",
    "global_source": "ec44be8e4d82634e30944739d10d3f80fbb9f6fee0883ec1f612690c38d90ab6",
    "global_json": "61cebd1cd67fcdc56de088855b1fc7b805d0f70f9f9b3029d4a61209d7a53944",
    "local_result": "2db55cb87ec1c01d537cdbc11010bc9ea740762c598108e4c2de0f3acca72cc8",
    "local_source": "69a5d7479a5df427cead76f82db31fe62a9190c28c967f699c846881634fb0f6",
    "local_json": "001212016553d006862e68edc4f780f37ca1476110b6e0aed3e987f52a43b5e3",
}

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
CANONICAL_BASE = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
))
REPRESENTATIVES = ((2, 5), (3, 11))
PRIMES = (1000003, 1000033)
EDGE_COUNT = 720
VERTEX_COUNT = 120
CELL_COLUMNS = 3600
DATA_COLUMNS = EDGE_COUNT + VERTEX_COUNT
TOTAL_COLUMNS = CELL_COLUMNS + DATA_COLUMNS

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def kernel(matrix):
    vectors = matrix.nullspace()
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(matrix.cols, 0)


def same_space(left, right):
    return bool(
        left.rows == right.rows
        and left.rank() == right.rank()
        and left.row_join(right).rank() == left.rank()
    )


# Exact Q(phi), represented by a+b*phi with phi^2=phi+1.
def golden(a=0, b=0):
    return (Fraction(a), Fraction(b))


def gadd(left, right):
    return (left[0] + right[0], left[1] + right[1])


def gmul(left, right):
    a, b = left
    c, d = right
    return (a * c + b * d, a * d + b * c + b * d)


def gscale(value, scalar):
    scalar = Fraction(scalar)
    return (scalar * value[0], scalar * value[1])


def gdot(left, right):
    result = golden()
    for a, b in zip(left, right):
        result = gadd(result, gmul(a, b))
    return result


def permutation_parity(permutation):
    return sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    ) % 2


def build_exact_complex():
    zero = golden()
    one = golden(1)
    half = golden(Fraction(1, 2))
    phi_half = golden(0, Fraction(1, 2))
    inverse_phi_half = golden(Fraction(-1, 2), Fraction(1, 2))
    vertices = set()
    for axis in range(4):
        for sign in (-1, 1):
            point = [zero] * 4
            point[axis] = gscale(one, sign)
            vertices.add(tuple(point))
    for signs in product((-1, 1), repeat=4):
        vertices.add(tuple(gscale(half, sign) for sign in signs))
    base = (phi_half, half, inverse_phi_half, zero)
    for permutation in permutations(range(4)):
        if permutation_parity(permutation):
            continue
        permuted = tuple(base[permutation[index]] for index in range(4))
        nonzero = tuple(
            index for index, value in enumerate(permuted) if value != zero
        )
        for signs in product((-1, 1), repeat=3):
            point = list(permuted)
            for index, sign in zip(nonzero, signs):
                point[index] = gscale(point[index], sign)
            vertices.add(tuple(point))
    vertices = tuple(sorted(vertices))

    adjacency = [set() for _ in vertices]
    edges = []
    for left, right in combinations(range(len(vertices)), 2):
        if gdot(vertices[left], vertices[right]) == phi_half:
            adjacency[left].add(right)
            adjacency[right].add(left)
            edges.append((left, right))
    triangles = []
    for left, right in edges:
        for third in adjacency[left] & adjacency[right]:
            if right < third:
                triangles.append((left, right, third))
    tetrahedra = []
    for a, b, c in triangles:
        for d in adjacency[a] & adjacency[b] & adjacency[c]:
            if c < d:
                tetrahedra.append((a, b, c, d))
    tetrahedra = tuple(sorted(tetrahedra))

    face_to_tetrahedra = {}
    edge_to_tetrahedra = {}
    for tetrahedron in tetrahedra:
        for face in combinations(tetrahedron, 3):
            face_to_tetrahedra.setdefault(tuple(sorted(face)), []).append(
                tetrahedron
            )
        for edge in combinations(tetrahedron, 2):
            edge_to_tetrahedra.setdefault(tuple(sorted(edge)), []).append(
                tetrahedron
            )
    return {
        "vertices": vertices,
        "edges": tuple(edges),
        "triangles": tuple(triangles),
        "tetrahedra": tetrahedra,
        "face_to_tetrahedra": face_to_tetrahedra,
        "edge_to_tetrahedra": edge_to_tetrahedra,
    }


def lorentz_basis(metric):
    result = []
    for a, b in ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)):
        generator = sp.zeros(4)
        generator[a, b] = 1
        generator[b, a] = -metric[a, a] / metric[b, b]
        result.append(generator)
    return tuple(result)


def local_poincare_coefficients(scale, lapse, basis, normal):
    if scale == 1:
        result = sp.zeros(10, 6)
        result[:3, :3] = sp.eye(3)
        result[6:9, 3:6] = sp.eye(3)
        return result
    result = sp.zeros(10, 6)
    result[:6, :6] = sp.eye(6)
    for column, generator in enumerate(basis):
        result[6:10, column] = (
            sp.Rational(lapse, scale - 1) * generator * normal
        )
    return result


def poincare_evaluation(points, basis):
    columns = [
        sp.Matrix.vstack(*(generator * point for point in points))
        for generator in basis
    ]
    for axis in range(4):
        direction = sp.eye(4)[:, axis]
        columns.append(sp.Matrix.vstack(*(direction for _ in points)))
    return sp.Matrix.hstack(*columns)


def coordinate_map(tetrahedron, canonical):
    return {
        vertex: canonical[index] for index, vertex in enumerate(tetrahedron)
    }


def affine_reflection(point, face_points, metric):
    anchor = face_points[0]
    tangents = sp.Matrix.hstack(*(value - anchor for value in face_points[1:]))
    normals = (tangents.T * metric).nullspace()
    if len(normals) != 1:
        raise RuntimeError("lateral face does not have a unique normal")
    normal = normals[0]
    denominator = (normal.T * metric * normal)[0]
    if denominator == 0:
        raise RuntimeError("lateral face normal is null")
    return sp.simplify(
        point
        - 2
        * normal
        * (normal.T * metric * (point - anchor))[0]
        / denominator
    )


def affine_map(domain, codomain):
    domain_h = sp.Matrix.hstack(
        *(point.col_join(sp.ones(1, 1)) for point in domain)
    )
    codomain_h = sp.Matrix.hstack(
        *(point.col_join(sp.ones(1, 1)) for point in codomain)
    )
    return sp.simplify(codomain_h * domain_h.inv())


def lateral_transition(source, target, shared, canonical, scale, lapse, metric):
    source_coordinates = coordinate_map(source, canonical)
    target_coordinates = coordinate_map(target, canonical)
    source_top = {
        vertex: scale * point + lapse * NORMAL
        for vertex, point in source_coordinates.items()
    }
    target_top = {
        vertex: scale * point + lapse * NORMAL
        for vertex, point in target_coordinates.items()
    }
    source_apex = next(vertex for vertex in source if vertex not in shared)
    target_apex = next(vertex for vertex in target if vertex not in shared)

    source_lower_face = tuple(source_coordinates[vertex] for vertex in shared)
    lateral_points = source_lower_face + (source_top[shared[0]],)
    reflected_lower = affine_reflection(
        source_coordinates[source_apex], lateral_points, metric
    )
    reflected_upper = affine_reflection(
        source_top[source_apex], lateral_points, metric
    )
    domain = (
        tuple(target_coordinates[vertex] for vertex in shared)
        + (target_top[shared[0]], target_coordinates[target_apex])
    )
    codomain = lateral_points + (reflected_lower,)
    transition = affine_map(domain, codomain)
    linear = transition[:4, :4]

    shared_mapping = all(
        sp.simplify(
            transition
            * target_coordinates[vertex].col_join(sp.ones(1, 1))
            - source_coordinates[vertex].col_join(sp.ones(1, 1))
        )
        == sp.zeros(5, 1)
        and sp.simplify(
            transition * target_top[vertex].col_join(sp.ones(1, 1))
            - source_top[vertex].col_join(sp.ones(1, 1))
        )
        == sp.zeros(5, 1)
        for vertex in shared
    )
    apex_mapping = bool(
        sp.simplify(
            transition * target_top[target_apex].col_join(sp.ones(1, 1))
            - reflected_upper.col_join(sp.ones(1, 1))
        )
        == sp.zeros(5, 1)
    )
    control = bool(
        transition[4, :] == sp.Matrix([[0, 0, 0, 0, 1]])
        and sp.simplify(linear.T * metric * linear) == metric
        and abs(linear.det()) == 1
        and shared_mapping
        and apex_mapping
    )
    return {
        "transition": transition,
        "linear": linear,
        "control": control,
        "source_coordinates": source_coordinates,
        "target_coordinates": target_coordinates,
        "source_top": source_top,
        "target_top": target_top,
    }


def selector(tetrahedron, shared):
    result = sp.zeros(12, 16)
    for shared_index, vertex in enumerate(shared):
        local_index = tetrahedron.index(vertex)
        result[
            4 * shared_index : 4 * shared_index + 4,
            4 * local_index : 4 * local_index + 4,
        ] = sp.eye(4)
    return result


def local_length_geometry(scale, lapse, metric, canonical):
    top = tuple(scale * point + lapse * NORMAL for point in canonical)
    edge_pairs = tuple(combinations(range(4), 2))
    jacobian = sp.zeros(10, 16)
    for row, (left, right) in enumerate(edge_pairs):
        covector = 2 * (top[left] - top[right]).T * metric
        jacobian[row, 4 * left : 4 * left + 4] = covector
        jacobian[row, 4 * right : 4 * right + 4] = -covector
    for local_index in range(4):
        covector = (
            2 * (top[local_index] - canonical[local_index]).T * metric
        )
        jacobian[
            6 + local_index,
            4 * local_index : 4 * local_index + 4,
        ] = covector

    pivot_columns = tuple(jacobian.rref()[1])
    if len(pivot_columns) != 10:
        raise RuntimeError("local natural-length Jacobian is not rank ten")
    square = jacobian[:, pivot_columns]
    right_inverse = sp.zeros(16, 10)
    inverse = square.inv()
    for row, column in enumerate(pivot_columns):
        right_inverse[column, :] = inverse[row, :]
    local_kernel = kernel(jacobian)
    basis = lorentz_basis(metric)
    poincare_image = (
        poincare_evaluation(top, basis)
        * local_poincare_coefficients(scale, lapse, basis, NORMAL)
    )
    graph = sp.Matrix(
        6, 10, lambda row, column: (row + 1) * (column + 1)
    )
    alternate = right_inverse + local_kernel * graph
    controls = bool(
        jacobian.rank() == 10
        and jacobian * right_inverse == sp.eye(10)
        and local_kernel.shape == (16, 6)
        and same_space(local_kernel, poincare_image)
        and jacobian * alternate == sp.eye(10)
    )
    return {
        "top": top,
        "edge_pairs": edge_pairs,
        "jacobian": jacobian,
        "kernel": local_kernel,
        "right_inverse": right_inverse,
        "alternate_right_inverse": alternate,
        "controls": controls,
        "pivot_columns": pivot_columns,
    }


def local_data_keys(tetrahedron, edge_index, corruption=None):
    result = []
    for row, (left, right) in enumerate(combinations(range(4), 2)):
        edge = tuple(sorted((tetrahedron[left], tetrahedron[right])))
        key = edge_index[edge]
        if corruption is not None and row == corruption:
            key = DATA_COLUMNS
        result.append(key)
    result.extend(EDGE_COUNT + vertex for vertex in tetrahedron)
    return tuple(result)


def independent_row_indices(matrix):
    return tuple(matrix.T.rref()[1])


def face_equations(
    source,
    target,
    shared,
    transition_data,
    local_geometry,
    right_inverse,
    edge_index,
    source_offset,
    target_offset,
    corruption=None,
):
    local_kernel = local_geometry["kernel"]
    source_selector = selector(source, shared)
    target_selector = selector(target, shared)
    linear_block = sp.diag(*([transition_data["linear"]] * 3))
    fixed = (source_selector * local_kernel).row_join(
        -linear_block * target_selector * local_kernel
    )

    basis = lorentz_basis(ETA if transition_data["linear"].T * ETA * transition_data["linear"] == ETA else -ETA)
    source_lower = tuple(
        transition_data["source_coordinates"][vertex] for vertex in shared
    )
    source_upper = tuple(
        transition_data["source_top"][vertex] for vertex in shared
    )
    lower_evaluation = poincare_evaluation(source_lower, basis)
    stabilizer = kernel(lower_evaluation)
    source_upper_evaluation = poincare_evaluation(source_upper, basis)
    connection_column = source_upper_evaluation * stabilizer
    annihilator = kernel(connection_column.T).T
    fixed_eliminated = annihilator * fixed

    source_response = annihilator * source_selector * right_inverse
    target_response = (
        -annihilator * linear_block * target_selector * right_inverse
    )
    source_keys = local_data_keys(source, edge_index)
    target_keys = local_data_keys(target, edge_index, corruption=corruption)
    data_keys = tuple(sorted(set(source_keys + target_keys)))
    data_positions = {key: index for index, key in enumerate(data_keys)}
    compact = fixed_eliminated.row_join(
        sp.zeros(fixed_eliminated.rows, len(data_keys))
    )
    for local_row, key in enumerate(source_keys):
        column = 12 + data_positions[key]
        compact[:, column] += source_response[:, local_row]
    for local_row, key in enumerate(target_keys):
        column = 12 + data_positions[key]
        compact[:, column] += target_response[:, local_row]

    row_indices = independent_row_indices(compact)
    rows = []
    for row_index in row_indices:
        row = {}
        for local_column in range(6):
            value = sp.cancel(compact[row_index, local_column])
            if value:
                row[source_offset + local_column] = value
            value = sp.cancel(compact[row_index, 6 + local_column])
            if value:
                row[target_offset + local_column] = value
        for local_column, key in enumerate(data_keys):
            value = sp.cancel(compact[row_index, 12 + local_column])
            if value:
                row[CELL_COLUMNS + key] = value
        rows.append(row)

    control = bool(
        stabilizer.shape == (10, 1)
        and connection_column.shape == (12, 1)
        and connection_column.rank() == 1
        and annihilator.shape == (11, 12)
        and annihilator.rank() == 11
        and annihilator * connection_column == sp.zeros(11, 1)
        and fixed.row_join(-connection_column).rank() == 6
        and fixed_eliminated.rank() == 5
        and compact.rank() == len(row_indices)
    )
    return {
        "rows": rows,
        "control": control,
        "fixed_rank": fixed_eliminated.rank(),
        "augmented_rank": compact.rank(),
        "compact": compact,
        "data_keys": data_keys,
    }


def rational_mod(value, prime):
    value = sp.Rational(value)
    numerator = int(value.p) % prime
    denominator = int(value.q) % prime
    if denominator == 0:
        raise ZeroDivisionError(f"denominator vanishes modulo {prime}")
    return numerator * pow(denominator, -1, prime) % prime


def modular_rank(rows, columns, prime):
    pivots = {}
    maximum_width = 0
    for original in rows:
        row = {}
        for column, value in original.items():
            reduced = rational_mod(value, prime)
            if reduced:
                row[column] = reduced
        while row:
            existing = [column for column in row if column in pivots]
            if not existing:
                pivot = min(row)
                inverse = pow(row[pivot], -1, prime)
                row = {
                    column: value * inverse % prime
                    for column, value in row.items()
                    if value % prime
                }
                pivots[pivot] = row
                maximum_width = max(maximum_width, len(row))
                break
            pivot = min(existing)
            factor = row[pivot]
            pivot_row = pivots[pivot]
            for column, value in pivot_row.items():
                updated = (row.get(column, 0) - factor * value) % prime
                if updated:
                    row[column] = updated
                elif column in row:
                    del row[column]
        if len(pivots) == columns:
            break
    return len(pivots), maximum_width


def row_dot(row, vector):
    return sp.cancel(sum(value * vector.get(column, 0) for column, value in row.items()))


def solve_full_column(matrix, right_hand_side):
    pivot_rows = independent_row_indices(matrix)
    if len(pivot_rows) != matrix.cols:
        raise RuntimeError("matrix does not have full column rank")
    square = matrix[list(pivot_rows), :]
    result = square.inv() * right_hand_side[list(pivot_rows), :]
    if matrix * result != right_hand_side:
        raise RuntimeError("right hand side is outside the declared image")
    return result


def homothetic_vectors(local_geometry, right_inverse, tetrahedron_count):
    jacobian = local_geometry["jacobian"]
    local_kernel = local_geometry["kernel"]
    displacement_scale = sp.Matrix.vstack(*CANONICAL_BASE)
    displacement_lapse = sp.Matrix.vstack(*(NORMAL for _ in range(4)))
    vectors = []
    data_columns = []
    for name, displacement in (
        ("scale", displacement_scale),
        ("lapse", displacement_lapse),
    ):
        local_data = jacobian * displacement
        top_values = tuple(local_data[index] for index in range(6))
        strut_values = tuple(local_data[6 + index] for index in range(4))
        if len(set(top_values)) != 1 or len(set(strut_values)) != 1:
            raise RuntimeError("homothetic data are not globally constant")
        residual = displacement - right_inverse * local_data
        local_flex = solve_full_column(local_kernel, residual)
        vector = {}
        for tetrahedron_index in range(tetrahedron_count):
            for local_column in range(6):
                value = sp.cancel(local_flex[local_column])
                if value:
                    vector[6 * tetrahedron_index + local_column] = value
        for edge in range(EDGE_COUNT):
            value = sp.cancel(top_values[0])
            if value:
                vector[CELL_COLUMNS + edge] = value
        for vertex in range(VERTEX_COUNT):
            value = sp.cancel(strut_values[0])
            if value:
                vector[CELL_COLUMNS + EDGE_COUNT + vertex] = value
        vectors.append((name, vector))
        data_columns.append(
            sp.Matrix([top_values[0], strut_values[0]])
        )
    data_rank = sp.Matrix.hstack(*data_columns).rank()
    return vectors, data_rank


def build_global(
    complex_data,
    scale,
    lapse,
    metric,
    canonical,
    use_alternate=False,
    reverse_orientation=False,
):
    tetrahedra = complex_data["tetrahedra"]
    tetrahedron_index = {
        tetrahedron: index for index, tetrahedron in enumerate(tetrahedra)
    }
    faces = tuple(sorted(complex_data["face_to_tetrahedra"]))
    edge_index = {
        edge: index for index, edge in enumerate(sorted(complex_data["edges"]))
    }
    local_geometry = local_length_geometry(
        scale, lapse, metric, canonical
    )
    right_inverse = local_geometry[
        "alternate_right_inverse" if use_alternate else "right_inverse"
    ]
    rows = []
    transition_ok = True
    inverse_ok = True
    face_ok = True
    local_fixed_ranks = Counter()
    local_augmented_ranks = Counter()
    first_face_data = None

    for face_number, face in enumerate(faces):
        left, right = sorted(complex_data["face_to_tetrahedra"][face])
        source, target = (
            (right, left) if reverse_orientation else (left, right)
        )
        shared = tuple(sorted(face))
        forward = lateral_transition(
            source, target, shared, canonical, scale, lapse, metric
        )
        backward = lateral_transition(
            target, source, shared, canonical, scale, lapse, metric
        )
        transition_ok &= forward["control"] and backward["control"]
        inverse_ok &= bool(
            sp.simplify(forward["transition"] * backward["transition"])
            == sp.eye(5)
        )
        data = face_equations(
            source,
            target,
            shared,
            forward,
            local_geometry,
            right_inverse,
            edge_index,
            6 * tetrahedron_index[source],
            6 * tetrahedron_index[target],
        )
        face_ok &= data["control"]
        local_fixed_ranks[data["fixed_rank"]] += 1
        local_augmented_ranks[data["augmented_rank"]] += 1
        rows.extend(data["rows"])
        if face_number == 0:
            corrupted = face_equations(
                source,
                target,
                shared,
                forward,
                local_geometry,
                right_inverse,
                edge_index,
                6 * tetrahedron_index[source],
                6 * tetrahedron_index[target],
                corruption=0,
            )
            first_face_data = {
                "source_index": tetrahedron_index[source],
                "target_index": tetrahedron_index[target],
                "corrupted_rows": corrupted["rows"],
            }

    return {
        "rows": rows,
        "local_geometry": local_geometry,
        "transition_control": transition_ok,
        "inverse_control": inverse_ok,
        "face_control": face_ok,
        "local_fixed_ranks": dict(local_fixed_ranks),
        "local_augmented_ranks": dict(local_augmented_ranks),
        "first_face_data": first_face_data,
    }


def rank_record(built):
    fixed_rows = [
        {column: value for column, value in row.items() if column < CELL_COLUMNS}
        for row in built["rows"]
    ]
    fixed = {}
    augmented = {}
    widths = {}
    for prime in PRIMES:
        fixed_rank, fixed_width = modular_rank(
            fixed_rows, CELL_COLUMNS, prime
        )
        augmented_rank, augmented_width = modular_rank(
            built["rows"], TOTAL_COLUMNS, prime
        )
        fixed[str(prime)] = fixed_rank
        augmented[str(prime)] = augmented_rank
        widths[str(prime)] = {
            "fixed": fixed_width,
            "augmented": augmented_width,
        }
    return {
        "fixed_ranks": fixed,
        "augmented_ranks": augmented,
        "modular_nullities": {
            key: TOTAL_COLUMNS - value for key, value in augmented.items()
        },
        "maximum_elimination_width": widths,
    }


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "global_protocol": GLOBAL_PROTOCOL,
    "global_source": GLOBAL_SOURCE,
    "global_json": GLOBAL_JSON,
    "local_result": LOCAL_RESULT,
    "local_source": LOCAL_SOURCE,
    "local_json": LOCAL_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all canonical-data inputs have frozen provenance", provenance_ok, str(hashes))

global_artifact = json.loads(GLOBAL_JSON.read_text())
local_artifact = json.loads(LOCAL_JSON.read_text())
upstream_records = {
    (record["scale"], record["lapse"]): record["complete_modular_ranks"]
    for record in global_artifact["records"]
}
upstream_ok = bool(
    global_artifact["outcome"] == "GLOBAL_VARIABLE_FACE_MODULAR_DEFICIT_OPEN"
    and global_artifact["passed"] == global_artifact["tests"] == 12
    and all(
        upstream_records[representative]
        == {str(prime): CELL_COLUMNS for prime in PRIMES}
        for representative in REPRESENTATIVES
    )
    and local_artifact["outcome"] == "ONE_CONNECTION_COUPLED_RELATIVE_MODE"
    and local_artifact["passed"] == local_artifact["tests"] == 11
)
check("the local and complete variable-face controls persist", upstream_ok)

complex_data = build_exact_complex()
tetrahedra = complex_data["tetrahedra"]
dual = nx.Graph()
dual.add_nodes_from(range(len(tetrahedra)))
tetrahedron_index = {
    tetrahedron: index for index, tetrahedron in enumerate(tetrahedra)
}
for face, incident in complex_data["face_to_tetrahedra"].items():
    left, right = incident
    dual.add_edge(tetrahedron_index[left], tetrahedron_index[right])
f_vector = (
    len(complex_data["vertices"]),
    len(complex_data["edges"]),
    len(complex_data["triangles"]),
    len(tetrahedra),
)
incidence_ok = bool(
    f_vector == (120, 720, 1200, 600)
    and Counter(map(len, complex_data["face_to_tetrahedra"].values()))
    == Counter({2: 1200})
    and Counter(map(len, complex_data["edge_to_tetrahedra"].values()))
    == Counter({5: 720})
    and nx.is_connected(dual)
    and set(dict(dual.degree()).values()) == {4}
)
check("the exact 600-cell incidence is reconstructed", incidence_ok, f"f={f_vector}")

records = []
baseline_builds = {}
baseline_controls = True
fixed_rank_ok = True
homothetic_ok = True
corruption_ok = True
for scale, lapse in REPRESENTATIVES:
    print(f"building baseline (lambda,tau)=({scale},{lapse})", flush=True)
    built = build_global(
        complex_data, scale, lapse, ETA, CANONICAL_BASE
    )
    baseline_builds[(scale, lapse)] = built
    ranks = rank_record(built)
    vectors, data_rank = homothetic_vectors(
        built["local_geometry"],
        built["local_geometry"]["right_inverse"],
        len(tetrahedra),
    )
    residuals = {
        name: max(
            (abs(row_dot(row, vector)) for row in built["rows"]),
            default=0,
        )
        for name, vector in vectors
    }
    first_face = built["first_face_data"]
    corrupted_residuals = {
        name: max(
            (
                abs(row_dot(row, vector))
                for row in first_face["corrupted_rows"]
            ),
            default=0,
        )
        for name, vector in vectors
    }
    baseline_controls &= bool(
        built["local_geometry"]["controls"]
        and built["transition_control"]
        and built["inverse_control"]
        and built["face_control"]
        and built["local_fixed_ranks"] == {5: 1200}
    )
    fixed_rank_ok &= all(
        value == CELL_COLUMNS for value in ranks["fixed_ranks"].values()
    )
    homothetic_ok &= bool(
        data_rank == 2 and all(value == 0 for value in residuals.values())
    )
    corruption_ok &= any(value != 0 for value in corrupted_residuals.values())
    records.append({
        "scale": scale,
        "lapse": lapse,
        "local_fixed_rank_census": built["local_fixed_ranks"],
        "local_augmented_rank_census": built["local_augmented_ranks"],
        "homothetic_data_rank": data_rank,
        "homothetic_residuals": {key: str(value) for key, value in residuals.items()},
        "corrupted_first_face_residuals": {
            key: str(value) for key, value in corrupted_residuals.items()
        },
        **ranks,
    })

check(
    "all local Jacobians, Poincare kernels and forced face blocks pass",
    baseline_controls,
    str([(r["scale"], r["local_augmented_rank_census"]) for r in records]),
)
check("the old nonstatic fixed-data rank 3600 is reproduced", fixed_rank_ok)
check("two independent rational homothetic tangents survive exactly", homothetic_ok)
check("the private-edge corruption destroys a homothetic face control", corruption_ok)

alternate_records = []
alternate_ok = True
for scale, lapse in REPRESENTATIVES:
    print(f"building alternate graph (lambda,tau)=({scale},{lapse})", flush=True)
    built = build_global(
        complex_data,
        scale,
        lapse,
        ETA,
        CANONICAL_BASE,
        use_alternate=True,
    )
    ranks = rank_record(built)
    baseline = next(
        record for record in records
        if (record["scale"], record["lapse"]) == (scale, lapse)
    )
    vectors, data_rank = homothetic_vectors(
        built["local_geometry"],
        built["local_geometry"]["alternate_right_inverse"],
        len(tetrahedra),
    )
    residuals = {
        name: max(
            (abs(row_dot(row, vector)) for row in built["rows"]),
            default=0,
        )
        for name, vector in vectors
    }
    record_ok = bool(
        ranks["fixed_ranks"] == baseline["fixed_ranks"]
        and ranks["augmented_ranks"] == baseline["augmented_ranks"]
        and data_rank == 2
        and all(value == 0 for value in residuals.values())
    )
    alternate_ok &= record_ok
    alternate_records.append({
        "scale": scale,
        "lapse": lapse,
        "fixed_ranks": ranks["fixed_ranks"],
        "augmented_ranks": ranks["augmented_ranks"],
        "homothetic_residuals": {key: str(value) for key, value in residuals.items()},
    })
check("the augmented rank is independent of the local right-inverse graph", alternate_ok)


def convention_record(name, metric, canonical, reverse):
    print(f"building convention control {name}", flush=True)
    built = build_global(
        complex_data,
        2,
        5,
        metric,
        canonical,
        reverse_orientation=reverse,
    )
    ranks = rank_record(built)
    return {
        "name": name,
        "fixed_ranks": ranks["fixed_ranks"],
        "augmented_ranks": ranks["augmented_ranks"],
        "local_controls": bool(
            built["local_geometry"]["controls"]
            and built["transition_control"]
            and built["inverse_control"]
            and built["face_control"]
        ),
    }


odd_canonical = (
    CANONICAL_BASE[1],
    CANONICAL_BASE[0],
    CANONICAL_BASE[2],
    CANONICAL_BASE[3],
)
convention_records = [
    convention_record("reverse_faces", ETA, CANONICAL_BASE, True),
    convention_record("odd_relabelling", ETA, odd_canonical, False),
    convention_record("metric_sign", -ETA, CANONICAL_BASE, False),
]
baseline_first = records[0]
conventions_ok = all(
    record["local_controls"]
    and record["fixed_ranks"] == baseline_first["fixed_ranks"]
    and record["augmented_ranks"] == baseline_first["augmented_ranks"]
    for record in convention_records
)
check("orientation, relabelling and metric-sign attacks preserve all ranks", conventions_ok)

controls_ok = bool(
    provenance_ok
    and upstream_ok
    and incidence_ok
    and baseline_controls
    and fixed_rank_ok
    and homothetic_ok
    and corruption_ok
    and alternate_ok
    and conventions_ok
)
nullity_pairs = [
    tuple(record["modular_nullities"][str(prime)] for prime in PRIMES)
    for record in records
]
primes_agree = all(left == right for left, right in nullity_pairs)
only_homothetic = controls_ok and primes_agree and all(
    left == 2 for left, _ in nullity_pairs
)
intermediate = controls_ok and primes_agree and all(
    2 < left < DATA_COLUMNS for left, _ in nullity_pairs
)
full_modular = controls_ok and primes_agree and all(
    left == DATA_COLUMNS for left, _ in nullity_pairs
)

if not controls_ok:
    outcome = "CANONICAL_DATA_ADMISSIBILITY_CONTROL_FAILED"
elif only_homothetic:
    outcome = "CANONICAL_DATA_ONLY_HOMOTHETIC"
elif intermediate:
    outcome = "CANONICAL_DATA_INTERMEDIATE_MODULAR_OPEN"
elif full_modular:
    outcome = "CANONICAL_DATA_FULL_MODULAR_OPEN"
elif not primes_agree:
    outcome = "CANONICAL_DATA_PRIME_DISAGREEMENT_OPEN"
else:
    outcome = "CANONICAL_DATA_ADMISSIBILITY_OPEN"

allowed_outcomes = {
    "CANONICAL_DATA_ADMISSIBILITY_CONTROL_FAILED",
    "CANONICAL_DATA_ONLY_HOMOTHETIC",
    "CANONICAL_DATA_INTERMEDIATE_MODULAR_OPEN",
    "CANONICAL_DATA_FULL_MODULAR_OPEN",
    "CANONICAL_DATA_PRIME_DISAGREEMENT_OPEN",
    "CANONICAL_DATA_ADMISSIBILITY_OPEN",
}
check("the preregistered admissibility hierarchy assigns one outcome", outcome in allowed_outcomes, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "f_vector": list(f_vector),
    "representatives": [list(value) for value in REPRESENTATIVES],
    "primes": list(PRIMES),
    "column_counts": {
        "cell_flex": CELL_COLUMNS,
        "upper_edges": EDGE_COUNT,
        "struts": VERTEX_COUNT,
        "augmented": TOTAL_COLUMNS,
    },
    "records": records,
    "alternate_right_inverse_records": alternate_records,
    "convention_records": convention_records,
    "classification": {
        "rational_lower_bound_from_explicit_tangents": 2,
        "rational_exact_dimension": 2 if only_homothetic else "OPEN",
        "action_hessian": "NOT EVALUATED",
        "physical_tick_or_speed": "NOT EVALUATED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
for record in records:
    print(
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"fixed={record['fixed_ranks']}, "
        f"augmented={record['augmented_ranks']}, "
        f"nullities={record['modular_nullities']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

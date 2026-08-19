#!/usr/bin/env python3
"""Exact complete variable-face closure on the regular 600-cell."""

from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import networkx as nx
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_global_variable_face_closure.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_global_variable_face_closure_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_global_variable_face_closure_protocol.md"
LOCAL_RESULT = ROOT / "docs/gravity/gravity_600cell_variable_face_connection_result.md"
LOCAL_SOURCE = HERE / "verify_gravity_600cell_variable_face_connection.py"
LOCAL_JSON = HERE / "gravity_600cell_variable_face_connection.json"
LOCAL_ADV_SOURCE = HERE / "verify_gravity_600cell_variable_face_connection_adversarial.py"
LOCAL_ADV_JSON = HERE / "gravity_600cell_variable_face_connection_adversarial.json"
FROZEN_SOURCE = HERE / "verify_gravity_600cell_global_flex_holonomy_adversarial.py"
FROZEN_JSON = HERE / "gravity_600cell_global_flex_holonomy_adversarial.json"

PROTOCOL_COMMIT = "77e5370"
EXPECTED_HASHES = {
    "prior_art": "a76a28e8247e2fd1d0ea3536e6a345ba4a091ce931bc0e9c570f3234286b5014",
    "protocol": "ed79c6a15ade377ae09854b3cad3028eb1c0f43cf8e85789d46993fe25ed1b49",
    "local_result": "2db55cb87ec1c01d537cdbc11010bc9ea740762c598108e4c2de0f3acca72cc8",
    "local_source": "69a5d7479a5df427cead76f82db31fe62a9190c28c967f699c846881634fb0f6",
    "local_json": "001212016553d006862e68edc4f780f37ca1476110b6e0aed3e987f52a43b5e3",
    "local_adv_source": "9a3c6985eb4d833ef4ecc21f9e964577d102d5d2e7beac6a4163c4225faa5984",
    "local_adv_json": "c8c8c58711e5bf4e49c110e84518ddf643b75cc4377d05fb5f577003b8395466",
    "frozen_source": "54fa9775a2f14d708359167d3f8b81e03d985f24594b453f16028d9981d9be0d",
    "frozen_json": "f224fe123c882ccda97d4ca6ec67c9fd810d58ed8377c5afb457a1dec69f4b87",
}

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
CANONICAL_BASE = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
))
REPRESENTATIVES = ((1, 5), (2, 5), (3, 11))
PRIMES = (1000003, 1000033)

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


# Exact Q(phi), represented by a+b*phi and phi^2=phi+1.
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
    return sum(permutation[i] > permutation[j]
               for i in range(len(permutation))
               for j in range(i + 1, len(permutation))) % 2


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
        nonzero = tuple(index for index, value in enumerate(permuted)
                        if value != zero)
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
            face_to_tetrahedra.setdefault(tuple(sorted(face)), []).append(tetrahedron)
        for edge in combinations(tetrahedron, 2):
            edge_to_tetrahedra.setdefault(tuple(sorted(edge)), []).append(tetrahedron)
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


def local_kernel(scale, lapse, basis, normal):
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
    return {vertex: canonical[index]
            for index, vertex in enumerate(tetrahedron)}


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
        - 2 * normal * (normal.T * metric * (point - anchor))[0]
        / denominator
    )


def affine_map(domain, codomain):
    domain_h = sp.Matrix.hstack(*(
        point.col_join(sp.ones(1, 1)) for point in domain
    ))
    codomain_h = sp.Matrix.hstack(*(
        point.col_join(sp.ones(1, 1)) for point in codomain
    ))
    return sp.simplify(codomain_h * domain_h.inv())


def lateral_transition(source, target, shared, canonical, scale, lapse, metric):
    source_coordinates = coordinate_map(source, canonical)
    target_coordinates = coordinate_map(target, canonical)
    source_top = {vertex: scale * point + lapse * NORMAL
                  for vertex, point in source_coordinates.items()}
    target_top = {vertex: scale * point + lapse * NORMAL
                  for vertex, point in target_coordinates.items()}
    source_apex = next(vertex for vertex in source if vertex not in shared)
    target_apex = next(vertex for vertex in target if vertex not in shared)

    shared_face_source = tuple(source_coordinates[vertex] for vertex in shared)
    lateral_points_source = shared_face_source + (source_top[shared[0]],)
    reflected_lower = affine_reflection(
        source_coordinates[source_apex], lateral_points_source, metric
    )
    reflected_upper = affine_reflection(
        source_top[source_apex], lateral_points_source, metric
    )

    domain = (
        tuple(target_coordinates[vertex] for vertex in shared)
        + (target_top[shared[0]], target_coordinates[target_apex])
    )
    codomain = lateral_points_source + (reflected_lower,)
    transition = affine_map(domain, codomain)
    linear = transition[:4, :4]

    shared_mapping = all(
        sp.simplify(
            transition * target_coordinates[vertex].col_join(sp.ones(1, 1))
            - source_coordinates[vertex].col_join(sp.ones(1, 1))
        ) == sp.zeros(5, 1)
        and sp.simplify(
            transition * target_top[vertex].col_join(sp.ones(1, 1))
            - source_top[vertex].col_join(sp.ones(1, 1))
        ) == sp.zeros(5, 1)
        for vertex in shared
    )
    apex_mapping = bool(
        sp.simplify(
            transition * target_top[target_apex].col_join(sp.ones(1, 1))
            - reflected_upper.col_join(sp.ones(1, 1))
        ) == sp.zeros(5, 1)
    )
    control = bool(
        transition[4, :] == sp.Matrix([[0, 0, 0, 0, 1]])
        and sp.simplify(linear.T * metric * linear) == metric
        and abs(linear.det()) == 1
        and shared_mapping and apex_mapping
    )
    return {
        "transition": transition,
        "linear": linear,
        "control": control,
        "time_space_mixing": any(
            linear[index, 3] != 0 or linear[3, index] != 0
            for index in range(3)
        ),
        "source_coordinates": source_coordinates,
        "target_coordinates": target_coordinates,
        "source_top": source_top,
        "target_top": target_top,
    }


def face_reduced_block(source, target, shared, transition_data,
                       scale, lapse, metric, basis):
    local = local_kernel(scale, lapse, basis, NORMAL)
    source_upper = tuple(transition_data["source_top"][vertex]
                         for vertex in shared)
    target_upper = tuple(transition_data["target_top"][vertex]
                         for vertex in shared)
    source_lower = tuple(transition_data["source_coordinates"][vertex]
                         for vertex in shared)
    source_evaluation = poincare_evaluation(source_upper, basis)
    target_evaluation = poincare_evaluation(target_upper, basis)
    lower_evaluation = poincare_evaluation(source_lower, basis)
    stabilizer = kernel(lower_evaluation)
    linear_block = sp.diag(*([transition_data["linear"]] * 3))
    cell_pair = (source_evaluation * local).row_join(
        -linear_block * target_evaluation * local
    )
    raw = cell_pair.row_join(-source_evaluation * stabilizer)
    raw_kernel = kernel(raw)
    projected = raw_kernel[:12, :]
    annihilator = kernel(projected.T)
    reduced = annihilator.T
    fixed_kernel = kernel(cell_pair)
    control = bool(
        stabilizer.shape == (10, 1)
        and raw.shape == (12, 13)
        and raw.rank() == 6
        and raw_kernel.shape == (13, 7)
        and raw_kernel[12:13, :].rank() == 1
        and fixed_kernel.shape == (12, 6)
        and projected.rank() == 7
        and reduced.shape == (5, 12)
        and reduced.rank() == 5
        and same_space(kernel(reduced), projected)
    )
    return reduced, control


def matrix_rows(blocks, tetrahedron_index, face_order):
    rows = []
    for face_index in face_order:
        source, target, block = blocks[face_index]
        source_offset = 6 * tetrahedron_index[source]
        target_offset = 6 * tetrahedron_index[target]
        for row_index in range(5):
            row = {}
            for local_column in range(6):
                left_value = sp.cancel(block[row_index, local_column])
                right_value = sp.cancel(block[row_index, 6 + local_column])
                if left_value:
                    row[source_offset + local_column] = left_value
                if right_value:
                    row[target_offset + local_column] = right_value
            rows.append(row)
    return rows


def generic_control_rows(face_pairs, tetrahedron_count):
    rows = []
    for number, (source_index, target_index) in enumerate(face_pairs, start=1):
        source_offset = 6 * source_index
        target_offset = 6 * target_index
        for coordinate in range(5):
            row = {
                source_offset + coordinate: -number,
                source_offset + coordinate + 1: 1,
                target_offset + coordinate: number,
                target_offset + coordinate + 1: -1,
            }
            rows.append(row)
    return rows


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
        row = {
            column: rational_mod(value, prime)
            for column, value in original.items()
            if rational_mod(value, prime)
        }
        while row:
            existing = [column for column in row if column in pivots]
            if not existing:
                pivot = min(row)
                inverse = pow(row[pivot], -1, prime)
                row = {column: value * inverse % prime
                       for column, value in row.items() if value % prime}
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


def build_blocks(complex_data, scale, lapse, metric, canonical,
                 reverse_orientation=False):
    tetrahedra = complex_data["tetrahedra"]
    tetrahedron_index = {tetrahedron: index
                         for index, tetrahedron in enumerate(tetrahedra)}
    faces = tuple(sorted(complex_data["face_to_tetrahedra"]))
    blocks = {}
    all_transitions = True
    all_inverses = True
    all_local_blocks = True
    any_mixing = False
    transition_cache = {}
    basis = lorentz_basis(metric)

    for face_index, face in enumerate(faces):
        left, right = sorted(complex_data["face_to_tetrahedra"][face])
        source, target = ((right, left) if reverse_orientation else (left, right))
        shared = tuple(sorted(face))
        forward = lateral_transition(
            source, target, shared, canonical, scale, lapse, metric
        )
        backward = lateral_transition(
            target, source, shared, canonical, scale, lapse, metric
        )
        transition_cache[(source, target)] = forward["transition"]
        transition_cache[(target, source)] = backward["transition"]
        all_transitions &= forward["control"] and backward["control"]
        all_inverses &= bool(
            sp.simplify(forward["transition"] * backward["transition"])
            == sp.eye(5)
        )
        any_mixing |= forward["time_space_mixing"]
        block, local_control = face_reduced_block(
            source, target, shared, forward, scale, lapse, metric, basis
        )
        all_local_blocks &= local_control
        blocks[face_index] = (source, target, block)
    return {
        "blocks": blocks,
        "faces": faces,
        "tetrahedron_index": tetrahedron_index,
        "transition_control": all_transitions,
        "inverse_control": all_inverses,
        "local_block_control": all_local_blocks,
        "any_time_space_mixing": any_mixing,
    }


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "local_result": LOCAL_RESULT,
    "local_source": LOCAL_SOURCE,
    "local_json": LOCAL_JSON,
    "local_adv_source": LOCAL_ADV_SOURCE,
    "local_adv_json": LOCAL_ADV_JSON,
    "frozen_source": FROZEN_SOURCE,
    "frozen_json": FROZEN_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all complete variable-face inputs have frozen provenance",
      provenance_ok, str(hashes))

local_primary = json.loads(LOCAL_JSON.read_text())
local_adversarial = json.loads(LOCAL_ADV_JSON.read_text())
frozen = json.loads(FROZEN_JSON.read_text())
upstream_ok = bool(
    local_primary["outcome"] == "ONE_CONNECTION_COUPLED_RELATIVE_MODE"
    and local_primary["passed"] == local_primary["tests"] == 11
    and local_adversarial["outcome"] == "ADVERSARIAL_ONE_CONNECTION_MODE"
    and local_adversarial["passed"] == local_adversarial["tests"] == 9
    and frozen["outcome"] == "ADVERSARIAL_GLOBAL_FLEX_SEED_KILLED"
    and frozen["passed"] == frozen["tests"] == 11
)
check("all local correction and frozen-loop controls persist", upstream_ok)

complex_data = build_exact_complex()
tetrahedra = complex_data["tetrahedra"]
tetrahedron_index = {tetrahedron: index
                     for index, tetrahedron in enumerate(tetrahedra)}
faces = tuple(sorted(complex_data["face_to_tetrahedra"]))
dual = nx.Graph()
dual.add_nodes_from(range(len(tetrahedra)))
face_pairs = []
pair_to_face = {}
for face_index, face in enumerate(faces):
    left, right = sorted(complex_data["face_to_tetrahedra"][face])
    pair = (tetrahedron_index[left], tetrahedron_index[right])
    dual.add_edge(*pair, face_index=face_index)
    face_pairs.append(pair)
    pair_to_face[frozenset(pair)] = face_index

f_vector = (
    len(complex_data["vertices"]), len(complex_data["edges"]),
    len(complex_data["triangles"]), len(tetrahedra),
)
face_incidence = Counter(map(len, complex_data["face_to_tetrahedra"].values()))
edge_incidence = Counter(map(len, complex_data["edge_to_tetrahedra"].values()))
incidence_ok = bool(
    f_vector == (120, 720, 1200, 600)
    and face_incidence == Counter({2: 1200})
    and edge_incidence == Counter({5: 720})
    and nx.is_connected(dual)
    and set(dict(dual.degree()).values()) == {4}
    and nx.edge_connectivity(dual) == 4
)
check("the exact carrier has the required four-edge-connected dual",
      incidence_ok,
      f"f={f_vector}, face={dict(face_incidence)}, edge={dict(edge_incidence)}")

bfs_pairs = list(nx.bfs_edges(dual, source=0))
tree_face_order = [pair_to_face[frozenset(pair)] for pair in bfs_pairs]
all_face_order = list(range(len(faces)))

records = []
transition_controls = True
local_block_controls = True
tree_controls = True
physical_ranks = {}
baseline_builds = {}

for scale, lapse in REPRESENTATIVES:
    built = build_blocks(
        complex_data, scale, lapse, ETA, CANONICAL_BASE,
        reverse_orientation=False,
    )
    baseline_builds[(scale, lapse)] = built
    transition_ok = bool(
        built["transition_control"] and built["inverse_control"]
        and (scale == 1 or built["any_time_space_mixing"])
    )
    transition_controls &= transition_ok
    local_block_controls &= built["local_block_control"]
    rows = matrix_rows(built["blocks"], tetrahedron_index, all_face_order)
    tree_rows = matrix_rows(
        built["blocks"], tetrahedron_index, tree_face_order
    )
    prime_ranks = {}
    tree_ranks = {}
    widths = {}
    for prime in PRIMES:
        rank, width = modular_rank(rows, 3600, prime)
        tree_rank, tree_width = modular_rank(tree_rows, 3600, prime)
        prime_ranks[str(prime)] = rank
        tree_ranks[str(prime)] = tree_rank
        widths[str(prime)] = {
            "complete": width,
            "tree": tree_width,
        }
    physical_ranks[(scale, lapse)] = prime_ranks
    tree_ok = all(rank == 2995 for rank in tree_ranks.values())
    tree_controls &= tree_ok
    records.append({
        "scale": scale,
        "lapse": lapse,
        "time_space_mixing": built["any_time_space_mixing"],
        "complete_modular_ranks": prime_ranks,
        "tree_modular_ranks": tree_ranks,
        "maximum_elimination_width": widths,
    })

check("all 2400 full lateral transitions are inverse Lorentz isometries",
      transition_controls)
check("all 1200 isolated face blocks reproduce the 6-to-7 theorem",
      local_block_controls)
check("every physical spanning-tree restriction has rank/nullity 2995/605",
      tree_controls, str(records))

generic_rows = generic_control_rows(face_pairs, len(tetrahedra))
generic_ranks = {}
generic_widths = {}
for prime in PRIMES:
    rank, width = modular_rank(generic_rows, 3600, prime)
    generic_ranks[str(prime)] = rank
    generic_widths[str(prime)] = width
generic_ok = all(rank == 3594 for rank in generic_ranks.values())
check("the deterministic generic body-hinge control has nullity six",
      generic_ok, str(generic_ranks))

check("all complete physical modular ranks are recorded without fitting",
      all(set(ranks) == {str(value) for value in PRIMES}
          for ranks in physical_ranks.values()), str(physical_ranks))

# Convention attacks are rebuilt at the static representative.  Exact
# per-face controls plus equal complete ranks are the declared criterion.
static_rows = matrix_rows(
    baseline_builds[(1, 5)]["blocks"], tetrahedron_index, all_face_order
)
static_baseline = physical_ranks[(1, 5)]

reversed_build = build_blocks(
    complex_data, 1, 5, ETA, CANONICAL_BASE, reverse_orientation=True
)
reversed_rows = matrix_rows(
    reversed_build["blocks"], tetrahedron_index, all_face_order
)
reversed_ranks = {
    str(prime): modular_rank(reversed_rows, 3600, prime)[0]
    for prime in PRIMES
}
reverse_ok = bool(
    reversed_build["transition_control"]
    and reversed_build["inverse_control"]
    and reversed_build["local_block_control"]
    and reversed_ranks == static_baseline
)
check("reversing every face orientation preserves the complete rank",
      reverse_ok, str(reversed_ranks))

odd_canonical = (
    CANONICAL_BASE[1], CANONICAL_BASE[0],
    CANONICAL_BASE[2], CANONICAL_BASE[3],
)
odd_build = build_blocks(
    complex_data, 1, 5, ETA, odd_canonical, reverse_orientation=False
)
odd_rows = matrix_rows(odd_build["blocks"], tetrahedron_index, all_face_order)
odd_ranks = {
    str(prime): modular_rank(odd_rows, 3600, prime)[0]
    for prime in PRIMES
}
odd_ok = bool(
    odd_build["transition_control"] and odd_build["inverse_control"]
    and odd_build["local_block_control"]
    and odd_ranks == static_baseline
)
check("an odd canonical relabelling preserves the complete rank",
      odd_ok, str(odd_ranks))

sign_build = build_blocks(
    complex_data, 1, 5, -ETA, CANONICAL_BASE, reverse_orientation=False
)
sign_rows = matrix_rows(sign_build["blocks"], tetrahedron_index, all_face_order)
sign_ranks = {
    str(prime): modular_rank(sign_rows, 3600, prime)[0]
    for prime in PRIMES
}
sign_ok = bool(
    sign_build["transition_control"] and sign_build["inverse_control"]
    and sign_build["local_block_control"]
    and sign_ranks == static_baseline
)
check("metric-sign reversal preserves the complete static rank",
      sign_ok, str(sign_ranks))

controls_ok = bool(
    provenance_ok and upstream_ok and incidence_ok and transition_controls
    and local_block_controls and tree_controls and generic_ok
    and reverse_ok and odd_ok and sign_ok
)
full_rank = bool(
    controls_ok
    and all(rank == 3600
            for ranks in physical_ranks.values() for rank in ranks.values())
)
modular_deficit = bool(
    controls_ok
    and any(rank < 3600
            for ranks in physical_ranks.values() for rank in ranks.values())
)

if not controls_ok:
    outcome = "GLOBAL_VARIABLE_FACE_CONTROL_FAILED"
elif full_rank:
    outcome = "GLOBAL_VARIABLE_FACE_KERNEL_ZERO"
elif modular_deficit:
    outcome = "GLOBAL_VARIABLE_FACE_MODULAR_DEFICIT_OPEN"
else:
    outcome = "GLOBAL_VARIABLE_FACE_OPEN"

allowed = {
    "GLOBAL_VARIABLE_FACE_CONTROL_FAILED",
    "GLOBAL_VARIABLE_FACE_KERNEL_ZERO",
    "GLOBAL_VARIABLE_FACE_POSITIVE_KERNEL",
    "GLOBAL_VARIABLE_FACE_MODULAR_DEFICIT_OPEN",
    "GLOBAL_VARIABLE_FACE_OPEN",
}
check("the complete variable-face hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "f_vector": list(f_vector),
    "dual_edge_connectivity": nx.edge_connectivity(dual),
    "primes": list(PRIMES),
    "records": records,
    "generic_control_ranks": generic_ranks,
    "generic_control_widths": generic_widths,
    "reverse_static_ranks": reversed_ranks,
    "odd_relabelling_static_ranks": odd_ranks,
    "metric_sign_static_ranks": sign_ranks,
    "classification": {
        "complete_variable_face_kernel": (
            "ZERO OVER Q" if full_rank else "OPEN"
        ),
        "frozen_connection_correction": "INCLUDED",
        "finite_reconstruction": "NOT TESTED",
        "action_hessian_or_dynamics": "NOT TESTED",
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
        f"ranks={record['complete_modular_ranks']}, "
        f"tree={record['tree_modular_ranks']}"
    )
print("generic:", generic_ranks)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

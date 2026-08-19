#!/usr/bin/env python3
"""Independent static rotation/tangential-translation decomposition audit."""

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
OUTPUT = HERE / "gravity_600cell_static_vertex_gradient_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_static_vertex_gradient_adversarial_protocol.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_static_vertex_gradient_prior_art.md"
PRIMARY_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_static_vertex_gradient_protocol.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_static_vertex_gradient.py"
PRIMARY_JSON = HERE / "gravity_600cell_static_vertex_gradient.json"
GLOBAL_JSON = HERE / "gravity_600cell_global_variable_face_closure.json"

PROTOCOL_COMMIT = "07906b9"
EXPECTED_HASHES = {
    "protocol": "53c094c120411e0b63cafd9b8b5e4b60880c7f1339eabec25eed7f75ccdf2805",
    "prior_art": "9541383a435e069be13ed9c2175674036a9cef5cd4e17ee455f524bd1c1c6a7d",
    "primary_protocol": "c689202b94abfe9436bdc8f5db6f79fa00d6a7bc743d79bc3ae34d391890c17b",
    "primary_source": "6974ef8c85ee62b32daa2277ba221ee1cfa96f1c7cc7a92a8ec91fad576b124f",
    "primary_json": "ce018db5c66c78e89e4ca32360385955ea520b9ac8e42955b110d190432239c0",
    "global_json": "61cebd1cd67fcdc56de088855b1fc7b805d0f70f9f9b3029d4a61209d7a53944",
}

CANONICAL = tuple(sp.Matrix(point) for point in (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
))
ROTATIONS = tuple(
    matrix for matrix in (
        sp.Matrix(((0, 1, 0), (-1, 0, 0), (0, 0, 0))),
        sp.Matrix(((0, 0, 1), (0, 0, 0), (-1, 0, 0))),
        sp.Matrix(((0, 0, 0), (0, 0, 1), (0, -1, 0))),
    )
)
ROTATION_COORDINATES = sp.Matrix.hstack(*(
    value.reshape(9, 1) for value in ROTATIONS
))
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


def build_complex():
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
    for tetrahedron in tetrahedra:
        for face in combinations(tetrahedron, 3):
            face_to_tetrahedra.setdefault(tuple(sorted(face)), []).append(tetrahedron)
    return {
        "vertices": vertices,
        "edges": tuple(edges),
        "triangles": tuple(triangles),
        "tetrahedra": tetrahedra,
        "face_to_tetrahedra": face_to_tetrahedra,
    }


def coordinates(tetrahedron, canonical):
    return {vertex: canonical[index]
            for index, vertex in enumerate(tetrahedron)}


def reflect(point, face_points):
    anchor, second, third = face_points
    normal = (second - anchor).cross(third - anchor)
    return sp.simplify(
        point - 2 * normal * (normal.T * (point - anchor))[0]
        / (normal.T * normal)[0]
    )


def affine_map(domain, codomain):
    domain_h = sp.Matrix.hstack(*(
        point.col_join(sp.ones(1, 1)) for point in domain
    ))
    codomain_h = sp.Matrix.hstack(*(
        point.col_join(sp.ones(1, 1)) for point in codomain
    ))
    return sp.simplify(codomain_h * domain_h.inv())


def spatial_transition(source, target, shared, canonical):
    source_coordinates = coordinates(source, canonical)
    target_coordinates = coordinates(target, canonical)
    source_apex = next(vertex for vertex in source if vertex not in shared)
    target_apex = next(vertex for vertex in target if vertex not in shared)
    reflected = reflect(
        source_coordinates[source_apex],
        tuple(source_coordinates[vertex] for vertex in shared),
    )
    domain = tuple(target_coordinates[vertex] for vertex in target)
    codomain = tuple(
        source_coordinates[vertex] if vertex in shared else reflected
        for vertex in target
    )
    affine = affine_map(domain, codomain)
    linear = affine[:3, :3]
    control = bool(
        affine[3, :] == sp.Matrix([[0, 0, 0, 1]])
        and linear.T * linear == sp.eye(3)
        and abs(linear.det()) == 1
        and all(
            affine * target_coordinates[vertex].col_join(sp.ones(1, 1))
            == source_coordinates[vertex].col_join(sp.ones(1, 1))
            for vertex in shared
        )
    )
    return affine, linear, control, source_coordinates, target_coordinates


def rotation_adjoint(linear):
    result = sp.zeros(3)
    for column, generator in enumerate(ROTATIONS):
        transformed = sp.simplify(linear * generator * linear.inv())
        coordinates_vector, free = ROTATION_COORDINATES.gauss_jordan_solve(
            transformed.reshape(9, 1)
        )
        if free.rows:
            raise RuntimeError("ambiguous rotation coordinates")
        result[:, column] = coordinates_vector
    return result


def local_gradient(canonical):
    directions = sp.Matrix.hstack(*(
        canonical[index] - canonical[0] for index in range(1, 4)
    ))
    difference = sp.Matrix((
        (-1, 1, 0, 0),
        (-1, 0, 1, 0),
        (-1, 0, 0, 1),
    ))
    return sp.simplify(directions.T.inv() * difference)


def cell_gradient_map(tetrahedron, union, gradient):
    selection = sp.zeros(4, len(union))
    locations = {vertex: index for index, vertex in enumerate(union)}
    for local_index, vertex in enumerate(tetrahedron):
        selection[local_index, locations[vertex]] = 1
    return gradient * selection


def build_operators(complex_data, canonical):
    tetrahedra = complex_data["tetrahedra"]
    tetrahedron_index = {tetrahedron: index
                         for index, tetrahedron in enumerate(tetrahedra)}
    faces = tuple(sorted(complex_data["face_to_tetrahedra"]))
    rotation_rows = []
    translation_rows = []
    flat_rotation_rows = []
    face_rotation_rows = {}
    face_translation_rows = {}
    transition_control = True
    inverse_control = True
    face_gradient_control = True
    discontinuity_detected = False
    gradient = local_gradient(canonical)

    for face_index, face in enumerate(faces):
        source, target = sorted(complex_data["face_to_tetrahedra"][face])
        shared = tuple(sorted(face))
        forward = spatial_transition(source, target, shared, canonical)
        backward = spatial_transition(target, source, shared, canonical)
        affine, linear, control, source_coordinates, _ = forward
        transition_control &= control and backward[2]
        inverse_control &= bool(affine * backward[0] == sp.eye(4))
        adjoint = rotation_adjoint(linear)
        source_offset = 3 * tetrahedron_index[source]
        target_offset = 3 * tetrahedron_index[target]

        local_rotation_rows = []
        local_flat_rows = []
        for coordinate in range(3):
            row = {source_offset + coordinate: sp.Integer(1)}
            flat = {
                source_offset + coordinate: sp.Integer(1),
                target_offset + coordinate: sp.Integer(-1),
            }
            for target_coordinate in range(3):
                value = -adjoint[coordinate, target_coordinate]
                if value:
                    row[target_offset + target_coordinate] = value
            local_rotation_rows.append(row)
            local_flat_rows.append(flat)
        rotation_rows.extend(local_rotation_rows)
        flat_rotation_rows.extend(local_flat_rows)
        face_rotation_rows[face_index] = local_rotation_rows

        source_face = tuple(source_coordinates[vertex] for vertex in shared)
        tangents = (source_face[1] - source_face[0],
                    source_face[2] - source_face[0])
        local_translation_rows = []
        for tangent in tangents:
            row = {}
            for coordinate in range(3):
                if tangent[coordinate]:
                    row[source_offset + coordinate] = tangent[coordinate]
                transported = -(tangent.T * linear)[coordinate]
                if transported:
                    row[target_offset + coordinate] = transported
            local_translation_rows.append(row)
        translation_rows.extend(local_translation_rows)
        face_translation_rows[face_index] = local_translation_rows

        union = tuple(sorted(set(source) | set(target)))
        source_gradient = cell_gradient_map(source, union, gradient)
        target_gradient = cell_gradient_map(target, union, gradient)
        local_matrix = sp.zeros(2, 6)
        for row_index, row in enumerate(local_translation_rows):
            for coordinate in range(3):
                local_matrix[row_index, coordinate] = row.get(
                    source_offset + coordinate, 0
                )
                local_matrix[row_index, 3 + coordinate] = row.get(
                    target_offset + coordinate, 0
                )
        pair_gradient = source_gradient.col_join(target_gradient)
        face_gradient_control &= bool(
            local_matrix * pair_gradient == sp.zeros(2, len(union))
        )
        if face_index == 0:
            wrong_target = target_gradient.copy()
            first, second = union.index(shared[0]), union.index(shared[1])
            wrong_target[:, first], wrong_target[:, second] = (
                wrong_target[:, second], wrong_target[:, first]
            )
            discontinuity_detected = bool(
                local_matrix * source_gradient.col_join(wrong_target)
                != sp.zeros(2, len(union))
            )

    gradient_rows = []
    for tetrahedron in tetrahedra:
        local = cell_gradient_map(tetrahedron, tuple(range(120)), gradient)
        for coordinate in range(3):
            gradient_rows.append({
                vertex: sp.cancel(local[coordinate, vertex])
                for vertex in tetrahedron if local[coordinate, vertex] != 0
            })
    return {
        "tetrahedron_index": tetrahedron_index,
        "faces": faces,
        "rotation_rows": rotation_rows,
        "translation_rows": translation_rows,
        "flat_rotation_rows": flat_rotation_rows,
        "face_rotation_rows": face_rotation_rows,
        "face_translation_rows": face_translation_rows,
        "gradient_rows": gradient_rows,
        "transition_control": transition_control,
        "inverse_control": inverse_control,
        "face_gradient_control": face_gradient_control,
        "discontinuity_detected": discontinuity_detected,
    }


def rational_mod(value, prime):
    value = sp.Rational(value)
    numerator = int(value.p) % prime
    denominator = int(value.q) % prime
    if denominator == 0:
        raise ZeroDivisionError
    return numerator * pow(denominator, -1, prime) % prime


def modular_rank(rows, columns, prime):
    pivots = {}
    for original in rows:
        row = {}
        for column, value in original.items():
            residue = rational_mod(value, prime)
            if residue:
                row[column] = residue
        while row:
            existing = [column for column in row if column in pivots]
            if not existing:
                pivot = min(row)
                inverse = pow(row[pivot], -1, prime)
                row = {column: value * inverse % prime
                       for column, value in row.items() if value % prime}
                pivots[pivot] = row
                break
            pivot = min(existing)
            factor = row[pivot]
            for column, value in pivots[pivot].items():
                updated = (row.get(column, 0) - factor * value) % prime
                if updated:
                    row[column] = updated
                elif column in row:
                    del row[column]
        if len(pivots) == columns:
            break
    return len(pivots)


paths = {
    "protocol": PROTOCOL,
    "prior_art": PRIOR_ART,
    "primary_protocol": PRIMARY_PROTOCOL,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "global_json": GLOBAL_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all adversarial static inputs have frozen provenance",
      provenance_ok, str(hashes))

primary = json.loads(PRIMARY_JSON.read_text())
global_artifact = json.loads(GLOBAL_JSON.read_text())
upstream_ok = bool(
    primary["outcome"] == "STATIC_KERNEL_EXACTLY_VERTEX_GRADIENTS"
    and primary["passed"] == primary["tests"] == 10
    and global_artifact["outcome"] == "GLOBAL_VARIABLE_FACE_MODULAR_DEFICIT_OPEN"
    and global_artifact["passed"] == global_artifact["tests"] == 12
)
check("the primary gradient theorem and target-blind census persist", upstream_ok)

complex_data = build_complex()
incidence_ok = bool(
    len(complex_data["vertices"]) == 120
    and len(complex_data["edges"]) == 720
    and len(complex_data["triangles"]) == 1200
    and len(complex_data["tetrahedra"]) == 600
    and Counter(map(len, complex_data["face_to_tetrahedra"].values()))
    == Counter({2: 1200})
)
check("the independent carrier has exact 600-cell incidence", incidence_ok)

operators = build_operators(complex_data, CANONICAL)
geometry_ok = bool(
    operators["transition_control"] and operators["inverse_control"]
    and operators["face_gradient_control"]
    and operators["discontinuity_detected"]
)
check("all spatial transitions and tangential-gradient identities are exact",
      geometry_ok)

rotation_ranks = {
    str(prime): modular_rank(operators["rotation_rows"], 1800, prime)
    for prime in PRIMES
}
translation_ranks = {
    str(prime): modular_rank(operators["translation_rows"], 1800, prime)
    for prime in PRIMES
}
gradient_ranks = {
    str(prime): modular_rank(operators["gradient_rows"], 120, prime)
    for prime in PRIMES
}
sector_ok = bool(
    set(rotation_ranks.values()) == {1800}
    and set(translation_ranks.values()) == {1681}
    and set(gradient_ranks.values()) == {119}
)
check("curved rotations are full rank and translations have exactly 119 gradients",
      sector_ok,
      f"rotation={rotation_ranks}, translation={translation_ranks}, "
      f"gradient={gradient_ranks}")

flat_ranks = {
    str(prime): modular_rank(operators["flat_rotation_rows"], 1800, prime)
    for prime in PRIMES
}
flat_ok = set(flat_ranks.values()) == {1797}
check("the flat rotation control retains exactly three constant rotations",
      flat_ok, str(flat_ranks))

tetrahedra = complex_data["tetrahedra"]
tetrahedron_index = operators["tetrahedron_index"]
dual = nx.Graph()
pair_to_face = {}
for face_index, face in enumerate(operators["faces"]):
    left, right = sorted(complex_data["face_to_tetrahedra"][face])
    pair = (tetrahedron_index[left], tetrahedron_index[right])
    dual.add_edge(*pair)
    pair_to_face[frozenset(pair)] = face_index
tree_faces = [pair_to_face[frozenset(pair)]
              for pair in nx.bfs_edges(dual, source=0)]
tree_rotation_rows = [row for face in tree_faces
                      for row in operators["face_rotation_rows"][face]]
tree_translation_rows = [row for face in tree_faces
                         for row in operators["face_translation_rows"][face]]
tree_rotation_ranks = {
    str(prime): modular_rank(tree_rotation_rows, 1800, prime)
    for prime in PRIMES
}
tree_translation_ranks = {
    str(prime): modular_rank(tree_translation_rows, 1800, prime)
    for prime in PRIMES
}
tree_ok = bool(
    set(tree_rotation_ranks.values()) == {1797}
    and set(tree_translation_ranks.values()) == {1198}
)
check("the spanning-tree sectors reproduce combined nullity 605", tree_ok,
      f"rotation={tree_rotation_ranks}, translation={tree_translation_ranks}")

odd = (CANONICAL[1], CANONICAL[0], CANONICAL[2], CANONICAL[3])
odd_operators = build_operators(complex_data, odd)
odd_rotation_ranks = {
    str(prime): modular_rank(odd_operators["rotation_rows"], 1800, prime)
    for prime in PRIMES
}
odd_translation_ranks = {
    str(prime): modular_rank(odd_operators["translation_rows"], 1800, prime)
    for prime in PRIMES
}
odd_gradient_ranks = {
    str(prime): modular_rank(odd_operators["gradient_rows"], 120, prime)
    for prime in PRIMES
}
odd_ok = bool(
    odd_operators["transition_control"] and odd_operators["inverse_control"]
    and odd_operators["face_gradient_control"]
    and odd_rotation_ranks == rotation_ranks
    and odd_translation_ranks == translation_ranks
    and odd_gradient_ranks == gradient_ranks
)
check("odd canonical relabelling preserves both sectors and the gradient image",
      odd_ok)

controls_ok = bool(
    provenance_ok and upstream_ok and incidence_ok and geometry_ok
    and sector_ok and flat_ok and tree_ok and odd_ok
)
corroborated = bool(
    controls_ok
    and all(rotation_ranks[str(prime)] + translation_ranks[str(prime)] == 3481
            and translation_ranks[str(prime)] + gradient_ranks[str(prime)] == 1800
            for prime in PRIMES)
)

if not controls_ok:
    outcome = "ADVERSARIAL_STATIC_DECOMPOSITION_CONTROL_FAILED"
elif corroborated:
    outcome = "ADVERSARIAL_STATIC_KERNEL_IS_P1_GRADIENTS"
else:
    outcome = "ADVERSARIAL_STATIC_DECOMPOSITION_OPEN"

allowed = {
    "ADVERSARIAL_STATIC_DECOMPOSITION_CONTROL_FAILED",
    "ADVERSARIAL_STATIC_KERNEL_IS_P1_GRADIENTS",
    "ADVERSARIAL_STATIC_EXTRA_ROTATION_OR_TRANSLATION",
    "ADVERSARIAL_STATIC_DECOMPOSITION_OPEN",
}
check("the adversarial static hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "rotation_modular_ranks": rotation_ranks,
    "translation_modular_ranks": translation_ranks,
    "gradient_modular_ranks": gradient_ranks,
    "flat_rotation_modular_ranks": flat_ranks,
    "tree_rotation_modular_ranks": tree_rotation_ranks,
    "tree_translation_modular_ranks": tree_translation_ranks,
    "odd_rotation_modular_ranks": odd_rotation_ranks,
    "odd_translation_modular_ranks": odd_translation_ranks,
    "odd_gradient_modular_ranks": odd_gradient_ranks,
    "classification": {
        "static_rotation_sector": "ZERO KERNEL" if corroborated else "OPEN",
        "static_translation_sector": (
            "CONTINUOUS P1 GRADIENTS MOD CONSTANTS" if corroborated else "OPEN"
        ),
        "static_combined_kernel_dimension": 119 if corroborated else "OPEN",
        "physical_interpretation": "OPEN",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
print("rotation:", rotation_ranks)
print("translation:", translation_ranks)
print("gradient:", gradient_ranks)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

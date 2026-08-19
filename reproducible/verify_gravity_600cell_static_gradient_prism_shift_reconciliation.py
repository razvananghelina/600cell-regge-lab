#!/usr/bin/env python3
"""Exact coordinate reconciliation of static gradients and prism shifts."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_static_gradient_prism_shift_reconciliation.json"
PROTOCOL_COMMIT = "e2c633e"
INPUT_HASHES = {
    "docs/gravity/gravity_600cell_static_gradient_prism_shift_reconciliation_prior_art.md":
        "46134f34b396ecc7dd844c85a83a561da29c8db0896d33afea8865fecc00718e",
    "docs/gravity/gravity_600cell_static_gradient_prism_shift_reconciliation_protocol.md":
        "c6bc2e4ecbc781b21ec83dc39f377131e28e95e3b990f3dc38ba46df7998302c",
    "docs/gravity/gravity_600cell_prism_shift_gluing_protocol.md":
        "f76a2216f20ea00aa66c31891d81c1f72a4ee86fe519458b6527118a4fca2251",
    "docs/gravity/gravity_600cell_prism_shift_gluing_result.md":
        "2cefd7a24a6ac132da34cbe450210446b1b3bd0b39a4920dcc36176fc4a68e1a",
    "reproducible/verify_gravity_600cell_prism_shift_gluing.py":
        "0faa50e20f3efd89b8828426d83aba5d92401bc59e72a1091653761c4ab23519",
    "reproducible/gravity_600cell_prism_shift_gluing.json":
        "1ab6654ae57c83a49dd4f427154b891c0b8ae613631773ab6733a1227b9999fa",
    "docs/gravity/gravity_600cell_static_vertex_gradient_protocol.md":
        "c689202b94abfe9436bdc8f5db6f79fa00d6a7bc743d79bc3ae34d391890c17b",
    "docs/gravity/gravity_600cell_static_vertex_gradient_adversarial_protocol.md":
        "53c094c120411e0b63cafd9b8b5e4b60880c7f1339eabec25eed7f75ccdf2805",
    "reproducible/verify_gravity_600cell_static_vertex_gradient.py":
        "6974ef8c85ee62b32daa2277ba221ee1cfa96f1c7cc7a92a8ec91fad576b124f",
    "reproducible/gravity_600cell_static_vertex_gradient.json":
        "ce018db5c66c78e89e4ca32360385955ea520b9ac8e42955b110d190432239c0",
    "reproducible/verify_gravity_600cell_static_vertex_gradient_adversarial.py":
        "ac5774353d516ee33355e7850f1f32fc6a678d28d0979bb59141607405690c5f",
    "reproducible/gravity_600cell_static_vertex_gradient_adversarial.json":
        "682ede43b77ada62ce8f7badb3c1672468ad14636e4cf7581af3e0b1a1a92632",
}
CANONICAL = tuple(sp.Matrix(point) for point in (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
))
PRIMES = (1000003, 1000033)

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


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def golden(a=0, b=0):
    return (Fraction(a), Fraction(b))


def gadd(left, right):
    return (left[0] + right[0], left[1] + right[1])


def gmul(left, right):
    a, b = left
    c, d = right
    return (a*c + b*d, a*d + b*c + b*d)


def gscale(value, scalar):
    scalar = Fraction(scalar)
    return (scalar*value[0], scalar*value[1])


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


def build_complex():
    zero = golden()
    one = golden(1)
    half = golden(Fraction(1, 2))
    phi_half = golden(0, Fraction(1, 2))
    inverse_phi_half = golden(Fraction(-1, 2), Fraction(1, 2))
    vertices = set()
    for axis in range(4):
        for sign in (-1, 1):
            point = [zero]*4
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
    incidence = {}
    for tetrahedron in tetrahedra:
        for face in combinations(tetrahedron, 3):
            incidence.setdefault(tuple(sorted(face)), []).append(tetrahedron)
    return {
        "vertices": vertices,
        "edges": tuple(edges),
        "triangles": tuple(triangles),
        "tetrahedra": tetrahedra,
        "incidence": incidence,
    }


def coordinates(tetrahedron, canonical):
    return {vertex: canonical[index]
            for index, vertex in enumerate(tetrahedron)}


def reflect(point, face_points):
    anchor, second, third = face_points
    normal = (second-anchor).cross(third-anchor)
    return sp.simplify(
        point - 2*normal*(normal.T*(point-anchor))[0]
        / (normal.T*normal)[0]
    )


def affine_map(domain, codomain):
    domain_h = sp.Matrix.hstack(*(
        point.col_join(sp.ones(1, 1)) for point in domain
    ))
    codomain_h = sp.Matrix.hstack(*(
        point.col_join(sp.ones(1, 1)) for point in codomain
    ))
    return sp.simplify(codomain_h*domain_h.inv())


def spatial_transition(source, target, shared, canonical):
    source_coordinates = coordinates(source, canonical)
    target_coordinates = coordinates(target, canonical)
    source_apex = next(vertex for vertex in source if vertex not in shared)
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
    return affine, affine[:3, :3], source_coordinates, target_coordinates


def add_value(row, column, value):
    value = sp.cancel(value)
    if not value:
        return
    updated = sp.cancel(row.get(column, 0) + value)
    if updated:
        row[column] = updated
    elif column in row:
        del row[column]


def local_edge_row(tetrahedron, tetrahedron_index, left, right):
    reference = tetrahedron[0]
    positions = {vertex: index for index, vertex in enumerate(tetrahedron)}
    row = {}
    for vertex, coefficient in ((right, 1), (left, -1)):
        if vertex != reference:
            add_value(
                row,
                3*tetrahedron_index + positions[vertex]-1,
                coefficient,
            )
    return row


def old_face_rows(complex_data):
    tetrahedra = complex_data["tetrahedra"]
    index = {tetrahedron: position
             for position, tetrahedron in enumerate(tetrahedra)}
    rows = []
    for face in sorted(complex_data["incidence"]):
        source, target = sorted(complex_data["incidence"][face])
        root, other_one, other_two = face
        for endpoint in (other_one, other_two):
            row = local_edge_row(source, index[source], root, endpoint)
            target_row = local_edge_row(target, index[target], root, endpoint)
            for column, value in target_row.items():
                add_value(row, column, -value)
            rows.append(row)
    return rows


def potential_rows(complex_data):
    rows = []
    for tetrahedron in complex_data["tetrahedra"]:
        reference = tetrahedron[0]
        for vertex in tetrahedron[1:]:
            rows.append({reference: sp.Integer(-1), vertex: sp.Integer(1)})
    return rows


def left_block_product(block, rows):
    result = []
    for start in range(0, len(rows), block.cols):
        local = rows[start:start + block.cols]
        for output in range(block.rows):
            row = {}
            for inner in range(block.cols):
                for column, value in local[inner].items():
                    add_value(row, column, block[output, inner]*value)
            result.append(row)
    return result


def transform_columns(rows, block):
    result = []
    for source in rows:
        target = {}
        for column, value in source.items():
            cell, local = divmod(column, block.rows)
            for coordinate in range(block.cols):
                add_value(
                    target,
                    block.cols*cell + coordinate,
                    value*block[local, coordinate],
                )
        result.append(target)
    return result


def transform_columns_first_cell_corrupted(rows, block, corrupted):
    result = []
    for source in rows:
        target = {}
        for column, value in source.items():
            cell, local = divmod(column, block.rows)
            active = corrupted if cell == 0 else block
            for coordinate in range(active.cols):
                add_value(
                    target,
                    active.cols*cell + coordinate,
                    value*active[local, coordinate],
                )
        result.append(target)
    return result


def new_face_rows(complex_data, canonical, identity_target=False):
    tetrahedra = complex_data["tetrahedra"]
    index = {tetrahedron: position
             for position, tetrahedron in enumerate(tetrahedra)}
    rows = []
    transition_ok = True
    identity_target_used = False
    for face in sorted(complex_data["incidence"]):
        source, target = sorted(complex_data["incidence"][face])
        affine, linear, source_coordinates, target_coordinates = (
            spatial_transition(source, target, face, canonical)
        )
        reverse = spatial_transition(target, source, face, canonical)[0]
        transition_ok &= bool(
            affine[3, :] == sp.Matrix([[0, 0, 0, 1]])
            and linear.T*linear == sp.eye(3)
            and affine*reverse == sp.eye(4)
            and all(
                affine*target_coordinates[vertex].col_join(sp.ones(1, 1))
                == source_coordinates[vertex].col_join(sp.ones(1, 1))
                for vertex in face
            )
        )
        root, other_one, other_two = face
        tangents = tuple(
            source_coordinates[endpoint]-source_coordinates[root]
            for endpoint in (other_one, other_two)
        )
        nontrivial_tangent_transport = any(
            tangent.T*linear != tangent.T for tangent in tangents
        )
        if (identity_target and not identity_target_used
                and nontrivial_tangent_transport):
            linear = sp.eye(3)
            identity_target_used = True
        for tangent in tangents:
            row = {}
            for coordinate in range(3):
                add_value(
                    row, 3*index[source] + coordinate, tangent[coordinate]
                )
                add_value(
                    row, 3*index[target] + coordinate,
                    -(tangent.T*linear)[coordinate],
                )
            rows.append(row)
    return rows, transition_ok


def rational_mod(value, prime):
    value = sp.Rational(value)
    numerator = int(value.p) % prime
    denominator = int(value.q) % prime
    if denominator == 0:
        raise ZeroDivisionError
    return numerator*pow(denominator, -1, prime) % prime


def modular_rank(rows, columns, prime):
    pivots = {}
    for original in rows:
        row = {}
        for column, value in original.items():
            residue = rational_mod(value, prime)
            if residue:
                row[column] = residue
        while row:
            old_pivots = [column for column in row if column in pivots]
            if not old_pivots:
                pivot = min(row)
                inverse = pow(row[pivot], -1, prime)
                row = {column: value*inverse % prime
                       for column, value in row.items() if value % prime}
                pivots[pivot] = row
                break
            pivot = min(old_pivots)
            factor = row[pivot]
            for column, value in pivots[pivot].items():
                updated = (row.get(column, 0)-factor*value) % prime
                if updated:
                    row[column] = updated
                elif column in row:
                    del row[column]
        if len(pivots) == columns:
            break
    return len(pivots)


def run_construction(complex_data, canonical):
    difference = sp.Matrix.hstack(*(
        canonical[index]-canonical[0] for index in range(1, 4)
    ))
    q_block = difference.T
    gradient = q_block.inv()
    b_rows = potential_rows(complex_data)
    g_rows = left_block_product(gradient, b_rows)
    qg_rows = left_block_product(q_block, g_rows)
    old_rows = old_face_rows(complex_data)
    new_rows, transition_ok = new_face_rows(complex_data, canonical)
    transformed_old = transform_columns(old_rows, q_block)
    ranks = {
        str(prime): {
            "B": modular_rank(b_rows, 120, prime),
            "G": modular_rank(g_rows, 120, prime),
            "C_old": modular_rank(old_rows, 1800, prime),
            "C_new": modular_rank(new_rows, 1800, prime),
        }
        for prime in PRIMES
    }
    return {
        "Q": q_block,
        "B": b_rows,
        "G": g_rows,
        "QG": qg_rows,
        "C_old": old_rows,
        "C_new": new_rows,
        "C_old_Q": transformed_old,
        "transition_ok": transition_ok,
        "ranks": ranks,
    }


print("="*78)
print("STATIC GRADIENT / PRISM-SHIFT EXACT RECONCILIATION")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = actual_hashes == INPUT_HASHES
check("all reconciliation inputs have frozen provenance", provenance_ok)

prism = json.loads((HERE/"gravity_600cell_prism_shift_gluing.json").read_text())
primary = json.loads((HERE/"gravity_600cell_static_vertex_gradient.json").read_text())
adversarial = json.loads((HERE/"gravity_600cell_static_vertex_gradient_adversarial.json").read_text())
upstream_ok = bool(
    prism["tests"] == prism["passed"] == 15
    and prism["verdict"] == "GLOBAL_STATIC_LENGTH_DATA_UNDERDETERMINED"
    and prism["matching"]["kernel_dimension"] == 119
    and prism["matching"]["kernel_equals_vertex_potential_image"]
    and primary["tests"] == primary["passed"] == 10
    and primary["outcome"] == "STATIC_KERNEL_EXACTLY_VERTEX_GRADIENTS"
    and adversarial["tests"] == adversarial["passed"] == 9
    and adversarial["outcome"] == "ADVERSARIAL_STATIC_KERNEL_IS_P1_GRADIENTS"
)
check("all three frozen 119-dimensional input theorems persist", upstream_ok)

complex_data = build_complex()
incidence_ok = bool(
    len(complex_data["vertices"]) == 120
    and len(complex_data["edges"]) == 720
    and len(complex_data["triangles"]) == 1200
    and len(complex_data["tetrahedra"]) == 600
    and Counter(map(len, complex_data["incidence"].values()))
    == Counter({2: 1200})
)
check("the exact carrier is the closed regular 600-cell", incidence_ok)

construction = run_construction(complex_data, CANONICAL)
q_ok = bool(construction["Q"].det() == -16)
check("the local Cartesian-to-covector map is invertible with determinant -16", q_ok,
      str(construction["Q"]))

potential_identity = construction["QG"] == construction["B"]
check("the two global potential embeddings satisfy Q G = B entrywise", potential_identity)

face_identity = construction["C_new"] == construction["C_old_Q"]
check("all 2400 face rows satisfy C_new = C_old Q entrywise", face_identity)

rank_ok = all(record == {
    "B": 119, "G": 119, "C_old": 1681, "C_new": 1681,
} for record in construction["ranks"].values())
check("both embeddings and both face operators have the exact disclosed ranks",
      rank_ok, str(construction["ranks"]))

wrong_rows, _ = new_face_rows(complex_data, CANONICAL, identity_target=True)
wrong_transition_detected = wrong_rows != construction["C_old_Q"]
check("the identity-target control breaks the first nontrivial face intertwiner",
      wrong_transition_detected)

bad_q = construction["Q"].copy()
bad_q[:, 0] = -bad_q[:, 0]
bad_qg = left_block_product(bad_q, construction["G"][:3])
bad_potential_detected = bad_qg != construction["B"][:3]
bad_face_transformed = transform_columns_first_cell_corrupted(
    construction["C_old"], construction["Q"], bad_q
)
bad_face_detected = bad_face_transformed != construction["C_new"]
check("a one-axis frame corruption breaks potential and face identities",
      bad_potential_detected and bad_face_detected)

odd = (CANONICAL[1], CANONICAL[0], CANONICAL[2], CANONICAL[3])
odd_construction = run_construction(complex_data, odd)
odd_ok = bool(
    odd_construction["Q"].det() == 16
    and odd_construction["QG"] == odd_construction["B"]
    and odd_construction["C_new"] == odd_construction["C_old_Q"]
    and odd_construction["transition_ok"]
    and odd_construction["ranks"] == construction["ranks"]
)
check("odd canonical relabelling preserves the exact carrier intertwiner", odd_ok)

controls_ok = bool(
    provenance_ok and upstream_ok and incidence_ok and q_ok
    and construction["transition_ok"] and rank_ok
    and wrong_transition_detected and bad_potential_detected and bad_face_detected
    and odd_ok
)
if not controls_ok:
    outcome = "RECONCILIATION_CONTROL_FAILED"
elif potential_identity and face_identity:
    outcome = "STATIC_GRADIENT_IS_PRISM_SHIFT_EXACTLY"
elif rank_ok:
    outcome = "STATIC_GRADIENT_ONLY_ABSTRACTLY_ISOMORPHIC"
else:
    outcome = "STATIC_GRADIENT_PRISM_SHIFT_RECONCILIATION_REFUTED"

allowed = {
    "RECONCILIATION_CONTROL_FAILED",
    "STATIC_GRADIENT_IS_PRISM_SHIFT_EXACTLY",
    "STATIC_GRADIENT_ONLY_ABSTRACTLY_ISOMORPHIC",
    "STATIC_GRADIENT_PRISM_SHIFT_RECONCILIATION_REFUTED",
}
check("the preregistered hierarchy assigns exactly one reconciliation outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "carrier_f_vector": [120, 720, 1200, 600],
    "local_intertwiner": {
        "matrix": [[str(value) for value in construction["Q"].row(index)]
                   for index in range(3)],
        "determinant": str(construction["Q"].det()),
    },
    "exact_identities": {
        "Q_G_equals_B": potential_identity,
        "C_new_equals_C_old_Q": face_identity,
        "face_rows_checked": len(construction["C_new"]),
    },
    "modular_ranks": construction["ranks"],
    "classification": {
        "carrier_identity": (
            "DERIVED EXACT RECONCILIATION"
            if outcome == "STATIC_GRADIENT_IS_PRISM_SHIFT_EXACTLY" else "OPEN"
        ),
        "new_physics": "NONE; DUPLICATE PROJECT LABELS RECONCILED",
        "physical_interpretation": "INHERITS PRISM-SHIFT RESULTS",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-"*78)
print("OUTCOME:", outcome)
print("ranks:", construction["ranks"])
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

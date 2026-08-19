#!/usr/bin/env python3
"""Exact rational certificate for the 119 static vertex-gradient modes."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_static_vertex_gradient.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_static_vertex_gradient_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_static_vertex_gradient_protocol.md"
GLOBAL_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_global_variable_face_closure_protocol.md"
GLOBAL_SOURCE = HERE / "verify_gravity_600cell_global_variable_face_closure.py"
GLOBAL_JSON = HERE / "gravity_600cell_global_variable_face_closure.json"
LOCAL_RESULT = ROOT / "docs/gravity/gravity_600cell_variable_face_connection_result.md"

PROTOCOL_COMMIT = "eff9744"
EXPECTED_HASHES = {
    "prior_art": "9541383a435e069be13ed9c2175674036a9cef5cd4e17ee455f524bd1c1c6a7d",
    "protocol": "c689202b94abfe9436bdc8f5db6f79fa00d6a7bc743d79bc3ae34d391890c17b",
    "global_protocol": "ed79c6a15ade377ae09854b3cad3028eb1c0f43cf8e85789d46993fe25ed1b49",
    "global_source": "ec44be8e4d82634e30944739d10d3f80fbb9f6fee0883ec1f612690c38d90ab6",
    "global_json": "61cebd1cd67fcdc56de088855b1fc7b805d0f70f9f9b3029d4a61209d7a53944",
    "local_result": "2db55cb87ec1c01d537cdbc11010bc9ea740762c598108e4c2de0f3acca72cc8",
}

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
CANONICAL = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
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


def kernel(matrix):
    vectors = matrix.nullspace()
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(matrix.cols, 0)


def same_space(left, right):
    return bool(
        left.rows == right.rows
        and left.rank() == right.rank()
        and left.row_join(right).rank() == left.rank()
    )


# Independent exact Q(phi) carrier.
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


def lorentz_basis(metric):
    result = []
    for a, b in ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)):
        generator = sp.zeros(4)
        generator[a, b] = 1
        generator[b, a] = -metric[a, a] / metric[b, b]
        result.append(generator)
    return tuple(result)


def static_kernel():
    result = sp.zeros(10, 6)
    result[:3, :3] = sp.eye(3)
    result[6:9, 3:6] = sp.eye(3)
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
        raise RuntimeError("nonunique lateral normal")
    normal = normals[0]
    denominator = (normal.T * metric * normal)[0]
    return sp.simplify(
        point - 2 * normal * (normal.T * metric * (point - anchor))[0]
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


def transition(source, target, shared, canonical, metric):
    source_coordinates = coordinate_map(source, canonical)
    target_coordinates = coordinate_map(target, canonical)
    source_top = {vertex: point + 5 * NORMAL
                  for vertex, point in source_coordinates.items()}
    target_top = {vertex: point + 5 * NORMAL
                  for vertex, point in target_coordinates.items()}
    source_apex = next(vertex for vertex in source if vertex not in shared)
    target_apex = next(vertex for vertex in target if vertex not in shared)
    lateral = tuple(source_coordinates[vertex] for vertex in shared) + (
        source_top[shared[0]],
    )
    reflected_lower = affine_reflection(
        source_coordinates[source_apex], lateral, metric
    )
    reflected_upper = affine_reflection(source_top[source_apex], lateral, metric)
    domain = tuple(target_coordinates[vertex] for vertex in shared) + (
        target_top[shared[0]], target_coordinates[target_apex],
    )
    codomain = lateral + (reflected_lower,)
    affine = affine_map(domain, codomain)
    linear = affine[:4, :4]
    control = bool(
        affine[4, :] == sp.Matrix([[0, 0, 0, 0, 1]])
        and sp.simplify(linear.T * metric * linear) == metric
        and all(
            sp.simplify(
                affine * target_coordinates[vertex].col_join(sp.ones(1, 1))
                - source_coordinates[vertex].col_join(sp.ones(1, 1))
            ) == sp.zeros(5, 1)
            and sp.simplify(
                affine * target_top[vertex].col_join(sp.ones(1, 1))
                - source_top[vertex].col_join(sp.ones(1, 1))
            ) == sp.zeros(5, 1)
            for vertex in shared
        )
        and sp.simplify(
            affine * target_top[target_apex].col_join(sp.ones(1, 1))
            - reflected_upper.col_join(sp.ones(1, 1))
        ) == sp.zeros(5, 1)
    )
    return {
        "linear": linear,
        "control": control,
        "source_coordinates": source_coordinates,
        "target_coordinates": target_coordinates,
        "source_top": source_top,
        "target_top": target_top,
    }


def reduced_block(source, target, shared, data, metric):
    basis = lorentz_basis(metric)
    local = static_kernel()
    source_upper = tuple(data["source_top"][vertex] for vertex in shared)
    target_upper = tuple(data["target_top"][vertex] for vertex in shared)
    source_lower = tuple(data["source_coordinates"][vertex] for vertex in shared)
    source_evaluation = poincare_evaluation(source_upper, basis)
    target_evaluation = poincare_evaluation(target_upper, basis)
    stabilizer = kernel(poincare_evaluation(source_lower, basis))
    linear_block = sp.diag(*([data["linear"]] * 3))
    pair = (source_evaluation * local).row_join(
        -linear_block * target_evaluation * local
    )
    raw = pair.row_join(-source_evaluation * stabilizer)
    raw_kernel = kernel(raw)
    projected = raw_kernel[:12, :]
    block = kernel(projected.T).T
    control = bool(
        raw.rank() == 6 and raw_kernel.shape == (13, 7)
        and block.shape == (5, 12) and block.rank() == 5
        and same_space(kernel(block), projected)
    )
    return block, control


def local_gradient(canonical):
    spatial = tuple(point[:3, :] for point in canonical)
    difference = sp.Matrix((
        (-1, 1, 0, 0),
        (-1, 0, 1, 0),
        (-1, 0, 0, 1),
    ))
    directions = sp.Matrix.hstack(*(
        spatial[index] - spatial[0] for index in range(1, 4)
    ))
    return sp.simplify(directions.T.inv() * difference)


def cell_gradient_map(tetrahedron, vertex_columns, gradient):
    selection = sp.zeros(4, len(vertex_columns))
    column_index = {vertex: index for index, vertex in enumerate(vertex_columns)}
    for local_index, vertex in enumerate(tetrahedron):
        selection[local_index, column_index[vertex]] = 1
    return sp.zeros(3, len(vertex_columns)).col_join(gradient * selection)


def build_blocks_and_gradient(complex_data, canonical, metric):
    tetrahedra = complex_data["tetrahedra"]
    tetrahedron_index = {tetrahedron: index
                         for index, tetrahedron in enumerate(tetrahedra)}
    faces = tuple(sorted(complex_data["face_to_tetrahedra"]))
    gradient = local_gradient(canonical)
    blocks = {}
    face_identity = True
    transition_control = True
    block_control = True
    discontinuity_detected = False

    for face_index, face in enumerate(faces):
        source, target = sorted(complex_data["face_to_tetrahedra"][face])
        shared = tuple(sorted(face))
        data = transition(source, target, shared, canonical, metric)
        block, local_control = reduced_block(
            source, target, shared, data, metric
        )
        transition_control &= data["control"]
        block_control &= local_control
        union = tuple(sorted(set(source) | set(target)))
        source_map = cell_gradient_map(source, union, gradient)
        target_map = cell_gradient_map(target, union, gradient)
        pair_map = source_map.col_join(target_map)
        face_identity &= bool(block * pair_map == sp.zeros(5, len(union)))
        if face_index == 0:
            wrong_target = target_map.copy()
            shared_columns = [union.index(vertex) for vertex in shared[:2]]
            first, second = shared_columns
            wrong_target[:, first], wrong_target[:, second] = (
                wrong_target[:, second], wrong_target[:, first]
            )
            wrong_pair = source_map.col_join(wrong_target)
            discontinuity_detected = bool(
                block * wrong_pair != sp.zeros(5, len(union))
            )
        blocks[face_index] = (source, target, block)

    global_gradient_rows = []
    for tetrahedron in tetrahedra:
        local_map = cell_gradient_map(
            tetrahedron, tuple(range(120)), gradient
        )
        for row_index in range(6):
            row = {
                column: sp.cancel(local_map[row_index, column])
                for column in tetrahedron
                if local_map[row_index, column] != 0
            }
            global_gradient_rows.append(row)
    return {
        "blocks": blocks,
        "tetrahedron_index": tetrahedron_index,
        "gradient": gradient,
        "gradient_rows": global_gradient_rows,
        "face_identity": face_identity,
        "transition_control": transition_control,
        "block_control": block_control,
        "discontinuity_detected": discontinuity_detected,
    }


def closure_rows(blocks, tetrahedron_index):
    rows = []
    for face_index in sorted(blocks):
        source, target, block = blocks[face_index]
        source_offset = 6 * tetrahedron_index[source]
        target_offset = 6 * tetrahedron_index[target]
        for block_row in range(5):
            row = {}
            for coordinate in range(6):
                left = sp.cancel(block[block_row, coordinate])
                right = sp.cancel(block[block_row, 6 + coordinate])
                if left:
                    row[source_offset + coordinate] = left
                if right:
                    row[target_offset + coordinate] = right
            rows.append(row)
    return rows


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
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "global_protocol": GLOBAL_PROTOCOL,
    "global_source": GLOBAL_SOURCE,
    "global_json": GLOBAL_JSON,
    "local_result": LOCAL_RESULT,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all static vertex-gradient inputs have frozen provenance",
      provenance_ok, str(hashes))

global_artifact = json.loads(GLOBAL_JSON.read_text())
upstream_ok = bool(
    global_artifact["outcome"] == "GLOBAL_VARIABLE_FACE_MODULAR_DEFICIT_OPEN"
    and global_artifact["passed"] == global_artifact["tests"] == 12
    and all(
        record["complete_modular_ranks"]
        == {str(PRIMES[0]): 3481, str(PRIMES[1]): 3481}
        for record in global_artifact["records"] if record["scale"] == 1
    )
    and all(
        set(record["complete_modular_ranks"].values()) == {3600}
        for record in global_artifact["records"] if record["scale"] != 1
    )
)
check("the target-blind static deficit and expanding full ranks persist",
      upstream_ok)

complex_data = build_complex()
incidence_ok = bool(
    len(complex_data["vertices"]) == 120
    and len(complex_data["edges"]) == 720
    and len(complex_data["triangles"]) == 1200
    and len(complex_data["tetrahedra"]) == 600
    and Counter(map(len, complex_data["face_to_tetrahedra"].values()))
    == Counter({2: 1200})
)
check("the independent exact carrier has f=(120,720,1200,600)", incidence_ok)

gradient = local_gradient(CANONICAL)
directions = sp.Matrix.hstack(*(
    CANONICAL[index][:3, :] - CANONICAL[0][:3, :]
    for index in range(1, 4)
))
differences = sp.Matrix((
    (-1, 1, 0, 0),
    (-1, 0, 1, 0),
    (-1, 0, 0, 1),
))
local_gradient_ok = bool(
    gradient.shape == (3, 4) and gradient.rank() == 3
    and gradient * sp.ones(4, 1) == sp.zeros(3, 1)
    and directions.T * gradient == differences
)
check("the canonical local P1 gradient is exact and target-independent",
      local_gradient_ok, str(gradient))

baseline = build_blocks_and_gradient(complex_data, CANONICAL, ETA)
construction_ok = bool(
    baseline["transition_control"] and baseline["block_control"]
    and baseline["face_identity"] and baseline["discontinuity_detected"]
)
check("all 1200 rational face blocks annihilate exactly the continuous gradients",
      construction_ok)

constant_zero = all(not row for row in (
    {
        0: sum(row.values())
    } if sum(row.values()) else {}
    for row in baseline["gradient_rows"]
))
gradient_ranks = {
    str(prime): modular_rank(baseline["gradient_rows"], 120, prime)
    for prime in PRIMES
}
gradient_rank_ok = constant_zero and set(gradient_ranks.values()) == {119}
check("the global gradient image has exact rational rank 119",
      gradient_rank_ok, str(gradient_ranks))

rows = closure_rows(baseline["blocks"], baseline["tetrahedron_index"])
closure_ranks = {
    str(prime): modular_rank(rows, 3600, prime) for prime in PRIMES
}
closure_rank_ok = set(closure_ranks.values()) == {3481}
check("the static closure rank lower bound is independently reproduced",
      closure_rank_ok, str(closure_ranks))

odd = (CANONICAL[1], CANONICAL[0], CANONICAL[2], CANONICAL[3])
odd_data = build_blocks_and_gradient(complex_data, odd, ETA)
odd_gradient_ranks = {
    str(prime): modular_rank(odd_data["gradient_rows"], 120, prime)
    for prime in PRIMES
}
odd_rows = closure_rows(odd_data["blocks"], odd_data["tetrahedron_index"])
odd_closure_ranks = {
    str(prime): modular_rank(odd_rows, 3600, prime) for prime in PRIMES
}
odd_ok = bool(
    odd_data["transition_control"] and odd_data["block_control"]
    and odd_data["face_identity"]
    and set(odd_gradient_ranks.values()) == {119}
    and set(odd_closure_ranks.values()) == {3481}
)
check("odd canonical relabelling preserves inclusion and exact ranks",
      odd_ok,
      f"gradient={odd_gradient_ranks}, closure={odd_closure_ranks}")

sign_data = build_blocks_and_gradient(complex_data, CANONICAL, -ETA)
sign_ok = bool(
    sign_data["transition_control"] and sign_data["block_control"]
    and sign_data["face_identity"] and sign_data["discontinuity_detected"]
)
check("metric-sign reversal preserves every face-gradient identity", sign_ok)

controls_ok = bool(
    provenance_ok and upstream_ok and incidence_ok and local_gradient_ok
    and construction_ok and gradient_rank_ok and closure_rank_ok
    and odd_ok and sign_ok
)
exhausted = bool(
    controls_ok
    and min(closure_ranks.values()) + min(gradient_ranks.values()) == 3600
)

if not controls_ok:
    outcome = "STATIC_VERTEX_GRADIENT_CONTROL_FAILED"
elif exhausted:
    outcome = "STATIC_KERNEL_EXACTLY_VERTEX_GRADIENTS"
else:
    outcome = "STATIC_VERTEX_GRADIENT_OPEN"

allowed = {
    "STATIC_VERTEX_GRADIENT_CONTROL_FAILED",
    "STATIC_KERNEL_EXACTLY_VERTEX_GRADIENTS",
    "STATIC_VERTEX_GRADIENTS_PROPER_SUBSPACE",
    "STATIC_VERTEX_GRADIENT_HYPOTHESIS_REFUTED",
    "STATIC_VERTEX_GRADIENT_OPEN",
}
check("the static vertex-gradient hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "local_gradient_matrix": [list(map(str, gradient.row(index)))
                              for index in range(gradient.rows)],
    "gradient_modular_ranks": gradient_ranks,
    "closure_modular_ranks": closure_ranks,
    "odd_gradient_modular_ranks": odd_gradient_ranks,
    "odd_closure_modular_ranks": odd_closure_ranks,
    "classification": {
        "static_rational_kernel": "CONTINUOUS P1 VERTEX GRADIENTS MOD CONSTANTS"
        if exhausted else "OPEN",
        "static_kernel_dimension": 119 if exhausted else "OPEN",
        "physical_interpretation": "OPEN",
        "action_hessian_or_dynamics": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
print("gradient ranks:", gradient_ranks)
print("closure ranks:", closure_ranks)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

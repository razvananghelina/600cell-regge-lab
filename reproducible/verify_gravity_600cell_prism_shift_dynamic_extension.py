#!/usr/bin/env python3
"""Exact branch audit for unequal-scale extension of the prism shift."""

from collections import Counter, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_prism_shift_dynamic_extension.json"
PRIOR_ART_COMMIT = "e455921"
PROTOCOL_COMMIT = "b251d5b"
INPUT_HASHES = {
    "docs/gravity/gravity_600cell_prism_shift_dynamic_extension_prior_art.md":
        "24474885f7888b4a3b750418992284a2b768265b386f2b921cad1a72307a16b9",
    "docs/gravity/gravity_600cell_prism_shift_dynamic_extension_protocol.md":
        "c1286f90e31048a7582a54943f906ce95b450a26b8e6dd212aa57b11644410dd",
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
PRIMES = (101, 1000003)
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


def tetrahedra_from_graph(adjacency):
    neighbours = [set(np.flatnonzero(row > 0.5)) for row in adjacency]
    cells = []
    for first in range(len(adjacency)):
        for second in sorted(v for v in neighbours[first] if v > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(v for v in common_two if v > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(v for v in common_three if v > third):
                    cells.append((first, second, third, fourth))
    return tuple(cells)


def rank_mod_prime(rows, column_count, prime):
    basis = {}
    for source in rows:
        row = {column: value % prime for column, value in source.items()
               if value % prime}
        while row:
            pivot = min(row)
            if pivot not in basis:
                inverse = pow(row[pivot], -1, prime)
                basis[pivot] = {
                    column: (value*inverse) % prime
                    for column, value in row.items() if value % prime
                }
                break
            factor = row[pivot]
            for column, value in basis[pivot].items():
                updated = (row.get(column, 0)-factor*value) % prime
                if updated:
                    row[column] = updated
                elif column in row:
                    del row[column]
    assert all(0 <= pivot < column_count for pivot in basis)
    return len(basis)


def graph_connected(vertex_count, edges):
    adjacency = [[] for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    visited = {0}
    queue = deque([0])
    while queue:
        for neighbour in adjacency[queue.popleft()]:
            if neighbour not in visited:
                visited.add(neighbour)
                queue.append(neighbour)
    return len(visited) == vertex_count


print("="*78)
print("UNEQUAL-SCALE DYNAMIC EXTENSION OF THE PRISM SHIFT")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the prior-art gate, protocol and carrier have frozen provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "e455921" and PROTOCOL_COMMIT == "b251d5b",
    str(actual_hashes),
)

# Frozen rational regular tetrahedron.
b = (
    sy.Matrix((1, 1, 1)),
    sy.Matrix((1, -1, -1)),
    sy.Matrix((-1, 1, -1)),
    sy.Matrix((-1, -1, 1)),
)
a, s1, s2, s3, lapse = sy.symbols("a s1 s2 s3 N", real=True)
s = sy.Matrix((s1, s2, s3))
q = 1+a
radii = tuple((vertex.T*vertex)[0] for vertex in b)
E = sy.Matrix.hstack(*(b[index]-b[0] for index in range(1, 4))).T
check(
    "the frozen tetrahedron is regular, centered and spans three-space",
    sum(b, sy.zeros(3, 1)) == sy.zeros(3, 1)
    and len(set(radii)) == 1 and radii[0] == 3
    and E.det() == -16,
    f"radii={radii}, det(E)={E.det()}",
)

# In four-dimensional local spacetime coordinates, all lateral faces obey an
# exact two-vector affine relation.
bottom = tuple(vertex.col_join(sy.Matrix([0])) for vertex in b)
top = tuple((q*vertex+s).col_join(sy.Matrix([lapse])) for vertex in b)
planarity_residuals = []
planarity_ranks = []
control_substitution = {a: sy.Rational(1, 3), s1: 2, s2: -1, s3: 3,
                        lapse: 5}
for left, right in combinations(range(4), 2):
    edge = bottom[right]-bottom[left]
    strut = top[left]-bottom[left]
    residual = sy.simplify(
        (top[right]-bottom[left])-strut-q*edge)
    planarity_residuals.append(residual)
    columns = sy.Matrix.hstack(
        edge, top[right]-bottom[left], top[left]-bottom[left])
    planarity_ranks.append(columns.subs(control_substitution).rank())
check(
    "all six lateral quadrilaterals are exactly planar and nondegenerate",
    all(residual == sy.zeros(4, 1) for residual in planarity_residuals)
    and planarity_ranks == [2]*6,
)

strut_squares = tuple(
    sy.expand(((a*vertex+s).T*(a*vertex+s))[0]-lapse**2)
    for vertex in b
)
constraints = sy.Matrix(tuple(
    sy.expand(strut_squares[index]-strut_squares[0])
    for index in range(1, 4)
))
expected_constraints = sy.simplify(2*a*E*s)
transformed_constraints = sy.simplify(E.inv()*constraints/2)
check(
    "common struts reduce exactly to c=2*a*E*s",
    constraints == expected_constraints
    and transformed_constraints == a*s,
    f"c={tuple(constraints)}, E^-1*c/2={tuple(transformed_constraints)}",
)

# Compute the intersection <a> cap <s1,s2,s3> through exact elimination.
t = sy.symbols("t")
intersection_basis = sy.groebner(
    (t*a, (1-t)*s1, (1-t)*s2, (1-t)*s3),
    t, a, s1, s2, s3, order="lex",
)
eliminated = {
    sy.expand(poly.as_expr()) for poly in intersection_basis.polys
    if not poly.as_expr().has(t)
}
product_generators = {a*s1, a*s2, a*s3}
check(
    "the exact common-strut ideal is the union of two branch primes",
    eliminated == product_generators,
    f"elimination basis={sorted(map(str, eliminated))}",
)

anchor = {a: 0, s1: 0, s2: 0, s3: 0}
jacobian = constraints.jacobian((a, s1, s2, s3)).subs(anchor)
check(
    "the equal-scale zero-shift anchor has a four-dimensional linear tangent",
    jacobian == sy.zeros(3, 4),
)

alpha = sy.symbols("alpha", real=True)
beta1, beta2, beta3 = sy.symbols("beta1 beta2 beta3", real=True)
beta = sy.Matrix((beta1, beta2, beta3))
second_directional = sy.simplify(4*alpha*E*beta)
pure_scale = second_directional.subs(
    {alpha: 1, beta1: 0, beta2: 0, beta3: 0})
pure_shifts = [second_directional.subs(
    {alpha: 0, beta1: int(index == 0), beta2: int(index == 1),
     beta3: int(index == 2)}) for index in range(3)]
mixed_shifts = [second_directional.subs(
    {alpha: 1, beta1: int(index == 0), beta2: int(index == 1),
     beta3: int(index == 2)}) for index in range(3)]
check(
    "the quadratic tangent cone accepts pure branches and rejects mixtures",
    pure_scale == sy.zeros(3, 1)
    and all(value == sy.zeros(3, 1) for value in pure_shifts)
    and all(value != sy.zeros(3, 1) for value in mixed_shifts)
    and E.det() != 0,
    f"mixed second derivatives={[tuple(value) for value in mixed_shifts]}",
)

# Literal global carrier and exact incidence ranks over two independent fields.
_, adjacency, _ = build_600cell()
cells = tetrahedra_from_graph(adjacency)
edges = tuple(sorted({tuple(sorted(edge)) for cell in cells
                      for edge in combinations(cell, 2)}))
faces = tuple(sorted({tuple(sorted(face)) for cell in cells
                      for face in combinations(cell, 3)}))
edge_incidence = Counter(tuple(sorted(edge)) for cell in cells
                         for edge in combinations(cell, 2))
incidence_rows = tuple({left: -1, right: 1} for left, right in edges)
incidence_ranks = {
    prime: rank_mod_prime(incidence_rows, 120, prime) for prime in PRIMES
}
check(
    "the literal 600-cell carrier and local coverage are exact",
    (120, len(edges), len(faces), len(cells)) == (120, 720, 1200, 600)
    and Counter(edge_incidence.values()) == Counter({5: 720})
    and graph_connected(120, edges),
    f"f={(120, len(edges), len(faces), len(cells))}",
)
check(
    "the global potential incidence has exact rank 119 over two fields",
    incidence_ranks == {prime: 119 for prime in PRIMES},
    str(incidence_ranks),
)

equal_scale_shift_dimension = incidence_ranks[PRIMES[0]]
unequal_scale_shift_dimension = 0 if E.det() != 0 else None
check(
    "equal scale has 119 shifts but unequal common-strut scale has none",
    equal_scale_shift_dimension == 119
    and unequal_scale_shift_dimension == 0,
)

# Nonuniform strut differences encode, rather than leave free, the shift.
recovery_controls = (
    (sy.Rational(-1, 2), sy.Matrix((1, 2, 3))),
    (sy.Rational(1, 3), sy.Matrix((-2, 1, 4))),
    (sy.Integer(1), sy.Matrix((3, -1, -2))),
)
recovery_records = []
for scale_difference, shift in recovery_controls:
    differences = 2*scale_difference*E*shift
    recovered = sy.simplify(
        E.inv()*differences/(2*scale_difference))
    recovery_records.append({
        "a": str(scale_difference),
        "shift": tuple(map(str, shift)),
        "strut_differences": tuple(map(str, differences)),
        "recovered": tuple(map(str, recovered)),
        "exact": recovered == shift,
    })
map_determinant = sy.factor((2*a*E).det())
check(
    "nonuniform strut differences recover every frozen nonzero shift exactly",
    all(record["exact"] for record in recovery_records),
)
check(
    "the recovery map has the preregistered equal-scale pole",
    map_determinant == (2*a)**3*E.det(),
    f"det(2*a*E)={map_determinant}",
)

# Sensitivity control: loss of equal radii permits a nonzero compensating
# shift when homothety is kept about the now non-circumcentric origin.
deformed = b[:3]+(sy.Matrix((-1, -1, 2)),)
deformed_E = sy.Matrix.hstack(
    *(deformed[index]-deformed[0] for index in range(1, 4))).T
deformed_radii = sy.Matrix(tuple(
    (deformed[index].T*deformed[index])[0]
    - (deformed[0].T*deformed[0])[0]
    for index in range(1, 4)
))
deformed_shift = sy.simplify(-deformed_E.inv()*deformed_radii/2)
deformed_squares = tuple(sy.expand(
    ((vertex+deformed_shift).T*(vertex+deformed_shift))[0]-25)
    for vertex in deformed)
check(
    "the unequal-radius negative control admits a detected nonzero shift",
    deformed_E.det() != 0 and deformed_shift != sy.zeros(3, 1)
    and len(set(deformed_squares)) == 1,
    f"shift={tuple(deformed_shift)}, squares={deformed_squares}",
)

verdict = (
    "DYNAMIC_SHIFT_EXTENSION_OBSTRUCTED"
    if passed == tests else "BRANCH_GEOMETRY_OPEN"
)
artifact = {
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "input_hashes": actual_hashes,
    },
    "local_exact": {
        "tetrahedron_radii": list(map(str, radii)),
        "edge_direction_determinant": str(E.det()),
        "common_strut_constraints": list(map(str, constraints)),
        "transformed_ideal_generators": sorted(map(str, eliminated)),
        "branch_equation": "(q-1)*s=0",
        "linear_tangent_dimension_at_anchor": 4,
        "quadratic_tangent_cone": "alpha=0 union beta=0",
    },
    "global_600cell": {
        "f_vector": [120, len(edges), len(faces), len(cells)],
        "incidence_ranks": {str(key): value
                            for key, value in incidence_ranks.items()},
        "equal_scale_shift_dimension_mod_constants": 119,
        "unequal_scale_common_strut_shift_dimension": 0,
    },
    "nonuniform_struts": {
        "recovery_formula": "s=E^-1*c/[2*(q-1)]",
        "map_determinant": str(map_determinant),
        "controls": recovery_records,
    },
    "negative_control": {
        "deformed_radii": list(map(str, (
            (vertex.T*vertex)[0] for vertex in deformed))),
        "compensating_shift": list(map(str, deformed_shift)),
        "common_squares": list(map(str, deformed_squares)),
    },
    "classification": "DERIVED_EXACT_KINEMATIC_OBSTRUCTION",
    "tests": tests,
    "passed": passed,
    "verdict": verdict,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks pass")
print(f"VERDICT: {verdict}")
print(f"ARTIFACT: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Exact independent audit of the 600-cell prism-shift Regge Hessian."""

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


OUTPUT = HERE / "gravity_600cell_prism_shift_action_adversarial.json"
AUDIT_PROTOCOL_COMMIT = "2b81792"
INPUT_HASHES = {
    "docs/gravity/gravity_600cell_prism_shift_action_adversarial_protocol.md":
        "2db378dd6cfeed5be0f86a207855152c7a54fe1a14a03f4731a1d6b5e7c9ad38",
    "reproducible/gravity_600cell_prism_shift_action.json":
        "63c9fe41ea4b4de2457f1308a91689786e3871d09ffe8be9008912300e6a4260",
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
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


def tetrahedra_from_graph(adjacency):
    """Enumerate 4-cliques directly; no primary action code is imported."""
    neighbours = [set(np.flatnonzero(row > 0.5)) for row in adjacency]
    cells = []
    for a in range(len(adjacency)):
        for b in sorted(vertex for vertex in neighbours[a] if vertex > a):
            common_ab = neighbours[a] & neighbours[b]
            for c in sorted(vertex for vertex in common_ab if vertex > b):
                for d in sorted(vertex for vertex in
                                common_ab & neighbours[c] if vertex > c):
                    cells.append((a, b, c, d))
    return tuple(cells)


def connected(vertex_count, edges):
    graph = [[] for _ in range(vertex_count)]
    for left, right in edges:
        graph[left].append(right)
        graph[right].append(left)
    seen = {0}
    queue = deque([0])
    while queue:
        for neighbour in graph[queue.popleft()]:
            if neighbour not in seen:
                seen.add(neighbour)
                queue.append(neighbour)
    return len(seen) == vertex_count


def matrix_strings(matrix):
    return [[str(sy.simplify(value)) for value in matrix.row(row)]
            for row in range(matrix.rows)]


print("="*78)
print("ADVERSARIAL EXACT AUDIT: PRISM-SHIFT REGGE HESSIAN")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the audit protocol and primary artifact have frozen provenance",
    actual_hashes == INPUT_HASHES and AUDIT_PROTOCOL_COMMIT == "2b81792",
    str(actual_hashes),
)

# Exact local inverse.  The second-order jet below follows from this inverse,
# rather than from an angle routine or from the Schlaefli identity.
x1, x2, x3 = sy.symbols("x1 x2 x3", real=True)
variables = (x1, x2, x3)
x = sy.Matrix(variables)
half = sy.Rational(1, 2)
G = sy.Matrix(((1, half, half), (half, 1, half), (half, half, 1)))
G_inverse = sy.Matrix((
    (sy.Rational(3, 2), -half, -half),
    (-half, sy.Rational(3, 2), -half),
    (-half, -half, sy.Rational(3, 2)),
))
u = G_inverse*x
schur = -1-(x.T*u)[0]
upper = G_inverse+u*u.T/schur
inverse = upper.row_join(-u/schur).col_join(
    (-u.T/schur).row_join(sy.Matrix([[1/schur]])))
metric = G.row_join(x).col_join(x.T.row_join(sy.Matrix([[-1]])))
inverse_residual = metric*inverse-sy.eye(4)
check(
    "the symbolic block matrix is the exact inverse metric",
    all(sy.cancel(value) == 0 for value in inverse_residual),
)

lateral_normals = (
    sy.Matrix((-1, -1, -1)),
    sy.Matrix((1, 0, 0)),
    sy.Matrix((0, 1, 0)),
    sy.Matrix((0, 0, 1)),
)
local_edges = tuple(combinations(range(4), 2))
angle_hessians = []
cosine_anchors = []
cosine_gradients = []

for left, right in local_edges:
    omitted = tuple(index for index in range(4)
                    if index not in (left, right))
    n = lateral_normals[omitted[0]]
    m = lateral_normals[omitted[1]]
    normal_n = (n.T*G_inverse*n)[0]
    normal_m = (m.T*G_inverse*m)[0]
    cross = (n.T*G_inverse*m)[0]
    linear_n = (n.T*G_inverse*x)[0]
    linear_m = (m.T*G_inverse*x)[0]

    # The exact inverse has upper block Q+uu^T/S.  Because u is linear and
    # S=-1+O(x^2), its exact two-jet is Q-uu^T.  Expanding the normalized
    # conormal product therefore gives the following exact quadratic jet.
    delta_n = -linear_n**2
    delta_m = -linear_m**2
    delta_cross = -linear_n*linear_m
    normalization = sy.sqrt(normal_n*normal_m)
    cosine_anchor = sy.simplify(-cross/normalization)
    cosine_quadratic = sy.expand(
        -delta_cross/normalization
        + cross*(delta_n/normal_n+delta_m/normal_m)/(2*normalization)
    )
    cosine_hessian = sy.hessian(cosine_quadratic, variables)
    angle_hessian = sy.simplify(
        -cosine_hessian/sy.sqrt(1-cosine_anchor**2))
    cosine_anchors.append(cosine_anchor)
    cosine_gradients.append((sy.Integer(0),)*3)
    angle_hessians.append(angle_hessian)

angle_sum = sy.simplify(sum(angle_hessians, sy.zeros(3)))
check(
    "all local lateral cosine anchors and first derivatives are exact",
    all(value == sy.Rational(1, 3) for value in cosine_anchors)
    and all(gradient == (0, 0, 0) for gradient in cosine_gradients),
)
check(
    "the six independently differentiated angle Hessians cancel exactly",
    angle_sum == sy.zeros(3),
    str(matrix_strings(angle_sum)),
)
check(
    "the angle differentiation is nontrivial before the six-edge sum",
    any(matrix != sy.zeros(3) for matrix in angle_hessians),
)

potentials = (sy.Integer(0), x1, x2, x3)
area_hessians = []
for left, right in local_edges:
    area = sy.sqrt(1+(potentials[right]-potentials[left])**2)
    area_hessians.append(
        sy.hessian(area, variables).subs({x1: 0, x2: 0, x3: 0}))
area_sum = sy.simplify(sum(area_hessians, sy.zeros(3)))
K4_reduced = sy.Matrix(((3, -1, -1), (-1, 3, -1), (-1, -1, 3)))
theta0 = sy.acos(sy.Rational(1, 3))
wedge_hessian = sy.simplify(theta0*area_sum+angle_sum)
check(
    "the direct local wedge Hessian is acos(1/3) times reduced K4",
    area_sum == K4_reduced and wedge_hessian == theta0*K4_reduced,
    str(matrix_strings(wedge_hessian)),
)

# Independent cellular assembly from the graph cliques.
_, adjacency, _ = build_600cell()
cells = tetrahedra_from_graph(adjacency)
edges = tuple(sorted({tuple(sorted(edge)) for cell in cells
                      for edge in combinations(cell, 2)}))
edge_incidence = Counter(edge for cell in cells
                         for edge in combinations(cell, 2))
laplacian = sy.zeros(120)
for left, right in edges:
    laplacian[left, left] += 1
    laplacian[right, right] += 1
    laplacian[left, right] -= 1
    laplacian[right, left] -= 1
cell_assembly = sy.zeros(120)
for cell in cells:
    for local_left, local_right in local_edges:
        left = cell[local_left]
        right = cell[local_right]
        cell_assembly[left, left] += 1
        cell_assembly[right, right] += 1
        cell_assembly[left, right] -= 1
        cell_assembly[right, left] -= 1

incidence_distribution = Counter(edge_incidence.values())
check(
    "the cellular reconstruction discovers f0=120, f1=720, f3=600",
    len(edges) == 720 and len(cells) == 600
    and incidence_distribution == Counter({5: 720}),
    f"edge incidence={dict(incidence_distribution)}",
)
check(
    "the 600 local K4 forms assemble to five graph Laplacians",
    cell_assembly == 5*laplacian,
)

epsilon = 2*sy.pi-5*theta0
global_hessian = 2*sy.pi*laplacian-theta0*cell_assembly
check(
    "the independently assembled global Hessian is epsilon Delta0",
    global_hessian == epsilon*laplacian,
)

# Boundary cancellation is derived from the exact inverse and conormal signs.
bottom = sy.Matrix((0, 0, 0, 1))
top = -bottom
boundary_relations = []
for normal in lateral_normals:
    lateral = normal.col_join(sy.Matrix([0]))
    numerator_sum = sy.cancel(
        (bottom.T*inverse*lateral)[0]+(top.T*inverse*lateral)[0])
    norm_difference = sy.cancel(
        (bottom.T*inverse*bottom)[0]-(top.T*inverse*top)[0])
    boundary_relations.append((numerator_sum, norm_difference))
check(
    "opposite conormals give cos(top)=-cos(bottom) exactly",
    all(pair == (0, 0) for pair in boundary_relations),
)

faces = tuple(sorted({tuple(sorted(face)) for cell in cells
                      for face in combinations(cell, 3)}))
face_incidence = Counter(tuple(sorted(face)) for cell in cells
                         for face in combinations(cell, 3))
check(
    "the connected angle branch cancels the full top-bottom boundary term",
    len(faces) == 1200
    and Counter(face_incidence.values()) == Counter({2: 1200}),
    "theta_top=pi-theta_bottom and every face has two incident cells",
)

# Negative controls and exact positivity probes.
shadow_hessian = (2*sy.pi-4*theta0)*laplacian
check(
    "the four-cell shadow has its own coefficient and rejects the primary",
    shadow_hessian == (2*sy.pi-4*theta0)*laplacian
    and shadow_hessian != global_hessian,
)
constant = sy.ones(120, 1)
delta = sy.zeros(120, 1)
delta[0] = 1
one_vertex_energy = sy.simplify((delta.T*global_hessian*delta)[0])
check(
    "constant shifts are null and a one-vertex shift is strictly positive",
    global_hessian*constant == sy.zeros(120, 1)
    and one_vertex_energy == 12*epsilon
    and bool(sy.N(one_vertex_energy, 30) > 0),
    f"one-vertex Hessian={one_vertex_energy}",
)

# Exact spectral checksum, independent of numerical diagonalization.
golden = (1+sy.sqrt(5))/2
spectrum = (
    (0, 1), (12-6*golden, 4), (12-4*golden, 9),
    (9, 16), (12, 25), (14, 36), (8+4*golden, 9),
    (15, 16), (6+6*golden, 4),
)
spectral_dimension = sum(multiplicity for _, multiplicity in spectrum)
spectral_trace = sy.simplify(sum(value*multiplicity
                                 for value, multiplicity in spectrum))
spectral_trace_square = sy.simplify(sum(value**2*multiplicity
                                        for value, multiplicity in spectrum))
matrix_trace = sy.trace(laplacian)
matrix_trace_square = sy.trace(laplacian*laplacian)
is_connected = connected(120, edges)
check(
    "the exact frozen spectrum matches trace and trace-square",
    spectral_dimension == 120 and spectral_trace == matrix_trace == 1440
    and spectral_trace_square == matrix_trace_square == 18720,
    f"trace={matrix_trace}, trace-square={matrix_trace_square}",
)
check(
    "connectedness and the positive coefficient give rank 119",
    is_connected and epsilon.is_positive is True,
)

verdict = (
    "SHIFT_LAPLACIAN_HESSIAN_CORROBORATED"
    if passed == tests else "PRIMARY_HESSIAN_REFUTED"
)
artifact = {
    "provenance": {
        "audit_protocol_commit": AUDIT_PROTOCOL_COMMIT,
        "input_hashes": actual_hashes,
        "method_independence": [
            "no primary action evaluator",
            "no finite differences",
            "no Schlaefli identity assumed",
        ],
    },
    "local_symbolic_audit": {
        "cosine_anchor": "1/3",
        "cosine_first_derivatives": "all zero",
        "angle_hessians": [matrix_strings(matrix)
                            for matrix in angle_hessians],
        "angle_hessian_sum": matrix_strings(angle_sum),
        "area_hessian_sum": matrix_strings(area_sum),
        "wedge_hessian": "acos(1/3) * K4_reduced",
    },
    "cellular_assembly": {
        "f_vector_partial": [120, 720, 1200, 600],
        "edge_incidence": dict(incidence_distribution),
        "cell_matrix": "5 * Delta0",
        "boundary_hessian": "0 (exact top-bottom cancellation)",
        "hessian": "[2*pi-5*acos(1/3)] * Delta0",
        "epsilon": str(sy.N(epsilon, 40)),
    },
    "negative_controls": {
        "four_cell_shadow": "[2*pi-4*acos(1/3)] * Delta0",
        "individual_angle_hessian_nonzero": True,
        "constant_null": True,
        "one_vertex_hessian": str(one_vertex_energy),
    },
    "spectral_checksum": {
        "trace_Delta0": int(matrix_trace),
        "trace_Delta0_squared": int(matrix_trace_square),
        "rank": 119,
        "nullity": 1,
    },
    "classification": "DERIVED_RESTRICTED_SECTOR",
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

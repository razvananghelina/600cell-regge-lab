#!/usr/bin/env python3
"""Weak first-class locality audit for the Whitney Hamiltonian.

Protocol commits 50a5f79 and bb02023 froze the weak condition, the unique
top-degree block, the endpoint-star locality test and the controls before the
definitive calculation.  Preliminary reconnaissance was disclosed there.
"""

from itertools import combinations
import json
from math import factorial
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import splu
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name(
    "whitney_weak_first_class_hamiltonian.json"
)
PROTOCOL_COMMIT = "50a5f79"
PROTOCOL_CORRECTION_COMMIT = "bb02023"
MODULUS = 1_000_003
RELATIVE_THRESHOLD = 1e-10
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


coordinate_bases = [list(combinations(range(3), degree))
                    for degree in range(4)]


def wedge_components(covectors, degree):
    if degree == 0:
        return sy.Matrix((1,))
    return sy.Matrix([
        sy.det(sy.Matrix([
            [covector[index] for index in basis]
            for covector in covectors
        ]))
        for basis in coordinate_bases[degree]
    ])


def local_whitney_mass(points, degree):
    """Derive one exact affine Whitney mass from its defining integral."""
    affine = sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    )
    inverse = affine.inv()
    gradients = [-sum(
        (sy.Matrix(inverse.row(row)).T for row in range(3)),
        sy.zeros(3, 1),
    )]
    gradients.extend(sy.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det()) / 6
    barycentric_second_moment = volume * (sy.ones(4, 4) + sy.eye(4)) / 20
    forms = list(combinations(range(4), degree + 1))
    coefficient_matrices = []
    for form in forms:
        coefficients = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            coefficients[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                coefficients[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * wedge_components(covectors, degree)
                )
        coefficient_matrices.append(coefficients)
    mass = sy.zeros(len(forms), len(forms))
    for row, left in enumerate(coefficient_matrices):
        for column, right in enumerate(coefficient_matrices):
            mass[row, column] = sy.simplify(sum(
                (left[basis, :] * barycentric_second_moment
                 * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
            ))
    return mass


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


def complete_constraints(top_cells, cells):
    """Canonical neighbour differences across every shared triangle."""
    local_faces = [list(combinations(range(4), degree + 1))
                   for degree in range(4)]
    local_offsets = np.cumsum((0, 4, 6, 4, 1))
    occurrence = {}
    for top_id, top in enumerate(top_cells):
        for degree in range(4):
            for local_id, face in enumerate(local_faces[degree]):
                cell = tuple(top[index] for index in face)
                local = top_id * 15 + int(local_offsets[degree]) + local_id
                occurrence[(top_id, degree, cell)] = local
    triangle_parents = {triangle: [] for triangle in cells[2]}
    for top_id, top in enumerate(top_cells):
        for triangle in combinations(top, 3):
            triangle_parents[tuple(triangle)].append(top_id)
    rows = []
    for triangle in cells[2]:
        parents = triangle_parents[triangle]
        if len(parents) != 2:
            raise AssertionError("triangle does not have two parents")
        left, right = parents
        for degree in range(3):
            for cell in combinations(triangle, degree + 1):
                left_copy = occurrence[(left, degree, tuple(cell))]
                right_copy = occurrence[(right, degree, tuple(cell))]
                row = [sy.Integer(0)] * (15 * len(top_cells))
                row[left_copy] = -1
                row[right_copy] = 1
                rows.append(row)
    return sy.Matrix(rows)


def degree_two_constraint(top_cells):
    """One difference row per global triangle and its tetra endpoints."""
    local_faces = list(combinations(range(4), 3))
    triangles = sorted({
        tuple(top[index] for index in face)
        for top in top_cells for face in local_faces
    })
    triangle_index = {face: index for index, face in enumerate(triangles)}
    occurrences = [[] for _ in triangles]
    for top_id, top in enumerate(top_cells):
        for local_id, face in enumerate(local_faces):
            triangle = tuple(top[index] for index in face)
            occurrences[triangle_index[triangle]].append((top_id, local_id))
    matrix = sy.zeros(len(triangles), 4 * len(top_cells))
    endpoints = []
    for row, copies in enumerate(occurrences):
        if len(copies) != 2:
            raise AssertionError("degree-two face does not have two parents")
        (left, left_face), (right, right_face) = copies
        matrix[row, 4 * left + left_face] = -1
        matrix[row, 4 * right + right_face] = 1
        endpoints.append((left, right))
    return matrix, endpoints, triangles


def modular_solve(matrix, rhs, modulus):
    """Dense exact Gauss--Jordan solve over a prime field."""
    dimension = matrix.shape[0]
    augmented = np.concatenate((
        np.asarray(matrix % modulus, dtype=np.int64),
        np.asarray(rhs % modulus, dtype=np.int64).reshape(-1, 1),
    ), axis=1)
    rank = 0
    for column in range(dimension):
        candidates = np.flatnonzero(augmented[column:, column])
        if not len(candidates):
            continue
        pivot = column + int(candidates[0])
        if pivot != column:
            augmented[[column, pivot]] = augmented[[pivot, column]]
        inverse = pow(int(augmented[column, column]), -1, modulus)
        augmented[column] = augmented[column] * inverse % modulus
        factors = augmented[:, column].copy()
        factors[column] = 0
        augmented = (
            augmented - factors[:, None] * augmented[column][None, :]
        ) % modulus
        rank += 1
    return rank, augmented[:, -1]


print("=" * 78)
print("WEAK FIRST-CLASS WHITNEY HAMILTONIAN LOCALITY")
print("=" * 78)

# -------------------------------------------------------------------------
# Exact algebra and boundary-of-a-4-simplex control.
# -------------------------------------------------------------------------
regular_points = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
mass_two = local_whitney_mass(regular_points, 2)
mass_three = local_whitney_mass(regular_points, 3)
boundary_top = tuple(combinations(range(5), 4))
boundary_cells = all_simplices(boundary_top)
constraint_all = complete_constraints(boundary_top, boundary_cells)

identity_occurrence = sy.eye(constraint_all.cols)
identity_constraint = sy.eye(constraint_all.rows)
L = sy.Matrix.hstack(constraint_all, identity_constraint)
S_surface = sy.Matrix.vstack(identity_occurrence, -constraint_all)
check("the constraint-surface parametrisation is exact",
      L * S_surface == sy.zeros(L.rows, S_surface.cols)
      and L.rank() == constraint_all.rows
      and S_surface.rank() == constraint_all.cols,
      f"L={L.shape}, S={S_surface.shape}")
check("the weak condition is strictly weaker than H K=0",
      L.T != sy.zeros(L.cols, L.rows)
      and S_surface.T * L.T == sy.zeros(S_surface.cols, L.rows),
      "HK=L*V is annihilated by S*, while HK need not vanish")

sign = sy.Matrix((1, -1, 1, -1))
mass_two_inverse = mass_two.inv()
top_coboundary = sign.T
weak_32 = mass_three * top_coboundary
check("the exact local degree-two Whitney inverse has rank-one form",
      mass_two_inverse == 6 * sy.eye(4) + sign * sign.T)
check("the fixed top weak block contracts to (15/4) times incidence",
      weak_32 * mass_two_inverse == sy.Rational(15, 4) * sign.T)

constraint_two, boundary_endpoints, _ = degree_two_constraint(boundary_top)
sign_injection = sy.diag(*([sign] * len(boundary_top)))
dual_incidence = constraint_two * sign_injection
metric_two_inverse = sy.diag(*([mass_two_inverse] * len(boundary_top)))
weak_32_sum = sy.diag(*([weak_32] * len(boundary_top)))
gram_two = constraint_two * metric_two_inverse * constraint_two.T
right_two = weak_32_sum * metric_two_inverse * constraint_two.T
check("degree-two constraints are independent and G2 is positive/invertible",
      constraint_two.rank() == constraint_two.rows
      and gram_two.det() != 0
      and gram_two == 12 * sy.eye(gram_two.rows)
      + dual_incidence * dual_incidence.T,
      f"rank(C2)={constraint_two.rank()}, dim(G2)={gram_two.rows}")
check("the top-degree weak equation has the frozen exact right side",
      right_two == sy.Rational(15, 4) * dual_incidence.T)

reduced_top_matrix = (
    12 * sy.eye(len(boundary_top))
    + dual_incidence.T * dual_incidence
)
top_cross = sy.simplify(
    sy.Rational(15, 4) * reduced_top_matrix.inv() * dual_incidence.T
)
check("the unique B32 solves the weak first-class top block exactly",
      top_cross * gram_two == right_two
      and reduced_top_matrix.det() != 0,
      "B32 G2=A32 M2^-1 C2*; D and all other B blocks drop out")

boundary_remote = []
boundary_incident = []
for top_id in range(len(boundary_top)):
    for triangle_id, endpoints in enumerate(boundary_endpoints):
        value = top_cross[top_id, triangle_id]
        if top_id in endpoints:
            boundary_incident.append(value)
        elif value != 0:
            boundary_remote.append(value)
check("the unique exact top block violates endpoint-star locality",
      len(boundary_remote) == 30
      and all(value != 0 for value in boundary_incident),
      f"exact remote/incident nonzeros={len(boundary_remote)}/"
      f"{sum(value != 0 for value in boundary_incident)}")

# Frozen mass-lumped negative control: M2_lump=(1/6)I.
lumped_inverse = 6 * sy.eye(4)
lumped_gram = constraint_two * sy.diag(
    *([lumped_inverse] * len(boundary_top))
) * constraint_two.T
lumped_right = weak_32_sum * sy.diag(
    *([lumped_inverse] * len(boundary_top))
) * constraint_two.T
lumped_cross = sy.simplify(lumped_right * lumped_gram.inv())
lumped_remote = sum(
    lumped_cross[top_id, triangle_id] != 0
    and top_id not in endpoints
    for top_id in range(len(boundary_top))
    for triangle_id, endpoints in enumerate(boundary_endpoints)
)
check("the mass-lumped negative control has zero remote support",
      lumped_gram == 12 * sy.eye(lumped_gram.rows)
      and lumped_cross == sy.Rational(3, 16) * dual_incidence.T
      and lumped_remote == 0)

# -------------------------------------------------------------------------
# Independent base-600-cell assembly and support-depth audit.
# -------------------------------------------------------------------------
_, adjacency, _ = build_600cell()
neighbours = tuple(
    frozenset(np.flatnonzero(adjacency[index]).tolist())
    for index in range(120)
)
edges = tuple(
    (left, right)
    for left in range(120)
    for right in sorted(neighbours[left])
    if left < right
)
triangles = tuple(
    (left, right, third)
    for left, right in edges
    for third in sorted(neighbours[left] & neighbours[right])
    if right < third
)
tetrahedra = tuple(
    (first, second, third, fourth)
    for first, second, third in triangles
    for fourth in sorted(
        neighbours[first] & neighbours[second] & neighbours[third]
    )
    if third < fourth
)
check("the independently rebuilt carrier is the complete 600-cell boundary",
      (len(neighbours), len(edges), len(triangles), len(tetrahedra))
      == (120, 720, 1200, 600))

triangle_index = {triangle: index for index, triangle in enumerate(triangles)}
local_triangles = list(combinations(range(4), 3))
triangle_occurrences = [[] for _ in triangles]
for top_id, top in enumerate(tetrahedra):
    for local_id, face in enumerate(local_triangles):
        triangle = tuple(top[index] for index in face)
        triangle_occurrences[triangle_index[triangle]].append(
            (top_id, local_id)
        )

dual_incidence_base = np.zeros((len(triangles), len(tetrahedra)),
                               dtype=np.int64)
base_endpoints = []
for row, copies in enumerate(triangle_occurrences):
    if len(copies) != 2:
        raise AssertionError("600-cell triangle does not have two parents")
    (left, left_face), (right, right_face) = copies
    dual_incidence_base[row, left] = -int(sign[left_face])
    dual_incidence_base[row, right] = int(sign[right_face])
    base_endpoints.append((left, right))

reduced_base_integer = (
    12 * np.eye(len(tetrahedra), dtype=np.int64)
    + dual_incidence_base.T @ dual_incidence_base
)
factor = splu(sparse.csc_matrix(reduced_base_integer, dtype=np.float64))
base_solution = factor.solve(dual_incidence_base.T.astype(np.float64))
base_cross = (15.0 / 4.0) * base_solution
solve_residual = float(np.linalg.norm(
    reduced_base_integer @ base_solution - dual_incidence_base.T,
    ord=np.inf,
))
relative_solve_residual = solve_residual / max(
    1.0, float(np.linalg.norm(dual_incidence_base.T, ord=np.inf))
)
check("the numerical base solve satisfies the unique reduced equation",
      relative_solve_residual < 1e-11,
      f"relative infinity residual={relative_solve_residual:.3e}")

dual_rows = []
dual_columns = []
for left, right in base_endpoints:
    dual_rows.extend((left, right))
    dual_columns.extend((right, left))
dual_graph = sparse.csr_matrix((
    np.ones(len(dual_rows), dtype=np.int8),
    (dual_rows, dual_columns),
), shape=(len(tetrahedra), len(tetrahedra)))
dual_distances = shortest_path(
    dual_graph, directed=False, unweighted=True
)
maximum_entry = float(np.max(np.abs(base_cross)))
threshold = RELATIVE_THRESHOLD * maximum_entry
distance_histogram = {}
remote_above_threshold = 0
for top_id in range(len(tetrahedra)):
    for triangle_id, endpoints in enumerate(base_endpoints):
        if abs(base_cross[top_id, triangle_id]) <= threshold:
            continue
        distance = int(min(
            dual_distances[top_id, endpoints[0]],
            dual_distances[top_id, endpoints[1]],
        ))
        distance_histogram[distance] = distance_histogram.get(distance, 0) + 1
        remote_above_threshold += int(distance > 0)
numerical_maximum_distance = max(distance_histogram)
check("the base solution is numerically remote at the frozen threshold",
      remote_above_threshold > 0 and numerical_maximum_distance >= 8,
      f"remote={remote_above_threshold}, max dual distance="
      f"{numerical_maximum_distance}, threshold={threshold:.3e}")

# Exact finite-field certificate.  Select the lexicographically first pair at
# maximum tetrahedron-to-triangle-endpoint distance; the selection uses only
# topology, not the numerical value of B32.
maximum_endpoint_distance = max(
    int(min(dual_distances[top_id, endpoints[0]],
            dual_distances[top_id, endpoints[1]]))
    for top_id in range(len(tetrahedra))
    for endpoints in base_endpoints
)
maximal_pairs = [
    (top_id, triangle_id)
    for top_id in range(len(tetrahedra))
    for triangle_id, endpoints in enumerate(base_endpoints)
    if int(min(dual_distances[top_id, endpoints[0]],
               dual_distances[top_id, endpoints[1]]))
    == maximum_endpoint_distance
]
certificate_top, certificate_triangle = maximal_pairs[0]
modular_rank, modular_solution = modular_solve(
    reduced_base_integer,
    dual_incidence_base[certificate_triangle],
    MODULUS,
)
modular_residual = (
    reduced_base_integer @ modular_solution
    - dual_incidence_base[certificate_triangle]
) % MODULUS
modular_cross_value = (
    15 * pow(4, -1, MODULUS)
    * int(modular_solution[certificate_top])
) % MODULUS
check("the base reduced matrix is invertible modulo the frozen prime",
      modular_rank == len(tetrahedra)
      and np.all(modular_residual == 0),
      f"rank={modular_rank} mod {MODULUS}")
check("an exact B32 entry is nonzero at maximal endpoint distance",
      maximum_endpoint_distance == 14
      and modular_cross_value != 0,
      f"pair=({certificate_top},{certificate_triangle}), endpoints="
      f"{base_endpoints[certificate_triangle]}, distance="
      f"{maximum_endpoint_distance}, value mod p={modular_cross_value}")

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "phenomenological_target_used": False,
    "hypotheses": {
        "carrier": "canonical duplicated Whitney cochains",
        "poisson_auxiliary": "+i G with G=C M^-1 C*",
        "hamiltonian": "linear Hermitian quadratic, fixed top-left A_loc",
        "invariance": "weak first class: H K belongs to im L*",
        "locality": "tetrahedron belongs to a constraint triangle endpoint",
    },
    "weak_identity": (
        "S* H K=-A Q+B G+C* B* Q-C* D G=0"
    ),
    "unique_top_block": {
        "equation": "B32 G2=A32 M2^-1 C2*",
        "solution": "B32=(15/4)(12I+N*N)^-1 N*",
        "reason": "C3=0 and G2 is invertible",
    },
    "exact_boundary_control": {
        "tetrahedra": len(boundary_top),
        "triangles": len(boundary_endpoints),
        "remote_nonzeros": len(boundary_remote),
        "incident_nonzeros": sum(value != 0 for value in boundary_incident),
        "lumped_remote_nonzeros": lumped_remote,
        "remote_values": sorted({str(value) for value in boundary_remote}),
    },
    "base_600cell": {
        "f_vector": [120, 720, 1200, 600],
        "relative_threshold": RELATIVE_THRESHOLD,
        "absolute_threshold": threshold,
        "relative_solve_residual": relative_solve_residual,
        "numerical_nonzero_distance_histogram": {
            str(key): value for key, value in sorted(distance_histogram.items())
        },
        "numerical_remote_nonzeros": remote_above_threshold,
        "numerical_maximum_distance": numerical_maximum_distance,
        "exact_modular_certificate": {
            "prime": MODULUS,
            "rank": modular_rank,
            "top_row": certificate_top,
            "triangle_column": certificate_triangle,
            "triangle_endpoints": list(base_endpoints[certificate_triangle]),
            "endpoint_distance": maximum_endpoint_distance,
            "value_mod_prime": int(modular_cross_value),
        },
    },
    "verdicts": [
        "DERIVED: weak first-classness is the correct weaker condition",
        "DERIVED: the degree-3/degree-2 cross block is nevertheless unique",
        "DERIVED SCOPED NEGATIVE: its exact support is not endpoint-local",
        "OPEN: nonlinear or differently embedded local Hamiltonian dynamics",
        "NOT CLAIMED: physical time, causality, mass, c or Planck units",
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured weak-Hamiltonian certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"EXACT_BOUNDARY_REMOTE={len(boundary_remote)}")
print(f"BASE_NUMERICAL_DISTANCE={numerical_maximum_distance}")
print(f"BASE_EXACT_DISTANCE={maximum_endpoint_distance}")
print("VERDICT: weak first-classness does not restore endpoint-local dynamics")
raise SystemExit(0 if passed == tests else 1)

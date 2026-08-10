#!/usr/bin/env python3
"""Exact gate: do the currently certified actions generate the Hopf selector?

The protocol was frozen in commit e18c00b.  The chamber construction is
rebuilt directly so this audit does not inherit the expensive witness script's
conclusions by import.
"""

from itertools import combinations, product
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import sympy as sy


OUTPUT = Path(__file__).with_name("hopf_selector_action_gate.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def rotation_matrix(axis, angle):
    axis = axis/np.linalg.norm(axis)
    cross = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
    return (np.eye(3)+np.sin(angle)*cross
            +(1-np.cos(angle))*(cross@cross))


print("="*78)
print("HOPF SELECTOR: CERTIFIED-ACTION GATE")
print("="*78)

# ---------------------------------------------------------------------------
# Independent oriented-chamber reconstruction and A5 action.
# ---------------------------------------------------------------------------
phi = (1+np.sqrt(5.0))/2
vertices = []
for base in ((0, 1, phi), (1, phi, 0), (phi, 0, 1)):
    zero = base.index(0)
    nonzero = [index for index in range(3) if index != zero]
    for signs in product((-1, 1), repeat=2):
        vertex = list(base)
        for coordinate, sign in zip(nonzero, signs):
            vertex[coordinate] *= sign
        vertices.append(tuple(vertex))
vertices = np.asarray(sorted(set(vertices)))
edges = [(left, right) for left, right in combinations(range(12), 2)
         if abs(np.sum((vertices[left]-vertices[right])**2)-4) < 1e-8]
adjacency = [set() for _ in range(12)]
for left, right in edges:
    adjacency[left].add(right)
    adjacency[right].add(left)
faces = [(left, middle, right) for left, middle in edges
         for right in adjacency[left] & adjacency[middle]
         if middle < right]


def induced_vertex_permutation(matrix):
    moved = (matrix@vertices.T).T
    distances = ((moved[:, None, :]-vertices[None, :, :])**2).sum(axis=2)
    permutation = distances.argmin(axis=1)
    if distances[np.arange(12), permutation].max() >= 1e-12:
        return None
    return tuple(map(int, permutation))


rotation_candidates = []
for vertex in vertices:
    for multiple in range(1, 5):
        rotation_candidates.append(
            rotation_matrix(vertex, 2*np.pi*multiple/5))
for face in faces:
    centre = sum((vertices[index] for index in face), np.zeros(3))
    rotation_candidates.extend((rotation_matrix(centre, 2*np.pi/3),
                                rotation_matrix(centre, 4*np.pi/3)))
for edge in edges:
    rotation_candidates.append(
        rotation_matrix(vertices[edge[0]]+vertices[edge[1]], np.pi))
vertex_rotations = {tuple(range(12))}
for matrix in rotation_candidates:
    permutation = induced_vertex_permutation(matrix)
    if permutation is not None:
        vertex_rotations.add(permutation)

chambers = []
for face in faces:
    face_edges = [edge for edge in edges if set(edge).issubset(face)]
    for edge in face_edges:
        for vertex in edge:
            chambers.append((vertex, edge, face))
chambers = tuple(chambers)
chamber_index = {chamber: index for index, chamber in enumerate(chambers)}


def chamber_permutation(vertex_permutation):
    result = []
    for vertex, edge, face in chambers:
        image = (vertex_permutation[vertex],
                 tuple(sorted(vertex_permutation[index] for index in edge)),
                 tuple(sorted(vertex_permutation[index] for index in face)))
        result.append(chamber_index[image])
    return tuple(result)


chamber_rotations = tuple(chamber_permutation(permutation)
                          for permutation in vertex_rotations)
remaining = set(range(120))
orbits = []
while remaining:
    seed = min(remaining)
    orbit = {permutation[seed] for permutation in chamber_rotations}
    orbits.append(orbit)
    remaining -= orbit
reflection = chamber_permutation(induced_vertex_permutation(-np.eye(3)))
gamma = np.empty(120, dtype=np.int64)
gamma[list(orbits[0])] = 1
gamma[list(orbits[1])] = -1
plus = np.flatnonzero(gamma == 1)

chamber_edges = [(left, right) for left, right in combinations(range(120), 2)
                 if sum(chambers[left][coordinate]
                        == chambers[right][coordinate]
                        for coordinate in range(3)) == 2]
D = np.zeros((120, 120), dtype=np.int64)
for left, right in chamber_edges:
    D[left, right] = D[right, left] = 1
J = np.zeros((120, 120), dtype=np.int64)
J[np.asarray(reflection), np.arange(120)] = 1
Gamma = np.diag(gamma)
check("independent carrier is the fixed A5 chamber KO6 geometry",
      len(vertices) == 12 and len(edges) == 30 and len(faces) == 20
      and len(chambers) == 120 and len(vertex_rotations) == 60
      and sorted(map(len, orbits)) == [60, 60]
      and len(chamber_edges) == 180
      and np.array_equal(D.sum(axis=1), 3*np.ones(120, dtype=np.int64))
      and np.array_equal(J@D, D@J)
      and np.array_equal(J@Gamma, -Gamma@J))

# ---------------------------------------------------------------------------
# The committed B1 colouring and its exact A5 stabilizer.
# ---------------------------------------------------------------------------
CELL_LABELS = (
    2, 0, 2, 2, 1, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1,
    0, 3, 3, 3, 3, 1, 3, 2, 2, 3, 0, 3, 3, 0, 3, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 3, 3, 3, 3, 2,
    3, 1, 2, 1, 1, 3, 1, 1, 2, 3, 1, 1, 3,
)
labels = np.asarray(CELL_LABELS, dtype=np.int64)
label_by_chamber = {int(chamber): int(labels[index])
                    for index, chamber in enumerate(plus)}
stabilizer = [permutation for permutation in chamber_rotations
              if all(label_by_chamber[permutation[int(chamber)]]
                     == label_by_chamber[int(chamber)] for chamber in plus)]
colouring_orbit = {
    tuple(label_by_chamber[permutation[int(chamber)]] for chamber in plus)
    for permutation in chamber_rotations
}
capacities = tuple(int(np.count_nonzero(labels == index))
                   for index in range(4))
check("the B1 central-support colouring has trivial A5 stabilizer",
      capacities == (4, 25, 12, 19)
      and len(stabilizer) == 1 and len(colouring_orbit) == 60,
      f"capacities={capacities}, stabilizer={len(stabilizer)}, "
      f"orbit={len(colouring_orbit)}")

# ---------------------------------------------------------------------------
# Rebuild the representation and isolate its only M2-to-scalar fluctuation
# block.  Every calculation is integral after the chamber permutations.
# ---------------------------------------------------------------------------
factor_sizes = (2, 1, 1, 1)
cells = ((0, 1, 2), (1, 2, 25), (3, 1, 12), (2, 3, 19))
cell_vertices = tuple(plus[np.flatnonzero(labels == index)]
                      for index in range(4))


def representation(matrix2, lambda1, lambda2, lambda3):
    matrix2 = np.asarray(matrix2, dtype=np.complex128)
    scalars = (None, lambda1, lambda2, lambda3)
    answer = np.zeros((120, 120), dtype=np.complex128)
    for cell_index, (left, right, multiplicity) in enumerate(cells):
        occupied = cell_vertices[cell_index]
        if left == 0:
            plus_block = np.kron(matrix2, np.eye(multiplicity))
        else:
            plus_block = scalars[left]*np.eye(len(occupied))
        answer[np.ix_(occupied, occupied)] = plus_block
        reflected = np.asarray(reflection)[occupied]
        if right == 0:
            minus_block = np.kron(np.eye(factor_sizes[left]),
                                  np.kron(matrix2, np.eye(multiplicity)))
        else:
            minus_block = scalars[right]*np.eye(len(occupied))
        answer[np.ix_(reflected, reflected)] = minus_block
    return answer


zero2 = np.zeros((2, 2), dtype=np.int64)
matrix_units = []
for row in range(2):
    for column in range(2):
        unit = zero2.copy()
        unit[row, column] = 1
        matrix_units.append(unit)
algebra_basis = [representation(unit, 0, 0, 0)
                 for unit in matrix_units]
algebra_basis.extend((representation(zero2, 1, 0, 0),
                      representation(zero2, 0, 1, 0),
                      representation(zero2, 0, 0, 1)))

block_rows = cell_vertices[0]
block_columns = np.asarray(reflection)[cell_vertices[1]]
D_block = D[np.ix_(block_rows, block_columns)]
one_form_blocks = []
for left in algebra_basis:
    for right in algebra_basis:
        one_form = left@(D@right-right@D)
        block = np.rint(one_form[np.ix_(block_rows, block_columns)].real)
        one_form_blocks.append(block.astype(np.int64))
one_form_columns = np.column_stack(
    [block.reshape(-1) for block in one_form_blocks]
)
one_form_gram = one_form_columns.T@one_form_columns
one_form_rank = sy.Matrix(one_form_gram.tolist()).rank()

m2_images = []
for unit in matrix_units:
    m2_images.append(np.kron(unit, np.eye(2, dtype=np.int64))@D_block)
m2_columns = np.column_stack([image.reshape(-1) for image in m2_images])
m2_rank = sy.Matrix((m2_columns.T@m2_columns).tolist()).rank()
combined = np.column_stack((one_form_columns, m2_columns))
combined_rank = sy.Matrix((combined.T@combined).tolist()).rank()
check("the noncommutative Dirac block is full-row-rank with twelve incidences",
      D_block.shape == (4, 25)
      and sy.Matrix(D_block.tolist()).rank() == 4
      and np.count_nonzero(D_block) == 12)
check("its represented one-forms are exactly (M2 tensor I2)D_block",
      one_form_rank == m2_rank == combined_rank == 4,
      f"one-form rank={one_form_rank}, M2 image rank={m2_rank}, "
      f"union rank={combined_rank}")

# Coefficient space M2 under left multiplication.  In the matrix-unit basis,
# the two column subspaces are two fundamental doublets.  The full commutant
# is right M2, so geometry has not selected one multiplicity column.
def action_matrix(unit, side):
    columns = []
    for basis_element in matrix_units:
        image = unit@basis_element if side == "left" else basis_element@unit
        columns.append(image.reshape(-1))
    return sy.Matrix(np.column_stack(columns).tolist())


left_actions = [action_matrix(unit, "left") for unit in matrix_units]
right_actions = [action_matrix(unit, "right") for unit in matrix_units]
commuting_lr = all(left*right == right*left
                   for left in left_actions for right in right_actions)
right_flat = sy.Matrix.hstack(*(matrix.reshape(16, 1)
                                for matrix in right_actions))

unknowns = sy.symbols("c0:16")
candidate = sy.Matrix(4, 4, unknowns)
constraints = []
for left in left_actions:
    constraints.extend(list(candidate*left-left*candidate))
constraint_matrix, _ = sy.linear_eq_to_matrix(constraints, unknowns)
commutant_dimension = 16-constraint_matrix.rank()

column_one = sy.diag(1, 0, 1, 0)
column_two = sy.eye(4)-column_one
doublet_split = (
    column_one.rank() == column_two.rank() == 2
    and all(column_one*left == left*column_one
            and column_two*left == left*column_two
            for left in left_actions)
)
check("the M2 coefficient module is two left-SU2 fundamental columns",
      doublet_split and commuting_lr
      and right_flat.rank() == commutant_dimension == 4,
      f"left commutant dimension={commutant_dimension} (right M2)")

# ---------------------------------------------------------------------------
# Degree ceiling of the actually certified D^2/D^4 moment truncation.
# ---------------------------------------------------------------------------
n1, n2, n3 = sy.symbols("n1 n2 n3")
linear_family = sy.Matrix([[1+n1, n2+n3],
                           [n2-n3, 2-n1]])
moment2 = sy.expand(sy.trace(linear_family**2))
moment4 = sy.expand(sy.trace(linear_family**4))
moment6 = sy.expand(sy.trace(linear_family**6))
degree2 = sy.Poly(moment2, n1, n2, n3).total_degree()
degree4 = sy.Poly(moment4, n1, n2, n3).total_degree()
degree6 = sy.Poly(moment6, n1, n2, n3).total_degree()
check("linear fluctuations have degree ceilings 2/4/6 in the three moments",
      (degree2, degree4, degree6) == (2, 4, 6),
      f"total degrees={(degree2, degree4, degree6)}")

spectral_source = Path(__file__).with_name("verify_spectral_action.py").read_text()
certified_moment_snippets = (
    "c0 = N_total",
    "c1 = np.sum(all_evals_D2)",
    "c2 = 0.5 * np.sum(all_evals_D2**2)",
)
check("the registered finite spectral-action file stops at unfluctuated D^4",
      all(snippet in spectral_source for snippet in certified_moment_snippets)
      and "all_evals_D2**3" not in spectral_source
      and "D_A" not in spectral_source,
      "certified powers are D^0,D^2,D^4; no D_A or D^6 moment")

payload = {
    "protocol_commit": "e18c00b",
    "chamber": {
        "A5_order": len(chamber_rotations),
        "B1_capacities": list(capacities),
        "B1_stabilizer_order": len(stabilizer),
        "B1_colouring_orbit_size": len(colouring_orbit),
    },
    "noncommutative_one_form_block": {
        "shape": list(D_block.shape),
        "incidences": int(np.count_nonzero(D_block)),
        "D_block_rank": int(sy.Matrix(D_block.tolist()).rank()),
        "complex_one_form_dimension": int(one_form_rank),
        "module": "M2(C) = two left-SU2 fundamental columns",
        "left_action_commutant": "right M2(C)",
        "commutant_complex_dimension": int(commutant_dimension),
    },
    "degree_gate": {
        "linear_fluctuation_moment_degrees": [degree2, degree4, degree6],
        "currently_certified_moment_powers": [0, 2, 4],
        "sixth_order_available_in_current_file": False,
    },
    "verdict": {
        "geometry_selected_A5_order_parameter": False,
        "reason": "B1 central supports have trivial A5 stabilizer",
        "current_action_generates_S6_with_required_sign": False,
        "status": "DERIVED NEGATIVE for current certified construction",
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("action-gate JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: the valid chamber witness contains two M2 doublets, not one.")
print("DERIVED NEGATIVE: its embedding has trivial A5 stabilizer (orbit 60).")
print("DERIVED NEGATIVE: the certified unfluctuated D^0/D^2/D^4 moments cannot")
print("                  contain the homogeneous sixth-order Hopf anisotropy.")
print("OPEN: an A5-equivariant selected triple, D_A^6 coefficient and its sign.")
raise SystemExit(0 if passed == tests else 1)

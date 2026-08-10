#!/usr/bin/env python3
"""Exact audit of superselection from the already-defined Hopf dynamics.

The mechanisms and decision boundary were frozen in commit 0807c44 before
their fixed algebras were computed.  All load-bearing ranks, commutators and
spectral multiplicities are exact.
"""

from collections import Counter
import json
from pathlib import Path

import numpy as np
import sympy as sp

from verify_hopf_fibration_invariants import (
    build_2I,
    build_adjacency,
    build_fiber_adjacency,
    find_all_hopf_fibrations,
    find_vertex_index,
    quat_mult,
)


OUTPUT = Path(__file__).with_name("hopf_dynamical_superselection.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def trace_product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = result @ matrix
    return int(np.trace(result))


print("="*78)
print("EXACT EXISTING-DYNAMICS HOPF SUPERSELECTION AUDIT")
print("="*78)

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
fiber_adjacencies = [
    np.rint(build_fiber_adjacency(adjacency, fibration)).astype(np.int64)
    for fibration in fibrations
]
boxes = [6*fiber-adjacency for fiber in fiber_adjacencies]
basis = [boxes[index]-boxes[5] for index in range(5)]
gram = sp.Matrix([
    [trace_product(basis[row], basis[col]) for col in range(5)]
    for row in range(5)
])
check("the fixed geometry reconstructs the exact five-dimensional Box space",
      len(boxes) == 6 and np.array_equal(sum(boxes), np.zeros_like(adjacency))
      and gram.rank() == 5)


def conjugation_permutation(group_element):
    group_inverse = group_element.copy()
    group_inverse[1:] *= -1
    return tuple(find_vertex_index(
        vertices,
        quat_mult(quat_mult(group_element, vertex), group_inverse),
    ) for vertex in vertices)


vertex_actions = sorted(set(conjugation_permutation(element)
                            for element in vertices))
fibration_by_signature = {
    tuple(sorted(tuple(sorted(fiber)) for fiber in fibration)): index
    for index, fibration in enumerate(fibrations)
}
label_actions = []
for vertex_action in vertex_actions:
    label_action = []
    for fibration in fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fiber))
            for fiber in fibration
        ))
        label_action.append(fibration_by_signature[signature])
    label_actions.append(tuple(label_action))
label_actions = sorted(set(label_actions))

n = 6
permutation_matrices = []
for action in label_actions:
    permutation = sp.zeros(n)
    for source, target in enumerate(action):
        permutation[target, source] = 1
    permutation_matrices.append(permutation)
check("the derived effective action has 60 elements and is transitive",
      len(label_actions) == 60
      and {action[0] for action in label_actions} == set(range(n)))

# Exact Reynolds superoperator for conjugation on M6.  With column or row
# vectorization the permutation matrices are orthogonal, so kron(P,P) gives
# the same exact rank and fixed space.
twirl_numerator = sum((sp.kronecker_product(P, P)
                       for P in permutation_matrices), sp.zeros(n*n))
twirl = twirl_numerator/sp.Integer(len(permutation_matrices))
identity = sp.eye(n)
all_ones = sp.ones(n)
check("the exact A5 conjugation twirl has rank two",
      twirl.rank() == 2 and twirl*twirl == twirl,
      "fixed algebra dimension=2, not dim C(F)=6")


def vectorize(matrix):
    return sp.Matrix(list(matrix))


check("the twirl image is exactly span{I,J}",
      twirl*vectorize(identity) == vectorize(identity)
      and twirl*vectorize(all_ones) == vectorize(all_ones)
      and sp.Matrix.hstack(vectorize(identity), vectorize(all_ones)).rank()
      == twirl.rank())

label_projection = sp.zeros(n)
label_projection[0, 0] = 1
off_diagonal_unit = sp.zeros(n)
off_diagonal_unit[0, 1] = 1
twirled_label = sp.Matrix(twirl*vectorize(label_projection)).reshape(n, n)
twirled_transition = sp.Matrix(
    twirl*vectorize(off_diagonal_unit)
).reshape(n, n)
expected_transition = (all_ones-identity)/sp.Integer(30)
check("no individual label projection is fixed by the symmetry twirl",
      twirled_label == identity/sp.Integer(6)
      and twirled_label != label_projection)
check("symmetry averaging retains an off-diagonal transition channel",
      twirled_transition == expected_transition
      and twirled_transition != sp.zeros(n),
      "E_01 averages to (J-I)/30")

# The permutation average on H_F is projection onto the constant vector.
# Hence its restriction to the zero-sum Box module W is zero.  Equivariance
# and linearity of H_X then make the simultaneous orbit average of H_X zero.
representation_average = sum(permutation_matrices, sp.zeros(n))/sp.Integer(60)
check("the group average kills every nontrivial Box order-parameter direction",
      representation_average == all_ones/sp.Integer(6)
      and all(representation_average*sp.Matrix([
          1 if label == coordinate else (-1 if label == 5 else 0)
          for label in range(6)
      ]) == sp.zeros(6, 1) for coordinate in range(5)),
      "the only invariant label direction is the excluded constant mode")

# Complete Hessian H_X of Tr(X^3) in the six synthesis directions.
box_products = [[boxes[row]@boxes[col] for col in range(n)]
                for row in range(n)]


def label_hessian(X):
    return sp.Matrix(n, n, lambda row, col:
                     3*(int(np.sum(X*box_products[row][col].T))
                        + int(np.sum(X*box_products[col][row].T))))


basis_hessians = [label_hessian(direction) for direction in basis]
vertex_hessians = [label_hessian(box) for box in boxes]
off_diagonal_pairs = [(row, col) for row in range(n)
                      for col in range(row+1, n)]
off_diagonal_map = sp.Matrix([
    [basis_hessians[coordinate][row, col] for coordinate in range(5)]
    for row, col in off_diagonal_pairs
])
check("the off-diagonal Hessian map is injective on W",
      off_diagonal_map.rank() == 5 and len(off_diagonal_map.nullspace()) == 0,
      "only X=0 gives a label-diagonal full Hessian")
check("no normalized X can fix all six labels under Hessian evolution",
      off_diagonal_map.rank() == gram.rank() == 5,
      "Tr(X^2)=7200 excludes the unique diagonal case X=0")

# At each desired Box point, compute all label-projection commutators and the
# graph of nonzero off-diagonal Hessian entries.
all_commutators_nonzero = True
graph_component_counts = []
diagonal_commutant_dimensions = []
fixed_algebra_dimensions = []
eigenvalue_multiplicity_patterns = []
for hessian in vertex_hessians:
    for label in range(n):
        projection = sp.zeros(n)
        projection[label, label] = 1
        commutator = hessian*projection-projection*hessian
        if commutator == sp.zeros(n):
            all_commutators_nonzero = False

    unseen = set(range(n))
    components = []
    while unseen:
        component = {unseen.pop()}
        frontier = list(component)
        while frontier:
            row = frontier.pop()
            neighbors = {
                col for col in range(n)
                if col != row and hessian[row, col] != 0
            }
            new = neighbors & unseen
            unseen -= new
            component |= new
            frontier.extend(new)
        components.append(component)
    graph_component_counts.append(len(components))
    # A diagonal matrix commutes with H iff its entries are constant on each
    # connected component of the nonzero off-diagonal graph.
    diagonal_commutant_dimensions.append(len(components))

    multiplicities = sorted(hessian.eigenvals().values())
    eigenvalue_multiplicity_patterns.append(multiplicities)
    fixed_algebra_dimensions.append(sum(value*value
                                        for value in multiplicities))

check("all 36 desired-point label projections fail to commute with H_Box",
      all_commutators_nonzero)
check("the nonzero Hessian graph is connected at every desired Box point",
      graph_component_counts == [1]*6)
check("C(F) intersect {H_Box}' consists only of scalars",
      diagonal_commutant_dimensions == [1]*6,
      "time dephasing preserves no individual fibration projection")
check("the full time-dephasing fixed algebra has dimension ten, not six",
      fixed_algebra_dimensions == [10]*6
      and eigenvalue_multiplicity_patterns == [[1, 1, 2, 2]]*6,
      "eigenvalue multiplicities 1,1,2,2 give sum m^2=10")

payload = {
    "protocol_commit": "0807c44",
    "symmetry_twirl": {
        "group": "A5",
        "superoperator_rank": int(twirl.rank()),
        "fixed_algebra": "span{I,J} ~= C+C",
        "fixed_algebra_dimension": int(twirl.rank()),
        "label_algebra_dimension": 6,
        "individual_label_projections_fixed": False,
        "twirl_E00": "I/6",
        "twirl_E01": "(J-I)/30",
        "nonzero_order_parameter_survives_orbit_average": False,
    },
    "Hessian_dephasing": {
        "off_diagonal_map_rank_on_W": int(off_diagonal_map.rank()),
        "dimension_W": int(gram.rank()),
        "nonzero_diagonal_Hessian_exists": False,
        "desired_Box_points_checked": 6,
        "label_projection_commutators_nonzero": 36,
        "off_diagonal_graph_components": graph_component_counts,
        "diagonal_fixed_algebra_dimension": diagonal_commutant_dimensions,
        "full_fixed_algebra_dimension": fixed_algebra_dimensions,
        "eigenvalue_multiplicities": eigenvalue_multiplicity_patterns,
    },
    "verdict": (
        "DERIVED NEGATIVE: neither existing parameter-free averaging "
        "mechanism produces the six-label diagonal algebra.  The A5 twirl "
        "fixes span{I,J} and retains a collective off-diagonal channel.  The "
        "full cubic Hessian is nondiagonal for every normalized X; at all six "
        "Box points its intersection with C(F) is only the scalars."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact dynamical-superselection audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED NEGATIVE: A5 twirl fixes C+C, not the six-label C(F).")
print("DERIVED NEGATIVE: every normalized full Hessian mixes label sectors.")
print("DERIVED NEGATIVE: at each Box point only scalar diagonal observables")
print("                  survive Hessian-generated time dephasing.")
raise SystemExit(0 if passed == tests else 1)

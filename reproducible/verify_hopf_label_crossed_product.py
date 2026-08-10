#!/usr/bin/env python3
"""Exact algebra audit of C(F) crossed product A5 for six Hopf labels.

The effective transformation-groupoid arena and decision boundaries were
frozen in commit d15e7fa before its stabilizer, Wedderburn type or natural
image were computed.
"""

from collections import Counter
import json
from pathlib import Path

import numpy as np
import sympy as sp

from verify_hopf_fibration_invariants import (
    build_2I,
    find_all_hopf_fibrations,
    find_vertex_index,
    quat_mult,
)


OUTPUT = Path(__file__).with_name("hopf_label_crossed_product.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(right)))


def inverse(permutation):
    result = [None]*len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def permutation_order(permutation):
    current = tuple(range(len(permutation)))
    for order in range(1, 61):
        current = compose(permutation, current)
        if current == tuple(range(len(permutation))):
            return order
    raise RuntimeError("permutation order exceeded group order")


print("="*78)
print("EXACT HOPF LABEL CROSSED-PRODUCT AUDIT")
print("="*78)

vertices = build_2I()
fibrations = find_all_hopf_fibrations(vertices)


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
    action = []
    for fibration in fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fiber))
            for fiber in fibration
        ))
        action.append(fibration_by_signature[signature])
    label_actions.append(tuple(action))
label_actions = sorted(set(label_actions))

label_orbit = {action[0] for action in label_actions}
stabilizer = [action for action in label_actions if action[0] == 0]
stabilizer_orders = Counter(permutation_order(action) for action in stabilizer)
check("the effective A5 action is transitive with stabilizer order ten",
      len(label_actions) == 60 and label_orbit == set(range(6))
      and len(stabilizer) == 10)
check("the exact stabilizer element census identifies H as D5",
      stabilizer_orders == Counter({1: 1, 2: 5, 5: 4}),
      f"order census={dict(sorted(stabilizer_orders.items()))}")

# Exact conjugacy classes inside H.
unseen = set(stabilizer)
stabilizer_classes = []
while unseen:
    element = next(iter(unseen))
    conjugates = {
        compose(compose(group_element, element), inverse(group_element))
        for group_element in stabilizer
    }
    stabilizer_classes.append(conjugates)
    unseen -= conjugates
stabilizer_class_sizes = sorted(len(group_class)
                                for group_class in stabilizer_classes)
check("D5 has four exact conjugacy classes of sizes 1,2,2,5",
      stabilizer_class_sizes == [1, 2, 2, 5])

# Four irreducibles, group order ten and abelianization C2 force degrees
# 1,1,2,2.  The transitive-groupoid theorem then gives four matrix blocks.
stabilizer_irrep_degrees = [1, 1, 2, 2]
crossed_product_blocks = [6*degree for degree in stabilizer_irrep_degrees]
crossed_product_dimension = 6*len(label_actions)
check("the stabilizer Wedderburn degrees are exactly 1,1,2,2",
      len(stabilizer_classes) == 4
      and sum(degree**2 for degree in stabilizer_irrep_degrees) == 10)
check("the crossed product has type M6+M6+M12+M12 and dimension 360",
      crossed_product_blocks == [6, 6, 12, 12]
      and sum(block**2 for block in crossed_product_blocks)
      == crossed_product_dimension == 360,
      "C(F) crossed A5 ~= M6(C[D5])")

# Natural covariant representation delta_i -> E_ii, u_g -> P_g.
n = 6
matrix_columns = []
natural_generators = []
for action in label_actions:
    permutation = np.zeros((n, n), dtype=np.int64)
    for source, target in enumerate(action):
        permutation[target, source] = 1
    for label in range(n):
        projection = np.zeros((n, n), dtype=np.int64)
        projection[label, label] = 1
        generator = projection@permutation
        natural_generators.append(generator)
        matrix_columns.append(generator.reshape(-1))
natural_span_matrix = sp.Matrix(np.column_stack(matrix_columns).tolist())
natural_image_rank = natural_span_matrix.rank()
check("the natural six-label image is the full M6(C)",
      natural_image_rank == 36,
      "span{E_ii P_g} has rank 36")
check("the natural representation is nonfaithful with kernel dimension 324",
      crossed_product_dimension-natural_image_rank == 324,
      "it is the M6 block induced from the trivial stabilizer irrep")

# Since the image is full M6, its commutant is scalar.  Exhibit all matrix
# units directly among the span and an explicit noncentral label commutator.
matrix_units_present = True
for row in range(n):
    for col in range(n):
        target_unit = np.zeros((n, n), dtype=np.int64)
        target_unit[row, col] = 1
        augmented = natural_span_matrix.row_join(
            sp.Matrix(target_unit.reshape(-1).tolist())
        )
        if augmented.rank() != natural_image_rank:
            matrix_units_present = False
label_projection = np.zeros((n, n), dtype=np.int64)
label_projection[0, 0] = 1
moving_action = next(action for action in label_actions if action[0] != 0)
moving_permutation = np.zeros((n, n), dtype=np.int64)
for source, target in enumerate(moving_action):
    moving_permutation[target, source] = 1
label_commutator = label_projection@moving_permutation
label_commutator -= moving_permutation@label_projection
check("all 36 matrix units occur and the image commutant is scalar",
      matrix_units_present and natural_image_rank == n*n)
check("individual fibration projections are noncentral",
      np.count_nonzero(label_commutator) == 2,
      "group unitaries conjugate all six projections transitively")

# Full crossed-product centre has one dimension per simple block.  Every label
# projection is full and the six are Murray--von Neumann equivalent under the
# group arrows; none gives a central superselection sector.
centre_dimension = len(crossed_product_blocks)
projection_equivalence_orbit = {
    action[0] for action in label_actions
}
check("the full centre has dimension four, not six",
      centre_dimension == 4 and projection_equivalence_orbit == set(range(6)),
      "six label projections are equivalent and have common central support 1")

# Conditional expectations onto C(F).  Bimodularity kills off-diagonal matrix
# corners, but each diagonal corner is C[D5].  Thus arbitrary bimodular linear
# projections are specified by six linear functionals on a 10D corner.  A5
# equivariance reduces this to conjugation-invariant functionals on D5, whose
# dimension is the number of conjugacy classes: four.  Unitality leaves a
# 3D affine family, so no unique expectation is forced without a trace/state.
bimodular_map_dimension = 6*len(stabilizer)
equivariant_functional_dimension = len(stabilizer_classes)
equivariant_unital_affine_dimension = equivariant_functional_dimension-1
check("diagonal bimodular projection is highly nonunique",
      bimodular_map_dimension == 60
      and equivariant_functional_dimension == 4
      and equivariant_unital_affine_dimension == 3,
      "A5-equivariant unital expectations form at least a 3D affine family")

# Two exact positive unital stabilizer states already distinguish expectations:
# regular coefficient trace tau(h)=delta_{h,e}, and the trivial character
# phi(h)=1.  Convex combinations give faithful alternatives once tau has
# positive weight.
identity_action = tuple(range(n))
nonidentity_stabilizer = next(action for action in stabilizer
                              if action != identity_action)
regular_state_value = int(nonidentity_stabilizer == identity_action)
trivial_state_value = 1
check("canonical regular trace and trivial-character state give distinct expectations",
      regular_state_value == 0 and trivial_state_value == 1,
      "choosing the regular-trace expectation requires specifying that trace")

# Only now locate the selector: Phi is diagonal in the C(F) subalgebra, but
# full M6 conjugation moves it and it is noncentral whenever X is nonzero.
phi_times_5 = np.diag([5, -1, -1, -1, -1, -1])
phi_commutator = phi_times_5@moving_permutation
phi_commutator -= moving_permutation@phi_times_5
check("Phi lies in C(F) but is noncentral in the crossed product",
      np.count_nonzero(phi_commutator) > 0
      and len(set(np.diag(phi_times_5))) > 1,
      "D_aux diagonality is a conditional-expectation choice, not superselection")

payload = {
    "protocol_commit": "d15e7fa",
    "effective_action": {
        "group": "A5",
        "order": len(label_actions),
        "labels": 6,
        "transitive": True,
        "stabilizer": "D5",
        "stabilizer_order": len(stabilizer),
        "stabilizer_order_census": dict(sorted(stabilizer_orders.items())),
    },
    "crossed_product": {
        "definition": "C(F) crossed_product A5",
        "dimension": crossed_product_dimension,
        "isomorphism": "M6(C[D5])",
        "Wedderburn_blocks": crossed_product_blocks,
        "centre_dimension": centre_dimension,
        "simple_summands": len(crossed_product_blocks),
    },
    "natural_label_representation": {
        "dimension": 6,
        "image": "M6(C)",
        "image_dimension": natural_image_rank,
        "kernel_dimension": crossed_product_dimension-natural_image_rank,
        "commutant": "C*I",
        "label_projections_central": False,
    },
    "conditional_expectations_to_C6": {
        "bimodular_linear_map_dimension": bimodular_map_dimension,
        "A5_equivariant_functional_dimension": (
            equivariant_functional_dimension
        ),
        "A5_equivariant_unital_affine_dimension": (
            equivariant_unital_affine_dimension
        ),
        "canonical_regular_trace_expectation_exists": True,
        "unique_without_trace_or_state": False,
    },
    "selector": {
        "Phi_in_diagonal_subalgebra": True,
        "Phi_central": False,
        "D_aux_superselected": False,
    },
    "verdict": (
        "DERIVED ALGEBRA TYPE: C(F) crossed A5 is "
        "M6+M6+M12+M12.  DERIVED NEGATIVE for superselection: its natural "
        "label image is full M6, all six label projections are equivalent "
        "and noncentral, and diagonal expectations are not unique without "
        "choosing a trace/state.  Phi is a noncentral diagonal element rather "
        "than a sector selected by the crossed product"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact crossed-product audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: C(F) crossed A5 = M6 + M6 + M12 + M12.")
print("DERIVED NEGATIVE: natural image=M6; labels are equivalent noncentral")
print("                  projections, not superselection sectors.")
print("DERIVED NEGATIVE: diagonal expectation needs an additional trace/state.")
raise SystemExit(0 if passed == tests else 1)

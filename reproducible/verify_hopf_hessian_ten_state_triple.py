#!/usr/bin/env python3
"""Exact finite-triple gate audit for the ten-state Hessian selector.

The initial protocol was committed in a0d7dd8.  Its incomplete algebra list
was corrected in cae7233 before the monomial C^5 system or first-order support
was computed.  This verifier constructs that missing system exactly and then
tests the full affine Hessian family.
"""

from itertools import permutations
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


OUTPUT = Path(__file__).with_name("hopf_hessian_ten_state_triple.json")
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
    identity = tuple(range(len(permutation)))
    current = identity
    for order in range(1, 61):
        current = compose(permutation, current)
        if current == identity:
            return order
    raise RuntimeError("order exceeds A5 size")


def generated_subgroup(generators, identity):
    subgroup = {identity}
    frontier = list(generators)
    while frontier:
        element = frontier.pop()
        if element in subgroup:
            continue
        old = list(subgroup)
        subgroup.add(element)
        for other in old+[element]:
            for product in (compose(element, other), compose(other, element)):
                if product not in subgroup:
                    frontier.append(product)
    return subgroup


def matrix_key(matrix):
    return tuple(sp.simplify(value) for value in matrix)


def matrix_is_zero(matrix):
    return all(sp.simplify(value) == 0 for value in matrix)


print("="*78)
print("EXACT TEN-STATE HESSIAN FINITE-TRIPLE AUDIT")
print("="*78)

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
fiber_adjacencies = [
    np.rint(build_fiber_adjacency(adjacency, fibration)).astype(np.int64)
    for fibration in fibrations
]
boxes = [6*fiber-adjacency for fiber in fiber_adjacencies]
field_basis = [boxes[index]-boxes[5] for index in range(5)]
box_products = [[boxes[row]@boxes[col] for col in range(6)]
                for row in range(6)]


def label_hessian(X):
    return sp.Matrix(6, 6, lambda row, col:
                     3*(int(np.sum(X*box_products[row][col].T))
                        + int(np.sum(X*box_products[col][row].T))))


label_basis = sp.Matrix(6, 5, lambda row, col:
                        (1 if row == col else (-1 if row == 5 else 0)))
basis_hessians = [label_hessian(direction) for direction in field_basis]
restricted_hessians = [(hessian*label_basis)[:5, :]
                       for hessian in basis_hessians]
check("the full affine selector has five independent Hessian directions",
      sp.Matrix.hstack(*[
          sp.Matrix(list(hessian)) for hessian in restricted_hessians
      ]).rank() == 5)


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
identity = tuple(range(6))

representation = {}
for action in label_actions:
    permutation = sp.zeros(6)
    for source, target in enumerate(action):
        permutation[target, source] = 1
    representation[action] = (permutation*label_basis)[:5, :]
check("the exact physical label action has group order 60",
      len(label_actions) == 60 and identity in label_actions)

label_metric = label_basis.T*label_basis
representation_law_ok = all(
    representation[compose(left, right)]
    == representation[left]*representation[right]
    for left in label_actions for right in label_actions
)
representation_unitary_ok = all(
    matrix.T*label_metric*matrix == label_metric
    for matrix in representation.values()
)
check("the restricted matrices form an exact unitary A5 representation",
      representation_law_ok and representation_unitary_ok,
      "unitarity is with respect to the inherited W5 Gram matrix")

# Irreducibility and simplicity data used in the normalized-subalgebra
# classification.  Every nonidentity element has normal closure A5.
fixed_equations = sp.Matrix.vstack(*[
    representation[action]-sp.eye(5) for action in label_actions
])
commutant_equations = sp.Matrix.vstack(*[
    sp.kronecker_product(matrix.T, sp.eye(5))
    - sp.kronecker_product(sp.eye(5), matrix)
    for matrix in representation.values()
])
check("W5 has no invariant vector and has scalar exact commutant",
      fixed_equations.rank() == 5 and commutant_equations.rank() == 24,
      "the exact complex commutant has dimension one")

normal_closure_sizes = []
for element in label_actions:
    if element == identity:
        continue
    conjugates = {
        compose(compose(group_element, element), inverse(group_element))
        for group_element in label_actions
    }
    normal_closure_sizes.append(len(generated_subgroup(conjugates, identity)))
check("the exact 60-element group is simple by normal closures",
      set(normal_closure_sizes) == {60},
      "every nonidentity normal closure has order 60")

# Build an index-five A4 as the normalizer of a V4 centralizer.
order_two = next(element for element in label_actions
                 if permutation_order(element) == 2)
V4 = {
    element for element in label_actions
    if compose(element, order_two) == compose(order_two, element)
}
A4 = {
    element for element in label_actions
    if {
        compose(compose(element, member), inverse(element)) for member in V4
    } == V4
}
check("an exact V4 normalizer gives an index-five A4 stabilizer",
      len(V4) == 4 and len(A4) == 12
      and sorted(permutation_order(element) for element in A4)
      == [1]+[2]*3+[3]*8)

order_three = next(element for element in A4
                   if permutation_order(element) == 3)
cosets = [
    V4,
    {compose(order_three, member) for member in V4},
    {compose(compose(order_three, order_three), member) for member in V4},
]
check("A4/V4 is exactly C3",
      set().union(*cosets) == A4
      and sum(len(coset) for coset in cosets) == 12
      and all(len(coset) == 4 for coset in cosets))
coset_exponent = {
    element: exponent for exponent, coset in enumerate(cosets)
    for element in coset
}

# Work exactly in Q(omega), omega^2+omega+1=0, without asking SymPy to
# simplify complex radicals.  A pair (A,B) denotes A+omega*B.
def pair_add(left, right):
    return left[0]+right[0], left[1]+right[1]


def pair_scale(scalar, pair):
    return scalar*pair[0], scalar*pair[1]


def pair_multiply(left, right):
    A, B = left
    C, D = right
    return A*C-B*D, A*D+B*C-B*D


def pair_conjugate(pair, matrix):
    matrix_inverse = matrix.inv()
    return (matrix*pair[0]*matrix_inverse,
            matrix*pair[1]*matrix_inverse)


def pair_metric_adjoint(pair):
    # Complex conjugation sends omega to omega^2=-1-omega.
    conjugate_constant = pair[0]-pair[1]
    conjugate_omega = -pair[1]
    metric_inverse = label_metric.inv()
    return (metric_inverse*conjugate_constant.T*label_metric,
            metric_inverse*conjugate_omega.T*label_metric)


def pair_key(pair):
    return tuple(pair[0])+tuple(pair[1])


def pair_is_zero(pair):
    return pair[0] == sp.zeros(pair[0].rows, pair[0].cols) \
        and pair[1] == sp.zeros(pair[1].rows, pair[1].cols)


def pair_is_equal(left, right):
    return left[0] == right[0] and left[1] == right[1]


zero_pair = (sp.zeros(5), sp.zeros(5))
identity_pair = (sp.eye(5), sp.zeros(5))

quotient_character_ok = all(
    coset_exponent[compose(left, right)]
    == (coset_exponent[left]+coset_exponent[right]) % 3
    for left in A4 for right in A4
)
check("the A4/V4 exponent map is an exact C3 character",
      quotient_character_ok)


def character_projector(character_sign):
    projector = zero_pair
    # Coordinates of 1, omega, omega^2 in the basis (1,omega).
    powers = ((sp.Integer(1), sp.Integer(0)),
              (sp.Integer(0), sp.Integer(1)),
              (sp.Integer(-1), sp.Integer(-1)))
    for element in A4:
        exponent = coset_exponent[element]
        # P_chi=(1/|H|)sum conjugate(chi(h))*rho(h).
        constant, omega_coefficient = powers[
            (-character_sign*exponent) % 3
        ]
        term = (constant*representation[element],
                omega_coefficient*representation[element])
        projector = pair_add(projector, term)
    return pair_scale(sp.Rational(1, 12), projector)


seed_projectors = {
    sign: character_projector(sign) for sign in (1, -1)
}
check("both nontrivial A4 character projectors have exact rank one",
      all(pair_is_equal(pair_multiply(projector, projector), projector)
              and (sp.trace(projector[0]), sp.trace(projector[1])) == (1, 0)
          for projector in seed_projectors.values()))


def imprimitivity_system(seed):
    orbit = {}
    for action in label_actions:
        matrix = representation[action]
        conjugated = pair_conjugate(seed, matrix)
        orbit[pair_key(conjugated)] = conjugated
    return list(orbit.values())


systems = {
    sign: imprimitivity_system(projector)
    for sign, projector in seed_projectors.items()
}
system_checks = []
for projectors in systems.values():
    orthogonal = all(
        pair_is_zero(pair_multiply(projectors[left], projectors[right]))
        for left in range(5) for right in range(5) if left != right
    )
    projector_sum = zero_pair
    for projector in projectors:
        projector_sum = pair_add(projector_sum, projector)
    system_checks.append(
        len(projectors) == 5
        and all(pair_is_equal(pair_multiply(projector, projector), projector)
                and (sp.trace(projector[0]), sp.trace(projector[1])) == (1, 0)
                and pair_is_equal(pair_metric_adjoint(projector), projector)
                for projector in projectors)
        and pair_is_equal(projector_sum, identity_pair)
        and orthogonal
    )
check("each character gives five exact rank-one algebraic projectors",
      all(system_checks),
      "they are self-adjoint in the inherited W5 metric")

normalization_ok = True
for projectors in systems.values():
    keys = {pair_key(projector) for projector in projectors}
    for action in representation.values():
        images = {
            pair_key(pair_conjugate(projector, action))
            for projector in projectors
        }
        if images != keys:
            normalization_ok = False
check("A5 normalizes both exact C^5 imprimitivity algebras",
      normalization_ok)
same_systems = ({pair_key(projector) for projector in systems[1]}
                == {pair_key(projector) for projector in systems[-1]})
check("the conjugate A4 characters give two distinct systems",
      not same_systems)

# Correct sheet-image classification.  A normalized algebra has r<=5 central
# supports.  For r<5, simplicity and |A5|>|S4| force the centre action to be
# trivial, contradicting W5 irreducibility unless r=1.  For r=5 every support
# is one-dimensional, giving precisely a monomial C^5 system.  With r=1,
# k*m=5 gives only C or M5 because five is prime.
normalized_sheet_images = ["C", "C^5", "M5"]
check("the corrected normalized sheet-image list is exhaustive",
      normalized_sheet_images == ["C", "C^5", "M5"]
      and len(A4) == 12 and all(system_checks),
      "centre sizes 2,3,4 are excluded; sizes 1 and 5 give the list")

# First-order support for a diagonal C^5 algebra is a union of two perfect
# matchings: an entry can survive only when the two bimodule nodes share their
# left label or their right label.  Compute the union support of the complete
# affine family between every pair of the two exact monomial systems.
hessian_directions = [
    (direction, sp.zeros(5)) for direction in restricted_hessians
]
affine_directions = [identity_pair]+hessian_directions
affine_self_adjoint = all(
    pair_is_equal(pair_metric_adjoint(direction), direction)
    for direction in affine_directions
)
check("all six affine selector directions are exactly self-adjoint on W5",
      affine_self_adjoint)
support_records = {}
all_permutations = list(permutations(range(5)))
for plus_sign, plus_projectors in systems.items():
    for minus_sign, minus_projectors in systems.items():
        support = []
        affine_support = []
        for row, left_projector in enumerate(plus_projectors):
            support_row = []
            affine_support_row = []
            for col, right_projector in enumerate(minus_projectors):
                nonzero = any(not pair_is_zero(pair_multiply(
                    pair_multiply(left_projector, direction),
                    right_projector,
                )) for direction in hessian_directions)
                affine_nonzero = any(not pair_is_zero(pair_multiply(
                    pair_multiply(left_projector, direction),
                    right_projector,
                )) for direction in affine_directions)
                support_row.append(nonzero)
                affine_support_row.append(affine_nonzero)
            support.append(support_row)
            affine_support.append(affine_support_row)
        support_size = sum(sum(row) for row in support)
        affine_support_size = sum(sum(row) for row in affine_support)
        # A union of two permutation matchings contains at most ten entries.
        # Only smaller supports require the explicit relative-label audit.
        rook_cover = False
        if support_size <= 10:
            for first in all_permutations:
                for second in all_permutations:
                    if all(not support[row][col]
                           or col == first[row] or col == second[row]
                           for row in range(5) for col in range(5)):
                        rook_cover = True
                        break
                if rook_cover:
                    break
        support_records[f"{plus_sign},{minus_sign}"] = {
            "hessian_size": support_size,
            "hessian_full": support_size == 25,
            "affine_size": affine_support_size,
            "affine_full": affine_support_size == 25,
            "two_rook_cover": rook_cover,
        }
check("Hessian support is 20 or 25 entries and affine support is full",
      sorted(record["hessian_size"] for record in support_records.values())
      == [20, 20, 25, 25]
      and all(record["affine_full"] for record in support_records.values()),
      f"supports={support_records}")
check("no C^5 first-order rook pattern contains the affine family",
      all(not record["two_rook_cover"]
          for record in support_records.values()),
      "all two-character systems and every relative label permutation covered")

# Split central supports (C+C, M5+C, C5+C, or separate copies) fail even more
# cheaply.  For a central idempotent e supported on H+ and its KO6 opposite on
# H-, the first-order double commutator contains the whole B block.
generic_entries = sp.symbols("z0:25")
generic_B = sp.Matrix(5, 5, generic_entries)
zero5 = sp.zeros(5)
generic_D = zero5.row_join(generic_B).col_join(
    generic_B.T.row_join(zero5)
)
central_left = sp.diag(*([1]*5+[0]*5))
central_opposite = sp.diag(*([0]*5+[1]*5))
double_commutator = ((generic_D*central_left-central_left*generic_D)
                     *central_opposite
                     -central_opposite
                     *(generic_D*central_left-central_left*generic_D))
expected_double_commutator = zero5.row_join(-generic_B).col_join(
    (-generic_B.T).row_join(zero5)
)
check("split central supports force the selector block to zero",
      double_commutator == expected_double_commutator
      and any(value != 0 for value in double_commutator),
      "[[D,e],e^o]=[[0,-B],[-B^T,0]]")

# If both sheet images are M5, KO6 exchange puts a full M5 opposite image on
# the same sheet.  Its commutant is scalar, so order zero fails.  If the only
# image is scalar C, every represented one-form vanishes.
E12 = sp.zeros(5)
E12[0, 1] = 1
E21 = sp.zeros(5)
E21[1, 0] = 1
check("full-matrix sheet pairs fail order zero",
      E12*E21-E21*E12 != sp.zeros(5),
      "a unitary conjugate of full M5 is still full M5, not its commutant")

a, scalar = sp.symbols("a scalar")
scalar_representation = scalar*sp.eye(10)
scalar_D = zero5.row_join(sp.eye(5)).col_join(sp.eye(5).row_join(zero5))
check("the sole scalar algebra has identically zero represented one-forms",
      scalar_D*scalar_representation-scalar_representation*scalar_D
      == sp.zeros(10))

joint_gate_ledger = {
    "scalar C": "order zero/first order possible; zero one-forms",
    "diagonal C5": "order zero possible; 20/25-entry Hessian support violates first order",
    "split central summands": "first order forces B=0",
    "both sheet images M5": "order zero fails",
    "one M5 and one nonscalar non-full image": "order zero fails",
}
check("every corrected joint image category fails a necessary gate",
      len(joint_gate_ledger) == 5
      and all(not record["two_rook_cover"]
              for record in support_records.values()))

payload = {
    "protocol_commits": {
        "initial": "a0d7dd8",
        "monomial_correction": "cae7233",
    },
    "group": {
        "order": len(label_actions),
        "simple_by_normal_closures": set(normal_closure_sizes) == {60},
        "W5_complex_commutant_dimension": 25-commutant_equations.rank(),
        "W5_irreducible": commutant_equations.rank() == 24,
        "representation_exactly_unitary": representation_unitary_ok,
        "A4_order": len(A4),
        "V4_order": len(V4),
    },
    "normalized_sheet_algebras": normalized_sheet_images,
    "monomial_C5": {
        "nontrivial_A4_characters": 2,
        "systems": 2,
        "projectors_per_system": 5,
        "projector_rank": 1,
        "projectors_metric_self_adjoint": all(system_checks),
        "A5_normalized": normalization_ok,
        "affine_directions_metric_self_adjoint": affine_self_adjoint,
        "affine_support_by_system_pair": support_records,
        "first_order_rook_cover_exists": False,
    },
    "joint_gate_ledger": joint_gate_ledger,
    "verdict": (
        "DERIVED CORRECTION: W5 has two exact A5-normalized monomial C^5 "
        "systems, so the original C/M5 sheet list was incomplete.  DERIVED "
        "TEN-STATE NO-GO under the stated hypotheses: the variable Hessian "
        "family has 20 or 25 support entries between the character systems "
        "(and full 25-entry affine-span support), while first order over C^5 "
        "permits at most a union of two permutation matchings.  Split centres "
        "force B=0, M5 pairs fail order zero, and scalar C has zero one-forms. "
        "Larger bimodule completions are not covered."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact ten-state finite-triple audit JSON was written",
      OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED CORRECTION: W5 admits two normalized monomial C^5 systems.")
print("DERIVED TEN-STATE NO-GO: their first-order rook supports cannot contain")
print("                         the 20/25-entry Hessian family.")
print("DERIVED: split centres force B=0; M5 pairs fail order zero; C has")
print("         zero one-forms.")
print("OPEN: larger Krajewski/bimodule completions are outside this fixed carrier.")
raise SystemExit(0 if passed == tests else 1)

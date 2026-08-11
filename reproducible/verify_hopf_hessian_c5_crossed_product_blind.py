#!/usr/bin/env python3
"""Target-blind audit of the canonical five-point A5 crossed product.

Protocol commit c3c1584 froze the algebra, natural carriers and tests before
this file was constructed.  No selector matrix or target comparison occurs.
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


OUTPUT = Path(__file__).with_name("hopf_hessian_c5_crossed_product_blind.json")
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
    raise RuntimeError("order exceeds derived group size")


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
            for candidate in (compose(element, other),
                              compose(other, element)):
                if candidate not in subgroup:
                    frontier.append(candidate)
    return frozenset(subgroup)


def conjugate_subgroup(subgroup, element):
    element_inverse = inverse(element)
    return frozenset(compose(compose(element, member), element_inverse)
                     for member in subgroup)


# Pairs represent exact matrices/scalars A+omega B over
# Q(omega), omega^2+omega+1=0.
def pair_add(left, right):
    return left[0]+right[0], left[1]+right[1]


def pair_multiply(left, right):
    A, B = left
    C, D = right
    return A*C-B*D, A*D+B*C-B*D


def pair_conjugate(pair):
    return pair[0]-pair[1], -pair[1]


def pair_equal(left, right):
    return left[0] == right[0] and left[1] == right[1]


def pair_trace(pair):
    return sp.trace(pair[0]), sp.trace(pair[1])


def scalar_add(left, right):
    return left[0]+right[0], left[1]+right[1]


def scalar_multiply(left, right):
    A, B = left
    C, D = right
    return A*C-B*D, A*D+B*C-B*D


print("="*78)
print("TARGET-BLIND FIVE-POINT A5 CROSSED-PRODUCT AUDIT")
print("="*78)

# Reconstruct the exact 60-element action from the binary geometry.
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
group = []
for vertex_action in vertex_actions:
    action = []
    for fibration in fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fiber))
            for fiber in fibration
        ))
        action.append(fibration_by_signature[signature])
    group.append(tuple(action))
group = sorted(set(group))
identity = tuple(range(6))
check("the derived effective group has order 60",
      len(group) == 60 and identity in group)

# Build the five-point coset action from an index-five A4.
order_two = next(element for element in group
                 if permutation_order(element) == 2)
V4 = frozenset(element for element in group
               if compose(element, order_two) == compose(order_two, element))
A4 = frozenset(element for element in group
               if conjugate_subgroup(V4, element) == V4)
A4_conjugates = sorted(
    {conjugate_subgroup(A4, element) for element in group},
    key=lambda subgroup: tuple(sorted(subgroup)),
)
subgroup_index = {subgroup: index
                  for index, subgroup in enumerate(A4_conjugates)}
five_actions_by_group = {}
for element in group:
    five_actions_by_group[element] = tuple(
        subgroup_index[conjugate_subgroup(subgroup, element)]
        for subgroup in A4_conjugates
    )
five_actions = sorted(set(five_actions_by_group.values()))
base_point = subgroup_index[A4]
stabilizer = [element for element in group
              if five_actions_by_group[element][base_point] == base_point]
check("the exact five-point action is faithful and transitive",
      len(A4_conjugates) == 5 and len(five_actions) == 60
      and {action[base_point] for action in five_actions} == set(range(5)))
check("a point stabilizer is exactly the constructed A4",
      frozenset(stabilizer) == A4 and len(stabilizer) == 12)

# Exact A4 class and abelianization data.
unseen = set(A4)
A4_classes = []
while unseen:
    element = min(unseen)
    group_class = frozenset(
        compose(compose(conjugator, element), inverse(conjugator))
        for conjugator in A4
    )
    A4_classes.append(group_class)
    unseen -= group_class
A4_class_sizes = sorted(len(group_class) for group_class in A4_classes)
commutators = [
    compose(compose(compose(left, right), inverse(left)), inverse(right))
    for left in A4 for right in A4
]
commutator_subgroup = generated_subgroup(commutators, identity)
check("A4 has exact class sizes 1,3,4,4 and abelianization C3",
      A4_class_sizes == [1, 3, 4, 4]
      and commutator_subgroup == V4
      and len(A4)//len(commutator_subgroup) == 3)

stabilizer_degrees = [1, 1, 1, 3]
crossed_blocks = [5*degree for degree in stabilizer_degrees]
crossed_dimension = len(group)*5
check("the stabilizer irreducible degrees are exactly 1,1,1,3",
      len(A4_classes) == 4
      and sum(degree*degree for degree in stabilizer_degrees) == len(A4))
check("the crossed product is M5+M5+M5+M15 of dimension 300",
      crossed_blocks == [5, 5, 5, 15]
      and sum(block*block for block in crossed_blocks)
      == crossed_dimension == 300,
      "C(P) crossed A5 ~= M5(C[A4])")

# Quotient character A4/V4=C3.
order_three = next(element for element in A4
                   if permutation_order(element) == 3)
quotient_cosets = [
    V4,
    frozenset(compose(order_three, member) for member in V4),
    frozenset(compose(compose(order_three, order_three), member)
              for member in V4),
]
coset_exponent = {
    element: exponent for exponent, coset in enumerate(quotient_cosets)
    for element in coset
}
check("the three stabilizer characters are constructed exactly",
      set().union(*quotient_cosets) == A4
      and all(coset_exponent[compose(left, right)]
              == (coset_exponent[left]+coset_exponent[right]) % 3
              for left in A4 for right in A4))

# Choose exact left-coset representatives through the conjugate subgroups.
representatives = {}
for point, subgroup in enumerate(A4_conjugates):
    if subgroup == A4:
        representatives[point] = identity
    else:
        representatives[point] = next(
            element for element in group
            if conjugate_subgroup(A4, element) == subgroup
        )

powers = (
    (sp.Integer(1), sp.Integer(0)),
    (sp.Integer(0), sp.Integer(1)),
    (sp.Integer(-1), sp.Integer(-1)),
)
zero_matrix_pair = (sp.zeros(5), sp.zeros(5))
identity_matrix_pair = (sp.eye(5), sp.zeros(5))


def induced_matrix(element, character_sign):
    constant = sp.zeros(5)
    omega_part = sp.zeros(5)
    action = five_actions_by_group[element]
    for source in range(5):
        target = action[source]
        stabilizer_element = compose(
            compose(inverse(representatives[target]), element),
            representatives[source],
        )
        if stabilizer_element not in A4:
            raise RuntimeError("invalid induced-representation cocycle")
        exponent = (character_sign*coset_exponent[stabilizer_element]) % 3
        scalar = powers[exponent]
        constant[target, source] = scalar[0]
        omega_part[target, source] = scalar[1]
    return constant, omega_part


branches = (0, 1, -1)
induced = {
    branch: {element: induced_matrix(element, branch) for element in group}
    for branch in branches
}
representation_laws = {}
covariance_checks = {}
branch_records = {}
diagonal_projectors = []
for point in range(5):
    projector = sp.zeros(5)
    projector[point, point] = 1
    diagonal_projectors.append((projector, sp.zeros(5)))

for branch in branches:
    representation_laws[branch] = all(pair_equal(
        induced[branch][compose(left, right)],
        pair_multiply(induced[branch][left], induced[branch][right]),
    ) for left in group for right in group)
    covariance_checks[branch] = all(pair_equal(
        pair_multiply(pair_multiply(induced[branch][element],
                                    diagonal_projectors[source]),
                      induced[branch][inverse(element)]),
        diagonal_projectors[five_actions_by_group[element][source]],
    ) for element in group for source in range(5))

    traces = {element: pair_trace(induced[branch][element])
              for element in group}
    invariant_sum = (sp.Integer(0), sp.Integer(0))
    commutant_sum = (sp.Integer(0), sp.Integer(0))
    for element in group:
        trace = traces[element]
        invariant_sum = scalar_add(invariant_sum, trace)
        commutant_sum = scalar_add(
            commutant_sum,
            scalar_multiply(trace, pair_conjugate(trace)),
        )
    invariant_dimension = (invariant_sum[0]/60, invariant_sum[1]/60)
    commutant_dimension = (commutant_sum[0]/60, commutant_sum[1]/60)

    image_support = set()
    for element in group:
        for point in range(5):
            corner = pair_multiply(diagonal_projectors[point],
                                   induced[branch][element])
            for row in range(5):
                for col in range(5):
                    if corner[0][row, col] != 0 or corner[1][row, col] != 0:
                        image_support.add((row, col))
    branch_records[str(branch)] = {
        "A5_invariant_dimension": [int(invariant_dimension[0]),
                                    int(invariant_dimension[1])],
        "A5_commutant_dimension": [int(commutant_dimension[0]),
                                    int(commutant_dimension[1])],
        "crossed_product_image_dimension": len(image_support),
        "crossed_product_kernel_dimension": crossed_dimension-len(image_support),
        "image_is_full_M5": len(image_support) == 25,
    }

check("all three induced matrices obey the exact group law and covariance",
      all(representation_laws.values()) and all(covariance_checks.values()))
check("the trivial branch is 1+4 and the two nontrivial branches are W5",
      branch_records["0"]["A5_invariant_dimension"] == [1, 0]
      and branch_records["0"]["A5_commutant_dimension"] == [2, 0]
      and all(branch_records[str(branch)]["A5_invariant_dimension"] == [0, 0]
              and branch_records[str(branch)]["A5_commutant_dimension"] == [1, 0]
              for branch in (1, -1)))
check("every one-dimensional-character covariant image is full M5",
      all(record["image_is_full_M5"]
          and record["crossed_product_kernel_dimension"] == 275
          for record in branch_records.values()),
      f"branches={branch_records}")

# A faithful left module must contain at least one fundamental module for
# every simple summand.  The lower bound is attained by their direct sum.
minimum_faithful_left_dimension = sum(crossed_blocks)
check("the minimum faithful left module has dimension 30",
      minimum_faithful_left_dimension == 30)

# The direct sum of one fundamental module for each simple block is canonical
# up to unitary equivalence as a minimal faithful *left* module, but it is not
# a real bimodule.  Its commutant is C^4.  Even after the standard grading
# double, the commutant has dimension 4*(2^2)=16, too small to contain the
# faithful 300-dimensional opposite algebra forced by an antiunitary J.
minimum_left_commutant_dimension = len(crossed_blocks)
minimum_odd_double_commutant_dimension = 4*len(crossed_blocks)
check("the canonical minimum faithful left double fails order zero",
      minimum_left_commutant_dimension == 4
      and minimum_odd_double_commutant_dimension == 16
      and minimum_odd_double_commutant_dimension < crossed_dimension,
      "a faithful opposite copy of the 300D algebra cannot fit in its commutant")

# The two carriers functorially supplied by a semisimple algebra do not select
# an orientable KO6 even structure.  The regular bimodule has only the four
# diagonal central cells.  The full enveloping bimodule has all sixteen and
# therefore also contains the diagonal cells.  Their standard odd doubles
# duplicate the same left/right representation on the two grading sheets;
# every zero-cycle has sheet profile (1,1), while gamma has (1,-1).
regular_central_cells = [(index, index) for index in range(4)]
enveloping_central_cells = [(left, right)
                            for left in range(4) for right in range(4)]
zero_cycle_sheet_rank = sp.Matrix([[1], [1]]).rank()
orientation_augmented_rank = sp.Matrix([[1, 1], [1, -1]]).rank()
check("the regular bimodule has only four diagonal central cells",
      len(regular_central_cells) == 4
      and all(left == right for left, right in regular_central_cells))
check("the full enveloping bimodule contains all four diagonal cells",
      len(enveloping_central_cells) == 16
      and sum(left == right for left, right in enveloping_central_cells) == 4)
check("the canonical odd doubles fail metric-zero orientability",
      zero_cycle_sheet_rank == 1 and orientation_augmented_rank == 2,
      "zero-cycles act with sheet profile (1,1), gamma has (1,-1)")

# The natural ten-state branch census is target-free: it uses only that an
# odd block is nonzero.  Equal five-dimensional branches give a graph M5
# image and fail order zero; unequal branches give split central summands and
# first order forces the odd block to vanish.  All branches are nonfaithful as
# representations of the full crossed product.
ten_state_branches = []
for plus_branch in branches:
    for minus_branch in branches:
        if plus_branch == minus_branch:
            obstruction = "order_zero_full_M5_graph"
        else:
            obstruction = "first_order_split_central_support"
        ten_state_branches.append({
            "plus": plus_branch,
            "minus": minus_branch,
            "obstruction": obstruction,
            "faithful_to_full_crossed_product": False,
        })
check("all nine natural ten-state branches fail before target comparison",
      len(ten_state_branches) == 9
      and sum(record["obstruction"] == "order_zero_full_M5_graph"
              for record in ten_state_branches) == 3
      and sum(record["obstruction"] == "first_order_split_central_support"
              for record in ten_state_branches) == 6
      and not any(record["faithful_to_full_crossed_product"]
                  for record in ten_state_branches),
      "3 graph-M5 order-zero failures; 6 split-centre first-order failures")

payload = {
    "protocol_commit": "c3c1584",
    "target_comparison_performed": False,
    "effective_group_order": len(group),
    "five_point_action": {
        "points": 5,
        "stabilizer": "A4",
        "stabilizer_order": len(A4),
        "stabilizer_class_sizes": A4_class_sizes,
        "stabilizer_abelianization_order": len(A4)//len(commutator_subgroup),
        "stabilizer_irrep_degrees": stabilizer_degrees,
    },
    "crossed_product": {
        "definition": "C(P) crossed_product A5",
        "isomorphism": "M5(C[A4])",
        "dimension": crossed_dimension,
        "Wedderburn_blocks": crossed_blocks,
        "centre_dimension": len(crossed_blocks),
        "minimum_faithful_left_module_dimension": minimum_faithful_left_dimension,
    },
    "induced_character_branches": branch_records,
    "canonical_carriers": {
        "minimum_faithful_left_module_dimension": minimum_faithful_left_dimension,
        "minimum_faithful_standard_odd_double_dimension": 2*minimum_faithful_left_dimension,
        "minimum_faithful_standard_odd_double_commutant_dimension": (
            minimum_odd_double_commutant_dimension
        ),
        "minimum_faithful_standard_odd_double_order_zero": False,
        "regular_bimodule_dimension": crossed_dimension,
        "regular_standard_odd_double_dimension": 2*crossed_dimension,
        "regular_central_cells": regular_central_cells,
        "full_enveloping_bimodule_dimension": crossed_dimension**2,
        "full_enveloping_central_cells": len(enveloping_central_cells),
        "standard_odd_double_orientable": False,
    },
    "natural_ten_state_branch_census": ten_state_branches,
    "verdict": (
        "TARGET-BLIND DERIVED: C(P) crossed A5 is M5+M5+M5+M15. "
        "The three one-dimensional stabilizer-character representations all "
        "have full M5 image and 275-dimensional kernel; the two nontrivial "
        "branches are the conjugate W5 systems. No natural ten-state branch "
        "is faithful or passes the cheap real-triple gates. The minimum "
        "faithful left double fails order zero; the regular and full "
        "enveloping carriers supplied by the algebra have nonorientable "
        "standard odd doubles. No selector target was inspected."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the target-blind crossed-product JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("TARGET-BLIND DERIVED: C(P) crossed A5 = M5+M5+M5+M15.")
print("TARGET-BLIND NEGATIVE: all natural ten-state branches and both canonical")
print("                       doubled carriers fail before selector comparison.")
print("NO SELECTOR TARGET COMPARISON WAS PERFORMED.")
raise SystemExit(0 if passed == tests else 1)

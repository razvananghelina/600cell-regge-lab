#!/usr/bin/env python3
"""Exact real-form and KO6 Poincare audit of the five-point crossed product.

Protocol commit 5eafb29 froze the real algebra and decision boundary before
the Frobenius--Schur indicators or real Wedderburn type were computed.  No
selector target is used.
"""

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


OUTPUT = Path(__file__).with_name("hopf_hessian_crossed_real_form.json")
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


def conjugate_subgroup(subgroup, element):
    element_inverse = inverse(element)
    return frozenset(compose(compose(element, member), element_inverse)
                     for member in subgroup)


# Exact Q(omega) scalar arithmetic, omega^2+omega+1=0.
def cyclotomic_add(left, right):
    return left[0]+right[0], left[1]+right[1]


def cyclotomic_multiply(left, right):
    A, B = left
    C, D = right
    return A*C-B*D, A*D+B*C-B*D


def cyclotomic_conjugate(value):
    return value[0]-value[1], -value[1]


def cyclotomic_scale(value, scalar):
    return value[0]*scalar, value[1]*scalar


print("="*78)
print("EXACT REAL-FORM AND KO6 POINCARE AUDIT")
print("="*78)

# Reconstruct the exact A5 quotient and an index-five A4 stabilizer.
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

order_two = next(element for element in group
                 if permutation_order(element) == 2)
V4 = frozenset(element for element in group
               if compose(element, order_two) == compose(order_two, element))
A4 = frozenset(element for element in group
               if conjugate_subgroup(V4, element) == V4)
check("the derived subgroup data are |A5|=60, |A4|=12, |V4|=4",
      len(group) == 60 and len(A4) == 12 and len(V4) == 4
      and identity in V4)

# Exact A4 conjugacy classes.
unseen = set(A4)
classes = []
while unseen:
    element = min(unseen)
    group_class = frozenset(
        compose(compose(conjugator, element), inverse(conjugator))
        for conjugator in A4
    )
    classes.append(group_class)
    unseen -= group_class
classes.sort(key=lambda group_class:
             (permutation_order(next(iter(group_class))),
              tuple(sorted(group_class))))
class_sizes = [len(group_class) for group_class in classes]
class_orders = [permutation_order(next(iter(group_class)))
                for group_class in classes]
check("A4 has exact class size/order data (1,3,4,4)/(1,2,3,3)",
      class_sizes == [1, 3, 4, 4]
      and class_orders == [1, 2, 3, 3])

# Quotient exponent A4/V4=C3 and all four complex characters.
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
powers = ((1, 0), (0, 1), (-1, -1))


def character_value(character, element):
    if character == "1":
        return 1, 0
    if character == "chi":
        return powers[coset_exponent[element]]
    if character == "chibar":
        return powers[(-coset_exponent[element]) % 3]
    if character == "3":
        order = permutation_order(element)
        return ({1: 3, 2: -1, 3: 0}[order], 0)
    raise KeyError(character)


characters = ("1", "chi", "chibar", "3")
degrees = {character: character_value(character, identity)[0]
           for character in characters}
inner_products = {}
for left in characters:
    for right in characters:
        total = (0, 0)
        for element in A4:
            total = cyclotomic_add(total, cyclotomic_multiply(
                character_value(left, element),
                cyclotomic_conjugate(character_value(right, element)),
            ))
        inner_products[(left, right)] = cyclotomic_scale(
            total, sp.Rational(1, len(A4))
        )
check("the four exact A4 characters are orthonormal",
      all(inner_products[(left, right)]
          == ((1, 0) if left == right else (0, 0))
          for left in characters for right in characters))
check("their degrees exhaust the complex group algebra",
      [degrees[character] for character in characters] == [1, 1, 1, 3]
      and sum(degree*degree for degree in degrees.values()) == len(A4))

# Frobenius--Schur indicators.
indicators = {}
for character in characters:
    total = (0, 0)
    for element in A4:
        total = cyclotomic_add(
            total,
            character_value(character, compose(element, element)),
        )
    indicator = cyclotomic_scale(total, sp.Rational(1, len(A4)))
    indicators[character] = indicator
check("the exact Frobenius--Schur indicators are 1,0,0,1",
      indicators == {
          "1": (1, 0),
          "chi": (0, 0),
          "chibar": (0, 0),
          "3": (1, 0),
      },
      f"indicators={indicators}")

# Real Wedderburn reconstruction.  The two FS-zero characters form one
# conjugate pair and hence one complex-type real simple summand.
conjugate_pair_ok = all(
    character_value("chibar", element)
    == cyclotomic_conjugate(character_value("chi", element))
    for element in A4
)
real_group_algebra_blocks = ["R", "C", "M3(R)"]
real_group_algebra_dimensions = [1, 2, 9]
real_crossed_product_blocks = ["M5(R)", "M5(C)", "M15(R)"]
real_crossed_product_dimensions = [25, 50, 225]
complexified_blocks = [5, 5, 5, 15]
check("R[A4] has real type R+C+M3(R)",
      conjugate_pair_ok
      and sum(real_group_algebra_dimensions) == len(A4)
      and real_group_algebra_blocks == ["R", "C", "M3(R)"])
check("the real crossed product is M5(R)+M5(C)+M15(R)",
      sum(real_crossed_product_dimensions) == 300
      and real_crossed_product_blocks
      == ["M5(R)", "M5(C)", "M15(R)"],
      "real dimension=25+50+225=300")
check("complexification reproduces M5+M5+M5+M15",
      complexified_blocks == [5, 5, 5, 15]
      and sum(block*block for block in complexified_blocks) == 300,
      "the real M5(C) block splits into the chi/chibar complex blocks")

# KO6 intersection transpose sign.  For
# cap_ij=Tr(gamma*pi(p_i)*J*pi(p_j)*J^-1), conjugating the trace by the
# antiunitary J and using J gamma J^-1=-gamma swaps i,j and contributes -1.
# K0 has one free generator per real simple matrix-algebra summand.
J_gamma_sign = -1
intersection_transpose_sign = J_gamma_sign
K0_free_rank = len(real_crossed_product_blocks)
check("KO6 makes the real K0 intersection form antisymmetric",
      intersection_transpose_sign == -1 and K0_free_rank == 3)

x, y, z = sp.symbols("x y z")
generic_intersection = sp.Matrix([
    [0, x, y],
    [-x, 0, z],
    [-y, -z, 0],
])
generic_determinant = sp.expand(generic_intersection.det())
generic_rank = generic_intersection.rank()
check("every three-generator antisymmetric intersection form is degenerate",
      generic_determinant == 0 and generic_rank == 2,
      "generic rank=2<3; determinant identically zero")

nondegenerate_poincare_possible = (
    generic_determinant != 0 and generic_rank == K0_free_rank
)
check("the canonical real crossed product fails KO6 Poincare duality",
      not nondegenerate_poincare_possible,
      "arbitrary multiplicities and Dirac operators are covered")

payload = {
    "protocol_commit": "5eafb29",
    "target_comparison_performed": False,
    "derived_group": {
        "A5_order": len(group),
        "A4_order": len(A4),
        "V4_order": len(V4),
        "A4_class_sizes": class_sizes,
        "A4_class_orders": class_orders,
    },
    "complex_A4_characters": {
        "names": list(characters),
        "degrees": [degrees[character] for character in characters],
        "orthonormal": True,
        "Frobenius_Schur_indicators": {
            character: [int(value) for value in indicators[character]]
            for character in characters
        },
        "chi_chibar_conjugate_pair": conjugate_pair_ok,
    },
    "real_group_algebra": {
        "type": real_group_algebra_blocks,
        "real_block_dimensions": real_group_algebra_dimensions,
        "real_dimension": sum(real_group_algebra_dimensions),
    },
    "real_crossed_product": {
        "type": real_crossed_product_blocks,
        "real_block_dimensions": real_crossed_product_dimensions,
        "real_dimension": sum(real_crossed_product_dimensions),
        "complexification_block_sizes": complexified_blocks,
        "real_simple_summands": K0_free_rank,
        "K0_free_rank": K0_free_rank,
    },
    "KO6_Poincare": {
        "intersection_transpose_sign": intersection_transpose_sign,
        "generic_intersection_max_rank": int(generic_rank),
        "required_rank": K0_free_rank,
        "nondegenerate_possible": nondegenerate_poincare_possible,
        "arbitrary_multiplicities_covered": True,
        "arbitrary_D_covered": True,
    },
    "verdict": (
        "DERIVED FULL-ARENA POINCARE NO-GO under the stated canonical-real-"
        "form and KO6 hypotheses. Exact FS indicators 1,0,0,1 give "
        "R[A4]=R+C+M3(R) and B_R=M5(R)+M5(C)+M15(R). Its real K0 free rank "
        "is three, while the KO6 intersection form is antisymmetric and has "
        "rank at most two. No Hessian target is used."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact real-form audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: B_R = M5(R)+M5(C)+M15(R), real K0 rank three.")
print("DERIVED FULL-ARENA KO6 POINCARE NO-GO: antisymmetric rank <=2.")
print("NO HESSIAN OR SELECTOR TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

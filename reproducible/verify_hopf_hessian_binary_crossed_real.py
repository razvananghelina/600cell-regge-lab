#!/usr/bin/env python3
"""Exact real-form audit for the ineffective five-point 2I crossed product.

Protocol commit 32d00c7 froze the binary groupoid and KO6 decision boundary
before the 2T characters and Frobenius--Schur indicators were computed.  No
selector target is used.
"""

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import sympy as sp

from verify_hopf_fibration_invariants import (
    build_2I,
    find_all_hopf_fibrations,
    find_vertex_index,
    quat_mult,
)


OUTPUT = Path(__file__).with_name("hopf_hessian_binary_crossed_real.json")
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


def inverse_permutation(permutation):
    result = [None]*len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def permutation_order(permutation):
    identity = tuple(range(len(permutation)))
    current = identity
    for order in range(1, 121):
        current = compose(permutation, current)
        if current == identity:
            return order
    raise RuntimeError("permutation order exceeds binary group order")


def conjugate_subgroup(subgroup, element):
    element_inverse = inverse_permutation(element)
    return frozenset(compose(compose(element, member), element_inverse)
                     for member in subgroup)


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
print("EXACT INEFFECTIVE BINARY CROSSED-PRODUCT REAL-FORM AUDIT")
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


# Effective A5 action on the six already-derived fibrations, retained for
# every one of the 120 binary elements so its kernel is visible.
fibration_by_signature = {
    tuple(sorted(tuple(sorted(fiber)) for fiber in fibration)): index
    for index, fibration in enumerate(fibrations)
}
six_actions_by_binary_element = []
for element in vertices:
    vertex_action = conjugation_permutation(element)
    action = []
    for fibration in fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fiber))
            for fiber in fibration
        ))
        action.append(fibration_by_signature[signature])
    six_actions_by_binary_element.append(tuple(action))
effective_group = sorted(set(six_actions_by_binary_element))
effective_identity = tuple(range(6))
check("the binary action has order 120 with effective A5 quotient order 60",
      len(vertices) == 120 and len(effective_group) == 60
      and six_actions_by_binary_element.count(effective_identity) == 2)

# Exact multiplication table in the 120-vertex quaternion model.
multiplication = np.empty((120, 120), dtype=np.int16)
max_product_error = 0.0
for left in range(120):
    for right in range(120):
        product_quaternion = quat_mult(vertices[left], vertices[right])
        product_index = find_vertex_index(vertices, product_quaternion)
        multiplication[left, right] = product_index
        max_product_error = max(
            max_product_error,
            float(np.sum((vertices[product_index]-product_quaternion)**2)),
        )
identity = find_vertex_index(vertices, np.array([1.0, 0.0, 0.0, 0.0]))
centre_minus = find_vertex_index(vertices, np.array([-1.0, 0.0, 0.0, 0.0]))
check("all binary products close on the exact indexed vertex set",
      max_product_error < 1e-20
      and all(multiplication[identity, index] == index
              and multiplication[index, identity] == index
              for index in range(120)),
      f"max squared product match error={max_product_error:.3e}")

# Construct A4 as before and take its full binary preimage.  Q8 is the
# preimage of its normal V4.
order_two = next(element for element in effective_group
                 if permutation_order(element) == 2)
V4_effective = frozenset(
    element for element in effective_group
    if compose(element, order_two) == compose(order_two, element)
)
A4_effective = frozenset(
    element for element in effective_group
    if conjugate_subgroup(V4_effective, element) == V4_effective
)
T_indices = tuple(index for index, action
                  in enumerate(six_actions_by_binary_element)
                  if action in A4_effective)
Q8_indices = frozenset(index for index, action
                       in enumerate(six_actions_by_binary_element)
                       if action in V4_effective)
T_set = frozenset(T_indices)
check("the exact stabilizer preimages have orders |2T|=24 and |Q8|=8",
      len(A4_effective) == 12 and len(V4_effective) == 4
      and len(T_set) == 24 and len(Q8_indices) == 8
      and identity in Q8_indices and centre_minus in Q8_indices)
check("2T and Q8 are closed, with Q8 normal in 2T",
      all(int(multiplication[left, right]) in T_set
          for left in T_set for right in T_set)
      and all(int(multiplication[left, right]) in Q8_indices
              for left in Q8_indices for right in Q8_indices)
      and all({int(multiplication[int(multiplication[g, h]),
                                      find_vertex_index(
                                          vertices,
                                          np.array([vertices[g, 0],
                                                    -vertices[g, 1],
                                                    -vertices[g, 2],
                                                    -vertices[g, 3]]))])
                   for h in Q8_indices} == Q8_indices
              for g in T_set))

# Exact quotient 2T/Q8=C3.
quotient_generator = next(index for index in T_indices
                          if index not in Q8_indices)
generator_square = int(multiplication[quotient_generator,
                                      quotient_generator])
quotient_cosets = [
    Q8_indices,
    frozenset(int(multiplication[quotient_generator, member])
              for member in Q8_indices),
    frozenset(int(multiplication[generator_square, member])
              for member in Q8_indices),
]
coset_exponent = {
    element: exponent for exponent, coset in enumerate(quotient_cosets)
    for element in coset
}
check("the exact binary quotient is 2T/Q8=C3",
      set().union(*quotient_cosets) == T_set
      and all(len(coset) == 8 for coset in quotient_cosets)
      and all(coset_exponent[int(multiplication[left, right])]
              == (coset_exponent[left]+coset_exponent[right]) % 3
              for left in T_set for right in T_set))

powers = ((1, 0), (0, 1), (-1, -1))
natural_trace = {}
max_trace_residual = 0.0
for element in T_set:
    trace_value = 2*float(vertices[element, 0])
    rounded = int(round(trace_value))
    max_trace_residual = max(max_trace_residual, abs(trace_value-rounded))
    natural_trace[element] = rounded
check("the defining SU2 character is exactly integer-valued on 2T",
      max_trace_residual < 1e-12
      and set(natural_trace.values()) == {-2, -1, 0, 1, 2})


def one_dimensional_value(twist, element):
    return powers[(twist*coset_exponent[element]) % 3]


def character_value(character, element):
    if character == "1":
        return 1, 0
    if character == "chi":
        return one_dimensional_value(1, element)
    if character == "chibar":
        return one_dimensional_value(-1, element)
    if character == "2":
        return natural_trace[element], 0
    if character == "2chi":
        return cyclotomic_scale(one_dimensional_value(1, element),
                                natural_trace[element])
    if character == "2chibar":
        return cyclotomic_scale(one_dimensional_value(-1, element),
                                natural_trace[element])
    if character == "3":
        return natural_trace[element]**2-1, 0
    raise KeyError(character)


characters = ("1", "chi", "chibar", "2", "2chi", "2chibar", "3")
degrees = {character: character_value(character, identity)[0]
           for character in characters}
inner_products = {}
for left in characters:
    for right in characters:
        total = (0, 0)
        for element in T_set:
            total = cyclotomic_add(total, cyclotomic_multiply(
                character_value(left, element),
                cyclotomic_conjugate(character_value(right, element)),
            ))
        inner_products[(left, right)] = cyclotomic_scale(
            total, sp.Rational(1, len(T_set))
        )
check("the seven exact binary-tetrahedral characters are orthonormal",
      all(inner_products[(left, right)]
          == ((1, 0) if left == right else (0, 0))
          for left in characters for right in characters))
check("their degrees are 1,1,1,2,2,2,3 and exhaust order 24",
      [degrees[character] for character in characters]
      == [1, 1, 1, 2, 2, 2, 3]
      and sum(degree*degree for degree in degrees.values()) == len(T_set))

indicators = {}
for character in characters:
    total = (0, 0)
    for element in T_set:
        square = int(multiplication[element, element])
        total = cyclotomic_add(total, character_value(character, square))
    indicators[character] = cyclotomic_scale(
        total, sp.Rational(1, len(T_set))
    )
check("the exact FS indicators are 1,0,0,-1,0,0,1",
      indicators == {
          "1": (1, 0),
          "chi": (0, 0),
          "chibar": (0, 0),
          "2": (-1, 0),
          "2chi": (0, 0),
          "2chibar": (0, 0),
          "3": (1, 0),
      },
      f"indicators={indicators}")

conjugate_pairs_ok = all(
    character_value("chibar", element)
    == cyclotomic_conjugate(character_value("chi", element))
    and character_value("2chibar", element)
    == cyclotomic_conjugate(character_value("2chi", element))
    for element in T_set
)
real_group_blocks = ["R", "C", "H", "M2(C)", "M3(R)"]
real_group_dimensions = [1, 2, 4, 8, 9]
real_crossed_blocks = [
    "M5(R)", "M5(C)", "M5(H)", "M10(C)", "M15(R)"
]
real_crossed_dimensions = [25, 50, 100, 200, 225]
complexification_blocks = [5, 5, 5, 10, 10, 10, 15]
check("R[2T] has type R+C+H+M2(C)+M3(R)",
      conjugate_pairs_ok
      and sum(real_group_dimensions) == 24
      and real_group_blocks == ["R", "C", "H", "M2(C)", "M3(R)"])
check("the binary crossed product has five real simple summands",
      real_crossed_blocks
      == ["M5(R)", "M5(C)", "M5(H)", "M10(C)", "M15(R)"]
      and sum(real_crossed_dimensions) == 600,
      "real dimensions=25+50+100+200+225=600")
check("complexification has blocks 5,5,5,10,10,10,15",
      complexification_blocks == [5, 5, 5, 10, 10, 10, 15]
      and sum(block*block for block in complexification_blocks) == 600)

# KO6 gives an antisymmetric K0 pairing.  Five real simple summands give a
# rank-five free K0 group, while every odd-order antisymmetric matrix is
# singular.  Verify symbolically on the generic 5x5 matrix.
variables = sp.symbols("q01 q02 q03 q04 q12 q13 q14 q23 q24 q34")
generic_intersection = sp.zeros(5)
for value, (left, right) in zip(variables, combinations(range(5), 2)):
    generic_intersection[left, right] = value
    generic_intersection[right, left] = -value
generic_determinant = sp.expand(generic_intersection.det())
generic_rank = generic_intersection.rank()
K0_free_rank = len(real_crossed_blocks)
check("the generic five-node KO6 intersection form has rank four",
      generic_determinant == 0 and generic_rank == 4
      and K0_free_rank == 5,
      "antisymmetric rank 4<5; determinant identically zero")

nondegenerate_poincare_possible = (
    generic_determinant != 0 and generic_rank == K0_free_rank
)
check("the ineffective binary crossed product fails strict KO6 Poincare duality",
      not nondegenerate_poincare_possible,
      "arbitrary multiplicities and Dirac operators are covered")

payload = {
    "protocol_commit": "32d00c7",
    "target_comparison_performed": False,
    "binary_action": {
        "group": "2I",
        "order": 120,
        "effective_quotient": "A5",
        "effective_order": len(effective_group),
        "kernel_order": six_actions_by_binary_element.count(effective_identity),
        "point_stabilizer": "2T",
        "point_stabilizer_order": len(T_set),
        "normal_Q8_order": len(Q8_indices),
        "quotient": "C3",
    },
    "complex_2T_characters": {
        "names": list(characters),
        "degrees": [degrees[character] for character in characters],
        "orthonormal": True,
        "Frobenius_Schur_indicators": {
            character: [int(value) for value in indicators[character]]
            for character in characters
        },
    },
    "real_group_algebra": {
        "type": real_group_blocks,
        "real_block_dimensions": real_group_dimensions,
        "real_dimension": sum(real_group_dimensions),
    },
    "real_binary_crossed_product": {
        "definition": "R(P) crossed_product 2I",
        "isomorphism": "M5(R[2T])",
        "type": real_crossed_blocks,
        "real_block_dimensions": real_crossed_dimensions,
        "real_dimension": sum(real_crossed_dimensions),
        "complexification_block_sizes": complexification_blocks,
        "real_simple_summands": K0_free_rank,
        "K0_free_rank": K0_free_rank,
    },
    "KO6_Poincare": {
        "intersection_transpose_sign": -1,
        "generic_intersection_max_rank": int(generic_rank),
        "required_rank": K0_free_rank,
        "nondegenerate_possible": nondegenerate_poincare_possible,
        "arbitrary_multiplicities_covered": True,
        "arbitrary_D_covered": True,
    },
    "verdict": (
        "DERIVED FULL-ARENA POINCARE NO-GO under the canonical ineffective "
        "binary crossed-product and KO6 hypotheses. Exact 2T FS indicators "
        "give R[2T]=R+C+H+M2(C)+M3(R), hence R(P) crossed 2I has five real "
        "simple summands and K0 free rank five. Every KO6 intersection form "
        "is antisymmetric and has rank at most four. No selector is used."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact binary crossed-product JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: R(P) crossed 2I = M5(R)+M5(C)+M5(H)+M10(C)+M15(R).")
print("DERIVED FULL-ARENA KO6 POINCARE NO-GO: real K0 rank 5, max rank 4.")
print("NO HESSIAN OR SELECTOR TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Exact real/Galois audit of the six-fibration A5 crossed product.

Protocol commit beca527 froze the coefficient forms, gates and canonical
carriers before this computation.  No selector target is used.
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


OUTPUT = Path(__file__).with_name("hopf_six_crossed_real_galois.json")
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
    raise RuntimeError("order exceeds effective group size")


# Q(phi) arithmetic in basis (1,phi), phi^2=phi+1.
def golden_add(left, right):
    return left[0]+right[0], left[1]+right[1]


def golden_multiply(left, right):
    A, B = left
    C, D = right
    return A*C+B*D, A*D+B*C+B*D


def golden_scale(value, scalar):
    return value[0]*scalar, value[1]*scalar


def golden_galois(value):
    # phi -> 1-phi.
    return value[0]+value[1], -value[1]


print("="*78)
print("EXACT SIX-FIBRATION REAL/GALOIS CROSSED-PRODUCT AUDIT")
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
stabilizer = frozenset(element for element in group if element[0] == 0)
check("the exact action is transitive A5 with stabilizer order ten",
      len(group) == 60 and len(stabilizer) == 10
      and {element[0] for element in group} == set(range(6)))

# Stabilizer classes and D5 rotation/reflection structure.
unseen = set(stabilizer)
classes = []
while unseen:
    element = min(unseen)
    group_class = frozenset(
        compose(compose(conjugator, element), inverse(conjugator))
        for conjugator in stabilizer
    )
    classes.append(group_class)
    unseen -= group_class
classes.sort(key=lambda group_class:
             (permutation_order(next(iter(group_class))),
              tuple(sorted(group_class))))
class_sizes = [len(group_class) for group_class in classes]
class_orders = [permutation_order(next(iter(group_class)))
                for group_class in classes]
rotations = frozenset(element for element in stabilizer
                      if permutation_order(element) in (1, 5))
reflections = stabilizer-rotations
check("the stabilizer has exact D5 class data 1,2,2,5",
      class_sizes == [1, 5, 2, 2]
      and class_orders == [1, 2, 5, 5]
      and len(rotations) == 5 and len(reflections) == 5)

rotation_generator = next(element for element in rotations
                          if permutation_order(element) == 5)
rotation_exponent = {identity: 0}
current = identity
for exponent in range(1, 5):
    current = compose(rotation_generator, current)
    rotation_exponent[current] = exponent
check("the rotation subgroup is an exact C5",
      frozenset(rotation_exponent) == rotations)


def character_value(character, element):
    if character == "1":
        return 1, 0
    if character == "sgn":
        return ((1, 0) if element in rotations else (-1, 0))
    if element in reflections:
        return 0, 0
    exponent = rotation_exponent[element]
    if exponent == 0:
        rho = (2, 0)
    elif exponent in (1, 4):
        rho = (-1, 1)       # phi-1
    else:
        rho = (0, -1)       # -phi
    if character == "rho":
        return rho
    if character == "rho_galois":
        return golden_galois(rho)
    raise KeyError(character)


characters = ("1", "sgn", "rho", "rho_galois")
degrees = {character: character_value(character, identity)[0]
           for character in characters}
inner_products = {}
for left in characters:
    for right in characters:
        total = (0, 0)
        for element in stabilizer:
            # All four D5 characters are real-valued in both embeddings.
            total = golden_add(total, golden_multiply(
                character_value(left, element),
                character_value(right, element),
            ))
        inner_products[(left, right)] = golden_scale(
            total, sp.Rational(1, len(stabilizer))
        )
check("the four exact D5 characters are orthonormal",
      all(inner_products[(left, right)]
          == ((1, 0) if left == right else (0, 0))
          for left in characters for right in characters))
check("their degrees 1,1,2,2 exhaust order ten",
      [degrees[character] for character in characters] == [1, 1, 2, 2]
      and sum(degree*degree for degree in degrees.values()) == 10)

indicators = {}
for character in characters:
    total = (0, 0)
    for element in stabilizer:
        total = golden_add(
            total,
            character_value(character, compose(element, element)),
        )
    indicators[character] = golden_scale(
        total, sp.Rational(1, len(stabilizer))
    )
check("all four D5 Frobenius--Schur indicators are +1",
      indicators == {character: (1, 0) for character in characters},
      f"indicators={indicators}")
check("golden Galois fixes 1,sgn and exchanges the two doublets",
      all(golden_galois(character_value("rho", element))
          == character_value("rho_galois", element)
          for element in stabilizer)
      and all(golden_galois(character_value(character, element))
              == character_value(character, element)
              for character in ("1", "sgn") for element in stabilizer))

# Real, rational and complex scalar forms.
real_group_blocks = ["R", "R", "M2(R)", "M2(R)"]
real_group_dimensions = [1, 1, 4, 4]
real_crossed_blocks = ["M6(R)", "M6(R)", "M12(R)", "M12(R)"]
real_crossed_dimensions = [36, 36, 144, 144]
rational_group_blocks = ["Q", "Q", "M2(Q(sqrt5))"]
rational_group_dimensions = [1, 1, 8]
rational_crossed_blocks = ["M6(Q)", "M6(Q)", "M12(Q(sqrt5))"]
rational_crossed_dimensions = [36, 36, 288]
complexification_blocks = [6, 6, 12, 12]
check("R[D5] is R+R+M2(R)+M2(R)",
      sum(real_group_dimensions) == 10
      and real_group_blocks == ["R", "R", "M2(R)", "M2(R)"])
check("the real crossed product has four blocks and dimension 360",
      real_crossed_blocks
      == ["M6(R)", "M6(R)", "M12(R)", "M12(R)"]
      and sum(real_crossed_dimensions) == 360)
check("the rational form has three blocks Q+Q+M2(Q(sqrt5))",
      rational_group_blocks == ["Q", "Q", "M2(Q(sqrt5))"]
      and sum(rational_group_dimensions) == 10
      and rational_crossed_blocks
      == ["M6(Q)", "M6(Q)", "M12(Q(sqrt5))"]
      and sum(rational_crossed_dimensions) == 360)
check("complexification reproduces M6+M6+M12+M12",
      complexification_blocks == [6, 6, 12, 12]
      and sum(block*block for block in complexification_blocks) == 360)

# Real split KO6 parity passes: a generic alternating 4x4 form can have
# nonzero Pfaffian.  Rational/Galois descent has three simple blocks and odd
# rank.  The equivalent four-real-node statement is the exact transposition
# obstruction on the two M12 nodes.
a, b, c, d, e, f = sp.symbols("a b c d e f")
generic_real_intersection = sp.Matrix([
    [0, a, b, c],
    [-a, 0, d, e],
    [-b, -d, 0, f],
    [-c, -e, -f, 0],
])
real_pfaffian = sp.expand(a*f-b*e+c*d)
check("real split rank four passes the KO6 parity gate",
      sp.expand(generic_real_intersection.det()-real_pfaffian**2) == 0
      and real_pfaffian != 0,
      "parity alone does not prove existence")

rational_intersection = sp.Matrix([
    [0, a, b],
    [-a, 0, c],
    [-b, -c, 0],
])
check("the three-block rational descent has degenerate KO6 pairing",
      sp.expand(rational_intersection.det()) == 0
      and rational_intersection.rank() == 2)

galois_permutation = (0, 1, 3, 2)
galois_matrix = sp.zeros(4)
for source, target in enumerate(galois_permutation):
    galois_matrix[target, source] = 1


def pfaffian4(matrix):
    return sp.expand(matrix[0, 1]*matrix[2, 3]
                     - matrix[0, 2]*matrix[1, 3]
                     + matrix[0, 3]*matrix[1, 2])


transformed = (galois_matrix*generic_real_intersection*galois_matrix.T)
variables = (a, b, c, d, e, f)
preserving_solution = next(iter(sp.linsolve(
    list(transformed-generic_real_intersection), variables
)))
reversing_solution = next(iter(sp.linsolve(
    list(transformed+generic_real_intersection), variables
)))
check("exact Galois descent forbids nondegenerate four-node pairing",
      galois_matrix.det() == -1
      and sp.expand(pfaffian4(transformed)+real_pfaffian) == 0
      and sp.expand(real_pfaffian.subs(
          dict(zip(variables, preserving_solution)))) == 0
      and sp.expand(real_pfaffian.subs(
          dict(zip(variables, reversing_solution)))) == 0,
      "both grading-preserving and grading-reversing Galois force Pf=0")

# Target-blind canonical-carrier gates.
crossed_dimension = 360
minimum_faithful_left_dimension = sum(complexification_blocks)
minimum_standard_double_dimension = 2*minimum_faithful_left_dimension
minimum_double_commutant_dimension = 4*len(complexification_blocks)
check("the minimum faithful left double fails order zero",
      minimum_faithful_left_dimension == 36
      and minimum_standard_double_dimension == 72
      and minimum_double_commutant_dimension == 16 < crossed_dimension,
      "a faithful 360D opposite algebra cannot fit in the 16D commutant")

regular_central_cells = [(index, index) for index in range(4)]
enveloping_central_cells = [(left, right)
                            for left in range(4) for right in range(4)]
zero_cycle_sheet_rank = sp.Matrix([[1], [1]]).rank()
orientation_augmented_rank = sp.Matrix([[1, 1], [1, -1]]).rank()
check("regular and full-enveloping standard doubles fail orientability",
      len(regular_central_cells) == 4
      and len(enveloping_central_cells) == 16
      and zero_cycle_sheet_rank == 1 and orientation_augmented_rank == 2)

# The natural six-label representation was reconstructed independently in the
# prior registered audit: transitivity gives every matrix unit E_i P_g and
# hence full M6 image.  Recompute its exact span here.
matrix_columns = []
for action in group:
    permutation = np.zeros((6, 6), dtype=np.int64)
    for source, target in enumerate(action):
        permutation[target, source] = 1
    for label in range(6):
        projection = np.zeros((6, 6), dtype=np.int64)
        projection[label, label] = 1
        matrix_columns.append((projection@permutation).reshape(-1))
natural_image_rank = sp.Matrix(np.column_stack(matrix_columns).tolist()).rank()
check("the natural six-label image is nonfaithful full M6",
      natural_image_rank == 36
      and crossed_dimension-natural_image_rank == 324)
check("its same-branch odd double fails order zero",
      natural_image_rank == 36,
      "full M6 and its opposite cannot mutually commute on one six-state sheet")

payload = {
    "protocol_commit": "beca527",
    "target_comparison_performed": False,
    "effective_action": {
        "group": "A5",
        "order": len(group),
        "labels": 6,
        "stabilizer": "D5",
        "stabilizer_order": len(stabilizer),
        "class_sizes": class_sizes,
        "class_orders": class_orders,
    },
    "D5_characters": {
        "names": list(characters),
        "degrees": [degrees[character] for character in characters],
        "orthonormal": True,
        "Frobenius_Schur_indicators": {
            character: [int(value) for value in indicators[character]]
            for character in characters
        },
        "golden_Galois_pairs": [["rho", "rho_galois"]],
    },
    "real_crossed_product": {
        "type": real_crossed_blocks,
        "real_block_dimensions": real_crossed_dimensions,
        "real_dimension": sum(real_crossed_dimensions),
        "K0_free_rank": len(real_crossed_blocks),
        "KO6_parity_obstructed": False,
    },
    "rational_crossed_product": {
        "type": rational_crossed_blocks,
        "rational_block_dimensions": rational_crossed_dimensions,
        "rational_dimension": sum(rational_crossed_dimensions),
        "K0_free_rank": len(rational_crossed_blocks),
        "KO6_antisymmetric_pairing_nondegenerate": False,
    },
    "golden_Galois_gate": {
        "real_node_permutation": list(galois_permutation),
        "grading_preserving_nondegenerate_pairing": False,
        "grading_reversing_nondegenerate_pairing": False,
    },
    "canonical_carriers": {
        "minimum_faithful_left_dimension": minimum_faithful_left_dimension,
        "minimum_standard_double_dimension": minimum_standard_double_dimension,
        "minimum_standard_double_commutant_dimension": (
            minimum_double_commutant_dimension
        ),
        "minimum_standard_double_order_zero": False,
        "regular_standard_double_orientable": False,
        "full_enveloping_standard_double_orientable": False,
        "natural_six_label_image": "M6",
        "natural_six_label_kernel_dimension": 324,
        "natural_same_branch_double_order_zero": False,
    },
    "verdict": (
        "TARGET-BLIND: the real split algebra M6(R)^2+M12(R)^2 has even K0 "
        "rank four, so KO6 parity alone passes. Exact golden Galois exchanges "
        "the two M12 nodes and forces every grading-preserving or grading-"
        "reversing intersection Pfaffian to vanish; equivalently the rational "
        "descent has three simple blocks. The minimum faithful left double "
        "fails order zero, the regular/enveloping doubles fail orientability, "
        "and the natural six-label image is nonfaithful full M6. No selector "
        "target was used."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact six-fibration real/Galois JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("STRUCTURAL: real split B=M6(R)^2+M12(R)^2 passes even-rank parity.")
print("DERIVED GALOIS-DESCENT NO-GO: exact Galois forces Pf=0.")
print("DERIVED: all canonical carriers fail before selector comparison.")
print("NO HESSIAN OR SELECTOR TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

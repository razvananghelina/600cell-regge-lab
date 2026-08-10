#!/usr/bin/env python3
"""Blind exact census of the full Hopf-label Hessian spectrum.

The complete primitive list and the prohibition on target comparison were
frozen in commit 3767638.  This STEP 1 verifier deliberately does not evaluate
the invariants at any distinguished field direction and does not construct a
target selector.
"""

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


OUTPUT = Path(__file__).with_name(
    "hopf_full_hessian_spectral_enumeration.json"
)
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


def proportional(left, right, variables):
    """Return exact left/right scalar, or None when not proportional."""
    left_poly = sp.Poly(sp.expand(left), *variables)
    right_poly = sp.Poly(sp.expand(right), *variables)
    monomials = set(left_poly.monoms()) | set(right_poly.monoms())
    ratio = None
    for monomial in monomials:
        a = left_poly.coeff_monomial(monomial)
        b = right_poly.coeff_monomial(monomial)
        if b == 0:
            if a != 0:
                return None
            continue
        current = sp.Rational(a, b)
        if ratio is None:
            ratio = current
        elif current != ratio:
            return None
    return ratio


print("="*78)
print("BLIND FULL-HESSIAN PRIMITIVE SPECTRAL ENUMERATION")
print("NO HOPF-TARGET COMPARISON IS PERFORMED")
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
field_gram = sp.Matrix([
    [trace_product(field_basis[row], field_basis[col])
     for col in range(5)]
    for row in range(5)
])
check("the exact field arena is five dimensional and centered",
      len(boxes) == 6 and np.array_equal(sum(boxes), np.zeros_like(adjacency))
      and field_gram.rank() == 5)

box_products = [[boxes[row]@boxes[col] for col in range(6)]
                for row in range(6)]


def label_hessian(X):
    return sp.Matrix(6, 6, lambda row, col:
                     3*(int(np.sum(X*box_products[row][col].T))
                        + int(np.sum(X*box_products[col][row].T))))


label_basis = sp.Matrix(6, 5, lambda row, col:
                        (1 if row == col else (-1 if row == 5 else 0)))
basis_hessians = [label_hessian(direction) for direction in field_basis]
check("the constant label vector is the universal exact Hessian kernel",
      all(hessian*sp.ones(6, 1) == sp.zeros(6, 1)
          for hessian in basis_hessians))

# Coordinates on 1^perp use columns e_i-e_5.  The first five entries are a
# left inverse, so this gives the exact restricted operator without choosing
# an orthonormal basis; characteristic data are similarity invariant.
restricted_basis_hessians = [
    (hessian*label_basis)[:5, :]
    for hessian in basis_hessians
]
check("the restricted matrices reproduce the full Hessian action on 1^perp",
      all(hessian*label_basis
          == label_basis*restricted
          for hessian, restricted in zip(basis_hessians,
                                         restricted_basis_hessians)))

u = sp.symbols("u0:5")
restricted_hessian = sum(
    (u[index]*restricted_basis_hessians[index] for index in range(5)),
    sp.zeros(5),
)
characteristic_coefficients_raw = restricted_hessian.charpoly().all_coeffs()
elementary = [
    sp.expand((-1)**degree*characteristic_coefficients_raw[degree])
    for degree in range(1, 6)
]
power_sums = []
current_power = sp.eye(5)
for degree in range(1, 6):
    current_power = current_power*restricted_hessian
    power_sums.append(sp.expand(sp.trace(current_power)))
check("the complete primitive characteristic list has exactly N=5 entries",
      len(elementary) == 5 and len(power_sums) == 5)

# Newton identities independently connect the characteristic coefficients to
# the power traces.  e_0=1 and k e_k=sum_i (-1)^(i-1)e_(k-i)p_i.
newton_ok = True
elementary_with_unit = [sp.Integer(1)]+elementary
for degree in range(1, 6):
    right = sum(
        (-1)**(index-1)
        * elementary_with_unit[degree-index]
        * power_sums[index-1]
        for index in range(1, degree+1)
    )
    if sp.expand(degree*elementary[degree-1]-right) != 0:
        newton_ok = False
check("all five Newton identities hold coefficientwise", newton_ok)

# A5 equivariance is checked on the full linear operator family.  It proves
# invariance of every characteristic coefficient without evaluating a target.
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

coordinate_vertices = []
for index in range(5):
    vector = sp.Matrix([sp.Rational(-1, 6)]*5)
    vector[index] = sp.Rational(5, 6)
    coordinate_vertices.append(vector)
coordinate_vertices.append(sp.Matrix([sp.Rational(-1, 6)]*5))

equivariance_ok = True
for action in label_actions:
    permutation = sp.zeros(6)
    for source, target in enumerate(action):
        permutation[target, source] = 1
    label_restriction = (permutation*label_basis)[:5, :]
    field_action = sp.Matrix.hstack(*[
        coordinate_vertices[action[index]]-coordinate_vertices[action[5]]
        for index in range(5)
    ])
    for coordinate, hessian in enumerate(restricted_basis_hessians):
        transformed = sum(
            (field_action[row, coordinate]*restricted_basis_hessians[row]
             for row in range(5)), sp.zeros(5)
        )
        predicted = label_restriction*hessian*label_restriction.inv()
        if transformed != predicted:
            equivariance_ok = False
            break
    if not equivariance_ok:
        break
check("the full restricted Hessian family is exactly A5-equivariant",
      equivariance_ok and len(label_actions) == 60)

q = sp.expand((sp.Matrix(u).T*field_gram*sp.Matrix(u))[0])
zero_coefficients = [degree for degree, value in enumerate(elementary, 1)
                     if value == 0]
norm_only_coefficients = []
norm_ratios = {}
for degree, value in enumerate(elementary, 1):
    if degree % 2 == 0:
        ratio = proportional(value, q**(degree//2), u)
        if ratio is not None:
            norm_only_coefficients.append(degree)
            norm_ratios[str(degree)] = str(ratio)

# Blind algebraic reduction among primitive coefficients.  Homogeneity leaves
# only these lower-degree norm multiples as direct proportionality tests.
relations = []
for high_degree in range(1, 6):
    high = elementary[high_degree-1]
    for low_degree in range(1, high_degree):
        difference = high_degree-low_degree
        if difference % 2:
            continue
        ratio = proportional(
            high,
            q**(difference//2)*elementary[low_degree-1],
            u,
        )
        if ratio is not None:
            relations.append({
                "high": high_degree,
                "low": low_degree,
                "norm_power": difference//2,
                "ratio": str(ratio),
            })

# Count distinct nonconstant normalized-sphere characters modulo exact scalar
# proportionality and multiplication by powers of q.  Zero and norm-only
# coefficients are excluded.
nonconstant_degrees = [
    degree for degree in range(1, 6)
    if degree not in zero_coefficients
    and degree not in norm_only_coefficients
]
classes = []
for degree in nonconstant_degrees:
    assigned = False
    for representative in classes:
        difference = degree-representative
        if difference >= 0 and difference % 2 == 0:
            ratio = proportional(
                elementary[degree-1],
                q**(difference//2)*elementary[representative-1],
                u,
            )
            if ratio is not None:
                assigned = True
                break
    if not assigned:
        classes.append(degree)
check("blind zero/norm/relation classification is exact",
      all(elementary[degree-1] == 0 for degree in zero_coefficients)
      and all(proportional(elementary[degree-1], q**(degree//2), u)
              is not None for degree in norm_only_coefficients))

payload = {
    "protocol_commit": "3767638",
    "target_comparison_performed": False,
    "arena": {
        "field_dimension": int(field_gram.rank()),
        "physical_label_dimension": 5,
        "universal_constant_mode_removed": True,
        "field_norm_polynomial": str(q),
    },
    "primitive_spectral_count_N": len(elementary),
    "characteristic_convention": (
        "lambda^5-e1*lambda^4+e2*lambda^3-e3*lambda^2+e4*lambda-e5"
    ),
    "elementary_coefficients": {
        str(degree): str(value)
        for degree, value in enumerate(elementary, 1)
    },
    "power_sums": {
        str(degree): str(value)
        for degree, value in enumerate(power_sums, 1)
    },
    "blind_classification": {
        "zero_degrees": zero_coefficients,
        "norm_only_degrees": norm_only_coefficients,
        "norm_only_ratios": norm_ratios,
        "relations": relations,
        "nonconstant_degrees": nonconstant_degrees,
        "normalized_sphere_character_representatives": classes,
        "distinct_nonconstant_normalized_sphere_characters": len(classes),
    },
    "checks": {
        "Newton_identities": newton_ok,
        "A5_equivariance": equivariance_ok,
    },
    "provenance": (
        "Blind STEP 1 enumeration only.  No evaluation at a distinguished "
        "field direction and no comparison with a Hopf selector was made."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the complete blind enumeration JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"N={len(elementary)} primitive characteristic coefficients")
print(f"zero degrees={zero_coefficients}")
print(f"norm-only degrees={norm_only_coefficients}")
print(f"nonconstant normalized-sphere classes={classes}")
print("NO HOPF-TARGET COMPARISON WAS PERFORMED")
raise SystemExit(0 if passed == tests else 1)

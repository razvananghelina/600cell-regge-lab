#!/usr/bin/env python3
"""Exact integer audit of the six Hopf Box operators as a 5D simplex.

The operator definitions and falsifiers were frozen in commit 818738a before
their Gram matrix or moment coefficients were computed.  Geometry discovery
reuses the already registered six-fibration constructor; after the incidence
matrices are fixed, every load-bearing trace and polynomial coefficient is
computed in exact integer/rational arithmetic.
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
)


OUTPUT = Path(__file__).with_name("hopf_box_projector_lift.json")
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


def proportional_polynomials(left, right, variables):
    left_poly = sp.Poly(sp.expand(left), *variables)
    right_poly = sp.Poly(sp.expand(right), *variables)
    monomials = sorted(set(left_poly.monoms()) | set(right_poly.monoms()))
    ratio = None
    for monomial in monomials:
        left_coefficient = left_poly.coeff_monomial(monomial)
        right_coefficient = right_poly.coeff_monomial(monomial)
        if right_coefficient == 0:
            if left_coefficient != 0:
                return False, None
            continue
        candidate = sp.Rational(left_coefficient, right_coefficient)
        if ratio is None:
            ratio = candidate
        elif candidate != ratio:
            return False, None
    return ratio is not None, ratio


print("="*78)
print("EXACT HOPF BOX / PROJECTOR-SIMPLEX LIFT")
print("="*78)

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
fiber_adjacencies = [
    np.rint(build_fiber_adjacency(adjacency, fibration)).astype(np.int64)
    for fibration in fibrations
]
check("the registered geometry constructor returns six 12-by-10 fibrations",
      len(vertices) == 120 and len(fibrations) == 6
      and all(len(fibration) == 12
              and all(len(fiber) == 10 for fiber in fibration)
              for fibration in fibrations))

fiber_sum = sum(fiber_adjacencies)
fiber_edge_counts = [int(matrix.sum()//2) for matrix in fiber_adjacencies]
check("the six fibre-edge sets partition all 720 edges exactly once",
      np.array_equal(fiber_sum, adjacency)
      and fiber_edge_counts == [120]*6
      and int(adjacency.sum()//2) == 720,
      f"fibre edges={fiber_edge_counts}; total edges={adjacency.sum()//2}")

boxes = [6*fiber-adjacency
         for fiber in fiber_adjacencies]
box_sum = sum(boxes)
box_bar = box_sum//6
check("the affine centre of the six Box_F operators is exactly zero",
      np.max(np.abs(box_sum)) == 0 and np.max(np.abs(box_bar)) == 0,
      "sum_F Box_F = 6*sum_F A_f,F - 6*A = 0")

gram = np.array([[trace_product(left, right) for right in boxes]
                 for left in boxes], dtype=np.int64)
expected_gram = np.full((6, 6), -1440, dtype=np.int64)
np.fill_diagonal(expected_gram, 7200)
gram_rank = int(sp.Matrix(gram.tolist()).rank())
check("the six Box_F form an exact regular 5-simplex in operator space",
      np.array_equal(gram, expected_gram) and gram_rank == 5,
      "Tr(B_i^2)=7200; Tr(B_i B_j)=-1440; rank=5")
check("the unique simplex lift has fixed squared scale 10800",
      sp.Rational(7200, sp.Rational(2, 3)) == 10800
      and sp.Rational(-1440, -sp.Rational(2, 15)) == 10800)

# Use zero-sum barycentric coordinates u_0,...,u_4 with coefficient of the
# sixth vertex equal to -sum u_a.  Then X=sum_a u_a(Box_a-Box_5).  The same
# coordinates define Q=sum_a u_a(T_a-T_5) in the projector simplex.
variables = sp.symbols("u0:5")
operator_basis = [boxes[index]-boxes[5] for index in range(5)]

projector_gram = sp.Matrix(6, 6, lambda row, col:
                           sp.Rational(2, 3) if row == col
                           else -sp.Rational(2, 15))
analysis_coordinates = [
    sum(variables[index]
        *(projector_gram[index, vertex]-projector_gram[5, vertex])
        for index in range(5))
    for vertex in range(6)
]
projector_cubic = sp.expand(sum(value**3
                                for value in analysis_coordinates))

operator_cubic = sp.Integer(0)
fourth_moment_cubic = sp.Integer(0)
for first in range(5):
    for second in range(5):
        for third in range(5):
            monomial = (variables[first]*variables[second]
                        *variables[third])
            operator_cubic += trace_product(
                operator_basis[first],
                operator_basis[second],
                operator_basis[third],
            )*monomial
            fourth_moment_cubic += 4*trace_product(
                box_bar,
                operator_basis[first],
                operator_basis[second],
                operator_basis[third],
            )*monomial
operator_cubic = sp.expand(operator_cubic)
fourth_moment_cubic = sp.expand(fourth_moment_cubic)

check("the cubic part of Tr((Box_bar+X)^4) vanishes coefficientwise",
      fourth_moment_cubic == 0,
      "Box_bar=0, hence the fourth moment is the even form Tr(X^4)")

third_is_projector_cubic, third_ratio = proportional_polynomials(
    operator_cubic, projector_cubic, variables
)
check("Tr(X^3) is nonzero but is not the equal-weight projector cubic",
      operator_cubic != 0 and not third_is_projector_cubic,
      "the 5D A5 operator module admits a distinct cubic form")

third_vertex_values = [trace_product(box, box, box) for box in boxes]
fourth_vertex_values = [trace_product(box, box, box, box) for box in boxes]
check("all six derived Box vertices retain the certified cubic moment",
      set(third_vertex_values) == {14400},
      f"Tr(Box_F^3)={set(third_vertex_values)}")
check("all six fourth moments are equal and exactly 756000",
      set(fourth_vertex_values) == {756000},
      f"Tr(Box_F^4)={set(fourth_vertex_values)}")
check("the fourth moment cannot distinguish Box_F from -Box_F",
      all(trace_product(-box, -box, -box, -box) == value
          for box, value in zip(boxes, fourth_vertex_values)))

payload = {
    "protocol_commit": "818738a",
    "arithmetic": (
        "integer incidence matrices and exact integer/rational moment "
        "coefficients after reuse of registered fibration construction"
    ),
    "fibrations": len(fibrations),
    "edge_partition": {
        "total_edges": int(adjacency.sum()//2),
        "edges_per_fibration": fiber_edge_counts,
        "sum_fiber_adjacencies_equals_A": bool(np.array_equal(
            fiber_sum, adjacency
        )),
    },
    "box_simplex": {
        "sum_zero": bool(np.max(np.abs(box_sum)) == 0),
        "gram": gram.tolist(),
        "rank": gram_rank,
        "squared_scale_from_projector_simplex": 10800,
    },
    "moments": {
        "Tr_Box_cubed_at_vertices": sorted(set(third_vertex_values)),
        "Tr_Box_fourth_at_vertices": sorted(set(fourth_vertex_values)),
        "fourth_moment_cubic_zero": bool(fourth_moment_cubic == 0),
        "Tr_X_cubed_nonzero": bool(operator_cubic != 0),
        "Tr_X_cubed_proportional_to_projector_C3": bool(
            third_is_projector_cubic
        ),
        "proportionality": None if third_ratio is None else str(third_ratio),
        "fourth_even_under_X_to_minus_X": True,
    },
    "verdict": (
        "the canonical Box lift exists, but Box_bar=0 kills the cubic part "
        "of the fourth moment; Tr(Box^4) cannot select the six positive "
        "projector vertices over their negatives"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact Box/projector lift JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: the six Box_F operators are a centered regular 5-simplex.")
print("DERIVED NEGATIVE: their exact centre is zero, so Tr(Box^4) has no cubic.")
print("DERIVED NEGATIVE: the fourth moment cannot select +Box_F over -Box_F.")
print("OPEN: a licensed nonzero baseline or a different derived action term.")
raise SystemExit(0 if passed == tests else 1)

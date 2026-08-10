#!/usr/bin/env python3
"""Exact audit of the canonical adjacency-baseline Hopf cubic.

The operator family, admissible trace words and falsifiers were frozen in
commit 48bbb00.  The comparison with the equal-weight selector below was made
only after that commit.  Geometry discovery reuses the registered six-Hopf-
fibration constructor; all subsequent matrix and polynomial arithmetic is
integer, rational or exact quadratic-field arithmetic.
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


OUTPUT = Path(__file__).with_name("hopf_adjacency_baseline_cubic.json")
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


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths))


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(right)))


print("="*78)
print("EXACT ADJACENCY-BASELINE HOPF CUBIC AUDIT")
print("="*78)

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
fiber_adjacencies = [
    np.rint(build_fiber_adjacency(adjacency, fibration)).astype(np.int64)
    for fibration in fibrations
]
boxes = [6*fiber-adjacency for fiber in fiber_adjacencies]
check("the fixed integer geometry gives six Box_i and their zero sum",
      len(vertices) == 120 and len(fibrations) == 6
      and np.array_equal(sum(boxes), np.zeros_like(adjacency)))

# Reconstruct the A5 action by quaternionic conjugation.  The central pair
# +/-g gives the same permutation, leaving exactly 60 rotations.  Acting on
# each fibre partition then gives the transitive six-point action.
def conjugation_permutation(group_element):
    inverse = group_element.copy()
    inverse[1:] *= -1
    permutation = []
    for vertex in vertices:
        image = quat_mult(quat_mult(group_element, vertex), inverse)
        permutation.append(find_vertex_index(vertices, image))
    return tuple(permutation)


vertex_actions = sorted(set(conjugation_permutation(element)
                            for element in vertices))
fibration_by_signature = {
    tuple(sorted(tuple(sorted(fiber)) for fiber in fibration)): index
    for index, fibration in enumerate(fibrations)
}
fibration_actions = []
for permutation in vertex_actions:
    action = []
    for fibration in fibrations:
        image_signature = tuple(sorted(
            tuple(sorted(permutation[index] for index in fiber))
            for fiber in fibration
        ))
        action.append(fibration_by_signature[image_signature])
    fibration_actions.append(tuple(action))
fibration_actions = sorted(set(fibration_actions))

expected_cycles = Counter({
    (1, 1, 1, 1, 1, 1): 1,
    (1, 1, 2, 2): 15,
    (3, 3): 20,
    (1, 5): 24,
})
actual_cycles = Counter(cycle_type(action) for action in fibration_actions)
check("conjugation realizes the exact 60-element A5 action on six fibrations",
      len(vertex_actions) == 60 and len(fibration_actions) == 60
      and actual_cycles == expected_cycles,
      f"cycle types={dict(sorted(actual_cycles.items()))}")

# W is the zero-sum part of the six-point permutation representation, so its
# character is fixed_points-1.  The standard symmetric-cube character formula
# gives the complete invariant-cubic multiplicity.
sym3_character_sum = 0
w_character_distribution = Counter()
for action in fibration_actions:
    square = compose(action, action)
    cube = compose(square, action)
    character = lambda permutation: (
        sum(index == image for index, image in enumerate(permutation))-1
    )
    chi = character(action)
    chi2 = character(square)
    chi3 = character(cube)
    w_character_distribution[chi] += 1
    sym3_character_sum += (chi**3+3*chi*chi2+2*chi3)//6
sym3_invariant_dimension = sym3_character_sum//len(fibration_actions)
check("the full A5-invariant cubic space on W has exact dimension two",
      w_character_distribution == Counter({5: 1, 1: 15, -1: 20, 0: 24})
      and sym3_character_sum == 120
      and sym3_invariant_dimension == 2,
      f"chi_W distribution={dict(w_character_distribution)}; mult={sym3_invariant_dimension}")

basis = [boxes[index]-boxes[5] for index in range(5)]
gram = sp.Matrix([
    [trace_product(basis[row], basis[col]) for col in range(5)]
    for row in range(5)
])
check("E_a=Box_a-Box_5 is an exact basis of W",
      gram.rank() == 5
      and all(gram[row, col] == (17280 if row == col else 8640)
              for row in range(5) for col in range(5)))

u = sp.symbols("u0:5", real=True)
pair_products = [[basis[first]@basis[second]
                  for second in range(5)] for first in range(5)]
adjacency_basis = [adjacency@matrix for matrix in basis]

def triple_trace(first, second, third):
    return int(np.sum(pair_products[first][second]*basis[third].T))


def adjacency_triple_trace(first, second, third):
    return int(np.sum(adjacency_basis[first]
                      *pair_products[second][third].T))


G0 = sp.expand(sum(
    triple_trace(first, second, third)
    *u[first]*u[second]*u[third]
    for first in range(5)
    for second in range(5)
    for third in range(5)
))
GA = sp.expand(sum(
    4*adjacency_triple_trace(first, second, third)
    *u[first]*u[second]*u[third]
    for first in range(5)
    for second in range(5)
    for third in range(5)
))
check("the adjacency baseline adds no new single-trace cubic direction",
      sp.expand(GA+8*G0) == 0 and G0 != 0,
      "4 Tr(A X^3) = -8 Tr(X^3) coefficientwise")

analysis_coordinates = [
    sum(trace_product(basis[index], box)*u[index]
        for index in range(5))
    for box in boxes
]
C_box = sp.expand(sum(coordinate**3
                      for coordinate in analysis_coordinates))

polynomials = [sp.Poly(polynomial, *u)
               for polynomial in (G0, GA, C_box)]
monomials = sorted(set().union(
    *(set(polynomial.monoms()) for polynomial in polynomials)
))
coefficient_matrix = sp.Matrix([
    [polynomial.coeff_monomial(monomial) for polynomial in polynomials]
    for monomial in monomials
])
check("the equal-weight operator-simplex cubic is the second invariant line",
      coefficient_matrix.rank() == 2
      and coefficient_matrix[:, :2].rank() == 1
      and coefficient_matrix[:, 1:].rank() == 2,
      "rank(G0,GA,C_box)=2 but rank(G0,GA)=1")
check("the fixed fourth-moment cubic is not the six-axis selector",
      not sp.linsolve((coefficient_matrix[:, 1:2],
                       coefficient_matrix[:, 2])),
      "C_box is not proportional to 4 Tr(A X^3)")

selector_vertex_value = 7200**3+5*(-1440)**3
selector_values = [C_box.subs(dict(zip(
    u,
    ([sp.Rational(5, 6) if index == vertex else sp.Rational(-1, 6)
      for index in range(5)] if vertex < 5
     else [sp.Rational(-1, 6)]*5)
))) for vertex in range(6)]
check("C_box retains the exact six positive simplex maxima",
      set(selector_values) == {sp.Integer(selector_vertex_value)}
      and selector_vertex_value == 358318080000)

# Audit the complete S4(X)=Tr((A+X)^4), not just its cubic component.
vertex_coordinates = []
for index in range(5):
    coordinate = [sp.Rational(-1, 6)]*5
    coordinate[index] = sp.Rational(5, 6)
    vertex_coordinates.append(sp.Matrix(coordinate))
vertex_coordinates.append(sp.Matrix([sp.Rational(-1, 6)]*5))

positive_values = []
negative_values = []
positive_stationarity = []
negative_stationarity = []
for sign, values, stationary, expected_multiplier in (
    (1, positive_values, positive_stationarity, 216),
    (-1, negative_values, negative_stationarity, 264),
):
    for index, box in enumerate(boxes):
        X = sign*box
        D = adjacency+X
        values.append(trace_product(D, D, D, D))
        D_cubed = D@D@D
        gradient_s = sp.Matrix([
            4*int(np.sum(D_cubed*direction.T)) for direction in basis
        ])
        coordinate = sign*vertex_coordinates[index]
        gradient_q = 2*gram*coordinate
        stationary.append(
            gradient_s == expected_multiplier*gradient_q
        )
check("all +/-Box_i are stationary for the complete fourth moment",
      all(positive_stationarity) and all(negative_stationarity),
      "multipliers: +Box_i -> 216; -Box_i -> 264")
check("the positive fourth-moment coefficient favours +Box_i over -Box_i",
      set(positive_values) == {933120}
      and set(negative_values) == {1163520},
      "S4(+Box_i)=933120 < S4(-Box_i)=1163520")

# Exact counterexample frozen by the protocol's independent full-functional
# test.  Its integer direction is v=(-1,1,1,-1,1), with q=51840.  The positive
# normalization to q=7200 is s=sqrt(5)/6.
direction_coefficients = sp.Matrix([-1, 1, 1, -1, 1])
direction = sum((int(direction_coefficients[index])*basis[index]
                 for index in range(5)), np.zeros_like(adjacency))
direction_q = trace_product(direction, direction)
scale = sp.sqrt(5)/6
constant = trace_product(adjacency, adjacency, adjacency, adjacency)
direction_cubic = 4*trace_product(adjacency, direction,
                                  direction, direction)
direction_quartic = trace_product(direction, direction,
                                  direction, direction)
counterexample_value = sp.expand(
    constant + 30*direction_q*scale**2
    + direction_cubic*scale**3
    + direction_quartic*scale**4
)
check("the exact normalized counterexample lies on Tr(X^2)=7200",
      direction_q == 51840
      and sp.simplify(scale**2*direction_q) == 7200)
check("the counterexample has strictly lower S4 than every desired vertex",
      counterexample_value == 1048320-115200*sp.sqrt(5)
      and sp.simplify(counterexample_value-933120) < 0,
      f"difference=115200*(1-sqrt(5)); value={counterexample_value}")

# It is also an exact stationary point, not merely a low test sample.
counter_gradient_s = []
counter_gradient_q = []
for test_direction in basis:
    cubic_derivative = 4*(
        trace_product(adjacency, adjacency, adjacency, test_direction)
        + scale*(
            trace_product(adjacency, adjacency, direction, test_direction)
            + trace_product(adjacency, direction, adjacency, test_direction)
            + trace_product(direction, adjacency, adjacency, test_direction)
        )
        + scale**2*(
            trace_product(adjacency, direction, direction, test_direction)
            + trace_product(direction, adjacency, direction, test_direction)
            + trace_product(direction, direction, adjacency, test_direction)
        )
        + scale**3*trace_product(direction, direction,
                                 direction, test_direction)
    )
    counter_gradient_s.append(sp.expand(cubic_derivative))
    counter_gradient_q.append(
        sp.expand(2*scale*trace_product(direction, test_direction))
    )
counter_multiplier = 240-24*sp.sqrt(5)
check("the lower counterexample is itself an exact constrained stationary point",
      all(sp.simplify(left-counter_multiplier*right) == 0
          for left, right in zip(counter_gradient_s, counter_gradient_q)),
      f"multiplier={counter_multiplier}")

# In zero-sum six-vertex barycentric coordinates this point is
# (-1,1,1,-1,1,-1).  Its A5 orbit has size ten and does not contain its
# negative; the latter belongs to the opposite ten-point orbit.
counter_barycentric = (-1, 1, 1, -1, 1, -1)
counter_orbit = set()
for action in fibration_actions:
    image = [None]*6
    for source, target in enumerate(action):
        image[target] = counter_barycentric[source]
    counter_orbit.add(tuple(image))
check("the lower stationary witness generates a distinct ten-point A5 orbit",
      len(counter_orbit) == 10
      and tuple(-value for value in counter_barycentric) not in counter_orbit
      and len(fibration_actions)//len(counter_orbit) == 6,
      "orbit size=10; stabilizer order=6; negative lies in the other orbit")

payload = {
    "protocol_commit": "48bbb00",
    "geometry": {
        "A5_actions": len(fibration_actions),
        "six_point_cycle_types": {
            str(key): value for key, value in sorted(actual_cycles.items())
        },
        "W_dimension": 5,
        "Sym3_W_invariant_dimension": sym3_invariant_dimension,
    },
    "cubic_span": {
        "invariant_dimension": 2,
        "single_trace_rank": 1,
        "identity": "4*Tr(A*X^3)=-8*Tr(X^3)",
        "C_box_in_single_trace_span": False,
    },
    "full_fourth_moment": {
        "definition": "S4(X)=Tr((A+X)^4), Tr(X^2)=7200",
        "S4_positive_Box": 933120,
        "S4_negative_Box": 1163520,
        "desired_vertices_stationary": True,
        "counterexample_coefficients": ["-sqrt(5)/6", "sqrt(5)/6",
                                         "sqrt(5)/6", "-sqrt(5)/6",
                                         "sqrt(5)/6"],
        "counterexample_value": str(counterexample_value),
        "counterexample_stationary": True,
        "counterexample_multiplier": str(counter_multiplier),
        "counterexample_A5_orbit_size": len(counter_orbit),
        "counterexample_stabilizer_order": (
            len(fibration_actions)//len(counter_orbit)
        ),
    },
    "verdict": (
        "DERIVED NEGATIVE: the canonical adjacency baseline reproduces only "
        "the old single-trace cubic line, not the equal-weight selector; "
        "the complete fourth moment has an exact stationary point below all "
        "six desired Hopf vertices"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact adjacency-baseline audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: dim Sym^3(W*)^A5 = 2.")
print("DERIVED NEGATIVE: 4 Tr(A X^3) = -8 Tr(X^3), not C_box.")
print("DERIVED NEGATIVE: the complete S4 has a lower exact stationary point.")
print("KILL: the canonical adjacency-baseline fourth moment does not select")
print("      the six Hopf fibrations.")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Exact audit of the six Hopf projectors as a 5D simplex order parameter.

The protocol and falsifiers were frozen in commit 402de35.  The central idea
was already recognized before preregistration and is not presented as a blind
discovery.  All load-bearing identities use exact SymPy arithmetic.
"""

from itertools import combinations, product
import json
from math import comb
from pathlib import Path

import sympy as sp


OUTPUT = Path(__file__).with_name("hopf_projector_cubic.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def exact_zero(expression):
    return sp.simplify(sp.expand(expression)) == 0


def exact_zero_matrix(matrix):
    return all(exact_zero(value) for value in matrix)


def dot(left, right):
    return sum(a*b for a, b in zip(left, right))


def projector(vector):
    column = sp.Matrix(vector)
    return sp.simplify(column*column.T/dot(vector, vector))


def matrix_key(matrix):
    return tuple(sp.radsimp(sp.simplify(value)) for value in matrix)


def frobenius(left, right):
    return sp.simplify(sp.trace(left.T*right))


print("="*78)
print("EXACT FIVE-DIMENSIONAL HOPF PROJECTOR CUBIC")
print("="*78)

sqrt5 = sp.sqrt(5)
phi = (1+sqrt5)/2

vertices = []
for first, second in product((1, -1), repeat=2):
    vertices.extend((
        (sp.Integer(0), sp.Integer(first), sp.Integer(second)*phi),
        (sp.Integer(first), sp.Integer(second)*phi, sp.Integer(0)),
        (sp.Integer(first)*phi, sp.Integer(0), sp.Integer(second)),
    ))
vertices = sorted(set(vertices), key=lambda item: tuple(map(str, item)))

projector_by_key = {}
for vertex in vertices:
    candidate = projector(vertex)
    projector_by_key.setdefault(matrix_key(candidate), candidate)
projectors = tuple(projector_by_key.values())
check("the twelve vertices give exactly six unoriented rank-one projectors",
      len(vertices) == 12 and len(projectors) == 6
      and all(exact_zero_matrix(P*P-P)
              and exact_zero(sp.trace(P)-1) for P in projectors))

identity = sp.eye(3)
centered = tuple(sp.simplify(P-identity/3) for P in projectors)
center_sum = sp.zeros(3)
for tensor in centered:
    center_sum += tensor
gram = sp.Matrix([[frobenius(left, right) for right in centered]
                  for left in centered])
expected_gram = sp.Matrix(6, 6, lambda row, col:
                          sp.Rational(2, 3) if row == col
                          else -sp.Rational(2, 15))
check("the centered projectors have the exact regular-simplex Gram matrix",
      center_sum == sp.zeros(3) and gram == expected_gram,
      "norm^2=2/3; cross/norm=-1/5; sum T_i=0")

flat_columns = [sp.Matrix(tensor).reshape(9, 1) for tensor in centered]
span_rank = sp.Matrix.hstack(*flat_columns).rank()
check("the six centered projectors span all of Sym^2_0(R^3)",
      span_rank == 5,
      f"exact span rank={span_rank}")

sqrt2 = sp.sqrt(2)
sqrt6 = sp.sqrt(6)
orthonormal_basis = (
    sp.diag(1, -1, 0)/sqrt2,
    sp.diag(1, 1, -2)/sqrt6,
    sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])/sqrt2,
    sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])/sqrt2,
    sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])/sqrt2,
)
basis_gram = sp.Matrix([[frobenius(left, right)
                         for right in orthonormal_basis]
                        for left in orthonormal_basis])
coordinates = tuple(sp.Matrix([frobenius(tensor, basis)
                               for basis in orthonormal_basis])
                    for tensor in centered)
frame = sp.zeros(5)
for coordinate in coordinates:
    frame += coordinate*coordinate.T
check("the simplex frame operator is exactly (4/5) I_5",
      exact_zero_matrix(basis_gram-sp.eye(5))
      and exact_zero_matrix(frame-sp.Rational(4, 5)*sp.eye(5)))

# The analysis coordinates s_i=Tr(Q T_i) identify Sym^2_0 with the hyperplane
# sum_i s_i=0.  On ||Q||_F^2=||T_i||_F^2=2/3, the tight-frame identity fixes
# sum_i s_i^2=8/15.  Lagrange stationarity of sum_i s_i^3 then forces every
# s_i to be one of at most two roots.  Enumerating their multiplicity k is an
# exhaustive constrained-coordinate proof, not a numerical search.
radius_squared = sp.Rational(8, 15)
stationary = []
for multiplicity in range(1, 6):
    other_multiplicity = 6-multiplicity
    first = sp.sqrt(radius_squared*other_multiplicity
                    /(6*multiplicity))
    second = sp.simplify(-multiplicity*first/other_multiplicity)
    value = sp.simplify(multiplicity*first**3
                        + other_multiplicity*second**3)
    stationary.append({
        "multiplicity": multiplicity,
        "count": comb(6, multiplicity),
        "first": first,
        "second": second,
        "C3": value,
    })

stationary_values = [record["C3"] for record in stationary]
expected_values = [
    sp.Rational(64, 225),
    8*sp.sqrt(10)/225,
    sp.Integer(0),
    -8*sp.sqrt(10)/225,
    -sp.Rational(64, 225),
]
check("the constrained cubic has the exact five stationary value classes",
      stationary_values == expected_values
      and sum(record["count"] for record in stationary) == 62,
      "values="+str(expected_values)+"; stationary points=62")

maximum_coordinates = [sp.Rational(2, 3)] + [-sp.Rational(2, 15)]*5
minimum_coordinates = [-sp.Rational(2, 3)] + [sp.Rational(2, 15)]*5
check("only the six simplex vertices maximize and their negatives minimize",
      stationary[0]["count"] == 6 and stationary[-1]["count"] == 6
      and stationary[0]["first"] == maximum_coordinates[0]
      and stationary[0]["second"] == maximum_coordinates[1]
      and stationary[-1]["second"] == minimum_coordinates[0]
      and stationary[-1]["first"] == minimum_coordinates[1]
      and all(expected_values[0] > value
              for value in expected_values[1:])
      and all(expected_values[-1] < value
              for value in expected_values[:-1]),
      "C3 max=64/225; min=-64/225")

simplex_cubic_values = [
    sp.simplify(sum(frobenius(tensor, other)**3 for other in centered))
    for tensor in centered
]
negative_cubic_values = [
    sp.simplify(sum(frobenius(-tensor, other)**3 for other in centered))
    for tensor in centered
]
check("direct tensor evaluation realizes the extrema on +/-T_i",
      set(simplex_cubic_values) == {sp.Rational(64, 225)}
      and set(negative_cubic_values) == {-sp.Rational(64, 225)})

x, y, z = sp.symbols("x y z", real=True)
n = sp.Matrix([x, y, z])
r2 = sp.expand(n.dot(n))
Q = sp.expand(n*n.T-r2*identity/3)
linear_coordinates = [frobenius(Q, tensor) for tensor in centered]
C3 = sp.expand(sum(value**3 for value in linear_coordinates))
axis_quadratics = [sp.expand((n.T*P*n)[0]) for P in projectors]
S6 = sp.expand(sum(value**3 for value in axis_quadratics))
pullback_residual = sp.expand(C3-S6+sp.Rational(34, 45)*r2**3)
check("the exact pullback is C3(Q(n))=S6(n)-34/45*r^6",
      exact_zero(pullback_residual))

payload = {
    "protocol_commit": "402de35",
    "provenance": (
        "centered-projector and expected simplex idea declared before the "
        "exhaustive verifier; not a blind target discovery"
    ),
    "space": {
        "name": "Sym^2_0(R^3)",
        "real_dimension": 5,
        "centered_projector_span_rank": int(span_rank),
    },
    "simplex": {
        "vertices": len(centered),
        "norm_squared": "2/3",
        "cross_inner_product": "-2/15",
        "normalized_cross_inner_product": "-1/5",
        "frame_operator": "4/5 * I_5",
    },
    "cubic": {
        "definition": "C3(Q)=sum_i Tr(Q*T_i)^3",
        "sphere_Q_norm_squared": "2/3",
        "stationary_points": sum(record["count"] for record in stationary),
        "stationary_classes": [
            {
                "multiplicity": record["multiplicity"],
                "count": record["count"],
                "C3": str(record["C3"]),
            }
            for record in stationary
        ],
        "global_maximum": "64/225 on the six T_i",
        "global_minimum": "-64/225 on the six -T_i",
    },
    "pullback": "C3(n*n^T-|n|^2*I/3)=S6(n)-34/45*|n|^6",
    "interpretation_boundary": (
        "the canonical 5D cubic is derived; its occurrence and sign in a "
        "certified fluctuated action remain open"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact projector-cubic audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: the six centered Hopf projectors form a regular 5-simplex.")
print("DERIVED: their equal-weight cubic selects exactly +/- the simplex vertices.")
print("DERIVED: the old degree-six vector selector is this cubic pulled back by Q(n).")
print("OPEN: a certified linear Q fluctuation, the cubic coefficient, and its sign.")
raise SystemExit(0 if passed == tests else 1)

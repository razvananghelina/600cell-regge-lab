#!/usr/bin/env python3
"""Exact falsification of selection by the canonical Hopf--Box cubic.

The protocol and decision boundary were frozen in commit 2b754de before the
stationarity and level-set calculations.  Geometry is reused from the
registered six-fibration constructor.  All load-bearing calculations after
that construction use exact integer/rational polynomial arithmetic.
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


OUTPUT = Path(__file__).with_name("hopf_box_cubic_selection.json")
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


print("="*78)
print("EXACT HOPF BOX CUBIC SELECTION TEST")
print("="*78)

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
fiber_adjacencies = [
    np.rint(build_fiber_adjacency(adjacency, fibration)).astype(np.int64)
    for fibration in fibrations
]
boxes = [6*fiber-adjacency for fiber in fiber_adjacencies]
check("the fixed constructor supplies the six certified Hopf Box operators",
      len(vertices) == 120 and len(fibrations) == 6
      and all(matrix.shape == (120, 120) for matrix in boxes))

# Use the exact basis E_a=Box_a-Box_5 of the five-dimensional span W.
basis = [boxes[index]-boxes[5] for index in range(5)]
u = sp.symbols("u0:5", real=True)
gram = sp.Matrix([
    [trace_product(basis[row], basis[col]) for col in range(5)]
    for row in range(5)
])
expected_gram = sp.Matrix(5, 5,
                          lambda row, col: 17280 if row == col else 8640)
check("E_a=Box_a-Box_5 is an exact positive-definite basis of W",
      gram == expected_gram and gram.det() > 0,
      f"det(G)={gram.det()}")

triple = [[[
    trace_product(basis[first], basis[second], basis[third])
    for third in range(5)] for second in range(5)] for first in range(5)]
q = sp.expand(sum(gram[first, second]*u[first]*u[second]
                  for first in range(5) for second in range(5)))
f = sp.expand(sum(triple[first][second][third]
                  *u[first]*u[second]*u[third]
                  for first in range(5)
                  for second in range(5)
                  for third in range(5)))

q_reduced = sp.expand(q/17280)
expected_q_reduced = sp.expand(
    sum(value**2 for value in u)
    + sum(u[left]*u[right]
          for left in range(5) for right in range(left+1, 5))
)
check("q=Tr(X^2) is reconstructed coefficientwise in exact arithmetic",
      q_reduced == expected_q_reduced)

expected_f_reduced = -(
    u[0]**2*u[1] + u[0]**2*u[2] + u[0]*u[1]**2
    + 2*u[0]*u[1]*u[2] + 2*u[0]*u[1]*u[3]
    + u[0]*u[2]**2 + 2*u[0]*u[2]*u[4]
    + u[1]**2*u[3] + u[1]*u[3]**2 + 2*u[1]*u[3]*u[4]
    + u[2]**2*u[4] + 2*u[2]*u[3]*u[4] + u[2]*u[4]**2
    + u[3]**2*u[4] + u[3]*u[4]**2
)
check("f=Tr(X^3) is reconstructed coefficientwise in exact arithmetic",
      sp.expand(f/155520) == sp.expand(expected_f_reduced),
      f"nonzero monomials={len(sp.Poly(f, *u).terms())}")

# In this basis Box_i has coordinate 5/6 at i and -1/6 elsewhere;
# Box_5 has all five coordinates -1/6.
vertex_coordinates = []
for index in range(5):
    coordinate = [sp.Rational(-1, 6)]*5
    coordinate[index] = sp.Rational(5, 6)
    vertex_coordinates.append(sp.Matrix(coordinate))
vertex_coordinates.append(sp.Matrix([sp.Rational(-1, 6)]*5))

reconstructed = [sum((coordinate[index]*basis[index]
                      for index in range(5)),
                     np.zeros_like(adjacency, dtype=object))
                 for coordinate in vertex_coordinates]
check("the six rational coordinate vectors reconstruct exactly the six Box_i",
      all(np.array_equal(np.asarray(matrix, dtype=object), box)
              for matrix, box in zip(reconstructed, boxes)))

q_values = [sp.expand(q).subs(dict(zip(u, coordinate)))
            for coordinate in vertex_coordinates]
f_values = [sp.expand(f).subs(dict(zip(u, coordinate)))
            for coordinate in vertex_coordinates]
check("all six Box vertices lie on q=7200 and f=14400",
      set(q_values) == {sp.Integer(7200)}
      and set(f_values) == {sp.Integer(14400)})

gradient_q = sp.Matrix([sp.diff(q, value) for value in u])
gradient_f = sp.Matrix([sp.diff(f, value) for value in u])
stationarity_residuals = []
for coordinate in vertex_coordinates:
    substitution = dict(zip(u, coordinate))
    stationarity_residuals.append(
        (gradient_f-3*gradient_q).subs(substitution)
    )
check("all six vertices are exact stationary points with multiplier lambda=3",
      all(residual == sp.zeros(5, 1)
          for residual in stationarity_residuals))

# An independent projective line supplies an exact extra solution.  Put
# u(t)=(-1,2t-1,0,0,0).  A positive rescaling reaches q=7200.  The condition
# that the rescaled cubic equal +14400 is exactly 1800*f(t)^2=q(t)^3,
# together with f(t)>0.
t = sp.symbols("t", real=True)
line_substitution = {
    u[0]: -1,
    u[1]: 2*t-1,
    u[2]: 0,
    u[3]: 0,
    u[4]: 0,
}
q_line = sp.factor(q.subs(line_substitution))
f_line = sp.factor(f.subs(line_substitution))
P = sp.Poly(
    256*t**6-1152*t**5+1764*t**4-972*t**3
    -27*t**2+162*t-27,
    t,
    domain=sp.ZZ,
)
level_residual = sp.factor(1800*f_line**2-q_line**3)
check("the level equation on the rational line reduces to the declared P(t)",
      sp.expand(q_line-17280*(4*t**2-6*t+3)) == 0
      and sp.expand(f_line-311040*(t-1)*(2*t-1)) == 0
      and sp.expand(level_residual
                    + 1289945088000*P.as_expr()) == 0)

left = sp.Rational(1, 10)
right = sp.Rational(1, 4)
check("P changes sign and has exactly one real root alpha in (1/10,1/4)",
      P.eval(left) < 0 < P.eval(right)
      and P.count_roots(left, right) == 1,
      f"P(1/10)={P.eval(left)}; P(1/4)={P.eval(right)}")
check("the isolated root is simple and the cubic has the positive target sign",
      sp.gcd(P, P.diff()).degree() == 0
      and f_line.subs(t, left) > 0
      and f_line.subs(t, right) > 0
      and q_line.subs(t, left) > 0
      and q_line.subs(t, right) > 0,
      "on the whole interval t<1/4<1/2<1, so f(t)>0 and q(t)>0")

# If alpha is that algebraic root and s=sqrt(7200/q(alpha)), then P(alpha)=0
# implies q(sX)=7200 and f(sX)=+14400.  The following identity certifies the
# normalization without a floating-point RootOf evaluation.
check("positive normalization of the algebraic root gives f=+14400 exactly",
      sp.Integer(7200)**3/sp.Integer(1800) == sp.Integer(14400)**2)

# No Box vertex has its last three E-coordinates zero: their entries are
# +/-1/6 or 5/6.  Thus the algebraic point is genuinely additional.
check("the algebraic witness cannot be any of the six Box vertices",
      all(any(coordinate[index] != 0 for index in (2, 3, 4))
          for coordinate in vertex_coordinates))

# P is square-free, hence the line meets the homogeneous level hypersurface
# transversely at alpha.  Were df proportional to dq there, Euler homogeneity
# would force d(1800 f^2-q^3)=0 in every direction.  The transverse line
# derivative contradicts that.  Therefore (q,f) has rank two at the scaled
# witness and its common level set is locally a real 3-manifold in dim W=5.
check("the extra solution is regular and lies on a local 3D continuum",
      sp.gcd(P, P.diff()) == sp.Poly(1, t, domain=sp.ZZ)
      and gram.rank() == 5,
      "simple line root => d(level residual)!=0 => dq,df independent")

payload = {
    "protocol_commit": "2b754de",
    "arithmetic": (
        "integer incidence matrices followed by exact integer/rational "
        "polynomial arithmetic and Sturm root counting"
    ),
    "space_dimension": 5,
    "sphere": "Tr(X^2)=7200",
    "cubic_level": "Tr(X^3)=14400",
    "desired_vertices": {
        "count": 6,
        "stationary": True,
        "lagrange_multiplier": 3,
    },
    "extra_witness": {
        "line": "u=(-1,2*t-1,0,0,0)",
        "root_polynomial": str(P.as_expr()),
        "isolating_interval": ["1/10", "1/4"],
        "roots_in_interval": int(P.count_roots(left, right)),
        "root_simple": bool(sp.gcd(P, P.diff()).degree() == 0),
        "normalization": "s=sqrt(7200/q(alpha))",
        "local_level_set_dimension": 3,
    },
    "verdict": (
        "DERIVED NEGATIVE: the existing q and cubic constraints do not "
        "select the six Hopf Box vertices; an additional regular algebraic "
        "point and therefore a local real 3D continuum share their levels"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact cubic-selection audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: all six Hopf Box vertices are stationary with lambda=3.")
print("DERIVED NEGATIVE: their q and cubic levels also contain a regular")
print("                  algebraic point and a local real 3D continuum.")
print("KILL: Tr(X^2)=7200 and Tr(X^3)=14400 do not select six fibrations.")
raise SystemExit(0 if passed == tests else 1)

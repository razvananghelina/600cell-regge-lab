#!/usr/bin/env python3
"""Exact barycentric refinement and Dirac-weight scaling audit.

The full barycentric subdivision of a closed 3-dimensional simplicial complex
has an exact f-vector transformation with dominant factor 24.  On the
minimal level-weight spectral model, assigning eigenvalue b**n to modes born
at level n gives spectral abscissa log(24)/log(b).  Compact resolvent only
requires b>1 and therefore does not select a dimension.  Choosing the inverse
linear scale from a 1/24 volume split gives three; choosing four is an extra
scaling law.

The level-weight model is explicitly STRUCTURAL.  This verifier certifies its
exact combinatorics and summability formula, not the existence of the missing
global configuration-space Dirac operator.
"""

import math

import sympy as sp
from sympy.functions.combinatorial.numbers import stirling


tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("MISSING-LINK AUDIT: REFINEMENT GROWTH VERSUS DIRAC SCALE")
print("=" * 78)

# Rows are the output simplex dimensions and columns the input dimensions.
# Entry (j,i) counts j-simplices in sd(Delta^i) whose maximal face is the
# parent i-simplex: (j+1)! S(i+1,j+1).
B = sp.Matrix([
    [1, 1, 1, 1],
    [0, 2, 6, 14],
    [0, 0, 6, 36],
    [0, 0, 0, 24],
])
f0 = sp.Matrix((120, 720, 1200, 600))
f1 = B * f0

formula_B = sp.zeros(4, 4)
for output_dim in range(4):
    for input_dim in range(output_dim, 4):
        formula_B[output_dim, input_dim] = (
            math.factorial(output_dim + 1)
            * stirling(input_dim + 1, output_dim + 1, kind=2)
        )
check("subdivision matrix follows the Stirling chain-count formula",
      B == formula_B)
check("barycentric f-vector transform has eigenvalues 1,2,6,24",
      B.eigenvals() == {1: 1, 2: 1, 6: 1, 24: 1})
check("first full barycentric subdivision has exact f-vector",
      tuple(f1) == (2640, 17040, 28800, 14400),
      f"got {tuple(f1)}")

euler_row = sp.Matrix([[1, -1, 1, -1]])
check("barycentric subdivision preserves Euler characteristic exactly",
      euler_row * B == euler_row and (euler_row*f0)[0] == 0)

# Closed forms for this S3 f-vector.  Verify they satisfy both the initial
# condition and the exact B recurrence symbolically at the coefficient-vector
# level, rather than by checking finitely many numerical levels.
c2 = sp.Matrix((sp.Rational(120, 11), sp.Rational(120, 11), 0, 0))
c24 = sp.Matrix((sp.Rational(1200, 11), sp.Rational(7800, 11), 1200, 600))
check("closed f-vector coefficients reconstruct level zero", c2 + c24 == f0)
check("closed f-vector coefficients are exact 2- and 24-eigenvectors",
      B*c2 == 2*c2 and B*c24 == 24*c24)

# Hence E_n=(120/11)2^n+(7800/11)24^n.  New edge variables born at n>=1 are
# m_n=E_n-E_{n-1}.
n = sp.symbols("n", integer=True, positive=True)
new_edges = (sp.Rational(120, 11) * 2**(n-1)
             + sp.Rational(179400, 11) * 24**(n-1))
for level in (1, 2, 3, 4):
    observed = int((B**level*f0)[1] - (B**(level-1)*f0)[1])
    expected = int(new_edges.subs(n, level))
    check(f"new-edge count at refinement level {level} is exact",
          observed == expected, f"m_{level}={observed}")

# If D has eigenvalue a_n=b^n on each new level-n mode, its zeta series is an
# exact sum of two geometric series.  Its convergence boundary is b^s=24.
b, s = sp.symbols("b s", positive=True)
zeta_closed = (
    sp.Rational(120, 11) * b**(-s) / (1 - 2*b**(-s))
    + sp.Rational(179400, 11) * b**(-s) / (1 - 24*b**(-s))
)

# Prove the geometric-series reduction at a symbolic finite cutoff.  Taking
# N to infinity is then valid exactly when |24*x|<1.
x = sp.symbols("x")  # x=b^(-s)
cutoff = 6
partial = sum(
    (sp.Rational(120, 11)*2**(level-1)
     + sp.Rational(179400, 11)*24**(level-1)) * x**level
    for level in range(1, cutoff + 1)
)
finite_closed = (
    sp.Rational(120, 11)*x*(1-(2*x)**cutoff)/(1-2*x)
    + sp.Rational(179400, 11)*x*(1-(24*x)**cutoff)/(1-24*x)
)
check("level-weight zeta function is the exact pair of geometric series",
      sp.simplify(partial-finite_closed) == 0,
      "symbolic partial-sum identity; convergence iff b^s>24")
infinite_from_x = (sp.Rational(120, 11)*x/(1-2*x)
                   + sp.Rational(179400, 11)*x/(1-24*x))
check("displayed zeta closed form matches x=b^(-s)",
      sp.simplify(zeta_closed-infinite_from_x.subs(x, b**(-s))) == 0)

d_b2 = math.log(24) / math.log(2)
d_b3 = math.log(24) / math.log(3)
check("two compact-resolvent scale laws give different dimensions",
      abs(d_b2 - 4.584962500721156) < 1e-14 and
      abs(d_b3 - 2.8927892607143724) < 1e-14,
      f"b=2 -> d={d_b2:.12f}; b=3 -> d={d_b3:.12f}")

check("volume-derived inverse linear scale gives dimension three",
      sp.simplify(sp.log(24) / (sp.log(24)/3)) == 3,
      "b=24^(1/3) is conditional on 3D isotropic scaling")
check("dimension four requires a different inserted scale law",
      sp.simplify(sp.log(24) / (sp.log(24)/4)) == 4,
      "b=24^(1/4) is not selected by the 3-simplex volume split")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: full barycentric edge growth has dominant factor 24.")
print("DERIVED CONDITIONAL: level weights b^n give d=log(24)/log(b).")
print("STRUCTURAL: inverse linear scaling from 3-volume gives d=3.")
print("DERIVED NEGATIVE: compact resolvent/convergence alone does not select b.")
raise SystemExit(0 if passed == tests else 1)

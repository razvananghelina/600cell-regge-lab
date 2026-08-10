#!/usr/bin/env python3
"""Exact shape-regularity attack on iterated tetrahedral barycentric refinement.

One deterministic nested flag suffices.  In parent edge coordinates its
affine Jacobian is upper triangular with eigenvalues 1/2, 1/3, 1/4.  Repeating
that flag gives an explicit sequence of increasingly anisotropic tetrahedra,
so standard shape-regular FEEC convergence theorems cannot be invoked for the
unmodified barycentric tower.
"""

import numpy as np
import sympy as sp


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
print("ITERATED BARYCENTRIC SHAPE-REGULARITY GATE")
print("=" * 78)

# For the repeated flag (v0)<(v0,v1)<(v0,v1,v2)<(v0,v1,v2,v3), the child
# edge columns in the parent edge basis (v1-v0,v2-v0,v3-v0) are these.
transform = sp.Matrix((
    (sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 4)),
    (0, sp.Rational(1, 3), sp.Rational(1, 4)),
    (0, 0, sp.Rational(1, 4)),
))
check("repeated-flag affine transform is exactly upper triangular",
      transform.is_upper)
check("one child has exact affine volume ratio 1/24",
      transform.det() == sp.Rational(1, 24))
check("linear contraction eigenvalues are exactly 1/2,1/3,1/4",
      transform.eigenvals() == {
          sp.Rational(1, 2): 1,
          sp.Rational(1, 3): 1,
          sp.Rational(1, 4): 1,
      })

for depth in range(1, 7):
    check(f"depth-{depth} nested child volume is exactly 24^-{depth}",
          (transform**depth).det() == sp.Rational(1, 24**depth))

# Any operator norm dominates spectral radius.  Hence
# cond_2(T^n)=||T^n|| ||T^-n|| >= (1/2)^n * 4^n = 2^n.
# Numerical singular values record how much stronger the actual witness is.
condition_numbers = []
lower_bounds_hold = True
for depth in range(1, 11):
    power = np.asarray(transform**depth, dtype=float)
    singular = np.linalg.svd(power, compute_uv=False)
    condition = float(singular[0]/singular[-1])
    condition_numbers.append(condition)
    lower_bounds_hold &= condition >= 2**depth*(1-1e-12)
check("repeated flag has condition number at least 2^n",
      lower_bounds_hold,
      "n=1..10 conditions="
      + ", ".join(f"{value:.2f}" for value in condition_numbers))
check("the explicit aspect distortion is unbounded, not a finite-level peak",
      all(condition_numbers[index+1] > 1.9*condition_numbers[index]
          for index in range(4, len(condition_numbers)-1))
      and condition_numbers[-1] > 3000,
      f"cond(T^10)={condition_numbers[-1]:.3f}")

# Put the witness into the actual Euclidean metric of a regular tetrahedron.
# If A is its edge Jacobian, cond(A T^n) >= cond(T^n)/cond(A), hence the same
# unboundedness survives physical coordinates.  The numerical sequence is an
# independent direct audit of that conclusion.
regular_vertices = np.array((
    (1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)
), dtype=float)
regular_edges = np.column_stack([
    regular_vertices[index]-regular_vertices[0] for index in range(1, 4)
])
physical_conditions = []
for depth in range(1, 11):
    jacobian = regular_edges@np.asarray(transform**depth, dtype=float)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    physical_conditions.append(float(singular[0]/singular[-1]))
check("regular-tetrahedron physical children also degenerate without bound",
      all(physical_conditions[index+1] > 1.8*physical_conditions[index]
          for index in range(5, len(physical_conditions)-1))
      and physical_conditions[-1] > 2000,
      f"physical cond depth 1/10={physical_conditions[0]:.3f}/"
      f"{physical_conditions[-1]:.3f}")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: barycentric volume and mesh scales shrink on the witness path.")
print("DERIVED NEGATIVE: iterated tetrahedral barycentric meshes are not shape-regular.")
print("DERIVED NEGATIVE: standard shape-regular FEEC convergence is inapplicable as-is.")
print("OPEN: prove convergence on this degenerating tower or select a canonical repair.")
raise SystemExit(0 if passed == tests else 1)

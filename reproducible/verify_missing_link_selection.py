#!/usr/bin/env python3
"""Exact local certificate for the projective spectral-selection audit.

This verifier does not claim that a spectral action is derived.  It checks the
algebraic statement used in ``missing_link_audit.md``: on the regular
tetrahedron and its complete barycentric subdivision, every cylindrically
compatible cometric in the already-established family splits into a fixed
horizontal block and a rank-44 vertical block ``t I``.  Consequently a
spectral functional ``Tr f(D^2)`` depends on the free scale only through
``44 f(t)``.  The functional selects a finite positive ``t`` only after the
function ``f`` (and its scale) has supplied one.

All matrix identities are exact over Q.  The exponential derivative is a
symbolic sign argument, not a floating-point fit.
"""

from itertools import combinations

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
print("MISSING-LINK AUDIT: CAN THE SPECTRAL ACTION SELECT THE PROJECTIVE SCALE?")
print("=" * 78)

# Complete barycentric subdivision of one tetrahedron.  A nonempty subset of
# the four parent vertices is a fine vertex; comparable subsets give the 50
# barycentric edges used by the HD-600 connection-space calculation.
coarse_edges = list(combinations(range(4), 2))
fine_edges = [
    (small, large)
    for small in range(1, 16)
    for large in range(small + 1, 16)
    if small & large == small
]
fine_edge_index = {edge: index for index, edge in enumerate(fine_edges)}

check("regular tetrahedron has 6 coarse and 50 barycentric edges",
      (len(coarse_edges), len(fine_edges)) == (6, 50))

# Coarse edge holonomy is the composition of the two oriented half-edges.
A = sp.zeros(6, 50)
for row, (i, j) in enumerate(coarse_edges):
    midpoint = (1 << i) | (1 << j)
    A[row, fine_edge_index[(1 << i, midpoint)]] = 1
    A[row, fine_edge_index[(1 << j, midpoint)]] = -1

H = A.T * (A * A.T).inv()
Q = sp.eye(50) - H * A

check("coarse composition has full row rank", A.rank() == 6)
check("H is an exact right inverse", A * H == sp.eye(6))
check("Q is the exact vertical orthogonal projector",
      Q.T == Q and Q * Q == Q and A * Q == sp.zeros(6, 50))
check("the vertical refinement sector has rank 44", Q.rank() == 44)
check("horizontal/vertical orthogonality is independent of the coarse metric",
      Q * H == sp.zeros(50, 6) and H.T * Q == sp.zeros(6, 50))

# The result is independent of the positive coarse cometric.  Identity is a
# convenient exact representative for checking the block decomposition.
Kc = sp.eye(6)
Kh = H * Kc * H.T
check("horizontal and vertical blocks are exactly orthogonal",
      Kh * Q == sp.zeros(50, 50) and Q * Kh == sp.zeros(50, 50))
check("every t has the same coarse pullback",
      A * (Kh + Q) * A.T == Kc and
      A * (Kh + 2 * Q) * A.T == Kc)

# Orthogonality and Q^n=Q imply the general polynomial trace theorem.  Check
# the first three moments explicitly at two rational values; the proof for a
# polynomial is then linearity in its monomials.
for power in (1, 2, 3):
    horizontal_moment = sp.trace(Kh ** power)
    for t in (sp.Rational(1, 3), sp.Rational(5, 2)):
        observed = sp.trace((Kh + t * Q) ** power)
        expected = horizontal_moment + 44 * t ** power
        check(f"moment p={power}, t={t}: vertical contribution is 44*t^p",
              observed == expected)

# Three standard variational choices illustrate the exhaustive alternatives.
# Positive moments are monotone and run to the degenerate boundary t=0.
a, b, t, s = sp.symbols("a b t s", positive=True)
positive_moment_derivative = sp.diff(44 * (a * t + b * t**2), t)
check("positive quadratic/quartic moments have no positive critical t",
      positive_moment_derivative == 44 * (a + 2 * b * t))

# A symmetry-breaking polynomial has a finite critical point, but its value is
# precisely the freely chosen coefficient ratio a/(2b).
broken_derivative = sp.diff(44 * (b * t**2 - a * t), t)
critical_t = sp.solve(sp.Eq(broken_derivative, 0), t)
check("symmetry-breaking polynomial fixes only the input coefficient ratio",
      critical_t == [a / (2 * b)])
check("two admissible coefficient choices select two different scales",
      critical_t[0].subs({a: 2, b: 1}) == 1 and
      critical_t[0].subs({a: 4, b: 1}) == 2,
      "(a,b)=(2,1) gives t=1; (4,1) gives t=2")

# A positive heat cutoff is strictly decreasing in t and therefore has no
# finite stationary point.  SymPy supplies the exact derivative; positivity
# of r,s,exp(-st) fixes its sign.
heat_derivative = sp.diff(44 * sp.exp(-s * t), t)
check("heat spectral action is strictly decreasing for finite positive t",
      heat_derivative == -44 * s * sp.exp(-s * t))

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: Tr f(K_h+tQ) = Tr_horizontal f(K_h) + 44 f(t).")
print("DERIVED NEGATIVE: positive moments and heat cutoffs select no finite t.")
print("STRUCTURAL: a nonmonotone f can select t only through its input coefficients.")
print("OPEN: a target-independent action and scale derived before choosing t.")
raise SystemExit(0 if passed == tests else 1)

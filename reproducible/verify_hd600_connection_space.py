#!/usr/bin/env python3
"""Projective connection-space gate for the HD-600 route.

No particle algebra or phenomenological target occurs here.  The calculation
asks what edge subdivision and gauge symmetry determine before a Dirac
operator on the connection space is chosen.
"""

from collections import Counter

import sympy as sp


tests = passed = 0


def check(name, ok, detail=""):
    global tests, passed
    tests += 1
    ok = bool(ok)
    passed += int(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("HD-600 GATE 2: PROJECTIVE CONNECTION SPACE AND METRIC FREEDOM")
print("=" * 78)

# The 600-cell 1-skeleton is connected.  One group element per unoriented
# edge describes a graph connection; the reverse orientation is its inverse.
vertices = 120
edges = 720
cycle_rank = edges-vertices+1
check("600-cell graph has 601 independent cycles",
      cycle_rank == 601,
      f"E-V+1={cycle_rank}")

# For SU(2), the generic stabilizer of an irreducible graph connection is its
# discrete centre.  Based gauge fixing leaves d(E-V+1) dimensions before the
# final global conjugation.
su2_edge_dim = 3*edges
su2_gauge_dim = 3*vertices
su2_based_dim = 3*cycle_rank
su2_generic_quotient_dim = 3*(edges-vertices)
check("SU(2) graph connection and generic quotient dimensions are exact",
      (su2_edge_dim, su2_gauge_dim, su2_based_dim,
       su2_generic_quotient_dim) == (2160, 360, 1803, 1800),
      "edge=2160, gauge=360, based=1803, full generic quotient=1800")

u2_edge_dim = 4*edges
u2_gauge_dim = 4*vertices
u2_based_dim = 4*cycle_rank
u2_generic_quotient_dim = 4*(edges-vertices)+1
check("U(2) graph connection is a genuinely larger, unselected space",
      (u2_edge_dim, u2_gauge_dim, u2_based_dim,
       u2_generic_quotient_dim) == (2880, 480, 2404, 2401),
      "edge=2880, gauge=480, based=2404, full generic quotient=2401; "
      "constant central U(1) is the generic stabilizer")

# Exact finite-group control for the Haar push-forward theorem.  For every
# group, multiplication GxG->G has exactly |G| preimages of each element.
# Z/120 is used only as a mechanically exhaustive finite control; the compact
# statement follows from left/right invariance and uniqueness of normalized
# Haar measure, not from this example.
order = 120
product_counts = Counter((a+b) % order
                         for a in range(order) for b in range(order))
check("uniform measure is exactly projective under group multiplication",
      len(product_counts) == order and
      set(product_counts.values()) == {order},
      "finite control: every coarse element has exactly 120 fine preimages")

# Linearize one edge subdivision at the identity:
#
#       dm(X,Y) = X+Y.
#
# For each Lie-algebra component, the most general exchange-symmetric,
# Ad-invariant quadratic form is M=[[a,b],[b,a]].  Requiring the canonical
# horizontal lift (Z/2,Z/2) to have the coarse norm c*|Z|^2 imposes only
# a+b=2c.  The vertical eigenvalue a-b remains an arbitrary positive delta.
c, delta = sp.symbols("c delta", positive=True)
a = c + delta/2
b = c - delta/2
metric = sp.Matrix(((a, b), (b, a)))
exchange = sp.Matrix(((0, 1), (1, 0)))
horizontal_lift = sp.Matrix((sp.Rational(1, 2), sp.Rational(1, 2)))
vertical = sp.Matrix((1, -1))
dm = sp.Matrix(((1, 1),))

check("metric family is invariant under exchange of the two subedges",
      sp.simplify(exchange.T*metric*exchange-metric) == sp.zeros(2))
check("every metric in the family makes the horizontal lift isometric",
      sp.simplify((horizontal_lift.T*metric*horizontal_lift)[0]-c) == 0 and
      dm*horizontal_lift == sp.eye(1),
      "a+b=2c is the only horizontal constraint")
vertical_norm = sp.simplify((vertical.T*metric*vertical)[0])
check("the new vertical subdivision mode retains one free positive scale",
      vertical_norm == 2*delta,
      f"||(X,-X)||^2 coefficient={vertical_norm}")

# The inverse metric controls the principal symbol of the Laplacian.  A coarse
# covector pulls back to (p,p); its norm must be p^2/c independently of delta.
coarse_covector_pullback = sp.Matrix((1, 1))
pulled_covector_norm = sp.simplify(
    (coarse_covector_pullback.T*metric.inv()
     * coarse_covector_pullback)[0])
check("coarse Laplacian symbol intertwines for every vertical scale",
      pulled_covector_norm == 1/c,
      f"pulled-back principal symbol={pulled_covector_norm}")

# Give three explicit, inequivalent positive metrics with identical coarse
# pullback.  The product metric is only the middle choice delta=2c; setting
# the cross term b to zero is therefore an additional locality assumption.
examples = [sp.simplify(metric.subs({c: 1, delta: value}))
            for value in (1, 2, 3)]
determinants = [sp.det(item) for item in examples]
check("at least three inequivalent positive metrics pass all projective gates",
      len({tuple(item) for item in examples}) == 3 and
      all(value > 0 for value in determinants),
      f"delta=(1,2,3), determinants={determinants}")
check("the product metric is an extra choice, not a consequence",
      examples[1][0, 1] == 0 and
      examples[0][0, 1] != 0 and examples[2][0, 1] != 0,
      "b=0 occurs only at delta=2c")

# First-order consequence.  On the abelian maximal-torus/local linear model,
# the Hodge--Dirac square is the Laplacian with Fourier eigenvalue
# k^T M^{-1} k.  Coarse cylindrical modes have k=(n,n) and are independent
# of delta.  New subdivision modes have k=(n,-n) and retain delta.  The same
# freedom is present in the principal symbol for every Lie-algebra component
# of SU(2); nonabelian lower-order terms cannot remove it.
coarse_mode = sp.Matrix((1, 1))
new_vertical_mode = sp.Matrix((1, -1))
coarse_eigenvalue = sp.simplify(
    (coarse_mode.T*metric.inv()*coarse_mode)[0])
vertical_eigenvalue = sp.simplify(
    (new_vertical_mode.T*metric.inv()*new_vertical_mode)[0])
check("Hodge--Dirac square agrees on every inherited coarse mode",
      coarse_eigenvalue == 1/c,
      f"lambda_coarse={coarse_eigenvalue}")
check("Hodge--Dirac spectrum of every new vertical mode is unselected",
      vertical_eigenvalue == 2/delta,
      f"lambda_vertical={vertical_eigenvalue}")
vertical_examples = [vertical_eigenvalue.subs(delta, value)
                     for value in (1, 2, 3)]
check("projectively equivalent Dirac operators have different new spectra",
      vertical_examples == [2, 1, sp.Rational(2, 3)],
      f"delta=(1,2,3) gives lambda_vertical={vertical_examples}")

print("\n" + "=" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("=" * 78)
print("HAAR_MEASURE=DERIVED_PROJECTIVELY_CONSISTENT")
print("HORIZONTAL_METRIC_SCALE=DERIVED_FIXED")
print("VERTICAL_METRIC_SCALE=DERIVED_FREE_POSITIVE_PARAMETER")
print("PROJECTIVE_CANONICITY_SELECTS_DIRAC=DERIVED_NEGATIVE")
print("GEOMETRIC_VERTICAL_SCALE=OPEN")
print("SM_TARGET_COMPARISON=NOT_PERFORMED")

if passed != tests:
    raise SystemExit(1)

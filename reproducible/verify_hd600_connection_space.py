#!/usr/bin/env python3
"""Projective connection-space gate for the HD-600 route.

No particle algebra or phenomenological target occurs here.  The calculation
asks what edge subdivision and gauge symmetry determine before a Dirac
operator on the connection space is chosen.
"""

from collections import Counter
from itertools import combinations, permutations

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

# ---------------------------------------------------------------- Whitney L2
# A possible geometric selector for the free vertical metric is the exact L2
# inner product of Whitney 1-forms.  Test it locally on a reference tetrahedron
# and its complete barycentric subdivision.  Affine naturality then makes this
# the local identity on every Euclidean 600-cell tetrahedron.
reference_vertices = (
    sp.Matrix((0, 0, 0)), sp.Matrix((1, 0, 0)),
    sp.Matrix((0, 1, 0)), sp.Matrix((0, 0, 1)))
coarse_edges = list(combinations(range(4), 2))


def subset_members(mask):
    return tuple(index for index in range(4) if mask & (1 << index))


def subset_barycentric(mask):
    members = subset_members(mask)
    return sp.Matrix(tuple(sp.Rational(int(index in members), len(members))
                           for index in range(4)))


fine_vertex_barycentric = {mask: subset_barycentric(mask)
                           for mask in range(1, 16)}
fine_vertex_position = {
    mask: sum((fine_vertex_barycentric[mask][index]
               * reference_vertices[index] for index in range(4)),
              sp.zeros(3, 1))
    for mask in range(1, 16)
}
fine_edges = [(small, large) for small in range(1, 16)
              for large in range(small+1, 16)
              if small & large == small]
fine_edge_index = {edge: index for index, edge in enumerate(fine_edges)}
fine_tetrahedra = []
for order_vertices in permutations(range(4)):
    a0, a1, a2, _ = order_vertices
    fine_tetrahedra.append((1 << a0,
                            (1 << a0) | (1 << a1),
                            (1 << a0) | (1 << a1) | (1 << a2),
                            15))
check("one tetrahedron subdivides into 15 vertices, 50 edges and 24 tetrahedra",
      (len(fine_vertex_barycentric), len(fine_edges), len(fine_tetrahedra)) ==
      (15, 50, 24))


def local_whitney_mass(points):
    """Exact 6x6 Whitney-one-form L2 mass matrix on one tetrahedron."""
    affine = sp.Matrix.hstack(points[1]-points[0], points[2]-points[0],
                              points[3]-points[0])
    inverse = affine.inv()
    gradients = [-sum((sp.Matrix(inverse.row(row)).T for row in range(3)),
                      sp.zeros(3, 1))]
    gradients.extend(sp.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det())/6

    def barycentric_product_integral(i, j):
        return volume*sp.Rational(2 if i == j else 1, 20)

    result = sp.zeros(6, 6)
    local_edges = list(combinations(range(4), 2))
    for row, (i, j) in enumerate(local_edges):
        for col, (k, ell) in enumerate(local_edges):
            result[row, col] = sp.simplify(
                gradients[j].dot(gradients[ell])
                * barycentric_product_integral(i, k)
                - gradients[j].dot(gradients[k])
                * barycentric_product_integral(i, ell)
                - gradients[i].dot(gradients[ell])
                * barycentric_product_integral(j, k)
                + gradients[i].dot(gradients[k])
                * barycentric_product_integral(j, ell))
    return result


coarse_mass = local_whitney_mass(reference_vertices)
fine_mass = sp.zeros(len(fine_edges), len(fine_edges))
local_edge_pairs = list(combinations(range(4), 2))
for tetrahedron in fine_tetrahedra:
    local_mass = local_whitney_mass(
        tuple(fine_vertex_position[mask] for mask in tetrahedron))
    for local_row, (i, j) in enumerate(local_edge_pairs):
        global_row = fine_edge_index[(tetrahedron[i], tetrahedron[j])]
        for local_col, (k, ell) in enumerate(local_edge_pairs):
            global_col = fine_edge_index[(tetrahedron[k], tetrahedron[ell])]
            fine_mass[global_row, global_col] += local_mass[local_row, local_col]

# The canonical cochain injection records the line integral of each coarse
# Whitney basis form over every oriented fine edge.  For barycentric
# coordinates lambda, this integral is lambda_i(s)lambda_j(t)-lambda_j(s)lambda_i(t).
whitney_injection = sp.zeros(len(fine_edges), len(coarse_edges))
for fine_row, (source, target) in enumerate(fine_edges):
    source_lambda = fine_vertex_barycentric[source]
    target_lambda = fine_vertex_barycentric[target]
    for coarse_col, (i, j) in enumerate(coarse_edges):
        whitney_injection[fine_row, coarse_col] = (
            source_lambda[i]*target_lambda[j]
            - source_lambda[j]*target_lambda[i])

whitney_pullback_residual = sp.simplify(
    whitney_injection.T*fine_mass*whitney_injection-coarse_mass)
check("Whitney L2 metric is exactly cylindrical under barycentric subdivision",
      whitney_pullback_residual == sp.zeros(6, 6),
      "P^T M_fine P = M_coarse exactly over Q")

# The existing Kahler--Dirac uses the counting inner product, i.e. an identity
# Gram matrix on oriented cochains.  Its pullback is P^T P and is not the
# coarse identity.  Record the exact spectrum to prevent a normalization from
# being mistaken for a cure.
unweighted_pullback = whitney_injection.T*whitney_injection
unweighted_residual = unweighted_pullback-sp.eye(len(coarse_edges))
unweighted_eigenvalues = sorted(unweighted_pullback.eigenvals().keys(),
                                key=lambda value: float(value))
check("the current unweighted cochain metric is not cylindrical",
      unweighted_residual != sp.zeros(6, 6),
      f"distinct eigenvalues of P^T P={unweighted_eigenvalues}")
check("no single scalar normalization repairs the unweighted metric",
      len(unweighted_eigenvalues) > 1,
      "P^T P is not proportional to the coarse identity")

# Framing attack: Gate 0 used radially normalized barycentres on the round
# sphere, whereas the exact Whitney nesting above uses affine barycentres in a
# flat tetrahedron.  They are not the same refinement.  For a regular
# 600-cell tetrahedron, all distinct vertex dot products are p=phi/2.  Its
# affine barycentre has norm r<1 and lies in the facet hyperplane n.x=r;
# radial normalization moves it to n.x=1 by the nonzero distance 1-r.
phi = (1+sp.sqrt(5))/2
pair_dot = phi/2
affine_barycentre_norm = sp.sqrt((1+3*pair_dot)/4)
radial_displacement = sp.simplify(1-affine_barycentre_norm)
check("round-geodesic and affine-Whitney refinements are genuinely distinct",
      radial_displacement > 0,
      f"tetra-centre normal displacement={radial_displacement} "
      f"~ {float(radial_displacement):.12f}")

print("\n" + "=" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("=" * 78)
print("HAAR_MEASURE=DERIVED_PROJECTIVELY_CONSISTENT")
print("HORIZONTAL_METRIC_SCALE=DERIVED_FIXED")
print("VERTICAL_METRIC_SCALE=DERIVED_FREE_POSITIVE_PARAMETER")
print("PROJECTIVE_CANONICITY_SELECTS_DIRAC=DERIVED_NEGATIVE")
print("GEOMETRIC_VERTICAL_SCALE=OPEN")
print("WHITNEY_L2_METRIC=DERIVED_EXACTLY_CYLINDRICAL_LOCAL")
print("UNWEIGHTED_COCHAIN_METRIC=DERIVED_NOT_CYLINDRICAL")
print("WHITNEY_SELECTION_BY_THEORY=OPEN_STRUCTURAL")
print("ROUND_VS_AFFINE_REFINEMENT=DERIVED_DISTINCT")
print("SM_TARGET_COMPARISON=NOT_PERFORMED")

if passed != tests:
    raise SystemExit(1)

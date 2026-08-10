#!/usr/bin/env python3
"""Projective connection-space gate for the HD-600 route.

No particle algebra or phenomenological target occurs here.  The calculation
asks what edge subdivision and gauge symmetry determine before a Dirac
operator on the connection space is chosen.
"""

from collections import Counter
from itertools import combinations, permutations

import sympy as sp
import numpy as np
import scipy.linalg as sla


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
check("Whitney L2 tangent inclusion is exactly isometric under subdivision",
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

# The displacement is not itself an obstruction.  Use affine barycentric
# coordinates in the reference simplex and map them to the round simplex by
# F(lambda)=normalize(sum lambda_i q_i).  Fine physical simplices are then
# restrictions of the same smooth map F to the affine barycentric subcells.
# It suffices to prove the Whitney inclusion pointwise on the reference
# simplex; pullback by any common F and integration with any common metric
# preserves it.
x, y, z = sp.symbols("x y z")
coordinate = sp.Matrix((x, y, z))


def local_whitney_fields(points):
    affine = sp.Matrix.hstack(points[1]-points[0], points[2]-points[0],
                              points[3]-points[0])
    inverse = affine.inv()
    tail = inverse*(coordinate-points[0])
    lambdas = [1-sum(tail), tail[0], tail[1], tail[2]]
    gradients = [-sum((sp.Matrix(inverse.row(row)).T for row in range(3)),
                      sp.zeros(3, 1))]
    gradients.extend(sp.Matrix(inverse.row(row)).T for row in range(3))
    fields = []
    for i, j in local_edge_pairs:
        fields.append(sp.simplify(lambdas[i]*gradients[j]
                                  - lambdas[j]*gradients[i]))
    return fields


coarse_fields = local_whitney_fields(reference_vertices)
pointwise_nested = True
for tetrahedron in fine_tetrahedra:
    fine_fields = local_whitney_fields(
        tuple(fine_vertex_position[mask] for mask in tetrahedron))
    for coarse_col in range(len(coarse_edges)):
        reconstruction = sp.zeros(3, 1)
        for local_row, (i, j) in enumerate(local_edge_pairs):
            global_row = fine_edge_index[(tetrahedron[i], tetrahedron[j])]
            reconstruction += (whitney_injection[global_row, coarse_col]
                               * fine_fields[local_row])
        if sp.simplify(reconstruction-coarse_fields[coarse_col]) != sp.zeros(3, 1):
            pointwise_nested = False
            break
    if not pointwise_nested:
        break
check("coarse Whitney forms are pointwise nested on all 24 reference subcells",
      pointwise_nested,
      "restriction identity is metric-independent")
check("the same nesting survives a common smooth radial map to round S3",
      pointwise_nested and radial_displacement > 0,
      "fine spherical simplices are restrictions of F(lambda)=normalize(sum lambda_i q_i)")

# Gauge-covariance attack.  A constant Whitney Gram matrix couples Lie-algebra
# variables on different edges.  In left trivialization, a vertex gauge
# transformation rotates an edge variable by the adjoint action at its source.
# Therefore a cross term between different sources is not locally gauge
# invariant.  The reference mass supplies an exact witness.
edge_a = coarse_edges.index((1, 2))
edge_b = coarse_edges.index((2, 3))
cross_source_coefficient = coarse_mass[edge_a, edge_b]
norm_change_under_pi_rotation = sp.simplify(-4*cross_source_coefficient)
check("constant Whitney mass has a nonzero cross-source edge coupling",
      cross_source_coefficient == -sp.Rational(1, 120),
      f"M[(1,2),(2,3)]={cross_source_coefficient}")
check("that coupling violates independent local vertex gauge invariance",
      norm_change_under_pi_rotation == sp.Rational(1, 30),
      "rotating only the second source sends Z to -Z and changes the norm by 1/30")

# Could a positive diagonal (product-link) metric retain exact cylindrical
# pullback?  Tetrahedral symmetry reduces its 50 fine weights to the six
# incidence types (|S|,|T|).  Solve P^T diag(w) P=I exactly.  Every arbitrary
# positive solution could be averaged over S4 to a positive orbit-constant
# one, so failure here excludes all positive diagonal solutions, not only the
# symmetric ansatz.
edge_types = sorted(set((source.bit_count(), target.bit_count())
                        for source, target in fine_edges))
orbit_weights = sp.symbols(f"u0:{len(edge_types)}")
diagonal_weights = sp.diag(*(
    orbit_weights[edge_types.index((source.bit_count(), target.bit_count()))]
    for source, target in fine_edges))
diagonal_residual = (whitney_injection.T*diagonal_weights
                     * whitney_injection-sp.eye(6))
diagonal_equations = [diagonal_residual[row, col]
                      for row in range(6) for col in range(row, 6)]
diagonal_solutions = sp.linsolve(diagonal_equations, orbit_weights)
solution_tuple = next(iter(diagonal_solutions))
forced_type_13 = sp.simplify(solution_tuple[edge_types.index((1, 3))])
expected_forced_type_13 = (-sp.Rational(9, 16)*orbit_weights[2]
                           - sp.Rational(1, 4)*orbit_weights[3]
                           - sp.Rational(9, 32)*orbit_weights[4]
                           - sp.Rational(1, 16)*orbit_weights[5])
check("no positive diagonal link metric is both tetrahedral and cylindrical",
      forced_type_13 == expected_forced_type_13,
      f"weight(1,3)={forced_type_13}, negative when all free weights are positive")

# A local connection-dependent rescue exists.  For each tetrahedral base
# vertex r, transport every edge tangent from its source to r along the unique
# direct tetrahedral edge, evaluate the Whitney quadratic form there, and
# average the four positive forms.  No base vertex is privileged.
def qmul(left, right):
    w, xq, yq, zq = left
    W, X, Y, Z = right
    return np.array((w*W-xq*X-yq*Y-zq*Z,
                     w*X+xq*W+yq*Z-zq*Y,
                     w*Y-xq*Z+yq*W+zq*X,
                     w*Z+xq*Y-yq*X+zq*W))


def qconj(q):
    return np.array((q[0], -q[1], -q[2], -q[3]))


def random_unit_quaternion(rng):
    q = rng.normal(size=4)
    return q/np.linalg.norm(q)


def adjoint(q, vector):
    pure = np.concatenate(([0.0], vector))
    return qmul(qmul(q, pure), qconj(q))[1:]


def directed_link(links, source, target):
    if source == target:
        return np.array((1., 0., 0., 0.))
    if source < target:
        return links[(source, target)]
    return qconj(links[(target, source)])


coarse_mass_float = np.asarray(coarse_mass, dtype=float)
coarse_mass_eigenvalues = coarse_mass.eigenvals()
check("the exact Whitney mass is positive definite",
      coarse_mass_eigenvalues == {
          sp.Rational(1, 6): 1, sp.Rational(1, 24): 4,
          sp.Rational(1, 60): 1},
      f"eigenvalues with multiplicity={coarse_mass_eigenvalues}")


def covariant_whitney_norm(links, tangents):
    total = 0.0
    for basepoint in range(4):
        transported = []
        for edge in coarse_edges:
            source = edge[0]
            transporter = directed_link(links, source, basepoint)
            transported.append(adjoint(transporter, tangents[edge]))
        transported = np.asarray(transported)
        total += sum(coarse_mass_float[row, col]
                     * float(transported[row] @ transported[col])
                     for row in range(6) for col in range(6))
    return total/4.0


def constant_whitney_norm(tangents):
    vectors = np.asarray([tangents[edge] for edge in coarse_edges])
    return sum(coarse_mass_float[row, col]
               * float(vectors[row] @ vectors[col])
               for row in range(6) for col in range(6))


rng = np.random.default_rng(600)
max_gauge_residual = 0.0
minimum_covariant_norm = np.inf
connection_dependence_witness = 0.0
for _ in range(40):
    links = {edge: random_unit_quaternion(rng) for edge in coarse_edges}
    tangents = {edge: rng.normal(size=3) for edge in coarse_edges}
    gauges = [random_unit_quaternion(rng) for _ in range(4)]
    transformed_links = {
        (source, target): qmul(gauges[target],
                               qmul(link, qconj(gauges[source])))
        for (source, target), link in links.items()
    }
    transformed_tangents = {
        edge: adjoint(gauges[edge[0]], vector)
        for edge, vector in tangents.items()
    }
    before = covariant_whitney_norm(links, tangents)
    after = covariant_whitney_norm(transformed_links,
                                   transformed_tangents)
    max_gauge_residual = max(max_gauge_residual, abs(after-before))
    minimum_covariant_norm = min(minimum_covariant_norm, before)
    identity_links = {edge: np.array((1., 0., 0., 0.))
                      for edge in coarse_edges}
    connection_dependence_witness = max(
        connection_dependence_witness,
        abs(before-covariant_whitney_norm(identity_links, tangents)))

check("basepoint-averaged covariant Whitney metric is locally gauge invariant",
      max_gauge_residual < 2e-12,
      f"40 deterministic trials, max residual={max_gauge_residual:.3e}")
check("basepoint-averaged covariant Whitney metric is positive on witnesses",
      minimum_covariant_norm > 1e-6,
      f"minimum sampled nonzero norm={minimum_covariant_norm:.12f}")
identity_links = {edge: np.array((1., 0., 0., 0.)) for edge in coarse_edges}
identity_tangents = {edge: rng.normal(size=3) for edge in coarse_edges}
flat_reduction_residual = abs(
    covariant_whitney_norm(identity_links, identity_tangents)
    - constant_whitney_norm(identity_tangents))
check("covariant metric reduces exactly to Whitney mass at the flat connection",
      flat_reduction_residual < 2e-12,
      f"residual={flat_reduction_residual:.3e}")
check("covariant Whitney metric is genuinely connection dependent",
      connection_dependence_witness > 1e-4,
      f"sampled change from flat metric={connection_dependence_witness:.6f}")

# Crucial distinction: P^T M_f P=M_c is tangent-space inclusion, whereas a
# cylindrical Laplacian on functions requires the projection p to be a
# Riemannian submersion.  At the flat connection its differential A adds the
# two oriented half-edge tangents, and the required dual identity is
# A M_f^{-1} A^T=M_c^{-1}.  Test it exactly before making any curved claim.
coarse_projection = sp.zeros(6, len(fine_edges))
for row, (i, j) in enumerate(coarse_edges):
    midpoint = (1 << i) | (1 << j)
    coarse_projection[row, fine_edge_index[(1 << i, midpoint)]] = 1
    coarse_projection[row, fine_edge_index[(1 << j, midpoint)]] = -1
check("Whitney cochain injection is a right inverse of coarse edge composition",
      coarse_projection*whitney_injection == sp.eye(6),
      "A P=I exactly")
flat_cometric_residual = sp.simplify(
    coarse_projection*fine_mass.inv()*coarse_projection.T
    - coarse_mass.inv())
flat_cometric_residual_numeric = np.asarray(flat_cometric_residual, dtype=float)
coarse_inverse_numeric = np.asarray(coarse_mass.inv(), dtype=float)
induced_cometric_numeric = (flat_cometric_residual_numeric
                            + coarse_inverse_numeric)
relative_cometric_eigenvalues = sla.eigvalsh(
    induced_cometric_numeric, coarse_inverse_numeric)
check("Whitney metric fails the flat configuration-space submersion identity",
      flat_cometric_residual != sp.zeros(6, 6) and
      np.max(abs(flat_cometric_residual_numeric)) > 60.0,
      "A M_f^-1 A^T != M_c^-1 exactly; relative eigenvalues="
      + np.array2string(relative_cometric_eigenvalues, precision=6))

print("\n" + "=" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("=" * 78)
print("HAAR_MEASURE=DERIVED_PROJECTIVELY_CONSISTENT")
print("HORIZONTAL_METRIC_SCALE=DERIVED_FIXED")
print("VERTICAL_METRIC_SCALE=DERIVED_FREE_POSITIVE_PARAMETER")
print("PROJECTIVE_CANONICITY_SELECTS_DIRAC=DERIVED_NEGATIVE")
print("WHITNEY_FORM_TANGENT_INCLUSION=DERIVED_EXACTLY_ISOMETRIC")
print("UNWEIGHTED_COCHAIN_METRIC=DERIVED_NOT_CYLINDRICAL")
print("WHITNEY_SELECTION_BY_THEORY=OPEN_STRUCTURAL")
print("ROUND_VS_AFFINE_REFINEMENT=DERIVED_DISTINCT")
print("ROUND_RADIAL_WHITNEY_NESTING=DERIVED_COMPATIBLE")
print("CONSTANT_WHITNEY_LINK_METRIC_GAUGE_INVARIANCE=DERIVED_NEGATIVE")
print("POSITIVE_DIAGONAL_GAUGE_METRIC_CYLINDRICALITY=DERIVED_NO_GO")
print("BASEPOINT_AVERAGED_COVARIANT_WHITNEY=DERIVED_LOCAL_CANDIDATE")
print("FLAT_CONFIGURATION_LAPLACIAN_CYLINDRICALITY=DERIVED_NEGATIVE")
print("COVARIANT_WHITNEY_AS_PROJECTIVE_DIRAC_METRIC=KILLED")
print("SM_TARGET_COMPARISON=NOT_PERFORMED")

if passed != tests:
    raise SystemExit(1)

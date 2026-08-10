#!/usr/bin/env python3
"""Geometry and first representation gates for the HD route on the 600-cell.

The fixed 600-cell boundary is treated only as the round spatial slice S^3.
Two transports are kept sharply separate:

* the group-difference transport, which is a pure gauge control;
* the spin lift of Levi--Civita parallel transport along the short geodesic
  edges of the round S^3.

No Standard-Model target or matter character is used in this verifier.
"""

from collections import defaultdict
from itertools import permutations, product

import numpy as np


tests = passed = 0


def check(name, ok, detail=""):
    global tests, passed
    tests += 1
    ok = bool(ok)
    passed += int(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


def permutation_sign(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p))
                     for j in range(i + 1, len(p))) % 2 else 1


def qmul(a, b):
    w, x, y, z = a
    W, X, Y, Z = b
    return np.array((w*W-x*X-y*Y-z*Z,
                     w*X+x*W+y*Z-z*Y,
                     w*Y-x*Z+y*W+z*X,
                     w*Z+x*Y-y*X+z*W), dtype=float)


def qconj(q):
    return np.array((q[0], -q[1], -q[2], -q[3]), dtype=float)


def qsqrt_principal(q):
    """Principal square root of a unit quaternion away from -1."""
    den = np.sqrt(2.0 * (1.0 + q[0]))
    return np.array((1.0 + q[0], q[1], q[2], q[3])) / den


def qnorm(q):
    return float(np.linalg.norm(q))


def qcomm(a, b):
    return qmul(a, b) - qmul(b, a)


def build_600cell():
    sqrt5 = np.sqrt(5.0)
    phi = (1.0 + sqrt5) / 2.0
    vertices = set()
    for i in range(4):
        for sign in (-1.0, 1.0):
            q = [0.0] * 4
            q[i] = sign
            vertices.add(tuple(q))
    for signs in product((-0.5, 0.5), repeat=4):
        vertices.add(signs)
    base = [phi / 2.0, 0.5, 1.0 / (2.0 * phi), 0.0]
    even_perms = [p for p in permutations(range(4))
                  if permutation_sign(p) == 1]
    for p in even_perms:
        q = [base[p[i]] for i in range(4)]
        nz = [i for i, value in enumerate(q) if abs(value) > 1e-12]
        for signs in product((-1, 1), repeat=3):
            r = q[:]
            for i, sign in zip(nz, signs):
                r[i] *= sign
            vertices.add(tuple(round(value, 10) for value in r))
    V = np.asarray(sorted(vertices), dtype=float)
    dots = V @ V.T
    edges = [(i, j) for i in range(120) for j in range(i + 1, 120)
             if abs(dots[i, j] - phi / 2.0) < 1e-8]
    adjacency = defaultdict(set)
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    faces = []
    for i, j in edges:
        for k in adjacency[i] & adjacency[j]:
            if j < k:
                faces.append((i, j, k))
    tetrahedra = []
    for i, j, k in faces:
        for ell in adjacency[i] & adjacency[j] & adjacency[k]:
            if k < ell:
                tetrahedra.append((i, j, k, ell))
    return phi, V, edges, faces, tetrahedra


print("=" * 78)
print("HD-600 GATES 0-1: SPIN HOLONOMY AND NATURAL HD REPRESENTATION")
print("=" * 78)

phi, V, edges, faces, tetrahedra = build_600cell()
check("600-cell has f-vector (120,720,1200,600)",
      (len(V), len(edges), len(faces), len(tetrahedra)) ==
      (120, 720, 1200, 600),
      f"counts={(len(V), len(edges), len(faces), len(tetrahedra))}")
check("all vertices are unit quaternions",
      np.max(abs(np.sum(V*V, axis=1) - 1.0)) < 2e-10)
check("every edge is the unique short geodesic with cos(theta)=phi/2",
      max(abs(float(V[i] @ V[j]) - phi/2.0) for i, j in edges) < 2e-10)


def flat_link(i, j):
    """Column convention: inverse of the Wilson link q_i^{-1}q_j."""
    return qmul(qconj(V[j]), V[i])


def lc_left_between(qi, qj):
    """Spin Levi--Civita transport qi->qj in the global left frame."""
    relative = qmul(qconj(qi), qj)
    return qconj(qsqrt_principal(relative))


def lc_left_link(i, j):
    """Spin Levi--Civita transport i->j in the global left frame."""
    return lc_left_between(V[i], V[j])


def lc_right_link(i, j):
    """The same geometric transport written in the global right frame."""
    relative = qmul(V[j], qconj(V[i]))
    return qsqrt_principal(relative)


def face_holonomy(link, face):
    """Column transport around i->j->k->i, based at i."""
    i, j, k = face
    return qmul(link(k, i), qmul(link(j, k), link(i, j)))


def point_face_holonomy(qi, qj, qk):
    """Levi--Civita holonomy around three arbitrary short-geodesic points."""
    return qmul(lc_left_between(qk, qi),
                qmul(lc_left_between(qj, qk),
                     lc_left_between(qi, qj)))


flat_inverse_residual = max(
    qnorm(qmul(flat_link(i, j), flat_link(j, i)) - np.array((1., 0., 0., 0.)))
    for i, j in edges)
flat_face_residual = max(
    qnorm(face_holonomy(flat_link, face) - np.array((1., 0., 0., 0.)))
    for face in faces)
check("group-difference links reverse by inversion",
      flat_inverse_residual < 2e-9,
      f"max residual={flat_inverse_residual:.3e}")
check("group-difference transport is pure gauge on every face",
      flat_face_residual < 3e-9,
      f"max face residual={flat_face_residual:.3e}")

lc_unit_residual = max(abs(qnorm(lc_left_link(i, j)) - 1.0)
                       for i, j in edges)
lc_inverse_residual = max(
    qnorm(qmul(lc_left_link(i, j), lc_left_link(j, i))
          - np.array((1., 0., 0., 0.)))
    for i, j in edges)
check("Levi--Civita spin links are unit quaternions",
      lc_unit_residual < 2e-10,
      f"max residual={lc_unit_residual:.3e}")
check("Levi--Civita spin links reverse by inversion",
      lc_inverse_residual < 2e-9,
      f"max residual={lc_inverse_residual:.3e}")

# The left and right invariant orthonormal frames satisfy
# y_R=Ad_q(y_L).  Their spin-frame gauge transformation is q at vertex q,
# hence T_R(i,j)=q_j T_L(i,j) q_i^{-1}.
frame_gauge_residual = max(
    qnorm(lc_right_link(i, j)
          - qmul(V[j], qmul(lc_left_link(i, j), qconj(V[i]))))
    for i, j in edges)
check("left- and right-frame Levi--Civita links are exactly gauge-related",
      frame_gauge_residual < 3e-9,
      f"max residual={frame_gauge_residual:.3e}")

holonomies = np.asarray([face_holonomy(lc_left_link, face) for face in faces])
scalar_parts = holonomies[:, 0]
vector_norms = np.linalg.norm(holonomies[:, 1:], axis=1)
check("Levi--Civita face holonomy is nontrivial on all 1200 faces",
      np.min(vector_norms) > 1e-6,
      f"min |Im H_f|={np.min(vector_norms):.12f}")
check("all equilateral faces have one holonomy conjugacy class",
      np.ptp(scalar_parts) < 3e-9 and np.ptp(vector_norms) < 3e-9,
      f"Re(H)={np.mean(scalar_parts):.12f}, |Im(H)|={np.mean(vector_norms):.12f}")

# Each face is a spherical equilateral triangle of side pi/5.  The spherical
# cosine rule gives cos(alpha)=1/sqrt(5) for its interior angle and spherical
# excess E=3 alpha-pi.  Levi--Civita tangent holonomy rotates by E, so its
# spin lift has scalar/vector norms cos(E/2), sin(E/2).  Half-angle algebra
# gives the exact radicals below.
expected_scalar = np.sqrt((25.0 + 11.0*np.sqrt(5.0)) / 50.0)
expected_vector = np.sqrt((25.0 - 11.0*np.sqrt(5.0)) / 50.0)
check("face holonomy equals the exact spherical-excess spin lift",
      np.max(abs(scalar_parts-expected_scalar)) < 3e-9 and
      np.max(abs(vector_norms-expected_vector)) < 3e-9,
      "Re(H)^2=(25+11 sqrt(5))/50; |Im(H)|^2=(25-11 sqrt(5))/50")

# A nonidentity SU(2) element has no fixed spinor unless its scalar part is
# one.  Thus the round connection has no parallel spinor.  Equivalently, its
# links cannot define a flat local-system coboundary: the square around every
# triangular 2-cell contains H_f-I rather than zero.
fixed_spinor_determinants = 2.0 * (1.0 - scalar_parts)
check("round Levi--Civita transport has no parallel spinor",
      np.min(abs(fixed_spinor_determinants)) > 1e-6,
      f"min |det(H_f-I)|={np.min(abs(fixed_spinor_determinants)):.12f}")
check("curvature obstructs using these links as a flat cochain twist",
      np.min(np.linalg.norm(holonomies-np.array((1., 0., 0., 0.)), axis=1))
      > 1e-6,
      "d_A^2 contains H_f-I on every triangular face")

# Gauge-invariant cross-check: face traces agree in the right frame.
right_scalars = np.asarray([
    face_holonomy(lc_right_link, face)[0] for face in faces
])
check("face holonomy traces are independent of the invariant frame",
      np.max(abs(right_scalars - scalar_parts)) < 4e-9,
      f"max trace/2 residual={np.max(abs(right_scalars-scalar_parts)):.3e}")

# A non-scalar quaternion generates only a commutative complex algebra.  Two
# noncommuting face holonomies generate H over R and therefore M2(C) after
# complexification.  Find the deterministic first witness and certify rank 4.
witness = None
for a in range(len(holonomies)):
    for b in range(a + 1, len(holonomies)):
        if qnorm(qcomm(holonomies[a], holonomies[b])) > 1e-8:
            witness = (a, b)
            break
    if witness is not None:
        break
if witness is None:
    algebra_rank = 0
    commutator_norm = 0.0
else:
    a, b = witness
    h1, h2 = holonomies[a], holonomies[b]
    algebra_basis = np.vstack((np.array((1., 0., 0., 0.)), h1, h2,
                               qmul(h1, h2)))
    algebra_rank = int(np.linalg.matrix_rank(algebra_basis, tol=1e-10))
    commutator_norm = qnorm(qcomm(h1, h2))
check("two face holonomies are noncommuting",
      witness is not None,
      f"faces={witness}, commutator norm={commutator_norm:.12f}")
check("face holonomies generate H over R, hence M2(C) as a complex *-algebra",
      algebra_rank == 4,
      f"real quaternion-algebra rank={algebra_rank}")

# The motivating HD construction acts on S+S.  The representation selected by
# ordinary spin geometry is the diagonal double U -> diag(U,U).  Test the
# paper's state-separation gate directly: this algebra must not be confused
# with a four-dimensional irreducible representation chosen by hand.
def quaternion_matrix(q):
    w, x, y, z = q
    return np.array(((w+1j*x, y+1j*z),
                     (-y+1j*z, w-1j*x)), dtype=complex)


spin_matrices = [quaternion_matrix(h) for h in holonomies]
doubled_matrices = [np.kron(np.eye(2), u) for u in spin_matrices]
doubled_algebra_rank = int(np.linalg.matrix_rank(
    np.asarray([u.reshape(-1) for u in doubled_matrices]), tol=1e-9))
check("natural diagonal action on S+S remains M2(C), not M4(C)",
      doubled_algebra_rank == 4,
      f"complex matrix-algebra rank={doubled_algebra_rank} (M4 would be 16)")

# Compute the commutant of two noncommuting doubled holonomies by solving
# [X,A]=0 on the 16-dimensional complex matrix space.
commutator_columns = []
for row in range(4):
    for col in range(4):
        matrix_unit = np.zeros((4, 4), dtype=complex)
        matrix_unit[row, col] = 1.0
        commutator_columns.append(np.concatenate((
            (matrix_unit-doubled_matrices[a] @ matrix_unit
             @ doubled_matrices[a].conj().T).reshape(-1),
            (matrix_unit-doubled_matrices[b] @ matrix_unit
             @ doubled_matrices[b].conj().T).reshape(-1))))
commutator_system = np.asarray(commutator_columns).T
doubled_commutant_dimension = 16-int(np.linalg.matrix_rank(
    commutator_system, tol=1e-9))
check("natural doubled holonomy algebra has a four-dimensional commutant",
      doubled_commutant_dimension == 4,
      f"complex commutant dimension={doubled_commutant_dimension}")

# Two pure states with the same spinor and different copy labels have equal
# expectation value on every diagonal holonomy, so the algebra cannot
# separate them.  This is the concrete finite-carrier version of the
# Stone--Weierstrass obstruction stated in arXiv:2504.03391.
copy_one = np.array((1., 0., 0., 0.), dtype=complex)
copy_two = np.array((0., 0., 1., 0.), dtype=complex)
separation_residual = max(abs(
    np.vdot(copy_one, u @ copy_one)-np.vdot(copy_two, u @ copy_two))
    for u in doubled_matrices)
check("diagonal spin doubling fails to separate the two copy states",
      separation_residual < 2e-12,
      f"max expectation difference={separation_residual:.3e}")

# Levi--Civita spin transport is SU(2)-valued.  The U(1) connection required
# to enlarge it to the U(2) configuration space used in the 2025 HD paper is
# additional data, not hidden in the quaternion construction.
determinant_residual = max(abs(np.linalg.det(u)-1.0) for u in spin_matrices)
check("Levi--Civita links and holonomies contain no selected U(1) factor",
      determinant_residual < 3e-10,
      f"max |det(H)-1|={determinant_residual:.3e}")

# Orientation-preserving round isometries q -> a q b^{-1} change the left
# frame by the constant gauge b.  It is enough to test all a,b in the exact
# vertex group on a deterministic sample of links; the formula itself proves
# the general statement.
sample_edges = edges[::47]
symmetry_residual = 0.0
for a in V[::11]:
    for b in V[::13]:
        binv = qconj(b)
        for i, j in sample_edges:
            qi = qmul(a, qmul(V[i], binv))
            qj = qmul(a, qmul(V[j], binv))
            transformed = qconj(qsqrt_principal(qmul(qconj(qi), qj)))
            expected = qmul(b, qmul(lc_left_link(i, j), binv))
            symmetry_residual = max(symmetry_residual,
                                    qnorm(transformed - expected))
check("Levi--Civita links are Spin(4)-equivariant up to constant fiber gauge",
      symmetry_residual < 4e-9,
      f"sampled max residual={symmetry_residual:.3e}")

# First canonical refinement check.  Radially normalized cell barycentres
# give 120+720+1200+600=2640 points, the same carrier size as the oriented
# cochain/Kahler--Dirac construction.  The cover-relation (Hasse) graph has
# 2E+3F+4T=7440 edges.  This is deliberately distinguished from the full
# barycentric 1-skeleton, which also joins non-consecutive comparable cells.
def normalized_sum(indices):
    q = np.sum(V[list(indices)], axis=0)
    return q / qnorm(q)


edge_centres = [normalized_sum(edge) for edge in edges]
face_centres = [normalized_sum(face) for face in faces]
tetra_centres = [normalized_sum(tetra) for tetra in tetrahedra]
all_cell_centres = np.vstack((V, edge_centres, face_centres, tetra_centres))
distinct_cell_centres = len({tuple(np.round(q, 9)) for q in all_cell_centres})
hasse_edges = 2*len(edges) + 3*len(faces) + 4*len(tetrahedra)
check("normalized cell barycentres give the 2640-point Hasse carrier",
      len(all_cell_centres) == distinct_cell_centres == 2640 and
      hasse_edges == 7440,
      f"points={len(all_cell_centres)}, distinct={distinct_cell_centres}, "
      f"cover relations={hasse_edges}")

# Along a coarse geodesic, transport must factor through its normalized
# midpoint.  This is exact in the one-parameter subgroup, independently of
# the two-dimensional curvature test below.
edge_refinement_residual = 0.0
for i, j in edges:
    midpoint = normalized_sum((i, j))
    composite = qmul(lc_left_between(midpoint, V[j]),
                     lc_left_between(V[i], midpoint))
    edge_refinement_residual = max(
        edge_refinement_residual,
        qnorm(composite-lc_left_between(V[i], V[j])))
check("coarse edge transport factors through every geodesic midpoint",
      edge_refinement_residual < 3e-9,
      f"max residual={edge_refinement_residual:.3e}")

# Each spherical face is divided by its edge midpoints and normalized face
# centre into six congruent geodesic triangles.  Each therefore carries one
# sixth of the coarse spherical excess E and spin angle E/12.
subface_scalars = []
subface_vectors = []
for i, j, k in faces:
    qi, qj, qk = V[i], V[j], V[k]
    mij = normalized_sum((i, j))
    mjk = normalized_sum((j, k))
    mki = normalized_sum((k, i))
    centre = normalized_sum((i, j, k))
    subfaces = ((qi, mij, centre), (qi, centre, mki),
                (qj, mjk, centre), (qj, centre, mij),
                (qk, mki, centre), (qk, centre, mjk))
    for a, b, c in subfaces:
        h = point_face_holonomy(a, b, c)
        subface_scalars.append(h[0])
        subface_vectors.append(qnorm(h[1:]))

alpha = np.arccos(1.0/np.sqrt(5.0))
spherical_excess = 3.0*alpha - np.pi
expected_sub_scalar = np.cos(spherical_excess/12.0)
expected_sub_vector = np.sin(spherical_excess/12.0)
subface_scalars = np.asarray(subface_scalars)
subface_vectors = np.asarray(subface_vectors)
check("7200 refined faces carry exactly one sixth of the coarse curvature",
      np.max(abs(subface_scalars-expected_sub_scalar)) < 4e-9 and
      np.max(abs(subface_vectors-expected_sub_vector)) < 4e-9,
      f"Re(H_small)={np.mean(subface_scalars):.12f}, "
      f"|Im(H_small)|={np.mean(subface_vectors):.12f}")

print("\n" + "=" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("=" * 78)
print("CONTROL_FLAT=DERIVED_NEGATIVE")
print("LC_FACE_HOLONOMY=DERIVED_NONTRIVIAL")
print("LC_GENERATED_FIBER_ALGEBRA=M2(C)")
print("BARYCENTRIC_REFINEMENT=DERIVED_CONSISTENT_AT_LEVEL_ONE")
print("KAEHLER_DIRAC_FLAT_TWIST=DERIVED_NEGATIVE")
print("NATURAL_HD_SPIN_DOUBLING=DERIVED_NEGATIVE_FOR_STATE_SEPARATION")
print("U1_CONNECTION=OPEN")
print("SM_TARGET_COMPARISON=NOT_PERFORMED")

if passed != tests:
    raise SystemExit(1)

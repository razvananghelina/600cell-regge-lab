"""Exact low-mode vector sampling theorem on the icosahedral Hopf base.

Continuum 1-forms are sampled by oriented geodesic edge integration.  For
real spherical harmonics with l=1,2, the script verifies the gradient and
Hodge-rotated families, the discrete Hodge decomposition, a band-limited
Hodge-Laplacian intertwiner, and the residual inner-product normalization
freedom.

Numerical residuals certify finite identities; no refinement limit or
nonabelian bracket is asserted.
"""

import math
import os
import sys

import numpy as np
from numpy.linalg import eigh, matrix_rank, norm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0
TOL = 2.0e-9
failed = False
tests_run = 0
tests_passed = 0


def check(name, condition, detail=""):
    global failed, tests_run, tests_passed
    tests_run += 1
    if condition:
        tests_passed += 1
        print(f"  [PASS] {name}")
    else:
        failed = True
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0],
    ])


def find_idx(v, verts, tol=1.0e-6):
    dots = verts @ v
    idx = int(np.argmax(dots))
    return idx if dots[idx] > 1.0 - tol else -1


def find_fibers(verts):
    for generator in verts:
        if abs(generator[0] - PHI / 2.0) >= 1.0e-6:
            continue
        power = generator.copy()
        valid = True
        for k in range(2, 11):
            power = qmul(power, generator)
            if k == 5 and not np.allclose(power, [-1, 0, 0, 0], atol=1.0e-6):
                valid = False
                break
            if k == 10 and not np.allclose(power, [1, 0, 0, 0], atol=1.0e-6):
                valid = False
        if not valid:
            continue
        subgroup = []
        power = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_idx(power, verts))
            power = qmul(power, generator)
        fibers, used = [], set()
        for seed in range(len(verts)):
            if seed in used:
                continue
            fiber = []
            for element in subgroup:
                idx = find_idx(qmul(verts[seed], verts[element]), verts)
                if idx >= 0 and idx not in used:
                    fiber.append(idx)
                    used.add(idx)
            if len(fiber) == 10:
                fibers.append(fiber)
        if len(fibers) == 12:
            return fibers
    raise RuntimeError("No Hopf fibration found")


def base_complex(adj, fibers, xyz):
    vertex_to_fiber = {v: fi for fi, fiber in enumerate(fibers) for v in fiber}
    base_adj = np.zeros((12, 12))
    for i in range(len(adj)):
        for j in range(i + 1, len(adj)):
            if adj[i, j] > 0.5 and vertex_to_fiber[i] != vertex_to_fiber[j]:
                fi, fj = vertex_to_fiber[i], vertex_to_fiber[j]
                base_adj[fi, fj] = base_adj[fj, fi] = 1.0
    edges = [(i, j) for i in range(12) for j in range(i + 1, 12) if base_adj[i, j] > 0.5]
    edge_index = {edge: k for k, edge in enumerate(edges)}
    faces = []
    for i in range(12):
        for j in range(i + 1, 12):
            for k in range(j + 1, 12):
                if base_adj[i, j] and base_adj[i, k] and base_adj[j, k]:
                    face = (i, j, k)
                    if np.linalg.det(np.column_stack([xyz[i], xyz[j], xyz[k]])) < 0:
                        face = (i, k, j)
                    faces.append(face)
    d0 = np.zeros((30, 12))
    for e, (i, j) in enumerate(edges):
        d0[e, i], d0[e, j] = -1.0, 1.0
    d1 = np.zeros((20, 30))
    for f, (i, j, k) in enumerate(faces):
        for a, b in ((i, j), (j, k), (k, i)):
            lo, hi = min(a, b), max(a, b)
            d1[f, edge_index[(lo, hi)]] = 1.0 if a < b else -1.0
    return edges, d0, d1


def scalar_values(xyz):
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    return np.column_stack([x, y, z, x*y, y*z, z*x, x*x-y*y, 2*z*z-x*x-y*y])


def polynomial_gradients(r):
    x, y, z = r
    return np.array([
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [y, x, 0], [0, z, y], [z, 0, x],
        [2*x, -2*y, 0], [-2*x, -2*y, 4*z],
    ], dtype=float)


def curl_edge_integrals(xyz, edges):
    # Gauss-Legendre integration is exact to displayed tolerance for these
    # trigonometric polynomials of degree at most two.
    nodes, weights = np.polynomial.legendre.leggauss(16)
    sampled = np.zeros((len(edges), 8))
    for edge, (i, j) in enumerate(edges):
        u, v = xyz[i], xyz[j]
        theta = math.acos(float(np.clip(u @ v, -1.0, 1.0)))
        tangent = (v - math.cos(theta)*u) / math.sin(theta)
        for node, weight in zip(nodes, weights):
            t = 0.5*theta*(node + 1.0)
            r = math.cos(t)*u + math.sin(t)*tangent
            dr = -math.sin(t)*u + math.cos(t)*tangent
            grads = polynomial_gradients(r)
            # (*dY)(dr) with the orientation induced by the outward normal.
            sampled[edge] += 0.5*theta*weight*np.einsum("ij,j->i", np.cross(r, grads), dr)
    return sampled


print("Verifying vector sampling/intertwining on the Hopf-base icosahedron...")
verts, adjacency_600, _ = build_600cell()
fibers = find_fibers(verts)

# Obtain the canonical regular embedding from the scalar spectral theorem.
vertex_to_fiber = {v: fi for fi, fiber in enumerate(fibers) for v in fiber}
base_adj = np.zeros((12, 12))
for i in range(120):
    for j in range(i + 1, 120):
        if adjacency_600[i, j] > 0.5 and vertex_to_fiber[i] != vertex_to_fiber[j]:
            fi, fj = vertex_to_fiber[i], vertex_to_fiber[j]
            base_adj[fi, fj] = base_adj[fj, fi] = 1.0
lap0 = np.diag(base_adj.sum(axis=1)) - base_adj
_, scalar_evecs = eigh(lap0)
xyz = scalar_evecs[:, 1:4].copy()
xyz /= np.linalg.norm(xyz, axis=1)[:, None]

edges, d0, d1 = base_complex(adjacency_600, fibers, xyz)
edge_set = set(edges)
# The 2I left action factors through A5 on fibers.  Compute the orbit of one
# base edge to certify that an invariant local diagonal edge metric has one
# weight only.
edge_orbit = set()
for group_element in verts:
    fiber_perm = []
    for fiber in fibers:
        moved = find_idx(qmul(group_element, verts[fiber[0]]), verts)
        fiber_perm.append(vertex_to_fiber[moved])
    a, b = fiber_perm[edges[0][0]], fiber_perm[edges[0][1]]
    edge_orbit.add((min(a, b), max(a, b)))
values = scalar_values(xyz)
sample_exact = d0 @ values
sample_coexact = curl_edge_integrals(xyz, edges)
sampling = np.column_stack([sample_exact[:, :3], sample_exact[:, 3:],
                            sample_coexact[:, :3], sample_coexact[:, 3:]])

b_exact = d0 @ d0.T
c_coexact = d1.T @ d1
lap1 = b_exact + c_coexact

check("icosahedral complex has (V,E,F)=(12,30,20)", len(edges) == 30 and d1.shape == (20, 30))
check("A5 is transitive on all 30 base edges", edge_orbit == edge_set,
      f"orbit size={len(edge_orbit)}")
check("d1 d0 = 0", norm(d1 @ d0) < TOL)
check("d0 has rank 11 and kills the constant gauge amplitude",
      matrix_rank(d0, tol=1.0e-10) == 11 and norm(d0 @ np.ones(12)) < TOL)
check("vector sampling has rank 16=2(3+5)", matrix_rank(sampling, tol=1.0e-10) == 16)
check("gradient sampling is exact", norm(c_coexact @ sample_exact) < TOL)
check("curl sampling is coexact", norm(b_exact @ sample_coexact) < TOL,
      f"||B S_curl||={norm(b_exact @ sample_coexact):.3e}")
check("gradient and curl samples are orthogonal", norm(sample_exact.T @ sample_coexact) < TOL)

# Raw discrete eigenvalues for exact l=1,l=2 and coexact l=1,l=2.
raw_eigenvalues = [5.0-math.sqrt(5.0), 6.0, 3.0-math.sqrt(5.0), 2.0]
blocks = [sample_exact[:, :3], sample_exact[:, 3:], sample_coexact[:, :3], sample_coexact[:, 3:]]
for label, block, eigenvalue in zip(("exact l=1", "exact l=2", "coexact l=1", "coexact l=2"), blocks, raw_eigenvalues):
    check(f"{label} is a Delta_1 eigenspace", norm(lap1 @ block - eigenvalue*block) < TOL)

# The unique cubic p on the four sampled eigenvalues which maps both Hodge
# types to the continuum eigenvalues l(l+1).
sqrt5 = math.sqrt(5.0)
poly = np.array([42.0-24.0*sqrt5, -24.0+22.0*sqrt5,
                 3.0-6.0*sqrt5, sqrt5/2.0])
check("exact cubic maps the four raw eigenvalues to 2,6,2,6",
      np.allclose([sum(poly[k]*x**k for k in range(4)) for x in raw_eigenvalues],
                  [2.0, 6.0, 2.0, 6.0], atol=TOL))
corrected = poly[0]*np.eye(30) + poly[1]*lap1 + poly[2]*(lap1@lap1) + poly[3]*(lap1@lap1@lap1)
continuum_hodge = np.diag([2.0]*3 + [6.0]*5 + [2.0]*3 + [6.0]*5)
intertwiner_error = norm(corrected @ sampling - sampling @ continuum_hodge)
check("p(Delta_1) S = S Delta_Hodge on l=1,2", intertwiner_error < TOL,
      f"residual={intertwiner_error:.3e}, p={poly.tolist()}")

reconstruction = np.linalg.inv(sampling.T @ sampling) @ sampling.T
check("band-limited vector sampling has exact left inverse", norm(reconstruction @ sampling - np.eye(16)) < TOL)
projector = sampling @ reconstruction
check("resolved alias complement has dimension 14", abs(np.trace(np.eye(30)-projector)-14.0) < TOL)
check("alias splits as exact 3 plus coexact 11",
      matrix_rank((np.eye(30)-projector) @ d0, tol=1.0e-9) == 3
      and matrix_rank((np.eye(30)-projector) @ d1.T, tol=1.0e-9) == 11)

# Compare edge Euclidean Gram matrices with exact continuum L2 Gram matrices.
continuum_l1 = (8.0*math.pi/3.0)*np.eye(3)
continuum_l2 = np.diag([8.0*math.pi/5.0]*3 + [32.0*math.pi/5.0, 96.0*math.pi/5.0])
continuum_blocks = [continuum_l1, continuum_l2, continuum_l1, continuum_l2]
scale_factors = []
gram_residuals = []
for block, continuum_gram in zip(blocks, continuum_blocks):
    discrete_gram = block.T @ block
    scale = np.trace(discrete_gram) / np.trace(continuum_gram)
    scale_factors.append(scale)
    gram_residuals.append(norm(discrete_gram-scale*continuum_gram))
check("each irreducible block is conformally sampled", max(gram_residuals) < TOL,
      f"Gram residuals={[float(x) for x in gram_residuals]}")
spread = max(scale_factors)-min(scale_factors)
theta = math.acos(1.0/math.sqrt(5.0))
expected_scales = [(4.0+PHI**-4)/math.pi, 3.0/math.pi,
                   15.0*theta*theta/(4.0*math.pi), (4.0+PHI**-4)/math.pi]
check("the four norm factors match their closed forms",
      np.allclose(scale_factors, expected_scales, atol=TOL),
      f"computed={scale_factors}, expected={expected_scales}")
check("no single edge weight makes all four blocks isometric", spread > 1.0e-3,
      f"discrete/continuum factors={scale_factors}")

print("\nClassification:")
print("  DERIVED: rank-16 edge-integral sampling and exact band-limited Hodge intertwiner")
print("  DERIVED: four inequivalent discrete/continuum norm factors; no common isometry")
print("  OPEN: a canonical sector-dependent renormalization, gauge-kernel lift, and Lie bracket")
print(f"\nTOTAL: {tests_passed}/{tests_run} tests PASSED")
sys.exit(1 if failed else 0)

"""
exp617: Vector-harmonic structure on the Hopf-base icosahedron.

Goal:
  Extend exp616 from scalar harmonics on the Hopf base to discrete 1-forms.
  The natural question is not yet whether the 12 gauge modes ARE the full
  continuum Yang-Mills connection, but whether the base graph already carries
  the correct low-harmonic vector operator structure.

Main checks:
  1. Build the Hopf-base icosahedron graph (12 vertices, 30 edges, 20 faces).
  2. Form the discrete Hodge pieces on base edges:
       B = d0 d0^T   (exact part)
       C = d1^T d1   (coexact part)
       D = B + C     (edge Hodge Laplacian)
  3. Show that gradients of the l=1 and l=2 scalar harmonics from exp616 are
     exact eigen-1-forms of D with eigenvalues 5-sqrt(5) and 6.
  4. Show that curls d1^T(phi_face) built from the dual-face low harmonics are
     coexact eigen-1-forms of D with eigenvalues 3-sqrt(5) and 2.
  5. Conclude that the base edge Hodge Laplacian already carries a clean low
     exact/coexact vector-harmonic sector of dimensions 3+5 and 3+5.

Interpretation:
  This is not yet the nonabelian gauge bracket. It is the next layer of the
  discrete-to-continuum bridge: the Hopf base supports the correct low-mode
  vector operator structure, not only scalar harmonics.
"""

import math
import sys

import numpy as np
from numpy.linalg import eigh, matrix_rank, norm

sys.path.insert(0, ".")
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0
TOL = 1e-8


def qmul(p, q):
    return np.array(
        [
            p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
            p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
            p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
            p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0],
        ]
    )


def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = int(np.argmax(dots))
    return idx if dots[idx] > 1.0 - tol else -1


def find_hopf_fibration(verts):
    nv = len(verts)
    for i in range(nv):
        if abs(verts[i, 0] - PHI / 2.0) >= 1e-6:
            continue

        g = verts[i]
        power = g.copy()
        ok = True
        for k in range(2, 11):
            power = qmul(power, g)
            if k == 5 and not np.allclose(power, [-1, 0, 0, 0], atol=1e-6):
                ok = False
                break
            if k == 10 and not np.allclose(power, [1, 0, 0, 0], atol=1e-6):
                ok = False
        if not ok:
            continue

        subgroup = []
        power = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_idx(power, verts))
            power = qmul(power, g)

        used = set()
        fibers = []
        for s in range(nv):
            if s in used:
                continue
            fiber = []
            for si in subgroup:
                idx = find_idx(qmul(verts[s], verts[si]), verts)
                if idx >= 0 and idx not in used:
                    fiber.append(idx)
                    used.add(idx)
            if len(fiber) == 10:
                fibers.append(fiber)

        if len(fibers) == 12:
            return fibers

    raise RuntimeError("Could not find Hopf fibration")


def build_base_graph(adj, fibers):
    vertex_to_fiber = {}
    for fi, fiber in enumerate(fibers):
        for v in fiber:
            vertex_to_fiber[v] = fi

    weighted = np.zeros((12, 12), dtype=float)
    nv = adj.shape[0]
    for i in range(nv):
        for j in range(i + 1, nv):
            if adj[i, j] > 0.5:
                fi = vertex_to_fiber[i]
                fj = vertex_to_fiber[j]
                if fi != fj:
                    weighted[fi, fj] += 1.0
                    weighted[fj, fi] += 1.0

    adjacency = (weighted > 0).astype(float)
    laplacian = np.diag(np.sum(adjacency, axis=1)) - adjacency
    return adjacency, laplacian


def harmonic_basis_from_embedding(xyz):
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    return np.column_stack(
        [
            np.ones(len(x)),
            x,
            y,
            z,
            x * y,
            y * z,
            z * x,
            x * x - y * y,
            2.0 * z * z - x * x - y * y,
        ]
    )


def oriented_faces(adjacency, xyz):
    faces = []
    seen = set()
    for i in range(12):
        for j in range(i + 1, 12):
            if adjacency[i, j] < 0.5:
                continue
            for k in range(j + 1, 12):
                if adjacency[i, k] > 0.5 and adjacency[j, k] > 0.5:
                    key = (i, j, k)
                    if key in seen:
                        continue
                    seen.add(key)
                    det = np.linalg.det(np.column_stack([xyz[i], xyz[j], xyz[k]]))
                    if det >= 0.0:
                        faces.append((i, j, k))
                    else:
                        faces.append((i, k, j))
    return faces


def build_edge_operators(adjacency, faces):
    edges = []
    edge_to_idx = {}
    for i in range(12):
        for j in range(i + 1, 12):
            if adjacency[i, j] > 0.5:
                edge_to_idx[(i, j)] = len(edges)
                edges.append((i, j))

    e = len(edges)
    f = len(faces)

    d0 = np.zeros((e, 12))
    for idx, (i, j) in enumerate(edges):
        d0[idx, i] = -1.0
        d0[idx, j] = +1.0

    d1 = np.zeros((f, e))
    for f_idx, (i, j, k) in enumerate(faces):
        boundary = [(i, j, +1.0), (j, k, +1.0), (i, k, -1.0)]
        for a, b, sign in boundary:
            e_idx = edge_to_idx[(min(a, b), max(a, b))]
            d1[f_idx, e_idx] = sign if a < b else -sign

    return edges, d0, d1


def subspace_action(op, basis):
    gram = basis.T @ basis
    action = np.linalg.solve(gram, basis.T @ op @ basis)
    residual = norm(op @ basis - basis @ action)
    return action, residual


def main():
    print("=" * 72)
    print("exp617: vector harmonics on the Hopf-base icosahedron")
    print("=" * 72)

    verts, adj, _ = build_600cell()
    fibers = find_hopf_fibration(verts)
    adjacency, laplacian = build_base_graph(adj, fibers)

    evals, evecs = eigh(laplacian)
    xyz = evecs[:, 1:4]
    xyz = xyz / np.linalg.norm(xyz, axis=1)[:, None]
    scalar_basis = harmonic_basis_from_embedding(xyz)

    faces = oriented_faces(adjacency, xyz)
    edges, d0, d1 = build_edge_operators(adjacency, faces)
    face_centers = []
    for i, j, k in faces:
        center = xyz[i] + xyz[j] + xyz[k]
        center = center / norm(center)
        face_centers.append(center)
    face_centers = np.array(face_centers)
    face_basis = harmonic_basis_from_embedding(face_centers)

    print("\n[1] Base Hodge data")
    print(f"  Vertices: 12, edges: {len(edges)}, faces: {len(faces)}")
    print(f"  rank(d0) = {matrix_rank(d0):d}, rank(d1) = {matrix_rank(d1):d}")

    b_exact = d0 @ d0.T
    c_coexact = d1.T @ d1
    delta_1 = b_exact + c_coexact

    be = np.unique(np.round(np.sort(eigh(b_exact)[0]), 8), return_counts=True)
    ce = np.unique(np.round(np.sort(eigh(c_coexact)[0]), 8), return_counts=True)
    de = np.unique(np.round(np.sort(eigh(delta_1)[0]), 8), return_counts=True)
    print(f"  Spec(B) = {be}")
    print(f"  Spec(C) = {ce}")
    print(f"  Spec(Delta_1) = {de}")

    print("\n[2] Exact low vector harmonics from gradients")
    grad_l1 = d0 @ scalar_basis[:, 1:4]
    grad_l2 = d0 @ scalar_basis[:, 4:]

    for name, basis, target in (
        ("grad l=1", grad_l1, 5.0 - math.sqrt(5.0)),
        ("grad l=2", grad_l2, 6.0),
    ):
        action, residual = subspace_action(delta_1, basis)
        internal = np.linalg.eigvalsh((action + action.T) / 2.0)
        print(f"  {name}: internal eigenvalues {np.round(internal, 8)}")
        print(f"      residual ||Delta_1 B - B T|| = {residual:.3e}")
        print(f"      exactness residual ||C B|| = {norm(c_coexact @ basis):.3e}")
        print(f"      target eigenvalue = {target:.8f}")

    print("\n[3] Coexact low vector harmonics from face harmonics")
    curl_l1 = d1.T @ face_basis[:, 1:4]
    curl_l2 = d1.T @ face_basis[:, 4:]

    for name, basis, target in (
        ("coexact l=1", curl_l1, 3.0 - math.sqrt(5.0)),
        ("coexact l=2", curl_l2, 2.0),
    ):
        action, residual = subspace_action(delta_1, basis)
        internal = np.linalg.eigvalsh((action + action.T) / 2.0)
        print(f"  {name}: internal eigenvalues {np.round(internal, 8)}")
        print(f"      residual ||Delta_1 B - B T|| = {residual:.3e}")
        print(f"      coexactness residual ||B B|| = {norm(b_exact @ basis):.3e}")
        print(f"      target eigenvalue = {target:.8f}")

    print("\n[4] Low vector sector dimensions")
    print("  Exact low sector dimensions: 3 + 5 = 8")
    print("  Coexact low sector dimensions: 3 + 5 = 8")
    print("  Total resolved low vector sector on edges: 16 of 30")

    print("\nInterpretation")
    print("  The Hopf-base icosahedron supports not only the exact scalar low")
    print("  harmonics of S^2 up to l=2, but also a clean low exact/coexact")
    print("  vector-harmonic sector on base edges. This is the operator-level")
    print("  continuum bridge relevant to gauge fields. What remains open is")
    print("  the nonabelian completion: how these low 1-form sectors couple to")
    print("  the 12-dimensional gauge skeleton and to a continuum Lie bracket.")


if __name__ == "__main__":
    main()

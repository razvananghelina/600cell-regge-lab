"""
exp616: Gauge continuum map on the Hopf base.

Goal:
  Replace the failed search for an exact discrete Lie bracket by a cleaner
  continuum statement: the 12-dimensional gauge skeleton on the Hopf base is
  the exact low-harmonic truncation of a round S^2 up to l = 2, plus one
  3-dimensional alias sector.

Main checks:
  1. The 12 Hopf fibers define an unweighted icosahedron graph on the base.
  2. The 3-dimensional eigenspace of the smallest nonzero base Laplacian
     gives a canonical spectral embedding into S^2.
  3. Real spherical harmonics built from that embedding satisfy:
       l = 0  -> eigenvalue 0
       l = 1  -> eigenvalue 5 - sqrt(5)
       l = 2  -> eigenvalue 6
     exactly on the discrete graph.
  4. The 9-dimensional sampled harmonic space equals the spectral projector
     onto 1 + 3 + 5. The remaining 3-dimensional eigenspace is the first
     lattice alias sector.
  5. A quadratic polynomial in the graph Laplacian reproduces the continuum
     S^2 Casimir values 0, 2, 6 on the low-harmonic sector exactly.

Interpretation:
  The relevant discrete-to-continuum bridge for the gauge sector is not an
  exact discrete Lie bracket. It is the fact that the Hopf-base icosahedron
  already carries the exact low-l harmonic content of S^2 up to l = 2.
"""

from collections import defaultdict
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

    raise RuntimeError("Could not find a Hopf fibration with 12 fibers")


def build_base_graph(adj, fibers):
    nv = adj.shape[0]
    vertex_to_fiber = {}
    for fi, fiber in enumerate(fibers):
        for v in fiber:
            vertex_to_fiber[v] = fi

    weighted = np.zeros((12, 12), dtype=float)
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
    return weighted, adjacency, laplacian


def harmonic_basis_from_embedding(xyz):
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    return np.column_stack(
        [
            np.ones(12),
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


def projector_from_columns(mat):
    gram = mat.T @ mat
    return mat @ np.linalg.inv(gram) @ mat.T


def subspace_action(op, basis):
    gram = basis.T @ basis
    action = np.linalg.solve(gram, basis.T @ op @ basis)
    residual = norm(op @ basis - basis @ action)
    return action, residual


def main():
    print("=" * 72)
    print("exp616: gauge continuum map on the Hopf base")
    print("=" * 72)

    verts, adj, _ = build_600cell()
    fibers = find_hopf_fibration(verts)
    weighted, adjacency, laplacian = build_base_graph(adj, fibers)

    print("\n[1] Hopf-base graph")
    degrees = np.sum(adjacency, axis=1)
    weight_values = sorted(set(weighted[weighted > 0].astype(int).tolist()))
    print(f"  Degrees: {degrees.astype(int).tolist()}")
    print(f"  Nonzero edge weights between fibers: {weight_values}")
    print(f"  Total base edges: {int(np.sum(adjacency) // 2)}")

    evals, evecs = eigh(laplacian)
    print(f"  Base Laplacian spectrum: {np.round(evals, 8)}")

    lam1 = 5.0 - math.sqrt(5.0)
    lam2 = 6.0
    lam3 = 5.0 + math.sqrt(5.0)
    expect = np.array([0.0] + [lam1] * 3 + [lam2] * 5 + [lam3] * 3)
    print(f"  Matches icosahedron spectrum: {np.allclose(evals, expect, atol=1e-8)}")

    print("\n[2] Canonical spectral embedding into S^2")
    xyz = evecs[:, 1:4]
    row_norms = np.linalg.norm(xyz, axis=1)
    print(
        "  Row-norm spread before normalization: "
        f"{row_norms.min():.12f} .. {row_norms.max():.12f}"
    )
    xyz = xyz / row_norms[:, None]

    dots = []
    for i in range(12):
        for j in range(i + 1, 12):
            dots.append(np.dot(xyz[i], xyz[j]))
    unique_dots = np.unique(np.round(np.sort(dots), 8))
    print(f"  Distinct pairwise dot products: {unique_dots}")
    print("  Expected regular-icosahedron dots: {-1, -1/sqrt(5), +1/sqrt(5)}")

    print("\n[3] Sampled low harmonics on the base")
    harmonics = harmonic_basis_from_embedding(xyz)
    rank = matrix_rank(harmonics, tol=1e-10)
    print(f"  Rank of sampled l<=2 basis: {rank}")

    l0 = harmonics[:, :1]
    l1 = harmonics[:, 1:4]
    l2 = harmonics[:, 4:]

    for name, basis, target in (
        ("l=0", l0, 0.0),
        ("l=1", l1, lam1),
        ("l=2", l2, lam2),
    ):
        action, residual = subspace_action(laplacian, basis)
        internal = np.linalg.eigvalsh((action + action.T) / 2.0)
        print(f"  {name}: internal eigenvalues {np.round(internal, 8)}")
        print(f"       residual ||L B - B T|| = {residual:.3e}")
        print(f"       target eigenvalue = {target:.8f}")

    print("\n[4] Exact low-mode projector and alias sector")
    low_mask = (
        (np.abs(evals) < TOL)
        | (np.abs(evals - lam1) < TOL)
        | (np.abs(evals - lam2) < TOL)
    )
    low_spec = evecs[:, low_mask]
    proj_spec = low_spec @ low_spec.T
    proj_harm = projector_from_columns(harmonics)
    proj_diff = norm(proj_spec - proj_harm)
    alias_dim = int(np.sum(np.abs(evals - lam3) < TOL))
    print(f"  Projector difference ||P_spec - P_harm|| = {proj_diff:.3e}")
    print(f"  Alias-sector dimension (lambda = 5+sqrt(5)): {alias_dim}")

    print("\n[5] Exact discrete representative of the continuum S^2 Casimir")
    coeffs = np.linalg.solve(
        np.array([[lam1 * lam1, lam1], [lam2 * lam2, lam2]], dtype=float),
        np.array([2.0, 6.0], dtype=float),
    )
    alpha, beta = float(coeffs[0]), float(coeffs[1])
    casimir = alpha * (laplacian @ laplacian) + beta * laplacian
    print(f"  Polynomial: Delta_S2^disc = {alpha:.12f} L^2 + {beta:.12f} L")

    for name, basis, target in (("l=0", l0, 0.0), ("l=1", l1, 2.0), ("l=2", l2, 6.0)):
        residual = norm(casimir @ basis - target * basis)
        print(f"  Casimir residual on {name}: {residual:.3e}")

    print("\nInterpretation")
    print("  The Hopf-base icosahedron realizes the exact low-l scalar harmonics")
    print("  of S^2 up to l=2: 1 + 3 + 5 = 9 modes. The remaining 3-dimensional")
    print("  eigenspace is the first discrete alias sector. This gives a clean")
    print("  discrete-to-continuum map for the gauge skeleton without requiring")
    print("  an exact Lie bracket on the finite graph.")


if __name__ == "__main__":
    main()

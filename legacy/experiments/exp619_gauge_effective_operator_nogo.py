"""
exp619: No-go test for a naive effective gauge operator on Hopf-base amplitudes.

Motivation:
  After exp618, the 12 gauge modes factor exactly as

      (fixed alternating fiber polarization) x (12 base amplitudes).

  The obvious next hope is that compressing natural edge operators to this
  amplitude space will produce the Hopf-base Laplacian or another simple
  continuum gauge operator.

This experiment tests two minimal constructions:
  1. First-order compression:      U^T O U
  2. Naive second-order Schur map: -U^T V Box_1^+ V U

where U is the canonical vertical gauge basis and Box_1^+ is the Moore-Penrose
pseudoinverse on the orthogonal complement of the gauge kernel.

Conclusion:
  - At first order, all natural operators collapse to scalars on amplitudes.
  - At second order, nontrivial operators appear, but they are not simple
    functions of the Hopf-base Laplacian and do not preserve the low
    1+3+5+3' block decomposition cleanly.

Interpretation:
  The gauge continuum completion is not obtained by naive operator projection.
  The vertical factorization is exact, but the nonabelian/base dynamics require
  a subtler construction than straightforward compression.
"""

from collections import defaultdict
import math
import sys

import numpy as np
from numpy.linalg import eigh, lstsq, norm

sys.path.insert(0, ".")
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0
A1 = 5.0
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


def build_geometry():
    verts, adj, _ = build_600cell()
    nv = len(verts)

    adj_list = defaultdict(set)
    edges = []
    edge_to_idx = {}
    for i in range(nv):
        for j in range(i + 1, nv):
            if adj[i, j] > 0.5:
                adj_list[i].add(j)
                adj_list[j].add(i)
                edge_to_idx[(i, j)] = len(edges)
                edges.append((i, j))

    fibers = find_hopf_fibration(verts)
    vertex_to_fiber = {}
    for fi, fiber in enumerate(fibers):
        for v in fiber:
            vertex_to_fiber[v] = fi

    fiber_edge_lists = []
    for fiber in fibers:
        fiber_edges = []
        for k in range(10):
            u = fiber[k]
            v = fiber[(k + 1) % 10]
            fiber_edges.append(edge_to_idx[(min(u, v), max(u, v))])
        fiber_edge_lists.append(fiber_edges)

    vte = defaultdict(list)
    for e_idx, (i, j) in enumerate(edges):
        vte[i].append(e_idx)
        vte[j].append(e_idx)

    a_line = np.zeros((len(edges), len(edges)))
    for v in range(nv):
        inc = vte[v]
        for a in range(len(inc)):
            for b in range(a + 1, len(inc)):
                ea = inc[a]
                eb = inc[b]
                a_line[ea, eb] = 1.0
                a_line[eb, ea] = 1.0

    a_fiber = np.zeros((len(edges), len(edges)))
    for fiber_edges in fiber_edge_lists:
        for k in range(10):
            ea = fiber_edges[k]
            eb = fiber_edges[(k + 1) % 10]
            a_fiber[ea, eb] = 1.0
            a_fiber[eb, ea] = 1.0

    l_fiber = np.diag(np.sum(a_fiber, axis=1)) - a_fiber
    l_cross = np.diag(np.sum(a_line - a_fiber, axis=1)) - (a_line - a_fiber)
    box1 = l_cross - A1 * l_fiber

    triangles = []
    for i in range(nv):
        for j in adj_list[i]:
            if j > i:
                for k in adj_list[i] & adj_list[j]:
                    if k > j:
                        triangles.append((i, j, k))

    d0 = np.zeros((len(edges), nv))
    for e_idx, (i, j) in enumerate(edges):
        d0[e_idx, i] = -1.0
        d0[e_idx, j] = +1.0

    d1 = np.zeros((len(triangles), len(edges)))
    for f_idx, (i, j, k) in enumerate(triangles):
        d1[f_idx, edge_to_idx[(i, j)]] = +1.0
        d1[f_idx, edge_to_idx[(j, k)]] = +1.0
        d1[f_idx, edge_to_idx[(i, k)]] = -1.0

    b_exact = d0 @ d0.T
    c_coexact = d1.T @ d1
    delta_1 = b_exact + c_coexact

    weighted = np.zeros((12, 12))
    for i in range(nv):
        for j in range(i + 1, nv):
            if adj[i, j] > 0.5:
                fi = vertex_to_fiber[i]
                fj = vertex_to_fiber[j]
                if fi != fj:
                    weighted[fi, fj] += 1.0
                    weighted[fj, fi] += 1.0

    a0 = (weighted > 0).astype(float)
    l0 = np.diag(np.sum(a0, axis=1)) - a0

    alt = np.array([(-1.0) ** k for k in range(10)], dtype=float)
    alt /= norm(alt)
    u = np.zeros((len(edges), 12))
    for fi, fiber_edges in enumerate(fiber_edge_lists):
        u[fiber_edges, fi] = alt

    return box1, l_fiber, l_cross, b_exact, c_coexact, delta_1, u, a0, l0


def fit_to_base(op12, a0, l0):
    design = np.column_stack(
        [
            np.eye(12).reshape(-1),
            l0.reshape(-1),
            a0.reshape(-1),
        ]
    )
    coeffs, _, _, _ = lstsq(design, op12.reshape(-1), rcond=None)
    rec = coeffs[0] * np.eye(12) + coeffs[1] * l0 + coeffs[2] * a0
    return coeffs, norm(op12 - rec)


def block_data(op12, l0):
    vals, vecs = eigh(l0)
    targets = [0.0, 5.0 - math.sqrt(5.0), 6.0, 5.0 + math.sqrt(5.0)]
    projs = []
    for target in targets:
        mask = np.abs(vals - target) < TOL
        projs.append(vecs[:, mask] @ vecs[:, mask].T)

    off = 0.0
    spectra = []
    for p in projs:
        blk = p @ op12 @ p
        evals = np.linalg.eigvalsh((blk + blk.T) / 2.0)
        evals = evals[np.abs(evals) > 1e-10]
        spectra.append(np.round(evals, 8))
    for i, p in enumerate(projs):
        for j, q in enumerate(projs):
            if i != j:
                off = max(off, norm(p @ op12 @ q))
    return spectra, off, norm(op12 @ l0 - l0 @ op12)


def main():
    print("=" * 72)
    print("exp619: no-go for naive effective gauge operators on amplitudes")
    print("=" * 72)

    box1, l_fiber, l_cross, b_exact, c_coexact, delta_1, u, a0, l0 = build_geometry()

    print("\n[1] First-order compression to amplitudes")
    for name, op in (
        ("L_fiber", l_fiber),
        ("L_cross", l_cross),
        ("Box_1", box1),
        ("B", b_exact),
        ("C", c_coexact),
        ("Delta_1", delta_1),
    ):
        comp = u.T @ op @ u
        coeffs, err = fit_to_base(comp, a0, l0)
        evals = np.round(np.sort(eigh(comp)[0]), 8)
        print(f"  {name}: coeffs[I,L0,A0] = {np.round(coeffs, 8)}, fit err = {err:.3e}")
        print(f"         spectrum = {evals}")

    print("\n[2] Naive second-order Schur complement")
    vals, vecs = eigh(box1)
    inv = np.zeros_like(vals)
    inv[np.abs(vals) > 1e-8] = 1.0 / vals[np.abs(vals) > 1e-8]
    box_pinv = vecs @ np.diag(inv) @ vecs.T

    for name, v in (("B", b_exact), ("C", c_coexact), ("Delta_1", delta_1)):
        heff = -u.T @ v @ box_pinv @ v @ u
        coeffs, err = fit_to_base(heff, a0, l0)
        spectra, off, comm = block_data(heff, l0)
        print(f"  Heff[{name}]: coeffs[I,L0,A0] = {np.round(coeffs, 8)}, fit err = {err:.3e}")
        print(f"              ||[Heff,L0]|| = {comm:.3e}, max off-block = {off:.3e}")
        print(f"              block spectra = {spectra}")

    print("\nInterpretation")
    print("  Direct compression of natural edge operators to the amplitude space")
    print("  is completely scalar: the vertical gauge kernel carries no nontrivial")
    print("  Hopf-base dynamics at first order. A naive Schur-complement treatment")
    print("  does generate nontrivial amplitude operators, but they are neither")
    print("  simple functions of the base Laplacian nor cleanly diagonal in the")
    print("  1+3+5+3' decomposition. This is evidence against a simplistic KK")
    print("  reduction by projection alone.")


if __name__ == "__main__":
    main()

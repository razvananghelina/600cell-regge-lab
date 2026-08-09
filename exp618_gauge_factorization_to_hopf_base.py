"""
exp618: Exact factorization of the 12 gauge modes into
        (fixed fiber polarization) x (Hopf-base amplitudes).

Goal:
  Couple the edge-gauge result to the Hopf-base continuum map.

Key observation:
  The 12 gauge modes in ker(Box_1) are supported on fiber edges only.  The
  natural next question is whether they are arbitrary edge patterns on each
  decagon, or whether they factor through a single canonical fiber
  polarization.  If they do, then the gauge skeleton is canonically isomorphic
  to functions on the 12-point Hopf base, and exp616 becomes directly relevant.

Checks:
  1. Extract the 12 gauge modes from ker(Box_1) = rho_0 + 2 rho_5.
  2. Build the canonical "vertical" basis U consisting of the alternating
     mode on each Hopf fiber C10.
  3. Show that the gauge subspace equals span(U) exactly.
  4. Compute the amplitude matrix A = U^T gauge and verify it is orthogonal.
  5. Show that the induced A5 action on amplitudes matches the 12-point
     permutation representation on fibers.
  6. Diagonalize the Hopf-base Laplacian on that amplitude space, recovering
     1 + 3 + 5 + 3 with the same icosahedral spectrum as exp616.

Interpretation:
  The gauge modes are not 12 arbitrary edge patterns. They are a fixed
  vertical fiber polarization tensored with a 12-dimensional amplitude space
  on the Hopf base. This is the clean discrete analogue of a Kaluza-Klein
  vertical 1-form with base-dependent coefficients.
"""

from collections import defaultdict
import math
import sys

import numpy as np
from numpy.linalg import eigh, norm, svd

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


def build_edge_data(verts, adj):
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
    return edges, edge_to_idx, adj_list


def build_box1(edges, fibers, edge_to_idx, nv):
    ne = len(edges)
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

    a_line = np.zeros((ne, ne))
    for v in range(nv):
        inc = vte[v]
        for a in range(len(inc)):
            for b in range(a + 1, len(inc)):
                ea = inc[a]
                eb = inc[b]
                a_line[ea, eb] = 1.0
                a_line[eb, ea] = 1.0

    a_fiber = np.zeros((ne, ne))
    for fiber_edges in fiber_edge_lists:
        for k in range(10):
            ea = fiber_edges[k]
            eb = fiber_edges[(k + 1) % 10]
            a_fiber[ea, eb] = 1.0
            a_fiber[eb, ea] = 1.0

    l_fiber = np.diag(np.sum(a_fiber, axis=1)) - a_fiber
    l_cross = np.diag(np.sum(a_line - a_fiber, axis=1)) - (a_line - a_fiber)
    box1 = l_cross - A1 * l_fiber
    return box1, fiber_edge_lists, vertex_to_fiber


def build_gauge_basis(box1, verts, edges, edge_to_idx):
    evals, evecs = eigh(box1)
    ker = evecs[:, np.abs(evals) < 1e-7]

    all_edge_perms = []
    for g_idx in range(len(verts)):
        vp = np.array([find_idx(qmul(verts[g_idx], verts[i]), verts) for i in range(len(verts))])
        ep = np.array([edge_to_idx[(min(vp[i], vp[j]), max(vp[i], vp[j]))] for i, j in edges])
        all_edge_perms.append(ep)

    p0 = np.zeros((13, 13))
    for ep in all_edge_perms:
        for ki in range(13):
            gv = ker[:, ki][ep]
            for kj in range(13):
                p0[ki, kj] += np.dot(ker[:, kj], gv)
    p0 /= float(len(all_edge_perms))

    p5 = np.eye(13) - p0
    p5e, p5v = eigh(p5)
    gauge = ker @ p5v[:, p5e > 0.5]
    return gauge, all_edge_perms


def build_canonical_vertical_basis(ne, fiber_edge_lists):
    alt = np.array([(-1.0) ** k for k in range(10)], dtype=float)
    alt /= norm(alt)
    u = np.zeros((ne, 12))
    for fi, fiber_edges in enumerate(fiber_edge_lists):
        u[fiber_edges, fi] = alt
    return u


def build_base_graph(adj, vertex_to_fiber):
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


def main():
    print("=" * 72)
    print("exp618: gauge-factorization to the Hopf base")
    print("=" * 72)

    verts, adj, _ = build_600cell()
    nv = len(verts)
    edges, edge_to_idx, _ = build_edge_data(verts, adj)
    fibers = find_hopf_fibration(verts)
    box1, fiber_edge_lists, vertex_to_fiber = build_box1(edges, fibers, edge_to_idx, nv)
    gauge, all_edge_perms = build_gauge_basis(box1, verts, edges, edge_to_idx)
    vertical = build_canonical_vertical_basis(len(edges), fiber_edge_lists)

    print("\n[1] Exact factorization of the gauge subspace")
    singular_values = svd(vertical.T @ gauge, compute_uv=False)
    print(f"  Singular values between canonical vertical basis and gauge space:")
    print(f"  {np.round(singular_values, 10)}")
    print(f"  Subspace match: {np.allclose(singular_values, 1.0, atol=1e-10)}")

    amplitude = vertical.T @ gauge
    print("\n[2] Amplitude matrix on the 12 fibers")
    print(f"  ||A^T A - I|| = {norm(amplitude.T @ amplitude - np.eye(12)):.3e}")
    print(f"  ||A A^T - I|| = {norm(amplitude @ amplitude.T - np.eye(12)):.3e}")

    print("\n[3] Microscopic fiber pattern")
    pattern_cov = np.zeros((10, 10))
    for fiber_edges in fiber_edge_lists:
        block = gauge[fiber_edges, :]
        pattern_cov += block @ block.T
    we, ve = eigh(pattern_cov)
    top_pattern = ve[:, -1]
    top_pattern /= norm(top_pattern)
    alt = np.array([(-1.0) ** k for k in range(10)], dtype=float)
    alt /= norm(alt)
    overlap = abs(np.dot(top_pattern, alt))
    print(f"  Top fiber-pattern eigenvalues: {np.round(we, 8)}")
    print(f"  Overlap with alternating C10 mode: {overlap:.12f}")
    print(f"  Top pattern: {np.round(top_pattern, 4)}")

    print("\n[4] A5 action on amplitudes = permutation action on fibers")
    amp_perms = []
    fiber_perms = []
    for ep in all_edge_perms:
        amp_perm = vertical.T @ vertical[ep, :]
        amp_perms.append(amp_perm)

    for g_idx in range(len(verts)):
        vp = np.array([find_idx(qmul(verts[g_idx], verts[i]), verts) for i in range(len(verts))])
        fiber_perm = np.zeros((12, 12))
        for fi in range(12):
            image = vertex_to_fiber[vp[fibers[fi][0]]]
            fiber_perm[image, fi] = 1.0
        fiber_perms.append(fiber_perm)

    diffs = [norm(np.abs(amp_perms[g]) - fiber_perms[g].T) for g in range(len(verts))]
    print(f"  max_g || |P_amp(g)| - P_fiber(g)^T || = {max(diffs):.3e}")

    minus_one = None
    for g_idx in range(len(verts)):
        if np.allclose(verts[g_idx], [-1.0, 0.0, 0.0, 0.0], atol=1e-6):
            minus_one = g_idx
            break
    if minus_one is not None:
        print(
            "  Central element -1 acts on amplitudes as: "
            f"{np.unique(np.round(np.diag(amp_perms[minus_one]), 8))}"
        )

    print("\n[5] Base Laplacian on the amplitude space")
    _, laplacian = build_base_graph(adj, vertex_to_fiber)
    evals = np.sort(eigh(laplacian)[0])
    print(f"  Base Laplacian spectrum: {np.round(evals, 8)}")
    print("  Expected decomposition: 1 + 3 + 5 + 3")
    print(
        "  Multiplicities: "
        f"{[int(np.sum(np.abs(evals - x) < 1e-8)) for x in [0.0, 5.0-math.sqrt(5.0), 6.0, 5.0+math.sqrt(5.0)]]}"
    )

    print("\nInterpretation")
    print("  The 12 gauge modes are exactly a fixed alternating polarization on")
    print("  each Hopf fiber tensored with a 12-dimensional amplitude space on")
    print("  the Hopf base. This identifies the gauge skeleton canonically with")
    print("  functions on the 12-point icosahedral base, i.e. with the same")
    print("  amplitude space used in the scalar continuum map of exp616.")
    print("  The induced group action is a signed permutation lift of the fiber")
    print("  permutation action; the sign comes from the alternating fiber mode,")
    print("  with the central element -1 acting as -I on amplitudes.")
    print("  The remaining open problem is then no longer the existence of a")
    print("  base-space amplitude sector, but how to lift it from this vertical")
    print("  polarization picture to a full nonabelian continuum connection.")


if __name__ == "__main__":
    main()

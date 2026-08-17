"""
exp620: No-go for a canonical local signed transport on the Hopf base.

Motivation:
  exp618 shows that the gauge amplitudes carry a signed-permutation lift of the
  Hopf-base fiber action. The natural next question is whether this lift
  defines a canonical local Z2 connection on the icosahedral base graph.

Test:
  For each adjacent pair of fibers i -> j, collect all group elements of 2I
  whose induced signed permutation sends the i-th amplitude basis vector to
  the j-th one. Then ask:

    Does there exist a canonical sign choice already from the group data?

  In particular, even if one chooses the "minimal rotation" representatives
  (maximal scalar part w), is the sign fixed?

Result:
  No. For every tested neighboring pair, both signs occur among the minimal
  rotation lifts. The global signed lift exists, but it does not descend to a
  unique local edge sign on the Hopf base from group action alone.

Interpretation:
  The missing gauge connection is not hidden in a naive local sign assignment
  on the base graph. Additional geometric/dynamical structure is required.
"""

from collections import defaultdict
import math
import sys

import numpy as np
from numpy.linalg import norm

sys.path.insert(0, ".")
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0


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


def build_signed_amplitude_matrices():
    verts, adj, _ = build_600cell()
    nv = len(verts)

    edges = []
    edge_to_idx = {}
    for i in range(nv):
        for j in range(i + 1, nv):
            if adj[i, j] > 0.5:
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

    alt = np.array([(-1.0) ** k for k in range(10)], dtype=float)
    alt /= norm(alt)
    u = np.zeros((len(edges), 12))
    for fi, fiber_edges in enumerate(fiber_edge_lists):
        u[fiber_edges, fi] = alt

    signed = []
    fiber_adj = np.zeros((12, 12))
    for g_idx in range(nv):
        vp = np.array([find_idx(qmul(verts[g_idx], verts[i]), verts) for i in range(nv)])
        ep = np.array([edge_to_idx[(min(vp[i], vp[j]), max(vp[i], vp[j]))] for i, j in edges])
        signed.append(u.T @ u[ep, :])

    for i in range(nv):
        for j in range(i + 1, nv):
            if adj[i, j] > 0.5:
                fi = vertex_to_fiber[i]
                fj = vertex_to_fiber[j]
                if fi != fj:
                    fiber_adj[fi, fj] = 1.0
                    fiber_adj[fj, fi] = 1.0

    return verts, signed, fiber_adj


def main():
    print("=" * 72)
    print("exp620: no-go for canonical local signed transport")
    print("=" * 72)

    verts, signed, fiber_adj = build_signed_amplitude_matrices()

    print("\n[1] Neighbor-pair sign ambiguity")
    examples = []
    for src in range(12):
        for dst in range(12):
            if fiber_adj[src, dst] < 0.5:
                continue

            candidates = []
            for g_idx, mat in enumerate(signed):
                nz = np.where(np.abs(mat[:, src]) > 1e-8)[0]
                if len(nz) == 1 and nz[0] == dst:
                    candidates.append((g_idx, verts[g_idx, 0], mat[dst, src]))

            signs = sorted(set(int(np.sign(c[2])) for c in candidates))
            max_w = max(c[1] for c in candidates)
            minimal = [c for c in candidates if abs(c[1] - max_w) < 1e-10]
            min_signs = sorted(set(int(np.sign(c[2])) for c in minimal))

            examples.append((src, dst, len(candidates), signs, max_w, len(minimal), min_signs))

    for src, dst, count, signs, max_w, nmin, min_signs in examples[:10]:
        print(
            f"  {src}->{dst}: {count} lifts, signs={signs}, "
            f"max w={max_w:.10f}, minimal reps={nmin}, minimal signs={min_signs}"
        )

    both_signs = sum(1 for _, _, _, signs, _, _, _ in examples if signs == [-1, 1])
    both_signs_min = sum(1 for _, _, _, _, _, _, signs in examples if signs == [-1, 1])
    print(f"\n  Neighbor pairs with both signs among all lifts: {both_signs}/{len(examples)}")
    print(f"  Neighbor pairs with both signs already among minimal lifts: {both_signs_min}/{len(examples)}")

    print("\n[2] Central element")
    minus_one = None
    for g_idx, v in enumerate(verts):
        if np.allclose(v, [-1.0, 0.0, 0.0, 0.0], atol=1e-6):
            minus_one = g_idx
            break
    if minus_one is not None:
        diag = np.unique(np.round(np.diag(signed[minus_one]), 8))
        print(f"  -1 acts as diagonal values {diag}")

    print("\nInterpretation")
    print("  The signed lift on amplitudes is real and exact, but it is global.")
    print("  It does not descend to a unique local Z2 transport on Hopf-base")
    print("  edges from the group action alone: even the minimal-rotation lifts")
    print("  already come in both signs. Any local gauge connection therefore")
    print("  needs extra structure beyond the bare signed permutation lift.")


if __name__ == "__main__":
    main()

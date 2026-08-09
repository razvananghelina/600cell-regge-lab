"""
exp610: IR projector selection of the neutral scale.

Idea:
  If n=25 is the neutral vacuum scale, it should be selected by a genuine
  spectral functional, not only by a multiplicity observation.

  The cleanest candidate is the IR projector of the vertex adjacency:
      P_IR(t) = exp(-t A^2)
  whose t -> infinity limit is the spectral projector onto ker(A).

Since dim ker(A) = 25 exactly, the nontrivial IR sector selected by A is 25.
The experiment also checks how Box refines this neutral sector.
"""

from collections import defaultdict
import math
import sys

import numpy as np
from numpy.linalg import eigh

sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1.0 + math.sqrt(5.0)) / 2.0
A1 = 5.0


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
    idx = np.argmax(dots)
    return idx if dots[idx] > 1.0 - tol else -1


def find_fibration(verts):
    n = len(verts)
    target_w = PHI / 2.0
    for i in range(n):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            p = g.copy()
            ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1, 0, 0, 0], atol=1e-6):
                    ok = False
                    break
                if k == 10 and not np.allclose(p, [1, 0, 0, 0], atol=1e-6):
                    ok = False
            if not ok:
                continue
            used = set()
            fibers = []
            subgroup = []
            pp = np.array([1.0, 0, 0, 0])
            for _ in range(10):
                subgroup.append(find_idx(pp, verts))
                pp = qmul(pp, g)
            for s in range(n):
                if s in used:
                    continue
                fib = []
                for si in subgroup:
                    idx = find_idx(qmul(verts[s], verts[si]), verts)
                    if idx >= 0 and idx not in used:
                        fib.append(idx)
                        used.add(idx)
                if len(fib) == 10:
                    fibers.append(fib)
            if len(fibers) == 12:
                return fibers
    raise RuntimeError("Could not find Hopf fibration")


def spectral_blocks(adj):
    evals, evecs = eigh(adj)
    groups = defaultdict(list)
    for i, val in enumerate(np.round(evals, 10)):
        groups[float(val)].append(i)
    return evals, evecs, [(lam, evecs[:, inds]) for lam, inds in sorted(groups.items())]


def main():
    print("=" * 72)
    print("EXP610: IR PROJECTOR SELECTION OF THE NEUTRAL SCALE")
    print("=" * 72)

    verts, adj, lap = build_600cell()
    _ = lap
    evals_a, evecs_a, blocks = spectral_blocks(adj)

    print("\n[1] Adjacency IR selector")
    zero_mask = np.abs(evals_a) < 1e-8
    n_zero = int(np.sum(zero_mask))
    gap = float(np.min(np.abs(evals_a[~zero_mask])))
    print(f"  dim ker(A) = {n_zero}")
    print(f"  Spectral gap away from zero = {gap:.6f}")
    print(f"  Thus any threshold |lambda| < {gap:.6f} selects exactly {n_zero} modes.")

    print("\n[2] Heat-kernel projector P_IR(t) = exp(-t A^2)")
    print(f"  {'t':>8s} {'Tr exp(-tA^2)':>18s} {'excess over 25':>18s}")
    for t in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        tr = float(np.sum(np.exp(-t * evals_a**2)))
        print(f"  {t:8.3f} {tr:18.10f} {tr - 25.0:18.10e}")
    print("  Limit t -> infinity: Tr exp(-t A^2) -> dim ker(A) = 25 exactly.")

    print("\n[3] Blockwise IR contributions")
    print(f"  {'lambda_A':>12s} {'mult':>6s} {'contrib@t=1':>16s}")
    for lam, block in blocks:
        mult = block.shape[1]
        contrib = mult * math.exp(-(lam ** 2))
        print(f"  {lam:12.6f} {mult:6d} {contrib:16.10f}")
    print("  The only surviving nontrivial IR plateau is the neutral 25-dimensional block.")

    print("\n[4] Build Box and restrict to ker(A)")
    fibers = find_fibration(verts)
    a_fiber = np.zeros_like(adj)
    for fib in fibers:
        for i in fib:
            for j in fib:
                if i != j and adj[i, j] > 0.5:
                    a_fiber[i, j] = 1.0
    a_cross = adj - a_fiber
    l_fiber = np.diag(np.sum(a_fiber, axis=1)) - a_fiber
    l_cross = np.diag(np.sum(a_cross, axis=1)) - a_cross
    box = l_cross - A1 * l_fiber

    ker_basis = evecs_a[:, zero_mask]
    rf = ker_basis.T @ l_fiber @ ker_basis
    rb = ker_basis.T @ box @ ker_basis
    evals_f = np.round(np.linalg.eigvalsh(rf), 10)
    evals_b = np.round(np.linalg.eigvalsh(rb), 10)

    print("  L_fiber spectrum on ker(A):")
    print(f"    {evals_f}")
    print("  Box spectrum on ker(A):")
    print(f"    {evals_b}")

    mult_f = defaultdict(int)
    mult_b = defaultdict(int)
    for x in evals_f:
        mult_f[float(x)] += 1
    for x in evals_b:
        mult_b[float(x)] += 1

    print("\n[5] Exact form inside the neutral sector")
    print("  L_fiber|ker(A) = 0^(5)  + (3-phi)^(10) + (phi+2)^(10)")
    print("  Because Box = 12*I - 6*L_fiber on vertices:")
    print("  Box|ker(A) = 12^(5) + (6/phi)^(10) + (-6phi)^(10)")
    print(f"  Observed multiplicities L_fiber: {dict(mult_f)}")
    print(f"  Observed multiplicities Box: {dict(mult_b)}")

    print("\n[6] Verdict")
    print("  The adjacency operator has an exact IR projector whose selected rank is 25.")
    print("  This makes 25 a genuine spectral output, not just an RG-fitted exponent.")
    print("  Box does not change that rank; it refines the neutral sector into")
    print("  a 5 + 10 + 10 internal structure.")
    print("  What remains open is the final physical identification")
    print("  ker(A) IR rank 25  ->  electroweak rung m_e * phi^25.")

    print("\n" + "=" * 72)
    print("EXP610 COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()

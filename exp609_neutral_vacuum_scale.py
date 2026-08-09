"""
exp609: Neutral vacuum scale from adjacency blindness.

Question:
  Why should the electroweak rung be phi^25 rather than phi^24 or phi^26?

Hypothesis:
  n = 25 is selected because the adjacency operator A has a unique neutral
  isotypic block: the 5-dimensional irrep with vanishing generator character.
  In the regular 2I action on the 120 vertices, that block appears with
  multiplicity 5, hence total dimension 5^2 = 25.

This does not yet derive the physical lattice spacing. It tests whether 25
is the unique internally distinguished "adjacency-blind" scale candidate.
"""

from collections import defaultdict
import math
import sys

import numpy as np
from numpy.linalg import eigh, norm

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
    idx = np.argmax(dots)
    return idx if dots[idx] > 1.0 - tol else -1


def spectral_blocks(adj):
    evals, evecs = eigh(adj)
    groups = defaultdict(list)
    for i, val in enumerate(np.round(evals, 10)):
        groups[float(val)].append(i)
    out = []
    for val, inds in sorted(groups.items()):
        block = evecs[:, inds]
        out.append((val, block))
    return out


def projector(block):
    return block @ block.T


def main():
    print("=" * 72)
    print("EXP609: NEUTRAL VACUUM SCALE FROM ADJACENCY BLINDNESS")
    print("=" * 72)

    verts, adj, lap = build_600cell()
    n = len(verts)
    deg = int(round(np.sum(adj[0])))
    identity = find_idx(np.array([1.0, 0, 0, 0]), verts)
    if identity < 0:
        raise RuntimeError("Identity quaternion not found")

    generators = [j for j in range(n) if adj[identity, j] > 0.5]
    print("\n[1] Generator class")
    print(f"  Identity index: {identity}")
    print(f"  Number of generators: {len(generators)}")
    traces = [2.0 * verts[g, 0] for g in generators]
    print(f"  Generator traces: min={min(traces):.6f}, max={max(traces):.6f}")
    print(f"  Common trace = phi: {all(abs(t - PHI) < 1e-6 for t in traces)}")

    print("\n[2] Adjacency spectral blocks")
    blocks = spectral_blocks(adj)
    print(f"  Distinct adjacency eigenvalues: {len(blocks)}")
    print(f"  {'lambda_A':>12s} {'mult':>6s} {'d=sqrt(mult)':>12s}")
    for lam, block in blocks:
        mult = block.shape[1]
        d = round(math.sqrt(mult))
        print(f"  {lam:12.6f} {mult:6d} {d:12d}")

    print("\n[3] Left regular action and characters")
    perms = []
    for g_idx in range(n):
        vp = np.array([find_idx(qmul(verts[g_idx], verts[i]), verts) for i in range(n)])
        perms.append(vp)

    g0 = generators[0]
    print(f"  Using generator index {g0} as class representative.")
    print(f"  Character relation to test: lambda_A = |S| * chi(g) / dim(rho)")
    print(f"  Here |S| = {len(generators)}.\n")

    neutral_block = None
    neutral_dim = None
    neutral_mult = None
    all_rows = []
    for lam, block in blocks:
        mult = block.shape[1]
        d = int(round(math.sqrt(mult)))
        p = projector(block)
        chi_isotypic = np.trace(p[:, perms[g0]])
        chi_irrep = chi_isotypic / d
        lam_from_char = len(generators) * chi_irrep / d
        blindness = norm(adj @ block) / math.sqrt(mult)
        all_rows.append((lam, mult, d, chi_irrep, lam_from_char, blindness))
        if abs(lam) < 1e-8:
            neutral_block = block
            neutral_dim = d
            neutral_mult = mult

    print(f"  {'lambda_A':>12s} {'mult':>6s} {'dim':>6s} {'chi(g)':>10s} {'12*chi/d':>10s} {'||A||_blk':>10s}")
    for lam, mult, d, chi_irrep, lam_from_char, blindness in all_rows:
        print(
            f"  {lam:12.6f} {mult:6d} {d:6d} {chi_irrep:10.6f} "
            f"{lam_from_char:10.6f} {blindness:10.3e}"
        )

    print("\n[4] Neutral block")
    neutral_blocks = [row for row in all_rows if abs(row[0]) < 1e-8]
    print(f"  Number of neutral adjacency blocks: {len(neutral_blocks)}")
    if len(neutral_blocks) != 1:
        raise RuntimeError("Expected a unique neutral adjacency block")

    lam, mult, d, chi_irrep, lam_from_char, blindness = neutral_blocks[0]
    _ = neutral_block, neutral_dim, neutral_mult, lam, chi_irrep, lam_from_char, blindness
    print(f"  Neutral block multiplicity: {mult}")
    print(f"  Underlying irrep dimension: {d}")
    print(f"  Generator character chi(g): {chi_irrep:.6f}")
    print(f"  A annihilates the block: ||A||_blk = {blindness:.3e}")
    print(f"  Candidate neutral scale = mult = d^2 = {d*d}")

    print("\n[5] Why 25 and not 24 or 26?")
    nearby = [24, 25, 26]
    multiplicities = sorted(block.shape[1] for _, block in blocks)
    square_mults = sorted({int(round(math.sqrt(m))) ** 2 for m in multiplicities})
    for m in nearby:
        is_mult = m in multiplicities
        is_square_mult = m in square_mults
        is_neutral = m == d * d
        print(
            f"  n={m}: spectral multiplicity={is_mult}, square-isotypic={is_square_mult}, "
            f"neutral-block={is_neutral}"
        )

    print("\n[6] Interpretation")
    print("  The 25-dimensional adjacency kernel is unique.")
    print("  It comes from the unique 5D irrep with vanishing generator character.")
    print("  In the regular 2I action, isotypic multiplicity equals irrep dimension,")
    print("  so the neutral block has total size 5 x 5 = 25.")
    print("  This makes n=25 the unique adjacency-blind scale candidate.")
    print("  What is still open is the physical step from this neutral spectral block")
    print("  to the electroweak lattice spacing m_e * phi^25.")

    print("\n" + "=" * 72)
    print("EXP609 COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()

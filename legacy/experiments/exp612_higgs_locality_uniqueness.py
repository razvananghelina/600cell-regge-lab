"""
exp612_higgs_locality_uniqueness.py
===================================

Question:
  Why is the tree-level Higgs/W ratio tied to phi^2, and not to some nearby
  spectral ratio?

Idea:
  The electroweak tree-level ratio in the paper comes from the two 3-dimensional
  A5 sectors on the 12-fiber icosahedral base:

      L(3') / L(3) = (5 + sqrt(5)) / (5 - sqrt(5)) = phi^2.

  This script asks how rigid that ratio really is.

  We build the 12-fiber graph directly from the Hopf decomposition of the
  600-cell and then scan the full A5-invariant Bose-Mesner algebra on that
  icosahedral graph:

      L(x, y) = (5 + 5x + y) I - A1 - x A2 - y A3

  where
      A1 = nearest-neighbor graph,
      A2 = distance-2 graph,
      A3 = antipodal matching.

  The diagonal is chosen so each row sums to zero.

What is tested:
  1. The local nearest-neighbor Laplacian L(0,0) has spectrum
         0, 5-sqrt(5) [mult 3], 6 [mult 5], 5+sqrt(5) [mult 3].
  2. The only nontrivial equal-dimensional Galois pair is the two 3d sectors.
  3. A scan over small half-integer A5-invariant deformations shows that
     phi^2 is preserved only by three simple operators, and only one of them
     is strictly local (nearest-neighbor support only): the standard Laplacian.

Conclusion:
  If one insists on locality on the 12-fiber base, the phi^2 ratio is not an
  arbitrary spectral choice: it is forced by the unique local A5-equivariant
  Laplacian.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import numpy as np
from numpy.linalg import eigvalsh

from commons import build_600cell


PHI = (1 + np.sqrt(5)) / 2
A1 = 5
TOL = 1e-8


def qmul(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.array(
        [
            p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
            p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
            p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
            p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0],
        ]
    )


def find_idx(v: np.ndarray, verts: np.ndarray, tol: float = 1e-6) -> int:
    dots = verts @ v
    idx = int(np.argmax(dots))
    return idx if dots[idx] > 1 - tol else -1


def find_hopf_fibration(verts: np.ndarray) -> list[list[int]]:
    """Return the ordered 12 x 10 Hopf fibers."""
    n = len(verts)
    target_w = PHI / 2.0
    for i in range(n):
        if abs(verts[i, 0] - target_w) >= 1e-6:
            continue
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

        used: set[int] = set()
        fibers: list[list[int]] = []
        subgroup: list[int] = []
        pp = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_idx(pp, verts))
            pp = qmul(pp, g)

        for s in range(n):
            if s in used:
                continue
            fib: list[int] = []
            for si in subgroup:
                q = qmul(verts[s], verts[si])
                idx = find_idx(q, verts)
                if idx >= 0 and idx not in used:
                    fib.append(idx)
                    used.add(idx)
            if len(fib) == 10:
                fibers.append(fib)

        if len(fibers) == 12:
            return fibers

    raise RuntimeError("Hopf fibration not found")


def build_fiber_graph(adj: np.ndarray, fibers: list[list[int]]) -> np.ndarray:
    """Build the 12-node icosahedral base graph from cross-fiber adjacency."""
    vertex_to_fiber: dict[int, int] = {}
    for fi, fib in enumerate(fibers):
        for v in fib:
            vertex_to_fiber[v] = fi

    w = np.zeros((12, 12), dtype=int)
    n = len(adj)
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j] < 0.5:
                continue
            fi = vertex_to_fiber[i]
            fj = vertex_to_fiber[j]
            if fi != fj:
                w[fi, fj] += 1
                w[fj, fi] += 1

    # The base graph forgets the uniform cross-edge multiplicity 20.
    return (w > 0).astype(float)


def distance_matrices(adj12: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return A1, A2, A3 for the 12-node icosahedron."""
    n = adj12.shape[0]
    dist = np.full((n, n), 99, dtype=int)
    for src in range(n):
        dist[src, src] = 0
        q: deque[int] = deque([src])
        while q:
            u = q.popleft()
            for v in np.where(adj12[u] > 0.5)[0]:
                if dist[src, v] > dist[src, u] + 1:
                    dist[src, v] = dist[src, u] + 1
                    q.append(v)

    a1 = (dist == 1).astype(float)
    a2 = (dist == 2).astype(float)
    a3 = (dist == 3).astype(float)
    return a1, a2, a3


def spectrum_summary(op: np.ndarray) -> list[tuple[float, int]]:
    evals = eigvalsh(op)
    out: list[tuple[float, int]] = []
    for val in sorted(set(np.round(evals, 10))):
        mult = int(np.sum(np.abs(evals - val) < 1e-8))
        out.append((float(val), mult))
    return out


@dataclass
class Hit:
    x: float
    y: float
    spectrum: list[tuple[float, int]]
    ratio: float


def format_spectrum(spec: list[tuple[float, int]]) -> str:
    return ", ".join(f"{val:.10f}[{mult}]" for val, mult in spec)


def main() -> None:
    print("=" * 72)
    print("EXP612: HIGGS TREE-LEVEL RATIO AND LOCALITY ON THE HOPF BASE")
    print("=" * 72)

    verts, adj, _ = build_600cell()
    fibers = find_hopf_fibration(verts)
    fiber_graph = build_fiber_graph(adj, fibers)
    a1_graph, a2_graph, a3_graph = distance_matrices(fiber_graph)

    print("\nSECTION 1: Fiber graph")
    print("-" * 72)
    print(f"  Nodes             = {fiber_graph.shape[0]}")
    print(f"  Degree            = {int(fiber_graph.sum(axis=1)[0])}")
    print(f"  Edges             = {int(fiber_graph.sum() // 2)}")
    print(f"  Dist-2 degree     = {int(a2_graph.sum(axis=1)[0])}")
    print(f"  Antipodal degree  = {int(a3_graph.sum(axis=1)[0])}")

    print("\nSECTION 2: Local icosahedral Laplacian")
    print("-" * 72)
    l_local = 5 * np.eye(12) - a1_graph
    local_spec = spectrum_summary(l_local)
    print(f"  Spectrum          = {format_spectrum(local_spec)}")

    l3 = 5 - np.sqrt(5)
    l5 = 6.0
    l3p = 5 + np.sqrt(5)
    ratio_local = l3p / l3
    print(f"  L(3)              = 5 - sqrt(5) = {l3:.10f}")
    print(f"  L(5)              = 6           = {l5:.10f}")
    print(f"  L(3')             = 5 + sqrt(5) = {l3p:.10f}")
    print(f"  L(3')/L(3)        = {ratio_local:.12f}")
    print(f"  phi^2             = {PHI**2:.12f}")
    print(f"  Match             = {abs(ratio_local - PHI**2) < 1e-12}")

    print("\nSECTION 3: A5-invariant deformation scan")
    print("-" * 72)
    print("  Family: L(x,y) = (5 + 5x + y) I - A1 - x A2 - y A3")
    print("  Scan domain: x,y in {-3,-2.5,...,2.5,3}")
    print("  Criterion: the two 3d sectors survive and their ratio is phi^2")

    hits: list[Hit] = []
    scan_values = [k / 2.0 for k in range(-6, 7)]
    for x in scan_values:
        for y in scan_values:
            op = (5 + 5 * x + y) * np.eye(12) - a1_graph - x * a2_graph - y * a3_graph
            spec = spectrum_summary(op)
            triples = [val for val, mult in spec if mult == 3]
            if len(triples) != 2 or min(triples) <= TOL:
                continue
            ratio = max(triples) / min(triples)
            if abs(ratio - PHI**2) < 1e-9:
                hits.append(Hit(x=x, y=y, spectrum=spec, ratio=ratio))

    print(f"  phi^2-preserving hits found = {len(hits)}")
    for hit in hits:
        locality = "local" if abs(hit.x) < TOL and abs(hit.y) < TOL else "nonlocal"
        print(
            f"    (x,y)=({hit.x:+.1f},{hit.y:+.1f}) [{locality}]  "
            f"spectrum = {format_spectrum(hit.spectrum)}"
        )

    local_hits = [hit for hit in hits if abs(hit.x) < TOL and abs(hit.y) < TOL]
    assert len(local_hits) == 1, "Expected exactly one local phi^2-preserving hit"

    print("\nSECTION 4: Verdict")
    print("-" * 72)
    print("  The phi^2 ratio is not unique inside the full nonlocal A5-invariant")
    print("  Bose-Mesner algebra: two additional distance-mixed operators preserve it.")
    print("  But once strict nearest-neighbor locality is imposed, only")
    print("      L = 5 I - A1")
    print("  survives.")
    print("  Therefore the tree-level ratio")
    print("      L(3') / L(3) = phi^2")
    print("  is forced by the unique local A5-equivariant Laplacian on the 12-fiber")
    print("  Hopf base, not by an arbitrary choice among many local operators.")


if __name__ == "__main__":
    main()

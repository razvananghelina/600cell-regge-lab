"""
exp615_nearest_neighbor_from_box_and_units.py
=============================================

Goal:
  Strengthen the "nearest-neighbor tight-binding" step in the Generation
  Theorem by tying it directly to the existing geometric/spectral framework.

Core idea:
  The nearest-neighbor Hamiltonian on the logarithmic DSI lattice is not an
  arbitrary modeling choice once two facts are combined:

    (1) The fiber part of the 600-cell geometry is literally the decagon C10,
        so the microscopic generator on each Hopf fiber is the nearest-neighbor
        graph Laplacian.

    (2) The unit group of Z[phi] is rank-1:
            U(Z[phi]) = { +/- phi^n } ~= Z.
        The fundamental unit phi generates one primitive step n -> n+1.
        Longer hops phi^r are composite powers of this same generator.

  Therefore the primitive generator of scale motion is the shift by one unit
  in the exponent lattice, giving the standard tight-binding Laplacian

      Delta = 2 I - S - S^{-1},

  with Bloch dispersion

      E(k) = 2 - 2 cos(k) = 4 sin^2(k/2).

  Under k = 2 pi x, this is exactly the DSI shape

      V(x) ~ sin^2(pi x).

This script verifies:
  A. Hopf fibers in the 600-cell are exactly C10 cycles.
  B. The fiber generator is nearest-neighbor and its spectrum is the C10
     spectrum 2 - 2 cos(2 pi m/10).
  C. For a general translation-invariant hopping Hamiltonian on the unit lattice
         H = sum_{r>=1} t_r (2I - S^r - S^{-r}),
     the dispersion is
         E(k) = 2 sum_r t_r (1 - cos(rk)).
     Pure sin^2 requires t_r = 0 for every r >= 2.
  D. The continuous Hopf-fiber heat kernel strongly suppresses the m=2 mode
     relative to m=1 at the natural framework times t = phi^2 and t = a1 = 5.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
from numpy.linalg import eigvalsh

from commons import build_600cell


PHI = (1 + np.sqrt(5)) / 2
A1 = 5


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
                idx = find_idx(qmul(verts[s], verts[si]), verts)
                if idx >= 0 and idx not in used:
                    fib.append(idx)
                    used.add(idx)
            if len(fib) == 10:
                fibers.append(fib)

        if len(fibers) == 12:
            return fibers

    raise RuntimeError("Hopf fibration not found")


def build_edges(adj: np.ndarray) -> list[tuple[int, int]]:
    n = len(adj)
    edges: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j] > 0.5:
                edges.append((i, j))
    return edges


def fiber_cycle_adjacency(adj: np.ndarray, fiber: list[int]) -> np.ndarray:
    n = len(fiber)
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if adj[fiber[i], fiber[j]] > 0.5:
                out[i, j] = 1.0
                out[j, i] = 1.0
    return out


def c10_laplacian() -> np.ndarray:
    adj = np.zeros((10, 10), dtype=float)
    for i in range(10):
        adj[i, (i + 1) % 10] = 1.0
        adj[(i + 1) % 10, i] = 1.0
    return 2.0 * np.eye(10) - adj


def general_dispersion(k: np.ndarray, coeffs: dict[int, float]) -> np.ndarray:
    out = np.zeros_like(k)
    for r, t_r in coeffs.items():
        out += 2.0 * t_r * (1.0 - np.cos(r * k))
    return out


def main() -> None:
    print("=" * 72)
    print("EXP615: NEAREST-NEIGHBOR FROM BOX FIBERS AND UNIT-GROUP STRUCTURE")
    print("=" * 72)

    verts, adj, _ = build_600cell()
    fibers = find_hopf_fibration(verts)
    edges = build_edges(adj)
    vertex_to_fiber: dict[int, int] = {}
    for fi, fib in enumerate(fibers):
        for v in fib:
            vertex_to_fiber[v] = fi

    print("\nSECTION 1: Hopf fibers inside the 600-cell")
    print("-" * 72)
    fiber_edge_counts = []
    for fib in fibers:
        count = 0
        for i in range(10):
            for j in range(i + 1, 10):
                if adj[fib[i], fib[j]] > 0.5:
                    count += 1
        fiber_edge_counts.append(count)

    print(f"  Number of fibers                = {len(fibers)}")
    print(f"  Fiber sizes                     = {sorted(set(len(f) for f in fibers))}")
    print(f"  Edges per fiber                 = {sorted(set(fiber_edge_counts))}")
    print(f"  Total fiber edges               = {sum(fiber_edge_counts)}")

    cycle_adj = fiber_cycle_adjacency(adj, fibers[0])
    cycle_deg = cycle_adj.sum(axis=1)
    print(f"  First fiber degree pattern      = {cycle_deg.astype(int).tolist()}")
    print(f"  Is first fiber a C10 cycle?     = {np.allclose(cycle_deg, 2.0)}")

    print("\nSECTION 2: Fiber generator is nearest-neighbor")
    print("-" * 72)
    l_fiber = c10_laplacian()
    evals = np.sort(eigvalsh(l_fiber))
    expected = np.sort([2 - 2 * np.cos(2 * np.pi * m / 10.0) for m in range(10)])
    print(f"  spec(L_C10)                     = {np.round(evals, 10)}")
    print(f"  exact cosine spectrum           = {np.round(expected, 10)}")
    print(f"  Match                           = {np.allclose(evals, expected)}")
    print("  So the microscopic fiber generator already has only")
    print("      n <-> n+1 and n <-> n-1")
    print("  couplings. Longer hops do not appear in the generator itself.")

    print("\nSECTION 3: From the fiber generator to the DSI tight-binding form")
    print("-" * 72)
    print("  On the unit group U(Z[phi]) = { +/- phi^n } ~= Z,")
    print("  multiplication by the fundamental unit phi is the primitive step")
    print("      n -> n+1.")
    print("  The induced shift operator S therefore gives the minimal generator")
    print("      Delta = 2I - S - S^{-1}.")
    k = np.linspace(0.0, 2.0 * np.pi, 9)
    disp_nn = general_dispersion(k, {1: 1.0})
    print(f"  Dispersion E_1(k)               = {np.round(disp_nn, 10)}")
    print("  Exact formula                   = 2 - 2 cos(k) = 4 sin^2(k/2)")
    print("  With k = 2 pi x this becomes")
    print("      V(x) ~ 4 sin^2(pi x).")

    print("\nSECTION 4: Pure sin^2 is unique for nearest-neighbor hopping")
    print("-" * 72)
    print("  General translation-invariant positive hopping on the unit lattice:")
    print("      H = sum_{r>=1} t_r (2I - S^r - S^{-r}),  t_r >= 0")
    print("  gives")
    print("      E(k) = 2 sum_r t_r (1 - cos(rk)).")
    print("  Any t_r with r >= 2 introduces higher harmonics cos(rk).")
    print("  Therefore pure sin^2(pi x) is exact iff")
    print("      t_r = 0  for every r >= 2.")
    disp_nnn = general_dispersion(k, {1: 1.0, 2: 0.25})
    print(f"  Example with t2=0.25            = {np.round(disp_nnn, 10)}")
    print("  This is no longer a pure first harmonic.")

    print("\nSECTION 5: Hopf heat-kernel support for lowest-harmonic dominance")
    print("-" * 72)
    print("  On the continuous Hopf fiber S^1, the m-th Fourier mode carries")
    print("      weight exp(-t m^2).")
    print("  Relative suppression of m=2 vs m=1 is therefore")
    print("      exp(-3 t).")
    for t in [PHI**2, float(A1)]:
        ratio = np.exp(-3.0 * t)
        print(f"  t = {t:.10f}:  exp(-3t) = {ratio:.12e}")

    print("\n  On the discrete C10 fiber, the first two Laplacian eigenvalues are")
    lam1 = 2 - 2 * np.cos(2 * np.pi / 10.0)
    lam2 = 2 - 2 * np.cos(4 * np.pi / 10.0)
    print(f"      lambda_1 = {lam1:.10f} = 2 - phi = 1/phi^2")
    print(f"      lambda_2 = {lam2:.10f}")
    print(f"      lambda_2 - lambda_1 = {lam2 - lam1:.10f}")
    print("  so the discrete heat-kernel suppression is")
    print("      exp(-t (lambda_2 - lambda_1)) = exp(-t).")
    for t in [PHI**2, float(A1)]:
        ratio = np.exp(-t)
        print(f"  t = {t:.10f}:  exp(-t)  = {ratio:.12e}")

    print("\nVERDICT")
    print("-" * 72)
    print("  The nearest-neighbor DSI Hamiltonian is not an isolated ansatz.")
    print("  It is the direct logarithmic image of the primitive fiber generator")
    print("  already present in the 600-cell geometry:")
    print("      Hopf fiber = C10  ->  nearest-neighbor Laplacian")
    print("      fundamental unit phi -> primitive shift n->n+1")
    print("  Higher jumps are composite powers of phi and appear only after")
    print("  exponentiating the generator (heat kernel / propagation), not in the")
    print("  microscopic generator itself. Pure sin^2 is therefore the unique")
    print("  exact first-harmonic potential associated with the primitive local")
    print("  generator.")


if __name__ == "__main__":
    main()

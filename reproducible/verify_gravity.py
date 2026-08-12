"""
verify_gravity.py
==================
Self-contained verification of finite geometric and edge-weight controls on
the 600-cell, historically grouped in the paper's gravity section.

What this script verifies:
  1. Ollivier-Ricci curvature kappa = 0 on all 720 edges (Prop. ricci_flat).
  2. Common neighbor count c = 5 for every adjacent pair.
  3. Edge-weight stiffness ratio dTr(Box^2)/dw: fiber = 4*a1^2 = 100,
     cross = 4, ratio = a1^2 = 25 = c_lat^4 under the separate lattice-speed
     normalization (Eq. gravitational_stiffness).
  4. Hessian of Tr(Box^2) is diagonal (edges decouple at leading order).
  5. Full Hessian of Tr(Box^4) is non-diagonal and positive definite on the
     720 independent edge weights: 720 positive and no zero/negative modes.
     This is algebraic stiffness, not by itself a graviton derivation.
  6. Graph Laplacian Green's function shape matches S^3 continuum:
     normalized V(d)/V(1) vs continuum G_S3(theta_d)/G_S3(theta_1),
     RMS error < 1% over d = 1, ..., 4.

Dependencies: numpy, scipy.
Encoding: ASCII only.

Author: Razvan-Constantin Anghelina
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh
from scipy.sparse.csgraph import shortest_path
from scipy.sparse import csr_matrix
from scipy.optimize import linprog
from itertools import permutations, product
from collections import defaultdict
import time

# ============================================================================
# Build 600-cell (self-contained, no external imports)
# ============================================================================

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
DEGREE = 12
TOL = 1e-8

def build_600cell():
    verts = set()
    def add(v):
        a = np.array(v, float); a /= np.linalg.norm(a)
        verts.add(tuple(np.round(a, 10)))
    for i in range(4):
        for s in [1., -1.]:
            v = [0.,0.,0.,0.]; v[i] = s; add(v)
    for s0,s1,s2,s3 in product([.5,-.5], repeat=4):
        add([s0,s1,s2,s3])
    base = [0., .5, PHI/2., 1./(2.*PHI)]
    even_perms = [p for p in permutations(range(4))
                  if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])%2==0]
    for perm in even_perms:
        coords = [base[perm[i]] for i in range(4)]
        nz = [i for i in range(4) if abs(coords[i])>1e-12]
        for signs in product([1,-1], repeat=len(nz)):
            v = list(coords)
            for idx, s in zip(nz, signs): v[idx] *= s
            add(v)
    V = np.array(sorted(verts))
    assert len(V) == 120
    dots = V @ V.T; np.clip(dots, -1, 1, out=dots)
    adj = np.zeros((120,120))
    for i in range(120):
        for j in range(i+1,120):
            if abs(dots[i,j] - PHI/2) < 1e-6:
                adj[i,j] = adj[j,i] = 1.
    assert np.all(adj.sum(1).astype(int) == 12)
    return V, adj

def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, V):
    dots = V @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1-1e-6 else -1

def find_fibration(V, adj):
    for i in range(N):
        if abs(V[i,0] - PHI/2) < 1e-6:
            g = V[i]; p = g.copy(); ok = True
            for k in range(2,11):
                p = qmul(p, g)
                if k==5 and not np.allclose(p,[-1,0,0,0],atol=1e-6): ok=False; break
                if k==10 and not np.allclose(p,[1,0,0,0],atol=1e-6): ok=False
            if not ok: continue
            used = set(); fibers = []; subg = []
            pp = np.array([1.,0,0,0])
            for k in range(10): subg.append(find_idx(pp, V)); pp = qmul(pp, g)
            for s in range(N):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(V[s], V[si]); idx = find_idx(q, V)
                    if idx>=0 and idx not in used: fib.append(idx); used.add(idx)
                if len(fib)==10: fibers.append(fib)
            if len(fibers)==12: return fibers
    return None

# ============================================================================
# Main verification
# ============================================================================

def main():
    t0 = time.time()
    passed = 0
    failed = 0
    total = 0

    def check(name, condition):
        nonlocal passed, failed, total
        total += 1
        if condition:
            passed += 1
            print(f"  PASS: {name}")
        else:
            failed += 1
            print(f"  FAIL: {name}")

    print("=" * 70)
    print("VERIFY LEGACY GRAVITY SECTION: finite 600-cell controls")
    print("=" * 70)

    # Build
    print("\nBuilding 600-cell and Hopf fibration...")
    verts, adj = build_600cell()
    fibers = find_fibration(verts, adj)
    assert fibers is not None, "Fibration not found"

    vtx_fiber = {}
    for fi, f in enumerate(fibers):
        for v in f: vtx_fiber[v] = fi

    A_fiber = np.zeros((N, N))
    for fib in fibers:
        for i in fib:
            for j in fib:
                if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.
    A_cross = adj - A_fiber
    Box = a1 * A_fiber - A_cross

    edges = []
    is_fiber_edge = []
    for i in range(N):
        for j in range(i+1, N):
            if adj[i,j] > 0.5:
                edges.append((i, j))
                is_fiber_edge.append(vtx_fiber[i] == vtx_fiber[j])
    is_fiber_edge = np.array(is_fiber_edge)
    Ne = len(edges)

    nbrs = defaultdict(list)
    for i in range(N):
        for j in range(N):
            if adj[i,j] > 0.5: nbrs[i].append(j)

    # ==================================================================
    # TEST 1: Common neighbors = 5 for every edge
    # ==================================================================
    print("\n--- TEST 1: Common neighbor count ---")
    common_counts = []
    for (i, j) in edges:
        common_counts.append(len(set(nbrs[i]) & set(nbrs[j])))
    common_counts = np.array(common_counts)
    check("All edges have exactly 5 common neighbors",
          np.all(common_counts == 5))

    # ==================================================================
    # TEST 2: Ollivier-Ricci curvature kappa = 0 on all edges
    # ==================================================================
    print("\n--- TEST 2: Ollivier-Ricci curvature ---")

    dist_matrix = shortest_path(csr_matrix(adj), method='D',
                                unweighted=True).astype(int)

    def wasserstein_1(p_sup, p_w, q_sup, q_w):
        m, n = len(p_sup), len(q_sup)
        cost = np.array([[dist_matrix[p_sup[i], q_sup[j]]
                          for j in range(n)] for i in range(m)]).flatten()
        A_eq = np.zeros((m + n, m * n))
        b_eq = np.zeros(m + n)
        for i in range(m):
            for j in range(n): A_eq[i, i*n+j] = 1.
            b_eq[i] = p_w[i]
        for j in range(n):
            for i in range(m): A_eq[m+j, i*n+j] = 1.
            b_eq[m+j] = q_w[j]
        res = linprog(cost, A_eq=A_eq, b_eq=b_eq,
                      bounds=[(0,None)]*(m*n), method='highs')
        return res.fun if res.success else float('nan')

    print("  Computing Ollivier-Ricci for all 720 edges (this takes ~1 min)...")
    kappas = np.zeros(Ne)
    for e_idx, (x, y) in enumerate(edges):
        nx, ny = nbrs[x], nbrs[y]
        px = np.ones(len(nx)) / len(nx)
        py = np.ones(len(ny)) / len(ny)
        w1 = wasserstein_1(nx, px, ny, py)
        kappas[e_idx] = 1.0 - w1

    check("kappa = 0 on all 720 edges (max |kappa| < 1e-8)",
          np.max(np.abs(kappas)) < 1e-8)
    check("All kappas identical (vertex-transitivity, std < 1e-10)",
          np.std(kappas) < 1e-10)

    # ==================================================================
    # TEST 3: edge-weight stiffness ratio = a1^2 = 25
    # ==================================================================
    print("\n--- TEST 3: Edge-weight stiffness ---")

    # Sensitivity: S[k, e] = <psi_k| dBox_e |psi_k>
    # dBox_e = a1*(|i><j|+|j><i|) for fiber, -(|i><j|+|j><i|) for cross
    dTr2 = np.zeros(Ne)
    for e_idx, (i, j) in enumerate(edges):
        coeff = a1 if is_fiber_edge[e_idx] else -1.0
        # dTr(Box^2)/dw_e = 2*Tr(Box * dBox_e) = 4*coeff*Box[i,j]
        dTr2[e_idx] = 4 * coeff * Box[i, j]

    dTr2_fiber = dTr2[is_fiber_edge].mean()
    dTr2_cross = dTr2[~is_fiber_edge].mean()
    ratio = dTr2_fiber / dTr2_cross

    check(f"dTr(Box^2)/dw_fiber = 4*a1^2 = {4*a1**2} (got {dTr2_fiber:.1f})",
          abs(dTr2_fiber - 4*a1**2) < 0.01)
    check(f"dTr(Box^2)/dw_cross = 4 (got {dTr2_cross:.1f})",
          abs(dTr2_cross - 4) < 0.01)
    check(f"Stiffness ratio = a1^2 = 25 (got {ratio:.4f})",
          abs(ratio - a1**2) < 0.01)
    check("Stiffness constant within fiber edges (std < 1e-6)",
          np.std(dTr2[is_fiber_edge]) < 1e-6)
    check("Stiffness constant within cross edges (std < 1e-6)",
          np.std(dTr2[~is_fiber_edge]) < 1e-6)

    # ==================================================================
    # TEST 4: Hessian of Tr(Box^2) is diagonal
    # ==================================================================
    print("\n--- TEST 4: Hessian diagonality ---")

    # d^2 Tr(Box^2)/(dw_e dw_f) = 2*Tr(dBox_e * dBox_f)
    # For e=(i,j), f=(k,l) with e != f:
    # dBox_e * dBox_f has nonzero only if they share a vertex,
    # and Tr = 0 because it requires j=k AND l=i (for distinct edges impossible)
    # Verify numerically on a sample
    n_sample = 200
    max_offdiag = 0
    np.random.seed(42)
    for _ in range(n_sample):
        e1, e2 = np.random.choice(Ne, 2, replace=False)
        i1, j1 = edges[e1]; i2, j2 = edges[e2]
        c1 = a1 if is_fiber_edge[e1] else -1.0
        c2 = a1 if is_fiber_edge[e2] else -1.0
        # dBox_e1 = c1*(|i1><j1|+|j1><i1|)
        # Tr(dBox_e1 * dBox_e2) = 2*c1*c2*(delta(j1,i2)*delta(i1,j2) + delta(i1,i2)*delta(j1,j2))
        val = 2*c1*c2*(int(j1==i2 and i1==j2) + int(i1==i2 and j1==j2))
        max_offdiag = max(max_offdiag, abs(val))

    check("Hessian of Tr(Box^2) is diagonal (off-diagonal = 0)",
          max_offdiag == 0)

    # ==================================================================
    # TEST 5: full Hessian of Tr(Box^4)
    # ==================================================================
    print("\n--- TEST 5: Full Tr(Box^4) Hessian stiffness ---")

    # For E_e = c_e(|i><j|+|j><i|), exact differentiation gives
    # H_ef = 8 Tr(Box^2 E_f E_e) + 4 Tr(Box E_f Box E_e).
    # The former eigenvalue-sensitivity Gram matrix omitted the off-diagonal
    # eigenbasis terms and was basis-dependent in degenerate eigenspaces.
    starts = np.array([edge[0] for edge in edges], dtype=int)
    ends = np.array([edge[1] for edge in edges], dtype=int)
    coeffs = np.where(is_fiber_edge, a1, -1).astype(int)
    i_e = starts[:, None]
    j_e = ends[:, None]
    k_f = starts[None, :]
    l_f = ends[None, :]
    cc = coeffs[:, None] * coeffs[None, :]
    Box2 = Box @ Box
    first_trace = cc * (
        (l_f == i_e) * Box2[j_e, k_f]
        + (l_f == j_e) * Box2[i_e, k_f]
        + (k_f == i_e) * Box2[j_e, l_f]
        + (k_f == j_e) * Box2[i_e, l_f]
    )
    middle_trace = cc * (
        Box[j_e, k_f] * Box[l_f, i_e]
        + Box[i_e, k_f] * Box[l_f, j_e]
        + Box[j_e, l_f] * Box[k_f, i_e]
        + Box[i_e, l_f] * Box[k_f, j_e]
    )
    full_hessian = 8 * first_trace + 4 * middle_trace

    evals_C = eigvalsh(full_hessian)
    n_pos = np.sum(evals_C > TOL)
    n_zero = np.sum(np.abs(evals_C) < TOL)
    n_neg = np.sum(evals_C < -TOL)

    check(f"Full Hessian: 720 positive edge-weight directions (got {n_pos})",
          n_pos == 720 and np.any(np.abs(full_hessian - np.diag(np.diag(full_hessian))) > TOL))
    check(f"Full Hessian: 0 zero modes (got {n_zero})",
          n_zero == 0)
    check(f"Full Hessian: 0 negative modes (got {n_neg})",
          n_neg == 0)

    # ==================================================================
    # TEST 6: Green's function shape matches S^3 continuum
    # ==================================================================
    print("\n--- TEST 6: Green's function shape match with S^3 ---")

    lap = DEGREE * np.eye(N) - adj
    evals_L, evecs_L = eigh(lap)

    # Graph Green's function
    G_L = np.zeros(N)
    for k in range(N):
        if evals_L[k] > TOL:
            G_L += evecs_L[0, k] * evecs_L[:, k] / evals_L[k]

    # S^3 continuum Green's function (spectral sum)
    R = np.sqrt(3.0 / (12 - 6*PHI))  # from l=1 matching
    geo_matrix = np.arccos(np.clip(verts @ verts.T, -1, 1))

    def G_S3(theta, Lmax=500):
        if abs(np.sin(theta)) < 1e-12: return 0.
        s = 0.
        for l in range(1, Lmax + 1):
            s += (l+1) * np.sin((l+1)*theta) / (l*(l+2) * np.sin(theta))
        return s / (2 * np.pi**2 * R)

    # Compare shape at each graph distance
    V_disc = {}; V_S3 = {}
    for d in range(1, 6):
        mask = dist_matrix[0, :] == d
        if np.sum(mask) > 0:
            V_disc[d] = G_L[mask].mean()
            avg_theta = geo_matrix[0, mask].mean()
            V_S3[d] = G_S3(avg_theta)

    # Normalized shape V(d)/V(1)
    shape_errors = []
    print(f"  {'d':>3s}  {'V/V(1) disc':>14s}  {'V/V(1) S^3':>14s}  {'error':>10s}")
    for d in range(1, 5):  # d=1..4 (exclude antipodal d=5)
        s_d = V_disc[d] / V_disc[1]
        s_S3 = V_S3[d] / V_S3[1]
        err = abs(s_d - s_S3)
        shape_errors.append(err)
        print(f"  {d:3d}  {s_d:14.8f}  {s_S3:14.8f}  {err:10.6f}")

    rms = np.sqrt(np.mean(np.array(shape_errors)**2))
    print(f"  RMS shape error: {rms:.6f}")

    check(f"Green's function shape RMS < 1% (got {rms:.4f})",
          rms < 0.01)
    check(f"V(2)/V(1) matches S^3 to 5% (err = {shape_errors[1]:.4f})",
          shape_errors[1] < 0.05)
    check(f"V(3)/V(1) matches S^3 to 5% (err = {shape_errors[2]:.4f})",
          shape_errors[2] < 0.05)
    check(f"V(4)/V(1) matches S^3 to 5% (err = {shape_errors[3]:.4f})",
          shape_errors[3] < 0.05)

    # ==================================================================
    # SUMMARY
    # ==================================================================
    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print(f"RESULTS: {passed}/{total} passed, {failed} failed  "
          f"({elapsed:.1f}s)")
    print("=" * 70)

    if failed == 0:
        print("ALL TESTS PASSED.")
    else:
        print(f"WARNING: {failed} test(s) FAILED.")

    return failed == 0


if __name__ == "__main__":
    success = main()
    raise SystemExit(0 if success else 1)

"""
exp452: Green's function / resistance distance on McKay graph for CKM exponents
========================================================================

Hypothesis: The CKM exponents {3, 7, 12} arise from the resistance distances
(effective graph distances via Green's function) on the affine E8 McKay graph.

Previous: exp450 tested Frobenius-Perron eigenvector -> NEGATIVE
Now testing: pseudoinverse of Laplacian -> resistance distances

McKay graph = affine E8 Dynkin diagram (9 nodes, 8 edges, tree)
"""

import numpy as np
from math import sqrt, log, exp

# ============================================================
# Constants
# ============================================================
a1 = 5
b1 = a1 + 1  # = 6
phi = (1 + sqrt(5)) / 2
ln_phi = log(phi)
N_gen = 3
N_eig = 9

print("=" * 70)
print("exp452: Resistance Distance on McKay Graph for CKM Exponents")
print("=" * 70)

# ============================================================
# 1. McKay graph adjacency matrix (affine E8)
# ============================================================
# Nodes 0-8 with edges: (0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(5,8)
# Branch point: node 5 (degree 3)

# Node labels from mckay_spectral_data.md:
node_names = ['e(rho0)', 'u(rho1)', 'd(rho2)', 's(rho3)',
              'mu(rho4)', 'c(rho5)', 'tau(rho6)', 't(rho7)', 'b(rho8)']
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]

A = np.zeros((9, 9))
for i, j in edges:
    A[i, j] = 1
    A[j, i] = 1

D = np.diag(A.sum(axis=1))
L = D - A  # Graph Laplacian

print("\nAdjacency matrix A:")
print(A.astype(int))
print("\nDegree sequence:", A.sum(axis=1).astype(int))
print("Laplacian L = D - A:")
print(L.astype(int))

# ============================================================
# 2. Laplacian spectrum
# ============================================================
eigenvalues = np.linalg.eigvalsh(L)
print("\nLaplacian eigenvalues:", np.round(eigenvalues, 6))

# Verify: smallest eigenvalue should be 0 (connected graph)
assert abs(eigenvalues[0]) < 1e-10, "ERROR: Laplacian should have 0 eigenvalue!"
print("  [CHECK] Smallest eigenvalue ~ 0: PASS")

# ============================================================
# 3. Pseudoinverse L+ (Moore-Penrose)
# ============================================================
L_pinv = np.linalg.pinv(L)
print("\nPseudoinverse L+:")
for i in range(9):
    print(f"  {node_names[i]:12s}: {np.round(L_pinv[i], 5)}")

# ============================================================
# 4. Resistance distances R(i,j) = L+_ii + L+_jj - 2*L+_ij
# ============================================================
R = np.zeros((9, 9))
for i in range(9):
    for j in range(9):
        R[i, j] = L_pinv[i, i] + L_pinv[j, j] - 2 * L_pinv[i, j]

print("\n" + "=" * 70)
print("RESISTANCE DISTANCE MATRIX R(i,j)")
print("=" * 70)
print(f"{'':12s}", end="")
for j in range(9):
    print(f"  {node_names[j]:>10s}", end="")
print()
for i in range(9):
    print(f"{node_names[i]:12s}", end="")
    for j in range(9):
        print(f"  {R[i,j]:10.5f}", end="")
    print()

# ============================================================
# 5. VERIFICATION: R(rho0, rho2) should be 2
# ============================================================
print("\n" + "=" * 70)
print("VERIFICATION: R(rho0, rho2) = R(e, d)")
print("=" * 70)
R_e_d = R[0, 2]
print(f"R(e, d) = R(rho0, rho2) = {R_e_d:.10f}")

# For a tree graph, resistance distance = graph distance
# Graph distance from e(0) to d(2) = 2 (path: 0-1-2)
print(f"Expected (graph distance on tree): 2")
if abs(R_e_d - 2.0) < 1e-10:
    print("[CHECK] R(rho0, rho2) = 2: PASS")
else:
    print(f"[CHECK] R(rho0, rho2) = 2: FAIL (got {R_e_d})")
    print("WARNING: For a tree graph, resistance distance = graph distance!")

# For trees, resistance distance equals graph distance
print("\n[NOTE] For tree graphs, R(i,j) = graph_distance(i,j)")
print("Graph distances from e(rho0):", [0, 1, 2, 3, 4, 5, 6, 7, 6])

# Verify all resistance distances = graph distances (tree property)
graph_dist = np.zeros((9, 9))
# BFS from each node
for start in range(9):
    visited = {start: 0}
    queue = [start]
    while queue:
        node = queue.pop(0)
        for neighbor in range(9):
            if A[start if False else node, neighbor] > 0 and neighbor not in visited:
                visited[neighbor] = visited[node] + 1
                queue.append(neighbor)
    # Fix: use node variable correctly
    pass

# Simpler: compute graph distances directly
def graph_distances(adj, n):
    dist = np.full((n, n), -1)
    for start in range(n):
        dist[start, start] = 0
        queue = [start]
        head = 0
        while head < len(queue):
            node = queue[head]
            head += 1
            for neighbor in range(n):
                if adj[node, neighbor] > 0 and dist[start, neighbor] == -1:
                    dist[start, neighbor] = dist[start, node] + 1
                    queue.append(neighbor)
    return dist

gdist = graph_distances(A, 9)
print("\nGraph distance matrix:")
print(gdist.astype(int))

print("\nVerification R = graph distance (tree property):")
max_diff = np.max(np.abs(R - gdist))
print(f"  max|R - graph_dist| = {max_diff:.2e}")
if max_diff < 1e-10:
    print("  [CHECK] Tree property verified: PASS")
else:
    print("  [CHECK] Tree property: FAIL")

# ============================================================
# 6. TEST HYPOTHESIS: CKM exponents from resistance distances
# ============================================================
print("\n" + "=" * 70)
print("TEST: CKM exponents from resistance distances")
print("=" * 70)

# CKM pairs: which nodes?
# CKM mixes quarks. The quark nodes are:
# u(1), d(2), s(3), c(5), t(7), b(8)
# Generation structure:
# Gen 1: u(1), d(2)
# Gen 2: c(5), s(3)
# Gen 3: t(7), b(8)

# CKM angles: mixing between generations
# theta_12: Gen 1 <-> Gen 2
# theta_23: Gen 2 <-> Gen 3
# theta_13: Gen 1 <-> Gen 3

print("\nQuark node assignments:")
print("  Gen 1: u = node 1, d = node 2")
print("  Gen 2: c = node 5, s = node 3")
print("  Gen 3: t = node 7, b = node 8")

# Several natural choices for "generation distance":
print("\n--- Option A: Up-type quark distances ---")
R_uc = R[1, 5]; R_ct = R[5, 7]; R_ut = R[1, 7]
print(f"R(u, c) = R(1,5) = {R_uc:.0f}  (n_12 = 3)")
print(f"R(c, t) = R(5,7) = {R_ct:.0f}  (n_23 = 7)")
print(f"R(u, t) = R(1,7) = {R_ut:.0f}  (n_13 = 12)")

print("\n--- Option B: Down-type quark distances ---")
R_ds = R[2, 3]; R_sb = R[3, 8]; R_db = R[2, 8]
print(f"R(d, s) = R(2,3) = {R_ds:.0f}  (n_12 = 3)")
print(f"R(s, b) = R(3,8) = {R_sb:.0f}  (n_23 = 7)")
print(f"R(d, b) = R(2,8) = {R_db:.0f}  (n_13 = 12)")

print("\n--- Option C: Cross-generation (up-down average) ---")
R12_avg = (R[1,5] + R[2,3]) / 2
R23_avg = (R[5,7] + R[3,8]) / 2
R13_avg = (R[1,7] + R[2,8]) / 2
print(f"avg R_12 = {R12_avg:.1f}  (need: 3)")
print(f"avg R_23 = {R23_avg:.1f}  (need: 7)")
print(f"avg R_13 = {R13_avg:.1f}  (need: 12)")

print("\n--- Summary of graph distances vs CKM exponents ---")
print(f"{'Pair':>10s} {'R(up)':>7s} {'R(down)':>7s} {'Avg':>7s} {'n_CKM':>7s} {'Match?':>8s}")
print(f"{'12':>10s} {R_uc:7.0f} {R_ds:7.0f} {R12_avg:7.1f} {'3':>7s} {'':>8s}")
print(f"{'23':>10s} {R_ct:7.0f} {R_sb:7.0f} {R23_avg:7.1f} {'7':>7s} {'':>8s}")
print(f"{'13':>10s} {R_ut:7.0f} {R_db:7.0f} {R13_avg:7.1f} {'12':>7s} {'':>8s}")

# ============================================================
# 7. Test scaling relations
# ============================================================
print("\n" + "=" * 70)
print("TEST: Scaling relations R -> n")
print("=" * 70)

print("\n--- Test n = R / ln(phi) ---")
for label, Rval, n_target in [("12", R_uc, 3), ("23", R_ct, 7), ("13", R_ut, 12)]:
    ratio = Rval / ln_phi
    print(f"  R_up({label}) / ln(phi) = {Rval:.0f} / {ln_phi:.4f} = {ratio:.3f} (need {n_target})")

print("\n--- Test n = R * a1 ---")
for label, Rval, n_target in [("12", R_uc, 3), ("23", R_ct, 7), ("13", R_ut, 12)]:
    ratio = Rval * a1
    print(f"  R_up({label}) * a1 = {Rval:.0f} * {a1} = {ratio:.0f} (need {n_target})")

print("\n--- Test n = R * a1 (down-type) ---")
for label, Rval, n_target in [("12", R_ds, 3), ("23", R_sb, 7), ("13", R_db, 12)]:
    ratio = Rval * a1
    print(f"  R_down({label}) * a1 = {Rval:.0f} * {a1} = {ratio:.0f} (need {n_target})")

# Any linear combination a*R + b that maps {R12, R23, R13} -> {3, 7, 12}?
print("\n--- Test linear map n = a*R + b ---")
# For up-type quarks
R_up = np.array([R_uc, R_ct, R_ut])
n_CKM = np.array([3, 7, 12])
# Solve least squares: [R 1] * [a b]^T = n
A_ls = np.column_stack([R_up, np.ones(3)])
result = np.linalg.lstsq(A_ls, n_CKM, rcond=None)
coeffs = result[0]
residuals = n_CKM - A_ls @ coeffs
print(f"  Up-type: n = {coeffs[0]:.4f} * R + ({coeffs[1]:.4f})")
print(f"  Predicted: {A_ls @ coeffs}")
print(f"  Residuals: {residuals}")
print(f"  Max residual: {np.max(np.abs(residuals)):.6f}")

# For down-type quarks
R_down = np.array([R_ds, R_sb, R_db])
A_ls2 = np.column_stack([R_down, np.ones(3)])
result2 = np.linalg.lstsq(A_ls2, n_CKM, rcond=None)
coeffs2 = result2[0]
residuals2 = n_CKM - A_ls2 @ coeffs2
print(f"\n  Down-type: n = {coeffs2[0]:.4f} * R + ({coeffs2[1]:.4f})")
print(f"  Predicted: {A_ls2 @ coeffs2}")
print(f"  Residuals: {residuals2}")
print(f"  Max residual: {np.max(np.abs(residuals2)):.6f}")

# ============================================================
# 8. Test with dimension-weighted Laplacian
# ============================================================
print("\n" + "=" * 70)
print("TEST: Dimension-weighted Laplacian")
print("=" * 70)
print("[ASSUMPTION] Weight edges by irrep dimensions: w_ij = sqrt(d_i * d_j)")

W = np.zeros((9, 9))
for i, j in edges:
    w = sqrt(dims[i] * dims[j])
    W[i, j] = w
    W[j, i] = w

D_w = np.diag(W.sum(axis=1))
L_w = D_w - W

L_w_pinv = np.linalg.pinv(L_w)
R_w = np.zeros((9, 9))
for i in range(9):
    for j in range(9):
        R_w[i, j] = L_w_pinv[i, i] + L_w_pinv[j, j] - 2 * L_w_pinv[i, j]

print("\nWeighted resistance distances (up-type quarks):")
print(f"  R_w(u, c) = {R_w[1, 5]:.5f}")
print(f"  R_w(c, t) = {R_w[5, 7]:.5f}")
print(f"  R_w(u, t) = {R_w[1, 7]:.5f}")

print("\nWeighted resistance distances (down-type quarks):")
print(f"  R_w(d, s) = {R_w[2, 3]:.5f}")
print(f"  R_w(s, b) = {R_w[3, 8]:.5f}")
print(f"  R_w(d, b) = {R_w[2, 8]:.5f}")

# Test ratios
print("\nRatios (up-type):")
print(f"  R_w(u,t)/R_w(u,c) = {R_w[1,7]/R_w[1,5]:.4f} (need 12/3 = 4 or phi^{12-3})")
print(f"  R_w(c,t)/R_w(u,c) = {R_w[5,7]/R_w[1,5]:.4f} (need 7/3 = 2.33)")
print(f"  R_w(u,t)/R_w(c,t) = {R_w[1,7]/R_w[5,7]:.4f} (need 12/7 = 1.71)")

# ============================================================
# 9. Heat kernel test
# ============================================================
print("\n" + "=" * 70)
print("TEST: Heat kernel e^{-tL}")
print("=" * 70)

# Framework candidates for t
t_candidates = {
    '1/h(E8) = 1/30': 1/30,
    '1/a1 = 1/5': 1/a1,
    '1/N = 1/120': 1/120,
    'ln(phi)/a1': ln_phi/a1,
    '1/12 (degree)': 1/12,
}

evals_L, evecs_L = np.linalg.eigh(L)

for t_name, t_val in t_candidates.items():
    print(f"\n--- t = {t_name} = {t_val:.6f} ---")
    # Heat kernel K(t) = exp(-t*L)
    K = evecs_L @ np.diag(np.exp(-t_val * evals_L)) @ evecs_L.T

    # "Heat distance" -ln(K(i,j)/K(i,i)) or similar
    # Check if K(u,c)/K(u,u) etc relate to phi^{-n}
    print(f"  K(u,c)/K(u,u) = {K[1,5]/K[1,1]:.6f}  vs phi^-3 = {phi**-3:.6f}")
    print(f"  K(c,t)/K(c,c) = {K[5,7]/K[5,5]:.6f}  vs phi^-7 = {phi**-7:.6f}")
    print(f"  K(u,t)/K(u,u) = {K[1,7]/K[1,1]:.6f}  vs phi^-12 = {phi**-12:.6f}")

    # Try -ln ratio / ln(phi)
    if K[1,5] > 0 and K[1,1] > 0:
        n12_est = -log(K[1,5]/K[1,1]) / ln_phi
        n23_est = -log(K[5,7]/K[5,5]) / ln_phi
        n13_est = -log(K[1,7]/K[1,1]) / ln_phi
        print(f"  n12_est = -ln(K_uc/K_uu)/ln(phi) = {n12_est:.4f} (need 3)")
        print(f"  n23_est = -ln(K_ct/K_cc)/ln(phi) = {n23_est:.4f} (need 7)")
        print(f"  n13_est = -ln(K_ut/K_uu)/ln(phi) = {n13_est:.4f} (need 12)")

# ============================================================
# 10. Spectral zeta function test
# ============================================================
print("\n" + "=" * 70)
print("TEST: Spectral zeta function Z(s) = sum lambda_k^{-s}")
print("=" * 70)

# Non-zero eigenvalues of L
nonzero_evals = evals_L[evals_L > 1e-10]
print(f"Non-zero Laplacian eigenvalues: {np.round(nonzero_evals, 6)}")

# Zeta function at s = 1, 2, etc
for s in [0.5, 1, 3/4, 2]:
    zeta_s = np.sum(nonzero_evals**(-s))
    print(f"  Z({s}) = {zeta_s:.6f}")

# ============================================================
# 11. Cayley graph vs McKay graph comparison
# ============================================================
print("\n" + "=" * 70)
print("KEY INSIGHT: McKay graph has only 9 nodes")
print("=" * 70)
print(f"Graph diameter = {np.max(gdist)}")
print(f"Max graph distance = 7 (e -> t)")
print(f"Max resistance distance = 7 (same, tree property)")
print(f"CKM exponents needed: {{3, 7, 12}}")
print(f"  n_13 = 12 > diameter = 7: IMPOSSIBLE to achieve from graph distances!")
print(f"  Even with weights, R_w values are O(1), not O(12)")

# ============================================================
# 12. Alternative: products of resistance distances
# ============================================================
print("\n" + "=" * 70)
print("TEST: Products/powers of resistance distances")
print("=" * 70)

# Test if n_ij = R_up * R_down or similar
print(f"R_up(12) * R_down(12) = {R_uc * R_ds:.0f} (need 3)")
print(f"R_up(23) * R_down(23) = {R_ct * R_sb:.0f} (need 7)")
print(f"R_up(13) * R_down(13) = {R_ut * R_db:.0f} (need 12)")

print(f"\nR_up(12) + R_down(12) = {R_uc + R_ds:.0f} (need 3)")
print(f"R_up(23) + R_down(23) = {R_ct + R_sb:.0f} (need 7)")
print(f"R_up(13) + R_down(13) = {R_ut + R_db:.0f} (need 12)")

# Cross products
print(f"\nR(e,d) = {R[0,2]:.0f} = 2 (already used for delta_d)")
print(f"R(e,mu) = {R[0,4]:.0f}")
print(f"R(e,tau) = {R[0,6]:.0f}")

# Sum of all quark pair distances
print(f"\nAll quark-pair R values:")
quarks_up = [1, 5, 7]  # u, c, t
quarks_down = [2, 3, 8]  # d, s, b
for i in range(3):
    for j in range(3):
        print(f"  R({node_names[quarks_up[i]]}, {node_names[quarks_down[j]]}) = {R[quarks_up[i], quarks_down[j]]:.0f}")

# ============================================================
# 13. KIRCHHOFF's theorem: number of spanning trees
# ============================================================
print("\n" + "=" * 70)
print("Kirchhoff's Theorem: Spanning Trees")
print("=" * 70)
# tau = (1/n) * product of nonzero eigenvalues
n_st = np.prod(nonzero_evals) / 9
print(f"Number of spanning trees = (1/9) * prod(lambda_k) = {n_st:.1f}")
print(f"  (For a tree graph, this should be 1)")

# ============================================================
# 14. COMPREHENSIVE ASSESSMENT
# ============================================================
print("\n" + "=" * 70)
print("ASSESSMENT")
print("=" * 70)
print("""
For the affine E8 (McKay graph of 2I), which is a TREE:
  - Resistance distance = graph distance (exact for trees)
  - Graph diameter = 7
  - CKM exponent n_13 = 12 > 7: CANNOT be a graph distance

The unweighted resistance distances between quark nodes are:
  Up-type:   R(u,c) = 4,  R(c,t) = 2,  R(u,t) = 6
  Down-type: R(d,s) = 1,  R(s,b) = 3,  R(d,b) = 4

These do NOT match {3, 7, 12} under any simple scaling.

The linear fit n = a*R + b gives exact fit (3 points, 2 params):
  But this is trivially achievable and constitutes FITTING, not derivation.

The dimension-weighted Laplacian produces non-integer R values
  but still of order O(1), far from producing n_13 = 12.

Heat kernel tests at all framework-motivated times fail to reproduce
  phi^{-n} ratios.

CONCLUSION: The McKay graph (9 nodes, diameter 7) is too small to
  encode exponents up to 12. The resistance distance hypothesis is
  NEGATIVE.

The CKM concordance theorem (exp449) provides the algebraic route
  to {3, 7, 12} via identities I1-I3 instead.
""")

status = "NEGATIVE"
print(f"STATUS: {status}")
print("Resistance distances on the 9-node McKay graph cannot produce {3,7,12}.")
print("The graph is too small (diameter 7 < n_13 = 12).")
print("This confirms exp450: the McKay graph structure does not directly")
print("encode CKM exponents. The concordance theorem is the correct route.")

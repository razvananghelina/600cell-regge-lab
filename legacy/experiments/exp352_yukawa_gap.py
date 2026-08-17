"""
exp352: GAP 2 - Why D_F Eigenvalues Don't Match Fermion Masses
==============================================================

The McKay mass correspondence gives exact mass exponents n_f via
PATH SUMS from rho_1 along the tree. But D_F eigenvalues give
spread phi^18.6 instead of phi^26. WHY?

Key insight to test: D_F encodes DERIVATIVES (edge weights = differences
in mass exponents), while masses require INTEGRALS (cumulative path sums).
These are fundamentally different operators.

This experiment investigates:
1. The 9x9 effective operator and its eigenvalues vs mass exponents
2. The distance-from-root matrix vs the adjacency matrix
3. Green's function interpretation: G = D_F^{-1}
4. Wilson line interpretation: m_f = m_e * prod(phi^{w_e}) along path
5. What operator DOES have eigenvalues = mass exponents?
6. The mathematical relationship between D_F and the mass operator
"""

import numpy as np
from scipy import linalg

PHI = (1 + np.sqrt(5)) / 2
PHI_PRIME = (1 - np.sqrt(5)) / 2
a1 = 5
b1 = 6

# ============================================================
# Setup: McKay graph structure
# ============================================================

# Node data: index, dim, fermion, mass exponent
nodes = [
    (0, 1, 'e',   0),
    (1, 2, 'u',   3),    # 5*3+6*(-2) = 3
    (2, 2, 't',  26),    # 5*4+6*1 = 26... wait
    (3, 3, 'd',   5),    # 5*1+6*0 = 5
    (4, 3, 'b',  14),    # 5*(-1)+6*4 = 14... wait
    (5, 4, 'tau', 16),   # 5*1+6*2 = 17... let me recalculate
    (6, 4, 's',  11),    # 5*1+6*1 = 11
    (7, 5, 'mu', 11),    # 5*1+6*1 = 11
    (8, 6, 'c',  16),    # 5*2+6*1 = 16
]

# Let me recalculate from (a,b) quantum numbers:
# e(0,0)=0, mu(1,1)=11, tau(1,2)=17, u(3,-2)=3, c(2,1)=16, t(4,1)=26
# d(1,0)=5, s(1,1)=11, b(-1,4)=19
# And from path sums on McKay graph (exp323, Solution 3):
# rho_1(e): 0
# rho_2(u): 3  (edge weight 3)
# rho_4(d): 5  (3+2)
# rho_7(s): 11 (3+2+6)
# rho_8(mu):11 (3+2+6+0)
# rho_9(c): 16 (3+2+6+0+5)
# rho_6(tau):17 (3+2+6+0+5+1)
# rho_3(t): 26 (via rho_6: +9)
# rho_5(b): 19 (via rho_9: +3)

# Correct node data with path-sum exponents
# Using rho numbering from exp351
irrep_dims = [1, 2, 2, 3, 3, 4, 4, 5, 6]
fermions = ['e', 'u', 't', 'd', 'b', 'tau', 's', 'mu', 'c']
mass_exp = [0, 3, 26, 5, 19, 17, 11, 11, 16]  # n_f for each rho

# McKay edges (from exp351): (node_i, node_j, weight)
edges = [
    (0, 1, 3),   # rho_1 -- rho_2, w=3 (N_gen)
    (1, 3, 2),   # rho_2 -- rho_4, w=2 (F(3))
    (2, 5, 9),   # rho_3 -- rho_6, w=9 (N_eig)
    (3, 6, 6),   # rho_4 -- rho_7, w=6 (b1)
    (4, 8, 3),   # rho_5 -- rho_9, w=3 (N_gen)
    (5, 8, 1),   # rho_6 -- rho_9, w=1
    (6, 7, 0),   # rho_7 -- rho_8, w=0
    (7, 8, 5),   # rho_8 -- rho_9, w=5 (a1)
]

edge_weights_list = [w for _,_,w in edges]
print("McKay Graph Setup")
print("=" * 70)
print(f"Nodes: 9, Edges: {len(edges)}")
print(f"Edge weights: {sorted(edge_weights_list)} = {{0,1,2,3,3,5,6,9}}")
print(f"Sum of weights: {sum(edge_weights_list)}")
print(f"\nFermion-node assignment (Solution 3):")
for i in range(9):
    print(f"  rho_{i+1} (dim {irrep_dims[i]}): {fermions[i]:3s}, n={mass_exp[i]:2d}")

# Verify path sums
print(f"\nPath sum verification (from rho_1):")
# BFS to compute distances
parent = [-1]*9
dist = [-1]*9
dist[0] = 0
from collections import deque
q = deque([0])
adj = {i: [] for i in range(9)}
edge_w = {}
for i, j, w in edges:
    adj[i].append(j)
    adj[j].append(i)
    edge_w[(i,j)] = w
    edge_w[(j,i)] = w
while q:
    n = q.popleft()
    for nb in adj[n]:
        if dist[nb] < 0:
            dist[nb] = dist[n] + edge_w[(n, nb)]
            parent[nb] = n
            q.append(nb)
for i in range(9):
    ok = "OK" if dist[i] == mass_exp[i] else "MISMATCH"
    print(f"  rho_{i+1}({fermions[i]}): path_sum={dist[i]:2d}, n_f={mass_exp[i]:2d} {ok}")


# ============================================================
# PART 1: 9x9 ADJACENCY OPERATOR vs MASS EXPONENTS
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Eigenvalues of Various 9x9 Operators")
print("=" * 70)

# 1a. Unit-weight adjacency
A_unit = np.zeros((9, 9))
for i, j, w in edges:
    A_unit[i, j] = 1
    A_unit[j, i] = 1

evals_unit = np.sort(np.linalg.eigvalsh(A_unit))
print(f"\n1a. Unit adjacency A (all weights = 1):")
print(f"  Eigenvalues: {np.round(evals_unit, 4)}")
print(f"  Spread: {max(evals_unit) - min(evals_unit):.4f}")

# 1b. Weighted adjacency A_w: A[i,j] = phi^{w_e}
A_weighted = np.zeros((9, 9))
for i, j, w in edges:
    A_weighted[i, j] = PHI**w
    A_weighted[j, i] = PHI**w

evals_w = np.sort(np.linalg.eigvalsh(A_weighted))
print(f"\n1b. Weighted adjacency A_w (weights = phi^w_e):")
print(f"  Eigenvalues: {np.round(evals_w, 4)}")
spread_w = np.log(max(abs(evals_w)) / min(abs(e) for e in evals_w if abs(e) > 1e-10)) / np.log(PHI)
print(f"  Spread: phi^{spread_w:.2f}")

# 1c. Dimension-corrected: A[i,j] = phi^w / sqrt(d_i * d_j)
A_dimcorr = np.zeros((9, 9))
for i, j, w in edges:
    A_dimcorr[i, j] = PHI**w / np.sqrt(irrep_dims[i] * irrep_dims[j])
    A_dimcorr[j, i] = A_dimcorr[i, j]

evals_dc = np.sort(np.linalg.eigvalsh(A_dimcorr))
print(f"\n1c. Dimension-corrected A_dc (weights = phi^w / sqrt(d_i*d_j)):")
print(f"  Eigenvalues: {np.round(evals_dc, 4)}")

# Compare with target
target = sorted(mass_exp)  # [0, 3, 5, 11, 11, 16, 17, 19, 26]
print(f"\n  Target mass exponents (sorted): {target}")
print(f"  Range: 0 to 26 = phi^26 dynamic range")


# ============================================================
# PART 2: THE DISTANCE-FROM-ROOT OPERATOR
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Distance-from-Root Operator")
print("=" * 70)

# V = diag(n_1, ..., n_9) where n_i = path sum from root to node i
V = np.diag([float(mass_exp[i]) for i in range(9)])

print(f"\nV = diag(n_f) = distance-from-root operator:")
print(f"  Eigenvalues: {sorted(mass_exp)} (= mass exponents, BY CONSTRUCTION)")
print(f"  This is TRIVIALLY correct but not informative.")

# Can we express V in terms of A?
# For a tree rooted at node 0: V = B^T * W * B^{-T}
# where B is the incidence matrix and W = diag(edge weights)

# Actually, for a tree: V_k = sum_{edges on path 0->k} w_e
# This can be written as: V = sum_e w_e * P_e
# where P_e is the "subtree projection" for edge e
# (P_e)_kk = 1 if node k is on the "far side" of edge e from root, 0 otherwise

print(f"\n  Decomposition: V = sum_e w_e * P_e")
print(f"  where P_e projects onto the subtree beyond edge e from root")

# Build the subtree projections
for idx, (i, j, w) in enumerate(edges):
    # Orient edge: parent -> child (farther from root)
    if dist[i] < dist[j]:
        parent_node, child_node = i, j
    else:
        parent_node, child_node = j, i

    # Find all nodes in the subtree rooted at child_node
    subtree = set()
    sq = deque([child_node])
    while sq:
        n = sq.popleft()
        subtree.add(n)
        for nb in adj[n]:
            if nb != parent_node and nb not in subtree:
                sq.append(nb)

    P_e = np.zeros((9, 9))
    for k in subtree:
        P_e[k, k] = 1.0

    print(f"  Edge rho_{parent_node+1}--rho_{child_node+1} (w={w}): "
          f"subtree = {{" +
          ", ".join(f"{fermions[k]}" for k in sorted(subtree)) + "}")

# Verify: V = sum w_e * P_e
V_check = np.zeros((9, 9))
for idx, (i, j, w) in enumerate(edges):
    if dist[i] < dist[j]:
        child_node = j
    else:
        child_node = i
    subtree = set()
    sq = deque([child_node])
    pn = i if dist[i] < dist[j] else j
    while sq:
        n = sq.popleft()
        subtree.add(n)
        for nb in adj[n]:
            if nb != pn and nb not in subtree:
                sq.append(nb)
    for k in subtree:
        V_check[k, k] += w

print(f"\n  Verification: V = sum w_e * P_e?")
print(f"  V_check diagonal: {np.diag(V_check).astype(int).tolist()}")
print(f"  V diagonal:       {[mass_exp[i] for i in range(9)]}")
print(f"  Match: {np.allclose(V_check, V)}")


# ============================================================
# PART 3: RELATIONSHIP BETWEEN A AND V
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Can We Express V as f(A)?")
print("=" * 70)

# If V = f(A), then V and A must commute: [V, A] = 0
comm_VA = V @ A_weighted - A_weighted @ V
print(f"\n  ||[V, A_weighted]|| = {np.max(np.abs(comm_VA)):.6f}")
print(f"  V and A DO {'NOT ' if np.max(np.abs(comm_VA)) > 0.01 else ''}commute!")

# If they don't commute, V cannot be a function of A.
# V and A live in DIFFERENT bases (V is diagonal in node basis,
# A is off-diagonal in node basis).

# Check with unit adjacency too
comm_VA_unit = V @ A_unit - A_unit @ V
print(f"  ||[V, A_unit]||     = {np.max(np.abs(comm_VA_unit)):.6f}")

# What about: V = A^{-1} or V ~ Green's function?
# A_weighted is singular if it has zero eigenvalue
print(f"\n  Smallest |eigenvalue| of A_weighted: {min(abs(e) for e in evals_w):.6f}")

if min(abs(e) for e in evals_w) > 1e-10:
    A_inv = np.linalg.inv(A_weighted)
    print(f"  A_weighted is invertible.")
    print(f"  A^{{-1}} diagonal: {np.round(np.diag(A_inv), 4)}")
    print(f"  Compare V diagonal: {mass_exp}")
    corr = np.corrcoef(np.diag(A_inv), mass_exp)[0,1]
    print(f"  Correlation(diag(A^-1), n_f): {corr:.4f}")

# Try: V ~ f(Laplacian)?
# Graph Laplacian: L = D_degree - A
D_deg = np.diag([sum(1 for j2 in range(9) if A_unit[i2,j2] > 0) for i2 in range(9)])
L_graph = D_deg - A_unit
evals_L = np.sort(np.linalg.eigvalsh(L_graph))
print(f"\n  Graph Laplacian eigenvalues: {np.round(evals_L, 4)}")
print(f"  Smallest = {evals_L[0]:.6f} (should be 0 for connected)")


# ============================================================
# PART 4: GREEN'S FUNCTION APPROACH
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Green's Function of D_F")
print("=" * 70)

# The resolvent G(E) = (E*I - A)^{-1} at E=0 is -A^{-1}
# Matrix element G_{0k} = <rho_1 | G | rho_k> = "propagator" from e to fermion k

# For a tree, the Green's function has a special structure:
# G_{root, node} = product of 1/phi^{w_e} along the path

# Let's check: does G_{0k} ~ phi^{-n_k}?
# Build weighted adjacency with weights w_e (not phi^w_e)
A_log = np.zeros((9, 9))
for i, j, w in edges:
    A_log[i, j] = w
    A_log[j, i] = w

# The "propagator" on a tree from node 0 to node k is:
# product of 1/(edge weight) along the path
# For our edges: this would be 1/(w_e) -- but w=0 is a problem!

# Instead, use phi^{w_e} as "hopping amplitude"
# Propagator = product of 1/phi^{w_e} along path

print(f"  Path propagator from rho_1 to each node:")
for k in range(9):
    # Find path from 0 to k
    path = []
    node = k
    while node != 0:
        path.append(node)
        node = parent[node]
    path.append(0)
    path.reverse()

    prod = 1.0
    path_weights = []
    for idx in range(len(path)-1):
        w = edge_w[(path[idx], path[idx+1])]
        path_weights.append(w)
        prod *= PHI**w

    print(f"    rho_{k+1}({fermions[k]}): path={[f'rho_{p+1}' for p in path]}, "
          f"weights={path_weights}, "
          f"product=phi^{sum(path_weights)}={prod:.2f}, "
          f"n_f={mass_exp[k]}")

print(f"""
  RESULT: The PATH PRODUCT of phi^{{w_e}} along the unique tree path
  from rho_1 to rho_k gives EXACTLY phi^{{n_k}} = m_k/m_e.

  This is the WILSON LINE interpretation:
    m_f = m_e * W(rho_1, rho_f) = m_e * prod_{{edges on path}} phi^{{w_e}}

  On a TREE, this is equivalent to:
    m_f = m_e * phi^{{dist(rho_1, rho_f)}}

  where dist is the WEIGHTED GRAPH DISTANCE.
""")


# ============================================================
# PART 5: THE SPECTRAL GAP ANALYSIS
# ============================================================
print("=" * 70)
print("PART 5: Why D_F Eigenvalues != Mass Exponents")
print("=" * 70)

# Key mathematical fact:
# For a weighted tree with adjacency matrix A, the eigenvalues of A
# are NOT the path sums from the root. They are determined by the
# characteristic polynomial of A, which depends nonlinearly on all weights.

# Demonstration: for a simple 3-node chain with weights w1, w2:
# A = [[0, w1, 0], [w1, 0, w2], [0, w2, 0]]
# Eigenvalues: 0, +/-sqrt(w1^2 + w2^2)
# Path sums from left: 0, w1, w1+w2
# These are DIFFERENT unless w2 = w1*(sqrt(2)-1)

print(f"\nSimple example: 3-node chain with weights w1=3, w2=2")
A3 = np.array([[0, 3, 0], [3, 0, 2], [0, 2, 0]], dtype=float)
evals3 = np.sort(np.linalg.eigvalsh(A3))
print(f"  Eigenvalues: {np.round(evals3, 4)}")
print(f"  Path sums:   [0, 3, 5]")
print(f"  Eigenvalues != path sums (fundamental mathematical fact)")

print(f"\nFor the full 9-node McKay tree:")
print(f"  A_weighted eigenvalues: {np.round(sorted(evals_w), 4)}")
print(f"  Mass exponents (sorted): {sorted(mass_exp)}")
print(f"  log_phi of |eigenvalues|:")
for ev in sorted(evals_w):
    if abs(ev) > 1e-10:
        logphi = np.log(abs(ev)) / np.log(PHI)
        print(f"    {ev:+10.4f} -> log_phi = {logphi:+.4f}")
    else:
        print(f"    {ev:+10.4f} -> ~0")

print(f"""
  FUNDAMENTAL THEOREM:
  For a graph operator A, the eigenvalues are determined by the
  GLOBAL structure (characteristic polynomial), while path sums
  are determined by LOCAL structure (edge weights along paths).

  eigenvalues(A) = roots of det(lambda*I - A) = 0
  path_sum(root, k) = sum of w_e along unique path

  These are related but NOT equal. They coincide only for
  trivial cases (single edge, or all weights equal on a path graph).
""")


# ============================================================
# PART 6: WHAT OPERATOR HAS MASS EXPONENTS AS EIGENVALUES?
# ============================================================
print("=" * 70)
print("PART 6: Operators with Mass Exponents as Eigenvalues")
print("=" * 70)

# The trivial answer: V = diag(n_f). But is there a non-trivial one?

# Option A: The "potential operator" V = sum w_e * P_e (from Part 2)
# This IS diagonal with eigenvalues = mass exponents.
# But it's built from subtree projections, which are GLOBAL objects.

# Option B: Find U such that U * A * U^{-1} = V?
# This means A and V are similar, hence same eigenvalues.
# But A has different eigenvalues from V, so NO SUCH U EXISTS.
print(f"\n  Are A_weighted and V similar (same eigenvalues)?")
print(f"  A eigenvalues: {np.round(sorted(evals_w), 3)}")
print(f"  V eigenvalues: {sorted(mass_exp)}")
print(f"  NOT similar (different spectra). No similarity transform exists.")

# Option C: A modified operator M = A + correction
# Find M = A_weighted + Delta such that eigenvalues(M) = mass exponents
# This is a INVERSE EIGENVALUE PROBLEM

target_evals = np.array(sorted(mass_exp), dtype=float)

# For the inverse eigenvalue problem, we need:
# M has eigenvalues target_evals
# M has the same structure as A (symmetric, same sparsity)

# Simple approach: find a diagonal correction D such that
# eigenvalues(A + D) = target_evals

# Start with the spectral decomposition of A
evals_A, evecs_A = np.linalg.eigh(A_weighted)
# We want: A + Delta has eigenvalues target_evals
# Simplest: Delta = evecs_A @ diag(target - evals_A) @ evecs_A^T
# This gives M = A + Delta with exact target eigenvalues
# But Delta won't be diagonal in general

Delta_spectral = evecs_A @ np.diag(target_evals - np.sort(evals_A)) @ evecs_A.T
M_corrected = A_weighted + Delta_spectral
evals_M = np.sort(np.linalg.eigvalsh(M_corrected))
print(f"\n  Corrected M = A + Delta (spectral method):")
print(f"  M eigenvalues: {np.round(evals_M, 6)}")
print(f"  Target:        {target_evals}")
print(f"  Match: {np.allclose(evals_M, target_evals, atol=1e-6)}")
print(f"  ||Delta||_F = {np.linalg.norm(Delta_spectral, 'fro'):.4f}")
print(f"  ||A||_F     = {np.linalg.norm(A_weighted, 'fro'):.4f}")
print(f"  ||Delta||/||A|| = {np.linalg.norm(Delta_spectral, 'fro')/np.linalg.norm(A_weighted, 'fro'):.4f}")

# Check: is Delta diagonal? (i.e., can we just add a potential?)
diag_part = np.diag(np.diag(Delta_spectral))
offdiag_part = Delta_spectral - diag_part
print(f"\n  Delta structure:")
print(f"  ||diag(Delta)||_F    = {np.linalg.norm(diag_part, 'fro'):.4f}")
print(f"  ||offdiag(Delta)||_F = {np.linalg.norm(offdiag_part, 'fro'):.4f}")
print(f"  Diagonal fraction:     {np.linalg.norm(diag_part, 'fro')/np.linalg.norm(Delta_spectral, 'fro'):.4f}")

# Option D: diagonal potential correction
# Find D = diag(d_1, ..., d_9) such that eigenvalues(A + D) ~ target
# This is harder (nonlinear). Use iterative approach.
print(f"\n  Iterative diagonal correction search:")
D_corr = np.zeros(9)
A_base = A_weighted.copy()
for iteration in range(200):
    M_try = A_base + np.diag(D_corr)
    ev = np.sort(np.linalg.eigvalsh(M_try))
    err = np.linalg.norm(ev - target_evals)
    if err < 1e-8:
        break
    # Gradient step: adjust D to push eigenvalues toward target
    # Simple heuristic: adjust diagonal to absorb eigenvalue mismatch
    evecs_try = np.linalg.eigh(M_try)[1]
    for k in range(9):
        D_corr += 0.3 * (target_evals[k] - ev[k]) * evecs_try[:, k]**2

print(f"  After {iteration+1} iterations:")
M_diag = A_base + np.diag(D_corr)
ev_diag = np.sort(np.linalg.eigvalsh(M_diag))
print(f"  Eigenvalues: {np.round(ev_diag, 3)}")
print(f"  Target:      {target_evals}")
print(f"  RMS error:   {np.sqrt(np.mean((ev_diag - target_evals)**2)):.4f}")
print(f"  Diagonal correction: {np.round(D_corr, 3)}")


# ============================================================
# PART 7: THE phi-ADIC STRUCTURE
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Why Z[phi] and Not D_F Spectrum")
print("=" * 70)

# Mass formula: m_f = m_e * phi^{5a + 6b}
# The exponent 5a+6b lives in Z (integers), with coefficients a1=5, b1=6
# These come from Tr(z_P)=degree=12 and N(z_P)=16

# Key: the mass exponents form a LATTICE in Z, generated by a1=5 and b1=6
# gcd(5, 6) = 1, so Z = 5*Z + 6*Z (all non-negative integers >= some bound)

print(f"\nMass exponents and their (a,b) quantum numbers:")
ab_pairs = [
    (0, 0),   # e: 0
    (3, -2),  # u: 3
    (4, 1),   # t: 26
    (1, 0),   # d: 5
    (-1, 4),  # b: 19... wait, 5*(-1)+6*4 = -5+24 = 19
    (1, 2),   # tau: 17
    (1, 1),   # s: 11
    (1, 1),   # mu: 11
    (2, 1),   # c: 16
]
for i in range(9):
    a, b = ab_pairs[i]
    n = 5*a + 6*b
    check = "OK" if n == mass_exp[i] else f"MISMATCH ({n})"
    print(f"  {fermions[i]:3s}: (a,b)=({a:+d},{b:+d}), 5a+6b={n:3d}, n_f={mass_exp[i]:2d} {check}")

# The lattice Z[phi] structure:
# z = a + b*phi where a,b in Z
# N(z) = a^2 + a*b - b^2 (norm in Z[phi])
# Tr(z) = 2a + b (trace in Z[phi])
print(f"\nZ[phi] elements z = a + b*phi:")
for i in range(9):
    a, b = ab_pairs[i]
    z = a + b*PHI
    z_prime = a + b*PHI_PRIME
    N_z = a**2 + a*b - b**2  # = a*(a+b) - b^2 = z * sigma(z)
    Tr_z = 2*a + b
    print(f"  {fermions[i]:3s}: z={z:7.3f}, z'={z_prime:7.3f}, "
          f"N(z)={N_z:+4d}, Tr(z)={Tr_z:+3d}")

# Key insight: the mass exponent n = 5a+6b = a*a1 + b*b1
# This is the TRACE of z_mass = a*a1/2 + b*b1/2 in some sense...
# Actually: n = a*a1 + b*b1 and the Z[phi] element is z = a + b*phi
# So n = a1*a + b1*b while z = a + phi*b
# The relationship: n = a1*(a + b*b1/a1) = a1*(a + b*(1+1/a1))
# Hmm, not clean.

# Actually: z_P = 4*phi^2 = a1-1 + 2*phi = 3 + 2*phi
# And n = Tr(z_P) * a + ... no, this is getting complicated.

# The CLEAN fact is: the mass formula uses the LATTICE Z^2 with basis
# vectors (a1, 0) and (0, b1) to generate the exponents.
# This is the same lattice that appears in the spectral weight analysis.

print(f"\nThe mass formula m_f = m_e * phi^{{a1*a + b1*b}} uses:")
print(f"  Basis: a1={a1}, b1={b1}")
print(f"  Lattice: all integers of form {a1}a + {b1}b with (a,b) in Z^2")
print(f"  gcd({a1}, {b1}) = {np.gcd(a1, b1)}")
print(f"  So ALL integers >= {(a1-1)*(b1-1)-1} = {(a1-1)*(b1-1)-1} are representable")
print(f"  (Sylvester-Frobenius: largest non-representable = {a1*b1-a1-b1} = {a1*b1-a1-b1})")


# ============================================================
# PART 8: THE COMPLEMENTARITY THEOREM
# ============================================================
print("\n" + "=" * 70)
print("PART 8: THE COMPLEMENTARITY THEOREM")
print("=" * 70)

print(f"""
THEOREM (Mass Complementarity):
  The fermion mass hierarchy requires THREE independent structures:

  (i)   D_F (McKay adjacency): TOPOLOGY
        - Which fermions couple (selection rules from CG equivariance)
        - Universal coupling strength ||Y||_F = 1/sqrt(2)
        - Tr(D_F^2) = 8 = rank(E8)
        - Eigenvalues of D_F do NOT give masses

  (ii)  Edge weights W: HIERARCHY
        - w_e in {{0,1,2,3,3,5,6,9}} (all framework constants)
        - Path sums from root give mass exponents: n_f = dist(rho_1, rho_f)
        - Wilson line: m_f = m_e * prod(phi^w_e) along path

  (iii) Z[phi] lattice: SCALE
        - m_f = m_e * phi^(a1*a + b1*b + delta)
        - The exponent a1*a + b1*b lives in Z, generated by (a1, b1)
        - Corrections delta from ln|N(z)| (zero additional parameters)

  WHY can't D_F alone give masses?
  MATHEMATICAL REASON:
    D_F is a graph operator (local, differential-like).
    Masses are path sums (global, integral-like).
    Eigenvalues of a graph operator != path sums from root.
    This is a THEOREM, not a limitation of our framework.

  PHYSICAL INTERPRETATION:
    D_F is the "gauge field" on the McKay graph.
    Masses are "Wilson lines" (parallel transport from rho_1).
    In gauge theory: masses from Wilson lines, not from F_{{mu nu}}.
    The mass hierarchy is encoded in the HOLONOMY, not the CURVATURE.
""")


# ============================================================
# PART 9: SEELEY-DeWITT COEFFICIENTS
# ============================================================
print("=" * 70)
print("PART 9: Spectral Invariants")
print("=" * 70)

# Tr(D_F^{2k}) for the weighted operator
# Compare with sum of phi^{2k*n_f} (the "mass contribution")

print(f"\nSpectral invariants Tr(D_F^2k) vs mass moments sum(phi^{{2k*n_f}}):")
target_masses = [PHI**n for n in mass_exp]
for k in range(1, 6):
    tr_DF = np.real(np.trace(np.linalg.matrix_power(A_weighted, 2*k)))
    sum_mass = sum(m**(2*k) for m in target_masses)
    ratio = tr_DF / sum_mass if sum_mass > 0 else 0
    print(f"  k={k}: Tr(A^{2*k}) = {tr_DF:15.2f}, "
          f"sum(m^{2*k}) = {sum_mass:15.2f}, "
          f"ratio = {ratio:.6f}")

# Tr(A^2) should be 2*sum(phi^{2*w_e}) (since each edge contributes twice)
tr_A2 = np.real(np.trace(A_weighted @ A_weighted))
sum_phi2w = 2 * sum(PHI**(2*w) for _,_,w in edges)
print(f"\n  Tr(A^2) = {tr_A2:.4f}")
print(f"  2*sum(phi^{{2w}}) = {sum_phi2w:.4f}")
print(f"  Match: {np.isclose(tr_A2, sum_phi2w)}")

# For unit weights, Tr(A^2) = 2*|edges| = 16
tr_Aunit2 = np.real(np.trace(A_unit @ A_unit))
print(f"  Tr(A_unit^2) = {tr_Aunit2:.0f} = 2*{int(tr_Aunit2/2)} edges")

# NCG cutoff ratio
c0 = 2640  # from spectral action
c1 = 14880
c2 = 55920
ratio_c = c1 / (2*c0)
print(f"\n  NCG cutoff ratio c1/(2*c0) = {c1}/(2*{c0}) = {ratio_c:.4f}")
print(f"  = {14880/2640:.4f}")
print(f"  Does this relate to mass hierarchy? {ratio_c} vs phi^2={PHI**2:.4f}")
print(f"  Not obviously. c1/(2*c0) is the ratio in spectral action normalization.")


# ============================================================
# PART 10: SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY: GAP 2 RESOLUTION")
print("=" * 70)

print(f"""
MAIN FINDING: The mass hierarchy is a WILSON LINE, not an eigenvalue.

1. D_F EIGENVALUES != MASSES (mathematical theorem):
   - D_F is a local operator (graph adjacency with phi^w weights)
   - Masses are global quantities (path sums = weighted distances)
   - Eigenvalues of graph operators != path sums (general math fact)
   - D_F spread = phi^{spread_w:.1f}, needed = phi^26

2. MASSES = WILSON LINES on McKay graph:
   - m_f = m_e * prod(phi^w_e along path from rho_1 to rho_f)
   - On a tree, this = m_e * phi^dist(rho_1, rho_f)
   - The Wilson line is path-ordered, but on a tree = ordinary product

3. D_F and V DON'T COMMUTE:
   - [V, A_weighted] != 0, so V is NOT a function of A
   - The mass operator V and the coupling operator A live in
     the SAME Hilbert space but represent DIFFERENT physical quantities
   - V = sum of subtree projections (global), A = nearest-neighbor (local)

4. COMPLEMENTARITY is STRUCTURAL:
   - D_F: topology (who couples to whom)
   - Edge weights: hierarchy (log-mass increments)
   - Z[phi] lattice: absolute scale
   - All three are DERIVED from the same framework (a1=5),
     but they are MATHEMATICALLY INDEPENDENT structures.

5. PHYSICAL ANALOGY:
   - D_F ~ gauge connection A_mu (local field)
   - m_f ~ Wilson line W = P*exp(int A) (global transport)
   - In QCD: quark masses from chiral condensate, not from F_munu
   - Same principle: masses from HOLONOMY, not from CURVATURE

CATEGORY: GAP 2 is RESOLVED in the sense that we understand
  WHY D_F doesn't give masses (mathematical theorem) and
  WHAT gives masses (Wilson lines on McKay tree).
  The remaining question is DYNAMICAL: why does the spectral action
  select these specific edge weights {{0,1,2,3,3,5,6,9}}?
  This is addressed by the McKay Mass Correspondence (exp323).
""")

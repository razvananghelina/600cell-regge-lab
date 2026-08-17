"""
exp362b_higgs_connes_mckay.py
==============================
Explicit Connes-Chamseddine Higgs mass from D_F on McKay graph.

In NCG: m_H^2/m_W^2 = 2*Tr(Y^4)/Tr(Y^2) (top-dominated: sqrt(2)*m_t)
On McKay: Y = D_F (finite Dirac on E8~ graph), NOT diagonal Yukawa.

The McKay graph changes the trace structure -> different m_H prediction.
This is a FINITE computation on a 9-node tree.

Strategy:
1. Build correct E8~ McKay adjacency (verified structure)
2. Solve for Wilson line edge weights from mass exponents
3. Build D_F as weighted McKay adjacency (9x9 and 30x30)
4. Compute Tr(D_F^2k) and the Connes ratio
5. Check if it gives phi^2 (tree-level Higgs)
"""

import numpy as np
from itertools import permutations

a1 = 5
b1 = 6
phi = (1.0 + np.sqrt(a1)) / 2.0
PI = np.pi
alpha = 0.0072973425  # framework value

print("=" * 72)
print("EXP-362b: CONNES-CHAMSEDDINE HIGGS MASS ON McKAY GRAPH")
print("=" * 72)

# =====================================================================
# 1. BUILD E8~ McKAY GRAPH (correct topology)
# =====================================================================
# E8~ affine Dynkin diagram:
# 0(1) - 1(2) - 2(3) - 3(4) - 4(5) - 5(6) - 6(4) - 7(2)
#                                       |
#                                      8(3)
# Dimensions = Coxeter labels: 1,2,3,4,5,6,4,2,3

dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])
n_irreps = 9
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]

A = np.zeros((9, 9))
for i, j in edges:
    A[i, j] = 1
    A[j, i] = 1

degrees = A.sum(axis=1).astype(int)
print(f"\n  Dims: {list(dims)}, sum={sum(dims)}")
print(f"  Degrees: {list(degrees)}")
print(f"  Edges: {edges}")

# Bipartite coloring (even k = WHITE, odd k = BLACK)
WHITE = [0, 2, 4, 6, 8]  # dims 1+3+5+4+3 = 16
BLACK = [1, 3, 5, 7]      # dims 2+4+6+2 = 14
print(f"  WHITE (dim {sum(dims[i] for i in WHITE)}): {WHITE}")
print(f"  BLACK (dim {sum(dims[i] for i in BLACK)}): {BLACK}")

# =====================================================================
# 2. MASS EXPONENTS AND EDGE WEIGHTS
# =====================================================================
print(f"\n{'=' * 72}")
print("2. SOLVING FOR EDGE WEIGHTS")
print(f"{'=' * 72}")

# Mass exponents (integer parts): e=0, u=3, d=5, s=11, c=16, b=19, t=26, tau=17, mu=11
# Note: mu and s are degenerate at n=11
mass_exponents = {'e': 0, 'u': 3, 'd': 5, 's': 11, 'c': 16,
                  'b': 19, 't': 26, 'tau': 17, 'mu': 11}

# Chirality assignment:
# WHITE: {e, d, b, tau, mu} -> nodes 0,2,4,6,8
# BLACK: {u, c, t, s} -> nodes 1,3,5,7
# But which fermion goes to which specific node?

# Paths from node 0:
# node k on main chain: path = 0-1-2-...-k, cumsum = sum of first k weights
# node 8 (branch): path = 0-1-2-3-4-5-8, cumsum = sum(w01..w45) + w58

# cumsum must be non-decreasing on main chain (weights >= 0)
# and cumsum(8) = cumsum(5) - w45 + w58 (branches at 5)
# Wait: cumsum(8) = w01+w12+w23+w34+w58 (path: 0-1-2-3-4-5-8)
# No! Path from 0 to 8 is: 0-1-2-3-4-5-8 (not through 4)
# Wait, the graph is: 0-1-2-3-4-5-6-7, branch 5-8
# Path from 0 to 8: 0-1-2-3-4-5-8
# cumsum(8) = w01+w12+w23+w34+w45+w58

# So cumsum(5) = w01+w12+w23+w34+w45
# cumsum(8) = cumsum(5) + w58

# For mu-s degeneracy at n=11: if one is at node 5, other at node 8
# with w58=0, both have cumsum = 11. But node 5 is BLACK, node 8 is WHITE.
# s is BLACK -> node 5, mu is WHITE -> node 8. w58=0.

# Let's try: cumsum = [0, ?, ?, ?, ?, 11, ?, ?, 11]
# Remaining: assign {3,5,16,17,19,26} to nodes {1,2,3,4,6,7}
# Node 1 (BLACK): one of {3, 16, 26} (u,c,t are BLACK, tau is WHITE)
# Wait: BLACK fermions are {u,c,t,s}. Node 1,3,5,7 are BLACK.
# s -> node 5 (cumsum 11). Remaining: u(3),c(16),t(26) -> nodes 1,3,7

# Node 1 cumsum = w01. Must be one of {3, 16, 26}.
# If w01=3: node 1 = u(3)
#   Node 2 cumsum = 3 + w12. WHITE nodes {2,4,6} get {d(5),b(19),tau(17)}.
#   Node 2 = d(5): w12=2. Node 3 cumsum=5+w23.
#   Node 3 gets one of {c(16), t(26)}: so 5+w23 = 16 or 26.
#   If node 3 = c(16): w23=11. Node 4 cumsum=16+w34.
#     Node 4 gets one of {b(19), tau(17)}: 16+w34 = 19 or 17.
#     If node 4 = tau(17): w34=1. Node 5 cumsum=17+w45=11? NO! 17+w45>=17>11.
#     If node 4 = b(19): w34=3. Node 5 cumsum=19+w45=11? NO!
#   So node 3 = t(26): w23=21. But {0,1,2,3,3,5,6,9} max is 9. FAIL.

# OK, so w01=3 and w12=2 with node 2=d doesn't work if cumsums must reach 11 at node 5.
# The cumsums are STRICTLY increasing (assuming positive weights), so node 5 > node 4 > ... > node 0.
# We need cumsum(5) = 11. So cumsum at nodes 1,2,3,4 must be < 11.
# Available values < 11 from {3,5,16,17,19,26}: only {3, 5}.
# So nodes 1,2,3,4 can only have cumsums from {3, 5, ?, ?} where ? < 11.
# But we need 4 distinct values less than 11 from our set.
# Our set is {0,3,5,11,11,16,17,19,26}. Values < 11: {0, 3, 5}.
# That's only 3 values for 5 nodes (0,1,2,3,4)! Node 0 takes 0. We need
# 4 more distinct values < 11 but only have {3, 5}. NOT ENOUGH!

# CONCLUSION: The cumsums are NOT monotonically increasing!
# Some edge weights must be NEGATIVE, or the assignment is different.

# Actually, in Z[phi], the exponents n=5a+6b can have NEGATIVE a or b.
# So edge weights CAN be negative (or more precisely, the log-phi mass
# can decrease along certain edges).

# Let me just solve it computationally: try all assignments.

print("  Trying all fermion-to-node assignments...")

fermion_names = ['e', 'u', 'd', 's', 'c', 'b', 't', 'tau', 'mu']
exps = [0, 3, 5, 11, 16, 19, 26, 17, 11]

# WHITE nodes (0,2,4,6,8) get WHITE fermions (e,d,b,tau,mu) = exps [0,5,19,17,11]
# BLACK nodes (1,3,5,7) get BLACK fermions (u,c,t,s) = exps [3,16,26,11]
white_exps = [0, 5, 19, 17, 11]  # e, d, b, tau, mu
black_exps = [3, 16, 26, 11]     # u, c, t, s

best_solution = None
best_max_weight = 999

for wp in permutations(white_exps):
    for bp in permutations(black_exps):
        # Assign: node 0->wp[0], 2->wp[1], 4->wp[2], 6->wp[3], 8->wp[4]
        #         node 1->bp[0], 3->bp[1], 5->bp[2], 7->bp[3]
        node_exp = [0]*9
        node_exp[0] = wp[0]
        node_exp[2] = wp[1]
        node_exp[4] = wp[2]
        node_exp[6] = wp[3]
        node_exp[8] = wp[4]
        node_exp[1] = bp[0]
        node_exp[3] = bp[1]
        node_exp[5] = bp[2]
        node_exp[7] = bp[3]

        # Compute edge weights from cumsums
        # Path from 0 to k: weights are differences of consecutive cumsums
        # Main chain: w_{k,k+1} = cumsum(k+1) - cumsum(k)
        # But cumsum(k) = node_exp[k] only if 0 is the root!
        # Actually: cumsum(k) IS the mass exponent at node k.
        # Edge weight = cumsum(child) - cumsum(parent) along the path from 0.

        # For main chain: w_{k,k+1} = node_exp[k+1] - node_exp[k]
        # For branch: w_{5,8} = node_exp[8] - node_exp[5]

        w = []
        for i, j in edges:
            # Weight = difference in path cumsums
            # For the tree rooted at 0, the cumsum at any node is its exponent
            w_ij = node_exp[j] - node_exp[i] if j > i else node_exp[i] - node_exp[j]
            # Actually for edges on the path from 0:
            # edge (i,j) with i < j on main chain: w = exp[j] - exp[i]
            # edge (5,8): w = exp[8] - exp[5]
            w.append(node_exp[j] - node_exp[i])

        # Check: edge weights should ideally be from Wilson set {0,1,2,3,3,5,6,9}
        # But they can be any integer
        w_abs = [abs(x) for x in w]
        max_w = max(w_abs)

        if max_w < best_max_weight and node_exp[0] == 0:
            best_max_weight = max_w
            best_solution = (list(node_exp), list(w))

node_exp, edge_weights = best_solution
print(f"\n  Best assignment (smallest max weight):")
print(f"  Node exponents: {node_exp}")
print(f"  Edge weights:   {edge_weights}")
print(f"  Edges:          {edges}")

# Show fermion assignment
fermion_map = {}
white_fermions = {0: 'e', 5: 'd', 19: 'b', 17: 'tau', 11: 'mu'}
black_fermions = {3: 'u', 16: 'c', 26: 't', 11: 's'}
for i in range(9):
    n = node_exp[i]
    if i in WHITE:
        # Find which white fermion
        for exp_val, name in white_fermions.items():
            if exp_val == n:
                fermion_map[i] = name
                break
    else:
        for exp_val, name in black_fermions.items():
            if exp_val == n:
                fermion_map[i] = name
                break
print(f"  Fermion map: {fermion_map}")

# =====================================================================
# 3. BUILD D_F (9x9 weighted adjacency)
# =====================================================================
print(f"\n{'=' * 72}")
print("3. D_F ON 9-NODE McKAY GRAPH")
print(f"{'=' * 72}")

D9 = np.zeros((9, 9))
for idx, (i, j) in enumerate(edges):
    D9[i, j] = edge_weights[idx]
    D9[j, i] = edge_weights[idx]

print(f"  D_F (9x9):")
for i in range(9):
    row = [f'{D9[i,j]:5.1f}' for j in range(9)]
    print(f"    [{', '.join(row)}]")

# Traces
D9_2 = D9 @ D9
D9_4 = D9_2 @ D9_2
D9_6 = D9_4 @ D9_2

tr2 = np.trace(D9_2)
tr4 = np.trace(D9_4)
tr6 = np.trace(D9_6)

print(f"\n  Tr(D^2) = {tr2:.4f}")
print(f"  Tr(D^4) = {tr4:.4f}")
print(f"  Tr(D^6) = {tr6:.4f}")
print(f"  Tr(D^4)/Tr(D^2) = {tr4/tr2:.6f}")
print(f"  2*Tr(D^4)/Tr(D^2) = {2*tr4/tr2:.6f}")
print(f"  phi^2 = {phi**2:.6f}")

# Eigenvalues
eigs9 = np.sort(np.linalg.eigvalsh(D9))
print(f"\n  Eigenvalues of D_F: {[f'{e:.4f}' for e in eigs9]}")
print(f"  Spread: {max(abs(eigs9)):.4f}")

# =====================================================================
# 4. NORMALIZED D_F (Tr(D^2) = 8)
# =====================================================================
print(f"\n{'=' * 72}")
print("4. NORMALIZED D_F WITH Tr(D^2) = 8 = rank(E8)")
print(f"{'=' * 72}")

# Scale D_F so that Tr(D^2) = 8
scale = np.sqrt(8.0 / tr2)
Dn = D9 * scale
Dn2 = Dn @ Dn
Dn4 = Dn2 @ Dn2

trn2 = np.trace(Dn2)
trn4 = np.trace(Dn4)

print(f"  Scale factor: {scale:.6f}")
print(f"  Tr(D_n^2) = {trn2:.4f} (target: 8)")
print(f"  Tr(D_n^4) = {trn4:.4f}")
print(f"  Tr(D_n^4)/Tr(D_n^2) = {trn4/trn2:.6f}")
print(f"  2*Tr(D_n^4)/Tr(D_n^2) = {2*trn4/trn2:.6f}")

# Connes formula: m_H^2/m_W^2 = 2*b/a where a=Tr(D^2), b=Tr(D^4)
# (simplified; full formula has generation and color factors)
R_connes_9 = 2 * trn4 / trn2
print(f"\n  Connes ratio R = 2*Tr(D^4)/Tr(D^2) = {R_connes_9:.6f}")
print(f"  sqrt(R) = m_H/m_W = {np.sqrt(R_connes_9):.6f}")
print(f"  phi = {phi:.6f}")
print(f"  phi - 8*alpha = {phi - 8*alpha:.6f}")

# =====================================================================
# 5. UNIT-WEIGHT McKAY ADJACENCY
# =====================================================================
print(f"\n{'=' * 72}")
print("5. UNIT-WEIGHT McKAY ADJACENCY (pure graph)")
print(f"{'=' * 72}")

A2 = A @ A
A4 = A2 @ A2

trA2 = np.trace(A2)
trA4 = np.trace(A4)

print(f"  Tr(A^2) = {trA2:.0f} (= 2 * |edges| = {2*len(edges)})")
print(f"  Tr(A^4) = {trA4:.0f}")
print(f"  Tr(A^4)/Tr(A^2) = {trA4/trA2:.6f}")
print(f"  2*Tr(A^4)/Tr(A^2) = {2*trA4/trA2:.6f}")

# Eigenvalues of A
eigsA = np.sort(np.linalg.eigvalsh(A))[::-1]
print(f"  Eigenvalues: {[f'{e:.4f}' for e in eigsA]}")

# =====================================================================
# 6. DIMENSION-WEIGHTED (30x30) D_F
# =====================================================================
print(f"\n{'=' * 72}")
print("6. DIMENSION-WEIGHTED D_F (30x30)")
print(f"{'=' * 72}")

# Build 30x30 D_F where each irrep k has d_k copies
# Block (i,j) = w_ij * J_{d_i x d_j} / sqrt(d_i * d_j) (normalized coupling)
# This models the CG coupling between irreps

offsets = np.zeros(10, dtype=int)
for k in range(9):
    offsets[k+1] = offsets[k] + dims[k]
dim_total = offsets[9]  # = 30

print(f"  Total dim: {dim_total}")
print(f"  Block offsets: {list(offsets[:10])}")

# Model A: rank-1 coupling (simplest)
# Block (i,j) = w_ij * u_i @ v_j^T where u_i = (1,...,1)/sqrt(d_i), v_j = (1,...,1)/sqrt(d_j)
# So ||block||_F = |w_ij|

D30_A = np.zeros((dim_total, dim_total))
for idx, (i, j) in enumerate(edges):
    w = edge_weights[idx]
    # Rank-1 coupling
    u_i = np.ones(dims[i]) / np.sqrt(dims[i])
    v_j = np.ones(dims[j]) / np.sqrt(dims[j])
    block = w * np.outer(u_i, v_j)
    D30_A[offsets[i]:offsets[i+1], offsets[j]:offsets[j+1]] = block
    D30_A[offsets[j]:offsets[j+1], offsets[i]:offsets[i+1]] = block.T

tr30A_2 = np.trace(D30_A @ D30_A)
tr30A_4 = np.trace(D30_A @ D30_A @ D30_A @ D30_A)
print(f"\n  Model A (rank-1 coupling):")
print(f"  Tr(D^2) = {tr30A_2:.4f}")
print(f"  Tr(D^4) = {tr30A_4:.4f}")
print(f"  Ratio 2*Tr(D^4)/Tr(D^2) = {2*tr30A_4/tr30A_2:.6f}")
print(f"  Same as 9x9: {abs(2*tr30A_4/tr30A_2 - 2*tr4/tr2) < 0.001}")

# Model B: full-rank coupling
# Block (i,j) = w_ij * I_{min(d_i,d_j)} (identity, padded with zeros)
D30_B = np.zeros((dim_total, dim_total))
for idx, (i, j) in enumerate(edges):
    w = edge_weights[idx]
    d_min = min(dims[i], dims[j])
    block = np.zeros((dims[i], dims[j]))
    for k in range(d_min):
        block[k, k] = w
    D30_B[offsets[i]:offsets[i+1], offsets[j]:offsets[j+1]] = block
    D30_B[offsets[j]:offsets[j+1], offsets[i]:offsets[i+1]] = block.T

tr30B_2 = np.trace(D30_B @ D30_B)
tr30B_4 = np.trace(D30_B @ D30_B @ D30_B @ D30_B)
print(f"\n  Model B (full-rank coupling, identity blocks):")
print(f"  Tr(D^2) = {tr30B_2:.4f}")
print(f"  Tr(D^4) = {tr30B_4:.4f}")
print(f"  Ratio 2*Tr(D^4)/Tr(D^2) = {2*tr30B_4/tr30B_2:.6f}")

# Model C: dimension-weighted coupling
# Block (i,j) = w_ij * sqrt(d_i * d_j) / d_coupling
# This gives higher weight to larger representations
D30_C = np.zeros((dim_total, dim_total))
for idx, (i, j) in enumerate(edges):
    w = edge_weights[idx]
    # Each of d_min independent channels gets weight w
    # Total coupling strength = w * sqrt(min(d_i,d_j))
    d_min = min(dims[i], dims[j])
    block = np.zeros((dims[i], dims[j]))
    for k in range(d_min):
        block[k, k] = w / np.sqrt(d_min)  # normalize per channel
    D30_C[offsets[i]:offsets[i+1], offsets[j]:offsets[j+1]] = block
    D30_C[offsets[j]:offsets[j+1], offsets[i]:offsets[i+1]] = block.T

tr30C_2 = np.trace(D30_C @ D30_C)
tr30C_4 = np.trace(D30_C @ D30_C @ D30_C @ D30_C)
print(f"\n  Model C (normalized per channel):")
print(f"  Tr(D^2) = {tr30C_2:.4f}")
print(f"  Tr(D^4) = {tr30C_4:.4f}")
print(f"  Ratio 2*Tr(D^4)/Tr(D^2) = {2*tr30C_4/tr30C_2:.6f}")

# =====================================================================
# 7. CONNES FORMULA COMPARISON
# =====================================================================
print(f"\n{'=' * 72}")
print("7. CONNES FORMULA: m_H/m_W FROM VARIOUS D_F")
print(f"{'=' * 72}")

# Standard Connes: m_H^2/m_W^2 = 2*Tr(Y^4)/Tr(Y^2)
# where Y is Yukawa matrix (diagonal in standard model)

# Our D_F (weighted McKay adjacency) replaces Y
# The RATIO Tr(D^4)/Tr(D^2) is scale-independent

R_9x9 = tr4 / tr2
R_A = tr30A_4 / tr30A_2
R_B = tr30B_4 / tr30B_2
R_C = tr30C_4 / tr30C_2
R_unit = trA4 / trA2

target = phi**2  # tree-level
target_corr = (phi - 8*alpha)**2  # with correction

print(f"\n  Model               Tr(D^4)/Tr(D^2)  2*ratio    sqrt(2R)   phi    phi-8a")
for name, R in [("9x9 weighted", R_9x9),
                ("30x30 rank-1", R_A),
                ("30x30 full-rank", R_B),
                ("30x30 normalized", R_C),
                ("9x9 unit weight", R_unit)]:
    print(f"  {name:20s}  {R:12.4f}   {2*R:8.4f}  {np.sqrt(2*R):8.4f}  {phi:6.4f}  {phi-8*alpha:.4f}")

# =====================================================================
# 8. DIAGONAL YUKAWA (standard Connes for comparison)
# =====================================================================
print(f"\n{'=' * 72}")
print("8. STANDARD CONNES: DIAGONAL YUKAWA")
print(f"{'=' * 72}")

# Physical Yukawa couplings: y_f = sqrt(2) * m_f / v
# m_f = m_e * phi^n_f, v = 246 GeV
m_e = 0.51099895e-3  # GeV
v = 246.22  # GeV

# Color factors: Nc=3 for quarks, Nc=1 for leptons
Nc = {'e': 1, 'mu': 1, 'tau': 1, 'u': 3, 'd': 3, 'c': 3, 's': 3, 't': 3, 'b': 3}

a_std = 0  # Tr(Y^2) with color
b_std = 0  # Tr(Y^4) with color
for name, n in mass_exponents.items():
    y = np.sqrt(2) * m_e * phi**n / v
    nc = Nc[name]
    a_std += nc * y**2
    b_std += nc * y**4

print(f"  a = Tr(Y^2) = {a_std:.6e}")
print(f"  b = Tr(Y^4) = {b_std:.6e}")
print(f"  b/a = {b_std/a_std:.6e}")
print(f"  2*b/a = {2*b_std/a_std:.6e}")
print(f"  sqrt(2*b/a) = {np.sqrt(2*b_std/a_std):.6f}")
print(f"  This gives m_H/m_W = {np.sqrt(2*b_std/a_std):.4f} (standard Connes)")
print(f"  ~ sqrt(2)*y_t/g ~ sqrt(2)*m_t/m_W = {np.sqrt(2)*172.76/80.377:.4f}")

# =====================================================================
# 9. THE KEY: WHAT RATIO GIVES phi?
# =====================================================================
print(f"\n{'=' * 72}")
print("9. WHAT D_F GIVES phi?")
print(f"{'=' * 72}")

# We need 2*Tr(D^4)/Tr(D^2) = phi^2 = 2.618
# i.e., Tr(D^4)/Tr(D^2) = phi^2/2 = 1.309

# For the 9x9 unit adjacency: ratio = 3.0 (too high)
# For weighted: depends on weights

# What if D_F is not the adjacency but the LAPLACIAN?
L9 = np.diag(degrees.astype(float)) - A
L9_2 = L9 @ L9
L9_4 = L9_2 @ L9_2
trL2 = np.trace(L9_2)
trL4 = np.trace(L9_4)
R_lap = trL4 / trL2

print(f"  Laplacian: Tr(L^4)/Tr(L^2) = {R_lap:.6f}, sqrt(2R) = {np.sqrt(2*R_lap):.6f}")

# Normalized Laplacian
D_inv_sq = np.diag(1.0/np.sqrt(np.maximum(degrees, 1).astype(float)))
Ln = np.eye(9) - D_inv_sq @ A @ D_inv_sq
Ln2 = Ln @ Ln
Ln4 = Ln2 @ Ln2
R_nlap = np.trace(Ln4) / np.trace(Ln2)
print(f"  Norm. Laplacian: ratio = {R_nlap:.6f}, sqrt(2R) = {np.sqrt(2*R_nlap):.6f}")

# What about the Coxeter-weighted adjacency?
# Weight each edge (i,j) by sqrt(c_i * c_j) where c_k = dim_k (Coxeter label)
D_cox = np.zeros((9, 9))
for i, j in edges:
    w = np.sqrt(dims[i] * dims[j])
    D_cox[i, j] = w
    D_cox[j, i] = w

D_cox2 = D_cox @ D_cox
D_cox4 = D_cox2 @ D_cox2
R_cox = np.trace(D_cox4) / np.trace(D_cox2)
print(f"  Coxeter-weighted: ratio = {R_cox:.6f}, sqrt(2R) = {np.sqrt(2*R_cox):.6f}")

# What if we weight by Coxeter labels (dimensions)?
D_dim = np.zeros((9, 9))
for i, j in edges:
    D_dim[i, j] = dims[j]  # weight by target dimension
    D_dim[j, i] = dims[i]  # symmetric
# This is NOT symmetric! Let me use (dims[i]+dims[j])/2
D_dim2 = np.zeros((9, 9))
for i, j in edges:
    w = (dims[i] + dims[j]) / 2.0
    D_dim2[i, j] = w
    D_dim2[j, i] = w
D_d2_2 = D_dim2 @ D_dim2
D_d2_4 = D_d2_2 @ D_d2_2
R_dim2 = np.trace(D_d2_4) / np.trace(D_d2_2)
print(f"  Avg-dim weighted: ratio = {R_dim2:.6f}, sqrt(2R) = {np.sqrt(2*R_dim2):.6f}")

# =====================================================================
# 10. SCAN: WHAT SINGLE WEIGHT GIVES phi?
# =====================================================================
print(f"\n{'=' * 72}")
print("10. SCANNING FOR phi^2/2")
print(f"{'=' * 72}")

# For a tree with adjacency A, if we use D = w*A for scalar w:
# Tr(D^4)/Tr(D^2) = w^2 * Tr(A^4)/Tr(A^2) = w^2 * 3
# This is always 3*w^2, never phi^2/2 regardless of w.
# SO: uniform weighting cannot give phi.

# We need NON-UNIFORM weights.
# Target: Tr(D^4)/Tr(D^2) = phi^2/2 = 1.309017

target_ratio = phi**2 / 2
print(f"  Target ratio Tr(D^4)/Tr(D^2) = phi^2/2 = {target_ratio:.6f}")

# For a diagonal matrix D = diag(d_0, ..., d_8):
# Tr(D^2) = sum d_i^2, Tr(D^4) = sum d_i^4
# Ratio = sum d_i^4 / sum d_i^2
# This is >= (max d_i)^2 * min(d_i^2/d_i^2) ... hard to make small.
# Actually ratio >= max(d_i^2) for concentrated distributions.
# For equal d: ratio = d^2. So we need d = sqrt(phi^2/2) ~ 1.144.

# But D_F is OFF-diagonal (adjacency), not diagonal!
# For an off-diagonal matrix, the trace structure is different.

# Let me parametrize: D = sum_{edges} w_e * (|i><j| + |j><i|)
# Then: D^2 = sum_i (sum_{j~i} w_{ij}^2) |i><i| + cross terms at dist 2
# Tr(D^2) = 2 * sum_e w_e^2
# Tr(D^4) = sum_i (sum_{j~i} w_{ij}^2)^2 + sum_{dist-2 paths i-k-j} 2*w_{ik}^2*w_{kj}^2
# Wait, need to be more careful.

# For a tree: (D^2)_{ii} = sum_{j~i} w_{ij}^2
# (D^2)_{ij} for dist(i,j)=2: w_{ik}*w_{kj} (unique intermediate k)
# (D^4)_{ii} = sum_j (D^2)_{ij}^2
#            = (D^2)_{ii}^2 + sum_{dist-2 neighbors j} (D^2)_{ij}^2
#            = (sum_{j~i} w_{ij}^2)^2 + sum_{k~i} sum_{j~k,j!=i} (w_{ik}*w_{kj})^2

# So for the E8~ tree, Tr(D^4)/Tr(D^2) depends on the weight distribution.

# Let me try: w_e = phi^{-dist_from_higgs_edge} (exponential decay)
# Higgs edge = (0,1) (trivial to fundamental)

print(f"\n  Trying various weight schemes:")

def compute_ratio(weights_dict):
    """Given edge weights, compute Tr(D^4)/Tr(D^2)."""
    D = np.zeros((9, 9))
    for idx, (i, j) in enumerate(edges):
        w = weights_dict[idx]
        D[i, j] = w
        D[j, i] = w
    D2 = D @ D
    D4 = D2 @ D2
    t2 = np.trace(D2)
    t4 = np.trace(D4)
    if t2 == 0:
        return 0, 0, 0
    return t4/t2, t2, t4

schemes = {
    'unit': [1]*8,
    'distance from 0': [1/(d+1) for d in [0,1,2,3,4,5,6,5]],  # edge distances
    'phi^{-k}': [phi**(-k) for k in range(8)],
    'Coxeter prod': [np.sqrt(dims[i]*dims[j]) for i,j in edges],
    'dim min': [min(dims[i],dims[j]) for i,j in edges],
    'inverse dim': [1.0/max(dims[i],dims[j]) for i,j in edges],
}

print(f"  {'Scheme':25s} {'Tr4/Tr2':10s} {'2*ratio':10s} {'sqrt(2R)':10s}")
for name, weights in schemes.items():
    R, t2, t4 = compute_ratio(weights)
    print(f"  {name:25s} {R:10.4f} {2*R:10.4f} {np.sqrt(2*R):10.4f}")

# =====================================================================
# 11. HIGGS POTENTIAL FROM SPECTRAL ACTION
# =====================================================================
print(f"\n{'=' * 72}")
print("11. HIGGS POTENTIAL FROM SPECTRAL ACTION")
print(f"{'=' * 72}")

# In NCG, the Higgs potential arises from expanding Tr(f(D(h)^2))
# where D(h) = D + h*A_H, with A_H the Higgs direction.
# On the McKay graph, A_H is the coupling along the Higgs edge.

# Higgs edge: (0,1) connecting trivial to fundamental
# A_H = |0><1| + |1><0| (unit coupling along Higgs edge)

A_H = np.zeros((9, 9))
A_H[0, 1] = 1
A_H[1, 0] = 1

# D_F(h) = D_F_bg + h * A_H
# where D_F_bg is D_F with the Higgs edge removed (background)
D_bg = D9.copy()
D_bg[0, 1] = 0
D_bg[1, 0] = 0

# Expand in h: D(h) = D_bg + h*A_H
# D(h)^2 = D_bg^2 + h*(D_bg*A_H + A_H*D_bg) + h^2*A_H^2
# D(h)^4 = ... (expand to 4th order in h)

# Compute Tr(D(h)^2) and Tr(D(h)^4) as polynomials in h
# Use numerical differentiation for simplicity

hs = np.linspace(-5, 5, 1001)
tr2_h = []
tr4_h = []
for h in hs:
    Dh = D_bg + h * A_H
    Dh2 = Dh @ Dh
    Dh4 = Dh2 @ Dh2
    tr2_h.append(np.trace(Dh2))
    tr4_h.append(np.trace(Dh4))

tr2_h = np.array(tr2_h)
tr4_h = np.array(tr4_h)

# Fit polynomials
c2 = np.polyfit(hs, tr2_h, 4)  # Tr(D^2) is quadratic in h
c4 = np.polyfit(hs, tr4_h, 8)  # Tr(D^4) is quartic in h

print(f"  Tr(D(h)^2) polynomial coefficients (descending):")
for k, coef in enumerate(c2):
    if abs(coef) > 1e-6:
        print(f"    h^{len(c2)-1-k}: {coef:12.4f}")

print(f"\n  Tr(D(h)^4) polynomial coefficients (descending):")
for k, coef in enumerate(c4):
    if abs(coef) > 1e-6:
        print(f"    h^{len(c4)-1-k}: {coef:12.4f}")

# The Higgs potential: V(h) = -mu^2 * h^2 + lambda * h^4
# From spectral action: S ~ f_2*Tr(D^2)/Lambda^2 - f_0*Tr(D^4)/(2*Lambda^4) + ...
# V(h) = -f_2*[h^2 coef of Tr(D^2)]/Lambda^2 + f_0*[h^4 coef of Tr(D^4)]/(2*Lambda^4)

# Extract h^2 coefficient from Tr(D^2) and h^4 coefficient from Tr(D^4)
# Tr(D^2) = a0 + a2*h^2 (no linear or quartic terms for tree)
a2_coef = c2[-3]  # coefficient of h^2

# Tr(D^4) = b0 + b2*h^2 + b4*h^4
b2_coef = c4[-3] if len(c4) > 2 else 0  # h^2 coef
b4_coef = c4[-5] if len(c4) > 4 else 0  # h^4 coef

print(f"\n  Key coefficients:")
print(f"  Tr(D^2): h^2 coef = {a2_coef:.4f}")
print(f"  Tr(D^4): h^2 coef = {b2_coef:.4f}")
print(f"  Tr(D^4): h^4 coef = {b4_coef:.4f}")

# Higgs potential V(h) = mu^2*h^2 + lambda*h^4
# mu^2 = -(c_1/Lambda^2)*a2 + (c_2/(2*Lambda^4))*b2
# lambda = (c_2/(2*Lambda^4))*b4

# For SSB: mu^2 < 0 and lambda > 0
# VEV: v^2 = -mu^2/(2*lambda)
# m_H^2 = -2*mu^2 = 4*lambda*v^2
# m_W^2 = g^2*v^2/4

# The RATIO m_H^2/m_W^2 = 16*lambda/g^2 = 8*(-mu^2)/g^2/v^2 = ...
# This depends on the cutoff Lambda and spectral moments f_k.

# BUT: if we look at the ratio of the quartic to quadratic terms
# lambda/mu^2 = b4 / (some combination of a2 and b2)
# At large Lambda: mu^2 ~ c_1*a2/Lambda^2, lambda ~ c_2*b4/(2*Lambda^4)
# v^2 ~ c_1*a2*Lambda^2/(c_2*b4)
# m_H^2 ~ 2*c_1*a2/Lambda^2
# m_W^2 ~ g^2*c_1*a2*Lambda^2/(4*c_2*b4)
# m_H^2/m_W^2 ~ 8*c_2*b4/(g^2*c_1*Lambda^4) * Lambda^2
# ... depends on Lambda. Need self-consistency.

# The SCALE-INDEPENDENT prediction comes from:
# m_H^2/m_W^2 = 2*b4/a2 (if V comes purely from Tr(D^4))

if a2_coef != 0 and b4_coef != 0:
    R_potential = 2 * b4_coef / a2_coef
    print(f"\n  Scale-independent ratio 2*b4/a2 = {R_potential:.6f}")
    print(f"  sqrt(2*b4/a2) = {np.sqrt(abs(R_potential)):.6f}")
    print(f"  phi = {phi:.6f}")
    print(f"  phi - 8*alpha = {phi - 8*alpha:.6f}")

# =====================================================================
# SUMMARY
# =====================================================================
print(f"\n{'=' * 72}")
print("SUMMARY")
print(f"{'=' * 72}")
print(f"""
  The Connes formula m_H^2/m_W^2 = 2*Tr(D_F^4)/Tr(D_F^2) gives:

  Standard (diagonal Yukawa):  sqrt(2)*m_t/m_W ~ {np.sqrt(2)*172.76/80.377:.2f} (WRONG)
  Unit McKay adjacency:        sqrt(6) ~ {np.sqrt(6):.4f} (too high)
  Weighted McKay (Wilson):     sqrt(2R) = {np.sqrt(2*R_9x9):.4f}

  Target: phi = {phi:.4f} or phi-8*alpha = {phi-8*alpha:.4f}

  The McKay graph topology DOES change the ratio compared to
  the standard diagonal Yukawa, but the specific value depends
  on the edge weighting scheme.
""")

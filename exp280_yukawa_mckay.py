"""
EXP-280: YUKAWA MATRIX FROM McKAY GRAPH
=========================================
Open problem: Full off-diagonal D_F on McKay/E8 graph
Can we derive CKM mixing from D_F?

McKay node assignment (from exp264):
  0(e) - 1(u) - 2(d) - 3(s) - 4(c) - 5(b) - 6(t) - 7(tau)
                                 |
                                8(mu)
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
N_eig = 9; N_gen = 3
ALPHA = 1/137.035999084

print("="*72)
print("EXP-280: YUKAWA MATRIX FROM McKAY GRAPH")
print("="*72)

# === Part 1: E8 Dynkin diagram adjacency ===
print("\n--- Part 1: E8 Adjacency Matrix ---")

# E8 Dynkin diagram (standard Bourbaki numbering):
# 1 - 2 - 3 - 4 - 5 - 6 - 7
#                  |
#                  8
# But our McKay assignment maps fermions differently:
# Node: 0(e) - 1(u) - 2(d) - 3(s) - 4(c) - 5(b) - 6(t) - 7(tau)
#                                      |
#                                     8(mu)
# This is the AFFINE E8 (E8~) with the extra node 0

# Affine E8 adjacency (9 nodes, 0-indexed):
# Main chain: 0-1-2-3-4-5-6-7, branch: 4-8
adj_E8 = np.zeros((9, 9))
# Main chain
for i in range(7):
    adj_E8[i, i+1] = 1
    adj_E8[i+1, i] = 1
# Branch
adj_E8[4, 8] = 1
adj_E8[8, 4] = 1

labels = ['e', 'u', 'd', 's', 'c', 'b', 't', 'tau', 'mu']
n_exp = [0, 3, 5, 11, 16, 19, 26, 17, 11]  # mass exponents

print("  Affine E8 adjacency (fermion assignment):")
print("  " + "  ".join(f"{l:>4s}" for l in labels))
for i in range(9):
    row = "  ".join(f"{int(adj_E8[i,j]):4d}" for j in range(9))
    print(f"  {labels[i]:>4s}: {row}")

# === Part 2: D_F diagonal ===
print("\n--- Part 2: D_F Diagonal ---")
D_F = np.diag(n_exp)
print(f"  D_F = diag({n_exp})")
print(f"  Tr(D_F) = {sum(n_exp)} = sum of all exponents")
print(f"  Tr(D_F^2) = {sum(n**2 for n in n_exp)}")

# === Part 3: Off-diagonal D_F from E8 edges ===
print("\n--- Part 3: Off-Diagonal D_F ---")

# In NCG: D_F has off-diagonal entries = Yukawa couplings
# On a graph: natural choice is D_F_ij = f(n_i, n_j) * A_ij
# where A_ij is the adjacency matrix

# Attempt 1: D_F = diag(n) + w * A_E8
# where w is the off-diagonal coupling

# For CKM: the quark mass matrices are
# M_up = diag(m_u, m_c, m_t) in mass basis
# M_down = V_CKM * diag(m_d, m_s, m_b) * V_CKM^dag in up-mass basis

# The off-diagonal D_F entries generate mixing
# D_F_ij for quarks: connects u-d (W coupling), d-s (Cabibbo), etc.

# From E8 edges in our assignment:
# u-d: edge exists (nodes 1-2)
# d-s: edge exists (nodes 2-3)
# s-c: edge exists (nodes 3-4)
# c-b: edge exists (nodes 4-5)
# b-t: edge exists (nodes 5-6)
# u-s: NO edge (nodes 1,3)
# u-c: NO edge (nodes 1,4)
# etc.

print("  E8 edges between quarks:")
quarks_up = [('u', 1), ('c', 4), ('t', 6)]
quarks_down = [('d', 2), ('s', 3), ('b', 5)]

for qu_name, qu_idx in quarks_up:
    for qd_name, qd_idx in quarks_down:
        edge = int(adj_E8[qu_idx, qd_idx])
        dist = abs(qu_idx - qd_idx) if qd_idx != 5 else "branch"
        print(f"    {qu_name}-{qd_name}: edge={edge}, graph distance ~ {dist}")

# Graph distances on E8~
print("\n  Graph distances (up x down):")
# BFS for distances
def bfs_dist(adj, start):
    n = adj.shape[0]
    dist = [-1]*n
    dist[start] = 0
    queue = [start]
    while queue:
        node = queue.pop(0)
        for nb in range(n):
            if adj[node, nb] > 0 and dist[nb] == -1:
                dist[nb] = dist[node] + 1
                queue.append(nb)
    return dist

dists = {}
for i in range(9):
    dists[i] = bfs_dist(adj_E8, i)

print("      " + "  ".join(f"{labels[q[1]]:>4s}" for q in quarks_down))
for qu_name, qu_idx in quarks_up:
    row = "  ".join(f"{dists[qu_idx][qd_idx]:4d}" for _, qd_idx in quarks_down)
    print(f"  {qu_name:>4s}: {row}")

# === Part 4: CKM from graph distances ===
print("\n--- Part 4: CKM from Graph Distances ---")

# CKM |V_ij| ~ phi^(-distance_ij) heuristic
# Let's test this

import math

# Experimental CKM magnitudes
V_CKM_exp = np.array([
    [0.97373, 0.2243, 0.00382],  # V_ud, V_us, V_ub
    [0.221, 0.975, 0.0408],      # V_cd, V_cs, V_cb
    [0.0086, 0.0415, 0.99914]    # V_td, V_ts, V_tb
])

print("  |V_ij| experimental:")
up_labels = ['u', 'c', 't']
down_labels = ['d', 's', 'b']
for i, qu in enumerate(up_labels):
    row = "  ".join(f"{V_CKM_exp[i,j]:.5f}" for j in range(3))
    print(f"    {qu}: {row}")

print("\n  Graph distances (up x down):")
dist_matrix = np.zeros((3, 3), dtype=int)
for i, (_, qu_idx) in enumerate(quarks_up):
    for j, (_, qd_idx) in enumerate(quarks_down):
        dist_matrix[i, j] = dists[qu_idx][qd_idx]

for i, qu in enumerate(up_labels):
    row = "  ".join(f"{dist_matrix[i,j]:7d}" for j in range(3))
    print(f"    {qu}: {row}")

print("\n  phi^(-dist) prediction:")
V_pred = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        V_pred[i, j] = PHI**(-dist_matrix[i, j])

# Normalize rows to unit probability
for i in range(3):
    V_pred[i] /= np.sqrt(np.sum(V_pred[i]**2))

for i, qu in enumerate(up_labels):
    row = "  ".join(f"{V_pred[i,j]:.5f}" for j in range(3))
    print(f"    {qu}: {row}")

print("\n  Ratios |V_exp|/|V_pred|:")
for i, qu in enumerate(up_labels):
    row = "  ".join(f"{V_CKM_exp[i,j]/V_pred[i,j]:.4f}" if V_pred[i,j] > 0.001 else "   ---" for j in range(3))
    print(f"    {qu}: {row}")

# === Part 5: Froggatt-Nielsen on E8 ===
print("\n--- Part 5: Froggatt-Nielsen Charges from E8 ---")

# FN mechanism: |V_ij| ~ epsilon^(|q_i - q_j|) where epsilon ~ phi^(-1)
# q_i = FN charge = distance from some reference node

# Reference = junction node 4 (charm):
# q_u = dist(u, junction) = dist(1, 4) = 3
# q_c = dist(c, junction) = 0
# q_t = dist(t, junction) = 2
# q_d = dist(d, junction) = 2
# q_s = dist(s, junction) = 1
# q_b = dist(b, junction) = 1

junction = 4
print(f"  Junction node: {junction} ({labels[junction]})")
q_up = [dists[qi][junction] for _, qi in quarks_up]
q_down = [dists[qi][junction] for _, qi in quarks_down]
print(f"  FN charges (up):   {list(zip(up_labels, q_up))}")
print(f"  FN charges (down): {list(zip(down_labels, q_down))}")

print("\n  |V_ij| ~ phi^(-|q_up_i - q_down_j|):")
for i, qu in enumerate(up_labels):
    for j, qd in enumerate(down_labels):
        charge_diff = abs(q_up[i] - q_down[j])
        pred = PHI**(-charge_diff)
        exp_val = V_CKM_exp[i, j]
        print(f"    V_{qu}{qd}: |dq|={charge_diff}, phi^(-{charge_diff})={pred:.5f}, exp={exp_val:.5f}, ratio={exp_val/pred:.3f}")

# === Part 6: Off-diagonal mass matrix ===
print("\n--- Part 6: Constructing Mass Matrix ---")

# The up-type mass matrix in flavor basis:
# M_up_ij = phi^(n_i) * delta_ij + epsilon * phi^((n_i+n_j)/2) * A_ij
# where epsilon is the off-diagonal coupling from E8

# Up quarks: nodes 1(u), 4(c), 6(t), exponents 3, 16, 26
n_up = np.array([3, 16, 26], dtype=float)
n_down = np.array([5, 11, 19], dtype=float)

# E8 adjacency for up quarks:
# u(1)-c(4): distance 3, no direct edge
# u(1)-t(6): distance 5, no direct edge
# c(4)-t(6): distance 2, no direct edge
# Actually NO direct E8 edges connect up quarks to each other!
# Same for down quarks: d(2)-s(3) edge, s(3)-b(5) distance 2, d(2)-b(5) distance 3

print("  Up-quark E8 connections:")
for i, (n1, idx1) in enumerate(zip(up_labels, [1, 4, 6])):
    for j, (n2, idx2) in enumerate(zip(up_labels, [1, 4, 6])):
        if i < j:
            print(f"    {n1}-{n2}: edge={int(adj_E8[idx1,idx2])}, dist={dists[idx1][idx2]}")

print("\n  Down-quark E8 connections:")
for i, (n1, idx1) in enumerate(zip(down_labels, [2, 3, 5])):
    for j, (n2, idx2) in enumerate(zip(down_labels, [2, 3, 5])):
        if i < j:
            print(f"    {n1}-{n2}: edge={int(adj_E8[idx1,idx2])}, dist={dists[idx1][idx2]}")

# Down quarks: d(2)-s(3) have direct edge!
# So M_down has off-diagonal d-s entry

# Attempt: M_down = m_e * diag(phi^n) + epsilon * A_sub
# where A_sub is the E8 adjacency restricted to down quarks
A_down = np.zeros((3, 3))
# d(2)-s(3): edge
A_down[0, 1] = 1; A_down[1, 0] = 1
# s(3)-b(5): no direct edge (distance 2)
# d(2)-b(5): no direct edge (distance 3)

print(f"\n  A_down (E8 edges among d,s,b):")
print(f"  {A_down}")

# With phi^(-dist) off-diagonal:
A_down_dist = np.zeros((3, 3))
down_nodes = [2, 3, 5]
for i in range(3):
    for j in range(3):
        if i != j:
            d = dists[down_nodes[i]][down_nodes[j]]
            A_down_dist[i, j] = PHI**(-d)

print(f"\n  phi^(-dist) off-diagonal for down quarks:")
for i in range(3):
    row = "  ".join(f"{A_down_dist[i,j]:.5f}" for j in range(3))
    print(f"    {down_labels[i]}: {row}")

# Build full mass matrix: M = diag(phi^n) + A_dist
M_down_full = np.diag(PHI**n_down) + A_down_dist
print(f"\n  M_down = diag(phi^n) + phi^(-dist):")
for i in range(3):
    row = "  ".join(f"{M_down_full[i,j]:.5f}" for j in range(3))
    print(f"    {row}")

# Diagonalize
eigenvalues, V = np.linalg.eigh(M_down_full)
print(f"\n  Eigenvalues: {eigenvalues}")
print(f"  Expected phi^n: {PHI**5:.5f}, {PHI**11:.5f}, {PHI**19:.5f}")

# CKM-like mixing from V
print(f"\n  Mixing matrix |V|:")
for i in range(3):
    row = "  ".join(f"{abs(V[i,j]):.5f}" for j in range(3))
    print(f"    {row}")

# === Part 7: W-boson mediated mixing ===
print("\n--- Part 7: W-Mediated Mixing via E8 Paths ---")

# The CKM matrix comes from the MISMATCH between up and down mass bases
# On E8: up-type node i connects to down-type node j via paths
# V_ud: u(1)-d(2), distance 1, DIRECT edge -> V_ud ~ 1
# V_us: u(1)-s(3), distance 2, through d -> V_us ~ phi^(-2)
# V_ub: u(1)-b(5), distance 4, through d,s,c -> V_ub ~ phi^(-4)
# V_cd: c(4)-d(2), distance 2, through s -> V_cd ~ phi^(-2)
# V_cs: c(4)-s(3), distance 1, DIRECT edge -> V_cs ~ 1
# V_cb: c(4)-b(5), distance 1, DIRECT edge -> V_cb ~ phi^(-1)??
# V_td: t(6)-d(2), distance 4 -> V_td ~ phi^(-4)
# V_ts: t(6)-s(3), distance 3 -> V_ts ~ phi^(-3)
# V_tb: t(6)-b(5), distance 1, DIRECT edge -> V_tb ~ 1

# BUT: our CKM exponents are n_12=3, n_23=7, n_13=12
# Not matching simple distances

print("  CKM exponents vs E8 distances:")
ckm_pairs = [
    ('ud', 1, 2, 0, 0), ('us', 1, 3, 0, 1), ('ub', 1, 5, 0, 2),
    ('cd', 4, 2, 1, 0), ('cs', 4, 3, 1, 1), ('cb', 4, 5, 1, 2),
    ('td', 6, 2, 2, 0), ('ts', 6, 3, 2, 1), ('tb', 6, 5, 2, 2)
]

# Known CKM structure: theta_ij ~ phi^(-n_ij)
# n_12 = 3 (Cabibbo), n_23 = 7, n_13 = 12
n_ckm = {(0,1): 3, (1,2): 7, (0,2): 12}

for name, ni, nj, gi, gj in ckm_pairs:
    d = dists[ni][nj]
    exp_v = V_CKM_exp[gi, gj]
    pred_d = PHI**(-d)
    # Also check n_up - n_down difference
    n_diff = abs(n_up[gi] - n_down[gj])
    print(f"    V_{name}: dist={d}, phi^(-{d})={pred_d:.5f}, exp={exp_v:.5f}, |n_u-n_d|={n_diff}")

# === Part 8: The right formula ===
print("\n--- Part 8: CKM from Exponent Differences ---")
# theta_12 ~ phi^(-3): 3 = n_u = exponent of up quark
# theta_23 ~ phi^(-7): 7 = ??
# theta_13 ~ phi^(-12): 12 = ?? but 3+7 = 10 = 2*a1, and 3*12 = 36 = b1^2

# More precisely from paper:
# theta_C = arctan(phi^-3 * (1 - 2*alpha_s/(3*pi)))
# theta_23 = arctan(phi^-7 * (1 + 5*alpha_s/pi))
# theta_13 = arctan(phi^-12 * (1 + sin^2(tW)))

# The exponents 3, 7, 12 come from E8 leg dimensions:
# n_ij = N_eig * b_j - a1 * b_i - b1
# For (1,2): 9*1 - 5*0 - 6 = 3
# For (2,3): 9*2 - 5*1 - 6 = 7
# For (1,3): 9*2 - 5*0 - 6 = 12... wait that's same as (2,3)+something

# Let me check: n(b1,b2) = 9*b2 - 5*b1 - 6
# b1 = generation index of up quark (0,1,2)
# b2 = generation index of down quark (0,1,2)... no

# Actually from paper: n_12=3, n_23=7, n_13=12
# 3 = n_u exponent, 7 = n_tau - 2*a1 = 17-10, 12 = 2*b1 = 12
# Or: n_ij = N_eig * j - a1 * i - b1 where i,j are CKM indices
for i in range(3):
    for j in range(3):
        if i < j:
            n_ij = N_eig * j - a1 * i - b1
            print(f"    n_{i+1}{j+1} = 9*{j} - 5*{i} - 6 = {n_ij}")

# Wait, let me recalculate
N_eig = 9
print(f"\n  Formula n(i,j) = {N_eig}*j - {a1}*i - {b1}:")
for i in range(3):
    for j in range(3):
        n_ij = N_eig * j - a1 * i - b1
        print(f"    n({i},{j}) = {n_ij}")

# With indices starting from 0:
# n(0,1) = 9*1 - 0 - 6 = 3 (Cabibbo!)
# n(1,2) = 9*2 - 5 - 6 = 7 (correct!)
# n(0,2) = 9*2 - 0 - 6 = 12 (correct!)
# n(0,0) = -6, n(1,1) = 9-5-6 = -2, n(2,2) = 18-10-6 = 2
print(f"\n  CKM exponents: n_12={N_eig*1-a1*0-b1}, n_23={N_eig*2-a1*1-b1}, n_13={N_eig*2-a1*0-b1}")
print(f"  ALL CORRECT from n(i,j) = N_eig*j - a1*i - b1")

# How does this relate to E8 graph distance?
# n(i,j) = N_eig*j - a1*i - b1
# For the E8 assignment: up=(1,4,6), down=(2,3,5)
# dist(up_i, down_j) on E8:
print(f"\n  E8 distance vs CKM exponent:")
up_nodes = [1, 4, 6]
down_nodes_list = [2, 3, 5]
for i in range(3):
    for j in range(3):
        if i <= j and not (i == j):
            d_E8 = dists[up_nodes[i]][down_nodes_list[j]]
            n_ckm_val = N_eig * j - a1 * i - b1
            print(f"    ({i},{j}): E8 dist = {d_E8}, CKM n = {n_ckm_val}, ratio = {n_ckm_val/d_E8:.2f}")

# === Summary ===
print("\n" + "="*72)
print("SUMMARY")
print("="*72)
print(f"  1. D_F diagonal = {n_exp} on McKay E8~ (9 nodes)")
print(f"  2. Off-diagonal from E8 edges: only d-s have direct connection")
print(f"  3. phi^(-dist) heuristic: qualitative but not quantitative for CKM")
print(f"  4. CKM exponents n(i,j) = N_eig*j - a1*i - b1: EXACT formula")
print(f"     But this is ALGEBRAIC, not from D_F graph structure")
print(f"  5. Full off-diagonal D_F: NEEDS the Yukawa structure beyond adjacency")
print(f"  STATUS: D_F diagonal complete, off-diagonal requires more work")
print("="*72)

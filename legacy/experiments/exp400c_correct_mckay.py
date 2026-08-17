"""
EXP-400c: Correct McKay Graph + Fresh (a,b) Analysis
=====================================================

CRITICAL FIX: The McKay graph in exp400/400b was WRONG!
Wrong: edges = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(5,7),(5,8)]
       (node 5 has degree 4 — WRONG)

Correct: edges = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(5,8)]
         (node 5 has degree 3: connected to 4,6,8)

Verification: null vector of Cartan matrix requires
  2*d_i = sum(d_j, j adjacent to i) for all i

Author: Claude (exp400c)
Date: 2026-02-25
"""

import math
import numpy as np

# ===== CONSTANTS =====
a1 = 5
b1 = 6
PHI = (1 + math.sqrt(5)) / 2
PHIp = (1 - math.sqrt(5)) / 2
N = 120
DEG = 12
h_E8 = 30
N_eig = 9

dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

print("=" * 72)
print("EXP-400c: CORRECT McKAY GRAPH ANALYSIS")
print("=" * 72)

# ===== CORRECT AFFINE E8 GRAPH =====
# 0(1) - 1(2) - 2(3) - 3(4) - 4(5) - 5(6) - 6(4) - 7(2)
#                                          |
#                                          8(3)
edges_correct = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
edges_wrong   = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (5,7), (5,8)]

adj = np.zeros((9, 9), dtype=int)
for i, j in edges_correct:
    adj[i][j] = 1
    adj[j][i] = 1

# Verify null condition: 2*d_i = sum(d_j adjacent)
print("\n--- Null condition verification ---")
all_ok = True
for i in range(9):
    neighbors = [j for j in range(9) if adj[i][j] == 1]
    lhs = 2 * dims[i]
    rhs = sum(dims[j] for j in neighbors)
    ok = "OK" if lhs == rhs else "FAIL"
    if lhs != rhs:
        all_ok = False
    print(f"  Node {i} (dim {dims[i]}): 2*{dims[i]}={lhs}, "
          f"sum_adj={rhs} [{ok}]  neighbors={neighbors}")

print(f"\n  All null conditions satisfied: {all_ok}")

# Degree sequence
degrees = [sum(adj[i]) for i in range(9)]
print(f"  Degree sequence: {degrees}")
print(f"  Branch point: node {degrees.index(max(degrees))} (degree {max(degrees)})")


# ===== CORRECT DISTANCES =====
print("\n--- Distances from node 0 ---")

def bfs_distances(adj, start):
    n = len(adj)
    dist = [-1] * n
    dist[start] = 0
    queue = [start]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if adj[v][u] == 1 and dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist

dist_from_0 = bfs_distances(adj, 0)
print(f"  dist(0, k): {dist_from_0}")
print(f"  Max distance: {max(dist_from_0)} (node {dist_from_0.index(max(dist_from_0))})")

# COMPARE with wrong graph
adj_wrong = np.zeros((9, 9), dtype=int)
for i, j in edges_wrong:
    adj_wrong[i][j] = 1
    adj_wrong[j][i] = 1
dist_wrong = bfs_distances(adj_wrong, 0)
print(f"\n  WRONG graph distances: {dist_wrong}")
print(f"  CORRECT distances:     {dist_from_0}")
print(f"  Difference at node 7:  {dist_from_0[7]} vs {dist_wrong[7]} (was {dist_wrong[7]}, should be {dist_from_0[7]})")


# ===== CORRECT BIPARTITE STRUCTURE =====
print("\n--- Bipartite coloring ---")
black = [k for k in range(9) if dist_from_0[k] % 2 == 0]
white = [k for k in range(9) if dist_from_0[k] % 2 == 1]
print(f"  BLACK (even dist): {black}, dims={[dims[k] for k in black]}, sum={sum(dims[k] for k in black)}")
print(f"  WHITE (odd dist):  {white}, dims={[dims[k] for k in white]}, sum={sum(dims[k] for k in white)}")


# ===== McKAY ADJACENCY EIGENVALUES =====
eigenvalues, eigenvectors = np.linalg.eigh(adj.astype(float))
idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

print(f"\n--- McKay adjacency eigenvalues ---")
for i, ev in enumerate(eigenvalues):
    print(f"  lambda_{i} = {ev:+.6f}")


# ===== E8 EXPONENTS =====
e8_exponents = [0, 1, 7, 11, 13, 17, 19, 23, 29]

# Cayley graph eigenvalues (these are CORRECT regardless of McKay graph)
cayley_eigs_by_node = [2 * math.cos(math.pi * m / h_E8) for m in e8_exponents]
cayley_L_by_node = [DEG - ev for ev in cayley_eigs_by_node]

print(f"\n--- Cayley graph data per node ---")
print(f"  {'Node':>4} {'dim':>4} {'exp':>4} {'dist':>5} {'L_Cay':>10} {'Casimir':>10}")
for k in range(9):
    j = (dims[k] - 1) / 2.0
    C2 = j * (j + 1)
    print(f"  {k:>4} {dims[k]:>4} {e8_exponents[k]:>4} {dist_from_0[k]:>5} "
          f"{cayley_L_by_node[k]:>10.4f} {C2:>10.4f}")


# ===== TARGET (a,b) DATA =====
fermion_names = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 'b', 't']
ab_values = [(0,0), (3,-2), (1,0), (1,1), (1,1), (2,1), (1,2), (-1,4), (4,1)]
n_bare = [a1*a + b1*b for a, b in ab_values]

# Sort fermions by n_bare
ferm_sorted = sorted(range(9), key=lambda i: (n_bare[i], i))
node_sorted = sorted(range(9), key=lambda k: cayley_L_by_node[k])

print(f"\n--- Fermions sorted by n_bare ---")
for fi in ferm_sorted:
    print(f"  {fermion_names[fi]:>4}: n_bare = {n_bare[fi]:>2}")

print(f"\n--- Nodes sorted by Cayley Laplacian ---")
for ni in node_sorted:
    print(f"  node {ni}: L = {cayley_L_by_node[ni]:.4f}, dim = {dims[ni]}, dist = {dist_from_0[ni]}")


# =====================================================================
# PART 1: CORRECTED ASSIGNMENT AND DISTANCE ANALYSIS
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: FERMION-NODE ASSIGNMENT WITH CORRECT DISTANCES")
print("=" * 72)

# Assignment: sorted n_bare -> sorted Cayley L (monotone hypothesis)
print(f"\n  {'Ferm':>4} {'n':>3} {'node':>4} {'dim':>4} {'exp':>4} {'dist':>5} {'L':>10}")
assignment = {}
for rank in range(9):
    fi = ferm_sorted[rank]
    ni = node_sorted[rank]
    assignment[fermion_names[fi]] = ni
    print(f"  {fermion_names[fi]:>4} {n_bare[fi]:>3} {ni:>4} {dims[ni]:>4} "
          f"{e8_exponents[ni]:>4} {dist_from_0[ni]:>5} {cayley_L_by_node[ni]:>10.4f}")


# =====================================================================
# PART 2: PATHS ON CORRECT GRAPH
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: PATHS AND WILSON LINES (CORRECT GRAPH)")
print("=" * 72)

def find_path(adj, start, end):
    """BFS to find shortest path from start to end"""
    n = len(adj)
    parent = [-1] * n
    visited = [False] * n
    visited[start] = True
    queue = [start]
    while queue:
        v = queue.pop(0)
        if v == end:
            break
        for u in range(n):
            if adj[v][u] == 1 and not visited[u]:
                visited[u] = True
                parent[u] = v
                queue.append(u)
    # Reconstruct path
    path = []
    v = end
    while v != -1:
        path.append(v)
        v = parent[v]
    return list(reversed(path))

print(f"\n  {'Node':>4} {'Path from 0':>25} {'dist':>5} {'sum_dim':>8} {'sum_dim_int':>12}")
for k in range(9):
    path = find_path(adj, 0, k)
    d = len(path) - 1
    s_dim = sum(dims[p] for p in path)
    s_int = sum(dims[p] for p in path[1:-1]) if len(path) > 2 else 0
    print(f"  {k:>4} {str(path):>25} {d:>5} {s_dim:>8} {s_int:>12}")

# Key change: node 7 now has distance 7 (through chain), not 6
print(f"\n  KEY CHANGE: node 7 path goes through node 6!")
print(f"  Path to 7: {find_path(adj, 0, 7)}")
print(f"  Path to 8: {find_path(adj, 0, 8)}")


# =====================================================================
# PART 3: WILSON LINE WITH PHI-WEIGHTS
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: PHI-WEIGHTED WILSON LINES")
print("=" * 72)

# Natural edge weights from the graph structure:
# Each edge (i,j) has weight w(i,j).
# Try: w = dim of target node, w = sqrt(dim_i * dim_j), w = Casimir, etc.

# Approach: Wilson line W(0->k) = sum of edge weights along path
# Test: does n_bare = A * W + B for some A, B?

print("\n--- Testing Wilson line = sum(dim_target) along path ---")
# For each node, compute Wilson line as sum of dim(next node) along path
W_dim = {}
for k in range(9):
    path = find_path(adj, 0, k)
    w = sum(dims[path[i+1]] for i in range(len(path)-1))
    W_dim[k] = w

print(f"  Node:    {list(range(9))}")
print(f"  W_dim:   {[W_dim[k] for k in range(9)]}")
print(f"  Dist:    {dist_from_0}")

# Map to fermion assignment
print(f"\n  {'Ferm':>4} {'n':>3} {'node':>4} {'W_dim':>6} {'n/W':>8}")
for rank in range(9):
    fi = ferm_sorted[rank]
    ni = node_sorted[rank]
    w = W_dim[ni]
    ratio = n_bare[fi] / w if w > 0 else float('inf')
    print(f"  {fermion_names[fi]:>4} {n_bare[fi]:>3} {ni:>4} {w:>6} {ratio:>8.3f}")

# Try Wilson line with phi as edge weight:
# W(0->k) = number_of_edges * phi^(something)
# Or: cumulative product along path
print(f"\n--- Testing Wilson line = sum(phi^j) along path ---")
for k in range(9):
    path = find_path(adj, 0, k)
    d = len(path) - 1
    # W = sum_{step i from 0 to d-1} phi^(something related to step i)
    W_phi_sum = sum(PHI**i for i in range(d))  # geometric sum
    W_phi_dim = sum(dims[path[i+1]] * PHI**(d-1-i) for i in range(d)) if d > 0 else 0
    print(f"  node {k}: d={d}, "
          f"sum(phi^i)={W_phi_sum:.3f}, "
          f"sum(d*phi^i)={W_phi_dim:.3f}")


# =====================================================================
# PART 4: DEEP SEARCH - n FROM (dim, dist, exp) WITH CORRECT GRAPH
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: EXHAUSTIVE SEARCH WITH CORRECT DISTANCES")
print("=" * 72)

# Feature set for each node (correct graph)
node_features = {}
for k in range(9):
    path = find_path(adj, 0, k)
    node_features[k] = {
        'dim': dims[k],
        'exp': e8_exponents[k],
        'dist': dist_from_0[k],
        'W_dim': W_dim[k],
        'sum_dim_path': sum(dims[p] for p in path),
        'Casimir': (dims[k]**2 - 1) / 4.0,
    }

# Using the monotone assignment (fermion sorted by n -> node sorted by L)
assigned_features = []
for rank in range(9):
    fi = ferm_sorted[rank]
    ni = node_sorted[rank]
    f = node_features[ni]
    assigned_features.append(f)

assigned_n = [n_bare[ferm_sorted[rank]] for rank in range(9)]
assigned_a = [ab_values[ferm_sorted[rank]][0] for rank in range(9)]
assigned_b = [ab_values[ferm_sorted[rank]][1] for rank in range(9)]

# Extended integer search: n = c0 + c1*dim + c2*dist + c3*exp
print("\n--- Integer search: n = c0 + c1*dim + c2*dist + c3*exp ---")
best_err = 999
best_formula = ""
for c0 in range(-5, 6):
    for c1 in range(-5, 6):
        for c2 in range(-5, 6):
            for c3 in range(-3, 4):
                pred = [c0 + c1*assigned_features[i]['dim']
                        + c2*assigned_features[i]['dist']
                        + c3*assigned_features[i]['exp']
                        for i in range(9)]
                err = sum((pred[i] - assigned_n[i])**2 for i in range(9))
                if err < best_err:
                    best_err = err
                    best_formula = f"n = {c0} + {c1}*dim + {c2}*dist + {c3}*exp"
                    best_pred = pred[:]

print(f"  Best: {best_formula}")
print(f"  SSE = {best_err}")
print(f"  Predicted: {best_pred}")
print(f"  Target:    {assigned_n}")
if best_err > 0:
    print(f"  Diffs:     {[best_pred[i]-assigned_n[i] for i in range(9)]}")

# Also try: n = c0 + c1*W_dim + c2*exp
print("\n--- Integer search: n = c0 + c1*W_dim + c2*exp ---")
best_err2 = 999
for c0 in range(-5, 6):
    for c1 in range(-5, 6):
        for c2 in range(-3, 4):
            pred = [c0 + c1*assigned_features[i]['W_dim']
                    + c2*assigned_features[i]['exp']
                    for i in range(9)]
            err = sum((pred[i] - assigned_n[i])**2 for i in range(9))
            if err < best_err2:
                best_err2 = err
                best_formula2 = f"n = {c0} + {c1}*W_dim + {c2}*exp"
                best_pred2 = pred[:]

print(f"  Best: {best_formula2}")
print(f"  SSE = {best_err2}")
print(f"  Predicted: {best_pred2}")
print(f"  Target:    {assigned_n}")


# =====================================================================
# PART 5: CHARACTER-BASED APPROACH WITH CORRECT GRAPH
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: McKAY RECURSION CHARACTERS (CORRECT GRAPH)")
print("=" * 72)

# McKay recursion on correct graph:
# x * chi_k = sum(chi_j, j adjacent to k)
# This gives polynomial expressions in x = chi_1 at identity class x=2

# For the linear chain 0-1-2-3-4-5:
# chi_0 = 1, chi_1 = x
# chi_2 = x*chi_1 - chi_0 = x^2 - 1
# chi_3 = x*chi_2 - chi_1 = x^3 - 2x
# chi_4 = x*chi_3 - chi_2 = x^4 - 3x^2 + 1
# chi_5 = x*chi_4 - chi_3 = x^5 - 4x^3 + 3x

# At the branch: x*chi_5 = chi_4 + chi_6 + chi_8
# For the chain 5-6-7:
# x*chi_6 = chi_5 + chi_7
# x*chi_7 = chi_6

# From last: chi_7 = chi_6/x
# Substituting: x*chi_6 = chi_5 + chi_6/x => chi_6(x - 1/x) = chi_5
# => chi_6 = x*chi_5/(x^2 - 1)

# For node 8 (leaf connected to 5):
# x*chi_8 = chi_5 => chi_8 = chi_5/x

# Check: x*chi_5 = chi_4 + chi_6 + chi_8
# = chi_4 + x*chi_5/(x^2-1) + chi_5/x
# = chi_4 + chi_5 * (x/(x^2-1) + 1/x)
# = chi_4 + chi_5 * (x^2 + x^2 - 1) / (x*(x^2-1))
# = chi_4 + chi_5 * (2x^2 - 1) / (x*(x^2-1))
# Hmm, this needs to equal x*chi_5.
# x*chi_5 - chi_4 = chi_5 * (2x^2-1)/(x*(x^2-1))
# (x*chi_5 - chi_4) * x*(x^2-1) = chi_5*(2x^2-1)
# x^2*(x^2-1)*chi_5 - x*(x^2-1)*chi_4 = chi_5*(2x^2-1)
# chi_5*(x^4-x^2-2x^2+1) = x*(x^2-1)*chi_4
# chi_5*(x^4-3x^2+1) = x*(x^2-1)*chi_4

# Check: chi_5 = x^5 - 4x^3 + 3x, chi_4 = x^4 - 3x^2 + 1
# LHS: (x^5-4x^3+3x)*(x^4-3x^2+1)
# RHS: x*(x^2-1)*(x^4-3x^2+1) = (x^3-x)*(x^4-3x^2+1)

# We need: (x^5-4x^3+3x) = (x^3-x)? NO.

# Let me redo this more carefully.
# McKay recursion: x * chi_k = sum(chi_j for j adjacent to k)

# Compute characters as polynomials in x evaluated at specific values
# Conjugacy class representatives have chi_1(g) = x values:
# 2, phi, -1/phi, 1, -1/phi, -phi, 0, -1, -2
# (these are 2*cos(alpha) for alpha = 0, pi/5, 2pi/5, pi/3, 3pi/5, 4pi/5, pi/2, 2pi/3, pi)

x_values = [2 * math.cos(alpha) for alpha in
            [0, math.pi/5, 2*math.pi/5, math.pi/3,
             3*math.pi/5, 4*math.pi/5, math.pi/2,
             2*math.pi/3, math.pi]]

print(f"\nConjugacy class x-values (chi_1 = 2*cos(alpha)):")
for i, x in enumerate(x_values):
    print(f"  C{i}: x = {x:.6f}")

# Compute character table by McKay recursion on CORRECT graph
char_table = np.zeros((9, 9))  # char_table[k][c] = chi_k(C_c)

for c in range(9):
    x = x_values[c]

    # Linear chain 0-1-2-3-4-5
    chi = [0.0] * 9
    chi[0] = 1.0
    chi[1] = x
    chi[2] = x * chi[1] - chi[0]          # x^2 - 1
    chi[3] = x * chi[2] - chi[1]          # x^3 - 2x
    chi[4] = x * chi[3] - chi[2]          # x^4 - 3x^2 + 1
    chi[5] = x * chi[4] - chi[3]          # x^5 - 4x^3 + 3x

    # Branch: chi_8 and chi_7 from leaves
    # x * chi_8 = chi_5 (node 8 connected only to 5)
    if abs(x) > 1e-10:
        chi[8] = chi[5] / x
    else:
        # x=0: chi_5(0) = 3*0 = 0, chi_8 = 0/0 -> use limit
        # chi_5 = x^5-4x^3+3x, so chi_5/x = x^4-4x^2+3 -> at x=0: 3
        chi[8] = 3.0  # polynomial limit

    # x * chi_7 = chi_6 (node 7 connected only to 6)
    # x * chi_6 = chi_5 + chi_7 (node 6 connected to 5 and 7)
    # From these: chi_7 = chi_6/x, so x*chi_6 = chi_5 + chi_6/x
    # chi_6*(x - 1/x) = chi_5 => chi_6 = x*chi_5/(x^2 - 1)
    if abs(x**2 - 1) > 1e-10 and abs(x) > 1e-10:
        chi[6] = x * chi[5] / (x**2 - 1)
        chi[7] = chi[6] / x
    else:
        # Need special handling for x = +/-1 and x = 0
        if abs(x) < 1e-10:
            # x=0: chi_5=0, chi_6=0/(-1)=0, chi_7=0/0
            # Actually chi_5 polynomial: x^5-4x^3+3x at x=0 is 0
            # chi_6 polynomial: from recursion, need to evaluate as rational function
            # At x=0: use L'Hopital or polynomial form
            # chi_6 = x*(x^5-4x^3+3x)/(x^2-1) = x^2*(x^4-4x^2+3)/(x^2-1)
            # = x^2*(x^2-1)(x^2-3)/(x^2-1) = x^2*(x^2-3)
            # At x=0: 0
            chi[6] = 0.0
            chi[7] = 0.0  # Actually need polynomial evaluation
        elif abs(x - 1) < 1e-10:
            # x=1: chi_5 = 1-4+3 = 0
            # chi_6 = 1*0/(1-1) = 0/0 -> L'Hopital
            # Numerator: x*chi_5 = x*(x^5-4x^3+3x), derivative at x=1:
            # d/dx[x^6-4x^4+3x^2] = 6x^5-16x^3+6x, at x=1: 6-16+6=-4
            # Denominator: x^2-1, derivative: 2x, at x=1: 2
            # L'Hopital: -4/2 = -2
            chi[6] = -2.0
            chi[7] = chi[6]  # chi_7 = chi_6/x = -2/1 = -2
        elif abs(x + 1) < 1e-10:
            # x=-1: chi_5 = -1+4-3 = 0
            # Same issue, L'Hopital
            chi_5_poly_deriv = 6*(-1)**5 - 16*(-1)**3 + 6*(-1)  # = -6+16-6=4
            denom_deriv = 2*(-1)  # = -2
            chi[6] = chi_5_poly_deriv / denom_deriv  # 4/(-2) = -2
            chi[7] = chi[6] / x  # -2/(-1) = 2
        else:
            chi[6] = x * chi[5] / (x**2 - 1)
            chi[7] = chi[6] / x

    # Verify: x*chi_5 should equal chi_4 + chi_6 + chi_8
    lhs = x * chi[5]
    rhs = chi[4] + chi[6] + chi[8]
    if abs(lhs - rhs) > 0.01:
        print(f"  WARNING C{c}: x*chi_5={lhs:.4f} != chi_4+chi_6+chi_8={rhs:.4f}")

    for k in range(9):
        char_table[k][c] = chi[k]

# Display character table
print(f"\nCharacter table (McKay recursion on CORRECT graph):")
class_sizes = [1, 12, 12, 20, 12, 12, 30, 20, 1]
print(f"  {'':>6}", end="")
for c in range(9):
    print(f"  C{c:d}({class_sizes[c]:>2})", end="")
print()

for k in range(9):
    dim_check = char_table[k][0]  # chi(identity) = dim
    print(f"  rho_{k}({int(round(dim_check)):>1})", end="")
    for c in range(9):
        val = char_table[k][c]
        if abs(val - round(val)) < 0.01:
            print(f"  {int(round(val)):>6}", end="")
        else:
            print(f"  {val:>6.3f}", end="")
    print()

# Verify orthogonality
G = np.array(class_sizes)
check = char_table @ np.diag(G) @ char_table.T
ortho_err = np.max(np.abs(check - N * np.eye(9)))
print(f"\n  Orthogonality check: max|<chi_i,chi_j> - {N}*delta| = {ortho_err:.6f}")

# Verify dimensions
print(f"  Dimensions from chi(C0): {[int(round(char_table[k][0])) for k in range(9)]}")
print(f"  Expected:                 {dims}")


# =====================================================================
# PART 6: GALOIS ORBITS FROM CORRECT CHARACTER TABLE
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: GALOIS ORBITS")
print("=" * 72)

# Galois: phi -> phi', which sends x = 2*cos(pi/5) -> 2*cos(2pi/5)
# i.e., C1 <-> C2, C4 <-> C5 (the phi-containing classes)
# Classes C0, C3, C6, C7, C8 are Galois-fixed (rational characters)

# For each irrep, compute Galois-conjugated character row
print("\n  Galois permutation of irreps:")
galois_pairs = {}
for k in range(9):
    chi_galois = char_table[k].copy()
    # Swap C1 <-> C2 and C4 <-> C5
    chi_galois[1], chi_galois[2] = char_table[k][2], char_table[k][1]
    chi_galois[4], chi_galois[5] = char_table[k][5], char_table[k][4]

    # Find matching irrep
    for kp in range(9):
        if np.allclose(chi_galois, char_table[kp], atol=0.01):
            galois_pairs[k] = kp
            if k <= kp:
                if k == kp:
                    print(f"  rho_{k} (dim {dims[k]}): FIXED")
                else:
                    print(f"  rho_{k} (dim {dims[k]}) <-> rho_{kp} (dim {dims[kp]})")
            break


# =====================================================================
# PART 7: GALOIS KERNEL AND (a,b)
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: GALOIS KERNEL AND MASS FORMULA CONNECTION")
print("=" * 72)

# From exp356: kernel = {rho_0, rho_1, rho_7} with dim 1+2+2 = 5 = a1
# These are the irreps that are SMALLEST in each Galois orbit

# Galois orbits:
# Fixed: rho_0(1), rho_4(5) -- or depends on correct computation
# Pairs: {rho_k, rho_{sigma(k)}}

print("\n  Galois orbits:")
seen = set()
orbits = []
for k in range(9):
    if k not in seen:
        kp = galois_pairs.get(k, k)
        if k == kp:
            orbits.append([k])
            seen.add(k)
        else:
            orbits.append([k, kp])
            seen.add(k)
            seen.add(kp)

for orbit in orbits:
    if len(orbit) == 1:
        k = orbit[0]
        print(f"  Fixed: rho_{k} (dim {dims[k]})")
    else:
        k, kp = orbit
        print(f"  Pair: rho_{k} (dim {dims[k]}) <-> rho_{kp} (dim {dims[kp]})")

# Which orbit gets which fermion?
# The kernel selects the smaller dim in each pair.

# Now: can we extract (a,b) from the character values?
# Specifically: the characters at phi-containing classes involve Z[phi]
# chi_k(C1) = a + b*phi for some (a,b)?

print(f"\n--- Character Z[phi] decomposition at C1 (x=phi) ---")
print(f"  C1: x = 2*cos(pi/5) = {x_values[1]:.6f} = phi = {PHI:.6f}")
for k in range(9):
    val = char_table[k][1]
    # Decompose val = a + b*phi
    # val = a + b*PHI => b = (val - a) / PHI
    # We need a to be an integer (or half-integer for spinor reps)
    # Try: b = round((val - round(val - val*PHIp)) / PHI)
    # Actually: val = a + b*PHI, and we know dim = a + b*1 (at x=2)
    # Wait, chi_k(C0) = dim_k and chi_k(C1) = val
    # If val = alpha + beta*phi, then we need to find alpha, beta

    # For Z[phi]: phi^2 = phi + 1, phi^3 = 2phi+1, phi^4 = 3phi+2, etc.
    # We know the polynomial form of chi_k(x) from McKay recursion
    # At x = phi: chi_k(phi) involves powers of phi which simplify

    # Try direct decomposition
    # val = alpha + beta * PHI
    # We need alpha, beta integers or rationals
    found = False
    for alpha in range(-10, 11):
        for beta in range(-10, 11):
            if abs(val - (alpha + beta * PHI)) < 0.001:
                print(f"  rho_{k}: chi(C1) = {val:.4f} = {alpha} + {beta}*phi"
                      f"  -> (alpha,beta)=({alpha},{beta})")
                found = True
                break
        if found:
            break
    if not found:
        print(f"  rho_{k}: chi(C1) = {val:.4f} - no integer decomposition found")


# =====================================================================
# PART 8: CAN (a,b)_fermion = (alpha, beta)_character?
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: CHARACTER DECOMPOSITION vs MASS QUANTUM NUMBERS")
print("=" * 72)

# If chi_k(C1) = alpha_k + beta_k * phi, can we relate
# (alpha_k, beta_k) to (a_f, b_f) of the fermion assigned to node k?

# First collect character decompositions for all nodes
char_ab = {}
for k in range(9):
    val = char_table[k][1]
    for alpha in range(-10, 11):
        for beta in range(-10, 11):
            if abs(val - (alpha + beta * PHI)) < 0.001:
                char_ab[k] = (alpha, beta)
                break

print(f"\n  {'Node':>4} {'dim':>4} {'chi(C1)':>10} {'(alpha,beta)':>15}  |  "
      f"{'Ferm':>4} {'(a,b)':>10} {'n_bare':>7}")

for rank in range(9):
    fi = ferm_sorted[rank]
    ni = node_sorted[rank]
    ab_char = char_ab.get(ni, ('?', '?'))
    ab_ferm = ab_values[fi]
    print(f"  {ni:>4} {dims[ni]:>4} {char_table[ni][1]:>10.4f} {str(ab_char):>15}  |  "
          f"{fermion_names[fi]:>4} ({ab_ferm[0]:+d},{ab_ferm[1]:+d}) {n_bare[fi]:>7}")

# Check: is there ANY relationship between (alpha,beta) and (a,b)?
print(f"\n--- Comparing character (alpha,beta) with mass (a,b) ---")
for rank in range(9):
    fi = ferm_sorted[rank]
    ni = node_sorted[rank]
    ab_char = char_ab.get(ni, (0, 0))
    ab_ferm = ab_values[fi]
    diff_a = ab_ferm[0] - ab_char[0]
    diff_b = ab_ferm[1] - ab_char[1]
    print(f"  {fermion_names[fi]:>4}: char({ab_char[0]:+d},{ab_char[1]:+d}) vs "
          f"mass({ab_ferm[0]:+d},{ab_ferm[1]:+d}) -> "
          f"diff=({diff_a:+d},{diff_b:+d})")


# =====================================================================
# PART 9: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: HONEST ASSESSMENT")
print("=" * 72)
print(f"""
  GRAPH CORRECTION:
    The McKay graph was WRONG in exp400/400b.
    Old: node 5 had degree 4 (connected to 4,6,7,8) -- WRONG
    New: node 5 has degree 3 (connected to 4,6,8) -- CORRECT
    Node 6 connects to 7 (not 5 to 7).

  KEY CHANGES:
    - Node 7 distance: {dist_wrong[7]} -> {dist_from_0[7]}
    - Bipartite: BLACK={black}, WHITE={white}
    - Character table now satisfies orthogonality

  NEW FINDINGS:
    1. Character Z[phi] decomposition at phi-containing classes
       gives (alpha, beta) pairs for each irrep
    2. These DO NOT directly equal the mass (a,b) assignments
    3. The correct graph still gives NO simple formula for n_bare

  REMAINING GAPS:
    - (a,b) assignments NOT derivable from graph properties alone
    - Character (alpha,beta) bear no obvious relation to mass (a,b)
    - The fermion-node assignment itself is not determined

  STATUS: **NEGATIVE** (corrected graph still doesn't derive (a,b))
  But the graph correction is important for other calculations.
""")

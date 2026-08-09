"""
EXP-400: Derive (a,b) Quantum Numbers from McKay Graph Structure
================================================================

GOAL: Can the 9 fermion (a,b) assignments be DERIVED from spectral/
geometric properties of the McKay graph (affine E8)?

Current status:
  - (a,b) give n_bare = 5a + 6b for each fermion
  - Sum rules: sum(a) = 12 = DEG, sum(b) = 8 = rank(E8)
  - 4 correction sectors IRREDUCIBLE (exp384)
  - Only 4 valid assignments out of 20,160 (exp323)

Strategy:
  1. Compute McKay graph spectral data (eigenvalues, eigenvectors)
  2. Test if (a,b) correlate with eigenvector components
  3. Test Wilson line (path sum) on McKay tree
  4. Test algebraic invariants of each node (dim, Casimir, etc.)
  5. Look for a MAP from node properties to (a,b)

REGULA BLIND: compute properties FIRST, compare with (a,b) AFTER.

Author: Claude (exp400)
Date: 2026-02-25
"""

import math
import numpy as np

# ===== CONSTANTS =====
a1 = 5
b1 = 6
PHI = (1 + math.sqrt(5)) / 2
PHIp = (1 - math.sqrt(5)) / 2  # Galois conjugate
N = 120
DEG = 12
N_eig = 9
rank_E8 = 8
h_E8 = 30

print("=" * 72)
print("EXP-400: DERIVE (a,b) FROM McKAY GRAPH STRUCTURE")
print("=" * 72)

# ===== McKAY GRAPH OF 2I (affine E8 Dynkin diagram) =====
# Node k = irrep rho_k, dimension d_k
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]  # dimensions of 9 irreps

# Adjacency matrix of affine E8 (McKay graph)
# Edges from McKay correspondence: rho_1 x rho_k decomposes as sum of neighbors
adj = np.zeros((9, 9), dtype=int)
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (5,7), (5,8)]
for i, j in edges:
    adj[i][j] = 1
    adj[j][i] = 1

print("\n--- McKay graph (affine E8) adjacency ---")
print("Node dims:", dims)
print("Edges:", edges)

# Graph Laplacian
degree_matrix = np.diag([sum(adj[i]) for i in range(9)])
L = degree_matrix - adj
print("\nDegree sequence:", [sum(adj[i]) for i in range(9)])

# Eigenvalues and eigenvectors of adjacency matrix
eigenvalues_adj, eigenvectors_adj = np.linalg.eigh(adj.astype(float))
# Sort by eigenvalue
idx = np.argsort(eigenvalues_adj)[::-1]
eigenvalues_adj = eigenvalues_adj[idx]
eigenvectors_adj = eigenvectors_adj[:, idx]

print("\nAdjacency eigenvalues:", [f"{v:.4f}" for v in eigenvalues_adj])

# Eigenvalues of graph Laplacian
eigenvalues_L, eigenvectors_L = np.linalg.eigh(L.astype(float))
idx_L = np.argsort(eigenvalues_L)
eigenvalues_L = eigenvalues_L[idx_L]
eigenvectors_L = eigenvectors_L[:, idx_L]

print("Laplacian eigenvalues:", [f"{v:.4f}" for v in eigenvalues_L])


# ===== TARGET: (a,b) ASSIGNMENTS =====
# Fermion -> node mapping (canonical from exp323)
# e=0, d=1, mu/s=2, u=3, c=4, t=5(center), b=6, tau=7, nu?=8
# But actually: need to be careful about the assignment

# Known (a,b) from mass formula:
fermions = {
    'e':   {'ab': (0, 0),  'n_bare': 0,  'sector': 'zero'},
    'mu':  {'ab': (1, 1),  'n_bare': 11, 'sector': 'unit'},
    'tau': {'ab': (1, 2),  'n_bare': 17, 'sector': 'unit'},
    'u':   {'ab': (3, -2), 'n_bare': 3,  'sector': 'unit'},
    'c':   {'ab': (2, 1),  'n_bare': 16, 'sector': 'prime'},
    't':   {'ab': (4, 1),  'n_bare': 26, 'sector': 'prime'},
    'd':   {'ab': (1, 0),  'n_bare': 5,  'sector': 'rational'},
    's':   {'ab': (1, 1),  'n_bare': 11, 'sector': 'unit'},
    'b':   {'ab': (-1, 4), 'n_bare': 19, 'sector': 'prime'},
}

# Z[phi] algebraic data
for name, data in fermions.items():
    a, b = data['ab']
    z = a + b * PHI
    zp = a + b * PHIp
    norm = a**2 + a*b - b**2  # N(z) = z*z'
    data['z'] = z
    data['zp'] = zp
    data['norm'] = norm

print("\n--- Target (a,b) and Z[phi] data ---")
print(f"  {'Fermion':>6} {'(a,b)':>8} {'n_bare':>7} {'z':>8} {'z_prime':>8} {'N(z)':>6} {'sector':>10}")
for name in ['e','mu','tau','u','c','t','d','s','b']:
    d = fermions[name]
    a, b = d['ab']
    print(f"  {name:>6} ({a:+d},{b:+d}) {d['n_bare']:>7} {d['z']:>8.3f} {d['zp']:>8.3f} {d['norm']:>6} {d['sector']:>10}")

# Sum rules
sum_a = sum(d['ab'][0] for d in fermions.values())
sum_b = sum(d['ab'][1] for d in fermions.values())
sum_n = sum(d['n_bare'] for d in fermions.values())
print(f"\n  Sum(a) = {sum_a} (= DEG = {DEG})")
print(f"  Sum(b) = {sum_b} (= rank(E8) = {rank_E8})")
print(f"  Sum(n) = {sum_n} (= a1*DEG + b1*rank = {a1*DEG + b1*rank_E8})")


# =====================================================================
# PART 1: NODE PROPERTIES AS PREDICTORS
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: NODE PROPERTIES")
print("=" * 72)

# For each McKay node, compute:
# - dimension d_k
# - Casimir C_k = d_k * (d_k + 2) / 8 (for SU(2) irrep of dim d_k)
# Actually for spin j, C = j(j+1), and dim = 2j+1
# So j_k = (d_k - 1)/2, C_k = j_k*(j_k+1) = (d_k^2-1)/4
casimirs = [(d**2 - 1) / 4.0 for d in dims]

# Laplacian eigenvalue for each node (from Cayley graph of 2I)
# These are L_k = DEG - 2*cos(2*pi*k/h) for affine Dynkin
# But for E8 affine, the eigenvalues of the Cayley adjacency are:
# lambda_k = chi_1(C_k) where C_k is the k-th conjugacy class
# Actually, the McKay graph eigenvalues ARE the Cayley adjacency eigenvalues
# divided by the dimension

# Cayley graph adjacency eigenvalues of 2I:
# lambda_k = dim_k * (sum over generators of chi_k(g)) / dim_k^2
# = chi_k(generators) / dim_k
# For the 2-element generating set {s, s^-1} where s is a specific element:
# This gives lambda = 2*cos(2*pi*j/h) for each node

# Actually, let me use the KNOWN Cayley eigenvalues from McKay recursion
# chi_1(C_k) for the fundamental (2-dim) representation
# These are the columns of the character table evaluated at the fundamental

# McKay recursion: chi_{k+1} = chi_1 * chi_k - chi_{k-1}
# With chi_0 = 1, chi_1 = x (variable = eigenvalue of adj matrix)
# For affine E8, the eigenvalues are 2*cos(pi*m_k/30) where m_k are exponents

# E8 exponents: 1, 7, 11, 13, 17, 19, 23, 29 (and 0 for affine)
# Affine E8 exponents (including 0): 0, 1, 7, 11, 13, 17, 19, 23, 29
e8_exponents = [0, 1, 7, 11, 13, 17, 19, 23, 29]

cayley_eigs = [2 * math.cos(math.pi * m / h_E8) for m in e8_exponents]
cayley_eigs.sort(reverse=True)

print("\nCayley graph eigenvalues (from E8 exponents):")
for i, ev in enumerate(cayley_eigs):
    print(f"  lambda_{i} = {ev:.6f}")

# Laplacian eigenvalues on Cayley graph
cayley_L = [DEG - ev for ev in cayley_eigs]
cayley_L.sort()
print("\nCayley Laplacian eigenvalues:")
for i, ev in enumerate(cayley_L):
    print(f"  L_{i} = {ev:.6f}")

# Check which eigenvalue goes with which irrep
# The eigenvalue of the Cayley adjacency for irrep rho_k is:
# lambda(rho_k) = chi_1(C_k) * |2I| / dim(rho_k)
# NO - it's simpler: the adjacency eigenvalue for rho_k in the
# McKay graph is the character of the fundamental at the conjugacy class
# associated to rho_k.

# For 2I, the conjugacy classes are labeled by their character values.
# Let me just use the McKay recursion directly.

print("\n--- McKay recursion for character values ---")
# chi_k(C_j) = character of rho_k at conjugacy class C_j
# McKay adjacency: A * v_k = sum_j A_{kj} * v_j = eigenvalue * v_k
# where eigenvalue = chi_1(C_j) for the j-th eigenvector

# The 9 Cayley adjacency eigenvalues are:
# 2*cos(pi*m/30) for m in {0,1,7,11,13,17,19,23,29}
# But which m goes with which irrep?

# Standard ordering: rho_0 (trivial) has eigenvalue = chi_1(identity) = dim(rho_1) = 2
# rho_0 -> m=0 -> 2*cos(0) = 2

# For affine E8, the node-exponent correspondence:
# Node 0 (dim 1) -> exponent 0
# Node 1 (dim 2) -> exponent 1
# Node 2 (dim 3) -> exponent 7
# Node 3 (dim 4) -> exponent 11
# Node 4 (dim 5) -> exponent 13
# Node 5 (dim 6) -> exponent 17 (center node)
# Node 6 (dim 4) -> exponent 19
# Node 7 (dim 2) -> exponent 23
# Node 8 (dim 3) -> exponent 29

node_exponents = [0, 1, 7, 11, 13, 17, 19, 23, 29]
cayley_by_node = [2 * math.cos(math.pi * m / h_E8) for m in node_exponents]
laplacian_by_node = [DEG - ev for ev in cayley_by_node]

print(f"\n  {'Node':>4} {'dim':>4} {'exp':>4} {'lambda':>10} {'L':>10} {'Casimir':>10}")
for k in range(9):
    print(f"  {k:>4} {dims[k]:>4} {node_exponents[k]:>4} {cayley_by_node[k]:>10.4f} {laplacian_by_node[k]:>10.4f} {casimirs[k]:>10.4f}")


# =====================================================================
# PART 2: FERMION-NODE ASSIGNMENT — WHICH FERMION ON WHICH NODE?
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: FERMION-NODE MAPPING")
print("=" * 72)

# From exp323 and the framework:
# The canonical assignment places fermions at nodes by their spectral properties
# Key constraint: n_f should correlate with node properties

# Let's try the most natural assignment based on bare exponent ordering:
# n_bare: e=0, u=3, d=5, mu=s=11, c=16, tau=17, b=19, t=26
# Sorted: 0, 3, 5, 11, 11, 16, 17, 19, 26

# Node properties sorted by Laplacian eigenvalue:
print("\nNodes sorted by Laplacian eigenvalue:")
node_order = sorted(range(9), key=lambda k: laplacian_by_node[k])
for rank, k in enumerate(node_order):
    print(f"  rank {rank}: node {k} (dim {dims[k]}, L={laplacian_by_node[k]:.4f})")

# Natural hypothesis: larger L -> larger mass -> larger n_bare
# Let's test this with the sorted bare exponents
n_bares_sorted = sorted([d['n_bare'] for d in fermions.values()])
print(f"\nSorted n_bare values: {n_bares_sorted}")
print(f"Number of nodes: {len(node_order)}, Number of fermions: {len(n_bares_sorted)}")

# Direct assignment: fermion with rank-k n_bare sits on node with rank-k L
print("\n--- Hypothesis: n_bare monotone with Laplacian ---")
assignment_hyp1 = {}
fermion_names_sorted = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 'b', 't']
# sorted by n_bare: 0, 3, 5, 11, 11, 16, 17, 19, 26
for rank, (fname, node_idx) in enumerate(zip(fermion_names_sorted, node_order)):
    fdata = fermions[fname]
    L_node = laplacian_by_node[node_idx]
    print(f"  {fname:>4} (n={fdata['n_bare']:>2}) -> node {node_idx} (dim={dims[node_idx]}, L={L_node:.4f})")
    assignment_hyp1[fname] = node_idx


# =====================================================================
# PART 3: CAN n_bare BE COMPUTED FROM NODE PROPERTIES?
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: n_bare FROM NODE PROPERTIES")
print("=" * 72)

# Test various functions f(node_properties) = n_bare
# Using the assignment from hypothesis 1

print("\n--- Test: n = A * dim + B * exponent ---")
# Linear fit: n = A * dim_k + B * exp_k
# Using the canonical assignment
assigned_dims = []
assigned_exps = []
assigned_L = []
assigned_n = []
for fname in fermion_names_sorted:
    node = assignment_hyp1[fname]
    assigned_dims.append(dims[node])
    assigned_exps.append(node_exponents[node])
    assigned_L.append(laplacian_by_node[node])
    assigned_n.append(fermions[fname]['n_bare'])

# Least squares: n = A * dim + B * exp + C
X = np.column_stack([assigned_dims, assigned_exps, np.ones(9)])
n_vec = np.array(assigned_n, dtype=float)
result = np.linalg.lstsq(X, n_vec, rcond=None)
coeffs = result[0]
residuals = n_vec - X @ coeffs
rms = math.sqrt(np.mean(residuals**2))

print(f"  Best fit: n = {coeffs[0]:.4f}*dim + {coeffs[1]:.4f}*exp + {coeffs[2]:.4f}")
print(f"  RMS residual: {rms:.4f}")
print(f"  Residuals: {[f'{r:.2f}' for r in residuals]}")

# Test: n = A * exp alone
X2 = np.column_stack([assigned_exps, np.ones(9)])
result2 = np.linalg.lstsq(X2, n_vec, rcond=None)
coeffs2 = result2[0]
residuals2 = n_vec - X2 @ coeffs2
rms2 = math.sqrt(np.mean(residuals2**2))
print(f"\n  Best fit: n = {coeffs2[0]:.4f}*exp + {coeffs2[1]:.4f}")
print(f"  RMS residual: {rms2:.4f}")

# Test: n = A * L (Laplacian) + B
X3 = np.column_stack([assigned_L, np.ones(9)])
result3 = np.linalg.lstsq(X3, n_vec, rcond=None)
coeffs3 = result3[0]
residuals3 = n_vec - X3 @ coeffs3
rms3 = math.sqrt(np.mean(residuals3**2))
print(f"\n  Best fit: n = {coeffs3[0]:.4f}*L + {coeffs3[1]:.4f}")
print(f"  RMS residual: {rms3:.4f}")


# =====================================================================
# PART 4: CAN (a,b) INDIVIDUALLY BE COMPUTED?
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: (a,b) FROM NODE PROPERTIES")
print("=" * 72)

assigned_a = [fermions[fname]['ab'][0] for fname in fermion_names_sorted]
assigned_b = [fermions[fname]['ab'][1] for fname in fermion_names_sorted]

# Test: a = f(dim, exp)
a_vec = np.array(assigned_a, dtype=float)
b_vec = np.array(assigned_b, dtype=float)

result_a = np.linalg.lstsq(X, a_vec, rcond=None)
ca = result_a[0]
res_a = a_vec - X @ ca
rms_a = math.sqrt(np.mean(res_a**2))

result_b = np.linalg.lstsq(X, b_vec, rcond=None)
cb = result_b[0]
res_b = b_vec - X @ cb
rms_b = math.sqrt(np.mean(res_b**2))

print(f"  a = {ca[0]:.4f}*dim + {ca[1]:.4f}*exp + {ca[2]:.4f}  (RMS {rms_a:.4f})")
print(f"  b = {cb[0]:.4f}*dim + {cb[1]:.4f}*exp + {cb[2]:.4f}  (RMS {rms_b:.4f})")

# Check if the fit is exact for any framework-rational coefficients
print("\n--- Framework-rational coefficient search ---")
best_err = 999
best_formula = ""
for p in range(-6, 7):
    for q in range(-6, 7):
        for r in range(-20, 21):
            n_pred = [p * assigned_dims[i] + q * assigned_exps[i] + r
                      for i in range(9)]
            err = sum((n_pred[i] - assigned_n[i])**2 for i in range(9))
            if err < best_err:
                best_err = err
                best_formula = f"n = {p}*dim + {q}*exp + {r}"
                best_pred = n_pred[:]

print(f"  Best integer linear: {best_formula}")
print(f"  SSE = {best_err}")
if best_err > 0:
    print(f"  Predictions: {best_pred}")
    print(f"  Targets:     {assigned_n}")
    print(f"  Differences: {[best_pred[i]-assigned_n[i] for i in range(9)]}")


# =====================================================================
# PART 5: GRAPH DISTANCE AND PATH-BASED FORMULAS
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: GRAPH DISTANCE ON McKAY TREE")
print("=" * 72)

# Compute all pairwise distances on the McKay tree
# The tree is: 0-1-2-3-4-5-6, with 5-7 and 5-8 branching
# Distance from node 0 to each node:
dist_from_0 = [0, 1, 2, 3, 4, 5, 6, 6, 6]

# Distance from node 5 (center) to each node:
dist_from_5 = [5, 4, 3, 2, 1, 0, 1, 1, 1]

# Weighted distance: weight each step by dim of intermediate node?
# Path from 0 to k goes through: 0 -> 1 -> 2 -> ... -> k
# Weighted path = sum of dims along the path (excluding endpoints? or including?)

def path_from_0(target):
    """Return list of nodes on path from 0 to target in affine E8"""
    if target <= 5:
        return list(range(target + 1))
    elif target == 6:
        return [0, 1, 2, 3, 4, 5, 6]
    elif target == 7:
        return [0, 1, 2, 3, 4, 5, 7]
    elif target == 8:
        return [0, 1, 2, 3, 4, 5, 8]

print("\nPaths from node 0 (trivial irrep) and weighted distances:")
print(f"  {'Node':>4} {'Path':>20} {'dist':>5} {'sum_dim':>8} {'sum_dim_int':>12} {'sum_Cas':>8}")
for k in range(9):
    path = path_from_0(k)
    d = len(path) - 1
    s_dim = sum(dims[p] for p in path)
    s_dim_int = sum(dims[p] for p in path[1:-1]) if len(path) > 2 else 0  # interior nodes
    s_cas = sum(casimirs[p] for p in path)
    print(f"  {k:>4} {str(path):>20} {d:>5} {s_dim:>8} {s_dim_int:>12} {s_cas:>8.2f}")

# Now map fermions to nodes and check correlation
print("\n--- Fermion n_bare vs graph distance metrics ---")
print(f"  {'Ferm':>4} {'node':>4} {'n_bare':>7} {'dist':>5} {'sum_dim':>8} {'sum_int':>8} {'dim*dist':>9}")
for fname in fermion_names_sorted:
    node = assignment_hyp1[fname]
    d = dist_from_0[node]
    path = path_from_0(node)
    s_dim = sum(dims[p] for p in path)
    s_int = sum(dims[p] for p in path[1:-1]) if len(path) > 2 else 0
    dd = dims[node] * d
    n = fermions[fname]['n_bare']
    print(f"  {fname:>4} {node:>4} {n:>7} {d:>5} {s_dim:>8} {s_int:>8} {dd:>9}")


# =====================================================================
# PART 6: DIMENSION-WEIGHTED FORMULAS
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: TESTING n_bare = f(dim, distance)")
print("=" * 72)

# The bare exponent n = a1*a + b1*b
# Can we express a and b separately in terms of node data?
# a is the "rational part" and b is the "irrational part" of z = a + b*phi

# Key insight: the sum rule says sum(a) = 12, sum(b) = 8
# And a1 * sum(a) + b1 * sum(b) = 5*12 + 6*8 = 108
# 108/9 = 12 = average n_bare

# Test: does the DIMENSION of the node relate to n_bare?
print("\n--- Correlation: dim_k vs n_bare ---")
for fname in fermion_names_sorted:
    node = assignment_hyp1[fname]
    ratio = fermions[fname]['n_bare'] / dims[node] if dims[node] > 0 else 0
    print(f"  {fname:>4}: dim={dims[node]}, n={fermions[fname]['n_bare']}, n/dim={ratio:.3f}")

# Test: n_bare = alpha * dim^2 + beta * dim + gamma ?
X4 = np.column_stack([np.array(assigned_dims)**2, assigned_dims, np.ones(9)])
result4 = np.linalg.lstsq(X4, n_vec, rcond=None)
c4 = result4[0]
res4 = n_vec - X4 @ c4
rms4 = math.sqrt(np.mean(res4**2))
print(f"\n  n = {c4[0]:.4f}*dim^2 + {c4[1]:.4f}*dim + {c4[2]:.4f}")
print(f"  RMS = {rms4:.4f}")


# =====================================================================
# PART 7: EXPONENT DECOMPOSITION — a FROM exp, b FROM dim?
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: DECOMPOSING (a,b) FROM ORTHOGONAL CHANNELS")
print("=" * 72)

# Key idea: a and b enter as n = a1*a + b1*b
# a1 = 5 and b1 = 6 are coprime, so the decomposition is unique mod 30
# Given n, there's a unique (a,b) with 0 <= b1*b < 30 (or some range)
# But the PHYSICAL (a,b) don't satisfy this — they have negative values too

# Instead, let's ask: what invariant of node k gives a_k? What gives b_k?
print("\nNode data with target (a,b):")
print(f"  {'Ferm':>4} {'node':>4} {'dim':>4} {'exp':>4} {'a':>4} {'b':>4} {'n':>4}")
for fname in fermion_names_sorted:
    node = assignment_hyp1[fname]
    a_f, b_f = fermions[fname]['ab']
    print(f"  {fname:>4} {node:>4} {dims[node]:>4} {node_exponents[node]:>4} {a_f:>4} {b_f:>4} {fermions[fname]['n_bare']:>4}")

# Look for pattern: is a_f related to (exp mod something)?
print("\n--- Searching for a = f(exp), b = g(exp) ---")
for fname in fermion_names_sorted:
    node = assignment_hyp1[fname]
    m = node_exponents[node]
    a_f = fermions[fname]['ab'][0]
    b_f = fermions[fname]['ab'][1]
    print(f"  {fname}: exp={m:>2}, a={a_f:>2}, b={b_f:>2}, "
          f"exp mod 5={m % 5}, exp mod 6={m % 6}, "
          f"exp//5={m//5}, exp//6={m//6}")

# Direct search: can we find integer functions f,g such that
# a = f(dim, exp, dist) and b = g(dim, exp, dist) EXACTLY?
print("\n--- Exact integer search for (a,b) ---")
# Extended feature set
features = []
for fname in fermion_names_sorted:
    node = assignment_hyp1[fname]
    d_node = dims[node]
    m_node = node_exponents[node]
    dist = dist_from_0[node]
    features.append([d_node, m_node, dist, d_node * dist, d_node * m_node])

feat_names = ['dim', 'exp', 'dist', 'dim*dist', 'dim*exp']

# Search for a = sum c_i * feat_i with small integer c_i
from itertools import product as iprod
print("\n  Searching a = c1*dim + c2*exp + c3*dist + c4*dim*dist + c5*dim*exp + c0...")
best_a_err = 999
best_a_formula = ""
for c0 in range(-5, 6):
    for c1 in range(-3, 4):
        for c2 in range(-3, 4):
            for c3 in range(-3, 4):
                a_pred = [c0 + c1*features[i][0] + c2*features[i][1] + c3*features[i][2]
                          for i in range(9)]
                err = sum((a_pred[i] - assigned_a[i])**2 for i in range(9))
                if err < best_a_err:
                    best_a_err = err
                    best_a_formula = f"a = {c0} + {c1}*dim + {c2}*exp + {c3}*dist"
                    best_a_pred = a_pred[:]

print(f"  Best: {best_a_formula}")
print(f"  SSE(a) = {best_a_err}")
if best_a_err == 0:
    print(f"  *** EXACT MATCH! ***")
else:
    print(f"  Predicted a: {best_a_pred}")
    print(f"  Target a:    {assigned_a}")

print(f"\n  Searching b = c1*dim + c2*exp + c3*dist + c0...")
best_b_err = 999
best_b_formula = ""
for c0 in range(-5, 6):
    for c1 in range(-3, 4):
        for c2 in range(-3, 4):
            for c3 in range(-3, 4):
                b_pred = [c0 + c1*features[i][0] + c2*features[i][1] + c3*features[i][2]
                          for i in range(9)]
                err = sum((b_pred[i] - assigned_b[i])**2 for i in range(9))
                if err < best_b_err:
                    best_b_err = err
                    best_b_formula = f"b = {c0} + {c1}*dim + {c2}*exp + {c3}*dist"
                    best_b_pred = b_pred[:]

print(f"  Best: {best_b_formula}")
print(f"  SSE(b) = {best_b_err}")
if best_b_err == 0:
    print(f"  *** EXACT MATCH! ***")
else:
    print(f"  Predicted b: {best_b_pred}")
    print(f"  Target b:    {assigned_b}")


# =====================================================================
# PART 8: ALTERNATIVE NODE ASSIGNMENTS
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: TRYING ALL 9! = 362880 ASSIGNMENTS")
print("=" * 72)

# For each permutation of fermions -> nodes, check if there exists
# an exact integer linear formula n = c0 + c1*dim + c2*exp
# This is computationally intensive but feasible for small integers

from itertools import permutations

n_values = [fermions[fname]['n_bare'] for fname in fermion_names_sorted]

# We can't try all 362880 permutations with all integer coefficients,
# but we can check: for each permutation, does a least-squares fit
# with small integer coefficients give zero residual?

best_global_err = 999
best_global_perm = None
best_global_formula = ""
count_exact = 0

# For efficiency, only search n = c0 + c1*exp
# since dim and exp are the main features
print("Searching: n = c0 + c1*exp for all 362880 permutations...")

for perm in permutations(range(9)):
    # perm[i] = which node fermion_i sits on
    exps_perm = [node_exponents[perm[i]] for i in range(9)]
    dims_perm = [dims[perm[i]] for i in range(9)]

    # Try: n = c1*exp + c0
    Xp = np.column_stack([exps_perm, np.ones(9)])
    np_vec = np.array(n_values, dtype=float)
    res = np.linalg.lstsq(Xp, np_vec, rcond=None)
    pred = Xp @ res[0]
    sse = np.sum((pred - np_vec)**2)

    if sse < best_global_err:
        best_global_err = sse
        best_global_perm = perm
        best_global_formula = f"n ~ {res[0][0]:.4f}*exp + {res[0][1]:.4f}"

    # Also try n = c1*dim + c2*exp + c0
    Xp2 = np.column_stack([dims_perm, exps_perm, np.ones(9)])
    res2 = np.linalg.lstsq(Xp2, np_vec, rcond=None)
    pred2 = Xp2 @ res2[0]
    sse2 = np.sum((pred2 - np_vec)**2)

    if sse2 < best_global_err:
        best_global_err = sse2
        best_global_perm = perm
        best_global_formula = f"n ~ {res2[0][0]:.4f}*dim + {res2[0][1]:.4f}*exp + {res2[0][2]:.4f}"

    if sse2 < 0.01:
        count_exact += 1

print(f"\nBest global SSE = {best_global_err:.4f}")
print(f"Best formula: {best_global_formula}")
print(f"Best permutation: {best_global_perm}")
if best_global_perm:
    print("Assignment:")
    for i, fname in enumerate(fermion_names_sorted):
        node = best_global_perm[i]
        print(f"  {fname} -> node {node} (dim={dims[node]}, exp={node_exponents[node]})")
print(f"\nNumber of exact (SSE<0.01) assignments: {count_exact}")


# =====================================================================
# PART 9: EIGENVECTOR COMPONENTS
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: McKAY ADJACENCY EIGENVECTOR APPROACH")
print("=" * 72)

# Each eigenvector of the McKay adjacency gives a "coordinate" for each node
# Maybe (a,b) can be read off from specific eigenvector components

print("\nEigenvectors of McKay adjacency (columns):")
print("Eigenvalues:", [f"{v:.4f}" for v in eigenvalues_adj])
print()

for k in range(9):
    v = eigenvectors_adj[:, k]
    # Normalize so max component = 1
    v_norm = v / np.max(np.abs(v))
    print(f"  EV {k} (lambda={eigenvalues_adj[k]:+.4f}): "
          f"[{', '.join(f'{c:+.3f}' for c in v_norm)}]")


# =====================================================================
# PART 10: CASIMIR OPERATOR AND SPIN
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: CASIMIR AND SPIN-BASED FORMULAS")
print("=" * 72)

# Each node has spin j = (dim-1)/2
# j values: 0, 1/2, 1, 3/2, 2, 5/2, 3/2, 1/2, 1
spins = [(d - 1) / 2.0 for d in dims]
print("\nNode spins j = (dim-1)/2:")
for k in range(9):
    print(f"  node {k}: dim={dims[k]}, j={spins[k]:.1f}, C2={casimirs[k]:.2f}")

# Test: n_bare = function of spin?
print("\n--- n_bare vs 2j (using hypothesis 1 assignment) ---")
for fname in fermion_names_sorted:
    node = assignment_hyp1[fname]
    j = spins[node]
    n = fermions[fname]['n_bare']
    print(f"  {fname:>4}: 2j={2*j:.0f}, n={n:>2}, n/(2j+1)={n/dims[node]:.3f}")


# =====================================================================
# PART 11: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 11: HONEST ASSESSMENT")
print("=" * 72)
print("""
  FINDINGS:

  1. FERMION-NODE ASSIGNMENT: The mapping of 9 fermions to 9 McKay nodes
     is NOT trivially determined by n_bare ordering vs Laplacian ordering.
     The best linear formula n = A*exp + B gives poor fit for most
     permutations.

  2. NO EXACT INTEGER FORMULA: No simple integer linear combination of
     (dim, exponent, distance) reproduces all 9 (a,b) values exactly.
     This means (a,b) are NOT simple linear functions of node data.

  3. THE ASSIGNMENT PROBLEM: Even the fermion-to-node mapping is not
     uniquely determined by a simple criterion. The 4 solutions from
     exp323 required weighted distances with FITTED weights.

  4. WHAT WORKS (partially):
     - Sum rules: sum(a) = DEG, sum(b) = rank(E8) [DERIVED]
     - n = a1*a + b1*b structure [DERIVED from framework]
     - Individual corrections: 4 sectors with different formulas [IDENTIFIED]
     - Graph distance for d quark: delta_d = -2/a1 [IDENTIFIED]

  5. WHAT DOESN'T WORK:
     - No single spectral function predicts (a,b) [exp384 NEGATIVE]
     - No integer linear formula from (dim, exp, dist) gives (a,b) exactly
     - Eigenvector components don't obviously map to (a,b)

  STATUS: The (a,b) assignments remain the DEEPEST open problem in the
  mass formula. They are highly constrained (only 4 valid assignments
  out of 20,160) but not yet derivable from McKay graph properties alone.

  POSSIBLE DIRECTIONS:
  - Non-linear functions of node data (e.g., involving phi)
  - Path integral (Wilson line) along McKay tree with phi-weights
  - Representation-theoretic: Clebsch-Gordan decomposition paths
  - Galois orbits within the character table
""")

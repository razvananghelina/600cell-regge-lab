"""
EXP-323c: Alternative Metrics on CORRECT McKay Graph
=====================================================

Three investigations on the correct affine E8 McKay graph (legs 1,2,5):

1. SIGNED WEIGHTS: Check if negative weights are on Galois-crossing edges
2. EFFECTIVE RESISTANCE: R_eff(rho_1, rho_i) with Cayley-weighted Laplacian
3. GREEN FUNCTION: G(rho_1, rho_i) = <rho_1|(L+eps)^-1|rho_i>
4. HEAT KERNEL: K_t(rho_1, rho_i) = <rho_1|e^{-tL}|rho_i>
5. NONLINEAR PATH: product phi^{w_e}, or sum w_e * dim(intermediate)

Author: Computational exploration for R.-C. Anghelina's framework
"""

import numpy as np
from itertools import combinations, permutations
from collections import deque

# ============================================================
# Constants
# ============================================================
phi = (1 + np.sqrt(5)) / 2
phi_prime = 1 - phi
a1 = 5
b1 = 6
N_gen = 3
N_eig = 9
rank_E8 = 8
SQRT5 = np.sqrt(5)

# Irrep dimensions
DIMS = np.array([1, 2, 2, 3, 3, 4, 4, 5, 6])
# 0=rho1(1), 1=rho2(2), 2=rho3(2'), 3=rho4(3), 4=rho5(3'),
# 5=rho6(4), 6=rho7(4'), 7=rho8(5), 8=rho9(6)

TARGET = np.array([0, 3, 5, 11, 11, 16, 17, 19, 26])
N_TO_F = {0:'e', 3:'u', 5:'d', 11:'mu/s', 16:'c', 17:'tau', 19:'b', 26:'t'}

# CORRECT McKay graph edges
EDGES = [(0,1), (1,3), (3,6), (6,7), (7,8), (8,4), (8,5), (5,2)]
N_EDGES = 8

# Cayley Laplacian eigenvalues per irrep (DERIVED from character table)
chi_10a = {0: 1, 1: phi, 2: -1/phi, 3: phi, 4: -1/phi,
           5: -1, 6: 1, 7: 0, 8: -1}
LAP_EIGS = {}
for i in range(9):
    LAP_EIGS[i] = 12 - 12 * chi_10a[i] / DIMS[i]

# Galois classification of irreps
# phi-sector (Sym^k of fund): rho2(1), rho4(3)
# phi'-sector (Sym^k of Galois fund): rho3(2), rho5(4)
# self-conjugate: rho1(0), rho6(5), rho7(6), rho8(7), rho9(8)
GALOIS_SECTOR = {0: 'self', 1: 'phi', 2: 'phi_prime', 3: 'phi', 4: 'phi_prime',
                 5: 'self', 6: 'self', 7: 'self', 8: 'self'}
GALOIS_PAIRS = [(1, 2), (3, 4)]  # (rho2,rho3), (rho4,rho5)

print("=" * 72)
print("EXP-323c: ALTERNATIVE METRICS ON CORRECT McKAY GRAPH")
print("=" * 72)

# ============================================================
# PART 0: Graph structure verification
# ============================================================
print("\n--- CORRECT McKay graph (affine E8, legs 1,2,5) ---")
adj = {}
for a, b in EDGES:
    adj.setdefault(a, []).append(b)
    adj.setdefault(b, []).append(a)

for i in range(9):
    deg = len(adj.get(i, []))
    sec = GALOIS_SECTOR[i]
    role = " [BRANCH]" if deg == 3 else (" [ENDPOINT]" if deg == 1 else "")
    print(f"  rho_{i+1}(d={DIMS[i]}, {sec:>10}): deg {deg}{role}, lambda={LAP_EIGS[i]:.4f}")

# BFS from node 0
parent = {0: None}
parent_edge = {0: None}
queue = deque([0])
visited = {0}
while queue:
    v = queue.popleft()
    for e_idx, (a, b) in enumerate(EDGES):
        for u in [a, b]:
            other = b if u == a else a
            if other == v and u not in visited:
                visited.add(u)
                parent[u] = v
                parent_edge[u] = e_idx
                queue.append(u)

# Path matrix
remaining = list(range(1, 9))
P = np.zeros((8, 8))
for idx, node in enumerate(remaining):
    v = node
    while parent[v] is not None:
        P[idx, parent_edge[v]] = 1
        v = parent[v]
P_inv = np.linalg.inv(P)


# ============================================================
# PART 1: SIGNED WEIGHT ANALYSIS
# ============================================================
print("\n" + "=" * 72)
print("PART 1: SIGNED WEIGHT SOLUTION ANALYSIS")
print("=" * 72)

# The unique signed solution with |w| = {0,1,2,3,5,6,8,9}
signed_w = np.array([5, -2, 8, 0, 6, 9, -1, 3])
signed_assign = {0: 0, 1: 5, 2: 19, 3: 3, 4: 26, 5: 16, 6: 11, 7: 11, 8: 17}

print(f"\n  Signed weights: {signed_w.tolist()}")
print(f"  |weights|: {sorted(np.abs(signed_w))}")
print(f"  Target:    [0, 1, 2, 3, 5, 6, 8, 9]")

print(f"\n  Edge analysis (sign vs Galois sector):")
for e_idx, (a, b) in enumerate(EDGES):
    w = signed_w[e_idx]
    sec_a = GALOIS_SECTOR[a]
    sec_b = GALOIS_SECTOR[b]
    crossing = sec_a != sec_b
    cross_type = "CROSS" if crossing else "same"
    fa = N_TO_F.get(signed_assign[a], '?')
    fb = N_TO_F.get(signed_assign[b], '?')
    sign = "+" if w >= 0 else "-"
    print(f"  rho_{a+1}({sec_a:>5}) --[{w:+3d}]--> rho_{b+1}({sec_b:>5})  "
          f"[{cross_type:>5}]  {fa}({signed_assign[a]}) -> {fb}({signed_assign[b]})")

# Check if negative weights are on specific structural edges
neg_edges = [(a, b, e_idx) for e_idx, (a, b) in enumerate(EDGES) if signed_w[e_idx] < 0]
print(f"\n  Negative weight edges: {len(neg_edges)}")
for a, b, e_idx in neg_edges:
    print(f"    rho_{a+1}(d={DIMS[a]}) -- rho_{b+1}(d={DIMS[b]}): w={signed_w[e_idx]}")
    print(f"    Sectors: {GALOIS_SECTOR[a]} -> {GALOIS_SECTOR[b]}")
    # Check if this edge is "going backwards" (child has smaller exponent)
    na, nb = signed_assign[a], signed_assign[b]
    if parent[b] == a:
        direction = "parent->child"
        exponent_dir = "DECREASING" if nb < na else "increasing"
    elif parent[a] == b:
        direction = "child->parent"
        exponent_dir = "N/A"
    else:
        direction = "?"
        exponent_dir = "?"
    print(f"    Direction: {direction}, exponent: {exponent_dir} ({na}->{nb})")


# ============================================================
# PART 2: EFFECTIVE RESISTANCE
# ============================================================
print("\n" + "=" * 72)
print("PART 2: EFFECTIVE RESISTANCE WITH CAYLEY-WEIGHTED EDGES")
print("=" * 72)

def compute_eff_resistance(edge_weights):
    """Compute effective resistance from node 0 to all other nodes."""
    # Build weighted Laplacian: L_ij = -w_ij, L_ii = sum(w_ij for j~i)
    L = np.zeros((9, 9))
    for e_idx, (a, b) in enumerate(EDGES):
        w = edge_weights[e_idx]
        if w > 1e-15:
            L[a, b] -= w
            L[b, a] -= w
            L[a, a] += w
            L[b, b] += w

    # Pseudoinverse (L is singular with null space = constant vector)
    L_pinv = np.linalg.pinv(L)

    # R_eff(i,j) = L_pinv[i,i] + L_pinv[j,j] - 2*L_pinv[i,j]
    R = np.zeros(9)
    for i in range(9):
        R[i] = L_pinv[0, 0] + L_pinv[i, i] - 2 * L_pinv[0, i]
    return R

# Weight scheme A: unit weights
R_unit = compute_eff_resistance([1.0]*8)
print(f"\n  A. Unit weights:")
print(f"  {'Node':>8} {'R_eff':>10} {'n_target':>10} {'ratio':>10}")
for i in range(9):
    ratio = TARGET[i] / R_unit[i] if R_unit[i] > 1e-10 else 0
    print(f"  rho_{i+1}  {R_unit[i]:10.4f}  {TARGET[i]:10d}  {ratio:10.4f}")

# Weight scheme B: Cayley eigenvalue weights
# Edge weight = geometric mean of endpoint Cayley eigenvalues
cayley_w = []
for a, b in EDGES:
    la, lb = LAP_EIGS[a], LAP_EIGS[b]
    cayley_w.append(np.sqrt(max(la * lb, 0.001)))  # avoid zero
R_cayley = compute_eff_resistance(cayley_w)
print(f"\n  B. Cayley eigenvalue weights (geometric mean):")
print(f"  Edge weights: {[f'{w:.3f}' for w in cayley_w]}")
print(f"  {'Node':>8} {'R_eff':>10} {'n_target':>10} {'ratio':>10}")
for i in range(9):
    ratio = TARGET[i] / R_cayley[i] if R_cayley[i] > 1e-10 else 0
    print(f"  rho_{i+1}  {R_cayley[i]:10.4f}  {TARGET[i]:10d}  {ratio:10.4f}")

# Weight scheme C: dimension weights
dim_w = [np.sqrt(DIMS[a] * DIMS[b]) for a, b in EDGES]
R_dim = compute_eff_resistance(dim_w)
print(f"\n  C. Dimension weights (geometric mean of dims):")
print(f"  {'Node':>8} {'R_eff':>10} {'n_target':>10} {'ratio':>10}")
for i in range(9):
    ratio = TARGET[i] / R_dim[i] if R_dim[i] > 1e-10 else 0
    print(f"  rho_{i+1}  {R_dim[i]:10.4f}  {TARGET[i]:10d}  {ratio:10.4f}")

# Weight scheme D: 1/Cayley eigenvalue (conductance ~ mass)
inv_cayley_w = []
for a, b in EDGES:
    la, lb = LAP_EIGS[a], LAP_EIGS[b]
    mean = (la + lb) / 2
    inv_cayley_w.append(1.0 / max(mean, 0.1))
R_inv = compute_eff_resistance(inv_cayley_w)
print(f"\n  D. Inverse Cayley weights (1/mean lambda):")
print(f"  {'Node':>8} {'R_eff':>10} {'n_target':>10} {'ratio':>10}")
for i in range(9):
    ratio = TARGET[i] / R_inv[i] if R_inv[i] > 1e-10 else 0
    print(f"  rho_{i+1}  {R_inv[i]:10.4f}  {TARGET[i]:10d}  {ratio:10.4f}")


# ============================================================
# PART 3: GREEN FUNCTION
# ============================================================
print("\n" + "=" * 72)
print("PART 3: GREEN FUNCTION G(rho_1, rho_i)")
print("=" * 72)

def compute_green(edge_weights, epsilon=0.01):
    """G = (L + eps*I)^{-1}, return G[0, :] and G[:, 0]."""
    L = np.zeros((9, 9))
    for e_idx, (a, b) in enumerate(EDGES):
        w = edge_weights[e_idx]
        if w > 1e-15:
            L[a, b] -= w
            L[b, a] -= w
            L[a, a] += w
            L[b, b] += w
    G = np.linalg.inv(L + epsilon * np.eye(9))
    return G

# Several epsilon values
for eps_name, eps_val in [("0.01", 0.01), ("1/phi^2", 1/phi**2), ("1/phi^4", 1/phi**4)]:
    G = compute_green([1.0]*8, float(eps_val))
    g_vals = G[0, :]  # Green function from rho_1

    # Try: n_f ~ C * (-log(G[0,i]) / log(phi))
    print(f"\n  Green function with eps = {eps_name}:")
    print(f"  {'Node':>8} {'G(0,i)':>12} {'-ln(G)/ln(phi)':>15} {'n_target':>10}")
    for i in range(9):
        if g_vals[i] > 1e-15:
            log_val = -np.log(g_vals[i]) / np.log(phi)
        else:
            log_val = float('inf')
        print(f"  rho_{i+1}  {g_vals[i]:12.6f}  {log_val:15.4f}  {TARGET[i]:10d}")


# ============================================================
# PART 4: HEAT KERNEL
# ============================================================
print("\n" + "=" * 72)
print("PART 4: HEAT KERNEL K_t(rho_1, rho_i)")
print("=" * 72)

def heat_kernel(edge_weights, t):
    """K_t = e^{-tL}, return K[0, :]."""
    L = np.zeros((9, 9))
    for e_idx, (a, b) in enumerate(EDGES):
        w = edge_weights[e_idx]
        if w > 1e-15:
            L[a, b] -= w
            L[b, a] -= w
            L[a, a] += w
            L[b, b] += w

    evals, evecs = np.linalg.eigh(L)
    K = evecs @ np.diag(np.exp(-t * evals)) @ evecs.T
    return K

# Test various t values
for t_name, t_val in [("1/phi", 1/phi), ("phi", phi), ("1", 1.0),
                       ("phi^2", phi**2), ("1/phi^2", 1/phi**2),
                       ("a1", float(a1)), ("log(phi)", np.log(phi))]:
    K = heat_kernel([1.0]*8, float(t_val))
    k_vals = K[0, :]

    print(f"\n  t = {t_name} ({float(t_val):.4f}):")
    print(f"  {'Node':>8} {'K_t(0,i)':>12} {'-ln(K)/ln(phi)':>15} {'n_target':>10} {'error':>8}")

    log_targets = []
    for i in range(9):
        if k_vals[i] > 1e-15:
            log_val = -np.log(k_vals[i]) / np.log(phi)
        else:
            log_val = float('inf')
        err = log_val - TARGET[i] if log_val < 1e10 else float('inf')
        log_targets.append((log_val, TARGET[i], err))
        print(f"  rho_{i+1}  {k_vals[i]:12.6e}  {log_val:15.4f}  {TARGET[i]:10d}  {err:+8.2f}")

    # RMS error
    errors = [x[2] for x in log_targets if abs(x[2]) < 1e10]
    if errors:
        rms = np.sqrt(np.mean(np.array(errors)**2))
        print(f"  RMS error in log_phi: {rms:.4f}")

# Also try Cayley-weighted Laplacian heat kernel
print(f"\n  --- Heat kernel with Cayley-weighted edges ---")
cayley_w_hk = []
for a, b in EDGES:
    cayley_w_hk.append(np.sqrt(max(LAP_EIGS[a] * LAP_EIGS[b], 0.01)))

for t_name, t_val in [("1/phi", 1/phi), ("phi", phi), ("1", 1.0),
                       ("1/(2*phi^3)", 1/(2*phi**3))]:
    K = heat_kernel(cayley_w_hk, float(t_val))
    k_vals = K[0, :]

    print(f"\n  t = {t_name} ({float(t_val):.4f}), Cayley weights:")
    errors = []
    for i in range(9):
        if k_vals[i] > 1e-300:
            log_val = -np.log(k_vals[i]) / np.log(phi)
        else:
            log_val = float('inf')
        err = log_val - TARGET[i] if log_val < 1e10 else float('inf')
        errors.append(err)
        print(f"  rho_{i+1}  {k_vals[i]:12.6e}  log_phi={log_val:8.2f}  n={TARGET[i]:2d}  err={err:+6.2f}")

    valid = [e for e in errors if abs(e) < 1e10]
    if valid:
        rms = np.sqrt(np.mean(np.array(valid)**2))
        print(f"  RMS error: {rms:.4f}")


# ============================================================
# PART 5: NONLINEAR PATH METRICS
# ============================================================
print("\n" + "=" * 72)
print("PART 5: NONLINEAR PATH METRICS")
print("=" * 72)

# For each node, the path from rho_1 is unique (tree)
# Test nonlinear accumulation along the path

def get_path(node):
    """Return list of (edge_idx, a, b) from node 0 to given node."""
    path = []
    v = node
    while parent[v] is not None:
        e = parent_edge[v]
        a, b = EDGES[e]
        path.append((e, a, b, v))  # edge, endpoints, which is child
        v = parent[v]
    path.reverse()
    return path

# Test A: Product accumulation d = product(phi^{w_e}) = phi^{sum w_e}
# This is equivalent to the standard weighted distance, so same as before.
# But with dimension-dependent weights:

print(f"\n  A. Weighted path: sum of w_e * dim(child) / dim(parent)")
# Using the best non-negative solution (Sol 3): [3,2,6,0,5,3,1,9]
sol3_w = [3, 2, 6, 0, 5, 3, 1, 9]
sol3_assign = {0:0, 1:3, 2:26, 3:5, 4:19, 5:17, 6:11, 7:11, 8:16}

for node in range(9):
    path = get_path(node)
    total = 0
    for e, a, b, child in path:
        p = a if a != child else b  # parent
        total += sol3_w[e] * DIMS[child] / DIMS[p]
    print(f"  rho_{node+1}: dim-weighted dist = {total:.4f}, n = {sol3_assign[node]}")

# Test B: Product of (1 + w_e * dim(node)/dim(branch))
print(f"\n  B. Multiplicative: product(1 + w_e/b1) along path")
for node in range(9):
    path = get_path(node)
    product = 1.0
    for e, a, b, child in path:
        product *= (1 + sol3_w[e] / b1)
    log_val = np.log(product) / np.log(phi) if product > 0 else 0
    print(f"  rho_{node+1}: prod = {product:.4f}, log_phi = {log_val:.4f}, n = {sol3_assign[node]}")

# Test C: Sum w_e * (depth of edge)
print(f"\n  C. Depth-weighted: sum w_e * depth(edge)")
# Build depth using BFS order
depth = {0: 0}
bfs_q = deque([0])
bfs_vis = {0}
while bfs_q:
    v = bfs_q.popleft()
    for e_idx, (a, b) in enumerate(EDGES):
        for u in [a, b]:
            other = b if u == a else a
            if other == v and u not in bfs_vis:
                bfs_vis.add(u)
                depth[u] = depth[v] + 1
                bfs_q.append(u)

for node in range(9):
    path = get_path(node)
    total = 0
    for e, a, b, child in path:
        total += sol3_w[e] * depth[child]
    print(f"  rho_{node+1}: depth-weighted = {total}, n = {sol3_assign[node]}")


# ============================================================
# PART 6: OPTIMAL CAYLEY-BASED METRIC SEARCH
# ============================================================
print("\n" + "=" * 72)
print("PART 6: SYSTEMATIC METRIC SEARCH")
print("=" * 72)

# For each metric candidate, compute d(0, i) for all i and find best
# linear fit n_i = C * d(0, i) + D. Report RMS.

# Node properties available:
# dim[i], LAP_EIGS[i], depth[i]

# Build distance functions parametrized by node function f(i):
# d(0, node) = sum over edges on path of |f(child) - f(parent)| or f(child) etc.

# Try: n(node) = sum over path of f(edge_endpoints)
node_funcs = {
    'dim': lambda i: DIMS[i],
    'dim^2': lambda i: DIMS[i]**2,
    'lambda': lambda i: LAP_EIGS[i],
    'lambda/dim': lambda i: LAP_EIGS[i]/DIMS[i],
    'lambda*dim': lambda i: LAP_EIGS[i]*DIMS[i],
    'dim*(dim+1)/2': lambda i: DIMS[i]*(DIMS[i]+1)//2,
}

print(f"\n  Testing: d(0, node) = sum_edge f(child)")
for name, f in node_funcs.items():
    dists = np.zeros(9)
    for node in range(9):
        path = get_path(node)
        for e, a, b, child in path:
            dists[node] += f(child)

    # Linear fit: n = C * d + D
    valid = [(dists[i], TARGET[i]) for i in range(9)]
    d_arr = np.array([x[0] for x in valid])
    n_arr = np.array([x[1] for x in valid])
    A_mat = np.column_stack([d_arr, np.ones(9)])
    result = np.linalg.lstsq(A_mat, n_arr, rcond=None)
    C, D = result[0]
    pred = C * d_arr + D
    rms = np.sqrt(np.mean((pred - n_arr)**2))
    print(f"  f={name:>15s}: n = {C:.4f}*d + {D:.4f}, RMS = {rms:.4f}")

print(f"\n  Testing: d(0, node) = sum_edge |f(child) - f(parent)|")
for name, f in node_funcs.items():
    dists = np.zeros(9)
    for node in range(9):
        path = get_path(node)
        for e, a, b, child in path:
            p = a if a != child else b
            dists[node] += abs(f(child) - f(p))

    A_mat = np.column_stack([dists, np.ones(9)])
    result = np.linalg.lstsq(A_mat, TARGET.astype(float), rcond=None)
    C, D = result[0]
    pred = C * dists + D
    rms = np.sqrt(np.mean((pred - TARGET)**2))
    print(f"  f={name:>15s}: n = {C:.4f}*d + {D:.4f}, RMS = {rms:.4f}")

# Try: d(0, node) = sum_edge f(child) * g(parent)
print(f"\n  Testing: d(0, node) = sum_edge dim(child) * lambda(parent)")
dists = np.zeros(9)
for node in range(9):
    path = get_path(node)
    for e, a, b, child in path:
        p = a if a != child else b
        dists[node] += DIMS[child] * LAP_EIGS[p]

A_mat = np.column_stack([dists, np.ones(9)])
result = np.linalg.lstsq(A_mat, TARGET.astype(float), rcond=None)
C, D = result[0]
pred = C * dists + D
rms = np.sqrt(np.mean((pred - TARGET)**2))
print(f"  n = {C:.6f}*d + {D:.4f}, RMS = {rms:.4f}")


# ============================================================
# PART 7: COMBINATORIAL OPTIMIZATION - BEST EDGE WEIGHTS
# ============================================================
print("\n" + "=" * 72)
print("PART 7: OPTIMAL EDGE WEIGHTS FOR EXACT FIT")
print("=" * 72)

# For each of the 3 non-negative solutions, display the full analysis
non_degen = [3, 5, 16, 17, 19, 26]
all_sols = []

for degen_pair in combinations(range(8), 2):
    other_idx = [i for i in range(8) if i not in degen_pair]
    for perm in permutations(range(6)):
        t = np.zeros(8)
        t[degen_pair[0]] = 11
        t[degen_pair[1]] = 11
        for k in range(6):
            t[other_idx[k]] = non_degen[perm[k]]
        w = P_inv @ t
        if np.min(w) < -1e-10:
            continue
        w_r = np.round(w)
        if not np.allclose(w, w_r, atol=1e-10):
            continue
        w_int = w_r.astype(int).tolist()
        assignment = {0: 0}
        for idx in range(8):
            assignment[remaining[idx]] = int(t[idx])
        all_sols.append({'assign': assignment, 'w': w_int})

print(f"\n  Found {len(all_sols)} non-negative integer solutions")

# Analyze each
for s_idx, sol in enumerate(all_sols):
    w = sol['w']
    assign = sol['assign']
    w_set = sorted(set(w))
    fw = {0, 1, 2, 3, 5, 6, 8, 9}
    n_fw = sum(1 for x in w if x in fw)
    all_fw = all(x in fw for x in w)

    print(f"\n  Solution {s_idx+1}: weights = {sorted(w)}")
    print(f"    All framework constants: {all_fw} ({n_fw}/8)")
    print(f"    Assignment:")
    for node in range(9):
        n = assign[node]
        f = N_TO_F.get(n, '?')
        print(f"      rho_{node+1}(d={DIMS[node]}): {f}(n={n})")

    # Galois pair analysis
    print(f"    Galois pairs:")
    for a, b in GALOIS_PAIRS:
        na, nb = assign[a], assign[b]
        fa, fb = N_TO_F.get(na, '?'), N_TO_F.get(nb, '?')
        diff = abs(na - nb)
        print(f"      rho_{a+1}<->rho_{b+1}: {fa}({na})<->{fb}({nb}), diff={diff}")

    # Special: is the (rho_6, rho_7) non-Galois pair interesting?
    na6, na7 = assign[5], assign[6]
    f6, f7 = N_TO_F.get(na6, '?'), N_TO_F.get(na7, '?')
    print(f"      rho_6<->rho_7 (self-conj pair): {f6}({na6})<->{f7}({na7}), diff={abs(na6-na7)}")

    # Check if branch weights have framework meaning
    main_edges = [0, 1, 2, 3, 4]  # edges on main chain
    branch_edges = [5, 6, 7]  # edges from branch
    main_sum = sum(w[e] for e in main_edges)
    branch_sum = sum(w[e] for e in branch_edges)
    print(f"    Main chain sum: {main_sum} (= n at branch point)")
    print(f"    Branch sum: {branch_sum}")


# ============================================================
# PART 8: DIMENSION-WEIGHTED DISTANCE (KEY TEST)
# ============================================================
print("\n" + "=" * 72)
print("PART 8: DIMENSION-WEIGHTED DISTANCE")
print("=" * 72)
print("""
  Key idea: instead of d(0,i) = sum w_e,
  try d(0,i) = sum w_e * F(dim) for some function of endpoint dims.

  This preserves the tree-path structure but modulates by dimension.
""")

# For Solution 3 (the best one with all framework weights)
sol3_w = [3, 2, 6, 0, 5, 3, 1, 9]
sol3_assign = {0:0, 1:3, 2:26, 3:5, 4:19, 5:17, 6:11, 7:11, 8:16}

# The "problematic" edges are the branch ones:
# rho_9->rho_5: w=3 but we want distance contribution to be 19-16=3. OK!
# rho_9->rho_6: w=1 but we want distance contribution to be 17-16=1. OK!
# rho_6->rho_3: w=9 but we want distance contribution to be 26-17=9. OK!

# The issue is that w=3 appears twice (main: rho_1->rho_2, branch: rho_9->rho_5)
# What if the EFFECTIVE weight is w * dim_correction?

# Try: n_f = sum over path of w_e * (dim_child / dim_parent)
print(f"\n  Test: n = sum w_e * dim(child)/dim(parent)")
for node in range(9):
    path = get_path(node)
    total = 0
    for e, a, b, child in path:
        p = a if a != child else b
        total += sol3_w[e] * DIMS[child] / DIMS[p]
    actual = sol3_assign[node]
    print(f"  rho_{node+1}: predicted={total:.3f}, actual={actual}, err={total-actual:.3f}")

# Try: n_f = sum over path of w_e * (dim_child + dim_parent) / 2
print(f"\n  Test: n = sum w_e * (dim_c + dim_p) / (2 * dim_p)")
for node in range(9):
    path = get_path(node)
    total = 0
    for e, a, b, child in path:
        p = a if a != child else b
        total += sol3_w[e] * (DIMS[child] + DIMS[p]) / (2 * DIMS[p])
    actual = sol3_assign[node]
    print(f"  rho_{node+1}: predicted={total:.3f}, actual={actual}, err={total-actual:.3f}")


print(f"\n{'='*72}")
print("EXP-323c COMPLETE")
print(f"{'='*72}")

"""
exp387_unified_mass_attempt.py
================================
GOAL: Find a SINGLE unified formula for all 9 fermion mass corrections.

Previous results:
  exp384: f(A)_{0,k} FAILS (Green's fn, heat kernel, adjacency powers)
  exp384c: Fourier decomposition - coefficients not "nice"
  exp364: Reduced 4->2 prescriptions (Q+L), proved 2 is minimum

NEW APPROACHES in this experiment:
  1. Edge-additive Wilson line decomposition
  2. Galois pair sum/difference identities (NEW observation)
  3. Bilinear form B(spectral, algebraic) search
  4. Character PRODUCT formula
  5. Regularized height with Galois correction

BLIND: compute first, interpret after.

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from math import sqrt, log, pi

PHI = (1 + sqrt(5)) / 2
PHIc = -1/PHI
a1 = 5
b1 = 6
N_gen = 3
C = 4.0 / (a1**2 + 1)  # = 2/13
d_ST = 4  # = phi^3 + phi'^3

# =====================================================================
# SETUP: McKay graph, character table, fermion data
# =====================================================================
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
fermions = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
ab_vals = [(0,0), (3,-2), (1,0), (1,1), (1,1), (2,1), (1,2), (4,1), (-1,4)]
T3_vals = [-0.5, +0.5, -0.5, -0.5, -0.5, +0.5, -0.5, +0.5, -0.5]
Nc_vals = [1, 3, 3, 3, 1, 3, 1, 3, 3]
n_bare = [5*a + 6*b for a,b in ab_vals]
delta_known = [0.0, 0.0, -0.4000, -0.17627, 0.07922, 0.24759, -0.05519, 0.45299, -0.28003]
colors = ['W', 'B', 'W', 'B', 'W', 'B', 'W', 'B', 'W']

# Galois pairs on McKay graph
galois_pairs = [(1,7), (2,8), (3,6)]  # (u,t), (d,b), (s,tau)
galois_self = [0, 4, 5]  # e, mu, c

def zphi_norm(a, b):
    return a**2 + a*b - b**2

def zphi_val(a, b):
    return a + b * PHI

def galois_conj(a, b):
    return a + b * PHIc

norms = [zphi_norm(a, b) for a, b in ab_vals]
z_vals = [zphi_val(a, b) for a, b in ab_vals]
zp_vals = [galois_conj(a, b) for a, b in ab_vals]

# Build McKay adjacency and character table
A_adj = np.zeros((9,9))
for i,j in edges:
    A_adj[i,j] = A_adj[j,i] = 1

evals_raw, evecs_raw = np.linalg.eigh(A_adj)
idx = np.argsort(evals_raw)[::-1]
evals = evals_raw[idx]
evecs = evecs_raw[:, idx]
for j in range(9):
    if evecs[0,j] < 0:
        evecs[:,j] *= -1

class_sizes = [1, 12, 20, 12, 30, 12, 20, 12, 1]
N_group = 120
chi = np.zeros((9,9))
for c in range(9):
    norm = sqrt(N_group / class_sizes[c])
    for k in range(9):
        chi[k,c] = evecs[k,c] * norm

# Galois anti-symmetric G
COL_10A, COL_5A = 1, 7
G_vals_spectral = [(chi[k, COL_10A]/dims[k] - chi[k, COL_5A]/dims[k]) * 6 for k in range(9)]

# Graph distance from node 0
def bfs_dist(start):
    dist = [-1]*9
    dist[start] = 0
    q = [start]
    while q:
        n = q.pop(0)
        for i,j in edges:
            neighbor = j if i == n else (i if j == n else -1)
            if neighbor >= 0 and dist[neighbor] < 0:
                dist[neighbor] = dist[n] + 1
                q.append(neighbor)
    return dist

dist_from_0 = bfs_dist(0)

# Path from node 0 to each node
def path_to(target):
    if target == 0: return [0]
    parent = {0: None}
    q = [0]
    while q:
        n = q.pop(0)
        for i,j in edges:
            nb = j if i == n else (i if j == n else -1)
            if nb >= 0 and nb not in parent:
                parent[nb] = n
                q.append(nb)
                if nb == target:
                    path = [target]
                    while path[-1] != 0:
                        path.append(parent[path[-1]])
                    return list(reversed(path))
    return [0]

paths = [path_to(k) for k in range(9)]

print("=" * 78)
print("exp387: UNIFIED MASS CORRECTION - FRESH APPROACHES")
print("=" * 78)

# =====================================================================
# PART 1: EDGE-ADDITIVE WILSON LINE DECOMPOSITION
# =====================================================================
print("\n" + "=" * 78)
print("PART 1: Edge-Additive Wilson Line Decomposition")
print("=" * 78)

# Idea: delta_k = sum over edges on path from 0 to k of f(edge)
# Since the tree has unique paths, we can INVERT: compute f(edge) from delta values.

# Incremental correction at each edge: f(i,j) = delta_j - delta_i
# (where j is farther from root than i, delta_0 = 0)

# Order: compute delta at each node as cumulative sum of edge contributions
# edge_delta(i->j) = delta(j) - delta(i) where i is parent

edge_corrections = {}
for i,j in edges:
    # Determine which is parent (closer to root)
    if dist_from_0[i] < dist_from_0[j]:
        parent, child = i, j
    else:
        parent, child = j, i
    ec = delta_known[child] - delta_known[parent]
    edge_corrections[(parent, child)] = ec

print("\nEdge corrections (incremental delta at each edge):")
print(f"  {'Edge':>8s}  {'dims':>8s}  {'f(edge)':>10s}  {'color':>6s}")
print("-" * 40)
for (p,c), ec in sorted(edge_corrections.items()):
    dp, dc = dims[p], dims[c]
    col_change = f"{colors[p]}->{colors[c]}"
    print(f"  {p}->{c}      {dp}->{dc}      {ec:+10.5f}   {col_change}")

# Check: can we express f(edge) in terms of (d_parent, d_child, color)?
print("\nTest: f(edge) = alpha * (d_child - d_parent) + beta?")
d_diffs = []
f_edges = []
for (p,c), ec in sorted(edge_corrections.items()):
    d_diffs.append(dims[c] - dims[p])
    f_edges.append(ec)

d_diffs = np.array(d_diffs, dtype=float)
f_edges = np.array(f_edges)

# Linear fit: f = alpha * d_diff + beta
X = np.column_stack([d_diffs, np.ones(len(d_diffs))])
params = np.linalg.lstsq(X, f_edges, rcond=None)[0]
pred = X @ params
rms = sqrt(np.mean((pred - f_edges)**2))
print(f"  Linear fit: f = {params[0]:.4f} * (d_c - d_p) + {params[1]:.4f}")
print(f"  RMS = {rms:.5f}")

# Try: f depends on d_child, d_parent, AND bipartite direction
print("\n  Edge corrections with dim product d_p * d_c:")
for (p,c), ec in sorted(edge_corrections.items()):
    dp, dc = dims[p], dims[c]
    prod = dp * dc
    ratio = ec / prod if prod > 0 else 0
    ratio_C = ratio / C if abs(ratio) > 1e-10 else 0
    print(f"    {p}->{c}: f={ec:+.5f}, d_p*d_c={prod:3d}, f/(d*d)={ratio:+.6f}, "
          f"f/(C*d*d)={ratio_C:+.5f}")

# =====================================================================
# PART 2: GALOIS PAIR IDENTITIES (NEW OBSERVATION)
# =====================================================================
print("\n" + "=" * 78)
print("PART 2: Galois Pair Identities")
print("=" * 78)

print("\nGalois pairs: corrections and their sums/differences")
print(f"  {'Pair':>10s}  {'delta_k':>8s}  {'delta_k*':>8s}  {'Sum':>8s}  {'Diff':>8s}")
print("-" * 55)

for k1, k2 in galois_pairs:
    d1 = delta_known[k1]
    d2 = delta_known[k2]
    s = d1 + d2
    diff = d1 - d2
    print(f"  ({fermions[k1]:>3s},{fermions[k2]:>3s})  {d1:+8.5f}  {d2:+8.5f}  {s:+8.5f}  {diff:+8.5f}")

# Test: are the differences related?
d_db = delta_known[2] - delta_known[8]  # d - b
d_st = delta_known[3] - delta_known[6]  # s - tau
d_ut = delta_known[1] - delta_known[7]  # u - t

print(f"\n  delta(d) - delta(b)   = {d_db:+.6f}")
print(f"  delta(s) - delta(tau) = {d_st:+.6f}")
print(f"  delta(u) - delta(t)   = {d_ut:+.6f}")

# Test: is delta(d)-delta(b) = -N_gen/a1^2?
target = -N_gen / a1**2
print(f"\n  Test: -N_gen/a1^2 = -3/25 = {target:.6f}")
print(f"    d-b: {d_db:.6f}  (diff from -3/25: {abs(d_db - target):.6f}, {abs(d_db-target)/abs(target)*100:.3f}%)")
print(f"    s-tau: {d_st:.6f}  (diff from -3/25: {abs(d_st - target):.6f}, {abs(d_st-target)/abs(target)*100:.3f}%)")

# Test: are the SUMS related?
s_db = delta_known[2] + delta_known[8]  # d + b
s_st = delta_known[3] + delta_known[6]  # s + tau
s_ut = delta_known[1] + delta_known[7]  # u + t

print(f"\n  Galois pair SUMS:")
print(f"    delta(d) + delta(b)   = {s_db:+.6f}")
print(f"    delta(s) + delta(tau) = {s_st:+.6f}")
print(f"    delta(u) + delta(t)   = {s_ut:+.6f}")

# What are these?
print(f"\n    d+b = {s_db:.6f}. -C*ln(19)*(1+1/phi) = {-C*log(19)*(1+1/PHI):.6f}")
print(f"    s+tau = {s_st:.6f}")
print(f"    u+t = {s_ut:.6f} = C*ln(19) = {C*log(19):.6f}")

# Total sum over ALL fermions
total = sum(delta_known)
print(f"\n  Total sum of all corrections: {total:.6f}")
print(f"    -C = {-C:.6f}")
print(f"    -2*C/phi = {-2*C/PHI:.6f}")

# =====================================================================
# PART 3: ALL CORRECTIONS ARE C * D_k
# =====================================================================
print("\n" + "=" * 78)
print("PART 3: Normalized Corrections D_k = delta_k / C")
print("=" * 78)

print(f"\nAll corrections expressed as delta = C * D_k:")
print(f"  {'Name':>4s}  {'delta':>8s}  {'D_k':>10s}  {'Identification':>40s}")
print("-" * 72)

for k in range(9):
    name = fermions[k]
    delta = delta_known[k]
    Dk = delta / C if abs(delta) > 1e-10 else 0.0
    a, b = ab_vals[k]
    Nz = norms[k]

    # Identify D_k
    ident = ""
    if name == 'e':
        ident = "0 (trivial)"
    elif name == 'u':
        ident = "0 (G < 0 => delta=0)"
    elif name == 'd':
        # delta = -2/a1, D = -2/(a1*C) = -2*13/10 = -13/5
        ident = f"-{a1**2+1}/(2*a1) = -{a1**2+1}/{2*a1} = {-(a1**2+1)/(2*a1):.4f}"
    elif name == 's':
        ident = f"-N_gen/phi^2 = {-N_gen/PHI**2:.4f}"
    elif name == 'mu':
        ident = f"+phi^(3/2)/d_ST = {PHI**1.5/d_ST:.4f}"
    elif name == 'c':
        ident = f"+ln(a1) = ln({a1}) = {log(a1):.4f}"
    elif name == 'tau':
        ident = f"-phi^(3/4)/d_ST = {-PHI**0.75/d_ST:.4f}"
    elif name == 't':
        ident = f"+ln(F) = ln({a1*b1-a1-b1}) = {log(a1*b1-a1-b1):.4f}"
    elif name == 'b':
        ident = f"-ln(F)/phi = {-log(a1*b1-a1-b1)/PHI:.4f}"

    print(f"  {name:>4s}  {delta:+8.5f}  {Dk:+10.5f}  {ident:>40s}")

# =====================================================================
# PART 4: BILINEAR FORM B(spectral, algebraic)
# =====================================================================
print("\n" + "=" * 78)
print("PART 4: Bilinear Form B(spectral, algebraic)")
print("=" * 78)

print("\nIdea: delta_k = C * B(v_k, w_k) where:")
print("  v_k = spectral data at McKay node k")
print("  w_k = algebraic data from Z[phi] element z = a + b*phi")

# Spectral features for each node
spec_features = {}
for k in range(9):
    G_k = G_vals_spectral[k]
    dist_k = dist_from_0[k]
    d_k = dims[k]
    # Character ratios at rational classes
    x3a = chi[k, 2] / dims[k]
    x4a = chi[k, 4] / dims[k]
    x6a = chi[k, 6] / dims[k]
    x2a = chi[k, 8] / dims[k]
    spec_features[k] = [G_k, d_k, dist_k, x3a, x4a, x6a, x2a]

# Algebraic features for each fermion
alg_features = {}
for k in range(9):
    a, b = ab_vals[k]
    z = z_vals[k]
    zp = zp_vals[k]
    Nz = norms[k]
    absN = abs(Nz)
    lnN = log(absN) if absN > 0 else 0
    sign_z = 1 if z > 0 else (-1 if z < 0 else 0)
    alg_features[k] = [z, zp, Nz, lnN, sign_z, a, b]

spec_names = ['G', 'd_k', 'dist', 'x(3a)', 'x(4a)', 'x(6a)', 'x(2a)']
alg_names = ['z', 'zp', 'N(z)', 'ln|N|', 'sgn(z)', 'a', 'b']

# Build cross-product feature matrix (exclude e which is input)
n_spec = len(spec_names)
n_alg = len(alg_names)
n_cross = n_spec * n_alg

feature_matrix = []
targets = []
for k in range(1, 9):  # skip e
    sv = spec_features[k]
    av = alg_features[k]
    cross = []
    for s in sv:
        for a in av:
            cross.append(s * a)
    feature_matrix.append(cross)
    targets.append(delta_known[k] / C)  # D_k = delta_k / C

X_cross = np.array(feature_matrix)
y_cross = np.array(targets)

# Use regularized least squares (too many features for 8 data points)
# Instead, try SELECTED cross-terms
selected_pairs = [
    ('G', 'ln|N|'),       # BLACK quarks with |N|>1
    ('G', 'sgn(z)'),      # sign selection for u
    ('d_k', 'ln|N|'),     # dimension weighting
    ('dist', 'b'),        # graph distance * b quantum number
    ('x(6a)', 'zp'),      # character ratio * Galois conjugate
    ('x(4a)', 'N(z)'),    # character ratio * norm
]

print(f"\n  Testing selected bilinear terms:")
for sn, an in selected_pairs:
    si = spec_names.index(sn)
    ai = alg_names.index(an)
    vals = []
    for k in range(1, 9):
        sv = spec_features[k][si]
        av = alg_features[k][ai]
        vals.append(sv * av)
    vals = np.array(vals)
    if np.std(vals) > 1e-10:
        # Linear fit: D_k = alpha * (s*a) + beta
        X_pair = np.column_stack([vals, np.ones(8)])
        params_pair = np.linalg.lstsq(X_pair, y_cross, rcond=None)[0]
        pred_pair = X_pair @ params_pair
        rms_pair = sqrt(np.mean((pred_pair - y_cross)**2))
        corr = np.corrcoef(vals, y_cross)[0,1]
        print(f"    {sn:>6s} x {an:>6s}: corr = {corr:+.4f}, RMS = {rms_pair:.4f}")

# Try 2-feature and 3-feature combinations
print("\n  Systematic 2-feature search (bilinear terms):")
best_rms_2 = 10.0
best_pair_2 = None

for i1 in range(n_cross):
    for i2 in range(i1+1, n_cross):
        X_2 = X_cross[:, [i1, i2]]
        if np.linalg.matrix_rank(X_2) < 2:
            continue
        params_2 = np.linalg.lstsq(X_2, y_cross, rcond=None)[0]
        pred_2 = X_2 @ params_2
        rms_2 = sqrt(np.mean((pred_2 - y_cross)**2))
        if rms_2 < best_rms_2:
            best_rms_2 = rms_2
            best_pair_2 = (i1, i2)
            best_params_2 = params_2

if best_pair_2:
    i1, i2 = best_pair_2
    s1, a1_idx = divmod(i1, n_alg)
    s2, a2_idx = divmod(i2, n_alg)
    print(f"  Best 2-feature: ({spec_names[s1]}*{alg_names[a1_idx]}, "
          f"{spec_names[s2]}*{alg_names[a2_idx]})")
    print(f"    RMS = {best_rms_2:.5f}")
    print(f"    Coefficients: {best_params_2[0]:.6f}, {best_params_2[1]:.6f}")
    # Show predictions
    X_best = X_cross[:, list(best_pair_2)]
    pred_best = X_best @ best_params_2
    print(f"    {'Name':>4s}  {'D_true':>8s}  {'D_pred':>8s}  {'err':>8s}")
    for k in range(1, 9):
        print(f"    {fermions[k]:>4s}  {y_cross[k-1]:+8.4f}  {pred_best[k-1]:+8.4f}  "
              f"{abs(pred_best[k-1]-y_cross[k-1]):8.4f}")

print("\n  Systematic 3-feature search (top 5):")
results_3 = []
for i1 in range(n_cross):
    for i2 in range(i1+1, n_cross):
        for i3 in range(i2+1, n_cross):
            X_3 = X_cross[:, [i1, i2, i3]]
            if np.linalg.matrix_rank(X_3) < 3:
                continue
            params_3 = np.linalg.lstsq(X_3, y_cross, rcond=None)[0]
            pred_3 = X_3 @ params_3
            rms_3 = sqrt(np.mean((pred_3 - y_cross)**2))
            results_3.append((rms_3, i1, i2, i3, params_3))

results_3.sort()
for rank, (rms_val, i1, i2, i3, params) in enumerate(results_3[:5]):
    s1, a1_idx = divmod(i1, n_alg)
    s2, a2_idx = divmod(i2, n_alg)
    s3, a3_idx = divmod(i3, n_alg)
    feat_str = (f"{spec_names[s1]}*{alg_names[a1_idx]}, "
                f"{spec_names[s2]}*{alg_names[a2_idx]}, "
                f"{spec_names[s3]}*{alg_names[a3_idx]}")
    print(f"  #{rank+1}: RMS={rms_val:.5f}  [{feat_str}]")

# =====================================================================
# PART 5: CHARACTER PRODUCT FORMULA
# =====================================================================
print("\n" + "=" * 78)
print("PART 5: Character Product Formula")
print("=" * 78)

# What if delta_k involves chi_k(c) * chi_k(c') for specific class pairs?
# This is a BILINEAR form in the character table.

print("\n  Testing chi_k(c1)/d_k * chi_k(c2)/d_k for all class pairs:")
best_rms_chi = 10.0
best_chi_pair = None

class_names_list = ['1a', '10a', '3a', '10b', '4a', '5b', '6a', '5a', '2a']

for c1 in range(9):
    for c2 in range(c1, 9):
        vals = []
        for k in range(1, 9):
            x1 = chi[k, c1] / dims[k]
            x2 = chi[k, c2] / dims[k]
            vals.append(x1 * x2)
        vals = np.array(vals)
        if np.std(vals) < 1e-10:
            continue
        # Linear fit: D_k = alpha * x1*x2 + beta
        X_chi = np.column_stack([vals, np.ones(8)])
        params_chi = np.linalg.lstsq(X_chi, y_cross, rcond=None)[0]
        pred_chi = X_chi @ params_chi
        rms_chi = sqrt(np.mean((pred_chi - y_cross)**2))
        if rms_chi < best_rms_chi:
            best_rms_chi = rms_chi
            best_chi_pair = (c1, c2)
            best_chi_params = params_chi

if best_chi_pair:
    c1, c2 = best_chi_pair
    print(f"  Best: chi({class_names_list[c1]})/d * chi({class_names_list[c2]})/d")
    print(f"    RMS = {best_rms_chi:.5f}")
    print(f"    D_k = {best_chi_params[0]:.4f} * x({class_names_list[c1]})*x({class_names_list[c2]}) + {best_chi_params[1]:.4f}")

# =====================================================================
# PART 6: THE HONEST ATTEMPT - SINGLE FORMULA WITH Z[phi] REGULARIZATION
# =====================================================================
print("\n" + "=" * 78)
print("PART 6: Regularized Height Formula")
print("=" * 78)

# ALL corrections are C * D_k. The norms are:
#   |N| > 1: D involves ln|N|
#   |N| = 1: D involves algebraic numbers
#   |N| = 0: D = 0
#
# Can we define a "regularized log-norm" that handles |N|=1?
# h_reg(z) = ln|N(z)| + correction_for_units

# Idea: USE the Galois anti-symmetric G from the McKay graph as the
# regularization. For BLACK nodes, G != 0 and gives a correction.
# For WHITE nodes, G = 0 and we need something else.

# What if: h_eff = ln|N| + beta * G / phi^2  (G provides the unit-sector correction)?
# For primes: ln|N| dominates, G/phi^2 is small correction
# For units: ln|N| = 0, G/phi^2 provides the value
# For WHITE units: BOTH are zero => need ANOTHER term

print("\nTrying: D_k = ln|N| + beta * G_k/phi^2 + gamma * dist_k/a1")
print("  (with sign(T3) and 1/phi factors for quark/lepton split)")

# Build for quarks only first
quark_idx = [k for k in range(1,9) if Nc_vals[k] == 3]
lepton_idx = [k for k in range(1,9) if Nc_vals[k] == 1]

print("\n  QUARKS: delta = sign(T3) * C * h_eff * phi^{(sign(T3)-1)/2}")
print("  back-calculating h_eff from known delta:")
for k in quark_idx:
    T3 = T3_vals[k]
    delta = delta_known[k]
    if T3 > 0:
        h_eff = delta / C if abs(delta) > 1e-10 else 0
    else:
        h_eff = -delta * PHI / C if abs(delta) > 1e-10 else 0
    Nz = abs(norms[k])
    lnN = log(Nz) if Nz > 1 else 0
    G = G_vals_spectral[k]
    dist_k = dist_from_0[k]
    print(f"    {fermions[k]:>3s}: h_eff = {h_eff:+8.4f}, ln|N| = {lnN:6.4f}, "
          f"G = {G:+8.4f}, dist = {dist_k}, |N| = {Nz}")

# For quarks: try h_eff = alpha*ln|N| + beta*G/phi^2 + gamma*(dist/a1)^2
print("\n  Fitting: h_eff = a*ln|N| + b*G/phi^2 + c*(dist/a1)")
h_effs_q = []
features_q = []
for k in quark_idx:
    T3 = T3_vals[k]
    delta = delta_known[k]
    if T3 > 0:
        h = delta / C if abs(delta) > 1e-10 else 0
    else:
        h = -delta * PHI / C if abs(delta) > 1e-10 else 0
    h_effs_q.append(h)

    Nz = abs(norms[k])
    lnN = log(Nz) if Nz > 1 else 0
    G = G_vals_spectral[k]
    dist_k = dist_from_0[k]
    features_q.append([lnN, G/PHI**2, dist_k/a1])

X_q = np.array(features_q)
y_q = np.array(h_effs_q)
params_q = np.linalg.lstsq(X_q, y_q, rcond=None)[0]
pred_q = X_q @ params_q
rms_q = sqrt(np.mean((pred_q - y_q)**2))

print(f"  Coefficients: ln|N|: {params_q[0]:.4f}, G/phi^2: {params_q[1]:.4f}, dist/a1: {params_q[2]:.4f}")
print(f"  RMS = {rms_q:.5f}")

for i, k in enumerate(quark_idx):
    print(f"    {fermions[k]:>3s}: h_true = {y_q[i]:+8.4f}, h_pred = {pred_q[i]:+8.4f}, "
          f"err = {abs(pred_q[i]-y_q[i]):.5f}")

# =====================================================================
# PART 7: THE FUNDAMENTAL QUESTION - WHY 2 IS THE MINIMUM
# =====================================================================
print("\n" + "=" * 78)
print("PART 7: Why Quarks and Leptons CANNOT Be Unified")
print("=" * 78)

# Compute the mass correction MECHANISM for each fermion
print("\n  Correction mechanisms by type:")
print("  QUARKS:    delta ~ C * ln|N|  or  C * (algebraic/a1)")
print("  LEPTONS:   delta ~ C * phi^alpha * |z'|^(3/4)")
print()

# Key difference: quarks use LOGARITHMIC corrections, leptons use POWER-LAW
# Can we express leptons as logarithms?
for k in lepton_idx:
    delta = delta_known[k]
    Dk = delta / C
    # Can we write Dk = ln(x) for some algebraic x?
    x = np.exp(Dk) if Dk != 0 else 1
    print(f"  {fermions[k]:>3s}: D_k = {Dk:+.5f}, exp(D_k) = {x:.5f}")

print(f"\n  mu: exp(D_mu) = exp({delta_known[4]/C:.5f}) = {np.exp(delta_known[4]/C):.5f}")
print(f"      phi^(1/2) = {PHI**0.5:.5f}")
print(f"  tau: exp(D_tau) = exp({delta_known[6]/C:.5f}) = {np.exp(delta_known[6]/C):.5f}")

# They're NOT simple algebraic numbers => can't be expressed as ln(nice thing)
# But CAN we express quarks as power laws?
print(f"\n  Conversely, quark corrections as power laws:")
for k in quark_idx:
    if abs(delta_known[k]) > 1e-10:
        Dk = delta_known[k] / C
        # D = |z'|^alpha  => alpha = ln|D| / ln|z'|
        zp = abs(zp_vals[k])
        if zp > 1e-10 and zp != 1:
            alpha_exp = log(abs(Dk)) / log(zp) if abs(Dk) > 1e-10 else 0
            print(f"    {fermions[k]:>3s}: D = {Dk:+.4f}, |z'| = {zp:.4f}, "
                  f"alpha = ln|D|/ln|z'| = {alpha_exp:.4f}")

# =====================================================================
# PART 8: NEW IDEA - NORM OF THE CORRECTION VECTOR
# =====================================================================
print("\n" + "=" * 78)
print("PART 8: Galois Norm of the Correction")
print("=" * 78)

# For each correction delta, compute what its "Galois partner" delta' would be
# under phi -> -1/phi
print("\n  If delta involves phi, its Galois conjugate is obtained by phi -> phi':")
for k in range(9):
    name = fermions[k]
    delta = delta_known[k]
    # Under Galois: phi -> -1/phi
    # C is invariant (C = 4/(a1^2+1), no phi)
    # ln|N| is invariant (N is rational)
    # 1/phi -> -phi  (Galois)
    # phi^(3/4) -> |-1/phi|^(3/4) = phi^(-3/4) (for magnitude)

    # For quark prime sector: delta involves C*ln|N| and 1/phi
    # Galois of delta_b = -C*ln(19)/phi is: -C*ln(19)*(-phi) = C*ln(19)*phi
    # For lepton: c_ell*sign(z')*|z'|^(3/4) -> different

    if name in ('c', 't'):
        delta_gal = delta  # C*ln|N| is Galois invariant
    elif name == 'b':
        delta_gal = C * log(19) * PHI  # phi -> -(-phi) = phi in |...|
    elif name == 's':
        delta_gal = -N_gen * C * PHI**2  # 1/phi^2 -> phi^2
    elif name in ('mu', 'tau'):
        a, b = ab_vals[k]
        zp_gal = zphi_val(a, b)  # z' under Galois becomes z
        delta_gal = (C * PHI**3 / d_ST) * np.sign(zp_gal) * abs(zp_gal)**0.75 if name != 'e' else 0
    else:
        delta_gal = delta  # e and u are 0

    if abs(delta) > 1e-10 and abs(delta_gal) > 1e-10:
        norm_product = delta * delta_gal
        print(f"  {name:>3s}: delta = {delta:+.5f}, delta' = {delta_gal:+.5f}, "
              f"delta*delta' = {norm_product:+.6f}")

# =====================================================================
# PART 9: MASTER TABLE AND HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 78)
print("PART 9: MASTER TABLE - All Approaches Compared")
print("=" * 78)

# Current best formula (2 prescriptions)
def delta_current(k):
    """Current 2-prescription formula."""
    name = fermions[k]
    a, b = ab_vals[k]
    T3 = T3_vals[k]
    Nc = Nc_vals[k]
    z = z_vals[k]
    zp = zp_vals[k]
    Nz = norms[k]

    if name == 'e': return 0.0

    # Leptons
    if Nc == 1:
        c_ell = C * PHI**3 / d_ST
        return c_ell * np.sign(zp) * abs(zp)**0.75

    # Quarks
    absN = abs(Nz)
    if absN > 1:
        if T3 > 0:
            return C * log(absN)
        else:
            return -C * log(absN) / PHI
    if b == 0:
        return -2.0 / a1
    if T3 > 0:
        return 0.0
    else:
        return -N_gen * C / PHI**2

print(f"\n  {'Name':>4s}  {'n_bare':>6s}  {'delta':>8s}  {'m_pred':>12s}  {'m_PDG':>12s}  {'Error':>8s}")
print("-" * 60)

m_e = 0.51099895
pdg = {'e':0.511, 'mu':105.66, 'tau':1776.86,
       'u':2.16, 'c':1270., 't':172760., 'd':4.67, 's':93.4, 'b':4180.}
all_errs = []

for k in range(9):
    name = fermions[k]
    nb = n_bare[k]
    delta = delta_current(k)
    n_eff = nb + delta
    m_pred = m_e * PHI**n_eff if name != 'e' else m_e
    m_exp = pdg[name]
    err = (m_pred - m_exp) / m_exp * 100
    if name != 'e':
        all_errs.append(abs(err))
    print(f"  {name:>4s}  {nb:>6d}  {delta:+8.5f}  {m_pred:>12.4f}  {m_exp:>12.4f}  {err:+7.3f}%")

rms_total = sqrt(np.mean(np.array(all_errs)**2))
print(f"\n  RMS (8 fermions): {rms_total:.4f}%")

# =====================================================================
# PART 10: HONEST CONCLUSION
# =====================================================================
print("\n" + "=" * 78)
print("PART 10: HONEST CONCLUSION")
print("=" * 78)

print("""
QUESTION: Can the 2 prescriptions be reduced to 1?

ANSWER: NO, for the following reasons:

1. QUARKS use LOGARITHMIC corrections (ln|N(z)|).
   LEPTONS use POWER-LAW corrections (|z'|^{3/4}).
   These are FUNDAMENTALLY different function classes.
   One cannot be re-expressed as the other without
   introducing new parameters.

2. Within QUARKS, the 3 Z[phi] sectors (prime, rational, unit)
   correspond to different ARITHMETIC structures:
   - Prime: ideal norm (class number theory)
   - Rational: global structure (graph distance / a1)
   - Unit: generation-coupling (N_gen appears)
   These are mathematically distinct sectors of Z[phi].

3. The bipartite split (WHITE/BLACK) on the McKay graph is EXACT:
   - BLACK nodes have G != 0 (Galois anti-symmetric spectral data)
   - WHITE nodes have G = 0 identically
   No single spectral function f(A) can give all corrections.

4. The bilinear form search (spectral x algebraic) finds fits with
   3+ free parameters, but these are OVERFITTING 8 data points,
   not a derivation.

NEW OBSERVATIONS from this experiment:

A) GALOIS PAIR NEAR-IDENTITY:
   delta(d) - delta(b) = -0.1200
   delta(s) - delta(tau) = -0.1211
   Both are CLOSE to -N_gen/a1^2 = -3/25 = -0.1200
   Accuracy: 0% and 0.9% respectively.
   This is an APPROXIMATE pattern (0.014% deviation from
   ln(19) = 91*phi/50), NOT exact.
   STATUS: PATTERN (interesting, not derived)

B) EDGE-ADDITIVE DECOMPOSITION:
   The corrections CAN be decomposed as sums of per-edge
   contributions on the McKay tree. But the edge weights
   don't simplify to a recognizable function of (d_i, d_j).
   STATUS: STRUCTURAL (correct but not illuminating)

C) ALL corrections have the form delta = C * D_k where
   C = 2/13 = 2*sin^2(tW)/N_gen is UNIVERSAL.
   The norms |N| that appear are EXACTLY a1 = 5 and
   F(a1,b1) = 19 (the Frobenius number of a1 and b1).
   STATUS: DERIVED (C is framework constant, norms are a1-determined)

WHAT WOULD CONSTITUTE TRUE UNIFICATION:
   Compute the lattice self-energy Sigma_k on the McKay graph:
   Sigma_k = sum_j G(k,j) * |V_kj|^2
   where G is the full propagator and V is the gauge coupling.
   This would give a SINGLE formula with NO sectors, because
   the self-energy naturally depends on the graph position,
   the coupling structure, and the Z[phi] arithmetic.
   This is FIX 17 (OPEN) from the paper review.

VERDICT: The 2-prescription minimum (Q + L) is IRREDUCIBLE
with current methods. True unification requires a dynamical
calculation (lattice self-energy), not a kinematic formula.
""")

print("=" * 78)
print("EXPERIMENT COMPLETE")
print("=" * 78)

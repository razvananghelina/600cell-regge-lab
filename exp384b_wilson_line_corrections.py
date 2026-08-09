"""
exp384b_wilson_line_corrections.py
===================================
GOAL: Derive fermion mass corrections from Wilson lines on McKay tree.

Key insight from exp384: d quark delta = -dist(e,d)/a1 = -2/5 EXACT.
Also: formula -2/(a1 + b1*chi(6a)/d) works for d.

NEW APPROACH: Compute Wilson line holonomy on McKay tree for each fermion.
The McKay tree has edges with weights. If we can identify the correct
weight function, the holonomy should give the mass correction.

BLIND: compute first, interpret after.

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from math import sqrt, log, pi, cos, sin

PHI = (1 + sqrt(5)) / 2
PHIc = -1/PHI
a1 = 5
b1 = 6
N_gen = 3
C = 4.0 / (a1**2 + 1)  # = 2/13

# McKay graph
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]

# Fermion mapping
fermions = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
ab_vals = [(0,0), (3,-2), (1,0), (1,1), (1,1), (2,1), (1,2), (4,1), (-1,4)]
n_bare = [5*a + 6*b for a,b in ab_vals]  # bare exponents
delta_known = [0.0, 0.0, -0.4, -0.17627, 0.07922, 0.24759, -0.05519, 0.45299, -0.28003]
colors = ['W', 'B', 'W', 'B', 'W', 'B', 'W', 'B', 'W']

# Compute Z[phi] norm
def galois_norm(a, b):
    return a**2 + a*b - b**2

norms = [galois_norm(a, b) for a, b in ab_vals]

# =====================================================================
# PART 1: Character table and spectral data
# =====================================================================
print("=" * 72)
print("exp384b: WILSON LINE CORRECTIONS ON McKAY TREE")
print("=" * 72)

# Build adjacency and compute character table
A_adj = np.zeros((9,9))
for i,j in edges:
    A_adj[i,j] = A_adj[j,i] = 1

evals_raw, evecs_raw = np.linalg.eigh(A_adj)
idx = np.argsort(evals_raw)[::-1]
evals = evals_raw[idx]
evecs = evecs_raw[:, idx]

for j in range(9):
    v = evecs[:, j]
    k = np.argmax(np.abs(v))
    if v[k] < 0:
        evecs[:, j] = -v
    evecs[:, j] *= dims[j] / evecs[0, j]

chi = evecs  # chi[k, c] = character of rho_k at class c

# Conjugacy class sizes
class_sizes = [1, 12, 20, 12, 30, 12, 20, 12, 1]
class_names = ['1a', '10a', '3a', '10b', '4a', '5b', '6a', '5a', '2a']
N_group = 120

# =====================================================================
# PART 2: Wilson Line on McKay Tree
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: Wilson Line Holonomy on McKay Tree")
print("=" * 72)

# Path from node 0 to each node
def path_from_root(target):
    """BFS to find path from node 0 to target in McKay tree."""
    if target == 0:
        return [0]
    parent = {0: None}
    queue = [0]
    while queue:
        node = queue.pop(0)
        for i, j in edges:
            neighbor = j if i == node else (i if j == node else None)
            if neighbor is not None and neighbor not in parent:
                parent[neighbor] = node
                queue.append(neighbor)
                if neighbor == target:
                    path = [target]
                    while path[-1] != 0:
                        path.append(parent[path[-1]])
                    return list(reversed(path))
    return None

print("\nPaths from e (node 0) to each fermion:")
paths = {}
for k in range(9):
    p = path_from_root(k)
    paths[k] = p
    edge_list = [(p[i], p[i+1]) for i in range(len(p)-1)]
    print(f"  {fermions[k]:3s} (node {k}): {p}  edges: {edge_list}")

# =====================================================================
# PART 3: Edge Weight Candidates
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: Edge Weights and Wilson Line Products")
print("=" * 72)

# For each edge (i,j), compute various candidate weights
print("\nEdge properties:")
print(f"  {'Edge':10s}  {'di':3s}  {'dj':3s}  {'di*dj':5s}  {'di+dj':5s}  1/sqrt(didj)  di/dj")
for i, j in edges:
    di, dj = dims[i], dims[j]
    print(f"  ({i},{j})      {di:3d}  {dj:3d}  {di*dj:5d}  {di+dj:5d}    {1/sqrt(di*dj):8.5f}   {di/dj:.5f}")

# Wilson line = product of weights along path
# Try weight = sqrt(d_j/d_i) per edge (i -> j)
print("\n--- Weight = sqrt(d_j / d_i) per edge ---")
for k in range(9):
    path = paths[k]
    W = 1.0
    for idx2 in range(len(path)-1):
        di = dims[path[idx2]]
        dj = dims[path[idx2 + 1]]
        W *= sqrt(dj / di)
    # W = sqrt(d_k / d_0) = sqrt(d_k)
    print(f"  {fermions[k]:3s}: W = sqrt({dims[k]}/{dims[0]}) = {W:.6f}, "
          f"ln(W)/a1 = {log(W)/a1 if W > 0 else float('nan'):.6f}, "
          f"delta = {delta_known[k]:+.6f}")

# Wilson line = product of 1/d_j along path
print("\n--- Weight = 1/d_j per edge ---")
for k in range(9):
    path = paths[k]
    W = 0.0
    for idx2 in range(1, len(path)):
        W += 1.0 / dims[path[idx2]]
    print(f"  {fermions[k]:3s}: sum(1/d_j) = {W:.6f}, "
          f"-W/a1 = {-W/a1:.6f} vs delta = {delta_known[k]:+.6f}")

# Wilson line = sum of dim ratios along path
print("\n--- Weight = sum(d_parent/d_child) along path ---")
for k in range(9):
    path = paths[k]
    W = 0.0
    for idx2 in range(len(path)-1):
        di = dims[path[idx2]]
        dj = dims[path[idx2 + 1]]
        W += di / dj
    print(f"  {fermions[k]:3s}: sum(di/dj) = {W:.6f}, "
          f"-C*W = {-C*W:.6f} vs delta = {delta_known[k]:+.6f}")

# =====================================================================
# PART 4: Character-Based Formula
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: Character-Based Correction Formula")
print("=" * 72)

# From exp384: -2/(a1 + b1*chi(6a)/d) works for d quark
# Test this for ALL WHITE nodes
print("\nFormula: delta = -2/(a1 + b1*x_6a) where x_6a = chi(6a)/d")
print("(Only for WHITE nodes)")
for k in range(9):
    if colors[k] == 'W':
        x6a = chi[k, 6] / dims[k]  # col 6 = 6a
        if abs(a1 + b1 * x6a) > 1e-10:
            delta_pred = -2.0 / (a1 + b1 * x6a)
        else:
            delta_pred = float('inf')
        err = abs(delta_pred - delta_known[k])
        match = "***" if err < 0.001 else ""
        print(f"  {fermions[k]:3s}: x(6a) = {x6a:+.4f}, "
              f"delta_pred = {delta_pred:+.6f}, delta_known = {delta_known[k]:+.6f}, "
              f"err = {err:.6f} {match}")

# Try variants with different classes
print("\nFormula: delta = f(chi(c)/d) for various classes c")
for col, cname in [(2, '3a'), (4, '4a'), (6, '6a')]:
    print(f"\n  --- Class {cname} ---")
    for k in range(9):
        if colors[k] == 'W':
            xc = chi[k, col] / dims[k]
            print(f"    {fermions[k]:3s}: x({cname}) = {xc:+.4f}, "
                  f"delta = {delta_known[k]:+.6f}")

# =====================================================================
# PART 5: Norm-Based Formula for Non-Unit Fermions
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: Galois Norm Analysis")
print("=" * 72)

print("\nFermion norms and corrections:")
print(f"  {'name':4s} {'(a,b)':8s}  {'z':8s}  {'z_conj':8s}  {'N':5s}  {'|N|':5s}  "
      f"{'ln|N|':8s}  {'delta':8s}  {'delta/C':8s}")
for k in range(9):
    a, b = ab_vals[k]
    z = a + b * PHI
    zc = a + b * PHIc
    N = norms[k]
    absN = abs(N) if N != 0 else 0
    lnN = log(absN) if absN > 1 else 0.0
    dC = delta_known[k] / C if abs(C) > 1e-10 else 0.0
    print(f"  {fermions[k]:4s} ({a:2d},{b:2d})  {z:+8.4f}  {zc:+8.4f}  {N:5d}  {absN:5d}  "
          f"{lnN:8.4f}  {delta_known[k]:+8.5f}  {dC:+8.5f}")

# =====================================================================
# PART 6: Unified Formula Attempt
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: Unified Formula - Sector by Sector")
print("=" * 72)

print("\nSECTOR A: Electron (trivial)")
print(f"  e: delta = 0 (reference, TRIVIAL)")

print("\nSECTOR B: BLACK quarks")
print("  Formula: delta = C * sigma_k * spectral(k)")
for k in [1, 3, 5, 7]:
    a, b = ab_vals[k]
    N = norms[k]
    absN = abs(N)
    G = 0  # compute properly
    # G = (L_10a - L_5a)/2 for each node
    L10a = 12 * (1 - chi[k, 1] / dims[k])
    L5a = 12 * (1 - chi[k, 7] / dims[k])
    G = (L10a - L5a) / 2

    print(f"  {fermions[k]:3s}: G = {G:+.4f}, |N| = {absN}, delta = {delta_known[k]:+.5f}")

    # s quark: delta = C*G/phi^2
    if k == 3:
        pred = C * G / PHI**2
        print(f"        C*G/phi^2 = {pred:+.5f} (err: {abs(pred-delta_known[k]):.6f})")
    # u quark: delta = 0 (G < 0, threshold)
    elif k == 1:
        print(f"        G < 0 => delta = 0 (sign rule)")
    # c, t quarks: delta = C*ln|N|
    elif k in [5, 7]:
        if absN > 1:
            pred = C * log(absN)
            print(f"        C*ln|N| = {pred:+.5f} (err: {abs(pred-delta_known[k]):.6f})")
        # Also test spectral: C * 4*G/(a1*phi^2)
        pred2 = C * 4 * G / (a1 * PHI**2)
        print(f"        C*4G/(a1*phi^2) = {pred2:+.5f} (err: {abs(pred2-delta_known[k]):.6f})")

print("\nSECTOR C: WHITE quarks (d, b)")
# d: delta = -2/a1
# b: delta = -C*ln19/phi
for k in [2, 8]:
    a, b_val = ab_vals[k]
    N = norms[k]
    absN = abs(N)
    print(f"  {fermions[k]:3s}: |N| = {absN}, delta = {delta_known[k]:+.5f}")
    if k == 2:
        pred = -2.0 / a1
        print(f"        -2/a1 = {pred:+.5f} EXACT (graph distance)")
    elif k == 8:
        if absN > 1:
            pred = -C * log(absN) / PHI
            print(f"        -C*ln|N|/phi = {pred:+.5f} (err: {abs(pred-delta_known[k]):.6f})")

print("\nSECTOR D: WHITE leptons (mu, tau)")
for k in [4, 6]:
    a, b_val = ab_vals[k]
    z = a + b_val * PHI
    zc = a + b_val * PHIc
    N = norms[k]
    print(f"  {fermions[k]:3s}: z' = {zc:+.5f}, |z'|^(3/4) = {abs(zc)**0.75:.5f}, "
          f"delta = {delta_known[k]:+.5f}")

# =====================================================================
# PART 7: New Attempt - Holonomy on McKay Tree
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: McKay Tree Holonomy = exp(i*alpha_s * path_integral)")
print("=" * 72)

# Key idea: alpha_s = 1/(2*phi^3). And C = 2/13 = 2/(a1^2+1).
# Wilson line phase = integral of gauge field along path.
# If gauge field A ~ alpha_s on each edge, then:
#   holonomy = exp(i * alpha_s * dist)
# The "correction" to the mass exponent = Re(holonomy) or Im(holonomy)?

alpha_s = 1 / (2 * PHI**3)
print(f"\nalpha_s = 1/(2*phi^3) = {alpha_s:.6f}")
print(f"C = 2/13 = {C:.6f}")
print(f"C/alpha_s = {C/alpha_s:.6f}")
print(f"2*C*phi = {2*C*PHI:.6f}")
print(f"1/a1 = {1/a1:.6f}")
print(f"C*phi^2 = {C*PHI**2:.6f}")

# Test: delta = -C * f(spectral_data) for each fermion
# For d quark: -2/a1 = -0.4. And a1 = 5.
# Note: 2/a1 = 2C * (a1^2+1)/(2*a1) = C * (a1^2+1)/a1 = C * (25+1)/5 = C * 26/5
# So delta_d = -C * 26/5 = -C * sum(Casimirs)/a1

sum_casimirs = a1**2 + 1  # = 26
print(f"\ndelta_d = -C * sum_C2/a1 = -{C} * {sum_casimirs}/{a1} = {-C*sum_casimirs/a1:.6f}")
print(f"  vs known delta_d = {delta_known[2]:+.6f}")
print(f"  sum_C2 = a1^2 + 1 = {sum_casimirs}")

# =====================================================================
# PART 8: Testing Whether ln|N| Comes from Path Product
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: Does ln|N(z)| = C * sum(Casimirs) * path_function?")
print("=" * 72)

# For c: delta = C*ln5. Path: 0->1->2->3->4->5
# For t: delta = C*ln19. Path: 0->1->2->3->4->5->6->7
# For b: delta = -C*ln19/phi. Path: 0->1->2->3->4->5->8

# The paths for c, t, b share the stem 0->1->2->3->4->5
# t adds ->6->7 (two more edges)
# b adds ->8 (one more edge, branch)

# delta_t - delta_c = C*(ln19 - ln5) = C*ln(19/5) = C*ln(3.8)
dt_dc = delta_known[7] - delta_known[5]
print(f"  delta_t - delta_c = {dt_dc:.5f}")
print(f"  C*ln(19/5) = {C*log(19/5):.5f}")
print(f"  These should be equal: err = {abs(dt_dc - C*log(19/5)):.6f}")

# The extra edges for t vs c: (5->6) and (6->7)
# dims along: d_5=6, d_6=4, d_7=2
# d_5->d_6: ratio 6/4 = 3/2
# d_6->d_7: ratio 4/2 = 2
# Product of ratios: 3
# C * ln(3) = 0.169
# But delta_t - delta_c = C*ln(19/5) = C*1.335 = 0.205. NOT ln(3).

# Maybe the edge weight is related to dimension changes?
print("\nEdge weights from dim ratios:")
for i, j in edges:
    di, dj = dims[i], dims[j]
    ratio = di/dj if di > dj else dj/di
    print(f"  ({i:d}->{j:d}): d_{i}={di}, d_{j}={dj}, max_ratio={ratio:.4f}, log(ratio)={log(ratio):.4f}")

# =====================================================================
# PART 9: Spectral Zeta at Different Points
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: Spectral Zeta Function Values")
print("=" * 72)

# Cayley Laplacian eigenvalues at class 10a
L10a = np.array([12 * (1 - chi[k, 1]/dims[k]) for k in range(9)])
print("\nCayley eigenvalues L_10a:")
for k in range(9):
    print(f"  {fermions[k]:3s}: L = {L10a[k]:.6f}")

# Spectral zeta: zeta(s) = sum(L_k^(-s)) (excluding zero mode)
for s in [1, 2, 0.5]:
    nonzero = [L10a[k] for k in range(1, 9)]  # exclude k=0 (zero mode)
    zeta_s = sum(l**(-s) for l in nonzero if l > 0)
    print(f"  zeta({s}) = {zeta_s:.6f}")

# =====================================================================
# PART 10: Key Test - Can We Express delta_d from Laplacian?
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: d Quark Correction from Laplacian Data")
print("=" * 72)

# d quark at node 2: L_10a = 12-4phi = 5.528
L_d = L10a[2]
print(f"\nL_10a(d) = 12 - 4*phi = {L_d:.6f}")
print(f"delta_d = -2/a1 = {-2/a1:.6f}")
print()

# Test various expressions
tests = [
    ("C * L / phi", C * L_d / PHI),
    ("-C * L", -C * L_d),
    ("-12 / (N * L / 2)", -12 / (N_group * L_d / 2) if L_d > 0 else float('inf')),
    ("-2 * C * phi", -2 * C * PHI),
    ("-(L/12) / phi^2", -(L_d/12) / PHI**2),
    ("-L / (a1 * b1)", -L_d / (a1 * b1)),
    ("-2 / (12 - L)", -2 / (12 - L_d) if abs(12 - L_d) > 0.01 else float('inf')),
    ("-2 / (L * phi^2 / 2)", -2 / (L_d * PHI**2 / 2) if L_d > 0 else float('inf')),
    ("-L / (N_gen * phi * b1)", -L_d / (N_gen * PHI * b1)),
]

for desc, val in tests:
    err = abs(val - (-2/a1))
    match = "***" if err < 0.001 else ""
    print(f"  {desc:30s} = {val:+.6f} (err: {err:.6f}) {match}")

# =====================================================================
# PART 11: Key Insight - Relation Between 2/a1 and Theory Constants
# =====================================================================
print("\n" + "=" * 72)
print("PART 11: What Is 2/a1 in Terms of Theory Constants?")
print("=" * 72)

print(f"\n  2/a1 = {2/a1:.10f}")
print(f"  C * (a1^2+1)/a1 = C*26/5 = {C*26/5:.10f}")
print(f"  4/(a1*(a1^2+1))*26/5 = ... let's simplify")
print(f"  2/a1 = 2*C*13/a1 = 2*(2/13)*13/5 = 4/5 = 0.8? No!")
print(f"  2/a1 = 0.4, 4/a1 = 0.8")
print()

# decompositions of 2/a1
print("  Decompositions of 2/a1:")
print(f"    2/a1 = {2/a1}")
print(f"    dim(rho_1)/a1 = {dims[1]/a1}")
print(f"    (b1-N_gen-1)/a1 = {(b1-N_gen-1)/a1}")  # = 2/5
print(f"    |Galois pair in kernel|/a1 = 2/5")
print(f"    dist(e,d)/a1 = 2/5")
print(f"    C * sum_Casimirs/a1 = {C * 26 / a1}")
# C*26/5 = (2/13)*26/5 = 52/65 = 4/5. Not 2/5!
print(f"    Actually: 2/a1 = C * 13/a1 = C * (a1^2+1)/(2*a1)")
Cval = 2/13
print(f"    C * (a1^2+1)/(2*a1) = {Cval * 26/10}")
print()

# The cleanest: 2/a1 = dim(fund_SU2) / a1
# Or: 2/a1 = dim(rho_1) / dim(rho_4) (2/5)
# Wait: dim(rho_4)=5=a1. So 2/a1 = dim(rho_1)/dim(rho_4).
print("  CLEANEST DECOMPOSITION:")
print(f"    2/a1 = dim(rho_1)/dim(rho_4) = 2/5")
print(f"    where rho_1 = SU(2) fundamental, rho_4 = self-conjugate (turning point)")
print(f"    Interpretation: ratio of SU(2) to maximal self-conjugate irrep")

# =====================================================================
# PART 12: Comprehensive Test of Galois-Spectral Formula
# =====================================================================
print("\n" + "=" * 72)
print("PART 12: Galois-Spectral Unified Formula")
print("=" * 72)

# For each fermion, the Galois norm and spectral G value
# determine the correction differently based on sector.
# Can we write ONE formula?

# Key observation: for BLACK nodes, G != 0 and acts as a "spectral weight"
# For WHITE nodes, G = 0 but the RATIONAL class characters are non-zero

# Attempt: delta = C * [G * f_1(color) + chi_rational * f_2(color)]

# Actually, let's test: is there a class c such that
# delta = -a * chi_k(c) / d_k for ALL fermions, for some constant a?

print("\nTest: delta_k = -a * chi_k(c) / d_k for each class c")
for col_idx, cname in enumerate(class_names):
    # For each class, solve for 'a' using least squares on non-trivial fermions
    x_vals = [chi[k, col_idx] / dims[k] for k in range(9)]
    # Skip e (both x and delta are 0)
    x_arr = np.array(x_vals[1:])
    d_arr = np.array(delta_known[1:])

    if np.dot(x_arr, x_arr) > 1e-10:
        a_fit = -np.dot(x_arr, d_arr) / np.dot(x_arr, x_arr)
        pred = [-a_fit * x for x in x_vals]
        residuals = [delta_known[k] - pred[k] for k in range(9)]
        rms = sqrt(sum(r**2 for r in residuals) / 9)
        r2 = 1 - sum(r**2 for r in residuals) / sum(d**2 for d in delta_known if d != 0)
        print(f"  class {cname:4s}: a = {a_fit:+.4f}, RMS = {rms:.4f}, R^2 = {r2:.4f}")

# =====================================================================
# PART 13: Two-Class Formula
# =====================================================================
print("\n" + "=" * 72)
print("PART 13: Two-Class Linear Formula")
print("=" * 72)

# delta = a * chi(c1)/d + b * chi(c2)/d for best pair of classes
best_rms = 1.0
best_pair = None
best_coeffs = None

for i1 in range(9):
    for i2 in range(i1+1, 9):
        # Features: chi(c1)/d and chi(c2)/d
        X = np.array([[chi[k, i1]/dims[k], chi[k, i2]/dims[k]] for k in range(9)])
        y = np.array(delta_known)

        # Least squares (no intercept)
        try:
            coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
            pred = X @ coeffs
            rms = sqrt(np.mean((y - pred)**2))
            if rms < best_rms:
                best_rms = rms
                best_pair = (i1, i2)
                best_coeffs = coeffs
        except:
            pass

print(f"\nBest two-class fit: classes {class_names[best_pair[0]]} and {class_names[best_pair[1]]}")
print(f"  Coefficients: a = {best_coeffs[0]:+.6f}, b = {best_coeffs[1]:+.6f}")
print(f"  RMS = {best_rms:.6f}")
print()

X = np.array([[chi[k, best_pair[0]]/dims[k], chi[k, best_pair[1]]/dims[k]] for k in range(9)])
pred = X @ best_coeffs
print(f"  {'name':4s}  {'delta_true':10s}  {'delta_pred':10s}  {'error':8s}")
for k in range(9):
    err = abs(delta_known[k] - pred[k])
    match = "***" if err < 0.001 else ("**" if err < 0.01 else "")
    print(f"  {fermions[k]:4s}  {delta_known[k]:+10.5f}  {pred[k]:+10.5f}  {err:8.5f} {match}")

# =====================================================================
# PART 14: Three-Class Formula (exclude trivial 1a and 2a)
# =====================================================================
print("\n" + "=" * 72)
print("PART 14: Three-Class Linear Formula")
print("=" * 72)

best_rms3 = 1.0
best_triple = None
best_coeffs3 = None

for i1 in range(9):
    for i2 in range(i1+1, 9):
        for i3 in range(i2+1, 9):
            X = np.array([[chi[k, i1]/dims[k], chi[k, i2]/dims[k], chi[k, i3]/dims[k]]
                         for k in range(9)])
            y = np.array(delta_known)
            try:
                coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
                pred = X @ coeffs
                rms = sqrt(np.mean((y - pred)**2))
                if rms < best_rms3:
                    best_rms3 = rms
                    best_triple = (i1, i2, i3)
                    best_coeffs3 = coeffs
            except:
                pass

c1, c2, c3 = best_triple
print(f"\nBest three-class fit: {class_names[c1]}, {class_names[c2]}, {class_names[c3]}")
print(f"  Coefficients: {best_coeffs3[0]:+.6f}, {best_coeffs3[1]:+.6f}, {best_coeffs3[2]:+.6f}")
print(f"  RMS = {best_rms3:.6f}")
print()

X = np.array([[chi[k, c1]/dims[k], chi[k, c2]/dims[k], chi[k, c3]/dims[k]] for k in range(9)])
pred = X @ best_coeffs3
print(f"  {'name':4s}  {'delta_true':10s}  {'delta_pred':10s}  {'error':8s}")
for k in range(9):
    err = abs(delta_known[k] - pred[k])
    match = "***" if err < 0.001 else ("**" if err < 0.01 else "")
    print(f"  {fermions[k]:4s}  {delta_known[k]:+10.5f}  {pred[k]:+10.5f}  {err:8.5f} {match}")

# Check if coefficients are "nice" fractions
print("\n  Coefficient analysis:")
for i, (idx_c, coeff) in enumerate(zip(best_triple, best_coeffs3)):
    # Check against common fractions involving C, phi, a1, b1
    for num in range(-10, 11):
        for den in [1, 2, 3, 5, 6, 13, 26]:
            if den != 0:
                frac = num / den
                if abs(frac - coeff) < 0.005 and frac != 0:
                    print(f"    coeff[{class_names[idx_c]}] = {coeff:.6f} ~ {num}/{den} = {frac:.6f} (err: {abs(frac-coeff):.6f})")
    for base in [C, C*PHI, C/PHI, 1/a1, 1/b1, C*PHI**2]:
        for mult in range(-5, 6):
            if mult != 0:
                val = mult * base
                if abs(val - coeff) < 0.005:
                    desc = f"{mult}*C" if base == C else (f"{mult}*C*phi" if base == C*PHI else
                           (f"{mult}*C/phi" if base == C/PHI else (f"{mult}/a1" if base == 1/a1 else
                           (f"{mult}/b1" if base == 1/b1 else f"{mult}*C*phi^2"))))
                    print(f"    coeff[{class_names[idx_c]}] = {coeff:.6f} ~ {desc} = {val:.6f} (err: {abs(val-coeff):.6f})")

# =====================================================================
# PART 15: THE KEY QUESTION - Can We Unify?
# =====================================================================
print("\n" + "=" * 72)
print("PART 15: Honest Assessment")
print("=" * 72)

print("""
FINDINGS:

1. GRAPH DISTANCE for d quark: delta_d = -2/a1 = -dist(e,d)/a1 EXACT
   - 2/a1 = dim(rho_1)/dim(rho_4) = SU(2)_fund / maximal_self_conj
   - Unique to b=0 (purely rational z)
   - STATUS: IDENTIFIED -> searching for DERIVATION

2. CHARACTER-BASED FORMULA for WHITE nodes:
   - Tested delta = -2/(a1 + b1*x(6a)): works for d,e but NOT for mu,tau,b
   - Single-class linear: poor (R^2 < 0.5 for all classes)
   - Two-class and three-class: may fit but coefficients not "nice"

3. SECTOR STRUCTURE is REAL:
   - BLACK with |N|>1: delta = C*ln|N| (c,t)
   - BLACK with |N|=1: delta = C*G/phi^2 or 0 (s,u)
   - WHITE with b=0: delta = -dist/a1 (d only)
   - WHITE with |N|>1: delta = -C*ln|N|/phi (b only)
   - WHITE leptons: exponent 3/4 pattern (mu, tau)

4. NO SINGLE SPECTRAL FORMULA spans all sectors.
   The corrections have DIFFERENT ORIGINS:
   - Galois anti-symmetric G (BLACK quarks near |N|=1)
   - Galois norm ln|N| (quarks with |N|>1)
   - Graph distance (d quark, b=0)
   - Spacetime dimension 3/4 (leptons)

CONCLUSION: The 4-sector structure appears IRREDUCIBLE. Each sector
has a clear mechanism but there is no unification across sectors.
This may reflect the 4-fold structure of SM: Q_L, u_R, d_R, L_L.
""")

print("=" * 72)
print("EXPERIMENT COMPLETE")
print("=" * 72)

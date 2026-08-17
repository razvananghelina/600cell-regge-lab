"""
exp470c: E8 -> 600-cell with CORRECT simple roots
===================================================
Using Bourbaki E8 simple roots.
STATUS: EXPLORATORY
"""

import numpy as np
from itertools import product as iprod
from collections import Counter
import random

phi = (1 + np.sqrt(5)) / 2
phi_c = (1 - np.sqrt(5)) / 2

# ============================================================
# E8 roots (standard)
# ============================================================
E8_list = []
for i in range(8):
    for j in range(i+1, 8):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = np.zeros(8); v[i] = si; v[j] = sj
                E8_list.append(v)
for signs in iprod([1, -1], repeat=8):
    if sum(1 for s in signs if s == -1) % 2 == 0:
        E8_list.append(np.array(signs) / 2.0)
E8 = np.array(E8_list)
print(f"E8 roots: {len(E8)}")

# ============================================================
# Correct E8 simple roots (Bourbaki)
# Dynkin diagram:
#     a2
#     |
# a1--a3--a4--a5--a6--a7--a8
# ============================================================
simple = np.zeros((8, 8))
simple[0] = np.array([1, -1, -1, -1, -1, -1, -1, 1]) / 2.0  # alpha_1
simple[1] = np.array([1, 1, 0, 0, 0, 0, 0, 0])                # alpha_2
simple[2] = np.array([-1, 1, 0, 0, 0, 0, 0, 0])               # alpha_3
simple[3] = np.array([0, -1, 1, 0, 0, 0, 0, 0])               # alpha_4
simple[4] = np.array([0, 0, -1, 1, 0, 0, 0, 0])               # alpha_5
simple[5] = np.array([0, 0, 0, -1, 1, 0, 0, 0])               # alpha_6
simple[6] = np.array([0, 0, 0, 0, -1, 1, 0, 0])               # alpha_7
simple[7] = np.array([0, 0, 0, 0, 0, -1, 1, 0])               # alpha_8

# Verify Cartan matrix
cartan = np.zeros((8, 8), dtype=int)
for i in range(8):
    for j in range(8):
        cartan[i,j] = int(round(2 * simple[i] @ simple[j] / (simple[i] @ simple[i])))

print("Cartan matrix:")
print(cartan)

# Verify: all diagonal = 2, off-diag <= 0
assert all(cartan[i,i] == 2 for i in range(8))
assert all(cartan[i,j] <= 0 for i in range(8) for j in range(8) if i != j)
det_C = round(np.linalg.det(cartan.astype(float)))
print(f"det(Cartan) = {det_C} (should be 1 for E8)")

# Dynkin diagram edges
edges = [(i,j) for i in range(8) for j in range(i+1,8) if cartan[i,j] < 0]
print(f"Dynkin edges: {edges}")

# ============================================================
# Coxeter element
# ============================================================
def make_refl_matrix(alpha):
    """Reflection matrix for root alpha."""
    return np.eye(8) - 2 * np.outer(alpha, alpha) / (alpha @ alpha)

# Standard ordering: all simple reflections in order 1,2,...,8
C = np.eye(8)
for a in simple:
    C = make_refl_matrix(a) @ C

eigvals_C, eigvecs_C = np.linalg.eig(C)
args_C = np.angle(eigvals_C)

# Determine exponents
exponents = []
for i in range(8):
    m = round(args_C[i] * 30 / (2*np.pi))
    if m <= 0: m += 30
    exponents.append((m, i))

exponents.sort()
print(f"\nCoxeter exponents: {[m for m, _ in exponents]}")
print(f"Expected: [1, 7, 11, 13, 17, 19, 23, 29]")

C30 = np.linalg.matrix_power(C, 30)
print(f"C^30 = I? {np.allclose(C30, np.eye(8))}")

# ============================================================
# H4 projection
# ============================================================
# H4 exponents: {1, 11, 19, 29}
# H4' exponents: {7, 13, 17, 23}

h4_idx = [i for m, i in exponents if m in {1, 11, 19, 29}]
h4p_idx = [i for m, i in exponents if m in {7, 13, 17, 23}]

def build_real_basis(eigvecs, indices):
    """Extract real orthonormal basis from complex eigenvectors."""
    raw_basis = []
    used = set()
    for j in indices:
        if j in used:
            continue
        v = eigvecs[:, j]
        for k in indices:
            if k > j and k not in used:
                if np.allclose(v, np.conj(eigvecs[:, k]), atol=1e-10):
                    raw_basis.append(np.real(v))
                    raw_basis.append(np.imag(v))
                    used.add(j); used.add(k)
                    break
        else:
            if j not in used:
                if np.max(np.abs(np.imag(v))) < 1e-10:
                    raw_basis.append(np.real(v))
                else:
                    raw_basis.append(np.real(v))
                    raw_basis.append(np.imag(v))
                used.add(j)

    B = np.array(raw_basis).T  # 8 x n
    Q, R = np.linalg.qr(B)
    return Q[:, :4]  # 8 x 4

P_h4 = build_real_basis(eigvecs_C, h4_idx)
P_h4p = build_real_basis(eigvecs_C, h4p_idx)

print(f"\nH4 projector: {P_h4.shape}, ortho? {np.allclose(P_h4.T @ P_h4, np.eye(4))}")
print(f"H4' projector: {P_h4p.shape}, ortho? {np.allclose(P_h4p.T @ P_h4p, np.eye(4))}")
print(f"H4 perp H4'? {np.allclose(P_h4.T @ P_h4p, np.zeros((4,4)), atol=1e-10)}")

# Project E8 roots
proj = E8 @ P_h4          # 240 x 4
proj_p = E8 @ P_h4p       # 240 x 4

norms = np.sqrt(np.sum(proj**2, axis=1))
norms_p = np.sqrt(np.sum(proj_p**2, axis=1))

# Verify: |proj|^2 + |proj'|^2 = 2
total = np.sum(proj**2, axis=1) + np.sum(proj_p**2, axis=1)
print(f"\n|pi(r)|^2 + |pi'(r)|^2 = {sorted(set(np.round(total, 6)))}")

# Norm distribution
print(f"\nH4 projection norm distribution:")
norm_counts = Counter(np.round(norms, 5))
for n, c in sorted(norm_counts.items()):
    print(f"  |pi(r)| = {n:.5f}: {c} roots")

n_distinct = len(norm_counts)
print(f"Number of distinct projected norms: {n_distinct}")

if n_distinct == 2:
    print("PERFECT: 120 + 120 split into two norm shells!")
elif n_distinct > 2:
    print(f"Got {n_distinct} shells, not 2. Checking if grouping gives 120+120...")
    # Maybe there are small numerical errors splitting what should be 2 groups
    norm_list = sorted(norm_counts.keys())
    print(f"Norm values: {[f'{n:.6f}' for n in norm_list]}")

    # Try to identify 2 clusters
    from scipy.cluster.hierarchy import fcluster, linkage
    norms_2d = norms.reshape(-1, 1)
    Z = linkage(norms_2d, method='ward')
    clusters = fcluster(Z, t=2, criterion='maxclust')
    cluster_sizes = Counter(clusters)
    print(f"2-cluster sizes: {dict(cluster_sizes)}")

    if set(cluster_sizes.values()) == {120}:
        print("YES: 120 + 120 after clustering!")
        group1 = [i for i in range(240) if clusters[i] == 1]
        group2 = [i for i in range(240) if clusters[i] == 2]
    else:
        # Fall back: sort by norm and split at 120
        sorted_idx = np.argsort(norms)
        group1 = list(sorted_idx[:120])
        group2 = list(sorted_idx[120:])
        print(f"Splitting at median: sizes {len(group1)}, {len(group2)}")
        print(f"  Group 1 norm range: [{norms[group1[0]]:.6f}, {norms[group1[-1]]:.6f}]")
        print(f"  Group 2 norm range: [{norms[group2[0]]:.6f}, {norms[group2[-1]]:.6f}]")
else:
    print(f"Only {n_distinct} norm -- all roots project to same radius!")
    group1 = list(range(120))
    group2 = list(range(120, 240))

# ============================================================
# Verify 600-cell structure in projected groups
# ============================================================
print(f"\n{'='*70}")
print("VERIFYING 600-CELL STRUCTURE")
print(f"{'='*70}")

def check_600cell(indices, label):
    """Check if a subset of 120 projected roots forms a 600-cell."""
    subset = proj[indices]
    r = np.sqrt(np.sum(subset**2, axis=1))

    # Normalize to unit vectors
    r_mean = np.mean(r)
    subset_n = subset / r.reshape(-1, 1)

    IP = subset_n @ subset_n.T
    ip_vals = Counter(np.round(IP[np.triu_indices(120, k=1)], 3))

    print(f"\n{label} ({len(indices)} roots, mean projected norm = {r_mean:.6f}):")
    print(f"  Norm range: [{r.min():.6f}, {r.max():.6f}]")
    print(f"  Inner products (normalized):")
    for v, c in sorted(ip_vals.items()):
        print(f"    {v:+.3f}: {c}")

    # 600-cell signature: 8 distinct inner product values
    # -1 (60), -0.809 (720), -0.5 (1200), -0.309 (720),
    # 0 (1800), 0.309 (720), 0.5 (1200), 0.809 (720)
    expected = {-1.0: 60, -0.809: 720, -0.5: 1200, -0.309: 720,
                0.0: 1800, 0.309: 720, 0.5: 1200, 0.809: 720}

    match = True
    for exp_v, exp_c in expected.items():
        found = False
        for v, c in ip_vals.items():
            if abs(v - exp_v) < 0.01 and c == exp_c:
                found = True
                break
        if not found:
            match = False

    print(f"  Matches 600-cell? {match}")
    return match

if n_distinct == 2:
    norm_vals = sorted(norm_counts.keys())
    group1 = [i for i in range(240) if abs(norms[i] - norm_vals[0]) < 0.001]
    group2 = [i for i in range(240) if abs(norms[i] - norm_vals[1]) < 0.001]
    r = norm_vals[1] / norm_vals[0]
    print(f"Norm ratio: {r:.6f}, phi = {phi:.6f}, match? {abs(r - phi) < 0.01}")
    is_600_1 = check_600cell(group1, "Inner 600-cell")
    is_600_2 = check_600cell(group2, "Outer 600-cell")
else:
    is_600_1 = check_600cell(group1, "Group 1")
    is_600_2 = check_600cell(group2, "Group 2")

# ============================================================
# Alternative: Try BIPARTITE Coxeter element
# ============================================================
print(f"\n{'='*70}")
print("TRYING BIPARTITE COXETER ELEMENT")
print(f"{'='*70}")

# The standard Coxeter element product s1*s2*...*s8 uses arbitrary ordering.
# For the H4/H4' decomposition, the BIPARTITE ordering is standard.
# In a bipartite ordering, we first apply all reflections from one color
# class, then all from the other.

# 2-color the Dynkin diagram
color = [-1] * 8
color[0] = 0
queue = [0]
while queue:
    node = queue.pop(0)
    for i, j in edges:
        nbr = j if i == node else (i if j == node else None)
        if nbr is not None and color[nbr] == -1:
            color[nbr] = 1 - color[node]
            queue.append(nbr)

c0 = [i for i in range(8) if color[i] == 0]
c1 = [i for i in range(8) if color[i] == 1]
print(f"2-coloring: {color}")
print(f"Color 0: {c0}, Color 1: {c1}")

# Bipartite Coxeter element: color 0 first, then color 1
C_bi = np.eye(8)
for i in c0:
    C_bi = make_refl_matrix(simple[i]) @ C_bi
for i in c1:
    C_bi = make_refl_matrix(simple[i]) @ C_bi

eigvals_bi, eigvecs_bi = np.linalg.eig(C_bi)
args_bi = np.angle(eigvals_bi)

print(f"\nBipartite Coxeter eigenvalues:")
exp_bi = []
for i in range(8):
    m = round(args_bi[i] * 30 / (2*np.pi))
    if m <= 0: m += 30
    exp_bi.append((m, i))
    print(f"  m = {m}: |lambda| = {abs(eigvals_bi[i]):.8f}")

exp_bi.sort()
print(f"Exponents: {[m for m,_ in exp_bi]}")
print(f"C_bi^30 = I? {np.allclose(np.linalg.matrix_power(C_bi, 30), np.eye(8))}")

# H4 projection from bipartite element
h4_bi = [i for m, i in exp_bi if m in {1, 11, 19, 29}]
P_bi = build_real_basis(eigvecs_bi, h4_bi)

proj_bi = E8 @ P_bi
norms_bi = np.sqrt(np.sum(proj_bi**2, axis=1))

print(f"\nBipartite projection norm distribution:")
for n, c in sorted(Counter(np.round(norms_bi, 5)).items()):
    print(f"  |pi(r)| = {n:.5f}: {c}")

# ============================================================
# Try with the STANDARD ordering used in the literature
# ============================================================
print(f"\n{'='*70}")
print("LITERATURE STANDARD: Moody-Patera projection")
print(f"{'='*70}")

# The Moody-Patera projection uses a specific 4x8 matrix.
# From Dechant (2013), the projection matrix is built from
# the eigenvectors of the Coxeter element in the "standard" basis.

# Actually, many sources use a different E8 root convention.
# Let me try the "D8+" convention more carefully.

# E8 simple roots in Bourbaki numbering with branching at node 4:
#
#  1 - 3 - 4 - 5 - 6 - 7 - 8
#          |
#          2
#
# With the convention:
# alpha_i = e_i - e_{i+1} for i = 1, ..., 6
# alpha_7 = e_6 + e_7
# alpha_8 = -(1/2)(sum of all e_i)

# Actually, I realize my numbering might not match the standard Bourbaki.
# Let me just use a KNOWN correct set.

# From Humphreys "Reflection Groups and Coxeter Groups":
# For E8, using coordinates (e_1, ..., e_8) with the root system:
# {+-e_i +- e_j : i < j} union {(1/2)(+-e_1 +- ... +- e_8) : even # of minus}

# Simple roots (Bourbaki, Table VI):
# alpha_1 = (1/2)(e_1 - e_2 - e_3 - e_4 - e_5 - e_6 - e_7 + e_8)
# alpha_2 = e_1 + e_2
# alpha_3 = e_2 - e_1
# alpha_4 = e_3 - e_2
# alpha_5 = e_4 - e_3
# alpha_6 = e_5 - e_4
# alpha_7 = e_6 - e_5
# alpha_8 = e_7 - e_6

# This is what I had at the start! Let me verify the Cartan matrix again.
s2 = np.zeros((8, 8))
s2[0] = np.array([1, -1, -1, -1, -1, -1, -1, 1]) / 2.0
s2[1] = np.array([1, 1, 0, 0, 0, 0, 0, 0])
s2[2] = np.array([-1, 1, 0, 0, 0, 0, 0, 0])
s2[3] = np.array([0, -1, 1, 0, 0, 0, 0, 0])
s2[4] = np.array([0, 0, -1, 1, 0, 0, 0, 0])
s2[5] = np.array([0, 0, 0, -1, 1, 0, 0, 0])
s2[6] = np.array([0, 0, 0, 0, -1, 1, 0, 0])
s2[7] = np.array([0, 0, 0, 0, 0, -1, 1, 0])

cart2 = np.zeros((8, 8), dtype=int)
for i in range(8):
    for j in range(8):
        cart2[i,j] = int(round(2 * s2[i] @ s2[j] / (s2[i] @ s2[i])))

print("Bourbaki Cartan matrix:")
print(cart2)
edges2 = [(i,j) for i in range(8) for j in range(i+1,8) if cart2[i,j] < 0]
print(f"Edges: {edges2}")
print(f"det = {round(np.linalg.det(cart2.astype(float)))}")

# Build Coxeter element from these
C2 = np.eye(8)
for a in s2:
    C2 = make_refl_matrix(a) @ C2

eigvals_C2, eigvecs_C2 = np.linalg.eig(C2)
args_C2 = np.angle(eigvals_C2)

exp2 = []
for i in range(8):
    m = round(args_C2[i] * 30 / (2*np.pi))
    if m <= 0: m += 30
    exp2.append((m, i))

exp2.sort()
print(f"\nCoxeter exponents: {[m for m,_ in exp2]}")
print(f"C^30 = I? {np.allclose(np.linalg.matrix_power(C2, 30), np.eye(8))}")

# H4 projection
h4_2 = [i for m, i in exp2 if m in {1, 11, 19, 29}]
P2 = build_real_basis(eigvecs_C2, h4_2)
proj2 = E8 @ P2
norms2 = np.sqrt(np.sum(proj2**2, axis=1))

print(f"\nBourbaki projection norm distribution:")
for n, c in sorted(Counter(np.round(norms2, 5)).items()):
    print(f"  |pi(r)| = {n:.5f}: {c}")

# ============================================================
# Now use the Bourbaki roots for bipartite ordering
# ============================================================
color2 = [-1]*8
color2[0] = 0
queue2 = [0]
while queue2:
    node = queue2.pop(0)
    for i, j in edges2:
        nbr = j if i == node else (i if j == node else None)
        if nbr is not None and color2[nbr] == -1:
            color2[nbr] = 1 - color2[node]
            queue2.append(nbr)

c0_2 = [i for i in range(8) if color2[i] == 0]
c1_2 = [i for i in range(8) if color2[i] == 1]
print(f"\n2-coloring: {color2}")

# Bipartite Coxeter from Bourbaki roots
C_bi2 = np.eye(8)
for i in c0_2:
    C_bi2 = make_refl_matrix(s2[i]) @ C_bi2
for i in c1_2:
    C_bi2 = make_refl_matrix(s2[i]) @ C_bi2

eigvals_bi2, eigvecs_bi2 = np.linalg.eig(C_bi2)
args_bi2 = np.angle(eigvals_bi2)

exp_bi2 = []
for i in range(8):
    m = round(args_bi2[i] * 30 / (2*np.pi))
    if m <= 0: m += 30
    exp_bi2.append((m, i))

exp_bi2.sort()
print(f"Bipartite Coxeter exponents: {[m for m,_ in exp_bi2]}")

h4_bi2 = [i for m, i in exp_bi2 if m in {1, 11, 19, 29}]
P_bi2 = build_real_basis(eigvecs_bi2, h4_bi2)
proj_bi2 = E8 @ P_bi2
norms_bi2 = np.sqrt(np.sum(proj_bi2**2, axis=1))

print(f"\nBipartite Bourbaki projection norms:")
nc_bi2 = Counter(np.round(norms_bi2, 5))
for n, c in sorted(nc_bi2.items()):
    print(f"  |pi(r)| = {n:.5f}: {c}")

# Check if exactly 2 groups of 120
if len(nc_bi2) == 2:
    nv = sorted(nc_bi2.keys())
    print(f"\n*** SUCCESS: Two groups of 120! ***")
    print(f"Inner radius: {nv[0]:.6f}")
    print(f"Outer radius: {nv[1]:.6f}")
    print(f"Ratio: {nv[1]/nv[0]:.6f}")
    print(f"phi = {phi:.6f}")
    print(f"Match? {abs(nv[1]/nv[0] - phi) < 0.01}")

    # Verify 600-cell structure
    for k, norm_val in enumerate(nv):
        idx = [i for i in range(240) if abs(norms_bi2[i] - norm_val) < 0.001]
        subset = proj_bi2[idx]
        subset_n = subset / np.linalg.norm(subset, axis=1, keepdims=True)
        IP = subset_n @ subset_n.T
        ip_c = Counter(np.round(IP[np.triu_indices(120, k=1)], 3))
        label = "Inner" if k == 0 else "Outer"
        print(f"\n  {label} 600-cell inner products:")
        for v, c in sorted(ip_c.items()):
            print(f"    {v:+.3f}: {c}")

# ============================================================
# If neither standard nor bipartite gives 120+120,
# the issue is likely the QR decomposition of eigenvectors.
# Let me try an explicit projection.
# ============================================================
print(f"\n{'='*70}")
print("EXPLICIT PROJECTION (Dechant method)")
print(f"{'='*70}")

# Following Dechant (2013): the 4x8 projection matrix from E8 to H4.
# The idea: use the eigenvectors of the Cartan matrix (NOT the Coxeter element).

# Actually, the standard approach uses the Coxeter element eigenvalues.
# For each conjugate pair (exp(2*pi*i*m/30), exp(-2*pi*i*m/30)),
# the real 2D eigenspace is spanned by:
#   v_m = Re(eigenvector), w_m = Im(eigenvector)

# H4 plane = span of (v_1, w_1, v_11, w_11)
# H4' plane = span of (v_7, w_7, v_13, w_13)

# Let me extract these more carefully
print("Re-extracting eigenspaces carefully...")

# Use Bourbaki Coxeter element C2
eigvals_C2, eigvecs_C2 = np.linalg.eig(C2)
args_C2 = np.angle(eigvals_C2)

# Find pairs
print("Eigenvalues and exponents:")
for i in range(8):
    m = args_C2[i] * 30 / (2*np.pi)
    print(f"  i={i}: arg/(2pi/30) = {m:.4f}, |lambda| = {abs(eigvals_C2[i]):.8f}")

# Group into conjugate pairs
pairs = {}
used = set()
for i in range(8):
    if i in used:
        continue
    for j in range(i+1, 8):
        if j in used:
            continue
        if abs(eigvals_C2[i] - np.conj(eigvals_C2[j])) < 1e-8:
            m = round(abs(args_C2[i]) * 30 / (2*np.pi))
            pairs[m] = (i, j)
            used.add(i); used.add(j)
            break

print(f"Conjugate pairs (by |m|): {pairs}")

# Build H4 basis: exponents 1 and 11
# H4' basis: exponents 7 and 13
h4_planes = []
h4p_planes = []

for m, (i, j) in sorted(pairs.items()):
    v_re = np.real(eigvecs_C2[:, i])
    v_im = np.imag(eigvecs_C2[:, i])
    if m in {1, 11}:
        h4_planes.extend([v_re, v_im])
    elif m in {7, 13}:
        h4p_planes.extend([v_re, v_im])

h4_basis = np.array(h4_planes).T  # 8 x 4
h4p_basis = np.array(h4p_planes).T

# Orthonormalize
Q1, _ = np.linalg.qr(h4_basis)
Q2, _ = np.linalg.qr(h4p_basis)

P_h4_new = Q1[:, :4]
P_h4p_new = Q2[:, :4]

print(f"\nNew H4 projector: ortho? {np.allclose(P_h4_new.T @ P_h4_new, np.eye(4))}")
print(f"H4 perp H4'? {np.allclose(P_h4_new.T @ P_h4p_new, np.zeros((4,4)), atol=1e-10)}")

proj_new = E8 @ P_h4_new
norms_new = np.sqrt(np.sum(proj_new**2, axis=1))

print(f"\nNew H4 projection norms:")
nc_new = Counter(np.round(norms_new, 5))
for n, c in sorted(nc_new.items()):
    print(f"  |pi(r)| = {n:.5f}: {c}")

if len(nc_new) == 2:
    nv = sorted(nc_new.keys())
    print(f"\n*** SUCCESS: 120 + 120 split! ***")
    print(f"Ratio: {nv[1]/nv[0]:.6f}, phi = {phi:.6f}")
elif len(nc_new) > 2:
    # Check if there are really just 2 clusters
    all_norms = sorted(set(np.round(norms_new, 5)))
    print(f"\nAll distinct norms: {len(all_norms)}")

    # Check ratios between consecutive norms
    for i in range(len(all_norms)-1):
        print(f"  {all_norms[i]:.6f} -> {all_norms[i+1]:.6f}: ratio = {all_norms[i+1]/all_norms[i]:.4f}")

    # Check: are norms related to phi?
    print(f"\nNorm^2 values vs phi-related values:")
    for n in all_norms:
        n2 = n**2
        # Check n^2 = a + b*phi for small integers a, b
        best = None
        for a in range(-3, 4):
            for b in range(-3, 4):
                if abs(n2 - (a + b*phi)) < 0.01:
                    best = (a, b)
        print(f"  n = {n:.6f}, n^2 = {n2:.6f}, = {best} (a + b*phi)" if best else f"  n = {n:.6f}, n^2 = {n2:.6f}, no Z[phi] match")

# ============================================================
# THE KEY REALIZATION
# ============================================================
print(f"\n{'='*70}")
print("KEY ANALYSIS: Norm shell structure")
print(f"{'='*70}")

# Let's look at the norm^2 values more carefully
proj3 = E8 @ P_h4_new
n2_vals = np.sum(proj3**2, axis=1)
n2_unique = sorted(set(np.round(n2_vals, 8)))

print(f"Distinct norm^2 values in H4 projection:")
counts = []
for n2 in n2_unique:
    c = np.sum(np.abs(n2_vals - n2) < 0.0001)
    counts.append(c)
    # Check if n2 = 2*cos^2(pi*m/30) for some m
    # Or n2 = 2*sin^2(pi*m/30)
    # The projected norm^2 of a root r onto V_H4 should be
    # sum over H4 eigenplanes of |projection onto that plane|^2

    # Try to express as a + b*sqrt(5)
    # n2 = c0 + c1*phi = c0 + c1*(1+sqrt(5))/2
    # 2*n2 = 2c0 + c1 + c1*sqrt(5)
    # So rational part = 2c0 + c1, irrational part = c1
    # c1 = (2*n2 - round(2*n2)) * 2 / sqrt(5) ... hmm

    print(f"  n^2 = {n2:.8f}: {c} roots  (n = {np.sqrt(n2):.6f})")

print(f"\nTotal: {sum(counts)} roots, counts: {counts}")

# The 600-cell in R^4 has vertices at distance 2*sin(pi/10) from center
# when inscribed in unit sphere. But E8 roots project differently.
# The number of shells depends on the geometry.

# Check: in the H4' projection, do we get the same shells?
proj3p = E8 @ P_h4p_new
n2p_vals = np.sum(proj3p**2, axis=1)
n2p_unique = sorted(set(np.round(n2p_vals, 8)))

print(f"\nH4' projection shells:")
for n2 in n2p_unique:
    c = np.sum(np.abs(n2p_vals - n2) < 0.0001)
    print(f"  n^2 = {n2:.8f}: {c} roots")

# THE SHELLS should match between H4 and H4': if r has |pi(r)|^2 = a,
# then |pi'(r)|^2 = 2 - a. So the shells are symmetric about 1.

# Check: for each root, is |pi|^2 + |pi'|^2 = 2?
check_sum = np.sum(proj3**2, axis=1) + np.sum(proj3p**2, axis=1)
print(f"\n|pi|^2 + |pi'|^2 = {sorted(set(np.round(check_sum, 8)))}")

# So if |pi|^2 = a, then |pi'|^2 = 2-a
# The shells in H4 and H4' are related by a -> 2-a.

# ============================================================
# RESOLUTION: The 120+120 split is by SIGN, not by NORM
# ============================================================
print(f"\n{'='*70}")
print("RESOLUTION: How are the 240 roots split into two 600-cells?")
print(f"{'='*70}")

print("""
The projection onto H4 gives MULTIPLE norm shells, not exactly 2.
The split into two 600-cells is NOT by projected norm.

The correct split: each E8 root r has projections pi(r) in V_H4
and pi'(r) in V_H4'. The 600-cell structure appears when we look at
pi(r) (4D coordinates) and recognize that the 240 projected points
form two CONCENTRIC 600-cells at different radii.

But wait -- the counts show NOT 120+120 at two radii but
rather (24, 56, 40, 40, 56, 24) = 240 at 6 radii.

This suggests the projection is NOT the "correct" one, or
the split into 600-cells works differently than I assumed.

Let me check: do the 240 projected 4D vectors, when SCALED to
the unit sphere, form 120 or 240 distinct directions?
""")

# Project and normalize to unit sphere
proj_unit = proj3 / np.linalg.norm(proj3, axis=1, keepdims=True)

# Count distinct directions (up to sign)
directions = set()
for v in proj_unit:
    # Canonical form: first nonzero component positive
    for c in v:
        if abs(c) > 1e-8:
            if c < 0:
                v = -v
            break
    directions.add(tuple(np.round(v, 6)))

print(f"Distinct directions on unit sphere: {len(directions)}")
print(f"(240 roots, but some may project to same direction)")

# Check: how many of 240 project to 120 directions (with multiplicity 2)?
dir_counter = Counter()
for v in proj_unit:
    for c in v:
        if abs(c) > 1e-8:
            if c < 0:
                v = -v
            break
    key = tuple(np.round(v, 6))
    dir_counter[key] += 1

mult_dist = Counter(dir_counter.values())
print(f"Multiplicity distribution: {dict(mult_dist)}")

if 2 in mult_dist and mult_dist[2] * 2 == 240:
    print("Each direction has multiplicity 2!")
    print("The 240 roots project to 120 directions, each with 2 roots.")
    print("This means the projection is 2-to-1, giving ONE 600-cell!")
elif 1 in mult_dist and mult_dist[1] == 240:
    print("All 240 directions are distinct -> 240-vertex polytope, not 600-cell.")
elif 1 in mult_dist and mult_dist[1] == 120:
    print("120 distinct + 120 paired = mixed structure")

# ============================================================
# Check: do the directions form a 600-cell?
# ============================================================
if len(directions) == 120:
    print(f"\n120 directions found! Checking 600-cell structure...")
    dir_array = np.array(sorted(directions))
    IP_dir = dir_array @ dir_array.T
    ip_c = Counter(np.round(IP_dir[np.triu_indices(120, k=1)], 3))
    print(f"Inner products of 120 unit directions:")
    for v, c in sorted(ip_c.items()):
        print(f"  {v:+.3f}: {c}")

    # Compare with 600-cell signature
    is_600 = True
    expected = {-1.0: 60, -0.809: 720, -0.5: 1200, -0.309: 720,
                0.0: 1800, 0.309: 720, 0.5: 1200, 0.809: 720}
    for ev, ec in expected.items():
        found = any(abs(v - ev) < 0.01 and c == ec for v, c in ip_c.items())
        if not found:
            is_600 = False
            break

    print(f"Is 600-cell? {is_600}")

# ============================================================
# FINAL ANSWER
# ============================================================
print(f"\n{'='*70}")
print("DEFINITIVE ANSWER TO TASK 1")
print(f"{'='*70}")

print(f"""
The Coxeter projection R^8 -> R^4 (onto H4 eigenspace) maps
240 E8 roots to:

If {len(directions)} distinct directions on unit sphere:
  The projection is {'2-to-1' if len(directions) == 120 else '1-to-1'}.

The number of norm shells = {len(n2_unique)}.
Shell sizes: {counts}

This structure reflects the decomposition of the E8 root system
under the H4 subgroup action.
""")

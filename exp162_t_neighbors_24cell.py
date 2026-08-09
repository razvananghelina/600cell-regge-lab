"""
EXP-162: Do the 24 T-neighbors of each S-vertex form a 24-cell?
================================================================
From exp160: each vertex in one 600-cell copy has exactly 24 neighbors
in the other copy (T-neighbors). 24 is the number of vertices of a 24-cell.

Uses ICOSIAN embedding: v = a + b*phi decomposition.
S in R^8: (a+b, b). T in R^8: (a-b, -a). Gram = [[1,1],[1,2]].
"""

import numpy as np
from itertools import product, permutations
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-162: T-Neighbors and the 24-Cell (Icosian Method)")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell with (a, b) decomposition
# ============================================================
print("\n--- Step 1: Build 600-cell with Q(phi) decomposition ---")

# Each vertex v = (v0,v1,v2,v3) where v_k = a_k + b_k * phi
# Store as list of (a_vec, b_vec) pairs

vertices_ab = []  # list of (a_4vec, b_4vec)
vertices_r4 = []  # R^4 coordinates

def add_vertex(a_vec, b_vec):
    r4 = tuple(np.round(np.array(a_vec) + np.array(b_vec) * PHI, 10))
    a_t = tuple(a_vec)
    b_t = tuple(b_vec)
    vertices_ab.append((a_t, b_t))
    vertices_r4.append(r4)

# Type 1: (+-1, 0, 0, 0) permutations - 8 vertices
# a = (+-1,0,0,0), b = (0,0,0,0)
for i in range(4):
    for s in [1, -1]:
        a = [0,0,0,0]; a[i] = s
        add_vertex(a, [0,0,0,0])

# Type 2: (+-1/2, +-1/2, +-1/2, +-1/2) - 16 vertices
# a = (+-1/2,...), b = (0,...)
for signs in product([0.5, -0.5], repeat=4):
    add_vertex(list(signs), [0,0,0,0])

# Type 3: even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0) - 96 vertices
# phi/2 = 0 + (1/2)*phi  => a=0, b=1/2
# 1/2 = 1/2 + 0*phi      => a=1/2, b=0
# 1/(2phi) = -1/2 + (1/2)*phi  => a=-1/2, b=1/2
# 0 = 0 + 0*phi           => a=0, b=0

# (a, b) values for the base: (phi/2, 1/2, 1/(2phi), 0)
ab_base = [
    (0, 0.5),      # phi/2
    (0.5, 0),      # 1/2
    (-0.5, 0.5),   # 1/(2phi)
    (0, 0),         # 0
]

even_perms = []
for p in permutations(range(4)):
    inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
    if inv % 2 == 0:
        even_perms.append(p)

for perm in even_perms:
    # Permuted (a, b) values
    ab_perm = [ab_base[perm[i]] for i in range(4)]
    a_base = [ab_perm[i][0] for i in range(4)]
    b_base = [ab_perm[i][1] for i in range(4)]

    # Find non-zero positions (where a or b is nonzero)
    nz = [i for i in range(4) if a_base[i] != 0 or b_base[i] != 0]

    for signs in product([1, -1], repeat=len(nz)):
        a_vec = list(a_base)
        b_vec = list(b_base)
        for idx, s in zip(nz, signs):
            a_vec[idx] *= s
            b_vec[idx] *= s
        add_vertex(a_vec, b_vec)

# Remove duplicates
unique = {}
for i, (ab, r4) in enumerate(zip(vertices_ab, vertices_r4)):
    if r4 not in unique:
        unique[r4] = ab

vertices_r4_unique = list(unique.keys())
vertices_ab_unique = [unique[r4] for r4 in vertices_r4_unique]

N = len(vertices_r4_unique)
print(f"  |S| = {N}")
assert N == 120, f"Expected 120, got {N}"

# ============================================================
# Step 2: Embed S and T in R^8
# ============================================================
print("\n--- Step 2: Embed S and T in R^8 ---")

# S: (a+b, b) in R^8
# T: (a-b, -a) in R^8

S_r8 = np.zeros((N, 8))
T_r8 = np.zeros((N, 8))

for i, (a_vec, b_vec) in enumerate(vertices_ab_unique):
    a = np.array(a_vec)
    b = np.array(b_vec)
    S_r8[i, :4] = a + b
    S_r8[i, 4:] = b
    T_r8[i, :4] = a - b
    T_r8[i, 4:] = -a

# Check norms
S_norms_sq = np.sum(S_r8**2, axis=1)
T_norms_sq = np.sum(T_r8**2, axis=1)
print(f"  S norms^2: unique = {np.unique(np.round(S_norms_sq, 6))}")
print(f"  T norms^2: unique = {np.unique(np.round(T_norms_sq, 6))}")

# E8 roots have norm^2 = 2
expected_norm_sq = 2.0
S_ok = np.allclose(S_norms_sq, expected_norm_sq, atol=0.01)
T_ok = np.allclose(T_norms_sq, expected_norm_sq, atol=0.01)
print(f"  S all norm^2=2: {S_ok}")
print(f"  T all norm^2=2: {T_ok}")

if not S_ok:
    # May need to scale
    s_scale = np.sqrt(2 / np.mean(S_norms_sq))
    print(f"  S needs scaling by {s_scale:.4f}")
    S_r8 *= s_scale
    T_r8 *= np.sqrt(2 / np.mean(T_norms_sq))
    S_norms_sq = np.sum(S_r8**2, axis=1)
    T_norms_sq = np.sum(T_r8**2, axis=1)
    print(f"  After scaling: S norms^2 = {np.unique(np.round(S_norms_sq, 4))}")
    print(f"  After scaling: T norms^2 = {np.unique(np.round(T_norms_sq, 4))}")

# Check S-T overlap in R^8
overlap_count = 0
for i in range(N):
    for j in range(N):
        if np.allclose(S_r8[i], T_r8[j], atol=0.001):
            overlap_count += 1
print(f"  S cap T in R^8: {overlap_count} (should be 0)")

# ============================================================
# Step 3: Verify E8 structure
# ============================================================
print("\n--- Step 3: Verify E8 adjacency ---")

# Combine S and T
E8_all = np.vstack([S_r8, T_r8])  # 240 x 8

# Inner products
dots_all = E8_all @ E8_all.T

# S-S adjacency
SS_dots = S_r8 @ S_r8.T
SS_adj = (np.abs(SS_dots - 1.0) < 0.01).astype(int)
np.fill_diagonal(SS_adj, 0)
SS_degrees = SS_adj.sum(axis=1)
print(f"  S-S degree: {Counter(SS_degrees.tolist())}")

# T-T adjacency
TT_dots = T_r8 @ T_r8.T
TT_adj = (np.abs(TT_dots - 1.0) < 0.01).astype(int)
np.fill_diagonal(TT_adj, 0)
TT_degrees = TT_adj.sum(axis=1)
print(f"  T-T degree: {Counter(TT_degrees.tolist())}")

# S-T adjacency (cross-connections)
ST_dots = S_r8 @ T_r8.T
ST_adj = (np.abs(ST_dots - 1.0) < 0.01).astype(int)
ST_degrees = ST_adj.sum(axis=1)  # For each S vertex, how many T neighbors
TS_degrees = ST_adj.T.sum(axis=1)  # For each T vertex, how many S neighbors
print(f"  S-to-T degree: {Counter(ST_degrees.tolist())}")
print(f"  T-to-S degree: {Counter(TS_degrees.tolist())}")

total_degree = SS_degrees + ST_degrees
print(f"  Total S degree: {Counter(total_degree.tolist())}")

# ============================================================
# Step 4: Analyze T-neighbors for each S-vertex
# ============================================================
print("\n--- Step 4: T-neighbors analysis ---")

# Check if norms are correct (need norm^2 = 2 for E8)
if not np.allclose(S_norms_sq, 2.0, atol=0.01):
    print("  WARNING: norms not exactly 2, adjusting threshold...")
    # Use relative threshold
    ss_dot_vals = np.unique(np.round(SS_dots[np.triu_indices(N, k=1)], 3))
    st_dot_vals = np.unique(np.round(ST_dots.flatten(), 3))
    print(f"  S-S dot products: {sorted(ss_dot_vals[ss_dot_vals > 0.5])[-5:]}")
    print(f"  S-T dot products: {sorted(st_dot_vals[st_dot_vals > 0.5])[-5:]}")

# Get dot product spectrum
ss_off = SS_dots[np.triu_indices(N, k=1)]
st_off = ST_dots.flatten()
print(f"\n  S-S dot product spectrum: {Counter(np.round(ss_off, 4).tolist()).most_common(10)}")
print(f"  S-T dot product spectrum: {Counter(np.round(st_off, 4).tolist()).most_common(10)}")

# If standard E8 adjacency (dot=1) doesn't work, find the right threshold
# E8 neighbors are at dot = 1 (for norm sqrt(2) vectors)
# If our norms are different, scale accordingly

actual_norm_sq = np.mean(S_norms_sq)
neighbor_dot = actual_norm_sq / 2  # For E8, neighbor dot = norm^2/2 = 1
print(f"\n  Expected neighbor dot product: {neighbor_dot:.4f}")

# Rebuild adjacency with correct threshold
SS_adj2 = (np.abs(SS_dots - neighbor_dot) < 0.01).astype(int)
np.fill_diagonal(SS_adj2, 0)
ST_adj2 = (np.abs(ST_dots - neighbor_dot) < 0.01).astype(int)

SS_deg2 = SS_adj2.sum(axis=1)
ST_deg2 = ST_adj2.sum(axis=1)
print(f"  S-S degree (adj threshold {neighbor_dot:.3f}): {Counter(SS_deg2.tolist())}")
print(f"  S-T degree (adj threshold {neighbor_dot:.3f}): {Counter(ST_deg2.tolist())}")

# Use 600-cell adjacency (dot = phi/2 in R^4) as baseline
R4 = np.array(vertices_r4_unique)
R4_dots = R4 @ R4.T
R4_adj = (np.abs(R4_dots - PHI/2) < 0.01).astype(int)
np.fill_diagonal(R4_adj, 0)
R4_degrees = R4_adj.sum(axis=1)
print(f"  600-cell R^4 degree: {Counter(R4_degrees.tolist())}")

# ============================================================
# Step 5: Direct 24-cell check using R^4 structure
# ============================================================
print("\n--- Step 5: Direct 24-cell analysis ---")

# The 24-cell vertices in R^4: permutations of (+-1, +-1, 0, 0)
# This is the D4 root system (24 vertices)

# Standard 24-cell
cell24_verts = []
for i in range(4):
    for j in range(i+1, 4):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = [0,0,0,0]; v[i]=si; v[j]=sj
                cell24_verts.append(tuple(v))
cell24 = np.array(cell24_verts)
print(f"  Standard 24-cell: {len(cell24)} vertices")

# 24-cell adjacency: dot product = 1 (for norm sqrt(2) vertices)
c24_dots = cell24 @ cell24.T
c24_adj = (np.abs(c24_dots - 1.0) < 0.01).astype(int)
np.fill_diagonal(c24_adj, 0)
c24_degrees = c24_adj.sum(axis=1)
print(f"  24-cell degree: {np.unique(c24_degrees)} (should be 8)")
c24_edges = c24_adj.sum() // 2
print(f"  24-cell edges: {c24_edges} (should be 96)")

# 5 inscribed 24-cells in 600-cell (from exp154: cosets of 2T in 2I)
# Identify them by clustering vertices with specific dot products

# First: which 600-cell vertices are D4 roots (24-cell vertices)?
r4_norms_sq = np.sum(R4**2, axis=1)
# D4 roots have norm sqrt(2) in standard form, but 600-cell vertices have norm 1
# In 600-cell: the 24 vertices that form a 24-cell are the (+-1,0,0,0) permutations
# AND (+-1/2,...)  -- wait, no.

# The 5 inscribed 24-cells: each is an orbit of 2T (binary tetrahedral, order 24) in 2I
# One coset is the 24 vertices forming the standard 24-cell (vertices of {3,4,3})
# Coset 0: the "Type A + B" vertices (from exp154)

# Let's identify Type A (8 vertices: (+-1,0,0,0) perms)
# and Type B (16 vertices: (+-1/2)^4)
type_A_indices = []
type_B_indices = []
type_C_indices = []

for i, r4 in enumerate(vertices_r4_unique):
    coords = np.array(r4)
    n_nonzero = np.sum(np.abs(coords) > 0.01)
    if n_nonzero == 1 and np.isclose(np.max(np.abs(coords)), 1.0):
        type_A_indices.append(i)
    elif n_nonzero == 4 and np.allclose(np.abs(coords), 0.5, atol=0.01):
        type_B_indices.append(i)
    else:
        type_C_indices.append(i)

print(f"\n  Type A: {len(type_A_indices)}, Type B: {len(type_B_indices)}, Type C: {len(type_C_indices)}")

# Coset 0 = Type A + Type B = first inscribed 24-cell
coset0 = type_A_indices + type_B_indices
print(f"  Coset 0 (A+B): {len(coset0)} vertices")

# Check: does coset 0 form a 24-cell?
c0_coords = R4[coset0]
c0_dots = c0_coords @ c0_coords.T
# For 600-cell vertices (norm 1), 24-cell neighbors are at dot = 0.5
# Because original 24-cell has norm sqrt(2), neighbors at dot=1
# Scaled to norm 1: dot = 1/2
c0_adj = (np.abs(c0_dots - 0.5) < 0.01).astype(int)
np.fill_diagonal(c0_adj, 0)
c0_deg = c0_adj.sum(axis=1)
c0_edges = c0_adj.sum() // 2
print(f"  Coset 0 degree: {np.unique(c0_deg)} (24-cell has 8)")
print(f"  Coset 0 edges: {c0_edges} (24-cell has 96)")

# ============================================================
# Step 6: Find all 5 inscribed 24-cells
# ============================================================
print("\n--- Step 6: Find 5 inscribed 24-cells ---")

# Method: use graph coloring. Two vertices in same 24-cell iff
# they are at distance 2 in 600-cell or have specific dot product.

# 600-cell dot products: 1 (same), phi/2 (d=1), 1/2 (d=2 type), 0 (d=2 type),
# -1/(2phi) (d=3 type), -phi/2 (d=4), -1 (antipodal)

# Vertices in same 24-cell have R^4 dot products from the set:
# For norm-1 24-cell: {1, 0.5, 0, -0.5, -1}
# That means: same 24-cell vertices have dots in {1, 0.5, 0, -0.5, -1}

# Partition using: v and w in same coset iff dot(v,w) in {1, 0.5, 0, -0.5, -1}
# (i.e., dot(v,w) is rational, no phi component)

cosets = [set() for _ in range(5)]
assigned = [False] * N
cosets[0] = set(coset0)
for i in coset0:
    assigned[i] = True

coset_count = 1
for seed in range(N):
    if assigned[seed]:
        continue
    # Start new coset
    current_coset = {seed}
    assigned[seed] = True
    queue = [seed]
    while queue:
        u = queue.pop(0)
        for v in range(N):
            if assigned[v]:
                continue
            d = R4_dots[u, v]
            # Same coset if dot is "rational" (no phi component)
            # Rational dots: {1, 0.5, 0, -0.5, -1}
            if any(abs(d - val) < 0.01 for val in [1.0, 0.5, 0.0, -0.5, -1.0]):
                current_coset.add(v)
                assigned[v] = True
                queue.append(v)

    if len(current_coset) == 24:
        cosets[coset_count] = current_coset
        coset_count += 1
    else:
        print(f"  WARNING: coset from seed {seed} has {len(current_coset)} vertices")

print(f"  Found {coset_count} cosets of size 24")

# Verify each coset forms a 24-cell
for c_idx in range(coset_count):
    c_list = sorted(cosets[c_idx])
    c_coords = R4[c_list]
    c_dots = c_coords @ c_coords.T
    c_adj = (np.abs(c_dots - 0.5) < 0.01).astype(int)
    np.fill_diagonal(c_adj, 0)
    c_edges = c_adj.sum() // 2
    c_deg = np.unique(c_adj.sum(axis=1))
    print(f"  Coset {c_idx}: {len(c_list)} verts, {c_edges} edges, degree {c_deg}")

# ============================================================
# Step 7: Cross-coset adjacency in 600-cell
# ============================================================
print("\n--- Step 7: Cross-coset adjacency (600-cell) ---")

# In 600-cell, are there edges between vertices of same coset?
for c_idx in range(coset_count):
    c_list = sorted(cosets[c_idx])
    intra_edges = 0
    for i in c_list:
        for j in c_list:
            if i < j and R4_adj[i, j]:
                intra_edges += 1
    print(f"  Coset {c_idx}: {intra_edges} intra-coset 600-cell edges")

# Cross-coset 600-cell edges
for ci in range(coset_count):
    for cj in range(ci+1, coset_count):
        cross = 0
        for i in cosets[ci]:
            for j in cosets[cj]:
                if R4_adj[i, j]:
                    cross += 1
        print(f"  Coset {ci}-{cj}: {cross} cross 600-cell edges")

# ============================================================
# Step 8: S-to-T analysis via (a,b) decomposition
# ============================================================
print("\n--- Step 8: Cross-copy analysis via (a,b) decomposition ---")

# For each S-vertex with (a,b), find T-vertices (a'-b', -a') that are E8-adjacent
# E8 adjacency in R^8: dot product = norm^2/2

# First verify S norms
S_norms_check = np.sum(S_r8**2, axis=1)
unique_norms = np.unique(np.round(S_norms_check, 6))
print(f"  S R^8 norms^2: {unique_norms}")

# If all norms are the same but not 2, we can still compute
# E8 adjacency: dot = 1 (for norm sqrt(2))
# For general norm r: adjacency at dot = r^2/2

if len(unique_norms) == 1:
    norm_sq = unique_norms[0]
    adj_dot = norm_sq / 2
    print(f"  Adjacency dot product: {adj_dot:.4f}")

    ST_adj_v2 = (np.abs(ST_dots - adj_dot) < 0.01 * norm_sq).astype(int)
    ST_deg_v2 = ST_adj_v2.sum(axis=1)
    print(f"  S-to-T degree: {Counter(ST_deg_v2.tolist())}")

    # For each S vertex, get its T neighbors
    for s_idx in range(min(5, N)):
        t_neighbors = np.where(ST_adj_v2[s_idx] > 0)[0]
        if len(t_neighbors) > 0:
            # Get R^4 coordinates of these T vertices
            t_r4 = R4[t_neighbors]

            # Mutual R^4 dot products
            t_mutual_dots = t_r4 @ t_r4.T
            off_diag = []
            for i in range(len(t_neighbors)):
                for j in range(i+1, len(t_neighbors)):
                    off_diag.append(round(t_mutual_dots[i,j], 4))

            # Check if they form a 24-cell
            # 24-cell (norm 1): mutual dots in {0.5, 0, -0.5, -1}
            # Adjacency at dot = 0.5, edges = 96

            t_adj = (np.abs(t_mutual_dots - 0.5) < 0.01).astype(int)
            np.fill_diagonal(t_adj, 0)
            t_edges = t_adj.sum() // 2
            t_deg = t_adj.sum(axis=1)

            # Which cosets are they in?
            coset_membership = []
            for tn in t_neighbors:
                for ci in range(coset_count):
                    if tn in cosets[ci]:
                        coset_membership.append(ci)
                        break
                else:
                    coset_membership.append(-1)

            if s_idx < 3:
                print(f"\n  S-vertex {s_idx}: {len(t_neighbors)} T-neighbors")
                print(f"    R^4 dot spectrum: {Counter(off_diag).most_common(8)}")
                print(f"    Mutual edges (24-cell test): {t_edges} (need 96)")
                print(f"    Degrees: {Counter(t_deg.tolist())}")
                print(f"    Coset membership: {Counter(coset_membership)}")

                # Check: do they form a single coset?
                if len(set(coset_membership)) == 1:
                    print(f"    >> All in SAME coset ({coset_membership[0]})")
                else:
                    print(f"    >> Spread across cosets: {set(coset_membership)}")

else:
    print("  Multiple norm values - checking structure differently")
    # Just use whatever the dot product spectrum gives us
    st_flat = ST_dots.flatten()
    st_unique = np.unique(np.round(st_flat, 4))
    print(f"  S-T dot product unique values: {st_unique}")

    # Find the value closest to what would be "adjacency"
    # In E8, adjacency is at dot = 1 (norm sqrt(2)), or dot = norm^2/2
    for candidate in sorted(st_unique, reverse=True):
        if abs(candidate) > 0.01 and abs(candidate) < 2:
            count = np.sum(np.abs(ST_dots - candidate) < 0.005)
            per_vertex = count / N
            print(f"    dot={candidate:.4f}: count={count}, per S-vertex={per_vertex:.1f}")

# ============================================================
# Step 9: Alternative - build cross adjacency from scratch
# ============================================================
print("\n--- Step 9: Direct cross-adjacency computation ---")

# For each pair (s_i, t_j), compute their R^8 distance squared
# E8 neighbors are at distance^2 = 2 (for norm sqrt(2) roots at dot=1)
# dist^2 = |s-t|^2 = |s|^2 + |t|^2 - 2*dot(s,t) = 2 + 2 - 2*1 = 2

dist_sq_ST = np.zeros((N, N))
for i in range(N):
    diff = S_r8[i] - T_r8  # (N, 8)
    dist_sq_ST[i] = np.sum(diff**2, axis=1)

# Unique distance values
dist_unique = np.unique(np.round(dist_sq_ST.flatten(), 4))
print(f"  S-T distance^2 unique values ({len(dist_unique)}):")
for d2 in sorted(dist_unique)[:10]:
    count = np.sum(np.abs(dist_sq_ST - d2) < 0.005)
    per_vertex = count / N
    print(f"    d^2={d2:.4f}: count={count}, per vertex={per_vertex:.1f}")

# E8 neighbor distance: for norm^2=r^2 roots, neighbors at d^2 = 2*r^2 - 2*r^2/2 = r^2
# Actually: d^2 = |s-t|^2 = |s|^2 + |t|^2 - 2*dot = 2*r^2 - 2*adj_dot
# For E8 (r^2=2, adj_dot=1): d^2 = 4-2 = 2
# For our vectors: d^2 = 2*norm_sq - 2*adj_dot

# The SMALLEST nonzero d^2 corresponds to E8 nearest neighbors
min_d2 = sorted(dist_unique)[0] if dist_unique[0] > 0.01 else sorted(dist_unique)[1]
print(f"\n  Smallest S-T distance^2: {min_d2:.4f}")
cross_adj = (np.abs(dist_sq_ST - min_d2) < 0.005).astype(int)
cross_degrees = cross_adj.sum(axis=1)
print(f"  Cross-degree at d^2={min_d2:.3f}: {Counter(cross_degrees.tolist())}")

# If cross-degree is 24, check 24-cell structure
if 24 in cross_degrees.tolist():
    print("\n  >> Found 24 cross-neighbors! Checking 24-cell structure...")

    for s_idx in range(min(5, N)):
        t_nbrs = np.where(cross_adj[s_idx] > 0)[0]
        if len(t_nbrs) != 24:
            continue

        t_r4 = R4[t_nbrs]
        t_dots = t_r4 @ t_r4.T
        t_adj24 = (np.abs(t_dots - 0.5) < 0.01).astype(int)
        np.fill_diagonal(t_adj24, 0)
        t_edges = t_adj24.sum() // 2
        t_deg = np.unique(t_adj24.sum(axis=1))

        # Coset membership
        cmem = []
        for tn in t_nbrs:
            for ci in range(coset_count):
                if tn in cosets[ci]:
                    cmem.append(ci)
                    break

        # Which coset is the S vertex in?
        s_coset = -1
        for ci in range(coset_count):
            if s_idx in cosets[ci]:
                s_coset = ci
                break

        print(f"\n  S-vertex {s_idx} (coset {s_coset}): {len(t_nbrs)} T-neighbors")
        print(f"    Edges: {t_edges} (24-cell=96)")
        print(f"    Degrees: {t_deg} (24-cell=8)")
        print(f"    T-neighbor cosets: {Counter(cmem)}")

        if t_edges == 96:
            print(f"    >> CONFIRMED: 24-cell structure!")
        elif t_edges == 0:
            # Try different adjacency threshold
            off_diag = t_dots[np.triu_indices(24, k=1)]
            print(f"    Dot products: {Counter(np.round(off_diag, 3).tolist()).most_common(6)}")

            # Try all possible adjacency thresholds
            for threshold in sorted(set(np.round(off_diag, 3))):
                adj_test = (np.abs(t_dots - threshold) < 0.01).astype(int)
                np.fill_diagonal(adj_test, 0)
                edges = adj_test.sum() // 2
                if edges > 0:
                    deg = np.unique(adj_test.sum(axis=1))
                    print(f"    At dot={threshold:.3f}: {edges} edges, deg={deg}")

# ============================================================
# Step 10: Also check S-S neighborhood structure
# ============================================================
print("\n--- Step 10: S-S neighborhood in R^8 ---")

dist_sq_SS = np.zeros((N, N))
for i in range(N):
    diff = S_r8[i] - S_r8
    dist_sq_SS[i] = np.sum(diff**2, axis=1)

ss_unique = np.unique(np.round(dist_sq_SS[np.triu_indices(N, k=1)], 4))
print(f"  S-S distance^2 values ({len(ss_unique)}):")
for d2 in sorted(ss_unique)[:8]:
    count = np.sum(np.abs(dist_sq_SS - d2) < 0.005) - (N if d2 < 0.01 else 0)
    per_vertex = count / N
    print(f"    d^2={d2:.4f}: count={count}, per vertex={per_vertex:.1f}")

# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Questions:
1. Do the 24 T-neighbors form a 24-cell? -> check edges=96, degree=8
2. Do T-neighbors belong to a single coset (inscribed 24-cell)?
3. Is the cross-copy structure the D4 root system?

STATUS: See output above
""")

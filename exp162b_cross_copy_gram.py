"""
EXP-162b: Cross-copy structure via GRAM inner product
=====================================================
Fix of exp162: use Gram matrix directly instead of R^8 embedding.

For v = a + b*phi in 600-cell (|v|=1), the E8 Gram inner product is:
  G(v1, v2) = a1.a2 + a1.b2 + b1.a2 + 2*b1.b2   [for same copy S-S]

For T = phi' * S: coefficients (a-b, -a), so:
  G(S_i, T_j) = -(a_i.b_j + b_i.a_j + b_i.b_j)  [cross S-T]

Scale by 2 for E8 (root norm^2=2). E8 adjacency at Gram=1.
"""

import numpy as np
from itertools import product, permutations
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-162b: Cross-Copy Structure (Gram Method)")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell with (a, b) decomposition
# ============================================================
print("\n--- Step 1: 600-cell with (a,b) decomposition ---")

ab_list = []  # (a_4vec, b_4vec)
r4_list = []

def add_v(a, b):
    r4 = tuple(np.round(np.array(a) + np.array(b) * PHI, 10))
    ab_list.append((tuple(a), tuple(b)))
    r4_list.append(r4)

# Type A: (+-1, 0, 0, 0) - a=(+-1,0,0,0), b=0
for i in range(4):
    for s in [1, -1]:
        a = [0,0,0,0]; a[i] = s
        add_v(a, [0,0,0,0])

# Type B: (+-1/2)^4 - a=(+-1/2,...), b=0
for signs in product([0.5, -0.5], repeat=4):
    add_v(list(signs), [0,0,0,0])

# Type C: even perms of (phi/2, 1/2, 1/(2phi), 0)
# phi/2: a=0, b=1/2.  1/2: a=1/2, b=0.  1/(2phi): a=-1/2, b=1/2.  0: a=0, b=0.
ab_base = [(0, 0.5), (0.5, 0), (-0.5, 0.5), (0, 0)]

even_perms = []
for p in permutations(range(4)):
    inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
    if inv % 2 == 0:
        even_perms.append(p)

for perm in even_perms:
    ab_perm = [ab_base[perm[i]] for i in range(4)]
    a_base = [ab_perm[i][0] for i in range(4)]
    b_base = [ab_perm[i][1] for i in range(4)]
    nz = [i for i in range(4) if a_base[i] != 0 or b_base[i] != 0]
    for signs in product([1, -1], repeat=len(nz)):
        a = list(a_base); b = list(b_base)
        for idx, s in zip(nz, signs):
            a[idx] *= s; b[idx] *= s
        add_v(a, b)

# Deduplicate
unique = {}
for (ab, r4) in zip(ab_list, r4_list):
    if r4 not in unique:
        unique[r4] = ab

r4_keys = list(unique.keys())
ab_vals = [unique[k] for k in r4_keys]
N = len(r4_keys)
print(f"  |S| = {N}")
assert N == 120

R4 = np.array(r4_keys)
A = np.array([ab[0] for ab in ab_vals])  # (120, 4)
B = np.array([ab[1] for ab in ab_vals])  # (120, 4)

# Verify: R4 = A + B*phi
check = np.max(np.abs(R4 - (A + B * PHI)))
print(f"  R4 = A + B*phi check: max err = {check:.2e}")

# ============================================================
# Step 2: Gram inner products (scaled by 2 for E8)
# ============================================================
print("\n--- Step 2: Gram inner products ---")

# S-S Gram: 2*(a1.a2 + a1.b2 + b1.a2 + 2*b1.b2)
SS_gram = 2 * (A @ A.T + A @ B.T + B @ A.T + 2 * B @ B.T)

# S-S self (should be 2 for all)
ss_self = np.diag(SS_gram)
print(f"  S-S self Gram: {np.unique(np.round(ss_self, 4))}")

# S-S off-diagonal spectrum
ss_off = SS_gram[np.triu_indices(N, k=1)]
print(f"  S-S Gram spectrum: {Counter(np.round(ss_off, 4).tolist()).most_common(10)}")

# S-S adjacency (Gram = 1 for E8 neighbors)
SS_adj = (np.abs(SS_gram - 1.0) < 0.01).astype(int)
np.fill_diagonal(SS_adj, 0)
SS_deg = SS_adj.sum(axis=1)
print(f"  S-S E8 degree: {Counter(SS_deg.tolist())}")

# T coefficients: for phi'*v = (a-b) + (-a)*phi
A_T = A - B   # a_T = a - b
B_T = -A      # b_T = -a

# T-T Gram: 2*(aT1.aT2 + aT1.bT2 + bT1.aT2 + 2*bT1.bT2)
TT_gram = 2 * (A_T @ A_T.T + A_T @ B_T.T + B_T @ A_T.T + 2 * B_T @ B_T.T)
tt_self = np.diag(TT_gram)
print(f"  T-T self Gram: {np.unique(np.round(tt_self, 4))}")

TT_adj = (np.abs(TT_gram - 1.0) < 0.01).astype(int)
np.fill_diagonal(TT_adj, 0)
TT_deg = TT_adj.sum(axis=1)
print(f"  T-T E8 degree: {Counter(TT_deg.tolist())}")

# S-T Gram: 2*(a_S.a_T + a_S.b_T + b_S.a_T + 2*b_S.b_T)
ST_gram = 2 * (A @ A_T.T + A @ B_T.T + B @ A_T.T + 2 * B @ B_T.T)
st_off = ST_gram.flatten()
print(f"  S-T Gram spectrum: {Counter(np.round(st_off, 4).tolist()).most_common(10)}")

ST_adj = (np.abs(ST_gram - 1.0) < 0.01).astype(int)
ST_deg = ST_adj.sum(axis=1)
print(f"  S-to-T E8 degree: {Counter(ST_deg.tolist())}")

# Total E8 degree
total_deg = SS_deg + ST_deg
print(f"\n  Total E8 degree: {Counter(total_deg.tolist())}")
print(f"  Expected: 56 for all")

# ============================================================
# Step 3: Verify E8 structure
# ============================================================
print("\n--- Step 3: Verify E8 ---")

# Build full 240x240 Gram matrix
# Upper-left: SS, Upper-right: ST, Lower-left: ST^T, Lower-right: TT
full_gram = np.block([[SS_gram, ST_gram], [ST_gram.T, TT_gram]])
full_adj = (np.abs(full_gram - 1.0) < 0.01).astype(int)
np.fill_diagonal(full_adj, 0)
full_deg = full_adj.sum(axis=1)
print(f"  Full E8 degree: {Counter(full_deg.tolist())}")

total_edges = full_adj.sum() // 2
print(f"  Total edges: {total_edges} (E8 has 6720)")

# Check kissing number: each root has 56 neighbors at dot=1, self=2, and dot=-2 for antipode
gram_vals = full_gram[np.triu_indices(240, k=1)]
print(f"  Gram value spectrum: {Counter(np.round(gram_vals, 4).tolist()).most_common(10)}")

# ============================================================
# Step 4: T-neighbor analysis
# ============================================================
print("\n--- Step 4: T-neighbors of S-vertices ---")

# Classify vertices
type_indices = {'A': [], 'B': [], 'C': []}
for i in range(N):
    if np.max(np.abs(B[i])) < 0.01:  # b = 0
        if np.sum(np.abs(A[i]) > 0.01) == 1:
            type_indices['A'].append(i)
        else:
            type_indices['B'].append(i)
    else:
        type_indices['C'].append(i)

print(f"  Type A: {len(type_indices['A'])}, B: {len(type_indices['B'])}, C: {len(type_indices['C'])}")

# 600-cell R^4 adjacency
R4_dots = R4 @ R4.T
R4_adj = (np.abs(R4_dots - PHI/2) < 0.01).astype(int)
np.fill_diagonal(R4_adj, 0)

# For each S vertex, find T neighbors and analyze
for vtype in ['A', 'B', 'C']:
    indices = type_indices[vtype]
    if not indices:
        continue
    test_idx = indices[0]
    t_nbrs = np.where(ST_adj[test_idx] > 0)[0]
    print(f"\n  {vtype}-vertex {test_idx}: {len(t_nbrs)} T-neighbors")

    if len(t_nbrs) == 0:
        print(f"    No T-neighbors!")
        continue

    # R^4 coordinates of T-neighbors
    t_r4 = R4[t_nbrs]

    # Mutual R^4 dot products
    t_dots = t_r4 @ t_r4.T
    off = []
    for i in range(len(t_nbrs)):
        for j in range(i+1, len(t_nbrs)):
            off.append(round(t_dots[i,j], 4))
    print(f"    R^4 dot spectrum: {Counter(off).most_common(8)}")

    # Check 24-cell: adjacency at dot=0.5 for norm-1 vertices
    t_adj24 = (np.abs(t_dots - 0.5) < 0.01).astype(int)
    np.fill_diagonal(t_adj24, 0)
    t_edges = t_adj24.sum() // 2
    t_deg = t_adj24.sum(axis=1)
    print(f"    24-cell test: edges={t_edges} (need 96), degree={np.unique(t_deg)}")

    # Types of T-neighbors
    t_types = []
    for tn in t_nbrs:
        if np.max(np.abs(B[tn])) < 0.01:
            t_types.append('AB')
        else:
            t_types.append('C')
    print(f"    T-neighbor types: {Counter(t_types)}")

    # Are T-neighbors in the same inscribed 24-cell?
    # Two vertices in same 24-cell iff their R^4 dot is rational (no phi component)
    # Check: for each pair, is the dot rational?
    n_rational = 0
    n_irrational = 0
    for i in range(len(t_nbrs)):
        for j in range(i+1, len(t_nbrs)):
            d = t_dots[i,j]
            # Rational values: 1, 0.5, 0, -0.5, -1
            is_rat = any(abs(d - v) < 0.01 for v in [1, 0.5, 0, -0.5, -1])
            if is_rat:
                n_rational += 1
            else:
                n_irrational += 1
    total_pairs = len(t_nbrs) * (len(t_nbrs)-1) // 2
    print(f"    Rational dot pairs: {n_rational}/{total_pairs}, Irrational: {n_irrational}")

    # 600-cell adjacency among T-neighbors
    t_600adj = R4_adj[np.ix_(t_nbrs, t_nbrs)]
    t_600edges = t_600adj.sum() // 2
    t_600deg = t_600adj.sum(axis=1)
    print(f"    600-cell edges among T-nbrs: {t_600edges}")
    print(f"    600-cell degrees: {Counter(t_600deg.tolist())}")

# ============================================================
# Step 5: Comprehensive analysis across all vertices
# ============================================================
print("\n--- Step 5: Statistics over all vertices ---")

all_t_counts = ST_deg
all_s_in_t_counts = np.zeros(N, dtype=int)  # How many S-neighbors each T-vertex has
for j in range(N):
    all_s_in_t_counts[j] = ST_adj[:, j].sum()

print(f"  S-to-T degree distribution: {Counter(all_t_counts.tolist())}")
print(f"  T-to-S degree distribution: {Counter(all_s_in_t_counts.tolist())}")

# For ALL vertices with 24 T-neighbors
has_24 = np.where(all_t_counts == 24)[0]
print(f"\n  Vertices with exactly 24 T-neighbors: {len(has_24)}")

n_24cell = 0
n_not_24cell = 0
coset_patterns = []

for s_idx in has_24[:20]:  # Check first 20
    t_nbrs = np.where(ST_adj[s_idx] > 0)[0]
    t_r4 = R4[t_nbrs]
    t_dots = t_r4 @ t_r4.T
    t_adj24 = (np.abs(t_dots - 0.5) < 0.01).astype(int)
    np.fill_diagonal(t_adj24, 0)
    t_edges = t_adj24.sum() // 2

    # T-neighbor types
    n_ab = sum(1 for tn in t_nbrs if np.max(np.abs(B[tn])) < 0.01)

    if t_edges == 96:
        n_24cell += 1
    else:
        n_not_24cell += 1

    coset_patterns.append((t_edges, n_ab, 24 - n_ab))

print(f"  24-cell confirmations: {n_24cell}/{min(20, len(has_24))}")
print(f"  Non-24-cell: {n_not_24cell}/{min(20, len(has_24))}")

print(f"\n  Patterns (edges, n_AB, n_C):")
for p in sorted(set(coset_patterns)):
    count = coset_patterns.count(p)
    print(f"    edges={p[0]}, AB={p[1]}, C={p[2]}: {count} vertices")

# ============================================================
# Step 6: 5 inscribed 24-cells revisited
# ============================================================
print("\n--- Step 6: Find 5 inscribed 24-cells ---")

# Use the 2T coset method: the binary tetrahedral group 2T (order 24)
# is a subgroup of 2I (order 120).
# 2I / 2T has 5 cosets.

# 2T = {(+-1,0,0,0), (+-1/2)^4, (0,0,0,0) excluded}
# Actually 2T consists of the 24 vertices that form a 24-cell
# Coset 0 = 2T itself = Type A + Type B

# For the other cosets, use left multiplication by a Type C vertex
# Pick a Type C vertex as coset representative
c_rep = type_indices['C'][0]
print(f"  Coset representative: vertex {c_rep}, R^4 = {np.round(R4[c_rep], 4)}")

# Quaternion multiplication: for q1=(a,b,c,d) and q2=(e,f,g,h):
# q1*q2 = (ae-bf-cg-dh, af+be+ch-dg, ag-bh+ce+df, ah+bg-cf+de)
def quat_mult(q1, q2):
    a,b,c,d = q1
    e,f,g,h = q2
    return np.array([
        a*e - b*f - c*g - d*h,
        a*f + b*e + c*h - d*g,
        a*g - b*h + c*e + d*f,
        a*h + b*g - c*f + d*e
    ])

# Build coset 0 (2T)
coset0 = set(type_indices['A'] + type_indices['B'])
print(f"  Coset 0: {len(coset0)} vertices")

# Build other cosets by left-multiplying coset 0 by Type C vertices
cosets = [coset0]
used = set(coset0)

for c_idx in type_indices['C']:
    if c_idx in used:
        continue
    rep = R4[c_idx]

    # Multiply rep by each coset 0 vertex
    new_coset = set()
    for v_idx in coset0:
        prod = quat_mult(rep, R4[v_idx])
        # Find this vertex in R4
        dists = np.sum((R4 - prod)**2, axis=1)
        closest = np.argmin(dists)
        if dists[closest] < 0.001:
            new_coset.add(closest)

    if len(new_coset) == 24 and not new_coset & used:
        cosets.append(new_coset)
        used |= new_coset
        print(f"  Coset {len(cosets)-1}: {len(new_coset)} vertices (rep=vertex {c_idx})")

    if len(cosets) == 5:
        break

print(f"\n  Total cosets found: {len(cosets)}")
print(f"  Total vertices covered: {len(used)}/120")

# Verify each coset is a 24-cell
for ci, coset in enumerate(cosets):
    c_list = sorted(coset)
    c_r4 = R4[c_list]
    c_dots = c_r4 @ c_r4.T
    c_adj = (np.abs(c_dots - 0.5) < 0.01).astype(int)
    np.fill_diagonal(c_adj, 0)
    c_edges = c_adj.sum() // 2
    c_deg = np.unique(c_adj.sum(axis=1))

    # Intra-coset 600-cell edges
    intra_600 = R4_adj[np.ix_(c_list, c_list)].sum() // 2
    print(f"  Coset {ci}: 24-cell edges={c_edges}, 600-cell intra={intra_600}, degree={c_deg}")

# Cross-coset analysis
if len(cosets) == 5:
    print(f"\n  Cross-coset 600-cell edges:")
    for ci in range(5):
        for cj in range(ci+1, 5):
            cross = 0
            for i in cosets[ci]:
                for j in cosets[cj]:
                    if R4_adj[i, j]:
                        cross += 1
            print(f"    Coset {ci}-{cj}: {cross}")

# ============================================================
# Step 7: T-neighbors vs cosets
# ============================================================
print("\n--- Step 7: T-neighbors and coset structure ---")

if len(cosets) == 5:
    for s_idx in range(min(10, N)):
        if ST_deg[s_idx] != 24:
            continue
        t_nbrs = np.where(ST_adj[s_idx] > 0)[0]

        # Which coset is S vertex in?
        s_coset = -1
        for ci, coset in enumerate(cosets):
            if s_idx in coset:
                s_coset = ci
                break

        # Which cosets are T-neighbors in?
        t_coset_counts = Counter()
        for tn in t_nbrs:
            for ci, coset in enumerate(cosets):
                if tn in coset:
                    t_coset_counts[ci] += 1
                    break

        print(f"  S-vertex {s_idx} (coset {s_coset}): T-neighbors in cosets {dict(t_coset_counts)}")

# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Key questions:
1. Is the E8 structure correct? (degree 56, 6720 edges)
2. Do T-neighbors form 24-cells?
3. T-neighbors from same or different coset as source?
4. Complete 5-partite structure in 600-cell confirmed?

STATUS: Check output above
""")

"""
EXP-143: Edge-Based Gauge - Fix (1,2,9) -> (1,3,8)
=====================================================
Currently: per C vertex, edges split as 1 AC + 2 BC + 9 CC = 12
SM needs: 1 U(1) + 3 SU(2) + 8 SU(3) = 12

Mismatch: 2 BC vs 3 SU(2), 9 CC vs 8 SU(3)
Question: can we find a natural 1+3+8 split of the 12 edges?

Approaches:
1. Local CC subgraph structure: 9 verts with degrees [3,3,3,3,3,3,4,4,4]
   - The 3 degree-4 vertices are special. Can 1 CC edge become "SU(2)-like"?
2. Use A vs B gauge type to induce splitting
3. Spectral splitting of local neighborhood
4. Geometric analysis: which CC edges are "closer" to gauge vertices?
5. 24-cell (D4) root system decomposition on edges
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict

PHI = (1 + 5**0.5) / 2

def generate_600cell():
    phi = PHI
    vertices = set()
    for i in range(4):
        for s in [1, -1]:
            v = [0.0]*4; v[i] = float(s)
            vertices.add(tuple(round(x,10) for x in v))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    vertices.add((s0, s1, s2, s3))
    base_vals = [phi/2, 0.5, 1/(2*phi), 0.0]
    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
        if inv % 2 == 0:
            even_perms.append(p)
    for perm in even_perms:
        for s0 in [1,-1]:
            for s1 in [1,-1]:
                for s2 in [1,-1]:
                    signed = [s0*base_vals[0], s1*base_vals[1], s2*base_vals[2], base_vals[3]]
                    v = tuple(round(signed[perm.index(i)],10) for i in range(4))
                    vertices.add(v)
    return [list(v) for v in sorted(vertices)]

def classify_vertex(v):
    nz = sum(1 for x in v if abs(x) > 1e-9)
    if nz == 1 and any(abs(abs(x)-1.0) < 1e-9 for x in v):
        return 'A'
    if all(abs(abs(x)-0.5) < 1e-9 for x in v):
        return 'B'
    return 'C'

print("="*70)
print("EXP-143: EDGE-BASED GAUGE FIX (1,2,9) -> (1,3,8)")
print("="*70)

verts = generate_600cell()
N = len(verts)
vtypes = [classify_vertex(v) for v in verts]
coords = np.array(verts)

neighbors = defaultdict(set)
for i in range(N):
    for j in range(i+1, N):
        d2 = np.sum((coords[i] - coords[j])**2)
        if d2 < 1.0/PHI**2 + 1e-6:
            neighbors[i].add(j)
            neighbors[j].add(i)

A_verts = sorted([i for i in range(N) if vtypes[i] == 'A'])
B_verts = sorted([i for i in range(N) if vtypes[i] == 'B'])
C_verts = sorted([i for i in range(N) if vtypes[i] == 'C'])
print(f"A={len(A_verts)}, B={len(B_verts)}, C={len(C_verts)}")

# ============================================================
# ANALYSIS 1: Detailed local CC structure
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 1: LOCAL CC NEIGHBORHOOD STRUCTURE")
print("="*70)

C_set = set(C_verts)
A_set = set(A_verts)
B_set = set(B_verts)

# For each C vertex, get its 9 CC neighbors, 1 A neighbor, 2 B neighbors
for ci_idx, ci in enumerate(C_verts[:1]):  # Analyze first C vertex in detail
    cc_nbrs = sorted([j for j in neighbors[ci] if j in C_set])
    ac_nbrs = sorted([j for j in neighbors[ci] if j in A_set])
    bc_nbrs = sorted([j for j in neighbors[ci] if j in B_set])
    print(f"\nC vertex {ci}: {len(ac_nbrs)} A-nbrs, {len(bc_nbrs)} B-nbrs, {len(cc_nbrs)} CC-nbrs")

    # Local CC subgraph
    cc_map = {v: i for i, v in enumerate(cc_nbrs)}
    local_adj = np.zeros((9, 9))
    for a in range(9):
        for b in range(a+1, 9):
            if cc_nbrs[b] in neighbors[cc_nbrs[a]]:
                local_adj[a][b] = 1
                local_adj[b][a] = 1

    local_degs = np.sum(local_adj, axis=1).astype(int)
    print(f"  Local degrees: {sorted(local_degs)} (6 x deg-3, 3 x deg-4)")

    # Which CC neighbors connect to which gauge vertices?
    print("\n  CC neighbor gauge connections:")
    for idx, cn in enumerate(cc_nbrs):
        a_conn = [a for a in A_verts if a in neighbors[cn]]
        b_conn = [b for b in B_verts if b in neighbors[cn]]
        # Does cn share the SAME A neighbor as ci?
        shared_A = set(ac_nbrs) & set(a_conn)
        shared_B = set(bc_nbrs) & set(b_conn)
        deg = local_degs[idx]
        print(f"    CC[{idx}] (deg={deg}): A-shared={len(shared_A)}, B-shared={len(shared_B)}, "
              f"A-total={len(a_conn)}, B-total={len(b_conn)}")

# Do this analysis for ALL C vertices and check uniformity
print("\n--- Systematic analysis across all C vertices ---")

# For each C vertex, classify its 9 CC neighbors by:
# (local_degree, shared_A_count, shared_B_count)
pattern_stats = []
for ci in C_verts:
    cc_nbrs = sorted([j for j in neighbors[ci] if j in C_set])
    ac_nbrs = set([j for j in neighbors[ci] if j in A_set])
    bc_nbrs = set([j for j in neighbors[ci] if j in B_set])

    # Local adjacency
    local_adj_temp = np.zeros((9, 9))
    for a in range(9):
        for b in range(a+1, 9):
            if cc_nbrs[b] in neighbors[cc_nbrs[a]]:
                local_adj_temp[a][b] = 1
                local_adj_temp[b][a] = 1
    local_degs = np.sum(local_adj_temp, axis=1).astype(int)

    nbr_profiles = []
    for idx, cn in enumerate(cc_nbrs):
        a_conn = set([a for a in A_verts if a in neighbors[cn]])
        b_conn = set([b for b in B_verts if b in neighbors[cn]])
        shared_A = len(ac_nbrs & a_conn)
        shared_B = len(bc_nbrs & b_conn)
        nbr_profiles.append((int(local_degs[idx]), shared_A, shared_B))

    profile_counts = Counter(nbr_profiles)
    pattern_stats.append(tuple(sorted(profile_counts.items())))

pattern_dist = Counter(pattern_stats)
print(f"Unique neighborhood patterns: {len(pattern_dist)}")
for pat, cnt in pattern_dist.most_common(5):
    print(f"  x{cnt}: {dict(pat)}")

# ============================================================
# ANALYSIS 2: Splitting CC edges by shared gauge structure
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 2: CC EDGE CLASSIFICATION BY GAUGE OVERLAP")
print("="*70)

# For each CC edge (i,j), count:
# - shared A neighbors
# - shared B neighbors
# - total shared gauge neighbors
cc_edge_profiles = []
for ci in C_verts:
    for cj in neighbors[ci]:
        if cj in C_set and cj > ci:
            a_i = set(a for a in A_verts if a in neighbors[ci])
            a_j = set(a for a in A_verts if a in neighbors[cj])
            b_i = set(b for b in B_verts if b in neighbors[ci])
            b_j = set(b for b in B_verts if b in neighbors[cj])
            shared_A = len(a_i & a_j)
            shared_B = len(b_i & b_j)
            cc_edge_profiles.append((shared_A, shared_B))

edge_profile_counts = Counter(cc_edge_profiles)
print(f"\nCC edge profiles (shared_A, shared_B):")
for prof, cnt in sorted(edge_profile_counts.items()):
    total_shared = prof[0] + prof[1]
    print(f"  shared_A={prof[0]}, shared_B={prof[1]} (total={total_shared}): {cnt} edges")

# Per C vertex: how many of its 9 CC edges have each profile?
print("\n--- Per vertex distribution ---")
per_vertex_edge_profiles = []
for ci in C_verts:
    profiles_at_ci = []
    for cj in neighbors[ci]:
        if cj in C_set:
            a_i = set(a for a in A_verts if a in neighbors[ci])
            a_j = set(a for a in A_verts if a in neighbors[cj])
            b_i = set(b for b in B_verts if b in neighbors[ci])
            b_j = set(b for b in B_verts if b in neighbors[cj])
            shared_A = len(a_i & a_j)
            shared_B = len(b_i & b_j)
            profiles_at_ci.append((shared_A, shared_B))
    per_vertex_edge_profiles.append(Counter(profiles_at_ci))

# Check uniformity
pvep_tuples = [tuple(sorted(c.items())) for c in per_vertex_edge_profiles]
pvep_dist = Counter(pvep_tuples)
print(f"Unique per-vertex distributions: {len(pvep_dist)}")
for pat, cnt in pvep_dist.most_common(5):
    print(f"  x{cnt}: {dict(pat)}")

# ============================================================
# ANALYSIS 3: Can we split 9 CC = 8 + 1?
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 3: NATURAL 8+1 SPLIT OF CC EDGES")
print("="*70)

# From analysis 2, if there are CC edges with different gauge overlap,
# we might split them into "SU(3)-like" (8) and "extra" (1)

# Or: the 3 degree-4 CC neighbors (in local graph) vs 6 degree-3
# Maybe the edge to one of the degree-4 neighbors should be reclassified?

# Check: for each C vertex, do the 3 degree-4 CC neighbors have different
# gauge connections than the 6 degree-3 ones?
print("\n--- Degree-4 vs Degree-3 CC neighbors ---")
deg4_profiles = []
deg3_profiles = []

for ci in C_verts:
    cc_nbrs = sorted([j for j in neighbors[ci] if j in C_set])
    ac_nbrs = set([j for j in neighbors[ci] if j in A_set])
    bc_nbrs = set([j for j in neighbors[ci] if j in B_set])

    local_adj_temp = np.zeros((9, 9))
    for a in range(9):
        for b in range(a+1, 9):
            if cc_nbrs[b] in neighbors[cc_nbrs[a]]:
                local_adj_temp[a][b] = 1
                local_adj_temp[b][a] = 1
    local_degs = np.sum(local_adj_temp, axis=1).astype(int)

    for idx, cn in enumerate(cc_nbrs):
        a_conn = set(a for a in A_verts if a in neighbors[cn])
        b_conn = set(b for b in B_verts if b in neighbors[cn])
        shared_A = len(ac_nbrs & a_conn)
        shared_B = len(bc_nbrs & b_conn)
        if local_degs[idx] == 4:
            deg4_profiles.append((shared_A, shared_B))
        else:
            deg3_profiles.append((shared_A, shared_B))

print(f"Degree-4 CC neighbors (3 per vertex, {len(deg4_profiles)} total):")
print(f"  Profiles: {Counter(deg4_profiles)}")
print(f"Degree-3 CC neighbors (6 per vertex, {len(deg3_profiles)} total):")
print(f"  Profiles: {Counter(deg3_profiles)}")

# ============================================================
# ANALYSIS 4: BC edges - can we split 2 BC = 2 + 1 (adding 1 to SU(2))?
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 4: BC EDGE STRUCTURE")
print("="*70)

# Each C has 2 B neighbors. Can one of them be "more SU(2)-like"?
# Or: can we reassign 1 CC edge as "weak" to get 1+3+8?

# Check the B-B structure through C
print("B neighbor pairs per C vertex:")
b_pair_profiles = []
for ci in C_verts:
    b_nbrs = sorted([j for j in neighbors[ci] if j in B_set])
    assert len(b_nbrs) == 2
    b1, b2 = b_nbrs
    # Are b1 and b2 connected to each other?
    b1b2_connected = b2 in neighbors[b1]
    # How many C neighbors do b1 and b2 share (besides ci)?
    b1_c = set(c for c in neighbors[b1] if c in C_set)
    b2_c = set(c for c in neighbors[b2] if c in C_set)
    shared_c = len(b1_c & b2_c)
    b_pair_profiles.append((b1b2_connected, shared_c))

bp_counts = Counter(b_pair_profiles)
print(f"  Profiles (B1-B2 connected, shared C neighbors):")
for prof, cnt in sorted(bp_counts.items()):
    print(f"    connected={prof[0]}, shared_C={prof[1]}: {cnt}")

# ============================================================
# ANALYSIS 5: Triangle-based splitting
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 5: TRIANGLE-BASED EDGE CLASSIFICATION")
print("="*70)

# For each C vertex, count triangles involving each of its edges
# ACC triangles: edge AC shared with another C
# BCC triangles: edge BC shared with another C
# CCC triangles: edge CC shared with another C

for ci in C_verts[:1]:
    print(f"\nDetailed triangle analysis for C vertex {ci}:")
    all_nbrs = sorted(neighbors[ci])
    for nj in all_nbrs:
        nt = vtypes[nj]
        # Triangles through edge ci-nj
        shared = neighbors[ci] & neighbors[nj]
        shared_types = Counter(vtypes[s] for s in shared)
        n_tri = len(shared)
        print(f"  Edge {ci}-{nj} (type {nt}): {n_tri} triangles, shared types: {dict(shared_types)}")

# Systematic: for each C vertex, classify edges by triangle count
print("\n--- Triangle count per edge type ---")
edge_triangle_data = {'AC': [], 'BC': [], 'CC': []}
for ci in C_verts:
    for nj in neighbors[ci]:
        nt = vtypes[nj]
        shared = neighbors[ci] & neighbors[nj]
        n_tri = len(shared)
        if nt == 'A':
            edge_triangle_data['AC'].append(n_tri)
        elif nt == 'B':
            edge_triangle_data['BC'].append(n_tri)
        else:
            edge_triangle_data['CC'].append(n_tri)

for et in ['AC', 'BC', 'CC']:
    data = edge_triangle_data[et]
    print(f"  {et} edges: triangle counts = {Counter(data)}")

# ============================================================
# ANALYSIS 6: Spectral splitting of local neighborhood
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 6: SPECTRAL SPLITTING OF 12-NEIGHBORHOOD")
print("="*70)

# For a C vertex, build the FULL 12-neighborhood adjacency
# and find eigenspaces that naturally split into 1+3+8
ci = C_verts[0]
nbrs = sorted(neighbors[ci])
assert len(nbrs) == 12

# Build 12x12 local adjacency
nbr_map = {v: i for i, v in enumerate(nbrs)}
local_full = np.zeros((12, 12))
for a in range(12):
    for b in range(a+1, 12):
        if nbrs[b] in neighbors[nbrs[a]]:
            local_full[a][b] = 1
            local_full[b][a] = 1

local_types = [vtypes[n] for n in nbrs]
print(f"Local 12-neighborhood types: {Counter(local_types)}")

local_eigs, local_vecs = np.linalg.eigh(local_full)
idx_sort = np.argsort(local_eigs)[::-1]
local_eigs = local_eigs[idx_sort]
local_vecs = local_vecs[:, idx_sort]

from collections import OrderedDict
def group_eigs_with_vecs(eigs, vecs, tol=0.01):
    groups = []
    i = 0
    while i < len(eigs):
        val = eigs[i]
        mult = 1
        while i + mult < len(eigs) and abs(eigs[i+mult] - val) < tol:
            mult += 1
        groups.append((round(val, 4), mult, vecs[:, i:i+mult]))
        i += mult
    return groups

local_groups = group_eigs_with_vecs(local_eigs, local_vecs)
print(f"\nLocal 12-neighborhood spectrum:")
for val, mult, vecs in local_groups:
    # Weight on each type
    a_idx = [i for i, t in enumerate(local_types) if t == 'A']
    b_idx = [i for i, t in enumerate(local_types) if t == 'B']
    c_idx = [i for i, t in enumerate(local_types) if t == 'C']
    a_w = np.sum(vecs[a_idx, :]**2) if a_idx else 0
    b_w = np.sum(vecs[b_idx, :]**2) if b_idx else 0
    c_w = np.sum(vecs[c_idx, :]**2) if c_idx else 0
    total = a_w + b_w + c_w
    print(f"  {val:7.3f} x {mult}  A={a_w/total:.3f} B={b_w/total:.3f} C={c_w/total:.3f}")

# Check: are there eigenspaces with dimension 1, 3, 8?
print("\n  Eigenspace dimensions:", [mult for _, mult, _ in local_groups])
dims = [mult for _, mult, _ in local_groups]
# Check if any combination gives 1+3+8
for i in range(len(dims)):
    for j in range(i+1, len(dims)):
        for k in range(j+1, len(dims)):
            if sorted([dims[i], dims[j], dims[k]]) == [1, 3, 8]:
                print(f"  FOUND 1+3+8 split at eigenspaces {i},{j},{k}!")
                print(f"    Eigenvalues: {local_groups[i][0]}, {local_groups[j][0]}, {local_groups[k][0]}")

# Also check 1+3+8 from consecutive groups
for i in range(len(dims)):
    cumsum = 0
    parts = []
    for j in range(i, len(dims)):
        cumsum += dims[j]
        parts.append(dims[j])
        if cumsum in [1, 3, 4, 8, 11, 12]:
            pass  # possible partial match

# ============================================================
# ANALYSIS 7: D4 root system on edges
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 7: D4 ROOT STRUCTURE ON EDGES")
print("="*70)

# The 24 gauge vertices form the D4 root system (24-cell)
# D4 = SU(3) x SU(2) x U(1) decomposition:
# D4 roots can be split into A2 (SU(3)) + A1 (SU(2)) + U(1) subalgebras
# This gives 6 + 2 + 1 = 9 positive roots... no, that's wrong

# D4 has rank 4, 24 roots (12 positive + 12 negative)
# Decomposition D4 -> A2 + A1 + U(1):
# dim = 28 -> 8 + 3 + 1 = 12 (adjoint) + 16 (off-diagonal)
# Wait, D4 has dim 28 as Lie algebra, but 24 roots + 4 Cartan = 28

# The edge assignment should use the D4 REPRESENTATION theory
# Each edge from gauge vertex g to fermion c transforms under
# the representation of the gauge group at g

# But simpler: how many distinct gauge vertices connect to a given edge?
# From exp142: each CC edge has gauge charge = 6 (uniform)
# With 1 or 2 gauge vertices connecting to BOTH endpoints

# Let's look at gauge vertex PAIRS that share C-C edges
print("Gauge vertex pairs sharing CC edges:")
gauge_cc_pairs = defaultdict(int)
for ci in C_verts:
    for cj in neighbors[ci]:
        if cj in C_set and cj > ci:
            # Which gauge vertices connect to BOTH ci and cj?
            for g in A_verts + B_verts:
                if ci in neighbors[g] and cj in neighbors[g]:
                    gauge_cc_pairs[g] += 1

print("Per gauge vertex, CC edges fully covered:")
for g in A_verts + B_verts:
    gt = vtypes[g]
    cnt = gauge_cc_pairs[g]
    print(f"  {gt} vertex {g}: {cnt} CC edges")

# ============================================================
# ANALYSIS 8: Icosahedron vertex figure
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 8: ICOSAHEDRON VERTEX FIGURE AND EDGE TYPES")
print("="*70)

# The vertex figure of the 600-cell is an icosahedron (12 vertices)
# The 12 neighbors of a C vertex form an icosahedron
# In this icosahedron: 1A + 2B + 9C vertices

# The icosahedron has 30 edges
# Question: how do these 30 edges decompose by type?
ci = C_verts[0]
nbrs = sorted(neighbors[ci])
ico_edges = []
for a in range(12):
    for b in range(a+1, 12):
        if nbrs[b] in neighbors[nbrs[a]]:
            ico_edges.append((a, b))

print(f"Icosahedron edges: {len(ico_edges)} (expected 30)")
ico_types = [vtypes[n] for n in nbrs]

ico_edge_types = Counter()
for a, b in ico_edges:
    key = tuple(sorted([ico_types[a], ico_types[b]]))
    ico_edge_types[key] += 1

print("Icosahedron edge types:")
for key, cnt in sorted(ico_edge_types.items()):
    print(f"  {key[0]}-{key[1]}: {cnt}")

# The A5 symmetry of icosahedron: irreps 1, 3, 3', 4, 5
# How does the vertex type assignment break this symmetry?
print(f"\nVertex types in icosahedron: {Counter(ico_types)}")
print(f"  A at positions: {[i for i, t in enumerate(ico_types) if t == 'A']}")
print(f"  B at positions: {[i for i, t in enumerate(ico_types) if t == 'B']}")

# Adjacency between type groups in the icosahedron
a_icos = [i for i, t in enumerate(ico_types) if t == 'A']
b_icos = [i for i, t in enumerate(ico_types) if t == 'B']
c_icos = [i for i, t in enumerate(ico_types) if t == 'C']

print(f"\n  A-C edges in ico: {sum(1 for a,b in ico_edges if (a in a_icos and b in c_icos) or (b in a_icos and a in c_icos))}")
print(f"  B-C edges in ico: {sum(1 for a,b in ico_edges if (a in b_icos and b in c_icos) or (b in b_icos and a in c_icos))}")
print(f"  C-C edges in ico: {sum(1 for a,b in ico_edges if a in c_icos and b in c_icos)}")
print(f"  A-A edges in ico: {sum(1 for a,b in ico_edges if a in a_icos and b in a_icos)}")
print(f"  A-B edges in ico: {sum(1 for a,b in ico_edges if (a in a_icos and b in b_icos) or (b in a_icos and a in b_icos))}")
print(f"  B-B edges in ico: {sum(1 for a,b in ico_edges if a in b_icos and b in b_icos)}")

# The A vertex in the icosahedron has 5 neighbors (icosahedron degree)
# Similarly B vertices have 5 neighbors each
for idx in a_icos:
    a_nbrs = [b for a, b in ico_edges if a == idx] + [a for a, b in ico_edges if b == idx]
    a_nbr_types = Counter(ico_types[n] for n in a_nbrs)
    print(f"\n  A-vertex {idx} icosahedron neighbors: {dict(a_nbr_types)}")

for idx in b_icos:
    b_nbrs = [b for a, b in ico_edges if a == idx] + [a for a, b in ico_edges if b == idx]
    b_nbr_types = Counter(ico_types[n] for n in b_nbrs)
    print(f"  B-vertex {idx} icosahedron neighbors: {dict(b_nbr_types)}")

# ============================================================
# ANALYSIS 9: Refined split using icosahedron structure
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 9: REFINED EDGE SPLIT")
print("="*70)

# In icosahedron: 5 AC + 10 BC + 15 CC = 30 edges
# These correspond to edges of 600-cell between pairs of ci's neighbors
# NOT to edges FROM ci

# The actual edges from ci are:
# 1 to A, 2 to B, 9 to C = 12 edges (spokes of the icosahedron)
# The icosahedron edges (30) are between neighbors (rim)

# For the 12 spokes: can we classify them using the icosahedron geometry?
# E.g., the A-spoke and the 2 B-spokes partition the icosahedron

# Check: the 5 icosahedron neighbors of A in the icosahedron
# do they correspond to a "pentagonal" structure?
print("Icosahedron structure around the A vertex:")
a_idx = a_icos[0]
a_ico_nbrs = sorted([b for a, b in ico_edges if a == a_idx] + [a for a, b in ico_edges if b == a_idx])
print(f"  A at pos {a_idx}, its 5 ico-neighbors: {a_ico_nbrs}")
print(f"  Their types: {[ico_types[n] for n in a_ico_nbrs]}")

# The 5 neighbors of A in icosahedron form a pentagon
# Check connectivity among them
pent_adj = []
for i in range(5):
    for j in range(i+1, 5):
        if a_ico_nbrs[j] in [b for a, b in ico_edges if a == a_ico_nbrs[i]] + \
                            [a for a, b in ico_edges if b == a_ico_nbrs[i]]:
            pent_adj.append((i, j))
print(f"  Pentagon edges: {pent_adj} (expected 5)")

# The opposite pentagonal antiprism: the 5 vertices NOT adjacent to A (and not B)
# Wait, we have 1A + 2B + 9C = 12 vertices
# A's ico-neighbors: 5 vertices.
# Non-neighbors of A in ico: 12 - 1 - 5 = 6 vertices
non_a_nbrs = [i for i in range(12) if i != a_idx and i not in a_ico_nbrs]
print(f"  Non-neighbors of A: {non_a_nbrs}")
print(f"  Their types: {[ico_types[n] for n in non_a_nbrs]}")

# ============================================================
# ANALYSIS 10: The 1+2+9 -> 1+3+8 reassignment hypothesis
# ============================================================
print("\n" + "="*70)
print("ANALYSIS 10: REASSIGNMENT HYPOTHESIS")
print("="*70)

print("""
Current:  1 AC + 2 BC + 9 CC = 12
SM:       1 U(1) + 3 SU(2) + 8 SU(3) = 12

Hypothesis: Reclassify 1 CC edge as "weak" (SU(2)-type)
This gives: 1 + (2+1) + (9-1) = 1 + 3 + 8

Question: which CC edge should be reclassified?
Candidates:
  A) The CC edge to the degree-4 neighbor closest to B vertices
  B) The CC edge with the most shared B-gauge vertices
  C) The CC edge in the A5-pentagon that breaks the most symmetry
""")

# Test: among the 9 CC neighbors, which one is most "B-like"?
# Criterion: closest to B vertices in graph distance
from collections import deque

def bfs_dist_to_set(start, target_set, adj_dict):
    """BFS distance from start to nearest vertex in target_set"""
    visited = {start}
    queue = deque([(start, 0)])
    while queue:
        v, d = queue.popleft()
        for nb in adj_dict[v]:
            if nb in target_set:
                return d + 1
            if nb not in visited:
                visited.add(nb)
                queue.append((nb, d + 1))
    return -1

# For the first C vertex, find which CC neighbor is "closest" to B
ci = C_verts[0]
cc_nbrs = sorted([j for j in neighbors[ci] if j in C_set])
print(f"\nFor C vertex {ci}, CC neighbor analysis:")

for idx, cn in enumerate(cc_nbrs):
    # Count B-neighbors of cn
    b_nbrs_cn = [b for b in B_verts if b in neighbors[cn]]
    # Count shared B with ci
    b_ci = set(b for b in B_verts if b in neighbors[ci])
    shared_b = len(set(b_nbrs_cn) & b_ci)
    # Count total gauge overlap
    g_ci = set(g for g in A_verts + B_verts if g in neighbors[ci])
    g_cn = set(g for g in A_verts + B_verts if g in neighbors[cn])
    overlap = len(g_ci & g_cn)
    print(f"  CC[{idx}]: B-nbrs={len(b_nbrs_cn)}, shared-B={shared_b}, gauge-overlap={overlap}")

# Summary: is there a unique "special" CC edge per C vertex?
print("\n--- Looking for unique special CC edge ---")
# Check if exactly one CC neighbor per C has maximum gauge overlap
max_overlap_counts = []
for ci in C_verts:
    cc_nbrs = sorted([j for j in neighbors[ci] if j in C_set])
    overlaps = []
    for cn in cc_nbrs:
        g_ci = set(g for g in A_verts + B_verts if g in neighbors[ci])
        g_cn = set(g for g in A_verts + B_verts if g in neighbors[cn])
        overlaps.append(len(g_ci & g_cn))
    max_ov = max(overlaps)
    n_max = overlaps.count(max_ov)
    max_overlap_counts.append((max_ov, n_max))

moc_dist = Counter(max_overlap_counts)
print(f"(max_overlap, count_at_max): {dict(sorted(moc_dist.items()))}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"\nCurrent split: 1 AC + 2 BC + 9 CC = 12")
print(f"Target: 1 + 3 + 8 = 12")
print(f"Gap: need to move 1 edge from CC to 'weak' sector")

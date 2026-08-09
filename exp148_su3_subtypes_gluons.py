"""
EXP-148: SU(3) Sub-types and the 8 Gluons
============================================
From exp147: SU3 edges have 3 sub-types:
  (0A, 1B, 4C): 96 edges, gauge_access=1
  (0A, 2B, 3C): 96 edges, gauge_access=2
  (1A, 1B, 3C): 192 edges, gauge_access=2

Per fermion: 2+2+4 = 8 SU3 edges (ratio 1:1:2)
Total: 384 = 8 * 48

Questions:
1. Can the 3 sub-types map to gluon structure (adjoint SU(3))?
2. Can we find a 3-coloring of fermions from SU3 edges?
3. Do the sub-types form recognizable algebraic structures?
4. How do pure SU3 plaquettes (288) decompose across sub-types?
"""

import numpy as np
from itertools import permutations
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
print("EXP-148: SU(3) SUB-TYPES AND THE 8 GLUONS")
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

A_set = set(i for i in range(N) if vtypes[i] == 'A')
B_set = set(i for i in range(N) if vtypes[i] == 'B')
C_set = set(i for i in range(N) if vtypes[i] == 'C')
C_list = sorted(C_set)

all_edges = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            all_edges.append((i,j))

# Classify edges with fine types
def get_edge_fine_type(i, j):
    ti, tj = vtypes[i], vtypes[j]
    if (ti == 'A' and tj == 'C') or (ti == 'C' and tj == 'A'):
        return 'U1'
    if (ti == 'B' and tj == 'C') or (ti == 'C' and tj == 'B'):
        return 'SU2_ch'
    if ti == 'C' and tj == 'C':
        shared_A = len((A_set & neighbors[i]) & (A_set & neighbors[j]))
        shared_B = len((B_set & neighbors[i]) & (B_set & neighbors[j]))
        if shared_A > 0 and shared_B == 0:
            return 'SU2_n'
        # SU3 sub-types based on shared neighbor profile
        shared_C = len((C_set & neighbors[i]) & (C_set & neighbors[j]))
        return f'SU3_{shared_A}A{shared_B}B{shared_C}C'
    return 'NONE'

edge_fine = {}
for i, j in all_edges:
    ft = get_edge_fine_type(i, j)
    edge_fine[(i,j)] = ft
    edge_fine[(j,i)] = ft

fine_counts = Counter(edge_fine[(i,j)] for i,j in all_edges)
print("Fine edge classification:")
for ft, cnt in sorted(fine_counts.items()):
    print(f"  {ft}: {cnt}")

# ============================================================
# PART 1: Per-fermion SU3 sub-type decomposition
# ============================================================
print("\n" + "="*70)
print("PART 1: PER-FERMION SU3 EDGE DECOMPOSITION")
print("="*70)

# For each fermion, list its 8 SU3 edges by sub-type
per_fermion_su3 = {}
for ci in C_list:
    su3_edges = {}
    for j in neighbors[ci]:
        ft = edge_fine[(ci, j)]
        if ft.startswith('SU3'):
            su3_edges[j] = ft
    per_fermion_su3[ci] = su3_edges

# Check uniformity
su3_profiles = []
for ci in C_list:
    profile = Counter(per_fermion_su3[ci].values())
    su3_profiles.append(tuple(sorted(profile.items())))

profile_counts = Counter(su3_profiles)
print(f"Unique per-fermion SU3 profiles: {len(profile_counts)}")
for prof, cnt in profile_counts.items():
    print(f"  x{cnt}: {dict(prof)}")

# ============================================================
# PART 2: SU3 sub-graph structure
# ============================================================
print("\n" + "="*70)
print("PART 2: SU3 SUB-TYPE GRAPH STRUCTURE")
print("="*70)

# Build sub-graphs for each SU3 sub-type
su3_types = ['SU3_0A1B4C', 'SU3_0A2B3C', 'SU3_1A1B3C']
for st in su3_types:
    edges_st = [(i,j) for i,j in all_edges if edge_fine[(i,j)] == st]
    verts_st = set()
    for i,j in edges_st:
        verts_st.add(i)
        verts_st.add(j)

    # Vertex type distribution
    v_types = Counter(vtypes[v] for v in verts_st)

    # Degree distribution within this sub-type
    deg = defaultdict(int)
    for i,j in edges_st:
        deg[i] += 1
        deg[j] += 1

    c_degs = [deg[v] for v in verts_st if vtypes[v] == 'C']
    deg_dist = Counter(c_degs)

    print(f"\n{st} ({len(edges_st)} edges):")
    print(f"  Vertices: {len(verts_st)} ({dict(v_types)})")
    print(f"  C-vertex degrees: {dict(deg_dist)}")

    # Check if sub-graph is connected
    if verts_st:
        adj_st = defaultdict(set)
        for i,j in edges_st:
            adj_st[i].add(j)
            adj_st[j].add(i)
        # BFS
        start = next(iter(verts_st))
        visited = {start}
        queue = [start]
        while queue:
            v = queue.pop(0)
            for nb in adj_st[v]:
                if nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
        n_components = 1
        remaining = verts_st - visited
        while remaining:
            n_components += 1
            start = next(iter(remaining))
            visited2 = {start}
            queue = [start]
            while queue:
                v = queue.pop(0)
                for nb in adj_st[v]:
                    if nb not in visited2:
                        visited2.add(nb)
                        queue.append(nb)
            remaining -= visited2
        print(f"  Connected components: {n_components}")

# ============================================================
# PART 3: 3-coloring attempt from SU3 edges
# ============================================================
print("\n" + "="*70)
print("PART 3: 3-COLORING FROM SU3 STRUCTURE")
print("="*70)

# In QCD, 3 colors = fundamental rep of SU(3)
# Can we 3-color the 96 C-vertices using SU3 edge structure?

# Approach 1: Use the 3 sub-types to define a coloring
# Each fermion has 2+2+4 SU3 edges of types T1, T2, T3
# The T3 edges (192 total, 4 per fermion) are the most connected
# Try: chromatic number of the SU3 sub-graph

# Build full SU3 adjacency (CC edges that are SU3 type)
su3_adj = defaultdict(set)
for i,j in all_edges:
    if edge_fine[(i,j)].startswith('SU3'):
        su3_adj[i].add(j)
        su3_adj[j].add(i)

# Each C vertex has degree 8 in SU3 subgraph
su3_degrees = [len(su3_adj[c]) for c in C_list]
print(f"SU3 subgraph: {len(C_list)} vertices, degree {Counter(su3_degrees)}")

# Chromatic number lower bound: max clique
# Check triangles (3-cliques) in SU3 subgraph
su3_triangles = 0
for i in C_list:
    for j in su3_adj[i]:
        if j > i:
            for k in su3_adj[i] & su3_adj[j]:
                if k > j:
                    su3_triangles += 1
print(f"SU3 triangles: {su3_triangles}")

# 4-cliques in SU3 subgraph
su3_4cliques = 0
for i in C_list:
    ni = su3_adj[i]
    for j in ni:
        if j > i:
            nij = ni & su3_adj[j]
            for k in nij:
                if k > j:
                    nijk = nij & su3_adj[k]
                    for l in nijk:
                        if l > k:
                            su3_4cliques += 1
print(f"SU3 4-cliques: {su3_4cliques}")

# Greedy coloring
import random
random.seed(42)
def greedy_color(adj, vertices):
    colors = {}
    for v in vertices:
        used = {colors[nb] for nb in adj[v] if nb in colors}
        c = 0
        while c in used:
            c += 1
        colors[v] = c
    return colors

colors = greedy_color(su3_adj, C_list)
n_colors = max(colors.values()) + 1
color_dist = Counter(colors.values())
print(f"\nGreedy coloring: {n_colors} colors needed")
print(f"  Color distribution: {dict(color_dist)}")

# Try with random orderings
min_colors = n_colors
for trial in range(100):
    shuffled = list(C_list)
    random.shuffle(shuffled)
    cols = greedy_color(su3_adj, shuffled)
    nc = max(cols.values()) + 1
    if nc < min_colors:
        min_colors = nc
print(f"Best greedy coloring (100 trials): {min_colors} colors")

# ============================================================
# PART 4: Sub-type adjacency (which sub-types connect?)
# ============================================================
print("\n" + "="*70)
print("PART 4: SUB-TYPE ADJACENCY (WHICH TYPES CONNECT?)")
print("="*70)

# For each pair of SU3 edges that share a vertex,
# what are the sub-type combinations?
pair_counts = defaultdict(int)
for ci in C_list:
    su3_nbrs = [(j, edge_fine[(ci,j)]) for j in neighbors[ci]
                if edge_fine[(ci,j)].startswith('SU3')]
    for idx1 in range(len(su3_nbrs)):
        for idx2 in range(idx1+1, len(su3_nbrs)):
            j1, t1 = su3_nbrs[idx1]
            j2, t2 = su3_nbrs[idx2]
            pair = tuple(sorted([t1, t2]))
            # Are j1 and j2 also connected by an SU3 edge?
            if j2 in su3_adj[j1]:
                pair_counts[('triangle', pair)] += 1
            else:
                pair_counts[('open', pair)] += 1

print("SU3 edge pair types (at shared vertex):")
print("\nTriangles (both endpoints connected):")
for (kind, pair), cnt in sorted(pair_counts.items()):
    if kind == 'triangle':
        print(f"  {pair[0]} -- {pair[1]}: {cnt}")

print("\nOpen pairs (endpoints NOT connected):")
for (kind, pair), cnt in sorted(pair_counts.items()):
    if kind == 'open':
        print(f"  {pair[0]} -- {pair[1]}: {cnt}")

# ============================================================
# PART 5: Pure SU3 plaquettes by sub-type
# ============================================================
print("\n" + "="*70)
print("PART 5: PURE SU3 PLAQUETTES BY SUB-TYPE")
print("="*70)

triangles = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            for k in neighbors[i] & neighbors[j]:
                if k > j:
                    triangles.append((i,j,k))

# Find pure SU3 triangles
pure_su3_tris = []
for i, j, k in triangles:
    types = [edge_fine[(i,j)], edge_fine[(i,k)], edge_fine[(j,k)]]
    if all(t.startswith('SU3') for t in types):
        pure_su3_tris.append((i,j,k, tuple(sorted(types))))

print(f"Pure SU3 triangles: {len(pure_su3_tris)}")

# Decompose by sub-type combination
tri_type_counts = Counter(t[3] for t in pure_su3_tris)
print("\nPure SU3 triangle types:")
for combo, cnt in sorted(tri_type_counts.items()):
    labels = [c.replace('SU3_', '') for c in combo]
    print(f"  {' -- '.join(labels)}: {cnt}")

# ============================================================
# PART 6: Each fermion's 8 SU3 neighbors - structure
# ============================================================
print("\n" + "="*70)
print("PART 6: LOCAL SU3 NEIGHBORHOOD STRUCTURE")
print("="*70)

# For one fermion, show the complete SU3 neighborhood
ci = C_list[0]
print(f"Fermion {ci}, SU3 neighbors:")
su3_nbrs_ci = []
for j in sorted(neighbors[ci]):
    ft = edge_fine[(ci, j)]
    if ft.startswith('SU3'):
        # What A and B does this neighbor share with ci?
        shared_A = sorted((A_set & neighbors[ci]) & (A_set & neighbors[j]))
        shared_B = sorted((B_set & neighbors[ci]) & (B_set & neighbors[j]))
        su3_nbrs_ci.append((j, ft, shared_A, shared_B))
        print(f"  -> {j}: {ft}, shared_A={shared_A}, shared_B={shared_B}")

# Are the 8 SU3 neighbors connected among themselves?
print(f"\nConnections among the 8 SU3 neighbors of fermion {ci}:")
su3_nbr_set = set(j for j, _, _, _ in su3_nbrs_ci)
internal_edges = 0
internal_types = defaultdict(int)
for j1 in su3_nbr_set:
    for j2 in su3_nbr_set:
        if j2 > j1 and j2 in neighbors[j1]:
            ft = edge_fine[(j1, j2)]
            internal_types[ft] += 1
            internal_edges += 1

print(f"  Internal edges: {internal_edges}")
for ft, cnt in sorted(internal_types.items()):
    print(f"    {ft}: {cnt}")

# Check uniformity across all fermions
print("\nInternal SU3 connections (all fermions):")
all_internal = []
for ci in C_list:
    su3_nbrs = [j for j in neighbors[ci] if edge_fine[(ci,j)].startswith('SU3')]
    su3_set = set(su3_nbrs)
    internal = 0
    for j1 in su3_set:
        for j2 in su3_set:
            if j2 > j1 and j2 in neighbors[j1]:
                internal += 1
    all_internal.append(internal)

print(f"  Mean internal SU3 edges per fermion: {np.mean(all_internal):.2f}")
print(f"  Std: {np.std(all_internal):.4f}")
print(f"  Min: {min(all_internal)}, Max: {max(all_internal)}")

# ============================================================
# PART 7: Graph coloring of SU3 sub-types
# ============================================================
print("\n" + "="*70)
print("PART 7: SUB-TYPE GRAPH EIGENVALUES")
print("="*70)

# Build adjacency matrices for each sub-type on C vertices
c_idx = {c: i for i, c in enumerate(C_list)}
nc = len(C_list)

for st in su3_types:
    A_st = np.zeros((nc, nc))
    edges_st = [(i,j) for i,j in all_edges if edge_fine[(i,j)] == st]
    for i, j in edges_st:
        if i in C_set and j in C_set:
            A_st[c_idx[i], c_idx[j]] = 1
            A_st[c_idx[j], c_idx[i]] = 1

    eigs = sorted(np.linalg.eigvalsh(A_st), reverse=True)

    # Group eigenvalues
    groups = []
    idx = 0
    while idx < len(eigs):
        val = eigs[idx]
        mult = 1
        while idx + mult < len(eigs) and abs(eigs[idx+mult] - val) < 0.01:
            mult += 1
        groups.append((round(val, 4), mult))
        idx += mult

    print(f"\n{st} ({len(edges_st)} edges on {nc} C-vertices):")
    print(f"  Degree per C: {2*len(edges_st)/nc:.1f}")
    for val, mult in groups:
        # Check if a+b*phi
        best = None
        for a in range(-10, 15):
            for b in range(-10, 15):
                if abs(a + b*PHI - val) < 0.01:
                    best = (a, b)
                    break
            if best: break
        label = f" = {best[0]}+{best[1]}phi" if best else ""
        print(f"  {val:8.4f} x {mult}{label}")

# Combined SU3 adjacency
A_su3 = np.zeros((nc, nc))
for i,j in all_edges:
    if edge_fine[(i,j)].startswith('SU3'):
        if i in C_set and j in C_set:
            A_su3[c_idx[i], c_idx[j]] = 1
            A_su3[c_idx[j], c_idx[i]] = 1

eigs_full = sorted(np.linalg.eigvalsh(A_su3), reverse=True)
groups_full = []
idx = 0
while idx < len(eigs_full):
    val = eigs_full[idx]
    mult = 1
    while idx + mult < len(eigs_full) and abs(eigs_full[idx+mult] - val) < 0.01:
        mult += 1
    groups_full.append((round(val, 4), mult))
    idx += mult

print(f"\nFull SU3 subgraph (384 edges, degree 8):")
for val, mult in groups_full:
    best = None
    for a in range(-10, 15):
        for b in range(-10, 15):
            if abs(a + b*PHI - val) < 0.01:
                best = (a, b)
                break
        if best: break
    label = f" = {best[0]}+{best[1]}phi" if best else ""
    print(f"  {val:8.4f} x {mult}{label}")

# ============================================================
# PART 8: B-vertex structure in SU3 sub-types
# ============================================================
print("\n" + "="*70)
print("PART 8: WHICH B-VERTICES APPEAR IN SU3 SHARED NEIGHBORS?")
print("="*70)

# For SU3 edges, the shared B-vertices are important
# Do specific B-vertices correspond to specific sub-types?

b_in_subtype = defaultdict(lambda: defaultdict(int))
for i, j in all_edges:
    ft = edge_fine[(i,j)]
    if ft.startswith('SU3'):
        shared_B = (B_set & neighbors[i]) & (B_set & neighbors[j])
        for b in shared_B:
            b_in_subtype[ft][b] += 1

for st in su3_types:
    b_counts = b_in_subtype[st]
    if b_counts:
        vals = list(b_counts.values())
        print(f"\n{st}:")
        print(f"  B-vertices involved: {len(b_counts)} of {len(B_set)}")
        print(f"  Occurrences per B: {Counter(vals)}")

# Similarly for A-vertices
a_in_subtype = defaultdict(lambda: defaultdict(int))
for i, j in all_edges:
    ft = edge_fine[(i,j)]
    if ft.startswith('SU3'):
        shared_A = (A_set & neighbors[i]) & (A_set & neighbors[j])
        for a in shared_A:
            a_in_subtype[ft][a] += 1

for st in su3_types:
    a_counts = a_in_subtype[st]
    if a_counts:
        print(f"\n{st} A-involvement:")
        print(f"  A-vertices involved: {len(a_counts)} of {len(A_set)}")
        print(f"  Occurrences per A: {Counter(a_counts.values())}")

# ============================================================
# PART 9: The 384/8 = 48 question
# ============================================================
print("\n" + "="*70)
print("PART 9: CAN WE GROUP 384 SU3 EDGES INTO 8 GROUPS OF 48?")
print("="*70)

# 384 = 8 * 48. If each "gluon" corresponds to 48 edges...
# The 48 CC* edges already give us a reference (1 per fermion)

# Approach: can we partition the 384 SU3 edges into 8 orbits?
# Under some symmetry group acting on the 600-cell

# Each fermion has 8 SU3 edges. If we label them 1-8 consistently,
# we get 8 groups of 96/2 = 48 (since each edge is shared by 2 fermions)

# But: 96 fermions * 8 edges / 2 (each edge counted twice) = 384 ✓

# Can we consistently label the 8 SU3 edges of each fermion?
# Use the sub-type + gauge overlap as a label

print("Per-fermion SU3 edge labeling attempt:")
labels_per_fermion = {}
for ci in C_list:
    labels = []
    for j in sorted(neighbors[ci]):
        ft = edge_fine[(ci, j)]
        if ft.startswith('SU3'):
            shared_A = sorted((A_set & neighbors[ci]) & (A_set & neighbors[j]))
            shared_B = sorted((B_set & neighbors[ci]) & (B_set & neighbors[j]))
            # Create a label from (sub-type, which specific A/B)
            labels.append((ft, tuple(shared_A), tuple(shared_B), j))
    labels_per_fermion[ci] = labels

# Check: are the labels consistent across fermions?
# Group by (sub-type, number of A, number of B)
label_types = Counter()
for ci in C_list:
    for ft, sa, sb, j in labels_per_fermion[ci]:
        label_types[(ft, len(sa), len(sb))] += 1

print("Label type counts (each counted from both endpoints):")
for (ft, na, nb), cnt in sorted(label_types.items()):
    print(f"  {ft} with {na}A, {nb}B: {cnt} (= {cnt//2} edges x 2)")

# ============================================================
# PART 10: The 48 = |CC*| question
# ============================================================
print("\n" + "="*70)
print("PART 10: RELATIONSHIP 48 CC* EDGES AND SU3 SUB-TYPES")
print("="*70)

# 48 CC* edges (SU2_n) - each fermion has exactly 1
# 96 SU3 type-1 edges (0A,1B,4C) - each fermion has 2
# 96 SU3 type-2 edges (0A,2B,3C) - each fermion has 2
# 192 SU3 type-3 edges (1A,1B,3C) - each fermion has 4

# Does each CC* edge pair with specific SU3 edges?
print("For each fermion, does its CC* edge determine SU3 sub-type partners?")

# For each fermion, find its CC* neighbor and its SU3 neighbors
cc_star_partner = {}
for ci in C_list:
    for j in neighbors[ci]:
        if edge_fine[(ci, j)] == 'SU2_n':
            cc_star_partner[ci] = j

# Check: are the SU3 type-3 edges (1A,1B,3C) connected to the CC* partner?
type3_connected_to_ccstar = 0
type3_total = 0
for ci in C_list:
    ccstar = cc_star_partner.get(ci)
    if ccstar is None:
        continue
    for j in neighbors[ci]:
        ft = edge_fine[(ci, j)]
        if ft == 'SU3_1A1B3C':
            type3_total += 1
            if j in neighbors[ccstar]:
                type3_connected_to_ccstar += 1

print(f"SU3 type-3 edges connected to CC* partner: {type3_connected_to_ccstar}/{type3_total}")

# What about other types?
for st in su3_types:
    connected = 0
    total = 0
    for ci in C_list:
        ccstar = cc_star_partner.get(ci)
        if ccstar is None:
            continue
        for j in neighbors[ci]:
            if edge_fine[(ci, j)] == st:
                total += 1
                if j in neighbors[ccstar]:
                    connected += 1
    print(f"  {st}: {connected}/{total} connected to CC* partner")

# ============================================================
# PART 11: Adjacency between sub-type graphs
# ============================================================
print("\n" + "="*70)
print("PART 11: COMMUTATION OF SUB-TYPE ADJACENCIES")
print("="*70)

# Build the 3 adjacency matrices for sub-types
A_t1 = np.zeros((nc, nc))  # 0A1B4C
A_t2 = np.zeros((nc, nc))  # 0A2B3C
A_t3 = np.zeros((nc, nc))  # 1A1B3C

for i, j in all_edges:
    ft = edge_fine[(i,j)]
    if i in C_set and j in C_set:
        ci, cj = c_idx[i], c_idx[j]
        if ft == 'SU3_0A1B4C':
            A_t1[ci, cj] = A_t1[cj, ci] = 1
        elif ft == 'SU3_0A2B3C':
            A_t2[ci, cj] = A_t2[cj, ci] = 1
        elif ft == 'SU3_1A1B3C':
            A_t3[ci, cj] = A_t3[cj, ci] = 1

print("Sub-type adjacency degrees:")
print(f"  T1 (0A1B4C): degree = {A_t1.sum(axis=1)[0]:.0f}")
print(f"  T2 (0A2B3C): degree = {A_t2.sum(axis=1)[0]:.0f}")
print(f"  T3 (1A1B3C): degree = {A_t3.sum(axis=1)[0]:.0f}")

# Do they commute?
comm_12 = np.max(np.abs(A_t1 @ A_t2 - A_t2 @ A_t1))
comm_13 = np.max(np.abs(A_t1 @ A_t3 - A_t3 @ A_t1))
comm_23 = np.max(np.abs(A_t2 @ A_t3 - A_t3 @ A_t2))
print(f"\nCommutators:")
print(f"  [A_t1, A_t2] max = {comm_12:.6f} ({'COMMUTE' if comm_12 < 0.01 else 'DO NOT COMMUTE'})")
print(f"  [A_t1, A_t3] max = {comm_13:.6f} ({'COMMUTE' if comm_13 < 0.01 else 'DO NOT COMMUTE'})")
print(f"  [A_t2, A_t3] max = {comm_23:.6f} ({'COMMUTE' if comm_23 < 0.01 else 'DO NOT COMMUTE'})")

# Products
print("\nProducts A_ti @ A_tj diagonal (shared neighbors through sub-type paths):")
for name, prod in [("T1*T1", A_t1@A_t1), ("T2*T2", A_t2@A_t2),
                    ("T3*T3", A_t3@A_t3), ("T1*T2", A_t1@A_t2),
                    ("T1*T3", A_t1@A_t3), ("T2*T3", A_t2@A_t3)]:
    diag = np.diag(prod)
    off = prod - np.diag(diag)
    off_vals = sorted(set(np.round(off.flatten(), 2)))
    print(f"  {name}: diag={diag[0]:.0f}, off-diag values={off_vals[:10]}")

# Does A_t1 + A_t2 + A_t3 = A_su3?
diff = np.max(np.abs(A_t1 + A_t2 + A_t3 - A_su3))
print(f"\nA_t1 + A_t2 + A_t3 == A_SU3? max diff = {diff:.10f}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print("""
SU3 SUB-TYPE STRUCTURE:
  Type T1 (0A,1B,4C): 96 edges, degree 2, gauge_access=1 (B only)
  Type T2 (0A,2B,3C): 96 edges, degree 2, gauge_access=2 (B only)
  Type T3 (1A,1B,3C): 192 edges, degree 4, gauge_access=2 (A+B mixed)

  Per fermion: 2 T1 + 2 T2 + 4 T3 = 8 SU3 edges
  Ratio: 1:1:2 (same as degree ratio 2:2:4)

  384 = 8 * 48 where 48 = |CC*| edges

KEY FINDINGS:
""")

print("="*70)
print("END EXP-148")
print("="*70)

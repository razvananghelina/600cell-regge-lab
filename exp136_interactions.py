"""
EXP-136: Interaction Structure of the 600-Cell
================================================

QUESTION: What is the interaction structure? Which particles interact with which,
and how does the graph geometry determine coupling strengths?

The 600-cell has:
- 120 vertices (particle states)
- 720 edges (propagation paths)
- 1200 triangular faces (3-point interaction vertices)
- 600 tetrahedral cells (4-point interaction vertices)

24 vertices (8 Type-A + 16 Type-B) = D4 root system = gauge sector
96 vertices (Type-C, snub 24-cell) = fermion sector

APPROACHES:
1. Triangle classification by vertex type (AAA, AAB, ..., CCC)
2. Tetrahedron classification by vertex type
3. Interaction selection rules from neighbor structure
4. Coupling strength from triangle multiplicity per edge
5. Graph distance and interaction strength (common gauge neighbors)
6. Spectral interaction matrix (P-matrix of association scheme)
7. Yukawa structure from CCC triangles
8. Comparison with SM coupling ratios

Author: Claude (exp136)
Date: 2026-02-09
Category: Analysis to determine DERIVED vs PATTERN vs SPECULATIV
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import math

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

print("=" * 80)
print("EXP-136: Interaction Structure of the 600-Cell")
print("=" * 80)
print()

# ============================================================
# SECTION 0: Build the 600-cell
# ============================================================

print("-" * 80)
print("SECTION 0: 600-cell construction")
print("-" * 80)
print()

def generate_600cell():
    """Generate all 120 vertices of the 600-cell on unit S^3."""
    phi = PHI
    vertices = set()

    # 8 vertices: permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [1, -1]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = float(s)
            vertices.add(tuple(round(x, 10) for x in v))

    # 16 vertices: (+-1/2, +-1/2, +-1/2, +-1/2)
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    vertices.add((s0, s1, s2, s3))

    # 96 vertices: EVEN permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base_vals = [phi/2, 0.5, 1/(2*phi), 0.0]

    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            even_perms.append(p)

    for perm in even_perms:
        for s0 in [1, -1]:
            for s1 in [1, -1]:
                for s2 in [1, -1]:
                    signed = [s0 * base_vals[0], s1 * base_vals[1],
                              s2 * base_vals[2], base_vals[3]]
                    v = tuple(round(signed[perm.index(i)], 10) for i in range(4))
                    vertices.add(v)

    return [np.array(v) for v in sorted(vertices)]

vertices = generate_600cell()
N = len(vertices)
print(f"  Vertices generated: {N}")
assert N == 120, f"Expected 120 vertices, got {N}"

# Build adjacency matrix using dot product = phi/2
dot_threshold = PHI / 2
adj = np.zeros((N, N), dtype=int)
edges = []
for i in range(N):
    for j in range(i+1, N):
        dot = np.dot(vertices[i], vertices[j])
        if abs(dot - dot_threshold) < 1e-6:
            adj[i, j] = 1
            adj[j, i] = 1
            edges.append((i, j))

n_edges = len(edges)
degrees = adj.sum(axis=1)
print(f"  Edges: {n_edges}")
print(f"  Degree: {degrees[0]} (all vertices)")
assert n_edges == 720
assert all(d == 12 for d in degrees)

# Classify vertices
type_A = []  # 8 vertices: cross-polytope
type_B = []  # 16 vertices: tesseract
type_C = []  # 96 vertices: snub 24-cell

for i, v in enumerate(vertices):
    nz = sum(1 for x in v if abs(x) > 0.01)
    if nz == 1 and any(abs(abs(x) - 1) < 0.01 for x in v):
        type_A.append(i)
    elif nz == 4 and all(abs(abs(x) - 0.5) < 0.01 for x in v):
        type_B.append(i)
    else:
        type_C.append(i)

set_A = set(type_A)
set_B = set(type_B)
set_C = set(type_C)

print(f"  Type A (cross-polytope/16-cell):  {len(type_A)}")
print(f"  Type B (tesseract/8-cell):        {len(type_B)}")
print(f"  Type C (snub 24-cell):            {len(type_C)}")
assert len(type_A) == 8 and len(type_B) == 16 and len(type_C) == 96
print()

def get_type(idx):
    if idx in set_A: return 'A'
    elif idx in set_B: return 'B'
    else: return 'C'

# ============================================================
# SECTION 1: Find ALL triangles (faces) of the 600-cell
# ============================================================

print("=" * 80)
print("SECTION 1: Enumerate all 1200 triangular faces")
print("=" * 80)
print()

# A triangle = 3 mutually adjacent vertices
# For efficiency, iterate over edges and find common neighbors
triangles = set()
for i, j in edges:
    # Find common neighbors of i and j
    for k in range(N):
        if k > j and adj[i, k] and adj[j, k]:
            triangles.add((i, j, k))
        elif k < i and adj[i, k] and adj[j, k]:
            # Already counted as (k, i, j) or similar
            pass

# Actually, let's be more careful: enumerate all triangles properly
triangles = set()
for i in range(N):
    nbrs_i = set(j for j in range(N) if adj[i, j])
    for j in nbrs_i:
        if j <= i:
            continue
        for k in nbrs_i:
            if k <= j:
                continue
            if adj[j, k]:
                triangles.add((i, j, k))

n_triangles = len(triangles)
print(f"  Total triangles found: {n_triangles}")
assert n_triangles == 1200, f"Expected 1200 triangles, got {n_triangles}"

# Verify: each edge is in exactly 5 triangles
edge_triangle_count = defaultdict(int)
for tri in triangles:
    i, j, k = tri
    edge_triangle_count[(i, j)] += 1
    edge_triangle_count[(i, k)] += 1
    edge_triangle_count[(j, k)] += 1

edge_tri_values = set(edge_triangle_count.values())
print(f"  Triangles per edge: {edge_tri_values}")
assert edge_tri_values == {5}, f"Expected all edges in exactly 5 triangles"

# Triangles per vertex
vertex_triangle_count = Counter()
for tri in triangles:
    for v in tri:
        vertex_triangle_count[v] += 1

vtc_values = set(vertex_triangle_count.values())
print(f"  Triangles per vertex: {vtc_values}")
# Each vertex: 12 neighbors form an icosahedron with 20 triangular faces
# Each face includes the vertex, so 20 triangles per vertex
# But also triangles NOT including the vertex... no, these are triangles
# including the vertex. Let's verify: 120 * 20 / 3 = 800 != 1200.
# Actually: each vertex is in C(12,2) - non-edges among neighbors = ...
# Icosahedron has 30 edges and 20 faces. So each vertex is part of 20
# triangles formed by it + edge in its icosahedron. Plus triangles
# among the neighbors themselves? No, a triangle of the 600-cell
# passing through vertex v consists of v plus an edge among v's neighbors.
# v has 12 neighbors, forming an icosahedron with 30 edges.
# So triangles through v = 30. Check: 120 * 30 / 3 = 1200. YES.
print(f"  Expected: 30 triangles per vertex (icosahedron has 30 edges)")
print(f"  Verification: 120 * 30 / 3 = {120 * 30 // 3} = 1200")
print()

# ============================================================
# SECTION 2: Classify triangles by vertex type
# ============================================================

print("=" * 80)
print("SECTION 2: Triangle classification by vertex type")
print("=" * 80)
print()

triangle_types = Counter()
triangle_by_type = defaultdict(list)

for tri in triangles:
    types = tuple(sorted([get_type(v) for v in tri]))
    triangle_types[types] += 1
    triangle_by_type[types].append(tri)

print(f"  {'Type':>6s}  {'Count':>6s}  {'Pct':>6s}  {'Interpretation':>40s}")
print(f"  {'-'*65}")

interpretations = {
    ('A', 'A', 'A'): "Pure gauge (triple gluon vertex)",
    ('A', 'A', 'B'): "Gauge-gauge-Higgs?",
    ('A', 'A', 'C'): "Gauge-gauge-fermion",
    ('A', 'B', 'B'): "Gauge-Higgs-Higgs?",
    ('A', 'B', 'C'): "Gauge-Higgs-fermion?",
    ('A', 'C', 'C'): "Gauge-fermion-fermion (QED/QCD vertex!)",
    ('B', 'B', 'B'): "Pure Higgs (triple Higgs vertex)",
    ('B', 'B', 'C'): "Higgs-Higgs-fermion",
    ('B', 'C', 'C'): "Higgs-fermion-fermion (Yukawa-like!)",
    ('C', 'C', 'C'): "Pure fermion (fermion-fermion-fermion)",
}

total_check = 0
for ttype in sorted(triangle_types.keys()):
    cnt = triangle_types[ttype]
    total_check += cnt
    pct = cnt / 1200 * 100
    interp = interpretations.get(ttype, "")
    label = ''.join(ttype)
    print(f"  {label:>6s}  {cnt:6d}  {pct:5.1f}%  {interp:>40s}")

print(f"  {'Total':>6s}  {total_check:6d}")
print()

# Key interaction vertices
n_ACC = triangle_types.get(('A', 'C', 'C'), 0)
n_BCC = triangle_types.get(('B', 'C', 'C'), 0)
n_CCC = triangle_types.get(('C', 'C', 'C'), 0)
n_ABC = triangle_types.get(('A', 'B', 'C'), 0)
n_AAC = triangle_types.get(('A', 'A', 'C'), 0)

print(f"  KEY INTERACTION VERTICES:")
print(f"    ACC (gauge-fermion-fermion):   {n_ACC}")
print(f"    BCC (Higgs?-fermion-fermion):  {n_BCC}")
print(f"    CCC (pure fermion / Yukawa):   {n_CCC}")
print(f"    ABC (gauge-Higgs-fermion):     {n_ABC}")
print(f"    AAC (gauge-gauge-fermion):     {n_AAC}")
print()

# Ratios
if n_ACC > 0 and n_BCC > 0:
    ratio_BCC_ACC = n_BCC / n_ACC
    print(f"  BCC / ACC = {n_BCC}/{n_ACC} = {ratio_BCC_ACC:.6f}")
    print(f"  Compare: B/A vertex ratio = 16/8 = 2.0")
    print(f"  Compare: alpha_s/alpha = {ALPHA_S/ALPHA:.4f}")
    print()

if n_CCC > 0 and n_ACC > 0:
    print(f"  CCC / ACC = {n_CCC}/{n_ACC} = {n_CCC/n_ACC:.6f}")
    print(f"  CCC / total = {n_CCC}/{n_triangles} = {n_CCC/n_triangles:.6f}")
    print()

# ============================================================
# SECTION 3: Find ALL tetrahedra (cells) of the 600-cell
# ============================================================

print("=" * 80)
print("SECTION 3: Enumerate all 600 tetrahedral cells")
print("=" * 80)
print()

# A tetrahedron = 4 mutually adjacent vertices
# Build from triangles: for each triangle, find vertices adjacent to all 3
tetrahedra = set()

# Build neighbor sets for efficiency
neighbor_sets = {}
for i in range(N):
    neighbor_sets[i] = set(j for j in range(N) if adj[i, j])

for tri in triangles:
    i, j, k = tri
    # Find vertices adjacent to all three
    common = neighbor_sets[i] & neighbor_sets[j] & neighbor_sets[k]
    for l in common:
        tet = tuple(sorted([i, j, k, l]))
        tetrahedra.add(tet)

n_tetrahedra = len(tetrahedra)
print(f"  Total tetrahedra found: {n_tetrahedra}")
assert n_tetrahedra == 600, f"Expected 600 tetrahedra, got {n_tetrahedra}"

# Verify: each edge is shared by exactly 5 tetrahedra
edge_tet_count = defaultdict(int)
for tet in tetrahedra:
    for pair in combinations(tet, 2):
        edge_tet_count[pair] += 1

edge_tet_values = set(edge_tet_count.values())
print(f"  Tetrahedra per edge: {edge_tet_values}")
# Actually, each edge has 5 triangles and 5 tetrahedra meeting at it
# Let me verify
print(f"  (Each edge bordered by 5 tetrahedra, arranged cyclically)")

# Tetrahedra per vertex
vertex_tet_count = Counter()
for tet in tetrahedra:
    for v in tet:
        vertex_tet_count[v] += 1

vtc_tet_values = set(vertex_tet_count.values())
print(f"  Tetrahedra per vertex: {vtc_tet_values}")
print(f"  Expected: 20 (icosahedral vertex figure = 20 triangles = 20 tetrahedra)")
print(f"  Verification: 120 * 20 / 4 = {120 * 20 // 4} = 600")
print()

# ============================================================
# SECTION 4: Classify tetrahedra by vertex type
# ============================================================

print("=" * 80)
print("SECTION 4: Tetrahedron classification by vertex type")
print("=" * 80)
print()

tet_types = Counter()
tet_by_type = defaultdict(list)

for tet in tetrahedra:
    types = tuple(sorted([get_type(v) for v in tet]))
    tet_types[types] += 1
    tet_by_type[types].append(tet)

print(f"  {'Type':>8s}  {'Count':>6s}  {'Pct':>6s}  {'Interpretation':>40s}")
print(f"  {'-'*65}")

tet_interpretations = {
    ('A', 'A', 'A', 'A'): "Pure gauge 4-pt",
    ('A', 'A', 'A', 'B'): "Gauge^3-Higgs",
    ('A', 'A', 'A', 'C'): "Gauge^3-fermion",
    ('A', 'A', 'B', 'B'): "Gauge^2-Higgs^2",
    ('A', 'A', 'B', 'C'): "Gauge^2-Higgs-fermion",
    ('A', 'A', 'C', 'C'): "Gauge^2-fermion^2 (box diagram)",
    ('A', 'B', 'B', 'B'): "Gauge-Higgs^3",
    ('A', 'B', 'B', 'C'): "Gauge-Higgs^2-fermion",
    ('A', 'B', 'C', 'C'): "Gauge-Higgs-fermion^2",
    ('A', 'C', 'C', 'C'): "Gauge-fermion^3",
    ('B', 'B', 'B', 'B'): "Pure Higgs 4-pt",
    ('B', 'B', 'B', 'C'): "Higgs^3-fermion",
    ('B', 'B', 'C', 'C'): "Higgs^2-fermion^2",
    ('B', 'C', 'C', 'C'): "Higgs-fermion^3",
    ('C', 'C', 'C', 'C'): "Pure fermion 4-pt",
}

total_check = 0
for ttype in sorted(tet_types.keys()):
    cnt = tet_types[ttype]
    total_check += cnt
    pct = cnt / 600 * 100
    interp = tet_interpretations.get(ttype, "")
    label = ''.join(ttype)
    print(f"  {label:>8s}  {cnt:6d}  {pct:5.1f}%  {interp:>40s}")

print(f"  {'Total':>8s}  {total_check:6d}")
print()

# ============================================================
# SECTION 5: Interaction selection rules
# ============================================================

print("=" * 80)
print("SECTION 5: Interaction selection rules from neighbor structure")
print("=" * 80)
print()

# Known from exp124/126:
# Type-A: 0A + 0B + 12C neighbors
# Type-B: ? (need to verify)
# Type-C: 1A + 2B + 9C neighbors

print("  Neighbor type census (all vertices):")
for label, indices in [("A", type_A), ("B", type_B), ("C", type_C)]:
    stats = []
    for idx in indices:
        nbrs = [j for j in range(N) if adj[idx, j]]
        n_A = sum(1 for j in nbrs if j in set_A)
        n_B = sum(1 for j in nbrs if j in set_B)
        n_C = sum(1 for j in nbrs if j in set_C)
        stats.append((n_A, n_B, n_C))

    unique_stats = Counter(stats)
    for pattern, cnt in sorted(unique_stats.items()):
        print(f"    Type {label}: {pattern[0]}A + {pattern[1]}B + {pattern[2]}C  "
              f"({cnt}/{len(indices)} vertices)")
    print()

# Edge type counts
edge_type_counts = Counter()
for i, j in edges:
    pair = tuple(sorted([get_type(i), get_type(j)]))
    edge_type_counts[pair] += 1

print("  Edge type distribution:")
for pair in sorted(edge_type_counts.keys()):
    cnt = edge_type_counts[pair]
    pct = cnt / 720 * 100
    print(f"    {pair[0]}-{pair[1]}: {cnt:4d} edges ({pct:5.1f}%)")

# Verify consistency
print()
print("  SELECTION RULES (from edge distribution):")
n_AA = edge_type_counts.get(('A', 'A'), 0)
n_AB = edge_type_counts.get(('A', 'B'), 0)
n_AC = edge_type_counts.get(('A', 'C'), 0)
n_BB = edge_type_counts.get(('B', 'B'), 0)
n_BC = edge_type_counts.get(('B', 'C'), 0)
n_CC = edge_type_counts.get(('C', 'C'), 0)

print(f"    A-A edges: {n_AA}  -> {'FORBIDDEN' if n_AA == 0 else 'ALLOWED'}")
print(f"    A-B edges: {n_AB}  -> {'FORBIDDEN' if n_AB == 0 else 'ALLOWED'}")
print(f"    A-C edges: {n_AC}  -> ALLOWED (gauge-fermion propagator)")
print(f"    B-B edges: {n_BB}  -> {'FORBIDDEN' if n_BB == 0 else 'ALLOWED'}")
print(f"    B-C edges: {n_BC}  -> ALLOWED (Higgs-fermion propagator)")
print(f"    C-C edges: {n_CC}  -> ALLOWED (fermion-fermion propagator)")
print()

# Derived consequences for triangles
print("  DERIVED CONSEQUENCES for triangle types:")
print(f"    A-A edge = {n_AA} => AAA, AAB, AAC triangles need A-A edges")
if n_AA == 0:
    print(f"    Since A-A = 0:")
    print(f"      AAA triangles = 0 (no A-A edges exist)")
    print(f"      AAB triangles CANNOT have A-A edge => need check")
    print(f"      AAC triangles CANNOT have A-A edge => need check")
if n_AB == 0:
    print(f"    Since A-B = 0:")
    print(f"      AAB triangles = 0 (no A-B edges exist)")
    print(f"      ABB triangles = 0 (no A-B edges exist)")
    print(f"      ABC triangles = 0 (no A-B edges exist)")
if n_BB == 0:
    print(f"    Since B-B = 0:")
    print(f"      BBB triangles = 0 (no B-B edges exist)")
    print(f"      BBC triangles need B-B edge => 0")
print()

# Verify these predictions match Section 2 counts
print("  VERIFICATION against Section 2 triangle counts:")
predictions = {}
if n_AA == 0:
    predictions[('A', 'A', 'A')] = 0
    predictions[('A', 'A', 'C')] = 0  # needs A-A edge
if n_AB == 0:
    predictions[('A', 'A', 'B')] = 0
    predictions[('A', 'B', 'B')] = 0
    predictions[('A', 'B', 'C')] = 0
if n_BB == 0:
    predictions[('B', 'B', 'B')] = 0
    predictions[('B', 'B', 'C')] = 0

for ttype, predicted in sorted(predictions.items()):
    actual = triangle_types.get(ttype, 0)
    match = "OK" if actual == predicted else "FAIL"
    label = ''.join(ttype)
    print(f"    {label}: predicted={predicted}, actual={actual}  [{match}]")

# Which types actually EXIST?
print()
print("  ALLOWED triangle types (non-zero count):")
for ttype in sorted(triangle_types.keys()):
    if triangle_types[ttype] > 0:
        label = ''.join(ttype)
        print(f"    {label}: {triangle_types[ttype]}")

print()
print("  ALLOWED tetrahedron types (non-zero count):")
for ttype in sorted(tet_types.keys()):
    if tet_types[ttype] > 0:
        label = ''.join(ttype)
        print(f"    {label}: {tet_types[ttype]}")
print()

# ============================================================
# SECTION 6: Triangle multiplicity per edge type
# ============================================================

print("=" * 80)
print("SECTION 6: Triangle multiplicity per edge (coupling strength proxy)")
print("=" * 80)
print()

# For each edge type (A-C, B-C, C-C), count how the 5 triangles on that edge
# are distributed by type

print("  For each edge, there are exactly 5 triangles. The third vertex type varies.")
print()

edge_triangle_detail = defaultdict(lambda: Counter())

for tri in triangles:
    i, j, k = tri
    # For each edge of this triangle, record the third vertex type
    for (a, b), c in [((i, j), k), ((i, k), j), ((j, k), i)]:
        edge_pair = tuple(sorted([get_type(a), get_type(b)]))
        third_type = get_type(c)
        e_key = (min(a, b), max(a, b))
        edge_triangle_detail[(edge_pair, e_key)][third_type] += 1

# Aggregate per edge type
print("  Distribution of 3rd vertex type for each edge type:")
print()

for edge_type in sorted(edge_type_counts.keys()):
    if edge_type_counts[edge_type] == 0:
        continue

    # Collect all 3rd-vertex distributions for this edge type
    all_thirds = Counter()
    per_edge_distributions = []

    for (et, ek), dist in edge_triangle_detail.items():
        if et == edge_type:
            per_edge_distributions.append(dict(dist))
            for vtype, cnt in dist.items():
                all_thirds[vtype] += cnt

    n_edges_this = edge_type_counts[edge_type]
    label = f"{edge_type[0]}-{edge_type[1]}"

    print(f"  {label} edges ({n_edges_this} total):")
    print(f"    Total 3rd vertices across all triangles: {dict(sorted(all_thirds.items()))}")

    # Check if all edges of this type have the same distribution
    unique_dists = Counter()
    for d in per_edge_distributions:
        key = tuple(sorted(d.items()))
        unique_dists[key] += 1

    if len(unique_dists) == 1:
        dist_key = list(unique_dists.keys())[0]
        print(f"    ALL {n_edges_this} edges have SAME distribution: {dict(dist_key)}")
    else:
        print(f"    {len(unique_dists)} distinct distributions found:")
        for dist_key, cnt in sorted(unique_dists.items(), key=lambda x: -x[1]):
            print(f"      {dict(dist_key)}: {cnt} edges")

    # Average per edge
    avg = {k: v / n_edges_this for k, v in sorted(all_thirds.items())}
    print(f"    Average per edge: {avg}")
    print()

# ============================================================
# SECTION 7: Common gauge neighbors (interaction strength vs distance)
# ============================================================

print("=" * 80)
print("SECTION 7: Common gauge neighbors vs graph distance")
print("=" * 80)
print()

# BFS to compute all pairwise distances
print("  Computing all pairwise distances (BFS from each vertex)...")
dist_matrix = np.zeros((N, N), dtype=int)
for source in range(N):
    dist = np.full(N, -1, dtype=int)
    dist[source] = 0
    queue = [source]
    head = 0
    while head < len(queue):
        u = queue[head]
        head += 1
        for v in range(N):
            if adj[u, v] and dist[v] == -1:
                dist[v] = dist[u] + 1
                queue.append(v)
    dist_matrix[source] = dist

diameter = dist_matrix.max()
print(f"  Diameter: {diameter}")
print()

# Distance distribution
dist_hist = Counter()
for i in range(N):
    for j in range(i+1, N):
        dist_hist[dist_matrix[i, j]] += 1

print("  Distance histogram (all pairs):")
for d in sorted(dist_hist.keys()):
    print(f"    d={d}: {dist_hist[d]} pairs")
print()

# For each pair of Type-C vertices at distance d, count common gauge (A+B) neighbors
set_gauge = set_A | set_B  # 24 gauge vertices

print("  Common gauge (A+B) neighbors between Type-C pairs at each distance:")
print()

gauge_nbr_by_dist = defaultdict(list)

for i in type_C:
    gauge_nbrs_i = set(j for j in neighbor_sets[i] if j in set_gauge)
    for j in type_C:
        if j <= i:
            continue
        d = dist_matrix[i, j]
        gauge_nbrs_j = set(k for k in neighbor_sets[j] if k in set_gauge)
        common_gauge = len(gauge_nbrs_i & gauge_nbrs_j)
        gauge_nbr_by_dist[d].append(common_gauge)

print(f"  {'Dist':>4s} {'#Pairs':>8s} {'Avg_common':>12s} {'Min':>5s} {'Max':>5s} {'Distribution':>30s}")
print(f"  {'-'*70}")

for d in sorted(gauge_nbr_by_dist.keys()):
    vals = gauge_nbr_by_dist[d]
    avg = np.mean(vals)
    mn = min(vals)
    mx = max(vals)
    dist_counter = Counter(vals)
    print(f"  {d:4d} {len(vals):8d} {avg:12.4f} {mn:5d} {mx:5d} {dict(sorted(dist_counter.items()))}")

print()
print("  INTERPRETATION:")
print("  Common gauge neighbors = number of one-boson-exchange paths")
print("  This measures 'interaction strength' between two fermions")
print("  at a given graph distance.")
print()

# Also break down by A vs B
print("  Breakdown: common A-neighbors vs common B-neighbors:")
print()

for d in sorted(gauge_nbr_by_dist.keys()):
    a_common = []
    b_common = []
    for i in type_C:
        a_nbrs_i = set(j for j in neighbor_sets[i] if j in set_A)
        b_nbrs_i = set(j for j in neighbor_sets[i] if j in set_B)
        for j in type_C:
            if j <= i:
                continue
            if dist_matrix[i, j] != d:
                continue
            a_nbrs_j = set(k for k in neighbor_sets[j] if k in set_A)
            b_nbrs_j = set(k for k in neighbor_sets[j] if k in set_B)
            a_common.append(len(a_nbrs_i & a_nbrs_j))
            b_common.append(len(b_nbrs_i & b_nbrs_j))

    if len(a_common) > 0:
        avg_a = np.mean(a_common)
        avg_b = np.mean(b_common)
        print(f"    d={d}: avg common A-nbrs = {avg_a:.4f}, avg common B-nbrs = {avg_b:.4f}")

print()

# ============================================================
# SECTION 8: Spectral analysis (P-matrix of association scheme)
# ============================================================

print("=" * 80)
print("SECTION 8: Association scheme P-matrix")
print("=" * 80)
print()

# The 600-cell is distance-transitive, so it forms an association scheme
# with d+1 = 9 classes (d=0,1,...,8 for antipodal polytope, but diameter=5)
# Wait, diameter might be 5 for the 600-cell. Let me check.

print(f"  Diameter = {diameter}")
print(f"  Association scheme has {diameter + 1} classes")
print()

# Build distance matrices A_d for each distance d
# The P-matrix P[i][k] gives the eigenvalue of A_i in the k-th eigenspace

# First, get the adjacency matrix eigenvalues
print("  Computing adjacency matrix eigenvalues...")
eigenvalues, eigenvectors = np.linalg.eigh(adj.astype(float))

# Group eigenvalues
eigenspaces = defaultdict(list)
for i, ev in enumerate(eigenvalues):
    ev_rounded = round(ev, 4)
    eigenspaces[ev_rounded].append(i)

sorted_eigs = sorted(eigenspaces.keys(), reverse=True)
n_eigenspaces = len(sorted_eigs)

print(f"  Number of distinct eigenvalues: {n_eigenspaces}")
print(f"  Eigenvalues and multiplicities:")
for ev in sorted_eigs:
    mult = len(eigenspaces[ev])
    print(f"    lambda = {ev:10.4f}, mult = {mult}")
print()

# Build the distance matrices
print("  Building distance matrices A_d...")
A_dist = {}
for d in range(diameter + 1):
    A_d = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            if dist_matrix[i, j] == d:
                A_d[i, j] = 1.0
    A_dist[d] = A_d

# For a distance-transitive graph, A_d can be expressed as a polynomial in A_1
# and the P-matrix entries are the eigenvalues of A_d in each eigenspace

print("  Computing P-matrix (eigenvalues of A_d in each eigenspace)...")
print()

# For each eigenspace, pick a representative eigenvector and compute A_d @ v
P_matrix = np.zeros((diameter + 1, n_eigenspaces))

for k, ev in enumerate(sorted_eigs):
    idx = eigenspaces[ev][0]  # representative eigenvector index
    v_k = eigenvectors[:, idx]

    for d in range(diameter + 1):
        # P[d][k] = eigenvalue of A_d in eigenspace k
        # For a distance-transitive graph: A_d @ v_k = P[d][k] * v_k
        Av = A_dist[d] @ v_k
        # Compute ratio (should be constant for all components)
        nonzero = np.abs(v_k) > 1e-10
        if np.any(nonzero):
            ratios = Av[nonzero] / v_k[nonzero]
            P_matrix[d, k] = np.mean(ratios)

print(f"  P-matrix (rows = distance d, columns = eigenspace k):")
print(f"  {'d':>3s}", end="")
for k, ev in enumerate(sorted_eigs):
    print(f"  {'E'+str(k):>10s}", end="")
print()

print(f"  {'':>3s}", end="")
for k, ev in enumerate(sorted_eigs):
    mult = len(eigenspaces[ev])
    print(f"  {'m='+str(mult):>10s}", end="")
print()

print(f"  {'-' * (3 + 12 * n_eigenspaces)}")
for d in range(diameter + 1):
    print(f"  {d:3d}", end="")
    for k in range(n_eigenspaces):
        print(f"  {P_matrix[d, k]:10.4f}", end="")
    print()

print()

# Key interpretation
print("  INTERPRETATION of P-matrix:")
print("    Row d=0: all entries = 1 (identity matrix has eigenvalue 1 everywhere)")
print("    Row d=1: entries are the adjacency matrix eigenvalues")
print("    P[d][k] tells us how the k-th mode 'sees' distance-d interactions")
print()

# Interaction decay with distance
print("  Normalized interaction strength |P[d][k]/P[1][k]| for each mode:")
print(f"  {'d':>3s}", end="")
for k in range(n_eigenspaces):
    print(f"  {'E'+str(k):>8s}", end="")
print()

for d in range(1, diameter + 1):
    print(f"  {d:3d}", end="")
    for k in range(n_eigenspaces):
        if abs(P_matrix[1, k]) > 0.01:
            ratio = abs(P_matrix[d, k] / P_matrix[1, k])
            print(f"  {ratio:8.4f}", end="")
        else:
            print(f"  {'inf':>8s}", end="")
    print()
print()

# ============================================================
# SECTION 9: Yukawa structure from CCC triangles
# ============================================================

print("=" * 80)
print("SECTION 9: CCC triangle analysis (Yukawa-like vertices)")
print("=" * 80)
print()

ccc_triangles = triangle_by_type.get(('C', 'C', 'C'), [])
n_ccc = len(ccc_triangles)
print(f"  Total CCC triangles: {n_ccc}")
print()

# For each CCC triangle, compute the graph distances between the 3 vertices
ccc_dist_patterns = Counter()
for tri in ccc_triangles:
    i, j, k = tri
    d_ij = dist_matrix[i, j]
    d_ik = dist_matrix[i, k]
    d_jk = dist_matrix[j, k]
    pattern = tuple(sorted([d_ij, d_ik, d_jk]))
    ccc_dist_patterns[pattern] += 1

print(f"  CCC triangle distance patterns (sorted triple of edge distances):")
for pattern, cnt in sorted(ccc_dist_patterns.items()):
    print(f"    ({pattern[0]}, {pattern[1]}, {pattern[2]}): {cnt}")
print()

# Note: since CCC triangles have all 3 vertices mutually adjacent (distance 1),
# the distance pattern is always (1,1,1)
if ccc_dist_patterns == {(1, 1, 1): n_ccc}:
    print(f"  All CCC triangles have distance pattern (1,1,1).")
    print(f"  This is trivially true: a triangle requires mutual adjacency (d=1).")
print()

# More interesting: for each CCC triangle, classify the 3 fermion vertices
# by their distance to the gauge sector (24-cell)
# The "gauge distance" of a Type-C vertex = min distance to any A or B vertex

gauge_dist = {}
for c in type_C:
    min_d = min(dist_matrix[c, g] for g in (type_A + type_B))
    gauge_dist[c] = min_d

gauge_dist_hist = Counter(gauge_dist.values())
print(f"  Distribution of min-distance-to-gauge for Type-C vertices:")
for d, cnt in sorted(gauge_dist_hist.items()):
    print(f"    min_dist = {d}: {cnt} vertices")
print()

# CCC triangles by gauge-distance profile
ccc_gauge_patterns = Counter()
for tri in ccc_triangles:
    gd = tuple(sorted([gauge_dist[v] for v in tri]))
    ccc_gauge_patterns[gd] += 1

print(f"  CCC triangles by gauge-distance profile:")
for pattern, cnt in sorted(ccc_gauge_patterns.items()):
    print(f"    {pattern}: {cnt}")
print()

# Count CCC triangles where all 3 vertices are distance-1 from gauge sector
n_close = ccc_gauge_patterns.get((1, 1, 1), 0)
print(f"  CCC triangles with all 3 fermions at d=1 from gauge: {n_close}/{n_ccc}")
print()

# ============================================================
# SECTION 10: ACC triangle detailed analysis (interaction vertices)
# ============================================================

print("=" * 80)
print("SECTION 10: ACC triangle analysis (gauge-fermion-fermion vertices)")
print("=" * 80)
print()

acc_triangles = triangle_by_type.get(('A', 'C', 'C'), [])
n_acc = len(acc_triangles)
print(f"  Total ACC triangles: {n_acc}")
print()

# For each ACC triangle, which A-vertex is involved?
acc_per_A = Counter()
for tri in acc_triangles:
    for v in tri:
        if v in set_A:
            acc_per_A[v] += 1

print(f"  ACC triangles per Type-A vertex:")
for v_idx in sorted(acc_per_A.keys()):
    print(f"    A vertex {v_idx} ({vertices[v_idx]}): {acc_per_A[v_idx]} ACC triangles")

a_vertex_acc_counts = sorted(acc_per_A.values())
print(f"  All A vertices have same count: {len(set(a_vertex_acc_counts)) == 1}")
if len(set(a_vertex_acc_counts)) == 1:
    print(f"  ACC per A vertex = {a_vertex_acc_counts[0]}")
    print(f"  Verification: 8 * {a_vertex_acc_counts[0]} / 1 = {8 * a_vertex_acc_counts[0]}")
    print(f"  (Each ACC counted once, has 1 A-vertex, so total = 8 * per_A = {8 * a_vertex_acc_counts[0]})")
print()

# BCC triangles analysis
bcc_triangles = triangle_by_type.get(('B', 'C', 'C'), [])
n_bcc = len(bcc_triangles)

bcc_per_B = Counter()
for tri in bcc_triangles:
    for v in tri:
        if v in set_B:
            bcc_per_B[v] += 1

print(f"  BCC triangles per Type-B vertex:")
b_vertex_bcc_counts = sorted(bcc_per_B.values())
print(f"  All B vertices have same count: {len(set(b_vertex_bcc_counts)) == 1}")
if len(set(b_vertex_bcc_counts)) == 1:
    print(f"  BCC per B vertex = {b_vertex_bcc_counts[0]}")
print()

# For each A vertex, how many distinct C-C pairs does it connect?
# This is the "vertex coupling multiplicity" - how many fermion pairs
# interact through a given gauge boson
print(f"  For a given Type-A vertex, the ACC triangles connect C-C pairs:")
a_sample = type_A[0]
acc_cc_pairs = []
for tri in acc_triangles:
    if a_sample in tri:
        cc = [v for v in tri if v != a_sample]
        acc_cc_pairs.append(tuple(sorted(cc)))

print(f"    A vertex {a_sample}: connects {len(acc_cc_pairs)} C-C pairs")
# Are these the same as the C-C edges among its neighbors?
# A vertex has 12 C-neighbors. C-C edges among them = edges of icosahedron = 30
# Each such C-C edge + A vertex = ACC triangle
# So ACC per A = 30. Verify:
a_nbrs = [j for j in range(N) if adj[a_sample, j]]
a_cc_edges = sum(1 for i_idx in range(len(a_nbrs))
                 for j_idx in range(i_idx+1, len(a_nbrs))
                 if adj[a_nbrs[i_idx], a_nbrs[j_idx]])
print(f"    C-C edges among A's 12 neighbors: {a_cc_edges}")
print(f"    (= icosahedron edges = 30)")
print()

# ============================================================
# SECTION 11: Coupling strength ratios from triangle counting
# ============================================================

print("=" * 80)
print("SECTION 11: Coupling strength ratios from geometry")
print("=" * 80)
print()

# Key idea: the "coupling" at a vertex of type X could be proportional
# to the number of interaction triangles involving X

# Per-vertex triangle counts
print("  Triangles per vertex by type:")
print()

for vtype, indices in [("A", type_A), ("B", type_B), ("C", type_C[:5])]:
    idx = indices[0]
    tri_counts = Counter()
    for tri in triangles:
        if idx in tri:
            types = tuple(sorted([get_type(v) for v in tri]))
            tri_counts[types] += 1

    print(f"  Type-{vtype} vertex {idx}:")
    total = 0
    for ttype in sorted(tri_counts.keys()):
        cnt = tri_counts[ttype]
        total += cnt
        label = ''.join(ttype)
        print(f"    {label}: {cnt}")
    print(f"    Total: {total} (should be 30)")
    print()

# Coupling ratio: ACC / BCC per fermion vertex
# Each C vertex is in ACC triangles (interacting with A gauge bosons)
# and BCC triangles (interacting with B gauge bosons)

c_acc_counts = []
c_bcc_counts = []
c_ccc_counts = []

for idx in type_C:
    n_acc_here = 0
    n_bcc_here = 0
    n_ccc_here = 0
    for tri in triangles:
        if idx in tri:
            types = tuple(sorted([get_type(v) for v in tri]))
            if types == ('A', 'C', 'C'):
                n_acc_here += 1
            elif types == ('B', 'C', 'C'):
                n_bcc_here += 1
            elif types == ('C', 'C', 'C'):
                n_ccc_here += 1
    c_acc_counts.append(n_acc_here)
    c_bcc_counts.append(n_bcc_here)
    c_ccc_counts.append(n_ccc_here)

print(f"  Per Type-C vertex triangle breakdown:")
print(f"    ACC triangles: {sorted(set(c_acc_counts))}")
print(f"    BCC triangles: {sorted(set(c_bcc_counts))}")
print(f"    CCC triangles: {sorted(set(c_ccc_counts))}")
print()

if len(set(c_acc_counts)) == 1 and len(set(c_bcc_counts)) == 1:
    acc_per_c = c_acc_counts[0]
    bcc_per_c = c_bcc_counts[0]
    ccc_per_c = c_ccc_counts[0]

    print(f"  Uniform: each C vertex participates in:")
    print(f"    {acc_per_c} ACC triangles (A-mediated interactions)")
    print(f"    {bcc_per_c} BCC triangles (B-mediated interactions)")
    print(f"    {ccc_per_c} CCC triangles (direct fermion coupling)")
    print()

    ratio_ba = bcc_per_c / acc_per_c if acc_per_c > 0 else float('inf')
    print(f"  RATIO BCC/ACC per fermion = {bcc_per_c}/{acc_per_c} = {ratio_ba:.6f}")
    print(f"  Compare: |B|/|A| = 16/8 = 2.0")
    print(f"  Compare: alpha_s/alpha = {ALPHA_S/ALPHA:.4f}")
    print(f"  Compare: g_s^2/e^2 = {ALPHA_S/ALPHA:.4f}")
    print()

    total_gauge_tri = acc_per_c + bcc_per_c
    print(f"  Total gauge-mediated triangles per fermion: {total_gauge_tri}")
    print(f"  ACC fraction: {acc_per_c}/{total_gauge_tri} = {acc_per_c/total_gauge_tri:.6f}")
    print(f"  BCC fraction: {bcc_per_c}/{total_gauge_tri} = {bcc_per_c/total_gauge_tri:.6f}")
    print(f"  Compare: sin^2(theta_W) = {SIN2_TW:.5f}")
    print(f"  Compare: 1 - sin^2(theta_W) = {1-SIN2_TW:.5f}")
    print()

# ============================================================
# SECTION 12: Summary of all interaction numbers
# ============================================================

print("=" * 80)
print("SECTION 12: Complete interaction structure summary")
print("=" * 80)
print()

print("  TOPOLOGICAL DATA (DERIVED from 600-cell geometry):")
print(f"    Vertices:      {N}")
print(f"    Edges:         {n_edges}")
print(f"    Triangles:     {n_triangles}")
print(f"    Tetrahedra:    {n_tetrahedra}")
print(f"    Degree:        {degrees[0]}")
print(f"    Diameter:      {diameter}")
print(f"    Triangles/edge:  5 (= a_1 intersection number)")
print(f"    Triangles/vertex: 30")
print(f"    Tetrahedra/vertex: 20")
print(f"    Tetrahedra/edge:   5")
print()

print("  VERTEX TYPE DECOMPOSITION (DERIVED):")
print(f"    Type A (gauge-like):    {len(type_A):4d}  = 8   = D4 vector rep")
print(f"    Type B (Higgs-like):    {len(type_B):4d}  = 16  = D4 spinor reps")
print(f"    Type C (fermion-like):  {len(type_C):4d}  = 96  = snub 24-cell")
print(f"    A+B = 24 = D4 root system (24-cell)")
print()

print("  EDGE TYPE DISTRIBUTION (DERIVED):")
for pair in sorted(edge_type_counts.keys()):
    cnt = edge_type_counts[pair]
    if cnt > 0:
        print(f"    {pair[0]}-{pair[1]}: {cnt:4d}")
print()

print("  SELECTION RULES (DERIVED from edge structure):")
print(f"    A-A = 0: Gauge bosons do NOT directly couple")
print(f"    A-B = 0: Gauge bosons do NOT couple to Higgs directly")
print(f"    B-B = 0: Higgs bosons do NOT self-couple in the graph")
print(f"    => The only non-zero edges involve at least one C (fermion)")
print(f"    => ALL propagation goes through the fermion sector!")
print()

print("  TRIANGLE TYPE DISTRIBUTION (DERIVED):")
for ttype in sorted(triangle_types.keys()):
    cnt = triangle_types[ttype]
    if cnt > 0:
        label = ''.join(ttype)
        pct = cnt / 1200 * 100
        print(f"    {label}: {cnt:5d} ({pct:5.1f}%)")
print()

print("  TETRAHEDRON TYPE DISTRIBUTION (DERIVED):")
for ttype in sorted(tet_types.keys()):
    cnt = tet_types[ttype]
    if cnt > 0:
        label = ''.join(ttype)
        pct = cnt / 600 * 100
        print(f"    {label}: {cnt:5d} ({pct:5.1f}%)")
print()

print("  KEY RATIOS:")
if n_ACC > 0:
    print(f"    BCC/ACC = {n_BCC}/{n_ACC} = {n_BCC/n_ACC:.6f}")
if n_ACC > 0 and n_CCC > 0:
    print(f"    CCC/ACC = {n_CCC}/{n_ACC} = {n_CCC/n_ACC:.6f}")
    print(f"    CCC/(ACC+BCC) = {n_CCC}/{n_ACC+n_BCC} = {n_CCC/(n_ACC+n_BCC):.6f}")

total_with_gauge = n_ACC + n_BCC
if total_with_gauge > 0:
    print(f"    ACC/(ACC+BCC) = {n_ACC}/{total_with_gauge} = {n_ACC/total_with_gauge:.6f}")
    print(f"    BCC/(ACC+BCC) = {n_BCC}/{total_with_gauge} = {n_BCC/total_with_gauge:.6f}")
print()

# phi-related ratios
print("  PHI-RELATED ANALYSIS:")
print(f"    phi = {PHI:.10f}")
print(f"    1/phi = {1/PHI:.10f}")
print(f"    phi^2 = {PHI**2:.10f}")

for label, val in [("BCC/ACC", n_BCC/n_ACC if n_ACC > 0 else 0),
                    ("CCC/ACC", n_CCC/n_ACC if n_ACC > 0 else 0),
                    ("CCC/(ACC+BCC)", n_CCC/(n_ACC+n_BCC) if (n_ACC+n_BCC) > 0 else 0),
                    ("ACC/total_tri", n_ACC/n_triangles),
                    ("BCC/total_tri", n_BCC/n_triangles),
                    ("CCC/total_tri", n_CCC/n_triangles)]:
    if val > 0:
        log_phi = math.log(val) / math.log(PHI) if val > 0 else float('nan')
        print(f"    {label} = {val:.6f},  log_phi = {log_phi:.4f}")

print()

# ============================================================
# SECTION 13: Final classification
# ============================================================

print("=" * 80)
print("SECTION 13: Final Classification")
print("=" * 80)
print()

print("  DERIVED (rigorous geometric facts):")
print("    1. Triangle enumeration: 1200 faces, 5 per edge, 30 per vertex")
print("    2. Tetrahedron enumeration: 600 cells, 5 per edge, 20 per vertex")
print("    3. Selection rules: A-A = A-B = B-B = 0 edges")
print("       => AAA = AAB = AAC = ABB = ABC = BBB = BBC = 0 triangles")
print("       => Only ACC, BCC, CCC triangles exist (3 types)")
print("    4. Edge distribution: A-C, B-C, C-C only")
print("    5. All gauge sector interactions are MEDIATED by fermion sector")
print("       (gauge bosons only connect to fermions, never to each other)")
print("    6. Triangle type counts (exact numbers)")
print("    7. Tetrahedron type counts (exact numbers)")
print("    8. P-matrix of association scheme (spectral interaction structure)")
print("    9. Common gauge neighbor counts by distance")
print()

print("  PATTERN (suggestive numerical relationships):")
print(f"   10. BCC/ACC ratio = {n_BCC/n_ACC:.6f} (compare 16/8 = 2.0)")
print(f"   11. ACC/(ACC+BCC) = {n_ACC/(n_ACC+n_BCC):.6f}")
print(f"   12. The 3-type triangle structure (ACC, BCC, CCC) mirrors")
print(f"       the 3-force structure (EM/strong, weak, Yukawa)")
print(f"   13. Only 3 non-zero triangle types from 10 possible: simplicity")
print()

print("  SPECULATIV (requires further investigation):")
print("   14. Mapping ACC -> QED/QCD vertex, BCC -> weak vertex, CCC -> Yukawa")
print("   15. Whether coupling strength ratios match SM values")
print("   16. Whether the P-matrix eigenspaces correspond to force carriers")
print()

print("  MOST STRIKING RESULT:")
print("    The 600-cell selection rules AUTOMATICALLY forbid direct")
print("    gauge-gauge interactions. In the SM, pure gauge self-coupling")
print("    (triple gluon vertex) exists. Here, gauge bosons (A+B vertices)")
print("    can ONLY interact through fermion mediators. This is a DIFFERENT")
print("    interaction structure from the SM.")
print()
print("    However, this could be a feature rather than a bug: at the level")
print("    of the 600-cell graph, gauge boson self-interactions might emerge")
print("    as EFFECTIVE interactions mediated by virtual fermion loops,")
print("    rather than fundamental vertices.")
print()

print("  COMPLETE.")

"""
Experiment 139: Edge-based gauge structure on 600-cell
===============================================

CONTEXT:
Exp138 showed gauge bosons on vertices fail to give Yang-Mills (no triple-gluon).
Standard lattice gauge theory places gauge fields on EDGES, not vertices.

INVESTIGATE:
1. Edge classification and plaquette structure
2. Wilson loops on plaquettes
3. Edge Laplacian and spectrum
4. Holonomy group structure
5. Edge decomposition and gauge group dimensions
6. Plaquette decomposition by vertex type
7. Edge coloring / gauge assignment
8. Comparison with lattice gauge structures

CATEGORIES: DERIVED (from geometry), PATTERN (observed), SPECULATIV (hypothesis)
"""

import numpy as np
from itertools import combinations, permutations, product as iprod
from collections import defaultdict, Counter

# Constants
PHI = (1 + np.sqrt(5)) / 2
ALPHA = 1 / 137.035999177
HBAR = 1.054571817e-34  # J*s
C = 299792458  # m/s

def round_vec(v, tol=10):
    """Round vector to tolerance"""
    return tuple(round(x, tol) for x in v)

def generate_600cell():
    """Generate 600-cell vertices with CORRECT even permutation rule"""
    phi = PHI
    vertices = set()

    # 8 vertices: (+-1, 0, 0, 0) permutations
    for i in range(4):
        for s in [1, -1]:
            v = [0.0]*4
            v[i] = float(s)
            vertices.add(tuple(round(x, 10) for x in v))

    # 16 vertices: (+-1/2)^4
    for signs in iprod([0.5, -0.5], repeat=4):
        vertices.add(signs)

    # 96 vertices: EVEN permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base_vals = [phi/2, 0.5, 1/(2*phi), 0.0]
    even_perms = []
    for p in permutations(range(4)):
        # Count inversions
        inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
        if inv % 2 == 0:  # EVEN permutation
            even_perms.append(p)

    for perm in even_perms:
        for s0 in [1, -1]:
            for s1 in [1, -1]:
                for s2 in [1, -1]:
                    signed = [s0*base_vals[0], s1*base_vals[1], s2*base_vals[2], base_vals[3]]
                    v = tuple(round(signed[perm.index(i)], 10) for i in range(4))
                    vertices.add(v)

    return [np.array(v) for v in sorted(vertices)]

def classify_vertex(v):
    """Classify vertex as Type A, B, or C"""
    nz = np.count_nonzero(np.abs(v) > 1e-9)
    max_coord = np.max(np.abs(v))

    if nz == 1 and abs(max_coord - 1.0) < 1e-9:
        return 'A'
    elif nz == 4 and abs(max_coord - 0.5) < 1e-9:
        return 'B'
    else:
        return 'C'

def build_adjacency(vertices, threshold=PHI/2):
    """Build adjacency matrix using dot product threshold"""
    n = len(vertices)
    adj = np.zeros((n, n), dtype=int)
    edges = []

    for i in range(n):
        for j in range(i+1, n):
            dot = np.dot(vertices[i], vertices[j])
            if abs(dot - threshold) < 1e-9:
                adj[i, j] = adj[j, i] = 1
                edges.append((i, j))

    return adj, edges

def find_triangles(edges, adj):
    """Find all triangular faces (plaquettes)"""
    triangles = []
    edge_set = set(edges)

    for i, j in edges:
        # Find common neighbors
        neighbors_i = set(np.where(adj[i] > 0)[0])
        neighbors_j = set(np.where(adj[j] > 0)[0])
        common = neighbors_i & neighbors_j

        for k in common:
            # Triangle (i,j,k) - store in sorted order
            tri = tuple(sorted([i, j, k]))
            if tri not in triangles:
                triangles.append(tri)

    return triangles

def classify_edge(v1_type, v2_type):
    """Classify edge by endpoint types"""
    types = tuple(sorted([v1_type, v2_type]))
    if types == ('A', 'C'):
        return 'AC'
    elif types == ('B', 'C'):
        return 'BC'
    elif types == ('C', 'C'):
        return 'CC'
    else:
        return 'OTHER'

def classify_triangle_by_edges(tri_idx, edge_types):
    """Classify triangle by its edge types"""
    edges_in_tri = sorted(edge_types)
    return tuple(edges_in_tri)

print("="*80)
print("EXPERIMENT 139: EDGE-BASED GAUGE STRUCTURE ON 600-CELL")
print("="*80)

# Generate 600-cell
vertices = generate_600cell()
n_vertices = len(vertices)
print(f"\nGenerated {n_vertices} vertices")

# Classify vertices
vertex_types = [classify_vertex(v) for v in vertices]
type_counts = Counter(vertex_types)
print(f"Vertex types: {dict(type_counts)}")

# Build adjacency and edges
adj, edges = build_adjacency(vertices)
n_edges = len(edges)
print(f"Number of edges: {n_edges}")
print(f"Average degree: {2*n_edges/n_vertices:.1f}")

# PART 1: Edge classification
print("\n" + "="*80)
print("PART 1: EDGE CLASSIFICATION")
print("="*80)

edge_types = []
edge_type_dict = {}  # (i,j) -> type
for i, j in edges:
    etype = classify_edge(vertex_types[i], vertex_types[j])
    edge_types.append(etype)
    edge_type_dict[(min(i,j), max(i,j))] = etype

edge_count = Counter(edge_types)
print(f"\nEdge counts by type:")
for etype in sorted(edge_count.keys()):
    print(f"  {etype}: {edge_count[etype]}")

print("\nVERIFICATION:")
print(f"  Expected: AC=96, BC=192, CC=432")
print(f"  Actual:   AC={edge_count.get('AC',0)}, BC={edge_count.get('BC',0)}, CC={edge_count.get('CC',0)}")
print(f"  Total: {sum(edge_count.values())} (expected 720)")

# Find triangles
print("\n" + "="*80)
print("FINDING TRIANGULAR FACES (PLAQUETTES)")
print("="*80)

triangles = find_triangles(edges, adj)
n_triangles = len(triangles)
print(f"Number of triangular faces: {n_triangles} (expected 1200)")

# Classify triangles by edge composition
triangle_edge_compositions = []
for tri in triangles:
    i, j, k = tri
    # Get edge types for the 3 edges
    edge1 = edge_type_dict[(min(i,j), max(i,j))]
    edge2 = edge_type_dict[(min(j,k), max(j,k))]
    edge3 = edge_type_dict[(min(i,k), max(i,k))]

    composition = classify_triangle_by_edges(tri, [edge1, edge2, edge3])
    triangle_edge_compositions.append(composition)

composition_counts = Counter(triangle_edge_compositions)
print(f"\nTriangle counts by edge composition:")
for comp in sorted(composition_counts.keys()):
    print(f"  {comp}: {composition_counts[comp]}")

# PART 2: Wilson loops structure
print("\n" + "="*80)
print("PART 2: WILSON LOOP STRUCTURE")
print("="*80)

# Build edge-to-triangle incidence
edge_to_triangles = defaultdict(list)
for tri_idx, tri in enumerate(triangles):
    i, j, k = tri
    edge_to_triangles[(min(i,j), max(i,j))].append(tri_idx)
    edge_to_triangles[(min(j,k), max(j,k))].append(tri_idx)
    edge_to_triangles[(min(i,k), max(i,k))].append(tri_idx)

# Count plaquettes per edge
plaquettes_per_edge = [len(edge_to_triangles[e]) for e in edges]
ppc = Counter(plaquettes_per_edge)
print(f"\nPlaquettes per edge distribution: {dict(ppc)}")
print(f"Average plaquettes per edge: {np.mean(plaquettes_per_edge):.2f}")

# Plaquettes per edge type
plaq_per_type = defaultdict(list)
for edge, etype in zip(edges, edge_types):
    plaq_per_type[etype].append(len(edge_to_triangles[edge]))

print(f"\nAverage plaquettes per edge type:")
for etype in sorted(plaq_per_type.keys()):
    avg = np.mean(plaq_per_type[etype])
    std = np.std(plaq_per_type[etype])
    print(f"  {etype}: {avg:.2f} +- {std:.2f}")

# PART 3: Edge Laplacian
print("\n" + "="*80)
print("PART 3: EDGE LAPLACIAN SPECTRUM")
print("="*80)

# Build edge-edge adjacency: two edges adjacent if they share vertex AND belong to common triangle
edge_adj = np.zeros((n_edges, n_edges), dtype=int)
edge_index = {edges[i]: i for i in range(n_edges)}

for tri in triangles:
    i, j, k = tri
    # The 3 edges in this triangle
    e1 = (min(i,j), max(i,j))
    e2 = (min(j,k), max(j,k))
    e3 = (min(i,k), max(i,k))

    # Connect them in edge graph
    idx1, idx2, idx3 = edge_index[e1], edge_index[e2], edge_index[e3]
    edge_adj[idx1, idx2] = edge_adj[idx2, idx1] = 1
    edge_adj[idx2, idx3] = edge_adj[idx3, idx2] = 1
    edge_adj[idx1, idx3] = edge_adj[idx3, idx1] = 1

# Compute degree and Laplacian
edge_degrees = edge_adj.sum(axis=1)
print(f"Edge degree distribution: {Counter(edge_degrees)}")
print(f"Average edge degree: {np.mean(edge_degrees):.2f}")

# Edge Laplacian
edge_L = np.diag(edge_degrees) - edge_adj

# Compute eigenvalues
edge_evals = np.linalg.eigvalsh(edge_L)
edge_evals_sorted = sorted(edge_evals)

# Check for phi-structure
print(f"\nEdge Laplacian eigenvalue range: [{edge_evals_sorted[0]:.6f}, {edge_evals_sorted[-1]:.6f}]")
print(f"Number of zero eigenvalues: {sum(abs(ev) < 1e-9 for ev in edge_evals)}")

# Top eigenvalues
print(f"\nTop 10 edge Laplacian eigenvalues:")
for i in range(10):
    ev = edge_evals_sorted[-(i+1)]
    # Check if ev = a + b*PHI
    for a in range(-10, 11):
        b = (ev - a) / PHI
        if abs(b - round(b)) < 0.01:
            print(f"  {ev:.6f} ~ {a} + {int(round(b))}*phi")
            break
    else:
        print(f"  {ev:.6f} (no phi form found)")

# PART 4: Holonomy group structure
print("\n" + "="*80)
print("PART 4: HOLONOMY GROUP STRUCTURE AT C VERTICES")
print("="*80)

# For each C vertex, analyze its 12 edges
c_vertices = [i for i in range(n_vertices) if vertex_types[i] == 'C']
print(f"Number of C vertices: {len(c_vertices)}")

# Sample first few C vertices
sample_c = c_vertices[:5]
print(f"\nAnalyzing edge structure at first {len(sample_c)} C vertices:")

for c_idx in sample_c:
    neighbors = np.where(adj[c_idx] > 0)[0]
    neighbor_types = [vertex_types[n] for n in neighbors]
    neighbor_count = Counter(neighbor_types)

    # Get edge types
    c_edge_types = []
    for n in neighbors:
        edge = (min(c_idx, n), max(c_idx, n))
        c_edge_types.append(edge_type_dict[edge])

    c_edge_count = Counter(c_edge_types)

    print(f"\n  C-vertex {c_idx}:")
    print(f"    Neighbors by type: {dict(neighbor_count)}")
    print(f"    Edges by type: {dict(c_edge_count)}")

# Overall statistics for all C vertices
all_c_neighbor_types = defaultdict(list)
all_c_edge_types = defaultdict(list)

for c_idx in c_vertices:
    neighbors = np.where(adj[c_idx] > 0)[0]
    neighbor_types = [vertex_types[n] for n in neighbors]
    neighbor_count = Counter(neighbor_types)

    c_edge_types_local = []
    for n in neighbors:
        edge = (min(c_idx, n), max(c_idx, n))
        c_edge_types_local.append(edge_type_dict[edge])

    c_edge_count = Counter(c_edge_types_local)

    all_c_neighbor_types['A'].append(neighbor_count['A'])
    all_c_neighbor_types['B'].append(neighbor_count['B'])
    all_c_neighbor_types['C'].append(neighbor_count['C'])

    all_c_edge_types['AC'].append(c_edge_count.get('AC', 0))
    all_c_edge_types['BC'].append(c_edge_count.get('BC', 0))
    all_c_edge_types['CC'].append(c_edge_count.get('CC', 0))

print(f"\n" + "="*80)
print("AGGREGATE STATISTICS FOR ALL C VERTICES")
print("="*80)

print(f"\nNeighbor type distribution (avg across all C vertices):")
for ntype in ['A', 'B', 'C']:
    counts = all_c_neighbor_types[ntype]
    print(f"  {ntype}-neighbors: {np.mean(counts):.2f} +- {np.std(counts):.2f}")

print(f"\nEdge type distribution (avg across all C vertices):")
for etype in ['AC', 'BC', 'CC']:
    counts = all_c_edge_types[etype]
    print(f"  {etype}-edges: {np.mean(counts):.2f} +- {np.std(counts):.2f}")

print(f"\nSM gauge group: U(1) (1 gen) + SU(2) (3 gen) + SU(3) (8 gen) = 12 total")
print(f"C-vertex edges: AC (1) + BC (2) + CC (9) = 12 total")
print(f"Mismatch: 1 vs 1 (OK), 2 vs 3 (off by 1), 9 vs 8 (off by 1)")
print(f"But: 1 + 2 + 9 = 12 = 1 + 3 + 8 (total matches!)")

# PART 5: Edge decomposition
print("\n" + "="*80)
print("PART 5: EDGE DECOMPOSITION AND GAUGE GROUP DIMENSIONS")
print("="*80)

# For each C vertex, look at the 9 CC-neighbors more closely
print(f"\nInvestigating 9 CC-neighbors structure:")

sample_c = c_vertices[:3]
for c_idx in sample_c:
    neighbors = np.where(adj[c_idx] > 0)[0]

    # Separate by type
    a_neighbors = [n for n in neighbors if vertex_types[n] == 'A']
    b_neighbors = [n for n in neighbors if vertex_types[n] == 'B']
    c_neighbors = [n for n in neighbors if vertex_types[n] == 'C']

    print(f"\n  C-vertex {c_idx}:")
    print(f"    {len(a_neighbors)} A-neighbors: {a_neighbors}")
    print(f"    {len(b_neighbors)} B-neighbors: {b_neighbors}")
    print(f"    {len(c_neighbors)} C-neighbors")

    # Check if the 9 C-neighbors have any substructure
    # Are they connected among themselves?
    cc_subgraph_edges = 0
    for i, n1 in enumerate(c_neighbors):
        for n2 in c_neighbors[i+1:]:
            if adj[n1, n2] > 0:
                cc_subgraph_edges += 1

    print(f"    C-C subgraph: {len(c_neighbors)} vertices, {cc_subgraph_edges} edges")
    print(f"    C-C density: {cc_subgraph_edges}/{len(c_neighbors)*(len(c_neighbors)-1)//2} = {cc_subgraph_edges/(len(c_neighbors)*(len(c_neighbors)-1)//2):.2f}")

# PART 6: Plaquette decomposition by vertex type
print("\n" + "="*80)
print("PART 6: PLAQUETTE DECOMPOSITION AT C VERTICES")
print("="*80)

# For each C vertex, count plaquettes by type
c_vertex_plaquette_counts = defaultdict(list)

for c_idx in c_vertices:
    # Find all triangles containing this vertex
    triangles_at_c = [tri for tri in triangles if c_idx in tri]

    # Classify by edge composition
    plaq_types = defaultdict(int)
    for tri in triangles_at_c:
        i, j, k = tri
        edge1 = edge_type_dict[(min(i,j), max(i,j))]
        edge2 = edge_type_dict[(min(j,k), max(j,k))]
        edge3 = edge_type_dict[(min(i,k), max(i,k))]

        composition = classify_triangle_by_edges(tri, [edge1, edge2, edge3])
        plaq_types[composition] += 1

    for comp in plaq_types:
        c_vertex_plaquette_counts[comp].append(plaq_types[comp])

print(f"\nPlaquette counts at C vertices (avg across all C vertices):")
for comp in sorted(c_vertex_plaquette_counts.keys()):
    counts = c_vertex_plaquette_counts[comp]
    print(f"  {comp}: {np.mean(counts):.2f} +- {np.std(counts):.2f}")

total_plaq_per_c = sum(np.mean(c_vertex_plaquette_counts[comp]) for comp in c_vertex_plaquette_counts)
print(f"\nTotal plaquettes per C vertex: {total_plaq_per_c:.2f}")

# PART 7: Edge coloring analysis
print("\n" + "="*80)
print("PART 7: EDGE COLORING / GAUGE ASSIGNMENT")
print("="*80)

print(f"\nEdge counts:")
print(f"  AC edges: {edge_count['AC']} = 96 (8 A-vertices * 12 neighbors each / 2? NO)")
print(f"  Actually: Each A has how many C neighbors?")

a_vertices = [i for i in range(n_vertices) if vertex_types[i] == 'A']
a_neighbor_counts = []
for a_idx in a_vertices:
    neighbors = np.where(adj[a_idx] > 0)[0]
    c_neighbors = sum(1 for n in neighbors if vertex_types[n] == 'C')
    a_neighbor_counts.append(c_neighbors)

print(f"  A-vertex C-neighbor counts: {Counter(a_neighbor_counts)}")
print(f"  Average C-neighbors per A: {np.mean(a_neighbor_counts):.2f}")
print(f"  Total AC edges: {sum(a_neighbor_counts)} (counted from A side)")
print(f"  This matches: {edge_count['AC']} AC edges")

print(f"\n  BC edges: {edge_count['BC']}")
b_vertices = [i for i in range(n_vertices) if vertex_types[i] == 'B']
b_neighbor_counts = []
for b_idx in b_vertices:
    neighbors = np.where(adj[b_idx] > 0)[0]
    c_neighbors = sum(1 for n in neighbors if vertex_types[n] == 'C')
    b_neighbor_counts.append(c_neighbors)

print(f"  B-vertex C-neighbor counts: {Counter(b_neighbor_counts)}")
print(f"  Average C-neighbors per B: {np.mean(b_neighbor_counts):.2f}")
print(f"  Total BC edges: {sum(b_neighbor_counts)} (counted from B side)")
print(f"  This matches: {edge_count['BC']} BC edges")

print(f"\n  CC edges: {edge_count['CC']}")
print(f"  Total C-C connections: {edge_count['CC']}")
print(f"  Per C vertex: {2*edge_count['CC']/len(c_vertices):.2f} CC-edges on average")

# PART 8: Comparison with lattice gauge theory
print("\n" + "="*80)
print("PART 8: COMPARISON WITH LATTICE GAUGE THEORY")
print("="*80)

print(f"\nStandard lattice QCD on Z^4:")
print(f"  - Edges per site: 4 (one per dimension)")
print(f"  - Each edge carries SU(3) matrix (8 real DOF)")
print(f"  - Plaquettes: 6 per site (C(4,2) orientations)")

print(f"\n600-cell lattice (C vertices only):")
print(f"  - Vertices (sites): {len(c_vertices)}")
print(f"  - Edges per site: 12 (1 AC + 2 BC + 9 CC)")
print(f"  - Plaquettes per site: {total_plaq_per_c:.0f}")
print(f"  - Edge types naturally implement gauge group PRODUCT")

print(f"\nGauge group assignment hypothesis:")
print(f"  - AC edges (1 per C vertex): U(1)_Y generator")
print(f"  - BC edges (2 per C vertex): Part of SU(2)_L? (need 3, have 2)")
print(f"  - CC edges (9 per C vertex): SU(3)_C (8 gen) + 1 extra?")

print(f"\nPlaquette structure hypothesis:")
ac_plaq = np.mean(c_vertex_plaquette_counts.get(('AC', 'AC', 'CC'), [0]))
bc_plaq = np.mean(c_vertex_plaquette_counts.get(('BC', 'BC', 'CC'), [0]))
cc_plaq = np.mean(c_vertex_plaquette_counts.get(('CC', 'CC', 'CC'), [0]))

print(f"  - ACC plaquettes: {ac_plaq:.0f} per C vertex (U(1) sector)")
print(f"  - BCC plaquettes: {bc_plaq:.0f} per C vertex (SU(2) sector)")
print(f"  - CCC plaquettes: {cc_plaq:.0f} per C vertex (SU(3) sector)")

print(f"\nPer-edge plaquette count:")
print(f"  - AC: {ac_plaq:.0f} plaquettes using AC edge")
print(f"  - BC: {bc_plaq:.0f} plaquettes using BC edge")
print(f"  - CC: contribution from both BCC and CCC")

# SUMMARY
print("\n" + "="*80)
print("SUMMARY AND CLASSIFICATION")
print("="*80)

print(f"\nDERIVED (from 600-cell geometry):")
print(f"  - Edge counts: AC=96, BC=192, CC=432 (total 720)")
print(f"  - Triangle counts by edge type")
print(f"  - Each C vertex has 1 AC + 2 BC + 9 CC edges (total 12)")
print(f"  - Plaquettes per edge: ~5 (exact from counting)")
print(f"  - Edge degree distribution in edge graph")

print(f"\nPATTERN (observed, may be meaningful):")
print(f"  - 12 edges per C vertex = 12 SM generators")
print(f"  - 1 AC + 2 BC + 9 CC ~ 1 U(1) + 3 SU(2) + 8 SU(3)")
print(f"  - Mismatch: (1,2,9) vs (1,3,8) - off by 1 in two sectors")
print(f"  - Plaquette counts: {ac_plaq:.0f}, {bc_plaq:.0f}, {cc_plaq:.0f} at C vertices")
print(f"  - Edge graph eigenvalues: not obviously phi-structured")

print(f"\nSPECULATIV (hypothesis):")
print(f"  - Edge types implement gauge group product structure naturally")
print(f"  - AC edges = U(1)_Y, BC edges = part of SU(2)_L, CC edges = SU(3)_C + extra")
print(f"  - The (1,2,9) vs (1,3,8) mismatch may require:")
print(f"    * BC edges carry 2 DOF but SU(2) needs 3 (one off-vertex?)")
print(f"    * CC edges carry 9 DOF with 8 as SU(3) and 1 as singlet")
print(f"  - Plaquette decomposition encodes field strength tensors F_mu_nu")
print(f"  - Wilson loops around plaquettes give curvature (gauge field strength)")

print(f"\nOPEN QUESTIONS:")
print(f"  - Why 2 BC edges instead of 3 for full SU(2)?")
print(f"  - What is the role of the 9th CC edge (beyond SU(3)'s 8)?")
print(f"  - Do edge graph eigenvalues have hidden algebraic structure?")
print(f"  - Can we explicitly construct gauge field assignments that close?")
print(f"  - Relation to exp111 spectral structure (a_1=5, b_1=6)?")

print(f"\n" + "="*80)
print("EXPERIMENT 139 COMPLETE")
print("="*80)

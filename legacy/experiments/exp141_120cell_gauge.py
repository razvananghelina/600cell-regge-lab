"""
EXP-141: Gauge bosons on 120-cell (dual of 600-cell)
=====================================================
Alternative 3: Place gauge bosons on 120-cell vertices (= 600-cell tetrahedra).
120-cell: V=600, E=1200, degree=4, diameter=15.
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict, deque

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
print("EXP-141: GAUGE BOSONS ON THE 120-CELL (DUAL OF 600-CELL)")
print("="*70)

# Step 1: Build 600-cell
print("\n--- STEP 1: Build 600-cell ---")
verts = generate_600cell()
N = len(verts)
print(f"600-cell vertices: {N}")
assert N == 120

vtypes = [classify_vertex(v) for v in verts]
type_counts = Counter(vtypes)
print(f"Type counts: A={type_counts['A']}, B={type_counts['B']}, C={type_counts['C']}")

coords = np.array(verts)
dist_sq = np.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=2)
edge_len_sq = 1.0 / PHI**2
adj = dist_sq < edge_len_sq + 1e-6
np.fill_diagonal(adj, False)

neighbors = defaultdict(set)
edges_600 = set()
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j]:
            edges_600.add((i,j))
            neighbors[i].add(j)
            neighbors[j].add(i)

print(f"600-cell edges: {len(edges_600)}")
assert len(edges_600) == 720

# Triangles
print("Finding triangles...")
triangles = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            for k in neighbors[i] & neighbors[j]:
                if k > j:
                    triangles.append((i, j, k))
print(f"Triangles: {len(triangles)}")
assert len(triangles) == 1200

# Tetrahedra
print("Finding tetrahedra...")
tetrahedra = []
for tri in triangles:
    i, j, k = tri
    common = neighbors[i] & neighbors[j] & neighbors[k]
    for l in common:
        if l > k:
            tetrahedra.append((i, j, k, l))
print(f"Tetrahedra: {len(tetrahedra)}")
assert len(tetrahedra) == 600

# Step 2: Build 120-cell dual
print("\n--- STEP 2: Build 120-cell dual graph ---")
tri_to_tet = defaultdict(list)
for t_idx, tet in enumerate(tetrahedra):
    for face in combinations(tet, 3):
        key = tuple(sorted(face))
        tri_to_tet[key].append(t_idx)

adj_120 = defaultdict(set)
for face_key, tet_list in tri_to_tet.items():
    if len(tet_list) == 2:
        t1, t2 = tet_list
        adj_120[t1].add(t2)
        adj_120[t2].add(t1)

edges_120 = set()
for t1 in adj_120:
    for t2 in adj_120[t1]:
        if t1 < t2:
            edges_120.add((t1, t2))

degrees_120 = [len(adj_120[i]) for i in range(600)]
V120 = 600
E120 = len(edges_120)
print(f"120-cell: V={V120}, E={E120}, degree={sorted(set(degrees_120))}")
assert E120 == 1200 and set(degrees_120) == {4}
print(">> 120-cell construction VERIFIED")

# Step 3: Classify 120-cell vertices by tet type
print("\n--- STEP 3: Classify by tetrahedron type ---")
tet_types = []
for tet in tetrahedra:
    tt = ''.join(sorted([vtypes[v] for v in tet]))
    tet_types.append(tt)

tet_type_counts = Counter(tet_types)
TYPE_NAMES = sorted(tet_type_counts.keys())
print("Tetrahedron types (= 120-cell vertex types):")
for tt in TYPE_NAMES:
    print(f"  {tt}: {tet_type_counts[tt]}")
print(f"  Total: {sum(tet_type_counts.values())}")

# Step 4: Edge structure by type
print("\n--- STEP 4: Edge structure by type ---")
edge_type_counts = Counter()
for t1, t2 in edges_120:
    key = tuple(sorted([tet_types[t1], tet_types[t2]]))
    edge_type_counts[key] += 1

print("Edge counts:")
for key in sorted(edge_type_counts.keys()):
    print(f"  {key[0]}-{key[1]}: {edge_type_counts[key]}")

print("\nSelection rules (forbidden edges):")
all_pairs = [(a,b) for a in TYPE_NAMES for b in TYPE_NAMES if a <= b]
forbidden = [p for p in all_pairs if edge_type_counts.get(p, 0) == 0]
if forbidden:
    for p in forbidden:
        print(f"  {p[0]}-{p[1]} = 0 (FORBIDDEN)")
else:
    print("  None - all type pairs have edges")

# Step 5: Neighbor type distribution
print("\n--- STEP 5: Neighbor type distribution ---")
for ttype in TYPE_NAMES:
    indices = [i for i in range(600) if tet_types[i] == ttype]
    patterns = []
    for i in indices:
        local = Counter()
        for nb in adj_120[i]:
            local[tet_types[nb]] += 1
        patterns.append(tuple(sorted(local.items())))
    unique_pats = set(patterns)
    print(f"\nType {ttype} ({len(indices)} verts):")
    print(f"  Unique neighbor patterns: {len(unique_pats)}")
    for pat in sorted(unique_pats):
        cnt = patterns.count(pat)
        print(f"    {dict(pat)} x {cnt}")

# Step 6: Triangles on 120-cell
print("\n--- STEP 6: Triangles on 120-cell ---")
triangles_120 = []
for i in range(600):
    for j in adj_120[i]:
        if j > i:
            for k in adj_120[i] & adj_120[j]:
                if k > j:
                    triangles_120.append((i,j,k))

print(f"120-cell triangles: {len(triangles_120)}")
tri_type_counts = Counter()
for tri in triangles_120:
    tt = tuple(sorted([tet_types[v] for v in tri]))
    tri_type_counts[tt] += 1
print("Triangle types:")
for tt, cnt in sorted(tri_type_counts.items()):
    print(f"  {'-'.join(tt)}: {cnt}")

# Step 7: Shared neighbors (self-interaction test)
print("\n--- STEP 7: Shared neighbors ---")
for et in sorted(edge_type_counts.keys()):
    shared_vals = []
    for i, j in edges_120:
        if tuple(sorted([tet_types[i], tet_types[j]])) == et:
            shared_vals.append(len(adj_120[i] & adj_120[j]))
    if shared_vals:
        print(f"  {et[0]}-{et[1]}: shared = {dict(sorted(Counter(shared_vals).items()))}")

# Step 8: 4-cliques
print("\n--- STEP 8: 4-cliques in 120-cell ---")
cliques4 = []
for tri in triangles_120:
    t0, t1, t2 = tri
    common = adj_120[t0] & adj_120[t1] & adj_120[t2]
    for t3 in common:
        if t3 > t2:
            cliques4.append((t0, t1, t2, t3))
print(f"4-cliques: {len(cliques4)}")
if cliques4:
    cl4_types = Counter()
    for cl in cliques4:
        tt = tuple(sorted([tet_types[v] for v in cl]))
        cl4_types[tt] += 1
    for tt, cnt in sorted(cl4_types.items()):
        print(f"  {'-'.join(tt)}: {cnt}")

# Step 9: Type subgraph analysis
print("\n--- STEP 9: Type subgraph analysis ---")
for ttype in TYPE_NAMES:
    indices = [i for i in range(600) if tet_types[i] == ttype]
    idx_set = set(indices)
    int_edges = sum(1 for i in indices for j in adj_120[i] if j in idx_set) // 2
    int_tri = sum(1 for tri in triangles_120 if all(tet_types[v] == ttype for v in tri))
    print(f"\nType {ttype}: {len(indices)} verts, {int_edges} edges, {int_tri} triangles")
    if len(indices) > 0:
        print(f"  Avg internal degree: {2*int_edges/len(indices):.2f}")
    # Connected components
    if int_edges > 0:
        local_adj = defaultdict(set)
        for i in indices:
            for j in adj_120[i]:
                if j in idx_set:
                    local_adj[i].add(j)
        visited = set()
        n_comp = 0
        for start in indices:
            if start in visited:
                continue
            n_comp += 1
            stack = [start]
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                for nb in local_adj[node]:
                    if nb not in visited:
                        stack.append(nb)
        print(f"  Connected components: {n_comp}")

# Step 10: Spectral analysis
print("\n--- STEP 10: Spectral analysis ---")

# Full 120-cell spectrum
A120 = np.zeros((600, 600))
for t1, t2 in edges_120:
    A120[t1, t2] = 1
    A120[t2, t1] = 1

eigvals = sorted(np.linalg.eigvalsh(A120), reverse=True)

def group_eigenvalues(eigs, tol=1e-4):
    groups = []
    i = 0
    while i < len(eigs):
        val = eigs[i]
        mult = 1
        while i + mult < len(eigs) and abs(eigs[i+mult] - val) < tol:
            mult += 1
        groups.append((round(val, 4), mult))
        i += mult
    return groups

unique = group_eigenvalues(eigvals)
print(f"\nFull 120-cell: {len(unique)} distinct eigenvalues")
print("Top 10:")
for val, mult in unique[:10]:
    print(f"  {val:8.4f}  x {mult}")
print("Bottom 5:")
for val, mult in unique[-5:]:
    print(f"  {val:8.4f}  x {mult}")

# Check phi structure
print("\nPhi-based check (a + b*phi):")
for val, mult in unique[:15]:
    best_err = 999
    best_ab = None
    for a in range(-10, 11):
        for b in range(-10, 11):
            err = abs(val - (a + b*PHI))
            if err < best_err:
                best_err = err
                best_ab = (a, b)
    if best_err < 0.01:
        a, b = best_ab
        print(f"  {val:8.4f} = {a} + {b}*phi  (m={mult})")
    else:
        print(f"  {val:8.4f} NOT a+b*phi  (m={mult}, best_err={best_err:.3f})")

# Multiplicities that are perfect squares
print("\nPerfect square multiplicities:")
for val, mult in unique:
    sq = int(round(mult**0.5))
    if sq*sq == mult:
        print(f"  {val:8.4f} m={mult}={sq}^2")

# Type subgraph spectra
for ttype in TYPE_NAMES:
    indices = sorted([i for i in range(600) if tet_types[i] == ttype])
    n_v = len(indices)
    if n_v < 2 or n_v > 400:
        continue
    idx = np.array(indices)
    sub_A = A120[np.ix_(idx, idx)]
    sub_eigs = sorted(np.linalg.eigvalsh(sub_A), reverse=True)
    sub_unique = group_eigenvalues(sub_eigs)
    print(f"\nType {ttype} subgraph ({n_v} verts):")
    print(f"  Top 5: {sub_unique[:5]}")
    print(f"  Bottom 3: {sub_unique[-3:]}")

# Step 11: Diameter
print("\n--- STEP 11: Diameter ---")
def bfs_dist(start):
    dist = [-1]*600
    dist[start] = 0
    q = deque([start])
    while q:
        v = q.popleft()
        for nb in adj_120[v]:
            if dist[nb] == -1:
                dist[nb] = dist[v] + 1
                q.append(nb)
    return dist

d0 = bfs_dist(0)
print(f"Max distance from vertex 0: {max(d0)}")
dist_hist = Counter(d0)
print("Distance distribution:")
for d in sorted(dist_hist.keys()):
    print(f"  d={d:2d}: {dist_hist[d]:4d}")

# Step 12: Connectivity matrix
print("\n--- STEP 12: Inter-type connectivity ---")
print("Edge matrix:")
header = "       " + "  ".join(f"{t:>6s}" for t in TYPE_NAMES)
print(header)
for t1 in TYPE_NAMES:
    row = f"  {t1}: "
    for t2 in TYPE_NAMES:
        key = tuple(sorted([t1, t2]))
        c = edge_type_counts.get(key, 0)
        row += f"  {c:6d}"
    print(row)

print("\nEdge density (actual/possible):")
for t1 in TYPE_NAMES:
    row = f"  {t1}: "
    for t2 in TYPE_NAMES:
        n1 = tet_type_counts[t1]
        n2 = tet_type_counts[t2]
        key = tuple(sorted([t1, t2]))
        c = edge_type_counts.get(key, 0)
        possible = n1*(n1-1)/2 if t1 == t2 else n1*n2
        d = c/possible if possible > 0 else 0
        row += f"  {d:.5f}"
    print(row)

# Summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"\n120-cell: V=600, E=1200, deg=4, diameter={max(d0)}")
print(f"Vertex types: {dict(sorted(tet_type_counts.items()))}")
print(f"Triangles: {len(triangles_120)}")
print(f"4-cliques: {len(cliques4)}")
print(f"\nForbidden edges: {len(forbidden)}")
for p in forbidden:
    print(f"  {p[0]}-{p[1]}")
print(f"\nKey comparisons:")
print(f"  600-cell: deg=12=dim(SM), 3 types (A,B,C), strong selection rules")
print(f"  120-cell: deg=4=rank(SM), 3 types (ACCC,BCCC,CCCC)")

print("\n" + "="*70)
print("END EXP-141")
print("="*70)

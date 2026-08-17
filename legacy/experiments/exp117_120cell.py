"""
EXP-117: 120-Cell Analysis (Dual of 600-Cell)
==============================================

The 120-cell {5,3,3} is the dual of the 600-cell {3,3,5}.
- 600 vertices, 1200 edges, 720 pentagonal faces, 120 dodecahedral cells
- Vertex degree: 4 (vertex figure = tetrahedron {3,3})
- Each vertex corresponds to a tetrahedron of the 600-cell

CONSTRUCTION:
Build 600-cell -> find all 600 tetrahedral cells -> dual graph
(two 120-cell vertices adjacent iff corresponding tetrahedra share a face)

QUESTIONS:
1. Are eigenvalues still a+b*phi?
2. What are the multiplicities? Still n^2?
3. What's the diameter?
4. What are the intersection numbers?
5. Does the dual structure illuminate particle physics?

Author: Claude (exp117)
Date: 2026-02-08
"""

import numpy as np
from itertools import permutations
from collections import Counter, deque

PHI = (1 + np.sqrt(5)) / 2

print("=" * 80)
print("EXP-117: 120-Cell Analysis (Dual of 600-Cell)")
print("=" * 80)
print()

# ============================================================
# SECTION 1: Build the 600-cell
# ============================================================

print("-" * 80)
print("SECTION 1: Building the 600-cell")
print("-" * 80)
print()

def build_600cell():
    """Build 600-cell vertices and adjacency."""
    phi = PHI
    vertices = []

    # 8 vertices: Type A (cross-polytope)
    for i in range(4):
        for s in [1, -1]:
            v = [0, 0, 0, 0]
            v[i] = s
            vertices.append(v)

    # 16 vertices: Type B (tesseract)
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    vertices.append([s0, s1, s2, s3])

    # 96 vertices: Type C (snub 24-cell)
    base_vals = [phi/2, 0.5, 1/(2*phi), 0]
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                signed = [s0*base_vals[0], s1*base_vals[1],
                          s2*base_vals[2], base_vals[3]]
                for p in permutations(range(4)):
                    inv = sum(1 for i in range(4)
                              for j in range(i+1, 4) if p[i] > p[j])
                    if inv % 2 == 0:
                        v = [signed[p[i]] for i in range(4)]
                        v_rounded = [round(x, 10) for x in v]
                        if v_rounded not in vertices:
                            vertices.append(v_rounded)

    N = len(vertices)
    verts = [np.array(v) for v in vertices]

    edge_threshold = 1.0 / PHI + 0.01
    adj = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(i+1, N):
            d = np.linalg.norm(verts[i] - verts[j])
            if d < edge_threshold:
                adj[i, j] = 1
                adj[j, i] = 1

    return vertices, verts, adj

vertices_600, verts_600, adj_600 = build_600cell()
N_600 = len(vertices_600)
E_600 = adj_600.sum() // 2
deg_600 = adj_600.sum(axis=1)

print(f"  600-cell: V={N_600}, E={E_600}, degree={deg_600[0]}")

# ============================================================
# SECTION 2: Find all tetrahedral cells
# ============================================================

print()
print("-" * 80)
print("SECTION 2: Finding all 600 tetrahedral cells")
print("-" * 80)
print()

neighbors = [set(j for j in range(N_600) if adj_600[i, j]) for i in range(N_600)]

tetrahedra = []
for i in range(N_600):
    for j in neighbors[i]:
        if j <= i:
            continue
        common_ij = neighbors[i] & neighbors[j]
        for k in common_ij:
            if k <= j:
                continue
            common_ijk = common_ij & neighbors[k]
            for l in common_ijk:
                if l <= k:
                    continue
                tetrahedra.append((i, j, k, l))

print(f"  Found {len(tetrahedra)} tetrahedra (expected 600)")

# ============================================================
# SECTION 3: Build 120-cell dual graph
# ============================================================

print()
print("-" * 80)
print("SECTION 3: Building 120-cell dual graph")
print("-" * 80)
print()

N_120 = len(tetrahedra)

# Build face-to-tetrahedron mapping
face_to_tet = {}
for idx, tet in enumerate(tetrahedra):
    for skip in range(4):
        face = frozenset(tet[i] for i in range(4) if i != skip)
        if face in face_to_tet:
            face_to_tet[face].append(idx)
        else:
            face_to_tet[face] = [idx]

# Build adjacency for 120-cell
adj_120 = np.zeros((N_120, N_120), dtype=int)
for face, tets in face_to_tet.items():
    if len(tets) == 2:
        i, j = tets
        adj_120[i, j] = 1
        adj_120[j, i] = 1

E_120 = adj_120.sum() // 2
deg_120 = adj_120.sum(axis=1)

n_faces_1 = sum(1 for t in face_to_tet.values() if len(t) == 1)
n_faces_2 = sum(1 for t in face_to_tet.values() if len(t) == 2)

print(f"  120-cell dual graph:")
print(f"    Vertices: {N_120}")
print(f"    Edges: {E_120}")
print(f"    Vertex degree: min={deg_120.min()}, max={deg_120.max()}")
print(f"    Total faces: {len(face_to_tet)} (boundary: {n_faces_1}, shared: {n_faces_2})")
print()
print(f"  Expected: V=600, E=1200, degree=4")
print(f"  Match: {'YES' if N_120 == 600 and E_120 == 1200 and np.all(deg_120 == 4) else 'NO'}")

# ============================================================
# SECTION 4: Distance distribution and diameter
# ============================================================

print()
print("-" * 80)
print("SECTION 4: Distance distribution and diameter")
print("-" * 80)
print()

def bfs_distances(adj, start):
    N = adj.shape[0]
    dist = [-1] * N
    dist[start] = 0
    queue = deque([start])
    while queue:
        v = queue.popleft()
        for u in range(N):
            if adj[v, u] and dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist

dist_0 = bfs_distances(adj_120, 0)
max_dist = max(dist_0)
dist_counts = Counter(dist_0)

print(f"  Distance distribution from vertex 0:")
print(f"  {'d':>4s} {'count':>6s} {'cum':>6s}")
print(f"  {'-'*20}")
cum = 0
for d in range(max_dist + 1):
    c = dist_counts[d]
    cum += c
    print(f"  {d:4d} {c:6d} {cum:6d}")

print(f"\n  Diameter: {max_dist}")

# Verify from a few vertices
diameters = [max(bfs_distances(adj_120, s)) for s in [0, 100, 300, 500]]
print(f"  Diameter from vertices [0,100,300,500]: {diameters}")
print(f"  All same: {len(set(diameters)) == 1}")

# Compare with 600-cell
dist_600_0 = bfs_distances(adj_600, 0)
dist_600_counts = Counter(dist_600_0)

print(f"\n  For reference - 600-cell distances from vertex 0:")
for d in sorted(dist_600_counts.keys()):
    print(f"    d={d}: {dist_600_counts[d]} vertices")

# ============================================================
# SECTION 5: Spectral analysis
# ============================================================

print()
print("-" * 80)
print("SECTION 5: Spectral analysis of 120-cell")
print("-" * 80)
print()

print("  Computing eigenvalues of 600x600 adjacency matrix...")
evals_120 = np.linalg.eigvalsh(adj_120.astype(float))
evals_120 = sorted(evals_120, reverse=True)

# Find unique eigenvalues with multiplicities
unique_evals = []
tol = 0.001
i = 0
while i < len(evals_120):
    val = evals_120[i]
    mult = 1
    while i + mult < len(evals_120) and abs(evals_120[i + mult] - val) < tol:
        mult += 1
    avg_val = sum(evals_120[i:i + mult]) / mult
    unique_evals.append((avg_val, mult))
    i += mult

print(f"  Number of distinct eigenvalues: {len(unique_evals)}")
print(f"  (600-cell has 9 distinct eigenvalues)")
print()

print(f"  {'theta':>12s} {'mult':>6s} {'n=sqrt':>7s} {'n^2?':>5s} {'a+b*phi':>15s}")
print(f"  {'-'*50}")

n_phi_form = 0
n_sq_mult = 0

for val, mult in unique_evals:
    sqrt_m = np.sqrt(mult)
    is_sq = abs(sqrt_m - round(sqrt_m)) < 0.01
    if is_sq:
        n_sq_mult += 1

    # Check if val = a + b*phi for small integers a, b
    best_match = ""
    for b in range(-20, 21):
        a_float = val - b * PHI
        a_int = round(a_float)
        if abs(a_float - a_int) < 0.005:
            check = a_int + b * PHI
            if abs(check - val) < 0.005:
                sign = "+" if b >= 0 else ""
                best_match = f"{a_int}{sign}{b}*phi"
                n_phi_form += 1
                break

    sq_note = f"{int(round(sqrt_m))}^2" if is_sq else ""
    print(f"  {val:12.4f} {mult:6d} {sqrt_m:7.2f} {sq_note:>5s} {best_match:>15s}")

print()
print(f"  Eigenvalues of form a+b*phi: {n_phi_form}/{len(unique_evals)}")
print(f"  Multiplicities that are n^2: {n_sq_mult}/{len(unique_evals)}")

# ============================================================
# SECTION 6: Intersection numbers
# ============================================================

print()
print("-" * 80)
print("SECTION 6: Intersection numbers")
print("-" * 80)
print()

# Build full distance matrix (needed for intersection numbers)
print("  Computing full distance matrix (600x600)...")
all_dists = np.zeros((N_120, N_120), dtype=int)
for i in range(N_120):
    d = bfs_distances(adj_120, i)
    all_dists[i] = d
print("  Done.")
print()

# For each edge (d=1), count how neighbors of j distribute by distance to i
edges_list = [(i, j) for i in range(N_120) for j in range(i + 1, N_120)
              if adj_120[i, j]]

# Check all edges
profile_counter = Counter()
all_profiles = []
for i, j in edges_list:
    profile = {}
    for k in range(N_120):
        if adj_120[j, k] and k != i:
            d_ik = all_dists[i, k]
            profile[d_ik] = profile.get(d_ik, 0) + 1
    p_tuple = tuple(sorted(profile.items()))
    profile_counter[p_tuple] += 1
    all_profiles.append(p_tuple)

print(f"  Checked all {len(edges_list)} edges:")
print(f"  Number of distinct intersection profiles: {len(profile_counter)}")
print()

if len(profile_counter) == 1:
    print("  ALL edges have SAME profile -> association scheme for d=1!")
    profile = dict(list(profile_counter.keys())[0])
    for d in sorted(profile.keys()):
        print(f"    p(1,{d}) = {profile[d]}")

    a1 = profile.get(1, 0)
    b1 = profile.get(2, 0)
    print(f"\n  a_1 = {a1} (common neighbors per edge)")
    print(f"  b_1 = {b1} (next-shell neighbors per edge)")
    print(f"  a_1 + b_1 + 1 = {a1 + b1 + 1} (should = degree = {deg_120[0]})")
else:
    for profile, count in profile_counter.most_common(5):
        print(f"  Profile {dict(profile)}: {count} edges")

# Also check: is it distance-regular?
print()
print("  Checking distance-regularity...")

is_dist_regular = True
for d_check in range(1, min(max_dist + 1, 6)):
    profiles_d = Counter()
    pair_count = 0
    for i in range(min(N_120, 100)):
        for j in range(N_120):
            if all_dists[i, j] == d_check:
                profile = {}
                for k in range(N_120):
                    if adj_120[j, k]:
                        d_ik = all_dists[i, k]
                        profile[d_ik] = profile.get(d_ik, 0) + 1
                profiles_d[tuple(sorted(profile.items()))] += 1
                pair_count += 1

    n_distinct = len(profiles_d)
    if n_distinct == 1:
        p = dict(list(profiles_d.keys())[0])
        print(f"    d={d_check}: REGULAR (profile: {p})")
    else:
        print(f"    d={d_check}: NOT regular ({n_distinct} distinct profiles)")
        is_dist_regular = False

print(f"\n  Distance-regular: {'YES' if is_dist_regular else 'Partial (checked d=1..5)'}")

# ============================================================
# SECTION 7: Comparison table
# ============================================================

print()
print("-" * 80)
print("SECTION 7: 600-cell vs 120-cell comparison")
print("-" * 80)
print()

# 600-cell spectral data (from exp110)
evals_600_data = [
    (12, 1), (6*PHI, 4), (4*PHI, 9), (3, 16), (0, 25),
    (-2, 36), (4 - 4*PHI, 9), (-3, 16), (6 - 6*PHI, 4)
]

print(f"  {'Property':30s} {'600-cell':>15s} {'120-cell':>15s} {'Relation':>15s}")
print(f"  {'-'*78}")
print(f"  {'Schlafli symbol':30s} {'{{3,3,5}}':>15s} {'{{5,3,3}}':>15s} {'dual':>15s}")
print(f"  {'Vertices':30s} {120:>15d} {N_120:>15d} {'cells<->verts':>15s}")
print(f"  {'Edges':30s} {720:>15d} {E_120:>15d} {'faces<->edges':>15s}")
print(f"  {'Faces':30s} {1200:>15d} {'720':>15s} {'edges<->faces':>15s}")
print(f"  {'Cells':30s} {600:>15d} {'120':>15s} {'verts<->cells':>15s}")
print(f"  {'Vertex degree':30s} {12:>15d} {deg_120[0]:>15d} {f'{12}/{deg_120[0]}={12//deg_120[0] if deg_120[0]>0 else 0}':>15s}")
print(f"  {'Diameter':30s} {'5':>15s} {f'{max_dist}':>15s} {f'ratio={max_dist/5:.1f}':>15s}")
print(f"  {'Distinct eigenvalues':30s} {'9':>15s} {f'{len(unique_evals)}':>15s}")

# Eigenvalue comparison
print()
print(f"  600-cell eigenvalues (9 total):")
for val, mult in evals_600_data:
    print(f"    theta = {val:8.4f}, mult = {mult:4d} = {int(round(np.sqrt(mult)))}^2")

print()
print(f"  120-cell eigenvalues ({len(unique_evals)} total):")
for val, mult in unique_evals:
    sqrt_m = np.sqrt(mult)
    sq_str = f" = {int(round(sqrt_m))}^2" if abs(sqrt_m - round(sqrt_m)) < 0.01 else ""
    print(f"    theta = {val:8.4f}, mult = {mult:4d}{sq_str}")

# ============================================================
# SECTION 8: Physical interpretation
# ============================================================

print()
print("-" * 80)
print("SECTION 8: Physical interpretation")
print("-" * 80)
print()

# Duality: vertices <-> cells means particles <-> interactions?
print("  DUALITY INTERPRETATION:")
print()
print("  600-cell: 120 vertices = particles?")
print("  120-cell: 600 vertices = interactions/cells?")
print()

# Check: degree 4 in 120-cell. What does this mean?
print("  120-cell vertex degree = 4:")
print("    Each tetrahedron shares faces with exactly 4 neighbors")
print("    4 = dimension of ambient space R^4")
print("    4 = rank of SM gauge group = rank(SU(3)xSU(2)xU(1))")
print()

# Ratio V(120-cell)/V(600-cell) = 600/120 = 5
print(f"  V(120-cell)/V(600-cell) = {N_120}/{N_600} = {N_120/N_600}")
print(f"  This is 5 = diameter of 600-cell = a_1 intersection number")
print()

# Check if 120-cell eigenvalue ratios give coupling constants
if len(unique_evals) >= 4:
    print("  Eigenvalue ratios (120-cell):")
    theta_max = unique_evals[0][0]
    for i, (val, mult) in enumerate(unique_evals[:8]):
        ratio = theta_max / val if val != 0 else float('inf')
        print(f"    theta_{i}/theta_0 = {val:.4f}/{theta_max:.4f} = {val/theta_max:.6f}")
    print()

# ============================================================
# SECTION 9: Summary
# ============================================================

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

print(f"  120-cell successfully constructed from 600-cell dual")
print(f"  V={N_120}, E={E_120}, degree={deg_120[0]}, diameter={max_dist}")
print()
print(f"  SPECTRAL PROPERTIES:")
print(f"    Distinct eigenvalues: {len(unique_evals)}")
print(f"    All a+b*phi form: {'YES' if n_phi_form == len(unique_evals) else f'NO ({n_phi_form}/{len(unique_evals)})'}")
print(f"    All mult = n^2:   {'YES' if n_sq_mult == len(unique_evals) else f'NO ({n_sq_mult}/{len(unique_evals)})'}")
print()

if len(profile_counter) == 1:
    profile = dict(list(profile_counter.keys())[0])
    a1 = profile.get(1, 0)
    b1 = profile.get(2, 0)
    print(f"  INTERSECTION NUMBERS:")
    print(f"    a_1 = {a1} (common neighbors per edge)")
    print(f"    b_1 = {b1} (next-shell per edge)")
    print(f"    600-cell: a_1=5, b_1=6")
    print()

print(f"  KEY COMPARISONS:")
print(f"    600-cell diameter = 5, 120-cell diameter = {max_dist}")
print(f"    600-cell degree = 12, 120-cell degree = {deg_120[0]}")
print(f"    Vertex ratio: 600/120 = 5")
print()

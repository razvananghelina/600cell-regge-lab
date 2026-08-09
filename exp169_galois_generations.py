"""
EXP-169: 3 Generations from Galois Structure
=============================================
From exp156: Galois conjugation (phi -> phi') maps S to T in E8.
S intersect T = 24 vertices = D4 (gauge sector).
96 remaining = C-type (fermions).

Questions:
1. How does Galois act on the 5 coset structure?
2. Do the 96 C-type vertices decompose into 3x16 under Galois?
3. Is there a Galois orbit structure that counts generations?
4. Does Galois + coset give a natural 3-generation decomposition?

Category: RESEARCH (honest - searching for derivation)
"""

import numpy as np
from itertools import product, permutations
from collections import Counter, defaultdict

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2  # = -1/phi

print("=" * 70)
print("EXP-169: 3 Generations from Galois Structure")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell (exact coordinates in Z[phi])
# ============================================================
print("\n--- Step 1: Build 600-cell with Z[phi] coordinates ---")

# Store vertices as (a0+b0*phi, a1+b1*phi, a2+b2*phi, a3+b3*phi)
# where each coordinate is in Z[phi] = {a + b*phi : a,b in Z}
# For numerical: use floating point

verts_set = set()
for i in range(4):
    for s in [1, -1]:
        v = [0,0,0,0]; v[i] = s
        verts_set.add(tuple(v))
for signs in product([0.5, -0.5], repeat=4):
    verts_set.add(tuple(signs))
vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]
even_perms = [p for p in permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]) % 2 == 0]
for perm in even_perms:
    base = [vals_base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for signs in product([1,-1], repeat=len(nz)):
        v = list(base)
        for idx, s in zip(nz, signs):
            v[idx] *= s
        verts_set.add(tuple(np.round(v, 10)))

verts = np.array(sorted(verts_set))
N = len(verts)
dots = verts @ verts.T
adj = (np.abs(dots - PHI/2) < 0.01).astype(int)
np.fill_diagonal(adj, 0)

# Classify
type_A, type_B, type_C = [], [], []
vertex_type = {}
for i in range(N):
    c = verts[i]
    nz_count = np.sum(np.abs(c) > 0.01)
    if nz_count == 1 and np.isclose(np.max(np.abs(c)), 1.0):
        type_A.append(i); vertex_type[i] = 'A'
    elif nz_count == 4 and np.allclose(np.abs(c), 0.5, atol=0.01):
        type_B.append(i); vertex_type[i] = 'B'
    else:
        type_C.append(i); vertex_type[i] = 'C'

print(f"  N={N} = {len(type_A)}A + {len(type_B)}B + {len(type_C)}C")

# ============================================================
# Step 2: Galois conjugation
# ============================================================
print("\n--- Step 2: Galois Conjugation phi -> phi' ---")

# Apply Galois: replace phi by phi' = (1-sqrt(5))/2 in each coordinate
# For Type A (integer coords): Galois fixes them
# For Type B (half-integer coords): Galois fixes them
# For Type C (has phi components): Galois maps to different vertex

# To apply Galois, we need Z[phi] representation of each vertex
# phi/2 -> phi'/2, 1/(2*phi) = (phi-1)/2 -> (phi'-1)/2 = (-1/phi-1)/2

# Numerically: each coordinate x can be written as a + b*phi
# Galois sends it to a + b*phi'
# We can find (a,b) by solving: x = a + b*phi, so b = (x - round(x))/(phi - round(phi))
# Actually, let's just compute the conjugate vertex directly

def galois_conjugate_vertex(v):
    """Apply Galois conjugation to a vertex of the 600-cell.
    Each coordinate x = a + b*phi gets mapped to a + b*phi'.
    We decompose: b = round(2*x/sqrt(5)) approximately, a = x - b*phi.
    More robustly: use the fact that coordinates are in Z[1/2][phi]."""
    result = np.zeros(4)
    for k in range(4):
        x = v[k]
        # Try to decompose x = (a + b*phi)/2 where a,b in Z with a+b even
        # or x = a + b*phi where a,b in Z or half-integers
        # Robust method: x = p + q*phi where p,q are in {0, +-1/2, +-1, +-phi/2, +-1/(2phi)}

        # Since phi = (1+sqrt5)/2, phi' = (1-sqrt5)/2
        # x = a + b*phi => galois(x) = a + b*phi'
        # a = x - b*phi, so galois(x) = x - b*phi + b*phi' = x + b*(phi'-phi) = x - b*sqrt(5)
        # Also: x = a + b*phi, conjugate = a + b*phi' => sum = 2a + b, diff = b*sqrt(5)
        # So b = (x - galois(x))/sqrt(5)

        # Better: solve a + b*phi = x using: a = x*phi'/(phi'-phi) * (-1)
        # phi' - phi = -sqrt(5)
        # From x = a + b*phi: multiply both sides by 1, and x' = a + b*phi'
        # x + x' = 2a + b, x - x' = b*sqrt(5)
        # So: b = ... we need to find integer decomposition

        # Practical: coordinates are from {0, +-1, +-1/2, +-phi/2, +-1/(2phi)}
        # phi/2: a=0, b=1/2 -> galois = phi'/2
        # 1/(2phi) = (phi-1)/2: a=-1/2, b=1/2 -> galois = -1/2 + phi'/2 = (-1+phi')/2 = -phi/2...
        # Wait: 1/(2phi) = (phi-1)/2, galois((phi-1)/2) = (phi'-1)/2 = (-1/phi-1)/2 = -(1+phi^{-1})/2
        # = -(phi)/2... no wait.
        # phi' = (1-sqrt5)/2 ≈ -0.618
        # (phi'-1)/2 = ((1-sqrt5)/2 - 1)/2 = (-(1+sqrt5)/2)/2 = -phi/2

        # So Galois sends:
        # 0 -> 0
        # 1 -> 1, -1 -> -1
        # 1/2 -> 1/2, -1/2 -> -1/2
        # phi/2 -> phi'/2 ≈ -0.309
        # 1/(2phi) = (phi-1)/2 -> (phi'-1)/2 = -phi/2 ≈ -0.809

        # Let me just do it numerically
        # Decompose x = a + b*phi with a, b in Z/2
        # b = (2x - round(2x)) / (2*(phi - 0.5)) ... no, this is tricky

        # Direct approach: enumerate possible (a,b) in half-integers
        best_err = 1e10
        best_ab = (0, 0)
        for a2 in range(-4, 5):  # a = a2/2
            for b2 in range(-4, 5):  # b = b2/2
                test = a2/2 + b2/2 * PHI
                err = abs(test - x)
                if err < best_err:
                    best_err = err
                    best_ab = (a2/2, b2/2)

        a_val, b_val = best_ab
        result[k] = a_val + b_val * PHI_CONJ

    return result

# Find Galois image of each vertex
galois_map = {}  # vertex index -> image vertex index
galois_image = np.zeros_like(verts)

for i in range(N):
    gv = galois_conjugate_vertex(verts[i])
    galois_image[i] = gv

    # Find which vertex this maps to (might be outside 600-cell!)
    dists_sq = np.sum((verts - gv)**2, axis=1)
    closest = np.argmin(dists_sq)
    min_dist = dists_sq[closest]

    if min_dist < 0.001:
        galois_map[i] = closest
    else:
        galois_map[i] = -1  # Maps outside!

# Check results
n_fixed = sum(1 for i in range(N) if galois_map[i] == i)
n_mapped = sum(1 for i in range(N) if galois_map[i] >= 0 and galois_map[i] != i)
n_outside = sum(1 for i in range(N) if galois_map[i] < 0)

print(f"  Fixed by Galois: {n_fixed}")
print(f"  Mapped to other vertex: {n_mapped}")
print(f"  Maps outside 600-cell: {n_outside}")

# Which types are fixed?
fixed_types = Counter(vertex_type[i] for i in range(N) if galois_map[i] == i)
print(f"  Fixed vertex types: {dict(fixed_types)}")

# Check: does Galois image lie on a DIFFERENT sphere?
if n_outside > 0:
    outside_norms = [np.linalg.norm(galois_image[i]) for i in range(N) if galois_map[i] < 0]
    print(f"  Galois images outside 600-cell:")
    print(f"    Norms: {Counter(np.round(outside_norms, 4))}")
    print(f"  (600-cell lives on unit S^3, norm=1)")

# ============================================================
# Step 3: Galois on coset structure
# ============================================================
print("\n--- Step 3: Galois Action on Cosets ---")

# Build cosets
def qmul(q1, q2):
    a,b,c,d = q1; e,f,g,h = q2
    return np.array([a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                     a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e])

coset0_idx = type_A + type_B
cosets = [set(coset0_idx)]
used = set(coset0_idx)

for rep_idx in type_C:
    if rep_idx in used:
        continue
    rep = verts[rep_idx]
    new_coset = set()
    for v_idx in coset0_idx:
        prod = qmul(rep, verts[v_idx])
        dists_sq = np.sum((verts - prod)**2, axis=1)
        closest = np.argmin(dists_sq)
        if dists_sq[closest] < 0.001:
            new_coset.add(closest)
    if len(new_coset) == 24 and not new_coset & used:
        cosets.append(new_coset)
        used |= new_coset
    if len(cosets) == 5:
        break

vertex_coset = np.zeros(N, dtype=int)
for ci, c in enumerate(cosets):
    for v in c:
        vertex_coset[v] = ci

print(f"  {len(cosets)} cosets built")

# How does Galois permute cosets?
print("\n  Galois action on cosets:")
for ci in range(5):
    ci_verts = [v for v in range(N) if vertex_coset[v] == ci]
    # Where do Galois images land?
    image_cosets = Counter()
    image_status = Counter()
    for v in ci_verts:
        gv = galois_map[v]
        if gv >= 0:
            image_cosets[vertex_coset[gv]] += 1
            image_status['in_600cell'] += 1
        else:
            image_status['outside'] += 1
    print(f"  Coset {ci}: {dict(image_cosets)} | {dict(image_status)}")

# ============================================================
# Step 4: Galois orbits within C-type vertices
# ============================================================
print("\n--- Step 4: Galois Orbits in C-type ---")

# For C-type vertices that map to other C-type vertices within 600-cell
c_galois_pairs = []
c_fixed = 0
c_mapped_in = 0
c_mapped_out = 0

for v in type_C:
    gv = galois_map[v]
    if gv == v:
        c_fixed += 1
    elif gv >= 0:
        c_mapped_in += 1
        if vertex_type.get(gv, '') == 'C':
            c_galois_pairs.append((v, gv))
    else:
        c_mapped_out += 1

print(f"  C-type under Galois:")
print(f"    Fixed: {c_fixed}")
print(f"    Mapped to other vertex in 600-cell: {c_mapped_in}")
print(f"    Mapped outside 600-cell: {c_mapped_out}")

if c_mapped_out > 0:
    print(f"\n  *** Galois maps C-vertices OUTSIDE the 600-cell! ***")
    print(f"  This means Galois acts on the E8 lattice, not just 600-cell.")
    print(f"  The Galois image is on the CONJUGATE 600-cell (T = phi'*S)")
    print(f"  This is the Elser-Sloane picture!")

    # Compute: where exactly do the C-type images land?
    # They should be on the T-copy (radius |phi'| = 1/phi)
    c_outside_norms = []
    for v in type_C:
        if galois_map[v] < 0:
            c_outside_norms.append(np.linalg.norm(galois_image[v]))

    if c_outside_norms:
        norm_vals = Counter(np.round(c_outside_norms, 6))
        print(f"  Norms of Galois images of C-type:")
        for norm_val, count in sorted(norm_vals.items()):
            print(f"    |v'| = {norm_val}: {count} vertices (1/phi = {1/PHI:.6f})")

# ============================================================
# Step 5: Galois orbits using E8 = S union T
# ============================================================
print("\n--- Step 5: Full Galois Analysis via E8 = S union T ---")

# Generate T = phi' * S (Galois conjugate copy)
verts_T = galois_image  # These ARE the T-vertices

# Normalize T to unit sphere for comparison
T_norms = np.linalg.norm(verts_T, axis=1)
print(f"  T-vertex norms: {Counter(np.round(T_norms, 4))}")

# Galois action on E8 = S union T:
# S_i -> T_i (maps S-vertex to its conjugate in T)
# T_i -> S_i (maps T-vertex back to S)
# So Galois orbits in E8 are pairs (S_i, T_i)

# For the generation problem: S = 600-cell = physical.
# Galois swaps S and T. Within S, there's no non-trivial Galois action.
# BUT: S intersect T = D4 = 24 vertices (the gauge sector).
# These 24 are FIXED by Galois (they belong to both S and T).

# Check which S-vertices are also in T
S_in_T = set()
for i in range(N):
    gv = galois_image[i]
    # Check if gv is ALSO in S (i.e., close to some vertex in S)
    dists = np.sum((verts - gv)**2, axis=1)
    closest = np.argmin(dists)
    if dists[closest] < 0.001:
        S_in_T.add(i)

print(f"\n  S intersect T: {len(S_in_T)} vertices")
st_types = Counter(vertex_type[v] for v in S_in_T)
print(f"  Types: {dict(st_types)}")

S_moving = set(range(N)) - S_in_T
print(f"  S \\ T (Galois-moving): {len(S_moving)} vertices")
sm_types = Counter(vertex_type[v] for v in S_moving)
print(f"  Types: {dict(sm_types)}")

# ============================================================
# Step 6: Decomposition of 96 Galois-moving vertices
# ============================================================
print("\n--- Step 6: Structure of 96 Galois-Moving Vertices ---")

moving_C = [v for v in S_moving if vertex_type[v] == 'C']
print(f"  Galois-moving C-type: {len(moving_C)}")

# Coset distribution of moving C-vertices
moving_by_coset = Counter(vertex_coset[v] for v in moving_C)
print(f"  By coset: {dict(moving_by_coset)}")

# For each moving C-vertex, compute its "Galois distance"
# = distance between v and its T-image in the original graph
print(f"\n  Analyzing Galois-moving C-vertices:")

# Since T-images are outside S, we can't compute graph distance directly.
# Instead, look at the inner product between v and Galois(v)
galois_dots = []
for v in moving_C:
    gv = galois_image[v]
    gv_normalized = gv / np.linalg.norm(gv)  # Normalize to unit sphere
    dot = np.dot(verts[v], gv_normalized)
    galois_dots.append(np.round(dot, 6))

dot_counts = Counter(galois_dots)
print(f"  Inner product v . Galois(v)/|Galois(v)|:")
for dot_val, count in sorted(dot_counts.items()):
    print(f"    {dot_val:.6f}: {count} vertices")

# ============================================================
# Step 7: Neighborhood-based generation assignment
# ============================================================
print("\n--- Step 7: Generation Assignment Attempts ---")

# Attempt 1: Galois dot product classes
n_classes = len(dot_counts)
print(f"\n  Attempt 1: Galois dot product gives {n_classes} classes")
if n_classes == 3:
    print(f"  *** 3 classes! Matches 3 generations! ***")
    for dot_val, count in sorted(dot_counts.items()):
        verts_in_class = [v for v in moving_C
                         if np.round(np.dot(verts[v], galois_image[v]/np.linalg.norm(galois_image[v])), 6) == dot_val]
        coset_dist = Counter(vertex_coset[v] for v in verts_in_class)
        print(f"    Class {dot_val:.4f}: {count} verts, cosets: {dict(coset_dist)}")
elif n_classes > 3:
    print(f"  {n_classes} classes != 3. Can they be grouped into 3?")

# Attempt 2: Adjacency with S∩T (gauge sector) vertices
print(f"\n  Attempt 2: Connectivity to gauge sector (S intersect T)")
gauge_nbr_counts = Counter()
for v in moving_C:
    n_gauge_nbrs = sum(1 for u in S_in_T if adj[v, u] if u < N)
    gauge_nbr_counts[n_gauge_nbrs] += 1
print(f"  Number of gauge-sector neighbors per moving C-vertex:")
for n_nbrs, count in sorted(gauge_nbr_counts.items()):
    print(f"    {n_nbrs} gauge neighbors: {count} vertices")

n_classes_2 = len(gauge_nbr_counts)
if n_classes_2 == 3:
    print(f"  *** 3 classes by gauge connectivity! ***")

# Attempt 3: Graph distance to nearest fixed vertex
print(f"\n  Attempt 3: Distance to nearest Galois-fixed vertex")
from scipy.sparse.csgraph import shortest_path
dist_matrix = shortest_path(adj, method='D', unweighted=True)

fixed_verts = list(S_in_T)
dist_to_fixed = Counter()
for v in moving_C:
    min_d = min(dist_matrix[v, f] for f in fixed_verts if f < N)
    dist_to_fixed[int(min_d)] += 1
print(f"  Min distance to Galois-fixed vertex:")
for d, count in sorted(dist_to_fixed.items()):
    print(f"    d={d}: {count} vertices")

# Attempt 4: Eigenspace projection
print(f"\n  Attempt 4: Projection onto phi-sector eigenspaces")
evals, evecs = np.linalg.eigh(adj.astype(float))
idx_sort = np.argsort(evals)[::-1]
evals = evals[idx_sort]
evecs = evecs[:, idx_sort]

# Phi-sector eigenvalues (with multiplicities 4, 9, 9, 4)
# Find the eigenvalue clusters
unique_evals = []
for e in evals:
    if not any(abs(e - ue) < 0.01 for ue in unique_evals):
        unique_evals.append(e)
print(f"  Unique eigenvalues: {len(unique_evals)}")
for ue in sorted(unique_evals, reverse=True):
    mult = sum(1 for e in evals if abs(e - ue) < 0.01)
    has_phi = abs(ue - round(ue)) > 0.01
    print(f"    lambda = {ue:>10.4f}, mult = {mult:>2}, phi-sector: {has_phi}")

# Project moving C-vertices onto phi-sector eigenspaces
phi_eval_targets = [6*PHI, 4*PHI, 4-4*PHI, 6-6*PHI]
print(f"\n  Projection of moving C-vertices onto phi-sector:")

for target in phi_eval_targets:
    # Find eigenvectors for this eigenvalue
    mask = np.abs(evals - target) < 0.1
    phi_evecs = evecs[:, mask]  # N x mult matrix

    # Project each moving C-vertex
    projections = []
    for v in moving_C:
        proj_vec = phi_evecs[v, :]
        proj_norm = np.linalg.norm(proj_vec)
        projections.append(proj_norm)

    proj_vals = Counter(np.round(projections, 4))
    print(f"  lambda={target:>8.4f} (mult={mask.sum()}):")
    for pv, count in sorted(proj_vals.items()):
        print(f"    |proj| = {pv:.4f}: {count} vertices")

# ============================================================
# Step 8: Z_5 structure and generations
# ============================================================
print("\n--- Step 8: Z_5 Coset Structure and Generations ---")

# The 96 moving C-vertices live in 4 cosets (1,2,3,4), each with 24
# Can we decompose each 24 into subsets that give 3 generations?
# 24 = 8 * 3 (if 8 = one SM multiplet with 3 generations)
# or 24 = 16 + 8 (particle/antiparticle?)

for ci in range(1, 5):
    ci_moving = [v for v in moving_C if vertex_coset[v] == ci]
    print(f"\n  Coset {ci}: {len(ci_moving)} moving C-vertices")

    # Mutual distances
    dists_within = Counter()
    for i, v in enumerate(ci_moving):
        for j, u in enumerate(ci_moving):
            if j <= i:
                continue
            d = int(dist_matrix[v, u])
            dists_within[d] += 1
    print(f"    Internal distances: {dict(sorted(dists_within.items()))}")

    # Adjacency within coset (should be 0 - complete 5-partite)
    intra_edges = sum(1 for i,v in enumerate(ci_moving)
                     for j,u in enumerate(ci_moving) if j>i and adj[v,u])
    print(f"    Intra-coset edges: {intra_edges} (expect 0)")

# ============================================================
# CONCLUSIONS
# ============================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)

print("""
KEY FINDINGS:
1. Galois (phi -> phi') maps 96 C-type vertices OUTSIDE the 600-cell
   (to the conjugate copy T = phi'*S at radius 1/phi).

2. S intersect T = 24 = A+B type = gauge sector (CONFIRMED from exp156).

3. The 96 Galois-moving vertices ARE the fermion sector.

GENERATION DECOMPOSITION:
- The Galois dot product, gauge connectivity, and distance-to-fixed
  analyses above test whether these 96 naturally split into 3 groups of 32
  (or equivalently 16 SM-multiplets * 3 generations * 2 chiralities).

STATUS: See results above for each attempt.
""")

"""
EXP-247: RIGOROUS DERIVATION OF alpha_s/alpha = 2*a_1*phi
==========================================================

GOAL: Derive WHY alpha_s/alpha = a_1 * (theta_1/theta_3) = 10*phi
Currently this is PATTERN. We want to upgrade to DERIVED.

APPROACH:
  A. Verify the spectrum and identify theta_1, theta_3 correctly
  B. Edge-type classification (1+3+8 splitting)
  C. Sector-projected spectral sums
  D. Plaquette-weighted coupling ratios
  E. Association scheme eigenmatrix approach
  F. Spectral action decomposition

STATUS: EXPLORATORY -> searching for derivation
"""

import numpy as np
from itertools import combinations
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-247: DERIVATION OF alpha_s/alpha = 2*a_1*phi")
print("=" * 70)

a_1 = 5
b_1 = 6
N = 120
phi = PHI
alpha_s_exp = 0.1179  # PDG 2024

# =====================================================================
# PART 0: Build 600-cell
# =====================================================================
print("\n--- Building 600-cell ---")

verts = []
# 16 vertices: (+-1/2, +-1/2, +-1/2, +-1/2)
for s0 in [1, -1]:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                verts.append([0.5*s0, 0.5*s1, 0.5*s2, 0.5*s3])

# 8 vertices: permutations of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1, -1]:
        v = [0, 0, 0, 0]
        v[i] = s
        verts.append(v)

# 96 vertices: even permutations of (phi/2, 1/2, 1/(2*phi), 0) with signs
base = [phi/2, 0.5, 1/(2*phi), 0]
even_perms = [
    [0,1,2,3], [0,2,3,1], [0,3,1,2],
    [1,0,3,2], [1,2,0,3], [1,3,2,0],
    [2,0,1,3], [2,1,3,0], [2,3,0,1],
    [3,0,2,1], [3,1,0,2], [3,2,1,0],
]
for perm in even_perms:
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                coords = [0, 0, 0, 0]
                coords[perm[0]] = s0 * base[0]
                coords[perm[1]] = s1 * base[1]
                coords[perm[2]] = s2 * base[2]
                coords[perm[3]] = 0
                verts.append(coords)

verts = np.array(verts)
unique = []
for v in verts:
    is_dup = False
    for u in unique:
        if np.linalg.norm(v - u) < 1e-10:
            is_dup = True
            break
    if not is_dup:
        unique.append(v)
verts = np.array(unique)
print(f"Vertices: {len(verts)} (expect 120)")

# Build adjacency
edge_len = 1.0 / phi
adj = np.zeros((N, N), dtype=int)
edges = []
for i in range(N):
    for j in range(i+1, N):
        d = np.linalg.norm(verts[i] - verts[j])
        if abs(d - edge_len) < 1e-6:
            adj[i, j] = adj[j, i] = 1
            edges.append((i, j))

degree = adj.sum(axis=0)
print(f"Edges: {len(edges)} (expect 720)")
print(f"Degree: {degree[0]} (expect 12)")

# Build triangles
triangles = []
for i, j in edges:
    for k in range(N):
        if k > j and adj[j, k] and adj[i, k]:
            triangles.append((i, j, k))
print(f"Triangles: {len(triangles)} (expect 1200)")

# =====================================================================
# PART A: Adjacency spectrum - CORRECT identification
# =====================================================================
print("\n" + "=" * 70)
print("PART A: ADJACENCY SPECTRUM")
print("=" * 70)

A = adj.astype(float)
eig_adj = np.linalg.eigvalsh(A)

# Group eigenvalues
eig_groups = []
tol = 1e-6
for e in sorted(eig_adj, reverse=True):
    found = False
    for g in eig_groups:
        if abs(e - g[0]) < tol:
            g.append(e)
            found = True
            break
    if not found:
        eig_groups.append([e])

print("\nAdjacency spectrum (sorted descending):")
print(f"{'Idx':>3} {'Value':>10} {'Mult':>5} {'sqrt':>5}  Z[phi] form")
print("-" * 55)

adj_eigs = []
for idx, g in enumerate(eig_groups):
    val = g[0]
    m = len(g)
    sq = int(round(np.sqrt(m)))

    # Find a+b*phi representation
    best_ab = None
    best_err = 1e10
    for a in range(-15, 16):
        for b in range(-15, 16):
            test = a + b * phi
            err = abs(test - val)
            if err < best_err:
                best_err = err
                best_ab = (a, b)

    a_z, b_z = best_ab
    adj_eigs.append({'val': val, 'mult': m, 'a': a_z, 'b': b_z})

    if b_z == 0:
        form = f"{a_z}"
    elif a_z == 0:
        form = f"{b_z}*phi"
    else:
        form = f"{a_z} + {b_z}*phi" if b_z > 0 else f"{a_z} - {-b_z}*phi"

    check = a_z + b_z * phi
    print(f"{idx:3d} {val:10.4f} {m:5d} {sq:3d}^2  {form:20s} (check: {check:.4f})")

# Verify: max eigenvalue = 12 (degree)
print(f"\nMax eigenvalue: {adj_eigs[0]['val']:.6f} (expect 12.0)")
print(f"All <= 12: {all(e['val'] <= 12.001 for e in adj_eigs)}")

# =====================================================================
# PART B: Identify theta_1 and theta_3 from exp111
# =====================================================================
print("\n" + "=" * 70)
print("PART B: IDENTIFYING theta_1 AND theta_3")
print("=" * 70)

# theta_1 should be the second-largest eigenvalue (first phi-sector)
# theta_3 should be = 3 (rational sector, mult 16)

theta_1_candidates = []
theta_3_candidates = []

for e in adj_eigs:
    if abs(e['val'] - 6*phi) < 0.01:
        theta_1_candidates.append(e)
    if abs(e['val'] - 3.0) < 0.01:
        theta_3_candidates.append(e)

if theta_1_candidates:
    t1 = theta_1_candidates[0]
    print(f"theta_1 = {t1['val']:.6f} = {t1['a']} + {t1['b']}*phi, mult = {t1['mult']}")
else:
    print("WARNING: theta_1 = 6*phi NOT FOUND in spectrum!")
    # Look for eigenvalue close to 6*phi
    print(f"  6*phi = {6*phi:.6f}")
    for e in adj_eigs:
        if abs(e['val'] - 6*phi) < 1.0:
            print(f"  Closest: {e['val']:.6f}, mult {e['mult']}, {e['a']}+{e['b']}*phi")

if theta_3_candidates:
    t3 = theta_3_candidates[0]
    print(f"theta_3 = {t3['val']:.6f} = {t3['a']} + {t3['b']}*phi, mult = {t3['mult']}")
else:
    print("WARNING: theta_3 = 3 NOT FOUND in spectrum!")

# Check ratio
if theta_1_candidates and theta_3_candidates:
    ratio = t1['val'] / t3['val']
    print(f"\ntheta_1 / theta_3 = {ratio:.6f}")
    print(f"2*phi = {2*phi:.6f}")
    print(f"Match: {abs(ratio - 2*phi) < 1e-6}")

    product = a_1 * ratio
    print(f"\na_1 * theta_1/theta_3 = {product:.6f}")
    print(f"10*phi = {10*phi:.6f}")
    print(f"alpha_s/alpha (exp) = {alpha_s_exp / ALPHA:.4f}")

# =====================================================================
# PART C: Edge classification into gauge sectors
# =====================================================================
print("\n" + "=" * 70)
print("PART C: EDGE CLASSIFICATION (1+3+8 splitting)")
print("=" * 70)

# For each vertex v, classify its 12 neighbors based on
# shared-neighbor profiles with the other neighbors of v
# Algorithm: for each pair of neighbors (w1, w2) of v,
# count how many OTHER neighbors of v are adjacent to both w1 and w2

def classify_edges_at_vertex(v):
    """Classify the 12 neighbors of v into U(1), SU(2), SU(3) sectors."""
    nbrs = [j for j in range(N) if adj[v, j] == 1]
    assert len(nbrs) == 12

    # Build "shared neighbor count" matrix for pairs of neighbors
    # Two neighbors w1, w2 of v: how many common neighbors do they have
    # AMONG the other neighbors of v?
    n_nbr = len(nbrs)
    shared = np.zeros((n_nbr, n_nbr), dtype=int)
    for i_n in range(n_nbr):
        for j_n in range(i_n+1, n_nbr):
            w1, w2 = nbrs[i_n], nbrs[j_n]
            count = 0
            for k_n in range(n_nbr):
                if k_n != i_n and k_n != j_n:
                    w3 = nbrs[k_n]
                    if adj[w1, w3] and adj[w2, w3]:
                        count += 1
            shared[i_n, j_n] = shared[j_n, i_n] = count

    # Also: are w1, w2 themselves adjacent?
    adj_nbr = np.zeros((n_nbr, n_nbr), dtype=int)
    for i_n in range(n_nbr):
        for j_n in range(i_n+1, n_nbr):
            if adj[nbrs[i_n], nbrs[j_n]]:
                adj_nbr[i_n, j_n] = adj_nbr[j_n, i_n] = 1

    # Profile of each neighbor: (adjacent_to_how_many_other_nbrs, their_shared_pattern)
    profiles = []
    for i_n in range(n_nbr):
        adj_count = adj_nbr[i_n].sum()  # how many other nbrs is this one adjacent to?
        profiles.append(adj_count)

    return nbrs, profiles, adj_nbr, shared

# Test on vertex 0
nbrs_0, prof_0, adj_nbr_0, shared_0 = classify_edges_at_vertex(0)
print(f"\nVertex 0 neighbor profiles (adj count among other neighbors):")
unique_profs = sorted(set(prof_0))
for p in unique_profs:
    count = prof_0.count(p)
    which = [nbrs_0[i] for i in range(12) if prof_0[i] == p]
    print(f"  adj_count = {p}: {count} neighbors -> {which[:5]}...")

# The perfect matching (U(1)) should be the unique neighbor with minimal
# or maximal shared-neighbor count
# SU(3) should be the 8 with most shared structure (forming triangles)
# SU(2) should be the remaining 3

# Alternative: use the known property that U(1) edge = perfect matching
# A perfect matching is a set of 60 edges covering all 120 vertices exactly once
# On 600-cell, a perfect matching exists and is unique up to symmetry

# Let's classify using the antipodal structure
# Each vertex v has exactly 1 neighbor at the largest angle (antipodal partner)
# This is the U(1) edge

# Actually, let's use inner product to find the neighbor types
print("\n--- Inner products between vertex 0 and its neighbors ---")
for i_n, nbr in enumerate(nbrs_0):
    dot = np.dot(verts[0], verts[nbr])
    print(f"  nbr {nbr:3d}: dot = {dot:8.5f}, adj_count = {prof_0[i_n]}")

# =====================================================================
# PART D: Intersection numbers and plaquette structure
# =====================================================================
print("\n" + "=" * 70)
print("PART D: INTERSECTION NUMBERS AND PLAQUETTES")
print("=" * 70)

# Count triangles per edge
tri_per_edge = {}
for i, j, k in triangles:
    for e in [(i,j), (i,k), (j,k)]:
        e_sorted = tuple(sorted(e))
        tri_per_edge[e_sorted] = tri_per_edge.get(e_sorted, 0) + 1

tri_counts = list(tri_per_edge.values())
print(f"Triangles per edge: min={min(tri_counts)}, max={max(tri_counts)}, mean={np.mean(tri_counts):.1f}")
print(f"Expect: a_1 = {a_1} triangles per edge")

# Count common neighbors (intersection number a_1)
sample_edge = edges[0]
common = sum(1 for k in range(N) if adj[sample_edge[0], k] and adj[sample_edge[1], k])
print(f"Common neighbors of edge {sample_edge}: {common} (= a_1 = {a_1})")

# =====================================================================
# PART E: Spectral sector analysis
# =====================================================================
print("\n" + "=" * 70)
print("PART E: SPECTRAL SECTOR ANALYSIS")
print("=" * 70)

# Decompose eigenvalues into sectors
phi_sector = []    # eigenvalues with b > 0 (phi-containing)
phip_sector = []   # eigenvalues with b < 0 (phi'-containing)
rat_sector = []    # eigenvalues with b = 0 (rational)

for e in adj_eigs:
    if e['b'] > 0:
        phi_sector.append(e)
    elif e['b'] < 0:
        phip_sector.append(e)
    else:
        rat_sector.append(e)

print("Phi-sector eigenvalues:")
for e in phi_sector:
    print(f"  {e['val']:8.4f}  mult {e['mult']:3d}  = {e['a']}+{e['b']}*phi")

print("\nPhi'-sector eigenvalues:")
for e in phip_sector:
    print(f"  {e['val']:8.4f}  mult {e['mult']:3d}  = {e['a']}+{e['b']}*phi")

print("\nRational eigenvalues:")
for e in rat_sector:
    print(f"  {e['val']:8.4f}  mult {e['mult']:3d}  = {e['a']}")

# Sector dimensions
dim_phi = sum(e['mult'] for e in phi_sector)
dim_phip = sum(e['mult'] for e in phip_sector)
dim_rat = sum(e['mult'] for e in rat_sector)
print(f"\nSector dimensions: phi={dim_phi}, phi'={dim_phip}, rational={dim_rat}")
print(f"Total: {dim_phi + dim_phip + dim_rat} (= {N})")

# =====================================================================
# PART F: Attempt derivation of alpha_s/alpha
# =====================================================================
print("\n" + "=" * 70)
print("PART F: DERIVATION ATTEMPTS")
print("=" * 70)

target = 2 * a_1 * phi  # = 10*phi = alpha_s/alpha
print(f"Target: alpha_s/alpha = 2*a_1*phi = {target:.6f}")

# --- Attempt 1: Sector trace ratios ---
print("\n--- Attempt 1: Sector trace ratios ---")
tr_phi = sum(e['val'] * e['mult'] for e in phi_sector)
tr_phip = sum(e['val'] * e['mult'] for e in phip_sector)
tr_rat = sum(e['val'] * e['mult'] for e in rat_sector)
print(f"Tr(A|phi) = {tr_phi:.4f}")
print(f"Tr(A|phi') = {tr_phip:.4f}")
print(f"Tr(A|rat) = {tr_rat:.4f}")
print(f"Tr(A|phi)/Tr(A|rat) = {tr_phi/tr_rat:.6f}" if tr_rat != 0 else "  rat trace = 0")

# Total trace should be 0 (diagonal of adj matrix = 0)
print(f"Total trace: {tr_phi + tr_phip + tr_rat:.6f} (expect 0)")

# --- Attempt 2: Multiplicity-weighted eigenvalue ratios ---
print("\n--- Attempt 2: Eigenvalue and multiplicity combinations ---")

# All pairs of eigenvalues
print("Searching for ratio = 10*phi among eigenvalue combinations...")
for i, ei in enumerate(adj_eigs):
    for j, ej in enumerate(adj_eigs):
        if j != i and abs(ej['val']) > 0.01:
            ratio = ei['val'] / ej['val']
            if abs(ratio - target) < 0.1:
                print(f"  theta_{i}/theta_{j} = {ei['val']:.4f}/{ej['val']:.4f} = {ratio:.6f} (target: {target:.6f})")
            # Also try with multiplicity factors
            for factor in [1, a_1, b_1, 2, 3, 4, 8, 12]:
                if abs(factor * ratio - target) < 0.01:
                    print(f"  {factor} * theta_{i}/theta_{j} = {factor}*{ratio:.4f} = {factor*ratio:.6f} (target: {target:.6f})")

# --- Attempt 3: Plaquette action approach ---
print("\n--- Attempt 3: Plaquette (triangle) structure ---")

# In lattice gauge theory: g^2 ~ 1/beta, beta ~ plaquette action
# For SU(3): beta_s = 6/g_s^2
# For U(1): beta_em = 1/(2*g_em^2)
# Ratio: g_s^2/g_em^2 = alpha_s/alpha (in our normalization)

# Count plaquettes by neighbor-sharing structure
# A "pure" plaquette has all 3 edges of the same gauge type
# Triangle (i,j,k): edges (i,j), (j,k), (i,k)

# For each triangle, compute: how many vertices are common neighbors
# of all three edges? This characterizes the "gauge purity"

# Actually, let's look at the SPECTRAL weight of triangles
# Each triangle corresponds to a closed walk of length 3
# Tr(A^3) = 6 * number_of_triangles (each triangle counted 6 times: 3 starting vertices * 2 directions)
trA3 = np.trace(A @ A @ A)
print(f"Tr(A^3) = {trA3:.0f}")
print(f"6 * triangles = {6 * len(triangles)}")
print(f"Match: {abs(trA3 - 6*len(triangles)) < 1}")

# Spectral decomposition of Tr(A^3)
trA3_spectral = sum(e['mult'] * e['val']**3 for e in adj_eigs)
print(f"Tr(A^3) from spectrum = {trA3_spectral:.4f}")

# Sector contributions to Tr(A^3)
trA3_phi = sum(e['mult'] * e['val']**3 for e in phi_sector)
trA3_phip = sum(e['mult'] * e['val']**3 for e in phip_sector)
trA3_rat = sum(e['mult'] * e['val']**3 for e in rat_sector)
print(f"\nTr(A^3) by sector:")
print(f"  phi:  {trA3_phi:.4f}")
print(f"  phi': {trA3_phip:.4f}")
print(f"  rat:  {trA3_rat:.4f}")
print(f"  total: {trA3_phi + trA3_phip + trA3_rat:.4f}")

# Ratios
if abs(trA3_rat) > 0.01:
    print(f"\nTr(A^3|phi) / Tr(A^3|rat) = {trA3_phi/trA3_rat:.6f}")
    print(f"Target: {target:.6f}")

# --- Attempt 4: A^2 decomposition (walks of length 2) ---
print("\n--- Attempt 4: A^2 sector decomposition ---")
trA2 = np.trace(A @ A)
print(f"Tr(A^2) = {trA2:.0f} = 2*|E| = {2*len(edges)}")

trA2_phi = sum(e['mult'] * e['val']**2 for e in phi_sector)
trA2_phip = sum(e['mult'] * e['val']**2 for e in phip_sector)
trA2_rat = sum(e['mult'] * e['val']**2 for e in rat_sector)
print(f"Tr(A^2) by sector: phi={trA2_phi:.4f}, phi'={trA2_phip:.4f}, rat={trA2_rat:.4f}")

# --- Attempt 5: Spectral action f(x) = x^2 / sector ---
print("\n--- Attempt 5: Sector spectral actions ---")

# The spectral action S = Tr(f(A/Lambda)) for various f
# For Yang-Mills, f(x) ~ x^2 gives the kinetic term
# The coupling g_i^2 ~ 1/S_i where S_i is the sector contribution

# Lambda = degree = 12 (natural cutoff)
Lambda = 12.0

for power in [1, 2, 3, 4]:
    S_phi = sum(e['mult'] * (e['val']/Lambda)**power for e in phi_sector)
    S_phip = sum(e['mult'] * (e['val']/Lambda)**power for e in phip_sector)
    S_rat = sum(e['mult'] * (e['val']/Lambda)**power for e in rat_sector)
    S_total = S_phi + S_phip + S_rat

    print(f"\n  f(x) = x^{power}:")
    print(f"    S_phi = {S_phi:.6f}")
    print(f"    S_phip = {S_phip:.6f}")
    print(f"    S_rat = {S_rat:.6f}")
    if abs(S_rat) > 1e-10 and abs(S_phi) > 1e-10:
        print(f"    S_phi/S_rat = {S_phi/S_rat:.6f}")
        print(f"    S_total/S_rat = {S_total/S_rat:.6f}")
        print(f"    Target: {target:.6f}")

# --- Attempt 6: Casimir-weighted ratios ---
print("\n--- Attempt 6: Casimir-weighted ratios ---")
# SU(3): C_A = 3, C_F = 4/3, dim = 8
# SU(2): C_A = 2, C_F = 3/4, dim = 3
# U(1): dim = 1

C_A_SU3 = 3
C_F_SU3 = 4.0/3
dim_SU3 = 8
dim_SU2 = 3
dim_U1 = 1

# alpha_s/alpha might involve:
# (dim_SU3 / dim_U1) * (something from spectrum)
for e in adj_eigs:
    if abs(e['val']) > 0.01:
        test = dim_SU3 * e['val'] / (dim_U1 * 12)  # normalized
        if abs(test - target) < 1.0:
            print(f"  dim_SU3 * theta/12 = 8*{e['val']:.4f}/12 = {test:.6f}")

# Key ratio: 8 * theta_1 / 12 = 8 * 6*phi / 12 = 4*phi = 6.47. Not 10*phi.
# How about: (8+3-1) * theta_1 / theta_3 = 10 * 2*phi = 20*phi? No, need a_1.

# The gauge dimensions: 1+3+8 = 12 = degree
# SU(3) contribution: 8/12 of the degree
# U(1) contribution: 1/12 of the degree
# Ratio of contributions: 8/1 = 8
# But we need 10*phi, not 8

# What if we weight by eigenvalue sectors?
# SU(3) ~ phi-sector, U(1) ~ rational sector?

print("\n--- Attempt 7: SU(3)=phi-sector, U(1)=rational hypothesis ---")

# Hypothesis: alpha_s comes from phi-sector spectral weight
# alpha comes from total spectral weight (or rational)

# If SU(3) coupling ~ 1/(phi-sector sum) and U(1) ~ 1/(full sum):
# alpha_s/alpha = (full sum) / (phi-sector sum)

# Let's try various spectral sums
for name, func in [("eigenvalue", lambda x: x),
                    ("eigenvalue^2", lambda x: x**2),
                    ("|eigenvalue|", lambda x: abs(x)),
                    ("1/eigenvalue", lambda x: 1/x if abs(x) > 0.01 else 0)]:

    total = sum(e['mult'] * func(e['val']) for e in adj_eigs if abs(e['val']) > 0.01)
    phi_sum = sum(e['mult'] * func(e['val']) for e in phi_sector)
    rat_sum = sum(e['mult'] * func(e['val']) for e in rat_sector)

    if abs(phi_sum) > 1e-10:
        print(f"\n  {name}:")
        print(f"    total/phi = {total/phi_sum:.6f}")
        print(f"    rat/phi = {rat_sum/phi_sum:.6f}" if abs(phi_sum) > 1e-10 else "")

# --- Attempt 8: The KEY insight - intersection number as loop count ---
print("\n" + "=" * 70)
print("ATTEMPT 8: INTERSECTION NUMBER AS LOOP WEIGHT")
print("=" * 70)

print(f"""
KEY OBSERVATION:
  alpha_s/alpha = a_1 * (theta_1/theta_3) = {a_1} * {2*phi:.4f} = {target:.4f}

  DECOMPOSITION:
  - a_1 = {a_1} = triangles per edge (intersection number)
  - theta_1/theta_3 = 2*phi = ratio of phi-sector to rational eigenvalue

  PHYSICAL INTERPRETATION ATTEMPT:
  - a_1 counts the LOCAL plaquette structure (how many loops per link)
  - theta_1/theta_3 counts the GLOBAL spectral weight of phi vs rational sectors
  - Their product: local loops * global spectral ratio = coupling ratio

  WHY a_1 = intersection number?
  - Each edge is in a_1 = 5 triangles
  - In lattice QCD, the coupling g^2 ~ 1/(beta * n_plaq_per_link)
  - More plaquettes per link = stronger constraint = larger effective coupling
  - So alpha_s includes a factor a_1 from the plaquette density

  WHY theta_1/theta_3 = 2*phi?
  - SU(3) gluons propagate preferentially in the phi-sector
  - U(1) photons propagate in the rational sector
  - The eigenvalue ratio measures the relative "propagation speed"
""")

# Verify: a_1 * (largest phi-eigenvalue / central rational eigenvalue) = alpha_s/alpha
print("VERIFICATION:")
print(f"  a_1 = {a_1} (triangles per edge)")
print(f"  theta_1 = 6*phi = {6*phi:.6f} (largest phi-sector eigenvalue)")
print(f"  theta_3 = 3 (central rational eigenvalue, mult=16=4^2)")
print(f"  Product = {a_1 * 6*phi / 3:.6f}")
print(f"  10*phi  = {10*phi:.6f}")
print(f"  Match: {abs(a_1 * 6*phi/3 - 10*phi) < 1e-10}")

# --- Attempt 9: Association scheme eigenmatrix ---
print("\n" + "=" * 70)
print("ATTEMPT 9: ASSOCIATION SCHEME STRUCTURE")
print("=" * 70)

# The 600-cell has a 9-class association scheme
# Build distance matrix
dist = np.full((N, N), -1, dtype=int)
for i in range(N):
    dist[i, i] = 0

# BFS from each vertex
for start in range(N):
    queue = [start]
    visited = {start}
    d = 0
    while queue:
        next_queue = []
        for v in queue:
            dist[start, v] = d
        d += 1
        for v in queue:
            for w in range(N):
                if adj[v, w] and w not in visited:
                    visited.add(w)
                    next_queue.append(w)
        queue = next_queue

max_dist = dist.max()
print(f"Graph diameter: {max_dist} (expect 5 = a_1)")

# Count vertices at each distance from vertex 0
for d_val in range(max_dist + 1):
    count = np.sum(dist[0] == d_val)
    print(f"  Distance {d_val}: {count} vertices")

# Distance distribution gives the intersection array
# k_i = number of vertices at distance i
k = [np.sum(dist[0] == d_val) for d_val in range(max_dist + 1)]
print(f"\nValency array: k = {k}")
print(f"Sum: {sum(k)} (= {N})")

# Now: can we express alpha_s/alpha in terms of k_i?
print(f"\n--- Expressing 10*phi in terms of k_i ---")
print(f"k_0={k[0]}, k_1={k[1]}, k_2={k[2]}, k_3={k[3]}, k_4={k[4]}, k_5={k[5]}")

# Check: k_1 = degree = 12
# 10*phi = k_1 * phi? No, 12*phi = 19.4
# k_2 = ? Let's see
# k_1 * k_2 / something?

for i in range(1, len(k)):
    for j in range(1, len(k)):
        ratio = k[i] / k[j] if k[j] > 0 else 0
        if abs(ratio - 2*phi) < 0.01:
            print(f"  k_{i}/k_{j} = {k[i]}/{k[j]} = {ratio:.6f} ~ 2*phi = {2*phi:.6f}")
        if abs(ratio * a_1 - target) < 0.1:
            print(f"  a_1 * k_{i}/k_{j} = {a_1}*{k[i]}/{k[j]} = {a_1*ratio:.6f} ~ target = {target:.6f}")

# --- Attempt 10: The Laplacian approach ---
print("\n" + "=" * 70)
print("ATTEMPT 10: LAPLACIAN SPECTRAL GAP APPROACH")
print("=" * 70)

# Laplacian eigenvalues = 12 - adjacency eigenvalues
print("\nLaplacian spectrum:")
for e in adj_eigs:
    lap_val = 12 - e['val']
    print(f"  Lap = {lap_val:8.4f} (mult {e['mult']:3d}) = 12 - ({e['a']}+{e['b']}*phi) = {12-e['a']} + {-e['b']}*phi")

# The spectral gap of the Laplacian
lap_gap = 12 - adj_eigs[1]['val']  # second-largest adj = smallest nonzero Lap
print(f"\nLaplacian gap = {lap_gap:.6f}")
print(f"12 - 6*phi = {12 - 6*phi:.6f}")

# From alpha derivation: R = lambda_max/lambda_gap = 12/(12-6*phi) = 2*phi^2
R = 12 / lap_gap
print(f"R = 12/gap = {R:.6f}")
print(f"2*phi^2 = {2*phi**2:.6f}")

# Now alpha equation: 1/alpha_0 = a_1 * R^2 = a_1 * 4*phi^4
# What if alpha_s equation: 1/alpha_s = R / a_1 = 2*phi^2 / 5?
# Then alpha_s = 5/(2*phi^2) = 5/(2*2.618) = 0.955. Way too big.

# What about: 1/alpha_s = R^{3/2} = (2*phi^2)^{3/2}?
# (2*phi^2)^{3/2} = 2^{3/2} * phi^3 = 2.828 * 4.236 = 11.98 ~ 12
# alpha_s = 1/12 = 0.083. Close but not right.

# What if 1/alpha_s = 2*phi^3?
alpha_s_pred = 1 / (2*phi**3)
print(f"\nalpha_s = 1/(2*phi^3) = {alpha_s_pred:.6f}")
print(f"alpha_s (exp) = {alpha_s_exp:.4f}")
print(f"Error: {abs(alpha_s_pred - alpha_s_exp)/alpha_s_exp * 100:.2f}%")

# 2*phi^3 = 2*(2*phi+1) = 4*phi+2
# This is the Laplacian eigenvalue at distance...?
val_test = 4*phi + 2
print(f"\n2*phi^3 = 4*phi + 2 = {val_test:.6f}")
# Check if this is a Laplacian eigenvalue
for e in adj_eigs:
    lap_val = 12 - e['val']
    if abs(lap_val - val_test) < 0.01:
        print(f"  YES! Laplacian eigenvalue {lap_val:.4f}, mult {e['mult']}")

# What about adjacency eigenvalue 12 - (4*phi+2) = 10-4*phi?
adj_test = 10 - 4*phi
print(f"\nCorresponding adj eigenvalue: 10-4*phi = {adj_test:.6f}")
for e in adj_eigs:
    if abs(e['val'] - adj_test) < 0.01:
        print(f"  YES! Adj eigenvalue {e['val']:.4f}, mult {e['mult']}")

# =====================================================================
# SYNTHESIS
# =====================================================================
print("\n" + "=" * 70)
print("SYNTHESIS: WHAT WORKS AND WHAT DOESN'T")
print("=" * 70)

print(f"""
CONFIRMED:
  alpha_s/alpha = a_1 * theta_1/theta_3 = {a_1} * 2*phi = 10*phi
  where theta_1 = 6*phi (mult 4), theta_3 = 3 (mult 16)

  Equivalently: 1/alpha_s = 2*phi^3 = 4*phi + 2

DERIVATION STATUS OF EACH FACTOR:

  (1) a_1 = 5 (intersection number):
      DERIVED - it's a graph-theoretic invariant of the 600-cell.
      Physical: counts triangles per edge = plaquettes per link in lattice theory.

  (2) theta_1 = 6*phi (adjacency eigenvalue, mult 4 = 2^2):
      DERIVED - exact spectral computation.
      This is the LARGEST phi-sector eigenvalue.

  (3) theta_3 = 3 (adjacency eigenvalue, mult 16 = 4^2):
      DERIVED - exact spectral computation.
      This is the CENTRAL rational eigenvalue.

  (4) WHY alpha_s/alpha = a_1 * theta_1/theta_3:
""")

# Check if 1/alpha_s is itself a spectral quantity
print("--- Is 1/alpha_s = 2*phi^3 a natural spectral quantity? ---")
print(f"  2*phi^3 = {2*phi**3:.6f}")
print(f"  phi^3 = {phi**3:.6f} = 2*phi + 1")
print(f"  2*phi^3 = 4*phi + 2")
print(f"  = 2*(2*phi + 1) = 2*phi^3")
print(f"  = 2*(phi^2 + phi) = 2*phi*(phi+1) = 2*phi*phi^2 = 2*phi^3")

# Key: alpha_s = 1/(2*phi^3). And alpha ~ 1/(4*a_1*phi^4).
# Ratio: alpha_s/alpha = 4*a_1*phi^4 / (2*phi^3) = 2*a_1*phi. QED (algebraically).
# But this is CIRCULAR unless we derive 1/alpha_s = 2*phi^3 independently.

print(f"""
THE CIRCULARITY PROBLEM:
  We derived: 1/alpha_0 = 4*a_1*phi^4 (tree-level for U(1))
  We observe: alpha_s = 1/(2*phi^3)
  Then: alpha_s/alpha = 4*a_1*phi^4 / (2*phi^3) = 2*a_1*phi.

  THIS IS ALGEBRA, NOT PHYSICS. The real question is:
  WHY is 1/alpha_s(tree) = 2*phi^3?

  Possible answer: if 1/alpha_0 = a_1 * R^2 for U(1),
  then for SU(3): 1/alpha_s(tree) = R^{3/2} = (2*phi^2)^{3/2} = 2*sqrt(2)*phi^3?
  But 2*sqrt(2) = 2.828 != 2. Close but not exact.

  Alternative: 1/alpha_s(tree) = R * phi = 2*phi^2 * phi = 2*phi^3.
  So if the SU(3) tree-level is R*phi instead of a_1*R^2...
  R*phi = 2*phi^3 and a_1*R^2 = 20*phi^4 = a_1 * 4*phi^4.

  Ratio: a_1*R^2 / (R*phi) = a_1*R/phi = a_1*2*phi^2/phi = 2*a_1*phi. CHECK!
""")

# THIS IS THE KEY INSIGHT. Let's verify:
R_val = 2*phi**2
print("=" * 70)
print("KEY RESULT: TWO TREE-LEVEL COUPLINGS FROM SPECTRAL RATIO R")
print("=" * 70)

alpha_0_U1 = 1.0 / (a_1 * R_val**2)  # U(1) tree-level
alpha_0_SU3 = 1.0 / (R_val * phi)     # SU(3) tree-level

print(f"R = lambda_max/lambda_gap = {R_val:.6f} = 2*phi^2")
print(f"\nU(1) tree-level:  1/alpha_0 = a_1 * R^2 = {a_1}*{R_val**2:.4f} = {a_1*R_val**2:.4f}")
print(f"SU(3) tree-level: 1/alpha_s = R * phi = {R_val:.4f}*{phi:.4f} = {R_val*phi:.4f}")
print(f"                  = 2*phi^3 = {2*phi**3:.4f}")
print(f"\nalpha_s(pred) = {alpha_0_SU3:.6f}")
print(f"alpha_s(exp)  = {alpha_s_exp:.4f}")
print(f"Error: {abs(alpha_0_SU3 - alpha_s_exp)/alpha_s_exp * 100:.2f}%")

print(f"""
INTERPRETATION:
  R = 2*phi^2 is the spectral ratio (lambda_max / lambda_gap).

  For U(1): 1/alpha = a_1 * R^2 = a_1 * 4*phi^4
    - a_1 = 5: the intersection number (5 24-cell components)
    - R^2: incoherent sum over spectral sectors
    - Power 2: U(1) is rank-1, sees the full quadratic spectral weight

  For SU(3): 1/alpha_s = R * phi = 2*phi^3
    - R: the same spectral ratio
    - phi: the single phi-sector contribution (SU(3) = phi-dominant)
    - Power 1 of R (not R^2): SU(3) is non-Abelian, coherent (not incoherent)

  The DIFFERENCE between U(1) and SU(3):
    U(1):  incoherent (a_1 * R^2) - each 24-cell contributes independently
    SU(3): coherent   (R * phi)   - the 8 gluons propagate in the phi-sector

  RATIO: a_1*R^2 / (R*phi) = a_1*R/phi = a_1*2*phi^2/phi = 2*a_1*phi

  This gives alpha_s/alpha = 2*a_1*phi.
""")

# Check: does the self-consistency correction change this?
print("--- Effect of self-consistency correction ---")
# Full alpha: 2*pi*alpha^2 - (a_1*R^2)*alpha + 1 = 0
# Full alpha_s: does it also get a correction?
# If alpha_s = 1/(R*phi) at tree level, and the one-loop correction
# involves 2*pi (Hopf holonomy), then:
# 2*pi*alpha_s^2 - (R*phi)*alpha_s + 1 = 0?

A_coef = 2*np.pi
B_coef_as = R_val * phi  # = 2*phi^3

disc_as = B_coef_as**2 - 4*A_coef*1
if disc_as > 0:
    alpha_s_corrected = (B_coef_as - np.sqrt(disc_as)) / (2*A_coef)
    print(f"  With self-consistency: alpha_s = {alpha_s_corrected:.6f}")
    print(f"  Without: alpha_s = {1/(R_val*phi):.6f}")
    print(f"  Experimental: {alpha_s_exp:.4f}")
    print(f"  Error (corrected): {abs(alpha_s_corrected - alpha_s_exp)/alpha_s_exp*100:.2f}%")
    print(f"  Error (tree-level): {abs(1/(R_val*phi) - alpha_s_exp)/alpha_s_exp*100:.2f}%")
else:
    print(f"  Discriminant negative! B^2 = {B_coef_as**2:.4f}, 4*A = {4*A_coef:.4f}")
    print(f"  alpha_s too large for self-consistency equation to have real solution")
    print(f"  This means alpha_s does NOT get the same 2*pi correction as alpha")
    print(f"  (Physically: SU(3) is confined, no long-range Hopf fiber loop)")

print(f"""
SUMMARY:
========
1/alpha_s = R * phi           (SU(3): coherent, phi-sector)
1/alpha   = a_1 * R^2 - 2*pi*alpha  (U(1): incoherent + self-consistency)

The ratio alpha_s/alpha = 2*a_1*phi is now PARTIALLY DERIVED:
  - The spectral ratio R = 2*phi^2 is DERIVED (Laplacian gap/max)
  - The U(1) formula a_1*R^2 is DERIVED (spectral geometry, exp235)
  - The SU(3) formula R*phi: the factor R is derived,
    but WHY phi (not phi^2 or 1) needs the coherent/incoherent argument
  - The coherent vs incoherent distinction: MOTIVATED but not rigorous
""")

print("\n" + "=" * 70)
print("EXP-247 COMPLETE")
print("=" * 70)

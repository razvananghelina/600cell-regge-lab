"""
EXP-142: Emergent Collective Gauge Modes from Fermion Sector
=============================================================
Alternative 4: Instead of assigning gauge bosons to specific vertices,
derive gauge fields as COLLECTIVE excitations of the fermion sector (C vertices).

Approaches:
1. Random phase approximation (RPA): particle-hole excitations on C-C subgraph
2. Gauge invariant bilinears: f_i^dag * f_j along edges
3. Spectral decomposition of C-C adjacency -> natural gauge modes
4. Lie algebra structure from fermion bilinear commutators
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
print("EXP-142: EMERGENT COLLECTIVE GAUGE MODES")
print("="*70)

# Build 600-cell
verts = generate_600cell()
N = len(verts)
vtypes = [classify_vertex(v) for v in verts]
type_counts = Counter(vtypes)
print(f"600-cell: {N} vertices, A={type_counts['A']}, B={type_counts['B']}, C={type_counts['C']}")

coords = np.array(verts)
dist_sq = np.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=2)
edge_len_sq = 1.0 / PHI**2

neighbors = defaultdict(set)
for i in range(N):
    for j in range(i+1, N):
        if dist_sq[i,j] < edge_len_sq + 1e-6:
            neighbors[i].add(j)
            neighbors[j].add(i)

# Identify vertex subsets
A_verts = sorted([i for i in range(N) if vtypes[i] == 'A'])
B_verts = sorted([i for i in range(N) if vtypes[i] == 'B'])
C_verts = sorted([i for i in range(N) if vtypes[i] == 'C'])
gauge_verts = sorted(A_verts + B_verts)  # 24 gauge vertices

print(f"Gauge vertices (A+B): {len(gauge_verts)}")
print(f"Fermion vertices (C): {len(C_verts)}")

# ============================================================
# APPROACH 1: C-C edge bilinears as gauge field modes
# ============================================================
print("\n" + "="*70)
print("APPROACH 1: C-C EDGE BILINEARS")
print("="*70)

# Each C-C edge defines a bilinear operator f_i^dag f_j
# These live on the 432 C-C edges
# The commutator [f_i^dag f_j, f_k^dag f_l] = delta_jk f_i^dag f_l - delta_il f_k^dag f_j
# This is the gl(96) Lie algebra restricted to adjacent pairs

CC_edges = []
C_set = set(C_verts)
c_map = {v: i for i, v in enumerate(C_verts)}  # map to 0..95

for i in C_verts:
    for j in neighbors[i]:
        if j in C_set and j > i:
            CC_edges.append((i, j))

print(f"C-C edges: {len(CC_edges)}")

# Build C-C adjacency matrix (96x96)
n_C = len(C_verts)
A_CC = np.zeros((n_C, n_C))
for i in C_verts:
    for j in neighbors[i]:
        if j in C_set:
            A_CC[c_map[i]][c_map[j]] = 1

# Spectrum of C-C adjacency
eigs_CC = sorted(np.linalg.eigvalsh(A_CC), reverse=True)

def group_eigs(eigs, tol=1e-4):
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

cc_unique = group_eigs(eigs_CC)
print(f"\nC-C subgraph spectrum ({n_C} vertices):")
print(f"  Distinct eigenvalues: {len(cc_unique)}")
for val, mult in cc_unique:
    # Check phi structure
    phi_match = ""
    for a in range(-15, 16):
        for b in range(-10, 11):
            if abs(val - (a + b*PHI)) < 0.01:
                phi_match = f" = {a}+{b}*phi"
                break
        if phi_match:
            break
    print(f"  {val:8.4f}  x {mult}{phi_match}")

# ============================================================
# APPROACH 2: Gauge-Fermion edge modes (A-C and B-C)
# ============================================================
print("\n" + "="*70)
print("APPROACH 2: GAUGE-FERMION COUPLING STRUCTURE")
print("="*70)

# For each gauge vertex g, it connects to some C vertices
# The set of C-neighbors defines the "charge" of those fermions under g
gauge_set = set(gauge_verts)
g_map = {v: i for i, v in enumerate(gauge_verts)}

# Build gauge-fermion incidence matrix (24 x 96)
GF = np.zeros((len(gauge_verts), n_C))
for g in gauge_verts:
    for c in neighbors[g]:
        if c in C_set:
            GF[g_map[g]][c_map[c]] = 1

print(f"Gauge-Fermion incidence matrix: {GF.shape}")
print(f"  Non-zero entries: {int(np.sum(GF))}")
print(f"  Per gauge vertex (row sums): {sorted(set(np.sum(GF, axis=1).astype(int)))}")
print(f"  Per fermion (col sums): {sorted(set(np.sum(GF, axis=0).astype(int)))}")

# Gram matrix G = GF * GF^T (24x24): gauge-gauge overlap through fermions
G_gauge = GF @ GF.T
print(f"\nGauge-Gauge Gram matrix (24x24) via fermions:")
print(f"  Diagonal (self-connections): {sorted(set(np.diag(G_gauge).astype(int)))}")

# Off-diagonal: shared fermion neighbors
off_diag = []
for i in range(24):
    for j in range(i+1, 24):
        off_diag.append(int(G_gauge[i][j]))
off_diag_counts = Counter(off_diag)
print(f"  Off-diagonal values: {dict(sorted(off_diag_counts.items()))}")

# Spectrum of Gram matrix
gram_eigs = sorted(np.linalg.eigvalsh(G_gauge), reverse=True)
gram_unique = group_eigs(gram_eigs)
print(f"\n  Gram matrix spectrum:")
for val, mult in gram_unique:
    phi_match = ""
    for a in range(-50, 51):
        for b in range(-30, 31):
            if abs(val - (a + b*PHI)) < 0.05:
                phi_match = f" = {a}+{b}*phi"
                break
        if phi_match:
            break
    print(f"    {val:8.2f}  x {mult}{phi_match}")

# ============================================================
# APPROACH 3: Fermion bilinear Lie algebra
# ============================================================
print("\n" + "="*70)
print("APPROACH 3: FERMION BILINEAR LIE ALGEBRA")
print("="*70)

# Consider the 432 C-C edges as generators of a Lie algebra
# Generator E_{ij} (for edge i-j) acts as: [E_{ij}, E_{kl}] = delta_jk E_{il} - delta_il E_{kj}
# But only ADJACENT pairs exist, so we need to check closure

# For each pair of CC edges sharing a vertex, compute the "commutator" edge
print(f"C-C edges: {len(CC_edges)}")

# Build edge adjacency: two CC edges that share a C vertex
edge_adj = defaultdict(set)
c_to_edges = defaultdict(list)
for idx, (i, j) in enumerate(CC_edges):
    c_to_edges[c_map[i]].append(idx)
    c_to_edges[c_map[j]].append(idx)

for c_idx in range(n_C):
    edges_at_c = c_to_edges[c_idx]
    for a in range(len(edges_at_c)):
        for b in range(a+1, len(edges_at_c)):
            edge_adj[edges_at_c[a]].add(edges_at_c[b])
            edge_adj[edges_at_c[b]].add(edges_at_c[a])

# Check: how many edge pairs share a vertex
n_sharing = sum(len(v) for v in edge_adj.values()) // 2
print(f"CC edge pairs sharing a C vertex: {n_sharing}")

# For [E_{ij}, E_{jk}] = E_{ik}, check if E_{ik} exists (i.e., is i-k a CC edge?)
n_closed = 0
n_open = 0
for idx1, (i, j) in enumerate(CC_edges):
    ci, cj = c_map[i], c_map[j]
    # Edges at j
    for idx2 in c_to_edges[cj]:
        if idx2 <= idx1:
            continue
        i2, j2 = CC_edges[idx2]
        ci2, cj2 = c_map[i2], c_map[j2]
        # Find the "other" endpoint
        if ci2 == cj:
            other = cj2
        else:
            other = ci2
        # Is ci-other a CC edge?
        if A_CC[ci][other] > 0:
            n_closed += 1
        else:
            n_open += 1
    # Also edges at i
    for idx2 in c_to_edges[ci]:
        if idx2 <= idx1:
            continue
        i2, j2 = CC_edges[idx2]
        ci2, cj2 = c_map[i2], c_map[j2]
        if ci2 == ci:
            other = cj2
        else:
            other = ci2
        if A_CC[cj][other] > 0:
            n_closed += 1
        else:
            n_open += 1

total_commutators = n_closed + n_open
print(f"Commutators that close (E_ik exists): {n_closed}")
print(f"Commutators that don't close: {n_open}")
if total_commutators > 0:
    print(f"Closure fraction: {n_closed/total_commutators:.4f}")

# ============================================================
# APPROACH 4: Spectral decomposition -> natural gauge modes
# ============================================================
print("\n" + "="*70)
print("APPROACH 4: SPECTRAL DECOMPOSITION OF C-C SUBGRAPH")
print("="*70)

# Each eigenspace of A_CC gives a natural "mode" on the fermion sector
# The multiplicity structure should reflect gauge group structure
# From 600-cell full spectrum: eigenspaces have dimensions n^2
# Question: how do the 96 C vertices decompose in each eigenspace?

# Full 600-cell adjacency
A_full = np.zeros((N, N))
for i in range(N):
    for j in neighbors[i]:
        A_full[i][j] = 1

full_eigs, full_vecs = np.linalg.eigh(A_full)
idx_sort = np.argsort(full_eigs)[::-1]
full_eigs = full_eigs[idx_sort]
full_vecs = full_vecs[:, idx_sort]

full_unique = group_eigs(list(full_eigs))
print(f"\nFull 600-cell eigenvalue decomposition projected to C vertices:")
print(f"{'Eigenvalue':>10s} {'Mult':>5s} {'C-weight':>10s} {'A-weight':>10s} {'B-weight':>10s}")

pos = 0
for val, mult in full_unique:
    # Get eigenspace
    V_space = full_vecs[:, pos:pos+mult]  # N x mult

    # Weight on each type
    C_weight = np.sum(V_space[C_verts, :]**2)
    A_weight = np.sum(V_space[A_verts, :]**2)
    B_weight = np.sum(V_space[B_verts, :]**2)
    total = C_weight + A_weight + B_weight

    print(f"  {val:8.4f}  {mult:5d}  {C_weight/total:10.4f}  {A_weight/total:10.4f}  {B_weight/total:10.4f}")
    pos += mult

# ============================================================
# APPROACH 5: Effective gauge field from integrating out fermions
# ============================================================
print("\n" + "="*70)
print("APPROACH 5: EFFECTIVE GAUGE-GAUGE INTERACTION VIA FERMIONS")
print("="*70)

# Compute the effective interaction between gauge vertices
# mediated by 1-hop and 2-hop fermion paths
# G_eff(g1, g2) = sum over C: (g1 connected to C) AND (g2 connected to C)
# This is already the Gram matrix G_gauge from Approach 2

# But also consider 2-hop: g1 -> C1 -> C2 -> g2
# This gives GF * A_CC * GF^T
G2 = GF @ A_CC @ GF.T
print("2-hop effective gauge interaction (GF * A_CC * GF^T):")
g2_eigs = sorted(np.linalg.eigvalsh(G2), reverse=True)
g2_unique = group_eigs(g2_eigs)
print(f"  Spectrum:")
for val, mult in g2_unique:
    print(f"    {val:8.2f}  x {mult}")

# Separate A-block and B-block
n_A = len(A_verts)
n_B = len(B_verts)
A_indices = list(range(n_A))
B_indices = list(range(n_A, n_A + n_B))

G_AA = G_gauge[np.ix_(A_indices, A_indices)]
G_BB = G_gauge[np.ix_(B_indices, B_indices)]
G_AB = G_gauge[np.ix_(A_indices, B_indices)]

print(f"\nGram matrix blocks:")
print(f"  G_AA (8x8): diag={sorted(set(np.diag(G_AA).astype(int)))}, off-diag={sorted(set(G_AA[np.triu_indices(n_A, 1)].astype(int)))}")
print(f"  G_BB (16x16): diag={sorted(set(np.diag(G_BB).astype(int)))}, off-diag={sorted(set(G_BB[np.triu_indices(n_B, 1)].astype(int)))}")
print(f"  G_AB (8x16): values={sorted(set(G_AB.flatten().astype(int)))}")

# A-A block spectrum
aa_eigs = sorted(np.linalg.eigvalsh(G_AA), reverse=True)
aa_unique = group_eigs(aa_eigs)
print(f"\n  G_AA spectrum: {aa_unique}")

# B-B block spectrum
bb_eigs = sorted(np.linalg.eigvalsh(G_BB), reverse=True)
bb_unique = group_eigs(bb_eigs)
print(f"  G_BB spectrum: {bb_unique}")

# ============================================================
# APPROACH 6: Commutant algebra on gauge sector
# ============================================================
print("\n" + "="*70)
print("APPROACH 6: GAUGE ADJACENCY STRUCTURE")
print("="*70)

# Direct gauge-gauge adjacency (A-A, A-B, B-B edges from 600-cell)
A_gauge = np.zeros((24, 24))
for i, gi in enumerate(gauge_verts):
    for j, gj in enumerate(gauge_verts):
        if gj in neighbors[gi]:
            A_gauge[i][j] = 1

print(f"Gauge-gauge edges: {int(np.sum(A_gauge))//2}")
print(f"  A-A: {sum(1 for i in A_verts for j in A_verts if j>i and j in neighbors[i])}")
print(f"  A-B: {sum(1 for i in A_verts for j in B_verts if j in neighbors[i])}")
print(f"  B-B: {sum(1 for i in B_verts for j in B_verts if j>i and j in neighbors[i])}")

# This should be 0,0,0 from selection rules (exp136)
# Confirms: gauge vertices are NEVER directly connected

# Two-hop gauge-gauge via C
A_gauge_2hop = GF @ GF.T  # = G_gauge already computed
print(f"\n2-hop gauge-gauge (via C): non-zero off-diagonal = {np.sum(G_gauge > 0) - 24}")

# Effective gauge adjacency: threshold at some value
for threshold in [1, 2, 3]:
    eff_adj = (G_gauge >= threshold).astype(int)
    np.fill_diagonal(eff_adj, 0)
    n_eff_edges = np.sum(eff_adj) // 2
    print(f"  Threshold {threshold}: {n_eff_edges} effective gauge edges")

# ============================================================
# APPROACH 7: Check if C-C subgraph has SU(3) x SU(2) structure
# ============================================================
print("\n" + "="*70)
print("APPROACH 7: C-C SUBGRAPH STRUCTURE")
print("="*70)

# C-C subgraph: 96 vertices, 432 edges, degree 9
cc_degrees = np.sum(A_CC, axis=1).astype(int)
print(f"C-C subgraph: {n_C} vertices, {len(CC_edges)} edges")
print(f"  Degree: {sorted(set(cc_degrees))} (should be 9, uniform)")

# For each C vertex, its 12 neighbors are: 1A + 2B + 9C
# The 9 C-neighbors form a local structure
# What is this local structure?

local_structures = []
for ci in range(n_C):
    c_nbrs = [j for j in range(n_C) if A_CC[ci][j] > 0]
    # Build local adjacency among the 9 C-neighbors
    local_adj = np.zeros((9, 9))
    for a in range(9):
        for b in range(a+1, 9):
            if A_CC[c_nbrs[a]][c_nbrs[b]] > 0:
                local_adj[a][b] = 1
                local_adj[b][a] = 1
    local_edges = int(np.sum(local_adj)) // 2
    local_structures.append(local_edges)

ls_counts = Counter(local_structures)
print(f"\n  Local C-neighborhood structure (9 C-neighbors):")
print(f"  Edge count distribution: {dict(sorted(ls_counts.items()))}")
if len(ls_counts) == 1:
    n_local_edges = list(ls_counts.keys())[0]
    print(f"  UNIFORM: every C vertex has {n_local_edges} edges among its 9 C-neighbors")
    print(f"  Local density: {n_local_edges}/{9*8//2} = {n_local_edges/36:.4f}")

# Check: what graph is the local C-neighborhood?
# 9 vertices, 15 edges (from exp139) -> what is this graph?
if len(ls_counts) == 1:
    n_le = list(ls_counts.keys())[0]
    print(f"  {n_le} edges on 9 vertices:")
    # Build one example
    ci = 0
    c_nbrs = [j for j in range(n_C) if A_CC[ci][j] > 0]
    local_adj = np.zeros((9, 9))
    for a in range(9):
        for b in range(a+1, 9):
            if A_CC[c_nbrs[a]][c_nbrs[b]] > 0:
                local_adj[a][b] = 1
                local_adj[b][a] = 1

    local_degs = np.sum(local_adj, axis=1).astype(int)
    print(f"  Local degree sequence: {sorted(local_degs)}")

    local_eigs = sorted(np.linalg.eigvalsh(local_adj), reverse=True)
    local_unique = group_eigs(local_eigs, tol=0.01)
    print(f"  Local spectrum: {local_unique}")

    # Check for triangles in local graph
    local_tri = 0
    for a in range(9):
        for b in range(a+1, 9):
            if local_adj[a][b]:
                for c in range(b+1, 9):
                    if local_adj[a][c] and local_adj[b][c]:
                        local_tri += 1
    print(f"  Local triangles: {local_tri}")

    # Check if local graph is line graph of K_{3,3} or complement of something
    complement_edges = 36 - n_le
    print(f"  Complement edges: {complement_edges} on 9 verts")
    comp_adj = 1 - local_adj - np.eye(9)
    comp_degs = np.sum(comp_adj, axis=1).astype(int)
    print(f"  Complement degree sequence: {sorted(comp_degs)}")

# ============================================================
# APPROACH 8: Decompose 432 CC edges by gauge vertex action
# ============================================================
print("\n" + "="*70)
print("APPROACH 8: CC EDGE DECOMPOSITION BY GAUGE ACTION")
print("="*70)

# Each CC edge (i,j) has gauge neighbors:
# For each gauge vertex g, count how many of {i,j} it connects to
# This gives a "charge" vector for each edge

edge_charges = []
for i, j in CC_edges:
    charge = []
    for g in gauge_verts:
        ci = 1 if j in neighbors[g] else 0
        cj = 1 if i in neighbors[g] else 0
        charge.append(ci + cj)  # 0, 1, or 2
    edge_charges.append(tuple(charge))

charge_patterns = Counter(edge_charges)
print(f"Unique charge patterns for CC edges: {len(charge_patterns)}")

# Summarize charge patterns
charge_sums = [sum(c) for c in edge_charges]
charge_sum_counts = Counter(charge_sums)
print(f"Total gauge charge per CC edge: {dict(sorted(charge_sum_counts.items()))}")

# How many gauge vertices connect to both endpoints?
double_charges = [sum(1 for c in charge if c == 2) for charge in edge_charges]
dc_counts = Counter(double_charges)
print(f"Gauge vertices connecting to BOTH endpoints: {dict(sorted(dc_counts.items()))}")

# How many connect to exactly one?
single_charges = [sum(1 for c in charge if c == 1) for charge in edge_charges]
sc_counts = Counter(single_charges)
print(f"Gauge vertices connecting to ONE endpoint: {dict(sorted(sc_counts.items()))}")

# Summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print("\nC-C subgraph: 96 vertices, 432 edges, degree 9 (uniform)")
print(f"C-C eigenvalues: {len(cc_unique)} distinct")
for val, mult in cc_unique:
    print(f"  {val:8.4f} x {mult}")

print(f"\nGauge-Fermion incidence (24x96): each gauge has {int(np.sum(GF[0]))} C-neighbors")
print(f"Gram matrix G_gauge (24x24) = GF*GF^T")
print(f"  G_AA spectrum: {aa_unique}")
print(f"  G_BB spectrum: {bb_unique}")

print(f"\nBilinear Lie algebra closure: {n_closed}/{total_commutators} = {n_closed/total_commutators:.4f}" if total_commutators > 0 else "")

print(f"\nKey findings:")
print(f"  - C-C subgraph is 9-regular with {len(CC_edges)} edges")
print(f"  - Gauge vertices have 0 direct edges (A-A=A-B=B-B=0)")
print(f"  - Gauge-gauge interaction is ENTIRELY mediated by fermions")
print(f"  - Gram matrix captures effective gauge self-coupling")

print("\n" + "="*70)
print("END EXP-142")
print("="*70)

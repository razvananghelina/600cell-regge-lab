#!/usr/bin/env python3
"""
EXP-203: BLOCK-SPIN RENORMALIZATION GROUP ON THE 600-CELL
===========================================================
The 600-cell has 120 vertices: 8 type A, 16 type B, 96 type C.
The 24 A+B vertices form a 24-cell (D4 root system = gauge sector).
The 96 C vertices are fermions (96 = 16 x 3 x 2).

On the bare graph, vertex-transitivity kills differential running.
But there is a NATURAL coarse-graining: the 5 inscribed 24-cells form
a hierarchy:
  Level 0: full 600-cell (120 vertices, 720 edges)
  Level 1: one 24-cell (24 vertices, 96 edges) = gauge sector
  Level 2: one tesseract inside the 24-cell (8 vertices from type A)

This gives a natural "renormalization group" flow by integrating out
fermion vertices.

EXPERIMENT PARTS:
  A. Build 600-cell, classify vertices, classify edges
  B. Effective gauge graph via decimation of C vertices (path counting)
  C. Schur complement: L_eff = L_AB - L_AC * pinv(L_CC) * L_CA
  D. Soliton background: repeat Schur complement with modified weights
  E. RG flow: spectral gap and effective couplings at each scale
  F. 50 soliton backgrounds: statistics on L_eff spectrum
  G. 5 inscribed 24-cells: construction and analysis

STATUS: EXPLORATORY
Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
import math
from itertools import permutations as iperms
from collections import deque, Counter

# ===========================================================================
# CONSTANTS
# ===========================================================================
PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

print("=" * 72)
print("EXP-203: BLOCK-SPIN RENORMALIZATION GROUP ON THE 600-CELL")
print("=" * 72)
print()

# ===========================================================================
# PART A: BUILD 600-CELL, CLASSIFY VERTICES AND EDGES
# ===========================================================================
print("=" * 72)
print("PART A: BUILD 600-CELL AND CLASSIFY")
print("=" * 72)
print()

def build_600cell_vertices():
    """Build the 120 vertices of the 600-cell on the unit 3-sphere."""
    verts = set()
    # Type A: 8 vertices - permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts.add(tuple(v))
    # Type B: 16 vertices - (+-1/2, +-1/2, +-1/2, +-1/2)
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    verts.add((s0, s1, s2, s3))
    # Type C: 96 vertices - even perms of (0, +-1/2, +-phi/2, +-1/(2*phi))
    half = 0.5
    phi_half = PHI / 2.0
    iphi_half = 1.0 / (2.0 * PHI)
    base = [0.0, half, phi_half, iphi_half]
    even_perms = []
    for p in iperms(range(4)):
        inv = 0
        for a in range(4):
            for b in range(a + 1, 4):
                if p[a] > p[b]:
                    inv += 1
        if inv % 2 == 0:
            even_perms.append(p)
    for p in even_perms:
        coords_unsigned = [base[p[0]], base[p[1]], base[p[2]], base[p[3]]]
        nonzero_positions = [i for i in range(4) if abs(coords_unsigned[i]) > 1e-12]
        n_nonzero = len(nonzero_positions)
        for sign_bits in range(1 << n_nonzero):
            v = list(coords_unsigned)
            for k, pos in enumerate(nonzero_positions):
                if (sign_bits >> k) & 1:
                    v[pos] = -v[pos]
            verts.add(tuple(round(x, 12) for x in v))
    return sorted(verts)

vertices = build_600cell_vertices()
N = len(vertices)
print("  Vertices: %d" % N)

# Classify vertex types
def classify_vertex(v):
    coords = [abs(x) for x in v]
    coords_sorted = sorted(coords, reverse=True)
    if abs(coords_sorted[0] - 1.0) < 1e-10 and coords_sorted[1] < 1e-10:
        return 'A'
    if all(abs(c - 0.5) < 1e-10 for c in coords):
        return 'B'
    return 'C'

vtypes = [classify_vertex(v) for v in vertices]
nA = vtypes.count('A')
nB = vtypes.count('B')
nC = vtypes.count('C')
print("  Type A: %d, Type B: %d, Type C: %d" % (nA, nB, nC))

# Index sets
idx_A = [i for i in range(N) if vtypes[i] == 'A']
idx_B = [i for i in range(N) if vtypes[i] == 'B']
idx_C = [i for i in range(N) if vtypes[i] == 'C']
idx_AB = sorted(idx_A + idx_B)
n_AB = len(idx_AB)
n_C = len(idx_C)
print("  A+B vertices (gauge sector): %d" % n_AB)
print("  C vertices (fermion sector): %d" % n_C)

# Build adjacency
varray = np.array(vertices)
dot_matrix = varray @ varray.T
threshold = PHI / 2.0
adj = (np.abs(dot_matrix - threshold) < 1e-8).astype(np.float64)
np.fill_diagonal(adj, 0)

degrees = adj.sum(axis=1).astype(int)
total_edges = int(adj.sum()) // 2
print("  Degree: %d (all vertices), Total edges: %d" % (degrees[0], total_edges))

# Build neighbor lists
neighbors = {}
for i in range(N):
    neighbors[i] = [j for j in range(N) if adj[i, j] > 0.5]

# Edge list
all_edges = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            all_edges.append((i, j))
print("  Edge list: %d edges" % len(all_edges))

# BFS all-pairs shortest paths
print("  Computing all-pairs BFS distances...")
dist_matrix = np.full((N, N), -1, dtype=int)
for start in range(N):
    queue = deque([start])
    dist_matrix[start, start] = 0
    while queue:
        v = queue.popleft()
        for w in neighbors[v]:
            if dist_matrix[start, w] == -1:
                dist_matrix[start, w] = dist_matrix[start, v] + 1
                queue.append(w)
print("  Diameter: %d" % dist_matrix.max())

# Sets for classification
A_set = set(idx_A)
B_set = set(idx_B)
C_set = set(idx_C)
nbr_sets = {i: set(neighbors[i]) for i in range(N)}

# Classify edges by gauge type
def classify_edge_gauge(i, j):
    ti, tj = vtypes[i], vtypes[j]
    if (ti == 'A' and tj == 'C') or (ti == 'C' and tj == 'A'):
        return 'U1'
    if (ti == 'B' and tj == 'C') or (ti == 'C' and tj == 'B'):
        return 'SU2'
    if ti == 'C' and tj == 'C':
        common = nbr_sets[i] & nbr_sets[j]
        shared_A = len(common & A_set)
        shared_B = len(common & B_set)
        if shared_A > 0 and shared_B == 0:
            return 'SU2'
        return 'SU3'
    # A-A, A-B, B-B edges (within 24-cell)
    if ti in ('A', 'B') and tj in ('A', 'B'):
        return 'gauge_internal'
    return 'SU3'

edge_gauge = {}
gauge_counts = Counter()
for i, j in all_edges:
    g = classify_edge_gauge(i, j)
    edge_gauge[(i, j)] = g
    edge_gauge[(j, i)] = g
    gauge_counts[g] += 1

print()
print("  Edge gauge classification:")
for g in ['U1', 'SU2', 'SU3', 'gauge_internal']:
    if gauge_counts[g] > 0:
        print("    %s: %d edges" % (g, gauge_counts[g]))

# Classify edges by vertex-type pairs
pair_counts = Counter()
for i, j in all_edges:
    ti, tj = vtypes[i], vtypes[j]
    pair = tuple(sorted([ti, tj]))
    pair_counts[pair] += 1
print()
print("  Edge type-pair counts:")
for pair, cnt in sorted(pair_counts.items()):
    print("    %s-%s: %d edges" % (pair[0], pair[1], cnt))

# Subgraph statistics
# AB subgraph
ab_edges = [(i, j) for i, j in all_edges if vtypes[i] in ('A', 'B') and vtypes[j] in ('A', 'B')]
print()
print("  24-cell (A+B) subgraph: %d vertices, %d edges" % (n_AB, len(ab_edges)))
# Check degrees within AB subgraph
ab_adj = adj[np.ix_(idx_AB, idx_AB)]
ab_degrees = ab_adj.sum(axis=1).astype(int)
print("  AB subgraph degrees: min=%d, max=%d, mean=%.1f" % (
    ab_degrees.min(), ab_degrees.max(), ab_degrees.mean()))

# AC and BC edges
ac_edges = [(i, j) for i, j in all_edges if
            (vtypes[i] == 'A' and vtypes[j] == 'C') or
            (vtypes[i] == 'C' and vtypes[j] == 'A')]
bc_edges = [(i, j) for i, j in all_edges if
            (vtypes[i] == 'B' and vtypes[j] == 'C') or
            (vtypes[i] == 'C' and vtypes[j] == 'B')]
cc_edges = [(i, j) for i, j in all_edges if vtypes[i] == 'C' and vtypes[j] == 'C']
print("  AC edges: %d, BC edges: %d, CC edges: %d" % (len(ac_edges), len(bc_edges), len(cc_edges)))
print("  Total check: %d + %d + %d + %d = %d (should be %d)" % (
    len(ab_edges), len(ac_edges), len(bc_edges), len(cc_edges),
    len(ab_edges) + len(ac_edges) + len(bc_edges) + len(cc_edges),
    total_edges))
print()


# ===========================================================================
# PART B: EFFECTIVE GAUGE GRAPH VIA DECIMATION OF C VERTICES
# ===========================================================================
print("=" * 72)
print("PART B: EFFECTIVE GAUGE GRAPH BY DECIMATING C VERTICES")
print("=" * 72)
print()

# For each pair of AB vertices, count 2-step paths through C
print("STEP B1: 2-step paths through C (shared C-neighbors)")
print("-" * 60)

# Map AB indices to local indices
ab_local = {glob: loc for loc, glob in enumerate(idx_AB)}

# J_eff[i_local, j_local] = number of shared C-neighbors
J_eff_2step = np.zeros((n_AB, n_AB))
for i_loc, i_glob in enumerate(idx_AB):
    c_nbrs_i = nbr_sets[i_glob] & C_set
    for j_loc, j_glob in enumerate(idx_AB):
        if j_loc <= i_loc:
            continue
        c_nbrs_j = nbr_sets[j_glob] & C_set
        shared_c = len(c_nbrs_i & c_nbrs_j)
        J_eff_2step[i_loc, j_loc] = shared_c
        J_eff_2step[j_loc, i_loc] = shared_c

print("  J_eff (2-step) matrix statistics:")
print("    Max = %d" % int(J_eff_2step.max()))
print("    Non-zero entries: %d (out of %d pairs)" % (
    int(np.count_nonzero(J_eff_2step)) // 2, n_AB * (n_AB - 1) // 2))

# How many pairs have direct AB edge vs path through C?
direct_ab = np.zeros((n_AB, n_AB))
for i_loc, i_glob in enumerate(idx_AB):
    for j_loc, j_glob in enumerate(idx_AB):
        if j_loc > i_loc and adj[i_glob, j_glob] > 0.5:
            direct_ab[i_loc, j_loc] = 1
            direct_ab[j_loc, i_loc] = 1

n_direct = int(direct_ab.sum()) // 2
n_mediated = int((J_eff_2step > 0).sum()) // 2
print("    Direct AB edges: %d" % n_direct)
print("    C-mediated connections (2-step): %d" % n_mediated)

# Distribution of 2-step coupling strengths
vals_2step = []
for i_loc in range(n_AB):
    for j_loc in range(i_loc + 1, n_AB):
        if J_eff_2step[i_loc, j_loc] > 0:
            vals_2step.append(int(J_eff_2step[i_loc, j_loc]))
if vals_2step:
    c_2step = Counter(vals_2step)
    print("    Distribution of shared-C-neighbor counts:")
    for k in sorted(c_2step.keys()):
        print("      %d shared C-nbrs: %d pairs" % (k, c_2step[k]))

# Break down by A-A, A-B, B-B type
print()
print("  2-step coupling by vertex-type pair:")
for pair_label, pair_list in [("A-A", [(i, j) for i in idx_A for j in idx_A if j > i]),
                               ("A-B", [(i, j) for i in idx_A for j in idx_B]),
                               ("B-B", [(i, j) for i in idx_B for j in idx_B if j > i])]:
    vals = []
    for i_glob, j_glob in pair_list:
        i_loc = ab_local[i_glob]
        j_loc = ab_local[j_glob]
        v = J_eff_2step[i_loc, j_loc]
        if v > 0:
            vals.append(int(v))
    if vals:
        print("    %s: %d connected pairs, mean=%.2f, values=%s" % (
            pair_label, len(vals), np.mean(vals), sorted(set(vals))))
    else:
        print("    %s: no 2-step connections through C" % pair_label)

# 3-step and 4-step paths through C subgraph
print()
print("STEP B2: Higher-order paths through C subgraph (Green's function)")
print("-" * 60)

# Build adjacency restricted to C vertices
c_local = {glob: loc for loc, glob in enumerate(idx_C)}
adj_CC = np.zeros((n_C, n_C))
for i_loc, i_glob in enumerate(idx_C):
    for j_glob in neighbors[i_glob]:
        if vtypes[j_glob] == 'C':
            j_loc = c_local[j_glob]
            adj_CC[i_loc, j_loc] = 1.0

# Build coupling matrices: L_AC (from AB to C) etc.
# L_AC[i_AB_local, j_C_local] = adj[i_AB_global, j_C_global]
L_AC_matrix = np.zeros((n_AB, n_C))
for i_loc, i_glob in enumerate(idx_AB):
    for j_glob in neighbors[i_glob]:
        if vtypes[j_glob] == 'C':
            j_loc = c_local[j_glob]
            L_AC_matrix[i_loc, j_loc] = 1.0

print("  C subgraph: %d vertices" % n_C)
print("  CC edges: %d" % (int(adj_CC.sum()) // 2))
print("  L_AC shape: %s, non-zeros: %d" % (str(L_AC_matrix.shape), int(np.count_nonzero(L_AC_matrix))))

# Powers of adj_CC give path counts
# 3-step: AB -> C -> C -> AB = L_AC @ adj_CC @ L_AC.T
J_eff_3step = L_AC_matrix @ adj_CC @ L_AC_matrix.T
n_connected_3 = int((J_eff_3step > 0).sum()) // 2
print("  3-step (AB->C->C->AB) connections: %d pairs" % n_connected_3)
print("    Max coupling: %d" % int(J_eff_3step.max()))

# 4-step: AB -> C -> C -> C -> AB = L_AC @ adj_CC^2 @ L_AC.T
adj_CC2 = adj_CC @ adj_CC
J_eff_4step = L_AC_matrix @ adj_CC2 @ L_AC_matrix.T
n_connected_4 = int((J_eff_4step > 0).sum()) // 2
print("  4-step (AB->C->C->C->AB) connections: %d pairs" % n_connected_4)
print("    Max coupling: %d" % int(J_eff_4step.max()))
print()


# ===========================================================================
# PART C: SCHUR COMPLEMENT (EXACT EFFECTIVE ACTION)
# ===========================================================================
print("=" * 72)
print("PART C: SCHUR COMPLEMENT -- EXACT EFFECTIVE ACTION")
print("=" * 72)
print()

print("Constructing graph Laplacian and Schur complement...")
print("  L_eff = L_AB - L_AC * pinv(L_CC) * L_CA")
print("-" * 60)

# Full graph Laplacian L = D - A
L_full = np.diag(adj.sum(axis=1)) - adj

# Extract blocks using the index ordering [AB | C]
# Reorder so AB indices come first, then C indices
all_order = idx_AB + idx_C
perm = np.array(all_order)

L_reord = L_full[np.ix_(perm, perm)]

# Blocks
L_AB_block = L_reord[:n_AB, :n_AB]
L_AC_block = L_reord[:n_AB, n_AB:]
L_CA_block = L_reord[n_AB:, :n_AB]
L_CC_block = L_reord[n_AB:, n_AB:]

print("  L_AB block: %s" % str(L_AB_block.shape))
print("  L_AC block: %s" % str(L_AC_block.shape))
print("  L_CC block: %s" % str(L_CC_block.shape))

# Check L_CC properties
evals_CC = np.linalg.eigvalsh(L_CC_block)
print()
print("  L_CC eigenvalues: min=%.6f, max=%.6f" % (evals_CC.min(), evals_CC.max()))
print("  L_CC zero eigenvalues (< 1e-10): %d" % np.sum(np.abs(evals_CC) < 1e-10))
print("  L_CC smallest 5: %s" % str(np.round(evals_CC[:5], 6)))

# Schur complement using pseudoinverse (L_CC may have zero eigenvalue)
print()
print("  Computing pseudoinverse of L_CC (%dx%d)..." % L_CC_block.shape)
L_CC_pinv = np.linalg.pinv(L_CC_block, rcond=1e-10)
print("  Done. Residual |L_CC @ pinv(L_CC) @ L_CC - L_CC| = %.2e" % (
    np.linalg.norm(L_CC_block @ L_CC_pinv @ L_CC_block - L_CC_block)))

# Schur complement
L_eff = L_AB_block - L_AC_block @ L_CC_pinv @ L_CA_block

print()
print("  L_eff (Schur complement) shape: %s" % str(L_eff.shape))

# Eigenvalues
evals_eff = np.linalg.eigvalsh(L_eff)
evals_bare = np.linalg.eigvalsh(L_AB_block)

print()
print("  BARE L_AB eigenvalues (just 24-cell direct edges):")
for k, ev in enumerate(evals_bare):
    print("    [%2d] %.6f" % (k, ev))

print()
print("  EFFECTIVE L_eff eigenvalues (after integrating out C):")
for k, ev in enumerate(evals_eff):
    print("    [%2d] %.6f" % (k, ev))

# Spectral gap comparison
gap_bare = evals_bare[1] if abs(evals_bare[0]) < 1e-10 else evals_bare[0]
gap_eff = evals_eff[1] if abs(evals_eff[0]) < 1e-10 else evals_eff[0]
print()
print("  Spectral gap (bare):      %.6f" % gap_bare)
print("  Spectral gap (effective): %.6f" % gap_eff)
print("  Ratio eff/bare:           %.6f" % (gap_eff / gap_bare if gap_bare > 0 else 0))

# Compare with full graph Laplacian eigenvalues
evals_full = np.linalg.eigvalsh(L_full)
print()
print("  FULL graph Laplacian first 10 eigenvalues:")
for k in range(min(10, len(evals_full))):
    print("    [%2d] %.6f" % (k, evals_full[k]))

# Effective adjacency: the off-diagonal of -L_eff
# Effective coupling between AB vertices
L_eff_offdiag = -L_eff.copy()
np.fill_diagonal(L_eff_offdiag, 0)
print()
print("  Effective coupling (off-diag of -L_eff):")
print("    Max = %.6f" % L_eff_offdiag.max())
print("    Min = %.6f" % L_eff_offdiag.min())
print("    Mean (non-zero) = %.6f" % (
    np.mean(L_eff_offdiag[L_eff_offdiag > 1e-10]) if np.any(L_eff_offdiag > 1e-10) else 0))

# Classify effective couplings by vertex-type pair within AB
print()
print("  Effective coupling by vertex type in AB:")
for pair_label, i_list, j_list in [
    ("A-A", idx_A, idx_A),
    ("A-B", idx_A, idx_B),
    ("B-B", idx_B, idx_B)
]:
    vals = []
    for i_g in i_list:
        for j_g in j_list:
            if i_g == j_g:
                continue
            i_loc = ab_local[i_g]
            j_loc = ab_local[j_g]
            v = L_eff_offdiag[i_loc, j_loc]
            if abs(v) > 1e-10:
                vals.append(v)
    if vals:
        print("    %s: %d couplings, mean=%.6f, std=%.6f, range=[%.6f, %.6f]" % (
            pair_label, len(vals), np.mean(vals), np.std(vals), min(vals), max(vals)))
    else:
        print("    %s: no effective couplings" % pair_label)
print()


# ===========================================================================
# PART D: SOLITON BACKGROUND -- SCHUR COMPLEMENT WITH MODIFIED WEIGHTS
# ===========================================================================
print("=" * 72)
print("PART D: SCHUR COMPLEMENT WITH SOLITON BACKGROUND")
print("=" * 72)
print()

# Proper 5-coloring (ground state)
def find_proper_coloring(nbr_list, n_v, n_colors=5, max_attempts=500):
    for seed in range(max_attempts):
        rng_c = np.random.default_rng(seed)
        order = list(range(n_v))
        rng_c.shuffle(order)
        color = [-1] * n_v
        ok = True
        for v in order:
            used = set()
            for w in nbr_list[v]:
                if color[w] >= 0:
                    used.add(color[w])
            assigned = False
            for c in range(n_colors):
                if c not in used:
                    color[v] = c
                    assigned = True
                    break
            if not assigned:
                ok = False
                break
        if ok:
            return color
    # Fallback: Metropolis annealing
    print("  Greedy failed, using Metropolis annealing...")
    best_color = None
    best_E = 9999
    for trial in range(20):
        rng_a = np.random.default_rng(10000 + trial * 77)
        color = [int(rng_a.integers(0, n_colors)) for _ in range(n_v)]
        T = 2.0
        for step in range(500000):
            T = max(2.0 * math.exp(-step / 20000.0), 0.001)
            v = int(rng_a.integers(0, n_v))
            old_c = color[v]
            new_c = int(rng_a.integers(0, n_colors))
            if new_c == old_c:
                continue
            dE = 0
            for w in nbr_list[v]:
                if color[w] == old_c:
                    dE -= 1
                if color[w] == new_c:
                    dE += 1
            if dE <= 0:
                color[v] = new_c
            elif rng_a.random() < math.exp(-dE / T):
                color[v] = new_c
        E = sum(1 for v2 in range(n_v) for w2 in nbr_list[v2]
                if w2 > v2 and color[v2] == color[w2])
        if E < best_E:
            best_E = E
            best_color = list(color)
        if E == 0:
            return color
    return best_color

print("Finding ground state (proper 5-coloring)...")
ground_coloring = find_proper_coloring(neighbors, N, 5)
n_mono_ground = sum(1 for i, j in all_edges if ground_coloring[i] == ground_coloring[j])
print("  Ground state monochromatic edges: %d" % n_mono_ground)
print()

# Kibble-Zurek cooling
def kibble_zurek_cooling(nbr_list, n_v, n_colors, T0, tau, n_steps, rng):
    coloring = [int(rng.integers(0, n_colors)) for _ in range(n_v)]
    for step in range(n_steps):
        t_frac = step / max(n_steps - 1, 1)
        T = T0 * math.exp(-t_frac * n_steps / tau)
        if T < 1e-10:
            T = 1e-10
        v = int(rng.integers(0, n_v))
        old_c = coloring[v]
        new_c = int(rng.integers(0, n_colors))
        if new_c == old_c:
            continue
        dE = 0
        for j in nbr_list[v]:
            if coloring[j] == old_c:
                dE -= 1
            if coloring[j] == new_c:
                dE += 1
        if dE <= 0:
            coloring[v] = new_c
        else:
            if rng.random() < math.exp(-dE / T):
                coloring[v] = new_c
    return coloring

def count_frustrated_edges(coloring, edge_list):
    return sum(1 for i, j in edge_list if coloring[i] == coloring[j])

# Build Laplacian with soliton-modified weights
def build_soliton_laplacian(coloring, adj_matrix, n_v, epsilon=0.5):
    """Build Laplacian where frustrated edges have weight epsilon (< 1)."""
    W = np.zeros((n_v, n_v))
    for i in range(n_v):
        for j in range(i + 1, n_v):
            if adj_matrix[i, j] > 0.5:
                if coloring[i] == coloring[j]:
                    # Frustrated edge
                    W[i, j] = epsilon
                    W[j, i] = epsilon
                else:
                    W[i, j] = 1.0
                    W[j, i] = 1.0
    L = np.diag(W.sum(axis=1)) - W
    return L, W

def schur_complement_from_laplacian(L, idx_keep, idx_integrate):
    """Compute Schur complement: L_keep - L_ki * pinv(L_ii) * L_ik."""
    n_k = len(idx_keep)
    n_i = len(idx_integrate)
    L_kk = L[np.ix_(idx_keep, idx_keep)]
    L_ki = L[np.ix_(idx_keep, idx_integrate)]
    L_ik = L[np.ix_(idx_integrate, idx_keep)]
    L_ii = L[np.ix_(idx_integrate, idx_integrate)]
    L_ii_pinv = np.linalg.pinv(L_ii, rcond=1e-10)
    L_eff = L_kk - L_ki @ L_ii_pinv @ L_ik
    return L_eff

# Generate one soliton background
print("Generating one soliton background (Kibble-Zurek)...")
rng_sol = np.random.default_rng(42)
sol_coloring = kibble_zurek_cooling(neighbors, N, 5, 5.0, 5000, 100000, rng_sol)
n_frust = count_frustrated_edges(sol_coloring, all_edges)
print("  Frustrated edges: %d / %d (fraction = %.4f)" % (n_frust, total_edges, n_frust / total_edges))
print()

# Schur complement with soliton background
print("Computing Schur complement with soliton background...")
for eps_val in [0.1, 0.3, 0.5, 0.8]:
    L_sol, W_sol = build_soliton_laplacian(sol_coloring, adj, N, epsilon=eps_val)
    L_eff_sol = schur_complement_from_laplacian(L_sol, idx_AB, idx_C)
    evals_sol = np.linalg.eigvalsh(L_eff_sol)

    print("  epsilon = %.1f:" % eps_val)
    print("    L_eff eigenvalues (first 8): %s" % str(np.round(evals_sol[:8], 4)))
    gap_sol = evals_sol[1] if abs(evals_sol[0]) < 1e-8 else evals_sol[0]
    print("    Spectral gap: %.6f (bare: %.6f, ratio: %.4f)" % (
        gap_sol, gap_bare, gap_sol / gap_bare if gap_bare > 0 else 0))

# Compare soliton vs ground state
print()
print("  Ground state Schur complement (epsilon=0.5, all edges weight=1):")
L_eff_ground = schur_complement_from_laplacian(L_full, idx_AB, idx_C)
evals_ground = np.linalg.eigvalsh(L_eff_ground)
print("    Eigenvalues: %s" % str(np.round(evals_ground[:8], 4)))

L_sol_05, _ = build_soliton_laplacian(sol_coloring, adj, N, epsilon=0.5)
L_eff_sol_05 = schur_complement_from_laplacian(L_sol_05, idx_AB, idx_C)
evals_sol_05 = np.linalg.eigvalsh(L_eff_sol_05)
print("  Soliton Schur complement (epsilon=0.5):")
print("    Eigenvalues: %s" % str(np.round(evals_sol_05[:8], 4)))
print("    Diff from ground: %s" % str(np.round(evals_sol_05[:8] - evals_ground[:8], 4)))

# Per-gauge-sector analysis
print()
print("  Effective coupling per gauge channel (soliton, eps=0.5):")
L_eff_sol_offdiag = -L_eff_sol_05.copy()
np.fill_diagonal(L_eff_sol_offdiag, 0)

for pair_label, i_list, j_list in [
    ("A-A", idx_A, idx_A),
    ("A-B", idx_A, idx_B),
    ("B-B", idx_B, idx_B)
]:
    vals = []
    for i_g in i_list:
        for j_g in j_list:
            if i_g == j_g:
                continue
            i_loc = ab_local[i_g]
            j_loc = ab_local[j_g]
            v = L_eff_sol_offdiag[i_loc, j_loc]
            if abs(v) > 1e-10:
                vals.append(v)
    if vals:
        print("    %s: mean=%.6f, std=%.6f" % (pair_label, np.mean(vals), np.std(vals)))
    else:
        print("    %s: no couplings" % pair_label)
print()


# ===========================================================================
# PART E: RG FLOW -- SPECTRAL GAP AND COUPLINGS AT EACH SCALE
# ===========================================================================
print("=" * 72)
print("PART E: RG FLOW ACROSS SCALES")
print("=" * 72)
print()

print("Three scales:")
print("  Level 0: Full 600-cell (120 vertices)")
print("  Level 1: 24-cell = A+B (24 vertices, C integrated out)")
print("  Level 2: Tesseract = A only (8 vertices, B+C integrated out)")
print("-" * 60)

# Level 0: Full graph
evals_L0 = np.linalg.eigvalsh(L_full)
gap_L0 = evals_L0[1]
print()
print("Level 0 (full 600-cell, 120 vertices):")
print("  Spectral gap: %.6f" % gap_L0)
print("  First 5 eigenvalues: %s" % str(np.round(evals_L0[:5], 6)))
print("  Average edge weight: 1.0 (uniform)")

# Level 1: A+B with C integrated out
L_eff_L1 = schur_complement_from_laplacian(L_full, idx_AB, idx_C)
evals_L1 = np.linalg.eigvalsh(L_eff_L1)
gap_L1 = evals_L1[1] if abs(evals_L1[0]) < 1e-8 else evals_L1[0]
offdiag_L1 = -L_eff_L1.copy()
np.fill_diagonal(offdiag_L1, 0)
avg_weight_L1 = np.mean(offdiag_L1[offdiag_L1 > 1e-10]) if np.any(offdiag_L1 > 1e-10) else 0
print()
print("Level 1 (24-cell, 24 vertices, C integrated out):")
print("  Spectral gap: %.6f" % gap_L1)
print("  First 5 eigenvalues: %s" % str(np.round(evals_L1[:5], 6)))
print("  Average effective edge weight: %.6f" % avg_weight_L1)

# Level 2: A only, with B+C integrated out
idx_BC = sorted(idx_B + idx_C)
L_eff_L2 = schur_complement_from_laplacian(L_full, idx_A, idx_BC)
evals_L2 = np.linalg.eigvalsh(L_eff_L2)
gap_L2 = evals_L2[1] if abs(evals_L2[0]) < 1e-8 else evals_L2[0]
offdiag_L2 = -L_eff_L2.copy()
np.fill_diagonal(offdiag_L2, 0)
avg_weight_L2 = np.mean(offdiag_L2[offdiag_L2 > 1e-10]) if np.any(offdiag_L2 > 1e-10) else 0
print()
print("Level 2 (tesseract, 8 vertices = type A, B+C integrated out):")
print("  Spectral gap: %.6f" % gap_L2)
print("  First 5 eigenvalues: %s" % str(np.round(evals_L2[:5], 6)))
print("  Average effective edge weight: %.6f" % avg_weight_L2)

# Also do Level 2b: A only, with ONLY B integrated out from 24-cell
print()
print("Level 2b: A from 24-cell subgraph (integrating out B from bare 24-cell):")
L_24 = L_AB_block  # This is the 24-cell Laplacian block
# Reorder within 24-cell: first A, then B
# Local indices of A within AB
a_in_ab = [ab_local[g] for g in idx_A]
b_in_ab = [ab_local[g] for g in idx_B]
L_eff_L2b = schur_complement_from_laplacian(L_24, a_in_ab, b_in_ab)
evals_L2b = np.linalg.eigvalsh(L_eff_L2b)
gap_L2b = evals_L2b[1] if abs(evals_L2b[0]) < 1e-8 else evals_L2b[0]
print("  Spectral gap: %.6f" % gap_L2b)
print("  Eigenvalues: %s" % str(np.round(evals_L2b, 6)))

# Summary table
print()
print("=" * 60)
print("RG FLOW SUMMARY")
print("=" * 60)
print("  Scale   | Vertices | Spect. Gap | Avg Weight | Gap Ratio")
print("  --------|----------|------------|------------|----------")
print("  Full    | %6d   | %10.6f | %10.6f | %9.4f" % (N, gap_L0, 1.0, 1.0))
print("  24-cell | %6d   | %10.6f | %10.6f | %9.4f" % (n_AB, gap_L1, avg_weight_L1, gap_L1 / gap_L0))
print("  Tess.   | %6d   | %10.6f | %10.6f | %9.4f" % (len(idx_A), gap_L2, avg_weight_L2, gap_L2 / gap_L0))
print()

# With soliton
print("RG FLOW WITH SOLITON (epsilon=0.5):")
L_sol_full, _ = build_soliton_laplacian(sol_coloring, adj, N, epsilon=0.5)
evals_sol_L0 = np.linalg.eigvalsh(L_sol_full)
gap_sol_L0 = evals_sol_L0[1]

L_eff_sol_L1 = schur_complement_from_laplacian(L_sol_full, idx_AB, idx_C)
evals_sol_L1 = np.linalg.eigvalsh(L_eff_sol_L1)
gap_sol_L1 = evals_sol_L1[1] if abs(evals_sol_L1[0]) < 1e-8 else evals_sol_L1[0]

L_eff_sol_L2 = schur_complement_from_laplacian(L_sol_full, idx_A, idx_BC)
evals_sol_L2 = np.linalg.eigvalsh(L_eff_sol_L2)
gap_sol_L2 = evals_sol_L2[1] if abs(evals_sol_L2[0]) < 1e-8 else evals_sol_L2[0]

print("  Scale   | Spect. Gap | Gap/Ground | Gap Ratio")
print("  --------|------------|------------|----------")
print("  Full    | %10.6f | %10.4f   | %9.4f" % (gap_sol_L0, gap_sol_L0 / gap_L0, 1.0))
print("  24-cell | %10.6f | %10.4f   | %9.4f" % (gap_sol_L1, gap_sol_L1 / gap_L1, gap_sol_L1 / gap_sol_L0))
print("  Tess.   | %10.6f | %10.4f   | %9.4f" % (gap_sol_L2, gap_sol_L2 / gap_L2, gap_sol_L2 / gap_sol_L0))
print()


# ===========================================================================
# PART F: 50 SOLITON BACKGROUNDS -- STATISTICS
# ===========================================================================
print("=" * 72)
print("PART F: 50 SOLITON BACKGROUNDS -- L_eff SPECTRUM STATISTICS")
print("=" * 72)
print()

N_REAL = 50
n_steps_kz = 100000
tau_kz = 5000
T0_kz = 5.0

print("Generating %d Kibble-Zurek realizations (tau=%d, steps=%d)..." % (N_REAL, tau_kz, n_steps_kz))

all_gaps_L1 = []
all_gaps_L2 = []
all_frust_fracs = []
all_evals_L1 = []

for r in range(N_REAL):
    rng_r = np.random.default_rng(2000 + r * 53)
    col_r = kibble_zurek_cooling(neighbors, N, 5, T0_kz, tau_kz, n_steps_kz, rng_r)
    n_fr = count_frustrated_edges(col_r, all_edges)
    frac = n_fr / total_edges
    all_frust_fracs.append(frac)

    L_r, _ = build_soliton_laplacian(col_r, adj, N, epsilon=0.5)
    L_eff_r_L1 = schur_complement_from_laplacian(L_r, idx_AB, idx_C)
    ev_r = np.linalg.eigvalsh(L_eff_r_L1)
    all_evals_L1.append(ev_r)
    g1 = ev_r[1] if abs(ev_r[0]) < 1e-8 else ev_r[0]
    all_gaps_L1.append(g1)

    L_eff_r_L2 = schur_complement_from_laplacian(L_r, idx_A, idx_BC)
    ev_r2 = np.linalg.eigvalsh(L_eff_r_L2)
    g2 = ev_r2[1] if abs(ev_r2[0]) < 1e-8 else ev_r2[0]
    all_gaps_L2.append(g2)

    if (r + 1) % 10 == 0:
        print("  Done %d/%d realizations" % (r + 1, N_REAL))

all_gaps_L1 = np.array(all_gaps_L1)
all_gaps_L2 = np.array(all_gaps_L2)
all_frust_fracs = np.array(all_frust_fracs)
all_evals_L1 = np.array(all_evals_L1)

print()
print("Results over %d realizations:" % N_REAL)
print("  Frustration fraction: %.4f +/- %.4f (range [%.4f, %.4f])" % (
    np.mean(all_frust_fracs), np.std(all_frust_fracs),
    np.min(all_frust_fracs), np.max(all_frust_fracs)))
print()
print("  L_eff spectral gap (24-cell scale):")
print("    Mean: %.6f +/- %.6f" % (np.mean(all_gaps_L1), np.std(all_gaps_L1)))
print("    Ground state value: %.6f" % gap_L1)
print("    Mean/Ground ratio: %.4f" % (np.mean(all_gaps_L1) / gap_L1))
print()
print("  L_eff spectral gap (tesseract scale):")
print("    Mean: %.6f +/- %.6f" % (np.mean(all_gaps_L2), np.std(all_gaps_L2)))
print("    Ground state value: %.6f" % gap_L2)
print("    Mean/Ground ratio: %.4f" % (np.mean(all_gaps_L2) / gap_L2))

# Correlation between frustration fraction and spectral gap
corr_frust_gap1 = np.corrcoef(all_frust_fracs, all_gaps_L1)[0, 1]
corr_frust_gap2 = np.corrcoef(all_frust_fracs, all_gaps_L2)[0, 1]
print()
print("  Correlation (frustration fraction vs gap):")
print("    24-cell: r = %.4f" % corr_frust_gap1)
print("    Tesseract: r = %.4f" % corr_frust_gap2)

# Ratio of gaps across scales
ratios_L2_L1 = all_gaps_L2 / all_gaps_L1
print()
print("  Ratio gap(tess)/gap(24-cell) per realization:")
print("    Mean: %.4f +/- %.4f" % (np.mean(ratios_L2_L1), np.std(ratios_L2_L1)))
print("    Ground state: %.4f" % (gap_L2 / gap_L1))

# Eigenvalue degeneracy breaking
print()
print("  Eigenvalue spread (std across realizations per eigenvalue):")
evals_std = np.std(all_evals_L1, axis=0)
evals_mean = np.mean(all_evals_L1, axis=0)
for k in range(min(n_AB, 10)):
    print("    eval[%2d]: mean=%.4f, std=%.4f, rel_std=%.4f" % (
        k, evals_mean[k], evals_std[k],
        evals_std[k] / abs(evals_mean[k]) if abs(evals_mean[k]) > 1e-10 else 0))
print()


# ===========================================================================
# PART G: 5 INSCRIBED 24-CELLS
# ===========================================================================
print("=" * 72)
print("PART G: THE 5 INSCRIBED 24-CELLS")
print("=" * 72)
print()

print("The 600-cell has 5 inscribed 24-cells, each with 24 vertices.")
print("The proper 5-coloring partitions vertices into 5 cosets of 24.")
print("-" * 60)

# Partition vertices by color
cosets = {c: [] for c in range(5)}
for i in range(N):
    cosets[ground_coloring[i]].append(i)

print()
for c in range(5):
    n_c = len(cosets[c])
    types_c = Counter(vtypes[i] for i in cosets[c])
    print("  Coset %d: %d vertices (A=%d, B=%d, C=%d)" % (
        c, n_c, types_c.get('A', 0), types_c.get('B', 0), types_c.get('C', 0)))

# Each coset should form a 24-cell: check internal edges (should be 0,
# since proper coloring means no monochromatic edges)
print()
print("Checking coset subgraphs (should have 0 internal edges by coloring):")
for c in range(5):
    verts_c = cosets[c]
    n_internal = 0
    for i in verts_c:
        for j in verts_c:
            if j > i and adj[i, j] > 0.5:
                n_internal += 1
    print("  Coset %d: %d internal edges" % (c, n_internal))

# The 5 inscribed 24-cells are NOT the color cosets.
# The cosets are INDEPENDENT SETS (no internal edges).
# The 24-cells are subgraphs with 96 edges each.
# A different approach: use the algebraic structure.
print()
print("NOTE: The cosets from proper 5-coloring are independent sets.")
print("      The inscribed 24-cells are NOT the same as color cosets.")
print("      The 24-cells overlap with the coloring structure differently.")
print()

# The A+B vertices (8+16=24) form ONE inscribed 24-cell.
# Let us verify: does A+B form a 24-cell?
print("Checking A+B subgraph as a 24-cell candidate:")
ab_verts = idx_AB
ab_internal_edges = 0
for i in ab_verts:
    for j in ab_verts:
        if j > i and adj[i, j] > 0.5:
            ab_internal_edges += 1
print("  A+B: %d vertices, %d edges" % (len(ab_verts), ab_internal_edges))
print("  Expected for 24-cell: 24 vertices, 96 edges")

# Eigenvalues of the A+B subgraph
adj_ab_sub = adj[np.ix_(idx_AB, idx_AB)]
L_ab_sub = np.diag(adj_ab_sub.sum(axis=1)) - adj_ab_sub
evals_ab_sub = np.linalg.eigvalsh(L_ab_sub)
print("  A+B subgraph Laplacian eigenvalues:")
# Group by unique values
evals_rounded = np.round(evals_ab_sub, 4)
eval_counts = Counter(evals_rounded)
for ev in sorted(eval_counts.keys()):
    print("    %.4f (multiplicity %d)" % (ev, eval_counts[ev]))

# Check 24-cell eigenvalues (the standard 24-cell has specific spectrum)
print()
print("  Standard 24-cell adjacency eigenvalues should be:")
print("    8 (mult 1), 4 (mult 3), 2 (mult 3), 0 (mult 8), -2 (mult 3), -4 (mult 6)")
adj_ab_evals = np.linalg.eigvalsh(adj_ab_sub)
adj_ab_rounded = np.round(adj_ab_evals, 4)
adj_ab_counts = Counter(adj_ab_rounded)
print("  Actual A+B adjacency eigenvalues:")
for ev in sorted(adj_ab_counts.keys(), reverse=True):
    print("    %.4f (multiplicity %d)" % (ev, adj_ab_counts[ev]))

# Find the other 4 inscribed 24-cells using quaternion multiplication
# Each left coset of the binary tetrahedral group in the binary icosahedral
# group gives a 24-cell.
print()
print("Finding 5 inscribed 24-cells by distance-2 shells from A+B...")
print("  (Alternative approach: find 24-vertex subgraphs isomorphic to 24-cell)")
print()

# Method: The 5 inscribed 24-cells can be found via the H4 -> D4 branching.
# One 24-cell = A+B vertices. For the others, we use the symmetry.
# The proper 5-coloring assigns each vertex to one of 5 cosets.
# The complement approach: take all vertices AT distance 2 from a given vertex.

# Instead, let's check which vertices have the same "norm pattern"
# as A+B but rotated. The 5 inscribed 24-cells are related by phi rotations.

# A simpler check: how many vertices at each graph distance from the A+B set?
dist_from_AB = {}
for v in range(N):
    min_d = min(dist_matrix[v, a] for a in idx_AB if v != a) if v not in set(idx_AB) else 0
    dist_from_AB[v] = min_d

dist_hist = Counter(dist_from_AB.values())
print("  Distance from A+B set:")
for d in sorted(dist_hist.keys()):
    print("    d=%d: %d vertices" % (d, dist_hist[d]))

# For each C vertex, what is its BFS distance to the nearest A vertex and B vertex?
print()
print("  C vertex distances to A and B subsets:")
c_dists_to_A = []
c_dists_to_B = []
for c in idx_C:
    dA = min(dist_matrix[c, a] for a in idx_A)
    dB = min(dist_matrix[c, b] for b in idx_B)
    c_dists_to_A.append(dA)
    c_dists_to_B.append(dB)
cnt_dA = Counter(c_dists_to_A)
cnt_dB = Counter(c_dists_to_B)
print("    dist to A: %s" % dict(sorted(cnt_dA.items())))
print("    dist to B: %s" % dict(sorted(cnt_dB.items())))

# Summary of degree within each sector
print()
print("  Connectivity summary (how many neighbors of each type per vertex):")
for vt in ['A', 'B', 'C']:
    vt_idx = [i for i in range(N) if vtypes[i] == vt]
    nbr_A = [sum(1 for j in neighbors[i] if vtypes[j] == 'A') for i in vt_idx]
    nbr_B = [sum(1 for j in neighbors[i] if vtypes[j] == 'B') for i in vt_idx]
    nbr_C = [sum(1 for j in neighbors[i] if vtypes[j] == 'C') for i in vt_idx]
    print("    Type %s: <nbr_A>=%.1f, <nbr_B>=%.1f, <nbr_C>=%.1f, total=%.1f" % (
        vt, np.mean(nbr_A), np.mean(nbr_B), np.mean(nbr_C),
        np.mean(nbr_A) + np.mean(nbr_B) + np.mean(nbr_C)))


# ===========================================================================
# FINAL SUMMARY
# ===========================================================================
print()
print("=" * 72)
print("FINAL SUMMARY: EXP-203")
print("=" * 72)
print()
print("1. 600-cell structure:")
print("   120 vertices (A=8, B=16, C=96), 720 edges, degree=12")
print("   Edge classification: U1=%d, SU2=%d, SU3=%d, internal=%d" % (
    gauge_counts['U1'], gauge_counts['SU2'], gauge_counts['SU3'], gauge_counts.get('gauge_internal', 0)))
print()
print("2. Decimation of C (fermion) vertices:")
print("   2-step paths through C: %d AB pairs connected" % n_mediated)
print("   3-step paths: %d pairs connected" % n_connected_3)
print()
print("3. Schur complement (exact effective action):")
print("   Bare 24-cell spectral gap:  %.6f" % gap_bare)
print("   Effective (C decimated):    %.6f" % gap_L1)
print("   Effective (B+C decimated):  %.6f" % gap_L2)
print()
print("4. RG flow (ground state):")
print("   Full -> 24-cell gap ratio:  %.4f" % (gap_L1 / gap_L0))
print("   Full -> tesseract gap ratio: %.4f" % (gap_L2 / gap_L0))
print()
print("5. Soliton statistics (%d realizations):" % N_REAL)
print("   Mean frustration fraction: %.4f +/- %.4f" % (
    np.mean(all_frust_fracs), np.std(all_frust_fracs)))
print("   Gap(24-cell) / Ground: %.4f +/- %.4f" % (
    np.mean(all_gaps_L1) / gap_L1, np.std(all_gaps_L1) / gap_L1))
print("   Gap correlation with frustration: r=%.4f (24-cell), r=%.4f (tess)" % (
    corr_frust_gap1, corr_frust_gap2))
print()
print("6. A+B subgraph = 24-cell:")
print("   %d vertices, %d edges" % (n_AB, ab_internal_edges))
print()
print("STATUS: EXPLORATORY")
print("KEY QUESTION: Does Schur complement generate differential running")
print("between gauge sectors even though bare graph is vertex-transitive?")
print()
print("DONE.")

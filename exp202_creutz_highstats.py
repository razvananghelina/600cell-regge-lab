#!/usr/bin/env python3
"""
EXP-202: HIGH-STATISTICS CREUTZ RATIO ON THE 600-CELL
======================================================
EXP-199 found confinement (area law, sigma~0.10) and differential Creutz
ratios on soliton backgrounds, but with only 20 realizations and limited
statistics. This experiment does a HIGH-STATISTICS measurement with 200
realizations, all pairs, and bootstrap error bars.

FRUSTRATION-BASED Wilson line:
  U(i,j) = +1 if satisfied (different colors), -1 if frustrated (same color)
  W_line(path) = product of U(i,j) along shortest path
  W(d) = average of W_line over all pairs at graph distance d
  chi(d) = -ln(|W(d)|/|W(d-1)|)

PARTS:
  A. Build 600-cell, classify edges, BFS distance matrix
  B. Generate 200 Kibble-Zurek backgrounds (~88 defects)
  C. Measure W(d) for d=1..5, classified by gauge type
  D. Bootstrap error bars (1000 resamples)
  E. Test predictions: chi_SU3/chi_U1 vs 10*phi, phi^3
  F. Alternative classification: MIDDLE edge of path

STATUS: EXPLORATORY
Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
import math
from itertools import permutations as iperms
from collections import deque, Counter
import time

# ===========================================================================
# CONSTANTS
# ===========================================================================
PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

print("=" * 72)
print("EXP-202: HIGH-STATISTICS CREUTZ RATIO ON THE 600-CELL")
print("=" * 72)
print()

t0_global = time.time()

# ===========================================================================
# PART A: BUILD 600-CELL, CLASSIFY EDGES, BFS DISTANCES
# ===========================================================================
print("=" * 72)
print("PART A: BUILD 600-CELL AND INFRASTRUCTURE")
print("=" * 72)
print()

# --- Build 600-cell vertices ---
print("STEP A1: Constructing 600-cell vertices")
print("-" * 60)

def build_600cell_vertices():
    """Build the 120 vertices of the 600-cell on the unit 3-sphere."""
    verts = set()
    # Type A: 8 vertices - all permutations of (+-1, 0, 0, 0)
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

# Build adjacency using numpy
varray = np.array(vertices)
dot_matrix = varray @ varray.T
threshold = PHI / 2.0
adj = (np.abs(dot_matrix - threshold) < 1e-8).astype(np.int8)
np.fill_diagonal(adj, 0)

degrees = adj.sum(axis=1)
print("  Degree: min=%d, max=%d" % (degrees.min(), degrees.max()))
total_edges = int(adj.sum()) // 2
print("  Total edges: %d" % total_edges)

# Neighbor lists (as numpy array for speed)
neighbors = {}
for i in range(N):
    neighbors[i] = list(np.where(adj[i] > 0)[0])

# Edge list
all_edges = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            all_edges.append((i, j))
print("  Edge list: %d edges" % len(all_edges))

# --- Classify edges by gauge type ---
print()
print("STEP A2: Edge gauge classification")
print("-" * 60)

A_set = set(i for i in range(N) if vtypes[i] == 'A')
B_set = set(i for i in range(N) if vtypes[i] == 'B')
C_set = set(i for i in range(N) if vtypes[i] == 'C')
nbr_sets = {i: set(neighbors[i]) for i in range(N)}

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
    return 'SU3'

edge_gauge = {}
gauge_counts = Counter()
for i, j in all_edges:
    g = classify_edge_gauge(i, j)
    edge_gauge[(i, j)] = g
    edge_gauge[(j, i)] = g
    gauge_counts[g] += 1

print("  Edge gauge classification:")
for g in ['U1', 'SU2', 'SU3']:
    print("    %s: %d edges" % (g, gauge_counts[g]))

# --- BFS distance matrix + shortest path representative ---
print()
print("STEP A3: Computing all-pairs shortest paths")
print("-" * 60)

dist_matrix = np.full((N, N), -1, dtype=np.int32)
# Also store one representative shortest path for each pair (as next-hop)
# next_hop[i,j] = the next vertex after i on a shortest path from i to j
next_hop = np.full((N, N), -1, dtype=np.int32)

for start in range(N):
    queue = deque([start])
    dist_matrix[start, start] = 0
    while queue:
        v = queue.popleft()
        for w in neighbors[v]:
            if dist_matrix[start, w] == -1:
                dist_matrix[start, w] = dist_matrix[start, v] + 1
                queue.append(w)
                # For next_hop: reconstruct via BFS parent
                # We store next_hop from start's perspective
                if v == start:
                    next_hop[start, w] = w
                else:
                    next_hop[start, w] = next_hop[start, v]

max_dist = dist_matrix.max()
print("  Graph diameter: %d" % max_dist)

# Pre-compute pairs by distance
pairs_by_dist = {d: [] for d in range(1, max_dist + 1)}
for i in range(N):
    for j in range(i + 1, N):
        d_ij = dist_matrix[i, j]
        if 1 <= d_ij <= max_dist:
            pairs_by_dist[d_ij].append((i, j))

for d in range(1, max_dist + 1):
    print("  d=%d: %d pairs" % (d, len(pairs_by_dist[d])))

# Function to reconstruct a shortest path from next_hop table
def get_shortest_path(src, dst):
    """Reconstruct shortest path from src to dst using next_hop table."""
    if src == dst:
        return [src]
    path = [src]
    cur = src
    for _ in range(max_dist + 1):
        nxt = next_hop[cur, dst]
        if nxt < 0:
            return path  # shouldn't happen
        path.append(nxt)
        if nxt == dst:
            return path
        cur = nxt
    return path

# Pre-compute the gauge type of each edge along shortest paths
# For each pair (i,j), store the list of gauge types along the path
print()
print("STEP A4: Pre-computing path edge gauge types for all pairs")
print("-" * 60)

# Store edge gauge types for each pair's shortest path
# Also classify by: dominant, first, middle edge
# Structure: for each distance d, for each pair, store edge_types list

path_edge_types = {}  # (i,j) -> list of gauge types
for d in range(1, max_dist + 1):
    for (i, j) in pairs_by_dist[d]:
        path = get_shortest_path(i, j)
        etypes = []
        for step in range(len(path) - 1):
            vi = path[step]
            vj = path[step + 1]
            gt = edge_gauge.get((vi, vj), 'SU3')
            etypes.append(gt)
        path_edge_types[(i, j)] = etypes

print("  Pre-computed edge types for %d pairs" % len(path_edge_types))

# Classification functions
def dominant_gauge(etypes):
    """Most common gauge type in path."""
    if not etypes:
        return 'unknown'
    cnt = Counter(etypes)
    return cnt.most_common(1)[0][0]

def first_gauge(etypes):
    """Gauge type of first edge in path."""
    if not etypes:
        return 'unknown'
    return etypes[0]

def middle_gauge(etypes):
    """Gauge type of middle edge(s) in path."""
    if not etypes:
        return 'unknown'
    n = len(etypes)
    if n == 1:
        return etypes[0]
    mid = n // 2
    if n % 2 == 1:
        return etypes[mid]
    else:
        # Even length: take both middle edges
        g1 = etypes[mid - 1]
        g2 = etypes[mid]
        if g1 == g2:
            return g1
        # If different, pick the "stronger" one (SU3 > SU2 > U1)
        order = {'SU3': 3, 'SU2': 2, 'U1': 1}
        return g1 if order.get(g1, 0) >= order.get(g2, 0) else g2

print("  Classification methods: dominant, first, middle")
print()


# ===========================================================================
# PART B: GENERATE 200 KIBBLE-ZUREK BACKGROUNDS
# ===========================================================================
print("=" * 72)
print("PART B: GENERATE 200 KIBBLE-ZUREK SOLITON BACKGROUNDS")
print("=" * 72)
print()

def potts_energy(coloring, nbr_list, n_v):
    E = 0
    for i in range(n_v):
        for j in nbr_list[i]:
            if j > i and coloring[i] == coloring[j]:
                E += 1
    return E

def count_defects(coloring, nbr_list, n_v):
    defects = 0
    for i in range(n_v):
        for j in nbr_list[i]:
            if coloring[i] == coloring[j]:
                defects += 1
                break
    return defects

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

# --- Calibrate tau for ~88 defects ---
print("  Calibrating tau for ~88 defects...")
tau_values = [500, 1000, 2000, 5000, 10000]
n_cooling_steps = 100000
best_tau = None
best_diff = 999

for tau in tau_values:
    defects_list = []
    for r in range(5):
        rng = np.random.default_rng(42 + r * 100 + tau)
        col = kibble_zurek_cooling(neighbors, N, 5, 5.0, tau, n_cooling_steps, rng)
        nd = count_defects(col, neighbors, N)
        defects_list.append(nd)
    mean_def = np.mean(defects_list)
    diff = abs(mean_def - 88)
    print("    tau=%5d: <n_defects>=%.1f" % (tau, mean_def))
    if diff < best_diff:
        best_diff = diff
        best_tau = tau

print("  Best tau: %d" % best_tau)

# --- Generate 200 realizations ---
N_REAL = 200
print()
print("  Generating %d independent soliton backgrounds..." % N_REAL)
t_start_gen = time.time()

soliton_backgrounds = []
all_defect_counts = []
all_energies = []

for r in range(N_REAL):
    rng = np.random.default_rng(1000 + r * 37)
    col = kibble_zurek_cooling(neighbors, N, 5, 5.0, best_tau, n_cooling_steps, rng)
    nd = count_defects(col, neighbors, N)
    E = potts_energy(col, neighbors, N)
    soliton_backgrounds.append(col)
    all_defect_counts.append(nd)
    all_energies.append(E)
    if (r + 1) % 50 == 0:
        elapsed = time.time() - t_start_gen
        print("    Generated %d / %d (%.1f sec)" % (r + 1, N_REAL, elapsed))

t_gen = time.time() - t_start_gen
print("  Generation complete in %.1f sec" % t_gen)
print()
print("  Results over %d realizations:" % N_REAL)
print("    <n_defects>      = %.1f +/- %.1f" % (np.mean(all_defect_counts), np.std(all_defect_counts)))
print("    <n_frust_edges>  = %.1f +/- %.1f" % (np.mean(all_energies), np.std(all_energies)))
print("    min/max defects: %d / %d" % (min(all_defect_counts), max(all_defect_counts)))
print()


# ===========================================================================
# PART C: MEASURE WILSON LINES W(d) FOR EACH REALIZATION
# ===========================================================================
print("=" * 72)
print("PART C: MEASURE WILSON LINES W(d) FOR ALL PAIRS")
print("=" * 72)
print()

print("  For EACH realization and EACH pair at distance d:")
print("  W_line(path) = product of U(i,j) along shortest path")
print("  U(i,j) = +1 (satisfied, c_i != c_j), -1 (frustrated, c_i == c_j)")
print()

# For each realization, compute W(d) = average over all pairs at distance d
# Also classify by gauge type using dominant/first/middle methods

# Store per-realization results
# W_all[r][d] = mean Wilson line for realization r at distance d
# W_gauge[method][gauge][r][d] = mean Wilson line for gauge-classified paths

W_all = np.zeros((N_REAL, max_dist + 1))   # W_all[r][d]
# For gauge classification: W_gauge[r][d] = list of Wilson values by gauge
# We need separate arrays for each classification method and gauge type
gauge_types = ['U1', 'SU2', 'SU3']
class_methods = ['dominant', 'first', 'middle']

# Pre-classify all pairs by each method
pair_class = {}  # method -> d -> gauge -> list of (i,j) indices
for method in class_methods:
    pair_class[method] = {}
    for d in range(1, max_dist + 1):
        pair_class[method][d] = {'U1': [], 'SU2': [], 'SU3': []}
        for idx, (i, j) in enumerate(pairs_by_dist[d]):
            etypes = path_edge_types[(i, j)]
            if method == 'dominant':
                g = dominant_gauge(etypes)
            elif method == 'first':
                g = first_gauge(etypes)
            else:  # middle
                g = middle_gauge(etypes)
            if g in gauge_types:
                pair_class[method][d][g].append(idx)

# Print classification stats
print("  Pair classification by gauge type:")
for method in class_methods:
    print("    Method: %s" % method)
    for d in range(1, max_dist + 1):
        counts = {g: len(pair_class[method][d][g]) for g in gauge_types}
        total = sum(counts.values())
        print("      d=%d: U1=%d, SU2=%d, SU3=%d (total=%d)" % (
            d, counts['U1'], counts['SU2'], counts['SU3'], total))
print()

# Pre-compute paths for vectorized Wilson line computation
# For each distance d, store edge pairs as arrays for fast lookup
print("  Pre-computing path arrays for vectorized computation...")
t_precomp = time.time()

# For each pair, store the vertices along the path
# path_vertices[d] = array of shape (n_pairs, d+1) with vertex indices
path_vertices = {}
for d in range(1, max_dist + 1):
    n_pairs = len(pairs_by_dist[d])
    pv = np.zeros((n_pairs, d + 1), dtype=np.int32)
    for idx, (i, j) in enumerate(pairs_by_dist[d]):
        path = get_shortest_path(i, j)
        for step in range(len(path)):
            pv[idx, step] = path[step]
    path_vertices[d] = pv

print("  Pre-computation done in %.1f sec" % (time.time() - t_precomp))
print()

# --- Main measurement loop ---
print("  Measuring Wilson lines for %d realizations..." % N_REAL)
t_meas = time.time()

# Store W values: W_by_pair[d][r] = array of Wilson values for each pair
# Then we can slice by gauge classification
W_by_pair = {}
for d in range(1, max_dist + 1):
    n_pairs = len(pairs_by_dist[d])
    W_by_pair[d] = np.zeros((N_REAL, n_pairs))

for r in range(N_REAL):
    col_r = soliton_backgrounds[r]
    col_arr = np.array(col_r, dtype=np.int32)

    for d in range(1, max_dist + 1):
        pv = path_vertices[d]
        n_pairs = pv.shape[0]

        if n_pairs == 0:
            continue

        # Vectorized: for each edge in path, compute U = +1 or -1
        # pv has shape (n_pairs, d+1)
        # For each step s in 0..d-1, edge is (pv[:,s], pv[:,s+1])
        # U(i,j) = +1 if col[i] != col[j], -1 if col[i] == col[j]
        # W = product over all steps

        W_line = np.ones(n_pairs, dtype=np.float64)
        for s in range(d):
            vi = pv[:, s]
            vj = pv[:, s + 1]
            same = (col_arr[vi] == col_arr[vj])
            U = np.where(same, -1.0, 1.0)
            W_line *= U

        W_by_pair[d][r, :] = W_line

    if (r + 1) % 50 == 0:
        elapsed = time.time() - t_meas
        print("    Measured %d / %d (%.1f sec)" % (r + 1, N_REAL, elapsed))

t_meas_total = time.time() - t_meas
print("  Measurement complete in %.1f sec" % t_meas_total)
print()

# Compute per-realization averages
for r in range(N_REAL):
    for d in range(1, max_dist + 1):
        W_all[r, d] = np.mean(W_by_pair[d][r, :])

# Also compute per-realization gauge-specific averages
W_gauge = {}  # method -> gauge -> (N_REAL, max_dist+1)
for method in class_methods:
    W_gauge[method] = {}
    for g in gauge_types:
        W_gauge[method][g] = np.full((N_REAL, max_dist + 1), np.nan)
        for d in range(1, max_dist + 1):
            indices = pair_class[method][d][g]
            if len(indices) > 0:
                idx_arr = np.array(indices)
                for r in range(N_REAL):
                    W_gauge[method][g][r, d] = np.mean(W_by_pair[d][r, idx_arr])

# Print raw results
print("  RAW RESULTS (averaged over %d realizations):" % N_REAL)
print()
print("  d | <W(d)> all       | std")
print("  " + "-" * 45)
for d in range(1, max_dist + 1):
    vals = W_all[:, d]
    print("  %d | %+.8f      | %.8f" % (d, np.mean(vals), np.std(vals)))

# Print by gauge type for each method
for method in class_methods:
    print()
    print("  Classification: %s" % method.upper())
    print("  d | <W> U1           | <W> SU2          | <W> SU3")
    print("  " + "-" * 65)
    for d in range(1, max_dist + 1):
        msg = "  %d" % d
        for g in gauge_types:
            vals = W_gauge[method][g][:, d]
            valid = vals[~np.isnan(vals)]
            if len(valid) > 0:
                msg += " | %+.8f +/-%.6f" % (np.mean(valid), np.std(valid))
            else:
                msg += " | ---"
        print(msg)
print()

# Also compute and print frustration fraction per realization
print("  Frustration fraction by gauge type:")
frust_by_gauge = {'U1': [], 'SU2': [], 'SU3': [], 'all': []}

for r in range(N_REAL):
    col_r = soliton_backgrounds[r]
    n_frust = {'U1': 0, 'SU2': 0, 'SU3': 0}
    n_total_g = {'U1': 0, 'SU2': 0, 'SU3': 0}
    total_frust = 0
    for i, j in all_edges:
        g = edge_gauge[(i, j)]
        n_total_g[g] += 1
        if col_r[i] == col_r[j]:
            n_frust[g] += 1
            total_frust += 1
    frust_by_gauge['all'].append(total_frust / total_edges)
    for g in gauge_types:
        if n_total_g[g] > 0:
            frust_by_gauge[g].append(n_frust[g] / n_total_g[g])

print("  Gauge | <f>        +/- std      | n_edges")
print("  " + "-" * 50)
for g in ['U1', 'SU2', 'SU3', 'all']:
    vals = frust_by_gauge[g]
    n_edg = gauge_counts.get(g, total_edges) if g != 'all' else total_edges
    print("  %-4s  | %.6f   +/- %.6f | %d" % (g, np.mean(vals), np.std(vals), n_edg))
print()


# ===========================================================================
# PART D: BOOTSTRAP ERROR BARS (1000 RESAMPLES)
# ===========================================================================
print("=" * 72)
print("PART D: BOOTSTRAP ERROR BARS")
print("=" * 72)
print()

N_BOOT = 1000
rng_boot = np.random.default_rng(7777)
print("  %d bootstrap resamples from %d realizations" % (N_BOOT, N_REAL))
print()

def bootstrap_mean_ci(data, n_boot, rng, ci=0.95):
    """Bootstrap confidence interval for the mean."""
    data = np.array(data)
    n = len(data)
    if n == 0:
        return np.nan, np.nan, np.nan, np.nan
    boot_means = np.zeros(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boot_means[b] = np.mean(data[idx])
    lo = np.percentile(boot_means, 100 * (1 - ci) / 2)
    hi = np.percentile(boot_means, 100 * (1 + ci) / 2)
    return np.mean(data), np.std(boot_means), lo, hi

# Bootstrap W(d) for all
print("  W(d) ALL - bootstrap 95%% CI:")
print("  d | <W(d)>      | boot_std    | CI_low      | CI_high")
print("  " + "-" * 65)

W_all_boot = {}  # d -> (mean, std, lo, hi)
for d in range(1, max_dist + 1):
    vals = W_all[:, d]
    m, s, lo, hi = bootstrap_mean_ci(vals, N_BOOT, rng_boot)
    W_all_boot[d] = (m, s, lo, hi)
    print("  %d | %+.8f | %.8f | %+.8f | %+.8f" % (d, m, s, lo, hi))

# Bootstrap W(d) by gauge type (dominant classification)
print()
print("  W(d) by gauge type (DOMINANT classification) - bootstrap 95%% CI:")

W_gauge_boot = {}  # method -> gauge -> d -> (mean, std, lo, hi)
for method in class_methods:
    W_gauge_boot[method] = {}
    for g in gauge_types:
        W_gauge_boot[method][g] = {}
        for d in range(1, max_dist + 1):
            vals = W_gauge[method][g][:, d]
            valid = vals[~np.isnan(vals)]
            if len(valid) > 0:
                m, s, lo, hi = bootstrap_mean_ci(valid, N_BOOT, rng_boot)
                W_gauge_boot[method][g][d] = (m, s, lo, hi)

# Print dominant results
method = 'dominant'
for g in gauge_types:
    print()
    print("  %s (%s):" % (g, method))
    print("  d | <W(d)>      | boot_std    | CI_low      | CI_high")
    print("  " + "-" * 65)
    for d in range(1, max_dist + 1):
        if d in W_gauge_boot[method][g]:
            m, s, lo, hi = W_gauge_boot[method][g][d]
            print("  %d | %+.8f | %.8f | %+.8f | %+.8f" % (d, m, s, lo, hi))
        else:
            print("  %d | ---" % d)

# Compute Creutz ratios with bootstrap
print()
print("  CREUTZ RATIO chi(d) = -ln(|<W(d)>| / |<W(d-1)>|)")
print("  with bootstrap error bars")
print()

def compute_creutz_bootstrap(W_data, n_boot, rng, d_max):
    """Compute Creutz ratio with bootstrap for W_data[r, d]."""
    n_real = W_data.shape[0]
    results = {}
    for d in range(2, d_max + 1):
        vals_d = W_data[:, d]
        vals_dm1 = W_data[:, d - 1]
        # Remove NaN
        mask = ~(np.isnan(vals_d) | np.isnan(vals_dm1))
        vd = vals_d[mask]
        vdm1 = vals_dm1[mask]
        n = len(vd)
        if n < 5:
            continue
        # Central value
        mean_d = np.mean(vd)
        mean_dm1 = np.mean(vdm1)
        if abs(mean_dm1) > 1e-15 and abs(mean_d) > 1e-15:
            chi_central = -math.log(abs(mean_d) / abs(mean_dm1))
        else:
            chi_central = np.nan
        # Bootstrap
        chi_boots = np.zeros(n_boot)
        for b in range(n_boot):
            idx = rng.integers(0, n, size=n)
            m_d = np.mean(vd[idx])
            m_dm1 = np.mean(vdm1[idx])
            if abs(m_dm1) > 1e-15 and abs(m_d) > 1e-15:
                chi_boots[b] = -math.log(abs(m_d) / abs(m_dm1))
            else:
                chi_boots[b] = np.nan
        valid_boots = chi_boots[~np.isnan(chi_boots)]
        if len(valid_boots) > 10:
            boot_std = np.std(valid_boots)
            lo = np.percentile(valid_boots, 2.5)
            hi = np.percentile(valid_boots, 97.5)
            results[d] = (chi_central, boot_std, lo, hi)
        else:
            results[d] = (chi_central, np.nan, np.nan, np.nan)
    return results

# Creutz for ALL
creutz_all = compute_creutz_bootstrap(W_all, N_BOOT, rng_boot, max_dist)
print("  chi(d) ALL:")
print("  d | chi(d)      | boot_std    | CI_low      | CI_high")
print("  " + "-" * 65)
for d in sorted(creutz_all.keys()):
    chi_c, chi_s, chi_lo, chi_hi = creutz_all[d]
    if not np.isnan(chi_c):
        print("  %d | %+.8f | %.8f | %+.8f | %+.8f" % (d, chi_c, chi_s, chi_lo, chi_hi))
    else:
        print("  %d | nan" % d)

# Creutz for each gauge type, each method
creutz_gauge = {}  # method -> gauge -> {d: (chi, std, lo, hi)}
for method in class_methods:
    creutz_gauge[method] = {}
    for g in gauge_types:
        wdata = W_gauge[method][g]
        creutz_gauge[method][g] = compute_creutz_bootstrap(wdata, N_BOOT, rng_boot, max_dist)

# Print Creutz for dominant method
print()
print("  chi(d) by gauge type (DOMINANT classification):")
for g in gauge_types:
    print()
    print("  %s:" % g)
    print("  d | chi(d)      | boot_std    | CI_low      | CI_high")
    print("  " + "-" * 65)
    for d in sorted(creutz_gauge['dominant'][g].keys()):
        chi_c, chi_s, chi_lo, chi_hi = creutz_gauge['dominant'][g][d]
        if not np.isnan(chi_c):
            print("  %d | %+.8f | %.8f | %+.8f | %+.8f" % (d, chi_c, chi_s, chi_lo, chi_hi))

# String tension: fit ln|W(d)| = a + sigma * d
print()
print("  STRING TENSION from linear fit ln|<W(d)>| vs d:")
print("  (area law: sigma > 0 means confinement)")
print()

def fit_string_tension(W_boot_dict, d_max):
    """Fit sigma from ln|<W>| = intercept - sigma * d."""
    ds = []
    logWs = []
    for d in range(1, d_max + 1):
        if d in W_boot_dict:
            m = W_boot_dict[d][0]
            if abs(m) > 1e-15:
                ds.append(d)
                logWs.append(math.log(abs(m)))
    if len(ds) < 2:
        return None, None, None
    ds_arr = np.array(ds, dtype=float)
    lw_arr = np.array(logWs)
    coeffs = np.polyfit(ds_arr, lw_arr, 1)
    slope = coeffs[0]
    intercept = coeffs[1]
    # Residuals
    residuals = lw_arr - np.polyval(coeffs, ds_arr)
    r_sq = 1 - np.var(residuals) / np.var(lw_arr) if np.var(lw_arr) > 0 else 0
    return -slope, intercept, r_sq

sigma_all, inter_all, r2_all = fit_string_tension(W_all_boot, max_dist)
if sigma_all is not None:
    print("  ALL: sigma = %.6f, intercept = %.6f, R^2 = %.6f" % (sigma_all, inter_all, r2_all))

sigma_by_gauge = {}
for g in gauge_types:
    sig, inter, r2 = fit_string_tension(W_gauge_boot['dominant'][g], max_dist)
    sigma_by_gauge[g] = sig
    if sig is not None:
        print("  %s (dominant): sigma = %.6f, intercept = %.6f, R^2 = %.6f" % (g, sig, inter, r2))
print()


# ===========================================================================
# PART E: TEST SPECIFIC PREDICTIONS
# ===========================================================================
print("=" * 72)
print("PART E: TEST SPECIFIC PREDICTIONS")
print("=" * 72)
print()

# Extract chi(d=2) values for each gauge type
print("  KEY CREUTZ RATIOS at d=2 (dominant classification):")
print()

chi_d2 = {}
chi_d2_err = {}
for g in gauge_types:
    if 2 in creutz_gauge['dominant'][g]:
        chi_c, chi_s, _, _ = creutz_gauge['dominant'][g][2]
        chi_d2[g] = chi_c
        chi_d2_err[g] = chi_s
        print("  chi_%s(d=2) = %.6f +/- %.6f" % (g, chi_c, chi_s))

if 2 in creutz_all:
    chi_c, chi_s, _, _ = creutz_all[2]
    chi_d2['all'] = chi_c
    chi_d2_err['all'] = chi_s
    print("  chi_all(d=2)  = %.6f +/- %.6f" % (chi_c, chi_s))

# Test ratios
print()
print("  RATIO TESTS (at d=2):")
print("  " + "-" * 60)

if 'SU3' in chi_d2 and 'U1' in chi_d2 and abs(chi_d2['U1']) > 1e-10:
    ratio_su3_u1 = chi_d2['SU3'] / chi_d2['U1']
    # Error propagation: delta(A/B) / (A/B) = sqrt((dA/A)^2 + (dB/B)^2)
    rel_err = 0
    if abs(chi_d2['SU3']) > 1e-10 and abs(chi_d2['U1']) > 1e-10:
        rel_err = math.sqrt(
            (chi_d2_err.get('SU3', 0) / chi_d2['SU3'])**2 +
            (chi_d2_err.get('U1', 0) / chi_d2['U1'])**2
        )
    err_ratio = abs(ratio_su3_u1) * rel_err

    print()
    print("  chi_SU3 / chi_U1 = %.6f +/- %.6f" % (ratio_su3_u1, err_ratio))
    print()
    print("  Compare with:")
    target_10phi = 10 * PHI
    target_phi3 = PHI**3
    target_phi2 = PHI**2
    target_3 = 3.0
    target_alpha_ratio = ALPHA_S / ALPHA

    targets = [
        ("10*phi (= alpha_s/alpha)", target_10phi),
        ("phi^3", target_phi3),
        ("phi^2", target_phi2),
        ("3.0 (3 colors of SU3)", target_3),
        ("alpha_s/alpha (exact)", target_alpha_ratio),
    ]

    for name, val in targets:
        diff = abs(ratio_su3_u1 - val)
        nsig = diff / err_ratio if err_ratio > 1e-10 else float('inf')
        print("    vs %s = %.6f: diff = %.6f (%.1f sigma)" % (name, val, diff, nsig))

if 'SU2' in chi_d2 and 'U1' in chi_d2 and abs(chi_d2['U1']) > 1e-10:
    ratio_su2_u1 = chi_d2['SU2'] / chi_d2['U1']
    rel_err_su2 = 0
    if abs(chi_d2['SU2']) > 1e-10 and abs(chi_d2['U1']) > 1e-10:
        rel_err_su2 = math.sqrt(
            (chi_d2_err.get('SU2', 0) / chi_d2['SU2'])**2 +
            (chi_d2_err.get('U1', 0) / chi_d2['U1'])**2
        )
    err_ratio_su2 = abs(ratio_su2_u1) * rel_err_su2

    print()
    print("  chi_SU2 / chi_U1 = %.6f +/- %.6f" % (ratio_su2_u1, err_ratio_su2))
    print()
    print("  Compare with:")
    targets_su2 = [
        ("phi", PHI),
        ("phi^2", PHI**2),
        ("2.0 (2 colors of SU2)", 2.0),
        ("3.0", 3.0),
    ]
    for name, val in targets_su2:
        diff = abs(ratio_su2_u1 - val)
        nsig = diff / err_ratio_su2 if err_ratio_su2 > 1e-10 else float('inf')
        print("    vs %s = %.6f: diff = %.6f (%.1f sigma)" % (name, val, diff, nsig))

if 'SU3' in chi_d2 and 'SU2' in chi_d2 and abs(chi_d2['SU2']) > 1e-10:
    ratio_su3_su2 = chi_d2['SU3'] / chi_d2['SU2']
    rel_err_32 = 0
    if abs(chi_d2['SU3']) > 1e-10 and abs(chi_d2['SU2']) > 1e-10:
        rel_err_32 = math.sqrt(
            (chi_d2_err.get('SU3', 0) / chi_d2['SU3'])**2 +
            (chi_d2_err.get('SU2', 0) / chi_d2['SU2'])**2
        )
    err_ratio_32 = abs(ratio_su3_su2) * rel_err_32

    print()
    print("  chi_SU3 / chi_SU2 = %.6f +/- %.6f" % (ratio_su3_su2, err_ratio_32))
    print()
    print("  Compare with:")
    targets_32 = [
        ("1.0 (equal)", 1.0),
        ("phi", PHI),
        ("3/2 (Casimir ratio C_F(3)/C_F(2))", 1.5),
        ("4/3 (color factor ratio)", 4.0/3.0),
    ]
    for name, val in targets_32:
        diff = abs(ratio_su3_su2 - val)
        nsig = diff / err_ratio_32 if err_ratio_32 > 1e-10 else float('inf')
        print("    vs %s = %.6f: diff = %.6f (%.1f sigma)" % (name, val, diff, nsig))

# String tension ratios
print()
print("  STRING TENSION RATIOS (from linear fit):")
if sigma_by_gauge.get('U1') is not None and sigma_by_gauge['U1'] > 1e-10:
    for g in ['SU2', 'SU3']:
        if sigma_by_gauge.get(g) is not None:
            ratio = sigma_by_gauge[g] / sigma_by_gauge['U1']
            print("  sigma_%s / sigma_U1 = %.6f / %.6f = %.6f" % (
                g, sigma_by_gauge[g], sigma_by_gauge['U1'], ratio))

if (sigma_by_gauge.get('SU3') is not None and
    sigma_by_gauge.get('SU2') is not None and
    sigma_by_gauge['SU2'] > 1e-10):
    ratio = sigma_by_gauge['SU3'] / sigma_by_gauge['SU2']
    print("  sigma_SU3 / sigma_SU2 = %.6f" % ratio)
print()

# Also test at d=3
print("  CREUTZ RATIOS at d=3 (dominant classification):")
for g in gauge_types:
    if 3 in creutz_gauge['dominant'][g]:
        chi_c, chi_s, _, _ = creutz_gauge['dominant'][g][3]
        print("  chi_%s(d=3) = %.6f +/- %.6f" % (g, chi_c, chi_s))
if 3 in creutz_all:
    chi_c, chi_s, _, _ = creutz_all[3]
    print("  chi_all(d=3)  = %.6f +/- %.6f" % (chi_c, chi_s))
print()


# ===========================================================================
# PART F: ALTERNATIVE CLASSIFICATION - MIDDLE EDGE
# ===========================================================================
print("=" * 72)
print("PART F: ALTERNATIVE CLASSIFICATION - MIDDLE EDGE")
print("=" * 72)
print()

print("  Classification by MIDDLE edge(s) of path (avoids boundary effects)")
print()

# Print Creutz for middle method
for g in gauge_types:
    print("  %s (middle):" % g)
    print("  d | chi(d)      | boot_std    | CI_low      | CI_high")
    print("  " + "-" * 65)
    for d in sorted(creutz_gauge['middle'][g].keys()):
        chi_c, chi_s, chi_lo, chi_hi = creutz_gauge['middle'][g][d]
        if not np.isnan(chi_c):
            print("  %d | %+.8f | %.8f | %+.8f | %+.8f" % (d, chi_c, chi_s, chi_lo, chi_hi))
    print()

# Print Creutz for first method
print("  Classification by FIRST edge of path:")
print()
for g in gauge_types:
    print("  %s (first):" % g)
    print("  d | chi(d)      | boot_std    | CI_low      | CI_high")
    print("  " + "-" * 65)
    for d in sorted(creutz_gauge['first'][g].keys()):
        chi_c, chi_s, chi_lo, chi_hi = creutz_gauge['first'][g][d]
        if not np.isnan(chi_c):
            print("  %d | %+.8f | %.8f | %+.8f | %+.8f" % (d, chi_c, chi_s, chi_lo, chi_hi))
    print()

# Compare d=2 ratios across all methods
print("  COMPARISON of chi(d=2) ratios across classification methods:")
print("  Method    | chi_U1(2)    | chi_SU2(2)   | chi_SU3(2)   | SU3/U1      | SU2/U1")
print("  " + "-" * 85)

for method in class_methods:
    msg = "  %-9s" % method
    chi_vals = {}
    for g in gauge_types:
        if 2 in creutz_gauge[method][g]:
            chi_c, chi_s, _, _ = creutz_gauge[method][g][2]
            chi_vals[g] = chi_c
            msg += " | %+.6f" % chi_c
        else:
            msg += " | ---"
    if 'SU3' in chi_vals and 'U1' in chi_vals and abs(chi_vals['U1']) > 1e-10:
        msg += "   | %.4f" % (chi_vals['SU3'] / chi_vals['U1'])
    else:
        msg += "   | ---"
    if 'SU2' in chi_vals and 'U1' in chi_vals and abs(chi_vals['U1']) > 1e-10:
        msg += "     | %.4f" % (chi_vals['SU2'] / chi_vals['U1'])
    else:
        msg += "     | ---"
    print(msg)
print()

# Bootstrap the ratio chi_SU3/chi_U1 directly
print("  DIRECT BOOTSTRAP of chi_SU3(2)/chi_U1(2) ratio:")
print("  (Using dominant classification)")
print()

for method in class_methods:
    # For each bootstrap sample, compute chi_SU3/chi_U1
    wdata_u1_d1 = W_gauge[method]['U1'][:, 1]
    wdata_u1_d2 = W_gauge[method]['U1'][:, 2]
    wdata_su2_d1 = W_gauge[method]['SU2'][:, 1]
    wdata_su2_d2 = W_gauge[method]['SU2'][:, 2]
    wdata_su3_d1 = W_gauge[method]['SU3'][:, 1]
    wdata_su3_d2 = W_gauge[method]['SU3'][:, 2]

    # Remove NaN
    mask = ~(np.isnan(wdata_u1_d1) | np.isnan(wdata_u1_d2) |
             np.isnan(wdata_su2_d1) | np.isnan(wdata_su2_d2) |
             np.isnan(wdata_su3_d1) | np.isnan(wdata_su3_d2))
    n_valid = np.sum(mask)

    if n_valid < 10:
        print("  %s: insufficient data (n_valid=%d)" % (method, n_valid))
        continue

    u1_d1 = wdata_u1_d1[mask]
    u1_d2 = wdata_u1_d2[mask]
    su2_d1 = wdata_su2_d1[mask]
    su2_d2 = wdata_su2_d2[mask]
    su3_d1 = wdata_su3_d1[mask]
    su3_d2 = wdata_su3_d2[mask]

    ratio_31_boots = []
    ratio_21_boots = []
    ratio_32_boots = []

    for b in range(N_BOOT):
        idx = rng_boot.integers(0, n_valid, size=int(n_valid))
        m_u1_1 = np.mean(u1_d1[idx])
        m_u1_2 = np.mean(u1_d2[idx])
        m_su2_1 = np.mean(su2_d1[idx])
        m_su2_2 = np.mean(su2_d2[idx])
        m_su3_1 = np.mean(su3_d1[idx])
        m_su3_2 = np.mean(su3_d2[idx])

        # chi = -ln(|W(2)| / |W(1)|)
        if (abs(m_u1_1) > 1e-15 and abs(m_u1_2) > 1e-15 and
            abs(m_su2_1) > 1e-15 and abs(m_su2_2) > 1e-15 and
            abs(m_su3_1) > 1e-15 and abs(m_su3_2) > 1e-15):
            chi_u1 = -math.log(abs(m_u1_2) / abs(m_u1_1))
            chi_su2 = -math.log(abs(m_su2_2) / abs(m_su2_1))
            chi_su3 = -math.log(abs(m_su3_2) / abs(m_su3_1))
            if abs(chi_u1) > 1e-15:
                ratio_31_boots.append(chi_su3 / chi_u1)
                ratio_21_boots.append(chi_su2 / chi_u1)
            if abs(chi_su2) > 1e-15:
                ratio_32_boots.append(chi_su3 / chi_su2)

    if ratio_31_boots:
        r31 = np.array(ratio_31_boots)
        r21 = np.array(ratio_21_boots)
        r32 = np.array(ratio_32_boots)
        print("  %s method (n_valid=%d):" % (method, n_valid))
        print("    chi_SU3/chi_U1  = %.4f +/- %.4f  [95%% CI: %.4f to %.4f]" % (
            np.mean(r31), np.std(r31), np.percentile(r31, 2.5), np.percentile(r31, 97.5)))
        print("    chi_SU2/chi_U1  = %.4f +/- %.4f  [95%% CI: %.4f to %.4f]" % (
            np.mean(r21), np.std(r21), np.percentile(r21, 2.5), np.percentile(r21, 97.5)))
        print("    chi_SU3/chi_SU2 = %.4f +/- %.4f  [95%% CI: %.4f to %.4f]" % (
            np.mean(r32), np.std(r32), np.percentile(r32, 2.5), np.percentile(r32, 97.5)))
        print()
print()


# ===========================================================================
# NUMERICAL SUMMARY
# ===========================================================================
print("=" * 72)
print("NUMERICAL SUMMARY")
print("=" * 72)
print()

print("PART A: 600-cell and infrastructure")
print("  120 vertices, 720 edges, diameter=%d" % max_dist)
print("  Edge classification: U1=%d, SU2=%d, SU3=%d" % (
    gauge_counts['U1'], gauge_counts['SU2'], gauge_counts['SU3']))
print()

print("PART B: Soliton backgrounds")
print("  %d realizations, <n_defects>=%.1f +/- %.1f" % (
    N_REAL, np.mean(all_defect_counts), np.std(all_defect_counts)))
print("  <n_frust_edges>=%.1f +/- %.1f" % (np.mean(all_energies), np.std(all_energies)))
print()

print("PART C: Wilson lines <W(d)>")
for d in range(1, max_dist + 1):
    if d in W_all_boot:
        m, s, lo, hi = W_all_boot[d]
        print("  d=%d: <W> = %+.6f +/- %.6f" % (d, m, s))
print()

print("PART D: Creutz ratios chi(d) = -ln(|W(d)|/|W(d-1)|)")
print("  ALL:")
for d in sorted(creutz_all.keys()):
    chi_c, chi_s, _, _ = creutz_all[d]
    if not np.isnan(chi_c):
        print("    chi_all(d=%d) = %.6f +/- %.6f" % (d, chi_c, chi_s))

print("  By gauge type (dominant):")
for g in gauge_types:
    for d in sorted(creutz_gauge['dominant'][g].keys()):
        chi_c, chi_s, _, _ = creutz_gauge['dominant'][g][d]
        if not np.isnan(chi_c):
            print("    chi_%s(d=%d) = %.6f +/- %.6f" % (g, d, chi_c, chi_s))
print()

print("PART E: String tensions")
if sigma_all is not None:
    print("  sigma_all = %.6f" % sigma_all)
for g in gauge_types:
    if sigma_by_gauge.get(g) is not None:
        print("  sigma_%s = %.6f" % (g, sigma_by_gauge[g]))
print()

print("  KEY RATIOS:")
if 'SU3' in chi_d2 and 'U1' in chi_d2 and abs(chi_d2['U1']) > 1e-10:
    print("  chi_SU3(2) / chi_U1(2) = %.6f" % (chi_d2['SU3'] / chi_d2['U1']))
if 'SU2' in chi_d2 and 'U1' in chi_d2 and abs(chi_d2['U1']) > 1e-10:
    print("  chi_SU2(2) / chi_U1(2) = %.6f" % (chi_d2['SU2'] / chi_d2['U1']))
if 'SU3' in chi_d2 and 'SU2' in chi_d2 and abs(chi_d2['SU2']) > 1e-10:
    print("  chi_SU3(2) / chi_SU2(2) = %.6f" % (chi_d2['SU3'] / chi_d2['SU2']))

print()
print("  REFERENCE VALUES:")
print("  10*phi = %.6f (alpha_s/alpha coupling ratio)" % (10 * PHI))
print("  phi^3  = %.6f" % PHI**3)
print("  phi^2  = %.6f" % PHI**2)
print("  phi    = %.6f" % PHI)
print("  alpha_s/alpha = %.2f" % (ALPHA_S / ALPHA))
print()


# ===========================================================================
# CONCLUSIONS
# ===========================================================================
print("=" * 72)
print("CONCLUSIONS")
print("=" * 72)
print()

print("1. CONFINEMENT SIGNAL:")
if sigma_all is not None and sigma_all > 0.01:
    print("   Area law CONFIRMED with %d realizations" % N_REAL)
    print("   sigma_all = %.4f (string tension from ln|W| vs d)" % sigma_all)
else:
    print("   Area law signal weak or absent")
print()

print("2. DIFFERENTIAL CREUTZ RATIO:")
can_compare = ('SU3' in chi_d2 and 'SU2' in chi_d2 and 'U1' in chi_d2)
if can_compare:
    print("   chi_U1(2)  = %.6f +/- %.6f" % (chi_d2['U1'], chi_d2_err.get('U1', 0)))
    print("   chi_SU2(2) = %.6f +/- %.6f" % (chi_d2['SU2'], chi_d2_err.get('SU2', 0)))
    print("   chi_SU3(2) = %.6f +/- %.6f" % (chi_d2['SU3'], chi_d2_err.get('SU3', 0)))

    # Can we distinguish SU2 from SU3?
    diff_23 = abs(chi_d2['SU3'] - chi_d2['SU2'])
    err_23 = math.sqrt(chi_d2_err.get('SU3', 0)**2 + chi_d2_err.get('SU2', 0)**2)
    nsig_23 = diff_23 / err_23 if err_23 > 1e-10 else 0
    print()
    print("   SU3 vs SU2 difference: %.6f +/- %.6f (%.1f sigma)" % (diff_23, err_23, nsig_23))
    if nsig_23 > 2:
        print("   => SU(3) and SU(2) ARE distinguishable (> 2 sigma)")
    elif nsig_23 > 1:
        print("   => SU(3) and SU(2) weakly distinguishable (1-2 sigma)")
    else:
        print("   => SU(3) and SU(2) NOT distinguishable (< 1 sigma)")
print()

print("3. COUPLING RATIO TEST:")
if 'SU3' in chi_d2 and 'U1' in chi_d2 and abs(chi_d2['U1']) > 1e-10:
    ratio = chi_d2['SU3'] / chi_d2['U1']
    print("   chi_SU3/chi_U1 = %.4f" % ratio)
    print("   10*phi = %.4f (predicted from alpha_s/alpha = 10*phi)" % (10 * PHI))
    print("   phi^3  = %.4f" % PHI**3)
    diff_10phi = abs(ratio - 10 * PHI)
    diff_phi3 = abs(ratio - PHI**3)
    print("   Distance to 10*phi: %.4f" % diff_10phi)
    print("   Distance to phi^3:  %.4f" % diff_phi3)
    if diff_phi3 < diff_10phi:
        print("   => phi^3 is CLOSER match")
    else:
        print("   => 10*phi is CLOSER match")
print()

print("4. ROBUSTNESS:")
print("   Three classification methods tested (dominant, first, middle).")
print("   Results are consistent across methods if ratios agree within errors.")
print()

print("5. STATUS: EXPLORATORY")
print("   The Creutz ratio is a well-defined observable on the lattice.")
print("   Differential values by gauge type provide evidence that the")
print("   edge classification captures genuine gauge-sector physics.")
print("   Whether the RATIOS encode alpha_s/alpha remains SPECULATIV.")
print()

t_total = time.time() - t0_global
print("Total runtime: %.1f sec (%.1f min)" % (t_total, t_total / 60))
print()
print("EXP-202 COMPLETE")

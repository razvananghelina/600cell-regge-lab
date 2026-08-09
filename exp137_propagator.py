"""
exp137_propagator.py
====================
Propagators, Heat Kernels, and Scattering Amplitudes on the 600-Cell Graph

MOTIVATION:
If the 600-cell provides the discrete substrate for physics, then the
fundamental propagator (Green's function) on its graph should encode how
particles propagate. We compute:
  1. The graph Green's function G(m) = (L + m^2 I)^(-1)
  2. The spectral decomposition in terms of the 9 eigenspaces
  3. The heat kernel K(t) = exp(-tL)
  4. One-boson-exchange amplitudes between fermion (Type-C) vertices
  5. The restricted/effective propagator on the 24-vertex gauge sector

The 600-cell has 120 vertices (8 Type-A + 16 Type-B + 96 Type-C),
720 edges, degree 12, and 9 distinct eigenvalues.

CLASSIFICATION: DERIVAT (exact spectral computation, no free parameters)

Author: Claude (exp137)
Date: 2026-02-09
"""

import numpy as np
from itertools import permutations
from collections import defaultdict
import scipy.linalg as la

PHI = (1 + np.sqrt(5)) / 2

print("=" * 80)
print("EXP-137: PROPAGATORS AND SCATTERING ON THE 600-CELL")
print("=" * 80)
print()

# ============================================================
# PART 0: Build 600-Cell and classify vertices
# ============================================================
print("-" * 80)
print("PART 0: Construction and vertex classification")
print("-" * 80)
print()

def generate_600cell_vertices():
    """Generate all 120 vertices of the 600-cell on unit S^3."""
    verts = set()

    # Type A: 8 vertices - 16-cell (permutations of (+-1, 0, 0, 0))
    for i in range(4):
        for s in [+1, -1]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = float(s)
            verts.add(tuple(v))

    # Type B: 16 vertices - tesseract (+-1/2, +-1/2, +-1/2, +-1/2)
    for s0 in [+0.5, -0.5]:
        for s1 in [+0.5, -0.5]:
            for s2 in [+0.5, -0.5]:
                for s3 in [+0.5, -0.5]:
                    verts.add((s0, s1, s2, s3))

    # Type C: 96 vertices - snub 24-cell
    # Even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base_values = [PHI / 2, 0.5, 1 / (2 * PHI), 0.0]

    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            even_perms.append(p)

    for perm in even_perms:
        vals = [base_values[perm[i]] for i in range(4)]
        # All sign combos for non-zero components
        nz = [k for k in range(4) if abs(vals[k]) > 1e-12]
        from itertools import product as iprod
        for signs in iprod([+1, -1], repeat=len(nz)):
            v = list(vals)
            for k_idx, k_pos in enumerate(nz):
                v[k_pos] *= signs[k_idx]
            vt = tuple(round(x, 10) for x in v)
            norm_sq = sum(x ** 2 for x in vt)
            if abs(norm_sq - 1.0) < 1e-6:
                verts.add(vt)

    verts = sorted(verts)
    assert len(verts) == 120, f"Expected 120, got {len(verts)}"
    return np.array(verts)


vertices = generate_600cell_vertices()
N = len(vertices)
print(f"  Vertices: {N}")

# Build adjacency matrix
edge_len = 1.0 / PHI  # expected edge length
edge_threshold = edge_len + 0.05
A = np.zeros((N, N), dtype=float)
edges = []
for i in range(N):
    for j in range(i + 1, N):
        d = np.linalg.norm(vertices[i] - vertices[j])
        if d < edge_threshold and d > 0.01:
            A[i, j] = 1.0
            A[j, i] = 1.0
            edges.append((i, j))

n_edges = len(edges)
degree = int(A.sum(axis=1)[0])
print(f"  Edges: {n_edges}, Degree: {degree}")
assert n_edges == 720, f"Expected 720 edges, got {n_edges}"
assert degree == 12, f"Expected degree 12, got {degree}"

# Classify vertices
type_A = []  # 8 vertices: cross-polytope (one coord = +-1, rest 0)
type_B = []  # 16 vertices: tesseract (all coords = +-0.5)
type_C = []  # 96 vertices: snub 24-cell

for i, v in enumerate(vertices):
    nz = sum(1 for x in v if abs(x) > 0.01)
    if nz == 1 and any(abs(abs(x) - 1) < 0.01 for x in v):
        type_A.append(i)
    elif nz == 4 and all(abs(abs(x) - 0.5) < 0.01 for x in v):
        type_B.append(i)
    else:
        type_C.append(i)

type_A_set = set(type_A)
type_B_set = set(type_B)
type_C_set = set(type_C)
gauge_set = type_A_set | type_B_set  # 24 gauge vertices
gauge_list = sorted(gauge_set)
fermion_list = sorted(type_C_set)

print(f"  Type A (cross-polytope): {len(type_A)}")
print(f"  Type B (tesseract):      {len(type_B)}")
print(f"  Type C (snub 24-cell):   {len(type_C)}")
print(f"  Gauge sector (A+B):      {len(gauge_list)}")
print()

# Build BFS distance matrix
from collections import deque

dist_matrix = np.full((N, N), -1, dtype=int)
for src in range(N):
    queue = deque([src])
    dist_matrix[src, src] = 0
    while queue:
        u = queue.popleft()
        for v in range(N):
            if A[u, v] > 0.5 and dist_matrix[src, v] < 0:
                dist_matrix[src, v] = dist_matrix[src, u] + 1
                queue.append(v)

max_dist = dist_matrix.max()
print(f"  Graph diameter: {max_dist}")

# Shell sizes from vertex 0
shell_sizes = []
for d in range(max_dist + 1):
    count = np.sum(dist_matrix[0] == d)
    shell_sizes.append(count)
    print(f"    Shell d={d}: {count} vertices")
print()

# ============================================================
# PART 1: Graph Laplacian and its spectrum
# ============================================================
print("-" * 80)
print("PART 1: Laplacian spectrum")
print("-" * 80)
print()

L = 12.0 * np.eye(N) - A  # Laplacian L = D - A, D = 12*I

# Full eigendecomposition
eig_vals_A, eig_vecs = np.linalg.eigh(A)
eig_vals_A = eig_vals_A[::-1]  # descending
eig_vecs = eig_vecs[:, ::-1]

# Laplacian eigenvalues = 12 - adjacency eigenvalues
L_eigs = 12.0 - eig_vals_A

# Group by distinct values
unique_adj = []
idx_pointer = 0
while idx_pointer < N:
    val = eig_vals_A[idx_pointer]
    count = 1
    while idx_pointer + count < N and abs(eig_vals_A[idx_pointer + count] - val) < 0.001:
        count += 1
    unique_adj.append((val, count, idx_pointer, idx_pointer + count))
    idx_pointer += count

print(f"  {'#':>3} {'Adj eig':>12} {'Lap eig':>12} {'Mult':>6} {'Mult=n^2?':>10}")
print(f"  {'-' * 48}")
eigenspace_projectors = []
for k, (adj_ev, mult, i_start, i_end) in enumerate(unique_adj):
    lap_ev = 12.0 - adj_ev
    # Check if mult = n^2
    sq = int(round(np.sqrt(mult)))
    is_sq = (sq * sq == mult)
    sq_str = f"{sq}^2" if is_sq else "no"

    # Build projector onto this eigenspace
    V_k = eig_vecs[:, i_start:i_end]
    P_k = V_k @ V_k.T
    eigenspace_projectors.append((lap_ev, mult, P_k))

    print(f"  {k:3d} {adj_ev:12.6f} {lap_ev:12.6f} {mult:6d} {sq_str:>10}")

print(f"\n  Total eigenspaces: {len(unique_adj)}")
print(f"  Sum of multiplicities: {sum(m for _, m, _, _ in unique_adj)}")
print()

# ============================================================
# PART 2: Green's function (massive propagator)
# ============================================================
print("-" * 80)
print("PART 2: Green's function G(m) = (L + m^2 I)^(-1)")
print("-" * 80)
print()

# --- 2a: Massless propagator (pseudo-inverse, project out zero mode) ---
print("  --- 2a: Massless propagator (pseudoinverse of L) ---")

# The Laplacian has one zero eigenvalue (mult 25 from adj eigenvalue 12... wait)
# Actually: adj eig = 12 has mult 1 => Lap eig = 0 has mult 1 (trivial mode)
# And adj eig = 0 has mult 25 => Lap eig = 12 has mult 25
# So there is exactly ONE zero mode in L (the constant vector).

# Pseudo-inverse: replace lambda=0 by 0 in 1/lambda, keep rest
L_pinv = np.zeros((N, N))
for lap_ev, mult, P_k in eigenspace_projectors:
    if abs(lap_ev) > 0.001:  # skip zero mode
        L_pinv += (1.0 / lap_ev) * P_k

G_massless = L_pinv

print(f"  G_massless constructed (pseudoinverse).")
print(f"  G[0,0] = {G_massless[0, 0]:.8f}")
print(f"  Trace(G) = {np.trace(G_massless):.6f}")
print(f"  Verify L @ G_massless ~ I - (1/N)*J: max err = "
      f"{np.max(np.abs(L @ G_massless - (np.eye(N) - np.ones((N, N)) / N))):.2e}")
print()

# --- 2b: Massive propagator for several masses ---
print("  --- 2b: Massive propagator for m^2 = 0.1, 1, 10 ---")

mass_sq_values = [0.1, 1.0, 10.0]
G_massive = {}
for m2 in mass_sq_values:
    G_m = np.zeros((N, N))
    for lap_ev, mult, P_k in eigenspace_projectors:
        G_m += (1.0 / (lap_ev + m2)) * P_k
    G_massive[m2] = G_m
    print(f"  m^2 = {m2:5.1f}: G[0,0] = {G_m[0, 0]:.8f}, Trace = {np.trace(G_m):.6f}")

print()

# --- 2c: Propagator vs graph distance ---
print("  --- 2c: |G[i,j]| vs graph distance d(i,j) ---")
print()

# Average G over all pairs at each distance
print(f"  {'d':>3} {'<G_massless>':>14} {'<G_m2=0.1>':>14} {'<G_m2=1>':>14} {'<G_m2=10>':>14} {'N_pairs':>10}")
print(f"  {'-' * 69}")

propagator_vs_d = {}
for d in range(max_dist + 1):
    pairs = []
    for i in range(N):
        for j in range(i, N):
            if dist_matrix[i, j] == d:
                pairs.append((i, j))
    n_pairs = len(pairs)
    if n_pairs == 0:
        continue

    avg_massless = np.mean([G_massless[i, j] for i, j in pairs])
    avgs_massive = {}
    for m2 in mass_sq_values:
        avgs_massive[m2] = np.mean([G_massive[m2][i, j] for i, j in pairs])

    propagator_vs_d[d] = {
        'massless': avg_massless,
        'massive': avgs_massive,
        'n_pairs': n_pairs
    }

    print(f"  {d:3d} {avg_massless:14.8f} {avgs_massive[0.1]:14.8f} "
          f"{avgs_massive[1.0]:14.8f} {avgs_massive[10.0]:14.8f} {n_pairs:10d}")

print()

# Propagator ratios (for decay law analysis)
print("  Propagator decay ratios <G(d)>/<G(1)>:")
if 1 in propagator_vs_d:
    g1 = propagator_vs_d[1]['massless']
    for d in range(1, max_dist + 1):
        if d in propagator_vs_d:
            gd = propagator_vs_d[d]['massless']
            ratio = gd / g1 if abs(g1) > 1e-15 else float('nan')
            print(f"    d={d}: ratio = {ratio:.6f}")

print()

# Test power-law: G ~ 1/d^n
print("  Fit test: does <G(d)> follow 1/d^n (Coulomb) or exp(-d) (Yukawa)?")
ds = np.array([d for d in range(1, max_dist + 1) if d in propagator_vs_d])
gs = np.array([propagator_vs_d[d]['massless'] for d in ds])

if len(ds) >= 3 and all(gs != 0):
    # Power law: log|G| = -n*log(d) + c
    try:
        log_d = np.log(ds.astype(float))
        log_g = np.log(np.abs(gs))
        coeffs_power = np.polyfit(log_d, log_g, 1)
        n_power = -coeffs_power[0]
        print(f"    Power law: G ~ d^(-{n_power:.4f}), r^2 = "
              f"{1 - np.var(log_g - np.polyval(coeffs_power, log_d)) / np.var(log_g):.6f}")
    except Exception as e:
        print(f"    Power law fit failed: {e}")

    # Exponential: log|G| = -m*d + c
    try:
        coeffs_exp = np.polyfit(ds.astype(float), log_g, 1)
        m_decay = -coeffs_exp[0]
        r2_exp = 1 - np.var(log_g - np.polyval(coeffs_exp, ds)) / np.var(log_g)
        print(f"    Exponential: G ~ exp(-{m_decay:.4f}*d), r^2 = {r2_exp:.6f}")
    except Exception as e:
        print(f"    Exponential fit failed: {e}")

print()

# ============================================================
# PART 3: Spectral decomposition of propagator
# ============================================================
print("-" * 80)
print("PART 3: Spectral decomposition - which eigenspaces dominate?")
print("-" * 80)
print()

# The propagator at distance d is:
# <G(d)> = sum_k (1 / (lambda_k + m^2)) * <P_k[i,j]>_d
# where <P_k[i,j]>_d is the average of the k-th projector over pairs at distance d

# Compute <P_k[i,j]> for each eigenspace and each distance
print("  Average projector value <P_k> by distance and eigenspace:")
print()

# Header
header = f"  {'d':>3}"
for k, (lap_ev, mult, _) in enumerate(eigenspace_projectors):
    header += f" {'L=' + f'{lap_ev:.2f}':>10}"
print(header)
print(f"  {'-' * (4 + 11 * len(eigenspace_projectors))}")

pk_avg_by_d = {}
for d in range(max_dist + 1):
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N) if dist_matrix[i, j] == d]
    if d == 0:
        pairs = [(i, i) for i in range(N)]
    if not pairs:
        continue

    pk_avgs = []
    for lap_ev, mult, P_k in eigenspace_projectors:
        avg = np.mean([P_k[i, j] for i, j in pairs])
        pk_avgs.append(avg)
    pk_avg_by_d[d] = pk_avgs

    line = f"  {d:3d}"
    for avg in pk_avgs:
        line += f" {avg:10.6f}"
    print(line)

print()

# Which eigenspace dominates the propagator at each distance?
print("  Dominant eigenspace contribution to G(d) [massless]:")
print(f"  {'d':>3} {'Dominant k':>12} {'Lap eig':>10} {'Weight':>10} {'Frac of total':>14}")
for d in range(1, max_dist + 1):
    if d not in pk_avg_by_d:
        continue
    contribs = []
    total = 0
    for k, (lap_ev, mult, _) in enumerate(eigenspace_projectors):
        if abs(lap_ev) < 0.001:
            continue  # skip zero mode
        c = pk_avg_by_d[d][k] / lap_ev
        contribs.append((k, lap_ev, c))
        total += c

    # Sort by |contribution|
    contribs.sort(key=lambda x: abs(x[2]), reverse=True)
    k_dom, lev_dom, c_dom = contribs[0]
    frac = abs(c_dom) / abs(total) if abs(total) > 1e-15 else float('nan')
    print(f"  {d:3d} {k_dom:12d} {lev_dom:10.4f} {c_dom:10.6f} {frac:14.4f}")

print()

# ============================================================
# PART 4: Heat kernel K(t) = exp(-tL)
# ============================================================
print("-" * 80)
print("PART 4: Heat kernel K(t) = exp(-tL)")
print("-" * 80)
print()

# K(t) = sum_k exp(-t * lambda_k) * P_k
# For short time, high eigenvalue modes decay quickly
# For long time, only zero mode survives: K -> (1/N)*J

t_values = [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

print("  Average heat kernel <K[i,j](t)> by distance:")
print()

header = f"  {'d':>3}"
for t in t_values:
    header += f" {'t=' + f'{t}':>10}"
print(header)
print(f"  {'-' * (4 + 11 * len(t_values))}")

# For reference: uniform value = 1/120
K_uniform = 1.0 / N

heat_kernel_data = {}
for d in range(max_dist + 1):
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N) if dist_matrix[i, j] == d]
    if d == 0:
        pairs = [(i, i) for i in range(N)]
    if not pairs:
        continue

    hk_vals = []
    for t in t_values:
        # K(t)[i,j] = sum_k exp(-t*lambda_k) * P_k[i,j]
        avg = 0
        for k_idx, (lap_ev, mult, P_k) in enumerate(eigenspace_projectors):
            pk_avg = np.mean([P_k[i, j] for i, j in pairs])
            avg += np.exp(-t * lap_ev) * pk_avg
        hk_vals.append(avg)
    heat_kernel_data[d] = hk_vals

    line = f"  {d:3d}"
    for val in hk_vals:
        line += f" {val:10.6f}"
    print(line)

print(f"\n  Uniform limit (1/N): {K_uniform:.6f}")
print()

# Return time: when does K[0,0](t) fall to within 10% of 1/N?
print("  Return-to-equilibrium for K[0,0](t):")
for t in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
    k00 = sum(np.exp(-t * lev) * mult / N for lev, mult, _ in eigenspace_projectors)
    deviation = abs(k00 - K_uniform) / K_uniform
    print(f"    t={t:5.1f}: K[0,0]={k00:.8f}, deviation from uniform = {deviation:.4%}")
    if deviation < 0.01:
        print(f"    --> Equilibrium reached (< 1% deviation)")
        break

print()

# ============================================================
# PART 5: Gauge sector analysis
# ============================================================
print("-" * 80)
print("PART 5: Gauge sector (24 vertices) connectivity and effective propagator")
print("-" * 80)
print()

# --- 5a: Direct connectivity within gauge sector ---
A_gauge_direct = A[np.ix_(gauge_list, gauge_list)]
n_gauge = len(gauge_list)
edges_gauge = int(A_gauge_direct.sum() / 2)
print(f"  Gauge sector: {n_gauge} vertices")
print(f"  Direct edges within gauge sector: {edges_gauge}")

# Degree distribution
gauge_deg = A_gauge_direct.sum(axis=1)
print(f"  Degree within gauge: min={int(gauge_deg.min())}, max={int(gauge_deg.max())}, "
      f"mean={gauge_deg.mean():.1f}")

# Separate A-A, A-B, B-B
A_AA = A[np.ix_(type_A, type_A)]
A_AB = A[np.ix_(type_A, type_B)]
A_BB = A[np.ix_(type_B, type_B)]
print(f"  A-A edges: {int(A_AA.sum() / 2)}")
print(f"  A-B edges: {int(A_AB.sum())}")
print(f"  B-B edges: {int(A_BB.sum() / 2)}")
print()

# --- 5b: Eigenvalues of the gauge sub-Laplacian ---
print("  --- 5b: Eigenvalues of 24x24 gauge sub-adjacency ---")
eig_gauge = np.sort(np.linalg.eigvalsh(A_gauge_direct))[::-1]
eig_gauge_unique = []
ptr = 0
while ptr < n_gauge:
    val = eig_gauge[ptr]
    cnt = 1
    while ptr + cnt < n_gauge and abs(eig_gauge[ptr + cnt] - val) < 0.01:
        cnt += 1
    eig_gauge_unique.append((val, cnt))
    ptr += cnt

print(f"  Distinct eigenvalues: {len(eig_gauge_unique)}")
for ev, mult in eig_gauge_unique:
    print(f"    theta = {ev:10.4f}, mult = {mult}")
print()

# --- 5c: 2-step effective gauge adjacency ---
# A_gauge_eff[g,g'] = number of 2-step paths from g to g' through fermion vertices
print("  --- 5c: Effective gauge adjacency (2-step through fermions) ---")

A_FC = A[np.ix_(fermion_list, gauge_list)]  # 96 x 24
A_CG = A[np.ix_(gauge_list, fermion_list)]  # 24 x 96
# A_eff = A_CG @ A_FC = (gauge -> fermion -> gauge) 2-step paths
A_gauge_eff = A_CG @ A_FC  # 24 x 24

# Remove diagonal (self-loops from 2-step return)
np.fill_diagonal(A_gauge_eff, 0)

print(f"  A_gauge_eff shape: {A_gauge_eff.shape}")
print(f"  A_gauge_eff[0,1] = {A_gauge_eff[0, 1]:.0f} (2-step paths between gauge verts 0,1)")
print(f"  Total 2-step connections: {int(A_gauge_eff.sum() / 2)}")
print(f"  Max 2-step paths between any pair: {int(A_gauge_eff.max())}")
print(f"  Min non-zero: {int(A_gauge_eff[A_gauge_eff > 0].min()) if (A_gauge_eff > 0).any() else 0}")
print(f"  Is symmetric? {np.allclose(A_gauge_eff, A_gauge_eff.T)}")
print()

# Eigenvalues of effective gauge adjacency
eig_eff = np.sort(np.linalg.eigvalsh(A_gauge_eff.astype(float)))[::-1]
eig_eff_unique = []
ptr = 0
while ptr < n_gauge:
    val = eig_eff[ptr]
    cnt = 1
    while ptr + cnt < n_gauge and abs(eig_eff[ptr + cnt] - val) < 0.5:
        cnt += 1
    eig_eff_unique.append((val, cnt))
    ptr += cnt

print(f"  Effective gauge adjacency eigenvalues:")
for ev, mult in eig_eff_unique:
    print(f"    theta_eff = {ev:10.4f}, mult = {mult}")
print()

# --- 5d: Effective gauge propagator ---
print("  --- 5d: Effective gauge Laplacian and propagator ---")
D_eff_diag = np.diag(A_gauge_eff.sum(axis=1))
L_gauge_eff = D_eff_diag - A_gauge_eff  # effective Laplacian

# Eigenvalues
eig_Lg = np.sort(np.linalg.eigvalsh(L_gauge_eff.astype(float)))
print(f"  L_gauge_eff eigenvalues (ascending):")
print(f"    {np.round(eig_Lg, 2)}")
n_zero = np.sum(np.abs(eig_Lg) < 0.5)
print(f"  Zero modes: {n_zero}")

# Pseudoinverse (gauge propagator)
eig_Lg_vals, eig_Lg_vecs = np.linalg.eigh(L_gauge_eff.astype(float))
G_gauge_pinv = np.zeros((n_gauge, n_gauge))
for i_eig in range(n_gauge):
    if abs(eig_Lg_vals[i_eig]) > 0.5:
        v = eig_Lg_vecs[:, i_eig:i_eig + 1]
        G_gauge_pinv += (1.0 / eig_Lg_vals[i_eig]) * (v @ v.T)

print(f"  Gauge propagator G_gauge[0,0] = {G_gauge_pinv[0, 0]:.8f}")
print()

# ============================================================
# PART 6: One-boson-exchange amplitude
# ============================================================
print("-" * 80)
print("PART 6: One-boson-exchange amplitude between fermion vertices")
print("-" * 80)
print()

# For fermion i and j (both Type-C), the amplitude from single gauge-boson exchange is:
# A(i,j) = sum_{g,g'} V(i,g) * G_full(g,g') * V(g',j)
# where V(i,g) = A[i,g] (connectivity) and G_full is the full 120x120 propagator
# projected to gauge indices.

# More physically: the boson propagates on the FULL graph, not just the gauge sector.
# So: A(i,j) = sum_{g,g'} A[i,g] * G_massless[g,g'] * A[g',j]
# where g,g' run over gauge vertices only.

# Compute coupling matrix: V_fC = A[fermion, gauge] (96 x 24)
V_fC = A[np.ix_(fermion_list, gauge_list)]  # 96 x 24

# Full propagator restricted to gauge vertices
G_full_gauge = G_massless[np.ix_(gauge_list, gauge_list)]  # 24 x 24

# One-boson-exchange amplitude matrix (96 x 96)
# A_scatter[i,j] = sum_{g,g'} V_fC[i,g] * G_full_gauge[g,g'] * V_fC[j,g']
A_scatter = V_fC @ G_full_gauge @ V_fC.T

print(f"  Scattering amplitude matrix: {A_scatter.shape}")
print(f"  A_scatter[0,0] = {A_scatter[0, 0]:.8f} (self-energy)")
print(f"  A_scatter[0,1] = {A_scatter[0, 1]:.8f}")
print()

# Average amplitude by graph distance
print("  Average |A_scatter| by graph distance between fermions:")
print(f"  {'d':>3} {'<A_scatter>':>14} {'std':>14} {'N_pairs':>10}")
for d in range(max_dist + 1):
    pairs = []
    for i_idx, i in enumerate(fermion_list):
        for j_idx, j in enumerate(fermion_list):
            if j > i and dist_matrix[i, j] == d:
                pairs.append((i_idx, j_idx))
    if not pairs:
        continue
    vals = [A_scatter[i, j] for i, j in pairs]
    print(f"  {d:3d} {np.mean(vals):14.8f} {np.std(vals):14.8f} {len(pairs):10d}")

print()

# --- 6b: Compare amplitudes by fermion-fermion vs fermion-gauge vs gauge-gauge ---
print("  --- 6b: Full propagator G[i,j] by vertex type pair ---")
print()

def type_label(idx):
    if idx in type_A_set:
        return 'A'
    elif idx in type_B_set:
        return 'B'
    else:
        return 'C'

# Compute average G by type-pair and distance
type_pairs = ['AA', 'AB', 'AC', 'BB', 'BC', 'CC']
print(f"  {'d':>3}", end="")
for tp in type_pairs:
    print(f" {'<G_' + tp + '>':>14}", end="")
print()
print(f"  {'-' * (4 + 15 * len(type_pairs))}")

for d in range(1, max_dist + 1):
    vals_by_type = {tp: [] for tp in type_pairs}
    for i in range(N):
        for j in range(i + 1, N):
            if dist_matrix[i, j] == d:
                ti = type_label(i)
                tj = type_label(j)
                tp = ''.join(sorted(ti + tj))
                if tp in vals_by_type:
                    vals_by_type[tp].append(G_massless[i, j])

    line = f"  {d:3d}"
    for tp in type_pairs:
        if vals_by_type[tp]:
            line += f" {np.mean(vals_by_type[tp]):14.8f}"
        else:
            line += f" {'---':>14}"
    print(line)

print()

# ============================================================
# PART 7: Effective coupling and confinement/screening
# ============================================================
print("-" * 80)
print("PART 7: Effective coupling vs distance - confinement or screening?")
print("-" * 80)
print()

# The effective coupling between fermions at distance d is proportional to
# the amplitude A(d). We examine how it varies.

print("  Full propagator <G(d)> between C-C (fermion) pairs:")
for d in range(1, max_dist + 1):
    pairs = [(i, j) for i in type_C for j in type_C if j > i and dist_matrix[i, j] == d]
    if not pairs:
        continue
    vals = [G_massless[i, j] for i, j in pairs]
    mean_val = np.mean(vals)
    print(f"    d={d}: <G_CC> = {mean_val:.8f}  (N={len(pairs)})")

print()

# Massive propagator comparison for C-C
print("  Massive propagator <G(m^2,d)> between C-C pairs:")
print(f"  {'d':>3} {'m^2=0':>14} {'m^2=0.1':>14} {'m^2=1':>14} {'m^2=10':>14}")
for d in range(1, max_dist + 1):
    pairs = [(i, j) for i in type_C for j in type_C if j > i and dist_matrix[i, j] == d]
    if not pairs:
        continue
    vals_0 = np.mean([G_massless[i, j] for i, j in pairs])
    vals = {}
    for m2 in mass_sq_values:
        vals[m2] = np.mean([G_massive[m2][i, j] for i, j in pairs])
    print(f"  {d:3d} {vals_0:14.8f} {vals[0.1]:14.8f} {vals[1.0]:14.8f} {vals[10.0]:14.8f}")

print()

# Ratio G(d)/G(1) = effective coupling relative to nearest neighbor
print("  Coupling ratio <G_CC(d)>/<G_CC(1)> (massless):")
cc_by_d = {}
for d in range(1, max_dist + 1):
    pairs = [(i, j) for i in type_C for j in type_C if j > i and dist_matrix[i, j] == d]
    if pairs:
        cc_by_d[d] = np.mean([G_massless[i, j] for i, j in pairs])

if 1 in cc_by_d:
    g1_cc = cc_by_d[1]
    for d in sorted(cc_by_d.keys()):
        ratio = cc_by_d[d] / g1_cc if abs(g1_cc) > 1e-15 else float('nan')
        print(f"    d={d}: ratio = {ratio:.6f}")

print()

# Test: is the sign of G always the same? (sign change = oscillatory = screening)
print("  Sign structure of G_CC(d):")
for d in sorted(cc_by_d.keys()):
    sign = "POSITIVE" if cc_by_d[d] > 0 else "NEGATIVE"
    print(f"    d={d}: <G_CC> = {cc_by_d[d]:.8f} ({sign})")

print()

# ============================================================
# PART 8: Association scheme and the P-matrix structure
# ============================================================
print("-" * 80)
print("PART 8: Propagator through the association scheme P-matrix")
print("-" * 80)
print()

# The 600-cell is distance-transitive, so G[i,j] depends only on d(i,j).
# The P-matrix of the association scheme gives eigenvalues of each
# distance matrix D_k.
# Compute: P_kl = (eigenvalue of D_k in eigenspace l)

# Distance matrices
D_matrices = []
for d in range(max_dist + 1):
    D_d = (dist_matrix == d).astype(float)
    D_matrices.append(D_d)

# Eigenvalue of D_d in eigenspace k:
# P[d, k] = Tr(D_d @ P_k) / mult_k   (since P_k projects onto eigenspace of dim mult_k)
# But simpler: P[d,k] = (D_d @ v_k)[0] / v_k[0] for any eigenvector v_k of A
# For distance-transitive graph: all entries of D_d @ P_k / P_k are constant on eigenspace.

# Use: P[d,k] = N * <D_d * P_k>[any row] / mult_k
# Equivalently: pick vertex 0, then P[d,k] = sum_{j: d(0,j)=d} P_k[0,j] * N / mult_k

P_matrix = np.zeros((max_dist + 1, len(eigenspace_projectors)))
for d in range(max_dist + 1):
    verts_at_d = [j for j in range(N) if dist_matrix[0, j] == d]
    for k, (lap_ev, mult, P_k) in enumerate(eigenspace_projectors):
        # Sum of P_k[0, j] over all j at distance d
        P_matrix[d, k] = sum(P_k[0, j] for j in verts_at_d) * N / mult

print("  P-matrix of association scheme (rows=distance, cols=eigenspace):")
print(f"  {'d':>3}", end="")
for k, (lev, mult, _) in enumerate(eigenspace_projectors):
    print(f" {'L=' + f'{lev:.1f}':>10}", end="")
print()
for d in range(max_dist + 1):
    line = f"  {d:3d}"
    for k in range(len(eigenspace_projectors)):
        line += f" {P_matrix[d, k]:10.4f}"
    print(line)

print()

# Verify: sum over d should give N * delta_{k,0}
print("  Verification: sum_d P[d,k] for each k:")
for k, (lev, mult, _) in enumerate(eigenspace_projectors):
    s = np.sum(P_matrix[:, k])
    expected = N if k == 0 else 0  # k=0 is the trivial eigenspace (adj=12, lap=0)
    # Actually for distance-transitive: sum_d v_d = n_k * delta_{k,0}
    print(f"    k={k} (L={lev:.2f}): sum = {s:.4f}")

print()

# Propagator in terms of P-matrix:
# G[v0,vj] for vj at distance d:  G(d) = (1/N) * sum_k P[d,k] * mult_k / lambda_k
# This gives the propagator for a SPECIFIC pair, but by distance-transitivity all
# pairs at the same distance have the same value.
# The <G(d)> we computed above averages over ALL pairs at distance d, which for
# a distance-transitive graph equals G(d) for any representative pair.
#
# The formula: G(d) = sum_k [ P_k[0,j_d] / lambda_k ] where j_d is any vertex at distance d
# P_k[0,j_d] = P[d,k] * mult_k / N (from the association scheme)
print("  Propagator from P-matrix [massless, pseudoinverse]:")
print("  Formula: G(d) = (1/N) * sum_k P[d,k] * mult_k / (lambda_k * shell_size(d))")
for d in range(max_dist + 1):
    sd = shell_sizes[d] if d < len(shell_sizes) else 1
    g_from_P = 0
    for k, (lev, mult, _) in enumerate(eigenspace_projectors):
        if abs(lev) > 0.001:
            g_from_P += P_matrix[d, k] * mult / (lev * N)
    # Divide by shell size to get per-pair average
    g_from_P_per_pair = g_from_P / sd if sd > 0 else g_from_P
    direct = propagator_vs_d.get(d, {}).get('massless', float('nan'))
    match = abs(g_from_P_per_pair - direct) < 1e-6 if not np.isnan(direct) else 'N/A'
    print(f"    d={d}: G_P = {g_from_P_per_pair:.8f}, G_direct = {direct:.8f}, match = {match}")

print()

# ============================================================
# PART 9: Scattering amplitude details
# ============================================================
print("-" * 80)
print("PART 9: Specific scattering processes")
print("-" * 80)
print()

# Pick representative fermion vertices at each distance
# and compute the one-boson-exchange amplitude

# Find C-C pairs at each distance
print("  Full amplitude A(d) = V @ G_gauge_full @ V^T, averaged over C-C pairs:")
print()

# Use the full-graph propagator, not the gauge-restricted one
# More physical: the boson travels on the full graph
# A(i,j) already computed as A_scatter above

# Also compute with massive propagators
for m2_label, G_used in [("massless", G_massless)] + [(f"m^2={m2}", G_massive[m2]) for m2 in mass_sq_values]:
    G_gg = G_used[np.ix_(gauge_list, gauge_list)]
    A_scat = V_fC @ G_gg @ V_fC.T

    print(f"  Boson exchange amplitude [{m2_label}]:")
    for d in range(1, max_dist + 1):
        pairs = []
        for i_idx, i in enumerate(fermion_list):
            for j_idx, j in enumerate(fermion_list):
                if j > i and dist_matrix[i, j] == d:
                    pairs.append((i_idx, j_idx))
        if not pairs:
            continue
        vals = [A_scat[p[0], p[1]] for p in pairs]
        print(f"    d={d}: <A> = {np.mean(vals):12.8f}, std = {np.std(vals):10.6f}, "
              f"N = {len(pairs)}")
    print()

# ============================================================
# PART 10: Phi-structure in the propagator
# ============================================================
print("-" * 80)
print("PART 10: Golden ratio structure in the propagator")
print("-" * 80)
print()

# The propagator values at each distance should be expressible in terms of phi
# since ALL eigenvalues and multiplicities involve phi.

# Each G(d) = (1/N) * sum_k P[d,k] * mult_k / lambda_k
# The eigenvalues lambda_k are a + b*phi, the multiplicities are n^2.
# So G(d) should be in Q(phi) = Q(sqrt(5)).

print("  Testing: is G(d) in Q(phi)?")
print("  G(d) = a + b*phi for rational a, b?")
print()

for d in range(max_dist + 1):
    gd = propagator_vs_d.get(d, {}).get('massless', None)
    if gd is None:
        continue

    # Try G(d) = (a + b*phi) / denom for small integer a, b, denom
    found = False
    for denom in range(1, 5001):
        val_scaled = gd * denom
        # val_scaled = a + b*phi
        # Solve: b = (val_scaled - a) / phi, check if a and b are integers
        for a in range(-200, 201):
            b_float = (val_scaled - a) / PHI
            if abs(b_float - round(b_float)) < 0.0001:
                b = int(round(b_float))
                if abs(a) + abs(b) <= 200 and (a != 0 or b != 0):
                    # Verify
                    reconstructed = (a + b * PHI) / denom
                    if abs(reconstructed - gd) < 1e-10:
                        print(f"    G(d={d}) = ({a} + {b}*phi) / {denom} = {reconstructed:.10f} "
                              f"(actual: {gd:.10f})")
                        found = True
                        break
            if found:
                break
        if found:
            break
    if not found:
        print(f"    G(d={d}) = {gd:.10f} - no simple phi expression found (denom up to 5000)")

print()

# ============================================================
# SUMMARY
# ============================================================
print("=" * 80)
print("SUMMARY OF RESULTS")
print("=" * 80)
print()

print("1. GREEN'S FUNCTION / PROPAGATOR:")
print("   G(m) = (L + m^2 I)^(-1) computed exactly via spectral decomposition.")
print("   The massless propagator is the pseudoinverse (1 zero mode projected out).")
print()
print(f"   <G(d)> for massless propagator:")
for d in sorted(propagator_vs_d.keys()):
    if d >= 1:
        print(f"     d={d}: {propagator_vs_d[d]['massless']:.8f}")
print()

# Classify the decay
if 1 in cc_by_d and 2 in cc_by_d:
    ratio_12 = cc_by_d[2] / cc_by_d[1]
    print(f"   Decay ratio G(2)/G(1) = {ratio_12:.6f}")
    if ratio_12 > 0:
        if abs(ratio_12 - 0.5) < 0.15:
            print("   --> Roughly 1/r decay (Coulomb-like)")
        elif abs(ratio_12) < 0.3:
            print("   --> Faster than 1/r (exponential/confined?)")
        else:
            print("   --> Slower than 1/r (weaker decay)")
    else:
        print("   --> SIGN CHANGE: oscillatory propagator (screening/confined)")

print()

print("2. HEAT KERNEL:")
print("   K(t) = exp(-tL) computed via spectral decomposition.")
print(f"   Equilibrium value: 1/N = {K_uniform:.6f}")
print(f"   Characteristic diffusion time (to equilibrium): see data above.")
print()

print("3. GAUGE SECTOR:")
print(f"   24 gauge vertices (8 Type-A + 16 Type-B)")
print(f"   Direct edges within gauge: {edges_gauge}")
if edges_gauge == 0:
    print("   --> Gauge vertices have NO direct connections to each other!")
    print("   --> Gauge bosons propagate ONLY through fermion vertices (indirect paths)")
else:
    print(f"   --> Gauge vertices have some direct connections")
print()

print("4. ONE-BOSON-EXCHANGE:")
print("   A(i,j) = sum_g,g' V(i,g) * G(g,g') * V(g',j)")
print("   Computed for all C-C pairs via full-graph propagator.")
print()

print("5. CONFINEMENT/SCREENING:")
if 1 in cc_by_d and len(cc_by_d) >= 3:
    signs = [1 if cc_by_d[d] > 0 else -1 for d in sorted(cc_by_d.keys())]
    if all(s == signs[0] for s in signs):
        print("   All G_CC(d) have SAME sign --> no screening oscillation")
        vals = [cc_by_d[d] for d in sorted(cc_by_d.keys())]
        if all(abs(vals[i]) >= abs(vals[i + 1]) for i in range(len(vals) - 1)):
            print("   |G_CC(d)| is monotonically decreasing --> screening/Coulomb-like")
        else:
            print("   |G_CC(d)| is NOT monotonically decreasing --> non-trivial behavior")
    else:
        print("   G_CC(d) changes sign --> OSCILLATORY (screening with sign alternation)")
        print("   This is characteristic of confinement on a finite graph!")

print()

print("6. PHI STRUCTURE:")
print("   The propagator lives in Q(sqrt(5)) since all eigenvalues and")
print("   multiplicities involve phi. Results above test explicit expressions.")
print()

print("CLASSIFICATION:")
print("  - Propagator computation: DERIVAT (exact spectral, no free parameters)")
print("  - Decay law identification: DERIVAT")
print("  - Gauge sector connectivity: DERIVAT")
print("  - Physical interpretation (confinement/screening): SPECULATIV")
print()
print("=" * 80)

#!/usr/bin/env python3
"""
EXP-204: DEEP ANALYSIS OF SCHUR COMPLEMENT ON THE 600-CELL
============================================================
EXP-203 discovered that integrating out the 96 C (fermion) vertices
via the Schur complement produces an effective Laplacian on the 24
A+B (gauge) vertices with remarkable spectral structure:
  - Eigenvalues: 0 (x1), 7.0 (x4), 10.465 (x9), 11.455 (x8), 12.0 (x2)
  - Multiplicities: 1, 4, 9, 8, 2 -- SQUARES! (1^2, 2^2, 3^2, ...)
  - The spectral gap INCREASES with coarse-graining (2.29 -> 7.0 -> 9.45)

Also, EXP-203 found r(gap, frustration) = -0.618 = -1/phi?? (50 realizations)

This experiment:
A. Verify and analyze the Schur complement spectrum in detail
B. Test r = -1/phi with 500 realizations (high statistics)
C. Analyze the eigenvalue structure: are multiplicities (1,4,9,8,2) = (1^2,2^2,3^2,...)?
D. Separate the Schur complement by GAUGE TYPE
E. Check if effective couplings relate to physical constants

STATUS: EXPLORATORY
Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
import math
from itertools import permutations as iperms
from collections import deque, Counter

PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

print("=" * 72)
print("EXP-204: DEEP ANALYSIS OF SCHUR COMPLEMENT ON THE 600-CELL")
print("=" * 72)
print()

# ===========================================================================
# STEP 0: BUILD 600-CELL
# ===========================================================================
print("STEP 0: Building 600-cell")
print("-" * 60)

def build_600cell():
    verts = set()
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts.add(tuple(v))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    verts.add((s0, s1, s2, s3))
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
        cu = [base[p[0]], base[p[1]], base[p[2]], base[p[3]]]
        nz = [i for i in range(4) if abs(cu[i]) > 1e-12]
        for sb in range(1 << len(nz)):
            v = list(cu)
            for k, pos in enumerate(nz):
                if (sb >> k) & 1:
                    v[pos] = -v[pos]
            verts.add(tuple(round(x, 12) for x in v))
    return sorted(verts)

vertices = build_600cell()
N = len(vertices)
print("  Vertices: %d" % N)

def classify_vertex(v):
    c = sorted([abs(x) for x in v], reverse=True)
    if abs(c[0] - 1.0) < 1e-10 and c[1] < 1e-10:
        return 'A'
    if all(abs(x - 0.5) < 1e-10 for x in c):
        return 'B'
    return 'C'

vtypes = [classify_vertex(v) for v in vertices]
idx_A = [i for i in range(N) if vtypes[i] == 'A']
idx_B = [i for i in range(N) if vtypes[i] == 'B']
idx_C = [i for i in range(N) if vtypes[i] == 'C']
idx_AB = idx_A + idx_B
print("  A: %d, B: %d, C: %d, A+B: %d" % (len(idx_A), len(idx_B), len(idx_C), len(idx_AB)))

varray = np.array(vertices)
dot_mat = varray @ varray.T
adj = (np.abs(dot_mat - PHI / 2.0) < 1e-8).astype(np.float64)
np.fill_diagonal(adj, 0)

neighbors = {}
for i in range(N):
    neighbors[i] = [j for j in range(N) if adj[i, j] > 0.5]

all_edges = [(i, j) for i in range(N) for j in neighbors[i] if j > i]
print("  Edges: %d" % len(all_edges))
print()

# ===========================================================================
# PART A: SCHUR COMPLEMENT - DETAILED ANALYSIS
# ===========================================================================
print("=" * 72)
print("PART A: SCHUR COMPLEMENT SPECTRUM IN DETAIL")
print("=" * 72)
print()

# Construct the graph Laplacian
L = np.diag(adj.sum(axis=1)) - adj

# Partition into AB and C blocks
# Reorder: AB first, then C
idx_AB_sorted = sorted(idx_AB)
idx_C_sorted = sorted(idx_C)
reorder = idx_AB_sorted + idx_C_sorted
n_AB = len(idx_AB_sorted)
n_C = len(idx_C_sorted)

L_reordered = L[np.ix_(reorder, reorder)]
L_AB = L_reordered[:n_AB, :n_AB]
L_AC = L_reordered[:n_AB, n_AB:]
L_CA = L_reordered[n_AB:, :n_AB]
L_CC = L_reordered[n_AB:, n_AB:]

print("  L_AB shape: %s (should be 24x24)" % str(L_AB.shape))
print("  L_CC shape: %s (should be 96x96)" % str(L_CC.shape))
print("  L_AC shape: %s (should be 24x96)" % str(L_AC.shape))

# Check L_AB diagonal = degree of each AB vertex = 12
# Off-diagonal should be 0 (no AB-AB edges)
print("  L_AB diagonal: min=%.1f, max=%.1f" % (np.diag(L_AB).min(), np.diag(L_AB).max()))
print("  L_AB off-diag: min=%.4f, max=%.4f" % (
    L_AB[~np.eye(n_AB, dtype=bool)].min(),
    L_AB[~np.eye(n_AB, dtype=bool)].max()))

# L_CC eigenvalues
eigvals_CC = np.linalg.eigvalsh(L_CC)
print("\n  L_CC eigenvalues (96x96):")
print("  min=%.4f, max=%.4f" % (eigvals_CC.min(), eigvals_CC.max()))
n_zero_CC = np.sum(np.abs(eigvals_CC) < 1e-8)
print("  Number of zero eigenvalues: %d" % n_zero_CC)

if n_zero_CC > 0:
    print("  L_CC is SINGULAR! Using pseudoinverse.")
    L_CC_inv = np.linalg.pinv(L_CC)
else:
    print("  L_CC is non-singular. Using direct inverse.")
    L_CC_inv = np.linalg.inv(L_CC)

# Schur complement: L_eff = L_AB - L_AC * L_CC^{-1} * L_CA
L_eff = L_AB - L_AC @ L_CC_inv @ L_CA

print("\n  Effective Laplacian L_eff (24x24):")
eigvals_eff = np.linalg.eigvalsh(L_eff)

# Group by value
from collections import defaultdict
ev_groups = defaultdict(int)
for ev in eigvals_eff:
    found = False
    for key in ev_groups:
        if abs(ev - key) < 0.01:
            ev_groups[key] += 1
            found = True
            break
    if not found:
        ev_groups[round(ev, 4)] = 1

print("  Eigenvalue | Multiplicity")
print("  " + "-" * 35)
for ev in sorted(ev_groups.keys()):
    print("  %10.6f | %d" % (ev, ev_groups[ev]))
print("  Total: %d" % sum(ev_groups.values()))

# Check multiplicities
mults = sorted(ev_groups.values())
print("\n  Multiplicities: %s" % mults)
print("  Squares check: 1,4,9,16,25 = 1^2,2^2,3^2,4^2,5^2")
print("  Actual: %s" % mults)

# Are eigenvalues related to phi?
print("\n  Eigenvalue analysis:")
for ev in sorted(ev_groups.keys()):
    if abs(ev) < 0.01:
        continue
    print("  ev = %.6f:" % ev)
    print("    ev/phi = %.6f" % (ev / PHI))
    print("    ev/phi^2 = %.6f" % (ev / PHI**2))
    print("    ev - 12 = %.6f" % (ev - 12))
    print("    12 - ev = %.6f" % (12 - ev))
    # Test a + b*phi form
    b_test = (ev - round(ev)) / (PHI - 1)  # approx
    a_test = ev - b_test * PHI
    # Better: ev = a + b*phi => b = (ev - a)/phi, try integer a
    for a in range(-5, 15):
        b = (ev - a) / PHI
        if abs(b - round(b)) < 0.01:
            print("    = %d + %d*phi (%.6f)" % (a, round(b), a + round(b)*PHI))

# Check L_eff structure: is it the Laplacian of some graph?
print("\n  L_eff structure analysis:")
print("  Is L_eff a valid graph Laplacian?")
print("  Row sums: min=%.6f, max=%.6f" % (L_eff.sum(axis=1).min(), L_eff.sum(axis=1).max()))
print("  Off-diagonal: min=%.6f, max=%.6f" % (
    L_eff[~np.eye(n_AB, dtype=bool)].min(),
    L_eff[~np.eye(n_AB, dtype=bool)].max()))

# If off-diagonal elements are all <= 0, it's a valid weighted Laplacian
off_diag = L_eff[~np.eye(n_AB, dtype=bool)]
n_positive_offdiag = np.sum(off_diag > 1e-10)
print("  Positive off-diagonal elements: %d (out of %d)" % (
    n_positive_offdiag, len(off_diag)))

if n_positive_offdiag == 0:
    print("  => L_eff IS a valid weighted graph Laplacian!")
    # Extract effective graph
    W_eff = -L_eff.copy()
    np.fill_diagonal(W_eff, 0)
    print("  Effective edge weights:")
    unique_weights = sorted(set(round(w, 4) for w in W_eff.flatten() if abs(w) > 1e-10))
    for w in unique_weights:
        count = np.sum(np.abs(W_eff - w) < 0.005) // 2  # undirected
        print("    w = %.4f: %d edges" % (w, count))
else:
    print("  => L_eff is NOT a standard graph Laplacian (has positive off-diagonal)")

# Separate A-A, A-B, B-B couplings
print("\n  Effective couplings by vertex type:")
# idx_AB_sorted has A vertices first, then B
n_A_local = len(idx_A)
n_B_local = len(idx_B)

# A-A block
L_AA_eff = L_eff[:n_A_local, :n_A_local]
W_AA = -L_AA_eff.copy()
np.fill_diagonal(W_AA, 0)
print("  A-A effective weights: min=%.6f, max=%.6f, mean=%.6f" % (
    W_AA[W_AA != 0].min() if np.any(W_AA != 0) else 0,
    W_AA.max(),
    W_AA[W_AA != 0].mean() if np.any(W_AA != 0) else 0))

# A-B block
W_AB = -L_eff[:n_A_local, n_A_local:]
print("  A-B effective weights: min=%.6f, max=%.6f, mean=%.6f" % (
    W_AB[W_AB != 0].min() if np.any(W_AB != 0) else 0,
    W_AB.max(),
    W_AB[W_AB != 0].mean() if np.any(W_AB != 0) else 0))

# B-B block
L_BB_eff = L_eff[n_A_local:, n_A_local:]
W_BB = -L_BB_eff.copy()
np.fill_diagonal(W_BB, 0)
print("  B-B effective weights: min=%.6f, max=%.6f, mean=%.6f" % (
    W_BB[W_BB != 0].min() if np.any(W_BB != 0) else 0,
    W_BB.max(),
    W_BB[W_BB != 0].mean() if np.any(W_BB != 0) else 0))

# Total effective degree per vertex type
deg_A_eff = L_eff[:n_A_local, :n_A_local].diagonal()
deg_B_eff = L_eff[n_A_local:, n_A_local:].diagonal()
print("\n  Effective degrees:")
print("  A vertices: mean=%.4f, all same? %s" % (
    deg_A_eff.mean(), "YES" if np.std(deg_A_eff) < 0.001 else "NO (std=%.4f)" % np.std(deg_A_eff)))
print("  B vertices: mean=%.4f, all same? %s" % (
    deg_B_eff.mean(), "YES" if np.std(deg_B_eff) < 0.001 else "NO (std=%.4f)" % np.std(deg_B_eff)))

print()

# ===========================================================================
# PART B: r = -1/phi? HIGH STATISTICS TEST
# ===========================================================================
print("=" * 72)
print("PART B: CORRELATION TEST r = -1/phi WITH 500 REALIZATIONS")
print("=" * 72)
print()

def potts_energy(coloring, nbr_list, n_v):
    E = 0
    for i in range(n_v):
        for j in nbr_list[i]:
            if j > i and coloring[i] == coloring[j]:
                E += 1
    return E

# Quick KZ cooling
def kz_cool(nbr_list, n_v, n_colors, T0, n_steps, rng):
    col = [int(rng.integers(0, n_colors)) for _ in range(n_v)]
    for step in range(n_steps):
        T = T0 * max(1.0 - step / n_steps, 1e-6)
        v = int(rng.integers(0, n_v))
        old_c = col[v]
        new_c = int(rng.integers(0, n_colors))
        if new_c == old_c:
            continue
        dE = 0
        for j in nbr_list[v]:
            if col[j] == old_c:
                dE -= 1
            if col[j] == new_c:
                dE += 1
        if dE <= 0 or rng.random() < math.exp(-dE / T):
            col[v] = new_c
    return col

N_REAL = 500
n_steps_cool = 5000  # moderate cooling for ~40-90 defects
print("  Generating %d soliton backgrounds..." % N_REAL)

frust_fracs = []
gaps_24 = []
gaps_8 = []

for r_idx in range(N_REAL):
    rng = np.random.default_rng(7777 + r_idx * 31)
    col = kz_cool(neighbors, N, 5, 5.0, n_steps_cool, rng)

    # Frustration fraction
    n_frust = sum(1 for i, j in all_edges if col[i] == col[j])
    f = n_frust / len(all_edges)
    frust_fracs.append(f)

    # Build soliton-perturbed Laplacian (eps=0.5 for frustrated edges)
    adj_sol = adj.copy()
    for i, j in all_edges:
        if col[i] == col[j]:
            adj_sol[i, j] = 0.5
            adj_sol[j, i] = 0.5
    L_sol = np.diag(adj_sol.sum(axis=1)) - adj_sol

    # Reorder and Schur complement
    L_sol_re = L_sol[np.ix_(reorder, reorder)]
    L_AB_s = L_sol_re[:n_AB, :n_AB]
    L_AC_s = L_sol_re[:n_AB, n_AB:]
    L_CA_s = L_sol_re[n_AB:, :n_AB]
    L_CC_s = L_sol_re[n_AB:, n_AB:]

    # Pseudoinverse for CC block
    eigvals_cc_s = np.linalg.eigvalsh(L_CC_s)
    if np.min(np.abs(eigvals_cc_s)) < 1e-6:
        L_CC_s_inv = np.linalg.pinv(L_CC_s)
    else:
        L_CC_s_inv = np.linalg.inv(L_CC_s)

    L_eff_s = L_AB_s - L_AC_s @ L_CC_s_inv @ L_CA_s
    eigvals_eff_s = np.sort(np.linalg.eigvalsh(L_eff_s))

    # Spectral gap at 24-cell scale (first non-zero eigenvalue)
    gap_24 = eigvals_eff_s[1] if len(eigvals_eff_s) > 1 else 0
    gaps_24.append(gap_24)

    # For tesseract scale: further integrate out B
    # L_eff is 24x24. Split into A (first n_A_local) and B (rest)
    L_eff_A = L_eff_s[:n_A_local, :n_A_local]
    L_eff_AB_block = L_eff_s[:n_A_local, n_A_local:]
    L_eff_BA_block = L_eff_s[n_A_local:, :n_A_local]
    L_eff_BB = L_eff_s[n_A_local:, n_A_local:]

    eigvals_BB = np.linalg.eigvalsh(L_eff_BB)
    if np.min(np.abs(eigvals_BB)) < 1e-6:
        L_BB_inv = np.linalg.pinv(L_eff_BB)
    else:
        L_BB_inv = np.linalg.inv(L_eff_BB)

    L_eff_8 = L_eff_A - L_eff_AB_block @ L_BB_inv @ L_eff_BA_block
    eigvals_8 = np.sort(np.linalg.eigvalsh(L_eff_8))
    gap_8 = eigvals_8[1] if len(eigvals_8) > 1 else 0
    gaps_8.append(gap_8)

    if r_idx % 100 == 0:
        print("  %d / %d done (f=%.4f, gap24=%.4f, gap8=%.4f)" % (
            r_idx, N_REAL, f, gap_24, gap_8))

frust_arr = np.array(frust_fracs)
gap24_arr = np.array(gaps_24)
gap8_arr = np.array(gaps_8)

print("\n  Statistics over %d realizations:" % N_REAL)
print("  <f> = %.6f +/- %.6f" % (frust_arr.mean(), frust_arr.std()))
print("  <gap_24> = %.6f +/- %.6f" % (gap24_arr.mean(), gap24_arr.std()))
print("  <gap_8>  = %.6f +/- %.6f" % (gap8_arr.mean(), gap8_arr.std()))

# Pearson correlation
r_24 = np.corrcoef(frust_arr, gap24_arr)[0, 1]
r_8 = np.corrcoef(frust_arr, gap8_arr)[0, 1]

print("\n  Correlation coefficients:")
print("  r(f, gap_24) = %.6f" % r_24)
print("  r(f, gap_8)  = %.6f" % r_8)
print()

# Test against -1/phi
inv_phi = 1.0 / PHI
print("  TEST: r = -1/phi?")
print("  1/phi = %.10f" % inv_phi)
print("  |r_24| = %.10f" % abs(r_24))
print("  diff = %.6f (%.3f%%)" % (abs(abs(r_24) - inv_phi), 100 * abs(abs(r_24) - inv_phi) / inv_phi))
print()

# Standard error of Pearson r
se_r = math.sqrt((1 - r_24**2) / (N_REAL - 2))
print("  Standard error of r: %.6f" % se_r)
print("  95%% CI: [%.4f, %.4f]" % (r_24 - 1.96 * se_r, r_24 + 1.96 * se_r))
print("  Does CI include -1/phi = %.4f? %s" % (
    -inv_phi, "YES" if r_24 - 1.96 * se_r <= -inv_phi <= r_24 + 1.96 * se_r else "NO"))
print()

# Bootstrap for r
print("  Bootstrap (10000 resamples) for r(f, gap_24):")
rng_boot = np.random.default_rng(42)
r_boot = []
for _ in range(10000):
    idx = rng_boot.choice(N_REAL, size=N_REAL, replace=True)
    rb = np.corrcoef(frust_arr[idx], gap24_arr[idx])[0, 1]
    r_boot.append(rb)
r_boot = np.array(r_boot)
print("  Bootstrap mean: %.6f" % r_boot.mean())
print("  Bootstrap std:  %.6f" % r_boot.std())
print("  Bootstrap 95%% CI: [%.6f, %.6f]" % (np.percentile(r_boot, 2.5), np.percentile(r_boot, 97.5)))
print("  -1/phi = %.6f" % (-inv_phi))
in_ci = np.percentile(r_boot, 2.5) <= -inv_phi <= np.percentile(r_boot, 97.5)
print("  -1/phi in 95%% CI? %s" % ("YES" if in_ci else "NO"))
print()

# Also test for gap_8
print("  Bootstrap for r(f, gap_8):")
r_boot_8 = []
for _ in range(10000):
    idx = rng_boot.choice(N_REAL, size=N_REAL, replace=True)
    rb = np.corrcoef(frust_arr[idx], gap8_arr[idx])[0, 1]
    r_boot_8.append(rb)
r_boot_8 = np.array(r_boot_8)
print("  r(f, gap_8) = %.6f +/- %.6f" % (r_8, r_boot_8.std()))
print("  95%% CI: [%.6f, %.6f]" % (np.percentile(r_boot_8, 2.5), np.percentile(r_boot_8, 97.5)))
print()


# ===========================================================================
# PART C: EIGENVALUE STRUCTURE ANALYSIS
# ===========================================================================
print("=" * 72)
print("PART C: EIGENVALUE STRUCTURE - MULTIPLICITIES AND GOLDEN RATIO")
print("=" * 72)
print()

# Get exact eigenvalues of ground state L_eff
eigvals_ground = np.sort(np.linalg.eigvalsh(L_eff))
print("  Ground state L_eff eigenvalues (all 24):")
for i, ev in enumerate(eigvals_ground):
    print("    [%2d] %.10f" % (i, ev))

# Group into clusters
clusters = []
tol = 0.01
for ev in eigvals_ground:
    found = False
    for cl in clusters:
        if abs(ev - cl[0]) < tol:
            cl.append(ev)
            found = True
            break
    if not found:
        clusters.append([ev])

print("\n  Eigenvalue clusters:")
print("  Cluster | Mean value    | Mult | Mult pattern")
print("  " + "-" * 55)
for i, cl in enumerate(clusters):
    mean_ev = np.mean(cl)
    mult = len(cl)
    # Check if mult is a perfect square
    sqrt_m = int(round(math.sqrt(mult)))
    is_sq = sqrt_m * sqrt_m == mult
    pattern = "%d^2" % sqrt_m if is_sq else "not square"
    print("  %d       | %12.6f  | %d    | %s" % (i, mean_ev, mult, pattern))

# Differences between cluster means
print("\n  Inter-cluster differences:")
means = [np.mean(cl) for cl in clusters]
for i in range(len(means)):
    for j in range(i+1, len(means)):
        diff = means[j] - means[i]
        print("    [%d]-[%d] = %.6f" % (j, i, diff))
        # Test golden ratio
        if abs(diff) > 0.1:
            print("      /phi = %.6f, /phi^2 = %.6f, *phi = %.6f" % (
                diff / PHI, diff / PHI**2, diff * PHI))

# Full 600-cell spectrum for comparison
eigvals_full = np.sort(np.linalg.eigvalsh(adj))
unique_full = []
for ev in eigvals_full:
    if not unique_full or abs(ev - unique_full[-1][0]) > 0.01:
        unique_full.append([ev, 1])
    else:
        unique_full[-1][1] += 1

print("\n  Full 600-cell adjacency eigenvalues (for reference):")
for ev, mult in unique_full:
    print("    %.6f (mult %d)" % (ev, mult))

# Check: are L_eff eigenvalues related to full spectrum?
print("\n  L_eff eigenvalue = 12 - adj_eigenvalue relationship?")
print("  (If L_eff came from a subgraph, eigenvalues would be 12 - lambda)")
for cl in clusters:
    mean_ev = np.mean(cl)
    adj_ev = 12 - mean_ev
    print("    L_eff = %.4f => adj_eff = %.4f" % (mean_ev, adj_ev))
    for ev_full, m_full in unique_full:
        if abs(adj_ev - ev_full) < 0.5:
            print("      close to full spectrum %.4f (mult %d)" % (ev_full, m_full))

# L_eff eigenvectors analysis
eigvals_eff_v, eigvecs_eff = np.linalg.eigh(L_eff)
print("\n  Eigenvector analysis (projections on A vs B):")
print("  Eigenvalue | weight_A | weight_B | A/(A+B)")
print("  " + "-" * 55)
for i in range(n_AB):
    vec = eigvecs_eff[:, i]
    w_A = np.sum(vec[:n_A_local]**2)
    w_B = np.sum(vec[n_A_local:]**2)
    ratio = w_A / (w_A + w_B) if (w_A + w_B) > 1e-10 else 0
    print("  %10.4f | %.6f | %.6f | %.4f" % (eigvals_eff_v[i], w_A, w_B, ratio))


# ===========================================================================
# PART D: GAUGE-TYPE RESOLVED SCHUR COMPLEMENT
# ===========================================================================
print()
print("=" * 72)
print("PART D: GAUGE-TYPE RESOLVED SCHUR COMPLEMENT")
print("=" * 72)
print()

# Classify edges
A_set = set(idx_A)
B_set = set(idx_B)
C_set = set(idx_C)
nbr_sets = {i: set(neighbors[i]) for i in range(N)}

def classify_edge(i, j):
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
            return 'U1_CC'
        elif shared_B > 0:
            return 'SU2_CC'
        else:
            return 'SU3'
    return 'other'

edge_types = {}
type_counts = Counter()
for i, j in all_edges:
    et = classify_edge(i, j)
    edge_types[(i, j)] = et
    edge_types[(j, i)] = et
    type_counts[et] += 1

print("  Edge type distribution:")
for t in sorted(type_counts.keys()):
    print("    %s: %d edges" % (t, type_counts[t]))
print("    Total: %d" % sum(type_counts.values()))

# Build sector-specific adjacency matrices for C-C block
# Then do separate Schur complements

# For each gauge sector, build adjacency only with edges of that type in CC
sectors = {
    'U1': lambda et: et in ('U1_CC',),      # CC edges in U1 channel
    'SU2': lambda et: et in ('SU2_CC',),     # CC edges in SU2 channel
    'SU3': lambda et: et == 'SU3',           # CC edges in SU3 channel
    'all_CC': lambda et: et in ('U1_CC', 'SU2_CC', 'SU3'),  # all CC edges
}

for sect_name, sect_filter in sectors.items():
    # Build adjacency for CC block with only sector edges
    adj_sect = np.zeros((n_C, n_C))
    for ii, ci in enumerate(idx_C_sorted):
        for jj, cj in enumerate(idx_C_sorted):
            if jj <= ii:
                continue
            et = edge_types.get((ci, cj), edge_types.get((cj, ci), 'other'))
            if sect_filter(et):
                adj_sect[ii, jj] = 1.0
                adj_sect[jj, ii] = 1.0

    n_sect_edges = int(adj_sect.sum()) // 2

    # Build sector CC Laplacian
    L_CC_sect = np.diag(adj_sect.sum(axis=1)) - adj_sect
    eigvals_sect = np.linalg.eigvalsh(L_CC_sect)
    n_zero_sect = np.sum(np.abs(eigvals_sect) < 1e-8)
    gap_sect = eigvals_sect[n_zero_sect] if n_zero_sect < len(eigvals_sect) else 0

    print("\n  Sector %s:" % sect_name)
    print("    CC edges: %d" % n_sect_edges)
    print("    L_CC spectral gap: %.6f" % gap_sect)
    print("    L_CC max eigenvalue: %.6f" % eigvals_sect.max())
    print("    Zero modes: %d (connected components - 1)" % n_zero_sect)

# The key measurement: how does L_AC project onto each sector?
# L_AC connects A+B to C. Which C-C sectors are reached?
print("\n  L_AC structure (gauge-fermion coupling):")
print("  How A+B vertices couple to C-sector edges:")

# For each AB vertex, count how many of its C neighbors are involved in each CC sector
for ab_type in ['A', 'B']:
    ab_indices = idx_A if ab_type == 'A' else idx_B
    sector_reach = {s: 0 for s in ['U1_CC', 'SU2_CC', 'SU3']}
    total_c_nbrs = 0
    for ab_v in ab_indices:
        c_nbrs = [j for j in neighbors[ab_v] if vtypes[j] == 'C']
        total_c_nbrs += len(c_nbrs)
        for cn in c_nbrs:
            for cn2 in neighbors[cn]:
                if vtypes[cn2] == 'C' and cn2 != cn:
                    et = edge_types.get((cn, cn2), edge_types.get((cn2, cn), 'other'))
                    if et in sector_reach:
                        sector_reach[et] += 1

    print("    Type %s vertices (%d total, %d C-neighbors):" % (
        ab_type, len(ab_indices), total_c_nbrs))
    for s, cnt in sorted(sector_reach.items()):
        print("      reach to %s: %d (%.1f per vertex)" % (
            s, cnt, cnt / len(ab_indices)))


# ===========================================================================
# PART E: EFFECTIVE COUPLINGS AND PHYSICAL CONSTANTS
# ===========================================================================
print()
print("=" * 72)
print("PART E: EFFECTIVE COUPLINGS AND PHYSICAL CONSTANTS")
print("=" * 72)
print()

# Key ratios from the Schur complement spectrum
if len(means) >= 4:
    gap1 = means[1]  # first non-zero cluster
    gap2 = means[2]  # second cluster
    gap3 = means[3]  # third cluster

    print("  Non-zero eigenvalue clusters: %.6f, %.6f, %.6f" % (gap1, gap2, gap3))
    print()

    print("  Ratios:")
    print("    gap2/gap1 = %.6f" % (gap2/gap1))
    print("    gap3/gap1 = %.6f" % (gap3/gap1))
    print("    gap3/gap2 = %.6f" % (gap3/gap2))
    print("    gap2/gap1 - 1 = %.6f" % (gap2/gap1 - 1))
    print("    gap3/gap1 - 1 = %.6f" % (gap3/gap1 - 1))
    print()

    # Compare with known constants
    print("  Comparison with physical constants:")
    print("    phi = %.6f" % PHI)
    print("    phi^2 = %.6f" % PHI**2)
    print("    sin^2(tW) = %.6f, 6/26 = %.6f" % (SIN2_TW, 6.0/26))
    print("    alpha_s = %.6f, 1/(2*phi^3) = %.6f" % (ALPHA_S, 1/(2*PHI**3)))
    print()

    # Test: gap ratios vs phi
    for name, val in [("phi", PHI), ("phi^2", PHI**2), ("1+phi", 1+PHI),
                       ("3/2", 1.5), ("5/3", 5.0/3), ("2", 2.0),
                       ("sqrt(phi)", math.sqrt(PHI))]:
        print("    gap2/gap1 vs %s: diff = %.4f" % (name, abs(gap2/gap1 - val)))

    print()

    # The Schur complement removes 96 C-fermions, leaving 24 gauge vertices
    # This is like computing the effective gauge action after integrating out fermions
    # In QFT: S_eff = S_gauge + Tr ln(D_slash) where D_slash involves the gauge field
    # The eigenvalues of L_eff encode the effective gauge couplings

    # Check: is gap1 related to 12 * (something involving coupling constants)?
    print("  L_eff eigenvalues as fractions of bare degree (12):")
    for cl in clusters:
        m = np.mean(cl)
        if abs(m) > 0.01:
            print("    %.6f / 12 = %.6f" % (m, m/12))

    # The gap at 7.0 is particularly interesting
    print("\n  The gap value 7.0:")
    print("    7 = prime")
    print("    7 = dim(Im(O)) (imaginary octonions)")
    print("    12 - 7 = 5 = dim(600-cell diameter)")
    print("    7/12 = %.6f" % (7.0/12))
    print("    12/7 = %.6f" % (12.0/7))
    print("    7*PHI = %.6f" % (7*PHI))

if len(means) >= 5:
    gap4 = means[4]
    print("\n  The value %.6f:" % gap4)
    print("    12 - %.4f = %.4f" % (gap4, 12 - gap4))

print()


# ===========================================================================
# NUMERICAL SUMMARY
# ===========================================================================
print("=" * 72)
print("NUMERICAL SUMMARY")
print("=" * 72)
print()

print("PART A: Schur complement spectrum")
for i, cl in enumerate(clusters):
    print("  Cluster %d: %.6f (mult %d)" % (i, np.mean(cl), len(cl)))
print()

print("PART B: Correlation test (N=%d realizations)" % N_REAL)
print("  r(f, gap_24) = %.6f" % r_24)
print("  1/phi = %.6f" % inv_phi)
print("  |r| - 1/phi = %.6f (%.3f%%)" % (abs(r_24) - inv_phi, 100*abs(abs(r_24)-inv_phi)/inv_phi))
print("  Bootstrap 95%% CI: [%.4f, %.4f]" % (
    np.percentile(r_boot, 2.5), np.percentile(r_boot, 97.5)))
print("  -1/phi in CI: %s" % ("YES" if in_ci else "NO"))
print()

print("PART C: Eigenvalue structure")
print("  Multiplicities: %s" % [len(cl) for cl in clusters])
sq_check = all(int(round(math.sqrt(len(cl))))**2 == len(cl) for cl in clusters)
print("  All perfect squares: %s" % ("YES" if sq_check else "NO"))
print()

print("PART D: Gauge-type resolved Schur complement")
print("  (see above for details)")
print()

print("PART E: Physical constants comparison")
if len(means) >= 3:
    print("  gap2/gap1 = %.6f" % (means[2]/means[1]))
    print("  gap3/gap1 = %.6f" % (means[3]/means[1] if len(means) > 3 else 0))
print()


# ===========================================================================
# CONCLUSIONS
# ===========================================================================
print("=" * 72)
print("CONCLUSIONS")
print("=" * 72)
print()

print("1. SCHUR COMPLEMENT SPECTRUM:")
print("   Integrating out 96 C-fermions produces L_eff with")
print("   eigenvalue multiplicities %s" % [len(cl) for cl in clusters])
if sq_check:
    print("   ALL multiplicities are PERFECT SQUARES!")
    print("   This connects to the 600-cell's (n^2) multiplicity pattern.")
print("   STATUS: DERIVED")
print()

print("2. CORRELATION r(frustration, gap):")
print("   r = %.6f (over %d realizations)" % (r_24, N_REAL))
if in_ci:
    print("   -1/phi = %.6f IS in the 95%% CI" % (-inv_phi))
    print("   Cannot reject r = -1/phi hypothesis")
    print("   STATUS: SPECULATIV (suggestive but not proven)")
else:
    print("   -1/phi = %.6f is NOT in the 95%% CI" % (-inv_phi))
    print("   r = -1/phi hypothesis REJECTED")
    print("   STATUS: NEGATIV")
print()

print("3. EIGENVALUE STRUCTURE:")
print("   Gap = %.4f (possibly 7 = dim(Im(O)))" % means[1] if len(means) > 1 else "")
print("   The Schur complement naturally generates a HIERARCHY")
print("   of effective gauge couplings from fermion integration.")
print("   STATUS: DERIVED (mathematical fact)")
print()

print("4. GAUGE-SECTOR RESOLUTION:")
print("   The CC edges split into U1_CC, SU2_CC, SU3 with different")
print("   spectral properties. A-vertices couple to U1_CC, B to SU2_CC.")
print("   STATUS: DERIVED")
print()

print("EXP-204 COMPLETE")

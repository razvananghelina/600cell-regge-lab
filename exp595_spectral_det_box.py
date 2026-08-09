#!/usr/bin/env python3
"""
EXP-595: SPECTRAL DETERMINANT AL BOX — Iese 57?
==================================================

Paper: log det'(L/N) / ln(alpha) = 56.97 (near-miss, NEGATIVE cu L).
Dar L nu e Box. Box are structura Hopf, Galois, gauge/KK.

INTREBARE: log det'(Box) / ln(alpha) = 57? Sau vreo varianta?

Calculam cu TOATE 120 eigenvalues exacte (nu medii per irrep).
"""
import numpy as np
from numpy.linalg import eigh
import sys
sys.path.insert(0, ".")
from commons import build_600cell
from commons.constants import PHI, SQRT5, a1, b1, N, alpha_em

print("=" * 80)
print("EXP-595: SPECTRAL DETERMINANT AL BOX")
print("=" * 80)

# Build 600-cell + Hopf + Box
verts, adj, lap = build_600cell()

def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

fibers = None
for i in range(N):
    if abs(verts[i, 0] - PHI/2) < 1e-6:
        g = verts[i]; p = g.copy(); ok = True
        for k in range(2, 11):
            p = qmul(p, g)
            if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok=False; break
            if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok=False
        if not ok: continue
        used = set(); fibers_tmp = []; subg = []
        pp = np.array([1.0,0,0,0])
        for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
        for s in range(N):
            if s in used: continue
            fib = []
            for si in subg:
                q = qmul(verts[s], verts[si]); idx = find_idx(q, verts)
                if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
            if len(fib) == 10: fibers_tmp.append(fib)
        if len(fibers_tmp) == 12:
            fibers = fibers_tmp
            break

A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.0

Box = b1 * A_fiber - adj

# Get ALL 120 eigenvalues of Box, L, and Box^2
evals_Box = np.sort(np.linalg.eigvalsh(Box))
evals_L = np.sort(np.linalg.eigvalsh(lap))

ln_alpha = np.log(alpha_em)  # negative! ln(1/137) < 0

# =====================================================================
# PART 1: Eigenvalue catalog
# =====================================================================
print("\n" + "=" * 80)
print("PART 1: EIGENVALUE CATALOG")
print("=" * 80)

# Box eigenvalues — group by value
from collections import Counter
box_rounded = Counter(np.round(evals_Box, 4))
print(f"\n  Box eigenvalues (distinct, with multiplicities):")
for val, mult in sorted(box_rounded.items()):
    print(f"    lambda = {val:>10.4f}  (x{mult})")

n_zero_box = sum(1 for e in evals_Box if abs(e) < 0.001)
n_nonzero_box = N - n_zero_box
print(f"\n  Zero eigenvalues: {n_zero_box}")
print(f"  Nonzero eigenvalues: {n_nonzero_box}")

# L eigenvalues
L_rounded = Counter(np.round(evals_L, 4))
print(f"\n  Laplacian eigenvalues (distinct):")
for val, mult in sorted(L_rounded.items()):
    print(f"    lambda = {val:>10.4f}  (x{mult})")

n_zero_L = sum(1 for e in evals_L if abs(e) < 0.001)

# =====================================================================
# PART 2: Spectral determinants
# =====================================================================
print("\n" + "=" * 80)
print("PART 2: SPECTRAL DETERMINANTS")
print("=" * 80)

# det'(M) = product of nonzero eigenvalues
# log det'(M) = sum of log|lambda_k| for nonzero lambda_k

# --- Box ---
nonzero_box = [e for e in evals_Box if abs(e) > 0.001]
log_det_box = sum(np.log(abs(e)) for e in nonzero_box)
log_det_box2 = sum(np.log(e**2) for e in nonzero_box)  # = 2 * log_det_|box|

# --- L ---
nonzero_L = [e for e in evals_L if abs(e) > 0.001]
log_det_L = sum(np.log(e) for e in nonzero_L)  # L eigenvalues are all positive

# --- Box^2 ---
log_det_box_sq = log_det_box2  # = sum log(lambda^2) = 2*sum log|lambda|

print(f"  log det'(|Box|) = {log_det_box:.10f}")
print(f"  log det'(Box^2) = {log_det_box2:.10f}")
print(f"  log det'(L)     = {log_det_L:.10f}")
print(f"  ln(alpha)       = {ln_alpha:.10f}")

# =====================================================================
# PART 3: Exponent search — divide by ln(alpha)
# =====================================================================
print("\n" + "=" * 80)
print("PART 3: EXPONENT z = log det' / ln(alpha)")
print("=" * 80)

target = 57 - 1/(2*PHI**3)  # = 56.882
print(f"  Target exponent: 57 - alpha_s = {target:.4f}")
print(f"  Near target: 57.0000")
print()

# Try various normalizations
quantities = {
    "log det'(|Box|)": log_det_box,
    "log det'(Box^2)": log_det_box2,
    "log det'(L)":     log_det_L,
    "log det'(|Box|) - log det'(L)": log_det_box - log_det_L,
    "log det'(L) - log det'(|Box|)": log_det_L - log_det_box,
    "log det'(|Box|) / 2": log_det_box / 2,
    "log det'(L) / 2": log_det_L / 2,
}

# Also try normalized versions
quantities["log det'(L/N)"] = log_det_L - n_nonzero_box * np.log(N)  # hmm wrong
quantities["log det'(L) - 119*ln(N)"] = log_det_L - (N-1)*np.log(N)
quantities["log det'(L) - 119*ln(12)"] = log_det_L - (N-1)*np.log(12)

# The paper's version: log det'(L/N) = log det'(L) - (N-1)*log(N)
log_det_L_over_N = log_det_L - (N-n_zero_L)*np.log(N)
quantities["log det'(L/N) [paper]"] = log_det_L_over_N

# Box normalized
log_det_box_over_N = log_det_box - n_nonzero_box*np.log(N)
quantities["log det'(|Box|/N)"] = log_det_box_over_N

# Box over various normalizations
for norm_name, norm_val in [
    ("a1", a1), ("b1", b1), ("degree", 12), ("N/b1", N/b1),
    ("phi^4", PHI**4), ("2*pi", 2*np.pi),
]:
    quantities[f"log det'(|Box|/{norm_name})"] = log_det_box - n_nonzero_box*np.log(norm_val)

print(f"  {'Quantity':>45s} {'Value':>12s} {'/ ln(a)':>12s} {'vs 57':>8s}")
print(f"  {'-'*82}")
for name, val in quantities.items():
    z = val / ln_alpha
    diff = z - 57
    flag = " ***" if abs(diff) < 0.5 else " **" if abs(diff) < 2 else ""
    print(f"  {name:>45s} {val:>12.4f} {z:>12.4f} {diff:>+8.4f}{flag}")

# =====================================================================
# PART 4: Galois-physical sector only
# =====================================================================
print("\n" + "=" * 80)
print("PART 4: SPECTRAL DETERMINANT — DOAR SECTORUL FIZIC")
print("=" * 80)

# Physical sector: rho_0, rho_1, rho_2, rho_3, rho_4, rho_5
# Dark sector: rho_6, rho_7, rho_8

# Identify eigenvalues per sector using simultaneous diag
M_combo = lap + np.e * Box
evals_M, evecs_M = eigh(M_combo)
L_vals = np.array([evecs_M[:,i] @ lap @ evecs_M[:,i] for i in range(N)])
mu_vals = np.array([evecs_M[:,i] @ Box @ evecs_M[:,i] for i in range(N)])

L_known = {
    0: 0.0, 1: 12-6*PHI, 2: 10-2*SQRT5, 3: 9.0,
    4: 12.0, 5: 14.0, 6: 10+2*SQRT5, 7: 15.0, 8: 6+6*PHI
}

phys_irreps = {0, 1, 2, 3, 4, 5}
dark_irreps = {6, 7, 8}

phys_evals = []
dark_evals = []

for idx in range(N):
    # Which irrep?
    best_rho = min(range(9), key=lambda r: abs(L_vals[idx] - L_known[r]))
    mu = mu_vals[idx]
    if best_rho in phys_irreps:
        phys_evals.append(mu)
    else:
        dark_evals.append(mu)

n_phys = len(phys_evals)
n_dark = len(dark_evals)
print(f"  Physical sector modes: {n_phys}")
print(f"  Dark sector modes: {n_dark}")
print(f"  Total: {n_phys + n_dark}")

# Multiplicities check
phys_mults = {0:1, 1:4, 2:9, 3:16, 4:25, 5:36}
dark_mults = {6:9, 7:16, 8:4}
print(f"  Expected physical: {sum(phys_mults.values())} = 1+4+9+16+25+36 = 91")
print(f"  Expected dark: {sum(dark_mults.values())} = 9+16+4 = 29")

# Spectral determinant per sector
phys_nonzero = [e for e in phys_evals if abs(e) > 0.001]
dark_nonzero = [e for e in dark_evals if abs(e) > 0.001]

log_det_phys = sum(np.log(abs(e)) for e in phys_nonzero)
log_det_dark = sum(np.log(abs(e)) for e in dark_nonzero)

print(f"\n  log det'(|Box_phys|) = {log_det_phys:.6f}  ({len(phys_nonzero)} nonzero modes)")
print(f"  log det'(|Box_dark|) = {log_det_dark:.6f}  ({len(dark_nonzero)} nonzero modes)")
print(f"  Difference = {log_det_phys - log_det_dark:.6f}")
print(f"  Sum = {log_det_phys + log_det_dark:.6f} (= log det'(|Box|) = {log_det_box:.6f})")

# Divide by ln(alpha)
z_phys = log_det_phys / ln_alpha
z_dark = log_det_dark / ln_alpha
z_diff = (log_det_phys - log_det_dark) / ln_alpha
z_ratio = log_det_phys / log_det_dark if abs(log_det_dark) > 0.01 else 0

print(f"\n  Exponents (/ ln(alpha)):")
print(f"    Physical: {z_phys:.4f}")
print(f"    Dark:     {z_dark:.4f}")
print(f"    Diff:     {z_diff:.4f}  (vs 57: {z_diff-57:+.4f})")
print(f"    Ratio phys/dark: {z_ratio:.6f}")

# =====================================================================
# PART 5: Gauge sector determinant
# =====================================================================
print("\n" + "=" * 80)
print("PART 5: SPECTRAL DETERMINANT — GAUGE SECTOR")
print("=" * 80)

# Project Box to gauge sector
P0 = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            P0[i, j] = 1.0 / len(fib)

Box_gauge = P0 @ Box @ P0
evals_gauge = np.sort(np.linalg.eigvalsh(Box_gauge))
gauge_nonzero = [e for e in evals_gauge if abs(e) > 0.001]

log_det_gauge = sum(np.log(abs(e)) for e in gauge_nonzero)
z_gauge = log_det_gauge / ln_alpha

# KK sector
P_KK = np.eye(N) - P0
Box_KK = P_KK @ Box @ P_KK
evals_KK = np.sort(np.linalg.eigvalsh(Box_KK))
KK_nonzero = [e for e in evals_KK if abs(e) > 0.001]
log_det_KK = sum(np.log(abs(e)) for e in KK_nonzero)
z_KK = log_det_KK / ln_alpha

print(f"  Gauge sector: {len(gauge_nonzero)} nonzero, log det' = {log_det_gauge:.6f}, z = {z_gauge:.4f}")
print(f"  KK sector:    {len(KK_nonzero)} nonzero, log det' = {log_det_KK:.6f}, z = {z_KK:.4f}")
print(f"  Gauge z vs 57: {z_gauge - 57:+.4f}")

# =====================================================================
# PART 6: Combinatii de sectoare
# =====================================================================
print("\n" + "=" * 80)
print("PART 6: COMBINATII INTRE SECTOARE")
print("=" * 80)

# Try: z_CC = some combination of z_phys, z_dark, z_gauge, z_KK
combos = {
    "z_phys": z_phys,
    "z_dark": z_dark,
    "z_gauge": z_gauge,
    "z_KK": z_KK,
    "z_phys - z_dark": z_phys - z_dark,
    "z_phys + z_dark": z_phys + z_dark,
    "z_gauge + z_KK": z_gauge + z_KK,
    "z_gauge - z_KK": z_gauge - z_KK,
    "z_phys - z_gauge": z_phys - z_gauge,
    "z_phys / a1": z_phys / a1,
    "z_dark * a1": z_dark * a1,
    "(z_phys+z_dark)/2": (z_phys+z_dark)/2,
    "z_phys - z_dark - z_gauge": z_phys - z_dark - z_gauge,
    "z_KK / a1": z_KK / a1,
    "z_gauge * (a1-1)": z_gauge * (a1-1),
    "-z_gauge": -z_gauge,
    "-z_KK": -z_KK,
    "-z_phys": -z_phys,
}

print(f"\n  {'Combination':>35s} {'Value':>12s} {'vs 57':>10s}")
print(f"  {'-'*60}")
for name, val in sorted(combos.items(), key=lambda x: abs(x[1]-57)):
    diff = val - 57
    flag = " <-- CLOSE!" if abs(diff) < 1 else ""
    print(f"  {name:>35s} {val:>12.4f} {diff:>+10.4f}{flag}")

# =====================================================================
# PART 7: Paper's near-miss recalculated
# =====================================================================
print("\n" + "=" * 80)
print("PART 7: RECALCULARE NEAR-MISS DIN PAPER")
print("=" * 80)

# Paper claims: log det'(L/N) / ln(alpha) ~ 56.97
# Let me compute this exactly
# det'(L/N) = det'(L) / N^(N-1) since L has one zero eigenvalue
# log det'(L/N) = log det'(L) - (N-1)*log(N)
log_det_L_normalized = log_det_L - (N-1)*np.log(N)
z_paper = log_det_L_normalized / ln_alpha

print(f"  log det'(L) = {log_det_L:.10f}")
print(f"  (N-1)*ln(N) = {(N-1)*np.log(N):.10f}")
print(f"  log det'(L/N) = {log_det_L_normalized:.10f}")
print(f"  z = log det'(L/N) / ln(alpha) = {z_paper:.6f}")
print(f"  Paper claims: 56.97")
print(f"  Target: 57.0000")
print(f"  Difference: {z_paper - 57:+.6f}")

# Now try same with Box
# Box has 9 zero eigenvalues (kernel)
n_zero_box_actual = sum(1 for e in evals_Box if abs(e) < 0.001)
log_det_box_normalized = log_det_box - n_nonzero_box*np.log(N)
z_box_normalized = log_det_box_normalized / ln_alpha

print(f"\n  Same normalization with Box:")
print(f"  log det'(|Box|/N) = {log_det_box_normalized:.10f}")
print(f"  z = {z_box_normalized:.6f}")
print(f"  Difference from 57: {z_box_normalized - 57:+.6f}")

# =====================================================================
# PART 8: Galois-invariant determinant
# =====================================================================
print("\n" + "=" * 80)
print("PART 8: GALOIS-INVARIANT SPECTRAL DETERMINANT")
print("=" * 80)

# The Galois-invariant quantity: det'(Box) * det'(sigma(Box))
# Since Box eigenvalues split into Galois pairs + fixed:
# Fixed: rho_0(0), rho_3(9), rho_4(12), rho_5(14)... wait, these are L eigenvalues
# The Galois-invariant determinant uses N(lambda) = lambda * sigma(lambda)
# For Galois pairs: product is rational (Galois norm)

# Group eigenvalues by Galois type
# Physical (non-dark) eigenvalues: all those from rho_0..rho_5
# Dark eigenvalues: rho_6, rho_7, rho_8
# The Galois norm: product of physical_eigenvalue * dark_partner_eigenvalue

# Actually for the spectral determinant:
# The Galois-invariant log det should be:
# (1/2) * (log det'_phys + log det'_dark) = log det'_full / 2 (for paired modes)
# Plus the fixed modes

# Let me compute det'(Box^2) which is automatically Galois-invariant
# since lambda^2 is always real and Galois-invariant for paired eigenvalues
log_det_box_sq = 2 * log_det_box  # = sum 2*log|lambda| = sum log(lambda^2)
z_box_sq = log_det_box_sq / ln_alpha

print(f"  log det'(Box^2) = {log_det_box_sq:.6f}")
print(f"  z = log det'(Box^2) / ln(alpha) = {z_box_sq:.4f}")
print(f"  z/2 = {z_box_sq/2:.4f}")

# Try Box^2 normalized
log_det_box_sq_norm = log_det_box_sq - n_nonzero_box*np.log(N**2)
z_box_sq_norm = log_det_box_sq_norm / ln_alpha

print(f"\n  log det'(Box^2/N^2) = {log_det_box_sq_norm:.6f}")
print(f"  z = {z_box_sq_norm:.4f}")
print(f"  vs 57: {z_box_sq_norm - 57:+.4f}")

# =====================================================================
# VERDICT
# =====================================================================
print("\n" + "=" * 80)
print("V E R D I C T")
print("=" * 80)

# Find the closest to 57 among all quantities tested
all_results = []
for name, val in quantities.items():
    z = val / ln_alpha
    all_results.append((name, z, abs(z-57)))
all_results.extend([
    ("z_phys - z_dark", z_phys - z_dark, abs(z_phys-z_dark-57)),
    ("z_gauge", z_gauge, abs(z_gauge-57)),
    ("z_paper (L/N)", z_paper, abs(z_paper-57)),
    ("z_box_sq_norm", z_box_sq_norm, abs(z_box_sq_norm-57)),
])

all_results.sort(key=lambda x: x[2])
print(f"\n  Top 5 closest to 57:")
for name, z, diff in all_results[:5]:
    print(f"    {name:>45s}: z = {z:.6f}, |z-57| = {diff:.6f}")

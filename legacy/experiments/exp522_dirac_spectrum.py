"""
EXP-522: Dirac Operator on the Poincaré Homology Sphere S^3/2I
================================================================
The spectral action Tr(f(D/Lambda)) uses the DIRAC operator, not the
graph Laplacian. The Dirac spectrum on S^3/2I is computed from the
McKay correspondence: decompose SU(2) reps restricted to 2I.

Dirac eigenvalues on S^3: ±(k+3/2), k=0,1,2,...
Multiplicity on S^3/2I: n_{rho_1, k+1} (times fundamental of 2I
appears in V_{k+1}|_2I, for spinor equivariance).

McKay recursion:
  n_{i, k+1} = sum_j A_{ij} * n_{j, k} - n_{i, k-1}
  where A = adjacency matrix of extended E8 Dynkin diagram.
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, SQRT5, a1, b1, N, h, degree, d_ST, rank_E8, dim_E8,
                      alpha_em, inv_alpha, sin2_tW, Omega_DM, N_gen)

print("=" * 70)
print("EXP-522: DIRAC SPECTRUM ON S^3/2I (POINCARÉ HOMOLOGY SPHERE)")
print("=" * 70)

# ============================================================
# STEP 1: McKAY ADJACENCY MATRIX (EXTENDED E8)
# ============================================================
print(f"\n{'='*70}")
print("STEP 1: McKAY GRAPH = EXTENDED E8 DYNKIN DIAGRAM")
print(f"{'='*70}")

# 9 irreps of 2I: dims 1,2,3,4,5,6,4,2,3
irrep_dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
n_irreps = len(irrep_dims)

# McKay graph: rho_1 ⊗ rho_i = sum_j A_{ij} rho_j
# Graph: 0-1-2-3-4-5-6-7, plus 5-8 (branch at rho_5)
A = np.zeros((n_irreps, n_irreps), dtype=int)
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
for i, j in edges:
    A[i,j] = 1
    A[j,i] = 1

print(f"  Irrep dims: {irrep_dims}")
print(f"  McKay graph (extended E8):")
print(f"    rho0(1)-rho1(2)-rho2(3)-rho3(4)-rho4(5)-rho5(6)-rho6(4)-rho7(2)")
print(f"                                              |")
print(f"                                            rho8(3)")

# Verify: 2*dim(rho_i) = sum_j A_{ij} * dim(rho_j)
print(f"\n  McKay verification: 2*dim(i) = sum_j A_ij*dim(j)?")
for i in range(n_irreps):
    lhs = 2 * irrep_dims[i]
    rhs = sum(A[i,j] * irrep_dims[j] for j in range(n_irreps))
    ok = "OK" if lhs == rhs else "FAIL"
    print(f"    rho_{i}: 2*{irrep_dims[i]} = {lhs}, sum = {rhs}  {ok}")

# ============================================================
# STEP 2: McKAY RECURSION - DECOMPOSE V_n|_{2I}
# ============================================================
print(f"\n{'='*70}")
print("STEP 2: DECOMPOSITION OF SU(2) REPS RESTRICTED TO 2I")
print(f"{'='*70}")

N_max = 250  # compute V_1 through V_{N_max}

# n[i][k] = multiplicity of rho_i in V_k|_{2I}
# V_k is the k-dimensional SU(2) representation (spin (k-1)/2)
decomp = np.zeros((n_irreps, N_max + 1), dtype=int)

# Initial conditions
decomp[0, 1] = 1  # V_1 = trivial = rho_0
decomp[1, 2] = 1  # V_2 = fundamental = rho_1

# Recursion: n_{i, k+1} = sum_j A_{ij} * n_{j, k} - n_{i, k-1}
for k in range(2, N_max + 1):
    for i in range(n_irreps):
        decomp[i, k] = sum(A[i,j] * decomp[j, k-1] for j in range(n_irreps)) - decomp[i, k-2]
        # Multiplicities should be non-negative
        if decomp[i, k] < 0:
            print(f"  WARNING: negative multiplicity n_{i},{k} = {decomp[i,k]}")

# Verify dimensions: sum_i n_{i,k} * dim(rho_i) = k
print(f"  Dimension check (first 15):")
for k in range(1, 16):
    dim_check = sum(decomp[i, k] * irrep_dims[i] for i in range(n_irreps))
    reps = [f"{decomp[i,k]}*rho{i}" for i in range(n_irreps) if decomp[i,k] > 0]
    print(f"    V_{k:>3} = {' + '.join(reps):>40}  dim={dim_check:>3} (expect {k})")

# ============================================================
# STEP 3: DIRAC SPECTRUM ON S^3/2I
# ============================================================
print(f"\n{'='*70}")
print("STEP 3: DIRAC EIGENVALUES ON POINCARÉ HOMOLOGY SPHERE")
print(f"{'='*70}")

# Eigenvalue +(k+3/2) has multiplicity = n_{1, k+1}
# (times rho_1 appears in V_{k+1}|_{2I}, for spinor equivariance)

# Also eigenvalue -(k+3/2) with possibly different multiplicity
# For negative eigenvalues: m_k^- = n_{1, k+2} (shifted by 1)
# Actually, the formula for negative eigenvalues might be:
# m_k^- = n_{1, k+1} as well (same formula, symmetric spectrum)
# OR it could involve a different irrep.

# Let's first compute the positive spectrum.
print(f"\n  Positive Dirac eigenvalues (first appearances of rho_1):")
print(f"  {'k':>4} {'lambda':>10} {'m_k=n(1,k+1)':>14} {'V_{k+1} content':>40}")

dirac_pos = []  # (eigenvalue, multiplicity)
for k in range(N_max - 1):
    m_k = decomp[1, k + 1]  # n_{rho_1, k+1}
    if m_k > 0:
        lam = k + 1.5
        dirac_pos.append((lam, m_k))
        if len(dirac_pos) <= 25:
            reps = [f"{decomp[i,k+1]}*rho{i}" for i in range(n_irreps) if decomp[i,k+1] > 0]
            print(f"  {k:>4} {lam:>10.1f} {m_k:>14} {(' + '.join(reps))[:40]:>40}")

print(f"\n  Total nonzero eigenvalues found: {len(dirac_pos)} (out of {N_max-1})")
print(f"  First eigenvalue: {dirac_pos[0][0]}")
print(f"  First 10 eigenvalues: {[e[0] for e in dirac_pos[:10]]}")
print(f"  First 10 multiplicities: {[e[1] for e in dirac_pos[:10]]}")

# ============================================================
# STEP 4: PATTERN ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("STEP 4: PATTERN IN DIRAC SPECTRUM")
print(f"{'='*70}")

# The eigenvalues where rho_1 appears
k_values = [k for k in range(N_max-1) if decomp[1, k+1] > 0]

# Differences between consecutive k values
if len(k_values) > 1:
    diffs = [k_values[i+1] - k_values[i] for i in range(len(k_values)-1)]
    print(f"  k values where m_k > 0: {k_values[:20]}")
    print(f"  Differences: {diffs[:20]}")
    print(f"  h (Coxeter number) = {h}")
    print(f"  Period? {set(diffs[:20]) if len(diffs)>0 else 'N/A'}")

# The actual eigenvalues
eigenvalues_pos = [k + 1.5 for k in k_values]
print(f"\n  Eigenvalues: {eigenvalues_pos[:15]}")

# Check for framework structure
print(f"\n  First eigenvalue: {dirac_pos[0][0]}")
print(f"    = {dirac_pos[0][0]:.1f}")
if len(dirac_pos) > 1:
    print(f"  Second eigenvalue: {dirac_pos[1][0]}")
    print(f"  Ratio: {dirac_pos[1][0]/dirac_pos[0][0]:.6f}")
    if len(dirac_pos) > 2:
        print(f"  Third eigenvalue: {dirac_pos[2][0]}")

# ============================================================
# STEP 5: GALOIS STRUCTURE OF DIRAC SPECTRUM
# ============================================================
print(f"\n{'='*70}")
print("STEP 5: GALOIS STRUCTURE - WHICH IRREPS ARE BROKEN?")
print(f"{'='*70}")

# The Galois-broken irreps of 2I are: rho_1, rho_2, rho_7, rho_8
# rho_1 <-> rho_7 (both dim 2, characters involve phi)
# rho_2 <-> rho_8 (both dim 3, characters involve sqrt5)

galois_broken = {1, 2, 7, 8}
galois_pairs = [(1, 7), (2, 8)]

# For the Dirac spectrum, we used rho_1. What about rho_7?
print(f"  Comparing rho_1 vs rho_7 (Galois partners):")
print(f"  {'k':>4} {'n(rho1,k)':>10} {'n(rho7,k)':>10} {'diff':>6}")
for k in range(1, min(40, N_max+1)):
    n1 = decomp[1, k]
    n7 = decomp[7, k]
    if n1 > 0 or n7 > 0:
        print(f"  {k:>4} {n1:>10} {n7:>10} {n1-n7:>6}")

# Also check rho_2 vs rho_8
print(f"\n  Comparing rho_2 vs rho_8 (Galois partners):")
print(f"  {'k':>4} {'n(rho2,k)':>10} {'n(rho8,k)':>10} {'diff':>6}")
for k in range(1, min(40, N_max+1)):
    n2 = decomp[2, k]
    n8 = decomp[8, k]
    if n2 > 0 or n8 > 0:
        print(f"  {k:>4} {n2:>10} {n8:>10} {n2-n8:>6}")

# ============================================================
# STEP 6: SPECTRAL ZETA OF THE DIRAC OPERATOR
# ============================================================
print(f"\n{'='*70}")
print("STEP 6: DIRAC SPECTRAL ZETA FUNCTION")
print(f"{'='*70}")

# zeta_D(s) = sum_{lambda>0} m(lambda) * lambda^{-s}
# Only include eigenvalues where m > 0

def dirac_zeta(s, spectrum):
    """Compute sum m_k * lambda_k^{-s}."""
    return sum(m * lam**(-s) for lam, m in spectrum)

# Also define zeta for |D| (both signs): 2 * zeta_D(s) if spectrum is symmetric

print(f"  Using {len(dirac_pos)} positive eigenvalues up to lambda = {dirac_pos[-1][0]:.1f}")
print(f"  Total modes: {sum(m for _,m in dirac_pos)}")

print(f"\n  {'s':>8} {'zeta_D(s)':>14} {'zeta_D^2(s)':>14}")
for s in [1, 2, 3, 4, 5, PHI, 1/PHI, PHI**2, 2*PHI, np.pi]:
    zd = dirac_zeta(s, dirac_pos)
    zd2 = dirac_zeta(s/2, dirac_pos)  # zeta of D^2
    print(f"  {s:>8.4f} {zd:>14.6f} {zd2:>14.6f}")

# ============================================================
# STEP 7: SEARCH FOR FRAMEWORK CONSTANTS
# ============================================================
print(f"\n{'='*70}")
print("STEP 7: SEARCHING DIRAC ZETA FOR FRAMEWORK CONSTANTS")
print(f"{'='*70}")

targets = {
    '1/alpha': inv_alpha,
    'alpha': alpha_em,
    '4a1*phi^4': 4*a1*PHI**4,
    'sin2tW': sin2_tW,
    '7-phi': 7-PHI,
    'b1': float(b1),
    'a1': float(a1),
    'degree': float(degree),
    'phi': PHI,
    'phi^2': PHI**2,
    'phi^3': PHI**3,
    'phi^4': PHI**4,
    '2sqrt5': 2*SQRT5,
    'sqrt5': SQRT5,
    '2': 2.0,
    'pi': np.pi,
    'N': float(N),
    'h': float(h),
}

best = {}
for s_test in np.arange(0.1, 12.01, 0.001):
    zd = dirac_zeta(s_test, dirac_pos)
    if zd < 1e-15:
        continue
    for tname, tval in targets.items():
        if tval < 1e-6:
            continue
        err = abs(zd / tval - 1)
        if err < 0.005:
            key = tname
            if key not in best or err < best[key][0]:
                best[key] = (err, s_test, zd, tval)

print(f"\n  {'target':>12} {'s*':>10} {'zeta_D(s*)':>14} {'value':>12} {'err%':>10}")
print(f"  {'-'*12} {'-'*10} {'-'*14} {'-'*12} {'-'*10}")
for tname in sorted(best, key=lambda k: best[k][0]):
    err, s, zd, tv = best[tname]
    print(f"  {tname:>12} {s:>10.4f} {zd:>14.6f} {tv:>12.6f} {err*100:>10.4f}%")

# ============================================================
# STEP 8: HEAT KERNEL OF DIRAC OPERATOR
# ============================================================
print(f"\n{'='*70}")
print("STEP 8: DIRAC HEAT KERNEL Tr(exp(-t*D^2))")
print(f"{'='*70}")

def dirac_heat(t, spectrum):
    """Tr(exp(-t*D^2)) = sum m_k * exp(-t * lambda_k^2)."""
    return sum(m * np.exp(-t * lam**2) for lam, m in spectrum)

# The spectral action f(D^2/Lambda^2) at scale Lambda evaluates the heat kernel
# Tr(f(D^2/Lambda^2)) ~ sum f_n * a_n(D^2) * Lambda^{3-n}

# For the cutoff function f(x) = Theta(1-x), f(D^2/Lambda^2) counts eigenvalues < Lambda^2
# = N(Lambda) = number of eigenvalues with |lambda| < Lambda

print(f"  Counting function N(Lambda) = eigenvalues below Lambda:")
print(f"  {'Lambda':>10} {'N(Lambda)':>10} {'N/Lambda^3':>12}")
for Lambda in [5, 10, 20, 30, 50, 100, 150]:
    count = sum(m for lam, m in dirac_pos if lam < Lambda)
    vol_ratio = count / Lambda**3 if Lambda > 0 else 0
    print(f"  {Lambda:>10} {count:>10} {vol_ratio:>12.6f}")

# The asymptotic formula: N(Lambda) ~ vol(S^3/2I)/(6*pi^2) * Lambda^3
# vol(S^3/2I) = vol(S^3)/|2I| = 2*pi^2/120 = pi^2/60
expected_coeff = np.pi**2 / (6 * np.pi**2 * 60)  # = 1/(360)
print(f"\n  Expected coefficient: 1/360 = {1/360:.6f}")
print(f"  (From vol(S^3/2I)/(6*pi^2) = (2*pi^2/120)/(6*pi^2) = 1/360)")

# ============================================================
# STEP 9: SEELEY-DEWITT COEFFICIENTS FROM DIRAC
# ============================================================
print(f"\n{'='*70}")
print("STEP 9: SPECTRAL ACTION COEFFICIENTS")
print(f"{'='*70}")

# In the spectral action Tr(f(D/Lambda)):
# The heat kernel expansion gives:
# Tr(exp(-t*D^2)) ~ (4*pi*t)^{-3/2} * (a_0 + a_2*t + a_4*t^2 + ...)
# for t -> 0+

# a_0 = vol(M) / (4*pi)^{3/2} ... proportional to volume
# a_2 involves scalar curvature R
# a_4 involves R^2, Ricci^2, Weyl^2

# For S^3/2I with radius r:
# vol = 2*pi^2*r^3/120
# R = 6/r^2 (scalar curvature of S^3 of radius r)

# We can extract the coefficients numerically from the spectrum:
# a_0 = lim_{t->0} (4*pi*t)^{3/2} * Tr(exp(-t*D^2))

# Using the spectrum, check small-t behavior
print(f"  Small-t expansion of (4*pi*t)^(3/2) * Tr(exp(-t*D^2)):")
print(f"  {'t':>12} {'Tr':>14} {'(4pi*t)^3/2*Tr':>18} {'expected a0':>14}")

# Include both + and - eigenvalues (double the positive spectrum)
# Actually for Tr(exp(-t*D^2)), both signs contribute the same since D^2 is the same
# So: Tr(exp(-t*D^2)) = 2 * sum_positive m_k * exp(-t*lambda_k^2)

a0_expected = 2 * np.pi**2 / (120 * (4*np.pi)**1.5)  # rough estimate
for t in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]:
    tr = 2 * dirac_heat(t, dirac_pos)  # factor 2 for both signs
    coeff = (4*np.pi*t)**1.5 * tr
    print(f"  {t:>12.4f} {tr:>14.4f} {coeff:>18.6f} {a0_expected:>14.6f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  DIRAC SPECTRUM ON S^3/2I (POINCARÉ HOMOLOGY SPHERE):

  Computed via McKay recursion on extended E8 Dynkin diagram.
  Eigenvalue +(k+3/2) has multiplicity n_{{rho_1, k+1}}.

  First eigenvalue: {dirac_pos[0][0]} (mult {dirac_pos[0][1]})
  Spectrum gap structure: eigenvalues at k = {k_values[:10]}
  {'Periodic with period ' + str(diffs[0]) if len(set(diffs[:5]))==1 else 'Pattern: ' + str(diffs[:10])}

  This is the CORRECT operator for the spectral action.
  The graph Laplacian (exp520-521) was the wrong object.
""")

print("=" * 70)
print("EXP-522 COMPLETE")
print("=" * 70)

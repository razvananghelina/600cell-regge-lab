"""
exp338b: One-Loop VEV Profile - Single Parameter Search
========================================================
From exp338 we learned:
- CG selection rule: only "even" irreps couple to fermion densities
- Lepton corrections = 0 automatically (rho_0 uniform density)
- Overlap matrix rank = 8 with all modes

Now: construct the ONE-PARAMETER VEV profile from the Higgs propagator
on the 600-cell and see if it reproduces the mass corrections.

Model: delta_k = T3_k * alpha_s * Nc * F_k(r)
where F_k(r) = CG-weighted sum over Laplacian propagators
and r = kappa/mu^2 is the SINGLE free parameter.
"""
import numpy as np
from itertools import product, permutations

PHI = (1 + np.sqrt(5)) / 2
ALPHA_S = 1 / (2 * PHI**3)

# ================================================================
# CG table (verified in exp338)
# ================================================================
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

# CG_self[k, m] = multiplicity of rho_m in rho_k x rho_k
CG_self = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, 0],  # rho_0^2
    [1, 0, 1, 0, 0, 0, 0, 0, 0],  # rho_1^2
    [1, 0, 1, 0, 1, 0, 0, 0, 0],  # rho_2^2
    [1, 0, 1, 0, 1, 0, 1, 0, 1],  # rho_3^2
    [1, 0, 1, 0, 2, 0, 2, 0, 1],  # rho_4^2
    [1, 0, 2, 0, 3, 0, 2, 0, 2],  # rho_5^2
    [1, 0, 1, 0, 1, 0, 1, 0, 1],  # rho_6^2
    [1, 0, 0, 0, 0, 0, 0, 0, 1],  # rho_7^2
    [1, 0, 0, 0, 1, 0, 0, 0, 1],  # rho_8^2
])

# Cayley Laplacian eigenvalues (verified)
char_C1 = [1, PHI, PHI, 1, 0, -1, -1, -1/PHI, -1/PHI]
cayley_eigs = np.array([12 * (1 - char_C1[k] / dims[k]) for k in range(9)])

# Only "even" irreps appear in CG: rho_0, rho_2, rho_4, rho_6, rho_8
even_irreps = [0, 2, 4, 6, 8]
even_eigs = [cayley_eigs[m] for m in even_irreps]
even_dims = [dims[m] for m in even_irreps]

print("="*72)
print("exp338b: One-Loop VEV Profile - Single Parameter")
print("="*72)
print(f"\nalpha_s = 1/(2*phi^3) = {ALPHA_S:.6f}")
print(f"\nCayley eigenvalues for 'even' irreps:")
for m in even_irreps:
    print(f"  rho_{m}: lambda = {cayley_eigs[m]:.4f}, dim = {dims[m]}")

# ================================================================
# MODEL A: CG-weighted propagator sum
# ================================================================
print("\n" + "="*72)
print("MODEL A: CG-Weighted Propagator (Analytical)")
print("="*72)

# F_k(r) = sum_{m in even, m!=0} C^m_{kk} * d_m / (r * lambda_m + 2)
# delta_k = scale * T3_k * F_k(r)

fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
fermion_T3 = [-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5]
fermion_Nc = [0, 0, 0, 1, 1, 1, 1, 1, 1]  # 0 for leptons (no QCD)
delta_exp = np.array([0.000, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278])

def compute_Fk(r):
    """CG-weighted propagator sum for each fermion."""
    F = np.zeros(9)
    for k in range(9):
        for m in even_irreps:
            if m == 0:
                continue  # skip trivial (this is the uniform part)
            F[k] += CG_self[k, m] * dims[m] / (r * cayley_eigs[m] + 2)
    return F

# Scan r
print(f"\n{'r':>8} {'scale':>8} {'RMS':>8} | predicted deltas")
print("-"*90)

best_rms = 999
best_r = 0
best_scale = 0

for r in np.concatenate([np.linspace(0.001, 1, 200), np.linspace(1, 20, 200)]):
    F = compute_Fk(r)

    # Model: delta_k = scale * T3_k * Nc_k * F_k
    # For leptons: Nc = 0 => delta = 0
    # For quarks: Nc = 1 (relative)
    pred_raw = np.array([fermion_T3[k] * fermion_Nc[k] * F[k] for k in range(9)])

    # Find optimal scale
    if np.sum(pred_raw**2) > 1e-20:
        scale = np.sum(delta_exp * pred_raw) / np.sum(pred_raw**2)
        pred = scale * pred_raw
        rms = np.sqrt(np.mean((delta_exp - pred)**2))

        if rms < best_rms:
            best_rms = rms
            best_r = r
            best_scale = scale

# Print best
r = best_r
F = compute_Fk(r)
pred_raw = np.array([fermion_T3[k] * fermion_Nc[k] * F[k] for k in range(9)])
scale = best_scale
pred = scale * pred_raw

print(f"\nBest fit: r = {best_r:.4f}, scale = {best_scale:.4f}, RMS = {best_rms:.4f}")
print(f"\n{'fermion':>8} {'T3':>5} {'Nc':>4} {'F_k':>8} {'pred':>8} {'exp':>8} {'err':>8}")
print("-"*55)
for k in range(9):
    print(f"{fermion_names[k]:>8} {fermion_T3[k]:>+5.1f} {fermion_Nc[k]:>4} {F[k]:>8.4f} "
          f"{pred[k]:>+8.4f} {delta_exp[k]:>+8.4f} {pred[k]-delta_exp[k]:>+8.4f}")


# ================================================================
# MODEL B: T3 * Nc * F_k + generation correction
# ================================================================
print("\n" + "="*72)
print("MODEL B: T3*Nc*F + generation-dependent correction")
print("="*72)

# Generation index
gen = [0, 1, 2, 0, 1, 2, 0, 1, 2]

# Model: delta_k = a * T3 * Nc * F_k(r) + b * Nc * g_gen(k) + c * Nc
# where g_gen encodes generation dependence from CG

# The generation dependence comes from: within each T3 sector,
# heavier generations have more CG channels.
# Let's define g_k = sum_m C^m_{kk} (total CG multiplicity - 1)
g_total = np.array([sum(CG_self[k, m] for m in even_irreps if m != 0) for k in range(9)])

print(f"\nTotal non-trivial CG multiplicity per fermion:")
for k in range(9):
    print(f"  rho_{k} ({fermion_names[k]}): {g_total[k]}")

# Scan r with 3-parameter fit: a*T3*Nc*F + b*Nc + c
best_rms_B = 999
best_r_B = 0

for r in np.concatenate([np.linspace(0.001, 1, 200), np.linspace(1, 20, 200)]):
    F = compute_Fk(r)

    # Feature matrix
    X = np.column_stack([
        np.array([fermion_T3[k] * fermion_Nc[k] * F[k] for k in range(9)]),
        np.array([fermion_Nc[k] * g_total[k] for k in range(9)]),
        np.ones(9),
    ])

    coeffs, _, _, _ = np.linalg.lstsq(X, delta_exp, rcond=None)
    pred = X @ coeffs
    rms = np.sqrt(np.mean((delta_exp - pred)**2))

    if rms < best_rms_B:
        best_rms_B = rms
        best_r_B = r
        best_coeffs_B = coeffs

r = best_r_B
F = compute_Fk(r)
X = np.column_stack([
    np.array([fermion_T3[k] * fermion_Nc[k] * F[k] for k in range(9)]),
    np.array([fermion_Nc[k] * g_total[k] for k in range(9)]),
    np.ones(9),
])
pred = X @ best_coeffs_B

print(f"\nBest fit: r={best_r_B:.4f}, coeffs={np.round(best_coeffs_B, 4)}, RMS={best_rms_B:.4f}")
for k in range(9):
    print(f"  {fermion_names[k]:>4}: pred={pred[k]:>+.4f}, exp={delta_exp[k]:>+.4f}, err={pred[k]-delta_exp[k]:>+.4f}")


# ================================================================
# MODEL C: Direct CG weights without T3 assumption
# ================================================================
print("\n" + "="*72)
print("MODEL C: CG-propagator scan without T3 assumption")
print("="*72)

# For each fermion, compute the FULL CG-weighted propagator
# without assuming the T3 sign. Let the sign emerge from the fit.

# delta_k = sum_m c_m * C^m_{kk} * d_m / (r * lambda_m + 2)
# where c_m are mode-specific amplitudes (4 free params: m=2,4,6,8)

best_rms_C = 999
best_r_C = 0

for r in np.concatenate([np.linspace(0.001, 1, 200), np.linspace(1, 50, 200)]):
    # Feature matrix: each even mode (m=2,4,6,8) is a feature
    X_C = np.zeros((9, 4))
    for idx, m in enumerate([2, 4, 6, 8]):
        lam = cayley_eigs[m]
        prop = 1.0 / (r * lam + 2)
        for k in range(9):
            X_C[k, idx] = CG_self[k, m] * dims[m] * prop

    coeffs_C, _, _, _ = np.linalg.lstsq(X_C, delta_exp, rcond=None)
    pred_C = X_C @ coeffs_C
    rms_C = np.sqrt(np.mean((delta_exp - pred_C)**2))

    if rms_C < best_rms_C:
        best_rms_C = rms_C
        best_r_C = r
        best_coeffs_C = coeffs_C.copy()

r = best_r_C
X_C = np.zeros((9, 4))
for idx, m in enumerate([2, 4, 6, 8]):
    lam = cayley_eigs[m]
    prop = 1.0 / (r * lam + 2)
    for k in range(9):
        X_C[k, idx] = CG_self[k, m] * dims[m] * prop

pred_C = X_C @ best_coeffs_C

print(f"\nBest fit: r={best_r_C:.4f}, RMS={best_rms_C:.4f}")
print(f"Mode amplitudes (c_2, c_4, c_6, c_8): {np.round(best_coeffs_C, 4)}")
print(f"\n{'fermion':>8} {'pred':>8} {'exp':>8} {'err':>8}")
print("-"*40)
for k in range(9):
    print(f"{fermion_names[k]:>8} {pred_C[k]:>+8.4f} {delta_exp[k]:>+8.4f} {pred_C[k]-delta_exp[k]:>+8.4f}")

# What is the rank of X_C?
rank_C = np.linalg.matrix_rank(X_C, tol=1e-8)
print(f"\nRank of CG feature matrix: {rank_C}")
print("(4 modes for 9 fermions => 4 free params, but some rows are identical)")

# Check: which fermions have identical CG patterns?
print("\nCG patterns (for even m=2,4,6,8):")
for k in range(9):
    pattern = tuple(CG_self[k, m] for m in [2, 4, 6, 8])
    print(f"  rho_{k} ({fermion_names[k]}): {pattern}")

# Note: rho_3 and rho_6 have IDENTICAL CG patterns => same prediction
# rho_0 has all zeros => always predicts 0
# rho_1 and rho_2 differ (rho_1: (1,0,0,0), rho_2: (1,1,0,0))


# ================================================================
# MODEL D: PHYSICAL one-loop model with alpha_s
# ================================================================
print("\n" + "="*72)
print("MODEL D: Physical One-Loop with alpha_s")
print("="*72)

# Physical picture:
# The QCD correction to quark mass on the 600-cell is:
#
# delta m_q / m_q = (alpha_s / pi) * C_F * sum_m D_m(q)
#
# where D_m(q) depends on the CG coupling of quark q to mode m
# and the gluon propagator at the Cayley eigenvalue lambda_m.
#
# C_F = 4/3, alpha_s = 1/(2*phi^3)
#
# The correction to the EXPONENT is:
# delta_n = delta m / (m * ln(phi)) = (alpha_s * C_F / (pi * ln(phi))) * sum_m D_m

# alpha_s * C_F / (pi * ln(phi))
prefactor = ALPHA_S * (4.0/3) / (np.pi * np.log(PHI))
print(f"alpha_s * C_F / (pi * ln(phi)) = {prefactor:.6f}")

# For T3-dependent correction:
# delta_n = 2 * T3 * prefactor * F_k(r)
# The factor of 2 comes from the SU(2)_L structure

# But we need to account for the fact that the Higgs doublet
# couples differently to up and down quarks.
# Specifically, the VEV is in the T3 = -1/2 component,
# so down-type quarks couple directly, while up-type couple
# through the conjugate tilde{Phi}.
#
# The one-loop vertex correction has opposite sign for the two.

print(f"\nModel: delta_k = 2*T3*Nc*prefactor * F_k(r)")
print(f"  where F_k(r) = sum_m C^m_kk * d_m / (r*lambda_m + 2)")

# Find optimal r
best_rms_D = 999
for r in np.logspace(-3, 2, 2000):
    F = compute_Fk(r)
    pred = np.array([2 * fermion_T3[k] * fermion_Nc[k] * prefactor * F[k] for k in range(9)])

    # Allow ONE overall scale parameter
    if np.sum(pred**2) > 1e-20:
        s = np.sum(delta_exp * pred) / np.sum(pred**2)
        pred_s = s * pred
        rms = np.sqrt(np.mean((delta_exp - pred_s)**2))
        if rms < best_rms_D:
            best_rms_D = rms
            best_r_D = r
            best_s_D = s

r = best_r_D
F = compute_Fk(r)
pred = np.array([2 * fermion_T3[k] * fermion_Nc[k] * prefactor * F[k] for k in range(9)])
pred_s = best_s_D * pred
print(f"\nBest r = {best_r_D:.4f}, overall scale = {best_s_D:.4f}")
print(f"RMS = {best_rms_D:.4f}")
print(f"\n{'fermion':>8} {'pred':>8} {'exp':>8} {'err':>8}")
print("-"*40)
for k in range(9):
    print(f"{fermion_names[k]:>8} {pred_s[k]:>+8.4f} {delta_exp[k]:>+8.4f} {pred_s[k]-delta_exp[k]:>+8.4f}")


# ================================================================
# MODEL E: Allow different assignments of fermions to irreps
# ================================================================
print("\n" + "="*72)
print("MODEL E: Scan ALL fermion-irrep assignments")
print("="*72)

# The McKay assignment might not be optimal. Try all permutations.
# But 9! = 362880 is too many. Constrain:
# - e must go to rho_0 (only irrep with delta=0 automatically)
# - leptons go to {rho_0, x, y} where x,y from remaining
# - quarks go to remaining

# Actually, let's use Model C (4-param CG fit) and scan assignments.
# For each permutation of fermion-to-irrep:

# To reduce search: fix e -> rho_0 (forced by delta_e = 0)
# and try all 8! = 40320 permutations for the rest.
# This is feasible.

print("\nScanning 8! = 40320 assignments (e fixed at rho_0)...")

remaining_fermions = list(range(1, 9))  # mu, tau, u, c, t, d, s, b
remaining_irreps = list(range(1, 9))    # rho_1 through rho_8

best_rms_E = 999
best_perm_E = None
count = 0

for perm in permutations(remaining_irreps):
    count += 1
    # Assignment: fermion i -> irrep perm[i-1]
    assignment = [0] + list(perm)  # e -> rho_0, then permutation

    # Build CG feature matrix for this assignment
    X_E = np.zeros((9, 4))
    for idx, m in enumerate([2, 4, 6, 8]):
        for f in range(9):
            k = assignment[f]
            X_E[f, idx] = CG_self[k, m] * dims[m]

    # Fit
    coeffs_E, _, _, _ = np.linalg.lstsq(X_E, delta_exp, rcond=None)
    pred_E = X_E @ coeffs_E
    rms_E = np.sqrt(np.mean((delta_exp - pred_E)**2))

    if rms_E < best_rms_E:
        best_rms_E = rms_E
        best_perm_E = assignment
        best_pred_E = pred_E.copy()
        best_coeffs_E = coeffs_E.copy()

print(f"\nBest assignment found:")
print(f"  RMS = {best_rms_E:.4f}")
print(f"  Assignment: ", end="")
for f in range(9):
    print(f"{fermion_names[f]}->rho_{best_perm_E[f]} ", end="")
print()

print(f"\n  Mode amplitudes: {np.round(best_coeffs_E, 4)}")
print(f"\n  {'fermion':>8} {'irrep':>6} {'pred':>8} {'exp':>8} {'err':>8}")
print(f"  {'-'*45}")
for f in range(9):
    k = best_perm_E[f]
    print(f"  {fermion_names[f]:>8} rho_{k:>2} {best_pred_E[f]:>+8.4f} {delta_exp[f]:>+8.4f} {best_pred_E[f]-delta_exp[f]:>+8.4f}")

# CG patterns for best assignment
print(f"\n  CG patterns in best assignment:")
for f in range(9):
    k = best_perm_E[f]
    pattern = tuple(CG_self[k, m] for m in [2, 4, 6, 8])
    print(f"  {fermion_names[f]:>4} -> rho_{k} (dim {dims[k]}): CG = {pattern}")


# ================================================================
# MODEL F: T3-signed model with assignment scan
# ================================================================
print("\n" + "="*72)
print("MODEL F: T3-signed CG model + assignment scan")
print("="*72)

# Model: delta_f = T3_f * Nc_f * sum_m c_m * C^m_{kk} * d_m
# Scan over assignments AND mode amplitudes

best_rms_F = 999
best_perm_F = None

for perm in permutations(remaining_irreps):
    assignment = [0] + list(perm)

    # Build feature matrix with T3 * Nc weighting
    X_F = np.zeros((9, 4))
    for idx, m_idx in enumerate([2, 4, 6, 8]):
        for f in range(9):
            k = assignment[f]
            X_F[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]

    coeffs_F, _, _, _ = np.linalg.lstsq(X_F, delta_exp, rcond=None)
    pred_F = X_F @ coeffs_F
    rms_F = np.sqrt(np.mean((delta_exp - pred_F)**2))

    if rms_F < best_rms_F:
        best_rms_F = rms_F
        best_perm_F = assignment
        best_pred_F = pred_F.copy()
        best_coeffs_F = coeffs_F.copy()

print(f"\nBest T3-signed assignment:")
print(f"  RMS = {best_rms_F:.4f}")
print(f"  Assignment: ", end="")
for f in range(9):
    print(f"{fermion_names[f]}->rho_{best_perm_F[f]} ", end="")
print()
print(f"  Mode amplitudes: {np.round(best_coeffs_F, 4)}")
print(f"\n  {'fermion':>8} {'irrep':>6} {'T3':>5} {'Nc':>4} {'pred':>8} {'exp':>8} {'err':>8}")
print(f"  {'-'*55}")
for f in range(9):
    k = best_perm_F[f]
    print(f"  {fermion_names[f]:>8} rho_{k:>2} {fermion_T3[f]:>+5.1f} {fermion_Nc[f]:>4} "
          f"{best_pred_F[f]:>+8.4f} {delta_exp[f]:>+8.4f} {best_pred_F[f]-delta_exp[f]:>+8.4f}")

# Check: do the mode amplitudes have the right sign for QCD?
# If this is a QCD correction, all c_m should be positive (since the
# QCD correction always increases mass magnitude).
print(f"\n  Mode amplitudes: {np.round(best_coeffs_F, 4)}")
print(f"  (For QCD interpretation: should all be positive)")


# ================================================================
# MODEL G: Minimal model - T3 * alpha_s * W_k with NO free scale
# ================================================================
print("\n" + "="*72)
print("MODEL G: Zero free parameters (framework prediction)")
print("="*72)

# The most constrained model: ALL quantities from framework.
# delta_k = T3 * Nc * (alpha_s * C_F / (pi * ln(phi))) * F_k
#
# F_k = sum_{m=2,4,6,8} C^m_kk * d_m / lambda_m
# (propagator with M^2 = 0, i.e., massless exchange)

# With framework alpha_s and Cayley eigenvalues:
F_massless = np.zeros(9)
for k in range(9):
    for m in even_irreps:
        if m == 0:
            continue
        F_massless[k] += CG_self[k, m] * dims[m] / cayley_eigs[m]

pred_G = np.array([2 * fermion_T3[k] * fermion_Nc[k] * prefactor * F_massless[k] for k in range(9)])

rms_G = np.sqrt(np.mean((delta_exp - pred_G)**2))
print(f"\nZero-parameter prediction (r=0 limit, pure framework):")
print(f"  RMS = {rms_G:.4f}")
print(f"\n{'fermion':>8} {'pred':>8} {'exp':>8} {'err':>8}")
print("-"*40)
for k in range(9):
    print(f"{fermion_names[k]:>8} {pred_G[k]:>+8.4f} {delta_exp[k]:>+8.4f} {pred_G[k]-delta_exp[k]:>+8.4f}")

# Also try: delta = T3 * alpha_s (universal)
pred_simple = np.array([fermion_T3[k] * fermion_Nc[k] * ALPHA_S / np.log(PHI) for k in range(9)])
rms_simple = np.sqrt(np.mean((delta_exp - pred_simple)**2))
print(f"\nSimplest: delta = T3*Nc*alpha_s/ln(phi) = T3*Nc*{ALPHA_S/np.log(PHI):.4f}")
print(f"  RMS = {rms_simple:.4f}")
for k in range(9):
    print(f"  {fermion_names[k]:>4}: {pred_simple[k]:>+.4f} vs {delta_exp[k]:>+.4f}")

# Best: delta = T3 * x for quarks (1 param)
# x_opt minimizes sum (T3*x*Nc - delta)^2 for quarks only
quark_idx = [3, 4, 5, 6, 7, 8]
T3_q = np.array([fermion_T3[k] for k in quark_idx])
delta_q = np.array([delta_exp[k] for k in quark_idx])
x_opt = np.sum(T3_q * delta_q) / np.sum(T3_q**2)
pred_T3 = np.zeros(9)
for k in quark_idx:
    pred_T3[k] = fermion_T3[k] * x_opt
rms_T3 = np.sqrt(np.mean((delta_exp - pred_T3)**2))
print(f"\nBest T3-only: delta = T3*{x_opt:.4f} for quarks")
print(f"  RMS = {rms_T3:.4f}")
for k in range(9):
    print(f"  {fermion_names[k]:>4}: {pred_T3[k]:>+.4f} vs {delta_exp[k]:>+.4f}")

# Also: T3 * (a + b*gen) for quarks (2 params)
X_T3g = np.zeros((6, 2))
for i, k in enumerate(quark_idx):
    X_T3g[i, 0] = fermion_T3[k]
    X_T3g[i, 1] = fermion_T3[k] * gen[k]
c_T3g, _, _, _ = np.linalg.lstsq(X_T3g, delta_q, rcond=None)
pred_T3g = np.zeros(9)
for k in quark_idx:
    pred_T3g[k] = fermion_T3[k] * (c_T3g[0] + c_T3g[1] * gen[k])
rms_T3g = np.sqrt(np.mean((delta_exp - pred_T3g)**2))
print(f"\nT3*(a + b*gen): a={c_T3g[0]:.4f}, b={c_T3g[1]:.4f}")
print(f"  RMS = {rms_T3g:.4f}")
for k in range(9):
    print(f"  {fermion_names[k]:>4}: {pred_T3g[k]:>+.4f} vs {delta_exp[k]:>+.4f}")


# ================================================================
# SUMMARY
# ================================================================
print("\n" + "="*72)
print("SUMMARY OF ALL MODELS")
print("="*72)
print(f"""
{'Model':>30} {'Params':>7} {'RMS':>8} {'Notes':>30}
{'-'*80}
{'G: T3*alpha_s/ln(phi)':>30} {'0':>7} {rms_simple:>8.4f} {'Universal T3*alpha_s':>30}
{'G: T3*Nc*pref*F(r=0)':>30} {'0':>7} {rms_G:>8.4f} {'CG-weighted, zero-param':>30}
{'T3*x (quarks only)':>30} {'1':>7} {rms_T3:>8.4f} {'Best single param':>30}
{'A: T3*Nc*F(r)':>30} {'2':>7} {best_rms:>8.4f} {'CG propagator + scale':>30}
{'D: Physical one-loop':>30} {'2':>7} {best_rms_D:>8.4f} {'With alpha_s prefactor':>30}
{'T3*(a+b*gen)':>30} {'2':>7} {rms_T3g:>8.4f} {'T3 + generation':>30}
{'B: T3*Nc*F + Nc*g + c':>30} {'4':>7} {best_rms_B:>8.4f} {'CG + generation + const':>30}
{'C: 4 mode amplitudes':>30} {'4':>7} {best_rms_C:>8.4f} {'Free CG mode amplitudes':>30}
{'E: C + reassign':>30} {'4':>7} {best_rms_E:>8.4f} {'Best reassignment, 4 modes':>30}
{'F: T3*C + reassign':>30} {'4':>7} {best_rms_F:>8.4f} {'T3-signed + reassignment':>30}
{'(Bare formula)':>30} {'0':>7} {'0.2469':>8} {'Reference: delta RMS':>30}
""")
print("Key: lower RMS = better fit. Bare formula has delta RMS = 0.247")
print("     A 'good' model should have RMS < 0.10 with few parameters.")

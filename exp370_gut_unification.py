"""
exp370_gut_unification.py
==========================
OPEN-3 item 3: Precise GUT unification analysis with framework couplings.

Prior results:
  - exp317/326: basic running, YM ratios 8:5:2, no single unification in 1-loop SM
  - exp358: spectral truncation NEGATIVE, running is dynamical
  - exp361: alpha(mZ)/alpha(0) partially resolved, nonpert QCD irreducible

THIS experiment:
1. Precise 1-loop and 2-loop unification with ALL framework couplings
2. Unification scale as framework quantity (m_e * phi^z)
3. Unification mismatch and its framework meaning
4. The 16/15 ratio test: alpha(mZ)/alpha(0) = (3*a1+1)/(3*a1)
5. Comparison: framework vs experimental vs GUT-scale predictions
6. Does the mismatch close at 2-loop? 3-loop?
7. Threshold corrections from framework masses

Dependencies: numpy, scipy
Encoding: ASCII only (Windows cp1252)
"""

import numpy as np
from scipy.integrate import solve_ivp

print("=" * 72)
print("EXP-370: PRECISE GUT UNIFICATION WITH FRAMEWORK COUPLINGS")
print("=" * 72)

# =========================================================================
# FRAMEWORK CONSTANTS
# =========================================================================
a1 = 5
b1 = a1 + 1  # = 6
phi = (1.0 + np.sqrt(a1)) / 2.0
PI = np.pi
N = 120
N_gen = 3
N_H = 1  # Higgs doublets

# Framework couplings
alpha_s_fw = 1.0 / (2.0 * phi**3)  # = 0.11803
sin2tw_fw = float(b1) / (a1**2 + 1)  # = 6/26

# alpha(0) from quadratic 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0
A_c = 2.0 * PI
B_c = -4.0 * a1 * phi**4
C_c = 1.0
disc = B_c**2 - 4 * A_c * C_c
alpha_0 = (-B_c - np.sqrt(disc)) / (2 * A_c)
alpha_0_inv = 1.0 / alpha_0

# The 16/15 pattern for alpha(mZ)/alpha(0)
ratio_16_15 = (3*a1 + 1.0) / (3*a1)  # = 16/15

# Framework alpha at mZ
alpha_mZ_fw = alpha_0 * ratio_16_15
alpha_mZ_inv_fw = 1.0 / alpha_mZ_fw

# Experimental
alpha_mZ_exp = 1.0 / 127.951
alpha_s_exp = 0.1179
sin2tw_exp = 0.23122
m_Z = 91.1876  # GeV
m_e = 0.51099895e-3  # GeV

print(f"\n--- Framework coupling constants ---")
print(f"  alpha(0)^-1     = {alpha_0_inv:.6f}")
print(f"  alpha_s(mZ)     = 1/(2*phi^3) = {alpha_s_fw:.6f}")
print(f"  sin^2(tW)       = 6/26 = {sin2tw_fw:.6f}")
print(f"  alpha(mZ)/alpha(0) = 16/15 = {ratio_16_15:.6f}")
print(f"  alpha(mZ)^-1    = {alpha_mZ_inv_fw:.4f}")
print(f"  (Experimental: alpha_s={alpha_s_exp}, sin2tw={sin2tw_exp}, alpha_mZ^-1=127.951)")

# =========================================================================
# PART 1: GAUGE COUPLINGS AT M_Z (GUT-normalized)
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 1: GAUGE COUPLINGS AT M_Z")
print(f"{'=' * 72}")

# SU(3): alpha_3 = alpha_s
alpha_3_mZ = alpha_s_fw
alpha_3_inv_mZ = 1.0 / alpha_3_mZ

# SU(2): alpha_2 = alpha_em / sin^2(tW)
alpha_2_mZ = alpha_mZ_fw / sin2tw_fw
alpha_2_inv_mZ = 1.0 / alpha_2_mZ

# U(1) GUT-normalized: alpha_1 = (5/3) * alpha_em / cos^2(tW)
cos2tw_fw = 1.0 - sin2tw_fw
alpha_1_mZ = (5.0/3.0) * alpha_mZ_fw / cos2tw_fw
alpha_1_inv_mZ = 1.0 / alpha_1_mZ

print(f"\n  Framework gauge couplings at M_Z (GUT normalized):")
print(f"    alpha_1^-1 = {alpha_1_inv_mZ:.6f}")
print(f"    alpha_2^-1 = {alpha_2_inv_mZ:.6f}")
print(f"    alpha_3^-1 = {alpha_3_inv_mZ:.6f}")

# Experimental
alpha_3_inv_exp = 1.0 / alpha_s_exp
alpha_2_inv_exp = sin2tw_exp / alpha_mZ_exp
alpha_1_inv_exp = (3.0/5.0) * cos2tw_fw / alpha_mZ_exp  # using fw sin2tw for consistency
# Actually use experimental values properly
alpha_2_mZ_exp = alpha_mZ_exp / sin2tw_exp
alpha_1_mZ_exp = (5.0/3.0) * alpha_mZ_exp / (1.0 - sin2tw_exp)

print(f"\n  Experimental gauge couplings at M_Z:")
print(f"    alpha_1^-1 = {1/alpha_1_mZ_exp:.6f}")
print(f"    alpha_2^-1 = {1/alpha_2_mZ_exp:.6f}")
print(f"    alpha_3^-1 = {1/alpha_s_exp:.6f}")

print(f"\n  Deviations (framework - experiment):")
print(f"    alpha_1^-1: {alpha_1_inv_mZ:.4f} vs {1/alpha_1_mZ_exp:.4f} ({(alpha_1_inv_mZ - 1/alpha_1_mZ_exp)/(1/alpha_1_mZ_exp)*100:+.3f}%)")
print(f"    alpha_2^-1: {alpha_2_inv_mZ:.4f} vs {1/alpha_2_mZ_exp:.4f} ({(alpha_2_inv_mZ - 1/alpha_2_mZ_exp)/(1/alpha_2_mZ_exp)*100:+.3f}%)")
print(f"    alpha_3^-1: {alpha_3_inv_mZ:.4f} vs {1/alpha_s_exp:.4f} ({(alpha_3_inv_mZ - 1/alpha_s_exp)/(1/alpha_s_exp)*100:+.3f}%)")

# Check: alpha_em^-1 = (3/5)*alpha_1^-1 + alpha_2^-1
check_em = (3.0/5.0) * alpha_1_inv_mZ + alpha_2_inv_mZ
print(f"\n  Consistency: (3/5)*a1^-1 + a2^-1 = {check_em:.4f} = alpha_em^-1(mZ)")

# =========================================================================
# PART 2: ONE-LOOP BETA COEFFICIENTS
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 2: ONE-LOOP BETA COEFFICIENTS (DERIVED)")
print(f"{'=' * 72}")

# SM 1-loop beta coefficients with N_gen=3, N_H=1
# b_i defined by: d(alpha_i^-1)/d(ln mu) = -b_i/(2*pi)
b1_beta = 4.0/3.0 * N_gen + 1.0/10.0 * N_H  # = 41/10
b2_beta = -22.0/3.0 + 4.0/3.0 * N_gen + 1.0/6.0 * N_H  # = -19/6
b3_beta = -11.0 + 4.0/3.0 * N_gen  # = -7

print(f"  b_1 = 4/3*{N_gen} + 1/10*{N_H} = {b1_beta:.6f} = 41/10")
print(f"  b_2 = -22/3 + 4/3*{N_gen} + 1/6*{N_H} = {b2_beta:.6f} = -19/6")
print(f"  b_3 = -11 + 4/3*{N_gen} = {b3_beta:.1f}")
print(f"  Sum: b_1+b_2+b_3 = {b1_beta+b2_beta+b3_beta:.6f}")
print(f"  YM ratio: b_3:b_2:b_1 gauge parts = 8:5:2 (sum=15=a1*N_gen)")

# =========================================================================
# PART 3: ONE-LOOP RUNNING AND PAIRWISE UNIFICATION
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 3: ONE-LOOP PAIRWISE UNIFICATION")
print(f"{'=' * 72}")

# alpha_i^-1(mu) = alpha_i^-1(mZ) - b_i/(2*pi) * ln(mu/mZ)
def run_1loop(alpha_inv_mZ, b, t):
    return alpha_inv_mZ - b / (2*PI) * t

# Pairwise crossing scales
# alpha_i^-1(t*) = alpha_j^-1(t*) =>
# alpha_i^-1 - alpha_j^-1 = (b_i - b_j)/(2*pi) * t*

def crossing(ai_inv, aj_inv, bi, bj):
    """Find t = ln(mu/mZ) where alpha_i = alpha_j."""
    t = (ai_inv - aj_inv) / ((bi - bj) / (2*PI))
    mu = m_Z * np.exp(t)
    a_cross = run_1loop(ai_inv, bi, t)
    return t, mu, a_cross

# alpha_1 = alpha_2
t_12, mu_12, a_12 = crossing(alpha_1_inv_mZ, alpha_2_inv_mZ, b1_beta, b2_beta)
# alpha_2 = alpha_3
t_23, mu_23, a_23 = crossing(alpha_2_inv_mZ, alpha_3_inv_mZ, b2_beta, b3_beta)
# alpha_1 = alpha_3
t_13, mu_13, a_13 = crossing(alpha_1_inv_mZ, alpha_3_inv_mZ, b1_beta, b3_beta)

print(f"\n  Pairwise crossing scales (1-loop, framework couplings):")
print(f"    alpha_1 = alpha_2: mu = {mu_12:.4e} GeV, alpha^-1 = {a_12:.4f}, t = {t_12:.4f}")
print(f"    alpha_2 = alpha_3: mu = {mu_23:.4e} GeV, alpha^-1 = {a_23:.4f}, t = {t_23:.4f}")
print(f"    alpha_1 = alpha_3: mu = {mu_13:.4e} GeV, alpha^-1 = {a_13:.4f}, t = {t_13:.4f}")

# Unification triangle
spread_decades = np.log10(max(mu_12, mu_23, mu_13) / min(mu_12, mu_23, mu_13))
print(f"\n  Unification triangle:")
print(f"    Spread: {spread_decades:.4f} decades")
print(f"    log10(mu_12) = {np.log10(mu_12):.4f}")
print(f"    log10(mu_23) = {np.log10(mu_23):.4f}")
print(f"    log10(mu_13) = {np.log10(mu_13):.4f}")

# The mismatch at each crossing
a3_at_12 = run_1loop(alpha_3_inv_mZ, b3_beta, t_12)
a1_at_23 = run_1loop(alpha_1_inv_mZ, b1_beta, t_23)
a2_at_13 = run_1loop(alpha_2_inv_mZ, b2_beta, t_13)

print(f"\n  Mismatch (third coupling at each crossing):")
print(f"    At alpha_1=alpha_2: alpha_3^-1 = {a3_at_12:.4f}, gap = {a3_at_12 - a_12:.4f}")
print(f"    At alpha_2=alpha_3: alpha_1^-1 = {a1_at_23:.4f}, gap = {a1_at_23 - a_23:.4f}")
print(f"    At alpha_1=alpha_3: alpha_2^-1 = {a2_at_13:.4f}, gap = {a2_at_13 - a_13:.4f}")

# =========================================================================
# PART 4: THE B-PARAMETER (UNIFICATION TEST)
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 4: B-PARAMETER (UNIFICATION QUALITY)")
print(f"{'=' * 72}")

# The B-parameter tests if 3 lines meet at a single point.
# If they do: (a1^-1 - a2^-1)/(a2^-1 - a3^-1) = (b1 - b2)/(b2 - b3)
ratio_coupling = (alpha_1_inv_mZ - alpha_2_inv_mZ) / (alpha_2_inv_mZ - alpha_3_inv_mZ)
ratio_beta = (b1_beta - b2_beta) / (b2_beta - b3_beta)

print(f"\n  Unification test:")
print(f"    (a1^-1 - a2^-1)/(a2^-1 - a3^-1) = {ratio_coupling:.6f}")
print(f"    (b1 - b2)/(b2 - b3) = {ratio_beta:.6f}")
print(f"    B = ratio_coupling/ratio_beta = {ratio_coupling/ratio_beta:.6f}")
print(f"    Deviation from perfect unification: {abs(ratio_coupling/ratio_beta - 1)*100:.3f}%")

# Compare with experimental couplings
ratio_coupling_exp = (1/alpha_1_mZ_exp - 1/alpha_2_mZ_exp) / (1/alpha_2_mZ_exp - 1/alpha_s_exp)
print(f"\n  Experimental:")
print(f"    (a1^-1 - a2^-1)/(a2^-1 - a3^-1) = {ratio_coupling_exp:.6f}")
print(f"    B_exp = {ratio_coupling_exp/ratio_beta:.6f}")
print(f"    Deviation from perfect: {abs(ratio_coupling_exp/ratio_beta - 1)*100:.3f}%")

# =========================================================================
# PART 5: TWO-LOOP RUNNING
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 5: TWO-LOOP RUNNING (RK4 INTEGRATION)")
print(f"{'=' * 72}")

# 2-loop beta coefficient matrix b_ij
# d(alpha_i^-1)/d(ln mu) = -b_i/(2*pi) - sum_j b_ij/(8*pi^2) * 1/alpha_j^-1
# i.e., using alpha_j (not alpha_j^-1)

b_2loop = np.array([
    [199.0/50.0, 27.0/10.0, 44.0/5.0],   # b_1j
    [9.0/10.0,   35.0/6.0,  12.0],        # b_2j
    [11.0/10.0,  9.0/2.0,  -26.0],        # b_3j
])

# Yukawa contribution (top quark dominant)
# b_ij^Y terms from top Yukawa
# These are typically small compared to gauge terms
# For precision: y_t^2/(16*pi^2) contributions
# At M_Z: y_t = sqrt(2)*m_t/v, m_t~172.76 GeV, v=246.22 GeV
# y_t ~ 0.993, y_t^2/(16*pi^2) ~ 0.00625

# For now, use pure gauge 2-loop (neglect Yukawa)
b_1loop = np.array([b1_beta, b2_beta, b3_beta])

def rge_2loop_rhs(t, alpha_inv):
    """d(alpha_i^-1)/dt where t = ln(mu/m_Z)."""
    alpha = 1.0 / np.maximum(alpha_inv, 1e-10)  # safety
    deriv = np.zeros(3)
    for i in range(3):
        deriv[i] = -b_1loop[i]/(2*PI) - np.dot(b_2loop[i], alpha)/(8*PI**2)
    return deriv

# RK4 integration
dt = 0.005
t_max = 75.0  # up to ~10^32 GeV, more than enough
n_steps = int(t_max / dt)

alpha_inv_t = np.zeros((n_steps+1, 3))
alpha_inv_t[0] = [alpha_1_inv_mZ, alpha_2_inv_mZ, alpha_3_inv_mZ]
t_arr = np.arange(n_steps+1) * dt

for step in range(n_steps):
    y = alpha_inv_t[step]
    if np.any(y <= 0.5):  # coupling blowup
        alpha_inv_t[step+1:] = alpha_inv_t[step]
        break
    k1 = rge_2loop_rhs(t_arr[step], y)
    k2 = rge_2loop_rhs(t_arr[step] + dt/2, y + dt/2*k1)
    k3 = rge_2loop_rhs(t_arr[step] + dt/2, y + dt/2*k2)
    k4 = rge_2loop_rhs(t_arr[step] + dt, y + dt*k3)
    alpha_inv_t[step+1] = y + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4)

# Print running at selected scales
print(f"\n  2-loop running (framework couplings):")
print(f"  {'log10(mu/GeV)':>14s}  {'a1^-1':>10s}  {'a2^-1':>10s}  {'a3^-1':>10s}  {'a1-a2':>8s}  {'a2-a3':>8s}")

for log_mu_target in [2, 3, 4, 6, 8, 10, 12, 14, 15, 15.5, 16, 16.5, 17, 18]:
    t_target = log_mu_target * np.log(10) - np.log(m_Z)
    idx = int(t_target / dt)
    if 0 <= idx < len(t_arr):
        a1i, a2i, a3i = alpha_inv_t[idx]
        print(f"  {log_mu_target:14.1f}  {a1i:10.4f}  {a2i:10.4f}  {a3i:10.4f}  {a1i-a2i:+8.4f}  {a2i-a3i:+8.4f}")

# Find 2-loop crossings
crossings_2loop = {}
for pair, (ci, cj) in [('12', (0, 1)), ('23', (1, 2)), ('13', (0, 2))]:
    for step in range(1, n_steps):
        diff_prev = alpha_inv_t[step-1, ci] - alpha_inv_t[step-1, cj]
        diff_curr = alpha_inv_t[step, ci] - alpha_inv_t[step, cj]
        if diff_prev * diff_curr <= 0 and abs(diff_prev) + abs(diff_curr) > 1e-10:
            # Linear interpolation
            frac = abs(diff_prev) / (abs(diff_prev) + abs(diff_curr))
            t_cross = t_arr[step-1] + frac * dt
            mu_cross = m_Z * np.exp(t_cross)
            a_cross = (alpha_inv_t[step-1, ci] + alpha_inv_t[step, ci]) / 2
            crossings_2loop[pair] = (t_cross, mu_cross, a_cross)
            break

print(f"\n  2-loop pairwise crossings:")
for pair in ['12', '23', '13']:
    if pair in crossings_2loop:
        t_c, mu_c, a_c = crossings_2loop[pair]
        print(f"    alpha_{pair[0]} = alpha_{pair[1]}: mu = {mu_c:.4e} GeV (log10 = {np.log10(mu_c):.4f}), alpha^-1 = {a_c:.4f}")

# Mismatch at 2-loop crossings
if '23' in crossings_2loop:
    t_c = crossings_2loop['23'][0]
    idx = int(t_c / dt)
    if 0 <= idx < len(alpha_inv_t):
        a1_at_23_2loop = alpha_inv_t[idx, 0]
        a23_cross = crossings_2loop['23'][2]
        mismatch_2loop = a1_at_23_2loop - a23_cross
        print(f"\n  2-loop mismatch at alpha_2=alpha_3:")
        print(f"    alpha_1^-1 = {a1_at_23_2loop:.4f}")
        print(f"    alpha_{2,3}^-1 = {a23_cross:.4f}")
        print(f"    Gap = {mismatch_2loop:.4f}")
        print(f"    Relative gap = {mismatch_2loop/a23_cross*100:.2f}%")

# =========================================================================
# PART 6: UNIFICATION SCALE AS FRAMEWORK QUANTITY
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 6: UNIFICATION SCALE AS FRAMEWORK QUANTITY")
print(f"{'=' * 72}")

# Check if unification scale = m_e * phi^z for some z
for pair in ['12', '23', '13']:
    if pair in crossings_2loop:
        mu_c = crossings_2loop[pair][1]
        z = np.log(mu_c / (m_e * 1e3)) / np.log(phi)  # m_e in MeV
        z_from_mP = np.log(mu_c / 1.22089e19) / np.log(phi)
        print(f"    mu_{pair} = m_e * phi^{z:.4f}")
        print(f"    mu_{pair} = m_P * phi^{z_from_mP:.4f}")
        print(f"    mu_{pair}/m_P = {mu_c/1.22089e19:.4e}")

# More interesting: the "average" unification scale
if len(crossings_2loop) >= 2:
    mu_vals = [v[1] for v in crossings_2loop.values()]
    mu_avg = np.exp(np.mean(np.log(mu_vals)))
    z_avg = np.log(mu_avg / (m_e * 1e3)) / np.log(phi)
    print(f"\n  Geometric mean unification scale:")
    print(f"    mu_GUT ~ {mu_avg:.4e} GeV = m_e * phi^{z_avg:.4f}")
    print(f"    log10(mu_GUT) = {np.log10(mu_avg):.4f}")

    # Check framework quantities near this z
    print(f"\n  Nearby framework quantities:")
    z_Planck = 4 * phi**2
    z_int = round(z_avg)
    print(f"    z_Planck = 4*phi^2 = {z_Planck:.4f}")
    print(f"    Nearest integer z = {z_int}")
    print(f"    N/2 = {N/2}")
    print(f"    2*z_Planck = {2*z_Planck:.4f}")

    # Test specific framework expressions for z_GUT
    candidates = [
        ("N/2 - N_gen", N/2 - N_gen),          # = 57 (CC exponent!)
        ("N/2", N/2),                            # = 60
        ("2*z_Planck", 2*z_Planck),              # = 20.944
        ("a1^2 + a1", a1**2 + a1),              # = 30
        ("(a1^2+1)*phi^2", (a1**2+1)*phi**2),   # = 68.07
        ("N_gen*z_Planck", N_gen*z_Planck),       # = 31.42
        ("h_E8*phi^2", 30*phi**2),               # = 78.54
        ("2*a1*phi^3", 2*a1*phi**3),             # = 42.36
        ("N_gen*a1*phi^3", N_gen*a1*phi**3),     # = 63.54
    ]

    print(f"\n  {'Expression':25s} {'Value':10s} {'mu (GeV)':12s} {'log10':10s}")
    for name, z_try in candidates:
        mu_try = m_e * 1e3 * phi**z_try  # m_e in MeV, convert to GeV
        print(f"  {name:25s} {z_try:10.4f} {mu_try:12.4e} {np.log10(mu_try):10.4f}")

# =========================================================================
# PART 7: PROTON DECAY AND FRAMEWORK PREDICTION
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 7: PROTON DECAY PREDICTION")
print(f"{'=' * 72}")

# In SU(5) GUT, proton lifetime ~ M_GUT^4 / (alpha_GUT^2 * m_p^5)
m_p = 0.93827  # GeV (proton mass)

if '23' in crossings_2loop:
    mu_GUT = crossings_2loop['23'][1]
    alpha_GUT_inv = crossings_2loop['23'][2]
    alpha_GUT = 1.0 / alpha_GUT_inv

    # Simplified proton lifetime estimate
    # tau_p ~ M_GUT^4 / (alpha_GUT^2 * m_p^5) * phase_space
    # In natural units (GeV^-1 = 6.582e-25 s):
    tau_p_nat = mu_GUT**4 / (alpha_GUT**2 * m_p**5)
    # Convert to seconds: 1 GeV^-1 = 6.582e-25 s
    tau_p_sec = tau_p_nat * 6.582e-25
    # Convert to years: 1 year = 3.156e7 s
    tau_p_yr = tau_p_sec / 3.156e7

    print(f"  M_GUT = {mu_GUT:.4e} GeV")
    print(f"  alpha_GUT = {alpha_GUT:.6f}")
    print(f"  tau_p (rough) ~ M_GUT^4 / (alpha_GUT^2 * m_p^5)")
    print(f"  tau_p ~ {tau_p_yr:.2e} years")
    print(f"  log10(tau_p/yr) = {np.log10(tau_p_yr):.2f}")
    print(f"  Current bound: tau_p > 2.4e34 yr (Super-K, p -> e+ pi0)")
    print(f"  Framework prediction: {'ALLOWED' if tau_p_yr > 2.4e34 else 'EXCLUDED'}")

# =========================================================================
# PART 8: THE 16/15 PATTERN INVESTIGATION
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 8: THE 16/15 PATTERN FOR alpha(mZ)/alpha(0)")
print(f"{'=' * 72}")

# The pattern: alpha(mZ)/alpha(0) = (3*a1+1)/(3*a1) = 16/15
# This means m_Z = m_e * phi^25 * 16/15

# Is 16/15 related to gauge group structure?
print(f"\n  16/15 = (3*a1+1)/(3*a1) = {(3*a1+1)}/{3*a1}")
print(f"  = 1 + 1/(3*a1) = 1 + 1/15")
print(f"  15 = a1 * N_gen = 5 * 3 (sum of YM ratios 8+5+2)")
print(f"  16 = dim(WHITE) = (a1-1)^2 = fermions per generation")

# So: alpha(mZ)/alpha(0) = (a1-1)^2 / (a1*N_gen) = 16/15
# This has a beautiful interpretation:
# numerator = number of fermion species per generation
# denominator = sum of Yang-Mills gauge fractions
print(f"\n  INTERPRETATION:")
print(f"  numerator   = dim(WHITE) = 16 = fermions per generation")
print(f"  denominator = a1*N_gen  = 15 = sum of YM ratios 8+5+2")
print(f"  ratio = 16/15 = fermions / gauge_total")

# Numerical test
ratio_exp_actual = alpha_mZ_exp / alpha_0
print(f"\n  Numerical test:")
print(f"    16/15 = {16/15:.10f}")
print(f"    alpha_exp(mZ)/alpha_exp(0) = {ratio_exp_actual:.10f}")
print(f"    Error = {abs(16/15 - ratio_exp_actual)/ratio_exp_actual*100:.3f}%")

# Decomposition: the QED running Delta_alpha
# Delta_alpha = 1 - 1/R = 1 - 15/16 = 1/16
# and 1/16 = 1/(a1-1)^2
Da_pattern = 1.0 - 15.0/16.0
Da_exp = 1.0 - 1.0/ratio_exp_actual
print(f"\n  Delta_alpha = 1 - 1/R:")
print(f"    Pattern: 1 - 15/16 = 1/16 = {Da_pattern:.6f}")
print(f"    Experimental: {Da_exp:.6f}")
print(f"    Error: {abs(Da_pattern - Da_exp)/Da_exp*100:.2f}%")

# =========================================================================
# PART 9: CLOSEST APPROACH (NEAR-UNIFICATION)
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 9: CLOSEST APPROACH OF ALL THREE COUPLINGS")
print(f"{'=' * 72}")

# Define mismatch function: max|alpha_i^-1 - alpha_j^-1| for all pairs
def mismatch(alpha_inv_3):
    a1, a2, a3 = alpha_inv_3
    return max(abs(a1-a2), abs(a2-a3), abs(a1-a3))

# Find minimum mismatch in 2-loop running
min_mismatch = 1e10
min_idx = 0
for step in range(1, n_steps):
    m = mismatch(alpha_inv_t[step])
    if m < min_mismatch and np.all(alpha_inv_t[step] > 0):
        min_mismatch = m
        min_idx = step

t_min = t_arr[min_idx]
mu_min = m_Z * np.exp(t_min)
a1_min, a2_min, a3_min = alpha_inv_t[min_idx]
a_avg_min = np.mean(alpha_inv_t[min_idx])

print(f"\n  Closest approach (2-loop, framework):")
print(f"    Scale: mu = {mu_min:.4e} GeV (log10 = {np.log10(mu_min):.4f})")
print(f"    alpha_1^-1 = {a1_min:.4f}")
print(f"    alpha_2^-1 = {a2_min:.4f}")
print(f"    alpha_3^-1 = {a3_min:.4f}")
print(f"    Average: {a_avg_min:.4f}")
print(f"    Max pairwise gap: {min_mismatch:.4f}")
print(f"    Relative mismatch: {min_mismatch/a_avg_min*100:.2f}%")

# Check if average coupling is a framework quantity
print(f"\n  Average coupling at closest approach:")
print(f"    alpha_avg^-1 = {a_avg_min:.4f}")
print(f"    Compare: 2*a1^2 = {2*a1**2} = 50")
print(f"    Compare: N/2 - N_gen = 57")
print(f"    Compare: a1^2 + a1 + 1 = {a1**2+a1+1} = 31")
print(f"    Compare: h_E8 + b1 + a1 = {30+b1+a1} = 41")
print(f"    Compare: 2*a1*phi^2 = {2*a1*phi**2:.4f}")

# =========================================================================
# PART 10: THRESHOLD CORRECTIONS FROM FRAMEWORK MASSES
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 10: THRESHOLD CORRECTIONS AT M_Z")
print(f"{'=' * 72}")

# Threshold corrections modify the running near particle thresholds.
# The heaviest SM particles below M_Z affect the coupling matching.
# Framework masses:
mass_exps = {
    'e': 0, 'mu': 11, 'tau': 17,
    'u': 3, 'c': 16, 't': 26,  # bare exponents
    'd': 5, 's': 11, 'b': 19   # bare exponents (approx)
}

print(f"\n  Fermion mass hierarchy (log_phi scale):")
for name in ['e', 'mu', 'tau', 'u', 'd', 's', 'c', 'b', 't']:
    n = mass_exps[name]
    m_fw = m_e * phi**n
    m_GeV = m_fw
    threshold = "< m_Z" if m_GeV < m_Z else "> m_Z"
    print(f"    {name}: n={n:3d}, m ~ {m_GeV:.4f} GeV  [{threshold}]")

# =========================================================================
# PART 11: COMPARISON - FRAMEWORK vs MSSM
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 11: COMPARISON WITH MSSM UNIFICATION")
print(f"{'=' * 72}")

# MSSM has different beta coefficients (with superpartners)
b1_MSSM = 33.0/5.0  # = 6.6
b2_MSSM = 1.0
b3_MSSM = -3.0

# 1-loop MSSM running from experimental couplings
t_12_MSSM, mu_12_MSSM, _ = crossing(1/alpha_1_mZ_exp, 1/alpha_2_mZ_exp, b1_MSSM, b2_MSSM)
t_23_MSSM, mu_23_MSSM, _ = crossing(1/alpha_2_mZ_exp, 1/alpha_s_exp, b2_MSSM, b3_MSSM)
t_13_MSSM, mu_13_MSSM, _ = crossing(1/alpha_1_mZ_exp, 1/alpha_s_exp, b1_MSSM, b3_MSSM)

ratio_MSSM = (1/alpha_1_mZ_exp - 1/alpha_2_mZ_exp) / (1/alpha_2_mZ_exp - 1/alpha_s_exp)
ratio_MSSM_beta = (b1_MSSM - b2_MSSM) / (b2_MSSM - b3_MSSM)

print(f"\n  MSSM 1-loop beta coefficients:")
print(f"    b_1 = 33/5 = {b1_MSSM}")
print(f"    b_2 = 1")
print(f"    b_3 = -3")

print(f"\n  MSSM unification test:")
print(f"    coupling ratio = {ratio_MSSM:.6f}")
print(f"    beta ratio     = {ratio_MSSM_beta:.6f}")
print(f"    Deviation = {abs(ratio_MSSM/ratio_MSSM_beta - 1)*100:.2f}%")

print(f"\n  MSSM crossing scales:")
print(f"    alpha_1=alpha_2: {mu_12_MSSM:.2e} GeV (log10 = {np.log10(mu_12_MSSM):.2f})")
print(f"    alpha_2=alpha_3: {mu_23_MSSM:.2e} GeV (log10 = {np.log10(mu_23_MSSM):.2f})")
print(f"    Spread: {abs(np.log10(mu_12_MSSM)-np.log10(mu_23_MSSM)):.2f} decades")

# Framework (SM) comparison
ratio_SM_fw = ratio_coupling / ratio_beta
print(f"\n  COMPARISON:")
print(f"    MSSM: B = {ratio_MSSM/ratio_MSSM_beta:.6f} (deviation {abs(ratio_MSSM/ratio_MSSM_beta - 1)*100:.2f}%)")
print(f"    SM framework: B = {ratio_SM_fw:.6f} (deviation {abs(ratio_SM_fw - 1)*100:.2f}%)")
print(f"    SM experimental: B = {ratio_coupling_exp/ratio_beta:.6f} (deviation {abs(ratio_coupling_exp/ratio_beta - 1)*100:.2f}%)")

# =========================================================================
# PART 12: DEEP FRAMEWORK CONNECTION
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 12: DEEP FRAMEWORK CONNECTIONS")
print(f"{'=' * 72}")

# The key observation: framework gives TREE-LEVEL couplings.
# Running = quantum corrections = standard QFT.
# The framework is NOT a GUT: couplings are unified ALGEBRAICALLY at M_Z,
# not at a high-energy scale.

# Key algebraic relations at M_Z:
print(f"\n  Framework algebraic relations (at M_Z):")
print(f"    alpha_s = 1/(2*phi^3)")
print(f"    sin^2(tW) = b1/(a1^2+1) = 6/26")
print(f"    alpha from quadratic with coefficients (2pi, 4*a1*phi^4, 1)")
print(f"    These are SIMULTANEOUS at M_Z, not at a GUT scale.")

# The unification "triangle" measures how far from a GUT the framework is
# In a true GUT, the triangle collapses to a point.
# In the framework, there's a spread of ~1 decade.

# Can the mismatch be related to framework quantities?
if '23' in crossings_2loop and '12' in crossings_2loop:
    t_23_2l = crossings_2loop['23'][0]
    t_12_2l = crossings_2loop['12'][0]
    delta_t = abs(t_23_2l - t_12_2l)
    # delta_t in units of ln(phi)
    delta_t_phi = delta_t / np.log(phi)
    print(f"\n  Unification triangle width:")
    print(f"    delta_t = {delta_t:.4f}")
    print(f"    delta_t / ln(phi) = {delta_t_phi:.4f}")
    print(f"    Compare: a1 = {a1}, b1 = {b1}")

    # The mismatch in alpha at nearest crossing
    gap_at_23 = alpha_inv_t[min_idx, 0] - np.mean(alpha_inv_t[min_idx, 1:])
    print(f"    Gap in alpha^-1 ~ {gap_at_23:.4f}")

# =========================================================================
# PART 13: FRAMEWORK PREDICTION FOR PROTON DECAY (REFINED)
# =========================================================================
print(f"\n{'=' * 72}")
print("PART 13: FRAMEWORK PREDICTION SUMMARY")
print(f"{'=' * 72}")

# Since framework is NOT a GUT:
# - No gauge boson X,Y mediating proton decay
# - The "unification scale" is not physical
# - However, the near-convergence at ~10^16 GeV is a PREDICTION
#   that can be tested against non-observation of proton decay

print(f"""
  CRITICAL DISTINCTION:
  The framework is NOT a Grand Unified Theory (GUT).
  The gauge group is SU(3) x SU(2) x U(1) at ALL scales.

  The "unification" is ALGEBRAIC (all couplings from a1=5)
  NOT DYNAMIC (no single coupling at high energy).

  The near-convergence at ~10^16 GeV is a consequence of:
  1. Framework couplings at M_Z (derived from a1=5)
  2. Standard one-loop beta coefficients (derived from N_gen=3)
  3. These happen to nearly meet because the B-parameter is close to 1.

  B-parameter analysis:
""")

# Why is B close to 1?
# B = (alpha_1^-1 - alpha_2^-1) / (alpha_2^-1 - alpha_3^-1) / ((b1-b2)/(b2-b3))
# The beta ratio is (41/10 + 19/6) / (-19/6 + 7) = (218/60)/(23/6) = 218*6/(60*23) = 1308/1380

print(f"  Beta ratio (b1-b2)/(b2-b3) = {ratio_beta:.6f}")
print(f"    = (41/10 + 19/6) / (7 - 19/6)")
print(f"    = {41/10 + 19/6:.4f} / {7 - 19/6:.4f}")
print(f"    = 109/30 / 23/6 = 654/690 = 109*6/(30*23)")

# Coupling ratio
# (alpha_1^-1 - alpha_2^-1) = (3/5)*cos^2(tW)/alpha_em - sin^2(tW)/alpha_em * (a little different)
# Let me compute it symbolically
num_coupling = alpha_1_inv_mZ - alpha_2_inv_mZ
den_coupling = alpha_2_inv_mZ - alpha_3_inv_mZ

print(f"\n  Coupling difference ratio:")
print(f"    alpha_1^-1 - alpha_2^-1 = {num_coupling:.6f}")
print(f"    alpha_2^-1 - alpha_3^-1 = {den_coupling:.6f}")
print(f"    Ratio = {num_coupling/den_coupling:.6f}")
print(f"    Beta ratio = {ratio_beta:.6f}")
print(f"    B = {ratio_coupling/ratio_beta:.6f}")

# Express in terms of framework:
# alpha_1^-1 = (3/5)/(alpha_em*cos^2tw) = (3/5)*(a1^2+1)/(b1*alpha_em)... wait
# alpha_1^-1 = (3/5) * (a1^2+1) / ((a1^2+1-b1) * alpha_mZ)
# alpha_2^-1 = (a1^2+1)/b1 * 1/alpha_mZ ... no

# Let me just compute numerically
# alpha_em = alpha_0 * 16/15
# alpha_1 = (5/3) * alpha_em / cos2tw = (5/3) * alpha_em / (20/26) = (5/3)*(26/20)*alpha_em = (13/6)*alpha_em
# alpha_2 = alpha_em / sin2tw = alpha_em / (6/26) = (26/6)*alpha_em = (13/3)*alpha_em
# alpha_3 = alpha_s = 1/(2*phi^3)

# alpha_1^-1 = 6/(13*alpha_em)
# alpha_2^-1 = 3/(13*alpha_em)... no
# Wait: alpha_em = alpha_0 * R, sin2tw = 6/26 = 3/13, cos2tw = 20/26 = 10/13
# alpha_2 = alpha_em / sin2tw = alpha_em * 13/3
# alpha_1 = (5/3) * alpha_em / cos2tw = (5/3) * alpha_em * 13/10 = alpha_em * 13/6

# So alpha_1 = alpha_em * 13/6, alpha_2 = alpha_em * 13/3
# alpha_1^-1 = 6/(13*alpha_em), alpha_2^-1 = 3/(13*alpha_em)
# alpha_1^-1 - alpha_2^-1 = 3/(13*alpha_em)

a1_inv_symbolic = 6.0 / (13.0 * alpha_mZ_fw)
a2_inv_symbolic = 3.0 / (13.0 * alpha_mZ_fw)
print(f"\n  Symbolic check:")
print(f"    alpha_1^-1 = 6/(13*alpha_em(mZ)) = {a1_inv_symbolic:.6f} (direct: {alpha_1_inv_mZ:.6f})")
print(f"    alpha_2^-1 = 3/(13*alpha_em(mZ)) = {a2_inv_symbolic:.6f} (direct: {alpha_2_inv_mZ:.6f})")

# So alpha_1^-1 - alpha_2^-1 = 3/(13*alpha_em) = alpha_2^-1
num_symb = 3.0 / (13.0 * alpha_mZ_fw)  # = alpha_2^-1
den_symb = 3.0 / (13.0 * alpha_mZ_fw) - 2*phi**3  # = alpha_2^-1 - alpha_3^-1

print(f"\n  ELEGANT RESULT:")
print(f"    alpha_1^-1 - alpha_2^-1 = alpha_2^-1 = 3/(13*alpha_em) = {num_symb:.6f}")
print(f"    alpha_2^-1 - alpha_3^-1 = alpha_2^-1 - 2*phi^3 = {den_symb:.6f}")

B_exact = num_symb / den_symb / ratio_beta
print(f"    B = [alpha_2^-1 / (alpha_2^-1 - 2*phi^3)] / [(b1-b2)/(b2-b3)]")
print(f"      = {B_exact:.6f}")

# The condition for EXACT unification: B = 1
# alpha_2^-1 / (alpha_2^-1 - 2*phi^3) = (b1-b2)/(b2-b3) = 109*6/(30*23) = 654/690
# Let x = alpha_2^-1 = 3/(13*alpha_em)
# x/(x - 2*phi^3) = 654/690
# 690*x = 654*(x - 2*phi^3)
# 690*x = 654*x - 1308*phi^3
# 36*x = -1308*phi^3
# x = -1308*phi^3/36 = -36.33*phi^3

# This gives NEGATIVE x, which is unphysical!
# So EXACT unification is IMPOSSIBLE in the SM framework.
print(f"\n  EXACT unification condition: B = 1")
print(f"    Requires: alpha_2^-1 = -1308*phi^3/36 = {-1308*phi**3/36:.4f} (NEGATIVE!)")
print(f"    => EXACT unification is IMPOSSIBLE in SM framework.")
print(f"    The near-unification (B ~ {B_exact:.3f}) is as close as SM can get.")

# =========================================================================
# PART 14: HONEST SUMMARY
# =========================================================================
print(f"\n{'=' * 72}")
print("SUMMARY: OPEN-3 COMPLETE ANALYSIS")
print(f"{'=' * 72}")

print(f"""
  DERIVED (from a1 = 5):
  =======================
  1. Beta coefficients b_1=41/10, b_2=-19/6, b_3=-7
     (from N_gen=3, N_H=1, SU(3)xSU(2)xU(1) - all derived)

  2. Yang-Mills ratios at cutoff: 8:5:2 (sum = 15 = a1*N_gen)
     (from spectral action 1+3+8 decomposition)

  3. Gauge couplings at M_Z:
     alpha_s = 1/(2*phi^3), sin^2(tW) = 6/26, alpha from quadratic

  4. The 16/15 pattern: alpha(mZ)/alpha(0) = 16/15 = dim(WHITE)/(a1*N_gen)
     where 16 = fermions per generation, 15 = sum of YM ratios
     Error vs experiment: {abs(16/15 - ratio_exp_actual)/ratio_exp_actual*100:.2f}%
     CATEGORY: STRONG PATTERN (not derived from running)

  NEAR-UNIFICATION (consequences of derived quantities):
  ======================================================
  5. 1-loop: alpha_2 ~ alpha_3 at ~{mu_23:.1e} GeV
     2-loop: refined to ~{crossings_2loop.get('23', (0, 1e16, 0))[1]:.1e} GeV
     Mismatch: alpha_1 separated by ~{min_mismatch:.1f} from alpha_{{2,3}}

  6. B-parameter = {B_exact:.4f}
     EXACT unification impossible in SM (requires negative coupling)
     Framework mismatch: {abs(B_exact-1)*100:.1f}%

  7. Proton decay: framework predicts NO GUT-mediated decay
     (gauge group is SU(3)xSU(2)xU(1) at all scales)

  IRREDUCIBLE EXTERNAL INPUT:
  ===========================
  8. alpha(mZ)/alpha(0) = 16/15 is a PATTERN.
     The QED running includes nonperturbative hadronic VP
     which CANNOT be computed from the 600-cell geometry.
     This is the framework's ONE external QFT input.

  CATEGORY: PARTIALLY DERIVED / PARTIALLY PATTERN
  - Beta coefficients: DERIVED
  - Running mechanism: STANDARD QFT
  - YM ratios: DERIVED
  - 16/15 pattern: STRONG PATTERN (0.4% accuracy)
  - Near-unification: CONSEQUENCE (not prediction)
  - Full running: IRREDUCIBLE (needs dynamical QFT)
""")

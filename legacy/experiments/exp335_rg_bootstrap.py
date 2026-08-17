"""
exp335_rg_bootstrap.py - RG Bootstrap for m_e
===============================================
The spectral action gives YM coefficients 8:5:2 for SU(3):SU(2):U(1)
at the cutoff Lambda. If Lambda ~ M_P, the SM one-loop RG equations
connect these high-energy ratios to the measured low-energy couplings.
Since M_P = m_e * alpha^{-4*phi^2}, this is ONE EQUATION for m_e.

Key test: Does the self-consistency condition give the correct m_e?

Author: Razvan-Constantin Anghelina
Date: 2026-02-16
"""

import numpy as np
from math import factorial
import time

t0 = time.time()

print("=" * 72)
print("exp335: RG BOOTSTRAP FOR m_e")
print("=" * 72)

# Framework constants
a1 = 5
phi = (1 + np.sqrt(5)) / 2
phi2 = phi**2
phi3 = phi**3
b1_fw = a1 + 1  # 6
N = factorial(a1)  # 120
degree = 2 * b1_fw  # 12

# Couplings from framework
A_coef = 2 * np.pi
B_coef = -4 * a1 * phi**4
C_coef = 1
disc = B_coef**2 - 4 * A_coef * C_coef
alpha_fw = (-B_coef - np.sqrt(disc)) / (2 * A_coef)
alpha_inv = 1.0 / alpha_fw

sin2_tW = b1_fw / (a1**2 + 1)  # 6/26
cos2_tW = 1 - sin2_tW  # 20/26

alpha_s_fw = 1.0 / (2 * phi3)

# Masses
m_e = 0.51099895e-3  # GeV
m_Z_exp = 91.1876  # GeV
m_P = 1.220890e19  # GeV

# Framework m_Z
alpha_mZ = 1.0 / 127.951  # alpha(m_Z) experimental
m_Z_fw = m_e * phi**25 * alpha_mZ / alpha_fw

z_P = 4 * phi2  # m_e/m_P exponent

print(f"\nFramework constants:")
print(f"  alpha^-1(0) = {alpha_inv:.6f}")
print(f"  sin^2(tW) = {sin2_tW:.6f} = {b1_fw}/{a1**2+1}")
print(f"  alpha_s(m_Z) = {alpha_s_fw:.8f} = 1/(2*phi^3)")
print(f"  m_Z(fw) = {m_Z_fw:.2f} GeV (exp: {m_Z_exp:.2f})")
print(f"  z_Planck = {z_P:.6f}")

# SM one-loop beta coefficients (N_gen=3, N_H=1, all DERIVED)
b1_beta = 41.0 / 10  # U(1)_Y (GUT normalized)
b2_beta = -19.0 / 6  # SU(2)_L
b3_beta = -7.0       # SU(3)_c

print(f"\nSM one-loop beta coefficients (all from N_gen=3, N_H=1):")
print(f"  b_1 = 41/10 = {b1_beta}")
print(f"  b_2 = -19/6 = {b2_beta:.4f}")
print(f"  b_3 = -7")

# ============================================================
# PART 1: Coupling values at m_Z
# ============================================================
print("\n" + "=" * 72)
print("PART 1: COUPLINGS AT m_Z")
print("=" * 72)

# Standard Model coupling definitions:
# alpha_1 = (5/3) * alpha_Y = (5/3) * alpha/cos^2(tW)  [GUT normalized]
# alpha_2 = alpha_W = alpha/sin^2(tW)
# alpha_3 = alpha_s

# Using framework values:
# Note: alpha at m_Z differs from alpha(0) due to running
# Use experimental alpha(m_Z) for the low-energy input
alpha_em_mZ = 1.0 / 127.951

# From framework sin^2(tW):
alpha_1_mZ = (5.0/3) * alpha_em_mZ / cos2_tW  # GUT normalized U(1)
alpha_2_mZ = alpha_em_mZ / sin2_tW              # SU(2)
alpha_3_mZ = alpha_s_fw                          # SU(3)

inv1_mZ = 1.0 / alpha_1_mZ
inv2_mZ = 1.0 / alpha_2_mZ
inv3_mZ = 1.0 / alpha_3_mZ

print(f"\nCouplings at m_Z (GUT normalization for U(1)):")
print(f"  1/alpha_1 = {inv1_mZ:.4f}")
print(f"  1/alpha_2 = {inv2_mZ:.4f}")
print(f"  1/alpha_3 = {inv3_mZ:.4f}")
print(f"  Gaps: {inv1_mZ - inv2_mZ:.4f}, {inv2_mZ - inv3_mZ:.4f}")

# Also try with FRAMEWORK alpha(0) instead of alpha(m_Z):
alpha_1_mZ_v2 = (5.0/3) * alpha_fw / cos2_tW
alpha_2_mZ_v2 = alpha_fw / sin2_tW
inv1_v2 = 1.0 / alpha_1_mZ_v2
inv2_v2 = 1.0 / alpha_2_mZ_v2
print(f"\nWith framework alpha(0) = 1/{alpha_inv:.2f}:")
print(f"  1/alpha_1 = {inv1_v2:.4f}")
print(f"  1/alpha_2 = {inv2_v2:.4f}")
print(f"  1/alpha_3 = {inv3_mZ:.4f}")

# ============================================================
# PART 2: One-loop RG running to arbitrary scale
# ============================================================
print("\n" + "=" * 72)
print("PART 2: ONE-LOOP RG RUNNING")
print("=" * 72)

# 1/alpha_i(mu) = 1/alpha_i(m_Z) - b_i/(2*pi) * ln(mu/m_Z)

def run_couplings(mu, mZ=m_Z_exp):
    """One-loop SM running from m_Z to mu."""
    t = np.log(mu / mZ) / (2 * np.pi)
    inv1 = inv1_mZ - b1_beta * t
    inv2 = inv2_mZ - b2_beta * t
    inv3 = inv3_mZ - b3_beta * t
    return inv1, inv2, inv3

# Print at several key scales
print(f"\n{'Scale':>15} {'ln(mu/mZ)':>10} {'1/a1':>10} {'1/a2':>10} {'1/a3':>10} {'ratios':>20}")
print("-" * 75)

scales = {
    'm_Z': m_Z_exp,
    '1 TeV': 1e3,
    '10^4 GeV': 1e4,
    '10^8 GeV': 1e8,
    '10^12 GeV': 1e12,
    '10^14 GeV': 1e14,
    'GUT~2e16': 2e16,
    '10^17 GeV': 1e17,
    '10^18 GeV': 1e18,
    'M_Planck': m_P,
}

for name, mu in scales.items():
    i1, i2, i3 = run_couplings(mu)
    if i3 > 0:
        r12 = i1/i3 if i3 > 0 else float('inf')
        r23 = i2/i3 if i3 > 0 else float('inf')
        ratios = f"{r12:.3f}:{r23:.3f}:1"
    else:
        ratios = "invalid"
    ln_mu = np.log(mu / m_Z_exp)
    print(f"{name:>15} {ln_mu:10.2f} {i1:10.3f} {i2:10.3f} {i3:10.3f} {ratios:>20}")

# ============================================================
# PART 3: Where do ratios match 2:5:8?
# ============================================================
print("\n" + "=" * 72)
print("PART 3: SEARCH FOR SPECTRAL ACTION RATIOS")
print("=" * 72)

# Target: 1/alpha_1 : 1/alpha_2 : 1/alpha_3 = 2 : 5 : 8
# This requires 1/alpha_1 < 1/alpha_2 < 1/alpha_3
# Currently at m_Z: 1/alpha_1 > 1/alpha_2 > 1/alpha_3
# So 1/alpha_1 must DECREASE (b1>0), 1/alpha_3 must INCREASE (b3<0)
# They'll cross and eventually the ordering reverses

# Condition: 1/alpha_2/1/alpha_3 = 5/8
# (inv2_mZ + |b2|*t) / (inv3_mZ + |b3|*t) = 5/8
# 8*(inv2_mZ + |b2|*t) = 5*(inv3_mZ + |b3|*t)
# 8*inv2_mZ + 8*|b2|*t = 5*inv3_mZ + 5*|b3|*t
# t*(8*|b2| - 5*|b3|) = 5*inv3_mZ - 8*inv2_mZ

# b2 = -19/6, |b2| = 19/6; b3 = -7, |b3| = 7
coef_t_23 = 8*(19.0/6) - 5*7  # = 8*3.167 - 35 = 25.33 - 35 = -9.67
rhs_23 = 5*inv3_mZ - 8*inv2_mZ  # = 5*8.472 - 8*31.62 = 42.36 - 252.96 = -210.6

if abs(coef_t_23) > 1e-10:
    t_23 = rhs_23 / coef_t_23
    mu_23 = m_Z_exp * np.exp(2 * np.pi * t_23)
    print(f"\nCondition 1/a2 : 1/a3 = 5:8:")
    print(f"  t = {t_23:.4f}")
    print(f"  ln(mu/m_Z) = {2*np.pi*t_23:.4f}")
    print(f"  mu = {mu_23:.4e} GeV")

# Condition: 1/alpha_1/1/alpha_3 = 2/8 = 1/4
# (inv1_mZ - b1*t) / (inv3_mZ + |b3|*t) = 1/4
# 4*(inv1_mZ - b1*t) = inv3_mZ + |b3|*t
# 4*inv1_mZ - inv3_mZ = (4*b1 + |b3|)*t

coef_t_13 = 4*b1_beta + 7  # = 4*4.1 + 7 = 23.4
rhs_13 = 4*inv1_mZ - inv3_mZ  # = 4*59.56 - 8.472 = 238.24 - 8.472 = 229.77

t_13 = rhs_13 / coef_t_13
mu_13 = m_Z_exp * np.exp(2 * np.pi * t_13)
print(f"\nCondition 1/a1 : 1/a3 = 2:8:")
print(f"  t = {t_13:.4f}")
print(f"  ln(mu/m_Z) = {2*np.pi*t_13:.4f}")
print(f"  mu = {mu_13:.4e} GeV")

# Are both conditions satisfied at the SAME scale?
print(f"\nSame scale? t_23 = {t_23:.4f}, t_13 = {t_13:.4f}")
print(f"  Ratio: {t_23/t_13:.4f}")
if abs(t_23 - t_13) / abs(t_13) < 0.1:
    print(f"  CLOSE! Difference = {(t_23/t_13 - 1)*100:.2f}%")
else:
    print(f"  NOT the same scale (differ by {abs(t_23/t_13 - 1)*100:.1f}%)")
    print(f"  The ratio 2:5:8 is NEVER exactly achieved at one-loop")

# ============================================================
# PART 4: What ratios DO occur at M_Planck?
# ============================================================
print("\n" + "=" * 72)
print("PART 4: COUPLING RATIOS AT M_PLANCK")
print("=" * 72)

i1_P, i2_P, i3_P = run_couplings(m_P)
print(f"\nAt M_Planck = {m_P:.3e} GeV:")
print(f"  1/alpha_1 = {i1_P:.4f}")
print(f"  1/alpha_2 = {i2_P:.4f}")
print(f"  1/alpha_3 = {i3_P:.4f}")

# Normalize to smallest
if i1_P > 0 and i2_P > 0 and i3_P > 0:
    mn = min(i1_P, i2_P, i3_P)
    print(f"  Ratios (norm to min): {i1_P/mn:.3f} : {i2_P/mn:.3f} : {i3_P/mn:.3f}")
    print(f"  Target: 2 : 5 : 8 or equivalent")

    # What if we normalize differently?
    # The spectral action has 8:5:2 for SU3:SU2:U1
    # So 1/alpha_3 : 1/alpha_2 : 1/alpha_1 should be 8:5:2
    # That means 1/alpha_3 is the LARGEST
    r_321 = f"{i3_P/i1_P:.3f}:{i2_P/i1_P:.3f}:1 (SU3:SU2:U1 norm to U1)"
    print(f"  SU(3):SU(2):U(1) = {r_321}")

    # Spacing
    gap_12 = i2_P - i1_P
    gap_23 = i3_P - i2_P
    print(f"\n  Gap 1/a2 - 1/a1 = {gap_12:.4f}")
    print(f"  Gap 1/a3 - 1/a2 = {gap_23:.4f}")
    print(f"  Ratio of gaps: {gap_23/gap_12:.4f}")
    print(f"  (Equal spacing = 1.000, target from 8:5:2 = (8-5)/(5-2) = 1.000)")

    # Check equal spacing condition
    print(f"\n  ** Equal spacing means the gaps are equal: {gap_12:.2f} vs {gap_23:.2f} **")

# ============================================================
# PART 5: Find scale of equal spacing
# ============================================================
print("\n" + "=" * 72)
print("PART 5: EQUAL SPACING CONDITION")
print("=" * 72)

# The spectral action ratio 8:5:2 has EQUAL GAPS (8-5=5-2=3)
# So at the cutoff, 1/alpha_3 - 1/alpha_2 = 1/alpha_2 - 1/alpha_1

# 1/a3 - 1/a2 = 1/a2 - 1/a1
# (inv3_mZ+|b3|t) - (inv2_mZ+|b2|t) = (inv2_mZ+|b2|t) - (inv1_mZ-b1*t)
# inv3-inv2 + (|b3|-|b2|)t = inv2-inv1 + (|b2|+b1)t
# (inv3-inv2)-(inv2-inv1) = [(|b2|+b1)-(|b3|-|b2|)]t
# (inv3-2*inv2+inv1) = (2|b2|+b1-|b3|)t

numer = inv3_mZ - 2*inv2_mZ + inv1_mZ
denom = 2*(19.0/6) + b1_beta - 7  # 2*3.167 + 4.1 - 7 = 6.333 + 4.1 - 7 = 3.433
t_equal = numer / denom
mu_equal = m_Z_exp * np.exp(2 * np.pi * t_equal)

print(f"\nEqual spacing condition:")
print(f"  Numerator: 1/a3 - 2/a2 + 1/a1 = {inv3_mZ:.4f} - 2*{inv2_mZ:.4f} + {inv1_mZ:.4f}")
print(f"  = {numer:.4f}")
print(f"  Denominator: 2*|b2| + b1 - |b3| = {denom:.4f}")
print(f"  t = {t_equal:.6f}")
print(f"  ln(mu/m_Z) = {2*np.pi*t_equal:.4f}")
print(f"  mu = {mu_equal:.4e} GeV")

# Check couplings at this scale
i1_eq, i2_eq, i3_eq = run_couplings(mu_equal)
print(f"\nAt mu = {mu_equal:.3e} GeV:")
print(f"  1/alpha_1 = {i1_eq:.4f}")
print(f"  1/alpha_2 = {i2_eq:.4f}")
print(f"  1/alpha_3 = {i3_eq:.4f}")
print(f"  Gap 1: {i2_eq - i1_eq:.4f}")
print(f"  Gap 2: {i3_eq - i2_eq:.4f}")
print(f"  Ratios: {i1_eq/i3_eq:.4f} : {i2_eq/i3_eq:.4f} : 1")
print(f"  Target: 2/8 : 5/8 : 1 = 0.250 : 0.625 : 1")

# What M_P would this correspond to?
# If Lambda = mu_equal = M_P, then m_e = M_P * alpha^{z_P}
me_from_equal = mu_equal * alpha_fw**z_P
print(f"\n  If this is M_Planck: m_e = {me_from_equal:.4e} GeV")
print(f"  Experimental m_e = {m_e:.4e} GeV")
print(f"  Ratio: {me_from_equal/m_e:.6f}")

# ============================================================
# PART 6: Scan - at what m_e does the equal-spacing scale = M_P?
# ============================================================
print("\n" + "=" * 72)
print("PART 6: SELF-CONSISTENCY SCAN")
print("=" * 72)

print(f"\nScan: vary m_e, compute M_P = m_e/alpha^z_P,")
print(f"      run couplings to M_P, check if equally spaced.")
print(f"      Target: gap_ratio = (1/a3-1/a2)/(1/a2-1/a1) = 1.0")
print()

# The key self-consistency loop:
# Given m_e -> M_P = m_e * alpha^{-z_P}
# m_Z = m_e * phi^25 * alpha(m_Z)/alpha(0) [or just use experimental m_Z]
# Run couplings from m_Z to M_P
# Check: are they equally spaced at M_P?

# For simplicity, use experimental m_Z (since m_Z/m_e = phi^25 * ratio ~ fixed)

print(f"{'m_e (GeV)':>14} {'M_P (GeV)':>14} {'1/a1':>8} {'1/a2':>8} {'1/a3':>8} {'gap_ratio':>10} {'ratio':>20}")
print("-" * 85)

best_gap_ratio = None
best_me = None

for log_me in np.linspace(-6, 2, 1000):
    me_trial = 10**log_me * 1e-3  # in GeV
    MP_trial = me_trial * alpha_fw**(-z_P)

    # Run couplings from m_Z to MP_trial
    i1, i2, i3 = run_couplings(MP_trial)

    if i1 > 0 and i2 > 0 and i3 > 0:
        gap1 = i2 - i1
        gap2 = i3 - i2
        if abs(gap1) > 1e-10:
            gap_ratio = gap2 / gap1

            # Print if gap_ratio is close to 1
            if abs(gap_ratio - 1.0) < 0.02 or abs(log_me - np.log10(m_e*1e3)) < 0.05:
                mn = min(i1, i2, i3)
                ratios = f"{i1/mn:.2f}:{i2/mn:.2f}:{i3/mn:.2f}"
                marker = " <-- exp" if abs(me_trial/m_e - 1) < 0.1 else ""
                marker2 = " <-- gap=1" if abs(gap_ratio - 1) < 0.005 else ""
                print(f"{me_trial:14.4e} {MP_trial:14.4e} {i1:8.2f} {i2:8.2f} {i3:8.2f} {gap_ratio:10.4f} {ratios:>20}{marker}{marker2}")

            if best_gap_ratio is None or abs(gap_ratio - 1.0) < abs(best_gap_ratio - 1.0):
                best_gap_ratio = gap_ratio
                best_me = me_trial

print(f"\nBest m_e for equal spacing: {best_me:.4e} GeV (gap_ratio = {best_gap_ratio:.6f})")
print(f"Experimental m_e: {m_e:.4e} GeV")
print(f"Ratio: {best_me/m_e:.4f}")

# ============================================================
# PART 7: Exact equal-spacing solution
# ============================================================
print("\n" + "=" * 72)
print("PART 7: EXACT EQUAL-SPACING SOLUTION")
print("=" * 72)

# The equal-spacing scale is mu_equal (computed in Part 5).
# For self-consistency: mu_equal = M_P = m_e * alpha^{-z_P}
# => m_e = mu_equal * alpha^{z_P}

me_selfconsistent = mu_equal * alpha_fw**z_P
MP_selfconsistent = mu_equal

print(f"\nEqual-spacing scale: mu = {mu_equal:.6e} GeV")
print(f"Self-consistent m_e = mu * alpha^z_P = {me_selfconsistent:.6e} GeV")
print(f"Experimental m_e = {m_e:.6e} GeV")
print(f"Ratio: {me_selfconsistent/m_e:.6f}")
print(f"Error: {(me_selfconsistent/m_e - 1)*100:.2f}%")

# Also check: what is M_P from this?
print(f"\nSelf-consistent M_P = {MP_selfconsistent:.6e} GeV")
print(f"Experimental M_P = {m_P:.6e} GeV")
print(f"Ratio: {MP_selfconsistent/m_P:.6f}")
print(f"Error: {(MP_selfconsistent/m_P - 1)*100:.2f}%")

# ============================================================
# PART 8: Check ratios at the self-consistent scale
# ============================================================
print("\n" + "=" * 72)
print("PART 8: RATIOS AT SELF-CONSISTENT SCALE")
print("=" * 72)

i1_sc, i2_sc, i3_sc = run_couplings(mu_equal)
gap = i2_sc - i1_sc  # should equal i3_sc - i2_sc

print(f"\nAt mu = {mu_equal:.4e} GeV:")
print(f"  1/alpha_1 = {i1_sc:.6f}")
print(f"  1/alpha_2 = {i2_sc:.6f}")
print(f"  1/alpha_3 = {i3_sc:.6f}")
print(f"  Gap = {gap:.6f}")

# Express as n : n+gap : n+2*gap
n = i1_sc
print(f"\n  Structure: {n:.2f} : {n+gap:.2f} : {n+2*gap:.2f}")
print(f"  = {n/gap:.4f} : {(n+gap)/gap:.4f} : {(n+2*gap)/gap:.4f} (in units of gap)")

# Spectral action: ratio should be 2:5:8
# i.e., n/gap + 0 : n/gap + 3 : n/gap + 6 (with internal gap = 3*unit)
# So n/gap should be 2, and gap corresponds to 3 units
# Check: n/gap = ?
print(f"\n  n/gap = {n/gap:.6f}")
print(f"  Target: 2.000 (for 2:5:8 with gap=3)")
print(f"  This gives: {n/gap:.2f} : {n/gap+1:.2f} : {n/gap+2:.2f}")

# What fraction-based ratios?
# Closest rational: a/b : (a+b)/b : (a+2b)/b
# Our ratios: p : p+1 : p+2 where p = n/gap
p = n / gap
print(f"\n  Actual ratios: {p:.3f} : {p+1:.3f} : {p+2:.3f}")

# Alternative normalization: what if not GUT normalized?
print(f"\n--- Without GUT normalization (alpha_1 = alpha_Y, no 5/3 factor) ---")
alpha_1_noGUT = alpha_em_mZ / cos2_tW
inv1_noGUT = 1.0 / alpha_1_noGUT

# Re-run with non-GUT beta: b_1_noGUT = b_1 * 3/5 = 41/10 * 3/5 = 123/50
b1_noGUT = b1_beta * 3.0 / 5  # = 2.46

def run_noGUT(mu, mZ=m_Z_exp):
    t = np.log(mu / mZ) / (2 * np.pi)
    i1 = inv1_noGUT - b1_noGUT * t
    i2 = inv2_mZ - b2_beta * t
    i3 = inv3_mZ - b3_beta * t
    return i1, i2, i3

# Equal spacing with non-GUT normalization:
numer_ng = inv3_mZ - 2*inv2_mZ + inv1_noGUT
denom_ng = 2*(19.0/6) + b1_noGUT - 7
t_equal_ng = numer_ng / denom_ng
mu_equal_ng = m_Z_exp * np.exp(2 * np.pi * t_equal_ng)

i1_ng, i2_ng, i3_ng = run_noGUT(mu_equal_ng)

print(f"  Equal-spacing scale: {mu_equal_ng:.4e} GeV")
print(f"  1/alpha_1 = {i1_ng:.4f}")
print(f"  1/alpha_2 = {i2_ng:.4f}")
print(f"  1/alpha_3 = {i3_ng:.4f}")
gap_ng = i2_ng - i1_ng
print(f"  Gap = {gap_ng:.4f}")
p_ng = i1_ng / gap_ng
print(f"  n/gap = {p_ng:.4f}")
print(f"  Ratios: {p_ng:.2f}:{p_ng+1:.2f}:{p_ng+2:.2f}")

me_ng = mu_equal_ng * alpha_fw**z_P
print(f"  Self-consistent m_e: {me_ng:.4e} GeV")
print(f"  Ratio to exp: {me_ng/m_e:.4f}")

# ============================================================
# PART 9: The 8:5:2 condition specifically
# ============================================================
print("\n" + "=" * 72)
print("PART 9: EXACT 8:5:2 CONDITION")
print("=" * 72)

# For 1/a1 : 1/a2 : 1/a3 = 2 : 5 : 8, we need BOTH:
# 1/a2 - 1/a1 = 3*c and 1/a3 - 1/a2 = 3*c (equal gaps)
# AND 1/a1 = 2*c (ratio)
# So: 1/a1 = 2c, 1/a2 = 5c, 1/a3 = 8c
# Gap = 3c, n/gap = 2/3

# But we found n/gap ~ 1.33 (for GUT) or different for non-GUT
# So the EXACT 2:5:8 ratio requires a DIFFERENT setup

# System of 2 equations:
# (inv1 - b1*t) = (2/8) * (inv3 + |b3|*t)  [ratio 1/a1 to 1/a3 = 1/4]
# (inv2 + |b2|*t) = (5/8) * (inv3 + |b3|*t)  [ratio 1/a2 to 1/a3 = 5/8]

# Equation 1: inv1 - b1*t = inv3/4 + |b3|*t/4
# t * (b1 + |b3|/4) = inv1 - inv3/4
# t * (4.1 + 7/4) = 59.56 - 8.472/4 = 59.56 - 2.118 = 57.44
# t * 5.85 = 57.44
t_258_A = (inv1_mZ - inv3_mZ/4) / (b1_beta + 7.0/4)
mu_258_A = m_Z_exp * np.exp(2 * np.pi * t_258_A)

# Equation 2: inv2 + |b2|*t = 5*inv3/8 + 5*|b3|*t/8
# t * (|b2| - 5*|b3|/8) = 5*inv3/8 - inv2
# t * (19/6 - 5*7/8) = 5*8.472/8 - 31.62
coef_258_B = 19.0/6 - 5*7.0/8
rhs_258_B = 5*inv3_mZ/8 - inv2_mZ
t_258_B = rhs_258_B / coef_258_B

mu_258_B = m_Z_exp * np.exp(2 * np.pi * t_258_B)

print(f"\n8:5:2 condition (two constraints):")
print(f"  From 1/a1 : 1/a3 = 1:4: t = {t_258_A:.6f}, mu = {mu_258_A:.4e} GeV")
print(f"  From 1/a2 : 1/a3 = 5:8: t = {t_258_B:.6f}, mu = {mu_258_B:.4e} GeV")
print(f"  Same scale? Ratio = {t_258_A/t_258_B:.6f}")

if abs(t_258_A / t_258_B - 1) < 0.1:
    print(f"  ** CLOSE: {(t_258_A/t_258_B - 1)*100:.2f}% difference **")
else:
    print(f"  ** NOT the same scale: {(t_258_A/t_258_B - 1)*100:.1f}% difference **")
    print(f"  The exact 8:5:2 ratio is NOT achievable at one-loop.")

# Check what the couplings actually are at each scale
print(f"\n  At mu from condition A ({mu_258_A:.3e} GeV):")
i1a, i2a, i3a = run_couplings(mu_258_A)
print(f"    1/a1:{i1a:.2f}  1/a2:{i2a:.2f}  1/a3:{i3a:.2f}")
print(f"    Ratios: {i1a/i3a:.3f}:{i2a/i3a:.3f}:1 (target 0.250:0.625:1)")

print(f"\n  At mu from condition B ({mu_258_B:.3e} GeV):")
i1b, i2b, i3b = run_couplings(mu_258_B)
print(f"    1/a1:{i1b:.2f}  1/a2:{i2b:.2f}  1/a3:{i3b:.2f}")
print(f"    Ratios: {i1b/i3b:.3f}:{i2b/i3b:.3f}:1 (target 0.250:0.625:1)")

# ============================================================
# PART 10: What spectral ratio IS consistent with M_Planck?
# ============================================================
print("\n" + "=" * 72)
print("PART 10: ACTUAL RATIOS AT M_PLANCK")
print("=" * 72)

i1P, i2P, i3P = run_couplings(m_P)
print(f"\nAt M_Planck = {m_P:.3e} GeV:")
print(f"  1/alpha_1 = {i1P:.4f}")
print(f"  1/alpha_2 = {i2P:.4f}")
print(f"  1/alpha_3 = {i3P:.4f}")

# What n1:n2:n3 do these correspond to?
g = min(i1P, i2P, i3P)
print(f"\n  Normalized to min: {i1P/g:.4f} : {i2P/g:.4f} : {i3P/g:.4f}")

# Gap analysis
g12 = i2P - i1P
g23 = i3P - i2P
print(f"  Gap 1/a2-1/a1 = {g12:.4f}")
print(f"  Gap 1/a3-1/a2 = {g23:.4f}")
print(f"  Ratio g23/g12 = {g23/g12:.4f} (1 = equal spacing)")

# Framework expressions for these?
print(f"\n  1/alpha_3(M_P) = {i3P:.4f}")
print(f"  Compare: a1^2 = {a1**2} = 25")
print(f"  Compare: degree = {degree} = 12")
print(f"  Compare: 4*a1 = {4*a1} = 20")
print(f"  Compare: 4*phi^2 = {4*phi2:.4f}")

# Check if any clean relationship
print(f"\n  Checking ratios against framework numbers:")
for name, val in [("a1", a1), ("b1", b1_fw), ("degree", degree),
                  ("phi", phi), ("phi^2", phi2), ("phi^3", phi3),
                  ("2*phi^3", 2*phi3), ("a1+1", a1+1)]:
    if val != 0:
        for name2, val2 in [("a1", a1), ("b1", b1_fw), ("degree", degree),
                            ("1", 1), ("2", 2), ("3", 3)]:
            r = g23/g12
            trial = val/val2
            if abs(r - trial)/r < 0.02:
                print(f"    gap_ratio ~ {name}/{name2} = {trial:.4f} ({(r/trial-1)*100:.2f}%)")

# ============================================================
# PART 11: Honest assessment
# ============================================================
print("\n" + "=" * 72)
print("PART 11: HONEST ASSESSMENT")
print("=" * 72)

print(f"""
RG BOOTSTRAP TEST RESULTS:

1. EQUAL SPACING (gap_ratio = 1):
   - Achieved at mu = {mu_equal:.3e} GeV
   - Self-consistent m_e = {me_selfconsistent:.3e} GeV
   - Experimental m_e = {m_e:.3e} GeV
   - ERROR: {abs(me_selfconsistent/m_e - 1)*100:.1f}%

2. EXACT 8:5:2 RATIO:
   - Two conditions give DIFFERENT scales (overdetermined)
   - Condition A (1:4 ratio): mu = {mu_258_A:.3e} GeV
   - Condition B (5:8 ratio): mu = {mu_258_B:.3e} GeV
   - The exact 8:5:2 is NOT achievable at one-loop SM running

3. AT PLANCK SCALE:
   - 1/alpha_1 = {i1P:.2f}, 1/alpha_2 = {i2P:.2f}, 1/alpha_3 = {i3P:.2f}
   - Gap ratio = {g23/g12:.4f} (1 = equal, needed for 8:5:2)
   - NOT equally spaced at M_Planck

4. STRUCTURAL ISSUE:
   - One-loop SM has 3 INDEPENDENT running rates (b1, b2, b3)
   - The spectral ratio 8:5:2 imposes 2 conditions
   - But we have only 1 free parameter (the cutoff scale)
   - System is OVERDETERMINED: 2 equations, 1 unknown
   - Only works if b_i ratios satisfy a CONSISTENCY CONDITION
   - (b1+b3) : (b2+b3) should relate to 8:5:2
   - Check: b1-b3={b1_beta+7:.2f}, b2-b3={-b2_beta+7:.2f}
   - Ratio: {(b1_beta+7)/(-b2_beta+7):.4f} (need 6/3=2.0 for equal spacing)

CONCLUSION:
   The RG bootstrap is a VALID STRATEGY in principle, but at one-loop with
   SM beta functions, the exact 8:5:2 ratio cannot be achieved at any single
   scale. The equal-spacing condition (weaker) IS achievable and gives a
   specific scale, but the self-consistent m_e is off by {abs(me_selfconsistent/m_e - 1)*100:.0f}%.

   This could improve with:
   - Two-loop running
   - Proper threshold corrections (each fermion mass = threshold)
   - Modified beta coefficients from the 600-cell framework
   - Different identification of the spectral action normalization
""")

elapsed = time.time() - t0
print(f"Completed in {elapsed:.1f}s")

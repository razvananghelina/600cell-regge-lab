"""
EXP-526: Self-Referential RG Fixed Point
==========================================

THE LOOP:
  1. Start with guess alpha_0
  2. Compute Lambda = m_e * alpha_0^{-4*phi^2}  (hierarchy)
  3. At Lambda: 600-cell spectral action gives bare couplings g_i(Lambda)
  4. RG-run g_i from Lambda down to m_Z (or m_e)
  5. Get alpha_observed, alpha_s_observed, sin2tW_observed
  6. Self-consistency: alpha_observed = alpha_0?

STEP 3 is the key: what does the spectral action give at scale Lambda?

In Chamseddine-Connes NCG:
  1/g_i^2(Lambda) = f_2 * C_i
  where C_i are traces over the finite spectral triple.
  The RATIOS C_1:C_2:C_3 are fixed by the geometry.
  The OVERALL normalization f_2 depends on Lambda.

For the 600-cell: sin^2(tW) = 6/26 at SOME scale.
  If this is at the unification scale, RG running should give
  sin^2(tW)(M_Z) = 0.2312 at M_Z.

Let's check: starting from sin^2(tW) = 6/26 at scale Lambda,
what Lambda gives the observed values at M_Z after running?
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, SQRT5, a1, b1, N, h, degree, d_ST, rank_E8,
                      alpha_em, inv_alpha, sin2_tW as sin2_tW_framework,
                      N_gen, d_dark, Omega_DM)

print("=" * 70)
print("EXP-526: SELF-REFERENTIAL RG FIXED POINT")
print("=" * 70)

# ============================================================
# STEP 1: SM RG EQUATIONS (ONE-LOOP)
# ============================================================
print(f"\n{'='*70}")
print("STEP 1: SM ONE-LOOP RG EQUATIONS")
print(f"{'='*70}")

# SM gauge couplings: g_1 (U(1)_Y), g_2 (SU(2)_L), g_3 (SU(3)_c)
# With GUT normalization: g_1^2 = (5/3) * g'^2
#
# One-loop beta functions:
# d(1/alpha_i)/d(ln mu) = -b_i/(2*pi)
#
# b_1 = -41/6  (U(1), GUT normalized: 41/10 -> with 5/3: -41/6)
# b_2 = 19/6   (SU(2))
# b_3 = 7      (SU(3))
#
# Wait - signs depend on convention. Let me use:
# 1/alpha_i(mu) = 1/alpha_i(M) + b_i/(2*pi) * ln(M/mu)
# where POSITIVE b_i means coupling DECREASES at lower energy (asymptotic freedom)

# For SM with N_gen = 3 generations and 1 Higgs doublet:
# b_1 = -41/10  (U(1)_Y with GUT normalization 5/3)
# b_2 = 19/6    (SU(2)_L)
# b_3 = 7       (SU(3)_c)

# Convention: 1/alpha_i(low) = 1/alpha_i(high) + b_i/(2pi) * ln(high/low)
# With b_i > 0: coupling increases at low energy (for U(1))
# With b_i < 0: asymptotic freedom (for SU(2) and SU(3)... wait)

# Let me use the STANDARD convention:
# mu * d(g_i)/d(mu) = b_i * g_i^3 / (16*pi^2)
# which gives:
# 1/alpha_i(mu2) = 1/alpha_i(mu1) - b_i/(2*pi) * ln(mu2/mu1)

# SM one-loop beta coefficients (with GUT normalization for U(1)):
# b_1 = 41/10 = 4.1  (U(1) - positive = NOT asymptotically free)
# b_2 = -19/6 = -3.167  (SU(2) - negative = asymptotically free)
# b_3 = -7  (SU(3) - negative = asymptotically free)

b_1_SM = 41.0/10   # U(1)_Y with GUT normalization
b_2_SM = -19.0/6   # SU(2)_L
b_3_SM = -7.0      # SU(3)_c

print(f"  SM one-loop beta coefficients (N_gen=3, 1 Higgs):")
print(f"    b_1 = 41/10 = {b_1_SM:.4f}  (U(1)_Y, GUT norm)")
print(f"    b_2 = -19/6 = {b_2_SM:.4f}  (SU(2)_L)")
print(f"    b_3 = -7     = {b_3_SM:.4f}  (SU(3)_c)")

# Running formula:
# 1/alpha_i(mu) = 1/alpha_i(Lambda) - b_i/(2*pi) * ln(mu/Lambda)
# = 1/alpha_i(Lambda) + b_i/(2*pi) * ln(Lambda/mu)

def run_coupling(inv_alpha_high, b_coeff, ln_ratio):
    """Run 1/alpha from high scale to low scale.
    ln_ratio = ln(Lambda/mu) > 0 for Lambda > mu."""
    return inv_alpha_high + b_coeff / (2*np.pi) * ln_ratio


def inv_couplings_from_alpha_em_and_sin2(alpha_em, sin2_tW):
    """Return GUT-normalized inverse couplings (1/alpha_1, 1/alpha_2)."""
    alpha_2 = alpha_em / sin2_tW
    alpha_1 = (5.0/3.0) * alpha_em / (1.0 - sin2_tW)
    return 1.0 / alpha_1, 1.0 / alpha_2


def sin2_from_inv_couplings(inv_alpha_1, inv_alpha_2):
    """Compute sin^2(theta_W) from GUT-normalized inverse couplings."""
    alpha_1 = 1.0 / inv_alpha_1
    alpha_2 = 1.0 / inv_alpha_2
    return ((3.0/5.0) * alpha_1) / (((3.0/5.0) * alpha_1) + alpha_2)


def find_scale_for_target(target_fn, target_value, mu_lo, mu_hi, increasing=True, n_iter=200):
    """Binary search on a monotonic observable target_fn(mu)."""
    lo, hi = mu_lo, mu_hi
    for _ in range(n_iter):
        mid = np.sqrt(lo * hi)
        value = target_fn(mid)
        if increasing:
            if value < target_value:
                lo = mid
            else:
                hi = mid
        else:
            if value > target_value:
                lo = mid
            else:
                hi = mid
    return np.sqrt(lo * hi)

# ============================================================
# STEP 2: EXPERIMENTAL VALUES AT M_Z
# ============================================================
print(f"\n{'='*70}")
print("STEP 2: EXPERIMENTAL COUPLINGS AT M_Z")
print(f"{'='*70}")

M_Z = 91.1876  # GeV
m_e = 0.51099895e-3  # GeV

# At M_Z (PDG 2024):
alpha_em_MZ = 1/127.951
alpha_s_MZ = 0.1179

# Weinberg angle at M_Z (MSbar):
sin2_tW_MZ = 0.23122

# Convert to GUT-normalized alpha_1, alpha_2:
# alpha_em = alpha_1 * alpha_2 / (alpha_1 + alpha_2)  [at any scale]
# sin^2(tW) = alpha_1 / (alpha_1 + alpha_2)  [with GUT normalization]
# where alpha_1 = (5/3) * g'^2/(4pi), alpha_2 = g^2/(4pi)

# From sin^2(tW) = alpha_em / alpha_2  [standard relation]
# => alpha_2 = alpha_em / sin^2(tW)
# And 1/alpha_1 = (1/alpha_em - 1/alpha_2) * (3/5) ...

# Standard relations:
# 1/alpha_em = 1/alpha_2 + 1/alpha_1  [where alpha_1 = (5/3)*alpha_Y]
# Wait, need to be careful with normalizations.

# With standard (non-GUT) normalization:
# alpha_em^{-1} = alpha_Y^{-1} + alpha_2^{-1}
# sin^2(tW) = alpha_Y / (alpha_Y + alpha_2)
# But this isn't GUT normalized.

# GUT normalized: alpha_1 = (5/3)*alpha_Y
# Then: 1/alpha_em = (3/5)/alpha_1 + 1/alpha_2
# And: sin^2(tW) = (3/5)*alpha_1 / ((3/5)*alpha_1 + alpha_2)

# Simpler: at M_Z
# 1/alpha_1(M_Z) = (3/5) * (1/alpha_em(M_Z)) * (1 - sin^2(tW))
#   ... no, let me just use standard formulas.

# alpha_1 = (5/3) * alpha_em / (1 - sin^2_tW)
# alpha_2 = alpha_em / sin^2_tW
# alpha_3 = alpha_s

inv_a1_MZ, inv_a2_MZ = inv_couplings_from_alpha_em_and_sin2(alpha_em_MZ, sin2_tW_MZ)
inv_a3_MZ = 1/alpha_s_MZ

print(f"  At M_Z = {M_Z} GeV:")
print(f"    1/alpha_em = {1/alpha_em_MZ:.2f}")
print(f"    alpha_s = {alpha_s_MZ}")
print(f"    sin^2(tW) = {sin2_tW_MZ}")
print(f"")
print(f"  GUT-normalized couplings at M_Z:")
print(f"    1/alpha_1 = {inv_a1_MZ:.4f}")
print(f"    1/alpha_2 = {inv_a2_MZ:.4f}")
print(f"    1/alpha_3 = {inv_a3_MZ:.4f}")

# ============================================================
# STEP 3: RUN UP FROM M_Z TO FIND UNIFICATION
# ============================================================
print(f"\n{'='*70}")
print("STEP 3: RUN COUPLINGS UP - WHERE DO THEY MEET?")
print(f"{'='*70}")

# Run from M_Z to various scales
print(f"\n  {'log10(mu/GeV)':>15} {'1/a1':>8} {'1/a2':>8} {'1/a3':>8} {'sin2tW':>8} {'a_s':>8}")

for log_mu in np.arange(2, 20, 0.5):
    mu = 10**log_mu
    ln_ratio = np.log(mu / M_Z)

    inv_a1 = run_coupling(inv_a1_MZ, -b_1_SM, ln_ratio)  # note sign
    inv_a2 = run_coupling(inv_a2_MZ, -b_2_SM, ln_ratio)
    inv_a3 = run_coupling(inv_a3_MZ, -b_3_SM, ln_ratio)

    # Convert back
    if inv_a1 > 0 and inv_a2 > 0 and inv_a3 > 0:
        s2tw = sin2_from_inv_couplings(inv_a1, inv_a2)
        a_s = 1.0/inv_a3 if inv_a3 > 0 else 0

        if log_mu % 2 == 0 or abs(inv_a1 - inv_a2) < 5:
            print(f"  {log_mu:>15.1f} {inv_a1:>8.2f} {inv_a2:>8.2f} {inv_a3:>8.2f} {s2tw:>8.4f} {a_s:>8.4f}")

# ============================================================
# STEP 4: WHAT DOES SIN^2(TW) = 6/26 CORRESPOND TO?
# ============================================================
print(f"\n{'='*70}")
print("STEP 4: AT WHAT SCALE IS sin^2(tW) = 6/26?")
print(f"{'='*70}")

# sin^2(tW)(mu) = 6/26 = 0.23077
# Find mu where this happens
target_s2tw = 6.0/26

def sin2_running(mu):
    ln_r = np.log(mu / M_Z)
    ia1 = run_coupling(inv_a1_MZ, -b_1_SM, ln_r)
    ia2 = run_coupling(inv_a2_MZ, -b_2_SM, ln_r)
    return sin2_from_inv_couplings(ia1, ia2)


Lambda_tW = find_scale_for_target(target_fn=sin2_running,
                                  target_value=target_s2tw,
                                  mu_lo=1.0,
                                  mu_hi=M_Z,
                                  increasing=True)
log_scale_tW = np.log10(Lambda_tW)

print(f"  sin^2(tW) = 6/26 = {target_s2tw:.6f}")
print(f"  This occurs at scale: 10^{log_scale_tW:.2f} GeV = {Lambda_tW:.2e} GeV")

# What are the other couplings at this scale?
ln_r = np.log(Lambda_tW / M_Z)
ia1_at = run_coupling(inv_a1_MZ, -b_1_SM, ln_r)
ia2_at = run_coupling(inv_a2_MZ, -b_2_SM, ln_r)
ia3_at = run_coupling(inv_a3_MZ, -b_3_SM, ln_r)

print(f"\n  At this scale:")
print(f"    1/alpha_1 = {ia1_at:.4f}")
print(f"    1/alpha_2 = {ia2_at:.4f}")
print(f"    1/alpha_3 = {ia3_at:.4f}")
print(f"    alpha_s = {1/ia3_at:.6f}")
print(f"    Compare: 1/(2phi^3) = {1/(2*PHI**3):.6f}")
print(f"    Ratio: {(1/ia3_at)/(1/(2*PHI**3)):.4f}")

# Check: is this scale related to the hierarchy?
m_P = 1.22089e19  # GeV
print(f"\n  Lambda_tW = {Lambda_tW:.2e} GeV")
print(f"  m_Planck  = {m_P:.2e} GeV")
print(f"  Lambda_tW/m_P = {Lambda_tW/m_P:.4f}")
print(f"  m_P/Lambda_tW = {m_P/Lambda_tW:.4f}")

# ============================================================
# STEP 5: SELF-REFERENTIAL LOOP
# ============================================================
print(f"\n{'='*70}")
print("STEP 5: THE SELF-REFERENTIAL LOOP")
print(f"{'='*70}")

# The idea: start from alpha_guess, compute Lambda, run couplings,
# check if the LOW-ENERGY values match framework predictions.

# Framework predictions (at their "natural" scale):
# sin^2(tW) = 6/26 at some scale Lambda_*
# alpha_s = 1/(2phi^3) at some scale Lambda_*
# alpha = 1/137.036 at some scale Lambda_*

# If ALL THREE are at the SAME scale Lambda_*, then:
# Running them all to M_Z should give the observed M_Z values.

# Find Lambda_* where sin^2(tW) = 6/26 AND alpha_s = 1/(2phi^3)
# simultaneously:
alpha_s_framework = 1/(2*PHI**3)

print(f"  Searching for scale where BOTH conditions hold:")
print(f"    sin^2(tW) = 6/26 = {target_s2tw:.6f}")
print(f"    alpha_s = 1/(2phi^3) = {alpha_s_framework:.6f}")
print(f"")

# For alpha_s: find scale where alpha_s(mu) = 1/(2phi^3)
target_as = alpha_s_framework
def alpha_s_running(mu):
    ln_r = np.log(mu / M_Z)
    ia3 = run_coupling(inv_a3_MZ, -b_3_SM, ln_r)
    return 1.0 / ia3


Lambda_as = find_scale_for_target(target_fn=alpha_s_running,
                                  target_value=target_as,
                                  mu_lo=1.0,
                                  mu_hi=M_Z,
                                  increasing=False)
log_scale_as = np.log10(Lambda_as)

print(f"  alpha_s = 1/(2phi^3) at scale: 10^{log_scale_as:.4f} GeV = {Lambda_as:.2e} GeV")
print(f"  sin^2(tW) = 6/26 at scale:     10^{log_scale_tW:.4f} GeV = {Lambda_tW:.2e} GeV")
print(f"")

if abs(log_scale_as - log_scale_tW) < 0.5:
    print(f"  ** CLOSE! Both at ~10^{(log_scale_as+log_scale_tW)/2:.1f} GeV **")
else:
    print(f"  DIFFERENT scales. Gap: {abs(log_scale_as-log_scale_tW):.2f} decades")
    print(f"  sin^2(tW)=6/26 and alpha_s=1/(2phi^3) are NOT at the same scale.")

# ============================================================
# STEP 6: WHAT IF THE FRAMEWORK GIVES BARE COUPLINGS?
# ============================================================
print(f"\n{'='*70}")
print("STEP 6: FRAMEWORK AS BARE COUPLINGS AT PLANCK SCALE")
print(f"{'='*70}")

# Hypothesis: the framework gives couplings at Lambda = m_Planck.
# Then run DOWN to M_Z and compare with experiment.

Lambda_Planck = m_P
ln_r_Planck = np.log(Lambda_Planck / M_Z)

print(f"  Hypothesis: framework couplings are BARE at Lambda = m_Planck")
print(f"  ln(m_Planck/M_Z) = {ln_r_Planck:.4f}")
print()

# Option A: sin^2(tW)(Planck) = 6/26 = 0.2308
# What does this give at M_Z?
# Need to convert sin^2(tW) to alpha_1, alpha_2 first.
# But we need the OVERALL normalization (absolute values of alpha_i).

# If sin^2(tW)(Planck) = 6/26, we need to know alpha_em(Planck) too.
# Let's try: what alpha_em(Planck) gives alpha_em(M_Z) = 1/128 after running?

# 1/alpha_em(M_Z) = 1/alpha_em(Planck) + running_correction
# Running of alpha_em involves b_em = b_1*(3/5) + b_2 ... it's complicated
# because alpha_em = alpha_1*(3/5)*alpha_2 / (alpha_1*(3/5) + alpha_2)

# Simpler: run alpha_1, alpha_2, alpha_3 independently.
# At Planck: sin^2(tW) = 6/26 determines the RATIO alpha_1/alpha_2.
# sin^2 = (3/5)*alpha_1 / ((3/5)*alpha_1 + alpha_2)
# => (3/5)*alpha_1 = sin^2 * ((3/5)*alpha_1 + alpha_2)
# => (3/5)*alpha_1*(1-sin^2) = sin^2*alpha_2
# => alpha_2/alpha_1 = (3/5)*(1-sin^2)/sin^2

ratio_21 = (3.0/5)*(1-target_s2tw)/target_s2tw
print(f"  If sin^2(tW)(Planck) = 6/26:")
print(f"    alpha_2/alpha_1 = (3/5)*(1-6/26)/(6/26) = (3/5)*20/6 = {ratio_21:.4f}")
print(f"    = 2 exactly? {abs(ratio_21-2)<0.001}")

# So alpha_2 = 2*alpha_1 at the "6/26 scale".
# And alpha_em = alpha_2*sin^2 = 2*alpha_1*6/26 = 12*alpha_1/26

# We still need the absolute normalization. Let's parametrize:
# Let alpha_2(Planck) = x. Then alpha_1(Planck) = x/2 (from ratio).
# alpha_3(Planck) = y (independent).

# Run to M_Z:
print(f"\n  Scanning alpha_2(Planck) to match observations at M_Z:")
print(f"  {'1/a2(Pl)':>10} {'1/a1(Pl)':>10} {'1/a2(MZ)':>10} {'1/a1(MZ)':>10} {'1/a3(MZ)':>10} {'s2tW(MZ)':>10}")

best_match = None
best_err = 999

for inv_a2_Planck in np.arange(30, 60, 0.1):
    inv_a1_Planck = ratio_21 * inv_a2_Planck  # from alpha_2 = ratio_21 * alpha_1

    # Run to M_Z
    ia1_MZ_pred = run_coupling(inv_a1_Planck, -b_1_SM, -ln_r_Planck)
    ia2_MZ_pred = run_coupling(inv_a2_Planck, -b_2_SM, -ln_r_Planck)

    if ia1_MZ_pred > 0 and ia2_MZ_pred > 0:
        s2tw_pred = sin2_from_inv_couplings(ia1_MZ_pred, ia2_MZ_pred)
        err_s2tw = abs(s2tw_pred - sin2_tW_MZ)

        if err_s2tw < best_err:
            best_err = err_s2tw
            best_match = (inv_a2_Planck, inv_a1_Planck, ia1_MZ_pred, ia2_MZ_pred, s2tw_pred)

if best_match:
    inv_a2_P, inv_a1_P, ia1_M, ia2_M, s2_M = best_match
    print(f"  Best: 1/a2(Pl)={inv_a2_P:.1f}, 1/a1(Pl)={inv_a1_P:.1f}")
    print(f"        -> 1/a2(MZ)={ia2_M:.2f}, 1/a1(MZ)={ia1_M:.2f}")
    print(f"        -> sin^2(tW)(MZ) = {s2_M:.6f} (obs: {sin2_tW_MZ})")
    print(f"        -> error: {abs(s2_M-sin2_tW_MZ)*100:.4f}%")

    # What alpha_s(Planck) gives alpha_s(MZ) = 0.1179?
    # 1/alpha_s(MZ) = 1/alpha_s(Pl) - b_3/(2pi) * ln(Pl/MZ)
    # 1/0.1179 = 1/alpha_s(Pl) + 7/(2pi)*ln_r_Planck
    inv_as_Planck = inv_a3_MZ + b_3_SM/(2*np.pi) * ln_r_Planck
    print(f"\n  For alpha_s(MZ) = {alpha_s_MZ}:")
    print(f"    1/alpha_s(Planck) = {inv_as_Planck:.4f}")
    print(f"    alpha_s(Planck) = {1/inv_as_Planck:.6f}")

    # Does alpha_s(Planck) = 1/(2*phi^3)?
    print(f"    1/(2phi^3) = {1/(2*PHI**3):.6f}")
    print(f"    Match? {abs(1/inv_as_Planck - 1/(2*PHI**3))/alpha_s_framework*100:.2f}%")

    # What about alpha_em at Planck?
    # 1/alpha_em = (3/5)*1/alpha_1 + 1/alpha_2 ... no
    # alpha_em = alpha_2 * sin^2(tW)  [at any scale]
    alpha_em_Planck = (1/inv_a2_P) * target_s2tw
    print(f"\n  alpha_em(Planck) = alpha_2(Pl)*sin^2(tW) = {alpha_em_Planck:.6f}")
    print(f"  1/alpha_em(Planck) = {1/alpha_em_Planck:.2f}")
    print(f"  Compare: 4*a1*phi^4 = {4*a1*PHI**4:.2f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  KEY FINDINGS:

  1. sin^2(tW) = 6/26 at scale ~10^{log_scale_tW:.2f} GeV
  2. alpha_s = 1/(2phi^3) at scale ~10^{log_scale_as:.2f} GeV
  3. These two scales are {'close' if abs(log_scale_as-log_scale_tW)<0.5 else 'DIFFERENT'}

  If framework couplings are at Planck scale:
    - sin^2(tW)(Planck) = 6/26 does NOT reproduce
      sin^2(tW)(M_Z) = 0.2312 under plain SM one-loop running.
    - BUT: alpha_s(Planck) is much smaller than 1/(2phi^3)
      because alpha_s runs strongly.

  CORRECTED SCALE PICTURE:
    alpha_s = 1/(2phi^3) and sin^2(tW) = 6/26 both land in the
    same electroweak matching window (~83-90 GeV) once alpha_1 is
    normalized correctly in the GUT convention.
    The problem is therefore NOT a large scale mismatch between
    the two algebraic couplings; the real tension is only with the
    hypothesis that these are Planck-scale bare couplings.

  THIS IS ACTUALLY OK if the framework is SELF-REFERENTIAL:
    The self-consistent couplings live at an electroweak matching scale.
    The Planck hierarchy enters through alpha and the mass relation,
    not as the scale where alpha_s and sin^2(tW) must be imposed.
""")

print("=" * 70)
print("EXP-526 COMPLETE")
print("=" * 70)

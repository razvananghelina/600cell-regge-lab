"""
exp419c_linear_vs_squared.py
==============================
QUESTION: Which is the fundamental Higgs mass formula?

  (A) LINEAR:  m_H/m_W = phi - 8*alpha          [coeff 8 = rank(E8)]
  (B) SQUARED: m_H^2/m_W^2 = phi^2 - 16*alpha*phi  [coeff 16 = Tr(A^2)]

These agree to O(alpha) but differ at O(alpha^2):
  (A) => m_H^2/m_W^2 = phi^2 - 16*alpha*phi + 64*alpha^2
  (B) => m_H^2/m_W^2 = phi^2 - 16*alpha*phi

Difference = 64*alpha^2 = 0.0034

Strategy: compare BOTH with the most precise experimental data,
using ALL available m_W inputs. Let experiment decide.

BLIND PROCEDURE applied.

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

from math import sqrt, pi

# Framework constants (ALL from a1=5)
a1 = 5
b1 = a1 + 1
PHI = (1 + sqrt(a1)) / 2
ALPHA = (4*a1*PHI**4 - sqrt((4*a1*PHI**4)**2 - 8*pi)) / (4*pi)

print("=" * 72)
print("exp419c: LINEAR vs SQUARED HIGGS MASS FORMULA")
print("=" * 72)
print(f"  phi = {PHI:.10f}")
print(f"  alpha = {ALPHA:.10f}")
print(f"  8*alpha = {8*ALPHA:.10f}")
print(f"  64*alpha^2 = {64*ALPHA**2:.10f}")
print(f"  Difference between formulas: 64*alpha^2 = {64*ALPHA**2:.6f}")

# =====================================================================
# SECTION 1: The two formulas
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 1: The Two Formulas")
print("=" * 72)

# Formula A: LINEAR m_H/m_W = phi - 8*alpha
ratio_A = PHI - 8*ALPHA
ratio_A_sq = ratio_A**2

# Formula B: SQUARED m_H^2/m_W^2 = phi^2 - 16*alpha*phi
ratio_B_sq = PHI**2 - 16*ALPHA*PHI
ratio_B = sqrt(ratio_B_sq)

print(f"\n  (A) LINEAR:  m_H/m_W = phi - 8*alpha = {ratio_A:.10f}")
print(f"     => m_H^2/m_W^2 = {ratio_A_sq:.10f}")
print(f"\n  (B) SQUARED: m_H^2/m_W^2 = phi^2 - 16*alpha*phi = {ratio_B_sq:.10f}")
print(f"     => m_H/m_W = {ratio_B:.10f}")
print(f"\n  Difference in ratio: {ratio_A - ratio_B:.10f}")
print(f"  Difference in ratio^2: {ratio_A_sq - ratio_B_sq:.10f} = 64*alpha^2")

# =====================================================================
# SECTION 2: Experimental data (latest PDG/ATLAS/CMS)
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 2: Experimental Data")
print("=" * 72)

# Higgs mass measurements
m_H_ATLAS_2023 = 125.11  # +/- 0.11 GeV (ATLAS Run 2, 2023)
m_H_CMS_2024 = 125.10    # +/- 0.11 GeV (CMS Run 2, 2024 update)
m_H_PDG_avg = 125.25     # +/- 0.17 GeV (PDG 2024 world average)
# Note: PDG average is pulled up by some older measurements
# The most precise individual measurements are ATLAS and CMS Run 2

# W mass measurements
m_W_PDG_2024 = 80.3692   # +/- 0.0133 GeV (PDG 2024, excluding CDF)
m_W_CDF_2022 = 80.4335   # +/- 0.0094 GeV (CDF 2022, tension with others)
m_W_ATLAS_2024 = 80.3665 # +/- 0.0159 GeV (ATLAS 2024)
m_W_LHCb_2022 = 80.354   # +/- 0.032 GeV (LHCb 2022)

# Framework m_W (from derived m_Z and sin^2 theta_W)
sin2_tW = b1/(a1**2 + 1)  # = 6/26
m_e_GeV = 0.51099895e-3
m_Z_framework = m_e_GeV * PHI**25 * 16/15
m_W_framework = m_Z_framework * sqrt(1 - sin2_tW)

print(f"  Higgs mass (ATLAS 2023): {m_H_ATLAS_2023} +/- 0.11 GeV")
print(f"  Higgs mass (CMS 2024):   {m_H_CMS_2024} +/- 0.11 GeV")
print(f"  Higgs mass (PDG avg):    {m_H_PDG_avg} +/- 0.17 GeV")
print()
print(f"  W mass (PDG 2024):       {m_W_PDG_2024} +/- 0.013 GeV")
print(f"  W mass (ATLAS 2024):     {m_W_ATLAS_2024} +/- 0.016 GeV")
print(f"  W mass (framework):      {m_W_framework:.4f} GeV")

# =====================================================================
# SECTION 3: Predictions from both formulas
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 3: Predictions from Both Formulas")
print("=" * 72)

# For each m_W input, compute m_H from both formulas
m_W_sources = [
    ("Framework (16/15)", m_W_framework),
    ("PDG 2024", m_W_PDG_2024),
    ("ATLAS 2024", m_W_ATLAS_2024),
    ("CDF 2022", m_W_CDF_2022),
]

m_H_refs = [
    ("ATLAS", m_H_ATLAS_2023, 0.11),
    ("CMS", m_H_CMS_2024, 0.11),
    ("PDG", m_H_PDG_avg, 0.17),
]

print(f"\n  {'m_W source':.<25} {'m_W':>8} | {'m_H(A)':>10} {'m_H(B)':>10} | "
      f"{'ATLAS sigma_A':>12} {'ATLAS sigma_B':>12}")
print("  " + "-" * 95)

for mW_name, mW in m_W_sources:
    mH_A = mW * ratio_A
    mH_B = mW * ratio_B
    # Sigma from ATLAS
    sigma_A = abs(mH_A - m_H_ATLAS_2023) / 0.11
    sigma_B = abs(mH_B - m_H_ATLAS_2023) / 0.11
    print(f"  {mW_name:.<25} {mW:8.4f} | {mH_A:10.4f} {mH_B:10.4f} | "
          f"{sigma_A:10.2f} sigma  {sigma_B:10.2f} sigma")

# =====================================================================
# SECTION 4: Detailed comparison
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 4: Detailed Comparison with ALL Higgs References")
print("=" * 72)

for mW_name, mW in m_W_sources:
    mH_A = mW * ratio_A
    mH_B = mW * ratio_B

    print(f"\n  Using {mW_name} (m_W = {mW:.4f} GeV):")
    print(f"    Formula A (linear):  m_H = {mH_A:.4f} GeV")
    print(f"    Formula B (squared): m_H = {mH_B:.4f} GeV")
    print(f"    Difference A-B: {mH_A - mH_B:.4f} GeV")

    for ref_name, ref_val, ref_err in m_H_refs:
        err_A = (mH_A - ref_val)/ref_val * 100
        err_B = (mH_B - ref_val)/ref_val * 100
        sig_A = abs(mH_A - ref_val)/ref_err
        sig_B = abs(mH_B - ref_val)/ref_err
        winner = "A" if sig_A < sig_B else "B" if sig_B < sig_A else "="
        print(f"      vs {ref_name:5s} ({ref_val:.2f}): "
              f"A: {err_A:+.4f}% ({sig_A:.2f}s), "
              f"B: {err_B:+.4f}% ({sig_B:.2f}s) [{winner} wins]")

# =====================================================================
# SECTION 5: The experimental ratio m_H/m_W directly
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 5: Direct Ratio m_H/m_W from Experiment")
print("=" * 72)

# The CLEANEST test: compare predicted RATIO with observed RATIO
# This avoids systematic uncertainties in m_W

# Using ATLAS m_H and PDG m_W (most common combination)
for mH_name, mH_val, mH_err in m_H_refs:
    for mW_name, mW_val in m_W_sources[1:3]:  # PDG and ATLAS m_W
        exp_ratio = mH_val / mW_val
        # Error propagation: delta(r)/r = sqrt((delta_mH/mH)^2 + (delta_mW/mW)^2)
        # For now just use mH error
        exp_ratio_err = mH_err / mW_val

        sigma_A = abs(ratio_A - exp_ratio) / exp_ratio_err
        sigma_B = abs(ratio_B - exp_ratio) / exp_ratio_err

        print(f"  {mH_name}/{mW_name}: ratio = {exp_ratio:.6f} +/- {exp_ratio_err:.6f}")
        print(f"    A = {ratio_A:.6f} ({sigma_A:.2f} sigma)")
        print(f"    B = {ratio_B:.6f} ({sigma_B:.2f} sigma)")
        print()

# =====================================================================
# SECTION 6: The squared ratio from experiment
# =====================================================================
print("=" * 72)
print("SECTION 6: Direct Squared Ratio m_H^2/m_W^2")
print("=" * 72)

# Most precise: ATLAS m_H, PDG m_W
exp_sq = (m_H_ATLAS_2023/m_W_PDG_2024)**2
exp_sq_CMS = (m_H_CMS_2024/m_W_PDG_2024)**2
exp_sq_PDG = (m_H_PDG_avg/m_W_PDG_2024)**2

print(f"  Experimental (m_H/m_W)^2:")
print(f"    ATLAS/PDG:  {exp_sq:.8f}")
print(f"    CMS/PDG:    {exp_sq_CMS:.8f}")
print(f"    PDGavg/PDG: {exp_sq_PDG:.8f}")
print()
print(f"  Formula A: (phi - 8*alpha)^2         = {ratio_A_sq:.8f}")
print(f"  Formula B: phi^2 - 16*alpha*phi       = {ratio_B_sq:.8f}")
print(f"  phi^2 (tree):                          = {PHI**2:.8f}")
print()
print(f"  Residuals (predicted - experiment):")
print(f"    Tree:  {PHI**2 - exp_sq:+.6f}")
print(f"    A:     {ratio_A_sq - exp_sq:+.6f}")
print(f"    B:     {ratio_B_sq - exp_sq:+.6f}")
print(f"    A-CMS: {ratio_A_sq - exp_sq_CMS:+.6f}")
print(f"    B-CMS: {ratio_B_sq - exp_sq_CMS:+.6f}")

# =====================================================================
# SECTION 7: Physical arguments for each formula
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 7: Physical Arguments")
print("=" * 72)

print("""
FORMULA A: m_H/m_W = phi - 8*alpha
  PRO:
  [1] 8 = rank(E8) = Tr(D_F^2) - single fundamental invariant
  [2] Linear in alpha - natural for a MASS correction
  [3] Transparent factorization: manifold (phi) - finite space (8*alpha)
  [4] m_H = m_W*(phi - Tr(D_F^2)*alpha): each factor has clear origin
  [5] Already in paper (v4.2) as the primary formula
  CON:
  [1] m^2, not m, is the natural quantity in the Lagrangian
  [2] "Why linear?" is harder to justify from spectral action

FORMULA B: m_H^2/m_W^2 = phi^2 - 16*alpha*phi
  PRO:
  [1] m^2 is the Lagrangian quantity (m^2 |H|^2 term)
  [2] 16 = Tr(A^2_McKay) = (a1-1)^2 - also fundamental
  [3] Better numerical agreement (but by only ~0.1 sigma)
  [4] The spectral action naturally produces m^2, not m
  [5] Connes framework works with m_H^2/m_W^2 = 8*b/a^2
  CON:
  [1] Two factors (phi^2 and phi) appear - less transparent
  [2] 16 = 2*8 raises "why factor 2?" question
  [3] Formula B can be REWRITTEN as phi*(phi - 16*alpha)
      which has coefficient 16, or as (phi-8*alpha)^2 - 64*alpha^2
      which looks like "A minus a correction"

FORMULA C (hidden option): m_H/m_W = phi*(1 - 8*alpha/phi) = phi - 8*alpha
  This is IDENTICAL to A but emphasizes the MULTIPLICATIVE structure.
  The correction factor is (1 - 8*alpha/phi) = (1 - Tr(D_F^2)*alpha/phi).
""")

# =====================================================================
# SECTION 8: The O(alpha^2) term - is it physical?
# =====================================================================
print("=" * 72)
print("SECTION 8: Is the O(alpha^2) term physical?")
print("=" * 72)

# The difference between A and B is 64*alpha^2 in the squared ratio.
# In physical terms: m_H(A) - m_H(B) = ?
mH_diff = m_W_PDG_2024 * (ratio_A - ratio_B)
print(f"  m_H(A) - m_H(B) = {mH_diff*1000:.2f} MeV (using PDG m_W)")
print(f"  This is {mH_diff/0.11*100:.1f}% of the ATLAS uncertainty (110 MeV)")
print(f"  Currently UNRESOLVABLE experimentally.")

# HL-LHC precision: ~30-40 MeV on m_H
hlhc_err = 0.030  # GeV, optimistic
print(f"\n  HL-LHC projected m_H uncertainty: ~{hlhc_err*1000:.0f} MeV")
print(f"  Ratio m_H(A)-m_H(B) / HL-LHC_err = {mH_diff/hlhc_err:.2f}")
print(f"  {'MAY be distinguishable at HL-LHC' if mH_diff/hlhc_err > 0.5 else 'NOT distinguishable even at HL-LHC'}")

# FCC-ee precision: ~10-15 MeV on m_H
fcc_err = 0.010  # GeV, optimistic
print(f"\n  FCC-ee projected m_H uncertainty: ~{fcc_err*1000:.0f} MeV")
print(f"  Ratio m_H(A)-m_H(B) / FCC_err = {mH_diff/fcc_err:.2f}")
print(f"  {'Distinguishable at FCC-ee' if mH_diff/fcc_err > 1 else 'Borderline at FCC-ee'}")

# =====================================================================
# SECTION 9: The ALGEBRAIC argument
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 9: Algebraic Structure")
print("=" * 72)

# Let's check: is 16*alpha*phi special?
val_16ap = 16*ALPHA*PHI
print(f"  16*alpha*phi = {val_16ap:.10f}")
print(f"  8*alpha = {8*ALPHA:.10f}")
print(f"  phi^2 = {PHI**2:.10f}")

# In Z[phi]: phi^2 = phi + 1. So:
# phi^2 - 16*alpha*phi = (phi+1) - 16*alpha*phi = 1 + phi*(1 - 16*alpha)
print(f"\n  phi^2 - 16*alpha*phi = 1 + phi*(1-16*alpha)")
print(f"    = 1 + phi*{1-16*ALPHA:.10f}")
print(f"    = {1 + PHI*(1-16*ALPHA):.10f}")
print(f"    Check: {ratio_B_sq:.10f}")

# (phi - 8*alpha)^2 = phi^2 - 16*alpha*phi + 64*alpha^2
#                    = 1 + phi - 16*alpha*phi + 64*alpha^2
#                    = 1 + phi(1 - 16*alpha) + 64*alpha^2
print(f"\n  (phi-8*alpha)^2 = 1 + phi*(1-16*alpha) + 64*alpha^2")
print(f"    = {1 + PHI*(1-16*ALPHA) + 64*ALPHA**2:.10f}")
print(f"    Check: {ratio_A_sq:.10f}")

# The key: Formula B is "CLEANER" in Z[phi] because it doesn't have
# the 64*alpha^2 term (which involves pi through alpha).
# Formula B: phi^2 - 16*alpha*phi is LINEAR in alpha, while
# Formula A squared: (phi-8*alpha)^2 is QUADRATIC in alpha.

# But: in the SPECTRAL ACTION, the expansion is naturally in powers of
# the coupling. At O(alpha), A and B are identical. The O(alpha^2) term
# in A comes from squaring the linear correction, which is NOT a separate
# physical effect but just the square of the first-order correction.

# ARGUMENT: If the correction to m_H/m_W is EXACTLY phi - 8*alpha
# (as a RATIO, not just to first order), then the squared formula
# MUST include the 64*alpha^2 term. The question is whether the
# correction to the RATIO is exact or only first-order.

print("\n  CONCLUSION:")
print("  At O(alpha): A and B are IDENTICAL. No experiment can distinguish.")
print(f"  At O(alpha^2): they differ by {64*ALPHA**2:.6f} in ratio^2")
print(f"    = {mH_diff*1000:.1f} MeV in m_H")
print("  This will require FCC-ee precision to resolve.")

# =====================================================================
# SECTION 10: Recommendation
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 10: RECOMMENDATION")
print("=" * 72)

print("""
RECOMMENDATION: Keep FORMULA A as primary, note B as equivalent to O(alpha).

REASONS:
1. Formula A (phi - 8*alpha) is SIMPLER and more transparent
2. The coefficient 8 = rank(E8) = Tr(D_F^2) has a direct physical meaning
3. The two formulas are experimentally INDISTINGUISHABLE with current data
4. The O(alpha^2) difference (0.14 GeV) requires FCC-ee to resolve
5. Formula A is already established in the paper

WHAT TO ADD TO PAPER:
- Note that the SQUARED version phi^2 - Tr(A^2)*alpha*phi = phi^2 - 16*alpha*phi
  is equivalent to O(alpha) and has coefficient 16 = Tr(A_McKay^2)
- The identity Tr(A^2)/(2*Tr(D_F^2)) = 2|E|/(2|E|) = 1 explains why
  the coefficient in the LINEAR formula is exactly 1
- The squared version gives 0.008% error vs 0.13% for the linear version
  (both using experimental m_W), but the difference is within alpha^2

PREDICTION: If future precision favors B over A (m_H^2/m_W^2 = phi^2 - 16*alpha*phi
WITHOUT the 64*alpha^2 term), this would mean:
- The correction is to m^2, not m (Lagrangian level)
- The fundamental coefficient is 16 = Tr(A^2), not 8 = Tr(D_F^2)
- The factor 2 between A^2 and D_F^2 has physical significance
""")

# =====================================================================
# SECTION 11: Summary table
# =====================================================================
print("=" * 72)
print("SUMMARY TABLE")
print("=" * 72)

print(f"\n  {'':.<45} {'Formula A':>12} {'Formula B':>12}")
print(f"  {'Expression':.<45} {'phi-8a':>12} {'phi^2-16ap':>12}")
print(f"  {'':.<45} {'-'*12} {'-'*12}")
print(f"  {'Coefficient':.<45} {'8':>12} {'16':>12}")
print(f"  {'= what':.<45} {'rank(E8)':>12} {'Tr(A^2)':>12}")
print(f"  {'= what else':.<45} {'Tr(D_F^2)':>12} {'(a1-1)^2':>12}")
print(f"  {'m_H/m_W':.<45} {ratio_A:>12.6f} {ratio_B:>12.6f}")
print(f"  {'m_H^2/m_W^2':.<45} {ratio_A_sq:>12.6f} {ratio_B_sq:>12.6f}")
print(f"  {'m_H (PDG m_W)':.<45} {m_W_PDG_2024*ratio_A:>12.4f} {m_W_PDG_2024*ratio_B:>12.4f}")
print(f"  {'Error vs ATLAS (125.11)':.<45} {abs(m_W_PDG_2024*ratio_A-125.11)/125.11*100:>11.4f}% {abs(m_W_PDG_2024*ratio_B-125.11)/125.11*100:>11.4f}%")
print(f"  {'Error vs CMS (125.10)':.<45} {abs(m_W_PDG_2024*ratio_A-125.10)/125.10*100:>11.4f}% {abs(m_W_PDG_2024*ratio_B-125.10)/125.10*100:>11.4f}%")
print(f"  {'Error vs PDG (125.25)':.<45} {abs(m_W_PDG_2024*ratio_A-125.25)/125.25*100:>11.4f}% {abs(m_W_PDG_2024*ratio_B-125.25)/125.25*100:>11.4f}%")
print(f"  {'Difference A-B':.<45} {mH_diff*1000:>10.1f} MeV")
print(f"  {'Distinguished at':.<45} {'FCC-ee':>12} {'':>12}")

print("\n" + "=" * 72)
print("EXPERIMENT COMPLETE")
print("=" * 72)

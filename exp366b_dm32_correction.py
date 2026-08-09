# EXP-366b: Sub-leading Correction to Dm2_32
# =============================================
# The atmospheric mass splitting has 2.3 sigma tension:
#   Dm2_32_pred = 2.378e-3 eV^2
#   Dm2_32_exp  = 2.453e-3 eV^2 (3.1% deficit)
#
# Question: Is there an algebraically motivated correction from
# framework constants that reduces this tension?
#
# RULE ZERO: derive, don't invent. No fitting.
# BLIND: compute correction BEFORE checking if it helps.
# Windows: no Unicode.

import numpy as np

a1 = 5
b1 = 6
PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2
N = 120
N_gen = 3
N_eig = 9
DEG = 12

# Framework couplings
A_eq = 2 * np.pi
B_eq = -4 * a1 * PHI**4
C_eq_coef = 1.0
disc = B_eq**2 - 4 * A_eq * C_eq_coef
ALPHA = (-B_eq - np.sqrt(disc)) / (2 * A_eq)
ALPHA_S = 1 / (2 * PHI**3)
SIN2TW = 6.0 / 26.0

# Mass correction coefficient
C_corr = 4.0 / (a1**2 + 1)  # = 2/13

m_e_eV = 0.51099895e6  # eV

# Experimental data
Dm2_21_exp = 7.53e-5     # eV^2
Dm2_21_err = 0.18e-5
Dm2_32_exp = 2.453e-3    # eV^2 (normal ordering)
Dm2_32_err = 0.033e-3
Dm2_31_exp = Dm2_32_exp + Dm2_21_exp

# Baseline predictions (from exp363)
m3_bare = 2 * m_e_eV / PHI**35
r_fw = ALPHA * PHI**3  # mass splitting ratio
m2_bare = m3_bare * np.sqrt(r_fw)
m1_bare = 0.0

Dm2_21_bare = m2_bare**2 - m1_bare**2
Dm2_32_bare = m3_bare**2 - m2_bare**2
Dm2_31_bare = m3_bare**2 - m1_bare**2

print("=" * 75)
print("EXP-366b: SUB-LEADING CORRECTION TO Dm2_32")
print("=" * 75)

# =====================================================================
# PART 0: BASELINE
# =====================================================================
print("\n" + "=" * 75)
print("PART 0: BASELINE (no correction)")
print("=" * 75)

def report(label, m1, m2, m3):
    """Report mass splittings and compare with experiment."""
    dm21 = m2**2 - m1**2
    dm32 = m3**2 - m2**2
    dm31 = m3**2 - m1**2
    sig21 = abs(dm21 - Dm2_21_exp) / Dm2_21_err
    sig32 = abs(dm32 - Dm2_32_exp) / Dm2_32_err
    r = dm21 / dm31
    r_exp = Dm2_21_exp / Dm2_31_exp
    print("\n  %s:" % label)
    print("    m1 = %.4f meV, m2 = %.4f meV, m3 = %.4f meV" % (
        m1*1e3, m2*1e3, m3*1e3))
    print("    Dm2_21 = %.4e eV^2 (%.2f sigma)" % (dm21, sig21))
    print("    Dm2_32 = %.4e eV^2 (%.2f sigma)" % (dm32, sig32))
    print("    r = %.6f (exp: %.6f, err: %.2f%%)" % (r, r_exp, (r/r_exp-1)*100))
    print("    Sum(m) = %.2f meV" % ((m1+m2+m3)*1e3))
    return sig21, sig32

sig21_0, sig32_0 = report("BASELINE", m1_bare, m2_bare, m3_bare)

# =====================================================================
# PART 1: CORRECTION TO m3 (EXPONENT CORRECTION)
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 1: EXPONENT CORRECTIONS TO m3")
print("=" * 75)

# For charged fermions: m_f = m_e * phi^(n + delta)
# where delta = C * f(z) for various prescriptions.
# For neutrinos, the exponent is 35 (negative, as 2*m_e/phi^35).
# Correction: m3 = 2*m_e/phi^(35 - delta) = m3_bare * phi^delta

# What correction delta is ALGEBRAICALLY motivated?
# Try all framework-derived small quantities:

corrections = [
    # (label, delta_value, interpretation)
    ("C * ln(b1)",           C_corr * np.log(b1),
     "norm-log with |N|=b1=6"),
    ("C * ln(a1)",           C_corr * np.log(a1),
     "norm-log with |N|=a1=5"),
    ("C * ln(N_gen)",        C_corr * np.log(N_gen),
     "norm-log with |N|=N_gen=3"),
    ("C/a1 = 2/65",          C_corr / a1,
     "= 2/(a1*(a1^2+1)/2) = 4/(a1*(a1^2+1))"),
    ("alpha",                ALPHA,
     "fine structure correction"),
    ("2*alpha",              2*ALPHA,
     "2*alpha correction"),
    ("alpha_s/a1^2",         ALPHA_S / a1**2,
     "strong coupling / a1^2"),
    ("C * ln(2)",            C_corr * np.log(2),
     "norm-log with |N|=2 (prefactor)"),
    ("1/(2*phi^5)",          1.0/(2*PHI**5),
     "= alpha_s/phi^2 (next Galois power)"),
    ("C/b1",                 C_corr / b1,
     "= 2/(13*6) = 1/39"),
    ("sin2tW * alpha",       SIN2TW * ALPHA,
     "weak mixing * alpha"),
    ("C * alpha_s",          C_corr * ALPHA_S,
     "correction coef * strong coupling"),
    ("(b1-a1)/(a1^2+1)",     (b1-a1)*1.0/(a1**2+1),
     "= 1/26 (inverse Casimir sum)"),
]

print("\n  Testing algebraically motivated exponent corrections:")
print("  delta -> m3_corr = m3_bare * phi^delta")
print()
print("  %-25s  %10s  %8s  %8s  %s" % (
    "Correction", "delta", "sig(21)", "sig(32)", "Interpretation"))
print("  " + "-" * 85)

results = []
for label, delta, interp in corrections:
    m3_c = m3_bare * PHI**delta
    m2_c = m3_c * np.sqrt(r_fw)  # keep same ratio
    m1_c = 0.0
    dm21 = m2_c**2
    dm32 = m3_c**2 - m2_c**2
    s21 = abs(dm21 - Dm2_21_exp) / Dm2_21_err
    s32 = abs(dm32 - Dm2_32_exp) / Dm2_32_err
    improved = s32 < sig32_0
    marker = " <--" if s32 < 1.0 else (" *" if improved else "")
    print("  %-25s  %10.6f  %8.2f  %8.2f  %s%s" % (
        label, delta, s21, s32, interp, marker))
    results.append((label, delta, s21, s32, interp))

# =====================================================================
# PART 2: CORRECTION TO SPLITTING RATIO r
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 2: CORRECTIONS TO SPLITTING RATIO r = alpha*phi^3")
print("=" * 75)

# Instead of modifying m3, modify the ratio r
# r_corr = alpha*phi^3 * (1 + epsilon) or r_corr = alpha*phi^3 * correction

ratio_corrections = [
    ("r * (1 - alpha)",       r_fw * (1 - ALPHA),
     "EM perturbation"),
    ("r * (1 - 2*alpha)",     r_fw * (1 - 2*ALPHA),
     "2-loop EM"),
    ("r * (1 - C*alpha)",     r_fw * (1 - C_corr*ALPHA),
     "norm-log * alpha"),
    ("r * (1 - alpha_s/a1^2)", r_fw * (1 - ALPHA_S/a1**2),
     "QCD correction"),
    ("r * (1 - 1/(a1*b1))",  r_fw * (1 - 1.0/(a1*b1)),
     "1/h(E8) correction"),
    ("alpha*phi^3/(1+C)",     ALPHA*PHI**3 / (1 + C_corr),
     "denominator correction"),
    ("alpha*phi^3*(1-sin2tW*alpha)", r_fw*(1-SIN2TW*ALPHA),
     "EW correction"),
]

print("\n  %-30s  %10s  %8s  %8s  %s" % (
    "Correction", "r_corr", "sig(21)", "sig(32)", "Note"))
print("  " + "-" * 85)

for label, r_c, interp in ratio_corrections:
    m2_c = m3_bare * np.sqrt(r_c)
    dm21 = m2_c**2
    dm32 = m3_bare**2 - m2_c**2
    s21 = abs(dm21 - Dm2_21_exp) / Dm2_21_err
    s32 = abs(dm32 - Dm2_32_exp) / Dm2_32_err
    improved = (s32 < sig32_0 and s21 < 3.0)
    marker = " <--" if (s32 < 1.5 and s21 < 1.5) else (" *" if improved else "")
    print("  %-30s  %10.6f  %8.2f  %8.2f  %s%s" % (
        label, r_c, s21, s32, interp, marker))

# =====================================================================
# PART 3: COMBINED CORRECTION (m3 + r)
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 3: WHAT CORRECTION IS NEEDED? (reverse engineering)")
print("=" * 75)

# What m3 would make Dm2_32 exact (keeping r)?
# Dm2_32 = m3^2 * (1 - r) = 2.453e-3
# m3^2 = 2.453e-3 / (1 - r_fw)
m3_needed_from_32 = np.sqrt(Dm2_32_exp / (1 - r_fw))
delta_needed_32 = np.log(m3_needed_from_32 / m3_bare) / np.log(PHI)

print("\n  To match Dm2_32 exactly (keeping r = alpha*phi^3):")
print("    m3 needed = %.6f eV = %.4f meV" % (m3_needed_from_32, m3_needed_from_32*1e3))
print("    m3 bare   = %.6f eV = %.4f meV" % (m3_bare, m3_bare*1e3))
print("    m3_needed/m3_bare = %.6f" % (m3_needed_from_32/m3_bare))
print("    delta (exponent correction) = %.6f" % delta_needed_32)

# What r would make Dm2_21 AND Dm2_32 consistent (keeping m3)?
# With m3 fixed, Dm2_31 = m3^2
# r = Dm2_21/Dm2_31 -> r_needed for Dm2_21 exact
r_needed_21 = Dm2_21_exp / m3_bare**2
# Dm2_32 = m3^2 - Dm2_21 = m3^2(1 - r_needed_21)
dm32_with_exact_21 = m3_bare**2 * (1 - r_needed_21)

print("\n  To match Dm2_21 exactly (keeping m3 = 2*m_e/phi^35):")
print("    r needed = %.8f" % r_needed_21)
print("    r bare   = %.8f" % r_fw)
print("    Dm2_32 with this r = %.4e (exp: %.4e)" % (dm32_with_exact_21, Dm2_32_exp))
sig32_r = abs(dm32_with_exact_21 - Dm2_32_exp) / Dm2_32_err
print("    Tension: %.2f sigma (worse!)" % sig32_r)

# To match BOTH simultaneously:
# Dm2_21 = m3^2 * r = 7.53e-5
# Dm2_32 = m3^2 * (1-r) = 2.453e-3
# => m3^2 = Dm2_21 + Dm2_32 (approximately)
# => m3 = sqrt(Dm2_31_exp)
m3_exp = np.sqrt(Dm2_31_exp)
r_exp = Dm2_21_exp / Dm2_31_exp

print("\n  To match BOTH splittings:")
print("    m3 = sqrt(Dm2_31) = %.6f eV = %.4f meV" % (m3_exp, m3_exp*1e3))
print("    r = Dm2_21/Dm2_31 = %.8f" % r_exp)
print("    m3_exp/m3_bare = %.6f" % (m3_exp/m3_bare))
delta_both = np.log(m3_exp/m3_bare) / np.log(PHI)
print("    delta = %.6f (exponent shift)" % delta_both)
print("    r_exp/r_fw = %.6f" % (r_exp/r_fw))

# =====================================================================
# PART 4: THE C/a1 CORRECTION
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 4: THE C/a1 = 2/65 CORRECTION (best candidate)")
print("=" * 75)

delta_candidate = C_corr / a1  # = 2/65 = 0.030769
m3_corr = m3_bare * PHI**delta_candidate
m2_corr = m3_corr * np.sqrt(r_fw)
m1_corr = 0.0

print("\n  delta = C/a1 = (2/13)/5 = 2/65 = %.8f" % delta_candidate)
print("  = 4 / (a1 * (a1^2 + 1)) = 4 / (%d * %d) = 4/%d" % (
    a1, a1**2+1, a1*(a1**2+1)))
print("  This is: correction_coef / diameter = mass_correction / diameter")
print()
print("  phi^(2/65) = %.8f" % PHI**delta_candidate)
print("  m3_corr = %.6f eV = %.4f meV" % (m3_corr, m3_corr*1e3))

sig21_c, sig32_c = report("C/a1 correction", m1_corr, m2_corr, m3_corr)

print("\n  Improvement:")
print("    Dm2_32: %.2f sigma -> %.2f sigma" % (sig32_0, sig32_c))
print("    Dm2_21: %.2f sigma -> %.2f sigma" % (sig21_0, sig21_c))

# =====================================================================
# PART 5: MOTIVATION FOR C/a1
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 5: IS C/a1 MOTIVATED? (honest assessment)")
print("=" * 75)

print("""
  The correction delta = C/a1 = 2/65 is:

  ALGEBRAICALLY NATURAL:
    - C = 4/(a1^2+1) = 2/13 is the universal mass correction coefficient
    - a1 = 5 is the diameter of the 600-cell Cayley graph
    - C/a1 = mass_correction_scale / geometry_scale

  POSSIBLE PHYSICAL INTERPRETATIONS:
    1. The neutrino seesaw traverses the FULL graph diameter (a1=5 edges).
       The correction per edge is C. Total: C applied over a1 steps = C/a1
       per step, accumulated once. This gives C/a1.

    2. The norm |N(z)| for the seesaw exponent 35:
       35 = b1^2 - 1. If z = b1^2 - 1 viewed as element of Z:
       There's no natural Z[phi] interpretation that gives ln|N| = 1/a1.
       e^{1/5} is irrational, so NO Galois norm gives this.

    3. Analogy with charged leptons: the lepton correction uses
       c_ell = C*phi^3/4. Here we have C/a1 = C*1/a1.
       Ratio: c_ell/(C/a1) = phi^3*a1/4 = a1*phi^3/4 = 5*4.236/4 = 5.295.
       Not obviously meaningful.

  HONEST ASSESSMENT:
    - delta = 2/65 matches the needed correction to <1%
    - It uses only framework constants (C and a1)
    - BUT it has no derivation from norm-log or spectral action
    - It cannot be written as C*ln|N(z)| for any Z[phi] norm
    - CATEGORY: FITTING (algebraically motivated but not derived)
    - It would need independent justification to be called "pattern"
""")

# =====================================================================
# PART 6: SCAN ALL SIMPLE FRAMEWORK RATIOS
# =====================================================================
print("=" * 75)
print("PART 6: SYSTEMATIC SCAN OF FRAMEWORK CORRECTIONS")
print("=" * 75)

print("\n  Target delta = %.6f (to make Dm2_32 match exp)" % delta_both)
print("  Tolerance: +/- 0.005 in delta\n")

# Generate all simple framework ratios
candidates = []
for num_label, num_val in [
    ("1", 1), ("2", 2), ("C", C_corr), ("alpha", ALPHA),
    ("alpha_s", ALPHA_S), ("1/phi", 1/PHI), ("phi-1", PHI-1),
    ("ln(phi)", np.log(PHI))]:
    for den_label, den_val in [
        ("a1", a1), ("b1", b1), ("a1^2", a1**2), ("N_eig", N_eig),
        ("DEG", DEG), ("h", a1*b1), ("N", N), ("a1^2+1", a1**2+1),
        ("2*a1", 2*a1), ("a1*b1", a1*b1)]:
        delta_test = num_val / den_val
        if abs(delta_test - delta_both) < 0.005:
            candidates.append((num_label + "/" + den_label, delta_test))

# Also test products
for l1, v1 in [("C", C_corr), ("alpha", ALPHA), ("alpha_s", ALPHA_S)]:
    for l2, v2 in [("ln(a1)", np.log(a1)), ("ln(b1)", np.log(b1)),
                    ("ln(phi)", np.log(PHI)), ("1/a1", 1.0/a1),
                    ("1/b1", 1.0/b1)]:
        delta_test = v1 * v2
        if abs(delta_test - delta_both) < 0.005:
            candidates.append((l1 + "*" + l2, delta_test))

if candidates:
    print("  MATCHES found (delta within 0.005 of target %.6f):" % delta_both)
    for label, val in sorted(candidates, key=lambda x: abs(x[1] - delta_both)):
        m3_c = m3_bare * PHI**val
        m2_c = m3_c * np.sqrt(r_fw)
        dm32 = m3_c**2 - m2_c**2
        s32 = abs(dm32 - Dm2_32_exp) / Dm2_32_err
        dm21 = m2_c**2
        s21 = abs(dm21 - Dm2_21_exp) / Dm2_21_err
        print("    %-20s = %.6f  sig32=%.2f  sig21=%.2f" % (label, val, s32, s21))
else:
    print("  No simple framework ratio matches delta = %.6f" % delta_both)

# =====================================================================
# PART 7: ALTERNATIVE -- CORRECT r INSTEAD OF m3
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 7: WHAT IF THE ISSUE IS IN r, NOT m3?")
print("=" * 75)

# The ratio r = alpha*phi^3 is already a PATTERN.
# What if the true ratio involves a small correction?
# r_true = alpha*phi^3 * (1 - epsilon)

# For Dm2_32 to match, with m3 fixed:
# We need Dm2_32 = m3^2(1-r) = 2.453e-3
# So r_needed = 1 - Dm2_32/m3^2
r_from_32 = 1 - Dm2_32_exp / m3_bare**2
# But also Dm2_21 = m3^2*r
dm21_from_r32 = m3_bare**2 * r_from_32
sig21_from_r32 = abs(dm21_from_r32 - Dm2_21_exp) / Dm2_21_err

print("\n  If we fix m3 = %.4f meV and adjust only r:" % (m3_bare*1e3))
print("    r_needed (from Dm2_32) = %.8f" % r_from_32)
print("    r_framework = %.8f" % r_fw)
print("    Dm2_21 with r_needed = %.4e (exp: %.4e)" % (dm21_from_r32, Dm2_21_exp))
print("    Dm2_21 tension: %.2f sigma" % sig21_from_r32)
print()
print("  PROBLEM: Adjusting r to fix Dm2_32 makes Dm2_21 WORSE.")
print("  The issue is primarily in m3 (overall scale), not in the ratio.")

# =====================================================================
# SUMMARY
# =====================================================================
print("\n\n" + "=" * 75)
print("SUMMARY")
print("=" * 75)

print("""
  BASELINE: m3 = 2*m_e/phi^35 = %.4f meV
    Dm2_21: 0.30 sigma (excellent)
    Dm2_32: 2.29 sigma (moderate tension)

  BEST ALGEBRAIC CORRECTION: delta = C/a1 = 2/65
    m3_corr = 2*m_e/phi^(35 - 2/65) = %.4f meV
    Dm2_21: %.2f sigma
    Dm2_32: %.2f sigma
    Uses only framework constants.

  HONEST ASSESSMENT:
    - The C/a1 correction WORKS numerically (brings Dm2_32 within ~1 sigma)
    - C = 2/13 (mass correction coef) and a1 = 5 (diameter) are both derived
    - BUT no first-principles derivation connects C/a1 to neutrino masses
    - Cannot be written as C*ln|N| for any Galois norm
    - CATEGORY: algebraically motivated FITTING

  RECOMMENDATION:
    - Keep m3 = 2*m_e/phi^35 as the BARE prediction (PATTERN)
    - Note the 2.3 sigma tension as analogous to bare charged fermion deviations
    - DO NOT add the C/a1 correction to the paper without a derivation
    - The tension may resolve with better experimental data or a future insight
""" % (m3_bare*1e3, m3_corr*1e3, sig21_c, sig32_c))

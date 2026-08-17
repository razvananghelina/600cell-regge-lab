# EXP-366: Algebraic Structure of alpha(mZ)/alpha(0)
# =====================================================
# QUESTION: Is the ratio 16/15 = (3*a1+1)/(3*a1) DERIVABLE from the
# 1-loop QED sum with pure framework masses?
#
# Previous work (exp361): 1-loop gives 1.0628, exp = 1.0631, 16/15 = 1.0667.
# Light quarks (u,d,s) are nonperturbative. Can we isolate the perturbative
# part and check if it has an algebraic expression in framework constants?
#
# RULE ZERO: derive, don't invent. Categorize honestly.
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

# Framework alpha from quadratic equation
A_eq = 2 * np.pi
B_eq = -4 * a1 * PHI**4
C_eq = 1.0
disc = B_eq**2 - 4 * A_eq * C_eq
ALPHA = (-B_eq - np.sqrt(disc)) / (2 * A_eq)
ALPHA_S = 1 / (2 * PHI**3)

m_e = 0.51099895  # MeV
m_Z = 91187.6     # MeV (PDG)

# Correction coefficient
C_corr = 4.0 / (a1**2 + 1)  # = 2/13

print("=" * 75)
print("EXP-366: ALGEBRAIC STRUCTURE OF alpha(mZ)/alpha(0)")
print("=" * 75)

# =====================================================================
# PART 1: FRAMEWORK FERMION MASSES (BARE)
# =====================================================================
print("\n" + "=" * 75)
print("PART 1: FRAMEWORK FERMION MASSES")
print("=" * 75)

# (name, a, b, Nc, Q, Q^2)
fermions = [
    ("e",   0,  0, 1, -1.0, 1.0),
    ("mu",  1,  1, 1, -1.0, 1.0),
    ("tau", 1,  2, 1, -1.0, 1.0),
    ("u",   3, -2, 3,  2.0/3, 4.0/9),
    ("d",   1,  0, 3, -1.0/3, 1.0/9),
    ("c",   2,  1, 3,  2.0/3, 4.0/9),
    ("s",   1,  1, 3, -1.0/3, 1.0/9),
    ("b",  -1,  4, 3, -1.0/3, 1.0/9),
    ("t",   4,  1, 3,  2.0/3, 4.0/9),
]

print("\n  %-5s  n   m_bare/MeV    Nc  Q^2    Nc*Q^2" % "Name")
print("  " + "-" * 55)

S_Q_total = 0  # total sum Nc*Q^2
sum_perturbative = 0  # sum Nc*Q^2 * ln(mZ^2/mf^2) for active fermions

for name, a, b, Nc, Q, Q2 in fermions:
    n = a1 * a + b1 * b
    m_bare = m_e * PHI**n
    NcQ2 = Nc * Q2
    S_Q_total += NcQ2

    active = m_bare < m_Z
    if active:
        ln_term = np.log(m_Z**2 / m_bare**2)
        contribution = NcQ2 * ln_term
        sum_perturbative += contribution
    else:
        ln_term = 0
        contribution = 0

    marker = "" if active else " (INACTIVE, m > mZ)"
    print("  %-5s %3d  %12.4f  %d  %4.2f   %5.2f%s" % (
        name, n, m_bare, Nc, Q2, NcQ2, marker))

print("\n  S_Q = sum(Nc*Q^2) for active = %.4f" % S_Q_total)
# Note: top quark is inactive (m_t > m_Z in bare framework)

# 1-loop running
Delta_alpha = ALPHA / (3 * np.pi) * sum_perturbative
ratio_1loop = 1.0 / (1.0 - Delta_alpha)
alpha_inv_mZ = 1.0 / (ALPHA * ratio_1loop)

print("\n--- 1-loop QED Running (BARE masses) ---")
print("  Delta_alpha = alpha/(3*pi) * sum = %.6f" % Delta_alpha)
print("  alpha(mZ)/alpha(0) = 1/(1-Da) = %.6f" % ratio_1loop)
print("  1/alpha(mZ) = %.4f" % alpha_inv_mZ)
print("  Experimental: 1/alpha(mZ) = 127.944")
print("  16/15 = %.6f" % (16.0/15))
print("  (3*a1+1)/(3*a1) = %d/%d = %.6f" % (3*a1+1, 3*a1, (3*a1+1.0)/(3*a1)))

# =====================================================================
# PART 2: ALGEBRAIC DECOMPOSITION OF THE SUM
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 2: ALGEBRAIC DECOMPOSITION")
print("=" * 75)

# For framework masses: m_f = m_e * phi^n_f, m_Z = m_e * phi^25 * 16/15
# ln(mZ/mf) = ln(phi^(25-n_f) * 16/15) = (25-n_f)*ln(phi) + ln(16/15)
# ln(mZ^2/mf^2) = 2*(25-n_f)*ln(phi) + 2*ln(16/15)

# So the 1-loop sum decomposes into:
# sum = 2*ln(phi) * sum_f[Nc*Q^2*(25-n_f)] + 2*ln(16/15) * sum_f[Nc*Q^2]

ln_phi = np.log(PHI)
ln_ratio = np.log(16.0/15)

S_weighted = 0  # sum Nc*Q^2*(25-n_f) for active fermions
S_charge = 0    # sum Nc*Q^2 for active fermions

print("\n  Decomposition: sum = 2*ln(phi)*S_w + 2*ln(16/15)*S_Q")
print("\n  %-5s   n   25-n   Nc*Q^2  Nc*Q^2*(25-n)" % "Name")
print("  " + "-" * 50)

for name, a, b, Nc, Q, Q2 in fermions:
    n = a1 * a + b1 * b
    m_bare = m_e * PHI**n
    NcQ2 = Nc * Q2

    if m_bare < m_Z:
        diff = 25 - n
        weighted = NcQ2 * diff
        S_weighted += weighted
        S_charge += NcQ2
        print("  %-5s  %3d   %3d   %5.2f   %8.2f" % (
            name, n, diff, NcQ2, weighted))

print("\n  S_weighted = sum Nc*Q^2*(25-n_f) = %.4f" % S_weighted)
print("  S_charge   = sum Nc*Q^2           = %.4f" % S_charge)

# Framework exponents
exponents = {}
for name, a, b, Nc, Q, Q2 in fermions:
    n = a1 * a + b1 * b
    exponents[name] = n

# Print the algebraic structure
print("\n--- Algebraic Structure ---")
print("  Lepton exponents: e=%d, mu=%d, tau=%d" % (
    exponents['e'], exponents['mu'], exponents['tau']))
print("  Up-type exponents: u=%d, c=%d, t=%d" % (
    exponents['u'], exponents['c'], exponents['t']))
print("  Down-type exponents: d=%d, s=%d, b=%d" % (
    exponents['d'], exponents['s'], exponents['b']))

# Compute S_weighted analytically
# Leptons: Nc=1, Q^2=1, n_e=0, n_mu=11, n_tau=17
# (25-0) + (25-11) + (25-17) = 25 + 14 + 8 = 47
S_lep = (25 - 0) + (25 - 11) + (25 - 17)

# Up quarks (Nc=3, Q^2=4/9): u(n=3), c(n=16), t(n=26)
# t is INACTIVE if bare mass < m_Z. Bare m_t = m_e*phi^26 = 138.7 GeV > m_Z.
# So only u and c: 3*(4/9)*[(25-3) + (25-16)] = (4/3)*(22+9) = (4/3)*31
S_up = 3 * (4.0/9) * ((25 - exponents['u']) + (25 - exponents['c']))
# Note: check if t active
m_t_bare = m_e * PHI**26
t_active = m_t_bare < m_Z
print("\n  m_t(bare) = %.2f MeV, m_Z = %.2f MeV" % (m_t_bare, m_Z))
print("  Top quark active? %s" % t_active)

if t_active:
    S_up += 3 * (4.0/9) * (25 - exponents['t'])

# Down quarks (Nc=3, Q^2=1/9): d(n=5), s(n=11), b(n=18)
S_down = 3 * (1.0/9) * ((25 - exponents['d']) + (25 - exponents['s']) + (25 - exponents['b']))

print("\n  S_lep (exact integers) = %d" % S_lep)
print("  S_up  = 3*(4/9)*sum = %.4f" % S_up)
print("  S_down = 3*(1/9)*sum = %.4f" % S_down)
print("  S_weighted = S_lep + S_up + S_down = %.4f" % (S_lep + S_up + S_down))
print("  Check: %.4f (from loop above)" % S_weighted)

# More explicitly:
# Up: (25-3) + (25-16) = 22 + 9 = 31. S_up = (4/3)*31 = 124/3
# Down: (25-5) + (25-11) + (25-18) = 20 + 14 + 7 = 41. S_down = (1/3)*41 = 41/3
# Total: 47 + 124/3 + 41/3 = 47 + 165/3 = 47 + 55 = 102

up_sum = (25 - exponents['u']) + (25 - exponents['c'])
if t_active:
    up_sum += (25 - exponents['t'])
down_sum = (25 - exponents['d']) + (25 - exponents['s']) + (25 - exponents['b'])

print("\n  Lepton (25-n) sum: %d" % S_lep)
print("  Up quark (25-n) sum: %d  (active: u, c%s)" % (
    up_sum, ", t" if t_active else ""))
print("  Down quark (25-n) sum: %d" % down_sum)
print("  Total S_weighted = %d + (4/3)*%d + (1/3)*%d = %.4f" % (
    S_lep, up_sum, down_sum, S_lep + (4.0/3)*up_sum + (1.0/3)*down_sum))

# S_charge for active fermions
# Leptons: 3*1 = 3
# Up (u,c): 2*3*(4/9) = 8/3
# Down (d,s,b): 3*3*(1/9) = 1
S_charge_exact = 3 + 2*(4.0/3) + 3*(1.0/3)
print("\n  S_charge = 3 + 8/3 + 1 = %.4f = %d/3" % (S_charge_exact, int(round(S_charge_exact*3))))

# =====================================================================
# PART 3: IS 16/15 ALGEBRAICALLY EXACT?
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 3: DOES THE 1-LOOP SUM GIVE 16/15?")
print("=" * 75)

# alpha(mZ)/alpha(0) = 1/(1 - Da)
# Da = alpha/(3*pi) * [2*ln(phi)*S_w + 2*ln(16/15)*S_Q]
# For 16/15 to be exact: Da = 1 - 15/16 = 1/16

Da_needed = 1.0 - 15.0/16  # = 1/16 = 0.0625
Da_computed = ALPHA / (3 * np.pi) * (2 * ln_phi * S_weighted + 2 * ln_ratio * S_charge)

print("\n  Delta_alpha needed for ratio = 16/15: %.6f (= 1/16)" % Da_needed)
print("  Delta_alpha from 1-loop sum:         %.6f" % Da_computed)
print("  Difference:                          %.6f" % (Da_computed - Da_needed))
print("  Ratio Da_computed/Da_needed:          %.6f" % (Da_computed/Da_needed))

# What sum_perturbative would give Da = 1/16?
sum_needed = Da_needed * 3 * np.pi / ALPHA
print("\n  For Da = 1/16, need sum_pert = %.4f" % sum_needed)
print("  Actual sum_pert = %.4f" % sum_perturbative)
print("  Deficit: %.4f (%.2f%%)" % (sum_needed - sum_perturbative,
    (sum_needed - sum_perturbative)/sum_needed * 100))

# =====================================================================
# PART 4: ALGEBRAIC IDENTITY SEARCH
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 4: SEARCHING FOR ALGEBRAIC IDENTITIES")
print("=" * 75)

# The 1-loop formula without the ln(16/15) correction (pure phi-based):
# Da_phi = alpha/(3*pi) * 2*ln(phi) * S_weighted
Da_phi = ALPHA / (3 * np.pi) * 2 * ln_phi * S_weighted
ratio_phi = 1.0 / (1.0 - Da_phi)

print("\n  Da (phi-only, no 16/15 correction) = %.6f" % Da_phi)
print("  Ratio (phi-only) = %.6f" % ratio_phi)

# Check if S_weighted has a clean value
# S_weighted = 47 + 124/3 + 41/3 = 47 + 55 = 102 (if t inactive)
# Exact: 141/3 + 124/3 + 41/3 = 306/3 = 102?
# Wait: 47 = 141/3? No, 47*3 = 141.
# Total = 47 + 124/3 + 41/3 = (141 + 124 + 41)/3 = 306/3 = 102
S_w_rational = S_lep + (4.0/3)*up_sum + (1.0/3)*down_sum
print("\n  S_weighted = %d + (4/3)*%d + (1/3)*%d" % (S_lep, up_sum, down_sum))
S_w_num = 3*S_lep + 4*up_sum + down_sum
print("  = (%d + %d + %d)/3 = %d/3" % (3*S_lep, 4*up_sum, down_sum, S_w_num))
if S_w_num % 3 == 0:
    print("  = %d (exact integer!)" % (S_w_num // 3))

# The condition for ratio = 16/15:
# 1/(1 - alpha*2*ln(phi)*S_w/(3*pi)) = 16/15
# 1 - 15/16 = alpha*2*ln(phi)*S_w/(3*pi)
# S_w = 3*pi / (16*2*alpha*ln(phi))
S_w_needed = 3 * np.pi / (16 * 2 * ALPHA * ln_phi)
print("\n  For EXACT ratio=16/15 (ignoring ln(16/15) term):")
print("  Need S_w = 3*pi / (32*alpha*ln(phi)) = %.4f" % S_w_needed)
print("  Have S_w = %d" % int(round(S_w_rational)))
print("  Ratio: %.6f" % (S_w_rational / S_w_needed))

# Check: is 3*pi/(32*alpha*ln(phi)) close to an integer or simple fraction?
print("\n  3*pi/(32*alpha*ln(phi)) = %.6f" % S_w_needed)
print("  Nearby integers: %d, %d" % (int(S_w_needed), int(S_w_needed)+1))

# Try different denominators
for den in range(1, 20):
    val = S_w_needed * den
    nearest = round(val)
    if abs(val - nearest) < 0.02:
        print("  %d * S_needed = %.4f ~ %d/%d" % (den, val, nearest, den))

# =====================================================================
# PART 5: CONTRIBUTION DECOMPOSITION BY SECTOR
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 5: SECTOR DECOMPOSITION")
print("=" * 75)

# Compute individual sector contributions to Da
sectors = {
    'leptons':    [],
    'up_quarks':  [],
    'down_quarks': [],
}

for name, a, b, Nc, Q, Q2 in fermions:
    n = a1 * a + b1 * b
    m_bare = m_e * PHI**n
    if m_bare >= m_Z:
        continue
    ln_term = np.log(m_Z**2 / m_bare**2)
    Da_f = ALPHA / (3 * np.pi) * Nc * Q2 * ln_term

    if Nc == 1:
        sectors['leptons'].append((name, n, Da_f))
    elif Q2 > 0.2:  # up-type
        sectors['up_quarks'].append((name, n, Da_f))
    else:
        sectors['down_quarks'].append((name, n, Da_f))

print("\n  %-10s %-5s  %4s  Delta_alpha_f  Cumulative" % ("Sector", "Name", "n"))
print("  " + "-" * 55)

cumulative = 0
for sector_name in ['leptons', 'up_quarks', 'down_quarks']:
    sector_total = 0
    for name, n, Da_f in sectors[sector_name]:
        cumulative += Da_f
        sector_total += Da_f
        print("  %-10s %-5s  %4d  %.6f    %.6f" % (
            sector_name, name, n, Da_f, cumulative))
    print("  %-10s TOTAL        %.6f" % (sector_name, sector_total))
    print()

print("  Total Delta_alpha = %.6f" % cumulative)
print("  Ratio = 1/(1-Da) = %.6f" % (1.0/(1.0 - cumulative)))
print("  16/15 = %.6f" % (16.0/15))
print("  Experimental = 1.0631")

# =====================================================================
# PART 6: THE TOP QUARK THRESHOLD
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 6: TOP QUARK THRESHOLD ANALYSIS")
print("=" * 75)

# If we include top (even though m_t_bare > m_Z, it contributes above threshold)
m_t_bare = m_e * PHI**26
print("\n  m_t(bare) = %.2f GeV (framework)" % (m_t_bare/1000))
print("  m_t(PDG)  = 172.76 GeV")
print("  m_Z       = 91.19 GeV")

# With top quark:
sum_with_top = sum_perturbative
if not t_active:
    ln_top = np.log(m_Z**2 / m_t_bare**2)  # negative!
    Da_top = ALPHA / (3 * np.pi) * 3 * (4.0/9) * ln_top
    sum_with_top += 3 * (4.0/9) * ln_top
    print("\n  Top contribution (above threshold, negative): %.6f" % Da_top)
    Da_with_top = ALPHA / (3 * np.pi) * sum_with_top
    ratio_with_top = 1.0 / (1.0 - Da_with_top)
    print("  Da (including top) = %.6f" % Da_with_top)
    print("  Ratio (with top) = %.6f" % ratio_with_top)
    print("  Still 16/15? Diff = %.6f" % (ratio_with_top - 16.0/15))

# =====================================================================
# PART 7: CORRECTED MASSES
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 7: WITH NORM-LOG CORRECTED MASSES")
print("=" * 75)

# Norm-log correction: delta = C * ln|N(z)| for quarks, C*|z'|^{3/4} for leptons
# Let's compute with the best-fit corrected masses from the framework

# PDG reference masses (what the running is measured against)
pdg_masses = {
    'e':   0.51099895,
    'mu':  105.6584,
    'tau': 1776.86,
    'u':   2.16,
    'c':   1270.0,
    't':   172760.0,
    'd':   4.67,
    's':   93.4,
    'b':   4180.0,
}

print("\n  Using PDG masses (as sanity check for 1-loop):")
sum_pdg = 0
for name, a, b, Nc, Q, Q2 in fermions:
    m_pdg = pdg_masses[name]
    if m_pdg >= m_Z:
        continue
    ln_term = np.log(m_Z**2 / m_pdg**2)
    sum_pdg += Nc * Q2 * ln_term

Da_pdg = ALPHA / (3 * np.pi) * sum_pdg
ratio_pdg = 1.0 / (1.0 - Da_pdg)
print("  Da (PDG masses) = %.6f" % Da_pdg)
print("  Ratio (PDG) = %.6f" % ratio_pdg)
print("  Experimental: 1.0631 (includes hadronic VP)")
print("  Note: 1-loop with PDG masses UNDERESTIMATES (missing hadronic VP)")

# =====================================================================
# PART 8: FRAMEWORK EXPRESSION FOR 16/15
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 8: FRAMEWORK EXPRESSION FOR 16/15")
print("=" * 75)

# 16/15 = (3*a1+1)/(3*a1) = 1 + 1/(3*a1) = 1 + 1/15
# Also: 16 = dim(WHITE) = fermions per generation
#        15 = a1*N_gen = 5*3 = sum(YM beta coefficients 8+5+2)
# So: alpha(mZ)/alpha(0) = dim(WHITE)/(a1*N_gen) ???

print("\n  16/15 = (3*a1+1)/(3*a1) = %d/%d" % (3*a1+1, 3*a1))
print("  16 = dim(WHITE) = (a1-1)^2 = fermions per SM generation")
print("  15 = a1 * N_gen = sum of YM beta coefs (8+5+2)")
print("  Correction = 1/(3*a1) = 1/%d" % (3*a1))
print()
print("  Physical interpretation:")
print("  - Numerator 16: number of chiral fermions contributing to running")
print("  - Denominator 15: gauge beta function sum (topological)")
print("  - 16/15 > 1: fermions dominate over gauge sector")
print()

# Check: 1/(3*a1) = 1/15 = 0.0667
# Da = 1/16 = 0.0625. Close but not equal.
# 1 - 15/16 = 1/16. If ratio = 16/15, Da = 1 - 15/16 = 1/16.
print("  If ratio = 16/15 exact:")
print("    Da = 1 - 15/16 = 1/16 = %.6f" % (1.0/16))
print("    Da = alpha * 2*ln(phi)*S_w / (3*pi)")
print("    => S_w = 3*pi/(32*alpha*ln(phi)) = %.4f" % S_w_needed)
print("    Actual S_w = %d" % int(round(S_w_rational)))
print()
ratio_actual = S_w_rational / S_w_needed
print("  S_w(actual)/S_w(needed) = %.6f" % ratio_actual)
print("  This ratio must be ~1 for 16/15 to be algebraically exact.")
print("  Deviation: %.2f%%" % ((ratio_actual - 1) * 100))

# =====================================================================
# SUMMARY
# =====================================================================
print("\n\n" + "=" * 75)
print("SUMMARY")
print("=" * 75)

print("""
  1-LOOP RESULTS:
    Framework bare masses:   alpha(mZ)/alpha(0) = %.6f
    PDG masses:              alpha(mZ)/alpha(0) = %.6f
    16/15 = (3a1+1)/(3a1):                      = %.6f
    Experimental:                                = 1.0631

  ALGEBRAIC STRUCTURE:
    S_weighted = %d (exact integer from framework exponents)
    S_charge = %d/3
    Sum decomposes as: leptons(%d) + up(4/3*%d) + down(1/3*%d)

  VERDICT:
    - The 1-loop sum S_w = %d is an exact integer (framework-derived).
    - But S_w = %d != 3*pi/(32*alpha*ln(phi)) = %.1f.
    - The ratio involves TRANSCENDENTAL quantities (pi, ln(phi)).
    - Therefore 16/15 is NOT algebraically exact from 1-loop QED.
    - It remains a PATTERN: (3a1+1)/(3a1) with framework meaning
      but no proof from the running sum.

  CATEGORY: PATTERN (16/15 motivated, not derived from 1-loop)
""" % (ratio_1loop, ratio_pdg, 16.0/15,
       int(round(S_w_rational)), int(round(S_charge*3)),
       S_lep, up_sum, down_sum,
       int(round(S_w_rational)), int(round(S_w_rational)), S_w_needed))

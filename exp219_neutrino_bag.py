# EXP-219: Neutrino Mass from alpha^2 and Soliton Bag Geometry
# =============================================================
# Proposal from masa_neutrinilor.md (user's idea):
#   m_nu = m_e * alpha^2 * phi^(-(13+k))
# where:
#   alpha^2 = Majorana reflection (order-2 process)
#   13 = soliton bag vertices (derived EXP-178)
#   k = PMNS exponents {0, 3, 4}
#
# RULE ZERO: verify ALL claims numerically. No rounding tricks.

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 1/137.035999084
ALPHA_S = 0.1179
SIN2TW = 6.0/26.0

m_e = 0.51099895e6  # eV

# Neutrino data (PDG 2024, normal hierarchy)
Dm2_21 = 7.53e-5   # eV^2
Dm2_31 = 2.455e-3   # eV^2
m3_exp = np.sqrt(Dm2_31)  # 0.04955 eV
m2_exp = np.sqrt(Dm2_21)  # 0.008678 eV

# Soliton bag
N_bag = 13  # vertices in soliton bag (DERIVED, EXP-178)

print("="*70)
print("EXP-219: NEUTRINO MASS FROM alpha^2 + SOLITON BAG")
print("="*70)

# ===================================================================
print("\n--- BASE FORMULA: m_nu = m_e * alpha^2 * phi^(-(13+k)) ---")
# ===================================================================
base = m_e * ALPHA**2
print(f"m_e * alpha^2 = {base:.4f} eV")
print(f"phi^13 = {PHI**13:.3f}")
print()

# User's claim: k = {0, 3, 4} = "PMNS exponents"
# BUT: actual PMNS exponents from EXP-085 are {0, 1, 4}
# k=3 is the CKM theta_12 (Cabibbo) exponent, NOT PMNS!

print("--- User's proposed k = {0, 3, 4} ---")
k_user = [0, 3, 4]
names = ["nu_3 (heaviest)", "nu_2 (middle)", "nu_1 (lightest)"]
m_exp_list = [m3_exp, m2_exp, 0.0]  # m1 ~ 0 in NH
m_pred_user = []

for i, (name, k, m_exp) in enumerate(zip(names, k_user, m_exp_list)):
    n = N_bag + k
    m_pred = base / PHI**n
    m_pred_user.append(m_pred)
    err_str = f"{abs(m_pred-m_exp)/m_exp*100:.1f}%" if m_exp > 0 else "N/A"
    print(f"  {name}: k={k}, n={n}, m = {m_pred:.5f} eV (exp: {m_exp:.5f}, err: {err_str})")

m_sum_user = sum(m_pred_user)
print(f"\n  Sum = {m_sum_user:.4f} eV (Planck limit: < 0.12 eV) {'OK' if m_sum_user < 0.12 else 'VIOLATED!'}")

# ===================================================================
print("\n--- CORRECT PMNS exponents k = {0, 1, 4} ---")
# ===================================================================
k_pmns = [0, 1, 4]
m_pred_pmns = []

for i, (name, k, m_exp) in enumerate(zip(names, k_pmns, m_exp_list)):
    n = N_bag + k
    m_pred = base / PHI**n
    m_pred_pmns.append(m_pred)
    err_str = f"{abs(m_pred-m_exp)/m_exp*100:.1f}%" if m_exp > 0 else "N/A"
    print(f"  {name}: k={k}, n={n}, m = {m_pred:.5f} eV (exp: {m_exp:.5f}, err: {err_str})")

m_sum_pmns = sum(m_pred_pmns)
print(f"\n  Sum = {m_sum_pmns:.4f} eV {'OK' if m_sum_pmns < 0.12 else 'VIOLATED!'}")

# ===================================================================
print("\n\n--- CRITICAL ISSUE: k=3 is NOT a PMNS exponent ---")
# ===================================================================
print("""
PMNS exponents (EXP-085, EXP-206):
  theta_23: n = 0 (atmospheric, maximal mixing)
  theta_12: n = 1 (solar, golden ratio mixing)
  theta_13: n = 4 (reactor, space dimension)

CKM exponents (EXP-191):
  theta_12: n = 3 (Cabibbo)
  theta_23: n = 7
  theta_13: n = 12

k=3 is CKM theta_12, NOT PMNS!
Using {0, 3, 4} mixes CKM and PMNS exponents.
""")

# ===================================================================
print("--- MASS-ANGLE CONNECTION (user's claim) ---")
# ===================================================================
print("If k = {0, 3, 4}, then n = {13, 16, 17}:")
print("  n(nu_1) - n(nu_2) = 17 - 16 = 1 --> PMNS theta_12 exponent")
print("  n(nu_1) - n(nu_3) = 17 - 13 = 4 --> PMNS theta_13 exponent")
print("  n(nu_2) - n(nu_3) = 16 - 13 = 3 --> CKM theta_12 exponent (!)")
print()
print("This is interesting but INCONSISTENT:")
print("  The input k-values are {0, 3, 4} (mixing CKM and PMNS)")
print("  The output differences are {1, 4, 3} (also mixing CKM and PMNS)")
print("  It works precisely BECAUSE 3 = CKM_theta12 = 3*PMNS_theta12")
print("  This is just the factor-3 rule in disguise, not a new result.")

# ===================================================================
print("\n\n--- COMPARISON: All neutrino mass proposals ---")
# ===================================================================
print(f"{'Formula':<45s} {'m3 (eV)':>10s} {'err_m3':>8s} {'m2 (eV)':>10s} {'err_m2':>8s} {'Sum':>8s}")
print("-"*90)

# 1. User's proposal: m_e * alpha^2 * phi^(-(13+k))
m3_a = base/PHI**13; m2_a = base/PHI**16; m1_a = base/PHI**17
print(f"{'m_e*a^2*phi^-(13+k), k={0,3,4}':<45s} {m3_a:>10.5f} {abs(m3_a-m3_exp)/m3_exp*100:>7.1f}% {m2_a:>10.5f} {abs(m2_a-m2_exp)/m2_exp*100:>7.1f}% {m3_a+m2_a+m1_a:>7.3f}")

# 2. EXP-218: m_e * (alpha/phi)^3
m3_b = m_e*(ALPHA/PHI)**3; m2_b = m3_b*2*np.sqrt(ALPHA)
print(f"{'m_e*(a/phi)^3, ratio=2*sqrt(a)':<45s} {m3_b:>10.5f} {abs(m3_b-m3_exp)/m3_exp*100:>7.1f}% {m2_b:>10.5f} {abs(m2_b-m2_exp)/m2_exp*100:>7.1f}% {m3_b+m2_b:>7.3f}")

# 3. Previous masa_neutrinilor: m_e * phi^(-34)
m3_c = m_e/PHI**34
print(f"{'m_e*phi^-34 (old n=h+d)':<45s} {m3_c:>10.5f} {abs(m3_c-m3_exp)/m3_exp*100:>7.1f}% {'---':>10s} {'---':>8s} {'---':>8s}")

# 4. With correct PMNS k={0,1,4}
m3_d = base/PHI**13; m2_d = base/PHI**14; m1_d = base/PHI**17
print(f"{'m_e*a^2*phi^-(13+k), k={0,1,4} PMNS':<45s} {m3_d:>10.5f} {abs(m3_d-m3_exp)/m3_exp*100:>7.1f}% {m2_d:>10.5f} {abs(m2_d-m2_exp)/m2_exp*100:>7.1f}% {m3_d+m2_d+m1_d:>7.3f}")

# ===================================================================
print("\n\n--- WHAT'S GENUINELY GOOD IN THIS PROPOSAL ---")
# ===================================================================
print("""
POSITIV:
1. |Bag| = 13 as base exponent: TOPOLOGICAL, DERIVED (EXP-178)
   Much better than 30+4=34 (numerology)
2. alpha^2 suppression: Has physical motivation
3. Sum < 0.12 eV: Cosmologically consistent
4. Mass-angle differences {1, 4}: Genuinely interesting pattern

NEGATIV:
1. k = {0, 3, 4}: NOT PMNS exponents! k=3 is CKM theta_12
   With correct PMNS {0, 1, 4}: m2 error = 272%
2. m2 error = 42% even with user's {0, 3, 4}: Too large
3. alpha^2 = "Majorana reflection": Not derived.
   Majorana mass doesn't require EM interaction.
   Why alpha^2 and not G_F or alpha_s?
4. The mass-angle connection is the factor-3 rule repackaged
""")

# ===================================================================
print("\n--- CAN WE SALVAGE THE GOOD PARTS? ---")
# ===================================================================
print("The core insight: m_nu ~ m_e * (coupling)^p * phi^(-N_bag)")
print("Let's find what coupling^p * phi^(-13) gives m3:\n")

# m3 = m_e * X * phi^(-13)
# X = m3 * phi^13 / m_e
X_needed = m3_exp * PHI**13 / m_e
print(f"Need: X = m3*phi^13/m_e = {X_needed:.6e}")
print(f"Compare:")
print(f"  alpha^2         = {ALPHA**2:.6e} (ratio: {X_needed/ALPHA**2:.4f})")
print(f"  alpha^2*phi     = {ALPHA**2*PHI:.6e} (ratio: {X_needed/(ALPHA**2*PHI):.4f})")
print(f"  alpha*alpha_s   = {ALPHA*ALPHA_S:.6e} (ratio: {X_needed/(ALPHA*ALPHA_S):.4f})")
print(f"  alpha_s^2/(4pi) = {ALPHA_S**2/(4*np.pi):.6e} (ratio: {X_needed/(ALPHA_S**2/(4*np.pi)):.4f})")
print(f"  alpha^2*(1+sin2tW/4) = {ALPHA**2*(1+SIN2TW/4):.6e} (ratio: {X_needed/(ALPHA**2*(1+SIN2TW/4)):.4f})")

# Best: alpha^2 needs a factor of ~0.95 to match.
# Or equivalently, phi^(-13) is slightly too small.

# What if we use phi^(-12.88) instead of phi^(-13)?
n_best = np.log(base/m3_exp) / np.log(PHI)
print(f"\nBest continuous exponent: n = {n_best:.3f}")
print(f"  Deficit from 13: {13 - n_best:.3f}")
print(f"  This is NOT a clean number.")

# What about m_e * alpha^2 * phi^(-N_bag) with N_bag = 12 + correction?
# E_bag = 6.0 (derived). What if n = 2*E_bag = 12? Then:
m3_2Ebag = base / PHI**12
print(f"\nm_e * alpha^2 * phi^(-2*E_bag) = m_e*a^2*phi^-12 = {m3_2Ebag:.5f} eV")
print(f"  Error vs m3: {abs(m3_2Ebag-m3_exp)/m3_exp*100:.1f}%")

# Or n = deg = 12?
m3_deg = base / PHI**12
print(f"\nm_e * alpha^2 * phi^(-deg) = m_e*a^2*phi^-12 = {m3_deg:.5f} eV")
print(f"  Error vs m3: {abs(m3_deg-m3_exp)/m3_exp*100:.1f}%")
# That's 70% off. Too much.

# ===================================================================
print("\n\n--- BEST HYBRID: alpha^2 * phi^-13 for m3, find m2 ---")
# ===================================================================
# If m3 ~ m_e * alpha^2 / phi^13 (5.4% error):
# What k gives m2?
k_m2 = np.log(base/(m2_exp * PHI**13)) / np.log(PHI)
# m2 = base / phi^(13+k) => phi^k = base/(m2*phi^13) => k = log(base/(m2*phi^13))/log(phi)
k_m2_v2 = np.log(m3_a / m2_exp) / np.log(PHI)
print(f"For m2 = {m2_exp:.5f} eV:")
print(f"  Need k = {k_m2_v2:.3f}")
print(f"  Nearest integer: {round(k_m2_v2)}")
print(f"  Nearest: k=4 gives phi^4={PHI**4:.3f}, m2={m3_a/PHI**4:.5f} (err: {abs(m3_a/PHI**4-m2_exp)/m2_exp*100:.1f}%)")
print(f"  Nearest: k=3 gives phi^3={PHI**3:.3f}, m2={m3_a/PHI**3:.5f} (err: {abs(m3_a/PHI**3-m2_exp)/m2_exp*100:.1f}%)")

# m2/m3 ratio
r = m2_exp / m3_exp
print(f"\nm2/m3 ratio: {r:.5f}")
print(f"  phi^(-3.72) = {PHI**(-3.72):.5f}")
print(f"  1/phi^4 = {1/PHI**4:.5f} (err: {abs(1/PHI**4-r)/r*100:.1f}%)")
print(f"  1/phi^3 = {1/PHI**3:.5f} (err: {abs(1/PHI**3-r)/r*100:.1f}%)")
print(f"  2*sqrt(alpha) = {2*np.sqrt(ALPHA):.5f} (err: {abs(2*np.sqrt(ALPHA)-r)/r*100:.1f}%)")

print("\n" + "="*70)
print("FINAL VERDICT")
print("="*70)
print(f"""
The proposal m_nu = m_e * alpha^2 * phi^-(13+k) is:

BETTER than m_e*phi^-34 (old proposal):
  - Uses derived quantity |Bag|=13 instead of numerological 34
  - m3 error: 5.4% vs 19%
  - Respects cosmological bound

COMPARABLE to m_e*(alpha/phi)^3 (EXP-218):
  - m3 error: 5.4% vs 5.4% (essentially same!)
  - But EXP-218 has cleaner formula (one expression)

PROBLEMATIC for m2:
  - k=3 is NOT a PMNS exponent (it's CKM Cabibbo)
  - m2 error: 42% (user's k) or 272% (correct PMNS k)
  - The mass-angle "unification" is the factor-3 rule repackaged

STATUS: PATTERN for m3, NEGATIVE for m2 hierarchy.
The alpha^2 * phi^-13 combination for m3 is a decent pattern,
but the generation splitting (k values) doesn't work.
""")

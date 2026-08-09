"""
exp311_cc_refined.py
=====================
REFINED COSMOLOGICAL CONSTANT FORMULAS

From exp310_cosmological_constant_v2: Lambda_P ~ alpha^57, factor 1.793 off.
KEY INSIGHT: ln(1.793)/ln(alpha) ~ -alpha_s = -1/(2*phi^3)!

This suggests: z_Lambda = 57 - alpha_s = (N/2 - N_gen) - 1/(2*phi^3)

Tests:
  1. z = 57 - alpha_s (= 58.5 - phi in exact form)
  2. z = 57 - alpha_s with Galois correction
  3. Galois vacuum energy cancellation
  4. Seesaw-type CC from SM x Galois sectors
  5. Various algebraic alternatives

CATEGORY: EXPLORATORY
"""

import math

# === FRAMEWORK CONSTANTS ===
a1 = 5
b1 = a1 + 1  # = 6
PHI = (1 + math.sqrt(a1)) / 2
phi_prime = (1 - math.sqrt(a1)) / 2  # = -1/PHI
N = 120  # = a1!
N_gen = 3
N_eig = 9
rank_E8 = 8
h_E8 = a1 * b1  # = 30
degree = 12

# Couplings
alpha = 7.2973525693e-3
alpha_s = 1 / (2 * PHI**3)  # framework: 0.117851...
sin2_tW = b1 / (a1**2 + 1)  # = 6/26

# Masses
m_e_MeV = 0.51099895
m_e_GeV = m_e_MeV * 1e-3
M_P_GeV = 1.22089e19

# Planck units & CC observed
HBAR = 1.054571817e-34
C = 299792458
G = 6.67430e-11
l_P = math.sqrt(HBAR * G / C**3)
E_CHARGE = 1.602176634e-19

H_0_SI = 67.36 * 1000 / 3.0857e22
Omega_Lambda = 0.6847
Lambda_m2 = 3 * Omega_Lambda * H_0_SI**2 / C**2
Lambda_Planck = Lambda_m2 * l_P**2

z_Planck = 4 * PHI**2
n_nu = 35

# Exact exponent from data
z_exact = math.log(Lambda_Planck) / math.log(alpha)

# Lambda^{1/4} in GeV
hbar_c = HBAR * C
GeV_to_J = 1e9 * E_CHARGE
GeV_inv_to_m = hbar_c / GeV_to_J
rho_crit = 3 * H_0_SI**2 / (8 * math.pi * G)
rho_Lambda_SI = Omega_Lambda * rho_crit
rho_Lambda_J = rho_Lambda_SI * C**2
rho_Lambda_GeV4 = rho_Lambda_J / GeV_to_J * GeV_inv_to_m**3
Lambda_quarter_GeV = rho_Lambda_GeV4**0.25
Lambda_quarter_meV = Lambda_quarter_GeV * 1e12

print("=" * 72)
print("EXP311: REFINED COSMOLOGICAL CONSTANT FORMULAS")
print("=" * 72)

print(f"\n--- Observed values ---")
print(f"Lambda_P = {Lambda_Planck:.6e}")
print(f"z_exact = ln(Lambda_P)/ln(alpha) = {z_exact:.6f}")
print(f"Lambda^{{1/4}} = {Lambda_quarter_meV:.4f} meV")

# =====================================================================
# SECTION 1: THE z = 57 - alpha_s FORMULA
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 1: z_Lambda = 57 - alpha_s")
print("=" * 72)

z_57 = 57.0
z_57_minus_as = 57.0 - alpha_s

print(f"\nalpha_s (framework) = 1/(2*phi^3) = {alpha_s:.10f}")
print(f"z = 57 = N/2 - N_gen     -> alpha^57 = {alpha**57:.6e}")
print(f"z = 57 - alpha_s = {z_57_minus_as:.10f}")
print(f"alpha^(57 - alpha_s) = {alpha**z_57_minus_as:.6e}")
print(f"Lambda_P observed    = {Lambda_Planck:.6e}")
print(f"")
print(f"Ratio Lambda_P / alpha^(57-alpha_s) = {Lambda_Planck / alpha**z_57_minus_as:.6f}")
err_pct = (alpha**z_57_minus_as - Lambda_Planck) / Lambda_Planck * 100
print(f"Error: {err_pct:+.4f}%")
print(f"")
print(f"z_exact - (57 - alpha_s) = {z_exact - z_57_minus_as:.6f}")
print(f"  (positive = formula gives z too large, Lambda too small)")

# What is 57 - alpha_s in exact form?
# alpha_s = 1/(2*phi^3) = 1/(4*phi+2)
# Rationalized: (4*phi'+2) / ((4*phi+2)(4*phi'+2))
# (4phi+2)(4phi'+2) = 16*phi*phi'+8*(phi+phi')+4 = -16+8+4 = -4
# So 1/(4phi+2) = (4phi'+2)/(-4) = -(phi'+1/2) = -(1-phi+1/2) = phi - 3/2
# alpha_s = phi - 3/2

alpha_s_exact = PHI - 1.5
print(f"\nExact form: alpha_s = 1/(2*phi^3) = phi - 3/2")
print(f"  Check: phi - 3/2 = {alpha_s_exact:.10f}")
print(f"  alpha_s = {alpha_s:.10f}")
print(f"  Match: {abs(alpha_s - alpha_s_exact):.2e}")

print(f"\nTherefore: z_Lambda = 57 - (phi - 3/2) = 58.5 - phi = (117 - 2*phi)/2")
z_formula = 58.5 - PHI
print(f"  z = 58.5 - phi = {z_formula:.10f}")
print(f"  z_exact =        {z_exact:.10f}")
print(f"  Difference: {z_exact - z_formula:.6f}")

# Observational uncertainty
H0_err = 0.54  # km/s/Mpc
Omega_err = 0.0073
# Propagate to z
H_0_lo = (67.36 - H0_err) * 1000 / 3.0857e22
H_0_hi = (67.36 + H0_err) * 1000 / 3.0857e22
Omega_lo = Omega_Lambda - Omega_err
Omega_hi = Omega_Lambda + Omega_err
Lambda_lo = 3 * Omega_lo * H_0_lo**2 / C**2 * l_P**2
Lambda_hi = 3 * Omega_hi * H_0_hi**2 / C**2 * l_P**2
z_lo = math.log(Lambda_hi) / math.log(alpha)  # larger Lambda -> smaller |z|
z_hi = math.log(Lambda_lo) / math.log(alpha)  # smaller Lambda -> larger |z|
z_err = abs(z_hi - z_lo) / 2

print(f"\n--- Observational uncertainty ---")
print(f"z_exact = {z_exact:.4f} +/- {z_err:.4f}")
print(f"z_formula = {z_formula:.4f}")
print(f"Deviation: {abs(z_exact - z_formula)/z_err:.2f} sigma")
print(f"57 deviation: {abs(z_exact - 57)/z_err:.2f} sigma")

# =====================================================================
# SECTION 2: PHYSICAL INTERPRETATION
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 2: PHYSICAL INTERPRETATION")
print("=" * 72)

print(f"""
Formula: Lambda_P = alpha^(N/2 - N_gen - alpha_s)
       = alpha^(57 - 1/(2*phi^3))
       = alpha^(58.5 - phi)

Decomposition:
  N/2 = 60 = half the 600-cell vertices (geometric)
  N_gen = 3 = number of generations (Nyquist on icosahedron)
  alpha_s = 1/(2*phi^3) = QCD coupling (one-loop correction)

Interpretation:
  - The INTEGER part: 57 = N/2 - N_gen = 60 - 3
    This is the "tree-level" CC exponent from the combinatorics
    of the 600-cell minus the fermionic (generation) correction.

  - The FRACTIONAL correction: -alpha_s = -1/(2*phi^3)
    This is a one-loop QCD correction. The strong coupling
    modifies the vacuum energy by alpha_s in the exponent.

  - Equivalently: z = 58.5 - phi = (117/2) - phi
    117 = 2*57 + 3 = 2*(N/2 - N_gen) + N_gen = N - N_gen = 117
    So z = (N - N_gen)/2 - phi = (120 - 3)/2 - phi

  NICE FORM: z = (N - N_gen)/2 - phi = 117/2 - phi
""")

# Check: (N - N_gen)/2 - phi
z_nice = (N - N_gen) / 2.0 - PHI
print(f"z = (N - N_gen)/2 - phi = ({N} - {N_gen})/2 - phi = {z_nice:.10f}")
print(f"z_exact =                  {z_exact:.10f}")
print(f"Match: {abs(z_nice - z_exact):.6f}")
print(f"Deviation: {abs(z_nice - z_exact)/z_err:.2f} sigma")

# =====================================================================
# SECTION 3: COMPARISON OF CANDIDATE FORMULAS
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 3: SYSTEMATIC COMPARISON")
print("=" * 72)

candidates = []

# The new formula
candidates.append(("57 - alpha_s = (N-N_gen)/2 - phi", z_57_minus_as))

# Pure integers
candidates.append(("57 = N/2 - N_gen", 57.0))
candidates.append(("56 = 2*h_E8 - N_gen - 1", 56.0))

# Z[phi] elements close to z
candidates.append(("31 + 16*phi", 31 + 16*PHI))
candidates.append(("86 - 18*phi", 86 - 18*PHI))
candidates.append(("49 + 5*phi (5a+6b form, a=49/5?)", 49 + 5*PHI))

# Framework expressions
candidates.append(("4*z_Planck + N_gen*phi^3", 4*z_Planck + N_gen*PHI**3))
candidates.append(("a1*z_Planck + b1 + alpha_s", a1*z_Planck + b1 + alpha_s))
candidates.append(("N/2 - N_gen - sin2_tW", N/2 - N_gen - sin2_tW))
candidates.append(("N/2 - N_gen*phi", N/2 - N_gen*PHI))
candidates.append(("rank_E8 * (b1 + alpha_s + 1)", rank_E8*(b1 + alpha_s + 1)))
candidates.append(("2*n_nu - degree - alpha_s", 2*n_nu - degree - alpha_s))
candidates.append(("N/2 - pi", N/2 - math.pi))
candidates.append(("57 - 2*alpha_s", 57 - 2*alpha_s))
candidates.append(("57 - alpha_s/2", 57 - alpha_s/2))
candidates.append(("4*n_nu/phi^2 + 1", 4*n_nu/PHI**2 + 1))

# phi-based
candidates.append(("phi^(a1+3) - phi^2", PHI**(a1+3) - PHI**2))
candidates.append(("N_eig*b1 + N_gen + alpha_s", N_eig*b1 + N_gen + alpha_s))

# Sort by closeness to z_exact
candidates.sort(key=lambda x: abs(x[1] - z_exact))

print(f"\nTarget: z_exact = {z_exact:.6f}")
print(f"Obs uncertainty: +/- {z_err:.4f}")
print(f"\n{'Formula':<45} | {'z_pred':>10} | {'z-z_ex':>10} | {'sigma':>8} | {'Lam_P':>12}")
print("-" * 100)
for name, z_pred in candidates:
    dz = z_pred - z_exact
    sigma = abs(dz) / z_err if z_err > 0 else float('inf')
    lam_p = alpha**z_pred
    marker = " ***" if sigma < 1 else (" **" if sigma < 3 else (" *" if sigma < 10 else ""))
    print(f"{name:<45} | {z_pred:10.6f} | {dz:+10.6f} | {sigma:8.2f} | {lam_p:12.4e}{marker}")

# =====================================================================
# SECTION 4: THE z = 57 - alpha_s FORMULA IN DETAIL
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 4: DETAILED ANALYSIS OF z = 57 - alpha_s")
print("=" * 72)

# Lambda_P = alpha^(57 - alpha_s) in terms of framework quantities
# alpha = solution of 2*pi*x^2 - 20*phi^4*x + 1 = 0
# alpha_s = 1/(2*phi^3)
# 57 = N/2 - N_gen

# Express as: Lambda_P = alpha^57 * alpha^(-alpha_s)
# = alpha^57 * exp(-alpha_s * ln(alpha))
# = alpha^57 * exp(alpha_s * |ln(alpha)|)

print(f"\nLambda_P = alpha^57 * alpha^(-alpha_s)")
print(f"  alpha^57 = {alpha**57:.6e}")
print(f"  alpha^(-alpha_s) = {alpha**(-alpha_s):.6f}")
print(f"  Product = {alpha**57 * alpha**(-alpha_s):.6e}")
print(f"  Observed = {Lambda_Planck:.6e}")
print(f"  Ratio: {Lambda_Planck / (alpha**57 * alpha**(-alpha_s)):.6f}")

# The correction factor alpha^{-alpha_s}
corr = alpha**(-alpha_s)
print(f"\nCorrection factor: alpha^{{-alpha_s}} = {corr:.6f}")
print(f"  ln(correction) = -alpha_s * ln(alpha) = {-alpha_s * math.log(alpha):.6f}")
print(f"  = alpha_s * |ln(alpha)| = {alpha_s * abs(math.log(alpha)):.6f}")

# All quantities from a1 = 5
print(f"\n--- Everything from a1 = {a1} ---")
print(f"  phi = (1+sqrt(a1))/2 = {PHI:.10f}")
print(f"  N = a1! = {N}")
print(f"  N_gen = 3 (Nyquist)")
print(f"  alpha_s = 1/(2*phi^3) = {alpha_s:.10f}")
print(f"  alpha from 2*pi*x^2 - 4*a1*phi^4*x + 1 = 0")
print(f"  z = N/2 - N_gen - alpha_s = 60 - 3 - {alpha_s:.6f} = {z_57_minus_as:.6f}")
print(f"  Lambda_P = alpha^z = {alpha**z_57_minus_as:.6e}")
print(f"  Observed = {Lambda_Planck:.6e}")

# =====================================================================
# SECTION 5: GALOIS VACUUM ENERGY CANCELLATION
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 5: GALOIS VACUUM ENERGY CANCELLATION")
print("=" * 72)

print(f"""
The SM and Galois sectors both contribute to vacuum energy.
If the total CC is a DIFFERENCE (partial cancellation):

  rho_Lambda = rho_vac(SM) - rho_vac(Galois)

In the framework, SM uses phi, Galois uses phi' = -1/phi.
The spectral action S = Tr f(D/Lambda) sums eigenvalues.

SM eigenvalues scale as phi^n (growing).
Galois eigenvalues scale as phi'^n = (-1)^n/phi^n (shrinking).

For the vacuum energy: Tr(D^4) ~ sum phi^(4n)  (SM)
                                 ~ sum phi'^(4n) (Galois)
Since phi'^4 = phi^(-4) > 0 (even power), Galois terms are POSITIVE.

No cancellation from sign: both sectors contribute positively.
However, the MAGNITUDE differs enormously.
""")

# Compute ratio of Galois to SM vacuum energies
# For each fermion with exponent n:
# SM contributes ~ phi^{4n}, Galois ~ phi^{-4n}
# Ratio for single fermion: phi^{-4n} / phi^{4n} = phi^{-8n}
print(f"Ratio rho_Galois / rho_SM per fermion:")
fermions = {
    'e':   (0, 0),  'mu':  (1, 1),  'tau': (1, 2),
    'u':   (3, -2), 'c':   (2, 1),  't':   (4, 1),
    'd':   (1, 0),  's':   (1, 1),  'b':   (-1, 4),
}

total_ratio = 0
for name, (a, b) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
    n = 5*a + 6*b
    ratio_n = PHI**(-8*n) if n != 0 else 1.0
    total_ratio += ratio_n
    print(f"  {name:>5}: n={n:3d}, phi^{{-8n}} = {ratio_n:.4e}")

print(f"\nSum of ratios: {total_ratio:.6e}")
print(f"Dominated by electron (n=0): ratio = 1")
print(f"Second: down quark (n=5): ratio = {PHI**(-40):.4e}")
print(f"")
print(f"Since electron has n=0, its SM and Galois contributions are EQUAL.")
print(f"The difference comes from other fermions where SM >> Galois.")
print(f"So rho_vac(SM) - rho_vac(Galois) ~ rho_vac(SM) * (1 - tiny)")
print(f"NO significant cancellation from Galois sector.")

# =====================================================================
# SECTION 6: LOOP-CORRECTED CC
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 6: LOOP CORRECTION INTERPRETATION")
print("=" * 72)

print(f"""
If Lambda_P = alpha^(N/2 - N_gen - alpha_s), the alpha_s correction
could represent a QCD vacuum energy contribution.

In QFT, the vacuum energy gets corrections from all interactions:
  rho_vac = rho_0 + rho_QCD + rho_EM + rho_weak + ...

The QCD contribution is proportional to alpha_s.
In the 600-cell framework, this appears as:
  z_CC = z_CC^(tree) + delta_z_QCD
       = (N/2 - N_gen) + (-alpha_s)
       = 57 - 0.11785

The SIGN: alpha_s contribution REDUCES the exponent (makes Lambda larger).
This is physically correct: QCD vacuum condensate contributes POSITIVE
vacuum energy, partially canceling the negative tree-level contribution.
""")

# Test: is the correction exactly alpha_s or could it be different?
# The exact delta from 57 is:
delta_exact = 57 - z_exact
print(f"Exact correction needed: 57 - z_exact = {delta_exact:.6f}")
print(f"alpha_s = {alpha_s:.6f}")
print(f"Difference: {delta_exact - alpha_s:.6f}")

# What is this difference in framework terms?
diff = delta_exact - alpha_s
print(f"\nResidual after alpha_s: {diff:.6f}")
print(f"  = {diff/alpha:.4f} * alpha")
print(f"  = {diff/alpha**2:.4f} * alpha^2")
print(f"  = {diff/alpha_s:.6f} * alpha_s")
print(f"  = {diff/(alpha*alpha_s):.4f} * alpha*alpha_s")

# Could it be alpha_s + alpha?
z_as_plus_a = 57 - alpha_s - alpha
print(f"\nTest z = 57 - alpha_s - alpha:")
print(f"  z = {z_as_plus_a:.6f}")
print(f"  z_exact = {z_exact:.6f}")
print(f"  Difference: {z_as_plus_a - z_exact:.6f}")
print(f"  sigma: {abs(z_as_plus_a - z_exact)/z_err:.2f}")

# Or 57 - alpha_s - alpha * sin^2(tW)?
z_test = 57 - alpha_s - alpha * sin2_tW
print(f"\nTest z = 57 - alpha_s - alpha*sin^2(tW):")
print(f"  z = {z_test:.6f}")
print(f"  Difference: {z_test - z_exact:.6f}")
print(f"  sigma: {abs(z_test - z_exact)/z_err:.2f}")

# Or 57 - alpha_s + alpha_s^2?
z_test2 = 57 - alpha_s + alpha_s**2
print(f"\nTest z = 57 - alpha_s + alpha_s^2 (QCD 2-loop):")
print(f"  z = {z_test2:.6f}")
print(f"  Difference: {z_test2 - z_exact:.6f}")
print(f"  sigma: {abs(z_test2 - z_exact)/z_err:.2f}")

# Systematic search of corrections
print(f"\n--- Systematic correction search ---")
corrections = {
    "alpha_s only": alpha_s,
    "alpha_s + alpha": alpha_s + alpha,
    "alpha_s + alpha*sin2_tW": alpha_s + alpha * sin2_tW,
    "alpha_s - alpha_s^2": alpha_s - alpha_s**2,
    "alpha_s + alpha_s^2": alpha_s + alpha_s**2,
    "alpha_s*(1+alpha_s)": alpha_s * (1 + alpha_s),
    "alpha_s/(1-alpha_s)": alpha_s / (1 - alpha_s),
    "alpha_s + alpha/(2*pi)": alpha_s + alpha / (2*math.pi),
    "alpha_s + alpha*N_gen": alpha_s + alpha * N_gen,
    "phi - 3/2 (exact alpha_s)": PHI - 1.5,
    "phi - 3/2 + alpha^2": PHI - 1.5 + alpha**2,
}

print(f"\nTarget correction: {delta_exact:.8f}")
print(f"\n{'Correction formula':<35} | {'Value':>12} | {'delta':>12} | {'sigma':>8}")
print("-" * 75)
for name, val in sorted(corrections.items(), key=lambda x: abs(x[1] - delta_exact)):
    z_test = 57 - val
    delta = val - delta_exact
    sigma = abs(z_test - z_exact) / z_err
    marker = " ***" if sigma < 1 else (" **" if sigma < 3 else "")
    print(f"{name:<35} | {val:12.8f} | {delta:+12.8f} | {sigma:8.2f}{marker}")

# =====================================================================
# SECTION 7: THE EXACT Z[phi] FORMULA
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 7: EXACT ALGEBRAIC FORMULA")
print("=" * 72)

# z = 58.5 - phi = (117 - 2*phi) / 2
# But 117/2 is not in Z[phi]. Can we write this differently?
# z = 57 - (phi - 3/2) where phi - 3/2 = alpha_s
# This involves half-integers.

# Alternative: 2*z = 117 - 2*phi = (117, -2) in Z[phi]
# Tr(2z) = 2*117 + (-2) = 232
# N(2z) = 117^2 + 117*(-2) - (-2)^2 = 13689 - 234 - 4 = 13451
# 13451 = ? Check: 13451/a1 = 2690.2, 13451/N = 112.09
two_z = 117 - 2*PHI
print(f"\n2*z_Lambda = 117 - 2*phi = {two_z:.6f}")
print(f"In Z[phi]: (117, -2)")
print(f"Tr(2z) = {2*117 + (-2)} = 232")
print(f"N(2z) = {117**2 + 117*(-2) - (-2)**2}")
print(f"2*z_exact = {2*z_exact:.6f}")
print(f"Difference: {abs(two_z - 2*z_exact):.6f}")

# Is 117 interesting?
print(f"\n117 = N - N_gen = {N} - {N_gen}")
print(f"117 = 9 * 13 = N_eig * 13")
print(f"117 = 3 * 39 = N_gen * 39")
print(f"13 = degree + 1 = {degree+1}")
print(f"39 = 3 * 13")
print(f"So 117 = N_eig * (degree + 1) = 9 * 13")

# Better: z = (N - N_gen) / 2 - phi
# Physical: (total vertices minus fermionic modes) / 2, corrected by phi
print(f"\nNICE FORM: z_Lambda = (N - N_gen)/2 - phi")
print(f"  = (a1! - N_gen)/2 - phi")
print(f"  = (120 - 3)/2 - phi")
print(f"  = 58.5 - phi")
print(f"  = {(N - N_gen)/2 - PHI:.10f}")
print(f"  z_exact = {z_exact:.10f}")

# Alternative: z = N/2 - N_gen - 1/(2*phi^3)
# = 60 - 3 - 1/(4*phi + 2)
# = 57 - 1/(4*phi + 2)
print(f"\nEQUIVALENT: z_Lambda = N/2 - N_gen - 1/(2*phi^3)")
print(f"  = 60 - 3 - 1/{2*PHI**3:.6f}")
print(f"  = {N/2 - N_gen - 1/(2*PHI**3):.10f}")

# =====================================================================
# SECTION 8: NEUTRINO-CC CONNECTION REFINED
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 8: NEUTRINO-CC CONNECTION REFINED")
print("=" * 72)

m_3_GeV = 2 * m_e_GeV / PHI**n_nu
m_3_meV = m_3_GeV * 1e12

print(f"\nm_3 = 2*m_e/phi^35 = {m_3_meV:.4f} meV")
print(f"Lambda^{{1/4}} (obs) = {Lambda_quarter_meV:.4f} meV")
print(f"Ratio m_3/Lambda^{{1/4}} = {m_3_meV/Lambda_quarter_meV:.4f}")

# From exp310: Lambda^{1/4} ~ m_3 / (2*z_Planck)
# z_Planck = 4*phi^2 ~ 10.472
# m_3/(2*z_P) = 49.53 / 20.944 = 2.365 meV
# Observed: 2.24 meV
# Error: 5.6%

divisors = {
    "2*z_Planck = 8*phi^2": 2 * z_Planck,
    "4*a1 = 20": 4.0 * a1,
    "m_3/Lambda^{1/4} exact": m_3_meV / Lambda_quarter_meV,
    "4*a1 + alpha_s": 4*a1 + alpha_s,
    "4*a1 + 2*alpha_s": 4*a1 + 2*alpha_s,
    "4*a1*phi/phi' (Galois)": abs(4*a1*PHI/phi_prime),
    "N_eig * phi^2": N_eig * PHI**2,
    "a1*(phi^2+N_gen)": a1*(PHI**2 + N_gen),
    "4*a1*(1+alpha_s/a1)": 4*a1*(1+alpha_s/a1),
    "8*phi^2 + 1": 8*PHI**2 + 1,
    "degree*phi + N_gen": degree*PHI + N_gen,
    "2*degree*phi^0.5": 2*degree*PHI**0.5,
}

print(f"\n{'Divisor formula':<35} | {'Value':>10} | {'L^1/4 pred':>10} | {'Error%':>8}")
print("-" * 75)
for name, div in sorted(divisors.items(), key=lambda x: abs(m_3_meV/x[1] - Lambda_quarter_meV)):
    if div > 0:
        pred = m_3_meV / div
        err = (pred - Lambda_quarter_meV) / Lambda_quarter_meV * 100
        print(f"{name:<35} | {div:10.4f} | {pred:10.4f} | {err:+8.4f}%")

# =====================================================================
# SECTION 9: CONNECTING z_CC TO NEUTRINO FORMULA
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 9: CC FROM m_e/M_P HIERARCHY")
print("=" * 72)

# m_e/M_P = alpha^{z_Planck} = alpha^{4*phi^2}
# Lambda_P = alpha^{z_CC}
# z_CC / z_Planck = ?
ratio_z = z_exact / z_Planck
print(f"\nz_CC = {z_exact:.6f}")
print(f"z_Planck = {z_Planck:.6f}")
print(f"z_CC / z_Planck = {ratio_z:.6f}")

# Is this ratio something?
print(f"\nz_CC/z_P candidates:")
print(f"  a1 + phi/2 = {a1 + PHI/2:.6f} (ratio {(a1+PHI/2):.6f} vs {ratio_z:.6f})")
print(f"  a1 + 1/(2*phi) = {a1 + 1/(2*PHI):.6f}")
print(f"  b1 - phi/N_gen = {b1 - PHI/N_gen:.6f}")
print(f"  a1 + alpha_s = {a1 + alpha_s:.6f}")

# z_CC = z_P * (a1 + alpha_s)?
z_test_ratio = z_Planck * (a1 + alpha_s)
print(f"\nTest: z_CC = z_P * (a1 + alpha_s):")
print(f"  = {z_Planck:.4f} * {a1+alpha_s:.6f} = {z_test_ratio:.6f}")
print(f"  z_exact = {z_exact:.6f}")
print(f"  Error: {abs(z_test_ratio - z_exact):.4f} ({abs(z_test_ratio - z_exact)/z_err:.1f} sigma)")

# z_CC = (z_P + 1) * a1 + phi?
z_test2 = (z_Planck + 1) * a1 + PHI
print(f"\nTest: z_CC = (z_P+1)*a1 + phi = {z_test2:.6f} (off by {z_test2-z_exact:.3f})")

# z_CC = z_P^2 / phi^2?
z_test3 = z_Planck**2 / PHI**2
print(f"Test: z_CC = z_P^2/phi^2 = {z_test3:.6f} (off by {z_test3-z_exact:.3f})")

# =====================================================================
# SECTION 10: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 10: HONEST ASSESSMENT")
print("=" * 72)

# Final numbers
z_best_formula = (N - N_gen)/2.0 - PHI  # = 58.5 - phi = 57 - alpha_s
lam_best = alpha**z_best_formula
lam_obs = Lambda_Planck
lam_err = abs(lam_best - lam_obs) / lam_obs * 100
z_dev = abs(z_best_formula - z_exact) / z_err

print(f"""
BEST FORMULA: Lambda_P = alpha^((N - N_gen)/2 - phi)
            = alpha^(N/2 - N_gen - alpha_s)
            = alpha^(57 - 1/(2*phi^3))

Numerical values:
  z_formula = {z_best_formula:.6f}
  z_exact   = {z_exact:.6f}
  z_err     = +/- {z_err:.4f}
  Deviation = {z_dev:.2f} sigma

  Lambda_pred = {lam_best:.6e}
  Lambda_obs  = {lam_obs:.6e}
  Error on Lambda: {lam_err:.2f}%

CATEGORIZATION:
  - The EXPONENT 57 = N/2 - N_gen: PATTERN (motivated but not derived)
  - The CORRECTION -alpha_s: PATTERN (QCD loop correction, plausible)
  - Combined formula z = 57 - alpha_s: PATTERN
  - Match to observation: {z_dev:.1f} sigma (EXCELLENT fit)

WHAT'S MISSING FOR DERIVATION:
  1. No mechanism from spectral action produces exponent 57
  2. No computation shows why QCD corrects CC at this level
  3. The formula is discovered empirically, not derived

HOWEVER:
  - Uses ONLY quantities already in the framework (N, N_gen, alpha_s)
  - No new parameters
  - Error is within observational uncertainty
  - Physical interpretation is plausible (tree + 1-loop)

UPGRADE: from "factor 1.8 off" (exp310) to "within 1 sigma" (this exp)
CATEGORY: STRONG PATTERN (but NOT DERIVED)
""")

print("=" * 72)
print("EXP311 COMPLETE")
print("=" * 72)

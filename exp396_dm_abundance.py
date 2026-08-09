"""
EXP-396: Dark Matter Abundance Ratio — Derivation Attempt
==========================================================

GOAL: Derive Omega_DM/Omega_b = 7 - phi = 5.382 (obs: 5.36, 0.41%)
from first principles of the 600-cell framework.

KNOWN (DERIVED):
  - Galois sector IS the dark sector (alpha' complex, alpha_s' < 0)
  - 9 stable fermion species with masses m_f' = m_e^2/m_f
  - Gravitational freeze-in is the UNIQUE viable production mechanism
  - Delta N_eff ~ 0 automatically

KNOWN (PATTERN):
  - Omega_DM/Omega_b = 7 - phi = (13-sqrt(5))/2 = 5.382 (0.41%)
  - T_reh ~ 1.4e17 GeV required but NOT derived

NEW APPROACHES:
  1. Algebraic decomposition of 7-phi in framework constants
  2. Galois mass sum ratios
  3. Spectral invariant search
  4. Gravitational freeze-in with spectral action T_reh
  5. Partition function argument

Author: Claude (exp396)
Date: 2026-02-25
"""

import math

# ===== FRAMEWORK CONSTANTS =====
a1 = 5
b1 = 6
PHI = (1 + math.sqrt(5)) / 2
PHI_C = (1 - math.sqrt(5)) / 2  # phi' = -1/phi
N = 120
N_eig = 9
N_gen = 3
DEG = 12
rank_E8 = 8
h_E8 = 30

A_eq = 2 * math.pi
B_eq = -4 * a1 * PHI**4
C_eq = 1.0
disc = B_eq**2 - 4 * A_eq * C_eq
ALPHA = (-B_eq - math.sqrt(disc)) / (2 * A_eq)
ALPHA_S = 1 / (2 * PHI**3)
sin2_tW = b1 / (a1**2 + 1)  # = 6/26

m_e_eV = 0.51099895e6  # eV

# Experimental
Omega_ratio_exp = 5.36  # Planck 2018: Omega_DM/Omega_b
eta_B_exp = 6.1e-10  # baryon-to-photon ratio
m_proton = 938.272e6  # eV

print("=" * 72)
print("EXP-396: DARK MATTER ABUNDANCE RATIO")
print("=" * 72)

# Target
target = 7 - PHI
print(f"\n  Target: 7 - phi = {target:.6f}")
print(f"  Experiment: {Omega_ratio_exp:.2f} +/- 0.05")
print(f"  Error: {abs(target-Omega_ratio_exp)/Omega_ratio_exp*100:.2f}%")


# =====================================================================
# PART 1: ALGEBRAIC DECOMPOSITIONS OF 7-phi
# =====================================================================
print("\n" + "-" * 72)
print("PART 1: Algebraic decompositions of 7 - phi")
print("-" * 72)

decomps = [
    ("7 - phi", 7 - PHI),
    ("b1 + 1 - phi", b1 + 1 - PHI),
    ("a1 + 2 - phi", a1 + 2 - PHI),
    ("b1 + phi'^2", b1 + PHI_C**2),
    ("a1 + phi' + 2", a1 + PHI_C + 2),
    ("(a1^2+1)/2 - phi + 1/2", (a1**2+1)/2 - PHI + 0.5),
    ("(13-sqrt(5))/2", (13 - math.sqrt(5))/2),
]

print(f"\n  {'Expression':>30} {'Value':>12} {'Match':>8}")
print(f"  {'-'*55}")
for name, val in decomps:
    match = "YES" if abs(val - target) < 1e-10 else "no"
    print(f"  {name:>30} {val:>12.6f} {match:>8}")

# Z[phi] properties
trace = 2*7 - (PHI + PHI_C)  # = 14 - 1 = 13
norm = (7 - PHI) * (7 - PHI_C)  # = (7-phi)(7+1/phi) = 49-7+... let me compute
norm_val = (7 - PHI) * (7 - PHI_C)
print(f"\n  Z[phi] properties of 7-phi:")
print(f"    Trace = (7-phi) + (7-phi') = {trace:.0f}")
print(f"    Norm  = (7-phi)(7-phi') = {norm_val:.0f}")
print(f"    13 = (a1^2+1)/2 [Fibonacci F_7]")
print(f"    41 = prime [= a1^2 + (a1-1)^2 = 25+16]")


# =====================================================================
# PART 2: GALOIS MASS SUMS
# =====================================================================
print("\n" + "-" * 72)
print("PART 2: Galois mass sums and ratios")
print("-" * 72)

# SM fermion exponents n_f = 5*a + 6*b
# (a,b): e(0,0), mu(1,1), tau(1,2), u(3,-2), c(2,1), t(4,1), d(1,0), s(1,1), b(-1,4)
n_f = {
    'e': 0, 'mu': 11, 'tau': 17,
    'u': 3, 'c': 16, 't': 26,
    'd': 5, 's': 11, 'b': 19
}

print(f"\n  SM and Galois masses:")
print(f"  {'Fermion':>8} {'n_f':>4} {'m_SM (MeV)':>12} {'m_Galois (MeV)':>14}")
print(f"  {'-'*42}")

sum_mf = 0  # sum of SM masses
sum_mfp = 0  # sum of Galois masses
for name in ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']:
    nf = n_f[name]
    m_sm = m_e_eV * PHI**nf / 1e6  # MeV
    m_gal = m_e_eV / PHI**nf / 1e6  # MeV (= m_e^2/m_f in natural units)
    sum_mf += m_sm
    sum_mfp += m_gal
    print(f"  {name:>8} {nf:>4} {m_sm:>12.4f} {m_gal:>14.6f}")

print(f"\n  Sum SM: {sum_mf:.2f} MeV")
print(f"  Sum Galois: {sum_mfp:.6f} MeV = {sum_mfp*1e3:.3f} keV")
print(f"  Ratio sum_Galois/m_e: {sum_mfp*1e6/m_e_eV:.6f}")

# Galois mass sum in units of m_e
gal_sum_over_me = sum(1/PHI**n_f[f] for f in n_f)
print(f"  sum(1/phi^n_f) = {gal_sum_over_me:.6f}")

# Is 7-phi related to any mass ratio?
print(f"\n  Testing mass ratios:")
print(f"    sum(m_f')/m_e = {gal_sum_over_me:.6f}")
print(f"    7-phi = {target:.6f}")
print(f"    Ratio: {gal_sum_over_me/target:.6f}")
print(f"    NOT equal.")

# What about sum(m_f')/m_f' for heaviest?
m_e_prime = m_e_eV  # n=0, so m_e' = m_e
print(f"    sum(m_f')/m_e' = {gal_sum_over_me:.6f} (same as above, e'=e)")


# =====================================================================
# PART 3: SPECTRAL INVARIANT SEARCH
# =====================================================================
print("\n" + "-" * 72)
print("PART 3: Spectral invariants near 5.38")
print("-" * 72)

# Try various framework expressions
candidates = []

# Simple expressions
candidates.append(("7 - phi", 7 - PHI))
candidates.append(("b1 - 1/phi", b1 - 1/PHI))
candidates.append(("a1 + phi'**2", a1 + PHI_C**2))
candidates.append(("N_eig - PHI**3", N_eig - PHI**3))  # 9 - 4.236 = 4.764
candidates.append(("a1*phi - PHI**2", a1*PHI - PHI**2))  # 8.09 - 2.618 = 5.472
candidates.append(("DEG / (a1*ALPHA_S)", DEG / (a1*ALPHA_S)))  # 12/(5*0.118) = 20.3
candidates.append(("1/sin2_tW - b1/a1", 1/sin2_tW - b1/a1))  # 26/6 - 6/5 = 3.133
candidates.append(("h_E8/b1 + phi'", h_E8/b1 + PHI_C))  # 5-0.618 = 4.382
candidates.append(("N/(a1*a1-1)", N/((a1*a1-1))))  # 120/24 = 5.0
candidates.append(("N/(a1*(a1-1)+2)", N/(a1*(a1-1)+2)))  # 120/22 = 5.45
candidates.append(("DEG*ALPHA_S/ALPHA", DEG*ALPHA_S/ALPHA))  # 12*0.118/0.0073 = 194
candidates.append(("N_gen*phi + 1/(N_gen*phi)", N_gen*PHI + 1/(N_gen*PHI)))  # 5.06
candidates.append(("(a1^2+1)/(a1-phi)", (a1**2+1)/(a1-PHI)))  # 26/3.382 = 7.69

# More targeted: expressions involving b1, a1, phi
candidates.append(("a1 + 1/PHI**2", a1 + 1/PHI**2))  # 5+0.382 = 5.382!
candidates.append(("a1 + phi - 1 - 1/phi", a1 + PHI - 1 - 1/PHI))  # 5+0 = 5.0

print(f"\n  {'Expression':>35} {'Value':>10} {'Error%':>8}")
print(f"  {'-'*58}")
for name, val in sorted(candidates, key=lambda x: abs(x[1]-target)):
    err = (val - target)/target*100
    marker = " <--" if abs(err) < 0.1 else ""
    print(f"  {name:>35} {val:>10.4f} {err:>8.2f}%{marker}")

# Key finding: a1 + 1/phi^2 = a1 + phi'^2 = 7-phi EXACTLY
print(f"\n  KEY: a1 + 1/phi^2 = a1 + (2-phi) = a1+2-phi = 7-phi EXACTLY")
print(f"       = {a1 + 1/PHI**2:.10f}")
print(f"       Because 1/phi^2 = 2-phi = phi'^2 (Galois conjugate of phi^2=phi+1)")


# =====================================================================
# PART 4: GRAVITATIONAL FREEZE-IN
# =====================================================================
print("\n" + "-" * 72)
print("PART 4: Gravitational freeze-in calculation")
print("-" * 72)

# The gravitational freeze-in yield per species (Ema, Sala, Takahashi 2019):
# Y_DM = (45 * zeta(3) / (2 * pi^4)) * (T_reh / M_Pl) * sum_i g_i
# where g_i = species-dependent factor

# For spin-1/2 fermions produced gravitationally:
# Omega_DM h^2 = sum_i (m_i * Y_i * s_0) / rho_crit
# where s_0 = 2891.2 cm^-3, rho_crit/h^2 = 1.054e4 eV/cm^3

# Simpler: Omega_DM/Omega_b = (sum m_f' * n_f') / (m_p * n_b)
# For gravitational production: all species have same n/s (massless production)
# n_DM_per_species / s = C * (T_reh/M_Pl)^3
# n_b / s = eta_B (given)

# So: Omega_DM/Omega_b = N_DM_species * (C/eta_B) * (T_reh/M_Pl)^3 * <m_DM> / m_p

# The question is: can we determine T_reh from the framework?

# Approach: T_reh from Starobinsky inflation (spectral action)
# The spectral action gives R + R^2/(6*M^2) where M is the scalaron mass.
# Reheating temperature: T_reh ~ 0.55 * (M * M_Pl)^(1/2) (for Starobinsky)
# M_scalaron is related to test function moments... but exp392 showed this is NEGATIVE.

print(f"  Gravitational freeze-in approach:")
print(f"  Omega_DM/Omega_b = (sum m_f') / m_p * N_species * C*(T_reh/M_Pl)^3 / eta_B")
print(f"")
print(f"  PROBLEM: T_reh depends on test function f_4 (NOT derived, exp392 NEGATIVE)")
print(f"")
print(f"  Let's invert: what T_reh is REQUIRED for 7-phi?")

# Required T_reh
# For gravitational freeze-in of N_s species of spin-1/2:
# Y_DM ~ N_s * (T_reh/M_Pl)^3 * 135*zeta(3)/(16*pi^5*g_star)
# with g_star ~ 106.75 (SM)

M_Pl = 1.22e19 * 1e9  # eV (Planck mass)
g_star = 106.75
zeta3 = 1.20206

# Yield coefficient
C_yield = 135 * zeta3 / (16 * math.pi**5 * g_star)
# This gives Y per species ~ C_yield * (T_reh/M_Pl)^3
print(f"  Yield coefficient C = {C_yield:.4e}")

# From Omega_DM/Omega_b = 7-phi:
# (sum m_f') * N_species_eff * C * (T_reh/M_Pl)^3 / (m_p * eta_B) = 7-phi
# N_species_eff accounts for degrees of freedom (9 fermions * 2 spin * 2 particle/anti)

N_species_eff = 9 * 2 * 2  # 9 species, 2 spins, particle+antiparticle = 36
sum_m_galois_eV = sum_mfp * 1e6  # eV

T_reh_required = M_Pl * (target * m_proton * eta_B_exp /
                  (sum_m_galois_eV * N_species_eff * C_yield))**(1/3)

print(f"  N_species_eff = {N_species_eff} (9 Galois fermions x 2 spin x 2 part/anti)")
print(f"  sum(m_f') = {sum_m_galois_eV:.2f} eV = {sum_mfp*1e3:.3f} keV")
print(f"  Required T_reh = {T_reh_required:.3e} eV = {T_reh_required/1e9:.3e} GeV")
print(f"  T_reh / M_Pl = {T_reh_required/M_Pl:.3e}")

# Express in framework units
z_reh = math.log(T_reh_required / m_e_eV) / math.log(PHI)
print(f"  T_reh = m_e * phi^{z_reh:.2f}")

# Is z_reh a nice number?
print(f"\n  Testing if z_reh is a framework integer:")
for z in range(60, 120):
    T = m_e_eV * PHI**z
    omega = sum_m_galois_eV * N_species_eff * C_yield * (T/M_Pl)**3 / (m_proton * eta_B_exp)
    if abs(omega - target)/target < 0.05:
        print(f"    z={z}: T_reh = m_e*phi^{z} = {T:.2e} eV, Omega_ratio = {omega:.3f} "
              f"(err {abs(omega-target)/target*100:.1f}%)")


# =====================================================================
# PART 5: ALTERNATIVE — DIRECT ALGEBRAIC ARGUMENT
# =====================================================================
print("\n" + "-" * 72)
print("PART 5: Direct algebraic argument (speculative)")
print("-" * 72)

print(f"""
  Can we derive Omega_DM/Omega_b = 7-phi WITHOUT knowing T_reh?

  IDEA 1: Galois symmetry argument.
  In the spectral action, the Galois involution sigma: phi -> phi'
  maps the SM sector to the dark sector. If the spectral action
  treats both sectors equally at the fundamental level, the ratio
  of energy densities should be determined by the Galois structure.

  The SM baryon mass is dominated by QCD (m_p ~ Lambda_QCD).
  The Galois sector has NO QCD (alpha_s' < 0), so masses are bare.

  The ratio might be:
    Omega_DM/Omega_b ~ (sum m_f') / (eta_B * m_p) * (n_DM/n_baryon)

  IDEA 2: Counting argument.
  If each Galois species has the SAME abundance as baryons (from
  gravitational production at equal rate), then:
    Omega_DM/Omega_b = sum(m_f') / m_p

  sum(m_f') = {sum_m_galois_eV:.2f} eV
  m_p = {m_proton:.2f} eV
  Ratio = {sum_m_galois_eV/m_proton:.6f}  -> WAY TOO SMALL (10^-6)

  This fails because baryons have an asymmetry (eta_B ~ 6e-10)
  while DM might not.

  IDEA 3: Temperature ratio.
  The spectral action on the product M^4 x F gives TWO sectors.
  If reheating populates both equally, but the Galois sector
  decouples immediately (gravity-only), it cools independently.
  The ratio of temperatures affects abundances through:
    Omega_DM/Omega_b ~ (T_dark/T_SM)^3 * sum(m_f') / (eta_B * m_p)
""")


# =====================================================================
# PART 6: SPECTRAL ACTION APPROACH TO T_reh
# =====================================================================
print("-" * 72)
print("PART 6: Can spectral action determine T_reh?")
print("-" * 72)

# Starobinsky inflation from spectral action: S ~ integral (R + R^2/(6M^2))
# After inflation, reheating: T_reh ~ C_rh * sqrt(Gamma * M_Pl)
# where Gamma is the inflaton decay rate.
# For Starobinsky: Gamma ~ M^3 / M_Pl^2
# T_reh ~ (M^3/M_Pl)^(1/2) * some coefficient

# The scalaron mass in spectral action:
# M^2 = f_0 * c_2 / (pi^2 * f_2) where c_2, f_0, f_2 are spectral moments
# exp392 showed M cannot be determined (test function problem).

# But what if M is related to m_e and phi?
# The natural scale in the framework is m_e and M_Pl.
# We know m_e/M_Pl = alpha^(4*phi^2)

# If M_scalaron = m_e * phi^z for some z:
print(f"  If M_scalaron = m_e * phi^z:")
print(f"  T_reh ~ (M^3/M_Pl)^(1/2)")
print(f"  Required: z such that T_reh gives Omega_ratio = 7-phi")
print(f"")

# From Part 4, required T_reh ~ 1e17 GeV = 1e26 eV
# T_reh^2 ~ M^3/M_Pl
# M^3 ~ T_reh^2 * M_Pl ~ (1e26)^2 * 1.22e28 ~ 1.22e80 eV^3
# M ~ (1.22e80)^(1/3) ~ 5e26 eV ~ 5e17 GeV

# This is a very high scale. Let's see what phi exponent it corresponds to:
M_required_approx = (T_reh_required**2 * M_Pl)**(1/3)
z_M = math.log(M_required_approx / m_e_eV) / math.log(PHI)
print(f"  Required M_scalaron ~ {M_required_approx:.2e} eV = {M_required_approx/1e9:.2e} GeV")
print(f"  z_M = {z_M:.1f}")
print(f"  Closest integers: {int(z_M)}, {int(z_M)+1}")
print(f"  NOT a clean framework number.")


# =====================================================================
# PART 7: THE BARYOGENESIS CONNECTION
# =====================================================================
print("\n" + "-" * 72)
print("PART 7: Can eta_B bridge the gap?")
print("-" * 72)

# eta_B from the framework (paper eq in sec:baryo):
# eta_B ~ alpha_W^4 * J_CP * |eta_A|
# where alpha_W = alpha/sin^2(tW), J = 3.12e-5, |eta_A| = 35

alpha_W = ALPHA / sin2_tW
J_CP = 3.12e-5  # Jarlskog invariant (framework value)
eta_A = 35

eta_B_framework = alpha_W**4 * J_CP * eta_A
print(f"  Framework baryogenesis:")
print(f"    alpha_W = alpha/sin^2(tW) = {alpha_W:.6f}")
print(f"    J_CP = {J_CP:.2e}")
print(f"    |eta_A| = {eta_A}")
print(f"    eta_B ~ alpha_W^4 * J * |eta| = {eta_B_framework:.2e}")
print(f"    Experiment: {eta_B_exp:.1e}")
print(f"    Ratio: {eta_B_framework/eta_B_exp:.2f}")

print(f"""
  The framework eta_B is ~1.8x too large. But if we USE the framework
  eta_B instead of the experimental one, the required T_reh changes.
""")

T_reh_framework = M_Pl * (target * m_proton * eta_B_framework /
                  (sum_m_galois_eV * N_species_eff * C_yield))**(1/3)
z_reh_fw = math.log(T_reh_framework / m_e_eV) / math.log(PHI)
print(f"  With framework eta_B:")
print(f"    Required T_reh = {T_reh_framework:.3e} eV = {T_reh_framework/1e9:.3e} GeV")
print(f"    T_reh = m_e * phi^{z_reh_fw:.2f}")


# =====================================================================
# PART 8: THE KEY IDENTITY — WHY 7?
# =====================================================================
print("\n" + "-" * 72)
print("PART 8: Why 7? Structural meaning")
print("-" * 72)

print(f"""
  In 7 - phi, the number 7 = b1 + 1 = a1 + 2.

  Appearances of 7 in the framework:
    7 = b1 + 1 = number of vertices in affine E8 main leg + 1
    7 = a1 + 2 = diameter + 2
    7 = dim(rho_5)/dim(rho_8) + 1 = 6/2 + 1 ... no, that's 4
    7 = inert prime in Z[phi] (7 mod 5 = 2)
    7 = C(b1+1, N_gen) identity: C(7,3) = 35 = n_seesaw!
    7 = number of non-zero Galois orbit classes

  Connection to 35:
    C(7, 3) = 35 = n_seesaw
    So 7 = the number from which the seesaw exponent is derived
    combinatorially.

  The combination 7-phi = b1+phi'^2:
    - b1 = number of independent kernel generators
    - phi'^2 = Galois conjugate of phi^2
    - b1 + phi'^2 = "kernel + Galois correction"

  Speculative interpretation:
    The DM abundance is set by:
    - b1 contributions from the kernel structure (integer part)
    - phi'^2 correction from the Galois sector (irrational part)
    This gives a1+2-phi, which is as "intrinsic" as 7-phi can be.
""")

# Test: is 7-phi = something more fundamental?
# 7-phi = (a1^2+1)/(a1-1) - (a1^2+1)/(a1^2+1) * phi ... no
# 7-phi = (b1+1) - phi = (N_eig-2) - phi ... 9-2=7, yes
# 7-phi = (N_eig - 2) - phi. Why N_eig-2? Because 2 = dim of kernel Galois pair.

print(f"  Alternative: 7-phi = (N_eig - 2) - phi")
print(f"    N_eig - 2 = {N_eig-2} (total irreps minus kernel Galois pair dimension)")
print(f"    (N_eig-2) - phi = {N_eig-2-PHI:.6f} = 7-phi: {abs(N_eig-2-PHI-(7-PHI))<1e-10}")


# =====================================================================
# PART 9: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: HONEST ASSESSMENT")
print("=" * 72)

print(f"""
  RESULT OF THIS EXPERIMENT:

  DERIVED:
    - 7-phi = b1+phi'^2 = a1+2-phi has clean algebraic meaning
    - Trace(7-phi) = 13 = (a1^2+1)/2  (connects to Weinberg angle)
    - Norm(7-phi) = 41 = a1^2+(a1-1)^2  (a prime)
    - C(7,3) = 35 = n_seesaw (combinatorial connection)

  NOT DERIVED:
    - WHY Omega_DM/Omega_b equals 7-phi (no dynamical mechanism found)
    - T_reh cannot be determined from spectral action (test function problem)
    - Gravitational freeze-in requires T_reh as input
    - No spectral invariant of the Cayley graph gives 5.38 directly

  FAILED APPROACHES:
    - Mass sum ratio: sum(m_f')/m_p ~ 10^-6 (way off)
    - Equal abundance: would give ~10^-6 (no eta_B cancellation)
    - Spectral action T_reh: BLOCKED by test function problem
    - Counting argument: no simple counting gives 5.38

  CONCLUSION:
    Omega_DM/Omega_b = 7 - phi remains a **PATTERN**.
    The algebraic structure is compelling (trace=13, norm=41, C(7,3)=35)
    but the dynamical origin is UNKNOWN.

    The main blocker is the **test function problem**: the spectral action
    contains undetermined moments f_0, f_2, f_4 that set the scale of
    reheating, and no mechanism within the framework fixes these.

    This is the SAME blocker as for the scalaron mass (exp392 NEGATIVE).

  STATUS: OPEN (pattern, not derived). Honest assessment per REGULA ZERO.
""")


# =====================================================================
# PART 10: WHAT WOULD A DERIVATION LOOK LIKE?
# =====================================================================
print("-" * 72)
print("PART 10: What would close this gap?")
print("-" * 72)

print(f"""
  To derive Omega_DM/Omega_b = 7-phi, one needs EITHER:

  (A) A way to determine T_reh from spectral data.
      This requires solving the test function problem: finding f_4
      from the 600-cell geometry. Currently BLOCKED (exp392).

  (B) A T_reh-independent argument.
      For example: a topological/spectral invariant that directly
      gives the ratio of DM to baryon energy density, bypassing
      the production mechanism entirely. No such invariant is known.

  (C) A Galois symmetry argument.
      If the Galois involution exactly determines the sector ratio
      from the spectral action partition function, the ratio could
      be computed without T_reh. This requires understanding how
      the Galois automorphism acts on the PHYSICAL Hilbert space,
      which is currently OPEN.

  The most promising direction appears to be (C): the Galois
  structure already derives darkness, mass spectrum, and stability.
  If it also determines abundances, the mechanism would be purely
  algebraic (no thermal history needed).

  KEY OPEN QUESTION: Does the spectral action Z = Tr(e^(-D^2/Lambda^2))
  factorize under Galois as Z = Z_SM * Z_dark, and if so, what is
  the ratio Z_dark/Z_SM?
""")

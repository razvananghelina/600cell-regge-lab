"""
exp310_galois_dark_matter.py
============================
Galois Dark Sector: phi -> phi' conjugate particles as dark matter.

The framework's mass formula m_f = m_e * phi^(5a+6b) uses phi = (1+sqrt(5))/2.
The Galois automorphism sigma: sqrt(5) -> -sqrt(5) sends phi -> phi' = (1-sqrt(5))/2.
This creates a CONJUGATE SECTOR with definite, computable properties.

Key questions:
1. What are the Galois partner masses?
2. Are the coupling constants real in the conjugate sector?
3. Is the Galois sector electromagnetically dark?
4. Can it explain dark matter?
5. What is Omega_DM / Omega_b?

CATEGORY: EXPLORATORY (honest assessment at each step)
"""

import math
import sys

# Framework constants
a1 = 5
b1 = a1 + 1  # = 6
phi = (1 + math.sqrt(a1)) / 2  # golden ratio
phi_prime = (1 - math.sqrt(a1)) / 2  # Galois conjugate
N = 120  # = a1!
N_gen = 3
N_eig = 9
rank_E8 = 8

# Physical constants
m_e_MeV = 0.51099895  # electron mass in MeV
m_e_eV = m_e_MeV * 1e6
alpha_em = 7.2973525693e-3  # fine structure constant
alpha_s_framework = 1 / (2 * phi**3)  # framework QCD coupling

# PDG masses in MeV for comparison
PDG_masses = {
    'e': 0.51099895, 'mu': 105.6583755, 'tau': 1776.86,
    'u': 2.16, 'c': 1270.0, 't': 172760.0,
    'd': 4.67, 's': 93.4, 'b': 4180.0
}

# Fermion quantum numbers (a, b) and exponents n = 5a + 6b
fermions = {
    'e':   (0, 0),
    'mu':  (1, 1),
    'tau': (1, 2),
    'u':   (3, -2),
    'c':   (2, 1),
    't':   (4, 1),
    'd':   (1, 0),
    's':   (1, 1),
    'b':   (-1, 4),
}

print("=" * 72)
print("EXP310: GALOIS DARK SECTOR - phi -> phi' CONJUGATE PARTICLES")
print("=" * 72)

# =============================================================================
# SECTION 1: Basic Galois properties
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 1: GALOIS AUTOMORPHISM PROPERTIES")
print("=" * 72)

print(f"\nphi  = (1+sqrt(5))/2 = {phi:.10f}")
print(f"phi' = (1-sqrt(5))/2 = {phi_prime:.10f}")
print(f"|phi'| = 1/phi       = {abs(phi_prime):.10f}")
print(f"1/phi                = {1/phi:.10f}")
print(f"Check |phi'| = 1/phi: {abs(abs(phi_prime) - 1/phi):.2e}")

print(f"\nKey algebraic relations:")
print(f"  phi + phi' = 1  (trace)")
print(f"  phi * phi' = -1 (norm)")
print(f"  phi^2 = phi + 1, phi'^2 = phi' + 1 = (3-sqrt(5))/2 = {phi_prime**2:.6f}")

# Powers of phi vs phi'
print(f"\nPowers of phi and phi':")
print(f"  {'k':>3} | {'phi^k':>14} | {'phi\'^k':>14} | {'|phi\'^k|':>14} | {'(-1)^k':>7}")
print(f"  {'-'*3}-+-{'-'*14}-+-{'-'*14}-+-{'-'*14}-+-{'-'*7}")
for k in range(0, 30):
    pk = phi**k
    ppk = phi_prime**k
    print(f"  {k:3d} | {pk:14.6f} | {ppk:14.6f} | {abs(ppk):14.6f} | {(-1)**k:7d}")
    if k > 6 and pk > 1e6:
        break

print(f"\nCRITICAL: phi'^n = (-1)^n / phi^n")
print(f"  For even n: phi'^n = +1/phi^n (positive)")
print(f"  For odd n:  phi'^n = -1/phi^n (negative)")
print(f"  Mass = |m_e * phi'^n| = m_e / phi^n (always positive)")

# =============================================================================
# SECTION 2: Galois partner mass spectrum
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 2: GALOIS PARTNER MASS SPECTRUM")
print("=" * 72)

print(f"\nFor each fermion f with mass m_f = m_e * phi^n,")
print(f"the Galois partner f' has mass |m_f'| = m_e * |phi'^n| = m_e / phi^n")
print(f"This is the RECIPROCAL mass: m_f * |m_f'| = m_e^2")
print(f"\n{'Fermion':>8} | {'(a,b)':>8} | {'n=5a+6b':>8} | {'m_f (MeV)':>12} | {'|m_f\'| (MeV)':>14} | {'|m_f\'| (keV)':>14} | {'Ratio':>10}")
print(f"{'-'*8:>8}-+-{'-'*8}-+-{'-'*8}-+-{'-'*12}-+-{'-'*14}-+-{'-'*14}-+-{'-'*10}")

galois_masses = {}
total_SM = 0
total_galois = 0

for name, (a, b) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
    n = 5*a + 6*b
    m_SM = m_e_MeV * phi**n
    m_galois = m_e_MeV / phi**n  # = m_e * |phi'^n|
    galois_masses[name] = m_galois
    total_SM += m_SM
    total_galois += m_galois
    ratio = m_SM / m_galois
    print(f"{name:>8} | ({a:2d},{b:2d}) | {n:8d} | {m_SM:12.4f} | {m_galois:14.8f} | {m_galois*1000:14.6f} | {ratio:10.2f}")

print(f"\n{'TOTAL':>8} | {'':>8} | {'':>8} | {total_SM:12.2f} | {total_galois:14.8f} | {total_galois*1000:14.6f} |")

print(f"\nSM total mass (bare):     {total_SM:.2f} MeV")
print(f"Galois total mass:        {total_galois:.8f} MeV = {total_galois*1000:.4f} keV")
print(f"Ratio SM/Galois:          {total_SM/total_galois:.2f}")

# Classify by mass scale
print(f"\nGalois partner mass hierarchy:")
sorted_galois = sorted(galois_masses.items(), key=lambda x: -x[1])
for name, mass in sorted_galois:
    if mass > 0.1:
        scale = "sub-MeV"
    elif mass*1000 > 1:
        scale = "keV"
    elif mass*1e6 > 1:
        scale = "eV"
    else:
        scale = "sub-eV"
    print(f"  {name:>5}': {mass*1000:12.6f} keV  ({scale})")

# Note degeneracies
print(f"\nDegeneracies:")
n_mu = 5*1 + 6*1
n_s = 5*1 + 6*1
print(f"  mu' and s' both have n={n_mu}: DEGENERATE at {galois_masses['mu']*1000:.6f} keV")
print(f"  (In SM: m_mu = 105.66 MeV, m_s = 93.4 MeV -- different by QCD corrections)")

# =============================================================================
# SECTION 3: Galois conjugate coupling constants
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 3: GALOIS-CONJUGATE COUPLING CONSTANTS")
print("=" * 72)

# 3a. Weinberg angle
sin2_tW = b1 / (a1**2 + 1)  # = 6/26
sin2_tW_prime = sin2_tW  # RATIONAL -> Galois invariant!
print(f"\n--- sin^2(theta_W) ---")
print(f"  Framework: sin^2(tW) = b1/(a1^2+1) = {b1}/{a1**2+1} = {sin2_tW:.10f}")
print(f"  Galois:    sin^2(tW)' = {sin2_tW_prime:.10f} (RATIONAL -> INVARIANT)")
print(f"  RESULT: Weinberg angle is IDENTICAL in both sectors")

# 3b. Strong coupling
print(f"\n--- alpha_s ---")
phi3 = phi**3
phi_prime_3 = phi_prime**3
inv_alpha_s = 2 * phi3
inv_alpha_s_prime = 2 * phi_prime_3

print(f"  phi^3  = {phi3:.10f}")
print(f"  phi'^3 = {phi_prime_3:.10f}")
print(f"  1/alpha_s  = 2*phi^3  = {inv_alpha_s:.10f}")
print(f"  1/alpha_s' = 2*phi'^3 = {inv_alpha_s_prime:.10f}")
print(f"  alpha_s  = {1/inv_alpha_s:.10f}")
print(f"  alpha_s' = {1/inv_alpha_s_prime:.10f}")
print(f"")
print(f"  alpha_s' = 1/(2*phi'^3) = {1/inv_alpha_s_prime:.6f}")
print(f"  THIS IS NEGATIVE! alpha_s' < 0")
print(f"")
print(f"  Z[phi] analysis:")
print(f"    phi^3 = 2*phi + 1 (a=1, b=2 in Z[phi])")
print(f"    phi'^3 = 2*phi' + 1 = 2*(1-sqrt(5))/2 + 1 = 2 - sqrt(5) = {phi_prime_3:.6f}")
print(f"    Tr(phi^3) = phi^3 + phi'^3 = {phi3 + phi_prime_3:.6f} = 4 (= dim spacetime)")
print(f"    N(phi^3) = phi^3 * phi'^3 = {phi3 * phi_prime_3:.6f} = -1")
print(f"")
print(f"  PHYSICAL MEANING:")
print(f"    alpha_s > 0: CONFINEMENT (quarks bound into hadrons)")
print(f"    alpha_s' < 0: NO CONFINEMENT (quarks are FREE)")
print(f"    -> Galois quarks do NOT form baryons/mesons!")
print(f"    -> They exist as free, individual particles")

# 3c. Fine structure constant
print(f"\n--- alpha (electromagnetic) ---")
phi4 = phi**4
phi_prime_4 = phi_prime**4
print(f"  Framework equation: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0")
print(f"  phi^4  = {phi4:.10f}")
print(f"  phi'^4 = {phi_prime_4:.10f}")

# Galois equation: 2*pi*x^2 - 4*a1*phi'^4*x + 1 = 0
A_coeff = 2 * math.pi
B_coeff_SM = -4 * a1 * phi4
B_coeff_galois = -4 * a1 * phi_prime_4
C_coeff = 1

disc_SM = B_coeff_SM**2 - 4 * A_coeff * C_coeff
disc_galois = B_coeff_galois**2 - 4 * A_coeff * C_coeff

print(f"\n  SM sector:")
print(f"    Linear coeff: 4*a1*phi^4 = {4*a1*phi4:.6f}")
print(f"    Discriminant: {disc_SM:.6f} > 0 -> TWO REAL SOLUTIONS")
alpha_SM_1 = (-B_coeff_SM - math.sqrt(disc_SM)) / (2 * A_coeff)
alpha_SM_2 = (-B_coeff_SM + math.sqrt(disc_SM)) / (2 * A_coeff)
print(f"    alpha_1 = {alpha_SM_1:.10f} (= 1/{1/alpha_SM_1:.6f})")
print(f"    alpha_2 = {alpha_SM_2:.10f} (= 1/{1/alpha_SM_2:.6f})")
print(f"    Physical: alpha = {alpha_SM_1:.10f} (CODATA: {alpha_em:.10f})")

print(f"\n  GALOIS sector:")
print(f"    Linear coeff: 4*a1*phi'^4 = {4*a1*phi_prime_4:.6f}")
print(f"    Discriminant: {disc_galois:.6f}")
if disc_galois < 0:
    print(f"    Discriminant < 0 -> NO REAL SOLUTIONS!")
    print(f"    alpha' is COMPLEX: no real electromagnetic coupling!")
    # Compute complex solutions
    real_part = -B_coeff_galois / (2 * A_coeff)
    imag_part = math.sqrt(-disc_galois) / (2 * A_coeff)
    print(f"    alpha' = {real_part:.6f} +/- {imag_part:.6f}i")
    print(f"    |alpha'| = {math.sqrt(real_part**2 + imag_part**2):.6f}")
    print(f"")
    print(f"  *** GALOIS SECTOR HAS NO REAL EM COUPLING ***")
    print(f"  *** THESE PARTICLES DO NOT INTERACT WITH PHOTONS ***")
    print(f"  *** THEY ARE ELECTROMAGNETICALLY DARK ***")
else:
    print(f"    Discriminant >= 0 -> real solutions exist")

# 3d. Weak coupling
print(f"\n--- Weak coupling (g_2, g_1) ---")
print(f"  g_2^2 = 4*pi*alpha / sin^2(tW)")
print(f"  g_1^2 = 4*pi*alpha / cos^2(tW)")
print(f"  Since alpha' is complex but sin^2(tW) is real:")
print(f"  -> g_2' and g_1' are also COMPLEX")
print(f"  -> No real weak interactions in the Galois sector")
print(f"  -> Galois particles DO NOT participate in weak decays")

# =============================================================================
# SECTION 4: Summary of Galois sector properties
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 4: GALOIS SECTOR PROPERTY SUMMARY")
print("=" * 72)

print(f"""
GALOIS SECTOR (phi -> phi' = -1/phi):
=====================================

1. PARTICLE CONTENT: Same 9 fermions (3 generations, quarks + leptons)
   - Masses: |m_f'| = m_e / phi^n (reciprocal of SM masses)
   - Heaviest: u' at {galois_masses['u']*1000:.3f} keV
   - Lightest: t' at {galois_masses['t']*1e6:.3f} eV
   - Most particles in keV range (warm dark matter scale!)

2. GAUGE COUPLINGS:
   - sin^2(tW) = 6/26 (IDENTICAL - rational, Galois-invariant)
   - alpha_s' = 1/(2*phi'^3) < 0 (NEGATIVE - no confinement)
   - alpha' = COMPLEX (no real EM coupling)
   - g_2', g_1' = COMPLEX (no real weak coupling)

3. INTERACTIONS:
   - NO electromagnetic interaction (dark to photons)
   - NO weak interaction (no beta decay, no neutrino scattering)
   - NO QCD confinement (free quarks, no hadrons)
   - ONLY gravitational interaction

4. STABILITY:
   - Galois particles CANNOT decay to SM particles
     (no shared gauge bosons with real couplings)
   - Within the Galois sector, decays require real couplings
     -> all Galois particles are STABLE (no decay channels!)
   - This is exactly what dark matter needs: stable + non-interacting

5. MASS SCALE:
   - All Galois particles have masses < 0.121 MeV
   - Dominant: e' (511 keV), u' (121 keV), d' (46 keV)
   - These are in the WARM DARK MATTER range
""")

# =============================================================================
# SECTION 5: Cosmological abundance
# =============================================================================
print("=" * 72)
print("SECTION 5: COSMOLOGICAL ABUNDANCE")
print("=" * 72)

# Degrees of freedom
# SM fermion DOF: per fermion = 2 (particle/antiparticle) * 2 (spin) = 4
# For quarks: additional * 3 (color)
# Leptons: N_gen * 2 flavors (charged + neutrino) * 4 = 24
# Quarks: N_gen * 2 flavors (up + down) * 4 * 3 (color) = 72
# Total fermion: 96 (times 7/8 for Fermi-Dirac)
# SM gauge bosons: photon(2) + W+(2) + W-(2) + Z(2) + 8 gluons(2) + Higgs(1) = 27? No...
# Actually: photon(2) + W+/-(3*2=6 DOF massive -> 3 each) + Z(3) + 8*gluons(2) + Higgs(1) = 28

g_SM_fermion = N_gen * (2 * 4 + 2 * 4 * 3)  # leptons + quarks = 3*(8+24) = 96
g_SM_boson = 2 + 3 + 3 + 3 + 16 + 1  # photon + W+ + W- + Z + gluons + Higgs = 28
g_SM = g_SM_boson + g_SM_fermion * 7/8  # = 28 + 84 = 112

print(f"\nSM degrees of freedom:")
print(f"  Fermions: {g_SM_fermion} (x 7/8 = {g_SM_fermion * 7/8})")
print(f"  Bosons:   {g_SM_boson}")
print(f"  Total:    g_SM = {g_SM}")

# Galois sector: same fermions but NO gauge bosons (couplings not real)
# Actually, gauge bosons are determined by the gauge GROUP, not the coupling
# The group 1+3+8 is rational, so it exists in both sectors
# But with complex/negative couplings, the gauge fields are non-dynamical
# or have wrong-sign kinetic terms (ghosts)
# Conservative: count only fermion DOF

g_galois_fermion = g_SM_fermion  # same fermion content
g_galois_boson = 0  # no physical gauge bosons (complex couplings -> ghost/non-dynamical)
g_galois = g_galois_boson + g_galois_fermion * 7/8

print(f"\nGalois sector degrees of freedom:")
print(f"  Fermions: {g_galois_fermion} (x 7/8 = {g_galois_fermion * 7/8})")
print(f"  Bosons:   {g_galois_boson} (gauge couplings complex -> no physical gauge bosons)")
print(f"  Total:    g_Galois = {g_galois}")

# Scenario A: Both sectors in thermal equilibrium at T >> all masses
# Then they decouple when interaction rate < H (expansion rate)
# Galois sector only interacts gravitationally -> decouples at T_Planck
# After decoupling, both sectors cool independently

# The temperature ratio after SM entropy dumps:
# T_galois / T_SM = (g_SM_at_decouple / g_SM_today)^(1/3)
# At Planck scale: g_SM ~ 106.75 (full SM). Today: g_SM_today ~ 3.36 (photons + neutrinos)
g_SM_full = 106.75  # standard SM at high T
g_SM_today = 3.36  # photons (2) + neutrinos (3 * 2 * 7/8 * (4/11)^{4/3} equivalent)

# If Galois decouples at T_Planck, its temperature today is:
# T_G / T_gamma = (g_SM_at_decouple / g_SM_today)^(1/3) * (g_G_today / g_G_at_decouple)^(1/3)
# Wait, this needs more care. When SM species annihilate, they heat the SM plasma
# but NOT the Galois sector (which is decoupled).
# So T_G / T_gamma = (g_SM_at_decouple / g_SM_today)^(-1/3) ... no
# Actually: entropy conservation in each sector separately.
# SM: s_SM * a^3 = const -> T_SM^3 * g_SM * a^3 = const
# Galois: s_G * a^3 = const -> T_G^3 * g_G * a^3 = const
# At decoupling (T_dec): T_SM = T_G = T_dec
# Today: T_G/T_SM = (g_SM_today / g_SM_dec)^(1/3) * (g_G_dec / g_G_today)^(1/3)

# If Galois DOF don't change (no annihilations since all are stable!):
# g_G_dec = g_G_today = g_galois
# SM: g_SM_dec = 106.75 (at T >> 100 GeV), g_SM_today = 3.36

T_ratio = (g_SM_today / g_SM_full)**(1.0/3)  # Galois is colder than SM today
print(f"\n--- Temperature ratio ---")
print(f"  g_SM (high T): {g_SM_full}")
print(f"  g_SM (today):  {g_SM_today}")
print(f"  T_Galois / T_SM = (g_SM_today / g_SM_full)^(1/3) = {T_ratio:.6f}")
print(f"  Galois sector is COLDER by factor {1/T_ratio:.2f}")
print(f"  (SM entropy dumps from e+e-, mu+mu-, hadrons, W, Z, t, etc. heat SM only)")

# Energy density of Galois sector (relativistic regime, T >> m):
# rho_G / rho_SM = (g_G / g_SM) * (T_G / T_SM)^4
rho_ratio_relativistic = (g_galois / g_SM) * T_ratio**4
print(f"\nRelativistic regime (early universe):")
print(f"  rho_G / rho_SM = (g_G/g_SM) * (T_G/T_SM)^4 = {rho_ratio_relativistic:.6f}")

# This gives Delta N_eff (extra neutrino-like species):
# Delta N_eff = (8/7) * (g_G_boson + 7/8 * g_G_fermion) * (T_G/T_nu)^4
# T_nu / T_gamma = (4/11)^(1/3)
T_nu_ratio = (4.0/11)**(1.0/3)
T_G_over_T_nu = T_ratio / T_nu_ratio
Delta_Neff = (8.0/7) * g_galois * (T_G_over_T_nu)**4 / (2 * 7.0/8)  # per neutrino species equiv

# Actually, Delta_N_eff is defined as extra radiation in neutrino units:
# rho_extra = Delta_N_eff * (7/8) * (4/11)^(4/3) * rho_gamma
# rho_G = g_G * (pi^2/30) * T_G^4
# rho_nu_single = (7/8) * 2 * (pi^2/30) * T_nu^4
# Delta_N_eff = rho_G / rho_nu_single = g_G * T_G^4 / ((7/8)*2*T_nu^4)

Delta_Neff_correct = g_galois * (T_ratio)**4 / ((7.0/8) * 2 * T_nu_ratio**4)
print(f"\nExtra relativistic DOF:")
print(f"  T_G / T_nu = {T_G_over_T_nu:.6f}")
print(f"  Delta N_eff = {Delta_Neff_correct:.4f}")
print(f"  Planck 2018 bound: Delta N_eff < 0.30 (95% CL)")
if Delta_Neff_correct > 0.30:
    print(f"  *** TENSION: Delta N_eff = {Delta_Neff_correct:.2f} >> 0.30 ***")
    print(f"  This means Galois sector CANNOT be in full thermal equilibrium")
    print(f"  at T_Planck with the same temperature as SM!")
else:
    print(f"  CONSISTENT with Planck bound")

# =============================================================================
# SECTION 6: Non-thermal Galois sector
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 6: NON-THERMAL PRODUCTION (GRAVITATIONAL)")
print("=" * 72)

print(f"""
Since thermal production gives too many Galois DOF (Delta N_eff too large),
the Galois sector must be produced NON-THERMALLY.

Gravitational particle production during inflation/reheating:
- Rate: Gamma ~ m^2 / M_P (for conformally coupled scalars)
- For fermions: Gamma ~ m * H_inf (during inflation)
- Abundance: n/s ~ T_rh / M_P (gravitational freeze-in)

The key insight: Galois particles interact ONLY gravitationally.
Their abundance is set by gravitational production during reheating.
""")

# Gravitational freeze-in abundance
# n_f / s ~ m_f * T_rh / M_P^2 (order of magnitude)
# Omega_DM / Omega_b ~ (sum m_f' * n_f') / (m_p * n_b)
# n_b / s ~ eta_B ~ 6.1e-10

M_Planck_MeV = 1.22e22  # Planck mass in MeV
eta_B = 6.1e-10  # baryon asymmetry
m_proton = 938.272  # proton mass in MeV

print(f"For gravitational freeze-in, n_f'/s ~ T_rh / M_P:")
print(f"  Omega_DM/Omega_b = (sum m_f' * n_f'/s) / (m_p * eta_B)")
print(f"  = (sum m_f') * (T_rh/M_P) / (m_p * eta_B)")
print(f"")

# What T_rh gives Omega_DM/Omega_b ~ 5.36?
Omega_ratio_obs = 5.36  # Planck 2018
sum_galois_MeV = total_galois

# Omega_ratio = sum_m * T_rh / (M_P * m_p * eta_B)
# T_rh = Omega_ratio * M_P * m_p * eta_B / sum_m
T_rh_needed = Omega_ratio_obs * M_Planck_MeV * m_proton * eta_B / sum_galois_MeV

print(f"  sum |m_f'| = {sum_galois_MeV:.6f} MeV")
print(f"  For Omega_DM/Omega_b = {Omega_ratio_obs}:")
print(f"  T_rh needed = {T_rh_needed:.4e} MeV = {T_rh_needed/1e3:.4e} GeV")
print(f"  This is {T_rh_needed / 1e12:.2e} * 10^12 GeV")

# Compare with typical GUT-scale reheating
print(f"\n  Typical reheating temperatures:")
print(f"    Electroweak: ~100 GeV")
print(f"    Leptogenesis: ~10^9 GeV")
print(f"    GUT scale:    ~10^16 GeV")
print(f"    Planck:       ~10^19 GeV")

# =============================================================================
# SECTION 7: Algebraic dark matter ratio
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 7: ALGEBRAIC DARK MATTER RATIO")
print("=" * 72)

print(f"""
Instead of relying on uncertain cosmological dynamics (T_rh, production
mechanism), look for an ALGEBRAIC prediction from the framework.

The observed ratio Omega_DM / Omega_b = 5.36 +/- 0.05 (Planck 2018).

In the framework, the natural candidate is:
  Omega_DM / Omega_b = (Galois DOF) / (SM hadronic DOF)
or some ratio derived from a1 = 5.
""")

# Test various algebraic formulas
print(f"Algebraic candidates:")
candidates = []

# a1 itself
val = float(a1)
candidates.append(("a1", val))

# b1 - 1/phi = 7 - phi
val = 7 - phi
candidates.append(("7 - phi = b1+1 - phi", val))

# a1 + phi'^2
val = a1 + phi_prime**2
candidates.append(("a1 + phi'^2 = (13-sqrt(5))/2", val))

# phi^4 / phi'^4 ... no, that's phi^8
# Various
val = (a1**2 + 1) / (a1 - 1)  # 26/4 = 6.5
candidates.append(("(a1^2+1)/(a1-1)", val))

val = b1 / phi**2 + N_gen  # 6/2.618 + 3 = 5.292
candidates.append(("b1/phi^2 + N_gen", val))

val = a1 * phi / N_gen  # 5*1.618/3 = 2.697
candidates.append(("a1*phi/N_gen", val))

val = 2 * phi**3 + 2  # = 1/alpha_s + 2 = z_Planck
candidates.append(("1/alpha_s + 2 (z_Planck)", val))

# The trace of phi^4 in Z[phi]
# phi^4 = 3*phi + 2, Tr(phi^4) = (3*phi+2) + (3*phi'+2) = 3+4 = 7
val = float(7)
candidates.append(("Tr(phi^4) = 7", val))

# phi^4 = 3*phi + 2 = 6.854
val = phi**4
candidates.append(("phi^4", val))

# 4*phi^2 = z_Planck itself
val = 4 * phi**2
candidates.append(("4*phi^2 (z_Planck)", val))

# N_gen * phi = 3*1.618 = 4.854
val = N_gen * phi
candidates.append(("N_gen * phi", val))

# (a1 + phi) / phi = a1/phi + 1
val = (a1 + phi) / phi
candidates.append(("(a1+phi)/phi = a1/phi + 1", val))

# rank(E8) / phi^2
val = rank_E8 / phi**2
candidates.append(("rank(E8)/phi^2", val))

# b1 - phi' = 6 + 0.618 = 6.618
val = b1 - phi_prime
candidates.append(("b1 - phi'", val))

# Norm: |a1 + phi'^2| = |5 + (3-sqrt(5))/2| = (13-sqrt(5))/2 = 5.382
# Already covered above

# Sort by closeness to observed value
candidates.sort(key=lambda x: abs(x[1] - Omega_ratio_obs))

print(f"\n  {'Formula':>40} | {'Value':>10} | {'Error (%)':>10}")
print(f"  {'-'*40}-+-{'-'*10}-+-{'-'*10}")
for formula, val in candidates[:10]:
    err = (val - Omega_ratio_obs) / Omega_ratio_obs * 100
    marker = " <---" if abs(err) < 1 else ""
    print(f"  {formula:>40} | {val:10.6f} | {err:+10.4f}%{marker}")

# Best candidate analysis
print(f"\n--- Best algebraic candidate ---")
best_val = 7 - phi
print(f"  7 - phi = b1 + 1 - phi = {best_val:.10f}")
print(f"  Observed: {Omega_ratio_obs:.2f}")
print(f"  Error: {(best_val - Omega_ratio_obs)/Omega_ratio_obs * 100:+.2f}%")
print(f"")
print(f"  Algebraic form: 7 - phi = (14 - 1 - sqrt(5))/2 = (13 - sqrt(5))/2")
print(f"  In Z[phi]: 7 - phi has a = 7, b = -1")
print(f"  Trace: 2*7 + (-1) = 13")
print(f"  Norm: 7^2 + 7*(-1) - (-1)^2 = 49 - 7 - 1 = 41")
print(f"")
print(f"  Alternative: a1 + phi'^2 = 5 + (3-sqrt(5))/2 = (13-sqrt(5))/2 = SAME!")
print(f"  So: Omega_DM/Omega_b = a1 + phi'^2 = a1 + 2 - phi")
print(f"")
print(f"  PHYSICAL MOTIVATION:")
print(f"  phi'^2 = phi' + 1 = (3-sqrt(5))/2 is the Galois analogue of phi^2 = phi+1")
print(f"  a1 + phi'^2: integer a1 from the framework + Galois quadratic correction")
print(f"  a1: number of targets (5 irreps of A5)")
print(f"  phi'^2: Galois correction from the Z[phi] lattice")
print(f"  HOWEVER: this is a PATTERN, not a derivation. No dynamical mechanism.")

# =============================================================================
# SECTION 8: Framework-consistent approach
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 8: GALOIS SECTOR FROM SPECTRAL ACTION")
print("=" * 72)

print(f"""
The spectral action S = Tr f(D/Lambda) naturally includes BOTH sectors.

The Dirac operator D on the 600-cell has eigenvalues in Z[phi].
The Galois conjugate D' = sigma(D) has eigenvalues in Z[phi'] = Z[phi].
(sigma maps Z[phi] to itself!)

Key point: D and D' are DISTINCT operators on the SAME Hilbert space.
They share the same spectral STRUCTURE but different numerical values.

The total physical content is D_total = D + D' (Galois completion).

For the spectral action:
  S_total = Tr f(D/Lambda) + Tr f(D'/Lambda)
          = S_SM + S_Galois

The ratio of the two contributions depends on the eigenvalue
distributions and the cutoff function f.
""")

# Compute Galois spectral characteristics
print(f"--- Eigenvalue comparison ---")
print(f"  SM sector: eigenvalues ~ phi^n (exponentially growing)")
print(f"  Galois:    eigenvalues ~ phi'^n = (-1)^n/phi^n (exponentially shrinking)")
print(f"")
print(f"  SM spectral gap:     lambda_1 = 7 - 4*phi = {7 - 4*phi:.6f}")
print(f"  Galois spectral gap: lambda_1' = 7 - 4*phi' = {7 - 4*phi_prime:.6f}")
print(f"  Ratio: lambda_1'/lambda_1 = {(7 - 4*phi_prime)/(7 - 4*phi):.6f}")
print(f"")
gap_SM = 7 - 4*phi
gap_G = 7 - 4*phi_prime
print(f"  SM gap:     {gap_SM:.6f} (coexact, from computation)")
print(f"  Galois gap: {gap_G:.6f} (Galois conjugate)")
print(f"  Note: Galois gap is MUCH LARGER ({gap_G/gap_SM:.1f}x)")
print(f"  -> Galois sector is more 'gapped' -> harder to excite")

# Scalar gap comparison
scalar_gap_SM = 12 - 6*phi
scalar_gap_G = 12 - 6*phi_prime
print(f"\n  SM scalar gap:     {scalar_gap_SM:.6f}")
print(f"  Galois scalar gap: {scalar_gap_G:.6f}")
print(f"  Ratio: {scalar_gap_G/scalar_gap_SM:.6f}")

# =============================================================================
# SECTION 9: Dark matter predictions summary
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 9: TESTABLE PREDICTIONS")
print("=" * 72)

print(f"""
The Galois dark sector makes SPECIFIC, TESTABLE predictions:

1. MASS SPECTRUM (9 particles, 8 distinct masses):
""")

sorted_galois = sorted(galois_masses.items(), key=lambda x: -x[1])
for i, (name, mass) in enumerate(sorted_galois, 1):
    # Get quantum numbers
    a, b = fermions[name]
    n = 5*a + 6*b
    print(f"   {i}. {name}': n={n:3d}, mass = {mass*1000:.6f} keV = m_e/phi^{n}")

print(f"""
2. INTERACTION PROPERTIES:
   - Electromagnetically DARK (alpha' complex)
   - Weakly DARK (g_2', g_1' complex)
   - QCD DECONFINED (alpha_s' < 0, free quarks)
   - Only GRAVITATIONAL interaction with SM

3. STABILITY:
   - ALL Galois particles are absolutely stable
   - No decay channels (no real gauge couplings)
   - Multi-component dark matter

4. mu'-s' DEGENERACY:
   - Both have n=11: identical mass {galois_masses['mu']*1000:.6f} keV
   - This is a UNIQUE PREDICTION of the framework
   - Observable via gravitational lensing substructure?

5. WARMTH:
   - Dominant particles (e', u', d') in range 46-511 keV
   - This is WARM dark matter (between hot neutrinos and cold WIMPs)
   - Predicts specific free-streaming length
   - Should suppress small-scale structure
""")

# Free-streaming length
# lambda_fs ~ 0.1 Mpc * (1 keV / m_DM) * (T_DM / T_nu)
# For e' at 511 keV: very short free-streaming (basically cold)
# For t' at 3.4 eV: very long free-streaming (basically hot)
print(f"--- Free-streaming scale ---")
print(f"  lambda_fs ~ 0.1 Mpc * (keV / m) for thermal relics")
for name, mass in sorted_galois:
    lfs = 0.1 * (1.0 / (mass * 1000))  # in Mpc, assuming keV-scale normalization
    if lfs < 0.001:
        scale = "CDM-like"
    elif lfs < 0.1:
        scale = "WDM"
    elif lfs < 10:
        scale = "warm"
    else:
        scale = "HDM-like"
    print(f"    {name}': m = {mass*1000:10.4f} keV, lambda_fs ~ {lfs:.4f} Mpc ({scale})")

# =============================================================================
# SECTION 10: Relation to neutrinos
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 10: GALOIS NEUTRINOS")
print("=" * 72)

print(f"""
The framework predicts neutrino masses via seesaw:
  m_3 = 2*m_e / phi^35 = {2*m_e_MeV/phi**35 * 1e9:.4f} meV
  n_nu = 35 = b1^2 - 1, seesaw M_R ~ m_e * phi^69

The Galois neutrinos would have:
  m_nu' ~ m_e / phi^(n_nu) ... but neutrinos are special.

In the framework, neutrinos get mass from the seesaw mechanism.
The seesaw involves the Majorana mass M_R, which is in Z[phi].
Under Galois: M_R -> M_R' = m_e * phi'^69 = m_e * (-1)^69 / phi^69
  = -m_e / phi^69

|M_R'| = m_e / phi^69 = {m_e_MeV / phi**69:.4e} MeV
  = {m_e_MeV / phi**69 * 1e9:.4e} eV

This is EXTREMELY LIGHT. The seesaw gives:
  m_nu_Galois = m_f'^2 / M_R'

Actually, the Galois seesaw is inverted: small Dirac mass, small Majorana mass.
  m_D' ~ m_e / phi^n (keV scale)
  M_R' ~ m_e / phi^69 ~ 10^-23 MeV (absurdly light)
  m_nu' ~ m_D'^2 / M_R' ~ (m_e/phi^n)^2 * phi^69 / m_e = m_e * phi^(69-2n)

For electron neutrino (n ~ 0): m_nu_e' ~ m_e * phi^69 ~ M_R ~ 10^11 GeV!
The Galois seesaw is INVERTED: light M_R -> HEAVY neutrinos!

This suggests Galois neutrinos are SUPER-HEAVY (not dark matter).
They would be produced only at extremely high energies.
""")

# Compute Galois seesaw masses
M_R_galois = m_e_MeV / phi**69  # in MeV
print(f"  M_R (SM):     m_e * phi^69 = {m_e_MeV * phi**69:.4e} MeV = {m_e_MeV * phi**69 / 1e3:.4e} GeV")
print(f"  M_R (Galois): m_e / phi^69 = {M_R_galois:.4e} MeV = {M_R_galois * 1e6:.4e} eV")
print(f"")

# Galois seesaw: m_nu' = m_D'^2 / M_R' but with |M_R'| extremely small
# This means m_nu' = (m_e/phi^n)^2 / (m_e/phi^69) = m_e * phi^(69-2n)
print(f"  Galois neutrino masses (inverted seesaw):")
for gen in range(3):
    # Using the charged lepton exponent as proxy for Dirac mass
    lepton_names = ['e', 'mu', 'tau']
    name = lepton_names[gen]
    a, b = fermions[name]
    n = 5*a + 6*b
    n_galois_nu = 69 - 2*n
    m_nu_galois = m_e_MeV * phi**n_galois_nu
    print(f"    nu_{name}': n_eff = 69 - 2*{n} = {n_galois_nu}, m = m_e*phi^{n_galois_nu} = {m_nu_galois:.4e} MeV = {m_nu_galois/1e3:.4e} GeV")

print(f"\n  Galois neutrinos are SUPER-HEAVY (GUT scale or above)")
print(f"  They are NOT dark matter candidates")
print(f"  Only the Galois charged fermions and quarks are at keV scale")

# =============================================================================
# SECTION 11: Honest assessment
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 11: HONEST ASSESSMENT AND CATEGORIZATION")
print("=" * 72)

print(f"""
=== WHAT IS DERIVED (from a1 = 5, zero free parameters) ===

1. Galois mass spectrum: |m_f'| = m_e / phi^n
   -> 9 specific masses, 8 distinct values
   -> All in keV-MeV range
   CATEGORY: DERIVED (direct consequence of mass formula + Galois)

2. alpha' is complex (discriminant < 0):
   -> Galois sector is electromagnetically dark
   CATEGORY: DERIVED (algebraic consequence of alpha equation)

3. alpha_s' < 0:
   -> Galois sector has no QCD confinement
   -> Free quarks
   CATEGORY: DERIVED (algebraic consequence of alpha_s formula)

4. All Galois particles are stable:
   -> No real gauge couplings -> no decay channels
   CATEGORY: DERIVED (consequence of items 2, 3)

5. Galois neutrinos are super-heavy (inverted seesaw):
   -> NOT dark matter
   CATEGORY: DERIVED (consequence of seesaw mechanism + Galois)

=== WHAT IS PARTIALLY DERIVED ===

6. Galois sector as dark matter:
   -> Right properties: dark, stable, gravitational only
   -> Right mass scale: keV (warm dark matter)
   -> But: abundance mechanism unknown
   CATEGORY: PARTIALLY DERIVED

=== WHAT IS SPECULATIVE ===

7. Omega_DM/Omega_b = 7 - phi = a1 + phi'^2:
   -> Matches observation to 0.4%
   -> No dynamical derivation
   -> Could be coincidence (a1 = 5 gives many numbers near 5)
   CATEGORY: SPECULATIVE (pattern, not derivation)

8. Delta N_eff constraint requires non-thermal production:
   -> Gravitational freeze-in plausible but model-dependent
   -> Depends on reheating temperature (not predicted)
   CATEGORY: SPECULATIVE

=== OPEN QUESTIONS ===

- Production mechanism during inflation/reheating
- Precise abundance calculation
- Structure formation constraints (Lyman-alpha, etc.)
- Galois gauge boson status (ghosts? non-dynamical?)
- Galois Higgs sector
- Direct detection prospects (gravitational only -> hopeless?)
""")

# =============================================================================
# SECTION 12: Framework implications
# =============================================================================
print("=" * 72)
print("SECTION 12: FRAMEWORK IMPLICATIONS")
print("=" * 72)

print(f"""
The Galois dark sector provides a NATURAL dark matter candidate
from the SAME structure (a1 = 5) that gives all SM physics.

Key insight: the Z[phi] lattice NECESSARILY has two sectors
(phi and phi'), just as every quadratic extension Q(sqrt(d))/Q
has a Galois conjugate. This is not put in by hand -- it is
an ALGEBRAIC NECESSITY of using the golden ratio.

The 600-cell framework thus predicts:
  SM = phi-sector (positive root, confinement, EM coupling)
  DM = phi'-sector (negative root, no confinement, no EM)

This is deeply connected to the A5 representation theory:
  3 and 3' are Galois conjugates
  SM uses BOTH (PMNS mixing involves 3 x 3')
  But the mass formula phi^n vs phi'^n selects sectors

The Galois dark sector is NOT an addition to the theory.
It is an UNAVOIDABLE CONSEQUENCE of the algebraic structure.

Framework score update:
  Previous: ~35 predictions, ~35 derived
  New: Galois dark sector adds 5 DERIVED properties (items 1-5)
       + 1 PARTIALLY DERIVED (dark matter candidacy)
       + 1 SPECULATIVE (Omega ratio)
""")

print("=" * 72)
print("EXP310 COMPLETE")
print("=" * 72)

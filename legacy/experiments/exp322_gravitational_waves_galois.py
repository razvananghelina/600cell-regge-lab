"""
EXP-322: Gravitational Wave Signature from Galois Dark Sector
==============================================================

Framework: One-parameter theory (a_1 = 5)
Question:  Does the Galois dark sector produce a detectable stochastic
           gravitational wave background?

The Galois sector has particles with masses in the keV range:
  e'=511keV, u'=121keV, d'=46keV, mu'/s'=2.57keV, ...

If a first-order phase transition (PT) occurred at T ~ keV, it would
produce a stochastic GW background. We investigate:
  1. Does the framework predict such a PT? (Physics arguments)
  2. What frequency/amplitude would the GW have? (Parametric calc)
  3. Which detectors could see it? (LISA, ET, PTA comparison)

GW formulas: Caprini et al. (2016) [1512.06239], (2020) [1910.13125]

SPOILER: The Galois sector has NO gauge interactions (alpha' complex,
alpha_s' < 0) and was NEVER thermalized. Three independent arguments
show no first-order PT occurs. Even parametrically, the signal falls
in the nHz band (PTA), NOT mHz (LISA) or Hz (Einstein Telescope).

Author: Computational exploration for R.-C. Anghelina's framework
Date:   15 Feb 2026
"""

import numpy as np

print("=" * 78)
print("EXP 322: GRAVITATIONAL WAVE SIGNATURE FROM GALOIS DARK SECTOR")
print("=" * 78)

# ==============================================================
# CONSTANTS
# ==============================================================
a1 = 5
b1 = a1 + 1   # 6
PHI = (1 + np.sqrt(a1)) / 2        # 1.6180339887...
PHI_prime = -1.0 / PHI             # -0.6180339887...
N = 2 * a1 * (a1 - 1)              # 120
N_gen = 3

# Physical constants
m_e_keV = 510.99895
m_e_GeV = m_e_keV * 1e-6           # 5.11e-4 GeV
M_P_GeV = 1.22089e19               # Planck mass in GeV
M_P_red = M_P_GeV / np.sqrt(8*np.pi)  # 2.435e18 GeV

# Cosmological
T_CMB_GeV = 2.349e-13              # 2.725 K in GeV
H_100 = 3.2408e-18                 # 100 km/s/Mpc in Hz
h_hub = 0.674
H_0_Hz = h_hub * H_100             # 2.18e-18 Hz

# Conversions
GeV_to_Hz = 1.519e24               # 1 GeV = 1.519e24 Hz

# ==============================================================
# PART 1: GALOIS MASS SPECTRUM
# ==============================================================
print("\n" + "=" * 78)
print("PART 1: GALOIS MASS SPECTRUM")
print("=" * 78)

# (name, a, b) quantum numbers
particles = [
    ("e",   0,  0),
    ("mu",  1,  1),
    ("tau", 1,  2),
    ("u",   3, -2),
    ("c",   2,  1),
    ("t",   4,  1),
    ("d",   1,  0),
    ("s",   1,  1),
    ("b",  -1,  4),
]

print(f"\nMass duality: m_f' = m_e^2 / m_f = m_e / phi^n,  n = 5a + 6b")
print("-" * 65)
print(f"  {'Particle':<8} {'(a,b)':<10} {'n':<5} {'m_SM (MeV)':<14} {'m_Galois (keV)':<14}")
print("-" * 65)

galois_data = []
for name, a, b in particles:
    n = 5*a + 6*b
    m_sm_MeV = m_e_keV * 1e-3 * PHI**n   # SM mass in MeV
    m_gal_keV = m_e_keV / PHI**n if n != 0 else m_e_keV
    m_gal_GeV = m_gal_keV * 1e-6
    galois_data.append((name + "'", n, m_gal_keV, m_gal_GeV))
    print(f"  {name + chr(39):<7} ({a:>2},{b:>2})    {n:<5} {m_sm_MeV:<14.4f} {m_gal_keV:<14.4f}")

galois_data.sort(key=lambda x: -x[2])

print(f"\nNote: mu' and s' are DEGENERATE (both n=11)")
print(f"Heaviest Galois: e' = {galois_data[0][2]:.2f} keV (= m_e, by n=0)")
print(f"Lightest Galois: t' = {galois_data[-1][2]:.4f} keV = {galois_data[-1][2]*1e3:.2f} eV")
print(f"\nCandidate PT temperatures: T_* ~ mass of heaviest particles")
print(f"  e':      T_* ~ 511 keV   (epoch: e+e- annihilation / BBN)")
print(f"  u':      T_* ~ 121 keV")
print(f"  d':      T_* ~  46 keV")
print(f"  mu'/s':  T_* ~ 2.6 keV   (epoch: well after BBN)")

# ==============================================================
# PART 2: PHASE TRANSITION ANALYSIS - THREE NO-GO ARGUMENTS
# ==============================================================
print("\n" + "=" * 78)
print("PART 2: DOES THE FRAMEWORK PREDICT A FIRST-ORDER PHASE TRANSITION?")
print("=" * 78)

# Check Galois couplings
disc_sm = (4*a1*PHI**4)**2 - 8*np.pi
disc_gal = (4*a1*PHI_prime**4)**2 - 8*np.pi
alpha_s_gal = 1.0 / (2 * PHI_prime**3)

print(f"""
=== ARGUMENT 1: No gauge interactions ===

  SM discriminant (alpha equation):    {disc_sm:.4f} > 0  => alpha REAL
  Galois discriminant:                 {disc_gal:.4f} < 0  => alpha' COMPLEX
  Galois alpha_s' = 1/(2*phi'^3) =    {alpha_s_gal:.4f} < 0

  Without gauge interactions:
    - No EW symmetry breaking (no Higgs mechanism without real alpha)
    - No QCD confinement (alpha_s' < 0 => no asymptotic freedom)
    - No order parameter that can change discontinuously
  => No first-order PT from gauge sector

=== ARGUMENT 2: Never thermalized ===

  Galois particles produced by gravitational freeze-in (exp318):
    sigma(grav) ~ s / M_P^4  (graviton exchange)
    Y_DM << Y_equilibrium at ALL temperatures""")

# Quantitative: interaction rate vs Hubble
T_star_GeV = m_e_GeV              # 511 keV
g_star = 10.75                     # at T ~ 500 keV
H_at_T = 1.66 * np.sqrt(g_star) * T_star_GeV**2 / M_P_GeV
sigma_grav = T_star_GeV**2 / M_P_GeV**4   # ~ s/M_P^4 at T~m
Y_freeze_in = 1e-7                 # rough yield from exp318
s_entropy = (2*np.pi**2/45) * g_star * T_star_GeV**3
n_galois = Y_freeze_in * s_entropy
Gamma_grav = n_galois * sigma_grav  # rate (natural units)

print(f"  Quantitative check at T = 511 keV:")
print(f"    Hubble rate:     H = {H_at_T:.2e} GeV")
print(f"    Galois density:  n = Y*s = {n_galois:.2e} GeV^3")
print(f"    Grav. cross sec: sigma = {sigma_grav:.2e} GeV^-2")
print(f"    Interaction rate: Gamma = {Gamma_grav:.2e} GeV")
print(f"    Gamma / H = {Gamma_grav/H_at_T:.2e}")
print(f"    => Galois particles NEVER interact (Gamma/H ~ 10^{np.log10(max(Gamma_grav/H_at_T, 1e-300)):.0f})")

print(f"""
  Phase transitions require thermal equilibrium + interactions.
  The Galois sector has neither.
  => No PT from thermal dynamics

=== ARGUMENT 3: No symmetry breaking ===

  SM masses: Higgs VEV breaks EW symmetry dynamically at T ~ 160 GeV
  Galois masses: m_f' = m_e / phi^n from spectral geometry of 600-cell
    - Masses are GEOMETRIC (built into the polytope spectrum)
    - No temperature-dependent potential
    - Galois symmetry sigma: phi -> phi' ALREADY BROKEN in geometry
    - phi = 1.618... != phi' = -0.618... (fixed, not dynamical)
  => No symmetry restoration at high T, no PT

CONCLUSION: THREE independent arguments => NO first-order PT predicted.
""")

# ==============================================================
# PART 3: PARAMETRIC GW SPECTRUM (IF PT EXISTED)
# ==============================================================
print("=" * 78)
print("PART 3: PARAMETRIC GW SPECTRUM (hypothetical first-order PT)")
print("=" * 78)

print("""
Despite the no-go arguments, we compute the GW spectrum parametrically
to show the frequency range and amplitude. This is useful because:
  (a) Shows which detector band is relevant (spoiler: PTA, not LISA/ET)
  (b) Quantifies how extreme parameters would need to be
  (c) Provides template for comparison with real data

GW from first-order PT has three contributions:
  1. Sound waves in plasma (usually dominant)
  2. Bubble wall collisions (envelope approximation)
  3. MHD turbulence (subdominant)

Formulas: Caprini et al. (2016) JCAP 1604, 001 [1512.06239]
Parameters:
  T_*      = temperature at transition
  alpha    = vacuum energy / radiation energy (strength)
  beta/H_* = inverse duration in Hubble units
  v_w      = bubble wall velocity (we use v_w = 1, detonation limit)
""")


def gw_sound_waves(T_GeV, alpha, beta_H, g_s=10.75, v_w=1.0):
    """Peak frequency and amplitude from sound waves."""
    kappa_v = alpha / (0.73 + 0.083 * np.sqrt(alpha) + alpha)
    f_sw = 1.9e-5 * (1.0/v_w) * beta_H * (T_GeV/100.0) * (g_s/100.0)**(1.0/6)
    h2_Omega = (2.65e-6 * (1.0/beta_H)
                * (kappa_v * alpha / (1 + alpha))**2
                * (100.0/g_s)**(1.0/3) * v_w)
    return f_sw, h2_Omega


def gw_bubbles(T_GeV, alpha, beta_H, g_s=10.75, v_w=1.0):
    """Peak frequency and amplitude from bubble collisions."""
    kappa = 1.0 / (1 + 0.715*alpha) * (0.715*alpha + 4.0/27*np.sqrt(3*alpha/2))
    f_env = (1.65e-5 * (0.62/(1.8 - 0.1*v_w + v_w**2))
             * beta_H * (T_GeV/100.0) * (g_s/100.0)**(1.0/6))
    h2_Omega = (1.67e-5 * (1.0/beta_H)**2
                * (kappa * alpha / (1 + alpha))**2
                * (100.0/g_s)**(1.0/3)
                * (0.11 * v_w**3 / (0.42 + v_w**2)))
    return f_env, h2_Omega


def gw_turbulence(T_GeV, alpha, beta_H, g_s=10.75, v_w=1.0, eps=0.05):
    """Peak frequency and amplitude from MHD turbulence."""
    kappa_v = alpha / (0.73 + 0.083 * np.sqrt(alpha) + alpha)
    kappa_t = eps * kappa_v
    f_turb = 2.7e-5 * (1.0/v_w) * beta_H * (T_GeV/100.0) * (g_s/100.0)**(1.0/6)
    h2_Omega = (3.35e-4 * (1.0/beta_H)
                * (kappa_t * alpha / (1 + alpha))**(3.0/2)
                * (100.0/g_s)**(1.0/3) * v_w)
    return f_turb, h2_Omega


# Scan over parameters
T_star = m_e_GeV   # 511 keV (heaviest Galois particle)
g_s = 10.75

print(f"Reference: T_* = {T_star*1e6:.0f} keV (= m_e, heaviest Galois mass)")
print(f"           g_* = {g_s} (at T ~ 500 keV)")
print(f"           v_w = 1.0 (detonation)")

alphas = [0.01, 0.1, 0.5, 1.0]
betas = [10, 100, 1000, 10000]

print(f"\nSound waves (dominant contribution):")
print(f"  {'alpha':<7} {'beta/H*':<9} {'f_peak (Hz)':<14} {'h2*Omega':<14} {'Band':<12}")
print("-" * 60)

for alpha_val in alphas:
    for beta_val in betas:
        f_sw, Om_sw = gw_sound_waves(T_star, alpha_val, beta_val, g_s)
        if f_sw < 1e-9:
            band = "sub-PTA"
        elif f_sw < 1e-7:
            band = "PTA (nHz)"
        elif f_sw < 1e-4:
            band = "muHz gap"
        elif f_sw < 1e-1:
            band = "LISA (mHz)"
        else:
            band = "HIGH"
        print(f"  {alpha_val:<7.2f} {beta_val:<9} {f_sw:<14.2e} {Om_sw:<14.2e} {band}")
    print()

# ==============================================================
# PART 4: DETECTOR SENSITIVITY COMPARISON
# ==============================================================
print("=" * 78)
print("PART 4: DETECTOR SENSITIVITY BANDS")
print("=" * 78)

detectors = [
    ("PTA (NANOGrav)",   1e-9,  1e-7,  1e-9),
    ("mu-Ares (prop.)",  1e-7,  1e-4,  1e-15),
    ("LISA",             1e-4,  1e-1,  1e-12),
    ("DECIGO/BBO",       1e-2,  1e1,   1e-15),
    ("Einstein Tel.",    1e0,   1e4,   1e-13),
]

print(f"\n  {'Detector':<18} {'f_min (Hz)':<12} {'f_max (Hz)':<12} {'h2*Omega sens':<14}")
print("-" * 60)
for name, fmin, fmax, sens in detectors:
    print(f"  {name:<18} {fmin:<12.0e} {fmax:<12.0e} {sens:<14.0e}")

# What beta/H_* is needed to reach each band?
factor = 1.9e-5 * (T_star/100.0) * (g_s/100.0)**(1.0/6)
print(f"\n  Frequency prefactor for T_* = 511 keV: {factor:.2e} Hz * beta/H_*")
print(f"\n  beta/H_* needed to reach each detector (with v_w=1):")
for name, fmin, fmax, sens in detectors:
    f_mid = np.sqrt(fmin * fmax)
    beta_needed = f_mid / factor
    print(f"    {name:<18}: beta/H_* ~ {beta_needed:.1e}  (for f_mid = {f_mid:.1e} Hz)")

print(f"""
KEY INSIGHT:
  For T_* = 511 keV, the sound wave peak frequency is:
    f_sw = {factor:.2e} * (beta/H_*) Hz

  To reach LISA band  (f ~ 1 mHz):  need beta/H_* ~ {1e-3/factor:.0e}
  To reach ET band    (f ~ 10 Hz):  need beta/H_* ~ {10/factor:.0e}
  To reach PTA band   (f ~ 10 nHz): need beta/H_* ~ {1e-8/factor:.0e}

  => LISA and Einstein Telescope are INACCESSIBLE for keV-scale PT
  => Only PTA is potentially relevant (beta/H_* ~ {1e-8/factor:.0f})
""")

# ==============================================================
# PART 5: BEST-CASE SCENARIO vs PTA / NANOGrav
# ==============================================================
print("=" * 78)
print("PART 5: BEST-CASE SCENARIO AND NANOGrav COMPARISON")
print("=" * 78)

# NANOGrav 15yr signal
A_nano = 2.4e-15       # characteristic strain amplitude
f_yr = 1.0 / (365.25 * 24 * 3600)   # 1/year in Hz = 3.17e-8 Hz
# h2*Omega_GW(f) = (2*pi^2/3) * f^2 * h_c^2 / (H_0/h)^2
# where H_0/h = 100 km/s/Mpc = 3.24e-18 Hz
H_100_val = 3.2408e-18
h2_Omega_nano = (2*np.pi**2/3) * f_yr**2 * A_nano**2 / H_100_val**2

print(f"\nNANOGrav 15yr stochastic GW background:")
print(f"  Amplitude: A = {A_nano:.1e} at f_ref = 1/yr = {f_yr:.2e} Hz")
print(f"  h2*Omega(1/yr) = {h2_Omega_nano:.2e}")
print(f"  Interpretation: likely supermassive BH binaries, possibly cosmological")

print(f"\nBest-case Galois scenario: alpha = 1.0 (very strong PT), v_w = 1")
print(f"  T_* = 511 keV, g_* = {g_s}")
print(f"\n  {'beta/H*':<10} {'f_sw (Hz)':<14} {'h2*Omega_sw':<14} {'vs NANOGrav':<14} {'Detectable?'}")
print("-" * 70)

for beta_val in [10, 50, 100, 500, 1000, 5000, 10000]:
    f_sw, Om_sw = gw_sound_waves(T_star, 1.0, beta_val, g_s)
    ratio_nano = Om_sw / h2_Omega_nano if h2_Omega_nano > 0 else 0
    # Check PTA detectability: need f in [1e-9, 1e-7] AND Omega > 1e-9
    in_PTA = (1e-9 < f_sw < 1e-7)
    above_sens = (Om_sw > 1e-9)
    if in_PTA and above_sens:
        detect = "YES (PTA)"
    elif in_PTA:
        detect = "in band, weak"
    else:
        detect = "NO (out of band)"
    print(f"  {beta_val:<10} {f_sw:<14.2e} {Om_sw:<14.2e} {ratio_nano:<14.2e} {detect}")

# Specific case: beta/H_* that puts f_sw at 10 nHz
f_target = 1e-8   # 10 nHz
beta_for_PTA = f_target / factor
f_sw_pta, Om_sw_pta = gw_sound_waves(T_star, 1.0, beta_for_PTA, g_s)
f_env_pta, Om_env_pta = gw_bubbles(T_star, 1.0, beta_for_PTA, g_s)
f_turb_pta, Om_turb_pta = gw_turbulence(T_star, 1.0, beta_for_PTA, g_s)

print(f"\n  Tuned to NANOGrav frequency (f = 10 nHz):")
print(f"    beta/H_* = {beta_for_PTA:.0f}")
print(f"    Sound waves:  f = {f_sw_pta:.2e} Hz, h2*Omega = {Om_sw_pta:.2e}")
print(f"    Bubbles:      f = {f_env_pta:.2e} Hz, h2*Omega = {Om_env_pta:.2e}")
print(f"    Turbulence:   f = {f_turb_pta:.2e} Hz, h2*Omega = {Om_turb_pta:.2e}")
print(f"    NANOGrav:     f = {f_yr:.2e} Hz, h2*Omega = {h2_Omega_nano:.2e}")
print(f"    Total/NANOGrav = {(Om_sw_pta+Om_env_pta+Om_turb_pta)/h2_Omega_nano:.2f}")

# ==============================================================
# PART 6: SCALING ANALYSIS - WHY keV PTs ARE HARD TO DETECT
# ==============================================================
print("\n" + "=" * 78)
print("PART 6: SCALING ANALYSIS")
print("=" * 78)

# Compare with EW-scale PT (the standard benchmark)
T_ew = 100.0  # GeV
T_gal = T_star  # 5.11e-4 GeV
ratio_T = T_ew / T_gal

print(f"""
  Comparison: EW-scale PT (T=100 GeV) vs Galois PT (T=511 keV)
  Temperature ratio: T_EW / T_Galois = {ratio_T:.0f}

  At SAME alpha and beta/H_*:
    f_peak ratio:    T_EW/T_Gal = {ratio_T:.0f}x  (frequency scales with T_*)
    Omega ratio:     1           (amplitude independent of T_*)

  EW-scale PT at beta/H_*=100 gives f ~ 1 mHz (LISA band).
  Galois PT at same beta/H_* gives f ~ {ratio_T:.0f}x lower = {1e-3/ratio_T:.1e} Hz

  To push Galois signal UP to LISA band:
    Need beta/H_* ~ {ratio_T:.0f}x larger
    But h2*Omega(sound waves) ~ 1/beta
    => Amplitude {ratio_T:.0f}x WEAKER ({np.log10(ratio_T):.1f} orders of magnitude)

  This is the FUNDAMENTAL OBSTRUCTION:
    Low T_* => low frequency AND low amplitude (when adjusted to higher f)
""")

# At what T_* do we enter LISA band naturally?
# f_sw = 1.9e-5 * beta * T/100 * (g/100)^{1/6}
# For f = 1 mHz, beta = 100, g = 100:
# T/100 = 1e-3 / (1.9e-5 * 100) = 0.53 => T ~ 53 GeV
print(f"  Minimum T_* for LISA detectability (beta/H_*=100, alpha=0.3):")
T_for_LISA = 1e-3 / (1.9e-5 * 100 * (10.75/100)**(1.0/6))
print(f"    T_* > {T_for_LISA:.0f} GeV")
print(f"    Galois T_* = {T_gal*1e3:.1f} MeV = {T_for_LISA/T_gal:.0f}x too low")

print(f"\n  Minimum T_* for Einstein Telescope (beta/H_*=100):")
T_for_ET = 10.0 / (1.9e-5 * 100 * (100/100)**(1.0/6))
print(f"    T_* > {T_for_ET:.0e} GeV")
print(f"    This is {T_for_ET/T_gal:.0e}x above Galois scale")

# ==============================================================
# PART 7: SCAN OVER ALL GALOIS MASS THRESHOLDS
# ==============================================================
print("\n" + "=" * 78)
print("PART 7: ALL GALOIS MASS THRESHOLDS")
print("=" * 78)

print(f"\nBest case (alpha=1, beta/H_*=100) for each Galois mass threshold:")
print(f"  {'Particle':<8} {'T_* (keV)':<12} {'g_*':<6} {'f_sw (Hz)':<14} {'h2*Omega':<14} {'PTA?'}")
print("-" * 65)

# g_* depends on temperature
def g_star_at_T(T_keV):
    """Approximate g_* at temperature T in keV."""
    T_MeV = T_keV * 1e-3
    if T_MeV > 200:     return 10.75  # gamma + e + 3*nu
    elif T_MeV > 1:     return 10.75
    elif T_MeV > 0.5:   return 10.75  # above e+e- annihilation
    elif T_MeV > 0.1:   return 7.25   # gamma + 3*nu (after e+e- ann.)
    else:                return 3.36   # gamma + nu (decoupled)

for name, n, m_keV, m_GeV in galois_data:
    g_s_local = g_star_at_T(m_keV)
    f_sw, Om_sw = gw_sound_waves(m_GeV, 1.0, 100, g_s_local)
    in_pta = "YES" if (1e-9 < f_sw < 1e-7 and Om_sw > 1e-9) else "no"
    print(f"  {name:<8} {m_keV:<12.4f} {g_s_local:<6.2f} {f_sw:<14.2e} {Om_sw:<14.2e} {in_pta}")

print(f"\nOnly e' (511 keV) and u' (121 keV) have f_peak in or near PTA band.")
print(f"Lighter particles produce GW below any detector band.")

# ==============================================================
# PART 8: ALTERNATIVE GW SOURCES
# ==============================================================
print("\n" + "=" * 78)
print("PART 8: ALTERNATIVE GW SOURCES FROM GALOIS SECTOR")
print("=" * 78)

print("""
(a) SECOND-ORDER INDUCED GW
    Large density perturbations at second order source tensor modes.
    Requires enhanced power spectrum at k ~ 1/lambda_Galois.
    The framework does NOT predict enhanced perturbations at any scale.
    Standard nearly scale-invariant spectrum (n_s=0.965) gives negligible signal.
    => NOT PREDICTED

(b) GRAVITATIONAL PRODUCTION AT REHEATING
    SM + SM -> Galois + Galois via graviton exchange during reheating.
    Individual quantum scatterings do NOT produce classical GW.
    Collective backreaction is suppressed by Y_DM / Y_eq << 1.
    The gravitational production is smooth, no violent dynamics.
    => NEGLIGIBLE

(c) DOMAIN WALLS FROM GALOIS Z_2
    Galois group Gal(Q(phi)/Q) = Z_2 gives sigma: phi <-> phi'.
    Domain walls require DEGENERATE vacua populated in different regions.
    But the geometric structure ALREADY selects phi = 1.618...
    Both SM and Galois sectors coexist in same spacetime, not in domains.
    No domain wall formation.
    => NOT APPLICABLE

(d) COSMIC STRINGS
    Require broken U(1) symmetry.
    Galois sector has no U(1)_EM (alpha' complex, no real gauge group).
    => NOT APPLICABLE

(e) PREHEATING / PARAMETRIC RESONANCE
    Could produce GW if inflaton couples to Galois sector.
    But Galois particles couple ONLY gravitationally to everything.
    Gravitational coupling to inflaton: delta(rho)/rho ~ (H/M_P)^2 ~ 10^{-10}.
    => NEGLIGIBLE
""")

# ==============================================================
# PART 9: COMPARISON WITH KNOWN SIGNALS
# ==============================================================
print("=" * 78)
print("PART 9: NANOGrav 15yr SIGNAL - COULD IT BE GALOIS?")
print("=" * 78)

print(f"""
NANOGrav (2023) detected a stochastic GW background:
  Frequency: f ~ 2-90 nHz
  Amplitude: A = (2.4 +/- 0.7) x 10^-15 at f = 1/yr
  h2*Omega ~ {h2_Omega_nano:.1e} at 1/yr = {f_yr:.2e} Hz
  Spectral index consistent with f^(2/3) (supermassive BH binaries)

Could the Galois sector contribute?""")

# For Galois PT at T_* = 511 keV to match NANOGrav frequency:
# need beta/H_* ~ 150 (from Part 5)
# At alpha = 1 (extreme), this gives h2*Omega ~ few x 10^-9
# NANOGrav signal is h2*Omega ~ 4e-9
# So parametrically, the AMPLITUDE could match!

print(f"  IF Galois had a strong PT (alpha=1) at T_*=511 keV:")
print(f"    At beta/H_*={beta_for_PTA:.0f}: f_peak={f_sw_pta:.1e} Hz")
print(f"    Amplitude: h2*Omega = {Om_sw_pta:.1e}")
print(f"    NANOGrav:  h2*Omega = {h2_Omega_nano:.1e}")
print(f"    Ratio: {Om_sw_pta/h2_Omega_nano:.1f}")
print(f"""
  Parametrically, a STRONG Galois PT could match NANOGrav amplitude!
  But spectral shape would differ from SMBHB power law.

  HOWEVER: The framework predicts NO first-order PT (Part 2).
  => The NANOGrav signal is NOT from the Galois dark sector.
  => Most likely supermassive BH binaries (consistent with f^{2/3} spectrum).
""")

# ==============================================================
# PART 10: SPECTRAL SHAPE (for completeness)
# ==============================================================
print("=" * 78)
print("PART 10: GW SPECTRAL SHAPE (hypothetical, alpha=1, beta/H_*=150)")
print("=" * 78)

def sw_spectrum(f, f_peak, h2_Omega_peak):
    """Sound wave spectral shape (Caprini et al. 2016 Eq. 36)."""
    x = f / f_peak
    shape = x**3 * (7.0 / (4.0 + 3.0 * x**2))**(7.0/2)
    return h2_Omega_peak * shape

f_peak_ref, Om_peak_ref = gw_sound_waves(T_star, 1.0, beta_for_PTA, g_s)

# Sample frequencies from 0.1 nHz to 1 muHz
log_f = np.linspace(-10, -6, 30)
freqs = 10**log_f

print(f"\n  f_peak = {f_peak_ref:.2e} Hz, h2*Omega_peak = {Om_peak_ref:.2e}")
print(f"\n  {'f (Hz)':<14} {'h2*Omega':<14} {'Notes'}")
print("-" * 50)

for f_val in freqs:
    Om_val = sw_spectrum(f_val, f_peak_ref, Om_peak_ref)
    notes = ""
    if abs(np.log10(f_val) - np.log10(f_peak_ref)) < 0.2:
        notes = "<-- peak"
    if 1e-9 < f_val < 1e-7:
        notes += " [PTA band]"
    if Om_val > 1e-20:  # only print non-negligible
        print(f"  {f_val:<14.2e} {Om_val:<14.2e} {notes}")

# ==============================================================
# FINAL SUMMARY
# ==============================================================
print("\n" + "=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print(f"""
QUESTION: Does the Galois dark sector produce detectable GW?

ANSWER: NO. Three independent DERIVED arguments:

  1. NO GAUGE INTERACTIONS
     alpha' = complex, alpha_s' = {alpha_s_gal:.4f} < 0
     => No symmetry breaking, no confinement, no PT mechanism

  2. NEVER THERMALIZED
     Gamma/H = {Gamma_grav/H_at_T:.0e} at T = 511 keV
     => No thermal bath, no collective dynamics

  3. GEOMETRIC MASSES (no Higgs mechanism)
     m_f' = m_e/phi^n from 600-cell spectrum
     => No order parameter, no restoration at high T

FREQUENCY ANALYSIS (even if PT existed):
  For T_* = 511 keV (heaviest Galois particle):
    - LISA (mHz):        need beta/H_* ~ 10^7   => INACCESSIBLE
    - Einstein Tel (Hz): need beta/H_* ~ 10^12  => INACCESSIBLE
    - PTA (nHz):         need beta/H_* ~ 150    => marginally possible
                         but requires alpha > 0.5 (extreme)

PREDICTION:  Galois sector produces ZERO GW signal.
             LISA: NO (wrong frequency band by 5 orders of magnitude)
             ET:   NO (wrong frequency band by 10+ orders)
             PTA:  NO (right band but no PT mechanism)

CATEGORY: DERIVED (no-go result, three independent arguments)

NOTE ON DETECTOR CHOICE:
  The user suggested LISA / Einstein Telescope detection.
  For T_* ~ keV, the GW frequency is in the nHz range (PTA band).
  LISA sees mHz (PT at T ~ 10-1000 GeV).
  ET sees Hz-kHz (PT at T ~ 10^6-10^9 GeV or compact binaries).
  The correct detector for keV-scale physics would be PTA.

BROADER IMPLICATION:
  The Galois dark sector is MAXIMALLY DARK:
    - No electromagnetic coupling  (alpha' complex)
    - No strong interaction        (alpha_s' < 0)
    - No phase transitions         (no symmetry breaking)
    - No gravitational waves       (no violent dynamics)
    - ONLY gravitational coupling  (through spacetime curvature)
  This maximal darkness is not assumed but DERIVED from the framework.
""")

print("=" * 78)
print("EXP 322 COMPLETE")
print("=" * 78)

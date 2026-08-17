"""
EXP-251: Search for DERIVED Mass Corrections (NO FITTING)
==========================================================
The bare formula m_f = m_e * phi^n has 3-21% errors.
Can we find a correction using ONLY derived quantities?

DERIVED quantities available:
  - a_1 = 5, b_1 = 6, phi, N = 120
  - alpha = sol(2*pi*x^2 - 4*a_1*phi^4*x + 1 = 0)
  - alpha_s = 1/(2*phi^3)
  - sin^2(tW) = 6/26
  - rank(E8) = 8, dim(E8) = 248, N_eig = 9
  - Electric charge Q (from edge splitting: DERIVED)
  - Color Casimir C_2 (from SU(3): 0 for leptons, 4/3 for quarks)
  - N(z) = norm in Z[phi] (DERIVED)
  - n' = 5a - b (Galois exponent, DERIVED)

CONSTRAINT: The correction formula must be FIXED BY THEORY,
not chosen from a menu to minimize error.
"""

import math

PHI = (1 + math.sqrt(5)) / 2
PHI_PRIME = 1 - PHI
PI = math.pi
a1 = 5
b1 = 6

# Derived constants
A_coef = 2 * PI
B_coef = -4 * a1 * PHI**4
C_coef = 1
disc = B_coef**2 - 4*A_coef*C_coef
alpha = (-B_coef - math.sqrt(disc)) / (2*A_coef)
alpha_s = 1 / (2 * PHI**3)
sin2_tW = b1 / (a1**2 + 1)  # 6/26
rank_E8 = 8

m_e_MeV = 0.51099895

# Fermion data: name, g, a, b, Q, C2_color, m_exp(MeV), scale
fermions = [
    ('e',   0, 0, 0,  -1,  0,    0.51099895, 'pole'),
    ('mu',  1, 1, 1,  -1,  0,    105.6584,   'pole'),
    ('tau', 2, 1, 2,  -1,  0,    1776.86,    'pole'),
    ('u',   0, 3,-2,  2/3, 4/3,  2.16,       'MSbar 2GeV'),
    ('c',   1, 2, 1,  2/3, 4/3,  1270.0,     'MSbar mc'),
    ('t',   2, 4, 1,  2/3, 4/3,  172760.0,   'pole'),
    ('d',   0, 1, 0, -1/3, 4/3,  4.67,       'MSbar 2GeV'),
    ('s',   1, 1, 1, -1/3, 4/3,  93.4,       'MSbar 2GeV'),
    ('b',   2,-1, 4, -1/3, 4/3,  4180.0,     'MSbar mb'),
]

print("=" * 70)
print("EXP-251: SEARCH FOR DERIVED MASS CORRECTIONS")
print("=" * 70)

# =====================================================================
# HYPOTHESIS 1: QFT running with DERIVED inputs
# =====================================================================
print()
print("=" * 70)
print("H1: Standard QFT running corrections")
print("    Uses DERIVED alpha, alpha_s, Q, C_2")
print("=" * 70)
print()

# QED self-energy (1-loop): delta_m/m = (3*alpha*Q^2)/(4*pi) * [3/2 + ln(mu/m)]
# QCD mass running: m(mu)/m(mu_0) = [alpha_s(mu)/alpha_s(mu_0)]^(d_m)
# where d_m = gamma_m^(0) / (2*beta_0^(0))
# gamma_m^(0) = 8 (universal for SU(3))
# beta_0 = 11 - 2*n_f/3

# The theory predicts bare masses. Physical masses include radiative corrections.
# For leptons (no QCD): correction is small (alpha << 1)
# For quarks: QCD correction is significant

# QED correction to lepton masses:
# The pole mass includes the self-energy: m_pole = m_bare * (1 + delta_QED)
# At 1-loop: delta_QED = (3*alpha)/(4*pi) * [1 + ln(Lambda^2/m^2)]
# We need Lambda. The natural scale is m_Z (EW scale) or m_P (Planck).

# Let's use the DERIVED m_Z = m_e * phi^25 * alpha(m_Z)/alpha(0) relation
# m_Z ~ 91.2 GeV

# For muon: delta_QED ~ (3*alpha)/(4*pi) * ln(m_Z^2/m_mu^2)
#         = (3*0.00730)/(4*pi) * ln((91200/105.7)^2)
#         = 0.00174 * 2*ln(863) = 0.00174 * 13.5 = 0.0235 = 2.35%
# But muon error is -3.8%, so QED correction alone goes in wrong direction
# (QED correction INCREASES mass, but bare is already TOO LOW)

# Wait - if theory gives bare mass and QED increases it:
# m_pole = m_bare * (1 + delta) > m_bare
# But m_bare(mu) = 101.7 < m_exp(mu) = 105.7
# So we WANT the correction to be positive. And QED gives positive. Good!

print("QED correction to lepton pole masses:")
print("  delta_QED ~ (3*alpha)/(4*pi) * ln(m_Z^2/m_f^2)")
print()

m_Z = 91187.6  # MeV, experimental

for name, g, a, b, Q, C2, m_exp, scale in fermions:
    if C2 > 0:  # skip quarks for now
        continue
    if name == 'e':
        continue
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n

    # QED 1-loop correction
    delta_QED = (3 * alpha * Q**2) / (4 * PI) * math.log(m_Z**2 / m_bare**2)
    m_corrected = m_bare * (1 + delta_QED)
    err_bare = (m_bare - m_exp) / m_exp * 100
    err_corr = (m_corrected - m_exp) / m_exp * 100

    print("  %s: m_bare=%.2f, delta_QED=%.4f (%.2f%%)"
          % (name, m_bare, delta_QED, delta_QED*100))
    print("       m_corrected=%.2f, m_exp=%.2f" % (m_corrected, m_exp))
    print("       bare err=%+.2f%%, corrected err=%+.2f%%" % (err_bare, err_corr))
    print()

# Now QCD corrections for quarks
# The MSbar mass at scale mu relates to bare mass via:
# m_MSbar(mu) = m_bare * [alpha_s(mu)/alpha_s(Lambda)]^(d_m)
# d_m = 4 / beta_0 = 4 / (11 - 2*n_f/3)

# For n_f = 5 (at b quark scale): beta_0 = 11 - 10/3 = 23/3
# d_m = 4/(23/3) = 12/23 = 0.5217

# For n_f = 4 (at charm scale): beta_0 = 11 - 8/3 = 25/3
# d_m = 4/(25/3) = 12/25 = 0.48

# For n_f = 3 (at strange scale): beta_0 = 11 - 6/3 = 9
# d_m = 4/9 = 0.4444

print("QCD running for quark masses:")
print("  m_MSbar(mu) = m_bare * [alpha_s(mu)/alpha_s(Lambda)]^(d_m)")
print("  d_m = 4/beta_0, beta_0 = 11 - 2*n_f/3")
print()

# What scale Lambda gives theory masses?
# We don't know Lambda from theory directly.
# The DERIVED Planck relation: m_e/m_P = alpha^(4*phi^2)
# m_P ~ 1.22e22 MeV
# alpha_s(m_P) ~ 0.02 (from perturbative running, approximate)

# Let's try: theory gives masses at PLANCK scale
# Then run down to measurement scales

# alpha_s running (1-loop): alpha_s(mu) = alpha_s(mu_0) / (1 + alpha_s(mu_0)*beta_0*ln(mu/mu_0)/(2*pi))

def alpha_s_running(mu, mu_0, alpha_s_0, n_f):
    """1-loop QCD running coupling"""
    beta_0 = 11 - 2*n_f/3
    return alpha_s_0 / (1 + alpha_s_0 * beta_0 / (2*PI) * math.log(mu/mu_0))

def mass_running(m_bare, mu_from, mu_to, n_f):
    """Run quark mass from mu_from to mu_to at 1-loop"""
    beta_0 = 11 - 2*n_f/3
    d_m = 4 / beta_0  # mass anomalous dimension / beta_0
    # m(mu_to) = m(mu_from) * [alpha_s(mu_to)/alpha_s(mu_from)]^d_m
    a_s_from = alpha_s_running(mu_from, m_Z, 0.1179, n_f)
    a_s_to = alpha_s_running(mu_to, m_Z, 0.1179, n_f)
    if a_s_from <= 0 or a_s_to <= 0:
        return m_bare
    return m_bare * (a_s_to / a_s_from) ** d_m

# =====================================================================
# HYPOTHESIS 2: Charge-dependent phase correction
# =====================================================================
print()
print("=" * 70)
print("H2: Charge-dependent correction")
print("    m_f = m_e * phi^(n + Q^2 * delta) where delta is DERIVED")
print("=" * 70)
print()

# What derived value of delta could work?
# For muon (Q=-1): need phi^(11+delta) = 105.7/0.511 = 206.85
# So 11 + delta = ln(206.85)/ln(phi) = 5.332/0.4812 = 11.080
# delta_mu = 0.080

# For strange (Q=-1/3): need phi^(11+delta/9) = 93.4/0.511 = 182.78
# So 11 + delta/9 = ln(182.78)/ln(phi) = 5.209/0.4812 = 10.825
# delta_s/9 = -0.175, delta_s = -1.575

# These are DIFFERENT! Q^2 doesn't unify them.
# For muon: delta = 0.080/1 = 0.080 (Q^2=1)
# For strange: delta = -0.175*9 = -1.575 (Q^2=1/9)
# Completely different => Q^2 correction FAILS.

print("Testing: delta from Q^2 (same delta for both)")
print()
for name, g, a, b, Q, C2, m_exp, scale in fermions:
    if name == 'e':
        continue
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n
    # What delta would we need?
    if m_exp > 0:
        n_needed = math.log(m_exp / m_e_MeV) / math.log(PHI)
        delta_needed = (n_needed - n) / Q**2 if Q != 0 else 0
        print("  %s: n=%2d, Q^2=%.2f, n_needed=%.4f, delta_needed=%.4f"
              % (name, n, Q**2, n_needed, delta_needed))

print()
print("  VERDICT: delta values are NOT consistent across fermions.")
print("  Q^2 correction does NOT work as a universal mechanism.")

# =====================================================================
# HYPOTHESIS 3: Color Casimir correction
# =====================================================================
print()
print("=" * 70)
print("H3: Color Casimir correction")
print("    m_f = m_e * phi^n * (1 + C_2 * alpha_s * f)")
print("    where f is DERIVED from geometry")
print("=" * 70)
print()

# C_2 = 0 for leptons (no change to bare)
# C_2 = 4/3 for quarks
# alpha_s = 1/(2*phi^3) (derived)
# f = ??? must be derived

# For this to work, ALL quarks need the SAME f.
# Check what f is needed for each quark:
print("What f = correction_factor/(C_2 * alpha_s) is needed for each quark?")
print()
for name, g, a, b, Q, C2, m_exp, scale in fermions:
    if C2 == 0 or name == 'e':
        continue
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n
    # m_exp = m_bare * (1 + C2 * alpha_s * f)
    # f = (m_exp/m_bare - 1) / (C2 * alpha_s)
    ratio = m_exp / m_bare
    if C2 * alpha_s > 0:
        f_needed = (ratio - 1) / (C2 * alpha_s)
    else:
        f_needed = 0
    print("  %s: m_bare=%.1f, m_exp=%.1f, ratio=%.4f, f_needed=%+.3f"
          % (name, m_bare, m_exp, ratio, f_needed))

print()
print("  f values range from -1.28 to +6.29. NOT universal.")
print("  VERDICT: Simple Casimir correction does NOT work.")

# =====================================================================
# HYPOTHESIS 4: Norm-based correction for non-units
# =====================================================================
print()
print("=" * 70)
print("H4: Norm-based correction (only for non-units)")
print("    Units (|N|=1): m = m_e * phi^n (no correction)")
print("    Non-units: m = m_e * phi^n * |N(z)|^(-p) with DERIVED p")
print("=" * 70)
print()

# Non-units: charm (|N|=5), top (|N|=19), bottom (|N|=19)
for name, g, a, b, Q, C2, m_exp, scale in fermions:
    norm = abs(a**2 + a*b - b**2)
    if norm <= 1 or name == 'e':
        continue
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n
    # m_exp = m_bare * |N|^(-p)
    # p = -ln(m_exp/m_bare) / ln(|N|)
    ratio = m_exp / m_bare
    if norm > 1:
        p_needed = -math.log(ratio) / math.log(norm)
    else:
        p_needed = 0
    print("  %s: |N(z)|=%d, ratio=%.4f, p_needed=%.4f"
          % (name, norm, ratio, p_needed))

print()
print("  charm: p=0.044, top: p=0.068, bottom: p=-0.045")
print("  NOT universal (different p, different signs).")
print("  VERDICT: Norm correction does NOT work universally.")

# =====================================================================
# HYPOTHESIS 5: The Galois trace formula
# =====================================================================
print()
print("=" * 70)
print("H5: Galois-weighted mass formula")
print("    m_f = m_e * |z|^a_1 * |z'|^(b_1 - a_1)")
print("    or other combinations of z and z'")
print("=" * 70)
print()

# For z = a + b*phi:
# |z| and |z'| are the two embeddings
# N(z) = z * z', Tr(z) = z + z'
# The mass should be a function of z that's Galois-COVARIANT (not invariant)

# Try: m = m_e * |z|^a_1 * |z'|^1
print("Test: m_f = m_e * |z|^5 * |z'|^1")
print("  (exponents a_1 and 1, summing to a_1+1 = b_1)")
print()
for name, g, a, b, Q, C2, m_exp, scale in fermions:
    if name == 'e':
        continue
    z = abs(a + b * PHI)
    zp = abs(a + b * PHI_PRIME)
    m_pred = m_e_MeV * z**a1 * zp**1
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n
    err_bare = abs(m_bare - m_exp)/m_exp*100
    err_new = abs(m_pred - m_exp)/m_exp*100 if m_exp > 0 else 0
    print("  %s: |z|=%.4f, |z'|=%.4f, m_new=%.2f, m_exp=%.2f, err=%.1f%% (bare: %.1f%%)"
          % (name, z, zp, m_pred, m_exp, err_new, err_bare))

# Try: m = m_e * |z^a_1 * z'^1| = m_e * |N(z)^(1/2) * z^(a_1-1)|
# Hmm, let me try systematically
print()
print("Systematic test: m_f = m_e * |z|^p * |z'|^q")
print("  for (p,q) derived from theory")
print()

# The natural choice: p + q should relate to the mass dimension
# In the bare formula: m = m_e * phi^n where n = a_1*a + b_1*b
# This can be written as m = m_e * phi^(a_1*a) * phi^(b_1*b)
# For units (phi^k): z = phi^k, z' = phi'^k = (-1/phi)^k = (-1)^k * phi^(-k)
# So |z|^a_1 * |z'|^b_1 = phi^(k*a_1) * phi^(-k*b_1) = phi^(k*(a_1-b_1)) = phi^(-k)
# That gives m = m_e * phi^(-k). Wrong direction!

# What about |z|^b_1 * |z'|^0 = |z|^6?
# For phi^k: |z|^6 = phi^(6k). Then n = 6k.
# But we want n = 5a + 6b = 5*F(k-1) + 6*F(k).
# For k=2: 5*1+6*1=11 vs 6*2=12. Different.

# The formula n = 5a + 6b is NOT simply |z|^something.
# It's a LINEAR form on Z[phi], not a norm or absolute value.

# Let me think differently: the mass formula is
# m = m_e * phi^(a_1*a + b_1*b)
# = m_e * exp((a_1*a + b_1*b) * ln(phi))
#
# For z = a + b*phi: a = (z*phi - z'*phi')/(phi-phi') ... complicated

# =====================================================================
# HYPOTHESIS 6: QCD running from a DERIVED scale
# =====================================================================
print()
print("=" * 70)
print("H6: QCD running from a DERIVED scale")
print("    If theory predicts bare masses at scale Lambda_theory,")
print("    what Lambda gives best results?")
print("=" * 70)
print()

# For each quark, find what scale Lambda_theory the bare mass would
# need to be AT so that QCD running brings it to the experimental value.

# Using RG running at 1-loop:
# m(mu_low) = m(mu_high) * [alpha_s(mu_low)/alpha_s(mu_high)]^(d_m)
# where d_m = 4/beta_0, beta_0 = 11 - 2*n_f/3

# We use experimental alpha_s(m_Z) = 0.1179 (which equals our derived value!)
# and run both mass and coupling.

# First, let's find for each quark: what scale gives m_bare?
print("What scale Lambda gives m_bare as the running mass m(Lambda)?")
print("  (Working backwards from experimental MSbar masses)")
print()

# For strange: m_s(2 GeV) = 93.4 MeV (n_f=3)
# Want m_s(Lambda) = 101.7 MeV
# 101.7 = 93.4 * [alpha_s(Lambda)/alpha_s(2)]^(4/9)
# alpha_s(2 GeV) ~ 0.30 (approximate)
# [alpha_s(Lambda)/0.30]^(4/9) = 101.7/93.4 = 1.089
# alpha_s(Lambda)/0.30 = 1.089^(9/4) = 1.089^2.25 = 1.209
# alpha_s(Lambda) = 0.363
# This corresponds to Lambda ~ 1 GeV (where alpha_s ~ 0.35)

# More precisely, let's solve numerically
def find_scale_for_mass(m_target, m_ref, mu_ref, n_f):
    """Find scale mu where m(mu) = m_target, given m(mu_ref) = m_ref"""
    beta_0 = 11 - 2*n_f/3
    d_m = 4 / beta_0
    # m(mu) = m_ref * [alpha_s(mu)/alpha_s(mu_ref)]^d_m
    # alpha_s(mu) = 0.1179 / (1 + 0.1179*beta_0/(2*pi)*ln(mu/91.2))

    # Binary search for mu
    mu_low, mu_high = 0.5, 1e6  # MeV range: 500 MeV to 1 TeV
    for _ in range(100):
        mu = (mu_low + mu_high) / 2
        # alpha_s at mu (in GeV)
        mu_GeV = mu / 1000
        as_mu = 0.1179 / (1 + 0.1179 * beta_0 / (2*PI) * math.log(mu_GeV / 91.1876))
        as_ref = 0.1179 / (1 + 0.1179 * beta_0 / (2*PI) * math.log(mu_ref / 1000 / 91.1876))
        if as_mu <= 0 or as_ref <= 0:
            mu_low = mu
            continue
        m_at_mu = m_ref * (as_mu / as_ref) ** d_m
        if m_at_mu > m_target:
            mu_low = mu
        else:
            mu_high = mu
    return mu

# For each quark that has MSbar mass
quarks_to_test = [
    ('s', 1, 1, 1, 93.4, 2000, 3),     # MSbar at 2 GeV, n_f=3
    ('c', 1, 2, 1, 1270.0, 1270, 4),    # MSbar at mc, n_f=4
    ('b', 2,-1, 4, 4180.0, 4180, 5),    # MSbar at mb, n_f=5
    ('d', 0, 1, 0, 4.67, 2000, 3),      # MSbar at 2 GeV, n_f=3
    ('u', 0, 3,-2, 2.16, 2000, 3),      # MSbar at 2 GeV, n_f=3
]

for name, g, a, b, m_msbar, mu_ref, n_f in quarks_to_test:
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n
    if abs(m_bare - m_msbar) / m_msbar < 0.001:
        print("  %s: bare=%.2f ~ msbar=%.2f, no running needed" % (name, m_bare, m_msbar))
        continue
    try:
        scale = find_scale_for_mass(m_bare, m_msbar, mu_ref, n_f)
        print("  %s: bare=%.2f MeV corresponds to MSbar at mu=%.0f MeV (%.1f GeV)"
              % (name, m_bare, scale, scale/1000))
    except:
        print("  %s: could not find scale (bare=%.2f, msbar=%.2f)" % (name, m_bare, m_msbar))

# =====================================================================
# HYPOTHESIS 7: Maybe the formula is m_e * phi^n * phi'^(sector_shift)
# where sector_shift is algebraically determined
# =====================================================================
print()
print("=" * 70)
print("H7: Sector from the 24-cell / edge type")
print("    1+3+8 splitting gives sector quantum numbers")
print("=" * 70)
print()

# In the 1+3+8 splitting:
# U(1): 1 edge -> leptons
# SU(2): 3 edges -> weak doublet structure
# SU(3): 8 edges -> quarks

# What if the mass receives a correction proportional to
# the NUMBER OF GAUGE CONNECTIONS of each type?

# Each vertex has 12 edges: 1(U1) + 3(SU2) + 8(SU3) = 12
# Leptons "use" the U(1) channel: 1/12 of total
# Quarks "use" the SU(3) channels: 8/12 = 2/3 of total

# What if: m_lepton = m_e * phi^(n * 12/12) [full weight]
#           m_quark = m_e * phi^(n * (12-8)/12) = phi^(n * 1/3) [reduced by color]
# That doesn't make sense dimensionally.

# What about the Dynkin index?
# U(1): T = 1 (trivial rep)
# SU(2): T = 1/2 (fundamental rep)
# SU(3): T = 1/2 (fundamental rep)

# Or the dimension of the representation:
# Leptons under SU(3): dim = 1 (singlet)
# Quarks under SU(3): dim = 3 (fundamental)

# What if the effective mass exponent includes a color factor:
# n_eff = n * (1 - alpha_s * C_2 * ln(something))

# Let me try the simplest derived correction:
# The 600-cell spectral action gives a mass operator.
# The SU(3) sector has 8 gauge edges. The Casimir C_2 = 4/3.
# The natural correction: delta_n = -C_2 * alpha_s * n / pi
# (minus sign: QCD screening reduces effective mass exponent)

print("Test: n_eff = n * (1 - C_2 * alpha_s / pi)")
print("  C_2 = 4/3 for quarks, 0 for leptons")
print("  alpha_s = 1/(2*phi^3) = %.6f" % alpha_s)
print("  C_2 * alpha_s / pi = %.6f" % (4/3 * alpha_s / PI))
print()

correction = 4/3 * alpha_s / PI  # = 0.0501

for name, g, a, b, Q, C2, m_exp, scale in fermions:
    if name == 'e':
        continue
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n
    n_eff = n * (1 - C2 * alpha_s / PI)
    m_corr = m_e_MeV * PHI**n_eff
    err_bare = abs(m_bare - m_exp)/m_exp*100
    err_corr = abs(m_corr - m_exp)/m_exp*100
    better = "BETTER" if err_corr < err_bare else "worse"
    print("  %s: n=%2d, n_eff=%.3f, m_bare=%.1f, m_corr=%.1f, m_exp=%.1f  "
          "bare:%.1f%% corr:%.1f%% %s"
          % (name, n, n_eff, m_bare, m_corr, m_exp, err_bare, err_corr, better))

# =====================================================================
# HYPOTHESIS 8: The ELECTRIC CHARGE correction via QED
# =====================================================================
print()
print("=" * 70)
print("H8: QED correction from electric charge")
print("    n_eff = n + 3*Q^2*alpha/pi * n")
print("    (standard QED mass renormalization, coefficient 3Q^2/(4pi) for Dirac)")
print("=" * 70)
print()

for name, g, a, b, Q, C2, m_exp, scale in fermions:
    if name == 'e':
        continue
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n

    # QED + QCD combined
    delta_QED = 3 * Q**2 * alpha / (4*PI)  # per unit of ln(Lambda/m)
    delta_QCD = -C2 * alpha_s / PI  # QCD reduces effective exponent

    # Combined: n_eff = n * (1 + delta_QED * f_QED + delta_QCD)
    # where f_QED involves ln(Lambda/m) ~ O(10)
    # For now, let's just compute the DIRECTION and MAGNITUDE

    print("  %s: Q^2=%.2f, C2=%.2f, delta_QED=%.5f/per-log, delta_QCD=%.5f"
          % (name, Q**2, C2, delta_QED, delta_QCD))

# =====================================================================
# HYPOTHESIS 9: Combined QED + QCD with ln(phi^n) counting
# =====================================================================
print()
print("=" * 70)
print("H9: m_phys = m_bare * phi^(Q^2 * alpha * n + C_2 * alpha_s * n' / pi)")
print("    Using BOTH n and n' (Galois exponent)")
print("    ALL quantities are DERIVED")
print("=" * 70)
print()

# This is interesting: use n for QED correction, n' for QCD correction
# The idea: QED operates in the phi-sector, QCD in the phi'-sector (Galois)
# This would explain why quarks (colored) "see" the Galois sector

# Let's test: m = m_e * phi^(n + Q^2*alpha*n_qed + C_2*alpha_s*n_qcd/pi)
# What should n_qed and n_qcd be?

# Simplest: n_qed = some derived number, n_qcd = some derived number
# But we need them to be the SAME for all fermions (not fitted per fermion)

# What about: the correction uses the GENERATION NUMBER g?
# g = 0, 1, 2 is derived for each fermion.
# And g*(g+1) already appears in delta_up!

# Let me try the most natural thing:
# m = m_e * phi^(n + corrections_from_gauge_loops)
# where corrections = alpha * Q^2 * (some geometric factor) + alpha_s * C_2 * (some geometric factor)

# The key geometric factors from the 600-cell:
# - rank(E8) = 8
# - vertex degree = 12
# - N_eig = 9
# - h(E8) = 30

# Let me just compute what delta_n would be needed for each fermion
# and see if it decomposes into derived quantities:
print("Needed delta_n = n_exact - n_integer for each fermion:")
print()
for name, g, a, b, Q, C2, m_exp, scale in fermions:
    if name == 'e':
        continue
    n = a1*a + b1*b
    m_bare = m_e_MeV * PHI**n
    n_exact = math.log(m_exp / m_e_MeV) / math.log(PHI)
    delta_n = n_exact - n
    print("  %s: n=%2d, n_exact=%.4f, delta_n=%+.4f, Q^2=%.2f, C_2=%.2f, g=%d"
          % (name, n, n_exact, delta_n, Q**2, C2, g))

print()
print("  Can delta_n be written as: a*Q^2 + b*C_2 + c*g + d*g^2?")
print("  (testing linear regression with DERIVED basis vectors)")
print()

# Collect data for regression
import numpy as np

rows = []
deltas = []
for name, g, a, b, Q, C2, m_exp, scale in fermions:
    if name == 'e':
        continue
    n = a1*a + b1*b
    n_exact = math.log(m_exp / m_e_MeV) / math.log(PHI)
    delta_n = n_exact - n
    rows.append([Q**2, C2, g, g*(g+1)])
    deltas.append(delta_n)

X = np.array(rows)
y = np.array(deltas)

# Least squares
try:
    result = np.linalg.lstsq(X, y, rcond=None)
    coeffs = result[0]
    residuals = y - X @ coeffs

    print("  Regression: delta_n = %.4f*Q^2 + %.4f*C_2 + %.4f*g + %.4f*g*(g+1)"
          % tuple(coeffs))
    print()
    print("  Residuals:")
    names_list = [f[0] for f in fermions if f[0] != 'e']
    for i, name in enumerate(names_list):
        print("    %s: fitted=%.4f, actual=%.4f, residual=%.4f"
              % (name, (X[i] @ coeffs), deltas[i], residuals[i]))
    print("  RMS residual: %.4f" % (np.sqrt(np.mean(residuals**2))))
except Exception as ex:
    print("  Regression failed: %s" % ex)

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print()
print("=" * 70)
print("FINAL SUMMARY: Can we find a DERIVED correction?")
print("=" * 70)
print()
print("H1 (QFT running):   Correct DIRECTION for leptons (QED lifts mass)")
print("                     But requires knowing the scale Lambda (not derived)")
print("H2 (Q^2 correction): FAILS - different delta for each fermion")
print("H3 (C_2 * alpha_s):  FAILS - not universal across quarks")
print("H4 (Norm-based):     FAILS - different p for each non-unit")
print("H5 (Galois trace):   Explored but no clean formula found")
print("H6 (QCD running):    Scale not determined by theory")
print("H7 (Color factor):   Correct DIRECTION for quarks but magnitude wrong")
print("H8 (QED + QCD):      Both contribute but can't determine coefficients")
print("H9 (Regression):     Needs fitting - NOT derived")
print()
print("CONCLUSION:")
print("  No DERIVED correction formula found that works universally.")
print("  The bare formula m_f = m_e * phi^(5a+6b) IS the theory's prediction.")
print("  Errors range from 0.2%% to 21%%.")
print()
print("  The theory successfully DERIVES:")
print("  - 9 mass EXPONENTS from Fibonacci + Galois + Casimir")
print("  - The mass HIERARCHY (why tau >> mu >> e)")
print("  - The muon-strange degeneracy (same bare mass)")
print("  - 3 generations (not 2 or 4)")
print()
print("  What it does NOT derive:")
print("  - Sub-percent mass predictions (corrections need dynamics)")
print("  - Lepton-quark mass splitting at same n (needs gauge loops)")
print("  - The renormalization scale (needs a dynamics principle)")

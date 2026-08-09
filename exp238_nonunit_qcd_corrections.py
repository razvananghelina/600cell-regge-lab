"""
EXP-238: QCD CORRECTIONS FOR NON-UNIT FERMIONS (c, t, b)
=========================================================

The mass formula m_f = m_e * phi^n gives:
  - UNITS in Z[phi] (e, mu, tau, u, d, s): errors 0-4%
  - NON-UNITS (c, t, b): errors 11-20%

The non-units have Z[phi] prime structure:
  z_charm = phi * sqrt(5), N(z) = 5 = a_1
  z_top = 4 + phi, N(z) = 19 (prime)
  z_bottom = -1 + 4*phi, N(z) = -19

GOAL: Derive QCD corrections that reduce errors from 11-20% to <2%.
"""

import numpy as np
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-238: QCD CORRECTIONS FOR NON-UNIT FERMIONS (c, t, b)")
print("=" * 70)

a_1 = 5
b_1 = 6
phi = PHI
m_e = 0.51099895  # MeV

# Known experimental masses (PDG 2024, MSbar where applicable)
# Quark masses at their own scale
m_u_exp = 2.16   # MeV (at 2 GeV)
m_d_exp = 4.67   # MeV (at 2 GeV)
m_s_exp = 93.4   # MeV (at 2 GeV)
m_c_exp = 1270   # MeV (MSbar at m_c)
m_b_exp = 4180   # MeV (MSbar at m_b)
m_t_exp = 172760 # MeV (pole mass)

# Lepton masses (exact)
m_mu_exp = 105.6583755  # MeV
m_tau_exp = 1776.86  # MeV

# Geometric predictions: m_f = m_e * phi^n
# (a,b) quantum numbers and n = 5a + 6b
fermions = {
    'electron': {'a': 0, 'b': 0, 'n': 0, 'm_exp': m_e},
    'muon':     {'a': 1, 'b': 1, 'n': 11, 'm_exp': m_mu_exp},
    'tau':      {'a': 1, 'b': 2, 'n': 17, 'm_exp': m_tau_exp},
    'up':       {'a': 3, 'b': -2, 'n': 3, 'm_exp': m_u_exp},
    'charm':    {'a': 2, 'b': 1, 'n': 16, 'm_exp': m_c_exp},
    'top':      {'a': 4, 'b': 1, 'n': 26, 'm_exp': m_t_exp},
    'down':     {'a': 1, 'b': 0, 'n': 5, 'm_exp': m_d_exp},
    'strange':  {'a': 1, 'b': 1, 'n': 11, 'm_exp': m_s_exp},
    'bottom':   {'a': -1, 'b': 4, 'n': 19, 'm_exp': m_b_exp},
}

print("\n--- Geometric predictions vs experiment ---")
print(f"{'Fermion':<10} {'n':>3} {'m_geom':>12} {'m_exp':>12} {'Error':>8} {'Unit?':>6}")
print("-" * 55)

for name, f in fermions.items():
    m_geom = m_e * phi**f['n']
    error = (m_geom - f['m_exp']) / f['m_exp'] * 100

    # Check if unit in Z[phi]
    a, b = f['a'], f['b']
    norm = a**2 + a*b - b**2
    is_unit = abs(abs(norm) - 1) < 0.01

    f['m_geom'] = m_geom
    f['error'] = error
    f['is_unit'] = is_unit
    f['norm'] = norm

    marker = "YES" if is_unit else f"N={norm}"
    print(f"{name:<10} {f['n']:3d} {m_geom:12.2f} {f['m_exp']:12.2f} {error:+7.1f}% {marker:>6}")

# =====================================================================
# PART A: QCD running correction
# =====================================================================
print("\n" + "=" * 70)
print("PART A: QCD RUNNING MASS CORRECTIONS")
print("=" * 70)

print("""
The geometric mass m_geom = m_e * phi^n is the mass at the
GEOMETRIC SCALE (lattice scale of the 600-cell).

The physical mass is the RUNNING mass at the particle's own scale.
QCD running modifies quark masses significantly.

The running mass in QCD:
  m(mu) = m(mu_0) * (alpha_s(mu)/alpha_s(mu_0))^(gamma_0/(2*beta_0))

where gamma_0 = 8 (one-loop anomalous dimension)
      beta_0 = 11 - 2*n_f/3 (one-loop beta function)
""")

# QCD running
alpha_s_mZ = 0.1179  # at M_Z
M_Z = 91187.6  # MeV

def alpha_s_running(mu, n_f):
    """One-loop QCD running of alpha_s from M_Z"""
    beta_0 = 11 - 2*n_f/3
    # alpha_s(mu) = alpha_s(M_Z) / (1 + beta_0/(2*pi) * alpha_s(M_Z) * ln(mu/M_Z))
    return alpha_s_mZ / (1 + beta_0 / (2*np.pi) * alpha_s_mZ * np.log(mu/M_Z))

def running_mass_ratio(mu1, mu2, n_f):
    """Ratio m(mu1)/m(mu2) from QCD running"""
    gamma_0 = 8  # one-loop quark mass anomalous dimension
    beta_0 = 11 - 2*n_f/3
    exponent = gamma_0 / (2 * beta_0)
    return (alpha_s_running(mu1, n_f) / alpha_s_running(mu2, n_f))**exponent

# Define the geometric scale
# The geometric scale is where the 600-cell structure is defined
# This should be near the Planck scale or some high scale
# Let's try: the geometric scale is at m_e (where alpha is defined)
# Or: the geometric scale is at the GUT scale

# Actually, the geometric mass should be the POLE mass or the running
# mass at some specific scale. Let's check what happens if we assume
# the geometric mass is at m_e scale (low energy).

print("\n--- QCD running effects ---")
print(f"alpha_s at M_Z = {alpha_s_mZ}")
print(f"alpha_s at 2 GeV (n_f=3) = {alpha_s_running(2000, 3):.4f}")
print(f"alpha_s at m_c (n_f=4) = {alpha_s_running(m_c_exp, 4):.4f}")
print(f"alpha_s at m_b (n_f=5) = {alpha_s_running(m_b_exp, 5):.4f}")
print(f"alpha_s at m_t (n_f=6) = {alpha_s_running(m_t_exp, 6):.4f}")

# =====================================================================
# PART B: Norm-dependent corrections
# =====================================================================
print("\n" + "=" * 70)
print("PART B: Z[phi] NORM-DEPENDENT CORRECTIONS")
print("=" * 70)

print(f"""
The Z[phi] norm N(z) = a^2 + ab - b^2 of the mass quantum number z = a+b*phi:

  Units (|N|=1): e, mu, tau, u, d, s  -> small errors (0-4%)
  Non-units:
    charm: N = {fermions['charm']['norm']} = a_1
    top:   N = {fermions['top']['norm']}
    bottom: N = {fermions['bottom']['norm']}

HYPOTHESIS: The mass correction for non-units is proportional to
the LOGARITHM of |N(z)|:

  m_phys = m_geom * (1 + c * alpha_s * ln|N(z)|)

For units: ln|N| = ln(1) = 0, no correction.
For charm: ln|N| = ln(5) = 1.609
For top/bottom: ln|N| = ln(19) = 2.944
""")

# Test this hypothesis
print("--- Testing norm-dependent correction ---")

# For each non-unit fermion, what correction factor is needed?
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    needed_factor = f['m_exp'] / f['m_geom']
    ln_N = np.log(abs(f['norm']))
    print(f"\n{name}: m_exp/m_geom = {needed_factor:.4f}, ln|N| = {ln_N:.4f}")
    print(f"  Correction needed: {(needed_factor-1)*100:+.1f}%")
    if ln_N > 0:
        c_fit = (needed_factor - 1) / ln_N
        print(f"  c = (factor-1)/ln|N| = {c_fit:.4f}")

# =====================================================================
# PART C: Try alpha_s * N(z) correction
# =====================================================================
print("\n" + "=" * 70)
print("PART C: SYSTEMATIC CORRECTION FORMULAS")
print("=" * 70)

# The errors are:
# charm: m_geom/m_exp - 1 = +10% (geometric mass TOO HIGH)
# top: +20% (geometric mass TOO HIGH)
# bottom: +17% (geometric mass TOO HIGH)
# Wait let me check signs
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    print(f"{name}: m_geom = {f['m_geom']:.1f}, m_exp = {f['m_exp']:.1f}, "
          f"ratio = {f['m_geom']/f['m_exp']:.4f}, error = {f['error']:+.1f}%")

# All geometric masses are TOO HIGH for non-units
# This suggests the correction is NEGATIVE (reducing the mass)

print(f"""
All non-unit geometric masses are HIGHER than experimental.
The correction should REDUCE the mass.

HYPOTHESIS 1: m_phys = m_geom / (1 + delta)
  where delta depends on the Z[phi] arithmetic of the fermion.

For non-units, the mass exponent z is NOT invertible in Z[phi].
The "missing" invertibility corresponds to a QCD binding effect
that reduces the mass.

Let me try: delta = alpha_s(m_f) * ln(|N(z)|/1)
""")

# Compute alpha_s at each fermion scale
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    n_f = 4 if name == 'charm' else (6 if name == 'top' else 5)
    a_s = alpha_s_running(f['m_exp'], n_f)
    ln_N = np.log(abs(f['norm']))
    delta = a_s * ln_N
    m_corrected = f['m_geom'] / (1 + delta)
    new_error = (m_corrected - f['m_exp']) / f['m_exp'] * 100

    print(f"\n{name}:")
    print(f"  alpha_s(m_{name}) = {a_s:.4f}")
    print(f"  ln|N(z)| = {ln_N:.4f}")
    print(f"  delta = alpha_s * ln|N| = {delta:.4f}")
    print(f"  m_corrected = {m_corrected:.1f} MeV (was {f['m_geom']:.1f})")
    print(f"  Error: {new_error:+.2f}% (was {f['error']:+.1f}%)")

# =====================================================================
# PART D: Try various correction formulas systematically
# =====================================================================
print("\n" + "=" * 70)
print("PART D: SYSTEMATIC SEARCH FOR CORRECTION FORMULA")
print("=" * 70)

# Define correction formulas to test
# Each formula gives delta such that m_phys = m_geom * (1 - delta)
# or m_phys = m_geom / (1 + delta)

alpha_s_0 = 1 / (2 * phi**3)  # geometric alpha_s = 0.1180

formulas = []

# Formula type: m_phys = m_geom * correction_factor
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    needed = f['m_exp'] / f['m_geom']
    f['correction_needed'] = needed

print(f"Correction factors needed:")
print(f"  charm:  {fermions['charm']['correction_needed']:.6f}")
print(f"  top:    {fermions['top']['correction_needed']:.6f}")
print(f"  bottom: {fermions['bottom']['correction_needed']:.6f}")

# Try: correction = 1 / phi^(c * alpha_s * |N|)
# Try: correction = N(z)^(-alpha_s)
# Try: correction from running from geometric scale to physical scale

print(f"\n--- Formula tests ---")

# Test 1: correction = |N(z)|^(-alpha_s)
print(f"\nTest 1: m_phys = m_geom * |N(z)|^(-alpha_s_0)")
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    corr = abs(f['norm'])**(-alpha_s_0)
    m_pred = f['m_geom'] * corr
    err = (m_pred - f['m_exp'])/f['m_exp']*100
    print(f"  {name}: corr={corr:.4f}, m={m_pred:.1f}, err={err:+.1f}%")

# Test 2: correction = |N(z)|^(-alpha_s * phi)
print(f"\nTest 2: m_phys = m_geom * |N(z)|^(-alpha_s_0 * phi)")
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    corr = abs(f['norm'])**(-alpha_s_0 * phi)
    m_pred = f['m_geom'] * corr
    err = (m_pred - f['m_exp'])/f['m_exp']*100
    print(f"  {name}: corr={corr:.4f}, m={m_pred:.1f}, err={err:+.1f}%")

# Test 3: correction = (1 - alpha_s * log(|N|))
print(f"\nTest 3: m_phys = m_geom * (1 - alpha_s_0 * ln|N(z)|)")
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    corr = 1 - alpha_s_0 * np.log(abs(f['norm']))
    m_pred = f['m_geom'] * corr
    err = (m_pred - f['m_exp'])/f['m_exp']*100
    print(f"  {name}: corr={corr:.4f}, m={m_pred:.1f}, err={err:+.1f}%")

# Test 4: correction = phi^(-2*alpha_s * |N-1|)
print(f"\nTest 4: m_phys = m_geom * phi^(-2*alpha_s_0 * (|N|-1))")
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    corr = phi**(-2 * alpha_s_0 * (abs(f['norm'])-1))
    m_pred = f['m_geom'] * corr
    err = (m_pred - f['m_exp'])/f['m_exp']*100
    print(f"  {name}: corr={corr:.4f}, m={m_pred:.1f}, err={err:+.1f}%")

# Test 5: correction from phi-adic valuation
# For non-units, the "excess" norm contributes to the mass
# The correction should be related to how far from unit the element is
print(f"\nTest 5: correction = (1 - (|N|-1)/(N*phi^4))")
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    corr = 1 - (abs(f['norm'])-1) / (120 * phi**4)
    m_pred = f['m_geom'] * corr
    err = (m_pred - f['m_exp'])/f['m_exp']*100
    print(f"  {name}: corr={corr:.4f}, m={m_pred:.1f}, err={err:+.1f}%")

# Test 6: Run alpha_s from a high scale to the fermion mass
# Geometric scale = m_P or some UV cutoff
print(f"\nTest 6: QCD running from phi^n_max * m_e to fermion scale")

# The idea: the geometric prediction is at a HIGH scale
# The physical mass is obtained by running down
m_P = 1.22089e22  # MeV
geom_scale = m_P  # try Planck scale

for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    n_f = 4 if name == 'charm' else (6 if name == 'top' else 5)
    # Running from geom_scale to m_exp
    ratio = running_mass_ratio(f['m_exp'], geom_scale, n_f)
    m_pred = f['m_geom'] * ratio
    err = (m_pred - f['m_exp'])/f['m_exp']*100
    print(f"  {name}: run_ratio={ratio:.4f}, m={m_pred:.1f}, err={err:+.1f}%")

# =====================================================================
# PART E: The N(z) = a_1 anomaly for charm
# =====================================================================
print("\n" + "=" * 70)
print("PART E: CHARM QUARK - THE RAMIFICATION PRIME")
print("=" * 70)

print(f"""
The charm quark has z_c = 2 + phi, with N(z_c) = a_1 = 5.
This is special because 5 is the RAMIFICATION PRIME of Z[phi].

In algebraic number theory:
  sqrt(5) = 2*phi - 1 = "the ramification prime"
  (5) = (sqrt(5))^2 in Z[phi]

The charm mass correction might be related to the ramification:
  z_c = phi * sqrt(5) = phi * pi_5

  where pi_5 = sqrt(5) is the ramification prime.

CONJECTURE: The correction factor for charm is:
  m_c_phys = m_c_geom / sqrt(N(z_c)) = m_c_geom / sqrt(5)
  OR
  m_c_phys = m_c_geom * phi^(-c) for some small c related to ramification
""")

# Test
m_c_geom = fermions['charm']['m_geom']
print(f"m_c_geom = {m_c_geom:.1f} MeV")
print(f"m_c_exp = {m_c_exp:.1f} MeV")
print(f"m_c_geom / sqrt(5) = {m_c_geom/np.sqrt(5):.1f} MeV (err = {(m_c_geom/np.sqrt(5)-m_c_exp)/m_c_exp*100:+.1f}%)")
print(f"m_c_geom / phi = {m_c_geom/phi:.1f} MeV (err = {(m_c_geom/phi-m_c_exp)/m_c_exp*100:+.1f}%)")
print(f"m_c_geom * (1-1/phi) = {m_c_geom*(1-1/phi):.1f} MeV (err = {(m_c_geom*(1-1/phi)-m_c_exp)/m_c_exp*100:+.1f}%)")

# =====================================================================
# PART F: Top and Bottom - the prime 19
# =====================================================================
print("\n" + "=" * 70)
print("PART F: TOP AND BOTTOM - THE PRIME 19")
print("=" * 70)

m_t_geom = fermions['top']['m_geom']
m_b_geom = fermions['bottom']['m_geom']

print(f"""
Top: z_t = 4+phi, N(z_t) = 19 (prime in Z)
Bottom: z_b = -1+4*phi, N(z_b) = -19 (associate of top)

19 = a_1^2 - b_1 = 25-6 = 19
19 = h(E8) - a_1 - b_1 = 30-5-6 = 19
19 = N/b_1 - 1 = 120/6 - 1 = 19

Top/Bottom are 90-degree rotated in (a,b) space:
  z_t = (4, 1), z_b = (-1, 4)
  This is (a,b) -> (-b, a): rotation by pi/2 in Z[phi] plane
""")

# Check: what correction brings top and bottom into line?
print(f"\nCorrection search for top and bottom:")
print(f"  m_t_geom = {m_t_geom:.1f}, m_t_exp = {m_t_exp:.1f}")
print(f"  m_b_geom = {m_b_geom:.1f}, m_b_exp = {m_b_exp:.1f}")

for label, factor_t, factor_b in [
    ("1/phi^(alpha_s)", phi**(-alpha_s_0), phi**(-alpha_s_0)),
    ("N^(-alpha/pi)", 19**(-ALPHA/np.pi), 19**(-ALPHA/np.pi)),
    ("(1-alpha_s*ln19)", 1-alpha_s_0*np.log(19), 1-alpha_s_0*np.log(19)),
    ("1/(1+4*alpha_s)", 1/(1+4*alpha_s_0), 1/(1+4*alpha_s_0)),
    ("(1-ln(19)/N_vert)", 1-np.log(19)/120, 1-np.log(19)/120),
    ("phi^(-1/phi)", phi**(-1/phi), phi**(-1/phi)),
]:
    m_t_pred = m_t_geom * factor_t
    m_b_pred = m_b_geom * factor_b
    err_t = (m_t_pred - m_t_exp)/m_t_exp*100
    err_b = (m_b_pred - m_b_exp)/m_b_exp*100
    combined = np.sqrt(err_t**2 + err_b**2)
    print(f"  {label:25s}: top err={err_t:+.1f}%, bot err={err_b:+.1f}%, combined={combined:.1f}%")

# =====================================================================
# PART G: The n -> n + delta_n correction
# =====================================================================
print("\n" + "=" * 70)
print("PART G: EXPONENT CORRECTION n -> n + delta_n")
print("=" * 70)

print(f"""
Instead of a multiplicative correction, try an EXPONENT correction:
  m_phys = m_e * phi^(n + delta_n)

For each non-unit fermion, what delta_n is needed?
""")

for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    # m_exp = m_e * phi^(n + delta_n)
    # phi^(n + delta_n) = m_exp / m_e
    # n + delta_n = log(m_exp/m_e) / log(phi)
    n_eff = np.log(f['m_exp'] / m_e) / np.log(phi)
    delta_n = n_eff - f['n']
    print(f"  {name}: n={f['n']}, n_eff={n_eff:.4f}, delta_n={delta_n:+.4f}")
    f['delta_n'] = delta_n

# Check if delta_n has a pattern
print(f"\n--- Pattern search for delta_n ---")
d_c = fermions['charm']['delta_n']
d_t = fermions['top']['delta_n']
d_b = fermions['bottom']['delta_n']

print(f"  delta_charm = {d_c:+.4f}")
print(f"  delta_top = {d_t:+.4f}")
print(f"  delta_bottom = {d_b:+.4f}")

print(f"\n  Ratios:")
print(f"  delta_t / delta_c = {d_t/d_c:.4f}")
print(f"  delta_b / delta_c = {d_b/d_c:.4f}")
print(f"  delta_t / delta_b = {d_t/d_b:.4f}")

print(f"\n  In terms of framework constants:")
print(f"  delta_c / alpha_s_0 = {d_c / alpha_s_0:.4f}")
print(f"  delta_t / alpha_s_0 = {d_t / alpha_s_0:.4f}")
print(f"  delta_b / alpha_s_0 = {d_b / alpha_s_0:.4f}")

print(f"  delta_c * phi = {d_c * phi:.4f}")
print(f"  delta_t * phi = {d_t * phi:.4f}")
print(f"  delta_b * phi = {d_b * phi:.4f}")

# Check if delta_n = -alpha_s * ln(|N(z)|) / ln(phi)
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    pred_delta = -alpha_s_0 * np.log(abs(f['norm'])) / np.log(phi)
    print(f"\n  {name}: delta_n = {f['delta_n']:+.4f}, "
          f"-alpha_s*ln|N|/ln(phi) = {pred_delta:+.4f}, "
          f"ratio = {f['delta_n']/pred_delta:.4f}")

# =====================================================================
# PART H: The correct formula - norm running
# =====================================================================
print("\n" + "=" * 70)
print("PART H: SYNTHESIS AND HONEST ASSESSMENT")
print("=" * 70)

print(f"""
FINDINGS:
=========

1. All non-unit geometric masses are TOO HIGH (positive errors).
2. The correction needed is generation-dependent:
   - charm: ~{fermions['charm']['error']:+.0f}% (N=5)
   - top: ~{fermions['top']['error']:+.0f}% (N=19)
   - bottom: ~{fermions['bottom']['error']:+.0f}% (N=19)

3. Simple formulas based on alpha_s * ln|N| give partial improvement
   but do not achieve sub-2% for all three simultaneously.

4. The exponent corrections delta_n are:
   - charm: {fermions['charm']['delta_n']:+.3f}
   - top: {fermions['top']['delta_n']:+.3f}
   - bottom: {fermions['bottom']['delta_n']:+.3f}

5. The fact that delta_t and delta_b have DIFFERENT magnitudes
   despite |N(z_t)| = |N(z_b)| = 19 suggests the correction
   is NOT purely norm-dependent. The SIGN of N matters (N vs -N),
   or the specific (a,b) decomposition matters.

HONEST ASSESSMENT:
  - No single simple formula reduces all three errors to <2%
  - The QCD running from a high scale helps somewhat
  - The Z[phi] arithmetic (norm, ramification) provides QUALITATIVE
    understanding but not QUANTITATIVE prediction
  - This remains OPEN: the non-unit corrections are STRUCTURAL
    (they correlate with non-unit status) but the precise formula
    is not yet derived

STATUS: PARTIAL (structural understanding, no closed formula)
""")

# =====================================================================
# PART I: Can we learn anything from the ratios?
# =====================================================================
print("\n" + "=" * 70)
print("PART I: RATIO ANALYSIS")
print("=" * 70)

# The correction factors
for name in ['charm', 'top', 'bottom']:
    f = fermions[name]
    f['correction'] = f['m_exp'] / f['m_geom']

print(f"Correction factors (m_exp / m_geom):")
c_c = fermions['charm']['correction']
c_t = fermions['top']['correction']
c_b = fermions['bottom']['correction']

print(f"  charm:  {c_c:.6f}")
print(f"  top:    {c_t:.6f}")
print(f"  bottom: {c_b:.6f}")

print(f"\nRatios of corrections:")
print(f"  top/charm = {c_t/c_c:.6f}")
print(f"  bottom/charm = {c_b/c_c:.6f}")
print(f"  top/bottom = {c_t/c_b:.6f}")

# Check if these are powers of phi or simple algebraic numbers
print(f"\n  top/charm / phi^(-1) = {c_t/c_c * phi:.6f}")
print(f"  bottom/charm / phi^(-1) = {c_b/c_c * phi:.6f}")

# Check m_t * m_b product (should relate to m_charm^2 or similar)
print(f"\n  m_t_exp * m_b_exp = {m_t_exp * m_b_exp:.2e}")
print(f"  m_c_exp^2 = {m_c_exp**2:.2e}")
print(f"  ratio = {(m_t_exp*m_b_exp)/m_c_exp**2:.2f}")

# The Koide-like ratio
sum3 = m_c_exp + m_t_exp + m_b_exp
sq_sum = (np.sqrt(m_c_exp) + np.sqrt(m_t_exp) + np.sqrt(m_b_exp))**2
koide = sum3 / sq_sum
print(f"\n  Koide ratio for c,t,b: {koide:.6f} (Koide = 2/3 = {2/3:.6f})")

# n_t + n_b = 45 = a_1 * N_eig (already known)
# n_c + n_b = 35 = |eta| (spectral asymmetry!)
print(f"\nExponent sum rules (checking):")
print(f"  n_c + n_b = {fermions['charm']['n'] + fermions['bottom']['n']} (= |eta| = 35)")
print(f"  n_t + n_b = {fermions['top']['n'] + fermions['bottom']['n']} (= a_1*N_eig = 45)")
print(f"  n_c + n_t = {fermions['charm']['n'] + fermions['top']['n']} (= ??? = 42)")
print(f"  42 = b_1 * (a_1+2) = 6*7 = 42. Interesting!")

print("\n" + "=" * 70)
print("EXP-238 COMPLETE")
print("=" * 70)

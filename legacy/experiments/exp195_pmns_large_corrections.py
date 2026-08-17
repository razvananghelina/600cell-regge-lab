"""
EXP-195: PMNS LARGE CORRECTIONS - THETA_12 AND THETA_23
========================================================

CONTEXT: The CKM angles have small corrections (perturbative QCD/EW).
But PMNS angles have LARGE corrections that can't be perturbative:

Bare PMNS values (phi-tower):
  theta_12 = arctan(phi^-1) = 31.72 deg (exp: 33.41 deg, needs 5.3% correction)
  theta_23 = arctan(phi^0)  = 45.0 deg  (exp: 49.0 deg,  needs ~9% correction)
  theta_13 = arctan(phi^-4) = 8.13 deg  (exp: 8.54 deg,  needs 5.0% correction)

theta_13 matched to 0.11% with 3*alpha_s/(4*pi) correction (exp192).
theta_12 and theta_23 need corrections of order 5-9%, far too large for alpha/pi.

PARTS:
  A: Non-perturbative corrections to theta_12
  B: Non-perturbative corrections to theta_23
  C: Symmetry-based approach (GRM, TBM, BM)
  D: Mass-dependent corrections (charged lepton corrections)
  E: Connection between CKM and PMNS
  F: The "exact" PMNS angles - algebraic expressions in Q(sqrt(5))

Author: Theory of Everything - 600-cell project
Date: 2026-02-10
"""

import math

# =============================================================================
# CONSTANTS
# =============================================================================
PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3      # fine structure constant
ALPHA_S = 0.1179             # strong coupling at M_Z
SIN2_TW = 0.23122            # sin^2(theta_W)
PI = math.pi

# Lepton masses (MeV)
M_E = 0.51099895             # electron
M_MU = 105.6583755           # muon
M_TAU = 1776.86              # tau

# Mass ratios
R_E_MU = M_E / M_MU          # 0.004836
R_MU_TAU = M_MU / M_TAU      # 0.05946
R_E_TAU = M_E / M_TAU        # 0.000288

# PMNS experimental angles (NuFIT 5.2, NO, 2022)
THETA12_EXP_DEG = 33.41       # +/- 0.75
THETA23_EXP_DEG = 49.0        # +/- 1.3 (Normal Ordering)
THETA13_EXP_DEG = 8.54        # +/- 0.15
DELTA_CP_EXP = 197.0          # +/- 25 deg

THETA12_EXP = math.radians(THETA12_EXP_DEG)
THETA23_EXP = math.radians(THETA23_EXP_DEG)
THETA13_EXP = math.radians(THETA13_EXP_DEG)

# Bare PMNS from phi-tower
THETA12_BARE = math.atan(PHI**(-1))   # arctan(1/phi) = 31.717 deg
THETA23_BARE = math.atan(PHI**(0))    # arctan(1) = 45.0 deg
THETA13_BARE = math.atan(PHI**(-4))   # arctan(phi^-4) = 8.127 deg

# CKM for reference
THETA_C = math.radians(12.96)         # Cabibbo angle

# Special symmetry mixing angles
TBM_12 = math.atan(1/math.sqrt(2))    # tribimaximal: 35.264 deg
TBM_23 = PI/4                          # 45 deg
BM_12 = PI/4                           # bimaximal: 45 deg
GRM_12 = math.atan(1/PHI)             # golden ratio: 31.717 deg = our bare


def pct_err(theory, experiment):
    if experiment == 0:
        return float('inf')
    return abs(theory - experiment) / abs(experiment) * 100


def deg(rad):
    return math.degrees(rad)


# =============================================================================
print("=" * 78)
print("EXP-195: PMNS LARGE CORRECTIONS - THETA_12 AND THETA_23")
print("=" * 78)

print("\n--- BARE PMNS ANGLES ---")
for name, bare, exp_deg in [("theta_12", THETA12_BARE, THETA12_EXP_DEG),
                             ("theta_23", THETA23_BARE, THETA23_EXP_DEG),
                             ("theta_13", THETA13_BARE, THETA13_EXP_DEG)]:
    err = pct_err(deg(bare), exp_deg)
    print(f"  {name}: bare = {deg(bare):.4f} deg, exp = {exp_deg:.2f} deg, err = {err:.2f}%")


# =============================================================================
# PART A: NON-PERTURBATIVE CORRECTIONS TO THETA_12
# =============================================================================
print("\n" + "=" * 78)
print("PART A: NON-PERTURBATIVE CORRECTIONS TO THETA_12 (PMNS)")
print("=" * 78)

# 1. What multiplicative correction C do we need?
# arctan(1/phi * (1+C)) = 33.41 deg
tan_exp12 = math.tan(THETA12_EXP)
tan_bare12 = math.tan(THETA12_BARE)
C_needed_12 = tan_exp12 / tan_bare12 - 1
print(f"\n  tan(33.41) = {tan_exp12:.6f}")
print(f"  tan(31.72) = 1/phi = {tan_bare12:.6f}")
print(f"  Multiplicative C needed: C = {C_needed_12:.6f} ({C_needed_12*100:.3f}%)")

# 2. What additive correction delta do we need? (in degrees and radians)
delta_needed_12 = THETA12_EXP - THETA12_BARE
print(f"\n  Additive delta needed: {deg(delta_needed_12):.4f} deg = {delta_needed_12:.6f} rad")

# 3. What epsilon in arctan(1/phi + epsilon) gives theta_12?
# tan(theta_12) = 1/phi + epsilon => epsilon = tan(33.41) - 1/phi
eps_needed_12 = tan_exp12 - 1/PHI
print(f"  Epsilon in arctan(1/phi + eps): eps = {eps_needed_12:.6f}")

# --- Test C candidates ---
print(f"\n  --- Test: arctan(1/phi * (1+C)) ---")
print(f"  C needed = {C_needed_12:.6f}")
C_candidates = {
    'sin2(tW)': SIN2_TW,
    'sin2(tW)/2': SIN2_TW/2,
    'sin2(tW)/3': SIN2_TW/3,
    'sin2(tW)/pi': SIN2_TW/PI,
    '1/phi': 1/PHI,
    '1/phi^2': PHI**(-2),
    '1/phi^3': PHI**(-3),
    '1/phi^4': PHI**(-4),
    '1/phi^5': PHI**(-5),
    'm_e/m_mu': R_E_MU,
    'm_mu/m_tau': R_MU_TAU,
    'alpha*phi': ALPHA*PHI,
    'alpha*phi^2': ALPHA*PHI**2,
    'alpha*phi^3': ALPHA*PHI**3,
    'alpha*phi^4': ALPHA*PHI**4,
    'alpha*phi^5': ALPHA*PHI**5,
    'sqrt(alpha)': math.sqrt(ALPHA),
    'alpha_s/pi': ALPHA_S/PI,
    'alpha_s/(2*pi)': ALPHA_S/(2*PI),
    '3*alpha_s/(4*pi)': 3*ALPHA_S/(4*PI),
    'alpha_s/3': ALPHA_S/3,
    'alpha/pi': ALPHA/PI,
    '2*alpha': 2*ALPHA,
    'alpha_s*phi/(4*pi)': ALPHA_S*PHI/(4*PI),
    'alpha_s^2': ALPHA_S**2,
    'sqrt(alpha_s)/pi': math.sqrt(ALPHA_S)/PI,
    'sqrt(alpha_s)/(2*pi)': math.sqrt(ALPHA_S)/(2*PI),
    'alpha_s/(pi*phi)': ALPHA_S/(PI*PHI),
    'alpha_s*phi/pi': ALPHA_S*PHI/PI,
    '(m_mu/m_tau)/phi': R_MU_TAU/PHI,
    '(m_mu/m_tau)*phi': R_MU_TAU*PHI,
    'sqrt(m_e/m_mu)': math.sqrt(R_E_MU),
    'sqrt(m_mu/m_tau)': math.sqrt(R_MU_TAU),
    '(m_e/m_mu)^(1/3)': R_E_MU**(1/3),
    '(m_mu/m_tau)^(2/3)': R_MU_TAU**(2/3),
    '1/(4*pi)': 1/(4*PI),
    '1/(2*pi)': 1/(2*PI),
    'alpha_s*sin2_tw': ALPHA_S*SIN2_TW,
}

results = []
for name, val in C_candidates.items():
    theta = math.atan(1/PHI * (1 + val))
    err = pct_err(theta, THETA12_EXP)
    ratio = val / C_needed_12
    results.append((name, val, theta, err, ratio))

results.sort(key=lambda x: x[3])
print(f"  {'Expression':<30} {'C value':<12} {'theta(deg)':<12} {'Err%':<10} {'C/C_need':<10}")
print(f"  {'-'*74}")
for name, val, theta, err, ratio in results[:20]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<30} {val:<12.6f} {deg(theta):<12.4f} {err:<10.4f} {ratio:<10.4f}{marker}")

# --- Test additive corrections ---
print(f"\n  --- Test: arctan(1/phi) + delta ---")
print(f"  delta needed = {delta_needed_12:.6f} rad = {deg(delta_needed_12):.4f} deg")
add_candidates = {
    'sqrt(m_e/m_mu) rad': math.sqrt(R_E_MU),
    'm_e/m_mu rad': R_E_MU,
    '(m_e/m_mu)^(1/3) rad': R_E_MU**(1/3),
    '(m_mu/m_tau)^(1/2) rad': math.sqrt(R_MU_TAU),
    'm_mu/m_tau rad': R_MU_TAU,
    'alpha_s/pi rad': ALPHA_S/PI,
    'alpha_s/(2*pi) rad': ALPHA_S/(2*PI),
    '3*alpha_s/(4*pi) rad': 3*ALPHA_S/(4*PI),
    'sin2_tw/pi rad': SIN2_TW/PI,
    'sin2_tw/(2*pi) rad': SIN2_TW/(2*PI),
    '1/phi^4 rad': PHI**(-4),
    '1/(phi^3*pi) rad': 1/(PHI**3*PI),
    '1/(phi^2*pi) rad': 1/(PHI**2*PI),
    'alpha*phi^3 rad': ALPHA*PHI**3,
    'alpha*phi^4 rad': ALPHA*PHI**4,
    'sqrt(alpha) rad': math.sqrt(ALPHA),
    'alpha_s*phi/(3*pi) rad': ALPHA_S*PHI/(3*PI),
}

results_add = []
for name, val in add_candidates.items():
    theta = THETA12_BARE + val
    err = pct_err(theta, THETA12_EXP)
    ratio = val / delta_needed_12
    results_add.append((name, val, theta, err, ratio))

results_add.sort(key=lambda x: x[3])
print(f"  {'Expression':<30} {'delta(rad)':<12} {'theta(deg)':<12} {'Err%':<10} {'d/d_need':<10}")
print(f"  {'-'*74}")
for name, val, theta, err, ratio in results_add[:15]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<30} {val:<12.6f} {deg(theta):<12.4f} {err:<10.4f} {ratio:<10.4f}{marker}")

# --- Test arctan(1/phi + epsilon) ---
print(f"\n  --- Test: arctan(1/phi + epsilon) ---")
print(f"  epsilon needed = {eps_needed_12:.6f}")
eps_candidates = {
    'alpha': ALPHA,
    'alpha*phi': ALPHA*PHI,
    'alpha*phi^2': ALPHA*PHI**2,
    'alpha*phi^3': ALPHA*PHI**3,
    'alpha_s/pi': ALPHA_S/PI,
    'alpha_s/(2*pi)': ALPHA_S/(2*PI),
    'alpha_s/(3*pi)': ALPHA_S/(3*PI),
    '3*alpha_s/(4*pi)': 3*ALPHA_S/(4*PI),
    'm_e/m_mu': R_E_MU,
    'sqrt(m_e/m_mu)': math.sqrt(R_E_MU),
    'm_mu/m_tau': R_MU_TAU,
    'sqrt(m_mu/m_tau)': math.sqrt(R_MU_TAU),
    '(m_e/m_tau)^(1/3)': R_E_TAU**(1/3),
    '1/phi^3': PHI**(-3),
    '1/phi^4': PHI**(-4),
    'sin2_tw/pi': SIN2_TW/PI,
    'alpha_s*phi^2/(4*pi)': ALPHA_S*PHI**2/(4*PI),
}

results_eps = []
for name, val in eps_candidates.items():
    theta = math.atan(1/PHI + val)
    err = pct_err(theta, THETA12_EXP)
    ratio = val / eps_needed_12
    results_eps.append((name, val, theta, err, ratio))

results_eps.sort(key=lambda x: x[3])
print(f"  {'Expression':<30} {'epsilon':<12} {'theta(deg)':<12} {'Err%':<10} {'e/e_need':<10}")
print(f"  {'-'*74}")
for name, val, theta, err, ratio in results_eps[:15]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<30} {val:<12.6f} {deg(theta):<12.4f} {err:<10.4f} {ratio:<10.4f}{marker}")

# --- TBM with corrections ---
print(f"\n  --- Tribimaximal starting point ---")
print(f"  TBM theta_12 = arctan(1/sqrt(2)) = {deg(TBM_12):.4f} deg")
print(f"  Exp theta_12 = {THETA12_EXP_DEG:.2f} deg")
delta_tbm = THETA12_EXP - TBM_12
print(f"  Correction needed (TBM -> exp): {deg(delta_tbm):.4f} deg = {delta_tbm:.6f} rad")
print(f"  Direction: DOWNWARD (TBM is above experiment)")
print(f"  |delta|/TBM = {abs(delta_tbm)/TBM_12*100:.2f}%")

# Can TBM correction be related to theta_C?
print(f"\n  TBM_12 - theta_C = {deg(TBM_12) - deg(THETA_C):.4f} deg")
print(f"  TBM_12 - exp_12 = {deg(TBM_12) - THETA12_EXP_DEG:.4f} deg")
print(f"  Is TBM_12 - exp_12 ~ theta_C/7? {deg(THETA_C)/7:.4f} deg")
print(f"  theta_C/7 = {deg(THETA_C)/7:.4f} vs {abs(deg(TBM_12) - THETA12_EXP_DEG):.4f}")

# GRM starting point
print(f"\n  --- Golden Ratio Mixing (GRM) starting point ---")
print(f"  GRM theta_12 = arctan(1/phi) = {deg(GRM_12):.4f} deg")
print(f"  GRM is BELOW experiment by {THETA12_EXP_DEG - deg(GRM_12):.4f} deg")


# =============================================================================
# PART B: NON-PERTURBATIVE CORRECTIONS TO THETA_23
# =============================================================================
print("\n" + "=" * 78)
print("PART B: NON-PERTURBATIVE CORRECTIONS TO THETA_23 (PMNS)")
print("=" * 78)

# theta_23 bare = pi/4 = 45 deg, exp = 49.0 deg
delta_23 = THETA23_EXP - THETA23_BARE  # positive: need to INCREASE
print(f"\n  Bare theta_23 = pi/4 = 45.0 deg (maximal mixing)")
print(f"  Exp theta_23 = {THETA23_EXP_DEG:.1f} +/- 1.3 deg")
print(f"  Deviation from maximality: delta = {deg(delta_23):.4f} deg = {delta_23:.6f} rad")
print(f"  delta/theta_23(bare) = {delta_23/THETA23_BARE*100:.2f}%")

# Multiplicative correction
tan_exp23 = math.tan(THETA23_EXP)
tan_bare23 = math.tan(THETA23_BARE)
C_needed_23 = tan_exp23 / tan_bare23 - 1
print(f"\n  tan(49.0) = {tan_exp23:.6f}")
print(f"  tan(45.0) = {tan_bare23:.6f}")
print(f"  Multiplicative C needed: C = {C_needed_23:.6f} ({C_needed_23*100:.3f}%)")

# --- Test: theta_23 = pi/4 + delta for various delta ---
print(f"\n  --- Test: theta_23 = pi/4 + delta ---")
print(f"  delta needed = {delta_23:.6f} rad = {deg(delta_23):.4f} deg")

delta23_candidates = {
    'arctan(m_mu/m_tau)': math.atan(R_MU_TAU),
    'm_mu/m_tau rad': R_MU_TAU,
    'sqrt(m_mu/m_tau) rad': math.sqrt(R_MU_TAU),
    '(m_mu/m_tau)^2 rad': R_MU_TAU**2,
    'theta_C/3 rad': THETA_C/3,
    'theta_C/phi^3': THETA_C/PHI**3,
    'theta_13(bare) / 2': THETA13_BARE/2,
    'theta_13(exp) / 2': THETA13_EXP/2,
    'theta_13(bare)': THETA13_BARE,
    'arctan(phi^-4)': math.atan(PHI**(-4)),
    'phi^(-4) rad': PHI**(-4),
    'alpha_s rad': ALPHA_S,
    'alpha_s/phi rad': ALPHA_S/PHI,
    'alpha_s*phi rad': ALPHA_S*PHI,
    '1/phi^2 rad': PHI**(-2),
    '1/phi^3 rad': PHI**(-3),
    'sin2_tw/pi rad': SIN2_TW/PI,
    'sin2_tw/3 rad': SIN2_TW/3,
    'sqrt(alpha_s) rad': math.sqrt(ALPHA_S),
    'alpha_s/(phi*pi) rad': ALPHA_S/(PHI*PI),
    'sqrt(m_e/m_tau) rad': math.sqrt(R_E_TAU),
    '(m_e/m_tau)^(1/3) rad': R_E_TAU**(1/3),
    'm_e/m_mu rad': R_E_MU,
    'pi/phi^5 rad': PI/PHI**5,
    'pi/(5*phi^2) rad': PI/(5*PHI**2),
    'pi*m_mu/(6*m_tau) rad': PI*R_MU_TAU/6,
}

results_23 = []
for name, val in delta23_candidates.items():
    theta = PI/4 + val
    err = pct_err(theta, THETA23_EXP)
    ratio = val / delta_23
    results_23.append((name, val, theta, err, ratio))

results_23.sort(key=lambda x: x[3])
print(f"  {'Expression':<30} {'delta(rad)':<12} {'theta(deg)':<12} {'Err%':<10} {'d/d_need':<10}")
print(f"  {'-'*74}")
for name, val, theta, err, ratio in results_23[:20]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<30} {val:<12.6f} {deg(theta):<12.4f} {err:<10.4f} {ratio:<10.4f}{marker}")

# --- pi/4 + arctan(m_mu/m_tau) ---
theta_23_mass = PI/4 + math.atan(R_MU_TAU)
err_mass = pct_err(theta_23_mass, THETA23_EXP)
print(f"\n  SPECIAL: pi/4 + arctan(m_mu/m_tau) = {deg(theta_23_mass):.4f} deg")
print(f"  Exp = {THETA23_EXP_DEG:.1f} deg, err = {err_mass:.2f}%")
print(f"  Difference from exp: {THETA23_EXP_DEG - deg(theta_23_mass):.4f} deg")

# --- arctan(1 + f(...)) ---
print(f"\n  --- Test: arctan(1 + f(alpha_s, phi, alpha)) ---")
f_candidates = {
    'm_mu/m_tau': R_MU_TAU,
    'sqrt(m_mu/m_tau)': math.sqrt(R_MU_TAU),
    'alpha_s': ALPHA_S,
    'alpha_s/phi': ALPHA_S/PHI,
    'alpha_s*phi': ALPHA_S*PHI,
    '2*alpha_s/(3*pi)': 2*ALPHA_S/(3*PI),
    'sin2_tw': SIN2_TW,
    'sin2_tw/phi': SIN2_TW/PHI,
    '1/phi^2': PHI**(-2),
    '1/phi^3': PHI**(-3),
    'alpha_s + alpha': ALPHA_S + ALPHA,
}

print(f"  C needed (= tan(49)/1 - 1) = {C_needed_23:.6f}")
print(f"\n  {'Expression':<30} {'f value':<12} {'theta(deg)':<12} {'Err%':<10} {'f/C_need':<10}")
print(f"  {'-'*74}")
results_f = []
for name, val in f_candidates.items():
    theta = math.atan(1 + val)
    err = pct_err(theta, THETA23_EXP)
    ratio = val / C_needed_23
    results_f.append((name, val, theta, err, ratio))

results_f.sort(key=lambda x: x[3])
for name, val, theta, err, ratio in results_f[:15]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<30} {val:<12.6f} {deg(theta):<12.4f} {err:<10.4f} {ratio:<10.4f}{marker}")

# --- arctan(phi^eps) for non-integer eps ---
print(f"\n  --- Test: arctan(phi^epsilon) for non-integer epsilon ---")
# arctan(phi^eps) = 49 deg => phi^eps = tan(49) => eps = ln(tan(49))/ln(phi)
eps_effective = math.log(tan_exp23) / math.log(PHI)
print(f"  phi^eps = tan(49.0) = {tan_exp23:.6f}")
print(f"  eps_effective = {eps_effective:.6f}")
print(f"  Bare: eps = 0. So correction = {eps_effective:.6f}")

# What is eps_effective close to?
eps_tests = {
    'm_mu/m_tau': R_MU_TAU,
    'sqrt(m_mu/m_tau)': math.sqrt(R_MU_TAU),
    'alpha_s': ALPHA_S,
    'alpha_s/2': ALPHA_S/2,
    'alpha_s/3': ALPHA_S/3,
    '2*alpha_s/(3*pi)': 2*ALPHA_S/(3*PI),
    'alpha_s/pi': ALPHA_S/PI,
    '3*alpha_s/(4*pi)': 3*ALPHA_S/(4*PI),
    'sin2_tw/phi^2': SIN2_TW/PHI**2,
    '1/phi^3': PHI**(-3),
    '1/phi^4': PHI**(-4),
    '1/(2*pi)': 1/(2*PI),
    'sqrt(alpha)': math.sqrt(ALPHA),
    '1/phi^2 - 1/phi^3': PHI**(-2) - PHI**(-3),
}

results_eff = []
for name, val in eps_tests.items():
    ratio = val / eps_effective
    theta_test = math.atan(PHI**val)
    err = pct_err(theta_test, THETA23_EXP)
    results_eff.append((name, val, theta_test, err, ratio))

results_eff.sort(key=lambda x: x[3])
print(f"  {'Expression':<30} {'epsilon':<12} {'theta(deg)':<12} {'Err%':<10} {'e/e_need':<10}")
print(f"  {'-'*74}")
for name, val, theta_test, err, ratio in results_eff[:15]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<30} {val:<12.6f} {deg(theta_test):<12.4f} {err:<10.4f} {ratio:<10.4f}{marker}")

# --- pi/4 + m_mu/(m_tau * phi^k) ---
print(f"\n  --- Test: pi/4 + m_mu/(m_tau * phi^k) ---")
for k in range(-3, 6):
    val = R_MU_TAU / PHI**k
    theta = PI/4 + val
    err = pct_err(theta, THETA23_EXP)
    marker = " <-- BEST" if err < 1.0 else (" <-- GOOD" if err < 3 else "")
    print(f"    k={k:3d}: delta = {val:.6f} rad, theta = {deg(theta):.4f} deg, err = {err:.3f}%{marker}")


# =============================================================================
# PART C: SYMMETRY-BASED APPROACH
# =============================================================================
print("\n" + "=" * 78)
print("PART C: SYMMETRY-BASED APPROACH")
print("=" * 78)

print("\n  --- Standard symmetry patterns ---")
patterns = {
    'Tribimaximal (TBM)':    (math.atan(1/math.sqrt(2)), PI/4, 0),
    'Bimaximal (BM)':        (PI/4, PI/4, 0),
    'Golden Ratio (GRM)':    (math.atan(1/PHI), PI/4, 0),
    'Hexagonal':             (PI/6, PI/4, 0),
    'Democratic':            (PI/4, PI/4, PI/4),
}

print(f"  {'Pattern':<25} {'th12(deg)':<12} {'th23(deg)':<12} {'th13(deg)':<12}")
print(f"  {'-'*61}")
for name, (t12, t23, t13) in patterns.items():
    print(f"  {name:<25} {deg(t12):<12.4f} {deg(t23):<12.4f} {deg(t13):<12.4f}")
print(f"  {'Experiment':<25} {THETA12_EXP_DEG:<12.2f} {THETA23_EXP_DEG:<12.1f} {THETA13_EXP_DEG:<12.2f}")

# GRM is EXACTLY our bare PMNS!
print(f"\n  NOTE: Golden Ratio Mixing = our bare PMNS from 600-cell!")
print(f"  This means the 600-cell NATURALLY gives GRM symmetry.")

# Corrections FROM GRM
delta_12_grm = THETA12_EXP - GRM_12
delta_23_grm = THETA23_EXP - PI/4
delta_13_grm = THETA13_EXP - 0  # GRM has theta_13 = 0, but our bare has arctan(phi^-4)
print(f"\n  Corrections from GRM to experiment:")
print(f"    delta_12 = +{deg(delta_12_grm):.4f} deg = +{delta_12_grm:.6f} rad")
print(f"    delta_23 = +{deg(delta_23_grm):.4f} deg = +{delta_23_grm:.6f} rad")
print(f"    delta_13 = +{THETA13_EXP_DEG:.4f} deg (from zero)")
print(f"    delta_13(from bare) = +{THETA13_EXP_DEG - deg(THETA13_BARE):.4f} deg")

# sin^2 analysis
sin2_12_grm = 1 / (1 + PHI**2)  # sin^2(arctan(1/phi)) = 1/(1+phi^2) = 1/(2+phi)
sin2_12_exp = math.sin(THETA12_EXP)**2
delta_sin2_12 = sin2_12_exp - sin2_12_grm
print(f"\n  --- sin^2(theta_12) analysis ---")
print(f"  sin^2(GRM) = 1/(1+phi^2) = 1/(2+phi) = {sin2_12_grm:.6f}")
print(f"  sin^2(exp) = sin^2(33.41) = {sin2_12_exp:.6f}")
print(f"  Delta(sin^2) = {delta_sin2_12:.6f}")

# Test what this delta matches
print(f"\n  Is Delta(sin^2) = {delta_sin2_12:.6f} related to known quantities?")
sin2_tests = {
    'alpha': ALPHA,
    'alpha/pi': ALPHA/PI,
    'alpha_s/pi': ALPHA_S/PI,
    'alpha_s/(4*pi)': ALPHA_S/(4*PI),
    'sin2_tw/10': SIN2_TW/10,
    'sin2_tw/8': SIN2_TW/8,
    'sin2_tw/pi^2': SIN2_TW/PI**2,
    'm_e/m_mu': R_E_MU,
    'sqrt(m_e/m_mu)/4': math.sqrt(R_E_MU)/4,
    'm_mu/m_tau': R_MU_TAU,
    'm_mu/(2*m_tau)': R_MU_TAU/2,
    'm_mu/(3*m_tau)': R_MU_TAU/3,
    'm_mu/(phi^2*m_tau)': R_MU_TAU/PHI**2,
    '1/phi^3': PHI**(-3),
    '1/phi^4': PHI**(-4),
    '1/(4*phi^2)': 1/(4*PHI**2),
    'alpha*phi^3': ALPHA*PHI**3,
    'alpha_s*alpha': ALPHA_S*ALPHA,
    '2/(3*phi^3)': 2/(3*PHI**3),
    '1/(2*phi^2*pi)': 1/(2*PHI**2*PI),
    '(2+phi)^(-2)': (2+PHI)**(-2),
    'phi^(-4)*pi/3': PHI**(-4)*PI/3,
}

results_sin2 = []
for name, val in sin2_tests.items():
    ratio = val / delta_sin2_12
    results_sin2.append((name, val, ratio))

results_sin2.sort(key=lambda x: abs(x[2] - 1.0))
print(f"  {'Expression':<30} {'Value':<14} {'Ratio to delta':<14}")
print(f"  {'-'*58}")
for name, val, ratio in results_sin2[:15]:
    marker = " <-- MATCH" if abs(ratio - 1.0) < 0.05 else (" <-- CLOSE" if abs(ratio-1.0) < 0.15 else "")
    print(f"  {name:<30} {val:<14.6f} {ratio:<14.4f}{marker}")

# Test exact algebraic: sin^2(theta_12) = 1/(2+phi) + correction
# Try: sin^2(theta_12) = (phi+1)/(phi^2+phi+3) or similar
# Or: sin^2(theta_12) = 2/(4+phi+1) etc
# Let's try: sin^2 = (a + b*phi) / (c + d*phi) for small integers
print(f"\n  --- Search: sin^2(theta_12) = (a + b*phi)/(c + d*phi) ---")
print(f"  target = {sin2_12_exp:.8f}")
best_alg = []
for a in range(-3, 6):
    for b in range(-3, 6):
        for c in range(1, 12):
            for d in range(-5, 6):
                denom = c + d*PHI
                if denom > 0.01:
                    val = (a + b*PHI) / denom
                    if 0.1 < val < 0.9:
                        err = abs(val - sin2_12_exp)
                        if err < 0.005:
                            best_alg.append((a, b, c, d, val, err))

best_alg.sort(key=lambda x: x[5])
print(f"  {'(a,b,c,d)':<20} {'Value':<14} {'Error':<12}")
print(f"  {'-'*46}")
for a, b, c, d, val, err in best_alg[:15]:
    marker = " <-- EXACT?" if err < 0.0001 else ""
    print(f"  ({a},{b},{c},{d})  => ({a}+{b}*phi)/({c}+{d}*phi)"
          f" = {val:<14.8f} err={err:.8f}{marker}")


# =============================================================================
# PART D: MASS-DEPENDENT CORRECTIONS
# =============================================================================
print("\n" + "=" * 78)
print("PART D: MASS-DEPENDENT CORRECTIONS (CHARGED LEPTON)")
print("=" * 78)

print(f"\n  Lepton mass ratios:")
print(f"    m_e/m_mu  = {R_E_MU:.6f}")
print(f"    m_mu/m_tau = {R_MU_TAU:.6f}")
print(f"    m_e/m_tau  = {R_E_TAU:.8f}")
print(f"    sqrt(m_e/m_mu) = {math.sqrt(R_E_MU):.6f}")
print(f"    sqrt(m_mu/m_tau) = {math.sqrt(R_MU_TAU):.6f}")
print(f"    sqrt(m_e/m_tau) = {math.sqrt(R_E_TAU):.6f}")
print(f"    (m_e/m_mu)^(1/3) = {R_E_MU**(1/3):.6f}")
print(f"    (m_mu/m_tau)^(1/3) = {R_MU_TAU**(1/3):.6f}")

# In Cabibbo-like models: theta_C ~ sqrt(m_d/m_s).
# Analogy for PMNS: theta_12 ~ sqrt(m_e/m_mu)?
print(f"\n  --- Cabibbo analogy: is PMNS theta_12 ~ sqrt(m_e/m_mu)? ---")
print(f"    sqrt(m_e/m_mu) = {math.sqrt(R_E_MU):.6f} rad = {deg(math.sqrt(R_E_MU)):.4f} deg")
print(f"    Exp theta_12 = {THETA12_EXP_DEG:.2f} deg")
print(f"    Not directly, but as CORRECTION:")
print(f"    arctan(1/phi) + sqrt(m_e/m_mu) = {deg(THETA12_BARE + math.sqrt(R_E_MU)):.4f} deg")
print(f"    (too large by {deg(THETA12_BARE + math.sqrt(R_E_MU)) - THETA12_EXP_DEG:.4f} deg)")

# More systematic: U_PMNS = U_nu^dag * U_e
# If U_nu = GRM and U_e provides corrections from charged lepton sector
print(f"\n  --- Model: U_PMNS = U_nu^dag * U_e ---")
print(f"  If U_nu gives GRM, then U_e gives charged lepton corrections")
print(f"  Typical estimates: theta_e ~ sqrt(m_e/m_mu) = {math.sqrt(R_E_MU):.6f} rad")
print(f"  Charged lepton rotation angle ~ {deg(math.sqrt(R_E_MU)):.4f} deg")

# Specific mass-based formulas for theta_12
print(f"\n  --- Mass-based formulas for theta_12 correction ---")
mass_formulas_12 = {
    'arctan(1/phi) + sqrt(m_e/m_mu)': THETA12_BARE + math.sqrt(R_E_MU),
    'arctan(1/phi) + (m_e/m_mu)^(1/3)': THETA12_BARE + R_E_MU**(1/3),
    'arctan(1/phi + sqrt(m_e/m_mu))': math.atan(1/PHI + math.sqrt(R_E_MU)),
    'arctan(1/phi + m_e/m_mu)': math.atan(1/PHI + R_E_MU),
    'arctan(1/phi * (1 + sqrt(m_e/m_mu)))': math.atan(1/PHI * (1 + math.sqrt(R_E_MU))),
    'arctan(1/phi * (1 + m_mu/m_tau))': math.atan(1/PHI * (1 + R_MU_TAU)),
    'arctan(1/phi) + m_mu/(3*m_tau)': THETA12_BARE + R_MU_TAU/3,
    'arctan(1/phi * (1+m_mu/(2*m_tau)))': math.atan(1/PHI * (1 + R_MU_TAU/2)),
    'arctan(1/phi * (1+sqrt(m_e/m_tau)))': math.atan(1/PHI * (1 + math.sqrt(R_E_TAU))),
    'arctan(1/phi + sqrt(m_e/m_tau))': math.atan(1/PHI + math.sqrt(R_E_TAU)),
    'arctan(1/phi) + (m_e/m_tau)^(1/3)': THETA12_BARE + R_E_TAU**(1/3),
    'arctan(1/phi * (1+(m_mu/m_tau)^2))': math.atan(1/PHI * (1 + R_MU_TAU**2)),
    'arctan(1/phi) + m_e^2/(m_mu*m_tau)': THETA12_BARE + M_E**2/(M_MU*M_TAU),
    'arctan((1+sqrt(m_e/m_mu))/phi)': math.atan((1+math.sqrt(R_E_MU))/PHI),
}

print(f"  {'Formula':<45} {'theta(deg)':<12} {'Err%':<10}")
print(f"  {'-'*67}")
results_mass12 = []
for name, val in mass_formulas_12.items():
    err = pct_err(val, THETA12_EXP)
    results_mass12.append((name, val, err))

results_mass12.sort(key=lambda x: x[2])
for name, val, err in results_mass12[:15]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<45} {deg(val):<12.4f} {err:<10.4f}{marker}")

# Mass-based formulas for theta_23
print(f"\n  --- Mass-based formulas for theta_23 correction ---")
mass_formulas_23 = {
    'pi/4 + arctan(m_mu/m_tau)': PI/4 + math.atan(R_MU_TAU),
    'pi/4 + m_mu/m_tau': PI/4 + R_MU_TAU,
    'pi/4 + sqrt(m_mu/m_tau)': PI/4 + math.sqrt(R_MU_TAU),
    'pi/4 + (m_mu/m_tau)^(2/3)': PI/4 + R_MU_TAU**(2/3),
    'pi/4 + (m_mu/m_tau)^(1/3)': PI/4 + R_MU_TAU**(1/3),
    'arctan(1 + m_mu/m_tau)': math.atan(1 + R_MU_TAU),
    'arctan(1 + sqrt(m_mu/m_tau))': math.atan(1 + math.sqrt(R_MU_TAU)),
    'arctan(1 + sqrt(m_e/m_mu))': math.atan(1 + math.sqrt(R_E_MU)),
    'pi/4 + m_mu/m_tau + sqrt(m_e/m_mu)': PI/4 + R_MU_TAU + math.sqrt(R_E_MU),
    'pi/4 + m_mu/(m_tau*phi)': PI/4 + R_MU_TAU/PHI,
    'pi/4 + m_mu*phi/m_tau': PI/4 + R_MU_TAU*PHI,
    'arctan(1 + 2*m_mu/(3*m_tau))': math.atan(1 + 2*R_MU_TAU/3),
    'pi/4 + (m_mu/m_tau)/phi^2': PI/4 + R_MU_TAU/PHI**2,
    'pi/4 + (m_mu/m_tau)*phi^2': PI/4 + R_MU_TAU*PHI**2,
    'arctan(1 + m_mu/m_tau + m_e/m_mu)': math.atan(1 + R_MU_TAU + R_E_MU),
}

print(f"  {'Formula':<45} {'theta(deg)':<12} {'Err%':<10}")
print(f"  {'-'*67}")
results_mass23 = []
for name, val in mass_formulas_23.items():
    err = pct_err(val, THETA23_EXP)
    results_mass23.append((name, val, err))

results_mass23.sort(key=lambda x: x[2])
for name, val, err in results_mass23[:15]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<45} {deg(val):<12.4f} {err:<10.4f}{marker}")

# Using exact mass ratios from the 600-cell formula
print(f"\n  --- Using 600-cell mass formula m_f = m_e * phi^(a*d + b*k) ---")
print(f"  600-cell masses: e(0,0), mu(1,1), tau(1,2)")
print(f"  mu:  m_e * phi^(1*5 + 1*6) = m_e * phi^11 = {M_E * PHI**11:.4f} MeV (exp: {M_MU:.4f})")
print(f"  tau: m_e * phi^(1*5 + 2*6) = m_e * phi^17 = {M_E * PHI**17:.4f} MeV (exp: {M_TAU:.4f})")
mu_600 = M_E * PHI**11
tau_600 = M_E * PHI**17
r_mu_tau_600 = mu_600 / tau_600
print(f"  m_mu/m_tau (600-cell) = phi^(11-17) = phi^-6 = {PHI**(-6):.6f}")
print(f"  m_mu/m_tau (exp) = {R_MU_TAU:.6f}")
print(f"  m_mu/m_tau (600) = {r_mu_tau_600:.6f}")

# Using phi^-6 instead of m_mu/m_tau
print(f"\n  Using phi^-6 as mass ratio:")
print(f"  pi/4 + phi^(-6) = {deg(PI/4 + PHI**(-6)):.4f} deg (err = {pct_err(PI/4 + PHI**(-6), THETA23_EXP):.3f}%)")
print(f"  arctan(1 + phi^(-6)) = {deg(math.atan(1 + PHI**(-6))):.4f} deg (err = {pct_err(math.atan(1 + PHI**(-6)), THETA23_EXP):.3f}%)")
print(f"  pi/4 + phi^(-6)/phi = pi/4 + phi^(-7) = {deg(PI/4 + PHI**(-7)):.4f} deg (err = {pct_err(PI/4 + PHI**(-7), THETA23_EXP):.3f}%)")


# =============================================================================
# PART E: CONNECTION BETWEEN CKM AND PMNS
# =============================================================================
print("\n" + "=" * 78)
print("PART E: CONNECTION BETWEEN CKM AND PMNS")
print("=" * 78)

# Factor-3 rule
print(f"\n  --- Factor-3 rule: n_CKM = 3 * n_PMNS ---")
print(f"  theta_12: CKM n=3, PMNS n=1.  3 = 3*1  CHECK")
print(f"  theta_13: CKM n=12, PMNS n=4. 12 = 3*4 CHECK")
print(f"  theta_23: CKM n=7, PMNS n=0.   7 != 3*0 FAIL")
print(f"  For theta_23: 7/3 = {7/3:.4f}")

# What if PMNS theta_23 has effective n = 7/3?
theta_23_7_3 = math.atan(PHI**(-7/3))
print(f"\n  If PMNS theta_23 uses n=7/3:")
print(f"  arctan(phi^(-7/3)) = {deg(theta_23_7_3):.4f} deg")
print(f"  Exp = {THETA23_EXP_DEG:.1f} deg, err = {pct_err(theta_23_7_3, THETA23_EXP):.2f}%")
print(f"  This is {deg(theta_23_7_3):.2f} deg, nowhere near 49 deg. REJECTED.")

# What effective n_PMNS gives theta_23?
# arctan(phi^(-n)) = 49 deg => phi^-n = tan(49) => n = -ln(tan(49))/ln(phi)
n_eff_23 = -math.log(tan_exp23) / math.log(PHI)
print(f"\n  Effective n for exp theta_23:")
print(f"  n_eff = -ln(tan(49))/ln(phi) = {n_eff_23:.6f}")
print(f"  This is NEGATIVE: {n_eff_23:.4f}")
print(f"  Meaning phi^(-n_eff) = phi^{-n_eff_23:.4f} = phi^{abs(n_eff_23):.4f}")
print(f"  (theta_23 > 45 deg means tan > 1, so exponent is positive)")

# Quark-lepton complementarity (QLC)
print(f"\n  --- Quark-lepton complementarity (QLC) ---")
print(f"  theta_C + theta_12(PMNS) = {deg(THETA_C):.2f} + {THETA12_EXP_DEG:.2f} = {deg(THETA_C) + THETA12_EXP_DEG:.2f} deg")
print(f"  Compare to 45 deg: off by {abs(deg(THETA_C) + THETA12_EXP_DEG - 45):.2f} deg")
print(f"")
print(f"  theta_23(CKM) + theta_23(PMNS) = {deg(math.asin(0.0408)):.2f} + {THETA23_EXP_DEG:.1f} = {deg(math.asin(0.0408)) + THETA23_EXP_DEG:.2f} deg")
print(f"  theta_13(CKM) + theta_13(PMNS) = {deg(math.asin(0.00382)):.3f} + {THETA13_EXP_DEG:.2f} = {deg(math.asin(0.00382)) + THETA13_EXP_DEG:.3f} deg")

# More precise: in bare phi-tower
print(f"\n  --- In bare phi-tower ---")
print(f"  Bare: arctan(phi^-3) + arctan(phi^-1) = {deg(math.atan(PHI**(-3)) + math.atan(PHI**(-1))):.4f} deg")
print(f"  Compare to 45 deg: off by {abs(deg(math.atan(PHI**(-3)) + math.atan(PHI**(-1))) - 45):.4f} deg")

# Test: arctan(phi^-a) + arctan(phi^-b) = pi/4? => phi^-a + phi^-b = 1
# This uses the identity arctan(x) + arctan(y) = pi/4 when xy + x + y = 1 (approximately)
# Actually arctan(x) + arctan(y) = arctan((x+y)/(1-xy)) (when xy < 1)
# = pi/4 when (x+y)/(1-xy) = 1, i.e., x+y = 1-xy, i.e., x+y+xy = 1
print(f"\n  Test: is phi^-3 + phi^-1 + phi^-3*phi^-1 = 1?")
x = PHI**(-3)
y = PHI**(-1)
print(f"  phi^-3 + phi^-1 + phi^-4 = {x + y + x*y:.6f}")
print(f"  Compare to 1: off by {abs(x + y + x*y - 1):.6f}")
print(f"  This uses arctan(a)+arctan(b)=pi/4 iff a+b+ab=1")

# What about the CORRECTED sum?
theta_C_bare = math.atan(PHI**(-3))
qlc_sum = theta_C_bare + THETA12_BARE
print(f"\n  Bare sum: arctan(phi^-3) + arctan(phi^-1) = {deg(qlc_sum):.4f} deg")
print(f"  Compare to pi/4 = 45.0 deg")
print(f"  The QLC sum BARE is {deg(qlc_sum):.4f} deg, off by {abs(deg(qlc_sum) - 45):.4f} deg")

# What if QLC is exact at the GUT scale and corrections break it?
print(f"\n  Experimental QLC: theta_C + theta_12(PMNS) = {deg(THETA_C) + THETA12_EXP_DEG:.2f} deg")
print(f"  If QLC is exact at GUT scale (sum = 45 deg),")
print(f"  then the 1.37 deg excess comes from RG running + corrections")

# Factor of 3 for all angles - modified version
print(f"\n  --- Modified factor-3 rule ---")
print(f"  CKM:  tan(theta_12)=phi^-3,  tan(theta_23)=phi^-7,  tan(theta_13)=phi^-12")
print(f"  PMNS: tan(theta_12)=phi^-1,  tan(theta_23)=phi^0,   tan(theta_13)=phi^-4")
print(f"  Exponent ratios: 3/1=3, 7/0=inf, 12/4=3")
print(f"  For theta_23: CKM_n / PMNS_n is undefined (PMNS_n=0)")
print(f"  Alternative: PMNS_n = CKM_n - 2*N_gen? (3-6<0, 7-6=1, 12-6=6)")
print(f"  Alternative: PMNS_n = CKM_n / N_gen? (3/3=1, 7/3=2.33, 12/3=4)")
print(f"  For theta_23: n_eff = 7/3 = 2.333")
print(f"  arctan(phi^(-2.333)) = {deg(math.atan(PHI**(-7/3))):.4f} deg (exp 49 - NO)")


# =============================================================================
# PART F: EXACT PMNS ANGLES IN Q(sqrt(5))
# =============================================================================
print("\n" + "=" * 78)
print("PART F: ARE EXPERIMENTAL PMNS ANGLES ALGEBRAIC IN Q(sqrt(5))?")
print("=" * 78)

print(f"\n  --- tan(theta) for experimental PMNS ---")
tan12 = math.tan(THETA12_EXP)
tan23 = math.tan(THETA23_EXP)
tan13 = math.tan(THETA13_EXP)
print(f"  tan(33.41) = {tan12:.8f}")
print(f"  tan(49.0)  = {tan23:.8f}")
print(f"  tan(8.54)  = {tan13:.8f}")

# Search: is tan(theta) = a + b*phi for small integers?
print(f"\n  --- Search: tan(theta) = a + b*phi for integers a,b (|a|,|b| <= 5) ---")
for label, tan_val in [("theta_12", tan12), ("theta_23", tan23), ("theta_13", tan13)]:
    print(f"\n  {label}: tan = {tan_val:.8f}")
    best = []
    for a in range(-5, 6):
        for b in range(-5, 6):
            val = a + b * PHI
            if val > 0:
                err = abs(val - tan_val)
                if err < 0.1:
                    best.append((a, b, val, err, pct_err(val, tan_val)))
    best.sort(key=lambda x: x[3])
    if best:
        print(f"  {'(a,b)':<12} {'a+b*phi':<14} {'Error':<12} {'Err%':<10}")
        print(f"  {'-'*48}")
        for a, b, val, err, perr in best[:5]:
            print(f"  ({a},{b}){'':<6} {val:<14.8f} {err:<12.8f} {perr:<10.4f}")
    else:
        print(f"  No match found for |a|,|b| <= 5 with error < 0.1")

# Search with fractions: tan(theta) = (a + b*phi) / (c + d*phi)
print(f"\n  --- Search: tan(theta) = (a + b*phi)/(c + d*phi) for small integers ---")
for label, tan_val in [("theta_12", tan12), ("theta_23", tan23), ("theta_13", tan13)]:
    print(f"\n  {label}: tan = {tan_val:.8f}")
    best_frac = []
    for a in range(-4, 5):
        for b in range(-4, 5):
            for c in range(-4, 5):
                for d in range(-4, 5):
                    denom = c + d*PHI
                    if abs(denom) > 0.01:
                        val = (a + b*PHI) / denom
                        if val > 0:
                            err = abs(val - tan_val)
                            if err < 0.01:
                                complexity = abs(a) + abs(b) + abs(c) + abs(d)
                                best_frac.append((a, b, c, d, val, err, complexity))
    best_frac.sort(key=lambda x: (x[5], x[6]))
    if best_frac:
        print(f"  {'(a,b,c,d)':<20} {'Value':<14} {'Error':<12}")
        print(f"  {'-'*46}")
        for a, b, c, d, val, err, comp in best_frac[:8]:
            marker = " <-- EXACT?" if err < 0.001 else ""
            print(f"  ({a},{b},{c},{d}){'':<8} {val:<14.8f} {err:<12.8f}{marker}")
    else:
        print(f"  No match found")

# What about phi-powers and combinations?
print(f"\n  --- More exotic: tan(theta) = phi^a * sqrt(b) + c ---")
for label, tan_val in [("theta_12", tan12), ("theta_23", tan23), ("theta_13", tan13)]:
    print(f"\n  {label}: tan = {tan_val:.8f}")
    exotic = []
    for a in range(-6, 7):
        for b in [1, 2, 3, 5, 6, 7, 10]:
            for c_int in range(-3, 4):
                val = PHI**a * math.sqrt(b) + c_int
                if val > 0:
                    err = abs(val - tan_val)
                    if err < 0.005:
                        exotic.append((a, b, c_int, val, err))
    exotic.sort(key=lambda x: x[4])
    if exotic:
        print(f"  {'phi^a * sqrt(b) + c':<30} {'Value':<14} {'Error':<12}")
        print(f"  {'-'*56}")
        for a, b, c_int, val, err in exotic[:5]:
            print(f"  phi^{a}*sqrt({b})+{c_int}{'':<14} {val:<14.8f} {err:<12.8f}")
    else:
        print(f"  No match found")


# =============================================================================
# PART G: COMBINED BEST RESULTS AND NEW IDEAS
# =============================================================================
print("\n" + "=" * 78)
print("PART G: COMBINED ANALYSIS AND BEST RESULTS")
print("=" * 78)

# Comprehensive table: for each PMNS angle, what is the best correction?
print(f"\n  --- COMPREHENSIVE BEST-FIT TABLE ---")

# theta_13: already solved (exp192)
theta13_corr = math.atan(PHI**(-4) * (1 + 3*ALPHA_S/(4*PI)))
err13 = pct_err(theta13_corr, THETA13_EXP)
print(f"\n  THETA_13: arctan(phi^-4 * (1 + 3*alpha_s/(4*pi)))")
print(f"    = {deg(theta13_corr):.4f} deg (exp: {THETA13_EXP_DEG:.2f} deg)")
print(f"    Error: {err13:.3f}%  [PATTERN - QCD loop with C_F factor]")

# theta_12: systematic scan of best formulas
print(f"\n  THETA_12: scanning combined corrections...")
best_12 = []

# Mass-based
formulas_12 = [
    ('arctan(1/phi*(1+m_mu/m_tau))', math.atan(1/PHI*(1+R_MU_TAU))),
    ('arctan(1/phi*(1+m_mu/(2*m_tau)))', math.atan(1/PHI*(1+R_MU_TAU/2))),
    ('arctan(1/phi + alpha_s/(2*pi))', math.atan(1/PHI + ALPHA_S/(2*PI))),
    ('arctan(1/phi + 3*as/(4*pi))', math.atan(1/PHI + 3*ALPHA_S/(4*PI))),
    ('arctan(1/phi*(1+sin2_tw/3))', math.atan(1/PHI*(1+SIN2_TW/3))),
    ('arctan(1/phi*(1+sin2_tw/pi))', math.atan(1/PHI*(1+SIN2_TW/PI))),
    ('arctan(1/phi*(1+alpha_s/(pi*phi)))', math.atan(1/PHI*(1+ALPHA_S/(PI*PHI)))),
    ('arctan(1/phi*(1+as*phi/pi))', math.atan(1/PHI*(1+ALPHA_S*PHI/PI))),
    ('arctan(1/phi*(1+1/(4*pi)))', math.atan(1/PHI*(1+1/(4*PI)))),
    ('arctan(1/phi*(1+1/phi^3))', math.atan(1/PHI*(1+PHI**(-3)))),
    ('arctan(1/phi + (m_e/m_tau)^(1/3))', math.atan(1/PHI + R_E_TAU**(1/3))),
    ('arctan(1/phi*(1+(m_mu/m_tau)^2))', math.atan(1/PHI*(1+R_MU_TAU**2))),
    ('arctan(1/phi + phi^(-5))', math.atan(1/PHI + PHI**(-5))),
    ('arctan((1+sqrt(m_e/m_mu))/phi)', math.atan((1+math.sqrt(R_E_MU))/PHI)),
    ('arctan(1/phi*(1+sqrt(m_e/m_tau)))', math.atan(1/PHI*(1+math.sqrt(R_E_TAU)))),
    ('arctan(1/phi*(1+as*sin2_tw))', math.atan(1/PHI*(1+ALPHA_S*SIN2_TW))),
    ('TBM - arctan(phi^-4)', TBM_12 - math.atan(PHI**(-4))),
    ('TBM - phi^(-5)', TBM_12 - PHI**(-5)),
    ('arctan(1/phi) + as/phi^4', THETA12_BARE + ALPHA_S/PHI**4),
    ('arctan(1/phi) + alpha*phi^5', THETA12_BARE + ALPHA*PHI**5),
]

for name, val in formulas_12:
    err = pct_err(val, THETA12_EXP)
    best_12.append((name, val, err))

best_12.sort(key=lambda x: x[2])
print(f"  {'Formula':<50} {'theta(deg)':<12} {'Err%':<10}")
print(f"  {'-'*72}")
for name, val, err in best_12[:20]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<50} {deg(val):<12.4f} {err:<10.4f}{marker}")

# theta_23: systematic scan of best formulas
print(f"\n  THETA_23: scanning combined corrections...")
best_23 = []
formulas_23 = [
    ('pi/4 + arctan(m_mu/m_tau)', PI/4 + math.atan(R_MU_TAU)),
    ('pi/4 + m_mu/m_tau', PI/4 + R_MU_TAU),
    ('pi/4 + m_mu*phi/m_tau', PI/4 + R_MU_TAU*PHI),
    ('pi/4 + m_mu/(m_tau*phi)', PI/4 + R_MU_TAU/PHI),
    ('pi/4 + m_mu/(m_tau*phi^2)', PI/4 + R_MU_TAU/PHI**2),
    ('arctan(1 + m_mu/m_tau)', math.atan(1 + R_MU_TAU)),
    ('arctan(1 + sqrt(m_mu/m_tau))', math.atan(1 + math.sqrt(R_MU_TAU))),
    ('pi/4 + sqrt(m_mu/m_tau)', PI/4 + math.sqrt(R_MU_TAU)),
    ('pi/4 + (m_mu/m_tau)^(1/3)', PI/4 + R_MU_TAU**(1/3)),
    ('pi/4 + alpha_s', PI/4 + ALPHA_S),
    ('pi/4 + alpha_s/phi', PI/4 + ALPHA_S/PHI),
    ('pi/4 + alpha_s*phi', PI/4 + ALPHA_S*PHI),
    ('pi/4 + 1/phi^3', PI/4 + PHI**(-3)),
    ('pi/4 + sin2_tw/3', PI/4 + SIN2_TW/3),
    ('arctan(1 + alpha_s/phi)', math.atan(1 + ALPHA_S/PHI)),
    ('arctan(1+m_mu/m_tau+m_e/m_mu)', math.atan(1+R_MU_TAU+R_E_MU)),
    ('pi/4 + phi^(-6)', PI/4 + PHI**(-6)),
    ('pi/4 + phi^(-6)*phi', PI/4 + PHI**(-5)),
    ('arctan(1+2*m_mu/(3*m_tau))', math.atan(1+2*R_MU_TAU/3)),
    ('pi/4+m_mu/m_tau+sqrt(m_e/m_mu)', PI/4 + R_MU_TAU + math.sqrt(R_E_MU)),
    ('pi/4 + arctan(m_mu/m_tau)+m_e/m_mu', PI/4 + math.atan(R_MU_TAU) + R_E_MU),
]

for name, val in formulas_23:
    err = pct_err(val, THETA23_EXP)
    best_23.append((name, val, err))

best_23.sort(key=lambda x: x[2])
print(f"  {'Formula':<50} {'theta(deg)':<12} {'Err%':<10}")
print(f"  {'-'*72}")
for name, val, err in best_23[:20]:
    marker = " <-- EXCELLENT" if err < 0.5 else (" <-- GOOD" if err < 2 else "")
    print(f"  {name:<50} {deg(val):<12.4f} {err:<10.4f}{marker}")

# Check if the theta_23 deviation relates to theta_13
print(f"\n  --- Is the theta_23 deviation related to theta_13? ---")
print(f"  delta_23 = {deg(delta_23):.4f} deg")
print(f"  theta_13(bare) = {deg(THETA13_BARE):.4f} deg")
print(f"  theta_13(exp)  = {THETA13_EXP_DEG:.2f} deg")
print(f"  delta_23 / theta_13(bare) = {deg(delta_23)/deg(THETA13_BARE):.4f}")
print(f"  delta_23 / theta_13(exp)  = {deg(delta_23)/THETA13_EXP_DEG:.4f}")
print(f"  theta_13(bare) / 2 = {deg(THETA13_BARE)/2:.4f} deg")
print(f"  Is delta_23 ~ theta_13/2? Off by {abs(deg(delta_23) - THETA13_EXP_DEG/2):.4f} deg")

# Seesaw-type: delta_23 ~ theta_13^2 / (theta_12)
ratio_seesaw = THETA13_EXP**2 / THETA12_EXP
print(f"  theta_13^2 / theta_12 = {deg(ratio_seesaw):.4f} deg vs delta_23 = {deg(delta_23):.4f} deg")
print(f"  Ratio: {deg(delta_23)/deg(ratio_seesaw):.4f}")


# =============================================================================
# PART H: COMBINED FORMULA ATTEMPT
# =============================================================================
print("\n" + "=" * 78)
print("PART H: COMBINED FORMULA STRUCTURE")
print("=" * 78)

print(f"""
  The 600-cell naturally gives GOLDEN RATIO MIXING (GRM):
    theta_12 = arctan(1/phi) = 31.72 deg   [from n=1 on phi-tower]
    theta_23 = arctan(1)     = 45.00 deg   [from n=0, maximal]
    theta_13 = arctan(phi^-4)= 8.13 deg    [from n=4 on phi-tower]

  The corrections needed are:
    delta_12 = +1.69 deg = +0.0295 rad  ({delta_12_grm/THETA12_BARE*100:.2f}% of bare)
    delta_23 = +4.00 deg = +0.0698 rad  ({delta_23/THETA23_BARE*100:.2f}% of bare)
    delta_13 = +0.41 deg = +0.0072 rad  ({(THETA13_EXP - THETA13_BARE)/THETA13_BARE*100:.2f}% of bare)
""")

# Classification of corrections by magnitude
print(f"  Corrections by magnitude:")
print(f"  theta_13: {(THETA13_EXP - THETA13_BARE)/THETA13_BARE*100:.2f}% -> 3*alpha_s/(4*pi) = {3*ALPHA_S/(4*PI)*100:.2f}% [PERTURBATIVE]")
print(f"  theta_12: {C_needed_12*100:.2f}% -> ???")
print(f"  theta_23: {C_needed_23*100:.2f}% -> ???")

# What if ALL corrections come from the same mechanism but with different coefficients?
# theta_13: C_F * alpha_s / (2*pi) * (1 factor) ~ 0.025
# theta_12: C_12 ~ 0.067 ~ 2.7 * 0.025
# theta_23: C_23 ~ 0.15 ~ 6.0 * 0.025
# Ratios: 1 : 2.7 : 6.0. Is this related to the phi-tower exponents? 4:1:0

C_13 = 3*ALPHA_S/(4*PI)  # ~ 0.0281
C_12 = C_needed_12
C_23 = C_needed_23
print(f"\n  If all corrections are multiples of C_13 = 3*as/(4*pi) = {C_13:.6f}:")
print(f"    C_12 / C_13 = {C_12/C_13:.4f}")
print(f"    C_23 / C_13 = {C_23/C_13:.4f}")
print(f"    Is C_12/C_13 ~ phi? phi = {PHI:.4f} (ratio = {C_12/C_13/PHI:.4f})")
print(f"    Is C_12/C_13 ~ phi^2? phi^2 = {PHI**2:.4f} (ratio = {C_12/C_13/PHI**2:.4f})")
print(f"    Is C_23/C_13 ~ phi^3? phi^3 = {PHI**3:.4f} (ratio = {C_23/C_13/PHI**3:.4f})")
print(f"    Is C_23/C_13 ~ 5? (ratio = {C_23/C_13/5:.4f})")
print(f"    Is C_23/C_13 ~ phi^2+1? {PHI**2+1:.4f} (ratio = {C_23/C_13/(PHI**2+1):.4f})")

# Pattern: C_ij = (something) * 3*alpha_s/(4*pi) * phi^(4-n)?
# n=4: phi^0 = 1 (theta_13, exact)
# n=1: phi^3 = 4.236 (theta_12, need 2.37, ratio = 0.56)
# n=0: phi^4 = 6.854 (theta_23, need 5.32, ratio = 0.78)
for n, label, C_val in [(4, "theta_13", C_13), (1, "theta_12", C_12), (0, "theta_23", C_23)]:
    phi_factor = PHI**(4-n)
    predicted = C_13 * phi_factor
    ratio = C_val / predicted
    print(f"\n  {label} (n={n}): C_13 * phi^(4-{n}) = {C_13} * {phi_factor:.4f} = {predicted:.6f}")
    print(f"    Actual C_{label[-2:]} = {C_val:.6f}, ratio = {ratio:.4f}")

# Alternative: is C proportional to (n_CKM - n_PMNS)?
# theta_13: n_CKM=12, n_PMNS=4, diff=8
# theta_12: n_CKM=3, n_PMNS=1, diff=2
# theta_23: n_CKM=7, n_PMNS=0, diff=7
print(f"\n  --- Alternative: C ~ (n_CKM - n_PMNS) * small ---")
for pmns_n, ckm_n, label, C_val in [(4,12,"theta_13",C_13), (1,3,"theta_12",C_12), (0,7,"theta_23",C_23)]:
    diff = ckm_n - pmns_n
    C_per_diff = C_val / diff if diff > 0 else float('inf')
    print(f"  {label}: n_CKM-n_PMNS = {diff}, C/{diff} = {C_per_diff:.6f}")

# =============================================================================
# PART I: EXACT MATCH ATTEMPTS WITH EXPERIMENTAL UNCERTAINTIES
# =============================================================================
print("\n" + "=" * 78)
print("PART I: WORKING WITHIN EXPERIMENTAL UNCERTAINTIES")
print("=" * 78)

# theta_23 has the largest uncertainty: 49.0 +/- 1.3 deg
# So values between 47.7 and 50.3 are acceptable
theta23_lo = math.radians(49.0 - 1.3)
theta23_hi = math.radians(49.0 + 1.3)
theta12_lo = math.radians(33.41 - 0.75)
theta12_hi = math.radians(33.41 + 0.75)
theta13_lo = math.radians(8.54 - 0.15)
theta13_hi = math.radians(8.54 + 0.15)

print(f"\n  Experimental ranges (1 sigma):")
print(f"  theta_12: {33.41-0.75:.2f} - {33.41+0.75:.2f} deg")
print(f"  theta_23: {49.0-1.3:.1f} - {50.3:.1f} deg")
print(f"  theta_13: {8.54-0.15:.2f} - {8.54+0.15:.2f} deg")

# For theta_23: can any nice formula fall within 1-sigma?
print(f"\n  --- theta_23: formulas within 1-sigma of 49.0 +/- 1.3 ---")
nice_23 = [
    ('pi/4 + arctan(m_mu/m_tau)', PI/4 + math.atan(R_MU_TAU)),
    ('pi/4 + m_mu/m_tau', PI/4 + R_MU_TAU),
    ('pi/4 + m_mu*phi/m_tau', PI/4 + R_MU_TAU*PHI),
    ('arctan(1+m_mu/m_tau)', math.atan(1+R_MU_TAU)),
    ('pi/4 + arctan(m_mu/m_tau)+m_e/m_mu', PI/4 + math.atan(R_MU_TAU) + R_E_MU),
    ('pi/4 + alpha_s/phi', PI/4 + ALPHA_S/PHI),
    ('arctan(1+alpha_s/phi)', math.atan(1+ALPHA_S/PHI)),
    ('pi/4 + phi^(-6)', PI/4 + PHI**(-6)),
    ('arctan(phi^0.29)', math.atan(PHI**eps_effective)),
    ('pi/4 + 1/phi^3', PI/4 + PHI**(-3)),
    ('arctan(1+sin2_tw/phi^2)', math.atan(1+SIN2_TW/PHI**2)),
    ('arctan(1+(m_mu/m_tau)^(2/3))', math.atan(1+R_MU_TAU**(2/3))),
]

print(f"  {'Formula':<50} {'theta(deg)':<12} {'Within 1sig?':<14} {'Err%':<10}")
print(f"  {'-'*86}")
for name, val in nice_23:
    err = pct_err(val, THETA23_EXP)
    within = "YES" if theta23_lo <= val <= theta23_hi else "no"
    marker = " <-- MATCH" if within == "YES" else ""
    print(f"  {name:<50} {deg(val):<12.4f} {within:<14}{err:<10.4f}{marker}")

# For theta_12: formulas within 1-sigma
print(f"\n  --- theta_12: formulas within 1-sigma of 33.41 +/- 0.75 ---")
nice_12 = [
    ('arctan(1/phi*(1+m_mu/m_tau))', math.atan(1/PHI*(1+R_MU_TAU))),
    ('arctan(1/phi + 3*as/(4*pi))', math.atan(1/PHI + 3*ALPHA_S/(4*PI))),
    ('arctan(1/phi*(1+sin2_tw/3))', math.atan(1/PHI*(1+SIN2_TW/3))),
    ('arctan(1/phi*(1+1/(4*pi)))', math.atan(1/PHI*(1+1/(4*PI)))),
    ('arctan(1/phi*(1+(m_mu/m_tau)^2))', math.atan(1/PHI*(1+R_MU_TAU**2))),
    ('arctan(1/phi + phi^(-5))', math.atan(1/PHI + PHI**(-5))),
    ('arctan(1/phi*(1+as*sin2_tw))', math.atan(1/PHI*(1+ALPHA_S*SIN2_TW))),
    ('arctan(1/phi*(1+sqrt(m_e/m_tau)))', math.atan(1/PHI*(1+math.sqrt(R_E_TAU)))),
    ('TBM - phi^(-5)', TBM_12 - PHI**(-5)),
    ('arctan(1/phi*(1+as/(pi*phi)))', math.atan(1/PHI*(1+ALPHA_S/(PI*PHI)))),
]

print(f"  {'Formula':<50} {'theta(deg)':<12} {'Within 1sig?':<14} {'Err%':<10}")
print(f"  {'-'*86}")
for name, val in nice_12:
    err = pct_err(val, THETA12_EXP)
    within = "YES" if theta12_lo <= val <= theta12_hi else "no"
    marker = " <-- MATCH" if within == "YES" else ""
    print(f"  {name:<50} {deg(val):<12.4f} {within:<14}{err:<10.4f}{marker}")


# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 78)
print("FINAL SUMMARY - EXP 195")
print("=" * 78)

print("""
+====================================================================+
|        PMNS LARGE CORRECTIONS - RESULTS TABLE                       |
+====================================================================+

BARE PMNS (Golden Ratio Mixing from 600-cell):
  theta_12 = arctan(1/phi)   = 31.72 deg (exp: 33.41, err: 5.1%)
  theta_23 = arctan(1) = pi/4 = 45.00 deg (exp: 49.0,  err: 8.2%)
  theta_13 = arctan(phi^-4)  =  8.13 deg (exp: 8.54,   err: 4.7%)

CORRECTED (exp192 result for theta_13):
  theta_13 = arctan(phi^-4 * (1 + 3*as/(4*pi)))
""")
print(f"           = {deg(theta13_corr):.4f} deg (exp: {THETA13_EXP_DEG:.2f}, err: {err13:.3f}%)")

print(f"""
KEY FINDING: CORRECTIONS SCALE WITH PHI-TOWER LEVEL
  The multiplicative correction C needed for each angle:
  theta_13 (n=4): C = 3*as/(4*pi) = {C_13:.6f} = {C_13*100:.3f}%
  theta_12 (n=1): C = {C_12:.6f} = {C_12*100:.3f}%
  theta_23 (n=0): C = {C_needed_23:.6f} = {C_needed_23*100:.3f}%

  The corrections GROW as n DECREASES:
  n=4 -> 2.8%   (perturbative)
  n=1 -> 6.7%   (borderline)
  n=0 -> 15.0%  (non-perturbative)

THETA_12 ANALYSIS:
""")

# Print best theta_12 results
best_12.sort(key=lambda x: x[2])
for name, val, err in best_12[:5]:
    status = "GOOD" if err < 2 else "OK"
    print(f"  {name:<50} {deg(val):.4f} deg  err={err:.3f}%  [{status}]")

print(f"""
THETA_23 ANALYSIS:
""")
best_23.sort(key=lambda x: x[2])
for name, val, err in best_23[:5]:
    status = "GOOD" if err < 2 else "OK"
    print(f"  {name:<50} {deg(val):.4f} deg  err={err:.3f}%  [{status}]")

print(f"""
SYMMETRY ANALYSIS:
  The 600-cell gives EXACTLY Golden Ratio Mixing (GRM).
  GRM: sin^2(theta_12) = 1/(1+phi^2) = 1/(2+phi) = {sin2_12_grm:.6f}
  Exp: sin^2(theta_12) = {sin2_12_exp:.6f}
  Delta(sin^2) = {delta_sin2_12:.6f}

QUARK-LEPTON COMPLEMENTARITY:
  theta_C + theta_12(PMNS) = {deg(THETA_C) + THETA12_EXP_DEG:.2f} deg (vs 45 deg, off by {abs(deg(THETA_C) + THETA12_EXP_DEG - 45):.2f} deg)
  Bare: arctan(phi^-3) + arctan(phi^-1) = {deg(math.atan(PHI**(-3)) + math.atan(PHI**(-1))):.4f} deg

FACTOR-3 RULE:
  theta_12: n_CKM=3, n_PMNS=1.  3/1=3  CHECK
  theta_13: n_CKM=12, n_PMNS=4. 12/4=3 CHECK
  theta_23: n_CKM=7, n_PMNS=0.  UNDEFINED (0 in denominator)

CLASSIFICATION:
  theta_13 correction: PATTERN (0.11%) - perturbative, 3*as/(4*pi)
  theta_12 correction: OPEN - multiple candidates at 1-3%, none compelling
  theta_23 correction: OPEN - large non-perturbative, mass-based best
  GRM = 600-cell bare:  DERIVED (structure from phi-tower)
  Factor-3 rule:         PATTERN (2 of 3 angles)
  QLC:                   PATTERN (approximate, 1.4 deg off from 45)

OPEN QUESTIONS:
  1. What mechanism gives 5-15% corrections to PMNS bare values?
  2. Is it charged lepton diagonalization (U_PMNS = U_nu^dag * U_e)?
  3. Is it RG running from GUT scale (but 1-loop too small)?
  4. Is there a non-perturbative 600-cell mechanism (soliton shifts)?
  5. Why does theta_23 deviate from maximal mixing by exactly ~4 deg?
  6. The correction C grows as 1/phi^n - is this a geometric series?
""")

print("=" * 78)
print("END EXP-195")
print("=" * 78)

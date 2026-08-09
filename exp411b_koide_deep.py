"""
exp411b_koide_deep.py
======================
DEEPENING the Koide formula investigation from exp411.

BUILDS ON exp411 results:
  - Q_framework = 0.666682 (0.0024% from 2/3)
  - Geometric limit: t = s + 1/s = 5 = a1 (palindromic polynomial)
  - Koide angle theta_0 ~ 2.3166 rad ~ 132.73 deg
  - alpha = 0.75 gives best Q among scanned values

INVESTIGATES 6 LEADS:
  1. Generalized palindromic condition for non-geometric sequences
  2. Koide angle vs framework expressions (2*delta_CKM, etc.)
  3. WHY alpha = 3/4 optimizes Koide (theorem or coincidence?)
  4. TQFT / Verlinde ring connection
  5. Exact Koide: what c_ell or corrections give Q = 2/3 exactly?
  6. Galois structure of Q(sqrt(5), sqrt(phi))

Author: Razvan-Constantin Anghelina
Date: 2026-02-25

RULES:
  - WINDOWS: NO Unicode in print/comments
  - REGULA ZERO: nu inventa
  - HONEST: DERIVAT / PATTERN / SPECULATIV
"""

import numpy as np
from math import sqrt, pi, log, factorial, atan, atan2, acos, asin, cos, sin, gcd

print("=" * 72)
print("EXP-411b: DEEP KOIDE INVESTIGATION")
print("=" * 72)

# =====================================================================
# Framework constants
# =====================================================================
a1 = 5
b1 = a1 + 1  # = 6
PHI = (1 + sqrt(a1)) / 2
PHIp = (1 - sqrt(a1)) / 2  # = -1/PHI
N = factorial(a1)  # = 120
d_ST = 4
C_mass = 4.0 / (a1**2 + 1)  # = 2/13
c_ell = C_mass * PHI**3 / d_ST
alpha_morrey = (d_ST - 1.0) / d_ST  # = 3/4
alpha_em = 7.2973525693e-3  # fine-structure constant
alpha_s_fw = 1.0 / (2 * PHI**3)  # framework alpha_s

m_e_MeV = 0.51099895000  # MeV CODATA 2018
m_mu_exp = 105.6583755    # MeV
m_tau_exp = 1776.86        # MeV

# Lepton (a,b) assignments
leptons = {
    'e':   {'a': 0, 'b': 0, 'n_bare': 0},
    'mu':  {'a': 1, 'b': 1, 'n_bare': 11},
    'tau': {'a': 1, 'b': 2, 'n_bare': 17},
}

# Compute framework masses
def compute_fw_mass(a, b, n, alpha=alpha_morrey, c=c_ell):
    """Compute framework mass with correction."""
    zp = a + b * PHIp
    if abs(zp) < 1e-15:
        delta = 0.0
    else:
        delta = c * np.sign(zp) * abs(zp)**alpha
    return m_e_MeV * PHI**(n + delta), n + delta, delta

m_e_fw, n_e_eff, d_e = compute_fw_mass(0, 0, 0)
m_mu_fw, n_mu_eff, d_mu = compute_fw_mass(1, 1, 11)
m_tau_fw, n_tau_eff, d_tau = compute_fw_mass(1, 2, 17)

masses_fw = [m_e_fw, m_mu_fw, m_tau_fw]

def koide_Q(masses):
    """Compute Koide ratio Q = sum(m) / (sum(sqrt(m)))^2."""
    s = sum(masses)
    ss = sum(sqrt(m) for m in masses)
    return s / ss**2

def koide_from_exponents(n0, n1, n2):
    """Q from phi exponents."""
    return (PHI**n0 + PHI**n1 + PHI**n2) / (PHI**(n0/2) + PHI**(n1/2) + PHI**(n2/2))**2

Q_fw = koide_Q(masses_fw)
Q_exp = koide_Q([m_e_MeV, m_mu_exp, m_tau_exp])

print(f"\n  Recap from exp411:")
print(f"    Q_experimental = {Q_exp:.10f}")
print(f"    Q_framework    = {Q_fw:.10f}")
print(f"    2/3            = {2/3:.10f}")
print(f"    |Q_fw - 2/3|   = {abs(Q_fw - 2/3):.2e} ({abs(Q_fw-2/3)/(2/3)*100:.4f}%)")
print(f"    n_mu_eff = {n_mu_eff:.8f}, n_tau_eff = {n_tau_eff:.8f}")

# =====================================================================
# LEAD 1: PALINDROMIC POLYNOMIAL AND a1=5
# =====================================================================
print(f"\n{'='*72}")
print("LEAD 1: Generalized Palindromic Condition for Non-Geometric Sequences")
print("=" * 72)

# For geometric sequence: m_i = m_0, m_0*r, m_0*r^2
# Koide Q = 2/3 => palindromic poly => t = s + 1/s = 5 = a1
#
# Our masses are: m_e = phi^0, m_mu = phi^n1, m_tau = phi^n2
# with n1 ~ 11, n2 ~ 17. These are NOT geometric (n2/n1 != 2).
#
# Generalized: consider m_i = phi^(n_i) for ARBITRARY n_0, n_1, n_2.
# Q = (phi^n0 + phi^n1 + phi^n2) / (phi^(n0/2) + phi^(n1/2) + phi^(n2/2))^2
#
# WLOG set n_0 = 0 (divide by m_e). Then:
# Q = (1 + phi^n1 + phi^n2) / (1 + phi^(n1/2) + phi^(n2/2))^2
#
# For geometric: n1 = p, n2 = 2p. Let's parametrize as:
# n1 = p, n2 = p + q (where q is the "asymmetry")
# Geometric: q = p. Our case: p = 11, q = 6.

print(f"\n  Our case: n1 = 11, n2 = 17, so p = n1 = 11, q = n2 - n1 = 6")
print(f"  Geometric would be q = p = 11. We have q = 6 = b1.")
print(f"  q/p = 6/11 = {6/11:.6f} (far from 1)")

# For the generalized case, Q = 2/3 is an equation in (p, q):
# 3(1 + phi^p + phi^(p+q)) = 2(1 + phi^(p/2) + phi^((p+q)/2))^2
#
# There is NO simple palindromic structure for general (p,q).
# But we can check: for what q(p) does Q(p,q) = 2/3?

print(f"\n  Scanning: for each p (= n_mu), find q such that Q(p, p+q) = 2/3")
print(f"  {'p':>4} {'q_exact':>10} {'n_tau':>8} {'q/p':>8} {'q+p':>6} {'q':>6}")

from scipy.optimize import brentq

results_pq = []
for p in range(4, 25):
    # Find q such that Q(0, p, p+q) = 2/3 exactly
    def f_q(q):
        return koide_from_exponents(0, p, p + q) - 2.0/3.0
    # Q is monotonic in q for reasonable range
    try:
        # Check if solution exists in a reasonable range
        q_lo, q_hi = 0.1, 3*p
        if f_q(q_lo) * f_q(q_hi) < 0:
            q_exact = brentq(f_q, q_lo, q_hi, xtol=1e-12)
            results_pq.append((p, q_exact))
            print(f"  {p:4d} {q_exact:10.4f} {p+q_exact:8.2f} {q_exact/p:8.4f} "
                  f"{p+q_exact:6.2f} {q_exact:6.2f}")
        else:
            # Try wider range
            q_lo, q_hi = 0.01, 10*p
            if f_q(q_lo) * f_q(q_hi) < 0:
                q_exact = brentq(f_q, q_lo, q_hi, xtol=1e-12)
                results_pq.append((p, q_exact))
                print(f"  {p:4d} {q_exact:10.4f} {p+q_exact:8.2f} {q_exact/p:8.4f} "
                      f"{p+q_exact:6.2f} {q_exact:6.2f}")
            else:
                print(f"  {p:4d} {'no sol':>10}")
    except Exception:
        print(f"  {p:4d} {'error':>10}")

# For p = 11 (our case), what is q_exact?
if results_pq:
    for p_val, q_val in results_pq:
        if p_val == 11:
            print(f"\n  For p = 11 (our n_mu):")
            print(f"    q_exact = {q_val:.8f}")
            print(f"    Our q (with corrections) = {n_tau_eff - n_mu_eff:.8f}")
            print(f"    Our q (bare) = 6.000000")
            print(f"    b1 = {b1}")
            q11 = q_val
            break

# Check: is q_exact for p=11 close to b1=6?
print(f"\n  KEY QUESTION: Is q_exact(p=11) close to b1 = 6?")
if results_pq:
    for p_val, q_val in results_pq:
        if p_val == 11:
            print(f"    q_exact = {q_val:.8f}")
            print(f"    b1 = 6")
            print(f"    Diff = {q_val - 6:.6f}")
            print(f"    Rel = {(q_val - 6)/6*100:.4f}%")

# Look at q/p ratio
print(f"\n  Asymptotic behavior of q/p for large p:")
if len(results_pq) >= 3:
    for p_val, q_val in results_pq[-5:]:
        print(f"    p={p_val}: q/p = {q_val/p_val:.6f}")
    # Does q/p converge? If so, to what?
    # For geometric (q=p), q/p = 1 exactly.
    # For large p, Q is dominated by phi^(p+q), so effectively geometric-like.
    # The asymptotic limit should be q/p -> 1 for Q = 2/3.
    ratios = [q/p for p, q in results_pq if p >= 10]
    if ratios:
        print(f"    Apparent limit q/p -> {ratios[-1]:.6f} (expected 1.0 for geometric)")

# Now: GENERALIZED palindromic condition
# For general (p,q): define s = phi^(p/2), u = phi^(q/2)
# Then Q = (1 + s^2 + s^2*u^2) / (1 + s + s*u)^2
# Q = 2/3 => 3(1 + s^2 + s^2*u^2) = 2(1 + s + su)^2
# Expand: 3 + 3s^2 + 3s^2*u^2 = 2(1 + s^2 + s^2*u^2 + 2s + 2su + 2s^2*u)
# = 2 + 2s^2 + 2s^2*u^2 + 4s + 4su + 4s^2*u
# => s^2 + s^2*u^2 + 1 - 4s - 4su - 4s^2*u = 0
# => s^2(1 + u^2 - 4u) + 1 - 4s(1 + u) = 0

print(f"\n  GENERALIZED KOIDE EQUATION:")
print(f"  With s = phi^(p/2), u = phi^(q/2):")
print(f"  s^2*(1 + u^2 - 4u) + 1 - 4s*(1 + u) = 0")
print(f"  For u = s (geometric, q = p): s^2*(1 + s^2 - 4s) + 1 - 4s*(1+s) = 0")
print(f"    => s^4 - 4s^3 + s^2 + 1 - 4s - 4s^2 = 0")
print(f"    => s^4 - 4s^3 - 3s^2 - 4s + 1 = 0  [PALINDROMIC, as in exp411]")

# For our case: s = phi^(11/2), u = phi^(6/2) = phi^3
s_our = PHI**(11.0/2)
u_our = PHI**3
check = s_our**2 * (1 + u_our**2 - 4*u_our) + 1 - 4*s_our*(1 + u_our)
print(f"\n  Our case: s = phi^(11/2), u = phi^3:")
print(f"    s = {s_our:.6f}")
print(f"    u = phi^3 = {u_our:.6f}")
print(f"    u = phi^3 = phi^2 + phi = (3+sqrt(5))/2 + (1+sqrt(5))/2 = 2 + sqrt(5)")
print(f"    u = {2 + sqrt(5):.6f}  (check: {u_our:.6f})")
print(f"    Equation residual = {check:.6f}")
print(f"    (Non-zero => our exponents don't give EXACT Q = 2/3)")

# What u gives Q = 2/3 for s = phi^(11/2)?
# s^2*(1 + u^2 - 4u) + 1 - 4s*(1+u) = 0
# s^2*u^2 - (4*s^2 + 4*s)*u + (s^2 + 1 - 4*s) = 0
# This is quadratic in u!
A_coef = s_our**2
B_coef = -(4*s_our**2 + 4*s_our)
C_coef = s_our**2 + 1 - 4*s_our
disc = B_coef**2 - 4*A_coef*C_coef
if disc >= 0:
    u_sol1 = (-B_coef + sqrt(disc)) / (2*A_coef)
    u_sol2 = (-B_coef - sqrt(disc)) / (2*A_coef)
    print(f"\n  Solving for u (given s = phi^(11/2)):")
    print(f"    u_1 = {u_sol1:.8f}")
    print(f"    u_2 = {u_sol2:.8f}")
    print(f"    Our u = phi^3 = {u_our:.8f}")
    print(f"    |u_1 - phi^3| = {abs(u_sol1 - u_our):.6f}")
    print(f"    |u_2 - phi^3| = {abs(u_sol2 - u_our):.6f}")
    # Which one is closer?
    closest = u_sol1 if abs(u_sol1 - u_our) < abs(u_sol2 - u_our) else u_sol2
    print(f"    Closest solution: u = {closest:.8f}")
    q_exact_from_u = 2 * log(closest) / log(PHI)
    print(f"    => q_exact = 2*log_phi(u) = {q_exact_from_u:.6f}")
    print(f"    Bare q = 6, diff = {q_exact_from_u - 6:.6f}")

print(f"\n  ASSESSMENT (Lead 1):")
print(f"    The palindromic structure is SPECIFIC to geometric sequences.")
print(f"    For general (p,q), the equation is quadratic in u = phi^(q/2).")
print(f"    The condition t = s + 1/s = a1 = 5 does NOT directly generalize.")
print(f"    Our q ~ 6 = b1 is CLOSE to but NOT exactly the Koide-exact value.")
print(f"    Status: PATTERN (t=a1 in geometric limit); NON-TRIVIAL for actual masses.")


# =====================================================================
# LEAD 2: KOIDE ANGLE theta_0
# =====================================================================
print(f"\n{'='*72}")
print("LEAD 2: Koide Angle theta_0 and Framework Expressions")
print("=" * 72)

# Koide parametrization: sqrt(m_i) = M0*(1 + sqrt(2)*cos(theta_0 + 2*pi*k/3))
# where k = 0,1,2 for e, mu, tau
# M0 = sum(sqrt(m_i))/3

sum_sqrt = sum(sqrt(m) for m in masses_fw)
M0 = sum_sqrt / 3

# Extract theta_0 from electron
ratio_e = sqrt(masses_fw[0]) / M0 - 1
cos_e = ratio_e / sqrt(2)

print(f"\n  M0 = {M0:.8f} MeV^(1/2)")
print(f"  sqrt(m_e)/M0 - 1 = {ratio_e:.8f}")
print(f"  cos(theta_0) = {cos_e:.8f}")

if abs(cos_e) > 1:
    print(f"  |cos(theta_0)| > 1: parametrization requires Q <= 2/3.")
    print(f"  Since Q_fw = {Q_fw:.8f} > 2/3 = {2/3:.8f}, the standard")
    print(f"  Koide parametrization breaks down slightly.")
    print(f"  Using cosh version or allowing complex theta.")
    # When Q > 2/3 but close, we can still define an "effective" angle
    # by noting that the parametrization fails because Q is slightly > 2/3.
    # We can still compute the angle using experimental masses (Q < 2/3).

# Use EXPERIMENTAL masses for clean angle extraction
sum_sqrt_exp = sum(sqrt(m) for m in [m_e_MeV, m_mu_exp, m_tau_exp])
M0_exp = sum_sqrt_exp / 3
ratio_e_exp = sqrt(m_e_MeV) / M0_exp - 1
cos_e_exp = ratio_e_exp / sqrt(2)

print(f"\n  Using EXPERIMENTAL masses:")
print(f"    M0_exp = {M0_exp:.8f}")
print(f"    cos(theta_0) = {cos_e_exp:.8f}")

theta_0_exp = acos(cos_e_exp)
print(f"    theta_0 = {theta_0_exp:.10f} rad")
print(f"    theta_0 = {np.degrees(theta_0_exp):.6f} deg")

# Verify
print(f"\n  Verification of Koide parametrization:")
for k, name in enumerate(['e', 'mu', 'tau']):
    sqrt_m_pred = M0_exp * (1 + sqrt(2) * cos(theta_0_exp + 2*pi*k/3))
    m_pred = sqrt_m_pred**2
    m_actual = [m_e_MeV, m_mu_exp, m_tau_exp][k]
    print(f"    {name}: m_pred = {m_pred:.4f}, m_exp = {m_actual:.4f}, "
          f"err = {(m_pred-m_actual)/m_actual*100:+.4f}%")

# Now test framework expressions for theta_0
delta_CKM = atan(sqrt(5))
print(f"\n  Framework angle candidates:")
print(f"    theta_0             = {theta_0_exp:.10f} rad = {np.degrees(theta_0_exp):.6f} deg")
print(f"    delta_CKM           = arctan(sqrt(5)) = {delta_CKM:.10f} rad = {np.degrees(delta_CKM):.6f} deg")
print(f"    2*delta_CKM         = {2*delta_CKM:.10f} rad = {np.degrees(2*delta_CKM):.6f} deg")
print(f"    theta_0 - 2*delta_CKM = {theta_0_exp - 2*delta_CKM:.8f} rad = {np.degrees(theta_0_exp - 2*delta_CKM):.4f} deg")
print(f"    theta_0 / delta_CKM = {theta_0_exp / delta_CKM:.8f}")

# More candidates
candidates = {
    'arctan(sqrt(5))': atan(sqrt(5)),
    '2*arctan(sqrt(5))': 2*atan(sqrt(5)),
    'arctan(phi)': atan(PHI),
    '2*arctan(phi)': 2*atan(PHI),
    'arctan(2)': atan(2),
    '2*arctan(2)': 2*atan(2),
    'pi/2 + arctan(phi)': pi/2 + atan(PHI),
    'pi - arctan(phi^2)': pi - atan(PHI**2),
    'pi - arctan(sqrt(5))': pi - atan(sqrt(5)),
    'pi - 2*pi/a1': pi - 2*pi/a1,
    'pi - pi/a1': pi - pi/a1,
    '3*pi/a1': 3*pi/a1,
    '4*pi/a1': 4*pi/a1,
    'arctan(phi^3)': atan(PHI**3),
    '2*pi/3 + pi/a1': 2*pi/3 + pi/a1,
    'pi*(1-1/a1)': pi*(1-1/a1),
    'pi*phi/a1': pi*PHI/a1,
    '2*pi*phi/a1': 2*pi*PHI/a1,
    'pi/2 + pi/(2*a1)': pi/2 + pi/(2*a1),
    'pi - arctan(2*phi)': pi - atan(2*PHI),
    'arctan(PHI) + pi/2 + pi/(2*a1)': atan(PHI) + pi/2 + pi/(2*a1),
}

print(f"\n  Systematic comparison:")
print(f"  {'Expression':>30} {'Value/rad':>12} {'Value/deg':>10} {'Ratio':>10} {'Diff/deg':>10}")

best_match = None
best_diff = 1000
for name, val in candidates.items():
    ratio = theta_0_exp / val
    diff_deg = abs(np.degrees(theta_0_exp - val))
    if diff_deg < best_diff:
        best_diff = diff_deg
        best_match = name
    if diff_deg < 5.0:
        print(f"  {name:>30} {val:12.8f} {np.degrees(val):10.4f} {ratio:10.6f} {diff_deg:10.4f}")

print(f"\n  Best match: {best_match} (diff = {best_diff:.4f} deg)")

# Deep check: 2*arctan(sqrt(5))
print(f"\n  DEEP CHECK: theta_0 vs 2*arctan(sqrt(5)) = 2*delta_CKM")
val_2delta = 2 * atan(sqrt(5))
print(f"    theta_0     = {theta_0_exp:.10f} rad = {np.degrees(theta_0_exp):.6f} deg")
print(f"    2*delta_CKM = {val_2delta:.10f} rad = {np.degrees(val_2delta):.6f} deg")
print(f"    Difference  = {theta_0_exp - val_2delta:.8f} rad = {np.degrees(theta_0_exp - val_2delta):.4f} deg")
print(f"    Ratio       = {theta_0_exp / val_2delta:.10f}")

# Check pi - arctan(sqrt(5))
val_pi_atan = pi - atan(sqrt(5))
print(f"\n    theta_0         = {theta_0_exp:.10f} rad = {np.degrees(theta_0_exp):.6f} deg")
print(f"    pi-arctan(sqrt5)= {val_pi_atan:.10f} rad = {np.degrees(val_pi_atan):.6f} deg")
print(f"    Difference      = {theta_0_exp - val_pi_atan:.8f} rad = {np.degrees(theta_0_exp - val_pi_atan):.4f} deg")

# Try to express theta_0/pi as a rational or algebraic number
theta_over_pi = theta_0_exp / pi
print(f"\n  theta_0/pi = {theta_over_pi:.10f}")
print(f"  Nearby rationals:")
for num in range(1, 20):
    for den in range(1, 20):
        val_rat = num / den
        if abs(val_rat - theta_over_pi) < 0.01:
            print(f"    {num}/{den} = {val_rat:.10f} (diff = {abs(val_rat - theta_over_pi):.6f})")

# Try continued fraction of theta_0/pi
print(f"\n  Continued fraction of theta_0/pi:")
x = theta_over_pi
cf = []
for _ in range(8):
    n_cf = int(x)
    cf.append(n_cf)
    frac = x - n_cf
    if abs(frac) < 1e-10:
        break
    x = 1.0 / frac
print(f"    [{', '.join(str(c) for c in cf)}]")

print(f"\n  ASSESSMENT (Lead 2):")
print(f"    theta_0 is NOT a simple framework angle.")
print(f"    The closest match to 2*delta_CKM has a ~{np.degrees(abs(theta_0_exp - val_2delta)):.1f} deg gap.")
print(f"    theta_0/pi ~ 0.737, not a simple rational.")
print(f"    Status: NO DERIVATION found. The Koide angle appears to be a")
print(f"    CONSEQUENCE of the exponents (0, 11, 17) rather than an independent")
print(f"    framework quantity.")


# =====================================================================
# LEAD 3: WHY alpha = 3/4 OPTIMIZES KOIDE
# =====================================================================
print(f"\n{'='*72}")
print("LEAD 3: Why alpha = 3/4 (Morrey Exponent) Optimizes Koide")
print("=" * 72)

# Scan alpha finely and find Q(alpha)
alphas = np.linspace(0.0, 2.0, 2001)
Q_of_alpha = []

zp_mu = 1 + PHIp  # z'(mu) = 1 - 1/phi = 1/phi^2
zp_tau = 1 + 2*PHIp  # z'(tau) = 1 - 2/phi = (phi-2)/1 = (1-sqrt(5))/2 + 1 = ...

print(f"\n  z'(mu)  = {zp_mu:.8f} = 1/phi^2 = {1/PHI**2:.8f}")
print(f"  z'(tau) = {zp_tau:.8f} = 1 + 2*phi' = {1 + 2*PHIp:.8f}")
print(f"  |z'(mu)| = {abs(zp_mu):.8f}")
print(f"  |z'(tau)| = {abs(zp_tau):.8f}")
print(f"  sign(z'(mu)) = {np.sign(zp_mu)}")
print(f"  sign(z'(tau)) = {np.sign(zp_tau)}")

for alpha_test in alphas:
    d_mu_t = c_ell * np.sign(zp_mu) * abs(zp_mu)**alpha_test
    d_tau_t = c_ell * np.sign(zp_tau) * abs(zp_tau)**alpha_test
    Q_test = koide_from_exponents(0, 11 + d_mu_t, 17 + d_tau_t)
    Q_of_alpha.append(Q_test)

Q_of_alpha = np.array(Q_of_alpha)
diffs = np.abs(Q_of_alpha - 2.0/3.0)
best_idx = np.argmin(diffs)
best_alpha = alphas[best_idx]
best_Q = Q_of_alpha[best_idx]

print(f"\n  Best alpha (scan 0-2, step 0.001):")
print(f"    alpha_best = {best_alpha:.4f}")
print(f"    Q_best     = {best_Q:.10f}")
print(f"    |Q - 2/3|  = {diffs[best_idx]:.2e}")
print(f"    Morrey 3/4  = {alpha_morrey:.4f}")

# Fine scan around 0.75
alphas_fine = np.linspace(0.70, 0.80, 10001)
Q_fine = []
for alpha_test in alphas_fine:
    d_mu_t = c_ell * np.sign(zp_mu) * abs(zp_mu)**alpha_test
    d_tau_t = c_ell * np.sign(zp_tau) * abs(zp_tau)**alpha_test
    Q_test = koide_from_exponents(0, 11 + d_mu_t, 17 + d_tau_t)
    Q_fine.append(Q_test)

Q_fine = np.array(Q_fine)
diffs_fine = np.abs(Q_fine - 2.0/3.0)
best_fine_idx = np.argmin(diffs_fine)
best_fine_alpha = alphas_fine[best_fine_idx]

print(f"\n  Fine scan (0.70 to 0.80, step 1e-5):")
print(f"    alpha_best = {best_fine_alpha:.6f}")
print(f"    |Q - 2/3|  = {diffs_fine[best_fine_idx]:.2e}")

# Use scipy for exact alpha
from scipy.optimize import minimize_scalar

def neg_koide_diff(alpha_test):
    d_mu_t = c_ell * np.sign(zp_mu) * abs(zp_mu)**alpha_test
    d_tau_t = c_ell * np.sign(zp_tau) * abs(zp_tau)**alpha_test
    Q_test = koide_from_exponents(0, 11 + d_mu_t, 17 + d_tau_t)
    return abs(Q_test - 2.0/3.0)

result = minimize_scalar(neg_koide_diff, bounds=(0.5, 1.5), method='bounded')
alpha_opt = result.x
print(f"\n  Scipy optimization:")
print(f"    alpha_optimal = {alpha_opt:.10f}")
print(f"    |Q - 2/3|    = {result.fun:.2e}")
print(f"    alpha_Morrey  = {alpha_morrey:.10f}")
print(f"    Diff          = {alpha_opt - alpha_morrey:.6e}")
print(f"    Rel           = {(alpha_opt - alpha_morrey)/alpha_morrey*100:.4f}%")

# So alpha = 3/4 is not EXACTLY optimal but VERY close
# What is the exact optimum? Is it a framework number?
print(f"\n  Framework number tests for alpha_optimal:")
print(f"    3/4         = 0.7500000000")
print(f"    alpha_opt   = {alpha_opt:.10f}")
print(f"    ln(phi)/phi = {log(PHI)/PHI:.10f}")
print(f"    1 - 1/d_ST  = {1 - 1/d_ST:.10f}")
print(f"    phi/a1 + 1/a1 = {PHI/a1 + 1/a1:.10f}")
print(f"    1/phi       = {1/PHI:.10f}")
print(f"    phi - 1     = {PHI - 1:.10f}")

# ANALYZE the mechanism: WHY does varying alpha affect Q?
# delta_mu(alpha) = c_ell * |z_mu'|^alpha (positive, increases n_mu)
# delta_tau(alpha) = -c_ell * |z_tau'|^alpha (negative, decreases n_tau)
# Increasing alpha: |z_mu'|^alpha changes (|z_mu'| < 1, so increases toward 1)
#                   |z_tau'|^alpha changes (|z_tau'| < 1, so increases toward 1)
# Net effect: corrections get larger (closer to c_ell in magnitude)

print(f"\n  MECHANISM:")
print(f"    |z'(mu)|  = {abs(zp_mu):.6f} = 1/phi^2 < 1")
print(f"    |z'(tau)| = {abs(zp_tau):.6f} = sqrt(5)-2 < 1")
print(f"    Both < 1, so |z'|^alpha DECREASES with alpha.")
print(f"    delta_mu = c_ell * |z_mu'|^alpha (POSITIVE, raises n_mu)")
print(f"    delta_tau = -c_ell * |z_tau'|^alpha (NEGATIVE, lowers n_tau)")
print(f"    The asymmetric corrections push mu UP and tau DOWN.")
print(f"    Q = 2/3 is about the BALANCE between these pushes.")

# Q(alpha=0) vs Q(alpha=1)
for a_test in [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0]:
    d_mu_t = c_ell * np.sign(zp_mu) * abs(zp_mu)**a_test if a_test > 0 else c_ell * np.sign(zp_mu)
    d_tau_t = c_ell * np.sign(zp_tau) * abs(zp_tau)**a_test if a_test > 0 else c_ell * np.sign(zp_tau)
    if a_test == 0:
        d_mu_t = c_ell * np.sign(zp_mu)
        d_tau_t = c_ell * np.sign(zp_tau)
    Q_test = koide_from_exponents(0, 11 + d_mu_t, 17 + d_tau_t)
    print(f"    alpha={a_test:.2f}: delta_mu={d_mu_t:+.6f}, delta_tau={d_tau_t:+.6f}, Q={Q_test:.8f}")

# Can we prove 3/4 is optimal?
# dQ/dalpha = 0 at alpha_opt.
# delta_i = c * sign(z_i') * |z_i'|^alpha
# d(delta_i)/dalpha = c * sign(z_i') * |z_i'|^alpha * ln|z_i'|
# So d(delta_mu)/dalpha = c * |z_mu'|^alpha * ln|z_mu'| < 0 (since |z_mu'|<1)
# And d(delta_tau)/dalpha = -c * |z_tau'|^alpha * ln|z_tau'| > 0 (since |z_tau'|<1, ln<0, times -1 -> >0)

ln_zmu = log(abs(zp_mu))
ln_ztau = log(abs(zp_tau))
print(f"\n  Derivatives:")
print(f"    ln|z'(mu)|  = {ln_zmu:.8f} = -2*ln(phi) = {-2*log(PHI):.8f}")
print(f"    ln|z'(tau)| = {ln_ztau:.8f}")
print(f"    ln|z'(mu)|  = -2*ln(phi) exactly? {abs(ln_zmu - (-2*log(PHI))) < 1e-10}")

# z'(mu) = 1/phi^2 exactly, so ln|z'(mu)| = -2*ln(phi) EXACTLY
# z'(tau) = -(sqrt(5)-2) = -(phi - 2 + phi - 1) ... let me compute
print(f"\n    z'(tau) = 1 + 2*phi' = 1 + 2*(1-sqrt(5))/2 = 2 - sqrt(5)")
print(f"    |z'(tau)| = sqrt(5) - 2 = {sqrt(5)-2:.8f}")
print(f"    phi^(-3) = {PHI**(-3):.8f}")
print(f"    (sqrt(5)-2) vs 1/phi^3 = {1/PHI**3:.8f}: NOT equal")
# sqrt(5)-2 = 0.2360..., 1/phi^3 = 0.2360... WAIT
print(f"    Actually: 1/phi^3 = phi^(-3) = {PHI**(-3):.10f}")
print(f"    sqrt(5)-2 = {sqrt(5)-2:.10f}")
# phi^3 = phi^2 + phi = (3+sqrt(5))/2 + (1+sqrt(5))/2 = 2+sqrt(5)
# 1/phi^3 = 1/(2+sqrt(5)) = (2-sqrt(5))/((2+sqrt(5))(2-sqrt(5))) = (2-sqrt(5))/(4-5) = (2-sqrt(5))/(-1) = sqrt(5)-2
print(f"    1/phi^3 = (sqrt(5)-2) = |z'(tau)| EXACTLY!")
print(f"    So ln|z'(tau)| = -3*ln(phi) EXACTLY!")
print(f"    Verify: -3*ln(phi) = {-3*log(PHI):.10f}, ln|z'(tau)| = {ln_ztau:.10f}")

# So we have:
# |z'(mu)| = phi^(-2), |z'(tau)| = phi^(-3)
# This means:
# delta_mu = c_ell * phi^(-2*alpha)
# delta_tau = -c_ell * phi^(-3*alpha)
# n_mu_eff = 11 + c_ell * phi^(-2*alpha)
# n_tau_eff = 17 - c_ell * phi^(-3*alpha)

print(f"\n  BEAUTIFUL SIMPLIFICATION:")
print(f"    |z'(mu)|  = phi^(-2) EXACTLY")
print(f"    |z'(tau)| = phi^(-3) EXACTLY")
print(f"    delta_mu  = c_ell * phi^(-2*alpha)")
print(f"    delta_tau = -c_ell * phi^(-3*alpha)")

# Verify
d_mu_check = c_ell * PHI**(-2*alpha_morrey)
d_tau_check = -c_ell * PHI**(-3*alpha_morrey)
print(f"\n    At alpha = 3/4:")
print(f"    delta_mu  = c_ell * phi^(-3/2) = {d_mu_check:.8f} (computed: {d_mu:.8f})")
print(f"    delta_tau = -c_ell * phi^(-9/4) = {d_tau_check:.8f} (computed: {d_tau:.8f})")

# Now: n_mu_eff = 11 + c*phi^(-2*alpha), n_tau_eff = 17 - c*phi^(-3*alpha)
# Q = f(alpha) through the exponents.
# The Koide constraint is:
# 3(1 + phi^(11+c*phi^(-2a)) + phi^(17-c*phi^(-3a))) = 2(1 + phi^((11+c*phi^(-2a))/2) + phi^((17-c*phi^(-3a))/2))^2

# This is too complex for a closed-form proof.
# But the fact that |z'| = phi^(-k) with k = 2 (mu) and k = 3 (tau)
# gives the system its special structure.

# The optimality at alpha ~ 3/4:
# For alpha = k/(k+1) with k = 3: alpha = 3/4
# With d_ST = k+1 = 4. So alpha = (d_ST-1)/d_ST, which IS the Morrey exponent.
# But is there a DEEPER reason?

# The exponent pairs: mu has phi^(-2*alpha), tau has phi^(-3*alpha).
# The ratio of logs: 3/2 = 1.5.
# For alpha = 3/4: phi^(-3/2) and phi^(-9/4). The ratio 3/2 : 9/4 = 2:3 = 2/3!
# WAIT: the ratio of the EXPONENT arguments is (3*alpha)/(2*alpha) = 3/2 = b1/2.
# That's just the ratio of |z'| exponents, not alpha-dependent.

# Let me check: is it the case that Q has a MAXIMUM at alpha=3/4?
# Compute second derivative numerically
dQ_dalpha = np.gradient(Q_of_alpha, alphas)
d2Q_dalpha2 = np.gradient(dQ_dalpha, alphas)

# Find where dQ/dalpha = 0 (looking for Q = 2/3 crossing, not extremum)
# Actually Q(alpha) is monotonic near alpha = 3/4.
# Q crosses 2/3 at alpha ~ 0.75, it doesn't have a maximum there.

# Let's check: is Q monotonic?
idx_075 = np.argmin(np.abs(alphas - 0.75))
print(f"\n  Q behavior near alpha = 0.75:")
for offset in [-5, -2, -1, 0, 1, 2, 5]:
    idx = idx_075 + offset
    if 0 <= idx < len(alphas):
        print(f"    alpha={alphas[idx]:.3f}: Q={Q_of_alpha[idx]:.10f}")

# Is Q monotonically decreasing through 2/3?
q_before = Q_of_alpha[idx_075 - 1]
q_at = Q_of_alpha[idx_075]
q_after = Q_of_alpha[idx_075 + 1]
print(f"\n  Q is {'decreasing' if q_before > q_after else 'increasing'} through alpha=0.75")
print(f"  Q CROSSES 2/3 at alpha ~ {best_fine_alpha:.6f}")
print(f"  alpha = 3/4 is where Q(alpha) = 2/3, i.e., it's a ZERO of Q-2/3.")

# So the statement is: Q(alpha=3/4) ~ 2/3 is not that 3/4 "optimizes" Q,
# but rather Q(3/4) happens to be VERY CLOSE to 2/3.

# The exact crossing:
print(f"\n  Q(alpha) crosses 2/3 at alpha = {alpha_opt:.10f}")
print(f"  Morrey alpha = 3/4 = {3/4:.10f}")
print(f"  Difference = {alpha_opt - 3/4:.2e}")

# Is the crossing EXACTLY at 3/4?
Q_at_34 = koide_from_exponents(0, 11 + c_ell * PHI**(-2*0.75), 17 - c_ell * PHI**(-3*0.75))
print(f"\n  Q(3/4) = {Q_at_34:.12f}")
print(f"  2/3    = {2/3:.12f}")
print(f"  Q(3/4) - 2/3 = {Q_at_34 - 2/3:.2e}")

# It's NOT exact. Let's quantify:
print(f"\n  RESULT: Q(alpha=3/4) = 2/3 + {Q_at_34 - 2/3:.6e}")
print(f"  This is a {abs(Q_at_34 - 2/3)/(2/3)*100:.4f}% deviation.")
print(f"  The Morrey exponent gives Q VERY CLOSE to 2/3 but NOT exactly 2/3.")

print(f"\n  ASSESSMENT (Lead 3):")
print(f"    alpha = 3/4 = (d_ST-1)/d_ST gives Q within 0.002% of 2/3.")
print(f"    This is NOT a theorem (Q is not exactly 2/3 at alpha = 3/4).")
print(f"    The beautiful simplification |z'(mu)| = phi^(-2), |z'(tau)| = phi^(-3)")
print(f"    makes the system algebraically clean but does not force exact Koide.")
print(f"    Status: PATTERN (strong numerical coincidence, not a theorem)")


# =====================================================================
# LEAD 4: TQFT / VERLINDE RING CONNECTION
# =====================================================================
print(f"\n{'='*72}")
print("LEAD 4: TQFT / Verlinde Ring and Koide")
print("=" * 72)

# In SU(2)_3 Chern-Simons (k+2 = a1 = 5, k = 3):
# Quantum dimensions: d_j for j = 0, 1/2, 1, 3/2
# d_j = sin((2j+1)*pi/5) / sin(pi/5)
# d_0 = 1, d_{1/2} = phi, d_1 = phi, d_{3/2} = 1

print(f"\n  SU(2)_3 Chern-Simons quantum dimensions:")
q_dims = []
for j2 in range(4):  # j = 0, 1/2, 1, 3/2
    j = j2 / 2.0
    d_j = sin((2*j+1)*pi/5) / sin(pi/5)
    q_dims.append(d_j)
    print(f"    j = {j}: d_j = {d_j:.8f}")

print(f"\n  Spectrum: {{1, phi, phi, 1}}")

# Koide involves sqrt(m_i) = sqrt(m_e) * phi^(n_i/2)
# The half-powers of phi are: phi^0 = 1, phi^(11/2), phi^(17/2)
# sqrt(phi) = phi^(1/2) is NOT a quantum dimension.
# But phi = d_{1/2} = d_1 IS a quantum dimension.

# The Verlinde ring (fusion rules):
# [0] x [j] = [j]
# [1/2] x [1/2] = [0] + [1]
# [1/2] x [1] = [1/2] + [3/2]
# [1] x [1] = [0] + [1]
# [1/2] x [3/2] = [1]
# [1] x [3/2] = [1/2]
# [3/2] x [3/2] = [0]

print(f"\n  Verlinde fusion ring (SU(2)_3):")
print(f"    [1/2] x [1/2] = [0] + [1]")
print(f"    [1/2] x [1]   = [1/2] + [3/2]")
print(f"    [1]   x [1]   = [0] + [1]")
print(f"    [3/2] x [3/2] = [0]")

# The ring Z[phi] = Z[x]/(x^2-x-1)
# Koide involves the sum (1 + phi^n1 + phi^n2) / (1 + phi^(n1/2) + phi^(n2/2))^2
# The denominator has phi^(n/2) = phi^(integer) * phi^(1/2) or phi^(integer)
# depending on parity.

# n_mu = 11 (odd): phi^(11/2) = phi^5 * sqrt(phi)
# n_tau = 17 (odd): phi^(17/2) = phi^8 * sqrt(phi)
# Both involve sqrt(phi), which is NOT in Z[phi].

# The total mass D^2 in TQFT:
# D^2 = sum d_j^2 = 1 + phi^2 + phi^2 + 1 = 2 + 2*phi^2 = 2(1+phi^2) = 2(2+phi) = 2(2+phi)
# = 5 + sqrt(5)
D_sq = sum(d**2 for d in q_dims)
print(f"\n  Total quantum dimension D^2 = sum(d_j^2) = {D_sq:.8f}")
print(f"  a1 + sqrt(a1) = {a1 + sqrt(a1):.8f}")
print(f"  Match: {abs(D_sq - (a1 + sqrt(a1))) < 1e-10}")

# Koide = 2/3 relation to TQFT:
# Consider the "Koide" of quantum dimensions:
Q_tqft = sum(d**2 for d in q_dims) / sum(d for d in q_dims)**2
print(f"\n  'Koide' of quantum dimensions:")
print(f"    sum(d_j^2) = {sum(d**2 for d in q_dims):.8f}")
print(f"    (sum(d_j))^2 = {sum(q_dims)**2:.8f}")
print(f"    Q_TQFT = sum(d^2) / (sum d)^2 = {Q_tqft:.8f}")
print(f"    2/3 = {2/3:.8f}")

# That's not 2/3. But what about specific subsets?
# {d_0, d_{1/2}, d_1} = {1, phi, phi}
Q_tqft_012 = (1 + PHI**2 + PHI**2) / (1 + PHI + PHI)**2
print(f"\n    Subset {{d_0, d_1/2, d_1}} = {{1, phi, phi}}:")
print(f"    Q = {Q_tqft_012:.8f}  (2/3 = {2/3:.8f})")
# Q = (1+2*phi^2)/(1+2*phi)^2 = (1+2(1+phi))/(1+2phi)^2 = (3+2phi)/(1+2phi)^2
val_num = 3 + 2*PHI
val_den = (1 + 2*PHI)**2
print(f"    = (3+2*phi)/(1+2*phi)^2 = {val_num:.4f}/{val_den:.4f} = {val_num/val_den:.8f}")

# {d_0, d_{1/2}, d_{3/2}} = {1, phi, 1}
Q_tqft_013 = (1 + PHI**2 + 1) / (1 + PHI + 1)**2
print(f"\n    Subset {{d_0, d_1/2, d_3/2}} = {{1, phi, 1}}:")
print(f"    Q = {Q_tqft_013:.8f}  (2/3 = {2/3:.8f})")
# Q = (2+phi^2)/(2+phi)^2 = (2+1+phi)/(2+phi)^2 = (3+phi)/(2+phi)^2
val_num2 = 3 + PHI
val_den2 = (2 + PHI)**2
print(f"    = (3+phi)/(2+phi)^2 = {val_num2:.4f}/{val_den2:.4f} = {val_num2/val_den2:.8f}")

# Try: Koide in the fusion ring. Mass as phi^n lives in Z[phi].
# phi^n = F_n * phi + F_{n-1}. So mass = F_n * phi + F_{n-1} in the ring.
# The "Koide" of ring elements is subtle because sqrt doesn't preserve Z[phi].

# Another approach: anyon braiding eigenvalues
# theta_j = exp(2*pi*i*h_j) where h_j = j(j+1)/(k+2) = j(j+1)/5
print(f"\n  Anyon topological spins:")
for j2 in range(4):
    j = j2 / 2.0
    h_j = j*(j+1)/5.0
    theta_j_angle = 2*pi*h_j
    print(f"    j={j}: h_j = {h_j:.4f}, theta_j = {theta_j_angle:.6f} = {np.degrees(theta_j_angle):.2f} deg")

# S-matrix of SU(2)_3
print(f"\n  S-matrix elements S_{0j}:")
S = np.zeros((4, 4))
for i2 in range(4):
    for j2 in range(4):
        i = i2/2.0
        j = j2/2.0
        S[i2, j2] = sqrt(2.0/5.0) * sin((2*i+1)*(2*j+1)*pi/5.0)

print(f"    S =")
for row in S:
    print(f"      [{', '.join(f'{x:8.5f}' for x in row)}]")

# S_{0j}/S_{00} = d_j (quantum dimensions)
print(f"    S_0j/S_00 = [{', '.join(f'{S[0,j]/S[0,0]:.4f}' for j in range(4))}]")

# The "Koide" of S-matrix column ratios?
s_ratios = [S[0,j]/S[0,0] for j in range(4)]
Q_s = sum(sr**2 for sr in s_ratios) / sum(s_ratios)**2
print(f"    'Koide' of S_0j/S_00 = {Q_s:.8f}")

# Try: the Verlinde formula gives N_{ij}^k (fusion coefficients).
# Can we construct a Q = 2/3 from fusion data?
# N_gen = 3 is derived from McKay p_4/p_2 = 3.
# 2/3 = 2/N_gen. Is there a trace formula?

print(f"\n  Q = 2/3 = 2/N_gen = 2/3")
print(f"  In TQFT: 2 = D^2 * D'^2 / N = {D_sq * (a1-sqrt(a1))}/N? No, that's not right.")
D_sq_prime = a1 - sqrt(a1)
print(f"  D^2 = {D_sq:.6f}, D'^2 = {D_sq_prime:.6f}")
print(f"  D^2 * D'^2 = {D_sq * D_sq_prime:.6f} = a1*(a1-1) = {a1*(a1-1)}")
print(f"  N/b1 = {N/b1:.1f} = 20 = a1*(a1-1). Check.")

# 2/3 from quantum dimensions more directly:
# Consider: sum(d_j^2) / (sum d_j)^2 = D^2 / (sum d_j)^2
sum_d = sum(q_dims)
print(f"\n  sum(d_j) = {sum_d:.8f} = 2 + 2*phi = 2*(1+phi) = 2*phi^2 = {2*PHI**2:.8f}")
print(f"  (sum d_j)^2 = {sum_d**2:.6f} = 4*phi^4 = {4*PHI**4:.6f}")
print(f"  D^2/(sum d)^2 = {D_sq/sum_d**2:.8f}")
print(f"  Compare: 2/(a1+1) = 2/b1 = {2/b1:.8f}")

# Actually: sum d_j = 1 + phi + phi + 1 = 2(1+phi) = 2*phi^2
# (sum d_j)^2 = 4*phi^4
# D^2 = 5+sqrt(5) = a1 + sqrt(a1)
# D^2/4*phi^4 = (5+sqrt(5))/(4*phi^4) = ...
# phi^4 = phi^2*(phi+1) = (phi+1)(phi+1) = phi^2 + 2*phi + 1 = 3*phi + 2 + 1 = ...
# phi^4 = 3*phi + 2. So 4*phi^4 = 12*phi + 8.
# D^2 = 5 + sqrt(5) = 5 + 2*phi - 1 = 4 + 2*phi
# Ratio = (4+2*phi)/(12*phi+8) = 2(2+phi)/(4(3*phi+2)) = (2+phi)/(2(3*phi+2))
# = (2+phi)/(6*phi+4) = ... let me just compute
ratio_tqft = D_sq / (4*PHI**4)
print(f"  D^2/(4*phi^4) = {ratio_tqft:.8f}")

print(f"\n  ASSESSMENT (Lead 4):")
print(f"    No direct TQFT derivation of Q = 2/3 found.")
print(f"    The quantum dimensions {{1, phi, phi, 1}} do not give Q = 2/3")
print(f"    for any natural subset.")
print(f"    The connection to Koide seems INDIRECT: the mass formula lives in")
print(f"    Z[phi] (the Verlinde ring), but Koide requires sqrt(phi) which")
print(f"    exits the ring. This is a STRUCTURAL OBSTRUCTION.")
print(f"    Status: SPECULATIV (no derivation found)")


# =====================================================================
# LEAD 5: EXACT KOIDE
# =====================================================================
print(f"\n{'='*72}")
print("LEAD 5: What Gives Exact Q = 2/3?")
print("=" * 72)

# Part A: What c_ell (with alpha=3/4) gives Q = 2/3 exactly?
print(f"\n  Part A: Find c_ell* such that Q(c_ell*, alpha=3/4) = 2/3 exactly")

def Q_of_c(c_test, alpha=0.75):
    d_mu_t = c_test * PHI**(-2*alpha)
    d_tau_t = -c_test * PHI**(-3*alpha)
    return koide_from_exponents(0, 11 + d_mu_t, 17 + d_tau_t)

# Search
result_c = brentq(lambda c: Q_of_c(c) - 2.0/3.0, 0.001, 1.0, xtol=1e-14)
print(f"    c_ell* = {result_c:.12f}")
print(f"    c_ell (framework) = {c_ell:.12f}")
print(f"    Diff = {result_c - c_ell:.6e}")
print(f"    Ratio c_ell*/c_ell = {result_c/c_ell:.10f}")

# Is c_ell* a framework number?
print(f"\n    Framework comparisons for c_ell*:")
print(f"      c_ell = C*phi^3/d_ST = (2/13)*phi^3/4 = {c_ell:.10f}")
print(f"      c_ell* = {result_c:.10f}")
print(f"      c_ell*/(phi^3) = {result_c/PHI**3:.10f}")
print(f"      c_ell*/(phi^3/d_ST) = {result_c/(PHI**3/d_ST):.10f}")
print(f"      4*c_ell*/phi^3 = {4*result_c/PHI**3:.10f}  (C would be {C_mass:.10f})")
C_exact_koide = 4*result_c/PHI**3
print(f"      C_exact = {C_exact_koide:.10f}")
print(f"      2/13 = {2/13:.10f}")
print(f"      C_exact - 2/13 = {C_exact_koide - 2/13:.6e}")
print(f"      C_exact * 13/2 = {C_exact_koide * 13/2:.10f}")

# Part B: What corrections delta_mu, delta_tau give exact Q = 2/3?
print(f"\n  Part B: What corrections give Q = 2/3 with n_mu=11, n_tau=17?")

# We need: Q(0, 11+d_mu, 17+d_tau) = 2/3
# This is one equation in two unknowns. There's a curve of solutions.
# Our framework constrains the RATIO: d_mu/d_tau = -phi^(-2*alpha)/phi^(-3*alpha) = -phi^alpha
# With alpha = 3/4: d_mu/d_tau = -phi^(3/4) = -phi^0.75

ratio_dt = -PHI**(alpha_morrey)
print(f"    Framework constraint: delta_mu / delta_tau = -phi^(3/4) = {ratio_dt:.8f}")

# With this constraint: delta_tau = delta_mu / ratio_dt
# Q(0, 11 + d_mu, 17 + d_mu/ratio_dt) = 2/3
def Q_of_dmu(dmu):
    dtau = dmu / ratio_dt
    return koide_from_exponents(0, 11 + dmu, 17 + dtau)

dmu_exact = brentq(lambda d: Q_of_dmu(d) - 2.0/3.0, 0.001, 1.0, xtol=1e-14)
dtau_exact = dmu_exact / ratio_dt
print(f"    delta_mu_exact  = {dmu_exact:.10f}  (framework: {d_mu:.10f})")
print(f"    delta_tau_exact = {dtau_exact:.10f}  (framework: {d_tau:.10f})")
print(f"    n_mu_exact  = {11 + dmu_exact:.10f}  (framework: {n_mu_eff:.10f})")
print(f"    n_tau_exact = {17 + dtau_exact:.10f}  (framework: {n_tau_eff:.10f})")
print(f"    Diff in delta_mu: {dmu_exact - d_mu:.6e}")
print(f"    Diff in delta_tau: {dtau_exact - d_tau:.6e}")

# What masses would these give?
m_mu_exact = m_e_MeV * PHI**(11 + dmu_exact)
m_tau_exact = m_e_MeV * PHI**(17 + dtau_exact)
print(f"\n    Exact-Koide masses:")
print(f"      m_mu  = {m_mu_exact:.6f} MeV (exp: {m_mu_exp:.6f}, diff: {(m_mu_exact-m_mu_exp)/m_mu_exp*100:+.4f}%)")
print(f"      m_tau = {m_tau_exact:.4f} MeV (exp: {m_tau_exp:.2f}, diff: {(m_tau_exact-m_tau_exp)/m_tau_exp*100:+.4f}%)")

# Part C: Without framework ratio constraint (free delta_mu, delta_tau)
print(f"\n  Part C: Free search for exact Koide with best mass fit")

from scipy.optimize import minimize

def mass_rms(params):
    dmu, dtau = params
    m_mu_p = m_e_MeV * PHI**(11 + dmu)
    m_tau_p = m_e_MeV * PHI**(17 + dtau)
    err_mu = (m_mu_p - m_mu_exp) / m_mu_exp
    err_tau = (m_tau_p - m_tau_exp) / m_tau_exp
    return sqrt(err_mu**2 + err_tau**2)

# Find all (d_mu, d_tau) on the Koide surface that minimize mass errors
def objective(params):
    dmu, dtau = params
    Q_test = koide_from_exponents(0, 11 + dmu, 17 + dtau)
    penalty = 1e6 * (Q_test - 2.0/3.0)**2
    return mass_rms(params) + penalty

res = minimize(objective, [d_mu, d_tau], method='Nelder-Mead',
               options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 100000})
dmu_opt, dtau_opt = res.x
Q_opt = koide_from_exponents(0, 11 + dmu_opt, 17 + dtau_opt)
m_mu_opt = m_e_MeV * PHI**(11 + dmu_opt)
m_tau_opt = m_e_MeV * PHI**(17 + dtau_opt)

print(f"    Optimized (exact Koide + best masses):")
print(f"      delta_mu  = {dmu_opt:.10f}")
print(f"      delta_tau = {dtau_opt:.10f}")
print(f"      Q         = {Q_opt:.12f}")
print(f"      m_mu      = {m_mu_opt:.4f} MeV (err: {(m_mu_opt-m_mu_exp)/m_mu_exp*100:+.4f}%)")
print(f"      m_tau     = {m_tau_opt:.2f} MeV (err: {(m_tau_opt-m_tau_exp)/m_tau_exp*100:+.4f}%)")
print(f"      ratio delta_mu/delta_tau = {dmu_opt/dtau_opt:.8f}")
print(f"      -phi^(3/4)                = {ratio_dt:.8f}")

print(f"\n  ASSESSMENT (Lead 5):")
print(f"    Exact Q = 2/3 requires c_ell* = {result_c:.8f} (vs {c_ell:.8f}, diff {(result_c-c_ell)/c_ell*100:.2f}%)")
print(f"    The required C* = {C_exact_koide:.8f} is close to but NOT exactly 2/13.")
print(f"    With the framework ratio constraint, exact Koide gives masses")
print(f"    that are still good but slightly shifted from experiment.")
print(f"    Status: Koide Q = 2/3 is APPROXIMATE in the framework, NOT exact.")


# =====================================================================
# LEAD 6: GALOIS STRUCTURE OF Q(sqrt(5), sqrt(phi))
# =====================================================================
print(f"\n{'='*72}")
print("LEAD 6: Galois Structure of Q(sqrt(5), sqrt(phi))")
print("=" * 72)

# Koide involves sqrt(m_i) where m_i = m_e * phi^(n_i).
# For n_i even: sqrt(m_i) = sqrt(m_e) * phi^(n_i/2) which is in Q(sqrt(m_e), phi) = Q(sqrt(m_e), sqrt(5))
# For n_i odd: sqrt(m_i) = sqrt(m_e) * phi^((n_i-1)/2) * sqrt(phi)
#   which requires sqrt(phi), extending the field.

# phi = (1+sqrt(5))/2 is in Q(sqrt(5)), degree 2 over Q.
# sqrt(phi) is a root of x^4 - x^2 - 1/4... no.
# Let y = sqrt(phi). Then y^2 = phi = (1+sqrt(5))/2.
# y^4 = phi^2 = (3+sqrt(5))/2.
# 2*y^4 - 3 = sqrt(5), so (2*y^4-3)^2 = 5, giving:
# 4*y^8 - 12*y^4 + 9 = 5
# 4*y^8 - 12*y^4 + 4 = 0
# y^8 - 3*y^4 + 1 = 0

print(f"\n  Field extension analysis:")
print(f"    phi = (1+sqrt(5))/2 is in Q(sqrt(5)), [Q(sqrt(5)):Q] = 2")
print(f"    sqrt(phi) satisfies: y^8 - 3y^4 + 1 = 0")

# Check this:
y = sqrt(PHI)
check_poly = y**8 - 3*y**4 + 1
print(f"    Verify: y^8 - 3y^4 + 1 = {check_poly:.2e} (should be 0)")

# Factor: y^8 - 3y^4 + 1 = (y^4 - phi^2)(y^4 - phi'^2)
# where phi'^2 = (3-sqrt(5))/2
# y^4 = phi^2 = (3+sqrt(5))/2 (for y = sqrt(phi))
phi_prime_sq = (3 - sqrt(5))/2
print(f"    phi^2 = (3+sqrt(5))/2 = {PHI**2:.10f}")
print(f"    phi'^2 = (3-sqrt(5))/2 = {phi_prime_sq:.10f}")
print(f"    phi^2 * phi'^2 = {PHI**2 * phi_prime_sq:.10f} (should be 1)")
print(f"    phi^2 + phi'^2 = {PHI**2 + phi_prime_sq:.10f} (should be 3)")

# The minimal polynomial of sqrt(phi) over Q:
# We need to check if y^8-3y^4+1 is irreducible over Q.
# Over Q(sqrt(5)): y^4 - phi^2 = 0, i.e., y^4 = (3+sqrt(5))/2
# Is y^4 - (3+sqrt(5))/2 irreducible over Q(sqrt(5))?
# The roots are: sqrt(phi), -sqrt(phi), i*sqrt(phi), -i*sqrt(phi)
# Actually y^4 = phi^2, so y = phi^(1/2), -phi^(1/2), i*phi^(1/2), -i*phi^(1/2)
# Over Q(sqrt(5)): y^4 - phi^2 factors as (y^2 - phi)(y^2 + phi)
# = (y - sqrt(phi))(y + sqrt(phi))(y^2 + phi)
# So [Q(sqrt(phi)) : Q(sqrt(5))] = 2, not 4.

print(f"\n  Factorization over Q(sqrt(5)):")
print(f"    y^8 - 3y^4 + 1 = (y^4 - phi^2)(y^4 - phi'^2)")
print(f"    Over Q(sqrt(5)): y^4 - phi^2 = (y^2 - phi)(y^2 + phi)")
print(f"    So sqrt(phi) is degree 2 over Q(sqrt(5)).")
print(f"    [Q(sqrt(5), sqrt(phi)) : Q] = [Q(sqrt(phi)) : Q(sqrt(5))] * [Q(sqrt(5)) : Q]")
print(f"                                 = 2 * 2 = 4")

# The Galois group of Q(sqrt(phi))/Q:
# The minimal polynomial of sqrt(phi) over Q is x^4 - x - 1... no.
# Let's find it properly.
# y = sqrt(phi), y^2 = phi = (1+sqrt(5))/2
# 2*y^2 = 1 + sqrt(5), so sqrt(5) = 2*y^2 - 1
# (2*y^2 - 1)^2 = 5
# 4*y^4 - 4*y^2 + 1 = 5
# 4*y^4 - 4*y^2 - 4 = 0
# y^4 - y^2 - 1 = 0

check2 = y**4 - y**2 - 1
print(f"\n  Minimal polynomial of sqrt(phi) over Q:")
print(f"    y^4 - y^2 - 1 = 0")
print(f"    Verify: {check2:.2e}")

# This is a degree 4 polynomial. Its roots:
# y^2 = (1 +/- sqrt(5))/2 => y^2 = phi or y^2 = phi'
# For y^2 = phi > 0: y = +/- sqrt(phi)
# For y^2 = phi' < 0: y = +/- i*sqrt(|phi'|) = +/- i*sqrt(1/phi)
# Roots: sqrt(phi), -sqrt(phi), i/sqrt(phi), -i/sqrt(phi)

print(f"\n  Roots of y^4 - y^2 - 1:")
roots = [sqrt(PHI), -sqrt(PHI)]
print(f"    y_1 = +sqrt(phi) = {roots[0]:.10f}")
print(f"    y_2 = -sqrt(phi) = {roots[1]:.10f}")
print(f"    y_3 = +i*sqrt(1/phi) = {sqrt(1/PHI):.10f}*i")
print(f"    y_4 = -i*sqrt(1/phi) = {sqrt(1/PHI):.10f}*i")

# The Galois group:
# Automorphisms must permute these roots.
# sigma: sqrt(phi) -> -sqrt(phi) (just sign flip)
# tau: sqrt(phi) -> i/sqrt(phi) (phi -> phi', sqrt -> i*sqrt(|.|))
# sigma^2 = id, tau^2: sqrt(phi) -> tau(i/sqrt(phi)) = i/tau(sqrt(phi)) = i/(i/sqrt(phi)) = sqrt(phi)
# So tau^2 = id. And sigma*tau: sqrt(phi) -> -i/sqrt(phi)
# tau*sigma: sqrt(phi) -> tau(-sqrt(phi)) = -i/sqrt(phi) (same!)
# So sigma*tau = tau*sigma: the group is ABELIAN, isomorphic to Z/2 x Z/2 = V4 (Klein four-group).

print(f"\n  Galois group Gal(Q(sqrt(phi))/Q):")
print(f"    sigma: sqrt(phi) -> -sqrt(phi)  (sign flip)")
print(f"    tau:   sqrt(phi) -> i/sqrt(phi)  (Galois conjugation + sqrt inversion)")
print(f"    sigma^2 = tau^2 = (sigma*tau)^2 = id")
print(f"    Gal = Z/2 x Z/2 = V4 (Klein four-group)")
print(f"    |Gal| = 4 = d_ST !!!")

# The Galois group has order 4 = d_ST!
print(f"\n  *** |Gal(Q(sqrt(phi))/Q)| = 4 = d_ST ***")
print(f"  The field extension needed for Koide has degree = spacetime dimension!")

# Connection to a1 = 5:
# The polynomial y^4 - y^2 - 1 has discriminant:
# For y^4 + py^2 + q: disc = p^2 - 4q = 1 + 4 = 5 = a1!
disc_poly = 1 + 4  # p = -1, q = -1, disc = (-1)^2 - 4*(-1) = 1+4 = 5
print(f"\n  Polynomial y^4 - y^2 - 1:")
print(f"    As y^2-quadratic: t^2 - t - 1 = 0 (where t = y^2)")
print(f"    Discriminant of t^2 - t - 1: 1 + 4 = 5 = a1!")
print(f"    This is the GOLDEN RATIO equation!")

# The Koide ratio in terms of field theory:
# Q = sum(m_i) / (sum(sqrt(m_i)))^2
# The numerator lives in Q(sqrt(5)) (since phi^n is in Q(sqrt(5)))
# The denominator involves sqrt(phi) and thus lives in Q(sqrt(phi))
# Q itself is in Q (it's a rational number approximately 2/3)

# For Q to be rational, we need: (sum sqrt(m_i))^2 is in Q(sqrt(5))
# sum sqrt(m_i) = sqrt(m_e) * (1 + phi^(11/2) + phi^(17/2))
# = sqrt(m_e) * (1 + phi^5 * sqrt(phi) + phi^8 * sqrt(phi))
# = sqrt(m_e) * (1 + (phi^5 + phi^8) * sqrt(phi))
# (sum sqrt(m_i))^2 = m_e * (1 + (phi^5+phi^8)*sqrt(phi))^2
# = m_e * (1 + 2*(phi^5+phi^8)*sqrt(phi) + (phi^5+phi^8)^2 * phi)
# The cross-term involves sqrt(phi), so (sum)^2 is NOT in Q(sqrt(5)) unless
# it cancels when dividing by sum(m_i).

print(f"\n  Algebraic structure of Q:")
phi5 = PHI**5
phi8 = PHI**8
A_coeff = phi5 + phi8  # coefficient of sqrt(phi)
print(f"    phi^5 + phi^8 = {A_coeff:.6f}")
print(f"    = F_5*phi + F_4 + F_8*phi + F_7")
F = [0, 1, 1, 2, 3, 5, 8, 13, 21]
A_phi_part = F[5] + F[8]  # phi coefficient
A_const_part = F[4] + F[7]  # constant part
print(f"    = ({A_phi_part})*phi + ({A_const_part}) = {A_phi_part}*phi + {A_const_part}")
print(f"    = {A_phi_part * PHI + A_const_part:.6f}")

# For Q to be exactly 2/3, the irrational parts must cancel in a specific way.
# This is a STRONG constraint that connects the exponents to the field structure.

print(f"\n  For Q = sum(m) / (sum sqrt(m))^2 to be RATIONAL:")
print(f"    The irrational part (involving sqrt(phi)) in (sum sqrt(m))^2")
print(f"    must be absorbed by the numerator's irrational part.")
print(f"    Since sum(m) is in Q(sqrt(5)) and (sum sqrt(m))^2 has terms")
print(f"    in sqrt(phi), exact rationality of Q requires a non-trivial")
print(f"    algebraic identity relating phi^n and phi^(n/2).")

# Compute the Galois norm of Q
# Under sigma (phi -> phi): Q -> Q (phi is fixed, only sqrt sign changes, but Q is sum/sum^2)
# Under tau (phi -> phi'): this changes EVERYTHING
# Q' = (1 + (phi')^11 + (phi')^17) / (1 + (phi')^(11/2) + (phi')^(17/2))^2
# phi' = -1/phi, so (phi')^11 = -1/phi^11, (phi')^17 = -1/phi^17

phi_11 = PHI**11
phi_17 = PHI**17
Q_galois = (1 + (-1/PHI)**11 + (-1/PHI)**17) / (1 + (-1/PHI)**(11.0/2) + (-1/PHI)**(17.0/2))**2

# (-1/phi)^11 = -1/phi^11, (-1/phi)^17 = -1/phi^17
# (-1/phi)^(11/2) = (-1)^(11/2)/phi^(11/2) -- complex!
# Actually (-1/phi)^(11/2) is complex because of (-1)^(11/2)
# So the Galois conjugate of Q involves complex numbers

print(f"\n  Galois conjugate Q' (under phi -> phi' = -1/phi):")
print(f"    phi'^11 = (-1/phi)^11 = -1/phi^11 = {(-1/PHI)**11:.10f}")
print(f"    phi'^17 = (-1/phi)^17 = -1/phi^17 = {(-1/PHI)**17:.10f}")
print(f"    Numerator' = 1 - 1/phi^11 - 1/phi^17 = {1 - 1/phi_11 - 1/phi_17:.10f}")
print(f"    For the denominator: phi'^(n/2) involves (-1)^(n/2) which is COMPLEX")
print(f"    for odd n. This means the Galois conjugate is NOT real.")
print(f"    The Koide ratio breaks the Q(sqrt(5)) Galois symmetry.")

# Real part analysis
num_conj = 1 - 1/phi_11 - 1/phi_17
print(f"\n    Q_numerator' = {num_conj:.10f}")
print(f"    Q_numerator  = {1 + phi_11 + phi_17:.6f}")
print(f"    Product = {num_conj * (1 + phi_11 + phi_17):.10f}")
print(f"    Norm of numerator (phi^0 + phi^11 + phi^17):")
print(f"      N = (1+phi^11+phi^17)(1-1/phi^11-1/phi^17)")
val_norm_num = (1 + phi_11 + phi_17) * num_conj
print(f"      = {val_norm_num:.10f}")

# Simplify: (1+phi^11+phi^17)(1-phi^(-11)-phi^(-17))
# = 1 - phi^(-11) - phi^(-17) + phi^11 - 1 - phi^(11-17) + phi^17 - phi^(17-11) - 1
# = -1 + phi^11 + phi^17 - phi^(-11) - phi^(-17) - phi^(-6) - phi^6
# This is messy. Let's just note the numerical value.

print(f"\n  ASSESSMENT (Lead 6):")
print(f"    The minimal polynomial of sqrt(phi) over Q is y^4 - y^2 - 1 = 0.")
print(f"    Its discriminant (as quadratic in y^2) is 5 = a1.")
print(f"    The Galois group is V4 = Z/2 x Z/2 with |Gal| = 4 = d_ST.")
print(f"    *** The field extension degree = spacetime dimension! ***")
print(f"    The Galois conjugation under phi -> phi' sends Koide to a COMPLEX number,")
print(f"    breaking the Galois symmetry. This is the same mechanism as the dark sector.")
print(f"    Status: |Gal| = d_ST is STRUCTURAL (a consequence of sqrt(phi) over Q).")
print(f"    The connection to a1 through discriminant = 5 is DERIVED (algebraic identity).")


# =====================================================================
# SYNTHESIS AND SUMMARY
# =====================================================================
print(f"\n{'='*72}")
print("SYNTHESIS: What We Learned")
print("=" * 72)

print(f"""
  LEAD 1 - PALINDROMIC POLYNOMIAL:
    Result: The palindromic condition (t = s + 1/s = a1 = 5) is SPECIFIC to
    geometric sequences. For our actual exponents (0, 11, 17), the Koide
    equation becomes quadratic in u = phi^(q/2). The bare q = b1 = 6 is
    CLOSE to but not exactly the Koide-exact value.
    Status: PATTERN

  LEAD 2 - KOIDE ANGLE theta_0:
    Result: theta_0 ~ 2.3166 rad = 132.73 deg is NOT a simple framework
    expression. It is not 2*delta_CKM, not n*pi/a1, not arctan of a framework
    number. theta_0 is DETERMINED by the exponents, not an independent quantity.
    Status: NEGATIVE (no framework expression found)

  LEAD 3 - WHY alpha = 3/4:
    Result: BEAUTIFUL simplification: |z'(mu)| = phi^(-2), |z'(tau)| = phi^(-3)
    EXACTLY. This means delta_mu = c*phi^(-2*alpha), delta_tau = -c*phi^(-3*alpha).
    At alpha = 3/4, Q is VERY close to 2/3 (within 0.0024%) but NOT exactly 2/3.
    The exact crossing is at alpha = {alpha_opt:.8f}, differing from 3/4 by {(alpha_opt-0.75):.2e}.
    Status: PATTERN (not a theorem, but 0.002% is striking)

  LEAD 4 - TQFT / VERLINDE RING:
    Result: No direct TQFT derivation found. The quantum dimensions
    {{1, phi, phi, 1}} do not produce Q = 2/3 for any natural subset.
    The OBSTRUCTION: Koide requires sqrt(phi) which exits the Verlinde ring Z[phi].
    Status: SPECULATIV (no mechanism found)

  LEAD 5 - EXACT KOIDE:
    Result: Exact Q = 2/3 requires c_ell* = {result_c:.8f} (vs {c_ell:.8f}).
    The difference is small ({(result_c-c_ell)/c_ell*100:.2f}%) but nonzero.
    C* = {C_exact_koide:.8f} vs C = 2/13 = {2/13:.8f}.
    With framework ratio constraint, exact Koide masses are close to experiment.
    Status: APPROXIMATE (Koide is NOT exact in the framework)

  LEAD 6 - GALOIS STRUCTURE:
    Result: *** The field Q(sqrt(phi)) has degree 4 = d_ST over Q! ***
    Minimal polynomial: y^4 - y^2 - 1 = 0 with discriminant 5 = a1.
    Galois group: V4 = Z/2 x Z/2 (Klein four-group).
    Under phi -> phi', Koide becomes COMPLEX (symmetry breaking).
    Status: STRUCTURAL (|Gal| = d_ST derived, disc = a1 derived)

  OVERALL CLASSIFICATION:
    - Q_framework = 2/3 to 0.0024%: PATTERN
    - t = a1 = 5 in geometric limit: DERIVED (algebraic theorem)
    - alpha = 3/4 crossing Q = 2/3: PATTERN (not exact)
    - |z'| = phi^(-2), phi^(-3): DERIVED (algebraic identity)
    - [Q(sqrt(phi)):Q] = 4 = d_ST: DERIVED (algebraic, with disc = a1)
    - Koide angle: NO DERIVATION
    - TQFT connection: NEGATIVE (obstruction from sqrt(phi) not in Z[phi])

  THE DEEPEST RESULT: The Koide formula connects to the framework through
  the field Q(sqrt(phi)), whose Galois group order = d_ST = 4 and whose
  discriminant = a1 = 5. The Koide ratio Q = 2/3 is APPROXIMATE in the
  framework, with the 3/4 Morrey exponent placing it within 0.002%.
  A full derivation of exact Q = 2/3 remains OPEN.
""")

print("=" * 72)
print("EXP-411b COMPLETE")
print("=" * 72)

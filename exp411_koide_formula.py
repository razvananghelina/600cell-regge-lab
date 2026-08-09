"""
exp411_koide_formula.py
========================
GOAL: Derive the Koide formula Q = 2/3 from the 600-cell mass framework.

The Koide formula (1981):
    Q = (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3

Experimentally: Q = 0.666661 +/- 0.000007  (vs 2/3 = 0.666667)
NO theoretical derivation exists in 40+ years.

Our framework: m_f = m_e * phi^(n_f + delta_f)
  n_e = 0,  n_mu = 11,  n_tau = 17
  delta_f = c_ell * sign(z') * |z'|^(3/4)  for leptons

PLAN:
  Part 1: Numerical check with experimental masses
  Part 2: Numerical check with framework masses (with corrections)
  Part 3: Algebraic analysis - can Q = 2/3 be DERIVED?
  Part 4: The role of the 3/4 exponent
  Part 5: Assessment

Author: Razvan-Constantin Anghelina
Date: 2026-02-25
"""

import numpy as np
from math import sqrt, pi, log, factorial

print("=" * 72)
print("EXP-411: THE KOIDE FORMULA FROM THE 600-CELL FRAMEWORK")
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

m_e_MeV = 0.51099895000  # MeV CODATA 2018
m_mu_exp = 105.6583755    # MeV
m_tau_exp = 1776.86        # MeV

# =====================================================================
# PART 1: KOIDE WITH EXPERIMENTAL MASSES
# =====================================================================
print(f"\n{'='*72}")
print("PART 1: Koide Formula with Experimental Masses")
print("=" * 72)

masses_exp = [m_e_MeV, m_mu_exp, m_tau_exp]

sum_m = sum(masses_exp)
sum_sqrt_m = sum(sqrt(m) for m in masses_exp)
Q_exp = sum_m / sum_sqrt_m**2

print(f"\n  Experimental masses:")
print(f"    m_e   = {m_e_MeV:.10f} MeV")
print(f"    m_mu  = {m_mu_exp:.7f} MeV")
print(f"    m_tau = {m_tau_exp:.2f} MeV")
print(f"\n  sum(m_i)        = {sum_m:.6f} MeV")
print(f"  sum(sqrt(m_i))  = {sum_sqrt_m:.6f} MeV^(1/2)")
print(f"  sum(sqrt(m_i))^2= {sum_sqrt_m**2:.6f} MeV")
print(f"\n  Q_exp = {Q_exp:.10f}")
print(f"  2/3   = {2/3:.10f}")
print(f"  Diff  = {abs(Q_exp - 2/3):.2e}")
print(f"  Rel   = {abs(Q_exp - 2/3)/(2/3)*100:.4f}%")


# =====================================================================
# PART 2: KOIDE WITH FRAMEWORK MASSES
# =====================================================================
print(f"\n{'='*72}")
print("PART 2: Koide Formula with Framework Masses")
print("=" * 72)

# Lepton (a,b) assignments
leptons = {
    'e':   {'a': 0, 'b': 0, 'n_bare': 0},
    'mu':  {'a': 1, 'b': 1, 'n_bare': 11},
    'tau': {'a': 1, 'b': 2, 'n_bare': 17},
}

print(f"\n  Framework parameters:")
print(f"    C = 4/(a1^2+1) = {C_mass:.8f}")
print(f"    c_ell = C*phi^3/d_ST = {c_ell:.8f}")
print(f"    alpha_Morrey = 3/4 = {alpha_morrey}")

framework_masses = {}
print(f"\n  {'Lepton':>6} {'n':>4} {'z_prime':>10} {'delta':>10} {'n_eff':>10} {'m/MeV':>12} {'m_exp/MeV':>12} {'err%':>8}")

for name in ['e', 'mu', 'tau']:
    data = leptons[name]
    a, b, n = data['a'], data['b'], data['n_bare']
    zp = a + b * PHIp

    if abs(zp) < 1e-10:
        delta = 0.0
    else:
        delta = c_ell * np.sign(zp) * abs(zp)**alpha_morrey

    n_eff = n + delta
    m_fw = m_e_MeV * PHI**n_eff
    framework_masses[name] = m_fw

    m_exp = [m_e_MeV, m_mu_exp, m_tau_exp][['e','mu','tau'].index(name)]
    err = (m_fw - m_exp) / m_exp * 100
    print(f"  {name:>6} {n:4d} {zp:10.6f} {delta:10.6f} {n_eff:10.6f} {m_fw:12.6f} {m_exp:12.6f} {err:8.4f}")

masses_fw = [framework_masses['e'], framework_masses['mu'], framework_masses['tau']]
sum_m_fw = sum(masses_fw)
sum_sqrt_fw = sum(sqrt(m) for m in masses_fw)
Q_fw = sum_m_fw / sum_sqrt_fw**2

print(f"\n  KOIDE with framework masses:")
print(f"    Q_fw  = {Q_fw:.10f}")
print(f"    2/3   = {2/3:.10f}")
print(f"    Diff  = {abs(Q_fw - 2/3):.2e}")
print(f"    Rel   = {abs(Q_fw - 2/3)/(2/3)*100:.4f}%")


# =====================================================================
# PART 3: KOIDE WITH BARE EXPONENTS (no corrections)
# =====================================================================
print(f"\n{'='*72}")
print("PART 3: Koide with Bare Exponents (delta = 0)")
print("=" * 72)

masses_bare = [m_e_MeV, m_e_MeV * PHI**11, m_e_MeV * PHI**17]
sum_m_bare = sum(masses_bare)
sum_sqrt_bare = sum(sqrt(m) for m in masses_bare)
Q_bare = sum_m_bare / sum_sqrt_bare**2

print(f"\n  Bare masses (no corrections):")
print(f"    m_e   = {masses_bare[0]:.6f} MeV")
print(f"    m_mu  = {masses_bare[1]:.3f} MeV (exp: {m_mu_exp:.3f})")
print(f"    m_tau = {masses_bare[2]:.2f} MeV (exp: {m_tau_exp:.2f})")
print(f"\n  Q_bare = {Q_bare:.10f}")
print(f"  2/3    = {2/3:.10f}")
print(f"  Diff   = {abs(Q_bare - 2/3):.2e}")
print(f"  Rel    = {abs(Q_bare - 2/3)/(2/3)*100:.4f}%")

print(f"\n  The CORRECTION brings Q from {Q_bare:.6f} to {Q_fw:.6f}")
print(f"  Movement toward 2/3: {abs(Q_bare - 2/3) - abs(Q_fw - 2/3):.2e}")
if abs(Q_fw - 2/3) < abs(Q_bare - 2/3):
    print(f"  --> Corrections IMPROVE Koide (move Q toward 2/3)")
else:
    print(f"  --> Corrections WORSEN Koide (move Q away from 2/3)")


# =====================================================================
# PART 4: ALGEBRAIC ANALYSIS
# =====================================================================
print(f"\n{'='*72}")
print("PART 4: Algebraic Structure of Koide in the Framework")
print("=" * 72)

# Write Q in terms of phi exponents
# m_f = m_e * phi^n_f, so m_f/m_e = phi^n_f
# Q = m_e * (1 + phi^n_mu + phi^n_tau) / (m_e * (1 + phi^(n_mu/2) + phi^(n_tau/2))^2)
#   = (1 + phi^n_mu + phi^n_tau) / (1 + phi^(n_mu/2) + phi^(n_tau/2))^2

# With n_mu = 11 + d_mu, n_tau = 17 + d_tau:
n_mu_eff = 11 + c_ell * abs(1 + PHIp)**alpha_morrey
n_tau_eff = 17 - c_ell * abs(1 + 2*PHIp)**alpha_morrey

print(f"\n  Effective exponents:")
print(f"    n_mu  = {n_mu_eff:.8f}")
print(f"    n_tau = {n_tau_eff:.8f}")

# Define function Q(n_mu, n_tau)
def koide_from_exponents(n1, n2):
    """Q from phi exponents (n_e = 0)."""
    return (1 + PHI**n1 + PHI**n2) / (1 + PHI**(n1/2) + PHI**(n2/2))**2

Q_func = koide_from_exponents(n_mu_eff, n_tau_eff)
print(f"    Q(n_mu, n_tau) = {Q_func:.10f}")

# What pair (n_mu, n_tau) gives EXACT Q = 2/3?
# Constraint: 3*(1 + phi^n1 + phi^n2) = 2*(1 + phi^(n1/2) + phi^(n2/2))^2
print(f"\n  KOIDE CONSTRAINT:")
print(f"    3*(1 + phi^n_mu + phi^n_tau) = 2*(1 + phi^(n_mu/2) + phi^(n_tau/2))^2")

lhs = 3 * (1 + PHI**n_mu_eff + PHI**n_tau_eff)
rhs = 2 * (1 + PHI**(n_mu_eff/2) + PHI**(n_tau_eff/2))**2
print(f"    LHS = {lhs:.6f}")
print(f"    RHS = {rhs:.6f}")
print(f"    LHS - RHS = {lhs - rhs:.6f}")
print(f"    (LHS-RHS)/LHS = {abs(lhs-rhs)/lhs*100:.4f}%")


# =====================================================================
# PART 5: KOIDE DECOMPOSITION - Geometric Meaning
# =====================================================================
print(f"\n{'='*72}")
print("PART 5: Geometric Decomposition of Koide")
print("=" * 72)

# The Koide formula has a beautiful geometric interpretation:
# Define the vector v = (sqrt(m_e), sqrt(m_mu), sqrt(m_tau))
# Then Q = |v|^2 / (v . (1,1,1))^2 * 3
# Wait, more precisely:
# Q = (v1^2 + v2^2 + v3^2) / (v1 + v2 + v3)^2
# This is related to the angle between v and (1,1,1):
# cos^2(theta) = (v.(1,1,1))^2 / (|v|^2 * 3) = (v1+v2+v3)^2 / (3*(v1^2+v2^2+v3^2))
# So: Q = 1/(3*cos^2(theta))
# Q = 2/3 means cos^2(theta) = 1/2, i.e., theta = pi/4 = 45 degrees!

v = np.array([sqrt(m) for m in masses_fw])
u = np.array([1.0, 1.0, 1.0])
cos_theta = np.dot(v, u) / (np.linalg.norm(v) * np.linalg.norm(u))
theta = np.arccos(cos_theta)

print(f"\n  GEOMETRIC INTERPRETATION:")
print(f"  Define v = (sqrt(m_e), sqrt(m_mu), sqrt(m_tau))")
print(f"  and u = (1, 1, 1)")
print(f"  Then Q = |v|^2 / (v.u)^2 = 1 / (3*cos^2(theta))")
print(f"  where theta = angle between v and u")
print(f"\n  Q = 2/3 iff cos^2(theta) = 1/2 iff theta = pi/4 = 45 deg")
print(f"\n  Framework result:")
print(f"    v = ({v[0]:.6f}, {v[1]:.6f}, {v[2]:.6f})")
print(f"    cos(theta) = {cos_theta:.8f}")
print(f"    cos^2(theta) = {cos_theta**2:.8f}")
print(f"    1/2 = 0.50000000")
print(f"    theta = {np.degrees(theta):.6f} deg")
print(f"    45 deg = 45.000000 deg")
print(f"    Diff = {abs(np.degrees(theta) - 45):.4f} deg")

# The question is: WHY should the mass vector make a 45-degree angle
# with the democratic vector (1,1,1)?
print(f"""
  KEY QUESTION: Why theta = 45 degrees?

  In our framework, the mass vector is:
    v = sqrt(m_e) * (1, phi^(n_mu/2), phi^(n_tau/2))
    = sqrt(m_e) * (phi^0, phi^5.5+..., phi^8.5+...)

  The components grow geometrically with ratio phi^(n/2).
  The "democratic direction" (1,1,1) is symmetric.
  The 45-degree angle between a phi-geometric sequence
  and the democratic vector is a NUMBER THEORY statement.
""")


# =====================================================================
# PART 6: CAN THE CORRECTIONS ENFORCE Q = 2/3?
# =====================================================================
print(f"{'='*72}")
print("PART 6: Do the Corrections Enforce Q = 2/3?")
print("=" * 72)

# Scan delta corrections and find Q as function of corrections
print(f"\n  Scanning: what correction alpha gives Q closest to 2/3?")
print(f"\n  {'alpha':>8} {'delta_mu':>10} {'delta_tau':>10} {'Q':>12} {'|Q-2/3|':>12}")

best_Q_alpha = 0
best_Q_diff = 100

for i in range(0, 201):
    alpha_test = i * 0.01
    zp_mu = 1 + PHIp  # = 1/phi^2
    zp_tau = 1 + 2*PHIp  # = 2 - sqrt(5)

    d_mu = c_ell * np.sign(zp_mu) * abs(zp_mu)**alpha_test
    d_tau = c_ell * np.sign(zp_tau) * abs(zp_tau)**alpha_test

    n1 = 11 + d_mu
    n2 = 17 + d_tau
    Q_test = koide_from_exponents(n1, n2)
    diff = abs(Q_test - 2/3)

    if diff < best_Q_diff:
        best_Q_diff = diff
        best_Q_alpha = alpha_test

    if i % 25 == 0 or abs(alpha_test - 0.75) < 0.001:
        print(f"  {alpha_test:8.2f} {d_mu:10.6f} {d_tau:10.6f} {Q_test:12.8f} {diff:12.2e}")

print(f"\n  Best: alpha = {best_Q_alpha:.2f}, |Q - 2/3| = {best_Q_diff:.2e}")
print(f"  Morrey alpha = 0.75: |Q - 2/3| = {abs(Q_fw - 2/3):.2e}")

# Fine scan near best
best_fine_alpha = best_Q_alpha
best_fine_diff = best_Q_diff
for i in range(-100, 101):
    alpha_test = best_Q_alpha + i * 0.0001
    if alpha_test < 0:
        continue
    zp_mu = 1 + PHIp
    zp_tau = 1 + 2*PHIp
    d_mu = c_ell * np.sign(zp_mu) * abs(zp_mu)**alpha_test
    d_tau = c_ell * np.sign(zp_tau) * abs(zp_tau)**alpha_test
    n1 = 11 + d_mu
    n2 = 17 + d_tau
    Q_test = koide_from_exponents(n1, n2)
    diff = abs(Q_test - 2/3)
    if diff < best_fine_diff:
        best_fine_diff = diff
        best_fine_alpha = alpha_test

print(f"  Fine-tuned: alpha = {best_fine_alpha:.4f}, |Q - 2/3| = {best_fine_diff:.2e}")


# =====================================================================
# PART 7: ALTERNATIVE KOIDE REPRESENTATIONS
# =====================================================================
print(f"\n{'='*72}")
print("PART 7: Koide in Alternative Forms")
print("=" * 72)

# Form 1: The "pole" form. Define r_i = sqrt(m_i/M) where M = (m_e+m_mu+m_tau)/3
M_avg = sum(masses_fw) / 3
r = [sqrt(m/M_avg) for m in masses_fw]
print(f"\n  Form 1: r_i = sqrt(m_i / M_avg)")
print(f"    M_avg = {M_avg:.4f} MeV")
for i, name in enumerate(['e', 'mu', 'tau']):
    print(f"    r_{name} = {r[i]:.8f}")
print(f"    sum(r_i^2) = {sum(ri**2 for ri in r):.8f}  (should be 3)")
print(f"    (sum r_i)^2 = {sum(r)**2:.8f}")
print(f"    Q = sum(r^2)/(sum r)^2 * 3 = ... no, Q = 3/sum(r_i)^2 * ... ")

# Form 2: Koide as trace condition
# Define the matrix M_ij = sqrt(m_i * m_j)
# Then Tr(M) = sum(m_i) and (Tr(sqrt(M)))^2 = (sum sqrt(m_i))^2 * ???
# Actually, the simplest form:
# Q = Tr(diag(m)) / (Tr(diag(sqrt(m))))^2

# Form 3: The "angle" parametrization (Koide's original)
# sqrt(m_i) = M0 * (1 + sqrt(2) * cos(theta_0 + 2*pi*i/3))
# with M0 and theta_0 as parameters.
# If Q = 2/3, this parametrization is exact.

# Solve for M0 and theta_0
# sum sqrt(m_i) = 3*M0 (since the cos terms sum to 0)
M0 = sum_sqrt_fw / 3
print(f"\n  Form 2: Koide parametrization")
print(f"    sqrt(m_i) = M0 * (1 + sqrt(2) * cos(theta_0 + 2*pi*i/3))")
print(f"    M0 = sum(sqrt(m_i))/3 = {M0:.8f}")

# For each lepton: sqrt(m_i)/M0 - 1 = sqrt(2) * cos(theta_0 + 2*pi*i/3)
for i, name in enumerate(['e', 'mu', 'tau']):
    ratio = sqrt(masses_fw[i]) / M0 - 1
    cos_val = ratio / sqrt(2)
    print(f"    {name}: sqrt(m)/M0 - 1 = {ratio:+.8f}, cos_val = {cos_val:+.8f}")

# Extract theta_0 from electron
ratio_e = sqrt(masses_fw[0]) / M0 - 1
cos_e = ratio_e / sqrt(2)
if abs(cos_e) <= 1:
    theta_0 = np.arccos(cos_e)
    print(f"\n    theta_0 = {theta_0:.8f} rad = {np.degrees(theta_0):.4f} deg")

    # Verify for mu and tau
    for i, name in enumerate(['mu', 'tau']):
        cos_pred = np.cos(theta_0 + 2*np.pi*(i+1)/3)
        sqrt_m_pred = M0 * (1 + sqrt(2) * cos_pred)
        m_pred = sqrt_m_pred**2
        print(f"    {name}: predicted m = {m_pred:.4f}, actual = {masses_fw[i+1]:.4f}, "
              f"diff = {abs(m_pred - masses_fw[i+1])/masses_fw[i+1]*100:.4f}%")

    # Is theta_0 a framework number?
    print(f"\n    theta_0 = {theta_0:.8f} rad")
    print(f"    theta_0/pi = {theta_0/pi:.8f}")
    print(f"    theta_0 * phi = {theta_0 * PHI:.8f}")
    print(f"    theta_0 * a1 = {theta_0 * a1:.8f}")
    print(f"    2/9 * pi = {2*pi/9:.8f}")
    print(f"    theta_0 / (2*pi/9) = {theta_0 / (2*pi/9):.8f}")
else:
    print(f"\n    cos_e = {cos_e:.8f} > 1, parametrization not applicable")


# =====================================================================
# PART 8: FRAMEWORK NUMBER ANALYSIS OF KOIDE ANGLE
# =====================================================================
print(f"\n{'='*72}")
print("PART 8: Is the Koide Angle a Framework Number?")
print("=" * 72)

# Try: is theta_0 related to 600-cell geometry?
# 600-cell has dihedral angle = arccos(-phi/sqrt(3)) = 164.48 deg
# Or: angles pi/5, 2*pi/5, 3*pi/5 from icosahedral symmetry
# Or: arctan(2) = 63.43 deg (related to icosahedron)

angles_framework = {
    'pi/5': pi/5,
    '2*pi/5': 2*pi/5,
    '3*pi/5': 3*pi/5,
    'pi/a1': pi/a1,
    '2*pi/a1': 2*pi/a1,
    'arctan(phi)': np.arctan(PHI),
    'arctan(2)': np.arctan(2),
    'arctan(sqrt(5))': np.arctan(sqrt(5)),
    'pi/3': pi/3,
    '2*pi/9': 2*pi/9,
    'pi/4': pi/4,
}

if abs(cos_e) <= 1:
    print(f"\n  theta_0 = {theta_0:.8f} rad = {np.degrees(theta_0):.4f} deg")
    print(f"\n  Comparison with framework angles:")
    print(f"  {'Angle':>20} {'Value/rad':>12} {'Value/deg':>10} {'Ratio':>10} {'Match?':>8}")
    for name, val in angles_framework.items():
        ratio = theta_0 / val
        match = "YES!" if abs(ratio - round(ratio)) < 0.02 else ""
        print(f"  {name:>20} {val:12.6f} {np.degrees(val):10.4f} {ratio:10.4f} {match:>8}")


# =====================================================================
# PART 9: EXTENDED KOIDE FOR QUARKS
# =====================================================================
print(f"\n{'='*72}")
print("PART 9: Extended Koide for Quark Sectors")
print("=" * 72)

# The original Koide is for charged leptons. Can we check quarks too?
# Up-type quarks: u, c, t
# Down-type quarks: d, s, b

# Quark masses (pole masses, MeV) - from PDG 2024
m_u = 2.16      # MeV (MS-bar at 2 GeV)
m_c = 1270.0    # MeV
m_t = 172760.0  # MeV (pole mass)

m_d = 4.67      # MeV
m_s = 93.4      # MeV
m_b = 4180.0    # MeV

quarks_up = [m_u, m_c, m_t]
quarks_down = [m_d, m_s, m_b]

def koide_Q(masses):
    s = sum(masses)
    ss = sum(sqrt(m) for m in masses)
    return s / ss**2

Q_up = koide_Q(quarks_up)
Q_down = koide_Q(quarks_down)

print(f"\n  Up-type quarks: Q = {Q_up:.6f}  (2/3 = {2/3:.6f}, diff = {abs(Q_up-2/3):.4f})")
print(f"  Down-type quarks: Q = {Q_down:.6f}  (2/3 = {2/3:.6f}, diff = {abs(Q_down-2/3):.4f})")
print(f"  Charged leptons: Q = {Q_fw:.6f}  (2/3 = {2/3:.6f}, diff = {abs(Q_fw-2/3):.6f})")

print(f"\n  NOTE: Koide works spectacularly for leptons, modestly for quarks.")
print(f"  This is CONSISTENT with our framework:")
print(f"  - Leptons are WHITE nodes (G_k=0, Morrey exponent 3/4)")
print(f"  - Quarks are mixed (different correction mechanisms)")
print(f"  - The Koide formula may be specific to the WHITE/Morrey sector")


# =====================================================================
# PART 10: THE KOIDE-PHI CONNECTION
# =====================================================================
print(f"\n{'='*72}")
print("PART 10: Deep Connection - Koide and the Golden Ratio")
print("=" * 72)

# Check: is there a Koide-like formula that uses phi explicitly?
# Try: sum(m_i) / sum(m_i^x)^(1/x) for different x
print(f"\n  Generalized Koide: Q(x) = sum(m^1) / [sum(m^x)]^(1/x)")
print(f"  x = 1/2 gives standard Koide.")
print(f"\n  {'x':>8} {'Q(x)':>12} {'2/3?':>10}")

for x in [0.25, 1.0/3, 0.4, 0.45, 0.5, 1/PHI, 0.55, 0.6, 2.0/3, 0.75, 1.0]:
    if x < 0.01:
        continue
    sum_mx = sum(m**x for m in masses_fw)
    Qx = sum(masses_fw) / sum_mx**(1.0/x)
    label = ""
    if abs(x - 0.5) < 0.001:
        label = " <-- standard Koide"
    if abs(x - 1/PHI) < 0.01:
        label = " <-- x = 1/phi"
    print(f"  {x:8.4f} {Qx:12.8f} {label}")

# Direct phi test: is Q related to phi?
print(f"\n  Framework number tests on Q:")
print(f"    Q = {Q_fw:.10f}")
print(f"    2/3 = {2/3:.10f}")
print(f"    Q * 3 = {Q_fw*3:.10f} (should be 2.0)")
print(f"    1 - Q = {1-Q_fw:.10f}")
print(f"    1/(3*Q) = {1/(3*Q_fw):.10f} (= cos^2(theta))")
print(f"    phi/a1 * 2 = {PHI/a1*2:.10f}")
print(f"    phi^2/(phi^2+2) = {PHI**2/(PHI**2+2):.10f}")


# =====================================================================
# PART 11: ANALYTIC KOIDE FROM PHI-SEQUENCE
# =====================================================================
print(f"\n{'='*72}")
print("PART 11: Analytic Koide for a Phi-Geometric Sequence")
print("=" * 72)

# For a GEOMETRIC sequence m_i = m_0 * r^i with i=0,1,2:
# Q = (1 + r + r^2) / (1 + sqrt(r) + r)^2
# Note: for large r, Q -> r^2 / r^2 = 1 (not 2/3)
# For r = 1: Q = 3/9 = 1/3
# Q = 2/3 requires solving: 3(1+r+r^2) = 2(1+sqrt(r)+r)^2

# Let s = sqrt(r). Then r = s^2:
# 3(1 + s^2 + s^4) = 2(1 + s + s^2)^2
# 3 + 3*s^2 + 3*s^4 = 2(1 + s^2 + s^4 + 2s + 2s^2 + 2s^3)
# 3 + 3*s^2 + 3*s^4 = 2 + 2s^4 + 4s + 6s^2 + 4s^3
# s^4 - 3*s^2 - 4*s^3 + 1 - 4*s = 0  ... hmm
# s^4 - 4*s^3 - 3*s^2 - 4*s + 1 = 0

# This is a PALINDROMIC polynomial! (coefficients: 1, -4, -3, -4, 1)
# For palindromic: if s is a root, so is 1/s
# Substitution: t = s + 1/s
# s^4 + 1 = (s^2+1)^2 - 2s^2 = (t^2-2)s^2 ... hmm
# s^2(t^2-2) - 4s(s^2+1) - 3s^2 = 0... this is getting messy
# Let me use: divide by s^2: s^2 + 1/s^2 - 4(s + 1/s) - 3 = 0
# With t = s + 1/s: t^2 - 2 - 4t - 3 = 0 => t^2 - 4t - 5 = 0
# t = (4 +/- sqrt(16+20))/2 = (4 +/- 6)/2
# t = 5 or t = -1

print(f"  For geometric sequence m_i = m_0 * r^i (i=0,1,2):")
print(f"  Koide Q = 2/3 requires (with s = sqrt(r)):")
print(f"    s^4 - 4s^3 - 3s^2 - 4s + 1 = 0  [PALINDROMIC!]")
print(f"\n  Substitution t = s + 1/s gives:")
print(f"    t^2 - 4t - 5 = 0")
print(f"    (t - 5)(t + 1) = 0")
print(f"    t = 5 = a1  or  t = -1 = -1")
print(f"\n  *** t = a1 = 5! ***")

print(f"\n  From t = s + 1/s = 5:")
print(f"    s^2 - 5s + 1 = 0")
print(f"    s = (5 +/- sqrt(21))/2")
s1 = (5 + sqrt(21)) / 2
s2 = (5 - sqrt(21)) / 2
print(f"    s = {s1:.8f} or s = {s2:.8f}")
print(f"    r = s^2 = {s1**2:.6f} or r = {s2**2:.6f}")

print(f"\n  From t = s + 1/s = -1:")
print(f"    s^2 + s + 1 = 0")
print(f"    s = exp(+/- 2*pi*i/3) [complex, cube roots of unity]")
print(f"    Not physical (negative masses)")

print(f"\n  RESULT: For a geometric sequence, Q = 2/3 iff")
print(f"  the ratio parameter satisfies s + 1/s = a1 = 5!")
print(f"  This is DIRECTLY connected to the framework integer a1!")

# But our masses are NOT a pure geometric sequence...
# m_e = phi^0, m_mu = phi^~11.08, m_tau = phi^~16.94
# The ratios are phi^11.08 and phi^5.86, not a geometric sequence.
# A geometric would have equal log-ratios.

# However, let's check: what is the "effective" geometric ratio?
log_ratios = [log(masses_fw[1]/masses_fw[0]) / log(PHI),
              log(masses_fw[2]/masses_fw[1]) / log(PHI)]
print(f"\n  Our mass ratios in phi-units:")
print(f"    log_phi(m_mu/m_e)  = {log_ratios[0]:.6f}")
print(f"    log_phi(m_tau/m_mu) = {log_ratios[1]:.6f}")
print(f"    Ratio of ratios    = {log_ratios[0]/log_ratios[1]:.6f}")
print(f"    (For geometric: ratio = 1.0)")

# So it's NOT geometric. The Koide formula is more subtle than
# a simple geometric sequence. The fact that t = a1 = 5 appears
# in the geometric LIMIT is striking but the actual mechanism
# must involve the specific exponents 0, 11, 17 and their corrections.


# =====================================================================
# PART 12: THE KEY IDENTITY: n_mu + n_tau = 28 = N/a1 + N_gen
# =====================================================================
print(f"\n{'='*72}")
print("PART 12: Numerology of the Exponents")
print("=" * 72)

n_mu_bare = 11
n_tau_bare = 17

print(f"\n  Bare exponents: n_mu = {n_mu_bare}, n_tau = {n_tau_bare}")
print(f"  Sum: {n_mu_bare + n_tau_bare}")
print(f"  Diff: {n_tau_bare - n_mu_bare}")
print(f"  Ratio: {n_tau_bare / n_mu_bare:.6f}")
print(f"\n  Framework decompositions:")
print(f"    n_mu + n_tau = {n_mu_bare + n_tau_bare}")
print(f"      = 4 * 7 = d_ST * 7")
print(f"      = N/a1 + N_gen + 1 = 24 + 3 + 1 = 28")
print(f"      = a1^2 + N_gen = 25 + 3 = 28")
print(f"    n_tau - n_mu = {n_tau_bare - n_mu_bare}")
print(f"      = b1 = {b1}")
print(f"      = a1 + 1")
print(f"    n_mu = 11 = A_0 = 2*a1 + 1 (Seeley-DeWitt)")
print(f"    n_tau = 17 = prime, = 12 + a1 = DEG + a1")

# The (a,b) quantum numbers:
# mu: (1,1), so n = 5*1 + 6*1 = 11. z = 1 + phi
# tau: (1,2), so n = 5*1 + 6*2 = 17. z = 1 + 2*phi
print(f"\n  From (a,b):")
print(f"    mu:  (1,1) -> n = 5*1 + 6*1 = 11")
print(f"    tau: (1,2) -> n = 5*1 + 6*2 = 17")
print(f"    Diff in b: 2-1 = 1, contributing b1*(2-1) = 6 to n_tau-n_mu = 6")
print(f"    (Same a=1 for both => same 5a contribution)")


# =====================================================================
# PART 13: KOIDE FROM SUM AND DIFFERENCE RULES
# =====================================================================
print(f"\n{'='*72}")
print("PART 13: Koide Constraint from Framework Rules")
print("=" * 72)

# The Koide formula Q = 2/3 is equivalent to:
# sum(m) = (2/3) * (sum(sqrt(m)))^2
# = (2/3) * (sum(sqrt(m)))^2

# In our framework: m_i = m_e * phi^n_i
# Q = (1 + phi^n1 + phi^n2) / (1 + phi^(n1/2) + phi^(n2/2))^2

# With n1 = 11, n2 = 17:
# The HALF-exponents are 5.5 and 8.5.
# phi^5.5 = phi^5 * phi^(1/2)
# phi^5 = (7 + 3*sqrt(5))/2 (Lucas number L_5 = 11, F_5 = 5 relation)
# phi^(1/2) = sqrt(phi)

# This involves sqrt(phi), which is NOT in Z[phi] but in Q(phi^{1/2})

# Let's try a different approach: can we express Q in terms of
# Fibonacci/Lucas numbers?
# phi^n = F_n * phi + F_{n-1} where F_n is Fibonacci
# phi^11 = F_11 * phi + F_10 = 89*phi + 55 = 199.005...
# phi^17 = F_17 * phi + F_16 = 1597*phi + 987 = 3571.00...

F = [0, 1]
for i in range(2, 20):
    F.append(F[-1] + F[-2])

L = [2, 1]
for i in range(2, 20):
    L.append(L[-1] + L[-2])

phi11 = F[11] * PHI + F[10]
phi17 = F[17] * PHI + F[16]
print(f"\n  phi^11 = F_11*phi + F_10 = {F[11]}*phi + {F[10]} = {phi11:.6f}")
print(f"  phi^17 = F_17*phi + F_16 = {F[17]}*phi + {F[16]} = {phi17:.6f}")
print(f"\n  Numerator of Q: 1 + phi^11 + phi^17")
print(f"  = 1 + ({F[11]}+{F[17]})*phi + ({F[10]}+{F[16]})")
print(f"  = (1+{F[10]}+{F[16]}) + ({F[11]}+{F[17]})*phi")
print(f"  = {1+F[10]+F[16]} + {F[11]+F[17]}*phi")

num_rational = 1 + F[10] + F[16]
num_irrational = F[11] + F[17]
print(f"  = {num_rational} + {num_irrational}*phi")
print(f"  = {num_rational + num_irrational * PHI:.6f}")

# The denominator involves phi^(11/2) and phi^(17/2) which are
# NOT expressible in terms of Fibonacci numbers (half-integer powers)
print(f"\n  Denominator involves phi^(11/2) and phi^(17/2)")
print(f"  These are NOT in Z[phi] -- they require sqrt(phi)")
print(f"  This makes an EXACT algebraic derivation difficult")


# =====================================================================
# FINAL SUMMARY
# =====================================================================
print(f"\n{'='*72}")
print("FINAL SUMMARY")
print("=" * 72)

print(f"""
  KOIDE FORMULA CHECK:
    Q_experimental = {Q_exp:.8f}  (from measured masses)
    Q_framework    = {Q_fw:.8f}  (from 600-cell masses)
    Q_bare         = {Q_bare:.8f}  (without corrections)
    2/3            = {2/3:.8f}

  PRECISION:
    |Q_exp - 2/3| = {abs(Q_exp - 2/3):.2e}  ({abs(Q_exp-2/3)/(2/3)*100:.4f}%)
    |Q_fw  - 2/3| = {abs(Q_fw - 2/3):.2e}  ({abs(Q_fw-2/3)/(2/3)*100:.4f}%)
    |Q_bare- 2/3| = {abs(Q_bare - 2/3):.2e}  ({abs(Q_bare-2/3)/(2/3)*100:.4f}%)

  KEY DISCOVERY:
    For a geometric sequence with Koide Q = 2/3,
    the parameter t = sqrt(r) + 1/sqrt(r) = a1 = 5.
    This connects Koide directly to our framework integer!

  GEOMETRIC INTERPRETATION:
    The mass vector makes angle {np.degrees(theta):.2f} deg with (1,1,1).
    Q = 2/3 requires theta = 45 deg exactly.

  STATUS:
    - Framework masses reproduce Koide to {abs(Q_fw-2/3)/(2/3)*100:.4f}%
    - The Koide-geometric connection gives t = a1 = 5 (DERIVED)
    - Full algebraic proof requires handling phi^(n/2) terms (OPEN)
    - The 3/4 Morrey correction IMPROVES Koide fit (from bare)

  CLASSIFICATION: PATTERN (strong connection to a1=5, not yet fully derived)
""")

print("=" * 72)
print("EXP-411 COMPLETE")
print("=" * 72)

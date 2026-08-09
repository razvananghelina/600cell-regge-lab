"""
exp343c: Electron Mass from Spectral Invariants
================================================
m_e/m_P = alpha^(4*phi^2) = alpha^(1/alpha_s + 2)

Previous attempts (exp327, exp334): 7 routes FAILED for dynamical mechanism.
Key insight from exp334: decomposition alpha^2 * alpha^(1/alpha_s) is ALGEBRAIC.

NEW APPROACHES:
1. Spectral determinant det'(L) of 600-cell Laplacian
2. Analytic torsion / Ray-Singer invariant
3. Heat kernel regularized trace
4. Spectral action Seeley-DeWitt coefficient ratios
5. Z[phi] algebraic identities for z_P = 4*phi^2

CATEGORY: EXPLORATORY (may be NEGATIVE like exp327)
"""

import numpy as np
from math import factorial
import time

t0 = time.time()

# Framework
a1 = 5
b1 = a1 + 1
PHI = (1 + np.sqrt(5)) / 2
phi_prime = (1 - np.sqrt(5)) / 2
N = factorial(a1)  # 120
N_gen = 3
N_eig = 9
rank_E8 = 8
h_E8 = a1 * b1  # 30
degree = 2 * b1  # 12

# Couplings
alpha_s = 1.0 / (2 * PHI**3)
A_c = 2 * np.pi
B_c = -4 * a1 * PHI**4
C_c = 1
disc = B_c**2 - 4 * A_c * C_c
alpha = (-B_c - np.sqrt(disc)) / (2 * A_c)

# Target
z_P = 4 * PHI**2  # = 1/alpha_s + 2 = 2*phi^3 + 2

m_e_GeV = 0.51099895e-3
m_P_GeV = 1.220890e19
ratio_exp = m_e_GeV / m_P_GeV
z_P_exp = np.log(ratio_exp) / np.log(alpha)

print("=" * 72)
print("EXP343c: ELECTRON MASS FROM SPECTRAL INVARIANTS")
print("=" * 72)

print(f"\n  Target: z_P = 4*phi^2 = {z_P:.10f}")
print(f"  m_e/m_P = alpha^z_P = {alpha**z_P:.6e}")
print(f"  Observed: m_e/m_P = {ratio_exp:.6e}")
print(f"  z_P from data: {z_P_exp:.6f}")
print(f"  Error: {abs(z_P - z_P_exp)/z_P_exp*100:.4f}%")

# ============================================================
# SECTION 1: ALGEBRAIC IDENTITIES FOR z_P
# ============================================================
print("\n" + "=" * 72)
print("SECTION 1: ALGEBRAIC IDENTITIES FOR z_P = 4*phi^2")
print("=" * 72)

print(f"""
  z_P = 4*phi^2 = 4*(phi+1) = 4*phi + 4
      = 2*phi^3 + 2 = 1/alpha_s + 2

  DECOMPOSITION:
    z_P = 1/alpha_s + 2 = 2*phi^3 + 2
    m_e/m_P = alpha^2 * alpha^(1/alpha_s) = alpha^2 * alpha^(2*phi^3)

  ALGEBRAIC PROPERTIES:
    z_P in Z[phi]: z_P = 4 + 4*phi (coordinates (4, 4) in Z^2)
    N(z_P) = 16 + 16 - 16 = 16 = (a1-1)^2 = 4^2
    Tr(z_P) = 8 + 4 = 12 = degree = vertex degree
""")

z_P_a, z_P_b = 4, 4
N_zP = z_P_a**2 + z_P_a*z_P_b - z_P_b**2
Tr_zP = 2*z_P_a + z_P_b

print(f"  z_P = ({z_P_a}, {z_P_b}) in Z[phi]")
print(f"  N(z_P) = {N_zP} = (a1-1)^2 = {(a1-1)**2}")
print(f"  Tr(z_P) = {Tr_zP} = degree = {degree}")
print(f"  Min polynomial: t^2 - {Tr_zP}*t + {N_zP} = 0")
print(f"  Roots: z_P = {z_P:.6f}, z_P' = {4 + 4*phi_prime:.6f}")
print(f"  z_P' = 4*phi'^2 = 4*(phi'+1) = 4*phi'+4 = {4+4*phi_prime:.6f}")

# KEY: z_P' = 4*phi'^2 = 4*(1-phi+1) = 4*(2-phi) ... no
# phi' = (1-sqrt(5))/2, phi'^2 = phi'+1 = (3-sqrt(5))/2
# z_P' = 4*(3-sqrt(5))/2 = 2*(3-sqrt(5)) = 6-2*sqrt(5) ~ 1.528
z_P_prime = 4 * phi_prime**2
print(f"  z_P' = {z_P_prime:.6f}")
print(f"  z_P * z_P' = N(z_P) = {z_P * z_P_prime:.6f} = {N_zP}")
print(f"  z_P + z_P' = Tr(z_P) = {z_P + z_P_prime:.6f} = {Tr_zP}")

# Is there a connection to the CC exponent?
z_CC = 57 - alpha_s
print(f"\n  CONNECTION TO CC:")
print(f"  z_CC = {z_CC:.6f}")
print(f"  z_P = {z_P:.6f}")
print(f"  z_CC / z_P = {z_CC / z_P:.6f}")
print(f"  z_CC - z_P = {z_CC - z_P:.6f}")
print(f"  z_CC * z_P = {z_CC * z_P:.6f}")
print(f"  (z_CC - z_P) / alpha_s = {(z_CC - z_P) / alpha_s:.4f}")

# z_CC / z_P ~ 5.43. Is this a framework number?
ratio_CC_P = z_CC / z_P
print(f"\n  z_CC/z_P = {ratio_CC_P:.6f}")
print(f"    a1 + alpha_s = {a1 + alpha_s:.6f}")
print(f"    Diff: {ratio_CC_P - (a1 + alpha_s):.6f} ({abs(ratio_CC_P-(a1+alpha_s))/(a1+alpha_s)*100:.3f}%)")
# Not a great match

# What about integer combinations?
print(f"\n  z_CC = c * z_P + d for (c,d):")
for c in range(1, 8):
    d = z_CC - c * z_P
    print(f"    c={c}: d = {d:.4f}, d in framework? ", end="")
    # Check if d is close to a framework number
    for name, val in [("a1", a1), ("b1", b1), ("N_gen", N_gen), ("degree", degree),
                       ("alpha_s", alpha_s), ("phi", PHI), ("phi^2", PHI**2),
                       ("0", 0), ("2", 2), ("4", 4)]:
        if abs(d - val) < 0.5:
            print(f" ~ {name} ({val}), diff={d-val:.4f}", end="")
    print()

# ============================================================
# SECTION 2: SPECTRAL DETERMINANT
# ============================================================
print("\n" + "=" * 72)
print("SECTION 2: SPECTRAL DETERMINANT OF 600-CELL LAPLACIAN")
print("=" * 72)

# Compute the 600-cell Laplacian eigenvalues
# Using character theory of 2I (from exp343a):
# Eigenvalue of adj matrix for irrep R_k: lambda_k = 12 * chi_k(12A) / d_k
# Laplacian eigenvalue: mu_k = 12 - lambda_k
# Multiplicity: d_k^2

sin_pi_10 = np.sin(np.pi/10)

def su2_char(d, theta):
    if abs(np.sin(theta/2)) < 1e-15:
        return d
    return np.sin(d * theta / 2) / np.sin(theta / 2)

theta_12A = np.pi / 5
dims_2I = [1, 2, 3, 4, 5, 6]

# Standard irreps (from SU(2))
chi_std = {}
for d in dims_2I:
    chi_std[d] = su2_char(d, theta_12A)

# Conjugate irreps
chi_R10 = su2_char(10, theta_12A)  # sin(10*pi/10)/sin(pi/10) = sin(pi)/... = 0
chi_R8 = su2_char(8, theta_12A)
chi_R7 = su2_char(7, theta_12A)

# R4' (dim 4): from j=9/2 -> R4' + R6
chi_R4p = chi_R10 - chi_std[6]
# R3' (dim 3): from j=3 -> R3' + R4'
chi_R3p = chi_R7 - chi_R4p
# R2' (dim 2): from j=7/2 -> R2' + R6
chi_R2p = chi_R8 - chi_std[6]

# All 9 irreps: (label, dim, chi_at_12A)
irreps = [
    ('R1', 1, chi_std[1]),
    ('R2', 2, chi_std[2]),
    ("R2'", 2, chi_R2p),
    ('R3', 3, chi_std[3]),
    ("R3'", 3, chi_R3p),
    ('R4', 4, chi_std[4]),
    ("R4'", 4, chi_R4p),
    ('R5', 5, chi_std[5]),
    ('R6', 6, chi_std[6]),
]

print(f"\n  600-cell Laplacian eigenvalues:")
print(f"  {'Irrep':>6} {'dim':>4} {'chi(12A)':>10} {'adj_eig':>10} {'lap_eig':>10} {'mult':>6}")

lap_eigs = []  # (eigenvalue, multiplicity)
for label, d_k, chi_val in irreps:
    adj_eig = 12 * chi_val / d_k
    lap_eig = 12 - adj_eig
    mult = d_k**2
    print(f"  {label:>6} {d_k:>4} {chi_val:>10.6f} {adj_eig:>10.6f} {lap_eig:>10.6f} {mult:>6}")
    lap_eigs.append((lap_eig, mult))

# Nonzero eigenvalues (exclude R1 which has lap_eig = 0)
nonzero = [(l, m) for l, m in lap_eigs if abs(l) > 1e-10]

# Spectral determinant: det'(L) = prod lambda_i^{m_i}
log_det = sum(m * np.log(abs(l)) for l, m in nonzero)
print(f"\n  log det'(L) = {log_det:.6f}")
print(f"  det'(L) = exp({log_det:.2f})")

# Can this be related to z_P?
print(f"\n  Ratios with framework numbers:")
print(f"  log det'(L) / z_P = {log_det / z_P:.6f}")
print(f"  log det'(L) / z_CC = {log_det / z_CC:.6f}")
print(f"  log det'(L) / N = {log_det / N:.6f}")
print(f"  log det'(L) / (N-1) = {log_det / (N-1):.6f}")
print(f"  log det'(L) / degree = {log_det / degree:.6f}")
print(f"  log det'(L) / ln(1/alpha) = {log_det / np.log(1/alpha):.6f}")

# The relation m_e/m_P = alpha^z_P means:
# ln(m_e/m_P) = z_P * ln(alpha) = -z_P * ln(1/alpha)
# Can log_det = -z_P * c for some c?
print(f"\n  log det'(L) / (-z_P) = {log_det / (-z_P):.6f}")
print(f"  ln(alpha) = {np.log(alpha):.6f}")
print(f"  Ratio: {(log_det / (-z_P)) / np.log(alpha):.6f}")

# ============================================================
# SECTION 3: HEAT KERNEL TRACE
# ============================================================
print("\n" + "=" * 72)
print("SECTION 3: HEAT KERNEL AND FUNCTIONAL DETERMINANT")
print("=" * 72)

# Spectral zeta function: zeta(s) = sum m_i * lambda_i^{-s}
def zeta_L(s):
    return sum(m * l**(-s) for l, m in nonzero if l > 0)

# Heat trace: K(t) = sum m_i * exp(-t*lambda_i)
def heat_trace(t):
    return sum(m * np.exp(-t*l) for l, m in lap_eigs)

# The heat trace has an asymptotic expansion:
# K(t) ~ sum_{k=0}^{infty} c_k * t^{k-n/2} for t -> 0
# where n is the "dimension" and c_k are Seeley-DeWitt coefficients.

# The functional determinant is:
# -zeta'(0) = log det'(L)
# zeta(0) = sum m_i (for nonzero eigenvalues) = N - 1 = 119

zeta_0 = sum(m for l, m in nonzero)
print(f"\n  zeta(0) = sum of multiplicities = {zeta_0} = N - 1 = {N-1}")

# Compute zeta'(0) numerically
ds = 1e-7
zeta_prime_0 = (zeta_L(ds) - zeta_L(-ds)) / (2*ds)
# More stable: use the series representation
# zeta'(s) = -sum m_i * lambda_i^{-s} * log(lambda_i)
zeta_prime_0_v2 = -sum(m * np.log(l) for l, m in nonzero if l > 0)

print(f"  zeta'(0) (finite diff) ~ {zeta_prime_0:.6f}")
print(f"  zeta'(0) (analytic) = {zeta_prime_0_v2:.6f}")
print(f"  -zeta'(0) = log det'(L) = {-zeta_prime_0_v2:.6f}")
print(f"  Matches direct: {abs(-zeta_prime_0_v2 - log_det) < 1e-6}")

# Special values of zeta
print(f"\n  Special zeta values:")
for s in [0.5, 1, 1.5, 2, 3, -1, -2]:
    try:
        z = zeta_L(s)
        print(f"    zeta({s:>5.1f}) = {z:>14.6f}")
    except Exception:
        print(f"    zeta({s:>5.1f}) = ERROR")

# Heat trace at special times
print(f"\n  Heat trace K(t):")
for t in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    K = heat_trace(t)
    print(f"    K({t:>5.2f}) = {K:>14.6f}")

# At t=0: K(0) = sum m_i = 120 (all eigenvalues including zero)
print(f"  K(0) = {sum(m for _, m in lap_eigs)} = N = {N}")
# At t->inf: K(inf) -> m_0 = 1 (only zero eigenvalue contributes)

# ============================================================
# SECTION 4: SEELEY-DEWITT COEFFICIENTS
# ============================================================
print("\n" + "=" * 72)
print("SECTION 4: SEELEY-DEWITT COEFFICIENTS")
print("=" * 72)

# From the spectral action on the 600-cell:
# c_0 = 2640, c_1 = 14880, c_2 = 55920 (from exp342)
c_0 = 2640
c_1 = 14880
c_2 = 55920

print(f"\n  Known Seeley-DeWitt coefficients:")
print(f"  c_0 = {c_0} = {c_0//N} * N = {c_0//N} * {N}")
print(f"  c_1 = {c_1}")
print(f"  c_2 = {c_2}")

# Ratios
print(f"\n  Ratios:")
print(f"  c_1/c_0 = {c_1/c_0:.6f} = {c_1//c_0}... wait")
# Check if rational
from fractions import Fraction
r10 = Fraction(c_1, c_0)
r21 = Fraction(c_2, c_1)
r20 = Fraction(c_2, c_0)
print(f"  c_1/c_0 = {r10} = {float(r10):.6f}")
print(f"  c_2/c_1 = {r21} = {float(r21):.6f}")
print(f"  c_2/c_0 = {r20} = {float(r20):.6f}")

# c_1/c_0 related to z_P?
print(f"\n  c_1/c_0 = {float(r10):.6f}")
print(f"  z_P/2 = {z_P/2:.6f} ... NO")
print(f"  c_1/(c_0 * alpha) = {c_1/(c_0*alpha):.4f}")
print(f"  But c_1/c_0 = {r10} in exact form")

# The VEV Lambda^2 = c_1/(2*c_0) from spectral action
Lambda2 = float(Fraction(c_1, 2*c_0))
print(f"\n  Higgs VEV^2: Lambda^2 = c_1/(2*c_0) = {Fraction(c_1, 2*c_0)} = {Lambda2:.6f}")

# m_H^2 = 2*c_0/c_2 * Lambda^2 in the spectral action
# Actually: m_H^2 = 8*lambda*v^2 where lambda = pi^2*a/(f0*b)
# This gets complicated. Let's focus on ratios.

# Can z_P emerge from the spectral action coefficients?
print(f"\n  Searching for z_P = {z_P:.6f} in spectral ratios:")
candidates_zP = []
for expr, val in [
    ("c_1/c_0", c_1/c_0),
    ("c_2/c_1", c_2/c_1),
    ("c_2/c_0", c_2/c_0),
    ("sqrt(c_1/c_0)", np.sqrt(c_1/c_0)),
    ("c_0/N", c_0/N),
    ("c_1/N", c_1/N),
    ("c_2/N", c_2/N),
    ("c_0/(N*degree)", c_0/(N*degree)),
    ("ln(c_0)", np.log(c_0)),
    ("ln(c_1)", np.log(c_1)),
    ("ln(c_2/c_0)", np.log(c_2/c_0)),
    ("c_0/c_1*degree^2", c_0/c_1*degree**2),
]:
    diff = abs(val - z_P)
    pct = diff / z_P * 100
    candidates_zP.append((expr, val, pct))

candidates_zP.sort(key=lambda x: x[2])
for expr, val, pct in candidates_zP:
    mark = " ***" if pct < 5 else ""
    print(f"    {expr:<25} = {val:>12.6f} ({pct:>8.3f}% off){mark}")

# ============================================================
# SECTION 5: PRODUCT FORMULAS AND z_P
# ============================================================
print("\n" + "=" * 72)
print("SECTION 5: PRODUCT FORMULAS FOR z_P")
print("=" * 72)

# z_P = 4*phi^2 = 4*(phi+1)
# z_P = 2*(2*phi+2) = 2*(phi^3+1)
# z_P = 2*phi^3 + 2

# Can z_P emerge from some product of eigenvalues?
# Product of distinct Cayley eigenvalues?
distinct_eigs = sorted(set(l for l, m in nonzero if l > 0))
print(f"\n  Distinct positive Laplacian eigenvalues:")
for l in distinct_eigs:
    print(f"    {l:.6f}")

prod_distinct = np.prod(distinct_eigs)
sum_distinct = sum(distinct_eigs)
print(f"\n  Product of distinct positive eigs: {prod_distinct:.6f}")
print(f"  Sum of distinct positive eigs: {sum_distinct:.6f}")
print(f"  z_P = {z_P:.6f}")
print(f"  sum/z_P = {sum_distinct/z_P:.6f}")
print(f"  prod/z_P = {prod_distinct/z_P:.6f}")

# Average eigenvalue
avg_eig = sum(l*m for l, m in nonzero) / sum(m for l, m in nonzero)
print(f"\n  Average nonzero eigenvalue: {avg_eig:.6f}")
print(f"  degree (=12): {degree}")
print(f"  avg/z_P = {avg_eig/z_P:.6f}")

# Trace of Laplacian (= sum of all eigenvalues = N * degree)
trace_L = sum(l*m for l, m in lap_eigs)
print(f"\n  Tr(L) = {trace_L:.1f} = N*degree = {N*degree}")

# ============================================================
# SECTION 6: THE FUNDAMENTAL ALGEBRAIC IDENTITY
# ============================================================
print("\n" + "=" * 72)
print("SECTION 6: ALGEBRAIC IDENTITY ANALYSIS")
print("=" * 72)

print(f"""
  The relation m_e/m_P = alpha^z_P is equivalent to:
    ln(m_e/m_P) = z_P * ln(alpha)
    ln(m_e/m_P) = (1/alpha_s + 2) * ln(alpha)

  In the framework, alpha is determined by:
    2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0

  And alpha_s = 1/(2*phi^3).

  So z_P = 2*phi^3 + 2 and we need:
    ln(m_e/m_P) = (2*phi^3 + 2) * ln(alpha)

  This is a single equation relating three quantities:
    m_e, m_P, alpha (where m_P = sqrt(hbar*c/G))

  Since alpha is DERIVED from a1=5 (via the quadratic),
  and z_P = 4*phi^2 is also from a1=5,
  the ONLY free quantity is one of {{m_e, m_P, G}}.

  Traditionally, G (Newton's constant) is considered fundamental.
  Then: m_e = m_P * alpha^z_P = sqrt(hbar*c/G) * alpha^(4*phi^2)
  This PREDICTS m_e given G (or equivalently, G given m_e).
""")

# Compute the prediction both ways
from math import sqrt
HBAR = 1.054571817e-34
C_light = 299792458
G_N = 6.67430e-11
E_CHARGE = 1.602176634e-19
GeV_to_kg = E_CHARGE * 1e9 / C_light**2

m_P_kg = sqrt(HBAR * C_light / G_N)
m_P_pred = m_P_kg / GeV_to_kg  # in GeV

m_e_pred_GeV = m_P_pred * alpha**z_P
m_e_exp_GeV = 0.51099895e-3

print(f"  Prediction: m_e = m_P * alpha^z_P")
print(f"    m_P = {m_P_pred:.6e} GeV")
print(f"    alpha = {alpha:.10e}")
print(f"    z_P = {z_P:.10f}")
print(f"    m_e(pred) = {m_e_pred_GeV:.6e} GeV")
print(f"    m_e(exp)  = {m_e_exp_GeV:.6e} GeV")
print(f"    Error: {abs(m_e_pred_GeV - m_e_exp_GeV)/m_e_exp_GeV*100:.4f}%")

# Or predict G from m_e:
# m_e/m_P = alpha^z_P => m_P = m_e/alpha^z_P
# m_P = sqrt(hbar*c/G) => G = hbar*c/m_P^2
m_P_from_me = m_e_exp_GeV / alpha**z_P
G_pred = HBAR * C_light / (m_P_from_me * GeV_to_kg)**2

print(f"\n  Inverse prediction: G from m_e")
print(f"    m_P = m_e/alpha^z_P = {m_P_from_me:.6e} GeV")
print(f"    G(pred) = {G_pred:.6e} m^3/(kg*s^2)")
print(f"    G(exp)  = {G_N:.6e}")
print(f"    Error: {abs(G_pred - G_N)/G_N*100:.4f}%")

# ============================================================
# SECTION 7: z_P FROM REPRESENTATION THEORY
# ============================================================
print("\n" + "=" * 72)
print("SECTION 7: z_P FROM 2I REPRESENTATION THEORY")
print("=" * 72)

# Can z_P = 4*phi^2 emerge from some representation-theoretic quantity?
dims_all = [1, 2, 2, 3, 3, 4, 4, 5, 6]

print(f"\n  2I irrep dimensions: {dims_all}")
print(f"  N_eig = {N_eig} irreps")

# Various sums
S1 = sum(dims_all)
S2 = sum(d**2 for d in dims_all)
S3 = sum(d**3 for d in dims_all)
S4 = sum(d**4 for d in dims_all)

print(f"  S1 = sum(d) = {S1} = h_E8 = {h_E8}")
print(f"  S2 = sum(d^2) = {S2} = N = {N}")
print(f"  S3 = sum(d^3) = {S3}")
print(f"  S4 = sum(d^4) = {S4}")

# Ratios
print(f"\n  Ratios:")
print(f"  S3/S2 = {S3/S2:.6f}")
print(f"  S3/S1 = {S3/S1:.6f}")
print(f"  S2/S1 = {S2/S1:.6f} = N/h = 120/30 = 4")
print(f"  S3/(S1*S2) = {S3/(S1*S2):.6f}")
print(f"  S4/S2 = {S4/S2:.6f}")
print(f"  z_P = {z_P:.6f}")

# Check: S2/S1 = 4 = a1-1. And z_P = 4*phi^2 = (S2/S1)*phi^2.
print(f"\n  KEY: S2/S1 = {S2//S1} = N/h = a1-1 = {a1-1}")
print(f"  And z_P = (S2/S1)*phi^2 = {a1-1}*phi^2 = {(a1-1)*PHI**2:.6f}")
print(f"  Matches z_P = {z_P:.6f}? {'YES' if abs((a1-1)*PHI**2 - z_P) < 1e-10 else 'NO'}")

print(f"""
  FORMULA: z_P = (N/h) * phi^2 = (|2I|/h(E8)) * phi^2

  Where:
    N = |2I| = 120 (order of binary icosahedral group)
    h = h(E8) = 30 (Coxeter number)
    phi^2 = phi + 1 (golden ratio squared)
    N/h = 4 = a1-1 (a framework number)

  So: z_P = (a1-1) * (phi+1) = (a1-1)*phi^2

  Equivalently:
    z_P = a1*phi^2 - phi^2 = a1*(phi+1) - (phi+1)
        = a1*phi + a1 - phi - 1
        = (a1-1)*phi + (a1-1) = (a1-1)*(phi+1)

  Or: z_P = 4*(phi+1) = 4*phi^2 (as before, with 4 = a1-1)

  But WHY (a1-1)*phi^2?
  N/h(E8) counts the number of "Coxeter orbits" on 2I.
  Each orbit has h=30 elements, and there are 120/30 = 4 orbits.
  phi^2 = phi+1 is the "Perron value" per orbit (from the golden ratio).

  INTERPRETATION: m_e/m_P = alpha^(orbits * phi^2)
  Each Coxeter orbit contributes phi^2 to the Planck exponent.
""")

# ============================================================
# SECTION 8: COMPARISON WITH CC
# ============================================================
print("\n" + "=" * 72)
print("SECTION 8: COMPARING z_P AND z_CC STRUCTURES")
print("=" * 72)

print(f"""
  z_P  = (a1-1) * phi^2            = 4*phi^2    = {z_P:.6f}
  z_CC = N_gen * (a1^2-a1-1) - alpha_s = 57 - alpha_s = {z_CC:.6f}

  Both come from a1 = 5. Can we unify them?

  z_P  = (a1-1) * phi^2 = (a1-1) * (phi+1)
  z_CC = N_gen * N(a1+(a1+1)*phi) - 1/(2*phi^3)
       = 3 * 19 - alpha_s

  Ratio: z_CC/z_P = {z_CC/z_P:.6f}
  Product: z_CC * z_P = {z_CC*z_P:.6f}

  Is z_CC/z_P special?
    z_CC/z_P = (57-alpha_s) / (4*phi^2)
    = (57 - phi + 3/2) / (4*phi+4)  [using alpha_s = phi - 3/2]
    = (58.5 - phi) / (4*phi + 4)

  In exact form:
    (117/2 - phi) / (4*phi + 4) = (117 - 2*phi) / (8*phi + 8)

  Product: z_CC * z_P = (58.5 - phi) * (4*phi+4) = (117-2*phi)(4*phi+4)/2
    = (117*4*phi + 117*4 - 8*phi^2 - 8*phi)/2
    = (468*phi + 468 - 8*phi - 8 - 8*phi)/2
    = (452*phi + 460)/2 = 226*phi + 230
""")

prod_zCzP = z_CC * z_P
print(f"  z_CC * z_P = {prod_zCzP:.6f}")
print(f"  226*phi + 230 = {226*PHI + 230:.6f}")
print(f"  Match: {abs(prod_zCzP - (226*PHI+230)) < 0.1}")

# N(z_CC * z_P in Z[phi]):
# z_CC*z_P = 230 + 226*phi
# N = 230^2 + 230*226 - 226^2 = 52900 + 51980 - 51076 = 53804
N_prod = 230**2 + 230*226 - 226**2
print(f"  N(z_CC*z_P) = {N_prod}")
print(f"  = {N_prod} = 4 * {N_prod//4} = 4 * {N_prod//4}")
# Factor
n = N_prod
factors = []
for p in range(2, 1000):
    while n % p == 0:
        factors.append(p)
        n //= p
if n > 1:
    factors.append(n)
print(f"  Factorization: {' * '.join(str(f) for f in factors)}")

# ============================================================
# SECTION 9: HONEST ASSESSMENT
# ============================================================
print("\n" + "=" * 72)
print("SECTION 9: HONEST ASSESSMENT")
print("=" * 72)

print(f"""
  RESULTS:

  1. z_P = (N/h) * phi^2 = (a1-1) * phi^2
     This gives a representation-theoretic interpretation:
     z_P counts "Coxeter orbits" (N/h = 4) weighted by phi^2.
     CATEGORY: STRUCTURAL OBSERVATION (not yet a mechanism)

  2. z_P has Z[phi] norm N(z_P) = (a1-1)^2 = 16 and trace Tr(z_P) = 12 = degree.
     These connections to framework numbers are NICE but not a derivation.

  3. The spectral determinant log det'(L) = {log_det:.2f}
     This is NOT simply related to z_P = {z_P:.6f}.
     NO spectral determinant mechanism for m_e.

  4. The Seeley-DeWitt coefficients do NOT directly give z_P.
     NO spectral action mechanism (confirming exp327 negative result).

  5. m_e prediction from m_e = m_P * alpha^z_P works to {abs(m_e_pred_GeV-m_e_exp_GeV)/m_e_exp_GeV*100:.2f}%.
     Equivalently, G prediction from G = hbar*c*(alpha^z_P / m_e)^2 works to {abs(G_pred-G_N)/G_N*100:.2f}%.

  CONCLUSION:
  The relation m_e/m_P = alpha^(4*phi^2) remains an ALGEBRAIC PATTERN.
  No spectral, determinantal, or heat-kernel mechanism produces z_P = 4*phi^2.
  The best interpretation is:
    z_P = (Coxeter orbits) * phi^2 = (N/h) * phi^2
  But WHY the Planck hierarchy is controlled by this quantity is OPEN.

  NEW INSIGHT: z_P = (a1-1)*phi^2 where:
    - a1-1 = N/h(E8) = number of Coxeter orbits = 4
    - phi^2 = phi+1 = Perron eigenvalue of the golden ratio

  CATEGORY: ALGEBRAIC PATTERN (unchanged from exp334)
  STATUS: OPEN (as predicted - may be fundamentally unfixable)
""")

elapsed = time.time() - t0
print(f"\nCompleted in {elapsed:.1f}s")

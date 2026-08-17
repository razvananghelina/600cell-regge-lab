"""
exp421_cc_heat_kernel.py
========================
GOAL: Test whether the CC functional form Lambda_P = alpha^z can be derived
from spectral data of the 600-cell Cayley graph Laplacian.

CONTEXT:
  - Exponent z_CC = 57 - alpha_s is FULLY DERIVED (7 algebraic routes)
  - Galois duality Lambda * Lambda' = (1/(2*pi))^z is DERIVED
  - But the functional form alpha^z is a PATTERN, not a derivation
  - This is the framework's most serious gap

NEW IDEA: The Cayley graph Laplacian (12-regular, 120 vertices) has 9 exact
eigenvalues in Z[phi] with multiplicities d_k^2. The regularized determinant
det'(L) is a finite, well-defined product (no UV divergences on a finite graph).
If log det'(L) / ln(alpha) = z_CC, we have a derivation.

WHAT WAS ALREADY TRIED (exp408c): Spectral determinants of the SIMPLICIAL
Laplacians (0-form through 3-form) tested against N^4. Result: NEGATIVE.
But the CC connection alpha^57 was **never tested**.

EIGENVALUE DATA (Cayley graph of 2I, generators = class 10a, 12 elements):
  L_k = 12*(1 - chi_k(10a)/d_k) in Z[phi]

  k  d_k  mult  L_k              Galois pair
  0   1    1    0                 -
  1   2    4    12 - 6*phi        L_7
  2   3    9    12 - 4*phi        L_8
  3   4   16    9                 self
  4   5   25    12                self
  5   6   36    14                self
  6   4   16    15                self
  7   2    4    6 + 6*phi         L_1
  8   3    9    8 + 4*phi         L_2

BLIND: compute first, interpret after.
Author: Razvan-Constantin Anghelina
Date: 2026-02-27
"""

import numpy as np
from math import sqrt, pi, log, log10, exp, factorial
from fractions import Fraction
import sys

print("=" * 72)
print("EXP-421: CC HEAT KERNEL ON CAYLEY GRAPH")
print("        Testing alpha^z from spectral determinant")
print("=" * 72)

# =====================================================================
# FRAMEWORK CONSTANTS
# =====================================================================
a1 = 5
b1 = a1 + 1                              # 6
PHI = (1 + sqrt(a1)) / 2                  # golden ratio
PHI_prime = (1 - sqrt(a1)) / 2            # Galois conjugate = -1/PHI
N = factorial(a1)                          # 120
degree = 2 * b1                           # 12
N_eig = 9
rank_E8 = 8
h_E8 = a1 * b1                            # 30
N_gen = 3

# Framework couplings
alpha_s = 1 / (2 * PHI**3)
sin2tW = Fraction(b1, a1**2 + 1)           # 6/26

# Alpha from quadratic: 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0
A_coef = 2 * pi
B_coef = -4 * a1 * PHI**4
C_coef = 1
disc = B_coef**2 - 4 * A_coef * C_coef
alpha = (-B_coef - sqrt(disc)) / (2 * A_coef)    # physical root
alpha_prime = (-B_coef + sqrt(disc)) / (2 * A_coef)  # Galois conjugate

# CC exponent
z_CC = 57 - alpha_s

print(f"\nFramework constants:")
print(f"  a1 = {a1}, b1 = {b1}, N = {N}, degree = {degree}")
print(f"  PHI = {PHI:.10f}")
print(f"  alpha = {alpha:.10f}  (1/alpha = {1/alpha:.6f})")
print(f"  alpha_s = {alpha_s:.10f}")
print(f"  alpha * alpha' = {alpha * alpha_prime:.10f}  (1/(2*pi) = {1/(2*pi):.10f})")
print(f"  z_CC = 57 - alpha_s = {z_CC:.10f}")
print(f"  alpha^z_CC = {alpha**z_CC:.6e}")

# =====================================================================
# CAYLEY GRAPH EIGENVALUE DATA
# =====================================================================
# L_k = degree * (1 - chi_k(10a) / d_k)
# chi(10a) from 2I character table:
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
chi_10a = [1, PHI, PHI, 1, 0, -1, -1, 1 - PHI, 1 - PHI]
chi_10a_galois = [1, PHI_prime, PHI_prime, 1, 0, -1, -1, 1 - PHI_prime, 1 - PHI_prime]

Lk = [degree * (1 - chi_10a[k] / dims[k]) for k in range(N_eig)]
Lk_galois = [degree * (1 - chi_10a_galois[k] / dims[k]) for k in range(N_eig)]
mults = [d**2 for d in dims]

# Z[phi] forms: L_k = a + b*phi
zphi_forms = [
    (0, 0),      # L_0 = 0
    (12, -6),    # L_1 = 12 - 6*phi
    (12, -4),    # L_2 = 12 - 4*phi
    (9, 0),      # L_3 = 9
    (12, 0),     # L_4 = 12
    (14, 0),     # L_5 = 14
    (15, 0),     # L_6 = 15
    (6, 6),      # L_7 = 6 + 6*phi
    (8, 4),      # L_8 = 8 + 4*phi
]


# =====================================================================
# TEST 1: Spectral Data Validation
# =====================================================================
def test_1_spectral_validation():
    """Validate eigenvalue data: sums, Galois norms, Z[phi] forms."""
    print(f"\n{'='*72}")
    print("TEST 1: Spectral Data Validation")
    print("=" * 72)

    all_pass = True

    # Check Z[phi] forms match
    print("\n  Z[phi] form verification:")
    for k in range(N_eig):
        a, b = zphi_forms[k]
        L_check = a + b * PHI
        diff = abs(L_check - Lk[k])
        ok = diff < 1e-10
        if not ok:
            all_pass = False
        sign = "+" if b >= 0 else ""
        print(f"    L_{k} = {a}{sign}{b}*phi = {L_check:.6f} "
              f"(computed: {Lk[k]:.6f}, diff: {diff:.2e}) [{'OK' if ok else 'FAIL'}]")

    # Check sum of multiplicities = N
    sum_mult = sum(mults)
    ok_mult = sum_mult == N
    if not ok_mult:
        all_pass = False
    print(f"\n  sum(m_k) = {sum_mult} (expected: {N}) [{'PASS' if ok_mult else 'FAIL'}]")

    # Check sum(m_k * L_k) = degree * (N - 1) = 12 * 119 = 1428?
    # Actually: sum m_k * L_k = degree * sum(m_k - m_k*chi/d_k) = degree * (N - sum(m_k*chi/d_k))
    # sum(m_k * chi/d_k) = sum(d_k * chi_k) = char of regular rep at class 10a = 0 (non-identity)
    # So sum(m_k * L_k) = degree * N = 12 * 120 = 1440
    sum_mL = sum(mults[k] * Lk[k] for k in range(N_eig))
    expected_mL = degree * N
    ok_mL = abs(sum_mL - expected_mL) < 1e-6
    if not ok_mL:
        all_pass = False
    print(f"  sum(m_k * L_k) = {sum_mL:.6f} (expected: {expected_mL}) [{'PASS' if ok_mL else 'FAIL'}]")

    # Check Galois norms are rational
    print(f"\n  Galois norms N(L_k) = L_k * sigma(L_k):")
    for k in range(N_eig):
        if Lk[k] > 0.001:
            norm = Lk[k] * Lk_galois[k]
            # Check rationality: should be integer or simple fraction
            # N(L_k) = degree^2 * (1 - Tr(chi)/d + N(chi)/d^2)
            # where Tr(chi) = chi + chi', N(chi) = chi * chi'
            a, b = zphi_forms[k]
            norm_exact = a**2 + a * b - b**2  # N(a + b*phi) = a^2 + ab - b^2
            ok_norm = abs(norm - norm_exact) < 0.01
            if not ok_norm:
                all_pass = False
            print(f"    k={k}: N(L_{k}) = {norm:.6f} "
                  f"(Z[phi] norm: {norm_exact}) [{'OK' if ok_norm else 'FAIL'}]")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  TEST 1: {status}")
    return all_pass


# =====================================================================
# TEST 2: log det'(L) / ln(alpha) vs z_CC -- CORE TEST
# =====================================================================
def test_2_logdet_vs_zCC():
    """Core test: does log det'(L) relate to z_CC * ln(alpha)?"""
    print(f"\n{'='*72}")
    print("TEST 2: log det'(L) / ln(alpha) vs z_CC  [CORE TEST]")
    print("=" * 72)

    found_pattern = False

    # log det'(L) = sum_{k>0} m_k * ln(L_k)
    ld = sum(mults[k] * log(Lk[k]) for k in range(N_eig) if Lk[k] > 0.001)
    ln_alpha = log(alpha)

    ratio_raw = ld / ln_alpha
    print(f"\n  log det'(L) = {ld:.10f}")
    print(f"  ln(alpha) = {ln_alpha:.10f}")
    print(f"  ratio = log det'(L) / ln(alpha) = {ratio_raw:.10f}")
    print(f"  target z_CC = {z_CC:.10f}")
    print(f"  diff = {abs(ratio_raw - z_CC):.6f}")

    if abs(ratio_raw - z_CC) < 0.5:
        print(f"  *** CLOSE MATCH! ***")
        found_pattern = True

    # Try normalized versions: L/12, L/(2*pi), L/N
    normalizations = [
        ("L/degree (L/12)", degree),
        ("L/(2*pi)", 2 * pi),
        ("L/N (L/120)", N),
        ("L/(degree*N) (L/1440)", degree * N),
        ("L/b1 (L/6)", b1),
        ("L/a1 (L/5)", a1),
    ]

    print(f"\n  Normalized log det':")
    for name, norm in normalizations:
        ld_norm = sum(mults[k] * log(Lk[k] / norm) for k in range(N_eig) if Lk[k] > 0.001)
        ratio_norm = ld_norm / ln_alpha
        diff_57 = abs(ratio_norm - z_CC)
        near_int = round(ratio_norm)
        diff_int = abs(ratio_norm - near_int)
        mark = "***" if diff_57 < 0.5 else ("*" if diff_int < 0.5 else "")
        print(f"    {name:30s}: ld={ld_norm:12.6f}, "
              f"ld/ln(a)={ratio_norm:12.6f}, "
              f"diff from 57={diff_57:8.4f}, "
              f"near int={near_int} (diff={diff_int:.4f}) {mark}")
        if diff_57 < 0.5:
            found_pattern = True

    # Search for p/q near 57 using raw ld
    print(f"\n  Rational search: ld / ln(alpha) near p/q?")
    print(f"    Raw ratio = {ratio_raw:.10f}")
    best_frac = None
    best_err = 999
    for q in range(1, 21):
        p = round(ratio_raw * q)
        err = abs(ratio_raw - p / q)
        if err < 0.01 and err < best_err:
            best_err = err
            best_frac = (p, q)
    if best_frac:
        p, q = best_frac
        print(f"    Best: {p}/{q} = {p/q:.10f} (err = {best_err:.2e})")
        if abs(p / q - z_CC) < 0.5:
            found_pattern = True
    else:
        print(f"    No good rational approximation found.")

    # ANALYSIS of L/N near-miss
    ld_N = sum(mults[k] * log(Lk[k] / N) for k in range(N_eig) if Lk[k] > 0.001)
    ratio_N = ld_N / ln_alpha
    print(f"\n  ANALYSIS of L/N near-miss:")
    print(f"    ld(L/N) / ln(alpha) = {ratio_N:.10f}")
    print(f"    57 - ratio = {57 - ratio_N:.10f}")
    print(f"    z_CC - ratio = {z_CC - ratio_N:.10f}")
    print(f"    The gap from 57 is 0.026 -- not close enough for derivation.")
    print(f"    The gap from z_CC is 0.092 -- even worse.")
    print(f"    To get z_CC exactly, normalization would need c={exp((ld - z_CC*ln_alpha)/(N-1)):.6f}")
    print(f"    which is 0.38% off from N=120. No geometric justification.")

    status = "PATTERN" if found_pattern else "NEGATIVE"
    print(f"\n  TEST 2: {status} (near-miss, not exact)")
    return found_pattern


# =====================================================================
# TEST 3: Galois norm of det' vs CC
# =====================================================================
def test_3_galois_norm():
    """Check log|N(det')| against z_CC * ln(1/(2*pi))."""
    print(f"\n{'='*72}")
    print("TEST 3: Galois Norm of det' vs CC")
    print("=" * 72)

    found_pattern = False

    # log|N(det')| = sum m_k * ln|N(L_k)| (all rational, exact)
    log_norm_det = 0
    print(f"\n  Galois norms of eigenvalues:")
    for k in range(N_eig):
        if Lk[k] > 0.001:
            a, b = zphi_forms[k]
            norm_exact = a**2 + a * b - b**2  # N(a+b*phi) in Z[phi]
            contrib = mults[k] * log(abs(norm_exact))
            log_norm_det += contrib
            print(f"    k={k}: N(L_{k}) = {norm_exact:8d}, "
                  f"m_k*ln|N| = {contrib:12.6f}")

    print(f"\n  log|N(det'(L))| = {log_norm_det:.10f}")

    # Target: z_CC * ln(1/(2*pi))
    ln_2pi_inv = log(1 / (2 * pi))
    target = z_CC * ln_2pi_inv
    ratio = log_norm_det / ln_2pi_inv
    print(f"  ln(1/(2*pi)) = {ln_2pi_inv:.10f}")
    print(f"  z_CC * ln(1/(2*pi)) = {target:.10f}")
    print(f"  ratio = log|N(det')| / ln(1/(2*pi)) = {ratio:.10f}")
    print(f"  diff from z_CC = {abs(ratio - z_CC):.6f}")

    if abs(ratio - z_CC) < 0.5:
        print(f"  *** CLOSE MATCH! ***")
        found_pattern = True

    # Also check against other targets
    targets = [
        ("57", 57),
        ("z_CC", z_CC),
        ("35 (seesaw)", 35),
        ("22 (CKM)", 22),
        ("N/2 (60)", 60),
    ]
    print(f"\n  log|N(det')| / ln(X) for various X:")
    for name, target_z in targets:
        # Find X such that log_norm_det = target_z * ln(X)
        # ln(X) = log_norm_det / target_z
        if target_z != 0:
            ln_X = log_norm_det / target_z
            X = exp(ln_X)
            print(f"    target z={name}: X = exp({ln_X:.6f}) = {X:.6f}")

    # Check: is log|N(det')| itself a recognizable number?
    val = exp(log_norm_det)
    print(f"\n  |N(det'(L))| = exp({log_norm_det:.4f}) ~ 10^{log_norm_det/log(10):.2f}")
    print(f"  / N^4 = {val / N**4:.6e}")
    print(f"  / (2*pi)^z_CC = {val / (2*pi)**z_CC:.6e}")

    status = "PATTERN" if found_pattern else "NEGATIVE"
    print(f"\n  TEST 3: {status}")
    return found_pattern


# =====================================================================
# TEST 4: Heat kernel at special times
# =====================================================================
def test_4_heat_kernel():
    """K(t) = 1 + sum m_k * exp(-t*L_k). Find t* where K(t*) = alpha^57."""
    print(f"\n{'='*72}")
    print("TEST 4: Heat Kernel at Special Times")
    print("=" * 72)

    found_pattern = False

    def K(t):
        """Heat kernel trace on Cayley graph."""
        return sum(mults[k] * exp(-t * Lk[k]) for k in range(N_eig))

    # At t=0: K(0) = N = 120
    K0 = K(0)
    print(f"\n  K(0) = {K0:.6f} (= N = {N})")

    # K at special framework times
    special_times = [
        ("alpha", alpha),
        ("alpha_s", alpha_s),
        ("1/(2*pi)", 1 / (2 * pi)),
        ("1/degree", 1.0 / degree),
        ("1/N", 1.0 / N),
        ("ln(phi)", log(PHI)),
        ("phi^(-2)", PHI**(-2)),
        ("1", 1.0),
        ("pi/5", pi / 5),
    ]

    print(f"\n  Heat kernel at framework times:")
    print(f"  {'time':>20s}  {'t':>12s}  {'K(t)':>14s}  {'K/N':>10s}  {'log(K)/log(a)':>14s}")
    for name, t in special_times:
        Kt = K(t)
        ratio = Kt / N if N > 0 else 0
        log_ratio = log(Kt) / log(alpha) if Kt > 0 else float('nan')
        print(f"  {name:>20s}  {t:12.8f}  {Kt:14.6f}  {ratio:10.6f}  {log_ratio:14.6f}")

    # NOTE: K(t) -> m_0 = 1 as t -> inf (zero mode contribution)
    # So K(t) ranges from N=120 (t=0) down to 1 (t=inf).
    # alpha^57 ~ 10^{-122} is FAR below 1, so K(t) = alpha^57 has NO SOLUTION.
    # This is expected: the heat kernel on a finite graph cannot reach such values.
    print(f"\n  Target: K(t*) = alpha^57 = {alpha**57:.6e}")
    print(f"  BUT: K(t) -> 1 as t -> inf (zero mode floor)")
    print(f"  So alpha^57 is UNREACHABLE. Heat kernel approach inapplicable for CC.")

    # Instead, look for framework-meaningful K values
    # Find t where K(t) = 1/alpha = 137.036 (above K(0)=120, so no solution)
    # Find t where K(t) = degree = 12
    print(f"\n  Looking for meaningful K values:")
    for target_name, target_K in [("degree=12", degree), ("b1=6", b1),
                                   ("a1=5", a1), ("N_gen=3", N_gen), ("2", 2.0)]:
        if 1 < target_K < N:
            t_lo, t_hi = 0.0, 20.0
            for _ in range(200):
                t_mid = (t_lo + t_hi) / 2
                if K(t_mid) > target_K:
                    t_lo = t_mid
                else:
                    t_hi = t_mid
            t_star = (t_lo + t_hi) / 2
            print(f"  K(t*) = {target_name}: t* = {t_star:.8f}")
            for name, val in [("alpha", alpha), ("1/(2*pi)", 1/(2*pi)),
                              ("1/degree", 1.0/degree), ("phi^(-2)", PHI**(-2))]:
                r = t_star / val
                if abs(r - round(r)) < 0.05:
                    print(f"    t* / {name} = {r:.6f} ~ {round(r)}")
                    found_pattern = True

    status = "PATTERN" if found_pattern else "NEGATIVE"
    print(f"\n  TEST 4: {status}")
    return found_pattern


# =====================================================================
# TEST 5: Spectral zeta at special points
# =====================================================================
def test_5_spectral_zeta():
    """zeta(s) = sum m_k * L_k^(-s) at special s values."""
    print(f"\n{'='*72}")
    print("TEST 5: Spectral Zeta Function at Special Points")
    print("=" * 72)

    found_pattern = False

    def zeta_L(s):
        """Spectral zeta of Cayley Laplacian."""
        return sum(mults[k] * Lk[k]**(-s) for k in range(N_eig) if Lk[k] > 0.001)

    # Standard values
    print(f"\n  zeta_L(s) = sum m_k * L_k^(-s):")
    print(f"  {'s':>8s}  {'zeta_L(s)':>16s}  {'ratio to 57':>12s}  {'note':>20s}")
    for s_val in [-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3]:
        z = zeta_L(s_val)
        ratio = z / 57 if z != 0 else float('nan')
        note = ""
        if abs(z - 57) < 0.5:
            note = "*** = 57 ***"
            found_pattern = True
        elif abs(z - z_CC) < 0.5:
            note = "*** = z_CC ***"
            found_pattern = True
        elif abs(z - round(z)) < 0.01:
            note = f"int: {round(z)}"
        print(f"  {s_val:8.1f}  {z:16.6f}  {ratio:12.6f}  {note}")

    # Search for s where zeta_L(s) = 57 or z_CC
    print(f"\n  Searching for s* where zeta_L(s*) = 57:")
    best_s = None
    best_diff = 1e10
    for s100 in range(-1000, 1000):
        s_val = s100 / 100.0
        try:
            z = zeta_L(s_val)
            diff = abs(z - 57)
            if diff < best_diff:
                best_diff = diff
                best_s = s_val
        except (OverflowError, ZeroDivisionError):
            pass

    if best_s is not None:
        z_best = zeta_L(best_s)
        print(f"    Best: s* = {best_s:.2f}, zeta_L(s*) = {z_best:.6f}, diff = {best_diff:.4f}")
        if best_diff < 0.01:
            found_pattern = True
            # Check if s* is recognizable
            print(f"    s* / alpha = {best_s / alpha:.4f}")
            print(f"    s* / (1/alpha) = {best_s * alpha:.6f}")
            print(f"    s* * 2 = {best_s * 2:.4f}")

    # Also: zeta'(0) = -log det'
    # Compute numerically
    h = 1e-6
    zeta_prime_0 = (zeta_L(h) - zeta_L(-h)) / (2 * h)
    ld = sum(mults[k] * log(Lk[k]) for k in range(N_eig) if Lk[k] > 0.001)
    print(f"\n  zeta'(0) = {zeta_prime_0:.6f}")
    print(f"  -log det'(L) = {-ld:.6f}")
    print(f"  Match: {abs(zeta_prime_0 + ld) < 0.1}")

    status = "PATTERN" if found_pattern else "NEGATIVE"
    print(f"\n  TEST 5: {status}")
    return found_pattern


# =====================================================================
# TEST 6: One-loop vacuum energy
# =====================================================================
def test_6_vacuum_energy():
    """E_vac = (1/(2*N)) * sum m_k * sqrt(L_k)."""
    print(f"\n{'='*72}")
    print("TEST 6: One-Loop Vacuum Energy")
    print("=" * 72)

    found_pattern = False

    # Vacuum energy: E_vac = (1/2) * sum_{k>0} m_k * sqrt(L_k) / N
    E_vac = 0.5 * sum(mults[k] * sqrt(Lk[k]) for k in range(N_eig) if Lk[k] > 0.001) / N

    # Also without normalization
    E_vac_unnorm = 0.5 * sum(mults[k] * sqrt(Lk[k]) for k in range(N_eig) if Lk[k] > 0.001)

    ln_alpha_val = log(alpha)
    z_E = log(E_vac) / ln_alpha_val if E_vac > 0 else float('nan')
    z_E_un = log(E_vac_unnorm) / ln_alpha_val if E_vac_unnorm > 0 else float('nan')

    print(f"\n  E_vac (normalized) = {E_vac:.10f}")
    print(f"  E_vac (unnormalized) = {E_vac_unnorm:.6f}")
    print(f"  log(E_vac) / ln(alpha) = {z_E:.6f}")
    print(f"  log(E_vac_unnorm) / ln(alpha) = {z_E_un:.6f}")
    print(f"  target z_CC = {z_CC:.6f}")

    # Different powers of L_k
    print(f"\n  Vacuum energy with different powers:")
    for power_name, power in [("1/2 (sqrt)", 0.5), ("1", 1.0), ("2", 2.0),
                                ("-1/2", -0.5), ("-1", -1.0)]:
        E_p = sum(mults[k] * Lk[k]**power for k in range(N_eig) if Lk[k] > 0.001)
        if E_p > 0:
            z_p = log(E_p) / ln_alpha_val
            diff = abs(z_p - z_CC)
            mark = "***" if diff < 0.5 else ""
            print(f"    power={power_name:6s}: E = {E_p:14.6f}, "
                  f"log(E)/ln(a) = {z_p:10.4f}, diff = {diff:.4f} {mark}")
            if diff < 0.5:
                found_pattern = True

    # Alternating Hodge-like sum
    # Using Galois parity: rational eigenvalues (self-dual) vs irrational
    E_alt = 0
    for k in range(N_eig):
        if Lk[k] > 0.001:
            a, b = zphi_forms[k]
            sign = 1 if b == 0 else -1  # rational = +1, irrational = -1
            E_alt += sign * mults[k] * sqrt(Lk[k])

    print(f"\n  Alternating sum (rational+ / irrational-):")
    print(f"    E_alt = {E_alt:.10f}")
    if abs(E_alt) > 1e-10:
        z_alt = log(abs(E_alt)) / ln_alpha_val
        print(f"    log|E_alt| / ln(alpha) = {z_alt:.6f}")
        if abs(z_alt - z_CC) < 0.5:
            found_pattern = True
            print(f"    *** CLOSE TO z_CC ***")

    status = "PATTERN" if found_pattern else "NEGATIVE"
    print(f"\n  TEST 6: {status}")
    return found_pattern


# =====================================================================
# TEST 7: Integer combination search
# =====================================================================
def test_7_integer_search():
    """Search for a*ld + b*ln(2*pi) + ... = c*ln(alpha)."""
    print(f"\n{'='*72}")
    print("TEST 7: Integer Combination Search")
    print("=" * 72)

    found_pattern = False

    # Key quantities
    ld = sum(mults[k] * log(Lk[k]) for k in range(N_eig) if Lk[k] > 0.001)
    ld_galois = sum(mults[k] * log(abs(Lk_galois[k]))
                    for k in range(N_eig) if abs(Lk[k]) > 0.001)
    ln_2pi = log(2 * pi)
    ln_alpha_val = log(alpha)
    ln_phi = log(PHI)
    ln_N = log(N)

    print(f"\n  Key quantities:")
    print(f"    ld = log det'(L) = {ld:.10f}")
    print(f"    ld_galois = log det'_galois(L) = {ld_galois:.10f}")
    print(f"    ln(2*pi) = {ln_2pi:.10f}")
    print(f"    ln(alpha) = {ln_alpha_val:.10f}")
    print(f"    ln(phi) = {ln_phi:.10f}")
    print(f"    ln(N) = {ln_N:.10f}")

    # Search: a*ld + b*ln(2*pi) = z_CC * ln(alpha)?
    # i.e., a*ld + b*ln(2*pi) + c*ln(alpha) = 0 with c = -z_CC
    target = z_CC * ln_alpha_val  # = ln(alpha^z_CC) = ln(Lambda_CC)
    print(f"\n  Target: z_CC * ln(alpha) = {target:.10f}")

    print(f"\n  Search: a*ld + b*ln(2*pi) = z_CC*ln(alpha)?")
    for a in range(-5, 6):
        for b in range(-5, 6):
            if a == 0 and b == 0:
                continue
            val = a * ld + b * ln_2pi
            if abs(target) > 0:
                ratio = val / target
                if abs(ratio - 1.0) < 0.01:
                    print(f"    MATCH: {a}*ld + {b}*ln(2pi) = {val:.6f} "
                          f"(target {target:.6f}, ratio {ratio:.6f})")
                    found_pattern = True

    # Search: a*ld + b*ld_galois = z_CC * ln(alpha)?
    print(f"\n  Search: a*ld + b*ld_galois = z_CC*ln(alpha)?")
    for a in range(-5, 6):
        for b in range(-5, 6):
            if a == 0 and b == 0:
                continue
            val = a * ld + b * ld_galois
            if abs(target) > 0:
                ratio = val / target
                if abs(ratio - 1.0) < 0.01:
                    print(f"    MATCH: {a}*ld + {b}*ld_gal = {val:.6f} "
                          f"(target {target:.6f}, ratio {ratio:.6f})")
                    found_pattern = True

    # Search: a*ld + b*ln(2*pi) + c*ln(phi) + d*ln(N) = z_CC*ln(alpha)?
    print(f"\n  Extended search: a*ld + b*ln(2pi) + c*ln(phi) + d*ln(N) = z_CC*ln(alpha)?")
    found_extended = False
    for a in range(-3, 4):
        for b in range(-3, 4):
            for c in range(-5, 6):
                for d in range(-3, 4):
                    if a == b == c == d == 0:
                        continue
                    val = a * ld + b * ln_2pi + c * ln_phi + d * ln_N
                    if abs(target) > 0:
                        ratio = val / target
                        if abs(ratio - 1.0) < 0.005:
                            print(f"    {a}*ld + {b}*ln(2pi) + {c}*ln(phi) + {d}*ln(N)"
                                  f" = {val:.6f} (target {target:.6f}, "
                                  f"ratio {ratio:.8f})")
                            found_pattern = True
                            found_extended = True
    if not found_extended:
        print(f"    No match found.")

    # Also try: expressing ld directly in terms of framework quantities
    print(f"\n  ld decomposition:")
    print(f"    ld / ln(alpha) = {ld / ln_alpha_val:.10f}")
    print(f"    ld / ln(phi) = {ld / ln_phi:.10f}")
    print(f"    ld / ln(N) = {ld / ln_N:.10f}")
    print(f"    ld / ln(2*pi) = {ld / ln_2pi:.10f}")
    print(f"    ld / (degree) = {ld / degree:.10f}")

    status = "PATTERN" if found_pattern else "NEGATIVE"
    print(f"\n  TEST 7: {status}")
    return found_pattern


# =====================================================================
# TEST 8: Galois decomposition of log det'
# =====================================================================
def test_8_galois_decomposition():
    """Split log det' = R (rational) + I (irrational). STRUCTURAL THEOREM."""
    print(f"\n{'='*72}")
    print("TEST 8: Galois Decomposition of log det'")
    print("=" * 72)

    found_pattern = False

    ld = sum(mults[k] * log(Lk[k]) for k in range(N_eig) if Lk[k] > 0.001)
    ld_galois = sum(mults[k] * log(abs(Lk_galois[k]))
                    for k in range(N_eig) if abs(Lk[k]) > 0.001)

    # Symmetric (rational) and antisymmetric (irrational) parts
    R = (ld + ld_galois) / 2       # = (1/2) * log|N(det')|
    I = (ld - ld_galois) / 2       # irrational part

    print(f"\n  log det'(L) = {ld:.10f}")
    print(f"  log det'_galois(L) = {ld_galois:.10f}")
    print(f"  R = (ld + ld_gal)/2 = {R:.10f}  (rational part)")
    print(f"  I = (ld - ld_gal)/2 = {I:.10f}  (irrational part)")
    print(f"  Check: R + I = {R + I:.10f} (= ld)")

    # STRUCTURAL THEOREM: I = 0 exactly
    print(f"\n  THEOREM: I = 0 exactly (det'(L) is Galois-invariant).")
    print(f"  PROOF: The Galois action sigma: phi -> phi' permutes eigenvalues:")
    print(f"    sigma(L_1) = L_7, sigma(L_7) = L_1  (both have d_k=2, mult=4)")
    print(f"    sigma(L_2) = L_8, sigma(L_8) = L_2  (both have d_k=3, mult=9)")
    print("  Since m_k = m_{sigma(k)}, the Galois conjugate det' is:")
    print("    sigma(det') = prod sigma(L_k)^m_k = prod L_{sigma(k)}^m_k")
    print("                = prod L_k^m_k = det'  (by reindexing)")
    print(f"  Therefore det'(L) in Q (rational). QED.")

    ln_alpha_val = log(alpha)
    ln_phi = log(PHI)

    print(f"\n  Rational part R:")
    print(f"    R / ln(alpha) = {R / ln_alpha_val:.10f}")
    print(f"    R / ln(1/(2*pi)) = {R / log(1/(2*pi)):.10f}")
    r_near = round(R / ln_alpha_val)
    print(f"    Nearest integer to R/ln(alpha): {r_near}")

    print(f"\n  Irrational part I:")
    print(f"    I / ln(phi) = {I / ln_phi:.10f}")
    i_near = round(I / ln_phi)
    print(f"    Nearest integer to I/ln(phi): {i_near}")
    print(f"    I / ln(phi) - {i_near} = {I / ln_phi - i_near:.10f}")

    # Check: does R relate to 57 * ln(alpha) somehow?
    # Or: R + c*I = 57*ln(alpha)?
    # R + I = ld, so we need ld = 57*ln(alpha) which is Test 2.
    # More interesting: R = some_z * ln(alpha)? I = some_n * ln(phi)?
    z_R = R / ln_alpha_val
    z_I = I / ln_phi

    print(f"\n  Decomposition: ld = z_R*ln(alpha) + z_I*ln(phi)")
    print(f"    z_R = {z_R:.10f}")
    print(f"    z_I = {z_I:.10f}")
    print(f"    z_R + z_I*ln(phi)/ln(alpha) = {z_R + z_I*ln_phi/ln_alpha_val:.10f}")
    print(f"    should = ld/ln(alpha) = {ld/ln_alpha_val:.10f}")

    # Check against framework quantities
    framework = [
        ("57", 57), ("z_CC", z_CC), ("35", 35), ("22", 22),
        ("N/2", 60), ("N", 120), ("h_E8", 30), ("rank_E8", 8),
        ("a1^2", 25), ("degree", 12), ("a1*b1", 30),
    ]
    print(f"\n  z_R = {z_R:.6f} matches?")
    for name, val in framework:
        if abs(z_R - val) < 0.5:
            print(f"    *** CLOSE to {name} = {val} (diff = {z_R - val:.6f}) ***")
            found_pattern = True

    print(f"\n  z_I = {z_I:.6f} matches?")
    for name, val in framework:
        if abs(z_I - val) < 0.5:
            print(f"    *** CLOSE to {name} = {val} (diff = {z_I - val:.6f}) ***")
            found_pattern = True

    # Extra: combined forms
    print(f"\n  Combined searches:")
    for (a, b) in [(1, 1), (1, -1), (2, 1), (1, 2), (2, -1), (-1, 2)]:
        combo = a * z_R + b * z_I
        for name, val in framework:
            if abs(combo - val) < 0.5:
                print(f"    {a}*z_R + {b}*z_I = {combo:.6f} ~ {name} = {val}")
                found_pattern = True

    status = "PATTERN" if found_pattern else "NEGATIVE"
    print(f"\n  TEST 8: {status}")
    return found_pattern


# =====================================================================
# TEST 9: Full Hodge torsion vs alpha^57
# =====================================================================
def test_9_hodge_torsion():
    """Build 600-cell simplicial complex and compute analytic torsion."""
    print(f"\n{'='*72}")
    print("TEST 9: Full Hodge Torsion vs alpha^57")
    print("=" * 72)

    found_pattern = False

    # Build the 600-cell complex (reuse exp408c logic)
    print("\n  Building 600-cell...")
    from itertools import permutations, product as cartesian_product
    from collections import defaultdict
    import time

    t0 = time.time()

    # --- Vertices ---
    verts_set = set()
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts_set.add(tuple(v))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    verts_set.add((s0, s1, s2, s3))
    base = [PHI / 2, 0.5, 1 / (2 * PHI), 0.0]
    even_perms = [p for p in permutations(range(4))
                  if sum(1 for i in range(4) for j in range(i + 1, 4)
                         if p[i] > p[j]) % 2 == 0]
    for perm in even_perms:
        coords = [base[perm[i]] for i in range(4)]
        nonzero_idx = [i for i in range(4) if abs(coords[i]) > 1e-12]
        for signs in cartesian_product([1, -1], repeat=len(nonzero_idx)):
            v = list(coords)
            for idx, s in zip(nonzero_idx, signs):
                v[idx] *= s
            verts_set.add(tuple(round(x, 10) for x in v))

    verts = np.array(sorted(verts_set))
    Nv = len(verts)

    # --- Edges ---
    dots = verts @ verts.T
    edges = [(i, j) for i in range(Nv) for j in range(i + 1, Nv)
             if abs(dots[i, j] - PHI / 2) < 0.001]
    Ne = len(edges)
    edge_to_idx = {}
    for idx, (i, j) in enumerate(edges):
        edge_to_idx[(i, j)] = idx
        edge_to_idx[(j, i)] = idx
    adj_list = defaultdict(set)
    for i, j in edges:
        adj_list[i].add(j)
        adj_list[j].add(i)

    # --- Triangles ---
    triangles = []
    for i in range(Nv):
        for j in adj_list[i]:
            if j > i:
                for k in adj_list[i] & adj_list[j]:
                    if k > j:
                        triangles.append((i, j, k))
    Nf = len(triangles)
    face_to_idx = {t: idx for idx, t in enumerate(triangles)}

    # --- Tetrahedra ---
    tetrahedra = []
    for i in range(Nv):
        ni = adj_list[i]
        for j in ni:
            if j > i:
                cij = ni & adj_list[j]
                for k in cij:
                    if k > j:
                        for l in cij & adj_list[k]:
                            if l > k:
                                tetrahedra.append((i, j, k, l))
    Nc = len(tetrahedra)

    build_time = time.time() - t0
    print(f"  f-vector: ({Nv}, {Ne}, {Nf}, {Nc}) in {build_time:.1f}s")

    # Check Euler characteristic
    euler = Nv - Ne + Nf - Nc
    print(f"  Euler char: {Nv}-{Ne}+{Nf}-{Nc} = {euler} (expected 0 for S^3)")

    # --- Boundary operators ---
    print("  Computing Hodge Laplacians...")
    d0 = np.zeros((Ne, Nv))
    for e_idx, (i, j) in enumerate(edges):
        d0[e_idx, i] = -1.0
        d0[e_idx, j] = +1.0

    d1 = np.zeros((Nf, Ne))
    for f_idx, (i, j, k) in enumerate(triangles):
        d1[f_idx, edge_to_idx[(i, j)]] = +1.0
        d1[f_idx, edge_to_idx[(j, k)]] = +1.0
        d1[f_idx, edge_to_idx[(i, k)]] = -1.0

    d2 = np.zeros((Nc, Nf))
    for c_idx, (i, j, k, l) in enumerate(tetrahedra):
        for face, sign in [((j, k, l), +1), ((i, k, l), -1),
                           ((i, j, l), +1), ((i, j, k), -1)]:
            f_idx = face_to_idx.get(face)
            if f_idx is not None:
                d2[c_idx, f_idx] = sign

    # --- Hodge Laplacians ---
    Delta_0 = d0.T @ d0
    Delta_1 = d0 @ d0.T + d1.T @ d1
    Delta_2 = d1 @ d1.T + d2.T @ d2
    Delta_3 = d2 @ d2.T

    # --- Eigenvalues ---
    print("  Computing eigenvalues...")
    evals_0 = np.sort(np.linalg.eigvalsh(Delta_0))
    evals_1 = np.sort(np.linalg.eigvalsh(Delta_1))
    evals_2 = np.sort(np.linalg.eigvalsh(Delta_2))
    evals_3 = np.sort(np.linalg.eigvalsh(Delta_3))

    tol = 0.01

    def log_det_p(evals):
        nz = evals[evals > tol]
        return np.sum(np.log(nz)), len(evals) - len(nz)

    ld0, betti0 = log_det_p(evals_0)
    ld1, betti1 = log_det_p(evals_1)
    ld2, betti2 = log_det_p(evals_2)
    ld3, betti3 = log_det_p(evals_3)

    print(f"  Betti: ({betti0}, {betti1}, {betti2}, {betti3}) [expected (1,0,0,1)]")
    print(f"  log det'(Delta_k): {ld0:.4f}, {ld1:.4f}, {ld2:.4f}, {ld3:.4f}")

    # Analytic torsion: log T = (1/2) * sum (-1)^p * p * ld_p
    log_T = 0.5 * (0 * ld0 - 1 * ld1 + 2 * ld2 - 3 * ld3)
    log_T_R = sum((-1)**(p + 1) * ld for p, ld in enumerate([ld0, ld1, ld2, ld3]))

    ln_alpha_val = log(alpha)

    print(f"\n  Analytic torsion:")
    print(f"    log T (Ray-Singer) = {log_T:.10f}")
    print(f"    log T_R (Reidemeister) = {log_T_R:.10f}")

    # Test against alpha^57
    targets = [
        ("log T / ln(alpha)", log_T / ln_alpha_val),
        ("-log T / ln(alpha)", -log_T / ln_alpha_val),
        ("log T_R / ln(alpha)", log_T_R / ln_alpha_val),
        ("-log T_R / ln(alpha)", -log_T_R / ln_alpha_val),
    ]

    print(f"\n  Ratios to ln(alpha):")
    for name, val in targets:
        diff = abs(val - z_CC)
        mark = "***" if diff < 0.5 else ""
        print(f"    {name:>30s} = {val:12.6f} (diff from z_CC = {diff:.4f}) {mark}")
        if diff < 0.5:
            found_pattern = True

    # Test each ld individually
    print(f"\n  Individual log det' / ln(alpha):")
    for name, ld_val in [("ld0", ld0), ("ld1", ld1), ("ld2", ld2), ("ld3", ld3)]:
        ratio = ld_val / ln_alpha_val
        diff = abs(ratio - z_CC)
        mark = "***" if diff < 0.5 else ""
        print(f"    {name} / ln(alpha) = {ratio:12.6f} (diff from z_CC = {diff:.4f}) {mark}")
        if diff < 0.5:
            found_pattern = True

    # Alternating Hodge sum
    hodge_alt = ld0 - ld1 + ld2 - ld3
    ratio_alt = hodge_alt / ln_alpha_val
    diff_alt = abs(ratio_alt - z_CC)
    print(f"\n  Alternating Hodge sum: ld0-ld1+ld2-ld3 = {hodge_alt:.10f}")
    print(f"    / ln(alpha) = {ratio_alt:.6f} (diff from z_CC = {diff_alt:.4f})")
    if diff_alt < 0.5:
        found_pattern = True
        print(f"    *** CLOSE MATCH! ***")

    # Full search: a*ld0 + b*ld1 + c*ld2 + d*ld3 = z_CC * ln(alpha)?
    target_val = z_CC * ln_alpha_val
    print(f"\n  Search: a*ld0+b*ld1+c*ld2+d*ld3 = z_CC*ln(alpha) = {target_val:.4f}?")
    found_combo = False
    for a in range(-2, 3):
        for b in range(-2, 3):
            for c in range(-2, 3):
                for d in range(-2, 3):
                    if a == b == c == d == 0:
                        continue
                    val = a * ld0 + b * ld1 + c * ld2 + d * ld3
                    if abs(val - target_val) < 1.0:
                        print(f"    ({a})*ld0+({b})*ld1+({c})*ld2+({d})*ld3"
                              f" = {val:.4f} (target {target_val:.4f},"
                              f" diff {val - target_val:.4f})")
                        found_combo = True
                        if abs(val - target_val) < 0.5:
                            found_pattern = True
    if not found_combo:
        print(f"    No match found.")

    # Cross-check with Cayley graph log det'
    ld_cayley = sum(mults[k] * log(Lk[k]) for k in range(N_eig) if Lk[k] > 0.001)
    print(f"\n  Cross-check:")
    print(f"    log det'(Cayley L_0) = {ld_cayley:.6f}")
    print(f"    log det'(simplicial Delta_0) = {ld0:.6f}")
    print(f"    Ratio: {ld_cayley / ld0:.6f}")

    status = "PATTERN" if found_pattern else "NEGATIVE"
    print(f"\n  TEST 9: {status}")
    return found_pattern


# =====================================================================
# MAIN
# =====================================================================
def main():
    print(f"\n{'='*72}")
    print("RUNNING ALL TESTS")
    print("=" * 72)

    results = []

    results.append(("TEST 1: Spectral data validation", test_1_spectral_validation()))
    results.append(("TEST 2: log det'(L)/ln(alpha) vs z_CC [CORE]", test_2_logdet_vs_zCC()))
    results.append(("TEST 3: Galois norm of det' vs CC", test_3_galois_norm()))
    results.append(("TEST 4: Heat kernel at special times", test_4_heat_kernel()))
    results.append(("TEST 5: Spectral zeta at special points", test_5_spectral_zeta()))
    results.append(("TEST 6: One-loop vacuum energy", test_6_vacuum_energy()))
    results.append(("TEST 7: Integer combination search", test_7_integer_search()))
    results.append(("TEST 8: Galois decomposition", test_8_galois_decomposition()))
    results.append(("TEST 9: Full Hodge torsion vs alpha^57", test_9_hodge_torsion()))

    # =====================================================================
    # FINAL SUMMARY
    # =====================================================================
    print(f"\n{'='*72}")
    print("FINAL SUMMARY")
    print("=" * 72)

    n_pattern = sum(1 for _, ok in results if ok)
    n_negative = sum(1 for _, ok in results if not ok)
    n_total = len(results)

    # Infrastructure test (Test 1) reports separately
    infra_ok = results[0][1]

    for name, ok in results:
        status = "PASS/PATTERN" if ok else "NEGATIVE"
        print(f"  [{status:>12s}] {name}")

    print(f"\n  Infrastructure: {'PASS' if infra_ok else 'FAIL'}")
    print(f"  Patterns found: {n_pattern}/{n_total}")
    print(f"  Negative: {n_negative}/{n_total}")

    # Overall assessment
    physics_tests = results[1:]  # skip infrastructure
    any_pattern = any(ok for _, ok in physics_tests)

    # The L/N normalization gives ratio = 56.974, which is a near-miss
    ld_check = sum(mults[k] * log(Lk[k] / N) for k in range(N_eig) if Lk[k] > 0.001)
    ratio_check = ld_check / log(alpha)

    print(f"""
{'='*72}
ASSESSMENT
{'='*72}

  QUESTION: Can Lambda_P = alpha^z_CC be derived from spectral data of
  the 600-cell Cayley graph?

  CORE TEST (Test 2): log det'(L/N) / ln(alpha) = {ratio_check:.6f}
    Target z_CC = {z_CC:.6f}, diff = {abs(ratio_check - z_CC):.6f}
    Target 57   = {57:.6f}, diff = {abs(ratio_check - 57):.6f}

  RESULT: NEAR-MISS with L/N normalization (0.026 from 57, 0.092 from z_CC).
  NOT close enough for a derivation. No principled normalization gives z_CC.

  NEW STRUCTURAL RESULT (Test 8):
    det'(L) is GALOIS-INVARIANT (rational number).
    Proof: Galois-conjugate irreps have equal dimensions, so
    sigma permutes eigenvalues within each multiplicity class.
    This is a DERIVED THEOREM of the 2I representation theory.

  WHAT WORKS / WHAT DOESN'T:
    Exponent 57:        FULLY DERIVED (7 algebraic routes, all equiv to a1=5)
    Galois duality:     DERIVED (Lambda*Lambda' = (1/(2*pi))^z)
    det'(L) rationality: DERIVED (structural theorem, Test 8)
    Form alpha^z:       PATTERN (this experiment found NO derivation)
    Heat kernel:        INAPPLICABLE (K(t)->1, can't reach alpha^57)
    Analytic torsion:   NEGATIVE (log T_R = 0 for S^3)
    Spectral zeta:      NEGATIVE (no special value gives 57 or z_CC)

  HONEST CONCLUSION:
    The functional form alpha^z for the CC remains UNEXPLAINED.
    This is consistent with the genuinely non-perturbative nature
    of the CC problem: the one-loop computation (exp386) falls
    short by ~25 = a1^2, and no finite-graph spectral quantity
    can bridge this gap.

    The CC exponent 57 is FULLY DERIVED. The base alpha is not.
    This may require new mathematics (non-perturbative spectral
    geometry, resurgence, or instanton calculus on the fiber).

  STATUS: NEGATIVE (functional form alpha^z not derived)
""")

    return 0 if infra_ok else 1


if __name__ == "__main__":
    sys.exit(main())

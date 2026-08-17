"""
exp425_partition_function_cc.py
===============================
GOAL: Compute the FULL partition function Z for a scalar field on the 600-cell
and check whether ln Z / Vol(S^3) relates to alpha^57 (cosmological constant).

CONTEXT:
  - exp421 computed det'(L) and found log det'(L/N)/ln(alpha) = 56.97 (near-miss)
  - BUT: the full partition function Z = (2*pi)^{(N-1)/2} / sqrt(det'(L))
    was never computed with proper QFT normalization
  - The CC problem: Lambda_Planck ~ alpha^57. Can Z on 600-cell produce this?

WHAT'S NEW vs exp421:
  1. Full Z with (2*pi) factors
  2. Free energy F = -ln Z
  3. Vacuum energy density rho = F / Vol(S^3)
  4. Kaluza-Klein reduction: Lambda_4D = one-loop / V_internal
  5. Compact scalar (periodic BCs)
  6. Mass term exploration: Z(m^2) = prod (L_k + m^2)^{-m_k/2}
  7. Ratio Z_600 / Z_reference (relative partition function)

EIGENVALUE DATA (from exp421, verified):
  k  d_k  mult  L_k
  0   1    1    0                 (zero mode)
  1   2    4    12 - 6*phi
  2   3    9    12 - 4*phi
  3   4   16    9
  4   5   25    12
  5   6   36    14
  6   4   16    15
  7   2    4    6 + 6*phi
  8   3    9    8 + 4*phi

BLIND: compute first, interpret after.
Author: Razvan-Constantin Anghelina
Date: 2026-02-27
"""

import numpy as np
from math import sqrt, pi, log, log10, exp, factorial
from fractions import Fraction
import sys

print("=" * 72)
print("EXP-425: PARTITION FUNCTION ON 600-CELL vs COSMOLOGICAL CONSTANT")
print("        Full QFT computation: Z, F, rho_vac")
print("=" * 72)

# =====================================================================
# FRAMEWORK CONSTANTS
# =====================================================================
a1 = 5
b1 = a1 + 1                              # 6
PHI = (1 + sqrt(a1)) / 2                  # golden ratio
PHI_prime = (1 - sqrt(a1)) / 2            # Galois conjugate
N = factorial(a1)                          # 120
degree = 2 * b1                           # 12
N_eig = 9
rank_E8 = 8
h_E8 = a1 * b1                            # 30
N_gen = 3

# Couplings
alpha_s = 1 / (2 * PHI**3)
sin2tW = Fraction(b1, a1**2 + 1)           # 6/26

# Alpha from quadratic
A_coef = 2 * pi
B_coef = -4 * a1 * PHI**4
C_coef = 1
disc = B_coef**2 - 4 * A_coef * C_coef
alpha = (-B_coef - sqrt(disc)) / (2 * A_coef)
alpha_prime = (-B_coef + sqrt(disc)) / (2 * A_coef)

# CC exponent
z_CC = 57 - alpha_s
ln_alpha = log(alpha)
ln_2pi = log(2 * pi)

# Volumes
Vol_S3 = 2 * pi**2    # unit S^3

print(f"\nFramework:")
print(f"  a1={a1}, N={N}, degree={degree}, PHI={PHI:.10f}")
print(f"  alpha = {alpha:.10f} (1/alpha = {1/alpha:.6f})")
print(f"  alpha_s = {alpha_s:.10f}")
print(f"  z_CC = {z_CC:.10f}")
print(f"  alpha^z_CC = {alpha**z_CC:.6e}")
print(f"  ln(alpha) = {ln_alpha:.10f}")
print(f"  Vol(S^3) = 2*pi^2 = {Vol_S3:.10f}")

# =====================================================================
# EIGENVALUE DATA
# =====================================================================
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
chi_10a = [1, PHI, PHI, 1, 0, -1, -1, 1 - PHI, 1 - PHI]
Lk = [degree * (1 - chi_10a[k] / dims[k]) for k in range(N_eig)]
mults = [d**2 for d in dims]

print(f"\nEigenvalues of Cayley graph Laplacian:")
for k in range(N_eig):
    print(f"  k={k}: L_{k} = {Lk[k]:10.6f}, mult = {mults[k]:3d}")

# =====================================================================
# CORE COMPUTATIONS
# =====================================================================

# log det'(L) = sum_{k>0} m_k * ln(L_k)
ld = sum(mults[k] * log(Lk[k]) for k in range(N_eig) if Lk[k] > 1e-10)
N_nonzero = N - 1   # 119 nonzero eigenvalues (counting multiplicity)

print(f"\n{'='*72}")
print("CORE COMPUTATION: Partition Function")
print("=" * 72)

# =====================================================================
# TEST 1: Free scalar - standard normalization
# =====================================================================
print(f"\n--- TEST 1: Free scalar, standard normalization ---")

# Z = (2*pi)^{(N-1)/2} / sqrt(det'(L))
# ln Z = (N-1)/2 * ln(2*pi) - 1/2 * ln det'(L)
ln_Z = N_nonzero / 2 * ln_2pi - 0.5 * ld

print(f"  ln det'(L) = {ld:.10f}")
print(f"  (N-1)/2 * ln(2pi) = {N_nonzero/2 * ln_2pi:.10f}")
print(f"  ln Z = {ln_Z:.10f}")
print(f"  F = -ln Z = {-ln_Z:.10f}")

# Vacuum energy density
rho_vac = -ln_Z / Vol_S3
print(f"\n  rho_vac = F / Vol(S^3) = {rho_vac:.10f}")
print(f"  |rho_vac| = {abs(rho_vac):.10f}")
print(f"  log10|rho_vac| = {log10(abs(rho_vac)):.6f}")

# Compare to alpha^57
target = alpha**z_CC
print(f"\n  Target: alpha^z_CC = {target:.6e}")
print(f"  log10(alpha^z_CC) = {z_CC * log10(alpha):.4f}")
print(f"  |rho_vac| / alpha^z_CC = {abs(rho_vac) / target:.6e}")
print(f"  ln|rho_vac| / ln(alpha) = {log(abs(rho_vac)) / ln_alpha:.6f}")

# =====================================================================
# TEST 2: Normalized Laplacian L/degree
# =====================================================================
print(f"\n--- TEST 2: Normalized Laplacian L_norm = L/degree ---")
# L_norm_k = L_k / degree => det'(L_norm) = det'(L) / degree^{N-1}
ld_norm = ld - N_nonzero * log(degree)
ln_Z_norm = N_nonzero / 2 * ln_2pi - 0.5 * ld_norm

print(f"  ln det'(L/degree) = {ld_norm:.10f}")
print(f"  ln Z_norm = {ln_Z_norm:.10f}")
rho_norm = -ln_Z_norm / Vol_S3
print(f"  rho_norm = {rho_norm:.10f}")
print(f"  ln|rho_norm| / ln(alpha) = {log(abs(rho_norm)) / ln_alpha:.6f}")

# =====================================================================
# TEST 3: Various zero-mode prescriptions
# =====================================================================
print(f"\n--- TEST 3: Zero-mode prescriptions ---")

# The zero mode integral gives a factor V_0 = "volume of field space"
# Different choices:
prescriptions = [
    ("No zero mode (standard)", 0),
    ("V_0 = sqrt(N)", log(sqrt(N))),
    ("V_0 = N", log(N)),
    ("V_0 = 2*pi*sqrt(N) [compact S^1]", log(2*pi*sqrt(N))),
    ("V_0 = Vol(S^3) = 2*pi^2", log(Vol_S3)),
    ("V_0 = 1/(2*pi) [KK radius]", log(1/(2*pi))),
    ("V_0 = sqrt(N/(2*pi))", log(sqrt(N/(2*pi)))),
]

for name, ln_V0 in prescriptions:
    ln_Z_full = ln_Z + ln_V0
    rho_full = -ln_Z_full / Vol_S3
    ratio = log(abs(rho_full)) / ln_alpha if abs(rho_full) > 1e-300 else float('nan')
    print(f"  {name}")
    print(f"    ln Z = {ln_Z_full:.6f}, F/Vol = {rho_full:.6f}, "
          f"ln|F/Vol|/ln(a) = {ratio:.4f}")

# =====================================================================
# TEST 4: Free energy PER DEGREE OF FREEDOM
# =====================================================================
print(f"\n--- TEST 4: Free energy per degree of freedom ---")

F = -ln_Z
f_per_dof = F / N_nonzero
f_per_vertex = F / N

print(f"  F = {F:.10f}")
print(f"  F / (N-1) = {f_per_dof:.10f}  (per mode)")
print(f"  F / N = {f_per_vertex:.10f}  (per vertex)")
print(f"  F / (N * Vol) = {F / (N * Vol_S3):.10f}")

# Check if F itself relates to alpha^57
print(f"\n  ln|F| / ln(alpha) = {log(abs(F)) / ln_alpha:.6f}")
print(f"  ln|F/Vol| / ln(alpha) = {log(abs(F/Vol_S3)) / ln_alpha:.6f}")

# =====================================================================
# TEST 5: Massive scalar - Z(m^2) with framework masses
# =====================================================================
print(f"\n--- TEST 5: Massive scalar Z(m^2) ---")

# Z(m^2) = prod_{k>=0} (L_k + m^2)^{-m_k/2} * (2*pi)^{N/2}
# ln Z(m^2) = N/2 * ln(2*pi) - 1/2 * sum m_k * ln(L_k + m^2)
# For m^2 > 0, ALL eigenvalues (including k=0) contribute

framework_masses_sq = [
    ("m^2 = 0 (massless)", 0),
    ("m^2 = alpha", alpha),
    ("m^2 = alpha_s", alpha_s),
    ("m^2 = 1/(2*pi)", 1/(2*pi)),
    ("m^2 = 1/PHI^2", 1/PHI**2),
    ("m^2 = 1/degree", 1.0/degree),
    ("m^2 = alpha^2", alpha**2),
    ("m^2 = PHI^2 [Higgs-like]", PHI**2),
    ("m^2 = 1/N", 1.0/N),
]

print(f"  {'Mass':35s}  {'ln Z':>12s}  {'F/Vol':>12s}  {'ln|F/V|/ln(a)':>14s}")
for name, m2 in framework_masses_sq:
    if m2 == 0:
        # Massless: skip zero mode
        ld_m = sum(mults[k] * log(Lk[k] + m2) for k in range(N_eig) if Lk[k] > 1e-10)
        ln_Z_m = N_nonzero / 2 * ln_2pi - 0.5 * ld_m
    else:
        # Massive: all modes contribute
        ld_m = sum(mults[k] * log(Lk[k] + m2) for k in range(N_eig))
        ln_Z_m = N / 2 * ln_2pi - 0.5 * ld_m
    F_m = -ln_Z_m
    rho_m = F_m / Vol_S3
    ratio_m = log(abs(rho_m)) / ln_alpha if abs(rho_m) > 1e-300 else float('nan')
    print(f"  {name:35s}  {ln_Z_m:12.4f}  {rho_m:12.4f}  {ratio_m:14.4f}")

# =====================================================================
# TEST 6: Ratio of partition functions (relative to flat space)
# =====================================================================
print(f"\n--- TEST 6: Relative partition function ---")

# On flat lattice (complete graph K_N), L_flat has eigenvalue N with mult N-1
# Or: regular lattice with the same degree
# Simple comparison: Z_600 / Z_free where Z_free has uniform eigenvalues = degree

ld_flat = N_nonzero * log(degree)   # all eigenvalues = degree
ln_Z_flat = N_nonzero / 2 * ln_2pi - 0.5 * ld_flat

ln_ratio = ln_Z - ln_Z_flat  # = -1/2 * (ld - ld_flat)
print(f"  ln Z_600 = {ln_Z:.6f}")
print(f"  ln Z_flat(deg={degree}) = {ln_Z_flat:.6f}")
print(f"  ln(Z_600/Z_flat) = {ln_ratio:.6f}")
print(f"  (Z_600/Z_flat)/Vol = {ln_ratio / Vol_S3:.6f}")
print(f"  ln|ratio/Vol| / ln(alpha) = {log(abs(ln_ratio / Vol_S3)) / ln_alpha:.4f}"
      if abs(ln_ratio / Vol_S3) > 1e-300 else "  N/A")

# Compare to det'(L) / degree^{N-1}
print(f"\n  ln det'(L) - (N-1)*ln(degree) = {ld - ld_flat:.6f}")
print(f"  = ln[det'(L)/degree^(N-1)]")

# =====================================================================
# TEST 7: KK vacuum energy (effective 4D CC from internal space)
# =====================================================================
print(f"\n--- TEST 7: Kaluza-Klein vacuum energy ---")

# In KK, the 4D CC gets a one-loop contribution:
# Lambda_4D = (1/(2*V_internal)) * sum_k ln(L_k) = (1/(2*V)) * ln det'(L)
# Or more precisely: rho_4D = (1/V_int) * sum_k (1/2) * omega_k
# where omega_k = sqrt(L_k) are the KK masses

# One-loop vacuum energy: E_vac = (1/2) * sum sqrt(L_k)
E_vac_sum = 0.5 * sum(mults[k] * sqrt(Lk[k]) for k in range(N_eig) if Lk[k] > 1e-10)

# 4D CC = E_vac / V_internal
rho_KK = E_vac_sum / Vol_S3

print(f"  E_vac = (1/2)*sum m_k*sqrt(L_k) = {E_vac_sum:.6f}")
print(f"  Lambda_4D = E_vac / Vol(S^3) = {rho_KK:.6f}")
print(f"  ln(Lambda_4D) / ln(alpha) = {log(rho_KK) / ln_alpha:.6f}")

# Try with different normalizations of the KK masses
# Physical KK mass: m_n^2 = L_k / R^2 where R = internal radius
# If R ~ l_Planck, then m_KK ~ M_Planck * sqrt(L_k/degree)
for R_name, R_val in [("R=1", 1.0), ("R=1/sqrt(degree)", 1/sqrt(degree)),
                       ("R=1/sqrt(2*pi)", 1/sqrt(2*pi)), ("R=PHI/sqrt(N)", PHI/sqrt(N))]:
    E_vac_R = 0.5 * sum(mults[k] * sqrt(Lk[k]) / R_val for k in range(N_eig) if Lk[k] > 1e-10)
    V_R = Vol_S3 * R_val**3
    rho_R = E_vac_R / V_R if V_R > 0 else float('inf')
    ratio_R = log(rho_R) / ln_alpha if rho_R > 0 else float('nan')
    print(f"  {R_name:25s}: E_vac={E_vac_R:10.4f}, V={V_R:10.4f}, "
          f"rho={rho_R:10.4f}, ln(rho)/ln(a)={ratio_R:.4f}")

# =====================================================================
# TEST 8: Zeta-regularized one-loop effective action
# =====================================================================
print(f"\n--- TEST 8: Zeta-regularized effective action ---")

# W = -1/2 * zeta'_L(0) = 1/2 * ln det'(L) (the effective action)
# This is the proper QFT regularized quantity

W = 0.5 * ld  # = 1/2 * ln det'(L)
print(f"  W = (1/2) * ln det'(L) = {W:.10f}")
print(f"  W / Vol(S^3) = {W / Vol_S3:.10f}")
print(f"  W / (N * Vol) = {W / (N * Vol_S3):.10f}")
print(f"  ln|W/Vol| / ln(alpha) = {log(abs(W / Vol_S3)) / ln_alpha:.6f}")

# With 2*pi normalization: W = (1/2)*ln det'(L/(2*pi))
W_2pi = 0.5 * (ld - N_nonzero * ln_2pi)
print(f"\n  W_{'{2pi}'} = (1/2)*ln det'(L/(2*pi)) = {W_2pi:.10f}")
print(f"  W_2pi / Vol = {W_2pi / Vol_S3:.10f}")
print(f"  ln|W_2pi/Vol| / ln(alpha) = {log(abs(W_2pi / Vol_S3)) / ln_alpha:.6f}")

# =====================================================================
# TEST 9: exp(-F/Vol) - the Boltzmann weight itself
# =====================================================================
print(f"\n--- TEST 9: Boltzmann weight exp(-F/Vol) ---")

# Maybe the CC isn't F/Vol but exp(-F/Vol)?
# Lambda ~ exp(-F/Vol) = exp(ln Z / Vol)
rho_exp = exp(ln_Z / Vol_S3) if abs(ln_Z / Vol_S3) < 500 else float('inf')
print(f"  ln Z / Vol = {ln_Z / Vol_S3:.10f}")
print(f"  exp(ln Z / Vol) = {rho_exp:.10e}")
print(f"  ln(exp(lnZ/V)) / ln(alpha) = {(ln_Z / Vol_S3) / ln_alpha:.6f}")

# Try: Lambda = alpha^{F/(Vol * ln(1/alpha))}
# i.e., the CC exponent = F/Vol / |ln alpha|
z_eff = (-ln_Z / Vol_S3) / abs(ln_alpha)
print(f"\n  z_eff = F / (Vol * |ln alpha|) = {z_eff:.10f}")
print(f"  Compare: z_CC = {z_CC:.10f}")
print(f"  Diff: {abs(z_eff - z_CC):.6f}")

# Try various normalization combos
print(f"\n  Systematic search: z = c * F / (Vol * |ln(alpha)|)")
for c_name, c_val in [("1", 1), ("2", 2), ("1/2", 0.5), ("pi", pi),
                        ("2*pi", 2*pi), ("N/2pi^2", N/(2*pi**2)),
                        ("1/pi", 1/pi), ("4*pi", 4*pi)]:
    z_try = c_val * (-ln_Z / Vol_S3) / abs(ln_alpha)
    diff = abs(z_try - z_CC)
    mark = "***" if diff < 0.5 else ("*" if diff < 2 else "")
    print(f"  c={c_name:10s}: z = {z_try:12.4f} (diff from z_CC: {diff:8.4f}) {mark}")

# =====================================================================
# TEST 10: The COMPLETE effective potential
# =====================================================================
print(f"\n--- TEST 10: Complete effective potential ---")

# V_eff = F/Vol = [-ln Z] / Vol
# But in QFT: V_eff = (1/Vol) * [1/2 * sum_k ln(L_k/(2*pi))]
# The (2*pi) comes from the proper-time regulator

# V_eff = (1/(2*Vol)) * sum_{k>0} m_k * ln(L_k / (2*pi))
V_eff_1 = sum(mults[k] * log(Lk[k] / (2*pi)) for k in range(N_eig) if Lk[k] > 1e-10) / (2 * Vol_S3)
print(f"  V_eff = (1/(2*Vol)) * sum m_k * ln(L_k/(2pi)) = {V_eff_1:.10f}")
print(f"  ln|V_eff| / ln(alpha) = {log(abs(V_eff_1)) / ln_alpha:.6f}")

# V_eff with L/(2*pi*degree)
V_eff_2 = sum(mults[k] * log(Lk[k] / (2*pi*degree)) for k in range(N_eig) if Lk[k] > 1e-10) / (2 * Vol_S3)
print(f"  V_eff [L/(2pi*deg)] = {V_eff_2:.10f}")
if abs(V_eff_2) > 1e-300:
    print(f"  ln|V_eff| / ln(alpha) = {log(abs(V_eff_2)) / ln_alpha:.6f}")

# =====================================================================
# TEST 11: Functional determinant ratios (key idea)
# =====================================================================
print(f"\n--- TEST 11: det'(L) / (2*pi)^{N_nonzero} ---")

# The partition function IS:
# Z = (2*pi)^{(N-1)/2} * det'(L)^{-1/2}
# So Z^2 = (2*pi)^{N-1} / det'(L)
# Or: det'(L) / (2*pi)^{N-1} = 1/Z^2

# log[det'(L) / (2pi)^{N-1}] = ld - (N-1)*ln(2pi) = -2*ln Z
val_11 = ld - N_nonzero * ln_2pi
print(f"  ln[det'(L)/(2pi)^119] = {val_11:.10f}")
print(f"  = -2 * ln Z = {-2 * ln_Z:.10f}")
print(f"  val / ln(alpha) = {val_11 / ln_alpha:.6f}")
print(f"  val / (2*Vol) = {val_11 / (2 * Vol_S3):.10f}")
print(f"  ln|val/(2*Vol)| / ln(alpha) = {log(abs(val_11 / (2*Vol_S3))) / ln_alpha:.6f}")

# What about det'(L/(2pi))^{1/(N-1)} = geometric mean of L_k/(2pi)?
geom_mean = exp(sum(mults[k] * log(Lk[k]/(2*pi)) for k in range(N_eig) if Lk[k] > 1e-10) / N_nonzero)
print(f"\n  Geometric mean of L_k/(2pi): {geom_mean:.10f}")
print(f"  ln(geom_mean) / ln(alpha) = {log(geom_mean) / ln_alpha:.6f}")
print(f"  (geom_mean)^{N_nonzero} / alpha^z_CC = exp({N_nonzero * log(geom_mean) - z_CC * ln_alpha:.6f})")

# =====================================================================
# TEST 12: Compact scalar on S^1 (discrete sum)
# =====================================================================
print(f"\n--- TEST 12: Compact scalar on S^1_{'{R}'} ---")

# For a compact scalar phi ~ phi + 2*pi*R:
# Z_compact = Z_noncompact * theta_3(0 | tau) where tau encodes the radius
# But on a finite graph, the compact partition function is:
# Z = sum_{n in Z^N} exp(-1/2 * R^2 * (n + phi_0)^T L (n + phi_0))
# For phi_0 = 0: Z = sum_{n in Z^N} exp(-R^2/2 * n^T L n)
# This is essentially a theta function of the lattice L.

# At small R (decompactification), Z_compact -> Z_noncompact * sqrt(det'(L)) * R^{N-1}
# At large R (compact limit), Z -> 1 (ground state dominates)

# The RATIO: Z_compact / Z_noncompact probes the topology

# For a compact boson of radius R = 1/sqrt(2*pi) (self-dual in 2D):
R_sd = 1 / sqrt(2*pi)
# Winding modes get mass m_w^2 = (2*pi*R)^2 * L_k
# At self-dual radius, momentum = winding

# Sum: Z_compact(R) / Z_gaussian = prod_{k>0} [theta_3(0, exp(-pi*R^2*L_k))]^{m_k}
# At R = 1: theta_3(0, exp(-pi*L_k))
print(f"  Compact scalar: Z_compact / Z_Gaussian at R=1:")
from math import cosh
ln_theta_ratio = 0
for k in range(N_eig):
    if Lk[k] > 1e-10:
        # theta_3(0, q) = 1 + 2*sum_{n=1}^inf q^{n^2}  where q = exp(-pi*R^2*L_k)
        # For large L_k, theta_3 ~ 1 (corrections exponentially small)
        q = exp(-pi * Lk[k])
        # theta_3 = 1 + 2*q + 2*q^4 + 2*q^9 + ...
        theta3 = 1.0
        for n in range(1, 20):
            theta3 += 2 * exp(-pi * Lk[k] * n**2)
        ln_contrib = mults[k] * log(theta3)
        ln_theta_ratio += ln_contrib
        if ln_contrib > 1e-15:
            print(f"    k={k}: L_k={Lk[k]:.4f}, q={q:.4e}, theta_3={theta3:.10f}, "
                  f"m_k*ln(theta_3)={ln_contrib:.4e}")

print(f"  Total: ln(Z_compact/Z_Gauss) = {ln_theta_ratio:.10e}")
print(f"  This is exponentially small: compact corrections negligible.")

# =====================================================================
# TEST 13: The key question - what COMBINATION gives alpha^57?
# =====================================================================
print(f"\n{'='*72}")
print("TEST 13: COMPREHENSIVE SEARCH - What gives alpha^57?")
print("=" * 72)

# We have these building blocks (all computable from the 600-cell):
blocks = {
    "ld": ld,                                    # ln det'(L)
    "ln_2pi": ln_2pi,                            # ln(2*pi)
    "ln_N": log(N),                              # ln(120)
    "ln_alpha": ln_alpha,                         # ln(alpha)
    "ln_phi": log(PHI),                          # ln(phi)
    "Vol": Vol_S3,                                # 2*pi^2
    "N-1": float(N_nonzero),                     # 119
    "N": float(N),                                # 120
}

target_val = z_CC * ln_alpha   # this is ln(alpha^z_CC) = ln(Lambda_CC)

print(f"\n  Target: z_CC * ln(alpha) = {target_val:.10f}")
print(f"  = ln(alpha^57 * alpha^(-alpha_s)) = ln({alpha**z_CC:.6e})")

print(f"\n  Building blocks:")
for name, val in blocks.items():
    print(f"    {name:10s} = {val:.10f}")

# Key combinations
print(f"\n  Key combinations tested:")
combos = [
    ("ln Z", ln_Z),
    ("-ln Z", -ln_Z),
    ("ln Z / Vol", ln_Z / Vol_S3),
    ("-ln Z / Vol", -ln_Z / Vol_S3),
    ("ld/2", ld/2),
    ("-ld/2", -ld/2),
    ("(ld - (N-1)*ln(2pi))/2", val_11/2),
    ("ln Z / N", ln_Z / N),
    ("ld / (2*Vol)", ld / (2*Vol_S3)),
    ("(ld - (N-1)*ln(N))", ld - N_nonzero*log(N)),
    ("(ld - (N-1)*ln(N))/2", (ld - N_nonzero*log(N))/2),
    ("ld / (2*N)", ld / (2*N)),
    ("ln Z * Vol", ln_Z * Vol_S3),
]

for name, val in combos:
    ratio = val / ln_alpha if ln_alpha != 0 else float('nan')
    diff = abs(ratio - z_CC)
    mark = " ***" if diff < 0.5 else (" *" if diff < 2 else "")
    print(f"  {name:40s} = {val:12.4f}, /ln(a) = {ratio:10.4f} "
          f"(diff from z_CC: {diff:8.4f}){mark}")

# =====================================================================
# TEST 14: The exp421 near-miss REVISITED with (2*pi) factors
# =====================================================================
print(f"\n{'='*72}")
print("TEST 14: exp421 near-miss REVISITED")
print("=" * 72)

# exp421 found: ln det'(L/N) / ln(alpha) = 56.97 (near 57)
ld_over_N = ld - N_nonzero * log(N)
ratio_old = ld_over_N / ln_alpha
print(f"\n  exp421 result: ln det'(L/N) / ln(alpha) = {ratio_old:.10f}")
print(f"  Gap from 57: {57 - ratio_old:.10f}")
print(f"  Gap from z_CC: {z_CC - ratio_old:.10f}")

# NEW: what if the normalization is L/(2*pi*N/degree) or similar?
# The gap is 57 - 56.97 = 0.026
# Need: (N-1) * ln(c/N) / ln(alpha) = 57 where c absorbs the correction
# => ln(c/N) = 57 * ln(alpha) / 119 + correction
# Actually: ld/ln(alpha) = raw ratio. ld - 119*ln(c) = 57*ln(alpha)
# => ln(c) = (ld - 57*ln(alpha)) / 119

ln_c_exact = (ld - 57 * ln_alpha) / N_nonzero
c_exact = exp(ln_c_exact)
print(f"\n  For EXACT 57: need ln(c) = {ln_c_exact:.10f}, c = {c_exact:.10f}")
print(f"  Compare: N = {N}, ln(N) = {log(N):.10f}")
print(f"  c/N = {c_exact/N:.10f}")
print(f"  Correction factor: c/N - 1 = {c_exact/N - 1:.6e}")

# What framework quantity is c?
print(f"\n  c = {c_exact:.10f}")
print(f"  N = {N}")
print(f"  N*exp(-alpha) = {N*exp(-alpha):.10f}")
print(f"  N*(1-alpha) = {N*(1-alpha):.10f}")
print(f"  N*(1-alpha/2) = {N*(1-alpha/2):.10f}")
print(f"  N/(1+alpha) = {N/(1+alpha):.10f}")
print(f"  N*exp(-alpha^2/(2*pi)) = {N*exp(-alpha**2/(2*pi)):.10f}")
print(f"  N - 1/2 = {N - 0.5}")

# What about z_CC instead of 57?
ln_c_zCC = (ld - z_CC * ln_alpha) / N_nonzero
c_zCC = exp(ln_c_zCC)
print(f"\n  For EXACT z_CC: c = {c_zCC:.10f}")
print(f"  c / N = {c_zCC/N:.10f}")
print(f"  Correction c/N - 1 = {c_zCC/N - 1:.6e}")

# =====================================================================
# TEST 15: Functional integral with MEASURE factor
# =====================================================================
print(f"\n{'='*72}")
print("TEST 15: Measure factor in functional integral")
print("=" * 72)

# In lattice QFT, the partition function includes a measure factor:
# Z = prod_i (a^{d/2}) * int prod_i dphi_i exp(-S)
# where a = lattice spacing, d = dimension of the manifold
# On S^3 (d=3), with "lattice spacing" a ~ (Vol/N)^{1/3}:
a_lattice = (Vol_S3 / N)**(1.0/3)
d_manifold = 3

print(f"  Lattice spacing: a = (Vol/N)^(1/3) = {a_lattice:.10f}")
print(f"  Measure factor per site: a^(d/2) = {a_lattice**(d_manifold/2):.10f}")
print(f"  Total measure: N * (d/2) * ln(a) = {N * d_manifold/2 * log(a_lattice):.10f}")

ln_Z_lattice = ln_Z + N * d_manifold / 2 * log(a_lattice)
F_lattice = -ln_Z_lattice
rho_lattice = F_lattice / Vol_S3

print(f"  ln Z_lattice = {ln_Z_lattice:.10f}")
print(f"  F_lattice / Vol = {rho_lattice:.10f}")
if abs(rho_lattice) > 1e-300:
    print(f"  ln|F_lat/Vol| / ln(alpha) = {log(abs(rho_lattice)) / ln_alpha:.6f}")

# =====================================================================
# TEST 16: N-dependence / large-N structure
# =====================================================================
print(f"\n--- TEST 16: Structure of ln Z ---")

# Decompose: ln Z = a_coeff * N + b_coeff * ln(N) + c_coeff
# From ln Z = (N-1)/2 * ln(2pi) - 1/2 * ld:
# At large N: the extensive part is ~N, the sub-extensive ~ln N

# Break down ld contribution by eigenvalue class:
print(f"  Eigenvalue contributions to ln det'(L):")
total_ld = 0
for k in range(N_eig):
    if Lk[k] > 1e-10:
        contrib = mults[k] * log(Lk[k])
        total_ld += contrib
        frac = contrib / ld * 100
        print(f"    k={k}: m_k={mults[k]:3d}, L_k={Lk[k]:8.4f}, "
              f"m*ln(L)={contrib:10.4f} ({frac:5.1f}%)")

print(f"  Total: {total_ld:.6f} (check: {ld:.6f})")

# The Galois-conjugate pairs:
print(f"\n  Galois pairs:")
pairs = [(1,7), (2,8)]
for (i,j) in pairs:
    c_i = mults[i] * log(Lk[i])
    c_j = mults[j] * log(Lk[j])
    print(f"    ({i},{j}): {c_i:.6f} + {c_j:.6f} = {c_i+c_j:.6f}")
    norm = Lk[i] * Lk[j]  # Galois norm
    print(f"    N(L) = L_{i}*L_{j} = {norm:.6f}")

# Rational eigenvalue contribution
rational_ld = sum(mults[k] * log(Lk[k]) for k in [3,4,5,6])
irrational_ld = sum(mults[k] * log(Lk[k]) for k in [1,2,7,8])
print(f"\n  Rational eigenvalue contribution: {rational_ld:.6f}")
print(f"  Irrational eigenvalue contribution: {irrational_ld:.6f}")
print(f"  Ratio irrational/rational: {irrational_ld/rational_ld:.6f}")


# =====================================================================
# FINAL SUMMARY
# =====================================================================
print(f"\n{'='*72}")
print("FINAL SUMMARY")
print("=" * 72)

print(f"""
  QUESTION: Can the partition function Z of a scalar field on the 600-cell
  produce the cosmological constant Lambda ~ alpha^57?

  CORE RESULTS:
    ln Z = (N-1)/2 * ln(2pi) - (1/2)*ln det'(L)
         = {ln_Z:.6f}

    F = -ln Z = {-ln_Z:.6f}
    F / Vol(S^3) = {-ln_Z / Vol_S3:.6f}

    ln|F/Vol| / ln(alpha) = {log(abs(-ln_Z / Vol_S3)) / ln_alpha:.4f}

  COMPARISON TO CC TARGET:
    alpha^z_CC = alpha^{{{z_CC:.4f}}} = {alpha**z_CC:.6e}
    F/Vol = {-ln_Z/Vol_S3:.6f}  (order 1, NOT order 10^{{-122}})

  WHY THIS DOESN'T WORK:
    The partition function Z of a FREE SCALAR on a FINITE GRAPH with
    N={N} vertices gives ln Z ~ O(N) ~ O(100). Dividing by Vol ~ O(10)
    gives F/Vol ~ O(10). The CC target alpha^57 ~ 10^{{-122}} is
    unreachable by ANY polynomial combination of O(1)-O(100) quantities.

    The CC suppression alpha^57 requires 57 powers of a SMALL number.
    A one-loop computation on 120 vertices produces at most O(N) = O(120)
    factors, all of order 1-15 (the eigenvalue range). No combination
    of 120 factors of order 10 can give 10^{{-122}}.

  STRUCTURAL INSIGHT:
    For the CC, we need a NON-PERTURBATIVE mechanism that produces
    exp(-1/alpha^n) ~ exp(-137^n) type suppression. The partition function
    Z ~ det'^{{-1/2}} is inherently one-loop/perturbative. The det'
    near-miss (56.97 vs 57) from exp421 remains the closest approach,
    but it's log det', NOT exp(-1/alpha) * det'.

  THE REAL PROBLEM:
    exp421 showed ln det'(L/N) / ln(alpha) = {ratio_old:.4f} (near 57).
    The new (2*pi) factors shift this but don't bridge the gap.
    The full partition function ADDS information (the (2*pi)^59.5 factor)
    but this just shifts ln Z by a constant, it doesn't change the
    power-law structure.

  CONFIRMED NEGATIVE: The partition function of a free scalar on the
  600-cell does NOT derive the CC magnitude alpha^57.
  The near-miss from exp421 (56.97 vs 57) is UNCHANGED by proper
  QFT normalization: (2*pi) factors contribute to ln Z but not to
  the ratio ln det'/ln(alpha).
""")

# Final status
print(f"  STATUS: NEGATIVE")
print(f"  (exp421 near-miss stands; proper normalization does not help)")
print(f"  (the CC functional form alpha^z remains UNEXPLAINED)")
print(f"\n{'='*72}")

sys.exit(0)

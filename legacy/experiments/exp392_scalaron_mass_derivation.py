"""
exp392_scalaron_mass_derivation.py
===================================
GOAL: Can the scalaron mass M be FULLY derived from 600-cell Seeley-DeWitt
coefficients, without external input (without A_s from CMB)?

Key idea: On a FINITE discrete space (600-cell), the test function f in
Tr(f(D^2/Lambda^2)) is no longer arbitrary -- the spectrum is discrete
(finite eigenvalues). If f = step function at cutoff, f_k are FIXED.

Framework: a1=5, b1=6, phi=(1+sqrt(5))/2, N=120
Seeley-DeWitt: A0=11, A1=62, A2=233=F_13, c_k = 2*N*A_k

BLIND: derive M first, predict A_s, compare with Planck 2018 AFTER.

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from math import sqrt, pi, log, factorial, gcd
from fractions import Fraction

print("=" * 72)
print("exp392: SCALARON MASS FROM 600-CELL SPECTRAL COEFFICIENTS")
print("=" * 72)

# =====================================================================
# PART 0: FRAMEWORK CONSTANTS
# =====================================================================
print("\nPART 0: Framework Constants from a1 = 5")
print("-" * 50)

a1 = 5
b1 = a1 + 1
PHI = (1 + sqrt(a1)) / 2
N = factorial(a1)  # 120
degree = 2 * (a1 + 1)  # 12
N_eig = 9
rank_E8 = N_eig - 1  # 8
h_E8 = a1 * b1  # 30
N_gen = 3

# Couplings
sin2_tW = b1 / (a1**2 + 1)  # 6/26
alpha_s_fw = 1 / (2 * PHI**3)
A_eq = 2 * pi
B_eq = -4 * a1 * PHI**4
C_eq = 1.0
disc = B_eq**2 - 4 * A_eq * C_eq
alpha = (-B_eq - sqrt(disc)) / (2 * A_eq)

# Physical constants
m_e_MeV = 0.51099895
m_e_GeV = m_e_MeV * 1e-3
M_Pl_GeV = 1.22089e19          # full Planck mass
M_Pl_red_GeV = M_Pl_GeV / sqrt(8 * pi)  # reduced Planck mass

# Seeley-DeWitt coefficients
A0 = 2 * a1 + 1                # 11
A1 = 2 * a1**2 + 2 * a1 + 2   # 62
A2 = 6 * a1**2 + 15 * a1 + 8  # 233
c0 = 2 * N * A0                # 2640
c1 = 2 * N * A1                # 14880
c2 = 2 * N * A2                # 55920

print(f"a1={a1}, b1={b1}, phi={PHI:.10f}, N={N}")
print(f"alpha = {alpha:.10e}, 1/alpha = {1/alpha:.6f}")
print(f"A0={A0}, A1={A1}, A2={A2}")
print(f"c0={c0}, c1={c1}, c2={c2}")


# =====================================================================
# PART 1: THE TEST FUNCTION PROBLEM ON FINITE SPECTRA
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: Test Function on Finite Discrete Spectrum")
print("=" * 72)

print("""
In the NCG spectral action S = Tr(f(D^2/Lambda^2)), the heat kernel
expansion gives:

  S ~ f_0 Lambda^4 a_0 + f_2 Lambda^2 a_2 + f_4 a_4 + ...

where the MOMENTS f_k depend on the test function f:
  f_0 = int_0^inf f(u) u du     (second moment)
  f_2 = int_0^inf f(u) du       (zeroth moment)
  f_4 = f(0)                    (value at origin)

On a CONTINUOUS manifold, f is arbitrary => f_k are FREE parameters.
On a FINITE discrete space (600-cell), the spectrum is discrete:
  D has exactly dim(H) = c0 = 2640 eigenvalues.

QUESTION: Does the discreteness of the spectrum fix f_k?
""")

# The eigenvalues of the Cayley Laplacian on 600-cell
# L_k = degree * (1 - chi_k(C10)/d_k) for irrep k of 2I
dims_2I = [1, 2, 3, 4, 5, 6, 4, 2, 3]
chi_C10 = [1, PHI, PHI, 1, 0, -1, -1, 1-PHI, 1-PHI]
Lk = [degree * (1 - chi_C10[k]/dims_2I[k]) for k in range(N_eig)]
mults = [d**2 for d in dims_2I]

print("Cayley Laplacian eigenvalues on 0-forms (from 2I characters):")
for k in range(N_eig):
    print(f"  L_{k} = {Lk[k]:10.6f}  (mult d_k^2 = {mults[k]:2d},"
          f"  dim = {dims_2I[k]})")

# These are eigenvalues of Lap_0 = D^2 on 0-forms
# The FULL D = d+d* on all forms (0,1,2,3-forms) has c0=2640 eigenvalues
# By Hodge duality: spectrum on p-forms = spectrum on (3-p)-forms

print(f"\nTotal eigenvalues: {sum(mults)} on 0-forms, {c0} on all forms")

# Spectral radius (maximum eigenvalue of D^2 on 0-forms)
L_max_0 = max(Lk)
L_min_nonzero = min(l for l in Lk if l > 0)
print(f"Max eigenvalue (0-forms): L_max = {L_max_0:.6f}")
print(f"Min nonzero eigenvalue: L_min = {L_min_nonzero:.6f}")
print(f"Spectral gap: L_min = {L_min_nonzero:.6f}")

# The D^2 eigenvalues on 1-forms are harder to compute exactly
# But we know: Tr(D^2) = c1 = 14880, Tr(D^4) = c2 = 55920

# For the spectral action: natural cutoff Lambda = sqrt(L_max)
# or Lambda^2 = degree = 12 (spectral radius of Cayley graph on 0-forms)

# Actually from exp383: Lambda^2 * R^2 = 35 where R^2 = phi^2/2
# => Lambda^2 = 70/phi^2 = 70/(phi+1) = 70/phi^2
R_sq = PHI**2 / 2
Lambda_sq_383 = 35 / R_sq
print(f"\nFrom exp383: R^2 = phi^2/2 = {R_sq:.6f}")
print(f"Lambda^2 * R^2 = 35 => Lambda^2 = {Lambda_sq_383:.6f}")
print(f"  = 70/phi^2 = {70/PHI**2:.6f}")


# =====================================================================
# PART 2: STEP FUNCTION TEST FUNCTION
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: Step Function f(x) = theta(1-x)")
print("=" * 72)

print("""
For the step (characteristic) function f(x) = theta(1-x):
  f_4 = f(0) = 1
  f_2 = int_0^inf theta(1-u) du = int_0^1 du = 1
  f_0 = int_0^inf theta(1-u) u du = int_0^1 u du = 1/2

Key ratio: f_2/f_4 = 1/1 = 1
""")

f4_step = 1.0
f2_step = 1.0
f0_step = 0.5

# Also try Gaussian and exponential
print("Other natural test functions:")
print(f"  Step:       f_4={f4_step}, f_2={f2_step}, f_0={f0_step},"
      f"  f_2/f_4 = {f2_step/f4_step}")

# Gaussian f(x) = exp(-x)
f4_gauss = 1.0
f2_gauss = 1.0   # int_0^inf e^(-u) du = 1
f0_gauss = 1.0   # int_0^inf e^(-u) u du = 1
print(f"  Gaussian:   f_4={f4_gauss}, f_2={f2_gauss}, f_0={f0_gauss},"
      f"  f_2/f_4 = {f2_gauss/f4_gauss}")

# Smooth cutoff f(x) = (1-x)^2 theta(1-x)
f4_smooth = 1.0
f2_smooth = 1.0/3    # int_0^1 (1-u)^2 du = 1/3
f0_smooth = 1.0/12   # int_0^1 (1-u)^2 u du = 1/12
print(f"  Smooth:     f_4={f4_smooth}, f_2={f2_smooth:.4f},"
      f"  f_0={f0_smooth:.4f}, f_2/f_4 = {f2_smooth/f4_smooth:.4f}")


# =====================================================================
# PART 3: SCALARON MASS FROM SPECTRAL ACTION
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: Scalaron Mass Computation")
print("=" * 72)

print("""
From the Chamseddine-Connes spectral action on M^4 x F:

  S = Tr(f(D^2/Lambda^2)) contains:

  (1) Einstein-Hilbert:  1/(16*pi*G) = f_2*Lambda^2*c1/(48*pi^2)
  (2) R^2 term:          alpha_R2 = f_4*c2/(320*pi^2)

  The Starobinsky action: S = (M_Pl^2/2) int [R + R^2/(6*M^2)] sqrt(g)

  Identifying: M^2 = M_Pl^2 / (12*alpha_R2)

  where M_Pl^2 = 16*pi/(G) = f_2*Lambda^2*c1/(3*pi)  [full Planck]
""")

# M^2/M_Pl^2 = 1/(12*alpha_R2/M_Pl^2)
# = M_Pl^2 / (12 * f_4*c2/(320*pi^2))
# = [f_2*Lambda^2*c1/(3*pi)] * 320*pi^2 / (12*f_4*c2)
# = (f_2*Lambda^2/f_4) * c1 * 320*pi / (36*c2)
# = (f_2/f_4) * Lambda^2 * (c1/c2) * 80*pi/9

# But M^2/M_Pl^2 = 1/(12*alpha_R2) * (16*pi*G)
# = (16*pi*G)/(12*f_4*c2/(320*pi^2))
# = (16*pi)/(12*f_4*c2/(320*pi^2)) * G
# = (16*pi*320*pi^2)/(12*f_4*c2) * G

# More carefully:
# M_Pl^2 = f_2*Lambda^2*c1/(3*pi)
# alpha_R2 = f_4*c2/(320*pi^2)
# M^2 = M_Pl^2/(12*alpha_R2) = [f_2*Lambda^2*c1/(3*pi)] / [12*f_4*c2/(320*pi^2)]
# = f_2*Lambda^2*c1/(3*pi) * 320*pi^2/(12*f_4*c2)
# = (f_2/f_4) * Lambda^2 * (c1/c2) * 320*pi/(36)
# = (f_2/f_4) * Lambda^2 * (A1/A2) * 80*pi/9

# KEY: M^2/M_Pl^2 depends ONLY on f_4 and c2 (Lambda^2 cancels in ratio):
# M^2/M_Pl^2 = [f_2*Lambda^2*c1/(3*pi)] / [12*f_4*c2/(320*pi^2)] / M_Pl^2
# Wait, M^2 = M_Pl^2/(12*alpha_R2), so
# M^2/M_Pl^2 = 1/(12*alpha_R2) = 320*pi^2/(12*f_4*c2) = 80*pi^2/(3*f_4*c2)

# WAIT! This is WRONG. alpha_R2 has dimensions. Let me redo.
# In the spectral action, the R^2 coefficient in the Lagrangian is:
# L = ... + alpha_R2 * R^2 + ...
# where alpha_R2 = f_4 * c2_fiber / (320*pi^2)  [dimensionless if f_4 dimless]
# And the EH term: L = 1/(2*kappa^2) * R + ...
# where 1/(2*kappa^2) = f_2*Lambda^2 * c1_fiber / (96*pi^2) [mass^2]

# Starobinsky: L = (1/2*kappa^2)[R + R^2/(6*M^2)]
# => R^2 coefficient = 1/(2*kappa^2) * 1/(6*M^2) = alpha_R2
# => M^2 = 1/(2*kappa^2) / (6*alpha_R2)
# => M^2 = [f_2*Lambda^2*c1/(96*pi^2)] / [6*f_4*c2/(320*pi^2)]
# = (f_2*Lambda^2*c1) * 320 / (96 * 6 * f_4 * c2)
# = (f_2/f_4) * Lambda^2 * (c1/c2) * 320/576
# = (f_2/f_4) * Lambda^2 * (c1/c2) * 5/9

print("Corrected derivation:")
print("  1/(2*kappa^2) = f_2*Lambda^2*c1_F/(96*pi^2)")
print("  alpha_R2 = f_4*c2_F/(320*pi^2)")
print("  M^2 = 1/(2*kappa^2*6*alpha_R2)")
print("       = (f_2/f_4)*Lambda^2*(c1/c2)*(320/(96*6))")
print("       = (f_2/f_4)*Lambda^2*(c1/c2)*(5/9)")
print()

# Also: M_Pl^2 = 1/kappa^2 = 2*f_2*Lambda^2*c1_F/(96*pi^2)
#              = f_2*Lambda^2*c1_F/(48*pi^2)
# => M^2/M_Pl^2 = M^2 * kappa^2
#                = [1/(12*alpha_R2)] * kappa^2 / kappa^2 ... hmm
# M^2/M_Pl^2 = M^2/(1/kappa^2) = M^2*kappa^2
# = [1/(2*kappa^2)] / (6*alpha_R2) * kappa^2
# = 1/(12*alpha_R2)
# = 320*pi^2/(12*f_4*c2_F)

# YES! Lambda^2 CANCELS in the ratio M^2/M_Pl^2:
ratio_M_Mpl_sq = 320 * pi**2 / (12 * f4_step * c2)
print(f"M^2/M_Pl^2 = 320*pi^2/(12*f_4*c2_F)")
print(f"  = {320*pi**2:.4f} / (12 * {f4_step} * {c2})")
print(f"  = {ratio_M_Mpl_sq:.6e}")
print(f"  M/M_Pl = {sqrt(ratio_M_Mpl_sq):.6e}")

# But WHICH M_Pl? This is the FULL Planck mass: M_Pl = sqrt(hbar*c/G)
# In Connes convention: M_Pl^2 = 1/kappa^2, where kappa^2 = 8*pi*G
# So M_Pl here = 1/sqrt(8*pi*G) = M_Pl_reduced

# Let me be explicit:
# kappa^2 = 8*pi*G => 1/kappa^2 = M_Pl_red^2 = M_Pl^2/(8*pi)
# The formula 1/(2*kappa^2) = f_2*Lambda^2*c1_F/(96*pi^2)
# => M_Pl_red^2 = f_2*Lambda^2*c1_F/(48*pi^2)

# And M^2/M_Pl_red^2 = 1/(12*alpha_R2) = 320*pi^2/(12*f_4*c2_F)

print(f"\nUsing REDUCED Planck mass convention (M_Pl_red = M_Pl/sqrt(8*pi)):")
print(f"  M_Pl_red = {M_Pl_red_GeV:.4e} GeV")
M_ratio_sq_red = 320 * pi**2 / (12 * f4_step * c2)
M_pred_red = M_Pl_red_GeV * sqrt(M_ratio_sq_red)
print(f"  M^2/M_Pl_red^2 = {M_ratio_sq_red:.6e}")
print(f"  M/M_Pl_red = {sqrt(M_ratio_sq_red):.6e}")
print(f"  M = {M_pred_red:.4e} GeV")

# Using FULL Planck mass
M_ratio_sq_full = 320 * pi**2 / (12 * f4_step * c2)
M_pred_full = M_Pl_red_GeV * sqrt(M_ratio_sq_full)
print(f"\n  M (step function, f_4=1) = {M_pred_full:.4e} GeV")

# Compare with CMB value
As_planck = 2.1e-9
for Ne_val in [60, 62]:
    M_cmb = M_Pl_red_GeV * sqrt(24 * pi**2 * As_planck) / Ne_val
    print(f"  M (from CMB, N_e={Ne_val}) = {M_cmb:.4e} GeV")

print(f"\n  Ratio M_pred/M_cmb(60) = {M_pred_full/(M_Pl_red_GeV * sqrt(24*pi**2*As_planck)/60):.2f}")


# =====================================================================
# PART 4: CAN A2/A0 FIX f_2?
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: Can A2/A0 = 233/11 Fix f_2?")
print("=" * 72)

print(f"""
The ratio A2/A0 = {A2}/{A0} = {A2/A0:.6f}

In the spectral action, the scalaron mass involves:
  M^2 = (f_2/f_4) * Lambda^2 * (c1/c2) * 5/9

The ratio M^2/M_Pl^2 = 320*pi^2/(12*f_4*c2_F) involves ONLY f_4 and c2.
Lambda^2 CANCELS.

This means A2/A0 affects the CC/EH ratio and EH/R^2 ratio separately,
but NOT the ratio M/M_Pl which depends only on f_4 and c2.

The question simplifies to: can we determine f_4?

For any smooth, positive, decreasing test function with f(0)=1:
  f_4 = f(0) = 1

This is UNIVERSAL for normalized test functions.
The variation comes from different normalization conventions.
""")

# The KEY result: M^2/M_Pl^2 = 320*pi^2/(12*c2_F) for f_4=1
# = 320*pi^2 / (12 * 55920)
# = 3158.27 / 671040
# = 4.706e-3

val = 320 * pi**2 / (12 * c2)
print(f"M^2/M_Pl_red^2 = 320*pi^2/(12*c2) = {val:.6e}")
print(f"  = {320*pi**2:.2f}/{12*c2}")
print(f"  M/M_Pl_red = {sqrt(val):.6e}")
print(f"  M = {M_Pl_red_GeV * sqrt(val):.4e} GeV")


# =====================================================================
# PART 5: WHAT f_4 IS NEEDED FOR THE OBSERVED M?
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: Required f_4 for Observed M (Reverse Engineering)")
print("=" * 72)

# From CMB: M ~ 2.86e13 GeV (for N_e=60)
M_cmb_60 = M_Pl_red_GeV * sqrt(24 * pi**2 * As_planck) / 60
M_cmb_62 = M_Pl_red_GeV * sqrt(24 * pi**2 * As_planck) / 62

print(f"Required: M^2/M_Pl_red^2 = 320*pi^2/(12*f_4*c2)")
print(f"  => f_4 = 320*pi^2/(12*c2*(M/M_Pl_red)^2)")

for Ne_val, M_cmb in [(60, M_cmb_60), (62, M_cmb_62)]:
    M_ratio = M_cmb / M_Pl_red_GeV
    f4_needed = 320 * pi**2 / (12 * c2 * M_ratio**2)
    print(f"\nN_e = {Ne_val}: M = {M_cmb:.4e} GeV, M/M_Pl_red = {M_ratio:.4e}")
    print(f"  f_4 needed = {f4_needed:.4e}")
    print(f"  log10(f_4) = {log(f4_needed)/log(10):.2f}")
    print(f"  This means f_4 ~ {f4_needed:.0e}, NOT O(1)")

print("""
CRITICAL FINDING: To reproduce the observed M ~ 3e13 GeV from CMB,
we need f_4 ~ 10^7. This is NOT natural for a test function with f(0)=1.

The standard convention is f_4 = f(0) = 1 for a normalized test function.
A large f_4 ~ 10^7 would require either:
  (a) An unnormalized test function (f(0) >> 1)
  (b) A very different convention for the spectral action
  (c) Additional physics not captured by the heat kernel expansion
""")


# =====================================================================
# PART 6: DISCRETE SPECTRAL ACTION (EXACT, NO EXPANSION)
# =====================================================================
print("=" * 72)
print("PART 6: Exact Discrete Spectral Action")
print("=" * 72)

print("""
On the FINITE 600-cell, we can compute the spectral action EXACTLY
(no asymptotic expansion needed). The trace is a finite sum:

  S(Lambda) = sum_i f(lambda_i^2 / Lambda^2)

where {lambda_i} are eigenvalues of D = d + d*.
For the step function: S(Lambda) = N(Lambda) = counting function.
""")

# Eigenvalues of D^2 on 0-forms (Laplacian eigenvalues)
# with multiplicities
print("0-form spectrum of D^2 = Laplacian:")
eigenvalues_0 = []
for k in range(N_eig):
    for _ in range(mults[k]):
        eigenvalues_0.append(Lk[k])
eigenvalues_0.sort()

print(f"  {len(eigenvalues_0)} eigenvalues from 0 to {max(eigenvalues_0):.4f}")

# By Hodge duality on S^3: spec(Lap_p) = spec(Lap_{3-p})
# So 3-form spectrum = 0-form spectrum (another 120 eigenvalues)
# 1-form and 2-form spectra are equal (Hodge) but harder to compute

# For the EXACT spectral action on all forms:
# We know c0=2640, c1=14880, c2=55920
# And we know the 0-form spectrum exactly.

# Compute exact moments on 0-forms:
Tr_D2_0 = sum(eigenvalues_0)  # = Tr(Lap_0) = N*degree = 1440
Tr_D4_0 = sum(l**2 for l in eigenvalues_0)
print(f"\n0-forms: Tr(D^2) = {Tr_D2_0:.4f} (should be {N*degree})")
print(f"         Tr(D^4) = {Tr_D4_0:.4f}")

# Full space = 0+1+2+3 forms. By Hodge: 0=3, 1=2
Tr_D2_full = c1  # = 14880
Tr_D4_full = c2  # = 55920

Tr_D2_1 = (c1 - 2 * Tr_D2_0) / 2
Tr_D4_1 = (c2 - 2 * Tr_D4_0) / 2
print(f"\n1-forms: Tr(D^2) = {Tr_D2_1:.4f}")
print(f"         Tr(D^4) = {Tr_D4_1:.4f}")

# Average eigenvalue on each form degree
N_vert = 120
N_edges = 720
N_triangles = 1200
N_tetra = 600
print(f"\nAverage D^2 eigenvalue:")
print(f"  0-forms: {Tr_D2_0/N_vert:.4f} (= degree = {degree})")
print(f"  1-forms: {Tr_D2_1/N_edges:.4f}")
print(f"  2-forms: {Tr_D2_1/N_triangles:.4f} (Hodge = 1-forms)")
print(f"  3-forms: {Tr_D2_0/N_tetra:.4f} (Hodge = 0-forms)")

# For the step function with Lambda^2 = L_max (spectral radius):
# S = N(L_max) = all eigenvalues included = dim(H) = c0
L_max = max(Lk)
print(f"\nSpectral radius on 0-forms: {L_max:.6f}")

# For cutoff at Lambda^2 = spectral_radius, all 0-form eigenvalues
# are included. But different cutoffs give different S.

# Count eigenvalues below various cutoffs (0-forms only)
print(f"\nCounting function N(Lambda^2) on 0-forms:")
for L_cut in sorted(set(Lk)):
    count = sum(1 for l in eigenvalues_0 if l <= L_cut + 1e-10)
    print(f"  Lambda^2 = {L_cut:10.4f}: N = {count:4d}/{N_vert}")


# =====================================================================
# PART 7: PREDICT A_s AND COMPARE
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: Predict A_s from Step Function")
print("=" * 72)

# With f_4 = 1 (step function):
# M^2/M_Pl_red^2 = 320*pi^2/(12*c2) = 4.706e-3
# M/M_Pl_red = 0.0686

# A_s = N_e^2 * M^2 / (24*pi^2*M_Pl_red^2)
# = N_e^2/(24*pi^2) * 320*pi^2/(12*c2)
# = N_e^2 * 320/(24*12*c2)
# = N_e^2 * 10/(9*c2)
# = N_e^2 * 10/(9*55920)

print("CMB scalar amplitude: A_s = N_e^2 * M^2 / (24*pi^2*M_Pl_red^2)")
print()

for Ne_val in [60, 62]:
    As_pred = Ne_val**2 * 320 * pi**2 / (24 * pi**2 * 12 * c2)
    As_pred_simplified = Ne_val**2 * 10 / (9 * c2)

    print(f"N_e = {Ne_val}:")
    print(f"  A_s = N_e^2 * 10/(9*c2)")
    print(f"      = {Ne_val}^2 * 10/(9*{c2})")
    print(f"      = {Ne_val**2 * 10}/{9 * c2}")
    print(f"      = {As_pred_simplified:.6e}")
    print(f"  Observed: A_s = {As_planck:.1e}")
    print(f"  Ratio pred/obs = {As_pred_simplified/As_planck:.2e}")
    print()

print("RESULT: The step function (f_4=1) predicts A_s ~ 0.07,")
print("which is 3.4e7 times LARGER than observed A_s = 2.1e-9.")
print("The predicted M is too LARGE by a factor of ~5800.")
print("This is the well-known hierarchy problem of inflation.")


# =====================================================================
# PART 8: FRAMEWORK-SPECIFIC TEST FUNCTION?
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: Can the 600-Cell Fix the Test Function?")
print("=" * 72)

print("""
The spectral action S = Tr(f(D^2/Lambda^2)) on M^4 x F uses:
  - D_total = D_M x 1 + gamma_5 x D_F
  - The test function f acts on the COMBINED spectrum

On the discrete 600-cell F, D_F has finitely many eigenvalues.
But D_M (on the continuous M^4) has a CONTINUOUS spectrum.
The expansion f_k comes from the continuous M^4 part.

CRITICAL INSIGHT: The discreteness of F does NOT fix f_k.
The f_k moments are properties of f acting on the CONTINUOUS
spectrum of D_M, not on D_F. The fiber D_F contributes only
through the traces c_k = Tr_F(D_F^{2k}).

Therefore: f_0, f_2, f_4 remain FREE even with a discrete fiber.
The 600-cell fixes c0, c1, c2 but NOT f_k.
""")

# However, let's check: is there a NATURAL value of f_4 in the framework?
print("--- Searching for natural f_4 from framework ---")
print()

# If M^2/M_Pl^2 = 320*pi^2/(12*f_4*c2), and we want M from CMB:
# f_4_needed = 320*pi^2/(12*c2*(M/M_Pl)^2)

# Can f_4 be a framework number?
M_cmb_60_ratio = M_cmb_60 / M_Pl_red_GeV
f4_needed_60 = 320 * pi**2 / (12 * c2 * M_cmb_60_ratio**2)

print(f"f_4 needed (N_e=60) = {f4_needed_60:.4e}")
print()

# Check various framework combinations
candidates = {
    'c0': c0,
    'c1': c1,
    'c2': c2,
    'c0/10': c0/10,
    'c1/10': c1/10,
    'N^2': N**2,
    'N^2/2': N**2/2,
    'c0^2/(4*pi^2)': c0**2/(4*pi**2),
    'A2^2/(4*pi^2)': A2**2/(4*pi**2),
    'N*A2': N*A2,
    'c2/(4*pi^2)': c2/(4*pi**2),
    'c2/(2*pi)': c2/(2*pi),
    'c1*A1/(4*pi^2)': c1*A1/(4*pi**2),
}

for name, val in candidates.items():
    if val > 0:
        ratio = f4_needed_60 / val
        if 0.5 < ratio < 2.0:
            print(f"  *** CLOSE: f_4 ~ {name} = {val:.2f},"
                  f" ratio = {ratio:.4f}")
        elif 0.1 < ratio < 10:
            print(f"  Near: f_4 ~ {name} = {val:.2f},"
                  f" ratio = {ratio:.4f}")


# =====================================================================
# PART 9: N_e FROM SLOW-ROLL WITH EXACT COEFFICIENTS
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: N_e = 60 from Slow-Roll")
print("=" * 72)

print("""
The slow-roll predictions are UNIVERSAL for Starobinsky:
  n_s = 1 - 2/N_e - 9/(2*N_e^2)
  r = 12/N_e^2

These depend on N_e ONLY, not on M or the spectral coefficients.
The 600-cell coefficients affect only M (and hence A_s), not n_s or r.
""")

ns_planck = 0.9649
ns_err = 0.0042

for Ne in [50, 55, 60, 62, 65]:
    ns = 1 - 2.0/Ne - 9.0/(2*Ne**2)
    r = 12.0/Ne**2
    sigma = (ns - ns_planck) / ns_err
    print(f"  N_e={Ne:3d}: n_s={ns:.6f} ({sigma:+.2f} sigma),"
          f"  r={r:.6f}")

print(f"\nN_e = 60 = N/2 = |A5| gives n_s = {1-2/60-9/7200:.6f}"
      f" ({(1-2/60-9/7200-ns_planck)/ns_err:+.2f} sigma)")
print(f"This is CONSISTENT but not a derivation of N_e.")
print(f"N_e depends on reheating physics, not on A_k values.")


# =====================================================================
# PART 10: ALTERNATIVE - CONNES f_4 FROM FIBER GEOMETRY
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: Alternative Approaches to Fix f_4")
print("=" * 72)

print("""
Three ideas to potentially fix f_4:

(A) SPECTRAL ZETA REGULARIZATION:
    On the finite space F, the spectral zeta function is:
    zeta_F(s) = sum_k d_k^2 * L_k^{-s}  (excluding zero mode)
    At s=0: zeta_F(0) gives the 'regularized dimension'
    Could this provide a natural normalization?
""")

# Compute spectral zeta on 0-forms
print("Spectral zeta on 0-forms: zeta(s) = sum d_k^2 * L_k^(-s)")
for s_val in [0.5, 1.0, 1.5, 2.0]:
    zeta_s = sum(mults[k] * Lk[k]**(-s_val)
                 for k in range(N_eig) if Lk[k] > 1e-10)
    print(f"  zeta({s_val}) = {zeta_s:.6f}")

# zeta(0) = sum d_k^2 (excluding zero mode) = sum mults - 1 = 120 - 1 = 119
zeta_0 = sum(mults[k] for k in range(N_eig) if Lk[k] > 1e-10)
print(f"  zeta(0) = {zeta_0} (= N - 1)")

print(f"""
(B) HEAT KERNEL AT t=1:
    The natural 'time' on the 600-cell is t=1 (one step).
    K(1) = Tr(exp(-D^2)) = sum d_k^2 * exp(-L_k)
""")

K1 = sum(mults[k] * np.exp(-Lk[k]) for k in range(N_eig))
print(f"  K(1) = Tr(exp(-D^2)) = {K1:.6f}")
print(f"  K(1)/c0 = {K1/N_vert:.6f} (fraction of states at low energy)")

print(f"""
(C) f_4 FROM FIBER TRACE:
    In the NCG product M^4 x F, one could argue that
    f_4 should be 'normalized' by the fiber:
    f_4_eff = f_4 * Tr_F(1) = f_4 * c0_F
    But this doesn't help (makes M even larger).
""")

print(f"""
(D) f FROM SPECTRAL PROJECTION:
    Choose f = projection onto eigenvalues below some scale.
    Natural scale: Lambda_F = degree = {degree}
    Then f_4 = 1, f_2 ~ Lambda_F, and M gets modified.
    But f_2/f_4 = Lambda_F = degree only shifts Lambda,
    not the ratio M/M_Pl.
""")


# =====================================================================
# PART 11: THE REAL CONSTRAINT - DIOPHANTINE
# =====================================================================
print("=" * 72)
print("PART 11: Diophantine Identity and M")
print("=" * 72)

# From exp381: 2*A1^2 + 1 = 3*A0*A2
print(f"Diophantine: 2*{A1}^2 + 1 = 3*{A0}*{A2}")
print(f"  2*3844 + 1 = {2*A1**2+1}, 3*11*233 = {3*A0*A2}")
print(f"  Check: {2*A1**2+1 == 3*A0*A2}")
print()

# This constrains: c1^2/(c0*c2) = A1^2/(A0*A2) = 3844/2563
# = 3/2 - 1/(2*A0*A2) = 3/2 - 1/5126
ratio_c = A1**2 / (A0 * A2)
print(f"c1^2/(c0*c2) = A1^2/(A0*A2) = {ratio_c:.10f}")
print(f"  = 3/2 - 1/{2*A0*A2} = 1.5 - {1/(2*A0*A2):.6e}")
print(f"  Near-geometric progression of c_k (ratio ~ sqrt(3/2))")

# If c_k formed a perfect GP: c1^2 = c0*c2 (ratio = 1)
# Actual: ratio = 1.4996... very close to 3/2
# This means: the spectral coefficients are ALMOST but not quite geometric
# The deviation 1/5126 encodes the unique structure of a1=5

# Does this constrain M? Only through c2:
print(f"\nM^2/M_Pl^2 = const/c2 = const/{c2}")
print(f"  = const/(2*N*A2) = const/(2*{N}*{A2})")
print(f"  The Diophantine constrains A0*A2 ~ (2/3)*A1^2,")
print(f"  but since M depends on c2 ALONE (not c0 or c1),")
print(f"  the Diophantine gives no additional constraint on M.")


# =====================================================================
# PART 12: f_4 FROM THE SPECTRAL ACTION PRINCIPLE ITSELF
# =====================================================================
print("\n" + "=" * 72)
print("PART 12: Is f_4 = f(0) Really Free?")
print("=" * 72)

print("""
In the Chamseddine-Connes framework, the test function f must be:
  - Smooth, positive, decreasing
  - f(0) > 0 (for the a_4 coefficient to contribute)

The convention f(0) = 1 is the standard normalization.
But physically, f represents the SPECTRAL DENSITY of states
above the cutoff Lambda. Its normalization is a priori free.

IN PRINCIPLE: if the spectral geometry completely determines
the physics, then f should be determined too. This is an
OPEN PROBLEM in NCG (the 'test function problem').

Connes (2008): "The test function f is the only remaining
freedom in the spectral action. It should ultimately be
determined by the dynamics."

For the 600-cell framework:
  - f_4 = f(0) is NOT determined
  - Therefore M_scalaron is NOT fully determined
  - Therefore A_s is NOT predicted

This is the SAME conclusion as exp391, now confirmed from
the discrete spectrum perspective.
""")


# =====================================================================
# PART 13: WHAT CAN BE DERIVED (POSITIVE RESULTS)
# =====================================================================
print("=" * 72)
print("PART 13: Positive Results - What IS Determined")
print("=" * 72)

print(f"""
Even though M is not derived, several quantities ARE fixed:

1. RATIO n_s/r: Both depend only on N_e.
   For any N_e: r/n_s_deviation = 12/N_e^2 / (2/N_e) = 6/N_e
   At N_e=60: r/(1-n_s) = 6/60 = 0.1 (PREDICTION)

2. SPECTRAL INDEX (if N_e=60): n_s = {1-2/60-9/7200:.6f}

3. TENSOR-TO-SCALAR ratio (if N_e=60): r = {12/3600:.6f}

4. RUNNING: dn_s/d(ln k) = -2/N_e^2 = {-2/3600:.2e}

5. The R^2 SHAPE of the potential (not its amplitude):
   V(phi) propto (1 - exp(-sqrt(2/3)*phi/M_Pl))^2
   This is DERIVED from the spectral action (any NCG).

6. VACUUM STABILITY: lambda_H > 0 to M_Pl (from exp390).
   This is SPECIFIC to the 600-cell (framework top mass).

7. LEPTOGENESIS VIABILITY: T_reh >> M_R for any M > 10^10 GeV.
""")


# =====================================================================
# PART 14: HONEST FINAL ASSESSMENT
# =====================================================================
print("=" * 72)
print("PART 14: Honest Final Assessment")
print("=" * 72)

print(f"""
QUESTION: Can M_scalaron be derived without A_s from CMB?

ANSWER: **NO.**

DETAILED REASONING:

1. The test function moments f_k (especially f_4 = f(0)) are FREE
   parameters in the NCG spectral action, even with a discrete fiber.
   The discreteness of the 600-cell fixes the FIBER traces c_k = 2*N*A_k
   but NOT the moments f_k which come from the continuous M^4 part.

2. The ratio M^2/M_Pl_red^2 = 320*pi^2/(12*f_4*c2_F) depends on f_4.
   Lambda^2 cancels, but f_4 does NOT cancel.
   Compare with Higgs (exp379): m_H^2/m_W^2 where f_0 DOES cancel.

3. For the step function (f_4=1):
   M_pred = {M_pred_full:.2e} GeV (too large by factor ~5800)
   A_s_pred = {60**2 * 10/(9*c2):.2e} (too large by factor ~3.4e7)

4. The A_k ratios (A2/A0 = 233/11, Diophantine 2*A1^2+1=3*A0*A2)
   constrain the SHAPE of the spectral action expansion but NOT
   the absolute scale M. A2/A0 does not fix M because M depends
   on c2 ALONE (c0 and c1 cancel in M^2/M_Pl^2).

5. N_e = 60 = N/2 = |A5| is NATURAL but not dynamically derived.
   Different N_e give different (but all consistent) predictions.

STATUS CLASSIFICATION:
  M_scalaron = NOT DERIVED (requires A_s, i.e. f_4 from CMB)
  n_s, r = STRUCTURAL (universal Starobinsky, from any NCG)
  R^2 potential = STRUCTURAL (spectral action on any M^4 x F)
  N_e = 60 = PATTERN (natural, not rigorous)
  Vacuum stability = DERIVED (specific to 600-cell, exp390)
  Leptogenesis = CONSISTENT (T_reh >> M_R)

NEGATIVE RESULT:
  [1] f_4 does NOT cancel in M^2/M_Pl^2 (unlike Higgs where f_0 cancels)
  [2] Discrete spectrum of F does NOT fix f_k (they come from M^4)
  [3] No framework combination gives the required f_4 ~ {f4_needed_60:.0e}
  [4] The Diophantine identity provides no additional constraint on M

WHAT REMAINS FREE:
  M_scalaron (equivalently A_s) is the ONLY free parameter in the
  inflationary sector. Everything else (n_s, r, potential shape,
  vacuum stability) is determined by N_e and the framework.

  This parallels the overall framework structure: m_e is the
  single free parameter for particle physics, and A_s (via f_4)
  is the single free parameter for inflation.

This is a GENUINE limitation of the spectral action approach,
not specific to the 600-cell. It applies to ANY NCG model.
The test function problem remains OPEN in noncommutative geometry.
""")

print("=" * 72)
print("EXPERIMENT COMPLETE")
print("=" * 72)

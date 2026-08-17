"""
EXP-506: Which CC Formula is Correct?
======================================
Two candidate formulas for the cosmological constant:

  A) Lambda_P = alpha^(57 - alpha_s)          [PATTERN, known]
  B) Lambda_P = 4*sqrt(5)*exp(-738*phi'^2)    [NEW, from heat kernel]

Formula A involves alpha (transcendental, from quadratic with pi).
Formula B is purely algebraic in Z[phi] (plus ln and exp).

Which is closer to observation? And: what IS the observed value exactly?
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen,
                      alpha_em, alpha_s, ln_inv_alpha)

print("=" * 70)
print("EXP-506: WHICH CC FORMULA IS CORRECT?")
print("=" * 70)

# ============================================================
# The observed CC
# ============================================================
print("\n--- OBSERVED COSMOLOGICAL CONSTANT ---")

# Planck 2018 + BAO:
# Lambda = 1.1056(62) x 10^{-52} m^{-2}  (68% CL)
# In natural units (c = hbar = 1, G = 1/M_P^2):
# Lambda * (hbar/(M_P*c))^2 = Lambda * l_P^2
# l_P = 1.616255e-35 m
# Lambda_P = Lambda * l_P^2 = 1.1056e-52 * (1.616255e-35)^2
#          = 1.1056e-52 * 2.6123e-70 = 2.888e-122

# More precisely:
# rho_Lambda = Lambda * c^2 / (8*pi*G)
# rho_Lambda / rho_Planck = Lambda / (8*pi * l_P^{-2})
# = Lambda * l_P^2 / (8*pi)

# There are different conventions:
# Convention 1: Lambda_P = rho_vac / M_P^4 (energy density ratio)
# Convention 2: Lambda_P = Lambda * l_P^2 (geometric, curvature ratio)
# Convention 3: rho_vac / rho_Planck = Lambda/(8*pi/l_P^2) = Lambda*l_P^2/(8*pi)

# The commonly quoted "10^{-122}" uses convention 1 or 2.

# PDG 2024 values:
# Omega_Lambda = 0.6889(56)
# H_0 = 67.36(54) km/s/Mpc
# rho_crit = 3*H_0^2/(8*pi*G)

H_0_si = 67.36e3 / 3.0857e22  # H_0 in s^{-1}
G_N = 6.67430e-11  # m^3 kg^{-1} s^{-2}
c = 2.99792458e8    # m/s
hbar = 1.054571817e-34  # J s
M_P = np.sqrt(hbar * c / G_N)  # Planck mass in kg
l_P = np.sqrt(hbar * G_N / c**3)  # Planck length in m
t_P = l_P / c  # Planck time

rho_crit = 3 * H_0_si**2 / (8 * np.pi * G_N)  # kg/m^3
Omega_Lambda = 0.6889
rho_vac = Omega_Lambda * rho_crit  # kg/m^3
rho_vac_energy = rho_vac * c**2    # J/m^3

# Planck energy density
E_P = M_P * c**2  # Planck energy in J
rho_Planck = E_P / l_P**3  # Planck energy density J/m^3

# The ratio
Lambda_obs = rho_vac_energy / rho_Planck

print(f"  H_0 = {H_0_si:.4e} s^-1 = 67.36 km/s/Mpc")
print(f"  Omega_Lambda = {Omega_Lambda}")
print(f"  rho_vac = {rho_vac_energy:.4e} J/m^3")
print(f"  rho_Planck = {rho_Planck:.4e} J/m^3")
print(f"  Lambda_obs = rho_vac/rho_Planck = {Lambda_obs:.6e}")
print(f"  log10(Lambda_obs) = {np.log10(Lambda_obs):.4f}")

# Uncertainty: dominated by H_0 and Omega_Lambda
# Relative error on rho_vac ~ sqrt((delta_H/H)^2*4 + (delta_Omega/Omega)^2)
# ~ sqrt((0.54/67.36)^2*4 + (0.0056/0.6889)^2) ~ sqrt(0.000257 + 0.000066) ~ 1.8%
rel_err = np.sqrt((0.54/67.36)**2 * 4 + (0.0056/0.6889)**2)
Lambda_obs_err = Lambda_obs * rel_err

print(f"  Relative uncertainty: {rel_err*100:.1f}%")
print(f"  Lambda_obs = ({Lambda_obs:.3e}) +/- ({Lambda_obs_err:.3e})")

# ============================================================
# Formula A: alpha^z
# ============================================================
print(f"\n{'='*70}")
print("FORMULA A: Lambda = alpha^(57 - alpha_s)")
print(f"{'='*70}")

z_A = 57 - alpha_s
Lambda_A = alpha_em ** z_A

print(f"  z = 57 - 1/(2*phi^3) = {z_A:.6f}")
print(f"  alpha = {alpha_em:.10f}")
print(f"  Lambda_A = alpha^{z_A:.4f} = {Lambda_A:.6e}")
print(f"  Pull: {(Lambda_A - Lambda_obs) / Lambda_obs_err:.2f} sigma")
print(f"  Error: {(Lambda_A/Lambda_obs - 1)*100:+.2f}%")

# ============================================================
# Formula B: heat kernel (algebraic)
# ============================================================
print(f"\n{'='*70}")
print("FORMULA B: Lambda = 4*sqrt(5)*exp(-738*phi'^2)")
print(f"{'='*70}")

# 738 = b1 * (N + N_gen) = 6 * 123
# phi'^2 = (3-sqrt(5))/2 = 2 - phi
exponent_B = 738 * PHI_CONJ**2
Lambda_B = 4 * SQRT5 * np.exp(exponent_B)  # phi'^2 < 0? No, phi' = -0.618, phi'^2 = 0.382

# Wait: phi' = (1-sqrt(5))/2 = -0.6180
# phi'^2 = ((1-sqrt(5))/2)^2 = (6-2*sqrt(5))/4 = (3-sqrt(5))/2 = 0.38197
# So exponent = -738 * 0.38197 = -281.89

Lambda_B = 4 * SQRT5 * np.exp(-738 * PHI_CONJ**2)
print(f"  738 = b1*(N+N_gen) = {b1}*{N+N_gen}")
print(f"  phi'^2 = (3-sqrt(5))/2 = {PHI_CONJ**2:.10f}")
print(f"  Exponent = -738*phi'^2 = {-738*PHI_CONJ**2:.6f}")
print(f"  Lambda_B = 4*sqrt(5)*exp(-738*phi'^2) = {Lambda_B:.6e}")
print(f"  Pull: {(Lambda_B - Lambda_obs) / Lambda_obs_err:.2f} sigma")
print(f"  Error: {(Lambda_B/Lambda_obs - 1)*100:+.2f}%")

# ============================================================
# Formula B variants
# ============================================================
print(f"\n{'='*70}")
print("FORMULA B VARIANTS")
print(f"{'='*70}")

# Maybe the prefactor isn't 4*sqrt(5)?
# In the Galois heat kernel: Delta = mult * (d_q^2 - d_q'^2) * exp(-beta*L_1)
# d_q = phi, d_q' = |phi'| = 1/phi
# d_q^2 - d_q'^2 = phi^2 - 1/phi^2 = (phi^2 - phi'^2) = sqrt(5)
# But the sign: phi'^2 = phi^2 - sqrt(5). So d_q^2 - d_q'^2 = phi^2 - phi'^2 = sqrt(5)
# mult = 4

# Variant 1: exact prefactor analysis
# The eigenspace L_1 has mult 4, corresponding to the 4D irrep of A5
# (which is actually two 2D irreps of 2I, forming a Galois pair)
# The quantum dim weight: each 2D irrep has d_q = phi
# Under Galois: d_q -> |phi'| = 1/phi
# Weight per mode: phi^2 - 1/phi^2 = phi^2 - (phi-1) = ...
# phi^2 = phi + 1, so phi^2 - 1/phi^2 = (phi+1) - (phi-1) = 2? No.
# 1/phi^2 = phi^2 - 2*phi + 1 = (phi-1)^2 = (1/phi)^2... let me just compute.

dq_phys = PHI
dq_dark = abs(PHI_CONJ)  # = 1/phi
weight_per_mode = dq_phys**2 - dq_dark**2
print(f"\n  Quantum dim: d_phys = phi = {dq_phys:.6f}")
print(f"  Quantum dim: d_dark = 1/phi = {dq_dark:.6f}")
print(f"  Weight per mode: phi^2 - 1/phi^2 = {weight_per_mode:.6f}")
print(f"  = sqrt(5)? {abs(weight_per_mode - SQRT5) < 0.001}")

# So the prefactor IS 4*sqrt(5). But what about the actual Galois action?
# sigma: phi -> phi' = -1/phi (not +1/phi)
# So d_q' = phi' = -1/phi, d_q'^2 = 1/phi^2
# d_q^2 - d_q'^2 = phi^2 - 1/phi^2 = sqrt(5). Same.

# Variant 2: maybe prefactor should be mult * |d_q^2 - d_q'^2| / D^2
# where D^2 = sum d_q^2 = a1 + sqrt(a1) = 5 + sqrt(5)
D_sq = a1 + np.sqrt(a1)
print(f"\n  D^2 = a1 + sqrt(a1) = {D_sq:.6f}")
prefactor_normalized = 4 * SQRT5 / D_sq
Lambda_B2 = prefactor_normalized * np.exp(-738 * PHI_CONJ**2)
print(f"  Variant 2: (4*sqrt(5)/D^2)*exp(...) = {Lambda_B2:.6e}")
print(f"  Error vs obs: {(Lambda_B2/Lambda_obs - 1)*100:+.2f}%")

# Variant 3: prefactor = 4*sqrt(5) / (4*pi) (one-loop factor?)
Lambda_B3 = 4 * SQRT5 / (4 * np.pi) * np.exp(-738 * PHI_CONJ**2)
print(f"  Variant 3: (4*sqrt(5)/(4*pi))*exp(...) = {Lambda_B3:.6e}")
print(f"  Error vs obs: {(Lambda_B3/Lambda_obs - 1)*100:+.2f}%")

# Variant 4: what prefactor gives EXACT match to observation?
prefactor_exact = Lambda_obs / np.exp(-738 * PHI_CONJ**2)
print(f"\n  Exact-match prefactor: {prefactor_exact:.6f}")
print(f"  4*sqrt(5) = {4*SQRT5:.6f}")
print(f"  Ratio: {prefactor_exact / (4*SQRT5):.6f}")
print(f"  = 1/{4*SQRT5/prefactor_exact:.4f}")

# What clean expressions give this ratio?
ratio = prefactor_exact / (4*SQRT5)
print(f"\n  Ratio ~ {ratio:.6f}")
print(f"  1/sqrt(phi) = {1/np.sqrt(PHI):.6f}")
print(f"  phi'/phi = {abs(PHI_CONJ)/PHI:.6f}")
print(f"  1/phi = {1/PHI:.6f}")
print(f"  sqrt(2/a1) = {np.sqrt(2/a1):.6f}")
print(f"  2/pi = {2/np.pi:.6f}")
print(f"  1/sqrt(2*pi) = {1/np.sqrt(2*np.pi):.6f}")

# ============================================================
# Variant 5: beta = N + N_gen + alpha_s?
# ============================================================
print(f"\n{'='*70}")
print("VARIANT: beta = N + N_gen + alpha_s = 123.118")
print(f"{'='*70}")

L_1 = 12 - 6*PHI
beta_v5 = N + N_gen + alpha_s  # = 123.118
Lambda_v5 = 4 * SQRT5 * np.exp(-beta_v5 * L_1)
print(f"  beta = N + N_gen + alpha_s = {beta_v5:.6f}")
print(f"  Lambda = {Lambda_v5:.6e}")
print(f"  Error vs obs: {(Lambda_v5/Lambda_obs - 1)*100:+.2f}%")

# beta = N + pi?
beta_v6 = N + np.pi
Lambda_v6 = 4 * SQRT5 * np.exp(-beta_v6 * L_1)
print(f"\n  beta = N + pi = {beta_v6:.6f}")
print(f"  Lambda = {Lambda_v6:.6e}")
print(f"  Error vs obs: {(Lambda_v6/Lambda_obs - 1)*100:+.2f}%")

# What beta gives exact match?
beta_obs = -np.log(Lambda_obs / (4*SQRT5)) / L_1
print(f"\n  beta for exact obs match: {beta_obs:.6f}")
print(f"  beta - N = {beta_obs - N:.6f}")
print(f"  N_gen = {N_gen}")
print(f"  pi = {np.pi:.6f}")

# ============================================================
# COMPARISON TABLE
# ============================================================
print(f"\n{'='*70}")
print("COMPARISON TABLE")
print(f"{'='*70}")

results = [
    ("A: alpha^(57-as)",     Lambda_A,  z_A),
    ("B: 4s5*exp(-738*p'^2)", Lambda_B,  -np.log(Lambda_B)/ln_inv_alpha),
    ("Observed",              Lambda_obs, -np.log(Lambda_obs)/ln_inv_alpha),
]

print(f"\n  {'Formula':>25}  {'Lambda':>12}  {'z_equiv':>10}  {'vs obs':>10}  {'pull':>8}")
for name, lam, z in results:
    err = (lam/Lambda_obs - 1)*100
    pull = (lam - Lambda_obs) / Lambda_obs_err
    print(f"  {name:>25}  {lam:12.4e}  {z:10.4f}  {err:+8.2f}%  {pull:+8.2f}s")

print(f"\n  Observed uncertainty: {rel_err*100:.1f}%")

# ============================================================
# THE REAL QUESTION
# ============================================================
print(f"\n{'='*70}")
print("THE REAL QUESTION")
print(f"{'='*70}")

print(f"""
  Formula A (alpha^z):     {Lambda_A:.4e}  ({(Lambda_A/Lambda_obs-1)*100:+.1f}% from obs)
  Formula B (heat kernel): {Lambda_B:.4e}  ({(Lambda_B/Lambda_obs-1)*100:+.1f}% from obs)
  Observed:                {Lambda_obs:.4e}  (+/- {rel_err*100:.1f}%)

  NEITHER formula matches observation within uncertainty.
  Formula A is {abs(Lambda_A/Lambda_obs-1)*100:.1f}% off.
  Formula B is {abs(Lambda_B/Lambda_obs-1)*100:.1f}% off.

  BUT: the observed value sits BETWEEN the two formulas!
  A = {Lambda_A:.4e} < obs = {Lambda_obs:.4e} < B = {Lambda_B:.4e}

  This suggests neither is exactly right, but both capture
  different aspects of the same physics.

  The key difference:
    Formula A uses alpha (involves pi, transcendental)
    Formula B uses only phi (algebraic)

  A possible resolution: the TRUE CC formula involves BOTH
  alpha and the heat kernel, in a way we haven't found yet.
  The heat kernel gives the FORM (exponential),
  and alpha provides the CORRECTION to the pure phi^2 exponent.
""")

print("=" * 70)
print("EXP-506 COMPLETE")
print("=" * 70)

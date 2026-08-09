"""
EXP-506b: CC Conventions - What Exactly Are We Comparing?
==========================================================
The CC comparison depends critically on convention.
Let's be PRECISE about what number we're trying to match.
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen,
                      alpha_em, alpha_s, ln_inv_alpha)

print("=" * 70)
print("EXP-506b: CC CONVENTIONS")
print("=" * 70)

# Physical constants (CODATA 2018)
G_N = 6.67430e-11       # m^3 kg^-1 s^-2
c = 2.99792458e8         # m/s
hbar = 1.054571817e-34   # J s
k_B = 1.380649e-23       # J/K

# Planck units
l_P = np.sqrt(hbar * G_N / c**3)   # 1.616255e-35 m
t_P = l_P / c                       # 5.391247e-44 s
M_P = np.sqrt(hbar * c / G_N)       # 2.176434e-8 kg
E_P = M_P * c**2                     # 1.956086e9 J
T_P = E_P / k_B                     # Planck temperature
rho_P = E_P / l_P**3                # Planck energy density

# Reduced Planck mass: M_Pl = M_P / sqrt(8*pi)
M_Pl_reduced = M_P / np.sqrt(8*np.pi)

print(f"  l_P = {l_P:.6e} m")
print(f"  M_P = {M_P:.6e} kg")
print(f"  M_Pl (reduced) = {M_Pl_reduced:.6e} kg")
print(f"  rho_P = E_P/l_P^3 = {rho_P:.6e} J/m^3")

# Cosmological parameters (Planck 2018)
H_0 = 67.36e3 / 3.0857e22  # s^-1
Omega_Lambda = 0.6889

# Critical density
rho_crit = 3 * H_0**2 * c**2 / (8 * np.pi * G_N)  # J/m^3
rho_vac = Omega_Lambda * rho_crit

print(f"\n  H_0 = {H_0:.4e} s^-1")
print(f"  rho_crit = {rho_crit:.6e} J/m^3")
print(f"  rho_vac = Omega_Lambda * rho_crit = {rho_vac:.6e} J/m^3")

# ============================================================
# Different conventions for "Lambda in Planck units"
# ============================================================
print(f"\n{'='*70}")
print("CONVENTIONS")
print(f"{'='*70}")

# Convention 1: rho_vac / rho_Planck (energy density ratio)
# rho_Planck = M_P * c^2 / l_P^3 = c^7 / (hbar * G^2)
conv1 = rho_vac / rho_P
print(f"\n  Convention 1: rho_vac / rho_Planck")
print(f"  = {conv1:.6e}")
print(f"  log10 = {np.log10(conv1):.2f}")

# Convention 2: Lambda_geometric * l_P^2
# Einstein: G_mu_nu + Lambda * g_mu_nu = 8*pi*G * T_mu_nu
# Lambda = 8*pi*G * rho_vac / c^4  (in m^-2)
Lambda_geom = 8 * np.pi * G_N * rho_vac / c**4
conv2 = Lambda_geom * l_P**2
print(f"\n  Convention 2: Lambda * l_P^2")
print(f"  Lambda = {Lambda_geom:.6e} m^-2")
print(f"  Lambda * l_P^2 = {conv2:.6e}")
print(f"  log10 = {np.log10(conv2):.2f}")

# Convention 3: rho_vac / (M_Pl_reduced^4 * c^... )
# In "natural units" with reduced Planck mass:
# rho_vac in units of M_Pl^4 (reduced)
E_Pl_reduced = M_Pl_reduced * c**2
rho_Pl_reduced = E_Pl_reduced**4 / (hbar**3 * c**3)
conv3 = rho_vac / rho_Pl_reduced
print(f"\n  Convention 3: rho_vac / M_Pl_reduced^4")
print(f"  = {conv3:.6e}")
print(f"  log10 = {np.log10(conv3):.2f}")

# Convention 4: rho_vac / (M_P^2 * c^2 / (8*pi*l_P^2))
# = 8*pi * rho_vac / rho_Planck
conv4 = 8 * np.pi * conv1
print(f"\n  Convention 4: 8*pi * rho_vac/rho_Planck")
print(f"  = {conv4:.6e}")
print(f"  log10 = {np.log10(conv4):.2f}")

# Convention 5: 3 * Omega_Lambda * H_0^2 / (8*pi) in Planck units
# = 3 * Omega_Lambda * (H_0 * t_P)^2 / (8*pi)
conv5 = 3 * Omega_Lambda * (H_0 * t_P)**2 / (8 * np.pi)
print(f"\n  Convention 5: 3*Omega_L*(H_0*t_P)^2 / (8*pi)")
print(f"  = {conv5:.6e}")
print(f"  log10 = {np.log10(conv5):.2f}")

# Summary
print(f"\n{'='*70}")
print("SUMMARY OF CONVENTIONS")
print(f"{'='*70}")

conventions = {
    "rho_vac / rho_Planck": conv1,
    "Lambda * l_P^2": conv2,
    "rho_vac / M_Pl_red^4": conv3,
    "8*pi * rho_vac/rho_P": conv4,
}

for name, val in conventions.items():
    z_equiv = -np.log(val) / ln_inv_alpha
    print(f"  {name:>25}: {val:.4e}  (10^{np.log10(val):.1f})  z_equiv = {z_equiv:.2f}")

# ============================================================
# Which convention was used in the paper?
# ============================================================
print(f"\n{'='*70}")
print("WHICH CONVENTION DOES THE PAPER USE?")
print(f"{'='*70}")

# The paper says: Lambda_P = alpha^(57-alpha_s), error 0.33%
# alpha^56.88 = 2.835e-122
# This matches convention 4 or something giving ~10^-122.

# Let me check: alpha^56.88 = 2.835e-122
target_alpha = alpha_em ** (57 - alpha_s)
print(f"\n  alpha^(57-alpha_s) = {target_alpha:.6e}")
print(f"  log10 = {np.log10(target_alpha):.4f}")

# Which convention matches?
print(f"\n  Convention matches:")
for name, val in conventions.items():
    ratio = target_alpha / val
    print(f"    {name}: ratio = {ratio:.4f} ({abs(ratio-1)*100:.1f}% if ~1)")

# The paper claims 0.33% error. Let's find which convention gives that.
for name, val in conventions.items():
    err = abs(target_alpha/val - 1) * 100
    if err < 5:
        print(f"\n  MATCH: {name} with {err:.2f}% error")

# Hmm, none match. The 2.888e-122 number must come from a specific definition.
# Let me try: rho_Lambda / (M_P^4 / (8*pi)^2) or similar normalization

# Actually: the standard "CC problem" number is:
# rho_Lambda / rho_Planck ~ 10^{-123}
# But many papers quote 10^{-122} because they use reduced Planck mass
# or different factors of 8*pi.

# Let me find the EXACT observed value that gives 0.33% error with alpha^56.88:
# alpha^56.88 = obs * (1 + 0.0033)
# obs = alpha^56.88 / 1.0033 = 2.826e-122
# or obs = alpha^56.88 * 1.0033 = 2.844e-122

obs_target = target_alpha / 1.0033
print(f"\n  For 0.33% error: obs = {obs_target:.4e}")
print(f"  This is {obs_target/conv1:.2f} * rho_vac/rho_Planck")
print(f"  = {obs_target/conv2:.2f} * Lambda*l_P^2")

# It's about 25 * rho_vac/rho_Planck
print(f"\n  Ratio: {target_alpha/conv1:.2f}")
print(f"  = 8*pi = {8*np.pi:.2f}? {abs(target_alpha/conv1 - 8*np.pi) < 1}")
print(f"  = a1^2 = 25? {abs(target_alpha/conv1 - 25) < 1}")

# THE PAPER USES rho_vac/rho_Planck WITH 8*pi FACTOR:
conv_paper = conv2  # Lambda*l_P^2
print(f"\n  Paper convention: Lambda*l_P^2 = {conv2:.4e}")
print(f"  alpha^56.88 / (Lambda*l_P^2) = {target_alpha/conv2:.4f}")

# Or: the paper uses Lambda in units where the Planck mass INCLUDES the 8pi:
# i.e., G = 1 (not 1/(8*pi*M_P^2) but 1/M_P^2)
# In that case rho_Planck is different by factor 8*pi
conv_paper2 = conv1 * (8*np.pi)**2
print(f"  alt: rho_vac * (8pi)^2 / rho_P = {conv_paper2:.4e}")
print(f"  ratio: {target_alpha/conv_paper2:.4f}")

# Let's just compute what the paper MUST be using:
# If alpha^56.88 = X and error is 0.33%, then X/obs ~ 1.003 or 0.997
# obs = X * 0.997 or X * 1.003

for sign, label in [(0.9967, "-0.33%"), (1.0033, "+0.33%")]:
    obs_calc = target_alpha / sign
    ratio_conv1 = obs_calc / conv1
    print(f"\n  If error = {label}: obs = {obs_calc:.4e}")
    print(f"    obs / conv1 = {ratio_conv1:.2f}")
    print(f"    obs / conv2 = {obs_calc/conv2:.2f}")

# ============================================================
# FINAL: Compare both formulas in SAME convention
# ============================================================
print(f"\n{'='*70}")
print("BOTH FORMULAS IN SAME CONVENTION")
print(f"{'='*70}")

# Use convention 2: Lambda * l_P^2 (geometric)
obs = conv2
print(f"\n  Using Lambda * l_P^2 = {obs:.4e}")

z_obs = -np.log(obs) / ln_inv_alpha
z_A = 57 - alpha_s
z_B = -np.log(4*SQRT5*np.exp(-738*PHI_CONJ**2)) / ln_inv_alpha

Lambda_A = alpha_em ** z_A
Lambda_B = 4 * SQRT5 * np.exp(-738 * PHI_CONJ**2)

print(f"  z_obs = {z_obs:.4f}")
print(f"  z_A = {z_A:.4f}, Lambda_A = {Lambda_A:.4e}, error = {(Lambda_A/obs-1)*100:+.1f}%")
print(f"  z_B = {z_B:.4f}, Lambda_B = {Lambda_B:.4e}, error = {(Lambda_B/obs-1)*100:+.1f}%")

# Convention 1: rho_vac / rho_Planck
obs1 = conv1
z_obs1 = -np.log(obs1) / ln_inv_alpha
print(f"\n  Using rho_vac/rho_Planck = {obs1:.4e}")
print(f"  z_obs = {z_obs1:.4f}")
print(f"  z_A = {z_A:.4f}, error = {(Lambda_A/obs1-1)*100:+.1f}%")
print(f"  z_B = {z_B:.4f}, error = {(Lambda_B/obs1-1)*100:+.1f}%")

print(f"\n  KEY: The '0.33%' claim in the paper requires knowing")
print(f"  WHICH exact convention they used for Lambda_P.")
print(f"  The factor between conventions is 8*pi ~ 25.")
print(f"  Both formulas A and B give 10^-122, not 10^-123.")
print(f"  The observed value is 10^-122.9 or 10^-123.0")
print(f"  depending on convention.")

print("\n" + "=" * 70)
print("EXP-506b COMPLETE")
print("=" * 70)

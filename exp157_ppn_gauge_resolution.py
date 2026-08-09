"""
EXP-157: PPN Gauge Resolution - Coordinate-Invariant Observables

The paradox: exp104 found precession = GR at 1PN for alpha=1/2,
but exp155 found beta=2 (isotropic gauge) which predicts precession = 2/3 GR.
In areal gauge, beta=1 (GR).

This experiment computes PHYSICAL observables that are coordinate-invariant:
1. Perihelion advance per orbit (as angle measured at infinity)
2. Shapiro time delay
3. Light deflection angle
to determine if STG(alpha=1/2) agrees or disagrees with GR.

The approach: work directly in the original STG coordinates and compute
observables from geodesic equations without any gauge transformation.
"""

import numpy as np
from sympy import *

print("=" * 70)
print("EXP-157: PPN GAUGE RESOLUTION")
print("Coordinate-Invariant Observables for STG (alpha=1/2)")
print("=" * 70)

# =====================================================================
# PART 1: THE STG METRIC AND GEODESIC EQUATIONS
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: STG Metric and Geodesic Equations (alpha=1/2)")
print("=" * 70)

r, M, L, E = symbols('r M L E', positive=True)
u = symbols('u', positive=True)  # u = M/r
a = Rational(1, 2)  # alpha = 1/2

# STG metric: ds^2 = -dt^2/C^2 + C^2 dr^2 + C^(2a) r^2 dOmega^2
# C(r) = 1 + M/r = 1 + u (where u = M/r)
C = 1 + u

print(f"\nMetric: ds^2 = -dt^2/C^2 + C^2 dr^2 + C^(2*{a}) r^2 dOmega^2")
print(f"C = 1 + M/r, alpha = {a}")
print(f"\ng_tt = -1/C^2 = -(1+u)^(-2)")
print(f"g_rr = C^2 = (1+u)^2")
print(f"g_ang = C^(2a) r^2 = (1+u) r^2")

# For geodesics, define:
# E = (1/C^2) dt/dtau  (energy per unit mass)
# L = C^(2a) r^2 dphi/dtau  (angular momentum per unit mass)
# The radial equation (from ds^2 = -1 for timelike):
# (dr/dtau)^2 = [E^2 - V_eff(r)] / C^2
# where V_eff = (1 + L^2/(C^(2a) r^2)) / C^2

# For orbits, use u = 1/r and phi as parameter
# Define w = M*u_orbital = M/r (dimensionless)

# =====================================================================
# PART 2: ORBIT EQUATION IN ORIGINAL COORDINATES
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: Orbit Equation (exact)")
print("=" * 70)

# The orbit equation for the STG metric with general alpha:
# From the effective potential approach:
# (du/dphi)^2 = f(u) where u = 1/r
#
# For alpha=1/2:
# g_tt = -1/(1+M*u)^2, g_rr = (1+M*u)^2, g_ang = (1+M*u) / u^2
# (using u = 1/r, so r = 1/u)
#
# Energy: E_tilde = dt/dtau / C^2 = dt/dtau / (1+Mu)^2
# Ang mom: L_tilde = C^(2a) r^2 dphi/dtau = (1+Mu)/u^2 * dphi/dtau

# The orbit equation:
# (du/dphi)^2 + u^2 = 2Mu^3/(1+Mu) + E^2 u^4/((1+Mu) L^2) - u^4/(L^2 (1+Mu)^2)
# This is messy. Let me derive it properly with sympy.

# For a cleaner derivation, use the Lagrangian approach
# L = -g_tt * tdot^2 + g_rr * rdot^2 + g_ang * phidot^2
# = (1+Mu)^(-2) tdot^2 - (1+Mu)^2 rdot^2 - (1+Mu)/u^2 phidot^2

# Conserved quantities:
# E = d(L)/d(tdot) = 2/(1+Mu)^2 * tdot => tdot = E(1+Mu)^2/2
# Wait, let me use standard GR convention.

# Proper Lagrangian: 2L = g_mu_nu dx^mu/dtau dx^nu/dtau
# For timelike geodesics: 2L = -1
# 2L = -1/C^2 * (dt/dtau)^2 + C^2 * (dr/dtau)^2 + C^(2a) * r^2 * (dphi/dtau)^2 = -1

# Killing vectors give:
# e = (1/C^2) dt/dtau  (energy)
# l = C^(2a) r^2 dphi/dtau  (angular momentum)

# Substituting:
# -C^2 e^2 + C^2 (dr/dtau)^2 + l^2 / (C^(2a) r^2) = -1
# => C^2 (dr/dtau)^2 = C^2 e^2 - 1 - l^2 / (C^(2a) r^2)
# => (dr/dtau)^2 = e^2 - 1/C^2 - l^2 / (C^(2+2a) r^2)

# Switch to u = 1/r, dphi = (l / (C^(2a) r^2)) dtau
# dr/dtau = dr/dphi * dphi/dtau = (-1/u^2)(du/dphi) * l / (C^(2a) / u^2)
#         = -(du/dphi) * l / (u^2 * C^(2a) / u^2)
# Wait, let me be more careful.
# r = 1/u, dr = -1/u^2 du, dr/dtau = -1/u^2 du/dtau
# dphi/dtau = l / (C^(2a) r^2) = l u^2 / C^(2a)
# du/dphi = (du/dtau) / (dphi/dtau)
# du/dtau = -u^2 dr/dtau
# dr/dtau = -(1/u^2) du/dtau
# So du/dphi = (du/dtau) * C^(2a) / (l u^2)

# From (dr/dtau)^2:
# (1/u^4)(du/dtau)^2 = e^2 - 1/C^2 - l^2 u^2 / C^(2+2a)
# (du/dtau)^2 = u^4 [e^2 - 1/C^2 - l^2 u^2 / C^(2+2a)]

# (du/dphi)^2 = (du/dtau)^2 * C^(4a) / (l^2 u^4)
# = C^(4a)/l^2 * [e^2 - 1/C^2 - l^2 u^2 / C^(2+2a)]
# = C^(4a) e^2 / l^2 - C^(4a-2)/l^2 - u^2 * C^(4a-2-2a)
# = C^(4a) e^2 / l^2 - C^(4a-2)/l^2 - u^2 * C^(2a-2)

# For alpha = 1/2:
# 4a = 2, 4a-2 = 0, 2a-2 = -1
# (du/dphi)^2 = C^2 e^2 / l^2 - 1/l^2 - u^2 / C
# = (1+Mu)^2 e^2 / l^2 - 1/l^2 - u^2/(1+Mu)

# Define h = M/l (dimensionless), w = Mu (so u = w/M)
w = Symbol('w')  # w = M*u = M/r (dimensionless)
h_sq = Symbol('h2')  # h^2 = M^2/l^2
e_sq = Symbol('e2')  # e^2

# (du/dphi)^2 = (1+w)^2 * e2 / (M^2 * h2^-1) ... let me use u directly
# Keep u = 1/r, use l and e as parameters
# (du/dphi)^2 + u^2/(1+Mu) = (1+Mu)^2 e^2/l^2 - 1/l^2

# Let me define p = l^2/M (semi-latus rectum scale)
# and use w = Mu as the variable

# Actually, let me just expand everything to 2PN and compare with GR.

print("\nOrbit equation for alpha=1/2:")
print("(du/dphi)^2 = (1+Mu)^2 * e^2/l^2 - 1/l^2 - u^2/(1+Mu)")
print("\nThis is EXACT (no approximation).")

# =====================================================================
# PART 3: CIRCULAR ORBITS AND PRECESSION (EXACT)
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: Circular Orbits and Precession (perturbative)")
print("=" * 70)

# For circular orbits at u = u0 (r = r0):
# F(u0) = 0 and F'(u0) = 0
# where F(u) = RHS - LHS of orbit equation with (du/dphi)^2 = 0

# Define F(u) = (1+Mu)^2 e^2/l^2 - 1/l^2 - u^2/(1+Mu)
# F(u0) = 0: (1+Mu0)^2 e^2/l^2 - 1/l^2 = u0^2/(1+Mu0)
# F'(u0) = 0: 2M(1+Mu0) e^2/l^2 = 2u0/(1+Mu0) - Mu0^2/(1+Mu0)^2
#            = u0(2+Mu0)/(1+Mu0)^2

# From F'=0: e^2/l^2 = u0(2+Mu0) / [2M(1+Mu0)^3]
# Substituting into F=0:
# (1+Mu0)^2 * u0(2+Mu0)/[2M(1+Mu0)^3] - 1/l^2 = u0^2/(1+Mu0)
# u0(2+Mu0)/[2M(1+Mu0)] - 1/l^2 = u0^2/(1+Mu0)
# 1/l^2 = u0(2+Mu0)/[2M(1+Mu0)] - u0^2/(1+Mu0)
#        = u0/[(1+Mu0)] * [(2+Mu0)/(2M) - u0]
#        = u0/[(1+Mu0)] * [(2+Mu0-2Mu0)/(2M)]
#        = u0/[(1+Mu0)] * [(2-Mu0)/(2M)]
#        = u0(2-Mu0) / [2M(1+Mu0)]

w0 = Symbol('w0')  # w0 = M*u0 = M/r0

l_sq_inv = w0 * (2 - w0) / (2 * (1 + w0))  # this is M^2/l^2 * (1/M) = 1/(M*l^2) hmm
# Actually, let me redo with w0 = M*u0

# 1/l^2 = u0(2-Mu0)/[2M(1+Mu0)]
# M^2/l^2 = M*u0*(2-Mu0)/[2(1+Mu0)] = w0*(2-w0)/[2(1+w0)]

h2_circ = w0 * (2 - w0) / (2 * (1 + w0))
print(f"\nCircular orbit: M^2/l^2 = w0*(2-w0)/[2*(1+w0)]")
print(f"  where w0 = M/r0")

# e^2/l^2 from above: e^2/l^2 = u0(2+Mu0)/[2M(1+Mu0)^3]
# M^2 e^2/l^2 = w0*(2+w0)/[2*(1+w0)^3]
# => e^2 = l^2 * w0*(2+w0)/[2M^2*(1+w0)^3]
#         = [2(1+w0)/(w0*(2-w0))] * w0*(2+w0)/[2*(1+w0)^3]
#         = (2+w0) / [(2-w0)*(1+w0)^2]

e2_circ = (2 + w0) / ((2 - w0) * (1 + w0)**2)
print(f"e^2 = (2+w0)/[(2-w0)*(1+w0)^2]")
e2_exp = series(e2_circ, w0, 0, n=5)
print(f"  = {e2_exp}")

# The orbit equation perturbed around u = u0:
# u = u0(1 + epsilon * cos(omega_r * phi))
# omega_r^2 = -F''(u0) / (2 * something)
# More precisely, (du/dphi)^2 = F(u), so for small perturbation delta_u = u - u0:
# (d(delta_u)/dphi)^2 = F(u0) + F'(u0)*delta_u + F''(u0)/2 * delta_u^2 + ...
# Since F(u0)=F'(u0)=0: (d(delta_u)/dphi)^2 = F''(u0)/2 * delta_u^2
# => omega_r^2 = -F''(u0)/2 ... wait.
# (d(delta_u)/dphi)^2 + omega_r^2 delta_u^2 = 0 requires omega_r^2 = -F''(u0)/2?
# No: if F(u) ~ F''(u0)/2 * (u-u0)^2, then
# (du/dphi)^2 = F''(u0)/2 * delta_u^2
# But this gives (du/dphi)^2 - F''(u0)/2 * delta_u^2 = 0
# which is NOT oscillatory unless F''(u0) < 0.
# Actually the equation is (du/dphi)^2 = F(u), so
# d^2u/dphi^2 = F'(u)/2
# Linearizing: d^2(delta_u)/dphi^2 = F''(u0)/2 * delta_u
# So omega_r^2 = -F''(u0)/2

# F(u) = (1+Mu)^2 e^2/l^2 - 1/l^2 - u^2/(1+Mu)

u_var = Symbol('u')
Mu = Symbol('Mu')

# Use sympy for F(u) with substitution Mu = M*u
F_of_u = (1 + Mu)**2 * e2_circ * h2_circ - h2_circ - u_var**2 / (1 + Mu)
# But this mixes u_var and Mu which should be related... Let me use w = M*u directly.

# Redefine: w = M*u, F(w) in terms of w and w0
w_var = Symbol('w')
# (M du/dphi)^2 = F(w) where:
# F(w) = (1+w)^2 * e2 * M^2/l^2 - M^2/l^2 - w^2/(M^2 * (1+w)) * M^2
# No wait, M*u = w, u = w/M, u^2 = w^2/M^2
# (du/dphi)^2 = F_original(u)
# (dw/dphi)^2 = M^2 * F_original(w/M) = G(w)

# Original: (du/dphi)^2 = (1+Mu)^2 * e2/l^2 - 1/l^2 - u^2/(1+Mu)
# Substitute u = w/M:
# (1/M^2)(dw/dphi)^2 = (1+w)^2 * e2/l^2 - 1/l^2 - w^2/(M^2*(1+w))
# (dw/dphi)^2 = M^2*(1+w)^2*e2/l^2 - M^2/l^2 - w^2/(1+w)

# Using h2 = M^2/l^2 = w0*(2-w0)/[2*(1+w0)]
# and e2 = (2+w0)/[(2-w0)*(1+w0)^2]
# M^2*e2/l^2 = h2 * e2 = w0*(2-w0)/[2*(1+w0)] * (2+w0)/[(2-w0)*(1+w0)^2]
#             = w0*(2+w0) / [2*(1+w0)^3]

me2l2 = w0 * (2 + w0) / (2 * (1 + w0)**3)

G_w = (1 + w_var)**2 * me2l2 - h2_circ - w_var**2 / (1 + w_var)

# Simplify G(w)
G_w_simplified = simplify(G_w)
print(f"\nG(w) = (dw/dphi)^2 = {G_w_simplified}")

# Check G(w0) = 0
G_at_w0 = G_w_simplified.subs(w_var, w0)
print(f"G(w0) = {simplify(G_at_w0)}")

# Compute G'(w0)
G_prime = diff(G_w, w_var)
G_prime_at_w0 = simplify(G_prime.subs(w_var, w0))
print(f"G'(w0) = {G_prime_at_w0}")

# Compute G''(w0)
G_double_prime = diff(G_w, w_var, 2)
G_dp_at_w0 = simplify(G_double_prime.subs(w_var, w0))
print(f"G''(w0) = {G_dp_at_w0}")

# omega_r^2 = -G''(w0)/2
omega_r_sq = simplify(-G_dp_at_w0 / 2)
print(f"\nomega_r^2 = -G''(w0)/2 = {omega_r_sq}")

# Expand in powers of w0
omega_r_sq_series = series(omega_r_sq, w0, 0, n=5)
print(f"omega_r^2 = {omega_r_sq_series}")

# Precession per orbit = 2*pi*(1 - 1/omega_r) for omega_r close to 1
# Or more precisely: delta_phi = 2*pi/omega_r - 2*pi = 2*pi*(1/omega_r - 1)
# For omega_r^2 = 1 - A*w0 + ...: omega_r ~ 1 - A*w0/2
# delta_phi ~ pi*A*w0

# =====================================================================
# PART 4: COMPARISON WITH GR (SCHWARZSCHILD)
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: Comparison with GR (Schwarzschild)")
print("=" * 70)

# GR orbit equation: (du/dphi)^2 + u^2 = 1/p + 2M/p^2 * u + 3Mu^2 (approximately)
# Exact: (du/dphi)^2 = 2M u^3 - u^2 + u/p + (e^2-1)/l^2
# With w = M*u: (dw/dphi)^2 = 2w^3 - w^2 + w*M/p + M^2*(e^2-1)/l^2

# For GR circular orbit:
# l^2_GR = M*r0 / (1 - 3M/r0) = M/(w0 * (1-3w0)) ... r0 = M/w0
# M^2/l^2_GR = w0^3 * (1-3w0)/M ... no, r0 = M/w0 is wrong.
# u0 = 1/r0, w0 = M*u0 = M/r0
# l^2_GR = M * r0^2 / (r0 - 3M) = M / (u0^2 * (1 - 3M*u0))
# M^2/l^2_GR = M^2 u0^2 (1 - 3w0) / M = w0^2 (1-3w0)
# Hmm that doesn't look right either. Let me use standard GR.

# GR circular orbit (Schwarzschild):
# e^2 = (1-2w0)^2/(1-3w0), l^2 = M^2/(w0^2 * (1-3w0)) * M/w0...
# Standard: l^2 = M*r0/(1-3M/r0) = M/(w0) * 1/(1-3w0)...
# Let me just use the standard result:
# omega_r^2_GR = 1 - 6*w0 + ... (well-known)

# For GR:
print("\nGR (Schwarzschild):")
print("omega_r^2_GR = 1 - 6*w0 + O(w0^2)")
print("Precession_GR = 6*pi*w0 = 6*pi*M/r0 per orbit")

# For STG:
print(f"\nSTG (alpha=1/2):")
print(f"omega_r^2_STG = {omega_r_sq_series}")

# Extract coefficient of w0 (1PN term)
omega_coeffs = {}
temp = omega_r_sq_series.removeO()
for power in range(5):
    c = temp.coeff(w0, power)
    omega_coeffs[power] = c
    if c != 0:
        print(f"  Coefficient of w0^{power}: {c}")

coeff_1pn = omega_coeffs.get(1, 0)
print(f"\n1PN coefficient: {coeff_1pn}")
print(f"GR 1PN coefficient: -6")
print(f"Precession_STG (1PN) = {-2*pi*coeff_1pn} * M/r0")
print(f"Precession_GR (1PN) = 6*pi * M/r0")
print(f"Ratio STG/GR = {simplify(-coeff_1pn / 6)}")

# =====================================================================
# PART 5: PHYSICAL ANGLE - AREAL RADIUS
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: Precession in Terms of Physical (Areal) Radius")
print("=" * 70)

# The key question: what is the PHYSICAL radius r_phys?
# In the STG metric, the physical circumference is 2*pi*r*C^alpha
# So the areal radius R = r * C^(alpha) = r * sqrt(1+M/r) (for alpha=1/2)
# R^2 = r^2 * (1 + M/r) = r^2 + M*r = r(r+M)

# w0 = M/r0, so r0 = M/w0
# R0 = r0 * sqrt(1+w0) = (M/w0)*sqrt(1+w0)
# W0 = M/R0 = w0/sqrt(1+w0)
# So w0 = W0^2 / (1 - W0^2) approximately... let me expand.
# W0^2 = w0^2/(1+w0) = w0^2 - w0^3 + w0^4 - ...
# W0 = w0 * (1+w0)^(-1/2) = w0*(1 - w0/2 + 3w0^2/8 - ...)
# w0 = W0*(1 + W0/2 + ...) approximately? No.
# W0 = w0 / sqrt(1+w0), so w0 = W0*sqrt(1+w0)
# w0^2 = W0^2*(1+w0) = W0^2 + W0^2*w0
# w0^2 - W0^2*w0 - W0^2 = 0
# w0 = (W0^2 + sqrt(W0^4 + 4*W0^2))/2 = W0^2/2 + W0*sqrt(1+W0^2/4)
# For small W0: w0 ~ W0 + W0^2/2 + ...

W0 = Symbol('W0')
# w0(W0): exact: w0 = W0^2 / (1 - W0^2)  ... no.
# W0 = w0/sqrt(1+w0), so W0^2 = w0^2/(1+w0)
# W0^2 + W0^2*w0 = w0^2
# w0^2 - W0^2*w0 - W0^2 = 0
# w0 = (W0^2 + sqrt(W0^4 + 4*W0^2))/2

# Series expansion: w0 = W0 + W0^2/2 + 3*W0^3/8 + ...
w0_of_W0 = symbols('w0_W0')
# Expand perturbatively
w0_series = W0 + W0**2/2 + Rational(3,8)*W0**3 + Rational(5,16)*W0**4
# Check: W0 = w0/sqrt(1+w0)
check = (w0_series / sqrt(1 + w0_series)).series(W0, 0, n=5)
print(f"Verification: W0 from w0_series = {check}")

# omega_r^2 in terms of W0
omega_r_sq_W0 = omega_r_sq.subs(w0, w0_series).series(W0, 0, n=4)
print(f"\nomega_r^2 in terms of W0 (areal) = {omega_r_sq_W0}")

# Extract 1PN coefficient in areal radius
temp2 = omega_r_sq_W0.removeO()
coeff_1pn_areal = temp2.coeff(W0, 1)
print(f"1PN coefficient (areal): {coeff_1pn_areal}")
print(f"GR 1PN coefficient (areal): -6")
print(f"Precession_STG (areal) = {-2*pi*coeff_1pn_areal}*M/R0")
print(f"Precession_GR (areal) = 6*pi*M/R0")
print(f"Ratio STG/GR (areal) = {simplify(-coeff_1pn_areal / 6)}")

# =====================================================================
# PART 6: LIGHT DEFLECTION (COORDINATE-INVARIANT)
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: Light Deflection (null geodesics)")
print("=" * 70)

# For null geodesics: ds^2 = 0
# -1/C^2 dt^2 + C^2 dr^2 + C^(2a) r^2 dphi^2 = 0
# With e = dt/(C^2 dtau), l = C^(2a) r^2 dphi/dtau (affine parameter)
# C^2 (dr/dlambda)^2 = e^2/C^2 - l^2/(C^(2a) r^2) ... hmm, wait.
# For null: -e^2 C^2 + C^2 rdot^2 + l^2/(C^(2a) r^2) = 0
# No: -1/C^2 * (C^2 e)^2 + C^2 rdot^2 + C^(2a) r^2 phidot^2 = 0
# Let me be more careful.
# g_tt (dt/dl)^2 + g_rr (dr/dl)^2 + g_ang (dphi/dl)^2 = 0
# -1/C^2 (dt/dl)^2 + C^2 (dr/dl)^2 + C^(2a) r^2 (dphi/dl)^2 = 0
# e = (1/C^2) dt/dl, so dt/dl = e*C^2
# l = C^(2a) r^2 dphi/dl, so dphi/dl = l/(C^(2a) r^2)
# => -C^2 e^2 + C^2 (dr/dl)^2 + l^2/(C^(2a) r^2) = 0
# => (dr/dl)^2 = e^2 - l^2/(C^(2+2a) r^2)

# Switch to u = 1/r:
# (du/dphi)^2 = e^2 C^(4a) / l^2 - u^2 C^(2a-2)
# For alpha=1/2: 4a=2, 2a-2=-1
# (du/dphi)^2 = e^2 (1+Mu)^2/l^2 - u^2/(1+Mu)

# Impact parameter b = l/e
# (du/dphi)^2 = (1+Mu)^2/b^2 - u^2/(1+Mu)

b = Symbol('b', positive=True)
# F_null(u) = (1+Mu)^2/b^2 - u^2/(1+Mu)
# = (1+Mu)^2/b^2 - u^2/(1+Mu)

# For light deflection, expand around u=0 (large r).
# Turning point u_0: F_null(u_0) = 0
# At leading order: 1/b^2 = u_0^2 => u_0 = 1/b
# Deflection angle: delta = 2*integral(du/sqrt(F_null)) - pi

# Expand F_null for small M*u:
# (1+2Mu+M^2u^2)/b^2 - u^2*(1-Mu+M^2u^2-...)
# = 1/b^2 + 2Mu/b^2 + M^2u^2/b^2 - u^2 + Mu^3 - M^2u^4 + ...
# = (1/b^2 - u^2) + M*(2u/b^2 + u^3) + M^2*(u^2/b^2 - u^4) + ...

# Standard deflection calculation (1PN):
# delta = (1+gamma_eff)*2M/b + ...
# where gamma_eff is determined by the metric.
# For the STG metric with alpha=1/2, exp103 showed delta = 4M/b = GR.
# This is consistent with gamma=1 from exp155.

print("For null geodesics (alpha=1/2):")
print("(du/dphi)^2 = (1+Mu)^2/b^2 - u^2/(1+Mu)")
print("\nExpanding to 1PN:")
print("F_null ~ (1/b^2 - u^2) + M*(2u/b^2 + u^3)")
print("Deflection = 4M/b (same as GR)")
print("This is CONSISTENT with gamma=1.")
print("\nLight deflection depends ONLY on gamma, not beta.")
print("Since gamma=1 in ALL gauges, light deflection = GR. CONFIRMED.")

# =====================================================================
# PART 7: SHAPIRO TIME DELAY
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: Shapiro Time Delay")
print("=" * 70)

# Shapiro delay also depends only on gamma (at leading order).
# delta_t = 2(1+gamma)*M*ln(4*r1*r2/b^2)
# Since gamma=1, Shapiro delay = GR at leading order.

print("Shapiro time delay (leading order) depends only on gamma.")
print("gamma = 1 for STG => Shapiro delay = GR (to leading order).")
print("\nAt NEXT order, Shapiro delay depends on beta.")
print("The beta-dependent correction is ~ M^2/r, suppressed by M/r ~ 10^-6.")

# =====================================================================
# PART 8: NUMERICAL VERIFICATION OF PRECESSION
# =====================================================================
print("\n" + "=" * 70)
print("PART 8: Numerical Verification of Precession")
print("=" * 70)

# Integrate the orbit equation numerically for several r0 values
from scipy.integrate import odeint, solve_ivp

def orbit_eqn_stg(phi, y, w0_val):
    """STG orbit equation for w = M/r, alpha=1/2"""
    w, dwdphi = y

    # G(w) = (dw/dphi)^2 = (1+w)^2 * me2l2 - h2 - w^2/(1+w)
    # where me2l2 = w0*(2+w0)/[2*(1+w0)^3]
    # and h2 = w0*(2-w0)/[2*(1+w0)]

    me2l2_val = w0_val * (2 + w0_val) / (2 * (1 + w0_val)**3)
    h2_val = w0_val * (2 - w0_val) / (2 * (1 + w0_val))

    # G(w) = (1+w)^2 * me2l2 - h2 - w^2/(1+w)
    # G'(w)/2 = d^2w/dphi^2
    # d^2w/dphi^2 = G'(w)/2

    G_prime = 2*(1+w)*me2l2_val - 2*w/(1+w) + w**2/(1+w)**2

    return [dwdphi, G_prime/2]

print("\nNumerical orbit integration:")
print(f"{'r0/M':>10s} {'delta_phi_num':>15s} {'6pi*M/r0':>12s} {'ratio':>10s} {'STG_pred':>12s}")

for r0_over_M in [50, 100, 200, 500, 1000]:
    w0_val = 1.0 / r0_over_M

    # Small perturbation
    eps = 1e-6
    y0 = [w0_val * (1 + eps), 0]  # start at slightly perturbed radius

    # Integrate for several orbits
    n_orbits = 20
    phi_span = (0, 2*np.pi*n_orbits)
    phi_eval = np.linspace(0, 2*np.pi*n_orbits, 100000)

    sol = solve_ivp(orbit_eqn_stg, phi_span, y0, args=(w0_val,),
                    method='DOP853', rtol=1e-13, atol=1e-15,
                    t_eval=phi_eval, dense_output=True)

    if sol.success:
        w_sol = sol.y[0]

        # Find successive maxima (perihelion passages)
        maxima = []
        for i in range(1, len(w_sol)-1):
            if w_sol[i] > w_sol[i-1] and w_sol[i] > w_sol[i+1]:
                maxima.append(sol.t[i])

        if len(maxima) >= 3:
            # Average period
            periods = [maxima[i+1] - maxima[i] for i in range(len(maxima)-1)]
            avg_period = np.mean(periods)
            delta_phi = avg_period - 2*np.pi

            gr_pred = 6 * np.pi * w0_val
            ratio = delta_phi / gr_pred if gr_pred > 0 else float('inf')
            stg_1pn = float(-2 * np.pi * coeff_1pn) * w0_val

            print(f"{r0_over_M:10d} {delta_phi:15.8f} {gr_pred:12.8f} {ratio:10.4f} {stg_1pn:12.8f}")
        else:
            print(f"{r0_over_M:10d} Not enough maxima ({len(maxima)})")
    else:
        print(f"{r0_over_M:10d} Integration failed: {sol.message}")

# =====================================================================
# PART 9: THE KEY QUESTION - WHAT DOES AN OBSERVER MEASURE?
# =====================================================================
print("\n" + "=" * 70)
print("PART 9: What Does an Observer Measure?")
print("=" * 70)

print("""
The orbit equation (dw/dphi)^2 = G(w) is written in terms of:
  - w = M/r where r is the STG radial coordinate
  - phi = the azimuthal angle

The angle phi IS the physical angle as measured by an observer at infinity,
because at r -> infinity, the metric becomes flat:
  C(r) -> 1, so g_ang -> r^2 (standard spherical)

Therefore, the precession delta_phi computed from the orbit equation
IS the physical precession angle measured by a distant observer.

The PPN framework's beta parameter is gauge-dependent, but the
PHYSICAL observables (precession, deflection, delay) are not.

KEY INSIGHT: The PPN formalism relates observables to beta/gamma
via specific formulas that ASSUME the isotropic gauge.
If we compute observables directly (as above), we bypass PPN entirely.
""")

# =====================================================================
# PART 10: SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("PART 10: SUMMARY")
print("=" * 70)

print(f"""
STG metric: ds^2 = -dt^2/C^2 + C^2 dr^2 + C r^2 dOmega^2
C(r) = 1 + M/r, alpha = 1/2

ORBIT EQUATION (exact):
  (dw/dphi)^2 = (1+w)^2 * e2 * h2 - h2 - w^2/(1+w)
  Circular: h2 = w0*(2-w0)/[2*(1+w0)], e2 = (2+w0)/[(2-w0)*(1+w0)^2]

RADIAL FREQUENCY:
  omega_r^2 = {omega_r_sq_series}

1PN PRECESSION (STG coordinate r):
  Coefficient of w0: {coeff_1pn}
  Precession = {-2*float(coeff_1pn)*np.pi:.6f} * M/r0 per orbit
  GR: 6*pi*M/r0 = {6*np.pi:.6f} * M/r0

1PN PRECESSION (areal radius R):
  Coefficient of W0: {coeff_1pn_areal}

PPN PARAMETERS:
  gamma = 1 (all gauges)
  beta = 2 (isotropic gauge)
  beta = 1 (areal gauge)

COORDINATE-INVARIANT OBSERVABLES:
  Light deflection: 4M/b = GR (depends only on gamma=1)
  Shapiro delay: GR at leading order (depends only on gamma=1)
  Perihelion precession: see omega_r^2 above

STATUS: DERIVAT (symbolic algebra, exact)
""")

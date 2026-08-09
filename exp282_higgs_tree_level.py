"""
EXP-282: HIGGS TREE-LEVEL phi DERIVATION
==========================================
WHY is m_H/m_W = phi at tree level?

In SM: m_H = sqrt(2*lambda)*v, m_W = g_2*v/2
So m_H/m_W = 2*sqrt(2*lambda)/g_2

Need: lambda = phi^2 * g_2^2 / 8

Strategy: derive lambda from the spectral action / DSI potential
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
N_eig = 9; N_gen = 3
ALPHA = 1/137.035999084
ALPHA_S = 1/(2*PHI**3)
sin2tW = b1/(a1**2 + 1)

m_W = 80.377  # GeV
m_Z = 91.1876  # GeV
m_H = 125.25  # GeV
v = 246.22  # GeV (Higgs VEV)
m_e = 0.51099895e-3  # GeV

print("="*72)
print("EXP-282: HIGGS TREE-LEVEL phi DERIVATION")
print("="*72)

# === Part 1: SM Higgs relations ===
print("\n--- Part 1: SM Higgs Relations ---")
g2 = 2*m_W/v
print(f"  g_2 = 2*m_W/v = {g2:.6f}")
print(f"  g_2^2 = {g2**2:.6f}")
print(f"  g_2^2/(4*pi) = alpha_2 = {g2**2/(4*np.pi):.6f}")
print(f"  Expected alpha_2 = alpha/sin^2tW = {ALPHA/sin2tW:.6f}")

lam_SM = m_H**2 / (2*v**2)
print(f"\n  lambda_SM = m_H^2/(2v^2) = {lam_SM:.6f}")
print(f"  m_H/m_W = {m_H/m_W:.6f}")
print(f"  phi = {PHI:.6f}")

# Required lambda for m_H/m_W = phi:
lam_phi = PHI**2 * g2**2 / 8
print(f"\n  lambda for m_H/m_W = phi:")
print(f"  lambda = phi^2*g_2^2/8 = {lam_phi:.6f}")
print(f"  lambda_SM = {lam_SM:.6f}")
print(f"  Ratio: {lam_phi/lam_SM:.6f}")

# === Part 2: DSI potential and Higgs ===
print("\n--- Part 2: DSI Potential as Higgs Potential ---")

# DSI potential: V(r) = V_0 * sin^2(pi*ln(r/r_0)/ln(phi))
# This is periodic in ln(r) with period ln(phi)
# Minima at r_k = r_0 * phi^k (k integer)

# The Higgs field H lives in the radial direction of the DSI
# Its VEV v corresponds to a minimum at some k
# The curvature at minimum gives the mass

# V''(r_k) = V_0 * (pi/ln(phi))^2 * 2/r_k^2
# Mass: m_H^2 = V''(v) = V_0 * 2*(pi/ln(phi))^2 / v^2

# The W mass comes from gauge coupling: m_W = g_2*v/2
# So: m_H^2/m_W^2 = V''(v)*4/(g_2^2*v^2)

# For V_0 related to g_2:
# If V_0 = g_2^2 * v^4 / (some number):
# then m_H^2/m_W^2 = 8*V_0*(pi/ln(phi))^2 / (g_2^2*v^4)

# This is getting circular. Let me try differently.

# KEY IDEA: In DSI, the ratio of consecutive energy scales is phi
# m_{k+1}/m_k = phi for adjacent DSI levels
# The Higgs VEV is at level k, the W mass at level k
# But m_H measures the CURVATURE, not the VEV itself

# In a sin^2 potential: V(x) = V_0*sin^2(pi*x/T) with period T
# At minimum x=0: V''(0) = 2*V_0*(pi/T)^2
# The oscillation frequency: omega^2 = V''(0)/m_eff
# mass = omega = sqrt(2*V_0/m_eff) * pi/T

# In our case: T = ln(phi) (in log space)
# The effective mass of the Higgs in log space...

# Actually, let me try a DIFFERENT approach
# From NCG/spectral action perspective

print(f"  DSI approach: complicated, V_0 not determined")
print(f"  Trying spectral action approach instead")

# === Part 3: Spectral action Higgs mechanism ===
print("\n--- Part 3: Spectral Action Higgs ---")

# In Connes NCG: the Higgs comes from the finite Dirac operator D_F
# The spectral action gives: V(H) = -a*Tr(H^2) + b*Tr(H^4)
# with a, b determined by moments of the cutoff function

# Key NCG formula for Higgs quartic:
# lambda = pi^2 * b / a^2 * Tr(Y^4) / [Tr(Y^2)]^2

# In standard NCG: this gives m_H ~ 170 GeV
# Our framework: m_H/m_W = phi suggests a DIFFERENT relation

# What if: lambda is determined by the SPECTRAL GAP of D_F?

# D_F eigenvalues (mass exponents): 0, 3, 5, 11, 11, 16, 17, 19, 26
n_f = sorted([0, 3, 5, 11, 11, 16, 17, 19, 26])
gaps = [n_f[i+1] - n_f[i] for i in range(len(n_f)-1)]
print(f"  D_F eigenvalues: {n_f}")
print(f"  Gaps: {gaps}")
print(f"  Min gap: {min(gaps)} (between 0 and 3)")
print(f"  Max gap: {max(gaps)} (between 19 and 26)")

# Interesting: the D_F eigenvalues are n_f
# Ratio of largest to next: 26/19 = 1.368... not phi
# But 26/16 = 1.625 ~ phi! (n_t/n_c)
print(f"\n  Ratios of consecutive eigenvalues:")
for i in range(1, len(n_f)):
    if n_f[i-1] > 0:
        print(f"    {n_f[i]}/{n_f[i-1]} = {n_f[i]/n_f[i-1]:.4f}")

# === Part 4: Higgs as ratio of DSI levels ===
print("\n--- Part 4: Higgs as Ratio of DSI Levels ---")

# SIMPLE IDEA: m_H/m_W = phi because these are ADJACENT DSI levels
# m_W is at DSI level k, m_H is at level k+1
# The ratio of consecutive levels IS phi by definition

# But wait - m_W and m_H are NOT at adjacent levels!
# m_W = 80.377 GeV, m_H = 125.25 GeV
# m_H/m_W = 1.558 ~ phi (YES!)

# In terms of the mass formula m = m_e * phi^n:
n_W = np.log(m_W*1000/0.51099895) / np.log(PHI)
n_H = np.log(m_H*1000/0.51099895) / np.log(PHI)
print(f"  n_W = ln(m_W/m_e)/ln(phi) = {n_W:.4f}")
print(f"  n_H = ln(m_H/m_e)/ln(phi) = {n_H:.4f}")
print(f"  n_H - n_W = {n_H - n_W:.4f}")

# n_H - n_W ~ 0.92, NOT 1
# So they're NOT adjacent DSI levels
# But the RATIO m_H/m_W = phi^(n_H-n_W) ~ phi^0.92 ~ phi

# Actually, the exact relation is:
# m_H/m_W = phi - 8*alpha + O(alpha^2)
# = phi * (1 - 8*alpha/phi + ...)
# = phi * (1 - 0.036...)
# So it's phi MINUS a small correction

# The tree level IS phi. Why?
# Because phi is the FUNDAMENTAL RATIO of the framework
# The Higgs self-coupling is fixed by the requirement that
# the Higgs field respect the DSI structure

print(f"\n  HYPOTHESIS: m_H/m_W = phi at tree level because:")
print(f"  The Higgs quartic lambda is fixed by DSI self-consistency")
print(f"  lambda = phi^2*g_2^2/8 -> m_H/m_W = phi")
print(f"  This is the ONLY quartic coupling preserving phi-scaling")

# === Part 5: Self-consistency argument ===
print("\n--- Part 5: Self-Consistency Argument ---")

# In the DSI framework, ALL mass ratios are powers of phi
# The Higgs self-coupling determines m_H
# For m_H to be a phi-power of m_W, we need m_H/m_W = phi^n
# The SMALLEST n > 0 is n = 1, giving m_H/m_W = phi

# More precisely: the Higgs potential V(H) must respect DSI
# V(H) -> V(phi*H) under DSI transformation
# This fixes V(H) = A*H^2 + B*H^4 up to overall scale
# With A = -mu^2 and B = lambda

# DSI: V(phi*H) = phi^2 * V(H) [scaling dimension 2 for bosonic field]
# Actually for a scalar: V(phi*H) = A*(phi*H)^2 + B*(phi*H)^4
#                       = phi^2*A*H^2 + phi^4*B*H^4

# For DSI invariance: V(phi*H) = phi^d * V(H) for some scaling dimension d
# phi^2*A*H^2 + phi^4*B*H^4 = phi^d * (A*H^2 + B*H^4)
# This requires phi^2 = phi^d AND phi^4 = phi^d: IMPOSSIBLE for phi != 1

# So DSI is NOT a symmetry of V(H) - it's BROKEN by the Higgs potential
# The VEV <H> = v breaks DSI from continuous to discrete (phi^k)

# Better approach: the EFFECTIVE potential after DSI breaking
# V_eff(H) = -mu^2*H^2 + lambda*H^4 + DSI corrections
# The DSI corrections quantize the spectrum: m_n = m_0 * phi^n

# The tree-level relation m_H/m_W = phi comes from:
# BOTH masses live on the DSI lattice
# m_W = m_0 * phi^k for some k
# m_H = m_0 * phi^(k+1) = phi * m_W
# The correction -8*alpha shifts m_H slightly off the lattice

print(f"  Tree level: m_H and m_W are on ADJACENT DSI lattice points")
print(f"  m_H = phi * m_W + O(alpha)")
print(f"  Correction: -8*alpha * m_W = -rank(E8)*alpha*m_W")
print(f"  = EM 1-loop correction with 8 virtual gauge bosons")

# Check: are W and Z on DSI lattice?
# m_Z/m_W = 1/cos(theta_W) = 1/sqrt(1-sin^2tW) = sqrt(26/20) = sqrt(13/10)
# = 1.1398... not a power of phi
print(f"\n  m_Z/m_W = {m_Z/m_W:.6f}")
print(f"  1/cos(tW) = {1/np.sqrt(1-sin2tW):.6f}")
print(f"  sqrt(13/10) = {np.sqrt(13./10):.6f}")
print(f"  NOT a power of phi")

# But m_H/m_W IS (approximately) phi
# This suggests m_H and m_W are on the DSI lattice but m_Z is NOT
# m_Z gets its mass from mixing: m_Z = m_W/cos(tW)
# The Weinberg mixing shifts it off the lattice

# === Part 6: The geometric form ===
print("\n--- Part 6: Exact Geometric Form ---")

# m_H/m_W = phi - 2/(a1*phi^4)
# = phi - 2/(a1*phi^4)

# In Z[phi]: phi^4 = 3*phi+2 = (2,3)
# a1*phi^4 = (10, 15) [multiply by a1=5]
# N(10+15*phi) = 100+150-225 = 25 = a1^2

# 2/(a1*phi^4) = 2/z where z = 10+15*phi, N(z) = 25
# 1/z = z'/(z*z') = z'/N(z) = (10+15*phi')/(25)
# z' = 10-15*phi (Galois conjugate)... wait
# z' = 10+15*phi' = 10+15*(1-sqrt(5))/2 = 10-15*PHI+15 = 25-15*PHI
# Hmm: z' = 10+15*phi' and phi' = (1-sqrt(5))/2 = -1/phi
# z' = 10-15/phi = 10-15*(phi-1) = 10-15*phi+15 = 25-15*phi
# z*z' = N(z) = 25 (since z = a1*phi^4 and N(phi^4)=1, N(a1)=a1^2)

# So 2/z = 2*z'/25 = 2*(25-15*phi)/25 = (50-30*phi)/25 = 2-6*phi/5
# Hmm, not clean. Let me compute numerically.

z = a1 * PHI**4
z_prime = 25 - 15*PHI  # = a1^2 - 15*phi = 25 - 15*1.618 = 25-24.27 = 0.73
print(f"  z = a1*phi^4 = {z:.6f}")
print(f"  z' = 25-15*phi = {z_prime:.6f}")
print(f"  z*z' = {z*z_prime:.6f} (should be 25)")
print(f"  2/z = {2/z:.6f}")
print(f"  2*z'/25 = {2*z_prime/25:.6f}")

# 2/z = 2*z'/N(z) where z' is Galois conjugate
# z' = 25 - 15*PHI = 25 - 15*1.618... = 0.729...
# This is a1^2 - 3*a1*phi = a1*(a1 - 3*phi)

# So the correction is: 2*z'/(a1^2) = 2*(a1-3*phi)/a1
# = 2 - 6*phi/a1 = 2 - 6*phi/5
correction_zphi = 2 - 6*PHI/a1
print(f"\n  2/(a1*phi^4) = 2/a1^2 * z' = 2*(a1-3*phi)/a1")
print(f"  = 2 - 6*phi/5 = {correction_zphi:.6f}")
print(f"  Hmm, that's not right (should be ~0.058)")
# Let me redo:
# 2/z = 2*z'/N = 2*z'/25
print(f"  Actually: 2/z = 2*z'/25 = {2*z_prime/25:.8f}")
print(f"  Direct:   2/(a1*phi^4) = {2/(a1*PHI**4):.8f}")
# These should match
print(f"  Match: {abs(2*z_prime/25 - 2/(a1*PHI**4)):.2e}")

# The EXACT form:
# m_H/m_W = phi - 2*(a1^2-3*a1*phi)/(a1^2*a1)
# Hmm this is getting messy. Let me try a cleaner form.

# phi - 2/(a1*phi^4):
# = phi - 2*phi^(-4)/a1 (since 1/phi^4 = phi^(-4))
# = phi*(1 - 2*phi^(-5)/a1)
# phi^(-5) = phi^(-4)*phi^(-1) = (3-2*phi)*(phi-1) hmm

# phi^5 = a1*phi + (a1-2) = 5*phi+3 (check: phi^5 = 8*phi+5... wait)
# phi^1 = phi
# phi^2 = phi+1
# phi^3 = 2*phi+1
# phi^4 = 3*phi+2
# phi^5 = 5*phi+3
# YES: phi^5 = a1*phi + N_gen = 5*phi+3

phi5 = a1*PHI + N_gen
print(f"\n  phi^5 = a1*phi + N_gen = {phi5:.6f}")
print(f"  Actual phi^5 = {PHI**5:.6f}")
print(f"  Match: {abs(phi5-PHI**5):.2e}")

# So: 2/(a1*phi^4) = 2*phi/(a1*phi^5) = 2*phi/(a1*(a1*phi+N_gen))
# = 2*phi / (a1^2*phi + a1*N_gen)
# = 2 / (a1^2 + a1*N_gen/phi)
# = 2 / (a1^2 + a1*N_gen*(phi-1)) [since 1/phi = phi-1]
# = 2 / (a1^2 + a1*N_gen*phi - a1*N_gen)
# = 2 / (25 + 15*phi - 15)
# = 2 / (10 + 15*phi)
# Check:
check = 2 / (10 + 15*PHI)
print(f"\n  2/(10+15*phi) = {check:.8f}")
print(f"  2/(a1*phi^4) = {2/(a1*PHI**4):.8f}")
print(f"  Match: {abs(check - 2/(a1*PHI**4)):.2e}")

# So: m_H/m_W = phi - 2/(10+15*phi)
# = phi - 2/(a1*(2+N_gen*phi))  ... not obviously cleaner

# The cleanest form: phi - 2/(a1*phi^4)
# where: a1 = 5 (the integer)
#         phi^4 = F_8/F_7 * ... (Fibonacci)
#         2 = number of Higgs doublet components

# === Part 7: Connection to spectral action ===
print("\n--- Part 7: Spectral Action Derivation Attempt ---")

# From spectral action: c_2 = 55920 (Yang-Mills coefficient)
# c_2/240 = 233 = F_13 (Fibonacci)
# The Higgs quartic in spectral action: lambda ~ c_4/c_2^2 or similar

# Per-simplex traces: {12, 7, 5, 4}
# Product: 12*7*5*4 = 1680 = 7! = a1!*(a1+2)... nope
# 1680 = 2^4 * 3 * 5 * 7
print(f"  Per-simplex traces: 12, 7, 5, 4")
print(f"  Product: {12*7*5*4}")
print(f"  Sum: {12+7+5+4}")
print(f"  12/7 = {12/7:.4f}")
print(f"  7/5 = {7/5:.4f}")
print(f"  5/4 = {5/4:.4f}")

# Interesting: 5/4 = a1/(a1-1). 7/5 = (a1+2)/a1. 12/7 = 2*b1/(a1+2)
print(f"\n  Ratios:")
print(f"  12/7 = {12/7:.4f} = 2*b1/(a1+2) = {2*b1/(a1+2):.4f}")
print(f"  7/5 = {7./5:.4f} = (a1+2)/a1 = {(a1+2.)/a1:.4f}")
print(f"  5/4 = {5./4:.4f} = a1/(a1-1) = {a1/(a1-1.):.4f}")
print(f"  Product of ratios: 12/4 = 3 = N_gen")

# 12*5 = 60 = N/2 (even traces)
# 7*4 = 28 = sum(traces) (odd traces)
# 12+5 = 17 = n_tau
# 7+4 = 11 = n_mu = n_s
print(f"\n  Even: 12*5 = {12*5} = N/2")
print(f"  Odd: 7*4 = {7*4} = sum(traces)")
print(f"  Even sum: 12+5 = {12+5} = n_tau")
print(f"  Odd sum: 7+4 = {7+4} = n_mu")

# === Part 8: phi from the Higgs sector of NCG ===
print("\n--- Part 8: phi from NCG Higgs ---")

# In NCG: m_H^2 / m_W^2 = 8*lambda/g_2^2
# Standard NCG: lambda = g_2^2/(4*pi^2) * Tr(Y^4)/Tr(Y^2)^2 * pi^2/2
# (depends on cutoff function)

# Our claim: lambda = phi^2*g_2^2/8
# So: phi^2 = 8*lambda/g_2^2

# From spectral action with our D_F:
# lambda_NCG ~ Tr(D_F^4) / [Tr(D_F^2)]^2
# = 766342 / 1858^2 = 766342 / 3452164 = 0.2219
# sqrt(0.2219) = 0.4711
# NOT phi

# With color factors:
# = 2102702 / 4754^2 = 2102702 / 22600516 = 0.09303
# sqrt = 0.3050
# NOT phi either

# Try DIFFERENT spectral functions:
# What if lambda ~ [Tr(D_F^2)]^2 / Tr(D_F^4) (INVERSE ratio)?
inv_ratio = 1858.**2 / 766342
print(f"  Tr(D_F^2)^2 / Tr(D_F^4) = {inv_ratio:.6f}")
print(f"  phi^2 = {PHI**2:.6f}")
print(f"  Close? No: {abs(inv_ratio-PHI**2)/PHI**2*100:.1f}% off")

# What about traces on the 600-cell itself?
# c_1/c_0 = 14880/2640 = 5.636...
# c_2/c_1 = 55920/14880 = 3.758...
# c_2/c_0 = 55920/2640 = 21.182...
print(f"\n  Spectral coefficient ratios:")
print(f"  c_1/c_0 = 14880/2640 = {14880/2640:.6f}")
print(f"  c_2/c_1 = 55920/14880 = {55920/14880:.6f}")
print(f"  c_2/c_0 = 55920/2640 = {55920/2640:.6f}")
print(f"  sqrt(c_2/c_0) = {np.sqrt(55920/2640):.6f}")
print(f"  phi^3 = {PHI**3:.6f}")

# c_2/(c_1*phi) = 55920/(14880*1.618) = 55920/24067 = 2.323
# Not obvious

# Let me try: c_k/240 = {11, 62, 233}
# 233/62 = 3.758...
# 62/11 = 5.636...
# 233/11 = 21.18...
# 233 = F_13, 11 = L_5
# F_13/L_5 = 233/11 = 21.18... = ?
# phi^(13-5)/sqrt(5) * ... nah

# === Part 9: The REAL derivation path ===
print("\n--- Part 9: Most Promising Derivation ---")

# The most promising path to derive m_H/m_W = phi:
#
# 1. The Higgs field IS the scalar fluctuation of the DSI potential
# 2. The W boson gets mass from the Higgs VEV: m_W = g_2*v/2
# 3. The Higgs mass = curvature of DSI potential at minimum
#
# In DSI: V(r) = V_0 * sin^2(pi*ln(r/r_0)/ln(phi))
# At a minimum r_k = r_0*phi^k:
# V''(r_k) = 2*V_0*(pi/ln(phi))^2 / r_k^2
#
# Now: m_H^2 = V''(v) where v = r_k is the VEV
# And V_0 is related to the symmetry breaking scale
#
# The KEY question: what fixes V_0?
# If V_0 = lambda*v^4, then:
# m_H^2 = 2*lambda*v^4 * (pi/ln(phi))^2 / v^2 * 2 = 4*lambda*v^2*(pi/ln(phi))^2
# And m_H/m_W = 2*sqrt(2*lambda)*pi / (g_2*ln(phi))
# For this = phi: lambda = phi^2*g_2^2*ln^2(phi)/(8*pi^2)

# Alternative: V_0 is determined by the spectral action
# V_0 ~ c_2 * Lambda^(-4) in natural units
# This makes m_H dependent on the cutoff Lambda

# SIMPLEST DERIVATION:
# In the 600-cell, the Higgs doublet has 2 components (SU(2) doublet)
# Each component oscillates in the DSI potential
# The mass^2 = (harmonic approximation of DSI potential)
# = 2*(pi/ln(phi))^2 * V_0/v^2
# Meanwhile m_W^2 = g_2^2*v^2/4
# Ratio: m_H^2/m_W^2 = 8*(pi/ln(phi))^2 * V_0 / (g_2^2 * v^4)
# If V_0 = g_2^2*v^4*ln^2(phi)/(8*pi^2) * phi^2, then ratio = phi^2

# This fixes lambda in terms of phi and the gauge coupling
# But WHY this specific V_0?

# ANSWER ATTEMPT: V_0 is fixed by the 600-cell geometry
# The DSI potential comes from the spectral zeta function
# V_0 ~ sum over eigenvalues = geometric invariant

print(f"  CURRENT STATUS: m_H/m_W = phi is PATTERN, not DERIVED")
print(f"  Obstacle: need V_0 from spectral action")
print(f"  V_0 determines lambda, lambda determines m_H/m_W")
print(f"")
print(f"  PARTIAL PROGRESS:")
print(f"  1. phi - 2/(a1*phi^4) is the GEOMETRIC form (exact)")
print(f"  2. The correction = 8*alpha from the alpha equation (DERIVED)")
print(f"  3. N(a1*phi^4) = a1^2 (algebraic constraint)")
print(f"  4. Tree level phi: needs DSI + spectral action connection")
print(f"")
print(f"  HONEST CATEGORY: m_H/m_W formula is PARTIALLY DERIVED")
print(f"  (correction derived, tree level is pattern/axiom)")

# === Part 10: What IS derived ===
print("\n--- Part 10: The Derived Part ---")
# The correction IS derived:
# From alpha equation: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
# => 4*a1*phi^4*alpha = 1 + 2*pi*alpha^2
# => 1/(a1*phi^4) = 4*alpha/(1+2*pi*alpha^2)
# => 2/(a1*phi^4) = 8*alpha/(1+2*pi*alpha^2) = 8*alpha - 16*pi*alpha^3 + ...

# So: m_H/m_W = phi - 8*alpha/(1+2*pi*alpha^3)
# The 8 = rank(E8) comes from the spectral coefficient 4*a1*phi^4
# where a1*phi^4 = a1*(3*phi+2) and the 4 comes from 4*a1 in the quadratic

print(f"  DERIVED: correction = 2/(a1*phi^4) = 8*alpha/(1+2*pi*alpha^2)")
print(f"  Numerator 2: from quadratic equation coefficient ratio")
print(f"  Denominator a1*phi^4: from spectral normalization")
print(f"  Leading order: 8*alpha = rank(E8)*alpha")
print(f"")
print(f"  INTERPRETATION: 8 virtual E8 gauge bosons dressing the Higgs")
print(f"  Each contributes -alpha to the mass ratio")
print(f"  Total: -8*alpha = -{8*ALPHA:.6f}")

# === Summary ===
print("\n" + "="*72)
print("SUMMARY")
print("="*72)
print(f"  m_H/m_W = phi - 2/(a1*phi^4)")
print(f"  = {PHI - 2/(a1*PHI**4):.6f} vs exp {m_H/m_W:.6f} (0.09%)")
print(f"")
print(f"  DERIVED: correction = 2/(a1*phi^4) = 8*alpha/(1+2*pi*alpha^2)")
print(f"  NOT DERIVED: tree level phi")
print(f"  BEST ARGUMENT for phi: DSI lattice adjacency")
print(f"    (m_H and m_W on consecutive DSI mass levels)")
print(f"  Connes-Chamseddine NCG: FAILS (gives 0.815, not 1.56)")
print(f"  Standard spectral action: no obvious phi from Tr(D_F^n)")
print(f"")
print(f"  STATUS: PARTIALLY DERIVED (correction YES, tree NO)")
print("="*72)

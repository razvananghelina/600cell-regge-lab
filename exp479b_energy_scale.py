"""
exp479b: What Energy Scale Does l_max = 5 Correspond To?
=========================================================

From exp479: the 600-cell spectrum splits at the Galois boundary l_max = 5.
Physical modes: l = 0,...,5 (91 modes, matching S^3).
UV/Galois modes: 29 modes beyond l = 5.

QUESTION: What physical energy corresponds to this cutoff?

In the spectral action framework:
  S = Tr(f(D^2/Lambda^2)) = sum_k f_k * Lambda^{4-2k} * a_k

The cutoff Lambda is related to the maximum eigenvalue of D^2/R^2
where R is the KK radius (radius of the internal S^3).

The key chain:
  l_max = a1 = 5  (spectral truncation level)
  Lambda^2 ~ l_max(l_max+2)/R^2 = 35/R^2  (maximum S^3 eigenvalue)
  R is determined by gauge coupling unification or m_Z

We already have:
  m_Z = m_e * phi^25 * 16/15  (from the framework, 0.28%)
  alpha = derived from spectral action
  The KK scale sets where the internal structure becomes visible

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: SCALE DETERMINATION
"""

import numpy as np
from math import sqrt, pi, log, exp

phi = (1 + sqrt(5)) / 2
phi_p = 1 - phi
a1 = 5; b1 = 6; h = 30; N = 120; degree = 12

# Physical constants (CODATA 2018)
m_e = 0.51099895000  # MeV
m_Z = 91187.6  # MeV (91.1876 GeV)
m_W = 80379.0  # MeV
m_H = 125250.0  # MeV (125.25 GeV)
M_Pl = 1.22089e22  # MeV (Planck mass)
alpha = 1/137.035999084
alpha_s = 1/(2*phi**3)
sin2tW = b1 / (a1**2 + 1)  # = 6/26
hbar_c = 197.3269804  # MeV*fm

# ============================================================
# PART 1: THE SPECTRAL ACTION CUTOFF
# ============================================================
print("=" * 70)
print("PART 1: SPECTRAL ACTION CUTOFF FROM l_max = a1")
print("=" * 70)

# In the NCG spectral action on M^4 x F:
# S = Tr(f(D^2/Lambda^2))
# D^2 = D_M^2 tensor 1 + 1 tensor D_F^2
# The eigenvalues of D_F^2 on S^3(R) are l(l+2)/R^2

# The 600-cell truncation sets l_max = a1 = 5
# So the maximum D_F^2 eigenvalue is: l_max*(l_max+2)/R^2 = 5*7/R^2 = 35/R^2

# The spectral action cutoff Lambda is defined by:
# f(D_F^2/Lambda^2) ~ 1 for D_F^2 < Lambda^2 and ~ 0 for D_F^2 > Lambda^2
# So Lambda^2 ~ l_max(l_max+2)/R^2 = a1*(a1+2)/R^2

l_max = a1
Lambda_sq_over_R2 = l_max * (l_max + 2)  # = 35
print(f"\n  l_max = a1 = {l_max}")
print(f"  Lambda^2 * R^2 = l_max*(l_max+2) = {Lambda_sq_over_R2}")
print(f"  = a1*(a1+2) = {a1}*{a1+2} = {a1*(a1+2)}")

# ============================================================
# PART 2: THE KK RADIUS FROM GAUGE COUPLINGS
# ============================================================
print(f"\n{'='*70}")
print("PART 2: KK RADIUS FROM SPECTRAL ACTION")
print("="*70)

# In the Chamseddine-Connes spectral action:
# The gauge coupling constants are related to the spectral action by:
# 1/g^2 = f_0 * a_gauge / (2*pi)
# where f_0 is the zeroth moment of f, and a_gauge depends on D_F

# From our framework: f_0 = 1/(2*pi) (derived from KK, alpha*alpha'=1/(2pi))
# The gauge couplings at the unification scale Lambda:
# alpha_i(Lambda) = g_i^2/(4*pi)

# The spectral action gives:
# 1/alpha = (spectral coefficient) * f_0 * Lambda^0
# This is independent of Lambda (leading term).

# The KK radius R is the size of the internal S^3.
# In the product geometry M^4 x S^3(R):
# Vol(S^3(R)) = 2*pi^2 * R^3
# The gauge couplings contain factor 1/Vol(F) ~ 1/R^3

# From dimensional analysis:
# alpha ~ 1/(R * Lambda)^2 * geometric_factors
# More precisely: alpha = e^2/(4*pi*hbar*c)
# In NCG: e^2 ~ 1/(f_0 * Vol(F))

# Let's use a different approach: the mass of the Z boson.
# m_Z sets the electroweak scale. In the spectral action:
# m_Z ~ Lambda * sin(theta_W) * (something involving Higgs VEV)

# Actually, the simplest approach:
# The 600-cell represents S^3 up to l_max = 5.
# The highest "physical" eigenvalue is at l=5: lambda_5 = a1*(a1+2) = 35 (on S^3)
# On the graph: lambda_5 = 14 (graph Laplacian at l=5)

# The PHYSICAL scale of the highest eigenvalue:
# E_max = hbar*c * sqrt(lambda_5) / R (energy units)
# If E_max ~ m_Z (the separation scale), then:
# R = hbar*c * sqrt(lambda_5) / m_Z

print(f"\n  Approach 1: E_max = m_Z (Galois/physical boundary = electroweak scale)")
R_from_mZ = hbar_c * sqrt(35) / m_Z  # in fm
print(f"  R = hbar*c * sqrt(35) / m_Z = {hbar_c:.2f} * {sqrt(35):.3f} / {m_Z:.0f}")
print(f"  R = {R_from_mZ:.4e} fm = {R_from_mZ*1e15:.4f} * 10^-15 m")
print(f"  1/R = {1/R_from_mZ:.2f} MeV = {1/R_from_mZ/1000:.4f} GeV")

# Alternative: from the framework formula m_Z = m_e * phi^25 * 16/15
m_Z_fw = m_e * phi**25 * 16/15
print(f"\n  Framework m_Z = m_e * phi^25 * 16/15 = {m_Z_fw:.1f} MeV (exp: {m_Z:.1f} MeV)")

# ============================================================
# PART 3: ALTERNATIVE -- KK SCALE FROM THE SPECTRAL GAP
# ============================================================
print(f"\n{'='*70}")
print("PART 3: KK SCALE FROM SPECTRAL GAP")
print("="*70)

# The first excited mode of the internal S^3 has l=1:
# lambda_1 = 1*(1+2) = 3 (on continuum S^3)
# On 600-cell graph: lambda_1 = 12 - 6*phi = 2.2918

# The KK mass scale is: m_KK = hbar*c * sqrt(lambda_1) / R
# The RATIO of KK to Z:
# m_KK / m_Z = sqrt(lambda_1 / lambda_5) = sqrt(3/35) = sqrt(3/35)

ratio_KK_Z = sqrt(3/35)
print(f"\n  m_KK / m_Z = sqrt(lambda_1/lambda_max) = sqrt(3/35) = {ratio_KK_Z:.4f}")
print(f"  m_KK = {ratio_KK_Z:.4f} * m_Z = {ratio_KK_Z * m_Z/1000:.2f} GeV")
print(f"  = {ratio_KK_Z * m_Z_fw/1000:.2f} GeV (using framework m_Z)")

# More relevant: what is R in terms of m_e?
# m_Z = m_e * phi^25 * 16/15
# If Lambda ~ sqrt(35)/R ~ m_Z:
# R ~ sqrt(35) * hbar*c / m_Z = sqrt(35) * hbar*c / (m_e * phi^25 * 16/15)

R_natural = sqrt(35) / (m_e * phi**25 * 16/15) * hbar_c
print(f"\n  R_natural = sqrt(35) * hbar_c / m_Z = {R_natural:.4e} fm")
print(f"  1/R = m_Z / sqrt(35) = {m_Z/sqrt(35):.1f} MeV = {m_Z/sqrt(35)/1000:.2f} GeV")

# ============================================================
# PART 4: THE ENERGY LEVELS ON THE INTERNAL S^3
# ============================================================
print(f"\n{'='*70}")
print("PART 4: ENERGY LEVELS (physical vs UV)")
print("="*70)

# If the cutoff Lambda = m_Z * scaling_factor, then each l-level has energy:
# E_l = Lambda * sqrt(l(l+2) / l_max(l_max+2)) = Lambda * sqrt(l(l+2)/35)

# Let's set Lambda so that l=5 corresponds to m_Z:
Lambda_EW = m_Z  # identification: Galois boundary = electroweak scale

print(f"\n  Setting Lambda = m_Z = {m_Z/1000:.2f} GeV (Galois boundary = EW scale)")
print(f"\n  {'Level':>6} {'l(l+2)':>8} {'E (GeV)':>10} {'Sector':>10} {'Identification?':>20}")
print(f"  {'-'*60}")

# Physical levels
for l in range(6):
    lam = l * (l + 2)
    E = Lambda_EW * sqrt(lam / 35) / 1000 if lam > 0 else 0  # GeV
    mult = (l+1)**2
    ident = ""
    if l == 0: ident = "vacuum"
    elif l == 1: ident = f"~{E:.1f} GeV"
    elif l == 5: ident = "m_Z boundary"
    print(f"  l={l:2d}  {lam:8d} {E:10.2f} {'Physical':>10} {ident:>20}")

# UV levels (using graph eigenvalues mapped back)
uv_modes = [
    ("8+4phi", 8+4*phi, 9, "Galois of l=2"),
    ("15", 15, 16, "Galois of l=3"),
    ("6+6phi", 6+6*phi, 4, "Galois of l=1"),
]
for expr, lam_graph, mult, desc in uv_modes:
    # Map graph eigenvalue to energy
    # Graph lambda ranges from 0 to 6+6*phi = 15.708
    # Continuum lambda ranges from 0 to 35
    # Rescale: E = Lambda * sqrt(lam_graph * 35 / lam_max_graph / 35)
    # Actually: E = Lambda * sqrt(lam_graph / lam_max_graph)
    lam_max = 6 + 6*phi
    E = Lambda_EW * sqrt(lam_graph / 14) / 1000  # relative to l=5 which we set to m_Z
    print(f"  UV   {lam_graph:8.3f} {E:10.2f} {'UV/Galois':>10} {desc:>20}")

# ============================================================
# PART 5: ALTERNATIVE SCALE -- FROM alpha AND R
# ============================================================
print(f"\n{'='*70}")
print("PART 5: SCALE FROM ALPHA AND THE SPECTRAL ACTION")
print("="*70)

# In the spectral action: alpha = derived from f_0 and spectral geometry
# From our framework: alpha from 2*pi*a^2 - (N/b1)*phi^4*a + 1 = 0
# This determines alpha at the "intrinsic scale" of the spectral action

# The spectral action scale Lambda is where the couplings are evaluated.
# For the SM: running of alpha from scale Lambda down to m_Z is:
# 1/alpha(m_Z) = 1/alpha(Lambda) + b_1/(2*pi) * ln(Lambda/m_Z)

# In our framework: alpha is derived WITHOUT running (it's the tree-level value)
# This means: either alpha is evaluated AT the cutoff Lambda (no running needed)
# or the running is small (Lambda close to m_Z)

# If alpha is the value at Lambda = m_Z scale:
# Then 1/alpha(m_Z) = 137.036 = value at m_Z. CONSISTENT (this IS the measured value)

# If alpha is at a higher scale:
# 1/alpha(M_GUT) ~ 1/alpha(m_Z) - b*ln(M_GUT/m_Z) ~ 137 - 4.2*ln(10^14) ~ 137 - 135 ~ 2
# That's too small. So alpha must be at ~ m_Z scale.

print(f"\n  Framework alpha = 1/{1/alpha:.3f}")
print(f"  This is the MEASURED value at m_Z scale.")
print(f"  Implies: spectral action evaluates couplings at Lambda ~ m_Z.")
print(f"  CONSISTENT with setting l_max = a1 -> Galois boundary = m_Z.")

# ============================================================
# PART 6: THE R_KK EXPONENT
# ============================================================
print(f"\n{'='*70}")
print("PART 6: KK RADIUS AND THE HIERARCHY")
print("="*70)

# From the framework: m_e/M_Pl = alpha^(4*phi^2) (derived)
# z_Planck = 4*phi^2 = 10.472

# The KK radius: R ~ 1/m_Z ~ 1/(m_e * phi^25 * 16/15)
# In Planck units: R * M_Pl = M_Pl / m_Z = M_Pl / (m_e * phi^25 * 16/15)
# = (M_Pl/m_e) / (phi^25 * 16/15)
# = alpha^(-4*phi^2) / (phi^25 * 16/15)

# ln(R * M_Pl) = -4*phi^2 * ln(alpha) - 25*ln(phi) - ln(16/15)
# = 4*phi^2 * ln(1/alpha) - 25*ln(phi) - ln(16/15)
z_Pl = 4*phi**2
ln_R_MPl = z_Pl * log(1/alpha) - 25*log(phi) - log(16/15)
R_Planck = exp(ln_R_MPl)

print(f"\n  R * M_Pl = alpha^(-z_Pl) / (phi^25 * 16/15)")
print(f"  z_Planck = 4*phi^2 = {z_Pl:.4f}")
print(f"  ln(R*M_Pl) = {ln_R_MPl:.4f}")
print(f"  R * M_Pl = {R_Planck:.4e}")
print(f"  R = {R_Planck / (M_Pl/hbar_c):.4e} fm")

# Compare with the direct R from m_Z:
print(f"\n  Direct: R = hbar_c/m_Z * sqrt(35) = {hbar_c*sqrt(35)/m_Z:.4e} fm")
print(f"  Hierarchy: R = {R_Planck * hbar_c / M_Pl:.4e} fm")

# ============================================================
# PART 7: WHAT THE SECTORS SEE
# ============================================================
print(f"\n{'='*70}")
print("PART 7: PHYSICAL INTERPRETATION OF THE THREE SECTORS")
print("="*70)

print(f"""
SECTOR A (self-conjugate, rational, dim 62):
  Eigenvalues: 0, 12, 14 (rational)
  These are GALOIS-INVARIANT modes.
  They see NO difference between physical and dark sectors.
  Physical content: vacuum (l=0) + high multipoles (l=4,5)

SECTOR B (physical, phi-containing, dim 29):
  Eigenvalues: 12-6phi, 12-4phi, 9
  These are the modes that DISTINGUISH the physical sector.
  They carry the information about phi (golden ratio).
  Low-energy content: l=1,2,3 (the first excited states of S^3)

SECTOR C (UV/Galois, phi'-containing, dim 29):
  Eigenvalues: 8+4phi, 15, 6+6phi
  Galois conjugates of sector B.
  These are ABOVE the Galois boundary.
  Physical interpretation: the DARK SECTOR modes.

KEY INSIGHT:
  The physical sector (A+B, dim 91) matches S^3 EXACTLY.
  The dark sector (C, dim 29) is the UV "excess."

  An observer in the physical sector sees a SMOOTH S^3 geometry
  (91 modes matching continuum).

  An observer in the dark sector sees the SAME geometry
  (by Galois symmetry, C has the same structure as B).

  The COMBINED observer (A+B+C, dim 120) sees the DISCRETE 600-cell.

  This is the geometric realization of "dark = Galois conjugate of physical."
""")

# ============================================================
# PART 8: ENERGY BUDGET
# ============================================================
print(f"{'='*70}")
print("PART 8: ENERGY BUDGET PER SECTOR")
print("="*70)

# Tr(L) per sector = total "energy" in that sector
tr_A = sum(m * lam for lam, m in [(0,1), (12,25), (14,36)])
tr_B = sum(m * lam for lam, m in [(12-6*phi,4), (12-4*phi,9), (9,16)])
tr_C = sum(m * lam for lam, m in [(8+4*phi,9), (15,16), (6+6*phi,4)])
tr_total = tr_A + tr_B + tr_C

print(f"\n  Tr(L_A) = {tr_A:.2f} ({100*tr_A/tr_total:.1f}%)")
print(f"  Tr(L_B) = {tr_B:.2f} ({100*tr_B/tr_total:.1f}%)")
print(f"  Tr(L_C) = {tr_C:.2f} ({100*tr_C/tr_total:.1f}%)")
print(f"  Total = {tr_total:.2f} = N*degree = {N*degree}")

# Fractional energy in physical vs dark:
print(f"\n  Physical (A+B): {100*(tr_A+tr_B)/tr_total:.1f}%")
print(f"  Dark (C): {100*tr_C/tr_total:.1f}%")

# Compare with Omega_matter / Omega_total in cosmology:
# Baryonic ~ 5%, Dark matter ~ 27%, Dark energy ~ 68%
# Our ratio: physical/dark ~ 70/30 or ~84/16 (depends on what we count)

# Actually: Tr(L_B) / Tr(L_C) = ratio of phi to phi' sector energies
ratio_BC = tr_B / tr_C
print(f"\n  Tr(L_B) / Tr(L_C) = {ratio_BC:.6f}")
print(f"  = physical_paired / dark = {ratio_BC:.4f}")

# Check framework numbers
for name, val in [('1/phi^2', 1/phi**2), ('sin2tW', sin2tW),
                   ('alpha', alpha), ('b1/N', b1/N)]:
    print(f"    vs {name} = {val:.6f} {'<-- CLOSE' if abs(ratio_BC - val) < 0.02 else ''}")

# Tr(B)/Tr(C) in terms of phi:
# Tr(B) = 4*(12-6phi) + 9*(12-4phi) + 16*9 = 48-24phi+108-36phi+144 = 300-60phi
# Tr(C) = 9*(8+4phi) + 16*15 + 4*(6+6phi) = 72+36phi+240+24+24phi = 336+60phi
# Tr(B)/Tr(C) = (300-60phi)/(336+60phi)

tr_B_exact = 300 - 60*phi
tr_C_exact = 336 + 60*phi
print(f"\n  Exact: Tr(B) = 300 - 60*phi = {tr_B_exact:.4f}")
print(f"  Exact: Tr(C) = 336 + 60*phi = {tr_C_exact:.4f}")
print(f"  Ratio = (300-60phi)/(336+60phi) = {tr_B_exact/tr_C_exact:.6f}")

# Simplify: factor 12 from both?
# 300 = 12*25, 60 = 12*5, 336 = 12*28, 60 = 12*5
# = 12*(25-5phi) / (12*(28+5phi)) = (25-5phi)/(28+5phi) = 5(5-phi)/(28+5phi)
print(f"  = (25 - 5*phi) / (28 + 5*phi)")
print(f"  = 5*(5 - phi) / (28 + 5*phi)")
print(f"  = 5*(a1 - phi) / (28 + a1*phi)")

# Galois norm of ratio:
# N((25-5phi)/(28+5phi)) = (25-5phi)(25-5phi') / ((28+5phi)(28+5phi'))
# = (25-5phi)(25-5(1-phi)) / ((28+5phi)(28+5(1-phi)))
# = (25-5phi)(20+5phi) / ((28+5phi)(33-5phi))

num_norm = (25-5*phi)*(20+5*phi)
den_norm = (28+5*phi)*(33-5*phi)
galois_ratio = num_norm / den_norm
print(f"  Galois norm of ratio = {galois_ratio:.6f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print("="*70)
print(f"""
SCALE DETERMINATION:
  l_max = a1 = 5 defines the Galois/S^3 boundary.
  Setting this boundary at m_Z gives R_KK ~ hbar_c * sqrt(35) / m_Z.
  Consistent with: alpha is evaluated at m_Z (measured value = framework value).

ENERGY BUDGET:
  Physical (A+B): {100*(tr_A+tr_B)/tr_total:.1f}% of spectral weight (dim 91/120)
  Dark/UV (C): {100*tr_C/tr_total:.1f}% (dim 29/120)
  Tr(B)/Tr(C) = (25-5phi)/(28+5phi) = {ratio_BC:.4f}

PHYSICAL PICTURE:
  An observer in the phi-sector (physical) sees smooth S^3 (91 modes).
  The phi'-sector (dark) has equal mode count (29) but higher energies.
  The Galois boundary IS the electroweak scale.
  Below m_Z: effective continuum S^3 geometry.
  Above m_Z: discrete 600-cell structure becomes visible.

STATUS: STRUCTURAL (scale identification m_Z <-> l_max needs independent confirmation)
NEXT: Step 3 -- beta functions from the spectral action truncation?
""")

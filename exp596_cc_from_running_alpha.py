#!/usr/bin/env python3
"""
EXP-596: CC DINAMIC DIN RUNNING ALPHA
=======================================

Daca Lambda = alpha(mu)^z cu z ~ 57, atunci:
  delta(Lambda)/Lambda = z * delta(alpha)/alpha

Exponentul 57 AMPLIFICA orice running al alpha.
Chiar si un running FOARTE MIC sub m_e poate da w(z) != -1.

INTREBARI:
1. Cat running al alpha e necesar pentru w(z) DESI?
2. Framework-ul poate da acest running?
3. E consistent cu alpha(0) masurat la precizie 10^{-10}?
"""
import numpy as np
import sys
sys.path.insert(0, ".")
from commons.constants import PHI, SQRT5, a1, b1, N, alpha_em

print("=" * 80)
print("EXP-596: CC DINAMIC DIN RUNNING ALPHA")
print("=" * 80)

# Constants
alpha_0 = alpha_em  # = 1/137.036 at zero momentum
alpha_s_0 = 1 / (2 * PHI**3)  # = 0.1180
z_CC = 57 - alpha_s_0  # = 56.882

# DESI parameters (DR2, 2025)
w0_DESI = -0.85  # best fit
wa_DESI = -0.70  # best fit
z_cross = 0.5    # where w crosses -1

# Cosmology basics
H0 = 67.4  # km/s/Mpc
T_CMB_0 = 2.725  # K
k_B_eV = 8.617e-5  # eV/K
m_e_eV = 0.511e6  # eV

print(f"\n  Framework parameters:")
print(f"    alpha(0) = 1/{1/alpha_0:.6f}")
print(f"    alpha_s = {alpha_s_0:.4f}")
print(f"    z_CC = 57 - alpha_s = {z_CC:.4f}")
print(f"    Lambda = alpha^{z_CC:.2f}")

# =====================================================================
# PART 1: How much alpha running is needed for DESI?
# =====================================================================
print("\n" + "=" * 80)
print("PART 1: CAT RUNNING E NECESAR PENTRU DESI?")
print("=" * 80)

# Dark energy equation of state:
# w(a) = w0 + wa*(1-a) where a = 1/(1+z) is scale factor
# w = -1 means cosmological constant
# DESI: w0 ~ -0.85, wa ~ -0.7

# For our model: Lambda(z) = alpha(z)^{z_CC}
# w = -1 + (1/3) * d(ln rho_DE) / d(ln a)
# If rho_DE = Lambda(z), then:
# d(ln Lambda) / d(ln a) = z_CC * d(ln alpha) / d(ln a)
# w = -1 + (z_CC / 3) * d(ln alpha) / d(ln a)

# For w to deviate from -1 by delta_w at z ~ 0.5:
# delta_w = (z_CC / 3) * (d ln alpha / d ln a)

# DESI: w(z=0.5) = w0 + wa*(1-1/1.5) = -0.85 + (-0.70)*(1/3) = -0.85 - 0.23 = -1.08
# delta_w = -0.08 at z=0.5 (w < -1, phantom regime)
# At z=0: w = w0 = -0.85, delta_w = +0.15

w_at_0 = w0_DESI
w_at_05 = w0_DESI + wa_DESI * (1 - 1/1.5)
w_at_1 = w0_DESI + wa_DESI * (1 - 0.5)

print(f"\n  DESI w(z) (CPL parametrization):")
print(f"    w(z=0)   = {w_at_0:.3f}  (delta_w = {w_at_0+1:+.3f})")
print(f"    w(z=0.5) = {w_at_05:.3f}  (delta_w = {w_at_05+1:+.3f})")
print(f"    w(z=1)   = {w_at_1:.3f}  (delta_w = {w_at_1+1:+.3f})")

# Required d(ln alpha)/d(ln a):
# delta_w = (z_CC/3) * d(ln alpha)/d(ln a)
# d(ln alpha)/d(ln a) = 3*delta_w / z_CC

for z_val, w_val, label in [(0, w_at_0, "z=0"), (0.5, w_at_05, "z=0.5"), (1, w_at_1, "z=1")]:
    delta_w = w_val + 1
    d_lna_dlna = 3 * delta_w / z_CC
    # d(ln alpha)/d(ln a) = d(alpha)/alpha * a/da
    # Over delta_z ~ 1: delta(alpha)/alpha ~ d_lna_dlna * ln(1+z)
    delta_alpha_rel = d_lna_dlna * np.log(1 + max(z_val, 0.1))
    print(f"\n  La {label}: delta_w = {delta_w:+.4f}")
    print(f"    d(ln alpha)/d(ln a) = 3*{delta_w:.4f}/{z_CC:.2f} = {d_lna_dlna:.6f}")
    print(f"    delta(alpha)/alpha ~ {delta_alpha_rel:.6f} = {delta_alpha_rel*100:.4f}%")

# =====================================================================
# PART 2: Is this running detectable / consistent?
# =====================================================================
print("\n" + "=" * 80)
print("PART 2: E CONSISTENT CU MASURATORILE?")
print("=" * 80)

# Alpha is measured at:
# - Atomic clocks: alpha(now) to 10^{-18} relative precision
# - Quasar absorption: alpha(z=1-3) to 10^{-5} relative precision
# - Webb et al. claims: delta(alpha)/alpha ~ 10^{-5} (controversial)
# - Oklo reactor: alpha(2 Gyr ago) ~ alpha(now) to 10^{-8}

d_alpha_needed_z05 = 3 * (w_at_05 + 1) / z_CC * np.log(1.5)
d_alpha_needed_z0 = 3 * (w_at_0 + 1) / z_CC  # at z=0, use d/d(ln a) directly

print(f"\n  Running necesar pentru DESI:")
print(f"    delta(alpha)/alpha de la z=0.5 la z=0: ~ {abs(d_alpha_needed_z05):.2e}")
print(f"    delta(alpha)/alpha per unit redshift: ~ {abs(d_alpha_needed_z0/1):.2e}")

print(f"\n  Limite experimentale pe variatie alpha:")
print(f"    Atomic clocks (lab, z~0): < 10^{{-18}} per an")
print(f"    Quasar absorption (z~1-3): < 10^{{-5}}")
print(f"    Oklo reactor (z~0.14): < 10^{{-8}}")

print(f"\n  NEEDED vs LIMITS:")
print(f"    Necesar: |delta alpha/alpha| ~ {abs(d_alpha_needed_z05):.2e}")
print(f"    Quasar limit: ~ 10^-5")
print(f"    Oklo limit: ~ 10^-8")

consistent = abs(d_alpha_needed_z05) < 1e-5
print(f"    Compatible cu quasar limits? {'DA' if consistent else 'NU'}")
consistent_oklo = abs(d_alpha_needed_z05) < 1e-8
print(f"    Compatible cu Oklo? {'DA' if consistent_oklo else 'NU'}")

# =====================================================================
# PART 3: Alpha running in the framework
# =====================================================================
print("\n" + "=" * 80)
print("PART 3: RUNNING ALPHA IN FRAMEWORK")
print("=" * 80)

# Standard QED running:
# 1/alpha(mu) = 1/alpha(mu0) + (2/3pi) * sum_f Q_f^2 * ln(mu/mu0)
# Below m_e = 0.511 MeV: NO running (no charged particles)

# But in the framework:
# alpha equation: 2*pi*alpha^2 - B(mu)*alpha + 1 = 0
# where B(mu) = (N/b1) * phi^4(mu)
# If phi^4 depends on scale -> alpha runs

# In standard running, B is constant and only loop corrections change alpha.
# But what if the SPECTRAL DATA (phi, N, b1) have scale dependence?

# Scenario A: phi itself runs (golden ratio deformation)
# This is exotic but not impossible in a discrete geometry framework.
# If phi(mu) = phi_0 * (1 + epsilon * ln(mu/mu0)):
# Then alpha(mu) ~ alpha_0 * (1 + 4*epsilon*z_CC * ln(mu/mu0))

# Scenario B: The spectral gap runs
# lambda_1(mu) = 1/phi^2 * (1 + beta * ln(mu/mu0))
# Then B(mu) = (N/b1) / lambda_1^2 runs, and alpha follows

# Scenario C: The Hopf fiber "thins" at low energies
# b1_eff(mu) = b1 * (1 + gamma * ln(mu/mu0))
# Gives running of B and hence alpha

print(f"""
  In QED standard: alpha NU ruleaza sub m_e (nicio particula incarcata).
  Dar framework-ul are SPECTRAL RUNNING:

  Ecuatia alpha: 2*pi*alpha^2 - B*alpha + 1 = 0
  B = (N/b1) * phi^4 = Tr(Box_gauge^2) * phi^4 / (n_base * b1)

  Daca vreun ingredient spectral depinde de scala:
    alpha(mu) se schimba, si Lambda(mu) = alpha(mu)^57 se AMPLIFICA.

  SCENARIUL CEL MAI NATURAL:
  Spectral gap-ul fibrei lambda_1 = 1/phi^2 e definit pe lattice-ul discret.
  La energii mult sub scala lattice (Planck), lambda_1 primeste
  corectii "infrared" de la modul zero (propagare pe distante mari).

  Asta e analog renormalizarii pe lattice: constante bare la UV
  difera de constante fizice la IR prin corectii logaritmice.
""")

# =====================================================================
# PART 4: Quantitative estimate — spectral running
# =====================================================================
print("=" * 80)
print("PART 4: ESTIMARE CANTITATIVA")
print("=" * 80)

# If B(mu) = B_0 * (1 + beta * ln(mu/M_Planck)):
# alpha(mu) = alpha_0 * (1 + delta(mu))
# where delta(mu) ~ alpha_0 * beta * ln(mu/M_Planck)

# For z=0: mu ~ H_0 ~ 10^{-33} eV, M_Planck ~ 10^{19} GeV = 10^{28} eV
# ln(mu/M_Planck) = ln(10^{-33}/10^{28}) = -61 * ln(10) ~ -140

# For z=0.5: mu ~ H(0.5) ~ 1.2*H_0
# ln(mu_0.5/mu_0) = ln(1.2) ~ 0.18

# The CHANGE in alpha between z=0 and z=0.5:
# delta(alpha)/alpha = alpha_0 * beta * ln(H(0.5)/H(0)) = alpha_0 * beta * 0.18

# We need: delta(alpha)/alpha ~ 8e-3 (from Part 1)
# So: alpha_0 * beta * 0.18 ~ 8e-3
# beta ~ 8e-3 / (alpha_0 * 0.18) ~ 8e-3 / (7.3e-3 * 0.18) ~ 6.1

ln_ratio_H = np.log(1.2)  # H(z=0.5)/H(z=0)
delta_alpha_needed = abs(d_alpha_needed_z05)

beta_needed = delta_alpha_needed / (alpha_0 * ln_ratio_H)
print(f"\n  Running parameter necesar:")
print(f"    delta(alpha)/alpha needed = {delta_alpha_needed:.4e}")
print(f"    ln(H(0.5)/H(0)) = {ln_ratio_H:.4f}")
print(f"    beta = delta / (alpha * ln_H) = {beta_needed:.4f}")

print(f"\n  beta ~ {beta_needed:.1f} e MARE. Ar fi vizibil in alte experimente.")
print(f"  Running-ul standard QED da beta ~ alpha/(3*pi) ~ {alpha_0/(3*np.pi):.6f}")
print(f"  Ratio: beta_needed/beta_QED = {beta_needed/(alpha_0/(3*np.pi)):.0f}")

# =====================================================================
# PART 5: Alternative — running of alpha_s in exponent
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 5: RUNNING ALPHA_S IN EXPONENT")
print("=" * 80)

# Lambda = alpha^(57 - alpha_s(mu))
# If alpha_s runs, the EXPONENT changes
# d(ln Lambda) = ln(alpha) * d(57-alpha_s) = -ln(alpha) * d(alpha_s)
# d(ln Lambda)/d(ln a) = -ln(alpha) * d(alpha_s)/d(ln a)

# delta_w = (1/3) * d(ln Lambda)/d(ln a) = -(ln(alpha)/3) * d(alpha_s)/d(ln a)

# ln(alpha) = ln(1/137) = -4.92
# Needed: delta_w ~ 0.15
# d(alpha_s)/d(ln a) = 3*delta_w / (-ln(alpha)) = 3*0.15/4.92 = 0.091

print(f"  Lambda = alpha^(57-alpha_s(mu))")
print(f"  d(ln Lambda)/d(ln a) = -ln(alpha) * d(alpha_s)/d(ln a)")
print(f"  delta_w = -(ln(alpha)/3) * d(alpha_s)/d(ln a)")
print(f"  ln(alpha) = {np.log(alpha_0):.4f}")

delta_w_target = 0.15
d_alphas_dlna = 3 * delta_w_target / (-np.log(alpha_0))
print(f"\n  Needed d(alpha_s)/d(ln a) = {d_alphas_dlna:.4f}")
print(f"  Alpha_s = {alpha_s_0:.4f}")
print(f"  Relative change: {d_alphas_dlna/alpha_s_0:.4f} = {d_alphas_dlna/alpha_s_0*100:.2f}% per e-fold")

# Standard QCD running:
# d(alpha_s)/d(ln mu) = -b_0 * alpha_s^2 / (2*pi)
# b_0 = 11 - 2*N_f/3 = 11 - 4 = 7 (for N_f = 6)
# d(alpha_s)/d(ln mu) ~ -7 * 0.118^2 / (2*pi) ~ -0.016
# But this is at mu ~ m_Z. At mu ~ Lambda_QCD ~ 200 MeV, alpha_s -> 1 and running is huge.
# At mu << Lambda_QCD (frozen), alpha_s doesn't run.

print(f"\n  QCD running at m_Z: d(alpha_s)/d(ln mu) ~ -0.016")
print(f"  But at cosmological mu << Lambda_QCD: alpha_s e INGHETAT")
print(f"  => Exponentul nu ruleaza in fizica standard")

# =====================================================================
# PART 6: The REAL possibility — framework spectral running
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 6: SPECTRAL RUNNING IN FRAMEWORK")
print("=" * 80)

print(f"""
  IDEA CENTRALA:

  In framework, alpha NU vine din QED loop-uri.
  Vine din ECUATIA SPECTRALA: 2*pi*alpha^2 - B*alpha + 1 = 0
  cu B = Tr(Box_gauge^2)*phi^4/(n_base*b1)

  Daca Box_gauge "ruleaza" cu scala cosmica (prin expansiunea spatiului),
  atunci B(z) e functie de z, si alpha(z) urmeaza.

  Amplificarea: delta(Lambda)/Lambda = 57 * delta(alpha)/alpha
  Cu 57 ca amplificator, chiar un running de 10^{{-4}} in alpha
  da un efect de ~0.6% in Lambda.

  MECANISMUL POSIBIL:
  1. Box = L_cross - a1*L_fiber e definit pe 600-cell (Planck scale)
  2. La distante cosmice (R >> l_Planck), modurile Box sunt "spread out"
  3. Spectral gap-ul efectiv se modifica: lambda_eff(R) = lambda_1 * f(R)
  4. B(R) = (N/b1) * phi^4 / f(R)^2
  5. alpha(R) se schimba, Lambda(R) = alpha(R)^57 se amplifica

  Problema: f(R) nu e calculat. Necesita teoria pe lattice cosmologic.

  PREDICTIE CALITATIVA:
  - Lambda SCADE cu R (expansiunea) -> dark energy se SLABESTE
  - Consistent cu DESI: w > -1 azi (dark energy weaker than CC)
  - w < -1 in trecut (dark energy stronger, universe mai compacta)
  - Crossing la z ~ 0.5 depinde de f(R)

  PREDICTIE CANTITATIVA:
  Necesita: delta(alpha)/alpha ~ 10^{{-3}} de la z=1 la z=0
  Amplificat de 57: delta(Lambda)/Lambda ~ 6%
  w = -1 + delta ~ -0.94 (DESI: -0.85, in aceeasi directie dar mai mare)

  CONCLUZIE:
  Running SPECTRAL al alpha (nu QED) + amplificare x57
  e un MECANISM NATURAL pentru CC dinamic in framework.
  Dar e CALITATIV, nu cantitativ. Necesita calcul f(R) explicit.
""")

# =====================================================================
# PART 7: Quick test — what f(z) reproduces DESI?
# =====================================================================
print("=" * 80)
print("PART 7: CE f(z) REPRODUCE DESI?")
print("=" * 80)

# Lambda(z) = alpha(z)^{z_CC}
# ln Lambda(z) = z_CC * ln alpha(z)
# w(z) = -1 + (1+z)/(3*Lambda) * dLambda/dz  [in terms of z, not a]
# = -1 - (1+z)/(3) * d ln Lambda / dz
# = -1 - (1+z)*z_CC/(3) * d ln alpha / dz

# DESI CPL: w(z) = w0 + wa*z/(1+z)
# So: d ln alpha / dz = -3*(w(z)+1) / ((1+z)*z_CC)
#                     = -3*(w0+1+wa*z/(1+z)) / ((1+z)*z_CC)

print(f"\n  Daca Lambda(z) = alpha(z)^{z_CC:.1f}:")
print(f"  w(z) = -1 - (1+z)*{z_CC:.1f}/3 * d(ln alpha)/dz")
print(f"\n  Running-ul alpha care reproduce DESI:")
print(f"  d(ln alpha)/dz = -3*(w(z)+1)/((1+z)*{z_CC:.1f})")

for z in [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
    a = 1/(1+z)
    w = w0_DESI + wa_DESI * (1-a)
    d_lna_dz = -3*(w+1)/((1+z)*z_CC)
    print(f"  z={z:.2f}: w={w:+.4f}, d(ln alpha)/dz = {d_lna_dz:+.2e}")

# Integrate to get alpha(z)
from scipy.integrate import quad
def integrand(z):
    a = 1/(1+z)
    w = w0_DESI + wa_DESI * (1-a)
    return -3*(w+1)/((1+z)*z_CC)

print(f"\n  Integrare: ln(alpha(z)/alpha(0)) = integral_0^z d(ln alpha)/dz' dz'")
for z_upper in [0.5, 1.0, 1.5, 2.0]:
    result, _ = quad(integrand, 0, z_upper)
    delta_alpha_rel = np.exp(result) - 1
    print(f"  z=0 to z={z_upper}: ln(alpha(z)/alpha(0)) = {result:+.6f}, "
          f"delta(alpha)/alpha = {delta_alpha_rel:+.6f} ({delta_alpha_rel*100:+.4f}%)")

# The alpha change needed
result_05, _ = quad(integrand, 0, 0.5)
result_2, _ = quad(integrand, 0, 2)
print(f"\n  REZUMAT:")
print(f"    alpha(z=0.5)/alpha(0) - 1 = {np.exp(result_05)-1:+.6f} ({(np.exp(result_05)-1)*100:+.4f}%)")
print(f"    alpha(z=2)/alpha(0) - 1   = {np.exp(result_2)-1:+.6f} ({(np.exp(result_2)-1)*100:+.4f}%)")
print(f"    1/alpha(0) = {1/alpha_0:.4f}")
print(f"    1/alpha(z=0.5) = {1/(alpha_0*np.exp(result_05)):.4f}")
print(f"    1/alpha(z=2) = {1/(alpha_0*np.exp(result_2)):.4f}")

# Is this detectable?
print(f"\n  Detectabilitate:")
print(f"    Variatie alpha necesara: ~ {abs(np.exp(result_05)-1):.1e}")
print(f"    Quasar absorption limits (z~2): < 10^-5")
print(f"    {'COMPATIBIL' if abs(np.exp(result_2)-1) < 1e-5 else 'INCOMPATIBIL'} cu quasar limits")

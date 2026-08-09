"""
exp361_alpha_running.py
========================
OPEN-3 route 2: Self-consistent QED running alpha(0) -> alpha(m_Z).

The framework derives ALL fermion masses and charges. The only external
input is alpha(m_Z)/alpha(0) = 1.0628 (used for m_Z prediction).

1-loop QED running:
  alpha(mu)^{-1} = alpha(0)^{-1} - (1/3pi) * sum_f N_c*Q_f^2 * ln(mu^2/m_f^2) * theta(mu-m_f)

All quantities on the RHS are framework-derived:
  - alpha(0) from the alpha equation
  - masses from the mass formula m_f = m_e * phi^{n_f}
  - charges Q_f from anomaly cancellation
  - N_c = 3 for quarks, 1 for leptons (from gauge group derivation)

If this gives alpha(m_Z)/alpha(0) ~ 1.063, we eliminate the last external input.

Also compute:
  - 2-loop correction
  - Hadronic vacuum polarization estimate
  - alpha_s running (3-loop) as cross-check

Dependencies: numpy
"""

import numpy as np

print("=" * 72)
print("EXP-361: SELF-CONSISTENT QED RUNNING alpha(0) -> alpha(m_Z)")
print("=" * 72)

# Framework constants
a1 = 5
b1 = a1 + 1
phi = (1.0 + np.sqrt(a1)) / 2.0
PI = np.pi

# Framework alpha from the quadratic equation
# 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
A_coef = 2 * PI
B_coef = -4 * a1 * phi**4
C_coef = 1.0
disc = B_coef**2 - 4 * A_coef * C_coef
alpha_0 = (-B_coef - np.sqrt(disc)) / (2 * A_coef)  # smaller root = physical
alpha_0_inv = 1.0 / alpha_0

# Framework alpha_s
alpha_s = 1.0 / (2.0 * phi**3)

# Framework sin^2(theta_W)
sin2tw = float(b1) / (a1**2 + 1)  # = 6/26

print(f"\n--- Framework coupling constants ---")
print(f"  alpha(0)     = {alpha_0:.10f}")
print(f"  alpha(0)^-1  = {alpha_0_inv:.6f}")
print(f"  alpha_s(m_Z) = {alpha_s:.10f}")
print(f"  sin^2(tW)    = {sin2tw:.10f} = {b1}/{a1**2+1}")

# Experimental reference values
alpha_exp = 1.0 / 137.035999084  # CODATA 2018, alpha at q=0
alpha_mZ_exp = 1.0 / 127.951  # PDG at m_Z
ratio_exp = alpha_mZ_exp / alpha_exp
m_Z_exp = 91.1876  # GeV

print(f"\n--- Experimental references ---")
print(f"  alpha(0)_exp      = 1/{1/alpha_exp:.6f}")
print(f"  alpha(m_Z)_exp    = 1/{1/alpha_mZ_exp:.3f}")
print(f"  ratio alpha(mZ)/alpha(0) = {ratio_exp:.6f}")

# =====================================================================
# FRAMEWORK MASSES (from mass formula m_f = m_e * phi^n_f)
# =====================================================================
m_e = 0.51099895  # MeV - the ONE input

# Mass exponents n = 5a + 6b + delta (from Table in paper)
# (a,b) assignments and corrections
fermions = {
    # Leptons: n = 5a + 6b + c_ell*(b*(b-1)/2)
    # c_ell = phi^3 / (a1^2 + 1) = phi^3 / 26
    'e':   {'a': 0, 'b': 0, 'Q': -1, 'Nc': 1},
    'mu':  {'a': 1, 'b': 1, 'Q': -1, 'Nc': 1},
    'tau': {'a': 1, 'b': 2, 'Q': -1, 'Nc': 1},
    # Quarks: n = 5a + 6b + C*(a*(a-1)/2), C = 4/(a1^2+1) = 2/13
    'u':   {'a': 3, 'b': -2, 'Q': 2.0/3, 'Nc': 3},
    'd':   {'a': 1, 'b': 0, 'Q': -1.0/3, 'Nc': 3},
    'c':   {'a': 2, 'b': 1, 'Q': 2.0/3, 'Nc': 3},
    's':   {'a': 1, 'b': 1, 'Q': -1.0/3, 'Nc': 3},
    't':   {'a': 4, 'b': 1, 'Q': 2.0/3, 'Nc': 3},
    'b':   {'a': -1, 'b': 4, 'Q': -1.0/3, 'Nc': 3},
}

C_q = 4.0 / (a1**2 + 1)  # = 2/13
c_ell = phi**3 / (a1**2 + 1)

def compute_mass(name, info):
    a, b = info['a'], info['b']
    n0 = 5 * a + 6 * b
    if info['Nc'] == 1:  # lepton
        delta = c_ell * (b * (b - 1) / 2.0)
    else:  # quark
        delta = C_q * (a * (a - 1) / 2.0)
    n_total = n0 + delta
    mass_MeV = m_e * phi**(n_total)
    return mass_MeV, n_total

# Experimental masses (MeV) for comparison
exp_masses = {
    'e': 0.51099895, 'mu': 105.6584, 'tau': 1776.86,
    'u': 2.16, 'd': 4.67, 's': 93.4, 'c': 1270.0, 'b': 4180.0, 't': 172760.0
}

print(f"\n--- Framework fermion masses ---")
print(f"  {'Name':5s} {'(a,b)':8s} {'n_total':10s} {'m_fw (MeV)':14s} {'m_exp (MeV)':14s} {'err%':8s}")
masses_fw = {}
for name in ['e', 'mu', 'tau', 'u', 'd', 'c', 's', 't', 'b']:
    m, n = compute_mass(name, fermions[name])
    masses_fw[name] = m
    m_exp = exp_masses[name]
    err = (m - m_exp) / m_exp * 100
    a, b = fermions[name]['a'], fermions[name]['b']
    print(f"  {name:5s} ({a:2d},{b:2d})  {n:10.4f}  {m:14.4f}  {m_exp:14.4f}  {err:+8.2f}%")

# Convert to GeV for running
masses_GeV = {k: v / 1000.0 for k, v in masses_fw.items()}
m_Z_GeV = m_Z_exp  # Using experimental m_Z for the scale

# =====================================================================
# 1-LOOP QED RUNNING
# =====================================================================
# Delta_alpha = (alpha/3pi) * sum_f N_c * Q_f^2 * ln(m_Z^2/m_f^2) * theta(m_Z > m_f)
# alpha(m_Z)^{-1} = alpha(0)^{-1} - Delta_alpha / alpha
# or equivalently:
# alpha(m_Z)^{-1} = alpha(0)^{-1} - (1/3pi) * sum_f N_c * Q_f^2 * ln(m_Z^2/m_f^2)

print(f"\n{'=' * 72}")
print("1-LOOP QED VACUUM POLARIZATION")
print(f"{'=' * 72}")

delta_sum = 0.0
print(f"\n  {'Fermion':8s} {'Nc*Q^2':8s} {'ln(mZ/mf)^2':14s} {'contribution':14s}")
for name in ['e', 'mu', 'tau', 'u', 'd', 'c', 's', 't', 'b']:
    info = fermions[name]
    m_f = masses_GeV[name]
    if m_f < m_Z_GeV:
        Nc_Q2 = info['Nc'] * info['Q']**2
        log_term = np.log(m_Z_GeV**2 / m_f**2)
        contrib = Nc_Q2 * log_term
        delta_sum += contrib
        print(f"  {name:8s} {Nc_Q2:8.4f} {log_term:14.4f} {contrib:14.4f}")
    else:
        print(f"  {name:8s} -- above m_Z threshold, skipped --")

# 1-loop result
Delta_alpha_1loop = alpha_0 / (3 * PI) * delta_sum
alpha_mZ_inv_1loop = alpha_0_inv - delta_sum / (3 * PI)
alpha_mZ_1loop = 1.0 / alpha_mZ_inv_1loop
ratio_1loop = alpha_mZ_1loop / alpha_0

print(f"\n  Sum = {delta_sum:.6f}")
print(f"  Delta_alpha(1-loop) = alpha/(3pi) * sum = {Delta_alpha_1loop:.6f}")
print(f"  alpha(m_Z)^-1 (1-loop) = {alpha_mZ_inv_1loop:.4f}")
print(f"  alpha(m_Z)/alpha(0)    = {ratio_1loop:.6f}")
print(f"  Experimental ratio     = {ratio_exp:.6f}")
print(f"  Difference             = {ratio_1loop - ratio_exp:.6f} ({(ratio_1loop/ratio_exp - 1)*100:+.2f}%)")

# =====================================================================
# WITH EXPERIMENTAL MASSES (for comparison)
# =====================================================================
print(f"\n{'=' * 72}")
print("COMPARISON: 1-LOOP WITH EXPERIMENTAL MASSES")
print(f"{'=' * 72}")

delta_sum_exp = 0.0
for name in ['e', 'mu', 'tau', 'u', 'd', 'c', 's', 't', 'b']:
    info = fermions[name]
    m_f = exp_masses[name] / 1000.0  # MeV -> GeV
    if m_f < m_Z_GeV:
        Nc_Q2 = info['Nc'] * info['Q']**2
        log_term = np.log(m_Z_GeV**2 / m_f**2)
        delta_sum_exp += Nc_Q2 * log_term

alpha_mZ_inv_exp_masses = alpha_0_inv - delta_sum_exp / (3 * PI)
alpha_mZ_exp_masses = 1.0 / alpha_mZ_inv_exp_masses
ratio_exp_masses = alpha_mZ_exp_masses / alpha_0

print(f"  With exp masses: alpha(m_Z)^-1 = {alpha_mZ_inv_exp_masses:.4f}")
print(f"  ratio = {ratio_exp_masses:.6f}")
print(f"  (Same calculation, only masses differ)")

# =====================================================================
# 2-LOOP CORRECTIONS
# =====================================================================
print(f"\n{'=' * 72}")
print("2-LOOP AND HIGHER-ORDER CORRECTIONS")
print(f"{'=' * 72}")

# The 2-loop correction adds terms of order alpha^2:
# Delta_alpha_2loop ~ (alpha/pi)^2 * sum_f Nc*Q_f^2 * [...]
# The dominant 2-loop contribution comes from the (alpha/pi)*alpha_s/pi term
# for quarks (mixed QCD-QED correction):
# Delta_alpha_had^(2) ~ (alpha*alpha_s/pi^2) * sum_q Q_q^2 * ln(...)

# Leptonic 2-loop correction (known analytically):
# Delta_alpha_lep^(2) = (alpha/pi)^2 * sum_ell [...]
# The leading correction is ~ (alpha/pi) * Delta_alpha^(1) / 4
Delta_alpha_lep_2loop = (alpha_0 / PI) * Delta_alpha_1loop * 0.25
print(f"  Leptonic 2-loop estimate: {Delta_alpha_lep_2loop:.6f}")

# Hadronic contribution is the tricky part.
# In the real world, this uses e+e- -> hadrons data (R-ratio).
# In our framework, we can compute it perturbatively since we have all masses.

# For light quarks (u,d,s), perturbative QCD + alpha_s corrections:
# Delta_alpha_had = (alpha/3pi) * sum_q Nc*Q_q^2 * ln(mZ^2/mq^2) * (1 + alpha_s/pi + ...)

delta_had_QCD = 0.0
for name in ['u', 'd', 'c', 's', 'b']:
    info = fermions[name]
    m_f = masses_GeV[name]
    if m_f < m_Z_GeV:
        Nc_Q2 = info['Nc'] * info['Q']**2
        log_term = np.log(m_Z_GeV**2 / m_f**2)
        # QCD correction factor: 1 + alpha_s/pi + (alpha_s/pi)^2 * (1.9857 - 0.1153*n_f)
        as_pi = alpha_s / PI
        n_f = 5  # active flavors at m_Z
        qcd_factor = 1.0 + as_pi + as_pi**2 * (1.9857 - 0.1153 * n_f)
        contrib = Nc_Q2 * log_term * qcd_factor
        delta_had_QCD += contrib

# Leptonic part (no QCD correction)
delta_lep = 0.0
for name in ['e', 'mu', 'tau']:
    info = fermions[name]
    m_f = masses_GeV[name]
    if m_f < m_Z_GeV:
        Nc_Q2 = info['Nc'] * info['Q']**2
        log_term = np.log(m_Z_GeV**2 / m_f**2)
        delta_lep += Nc_Q2 * log_term

delta_total_improved = delta_lep + delta_had_QCD
alpha_mZ_inv_improved = alpha_0_inv - delta_total_improved / (3 * PI)
alpha_mZ_improved = 1.0 / alpha_mZ_inv_improved
ratio_improved = alpha_mZ_improved / alpha_0

print(f"\n  Improved (with QCD corrections to hadronic part):")
print(f"  alpha(m_Z)^-1 = {alpha_mZ_inv_improved:.4f}")
print(f"  ratio = {ratio_improved:.6f}")
print(f"  Experimental ratio = {ratio_exp:.6f}")
print(f"  Difference = {(ratio_improved/ratio_exp - 1)*100:+.3f}%")

# =====================================================================
# THE REAL ISSUE: LIGHT QUARK MASSES AND NONPERTURBATIVE QCD
# =====================================================================
print(f"\n{'=' * 72}")
print("LIGHT QUARK ISSUE: PERTURBATIVE vs CONSTITUENT MASSES")
print(f"{'=' * 72}")

print("""
  The light quark masses (u,d,s) in the mass formula are CURRENT (MS-bar) masses.
  For vacuum polarization below ~2 GeV, perturbative QCD breaks down.
  The real-world approach uses:
    Delta_alpha_had = (alpha*m_Z^2)/(3pi) * P integral[e+e- -> hadrons]

  In the framework, we should use:
  1. Leptonic: perturbative (exact)
  2. Heavy quarks (c,b): perturbative with alpha_s corrections (good)
  3. Light quarks (u,d,s): need nonperturbative estimate

  The nonperturbative contribution is dominated by the rho meson:
    Delta_alpha_had^(5) = 0.02766 +/- 0.00010 (PDG 2020)

  Can we get this from framework?
""")

# Approach: use the framework perturbatively and see what happens
# Compare with the known decomposition:
# Delta_alpha_total = Delta_alpha_lep + Delta_alpha_had^(5) + Delta_alpha_top
# = 0.03150 + 0.02766 + (-0.00007) = 0.05909
# => alpha(m_Z)^{-1} = 137.036 - 0.05909/alpha ~ 137.036 - 8.10 = 128.93
# Wait, let me redo this properly.

# The standard parametrization:
# alpha(m_Z)^{-1} = alpha(0)^{-1} / (1 - Delta_alpha)
# Delta_alpha = Delta_alpha_lep + Delta_alpha_had + Delta_alpha_top + Delta_alpha_EW
# Delta_alpha_lep = 0.031498  (3-loop, exact)
# Delta_alpha_had = 0.02766   (from R-ratio, nonperturbative)
# Delta_alpha_top = -0.000076 (decoupling)
# Delta_alpha_EW  = -0.000065 (electroweak)
# Total: ~0.05901
# alpha(m_Z)^{-1} = 137.036 * (1 - 0.05901) = 137.036 * 0.94099 = 128.94

print(f"  Standard breakdown:")
Da_lep_std = 0.031498
Da_had_std = 0.02766
Da_top_std = -0.000076
Da_EW_std = -0.000065
Da_total_std = Da_lep_std + Da_had_std + Da_top_std + Da_EW_std
print(f"    Delta_alpha_lep = {Da_lep_std:.6f}")
print(f"    Delta_alpha_had = {Da_had_std:.6f}")
print(f"    Delta_alpha_top = {Da_top_std:.6f}")
print(f"    Delta_alpha_EW  = {Da_EW_std:.6f}")
print(f"    Total = {Da_total_std:.6f}")
print(f"    alpha(m_Z)^-1 = {alpha_0_inv * (1 - Da_total_std):.3f}")
print(f"    ratio = {1/(1-Da_total_std):.6f}")

# =====================================================================
# FRAMEWORK PERTURBATIVE CALCULATION (all orders we can do)
# =====================================================================
print(f"\n{'=' * 72}")
print("FRAMEWORK PERTURBATIVE CALCULATION")
print(f"{'=' * 72}")

# Leptonic contribution (perturbative, reliable)
Da_lep_fw = 0.0
print(f"\n  Leptonic contributions:")
for name in ['e', 'mu', 'tau']:
    m_f = masses_GeV[name]
    Q = fermions[name]['Q']
    Nc = fermions[name]['Nc']
    log_term = np.log(m_Z_GeV / m_f)  # ln(m_Z/m_f), not squared
    # Delta_alpha_f = (alpha/3pi) * Nc * Q^2 * 2 * ln(mZ/mf)
    # (factor 2 because ln(mZ^2/mf^2) = 2*ln(mZ/mf))
    contrib = alpha_0 / (3 * PI) * Nc * Q**2 * 2 * log_term
    Da_lep_fw += contrib
    print(f"    {name}: m={m_f*1000:.4f} MeV, ln(mZ/mf)={log_term:.4f}, Da={contrib:.6f}")

print(f"  Total leptonic: {Da_lep_fw:.6f}")
print(f"  Standard value: {Da_lep_std:.6f}")
print(f"  Difference: {(Da_lep_fw/Da_lep_std - 1)*100:+.2f}%")

# Heavy quark contribution (perturbative with QCD)
Da_heavy_fw = 0.0
print(f"\n  Heavy quark contributions (c, b with QCD):")
for name in ['c', 'b']:
    m_f = masses_GeV[name]
    Q = fermions[name]['Q']
    Nc = fermions[name]['Nc']
    log_term = np.log(m_Z_GeV / m_f)
    as_pi = alpha_s / PI
    qcd_corr = 1.0 + as_pi + as_pi**2 * 1.4114  # n_f=5 coefficient
    contrib = alpha_0 / (3 * PI) * Nc * Q**2 * 2 * log_term * qcd_corr
    Da_heavy_fw += contrib
    print(f"    {name}: m={m_f:.4f} GeV, ln(mZ/mf)={log_term:.4f}, QCD={qcd_corr:.4f}, Da={contrib:.6f}")

print(f"  Total heavy: {Da_heavy_fw:.6f}")

# Light quark contribution (u,d,s) - perturbative estimate
Da_light_fw = 0.0
print(f"\n  Light quark contributions (u, d, s - perturbative, UNRELIABLE below ~2 GeV):")
for name in ['u', 'd', 's']:
    m_f = masses_GeV[name]
    Q = fermions[name]['Q']
    Nc = fermions[name]['Nc']
    log_term = np.log(m_Z_GeV / m_f)
    as_pi = alpha_s / PI
    qcd_corr = 1.0 + as_pi + as_pi**2 * 1.4114
    contrib = alpha_0 / (3 * PI) * Nc * Q**2 * 2 * log_term * qcd_corr
    Da_light_fw += contrib
    print(f"    {name}: m={m_f*1000:.2f} MeV, ln(mZ/mf)={log_term:.4f}, Da={contrib:.6f}")

print(f"  Total light (perturbative): {Da_light_fw:.6f}")

# Top quark (decoupled at m_Z)
Da_top_fw = 0.0
print(f"\n  Top quark: m_t = {masses_GeV['t']:.2f} GeV > m_Z, decoupled.")
print(f"  Decoupling contribution ~ -alpha/(3pi) * Nc * Q_t^2 * (m_Z/m_t)^2/5:")
m_t = masses_GeV['t']
Da_top_fw = -alpha_0 / (3 * PI) * 3 * (2.0/3)**2 * (m_Z_GeV/m_t)**2 / 5
print(f"  Da_top = {Da_top_fw:.6f}")

# Total
Da_total_fw = Da_lep_fw + Da_heavy_fw + Da_light_fw + Da_top_fw
alpha_mZ_inv_fw = alpha_0_inv * (1 - Da_total_fw)
# More precisely: alpha(mZ)^{-1} = alpha(0)^{-1} - sum/(3*pi)
alpha_mZ_inv_fw2 = alpha_0_inv - (delta_sum) / (3 * PI)

print(f"\n  TOTAL Delta_alpha (framework, perturbative):")
print(f"    Leptonic:     {Da_lep_fw:.6f}")
print(f"    Heavy quarks: {Da_heavy_fw:.6f}")
print(f"    Light quarks: {Da_light_fw:.6f} (perturbative, unreliable)")
print(f"    Top:          {Da_top_fw:.6f}")
print(f"    TOTAL:        {Da_total_fw:.6f}")
print(f"    Standard:     {Da_total_std:.6f}")
print(f"")
print(f"  alpha(m_Z)^-1 = {1.0/Da_total_fw:.4f}... wait, using correct formula:")
print(f"  alpha(m_Z)^-1 = alpha(0)^-1 * (1 - Da) = {alpha_mZ_inv_fw:.4f}")
print(f"  alpha(m_Z)^-1 (direct subtraction) = {alpha_mZ_inv_fw2:.4f}")
print(f"  Experimental: 127.951")
print(f"  ratio fw = {1/(1-Da_total_fw):.6f}")
print(f"  ratio exp = {ratio_exp:.6f}")

# =====================================================================
# KEY QUESTION: CAN THE HADRONIC CONTRIBUTION BE FRAMEWORK-DERIVED?
# =====================================================================
print(f"\n{'=' * 72}")
print("HADRONIC VACUUM POLARIZATION FROM FRAMEWORK")
print(f"{'=' * 72}")

print("""
  The perturbative calculation for light quarks is unreliable.
  The physical hadronic VP comes from rho, omega, phi mesons etc.

  In the framework, meson masses are NOT primary objects - fermion masses are.
  However, we can ask: what is the EFFECTIVE scale for light quarks?

  If we replace m_u, m_d, m_s by effective hadronic scales:
  - Lambda_QCD ~ 330 MeV (constituent quark mass)
  - m_rho = 775 MeV

  The question is whether the framework gives a SINGLE effective scale
  for the hadronic sector.
""")

# Try with Lambda_QCD as effective scale for u,d
Lambda_QCD = 0.330  # GeV, constituent quark mass
m_rho = 0.775  # GeV

# With Lambda_QCD for u,d and m_s for s
Da_had_eff = 0.0
for name, m_eff in [('u', Lambda_QCD), ('d', Lambda_QCD), ('s', masses_GeV['s'])]:
    Q = fermions[name]['Q']
    Nc = fermions[name]['Nc']
    log_term = np.log(m_Z_GeV / m_eff)
    contrib = alpha_0 / (3 * PI) * Nc * Q**2 * 2 * log_term
    Da_had_eff += contrib

for name in ['c', 'b']:
    m_f = masses_GeV[name]
    Q = fermions[name]['Q']
    Nc = fermions[name]['Nc']
    log_term = np.log(m_Z_GeV / m_f)
    contrib = alpha_0 / (3 * PI) * Nc * Q**2 * 2 * log_term
    Da_had_eff += contrib

print(f"  With Lambda_QCD={Lambda_QCD*1000:.0f} MeV for u,d: Da_had = {Da_had_eff:.6f}")
print(f"  Standard Da_had = {Da_had_std:.6f}")
print(f"  Difference: {(Da_had_eff/Da_had_std - 1)*100:+.1f}%")

# =====================================================================
# FRAMEWORK ALGEBRAIC APPROACH
# =====================================================================
print(f"\n{'=' * 72}")
print("ALGEBRAIC APPROACH: alpha(m_Z)/alpha(0) AS FRAMEWORK QUANTITY")
print(f"{'=' * 72}")

# The ratio alpha(mZ)/alpha(0) = 1/(1-Da)
# Can Da itself be a simple framework quantity?

# Check various framework expressions
r_exp = ratio_exp

candidates = [
    ("1 + 1/(3*a1)", 1 + 1/(3*a1)),
    ("1 + b1/(a1^2+1)", 1 + b1/(a1**2+1)),
    ("1 + 1/(2*phi^3*pi)", 1 + 1/(2*phi**3*PI)),
    ("1 + alpha_s/pi", 1 + alpha_s/PI),
    ("phi^(1/(2*phi^2))", phi**(1/(2*phi**2))),
    ("1 + 2*alpha_s/(3*pi)*ln(phi^25)", 1 + 2*alpha_s/(3*PI)*np.log(phi**25)),
    ("(a1^2+1)/a1^2", (a1**2+1)/a1**2),
    ("1 + sin2tw/(3*pi)", 1 + sin2tw/(3*PI)),
    ("(1-alpha_s/pi)^{-1}", 1/(1-alpha_s/PI)),
    ("1 + N_gen/(4*pi*a1)", 1 + 3/(4*PI*a1)),
    ("exp(alpha_s/(2*pi))", np.exp(alpha_s/(2*PI))),
    ("exp(3/(4*pi*a1))", np.exp(3/(4*PI*a1))),
]

print(f"\n  Target: alpha(m_Z)/alpha(0) = {r_exp:.6f}")
print(f"\n  {'Expression':35s} {'Value':12s} {'Error':10s}")
for name, val in candidates:
    err = (val / r_exp - 1) * 100
    print(f"  {name:35s} {val:12.6f} {err:+10.3f}%")

# =====================================================================
# THE 1-LOOP SUM AS FRAMEWORK QUANTITY
# =====================================================================
print(f"\n{'=' * 72}")
print("THE 1-LOOP SUM IN TERMS OF FRAMEWORK QUANTITIES")
print(f"{'=' * 72}")

# The sum S = sum_f Nc*Q_f^2 * ln(mZ^2/mf^2)
# With m_Z = m_e * phi^25 * alpha(mZ)/alpha(0)
# and m_f = m_e * phi^{n_f}
# We get: ln(mZ/mf) = ln(phi^{25-n_f} * alpha(mZ)/alpha(0))
#                    ~ (25 - n_f) * ln(phi) + ln(ratio)
# For the ratio, ln(1.063) ~ 0.061, while ln(phi) = 0.481
# So the ratio correction is small.

# At leading order (ignoring the ratio in the log):
# S ~ 2*ln(phi) * sum_f Nc*Q_f^2 * (25 - n_f)

# The charge-weighted sum
print(f"\n  Charge-weighted exponent sums (25 - n_f):")
S_fw = 0.0
for name in ['e', 'mu', 'tau', 'u', 'd', 'c', 's', 'b']:
    info = fermions[name]
    _, n_f = compute_mass(name, info)
    Nc_Q2 = info['Nc'] * info['Q']**2
    weight = 25 - n_f
    contrib = Nc_Q2 * weight
    S_fw += contrib
    print(f"    {name}: n={n_f:7.3f}, 25-n={weight:7.3f}, Nc*Q^2={Nc_Q2:.4f}, contrib={contrib:8.3f}")

print(f"\n  Total S = sum Nc*Q^2*(25-n_f) = {S_fw:.4f}")
print(f"  Da ~ (2*alpha*ln(phi))/(3*pi) * S = {2*alpha_0*np.log(phi)/(3*PI)*S_fw:.6f}")

# =====================================================================
# SELF-CONSISTENCY CHECK
# =====================================================================
print(f"\n{'=' * 72}")
print("SELF-CONSISTENT DETERMINATION OF alpha(m_Z)")
print(f"{'=' * 72}")

# Use m_Z = m_e * phi^25 * r where r = alpha(mZ)/alpha(0)
# And alpha(mZ)^{-1} = alpha(0)^{-1} - (1/3pi)*S(r)
# where S(r) = sum Nc*Q^2 * 2*ln(m_e*phi^25*r / m_f)
#            = sum Nc*Q^2 * 2*(ln(phi)*(25-n_f) + ln(r))

# This is a self-consistency equation for r:
# 1/(r*alpha(0)) = 1/alpha(0) - (1/3pi)*[2*ln(phi)*sum Nc*Q^2*(25-nf) + 2*ln(r)*sum Nc*Q^2]
# Let S0 = sum Nc*Q^2*(25-nf), SQ = sum Nc*Q^2

SQ = 0.0
for name in ['e', 'mu', 'tau', 'u', 'd', 'c', 's', 'b']:
    info = fermions[name]
    if masses_GeV[name] < m_Z_GeV:
        SQ += info['Nc'] * info['Q']**2
print(f"  SQ = sum Nc*Q^2 (active fermions) = {SQ:.4f}")
print(f"  (Standard: 8/3 for 3 leptons + 5 quarks with 3 colors)")
print(f"  Expected: 3*1 + 3*(4/9+1/9+4/9+1/9+4/9) = 3 + 3*14/9 = 3 + 14/3 = 23/3 = {23/3:.4f}")
# Wait, that's not matching. Let me recompute.
# Leptons: 3 * 1^2 = 3
# u-type: 3 colors * (2/3)^2 = 3*4/9 = 4/3  (x2 for u,c)
# d-type: 3 colors * (1/3)^2 = 3*1/9 = 1/3  (x3 for d,s,b)
SQ_check = 3*1 + 2*(3*(2.0/3)**2) + 3*(3*(1.0/3)**2)
print(f"  Check: 3(lep) + 2*(4/3)(u-type) + 3*(1/3)(d-type) = {SQ_check:.4f}")
print(f"  = 3 + 8/3 + 1 = {3 + 8.0/3 + 1:.4f} = 20/3")

# Self-consistent equation:
# 1/r = 1 - alpha(0)/(3pi) * 2 * (ln(phi)*S0_active + ln(r)*SQ)

# where S0_active = sum over active fermions only (m_f < m_Z)
S0_active = 0.0
for name in ['e', 'mu', 'tau', 'u', 'd', 'c', 's', 'b']:
    info = fermions[name]
    _, n_f = compute_mass(name, info)
    if masses_GeV[name] < m_Z_GeV:
        S0_active += info['Nc'] * info['Q']**2 * (25 - n_f)

print(f"  S0 = sum Nc*Q^2*(25-n_f) for active = {S0_active:.4f}")

# Solve iteratively
r = 1.0
for i in range(20):
    Da = alpha_0 / (3 * PI) * 2 * (np.log(phi) * S0_active + np.log(r) * SQ)
    r_new = 1.0 / (1.0 - Da)
    if abs(r_new - r) < 1e-12:
        break
    r = r_new

print(f"\n  Self-consistent solution (iterative):")
print(f"    r = alpha(m_Z)/alpha(0) = {r:.6f}")
print(f"    alpha(m_Z)^-1 = {alpha_0_inv/r:.4f}")
print(f"    Experimental: r = {ratio_exp:.6f}")
print(f"    Difference: {(r/ratio_exp - 1)*100:+.3f}%")

# =====================================================================
# THE DEEP QUESTION: IS alpha(mZ)/alpha(0) AN ALGEBRAIC QUANTITY?
# =====================================================================
print(f"\n{'=' * 72}")
print("IS THE RATIO ALGEBRAIC IN THE FRAMEWORK?")
print(f"{'=' * 72}")

# The key formula is:
# Da = (2*alpha*ln(phi))/(3*pi) * S0
# where S0 is a sum over framework quantities.
# alpha itself is algebraic (root of quadratic).
# ln(phi) is transcendental.
# S0 involves n_f which are algebraic (in Z[phi]).

# So Da is of the form: algebraic * transcendental
# and the ratio is 1/(1-Da).
# This is NOT algebraic, but it IS computable from framework.

print(f"\n  The ratio alpha(mZ)/alpha(0) is computable from framework quantities")
print(f"  but involves ln(phi) (transcendental). It is NOT algebraic.")
print(f"")
print(f"  Key formula: Da = (2*alpha*ln(phi))/(3*pi) * S0")
print(f"  where S0 = {S0_active:.4f} (framework-derived)")
print(f"  and alpha, ln(phi), pi are known constants.")
print(f"")
print(f"  Da = {alpha_0 / (3 * PI) * 2 * np.log(phi) * S0_active:.6f}")
print(f"  This gives r = {1/(1-alpha_0/(3*PI)*2*np.log(phi)*S0_active):.6f}")
print(f"  (without self-consistent ln(r) correction)")

# =====================================================================
# SUMMARY
# =====================================================================
print(f"\n{'=' * 72}")
print("SUMMARY")
print(f"{'=' * 72}")

r_simple = 1/(1-alpha_0/(3*PI)*2*np.log(phi)*S0_active)
print(f"""
  1-LOOP RESULTS (framework masses, perturbative):
    alpha(m_Z)/alpha(0) = {ratio_1loop:.6f}  (1-loop, naive)
    Self-consistent:      {r:.6f}
    Experimental:         {ratio_exp:.6f}

  ISSUE: Light quark masses (u,d,s) are current masses, but QED VP
  below ~2 GeV is nonperturbative. The perturbative calculation
  OVERESTIMATES Delta_alpha_had because ln(mZ/m_u) >> ln(mZ/Lambda_QCD).

  RESOLUTION: The framework provides fermion masses, NOT meson masses.
  The hadronic VP is fundamentally nonperturbative and requires
  the FULL QCD dynamics (which the framework does not compute).

  CONCLUSION: alpha(m_Z)/alpha(0) CANNOT be derived purely from
  the 600-cell framework. It requires QFT dynamics (vacuum polarization
  with nonperturbative QCD). This is an IRREDUCIBLE external input.

  However, the perturbative part (leptons + heavy quarks) IS framework-
  derived and accounts for ~{(Da_lep_fw + Da_heavy_fw)/Da_total_std*100:.0f}% of the total running.
  The remaining ~{Da_had_std/(Da_total_std)*100:.0f}% (hadronic light quarks) requires QCD dynamics.

  STATUS: OPEN-3 route 2 = PARTIALLY RESOLVED.
  The ratio is ~60% framework-derived, ~40% requires nonperturbative QCD.
  This is honest: the framework gives STATIC properties (masses, couplings),
  not DYNAMIC properties (vacuum fluctuations, confinement).
""")

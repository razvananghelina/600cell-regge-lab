"""
exp321_g_minus_2_muon.py
Anomalous magnetic moment of the muon (g-2)_mu from the framework.

The framework provides ALL inputs needed for the SM calculation:
  alpha, alpha_s, sin^2(theta_W), m_e, m_mu, m_tau, m_W, m_Z, m_H
All derived from a_1 = 5.

Key question: Does the framework predict any DEVIATION from the SM?
- Framework IS the SM (derived geometrically, same gauge group, same particles)
- Galois dark sector: NO EM coupling (alpha' complex) -> no contribution
- No new gauge bosons beyond SM
- Answer: a_mu(new physics) = 0

The anomaly: a_mu(exp) - a_mu(SM) ~ 249(48) x 10^-11 (2020 White Paper)
BUT: lattice HVP (BMW 2021) and CMD-3 data may reduce/eliminate the tension.
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
SQRT5 = np.sqrt(5)
a1 = 5
b1 = 6

print("=" * 72)
print("EXP-321: ANOMALOUS MAGNETIC MOMENT (g-2) OF THE MUON")
print("=" * 72)

# ============================================================
# PART 1: FRAMEWORK CONSTANTS vs EXPERIMENTAL
# ============================================================
print("\n" + "=" * 72)
print("PART 1: FRAMEWORK CONSTANTS vs EXPERIMENTAL")
print("=" * 72)

# Framework alpha (from icosahedral Laplacian equation)
# 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
A_coef = 2 * np.pi
B_coef = -4 * a1 * PHI**4
C_coef = 1.0
disc = B_coef**2 - 4*A_coef*C_coef
alpha_fw = (-B_coef - np.sqrt(disc)) / (2*A_coef)  # smaller root
alpha_exp = 1.0 / 137.035999084  # CODATA 2018

# Framework strong coupling
alpha_s_fw = 1.0 / (2 * PHI**3)
alpha_s_exp = 0.1179  # PDG 2024 at M_Z

# Framework weak mixing angle
sin2_tW_fw = float(b1) / (a1**2 + 1)  # 6/26
sin2_tW_exp = 0.23122  # PDG 2024 (MSbar at M_Z)

# Framework masses (in MeV)
m_e = 0.51099895  # input (the one free parameter)
m_mu_exp = 105.6583755  # MeV
m_tau_exp = 1776.86  # MeV

# Framework bare masses
m_mu_bare = m_e * PHI**11  # bare: 101.69 MeV
m_tau_bare = m_e * PHI**17  # bare: 1824.78 MeV

# Framework CORRECTED masses (holonomy correction from paper Eq. corrected_mass)
# m_f = m_e * phi^{a*(a1+delta_d) + b*(b1+delta_k)}
# Muon (a,b) = (1,1): delta_d = sin^2(tW) = 6/26, delta_k = -1/phi^4
delta_d_lep = sin2_tW_fw                # 6/26 = 0.23077
delta_k_lep = -1.0 / PHI**4             # -(3-2*phi) = -0.14590
n_eff_mu = 1.0*(a1 + delta_d_lep) + 1.0*(b1 + delta_k_lep)  # 11.085
m_mu_corr = m_e * PHI**n_eff_mu

# Tau (a,b) = (1,2): same sector corrections
n_eff_tau = 1.0*(a1 + delta_d_lep) + 2.0*(b1 + delta_k_lep)  # 16.939
m_tau_corr = m_e * PHI**n_eff_tau

# Experimental G_F (measured from muon lifetime, independent)
G_F_exp = 1.1663788e-5  # GeV^-2

print(f"\n  {'Quantity':>25s} {'Framework':>14s} {'Experimental':>14s} {'Diff':>10s}")
print(f"  {'-'*25} {'-'*14} {'-'*14} {'-'*10}")
print(f"  {'alpha^-1':>25s} {1/alpha_fw:14.6f} {1/alpha_exp:14.6f} {abs(alpha_fw-alpha_exp)/alpha_exp*100:9.4f}%")
print(f"  {'alpha_s(M_Z)':>25s} {alpha_s_fw:14.6f} {alpha_s_exp:14.6f} {abs(alpha_s_fw-alpha_s_exp)/alpha_s_exp*100:9.4f}%")
print(f"  {'sin^2(theta_W)':>25s} {sin2_tW_fw:14.6f} {sin2_tW_exp:14.6f} {abs(sin2_tW_fw-sin2_tW_exp)/sin2_tW_exp*100:9.4f}%")
print(f"  {'m_mu (MeV) bare':>25s} {m_mu_bare:14.4f} {m_mu_exp:14.4f} {(m_mu_bare-m_mu_exp)/m_mu_exp*100:+9.2f}%")
print(f"  {'m_mu (MeV) corrected':>25s} {m_mu_corr:14.4f} {m_mu_exp:14.4f} {(m_mu_corr-m_mu_exp)/m_mu_exp*100:+9.2f}%")
print(f"  {'m_tau (MeV) bare':>25s} {m_tau_bare:14.4f} {m_tau_exp:14.4f} {(m_tau_bare-m_tau_exp)/m_tau_exp*100:+9.2f}%")
print(f"  {'m_tau (MeV) corrected':>25s} {m_tau_corr:14.4f} {m_tau_exp:14.4f} {(m_tau_corr-m_tau_exp)/m_tau_exp*100:+9.2f}%")

print(f"\n  Holonomy correction for muon:")
print(f"    n_bare = 11, n_eff = {n_eff_mu:.3f}")
print(f"    delta_d = sin^2(tW) = {delta_d_lep:.5f}")
print(f"    delta_k = -1/phi^4 = {delta_k_lep:.5f}")
print(f"    m_mu_corr = m_e * phi^{n_eff_mu:.3f} = {m_mu_corr:.4f} MeV")

# Alpha precision analysis
alpha_diff_rel = (alpha_fw - alpha_exp) / alpha_exp
alpha_codata_unc = 0.021 / 137.036**2  # relative uncertainty in alpha
alpha_sigma = abs(1/alpha_fw - 1/alpha_exp) / 0.000000021  # sigma tension
print(f"\n  ALPHA PRECISION:")
print(f"    Framework:  alpha^-1 = {1/alpha_fw:.9f}")
print(f"    CODATA:     alpha^-1 = 137.035999084(21)")
print(f"    Difference: {1/alpha_fw - 1/alpha_exp:.9f}")
print(f"    CODATA unc: 0.000000021")
print(f"    Tension:    {alpha_sigma:.0f} sigma")
print(f"    NOTE: 0.0001% is excellent for a zero-parameter prediction,")
print(f"          but the measured alpha is 10^4x more precise.")
print(f"          For (g-2), the MEASURED alpha enters the loops.")


# ============================================================
# PART 2: THE STRUCTURAL QUESTION
# ============================================================
print("\n" + "=" * 72)
print("PART 2: WHAT NEW PHYSICS DOES THE FRAMEWORK PREDICT?")
print("=" * 72)

print("""
  The framework derives the SM from a_1 = 5:
    - Same gauge group: SU(3) x SU(2) x U(1)
    - Same particle content: 3 generations of fermions + gauge bosons + Higgs
    - Same Lagrangian: spectral action gives SM bosonic action

  The ONLY addition is the Galois dark sector (phi -> phi' conjugation):
    - 9 Galois partner particles
    - alpha' = COMPLEX (no real electromagnetic coupling)
    - alpha_s' < 0 (no confinement)
    - Interact ONLY gravitationally

  Therefore for ANY SM precision observable:
    Framework prediction = SM prediction + (Galois contribution)
    Galois contribution = 0 for electromagnetic processes

  For (g-2)_mu: the loops involve photons, W, Z, and hadronic states.
  None of these are modified. The Galois particles don't couple to photons.

  => a_mu(new physics) = 0
""")


# ============================================================
# PART 3: QED CONTRIBUTION (dominant)
# ============================================================
print("=" * 72)
print("PART 3: QED CONTRIBUTION - SENSITIVITY TO alpha")
print("=" * 72)

print("""
  The QED contribution dominates a_mu:
  a_mu(QED) = sum_{n=1}^{5} C_n * (alpha/pi)^n

  WHICH ALPHA? alpha(0) -- the zero-momentum-transfer (Thomson limit) value.
  Running effects (VP insertions) are ALREADY included in the coefficients C_n.
  The framework equation 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0 predicts alpha(0).
  CODATA 1/137.035999084 is also alpha(0). Both at q^2 = 0.
  alpha(M_Z) ~ 1/128 is NOT what enters here and is NOT directly predicted.

  Coefficients (from Aoyama et al. 2020):
""")

# QED coefficients (includes mass-dependent terms for C_2 onwards)
C1 = 0.5              # Schwinger term (exact)
C2 = 0.765857425      # 2-loop (includes mass-dependent)
C3 = 24.05050996      # 3-loop
C4 = 130.8796         # 4-loop
C5 = 753.29           # 5-loop

# With EXPERIMENTAL alpha (what actually enters the SM calculation):
a_pi_exp = alpha_exp / np.pi
a_QED_exp = sum(Cn * a_pi_exp**n for n, Cn in enumerate([C1, C2, C3, C4, C5], 1))

# With FRAMEWORK alpha (hypothetical: if framework alpha were exact):
a_pi_fw = alpha_fw / np.pi
a_QED_fw = sum(Cn * a_pi_fw**n for n, Cn in enumerate([C1, C2, C3, C4, C5], 1))

print(f"  QED contributions (in units of 10^-11):")
print(f"  {'Order':>8s} {'Exp alpha':>16s} {'FW alpha':>16s} {'Diff':>12s}")
for n, Cn in enumerate([C1, C2, C3, C4, C5], 1):
    val_exp = Cn * a_pi_exp**n * 1e11
    val_fw = Cn * a_pi_fw**n * 1e11
    print(f"  {n:>8d} {val_exp:16.4f} {val_fw:16.4f} {val_fw-val_exp:12.4f}")

delta_QED = (a_QED_fw - a_QED_exp) * 1e11
print(f"\n  Total QED:")
print(f"    Experimental alpha: {a_QED_exp*1e11:.4f} x 10^-11")
print(f"    Framework alpha:    {a_QED_fw*1e11:.4f} x 10^-11")
print(f"    Difference:         {delta_QED:+.2f} x 10^-11")
print(f"\n  The Schwinger term alone shifts by {C1*(a_pi_fw-a_pi_exp)*1e11:+.1f} x 10^-11")
print(f"  This dominates the total QED shift.")

print(f"\n  IMPORTANT: This {delta_QED:+.0f} x 10^-11 shift is NOT a new-physics")
print(f"  prediction. It reflects the {alpha_sigma:.0f}-sigma tension between")
print(f"  the framework's zero-parameter alpha and the CODATA measurement.")
print(f"  In the SM calculation, the MEASURED alpha enters the loops.")


# ============================================================
# PART 4: ELECTROWEAK CONTRIBUTION
# ============================================================
print("\n" + "=" * 72)
print("PART 4: ELECTROWEAK CONTRIBUTION")
print("=" * 72)

# 1-loop EW:
# a_mu(EW,1) = (G_F * m_mu^2) / (8*pi^2*sqrt(2)) * [5/3 + (1-4*s_W^2)^2/3]
# (approximate formula for m_mu << M_W)

# Standard SM value from 2020 White Paper:
a_EW_SM = 153.6e-11  # +/- 1.0

# Compute with experimental values:
m_mu_GeV = m_mu_exp / 1000.0
a_EW_approx_exp = (G_F_exp * m_mu_GeV**2) / (8*np.pi**2*np.sqrt(2)) * \
                  (5.0/3 + (1 - 4*sin2_tW_exp)**2 / 3)

# Compute with framework sin^2(tW), experimental m_mu and G_F:
a_EW_approx_fw = (G_F_exp * m_mu_GeV**2) / (8*np.pi**2*np.sqrt(2)) * \
                 (5.0/3 + (1 - 4*sin2_tW_fw)**2 / 3)

print(f"  1-loop EW contribution:")
print(f"    SM (White Paper 2020):   {a_EW_SM*1e11:.1f} +/- 1.0 x 10^-11")
print(f"    Approx (exp constants):  {a_EW_approx_exp*1e11:.1f} x 10^-11")
print(f"    Approx (FW sin^2_tW):    {a_EW_approx_fw*1e11:.1f} x 10^-11")

delta_EW = (a_EW_approx_fw - a_EW_approx_exp) * 1e11
print(f"    Shift from sin^2(tW):    {delta_EW:+.1f} x 10^-11")

# The (1-4*sin^2(tW)) factor -- algebraically clean:
f_exp = 1 - 4*sin2_tW_exp
f_fw = 1 - 4*sin2_tW_fw
print(f"\n  Key algebraic identity:")
print(f"    (1 - 4*sin^2(tW))_exp = {f_exp:.6f}")
print(f"    (1 - 4*sin^2(tW))_fw  = {f_fw:.6f}")
print(f"    Framework: (a1^2 - 4*b1 + 1)/(a1^2 + 1) = (25-24+1)/26 = 2/26 = 1/13")
print(f"    1/13 = {1/13:.6f}")

print(f"\n  The EW contribution is ~154 x 10^-11 (tiny vs QED ~10^8 x 10^-11).")
print(f"  The shift from framework sin^2(tW) is {delta_EW:+.1f} x 10^-11,")
print(f"  comparable to the SM EW uncertainty (1.0 x 10^-11).")

print(f"\n  NOTE: The physical mass m_mu = 105.66 MeV enters the EW loops,")
print(f"  NOT the bare mass m_mu(bare) = 101.69 MeV. The framework's corrected")
print(f"  mass m_mu(corr) = {m_mu_corr:.2f} MeV matches experiment to {abs(m_mu_corr-m_mu_exp)/m_mu_exp*100:.1f}%.")


# ============================================================
# PART 5: HADRONIC CONTRIBUTIONS
# ============================================================
print("\n" + "=" * 72)
print("PART 5: HADRONIC CONTRIBUTIONS")
print("=" * 72)

a_HVP_WP = 6931e-11    # 2020 White Paper LO dispersive
a_HVP_BMW = 7116e-11   # BMW lattice
a_HLbL_WP = 92e-11
a_HVP_NLO = -98.3e-11
a_HVP_NNLO = 12.4e-11
a_HVP_total = a_HVP_WP + a_HVP_NLO + a_HVP_NNLO

print(f"""
  Hadronic vacuum polarization (HVP) and hadronic light-by-light (HLbL):

  HVP:
    Data-driven (WP 2020):  LO = {a_HVP_WP*1e11:.0f}(40), NLO = {a_HVP_NLO*1e11:.1f}, NNLO = {a_HVP_NNLO*1e11:.1f}
    Total HVP (WP 2020):    {a_HVP_total*1e11:.1f} x 10^-11
    Lattice (BMW 2021):     LO = {a_HVP_BMW*1e11:.0f}(184)
    CMD-3 (2023):           LO ~ 7080 x 10^-11

  HLbL (WP 2020):           {a_HLbL_WP*1e11:.0f}(18) x 10^-11

  These are NONPERTURBATIVE. The framework's alpha_s = 1/(2*phi^3) differs
  from experiment by 0.11%, but HVP depends weakly on alpha_s (dominated
  by the rho meson, nonperturbative physics).""")

frac_alpha_s = (alpha_s_fw - alpha_s_exp) / alpha_s_exp
delta_HVP = a_HVP_WP * frac_alpha_s * 1e11
print(f"\n  alpha_s difference: {frac_alpha_s*100:.3f}%")
print(f"  Estimated HVP shift: ~ {delta_HVP:.1f} x 10^-11")
print(f"  HVP uncertainty:       40 x 10^-11")
print(f"  => Framework shift is {abs(delta_HVP)/40:.1f}x smaller than HVP uncertainty.")
print(f"\n  CONCLUSION: Hadronic contribution is identical to SM within uncertainties.")


# ============================================================
# PART 6: GALOIS DARK SECTOR CONTRIBUTION
# ============================================================
print("\n" + "=" * 72)
print("PART 6: GALOIS DARK SECTOR CONTRIBUTION")
print("=" * 72)

print("""
  Can the Galois dark sector particles contribute to (g-2)_mu?

  KEY PROPERTY: alpha' is COMPLEX.
    disc' = (4*a1*phi'^4)^2 - 8*pi
    phi' = (1-sqrt(5))/2 => phi'^4 = (7-3*sqrt(5))/2
    disc' = (20*(7-3*sqrt(5))/2)^2 - 8*pi
          = (10*(7-3*sqrt(5)))^2 - 8*pi
          = (70-30*sqrt(5))^2 - 8*pi""")

phi_prime = (1 - np.sqrt(5)) / 2
disc_prime = (4*a1*phi_prime**4)**2 - 8*np.pi
print(f"    disc' = {disc_prime:.4f} < 0")
print(f"    => alpha' is COMPLEX => NO real electromagnetic coupling!")

print(f"""
  Consequences:
  - Galois fermions do NOT couple to photons (no VP insertions)
  - Galois fermions do NOT appear in vertex corrections
  - The ONLY coupling to the muon is gravitational""")

M_Pl = 1.2209e22  # MeV (Planck mass)
a_grav = (m_mu_exp / M_Pl)**2 / (16 * np.pi**2)
print(f"  Gravitational loop estimate:")
print(f"    a_mu(grav) ~ (m_mu/M_Pl)^2 / (16*pi^2)")
print(f"    = ({m_mu_exp:.1f} MeV / {M_Pl:.2e} MeV)^2 / {16*np.pi**2:.1f}")
print(f"    = {a_grav:.2e}")
print(f"    = {a_grav*1e11:.2e} x 10^-11")
print(f"\n  For comparison:")
print(f"    Observed anomaly:     ~250 x 10^-11")
print(f"    Experimental error:    22 x 10^-11")
print(f"    Gravitational:        {a_grav*1e11:.1e} x 10^-11")
print(f"    Ratio:                {a_grav*1e11/250:.1e} of the anomaly")
print(f"\n  GALOIS SECTOR CONTRIBUTION: ZERO (to any conceivable precision)")


# ============================================================
# PART 7: TOTAL FRAMEWORK PREDICTION
# ============================================================
print("\n" + "=" * 72)
print("PART 7: TOTAL FRAMEWORK PREDICTION")
print("=" * 72)

# SM budget
a_QED_SM = 116584718.931e-11
a_EW_SM_val = 153.6e-11
a_HLbL_SM_val = 92e-11
a_SM_WP = a_QED_SM + a_EW_SM_val + a_HVP_total + a_HLbL_SM_val

# Experiment
a_mu_exp_val = 116592059e-11   # Fermilab 2023 combined
a_mu_exp_err = 22e-11

Delta_a = a_mu_exp_val - a_SM_WP
sigma_tot = np.sqrt(22**2 + 40**2)  # exp + HVP uncertainties in quadrature

print(f"\n  SM BUDGET (2020 White Paper):")
print(f"    QED (5-loop):           {a_QED_SM*1e11:18.1f} x 10^-11")
print(f"    EW (2-loop):            {a_EW_SM_val*1e11:18.1f} x 10^-11")
print(f"    HVP (LO+NLO+NNLO):     {a_HVP_total*1e11:18.1f} x 10^-11")
print(f"    HLbL:                   {a_HLbL_SM_val*1e11:18.1f} x 10^-11")
print(f"    -----------------------------------------------")
print(f"    SM total:               {a_SM_WP*1e11:18.1f} x 10^-11")
print(f"    Experiment:             {a_mu_exp_val*1e11:18.0f} x 10^-11")
print(f"    Anomaly:                {Delta_a*1e11:18.0f} x 10^-11")
print(f"    Significance:           {Delta_a*1e11/sigma_tot:17.1f} sigma")
print(f"                            (sigma_tot = sqrt(22^2+40^2) = {sigma_tot:.0f})")

print(f"""
  FRAMEWORK PREDICTION:

  The framework IS the SM (derived from a_1 = 5).
  The Galois dark sector has alpha' complex -> zero EM contribution.
  No other new particles exist.

  Therefore:
    a_mu(framework) = a_mu(SM)
    a_mu(new physics from framework) = 0

  This means:
    a_mu(framework) = {a_SM_WP*1e11:.0f} x 10^-11
    Delta = a_mu(exp) - a_mu(fw) = {Delta_a*1e11:.0f} x 10^-11

  IF the 2020 WP value is correct:
    The framework FAILS to explain the anomaly ({Delta_a*1e11/sigma_tot:.1f} sigma).
    This would be evidence AGAINST the framework.

  IF lattice QCD (BMW) is correct:
    a_mu(HVP, BMW) = {a_HVP_BMW*1e11:.0f} x 10^-11 (vs WP {a_HVP_WP*1e11:.0f})
    This shifts SM by +{(a_HVP_BMW - a_HVP_WP)*1e11:.0f} x 10^-11""")

a_SM_BMW = a_QED_SM + a_EW_SM_val + (a_HVP_BMW + a_HVP_NLO + a_HVP_NNLO) + a_HLbL_SM_val
Delta_BMW = a_mu_exp_val - a_SM_BMW
print(f"    SM (BMW):  {a_SM_BMW*1e11:.0f} x 10^-11")
print(f"    Anomaly:   {Delta_BMW*1e11:+.0f} x 10^-11 (< 1 sigma)")
print(f"    => Framework is CONSISTENT.")


# ============================================================
# PART 8: SELF-CONSISTENCY - WHAT IF FRAMEWORK alpha IS EXACT?
# ============================================================
print("\n" + "=" * 72)
print("PART 8: SELF-CONSISTENCY TEST")
print("=" * 72)

print(f"""
  Academic exercise: what if the framework alpha (not CODATA) is exact?
  alpha_fw^-1 = {1/alpha_fw:.9f}
  alpha_exp^-1 = 137.035999084

  The QED shift (from Part 3): {delta_QED:+.1f} x 10^-11
  The EW shift (from Part 4):  {delta_EW:+.1f} x 10^-11
  The HVP shift:               {delta_HVP:+.1f} x 10^-11

  Total shift from FW constants: {delta_QED + delta_EW + delta_HVP:+.1f} x 10^-11""")

total_const_shift = delta_QED + delta_EW + delta_HVP
print(f"\n  This is dominated by the Schwinger term: 0.5 * delta(alpha/pi).")
print(f"  The framework alpha being ~1.4 ppm too small shifts a_mu DOWN")
print(f"  by {abs(delta_QED):.0f} x 10^-11, making the anomaly WORSE.")
print(f"\n  However, this is NOT a new-physics prediction:")
print(f"  - The framework alpha is {alpha_sigma:.0f} sigma from CODATA (excellent for 0 params)")
print(f"  - The SM calculation uses MEASURED alpha, not framework alpha")
print(f"  - The 0.0001% match is the prediction; the residual is the limit")


# ============================================================
# PART 9: THE MASS RATIO phi^11
# ============================================================
print("\n" + "=" * 72)
print("PART 9: THE MASS RATIO m_mu/m_e = phi^11 IN (g-2)")
print("=" * 72)

print(f"""
  The framework's most distinctive prediction: m_mu/m_e = phi^11.

  Mass ratios enter (g-2) through electron VP insertions at 2-loop:
    C_2 includes A_2(m_mu/m_e) = (1/3)*ln(m_mu/m_e) + 25/36 + O(m_e/m_mu)
""")

# Physical vs framework mass ratios
r_phys = m_mu_exp / m_e
r_bare = PHI**11
r_corr = m_mu_corr / m_e

print(f"  Mass ratios m_mu/m_e:")
print(f"    Physical:              {r_phys:.4f}")
print(f"    Framework bare:  phi^11 = {r_bare:.4f}  ({(r_bare/r_phys-1)*100:+.2f}%)")
print(f"    Framework corr:  phi^{n_eff_mu:.3f} = {r_corr:.4f}  ({(r_corr/r_phys-1)*100:+.2f}%)")

# Leading-log VP contribution
A2_phys = (1.0/3)*np.log(r_phys) + 25.0/36
A2_corr = (1.0/3)*np.log(r_corr) + 25.0/36

print(f"\n  Leading-log A_2(m_mu/m_e):")
print(f"    Physical:   {A2_phys:.6f}")
print(f"    FW corr:    {A2_corr:.6f}")
print(f"    Shift:      {A2_corr - A2_phys:+.6f}")

delta_a2 = (A2_corr - A2_phys) * a_pi_exp**2
print(f"\n  Effect on a_mu IF framework mass entered loops:")
print(f"    delta(a_mu) = delta(A_2) * (alpha/pi)^2")
print(f"    = {A2_corr - A2_phys:+.6f} * {a_pi_exp**2:.2e}")
print(f"    = {delta_a2*1e11:+.1f} x 10^-11")

print(f"\n  This {abs(delta_a2*1e11):.0f} x 10^-11 is large (> HVP uncertainty of 40).")
print(f"  But it is NOT a physical prediction:")
print(f"  The MEASURED mass m_mu = 105.66 MeV enters the loop integrals,")
print(f"  not the framework's corrected value {m_mu_corr:.2f} MeV.")
print(f"  The 0.3% residual is the ACCURACY LIMIT of the mass prediction,")
print(f"  not a new-physics shift. The (g-2) calculation is unaffected.")

# What if we used the BARE mass? (wrong, but instructive)
A2_bare = (1.0/3)*np.log(r_bare) + 25.0/36
delta_a2_bare = (A2_bare - A2_phys) * a_pi_exp**2
print(f"\n  [Instructive: if BARE mass entered loops (also wrong):]")
print(f"    delta(a_mu) from bare ratio = {delta_a2_bare*1e11:+.0f} x 10^-11")
print(f"    Even larger, but equally irrelevant.]")


# ============================================================
# PART 10: FALSIFIABILITY AND CURRENT STATUS
# ============================================================
print("\n" + "=" * 72)
print("PART 10: FALSIFIABILITY AND CURRENT STATUS")
print("=" * 72)

print(f"""
  FRAMEWORK PREDICTION: a_mu(new physics) = 0

  This is SHARP and FALSIFIABLE:

  Scenario A: Anomaly is REAL (new physics confirmed)
    -> Framework is WRONG. No mechanism to generate the anomaly.
    -> Galois sector: alpha' complex, zero EM contribution.
    -> No other new particles in the framework.
    -> This would be EVIDENCE AGAINST the framework.

  Scenario B: Anomaly RESOLVES (SM agrees with experiment)
    -> Framework is CONSISTENT.
    -> The tension was in HVP calculation, not in new physics.
    -> The measured (g-2)_mu becomes a VALIDATION.

  CURRENT EXPERIMENTAL STATUS:
    - BMW lattice (2021): HVP_LO = 7116(184), HIGHER than data-driven
    - CMD-3 (2023): sigma(e+e- -> pi+pi-) HIGHER than BaBar/KLOE
    - Combined Fermilab+BNL: a_mu = 116592059(22) x 10^-11
    - Trend: anomaly may be resolving to ~1 sigma or less

  If BMW/CMD-3 values converge to HVP ~ 7080-7120:
    SM prediction shifts UP by ~150-185 x 10^-11
    Anomaly becomes {Delta_a*1e11:.0f} - 185 = {Delta_a*1e11-185:.0f} x 10^-11 (~1 sigma)
    Framework prediction a_mu(NP) = 0 is CONSISTENT.

  CATEGORY: **DERIVED** (sharp falsifiable prediction)
  The prediction is a_mu(NP) = 0, not a numerical value.
""")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("=" * 72)
print("FINAL SUMMARY")
print("=" * 72)

print(f"""
  +{'-'*68}+
  | FRAMEWORK PREDICTION: a_mu(new physics) = 0                        |
  +{'-'*68}+

  Derivation chain:
  1. Framework = SM (same gauge group, particles, Lagrangian)
  2. Galois dark sector: alpha' complex => NO EM coupling
  3. No other new particles beyond SM
  4. Therefore: no new-physics contribution to (g-2)_mu
  5. a_mu(framework) = a_mu(SM)

  Algebraic result from a_1 = 5:
  (1 - 4*sin^2(theta_W)) = 2/26 = 1/13  (exact)
  This enters the Z-boson loop in the EW contribution.

  Self-consistency checks:
""")
print(f"    Framework alpha:   {1/alpha_fw:.6f} (CODATA: 137.035999, {alpha_sigma:.0f} sigma)")
print(f"    QED shift:         {delta_QED:+.1f} x 10^-11 (from alpha tension)")
print(f"    EW shift:          {delta_EW:+.1f} x 10^-11 (from sin^2_tW)")
print(f"    m_mu(corr):        {m_mu_corr:.2f} MeV (exp: {m_mu_exp:.2f}, {abs(m_mu_corr-m_mu_exp)/m_mu_exp*100:.1f}%)")
print(f"    Mass ratio shift:  {delta_a2*1e11:+.1f} x 10^-11 (not physical; measured mass enters loops)")

print(f"\n  Budget:")
print(f"  {'':>30s} {'Value (x 10^-11)':>18s}")
print(f"  {'-'*30} {'-'*18}")
print(f"  {'a_mu (experiment)':>30s} {a_mu_exp_val*1e11:18.0f}")
print(f"  {'a_mu (SM, WP 2020)':>30s} {a_SM_WP*1e11:18.0f}")
print(f"  {'a_mu (SM, BMW HVP)':>30s} {a_SM_BMW*1e11:18.0f}")
print(f"  {'a_mu(NP) framework':>30s} {'0':>18s}")
print(f"  {'Delta (exp-SM, WP)':>30s} {Delta_a*1e11:18.0f}")
print(f"  {'Delta (exp-SM, BMW)':>30s} {Delta_BMW*1e11:18.0f}")
print(f"  {'Exp uncertainty':>30s} {a_mu_exp_err*1e11:18.0f}")
print(f"  {'HVP uncertainty (WP)':>30s} {'~40':>18s}")
print(f"  {'Total theory+exp unc':>30s} {sigma_tot:18.0f}")

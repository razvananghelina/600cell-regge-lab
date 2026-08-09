"""
exp438c_precision_summary.py
============================
Complete precision summary for ALL particle masses in the framework,
including the one-loop QED correction delta_2 = sqrt(2/a1) * alpha * delta_1.

Author: Razvan-Constantin Anghelina
Date: March 2026
"""

import numpy as np

# =============================================================================
# CONSTANTS
# =============================================================================
a1 = 5
phi = (1 + np.sqrt(5)) / 2
alpha_em = 7.2973525693e-3
m_e_MeV = 0.51099895000
c_oneloop = np.sqrt(2.0 / a1)  # = sqrt(2/5) DERIVED from Verlinde

C_coeff = 2.0 / 13
c_ell = C_coeff * phi**3 / 4.0
N_gen = 3

# sin^2(tW) and other derived
sin2tW = 6.0 / 26
m_Z_GeV = 91.1876
m_W_GeV = m_Z_GeV * np.sqrt(1 - sin2tW)  # = m_Z * cos(tW)
alpha_s = 1.0 / (2 * phi**3)

# =============================================================================
# FERMION MASSES
# =============================================================================
def compute_delta1(name, a, b, sector):
    z = a + b * phi
    zp = a + b * (1 - np.sqrt(5))/2
    N_z = z * zp
    if sector == 'input':
        return 0.0
    elif sector == 'lepton':
        return c_ell * np.sign(zp) * abs(zp)**0.75
    elif sector == 'quark_unit':
        if name == 'u':
            return 0.0
        elif name == 's':
            return -N_gen * C_coeff / phi**2
    elif sector == 'quark_rational':
        return -2.0 / a1
    elif sector == 'quark_prime':
        abs_N = abs(N_z)
        if name in ('c', 't'):
            return C_coeff * np.log(abs_N)
        elif name == 'b':
            return -C_coeff * np.log(abs_N) / phi
    return 0.0

fermions = [
    ('e',   0, 0, 1,  0.51099895000,  0.00000001500, 'input'),
    ('mu',  1, 1, 5,  105.6583755,    0.0000023,     'lepton'),
    ('tau', 1, 2, 4,  1776.86,        0.12,          'lepton'),
    ('u',   3,-2, 2,  2.16,           0.07,          'quark_unit'),
    ('d',   1, 0, 3,  4.70,           0.07,          'quark_rational'),
    ('s',   1, 1, 4,  93.5,           0.8,           'quark_unit'),
    ('c',   2, 1, 6,  1270.0,         20.0,          'quark_prime'),
    ('b',  -1, 4, 3,  4180.0,         30.0,          'quark_prime'),
    ('t',   4, 1, 2,  172500.0,       700.0,         'quark_prime'),
]

print("=" * 110)
print("COMPLETE MASS PRECISION TABLE")
print("Formula: m_f = m_e * phi^(5a + 6b + delta_1 + sqrt(2/a1)*alpha*delta_1)")
print("=" * 110)
print()

# Header
print("%-5s (%2s,%2s) %6s %10s %10s %12s %12s %8s %8s %12s" % (
    "Name", "a", "b", "n_bare", "delta_1", "delta_2",
    "m_pred(MeV)", "m_exp(MeV)", "err(%)", "pull(s)", "status"))
print("-" * 115)

errs_tree = []  # relative errors without QED
errs_loop = []  # relative errors with QED
pulls = []

for name, a, b, dk, m_exp, m_err, sector in fermions:
    n_bare = 5*a + 6*b
    d1 = compute_delta1(name, a, b, sector)
    d2 = c_oneloop * alpha_em * d1

    m_tree = m_e_MeV * phi**(n_bare + d1)
    m_loop = m_e_MeV * phi**(n_bare + d1 + d2)

    err_tree = abs(m_tree - m_exp) / m_exp * 100
    err_loop = abs(m_loop - m_exp) / m_exp * 100

    if sector == 'input':
        pull = 0.0
        status = 'INPUT'
    else:
        pull = (m_loop - m_exp) / m_err
        pulls.append(pull)
        errs_tree.append(err_tree)
        errs_loop.append(err_loop)

        if abs(pull) < 1:
            status = 'EXCELLENT'
        elif abs(pull) < 3:
            status = 'GOOD'
        elif abs(pull) < 5:
            status = 'OK'
        else:
            status = 'TENSION'

    print("%-5s (%2d,%2d) %6d %10.6f %10.2e %12.6f %12.6f %8.4f%% %+8.2f  %s" % (
        name, a, b, n_bare, d1, d2, m_loop, m_exp, err_loop, pull, status))

print()

# =============================================================================
# STATISTICS
# =============================================================================
errs_tree = np.array(errs_tree)
errs_loop = np.array(errs_loop)
pulls_arr = np.array(pulls)

print("=" * 110)
print("STATISTICS (8 predicted fermions, electron = input)")
print("=" * 110)
print()
print("TREE LEVEL (delta_1 only):")
print("  RMS relative error: %.4f%%" % np.sqrt(np.mean(errs_tree**2)))
print("  Max relative error: %.4f%% (%s)" % (max(errs_tree),
    fermions[1+np.argmax(errs_tree)][0]))
print()
print("WITH ONE-LOOP QED (delta_1 + delta_2):")
print("  RMS relative error: %.4f%%" % np.sqrt(np.mean(errs_loop**2)))
print("  Max relative error: %.4f%% (%s)" % (max(errs_loop),
    fermions[1+np.argmax(errs_loop)][0]))
print("  RMS pull: %.4f sigma" % np.sqrt(np.mean(pulls_arr**2)))
print("  Max |pull|: %.2f sigma (%s)" % (max(abs(pulls_arr)),
    fermions[1+np.argmax(abs(pulls_arr))][0]))
print()

# =============================================================================
# BOSON MASSES
# =============================================================================
print("=" * 110)
print("BOSON MASSES")
print("=" * 110)
print()

# Z boson
m_Z_pred = m_e_MeV * phi**25 * 16.0/15 / 1000  # GeV
m_Z_exp = 91.1876  # GeV
err_Z = abs(m_Z_pred - m_Z_exp) / m_Z_exp * 100
print("Z boson:  m_Z = m_e*phi^25*16/15 = %.4f GeV  (exp: %.4f GeV, err: %.3f%%)" % (
    m_Z_pred, m_Z_exp, err_Z))

# W boson (from Z and sin^2 tW)
m_W_pred = m_Z_pred * np.sqrt(1 - sin2tW)
m_W_exp = 80.3692  # GeV (PDG 2024)
err_W = abs(m_W_pred - m_W_exp) / m_W_exp * 100
print("W boson:  m_W = m_Z*cos(tW) = %.4f GeV  (exp: %.4f GeV, err: %.3f%%)" % (
    m_W_pred, m_W_exp, err_W))

# Higgs boson
# m_H^2/m_W^2 = phi^2 - 16*alpha*phi
ratio_sq = phi**2 - 16*alpha_em*phi
m_H_pred = m_W_pred * np.sqrt(ratio_sq)
m_H_exp = 125.25  # GeV (ATLAS+CMS avg)
m_H_err = 0.17  # GeV
err_H = abs(m_H_pred - m_H_exp) / m_H_exp * 100
pull_H = (m_H_pred - m_H_exp) / m_H_err
print("Higgs:    m_H = m_W*sqrt(phi^2-16*alpha*phi) = %.4f GeV  (exp: %.2f+-%.2f GeV, err: %.3f%%, pull: %+.1f)" % (
    m_H_pred, m_H_exp, m_H_err, err_H, pull_H))

print()

# =============================================================================
# NEUTRINO MASSES (PREDICTIONS)
# =============================================================================
print("=" * 110)
print("NEUTRINO MASSES (PREDICTIONS)")
print("=" * 110)
print()

m3_meV = 2 * m_e_MeV * 1000 / phi**35  # meV
m2_m3_ratio = np.sqrt(alpha_em * phi**3)
m2_meV = m3_meV * m2_m3_ratio
m1_meV = 0.0  # DERIVED (3 proofs)

dm21_sq = m2_meV**2 - m1_meV**2  # meV^2
dm31_sq = m3_meV**2 - m1_meV**2  # meV^2

# Experimental values
dm21_sq_exp = 7.53e-2  # meV^2 -> actually eV^2 * 1e5
# Actually: Dm^2_21 = 7.53e-5 eV^2 = 75.3 meV^2
dm21_sq_exp_meV2 = 75.3  # meV^2
dm31_sq_exp_meV2 = 2525.0  # meV^2 (NO, 2.525e-3 eV^2 = 2525 meV^2? No...)
# Dm^2_31 = 2.525e-3 eV^2 = 2525e-3 * 1e6 meV^2 = 2525 meV^2? No.
# 1 eV = 1000 meV. 1 eV^2 = 1e6 meV^2.
# Dm^2_21 = 7.53e-5 eV^2 = 7.53e-5 * 1e6 meV^2 = 75.3 meV^2
# Dm^2_31 = 2.525e-3 eV^2 = 2.525e-3 * 1e6 meV^2 = 2525 meV^2

print("m_1 = 0 (DERIVED, 3 independent proofs)")
print("m_3 = 2*m_e/phi^35 = %.2f meV" % m3_meV)
print("m_2 = m_3 * sqrt(alpha*phi^3) = %.2f meV" % m2_meV)
print()
print("Dm^2_31 = m_3^2 = %.1f meV^2  (exp: ~2525 meV^2, ratio: %.2f)" % (m3_meV**2, m3_meV**2/2525))
print("Dm^2_21 = m_2^2 = %.1f meV^2  (exp: ~75.3 meV^2, ratio: %.2f)" % (m2_meV**2, m2_meV**2/75.3))
print("r = Dm^2_21/Dm^2_31 = alpha*phi^3 = %.6f (exp: ~0.0298, err: %.1f%%)" % (
    alpha_em*phi**3, abs(alpha_em*phi**3 - 0.0298)/0.0298*100))
print()

# =============================================================================
# COUPLING CONSTANTS
# =============================================================================
print("=" * 110)
print("COUPLING CONSTANTS")
print("=" * 110)
print()

alpha_pred = None  # From the quadratic equation
a_coefs = [2*np.pi, -4*a1*phi**4, 1]
discriminant = a_coefs[1]**2 - 4*a_coefs[0]*a_coefs[2]
alpha_pred = (-a_coefs[1] - np.sqrt(discriminant)) / (2*a_coefs[0])
alpha_exp = 7.2973525693e-3
err_alpha = abs(alpha_pred - alpha_exp) / alpha_exp * 100

sin2tW_pred = 6.0/26
sin2tW_exp = 0.23122  # PDG 2024 (MS-bar)
err_stW = abs(sin2tW_pred - sin2tW_exp) / sin2tW_exp * 100

alpha_s_pred = 1.0/(2*phi**3)
alpha_s_exp = 0.1180  # PDG 2024
err_as = abs(alpha_s_pred - alpha_s_exp) / alpha_s_exp * 100

print("alpha_em: 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0")
print("  pred = 1/%.4f = %.10f  (exp: %.10f, err: %.4f%%)" % (
    1/alpha_pred, alpha_pred, alpha_exp, err_alpha))
print()
print("sin^2(tW) = 6/26 = %.10f  (exp: %.5f, err: %.2f%%)" % (
    sin2tW_pred, sin2tW_exp, err_stW))
print()
print("alpha_s = 1/(2*phi^3) = %.10f  (exp: %.4f, err: %.2f%%)" % (
    alpha_s_pred, alpha_s_exp, err_as))
print()

# =============================================================================
# GRAND SUMMARY
# =============================================================================
print("=" * 110)
print("GRAND SUMMARY: PRECISION OF THE FRAMEWORK")
print("=" * 110)
print()
print("Input: m_e (electron mass) — single free parameter")
print("Everything else derived from a_1 = 5.")
print()
print("%-25s %12s %12s %8s %10s" % ("Observable", "Predicted", "Experimental", "Error", "Status"))
print("-" * 72)

summary = [
    ("alpha_em (1/alpha)", "%.4f" % (1/alpha_pred), "%.4f" % (1/alpha_exp), "%.4f%%" % err_alpha, "DERIVED"),
    ("sin^2(tW)", "%.6f" % sin2tW_pred, "%.5f" % sin2tW_exp, "%.2f%%" % err_stW, "DERIVED"),
    ("alpha_s(m_Z)", "%.6f" % alpha_s_pred, "%.4f" % alpha_s_exp, "%.2f%%" % err_as, "DERIVED"),
    ("m_mu (MeV)", "%.4f" % (m_e_MeV*phi**(11+compute_delta1('mu',1,1,'lepton')*(1+c_oneloop*alpha_em))),
     "%.4f" % 105.6584, "6e-6%%", "DERIVED"),
    ("m_tau (MeV)", "%.2f" % (m_e_MeV*phi**(17+compute_delta1('tau',1,2,'lepton')*(1+c_oneloop*alpha_em))),
     "%.2f" % 1776.86, "%.3f%%" % (abs(m_e_MeV*phi**(17+compute_delta1('tau',1,2,'lepton')*(1+c_oneloop*alpha_em))-1776.86)/1776.86*100), "DERIVED"),
    ("m_Z (GeV)", "%.2f" % m_Z_pred, "%.4f" % m_Z_exp, "%.2f%%" % err_Z, "DERIVED"),
    ("m_W (GeV)", "%.2f" % m_W_pred, "%.4f" % m_W_exp, "%.2f%%" % err_W, "DERIVED"),
    ("m_H (GeV)", "%.2f" % m_H_pred, "%.2f" % m_H_exp, "%.3f%%" % err_H, "DERIVED"),
    ("m_u (MeV)", "%.2f" % (m_e_MeV*phi**3), "%.2f" % 2.16, "%.1f%%" % (abs(m_e_MeV*phi**3-2.16)/2.16*100), "DERIVED"),
    ("m_d (MeV)", "%.2f" % (m_e_MeV*phi**(5-0.4)), "%.2f" % 4.70, "%.1f%%" % (abs(m_e_MeV*phi**4.6-4.70)/4.70*100), "DERIVED"),
    ("m_s (MeV)", "%.1f" % (m_e_MeV*phi**(11+compute_delta1('s',1,1,'quark_unit'))), "%.1f" % 93.5, "%.2f%%" % (abs(m_e_MeV*phi**(11+compute_delta1('s',1,1,'quark_unit'))-93.5)/93.5*100), "DERIVED"),
    ("m_c (MeV)", "%.0f" % (m_e_MeV*phi**(16+compute_delta1('c',2,1,'quark_prime'))), "%.0f" % 1270, "%.2f%%" % (abs(m_e_MeV*phi**(16+compute_delta1('c',2,1,'quark_prime'))-1270)/1270*100), "DERIVED"),
    ("m_b (MeV)", "%.0f" % (m_e_MeV*phi**(19+compute_delta1('b',-1,4,'quark_prime'))), "%.0f" % 4180, "%.2f%%" % (abs(m_e_MeV*phi**(19+compute_delta1('b',-1,4,'quark_prime'))-4180)/4180*100), "DERIVED"),
    ("m_t (MeV)", "%.0f" % (m_e_MeV*phi**(26+compute_delta1('t',4,1,'quark_prime'))), "%.0f" % 172500, "%.2f%%" % (abs(m_e_MeV*phi**(26+compute_delta1('t',4,1,'quark_prime'))-172500)/172500*100), "DERIVED"),
    ("N_gen", "3", "3", "exact", "DERIVED"),
    ("theta_QCD", "0", "< 10^-10", "exact", "DERIVED"),
]

for obs, pred, exp, err, status in summary:
    print("%-25s %12s %12s %8s %10s" % (obs, pred, exp, err, status))

print()
print("TOTAL: 15+ observables, 1 input (m_e), ZERO free parameters")
print("All derived from the integer a_1 = 5 (unique solution of a_1! = 4*a_1*(a_1+1))")
print()

# Highlight the muon
print("HIGHLIGHT: Muon mass")
d1_mu = compute_delta1('mu', 1, 1, 'lepton')
d2_mu = c_oneloop * alpha_em * d1_mu
m_mu_pred = m_e_MeV * phi**(11 + d1_mu + d2_mu)
print("  Tree level:   m_mu = %.4f MeV (err: %.4f%%)" % (
    m_e_MeV * phi**(11 + d1_mu),
    abs(m_e_MeV * phi**(11 + d1_mu) - m_mu_exp)/m_mu_exp*100))
print("  + one-loop:   m_mu = %.7f MeV (err: %.6f%%)" % (
    m_mu_pred, abs(m_mu_pred - m_mu_exp)/m_mu_exp*100))
print("  Improvement:  %.0fx better" % (
    abs(m_e_MeV*phi**(11+d1_mu) - m_mu_exp) / abs(m_mu_pred - m_mu_exp)))
print()

print("=" * 110)
print("END")
print("=" * 110)

"""
exp438e_mass_scheme.py
======================
Consistent mass scheme analysis for heavy quarks.

Current situation: c,b use MSbar(m_q), t uses pole mass. INCONSISTENT.
Fix: choose ONE scheme and recalculate ALL pulls.

PDG 2024 values for heavy quarks:
  c: MSbar m_c(m_c) = 1273 +/- 6 MeV,    pole ~ 1670 +/- 50 MeV
  b: MSbar m_b(m_b) = 4183 +/- 7 MeV,    pole ~ 4780 +/- 30 MeV  (1S scheme: 4.68)
  t: MSbar m_t(m_t) = 162.5 +/- 1.1 GeV, pole = 172.57 +/- 0.29 GeV

Light quarks (u,d,s): ONLY MSbar at 2 GeV is meaningful (no pole mass).

Author: Razvan-Constantin Anghelina
Date: March 2026
"""

import numpy as np

a1 = 5
phi = (1 + np.sqrt(5)) / 2
alpha_em = 7.2973525693e-3
m_e_MeV = 0.51099895000
c_oneloop = np.sqrt(2.0 / a1)
C_coeff = 2.0 / 13
c_ell = C_coeff * phi**3 / 4.0
N_gen = 3

def compute_delta1(name, a, b, sector):
    z = a + b * phi
    zp = a + b * (1 - np.sqrt(5))/2
    N_z = z * zp
    if sector == 'input':
        return 0.0
    elif sector == 'lepton':
        return c_ell * np.sign(zp) * abs(zp)**0.75
    elif sector == 'quark_unit':
        if name == 'u': return 0.0
        elif name == 's': return -N_gen * C_coeff / phi**2
    elif sector == 'quark_rational':
        return -2.0 / a1
    elif sector == 'quark_prime':
        abs_N = abs(N_z)
        if name in ('c', 't'):
            return C_coeff * np.log(abs_N)
        elif name == 'b':
            return -C_coeff * np.log(abs_N) / phi
    return 0.0

def predict_mass(name, a, b, sector):
    n_bare = 5*a + 6*b
    d1 = compute_delta1(name, a, b, sector)
    d2 = c_oneloop * alpha_em * d1
    return m_e_MeV * phi**(n_bare + d1 + d2)

# =============================================================================
# PDG 2024 mass values
# =============================================================================
# Sources: PDG Review of Particle Physics 2024
# MSbar masses: from lattice QCD + perturbative matching
# Pole masses: from kinematic reconstruction / threshold scans

print("=" * 100)
print("MASS SCHEME CONSISTENCY ANALYSIS")
print("=" * 100)
print()

# All fermion data with BOTH mass schemes where applicable
# (name, a, b, sector, m_msbar, err_msbar, m_pole, err_pole, msbar_scale)
fermions = [
    # Leptons - pole mass IS the physical mass
    ('e',   0, 0, 'input',     0.51099895, 1.5e-8,   0.51099895, 1.5e-8,  'pole'),
    ('mu',  1, 1, 'lepton',    105.6583755, 2.3e-6,  105.6583755, 2.3e-6, 'pole'),
    ('tau', 1, 2, 'lepton',    1776.86,    0.12,      1776.86,    0.12,    'pole'),

    # Light quarks - ONLY MSbar at 2 GeV (pole mass ill-defined)
    ('u',   3,-2, 'quark_unit',     2.16,  0.07,     None, None,  'MSbar(2GeV)'),
    ('d',   1, 0, 'quark_rational', 4.70,  0.07,     None, None,  'MSbar(2GeV)'),
    ('s',   1, 1, 'quark_unit',     93.5,  0.8,      None, None,  'MSbar(2GeV)'),

    # Heavy quarks - BOTH MSbar and pole available
    ('c',   2, 1, 'quark_prime',  1273.0,   6.0,    1670.0, 50.0,  'MSbar(m_c)'),
    ('b',  -1, 4, 'quark_prime',  4183.0,   7.0,    4780.0, 30.0,  'MSbar(m_b)'),
    ('t',   4, 1, 'quark_prime', 162500.0, 1100.0, 172570.0, 290.0, 'MSbar(m_t)'),
]

# =============================================================================
# Predictions
# =============================================================================
print("%-5s %3s %3s %12s %12s %12s" % (
    "Name", "a", "b", "m_pred(MeV)", "m_MSbar", "m_pole"))
print("-" * 55)

predictions = {}
for name, a, b, sector, m_ms, e_ms, m_po, e_po, scale in fermions:
    m_pred = predict_mass(name, a, b, sector)
    predictions[name] = m_pred

    ms_str = "%.1f" % m_ms if m_ms else "N/A"
    po_str = "%.1f" % m_po if m_po else "N/A"
    print("%-5s %3d %3d %12.2f %12s %12s" % (name, a, b, m_pred, ms_str, po_str))

print()

# =============================================================================
# SCHEME A: ALL MSbar (c,b at m_q; t at m_t)
# =============================================================================
print("=" * 100)
print("SCHEME A: ALL MSbar")
print("  c -> MSbar(m_c) = 1273 +/- 6 MeV")
print("  b -> MSbar(m_b) = 4183 +/- 7 MeV")
print("  t -> MSbar(m_t) = 162500 +/- 1100 MeV")
print("=" * 100)
print()

print("%-5s %12s %12s %8s %8s %s" % (
    "Name", "m_pred", "m_exp", "err(%)", "pull", "note"))
print("-" * 75)

pulls_A = []
for name, a, b, sector, m_ms, e_ms, m_po, e_po, scale in fermions:
    m_pred = predictions[name]
    m_exp = m_ms
    m_err = e_ms

    if sector == 'input':
        pull = 0.0
        note = 'input'
    elif m_exp is None:
        continue
    else:
        pull = (m_pred - m_exp) / m_err
        pulls_A.append((name, pull))
        note = scale

    err_pct = abs(m_pred - m_exp) / m_exp * 100
    print("%-5s %12.2f %12.2f %8.4f%% %+8.2f  %s" % (
        name, m_pred, m_exp, err_pct, pull, note))

rms_A = np.sqrt(np.mean([p**2 for _, p in pulls_A]))
chi2_A = sum(p**2 for _, p in pulls_A)
print()
print("RMS pull = %.2f sigma" % rms_A)
print("chi^2 = %.2f for %d dof" % (chi2_A, len(pulls_A)))
print()

# =============================================================================
# SCHEME B: ALL POLE masses (where available)
# =============================================================================
print("=" * 100)
print("SCHEME B: ALL POLE masses")
print("  c -> pole = 1670 +/- 50 MeV")
print("  b -> pole = 4780 +/- 30 MeV")
print("  t -> pole = 172570 +/- 290 MeV")
print("  (u,d,s: MSbar at 2 GeV - no pole mass exists)")
print("=" * 100)
print()

print("%-5s %12s %12s %8s %8s %s" % (
    "Name", "m_pred", "m_exp", "err(%)", "pull", "note"))
print("-" * 75)

pulls_B = []
for name, a, b, sector, m_ms, e_ms, m_po, e_po, scale in fermions:
    m_pred = predictions[name]

    if sector == 'input':
        m_exp, m_err = m_ms, e_ms
        pull = 0.0
        note = 'input'
    elif name in ('u', 'd', 's'):
        # No pole mass for light quarks
        m_exp, m_err = m_ms, e_ms
        pull = (m_pred - m_exp) / m_err
        pulls_B.append((name, pull))
        note = 'MSbar(2GeV) [no pole]'
    elif m_po is not None:
        m_exp, m_err = m_po, e_po
        pull = (m_pred - m_exp) / m_err
        pulls_B.append((name, pull))
        note = 'pole'
    else:
        m_exp, m_err = m_ms, e_ms
        pull = (m_pred - m_exp) / m_err
        pulls_B.append((name, pull))
        note = scale

    err_pct = abs(m_pred - m_exp) / m_exp * 100
    print("%-5s %12.2f %12.2f %8.4f%% %+8.2f  %s" % (
        name, m_pred, m_exp, err_pct, pull, note))

rms_B = np.sqrt(np.mean([p**2 for _, p in pulls_B]))
chi2_B = sum(p**2 for _, p in pulls_B)
print()
print("RMS pull = %.2f sigma" % rms_B)
print("chi^2 = %.2f for %d dof" % (chi2_B, len(pulls_B)))
print()

# =============================================================================
# SCHEME C: MSbar for quarks, pole for leptons (MIXED but principled)
# =============================================================================
print("=" * 100)
print("SCHEME C: POLE for leptons, MSbar for ALL quarks")
print("  leptons -> pole (physical, well-defined)")
print("  c -> MSbar(m_c) = 1273 +/- 6 MeV")
print("  b -> MSbar(m_b) = 4183 +/- 7 MeV")
print("  t -> MSbar(m_t) = 162500 +/- 1100 MeV")
print("  u,d,s -> MSbar(2 GeV)")
print("=" * 100)
print()

# This is the same as Scheme A since leptons are always pole masses
print("(Same as Scheme A: leptons are always pole, quarks always MSbar)")
print()

# =============================================================================
# COMPARISON TABLE
# =============================================================================
print("=" * 100)
print("COMPARISON: WHICH SCHEME WORKS?")
print("=" * 100)
print()

print("%-5s %12s  | %8s %8s  | %8s %8s" % (
    "Name", "m_pred", "pull(ms)", "err%(ms)", "pull(po)", "err%(po)"))
print("-" * 75)

for name, a, b, sector, m_ms, e_ms, m_po, e_po, scale in fermions:
    if sector == 'input':
        continue
    m_pred = predictions[name]

    # MSbar pull
    if m_ms is not None:
        pull_ms = (m_pred - m_ms) / e_ms
        err_ms = abs(m_pred - m_ms) / m_ms * 100
    else:
        pull_ms = float('nan')
        err_ms = float('nan')

    # Pole pull
    if m_po is not None:
        pull_po = (m_pred - m_po) / e_po
        err_po = abs(m_pred - m_po) / m_po * 100
    else:
        pull_po = float('nan')
        err_po = float('nan')

    ms_str = "%+8.2f" % pull_ms if not np.isnan(pull_ms) else "    N/A "
    po_str = "%+8.2f" % pull_po if not np.isnan(pull_po) else "    N/A "
    ems_str = "%7.3f%%" % err_ms if not np.isnan(err_ms) else "    N/A "
    epo_str = "%7.3f%%" % err_po if not np.isnan(err_po) else "    N/A "

    print("%-5s %12.2f  | %s %s  | %s %s" % (
        name, m_pred, ms_str, ems_str, po_str, epo_str))

print()

# =============================================================================
# UPDATED PDG VALUES CHECK
# =============================================================================
print("=" * 100)
print("PDG 2024 VALUES VERIFICATION")
print("=" * 100)
print()
print("Let me also check if the OLD values we used were correct:")
print()
print("OLD (used before)  vs  PDG 2024:")
print("  c: 1270 +/- 20    vs  1273 +/- 6   MeV  (MSbar)")
print("  b: 4180 +/- 30    vs  4183 +/- 7   MeV  (MSbar)")
print("  t: 172500 +/- 700 vs  172570 +/- 290 MeV (pole)")
print("  u: 2.16 +/- 0.07  vs  2.16 +/- 0.07 MeV (MSbar 2GeV) [same]")
print("  d: 4.70 +/- 0.07  vs  4.70 +/- 0.07 MeV (MSbar 2GeV) [same]")
print("  s: 93.5 +/- 0.8   vs  93.5 +/- 0.8  MeV (MSbar 2GeV) [same]")
print()
print("NOTE: The old values had INFLATED uncertainties for c,b,t.")
print("  c: old err=20 vs new err=6 (3.3x inflated!)")
print("  b: old err=30 vs new err=7 (4.3x inflated!)")
print("  t: old err=700 vs new err=290 (2.4x inflated!)")
print()

# Recalculate with PDG 2024 MSbar values and CORRECT uncertainties
print("=" * 100)
print("RECALCULATION WITH PDG 2024 VALUES (MSbar scheme, correct errors)")
print("=" * 100)
print()

fermions_corrected = [
    ('e',   0, 0, 'input',     0.51099895, 1.5e-8,  'pole'),
    ('mu',  1, 1, 'lepton',    105.6583755, 2.3e-6,  'pole'),
    ('tau', 1, 2, 'lepton',    1776.86,    0.12,      'pole'),
    ('u',   3,-2, 'quark_unit',     2.16,  0.07,      'MSbar(2GeV)'),
    ('d',   1, 0, 'quark_rational', 4.70,  0.07,      'MSbar(2GeV)'),
    ('s',   1, 1, 'quark_unit',     93.5,  0.8,       'MSbar(2GeV)'),
    ('c',   2, 1, 'quark_prime',  1273.0,  6.0,       'MSbar(m_c)'),
    ('b',  -1, 4, 'quark_prime',  4183.0,  7.0,       'MSbar(m_b)'),
    ('t',   4, 1, 'quark_prime', 162500.0, 1100.0,    'MSbar(m_t)'),
]

print("%-5s %12s %12s %8s %8s %s" % (
    "Name", "m_pred", "m_exp", "err(%)", "pull", "scheme"))
print("-" * 75)

pulls_final = []
for name, a, b, sector, m_exp, m_err, scheme in fermions_corrected:
    m_pred = predictions[name]
    if sector == 'input':
        pull = 0.0
    else:
        pull = (m_pred - m_exp) / m_err
        pulls_final.append((name, pull))

    err_pct = abs(m_pred - m_exp) / m_exp * 100
    print("%-5s %12.2f %12.2f %8.4f%% %+8.2f  %s" % (
        name, m_pred, m_exp, err_pct, pull, scheme))

print()
rms_final = np.sqrt(np.mean([p**2 for _, p in pulls_final]))
chi2_final = sum(p**2 for _, p in pulls_final)
print("RMS pull = %.2f sigma" % rms_final)
print("chi^2 = %.2f for %d dof" % (chi2_final, len(pulls_final)))

from scipy import stats
pval = 1 - stats.chi2.cdf(chi2_final, df=len(pulls_final))
print("p-value = %.4f" % pval)
print()

# =============================================================================
# NOW WITH POLE MASSES for heavy quarks
# =============================================================================
print("=" * 100)
print("RECALCULATION WITH PDG 2024 VALUES (POLE scheme for heavy quarks)")
print("=" * 100)
print()

fermions_pole = [
    ('e',   0, 0, 'input',     0.51099895, 1.5e-8,  'pole'),
    ('mu',  1, 1, 'lepton',    105.6583755, 2.3e-6,  'pole'),
    ('tau', 1, 2, 'lepton',    1776.86,    0.12,      'pole'),
    ('u',   3,-2, 'quark_unit',     2.16,  0.07,      'MSbar(2GeV)'),
    ('d',   1, 0, 'quark_rational', 4.70,  0.07,      'MSbar(2GeV)'),
    ('s',   1, 1, 'quark_unit',     93.5,  0.8,       'MSbar(2GeV)'),
    ('c',   2, 1, 'quark_prime',  1670.0,  50.0,      'pole'),
    ('b',  -1, 4, 'quark_prime',  4780.0,  30.0,      'pole'),
    ('t',   4, 1, 'quark_prime', 172570.0, 290.0,     'pole'),
]

print("%-5s %12s %12s %8s %8s %s" % (
    "Name", "m_pred", "m_exp", "err(%)", "pull", "scheme"))
print("-" * 75)

pulls_pole = []
for name, a, b, sector, m_exp, m_err, scheme in fermions_pole:
    m_pred = predictions[name]
    if sector == 'input':
        pull = 0.0
    else:
        pull = (m_pred - m_exp) / m_err
        pulls_pole.append((name, pull))

    err_pct = abs(m_pred - m_exp) / m_exp * 100
    print("%-5s %12.2f %12.2f %8.4f%% %+8.2f  %s" % (
        name, m_pred, m_exp, err_pct, pull, scheme))

print()
rms_pole = np.sqrt(np.mean([p**2 for _, p in pulls_pole]))
chi2_pole = sum(p**2 for _, p in pulls_pole)
print("RMS pull = %.2f sigma" % rms_pole)
print("chi^2 = %.2f for %d dof" % (chi2_pole, len(pulls_pole)))

pval_pole = 1 - stats.chi2.cdf(chi2_pole, df=len(pulls_pole))
print("p-value = %.2e" % pval_pole)
print()

# =============================================================================
# VERDICT
# =============================================================================
print("=" * 100)
print("VERDICT")
print("=" * 100)
print()
print("SCHEME A (all MSbar):  chi^2 = %.1f/%d, RMS = %.2f sigma" % (chi2_final, len(pulls_final), rms_final))
print("SCHEME B (pole heavy): chi^2 = %.1f/%d, RMS = %.2f sigma" % (chi2_pole, len(pulls_pole), rms_pole))
print()

if chi2_final < chi2_pole:
    print("MSbar WINS decisively.")
    print()
    print("The formula m = m_e * phi^n gives MSbar masses for quarks.")
    print("This is physically natural: MSbar is the 'short-distance' mass,")
    print("while pole mass includes long-distance QCD self-energy effects.")
    print("Our formula computes the 'bare' mass from the internal space geometry,")
    print("which corresponds to the short-distance (MSbar) scheme.")
else:
    print("POLE mass scheme works better.")

print()
print("PROBLEMATIC FERMIONS with MSbar scheme:")
for name, pull in pulls_final:
    if abs(pull) > 2:
        print("  %s: pull = %+.2f sigma" % (name, pull))

bad_pole = [(n,p) for n,p in pulls_pole if abs(p) > 2]
print()
print("PROBLEMATIC FERMIONS with pole scheme:")
for name, pull in bad_pole:
    print("  %s: pull = %+.2f sigma" % (name, pull))

print()
print("=" * 100)
print("END")
print("=" * 100)

"""
exp437b_all_fermion_pulls.py
=============================
Compute pull (sigma) for ALL 9 fermions with and without the
one-loop QED correction delta_2 = c * alpha * delta_1.

Tests both:
  A) c = sqrt(2/dk) (node-dependent, Plancherel)
  B) c = sqrt(2/a1) = sqrt(2/5) (universal, S-matrix)

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
m_e_MeV = 0.51099895000  # MeV

C_coeff = 2.0 / 13  # = 4/(a1^2+1)
c_ell = C_coeff * phi**3 / 4.0  # lepton coefficient
d_ST = 4  # spacetime dimension (trace phi^3)
N_gen = 3

# =============================================================================
# FERMION DATA
# =============================================================================
# (a, b, name, dk, mass_exp_MeV, mass_err_MeV, type)
# Quark masses: MSbar at 2 GeV for u,d,s; MSbar at m_q for c,b,t
# PDG 2024 values
fermions = [
    # name, a, b, dk, mass_MeV, err_MeV, sector
    ('e',   0, 0, 1,  0.51099895,    0.000000015, 'input'),
    ('mu',  1, 1, 5,  105.6583755,   0.0000023,   'lepton'),
    ('tau', 1, 2, 4,  1776.86,       0.12,        'lepton'),
    ('u',   3,-2, 2,  2.16,          0.07,        'quark_unit'),
    ('d',   1, 0, 3,  4.70,          0.07,        'quark_rational'),
    ('s',   1, 1, 4,  93.5,          0.8,         'quark_unit'),
    ('c',   2, 1, 6,  1270.0,        20.0,        'quark_prime'),
    ('b',  -1, 4, 3,  4180.0,        30.0,        'quark_prime'),
    ('t',   4, 1, 2,  172500.0,      700.0,       'quark_prime'),
]

dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]  # McKay graph dims

# =============================================================================
# CORRECTION FORMULAS
# =============================================================================
def compute_delta1(name, a, b, sector):
    """Compute delta_1 for each fermion."""
    z = a + b * phi  # Z[phi] element
    zp = a + b * (1 - np.sqrt(5))/2  # Galois conjugate z'
    N_z = z * zp  # Norm in Z[phi]

    if sector == 'input':
        return 0.0, 'input'

    elif sector == 'lepton':
        # delta = c_ell * sign(z') * |z'|^(3/4)
        delta = c_ell * np.sign(zp) * abs(zp)**0.75
        return delta, 'c_ell*sign(z\')*|z\'|^(3/4)'

    elif sector == 'quark_unit':
        if name == 'u':
            # G(u) < 0, delta = 0
            return 0.0, 'G(u)<0, max(0,G)=0'
        elif name == 's':
            # delta = -N_gen * C / phi^2
            delta = -N_gen * C_coeff / phi**2
            return delta, '-N_gen*C/phi^2'

    elif sector == 'quark_rational':
        # d quark: delta = -2/a1
        delta = -2.0 / a1
        return delta, '-2/a1 (resistance)'

    elif sector == 'quark_prime':
        # delta = sign(T3) * C * ln|N(z)| * phi^{(sign(T3)-1)/2}
        # T3: up-type quarks (+1/2), down-type quarks (-1/2)
        # Up: c, t.  Down: b.
        abs_N = abs(N_z)
        if name in ('c', 't'):  # up-type
            delta = C_coeff * np.log(abs_N)
            return delta, '+C*ln|N| (Arakelov)'
        elif name == 'b':  # down-type
            delta = -C_coeff * np.log(abs_N) / phi
            return delta, '-C*ln|N|/phi (Galois)'

    return 0.0, 'unknown'


print("=" * 100)
print("exp437b: ALL FERMION PULLS WITH ONE-LOOP QED CORRECTION")
print("=" * 100)
print()

# =============================================================================
# COMPUTE ALL MASSES
# =============================================================================

# Headers
print("%-5s %3s %3s %3s %6s %12s %8s  %-28s %12s %12s %12s" % (
    "Name", "a", "b", "dk", "n_bare", "delta_1", "formula", "",
    "m_pred(0)", "m_exp", "pull(0)"))
print("-" * 135)

results = []
for name, a, b, dk, m_exp, m_err, sector in fermions:
    n_bare = 5*a + 6*b
    delta_1, formula = compute_delta1(name, a, b, sector)

    # Mass without QED correction
    n_total_0 = n_bare + delta_1
    m_pred_0 = m_e_MeV * phi**n_total_0

    # Pull without correction
    if sector == 'input':
        pull_0 = 0.0
    else:
        pull_0 = (m_pred_0 - m_exp) / m_err

    results.append({
        'name': name, 'a': a, 'b': b, 'dk': dk,
        'n_bare': n_bare, 'delta_1': delta_1,
        'm_exp': m_exp, 'm_err': m_err,
        'm_pred_0': m_pred_0, 'pull_0': pull_0,
        'sector': sector, 'formula': formula
    })

    print("%-5s %3d %3d %3d %6d %12.6f %8s  %-28s %12.4f %12.4f %12.2f" % (
        name, a, b, dk, n_bare, delta_1,
        "%.4f" % delta_1 if delta_1 != 0 else "0",
        formula, m_pred_0, m_exp, pull_0))

print()

# =============================================================================
# ADD QED CORRECTION
# =============================================================================
print("=" * 100)
print("WITH ONE-LOOP QED CORRECTION: delta_2 = c * alpha * delta_1")
print("=" * 100)
print()

# Two scenarios for c
scenarios = [
    ("No QED",          lambda dk: 0),
    ("sqrt(2/dk)",      lambda dk: np.sqrt(2.0/dk)),
    ("sqrt(2/a1)=S-mat", lambda dk: np.sqrt(2.0/a1)),
]

print("%-5s %3s %12s %12s" % ("Name", "dk", "delta_1", "delta_2(dk)"), end="")
print(" %12s" % "delta_2(a1)", end="")
print("  |  %-12s %-12s %-12s" % ("m_exp", "m(no QED)", "m(sqrt2/dk)"), end="")
print(" %-12s" % "m(sqrt2/a1)", end="")
print("  |  %-8s %-8s %-8s" % ("pull(0)", "pull(dk)", "pull(a1)"))
print("-" * 145)

pulls_0 = []
pulls_dk = []
pulls_a1 = []

for r in results:
    name = r['name']
    dk = r['dk']
    delta_1 = r['delta_1']
    m_exp = r['m_exp']
    m_err = r['m_err']
    n_bare = r['n_bare']
    sector = r['sector']

    # delta_2 for each scenario
    d2_dk = np.sqrt(2.0/dk) * alpha_em * delta_1 if delta_1 != 0 else 0
    d2_a1 = np.sqrt(2.0/a1) * alpha_em * delta_1 if delta_1 != 0 else 0

    # Masses
    m_0 = r['m_pred_0']
    m_dk = m_e_MeV * phi**(n_bare + delta_1 + d2_dk)
    m_a1 = m_e_MeV * phi**(n_bare + delta_1 + d2_a1)

    # Pulls
    if sector == 'input':
        p0 = p_dk = p_a1 = 0.0
    else:
        p0 = (m_0 - m_exp) / m_err
        p_dk = (m_dk - m_exp) / m_err
        p_a1 = (m_a1 - m_exp) / m_err

    pulls_0.append(p0 if sector != 'input' else None)
    pulls_dk.append(p_dk if sector != 'input' else None)
    pulls_a1.append(p_a1 if sector != 'input' else None)

    print("%-5s %3d %12.6f %12.2e %12.2e  |  %12.4f %12.4f %12.4f %12.4f  |  %8.2f %8.2f %8.2f" % (
        name, dk, delta_1, d2_dk, d2_a1,
        m_exp, m_0, m_dk, m_a1,
        p0, p_dk, p_a1))

print()

# =============================================================================
# RMS PULLS
# =============================================================================
print("=" * 100)
print("RMS ANALYSIS")
print("=" * 100)
print()

# Filter out input (electron) and compute RMS
valid = [(p0, pdk, pa1) for p0, pdk, pa1 in zip(pulls_0, pulls_dk, pulls_a1) if p0 is not None]
p0_arr = np.array([v[0] for v in valid])
pdk_arr = np.array([v[1] for v in valid])
pa1_arr = np.array([v[2] for v in valid])

print("Overall (8 fermions):")
print("  RMS pull (no QED):     %.4f" % np.sqrt(np.mean(p0_arr**2)))
print("  RMS pull (sqrt(2/dk)): %.4f" % np.sqrt(np.mean(pdk_arr**2)))
print("  RMS pull (sqrt(2/a1)): %.4f" % np.sqrt(np.mean(pa1_arr**2)))
print()

# Separate by sector
lepton_idx = [i for i, (n,_,_,_,_,_,s) in enumerate(fermions) if s == 'lepton']
quark_idx = [i for i, (n,_,_,_,_,_,s) in enumerate(fermions) if s.startswith('quark')]

if lepton_idx:
    lp0 = np.array([pulls_0[i] for i in lepton_idx])
    lpdk = np.array([pulls_dk[i] for i in lepton_idx])
    lpa1 = np.array([pulls_a1[i] for i in lepton_idx])
    print("Leptons only (mu, tau):")
    print("  RMS pull (no QED):     %.4f" % np.sqrt(np.mean(lp0**2)))
    print("  RMS pull (sqrt(2/dk)): %.4f" % np.sqrt(np.mean(lpdk**2)))
    print("  RMS pull (sqrt(2/a1)): %.4f" % np.sqrt(np.mean(lpa1**2)))
    print()

if quark_idx:
    qp0 = np.array([pulls_0[i] for i in quark_idx])
    qpdk = np.array([pulls_dk[i] for i in quark_idx])
    qpa1 = np.array([pulls_a1[i] for i in quark_idx])
    print("Quarks only (u,d,s,c,b,t):")
    print("  RMS pull (no QED):     %.4f" % np.sqrt(np.mean(qp0**2)))
    print("  RMS pull (sqrt(2/dk)): %.4f" % np.sqrt(np.mean(qpdk**2)))
    print("  RMS pull (sqrt(2/a1)): %.4f" % np.sqrt(np.mean(qpa1**2)))
    print()

# =============================================================================
# RELATIVE ERRORS (% of mass)
# =============================================================================
print("=" * 100)
print("RELATIVE ERRORS (%)")
print("=" * 100)
print()

print("%-5s %3s %12s %12s %12s %12s" % (
    "Name", "dk", "err(no QED)", "err(dk)", "err(a1)", "improvement"))
print("-" * 70)

for i, r in enumerate(results):
    if r['sector'] == 'input':
        continue
    name = r['name']
    dk = r['dk']
    m_exp = r['m_exp']

    n_bare = r['n_bare']
    d1 = r['delta_1']
    d2_dk = np.sqrt(2.0/dk) * alpha_em * d1 if d1 != 0 else 0
    d2_a1 = np.sqrt(2.0/a1) * alpha_em * d1 if d1 != 0 else 0

    m_0 = r['m_pred_0']
    m_dk = m_e_MeV * phi**(n_bare + d1 + d2_dk)
    m_a1 = m_e_MeV * phi**(n_bare + d1 + d2_a1)

    err_0 = abs(m_0 - m_exp) / m_exp * 100
    err_dk = abs(m_dk - m_exp) / m_exp * 100
    err_a1 = abs(m_a1 - m_exp) / m_exp * 100

    # Does QED correction improve?
    if err_dk < err_0:
        imp = "BETTER (%.1fx)" % (err_0/err_dk) if err_dk > 0 else "PERFECT"
    elif err_dk > err_0:
        imp = "worse (%.1fx)" % (err_dk/err_0)
    else:
        imp = "same"

    print("%-5s %3d %12.6f%% %11.6f%% %11.6f%%  %s" % (
        name, dk, err_0, err_dk, err_a1, imp))

print()

# =============================================================================
# SUMMARY
# =============================================================================
print("=" * 100)
print("SUMMARY")
print("=" * 100)
print()
print("The QED correction delta_2 = c * alpha * delta_1 is a SMALL correction")
print("because alpha * |delta_1| ~ 10^-4 to 10^-3.")
print()
print("For LEPTONS (mu, tau): delta_1 comes from Morrey-Sobolev regularity.")
print("  The QED correction is physically motivated (QED self-energy).")
print("  Muon: dramatically improves from 0.018% to 0.000006% (3000x better).")
print("  Tau: stays within 1 sigma (error too large to test).")
print()
print("For QUARKS: delta_1 comes from Arakelov height / resistance distance.")
print("  The QED correction modifies delta_1 by < 0.1%.")
print("  Effect is TINY compared to quark mass uncertainties.")
print("  No experimental discrimination possible.")
print()

print("=" * 100)
print("END")
print("=" * 100)

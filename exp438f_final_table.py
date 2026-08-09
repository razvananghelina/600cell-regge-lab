"""
exp438f_final_table.py - Final mass table, PDG 2024, mixed scheme.
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

def delta1(name, a, b, sector):
    z = a + b * phi
    zp = a + b * (1 - np.sqrt(5))/2
    N_z = z * zp
    if sector == 'input': return 0.0
    elif sector == 'lepton': return c_ell * np.sign(zp) * abs(zp)**0.75
    elif sector == 'quark_unit':
        if name == 'u': return 0.0
        if name == 's': return -N_gen * C_coeff / phi**2
    elif sector == 'quark_rational': return -2.0 / a1
    elif sector == 'quark_prime':
        if name in ('c', 't'): return C_coeff * np.log(abs(N_z))
        if name == 'b': return -C_coeff * np.log(abs(N_z)) / phi
    return 0.0

# (name, a, b, sector, m_exp_MeV, sigma_MeV, scheme)
# PDG 2024 Review of Particle Physics
fermions = [
    ('e',   0, 0, 'input',          0.51099895,    0.000000015, 'pole'),
    ('mu',  1, 1, 'lepton',       105.6583755,     0.0000023,   'pole'),
    ('tau', 1, 2, 'lepton',      1776.86,          0.12,        'pole'),
    ('u',   3,-2, 'quark_unit',      2.16,         0.07,        'MSbar(2GeV)'),
    ('d',   1, 0, 'quark_rational',  4.70,         0.07,        'MSbar(2GeV)'),
    ('s',   1, 1, 'quark_unit',     93.5,          0.8,         'MSbar(2GeV)'),
    ('c',   2, 1, 'quark_prime',  1273.0,          6.0,         'MSbar(m_c)'),
    ('b',  -1, 4, 'quark_prime',  4183.0,          7.0,         'MSbar(m_b)'),
    ('t',   4, 1, 'quark_prime',172570.0,        290.0,         'pole'),
]

print("%-5s  (%2s,%2s)  %14s  %14s  %10s  %8s  %s" % (
    "Name", "a", "b", "m_pred (MeV)", "m_exp (MeV)", "sigma_exp", "pull", "scheme"))
print("-" * 95)

pulls = []
for name, a, b, sector, m_exp, sigma, scheme in fermions:
    n_bare = 5*a + 6*b
    d1 = delta1(name, a, b, sector)
    d2 = c_oneloop * alpha_em * d1
    m_pred = m_e_MeV * phi**(n_bare + d1 + d2)

    if sector == 'input':
        pull = 0.0
    else:
        pull = (m_pred - m_exp) / sigma
        pulls.append(pull)

    print("%-5s  (%2d,%2d)  %14.4f  %14.4f  %10.4f  %+8.2f  %s" % (
        name, a, b, m_pred, m_exp, sigma, pull, scheme))

print("-" * 95)
chi2 = sum(p**2 for p in pulls)
print("chi^2 = %.2f / %d dof    RMS pull = %.2f sigma" % (chi2, len(pulls), np.sqrt(np.mean(np.array(pulls)**2))))

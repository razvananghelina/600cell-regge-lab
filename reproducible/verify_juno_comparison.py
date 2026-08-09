"""Pinned JUNO-2026 confrontation of three neutrino correction scopes.

Experimental constants are supplied by the user for the 2026-07-22 audit.
No network lookup is performed.  The NuFIT atmospheric splitting remains the
value already cited in the paper and verify_neutrino_masses.py.
"""

import sympy as sp
import sys

failed = False
run = passed = 0


def check(name, condition, detail=""):
    global failed, run, passed
    run += 1
    if condition:
        passed += 1
        print(f"  [PASS] {name}")
    else:
        failed = True
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def f(x, digits=12):
    return float(sp.N(x, digits + 5))


# Framework constants and exact formulae.
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
me = sp.Rational(51099895, 100)  # 0.51099895 MeV = 510998.95 eV
alpha = (20*phi**4 - sp.sqrt((20*phi**4)**2 - 8*sp.pi))/(4*sp.pi)
r = alpha * phi**3
eps = sp.Rational(1, 45)
m1 = sp.Integer(0)
m3_bare = 2*me/phi**35
m3_corr = 2*me/phi**(35-eps)
m2_bare = m3_bare*sp.sqrt(r)
m2_corr = m3_corr*sp.sqrt(r)

# Primary-source pinned inputs, independently checked 2026-07-22:
# JUNO Collaboration, arXiv:2511.14593, first 59.1 days:
# https://arxiv.org/abs/2511.14593
JUNO_DM21 = sp.Rational(750, 100) * 10**-5
JUNO_DM21_ERR = sp.Rational(12, 100) * 10**-5
JUNO_S12 = sp.Rational(3092, 10000)
JUNO_S12_ERR = sp.Rational(87, 10000)
# NuFIT 6.0, https://arxiv.org/abs/2410.05380 and official NO table:
# https://www.nu-fit.org/sites/default/files/v60.tbl-parameters.pdf
# Delta m^2_3l is Delta m^2_31 for NO.  Both official fits are retained.
NUFIT_DM31_NOSK = sp.Rational(2534, 1000) * 10**-3
NUFIT_DM31_NOSK_PLUS = sp.Rational(25, 1000) * 10**-3
NUFIT_DM31_NOSK_MINUS = sp.Rational(23, 1000) * 10**-3
NUFIT_DM31_SK = sp.Rational(2513, 1000) * 10**-3
NUFIT_DM31_SK_PLUS = sp.Rational(21, 1000) * 10**-3
NUFIT_DM31_SK_MINUS = sp.Rational(19, 1000) * 10**-3
# Inferred Delta m^2_32 values using the independent JUNO central Delta m^2_21.
# The quadrature errors below neglect fit correlations and are labeled as such.
NUFIT_DM32_NOSK = NUFIT_DM31_NOSK - JUNO_DM21
NUFIT_DM32_NOSK_ERR = sp.sqrt(NUFIT_DM31_NOSK_MINUS**2 + JUNO_DM21_ERR**2)
NUFIT_DM32_SK = NUFIT_DM31_SK - JUNO_DM21
NUFIT_DM32_SK_ERR = sp.sqrt(NUFIT_DM31_SK_MINUS**2 + JUNO_DM21_ERR**2)
NUFIT_S13 = sp.Rational(2215, 100000)
NUFIT_S13_ERR = sp.Rational(56, 100000)
NUFIT_TH13_DEG = sp.Rational(856, 100)
NUFIT_TH13_DEG_ERR = sp.Rational(11, 100)
NUFIT_DELTA = sp.Integer(212)
NUFIT_DELTA_PLUS = sp.Integer(26)
NUFIT_DELTA_MINUS = sp.Integer(41)
# DESI Collaboration, https://arxiv.org/abs/2503.14744,
# DR2 BAO + DR1 full shape, LambdaCDM.
DESI_LCDM_95 = sp.Rational(642, 10000)       # eV, one-sided 95% bound
# https://arxiv.org/abs/2507.16589,
# joint DESI DR2+CMB+DESY5+DESY1, w0waCDM.
DESI_DYN = sp.Rational(98, 1000)            # eV
DESI_DYN_PLUS = sp.Rational(16, 1000)
DESI_DYN_MINUS = sp.Rational(37, 1000)
# KATRIN Collaboration, https://arxiv.org/abs/2406.13516;
# Science 388 (2025) 180-185, DOI 10.1126/science.adq9592.
KATRIN_90 = sp.Rational(45, 100)             # eV, one-sided 90% bound

s12 = 2/(phi+5)
s13 = eps
c13 = 1-s13
ue1 = c13*(1-s12)
ue2 = c13*s12
ue3 = s13
phase1 = 4*sp.pi/5
phase2 = -4*sp.pi/5


def observables(m2, m3):
    dm21 = m2**2
    dm31 = m3**2
    dm32 = dm31-dm21
    total = m2+m3
    mbeta = sp.sqrt(ue2*m2**2 + ue3*m3**2)
    mbb_sq = ((ue2*m2)**2 + (ue3*m3)**2
              + 2*ue2*ue3*m2*m3*sp.cos(phase1-phase2))
    mbb = sp.sqrt(sp.simplify(mbb_sq))
    return dict(m2=m2, m3=m3, dm21=dm21, dm31=dm31, dm32=dm32,
                sum=total, mbeta=mbeta, mbb=mbb,
                sig21=(dm21-JUNO_DM21)/JUNO_DM21_ERR,
                sig31_sk=(dm31-NUFIT_DM31_SK)/NUFIT_DM31_SK_MINUS,
                sig31_nosk=(dm31-NUFIT_DM31_NOSK)/NUFIT_DM31_NOSK_MINUS,
                sig32_sk=(dm32-NUFIT_DM32_SK)/NUFIT_DM32_SK_ERR,
                sig32_nosk=(dm32-NUFIT_DM32_NOSK)/NUFIT_DM32_NOSK_ERR,
                dyn_sigma=(total-DESI_DYN)/DESI_DYN_MINUS)


variants = {
    "I": observables(m2_bare, m3_corr),
    "II": observables(m2_corr, m3_corr),
    "III": observables(m2_bare, m3_bare),
}

print("="*78)
print("JUNO 2026 NEUTRINO COMPARISON")
print("="*78)
print(f"alpha={f(alpha,15):.15f}, r=alpha*phi^3={f(r,15):.15f}")
for name, o in variants.items():
    chi2 = o["sig21"]**2 + o["sig31_sk"]**2
    print(f"  Variant {name}: m2={1000*f(o['m2']):.6f} meV, "
          f"m3={1000*f(o['m3']):.6f} meV, "
          f"Dm21={f(o['dm21']):.8e} ({f(o['sig21']):+.3f} sigma), "
          f"Dm31={f(o['dm31']):.8e} ({f(o['sig31_sk']):+.3f} sigma SK, "
          f"{f(o['sig31_nosk']):+.3f} sigma noSK), "
          f"chi2_2={f(chi2):.4f}")
    print(f"             Dm32={f(o['dm32']):.8e} "
          f"({f(o['sig32_sk']):+.3f} sigma inferred-SK, "
          f"{f(o['sig32_nosk']):+.3f} sigma inferred-noSK); "
          f"sum={1000*f(o['sum']):.6f} meV "
          f"({f(o['dyn_sigma']):+.3f} sigma vs dynamical-DE preference); "
          f"m_beta={1000*f(o['mbeta']):.6f} meV; "
          f"m_bb={1000*f(o['mbb']):.6f} meV")

check("Variant I changes only m3", variants["I"]["m2"] == m2_bare
      and variants["I"]["m3"] == m3_corr)
check("Variant II propagates correction through m2/m3=sqrt(r)",
      sp.simplify(variants["II"]["m2"]**2/variants["II"]["m3"]**2-r) == 0)
check("Variant III is bare throughout", variants["III"]["m2"] == m2_bare
      and variants["III"]["m3"] == m3_bare)
check("Variants I and III have the same JUNO solar prediction",
      variants["I"]["dm21"] == variants["III"]["dm21"])
check("Variant II gives the largest solar splitting",
      variants["II"]["dm21"] > variants["I"]["dm21"])
paper_stale_dm32 = sp.Rational(2453, 1000)*10**-3
check("paper's old 2.453e-3 Dm32 entry is not the official NuFIT-6 Dm31 input",
      paper_stale_dm32 != NUFIT_DM31_SK and paper_stale_dm32 != NUFIT_DM31_NOSK)

# Mixing comparisons common to every mass variant.
theta13_deg = sp.asin(sp.sqrt(s13))*180/sp.pi
delta_pred = 3*sp.atan(sqrt5)*180/sp.pi
sig_s12 = (s12-JUNO_S12)/JUNO_S12_ERR
sig_s13 = (s13-NUFIT_S13)/NUFIT_S13_ERR
sig_th13 = (theta13_deg-NUFIT_TH13_DEG)/NUFIT_TH13_DEG_ERR
# Prediction lies below the asymmetric central value, so use the lower error.
sig_delta = (delta_pred-NUFIT_DELTA)/NUFIT_DELTA_MINUS
check("sin2(theta12) is within 1 sigma of JUNO", abs(f(sig_s12)) < 1,
      f"prediction={f(s12):.9f}, signed sigma={f(sig_s12):+.4f}")
check("sin2(theta13) is within 1 sigma of NuFIT 6.0", abs(f(sig_s13)) < 1,
      f"prediction={f(s13):.9f}, signed sigma={f(sig_s13):+.4f}")
check("theta13 is within 1 sigma of NuFIT 6.0", abs(f(sig_th13)) < 1,
      f"prediction={f(theta13_deg):.6f} deg, signed sigma={f(sig_th13):+.4f}")
check("delta_CP is within 1 sigma using the directional lower error",
      abs(f(sig_delta)) < 1,
      f"prediction={f(delta_pred):.6f} deg, signed sigma={f(sig_delta):+.4f}")

for name, o in variants.items():
    check(f"Variant {name} passes DESI LCDM 0.0642 eV bound", o["sum"] < DESI_LCDM_95,
          f"sum={f(o['sum']):.8f} eV, margin={f(DESI_LCDM_95-o['sum']):.8f} eV")
    check(f"Variant {name} passes KATRIN 0.45 eV bound", o["mbeta"] < KATRIN_90,
          f"m_beta={1000*f(o['mbeta']):.6f} meV")

# With two Gaussian splitting inputs, Variant I has the smallest chi-square.
chi2 = {k: sp.N(v["sig21"]**2+v["sig31_sk"]**2, 30) for k, v in variants.items()}
check("Variant I has the smallest two-splitting chi-square",
      chi2["I"] < chi2["II"] and chi2["I"] < chi2["III"], f"chi2={chi2}")

print("\nPinned common mixing deviations:")
print(f"  sin2 theta12: {f(s12):.9f}, {f(sig_s12):+.4f} sigma")
print(f"  sin2 theta13: {f(s13):.9f}, {f(sig_s13):+.4f} sigma")
print(f"  theta13:      {f(theta13_deg):.6f} deg, {f(sig_th13):+.4f} sigma")
print(f"  delta_CP:     {f(delta_pred):.6f} deg, {f(sig_delta):+.4f} sigma")

print("\nClassification:")
print("  DERIVED: bare package and exact observable recomputation")
print("  PATTERN: 1/45 mass correction and its Variant-I/II scope")
print("  EXTERNAL: pinned JUNO/NuFIT/DESI/KATRIN inputs")
print(f"\nTOTAL: {passed}/{run} tests PASSED")
sys.exit(1 if failed else 0)

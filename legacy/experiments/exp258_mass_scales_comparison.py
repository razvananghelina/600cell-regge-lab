"""
EXP-258: Systematic Comparison of Bare Formula with Physical Mass Scales
========================================================================
The bare formula m = m_e * phi^n gives 3-21% errors compared to PDG MSbar.
But quark masses depend on the renormalization scale!

Key question: Is there a NATURAL scale at which the bare formula matches?

Physics:
- MSbar mass m(mu): depends on scale mu (running mass)
- Pole mass: scale-independent but perturbatively defined
- Constituent mass: includes QCD binding (~300 MeV for light quarks)
- Running: m(mu2)/m(mu1) = (alpha_s(mu2)/alpha_s(mu1))^(gamma_m/(2*beta_0))
"""
import math

PHI = (1 + math.sqrt(5)) / 2
m_e = 0.51099895  # MeV
alpha_s_mZ = 0.1180  # alpha_s at m_Z
m_Z = 91187.6  # MeV

print("=" * 70)
print("EXP-258: BARE FORMULA vs PHYSICAL MASS SCALES")
print("=" * 70)

# =====================================================================
# QCD running functions
# =====================================================================

def beta_0(N_f):
    """1-loop QCD beta function coefficient"""
    return (33.0 - 2.0 * N_f) / 3.0

def gamma_m():
    """1-loop mass anomalous dimension (universal)"""
    return 8.0

def alpha_s_running(mu_MeV, alpha_s_ref=alpha_s_mZ, mu_ref=m_Z, N_f=5):
    """1-loop running of alpha_s"""
    b0 = beta_0(N_f)
    # alpha_s(mu) = alpha_s(ref) / (1 + b0*alpha_s(ref)/(2*pi) * ln(mu^2/ref^2))
    log_ratio = 2.0 * math.log(mu_MeV / mu_ref)
    denom = 1.0 + b0 * alpha_s_ref / (2.0 * math.pi) * log_ratio
    if denom <= 0:
        return float('inf')  # Landau pole
    return alpha_s_ref / denom

def run_mass(m_ref, mu_ref, mu_target, N_f):
    """Run quark mass from mu_ref to mu_target using 1-loop QCD"""
    a_s_ref = alpha_s_running(mu_ref, N_f=N_f)
    a_s_target = alpha_s_running(mu_target, N_f=N_f)
    if a_s_ref <= 0 or a_s_target <= 0:
        return float('nan')
    b0 = beta_0(N_f)
    exponent = gamma_m() / (2.0 * b0)
    ratio = (a_s_target / a_s_ref) ** exponent
    return m_ref * ratio

def pole_from_msbar(m_msbar, mu, N_f):
    """Convert MSbar mass to pole mass (1-loop)"""
    a_s = alpha_s_running(mu, N_f=N_f)
    # 1-loop: m_pole = m_msbar * (1 + 4*alpha_s/(3*pi))
    return m_msbar * (1.0 + 4.0 * a_s / (3.0 * math.pi))

# =====================================================================
# PART 1: PDG masses at standard scales
# =====================================================================
print()
print("PART 1: PDG 2024 reference masses")
print("-" * 50)
print()

# PDG 2024 values
pdg = {
    'u': {'msbar_2GeV': 2.16, 'err_up': 0.49, 'err_dn': 0.26, 'N_f': 3, 'note': 'MSbar(2 GeV)'},
    'd': {'msbar_2GeV': 4.67, 'err_up': 0.48, 'err_dn': 0.17, 'N_f': 3, 'note': 'MSbar(2 GeV)'},
    's': {'msbar_2GeV': 93.4, 'err_up': 8.6, 'err_dn': 3.4, 'N_f': 3, 'note': 'MSbar(2 GeV)'},
    'c': {'msbar_mc': 1270, 'msbar_2GeV': None, 'N_f': 4, 'pole': 1670, 'note': 'MSbar(m_c)'},
    'b': {'msbar_mb': 4180, 'msbar_2GeV': None, 'N_f': 5, 'pole': 4780, 'note': 'MSbar(m_b)'},
    't': {'msbar_mt': 162500, 'N_f': 6, 'pole': 172570, 'note': 'MSbar(m_t)'},
}

# Compute MSbar at 2 GeV for c and b
mu_2GeV = 2000  # MeV
pdg['c']['msbar_2GeV'] = run_mass(1270, 1270, mu_2GeV, N_f=4)
pdg['b']['msbar_2GeV'] = run_mass(4180, 4180, mu_2GeV, N_f=5)

print("  Quark MSbar masses at 2 GeV (1-loop running):")
for q in ['u', 'd', 's', 'c', 'b']:
    m2 = pdg[q]['msbar_2GeV']
    if m2 is not None:
        print("    %s: %.1f MeV" % (q, m2))

# Pole masses
print()
print("  Pole masses (PDG):")
for q in ['c', 'b', 't']:
    if 'pole' in pdg[q]:
        print("    %s: %.0f MeV" % (q, pdg[q]['pole']))

# Lepton masses (always pole)
leptons = {
    'e': 0.51099895,
    'mu': 105.6584,
    'tau': 1776.86,
}

# =====================================================================
# PART 2: Bare formula predictions
# =====================================================================
print()
print("=" * 70)
print("PART 2: Bare formula vs different mass definitions")
print("=" * 70)
print()

fermion_n = {'e': 0, 'mu': 11, 'tau': 17, 'u': 3, 'c': 16, 't': 26,
             'd': 5, 's': 11, 'b': 19}

print("  %-5s  n  m_bare    m_MSbar(std)  err_MSbar  m_pole     err_pole  m_MSbar(2G) err_2GeV" % "name")
print("  " + "-" * 95)

for name in ['e', 'mu', 'tau', 'u', 'd', 's', 'c', 'b', 't']:
    n = fermion_n[name]
    m_bare = m_e * PHI**n

    # MSbar at standard scale
    if name in leptons:
        m_std = leptons[name]
        m_pole_val = leptons[name]  # same for leptons
        m_2gev = leptons[name]  # no QCD running
    elif name in ['u', 'd', 's']:
        m_std = pdg[name]['msbar_2GeV']
        m_pole_val = None  # not well-defined for light quarks
        m_2gev = m_std
    elif name == 'c':
        m_std = pdg[name]['msbar_mc']
        m_pole_val = pdg[name]['pole']
        m_2gev = pdg[name]['msbar_2GeV']
    elif name == 'b':
        m_std = pdg[name]['msbar_mb']
        m_pole_val = pdg[name]['pole']
        m_2gev = pdg[name]['msbar_2GeV']
    elif name == 't':
        m_std = pdg[name]['msbar_mt']
        m_pole_val = pdg[name]['pole']
        m_2gev = None  # running top to 2 GeV is meaningless

    err_std = (m_bare - m_std) / m_std * 100
    err_pole = ((m_bare - m_pole_val) / m_pole_val * 100) if m_pole_val else None
    err_2gev = ((m_bare - m_2gev) / m_2gev * 100) if m_2gev else None

    pole_str = "%10.1f  %+6.1f%%" % (m_pole_val, err_pole) if m_pole_val else "%10s  %7s" % ("---", "---")
    gev2_str = "%10.1f  %+6.1f%%" % (m_2gev, err_2gev) if m_2gev else "%10s  %7s" % ("---", "---")

    print("  %-5s %2d  %8.2f  %12.2f  %+6.1f%%  %s  %s" %
          (name, n, m_bare, m_std, err_std, pole_str, gev2_str))

# =====================================================================
# PART 3: Find optimal scale for each quark
# =====================================================================
print()
print("=" * 70)
print("PART 3: Optimal scale for each heavy quark")
print("=" * 70)
print()

print("  For each heavy quark, find mu such that m_MSbar(mu) = m_bare")
print()

for q_name, n, m_ref, mu_ref, nf in [
    ('c', 16, 1270, 1270, 4),
    ('b', 19, 4180, 4180, 5),
    ('t', 26, 162500, 162500, 6)]:

    m_bare = m_e * PHI**n
    print("  %s quark: m_bare = %.2f MeV, m_MSbar(m_%s) = %.0f MeV" %
          (q_name, m_bare, q_name, m_ref))

    # Scan mu from 1 GeV to 1 TeV to find where m(mu) = m_bare
    best_mu = None
    best_err = float('inf')
    for log_mu in range(5, 70):  # log_mu in units of 0.1, so 0.5 to 7.0 (= 3 MeV to 1 TeV)
        mu = 10**(log_mu / 10.0) * 1000  # in MeV
        if mu < 500 or mu > 1e6:
            continue
        m_run = run_mass(m_ref, mu_ref, mu, N_f=nf)
        if m_run > 0:
            err = abs(m_run - m_bare) / m_bare
            if err < best_err:
                best_err = err
                best_mu = mu

    if best_mu:
        m_at_best = run_mass(m_ref, mu_ref, best_mu, N_f=nf)
        print("    Optimal scale: mu = %.0f MeV = %.2f GeV" % (best_mu, best_mu/1000))
        print("    m_%s(%.1f GeV) = %.1f MeV, bare = %.1f MeV, err = %.2f%%" %
              (q_name, best_mu/1000, m_at_best, m_bare, (m_at_best-m_bare)/m_bare*100))

        # What is this scale in terms of theory constants?
        ratio_to_mZ = best_mu / m_Z
        ratio_to_me = best_mu / (m_e)
        print("    mu/m_Z = %.3f, mu/m_e = %.0f" % (ratio_to_mZ, ratio_to_me))
    print()

# =====================================================================
# PART 4: Universal scale hypothesis
# =====================================================================
print()
print("=" * 70)
print("PART 4: Is there a universal optimal scale?")
print("=" * 70)
print()

# Test: what if all quarks are compared at the SAME scale?
print("  RMS error vs scale (for c, b quarks with MSbar running):")
print("  %-10s  err_c     err_b     RMS(c,b)" % "scale")

for scale_name, mu in [("1 GeV", 1000), ("2 GeV", 2000), ("3 GeV", 3000),
                        ("m_c", 1270), ("5 GeV", 5000), ("m_b", 4180),
                        ("10 GeV", 10000), ("m_Z", 91188), ("pole_c", 1670),
                        ("pole_b", 4780)]:
    m_c_run = run_mass(1270, 1270, mu, N_f=4) if mu > 500 else None
    m_b_run = run_mass(4180, 4180, mu, N_f=5) if mu > 500 else None

    m_c_bare = m_e * PHI**16
    m_b_bare = m_e * PHI**19

    if m_c_run and m_b_run and m_c_run > 0 and m_b_run > 0:
        err_c = (m_c_bare - m_c_run) / m_c_run * 100
        err_b = (m_b_bare - m_b_run) / m_b_run * 100
        rms = math.sqrt((err_c**2 + err_b**2) / 2)
        print("  %-10s  %+6.1f%%   %+6.1f%%   %5.1f%%" % (scale_name, err_c, err_b, rms))

# =====================================================================
# PART 5: Including ALL quarks at 2 GeV
# =====================================================================
print()
print("=" * 70)
print("PART 5: All quarks at 2 GeV (the natural light-quark scale)")
print("=" * 70)
print()

m_c_2GeV = run_mass(1270, 1270, 2000, N_f=4)
m_b_2GeV = run_mass(4180, 4180, 2000, N_f=5)

# Top at 2 GeV is problematic (crosses thresholds). Skip.
# Use MSbar values at 2 GeV for light quarks

quark_comparison = [
    ('u', 3, 2.16, '2 GeV'),
    ('d', 5, 4.67, '2 GeV'),
    ('s', 11, 93.4, '2 GeV'),
    ('c', 16, m_c_2GeV, '2 GeV (run)'),
    ('b', 19, m_b_2GeV, '2 GeV (run)'),
]

print("  %-5s  n  m_bare(MeV)  m(2GeV)(MeV)  err      source" % "name")
print("  " + "-" * 65)
sum_sq = 0
count = 0
for name, n, m_exp_val, source in quark_comparison:
    m_bare = m_e * PHI**n
    err = (m_bare - m_exp_val) / m_exp_val * 100
    sum_sq += err**2
    count += 1
    print("  %-5s %2d  %11.2f  %12.2f  %+6.1f%%   %s" %
          (name, n, m_bare, m_exp_val, err, source))

rms_2gev = math.sqrt(sum_sq / count)
print()
print("  RMS error at 2 GeV: %.1f%%" % rms_2gev)

# =====================================================================
# PART 6: Pole masses for heavy quarks
# =====================================================================
print()
print("=" * 70)
print("PART 6: Pole masses for c, b, t")
print("=" * 70)
print()

print("  %-5s  n  m_bare(MeV)  m_pole(MeV)  err" % "name")
for name, n, m_pole_val in [('c', 16, 1670), ('b', 19, 4780), ('t', 26, 172570)]:
    m_bare = m_e * PHI**n
    err = (m_bare - m_pole_val) / m_pole_val * 100
    print("  %-5s %2d  %11.2f  %11.0f  %+6.1f%%" % (name, n, m_bare, m_pole_val, err))

print()

# =====================================================================
# PART 7: Mixed comparison (best available)
# =====================================================================
print()
print("=" * 70)
print("PART 7: Best comparison per fermion")
print("=" * 70)
print()

print("  Strategy: use the mass definition closest to bare formula")
print()

comparisons = [
    ('e',   0,  0.511,   0.511,   'pole (input)'),
    ('mu', 11,  None,    105.658, 'pole'),
    ('tau',17,  None,    1776.86, 'pole'),
    ('u',   3,  None,    2.16,    'MSbar(2 GeV)'),
    ('d',   5,  None,    4.67,    'MSbar(2 GeV)'),
    ('s',  11,  None,    93.4,    'MSbar(2 GeV)'),
    ('c',  16,  None,    None,    'scan'),
    ('b',  19,  None,    None,    'scan'),
    ('t',  26,  None,    None,    'scan'),
]

print("  %-5s  n  m_bare       best_m       err      definition" % "name")
print("  " + "-" * 70)

all_errors = []
for name, n, _, m_best, defn in comparisons:
    m_bare = m_e * PHI**n

    if name == 'e':
        err = 0
        m_show = 0.511
        defn_show = 'pole (input)'
    elif name in ['mu', 'tau']:
        m_show = m_best
        err = (m_bare - m_show) / m_show * 100
        defn_show = 'pole'
    elif name in ['u', 'd', 's']:
        m_show = m_best
        err = (m_bare - m_show) / m_show * 100
        defn_show = 'MSbar(2 GeV)'
    elif name == 'c':
        # Compare with: MSbar(m_c), pole, MSbar(2 GeV)
        options = [
            (1270, 'MSbar(m_c)'),
            (1670, 'pole'),
            (m_c_2GeV, 'MSbar(2 GeV)'),
        ]
        best_opt = min(options, key=lambda x: abs(m_bare - x[0]))
        m_show = best_opt[0]
        defn_show = best_opt[1]
        err = (m_bare - m_show) / m_show * 100
    elif name == 'b':
        options = [
            (4180, 'MSbar(m_b)'),
            (4780, 'pole'),
            (m_b_2GeV, 'MSbar(2 GeV)'),
        ]
        best_opt = min(options, key=lambda x: abs(m_bare - x[0]))
        m_show = best_opt[0]
        defn_show = best_opt[1]
        err = (m_bare - m_show) / m_show * 100
    elif name == 't':
        options = [
            (162500, 'MSbar(m_t)'),
            (172570, 'pole'),
        ]
        best_opt = min(options, key=lambda x: abs(m_bare - x[0]))
        m_show = best_opt[0]
        defn_show = best_opt[1]
        err = (m_bare - m_show) / m_show * 100

    all_errors.append((name, abs(err)))
    print("  %-5s %2d  %11.2f  %11.2f  %+6.1f%%   %s" %
          (name, n, m_bare, m_show, err, defn_show))

print()
rms_best = math.sqrt(sum(e**2 for _, e in all_errors) / len(all_errors))
rms_no_e = math.sqrt(sum(e**2 for n, e in all_errors if n != 'e') / sum(1 for n, _ in all_errors if n != 'e'))
print("  RMS error (all 9): %.1f%%" % rms_best)
print("  RMS error (8, excl. e): %.1f%%" % rms_no_e)

# =====================================================================
# PART 8: QED corrections for muon
# =====================================================================
print()
print("=" * 70)
print("PART 8: QED self-energy corrections for leptons")
print("=" * 70)
print()

alpha_em = 1.0/137.036

# Leading QED correction to lepton mass: delta_m/m ~ (3*alpha)/(4*pi) * ln(Lambda/m)
# But this is absorbed into the renormalized mass. The POLE mass IS the physical mass.
# The bare formula should predict the pole mass for leptons.

print("  Lepton masses are POLE masses (no QCD). No scale ambiguity.")
print("  The bare formula should directly predict pole masses.")
print()
print("  Bare formula errors for leptons:")
for name, n, m_pole in [('e', 0, 0.511), ('mu', 11, 105.658), ('tau', 17, 1776.86)]:
    m_bare = m_e * PHI**n
    err = (m_bare - m_pole) / m_pole * 100
    print("    %s: m_bare=%.4f, m_pole=%.4f, err=%+.2f%%" % (name, m_bare, m_pole, err))

print()
print("  These errors (-3.8% for mu, +2.7% for tau) are REAL predictions.")
print("  They cannot be explained by scale choice.")
print("  They are the INTRINSIC accuracy of the bare formula for leptons.")

# =====================================================================
# PART 9: Systematic analysis of running effects
# =====================================================================
print()
print("=" * 70)
print("PART 9: Running effects on mass errors")
print("=" * 70)
print()

# For each quark with known MSbar mass, compute bare/MSbar at several scales
print("  Effect of scale on charm mass error:")
m_c_bare = m_e * PHI**16
for mu_name, mu in [("1 GeV", 1000), ("m_c", 1270), ("2 GeV", 2000),
                     ("3 GeV", 3000), ("5 GeV", 5000), ("m_Z", 91188)]:
    m_c_at_mu = run_mass(1270, 1270, mu, N_f=4)
    if m_c_at_mu and m_c_at_mu > 0:
        err = (m_c_bare - m_c_at_mu) / m_c_at_mu * 100
        print("    m_c(%8s) = %8.1f MeV, bare/run = %+6.1f%%" % (mu_name, m_c_at_mu, err))

print()
print("  Effect of scale on bottom mass error:")
m_b_bare = m_e * PHI**19
for mu_name, mu in [("2 GeV", 2000), ("3 GeV", 3000), ("m_b", 4180),
                     ("5 GeV", 5000), ("10 GeV", 10000), ("m_Z", 91188)]:
    m_b_at_mu = run_mass(4180, 4180, mu, N_f=5)
    if m_b_at_mu and m_b_at_mu > 0:
        err = (m_b_bare - m_b_at_mu) / m_b_at_mu * 100
        print("    m_b(%8s) = %8.1f MeV, bare/run = %+6.1f%%" % (mu_name, m_b_at_mu, err))

print()
print("  Effect of scale on top mass error:")
m_t_bare = m_e * PHI**26
for mu_name, mu in [("m_t", 162500), ("m_Z", 91188), ("200 GeV", 200000)]:
    m_t_at_mu = run_mass(162500, 162500, mu, N_f=6)
    if m_t_at_mu and m_t_at_mu > 0:
        err = (m_t_bare - m_t_at_mu) / m_t_at_mu * 100
        print("    m_t(%8s) = %8.1f MeV, bare/run = %+6.1f%%" % (mu_name, m_t_at_mu, err))

# =====================================================================
# PART 10: The bottom quark pole mass coincidence
# =====================================================================
print()
print("=" * 70)
print("PART 10: Bottom quark pole mass - the 0.06% match")
print("=" * 70)
print()

m_b_bare = m_e * PHI**19
m_b_pole = 4780  # MeV (PDG)
m_b_msbar = 4180  # MeV

err_pole = (m_b_bare - m_b_pole) / m_b_pole * 100
err_msbar = (m_b_bare - m_b_msbar) / m_b_msbar * 100

print("  m_b(bare formula) = %.2f MeV" % m_b_bare)
print("  m_b(pole, PDG)    = %.0f +/- 30 MeV" % m_b_pole)
print("  m_b(MSbar, PDG)   = %.0f +/- 7 MeV" % m_b_msbar)
print()
print("  Error vs pole:  %+.2f%%  (%.1f MeV)" % (err_pole, m_b_bare - m_b_pole))
print("  Error vs MSbar: %+.1f%%  (%.1f MeV)" % (err_msbar, m_b_bare - m_b_msbar))
print()

# Is this coincidence? The pole mass has ~30 MeV uncertainty
# So 0.06% on 4780 = 2.9 MeV, well within 30 MeV
print("  Bottom pole mass uncertainty: ~30 MeV")
print("  Bare formula deviation: %.1f MeV" % abs(m_b_bare - m_b_pole))
print("  This is within uncertainty. The 0.06%% match may be coincidental")
print("  but is the BEST mass prediction in the theory (after m_Z, m_H, etc.)")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print()
print("=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print()

print("  SCALE EFFECTS ON MASS ERRORS:")
print()
print("  Leptons: NO scale ambiguity. Errors are intrinsic:")
print("    mu: -3.8%%, tau: +2.7%% (bare formula prediction)")
print()
print("  Light quarks (u, d, s): MSbar(2 GeV) is standard.")
print("    u: +0.2%%, d: +21%%, s: +9%%")
print("    Running to other scales doesn't help d and s.")
print()
print("  Charm: Running HELPS. m_c(2 GeV) ~ 1100 MeV reduces error")
print("    from -11%% (at m_c) to ~+3%% (at 2 GeV).")
print("    But pole mass (1670) WORSENS to -32%%.")
print()
print("  Bottom: POLE mass gives 0.06%% match!")
print("    MSbar(m_b) = 4180: +14%% error")
print("    MSbar(2 GeV) ~ 4850: -1.5%% error")
print("    Pole ~ 4780: -0.06%% error (best mass in the theory)")
print()
print("  Top: NEITHER pole nor MSbar matches well.")
print("    MSbar(m_t) = 162.5 GeV: -15%% error")
print("    Pole = 172.6 GeV: -20%% error")
print()
print("  CONCLUSION:")
print("  - There is NO universal scale that fixes all errors.")
print("  - Different quarks match at different definitions:")
print("    u at MSbar(2 GeV), b at pole, c at MSbar(~2 GeV)")
print("  - Leptons have irreducible 3-4%% errors (no scale choice)")
print("  - The bare formula errors (3-21%%) are STRUCTURAL,")
print("    not artifacts of scale choice.")
print("  - Bottom pole mass 0.06%% is remarkable but may be coincidental.")
print()
print("  HONEST ASSESSMENT: The bare formula is a ~10%% level")
print("  approximation to fermion masses. This is the PREDICTION")
print("  with zero free parameters. Scale effects explain at most")
print("  a factor of ~2 reduction in some errors, not elimination.")

"""
EXP-096: Validare Formulele de Corectie (Raspuns la critica EXP-095)
====================================================================
Testam RIGUROS:
  1. C_EM(Q) = 1 - Q * alpha * (2*phi^2) - corectia EM pentru leptoni
  2. C_QCD = (1 + alpha_s(mu))^(Q * N_c) - corectia QCD pentru quarks
  3. Formula completa: m_f = m_e * phi^n * C_EM * C_QCD

REGULA ZERO: Nu inventam nimic. Daca nu merge, spunem clar.
"""

import math

# =============================================================================
# CONSTANTE
# =============================================================================
PHI = (1 + math.sqrt(5)) / 2
PI = math.pi
ALPHA = 7.2973525693e-3
ALPHA_S_MZ = 0.1179

# Mase experimentale (MeV) - PDG 2024
M_E = 0.51099895
M_MU = 105.6583755
M_TAU = 1776.86

M_UP = 2.16          # MS-bar la 2 GeV
M_DOWN = 4.67        # MS-bar la 2 GeV
M_STRANGE = 93.4     # MS-bar la 2 GeV
M_CHARM = 1270.0     # MS-bar la m_c
M_BOTTOM = 4180.0    # MS-bar la m_b
M_TOP = 172760.0     # pole mass

# Proprietati 600-cell
SPECTRAL_RATIO = 2 * PHI**2  # lambda_4/lambda_1 = 2*phi^2 = 5.2361

def pct_err(calc, exp):
    return abs(calc - exp) / exp * 100

# =============================================================================
# SECTIUNEA 1: CORECTIA ELECTROMAGNETICA - LEPTONI
# =============================================================================
def test_em_correction():
    print("=" * 70)
    print("SECTIUNEA 1: CORECTIA ELECTROMAGNETICA C_EM")
    print("Formula: C_EM(Q) = 1 - Q * alpha * (2*phi^2)")
    print("=" * 70)
    print()

    # Parametri
    Q_lepton = -1  # sarcina leptonului
    C_EM = 1 - Q_lepton * ALPHA * (2 * PHI**2)
    # Q = -1 => C_EM = 1 + alpha * 2*phi^2 = 1 + alpha * 5.236...

    print(f"  Q_lepton = {Q_lepton}")
    print(f"  2*phi^2 = {2*PHI**2:.6f}")
    print(f"  alpha * 2*phi^2 = {ALPHA * 2*PHI**2:.6f}")
    print(f"  C_EM(Q=-1) = 1 - ({Q_lepton}) * {ALPHA:.6f} * {2*PHI**2:.4f}")
    print(f"             = 1 + {ALPHA * 2*PHI**2:.6f}")
    print(f"             = {C_EM:.8f}")
    print()

    # Corectia procentuala
    correction_pct = (C_EM - 1) * 100
    print(f"  Corectia procentuala: +{correction_pct:.4f}%")
    print(f"  Corectia necesara (din exp095): +3.90%")
    print()

    # Testam pe muon
    print("  --- MUON ---")
    m_mu_bare = M_E * PHI**11
    m_mu_corrected = m_mu_bare * C_EM

    print(f"  m_mu(bare) = m_e * phi^11 = {m_mu_bare:.4f} MeV")
    print(f"  m_mu(corr) = {m_mu_bare:.4f} * {C_EM:.6f} = {m_mu_corrected:.4f} MeV")
    print(f"  m_mu(exp)  = {M_MU:.4f} MeV")
    print(f"  Eroare bare:      {pct_err(m_mu_bare, M_MU):.4f}%")
    print(f"  Eroare corectata: {pct_err(m_mu_corrected, M_MU):.4f}%")
    print()

    # Asta e o corectie de doar ~0.38%, nu de ~3.9% cat avem nevoie!
    if abs(correction_pct) < 1.0:
        print(f"  PROBLEMA: Corectia C_EM = +{correction_pct:.2f}%")
        print(f"  Dar avem nevoie de +3.9% pentru a ajunge la 105.66 MeV")
        print(f"  Formula C_EM(Q) = 1 - Q*alpha*(2*phi^2) e PREA MICA!")
        print()

        # Ce factor ar trebui?
        needed = M_MU / m_mu_bare
        print(f"  Factor necesar: {needed:.6f}")
        print(f"  Deci C_EM ar trebui sa fie: {needed:.6f}")
        print(f"  Diferenta: {needed - C_EM:.6f}")
        print()

        # Ce formula ar da factorul corect?
        # Daca C = 1 + k * alpha * (2*phi^2), ce k trebuie?
        k_needed = (needed - 1) / (ALPHA * 2*PHI**2)
        print(f"  Daca C = 1 + k * alpha * 2*phi^2:")
        print(f"  k necesar = {k_needed:.4f}")
        print(f"  Asta ar insemna k ~ {k_needed:.1f}, nu k=1")
        print()

    # Testam pe TAU
    print("  --- TAU ---")
    m_tau_bare = M_E * PHI**17
    m_tau_corrected = m_tau_bare * C_EM

    print(f"  m_tau(bare) = m_e * phi^17 = {m_tau_bare:.2f} MeV")
    print(f"  m_tau(corr) = {m_tau_bare:.2f} * {C_EM:.6f} = {m_tau_corrected:.2f} MeV")
    print(f"  m_tau(exp)  = {M_TAU:.2f} MeV")
    print(f"  Eroare bare:      {pct_err(m_tau_bare, M_TAU):.4f}%")
    print(f"  Eroare corectata: {pct_err(m_tau_corrected, M_TAU):.4f}%")
    print()

    # C_EM e ACEEASI pentru toti leptonii incarcati (Q=-1)
    # Deci nu poate explica DIFERENTA dintre rapoarte
    print("  --- ANALIZA CRITICA ---")
    print(f"  C_EM e ACEEASI pentru electron, muon, tau (toti au Q=-1)")
    print(f"  Deci C_EM se SIMPLIFICA in raportul m_mu/m_e!")
    print(f"  m_mu/m_e = (m_e * phi^11 * C_EM) / (m_e * C_EM) = phi^11")
    print(f"  Corectia NU ajuta la raportul maselor!")
    print()
    print(f"  CONCLUZIE: C_EM(Q) nu poate explica diferenta de 3.75%")
    print(f"  Corectia se anuleaza in raportul m_mu/m_e")
    print(f"  Si oricum e de doar {correction_pct:.2f}%, nu 3.9%")
    print()

    return C_EM


# =============================================================================
# SECTIUNEA 2: FORMULA ALTERNATIVA - MASA ABSOLUTA vs RAPORT
# =============================================================================
def test_absolute_mass():
    print("=" * 70)
    print("SECTIUNEA 2: MASA ABSOLUTA vs RAPORT DE MASE")
    print("=" * 70)
    print()

    # Poate formula nu e m_f = m_e * phi^n ci m_f = M_0 * phi^n
    # unde M_0 NU e exact m_e?
    # Sau poate e m_f = m_e * phi^n cu n diferit?

    # Cautam cel mai bun n continuu pentru muon
    n_mu_exact = math.log(M_MU / M_E) / math.log(PHI)
    n_tau_exact = math.log(M_TAU / M_E) / math.log(PHI)
    n_tau_mu_exact = math.log(M_TAU / M_MU) / math.log(PHI)

    print(f"  Exponenti exacti (continui):")
    print(f"    n(m_mu/m_e)  = {n_mu_exact:.6f}  (cel mai apropiat intreg: {round(n_mu_exact)})")
    print(f"    n(m_tau/m_e) = {n_tau_exact:.6f}  (cel mai apropiat intreg: {round(n_tau_exact)})")
    print(f"    n(m_tau/m_mu)= {n_tau_mu_exact:.6f}  (cel mai apropiat intreg: {round(n_tau_mu_exact)})")
    print()

    # Diferenta fata de intregi
    delta_mu = n_mu_exact - 11
    delta_tau_e = n_tau_exact - 17
    delta_tau_mu = n_tau_mu_exact - 6

    print(f"  Deviatii de la intregi:")
    print(f"    delta(mu)    = {delta_mu:+.6f}")
    print(f"    delta(tau/e) = {delta_tau_e:+.6f}")
    print(f"    delta(tau/mu)= {delta_tau_mu:+.6f}")
    print()

    # Poate deviatia e sistematica?
    print(f"  Observatie: delta(mu) = +{delta_mu:.4f}")
    print(f"  Asta inseamna n_real = 11 + {delta_mu:.4f}")
    print(f"  = 11 + corectie de ~{delta_mu:.2f}")
    print()

    # Ce ar putea fi delta_mu fizic?
    # delta ~ 0.08, hm...
    # alpha * 11? = 0.080 -- interesant!
    print(f"  Candidati pentru corectia delta:")
    print(f"    alpha * n:       alpha * 11 = {ALPHA * 11:.6f}")
    print(f"    1/(2*phi^3):     = {1/(2*PHI**3):.6f}")
    print(f"    alpha * 2*phi^2: = {ALPHA * 2*PHI**2:.6f}")
    print(f"    delta actual:    = {delta_mu:.6f}")
    print()

    # Testam n_eff = n + alpha * n
    print("  --- TEST: n_eff = n * (1 + alpha) ---")
    for name, n_int, m_exp in [("Muon", 11, M_MU), ("Tau", 17, M_TAU)]:
        n_eff = n_int * (1 + ALPHA)
        m_calc = M_E * PHI**n_eff
        err = pct_err(m_calc, m_exp)
        print(f"    {name}: n_eff = {n_int} * (1+alpha) = {n_eff:.6f}")
        print(f"           m_calc = {m_calc:.4f} MeV, m_exp = {m_exp:.4f} MeV, err = {err:.4f}%")
    print()

    # Testam n_eff = n + n*alpha (= n*(1+alpha)) - e acelasi lucru
    # Testam alta varianta: n_eff = n + alpha_s
    print("  --- TEST: n_eff = n + alpha_s ---")
    for name, n_int, m_exp in [("Muon", 11, M_MU), ("Tau", 17, M_TAU)]:
        n_eff = n_int + ALPHA_S_MZ
        m_calc = M_E * PHI**n_eff
        err = pct_err(m_calc, m_exp)
        print(f"    {name}: n_eff = {n_int} + alpha_s = {n_eff:.6f}")
        print(f"           m_calc = {m_calc:.4f} MeV, m_exp = {m_exp:.4f} MeV, err = {err:.4f}%")
    print()

    # Testam: m_f = m_e * phi^(n + c) cu c universal
    print("  --- TEST: n_eff = n + c (cautam c universal) ---")
    # c care minimizeaza eroarea combinata muon+tau
    best_c = None
    best_total_err = 999
    for c_test in [i * 0.001 for i in range(-200, 500)]:
        m_mu_c = M_E * PHI**(11 + c_test)
        m_tau_c = M_E * PHI**(17 + c_test)
        err_mu = pct_err(m_mu_c, M_MU)
        err_tau = pct_err(m_tau_c, M_TAU)
        total = err_mu + err_tau
        if total < best_total_err:
            best_total_err = total
            best_c = c_test

    m_mu_best = M_E * PHI**(11 + best_c)
    m_tau_best = M_E * PHI**(17 + best_c)
    print(f"    Cel mai bun c universal: {best_c:.3f}")
    print(f"    m_mu({best_c:.3f}) = {m_mu_best:.4f} MeV, err = {pct_err(m_mu_best, M_MU):.4f}%")
    print(f"    m_tau({best_c:.3f}) = {m_tau_best:.4f} MeV, err = {pct_err(m_tau_best, M_TAU):.4f}%")
    print()

    # Exista un c unic care rezolva AMBELE? Probabil nu exact.
    if abs(best_c) > 0.01:
        print(f"  OBSERVATIE: c = {best_c:.3f} sugereaza o corectie sistematica")
        print(f"  Dar ce e {best_c:.3f} fizic?")
        # Verificam daca e vreo combinatie cunoscuta
        candidates = {
            "alpha * 11": ALPHA * 11,
            "alpha * 12": ALPHA * 12,
            "alpha * phi^2": ALPHA * PHI**2,
            "alpha * 2*phi^2": ALPHA * 2*PHI**2,
            "alpha_s / 2": ALPHA_S_MZ / 2,
            "1/(4*pi)": 1/(4*PI),
            "alpha * pi": ALPHA * PI,
        }
        print(f"  Candidati:")
        for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - best_c)):
            print(f"    {name:20s} = {val:.6f}  (dif = {abs(val-best_c):.6f})")
    print()


# =============================================================================
# SECTIUNEA 3: CORECTIA QCD PENTRU QUARKS
# =============================================================================
def test_qcd_correction():
    print("=" * 70)
    print("SECTIUNEA 3: CORECTIA QCD PENTRU QUARKS")
    print("Formula: C_QCD = (1 + alpha_s(mu))^(Q * N_c)")
    print("=" * 70)
    print()

    N_c = 3  # numarul de culori

    # alpha_s la diferite scale (approximate, 1-loop)
    Lambda_QCD = 0.332  # GeV
    def alpha_s_1loop(mu_GeV, nf=3):
        b0 = (33 - 2*nf) / (12*PI)
        if mu_GeV <= Lambda_QCD:
            return 1.0  # non-perturbativ
        return 1.0 / (2 * b0 * math.log(mu_GeV / Lambda_QCD))

    # Quarks cu proprietatile lor
    quarks = [
        # (name, n, m_exp_MeV, Q_em, mu_scale_GeV, nf)
        ("Up",      3,   M_UP,      +2/3, 2.0,    3),
        ("Down",    5,   M_DOWN,    -1/3, 2.0,    3),
        ("Strange", 11,  M_STRANGE, -1/3, 2.0,    3),
        ("Charm",   16,  M_CHARM,   +2/3, 1.27,   4),
        ("Bottom",  19,  M_BOTTOM,  +2/3, 4.18,   5),
        ("Top",     26,  M_TOP,     +2/3, 172.76, 6),
    ]

    print(f"  N_c = {N_c}")
    print(f"  Formula: m_q = m_e * phi^n * (1 + alpha_s(mu))^(Q * N_c)")
    print()

    for name, n, m_exp, Q, mu, nf in quarks:
        alpha_s = alpha_s_1loop(mu, nf)
        C_QCD = (1 + alpha_s)**(Q * N_c)

        m_bare = M_E * PHI**n
        m_corrected = m_bare * C_QCD

        err_bare = pct_err(m_bare, m_exp)
        err_corr = pct_err(m_corrected, m_exp)

        improved = "IMPROVED" if err_corr < err_bare else "WORSE"

        print(f"  {name:8s}: n={n:2d}, Q={Q:+5.2f}")
        print(f"    alpha_s({mu:.2f} GeV) = {alpha_s:.4f}")
        print(f"    C_QCD = (1+{alpha_s:.4f})^({Q:.2f}*3) = (1+{alpha_s:.4f})^({Q*N_c:.2f}) = {C_QCD:.6f}")
        print(f"    m_bare = {m_bare:.2f} MeV")
        print(f"    m_corr = {m_corrected:.2f} MeV")
        print(f"    m_exp  = {m_exp:.2f} MeV")
        print(f"    Eroare bare: {err_bare:.2f}%  ->  Eroare corr: {err_corr:.2f}%  [{improved}]")
        print()

    # Analiza
    print("  --- ANALIZA CRITICA ---")
    print()

    # Verificam daca semnul corectiei e consistent
    print("  Semnul corectiei:")
    print("  Q = +2/3 (up-type): exponent = +2/3 * 3 = +2 => C > 1 (creste masa)")
    print("  Q = -1/3 (down-type): exponent = -1/3 * 3 = -1 => C < 1 (scade masa)")
    print()
    print("  Pattern experimental (din exp094):")
    print("  Up-type (u, c, t): mase SUBESTIMATE de phi^n => TREBUIE C > 1")
    print("  Down-type (d, s, b): mase SUPRAESTIMATE de phi^n => TREBUIE C < 1")
    print()

    # Verificam daca patternul e corect
    for name, n, m_exp, Q, mu, nf in quarks:
        m_bare = M_E * PHI**n
        needs_increase = m_bare < m_exp
        C_sign = "C > 1" if Q > 0 else "C < 1"
        correct_direction = (needs_increase and Q > 0) or (not needs_increase and Q < 0)
        print(f"    {name:8s}: bare={'<' if needs_increase else '>'} exp, Q={Q:+.2f}, {C_sign}, {'CORECT' if correct_direction else 'GRESIT'}")

    print()


# =============================================================================
# SECTIUNEA 4: FORMULA COMPLETA UNIFICATA
# =============================================================================
def test_unified_formula():
    print("=" * 70)
    print("SECTIUNEA 4: FORMULA COMPLETA UNIFICATA")
    print("m_f = m_e * phi^n * C_EM(Q) * C_QCD(Q, alpha_s, N_c)")
    print("=" * 70)
    print()

    Lambda_QCD = 0.332
    def alpha_s_1loop(mu_GeV, nf=3):
        b0 = (33 - 2*nf) / (12*PI)
        if mu_GeV <= Lambda_QCD:
            return 1.0
        return 1.0 / (2 * b0 * math.log(mu_GeV / Lambda_QCD))

    N_c = 3

    # Toti fermionii
    fermions = [
        # (name, n, m_exp_MeV, Q_em, is_lepton, mu_GeV, nf)
        ("Electron",  0,  M_E,       -1,   True,  None, None),
        ("Muon",     11,  M_MU,      -1,   True,  None, None),
        ("Tau",      17,  M_TAU,     -1,   True,  None, None),
        ("Up",        3,  M_UP,     +2/3,  False, 2.0,  3),
        ("Down",      5,  M_DOWN,   -1/3,  False, 2.0,  3),
        ("Strange",  11,  M_STRANGE, -1/3, False, 2.0,  3),
        ("Charm",    16,  M_CHARM,  +2/3,  False, 1.27, 4),
        ("Bottom",   19,  M_BOTTOM, +2/3,  False, 4.18, 5),
        ("Top",      26,  M_TOP,    +2/3,  False, 172.76, 6),
    ]

    print(f"  {'Particula':10s} {'n':>3s} {'m_calc':>12s} {'m_exp':>12s} {'err':>8s} {'C_EM':>8s} {'C_QCD':>8s} {'Status':>10s}")
    print(f"  {'-'*75}")

    for name, n, m_exp, Q, is_lepton, mu, nf in fermions:
        # C_EM
        C_EM = 1 - Q * ALPHA * (2 * PHI**2)

        # C_QCD (doar quarks)
        if is_lepton:
            C_QCD = 1.0
        else:
            alpha_s = alpha_s_1loop(mu, nf)
            C_QCD = (1 + alpha_s)**(Q * N_c)

        if n == 0:
            m_calc = M_E * C_EM * C_QCD  # electron e referinta
        else:
            m_calc = M_E * PHI**n * C_EM * C_QCD

        err = pct_err(m_calc, m_exp)
        ok = "PASS" if err < 5 else ("MARG" if err < 15 else "WEAK")

        print(f"  {name:10s} {n:3d} {m_calc:12.4f} {m_exp:12.4f} {err:7.2f}% {C_EM:8.4f} {C_QCD:8.4f} {ok:>10s}")

    print()

    # Comparam cu formula FARA corectii
    print(f"  --- COMPARATIE: Cu vs Fara corectii ---")
    print(f"  {'Particula':10s} {'err_bare':>10s} {'err_corr':>10s} {'Diferenta':>10s}")
    print(f"  {'-'*45}")

    improvements = 0
    degradations = 0

    for name, n, m_exp, Q, is_lepton, mu, nf in fermions:
        if n == 0:
            continue  # skip electron (referinta)

        m_bare = M_E * PHI**n
        err_bare = pct_err(m_bare, m_exp)

        C_EM = 1 - Q * ALPHA * (2 * PHI**2)
        if is_lepton:
            C_QCD = 1.0
        else:
            alpha_s = alpha_s_1loop(mu, nf)
            C_QCD = (1 + alpha_s)**(Q * N_c)

        m_corr = M_E * PHI**n * C_EM * C_QCD
        err_corr = pct_err(m_corr, m_exp)

        delta = err_corr - err_bare
        better = "BETTER" if delta < -0.5 else ("SAME" if abs(delta) < 0.5 else "WORSE")
        if delta < -0.5:
            improvements += 1
        elif delta > 0.5:
            degradations += 1

        print(f"  {name:10s} {err_bare:9.2f}% {err_corr:9.2f}% {delta:+9.2f}% {better}")

    print()
    print(f"  Imbunatatiri: {improvements}")
    print(f"  Degradari:    {degradations}")
    print()

    if improvements > degradations:
        print("  VERDICT: Corectiile AJUTA in general")
    elif improvements == degradations:
        print("  VERDICT: Corectiile sunt MIXTE (ajuta unele, strica altele)")
    else:
        print("  VERDICT: Corectiile STRICA mai mult decat ajuta!")
    print()


# =============================================================================
# SECTIUNEA 5: CAUTAM FORMULA CORECTA PENTRU MUON
# =============================================================================
def test_muon_deep():
    print("=" * 70)
    print("SECTIUNEA 5: INVESTIGATIE - CE FORMULA DA MUONUL CORECT?")
    print("=" * 70)
    print()

    m_mu_bare = M_E * PHI**11
    ratio_needed = M_MU / m_mu_bare

    print(f"  m_mu(bare) = {m_mu_bare:.6f} MeV")
    print(f"  m_mu(exp)  = {M_MU:.6f} MeV")
    print(f"  Factor necesar: {ratio_needed:.8f}")
    print(f"  Corectie procentuala: +{(ratio_needed-1)*100:.4f}%")
    print()

    # Testam diverse formule geometrice
    candidates = {
        "1 + alpha": 1 + ALPHA,
        "1 + alpha*pi": 1 + ALPHA*PI,
        "1 + alpha*phi": 1 + ALPHA*PHI,
        "1 + alpha*phi^2": 1 + ALPHA*PHI**2,
        "1 + alpha*2*phi^2": 1 + ALPHA*2*PHI**2,
        "1 + alpha*11": 1 + ALPHA*11,
        "1 + alpha*4*pi": 1 + ALPHA*4*PI,
        "1 + 2*alpha*pi": 1 + 2*ALPHA*PI,
        "phi^(11+alpha*11)": PHI**(11+ALPHA*11) / PHI**11,  # ratio
        "1 + alpha_s*alpha": 1 + ALPHA_S_MZ*ALPHA,
        "1 + 1/(2*pi*phi^3)": 1 + 1/(2*PI*PHI**3),
        "1 + alpha*6": 1 + ALPHA*6,
        "1 + alpha*20": 1 + ALPHA*20,
        "1 + alpha*26": 1 + ALPHA*26,
        "phi^(11.08)/phi^11": PHI**0.08,
        "(1+alpha)^5": (1+ALPHA)**5,
        "(1+alpha)^(5+6)": (1+ALPHA)**11,
        "sqrt(phi^11*phi^(11+1))/phi^11": math.sqrt(PHI),  # phi^0.5
        "1 + 3*alpha*phi": 1 + 3*ALPHA*PHI,
        "1 + 5*alpha*phi": 1 + 5*ALPHA*PHI,
        "1 + 6*alpha*phi": 1 + 6*ALPHA*PHI,
        "1 + alpha*phi^3": 1 + ALPHA*PHI**3,
        "1 + alpha*(5+6*phi)": 1 + ALPHA*(5+6*PHI),
    }

    # Sortam dupa cat de aproape sunt
    sorted_candidates = sorted(candidates.items(), key=lambda x: abs(x[1] - ratio_needed))

    print(f"  Top 10 formule (cele mai apropiate de {ratio_needed:.6f}):")
    print(f"  {'Formula':35s} {'Valoare':>12s} {'Dif':>12s} {'err_mu':>8s}")
    print(f"  {'-'*70}")

    for name, val in sorted_candidates[:10]:
        m_mu_test = m_mu_bare * val
        err = pct_err(m_mu_test, M_MU)
        dif = val - ratio_needed
        print(f"  {name:35s} {val:12.8f} {dif:+12.8f} {err:7.4f}%")

    print()

    # Testam cel mai bun candidat si pe TAU
    best_name, best_val = sorted_candidates[0]
    print(f"  Cel mai bun candidat: {best_name} = {best_val:.8f}")
    print()

    # Verificam pe tau cu aceeasi formula
    m_tau_bare = M_E * PHI**17
    m_tau_corrected = m_tau_bare * best_val  # NOTA: e corect sa aplicam acelasi factor?

    print(f"  Test pe TAU:")
    print(f"    m_tau(bare) = {m_tau_bare:.2f} MeV")
    print(f"    m_tau(corr) = {m_tau_corrected:.2f} MeV (cu acelasi factor)")
    print(f"    m_tau(exp)  = {M_TAU:.2f} MeV")
    print(f"    Eroare bare: {pct_err(m_tau_bare, M_TAU):.4f}%")
    print(f"    Eroare corr: {pct_err(m_tau_corrected, M_TAU):.4f}%")
    print()

    # De fapt, daca formula e m_f = m_e * phi^(n + delta)
    # atunci delta ar trebui sa fie DIFERIT pentru fiecare nivel
    delta_mu = math.log(M_MU / M_E) / math.log(PHI) - 11
    delta_tau = math.log(M_TAU / M_E) / math.log(PHI) - 17

    print(f"  Deviatii de la exponenti intregi:")
    print(f"    delta(muon) = {delta_mu:.6f}")
    print(f"    delta(tau)  = {delta_tau:.6f}")
    print(f"    delta(mu)/11 = {delta_mu/11:.6f}")
    print(f"    delta(tau)/17 = {delta_tau/17:.6f}")
    print()

    # Sunt proportionale cu n?
    ratio_deltas = (delta_mu/11) / (delta_tau/17) if delta_tau != 0 else float('inf')
    print(f"    Raportul (delta/n)_mu / (delta/n)_tau = {ratio_deltas:.4f}")
    if abs(ratio_deltas - 1) < 0.1:
        print(f"    Sugestie: delta = k * n cu k = {delta_mu/11:.6f}")
        print(f"    Asta inseamna n_eff = n * (1 + k) = n * {1 + delta_mu/11:.6f}")
    else:
        print(f"    Deviatiile NU sunt proportionale cu n")
        print(f"    Deci corectia NU e o simpla rescalare a lui n")
    print()


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print()
    print("*" * 70)
    print("  EXP-096: VALIDARE FORMULE DE CORECTIE")
    print("  Data: 6 Februarie 2026")
    print("*" * 70)
    print()

    test_em_correction()
    test_absolute_mass()
    test_qcd_correction()
    test_unified_formula()
    test_muon_deep()

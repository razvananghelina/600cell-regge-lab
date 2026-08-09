"""
EXP-097: Cautarea Corectiei Simultane Muon + Tau
=================================================
Problema: phi^11 si phi^6 nu dau exact masele.
  - m_mu/m_e: phi^11 = 199.0 vs 206.77 (nevoie +3.9%)
  - m_tau/m_mu: phi^6 = 17.94 vs 16.82 (nevoie -6.2%)
  - m_tau/m_e: phi^17 = 3571 vs 3477 (nevoie -2.6%)

Semnele OPUSE inseamna ca nu exista un factor universal simplu.

STRATEGIE: Decomponem in d (diametru=5) si k (decagoane=6).
  n_mu    = d + k = 11
  n_tau/mu= k    = 6
  n_tau   = d + 2k = 17

Cautam corectii delta_d si delta_k din geometria 600-cell.
"""

import math
from itertools import product

PHI = (1 + math.sqrt(5)) / 2
PI = math.pi
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

# Mase experimentale (MeV)
M_E = 0.51099895
M_MU = 105.6583755
M_TAU = 1776.86

def pct_err(calc, exp):
    return abs(calc - exp) / exp * 100

# =============================================================================
# PARTEA 1: DECOMPOZIȚIA delta_d si delta_k
# =============================================================================
def decompose_corrections():
    print("=" * 70)
    print("PARTEA 1: DECOMPOZITIA IN delta_d SI delta_k")
    print("=" * 70)
    print()

    # Exponenti exacti
    n_mu_exact = math.log(M_MU / M_E) / math.log(PHI)   # m_mu/m_e
    n_tau_mu_exact = math.log(M_TAU / M_MU) / math.log(PHI)  # m_tau/m_mu
    n_tau_exact = math.log(M_TAU / M_E) / math.log(PHI)  # m_tau/m_e

    print(f"  Exponenti exacti continui:")
    print(f"    n(mu/e)    = {n_mu_exact:.6f}  (intreg: 11)")
    print(f"    n(tau/mu)  = {n_tau_mu_exact:.6f}  (intreg: 6)")
    print(f"    n(tau/e)   = {n_tau_exact:.6f}  (intreg: 17)")
    print(f"    Verificare: {n_mu_exact:.6f} + {n_tau_mu_exact:.6f} = {n_mu_exact + n_tau_mu_exact:.6f} vs {n_tau_exact:.6f}")
    print()

    # Decompozitie: n_mu = d + k, n_tau/mu = k, n_tau = d + 2k
    # d = 5, k = 6
    # delta: n_mu = (5+dd) + (6+dk) = 11 + dd + dk
    #        n_tau/mu = 6 + dk
    #        n_tau = (5+dd) + 2*(6+dk) = 17 + dd + 2*dk

    # Din n_tau/mu:
    delta_k = n_tau_mu_exact - 6
    # Din n_mu:
    delta_d = n_mu_exact - 11 - delta_k
    # Verificare n_tau:
    n_tau_check = 17 + delta_d + 2*delta_k

    print(f"  Corectii individuale:")
    print(f"    delta_d (diametru) = {delta_d:+.6f}")
    print(f"    delta_k (decagon)  = {delta_k:+.6f}")
    print()
    print(f"  Verificare:")
    print(f"    d_eff + k_eff = {5+delta_d:.6f} + {6+delta_k:.6f} = {11+delta_d+delta_k:.6f}  vs n_mu = {n_mu_exact:.6f}")
    print(f"    k_eff = {6+delta_k:.6f}  vs n_tau/mu = {n_tau_mu_exact:.6f}")
    print(f"    d_eff + 2*k_eff = {17+delta_d+2*delta_k:.6f}  vs n_tau = {n_tau_exact:.6f}")
    print()

    # Ce sunt delta_d si delta_k?
    print(f"  PROBLEMA: Ce sunt delta_d={delta_d:+.6f} si delta_k={delta_k:+.6f}?")
    print(f"    delta_d > 0: diametrul e putin MAI MARE decat 5")
    print(f"    delta_k < 0: contributia decagonala e putin MAI MICA decat 6")
    print()

    return delta_d, delta_k


# =============================================================================
# PARTEA 2: CAUTARE EXPRESII GEOMETRICE PENTRU delta_d si delta_k
# =============================================================================
def search_geometric_expressions(delta_d, delta_k):
    print("=" * 70)
    print("PARTEA 2: CAUTARE EXPRESII GEOMETRICE")
    print("=" * 70)
    print()

    # Biblioteca de expresii candidate bazate pe 600-cell
    expressions = {
        # Constante de cuplaj
        "alpha": ALPHA,
        "-alpha": -ALPHA,
        "alpha_s": ALPHA_S,
        "-alpha_s": -ALPHA_S,
        "sin2_tW": SIN2_TW,
        "-sin2_tW": -SIN2_TW,

        # Puteri ale lui phi
        "1/phi": 1/PHI,
        "-1/phi": -1/PHI,
        "1/phi^2": 1/PHI**2,
        "-1/phi^2": -1/PHI**2,
        "1/phi^3": 1/PHI**3,
        "-1/phi^3": -1/PHI**3,
        "1/phi^4": 1/PHI**4,
        "-1/phi^4": -1/PHI**4,
        "1/phi^5": 1/PHI**5,
        "-1/phi^5": -1/PHI**5,

        # Combinatii cu alpha
        "5*alpha": 5*ALPHA,
        "6*alpha": 6*ALPHA,
        "12*alpha": 12*ALPHA,
        "20*alpha": 20*ALPHA,
        "26*alpha": 26*ALPHA,
        "-5*alpha": -5*ALPHA,
        "-6*alpha": -6*ALPHA,
        "-12*alpha": -12*ALPHA,
        "-20*alpha": -20*ALPHA,

        # Combinatii cu phi si alpha
        "alpha*phi": ALPHA*PHI,
        "alpha*phi^2": ALPHA*PHI**2,
        "alpha*2*phi^2": ALPHA*2*PHI**2,
        "-alpha*phi": -ALPHA*PHI,
        "-alpha*phi^2": -ALPHA*PHI**2,

        # Combinatii cu alpha_s
        "alpha_s/pi": ALPHA_S/PI,
        "-alpha_s/pi": -ALPHA_S/PI,
        "alpha_s/phi": ALPHA_S/PHI,
        "-alpha_s/phi": -ALPHA_S/PHI,
        "alpha_s*alpha": ALPHA_S*ALPHA,

        # Rapoarte geometrice
        "1/12": 1/12,
        "-1/12": -1/12,
        "1/6": 1/6,
        "-1/6": -1/6,
        "1/5": 1/5,
        "-1/5": -1/5,
        "1/20": 1/20,
        "-1/20": -1/20,
        "1/26": 1/26,
        "-1/26": -1/26,
        "2/26": 2/26,
        "-2/26": -2/26,
        "6/120": 6/120,
        "-6/120": -6/120,

        # Trigonometrie pe 600-cell
        "sin(pi/5)": math.sin(PI/5),
        "-sin(pi/5)": -math.sin(PI/5),
        "cos(pi/5)": math.cos(PI/5),
        "1-cos(pi/5)": 1-math.cos(PI/5),
        "sin(pi/5)-1/2": math.sin(PI/5)-0.5,
        "sin(2*pi/5)-1": math.sin(2*PI/5)-1,
        "pi/5 - 0.6": PI/5 - 0.6,

        # Euler / topologice
        "chi=0": 0,

        # Spectrale
        "1/(2*phi^3)": 1/(2*PHI**3),
        "-1/(2*phi^3)": -1/(2*PHI**3),
        "alpha/(1-2*alpha)": ALPHA/(1-2*ALPHA),
        "alpha*pi/phi": ALPHA*PI/PHI,
        "-alpha*pi/phi": -ALPHA*PI/PHI,

        # Combinatii dimensionale
        "1/(4*pi)": 1/(4*PI),
        "-1/(4*pi)": -1/(4*PI),
        "alpha*4*pi": ALPHA*4*PI,
        "1/(2*pi*phi)": 1/(2*PI*PHI),
        "-1/(2*pi*phi)": -1/(2*PI*PHI),
    }

    # Cautam cele mai bune potriviri pentru delta_d
    print(f"  --- delta_d = {delta_d:+.6f} (corectia diametrului) ---")
    matches_d = sorted(expressions.items(), key=lambda x: abs(x[1] - delta_d))
    print(f"  {'Expresie':30s} {'Valoare':>12s} {'Dif':>12s} {'Err%':>8s}")
    print(f"  {'-'*65}")
    for name, val in matches_d[:15]:
        err = abs(val - delta_d) / abs(delta_d) * 100 if delta_d != 0 else 999
        print(f"  {name:30s} {val:+12.6f} {val-delta_d:+12.6f} {err:7.2f}%")

    print()

    # Cautam cele mai bune potriviri pentru delta_k
    print(f"  --- delta_k = {delta_k:+.6f} (corectia decagonala) ---")
    matches_k = sorted(expressions.items(), key=lambda x: abs(x[1] - delta_k))
    print(f"  {'Expresie':30s} {'Valoare':>12s} {'Dif':>12s} {'Err%':>8s}")
    print(f"  {'-'*65}")
    for name, val in matches_k[:15]:
        err = abs(val - delta_k) / abs(delta_k) * 100 if delta_k != 0 else 999
        print(f"  {name:30s} {val:+12.6f} {val-delta_k:+12.6f} {err:7.2f}%")

    print()

    return matches_d[:5], matches_k[:5]


# =============================================================================
# PARTEA 3: CAUTARE COMBINATORIE - PERECHI (delta_d, delta_k)
# =============================================================================
def combinatorial_search():
    print("=" * 70)
    print("PARTEA 3: CAUTARE COMBINATORIE - PERECHI (delta_d, delta_k)")
    print("=" * 70)
    print()

    # Construim expresii candidate
    building_blocks = {
        "alpha": ALPHA,
        "alpha_s": ALPHA_S,
        "1/phi": 1/PHI,
        "1/phi^2": 1/PHI**2,
        "1/phi^3": 1/PHI**3,
        "1/phi^4": 1/PHI**4,
        "1/phi^5": 1/PHI**5,
        "sin2_tW": SIN2_TW,
        "1/12": 1/12,
        "1/5": 1/5,
        "1/6": 1/6,
        "pi/5": PI/5,
        "sin(pi/5)": math.sin(PI/5),
        "alpha*phi": ALPHA*PHI,
        "alpha*phi^2": ALPHA*PHI**2,
    }

    # Multiplicatori
    multipliers = {
        "": 1,
        "-": -1,
        "2*": 2,
        "-2*": -2,
        "5*": 5,
        "-5*": -5,
        "6*": 6,
        "-6*": -6,
    }

    # Generare candidati
    candidates_d = {}
    candidates_k = {}
    for (m_name, m_val), (b_name, b_val) in product(multipliers.items(), building_blocks.items()):
        val = m_val * b_val
        name = f"{m_name}{b_name}"
        if -0.5 < val < 0.5:  # range rezonabil
            candidates_d[name] = val
            candidates_k[name] = val

    # Testam TOATE perechile
    n_mu_target = math.log(M_MU / M_E) / math.log(PHI)
    n_tau_mu_target = math.log(M_TAU / M_MU) / math.log(PHI)
    n_tau_target = math.log(M_TAU / M_E) / math.log(PHI)

    best_pairs = []

    for (d_name, d_val), (k_name, k_val) in product(candidates_d.items(), candidates_k.items()):
        d_eff = 5 + d_val
        k_eff = 6 + k_val

        n_mu = d_eff + k_eff
        n_tau_mu = k_eff
        n_tau = d_eff + 2 * k_eff

        m_mu_calc = M_E * PHI**n_mu
        m_tau_calc = M_E * PHI**n_tau

        err_mu = pct_err(m_mu_calc, M_MU)
        err_tau = pct_err(m_tau_calc, M_TAU)
        total_err = err_mu + err_tau

        if total_err < 2.0:  # sub 2% total
            best_pairs.append((d_name, k_name, d_val, k_val, err_mu, err_tau, total_err))

    best_pairs.sort(key=lambda x: x[6])

    print(f"  Perechi cu eroare totala < 2%:")
    print(f"  {'delta_d':25s} {'delta_k':25s} {'d_val':>8s} {'k_val':>8s} {'err_mu':>7s} {'err_tau':>7s} {'TOTAL':>7s}")
    print(f"  {'-'*105}")

    for d_name, k_name, d_val, k_val, err_mu, err_tau, total in best_pairs[:20]:
        print(f"  {d_name:25s} {k_name:25s} {d_val:+8.4f} {k_val:+8.4f} {err_mu:6.3f}% {err_tau:6.3f}% {total:6.3f}%")

    print(f"\n  Total perechi gasite: {len(best_pairs)}")
    print()

    if best_pairs:
        return best_pairs[0]
    return None


# =============================================================================
# PARTEA 4: ANALIZA FIZICA A CELEI MAI BUNE PERECHI
# =============================================================================
def analyze_best_pair(best):
    if best is None:
        print("  Nu s-a gasit nicio pereche buna!")
        return

    d_name, k_name, d_val, k_val, err_mu, err_tau, total = best

    print("=" * 70)
    print("PARTEA 4: ANALIZA CELEI MAI BUNE PERECHI")
    print("=" * 70)
    print()

    d_eff = 5 + d_val
    k_eff = 6 + k_val

    print(f"  delta_d = {d_name} = {d_val:+.6f}")
    print(f"  delta_k = {k_name} = {k_val:+.6f}")
    print()
    print(f"  d_eff = 5 + ({d_val:+.6f}) = {d_eff:.6f}")
    print(f"  k_eff = 6 + ({k_val:+.6f}) = {k_eff:.6f}")
    print()

    # Calcule complete
    n_mu = d_eff + k_eff
    n_tau_mu = k_eff
    n_tau = d_eff + 2 * k_eff

    m_mu_calc = M_E * PHI**n_mu
    m_tau_calc = M_E * PHI**n_tau
    ratio_tau_mu = PHI**n_tau_mu

    print(f"  Exponenti efecivi:")
    print(f"    n_mu = d_eff + k_eff = {n_mu:.6f}")
    print(f"    n_tau/mu = k_eff = {n_tau_mu:.6f}")
    print(f"    n_tau = d_eff + 2*k_eff = {n_tau:.6f}")
    print()

    print(f"  Mase:")
    print(f"    m_mu  = m_e * phi^{n_mu:.4f} = {m_mu_calc:.4f} MeV  (exp: {M_MU:.4f}, err: {pct_err(m_mu_calc, M_MU):.4f}%)")
    print(f"    m_tau = m_e * phi^{n_tau:.4f} = {m_tau_calc:.2f} MeV  (exp: {M_TAU:.2f}, err: {pct_err(m_tau_calc, M_TAU):.4f}%)")
    print(f"    m_tau/m_mu = phi^{n_tau_mu:.4f} = {ratio_tau_mu:.4f}  (exp: {M_TAU/M_MU:.4f}, err: {pct_err(ratio_tau_mu, M_TAU/M_MU):.4f}%)")
    print()

    # Interpretare fizica
    print(f"  INTERPRETARE FIZICA:")
    print(f"  d_eff = 5 + {d_name}")
    print(f"  k_eff = 6 + {k_name}")
    print()


# =============================================================================
# PARTEA 5: ABORDARE ALTERNATIVA - FORMULE MULTIPLICATIVE
# =============================================================================
def test_multiplicative():
    print("=" * 70)
    print("PARTEA 5: ABORDARE MULTIPLICATIVA")
    print("m_l = m_e * phi^n * f(n, geometrie)")
    print("=" * 70)
    print()

    # In loc de n_eff = n + delta, incercam:
    # m = m_e * phi^n * (1 + c_1)^p1 * (1 + c_2)^p2
    # unde p1, p2 depind de generatie

    # Muon: n=11 = 5 + 6, cu d=5 "traversari de diametru" si k=6 "bucle decagonale"
    # Tau: n=17 = 5 + 2*6, cu d=5 si k=12 (doua seturi de bucle)

    # IDEA: factorul e diferit pentru componentele d si k
    # m_mu = m_e * phi^5 * F_d * phi^6 * F_k
    # m_tau = m_e * phi^5 * F_d * phi^12 * F_k^2  (doua aplicari ale buclei)
    # sau m_tau = m_e * phi^5 * F_d * (phi^6 * F_k)^2

    # Ecuatii:
    # m_mu/m_e = phi^5 * F_d * phi^6 * F_k = phi^11 * F_d * F_k
    # m_tau/m_e = phi^5 * F_d * phi^12 * F_k^2 = phi^17 * F_d * F_k^2
    # m_tau/m_mu = phi^6 * F_k

    # Din m_tau/m_mu:
    ratio_tau_mu_exp = M_TAU / M_MU
    F_k = ratio_tau_mu_exp / PHI**6
    print(f"  Din m_tau/m_mu:")
    print(f"    F_k = (m_tau/m_mu) / phi^6 = {ratio_tau_mu_exp:.6f} / {PHI**6:.6f} = {F_k:.8f}")
    print()

    # Din m_mu/m_e:
    ratio_mu_e_exp = M_MU / M_E
    F_d = ratio_mu_e_exp / (PHI**11 * F_k)
    print(f"  Din m_mu/m_e:")
    print(f"    F_d = (m_mu/m_e) / (phi^11 * F_k) = {ratio_mu_e_exp:.4f} / ({PHI**11:.4f} * {F_k:.6f}) = {F_d:.8f}")
    print()

    # Verificare pe tau:
    m_tau_check = M_E * PHI**17 * F_d * F_k**2
    print(f"  Verificare m_tau = m_e * phi^17 * F_d * F_k^2:")
    print(f"    = {M_E:.6f} * {PHI**17:.2f} * {F_d:.6f} * {F_k**2:.6f}")
    print(f"    = {m_tau_check:.2f} MeV  (exp: {M_TAU:.2f}, err: {pct_err(m_tau_check, M_TAU):.6f}%)")
    print()

    # Ce sunt F_d si F_k?
    print(f"  VALORILE FACTORILOR:")
    print(f"    F_d = {F_d:.8f}  (corectie diametru: {(F_d-1)*100:+.4f}%)")
    print(f"    F_k = {F_k:.8f}  (corectie decagon: {(F_k-1)*100:+.4f}%)")
    print()

    # Cautam expresii
    print(f"  CAUTAM EXPRESII GEOMETRICE:")
    print()

    # F_d
    fd_candidates = {
        "1 + alpha*phi^3": 1 + ALPHA*PHI**3,
        "1 + alpha*5*phi": 1 + ALPHA*5*PHI,
        "1 + alpha*26": 1 + ALPHA*26,
        "1 + 5*alpha*phi": 1 + 5*ALPHA*PHI,
        "1 + alpha*4*pi": 1 + ALPHA*4*PI,
        "1 + sin(pi/5)/5": 1 + math.sin(PI/5)/5,
        "(phi^5+1)/phi^5": (PHI**5+1)/PHI**5,
        "1 + 1/(5*phi^3)": 1 + 1/(5*PHI**3),
        "1 + alpha*(5+6)": 1 + ALPHA*11,
        "1 + alpha*12": 1 + ALPHA*12,
        "1 + 2*alpha*pi": 1 + 2*ALPHA*PI,
        "1 + alpha*(6+20)": 1 + ALPHA*26,
        "phi^(1/5)": PHI**(1/5),
        "1 + alpha*20": 1 + ALPHA*20,
        "1 + alpha*2*phi^2": 1 + ALPHA*2*PHI**2,
    }

    print(f"  F_d = {F_d:.8f}")
    sorted_fd = sorted(fd_candidates.items(), key=lambda x: abs(x[1] - F_d))
    for name, val in sorted_fd[:10]:
        err = abs(val - F_d) / abs(F_d - 1) * 100 if F_d != 1 else 999
        # Testam pe ambele
        m_mu_test = M_E * PHI**11 * val * F_k  # F_k e exact din definitie
        m_tau_test = M_E * PHI**17 * val * F_k**2
        print(f"    {name:25s} = {val:.8f}  dif={(val-F_d):+.6f}  mu_err={pct_err(m_mu_test, M_MU):.3f}%  tau_err={pct_err(m_tau_test, M_TAU):.3f}%")
    print()

    # F_k
    fk_candidates = {
        "1 - alpha_s/pi": 1 - ALPHA_S/PI,
        "1 - 1/phi^4": 1 - 1/PHI**4,
        "1 - 6*alpha": 1 - 6*ALPHA,
        "1 - alpha*6*phi": 1 - ALPHA*6*PHI,
        "1 - alpha_s/2*phi": 1 - ALPHA_S/(2*PHI),
        "1 - 1/(6*phi)": 1 - 1/(6*PHI),
        "(phi^6-1)/phi^6": (PHI**6-1)/PHI**6,
        "1 - alpha*12": 1 - ALPHA*12,
        "1 - sin(pi/5)/6": 1 - math.sin(PI/5)/6,
        "phi^(-1/6)": PHI**(-1/6),
        "1 - 1/12": 1 - 1/12,
        "1 - alpha*2*phi^2": 1 - ALPHA*2*PHI**2,
        "cos(pi/5)/phi": math.cos(PI/5)/PHI,
        "1 - alpha*phi^3": 1 - ALPHA*PHI**3,
        "1 - 5*alpha": 1 - 5*ALPHA,
    }

    print(f"  F_k = {F_k:.8f}")
    sorted_fk = sorted(fk_candidates.items(), key=lambda x: abs(x[1] - F_k))
    for name, val in sorted_fk[:10]:
        # Testam cu F_d fix
        m_mu_test = M_E * PHI**11 * F_d * val
        m_tau_test = M_E * PHI**17 * F_d * val**2
        print(f"    {name:25s} = {val:.8f}  dif={(val-F_k):+.6f}  mu_err={pct_err(m_mu_test, M_MU):.3f}%  tau_err={pct_err(m_tau_test, M_TAU):.3f}%")
    print()

    # CAUTARE COMBINATA: gasim PERECHEA (F_d_expr, F_k_expr) care minimizeaza eroarea totala
    print("  --- CAUTARE COMBINATA: Cele mai bune perechi (F_d, F_k) ---")
    print()

    best_combo = []
    for (fd_name, fd_val), (fk_name, fk_val) in product(fd_candidates.items(), fk_candidates.items()):
        m_mu_test = M_E * PHI**11 * fd_val * fk_val
        m_tau_test = M_E * PHI**17 * fd_val * fk_val**2
        err_mu = pct_err(m_mu_test, M_MU)
        err_tau = pct_err(m_tau_test, M_TAU)
        total = err_mu + err_tau
        if total < 1.0:
            best_combo.append((fd_name, fk_name, fd_val, fk_val, err_mu, err_tau, total))

    best_combo.sort(key=lambda x: x[6])

    print(f"  Perechi cu eroare totala < 1%:")
    print(f"  {'F_d':25s} {'F_k':25s} {'err_mu':>8s} {'err_tau':>8s} {'TOTAL':>8s}")
    print(f"  {'-'*90}")
    for fd_name, fk_name, fd_val, fk_val, err_mu, err_tau, total in best_combo[:15]:
        print(f"  {fd_name:25s} {fk_name:25s} {err_mu:7.3f}% {err_tau:7.3f}% {total:7.3f}%")

    print(f"\n  Perechi sub 1%: {len(best_combo)}")
    print()

    if best_combo:
        return best_combo[0], F_d, F_k
    return None, F_d, F_k


# =============================================================================
# PARTEA 6: FORMULA FINALA SI VALIDARE COMPLETA
# =============================================================================
def final_validation(best_combo, F_d_exact, F_k_exact):
    print("=" * 70)
    print("PARTEA 6: FORMULA FINALA SI VALIDARE")
    print("=" * 70)
    print()

    if best_combo is None:
        print("  Nu s-a gasit nicio pereche sub 1%. Folosim valorile exacte.")
        fd_name, fk_name = "EXACT", "EXACT"
        fd_val, fk_val = F_d_exact, F_k_exact
    else:
        fd_name, fk_name, fd_val, fk_val = best_combo[0], best_combo[1], best_combo[2], best_combo[3]

    print(f"  Formula: m_lepton = m_e * phi^(5*p_d + 6*p_k) * F_d^p_d * F_k^p_k")
    print()
    print(f"  F_d ({fd_name}) = {fd_val:.8f}")
    print(f"  F_k ({fk_name}) = {fk_val:.8f}")
    print()

    # Puterile
    leptons = [
        ("Electron", 0, 0, M_E),
        ("Muon",     1, 1, M_MU),
        ("Tau",      1, 2, M_TAU),
    ]

    print(f"  {'Particula':10s} {'p_d':>3s} {'p_k':>3s} {'n_eff':>7s} {'m_calc':>12s} {'m_exp':>12s} {'err':>8s}")
    print(f"  {'-'*60}")

    for name, pd, pk, m_exp in leptons:
        n_phi = 5*pd + 6*pk
        m_calc = M_E * PHI**n_phi * fd_val**pd * fk_val**pk
        err = pct_err(m_calc, m_exp) if name != "Electron" else 0.0
        print(f"  {name:10s} {pd:3d} {pk:3d} {n_phi:7d} {m_calc:12.4f} {m_exp:12.4f} {err:7.4f}%")

    print()

    # Interpretare
    print(f"  INTERPRETARE:")
    print(f"    Electron (p_d=0, p_k=0): starea fundamentala (n=0)")
    print(f"    Muon (p_d=1, p_k=1):     o traversare de diametru + o bucla decagonala")
    print(f"    Tau (p_d=1, p_k=2):      o traversare de diametru + doua bucle decagonale")
    print()
    print(f"    F_d = corectia de faza acumulata la traversarea diametrului (5 pasi)")
    print(f"    F_k = corectia de faza acumulata pe o bucla decagonala (6 laturi)")
    print()

    # Rapoarte
    print(f"  RAPOARTE:")
    m_mu_calc = M_E * PHI**11 * fd_val * fk_val
    m_tau_calc = M_E * PHI**17 * fd_val * fk_val**2
    print(f"    m_mu/m_e  = {m_mu_calc/M_E:.4f}  (exp: {M_MU/M_E:.4f}, err: {pct_err(m_mu_calc/M_E, M_MU/M_E):.4f}%)")
    print(f"    m_tau/m_mu = {m_tau_calc/m_mu_calc:.4f}  (exp: {M_TAU/M_MU:.4f}, err: {pct_err(m_tau_calc/m_mu_calc, M_TAU/M_MU):.4f}%)")
    print(f"    m_tau/m_e  = {m_tau_calc/M_E:.2f}  (exp: {M_TAU/M_E:.2f}, err: {pct_err(m_tau_calc/M_E, M_TAU/M_E):.4f}%)")
    print()


# =============================================================================
# PARTEA 7: ABORDARE DIN PRIMA PRINCIPII - HOLONOMIE REALA
# =============================================================================
def holonomy_approach():
    print("=" * 70)
    print("PARTEA 7: HOLONOMIE PE S^3 - CALCUL EXACT")
    print("=" * 70)
    print()

    # Pe S^3, holonomia unui spinor pe un triunghi geodezic
    # cu arii A1, A2, A3 da o faza = A/2 (pentru spinori pe S^3)
    # Unde A e aria solida pe sfera

    # 600-cell inscris in S^3 cu raza R=1
    # Unghi muchie = pi/5 (EXACT)
    # Un decagon pe 600-cell = 10 muchii, formand un cerc mare

    edge_angle = PI/5  # 36 grade
    print(f"  Unghiul muchie 600-cell: pi/5 = {math.degrees(edge_angle):.1f} grade")
    print()

    # Diametrul: 5 pasi x pi/5 = pi
    # Asta e o traversare completa de la nord la sud pe S^3
    diameter_phase = 5 * edge_angle
    print(f"  Faza diametru: 5 * pi/5 = {diameter_phase/PI:.4f} * pi = {math.degrees(diameter_phase):.1f} grade")
    print(f"  = pi (antipod pe S^3)")
    print()

    # Decagon: 10 pasi, dar fiecare vertex vede 6 decagoane
    # O bucla completa pe un decagon = 10 * pi/5 = 2*pi
    decagon_phase = 10 * edge_angle
    print(f"  Faza decagon complet: 10 * pi/5 = {decagon_phase/PI:.4f} * pi = {math.degrees(decagon_phase):.1f} grade")
    print(f"  = 2*pi (bucla completa)")
    print()

    # Holonomia spinoriala pe S^3:
    # Un spinor transportat pe un cerc mare de S^3 acumuleaza faza -pi (semnul Berry)
    # Deci fiecare decagon (cerc mare) da o faza Berry de -pi

    # Masa ca exponentiala a holonomiei:
    # m = m_e * exp(sum of phases / fundamental_phase)
    # = m_e * exp(n * log(phi))
    # = m_e * phi^n

    # Dar phi = exp(log(phi)) = exp(0.4812...)
    # Si pi/5 = 0.6283...
    # Raportul: log(phi) / (pi/5) = 0.4812 / 0.6283 = 0.7659

    log_phi = math.log(PHI)
    print(f"  log(phi) = {log_phi:.6f}")
    print(f"  pi/5 = {edge_angle:.6f}")
    print(f"  Raportul log(phi)/(pi/5) = {log_phi/edge_angle:.6f}")
    print()

    # Idea: masa = m_e * exp(fase_totala * raport)
    # Muon: faza totala = pi (diametru) + 6*(????)
    # Trebuie sa intelegem ce fractiune din decagon contribuie

    # Holonomia exacta pe S^3 pentru un triunghi sferic cu laturi a, b, c
    # este Omega = A + B + C - pi, unde A, B, C sunt unghiurile
    # Pe 600-cell, triunghiurile sunt echilaterale cu unghi diedru

    # Unghi diedru tetrahedral pe 600-cell: arccos(1/4) ≈ 75.52 grade
    # Dar asta e unghiul dintre fetele 3D

    # Sa incercam o abordare mai directa:
    # Daca faza per pas pe diametru = log(phi) + delta_d_per_step
    # si faza per latura de decagon = log(phi) + delta_k_per_step

    # Muon: 5*(log_phi + eps_d) + 6*(log_phi + eps_k) = log(M_MU/M_E)
    # Tau: 5*(log_phi + eps_d) + 12*(log_phi + eps_k) = log(M_TAU/M_E)

    log_mu = math.log(M_MU / M_E)
    log_tau = math.log(M_TAU / M_E)

    # 5*eps_d + 6*eps_k = log_mu - 11*log_phi
    # 5*eps_d + 12*eps_k = log_tau - 17*log_phi

    rhs_mu = log_mu - 11 * log_phi
    rhs_tau = log_tau - 17 * log_phi

    # Sistem liniar:
    # 5*eps_d + 6*eps_k = rhs_mu
    # 5*eps_d + 12*eps_k = rhs_tau
    # Scazand: 6*eps_k = rhs_tau - rhs_mu
    eps_k = (rhs_tau - rhs_mu) / 6
    eps_d = (rhs_mu - 6*eps_k) / 5

    print(f"  SISTEM LINIAR:")
    print(f"    5*eps_d + 6*eps_k  = log(m_mu/m_e) - 11*log(phi) = {rhs_mu:+.8f}")
    print(f"    5*eps_d + 12*eps_k = log(m_tau/m_e) - 17*log(phi) = {rhs_tau:+.8f}")
    print()
    print(f"  SOLUTIE:")
    print(f"    eps_d = {eps_d:+.8f}  (corectie de faza per pas de diametru)")
    print(f"    eps_k = {eps_k:+.8f}  (corectie de faza per latura de decagon)")
    print()

    # Faza totala per pas
    phase_d = log_phi + eps_d
    phase_k = log_phi + eps_k

    print(f"  FAZA TOTALA PER PAS:")
    print(f"    Diametru: log(phi) + eps_d = {log_phi:.6f} + {eps_d:+.6f} = {phase_d:.6f}")
    print(f"    Decagon:  log(phi) + eps_k = {log_phi:.6f} + {eps_k:+.6f} = {phase_k:.6f}")
    print()

    # Ce sunt aceste faze?
    print(f"  CE SUNT ACESTE FAZE?")
    print(f"    Diametru: {phase_d:.6f} = ?")
    print(f"    Decagon:  {phase_k:.6f} = ?")
    print()

    # Candidati pentru phase_d
    print(f"  Candidati pentru faza de diametru ({phase_d:.6f}):")
    phase_candidates = {
        "log(phi)": log_phi,
        "pi/5 - log(2)": PI/5 - math.log(2),
        "pi/6": PI/6,
        "log(phi) + alpha": log_phi + ALPHA,
        "1/2": 0.5,
        "pi/5 * phi/2": PI/5 * PHI/2,
        "atan(1/phi)": math.atan(1/PHI),
        "log(phi+1)/2": math.log(PHI+1)/2,
        "log(phi) + alpha*pi/5": log_phi + ALPHA*PI/5,
        "pi/(5+1)": PI/6,
        "log(phi) + alpha^2*pi": log_phi + ALPHA**2*PI,
        "pi/(2*phi^2)": PI/(2*PHI**2),
        "asin(1/phi)/2": math.asin(1/PHI)/2,
        "log(phi)*(1+alpha)": log_phi*(1+ALPHA),
    }
    sorted_pd = sorted(phase_candidates.items(), key=lambda x: abs(x[1] - phase_d))
    for name, val in sorted_pd[:8]:
        print(f"      {name:30s} = {val:.8f}  dif = {val-phase_d:+.8f}")
    print()

    # Candidati pentru phase_k
    print(f"  Candidati pentru faza de decagon ({phase_k:.6f}):")
    phase_candidates_k = {
        "log(phi)": log_phi,
        "log(phi) - alpha": log_phi - ALPHA,
        "pi/5 - 1/2": PI/5 - 0.5,
        "log(phi)*phi/phi": log_phi,
        "pi/5 * (1-1/phi)": PI/5 * (1-1/PHI),
        "atan(1/2)": math.atan(0.5),
        "log(phi) - alpha_s/5": log_phi - ALPHA_S/5,
        "log(phi) - alpha*phi": log_phi - ALPHA*PHI,
        "2*log(phi) - log(phi+1/phi)": 2*log_phi - math.log(PHI + 1/PHI),
        "log(phi)*(1-alpha)": log_phi*(1-ALPHA),
        "pi/(2*phi^2+1)": PI/(2*PHI**2+1),
        "log(phi) - alpha*pi/5": log_phi - ALPHA*PI/5,
    }
    sorted_pk = sorted(phase_candidates_k.items(), key=lambda x: abs(x[1] - phase_k))
    for name, val in sorted_pk[:8]:
        print(f"      {name:30s} = {val:.8f}  dif = {val-phase_k:+.8f}")
    print()

    # Verificare finala
    print(f"  VERIFICARE FINALA:")
    m_mu_check = M_E * math.exp(5*phase_d + 6*phase_k)
    m_tau_check = M_E * math.exp(5*phase_d + 12*phase_k)
    print(f"    m_mu  = m_e * exp(5*phi_d + 6*phi_k)  = {m_mu_check:.4f} MeV  (exp: {M_MU:.4f})")
    print(f"    m_tau = m_e * exp(5*phi_d + 12*phi_k) = {m_tau_check:.2f} MeV  (exp: {M_TAU:.2f})")
    print()

    return eps_d, eps_k, phase_d, phase_k


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print()
    print("*" * 70)
    print("  EXP-097: CAUTAREA CORECTIEI SIMULTANE MUON + TAU")
    print("  Data: 6 Februarie 2026")
    print("*" * 70)
    print()

    delta_d, delta_k = decompose_corrections()
    search_geometric_expressions(delta_d, delta_k)
    best_pair = combinatorial_search()
    analyze_best_pair(best_pair)
    best_combo, F_d, F_k = test_multiplicative()
    final_validation(best_combo, F_d, F_k)
    holonomy_approach()

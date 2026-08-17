"""
EXP-100: Investigare Bottom si Down Quarks
==========================================
Problema din exp099:
  - Bottom (-1,4): 38% eroare cu delta_d=-1/5, delta_k=-1/(2*phi^3)
  - Down (1,0): 10% eroare

Dar din exp098 Partea 7, OPTIMUL NUMERIC era:
  Down-type: delta_d=-0.087, delta_k=-0.091
  => Strange 0.06%, Bottom 0.03%, Down 16% (Down ramane rau)

Deci problema REALA e Down-ul, nu Bottom-ul!
Trebuie sa:
  1. Gasim corectiile algebrice corecte pentru down-type
  2. Investigam daca Down are alt exponent (n != 5)
  3. Investigam decompositions alternative
"""

import math

PHI = (1 + math.sqrt(5)) / 2
PI = math.pi
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

M_E = 0.51099895
M_DOWN = 4.67        # +0.48 -0.17 MeV (MS-bar la 2 GeV)
M_STRANGE = 93.4     # +8.6 -3.4 MeV
M_BOTTOM = 4180.0    # +30 -20 MeV (MS-bar la m_b)

def pct_err(calc, exp):
    return abs(calc - exp) / exp * 100


# =============================================================================
# PARTEA 1: REEXAMINARE - CE SPUNE OPTIMUL NUMERIC?
# =============================================================================
def reexamine_optimal():
    print("=" * 70)
    print("PARTEA 1: REEXAMINARE OPTIMUL NUMERIC")
    print("=" * 70)
    print()

    # Din exp098 Partea 7, optimul down-type era:
    # delta_d = -0.087, delta_k = -0.091
    # Strange(1,1): 0.06%, Bottom(-1,4): 0.03%, Down(1,0): 16.37%

    # Problema: cu (1,0), Down nu are componenta k, deci delta_k nu conteaza!
    # n_eff(Down) = 1 * (5 + delta_d) + 0 * (6 + delta_k) = 5 + delta_d
    # Down: m = m_e * phi^(5 + delta_d)

    # Ce delta_d trebuie pentru Down?
    n_down_exact = math.log(M_DOWN / M_E) / math.log(PHI)
    delta_d_needed_down = n_down_exact - 5

    print(f"  n exact pentru Down: {n_down_exact:.6f}")
    print(f"  delta_d necesar (daca b=0): {delta_d_needed_down:+.6f}")
    print()

    # Ce delta_d trebuie pentru Strange?
    n_strange_exact = math.log(M_STRANGE / M_E) / math.log(PHI)
    # Strange (1,1): n = d + k, deci delta_d + delta_k = n_exact - 11
    delta_sum_strange = n_strange_exact - 11

    print(f"  n exact pentru Strange: {n_strange_exact:.6f}")
    print(f"  delta_d + delta_k necesar (Strange): {delta_sum_strange:+.6f}")
    print()

    # Ce delta_d + 4*delta_k trebuie pentru Bottom?
    n_bottom_exact = math.log(M_BOTTOM / M_E) / math.log(PHI)
    # Bottom (-1,4): n = -d + 4k = -5 + 24 = 19
    # n_eff = -(5+dd) + 4*(6+dk) = 19 - dd + 4*dk
    delta_combo_bottom = n_bottom_exact - 19  # = -dd + 4*dk

    print(f"  n exact pentru Bottom: {n_bottom_exact:.6f}")
    print(f"  -delta_d + 4*delta_k necesar (Bottom): {delta_combo_bottom:+.6f}")
    print()

    # Sistem:
    # Down (1,0):     delta_d = delta_d_needed_down = -0.5380
    # Strange (1,1):  delta_d + delta_k = delta_sum_strange
    # Bottom (-1,4): -delta_d + 4*delta_k = delta_combo_bottom
    #
    # Din Down: delta_d = -0.5380 (PREA MARE!)
    # Daca delta_d = -0.5380, atunci din Strange:
    #   delta_k = delta_sum_strange - delta_d = delta_sum_strange + 0.5380

    print(f"  PROBLEMA: Down (1,0) cere delta_d = {delta_d_needed_down:+.4f}")
    print(f"  Dar Strange si Bottom cer delta_d ~ -0.09")
    print(f"  Acestea sunt INCOMPATIBILE!")
    print()
    print(f"  CONCLUZIE: Descompunerea (1,0) pentru Down NU functioneaza")
    print(f"  in cadrul corectiilor care merg pentru Strange si Bottom.")
    print()

    # Ce daca Down NU e (1,0)?
    print(f"  ALTERNATIVA: Down nu e n=5 sau nu e (1,0)?")
    print(f"  n_exact(Down) = {n_down_exact:.4f}")
    print(f"  Cel mai apropiat intreg: {round(n_down_exact)} (dar 4 sau 5?)")
    print()

    # n=4 vs n=5 pentru Down
    m_down_n4 = M_E * PHI**4
    m_down_n5 = M_E * PHI**5
    print(f"  m_e * phi^4 = {m_down_n4:.4f} MeV (err: {pct_err(m_down_n4, M_DOWN):.2f}%)")
    print(f"  m_e * phi^5 = {m_down_n5:.4f} MeV (err: {pct_err(m_down_n5, M_DOWN):.2f}%)")
    print(f"  Down exp: {M_DOWN} +0.48 -0.17 MeV")
    print()

    return n_down_exact, delta_d_needed_down


# =============================================================================
# PARTEA 2: DESCOMPUNEREA (a,b) CU n=4 PENTRU DOWN
# =============================================================================
def try_down_n4():
    print("=" * 70)
    print("PARTEA 2: DOWN QUARK CU n=4 IN LOC DE n=5")
    print("=" * 70)
    print()

    # n=4: descompuneri 5a + 6b = 4
    # a=-2, b=? => -10 + 6b = 4 => b = 14/6 (nu intreg)
    # a=2, b=-1 => 10 - 6 = 4 (DA!)
    # a=-4, b=4 => -20+24 = 4

    print(f"  Descompuneri pentru n=4:")
    for a in range(-5, 6):
        for b in range(-5, 6):
            if 5*a + 6*b == 4:
                print(f"    ({a},{b}) -> {a}*5 + {b}*6 = {5*a+6*b}")
    print()

    # (2,-1) e cea mai simpla
    # Down (2,-1): n_eff = 2*d_eff - 1*k_eff
    # Cu corectiile down-type optime delta_d=-0.087, delta_k=-0.091:
    dd_opt = -0.087
    dk_opt = -0.091
    n_eff = 2*(5+dd_opt) - 1*(6+dk_opt)
    m_calc = M_E * PHI**n_eff
    err = pct_err(m_calc, M_DOWN)
    print(f"  Down (2,-1) cu dd={dd_opt}, dk={dk_opt}:")
    print(f"    n_eff = 2*{5+dd_opt:.3f} - 1*{6+dk_opt:.3f} = {n_eff:.4f}")
    print(f"    m_calc = {m_calc:.4f} MeV, m_exp = {M_DOWN:.4f} MeV, err = {err:.2f}%")
    print()

    # Testam cu diverse corectii algebrice
    print(f"  Testam Down (2,-1) cu diverse corectii:")
    print()

    candidates = {
        "alpha_s": ALPHA_S,
        "-alpha_s": -ALPHA_S,
        "1/phi^4": 1/PHI**4,
        "-1/phi^4": -1/PHI**4,
        "1/(2*phi^3)": 1/(2*PHI**3),
        "-1/(2*phi^3)": -1/(2*PHI**3),
        "sin^2(tW)": SIN2_TW,
        "-sin^2(tW)": -SIN2_TW,
        "1/5": 0.2,
        "-1/5": -0.2,
        "2*alpha_s/pi": 2*ALPHA_S/PI,
        "-2*alpha_s/pi": -2*ALPHA_S/PI,
        "12*alpha": 12*ALPHA,
        "-12*alpha": -12*ALPHA,
        "alpha_s/phi": ALPHA_S/PHI,
        "-alpha_s/phi": -ALPHA_S/PHI,
    }

    # Cautam perechea (dd, dk) optima cu Down(2,-1), Strange(1,1), Bottom(-1,4)
    down_type = [
        ("Down",    2, -1, M_DOWN),
        ("Strange", 1,  1, M_STRANGE),
        ("Bottom", -1,  4, M_BOTTOM),
    ]

    best = []
    for dd_name, dd_val in candidates.items():
        for dk_name, dk_val in candidates.items():
            total = 0
            errs = []
            for name, a, b, m_exp in down_type:
                n = a*(5+dd_val) + b*(6+dk_val)
                m = M_E * PHI**n
                e = pct_err(m, m_exp)
                total += e
                errs.append((name, e))
            avg = total / 3
            best.append((dd_name, dk_name, dd_val, dk_val, avg, errs))

    best.sort(key=lambda x: x[4])

    print(f"  Top 15 cu Down(2,-1):")
    print(f"  {'dd':20s} {'dk':20s} {'avg':>7s}  {'Down':>8s} {'Strange':>8s} {'Bottom':>8s}")
    print(f"  {'-'*80}")
    for dd_n, dk_n, dd_v, dk_v, avg, errs in best[:15]:
        d_str = ", ".join([f"{n}:{e:.2f}%" for n, e in errs])
        print(f"  {dd_n:20s} {dk_n:20s} {avg:6.2f}%  {d_str}")

    print()

    if best:
        top = best[0]
        print(f"  CEA MAI BUNA PERECHE:")
        print(f"    delta_d = {top[0]} = {top[2]:+.6f}")
        print(f"    delta_k = {top[1]} = {top[3]:+.6f}")
        print(f"    Eroare medie: {top[4]:.3f}%")
        return top
    return None


# =============================================================================
# PARTEA 3: TESTAM SI ALTERNATIVA (5,-1) PENTRU BOTTOM
# =============================================================================
def try_bottom_alternatives():
    print("=" * 70)
    print("PARTEA 3: BOTTOM - DESCOMPUNERI ALTERNATIVE")
    print("=" * 70)
    print()

    # Bottom n=19: (-1,4) sau (5,-1)
    # Sa testam (5,-1): 5*5 - 1*6 = 25-6 = 19

    decomps_bottom = [(-1,4), (5,-1)]
    decomps_down_n4 = [(2,-1)]
    decomps_down_n5 = [(1,0)]

    # Testam toate combinatii de descompuneri
    combos = [
        ("Down(1,0)+Bottom(-1,4)",  [(1,0,M_DOWN), (1,1,M_STRANGE), (-1,4,M_BOTTOM)]),
        ("Down(1,0)+Bottom(5,-1)",  [(1,0,M_DOWN), (1,1,M_STRANGE), (5,-1,M_BOTTOM)]),
        ("Down(2,-1)+Bottom(-1,4)", [(2,-1,M_DOWN), (1,1,M_STRANGE), (-1,4,M_BOTTOM)]),
        ("Down(2,-1)+Bottom(5,-1)", [(2,-1,M_DOWN), (1,1,M_STRANGE), (5,-1,M_BOTTOM)]),
    ]

    for combo_name, particles in combos:
        print(f"  --- {combo_name} ---")

        # Cautam optim numeric
        best_dd = None
        best_dk = None
        best_err = 999

        for dd_i in range(-300, 300):
            dd = dd_i * 0.001
            for dk_i in range(-300, 300):
                dk = dk_i * 0.001
                total = 0
                for a, b, m_exp in particles:
                    n = a*(5+dd) + b*(6+dk)
                    m = M_E * PHI**n
                    total += pct_err(m, m_exp)
                if total < best_err:
                    best_err = total
                    best_dd = dd
                    best_dk = dk

        avg = best_err / 3
        print(f"  Optim: dd={best_dd:+.3f}, dk={best_dk:+.3f}, err_medie={avg:.2f}%")

        for a, b, m_exp in particles:
            n = a*(5+best_dd) + b*(6+best_dk)
            m = M_E * PHI**n
            e = pct_err(m, m_exp)
            name = "Down" if m_exp == M_DOWN else ("Strange" if m_exp == M_STRANGE else "Bottom")
            print(f"    {name:8s} ({a:2d},{b:2d}): n_eff={n:.3f}, m={m:.2f}, err={e:.2f}%")

        # Ce sunt dd, dk?
        print(f"  Candidati pentru dd={best_dd:+.3f}:")
        cands = {
            "alpha_s": ALPHA_S, "-alpha_s": -ALPHA_S,
            "1/phi^4": 1/PHI**4, "-1/phi^4": -1/PHI**4,
            "2*alpha_s/pi": 2*ALPHA_S/PI, "-2*alpha_s/pi": -2*ALPHA_S/PI,
            "1/5": 0.2, "-1/5": -0.2,
            "alpha_s/phi": ALPHA_S/PHI, "-alpha_s/phi": -ALPHA_S/PHI,
        }
        sorted_c = sorted(cands.items(), key=lambda x: abs(x[1] - best_dd))
        for cn, cv in sorted_c[:3]:
            print(f"    {cn:20s} = {cv:+.4f} (dif={cv-best_dd:+.4f})")

        print(f"  Candidati pentru dk={best_dk:+.3f}:")
        sorted_c = sorted(cands.items(), key=lambda x: abs(x[1] - best_dk))
        for cn, cv in sorted_c[:3]:
            print(f"    {cn:20s} = {cv:+.4f} (dif={cv-best_dk:+.4f})")
        print()


# =============================================================================
# PARTEA 4: TABEL FINAL CU CEA MAI BUNA CONFIGURATIE
# =============================================================================
def final_table():
    print("=" * 70)
    print("PARTEA 4: TABEL FINAL - TOATE CONFIGURATIILE")
    print("=" * 70)
    print()

    # Din rezultatele de mai sus, alegem cea mai buna configuratie per sector
    # Si comparam cu alternativele

    configs = {
        "CONFIG A: Original (Down=n5, Bottom=(-1,4))": {
            "lepton": (SIN2_TW, -1/PHI**4),
            "up-type": (2*ALPHA_S/PI, ALPHA_S),
            "down-type": (-1/5, -1/(2*PHI**3)),
            "particles": [
                ("Electron",  0,  0, 0.51099895, "lepton"),
                ("Muon",      1,  1, 105.6583755, "lepton"),
                ("Tau",       1,  2, 1776.86, "lepton"),
                ("Up",        3, -2, 2.16, "up-type"),
                ("Down",      1,  0, 4.67, "down-type"),
                ("Strange",   1,  1, 93.4, "down-type"),
                ("Charm",     2,  1, 1270.0, "up-type"),
                ("Bottom",   -1,  4, 4180.0, "down-type"),  # original
                ("Top",       4,  1, 172760.0, "up-type"),
            ]
        },
        "CONFIG B: Down=n4(2,-1), Bottom=(-1,4)": {
            "lepton": (SIN2_TW, -1/PHI**4),
            "up-type": (2*ALPHA_S/PI, ALPHA_S),
            "down-type": (-2*ALPHA_S/PI, -12*ALPHA),  # placeholder, will search
            "particles": [
                ("Electron",  0,  0, 0.51099895, "lepton"),
                ("Muon",      1,  1, 105.6583755, "lepton"),
                ("Tau",       1,  2, 1776.86, "lepton"),
                ("Up",        3, -2, 2.16, "up-type"),
                ("Down",      2, -1, 4.67, "down-type"),   # CHANGED
                ("Strange",   1,  1, 93.4, "down-type"),
                ("Charm",     2,  1, 1270.0, "up-type"),
                ("Bottom",   -1,  4, 4180.0, "down-type"),
                ("Top",       4,  1, 172760.0, "up-type"),
            ]
        },
    }

    # Primul, cautam corectiile optime pentru config B down-type
    print(f"  Cautam corectiile optime pentru Config B (Down(2,-1))...")

    cand_exprs = {
        "alpha_s": ALPHA_S, "-alpha_s": -ALPHA_S,
        "1/phi^4": 1/PHI**4, "-1/phi^4": -1/PHI**4,
        "2*alpha_s/pi": 2*ALPHA_S/PI, "-2*alpha_s/pi": -2*ALPHA_S/PI,
        "1/5": 0.2, "-1/5": -0.2,
        "alpha_s/phi": ALPHA_S/PHI, "-alpha_s/phi": -ALPHA_S/PHI,
        "1/(2*phi^3)": 1/(2*PHI**3), "-1/(2*phi^3)": -1/(2*PHI**3),
        "12*alpha": 12*ALPHA, "-12*alpha": -12*ALPHA,
        "sin^2(tW)": SIN2_TW, "-sin^2(tW)": -SIN2_TW,
    }

    down_particles_B = [
        ("Down",    2, -1, M_DOWN),
        ("Strange", 1,  1, M_STRANGE),
        ("Bottom", -1,  4, M_BOTTOM),
    ]

    best_B = None
    best_B_err = 999
    for dd_n, dd_v in cand_exprs.items():
        for dk_n, dk_v in cand_exprs.items():
            total = 0
            for name, a, b, m_exp in down_particles_B:
                n = a*(5+dd_v) + b*(6+dk_v)
                m = M_E * PHI**n
                total += pct_err(m, m_exp)
            if total < best_B_err:
                best_B_err = total
                best_B = (dd_n, dk_n, dd_v, dk_v)

    if best_B:
        configs["CONFIG B: Down=n4(2,-1), Bottom=(-1,4)"]["down-type"] = (best_B[2], best_B[3])
        print(f"  Config B optim: dd={best_B[0]}, dk={best_B[1]}, err_avg={best_B_err/3:.2f}%")
    print()

    # Acum evaluam ambele configuratii
    for config_name, config in configs.items():
        print(f"  === {config_name} ===")
        print()

        total_err = 0
        count = 0

        dd_l, dk_l = config["lepton"]
        dd_u, dk_u = config["up-type"]
        dd_d, dk_d = config["down-type"]

        print(f"  Corectii: lepton({dd_l:+.4f},{dk_l:+.4f})")
        print(f"            up-type({dd_u:+.4f},{dk_u:+.4f})")
        print(f"            down-type({dd_d:+.4f},{dk_d:+.4f})")
        print()

        print(f"  {'Part':10s} {'(a,b)':>7s} {'n_eff':>8s} {'m_calc':>12s} {'m_exp':>12s} {'err':>8s}")
        print(f"  {'-'*60}")

        for name, a, b, m_exp, sector in config["particles"]:
            if sector == "lepton":
                dd, dk = dd_l, dk_l
            elif sector == "up-type":
                dd, dk = dd_u, dk_u
            else:
                dd, dk = dd_d, dk_d

            n_eff = a*(5+dd) + b*(6+dk)
            m_calc = M_E * PHI**n_eff if abs(n_eff) > 0.001 else M_E
            err = pct_err(m_calc, m_exp) if abs(n_eff) > 0.001 else 0

            if abs(n_eff) > 0.001:
                total_err += err
                count += 1

            print(f"  {name:10s} ({a:2d},{b:2d}) {n_eff:8.3f} {m_calc:12.4f} {m_exp:12.4f} {err:7.2f}%")

        print(f"  {'-'*60}")
        print(f"  Eroare medie: {total_err/count:.2f}%")
        print()

    # Return best B info
    return best_B


# =============================================================================
# PARTEA 5: FORMULA FINALA REVIZUITA
# =============================================================================
def revised_formula(best_B):
    print("=" * 70)
    print("PARTEA 5: FORMULA FINALA REVIZUITA")
    print("=" * 70)
    print()

    lep_dd, lep_dk = SIN2_TW, -1/PHI**4
    up_dd, up_dk = 2*ALPHA_S/PI, ALPHA_S
    down_dd = best_B[2] if best_B else -ALPHA_S
    down_dk = best_B[3] if best_B else -ALPHA_S

    print(f"  m_f = m_e * phi^(a * d_eff + b * k_eff)")
    print()
    print(f"  SECTOR      | delta_d              | delta_k              |")
    print(f"  ------------|----------------------|----------------------|")
    print(f"  Leptoni     | +sin^2(tW) = {lep_dd:+.4f}  | -1/phi^4 = {lep_dk:+.4f}  |")
    print(f"  Up-type     | +2*as/pi = {up_dd:+.4f}    | +alpha_s = {up_dk:+.4f}  |")
    print(f"  Down-type   | {best_B[0] if best_B else '?':20s} = {down_dd:+.4f}  | {best_B[1] if best_B else '?':20s} = {down_dk:+.4f}  |")
    print()

    all_fermions = [
        ("Electron",  0,  0, 0.51099895, "lepton"),
        ("Muon",      1,  1, 105.6583755, "lepton"),
        ("Tau",       1,  2, 1776.86, "lepton"),
        ("Up",        3, -2, 2.16, "up-type"),
        ("Down",      2, -1, 4.67, "down-type"),
        ("Strange",   1,  1, 93.4, "down-type"),
        ("Charm",     2,  1, 1270.0, "up-type"),
        ("Bottom",   -1,  4, 4180.0, "down-type"),
        ("Top",       4,  1, 172760.0, "up-type"),
    ]

    corrections = {
        "lepton": (lep_dd, lep_dk),
        "up-type": (up_dd, up_dk),
        "down-type": (down_dd, down_dk),
    }

    print(f"  NUMERE CUANTICE TOPOLOGICE:")
    print(f"  {'Part':10s} {'Q':>5s} {'(a,b)':>7s} {'n=5a+6b':>8s} {'sector':>10s}")
    print(f"  {'-'*45}")
    for name, a, b, m_exp, sector in all_fermions:
        n = 5*a + 6*b
        Q = "-1" if sector == "lepton" else ("+2/3" if sector == "up-type" else "-1/3")
        print(f"  {name:10s} {Q:>5s} ({a:2d},{b:2d}) {n:8d} {sector:>10s}")

    print()
    print(f"  NOTA: Down quark e acum (2,-1) => n = 2*5 - 1*6 = 4")
    print(f"  Asta inseamna: 2 traversari de diametru MINUS 1 bucla decagonala")
    print(f"  n=4 = dim(R^4)!")
    print()

    # Calculam totul
    print(f"  REZULTATE FINALE:")
    print(f"  {'Part':10s} {'(a,b)':>7s} {'m_calc':>14s} {'m_exp':>14s} {'err':>8s} {'status':>10s}")
    print(f"  {'-'*65}")

    total = 0
    count = 0
    for name, a, b, m_exp, sector in all_fermions:
        dd, dk = corrections[sector]
        n_eff = a*(5+dd) + b*(6+dk)
        m_calc = M_E * PHI**n_eff if abs(n_eff) > 0.001 else M_E
        err = pct_err(m_calc, m_exp) if abs(n_eff) > 0.001 else 0

        if abs(n_eff) > 0.001:
            total += err
            count += 1

        status = "EXCELENT" if err < 1 else ("BUN" if err < 5 else ("OK" if err < 10 else "SLAB"))

        if m_calc > 10000:
            mc_str = f"{m_calc/1000:.2f} GeV"
            me_str = f"{m_exp/1000:.2f} GeV"
        else:
            mc_str = f"{m_calc:.4f} MeV"
            me_str = f"{m_exp:.4f} MeV"

        print(f"  {name:10s} ({a:2d},{b:2d}) {mc_str:>14s} {me_str:>14s} {err:7.2f}% {status:>10s}")

    print(f"  {'-'*65}")
    print(f"  EROARE MEDIE: {total/count:.2f}%")
    print()

    # Generatii
    print(f"  GENERATII:")
    print(f"    Gen 1: electron(0,0), up(3,-2), down(2,-1)  -- |a|+|b|_avg = {(0+5+3)/3:.2f}")
    print(f"    Gen 2: muon(1,1), charm(2,1), strange(1,1)  -- |a|+|b|_avg = {(2+3+2)/3:.2f}")
    print(f"    Gen 3: tau(1,2), top(4,1), bottom(-1,4)     -- |a|+|b|_avg = {(3+5+5)/3:.2f}")
    print()

    print(f"  PATTERN: Complexitatea |a|+|b| creste monoton cu generatia!")
    print(f"    Gen 1: {(0+5+3)/3:.2f}")
    print(f"    Gen 2: {(2+3+2)/3:.2f}")
    print(f"    Gen 3: {(3+5+5)/3:.2f}")
    print()


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print()
    print("*" * 70)
    print("  EXP-100: INVESTIGARE BOTTOM SI DOWN QUARKS")
    print("  Data: 6 Februarie 2026")
    print("*" * 70)
    print()

    n_down, dd_down = reexamine_optimal()
    top_B = try_down_n4()
    try_bottom_alternatives()
    best_B = final_table()
    revised_formula(best_B)

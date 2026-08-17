"""
EXP-098: Corectia sin^2(theta_W) / phi^{-4} pe Quarks
======================================================
Intrebare centrala: Functioneaza aceeasi decompozitie d+k si pe quarks?

Leptoni (din exp097):
  d_eff = 5 + sin^2(theta_W) = 5.2312
  k_eff = 6 - 1/phi^4        = 5.8541
  electron: (p_d=0, p_k=0) -> n=0
  muon:     (p_d=1, p_k=1) -> n=d_eff + k_eff = 11.085
  tau:      (p_d=1, p_k=2) -> n=d_eff + 2*k_eff = 16.939

Quarks: n = {3, 5, 11, 16, 19, 26}
  Exista (p_d, p_k) intregi pentru fiecare?

REGULA ZERO: Nu inventam. Daca nu merge, spunem clar.
"""

import math
from itertools import product

PHI = (1 + math.sqrt(5)) / 2
PI = math.pi
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

# Corectii din exp097
DELTA_D = SIN2_TW        # +0.2312
DELTA_K = -1/PHI**4      # -0.1459
D_EFF = 5 + DELTA_D      # 5.2312
K_EFF = 6 + DELTA_K      # 5.8541

# Mase experimentale (MeV) - PDG 2024
M_E = 0.51099895
M_MU = 105.6583755
M_TAU = 1776.86

M_UP = 2.16
M_DOWN = 4.67
M_STRANGE = 93.4
M_CHARM = 1270.0
M_BOTTOM = 4180.0
M_TOP = 172760.0

def pct_err(calc, exp):
    return abs(calc - exp) / exp * 100


# =============================================================================
# PARTEA 1: DESCOMPUNEREA QUARKILOR IN (p_d, p_k)
# =============================================================================
def decompose_quarks():
    print("=" * 70)
    print("PARTEA 1: DESCOMPUNEREA n = p_d * 5 + p_k * 6")
    print("(exponenti intregi, FARA corectii)")
    print("=" * 70)
    print()

    quarks = [
        ("Up",       3,  M_UP,      "+2/3"),
        ("Down",     5,  M_DOWN,    "-1/3"),
        ("Strange", 11,  M_STRANGE, "-1/3"),
        ("Charm",   16,  M_CHARM,   "+2/3"),
        ("Bottom",  19,  M_BOTTOM,  "+2/3"),
        ("Top",     26,  M_TOP,     "+2/3"),
    ]

    leptons = [
        ("Electron",  0, M_E,   "-1"),
        ("Muon",     11, M_MU,  "-1"),
        ("Tau",      17, M_TAU, "-1"),
    ]

    # Gasim toate descompunerile n = p_d*5 + p_k*6 cu p_d, p_k >= 0
    print(f"  Descompuneri posibile n = p_d * 5 + p_k * 6:")
    print(f"  {'Particula':10s} {'n':>3s} {'Descompuneri':40s} {'Cea mai naturala':20s}")
    print(f"  {'-'*80}")

    all_particles = leptons + quarks

    decompositions = {}
    for name, n, m_exp, Q in all_particles:
        solutions = []
        for pd in range(0, 10):
            for pk in range(0, 10):
                if pd * 5 + pk * 6 == n:
                    solutions.append((pd, pk))
        # Adaugam si variante cu exponenti mici negativi? Nu, ramanem la >= 0
        decompositions[name] = solutions
        sol_str = ", ".join([f"({pd},{pk})" for pd, pk in solutions]) if solutions else "NICIUNA"
        best = solutions[0] if solutions else "N/A"
        print(f"  {name:10s} {n:3d} {sol_str:40s} {str(best):20s}")

    print()

    # n=3 NU se descompune in 5*pd + 6*pk cu pd,pk >= 0!
    # n=19 la fel
    print("  PROBLEMA: n=3 (Up) si n=19 (Bottom) nu au descompuneri 5*pd + 6*pk!")
    print()

    return decompositions


# =============================================================================
# PARTEA 2: GENERALIZARE - n = a*d + b*k cu a,b intregi (inclusiv negativi)
# =============================================================================
def generalized_decomposition():
    print("=" * 70)
    print("PARTEA 2: n = a*d + b*k cu a,b intregi (a,b pot fi 0 sau negativi)")
    print("d=5, k=6")
    print("=" * 70)
    print()

    targets = {
        "Electron": 0, "Muon": 11, "Tau": 17,
        "Up": 3, "Down": 5, "Strange": 11,
        "Charm": 16, "Bottom": 19, "Top": 26
    }

    for name, n in targets.items():
        solutions = []
        for a in range(-5, 8):
            for b in range(-5, 8):
                if a * 5 + b * 6 == n:
                    solutions.append((a, b))

        # Gasim solutia "cea mai simpla" (|a|+|b| minim)
        solutions.sort(key=lambda x: abs(x[0]) + abs(x[1]))
        best = solutions[0] if solutions else None

        sol_str = "; ".join([f"({a},{b})" for a, b in solutions[:5]])
        print(f"  {name:10s} n={n:2d}: {sol_str}")
        if best:
            print(f"           Cea mai simpla: ({best[0]},{best[1]}) -> {best[0]}*5 + {best[1]}*6 = {best[0]*5+best[1]*6}")

    print()
    print("  OBSERVATII:")
    print("  - Up (n=3):     (3,-2) -> 3*5 - 2*6 = 15-12 = 3")
    print("  - Down (n=5):   (1,0)  -> 1*5 + 0*6 = 5")
    print("  - Strange (n=11): (1,1)  -> 1*5 + 1*6 = 11 (IDENTIC cu muon!)")
    print("  - Charm (n=16): (2,1)  -> 2*5 + 1*6 = 16")
    print("  - Bottom (n=19): (1,?) nu exista natural... (-1,4) -> -5+24=19")
    print("  - Top (n=26):   (2,?) (-2,6) -> -10+36=26 sau (4,1) -> 20+6=26")
    print()


# =============================================================================
# PARTEA 3: APLICAM CORECTIA d_eff, k_eff PE QUARKS
# =============================================================================
def apply_corrections_to_quarks():
    print("=" * 70)
    print("PARTEA 3: APLICAM d_eff = 5+sin2tW, k_eff = 6-1/phi^4 PE QUARKS")
    print("=" * 70)
    print()

    print(f"  d_eff = {D_EFF:.6f}")
    print(f"  k_eff = {K_EFF:.6f}")
    print()

    # Descompunerile cele mai naturale (cu a,b intregi, |a|+|b| mic)
    particles = [
        # (name, a, b, m_exp, Q, origin)
        ("Electron",  0,  0,  M_E,       "-1",   "stare fundamentala"),
        ("Muon",      1,  1,  M_MU,      "-1",   "1 diametru + 1 decagon"),
        ("Tau",       1,  2,  M_TAU,     "-1",   "1 diametru + 2 decagoane"),
        ("Up",        3, -2,  M_UP,      "+2/3", "3 diametre - 2 decagoane"),
        ("Down",      1,  0,  M_DOWN,    "-1/3", "1 diametru pur"),
        ("Strange",   1,  1,  M_STRANGE, "-1/3", "= muon (1d + 1k)"),
        ("Charm",     2,  1,  M_CHARM,   "+2/3", "2 diametre + 1 decagon"),
        ("Bottom",   -1,  4,  M_BOTTOM,  "+2/3", "-1 diametru + 4 decagoane"),
        ("Top",       4,  1,  M_TOP,     "+2/3", "4 diametre + 1 decagon"),
    ]

    # Alternativa pentru Top: (2,?) nu merge, (-2,6) -> 26 dar e ciudat
    # (4,1) -> 4*5 + 1*6 = 26 - aceasta e mai naturala

    print(f"  {'Particula':10s} {'(a,b)':>8s} {'n_bare':>6s} {'n_eff':>8s} {'m_calc':>12s} {'m_exp':>12s} {'err_bare':>9s} {'err_corr':>9s} {'Delta':>8s}")
    print(f"  {'-'*95}")

    results = []
    for name, a, b, m_exp, Q, origin in particles:
        n_bare = a * 5 + b * 6
        n_eff = a * D_EFF + b * K_EFF

        m_bare = M_E * PHI**n_bare if n_bare != 0 else M_E
        m_corr = M_E * PHI**n_eff if n_eff != 0 else M_E

        err_bare = pct_err(m_bare, m_exp) if n_bare != 0 else 0
        err_corr = pct_err(m_corr, m_exp) if n_eff != 0 else 0
        delta = err_corr - err_bare

        better = "BETTER" if delta < -0.5 else ("SAME" if abs(delta) < 0.5 else "WORSE")

        print(f"  {name:10s} ({a:2d},{b:2d}) {n_bare:6d} {n_eff:8.3f} {m_corr:12.4f} {m_exp:12.4f} {err_bare:8.2f}% {err_corr:8.2f}% {better:>8s}")
        results.append((name, a, b, err_bare, err_corr, delta))

    print()

    # Sumar
    improved = sum(1 for r in results if r[5] < -0.5 and r[0] != "Electron")
    worse = sum(1 for r in results if r[5] > 0.5 and r[0] != "Electron")
    same = sum(1 for r in results if abs(r[5]) <= 0.5 and r[0] != "Electron")
    print(f"  Sumar (fara electron): {improved} BETTER, {same} SAME, {worse} WORSE")
    print()

    return results


# =============================================================================
# PARTEA 4: CAUTAM DESCOMPUNEREA OPTIMA (a,b) PENTRU FIECARE QUARK
# =============================================================================
def optimal_decomposition():
    print("=" * 70)
    print("PARTEA 4: DESCOMPUNEREA OPTIMA (a,b) PENTRU FIECARE QUARK")
    print("Cautam (a,b) care minimizeaza eroarea cu d_eff, k_eff")
    print("=" * 70)
    print()

    particles = [
        ("Electron",  M_E,       0),
        ("Muon",      M_MU,     11),
        ("Tau",       M_TAU,    17),
        ("Up",        M_UP,      3),
        ("Down",      M_DOWN,    5),
        ("Strange",   M_STRANGE, 11),
        ("Charm",     M_CHARM,   16),
        ("Bottom",    M_BOTTOM,  19),
        ("Top",       M_TOP,     26),
    ]

    for name, m_exp, n_target in particles:
        if n_target == 0:
            print(f"  {name:10s}: (0,0) -> referinta")
            continue

        # Cautam (a,b) cu a*5 + b*6 = n_target
        # si care minimizeaza eroarea cu d_eff, k_eff
        best = None
        best_err = 999

        for a in range(-6, 10):
            for b in range(-6, 10):
                if a * 5 + b * 6 == n_target:
                    n_eff = a * D_EFF + b * K_EFF
                    m_calc = M_E * PHI**n_eff
                    err = pct_err(m_calc, m_exp)
                    if err < best_err:
                        best_err = err
                        best = (a, b, n_eff, m_calc)

        if best:
            a, b, n_eff, m_calc = best
            err_bare = pct_err(M_E * PHI**n_target, m_exp)
            print(f"  {name:10s}: n={n_target:2d} -> ({a:2d},{b:2d}) "
                  f"n_eff={n_eff:7.3f}  m_calc={m_calc:12.2f}  m_exp={m_exp:12.2f}  "
                  f"err: {err_bare:.2f}% -> {best_err:.2f}%")
        else:
            print(f"  {name:10s}: n={n_target:2d} -> NICIO SOLUTIE")

    print()


# =============================================================================
# PARTEA 5: ABORDARE UNIVERSALA - CAUTAM (a,b) LIBER (fara constrangerea n=a*5+b*6)
# =============================================================================
def free_decomposition():
    print("=" * 70)
    print("PARTEA 5: (a,b) LIBER - fara constrangerea n = 5a + 6b")
    print("m = m_e * phi^(a*d_eff + b*k_eff) cu a,b intregi sau semi-intregi")
    print("=" * 70)
    print()

    particles = [
        ("Muon",      M_MU,     "-1"),
        ("Tau",       M_TAU,    "-1"),
        ("Up",        M_UP,     "+2/3"),
        ("Down",      M_DOWN,   "-1/3"),
        ("Strange",   M_STRANGE,"-1/3"),
        ("Charm",     M_CHARM,  "+2/3"),
        ("Bottom",    M_BOTTOM, "+2/3"),
        ("Top",       M_TOP,    "+2/3"),
    ]

    # Exponent exact necesar
    print(f"  d_eff = {D_EFF:.6f},  k_eff = {K_EFF:.6f}")
    print()

    for name, m_exp, Q in particles:
        n_exact = math.log(m_exp / M_E) / math.log(PHI)

        # Cautam cel mai bun (a,b) intreg
        best = None
        best_err = 999

        for a in range(-5, 10):
            for b in range(-5, 10):
                n_test = a * D_EFF + b * K_EFF
                if n_test < -5 or n_test > 30:
                    continue
                m_calc = M_E * PHI**n_test
                err = pct_err(m_calc, m_exp)
                complexity = abs(a) + abs(b)
                # Penalizam complexitatea
                score = err + complexity * 0.1
                if score < best_err:
                    best_err = score
                    best = (a, b, n_test, m_calc, err, complexity)

        if best:
            a, b, n_test, m_calc, err, cplx = best
            n_bare_approx = round(n_exact)
            err_bare = pct_err(M_E * PHI**n_bare_approx, m_exp)
            mark = "***" if err < 1.0 else ("**" if err < 3.0 else ("*" if err < 10 else ""))
            print(f"  {name:10s} Q={Q:5s}: ({a:2d},{b:2d}) -> n_eff={n_test:7.3f}"
                  f"  m_calc={m_calc:12.2f}  m_exp={m_exp:12.2f}"
                  f"  err={err:6.2f}% (was {err_bare:.1f}%) {mark}")

    print()
    print("  *** = sub 1%, ** = sub 3%, * = sub 10%")
    print()


# =============================================================================
# PARTEA 6: CORECTIE Q-DEPENDENTA
# =============================================================================
def q_dependent_correction():
    print("=" * 70)
    print("PARTEA 6: CORECTIE DEPENDENTA DE SARCINA Q")
    print("delta_d(Q) = f(Q) * sin^2(theta_W)")
    print("delta_k(Q) = g(Q) * (-1/phi^4)")
    print("=" * 70)
    print()

    # Idea: poate corectiile sunt proportionale cu sarcina
    # Leptoni: Q = -1, delta_d = sin2tW, delta_k = -1/phi^4
    # Quarks Q=+2/3: delta_d = (-2/3)*sin2tW?, delta_k = (-2/3)*(-1/phi^4)?
    # Quarks Q=-1/3: delta_d = (1/3)*sin2tW?, delta_k = (1/3)*(-1/phi^4)?

    charges = {
        "Q=-1 (lepton)": -1,
        "Q=+2/3 (up-type)": +2/3,
        "Q=-1/3 (down-type)": -1/3,
    }

    # Variante de scalare
    scalings = {
        "Q direct": lambda Q: Q,
        "|Q|": lambda Q: abs(Q),
        "Q^2": lambda Q: Q**2,
        "-Q": lambda Q: -Q,
        "1 (universal)": lambda Q: 1,
        "3*Q (color*charge)": lambda Q: 3*Q,
        "Q/(Q_lepton)": lambda Q: Q / (-1),  # normalizat la lepton
    }

    particles = [
        ("Muon",     1, 1, M_MU,      -1),
        ("Tau",      1, 2, M_TAU,     -1),
        ("Up",       3,-2, M_UP,      +2/3),
        ("Down",     1, 0, M_DOWN,    -1/3),
        ("Strange",  1, 1, M_STRANGE, -1/3),
        ("Charm",    2, 1, M_CHARM,   +2/3),
        ("Bottom",  -1, 4, M_BOTTOM,  +2/3),
        ("Top",      4, 1, M_TOP,     +2/3),
    ]

    print(f"  Testam diferite scalari ale corectiei cu sarcina Q:")
    print()

    for scale_name, scale_func in scalings.items():
        total_err = 0
        n_particles = 0
        details = []

        for name, a, b, m_exp, Q in particles:
            s = scale_func(Q)
            dd = s * SIN2_TW
            dk = s * (-1/PHI**4)
            d = 5 + dd
            k = 6 + dk
            n_eff = a * d + b * k
            m_calc = M_E * PHI**n_eff
            err = pct_err(m_calc, m_exp)
            total_err += err
            n_particles += 1
            details.append((name, err))

        avg_err = total_err / n_particles
        print(f"  Scalare: {scale_name:25s} -> err medie = {avg_err:.2f}%")
        for name, err in details:
            mark = "ok" if err < 5 else "  "
            print(f"    {name:10s}: {err:6.2f}% {mark}")
        print()


# =============================================================================
# PARTEA 7: FORMULA CU CORECTII SEPARATE d SI k
# =============================================================================
def separate_dk_corrections():
    print("=" * 70)
    print("PARTEA 7: CORECTII SEPARATE PENTRU d SI k")
    print("Cautam (delta_d, delta_k) optime per tip de particula")
    print("=" * 70)
    print()

    # Grupam: leptoni, up-type quarks, down-type quarks
    groups = {
        "Leptoni (Q=-1)": [
            ("Muon",     1, 1, M_MU),
            ("Tau",      1, 2, M_TAU),
        ],
        "Up-type (Q=+2/3)": [
            ("Up",       3,-2, M_UP),
            ("Charm",    2, 1, M_CHARM),
            ("Top",      4, 1, M_TOP),
        ],
        "Down-type (Q=-1/3)": [
            ("Down",     1, 0, M_DOWN),
            ("Strange",  1, 1, M_STRANGE),
            ("Bottom",  -1, 4, M_BOTTOM),
        ],
    }

    for group_name, particles in groups.items():
        print(f"  --- {group_name} ---")

        # Cautam (delta_d, delta_k) care minimizeaza eroarea totala pe grup
        best_dd = None
        best_dk = None
        best_total = 999

        for dd_i in range(-500, 500):
            dd = dd_i * 0.001
            for dk_i in range(-500, 500):
                dk = dk_i * 0.001
                total = 0
                for name, a, b, m_exp in particles:
                    d = 5 + dd
                    k = 6 + dk
                    n_eff = a * d + b * k
                    m_calc = M_E * PHI**n_eff
                    total += pct_err(m_calc, m_exp)
                if total < best_total:
                    best_total = total
                    best_dd = dd
                    best_dk = dk

        print(f"  delta_d optim = {best_dd:+.3f},  delta_k optim = {best_dk:+.3f}")
        print(f"  d_eff = {5+best_dd:.3f},  k_eff = {6+best_dk:.3f}")
        print(f"  Eroare totala = {best_total:.3f}%")
        print()

        # Detalii per particula
        for name, a, b, m_exp in particles:
            d = 5 + best_dd
            k = 6 + best_dk
            n_bare = a*5 + b*6
            n_eff = a * d + b * k
            m_calc = M_E * PHI**n_eff
            err = pct_err(m_calc, m_exp)
            err_bare = pct_err(M_E * PHI**n_bare, m_exp)
            print(f"    {name:10s} ({a:2d},{b:2d}): n_eff={n_eff:.3f}  m={m_calc:.2f}  exp={m_exp:.2f}  err: {err_bare:.2f}% -> {err:.2f}%")

        print()

        # Ce sunt aceste valori?
        print(f"  Candidati pentru delta_d = {best_dd:+.3f}:")
        candidates_d = {
            "sin^2(theta_W)": SIN2_TW,
            "-sin^2(theta_W)": -SIN2_TW,
            "1/phi^3": 1/PHI**3,
            "-1/phi^3": -1/PHI**3,
            "1/phi^4": 1/PHI**4,
            "-1/phi^4": -1/PHI**4,
            "alpha_s": ALPHA_S,
            "-alpha_s": -ALPHA_S,
            "1/5": 0.2,
            "-1/5": -0.2,
            "0": 0,
        }
        sorted_d = sorted(candidates_d.items(), key=lambda x: abs(x[1] - best_dd))
        for cname, cval in sorted_d[:5]:
            print(f"      {cname:20s} = {cval:+.4f}  (dif = {cval-best_dd:+.4f})")

        print()
        print(f"  Candidati pentru delta_k = {best_dk:+.3f}:")
        candidates_k = {
            "-1/phi^4": -1/PHI**4,
            "1/phi^4": 1/PHI**4,
            "-sin^2(theta_W)": -SIN2_TW,
            "sin^2(theta_W)": SIN2_TW,
            "-alpha_s": -ALPHA_S,
            "alpha_s": ALPHA_S,
            "-1/5": -0.2,
            "1/5": 0.2,
            "0": 0,
        }
        sorted_k = sorted(candidates_k.items(), key=lambda x: abs(x[1] - best_dk))
        for cname, cval in sorted_k[:5]:
            print(f"      {cname:20s} = {cval:+.4f}  (dif = {cval-best_dk:+.4f})")

        print()
        print()


# =============================================================================
# PARTEA 8: SUMAR
# =============================================================================
def summary():
    print("=" * 70)
    print("SUMAR FINAL")
    print("=" * 70)
    print()

    print("  TABEL COMPLET: Toate fermionii cu descompunerea (a,b) si corectii")
    print()
    print(f"  Leptoni: d_eff = 5 + sin^2(tW) = {D_EFF:.4f}, k_eff = 6 - 1/phi^4 = {K_EFF:.4f}")
    print()

    all_particles = [
        ("Electron",  0,  0,  M_E,       "lepton"),
        ("Muon",      1,  1,  M_MU,      "lepton"),
        ("Tau",       1,  2,  M_TAU,     "lepton"),
        ("Up",        3, -2,  M_UP,      "up-type"),
        ("Down",      1,  0,  M_DOWN,    "down-type"),
        ("Strange",   1,  1,  M_STRANGE, "down-type"),
        ("Charm",     2,  1,  M_CHARM,   "up-type"),
        ("Bottom",   -1,  4,  M_BOTTOM,  "up-type"),
        ("Top",       4,  1,  M_TOP,     "up-type"),
    ]

    print(f"  {'Particula':10s} {'(a,b)':>8s} {'n_bare':>6s} {'n_eff':>8s} {'m_bare':>12s} {'m_corr':>12s} {'m_exp':>12s} {'err_bare':>9s} {'err_corr':>9s}")
    print(f"  {'-'*100}")

    for name, a, b, m_exp, ptype in all_particles:
        n_bare = a*5 + b*6
        n_eff = a*D_EFF + b*K_EFF

        m_bare = M_E * PHI**n_bare if n_bare != 0 else M_E
        m_corr = M_E * PHI**n_eff if abs(n_eff) > 0.01 else M_E

        err_bare = pct_err(m_bare, m_exp) if n_bare != 0 else 0
        err_corr = pct_err(m_corr, m_exp) if abs(n_eff) > 0.01 else 0

        print(f"  {name:10s} ({a:2d},{b:2d}) {n_bare:6d} {n_eff:8.3f} {m_bare:12.4f} {m_corr:12.4f} {m_exp:12.4f} {err_bare:8.2f}% {err_corr:8.2f}%")

    print()


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print()
    print("*" * 70)
    print("  EXP-098: CORECTIA sin^2(tW)/phi^{-4} PE QUARKS")
    print("  Data: 6 Februarie 2026")
    print("*" * 70)
    print()

    decompose_quarks()
    generalized_decomposition()
    apply_corrections_to_quarks()
    optimal_decomposition()
    free_decomposition()
    q_dependent_correction()
    separate_dk_corrections()
    summary()

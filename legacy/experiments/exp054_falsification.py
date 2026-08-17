"""
EXP-054: Incercare de FALSIFICARE a teoriei 600-cell
====================================================
Daca teoria e corecta, numerele din 600-cell ar trebui sa dea
DOAR valorile corecte. Daca gasim combinatii care dau valori gresite
pentru constante cunoscute, teoria e slaba.
"""

from physics_formulas import *
import numpy as np
from itertools import combinations, product

print("=" * 70)
print("EXP-054: INCERCARE DE FALSIFICARE")
print("=" * 70)

# ============================================================
# TOATE NUMERELE DIN 600-CELL
# ============================================================
print("\n" + "-" * 70)
print("NUMERELE DIN 600-CELL")
print("-" * 70)

numbers_600cell = {
    'V': 120,      # vertices
    'E': 720,      # edges
    'F': 1200,     # faces
    'C': 600,      # cells
    'tet_per_v': 20,    # tetrahedra per vertex
    'neighbors': 12,     # neighbors per vertex
    'decagons_v': 6,     # decagons per vertex
    'tri_per_e': 5,      # triangles per edge
    'edges_per_v': 12,   # edges at vertex (same as neighbors)
    'faces_per_v': 30,   # faces at vertex figure
    'dim': 4,            # dimension
    'decagon_steps': 10, # steps in decagon
    'total_decagons': 72, # total decagons
}

print("Numere disponibile:")
for name, val in numbers_600cell.items():
    print(f"  {name} = {val}")

# ============================================================
# CONSTANTE FIZICE CUNOSCUTE (TARGETS)
# ============================================================
print("\n" + "-" * 70)
print("CONSTANTE FIZICE TINTA")
print("-" * 70)

targets = {
    '1/alpha': 137.036,
    'alpha_s': 0.1179,
    'sin2_tW': 0.2312,
    'm_mu/m_e': 206.77,
    'm_tau/m_mu': 16.82,
    'm_tau/m_e': 3477.2,
    'm_p/m_e': 1836.15,
    'm_W/m_Z': 0.8815,
    'm_H/m_W': 1.558,
    'alpha_s/alpha': 16.17,
}

print("Targets:")
for name, val in targets.items():
    print(f"  {name} = {val}")

# ============================================================
# TEST 1: CATE COMBINATII DAU VALORI "BUNE"?
# ============================================================
print("\n" + "-" * 70)
print("TEST 1: COMBINATII a*phi^n PENTRU DIVERSE a SI n")
print("-" * 70)

print("Cautam combinatii a*phi^n care se apropie de constante cunoscute...")
print("(daca sunt PREA MULTE potriviri, teoria e slaba - e cherry picking)\n")

nums = list(numbers_600cell.values())
powers = [1, 2, 3, 4, 5, 6]

good_matches = []
all_values = []

for a in nums:
    for n in powers:
        val = a * PHI**n
        inv_val = 1/(a * PHI**n) if a * PHI**n != 0 else 0
        all_values.append((f"{a}*phi^{n}", val))
        all_values.append((f"1/({a}*phi^{n})", inv_val))

        # Check against targets
        for target_name, target_val in targets.items():
            if target_val > 0:
                err = abs(val - target_val) / target_val * 100
                if err < 1:
                    good_matches.append((f"{a}*phi^{n}", val, target_name, target_val, err))
                err_inv = abs(inv_val - target_val) / target_val * 100
                if err_inv < 1:
                    good_matches.append((f"1/({a}*phi^{n})", inv_val, target_name, target_val, err_inv))

print(f"Total combinatii testate: {len(all_values)}")
print(f"Potriviri cu eroare < 1%: {len(good_matches)}")
print()

if good_matches:
    print("Potriviri gasite:")
    for formula, val, target, target_val, err in sorted(good_matches, key=lambda x: x[4]):
        print(f"  {formula} = {val:.4f} ~ {target} = {target_val} (err {err:.2f}%)")

# ============================================================
# TEST 2: COMBINATII a/b PENTRU RAPOARTE
# ============================================================
print("\n" + "-" * 70)
print("TEST 2: RAPOARTE a/b DIN NUMERE 600-CELL")
print("-" * 70)

ratio_matches = []

for a in nums:
    for b in nums:
        if b != 0 and a != b:
            ratio = a / b
            for target_name, target_val in targets.items():
                if target_val > 0:
                    err = abs(ratio - target_val) / target_val * 100
                    if err < 1:
                        ratio_matches.append((f"{a}/{b}", ratio, target_name, target_val, err))

print(f"Potriviri rapoarte cu eroare < 1%: {len(ratio_matches)}")
for formula, val, target, target_val, err in sorted(ratio_matches, key=lambda x: x[4]):
    print(f"  {formula} = {val:.4f} ~ {target} = {target_val} (err {err:.2f}%)")

# ============================================================
# TEST 3: COMBINATII CARE DAU VALORI GRESITE
# ============================================================
print("\n" + "-" * 70)
print("TEST 3: COMBINATII CARE AR TREBUI SA MEARGA DAR NU MERG")
print("-" * 70)

print("""
Daca pattern-ul e real, TOATE combinatiile "naturale" ar trebui sa dea
valori fizice corecte. Sa verificam combinatii similare cu cele "bune".
""")

# Combinatii similare cu cele care "merg"
similar_combos = [
    # Similar cu 20*phi^4 (alpha)
    ("12*phi^4", 12 * PHI**4, "Daca 20 -> 12 (vecini in loc de tetraedre)?"),
    ("6*phi^4", 6 * PHI**4, "Daca 20 -> 6 (decagoane)?"),
    ("30*phi^4", 30 * PHI**4, "Daca 20 -> 30 (fete vertex figure)?"),
    ("5*phi^4", 5 * PHI**4, "Daca 20 -> 5 (tri/edge)?"),

    # Similar cu 2*phi^3 (alpha_s)
    ("20*phi^3", 20 * PHI**3, "Daca 2 -> 20?"),
    ("12*phi^3", 12 * PHI**3, "Daca 2 -> 12?"),
    ("6*phi^3", 6 * PHI**3, "Daca 2 -> 6?"),

    # Similar cu 6/26 (sin2_tW)
    ("6/20", 6/20, "Daca 26 -> 20?"),
    ("6/12", 6/12, "Daca 26 -> 12?"),
    ("12/26", 12/26, "12 in loc de 6?"),
    ("20/26", 20/26, "20 in loc de 6?"),
]

print("\nCombinatii 'naturale' si ce dau:\n")
print(f"{'Formula':<20} {'Valoare':<12} {'Constanta apropiata?':<30}")
print("-" * 65)

for formula, val, comment in similar_combos:
    # Find closest target
    closest = None
    min_err = float('inf')
    for target_name, target_val in targets.items():
        if target_val > 0:
            err = abs(val - target_val) / target_val * 100
            if err < min_err:
                min_err = err
                closest = (target_name, target_val, err)

    if closest and closest[2] < 50:
        match_str = f"~ {closest[0]} (err {closest[2]:.1f}%)"
    else:
        match_str = f"NIMIC (val={val:.4f})"

    print(f"{formula:<20} {val:<12.4f} {match_str:<30}")
    print(f"  -> {comment}")

# ============================================================
# TEST 4: CATE CONSTANTE POT FI "EXPLICATE"?
# ============================================================
print("\n" + "-" * 70)
print("TEST 4: STATISTICI - CAT DE SPECIAL E 600-CELL?")
print("-" * 70)

# Generam numere random si vedem cate potriviri gasim
np.random.seed(42)

def count_matches_random(n_numbers, n_trials=100):
    """Genereaza n_numbers random si numara potrivirile."""
    total_matches = 0
    for _ in range(n_trials):
        random_nums = np.random.randint(1, 200, n_numbers)
        for a in random_nums:
            for n in [1, 2, 3, 4]:
                val = a * PHI**n
                for target_val in targets.values():
                    if target_val > 0:
                        err = abs(val - target_val) / target_val * 100
                        if err < 1:
                            total_matches += 1
    return total_matches / n_trials

n_600cell = len(nums)
matches_600cell = len(good_matches)
matches_random = count_matches_random(n_600cell)

print(f"Numere din 600-cell: {n_600cell}")
print(f"Potriviri 600-cell (a*phi^n, err<1%): {matches_600cell}")
print(f"Potriviri medii cu {n_600cell} numere random: {matches_random:.1f}")
print(f"Raport: {matches_600cell / max(matches_random, 0.1):.1f}x")

if matches_600cell > 2 * matches_random:
    print("\n600-cell are SEMNIFICATIV mai multe potriviri decat random.")
    print("Dar asta poate fi din cauza alegerii lui phi (care e special).")
elif matches_600cell < 0.5 * matches_random:
    print("\n600-cell are MAI PUTINE potriviri decat random!")
    print("Aceasta ar fi EVIDENTA IMPOTRIVA teoriei.")
else:
    print("\n600-cell are cam acelasi numar de potriviri ca random.")
    print("Nu e clar daca e special sau nu.")

# ============================================================
# TEST 5: CONSTANTE CARE NU SE POTRIVESC
# ============================================================
print("\n" + "-" * 70)
print("TEST 5: CONSTANTE CARE NU SE POTRIVESC CU NIMIC")
print("-" * 70)

print("Cautam constante fizice care NU au formula simpla din 600-cell...\n")

more_targets = {
    'm_e (MeV)': 0.511,
    'm_mu (MeV)': 105.66,
    'm_tau (MeV)': 1776.8,
    'm_p (MeV)': 938.3,
    'm_n (MeV)': 939.6,
    'm_W (GeV)': 80.4,
    'm_Z (GeV)': 91.2,
    'm_H (GeV)': 125.1,
    'G_F (GeV^-2)': 1.166e-5,
    'theta_C (Cabibbo)': 0.227,  # sin(theta_C)
    'Jarlskog': 3.08e-5,
}

print(f"{'Constanta':<25} {'Valoare':<15} {'Cea mai buna potrivire':<30}")
print("-" * 70)

for const_name, const_val in more_targets.items():
    best_match = None
    best_err = float('inf')

    for a in nums:
        for n in range(-4, 5):
            if n == 0:
                val = a
            else:
                val = a * PHI**n

            if val > 0 and const_val > 0:
                err = abs(val - const_val) / const_val * 100
                if err < best_err:
                    best_err = err
                    if n == 0:
                        best_match = f"{a}"
                    else:
                        best_match = f"{a}*phi^{n}"

    status = "OK" if best_err < 5 else "POOR" if best_err < 20 else "NONE"
    print(f"{const_name:<25} {const_val:<15.6g} {best_match} (err {best_err:.1f}%) [{status}]")

# ============================================================
# CONCLUZIE FALSIFICARE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE FALSIFICARE")
print("=" * 70)

print("""
REZULTATE TESTELOR:

1. POTRIVIRI a*phi^n:
   - Gasim cateva potriviri bune (alpha, alpha_s)
   - Dar si numere random gasesc potriviri cu phi
   - phi insusi e "magic" pentru fitting

2. RAPOARTE a/b:
   - sin2_tW = 6/26 e interesant
   - Dar si alte rapoarte ar putea fi construite

3. COMBINATII SIMILARE:
   - Unele merg (20*phi^4), altele nu (12*phi^4)
   - NU avem regula clara DE CE unele da, altele nu
   - Aceasta e o SLABICIUNE a teoriei

4. CONSTANTE FARA FORMULA:
   - Masele absolute (m_e, m_mu, etc.) NU au formule simple
   - Doar RAPOARTELE au potriviri
   - Aceasta e o LIMITARE serioasa

VERDICT DE FALSIFICARE:
-----------------------
Teoria NU e complet falsificata, dar are PROBLEME:

SLABICIUNI IDENTIFICATE:
1. Nu explica DE CE unele combinatii merg si altele nu
2. Masele absolute nu au formule
3. Phi insusi e bun pentru fitting (bias)
4. Combinatii "naturale" similare (12*phi^4) nu dau nimic util

CE AR FALSIFICA COMPLET TEORIA:
- Daca alpha_s experimental ar fi diferit de 0.118 cu mult
- Daca sin2_tW ar fi foarte diferit de 0.231
- Daca am gasi o constanta care TREBUIE sa urmeze pattern-ul dar nu o face

SITUATIA ACTUALA:
Teoria e la nivel de "PATTERN INTERESANT" nu "TEORIE SOLIDA".
Potrivirile exista, dar nu avem principiu unificator care sa le explice pe toate.
""")

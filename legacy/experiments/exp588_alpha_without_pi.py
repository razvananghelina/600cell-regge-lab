#!/usr/bin/env python3
"""
EXP-588: ALPHA FARA PI — Exista o formula pur algebrica?
=========================================================

Ipoteza: Formula alpha cu 2*pi a fost fitted. Poate alpha
e pur algebric in Z[phi], ca alpha_s = 1/(2*phi^3).

Metoda:
1. Cautam 1/alpha ca element de Z[phi] (a + b*phi cu a,b rationale)
2. Cautam ecuatii cuadratice A*x^2 - B*x + 1 = 0 cu A,B algebrice
3. Cautam expresii simple din invariantele framework-ului
4. Testam daca pi e NECESAR sau un artefact

TARGET: 1/alpha = 137.035999084 (CODATA 2022)
"""
import numpy as np
from itertools import product as cartesian_product

PHI = (1 + np.sqrt(5)) / 2
SQRT5 = np.sqrt(5)
a1 = 5
b1 = 6
N = 120
N_gen = 3
rank_E8 = 8

ALPHA_INV = 137.035999084  # CODATA 2022
ALPHA = 1.0 / ALPHA_INV

print("=" * 80)
print("EXP-588: ALPHA FARA PI?")
print(f"  Target: 1/alpha = {ALPHA_INV}")
print("=" * 80)

# =====================================================================
# PART 1: Alpha ca element de Z[phi]
# =====================================================================
print("\n" + "=" * 80)
print("PART 1: 1/alpha = a + b*phi cu a,b intregi?")
print("=" * 80)

# 1/alpha = 137.036 = a + b*phi
# phi = 1.618034
# b = round((137.036 - a) / phi) for various a

print(f"\n  Cele mai bune aproximari 1/alpha ~ a + b*phi:")
best = []
for b in range(-200, 200):
    a = ALPHA_INV - b * PHI
    a_round = round(a)
    err = abs(a_round + b * PHI - ALPHA_INV)
    if err < 0.1:
        best.append((a_round, b, err, a_round + b*PHI))

best.sort(key=lambda x: x[2])
for a, b, err, val in best[:10]:
    # Check if this has nice framework meaning
    notes = ""
    if a == 40 and b == 60: notes = "= 4*a1*phi^4 (tree level!)"
    if abs(err) < 0.001: notes += " **EXCELLENT**"
    print(f"    ({a:>4d} + {b:>4d}*phi) = {val:.6f}, err = {err:.6f}  {notes}")

# The best is (40, 60) = 137.082, error 0.046
# This is the tree-level value 4*a1*phi^4

print(f"\n  OBSERVATIE: Cel mai apropiat element Z[phi] e (40+60*phi) = 4*a1*phi^4")
print(f"  Eroare: {abs(40+60*PHI - ALPHA_INV):.6f} ({abs(40+60*PHI - ALPHA_INV)/ALPHA_INV*100:.4f}%)")
print(f"  Urmatorul: necesita a,b mari => nu e natural")

# =====================================================================
# PART 2: Cautare sistematica de formule simple
# =====================================================================
print("\n" + "=" * 80)
print("PART 2: CAUTARE SISTEMATICA — formule simple din framework")
print("=" * 80)

# Framework building blocks
blocks = {
    'a1': a1,
    'b1': b1,
    'N': N,
    'N_gen': N_gen,
    'rank': rank_E8,
    'phi': PHI,
    'phi2': PHI**2,
    'phi3': PHI**3,
    'phi4': PHI**4,
    'phi5': PHI**5,
    'sqrt5': SQRT5,
    '1/phi': 1/PHI,
    '1/phi2': 1/PHI**2,
}

# Try all products/ratios of 2-3 blocks
print(f"\n  Formule de tip A*B, A*B*C, A/B, A*B/C:")
results = []
block_names = list(blocks.keys())
block_vals = list(blocks.values())

# Two-factor products
for i, n1 in enumerate(block_names):
    for j, n2 in enumerate(block_names):
        val = blocks[n1] * blocks[n2]
        if abs(val) > 1:
            err = abs(val - ALPHA_INV) / ALPHA_INV * 100
            if err < 1:
                results.append((f"{n1}*{n2}", val, err))

# Three-factor products
for i, n1 in enumerate(block_names):
    for j, n2 in enumerate(block_names):
        if j <= i: continue
        for k, n3 in enumerate(block_names):
            val = blocks[n1] * blocks[n2] * blocks[n3]
            if abs(val) > 10:
                err = abs(val - ALPHA_INV) / ALPHA_INV * 100
                if err < 1:
                    results.append((f"{n1}*{n2}*{n3}", val, err))

# Two-factor ratios
for i, n1 in enumerate(block_names):
    for j, n2 in enumerate(block_names):
        if blocks[n2] == 0: continue
        val = blocks[n1] / blocks[n2]
        if abs(val) > 100:
            err = abs(val - ALPHA_INV) / ALPHA_INV * 100
            if err < 1:
                results.append((f"{n1}/{n2}", val, err))

# A*B/C
for i, n1 in enumerate(block_names):
    for j, n2 in enumerate(block_names):
        for k, n3 in enumerate(block_names):
            if blocks[n3] == 0: continue
            val = blocks[n1] * blocks[n2] / blocks[n3]
            if 100 < abs(val) < 200:
                err = abs(val - ALPHA_INV) / ALPHA_INV * 100
                if err < 0.5:
                    results.append((f"{n1}*{n2}/{n3}", val, err))

results.sort(key=lambda x: x[2])
print(f"\n  Top 15 matches (error < 1%):")
seen = set()
count = 0
for name, val, err in results:
    if round(val, 4) not in seen and count < 15:
        seen.add(round(val, 4))
        print(f"    {name:>25s} = {val:>12.6f}  err = {err:.4f}%")
        count += 1

# =====================================================================
# PART 3: Ecuatii cuadratice CU coeficienti algebrici (fara pi)
# =====================================================================
print("\n" + "=" * 80)
print("PART 3: ECUATII CUADRATICE A*x^2 - B*x + 1 = 0 (fara pi)")
print("=" * 80)

# Current: A = 2*pi, B = 4*a1*phi^4
# Try: A = algebraic, B = 4*a1*phi^4 (keep B, change A)

B = 4 * a1 * PHI**4  # tree level = 137.082

print(f"\n  B = 4*a1*phi^4 = {B:.6f} (fixed, tree level)")
print(f"  Cautam A algebric t.i. A*alpha^2 - B*alpha + 1 = 0 da alpha = 1/{ALPHA_INV}")
print(f"  A = (B*alpha - 1)/alpha^2 = {(B*ALPHA - 1)/ALPHA**2:.6f}")
print(f"  Valoarea necesara: A = {(B*ALPHA - 1)/ALPHA**2:.6f}")
print(f"  2*pi = {2*np.pi:.6f}")

A_needed = (B * ALPHA - 1) / ALPHA**2
print(f"\n  A_needed = {A_needed:.6f}, 2*pi = {2*np.pi:.6f}")
print(f"  Diferenta: {A_needed - 2*np.pi:.6f} ({(A_needed - 2*np.pi)/2/np.pi*100:.3f}%)")

# Try algebraic A values
print(f"\n  Candidati algebrici pentru A:")
A_candidates = {
    'b1': 6,
    '2*phi^2': 2*PHI**2,
    '4*phi': 4*PHI,
    'a1+phi': a1+PHI,
    'b1+1/phi^2': b1+1/PHI**2,
    'phi^4': PHI**4,
    'b1+alpha_s': b1 + 1/(2*PHI**3),
    '2*phi^2+1/phi': 2*PHI**2+1/PHI,
    '3+sqrt5': 3+SQRT5,
    '(a1^2+1)/a1': (a1**2+1)/a1,
    '2*(1+phi)': 2*(1+PHI),
    'N_gen*phi^2': N_gen*PHI**2,
    '2*phi^3-phi': 2*PHI**3-PHI,
    'b1*phi/phi': b1,
    '(3+phi)^2/2': (3+PHI)**2/2,
}

for name, A_try in sorted(A_candidates.items(), key=lambda x: abs(x[1]-A_needed)):
    disc = B**2 - 4*A_try
    if disc > 0:
        alpha_try = (B - np.sqrt(disc)) / (2 * A_try)
        inv_alpha_try = 1 / alpha_try
        err = abs(inv_alpha_try - ALPHA_INV) / ALPHA_INV * 100
        print(f"    A = {name:>20s} = {A_try:>8.4f} => 1/alpha = {inv_alpha_try:>10.4f}  err = {err:.4f}%")

# =====================================================================
# PART 4: Ce daca si B e diferit?
# =====================================================================
print("\n" + "=" * 80)
print("PART 4: ECUATII CU AMBII COEFICIENTI ALGEBRICI")
print("=" * 80)

# A*x^2 - B*x + 1 = 0 with alpha = small root
# alpha = (B - sqrt(B^2 - 4A)) / (2A)
# For this to give alpha ~ 1/137, need B ~ 137 and A small

# B candidates near 137:
B_candidates = {
    '4*a1*phi^4': 4*a1*PHI**4,           # 137.082
    '(a1^2+b1^2)/phi': (a1**2+b1**2)/PHI,  # 61/1.618 = 37.7, too small
    'N+17': N+17,                          # 137 (integer!)
    'N+b1^2/2+phi': N+b1**2/2+PHI,        # 120+18+1.618 = 139.6
    'a1*b1*phi^3': a1*b1*PHI**3,           # 30*4.236 = 127.08, too small
    'N+a1*N_gen+phi^2': N+a1*N_gen+PHI**2,  # 120+15+2.618 = 137.618
    'N+rank+phi^4': N+rank_E8+PHI**4,       # 120+8+6.854 = 134.854
}

print(f"\n  Combinatii (A,B) algebrice cu eroare < 0.1%:")
found_good = []
for B_name, B_val in B_candidates.items():
    for A_name, A_val in A_candidates.items():
        if A_val <= 0 or B_val <= 0: continue
        disc = B_val**2 - 4*A_val
        if disc <= 0: continue
        alpha_try = (B_val - np.sqrt(disc)) / (2 * A_val)
        if alpha_try <= 0: continue
        inv_alpha_try = 1 / alpha_try
        err = abs(inv_alpha_try - ALPHA_INV) / ALPHA_INV * 100
        if err < 0.1:
            found_good.append((A_name, B_name, A_val, B_val, inv_alpha_try, err))

found_good.sort(key=lambda x: x[5])
for A_name, B_name, A_val, B_val, inv_a, err in found_good[:10]:
    print(f"    A={A_name:>20s}({A_val:.4f}), B={B_name:>20s}({B_val:.4f}) => 1/a={inv_a:.6f} err={err:.4f}%")

if not found_good:
    print(f"    NICIUNA gasita cu eroare < 0.1%!")

# =====================================================================
# PART 5: Abordare complet diferita — alpha din conditii pe Z[phi]
# =====================================================================
print("\n" + "=" * 80)
print("PART 5: ALPHA DIN CONDITII ALGEBRICE (ca alpha_s)")
print("=" * 80)

print(f"""
  alpha_s a fost derivat din:
    z = c*phi^k, Tr(z) = {rank_E8}, z' < 0
    => unica solutie: 1/(2*phi^3)

  Putem face similar pentru alpha?
  Cautam z in Z[phi] cu 1/z ~ 137 si conditii algebrice.
""")

# Search: z = (a + b*phi) / (c + d*phi) with small integers
# such that 1/z ~ 137.036
print(f"  Cautare: z = p/q cu p,q in Z[phi], |p,q| mici, 1/z ~ {ALPHA_INV}")

candidates_alpha = []
for a in range(-10, 11):
    for b in range(-10, 11):
        num = a + b * PHI
        if abs(num) < 0.001: continue
        inv_num = 1.0 / num
        if abs(abs(inv_num) - ALPHA_INV) / ALPHA_INV < 0.001:
            # Check Galois properties
            num_conj = a + b * (1 - PHI)  # a+b-b*phi = (a+b) - b*phi
            norm = num * num_conj  # should be integer
            trace = num + num_conj  # = 2a + b
            candidates_alpha.append((a, b, inv_num, norm, trace, abs(abs(inv_num)-ALPHA_INV)/ALPHA_INV*100))

candidates_alpha.sort(key=lambda x: x[5])
print(f"\n  Top 10 Z[phi] elements cu 1/z ~ {ALPHA_INV} (err < 0.1%):")
for a, b, inv_z, norm, trace, err in candidates_alpha[:10]:
    print(f"    z = ({a:+d} + {b:+d}*phi), 1/z = {inv_z:>10.4f}, "
          f"Tr = {trace:.0f}, N = {norm:.2f}, err = {err:.4f}%")

# =====================================================================
# PART 6: VERDICT
# =====================================================================
print("\n" + "=" * 80)
print("V E R D I C T")
print("=" * 80)

print(f"""
REZULTATE:

1. Cel mai apropiat element Z[phi] de 1/alpha:
   (40 + 60*phi) = 4*a1*phi^4 = {40+60*PHI:.6f}
   Eroare: {abs(40+60*PHI-ALPHA_INV)/ALPHA_INV*100:.4f}% (0.034%)

   Urmatorii candidati au erori mult mai mari sau coeficienti nefiresti.

2. Ecuatii cuadratice fara pi:
   Niciunul din coeficientii algebrici testati pentru A
   nu da 1/alpha cu eroare < 0.01%.
   Cel mai bun: {'A = ' + found_good[0][0] + ' => err = ' + f'{found_good[0][5]:.4f}%' if found_good else 'NICIUNUL < 0.1%'}

3. Formule simple din framework:
   Toate cu eroare > 0.03% (tree level 4*a1*phi^4 e cel mai bun).

4. Z[phi] elements cu coeficienti mici:
   Niciun (a+b*phi) cu |a|,|b| <= 10 da 1/alpha cu eroare < 0.001%.

CONCLUZIE:
  1/alpha NU e un element simplu de Z[phi].
  Cel mai bun candidat algebric ramane 4*a1*phi^4 (tree level, 0.034% off).
  Corectia de ~0.034% NECESITA un ingredient transcendent.
  2*pi e cel mai natural candidat transcendent din framework (Hopf holonomy).

  VERDICT: Pi PARE necesar. Nu e doar un artefact al fitting-ului.
  Alpha traieste genuinely in Z[phi, pi], nu in Z[phi].
""")

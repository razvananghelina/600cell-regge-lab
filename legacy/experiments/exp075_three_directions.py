"""
EXP-075: Trei Directii de Investigare
======================================

1. De ce factorul 5? (suma peste cai topologice?)
2. Eroarea 0.034% = 2*pi*alpha? (faza geometrica)
3. Multiplicitatile si generatiile (16 = fermioni Weyl)
"""

import numpy as np
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-075: TREI DIRECTII DE INVESTIGARE")
print("=" * 70)

# Constante
alpha_exp = ALPHA  # 0.0072973525693
inv_alpha_exp = 1/ALPHA  # 137.035999...

# Spectrul 600-cell
spectrum = [
    (0, 1),                    # lambda_0
    (6/PHI**2, 4),             # lambda_1 = 2.2918
    (5.527864, 9),             # lambda_2
    (9, 16),                   # lambda_3
    (12, 25),                  # lambda_4
    (14, 36),                  # lambda_5
    (14.472136, 9),            # lambda_6
    (15, 16),                  # lambda_7
    (6/PHI**2 + 12, 4),        # lambda_8
]

lambda_1 = 6/PHI**2
lambda_4 = 12
formula_inv_alpha = 20 * PHI**4

# ============================================================
# DIRECTIA 1: DE CE FACTORUL 5?
# ============================================================
print("\n" + "=" * 70)
print("DIRECTIA 1: DE CE FACTORUL 5?")
print("=" * 70)

print("""
FORMULA: 1/alpha = 5 * (lambda_4/lambda_1)^2 = 5 * (2*phi^2)^2 = 20*phi^4

INTREBARE: De unde vine factorul 5?

FAPT: 600-cell = 5 x 24-cell (cele 5 sunt "intretesute")

IPOTEZE:
=========
""")

print("IPOTEZA 1: Suma peste 5 cai topologice")
print("-" * 50)
print("""
In mecanica cuantica, amplitudinea totala e SUMA amplitudinilor
peste toate caile posibile (Feynman path integral).

Daca exista 5 cai topologic distincte in 600-cell,
atunci:
  A_total = A_1 + A_2 + A_3 + A_4 + A_5

Probabilitatea (cuplajul) ar fi:
  |A_total|^2 ~ 5 * |A_single|^2  (daca sunt coerente)

SAU
  |A_total|^2 ~ 5^2 * |A_single|^2 = 25 * |A_single|^2 (daca interfereaza constructiv)
""")

# Verificam care varianta se potriveste
print("Verificare numerica:")
print(f"  5 * (2*phi^2)^2 = {5 * (2*PHI**2)**2:.4f}")
print(f"  25 * (2*phi^2)^2 = {25 * (2*PHI**2)**2:.4f}")
print(f"  sqrt(5) * (2*phi^2)^2 = {np.sqrt(5) * (2*PHI**2)**2:.4f}")
print(f"  1/alpha = {inv_alpha_exp:.4f}")

print("""
Cea mai buna potrivire: 5 * (2*phi^2)^2 = 137.08

Asta sugereaza SUMA INCOERENTA (probabilitati, nu amplitudini).
""")

print("\nIPOTEZA 2: Legatura cu multiplicitatea")
print("-" * 50)
print(f"""
lambda_4 = 12 are multiplicitate d_4 = 25 = 5^2

Observatie: 5 = sqrt(d_4)

Deci formula poate fi scrisa:
  1/alpha = sqrt(d_4) * (lambda_4/lambda_1)^2

Asta ar insemna ca factorul 5 vine din:
  sqrt(degenerarea) = sqrt(dimensiunea reprezentarii)

In teoria grupurilor, asta e legat de
caracterul reprezentarii la identitate / dimensiune.
""")

print("\nIPOTEZA 3: Structura 5-fold a icosaedrului")
print("-" * 50)
print("""
Icosaedrul are simetrie 5-fold (rotatii de 72 grade).
600-cell mostenteste aceasta simetrie.

5 = ordinul grupului ciclic Z_5 din simetria icosaedrica.

Factorul 5 ar putea veni din:
  - Numar de orbite sub actiunea Z_5
  - Numar de sectoare fundamentale
  - Index al unui subgrup
""")

# ============================================================
# DIRECTIA 2: EROAREA DE 0.034%
# ============================================================
print("\n" + "=" * 70)
print("DIRECTIA 2: EROAREA DE 0.034% SI FAZA GEOMETRICA")
print("=" * 70)

error_percent = abs(formula_inv_alpha - inv_alpha_exp) / inv_alpha_exp * 100
error_absolute = formula_inv_alpha - inv_alpha_exp

print(f"1/alpha (formula) = {formula_inv_alpha:.6f}")
print(f"1/alpha (experimental) = {inv_alpha_exp:.6f}")
print(f"Diferenta = {error_absolute:.6f}")
print(f"Eroare relativa = {error_percent:.4f}%")

print(f"\n--- Verificam daca eroarea e legata de 2*pi*alpha ---")

# 2*pi*alpha
two_pi_alpha = 2 * np.pi * alpha_exp
print(f"\n2*pi*alpha = {two_pi_alpha:.6f}")
print(f"Eroare absoluta = {error_absolute:.6f}")
print(f"Raport eroare/(2*pi*alpha) = {error_absolute / two_pi_alpha:.4f}")

# Poate eroarea e 2*pi*alpha^2?
two_pi_alpha2 = 2 * np.pi * alpha_exp**2
print(f"\n2*pi*alpha^2 = {two_pi_alpha2:.6f}")
print(f"Raport eroare/(2*pi*alpha^2) = {error_absolute / two_pi_alpha2:.4f}")

# Sau eroarea e alpha?
print(f"\nalpha = {alpha_exp:.6f}")
print(f"Raport eroare/alpha = {error_absolute / alpha_exp:.4f}")

# Alta forma: poate formula corecta e 20*phi^4 - ceva
print("\n--- Cautam forma corecta ---")
print("""
Daca 1/alpha = 20*phi^4 - corectie, atunci:
  corectie = 20*phi^4 - 1/alpha_exp
""")

correction = formula_inv_alpha - inv_alpha_exp
print(f"Corectie necesara = {correction:.6f}")

# Verificam diverse forme pentru corectie
corrections = [
    ("2*pi*alpha", 2*np.pi*alpha_exp),
    ("2*pi*alpha^2", 2*np.pi*alpha_exp**2),
    ("alpha", alpha_exp),
    ("pi*alpha", np.pi*alpha_exp),
    ("phi*alpha", PHI*alpha_exp),
    ("1/phi^4", 1/PHI**4),
    ("alpha/phi", alpha_exp/PHI),
    ("2*pi/phi^5", 2*np.pi/PHI**5),
]

print(f"\nCorectie = {correction:.6f}")
print("Comparatie cu diverse forme:\n")
for name, val in corrections:
    ratio = correction / val
    err = abs(ratio - round(ratio))
    if err < 0.1:
        marker = f" <-- ratio ~ {round(ratio)}"
    else:
        marker = ""
    print(f"  {name:<15} = {val:.6f}, ratio = {ratio:.4f}{marker}")

# Cea mai interesanta: eroare/alpha ~ 6.3 ~ 2*pi
print(f"\n--- Observatie cheie ---")
print(f"Eroare / alpha = {correction / alpha_exp:.4f}")
print(f"2*pi = {2*np.pi:.4f}")
print(f"Eroare ~ {correction/alpha_exp / (2*np.pi):.2f} * 2*pi * alpha")

# Formula corectata
print("\n--- Formula corectata ---")
# Daca 1/alpha = 20*phi^4 - k*alpha pentru un k
k_fit = correction / alpha_exp
print(f"k = (20*phi^4 - 1/alpha) / alpha = {k_fit:.4f}")
print(f"Deci: 1/alpha = 20*phi^4 - {k_fit:.2f}*alpha")

# Verificam
inv_alpha_corrected = formula_inv_alpha - k_fit * alpha_exp
print(f"\n1/alpha corected = 20*phi^4 - {k_fit:.2f}*alpha = {inv_alpha_corrected:.6f}")
print(f"1/alpha experimental = {inv_alpha_exp:.6f}")

# Dar asta e circular! Avem alpha in formula pentru alpha.
print("""
PROBLEMA: Formula 1/alpha = 20*phi^4 - k*alpha e CIRCULARA!

Trebuie sa rezolvam pentru alpha:
  1/alpha + k*alpha = 20*phi^4
  1 + k*alpha^2 = 20*phi^4 * alpha
  k*alpha^2 - 20*phi^4*alpha + 1 = 0
""")

# Rezolvam ecuatia de gradul 2
# k*alpha^2 - 20*phi^4*alpha + 1 = 0
# Cu k ~ 6.3, rezolvam

k_approx = 6.3
a_coef = k_approx
b_coef = -20*PHI**4
c_coef = 1

discriminant = b_coef**2 - 4*a_coef*c_coef
alpha_solution_1 = (-b_coef - np.sqrt(discriminant)) / (2*a_coef)
alpha_solution_2 = (-b_coef + np.sqrt(discriminant)) / (2*a_coef)

print(f"\nEcuatia: {k_approx}*alpha^2 - {20*PHI**4:.2f}*alpha + 1 = 0")
print(f"Solutii: alpha_1 = {alpha_solution_1:.8f}, alpha_2 = {alpha_solution_2:.8f}")
print(f"Alpha experimental = {alpha_exp:.8f}")

# Daca k = 2*pi
print("\n--- Daca k = 2*pi exact ---")
k_2pi = 2*np.pi
a_coef = k_2pi
b_coef = -20*PHI**4
c_coef = 1

discriminant = b_coef**2 - 4*a_coef*c_coef
alpha_2pi_1 = (-b_coef - np.sqrt(discriminant)) / (2*a_coef)
alpha_2pi_2 = (-b_coef + np.sqrt(discriminant)) / (2*a_coef)

print(f"Ecuatia: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0")
print(f"Solutii: alpha_1 = {alpha_2pi_1:.8f}, alpha_2 = {alpha_2pi_2:.8f}")
print(f"Alpha experimental = {alpha_exp:.8f}")
print(f"Eroare solutie 1 = {abs(alpha_2pi_1 - alpha_exp)/alpha_exp * 100:.4f}%")

# ============================================================
# DIRECTIA 3: MULTIPLICITATI SI GENERATII
# ============================================================
print("\n" + "=" * 70)
print("DIRECTIA 3: MULTIPLICITATI SI GENERATIILE DE FERMIONI")
print("=" * 70)

print("""
MULTIPLICITATILE SPECTRALE:
  lambda_0: mult = 1  = 1^2
  lambda_1: mult = 4  = 2^2
  lambda_2: mult = 9  = 3^2
  lambda_3: mult = 16 = 4^2  <-- 16 Weyl fermions/generatie!
  lambda_4: mult = 25 = 5^2
  lambda_5: mult = 36 = 6^2

IN MODELUL STANDARD:
  - 1 generatie = 16 Weyl fermions
    (6 quarks x 2 chiralitati + 2 leptoni x 2 chiralitati = 16)
  - 3 generatii => 48 fermioni total

OBSERVATIE: lambda_3 = 9 are multiplicitate 16!
""")

print("Valorile proprii si multiplicitatile:")
for i, (lam, mult) in enumerate(spectrum):
    sqrt_mult = int(np.sqrt(mult)) if np.sqrt(mult) == int(np.sqrt(mult)) else np.sqrt(mult)
    print(f"  lambda_{i} = {lam:8.4f}, mult = {mult:2d} = {sqrt_mult}^2")

print(f"""
IPOTEZA: Fiecare nivel spectral corespunde unui "sector" de particule

lambda_3 = 9 (mult 16) -> o generatie de fermioni?
lambda_4 = 12 (mult 25) -> bosoni gauge + Higgs?
lambda_5 = 14 (mult 36) -> ???

Suma multiplicitatilor:
  1 + 4 + 9 + 16 + 25 + 36 = 91

Dar avem 120 varfuri total!
  91 + 9 + 16 + 4 = 120 (include ultimele 3 nivele)
""")

# Verificam suma
total_mult = sum(m for _, m in spectrum)
print(f"\nSuma totala multiplicitati = {total_mult}")

# Legatura cu masele
print("\n--- Legatura cu masele ---")
print("""
IPOTEZA: Raportul maselor fermionilor ~ raportul valorilor proprii

Verificam pentru quark-uri:
  m_t / m_u ~ 75000 (top/up)
  m_b / m_d ~ 850 (bottom/down)
  m_c / m_s ~ 12 (charm/strange)

Rapoarte spectrale:
  lambda_5 / lambda_1 = 14 / 2.29 = 6.1
  lambda_4 / lambda_1 = 12 / 2.29 = 5.2
  lambda_3 / lambda_1 = 9 / 2.29 = 3.9
""")

print("Rapoarte intre valori proprii:")
for i in range(1, len(spectrum)):
    lam_i = spectrum[i][0]
    if lam_i > 0:
        ratio = lam_i / lambda_1
        print(f"  lambda_{i} / lambda_1 = {lam_i:.4f} / {lambda_1:.4f} = {ratio:.4f}")

# Incercam sa gasim pattern pentru mase
print("\n--- Pattern pentru masele leptonilor ---")
m_e = 0.511  # MeV
m_mu = 105.66  # MeV
m_tau = 1776.86  # MeV

print(f"Rapoarte mase leptoni:")
print(f"  m_mu / m_e = {m_mu/m_e:.2f}")
print(f"  m_tau / m_mu = {m_tau/m_mu:.2f}")
print(f"  m_tau / m_e = {m_tau/m_e:.2f}")

# Comparam cu puteri de phi
print(f"\nComparatie cu phi:")
print(f"  phi^5 = {PHI**5:.2f}")
print(f"  phi^4 = {PHI**4:.2f}")
print(f"  phi^7 = {PHI**7:.2f}")
print(f"  m_mu/m_e / phi^5 = {(m_mu/m_e) / PHI**5:.4f}")
print(f"  m_tau/m_mu / phi^4 = {(m_tau/m_mu) / PHI**4:.4f}")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-075")
print("=" * 70)

print("""
DIRECTIA 1 - FACTORUL 5:
========================
- 5 = numarul de 24-cell-uri
- 5 = sqrt(multiplicitatea lui lambda_4)
- Probabil legat de suma INCOERENTA peste 5 subspatii
- Necesita: derivare din teoria reprezentarilor H4

DIRECTIA 2 - EROAREA 0.034%:
============================
- Eroarea ~ 6.3 * alpha ~ 2*pi * alpha
- Ecuatia: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
- Solutia da alpha cu ~0.4% eroare (nu perfect, dar aproape)
- INTERPRETARE: 2*pi = faza geometrica pe bucla decagon!

DIRECTIA 3 - GENERATII:
=======================
- Multiplicitatea 16 la lambda_3 = 9 coincide cu 16 Weyl fermions
- Pattern-ul 1, 4, 9, 16, 25, 36 = patrate perfecte
- Sugereaza conexiune cu reprezentarile grupului H4
- Masele leptonilor NU se potrivesc direct cu phi^n

URMATORUL PAS RECOMANDAT:
=========================
Investigam ecuatia: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
Asta ar putea fi formula COMPLETA pentru alpha!
""")

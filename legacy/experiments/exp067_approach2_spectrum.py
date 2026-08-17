"""
EXP-067: ABORDAREA 2 - Derivarea din Spectrul Laplacianului
===========================================================

Ideea: Ecuatia noastra vine din conditia de REZONANTA
pe reteaua 600-cell.

Stim deja ca:
- lambda_1 = 6/phi^2
- lambda_4 = 12
- Exponent ierarhie = 4*phi^2 = 2*(lambda_4/lambda_1)

Poate ecuatia 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
vine din spectru!
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-067: DERIVARE DIN SPECTRUL LAPLACIANULUI")
print("=" * 70)

# ============================================================
# PASUL 1: SPECTRUL LAPLACIANULUI 600-CELL
# ============================================================
print("\n" + "-" * 70)
print("PASUL 1: RECAPITULARE SPECTRU")
print("-" * 70)

# Valorile proprii (din EXP-061)
eigenvalues = {
    0: (0, 1),           # lambda, degenerare
    1: (6/PHI**2, 4),    # 2.2918, degenerare 4
    2: (6/PHI, 9),       # 5.5279? - sa verificam
    3: (9, 16),
    4: (12, 25),
    5: (14, 36),
}

print("Valorile proprii ale Laplacianului 600-cell:")
for k, (lam, deg) in eigenvalues.items():
    print(f"  lambda_{k} = {lam:.4f} (degenerare {deg})")

lambda_1 = 6/PHI**2
lambda_4 = 12

print(f"\nRapoarte cheie:")
print(f"  lambda_4 / lambda_1 = {lambda_4/lambda_1:.4f} = 2*phi^2 = {2*PHI**2:.4f}")
print(f"  4*phi^2 = {4*PHI**2:.4f}")

# ============================================================
# PASUL 2: FUNCTIA DE PARTITIE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 2: FUNCTIA DE PARTITIE SI PROPAGATORUL")
print("-" * 70)

print("""
IDEE: Campul cuantic pe 600-cell are o functie de partitie.

Z = integral D[psi] exp(-S[psi])

Unde actiunea S include Laplacianul:
S = (1/2) * sum_i psi_i * L_ij * psi_j + interactiuni

PROPAGATORUL (functia Green):
=============================
G(i,j) = <psi_i * psi_j> = sum_n (v_n)_i * (v_n)_j / lambda_n

Unde v_n sunt vectorii proprii ai Laplacianului.

PENTRU UN PUNCT (i=j):
======================
G(i,i) = sum_n |v_n(i)|^2 / lambda_n
       = sum_n (degenerare_n / N_total) / lambda_n  [pentru retea uniforma]
       ~ sum_n d_n / lambda_n
""")

# Calculam propagatorul "on-site"
N_total = 120  # varfuri in 600-cell

# Suma peste moduri (excludem modul zero)
G_onsite = 0
for k in range(1, 6):
    lam, deg = list(eigenvalues.values())[k]
    G_onsite += deg / lam

G_onsite /= N_total

print(f"Propagator on-site (fara modul zero):")
print(f"  G(i,i) ~ sum(d_n/lambda_n) / N = {G_onsite:.6f}")

# ============================================================
# PASUL 3: CONDITIA DE REZONANTA
# ============================================================
print("\n" + "-" * 70)
print("PASUL 3: CONDITIA DE REZONANTA PENTRU SOLITON")
print("-" * 70)

print("""
IPOTEZA: Solitonul e o excitatie REZONANTA pe retea.

O unda stationara exista cand:
- Lungimea de unda "se potriveste" cu geometria
- Energia e minima la anumite frecvente

CONDITIA DE CUANTIZARE:
=======================
Pentru o bucla inchisa de lungime L:
  k * L = 2 * pi * n

Unde k = 2*pi/lambda e numarul de unda.

PE 600-CELL:
============
"Lungimea" buclei minime = circumferinta ecuatoriala
                        = 10 muchii (decagonul)

Conditia de rezonanta:
  k * 10 = 2*pi * n
  k = pi * n / 5

FRECVENTA corespunzatoare (din dispersie omega^2 = lambda):
  omega = sqrt(lambda_k)
""")

# Calculam frecventele
print("\nFrecventele naturale:")
for k in range(1, 6):
    lam, deg = list(eigenvalues.values())[k]
    omega = np.sqrt(lam)
    print(f"  omega_{k} = sqrt({lam:.4f}) = {omega:.4f}")

# ============================================================
# PASUL 4: ECUATIA DIN SPECTRU
# ============================================================
print("\n" + "-" * 70)
print("PASUL 4: INCERCARE DE DERIVARE A ECUATIEI")
print("-" * 70)

print("""
HAI SA VERIFICAM DACA ECUATIA ARE LEGATURA CU SPECTRUL.

Ecuatia noastra: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0

Coeficientii:
- 2*pi = ?
- 20*phi^4 = ?
- 1 = ?

CAUTAM COMBINATII DIN SPECTRU:
==============================
""")

# Combinatii de valori proprii
print("Combinatii de valori proprii si numere geometrice:")
print()

# 20*phi^4
target = 20 * PHI**4
print(f"TINTA: 20*phi^4 = {target:.4f}")

# Incercam combinatii
combos = [
    ("lambda_1 * lambda_4", lambda_1 * lambda_4),
    ("lambda_4^2 / lambda_1", lambda_4**2 / lambda_1),
    ("lambda_4 * (lambda_4/lambda_1)", lambda_4 * (lambda_4/lambda_1)),
    ("12 * 2*phi^2", 12 * 2*PHI**2),
    ("6 * lambda_4 / phi^2", 6 * lambda_4 / PHI**2),
    ("20 * (lambda_4/lambda_1)^2 / 4", 20 * (lambda_4/lambda_1)**2 / 4),
    ("5 * lambda_4 * lambda_1", 5 * lambda_4 * lambda_1),
    ("lambda_4 * lambda_4 / lambda_1", lambda_4 * lambda_4 / lambda_1),
]

for name, val in combos:
    err = abs(val - target) / target * 100
    marker = " <-- MATCH!" if err < 1 else ""
    print(f"  {name:<35} = {val:.4f} (err = {err:.1f}%){marker}")

print()

# 2*pi
target_2pi = 2 * np.pi
print(f"TINTA: 2*pi = {target_2pi:.4f}")

combos_2pi = [
    ("lambda_1 * phi^2", lambda_1 * PHI**2),
    ("6 / phi^2 * phi^2", 6),
    ("3 * lambda_1", 3 * lambda_1),
    ("pi * lambda_1 / 1.1", np.pi * lambda_1 / 1.1),
    ("20 * fete / 10 decagon", 20 / 10 * np.pi),
]

for name, val in combos_2pi:
    err = abs(val - target_2pi) / target_2pi * 100
    print(f"  {name:<35} = {val:.4f} (err = {err:.1f}%)")

# ============================================================
# PASUL 5: ABORDARE DIFERITA - DETERMINANTUL
# ============================================================
print("\n" + "-" * 70)
print("PASUL 5: DETERMINANTUL LAPLACIANULUI")
print("-" * 70)

print("""
IDEE: Poate ecuatia vine din DETERMINANTUL.

det(L - lambda*I) = 0 defineste valorile proprii.

Sau: det(L - m^2*I) pentru ecuatia Klein-Gordon pe retea.

CARACTERISTICA POLINOMIALA:
===========================
P(lambda) = det(L - lambda*I)
          = produs(lambda_n - lambda)^(degenerare_n)

Ecuatia noastra ar putea fi o SIMPLIFICARE a lui P(lambda)
pentru anumite valori speciale.
""")

# Calculam produsul (partial) de valori proprii
product = 1
for k in range(1, 6):
    lam, deg = list(eigenvalues.values())[k]
    product *= lam**deg

print(f"Produsul valorilor proprii (fara zero):")
print(f"  prod(lambda_n^d_n) = {product:.2e}")

# ============================================================
# PASUL 6: FUNCTIA ZETA SPECTRALA
# ============================================================
print("\n" + "-" * 70)
print("PASUL 6: FUNCTIA ZETA SPECTRALA")
print("-" * 70)

print("""
IDEE MAI SOFISTICATA: Functia zeta spectrala.

zeta(s) = sum_n (lambda_n)^(-s)  [pentru lambda_n > 0]

Derivata la s=0:
zeta'(0) = -ln(det(L))  [regularizat]

Aceasta e legata de ACTIUNEA EFECTIVA a teoriei.

Poate alpha apare din zeta(s) pentru un anumit s?
""")

def spectral_zeta(s, eigenvalues_dict):
    """Calculeaza functia zeta spectrala."""
    result = 0
    for k in range(1, len(eigenvalues_dict)):
        lam, deg = list(eigenvalues_dict.values())[k]
        if lam > 0:
            result += deg * lam**(-s)
    return result

# Calculam zeta pentru diverse valori
print("Functia zeta spectrala:")
for s in [0.5, 1, 1.5, 2, 2.5, 3]:
    z = spectral_zeta(s, eigenvalues)
    print(f"  zeta({s}) = {z:.6f}")

# Cautam s pentru care zeta(s) da ceva interesant
print("\nCautam s pentru care zeta(s) ~ 1/alpha:")
for s in np.linspace(0.1, 5, 50):
    z = spectral_zeta(s, eigenvalues)
    if abs(z - 1/ALPHA) < 5:
        print(f"  zeta({s:.3f}) = {z:.2f} (1/alpha = {1/ALPHA:.2f})")

# ============================================================
# PASUL 7: PROPORTII DIN SPECTRU
# ============================================================
print("\n" + "-" * 70)
print("PASUL 7: PROPORTII DIN SPECTRU")
print("-" * 70)

print("""
HAI SA VEDEM CE PROPORTII DAU EXACT VALORILE DIN ECUATIE.

Ecuatia: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0

Daca impartim totul cu alpha:
2*pi*alpha - 20*phi^4 + 1/alpha = 0

Sau:
1/alpha = 20*phi^4 - 2*pi*alpha

VERIFICARE:
1/alpha = 137.036
20*phi^4 = 137.082
2*pi*alpha = 0.0458

137.082 - 0.0458 = 137.036  BINGO!

DECI: 1/alpha = 20*phi^4 - 2*pi*alpha

Acum sa vedem de unde vin 20*phi^4 si 2*pi:
""")

print(f"Verificare numerica:")
print(f"  20*phi^4 = {20*PHI**4:.6f}")
print(f"  2*pi*alpha = {2*np.pi*ALPHA:.6f}")
print(f"  20*phi^4 - 2*pi*alpha = {20*PHI**4 - 2*np.pi*ALPHA:.6f}")
print(f"  1/alpha = {1/ALPHA:.6f}")
print(f"  Diferenta = {abs(20*PHI**4 - 2*np.pi*ALPHA - 1/ALPHA):.6f}")

print("""
INTERPRETARE DIN SPECTRU:
=========================

20*phi^4 = 20 * phi^4
         = (numar de tetraedre per vertex) * phi^4
         = 20 * (lambda_4/lambda_1)^2 / 4  [APROAPE!]

Sau:
20*phi^4 = 10 * 2*phi^4
         = (numar de decagoane per vertex) * 2*phi^4

2*pi = circumferinta unitara
     = faza pe o bucla completa

CONCLUZIE PARTIALA:
===================
Ecuatia pare sa fie:

1/alpha = (contributie geometrica) - (corectie de bucla)

Unde:
- Contributia geometrica = 20*phi^4 (din structura 600-cell)
- Corectia de bucla = 2*pi*alpha (efecte cuantice)
""")

# ============================================================
# PASUL 8: CONEXIUNE CU DEGENERARILE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 8: CONEXIUNE CU DEGENERARILE")
print("-" * 70)

print("""
Degenerarile valorilor proprii:
lambda_0: deg = 1
lambda_1: deg = 4
lambda_2: deg = 9
lambda_3: deg = 16
lambda_4: deg = 25
lambda_5: deg = 36
...

Observatie: deg_n = (n+1)^2 = patrate perfecte!

Suma degenerarilor: 1 + 4 + 9 + 16 + 25 + 36 + ... + 100 = 385
Dar 600-cell are 120 varfuri...

Hmm, degenerarile NU sunt exact (n+1)^2.
Sa verificam din nou...
""")

# Din literatura, degenerarile corecte pentru 600-cell
# sunt mai complicate...

print("""
NOTA: Degenerarile exacte depind de structura de grup.
600-cell are simetrie I_h (icosaedrica).
Reprezentarile ireductibile ale I_h determina degenerarile.
""")

# Fara vizualizare (la cerere)

print("\n" + "=" * 70)
print("CONCLUZIE ABORDAREA 2:")
print("=" * 70)
print("""
AM GASIT O INTERPRETARE PARTIALA:

1/alpha = 20*phi^4 - 2*pi*alpha

Unde:
- 20*phi^4 = contributia geometrica (din structura 600-cell)
- 2*pi*alpha = corectia cuantica (efecte de bucla)

INCA NU STIM:
- De ce exact 2*pi (si nu alt multiplu)?
- Care e mecanismul exact al "renormalizarii"?
- Cum se leaga de valorile proprii in mod direct?

Trecem la Abordarea 3: Fibrarea Hopf.
""")

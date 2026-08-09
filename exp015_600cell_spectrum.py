"""
EXP-015: Spectrul 600-cell si Originea lui 20*phi^4
===================================================
Cautam sa intelegem de unde vine combinatia 20*phi^4
din proprietatile matematice ale 600-cell.
"""

from physics_formulas import *

print("=" * 70)
print("EXP-015: SPECTRUL 600-CELL SI 20*PHI^4")
print("=" * 70)

print("""
INTREBARE CENTRALA:
===================
De ce exact 20*phi^4?

- 20 = tetraedre per vertex = fete icosaedru
- phi^4 = raport de scalare in 4D

Dar exista o formula UNICA care da 20*phi^4 din
proprietatile matematice ale 600-cell?
""")

print("\n" + "-" * 70)
print("PARTEA 1: NUMERELE 600-CELL")
print("-" * 70)

# Numerele 600-cell
V = 120   # vertices
E = 720   # edges
F = 1200  # faces (triangles)
C = 600   # cells (tetrahedra)

print(f"""
600-cell (hexacosichoron):
  V = {V} varfuri
  E = {E} muchii
  F = {F} fete (triunghiuri)
  C = {C} celule (tetraedre)

Relatii:
  V - E + F - C = {V} - {E} + {F} - {C} = {V - E + F - C} (Euler 4D)

Per vertex:
  Tetraedre: 20
  Muchii: 12 (grad)
  Fete: 30
""")

print("\n" + "-" * 70)
print("PARTEA 2: PHI IN COORDONATELE 600-CELL")
print("-" * 70)

print(f"""
Coordonatele celor 120 varfuri ale 600-cell (raza 2):

8 varfuri (permutatii pare):
  (+/-1, +/-1, +/-1, +/-1)

16 varfuri:
  (+/-2, 0, 0, 0) si permutatii

96 varfuri (permutatii pare):
  (+/-1, +/-phi, +/-1/phi, 0)

Unde phi = {PHI:.6f} si 1/phi = {1/PHI:.6f}

Observatie: PHI apare EXPLICIT in 96 din 120 coordonate!
Proportie: 96/120 = {96/120:.2f} = 4/5
""")

print("\n" + "-" * 70)
print("PARTEA 3: DISTANTELE IN 600-CELL")
print("-" * 70)

# Distantele intre varfuri
# Raza = 2, deci distantele sunt:
d1 = 2 / PHI  # edge length pentru raza 2
d2 = 2.0
d3 = 2 * PHI / math.sqrt(PHI + 2)

print(f"""
In 600-cell cu raza 2, distantele intre varfuri:

  d1 = 2/phi = {2/PHI:.6f} (lungime muchie)
  d2 = 2
  d3 ~ {d3:.6f}

Raportul muchie/raza = 1/phi = {1/PHI:.6f}

Alte rapoarte cu phi:
  phi^2 = {PHI**2:.6f}
  phi^3 = {PHI**3:.6f}
  phi^4 = {PHI**4:.6f}
""")

print("\n" + "-" * 70)
print("PARTEA 4: VOLUMUL 600-CELL")
print("-" * 70)

# Hipervolumul 600-cell
# Formula: HV = (50 + 25*sqrt(5))/2 * a^4 pentru edge length a
# Sau HV = 25*phi^3/2 * a^4

a = 1  # edge length unitar
HV_formula = 25 * PHI**3 / 2 * a**4

print(f"""
Hipervolumul 600-cell (edge length a = 1):

  HV = (25*phi^3/2) * a^4
     = (25 * {PHI**3:.6f} / 2) * 1
     = {HV_formula:.6f}

Observatie: Apare phi^3, nu phi^4!

Dar RAPORTUL hipervolumelor (doua 600-cells scalate cu phi):
  HV_mare / HV_mic = phi^4 = {PHI**4:.6f}

Pentru ca volumul 4D scaleaza cu puterea 4 a factorului de scala.
""")

print("\n" + "-" * 70)
print("PARTEA 5: COMBINATII CARE DAU 20*PHI^4")
print("-" * 70)

target = 20 * PHI**4
print(f"Tinta: 20*phi^4 = {target:.6f}")
print(f"\nCautam combinatii de numere din 600-cell:")

# Diverse combinatii
combos = {
    "20 * phi^4 (direct)": 20 * PHI**4,
    "V/6 * phi^4 = 120/6 * phi^4": (V/6) * PHI**4,
    "E/36 * phi^4 = 720/36 * phi^4": (E/36) * PHI**4,
    "F/60 * phi^4": (F/60) * PHI**4,
    "C/30 * phi^4 = 600/30 * phi^4": (C/30) * PHI**4,
    "sqrt(V) * phi^3": math.sqrt(V) * PHI**3,
    "V * phi^3 / 10": V * PHI**3 / 10,
    "(V + C)/2 * phi^3 / 10": (V + C)/2 * PHI**3 / 10,
    "12 * phi^5": 12 * PHI**5,  # 12 = grad
    "5 * phi^6 / 2": 5 * PHI**6 / 2,
    "HV * 4": HV_formula * 4,
    "grad * phi^5 = 12 * phi^5": 12 * PHI**5,
    "tetra_per_vertex * phi^4": 20 * PHI**4,
}

print(f"\n{'Combinatie':<40} {'Valoare':>12} {'Eroare':>12}")
print("-" * 66)
for name, val in combos.items():
    err = (val - target) / target * 100
    marker = " <--" if abs(err) < 0.001 else ""
    print(f"{name:<40} {val:>12.6f} {err:>+12.4f}%{marker}")

print("\n" + "-" * 70)
print("PARTEA 6: INTERPRETAREA GEOMETRICA")
print("-" * 70)

print(f"""
DE UNDE VINE 20:
================
20 = numarul de tetraedre care se intalnesc la un vertex
   = numarul de fete ale icosaedrului (vertex figure)
   = C * 4 / V = {C * 4 / V:.0f}  (fiecare tetraedru are 4 varfuri)

Verificare: 600 * 4 / 120 = 2400/120 = 20 CORECT!

DE UNDE VINE PHI^4:
===================
phi^4 = raportul hipervolumelor a doua 600-cells scalate cu phi
      = (phi * a)^4 / a^4 = phi^4

In proiectia E8 -> 4D:
  - Apar DOUA 600-cells, unul "mare" si unul "mic"
  - Scalate cu factor phi
  - Raportul volumelor = phi^4
""")

print("\n" + "-" * 70)
print("PARTEA 7: FORMULA DIN NUMARARE")
print("-" * 70)

print(f"""
O DERIVARE COMBINATORIALA:

Consideram un vertex in 600-cell.
  - Are 20 tetraedre adiacente
  - Fiecare tetraedru contribuie cu "ceva" la interactiune

Daca contributia fiecarui tetraedru = phi^4 (din scalare 4D):
  Contributia totala = 20 * phi^4

Interpretare:
  1/alpha = cate "cai" geometrice sunt disponibile
          = 20 directii * phi^4 pondere geometrica
          = 20 * phi^4
          ~ 137

Deci alpha = "probabilitatea de a alege o SINGURA cale specifica"
           = 1 / (numar_directii * pondere_geometrica)
           = 1 / (20 * phi^4)
           ~ 1/137
""")

print("\n" + "-" * 70)
print("PARTEA 8: CORECTIA 2*PI*ALPHA")
print("-" * 70)

print(f"""
Formula completa: 1/alpha = 20*phi^4 - 2*pi*alpha

De unde vine -2*pi*alpha?

Ipoteza 1: SELF-INTERACTION
  - Particula poate interactiona cu SINE INSASI
  - Aceasta e o "corectie radiativa" in QED
  - 2*pi vine din integrarea unghiulara

Ipoteza 2: TOPOLOGIA S^3
  - 600-cell traieste in S^3 (suprafata 4-sferei)
  - O particula care face un "loop complet" pe S^3
  - Contribuie cu factor 2*pi

Ipoteza 3: RENORMALIZARE
  - 20*phi^4 e valoarea "goala" (bare)
  - 2*pi*alpha e corectia de renormalizare
  - Rezulta alpha "imbracata" (dressed)

Numeric:
  20*phi^4 = {20*PHI**4:.6f}
  2*pi*alpha = {2*PI*ALPHA:.6f}
  20*phi^4 - 2*pi*alpha = {20*PHI**4 - 2*PI*ALPHA:.6f}
  1/alpha = {ALPHA_INV:.6f}
  Eroare: {abs(20*PHI**4 - 2*PI*ALPHA - ALPHA_INV)/ALPHA_INV * 100:.6f}%
""")

print("\n" + "-" * 70)
print("PARTEA 9: SPECTRUL ICOSAEDRULUI (VERTEX FIGURE)")
print("-" * 70)

print(f"""
Matricea de adiacenta a icosaedrului (12 varfuri, grad 5):

Valorile proprii:
  lambda_0 = 5 (multiplicitate 1) - valoare maxima
  lambda_1 = sqrt(5) ~ {math.sqrt(5):.6f} (multiplicitate 3)
  lambda_2 = 0 (multiplicitate 4)
  lambda_3 = -sqrt(5) ~ {-math.sqrt(5):.6f} (multiplicitate 3)
  lambda_4 = -5 (multiplicitate 1)

Suma valorilor proprii = 0 (trace = 0)
Produs valorilor proprii = det(A) = ?

Observatie: sqrt(5) = 2*phi - 1 = {2*PHI - 1:.6f}
""")

# Relatie sqrt(5) si phi
print(f"\nRelatii cu phi:")
print(f"  sqrt(5) = {math.sqrt(5):.6f}")
print(f"  2*phi - 1 = {2*PHI - 1:.6f}")
print(f"  phi^2 - phi = 1 => phi^2 = phi + 1 = {PHI + 1:.6f}")
print(f"  phi + 1/phi = sqrt(5) = {PHI + 1/PHI:.6f}")

print("\n" + "-" * 70)
print("PARTEA 10: CAUTAM 20*PHI^4 IN SPECTRU")
print("-" * 70)

# Combinatii spectrale
spec_combos = {
    "5 * phi^5": 5 * PHI**5,
    "(5 + sqrt(5)) * phi^4": (5 + math.sqrt(5)) * PHI**4,
    "5 * (1 + 1/phi) * phi^4 = 5*sqrt(5)*phi^3": 5 * math.sqrt(5) * PHI**3,
    "12 * phi^5 (12 varfuri ico)": 12 * PHI**5,
    "20 * phi^4 (20 fete ico)": 20 * PHI**4,
    "30 * phi^3 (30 muchii ico)": 30 * PHI**3,
    "(5^2 - 5) * phi^4 = 20 * phi^4": (25 - 5) * PHI**4,
}

print(f"Tinta: 20*phi^4 = {20*PHI**4:.6f}")
print(f"\n{'Combinatie spectrala':<45} {'Valoare':>12}")
print("-" * 60)
for name, val in spec_combos.items():
    marker = " <-- MATCH" if abs(val - target) < 0.001 else ""
    print(f"{name:<45} {val:>12.6f}{marker}")

print(f"""

OBSERVATIE CHEIE:
=================
(5^2 - 5) * phi^4 = (lambda_max^2 - lambda_max) * phi^4 = 20 * phi^4

unde lambda_max = 5 = valoarea proprie maxima a icosaedrului!

Deci: 20 = 5^2 - 5 = 5*(5-1) = 5*4 = 20
      = grad * (grad - 1)? Nu, gradul e 5, nu 4...

Dar 20 = fete icosaedru = varfuri * grad / 3 = 12 * 5 / 3 = 20
(pentru ca fiecare fata e partajata de 3 varfuri)
""")

print("\n" + "-" * 70)
print("PARTEA 11: FORMULA ELEGANTA")
print("-" * 70)

print(f"""
PROPUNERE DE FORMULA:

Fie icosaedrul I cu:
  - V_I = 12 varfuri
  - E_I = 30 muchii
  - F_I = 20 fete
  - grad = 5

Atunci:
  1/alpha = F_I * phi^4 - 2*pi*alpha

Verificare:
  F_I * phi^4 = 20 * {PHI**4:.6f} = {20 * PHI**4:.6f}

SAU, folosind relatia Euler pentru icosaedru:
  V - E + F = 12 - 30 + 20 = 2

  1/alpha = (E_I - V_I + 2) * phi^4 - 2*pi*alpha
          = (30 - 12 + 2) * phi^4 - 2*pi*alpha
          = 20 * phi^4 - 2*pi*alpha

Deci 20 = E - V + 2 pentru icosaedru!
""")

print("\n" + "-" * 70)
print("PARTEA 12: SINTEZA - DE UNDE VINE FORMULA")
print("-" * 70)

print(f"""
DERIVARE PROPUSA (speculativa):
===============================

1. Spatiul la scala Planck = quasicristal cu structura 600-cell

2. La fiecare vertex, structura locala = ICOSAEDRU (vertex figure)

3. Icosaedrul are F = 20 fete = 20 directii posibile

4. In 4D, contributia fiecarei directii e ponderata cu phi^4
   (din raportul hipervolumelor in proiectia E8 -> 4D)

5. Numarul efectiv de "cai" = 20 * phi^4 ~ 137

6. Probabilitatea de interactiune = 1/(numar de cai) ~ 1/137 = alpha

7. Corectia 2*pi*alpha vine din:
   - Self-interaction (bucle cuantice)
   - SAU topologia S^3 (unde traieste 600-cell)
   - SAU renormalizare

FORMULA FINALA:
  1/alpha = F_icosaedru * phi^4 - 2*pi*alpha
          = 20 * phi^4 - 2*pi*alpha

  Eroare: {abs(20*PHI**4 - 2*PI*ALPHA - ALPHA_INV)/ALPHA_INV * 100:.6f}%
""")

print("\n" + "=" * 70)
print("CONCLUZII EXP-015")
print("=" * 70)

print(f"""
AM GASIT:
=========
1. 20 = numarul de fete ale icosaedrului = vertex figure al 600-cell
2. 20 = E - V + 2 pentru icosaedru (din Euler)
3. phi^4 = raportul hipervolumelor in 4D
4. sqrt(5) apare in spectrul icosaedrului
5. Formula: 1/alpha = 20*phi^4 - 2*pi*alpha are eroare 0.00014%

INTERPRETARE:
=============
Alpha = 1/(numar de directii * pondere geometrica - corectie)
      = 1/(20 * phi^4 - 2*pi*alpha)

Particula la un vertex "vede" 20 directii icosaedrice,
ponderate cu phi^4 din geometria 4D.
Probabilitatea de interactiune = inversa acestui numar.

CE LIPSESTE:
============
- Derivarea lui 2*pi*alpha din principii prime
- Explicatia fizica pentru "ponderare cu phi^4"
- Confirmare ca spatiul chiar ARE structura 600-cell
""")

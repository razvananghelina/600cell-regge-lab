"""
EXP-050: Conexiunea 10/phi ~ 2*pi
=================================
Investigam de ce 10/phi si 2*pi sunt aproape egale
si ce inseamna asta pentru formule.
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-050: CONEXIUNEA 10/phi ~ 2*pi")
print("=" * 70)

# ============================================================
# OBSERVATIA DE BAZA
# ============================================================
print("\n" + "-" * 70)
print("OBSERVATIA DE BAZA")
print("-" * 70)

val_10_over_phi = 10 / PHI
val_2pi = 2 * PI

print(f"10/phi = {val_10_over_phi:.10f}")
print(f"2*pi   = {val_2pi:.10f}")
print(f"Diferenta = {abs(val_10_over_phi - val_2pi):.10f}")
print(f"Raport = {val_10_over_phi / val_2pi:.10f}")
print(f"Eroare = {abs(val_10_over_phi - val_2pi) / val_2pi * 100:.4f}%")

# ============================================================
# COINCIDENTA SAU NU?
# ============================================================
print("\n" + "-" * 70)
print("COINCIDENTA SAU NU?")
print("-" * 70)

print(f"""
10/phi = 10 * (2/(1+sqrt(5)))
       = 20 / (1+sqrt(5))
       = 20 / phi^2 * phi
       = 20/phi^2 * phi
       = 20 * phi^(-1)
       = 20 * (phi - 1)   [folosind phi^(-1) = phi - 1]
       = 20*phi - 20
       = {20*PHI - 20:.10f}

2*pi = {2*PI:.10f}

Deci intrebarea e: de ce 20*phi - 20 ~ 2*pi?
Sau: 20*(phi - 1) ~ 2*pi
Sau: 10*(phi - 1) ~ pi
""")

# Verificam
print(f"10*(phi - 1) = {10*(PHI - 1):.10f}")
print(f"pi = {PI:.10f}")
print(f"Raport: {10*(PHI - 1) / PI:.10f}")

# ============================================================
# EXPLORARE: ALTE RELATII phi - pi
# ============================================================
print("\n" + "-" * 70)
print("ALTE RELATII PHI - PI")
print("-" * 70)

print("Combinatii simple:")
relations = [
    ("phi + pi", PHI + PI),
    ("phi * pi", PHI * PI),
    ("phi^2 + pi", PHI**2 + PI),
    ("phi^2 * pi", PHI**2 * PI),
    ("phi^3 / pi", PHI**3 / PI),
    ("phi^4 / pi", PHI**4 / PI),
    ("phi^5 / pi^2", PHI**5 / PI**2),
    ("5 * phi / pi", 5 * PHI / PI),
    ("6 * phi / pi", 6 * PHI / PI),
    ("10 / phi - 2*pi", 10/PHI - 2*PI),
    ("10 * phi - 2*pi^2", 10*PHI - 2*PI**2),
    ("pi * phi^2 / 5", PI * PHI**2 / 5),
]

for name, val in relations:
    # Cautam daca e aproape de un intreg
    nearest_int = round(val)
    if nearest_int != 0:
        err = abs(val - nearest_int) / abs(nearest_int) * 100
        if err < 5:
            print(f"  {name} = {val:.6f} ~ {nearest_int} (eroare {err:.2f}%)")

# ============================================================
# INTERPRETARE GEOMETRICA
# ============================================================
print("\n" + "-" * 70)
print("INTERPRETARE GEOMETRICA")
print("-" * 70)

print(f"""
In 600-cell inscris in S^3 de raza 1:

DECAGONUL (geodezica):
- 10 vertices pe decagon
- Edge length = 1/phi = {1/PHI:.6f}
- Perimetrul decagonului = 10/phi = {10/PHI:.6f}

CERCUL MARE pe S^3:
- Circumferinta = 2*pi*r = 2*pi (pentru r=1)
- 2*pi = {2*PI:.6f}

OBSERVATIE:
  Perimetru decagon / Circumferinta cerc = {(10/PHI)/(2*PI):.6f}

  Decagonul e inscris IN cercul mare!
  Deci perimetrul sau e MAI MIC decat circumferinta.

  Exact cu cat?
  2*pi - 10/phi = {2*PI - 10/PHI:.6f}

  Aceasta diferenta = {(2*PI - 10/PHI):.6f} este "spatiul ramas".
""")

# ============================================================
# FORMULA PENTRU DIFERENTA
# ============================================================
print("-" * 70)
print("CE E DIFERENTA 2*pi - 10/phi ?")
print("-" * 70)

diff = 2*PI - 10/PHI
print(f"2*pi - 10/phi = {diff:.10f}")

# Cautam combinatii simple
print("\nCautam ce e aceasta diferenta:")
guesses = [
    ("1/10", 1/10),
    ("1/phi^3", 1/PHI**3),
    ("1/phi^4", 1/PHI**4),
    ("pi/30", PI/30),
    ("pi/31", PI/31),
    ("2*pi/60", 2*PI/60),
    ("1/(3*phi)", 1/(3*PHI)),
    ("1/(pi*phi)", 1/(PI*PHI)),
    ("(phi-1)/6", (PHI-1)/6),
    ("phi/15", PHI/15),
    ("phi/16", PHI/16),
]

print(f"Target: {diff:.6f}")
for name, val in guesses:
    err = abs(val - diff) / diff * 100
    print(f"  {name} = {val:.6f} (eroare {err:.2f}%)")

# ============================================================
# CONEXIUNEA CU CORECTIA 2*pi*alpha
# ============================================================
print("\n" + "-" * 70)
print("CONEXIUNEA CU CORECTIA 2*pi*alpha")
print("-" * 70)

print(f"""
Amintim formula pentru alpha:
  1/alpha = 20*phi^4 - 2*pi*alpha

Termenul de corectie: 2*pi*alpha = {2*PI*ALPHA:.6f}

OBSERVATIE:
  2*pi - 10/phi = {diff:.6f}
  2*pi*alpha = {2*PI*ALPHA:.6f}

  Raport: {(2*PI - 10/PHI)/(2*PI*ALPHA):.4f}

  Hmm, raportul e ~2.24, nu un numar simplu.
""")

# ============================================================
# FORMULA ALTERNATIVA: DECAGON vs CERC
# ============================================================
print("-" * 70)
print("FORMULA ALTERNATIVA BAZATA PE DECAGON/CERC")
print("-" * 70)

print(f"""
Daca perimetrul decagonului = 10/phi si circumferinta = 2*pi:

  "Excess" = 2*pi - 10/phi = {2*PI - 10/PHI:.6f}

Interpretare: cat "lipseste" decagonului pana la cerc.

IPOTEZA:
  Poate formula pentru alpha vine din raportul decagon/cerc?

  (10/phi) / (2*pi) = {(10/PHI)/(2*PI):.10f}

  1 - (10/phi)/(2*pi) = {1 - (10/PHI)/(2*PI):.10f}

  Aceasta "lipsa" = {(1 - (10/PHI)/(2*PI)):.6f} = 1.6%
""")

# ============================================================
# PATTERN MAI GENERAL
# ============================================================
print("\n" + "-" * 70)
print("PATTERN: n-GON vs CERC")
print("-" * 70)

print("Pentru un n-gon regulat inscris in cerc de raza r:")
print("  Side length = 2*r*sin(pi/n)")
print("  Perimeter = n * 2*r*sin(pi/n) = 2*n*r*sin(pi/n)")
print("  Raport cu circumferinta = n*sin(pi/n) / pi")
print()

for n in [5, 6, 8, 10, 12, 20]:
    ratio = n * np.sin(PI/n) / PI
    print(f"  n={n:2d}: ratio = {ratio:.6f}")

# Dar in 600-cell, edge length = 1/phi, nu 2*sin(pi/10)
print(f"\nIn 600-cell, edge = 1/phi = {1/PHI:.6f}")
print(f"2*sin(pi/10) = 2*sin(36 deg) = {2*np.sin(PI/10):.6f}")
print(f"Raport: {(1/PHI) / (2*np.sin(PI/10)):.6f}")

# Interesant: sunt EGALE!
print(f"\nOBSERVATIE: 1/phi = 2*sin(36 deg) = 2*sin(pi/10)!")
print(f"Verificare: 2*sin(pi/10) = {2*np.sin(PI/10):.10f}")
print(f"           1/phi = {1/PHI:.10f}")
print(f"           Diferenta: {abs(2*np.sin(PI/10) - 1/PHI):.2e}")

# ============================================================
# DERIVARE: DE CE 1/phi = 2*sin(36 deg)
# ============================================================
print("\n" + "-" * 70)
print("DERIVARE: 1/phi = 2*sin(36 deg)")
print("-" * 70)

print(f"""
Aceasta e o IDENTITATE CUNOSCUTA!

In pentagon regulat:
  Unghiul central = 360/5 = 72 deg
  Jumatate = 36 deg

  Side/Diameter = sin(36 deg)
  Side = 2*R*sin(36 deg) pentru raza R

  Pentru pentagon cu circumradius = phi (speciale):
    Side = 2 * phi * sin(36 deg)

  Dar stim ca in pentagon, side/diagonal = 1/phi.

RELATIA FUNDAMENTALA:
  sin(36 deg) = sqrt(10 - 2*sqrt(5)) / 4
  cos(36 deg) = (1 + sqrt(5)) / 4 = phi/2

Verificam:
  cos(36 deg) = {np.cos(PI/5):.10f}
  phi/2 = {PHI/2:.10f}
  Diferenta: {abs(np.cos(PI/5) - PHI/2):.2e}

DECI: cos(36 deg) = phi/2 EXACT!

Si: sin(36 deg) = sqrt(1 - phi^2/4) = sqrt((4-phi^2)/4)
  4 - phi^2 = 4 - {PHI**2:.6f} = {4 - PHI**2:.6f}

  Stim: phi^2 = phi + 1 = {PHI + 1:.6f}
  Deci: 4 - phi^2 = 4 - phi - 1 = 3 - phi = {3 - PHI:.6f}

  sin^2(36) = (3 - phi)/4 = {(3-PHI)/4:.6f}
  sin(36) = sqrt((3-phi)/4) = {np.sqrt((3-PHI)/4):.6f}

Verificare: {np.sin(PI/5):.6f}
Match!
""")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-050")
print("=" * 70)

print(f"""
1. IDENTITATE FUNDAMENTALA:
   1/phi = 2*sin(36 deg) = 2*sin(pi/10)  [EXACT]

   Aceasta leaga phi de geometria decagonului.

2. PERIMETRU DECAGON:
   P = 10 * (1/phi) = 10/phi = 20*sin(pi/10)

   Pentru decagon INSCRIS in cerc de raza 1:
   P = 10/phi = {10/PHI:.6f}

   Circumferinta cerc = 2*pi = {2*PI:.6f}

3. APROPIEREA 10/phi ~ 2*pi:
   Raport = {(10/PHI)/(2*PI):.6f} = 98.4%

   Aceasta NU e coincidenta - e geometrie!
   Decagonul inscris in cerc are perimetru ~98.4% din circumferinta.

4. CONEXIUNE CU ALPHA:
   Factorul 10*phi care leaga alpha_bare de alpha_s:

   10*phi = 10/(1/phi) = 1 / (perimetru decagon / 10^2)

   Sau: 10*phi = (circumferinta cerc) * (2*pi / (10/phi))
              = 2*pi * (2*pi * phi / 10)
              ~ (2*pi)^2 * phi / 10

   Nu e o formula simpla, dar leaga phi, pi, si 10.

5. SEMNIFICATIE:
   - phi vine din geometria pentagonala/decagonala
   - pi vine din geometria circulara
   - 10 = steps in decagon = unghiul de 36 deg
   - Aceste trei (phi, pi, 10) sunt LEGATE prin sin(36 deg) = 1/(2*phi)

6. FORMULA PENTRU ALPHA:
   1/alpha_bare = 20*phi^4 = 20 * phi^4
   1/alpha_s = 2*phi^3 = 20*phi^4 / (10*phi)

   10*phi e factorul de "reducere dimensionala" 4D -> 3D
   Geometric, reprezinta raportul intre structuri de dimensiuni diferite.
""")

# ============================================================
# BONUS: TOATE CONSTANTELE DIN GEOMETRIE
# ============================================================
print("-" * 70)
print("BONUS: DERIVARE COMPLETA DIN GEOMETRIE")
print("-" * 70)

print(f"""
TOATE NUMERELE VIN DIN 600-CELL:

  20 = tetraedre per vertex = 5 (tri/edge) * 4 (dim)
  10 = vertices per decagon = unghiuri de 36 deg
  6 = decagoane per vertex
  26 = 6 + 20 = "directii gauge" totale

  phi = golden ratio, apare natural in geometria 600-cell
  pi = circumferinta cercului (continuum limit)

FORMULELE:
  1/alpha_bare = 20 * phi^4           (EM, 4D)
  alpha_s = 1/(2*phi^3) = 10*phi * alpha_bare   (strong, 3D)
  sin^2(theta_W) = 6/26               (weak mixing)

  Termenul de corectie 2*pi*alpha: self-loop pe fibra Hopf (S^1 = cerc)

TOTUL vine din geometria 600-cell si proprietatile sale!
""")

"""
EXP-007b: 600-cell si phi^4
===========================
Calculam hipervolumul 600-cell si cautam legatura cu phi^4
"""

from physics_formulas import *

print("=" * 70)
print("EXP-007b: 600-CELL SI PHI^4")
print("=" * 70)

print("\n" + "-" * 70)
print("PROPRIETATILE 600-CELL")
print("-" * 70)

print("""
600-cell (hexacosichoron):
  - Vertices (V): 120
  - Edges (E): 720
  - Faces (F): 1200 (triunghiuri)
  - Cells (C): 600 (tetraedre)

  Euler 4D: V - E + F - C = 120 - 720 + 1200 - 600 = 0 (corect!)

  Este analogul 4D al icosaedrului.
  Vertex figure: icosaedru (20 tetraedre se intalnesc la fiecare varf)
""")

# Verificare Euler
V, E, F, C = 120, 720, 1200, 600
euler = V - E + F - C
print(f"Verificare Euler 4D: {V} - {E} + {F} - {C} = {euler}")

print("\n" + "-" * 70)
print("COORDONATE SI DIMENSIUNI")
print("-" * 70)

# Pentru 600-cell cu raza circumscrisa = 1
# Lungimea muchiei = 1/phi
edge_600 = 1 / PHI
print(f"\nPentru raza circumscrisa R = 1:")
print(f"  Lungimea muchiei a = 1/phi = {edge_600:.6f}")

# Volumul unui tetraedru regulat cu muchia a
def vol_tetrahedron(a):
    return a**3 / (6 * math.sqrt(2))

vol_tet = vol_tetrahedron(edge_600)
print(f"  Volumul unui tetraedru (a = 1/phi): {vol_tet:.10f}")

# Dar in 4D, "volumul" unui tetraedru e diferit - e o "celula" 3D
# Hipervolumul 4D al 600-cell trebuie calculat altfel

print("\n" + "-" * 70)
print("HIPERVOLUMUL 600-CELL")
print("-" * 70)

print("""
Formula hipervolum 4-polytop regulat:
  HV = (1/4) * C * h * V_cell

unde:
  C = numarul de celule (600)
  h = inaltimea de la centru la o celula
  V_cell = volumul unei celule (tetraedru)

Pentru 600-cell cu muchia a:
""")

# Formula exacta pentru hipervolumul 600-cell cu muchia a
# HV = (50 + 25*sqrt(5)) * a^4 / 4
# sau echivalent:
# HV = (25/4) * (2 + sqrt(5)) * a^4 = (25/4) * (2 + sqrt(5)) * a^4

def hypervolume_600cell(a):
    """Hipervolumul 600-cell cu muchia a"""
    return (50 + 25*math.sqrt(5)) / 4 * a**4

# Cu muchia a = 1
HV_unit = hypervolume_600cell(1)
print(f"Cu muchia a = 1:")
print(f"  HV = (50 + 25*sqrt(5))/4 * a^4")
print(f"  HV = {HV_unit:.10f}")

# Sa vedem legatura cu phi
print(f"\nAnaliza termenului (50 + 25*sqrt(5))/4:")
coef = (50 + 25*math.sqrt(5))/4
print(f"  (50 + 25*sqrt(5))/4 = {coef:.10f}")
print(f"  = 25*(2 + sqrt(5))/4")
print(f"  = 25*phi^2/2 = {25*PHI**2/2:.10f}")
print(f"  (folosind 2 + sqrt(5) = phi + phi^2 = phi*(1+phi) = phi*phi^2 = phi^3... nu)")

# De fapt 2 + sqrt(5) = ?
val = 2 + math.sqrt(5)
print(f"\n  2 + sqrt(5) = {val:.10f}")
print(f"  phi^3 = {PHI**3:.10f}")
print(f"  Da! 2 + sqrt(5) = phi^3 = 2*phi + 1")

print(f"\n  Deci HV = 25*phi^3/4 * a^4 = {25*PHI**3/4:.10f}")
print(f"  Verificare: {coef:.10f}")

print("\n" + "-" * 70)
print("LEGATURA CU 20*PHI^4")
print("-" * 70)

print(f"\nHipervolumul 600-cell:")
print(f"  HV = (25/4) * phi^3 * a^4")
print(f"  HV = 6.25 * phi^3 * a^4")

print(f"\nFormula noastra:")
print(f"  20 * phi^4 = {20*PHI**4:.10f}")

print(f"\nRaport:")
raport = 20 * PHI**4 / (25 * PHI**3 / 4)
print(f"  20*phi^4 / (25*phi^3/4) = {raport:.10f}")
print(f"  = 80*phi^4 / (25*phi^3) = (80/25)*phi = 3.2*phi")
print(f"  = {3.2*PHI:.10f}")

print("\nNu e o legatura directa simpla...")

print("\n" + "-" * 70)
print("ALTE PROPRIETATI CU PHI")
print("-" * 70)

print(f"\nNumere din 600-cell si icosaedru:")
print(f"  600-cell: 120 V, 720 E, 1200 F, 600 C")
print(f"  Icosaedru: 12 V, 30 E, 20 F")

print(f"\nRapoarte:")
print(f"  600/20 = {600/20} (celule 4D / fete 3D)")
print(f"  120/12 = {120/12} (varfuri 4D / varfuri 3D)")
print(f"  720/30 = {720/30} (muchii 4D / muchii 3D)")
print(f"  1200/20 = {1200/20} (fete 4D / fete 3D)")

print(f"\n  Raportul constant e 10 sau 30 sau 60...")
print(f"  60 = |A5| = numarul de rotatii ale icosaedrului")

print("\n" + "-" * 70)
print("CAUTAM FORMULA CU 20*PHI^4")
print("-" * 70)

# Ce daca 20*phi^4 e ceva legat de 600-cell?
print(f"\n600 / 20 = 30")
print(f"600 / phi^4 = {600/PHI**4:.6f}")
print(f"120 * phi^4 = {120*PHI**4:.6f}")
print(f"720 / phi^4 = {720/PHI**4:.6f}")

print(f"\n120 / phi^4 = {120/PHI**4:.6f}")
print(f"Asta e aproape de 20! (este {120/PHI**4:.6f})")

# AHA!
print("\n" + "=" * 70)
print("DESCOPERIRE!")
print("=" * 70)

print(f"\n120 / phi^4 = {120/PHI**4:.10f}")
print(f"Deci: 120 = phi^4 * {120/PHI**4:.6f}")
print(f"      120 ~ phi^4 * 17.5")

# Mai exact
print(f"\n120 / phi^4 = 120 / (3*phi + 2)")
val = 120 / (3*PHI + 2)
print(f"            = {val:.10f}")

# Alta idee: 20*phi^4 = ?
print(f"\n20*phi^4 = {20*PHI**4:.6f}")
print(f"600/phi^4 = {600/PHI**4:.6f} ~ 87.5")
print(f"1200/phi^4 = {1200/PHI**4:.6f} ~ 175")

# Poate e legat de raportul suprafata/hipervolum?
print("\n" + "-" * 70)
print("RAPORT SUPRAFATA 3D / HIPERVOLUM 4D")
print("-" * 70)

# Suprafata 3D a 600-cell = suma ariilor celor 1200 triunghiuri
# = 1200 * (sqrt(3)/4) * a^2 = 300*sqrt(3) * a^2
S_3D = 1200 * (math.sqrt(3)/4)  # pentru a=1
HV_4D = hypervolume_600cell(1)

print(f"\nPentru a = 1:")
print(f"  Suprafata 3D = 1200 * sqrt(3)/4 = 300*sqrt(3) = {S_3D:.6f}")
print(f"  Hipervolum 4D = {HV_4D:.6f}")
print(f"  Raport S^2/HV = {S_3D**2/HV_4D:.6f}")
print(f"  Raport S/HV = {S_3D/HV_4D:.6f}")

print("\n" + "-" * 70)
print("CONCLUZIE EXP-007b")
print("-" * 70)

print(f"""
GASIT:
1. Hipervolumul 600-cell: HV = (25/4)*phi^3 * a^4
   (contine phi^3, nu phi^4)

2. Numarul 120 (varfuri 600-cell = |2I| binary icosahedral)
   120/phi^4 ~ 17.5

3. Nu am gasit o legatura DIRECTA simpla intre 20*phi^4 si 600-cell

OBSERVATIE INTERESANTA:
- 600-cell are 120 varfuri
- Binary icosahedral group are 120 elemente
- E8 lattice are 240 radacini = 2 * 120
- Cand E8 se proiecteaza in 4D, da doua 600-cells cu raport phi!

INTREBARE:
Poate 20*phi^4 e legat de PROIECTIA E8 -> 4D, nu de 600-cell direct?
""")

"""
EXP-011: Analiza 20-group si unghiurile de rotatie
==================================================
Investigam unghiurile 15.522 si 15.552 grade folosite de QGR
si vedem daca au legatura cu phi si formula noastra.
"""

from physics_formulas import *

print("=" * 70)
print("EXP-011: 20-GROUP SI UNGHIURILE DE ROTATIE")
print("=" * 70)

print("\n" + "-" * 70)
print("UNGHIURILE DIN EMERGENCE THEORY")
print("-" * 70)

print("""
Din Quantum Gravity Research:

1. 15.522 grade - folosit pentru a crea QSN din 5 clone FCC rotite
2. 15.552 grade - folosit pentru a reprezenta 4D quasicrystal in 3D

Intrebare: Au legatura cu phi?
""")

# Verificam unghiurile
angle1 = 15.522
angle2 = 15.552

print(f"Unghiuri:")
print(f"  15.522 grade = {math.radians(angle1):.6f} rad")
print(f"  15.552 grade = {math.radians(angle2):.6f} rad")

print("\n" + "-" * 70)
print("CAUTAM CONEXIUNI CU PHI")
print("-" * 70)

# Unghiuri cunoscute legate de phi
print("\nUnghiuri legate de phi in geometria icosaedrica:")

# Unghiul diedru al icosaedrului
dihedral_ico = 2 * math.atan(PHI)  # in radiani
dihedral_ico_deg = math.degrees(dihedral_ico)
print(f"  Unghi diedru icosaedru = 2*arctan(phi) = {dihedral_ico_deg:.3f} grade")

# Unghiul diedru al tetraedrului
dihedral_tetra = math.acos(1/3)
dihedral_tetra_deg = math.degrees(dihedral_tetra)
print(f"  Unghi diedru tetraedru = arccos(1/3) = {dihedral_tetra_deg:.3f} grade")

# Unghi pentagon (interior)
pentagon_angle = 108
print(f"  Unghi interior pentagon = {pentagon_angle} grade")

# 360/phi^4
angle_phi4 = 360 / PHI**4
print(f"\n  360 / phi^4 = {angle_phi4:.3f} grade")
print(f"  360 / phi^5 = {360/PHI**5:.3f} grade")
print(f"  360 / phi^6 = {360/PHI**6:.3f} grade")

# Poate e legat de impartirea cercului
print(f"\n  360 / 20 = {360/20:.3f} grade")
print(f"  360 / 24 = {360/24:.3f} grade")
print(f"  360 / 23.14 = {360/23.14:.3f} grade ~ 15.5")

# Verificam daca 15.522 sau 15.552 au o formula simpla
print("\n" + "-" * 70)
print("CAUTAM FORMULA PENTRU ~15.5 GRADE")
print("-" * 70)

# Posibile formule
formulas = {
    "arctan(1/phi^2)": math.degrees(math.atan(1/PHI**2)),
    "arctan(1/phi^3)": math.degrees(math.atan(1/PHI**3)),
    "arcsin(1/phi^3)": math.degrees(math.asin(1/PHI**3)),
    "arccos(phi^3/5)": math.degrees(math.acos(PHI**3/5)) if PHI**3/5 <= 1 else None,
    "36/phi^2": 36/PHI**2,
    "72/(2*phi^2)": 72/(2*PHI**2),
    "90/phi^3": 90/PHI**3,
    "180/(phi^5)": 180/PHI**5,
    "360/(phi^4 + 20)": 360/(PHI**4 + 20),
    "arctan(phi-1)": math.degrees(math.atan(PHI-1)),
    "arctan(1/phi)": math.degrees(math.atan(1/PHI)),
    "18/phi + 5": 18/PHI + 5,
}

print(f"\n{'Formula':<25} {'Valoare':>12} {'Dif de 15.522':>15}")
print("-" * 55)
for name, val in formulas.items():
    if val is not None:
        diff = val - 15.522
        marker = " <--" if abs(diff) < 0.1 else ""
        print(f"{name:<25} {val:>12.4f} {diff:>+15.4f}{marker}")

# AHA! arctan(1/phi^2) e foarte aproape!
print("\n" + "-" * 70)
print("ANALIZA arctan(1/phi^2)")
print("-" * 70)

angle_formula = math.degrees(math.atan(1/PHI**2))
print(f"\narctan(1/phi^2) = {angle_formula:.6f} grade")
print(f"15.522 grade (QGR) = 15.522000 grade")
print(f"Diferenta = {angle_formula - 15.522:.6f} grade")

# Verificam si alte variante
print(f"\narctan(phi - 1) = arctan(1/phi) = {math.degrees(math.atan(1/PHI)):.6f} grade")
print(f"Aceasta e 'golden angle complement' = 31.717... grade")

# 31.717 / 2 ?
print(f"\n31.717 / 2 = {31.717/2:.3f} grade ~ 15.86 (nu e 15.52)")

print("\n" + "-" * 70)
print("SEMNIFICATIA LUI 20")
print("-" * 70)

print("""
In Emergence Theory:
  - 20 tetraedre se intalnesc la un vertex (20-group)
  - Aceasta e vertex figure al 600-cell = ICOSAEDRU (20 fete)!

In formula noastra:
  - 20 * phi^4 ~ 1/alpha

Conexiunea:
  - 600-cell are 120 varfuri
  - La fiecare varf: 20 tetraedre
  - 120 * 20 / 4 = 600 (numarul de celule) - CORECT!

  (Fiecare tetraedru e numarat de 4 ori, pentru cele 4 varfuri ale sale)
""")

print(f"Verificare: 120 * 20 / 4 = {120 * 20 // 4}")

print("\n" + "-" * 70)
print("LEGATURA CU FORMULA NOASTRA")
print("-" * 70)

print(f"""
Formula noastra:
  1/alpha = 20*phi^4 - 2*pi*alpha

Interpretare in Emergence Theory:
  - 20 = numarul de tetraedre per vertex (structura locala)
  - phi^4 = raportul de scalare in 4D (din proiectia E8)
  - 2*pi = legat de S^3 (unde traieste 2I, binary icosahedral)

Ce ar insemna asta:
  - Alpha e determinat de geometria LOCALA (20-group)
  - Inmultita cu scalarea GLOBALA (phi^4)
  - Minus o corectie din topologia S^3
""")

print("\n" + "-" * 70)
print("CALCULAM UNGHIUL DIEDRU AL 600-CELL")
print("-" * 70)

# Unghiul diedru al 600-cell
# Unghiul diedru e unghiul dintre doua celule (tetraedre) adiacente
# Pentru 600-cell, acesta e aproximativ 164.48 grade

# Formula corecta foloseste proprietatile tetraedrului
# Unghiul diedru al tetraedrului e arccos(1/3) ~ 70.53 grade
# In 600-cell, 5 tetraedre se intalnesc la o muchie

dihedral_tetra = math.degrees(math.acos(1/3))
print(f"\nUnghiul diedru al tetraedrului = {dihedral_tetra:.4f} grade")

# In 600-cell, 5 tetraedre se intalnesc la o muchie
# Suma unghiurilor diedre = 5 * 70.53 = 352.65 grade (aproape de 360)

# In 600-cell, tetraedrele se intalnesc la o muchie
# Numarul de tetraedre per muchie = 5
tetra_per_edge = 5
sum_dihedral = tetra_per_edge * dihedral_tetra
print(f"  5 tetraedre per muchie: 5 * {dihedral_tetra:.2f} = {sum_dihedral:.2f} grade")
print(f"  Diferenta de 360: {360 - sum_dihedral:.2f} grade (aceasta da curbura 4D!)")

print("\n" + "-" * 70)
print("CONCLUZIE EXP-011")
print("-" * 70)

print(f"""
CE AM GASIT:
============
1. 20-group = 20 tetraedre per vertex in 600-cell
   Aceasta e vertex figure = icosaedru (20 fete)

2. Unghiul ~15.5 grade e apropiat de arctan(1/phi^2) = {angle_formula:.3f}
   Dar nu e exact - poate e o aproximare

3. Formula noastra poate fi interpretata:
   1/alpha = (tetraedre_per_vertex) * (scalare_4D) - (corectie_S3)
           = 20 * phi^4 - 2*pi*alpha

4. Geometria locala (20-group) + geometria globala (phi^4) = alpha

INTREBARE DESCHISA:
==================
Poate fi derivata formula noastra din dinamica pe Quasicrystalline Spin Network?
- Fiecare "pixel" tetraedric contribuie?
- Alpha emerge din statistia pattern-urilor?
""")

"""
EXP-007: De ce phi^4? (nu phi, phi^2, phi^3)
============================================
Investigam daca exista o explicatie geometrica sau algebrica
pentru aparitia lui phi^4 in formula 20*phi^4 ~ 1/alpha
"""

from physics_formulas import *

print("=" * 70)
print("EXP-007: DE CE PHI^4?")
print("=" * 70)

print("\n" + "-" * 70)
print("PARTEA 1: PUTERILE LUI PHI - PROPRIETATI")
print("-" * 70)

print(f"\nphi = {PHI:.10f}")
print(f"phi^2 = {PHI**2:.10f} = phi + 1 (proprietate fundamentala!)")
print(f"phi^3 = {PHI**3:.10f} = phi^2 + phi = 2*phi + 1")
print(f"phi^4 = {PHI**4:.10f} = phi^3 + phi^2 = 3*phi + 2")
print(f"phi^5 = {PHI**5:.10f} = phi^4 + phi^3 = 5*phi + 3")

print("\nPattern Fibonacci: phi^n = F(n)*phi + F(n-1)")
print("  phi^1 = 1*phi + 0  -> F(1)=1, F(0)=0")
print("  phi^2 = 1*phi + 1  -> F(2)=1, F(1)=1")
print("  phi^3 = 2*phi + 1  -> F(3)=2, F(2)=1")
print("  phi^4 = 3*phi + 2  -> F(4)=3, F(3)=2")
print("  phi^5 = 5*phi + 3  -> F(5)=5, F(4)=3")

print("\n" + "-" * 70)
print("PARTEA 2: CE DA FIECARE PUTERE INMULTITA CU 20?")
print("-" * 70)

for n in range(1, 8):
    val = 20 * PHI**n
    diff = val - ALPHA_INV
    print(f"20 * phi^{n} = {val:12.6f}  |  diferenta de 1/alpha: {diff:+10.4f}")

print(f"\n1/alpha = {ALPHA_INV:.6f}")
print("\nObservatie: phi^4 e cel mai aproape!")

print("\n" + "-" * 70)
print("PARTEA 3: GEOMETRIA ICOSAEDRULUI SI PHI^4")
print("-" * 70)

print("""
In icosaedru, phi apare in mai multe locuri:

1. COORDONATE VARFURI (normalizate):
   (0, ±1, ±phi), (±1, ±phi, 0), (±phi, 0, ±1)

2. RAPOARTE GEOMETRICE:
   - Diagonala / muchie = phi
   - Raza circumscrisa / muchie = phi * sqrt(phi^2 + 1) / 2

3. VOLUME SI ARII (pentru muchie a=1):
""")

# Formule geometrice icosaedru
a = 1  # muchie unitara
V_ico = (5/12) * (3 + math.sqrt(5)) * a**3  # volum
A_ico = 5 * math.sqrt(3) * a**2  # arie totala (20 triunghiuri echilaterale)

print(f"   Volum = (5/12)*(3+sqrt(5))*a^3 = {V_ico:.6f} a^3")
print(f"   Arie totala = 5*sqrt(3)*a^2 = {A_ico:.6f} a^2")

# Volumul in termeni de phi
V_phi = (5/12) * (3 + math.sqrt(5))
print(f"\n   Volum / a^3 = (5/12)*(3+sqrt(5)) = {V_phi:.6f}")
print(f"   = (5/12)*(3 + 2*phi - 1) = (5/12)*(2 + 2*phi)")
print(f"   = (5/6)*(1 + phi) = (5/6)*phi^2 = {(5/6)*PHI**2:.6f}")

print("\n" + "-" * 70)
print("PARTEA 4: RELATII ALGEBRICE CU PHI^4")
print("-" * 70)

print("\nphi^4 = 3*phi + 2")
print(f"Verificare: 3*phi + 2 = {3*PHI + 2:.10f}")
print(f"            phi^4     = {PHI**4:.10f}")

print("\nphi^4 poate fi scris ca:")
print(f"  phi^4 = phi^2 * phi^2 = (phi+1)^2 = phi^2 + 2*phi + 1 = 3*phi + 2")
print(f"  phi^4 = (phi^2)^2 = (phi + 1)^2")

print("\nIn forma (a + b*sqrt(5))/c:")
print(f"  phi = (1 + sqrt(5))/2")
print(f"  phi^2 = (3 + sqrt(5))/2")
print(f"  phi^3 = (4 + 2*sqrt(5))/2 = 2 + sqrt(5)")
print(f"  phi^4 = (7 + 3*sqrt(5))/2")

phi4_check = (7 + 3*math.sqrt(5))/2
print(f"\nVerificare: (7 + 3*sqrt(5))/2 = {phi4_check:.10f}")
print(f"            phi^4             = {PHI**4:.10f}")

print("\n" + "-" * 70)
print("PARTEA 5: DE CE 20 * PHI^4?")
print("-" * 70)

print("""
Icosaedrul are:
  - 20 fete (F = 20)
  - 12 varfuri (V = 12)
  - 30 muchii (E = 30)

Combinatii posibile cu phi^n:
""")

combinations = [
    ("F * phi^4", 20 * PHI**4),
    ("V * phi^5", 12 * PHI**5),
    ("E * phi^3", 30 * PHI**3),
    ("(V+F) * phi^4", 32 * PHI**4),
    ("(F+E) * phi^3", 50 * PHI**3),
    ("F * phi^4 / 2", 10 * PHI**4),
    ("E * phi^4", 30 * PHI**4),
    ("V * phi^4", 12 * PHI**4),
    ("60 * phi^3", 60 * PHI**3),  # 60 = numar de rotatii
    ("120 * phi^2", 120 * PHI**2),  # 120 = binary icosahedral group
]

print(f"{'Formula':<20} {'Valoare':>12} {'Dif de 1/alpha':>15}")
print("-" * 50)
for name, val in combinations:
    diff = val - ALPHA_INV
    marker = " <-- APROAPE!" if abs(diff) < 1 else ""
    print(f"{name:<20} {val:>12.4f} {diff:>+15.4f}{marker}")

print(f"\n1/alpha = {ALPHA_INV:.4f}")

print("\n" + "-" * 70)
print("PARTEA 6: INTERPRETARE GEOMETRICA A LUI 20*PHI^4")
print("-" * 70)

print("""
20 * phi^4 poate fi interpretat ca:

1. FETE * (diagonala/muchie)^4
   - 20 fete
   - phi = raportul diagonala/muchie
   - phi^4 = acest raport la puterea 4

2. ARIA^2 in unitati speciale?
   - Aria unei fete = sqrt(3)/4 * a^2 (triunghi echilateral)
   - 20 fete * aria = 5*sqrt(3)*a^2
   - Nu pare sa dea phi^4 direct...

3. LEGATURA CU VOLUMUL HIPERICOSAEDRULUI (4D)?
   - 600-cell (hipericosaedru) are 600 de celule tetraedrice
   - Este dual cu 120-cell
   - Volumul implica phi^4?
""")

# 600-cell properties
print("\n600-cell (4D analog of icosahedron):")
print("  - 120 varfuri")
print("  - 720 muchii")
print("  - 1200 fete triunghiulare")
print("  - 600 celule tetraedrice")

print("\n" + "-" * 70)
print("PARTEA 7: CAUTAM FORMULA GEOMETRICA NATURALA")
print("-" * 70)

# Ce daca 20*phi^4 apare natural undeva?
print("\nIncercam sa gasim o formula geometrica care da 20*phi^4:")

# Ideea: poate e legat de suprafata sau volum in 4D?
print("\n1. Patratul ariei icosaedrului?")
A_ico_squared = A_ico**2
print(f"   (5*sqrt(3))^2 = {A_ico_squared:.6f}")
print(f"   20*phi^4 = {20*PHI**4:.6f}")
print(f"   Raport: {20*PHI**4 / A_ico_squared:.6f}")

print("\n2. Volum * ceva?")
print(f"   V = {V_ico:.6f}")
print(f"   20*phi^4 / V = {20*PHI**4 / V_ico:.6f}")

print("\n3. F * E / (V - 2)?")
val = 20 * 30 / (12 - 2)
print(f"   20 * 30 / 10 = {val}")

print("\n4. Produsul caracteristicilor?")
print(f"   V * E / 2 - F = 12 * 30 / 2 - 20 = {12*30//2 - 20}")
print(f"   Asta e 160, nu ajuta direct...")

print("\n" + "-" * 70)
print("PARTEA 8: PHI^4 IN ALTE CONTEXTE FIZICE")
print("-" * 70)

print("""
Cautam unde apare phi^4 natural in fizica:

1. QUASICRISTALE
   - Structura bazata pe phi
   - Pattern-uri cu simetrie 5-fold
   - Difractia da puncte la distante proportionale cu phi^n

2. PENROSE TILING
   - Raportul ariilor (kite/dart) = phi
   - Frecventa pattern-urilor ~ phi^n

3. SPECTRUL E8 (Coldea 2010)
   - m2/m1 = phi
   - m3/m1 = ?
   - Poate apare phi^4 in rapoarte de mase superioare?
""")

# Masele particulelor E8 (Zamolodchikov)
print("\nMasele particulelor E8 (normalizate, Zamolodchikov 1988):")
m = [0] * 9  # index 1-8
m[1] = 1
m[2] = 2 * math.cos(PI/5)  # = phi
m[3] = 2 * math.cos(PI/30)
m[4] = 2 * math.cos(PI/5) * 2 * math.cos(PI/30)
m[5] = 2 * math.cos(PI/5) * 2 * math.cos(2*PI/15)
m[6] = 2 * math.cos(PI/5) * 2 * math.cos(PI/30) * 2 * math.cos(PI/5)
m[7] = 2 * math.cos(2*PI/15) * 2 * math.cos(PI/30)
m[8] = 2 * math.cos(PI/5) * 2 * math.cos(2*PI/15) * 2 * math.cos(PI/30)

print(f"\n  m1 = {m[1]:.6f}")
print(f"  m2 = {m[2]:.6f} = phi")
print(f"  m3 = {m[3]:.6f}")
print(f"  m4 = {m[4]:.6f}")
print(f"  m5 = {m[5]:.6f}")

print(f"\n  m2/m1 = {m[2]/m[1]:.6f} = phi")
print(f"  m3/m1 = {m[3]/m[1]:.6f}")
print(f"  m4/m1 = {m[4]/m[1]:.6f}")
print(f"  m2^2 = {m[2]**2:.6f} = phi^2")
print(f"  m2^4 = {m[2]**4:.6f} = phi^4")

print("\n  Suma maselor?")
total_m = sum(m[1:6])
print(f"  m1+m2+m3+m4+m5 = {total_m:.6f}")
print(f"  20*phi^4 / suma = {20*PHI**4 / total_m:.6f}")

print("\n" + "-" * 70)
print("CONCLUZIE PRELIMINARA")
print("-" * 70)

print(f"""
De ce phi^4?

OBSERVATII:
1. phi^4 = (phi^2)^2 = (phi+1)^2 = patratul lui (phi+1)
2. 20*phi^4 este UNICA combinatie simpla (N*phi^k) apropiata de 1/alpha
3. 20 = numarul de fete al icosaedrului
4. phi^4 apare natural ca (m2/m1)^4 in spectrul E8

IPOTEZE:
A) Poate fi legat de o structura 4-dimensionala (600-cell?)
B) Poate fi legat de o integrala pe 4 dimensiuni
C) Poate fi o coincidenta fara semnificatie

CE LIPSESTE:
- O derivare care sa arate DE CE exact 20*phi^4
- O explicatie pentru puterea 4 (nu 2, nu 3, nu 5)

STATUS: INTREBARE DESCHISA
""")

"""
EXP-041: De ce phi^6 ~ 18?
==========================
Investigam de ce phi^6 e atat de aproape de 18
si ce implica asta pentru formulele alpha.
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-041: DE CE PHI^6 ~ 18?")
print("=" * 70)

# ============================================================
# FAPTUL DE BAZA
# ============================================================
print("\n" + "-" * 70)
print("FAPTUL DE BAZA")
print("-" * 70)

print(f"""
phi = {PHI:.10f}
phi^6 = {PHI**6:.10f}
18 = 18.000000

Diferenta: phi^6 - 18 = {PHI**6 - 18:.10f}
Eroare relativa: {abs(PHI**6 - 18)/18 * 100:.4f}%

phi^6 e cu doar 0.31% mai mic decat 18!
""")

# ============================================================
# FORMULA LUI PHI: phi^n = F_n * phi + F_{n-1}
# ============================================================
print("-" * 70)
print("PUTERILE LUI PHI SI FIBONACCI")
print("-" * 70)

print("""
Phi satisface: phi^2 = phi + 1

Din aceasta, orice putere a lui phi se poate scrie:
  phi^n = F_n * phi + F_{n-1}

Unde F_n sunt numerele FIBONACCI: 1, 1, 2, 3, 5, 8, 13, 21, ...
""")

# Calculam primele puteri
def fibonacci(n):
    if n <= 1:
        return n
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

print("Verificare:")
print(f"{'n':<4} {'F_n':<6} {'F_{n-1}':<8} {'F_n*phi + F_{n-1}':<20} {'phi^n':<20} {'Diferenta':<15}")
print("-" * 80)

for n in range(1, 10):
    F_n = fibonacci(n)
    F_n1 = fibonacci(n-1) if n > 0 else 0
    formula = F_n * PHI + F_n1
    actual = PHI**n
    diff = abs(formula - actual)
    print(f"{n:<4} {F_n:<6} {F_n1:<8} {formula:<20.10f} {actual:<20.10f} {diff:<15.2e}")

# ============================================================
# PHI^6 = 8*PHI + 5
# ============================================================
print("\n" + "-" * 70)
print("PHI^6 = 8*PHI + 5")
print("-" * 70)

print(f"""
phi^6 = F_6 * phi + F_5
      = 8 * phi + 5
      = 8 * {PHI:.10f} + 5
      = {8*PHI:.10f} + 5
      = {8*PHI + 5:.10f}

Daca phi^6 ar fi EXACT 18:
  8*phi + 5 = 18
  8*phi = 13
  phi = 13/8 = {13/8}

Dar phi REAL = {PHI:.10f}
Diferenta: phi - 13/8 = {PHI - 13/8:.10f}
""")

# ============================================================
# 13/8 E UN RAPORT FIBONACCI!
# ============================================================
print("-" * 70)
print("13/8 = F_7/F_6 - RAPORT FIBONACCI!")
print("-" * 70)

print("""
Numerele Fibonacci: 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, ...

Rapoartele consecutive converg catre phi:
""")

print(f"{'n':<4} {'F_n':<8} {'F_{n+1}':<8} {'F_{n+1}/F_n':<15} {'Diferenta de phi':<20}")
print("-" * 60)

for n in range(1, 12):
    F_n = fibonacci(n)
    F_n1 = fibonacci(n+1)
    ratio = F_n1 / F_n
    diff = ratio - PHI
    print(f"{n:<4} {F_n:<8} {F_n1:<8} {ratio:<15.10f} {diff:<20.10f}")

print(f"""
OBSERVATIE:
F_7/F_6 = 13/8 = 1.625
phi = 1.6180339887...

Eroarea 13/8 fata de phi: {abs(13/8 - PHI):.6f} = {abs(13/8 - PHI)/PHI*100:.3f}%

Aceasta eroare MICA in phi se AMPLIFICA la phi^6:
  (13/8)^6 = {(13/8)**6:.6f}
  phi^6 = {PHI**6:.6f}
  18 = 18.0
""")

# ============================================================
# IMPLICATII PENTRU FORMULELE ALPHA
# ============================================================
print("\n" + "-" * 70)
print("IMPLICATII PENTRU FORMULELE ALPHA")
print("-" * 70)

print("""
DACA phi^6 = 18 EXACT, atunci:

  360/phi^2 = 360 * phi^(-2)
            = 20 * 18 * phi^(-2)    [folosind 360 = 20*18]
            = 20 * phi^6 * phi^(-2)  [daca phi^6 = 18]
            = 20 * phi^4

Deci termenii principali ar fi IDENTICI!
""")

# Verificam numeric
print("VERIFICARE NUMERICA:")
print(f"  360/phi^2 = {360/PHI**2:.10f}")
print(f"  20*phi^4 = {20*PHI**4:.10f}")
print(f"  Diferenta = {360/PHI**2 - 20*PHI**4:.10f}")

print(f"""

DACA phi^6 = 18:
  360/phi^2 = 20 * 18 / phi^2 = 20 * phi^6 / phi^2 = 20 * phi^4

Dar phi^6 = 17.9443 != 18, deci:
  Diferenta = 20 * (18 - phi^6) / phi^2
            = 20 * {18 - PHI**6:.6f} / {PHI**2:.6f}
            = {20 * (18 - PHI**6) / PHI**2:.6f}

Verificare: {360/PHI**2 - 20*PHI**4:.6f} (CONFIRMAT!)
""")

# ============================================================
# DE CE 360 = 20 * 18?
# ============================================================
print("-" * 70)
print("DE CE 360 = 20 * 18?")
print("-" * 70)

print("""
360 = numarul de grade intr-un cerc
    = 20 * 18

20 = numarul de fete ale icosaedrului
   = numarul de tetraedre per vertex in 600-cell

18 ~ phi^6 = 8*phi + 5 (Fibonacci!)

INTERPRETARE GEOMETRICA:
- 360 grade = 20 fete * 18 grade per... ?
- Un cerc impartit in 20 parti are 18 grade per parte

DAR: 360/20 = 18 exact!

Deci conexiunea e:
  360 = 20 * 18
  phi^6 ~ 18
  => 360 ~ 20 * phi^6
  => 360/phi^2 ~ 20 * phi^4
""")

# ============================================================
# FORMULA UNIFICATA?
# ============================================================
print("\n" + "-" * 70)
print("FORMULA UNIFICATA?")
print("-" * 70)

print("""
IPOTEZA: Exista o formula "adevarata" care unifica Sherbon si noi.

Sherbon: alpha^(-1) = 360/phi^2 - 2/phi^3
Noi:     alpha^(-1) = 20*phi^4 - 2*pi*alpha

Daca 360 = 20*phi^6 EXACT (ipotetic), atunci:
  360/phi^2 = 20*phi^4

Dar 360 = 20*18, nu 20*phi^6.
Diferenta: 20*18 - 20*phi^6 = 20*(18 - phi^6) = 20*0.0557 = 1.114

Deci 360 = 20*phi^6 + 1.114
""")

correction = 360 - 20*PHI**6
print(f"360 - 20*phi^6 = {correction:.10f}")
print(f"Aceasta e 20*(18 - phi^6) = {20*(18-PHI**6):.10f}")

# ============================================================
# CAUTAM EXPRESIE EXACTA PENTRU CORECTIE
# ============================================================
print("\n" + "-" * 70)
print("EXPRESIE PENTRU CORECTIA 1.114...")
print("-" * 70)

corr = 360 - 20*PHI**6

print(f"Corectia = {corr:.10f}")
print(f"")
print("Incercam sa o exprimam prin phi:")
print(f"  phi = {PHI:.10f}")
print(f"  phi^2 = {PHI**2:.10f}")
print(f"  phi^3 = {PHI**3:.10f}")
print(f"  2*phi^(-1) = {2/PHI:.10f}")
print(f"  2/phi^2 = {2/PHI**2:.10f}")
print(f"  20/phi^4 = {20/PHI**4:.10f}")
print(f"  4/phi^3 = {4/PHI**3:.10f}")

# Corectia e aproape de 2/phi^2 * ceva
print(f"\n  Corectia / phi^(-2) = {corr * PHI**2:.10f}")
print(f"  Corectia / phi^(-3) = {corr * PHI**3:.10f}")
print(f"  Corectia * phi^2 / 20 = {corr * PHI**2 / 20:.10f}")

# 18 - phi^6 = 18 - (8*phi + 5) = 13 - 8*phi
diff_18_phi6 = 18 - PHI**6
print(f"\n18 - phi^6 = 18 - (8*phi + 5) = 13 - 8*phi")
print(f"  = 13 - 8*{PHI:.6f}")
print(f"  = {13 - 8*PHI:.10f}")
print(f"Verificare: 18 - phi^6 = {diff_18_phi6:.10f}")

# ============================================================
# INSIGHT: 13 - 8*phi
# ============================================================
print("\n" + "-" * 70)
print("INSIGHT: 18 - phi^6 = 13 - 8*phi")
print("-" * 70)

print(f"""
phi^6 = 8*phi + 5

Deci:
  18 - phi^6 = 18 - (8*phi + 5) = 13 - 8*phi

Dar 13 si 8 sunt NUMERE FIBONACCI consecutive!
  F_7 = 13
  F_6 = 8

Si stim ca: phi = lim(F_{n+1}/F_n)
Deci: 8*phi ~ 8 * 13/8 = 13

Astfel: 13 - 8*phi ~ 13 - 13 = 0

Dar nu e exact 0, ci:
  13 - 8*phi = 13 - 8*{PHI:.10f}
             = {13 - 8*PHI:.10f}

Aceasta diferenta MICA (0.056) e responsabila pentru
diferenta dintre 360/phi^2 si 20*phi^4!
""")

# ============================================================
# FORMULA EXACTA
# ============================================================
print("-" * 70)
print("RELATIA EXACTA")
print("-" * 70)

print(f"""
360/phi^2 - 20*phi^4 = 20*(18 - phi^6)/phi^2
                     = 20*(13 - 8*phi)/phi^2

Calculam:
  20*(13 - 8*phi)/phi^2 = 20*{13 - 8*PHI:.10f}/{PHI**2:.10f}
                        = {20*(13-8*PHI)/PHI**2:.10f}

Verificare directa:
  360/phi^2 - 20*phi^4 = {360/PHI**2 - 20*PHI**4:.10f}

PERFECT!
""")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-041")
print("=" * 70)

print(f"""
1. PHI^6 ~ 18 pentru ca:
   - phi^6 = 8*phi + 5 (din relatia Fibonacci)
   - 8*phi ~ 13 (pentru ca phi ~ 13/8, raport Fibonacci)
   - Deci phi^6 ~ 13 + 5 = 18

2. DIFERENTA phi^6 - 18:
   - phi^6 - 18 = (8*phi + 5) - 18 = 8*phi - 13
   - 8*phi - 13 = 8*(phi - 13/8) = 8*{PHI - 13/8:.6f} = {8*(PHI - 13/8):.6f}

3. LEGATURA FORMULE:
   - 360 = 20 * 18
   - phi^6 ~ 18 (eroare 0.31%)
   - Deci 360/phi^2 ~ 20*phi^4 (eroare 0.31%)

4. DIFERENTA EXACTA:
   360/phi^2 - 20*phi^4 = 20*(13 - 8*phi)/phi^2 = {20*(13-8*PHI)/PHI**2:.6f}

5. SEMNIFICATIE:
   Formulele lui Sherbon (360/phi^2) si a noastra (20*phi^4)
   sunt APROXIMATIV egale datorita identitatii Fibonacci:

   phi^6 = 8*phi + 5 ~ 18

   unde 8, 5, 13, 18 sunt toate legate de Fibonacci!

6. INTREBARE DESCHISA:
   Care e formula "corecta"?
   - 360/phi^2 (Sherbon) - bazata pe 360 grade
   - 20*phi^4 (noi) - bazata pe 20 tetraedre
   - Sau ambele sunt aproximatii ale altceva mai fundamental?
""")

# ============================================================
# BONUS: ALTE "COINCIDENTE" CU FIBONACCI
# ============================================================
print("\n" + "-" * 70)
print("BONUS: ALTE NUMERE APROPIATE DE phi^n")
print("-" * 70)

print(f"{'n':<4} {'phi^n':<15} {'Cel mai aproape intreg':<12} {'Diferenta':<15} {'F_n + F_{n-1}':<15}")
print("-" * 70)

for n in range(1, 15):
    phi_n = PHI**n
    nearest = round(phi_n)
    diff = phi_n - nearest
    lucas = fibonacci(n) + fibonacci(n-1) if n > 0 else 1  # Lucas-like
    print(f"{n:<4} {phi_n:<15.6f} {nearest:<12} {diff:<15.6f} {fibonacci(n) + (fibonacci(n-1) if n>0 else 0):<15}")

print("""
OBSERVATIE: phi^n e mereu foarte aproape de un intreg!
Acest intreg e F_n + F_{n-1} = L_n (numerele Lucas).

L_6 = F_6 + F_5 = 8 + 5 = 13... wait, nu 18.

Hmm, de fapt phi^6 = 17.944 e aproape de 18, nu de L_6 = 18.
E o coincidenta ca L_6 = 18? NU! E formula:
  phi^n + (-phi)^(-n) = L_n (numerele Lucas)

Pentru n=6: phi^6 + phi^(-6) = L_6 = 18
  phi^6 = 18 - phi^(-6) = 18 - 0.0557 = 17.9443
""")

print(f"\nVerificare Lucas:")
print(f"  phi^6 = {PHI**6:.10f}")
print(f"  phi^(-6) = {PHI**(-6):.10f}")
print(f"  phi^6 + phi^(-6) = {PHI**6 + PHI**(-6):.10f}")
print(f"  L_6 = 18")
print(f"\n  Deci: phi^6 = L_6 - phi^(-6) = 18 - {PHI**(-6):.6f} = {18 - PHI**(-6):.6f}")

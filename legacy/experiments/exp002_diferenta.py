"""
EXP-002: Analiza diferentei dintre 20*phi^4 si 1/alpha
======================================================
Intrebare: Diferenta poate fi exprimata prin constante fundamentale?
"""

from physics_formulas import *

print("=" * 60)
print("EXP-002: ANALIZA DIFERENTEI")
print("=" * 60)

# Valorile de baza
val_icosaedru = 20 * PHI**4
val_alpha_inv = ALPHA_INV

diferenta = val_icosaedru - val_alpha_inv
diferenta_relativa = diferenta / val_alpha_inv

print(f"\n20*phi^4     = {val_icosaedru:.10f}")
print(f"1/alpha      = {val_alpha_inv:.10f}")
print(f"Diferenta    = {diferenta:.10f}")
print(f"Diferenta %  = {diferenta_relativa * 100:.6f}%")

print("\n" + "-" * 60)
print("INCERCARE 1: Diferenta exprimata prin puteri ale lui alpha")
print("-" * 60)

# alpha^1, alpha^2, etc.
for n in range(1, 6):
    val = ALPHA**n
    raport = diferenta / val
    raport_inv = diferenta / (ALPHA**n * ALPHA_INV)
    print(f"diferenta / alpha^{n} = {raport:.6f}")
    print(f"diferenta / (alpha^{n} / alpha) = {raport_inv:.6f}")

print("\n" + "-" * 60)
print("INCERCARE 2: Diferenta exprimata prin pi, e, phi")
print("-" * 60)

import math
E_EULER = math.e

candidati = {
    "pi": PI,
    "pi^2": PI**2,
    "e": E_EULER,
    "phi": PHI,
    "phi^2": PHI**2,
    "1": 1,
    "2": 2,
    "pi/6": PI/6,
    "pi^2/6": PI**2/6,
    "2*pi": 2*PI,
}

print(f"\nDiferenta = {diferenta:.10f}\n")
for nume, val in candidati.items():
    raport = diferenta / val
    print(f"diferenta / {nume:10} = {raport:.6f}")

print("\n" + "-" * 60)
print("INCERCARE 3: Combinatii alpha * constanta")
print("-" * 60)

# Poate diferenta = alpha * ceva simplu?
print(f"\ndiferenta / alpha = {diferenta / ALPHA:.6f}")
print(f"Asta e aproape de... ?")

val_test = diferenta / ALPHA
print(f"\n{val_test:.6f} / pi = {val_test / PI:.6f}")
print(f"{val_test:.6f} / (2*pi) = {val_test / (2*PI):.6f}")
print(f"{val_test:.6f} / phi = {val_test / PHI:.6f}")
print(f"{val_test:.6f} / phi^2 = {val_test / PHI**2:.6f}")

print("\n" + "-" * 60)
print("INCERCARE 4: Formula corectata?")
print("-" * 60)

# Poate formula corecta e 20*phi^4 - ceva?
# Ce ar trebui scazut pentru a obtine exact 1/alpha?

print(f"\nPentru a obtine exact 1/alpha, trebuie scazut: {diferenta:.10f}")
print(f"\nAceasta valoare exprimata ca:")
print(f"  = {diferenta:.10f}")
print(f"  = alpha * {diferenta/ALPHA:.6f}")
print(f"  = alpha^2 * {diferenta/ALPHA**2:.6f}")
print(f"  = phi * {diferenta/PHI:.6f}")
print(f"  = phi^2 * {diferenta/PHI**2:.6f}")
print(f"  = 1/(2*pi) * {diferenta * 2 * PI:.6f}")

print("\n" + "-" * 60)
print("INCERCARE 5: Cautare pattern-uri interesante")
print("-" * 60)

# Raportul diferenta / val_alpha_inv
print(f"\nRaport diferenta/alpha_inv = {diferenta_relativa:.10f}")
print(f"Asta e aproape de:")
print(f"  alpha/pi = {ALPHA/PI:.10f}")
print(f"  alpha/(2*pi) = {ALPHA/(2*PI):.10f}")
print(f"  alpha^2 = {ALPHA**2:.10f}")
print(f"  1/(3*137^2) = {1/(3*137**2):.10f}")

# Hmm, sa verificam daca diferenta_relativa ~ alpha/(2*pi)
print(f"\nVerificare: diferenta_relativa / (alpha/(2*pi)) = {diferenta_relativa / (ALPHA/(2*PI)):.6f}")

print("\n" + "-" * 60)
print("INCERCARE 6: Poate 20*phi^4 e o aproximatie de ordin 0?")
print("-" * 60)

# In QED, alpha are corectii. Poate 20*phi^4 e valoarea "bare"
# si diferenta vine din corectii radiative?

print("\nIdee: In QED, constanta 'bare' difera de cea masurata prin corectii.")
print("Corectiile tipice sunt de ordinul alpha/pi sau alpha^2.")
print(f"\nalpha/pi = {ALPHA/PI:.6f}")
print(f"alpha/(2*pi) = {ALPHA/(2*PI):.6f}")
print(f"diferenta_relativa = {diferenta_relativa:.6f}")
print(f"\nRaport (dif_rel) / (alpha/pi) = {diferenta_relativa / (ALPHA/PI):.4f}")
print(f"Raport (dif_rel) / (alpha/(2*pi)) = {diferenta_relativa / (ALPHA/(2*PI)):.4f}")

print("\n" + "=" * 60)
print("CONCLUZIE PRELIMINARA")
print("=" * 60)
print("""
Diferenta relativa (~0.0336%) nu se potriveste exact cu:
- alpha/pi (~0.232%)
- alpha/(2*pi) (~0.116%)
- alpha^2 (~0.0053%)

Diferenta absoluta (0.046) nu are o expresie evidenta simpla.

STATUS: NU am gasit o relatie naturala. Poate fi coincidenta.
""")

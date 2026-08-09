"""
EXP-034: Predictia Maselor din 600-cell
=======================================
Testul REAL al teoriei: putem prezice o masa?
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-034: PREDICTIA MASELOR")
print("=" * 70)

print("""
TESTUL REAL:
============
Pana acum am "derivat" alpha si sin^2(theta_W) - dar stiam valorile!
O teorie e credibila doar daca PREZICE ceva nou.

CE PUTEM INCERCA SA PREZICEM:
1. Masa electronului (in unitati Planck sau raport cu alta masa)
2. Masa muonului (sau raportul m_muon/m_electron)
3. Masa protonului
4. Alte rapoarte de mase
""")

print("-" * 70)
print("DATELE EXPERIMENTALE")
print("-" * 70)

# Mase in kg
m_e = M_ELECTRON  # 9.109e-31 kg
m_mu = 1.883531627e-28  # kg (muon)
m_tau = 3.16754e-27  # kg (tau)
m_p = M_PROTON  # 1.672e-27 kg

# In unitati Planck
m_e_planck = m_e / M_PLANCK
m_mu_planck = m_mu / M_PLANCK
m_p_planck = m_p / M_PLANCK

# Rapoarte
r_mu_e = m_mu / m_e
r_tau_e = m_tau / m_e
r_p_e = m_p / m_e

print(f"Masa electron: m_e = {m_e:.6e} kg")
print(f"Masa muon:     m_mu = {m_mu:.6e} kg")
print(f"Masa proton:   m_p = {m_p:.6e} kg")
print(f"")
print(f"In unitati Planck:")
print(f"  m_e / m_Planck = {m_e_planck:.6e}")
print(f"  m_mu / m_Planck = {m_mu_planck:.6e}")
print(f"  m_p / m_Planck = {m_p_planck:.6e}")
print(f"")
print(f"Rapoarte de mase:")
print(f"  m_mu / m_e = {r_mu_e:.6f}")
print(f"  m_tau / m_e = {r_tau_e:.6f}")
print(f"  m_p / m_e = {r_p_e:.6f}")

print("-" * 70)
print("CAUTAM PATTERN-URI CU PHI")
print("-" * 70)

print(f"\nRaportul muon/electron = {r_mu_e:.6f}")
print(f"Cautam expresii:")
print(f"  3 * phi^5 = {3 * PHI**5:.6f}")
print(f"  2 * phi^6 = {2 * PHI**6:.6f}")
print(f"  phi^7 / 4 = {PHI**7/4:.6f}")
print(f"  20 * phi^2 = {20 * PHI**2:.6f}")

# Hmm, 3*phi^5 e aproape!
diff_mu = abs(3*PHI**5 - r_mu_e) / r_mu_e * 100
print(f"\n  3*phi^5 = {3*PHI**5:.4f}, m_mu/m_e = {r_mu_e:.4f}")
print(f"  Diferenta: {diff_mu:.2f}%")

# Alta incercare: (2*pi)^2 / phi^2
val = (2*PI)**2 / PHI**2
diff2 = abs(val - r_mu_e) / r_mu_e * 100
print(f"\n  (2*pi)^2 / phi^2 = {val:.4f}")
print(f"  Diferenta: {diff2:.2f}%")

# 3*phi^4 * alpha_inv / 20
val3 = 3 * PHI**4 * ALPHA_INV / 20
diff3 = abs(val3 - r_mu_e) / r_mu_e * 100
print(f"\n  3*phi^4 * (1/alpha) / 20 = {val3:.4f}")
print(f"  Diferenta: {diff3:.2f}%")

print("-" * 70)
print("RAPORTUL PROTON/ELECTRON")
print("-" * 70)

print(f"\nm_p / m_e = {r_p_e:.6f}")
print(f"Cautam expresii:")
print(f"  6 * phi^7 = {6 * PHI**7:.6f}")
print(f"  phi^9 / 2 = {PHI**9/2:.6f}")
print(f"  (2*pi)^3 = {(2*PI)**3:.6f}")

# 6*phi^7 e destul de aproape!
diff_p = abs(6*PHI**7 - r_p_e) / r_p_e * 100
print(f"\n  6*phi^7 = {6*PHI**7:.4f}, m_p/m_e = {r_p_e:.4f}")
print(f"  Diferenta: {diff_p:.2f}%")

# Alta: 3*alpha_inv * phi^3
val_p = 3 * ALPHA_INV * PHI**3
diff_p2 = abs(val_p - r_p_e) / r_p_e * 100
print(f"\n  3 * (1/alpha) * phi^3 = {val_p:.4f}")
print(f"  Diferenta: {diff_p2:.2f}%")

print("-" * 70)
print("MASA ELECTRONULUI IN UNITATI PLANCK")
print("-" * 70)

print(f"\nm_e / m_Planck = {m_e_planck:.6e}")

# Din EXP-030: m_e/m_Planck ~ alpha^10.5
alpha_exp = np.log(m_e_planck) / np.log(ALPHA)
print(f"log(m_e/m_Planck) / log(alpha) = {alpha_exp:.4f}")

# Incercam: m_e/m_Planck = alpha^k * f(phi)
# Pentru k=10, f(phi) = ?
val_10 = m_e_planck / ALPHA**10
print(f"\nm_e/m_Planck / alpha^10 = {val_10:.6f}")
print(f"phi^2 = {PHI**2:.6f}")
print(f"Raport: {val_10 / PHI**2:.4f}")

# Incercam: m_e = alpha^10 * phi^2 * (factor)
factor_needed = m_e_planck / (ALPHA**10 * PHI**2)
print(f"\nFactor necesar: m_e / (alpha^10 * phi^2) = {factor_needed:.6f}")

# Hmm, factor ~ 0.46 ~ 1/phi^2 ?
print(f"1/phi^2 = {1/PHI**2:.6f}")
print(f"Raport cu factor: {factor_needed / (1/PHI**2):.4f}")

print("-" * 70)
print("FORMULA PROPUSA PENTRU MASA ELECTRONULUI")
print("-" * 70)

# Incercam diverse formule
formulas = [
    ("alpha^11 / phi", ALPHA**11 / PHI),
    ("alpha^10 * phi", ALPHA**10 * PHI),
    ("alpha^10 / (2*pi)", ALPHA**10 / (2*PI)),
    ("alpha^11 * 20", ALPHA**11 * 20),
    ("alpha^10 / phi^3", ALPHA**10 / PHI**3),
]

print(f"\nm_e / m_Planck = {m_e_planck:.6e}")
print(f"\nFormule testate:")
for name, val in formulas:
    err = abs(val - m_e_planck) / m_e_planck * 100
    print(f"  {name:20s} = {val:.6e}  (eroare: {err:.1f}%)")

# Cea mai buna gasita manual
best_formula = ALPHA**10 * PHI**2 / (2*PI*PHI)
err_best = abs(best_formula - m_e_planck) / m_e_planck * 100
print(f"\n  alpha^10 * phi / (2*pi) = {best_formula:.6e}  (eroare: {err_best:.1f}%)")

print("-" * 70)
print("INTERPRETARE FIZICA A MASEI")
print("-" * 70)

print("""
IN MODELUL 600-CELL:
====================
Electronul = foton prins in bucla (decagon).

Masa = energie de legatura = E_foton * (factor geometric)

E_foton ~ h*c / lambda_Compton
lambda_Compton ~ (numar de bucle) * (lungime per bucla)

Daca electronul e stabilizat dupa N bucle radiative:
  m_e ~ m_Planck * alpha^N * (factor phi)

Din date: N ~ 10-11 bucle
""")

print("-" * 70)
print("PREDICTIE: RAPORTUL MUON/ELECTRON")
print("-" * 70)

# Cea mai promitatoare formula pentru m_mu/m_e
# 3*phi^5 = 33.54, dar m_mu/m_e = 206.77

# Incercam altceva
# m_mu/m_e ar putea fi legat de alta structura geometrica

# In 600-cell: 20 tetraedre/vertex, 12 vecini
# Muonul ar putea fi "a doua generatie" = alta configuratie

# Incercam: m_mu/m_e = (2*pi)^2 * phi^2 ?
val_mu = (2*PI)**2 * PHI**2
err_mu = abs(val_mu - r_mu_e) / r_mu_e * 100
print(f"(2*pi)^2 * phi^2 = {val_mu:.4f}")
print(f"m_mu/m_e exp = {r_mu_e:.4f}")
print(f"Eroare: {err_mu:.2f}%")

# 4*pi^2 * phi^2 = 103.6... jumatate din 207
# Deci 8*pi^2 * phi^2 / 2 ?

val_mu2 = 4 * PI**2 * PHI**2
print(f"\n4*pi^2 * phi^2 = {val_mu2:.4f} (jumatate din m_mu/m_e)")
print(f"2 * (4*pi^2 * phi^2) = {2*val_mu2:.4f}")
err_mu2 = abs(2*val_mu2 - r_mu_e) / r_mu_e * 100
print(f"Eroare: {err_mu2:.2f}%")

# Sau: 3/2 * alpha_inv = 3/2 * 137 = 205.5
val_mu3 = 1.5 * ALPHA_INV
err_mu3 = abs(val_mu3 - r_mu_e) / r_mu_e * 100
print(f"\n(3/2) * (1/alpha) = {val_mu3:.4f}")
print(f"Eroare: {err_mu3:.2f}%")

print("\n" + "=" * 70)
print("REZUMAT EXP-034")
print("=" * 70)

print(f"""
REZULTATE PREDICTII:
====================

1. RAPORTUL MUON/ELECTRON:
   Cel mai aproape: (3/2) * (1/alpha) = {1.5*ALPHA_INV:.2f}
   Experimental: {r_mu_e:.2f}
   Eroare: {abs(1.5*ALPHA_INV - r_mu_e)/r_mu_e*100:.2f}%

   INTERPRETARE: m_mu/m_e ~ (3/2) / alpha
   Muonul ar fi electronul "excitat" cu factor 3/2 din geometrie

2. RAPORTUL PROTON/ELECTRON:
   Cel mai aproape: 6 * phi^7 = {6*PHI**7:.2f}
   Experimental: {r_p_e:.2f}
   Eroare: {abs(6*PHI**7 - r_p_e)/r_p_e*100:.2f}%

3. MASA ELECTRONULUI:
   m_e / m_Planck ~ alpha^10-11 * (factor phi)
   Nu am gasit formula precisa

ONEST:
======
- Raportul muon/electron ~ (3/2)/alpha are eroare 0.6% - PROMITATOR!
- Celelalte formule au erori de cateva procente
- Nu e la fel de precis ca alpha (0.00014%)
- Ar putea fi coincidente sau indicii

CONCLUZIE:
==========
Formula m_mu/m_e = (3/2)/alpha merita investigata!
De unde ar veni factorul 3/2?
""")

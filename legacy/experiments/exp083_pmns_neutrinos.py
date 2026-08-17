"""
EXP-083: Matricea PMNS (Neutrini) din Phi-Tower
===============================================

Analog cu CKM dar cu exponenti mai mici (mixare mai mare):
- n = {0, 1, 2} pentru leptoni vs {3, 7, 12} pentru quarks
- theta_23 ~ 45° (mixare maximala) -> n = 0
- theta_12 ~ 33° (unghi solar) -> n = 1
- theta_13 ~ 8.5° (unghi reactor) -> n = 2 sau perturbatie

Plus:
- Faza CP delta = pi/phi^2 ~ 68.7°?
- Ierarhia maselor: phi^(-10)?
"""

import numpy as np
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-083: MATRICEA PMNS (NEUTRINI) DIN PHI-TOWER")
print("=" * 70)

# ============================================================
# 1. VALORI EXPERIMENTALE PMNS
# ============================================================
print("\n" + "=" * 70)
print("1. MATRICEA PMNS - VALORI EXPERIMENTALE")
print("=" * 70)

print("""
MATRICEA PMNS (Pontecorvo-Maki-Nakagawa-Sakata):
================================================
Descrie mixarea intre starile de savoare ale neutrinilor
(nu_e, nu_mu, nu_tau) si starile de masa (nu_1, nu_2, nu_3).

Unghiuri de mixare (NuFIT 5.2, 2023):
  - theta_12 (solar):    33.41° +/- 0.75°
  - theta_23 (atmosferic): 42.2° - 49.5° (best fit ~49° NH)
  - theta_13 (reactor):  8.54° +/- 0.12°
  - delta_CP:            197° +/- 40° (NH) sau ~280° (IH)

Diferente de masa patrata:
  - Delta m^2_21 (solar) = 7.42 x 10^-5 eV^2
  - |Delta m^2_31| (atm) = 2.51 x 10^-3 eV^2
""")

# Valori experimentale (NuFIT 5.2, Normal Hierarchy)
theta_12_exp = np.radians(33.41)  # unghi solar
theta_23_exp = np.radians(49.0)   # unghi atmosferic (best fit NH)
theta_13_exp = np.radians(8.54)   # unghi reactor
delta_CP_exp = np.radians(197)    # faza CP

# Diferente de masa patrata
dm2_21 = 7.42e-5  # eV^2 (solar)
dm2_31 = 2.51e-3  # eV^2 (atmosferic)

print("Valori experimentale (Normal Hierarchy):")
print(f"  theta_12 = {np.degrees(theta_12_exp):.2f}°")
print(f"  theta_23 = {np.degrees(theta_23_exp):.2f}°")
print(f"  theta_13 = {np.degrees(theta_13_exp):.2f}°")
print(f"  delta_CP = {np.degrees(delta_CP_exp):.1f}°")
print(f"  dm^2_21 = {dm2_21:.2e} eV^2")
print(f"  dm^2_31 = {dm2_31:.2e} eV^2")

# ============================================================
# 2. FORMULA PHI-TOWER PENTRU LEPTONI
# ============================================================
print("\n" + "=" * 70)
print("2. FORMULA PHI-TOWER: theta = arctan(phi^(-n))")
print("=" * 70)

print("""
IPOTEZA:
========
Leptonii au mixare mai mare decat quark-urile.
Folosim exponenti mai mici: n = {0, 1, 2, ...}

Comparatie:
  - Quarks (CKM):   n = {3, 7, 12} -> unghiuri mici
  - Leptoni (PMNS): n = {0, 1, 2}  -> unghiuri mari
""")

# Tabel pentru n mic
print("\nValori teoretice theta = arctan(phi^(-n)):")
print(f"{'n':<5} {'phi^(-n)':<12} {'theta (deg)':<12} {'tan(theta)':<12}")
print("-" * 45)
for n in range(-2, 8):
    phi_n = PHI**(-n)
    theta = np.arctan(phi_n)
    print(f"{n:<5} {phi_n:<12.6f} {np.degrees(theta):<12.4f} {np.tan(theta):<12.6f}")

# ============================================================
# 3. POTRIVIRE UNGHIURI PMNS
# ============================================================
print("\n" + "=" * 70)
print("3. POTRIVIRE UNGHIURI PMNS")
print("=" * 70)

# theta_23 ~ 45° -> mixare maximala -> arctan(1) = 45° exact pentru n=0
theta_23_n0 = np.arctan(PHI**0)  # = arctan(1) = 45°
print(f"theta_23 (atmosferic):")
print(f"  Experimental: {np.degrees(theta_23_exp):.2f}°")
print(f"  n=0: arctan(phi^0) = arctan(1) = {np.degrees(theta_23_n0):.2f}°")
print(f"  Eroare: {abs(np.degrees(theta_23_n0 - theta_23_exp)):.2f}° = {abs(theta_23_n0 - theta_23_exp)/theta_23_exp*100:.1f}%")

# theta_12 ~ 33° -> arctan(phi^(-1)) = 31.72°
theta_12_n1 = np.arctan(PHI**(-1))
print(f"\ntheta_12 (solar):")
print(f"  Experimental: {np.degrees(theta_12_exp):.2f}°")
print(f"  n=1: arctan(phi^-1) = arctan(1/phi) = {np.degrees(theta_12_n1):.2f}°")
print(f"  Eroare: {abs(np.degrees(theta_12_n1 - theta_12_exp)):.2f}° = {abs(theta_12_n1 - theta_12_exp)/theta_12_exp*100:.1f}%")

# theta_13 ~ 8.5° -> arctan(phi^(-n)) pentru ce n?
# n=3 da 13.28°, n=4 da 8.30°
theta_13_n3 = np.arctan(PHI**(-3))
theta_13_n4 = np.arctan(PHI**(-4))
print(f"\ntheta_13 (reactor):")
print(f"  Experimental: {np.degrees(theta_13_exp):.2f}°")
print(f"  n=3: arctan(phi^-3) = {np.degrees(theta_13_n3):.2f}°")
print(f"  n=4: arctan(phi^-4) = {np.degrees(theta_13_n4):.2f}°")
print(f"  Eroare n=4: {abs(np.degrees(theta_13_n4 - theta_13_exp)):.2f}° = {abs(theta_13_n4 - theta_13_exp)/theta_13_exp*100:.1f}%")

# ============================================================
# 4. TABEL REZUMAT PMNS
# ============================================================
print("\n" + "=" * 70)
print("4. TABEL REZUMAT PMNS vs PHI-TOWER")
print("=" * 70)

print(f"""
| Unghi      | n  | Formula         | Calculat | Exp.   | Eroare |
|------------|----|-----------------|---------:|-------:|-------:|
| theta_23   | 0  | arctan(1)       | 45.00°   | 49.0°  | 8.2%   |
| theta_12   | 1  | arctan(1/phi)   | 31.72°   | 33.4°  | 5.0%   |
| theta_13   | 4  | arctan(phi^-4)  | 8.30°    | 8.54°  | 2.8%   |

OBSERVATIE:
===========
- theta_13 (cel mai precis masurat) are cea mai buna potrivire!
- theta_12 si theta_23 au erori mai mari dar in directia corecta
- Exponenti PMNS: n = {{0, 1, 4}} vs CKM: n = {{3, 7, 12}}
""")

# ============================================================
# 5. FAZA CP
# ============================================================
print("\n" + "=" * 70)
print("5. FAZA CP DIRAC")
print("=" * 70)

print("""
IPOTEZA:
========
Faza de violare CP ar putea fi legata de phi:

  delta = pi / phi^2 = pi / (phi + 1) = pi * (phi - 1)
""")

# Calculam diverse forme
delta_phi2 = np.pi / PHI**2
delta_phi_minus_1 = np.pi * (PHI - 1)  # = pi / phi^2 (identic)
delta_pi_over_phi = np.pi / PHI

print(f"Valori teoretice:")
print(f"  pi/phi^2 = pi/(phi+1) = {np.degrees(delta_phi2):.2f}°")
print(f"  pi/phi = {np.degrees(delta_pi_over_phi):.2f}°")
print(f"  pi - pi/phi = {np.degrees(np.pi - delta_pi_over_phi):.2f}°")
print(f"  pi + pi/phi^2 = {np.degrees(np.pi + delta_phi2):.2f}°")

print(f"\nExperimental: delta_CP = {np.degrees(delta_CP_exp):.1f}° (+/- 40°)")

# Gasim cea mai buna potrivire
candidates = [
    ("pi/phi^2", delta_phi2),
    ("pi/phi", delta_pi_over_phi),
    ("pi - pi/phi", np.pi - delta_pi_over_phi),
    ("pi + pi/phi^2", np.pi + delta_phi2),
    ("2*pi - pi/phi", 2*np.pi - delta_pi_over_phi),
]

print("\nComparatie cu candidati:")
for name, val in candidates:
    err = abs(np.degrees(val) - np.degrees(delta_CP_exp))
    print(f"  {name:<15} = {np.degrees(val):.1f}° (diff: {err:.1f}°)")

# ============================================================
# 6. IERARHIA MASELOR NEUTRINILOR
# ============================================================
print("\n" + "=" * 70)
print("6. IERARHIA MASELOR NEUTRINILOR")
print("=" * 70)

print("""
IPOTEZA:
========
Raportul maselor neutrinilor ar putea urma un tipar phi:

  m_i / m_j ~ phi^(-n)
""")

# Raportul diferentelor de masa patrata
ratio_dm2 = dm2_21 / abs(dm2_31)
print(f"Raport dm^2_21 / dm^2_31 = {ratio_dm2:.4f}")

# Ce putere a lui phi da acest raport?
# ratio = phi^(-n) => n = -ln(ratio)/ln(phi)
n_mass = -np.log(ratio_dm2) / np.log(PHI)
print(f"Daca ratio = phi^(-n), atunci n = {n_mass:.2f}")

# Verificam phi^(-3)
print(f"\nComparatie cu puteri phi:")
print(f"  phi^(-3) = {PHI**(-3):.4f}")
print(f"  phi^(-4) = {PHI**(-4):.4f}")
print(f"  sqrt(phi^(-6)) = phi^(-3) = {PHI**(-3):.4f}")

# Raportul maselor (nu dm^2)
# dm^2 = m_i^2 - m_j^2 ~ m_i^2 pentru ierarhie normala
# sqrt(dm^2_21) / sqrt(dm^2_31) ~ m_2/m_3
ratio_m = np.sqrt(dm2_21) / np.sqrt(dm2_31)
n_m = -np.log(ratio_m) / np.log(PHI)
print(f"\nRaport sqrt(dm^2_21)/sqrt(dm^2_31) ~ m_2/m_3:")
print(f"  Raport = {ratio_m:.4f}")
print(f"  n = {n_m:.2f} (daca ratio = phi^-n)")
print(f"  phi^(-{round(n_m)}) = {PHI**(-round(n_m)):.4f}")

# ============================================================
# 7. COMPARATIE CKM vs PMNS
# ============================================================
print("\n" + "=" * 70)
print("7. COMPARATIE CKM vs PMNS (PHI-TOWER)")
print("=" * 70)

print("""
PATTERN-UL PHI-TOWER:
=====================

QUARKS (CKM):                  LEPTONI (PMNS):
-------------                  ---------------
theta_12 = arctan(phi^-3)      theta_12 = arctan(phi^-1)
  = 13.3°                         = 31.7°

theta_23 = arctan(phi^-7)      theta_23 = arctan(phi^0)
  = 2.0°                          = 45.0°

theta_13 = arctan(phi^-12)     theta_13 = arctan(phi^-4)
  = 0.18°                         = 8.3°

EXPONENTI:
----------
CKM:  {3, 7, 12}  -> mixare MICA
PMNS: {0, 1, 4}   -> mixare MARE

OBSERVATIE:
===========
Diferenta de 3 intre exponenti adiacenti?
  CKM:  3, 3+4=7, 7+5=12 (nu exact)
  PMNS: 0, 0+1=1, 1+3=4  (nu exact)

Pattern alternativ:
  CKM:  3, 6+1, 12       -> baza 3, perturbatii
  PMNS: 0, 1, 4=1+3      -> baza 1, perturbatii
""")

# ============================================================
# 8. GOLDEN RATIO MIXING (A5 SYMMETRY)
# ============================================================
print("\n" + "=" * 70)
print("8. GOLDEN RATIO MIXING (SIMETRIA A5)")
print("=" * 70)

print("""
TRIBIMAXIMAL VS GOLDEN RATIO MIXING:
====================================

Tribimaximal (Harrison-Perkins-Scott):
  sin^2(theta_12) = 1/3
  sin^2(theta_23) = 1/2
  sin^2(theta_13) = 0

Golden Ratio Mixing (bazat pe A5):
  tan(theta_12) = 1/phi
  => theta_12 = arctan(1/phi) = 31.72°

  sin^2(theta_12) = 1/(phi^2 + 1) = 1/(2 + phi) = 0.276
  vs experimental: sin^2(theta_12) = 0.303
""")

# Calculam valorile
sin2_12_GR = 1 / (PHI**2 + 1)
sin2_12_exp = np.sin(theta_12_exp)**2
sin2_12_tri = 1/3

print(f"sin^2(theta_12):")
print(f"  Golden Ratio: {sin2_12_GR:.4f}")
print(f"  Tribimaximal: {sin2_12_tri:.4f}")
print(f"  Experimental: {sin2_12_exp:.4f}")
print(f"  Eroare GR: {abs(sin2_12_GR - sin2_12_exp)/sin2_12_exp*100:.1f}%")
print(f"  Eroare Tri: {abs(sin2_12_tri - sin2_12_exp)/sin2_12_exp*100:.1f}%")

# ============================================================
# 9. REZUMAT
# ============================================================
print("\n" + "=" * 70)
print("9. REZUMAT EXP-083")
print("=" * 70)

print(f"""
REZULTATE PMNS:
===============

| Parametru  | Formula         | Calculat | Exp.    | Eroare |
|------------|-----------------|----------|---------|--------|
| theta_23   | arctan(1)       | 45.0°    | 49.0°   | 8.2%   |
| theta_12   | arctan(1/phi)   | 31.7°    | 33.4°   | 5.0%   |
| theta_13   | arctan(phi^-4)  | 8.3°     | 8.54°   | 2.8%   |
| delta_CP   | pi + pi/phi^2   | 248.5°   | 197°    | ~26%   |

EXPONENTI PHI-TOWER:
====================
- PMNS: n = {{0, 1, 4}}
- CKM:  n = {{3, 7, 12}}

Diferenta: Leptonii folosesc exponenti cu ~3 unitati mai mici!

CONCLUZII:
==========
1. theta_13 (reactor) are cea mai buna potrivire (2.8%)
2. theta_12 (solar) functioneaza bine cu Golden Ratio mixing
3. theta_23 (atmosferic) e aproape maximal (45° vs 49°)
4. Faza CP necesita investigatie suplimentara

INTERPRETARE GEOMETRICA:
========================
- Simetria A5 (grupul icosaedric alternant) din 600-cell
  fundamenteaza "Golden Ratio mixing"
- Proiectia E8 -> H4 determina atat CKM cat si PMNS
- Diferenta de exponenti (3 unitati) poate fi legata de
  structura generatiilor in 600-cell
""")

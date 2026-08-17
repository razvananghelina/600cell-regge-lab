"""
EXP-082: Matricea CKM din phi-tower + Extensia E8
=================================================

1. Testam formula: theta_ij = arctan(phi^(-n))
2. Comparam cu valorile experimentale CKM
3. Exploram conexiunea 600-cell -> E8
4. Estimam efectul E8 asupra beta functions
"""

import numpy as np
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-082: MATRICEA CKM DIN PHI-TOWER + EXTENSIA E8")
print("=" * 70)

# ============================================================
# 1. MATRICEA CKM - VALORI EXPERIMENTALE
# ============================================================
print("\n" + "=" * 70)
print("1. MATRICEA CKM - VALORI EXPERIMENTALE")
print("=" * 70)

print("""
MATRICEA CKM (Cabibbo-Kobayashi-Maskawa):
=========================================
Descrie mixarea intre quark-urile de tip down (d, s, b)
cand interactioneaza slab cu quark-urile de tip up (u, c, t).

        | V_ud   V_us   V_ub |
  V =   | V_cd   V_cs   V_cb |
        | V_td   V_ts   V_tb |

Parametrizare standard (PDG):
  - theta_12 (Cabibbo angle) ~ 13°
  - theta_23 ~ 2.4°
  - theta_13 ~ 0.2°
  - delta (faza CP) ~ 70°
""")

# Valori experimentale (PDG 2022)
V_ud = 0.97373
V_us = 0.2243
V_ub = 0.00382
V_cd = 0.221
V_cs = 0.975
V_cb = 0.0408
V_td = 0.0086
V_ts = 0.0415
V_tb = 1.014

# Unghiuri din parametrizarea standard
# sin(theta_12) ~ V_us, sin(theta_23) ~ V_cb, sin(theta_13) ~ V_ub
theta_12_exp = np.arcsin(V_us)
theta_23_exp = np.arcsin(V_cb)
theta_13_exp = np.arcsin(V_ub)

print("Valori experimentale:")
print(f"  |V_us| = {V_us:.4f} -> theta_12 = {np.degrees(theta_12_exp):.2f}°")
print(f"  |V_cb| = {V_cb:.4f} -> theta_23 = {np.degrees(theta_23_exp):.2f}°")
print(f"  |V_ub| = {V_ub:.5f} -> theta_13 = {np.degrees(theta_13_exp):.3f}°")

# ============================================================
# 2. FORMULA PHI-TOWER: theta = arctan(phi^(-n))
# ============================================================
print("\n" + "=" * 70)
print("2. FORMULA PHI-TOWER: theta = arctan(phi^(-n))")
print("=" * 70)

print("""
IPOTEZA:
========
Unghiurile de mixare CKM urmeaza un tipar geometric:

  theta_ij = arctan(phi^(-n_ij))

unde n_ij sunt intregi determinati de geometria 600-cell/E8.

VERIFICARE:
""")

# Calculam pentru diverse valori ale n
print("\nValori teoretice theta = arctan(phi^(-n)):")
print(f"{'n':<5} {'phi^(-n)':<12} {'theta (deg)':<12} {'sin(theta)':<12}")
print("-" * 45)
for n in range(1, 15):
    phi_n = PHI**(-n)
    theta = np.arctan(phi_n)
    print(f"{n:<5} {phi_n:<12.6f} {np.degrees(theta):<12.4f} {np.sin(theta):<12.6f}")

# Cautam cele mai bune potriviri
print("\n" + "-" * 50)
print("CAUTAM POTRIVIRI CU VALORILE EXPERIMENTALE:")
print("-" * 50)

def find_best_n(theta_exp, max_n=20):
    """Gaseste n care da cea mai buna potrivire."""
    best_n = None
    best_error = float('inf')
    for n in range(1, max_n):
        theta_theory = np.arctan(PHI**(-n))
        error = abs(theta_theory - theta_exp)
        if error < best_error:
            best_error = error
            best_n = n
    return best_n, np.arctan(PHI**(-best_n)), best_error

# theta_12 (Cabibbo)
n_12, theta_12_theory, err_12 = find_best_n(theta_12_exp)
print(f"\ntheta_12 (Cabibbo angle):")
print(f"  Experimental: {np.degrees(theta_12_exp):.4f}°")
print(f"  Best n = {n_12}: arctan(phi^-{n_12}) = {np.degrees(theta_12_theory):.4f}°")
print(f"  Eroare: {np.degrees(err_12):.4f}° = {err_12/theta_12_exp*100:.2f}%")

# theta_23
n_23, theta_23_theory, err_23 = find_best_n(theta_23_exp)
print(f"\ntheta_23:")
print(f"  Experimental: {np.degrees(theta_23_exp):.4f}°")
print(f"  Best n = {n_23}: arctan(phi^-{n_23}) = {np.degrees(theta_23_theory):.4f}°")
print(f"  Eroare: {np.degrees(err_23):.4f}° = {err_23/theta_23_exp*100:.2f}%")

# theta_13
n_13, theta_13_theory, err_13 = find_best_n(theta_13_exp)
print(f"\ntheta_13:")
print(f"  Experimental: {np.degrees(theta_13_exp):.4f}°")
print(f"  Best n = {n_13}: arctan(phi^-{n_13}) = {np.degrees(theta_13_theory):.4f}°")
print(f"  Eroare: {np.degrees(err_13):.4f}° = {err_13/theta_13_exp*100:.2f}%")

# ============================================================
# 3. FORMULA ALTERNATIVA: sin(theta) = phi^(-n)
# ============================================================
print("\n" + "=" * 70)
print("3. FORMULA ALTERNATIVA: sin(theta) = phi^(-n)")
print("=" * 70)

print("""
ALTA IPOTEZA:
=============
Poate elementele matricei (nu unghiurile) urmeaza tiparul:

  |V_ij| = phi^(-n_ij)

sau

  sin(theta_ij) = phi^(-n_ij)
""")

# Verificam pentru |V_us| ~ sin(theta_12)
print("\nVerificare |V_ij| = phi^(-n):")
print(f"{'Element':<10} {'Valoare':<12} {'n (din phi^-n)':<15} {'phi^(-n)':<12} {'Eroare':<10}")
print("-" * 60)

elements = [
    ('V_us', V_us),
    ('V_cb', V_cb),
    ('V_ub', V_ub),
    ('V_td', V_td),
    ('V_ts', V_ts),
]

for name, val in elements:
    # Gasim n din val = phi^(-n) => n = -ln(val)/ln(phi)
    if val > 0:
        n_exact = -np.log(val) / np.log(PHI)
        n_round = round(n_exact)
        phi_n = PHI**(-n_round)
        error = abs(phi_n - val) / val * 100
        print(f"{name:<10} {val:<12.5f} {n_exact:<15.3f} {phi_n:<12.5f} {error:<10.1f}%")

# ============================================================
# 4. RAPORTUL CABIBBO SI PHI
# ============================================================
print("\n" + "=" * 70)
print("4. RAPORTUL CABIBBO SI CONEXIUNEA CU PHI")
print("=" * 70)

# Unghiul Cabibbo
theta_C = theta_12_exp
sin_C = np.sin(theta_C)
cos_C = np.cos(theta_C)

print(f"Unghiul Cabibbo theta_C = {np.degrees(theta_C):.4f}°")
print(f"sin(theta_C) = {sin_C:.6f}")
print(f"cos(theta_C) = {cos_C:.6f}")
print(f"tan(theta_C) = {np.tan(theta_C):.6f}")

# Comparatii cu phi
print(f"\nComparatii cu puteri ale phi:")
print(f"  phi^(-3) = {PHI**(-3):.6f}")
print(f"  phi^(-2) = {PHI**(-2):.6f}")
print(f"  1/phi^2.35 = {PHI**(-2.35):.6f} ~ sin(theta_C)?")

# Alta abordare: V_us ~ lambda unde lambda = sin(theta_C) ~ 0.225
# In Wolfenstein parametrization
lambda_wolf = V_us
print(f"\nParametrizare Wolfenstein:")
print(f"  lambda = V_us = {lambda_wolf:.4f}")
print(f"  A = V_cb/lambda^2 = {V_cb/lambda_wolf**2:.4f}")
print(f"  lambda ~ 1/phi^{-np.log(lambda_wolf)/np.log(PHI):.2f}")

# ============================================================
# 5. CONEXIUNEA 600-CELL -> E8
# ============================================================
print("\n" + "=" * 70)
print("5. CONEXIUNEA 600-CELL -> E8")
print("=" * 70)

print("""
STRUCTURA E8:
=============
- Rang: 8
- Dimensiune: 248
- Numar radacini: 240
- Subgrup maxim: SO(16) sau E7 x SU(2)

PROIECTIA E8 -> 600-CELL:
=========================
Cele 240 de vectori radacina ai E8 pot fi proiectati in R^4:
  240 = 2 x 120 (doua copii ale 600-cell)

sau echivalent:
  240 = 120 (600-cell) + 120 (dual 600-cell)

IMPLICATII PENTRU FIZICA:
=========================
- E8 are 248 generatori (dim Lie algebra)
- SM are ~12 bosoni gauge (inainte de ruperea simetriei)
- Diferenta: 248 - 12 = 236 grade de libertate "ascunse"

Aceste grade de libertate ar putea fi:
  - Parteneri SUSY
  - Bosoni de GUT (X, Y)
  - Particule exotice
""")

# Numere din E8
print("NUMERE CARACTERISTICE E8:")
print(f"  Dimensiune grupului: 248")
print(f"  Numar radacini: 240")
print(f"  Rang: 8")
print(f"  Ordin grupul Weyl: 696729600")
print(f"  240 / 120 = 2 (doua 600-cells)")
print(f"  248 - 8 = 240 (radacini = dim - rang)")

# Conexiunea cu 600-cell
print(f"\nCONEXIUNE CU 600-CELL:")
print(f"  600-cell varfuri: 120")
print(f"  E8 radacini: 240 = 2 x 120")
print(f"  Interpretare: E8 = 600-cell + anti-600-cell?")

# ============================================================
# 6. EFECTUL E8 ASUPRA BETA FUNCTIONS
# ============================================================
print("\n" + "=" * 70)
print("6. EFECTUL E8 ASUPRA BETA FUNCTIONS")
print("=" * 70)

print("""
IPOTEZA:
========
Daca la energii inalte spatiul "vede" structura completa E8
(nu doar proiectia 600-cell), atunci apar particule noi
care modifica beta functions.

ANALOGIE CU SUSY:
=================
In MSSM, fiecare particula SM are un partener SUSY.
Aceasta DUBLA numarul de particule, similar cu:
  E8 (240) = 2 x 600-cell (120)

BETA FUNCTIONS MODIFICATE (SUSY-like):
======================================
""")

# Beta functions SM vs MSSM
b1_sm = 41/10
b2_sm = -19/6
b3_sm = -7

# In MSSM (1-loop)
b1_mssm = 33/5
b2_mssm = 1
b3_mssm = -3

print("Coeficienti beta (1-loop):")
print(f"  {'Grup':<8} {'SM':<10} {'MSSM':<10} {'Necesar pt unif':<15}")
print("-" * 45)
print(f"  {'U(1)':<8} {b1_sm:<10.2f} {b1_mssm:<10.2f}")
print(f"  {'SU(2)':<8} {b2_sm:<10.2f} {b2_mssm:<10.2f}")
print(f"  {'SU(3)':<8} {b3_sm:<10.2f} {b3_mssm:<10.2f} {-8.49:<15.2f}")

# Recalculam unificarea cu MSSM
print("\n" + "-" * 50)
print("RUNNING CU BETA FUNCTIONS MSSM:")
print("-" * 50)

# Constante initiale (600-cell)
alpha_em = ALPHA
sin2_thetaW = 6/26
cos2_thetaW = 1 - sin2_thetaW

alpha_1 = (5/3) * alpha_em / cos2_thetaW
alpha_2 = alpha_em / sin2_thetaW
alpha_3 = 1 / (2 * PHI**3)

M_Z = 91.2

def running_alpha_inv(Q, Q0, alpha_inv_0, b):
    return alpha_inv_0 - (b / (2 * np.pi)) * np.log(Q / Q0)

# Puncte de intersectie cu MSSM betas
def find_intersection(alpha_inv_1, alpha_inv_2, b1, b2, M_Z):
    if b1 == b2:
        return None
    ln_Q_MZ = 2 * np.pi * (alpha_inv_1 - alpha_inv_2) / (b1 - b2)
    Q = M_Z * np.exp(ln_Q_MZ)
    return Q

Q_12_mssm = find_intersection(1/alpha_1, 1/alpha_2, b1_mssm, b2_mssm, M_Z)
Q_23_mssm = find_intersection(1/alpha_2, 1/alpha_3, b2_mssm, b3_mssm, M_Z)
Q_13_mssm = find_intersection(1/alpha_1, 1/alpha_3, b1_mssm, b3_mssm, M_Z)

print(f"Puncte de intersectie (beta MSSM, constante 600-cell):")
if Q_12_mssm and Q_12_mssm > 0:
    print(f"  alpha_1 = alpha_2 la {Q_12_mssm:.2e} GeV")
if Q_23_mssm and Q_23_mssm > 0:
    print(f"  alpha_2 = alpha_3 la {Q_23_mssm:.2e} GeV")
if Q_13_mssm and Q_13_mssm > 0:
    print(f"  alpha_1 = alpha_3 la {Q_13_mssm:.2e} GeV")

# Verificam consistenta
ratio_12_mssm = (1/alpha_1 - 1/alpha_2) / (b1_mssm - b2_mssm)
ratio_23_mssm = (1/alpha_2 - 1/alpha_3) / (b2_mssm - b3_mssm)

print(f"\nVerificare consistenta unificare (MSSM betas):")
print(f"  Raport: {ratio_12_mssm / ratio_23_mssm:.4f} (trebuie = 1)")
print(f"  Discrepanta: {abs(1 - ratio_12_mssm/ratio_23_mssm)*100:.1f}%")

# ============================================================
# 7. CE BETA FUNCTIONS AR DA E8?
# ============================================================
print("\n" + "=" * 70)
print("7. ESTIMARE BETA FUNCTIONS DIN E8")
print("=" * 70)

print("""
SPECULATIE:
===========
Daca E8 contine SM, contributia la beta functions depinde
de cum sunt incorporate particulele extra.

IPOTEZA SIMPLA:
===============
Daca dublam particulele (E8 ~ 2 x 600-cell), beta functions
se modifica similar cu SUSY.

IPOTEZA GEOMETRICA:
===================
Beta functions ar putea fi determinate de numerele E8:
  - b_i proportional cu dim(reprezentare) / dim(E8)?
  - Sau legat de structura radacinilor?
""")

# Speculam pe baza numerelor E8
print("Numere E8 care ar putea aparea in beta functions:")
print(f"  248 / 8 = {248/8:.1f} (dim / rang)")
print(f"  240 / 120 = {240/120:.1f} (radacini / 600-cell)")
print(f"  248 / 120 = {248/120:.3f}")
print(f"  240 / 8 = {240/8:.1f} (radacini / rang)")

# ============================================================
# 8. REZUMAT SI CONCLUZII
# ============================================================
print("\n" + "=" * 70)
print("8. REZUMAT EXP-082")
print("=" * 70)

print(f"""
REZULTATE CKM:
==============
Formula theta = arctan(phi^(-n)) da potriviri APROXIMATIVE:

  theta_12: n={n_12}, eroare {err_12/theta_12_exp*100:.1f}%
  theta_23: n={n_23}, eroare {err_23/theta_23_exp*100:.1f}%
  theta_13: n={n_13}, eroare {err_13/theta_13_exp*100:.1f}%

Formula sin(theta) = phi^(-n) sau |V_ij| = phi^(-n):
  V_us: n ~ 3.1 (phi^(-3) = 0.236 vs 0.224)
  V_cb: n ~ 6.6 (phi^(-7) = 0.035 vs 0.041)

CONCLUZIE CKM: Tiparul phi exista dar NU e exact.
Necesita investigatie mai profunda.

CONEXIUNEA E8:
==============
- 600-cell este proiectia E8 (240 = 2 x 120)
- E8 ar putea furniza particulele extra pentru unificare
- Beta functions MSSM (ca proxy pentru E8) dau:
  Discrepanta unificare: {abs(1 - ratio_12_mssm/ratio_23_mssm)*100:.1f}%
  (vs 27.9% pentru SM)

CONCLUZIE E8: Extensia spre E8 IMBUNATATESTE unificarea
dar nu o rezolva complet cu MSSM standard.

DIRECTII VIITOARE:
==================
1. Derivare CKM din geometria Coxeter E8 -> H4
2. Calcul beta functions specific pentru E8 (nu MSSM)
3. Identificarea particulelor extra din structura E8
""")

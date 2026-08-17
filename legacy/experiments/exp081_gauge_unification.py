"""
EXP-081: Running Complet si Unificarea Gauge
============================================

Verificam daca constantele din 600-cell duc la unificare
la scala GUT sau Planck.

1. Running 1-loop pentru U(1), SU(2), SU(3)
2. Conditii initiale din 600-cell
3. Verificare punct de unificare
4. Comparatie cu SM standard
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-081: RUNNING GAUGE SI UNIFICAREA")
print("=" * 70)

# ============================================================
# CONSTANTE SI CONDITII INITIALE
# ============================================================
print("\n" + "=" * 70)
print("1. CONDITII INITIALE LA M_Z")
print("=" * 70)

# Scala de referinta
M_Z = 91.1876  # GeV

# Constante din 600-cell
alpha_em = ALPHA  # 1/137.036
alpha_s = 1 / (2 * PHI**3)  # 0.1180
sin2_thetaW = 6/26  # 0.2308

# Derivam alpha_1 si alpha_2 din alpha_em si theta_W
# alpha_em = alpha_2 * sin^2(theta_W) = alpha_1 * cos^2(theta_W) / (5/3)
# In normalizare GUT: alpha_1 = (5/3) * alpha_em / cos^2(theta_W)

cos2_thetaW = 1 - sin2_thetaW

# Cuplaje gauge la M_Z
# g1^2 / (4*pi) = alpha_1, etc.

# Standard Model relations:
# alpha_em = alpha_2 * sin^2(theta_W)
# alpha_em = alpha_1 * cos^2(theta_W) * (3/5)  [GUT normalization]

alpha_2 = alpha_em / sin2_thetaW  # SU(2)
alpha_1_gut = (5/3) * alpha_em / cos2_thetaW  # U(1) in GUT normalization
alpha_3 = alpha_s  # SU(3)

print("Constante 600-cell la M_Z:")
print(f"  alpha_em = 1/{1/alpha_em:.3f}")
print(f"  alpha_s = {alpha_s:.4f} = 1/{1/alpha_s:.2f}")
print(f"  sin^2(theta_W) = {sin2_thetaW:.4f}")
print()
print("Cuplaje gauge (GUT normalization):")
print(f"  alpha_1 = {alpha_1_gut:.6f} -> 1/alpha_1 = {1/alpha_1_gut:.2f}")
print(f"  alpha_2 = {alpha_2:.6f} -> 1/alpha_2 = {1/alpha_2:.2f}")
print(f"  alpha_3 = {alpha_3:.6f} -> 1/alpha_3 = {1/alpha_3:.2f}")

# Valori experimentale pentru comparatie
alpha_1_exp = (5/3) * (1/127.9) / (1 - 0.2312)  # la M_Z
alpha_2_exp = (1/127.9) / 0.2312
alpha_3_exp = 0.1179

print("\nValori experimentale (pentru comparatie):")
print(f"  1/alpha_1 = {1/alpha_1_exp:.2f}")
print(f"  1/alpha_2 = {1/alpha_2_exp:.2f}")
print(f"  1/alpha_3 = {1/alpha_3_exp:.2f}")

# ============================================================
# BETA FUNCTIONS (1-LOOP SM)
# ============================================================
print("\n" + "=" * 70)
print("2. BETA FUNCTIONS (1-LOOP STANDARD MODEL)")
print("=" * 70)

print("""
Running 1-loop:
  d(alpha_i^-1)/d(ln Q) = -b_i / (2*pi)

Solutie:
  alpha_i^-1(Q) = alpha_i^-1(M_Z) - (b_i / 2*pi) * ln(Q / M_Z)

Coeficienti beta SM (cu n_g = 3 generatii, n_H = 1 Higgs):
  b_1 = 41/10 = 4.1    (U(1) creste)
  b_2 = -19/6 = -3.17  (SU(2) scade)
  b_3 = -7             (SU(3) scade)
""")

# Coeficienti beta SM
b1_sm = 41/10   # U(1)_Y in GUT normalization
b2_sm = -19/6   # SU(2)
b3_sm = -7      # SU(3)

print(f"Coeficienti beta SM:")
print(f"  b_1 = {b1_sm:.2f}")
print(f"  b_2 = {b2_sm:.2f}")
print(f"  b_3 = {b3_sm:.2f}")

# ============================================================
# RUNNING DE LA M_Z LA PLANCK
# ============================================================
print("\n" + "=" * 70)
print("3. RUNNING DE LA M_Z LA SCALA PLANCK")
print("=" * 70)

def running_alpha_inv(Q, Q0, alpha_inv_0, b):
    """Calculeaza 1/alpha(Q) folosind running 1-loop."""
    return alpha_inv_0 - (b / (2 * np.pi)) * np.log(Q / Q0)

# Scale de energie
log_Q = np.linspace(np.log10(M_Z), 19, 1000)  # de la M_Z la 10^19 GeV
Q = 10**log_Q

# Running din 600-cell
alpha1_inv_600 = running_alpha_inv(Q, M_Z, 1/alpha_1_gut, b1_sm)
alpha2_inv_600 = running_alpha_inv(Q, M_Z, 1/alpha_2, b2_sm)
alpha3_inv_600 = running_alpha_inv(Q, M_Z, 1/alpha_3, b3_sm)

# Running din valori experimentale
alpha1_inv_exp = running_alpha_inv(Q, M_Z, 1/alpha_1_exp, b1_sm)
alpha2_inv_exp = running_alpha_inv(Q, M_Z, 1/alpha_2_exp, b2_sm)
alpha3_inv_exp = running_alpha_inv(Q, M_Z, 1/alpha_3_exp, b3_sm)

# Gasim punctele de intersectie
print("Valori la scale caracteristice (din 600-cell):")
scales = [1e3, 1e6, 1e10, 1e14, 1e16, 1e19]
print(f"{'Q (GeV)':<12} {'1/alpha_1':<12} {'1/alpha_2':<12} {'1/alpha_3':<12}")
print("-" * 50)
for scale in scales:
    a1 = running_alpha_inv(scale, M_Z, 1/alpha_1_gut, b1_sm)
    a2 = running_alpha_inv(scale, M_Z, 1/alpha_2, b2_sm)
    a3 = running_alpha_inv(scale, M_Z, 1/alpha_3, b3_sm)
    print(f"{scale:<12.0e} {a1:<12.2f} {a2:<12.2f} {a3:<12.2f}")

# ============================================================
# CAUTAM PUNCTUL DE UNIFICARE
# ============================================================
print("\n" + "=" * 70)
print("4. CAUTAM PUNCTUL DE UNIFICARE")
print("=" * 70)

# Intersectia alpha_1 si alpha_2
# 1/alpha_1(Q) = 1/alpha_2(Q)
# 1/alpha_1(M_Z) - b1*ln(Q/M_Z)/(2*pi) = 1/alpha_2(M_Z) - b2*ln(Q/M_Z)/(2*pi)
# (1/alpha_1 - 1/alpha_2) = (b1 - b2) * ln(Q/M_Z) / (2*pi)
# ln(Q/M_Z) = 2*pi * (1/alpha_1 - 1/alpha_2) / (b1 - b2)

def find_intersection(alpha_inv_1, alpha_inv_2, b1, b2, M_Z):
    """Gaseste scala Q unde alpha_1 = alpha_2."""
    if b1 == b2:
        return None
    ln_Q_MZ = 2 * np.pi * (alpha_inv_1 - alpha_inv_2) / (b1 - b2)
    Q = M_Z * np.exp(ln_Q_MZ)
    return Q

# Intersectii pentru 600-cell
Q_12 = find_intersection(1/alpha_1_gut, 1/alpha_2, b1_sm, b2_sm, M_Z)
Q_13 = find_intersection(1/alpha_1_gut, 1/alpha_3, b1_sm, b3_sm, M_Z)
Q_23 = find_intersection(1/alpha_2, 1/alpha_3, b2_sm, b3_sm, M_Z)

print("Puncte de intersectie (conditii initiale 600-cell):")
print(f"  alpha_1 = alpha_2 la Q_12 = {Q_12:.2e} GeV")
print(f"  alpha_1 = alpha_3 la Q_13 = {Q_13:.2e} GeV")
print(f"  alpha_2 = alpha_3 la Q_23 = {Q_23:.2e} GeV")

# Valoarea la intersectie
if Q_12 and Q_12 > 0:
    alpha_unif_12 = 1/running_alpha_inv(Q_12, M_Z, 1/alpha_1_gut, b1_sm)
    print(f"\n  La Q_12: alpha_unif = {alpha_unif_12:.4f} = 1/{1/alpha_unif_12:.1f}")

# Verificam cat de aproape sunt cele 3 linii
if Q_12 and Q_12 > 0:
    a1_at_Q12 = running_alpha_inv(Q_12, M_Z, 1/alpha_1_gut, b1_sm)
    a2_at_Q12 = running_alpha_inv(Q_12, M_Z, 1/alpha_2, b2_sm)
    a3_at_Q12 = running_alpha_inv(Q_12, M_Z, 1/alpha_3, b3_sm)
    print(f"\n  Verificare la Q_12 = {Q_12:.2e} GeV:")
    print(f"    1/alpha_1 = {a1_at_Q12:.2f}")
    print(f"    1/alpha_2 = {a2_at_Q12:.2f}")
    print(f"    1/alpha_3 = {a3_at_Q12:.2f}")
    print(f"    Discrepanta alpha_3: {abs(a3_at_Q12 - a1_at_Q12):.2f}")

# ============================================================
# COMPARATIE CU SM STANDARD
# ============================================================
print("\n" + "=" * 70)
print("5. COMPARATIE CU SM STANDARD")
print("=" * 70)

Q_12_exp = find_intersection(1/alpha_1_exp, 1/alpha_2_exp, b1_sm, b2_sm, M_Z)
Q_13_exp = find_intersection(1/alpha_1_exp, 1/alpha_3_exp, b1_sm, b3_sm, M_Z)
Q_23_exp = find_intersection(1/alpha_2_exp, 1/alpha_3_exp, b2_sm, b3_sm, M_Z)

print("Puncte de intersectie (valori experimentale):")
print(f"  alpha_1 = alpha_2 la Q_12 = {Q_12_exp:.2e} GeV")
print(f"  alpha_1 = alpha_3 la Q_13 = {Q_13_exp:.2e} GeV")
print(f"  alpha_2 = alpha_3 la Q_23 = {Q_23_exp:.2e} GeV")

print("""
OBSERVATIE:
===========
In SM standard, cele 3 linii NU se intalnesc intr-un singur punct!
Aceasta e motivul pentru care se cauta extensii (SUSY, etc.).
""")

# ============================================================
# CE AR TREBUI PENTRU UNIFICARE EXACTA?
# ============================================================
print("\n" + "=" * 70)
print("6. CONDITII PENTRU UNIFICARE EXACTA")
print("=" * 70)

print("""
Pentru unificare exacta la o scala Q_GUT, trebuie:
  1/alpha_1(Q_GUT) = 1/alpha_2(Q_GUT) = 1/alpha_3(Q_GUT) = 1/alpha_GUT

Din running:
  1/alpha_i(Q_GUT) = 1/alpha_i(M_Z) - (b_i / 2*pi) * ln(Q_GUT / M_Z)

Fie L = ln(Q_GUT / M_Z). Conditiile sunt:
  1/alpha_1 - b_1*L/(2*pi) = 1/alpha_2 - b_2*L/(2*pi) = 1/alpha_3 - b_3*L/(2*pi)

Din primele doua:
  L = 2*pi * (1/alpha_1 - 1/alpha_2) / (b_1 - b_2)

Din ultimele doua:
  L = 2*pi * (1/alpha_2 - 1/alpha_3) / (b_2 - b_3)

Pentru consistenta, trebuie:
  (1/alpha_1 - 1/alpha_2) / (b_1 - b_2) = (1/alpha_2 - 1/alpha_3) / (b_2 - b_3)
""")

# Verificam conditia de consistenta
ratio_12 = (1/alpha_1_gut - 1/alpha_2) / (b1_sm - b2_sm)
ratio_23 = (1/alpha_2 - 1/alpha_3) / (b2_sm - b3_sm)

print(f"Verificare consistenta (600-cell):")
print(f"  (1/a1 - 1/a2) / (b1 - b2) = {ratio_12:.4f}")
print(f"  (1/a2 - 1/a3) / (b2 - b3) = {ratio_23:.4f}")
print(f"  Raport: {ratio_12 / ratio_23:.4f} (trebuie = 1 pentru unificare exacta)")
print(f"  Discrepanta: {abs(1 - ratio_12/ratio_23)*100:.1f}%")

# Pentru SM experimental
ratio_12_exp = (1/alpha_1_exp - 1/alpha_2_exp) / (b1_sm - b2_sm)
ratio_23_exp = (1/alpha_2_exp - 1/alpha_3_exp) / (b2_sm - b3_sm)

print(f"\nVerificare consistenta (experimental):")
print(f"  Raport: {ratio_12_exp / ratio_23_exp:.4f}")
print(f"  Discrepanta: {abs(1 - ratio_12_exp/ratio_23_exp)*100:.1f}%")

# ============================================================
# 7. CE BETA FUNCTIONS AR DA UNIFICARE?
# ============================================================
print("\n" + "=" * 70)
print("7. CE BETA FUNCTIONS AR DA UNIFICARE EXACTA?")
print("=" * 70)

print("""
Daca pastram b_1 si b_2 din SM si cerem unificare,
ce b_3 ar trebui sa avem?

Din conditia de consistenta:
  (1/a2 - 1/a3) / (b2 - b3) = (1/a1 - 1/a2) / (b1 - b2)

Rezolvam pentru b_3:
  b3 = b2 - (1/a2 - 1/a3) * (b1 - b2) / (1/a1 - 1/a2)
""")

b3_needed = b2_sm - (1/alpha_2 - 1/alpha_3) * (b1_sm - b2_sm) / (1/alpha_1_gut - 1/alpha_2)
print(f"Pentru unificare cu constante 600-cell:")
print(f"  b_3 necesar = {b3_needed:.2f}")
print(f"  b_3 SM = {b3_sm:.2f}")
print(f"  Diferenta: {b3_needed - b3_sm:.2f}")

print("""
INTERPRETARE:
=============
b_3 = -7 in SM da running prea lent pentru alpha_3.
Pentru unificare, ar trebui b_3 mai negativ (running mai rapid)
sau noi particule care modifica beta functions.
""")

# ============================================================
# 8. ALTERNATIVA: SCALA DE UNIFICARE DIN 600-CELL
# ============================================================
print("\n" + "=" * 70)
print("8. SCALA DE UNIFICARE DIN GEOMETRIA 600-CELL?")
print("=" * 70)

print("""
IPOTEZA:
========
Daca 600-cell determina constantele, poate determina si
scala de unificare!

Candidati naturali:
  - M_Planck = 1.22 * 10^19 GeV
  - M_GUT tipic = 10^16 GeV
  - M_600 = M_Planck / phi^n pentru vreun n?
""")

M_Planck = 1.22e19  # GeV

# Incercam diverse puteri ale phi
print("Scale bazate pe phi:")
for n in range(1, 20):
    M_n = M_Planck / PHI**n
    if 1e14 < M_n < 1e18:
        print(f"  M_Planck / phi^{n} = {M_n:.2e} GeV")

# Scala din intersectia 1-2
print(f"\nScala din intersectia alpha_1 = alpha_2:")
print(f"  Q_12 = {Q_12:.2e} GeV")
print(f"  log10(Q_12) = {np.log10(Q_12):.2f}")
print(f"  Q_12 / M_Planck = {Q_12/M_Planck:.2e}")
print(f"  ln(M_Planck / Q_12) = {np.log(M_Planck/Q_12):.2f}")

# ============================================================
# 9. GRAFIC RUNNING
# ============================================================
print("\n" + "=" * 70)
print("9. GENERARE GRAFIC")
print("=" * 70)

plt.figure(figsize=(12, 8))

# Running 600-cell
plt.plot(log_Q, alpha1_inv_600, 'b-', linewidth=2, label=r'$\alpha_1^{-1}$ (600-cell)')
plt.plot(log_Q, alpha2_inv_600, 'g-', linewidth=2, label=r'$\alpha_2^{-1}$ (600-cell)')
plt.plot(log_Q, alpha3_inv_600, 'r-', linewidth=2, label=r'$\alpha_3^{-1}$ (600-cell)')

# Running experimental (dashed)
plt.plot(log_Q, alpha1_inv_exp, 'b--', linewidth=1, alpha=0.5, label=r'$\alpha_1^{-1}$ (exp)')
plt.plot(log_Q, alpha2_inv_exp, 'g--', linewidth=1, alpha=0.5, label=r'$\alpha_2^{-1}$ (exp)')
plt.plot(log_Q, alpha3_inv_exp, 'r--', linewidth=1, alpha=0.5, label=r'$\alpha_3^{-1}$ (exp)')

# Markere pentru scale importante
plt.axvline(x=np.log10(M_Z), color='gray', linestyle=':', alpha=0.5)
plt.axvline(x=16, color='purple', linestyle=':', alpha=0.5, label='GUT scale')
plt.axvline(x=19, color='black', linestyle=':', alpha=0.5, label='Planck scale')

plt.xlabel(r'$\log_{10}(Q/\mathrm{GeV})$', fontsize=12)
plt.ylabel(r'$\alpha_i^{-1}$', fontsize=12)
plt.title('Running al Constantelor de Cuplaj (1-loop SM)', fontsize=14)
plt.legend(loc='best')
plt.grid(True, alpha=0.3)
plt.xlim([2, 19])
plt.ylim([0, 70])

plt.savefig('gauge_unification.png', dpi=150, bbox_inches='tight')
print("Grafic salvat: gauge_unification.png")

# ============================================================
# 10. REZUMAT
# ============================================================
print("\n" + "=" * 70)
print("10. REZUMAT EXP-081")
print("=" * 70)

print(f"""
REZULTATE:
==========

1. CONDITII INITIALE (600-cell la M_Z):
   1/alpha_1 = {1/alpha_1_gut:.2f}
   1/alpha_2 = {1/alpha_2:.2f}
   1/alpha_3 = {1/alpha_3:.2f}

2. PUNCTE DE INTERSECTIE:
   alpha_1 = alpha_2 la {Q_12:.2e} GeV
   alpha_2 = alpha_3 la {Q_23:.2e} GeV
   alpha_1 = alpha_3 la {Q_13:.2e} GeV

3. UNIFICARE EXACTA?
   NU - cele 3 linii nu se intalnesc in acelasi punct
   Discrepanta: {abs(1 - ratio_12/ratio_23)*100:.1f}%
   (Similar cu SM standard: {abs(1 - ratio_12_exp/ratio_23_exp)*100:.1f}%)

4. CONCLUZIE:
   Constantele 600-cell dau running SIMILAR cu SM experimental.
   NU rezolva problema neunificarii din SM.

   Pentru unificare exacta ar trebui:
   - Beta functions modificate (noi particule)
   - SAU conditii initiale diferite
   - SAU geometrie 600-cell extinde natural spre E8?

5. OBSERVATIE POZITIVA:
   Diferentele intre 600-cell si SM experimental sunt MICI.
   Teoria e consistenta cu fizica cunoscuta.
""")

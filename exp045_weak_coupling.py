"""
EXP-045: Verificare 1/alpha_weak = 12*phi^2
===========================================
Verificam daca pattern-ul continua pentru forta slaba.
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-045: VERIFICARE PATTERN PENTRU WEAK COUPLING")
print("=" * 70)

# ============================================================
# PATTERN-UL OBSERVAT
# ============================================================
print("\n" + "-" * 70)
print("PATTERN-UL OBSERVAT")
print("-" * 70)

print(f"""
Am gasit:
  1/alpha_bare = 20*phi^4 = {20*PHI**4:.6f}  (electromagnetic, n=4)
  1/alpha_s = 2*phi^3 = {2*PHI**3:.6f}       (strong, n=3)

Daca pattern-ul continua:
  1/alpha_weak = a*phi^2                      (weak, n=2)

INTREBARE: Care e alpha_weak si ce e 'a'?
""")

# ============================================================
# DEFINIREA ALPHA_WEAK
# ============================================================
print("-" * 70)
print("CE E ALPHA_WEAK?")
print("-" * 70)

print("""
In Modelul Standard avem:
- U(1)_Y cu coupling g' (hypercharge)
- SU(2)_L cu coupling g (weak isospin)

Relatiile:
  e = g * sin(theta_W) = g' * cos(theta_W)
  alpha = e^2 / (4*pi)

Deci:
  alpha_1 = g'^2 / (4*pi) = alpha / cos^2(theta_W)  [U(1)]
  alpha_2 = g^2 / (4*pi) = alpha / sin^2(theta_W)   [SU(2)]

Cu sin^2(theta_W) = 0.2312 (experimental) sau 6/26 = 0.2308 (noi):
""")

# Valori experimentale
sin2_tW_exp = 0.23121  # experimental
sin2_tW_ours = 6/26    # formula noastra
cos2_tW_exp = 1 - sin2_tW_exp
cos2_tW_ours = 1 - sin2_tW_ours

# Alpha la scala M_Z (running)
alpha_MZ = 1/128.0  # alpha electromagnetic la M_Z
alpha_0 = ALPHA      # alpha la q^2 = 0

print(f"sin^2(theta_W) experimental = {sin2_tW_exp:.6f}")
print(f"sin^2(theta_W) noi (6/26) = {sin2_tW_ours:.6f}")
print(f"cos^2(theta_W) experimental = {cos2_tW_exp:.6f}")

# ============================================================
# CALCULUL COUPLINGS
# ============================================================
print("\n" + "-" * 70)
print("CALCULUL COUPLING-URILOR")
print("-" * 70)

# Folosind alpha la q=0 (137)
alpha_2_from_0 = alpha_0 / sin2_tW_ours  # SU(2) coupling
alpha_1_from_0 = alpha_0 / cos2_tW_ours  # U(1) coupling

print("Folosind alpha(0) = 1/137:")
print(f"  alpha_2 = alpha/sin^2(tW) = {alpha_2_from_0:.6f}")
print(f"  alpha_1 = alpha/cos^2(tW) = {alpha_1_from_0:.6f}")
print(f"  1/alpha_2 = {1/alpha_2_from_0:.6f}")
print(f"  1/alpha_1 = {1/alpha_1_from_0:.6f}")

# Folosind alpha la M_Z (128)
alpha_2_from_MZ = alpha_MZ / sin2_tW_exp
alpha_1_from_MZ = alpha_MZ / cos2_tW_exp

print(f"\nFolosind alpha(M_Z) = 1/128:")
print(f"  alpha_2 = alpha/sin^2(tW) = {alpha_2_from_MZ:.6f}")
print(f"  alpha_1 = alpha/cos^2(tW) = {alpha_1_from_MZ:.6f}")
print(f"  1/alpha_2 = {1/alpha_2_from_MZ:.6f}")
print(f"  1/alpha_1 = {1/alpha_1_from_MZ:.6f}")

# ============================================================
# VERIFICARE FORMULA 12*phi^2
# ============================================================
print("\n" + "-" * 70)
print("VERIFICARE: 1/alpha_weak = 12*phi^2 ?")
print("-" * 70)

formula_12phi2 = 12 * PHI**2

print(f"12*phi^2 = {formula_12phi2:.6f}")
print()
print("Comparatie cu diverse definitii ale alpha_weak:")
print(f"  1/alpha_2 (din alpha_0) = {1/alpha_2_from_0:.6f}")
print(f"  1/alpha_1 (din alpha_0) = {1/alpha_1_from_0:.6f}")
print(f"  1/alpha_2 (din alpha_MZ) = {1/alpha_2_from_MZ:.6f}")
print(f"  1/alpha_1 (din alpha_MZ) = {1/alpha_1_from_MZ:.6f}")

print(f"\nErori relative:")
print(f"  |12*phi^2 - 1/alpha_2(0)| / (1/alpha_2) = {abs(formula_12phi2 - 1/alpha_2_from_0)/(1/alpha_2_from_0)*100:.2f}%")
print(f"  |12*phi^2 - 1/alpha_2(MZ)| / (1/alpha_2) = {abs(formula_12phi2 - 1/alpha_2_from_MZ)/(1/alpha_2_from_MZ)*100:.2f}%")

# ============================================================
# CAUTAM COEFICIENTUL CORECT
# ============================================================
print("\n" + "-" * 70)
print("CAUTAM COEFICIENTUL CORECT")
print("-" * 70)

# Ce coeficient ar da 1/alpha_2?
coef_for_alpha2_0 = (1/alpha_2_from_0) / PHI**2
coef_for_alpha2_MZ = (1/alpha_2_from_MZ) / PHI**2

print(f"Pentru 1/alpha_2(0) = a*phi^2:")
print(f"  a = {coef_for_alpha2_0:.6f}")
print(f"  Cel mai aproape intreg: {round(coef_for_alpha2_0)}")

print(f"\nPentru 1/alpha_2(MZ) = a*phi^2:")
print(f"  a = {coef_for_alpha2_MZ:.6f}")
print(f"  Cel mai aproape intreg: {round(coef_for_alpha2_MZ)}")

# ============================================================
# TESTAM DIVERSE VALORI PENTRU 'a'
# ============================================================
print("\n" + "-" * 70)
print("TESTAM DIVERSE VALORI PENTRU COEFICIENT")
print("-" * 70)

target = 1/alpha_2_from_0  # 31.62

print(f"Target: 1/alpha_2 = {target:.6f}")
print()

for a in [10, 11, 12, 13, 14, 15, 20]:
    val = a * PHI**2
    err = abs(val - target) / target * 100
    print(f"  {a}*phi^2 = {val:.6f} (eroare {err:.2f}%)")

# ============================================================
# OBSERVATIE: 12 = VECINI PER VERTEX!
# ============================================================
print("\n" + "-" * 70)
print("OBSERVATIE: SEMNIFICATIA LUI 12")
print("-" * 70)

print(f"""
12 = numarul de VECINI per vertex in 600-cell!

In 600-cell:
- 20 tetraedre per vertex  -> 1/alpha_bare = 20*phi^4
- 12 vecini per vertex     -> 1/alpha_2 = 12*phi^2 ?
- 6 decagoane per vertex   -> sin^2(tW) = 6/26

VERIFICARE:
  12*phi^2 = {12*PHI**2:.6f}
  1/alpha_2 = {1/alpha_2_from_0:.6f}

  Eroare: {abs(12*PHI**2 - 1/alpha_2_from_0)/(1/alpha_2_from_0)*100:.2f}%

Hmm, eroarea e cam mare (0.6%). Nu e la fel de buna ca celelalte formule.
""")

# ============================================================
# ABORDARE ALTERNATIVA: FORMULA SELF-CONSISTENTA
# ============================================================
print("-" * 70)
print("ABORDARE ALTERNATIVA")
print("-" * 70)

print("""
Poate pattern-ul e diferit pentru weak.

ALTA IDEE:
Daca sin^2(theta_W) = 6/26, atunci:
  1/alpha_2 = 1/alpha * sin^2(theta_W)
            = (20*phi^4 - 2*pi*alpha) * (6/26)
""")

# Calculam
alpha_solved = (20*PHI**4 - np.sqrt(400*PHI**8 - 8*PI)) / (4*PI)
inv_alpha_solved = 1/alpha_solved
inv_alpha_2_derived = inv_alpha_solved * (6/26)

print(f"1/alpha (din formula noastra) = {inv_alpha_solved:.6f}")
print(f"sin^2(theta_W) = 6/26 = {6/26:.6f}")
print(f"1/alpha_2 = (1/alpha) * sin^2(tW) = {inv_alpha_2_derived:.6f}")

# Verificam daca e a*phi^2
coef = inv_alpha_2_derived / PHI**2
print(f"\nAceasta e a*phi^2 cu a = {coef:.6f}")
print(f"Cel mai aproape: {round(coef)} sau {6/26 * 20 * PHI**2:.4f}")

# ============================================================
# FORMULA DIRECTA DIN 600-CELL
# ============================================================
print("\n" + "-" * 70)
print("DERIVARE DIRECTA DIN 600-CELL")
print("-" * 70)

print(f"""
Daca:
  1/alpha = 20*phi^4 (cu corectie)
  sin^2(theta_W) = 6/26

Atunci:
  alpha_2 = alpha / sin^2(theta_W)
  1/alpha_2 = (1/alpha) * sin^2(theta_W)
            = 20*phi^4 * (6/26) - corectie
            = (120/26)*phi^4 - corectie
            = (60/13)*phi^4 - corectie

  60/13 = {60/13:.6f}
  (60/13)*phi^4 = {(60/13)*PHI**4:.6f}
""")

# ============================================================
# TABEL COMPLET DE PATTERN-URI
# ============================================================
print("\n" + "=" * 70)
print("TABEL COMPLET DE PATTERN-URI")
print("=" * 70)

print(f"""
Forta           | Formula principala      | Valoare    | Experimental | Eroare
----------------|-------------------------|------------|--------------|--------
Electromagnetic | 20*phi^4 - 2*pi*alpha   | 137.036    | 137.036      | 0.00014%
Strong (SU3)    | 1/(2*phi^3) = alpha_s   | 0.1180     | 0.1179       | 0.11%
Weak mixing     | 6/26 = sin^2(tW)        | 0.2308     | 0.2312       | 0.19%
Weak SU(2)      | ???                     | ???        | 31.62        | ???

COEFICIENTI DIN 600-CELL:
  20 = tetraedre per vertex
  12 = vecini per vertex
  6 = decagoane per vertex
  2 = ? (poate legat de SU(3) Cartan dimension)
""")

# ============================================================
# PATTERN FINAL
# ============================================================
print("-" * 70)
print("PATTERN FINAL PROPUS")
print("-" * 70)

print(f"""
IPOTEZA:
========
Couplings sunt determinate de geometria 600-cell:

  Electromagnetic: 1/alpha_bare = 20*phi^4
                   (20 = tetraedre/vertex, phi^4 = 4D scaling)

  Strong:          alpha_s = 1/(2*phi^3)
                   (2 = Cartan dim SU(3)?, phi^3 = ?)

  Weak mixing:     sin^2(theta_W) = 6/26
                   (6 = decagoane/vertex, 26 = 6 + 20)

  Weak SU(2):      1/alpha_2 = (1/alpha) * sin^2(theta_W)
                   (derivat din cele de mai sus)

VERIFICARE NUMERICA:
  1/alpha_2 = {1/alpha_2_from_0:.4f}

  Incercari:
    12*phi^2 = {12*PHI**2:.4f} (eroare {abs(12*PHI**2 - 1/alpha_2_from_0)/1/alpha_2_from_0*100:.2f}%)
    (120/26)*phi^2 = {(120/26)*PHI**2:.4f} (eroare {abs((120/26)*PHI**2 - 1/alpha_2_from_0)/1/alpha_2_from_0*100:.2f}%)
    10*pi = {10*PI:.4f} (eroare {abs(10*PI - 1/alpha_2_from_0)/1/alpha_2_from_0*100:.2f}%)
""")

# 120/26 * phi^2
test1 = (120/26) * PHI**2
print(f"  (120/26)*phi^2 = (varfuri/26)*phi^2 = {test1:.4f}")

# (20*6/26) * phi^2 = (120/26) * phi^2
test2 = (20 * 6 / 26) * PHI**2
print(f"  (20*6/26)*phi^2 = {test2:.4f}")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-045")
print("=" * 70)

print(f"""
1. PATTERN-UL NU E SIMPLU PENTRU WEAK:
   12*phi^2 = {12*PHI**2:.4f}
   1/alpha_2 = {1/alpha_2_from_0:.4f}
   Eroare: ~0.6% (mai mare decat celelalte formule)

2. FORMULA MAI BUNA:
   1/alpha_2 = (1/alpha) * sin^2(theta_W)
             = (20*phi^4 - 2*pi*alpha) * (6/26)

   Aceasta nu e o formula independenta, ci derivata din alpha si theta_W.

3. PATTERN ACTUALIZAT:

   Constanta      | Formula             | Coef | phi^n | Numar 600-cell
   ---------------|---------------------|------|-------|---------------
   1/alpha_bare   | 20*phi^4            | 20   | 4     | tetraedre/vertex
   alpha_s        | 1/(2*phi^3)         | 2    | 3     | ???
   sin^2(theta_W) | 6/26                | 6,26 | -     | decagoane, total

4. OBSERVATIE:
   - Pattern-ul a*phi^n functioneaza pentru alpha si alpha_s
   - Dar pentru weak, formula e derivata (sin^2 = 6/26)
   - Nu avem formula independenta pentru SU(2) coupling

5. INTREBARE DESCHISA:
   De ce coeficientul pentru alpha_s e 2, nu un numar din 600-cell?
   Poate 2 vine din rank(SU(3)) - 1 = 3 - 1 = 2?
   Sau 2 = dim(Cartan subalgebra SU(3))?
""")

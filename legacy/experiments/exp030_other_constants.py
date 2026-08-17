"""
EXP-030: Alte Constante din Geometria 600-cell
==============================================
Putem deriva si alte constante fundamentale?
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-030: ALTE CONSTANTE DIN 600-CELL")
print("=" * 70)

print("""
CONTEXT:
========
Am derivat alpha din 600-cell.
Ce alte constante ar trebui sa poata fi derivate?
""")

print("-" * 70)
print("1. UNGHIUL WEINBERG (sin^2 theta_W)")
print("-" * 70)

# Unghiul Weinberg (parametrul de mixing electroslab)
sin2_theta_W = 0.23121  # valoare experimentala la scala M_Z
theta_W = np.arcsin(np.sqrt(sin2_theta_W))

print(f"sin^2(theta_W) experimental = {sin2_theta_W:.5f}")
print(f"theta_W = {np.degrees(theta_W):.2f} grade")

# Cautam expresii in termeni de phi si geometrie
print(f"\nCautam pattern-uri:")
print(f"  1/4 = {0.25:.5f}")
print(f"  1/4 - 1/phi^4 = {0.25 - 1/PHI**4:.5f}")
print(f"  3/13 = {3/13:.5f}")
print(f"  1/phi^3 = {1/PHI**3:.5f}")

# Hmm, 1/phi^3 = 0.236... aproape de sin^2(theta_W)!
diff_phi3 = abs(1/PHI**3 - sin2_theta_W) / sin2_theta_W * 100
print(f"\n  1/phi^3 = {1/PHI**3:.5f}")
print(f"  sin^2(theta_W) = {sin2_theta_W:.5f}")
print(f"  Diferenta: {diff_phi3:.2f}%")

# Sau 1/(phi^2 + 2) ?
val = 1/(PHI**2 + 2)
diff = abs(val - sin2_theta_W) / sin2_theta_W * 100
print(f"\n  1/(phi^2 + 2) = {val:.5f}")
print(f"  Diferenta: {diff:.2f}%")

# 6/26 = 3/13?
val2 = 6/26
diff2 = abs(val2 - sin2_theta_W) / sin2_theta_W * 100
print(f"\n  6/26 = {val2:.5f}")
print(f"  Diferenta: {diff2:.2f}%")

print("-" * 70)
print("2. CONSTANTA DE CUPLAJ TARE (alpha_s)")
print("-" * 70)

# Alpha_s la scala M_Z
alpha_s = 0.1179  # valoare la M_Z

print(f"alpha_s(M_Z) = {alpha_s:.4f}")
print(f"alpha / alpha_s = {ALPHA / alpha_s:.6f}")

# Raportul alpha_s / alpha
ratio = alpha_s / ALPHA
print(f"alpha_s / alpha = {ratio:.2f}")
print(f"Aproape de 4*pi = {4*PI:.2f}")

# E8 are 240 radacini
# SU(3) (forta tare) e subgrup al E8
print(f"\nDimensiuni grupuri:")
print(f"  dim(E8) = 248")
print(f"  dim(SU(3)) = 8")
print(f"  248/8 = {248/8:.0f}")

print("-" * 70)
print("3. RAPORTUL MASELOR ELECTRON/PROTON")
print("-" * 70)

m_ratio = M_ELECTRON / M_PROTON
print(f"m_e / m_p = {m_ratio:.8f}")
print(f"1 / (m_p/m_e) = {m_ratio:.8f}")

# Cautam expresii
print(f"\nCautam pattern-uri:")
print(f"  1/(12*phi^7) = {1/(12*PHI**7):.8f}")
print(f"  1/(20*phi^6) = {1/(20*PHI**6):.8f}")
print(f"  alpha*phi = {ALPHA*PHI:.8f}")

# Alpha * phi^2 ?
val = ALPHA * PHI**2
print(f"  alpha * phi^2 = {val:.8f}")
diff = abs(val - m_ratio) / m_ratio * 100
print(f"  Diferenta de m_e/m_p: {diff:.2f}%")

# 1/(6*phi^6) ?
val2 = 1/(6*PHI**6)
diff2 = abs(val2 - m_ratio) / m_ratio * 100
print(f"  1/(6*phi^6) = {val2:.8f}, diferenta: {diff2:.2f}%")

print("-" * 70)
print("4. MASA ELECTRONULUI IN UNITATI PLANCK")
print("-" * 70)

m_e_planck = M_ELECTRON / M_PLANCK
print(f"m_e / m_Planck = {m_e_planck:.6e}")

# Expresii
print(f"\nCautam pattern-uri:")
print(f"  alpha^12 = {ALPHA**12:.6e}")
print(f"  exp(-1/alpha) = {np.exp(-1/ALPHA):.6e}")  # prea mic
print(f"  alpha^3 / (2*pi) = {ALPHA**3/(2*PI):.6e}")

# m_e / m_Planck ~ alpha^(ceva) ?
# log(m_e/m_Planck) / log(alpha) = ?
exponent = np.log(m_e_planck) / np.log(ALPHA)
print(f"  log(m_e/m_Planck) / log(alpha) = {exponent:.4f}")
print(f"  Deci m_e/m_Planck ~ alpha^{exponent:.1f}")

# Verificam alpha^10.5
val = ALPHA**10.5
print(f"  alpha^10.5 = {val:.6e}")
print(f"  m_e/m_Planck = {m_e_planck:.6e}")

print("-" * 70)
print("5. NUMERELE 600-CELL SI CONSTANTE")
print("-" * 70)

print("""
NUMERE CHEIE din 600-cell:
  V = 120 (varfuri)
  E = 720 (muchii)
  F = 1200 (fete)
  C = 600 (celule)

  20 = tetraedre/vertex
  12 = vecini/vertex
  5 = tetraedre/muchie
  72 = decagoane
  60 = perechi antipode
""")

# Verificam rapoarte
print(f"\nRapoarte interesante:")
print(f"  720/120 = {720/120:.0f} (muchii/varfuri)")
print(f"  1200/720 = {1200/720:.4f} (fete/muchii)")
print(f"  600/1200 = {600/1200:.4f} (celule/fete)")
print(f"  120 * 20 / 4 = {120*20/4:.0f} = 600 (celule)")

# 1/137 vs numere 600-cell
print(f"\n  1/137 = {1/137:.6f}")
print(f"  (V-1)/E/C = {(120-1)/720/600:.6f}")  # random
print(f"  F/(E*V) = {1200/(720*120):.6f}")

print("-" * 70)
print("6. NUMARUL 137 SI GEOMETRIE")
print("-" * 70)

# 137 e PRIM!
print("137 este NUMAR PRIM.")
print(f"\n137 in diverse baze:")
print(f"  137 = 128 + 8 + 1 = 2^7 + 2^3 + 2^0")
print(f"  137 in binar: {bin(137)}")

# Aproape de numere speciale
print(f"\n137 aproape de:")
print(f"  phi^10 / 2 = {PHI**10/2:.2f}")
print(f"  20*phi^4 = {20*PHI**4:.2f}")
print(f"  e^5 - 10 = {np.exp(5) - 10:.2f}")

print("-" * 70)
print("7. FORMULA PENTRU WEINBERG ANGLE")
print("-" * 70)

# Din geometria 600-cell, poate exista o relatie
# SU(2) x U(1) -> U(1)_EM
# sin^2(theta_W) = g'^2 / (g^2 + g'^2)

# Daca g si g' au expresii geometrice...
# U(1) = fibrele (decagoane) = 72
# SU(2) ~ S^3 = tot spatiul = 120 varfuri

# Raport?
ratio_geo = 72 / 120  # decagoane / varfuri
print(f"Decagoane / Varfuri = 72/120 = {ratio_geo:.4f}")
print(f"sin^2(theta_W) = {sin2_theta_W:.4f}")

# Sau (72*10) / (120*12) ?
ratio2 = (72*10) / (120*12)
print(f"(Dec*10)/(V*12) = {ratio2:.4f}")

# 3/8 (valoare teoretica la tree-level GUT)?
print(f"\n3/8 (GUT tree-level) = {3/8:.4f}")

print("\n" + "=" * 70)
print("REZUMAT EXP-030")
print("=" * 70)

print(f"""
REZULTATE:
==========

1. sin^2(theta_W):
   - 1/phi^3 = {1/PHI**3:.5f} (diferenta {abs(1/PHI**3 - sin2_theta_W)/sin2_theta_W*100:.1f}%)
   - Nu e o potrivire excelenta, dar e aproape

2. alpha_s / alpha ~ 16 (aproape de 4*pi)
   - Ar putea fi o relatie fundamentala

3. m_e / m_p:
   - Nu am gasit expresie simpla in phi

4. m_e / m_Planck ~ alpha^10.5
   - Sugereaza o relatie perturbativa de ordin inalt

CONCLUZIE:
==========
Alpha e singura constanta cu formula PRECISA in phi.
Celelalte constante ar putea necesita:
  - Dinamica completa (nu doar geometrie)
  - Corectii radiative
  - Alte structuri din E8 (nu doar 600-cell)
""")

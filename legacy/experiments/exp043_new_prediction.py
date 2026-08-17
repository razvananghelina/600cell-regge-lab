"""
EXP-043: Predictie Noua - Testul Crucial
========================================
Incercam sa prezicem o constanta pe care NU am verificat-o inca.
Aceasta e diferenta dintre fitting si teorie reala.
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-043: PREDICTIE NOUA - TESTUL CRUCIAL")
print("=" * 70)

print("""
CE AM VERIFICAT DEJA:
====================
1. alpha: 1/alpha = 20*phi^4 - 2*pi*alpha (eroare 0.00014%)
2. sin^2(theta_W) = 6/26 (eroare 0.19%)
3. m_mu/m_e = 8*pi^2*phi^2 (eroare 0.027%)

CE NU AM VERIFICAT:
==================
- m_tau/m_e sau m_tau/m_mu
- m_p/m_e (proton/electron)
- Unghiul Cabibbo
- alpha_s (strong coupling)
- Rapoarte de mase quark
""")

# ============================================================
# DATE EXPERIMENTALE (pentru verificare)
# ============================================================
print("-" * 70)
print("DATE EXPERIMENTALE")
print("-" * 70)

# Mase leptoni (MeV)
m_e = 0.51099895  # electron
m_mu = 105.6583755  # muon
m_tau = 1776.86  # tau

# Rapoarte
ratio_mu_e = m_mu / m_e  # 206.768
ratio_tau_e = m_tau / m_e  # 3477.23
ratio_tau_mu = m_tau / m_mu  # 16.817

print(f"Mase leptoni (MeV):")
print(f"  m_e = {m_e}")
print(f"  m_mu = {m_mu}")
print(f"  m_tau = {m_tau}")
print(f"\nRapoarte experimentale:")
print(f"  m_mu/m_e = {ratio_mu_e:.6f}")
print(f"  m_tau/m_e = {ratio_tau_e:.6f}")
print(f"  m_tau/m_mu = {ratio_tau_mu:.6f}")

# Alte constante
m_proton = 938.272  # MeV
ratio_p_e = m_proton / m_e  # 1836.15

print(f"\nMasa proton: {m_proton} MeV")
print(f"  m_p/m_e = {ratio_p_e:.6f}")

# Unghiul Cabibbo
theta_cabibbo = 13.04  # grade
sin_cabibbo = np.sin(np.radians(theta_cabibbo))
sin2_cabibbo = sin_cabibbo**2

print(f"\nUnghiul Cabibbo: {theta_cabibbo} grade")
print(f"  sin(theta_C) = {sin_cabibbo:.6f}")
print(f"  sin^2(theta_C) = {sin2_cabibbo:.6f}")

# Strong coupling
alpha_s_mZ = 0.1179  # la scala M_Z

print(f"\nStrong coupling la M_Z: alpha_s = {alpha_s_mZ}")

# ============================================================
# PREDICTIE 1: MASA TAU
# ============================================================
print("\n" + "-" * 70)
print("PREDICTIE 1: RAPORTUL m_tau/m_e")
print("-" * 70)

print("""
LOGICA:
- m_mu/m_e = 8*pi^2*phi^2 (generatia 2)
- m_tau/m_e = ? (generatia 3)

Daca pattern-ul continua, ar trebui sa avem un factor suplimentar
pentru a trece de la generatia 2 la generatia 3.

INCERCARI:
""")

# Diverse formule pentru tau/e
predictions_tau_e = [
    ("8*pi^2*phi^2 * phi^2 (factor phi^2)", 8*PI**2*PHI**2 * PHI**2),
    ("8*pi^2*phi^4 (phi^4 in loc de phi^2)", 8*PI**2*PHI**4),
    ("8*pi^3*phi^2 (pi^3 in loc de pi^2)", 8*PI**3*PHI**2),
    ("8*pi^2*phi^2 * tau/mu_exp", 8*PI**2*PHI**2 * ratio_tau_mu),
    ("20*phi^4 * phi^2", 20*PHI**4 * PHI**2),
    ("(8*pi^2*phi^2)^(3/2) / sqrt(phi)", (8*PI**2*PHI**2)**1.5 / np.sqrt(PHI)),
    ("8*pi^2*phi^2 + 20*phi^4", 8*PI**2*PHI**2 + 20*PHI**4),
    ("120*phi^4 (varfuri * phi^4)", 120*PHI**4),
    ("720/phi^2 (muchii / phi^2)", 720/PHI**2),
    ("600*phi^2 (celule * phi^2)", 600*PHI**2),
    ("1200/phi^2 (fete / phi^2)", 1200/PHI**2),
]

print(f"Target experimental: m_tau/m_e = {ratio_tau_e:.4f}")
print(f"\n{'Formula':<45} {'Valoare':<12} {'Eroare %':<10}")
print("-" * 70)

for name, val in predictions_tau_e:
    err = abs(val - ratio_tau_e) / ratio_tau_e * 100
    mark = " <-- BUNA!" if err < 1 else ""
    print(f"{name:<45} {val:<12.4f} {err:<10.2f}{mark}")

# ============================================================
# PREDICTIE 2: RAPORTUL m_tau/m_mu
# ============================================================
print("\n" + "-" * 70)
print("PREDICTIE 2: RAPORTUL m_tau/m_mu")
print("-" * 70)

predictions_tau_mu = [
    ("phi^5", PHI**5),
    ("phi^6 / phi", PHI**6 / PHI),
    ("3*phi^3", 3*PHI**3),
    ("4*phi^2", 4*PHI**2),
    ("2*pi*phi", 2*PI*PHI),
    ("pi^2/phi", PI**2/PHI),
    ("8*phi", 8*PHI),
    ("10*phi", 10*PHI),
    ("12/phi", 12/PHI),
    ("20/phi^2", 20/PHI**2),
    ("6*phi^2", 6*PHI**2),
    ("L_6 - 1 = 17", 17),
    ("L_7 / phi = 29/phi", 29/PHI),
    ("F_8 = 21", 21),
    ("5*pi", 5*PI),
    ("phi^4 + phi^2", PHI**4 + PHI**2),
]

print(f"Target experimental: m_tau/m_mu = {ratio_tau_mu:.4f}")
print(f"\n{'Formula':<30} {'Valoare':<12} {'Eroare %':<10}")
print("-" * 60)

best_tau_mu = None
best_err = 100

for name, val in predictions_tau_mu:
    err = abs(val - ratio_tau_mu) / ratio_tau_mu * 100
    if err < best_err:
        best_err = err
        best_tau_mu = (name, val, err)
    mark = " <-- BUNA!" if err < 2 else ""
    print(f"{name:<30} {val:<12.4f} {err:<10.2f}{mark}")

print(f"\nCea mai buna: {best_tau_mu[0]} = {best_tau_mu[1]:.4f} (eroare {best_tau_mu[2]:.2f}%)")

# ============================================================
# PREDICTIE 3: MASA PROTON
# ============================================================
print("\n" + "-" * 70)
print("PREDICTIE 3: RAPORTUL m_p/m_e")
print("-" * 70)

predictions_p_e = [
    ("6*pi^2*phi^4", 6*PI**2*PHI**4),
    ("8*pi^2*phi^3", 8*PI**2*PHI**3),
    ("20*phi^4 * alpha^(-1) / 10", 20*PHI**4 * ALPHA_INV / 10),
    ("1200*phi^2", 1200*PHI**2),
    ("720*phi^2", 720*PHI**2),
    ("600*phi^3 / pi", 600*PHI**3/PI),
    ("240*phi^4", 240*PHI**4),
    ("120*pi*phi^2", 120*PI*PHI**2),
    ("(m_mu/m_e)^2 / phi^4", (ratio_mu_e)**2 / PHI**4),
    ("alpha^(-1) * phi^5", ALPHA_INV * PHI**5),
    ("20*phi^6", 20*PHI**6),
    ("3 * 600 / phi", 3*600/PHI),
    ("1836 (fitting)", 1836),
]

print(f"Target experimental: m_p/m_e = {ratio_p_e:.4f}")
print(f"\n{'Formula':<40} {'Valoare':<12} {'Eroare %':<10}")
print("-" * 65)

best_p_e = None
best_err = 100

for name, val in predictions_p_e:
    err = abs(val - ratio_p_e) / ratio_p_e * 100
    if err < best_err:
        best_err = err
        best_p_e = (name, val, err)
    mark = " <-- BUNA!" if err < 2 else ""
    print(f"{name:<40} {val:<12.4f} {err:<10.2f}{mark}")

print(f"\nCea mai buna: {best_p_e[0]} = {best_p_e[1]:.4f} (eroare {best_p_e[2]:.2f}%)")

# ============================================================
# PREDICTIE 4: UNGHIUL CABIBBO
# ============================================================
print("\n" + "-" * 70)
print("PREDICTIE 4: UNGHIUL CABIBBO")
print("-" * 70)

predictions_cabibbo = [
    ("1/phi^3", 1/PHI**3),
    ("1/(2*phi^2)", 1/(2*PHI**2)),
    ("phi/8", PHI/8),
    ("1/5", 1/5),
    ("pi/(4*phi^3)", PI/(4*PHI**3)),
    ("1/phi^4", 1/PHI**4),
    ("alpha*phi", ALPHA*PHI),
    ("1/(3*phi^2)", 1/(3*PHI**2)),
    ("sin^2(thetaW)/phi", (6/26)/PHI),
    ("1/20 * phi", PHI/20),
    ("2/(3*pi)", 2/(3*PI)),
]

print(f"Target experimental: sin^2(theta_C) = {sin2_cabibbo:.6f}")
print(f"\n{'Formula':<30} {'Valoare':<12} {'Eroare %':<10}")
print("-" * 55)

best_cab = None
best_err = 100

for name, val in predictions_cabibbo:
    err = abs(val - sin2_cabibbo) / sin2_cabibbo * 100
    if err < best_err:
        best_err = err
        best_cab = (name, val, err)
    mark = " <-- BUNA!" if err < 5 else ""
    print(f"{name:<30} {val:<12.6f} {err:<10.2f}{mark}")

print(f"\nCea mai buna: {best_cab[0]} = {best_cab[1]:.6f} (eroare {best_cab[2]:.2f}%)")

# ============================================================
# PREDICTIE 5: STRONG COUPLING
# ============================================================
print("\n" + "-" * 70)
print("PREDICTIE 5: STRONG COUPLING alpha_s")
print("-" * 70)

predictions_alpha_s = [
    ("1/phi^4", 1/PHI**4),
    ("1/(2*phi^3)", 1/(2*PHI**3)),
    ("phi/20", PHI/20),
    ("1/(8*phi)", 1/(8*PHI)),
    ("pi*alpha", PI*ALPHA),
    ("6*alpha", 6*ALPHA),
    ("1/phi^3 - alpha", 1/PHI**3 - ALPHA),
    ("alpha * phi^3", ALPHA * PHI**3),
    ("1/(3*phi^2)", 1/(3*PHI**2)),
    ("2/20", 2/20),
    ("1/8", 1/8),
]

print(f"Target experimental: alpha_s(M_Z) = {alpha_s_mZ:.6f}")
print(f"\n{'Formula':<30} {'Valoare':<12} {'Eroare %':<10}")
print("-" * 55)

best_alpha_s = None
best_err = 100

for name, val in predictions_alpha_s:
    err = abs(val - alpha_s_mZ) / alpha_s_mZ * 100
    if err < best_err:
        best_err = err
        best_alpha_s = (name, val, err)
    mark = " <-- BUNA!" if err < 5 else ""
    print(f"{name:<30} {val:<12.6f} {err:<10.2f}{mark}")

print(f"\nCea mai buna: {best_alpha_s[0]} = {best_alpha_s[1]:.6f} (eroare {best_alpha_s[2]:.2f}%)")

# ============================================================
# ANALIZA: CE FUNCTIONEAZA?
# ============================================================
print("\n" + "=" * 70)
print("ANALIZA: CE FUNCTIONEAZA?")
print("=" * 70)

print(f"""
REZUMAT PREDICTII:

1. m_tau/m_mu: Cea mai buna e {best_tau_mu[0]} = {best_tau_mu[1]:.4f}
   Experimental: {ratio_tau_mu:.4f}
   Eroare: {best_tau_mu[2]:.2f}%

2. m_p/m_e: Cea mai buna e {best_p_e[0]} = {best_p_e[1]:.4f}
   Experimental: {ratio_p_e:.4f}
   Eroare: {best_p_e[2]:.2f}%

3. sin^2(theta_C): Cea mai buna e {best_cab[0]} = {best_cab[1]:.6f}
   Experimental: {sin2_cabibbo:.6f}
   Eroare: {best_cab[2]:.2f}%

4. alpha_s: Cea mai buna e {best_alpha_s[0]} = {best_alpha_s[1]:.6f}
   Experimental: {alpha_s_mZ:.6f}
   Eroare: {best_alpha_s[2]:.2f}%
""")

# ============================================================
# CAUTAM FORMULA PENTRU TAU/MU MAI SISTEMATIC
# ============================================================
print("-" * 70)
print("CAUTARE SISTEMATICA PENTRU m_tau/m_mu = 16.817")
print("-" * 70)

target = ratio_tau_mu
print(f"Target: {target:.6f}")
print()

# Cautam combinatii de forma a*phi^n + b*pi^m
found = []
for a in range(1, 25):
    for n in range(-4, 6):
        val = a * PHI**n
        err = abs(val - target) / target * 100
        if err < 3:
            found.append((f"{a}*phi^{n}", val, err))

for a in range(1, 15):
    for m in range(1, 4):
        val = a * PI**m / PHI
        err = abs(val - target) / target * 100
        if err < 3:
            found.append((f"{a}*pi^{m}/phi", val, err))

# Sortam dupa eroare
found.sort(key=lambda x: x[2])

print("Formule cu eroare < 3%:")
for name, val, err in found[:10]:
    print(f"  {name:<20} = {val:.4f} (eroare {err:.2f}%)")

# ============================================================
# FORMULA NOUA PENTRU TAU
# ============================================================
print("\n" + "-" * 70)
print("FORMULA CANDIDAT PENTRU m_tau/m_mu")
print("-" * 70)

# Observatie: 16.817 ~ 3*phi^3 + 4 sau 4*phi^2 + 6 sau ...
# Sau: 16.817 ~ L_7 / phi = 29/phi = 17.92... nu prea

# Incercam: phi^4 + phi^2 + 6
test = PHI**4 + PHI**2 + 6
print(f"phi^4 + phi^2 + 6 = {test:.4f} (target: {target:.4f}, eroare: {abs(test-target)/target*100:.2f}%)")

# 4*phi^2 + 6/phi
test = 4*PHI**2 + 6/PHI
print(f"4*phi^2 + 6/phi = {test:.4f} (eroare: {abs(test-target)/target*100:.2f}%)")

# 6*phi^2
test = 6*PHI**2
print(f"6*phi^2 = {test:.4f} (eroare: {abs(test-target)/target*100:.2f}%)")

# 2*pi*phi + 3
test = 2*PI*PHI + 3
print(f"2*pi*phi + 3 = {test:.4f} (eroare: {abs(test-target)/target*100:.2f}%)")

# Daca m_mu/m_e = 8*pi^2*phi^2, atunci m_tau/m_mu ar putea fi legat
# m_tau/m_e = (m_mu/m_e) * (m_tau/m_mu)
# ratio_tau_e = 8*pi^2*phi^2 * x
# x = ratio_tau_e / (8*pi^2*phi^2) = 3477.23 / 206.71 = 16.82

# Ce e 16.82 in termeni de phi?
x = ratio_tau_e / (8*PI**2*PHI**2)
print(f"\nm_tau/m_mu calculat = m_tau/m_e / (8*pi^2*phi^2) = {x:.6f}")
print(f"Experimental = {ratio_tau_mu:.6f}")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-043")
print("=" * 70)

print("""
REZULTATE:
==========
NU am gasit o formula simpla si convingatoare pentru:
- m_tau/m_mu (cea mai apropiata: 6*phi^2, eroare 2.6%)
- m_p/m_e (nicio formula buna)
- sin^2(theta_C) (nicio formula buna)
- alpha_s (cea mai apropiata: 1/phi^4 sau 2/20, erori ~5-15%)

INTERPRETARE:
=============
1. POSIBIL: Formula 6*phi^2 pentru m_tau/m_mu are potential
   - 6 = decagoane per vertex in 600-cell
   - phi^2 = factor geometric
   - Eroare 2.6% - nu e grozava dar nici rau

2. Celelalte constante (m_p/m_e, Cabibbo, alpha_s) nu par
   sa aiba formule simple cu phi si numere din 600-cell.

3. Aceasta ar putea insemna:
   a) Teoria noastra e incompleta
   b) Aceste constante vin din alta fizica (QCD, etc.)
   c) Formulele pentru alpha/Weinberg/muon sunt coincidente

PROPUNERE PREDICTIE:
====================
Daca acceptam 6*phi^2 pentru m_tau/m_mu:

  m_tau/m_mu = 6*phi^2 = 15.708 (predictie)
  Experimental = 16.817
  Eroare = 6.6%

Aceasta NU e o predictie buna. Eroarea e prea mare.

VERDICT:
========
Nu am reusit sa facem o PREDICTIE NOUA convingatoare.
Aceasta e o problema pentru teoria 600-cell.
""")

# ============================================================
# ULTIMA INCERCARE: COMBINATII
# ============================================================
print("\n" + "-" * 70)
print("ULTIMA INCERCARE: m_tau/m_e direct")
print("-" * 70)

# m_tau/m_e = 3477.23
target_tau_e = ratio_tau_e

# Incercam combinatii mai complexe
tests = [
    ("8*pi^2*phi^4", 8*PI**2*PHI**4),
    ("8*pi^3*phi^2", 8*PI**3*PHI**2),
    ("48*pi^2*phi^2", 48*PI**2*PHI**2),
    ("8*pi^2*phi^2 * 6*phi^2 / (2*pi)", 8*PI**2*PHI**2 * 6*PHI**2 / (2*PI)),
    ("20*phi^6 * pi", 20*PHI**6 * PI),
    ("600*phi^3", 600*PHI**3),
    ("1200*phi^2 / phi^(-1)", 1200*PHI**2 * PHI),
    ("120*phi^4 * phi^2 / 5", 120*PHI**4 * PHI**2 / 5),
    ("8*pi^2*phi^2 * phi^4 / phi^2", 8*PI**2*PHI**4),
]

print(f"Target: m_tau/m_e = {target_tau_e:.4f}")
print()
for name, val in tests:
    err = abs(val - target_tau_e) / target_tau_e * 100
    mark = " <--" if err < 2 else ""
    print(f"{name:<45} = {val:<12.4f} (eroare {err:.2f}%){mark}")

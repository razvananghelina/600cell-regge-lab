"""
EXP-057: Particule la Scale Multiple
====================================
Ideea: 600-cell e la scala Planck, dar particulele
pot fi rezonante care se intind pe MAI MULTI pixeli.

Fotonul: merge pixel cu pixel, nu rezoneaza
Electronul: bucla care se intinde pe N pixeli
Muonul: bucla mai mica (mai putini pixeli) = frecventa mai mare = masa mai mare
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-057: PARTICULE LA SCALE MULTIPLE")
print("=" * 70)

# ============================================================
# IDEEA DE BAZA
# ============================================================
print("\n" + "-" * 70)
print("IDEEA DE BAZA")
print("-" * 70)

print("""
FOTONUL:
========
- Excitatie care se PROPAGA
- Merge 1 pixel per tick
- NU rezoneaza (nu sta pe loc)
- Energia = h*f (depinde de frecventa undei)
- Masa = 0 (pentru ca nu e "prins")

PARTICULA MASIVA (electron, muon, etc.):
========================================
- Excitatie care REZONEAZA (sta pe loc sau se misca incet)
- E "prinsa" intr-o bucla
- Bucla poate fi de orice marime (nu doar 1 pixel)

MASA din E = mc^2:
- Daca particula e o rezonanta cu frecventa f:
  E = h * f
  m = E/c^2 = h*f/c^2

- Frecventa depinde de MARIMEA buclei:
  f = c / L  (unde L = circumferinta buclei)

- Deci:
  m = h / (c * L)

MASA MAI MARE = BUCLA MAI MICA!

Asta e contra-intuitiv dar corect:
- Electron (masa mica) = bucla MARE
- Muon (masa mare) = bucla MAI MICA
- Tau (masa si mai mare) = bucla SI MAI MICA
""")

# ============================================================
# CALCULUL MARIMII BUCLELOR
# ============================================================
print("-" * 70)
print("CALCULUL MARIMII BUCLELOR")
print("-" * 70)

# Constante
h = 6.626e-34  # J*s
c = 3e8        # m/s
m_e = 9.109e-31  # kg
m_mu = 1.883e-28  # kg
m_tau = 3.167e-27  # kg
l_planck = 1.616e-35  # m

print(f"Lungime Planck: l_P = {l_planck:.3e} m")
print(f"Masa electron: m_e = {m_e:.3e} kg")
print(f"Masa muon: m_mu = {m_mu:.3e} kg")
print(f"Masa tau: m_tau = {m_tau:.3e} kg")

# Din m = h / (c * L), scoatem L = h / (m * c)
# Aceasta e LUNGIMEA DE UNDA COMPTON

L_e = h / (m_e * c)
L_mu = h / (m_mu * c)
L_tau = h / (m_tau * c)

print(f"\nLungimi de unda Compton (= marimea buclei?):")
print(f"  L_e = h/(m_e*c) = {L_e:.3e} m")
print(f"  L_mu = h/(m_mu*c) = {L_mu:.3e} m")
print(f"  L_tau = h/(m_tau*c) = {L_tau:.3e} m")

# In unitati Planck
print(f"\nIn unitati Planck (l_P = 1):")
print(f"  L_e / l_P = {L_e / l_planck:.2e}")
print(f"  L_mu / l_P = {L_mu / l_planck:.2e}")
print(f"  L_tau / l_P = {L_tau / l_planck:.2e}")

# ============================================================
# INTERPRETARE
# ============================================================
print("\n" + "-" * 70)
print("INTERPRETARE")
print("-" * 70)

print(f"""
Electronul are bucla de ~{L_e/l_planck:.0e} pixeli Planck!
Muonul are bucla de ~{L_mu/l_planck:.0e} pixeli Planck!

Asta e ENORM - mult mai mare decat un singur 600-cell.

CONCLUZIE:
Daca 600-cell e la scala Planck, particulele sunt rezonante
care se intind pe MILIARDE de 600-cells!

Sau poate 600-cell nu e la scala Planck?
""")

# ============================================================
# ALTERNATIVA: 600-CELL LA ALTA SCALA
# ============================================================
print("-" * 70)
print("ALTERNATIVA: 600-CELL LA ALTA SCALA")
print("-" * 70)

print("""
Poate 600-cell nu e la scala Planck, ci la scala ELECTRONULUI.

Daca raza 600-cell = lungimea Compton a electronului:
  R_600cell = L_e = 2.43e-12 m

Atunci:
  - Electronul = rezonanta pe UN 600-cell
  - Muonul = rezonanta pe o FRACTIUNE din 600-cell
  - Sau alta interpretare...
""")

# Raportul
R_600cell_electron = L_e
print(f"Daca R_600cell = L_e = {R_600cell_electron:.3e} m:")
print(f"  Raport L_e/L_mu = {L_e/L_mu:.2f} = m_mu/m_e [corect!]")
print(f"  Raport L_e/L_tau = {L_e/L_tau:.2f} = m_tau/m_e [corect!]")

# ============================================================
# MODEL: REZONANTE PE DECAGON LA DIFERITE SCALE
# ============================================================
print("\n" + "-" * 70)
print("MODEL: REZONANTE PE DECAGON")
print("-" * 70)

print("""
Ipoteza: Toate leptonii sunt rezonante pe DECAGON,
dar la scale diferite (moduri diferite de vibratie).

Decagonul fundamental: 10 pasi, fiecare de lungime d.
Circumferinta: L = 10*d

Moduri de vibratie:
  Modul n: lungime de unda = L/n, frecventa f_n = n*c/L

Masa modului n:
  m_n = h*f_n/c^2 = h*n/(c*L) = n * m_1

Unde m_1 = masa modului fundamental.

DACA electron = modul 1, muon = modul n:
  m_mu/m_e = n
  n = 207

Deci muonul ar fi modul 207 pe acelasi decagon!
""")

# Verificare
print("Verificare:")
print(f"  m_mu/m_e = {m_mu/m_e:.2f}")
print(f"  m_tau/m_e = {m_tau/m_e:.2f}")
print(f"  m_tau/m_mu = {m_tau/m_mu:.2f}")

# ============================================================
# PROBLEMA: DE CE DOAR 3 LEPTONI?
# ============================================================
print("\n" + "-" * 70)
print("PROBLEMA: DE CE DOAR 3 LEPTONI?")
print("-" * 70)

print("""
Daca modurile sunt n = 1, 2, 3, ..., de ce avem doar 3 leptoni?

Electron: n = ?
Muon: n = ?
Tau: n = ?

Trebuie sa gasim ce e special la aceste 3 moduri.

IPOTEZA 1: Moduri stabile
- Doar anumite moduri sunt stabile (nu decad instant)
- Poate modurile prime? Sau Fibonacci?

IPOTEZA 2: Legate de geometria 600-cell
- Poate n = numere din 600-cell?
- n_e = 1, n_mu = ?, n_tau = ?

IPOTEZA 3: Raportul e important, nu valorile absolute
- m_mu/m_e ~ 207 ~ ?
- m_tau/m_mu ~ 17 ~ ?
""")

# ============================================================
# CAUTAM PATTERN-UL
# ============================================================
print("-" * 70)
print("CAUTAM PATTERN-UL IN RAPOARTE")
print("-" * 70)

ratio_mu_e = m_mu / m_e
ratio_tau_mu = m_tau / m_mu
ratio_tau_e = m_tau / m_e

print(f"m_mu/m_e = {ratio_mu_e:.4f}")
print(f"m_tau/m_mu = {ratio_tau_mu:.4f}")
print(f"m_tau/m_e = {ratio_tau_e:.4f}")

# Verificam daca sunt legate de phi
print(f"\nIn termeni de phi:")
print(f"  m_mu/m_e / phi^11 = {ratio_mu_e / PHI**11:.4f}")
print(f"  m_tau/m_mu / phi^6 = {ratio_tau_mu / PHI**6:.4f}")

# Verificam daca 207 si 17 au ceva special
print(f"\nFactorizari:")
print(f"  207 = 9 * 23 = 3^2 * 23")
print(f"  17 = prim")
print(f"  207 / 17 = {207/17:.2f} ~ 12.2")

# ============================================================
# IPOTEZA: MODURI PE STRUCTURI DIFERITE
# ============================================================
print("\n" + "-" * 70)
print("IPOTEZA: MODURI PE STRUCTURI DIFERITE")
print("-" * 70)

print("""
Poate fiecare lepton e pe o STRUCTURA diferita, nu doar mod diferit:

Electron: rezonanta pe DECAGON (10 pasi)
Muon: rezonanta pe TRIUNGHI (3 pasi)?
Tau: rezonanta pe MUCHIE (2 pasi)?

Sau invers (masa mai mare = structura mai mica):

Electron: rezonanta pe structura MARE
Muon: rezonanta pe structura MEDIE
Tau: rezonanta pe structura MICA

Sa verificam rapoartele:
""")

# Daca masa ~ 1/L, atunci m1/m2 = L2/L1
# Pentru decagon (10), triunghi (3), muchie (2):
print("Daca electron=decagon, muon=triunghi, tau=muchie:")
print(f"  m_mu/m_e = L_decagon/L_triunghi = 10/3 = {10/3:.2f}")
print(f"  Experimental: {ratio_mu_e:.2f}")
print(f"  FOARTE DIFERIT!")

print("\nDaca combinam structuri cu moduri:")
print("  Electron: decagon, mod 1 -> L_eff = 10")
print("  Muon: decagon, mod n -> L_eff = 10/n")
print("  m_mu/m_e = n")
print(f"  Deci n = {ratio_mu_e:.0f}")

# ============================================================
# MODEL CONCRET: STANDING WAVES
# ============================================================
print("\n" + "-" * 70)
print("MODEL CONCRET: UNDE STATIONARE")
print("-" * 70)

print("""
Imaginea finala:

1. SPATIUL e "tesut" din 600-cells la scala Planck
   (sau la alta scala fundamentala)

2. FOTONUL e o unda care se propaga liber
   - Nu e legat de nicio structura
   - Merge cu c (1 pixel/tick)

3. ELECTRONUL e o unda stationara pe o bucla MARE
   - Bucla are circumferinta L_e (Compton wavelength)
   - L_e = h/(m_e*c) = 2.43e-12 m
   - In pixeli Planck: ~10^23 pixeli
   - E o rezonanta de frecventa joasa

4. MUONUL e o unda stationara pe o bucla MAI MICA
   - L_mu = L_e / 207
   - Frecventa mai mare = masa mai mare

5. TAU e o unda stationara pe o bucla SI MAI MICA
   - L_tau = L_mu / 17
   - Frecventa si mai mare = masa si mai mare

INTREBARE CHEIE:
Ce determina marimile L_e, L_mu, L_tau?
De ce exact aceste valori si nu altele?
""")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-057")
print("=" * 70)

print(f"""
CE AM INTELES:

1. Masa vine din FRECVENTA rezonantei:
   m = h*f/c^2 = h/(c*L)

   Bucla mai mica = frecventa mai mare = masa mai mare

2. Electronul are bucla de ~10^23 pixeli Planck
   - E ENORM comparativ cu un 600-cell
   - 600-cell la scala Planck are 120 pixeli

3. Rapoartele de mase:
   m_mu/m_e = 207 = raport de circumferinte inverse
   m_tau/m_mu = 17 = raport de circumferinte inverse

4. PROBLEMA RAMASA:
   De ce 207 si 17?
   Ce determina aceste numere specifice?

   Posibil: sunt legate de moduri de vibratie permise
   pe geometria 600-cell (sau alta structura).

5. OBSERVATIE:
   207 ~ phi^11 / 0.96
   17 ~ phi^6 / 1.06

   Poate sunt puteri ale lui phi cu corectii mici?

DIRECTIE:
Trebuie inteles ce face ca anumite bucle sa fie STABILE
(sa nu decada) si altele nu. Poate e legat de:
- Simetria buclei
- Numarul de pasi (trebuie sa "inchida" corect)
- Interactiunea cu alte campuri (Higgs?)
""")

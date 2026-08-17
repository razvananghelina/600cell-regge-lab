"""
EXP-060: De la Energie la Masa - Mecanismul
===========================================
Ideea: La Big Bang, energia era atat de concentrata incat
nu putea sa se propage liber. S-a "blocat" in bucle,
creand particule.

E = mc^2 invers: E -> m cand energia e "prinsa"
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-060: DE LA ENERGIE LA MASA")
print("=" * 70)

# ============================================================
# IDEEA FUNDAMENTALA
# ============================================================
print("\n" + "-" * 70)
print("IDEEA FUNDAMENTALA")
print("-" * 70)

print("""
FOTONUL LIBER:
- Energie care se PROPAGA
- Merge de la pixel la pixel
- Nu sta pe loc
- Masa = 0 (pentru ca nu e localizat)

CE SE INTAMPLA LA DENSITATE MARE DE ENERGIE?

Imagineaza: 1000 de fotoni incearca sa treaca prin acelasi vertex.
- Nu pot trece toti simultan
- Se "lovesc" si interfereaza
- Unii raman "blocati" in bucle

BUCLA DE ENERGIE = PARTICULA CU MASA

E ca un ambuteiaj:
- Masini libere = fotoni (se misca)
- Ambuteiaj = particula (energie blocata)
- Ambuteiajul NU se misca (sau se misca incet)
- Dar contine energie!
""")

# ============================================================
# CONDITIA DE "BLOCARE"
# ============================================================
print("-" * 70)
print("CONDITIA DE BLOCARE")
print("-" * 70)

print("""
Cand se "blocheaza" energia intr-o bucla?

CONDITIE 1: Densitate critica
- Daca energia per pixel > E_critic, apar bucle
- E_critic ar putea fi legat de geometria 600-cell

CONDITIE 2: Rezonanta
- Bucla trebuie sa "se inchida" corect (faza = 2*pi*n)
- Nu orice bucla e stabila, doar cele rezonante

CONDITIE 3: Stabilitate
- Bucla trebuie sa nu se "destrame" imediat
- Anumite configuratii sunt stabile, altele nu

La Big Bang:
- Densitatea era ENORMA (mult peste E_critic)
- Toata energia era "blocata" in bucle
- Pe masura ce universul s-a extins, densitatea a scazut
- Unele bucle s-au "destins" si au devenit fotoni liberi
- Altele au ramas stabile = particulele pe care le vedem azi
""")

# ============================================================
# MODEL CANTITATIV SIMPLU
# ============================================================
print("-" * 70)
print("MODEL CANTITATIV")
print("-" * 70)

print("""
Fie:
- E = energia totala intr-o regiune
- V = volumul regiunii (in pixeli)
- rho = E/V = densitatea de energie

DENSITATEA CRITICA:
- Cand rho > rho_critic, apar particule
- rho_critic ar trebui sa fie legat de scala Planck

La scala Planck:
- E_Planck = sqrt(hbar * c^5 / G) = 1.22 * 10^19 GeV
- L_Planck = 1.62 * 10^-35 m
- rho_Planck = E_Planck / L_Planck^3 = densitatea maxima posibila
""")

# Calcule
E_planck = 1.22e19  # GeV
L_planck = 1.62e-35  # m
rho_planck = E_planck / L_planck**3  # GeV/m^3

print(f"E_Planck = {E_planck:.2e} GeV")
print(f"L_Planck = {L_planck:.2e} m")
print(f"rho_Planck = {rho_planck:.2e} GeV/m^3")

# ============================================================
# INTENSITATEA IN BUCLA
# ============================================================
print("\n" + "-" * 70)
print("INTENSITATEA IN BUCLA")
print("-" * 70)

print("""
Ideea ta: "intensitatea de fotoni intr-un loop"

Daca N fotoni sunt "prinsi" intr-o bucla de L pixeli:
- Energia totala: E = N * E_foton
- Dar fotonii interfereaza!

INTERFERENTA CONSTRUCTIVA:
- Daca fazele se aliniaza: E_total = N^2 * E_foton (intensitate)
- E ca in laser: N fotoni coerenti au intensitate N^2

PENTRU PARTICULA STABILA:
- Trebuie interferenta constructiva PERMANENTA
- Asta inseamna faza totala pe bucla = 2*pi*n

Daca bucla are L pasi, fiecare pas adauga faza delta:
  L * delta = 2*pi*n
  delta = 2*pi*n / L

Pentru decagon (L=10):
  delta = 2*pi*n / 10 = n * 36 grade
""")

# ============================================================
# MASA DIN ENERGIE BLOCATA
# ============================================================
print("-" * 70)
print("MASA DIN ENERGIE BLOCATA")
print("-" * 70)

print("""
Daca E e energia blocata in bucla:
  m = E / c^2

Energia depinde de:
1. Numarul de "fotoni" in bucla (N)
2. Frecventa fiecarui foton (f)
3. Cum interfereaza (constructiv sau nu)

IPOTEZA SIMPLA:
  E_bucla = h * f_bucla
  f_bucla = c / L_bucla

  m = h * f / c^2 = h / (c * L)

Deci: masa ~ 1 / marimea buclei

ELECTRON vs MUON:
  m_mu / m_e = L_e / L_mu = 207

Electronul e o bucla de 207x mai mare decat muonul!
(sau muonul are 207x mai multi "fotoni" in aceeasi bucla)
""")

# ============================================================
# DE CE ANUMITE MASE?
# ============================================================
print("-" * 70)
print("DE CE ANUMITE MASE?")
print("-" * 70)

print("""
De ce exista electronul cu masa lui specifica?
De ce nu particule cu orice masa?

RASPUNS POSIBIL: REZONANTA PE 600-CELL

Geometria 600-cell permite doar anumite bucle stabile.
Bucla trebuie sa:
1. Se inchida (sa revina la punctul initial)
2. Sa aiba faza corecta (interferenta constructiva)
3. Sa fie stabila (sa nu se destrame)

CONDITIA DE FAZA:
  Faza totala = k * L = 2*pi*n

  Unde k = numar de unda = 2*pi / lambda

  Deci: L = n * lambda

BUCLE PERMISE PE 600-CELL:
- Decagon: L = 10 pasi
- Triunghi: L = 3 pasi
- Pentagon: L = 5 pasi
- Etc.

Dar nu doar lungimea conteaza - si TOPOLOGIA!
- O bucla pe un decagon e diferita de o bucla pe alt decagon
- Pot fi 6 decagoane diferite prin fiecare vertex
- Combinatiile de bucle dau particule diferite
""")

# ============================================================
# MODELUL COMPLET
# ============================================================
print("\n" + "-" * 70)
print("MODELUL COMPLET: BIG BANG -> PARTICULE")
print("-" * 70)

print("""
SCENARIUL:

1. INAINTE DE BIG BANG (sau la t=0):
   - Toata energia "comprimata" in spatiu minim
   - Densitate >> rho_Planck
   - Spatiul insusi e deformat/comprimat

2. BIG BANG (t ~ t_Planck):
   - Energia incepe sa se "propage"
   - Dar densitatea e inca prea mare
   - Energia se "blocheaza" in bucle
   - Apar TOATE tipurile posibile de particule

3. EXPANSIUNE (t > t_Planck):
   - Universul se extinde
   - Densitatea scade
   - Particulele INSTABILE decad (bucle care se "destram")
   - Particulele STABILE raman (bucle rezonante)

4. ASTAZI:
   - Densitate mica
   - Doar particule stabile exista natural
   - Putem crea particule instabile in acceleratoare
     (concentram energie local)

PARTICULELE STABILE = BUCLE REZONANTE PE 600-CELL
- Electron: bucla stabila de un anumit tip
- Proton: combinatie de 3 bucle (quarks) legate
- Foton: energie care NU e in bucla (se propaga liber)
""")

# ============================================================
# CONEXIUNEA CU NUMERELE NOASTRE
# ============================================================
print("-" * 70)
print("CONEXIUNEA CU NUMERELE")
print("-" * 70)

print(f"""
Am gasit ca masele urmeaza phi^n cu n din 600-cell.

INTERPRETARE:
- n = numarul de "moduri" sau "quante" in bucla
- phi = factorul geometric fundamental
- phi^n = frecventa relativa a buclei

DE CE PHI?
- phi apare natural in geometria 600-cell
- Distanta intre vecini = 1/phi
- Unghiuri = multipli de 36 grade (legate de phi)
- phi^n e "eigenvalue" natural al geometriei

DE CE n = 0, 3, 5, 11, 17, ...?
- Acestea sunt buclele STABILE
- 5 si 6 sunt "blocuri de baza" (din geometrie)
- 3 ar putea fi legat de culoare (pentru quarks)
- Combinatiile stabile: 5+6=11, 11+6=17, etc.

INTREBARE DESCHISA:
Ce face ca n=11 sa fie stabil dar n=10 nu?
Trebuie sa intelegem STABILITATEA buclelor.
""")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-060")
print("=" * 70)

print("""
MODELUL PROPUS:

1. ENERGIA LIBERA = unde care se propaga (fotoni)
   - Masa = 0
   - Viteza = c

2. ENERGIA BLOCATA = unde statice in bucle (particule)
   - Masa = E/c^2 = h/(c*L)
   - Viteza < c

3. BIG BANG = energie la densitate extrema
   - Prea multa energie pentru a se propaga liber
   - Se formeaza bucle (particule)

4. GEOMETRIA 600-CELL determina:
   - Ce bucle sunt POSIBILE (topologie)
   - Ce bucle sunt STABILE (rezonanta)
   - Ce MASE au particulele (phi^n)

5. PARTICULELE SUNT "AMBUTEIAJE" DE ENERGIE
   - Electronul = ambuteiaj mic, stabil
   - Muonul = ambuteiaj mai "strans" (masa mai mare)
   - Fotonul = trafic liber (masa = 0)

RAMANE DE INTELES:
- DE CE anumite bucle sunt stabile (n = 0, 3, 5, 11, 17)?
- Cum se COMBINA buclele (protonul = 3 quarks)?
- Care e MECANISMUL exact de stabilitate?
""")

"""
EXP-055: Spatiu-timp Discret pe 600-cell
========================================
Imaginea:
- Fiecare vertex = 1 pixel (pozitie posibila)
- Lumina se misca cu 1 pixel/tick
- Particulele = pattern-uri stabile (rezonante)
- Undele = excitatii care se propaga

Intrebare: Cum devine unda particula si invers?
"""

from physics_formulas import *
import numpy as np

print("=" * 70)
print("EXP-055: SPATIU-TIMP DISCRET PE 600-CELL")
print("=" * 70)

# ============================================================
# STRUCTURA DE BAZA
# ============================================================
print("\n" + "-" * 70)
print("STRUCTURA DISCRETA")
print("-" * 70)

print("""
IMAGINEA FIZICA:
================

SPATIUL:
- 120 vertexuri = 120 "pixeli" posibili
- Fiecare vertex are 12 vecini
- Distanta intre vecini = 1 (unitate Planck)
- Raza totala = 1 (sau l_Planck in unitati reale)

TIMPUL:
- Discret: tick-uri de timp Planck
- t_Planck = 5.39 * 10^-44 secunde
- La fiecare tick, starea se updateaza

LUMINA:
- Se misca exact 1 pixel per tick
- c = 1 pixel / 1 tick = l_Planck / t_Planck
- Aceasta e DEFINITIA vitezei luminii in acest model

PARTICULE:
- Pattern-uri STABILE pe retea
- Trebuie sa "rezoneze" - sa se auto-sustina
- Miscarea = propagarea pattern-ului
""")

# ============================================================
# TIPURI DE EXCITATII
# ============================================================
print("-" * 70)
print("TIPURI DE EXCITATII PE 600-CELL")
print("-" * 70)

print("""
1. UNDA SIMPLA (foton):
   - Excitatie la un vertex
   - Se propaga la vecini la fiecare tick
   - Se "imprastie" in 12 directii
   - Amplitudinea scade cu 1/sqrt(12) per pas

   Dupa n pasi: amplitudine ~ (1/sqrt(12))^n = 12^(-n/2)

2. BUCLA INCHISA (particula stabila?):
   - Excitatie care revine la punctul initial
   - Cea mai scurta bucla = DECAGONUL (10 pasi)
   - Dupa 10 pasi, unda revine si interfereaza cu ea insasi

   Conditie de rezonanta: faza dupa 10 pasi = multiplu de 2*pi

3. PATTERN LOCALIZAT (particula masiva?):
   - Excitatie care NU se propaga (sau se propaga incet)
   - Trebuie sa fie "prinsa" in geometria locala
   - 20 tetraedre per vertex = 20 "cavitati" locale
""")

# ============================================================
# CONDITIA DE REZONANTA PE DECAGON
# ============================================================
print("-" * 70)
print("CONDITIA DE REZONANTA PE DECAGON")
print("-" * 70)

print("""
Decagonul = geodezica inchisa cu 10 pasi.

O unda care face o bucla completa acumuleaza faza:
  phi_total = k * L = k * 10 * (distanta per pas)

Unde k = numar de unda.

Pentru REZONANTA (unda stabila):
  phi_total = n * 2*pi  (n = intreg)

Deci:
  k * 10 * d = n * 2*pi
  k = n * 2*pi / (10*d)

Daca d = 1/phi (distanta dintre vecini in 600-cell):
  k = n * 2*pi * phi / 10

ENERGIE (E = hbar * omega = hbar * c * k):
  E_n = n * hbar * c * (2*pi*phi/10)
      = n * E_0

Unde E_0 = hbar * c * 2*pi*phi/10 = energie fundamentala
""")

# Calcul numeric
d_edge = 1/PHI  # distanta intre vecini
L_decagon = 10 * d_edge  # lungime decagon
print(f"Distanta intre vecini: d = 1/phi = {d_edge:.6f}")
print(f"Lungime decagon: L = 10/phi = {L_decagon:.6f}")
print(f"Circumferinta cerc (comparatie): 2*pi = {2*PI:.6f}")
print(f"Raport L/2*pi = {L_decagon/(2*PI):.6f} (decagonul e 98.4% din cerc)")

# ============================================================
# UNDE VS PARTICULE
# ============================================================
print("\n" + "-" * 70)
print("MECANISMUL UNDA -> PARTICULA")
print("-" * 70)

print("""
UNDA (delocalizata):
- Excitatie care se propaga liber
- Ocupa multi pixeli simultan (superpozitie)
- Energia distribuita pe tot spatiul
- Exemplu: foton in propagare

PARTICULA (localizata):
- Excitatie care sta "pe loc" (sau se misca lent)
- Ocupa o regiune mica (cativa pixeli)
- Energia concentrata local
- Exemplu: electron in repaus

TRANZITIA UNDA -> PARTICULA:
=============================
Cand o unda intalneste o "capcana geometrica" unde poate rezona.

Pe 600-cell, "capcanele" sunt:
1. DECAGOANE - bucle de 10 pasi
2. TETRAEDRE - cavitati locale (20 per vertex)
3. ICOSAEDRE - vertex figures (12 vecini formand icosaedru)

Daca unda are LUNGIMEA DE UNDA POTRIVITA, ramane "prinsa".
""")

# ============================================================
# CALCULUL ALPHA DIN ACEASTA IMAGINE
# ============================================================
print("-" * 70)
print("ALPHA DIN PERSPECTIVA DISCRETA")
print("-" * 70)

print("""
INTREBARE: Ce e alpha in acest model?

Alpha = probabilitatea ca doua particule sa interactioneze
        cand se intalnesc la acelasi vertex.

MODEL SIMPLU:
- Particula A e la vertex V
- Particula B ajunge la vertex V
- Probabilitatea de interactiune = ?

GEOMETRIE:
- Particula A "ocupa" o fractiune din spatiul local
- Spatiul local la V = 20 tetraedre
- Fiecare tetraedru are "volum" ~ phi^4 (in 4D)

PROBABILITATE:
- Daca A ocupa 1 tetraedru din 20: P = 1/20
- Dar A e "diluata" pe volum phi^4: P = 1/(20*phi^4)

Aceasta da: alpha = 1/(20*phi^4) = 0.00729...

CORECTIA 2*pi*alpha:
- Cand A interactioneaza cu sine insasi (self-energy)
- Unda face o bucla pe decagon (faza 2*pi)
- Probabilitate self-interaction = 2*pi * alpha
""")

# Verificare
alpha_model = 1/(20*PHI**4)
print(f"Alpha din model: {alpha_model:.10f}")
print(f"Alpha experimental: {ALPHA:.10f}")
print(f"Eroare: {abs(alpha_model - ALPHA)/ALPHA * 100:.4f}%")

# ============================================================
# DE CE 20 SI NU 12?
# ============================================================
print("\n" + "-" * 70)
print("DE CE 20 SI NU 12?")
print("-" * 70)

print("""
Intrebarea critica: de ce alpha depinde de 20 (tetraedre)
si nu de 12 (vecini)?

RASPUNS POSIBIL:

TETRAEDRELE (20) = VOLUMUL LOCAL
- Reprezinta "spatiul disponibil" pentru interactiune
- O particula poate fi in oricare din cele 20 tetraedre
- Interactiunea = overlap in acelasi tetraedru
- Probabilitate ~ 1/(numar de tetraedre)

VECINII (12) = DIRECTIILE DE MISCARE
- Reprezinta "unde poate merge" particula
- NU reprezinta "unde poate fi" particula
- Relevant pentru PROPAGARE, nu pentru INTERACTIUNE

DECAGOANELE (6) = BUCLE INCHISE
- Reprezinta modurile de "self-interaction"
- O particula poate face bucla pe oricare din cele 6 decagoane
- Relevant pentru STRUCTURA INTERNA, nu pentru interactiune externa

DECI:
- 20 apare in alpha pentru ca masoara VOLUMUL LOCAL
- 12 ar aparea in VITEZA sau PROPAGARE
- 6 apare in UNGHIUL WEINBERG (raport de structuri gauge)
""")

# ============================================================
# MASA DIN PERSPECTIVA DISCRETA
# ============================================================
print("-" * 70)
print("MASA = FRECVENTA DE REZONANTA")
print("-" * 70)

print("""
In mecanica cuantica: E = m*c^2 = hbar*omega

Daca particula e o "rezonanta" pe 600-cell:
  masa ~ frecventa de rezonanta

FRECVENTE POSIBILE:
1. Rezonanta pe decagon (10 pasi):
   omega_decagon = 2*pi / (10 * t_tick) = 2*pi / 10

2. Rezonanta pe tetraedru (bucla locala):
   omega_tet = 2*pi / (4 * t_tick) = 2*pi / 4 = pi/2

3. Rezonanta pe icosaedru (vertex figure):
   omega_ico = 2*pi / (12 * t_tick) = 2*pi / 12 = pi/6

RAPOARTE DE MASE:
  m1/m2 = omega1/omega2

Exemple:
  m_decagon / m_tet = (2*pi/10) / (pi/2) = 2/5
  m_decagon / m_ico = (2*pi/10) / (pi/6) = 6/5
""")

# ============================================================
# ELECTRON VS FOTON
# ============================================================
print("-" * 70)
print("ELECTRON VS FOTON IN MODEL DISCRET")
print("-" * 70)

print("""
FOTONUL:
========
- Unda care se propaga LIBER
- La fiecare tick, sare la un vecin
- Viteza = 1 pixel/tick = c
- Masa = 0 (nu rezoneaza local)
- Energie = hbar * k * c (depinde de lungimea de unda)

ELECTRONUL:
===========
- Pattern care REZONEAZA local
- "Prins" intr-o bucla (decagon sau alta structura)
- Se poate MISCA daca bucla se "rostogoleste"
- Masa = energie de rezonanta / c^2
- Viteza < c (trebuie timp sa "rostogoleasca" bucla)

TRANZITIE FOTON -> ELECTRON (pair creation):
- Fotonul ajunge la un vertex cu geometrie potrivita
- Energia fotonului = 2 * m_e * c^2 (sau mai mult)
- Fotonul se "sparge" in doua bucle (e+ si e-)
- Buclele rezoneaza in directii opuse (conservare impuls)

TRANZITIE ELECTRON -> FOTON (anihilare):
- Doua bucle opuse (e+ si e-) se intalnesc
- Buclele se "anuleaza" (interferenta destructiva)
- Energia se elibereaza ca unde propagante (fotoni)
""")

# ============================================================
# VIZUALIZARE CONCEPTUALA
# ============================================================
print("\n" + "-" * 70)
print("VIZUALIZARE CONCEPTUALA")
print("-" * 70)

print("""
Imagineaza-ti 600-cell ca o "plasa" 4D:

     [V]---[V]---[V]
      |\\   /|\\   /|
      | \\ / | \\ / |
     [V]---[V]---[V]
      | / \\ | / \\ |
      |/   \\|/   \\|
     [V]---[V]---[V]

Fiecare [V] = vertex = pixel
Fiecare --- sau / sau \\ = muchie = cale posibila

FOTON = "*" care sare de la V la V:

  Tick 0:  [*]---[V]---[V]
  Tick 1:  [V]---[*]---[V]
  Tick 2:  [V]---[V]---[*]

ELECTRON = pattern care se "invarte" local:

  Tick 0:  [*]---[V]          [V]---[*]
            |                   |
           [V]                 [V]

  Tick 1:  [V]---[*]          [*]---[V]
            |                   |
           [V]                 [V]

  (bucla care oscileaza intre configuratii)
""")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-055")
print("=" * 70)

print("""
IMAGINEA DISCRETA:
==================

1. SPATIUL = 600-cell cu 120 pixeli (vertexuri)
2. TIMPUL = tick-uri discrete (timp Planck)
3. LUMINA = excitatie care sare 1 pixel/tick
4. PARTICULE = pattern-uri rezonante (bucle stabile)

DE UNDE VIN NUMERELE:
=====================
- 20 (tetraedre) = "volum local" = apare in ALPHA
- 12 (vecini) = "directii de miscare" = apare in PROPAGARE
- 6 (decagoane) = "bucle inchise" = apare in STRUCTURA GAUGE
- 10 (pasi/decagon) = "circumferinta discreta" = apare in FAZA

MECANISMUL UNDA-PARTICULA:
==========================
- UNDA = excitatie delocalizata care se propaga
- PARTICULA = excitatie localizata care rezoneaza
- TRANZITIE = cand unda "se prinde" in geometrie (sau invers)

CE RAMANE DE CALCULAT:
======================
1. Care sunt frecventele exacte de rezonanta?
2. Ce pattern corespunde electronului? Muonului? Quark-ului?
3. Cum apar masele diferite din frecvente diferite?
4. Cum se propaga efectiv un "electron" (rostogolire bucla)?
""")

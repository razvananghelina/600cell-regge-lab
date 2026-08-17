"""
EXP-024: Originea Semnului MINUS
================================
Formula: 1/alpha = 20*phi^4 - 2*pi*alpha

De unde vine MINUS-ul?
"""

from physics_formulas import *

print("=" * 70)
print("EXP-024: DE UNDE VINE SEMNUL MINUS?")
print("=" * 70)

print("""
FORMULA NOASTRA:
================
1/alpha = 20*phi^4 - 2*pi*alpha

Echivalent:
1/alpha + 2*pi*alpha = 20*phi^4

INTREBARE: De ce MINUS (sau echivalent, de ce PLUS in forma a doua)?
""")

print("-" * 70)
print("ANALOGIE DIN FIZICA: SELF-ENERGY")
print("-" * 70)

print("""
In QED, electronul are o "self-energy" din interactiunea cu campul sau propriu.

MASA OBSERVATA vs MASA "GOALA":
  m_obs = m_bare + delta_m (corectie radiativa)

Corectia delta_m vine din:
  - Electronul emite un foton virtual
  - Fotonul se propaga
  - Electronul reabsoarbe fotonul

Aceasta corectie e NEGATIVA pentru energie (leaga electronul de campul sau).

ANALOGIE CU FORMULA NOASTRA:
  1/alpha_geometric = 20*phi^4 (termenul "gol")
  1/alpha_obs = 20*phi^4 - (corectie din self-loop)
""")

print("-" * 70)
print("INTERPRETARE FIZICA A TERMENILOR")
print("-" * 70)

print(f"""
TERMENUL 20*phi^4 = {20*PHI**4:.6f}:
===================================
Provenienta (din EXP-021, 022):
  - 5 = tetraedre per muchie (topologic)
  - 4*phi^4 = 1/lambda_1^2 = propagator din spectrul Laplacian
  - 20*phi^4 = probabilitate de "a nu interactiona" per traversare

TERMENUL 2*pi*alpha = {2*PI*ALPHA:.6f}:
=======================================
Provenienta (din EXP-023):
  - 2*pi = bucla completa (10 pasi x 36 grade)
  - alpha = probabilitatea de interactiune pe pas

COMBINATIA:
  - 20*phi^4 = "impedanta geometrica" a spatiului
  - 2*pi*alpha = "pierdere" din self-interaction in bucla
""")

print("-" * 70)
print("DE CE MINUS? (IPOTEZA 1: ECRANARE)")
print("-" * 70)

print("""
IN QED: ECRANARE (SCREENING)
============================

Sarcina electrica e "ecranata" de perechile virtuale electron-pozitron.

  q_obs = q_bare / (1 + corectie)

La distante mari (energie mica):
  - Vedem sarcina ecranata (mai mica)
  - alpha_obs < alpha_bare

Deci: 1/alpha_obs > 1/alpha_bare
      1/alpha_obs = 1/alpha_bare + (ceva pozitiv)

Dar noi avem:
      1/alpha_obs = 20*phi^4 - 2*pi*alpha

Deci termenul 2*pi*alpha REDUCE 1/alpha,
adica CRESTE alpha!

Asta e OPUS ecranarii normale.
Ar fi mai degraba ANTI-ECRANARE.
""")

print("-" * 70)
print("DE CE MINUS? (IPOTEZA 2: BINDING ENERGY)")
print("-" * 70)

print("""
ENERGIA DE LEGATURA E NEGATIVA
==============================

Cand doua particule sunt legate:
  E_total = E_individual - E_binding

Energia de legatura se SCADE.

Pentru electronul ca "foton prins in bucla":
  - Fara bucla: ar avea "impedanta" 20*phi^4
  - Cu bucla: pierde energie prin self-interaction

Pierderea = 2*pi*alpha

1/alpha = 20*phi^4 - (energie pierduta in bucla)

ACEASTA INTERPRETARE ARE SENS!
""")

print("-" * 70)
print("VERIFICARE NUMERICA")
print("-" * 70)

geometric_term = 20 * PHI**4
loop_term = 2 * PI * ALPHA
result = geometric_term - loop_term

print(f"Termen geometric: 20*phi^4 = {geometric_term:.10f}")
print(f"Termen bucla:     2*pi*alpha = {loop_term:.10f}")
print(f"Diferenta:        {geometric_term} - {loop_term} = {result:.10f}")
print(f"1/alpha CODATA:   {ALPHA_INV:.10f}")
print(f"Eroare:           {abs(result - ALPHA_INV)/ALPHA_INV * 100:.6f}%")

# Cat % din total e termenul de bucla?
fraction = loop_term / geometric_term * 100
print(f"\nTermenul bucla e {fraction:.4f}% din termenul geometric")

print("-" * 70)
print("IPOTEZA 3: INTERFERENTA DESTRUCTIVA")
print("-" * 70)

print("""
UNDE SI FAZE
============

Daca electronul e un "foton in bucla", atunci:
  - Fiecare "parcurgere" a buclei adauga o faza
  - Dupa o bucla completa: faza = 2*pi

Interferenta:
  - Unda directa (geometrica): amplitudine A_geo
  - Unda din self-loop: amplitudine A_loop

Daca fazele sunt OPUSE (pi diferenta):
  A_total = A_geo - A_loop

Probabilitatea:
  P = |A_total|^2 = |A_geo - A_loop|^2

La ordinul 1:
  P ~ P_geo - 2*Re(A_geo* A_loop)

Termenul de interferenta e NEGATIV daca fazele sunt opuse!
""")

print("-" * 70)
print("IPOTEZA 4: RENORMALIZARE")
print("-" * 70)

print("""
IN QFT: RUNNING COUPLING
========================

Constanta de cuplaj "ruleaza" cu scala de energie.

alpha(Q^2) = alpha(0) / (1 - (alpha(0)/3*pi) * ln(Q^2/m_e^2))

La energie Q = m_e (scala electronului):
  alpha(m_e^2) ~ alpha(0) * (1 + corectie)

Deci 1/alpha scade cand adaugam corectii.

ANALOGIE:
  1/alpha_obs = 1/alpha_bare - (corectii radiative)

  Cu bare = 20*phi^4 si corectie = 2*pi*alpha:
  1/alpha = 20*phi^4 - 2*pi*alpha
""")

print("-" * 70)
print("CEA MAI PLAUZIBILA EXPLICATIE")
print("-" * 70)

print("""
CONCLUZIE PROVIZORIE:
=====================

Semnul MINUS vine din faptul ca self-interaction-ul
REDUCE impedanta efectiva (creste cuplajul).

Analogii fizice:
1. Binding energy - energia de legatura e negativa
2. Interferenta - self-loop interfereaza destructiv cu propagarea directa
3. Renormalizare - corectiile radiative reduc 1/alpha

TOATE aceste mecanisme produc un MINUS!

Ce trebuie sa arate Lagrangianul:
- Termen "cinetic" care da 20*phi^4
- Termen de self-interaction care da -2*pi*alpha
- Minusul trebuie sa emerga din structura matematica
""")

print("-" * 70)
print("URMATORUL PAS: CE FEL DE LAGRANGIAN?")
print("-" * 70)

print("""
Din analiza minusului, stim ca Lagrangianul trebuie sa aiba:

L = L_cinetic + L_interaction

unde:
  L_cinetic -> da termenul 20*phi^4 (propagare pe 600-cell)
  L_interaction -> da termenul -2*pi*alpha (self-coupling)

Forma generala pentru un camp scalar pe graf:
  L = (1/2) * sum_ij psi_i * L_ij * psi_j - (lambda/4!) * sum_i psi_i^4

Trebuie sa:
1. Adaptam pentru 600-cell
2. Calculam contributiile
3. Verificam ca obtinem formula corecta
""")

print("\n" + "=" * 70)
print("REZUMAT EXP-024")
print("=" * 70)

print(f"""
SEMNUL MINUS E NATURAL!
=======================

Vine din faptul ca self-interaction-ul:
- Leaga particula de propriul camp (binding energy < 0)
- Produce interferenta destructiva
- Reduce "impedanta" geometrica

NUMERIC:
  20*phi^4 = {geometric_term:.6f} (termen geometric)
  2*pi*alpha = {loop_term:.6f} (corectie self-loop)
  Diferenta = {result:.6f}
  1/alpha = {ALPHA_INV:.6f}

Termenul de bucla e {fraction:.3f}% din total - o CORECTIE mica,
exact cum ne-am astepta de la o interactiune perturbativa!

URMATORUL PAS:
Scriem Lagrangianul explicit si verificam ca produce
acesti termeni cu semnele corecte.
""")

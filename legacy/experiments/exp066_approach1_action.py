"""
EXP-066: ABORDAREA 1 - Derivarea din Actiune (Lagrangian)
=========================================================

Scopul: Sa derivam ecuatia 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
pornind de la un Lagrangian pe reteaua 600-cell.

Metoda:
1. Definim campul pe varfurile 600-cell
2. Scriem Lagrangianul discret
3. Minimizam actiunea
4. Cautam conditia pentru soliton stabil
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-066: DERIVARE DIN ACTIUNE PE 600-CELL")
print("=" * 70)

# ============================================================
# PASUL 1: SETUP - CAMPUL PE RETEA
# ============================================================
print("\n" + "-" * 70)
print("PASUL 1: CAMPUL PE RETEAUA 600-CELL")
print("-" * 70)

print("""
DEFINITII:
==========
- psi_i = valoarea campului la varful i (i = 1, ..., 120)
- A_ij = conexiunea gauge pe muchia (i,j)
- L_ij = Laplacianul discret (matrice de adiacenta)

LAGRANGIANUL PE RETEA:
======================
L = T - V

Unde:
T = (1/2) * sum_i |dpsi_i/dt|^2       (energie cinetica)

V = V_gradient + V_potential + V_gauge

V_gradient = (1/2) * sum_<ij> |psi_i - e^(i*A_ij)*psi_j|^2
V_potential = sum_i [ (m^2/2)|psi_i|^2 + (lambda/4)|psi_i|^4 ]
V_gauge = (1/4) * sum_plachete F_munu^2

SIMPLIFICARE (pentru inceput):
==============================
Ignoram campul gauge (A_ij = 0) si termenul cinetic.
Cautam configuratii STATICE.
""")

# ============================================================
# PASUL 2: CONFIGURATIA SOLITON
# ============================================================
print("\n" + "-" * 70)
print("PASUL 2: ANSATZ PENTRU SOLITON")
print("-" * 70)

print("""
ANSATZ: Solitonul e localizat in jurul unui varf central.

psi_i = psi_0 * f(d_i)

Unde:
- psi_0 = amplitudinea centrala
- d_i = distanta (numar de pasi) de la varful central
- f(d) = functie de profil (descrescatoare)

PROFILE POSIBILE:
=================
1. Exponential: f(d) = exp(-d/xi)
2. Gaussian: f(d) = exp(-d^2/(2*xi^2))
3. Algebraic: f(d) = 1/(1 + (d/xi)^2)

Vom folosi profilul EXPONENTIAL pentru inceput.
""")

# ============================================================
# PASUL 3: ENERGIA CONFIGURATIEI
# ============================================================
print("\n" + "-" * 70)
print("PASUL 3: CALCULUL ENERGIEI")
print("-" * 70)

# Construim o retea simplificata (prima coaja = icosaedru)
# Varful central + 12 vecini + 20 vecini la distanta 2 + ...

# Numarul de varfuri la fiecare distanta de centru
# Pentru 600-cell: 1, 12, 20, 12, 30, 12, 20, 12, 1 (schematic)
shells = {
    0: 1,    # centru
    1: 12,   # prima coaja (icosaedru)
    2: 20,   # a doua coaja
    3: 12,   # a treia coaja
    4: 30,   # a patra coaja
    5: 12,   # ...
}

print("Structura de cochilii a 600-cell:")
for d, n in shells.items():
    print(f"  Distanta {d}: {n} varfuri")

# Energia gradientului pentru profil exponential
def gradient_energy(psi0, xi, z=12):
    """
    Energia de gradient pentru soliton exponential.

    E_grad = (1/2) * sum_<ij> |psi_i - psi_j|^2

    Pentru profil psi = psi0 * exp(-d/xi):
    La fiecare muchie de la d la d+1:
    |psi_d - psi_{d+1}|^2 = psi0^2 * exp(-2d/xi) * (1 - exp(-1/xi))^2
    """
    E = 0
    # Muchii de la centru (d=0) la prima coaja (d=1)
    # Sunt 12 astfel de muchii
    delta = 1 - np.exp(-1/xi)

    for d in range(5):
        n_edges = shells.get(d, 0) * z / 2  # aproximativ
        psi_d = psi0 * np.exp(-d/xi)
        contribution = n_edges * psi_d**2 * delta**2
        E += contribution

    return E

# Energia potentialului
def potential_energy(psi0, xi, m2, lam):
    """
    Energia potentialului V = (m^2/2)|psi|^2 + (lambda/4)|psi|^4

    Suma pe toate varfurile.
    """
    E = 0
    for d, n in shells.items():
        psi_d = psi0 * np.exp(-d/xi)
        E += n * (m2/2 * psi_d**2 + lam/4 * psi_d**4)
    return E

# Energia totala
def total_energy(psi0, xi, m2, lam, z=12):
    return gradient_energy(psi0, xi, z) + potential_energy(psi0, xi, m2, lam)

print("""
ENERGIA TOTALA:
===============
E = E_gradient + E_potential

E_gradient ~ psi0^2 * (1 - e^(-1/xi))^2 * N_muchii
E_potential ~ psi0^2 * m^2 * N_varfuri + psi0^4 * lambda * N_varfuri

Pentru soliton STABIL, minimizam E in raport cu psi0 si xi.
""")

# ============================================================
# PASUL 4: CONDITIA DE MINIMIZARE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 4: CONDITIA DE STATIONARITATE")
print("-" * 70)

print("""
Minimizam E(psi0, xi):

dE/d(psi0) = 0  =>  ecuatia pentru amplitudine
dE/d(xi) = 0    =>  ecuatia pentru marime

FORMA GENERALA:
===============
dE/d(psi0) = A * psi0 + B * psi0^3 = 0

Solutii:
1. psi0 = 0 (vacuum trivial)
2. psi0^2 = -A/B (soliton netrivial, daca A/B < 0)

Pentru soliton stabil, trebuie A < 0 (masa negativa = tachion)
si B > 0 (autointeractiune pozitiva).
""")

# ============================================================
# PASUL 5: INTRODUCEREA CONSTANTEI DE CUPLAJ
# ============================================================
print("\n" + "-" * 70)
print("PASUL 5: CONSTANTA DE CUPLAJ ALPHA")
print("-" * 70)

print("""
In teoria gauge U(1), constanta de cuplaj e apare in:

D_mu = partial_mu - i*e*A_mu

Lagrangianul gauge:
L_gauge = -(1/4) * F_munu * F^munu

Constanta de structura fina:
alpha = e^2 / (4*pi) [in unitati Gaussian/naturale]

PROBLEMA: Cum apare ECUATIA NOASTRA?
====================================

Ecuatia 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0 sugereaza:

- Termenul 2*pi*alpha^2: vine din energia gauge (F^2)?
- Termenul 20*phi^4*alpha: vine din cuplajul camp-gauge?
- Termenul 1: normalizare sau energie de vacuum?

HAI SA INCERCAM O DERIVARE:
===========================
""")

# ============================================================
# PASUL 6: DERIVARE EXPLICITA (INCERCARE)
# ============================================================
print("\n" + "-" * 70)
print("PASUL 6: INCERCARE DE DERIVARE")
print("-" * 70)

print("""
IPOTEZA: Energia solitonului are forma:

E = E_grad + E_gauge + E_potential

Cu:
E_grad = (1/2) * integral |D*psi|^2 ~ psi0^2 * (ceva geometric)
E_gauge = (1/4) * integral F^2 ~ alpha * (flux magnetic)^2
E_potential = integral V(|psi|) ~ psi0^4 * lambda

CONDITIA DE CUANTIZARE A FLUXULUI:
==================================
Pentru soliton topologic (vortex), fluxul magnetic e cuantizat:

Phi = n * (2*pi/e) = n * (2*pi/sqrt(4*pi*alpha))

Energia gauge:
E_gauge ~ alpha * Phi^2 ~ alpha * (2*pi)^2 / alpha = 4*pi^2 / alpha

Hmm, asta nu da exact ce vrem...

ALTA INCERCARE - Energia de legatura:
=====================================
Energia totala a solitonului (in unitati de m_Planck):

E_total = m_e = f(alpha, geometrie)

Din EXP-062, stim:
m_e / m_Planck = alpha^(4*phi^2)

Deci:
ln(m_e/m_P) = 4*phi^2 * ln(alpha)

Derivand in raport cu alpha:
d(ln E)/d(alpha) = 4*phi^2 / alpha

Asta e o relatie diferentiala, nu algebrica...
""")

# ============================================================
# PASUL 7: ABORDARE ALTERNATIVA - AUTO-CONSISTENTA
# ============================================================
print("\n" + "-" * 70)
print("PASUL 7: ECUATIA DE AUTO-CONSISTENTA")
print("-" * 70)

print("""
IDEE: Campul electromagnetic creeaza solitonul,
      iar solitonul sustine campul electromagnetic.

Ecuatia de auto-consistenta:
============================

1. Campul electric al solitonului:
   E ~ e / r^2 = sqrt(4*pi*alpha) / r^2

2. Energia campului:
   U_field = (1/2) * integral E^2 d^3r
           = (1/2) * (4*pi*alpha) * integral (1/r^4) * 4*pi*r^2 dr
           = 2*pi * alpha * integral dr/r^2

   Cu cutoff la r_min (marimea solitonului):
   U_field ~ 2*pi * alpha / r_min

3. Energia de masa a solitonului:
   U_mass = m_e = (ceva) * m_Planck * alpha^(4*phi^2)

4. CONDITIA: Energia totala e MINIMA

   d(U_field + U_mass) / d(r_min) = 0

Asta da o relatie intre alpha si r_min, dar nu exact ecuatia noastra...
""")

# ============================================================
# PASUL 8: CE AM GASIT SI CE LIPSESTE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 8: REZULTATE SI LIMITARI")
print("-" * 70)

print("""
CE AM REUSIT:
=============
1. Am scris Lagrangianul pe retea (forma generala)
2. Am definit ansatz-ul pentru soliton
3. Am identificat termenii de energie

CE NU AM REUSIT:
================
1. Nu am derivat EXACT ecuatia 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
2. Nu am aratat de unde vin coeficientii 2*pi si 20*phi^4
3. Nu am conectat explicit cu geometria 600-cell

POSIBILE MOTIVE:
================
1. Ecuatia nu vine dintr-un Lagrangian simplu
2. Trebuie incluse efecte CUANTICE (nu doar clasice)
3. Geometria 600-cell trebuie folosita mai explicit
4. Poate ecuatia vine din ALTA abordare (spectru, Hopf)

CONCLUZIE ABORDAREA 1:
======================
Abordarea clasica din Lagrangian nu a dat direct ecuatia.
Trebuie sa incercam alta metoda.
""")

# ============================================================
# VIZUALIZARE
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 1. Profilul solitonului
ax = axes[0, 0]
ax.set_title('Profilul Solitonului pe Retea', fontweight='bold')
d = np.arange(0, 6)
for xi in [0.5, 1, 2, 3]:
    psi = np.exp(-d/xi)
    ax.plot(d, psi, 'o-', label=f'xi = {xi}', markersize=10)
ax.set_xlabel('Distanta de la centru (pasi)')
ax.set_ylabel('psi / psi_0')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xticks(d)

# 2. Energia vs marime
ax = axes[0, 1]
ax.set_title('Energia Solitonului vs Marime', fontweight='bold')
xi_range = np.linspace(0.3, 5, 100)
psi0 = 1
m2 = -1  # masa negativa pentru soliton
lam = 0.5

E_grad = np.array([gradient_energy(psi0, xi) for xi in xi_range])
E_pot = np.array([potential_energy(psi0, xi, m2, lam) for xi in xi_range])
E_tot = E_grad + E_pot

ax.plot(xi_range, E_grad, 'b-', label='E_gradient')
ax.plot(xi_range, E_pot, 'r-', label='E_potential')
ax.plot(xi_range, E_tot, 'k-', linewidth=2, label='E_total')
ax.axhline(y=0, color='gray', linestyle='--')

# Gasim minimul
idx_min = np.argmin(E_tot)
ax.axvline(x=xi_range[idx_min], color='green', linestyle='--',
           label=f'xi_min = {xi_range[idx_min]:.2f}')

ax.set_xlabel('xi (marimea solitonului)')
ax.set_ylabel('Energie')
ax.legend()
ax.grid(True, alpha=0.3)

# 3. Termenii ecuatiei
ax = axes[1, 0]
ax.set_title('Termenii Ecuatiei vs Alpha', fontweight='bold')
alpha_range = np.linspace(0.001, 0.05, 100)

term1 = 2 * np.pi * alpha_range**2
term2 = 20 * PHI**4 * alpha_range
term3 = np.ones_like(alpha_range)

ax.plot(alpha_range, term1, 'b-', label=r'$2\pi\alpha^2$')
ax.plot(alpha_range, term2, 'r-', label=r'$20\phi^4\alpha$')
ax.plot(alpha_range, term3, 'g-', label='1')
ax.axvline(x=ALPHA, color='purple', linestyle='--', label=f'alpha = 1/137')

ax.set_xlabel('alpha')
ax.set_ylabel('Valoare')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 0.05)

# 4. Diagrama fluxului de lucru
ax = axes[1, 1]
ax.axis('off')
ax.set_title('Fluxul de Derivare (Abordarea 1)', fontweight='bold')

diagram = """
    +-----------------+
    |   LAGRANGIAN    |
    |  L = T - V      |
    +--------+--------+
             |
             v
    +--------+--------+
    |  ACTIUNE        |
    |  S = integral L |
    +--------+--------+
             |
             v
    +--------+--------+
    |  EULER-LAGRANGE |
    |  dS/d(psi) = 0  |
    +--------+--------+
             |
             v
    +--------+--------+
    |  ANSATZ SOLITON |
    |  psi = f(r)     |
    +--------+--------+
             |
             v
    +--------+--------+
    |  MINIMIZARE     |
    |  dE/d(param)=0  |
    +--------+--------+
             |
             v
    +--------+--------+
    |  ECUATIA?       |
    |  ???            |  <-- Nu am ajuns aici
    +--------+--------+
"""
ax.text(0.1, 0.95, diagram, transform=ax.transAxes,
       fontsize=10, verticalalignment='top',
       fontfamily='monospace',
       bbox=dict(boxstyle='round', facecolor='lightyellow'))

ax.text(0.6, 0.3,
       "PROBLEMA:\n"
       "Nu am reusit sa derivam\n"
       "ecuatia din Lagrangian.\n\n"
       "Poate trebuie:\n"
       "- Efecte cuantice\n"
       "- Geometrie explicita\n"
       "- Alta abordare",
       fontsize=11,
       bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))

plt.tight_layout()
plt.savefig('approach1_action.png', dpi=150, bbox_inches='tight')
print("\nSalvat: approach1_action.png")
plt.show()

print("\n" + "=" * 70)
print("CONCLUZIE: Abordarea 1 nu a dat rezultatul dorit.")
print("Trecem la Abordarea 2: Spectrul Laplacianului.")
print("=" * 70)

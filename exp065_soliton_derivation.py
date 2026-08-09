"""
EXP-065: Cum ar putea aparea un soliton din ecuatia alfa?
=========================================================

Ecuatia: 2*pi*alpha^2 - (20*phi^4)*alpha + 1 = 0

Intrebari:
1. Ce camp oscileaza?
2. Ce creeaza neliniaritatea care stabilizeaza solitonul?
3. De unde vine probabilitatea cuantica?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-065: DERIVAREA SOLITONULUI DIN GEOMETRIA 600-CELL")
print("=" * 70)

# ============================================================
# PASUL 1: ANALIZA ECUATIEI
# ============================================================
print("\n" + "-" * 70)
print("PASUL 1: ANALIZA ECUATIEI FUNDAMENTALE")
print("-" * 70)

print("""
Ecuatia noastra:
    2*pi*alpha^2 - (20*phi^4)*alpha + 1 = 0

Aceasta e o ecuatie CUADRATICA in alpha.
Are forma: a*x^2 + b*x + c = 0

Unde:
    a = 2*pi       (geometrie circulara - bucle?)
    b = -20*phi^4  (20 tetraedre, phi din geometrie)
    c = 1          (normalizare)

INTREBARE CHEIE: De unde ar putea veni aceasta ecuatie?
""")

a = 2 * np.pi
b = -20 * PHI**4
c = 1

# Solutiile
disc = b**2 - 4*a*c
alpha1 = (-b + np.sqrt(disc)) / (2*a)
alpha2 = (-b - np.sqrt(disc)) / (2*a)

print(f"Coeficienti: a = 2*pi = {a:.4f}")
print(f"             b = -20*phi^4 = {b:.4f}")
print(f"             c = 1")
print(f"\nDiscriminant: {disc:.4f}")
print(f"Solutii: alpha_1 = {alpha1:.6f} = 1/{1/alpha1:.2f}")
print(f"         alpha_2 = {alpha2:.6f} = 1/{1/alpha2:.2f}")
print(f"\nalpha fizic = {ALPHA:.6f} = 1/137.036")

# ============================================================
# PASUL 2: POSIBILE ORIGINI FIZICE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 2: POSIBILE ORIGINI FIZICE ALE ECUATIEI")
print("-" * 70)

print("""
IPOTEZA A: MINIMIZAREA ACTIUNII
===============================
Daca avem un Lagrangian de forma:
    L = (1/2)(dpsi/dt)^2 - V(psi)

Si potential:
    V(psi) = (1/2)*m^2*psi^2 + (lambda/4)*psi^4

Ecuatia de miscare (Klein-Gordon neliniara):
    d^2(psi)/dt^2 - nabla^2(psi) + m^2*psi + lambda*psi^3 = 0

Solitonul apare cand:
    - Dispersia (nabla^2) e balansata de neliniaritate (psi^3)
    - Rezulta o configuratie stabila, localizata

IPOTEZA B: CONDITIE DE AUTOCONSISTENTA
======================================
Campul creeaza geometria care il sustine:
    1. Campul psi deformeaza spatiul (curbura)
    2. Curbura afecteaza propagarea lui psi
    3. Starea stationara: geometria si campul sunt consistente

Ecuatia ar putea fi conditia de INCHIDERE A BUCLEI.

IPOTEZA C: REZONANTA PE 600-CELL
================================
Undele stationare pe 600-cell au frecvente discrete:
    omega_n = sqrt(lambda_n)

Solitonul = superozitie de moduri care interfereaza constructiv
in centru si destructiv la distanta.
""")

# ============================================================
# PASUL 3: MODEL SIMPLIFICAT - SOLITON 1D
# ============================================================
print("\n" + "-" * 70)
print("PASUL 3: MODEL SIMPLIFICAT - SOLITON 1D")
print("-" * 70)

print("""
Ecuatia solitonului (sine-Gordon):
    d^2(phi)/dx^2 - d^2(phi)/dt^2 = sin(phi)

Solutie soliton (kink):
    phi(x,t) = 4 * arctan(exp((x - v*t) / sqrt(1-v^2)))

In contextul nostru, phi ar putea fi FAZA campului
pe reteaua 600-cell.
""")

# Simulam solitonul
x = np.linspace(-10, 10, 500)
v = 0.5  # viteza

def soliton_kink(x, v=0):
    gamma = 1 / np.sqrt(1 - v**2) if v < 1 else 1
    return 4 * np.arctan(np.exp(gamma * x))

def soliton_envelope(x, width=1):
    """Anvelopa gaussiana a solitonului."""
    return np.exp(-x**2 / (2 * width**2))

# ============================================================
# PASUL 4: CONEXIUNE CU GEOMETRIA 600-CELL
# ============================================================
print("\n" + "-" * 70)
print("PASUL 4: CONEXIUNE CU GEOMETRIA 600-CELL")
print("-" * 70)

print(f"""
Din ecuatia noastra, putem extrage:

1. TERMENUL 2*pi:
   - Bucla inchisa (fotonul "prins")
   - Integrala pe un cerc = faza totala
   - Conditia de cuantizare: faza = 2*pi*n

2. TERMENUL 20*phi^4:
   - 20 = tetraedre in jurul fiecarui vertex
   - phi^4 = propagator in 4D? (fiecare phi = salt pe o muchie)
   - phi^4 ~ 6.85 = amplitudine de salt

3. TERMENUL 1:
   - Normalizare (probabilitatea totala = 1)
   - Sau: energia de repaus (m*c^2 = 1 in unitati naturale)

INTERPRETARE POSIBILA:
======================
Electronul = foton care se propaga in BUCLA inchisa
pe reteaua 600-cell, cu:
- 2*pi = faza acumulata pe bucla
- 20*phi^4 = amplitudine de propagare prin 20 tetraedre
- 1 = conditia de normalizare/inchidere

Ecuatia = conditia ca bucla sa fie STABILA (soliton).
""")

# ============================================================
# PASUL 5: DERIVARE DIN PRINCIPIUL VARIATIONAL
# ============================================================
print("\n" + "-" * 70)
print("PASUL 5: INCERCARE DE DERIVARE RIGUROASA")
print("-" * 70)

print("""
IPOTEZA: Actiunea pe 600-cell

S = integral(L * d^4x)

Unde Lagrangianul e:
L = (1/2)|D_mu * psi|^2 - (1/4)*F_munu*F^munu - V(|psi|)

Cu:
- psi = camp scalar complex (functia de unda)
- D_mu = derivata covarianta (include campul gauge A_mu)
- F_munu = tensorul electromagnetic
- V = potential

Pe reteaua 600-cell:
- D_mu -> diferente finite pe muchii
- Integrala -> suma pe celule
- Simetria I_h impune constrangeri

PROBLEMA: Cum obtinem ecuatia noastra din acest Lagrangian?

Raspuns posibil: Ecuatia apare din CONDITIA DE STATIONARITATE
a actiunii pentru configuratia soliton.
""")

# ============================================================
# PASUL 6: PROBABILITATEA CUANTICA
# ============================================================
print("\n" + "-" * 70)
print("PASUL 6: DE UNDE VINE PROBABILITATEA CUANTICA?")
print("-" * 70)

print("""
In mecanica cuantica standard:
    |psi(x)|^2 = probabilitatea de a gasi particula la x

In teoria noastra:
    Solitonul nu e localizat intr-un punct!
    E o configuratie EXTINSA pe reteaua 600-cell.

POSIBILA INTERPRETARE:
======================
1. "Pozitia" electronului = centrul de masa al solitonului
2. |psi|^2 = densitatea de energie/informatie a solitonului
3. Masurarea = interactiune care "colecteaza" solitonul

Nedeterminismul cuantic ar putea veni din:
- Fluctuatii ale retelei (la scala Planck)
- Interactiuni cu alte solitoni
- Limitarea fundamentala a informatiei accesibile

CHEIE: Solitonul are o STRUCTURA INTERNA
       Nu e un punct, dar la distante >> l_Planck
       PARE un punct (limita punctiforma).
""")

# ============================================================
# VIZUALIZARE
# ============================================================
print("\n" + "-" * 70)
print("VIZUALIZARE")
print("-" * 70)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# 1. Soliton kink
ax = axes[0, 0]
ax.set_title('Soliton Kink (Sine-Gordon)', fontweight='bold')
for v in [0, 0.3, 0.6, 0.9]:
    y = soliton_kink(x, v)
    ax.plot(x, y, label=f'v = {v}')
ax.set_xlabel('x')
ax.set_ylabel('phi(x)')
ax.legend()
ax.grid(True, alpha=0.3)
ax.axhline(y=np.pi, color='r', linestyle='--', alpha=0.5)
ax.text(5, np.pi+0.3, 'phi = pi', color='r')

# 2. Densitate de probabilitate
ax = axes[0, 1]
ax.set_title('Densitate de Probabilitate |psi|^2', fontweight='bold')
widths = [0.5, 1, 2, 4]
for w in widths:
    psi_sq = soliton_envelope(x, w)**2
    psi_sq = psi_sq / np.sum(psi_sq) / (x[1]-x[0])  # Normalizam
    ax.plot(x, psi_sq, label=f'width = {w}')
ax.set_xlabel('x')
ax.set_ylabel('|psi(x)|^2')
ax.legend()
ax.grid(True, alpha=0.3)

# 3. Faza pe bucla
ax = axes[0, 2]
ax.set_title('Faza pe Bucla Inchisa', fontweight='bold')
theta = np.linspace(0, 2*np.pi, 100)
r = 1 + 0.2 * np.sin(5*theta)  # Perturbatie pentagonala
ax.plot(r * np.cos(theta), r * np.sin(theta), 'b-', linewidth=2)
# Sageti pentru faza
for i in range(5):
    angle = i * 2*np.pi / 5
    ax.annotate('', xy=(1.3*np.cos(angle+0.3), 1.3*np.sin(angle+0.3)),
               xytext=(1.3*np.cos(angle), 1.3*np.sin(angle)),
               arrowprops=dict(arrowstyle='->', color='red', lw=2))
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect('equal')
ax.text(0, 0, r'$\oint d\phi = 2\pi n$', fontsize=14, ha='center')
ax.axis('off')

# 4. Potential V(psi)
ax = axes[1, 0]
ax.set_title('Potential V(psi) - "Mexican Hat"', fontweight='bold')
psi = np.linspace(-2, 2, 200)
# Potential phi^4 cu simetrie rupta
mu_sq = -1
lam = 0.5
V = -mu_sq/2 * psi**2 + lam/4 * psi**4
ax.plot(psi, V, 'b-', linewidth=2)
ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
ax.set_xlabel(r'$\psi$')
ax.set_ylabel(r'$V(\psi)$')
ax.grid(True, alpha=0.3)
# Marcam minimele
psi_min = np.sqrt(-mu_sq/lam)
ax.axvline(x=psi_min, color='r', linestyle='--', alpha=0.5)
ax.axvline(x=-psi_min, color='r', linestyle='--', alpha=0.5)
ax.text(psi_min, -0.3, 'vacuum', ha='center', color='r')

# 5. Structura solitonului pe retea
ax = axes[1, 1]
ax.set_title('Soliton pe Reteaua 600-cell\n(Vedere 2D simplificata)', fontweight='bold')

# Desenam o retea hexagonala
for i in range(-5, 6):
    for j in range(-5, 6):
        x_pos = i + 0.5*j
        y_pos = j * np.sqrt(3)/2
        # Intensitatea = distanta de la centru
        r = np.sqrt(x_pos**2 + y_pos**2)
        intensity = np.exp(-r**2 / 4)
        color = plt.cm.hot(intensity)
        circle = plt.Circle((x_pos, y_pos), 0.2, color=color, alpha=0.8)
        ax.add_patch(circle)

ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')

# Contur pentru "electron"
theta = np.linspace(0, 2*np.pi, 100)
ax.plot(1.5*np.cos(theta), 1.5*np.sin(theta), 'w--', linewidth=2, alpha=0.7)
ax.text(0, -4, '"Marimea" electronului ~ cateva l_Planck', ha='center', color='white',
       fontsize=10, bbox=dict(facecolor='black', alpha=0.5))

# 6. Ecuatia si interpretarea
ax = axes[1, 2]
ax.axis('off')
ax.set_title('Interpretarea Ecuatiei', fontweight='bold')

text = r"""
ECUATIA FUNDAMENTALA:
$2\pi\alpha^2 - 20\phi^4\alpha + 1 = 0$

INTERPRETARE CA SOLITON:
========================

$2\pi$ = faza totala pe bucla inchisa
       (fotonul face un tur complet)

$20\phi^4$ = amplitudine de propagare
           prin cele 20 tetraedre
           ($\phi^4$ = 4 salturi pe muchii)

$1$ = conditia de normalizare
    (probabilitate totala = 1)

REZULTAT:
=========
$\alpha = \frac{20\phi^4 - \sqrt{(20\phi^4)^2 - 8\pi}}{4\pi}$

$\alpha \approx 1/137.036$

Solitonul e STABIL cand aceasta
ecuatie e satisfacuta!
"""

ax.text(0.1, 0.95, text, transform=ax.transAxes,
       fontsize=11, verticalalignment='top',
       fontfamily='monospace',
       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()
plt.savefig('soliton_interpretation.png', dpi=150, bbox_inches='tight')
print("\nSalvat: soliton_interpretation.png")
plt.show()

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-065")
print("=" * 70)

print("""
CE AM INTELES:
==============
1. Ecuatia 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0 ar putea fi
   CONDITIA DE STABILITATE pentru un soliton pe 600-cell.

2. Solitonul = foton "prins" intr-o bucla inchisa,
   propagandu-se prin cele 20 de tetraedre.

3. Alpha = constanta de cuplaj care face bucla STABILA.

CE NU STIM INCA:
================
1. Cum se deriva EXACT aceasta ecuatie din un Lagrangian?
2. Care e campul fundamental (scalar, spinor, gauge)?
3. Cum apare masa electronului din aceasta configuratie?
4. De unde vine nedeterminismul cuantic?

PASI URMATORI:
==============
1. Incerca sa derivam ecuatia din actiunea pe 600-cell
2. Verifica daca spectrul Laplacianului da termenii corecti
3. Explora conexiunea cu fibrarea Hopf (bucle inchise pe S^3)
4. Compara cu alte teorii de solitoni (Skyrme, 't Hooft-Polyakov)

NOTA IMPORTANTA:
================
Aceasta e INCA o IPOTEZA. Nu am derivat riguros ecuatia.
Conexiunea geometrica e intriganta, dar trebuie verificata!
""")

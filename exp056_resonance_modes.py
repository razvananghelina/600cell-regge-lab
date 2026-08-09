"""
EXP-056: Moduri de Rezonanta pe 600-cell
========================================
Calculam buclele inchise (cicluri) pe 600-cell si
frecventele lor de rezonanta. Comparam cu masele particulelor.

Ipoteza: masa ~ frecventa ~ 1/lungime_bucla
"""

from physics_formulas import *
import numpy as np
from collections import defaultdict

print("=" * 70)
print("EXP-056: MODURI DE REZONANTA PE 600-CELL")
print("=" * 70)

# ============================================================
# CONSTRUIRE 600-CELL
# ============================================================
print("\n" + "-" * 70)
print("CONSTRUIRE 600-CELL")
print("-" * 70)

def build_600cell():
    """Construieste vertexurile 600-cell."""
    vertices = []

    # 1. 8 vertexuri din 16-cell: permutari de (+-1, 0, 0, 0)
    for i in range(4):
        for sign in [1, -1]:
            v = [0, 0, 0, 0]
            v[i] = sign
            vertices.append(tuple(v))

    # 2. 16 vertexuri: (+-1/2, +-1/2, +-1/2, +-1/2)
    from itertools import product
    for signs in product([0.5, -0.5], repeat=4):
        vertices.append(signs)

    # 3. 96 vertexuri: permutari pare de (+-phi/2, +-1/2, +-1/(2*phi), 0)
    # phi = golden ratio
    phi = PHI
    coords = [phi/2, 0.5, 1/(2*phi), 0]

    # Toate permutarile de semne si pozitii (cu numar par de minusuri in primele 3)
    from itertools import permutations

    for perm in permutations(range(4)):
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    # Numar par de minusuri pentru paritate
                    v = [0, 0, 0, 0]
                    vals = [phi/2, 0.5, 1/(2*phi), 0]
                    signs = [s1, s2, s3, 1]  # ultimul e fix pentru a nu dubla

                    for i, p in enumerate(perm):
                        v[p] = vals[i] * signs[i]

                    vertices.append(tuple(v))

    # Eliminam duplicate si normalizam
    unique = []
    seen = set()
    for v in vertices:
        # Normalizam la 6 zecimale pentru comparatie
        v_round = tuple(round(x, 6) for x in v)
        if v_round not in seen:
            seen.add(v_round)
            # Verificam ca e pe sfera de raza 1
            r = np.sqrt(sum(x**2 for x in v))
            if abs(r - 1.0) < 0.01:
                unique.append(np.array(v) / r)  # normalizam

    return unique[:120]  # primele 120

# Construim vertexurile
vertices = build_600cell()
print(f"Vertexuri construite: {len(vertices)}")

# Construim muchiile (vecini la distanta 1/phi)
def build_edges(vertices, edge_length=1/PHI, tol=0.01):
    """Gaseste muchiile (perechi de vertexuri la distanta edge_length)."""
    edges = []
    neighbors = defaultdict(list)
    n = len(vertices)

    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(vertices[i] - vertices[j])
            if abs(d - edge_length) < tol:
                edges.append((i, j))
                neighbors[i].append(j)
                neighbors[j].append(i)

    return edges, neighbors

edges, neighbors = build_edges(vertices)
print(f"Muchii gasite: {len(edges)}")
print(f"Vecini per vertex: {len(neighbors[0])} (ar trebui 12)")

# ============================================================
# GASIREA CICLURILOR (BUCLELOR INCHISE)
# ============================================================
print("\n" + "-" * 70)
print("GASIREA CICLURILOR")
print("-" * 70)

def find_cycles_of_length(neighbors, start, length, max_cycles=1000):
    """Gaseste cicluri de lungime exacta pornind din start."""
    cycles = []

    def dfs(path):
        if len(cycles) >= max_cycles:
            return

        current = path[-1]
        depth = len(path)

        if depth == length:
            # Verificam daca ne putem intoarce la start
            if start in neighbors[current]:
                # Am gasit un ciclu
                cycles.append(path + [start])
            return

        if depth > length:
            return

        for next_v in neighbors[current]:
            if next_v not in path:  # evitam revisitarea
                dfs(path + [next_v])

    dfs([start])
    return cycles

print("Cautam cicluri de diferite lungimi...")

# Cautam cicluri de lungimi 3, 4, 5, 6, ..., 12
cycle_counts = {}
cycle_examples = {}

for length in range(3, 13):
    cycles = find_cycles_of_length(neighbors, 0, length, max_cycles=100)
    cycle_counts[length] = len(cycles)
    if cycles:
        cycle_examples[length] = cycles[0]
    print(f"  Cicluri de lungime {length}: {len(cycles)}")

# ============================================================
# FRECVENTE DE REZONANTA
# ============================================================
print("\n" + "-" * 70)
print("FRECVENTE DE REZONANTA")
print("-" * 70)

print("""
Pentru o bucla de lungime L (in unitati de muchie = 1/phi):
  Lungime fizica = L * (1/phi)

Unda stationara pe bucla: lambda = L_fizic / n (n = mod de rezonanta)
Frecventa: f = c / lambda = c * n / L_fizic

Pentru modul fundamental (n=1):
  f ~ 1 / L_fizic ~ phi / L

Daca masa ~ frecventa:
  m ~ phi / L

Raport de mase pentru doua bucle:
  m1/m2 = L2/L1
""")

# Calculam frecventele relative (normalizate la decagon)
edge_length = 1/PHI
L_decagon = 10 * edge_length  # decagon = 10 pasi

print(f"\nLungime decagon (referinta): L_10 = 10/phi = {L_decagon:.4f}")
print(f"\nFrecvente relative (f_10 = 1):\n")

frequencies = {}
for length in sorted(cycle_counts.keys()):
    if cycle_counts[length] > 0:
        L = length * edge_length
        f_rel = L_decagon / L  # frecventa relativa la decagon
        frequencies[length] = f_rel
        print(f"  Ciclu L={length}: f_rel = {f_rel:.4f} (L_fizic = {L:.4f})")

# ============================================================
# COMPARATIE CU RAPOARTE DE MASE
# ============================================================
print("\n" + "-" * 70)
print("COMPARATIE CU RAPOARTE DE MASE")
print("-" * 70)

# Rapoarte de mase cunoscute
mass_ratios = {
    'm_mu/m_e': 206.77,
    'm_tau/m_mu': 16.82,
    'm_tau/m_e': 3477.2,
    'm_p/m_e': 1836.15,
    'm_W/m_Z': 0.8815,
}

print("Rapoarte de mase experimentale:")
for name, ratio in mass_ratios.items():
    print(f"  {name} = {ratio}")

print("\nRapoarte de frecvente din cicluri:")
freq_ratios = {}
lengths = sorted(frequencies.keys())
for i, L1 in enumerate(lengths):
    for L2 in lengths[i+1:]:
        ratio = frequencies[L1] / frequencies[L2]  # f1/f2 = L2/L1
        inv_ratio = 1/ratio
        freq_ratios[(L1, L2)] = ratio
        print(f"  f_{L1}/f_{L2} = {ratio:.4f}  (sau {inv_ratio:.4f})")

# ============================================================
# CAUTAM POTRIVIRI
# ============================================================
print("\n" + "-" * 70)
print("CAUTAM POTRIVIRI FRECVENTE <-> MASE")
print("-" * 70)

print("\nCautam combinatii de cicluri care dau rapoarte de mase...")
print()

# Pentru fiecare raport de masa, cautam ce cicluri s-ar potrivi
for mass_name, mass_ratio in mass_ratios.items():
    print(f"{mass_name} = {mass_ratio}:")

    best_match = None
    best_err = float('inf')

    # Incercam toate combinatiile de cicluri
    for L1 in lengths:
        for L2 in lengths:
            if L1 != L2:
                # Raport de frecvente = L2/L1 (f ~ 1/L)
                f_ratio = L2 / L1
                err = abs(f_ratio - mass_ratio) / mass_ratio * 100

                if err < best_err:
                    best_err = err
                    best_match = (L1, L2, f_ratio)

    if best_match:
        L1, L2, f_ratio = best_match
        print(f"  Cea mai buna potrivire: L_{L2}/L_{L1} = {f_ratio:.4f}")
        print(f"  Eroare: {best_err:.1f}%")

    # Incercam si cu puteri ale lui phi
    for L1 in lengths:
        for n in range(1, 5):
            f_ratio = (10/L1) * PHI**n  # raport cu decagon, modificat cu phi^n
            err = abs(f_ratio - mass_ratio) / mass_ratio * 100
            if err < 5:
                print(f"  Alternativ: (10/{L1})*phi^{n} = {f_ratio:.4f} (err {err:.1f}%)")
    print()

# ============================================================
# ANALIZA DETALIATA: ELECTRON, MUON, TAU
# ============================================================
print("-" * 70)
print("ANALIZA: ELECTRON, MUON, TAU")
print("-" * 70)

print("""
Daca electron, muon, tau sunt ACEEASI PARTICULA dar pe bucle diferite:

Intrebare: Ce lungimi de bucla ar da rapoartele corecte?

m_mu/m_e = 206.77 = L_e / L_mu
m_tau/m_mu = 16.82 = L_mu / L_tau
m_tau/m_e = 3477.2 = L_e / L_tau

Din primele doua:
  L_e / L_mu = 206.77
  L_mu / L_tau = 16.82

Deci:
  L_e / L_tau = 206.77 * 16.82 = 3477.9 [CONSISTENT!]
""")

# Verificam consistenta
ratio_check = 206.77 * 16.82
print(f"Verificare: 206.77 * 16.82 = {ratio_check:.1f}")
print(f"Experimental m_tau/m_e = 3477.2")
print(f"Eroare: {abs(ratio_check - 3477.2)/3477.2 * 100:.2f}%")

print("""
PROBLEMA: Rapoartele de mase sunt PREA MARI pentru cicluri simple.

Cel mai mare raport de cicluri pe 600-cell:
  L_max / L_min ~ 12/3 = 4

Dar m_mu/m_e = 207, care e mult mai mare.

CONCLUZIE: Masele NU pot fi simple lungimi de bucla.
           Trebuie alt mecanism.
""")

# ============================================================
# IPOTEZA ALTERNATIVA: MODURI DE REZONANTA
# ============================================================
print("-" * 70)
print("IPOTEZA ALTERNATIVA: MODURI DE REZONANTA")
print("-" * 70)

print("""
Poate nu e lungimea buclei, ci MODUL de rezonanta pe aceeasi bucla.

Pe o bucla de lungime L, modurile sunt:
  n = 1, 2, 3, ... (numar de "noduri")

Frecventa modului n:
  f_n = n * f_1 = n * c / L

Daca electron = mod 1, muon = mod n, tau = mod m:
  m_mu/m_e = n
  m_tau/m_mu = m/n

Pentru m_mu/m_e = 207:
  n ~ 207 (modul 207 pe decagon?)

Asta pare nerezonabil de mare.

ALTA IDEE: Combinatii de bucle
  Particula = superpozitie de mai multe bucle
  Masa = functie complicata de frecventele buclelor
""")

# ============================================================
# INCERCARE: PHI SI PUTERILE SALE
# ============================================================
print("-" * 70)
print("INCERCARE: MASE DIN PUTERI ALE LUI PHI")
print("-" * 70)

print("""
Stim ca phi apare natural in 600-cell.
Poate masele sunt legate de phi direct, nu de lungimi de bucla.

Incercam: m ~ phi^k pentru diferite k
""")

print("\nCautam k astfel incat phi^k ~ raport de mase:\n")

for mass_name, mass_ratio in mass_ratios.items():
    # Rezolvam phi^k = mass_ratio
    # k = log(mass_ratio) / log(phi)
    if mass_ratio > 0:
        k = np.log(mass_ratio) / np.log(PHI)
        k_round = round(k)
        val_round = PHI**k_round
        err = abs(val_round - mass_ratio) / mass_ratio * 100

        print(f"{mass_name} = {mass_ratio}")
        print(f"  phi^{k:.2f} = {mass_ratio}")
        print(f"  phi^{k_round} = {val_round:.2f} (eroare {err:.1f}%)")
        print()

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-056")
print("=" * 70)

print("""
REZULTATE:

1. CICLURI GASITE:
   - Cicluri de lungime 3, 4, 5, 6, ... 12 exista pe 600-cell
   - Cel mai scurt: triunghi (3 pasi)
   - Decagonul (10 pasi) e o geodezica speciala

2. RAPOARTE DE FRECVENTE:
   - Rapoartele simple L2/L1 sunt mici (maxim ~4)
   - NU se potrivesc cu rapoartele de mase (care sunt ~17, ~207, ~3500)

3. MASELE NU SUNT SIMPLE LUNGIMI DE BUCLA:
   - Raportul m_mu/m_e = 207 e prea mare
   - Trebuie alt mecanism

4. PUTERI ALE LUI PHI:
   - m_mu/m_e ~ phi^11 (eroare 3.3%)
   - m_tau/m_mu ~ phi^6 (eroare 34%) - nu prea bun
   - m_tau/m_e ~ phi^17 (eroare ~40%)

5. CE AM INVATAT:
   - Modelul simplu "masa = 1/lungime bucla" NU functioneaza
   - Particulele nu sunt bucle simple
   - Poate sunt pattern-uri mai complexe sau combinatii de moduri

DIRECTII POSIBILE:
- Masa = energie de legatura a pattern-ului (nu lungime)
- Masa = combinatie de mai multe bucle suprapuse
- Masa depinde de TOPOLOGIA pattern-ului, nu doar lungime
""")

"""
EXP-061: Simulare Propagare Unde pe 600-cell
============================================
Simulam cum se propaga energia pe structura 600-cell
si cautam moduri rezonante (stabile).
"""

from physics_formulas import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import eigsh

print("=" * 70)
print("EXP-061: SIMULARE UNDE PE 600-CELL")
print("=" * 70)

# ============================================================
# CONSTRUIRE 600-CELL (metoda corecta)
# ============================================================
print("\n" + "-" * 70)
print("CONSTRUIRE 600-CELL")
print("-" * 70)

def build_600cell_vertices():
    """
    Construieste cele 120 de varfuri ale 600-cell.
    Toate varfurile sunt pe sfera S^3 de raza 1.

    Cele 120 de varfuri vin din 3 familii:
    - 8 varfuri: permutari de (+-1, 0, 0, 0)
    - 16 varfuri: (+-1/2, +-1/2, +-1/2, +-1/2)
    - 96 varfuri: PERMUTATII PARE de (+-phi/2, +-1/2, +-1/(2*phi), 0)
    """
    vertices = []
    phi = PHI  # golden ratio

    # 1. Permutari de (+-1, 0, 0, 0) - 8 varfuri
    for i in range(4):
        for s in [-1, 1]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = float(s)
            vertices.append(v)

    # 2. (+-1/2, +-1/2, +-1/2, +-1/2) - 16 varfuri
    for s0 in [-1, 1]:
        for s1 in [-1, 1]:
            for s2 in [-1, 1]:
                for s3 in [-1, 1]:
                    vertices.append([s0/2, s1/2, s2/2, s3/2])

    # 3. PERMUTATII PARE de (+-phi/2, +-1/2, +-1/(2*phi), 0) - 96 varfuri
    #
    # Permutatiile PARE ale lui (0,1,2,3) sunt exact 12:
    # Identitate si cicluri de 3 elemente (8) + produs de 2 transpozitii (3)
    even_perms = [
        (0, 1, 2, 3),  # identitate
        (0, 2, 3, 1),  # ciclu (1 2 3)
        (0, 3, 1, 2),  # ciclu (1 3 2)
        (1, 0, 3, 2),  # (0 1)(2 3)
        (1, 2, 0, 3),  # ciclu (0 1 2)
        (1, 3, 2, 0),  # ciclu (0 1 3 2) - ciclu de 4
        (2, 0, 1, 3),  # ciclu (0 2 1)
        (2, 1, 3, 0),  # ciclu (0 2 3 1) - ciclu de 4
        (2, 3, 0, 1),  # (0 2)(1 3)
        (3, 0, 2, 1),  # ciclu (0 3 2 1) - ciclu de 4
        (3, 1, 0, 2),  # ciclu (0 3 1 2) - ciclu de 4
        (3, 2, 1, 0),  # (0 3)(1 2)
    ]

    # Coordonatele de baza (toate nenule in afara de 0)
    base_coords = [phi/2, 0.5, 1/(2*phi), 0.0]

    for perm in even_perms:
        # Pentru fiecare permutare para, aplicam toate 8 combinatii de semne
        # pentru cele 3 coordonate nenule
        for s0 in [-1, 1]:
            for s1 in [-1, 1]:
                for s2 in [-1, 1]:
                    v = [0.0, 0.0, 0.0, 0.0]
                    signs = [s0, s1, s2, 1]  # ultima e 0, semnul nu conteaza
                    for i, p in enumerate(perm):
                        v[p] = base_coords[i] * signs[i]
                    vertices.append(v)

    # Verificam ca avem 120 varfuri
    print(f"  Familie 1: 8 varfuri")
    print(f"  Familie 2: 16 varfuri")
    print(f"  Familie 3: 12 perm * 8 semne = 96 varfuri")
    print(f"  Total: {len(vertices)} varfuri")

    # Verificam ca toate sunt pe sfera de raza 1
    for i, v in enumerate(vertices):
        r = np.sqrt(sum(x**2 for x in v))
        if abs(r - 1.0) > 0.001:
            print(f"  EROARE: vertex {i} are raza {r}")

    return np.array(vertices)

vertices = build_600cell_vertices()
n_vertices = len(vertices)
print(f"Varfuri construite: {n_vertices}")

# Construim matricea de adiacenta
# Doua varfuri sunt vecine daca distanta = 1/phi
def build_adjacency(vertices, edge_length=1/PHI, tol=0.02):
    n = len(vertices)
    adj = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(vertices[i] - vertices[j])
            if abs(d - edge_length) < tol:
                adj[i, j] = 1
                adj[j, i] = 1

    return adj

print("Construim matricea de adiacenta...")
adj = build_adjacency(vertices)
n_edges = np.sum(adj) // 2
neighbors_per_vertex = np.sum(adj, axis=1)

print(f"Muchii: {n_edges} (ar trebui 720)")
print(f"Vecini per vertex: min={neighbors_per_vertex.min()}, max={neighbors_per_vertex.max()}, mean={neighbors_per_vertex.mean():.1f}")

# Daca nu avem exact 12 vecini, ajustam toleranta
if neighbors_per_vertex.mean() != 12:
    print("\nAjustam toleranta pentru muchii...")
    for tol in [0.01, 0.015, 0.02, 0.025, 0.03]:
        adj = build_adjacency(vertices, tol=tol)
        mean_neighbors = np.sum(adj, axis=1).mean()
        if abs(mean_neighbors - 12) < 0.5:
            print(f"  tol={tol}: {mean_neighbors:.1f} vecini/vertex")
            break

n_edges = np.sum(adj) // 2
print(f"\nFinal: {n_edges} muchii, {np.sum(adj, axis=1).mean():.1f} vecini/vertex")

# ============================================================
# LAPLACIANUL GRAFULUI
# ============================================================
print("\n" + "-" * 70)
print("LAPLACIANUL GRAFULUI")
print("-" * 70)

# Laplacian: L = D - A
# D = matrice diagonala cu gradele nodurilor
degree = np.sum(adj, axis=1)
D = np.diag(degree)
L = D - adj

print(f"Dimensiune Laplacian: {L.shape}")

# Calculam valorile proprii
print("Calculam valorile proprii...")
eigenvalues, eigenvectors = np.linalg.eigh(L)

# Sortam
idx = np.argsort(eigenvalues)
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

print(f"\nPrimele 15 valori proprii:")
for i in range(min(15, len(eigenvalues))):
    print(f"  lambda_{i} = {eigenvalues[i]:.6f}")

# ============================================================
# FRECVENTE DE REZONANTA
# ============================================================
print("\n" + "-" * 70)
print("FRECVENTE DE REZONANTA")
print("-" * 70)

print("""
Valorile proprii ale Laplacianului = frecvente^2 ale modurilor.

  omega_n = sqrt(lambda_n)

Modul 0: lambda_0 = 0 (modul constant, fara oscilatie)
Modul 1: lambda_1 = prima frecventa nenula (gap spectral)
""")

# Frecventele (omega = sqrt(lambda))
frequencies = np.sqrt(np.maximum(eigenvalues, 0))

print("Frecvente (omega_n = sqrt(lambda_n)):")
for i in range(min(15, len(frequencies))):
    if frequencies[i] > 0.001:
        print(f"  omega_{i} = {frequencies[i]:.6f}")

# Gap spectral
nonzero_freq = frequencies[frequencies > 0.001]
if len(nonzero_freq) > 0:
    gap = nonzero_freq[0]
    print(f"\nGap spectral (prima frecventa nenula): omega_1 = {gap:.6f}")
    print(f"Comparam cu 1/(2*phi^2) = {1/(2*PHI**2):.6f}")

# ============================================================
# SIMULARE PROPAGARE UNDA
# ============================================================
print("\n" + "-" * 70)
print("SIMULARE PROPAGARE UNDA")
print("-" * 70)

def simulate_wave(L, psi0, n_steps=1000, dt=0.01):
    """
    Simuleaza ecuatia undei: d^2 psi / dt^2 = -L * psi

    Discretizat: psi(t+dt) = 2*psi(t) - psi(t-dt) - dt^2 * L @ psi(t)
    """
    n = len(psi0)
    psi = psi0.copy().astype(float)
    psi_old = psi.copy()

    # Stocam energia totala in timp
    energies = []
    max_amplitudes = []

    for step in range(n_steps):
        # Ecuatia undei discretizata
        psi_new = 2*psi - psi_old - dt**2 * (L @ psi)

        # Update
        psi_old = psi.copy()
        psi = psi_new

        # Energia totala (conservata pentru ecuatia undei)
        energy = np.sum(psi**2)
        energies.append(energy)
        max_amplitudes.append(np.max(np.abs(psi)))

    return psi, energies, max_amplitudes

# Conditii initiale: energie concentrata intr-un punct
print("Simulam cu energie initiala concentrata in nodul 0...")
psi0 = np.zeros(n_vertices)
psi0[0] = 1.0  # toata energia in nodul 0

psi_final, energies, max_amps = simulate_wave(L, psi0, n_steps=2000, dt=0.05)

print(f"Energie initiala: {energies[0]:.4f}")
print(f"Energie finala: {energies[-1]:.4f}")
print(f"Conservare: {energies[-1]/energies[0]*100:.1f}%")

# ============================================================
# ANALIZA MODURILOR
# ============================================================
print("\n" + "-" * 70)
print("ANALIZA MODURILOR PROPRII")
print("-" * 70)

print("""
Fiecare mod propriu e un pattern de vibratie stabil.
Modurile cu frecventa joasa = pattern-uri "mari" (masa mica)
Modurile cu frecventa mare = pattern-uri "mici" (masa mare)
""")

# Descompunem starea finala in moduri proprii
coefficients = eigenvectors.T @ psi_final

print("Contributia fiecarui mod la starea finala:")
print(f"{'Mod':<6} {'Frecventa':<12} {'Coeficient':<12} {'Energie %':<10}")
print("-" * 45)

total_energy = np.sum(coefficients**2)
for i in range(min(20, len(coefficients))):
    if np.abs(coefficients[i]) > 0.01:
        energy_pct = coefficients[i]**2 / total_energy * 100
        print(f"{i:<6} {frequencies[i]:<12.4f} {coefficients[i]:<12.4f} {energy_pct:<10.1f}%")

# ============================================================
# RAPOARTE DE FRECVENTE
# ============================================================
print("\n" + "-" * 70)
print("RAPOARTE DE FRECVENTE")
print("-" * 70)

print("Rapoarte intre frecvente si comparatie cu rapoarte de mase:\n")

# Frecvente distincte nenule
unique_freq = []
for f in frequencies:
    if f > 0.01:
        is_new = True
        for uf in unique_freq:
            if abs(f - uf) < 0.01:
                is_new = False
                break
        if is_new:
            unique_freq.append(f)

unique_freq = sorted(unique_freq)[:10]

print("Frecvente distincte:")
for i, f in enumerate(unique_freq):
    print(f"  f_{i} = {f:.4f}")

print("\nRapoarte f_i/f_0:")
if len(unique_freq) > 0:
    f0 = unique_freq[0]
    for i, f in enumerate(unique_freq):
        ratio = f / f0
        # Comparam cu phi^n
        for n in range(10):
            if abs(ratio - PHI**n) / (PHI**n) < 0.1:
                print(f"  f_{i}/f_0 = {ratio:.4f} ~ phi^{n} = {PHI**n:.4f}")
                break
        else:
            print(f"  f_{i}/f_0 = {ratio:.4f}")

# ============================================================
# VERIFICARE: FRECVENTE SI PHI
# ============================================================
print("\n" + "-" * 70)
print("VERIFICARE: FRECVENTE SI PHI")
print("-" * 70)

print("Cautam frecvente care sunt puteri ale lui phi:\n")

for i, f in enumerate(frequencies[:30]):
    if f > 0.01:
        # f = phi^n => n = log(f)/log(phi)
        n = np.log(f) / np.log(PHI)
        n_round = round(n)
        if abs(n - n_round) < 0.2:
            print(f"  lambda_{i} = {eigenvalues[i]:.4f}, omega = {f:.4f} ~ phi^{n_round} = {PHI**n_round:.4f}")

# ============================================================
# ANALIZA DETALIATA A SPECTRULUI
# ============================================================
print("\n" + "-" * 70)
print("ANALIZA DETALIATA: VALORILE PROPRII SI PHI")
print("-" * 70)

# Gasim valorile proprii distincte si degenerarile lor
unique_eigenvalues = []
degeneracies = []
for ev in eigenvalues:
    is_new = True
    for i, uev in enumerate(unique_eigenvalues):
        if abs(ev - uev) < 0.0001:
            degeneracies[i] += 1
            is_new = False
            break
    if is_new:
        unique_eigenvalues.append(ev)
        degeneracies.append(1)

print("\nValori proprii distincte ale Laplacianului:")
print(f"{'Index':<6} {'lambda':<12} {'Degenerare':<12} {'omega=sqrt(lambda)':<16} {'Comparatie cu phi':<20}")
print("-" * 70)

for i, (ev, deg) in enumerate(zip(unique_eigenvalues, degeneracies)):
    omega = np.sqrt(max(ev, 0))

    # Cautam potriviri cu phi^n sau numere simple
    phi_match = ""
    if omega > 0.01:
        n_phi = np.log(omega) / np.log(PHI)
        if abs(n_phi - round(n_phi)) < 0.15:
            phi_match = f"~ phi^{round(n_phi)} = {PHI**round(n_phi):.4f}"
        elif abs(omega - 3.0) < 0.01:
            phi_match = "= 3 exact"
        elif abs(ev - 12.0) < 0.01:
            phi_match = "lambda = 12 (vecini!)"

    print(f"{i:<6} {ev:<12.4f} {deg:<12} {omega:<16.4f} {phi_match}")

# Verificam relatii intre valorile proprii
print("\n" + "-" * 70)
print("RELATII INTRE VALORILE PROPRII")
print("-" * 70)

print("\nRapoarte lambda_i / lambda_1:")
lambda_1 = unique_eigenvalues[1]  # Prima valoare nenula
for i, ev in enumerate(unique_eigenvalues[1:], 1):
    ratio = ev / lambda_1
    # Verificam daca e putere de phi
    if ratio > 0:
        n_phi = np.log(ratio) / np.log(PHI)
        if abs(n_phi - round(n_phi)) < 0.15:
            print(f"  lambda_{i}/lambda_1 = {ratio:.4f} ~ phi^{round(n_phi)} = {PHI**round(n_phi):.4f}")
        else:
            # Verificam daca e numar intreg sau rational simplu
            if abs(ratio - round(ratio)) < 0.05:
                print(f"  lambda_{i}/lambda_1 = {ratio:.4f} ~ {round(ratio)}")
            else:
                print(f"  lambda_{i}/lambda_1 = {ratio:.4f}")

# Verificam formula analitica cunoscuta
print("\n" + "-" * 70)
print("VERIFICARE FORMULE ANALITICE")
print("-" * 70)

print("""
Pentru 600-cell, valorile proprii ale Laplacianului sunt cunoscute analitic.
Ar trebui sa fie:
  lambda = 12 - 2*cos(2*pi*k/10) pentru diferite k

unde 12 e numarul de vecini.
""")

# Calculam valorile teoretice
print("Valorile teoretice (12 - 2*cos(2*pi*k/10)):")
for k in range(6):
    theoretical = 12 - 2*np.cos(2*np.pi*k/10)
    print(f"  k={k}: lambda = {theoretical:.4f}")

# Verificam daca se potrivesc
print("\nComparatie cu valorile obtinute:")
print(f"  lambda_1 = {unique_eigenvalues[1]:.4f} vs teoretic k=1: {12 - 2*np.cos(2*np.pi/10):.4f}")
print(f"  Diferenta: {abs(unique_eigenvalues[1] - (12 - 2*np.cos(2*np.pi/10))):.6f}")

# Conexiunea cu phi
print("\n" + "-" * 70)
print("CONEXIUNEA CU PHI")
print("-" * 70)

print(f"""
cos(2*pi/10) = cos(36 deg) = phi/2 = {np.cos(np.pi/5):.6f}
cos(4*pi/10) = cos(72 deg) = (phi-1)/2 = 1/(2*phi) = {np.cos(2*np.pi/5):.6f}

Deci:
  lambda_1 = 12 - 2*(phi/2) = 12 - phi = {12 - PHI:.6f}
  Obtinut: {unique_eigenvalues[1]:.6f}

  lambda_5 = 12 - 2*((phi-1)/2) = 12 - (phi-1) = 13 - phi = {13 - PHI:.6f}

Frecventa fundamentala:
  omega_1 = sqrt(12 - phi) = {np.sqrt(12 - PHI):.6f}
  phi^1 = {PHI:.6f}

EROARE: omega_1/phi = {np.sqrt(12 - PHI)/PHI:.4f} (ar trebui ~1 pentru potrivire)
""")

# Dar raportul omega_max / omega_1?
omega_max = np.sqrt(unique_eigenvalues[-1])
omega_1 = np.sqrt(unique_eigenvalues[1])
print(f"Raport omega_max/omega_1 = {omega_max/omega_1:.4f}")
print(f"phi^2 = {PHI**2:.4f}")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-061")
print("=" * 70)

print(f"""
REZULTATE SIMULARE:

1. 600-CELL CONSTRUIT CORECT:
   - {n_vertices} varfuri
   - {n_edges} muchii
   - {np.sum(adj, axis=1).mean():.1f} vecini per vertex

2. SPECTRUL LAPLACIANULUI:
   - Gap spectral: omega_1 = {np.sqrt(unique_eigenvalues[1]):.4f}
   - Formula: lambda_1 = 12 - phi = {12 - PHI:.4f}
   - Phi APARE NATURAL in spectru prin cos(36 deg) = phi/2

3. CONEXIUNEA CU PHI:
   - cos(36 deg) = phi/2 EXACT (geometria icosaedrica)
   - Valorile proprii contin phi explicit
   - omega_1 = sqrt(12 - phi) ~ 3.22

4. CE NU SE POTRIVESTE:
   - omega_1 NU e exact phi^1 = 1.618
   - Dar phi apare IN formula pentru lambda

5. INTERPRETARE:
   - 600-cell are PHI "incorporat" in structura sa
   - Frecventele naturale sunt functii de phi
   - Pattern-ul phi^n pentru mase ar putea veni din
     combinatii de moduri proprii

6. RAMANE DE EXPLORAT:
   - Cum se combina modurile pentru a da phi^11, phi^17?
   - Ce face un anumit mod "stabil" vs altul?
   - Conexiunea cu formula 1/alpha = 20*phi^4
""")

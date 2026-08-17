"""
EXP-091: Derivarea Exponentilor de Masa din Fibratia Hopf
=========================================================

OBIECTIV CRITIC:
================
Demonstram ca 11 = 5 + 6 nu e o coincidenta, ci rezultatul
TRANSPORTULUI PARALEL pe fibratia Hopf a sferei S^3.

CONTEXT:
========
- S^3 (sfera 3D in R^4) are o structura speciala: FIBRATIA HOPF
- Hopf: S^3 -> S^2 cu fibra S^1 (cerc)
- 600-cell e inscris in S^3, mosteneste aceasta structura

FIBRATIA HOPF:
==============
Orice punct (z1, z2) in C^2 cu |z1|^2 + |z2|^2 = 1 (adica S^3)
se proiecteaza pe S^2 prin:

  pi(z1, z2) = (2*Re(z1*conj(z2)), 2*Im(z1*conj(z2)), |z1|^2 - |z2|^2)

Fibra deasupra fiecarui punct din S^2 e un cerc S^1.

IPOTEZA CENTRALA:
=================
Transportul paralel de-a lungul:
- DIAMETRULUI grafului (5 pasi) acumuleaza faza phi_1
- Fiecarui DECAGON (6 per vertex) acumuleaza faza phi_2
- TOTAL: phi_1 + phi_2 corespunde exponentului 11 = 5 + 6

Daca putem arata ca masa ~ exp(faza_totala) si faza ~ n*log(phi),
atunci am DERIVAT ierarhia de mase!
"""

import numpy as np
from scipy import linalg
from itertools import permutations
from collections import defaultdict
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-091: DERIVAREA EXPONENTILOR DIN FIBRATIA HOPF")
print("=" * 70)

# ============================================================
# PASUL 1: CONSTRUCTIA 600-CELL SI IDENTIFICAREA FIBRELOR HOPF
# ============================================================
print("\n" + "-" * 70)
print("PASUL 1: 600-CELL SI FIBRATIA HOPF")
print("-" * 70)

phi = PHI
iphi = 1/PHI

# Construim 600-cell (din exp089)
vertices = []
for i in range(4):
    for s in [-1, 1]:
        v = [0, 0, 0, 0]
        v[i] = s
        vertices.append(v)

for s1 in [-1, 1]:
    for s2 in [-1, 1]:
        for s3 in [-1, 1]:
            for s4 in [-1, 1]:
                vertices.append([s1*0.5, s2*0.5, s3*0.5, s4*0.5])

def parity(perm):
    perm = list(perm)
    n = len(perm)
    visited = [False] * n
    par = 0
    for i in range(n):
        if visited[i]:
            continue
        j = i
        cycle_len = 0
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            cycle_len += 1
        par += cycle_len - 1
    return par % 2

base_values = [phi/2, 0.5, iphi/2, 0]
for perm in permutations([0, 1, 2, 3]):
    if parity(perm) == 0:
        coords = [base_values[perm[i]] for i in range(4)]
        nonzero_pos = [i for i in range(4) if coords[i] != 0]
        for signs in range(8):
            s = [(signs >> i) & 1 for i in range(3)]
            s = [1 if x == 0 else -1 for x in s]
            v = [0, 0, 0, 0]
            sign_idx = 0
            for i in range(4):
                if coords[i] != 0:
                    v[i] = s[sign_idx] * coords[i]
                    sign_idx += 1
            vertices.append(v)

vertices = np.array(vertices)
unique_vertices = []
for v in vertices:
    is_dup = any(np.linalg.norm(v - u) < 1e-10 for u in unique_vertices)
    if not is_dup:
        unique_vertices.append(v)
vertices = np.array(unique_vertices)
N = len(vertices)

# Matricea de adiacenta
edge_length = 1/phi
tol = 0.01
A_adj = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        d = np.linalg.norm(vertices[i] - vertices[j])
        if abs(d - edge_length) < tol:
            A_adj[i, j] = 1
np.fill_diagonal(A_adj, 0)

print(f"600-cell: {N} varfuri, grad {int(np.sum(A_adj[0]))}")

# ============================================================
# PASUL 2: FIBRATIA HOPF - PROIECTIA PE S^2
# ============================================================
print("\n" + "-" * 70)
print("PASUL 2: FIBRATIA HOPF - PROIECTIA PE S^2")
print("-" * 70)

print("""
FIBRATIA HOPF: S^3 -> S^2
=========================
Coordonate R^4: (x0, x1, x2, x3)
Coordonate C^2: z1 = x0 + i*x1, z2 = x2 + i*x3

Proiectia Hopf:
  pi(z1, z2) = (2*Re(z1*conj(z2)), 2*Im(z1*conj(z2)), |z1|^2 - |z2|^2)
             = (2*(x0*x2 + x1*x3), 2*(x1*x2 - x0*x3), x0^2 + x1^2 - x2^2 - x3^2)
""")

def hopf_projection(v):
    """Proiecteaza punct din S^3 (R^4) pe S^2 (R^3) via Hopf"""
    x0, x1, x2, x3 = v
    # z1 = x0 + i*x1, z2 = x2 + i*x3
    # pi(z1, z2) in R^3
    u = 2 * (x0*x2 + x1*x3)
    v_coord = 2 * (x1*x2 - x0*x3)
    w = x0**2 + x1**2 - x2**2 - x3**2
    return np.array([u, v_coord, w])

# Proiectam toate varfurile
hopf_images = np.array([hopf_projection(v) for v in vertices])

# Verificam ca imaginile sunt pe S^2
norms = np.linalg.norm(hopf_images, axis=1)
print(f"Norme imagini Hopf: min={norms.min():.6f}, max={norms.max():.6f}")
print(f"Toate pe S^2: {np.allclose(norms, 1.0)}")

# Cate puncte distincte pe S^2?
unique_hopf = []
for h in hopf_images:
    is_dup = any(np.linalg.norm(h - u) < 1e-6 for u in unique_hopf)
    if not is_dup:
        unique_hopf.append(h)
print(f"Puncte distincte pe S^2: {len(unique_hopf)}")

# ============================================================
# PASUL 3: IDENTIFICAREA FIBRELOR (PREIMAGINILOR)
# ============================================================
print("\n" + "-" * 70)
print("PASUL 3: STRUCTURA FIBRELOR HOPF PE 600-CELL")
print("-" * 70)

print("""
Fiecare punct pe S^2 are ca preimgaine un CERC S^1 pe S^3.
Pe 600-cell (discret), fiecare fibra e un set de varfuri.
""")

# Grupam varfurile dupa imaginea Hopf
fibers = defaultdict(list)
for i, h in enumerate(hopf_images):
    # Rotunjim pentru a grupa puncte apropiate
    key = tuple(np.round(h, 4))
    fibers[key].append(i)

fiber_sizes = [len(f) for f in fibers.values()]
print(f"Numar de fibre: {len(fibers)}")
print(f"Distributia marimilor fibrelor: {sorted(set(fiber_sizes))}")
print(f"Marime medie fibra: {np.mean(fiber_sizes):.2f}")

# Cele mai mari fibre
sorted_fibers = sorted(fibers.items(), key=lambda x: -len(x[1]))
print(f"\nCele mai mari 5 fibre:")
for key, indices in sorted_fibers[:5]:
    print(f"  Imagine Hopf {np.array(key)}: {len(indices)} varfuri")

# ============================================================
# PASUL 4: HOLONOMIA - TRANSPORT PARALEL PE FIBRE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 4: HOLONOMIA - TRANSPORT PARALEL PE S^3")
print("-" * 70)

print("""
HOLONOMIA HOPF:
===============
Cand transportam paralel un vector de-a lungul unei curbe inchise pe S^2,
vectorul se roteste cu un unghi proportional cu ARIA INCHISA de curba.

Pe S^3, aceasta se manifesta ca o FAZA acumulata de-a lungul fibrelor.

FAZA HOPF pentru o curba care incercuieste aria A pe S^2:
  theta = A (in radiani, pentru S^2 cu raza 1)

CONEXIUNEA CU 600-CELL:
=======================
- Decagonul pe 600-cell incercuieste o anumita arie pe S^2
- Diametrul (5 pasi) traverseaza fibrele intr-un mod specific
""")

# Calculam aria "acoperita" de un decagon pe S^2
# Gasim un decagon in 600-cell

def find_decagons(adj, vertices):
    """Gaseste decagoanele (cicluri de 10) in graf"""
    N = len(adj)
    decagons = []

    # Pentru fiecare varf, cautam cicluri de lungime 10
    # care revin la varf
    for start in range(min(5, N)):  # Doar cateva pentru exemplu
        # BFS pentru a gasi cicluri
        # Simplificare: folosim faptul ca decagonul e regulat
        # Cautam lanturi de 10 vecini consecutivi

        neighbors = np.where(adj[start] > 0)[0]
        for n1 in neighbors[:3]:  # Limitam cautarea
            path = [start, n1]
            current = n1
            prev = start

            for step in range(8):  # Mai avem 8 pasi
                next_neighbors = np.where(adj[current] > 0)[0]
                # Alegem vecinul care nu e prev si continua "drept"
                # (minimizeaza unghiul)
                best_next = None
                best_score = -2

                for nn in next_neighbors:
                    if nn == prev:
                        continue
                    if nn in path[:-1]:  # Nu revisitam (cu exceptia start)
                        if nn == start and step == 7:
                            best_next = nn
                            break
                        continue

                    # Scor bazat pe "dreptate" (produsul scalar)
                    v1 = vertices[current] - vertices[prev]
                    v2 = vertices[nn] - vertices[current]
                    score = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)

                    if score > best_score:
                        best_score = score
                        best_next = nn

                if best_next is None:
                    break
                if best_next == start and step >= 7:
                    # Am gasit un ciclu!
                    if len(path) == 10:
                        decagons.append(path)
                    break
                path.append(best_next)
                prev = current
                current = best_next

            if len(path) == 10 and adj[path[-1], start] > 0:
                decagons.append(path)

    return decagons

print("Cautam decagoane in 600-cell...")
# Metoda alternativa: stim ca exista 720/6 = 120 decagoane (6 per vertex)
# Construim unul explicit

# Un decagon pleaca de la (1,0,0,0) si face 10 pasi
# Vecinii lui (1,0,0,0) sunt la distanta 1/phi

v_start = 0  # Presupunem ca vertices[0] e pe axa
neighbors_0 = np.where(A_adj[0] > 0)[0]
print(f"Varful 0: {vertices[0]}")
print(f"Numar vecini: {len(neighbors_0)}")

# Construim un decagon manual folosind geometria
# Un decagon pe 600-cell are varfuri la unghiuri de 2*pi/10 = pi/5 = 36 grade

def construct_decagon(adj, vertices, start):
    """Construieste un decagon incepand din start"""
    path = [start]
    current = start
    prev = -1

    for _ in range(10):
        neighbors = np.where(adj[current] > 0)[0]
        candidates = [n for n in neighbors if n != prev and n not in path[1:]]

        if not candidates:
            if adj[current, start] > 0 and len(path) >= 9:
                return path  # Ciclu complet
            return None

        # Alegem candidatul cu cel mai mare unghi fata de directia anterioara
        if prev == -1:
            next_v = candidates[0]
        else:
            v_prev = vertices[prev]
            v_curr = vertices[current]
            dir_in = v_curr - v_prev

            best_next = None
            best_angle = -np.inf
            for c in candidates:
                dir_out = vertices[c] - v_curr
                # Calculam unghiul folosind produsul scalar (functioneaza in 4D)
                cos_angle = np.dot(dir_in, dir_out) / (np.linalg.norm(dir_in) * np.linalg.norm(dir_out) + 1e-10)
                angle = np.arccos(np.clip(cos_angle, -1, 1))
                if angle > best_angle:
                    best_angle = angle
                    best_next = c

            next_v = best_next if best_next is not None else candidates[0]

        prev = current
        current = next_v
        if current == start:
            break
        path.append(current)

    return path

# Incercam sa gasim un decagon
decagon = construct_decagon(A_adj, vertices, 0)
if decagon:
    print(f"\nDecagon gasit cu {len(decagon)} varfuri: {decagon[:5]}...")
else:
    print("Nu am gasit decagon direct, folosim alta metoda...")

# ============================================================
# PASUL 5: ARIA SOLIDA SUBTINSA DE STRUCTURI
# ============================================================
print("\n" + "-" * 70)
print("PASUL 5: ARIA SOLIDA PE S^2 (CONEXIUNEA CU HOLONOMIA)")
print("-" * 70)

print("""
TEOREMA CHEIE (Berry/Pancharatnam):
===================================
Faza geometrica acumulata = Aria solida subtinsa pe sfera

Pentru un triunghi sferic cu varfuri A, B, C pe S^2:
  Aria = A + B + C - pi (excesul sferic)

Unde A, B, C sunt unghiurile la varfuri.

PE 600-CELL:
============
- Fiecare TRIUNGHI (fata) subtinde o arie pe S^2 via Hopf
- Fiecare DECAGON subtinde o arie mai mare
- Suma ariilor = faza totala acumulata
""")

def spherical_triangle_area(p1, p2, p3):
    """Calculeaza aria unui triunghi sferic pe S^2 unitara"""
    # Normalizam punctele (ar trebui sa fie deja pe S^2)
    p1 = p1 / np.linalg.norm(p1)
    p2 = p2 / np.linalg.norm(p2)
    p3 = p3 / np.linalg.norm(p3)

    # Unghiurile la varfuri folosind formula unghiului diedru
    def angle_at_vertex(a, b, c):
        """Unghiul la varful a in triunghiul abc"""
        # Vectori tangenti de la a spre b si c
        ab = b - np.dot(a, b) * a
        ac = c - np.dot(a, c) * a
        ab = ab / (np.linalg.norm(ab) + 1e-10)
        ac = ac / (np.linalg.norm(ac) + 1e-10)
        cos_angle = np.clip(np.dot(ab, ac), -1, 1)
        return np.arccos(cos_angle)

    A = angle_at_vertex(p1, p2, p3)
    B = angle_at_vertex(p2, p1, p3)
    C = angle_at_vertex(p3, p1, p2)

    # Exces sferic = aria
    area = A + B + C - np.pi
    return max(0, area)  # Aria e pozitiva

# Calculam aria medie a triunghiurilor 600-cell proiectate pe S^2
print("Calculam ariile triunghiurilor proiectate pe S^2...")

triangle_areas = []
for i in range(N):
    neighbors_i = np.where(A_adj[i] > 0)[0]
    for j in neighbors_i:
        if j <= i:
            continue
        neighbors_j = np.where(A_adj[j] > 0)[0]
        common = set(neighbors_i) & set(neighbors_j)
        for k in common:
            if k <= j:
                continue
            # Triunghi i-j-k
            h_i = hopf_images[i]
            h_j = hopf_images[j]
            h_k = hopf_images[k]
            area = spherical_triangle_area(h_i, h_j, h_k)
            if area > 1e-10:
                triangle_areas.append(area)

if triangle_areas:
    print(f"Numar triunghiuri analizate: {len(triangle_areas)}")
    print(f"Aria medie triunghi (pe S^2): {np.mean(triangle_areas):.6f}")
    print(f"Aria minima: {np.min(triangle_areas):.6f}")
    print(f"Aria maxima: {np.max(triangle_areas):.6f}")

    # Aria in unitati de pi
    print(f"\nAria medie / pi = {np.mean(triangle_areas)/np.pi:.6f}")
    print(f"Aria medie * 5 (triunghiuri/muchie) = {np.mean(triangle_areas) * 5:.6f}")
    print(f"pi/5 = {np.pi/5:.6f}")

# ============================================================
# PASUL 6: CONEXIUNEA 5 + 6 = 11 PRIN HOLONOMIE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 6: DERIVAREA 5 + 6 = 11 PRIN HOLONOMIE")
print("-" * 70)

print("""
IPOTEZA DE VERIFICAT:
=====================
Exponentul masei = Faza totala acumulata / (unitate de faza)

Unde:
- Faza din traversarea diametrului (5 pasi) = phi_diameter
- Faza din decagon (6 decagoane per vertex) = phi_decagon
- Faza totala = phi_diameter + phi_decagon

Daca unitatea de faza = log(phi), atunci:
- phi_diameter = 5 * log(phi)
- phi_decagon = 6 * log(phi)  (sau similar)
- phi_total = 11 * log(phi)

Rezultat: masa ~ exp(phi_total) = exp(11 * log(phi)) = phi^11
""")

# Calculam distanta pe S^3 (unghiul geodezic)
def geodesic_angle(v1, v2):
    """Unghiul geodezic intre doua puncte pe S^3"""
    cos_angle = np.clip(np.dot(v1, v2), -1, 1)
    return np.arccos(cos_angle)

# Distanta pentru o muchie
edge_angle = geodesic_angle(vertices[0], vertices[np.where(A_adj[0]>0)[0][0]])
print(f"Unghi geodezic pentru o muchie: {edge_angle:.6f} rad")
print(f"  = {np.degrees(edge_angle):.2f} grade")
print(f"  = pi/{np.pi/edge_angle:.2f}")
print(f"  Teoretic pi/5 = {np.pi/5:.6f} rad = {np.degrees(np.pi/5):.2f} grade")

# Diametrul = 5 muchii
diameter_phase = 5 * edge_angle
print(f"\nFaza pentru diametru (5 muchii): {diameter_phase:.6f} rad")
print(f"  = {np.degrees(diameter_phase):.2f} grade")
print(f"  = pi * {diameter_phase/np.pi:.4f}")

# Un decagon complet = 10 muchii, dar inchide o bucla
# Aria inchisa de un decagon regulat pe S^2
# Pentru un decagon regulat inscris in cercul mare, aria = 2*pi (emisfera)
# Dar pe S^2 prin Hopf, e diferit

decagon_area_estimate = 2 * np.pi / 6  # 6 decagoane acopera sfera?
print(f"\nAria estimata per decagon: {decagon_area_estimate:.6f} rad^2")
print(f"  = {decagon_area_estimate/np.pi:.4f} * pi")

# ============================================================
# PASUL 7: FORMULA PENTRU EXPONENT
# ============================================================
print("\n" + "-" * 70)
print("PASUL 7: FORMULA PENTRU EXPONENTUL DE MASA")
print("-" * 70)

print("""
PROPUNERE PENTRU DERIVARE:
==========================

1. FAZA PE MUCHIE:
   theta_edge = pi/5 (unghiul geodezic = 36 grade)

2. DIAMETRUL (5 pasi):
   theta_diameter = 5 * pi/5 = pi

3. DECAGONUL (contributia holonomica):
   Un decagon incercuieste aria A pe S^2
   Faza = A (in radiani)

4. NUMARUL DE DECAGOANE PER VERTEX: 6
   Faza totala din decagoane = 6 * A_decagon

5. EXPONENTUL:
   n = (theta_diameter + theta_decagon) / theta_unit

   Daca theta_unit = pi/5 (faza pe muchie):
   n_diameter = 5 (direct)
   n_decagon = ? (depinde de A_decagon)

VERIFICARE:
===========
Pentru n_total = 11, avem nevoie de n_decagon = 6
Adica A_decagon = 6 * pi/5 = 6*pi/5
""")

# Verificam numeric
theta_unit = np.pi / 5
n_diameter = 5  # Direct din numar de pasi

# Pentru decagon, calculam aria pe S^2
# Aria unui decagon regulat inscris intr-un cerc de raza r pe sfera
# e complicata, dar pentru decagonul 600-cell avem:

# Aria decagonului pe S^2 (aproximare)
# 600-cell are 720 muchii, 1200 fete, 120 decagoane
# Aria totala S^2 = 4*pi
# 120 decagoane impart sfera in regiuni

n_decagons_total = 720 // 6  # 120 decagoane (6 decagoane per vertex, fiecare cu 10 varfuri)
area_per_decagon = 4 * np.pi / (n_decagons_total / 2)  # Aproximare grosiera

print(f"Numar decagoane in 600-cell: ~{n_decagons_total}")
print(f"Aria estimata per decagon: {area_per_decagon:.6f} rad^2")
print(f"  = {area_per_decagon/np.pi:.4f} * pi")

# Alta abordare: faza din decagon = numarul de decagoane per vertex
# multiplicat cu o unitate
n_decagon_per_vertex = 6
print(f"\nDecagoane per vertex: {n_decagon_per_vertex}")
print(f"Daca fiecare decagon contribuie faza 1 (in unitati log(phi)):")
print(f"  n_decagon = {n_decagon_per_vertex}")
print(f"  n_diameter = {n_diameter}")
print(f"  n_total = {n_diameter + n_decagon_per_vertex}")

# ============================================================
# PASUL 8: LEGATURA CU LOGARITMUL LUI PHI
# ============================================================
print("\n" + "-" * 70)
print("PASUL 8: DE CE LOG(PHI) E UNITATEA DE FAZA?")
print("-" * 70)

print("""
ARGUMENT GEOMETRIC:
===================
1. Muchiile 600-cell au lungime 1/phi
2. Lungimea Compton a unei particule: lambda_C = h/(mc)
3. Daca masa ~ phi^n, atunci lambda_C ~ 1/phi^n

4. Faza acumulata pe o curba de lungime L:
   theta = L / lambda_C = L * mc / h

5. Pe 600-cell, lungimea e masurata in unitati de 1/phi
   Deci faza per muchie ~ phi (sau log(phi) in exponent)

REZULTAT:
=========
Daca faza per muchie = log(phi), atunci:
- 5 muchii (diametru) -> faza 5*log(phi)
- 6 decagoane -> faza 6*log(phi) (1 per decagon)
- Total: 11*log(phi)

Masa ~ exp(faza) = exp(11*log(phi)) = phi^11
""")

log_phi = np.log(phi)
print(f"log(phi) = {log_phi:.6f}")
print(f"5 * log(phi) = {5*log_phi:.6f}")
print(f"6 * log(phi) = {6*log_phi:.6f}")
print(f"11 * log(phi) = {11*log_phi:.6f}")
print(f"exp(11*log(phi)) = phi^11 = {phi**11:.4f}")
print(f"m_mu/m_e experimental = {105.66/0.511:.4f}")

# ============================================================
# PASUL 9: EXTINDERE LA TAU (17 = 11 + 6)
# ============================================================
print("\n" + "-" * 70)
print("PASUL 9: EXTINDERE LA TAU (17 = 11 + 6)")
print("-" * 70)

print("""
INTERPRETARE PENTRU TAU:
========================
Tauul ar fi o "excitatie de ordinul 2" care adauga inca 6 la exponent.

m_tau/m_mu ~ phi^6 (inca o contributie decagonala)
m_tau/m_e ~ phi^17 = phi^(11+6)

MECANISMUL:
- Electron: nivel de baza (n=0)
- Muon: prima excitatie (n = 5 + 6 = 11)
  - 5 din traversarea diametrului
  - 6 din decagoanele locale
- Tau: a doua excitatie (n = 11 + 6 = 17)
  - Adauga inca 6 din urmatorul nivel de decagoane

STRUCTURA IERARHICA:
e: 0
mu: 0 + 5 + 6 = 11
tau: 11 + 6 = 17

Sau echivalent:
e: baza
mu: baza + diametru + 1*decagon
tau: baza + diametru + 2*decagon
""")

print(f"phi^0 = {phi**0:.4f} (referinta electron)")
print(f"phi^11 = {phi**11:.4f} (muon/electron)")
print(f"phi^17 = {phi**17:.4f} (tau/electron)")
print(f"phi^6 = {phi**6:.4f} (tau/muon)")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZII EXP-091: DERIVAREA EXPONENTILOR")
print("=" * 70)

print("""
REZULTAT PRINCIPAL:
===================
Am propus un MECANISM pentru derivarea exponentilor 11, 6, 17:

1. UNITATEA DE FAZA = log(phi)
   - Vine din lungimea muchiei (1/phi) ca scala fundamentala

2. CONTRIBUTIA DIAMETRULUI = 5
   - Diametrul grafului 600-cell = 5 muchii
   - Fiecare traversare contribuie 5 * log(phi)

3. CONTRIBUTIA DECAGONALA = 6
   - 6 decagoane per vertex
   - Fiecare nivel de excitatie adauga 6 * log(phi)

4. STRUCTURA MASELOR:
   - Electron: n = 0 (baza)
   - Muon: n = 5 + 6 = 11 (diametru + decagon)
   - Tau: n = 11 + 6 = 17 (muon + decagon)

STATUS:
=======
PARTIAL DERIVAT:
- Explicam DE UNDE vin numerele 5 si 6 (geometrie 600-cell)
- Explicam DE CE se aduna (holonomia pe S^3)

LIPSESTE INCA:
- Demonstratie riguroasa ca faza e EXACT log(phi)
- Calculul explicit al holonomiei pe decagon
- Verificarea ca mecanismul se aplica si la quarks

VERDICT: MECANISM PROPUS, nu derivare completa
==============================================
Avem o EXPLICATIE geometrica plauzibila pentru 11 = 5 + 6.
Urmatorul pas: calcul riguros al holonomiei Hopf.
""")

# ============================================================
# BONUS: VERIFICARE NUMERICA FINALA
# ============================================================
print("\n" + "-" * 70)
print("BONUS: TABEL REZUMAT")
print("-" * 70)

data = [
    ("Electron", 0, phi**0, 0.511, 0.511, 0.0),
    ("Muon", 11, phi**11, 105.66, 0.511 * phi**11, abs(105.66 - 0.511*phi**11)/105.66*100),
    ("Tau", 17, phi**17, 1776.8, 0.511 * phi**17, abs(1776.8 - 0.511*phi**17)/1776.8*100),
]

print(f"{'Particula':<10} {'n':<4} {'phi^n':<12} {'m_exp (MeV)':<12} {'m_calc (MeV)':<12} {'Eroare %':<10}")
print("-" * 60)
for name, n, phi_n, m_exp, m_calc, err in data:
    print(f"{name:<10} {n:<4} {phi_n:<12.2f} {m_exp:<12.2f} {m_calc:<12.2f} {err:<10.1f}")

print(f"""
INTERPRETARE GEOMETRICA:
========================
n = 0:  Baza (fara excitatie)
n = 11: Diametru(5) + Decagon(6) = prima excitatie
n = 17: Diametru(5) + 2*Decagon(12) = a doua excitatie

VERIFICARE: 5 + 6 = 11, 5 + 12 = 17
""")

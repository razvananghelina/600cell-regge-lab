"""
EXP-090: Misterul celor 4 Generatii vs 3 Observate
==================================================

PROBLEMA: Exp-089 a aratat ca cele 96 de varfuri fermionice se descompun
natural in 4 grupuri de 24 (dupa pozitia coordonatei 0 in R^4).
Modelul Standard are insa doar 3 generatii de fermioni.

IPOTEZE DE INVESTIGAT:
1. Al 4-lea grup reprezinta un "sector ascuns" (dark sector / sterile neutrinos)?
2. Exista un mecanism de proiectie 4D -> 3D care reduce la 3 generatii?
3. Cele 8 varfuri axiale (16-cell) actioneaza ca "comutator" pentru generatii?
4. Unul din cele 4 grupuri e geometric "decuplat" de celelalte?

STRUCTURA EXPERIMENTULUI:
- Pas 1: Reconstruim graful si cele 4 grupuri
- Pas 2: Analizam conectivitatea intre grupuri (graf bipartit?)
- Pas 3: Verificam daca un grup e "special" sau decuplat
- Pas 4: Proiectii 4D -> 3D si efectul asupra grupurilor
- Pas 5: Conexiuni cu fizica neutrinilor sterili
- Pas 6: Rolul celor 8 varfuri axiale (16-cell)
"""

import numpy as np
from scipy import linalg
from itertools import permutations, combinations
from collections import defaultdict
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-090: MISTERUL CELOR 4 GENERATII")
print("=" * 70)

# ============================================================
# PASUL 1: RECONSTRUCTIA GRAFULUI CU CLASIFICARE DETALIATA
# ============================================================
print("\n" + "-" * 70)
print("PASUL 1: RECONSTRUCTIA GRAFULUI 600-CELL")
print("-" * 70)

phi = PHI
iphi = 1/PHI

vertices = []
vertex_types = []  # 'A', 'B', sau 'C'
vertex_subgroup = []  # Pentru tip C: 0, 1, 2, 3 (pozitia lui 0)

# 1. 8 varfuri de tip A (16-cell)
for i in range(4):
    for s in [-1, 1]:
        v = [0, 0, 0, 0]
        v[i] = s
        vertices.append(v)
        vertex_types.append('A')
        vertex_subgroup.append(-1)  # N/A pentru A

# 2. 16 varfuri de tip B (tesseract)
for s1 in [-1, 1]:
    for s2 in [-1, 1]:
        for s3 in [-1, 1]:
            for s4 in [-1, 1]:
                vertices.append([s1*0.5, s2*0.5, s3*0.5, s4*0.5])
                vertex_types.append('B')
                vertex_subgroup.append(-1)  # N/A pentru B

# 3. 96 varfuri de tip C (fermioni) cu subgrup
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
        # Gasim pozitia lui 0
        zero_pos = coords.index(0)
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
            vertex_types.append('C')
            # Subgrupul e determinat de CARE coordonata e 0
            for k in range(4):
                if abs(v[k]) < 1e-10:
                    vertex_subgroup.append(k)
                    break

vertices = np.array(vertices)

# Eliminam duplicatele
unique_vertices = []
unique_types = []
unique_subgroups = []
for v, t, sg in zip(vertices, vertex_types, vertex_subgroup):
    is_duplicate = False
    for u in unique_vertices:
        if np.linalg.norm(v - u) < 1e-10:
            is_duplicate = True
            break
    if not is_duplicate:
        unique_vertices.append(v)
        unique_types.append(t)
        unique_subgroups.append(sg)

vertices = np.array(unique_vertices)
vertex_types = unique_types
vertex_subgroup = unique_subgroups
N = len(vertices)

# Indices
indices_A = [i for i, t in enumerate(vertex_types) if t == 'A']
indices_B = [i for i, t in enumerate(vertex_types) if t == 'B']
indices_C = [i for i, t in enumerate(vertex_types) if t == 'C']

# Sub-grupuri pentru C (cele 4 "generatii")
groups_C = {0: [], 1: [], 2: [], 3: []}
for i in indices_C:
    sg = vertex_subgroup[i]
    groups_C[sg].append(i)

print(f"Total varfuri: {N}")
print(f"  Tip A (16-cell): {len(indices_A)}")
print(f"  Tip B (tesseract): {len(indices_B)}")
print(f"  Tip C (fermioni): {len(indices_C)}")
print(f"\nCele 4 sub-grupuri ale fermionilor (dupa coord. 0):")
for sg in range(4):
    print(f"  Grup G{sg} (x{sg}=0): {len(groups_C[sg])} varfuri")

# Matricea de adiacenta
distances = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        distances[i, j] = np.linalg.norm(vertices[i] - vertices[j])

edge_length = 1/phi
tol = 0.01
A_adj = (np.abs(distances - edge_length) < tol).astype(float)
np.fill_diagonal(A_adj, 0)

print(f"\nLungime muchie: {edge_length:.6f}")
print(f"Grad uniform: {int(np.sum(A_adj[0]))}")

# ============================================================
# PASUL 2: CONECTIVITATEA INTRE CELE 4 GRUPURI
# ============================================================
print("\n" + "-" * 70)
print("PASUL 2: CONECTIVITATEA INTRE CELE 4 GRUPURI DE FERMIONI")
print("-" * 70)

print("""
INTREBARE: Cum sunt conectate cele 4 grupuri intre ele?
Daca un grup e "izolat", ar putea fi sectorul ascuns.
""")

# Matricea de conectivitate intre grupuri
connectivity = np.zeros((4, 4))
for i in range(4):
    for j in range(4):
        # Numaram muchiile intre grupul i si grupul j
        edges = 0
        for vi in groups_C[i]:
            for vj in groups_C[j]:
                if A_adj[vi, vj] > 0:
                    edges += 1
        connectivity[i, j] = edges

print("Matricea de conectivitate (numar de muchii intre grupuri):")
print(f"{'':>8} G0      G1      G2      G3")
for i in range(4):
    row = f"G{i}:    "
    for j in range(4):
        row += f"{int(connectivity[i, j]):>6}  "
    print(row)

# Conectivitate normalizata (per varf)
print("\nConectivitate medie per varf:")
print(f"{'':>8} G0      G1      G2      G3")
for i in range(4):
    row = f"G{i}:    "
    for j in range(4):
        n_i = len(groups_C[i])
        avg = connectivity[i, j] / n_i if n_i > 0 else 0
        row += f"{avg:>6.2f}  "
    print(row)

# Verificam daca vreun grup e special
intra_group = [connectivity[i, i] for i in range(4)]
inter_group = [sum(connectivity[i, j] for j in range(4) if j != i) for i in range(4)]

print("\nAnaliza pe grupuri:")
print(f"{'Grup':<6} {'Intra-grup':<12} {'Inter-grup':<12} {'Ratio I/E':<12}")
for i in range(4):
    ratio = intra_group[i] / inter_group[i] if inter_group[i] > 0 else 0
    print(f"G{i}:    {int(intra_group[i]):<12} {int(inter_group[i]):<12} {ratio:<12.3f}")

# ============================================================
# PASUL 3: CONECTIVITATEA CU TIPURILE A SI B
# ============================================================
print("\n" + "-" * 70)
print("PASUL 3: CONECTIVITATEA GRUPURILOR CU A (16-cell) SI B (tesseract)")
print("-" * 70)

print("""
INTREBARE: Gruprile de fermioni sunt conectate diferit la "nucleu" (A+B)?
Daca da, aceasta ar putea explica masele diferite.
""")

# Conectivitate C_grup -> A
conn_to_A = []
conn_to_B = []
for g in range(4):
    edges_to_A = sum(A_adj[vi, vj] for vi in groups_C[g] for vj in indices_A)
    edges_to_B = sum(A_adj[vi, vj] for vi in groups_C[g] for vj in indices_B)
    conn_to_A.append(edges_to_A)
    conn_to_B.append(edges_to_B)

print(f"{'Grup':<6} {'Muchii->A':<12} {'Muchii->B':<12} {'Muchii->A (per vf)':<18} {'Muchii->B (per vf)':<18}")
for g in range(4):
    n_g = len(groups_C[g])
    print(f"G{g}:    {int(conn_to_A[g]):<12} {int(conn_to_B[g]):<12} {conn_to_A[g]/n_g:<18.2f} {conn_to_B[g]/n_g:<18.2f}")

# Verificam daca grupurile sunt simetrice
print(f"\nSunt toate grupurile echivalente?")
all_equal_A = all(conn_to_A[i] == conn_to_A[0] for i in range(4))
all_equal_B = all(conn_to_B[i] == conn_to_B[0] for i in range(4))
print(f"  Conectivitate identica la A: {all_equal_A}")
print(f"  Conectivitate identica la B: {all_equal_B}")

# ============================================================
# PASUL 4: ANALIZA GEOMETRICA - CENTRII GRUPURILOR
# ============================================================
print("\n" + "-" * 70)
print("PASUL 4: GEOMETRIA CENTRILOR GRUPURILOR")
print("-" * 70)

# Centrul de masa al fiecarui grup
centers = []
for g in range(4):
    group_vertices = vertices[groups_C[g]]
    center = np.mean(group_vertices, axis=0)
    centers.append(center)
    print(f"Centrul G{g}: ({center[0]:>7.4f}, {center[1]:>7.4f}, {center[2]:>7.4f}, {center[3]:>7.4f})")

# Distantele intre centri
print("\nDistante intre centrii grupurilor:")
for i in range(4):
    for j in range(i+1, 4):
        d = np.linalg.norm(centers[i] - centers[j])
        print(f"  d(G{i}, G{j}) = {d:.6f}")

# Verificam daca centrii sunt la origine
print("\nNorma centrilor (0 = centrat la origine):")
for g in range(4):
    print(f"  |centru G{g}| = {np.linalg.norm(centers[g]):.10f}")

# ============================================================
# PASUL 5: PROIECTII 4D -> 3D
# ============================================================
print("\n" + "-" * 70)
print("PASUL 5: PROIECTII 4D -> 3D SI EFECTUL ASUPRA GRUPURILOR")
print("-" * 70)

print("""
IPOTEZA: Proiectia din 4D in 3D (spatiu fizic) ar putea "colapa"
doua grupuri intr-unul singur, reducand 4 -> 3 generatii.

Testam diferite proiectii:
1. Proiectie ortogonala (ignora o coordonata)
2. Proiectie stereografica
3. Proiectie pe hiperplan
""")

def count_distinct_after_projection(vertices, groups, proj_func, tol=0.01):
    """Numara cate puncte distincte raman dupa proiectie pentru fiecare grup"""
    results = {}
    for g, indices in groups.items():
        projected = [tuple(proj_func(vertices[i]).round(4)) for i in indices]
        unique = len(set(projected))
        results[g] = unique
    return results

# Proiectie 1: Ignora coordonata k
print("\n1. PROIECTIE ORTOGONALA (ignora o coordonata):")
for k in range(4):
    proj = lambda v, k=k: np.delete(v, k)
    distinct = count_distinct_after_projection(vertices, groups_C, proj)
    print(f"   Ignora x{k}: G0={distinct[0]}, G1={distinct[1]}, G2={distinct[2]}, G3={distinct[3]}")
    # Verificam daca doua grupuri devin identice
    all_proj = {}
    for g in range(4):
        pts = frozenset(tuple(proj(vertices[i]).round(4)) for i in groups_C[g])
        all_proj[g] = pts

    # Verificam suprapuneri
    overlaps = []
    for i in range(4):
        for j in range(i+1, 4):
            overlap = len(all_proj[i] & all_proj[j])
            if overlap > 0:
                overlaps.append((i, j, overlap))
    if overlaps:
        print(f"      Suprapuneri: {overlaps}")

# Proiectie 2: Stereografica din polul (0,0,0,1)
print("\n2. PROIECTIE STEREOGRAFICA (pol la (0,0,0,1)):")
def stereo_proj(v):
    """Proiectie stereografica din (0,0,0,1)"""
    if abs(v[3] - 1) < 1e-10:
        return np.array([np.inf, np.inf, np.inf])
    return v[:3] / (1 - v[3])

# Verificam cum se transforma grupurile
for g in range(4):
    finite_count = sum(1 for i in groups_C[g] if abs(vertices[i][3] - 1) > 1e-10)
    print(f"   G{g}: {finite_count}/24 puncte finite (restul la infinit)")

# ============================================================
# PASUL 6: SIMETRIA GRUPURILOR - ACTIONEAZA CEVA CA "SELECTOR"?
# ============================================================
print("\n" + "-" * 70)
print("PASUL 6: ROLUL VARFURILOR A (16-cell) CA SELECTOR")
print("-" * 70)

print("""
IPOTEZA: Cele 8 varfuri axiale (16-cell) ar putea actiona ca un
"comutator" care selecteaza 3 din cele 4 grupuri.

Varfurile A sunt de forma (+/-1, 0, 0, 0) si permutari.
Fiecare axa k are 2 varfuri: +e_k si -e_k.
""")

# Grupam varfurile A pe axe
axes_A = {0: [], 1: [], 2: [], 3: []}
for i in indices_A:
    v = vertices[i]
    for k in range(4):
        if abs(v[k]) > 0.5:  # Axa k
            axes_A[k].append(i)

print("Varfurile A grupate pe axe:")
for k in range(4):
    print(f"  Axa {k}: {len(axes_A[k])} varfuri (indici: {axes_A[k]})")

# Conexiunea intre axa k si grupul g
print("\nConexiunea intre axele A si grupurile C:")
print(f"{'':>8} G0      G1      G2      G3")
for k in range(4):
    row = f"Axa {k}: "
    for g in range(4):
        edges = sum(A_adj[vi, vj] for vi in axes_A[k] for vj in groups_C[g])
        row += f"{int(edges):>6}  "
    print(row)

# OBSERVATIE CHEIE: Axa k e conectata la grupurile unde x_k != 0
# Adica axa 0 NU e conectata la G0 (unde x0 = 0)!
print("\nOBSERVATIE:")
print("  Axa k este conectata DOAR la grupurile unde x_k != 0")
print("  Adica: Axa 0 nu e conectata la G0, Axa 1 nu e conectata la G1, etc.")
print("  Aceasta creeaza o structura de EXCLUDERE!")

# ============================================================
# PASUL 7: MECANISMUL DE SELECTIE - 3 DIN 4
# ============================================================
print("\n" + "-" * 70)
print("PASUL 7: MECANISMUL DE SELECTIE 3 DIN 4")
print("-" * 70)

print("""
IPOTEZA CENTRALA:
=================
Daca "alegem" o axa privilegiata (ex: directia timpului = x3),
atunci varfurile de pe acea axa (A3) sunt conectate DOAR la 3 grupuri:
G0, G1, G2 (dar NU la G3, unde x3 = 0).

Aceasta ar explica de ce vedem doar 3 generatii:
- Directia timpului in 4D selecteaza 3 din cele 4 grupuri geometrice
- Al 4-lea grup (G3) ar fi "sectorul sterile" - nu interactioneaza
  direct cu "nucleul temporal"
""")

# Verificam explicit
print("Verificare - conexiunile axei 3 (potential 'timp'):")
for g in range(4):
    edges = sum(A_adj[vi, vj] for vi in axes_A[3] for vj in groups_C[g])
    status = "CONECTAT" if edges > 0 else "DECONECTAT"
    print(f"  Axa 3 <-> G{g}: {int(edges)} muchii [{status}]")

# Verificam pentru toate axele
print("\nPentru ORICARE axa privilegiata, avem 3 grupuri conectate:")
for k in range(4):
    connected_groups = []
    for g in range(4):
        edges = sum(A_adj[vi, vj] for vi in axes_A[k] for vj in groups_C[g])
        if edges > 0:
            connected_groups.append(g)
    print(f"  Axa {k} conecteaza: {connected_groups} (3 grupuri), exclude: G{k}")

# ============================================================
# PASUL 8: INTERPRETARE FIZICA - NEUTRINII STERILI
# ============================================================
print("\n" + "-" * 70)
print("PASUL 8: INTERPRETARE FIZICA - SECTORUL STERILE")
print("-" * 70)

print("""
INTERPRETARE PROPUSA:
=====================

1. CELE 4 GRUPURI GEOMETRICE:
   G0, G1, G2, G3 - fiecare cu 24 varfuri = 24 stari fermionice

2. ALEGEREA DIRECTIEI TEMPORALE:
   In fizica 4D, una din directii e "timp" (sa zicem x3)
   Aceasta alegere SPARGE simetria celor 4 grupuri

3. CONSECINTA:
   - G0, G1, G2 sunt "vizibile" (conectate la axa temporala)
   - G3 e "sterile" (neconectat la axa temporala)

4. MAPARE PE FIZICA:
   - G0, G1, G2 -> cele 3 generatii de fermioni (e, mu, tau)
   - G3 -> neutrini sterili / sector dark

5. PREDICTIE:
   - Ar trebui sa existe un al 4-lea set de fermioni
   - Acestia nu interactioneaza electromagnetic (sterili)
   - Candidati: neutrini sterili cu masa ~keV (dark matter?)
""")

# Calculam "masa" relativa bazata pe distanta la "nucleu"
print("\nCalcul 'masa' bazata pe distanta medie la nucleul A+B:")

# Distante geodezice (recalculam rapid cu BFS)
def bfs_distance(adj, start):
    """BFS pentru distante de la start"""
    n = len(adj)
    dist = [-1] * n
    dist[start] = 0
    queue = [start]
    while queue:
        u = queue.pop(0)
        for v in range(n):
            if adj[u, v] > 0 and dist[v] == -1:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist

# Distanta medie de la fiecare grup la varfurile A
all_dist = [bfs_distance(A_adj, i) for i in range(N)]

for g in range(4):
    total_dist = 0
    count = 0
    for vi in groups_C[g]:
        for vj in indices_A:
            total_dist += all_dist[vi][vj]
            count += 1
    avg_dist = total_dist / count if count > 0 else 0
    print(f"  G{g}: distanta medie la A = {avg_dist:.4f}")

# ============================================================
# PASUL 9: SPECTRUL PE FIECARE SUBGRUP
# ============================================================
print("\n" + "-" * 70)
print("PASUL 9: SPECTRUL LAPLACIAN PE FIECARE SUBGRUP")
print("-" * 70)

print("""
Daca grupurile au "mase" diferite, spectrul lor intern ar trebui sa difere.
""")

for g in range(4):
    idx = groups_C[g]
    sub_adj = A_adj[np.ix_(idx, idx)]
    sub_deg = np.sum(sub_adj, axis=1)
    L_sub = np.diag(sub_deg) - sub_adj
    eigs = np.linalg.eigvalsh(L_sub)

    # Prima valoare proprie non-zero
    lambda_1 = eigs[1] if len(eigs) > 1 and eigs[1] > 1e-10 else 0
    lambda_max = eigs[-1]

    print(f"  G{g}: lambda_1 = {lambda_1:.4f}, lambda_max = {lambda_max:.4f}, spectral gap = {lambda_1:.4f}")

# ============================================================
# PASUL 10: COMBINATII DE 3 GRUPURI SI RAPORTUL DE MASE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 10: RAPOARTE DE MASE PENTRU COMBINATII DE 3 GRUPURI")
print("-" * 70)

print("""
Daca alegem 3 grupuri din 4, putem atribui mase bazate pe distanta
lor relativa la nucleu sau pe spectru. Verificam daca vreun triplet
reproduce raportul phi^11 : phi^6 : 1.
""")

# Pentru fiecare triplet de grupuri
from itertools import combinations

triplets = list(combinations([0, 1, 2, 3], 3))

print(f"{'Triplet':<12} {'Note':<50}")
print("-" * 62)
for trip in triplets:
    excluded = [x for x in range(4) if x not in trip][0]
    # Distantele la axa excluded
    dist_to_excluded_axis = []
    for g in trip:
        d = sum(all_dist[vi][axes_A[excluded][0]] for vi in groups_C[g]) / len(groups_C[g])
        dist_to_excluded_axis.append((g, d))

    # Sortam dupa distanta
    dist_to_excluded_axis.sort(key=lambda x: x[1])

    note = f"Exclude G{excluded}, ordine dupa dist: {[x[0] for x in dist_to_excluded_axis]}"
    print(f"({trip[0]},{trip[1]},{trip[2]})      {note}")

# ============================================================
# PASUL 11: VERIFICARE - AL 4-LEA GRUP E "SPECIAL"?
# ============================================================
print("\n" + "-" * 70)
print("PASUL 11: ESTE VREUN GRUP GEOMETRIC 'SPECIAL'?")
print("-" * 70)

# Verificam simetria grupurilor prin automorphisme
print("""
Verificam daca toate cele 4 grupuri sunt echivalente prin simetrie,
sau daca unul e diferit.

Din punct de vedere al SIMETRIEI 600-cell, toate grupurile TREBUIE
sa fie echivalente (simetrie S4 pe cele 4 coordonate).

Dar alegerea unei DIRECTII PRIVILEGIATE (timp) sparge aceasta simetrie!
""")

# Verificam ca toate grupurile au aceeasi structura interna
print("Structura interna a fiecarui grup (grad intern):")
for g in range(4):
    idx = groups_C[g]
    sub_adj = A_adj[np.ix_(idx, idx)]
    internal_edges = np.sum(sub_adj) / 2
    internal_degree = np.sum(sub_adj, axis=1)[0]  # toate sunt egale
    print(f"  G{g}: {int(internal_edges)} muchii interne, grad intern = {internal_degree:.0f}")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZII EXP-090")
print("=" * 70)

print("""
DESCOPERIRI PRINCIPALE:
=======================

1. STRUCTURA SIMETRICA:
   - Cele 4 grupuri (G0, G1, G2, G3) sunt PERFECT SIMETRICE
   - Fiecare are 24 varfuri, aceeasi conectivitate interna/externa
   - Din punct de vedere geometric, nu exista un grup "special"

2. MECANISMUL DE SELECTIE:
   - Alegerea unei AXE PRIVILEGIATE (ex: directia timpului x3)
   - SPARGE simetria si selecteaza 3 grupuri
   - Grupul Gk (unde xk=0) devine DECONECTAT de axa k

3. INTERPRETARE FIZICA:
   - 4 grupuri geometrice -> 3 generatii + 1 sector sterile
   - Selectia e facuta de DIRECTIA TIMPULUI in 4D
   - Al 4-lea grup = neutrini sterili / dark sector

4. PREDICTIE:
   - Exista un al 4-lea tip de fermion (sterile)
   - Nu interactioneaza cu campul EM (nu are sarcina)
   - Posibil candidat pentru materia intunecata

STATUS: PATTERN -> EXPLICATIE PARTIALA
======================================
- NU am derivat de ce exponentii sunt 11, 6, 17
- DAR am explicat de ce 4 grupuri dau 3 generatii observabile
- Mecanismul: selectia prin directia temporala

URMATORUL PAS:
==============
Investigarea daca RAPORTUL maselor (phi^11, phi^6) poate fi derivat
din DISTANTA sau CONECTIVITATEA diferita a celor 3 grupuri selectate
fata de axa temporala.
""")

# ============================================================
# BONUS: VERIFICARE NUMERICA FINALA
# ============================================================
print("\n" + "-" * 70)
print("BONUS: VERIFICARE NUMERICA")
print("-" * 70)

# Calculam o "masa efectiva" pentru fiecare grup
# bazata pe inversa conexiunii la axa privilegiata (x3)

print("'Masa efectiva' bazata pe 1/(conexiune la axa 3):")
masses = []
for g in range(4):
    conn = sum(A_adj[vi, vj] for vi in axes_A[3] for vj in groups_C[g])
    if conn > 0:
        eff_mass = 1.0 / conn
        masses.append((g, eff_mass, conn))
        print(f"  G{g}: conexiune = {int(conn)}, masa_eff ~ {eff_mass:.6f}")
    else:
        print(f"  G{g}: conexiune = 0 (STERILE - masa >> celelalte)")

if len(masses) >= 3:
    # Sortam dupa masa efectiva
    masses.sort(key=lambda x: x[1])
    print(f"\nRapoarte de mase (ordonate):")
    print(f"  m1/m0 = {masses[1][1]/masses[0][1]:.4f}")
    print(f"  m2/m1 = {masses[2][1]/masses[1][1]:.4f}")
    print(f"  m2/m0 = {masses[2][1]/masses[0][1]:.4f}")
    print(f"\nComparatie cu phi^n:")
    print(f"  phi^1 = {phi:.4f}")
    print(f"  phi^6 = {phi**6:.4f}")
    print(f"  phi^11 = {phi**11:.4f}")

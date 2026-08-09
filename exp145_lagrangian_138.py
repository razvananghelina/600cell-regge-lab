"""
EXP-145: Lattice Gauge Lagrangian with 1+3+8 Edge Split
=========================================================
Build explicit lattice gauge action on 600-cell using:
- 1+3+8 edge classification from exp143
- Gram matrix = 24-cell from exp144
- Wilson plaquette action decomposed by gauge sector

Goal: derive coupling constant ratios from geometry.
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict

PHI = (1 + 5**0.5) / 2
ALPHA = 1/137.036
ALPHA_S = 0.1180

def generate_600cell():
    phi = PHI
    vertices = set()
    for i in range(4):
        for s in [1, -1]:
            v = [0.0]*4; v[i] = float(s)
            vertices.add(tuple(round(x,10) for x in v))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    vertices.add((s0, s1, s2, s3))
    base_vals = [phi/2, 0.5, 1/(2*phi), 0.0]
    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
        if inv % 2 == 0:
            even_perms.append(p)
    for perm in even_perms:
        for s0 in [1,-1]:
            for s1 in [1,-1]:
                for s2 in [1,-1]:
                    signed = [s0*base_vals[0], s1*base_vals[1], s2*base_vals[2], base_vals[3]]
                    v = tuple(round(signed[perm.index(i)],10) for i in range(4))
                    vertices.add(v)
    return [list(v) for v in sorted(vertices)]

def classify_vertex(v):
    nz = sum(1 for x in v if abs(x) > 1e-9)
    if nz == 1 and any(abs(abs(x)-1.0) < 1e-9 for x in v):
        return 'A'
    if all(abs(abs(x)-0.5) < 1e-9 for x in v):
        return 'B'
    return 'C'

print("="*70)
print("EXP-145: LATTICE GAUGE LAGRANGIAN WITH 1+3+8 SPLIT")
print("="*70)

verts = generate_600cell()
N = len(verts)
vtypes = [classify_vertex(v) for v in verts]
coords = np.array(verts)

neighbors = defaultdict(set)
for i in range(N):
    for j in range(i+1, N):
        d2 = np.sum((coords[i] - coords[j])**2)
        if d2 < 1.0/PHI**2 + 1e-6:
            neighbors[i].add(j)
            neighbors[j].add(i)

A_set = set(i for i in range(N) if vtypes[i] == 'A')
B_set = set(i for i in range(N) if vtypes[i] == 'B')
C_set = set(i for i in range(N) if vtypes[i] == 'C')
gauge_set = A_set | B_set

# ============================================================
# PART 1: Classify all 720 edges in 1+3+8 basis
# ============================================================
print("\n" + "="*70)
print("PART 1: EDGE CLASSIFICATION (1+3+8)")
print("="*70)

# For each CC edge, determine if it's "special" (CC*) or regular
# CC* criterion: the two C endpoints share an A-gauge neighbor but zero B-gauge neighbors
all_edges = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            all_edges.append((i, j))

print(f"Total edges: {len(all_edges)}")

def get_gauge_type(i, j):
    """Classify edge (i,j) as U1, SU2_charged, SU2_neutral, or SU3"""
    ti, tj = vtypes[i], vtypes[j]
    if (ti == 'A' and tj == 'C') or (ti == 'C' and tj == 'A'):
        return 'U1'      # AC edge
    if (ti == 'B' and tj == 'C') or (ti == 'C' and tj == 'B'):
        return 'SU2_ch'  # BC edge (charged W)
    if ti == 'C' and tj == 'C':
        # Check if special: shared A > 0 and shared B = 0
        a_i = A_set & neighbors[i]
        a_j = A_set & neighbors[j]
        b_i = B_set & neighbors[i]
        b_j = B_set & neighbors[j]
        shared_A = len(a_i & a_j)
        shared_B = len(b_i & b_j)
        if shared_A > 0 and shared_B == 0:
            return 'SU2_n'  # CC* (neutral Z)
        else:
            return 'SU3'    # Regular CC
    # A-A, A-B, B-B (shouldn't exist)
    return 'NONE'

edge_types = {}
for i, j in all_edges:
    edge_types[(i,j)] = get_gauge_type(i, j)
    edge_types[(j,i)] = get_gauge_type(i, j)

type_counts = Counter(edge_types[(i,j)] for i,j in all_edges)
print("\nEdge classification:")
for t in ['U1', 'SU2_ch', 'SU2_n', 'SU3', 'NONE']:
    if t in type_counts:
        print(f"  {t}: {type_counts[t]}")

# Map to SM sectors
sector_map = {'U1': 'U1', 'SU2_ch': 'SU2', 'SU2_n': 'SU2', 'SU3': 'SU3'}
sector_counts = Counter()
for t, c in type_counts.items():
    sector_counts[sector_map.get(t, t)] += c

print("\nSM sector edge counts:")
for s in ['U1', 'SU2', 'SU3']:
    print(f"  {s}: {sector_counts[s]}")
print(f"  Total: {sum(sector_counts.values())}")

# Per fermion check
print("\nPer-fermion edge distribution:")
for ci in list(C_set)[:1]:
    u1 = sum(1 for j in neighbors[ci] if edge_types.get((ci,j)) == 'U1')
    su2c = sum(1 for j in neighbors[ci] if edge_types.get((ci,j)) == 'SU2_ch')
    su2n = sum(1 for j in neighbors[ci] if edge_types.get((ci,j)) == 'SU2_n')
    su3 = sum(1 for j in neighbors[ci] if edge_types.get((ci,j)) == 'SU3')
    print(f"  C vertex {ci}: U1={u1}, SU2_ch={su2c}, SU2_n={su2n}, SU3={su3}")
    print(f"  => 1 + 3 + 8 = {u1} + {su2c+su2n} + {su3}")

# ============================================================
# PART 2: Classify all 1200 triangles by gauge sector
# ============================================================
print("\n" + "="*70)
print("PART 2: TRIANGLE (PLAQUETTE) CLASSIFICATION")
print("="*70)

triangles = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            for k in neighbors[i] & neighbors[j]:
                if k > j:
                    triangles.append((i,j,k))

print(f"Total triangles: {len(triangles)}")

# Classify each triangle by the gauge types of its 3 edges
tri_gauge_types = []
for i, j, k in triangles:
    e1 = edge_types[(i,j)]
    e2 = edge_types[(i,k)]
    e3 = edge_types[(j,k)]
    sector_combo = tuple(sorted([sector_map.get(e1,e1), sector_map.get(e2,e2), sector_map.get(e3,e3)]))
    tri_gauge_types.append(sector_combo)

tri_sector_counts = Counter(tri_gauge_types)
print("\nTriangle types by gauge sector:")
for combo, cnt in sorted(tri_sector_counts.items()):
    label = '-'.join(combo)
    print(f"  {label}: {cnt}")

# Pure sector plaquettes (all 3 edges from same sector)
print("\nPure sector plaquettes:")
for sector in ['U1', 'SU2', 'SU3']:
    key = (sector, sector, sector)
    cnt = tri_sector_counts.get(key, 0)
    print(f"  {sector}-{sector}-{sector}: {cnt}")

# ============================================================
# PART 3: Wilson action weights
# ============================================================
print("\n" + "="*70)
print("PART 3: WILSON ACTION DECOMPOSITION")
print("="*70)

# In Wilson lattice gauge theory:
# S = sum_plaquettes beta_p * Re(Tr(U_p))
# where beta_p = 2N_c / g^2 for SU(N_c)
#
# On 600-cell: plaquettes = triangles
# For each SM sector, the effective beta is:
# beta_eff = (number of pure plaquettes in sector) * (coupling weight)

# Count plaquettes involving each sector
sector_plaquettes = defaultdict(int)
for combo in tri_gauge_types:
    for s in set(combo):
        sector_plaquettes[s] += 1

print("Plaquettes involving each sector:")
for s in ['U1', 'SU2', 'SU3']:
    print(f"  {s}: {sector_plaquettes[s]}")

# Weight per fermion: how many plaquettes does each fermion participate in per sector?
print("\nPer-fermion plaquette participation:")
for ci in list(C_set)[:1]:
    ci_tris = [(i,j,k) for idx, (i,j,k) in enumerate(triangles) if ci in (i,j,k)]
    per_sector = defaultdict(int)
    for i, j, k in ci_tris:
        combo = tuple(sorted([sector_map.get(edge_types[(i,j)], '?'),
                              sector_map.get(edge_types[(i,k)], '?'),
                              sector_map.get(edge_types[(j,k)], '?')]))
        for s in set(combo):
            per_sector[s] += 1
    print(f"  C vertex {ci}: {dict(per_sector)}")

# ============================================================
# PART 4: Gram matrix eigenspaces vs gauge sectors
# ============================================================
print("\n" + "="*70)
print("PART 4: GRAM EIGENSPACES AND GAUGE SECTORS")
print("="*70)

# Build Gram matrix
gauge_verts = sorted(A_set | B_set)
C_verts = sorted(C_set)
n_G = len(gauge_verts)
n_C = len(C_verts)
g_map = {v: i for i, v in enumerate(gauge_verts)}
c_map = {v: i for i, v in enumerate(C_verts)}

GF = np.zeros((n_G, n_C))
for gi, g in enumerate(gauge_verts):
    for c in neighbors[g]:
        if c in C_set:
            GF[gi][c_map[c]] = 1

G = GF @ GF.T
eigs_G, vecs_G = np.linalg.eigh(G)
idx = np.argsort(eigs_G)[::-1]
eigs_G = eigs_G[idx]
vecs_G = vecs_G[:, idx]

def group_eigs(eigs, tol=0.01):
    groups = []
    i = 0
    while i < len(eigs):
        val = eigs[i]
        mult = 1
        while i + mult < len(eigs) and abs(eigs[i+mult] - val) < tol:
            mult += 1
        groups.append((round(val, 2), mult, i))
        i += mult
    return groups

g_groups = group_eigs(eigs_G)
print("Gram eigenspaces:")
for val, mult, start in g_groups:
    print(f"  G={val:.1f} (m={mult})")

# For each eigenspace, check overlap with edge types
# An eigenvector of G lives on gauge vertices (24-dim)
# For each gauge vertex g, it has edges to C vertices classified as U1, SU2, SU3
# The eigenvector weight on g tells us how much that eigenspace "couples" to g's sector

print("\nEigenspace decomposition by vertex type and edge sector:")
A_idx = [g_map[a] for a in sorted(A_set)]
B_idx = [g_map[b] for b in sorted(B_set)]

for val, mult, start in g_groups:
    V = vecs_G[:, start:start+mult]  # 24 x mult

    # Weight on A vs B
    A_w = np.sum(V[A_idx, :]**2)
    B_w = np.sum(V[B_idx, :]**2)

    # For each gauge vertex in this eigenspace, count its edge types
    # weighted by eigenvector amplitude
    u1_coupling = 0
    su2_coupling = 0
    su3_coupling = 0

    for gi, g in enumerate(gauge_verts):
        weight = np.sum(V[gi, :]**2)  # total weight on this gauge vertex
        # Count edges from g by sector
        for c in neighbors[g]:
            if c in C_set:
                et = edge_types.get((g, c), edge_types.get((c, g), '?'))
                s = sector_map.get(et, '?')
                if s == 'U1':
                    u1_coupling += weight
                elif s == 'SU2':
                    su2_coupling += weight
                elif s == 'SU3':
                    su3_coupling += weight

    total_coupling = u1_coupling + su2_coupling + su3_coupling
    print(f"  G={val:.0f} (m={mult}): A={A_w:.3f} B={B_w:.3f} | "
          f"U1={u1_coupling/total_coupling:.3f} SU2={su2_coupling/total_coupling:.3f} "
          f"SU3={su3_coupling/total_coupling:.3f}")

# ============================================================
# PART 5: Coupling constant extraction
# ============================================================
print("\n" + "="*70)
print("PART 5: COUPLING CONSTANT RATIOS")
print("="*70)

# Method 1: From plaquette counting
# In lattice gauge theory: 1/g^2 ~ N_plaquettes / N_edges
# For each sector:
pure_U1 = tri_sector_counts.get(('U1','U1','U1'), 0)
pure_SU2 = tri_sector_counts.get(('SU2','SU2','SU2'), 0)
pure_SU3 = tri_sector_counts.get(('SU3','SU3','SU3'), 0)

n_U1 = sector_counts['U1']
n_SU2 = sector_counts['SU2']
n_SU3 = sector_counts['SU3']

print("Method 1: Plaquette/edge ratio")
print(f"  Pure plaquettes: U1={pure_U1}, SU2={pure_SU2}, SU3={pure_SU3}")
print(f"  Edges: U1={n_U1}, SU2={n_SU2}, SU3={n_SU3}")
for s, pp, ne in [('U1', pure_U1, n_U1), ('SU2', pure_SU2, n_SU2), ('SU3', pure_SU3, n_SU3)]:
    ratio = pp / ne if ne > 0 else 0
    print(f"  {s}: plaq/edge = {ratio:.4f}")

# Method 2: From Gram eigenvalues
print("\nMethod 2: Gram eigenvalue ratios")
gram_vals = [val for val, _, _ in g_groups if val > 0.01]
print(f"  Non-zero Gram eigenvalues: {gram_vals}")
if len(gram_vals) >= 2:
    print(f"  Ratios: {[gram_vals[0]/v for v in gram_vals]}")
    # 36/6 = 6, 36/24 = 1.5, 36/12 = 3
    # alpha_s/alpha = 10*phi = 16.18
    print(f"  36/6 = {36/6:.2f} (target alpha_s/alpha = {1/ALPHA_S * ALPHA:.2f} = 1/alpha_s * alpha)")
    print(f"  alpha_s/alpha = {ALPHA_S/ALPHA:.2f} = 10*phi = {10*PHI:.2f}")
    print(f"  1/alpha_s = {1/ALPHA_S:.2f}")
    print(f"  1/alpha = {1/ALPHA:.2f}")

# Method 3: From eigenvalue * multiplicity (= trace contribution)
print("\nMethod 3: Eigenvalue x multiplicity (trace contribution)")
for val, mult, _ in g_groups:
    print(f"  G={val:.0f} x {mult} = {val*mult:.0f}")
total_tr = sum(val*mult for val, mult, _ in g_groups)
print(f"  Total trace = {total_tr:.0f} (should be 24*12 = {24*12})")

# Method 4: Effective coupling from per-sector trace
print("\nMethod 4: Sector-specific effective coupling")
# The Gram matrix block structure:
# G_ij = number of shared C-neighbors between gauge vertices i and j
# For each sector, the "coupling strength" is proportional to
# the number of fermion loops (shared C vertices)

# Build sector-projected Gram matrices
# Project GF onto edges of each sector type
for sector in ['U1', 'SU2', 'SU3']:
    GF_sector = np.zeros((n_G, n_C))
    for gi, g in enumerate(gauge_verts):
        for c in neighbors[g]:
            if c in C_set:
                et = edge_types.get((g, c), edge_types.get((c, g), '?'))
                s = sector_map.get(et, '?')
                if s == sector:
                    GF_sector[gi][c_map[c]] = 1

    G_sector = GF_sector @ GF_sector.T
    tr = np.trace(G_sector)
    n_entries = int(np.sum(G_sector > 0.5))
    eigs_s = sorted(np.linalg.eigvalsh(G_sector), reverse=True)
    nz = [round(e, 2) for e in eigs_s if abs(e) > 0.01]
    print(f"\n  {sector} sector Gram:")
    print(f"    Trace = {tr:.0f}")
    print(f"    Non-zero eigenvalues: {nz[:10]}")
    if nz:
        print(f"    Top eigenvalue: {nz[0]:.2f}")
        print(f"    Sum eigenvalues: {sum(nz):.2f}")

# ============================================================
# PART 6: The lattice gauge Lagrangian
# ============================================================
print("\n" + "="*70)
print("PART 6: EXPLICIT LAGRANGIAN STRUCTURE")
print("="*70)

print("""
The lattice gauge Lagrangian on the 600-cell:

    S = S_U1 + S_SU2 + S_SU3 + S_fermion

where each gauge sector action is a Wilson plaquette sum:

    S_G = beta_G * sum_{triangles in G} Re(Tr(U_triangle))

with U_triangle = U_e1 * U_e2 * U_e3^dag (product of link variables)
""")

# Count mixed plaquettes
print("Plaquette decomposition:")
for combo, cnt in sorted(tri_sector_counts.items()):
    label = '-'.join(combo)
    purity = "PURE" if len(set(combo)) == 1 else "MIXED"
    print(f"  {label}: {cnt} ({purity})")

total_pure = pure_U1 + pure_SU2 + pure_SU3
total_mixed = len(triangles) - total_pure
print(f"\n  Pure plaquettes: {total_pure}")
print(f"  Mixed plaquettes: {total_mixed}")
print(f"  Pure fraction: {total_pure/len(triangles):.4f}")

# ============================================================
# PART 7: Ratio analysis - what gives alpha_s/alpha?
# ============================================================
print("\n" + "="*70)
print("PART 7: COUPLING RATIO ANALYSIS")
print("="*70)

# Target: alpha_s/alpha = 10*phi = 16.18
# Or equivalently: g_s^2/g_em^2 = alpha_s/alpha

# Various geometric ratios
print("Geometric ratios to test:")

ratios = {
    "n_SU3/n_U1": n_SU3/n_U1,
    "n_SU2/n_U1": n_SU2/n_U1,
    "n_SU3/n_SU2": n_SU3/n_SU2,
    "pure_SU3/pure_U1": pure_SU3/pure_U1 if pure_U1 > 0 else float('inf'),
    "pure_SU3/pure_SU2": pure_SU3/pure_SU2 if pure_SU2 > 0 else float('inf'),
    "sector_plaq_SU3/U1": sector_plaquettes['SU3']/sector_plaquettes['U1'],
    "sector_plaq_SU3/SU2": sector_plaquettes['SU3']/sector_plaquettes['SU2'],
    "Gram 36/6": 36/6,
    "Gram 36/12": 36/12,
    "Gram 24/6": 24/6,
    "8/1 (SU3 dim / U1 dim)": 8/1,
    "8*3/1 (from mults)": 8*3,
}

target = ALPHA_S / ALPHA
print(f"\nTarget: alpha_s/alpha = {target:.4f} = 10*phi = {10*PHI:.4f}")
print(f"        alpha/alpha_s = {ALPHA/ALPHA_S:.6f}")
print(f"        1/alpha = {1/ALPHA:.2f}")
print(f"        1/alpha_s = {1/ALPHA_S:.2f}")

for name, val in sorted(ratios.items(), key=lambda x: abs(x[1] - target)):
    err = abs(val - target) / target * 100
    print(f"  {name:35s} = {val:8.4f}  (err = {err:.1f}%)")

# Also check for formulas involving phi
print("\nPhi-based combinations:")
phi_combos = {
    "10*phi": 10*PHI,
    "8*phi^2": 8*PHI**2,
    "8*phi + 2": 8*PHI + 2,
    "n_SU3/n_U1 * phi": n_SU3/n_U1 * PHI,
    "Gram_36 * Gram_6 / Gram_12": 36*6/12,
    "(SU3 edges) / (U1 edges) / phi^2": n_SU3/n_U1/PHI**2,
    "sector_plaq_SU3 / sector_plaq_U1 / phi": sector_plaquettes['SU3']/sector_plaquettes['U1']/PHI,
}

for name, val in sorted(phi_combos.items(), key=lambda x: abs(x[1] - target)):
    err = abs(val - target) / target * 100
    print(f"  {name:45s} = {val:8.4f}  (err = {err:.1f}%)")

# ============================================================
# PART 8: The beta functions from plaquette structure
# ============================================================
print("\n" + "="*70)
print("PART 8: EFFECTIVE BETA FUNCTIONS")
print("="*70)

# In lattice gauge theory, the coupling runs with scale
# At one loop: 1/g^2(mu) = 1/g^2(mu_0) + b * ln(mu/mu_0)
# where b = beta_0 / (8*pi^2)
# For SU(N): beta_0 = (11*N - 2*N_f)/3
# For U(1): beta_0 = -2*N_f/3 (QED: no gluon contribution)

# On the 600-cell, the "scale" is the graph distance d
# The hop-invariant eigenvalue 6 (m=8) means the SU(3) coupling
# does NOT run with graph distance! This IS confinement.

print("SM beta functions (one-loop):")
print(f"  SU(3): beta_0 = (11*3 - 2*6)/3 = {(11*3 - 2*6)/3:.1f} (N_f=6 quarks)")
print(f"  SU(2): beta_0 = (11*2 - 2*6)/3 = {(11*2 - 2*6)/3:.1f} (N_f=6 doublets)")
print(f"  U(1):  beta_0 = -2*11/3 = {-2*11/3:.1f} (N_f=11 charged)")

print("\n600-cell effective beta (from Gram hop behavior):")
print("  Recall: G_n eigenvalues for each eigenspace:")

A_CC = np.zeros((n_C, n_C))
for c in C_verts:
    for nb in neighbors[c]:
        if nb in C_set:
            A_CC[c_map[c]][c_map[nb]] = 1

g0_vals = [36, 24, 12, 6, 0]
g0_mults = [1, 4, 9, 8, 2]

for n_hop in range(5):
    G_n = GF @ np.linalg.matrix_power(A_CC, n_hop) @ GF.T
    gn_eigs = sorted(np.linalg.eigvalsh(G_n), reverse=True)
    gn_groups = group_eigs(gn_eigs, tol=0.5)
    gn_vals = [v for v, _, _ in gn_groups]
    print(f"  n={n_hop}: eigenvalues = {gn_vals}")

print("\nEffective running (G_n / G_0):")
for n_hop in range(1, 5):
    G_n = GF @ np.linalg.matrix_power(A_CC, n_hop) @ GF.T
    gn_eigs = sorted(np.linalg.eigvalsh(G_n), reverse=True)
    gn_groups = group_eigs(gn_eigs, tol=0.5)
    gn_vals = [v for v, _, _ in gn_groups]
    if len(gn_vals) == len(g0_vals):
        running = []
        for i in range(len(g0_vals)):
            if abs(g0_vals[i]) > 0.01:
                running.append(round(gn_vals[i]/g0_vals[i], 4))
            else:
                running.append('inf')
        print(f"  n={n_hop}: {running}")

print("\nKey observation:")
print("  Eigenvalue 6 (m=8): hop-invariant => SU(3) coupling DOES NOT RUN")
print("  Eigenvalue 0 (m=2): stays zero => two modes completely decoupled")
print("  Eigenvalue 36 (m=1): grows as 9^n => strongly coupled mode")
print("  This is OPPOSITE to SM running (where SU(3) runs most)")
print("  But on finite graph, 'confinement' = no running = fixed coupling")

# ============================================================
# PART 9: Spectral action computation
# ============================================================
print("\n" + "="*70)
print("PART 9: SPECTRAL ACTION")
print("="*70)

# Connes-Chamseddine spectral action: S = Tr(f(D^2/Lambda^2))
# On a lattice: D is the graph Dirac operator, f is a cutoff function
# For the gauge sector: S_gauge = sum_k f(lambda_k) where lambda_k are
# Laplacian eigenvalues

# The 600-cell Laplacian L = 12*I - A has eigenvalues 12 - lambda_adj
# Gauge sector Laplacian: L_gauge = 8*I - A_eff (degree 8 for 24-cell)

L_eff = 8*np.eye(n_G) - (G - 12*np.eye(n_G))/3  # = 8I - A_eff
L_eigs = sorted(np.linalg.eigvalsh(L_eff))

print("Gauge Laplacian (L = 8I - A_eff) eigenvalues:")
l_groups = group_eigs(L_eigs, tol=0.01)
for val, mult, _ in l_groups:
    print(f"  {val:.2f} x {mult}")

# Spectral action for heat kernel: Tr(exp(-t*L))
print("\nHeat kernel Tr(exp(-t*L_gauge)):")
for t in [0.1, 0.5, 1.0, 2.0, 5.0]:
    hk = sum(mult * np.exp(-t * val) for val, mult, _ in l_groups)
    print(f"  t={t:.1f}: {hk:.4f}")

# The ratio of heat kernels at different t gives effective coupling evolution
print("\nHeat kernel ratio (gauge sector):")
for t in [0.5, 1.0, 2.0]:
    hk = sum(mult * np.exp(-t * val) for val, mult, _ in l_groups)
    hk0 = 24  # t=0 limit
    print(f"  t={t:.1f}: Tr/24 = {hk/24:.4f}")

# ============================================================
# PART 10: The 5/12 closure and Lagrangian
# ============================================================
print("\n" + "="*70)
print("PART 10: CLOSURE FRACTION AND YANG-MILLS STRUCTURE")
print("="*70)

# From exp142: bilinear Lie algebra closure = 5/12
# From exp143: local CC density = 5/12
# 5/12 = 1 - 7/12
# In SM: the non-Abelian fraction of the gauge group:
# dim(SU2) + dim(SU3) / dim(SM) = (3+8)/12 = 11/12
# Or: dim(SU3)/dim(SM) = 8/12 = 2/3
# Or: dim(SU2)/dim(SM) = 3/12 = 1/4

# The 5/12 could be: mixed sector fraction
mixed_frac = total_mixed / len(triangles)
pure_frac = total_pure / len(triangles)
print(f"Mixed plaquette fraction: {mixed_frac:.4f}")
print(f"Pure plaquette fraction: {pure_frac:.4f}")
print(f"5/12 = {5/12:.4f}")
print(f"7/12 = {7/12:.4f}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"""
EDGE CLASSIFICATION (1+3+8):
  U(1):  {type_counts.get('U1', 0)} edges (96 AC)
  SU(2): {type_counts.get('SU2_ch', 0) + type_counts.get('SU2_n', 0)} edges (192 BC + 48 CC*)
  SU(3): {type_counts.get('SU3', 0)} edges (384 CC_regular)
  Total: {sum(type_counts.values())} = 720

PLAQUETTES:
  Pure U(1):  {pure_U1}
  Pure SU(2): {pure_SU2}
  Pure SU(3): {pure_SU3}
  Mixed:      {total_mixed}
  Total:      {len(triangles)} = 1200

GRAM MATRIX:
  Spectrum: {{36, 24, 12, 6, 0}} x {{1, 4, 9, 8, 2}}
  A_eff = A_24cell = A^2/3

COUPLING RATIOS:
  SU3/U1 edges = {n_SU3/n_U1:.2f}
  SU3/SU2 edges = {n_SU3/n_SU2:.2f}
  Target alpha_s/alpha = {ALPHA_S/ALPHA:.2f} = 10*phi

HOP BEHAVIOR:
  Eigenvalue 6 (m=8): HOP-INVARIANT (SU(3) confinement?)
  Eigenvalue 36 (m=1): grows as 9^n
  Eigenvalue 0 (m=2): stays zero
""")

print("STATUS: MIXED")
print("  - Edge 1+3+8 classification: DERIVED")
print("  - Plaquette decomposition: DERIVED")
print("  - Pure vs mixed plaquettes: DERIVED")
print("  - Coupling ratio alpha_s/alpha from geometry: NOT YET DERIVED")
print("  - Hop-invariant eigenvalue = confinement: PATTERN")
print("  - Lagrangian structure: DERIVED (form), PATTERN (coupling values)")

print("\n" + "="*70)
print("END EXP-145")
print("="*70)

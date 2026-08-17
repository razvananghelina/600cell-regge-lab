"""
EXP-313c: Combine Hopf Fibration x Gauge Decomposition
=======================================================
GOAL: Cross the TWO independent edge classifications:
  1. Hopf: 3 classes x 240 edges (fiber = 0, 1/phi^2, 1/phi)
  2. Gauge: 1+3+8 per vertex (from A5 irrep projectors on 12 neighbors)

ALSO: Fix vertex classification to A=2 (poles), B=22 (other axis+half), C=96 (golden)
consistent with the paper's spectral action section.

KEY TEST: If the 3x3 contingency table (Hopf_class x gauge_sector) is NOT
proportional to the product of marginals, we have genuine anisotropy --
and that could produce the sector-dependent holonomy corrections.

Compute: <f>_U(1), <f>_SU(2), <f>_SU(3) -- mean fiber fraction per gauge sector.

NOTE: All print uses ASCII only (Windows cp1252).
"""

import numpy as np
import math
from itertools import product as iterproduct
from collections import Counter, defaultdict

PHI = (1 + math.sqrt(5)) / 2
PHI_P = (1 - math.sqrt(5)) / 2
SQRT5 = math.sqrt(5)
PI = math.pi
a1 = 5
b1 = 6
N = 120
K3P = "3'"

# Framework constants
ALPHA_S_FW = 1 / (2 * PHI**3)
SIN2_TW_FW = b1 / (a1**2 + 1)
_disc = (4*a1*PHI**4)**2 - 8*PI
ALPHA_FW = (4*a1*PHI**4 - math.sqrt(_disc)) / (4*PI)

print("=" * 72)
print("EXP-313c: HOPF FIBRATION x GAUGE DECOMPOSITION")
print("=" * 72)

# =====================================================================
# SECTION 1: Build 600-cell
# =====================================================================
print("\n--- SECTION 1: Build 600-cell ---")

def qmul(a, b):
    return np.array([
        a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
        a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0]
    ])

def qconj(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

def qinv(a):
    return qconj(a) / np.dot(a, a)

def generate_2I():
    elements = set()
    def add(q):
        q = tuple(round(x, 10) for x in q)
        elements.add(q)
        elements.add(tuple(-x for x in q))
    for i in range(4):
        v = [0,0,0,0]; v[i] = 1.0; add(v)
    for s0 in [0.5,-0.5]:
        for s1 in [0.5,-0.5]:
            for s2 in [0.5,-0.5]:
                for s3 in [0.5,-0.5]:
                    add([s0,s1,s2,s3])
    abs_vals = [0.0, 0.5, abs(PHI_P)/2, PHI/2]
    even_perms = [
        (0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
        (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0),
    ]
    for perm in even_perms:
        base = [abs_vals[perm[k]] for k in range(4)]
        non_zero = [i for i in range(4) if abs(base[i]) > 1e-10]
        for signs in iterproduct([1,-1], repeat=len(non_zero)):
            v = list(base)
            for idx, s in zip(non_zero, signs):
                v[idx] *= s
            add(v)
    return [np.array(q) for q in sorted(elements)]

elements_2I = generate_2I()
V = np.array(elements_2I)
n_vert = len(elements_2I)
print(f"  |2I| = {n_vert}")
assert n_vert == 120

# Adjacency
IP = V @ V.T
A_adj = np.zeros((n_vert, n_vert), dtype=int)
for i in range(n_vert):
    for j in range(i+1, n_vert):
        if abs(IP[i,j] - PHI/2) < 0.01:
            A_adj[i,j] = A_adj[j,i] = 1
n_edges = np.sum(A_adj) // 2
degree = int(np.sum(A_adj[0]))
print(f"  Edges = {n_edges}, Degree = {degree}")

def find_vertex_index(q, verts, tol=1e-6):
    dists = np.sum((verts - q)**2, axis=1)
    idx = np.argmin(dists)
    return idx if dists[idx] < tol else -1

# =====================================================================
# SECTION 2: Corrected vertex classification (A=2, B=22, C=96)
# =====================================================================
print("\n--- SECTION 2: Vertex classification (A=2, B=22, C=96) ---")

# A = 2 real-axis POLES: (1,0,0,0) and (-1,0,0,0)
# B = 22 remaining non-golden: 6 other axis + 16 half-integer
# C = 96 golden ratio vertices

def classify_vertex_ABC(q):
    """Corrected classification matching paper's spectral action."""
    abs_coords = [abs(x) for x in q]
    n_zero = sum(1 for x in abs_coords if x < 0.01)

    # Type A: real-axis poles (+/-1, 0, 0, 0) only
    if n_zero == 3:
        max_coord = max(abs_coords)
        if abs(max_coord - 1.0) < 0.01:
            # Check it's the q0 axis (real part)
            if abs(abs(q[0]) - 1.0) < 0.01:
                return 'A'
            else:
                return 'B'  # Other axis vertices like (0,1,0,0)

    # Half-integer vertices: all |coords| = 0.5
    n_half = sum(1 for x in abs_coords if abs(x - 0.5) < 0.01)
    if n_half == 4:
        return 'B'

    # Remaining are golden
    return 'C'

vtypes = [classify_vertex_ABC(V[i]) for i in range(n_vert)]
type_counts = Counter(vtypes)
print(f"  Vertex types: {dict(sorted(type_counts.items()))}")
print(f"  Expected: A=2, B=22, C=96")

# Verify A vertices
A_indices = [i for i in range(n_vert) if vtypes[i] == 'A']
for i in A_indices:
    print(f"    A vertex {i}: {V[i]}")

# Check: A vertices should be (1,0,0,0) and (-1,0,0,0)
assert type_counts['A'] == 2, f"Expected 2 A vertices, got {type_counts['A']}"
assert type_counts['B'] == 22, f"Expected 22 B vertices, got {type_counts['B']}"
assert type_counts['C'] == 96, f"Expected 96 C vertices, got {type_counts['C']}"
print("  Classification VERIFIED: A=2, B=22, C=96")

# Neighbor type distribution for each vertex type
print(f"\n  Neighbor type distribution:")
for vtype in ['A', 'B', 'C']:
    v_indices = [i for i in range(n_vert) if vtypes[i] == vtype]
    neighbor_type_counts = Counter()
    for i in v_indices:
        for j in range(n_vert):
            if A_adj[i, j]:
                neighbor_type_counts[vtypes[j]] += 1
    # Average per vertex
    n_v = len(v_indices)
    avg_neighbors = {t: neighbor_type_counts[t] / n_v for t in ['A', 'B', 'C']}
    print(f"    {vtype} ({n_v:3d} vertices): "
          f"<n_A>={avg_neighbors.get('A',0):.2f}, "
          f"<n_B>={avg_neighbors.get('B',0):.2f}, "
          f"<n_C>={avg_neighbors.get('C',0):.2f}, "
          f"total={sum(avg_neighbors.values()):.2f}")

# =====================================================================
# SECTION 3: Build edge list with vertex-type classification
# =====================================================================
print("\n--- SECTION 3: Edge classification by vertex type ---")

edge_list = []
edge_vtype = []  # (type_i, type_j) sorted
for i in range(n_vert):
    for j in range(i+1, n_vert):
        if A_adj[i, j]:
            edge_list.append((i, j))
            edge_vtype.append(tuple(sorted([vtypes[i], vtypes[j]])))

print(f"  Total edges: {len(edge_list)}")
edge_type_counts = Counter(edge_vtype)
print(f"\n  Edge types by vertex classification:")
for etype, cnt in sorted(edge_type_counts.items()):
    print(f"    {etype}: {cnt}")

# =====================================================================
# SECTION 4: Hopf fiber/base decomposition (from exp313b)
# =====================================================================
print("\n--- SECTION 4: Hopf fiber/base decomposition ---")

def hopf_fiber_tangent(q):
    """Fiber tangent at q (using i as generator)."""
    return np.array([-q[1], q[0], -q[3], q[2]])

fiber_fractions = []
fiber_weights = []
base_weights = []

for i_e, (i, j) in enumerate(edge_list):
    delta = V[j] - V[i]
    delta_tang = delta - np.dot(delta, V[i]) * V[i]
    norm_tang = np.linalg.norm(delta_tang)
    if norm_tang < 1e-10:
        fiber_weights.append(0.0)
        base_weights.append(0.0)
        fiber_fractions.append(0.0)
        continue

    delta_hat = delta_tang / norm_tang
    T = hopf_fiber_tangent(V[i])
    fc = abs(np.dot(delta_hat, T))
    bc = math.sqrt(max(0, 1 - fc**2))
    fiber_weights.append(fc)
    base_weights.append(bc)
    # Fiber fraction = fc / (fc + bc) if fc+bc > 0 else 0
    total = fc + bc
    fiber_fractions.append(fc / total if total > 0.001 else 0.0)

fiber_weights = np.array(fiber_weights)
base_weights = np.array(base_weights)
fiber_fractions = np.array(fiber_fractions)

# Classify into Hopf classes
fw_rounded = np.round(fiber_weights, 4)
distinct_fw = sorted(set(fw_rounded))
print(f"  Distinct fiber weights: {len(distinct_fw)}")
for val in distinct_fw:
    cnt = np.sum(np.abs(fw_rounded - val) < 0.001)
    # Compute fiber fraction for this class
    mask = np.abs(fw_rounded - val) < 0.001
    ff_vals = fiber_fractions[mask]
    print(f"    f_fiber = {val:.4f}: {cnt} edges, fiber_fraction = {ff_vals.mean():.6f}")

# Assign Hopf class label (0, 1, 2 for the 3 groups)
hopf_class = np.zeros(len(edge_list), dtype=int)
for i_e in range(len(edge_list)):
    if fw_rounded[i_e] < 0.001:
        hopf_class[i_e] = 0  # purely base
    elif abs(fw_rounded[i_e] - 0.5257) < 0.01:
        hopf_class[i_e] = 1  # mixed (1/phi^2)
    else:
        hopf_class[i_e] = 2  # fiber-dominant (1/phi)

hopf_class_names = {0: 'base(f=0)', 1: 'mixed(1/phi^2)', 2: 'fiber(1/phi)'}
for hc in range(3):
    cnt = np.sum(hopf_class == hc)
    print(f"    Hopf class {hc} [{hopf_class_names[hc]:16s}]: {cnt} edges")

# =====================================================================
# SECTION 5: A5 conjugation action & irrep projectors (from exp313)
# =====================================================================
print("\n--- SECTION 5: A5 irrep projectors on 12 neighbors ---")

# Find identity and its neighbors
idx_e = find_vertex_index(np.array([1.0, 0, 0, 0]), V)
neighbors_e = [j for j in range(n_vert) if A_adj[idx_e, j]]
print(f"  Identity: vertex {idx_e}, neighbors: {len(neighbors_e)}")

# Build A5 quotient
neg_partner = np.zeros(n_vert, dtype=int)
for i in range(n_vert):
    neg_partner[i] = find_vertex_index(-V[i], V)

a5_reps = []
used = set()
for i in range(n_vert):
    if i not in used:
        a5_reps.append(i)
        used.add(i)
        used.add(neg_partner[i])

# Conjugation permutations on 12 neighbors
conj_matrices = []
for g_idx in a5_reps:
    g = elements_2I[g_idx]
    g_inv = qinv(g)
    perm = np.zeros(12, dtype=int)
    for ni, h_idx in enumerate(neighbors_e):
        h = elements_2I[h_idx]
        ghg = qmul(qmul(g, h), g_inv)
        found = False
        for nj, k_idx in enumerate(neighbors_e):
            if np.sum((ghg - elements_2I[k_idx])**2) < 1e-8:
                perm[ni] = nj
                found = True
                break
        if not found:
            g2 = -g
            g2_inv = qinv(g2)
            ghg2 = qmul(qmul(g2, h), g2_inv)
            for nj, k_idx in enumerate(neighbors_e):
                if np.sum((ghg2 - elements_2I[k_idx])**2) < 1e-8:
                    perm[ni] = nj
                    found = True
                    break
        assert found
    P = np.zeros((12, 12))
    for i_p in range(12):
        P[perm[i_p], i_p] = 1
    conj_matrices.append(P)

# Classify A5 elements
def a5_order(g_idx):
    g = elements_2I[g_idx]
    power = np.array([1.0, 0, 0, 0])
    for k in range(1, 11):
        power = qmul(power, g)
        if abs(abs(power[0]) - 1) < 1e-8 and np.sum(power[1:]**2) < 1e-12:
            return k
    return -1

cls_indices = {'e': [], '2': [], '3': [], '5a': [], '5b': []}
for ci, g_idx in enumerate(a5_reps):
    o = a5_order(g_idx)
    if o == 1: cls_indices['e'].append(ci)
    elif o == 2: cls_indices['2'].append(ci)
    elif o == 3: cls_indices['3'].append(ci)
    elif o in (5, 10):
        chi2 = abs(2 * elements_2I[g_idx][0])
        if abs(chi2 - PHI) < 0.1: cls_indices['5a'].append(ci)
        elif abs(chi2 - abs(PHI_P)) < 0.1: cls_indices['5b'].append(ci)

# Build irrep projectors
chi_A5 = {
    '1':  {'e': 1, '2': 1,  '3': 1,  '5a': 1,    '5b': 1},
    '3':  {'e': 3, '2': -1, '3': 0,  '5a': PHI,  '5b': PHI_P},
    K3P:  {'e': 3, '2': -1, '3': 0,  '5a': PHI_P, '5b': PHI},
    '4':  {'e': 4, '2': 0,  '3': 1,  '5a': -1,   '5b': -1},
    '5':  {'e': 5, '2': 1,  '3': -1, '5a': 0,    '5b': 0},
}
class_sizes = {'e': 1, '2': 15, '3': 20, '5a': 12, '5b': 12}
order_A5 = 60

cls_avg_matrices = {}
for cls_name, cls_elts in cls_indices.items():
    cls_avg_matrices[cls_name] = sum(conj_matrices[ci] for ci in cls_elts) / len(cls_elts)

P_12 = {}
for R_name, R_chi in chi_A5.items():
    d_R = R_chi['e']
    P = np.zeros((12, 12))
    for cls_name in cls_indices:
        P += class_sizes[cls_name] * R_chi[cls_name] * cls_avg_matrices[cls_name]
    P *= d_R / order_A5
    P_12[R_name] = P

# Print edge weights (diagonal of projectors)
print(f"\n  Edge weights w_R for each of 12 neighbors of identity:")
print(f"  {'ni':>3s} {'w_1':>8s} {'w_3':>8s} {'w_3p':>8s} {'w_4':>8s} {'w_5':>8s} | {'type':>5s}")
for ni in range(12):
    h_idx = neighbors_e[ni]
    w = {R: P_12[R][ni, ni] for R in P_12}
    print(f"  {ni:3d} {w['1']:8.5f} {w['3']:8.5f} {w[K3P]:8.5f} {w['4']:8.5f} {w['5']:8.5f} | {vtypes[h_idx]:>5s}")

# Check: are all weights really isotropic?
w_vals = {R: [] for R in ['1', '3', K3P, '4', '5']}
for ni in range(12):
    for R in w_vals:
        w_vals[R].append(P_12[R][ni, ni])

print(f"\n  Isotropy check (std of diagonal elements):")
for R in ['1', '3', K3P, '4', '5']:
    vals = w_vals[R]
    print(f"    w_{R:2s}: mean={np.mean(vals):.6f}, std={np.std(vals):.2e}, "
          f"d_R/12={chi_A5[R]['e']/12:.6f}")

# =====================================================================
# SECTION 6: Assign gauge sector to each edge
# =====================================================================
print("\n--- SECTION 6: Gauge sector assignment ---")

# METHOD A: From A5 projector weights (dominant irrep per edge)
# Since w_R = d_R/12 for all edges (isotropic), this doesn't work.
# Every edge has the same decomposition.

# METHOD B: From vertex type (A, B, C)
# Map: edge (i,j) where vtypes are:
#   (A, C): A vertex is pole -> U(1) singlet direction (1 neighbor type)
#   (B, C): B vertex is intermediate -> SU(2) (3 neighbors)
#   (C, C): golden-golden -> SU(3) (8 neighbors)
# But this doesn't quite work because each C vertex may have different
# numbers of A, B, C neighbors.

# METHOD C: Use GLOBAL sub-adjacency matrices from exp313
# Build A_sub[R] for R in {1, 3, 3', 5}
# Then for each edge (i,j), read off A_sub[R][i,j]
# This IS the correct gauge weight per edge.

print("  Building global sub-adjacency matrices...")

A_sub = {}
for R in ['1', '3', K3P, '5']:
    A_sub[R] = np.zeros((n_vert, n_vert))

for i in range(n_vert):
    g = elements_2I[i]
    g_inv = qinv(g)
    for j in range(n_vert):
        if A_adj[i, j]:
            h = qmul(g_inv, elements_2I[j])
            h_idx = find_vertex_index(h, V)
            ni = -1
            if h_idx in neighbors_e:
                ni = neighbors_e.index(h_idx)
            elif neg_partner[h_idx] in neighbors_e:
                ni = neighbors_e.index(neg_partner[h_idx])
            else:
                h2 = qmul(-g_inv, elements_2I[j])
                h2_idx = find_vertex_index(h2, V)
                if h2_idx in neighbors_e:
                    ni = neighbors_e.index(h2_idx)
            if ni >= 0:
                for R in ['1', '3', K3P, '5']:
                    A_sub[R][i, j] = P_12[R][ni, ni]

# Verify
A_sum_check = sum(A_sub[R] for R in A_sub)
err = np.max(np.abs(A_sum_check - A_adj))
print(f"  ||A_1 + A_3 + A_3' + A_5 - A|| = {err:.2e}")

# Now for each edge in edge_list, extract the gauge weights
gauge_weights = {R: np.zeros(len(edge_list)) for R in ['1', '3', K3P, '5']}
for i_e, (i, j) in enumerate(edge_list):
    for R in gauge_weights:
        # Use average of both directions
        gauge_weights[R][i_e] = (A_sub[R][i, j] + A_sub[R][j, i]) / 2

# Verify gauge weights sum to 1 for each edge
gauge_sum = sum(gauge_weights[R] for R in gauge_weights)
print(f"  Gauge weight sum per edge: mean={gauge_sum.mean():.6f}, std={gauge_sum.std():.2e}")

# Map to gauge sectors: U(1) = 1, SU(2) = 3' (adjoint), SU(3) = 3+5
# From A5 rep theory: 12 = 1 + 3' + (3 + 5) = 1 + 3 + 8
gauge_sector_weights = {}
gauge_sector_weights['U1'] = gauge_weights['1']
gauge_sector_weights['SU2'] = gauge_weights[K3P]
gauge_sector_weights['SU3'] = gauge_weights['3'] + gauge_weights['5']

print(f"\n  Gauge sector weights per edge:")
for sector in ['U1', 'SU2', 'SU3']:
    vals = gauge_sector_weights[sector]
    print(f"    {sector:4s}: mean={vals.mean():.6f}, std={vals.std():.2e}")

# =====================================================================
# SECTION 7: CONTINGENCY TABLE - Hopf class x Gauge sector
# =====================================================================
print("\n--- SECTION 7: Contingency table (Hopf x Gauge) ---")

# Since gauge weights are continuous and ISOTROPIC (all edges have same w_R),
# the contingency table with continuous weights will be trivially proportional.
# But let's verify this explicitly and look for any deviation.

print(f"\n  APPROACH A: Mean gauge weight per Hopf class")
print(f"  (If isotropic, all Hopf classes have same gauge weights)")
print(f"  {'Hopf class':>20s} | {'w_U1':>10s} {'w_SU2':>10s} {'w_SU3':>10s} | {'n_edges':>7s}")
print(f"  {'-'*72}")

for hc in range(3):
    mask = hopf_class == hc
    n_hc = np.sum(mask)
    if n_hc > 0:
        w_u1 = gauge_sector_weights['U1'][mask].mean()
        w_su2 = gauge_sector_weights['SU2'][mask].mean()
        w_su3 = gauge_sector_weights['SU3'][mask].mean()
        print(f"  {hopf_class_names[hc]:>20s} | {w_u1:10.6f} {w_su2:10.6f} {w_su3:10.6f} | {n_hc:7d}")

print(f"\n  VERDICT: The gauge weights are ISOTROPIC from A5 projectors.")
print(f"  w_R = d_R/12 for ALL edges regardless of Hopf class.")
print(f"  This confirms exp313's negative result at the projector level.")

# =====================================================================
# SECTION 8: APPROACH B - Vertex-type edge classification x Hopf
# =====================================================================
print("\n\n--- SECTION 8: Vertex-type edges x Hopf (A=2, B=22, C=96) ---")

# For each edge, classify by BOTH:
#   - endpoint types: (A,B), (A,C), (B,B), (B,C), (C,C)
#   - Hopf class: 0, 1, 2

contingency = defaultdict(int)
fiber_by_vtype = defaultdict(list)

for i_e, (i, j) in enumerate(edge_list):
    et = tuple(sorted([vtypes[i], vtypes[j]]))
    hc = hopf_class[i_e]
    contingency[(et, hc)] += 1
    fiber_by_vtype[et].append(fiber_fractions[i_e])

# Print contingency table
vtypes_edges = sorted(set(edge_vtype))
print(f"\n  Contingency table: (vertex type) x (Hopf class)")
print(f"  {'Edge type':>12s} | {'base(f=0)':>10s} {'mixed':>10s} {'fiber':>10s} | {'total':>7s}")
print(f"  {'-'*60}")
for et in vtypes_edges:
    counts = [contingency.get((et, hc), 0) for hc in range(3)]
    total = sum(counts)
    print(f"  {str(et):>12s} | {counts[0]:10d} {counts[1]:10d} {counts[2]:10d} | {total:7d}")

# Marginals
col_totals = [sum(contingency.get((et, hc), 0) for et in vtypes_edges) for hc in range(3)]
print(f"  {'Total':>12s} | {col_totals[0]:10d} {col_totals[1]:10d} {col_totals[2]:10d} | {sum(col_totals):7d}")

# Mean fiber fraction per vertex-type edge
print(f"\n  Mean fiber fraction per vertex-type edge:")
for et in vtypes_edges:
    vals = fiber_by_vtype[et]
    if vals:
        print(f"    {str(et):12s}: <f> = {np.mean(vals):.6f} +- {np.std(vals):.6f} (n={len(vals)})")

# =====================================================================
# SECTION 9: Physical gauge sectors via vertex type
# =====================================================================
print("\n--- SECTION 9: Physical gauge sectors from vertex types ---")

# The PHYSICAL interpretation of vertex types:
# A vertices (2 poles) = distinguished by the Hopf fibration axis
# B vertices (22) = intermediate (axis + half-integer)
# C vertices (96) = generic golden vertices
#
# The edge types map to gauge sectors:
# Idea: Each C vertex has 12 neighbors. In the 1+3+8 decomposition:
#   1 = U(1): the edge pointing toward the closest pole (A vertex)
#   3 = SU(2): edges pointing toward B vertices
#   8 = SU(3): edges between C vertices
#
# This gives a DISCRETE gauge assignment based on geometry!

# Count neighbor types for C vertices
print(f"  Neighbor types for C vertices (96 golden):")
C_indices = [i for i in range(n_vert) if vtypes[i] == 'C']
neighbor_A_counts = []
neighbor_B_counts = []
neighbor_C_counts = []
for i in C_indices:
    nA = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'A')
    nB = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'B')
    nC = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'C')
    neighbor_A_counts.append(nA)
    neighbor_B_counts.append(nB)
    neighbor_C_counts.append(nC)

nA_vals = Counter(neighbor_A_counts)
nB_vals = Counter(neighbor_B_counts)
nC_vals = Counter(neighbor_C_counts)
print(f"    n_A distribution: {dict(sorted(nA_vals.items()))}")
print(f"    n_B distribution: {dict(sorted(nB_vals.items()))}")
print(f"    n_C distribution: {dict(sorted(nC_vals.items()))}")

# For A vertices (2 poles):
print(f"\n  Neighbor types for A vertices (2 poles):")
for i in A_indices:
    nA = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'A')
    nB = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'B')
    nC = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'C')
    print(f"    Vertex {i}: nA={nA}, nB={nB}, nC={nC} (total={nA+nB+nC})")

# For B vertices:
print(f"\n  Neighbor types for B vertices (22):")
B_indices = [i for i in range(n_vert) if vtypes[i] == 'B']
B_neighbor_type_dist = Counter()
for i in B_indices:
    nA = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'A')
    nB = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'B')
    nC = sum(1 for j in range(n_vert) if A_adj[i,j] and vtypes[j] == 'C')
    pattern = (nA, nB, nC)
    B_neighbor_type_dist[pattern] += 1

for pattern, cnt in sorted(B_neighbor_type_dist.items()):
    print(f"    (nA={pattern[0]}, nB={pattern[1]}, nC={pattern[2]}): {cnt} B vertices")

# =====================================================================
# SECTION 10: Gauge assignment from vertex types
# =====================================================================
print("\n--- SECTION 10: Gauge assignment from vertex types ---")

# Assign gauge type to each edge based on endpoint types:
# Approach: the gauge sector depends on the relationship of the edge
# to the Hopf fibration structure.
#
# Edge (A,*) -> U(1) sector (connected to pole)
# Edge (B,B) or (B,C) -> could be SU(2) or SU(3) sector
# Edge (C,C) -> likely SU(3) sector
#
# Count edges per type:
print(f"  Edge type counts:")
for et in sorted(edge_type_counts.keys()):
    pct = 100 * edge_type_counts[et] / n_edges
    print(f"    {str(et):12s}: {edge_type_counts[et]:4d} ({pct:.1f}%)")

# The 1+3+8 = 12 decomposition suggests per C vertex:
# 1 edge to U(1), 3 edges to SU(2), 8 edges to SU(3)
# Total from 96 C vertices: 96*1=96 U(1) edge-ends, 96*3=288 SU(2) edge-ends, 96*8=768 SU(3) edge-ends
# BUT each edge has 2 ends, so: 96/2=48 U(1) edges? No, it's per-vertex.
# The edges touching A: 2*12=24 edge-ends -> 24 edges (each A end counted once)
# If exactly 1 of 12 neighbors of each C vertex is A, then nA=96 edges of type (A,C) or (A,B)?
# But A has only 2 vertices with 12 edges each = 24 edges total.
# 96 C vertices with nA=0 or 1 each.

# Let's try a different approach: physical gauge type from DISTANCE to poles
# Idea: The Hopf fibration naturally defines a "latitude" on S^2 base.
# Vertices near the poles (small |Im(q)|/|Re(q)|) are in U(1) sector.

# =====================================================================
# SECTION 11: COMBINED Hopf x Vertex-type contingency
# =====================================================================
print("\n--- SECTION 11: Combined Hopf x Vertex-type (full table) ---")

# Build full contingency: for each edge, record (Hopf_class, vertex_edge_type)
# This is the KEY calculation the user requested.

# Map edge types to gauge sectors
def edge_to_gauge_sector(et):
    """Map vertex-type edge to gauge sector."""
    if 'A' in et:
        return 'U1'   # Connected to pole
    elif et == ('B', 'B'):
        return 'SU2'  # Between intermediate vertices
    elif et == ('B', 'C'):
        return 'SU2'  # Intermediate-golden
    elif et == ('C', 'C'):
        return 'SU3'  # Between golden vertices
    return 'unknown'

# Assign gauge sector to each edge
edge_gauge = []
for i_e in range(len(edge_list)):
    edge_gauge.append(edge_to_gauge_sector(edge_vtype[i_e]))

gauge_counts = Counter(edge_gauge)
print(f"  Gauge sector assignment (from vertex types):")
for sector, cnt in sorted(gauge_counts.items()):
    pct = 100 * cnt / n_edges
    print(f"    {sector:7s}: {cnt:4d} edges ({pct:.1f}%)")

# Check expected: U(1) ~ 1/12 of edges, SU(2) ~ 3/12, SU(3) ~ 8/12
print(f"  Expected (1:3:8): U1={n_edges/12:.0f}, SU2={3*n_edges/12:.0f}, SU3={8*n_edges/12:.0f}")

# NOW: the full contingency table Hopf x Gauge
print(f"\n  CONTINGENCY TABLE: Hopf class x Gauge sector")
print(f"  {'':>20s} | {'U1':>8s} {'SU2':>8s} {'SU3':>8s} | {'Total':>7s}")
print(f"  {'-'*60}")

contingency_hg = defaultdict(int)
for i_e in range(len(edge_list)):
    hc = hopf_class[i_e]
    gs = edge_gauge[i_e]
    contingency_hg[(hc, gs)] += 1

for hc in range(3):
    row = []
    for gs in ['U1', 'SU2', 'SU3']:
        row.append(contingency_hg.get((hc, gs), 0))
    total = sum(row)
    print(f"  {hopf_class_names[hc]:>20s} | {row[0]:8d} {row[1]:8d} {row[2]:8d} | {total:7d}")

# Column totals
col_tots = [sum(contingency_hg.get((hc, gs), 0) for hc in range(3)) for gs in ['U1', 'SU2', 'SU3']]
print(f"  {'Total':>20s} | {col_tots[0]:8d} {col_tots[1]:8d} {col_tots[2]:8d} | {sum(col_tots):7d}")

# Chi-squared independence test
print(f"\n  Independence test:")
total_edges = sum(col_tots)
row_tots = [sum(contingency_hg.get((hc, gs), 0) for gs in ['U1', 'SU2', 'SU3']) for hc in range(3)]
chi2 = 0
print(f"  {'':>20s} | {'E(U1)':>8s} {'E(SU2)':>8s} {'E(SU3)':>8s}")
print(f"  {'-'*50}")
for hc in range(3):
    expected = []
    for ci, gs in enumerate(['U1', 'SU2', 'SU3']):
        e_ij = row_tots[hc] * col_tots[ci] / total_edges
        o_ij = contingency_hg.get((hc, gs), 0)
        if e_ij > 0:
            chi2 += (o_ij - e_ij)**2 / e_ij
        expected.append(e_ij)
    print(f"  {hopf_class_names[hc]:>20s} | {expected[0]:8.1f} {expected[1]:8.1f} {expected[2]:8.1f}")

print(f"\n  Chi-squared = {chi2:.4f}")
print(f"  Degrees of freedom = (3-1)*(3-1) = 4")
# For df=4, chi2 > 9.49 means p < 0.05
print(f"  {'DEPENDENT (anisotropic)' if chi2 > 9.49 else 'INDEPENDENT (isotropic)'} at 95% level")
if chi2 > 13.28:
    print(f"  DEPENDENT at 99% level (chi2 > 13.28)")

# =====================================================================
# SECTION 12: Mean fiber fraction per gauge sector
# =====================================================================
print("\n--- SECTION 12: Mean fiber fraction per gauge sector ---")

# THIS IS THE KEY RESULT
fiber_by_gauge = defaultdict(list)
for i_e in range(len(edge_list)):
    gs = edge_gauge[i_e]
    fiber_by_gauge[gs].append(fiber_fractions[i_e])

print(f"  Mean fiber fraction <f> per gauge sector:")
print(f"  {'Sector':>8s} {'<f>':>10s} {'std':>10s} {'n':>6s}")
print(f"  {'-'*40}")
sector_means = {}
for gs in ['U1', 'SU2', 'SU3']:
    vals = fiber_by_gauge[gs]
    if vals:
        m = np.mean(vals)
        s = np.std(vals)
        sector_means[gs] = m
        print(f"  {gs:>8s} {m:10.6f} {s:10.6f} {len(vals):6d}")

# Check if the three means are different
if len(sector_means) == 3:
    f_u1 = sector_means['U1']
    f_su2 = sector_means['SU2']
    f_su3 = sector_means['SU3']
    print(f"\n  ANISOTROPY CHECK:")
    print(f"    <f>_U1  = {f_u1:.6f}")
    print(f"    <f>_SU2 = {f_su2:.6f}")
    print(f"    <f>_SU3 = {f_su3:.6f}")

    if abs(f_u1 - f_su2) > 0.001 or abs(f_su2 - f_su3) > 0.001:
        diff_12 = f_u1 - f_su2
        diff_13 = f_u1 - f_su3
        diff_23 = f_su2 - f_su3
        print(f"\n    <f>_U1  - <f>_SU2 = {diff_12:.6f}")
        print(f"    <f>_U1  - <f>_SU3 = {diff_13:.6f}")
        print(f"    <f>_SU2 - <f>_SU3 = {diff_23:.6f}")

        # Compare with known corrections
        print(f"\n  COMPARISON WITH KNOWN CORRECTIONS:")
        print(f"    Known: delta_d(lepton) = sin2tW = {SIN2_TW_FW:.6f}")
        print(f"    Known: delta_d(up)     = alpha_s*2/pi = {ALPHA_S_FW*2/PI:.6f}")
        print(f"    Known: delta_d(down)   = -1/a1 = {-1/a1:.6f}")
        print(f"    Known: delta_k(lepton) = -1/phi^4 = {-1/PHI**4:.6f}")
        print(f"    Known: delta_k(up)     = +alpha_s = {ALPHA_S_FW:.6f}")
        print(f"    Known: delta_k(down)   = -alpha_s = {-ALPHA_S_FW:.6f}")

# =====================================================================
# SECTION 13: Alternative gauge assignment (1+3+8 from neighbor count)
# =====================================================================
print("\n\n--- SECTION 13: Alternative gauge assignment ---")

# Instead of using vertex TYPE, use the actual NEIGHBOR PATTERN:
# For each C vertex, which of its 12 neighbors is the "1" (U(1) direction)?
# It should be the neighbor closest to a pole (highest |Re(q)|).
#
# For a generic C vertex g, its 12 neighbors h_i are at g*h where h runs
# over the 12 neighbors of identity. The "U(1) direction" is the one
# closest to the Hopf fiber axis.

print(f"  For each C vertex: assign gauge type to each of 12 edges")
print(f"  based on neighbor's Re(q) ranking (proxy for U(1) alignment)")

# For each vertex, rank its 12 neighbors by |Re(q)|
# The one with highest |Re(q)| is most aligned with the pole -> U(1)
# Next 3 -> SU(2), remaining 8 -> SU(3)

gauge_labels_global = {}  # (i,j) -> 'U1', 'SU2', 'SU3'

for i in range(n_vert):
    neighbors = [(j, abs(V[j][0])) for j in range(n_vert) if A_adj[i, j]]
    # Sort by |Re(q)| descending (closest to pole = U(1))
    neighbors_sorted = sorted(neighbors, key=lambda x: -x[1])

    for rank, (j, _) in enumerate(neighbors_sorted):
        edge_key = (min(i, j), max(i, j))
        if rank == 0:
            label = 'U1'
        elif rank < 4:
            label = 'SU2'
        else:
            label = 'SU3'

        # Each edge gets TWO labels (one from each endpoint)
        # Store as list to later check consistency
        if edge_key not in gauge_labels_global:
            gauge_labels_global[edge_key] = []
        gauge_labels_global[edge_key].append(label)

# Check consistency: does each edge get the same label from both endpoints?
consistent = 0
inconsistent = 0
edge_gauge_alt = {}
for edge_key, labels in gauge_labels_global.items():
    if len(labels) == 2 and labels[0] == labels[1]:
        consistent += 1
        edge_gauge_alt[edge_key] = labels[0]
    else:
        inconsistent += 1
        # Use "lower" gauge group (more restrictive assignment)
        # Priority: U1 > SU2 > SU3
        priority = {'U1': 0, 'SU2': 1, 'SU3': 2}
        edge_gauge_alt[edge_key] = min(labels, key=lambda x: priority[x])

print(f"  Consistent labels: {consistent}/{n_edges} ({100*consistent/n_edges:.1f}%)")
print(f"  Inconsistent labels: {inconsistent}/{n_edges}")

# Count gauge sectors with alternative assignment
gauge_alt_counts = Counter(edge_gauge_alt.values())
print(f"\n  Alternative gauge sector counts:")
for gs in ['U1', 'SU2', 'SU3']:
    cnt = gauge_alt_counts.get(gs, 0)
    print(f"    {gs:7s}: {cnt:4d} edges ({100*cnt/n_edges:.1f}%)")

# Contingency table with alternative gauge assignment
print(f"\n  CONTINGENCY TABLE (alternative): Hopf class x Gauge sector")
print(f"  {'':>20s} | {'U1':>8s} {'SU2':>8s} {'SU3':>8s} | {'Total':>7s}")
print(f"  {'-'*60}")

contingency_alt = defaultdict(int)
fiber_by_gauge_alt = defaultdict(list)

for i_e, (i, j) in enumerate(edge_list):
    edge_key = (min(i, j), max(i, j))
    hc = hopf_class[i_e]
    gs = edge_gauge_alt.get(edge_key, 'SU3')
    contingency_alt[(hc, gs)] += 1
    fiber_by_gauge_alt[gs].append(fiber_fractions[i_e])

for hc in range(3):
    row = []
    for gs in ['U1', 'SU2', 'SU3']:
        row.append(contingency_alt.get((hc, gs), 0))
    total = sum(row)
    print(f"  {hopf_class_names[hc]:>20s} | {row[0]:8d} {row[1]:8d} {row[2]:8d} | {total:7d}")

col_tots_alt = [sum(contingency_alt.get((hc, gs), 0) for hc in range(3)) for gs in ['U1', 'SU2', 'SU3']]
print(f"  {'Total':>20s} | {col_tots_alt[0]:8d} {col_tots_alt[1]:8d} {col_tots_alt[2]:8d} | {sum(col_tots_alt):7d}")

# Chi-squared for alternative
chi2_alt = 0
row_tots_alt = [sum(contingency_alt.get((hc, gs), 0) for gs in ['U1', 'SU2', 'SU3']) for hc in range(3)]
total_alt = sum(col_tots_alt)
for hc in range(3):
    for ci, gs in enumerate(['U1', 'SU2', 'SU3']):
        e_ij = row_tots_alt[hc] * col_tots_alt[ci] / total_alt if total_alt > 0 else 0
        o_ij = contingency_alt.get((hc, gs), 0)
        if e_ij > 0:
            chi2_alt += (o_ij - e_ij)**2 / e_ij

print(f"\n  Chi-squared (alt) = {chi2_alt:.4f}")
print(f"  {'DEPENDENT' if chi2_alt > 9.49 else 'INDEPENDENT'} at 95% level")

# Mean fiber fraction per alternative gauge sector
print(f"\n  Mean fiber fraction <f> per gauge sector (alternative):")
print(f"  {'Sector':>8s} {'<f>':>10s} {'std':>10s} {'n':>6s}")
print(f"  {'-'*40}")
sector_means_alt = {}
for gs in ['U1', 'SU2', 'SU3']:
    vals = fiber_by_gauge_alt[gs]
    if vals:
        m = np.mean(vals)
        s = np.std(vals)
        sector_means_alt[gs] = m
        print(f"  {gs:>8s} {m:10.6f} {s:10.6f} {len(vals):6d}")

# =====================================================================
# SECTION 14: Hopf fiber fraction vs geodesic distance to poles
# =====================================================================
print("\n\n--- SECTION 14: Fiber fraction vs distance to poles ---")

# For each edge, compute the geodesic distance of its midpoint to the
# nearest pole. This is a continuous version of the A/B/C classification.

# Distance of each vertex to the north pole (1,0,0,0)
north_pole = V[A_indices[0]]
vertex_dist_to_pole = np.array([
    np.arccos(np.clip(np.dot(V[i], north_pole), -1, 1))
    for i in range(n_vert)
])

# For each edge: midpoint distance to nearest pole
edge_midpoint_dist = []
for i, j in edge_list:
    # Midpoint on S^3 (normalize the average)
    mid = (V[i] + V[j])
    mid = mid / np.linalg.norm(mid)
    d_north = np.arccos(np.clip(np.dot(mid, north_pole), -1, 1))
    d_south = PI - d_north  # south pole is antipodal
    edge_midpoint_dist.append(min(d_north, d_south))

edge_midpoint_dist = np.array(edge_midpoint_dist)

# Correlation between fiber fraction and distance to pole
from numpy import corrcoef
corr = np.corrcoef(fiber_fractions, edge_midpoint_dist)[0, 1]
print(f"  Correlation(fiber_fraction, dist_to_pole) = {corr:.6f}")

# Bin edges by distance to pole and compute mean fiber fraction
n_bins = 8
dist_bins = np.linspace(0, PI/2, n_bins + 1)
print(f"\n  Fiber fraction vs distance to nearest pole:")
print(f"  {'Dist range':>20s} {'<f>':>10s} {'n':>6s}")
print(f"  {'-'*40}")
for b in range(n_bins):
    mask = (edge_midpoint_dist >= dist_bins[b]) & (edge_midpoint_dist < dist_bins[b+1])
    n_bin = np.sum(mask)
    if n_bin > 0:
        mf = fiber_fractions[mask].mean()
        print(f"  [{dist_bins[b]:.3f}, {dist_bins[b+1]:.3f}){' ':>5s} {mf:10.6f} {n_bin:6d}")

# =====================================================================
# SECTION 15: Fiber fraction by Hopf latitude
# =====================================================================
print("\n--- SECTION 15: Hopf map and base-point analysis ---")

def hopf_map(q):
    """Hopf map S^3 -> S^2 using i-fibration."""
    qi = qmul(q, np.array([0, 1, 0, 0]))
    qinv_q = np.array([q[0], -q[1], -q[2], -q[3]])
    result = qmul(qi, qinv_q)
    return result[1:4]

# Hopf image of each vertex
hopf_pts = np.array([hopf_map(V[i]) for i in range(n_vert)])

# "Latitude" on S^2 base = z-component of Hopf image
hopf_z = hopf_pts[:, 2]  # z-component (using third coord)
hopf_lat = np.arcsin(np.clip(hopf_z, -1, 1))  # latitude in [-pi/2, pi/2]

print(f"  Hopf latitude distribution:")
print(f"    Range: [{hopf_lat.min():.4f}, {hopf_lat.max():.4f}]")
print(f"    Mean: {hopf_lat.mean():.4f}")

# Edge midpoint Hopf latitude
edge_hopf_lat = []
for i, j in edge_list:
    mid = (V[i] + V[j])
    mid = mid / np.linalg.norm(mid)
    hpt = hopf_map(mid)
    edge_hopf_lat.append(np.arcsin(np.clip(hpt[2], -1, 1)))

edge_hopf_lat = np.array(edge_hopf_lat)

# Correlation between fiber fraction and Hopf latitude
corr_lat = np.corrcoef(fiber_fractions, np.abs(edge_hopf_lat))[0, 1]
print(f"\n  Correlation(fiber_fraction, |Hopf latitude|) = {corr_lat:.6f}")

# =====================================================================
# SECTION 16: Summary and assessment
# =====================================================================
print("\n" + "=" * 72)
print("SUMMARY: EXP-313c")
print("=" * 72)

print(f"""
VERTEX CLASSIFICATION (corrected):
  A = 2 poles: (+1,0,0,0) and (-1,0,0,0)
  B = 22 vertices: 6 other axis + 16 half-integer
  C = 96 golden ratio vertices
  VERIFIED: matches paper specification

HOPF x GAUGE CROSS-TABULATION:

  Method 1 (A5 projector weights):
    Gauge weights are ISOTROPIC: w_R = d_R/12 for ALL edges.
    No anisotropy between Hopf classes. CONFIRMED negative from exp313.

  Method 2 (vertex-type classification):""")

if contingency_hg:
    print(f"    Vertex-type edges x Hopf: chi2 = {chi2:.4f}")
    print(f"    {'ANISOTROPIC' if chi2 > 9.49 else 'ISOTROPIC'} at 95% level")

print(f"""
  Method 3 (rank-based 1+3+8 per vertex):
    Chi-squared = {chi2_alt:.4f}
    {'ANISOTROPIC' if chi2_alt > 9.49 else 'ISOTROPIC'} at 95% level""")

print(f"""
MEAN FIBER FRACTION PER GAUGE SECTOR:""")

if sector_means:
    for gs in ['U1', 'SU2', 'SU3']:
        if gs in sector_means:
            print(f"  Method 2: <f>_{gs} = {sector_means[gs]:.6f}")

print()
if sector_means_alt:
    for gs in ['U1', 'SU2', 'SU3']:
        if gs in sector_means_alt:
            print(f"  Method 3: <f>_{gs} = {sector_means_alt[gs]:.6f}")

print(f"""
CORRELATIONS:
  Fiber fraction vs dist to pole: r = {corr:.6f}
  Fiber fraction vs |Hopf latitude|: r = {corr_lat:.6f}
""")

print("DONE.")

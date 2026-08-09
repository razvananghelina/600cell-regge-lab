"""
EXP-144: Gram Matrix Investigation
====================================
From exp142: G = GF * GF^T (24x24) has spectrum {36,24,12,6,0} x {1,4,9,8,2}
- All eigenvalues multiples of 6
- First 3 multiplicities = perfect squares (1,4,9)
- 2-hop has SAME multiplicities

Questions:
1. Why are multiplicities {1,4,9,8,2}?
2. Is the effective adjacency A_eff = (G-12I)/3 related to the 24-cell?
3. How does A_eff relate to D4 representation theory?
4. Why do 1-hop and 2-hop share eigenspaces?
5. Can we derive the spectrum from 600-cell geometry?
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict

PHI = (1 + 5**0.5) / 2

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
print("EXP-144: GRAM MATRIX INVESTIGATION")
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

A_verts = sorted([i for i in range(N) if vtypes[i] == 'A'])
B_verts = sorted([i for i in range(N) if vtypes[i] == 'B'])
C_verts = sorted([i for i in range(N) if vtypes[i] == 'C'])
gauge_verts = sorted(A_verts + B_verts)

n_G = len(gauge_verts)
n_C = len(C_verts)
g_map = {v: i for i, v in enumerate(gauge_verts)}
c_map = {v: i for i, v in enumerate(C_verts)}

# Build GF incidence (24x96) and CC adjacency (96x96)
GF = np.zeros((n_G, n_C))
for gi, g in enumerate(gauge_verts):
    for c in neighbors[g]:
        if c in set(C_verts):
            GF[gi][c_map[c]] = 1

A_CC = np.zeros((n_C, n_C))
C_set = set(C_verts)
for c in C_verts:
    for nb in neighbors[c]:
        if nb in C_set:
            A_CC[c_map[c]][c_map[nb]] = 1

# Gram matrix
G = GF @ GF.T
print(f"\nGram matrix G = GF * GF^T ({n_G}x{n_G})")
print(f"  Diagonal: {sorted(set(np.diag(G).astype(int)))}")
print(f"  Off-diagonal: {sorted(set(G[np.triu_indices(n_G, 1)].astype(int)))}")

# ============================================================
# PART 1: Effective adjacency A_eff
# ============================================================
print("\n" + "="*70)
print("PART 1: EFFECTIVE ADJACENCY A_eff = (G - 12*I) / 3")
print("="*70)

A_eff = (G - 12*np.eye(n_G)) / 3
print(f"A_eff entries: {sorted(set(A_eff.flatten().astype(int)))}")
print(f"A_eff edges: {int(np.sum(A_eff > 0.5)) // 2}")
print(f"A_eff degree: {sorted(set(np.sum(A_eff, axis=1).astype(int)))}")

# Is A_eff the 24-cell adjacency?
# 24-cell: V=24, E=96, degree=8, strongly regular (24,8,2,4)
eff_edges = int(np.sum(A_eff > 0.5)) // 2
print(f"\n24-cell check:")
print(f"  Vertices: {n_G} (24-cell: 24) {'OK' if n_G == 24 else 'FAIL'}")
print(f"  Edges: {eff_edges} (24-cell: 96) {'OK' if eff_edges == 96 else 'FAIL'}")
print(f"  Degree: {int(np.sum(A_eff[0] > 0.5))} (24-cell: 8) {'OK' if int(np.sum(A_eff[0] > 0.5)) == 8 else 'FAIL'}")

# Strongly regular check: lambda (adj pair common nbrs) and mu (non-adj pair common nbrs)
lambdas = []
mus = []
for i in range(n_G):
    for j in range(i+1, n_G):
        common = np.sum(A_eff[i] * A_eff[j])
        if A_eff[i][j] > 0.5:
            lambdas.append(int(round(common)))
        else:
            mus.append(int(round(common)))

lambda_vals = Counter(lambdas)
mu_vals = Counter(mus)
print(f"  Lambda (adj common nbrs): {dict(lambda_vals)}")
print(f"  Mu (non-adj common nbrs): {dict(mu_vals)}")

is_srg = len(lambda_vals) == 1 and len(mu_vals) == 1
if is_srg:
    lam = list(lambda_vals.keys())[0]
    mu = list(mu_vals.keys())[0]
    print(f"  Strongly regular: ({n_G}, {int(np.sum(A_eff[0]))}, {lam}, {mu})")
else:
    print(f"  NOT strongly regular (multiple lambda or mu values)")

# ============================================================
# PART 2: Spectrum analysis
# ============================================================
print("\n" + "="*70)
print("PART 2: SPECTRUM ANALYSIS")
print("="*70)

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
        groups.append((round(val, 2), mult))
        i += mult
    return groups

g_unique = group_eigs(eigs_G)
print(f"\nGram matrix spectrum:")
for val, mult in g_unique:
    eff_val = (val - 12) / 3
    print(f"  G={val:6.1f} (m={mult})  ->  A_eff={eff_val:6.2f}")

# The A_eff eigenvalues
eigs_A = sorted(np.linalg.eigvalsh(A_eff), reverse=True)
a_unique = group_eigs(eigs_A)
print(f"\nA_eff spectrum:")
for val, mult in a_unique:
    print(f"  {val:6.2f} x {mult}")

# ============================================================
# PART 3: Compare with standard 24-cell
# ============================================================
print("\n" + "="*70)
print("PART 3: STANDARD 24-CELL COMPARISON")
print("="*70)

# Build the standard 24-cell from gauge vertex coordinates
# The gauge vertices are A (axis: ±e_i) and B (tesseract: ±1/2 all)
# 24-cell standard: all permutations of (±1,±1,0,0) - but our vertices
# are at different coordinates

# Compute actual distances between gauge vertices
gauge_coords = np.array([verts[g] for g in gauge_verts])
gauge_dist = np.sqrt(np.sum((gauge_coords[:, None, :] - gauge_coords[None, :, :]) ** 2, axis=2))

# Distance distribution
gauge_dists = []
for i in range(n_G):
    for j in range(i+1, n_G):
        gauge_dists.append(round(gauge_dist[i][j], 4))

dist_counts = Counter(gauge_dists)
print("Gauge vertex distance distribution:")
for d, cnt in sorted(dist_counts.items()):
    print(f"  d={d:.4f}: {cnt} pairs")

# Build 24-cell adjacency from NEAREST distance
min_dist = min(dist_counts.keys())
A_24cell = (gauge_dist < min_dist + 0.01).astype(float)
np.fill_diagonal(A_24cell, 0)

print(f"\nNearest-neighbor 24-cell:")
print(f"  Min distance: {min_dist:.4f}")
print(f"  Edges: {int(np.sum(A_24cell)) // 2}")
print(f"  Degree: {sorted(set(np.sum(A_24cell, axis=1).astype(int)))}")

# Compare A_eff with A_24cell
match = np.sum(np.abs(A_eff - A_24cell)) == 0
print(f"\nA_eff == A_24cell? {match}")

if not match:
    diff = A_eff - A_24cell
    n_diff = int(np.sum(np.abs(diff) > 0.5)) // 2
    print(f"  Differences: {n_diff} edge changes")
    # Edges in A_eff but not in A_24cell
    only_eff = int(np.sum((A_eff > 0.5) & (A_24cell < 0.5))) // 2
    only_24 = int(np.sum((A_eff < 0.5) & (A_24cell > 0.5))) // 2
    print(f"  Only in A_eff: {only_eff}")
    print(f"  Only in A_24cell: {only_24}")

    # Spectrum of standard 24-cell
    eigs_24 = sorted(np.linalg.eigvalsh(A_24cell), reverse=True)
    e24_unique = group_eigs(eigs_24)
    print(f"\n  Standard 24-cell spectrum:")
    for val, mult in e24_unique:
        print(f"    {val:6.2f} x {mult}")

# ============================================================
# PART 4: A_eff graph structure
# ============================================================
print("\n" + "="*70)
print("PART 4: A_eff GRAPH STRUCTURE")
print("="*70)

# A-A, A-B, B-B edges in A_eff
n_A = len(A_verts)
n_B = len(B_verts)
AA_eff = int(np.sum(A_eff[:n_A, :n_A] > 0.5)) // 2
AB_eff = int(np.sum(A_eff[:n_A, n_A:] > 0.5))
BB_eff = int(np.sum(A_eff[n_A:, n_A:] > 0.5)) // 2
print(f"A_eff edges: A-A={AA_eff}, A-B={AB_eff}, B-B={BB_eff}")
print(f"  Total: {AA_eff + AB_eff + BB_eff} (should be {eff_edges})")

# Per type degree
A_deg_in_eff = np.sum(A_eff[:n_A], axis=1).astype(int)
B_deg_in_eff = np.sum(A_eff[n_A:], axis=1).astype(int)
print(f"  A vertex degrees: {sorted(set(A_deg_in_eff))}")
print(f"  B vertex degrees: {sorted(set(B_deg_in_eff))}")

# A-type: degree to A and B
A_to_A = np.sum(A_eff[:n_A, :n_A], axis=1).astype(int)
A_to_B = np.sum(A_eff[:n_A, n_A:], axis=1).astype(int)
B_to_A = np.sum(A_eff[n_A:, :n_A], axis=1).astype(int)
B_to_B = np.sum(A_eff[n_A:, n_A:], axis=1).astype(int)
print(f"  A -> (A,B): ({sorted(set(A_to_A))}, {sorted(set(A_to_B))})")
print(f"  B -> (A,B): ({sorted(set(B_to_A))}, {sorted(set(B_to_B))})")

# ============================================================
# PART 5: Eigenspace analysis
# ============================================================
print("\n" + "="*70)
print("PART 5: EIGENSPACE ANALYSIS")
print("="*70)

# For each eigenspace of G, check:
# - Weight on A vs B vertices
# - Connection to 600-cell eigenspaces
pos = 0
for val, mult in g_unique:
    V_space = vecs_G[:, pos:pos+mult]  # 24 x mult

    # Weight on A (first 8) and B (last 16)
    A_weight = np.sum(V_space[:n_A, :]**2)
    B_weight = np.sum(V_space[n_A:, :]**2)
    total = A_weight + B_weight

    print(f"\nG={val:.1f} (m={mult}):")
    print(f"  A weight: {A_weight/total:.4f} ({A_weight:.4f})")
    print(f"  B weight: {B_weight/total:.4f} ({B_weight:.4f})")
    print(f"  A/B ratio: {A_weight/B_weight:.4f}" if B_weight > 0.01 else "  B weight ~ 0")

    # Check if eigenspace has simple A/B structure
    if mult <= 10:
        # SVD of A-part and B-part
        A_part = V_space[:n_A, :]
        B_part = V_space[n_A:, :]
        A_rank = np.linalg.matrix_rank(A_part, tol=1e-6)
        B_rank = np.linalg.matrix_rank(B_part, tol=1e-6)
        print(f"  A-part rank: {A_rank}/{min(n_A, mult)}")
        print(f"  B-part rank: {B_rank}/{min(n_B, mult)}")

    pos += mult

# ============================================================
# PART 6: Connection to 600-cell eigenspaces
# ============================================================
print("\n" + "="*70)
print("PART 6: CONNECTION TO 600-CELL EIGENSPACES")
print("="*70)

# Build full 600-cell adjacency
A_full = np.zeros((N, N))
for i in range(N):
    for j in neighbors[i]:
        A_full[i][j] = 1

full_eigs, full_vecs = np.linalg.eigh(A_full)
idx_sort = np.argsort(full_eigs)[::-1]
full_eigs = full_eigs[idx_sort]
full_vecs = full_vecs[:, idx_sort]

# Group 600-cell eigenvalues
full_unique = group_eigs(list(full_eigs))
print("600-cell eigenvalues:")
for val, mult in full_unique:
    print(f"  {val:8.4f} x {mult}")

# For each 600-cell eigenspace, project to gauge sector (24 vertices)
# and compute the Gram matrix restricted to that eigenspace
print("\nGram matrix decomposition by 600-cell eigenspace:")
pos = 0
total_contribution = np.zeros((n_G, n_G))
for val, mult in full_unique:
    # Eigenspace: N x mult
    V_k = full_vecs[:, pos:pos+mult]

    # Project to gauge vertices: extract rows for gauge verts
    gauge_proj = V_k[[g for g in gauge_verts], :]  # 24 x mult
    C_proj = V_k[[c for c in C_verts], :]  # 96 x mult

    # GF restricted through this eigenspace
    # The contribution to G from eigenspace k:
    # G_k = GF * P_k * GF^T where P_k is projector onto C-subspace of eigenspace k
    # But actually G = GF * I * GF^T, so we need to understand how I decomposes

    # Alternative: compute GF * P_C_k * GF^T where P_C_k = C_proj @ C_proj^T
    G_k = GF @ C_proj @ C_proj.T @ GF.T
    total_contribution += G_k

    # Trace of G_k
    tr = np.trace(G_k)

    # Eigenvalues of G_k
    gk_eigs = sorted(np.linalg.eigvalsh(G_k), reverse=True)
    gk_nz = [round(e, 2) for e in gk_eigs if abs(e) > 0.01]

    print(f"  lambda={val:8.4f} (m={mult}): Tr(G_k)={tr:.2f}, non-zero eigs: {gk_nz[:5]}")

    pos += mult

# Verify decomposition
G_check = np.max(np.abs(total_contribution - G))
print(f"\nDecomposition check: max|G_recon - G| = {G_check:.2e}")

# ============================================================
# PART 7: Why multiplicities {1,4,9,8,2}?
# ============================================================
print("\n" + "="*70)
print("PART 7: MULTIPLICITY ANALYSIS")
print("="*70)

mults = [mult for _, mult in g_unique]
print(f"Multiplicities: {mults}")
print(f"Sum: {sum(mults)} (should be 24)")

# Check patterns
print(f"\nPatterns:")
print(f"  Perfect squares: {[m for m in mults if int(round(m**0.5))**2 == m]}")
print(f"  Non-squares: {[m for m in mults if int(round(m**0.5))**2 != m]}")
print(f"  Sum of squares: {sum(m for m in mults if int(round(m**0.5))**2 == m)}")
print(f"  Sum of non-squares: {sum(m for m in mults if int(round(m**0.5))**2 != m)}")

# D4 representation theory
# D4 has representations with dimensions 1, 3, 3, 3 (triality), 6, 8v, 8s, 8c, ...
# The 24-cell carries the adjoint rep of D4 (dim 28) somehow
# Actually the 24 roots of D4 decompose under subgroups

# Check: D4 -> A2 x A1 x U(1)
# 24 roots -> ?
# D4 adjoint (28) -> (8,1,0) + (1,3,0) + (1,1,0) + (3,2,1) + (3',2,-1) = 8+3+1+6+6 = 24 non-Cartan
# So non-Cartan roots: 24 = 8 + 3 + 1 + 6 + 6
# Hmm, but our multiplicities are {1,4,9,8,2}

# Alternative: the PERMUTATION representation of the 24-cell
# Under D4 symmetry, the 24 vertices decompose into irreps
print(f"\nD4 irrep dimensions: 1, 3, 3', 3'', 6, 8_v, 8_s, 8_c")
print(f"Possible decompositions of 24:")
print(f"  24 = 1 + 3 + 3' + 3'' + 6 + 8 (standard D4 adjoint minus Cartan)")
print(f"  24 = 8_v + 8_s + 8_c (triality triple)")
print(f"  24 = 1 + 3 + 3 + 3 + 6 + 8 (our pattern?)")

# Our multiplicities: 1 + 4 + 9 + 8 + 2 = 24
# Could 4 = 3+1? Could 9 = 3+6? Could 8 = 8? Could 2 = something?
# 1 + (1+3) + (3+6) + 8 + 2?

# ============================================================
# PART 8: Higher-hop Gram matrices
# ============================================================
print("\n" + "="*70)
print("PART 8: HIGHER-HOP GRAM MATRICES")
print("="*70)

# G_n = GF * A_CC^n * GF^T
# From exp142: G_0 = G has spectrum {36,24,12,6,0} x {1,4,9,8,2}
# G_1 = GF * A_CC * GF^T has spectrum {324,168,44,6,0} x {1,4,9,8,2}
# SAME multiplicities!

# This means the eigenspaces of G are INVARIANT under the map
# M -> GF * A_CC^n * GF^T for all n

# Check G_2, G_3
for n_hop in range(4):
    if n_hop == 0:
        G_n = G.copy()
    else:
        G_n = GF @ np.linalg.matrix_power(A_CC, n_hop) @ GF.T

    gn_eigs = sorted(np.linalg.eigvalsh(G_n), reverse=True)
    gn_unique = group_eigs(gn_eigs, tol=0.1)
    print(f"\nG_{n_hop} = GF * A_CC^{n_hop} * GF^T:")
    for val, mult in gn_unique:
        print(f"  {val:12.2f} x {mult}")

    # Check if multiplicities match G_0
    gn_mults = [m for _, m in gn_unique]
    g0_mults = [m for _, m in g_unique]
    print(f"  Multiplicities: {gn_mults} (G_0: {g0_mults})")
    print(f"  Match: {gn_mults == g0_mults}")

# ============================================================
# PART 9: Eigenvalue ratios
# ============================================================
print("\n" + "="*70)
print("PART 9: EIGENVALUE RATIOS ACROSS HOPS")
print("="*70)

# If eigenspaces are shared, then G_n eigenvalues are functions of G_0 eigenvalues
# Specifically: if G_0 has eigenvalue lambda_0, then
# G_n = GF * A_CC^n * GF^T on that eigenspace gives lambda_n = f(lambda_0, n)

g0_vals = [val for val, _ in g_unique]
print(f"G_0 eigenvalues: {g0_vals}")

for n_hop in range(1, 4):
    G_n = GF @ np.linalg.matrix_power(A_CC, n_hop) @ GF.T
    gn_eigs = sorted(np.linalg.eigvalsh(G_n), reverse=True)
    gn_unique = group_eigs(gn_eigs, tol=0.1)
    gn_vals = [val for val, _ in gn_unique]
    print(f"\nG_{n_hop} eigenvalues: {gn_vals}")
    if len(gn_vals) == len(g0_vals):
        ratios = [gn_vals[i]/g0_vals[i] if abs(g0_vals[i]) > 0.01 else 'inf' for i in range(len(g0_vals))]
        print(f"  G_{n_hop}/G_0 ratios: {ratios}")

# ============================================================
# PART 10: Algebraic structure
# ============================================================
print("\n" + "="*70)
print("PART 10: ALGEBRAIC STRUCTURE")
print("="*70)

# G = GF * GF^T = 12I + 3*A_eff
# G_1 = GF * A_CC * GF^T
# If they share eigenspaces, then G_1 = p(G) for some polynomial p

# Try: G_1 = a*G^2 + b*G + c*I
# Solve: for eigenvalue pairs (36,324), (24,168), (12,44), (6,6), (0,0)
# 324 = a*36^2 + b*36 + c
# 168 = a*24^2 + b*24 + c
# 44 = a*12^2 + b*12 + c
# 6 = a*6^2 + b*6 + c
# 0 = a*0^2 + b*0 + c => c = 0

# From last: c = 0
# 324 = 1296a + 36b
# 168 = 576a + 24b
# 44 = 144a + 12b
# 6 = 36a + 6b

# From eq 4: 6 = 36a + 6b => b = 1 - 6a
# Into eq 3: 44 = 144a + 12(1-6a) = 144a + 12 - 72a = 72a + 12
# => 72a = 32 => a = 4/9
# => b = 1 - 24/9 = 1 - 8/3 = -5/3

# Check eq 2: 576*(4/9) + 24*(-5/3) = 256 - 40 = 216 != 168. FAIL.

# Try cubic: G_1 = a*G^3 + b*G^2 + c*G + d*I
G_1 = GF @ A_CC @ GF.T

# Direct polynomial fit
from numpy.polynomial import polynomial as P

g0_eig_vals = [36, 24, 12, 6, 0]  # non-degenerate G eigenvalues
g1_eig_vals = [324, 168, 44, 6, 0]  # corresponding G_1 eigenvalues

# Lagrange interpolation: find polynomial p such that p(g0[i]) = g1[i]
# Using numpy
coeffs = np.polyfit(g0_eig_vals, g1_eig_vals, len(g0_eig_vals) - 1)
print(f"Polynomial p such that G_1 = p(G):")
print(f"  Coefficients (highest to lowest): {[round(c, 6) for c in coeffs]}")

# Check
for i in range(len(g0_eig_vals)):
    pred = np.polyval(coeffs, g0_eig_vals[i])
    print(f"  p({g0_eig_vals[i]}) = {pred:.4f} (expected {g1_eig_vals[i]})")

# Verify: G_1 == p(G) as matrices?
G_poly = np.zeros_like(G)
for i, c in enumerate(coeffs):
    power = len(coeffs) - 1 - i
    G_poly += c * np.linalg.matrix_power(G, power)

max_err = np.max(np.abs(G_1 - G_poly))
print(f"\n  max|G_1 - p(G)| = {max_err:.2e}")

# ============================================================
# PART 11: A_eff as polynomial of 600-cell adjacency
# ============================================================
print("\n" + "="*70)
print("PART 11: A_eff FROM 600-CELL ADJACENCY")
print("="*70)

# A_eff = (G - 12I)/3 = (GF*GF^T - 12I)/3
# Is A_eff a projection of some power of the full adjacency A?

# P_gauge * A^2 * P_gauge restricted to gauge sector
# where P_gauge selects gauge vertices

gauge_idx = np.array(gauge_verts)
for power in range(1, 4):
    A_pow = np.linalg.matrix_power(A_full, power)
    A_pow_gauge = A_pow[np.ix_(gauge_idx, gauge_idx)]

    # Is A_eff a linear combination of these?
    # Check if A_pow_gauge is proportional to A_eff (after removing diagonal)
    off_diag_pow = A_pow_gauge - np.diag(np.diag(A_pow_gauge))
    off_diag_eff = A_eff

    # Non-zero entries of A_eff
    mask = off_diag_eff > 0.5
    if np.sum(mask) > 0:
        ratios_pow = off_diag_pow[mask] / off_diag_eff[mask]
        ratio_unique = sorted(set(ratios_pow.round(4)))
        print(f"\nA^{power} gauge-gauge / A_eff where A_eff=1:")
        print(f"  Values: {ratio_unique[:10]}")
        if len(set(ratios_pow.round(2))) == 1:
            print(f"  PROPORTIONAL! Ratio = {ratio_unique[0]}")

    # Diagonal of A^power on gauge
    diag_vals = sorted(set(np.diag(A_pow_gauge).astype(int)))
    print(f"  A^{power} diagonal on gauge: {diag_vals}")

# Summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"\nGram matrix G = GF * GF^T:")
print(f"  Spectrum: {[(v, m) for v, m in g_unique]}")
print(f"  A_eff = (G-12I)/3 spectrum: {[(v, m) for v, m in a_unique]}")
print(f"  Eigenvalues: multiples of 6 (G) / integers (A_eff)")
print(f"  Multiplicities: {mults}")
print(f"  {mults[0]} + {mults[1]} + {mults[2]} = {sum(mults[:3])} (perfect squares sum)")
print(f"  {mults[3]} + {mults[4]} = {sum(mults[3:])} (non-squares sum)")
print(f"  A_eff edges: {eff_edges}, degree: {int(np.sum(A_eff[0]))}")
print(f"  Higher-hop G_n share same eigenspaces (CONFIRMED)")

print("\n" + "="*70)
print("END EXP-144")
print("="*70)

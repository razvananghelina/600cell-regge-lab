"""
exp427d_4d_hessian.py -- 4D Regge Hessian on S^4 (double-cone over 600-cell)
=============================================================================
Compute the full 960x960 Hessian H = d^2 S_R / (dl_a dl_b) of the 4D Regge
action on the S^4 triangulation built in exp427.

S^4 triangulation:
  122 vertices: 120 equatorial (600-cell on S^3) + N pole + S pole
  960 edges:    720 equatorial + 120 N-radial + 120 S-radial
  2640 triangles (hinges): 1200 equatorial + 720 N-rad + 720 S-rad
  1200 four-simplices: ALL congruent (4 edges sqrt(2), 6 edges 1/phi)

Key insight: per-simplex Hessian is UNIVERSAL (computed once, reused 1200x).

4D Regge action: S = sum_triangles A_tri * epsilon_tri
               = sum_tri A_tri * (2*pi - sum_{simp} theta_tri^simp)
               = 2*pi * sum_tri A_tri  -  sum_simp sum_{10 tri} A_tri * theta_tri
               = 2*pi * sum_tri A_tri  -  sum_simp f_simp

So: H = 2*pi * H_area - H_f
"""

import numpy as np
import itertools
import time
from collections import defaultdict

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 7.2973525693e-3

print("=" * 72)
print("EXP-427d: 4D REGGE HESSIAN ON S^4 (DOUBLE-CONE OVER 600-CELL)")
print("=" * 72)

# ============================================================
# PART 1: INFRASTRUCTURE
# ============================================================
print("\n" + "=" * 72)
print("PART 1: INFRASTRUCTURE")
print("=" * 72)
t0 = time.time()

# --- Build 600-cell vertices on unit S^3 in R^4 ---
verts_4d = []
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0, 0, 0, 0]; v[i] = s
        verts_4d.append(v)
for signs in itertools.product([0.5, -0.5], repeat=4):
    verts_4d.append(list(signs))
vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]
even_perms = [p for p in itertools.permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]) % 2 == 0]
for perm in even_perms:
    base = [vals_base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for signs in itertools.product([1, -1], repeat=len(nz)):
        v = list(base)
        for k, idx in enumerate(nz):
            v[idx] *= signs[k]
        verts_4d.append(v)
verts_4d = np.array(verts_4d)

# Remove duplicates
unique = []
for i in range(len(verts_4d)):
    dup = False
    for j in unique:
        if np.linalg.norm(verts_4d[i] - verts_4d[j]) < 1e-10:
            dup = True; break
    if not dup:
        unique.append(i)
verts_4d = verts_4d[unique]
assert len(verts_4d) == 120, f"Got {len(verts_4d)} vertices"

# Build 600-cell edges (chord distance 1/phi)
target_len = 1.0 / PHI
eq_edges = []
adj = {i: set() for i in range(120)}
for i in range(120):
    for j in range(i+1, 120):
        if abs(np.linalg.norm(verts_4d[i] - verts_4d[j]) - target_len) < 1e-8:
            eq_edges.append((i, j))
            adj[i].add(j); adj[j].add(i)
assert len(eq_edges) == 720
eq_edge_set = {(min(i,j), max(i,j)) for i,j in eq_edges}
eq_edge_idx = {}
for idx, (i,j) in enumerate(eq_edges):
    eq_edge_idx[(i,j)] = idx; eq_edge_idx[(j,i)] = idx

# Build 600-cell tetrahedra (4-cliques)
cells = set()
for (i,j) in eq_edges:
    common = adj[i] & adj[j]
    for k in common:
        for l in common:
            if l > k and (min(k,l), max(k,l)) in eq_edge_set:
                cells.add(tuple(sorted([i,j,k,l])))
cells = sorted(cells)
assert len(cells) == 600

# Build 600-cell faces
eq_faces = set()
for cell in cells:
    for tri in itertools.combinations(cell, 3):
        eq_faces.add(tri)
eq_faces = sorted(eq_faces)
assert len(eq_faces) == 1200

# --- S^4 double-cone: embed in R^5, add N and S poles ---
idx_N, idx_S = 120, 121
verts_5d = np.zeros((122, 5))
verts_5d[:120, :4] = verts_4d
verts_5d[idx_N] = [0, 0, 0, 0, 1]
verts_5d[idx_S] = [0, 0, 0, 0, -1]

# Build 1200 four-simplices
four_simplices = []
for tet in cells:
    four_simplices.append(tuple(sorted((idx_N,) + tet)))  # North cone
    four_simplices.append(tuple(sorted((idx_S,) + tet)))  # South cone
assert len(four_simplices) == 1200

# --- Build the FULL edge list (960 edges) ---
all_edges_set = set()
for s4 in four_simplices:
    for e in itertools.combinations(s4, 2):
        all_edges_set.add(e)
all_edges = sorted(all_edges_set)
N_E = len(all_edges)
assert N_E == 960, f"Expected 960 edges, got {N_E}"

edge_to_idx = {}
for idx, e in enumerate(all_edges):
    edge_to_idx[e] = idx
    edge_to_idx[(e[1], e[0])] = idx

# Classify edges
n_eq_e, n_Nrad_e, n_Srad_e = 0, 0, 0
edge_type = []
for e in all_edges:
    if idx_N in e:
        n_Nrad_e += 1; edge_type.append('Nr')
    elif idx_S in e:
        n_Srad_e += 1; edge_type.append('Sr')
    else:
        n_eq_e += 1; edge_type.append('eq')
print(f"  Edges: {N_E} total = {n_eq_e} eq + {n_Nrad_e} Nr + {n_Srad_e} Sr")
assert n_eq_e == 720 and n_Nrad_e == 120 and n_Srad_e == 120

# --- Build FULL triangle list (2640 triangles) ---
all_tris_set = set()
for s4 in four_simplices:
    for t in itertools.combinations(s4, 3):
        all_tris_set.add(t)
all_tris = sorted(all_tris_set)
N_F = len(all_tris)
assert N_F == 2640, f"Expected 2640 triangles, got {N_F}"

tri_to_idx = {}
for idx, t in enumerate(all_tris):
    tri_to_idx[t] = idx

# --- Build simp_edges_list: 10 global edge indices per simplex ---
simp_edges_list = []
for s4 in four_simplices:
    sv = sorted(s4)
    te = []
    for pair in itertools.combinations(sv, 2):
        te.append(edge_to_idx[pair])
    simp_edges_list.append(te)

# --- Build tri_edges_list: 3 global edge indices per triangle ---
tri_edges_list = []
for t in all_tris:
    te = []
    for pair in itertools.combinations(t, 2):
        te.append(edge_to_idx[pair])
    tri_edges_list.append(te)

# --- Triangle-to-simplices map ---
tri_to_simp = [[] for _ in range(N_F)]
for si, s4 in enumerate(four_simplices):
    for t in itertools.combinations(sorted(s4), 3):
        fi = tri_to_idx[t]
        tri_to_simp[fi].append(si)

tri_valences = [len(x) for x in tri_to_simp]
print(f"  Triangle valences: min={min(tri_valences)}, max={max(tri_valences)}")

# --- Build FULL tetrahedra list ---
all_tets_set = set()
for s4 in four_simplices:
    for tt in itertools.combinations(s4, 4):
        all_tets_set.add(tuple(sorted(tt)))
N_T = len(all_tets_set)

chi = 122 - N_E + N_F - N_T + 1200
print(f"  V=122, E={N_E}, F={N_F}, T={N_T}, S4=1200")
print(f"  Euler char = {chi} (expect 2)")
assert chi == 2

# --- Boundary operators ---
print("  Building boundary operators d0, d1...")

d0 = np.zeros((N_E, 122))
for e_idx, (i, j) in enumerate(all_edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = 1.0

d1 = np.zeros((N_F, N_E))
for f_idx, (i, j, k) in enumerate(all_tris):
    e_ij = edge_to_idx[(i, j)]
    e_ik = edge_to_idx[(i, k)]
    e_jk = edge_to_idx[(j, k)]
    d1[f_idx, e_ij] = 1.0
    d1[f_idx, e_jk] = 1.0
    d1[f_idx, e_ik] = -1.0

dd = np.max(np.abs(d1 @ d0))
print(f"  |d1 @ d0| = {dd:.2e} (expect ~0)")
assert dd < 1e-10

# --- Hodge decomposition ---
print("  Computing Hodge decomposition...")
L_coexact = d1.T @ d1
evals_co, evecs_co = np.linalg.eigh(L_coexact)
coexact_mask = evals_co > 0.01
n_coexact = np.sum(coexact_mask)
V_co = evecs_co[:, coexact_mask]
lambda_co = evals_co[coexact_mask]
P_coexact = V_co @ V_co.T

L_exact = d0 @ d0.T
evals_ex, evecs_ex = np.linalg.eigh(L_exact)
exact_mask = evals_ex > 0.01
V_ex = evecs_ex[:, exact_mask]
P_exact = V_ex @ V_ex.T

n_exact = np.sum(exact_mask)
n_harmonic = N_E - n_coexact - n_exact

print(f"  Coexact dim: {n_coexact} (expect {N_E}-121={N_E-121})")
print(f"  Exact dim:   {n_exact} (expect 121)")
print(f"  Harmonic:    {n_harmonic} (expect 0, b_1(S^4)=0)")
print(f"  Tr(P_coexact): {np.trace(P_coexact):.1f}")
print(f"  Tr(P_exact):   {np.trace(P_exact):.1f}")

t1 = time.time()
print(f"  Infrastructure time: {t1-t0:.1f}s")


# ============================================================
# PART 2: PER-SIMPLEX 10x10 HESSIAN (UNIVERSAL)
# ============================================================
print("\n" + "=" * 72)
print("PART 2: PER-SIMPLEX 10x10 HESSIAN (UNIVERSAL)")
print("=" * 72)

# Determine reference edge lengths from first simplex
s0 = four_simplices[0]
sv0 = sorted(s0)
L0_10 = []
local_edge_types = []
EDGE_COMBOS = list(itertools.combinations(range(5), 2))
for pi, pair in enumerate(itertools.combinations(sv0, 2)):
    d = np.linalg.norm(verts_5d[pair[0]] - verts_5d[pair[1]])
    L0_10.append(d)
    if idx_N in pair or idx_S in pair:
        local_edge_types.append('rad')
    else:
        local_edge_types.append('eq')
L0_10 = np.array(L0_10)

radial_mask = np.array([t == 'rad' for t in local_edge_types])
equat_mask = ~radial_mask

print(f"  Reference simplex: {sv0}")
print(f"  Equatorial (idx {np.where(equat_mask)[0]}): "
      f"mean={np.mean(L0_10[equat_mask]):.8f}, std={np.std(L0_10[equat_mask]):.2e}")
print(f"  Radial (idx {np.where(radial_mask)[0]}): "
      f"mean={np.mean(L0_10[radial_mask]):.8f}, std={np.std(L0_10[radial_mask]):.2e}")
print(f"  Expected: eq=1/phi={1/PHI:.8f}, rad=sqrt(2)={np.sqrt(2):.8f}")


def dihedral_angle_4d(p0, p1, p2, p3, p4):
    """Dihedral angle at triangle (p0,p1,p2) in 4-simplex (p0,p1,p2,p3,p4)."""
    u1 = p1 - p0; u2 = p2 - p0; u3 = p3 - p0; u4 = p4 - p0
    e1 = u1 / np.linalg.norm(u1)
    u2o = u2 - np.dot(u2, e1) * e1
    e2 = u2o / np.linalg.norm(u2o)
    u3o = u3 - np.dot(u3, e1)*e1 - np.dot(u3, e2)*e2
    e3 = u3o / np.linalg.norm(u3o)
    u4o = u4 - np.dot(u4, e1)*e1 - np.dot(u4, e2)*e2 - np.dot(u4, e3)*e3
    e4 = u4o / np.linalg.norm(u4o)
    d3 = p3 - p0
    d3c = np.array([np.dot(d3, e3), np.dot(d3, e4)])
    d4 = p4 - p0
    d4c = np.array([np.dot(d4, e3), np.dot(d4, e4)])
    cos_a = np.dot(d3c, d4c) / (np.linalg.norm(d3c) * np.linalg.norm(d4c))
    return np.arccos(np.clip(cos_a, -1.0, 1.0))


def simplex_reconstruct_4d(L10):
    """Reconstruct 5 vertices in R^4 from 10 edge lengths via Cayley-Menger.
    L10 ordered by combinations(range(5), 2): (01,02,03,04,12,13,14,23,24,34)."""
    l01, l02, l03, l04, l12, l13, l14, l23, l24, l34 = L10

    v0 = np.zeros(4)
    v1 = np.array([l01, 0.0, 0.0, 0.0])

    x2 = (l01**2 + l02**2 - l12**2) / (2*l01)
    y2 = np.sqrt(max(l02**2 - x2**2, 0))
    v2 = np.array([x2, y2, 0.0, 0.0])

    x3 = (l01**2 + l03**2 - l13**2) / (2*l01)
    dot23 = (l02**2 + l03**2 - l23**2) / 2
    y3 = (dot23 - x2*x3) / y2 if y2 > 1e-15 else 0.0
    z3 = np.sqrt(max(l03**2 - x3**2 - y3**2, 0))
    v3 = np.array([x3, y3, z3, 0.0])

    x4 = (l01**2 + l04**2 - l14**2) / (2*l01)
    dot24 = (l02**2 + l04**2 - l24**2) / 2
    y4 = (dot24 - x2*x4) / y2 if y2 > 1e-15 else 0.0
    dot34 = (l03**2 + l04**2 - l34**2) / 2
    z4 = (dot34 - x3*x4 - y3*y4) / z3 if z3 > 1e-15 else 0.0
    w4 = np.sqrt(max(l04**2 - x4**2 - y4**2 - z4**2, 0))
    v4 = np.array([x4, y4, z4, w4])

    return np.array([v0, v1, v2, v3, v4])


def triangle_area(p0, p1, p2):
    """Area of triangle in R^d."""
    u = p1 - p0; v = p2 - p0
    return 0.5 * np.sqrt(max(np.dot(u,u)*np.dot(v,v) - np.dot(u,v)**2, 0))


TRI_COMBOS = list(itertools.combinations(range(5), 3))
TRI_OPPOSITES = []
for tri in TRI_COMBOS:
    TRI_OPPOSITES.append([k for k in range(5) if k not in tri])


def compute_f_simplex(L10):
    """Compute f = sum_{10 tri} A_tri * theta_tri for one 4-simplex."""
    verts = simplex_reconstruct_4d(L10)
    f = 0.0
    for tri, opp in zip(TRI_COMBOS, TRI_OPPOSITES):
        A = triangle_area(verts[tri[0]], verts[tri[1]], verts[tri[2]])
        theta = dihedral_angle_4d(verts[tri[0]], verts[tri[1]], verts[tri[2]],
                                  verts[opp[0]], verts[opp[1]])
        f += A * theta
    return f


# Compute per-simplex 10x10 Hessian by central finite differences
print("\n  Computing per-simplex Hessian Hf (10x10) by FD...")

h_fd = 1e-4  # Optimal for 2nd derivative FD: h ~ eps^(1/4) ~ 1e-4
f0 = compute_f_simplex(L0_10)
print(f"  f0 = {f0:.10f}")

Hf = np.zeros((10, 10))

fp = np.zeros(10)
fm = np.zeros(10)
for i in range(10):
    Lp = L0_10.copy(); Lp[i] += h_fd
    Lm = L0_10.copy(); Lm[i] -= h_fd
    fp[i] = compute_f_simplex(Lp)
    fm[i] = compute_f_simplex(Lm)

for i in range(10):
    Hf[i, i] = (fp[i] - 2*f0 + fm[i]) / h_fd**2

for i in range(10):
    for j in range(i+1, 10):
        Lpp = L0_10.copy(); Lpp[i] += h_fd; Lpp[j] += h_fd
        Lpm = L0_10.copy(); Lpm[i] += h_fd; Lpm[j] -= h_fd
        Lmp = L0_10.copy(); Lmp[i] -= h_fd; Lmp[j] += h_fd
        Lmm = L0_10.copy(); Lmm[i] -= h_fd; Lmm[j] -= h_fd
        Hf[i, j] = (compute_f_simplex(Lpp) - compute_f_simplex(Lpm)
                     - compute_f_simplex(Lmp) + compute_f_simplex(Lmm)) / (4*h_fd**2)
        Hf[j, i] = Hf[i, j]

print(f"  Hf symmetry: |Hf - Hf^T| = {np.max(np.abs(Hf - Hf.T)):.2e}")
print(f"  Hf range: [{Hf.min():.6f}, {Hf.max():.6f}]")
print(f"  Hf norm: {np.linalg.norm(Hf):.6f}")

# Check S4 symmetry: classify the 10x10 entries
eq_local = [i for i in range(10) if not radial_mask[i]]
rad_local = [i for i in range(10) if radial_mask[i]]
print(f"  Local equatorial edges: {eq_local}")
print(f"  Local radial edges:     {rad_local}")

def classify_pair(i, j):
    ti = 'r' if radial_mask[i] else 'e'
    tj = 'r' if radial_mask[j] else 'e'
    vi_set = set(EDGE_COMBOS[i])
    vj_set = set(EDGE_COMBOS[j])
    shared = len(vi_set & vj_set)
    if i == j:
        return f"{ti}{tj}_same"
    elif shared == 1:
        return f"{ti}{tj}_adj"
    else:
        return f"{ti}{tj}_opp"

classes = defaultdict(list)
for i in range(10):
    for j in range(i, 10):
        cls = classify_pair(i, j)
        classes[cls].append(Hf[i, j])

print(f"\n  Symmetry classes:")
for cls in sorted(classes.keys()):
    vals = classes[cls]
    print(f"    {cls:10s}: n={len(vals):2d}, mean={np.mean(vals):+.8f}, "
          f"spread={np.std(vals):.2e}")

n_classes = len(classes)
max_spread = max(np.std(v) for v in classes.values())
print(f"  Classes: {n_classes} (expect 7)")
print(f"  Max spread: {max_spread:.2e}")

t2 = time.time()
print(f"  Per-simplex time: {t2-t1:.1f}s")


# ============================================================
# PART 3: AREA HESSIANS
# ============================================================
print("\n" + "=" * 72)
print("PART 3: AREA HESSIANS (2 triangle types)")
print("=" * 72)

def heron_area(l1, l2, l3):
    s = (l1 + l2 + l3) / 2
    return np.sqrt(max(s * (s-l1) * (s-l2) * (s-l3), 0))

L_eq_tri = np.array([1/PHI, 1/PHI, 1/PHI])
L_rad_tri = np.array([1/PHI, np.sqrt(2), np.sqrt(2)])
A_eq = heron_area(*L_eq_tri)
A_rad = heron_area(*L_rad_tri)
print(f"  A_eq  = {A_eq:.8f}")
print(f"  A_rad = {A_rad:.8f}")

h_a = 1e-4  # Optimal for 2nd derivative FD

def area_hessian_3x3(L3):
    A0 = heron_area(*L3)
    H = np.zeros((3, 3))
    fp_a = np.zeros(3)
    fm_a = np.zeros(3)
    for i in range(3):
        Lp = L3.copy(); Lp[i] += h_a
        Lm = L3.copy(); Lm[i] -= h_a
        fp_a[i] = heron_area(*Lp)
        fm_a[i] = heron_area(*Lm)
    for i in range(3):
        H[i, i] = (fp_a[i] - 2*A0 + fm_a[i]) / h_a**2
    for i in range(3):
        for j in range(i+1, 3):
            Lpp = L3.copy(); Lpp[i] += h_a; Lpp[j] += h_a
            Lpm = L3.copy(); Lpm[i] += h_a; Lpm[j] -= h_a
            Lmp = L3.copy(); Lmp[i] -= h_a; Lmp[j] += h_a
            Lmm = L3.copy(); Lmm[i] -= h_a; Lmm[j] -= h_a
            H[i, j] = (heron_area(*Lpp) - heron_area(*Lpm) -
                        heron_area(*Lmp) + heron_area(*Lmm)) / (4*h_a**2)
            H[j, i] = H[i, j]
    return H

HA_eq = area_hessian_3x3(L_eq_tri)
HA_rad = area_hessian_3x3(L_rad_tri)

print(f"\n  Hess_A_eq (3x3):")
for i in range(3):
    print(f"    [{HA_eq[i,0]:+.8f} {HA_eq[i,1]:+.8f} {HA_eq[i,2]:+.8f}]")
print(f"  Hess_A_rad (3x3):")
for i in range(3):
    print(f"    [{HA_rad[i,0]:+.8f} {HA_rad[i,1]:+.8f} {HA_rad[i,2]:+.8f}]")

def get_tri_type(tri):
    if idx_N in tri or idx_S in tri:
        return 'rad'
    return 'eq'

# For radial triangles: sorted (eq_i, eq_j, pole), edges in combinations order:
#   (eq_i,eq_j)=1/phi, (eq_i,pole)=sqrt(2), (eq_j,pole)=sqrt(2)
# This matches L_rad_tri = [1/phi, sqrt(2), sqrt(2)].

t3 = time.time()
print(f"  Area Hessian time: {t3-t2:.1f}s")


# ============================================================
# PART 4: GLOBAL HESSIAN ASSEMBLY
# ============================================================
print("\n" + "=" * 72)
print("PART 4: GLOBAL HESSIAN ASSEMBLY (960x960)")
print("=" * 72)

# H_edge = 2*pi * H_area - H_f

# Assemble H_f
print("  Assembling H_f (1200 simplices)...")
H_f_global = np.zeros((N_E, N_E))
for si in range(1200):
    te = simp_edges_list[si]
    for i_loc in range(10):
        for j_loc in range(10):
            H_f_global[te[i_loc], te[j_loc]] += Hf[i_loc, j_loc]
H_f_global = 0.5 * (H_f_global + H_f_global.T)

# Assemble H_area
print("  Assembling H_area (2640 triangles)...")
H_area_global = np.zeros((N_E, N_E))
for fi in range(N_F):
    tri = all_tris[fi]
    te = tri_edges_list[fi]
    ttype = get_tri_type(tri)
    HA = HA_eq if ttype == 'eq' else HA_rad
    for i_loc in range(3):
        for j_loc in range(3):
            H_area_global[te[i_loc], te[j_loc]] += HA[i_loc, j_loc]
H_area_global = 0.5 * (H_area_global + H_area_global.T)

# Combine
H_edge = 2 * np.pi * H_area_global - H_f_global
H_edge = 0.5 * (H_edge + H_edge.T)

print(f"  |H_f - H_f^T|     = {np.max(np.abs(H_f_global - H_f_global.T)):.2e}")
print(f"  |H_area - H_area^T| = {np.max(np.abs(H_area_global - H_area_global.T)):.2e}")
print(f"  |H_edge - H_edge^T| = {np.max(np.abs(H_edge - H_edge.T)):.2e}")
print(f"  ||H_f||     = {np.linalg.norm(H_f_global):.4f}")
print(f"  ||H_area||  = {np.linalg.norm(H_area_global):.4f}")
print(f"  ||H_edge||  = {np.linalg.norm(H_edge):.4f}")

# --- FD VALIDATION (edge-length space) ---
print("\n  FD VALIDATION (perturbing edge LENGTHS directly)...")

# Build edge-to-simplices and edge-to-triangles maps
edge_to_simps = [[] for _ in range(N_E)]
for si_idx in range(1200):
    for e_idx in simp_edges_list[si_idx]:
        edge_to_simps[e_idx].append(si_idx)

edge_to_tris = [[] for _ in range(N_E)]
for fi_idx in range(N_F):
    for e_idx in tri_edges_list[fi_idx]:
        edge_to_tris[e_idx].append(fi_idx)

# Reference edge lengths
L0_all = np.array([np.linalg.norm(verts_5d[e[0]] - verts_5d[e[1]]) for e in all_edges])

def S_partial(L_vec, affected_simps, affected_tris):
    """Compute partial Regge action involving only specified simplices and triangles."""
    area_part = 0.0
    for fi_idx in affected_tris:
        te = tri_edges_list[fi_idx]
        area_part += heron_area(L_vec[te[0]], L_vec[te[1]], L_vec[te[2]])

    f_part = 0.0
    for si_idx in affected_simps:
        L10 = np.array([L_vec[e] for e in simp_edges_list[si_idx]])
        f_part += compute_f_simplex(L10)

    return 2 * np.pi * area_part - f_part

# Verify S0 from edge lengths
S0_edge = S_partial(L0_all, list(range(1200)), list(range(N_F)))
print(f"  S_Regge (from edge lengths): {S0_edge:.6f}")

# Pick test edge pairs
eq_indices = [i for i, t in enumerate(edge_type) if t == 'eq']
nr_indices = [i for i, t in enumerate(edge_type) if t == 'Nr']

test_pairs = []
# eq-eq: two equatorial edges sharing a simplex
for ea in eq_indices[:30]:
    for eb in eq_indices[:30]:
        if ea < eb:
            common_simps = set(edge_to_simps[ea]) & set(edge_to_simps[eb])
            if len(common_simps) > 0:
                test_pairs.append((ea, eb, 'eq-eq'))
                break
    if len(test_pairs) > 0:
        break
# Nr-Nr: two N-radial edges sharing a simplex
for ea in nr_indices[:10]:
    for eb in nr_indices[:10]:
        if ea < eb:
            common_simps = set(edge_to_simps[ea]) & set(edge_to_simps[eb])
            if len(common_simps) > 0:
                test_pairs.append((ea, eb, 'Nr-Nr'))
                break
    if len(test_pairs) > 1:
        break
# eq-Nr
for ea in eq_indices[:10]:
    for eb in nr_indices[:10]:
        common_simps = set(edge_to_simps[ea]) & set(edge_to_simps[eb])
        if len(common_simps) > 0:
            test_pairs.append((ea, eb, 'eq-Nr'))
            break
    if len(test_pairs) > 2:
        break
# diagonal
test_pairs.append((eq_indices[0], eq_indices[0], 'eq-diag'))

h_val = 1e-4  # Match per-simplex FD step
max_rel_err = 0.0
print(f"  Testing {len(test_pairs)} edge pairs...")
for ea, eb, label in test_pairs:
    # Find affected simplices and triangles (those containing ea OR eb)
    aff_simps = sorted(set(edge_to_simps[ea]) | set(edge_to_simps[eb]))
    aff_tris = sorted(set(edge_to_tris[ea]) | set(edge_to_tris[eb]))

    if ea == eb:
        # Diagonal: d^2S/dl_a^2
        Lp = L0_all.copy(); Lp[ea] += h_val
        Lm = L0_all.copy(); Lm[ea] -= h_val
        Sp = S_partial(Lp, aff_simps, aff_tris)
        S0p = S_partial(L0_all, aff_simps, aff_tris)
        Sm = S_partial(Lm, aff_simps, aff_tris)
        H_fd_val = (Sp - 2*S0p + Sm) / h_val**2
    else:
        # Off-diagonal: mixed partial
        Lpp = L0_all.copy(); Lpp[ea] += h_val; Lpp[eb] += h_val
        Lpm = L0_all.copy(); Lpm[ea] += h_val; Lpm[eb] -= h_val
        Lmp = L0_all.copy(); Lmp[ea] -= h_val; Lmp[eb] += h_val
        Lmm = L0_all.copy(); Lmm[ea] -= h_val; Lmm[eb] -= h_val
        H_fd_val = (S_partial(Lpp, aff_simps, aff_tris)
                    - S_partial(Lpm, aff_simps, aff_tris)
                    - S_partial(Lmp, aff_simps, aff_tris)
                    + S_partial(Lmm, aff_simps, aff_tris)) / (4*h_val**2)

    H_assembled = H_edge[ea, eb]
    err = abs(H_assembled - H_fd_val)
    denom = max(abs(H_fd_val), abs(H_assembled), 1e-15)
    rel = err / denom
    max_rel_err = max(max_rel_err, rel)
    print(f"    {label}: assembled={H_assembled:+.6f}, FD={H_fd_val:+.6f}, "
          f"err={err:.2e}, rel={rel:.2e}")

print(f"  Max relative error: {max_rel_err:.2e}")
if max_rel_err < 1e-3:
    print(f"  FD VALIDATION: PASSED")
elif max_rel_err < 0.1:
    print(f"  FD VALIDATION: MARGINAL")
else:
    print(f"  FD VALIDATION: CHECK NEEDED")

t4 = time.time()
print(f"  Assembly + validation time: {t4-t3:.1f}s")


# ============================================================
# PART 5: FULL SPECTRUM
# ============================================================
print("\n" + "=" * 72)
print("PART 5: FULL SPECTRUM (960x960)")
print("=" * 72)

evals_H = np.linalg.eigvalsh(H_edge)
n_pos = np.sum(evals_H > 0.01)
n_neg = np.sum(evals_H < -0.01)
n_zero_H = np.sum(np.abs(evals_H) < 0.01)

print(f"  Spectrum: {n_pos} positive, {n_neg} negative, {n_zero_H} near-zero")
print(f"  Range: [{evals_H[0]:.4f}, {evals_H[-1]:.4f}]")
print(f"  ||H_edge||: {np.linalg.norm(H_edge):.4f}")
print(f"  Tr(H_edge): {np.trace(H_edge):.4f}")

print(f"\n  Lowest 20 eigenvalues:")
for k in range(min(20, len(evals_H))):
    print(f"    {k:3d}: {evals_H[k]:+.6f}")

print(f"\n  Highest 20 eigenvalues:")
for k in range(max(0, len(evals_H)-20), len(evals_H)):
    print(f"    {k:3d}: {evals_H[k]:+.6f}")

# Group by multiplicity
groups = []
i = 0
while i < len(evals_H):
    val = evals_H[i]
    count = 1
    while i + count < len(evals_H) and abs(evals_H[i+count] - val) < 0.01:
        count += 1
    groups.append((val, count))
    i += count

print(f"\n  Distinct levels: {len(groups)}")
print(f"  {'level':>6s} {'eigenvalue':>14s} {'mult':>6s}")
for gi, (val, mult) in enumerate(groups):
    if gi < 25 or gi >= len(groups) - 10:
        print(f"  {gi:6d} {val:+14.4f} {mult:6d}")
    elif gi == 25:
        print(f"  {'...':>6s}")

print(f"\n  Multiplicities: {[g[1] for g in groups]}")

t5 = time.time()
print(f"  Spectrum time: {t5-t4:.1f}s")


# ============================================================
# PART 6: COEXACT PROJECTION & PHYSICAL SPECTRUM
# ============================================================
print("\n" + "=" * 72)
print("PART 6: COEXACT (PHYSICAL) SPECTRUM")
print("=" * 72)

H_co = P_coexact @ H_edge @ P_coexact
evals_Hco_full = np.linalg.eigvalsh(H_co)
evals_Hco = evals_Hco_full[np.abs(evals_Hco_full) > 0.01]

n_pos_co = np.sum(evals_Hco > 0.01)
n_neg_co = np.sum(evals_Hco < -0.01)
norm_Hco = np.linalg.norm(H_co)

print(f"  Coexact Hessian: {len(evals_Hco)} nonzero eigenvalues")
print(f"  Positive: {n_pos_co}, Negative: {n_neg_co}")
print(f"  ||H_co|| = {norm_Hco:.4f}")

if len(evals_Hco) > 0:
    print(f"  Range: [{evals_Hco.min():.6f}, {evals_Hco.max():.6f}]")

    if n_neg_co == 0:
        print(f"\n  STABILITY: Round S^4 IS stable (all physical eigenvalues > 0)")
    else:
        print(f"\n  STABILITY: Coexact Hessian INDEFINITE ({n_neg_co} negative modes)")
        print(f"  (Conformal factor problem in Euclidean gravity - EXPECTED)")

print(f"\n  Lowest 20 coexact eigenvalues:")
for k in range(min(20, len(evals_Hco))):
    print(f"    {k:3d}: {evals_Hco[k]:+.6f}")

print(f"\n  Highest 10 coexact eigenvalues:")
for k in range(max(0, len(evals_Hco)-10), len(evals_Hco)):
    print(f"    {k:3d}: {evals_Hco[k]:+.6f}")

# Group coexact by multiplicity
groups_co = []
sorted_Hco = np.sort(evals_Hco)
i = 0
while i < len(sorted_Hco):
    val = sorted_Hco[i]
    count = 1
    while i + count < len(sorted_Hco) and abs(sorted_Hco[i+count] - val) < 0.01:
        count += 1
    groups_co.append((val, count))
    i += count

print(f"\n  Coexact levels: {len(groups_co)}")
print(f"  {'level':>6s} {'eigenvalue':>14s} {'mult':>6s}")
for gi, (val, mult) in enumerate(groups_co):
    if gi < 20 or gi >= len(groups_co) - 5:
        print(f"  {gi:6d} {val:+14.6f} {mult:6d}")
    elif gi == 20:
        print(f"  {'...':>6s}")

# Exact sector
H_ex = P_exact @ H_edge @ P_exact
evals_Hex = np.linalg.eigvalsh(H_ex)
evals_Hex_nz = evals_Hex[np.abs(evals_Hex) > 0.01]
n_pos_ex = np.sum(evals_Hex_nz > 0.01)
n_neg_ex = np.sum(evals_Hex_nz < -0.01)
print(f"\n  Exact sector: {len(evals_Hex_nz)} nonzero evals")
print(f"  Positive: {n_pos_ex}, Negative: {n_neg_ex}")

ratio_HL = norm_Hco / np.linalg.norm(L_coexact) if np.linalg.norm(L_coexact) > 0 else 0
print(f"\n  ||H_co|| / ||L_coexact|| = {ratio_HL:.4f}")

print(f"\n  Continuum Lichnerowicz on S^4:")
print(f"    n=2: lambda=8  (mult 14)")
print(f"    n=3: lambda=16 (mult 30)")
print(f"    n=4: lambda=26 (mult 55)")
print(f"    n=5: lambda=38 (mult 91)")

t6 = time.time()
print(f"  Coexact time: {t6-t5:.1f}s")


# ============================================================
# PART 7: FRAMEWORK CONNECTIONS
# ============================================================
print("\n" + "=" * 72)
print("PART 7: FRAMEWORK CONNECTIONS")
print("=" * 72)

a1 = 5
N = 120

print(f"\n  Counts:")
print(f"    Full:    {n_pos} pos, {n_neg} neg, {n_zero_H} zero")
print(f"    Coexact: {n_pos_co} pos, {n_neg_co} neg")
print(f"    Exact:   {n_pos_ex} pos, {n_neg_ex} neg")

# Test framework numbers
framework_tests = {
    'a1^2+1=26': a1**2 + 1,
    'a1^2=25': a1**2,
    'a1*19=95': a1*19,
    '5*19=95': 5*19,
    'N-1=119': N-1,
    'N=120': N,
    'N+1=121': N+1,
}
print(f"\n  Framework number matches:")
all_counts = {'n_pos': n_pos, 'n_neg': n_neg, 'n_zero': n_zero_H,
              'n_pos_co': n_pos_co, 'n_neg_co': n_neg_co,
              'n_pos_ex': n_pos_ex, 'n_neg_ex': n_neg_ex}
for label, val in framework_tests.items():
    matches = [f"{k}={v}" for k, v in all_counts.items() if v == val]
    if matches:
        print(f"    {label} -> {', '.join(matches)}")

# Multiplicities
mults = [g[1] for g in groups]
sq_mults = [m for m in mults if int(round(np.sqrt(m)))**2 == m]
print(f"\n  Full multiplicities: {mults}")
print(f"  Perfect squares: {sq_mults}")
mults_co = [g[1] for g in groups_co]
print(f"  Coexact multiplicities: {mults_co}")

# Traces
Tr_H = np.trace(H_edge)
Tr_Hco = np.trace(H_co)
Tr_Hex = np.trace(H_ex)
print(f"\n  Traces:")
print(f"    Tr(H)    = {Tr_H:.4f}")
print(f"    Tr(H_co) = {Tr_Hco:.4f}")
print(f"    Tr(H_ex) = {Tr_Hex:.4f}")

Tr_H2 = np.trace(H_edge @ H_edge)
Tr_Hco2 = np.trace(H_co @ H_co)
print(f"    Tr(H^2)    = {Tr_H2:.4f}")
print(f"    Tr(H_co^2) = {Tr_Hco2:.4f}")

# Perturbation theory
print(f"\n  Coupling: ||H_co||/||L_co|| = {ratio_HL:.4f}")
if ratio_HL < 0.1:
    pt_status = "WEAKLY COUPLED"
elif ratio_HL < 1.0:
    pt_status = "MARGINAL"
else:
    pt_status = "STRONGLY COUPLED"
print(f"  -> {pt_status}")

# Mass gaps
if len(evals_Hco) > 0:
    mass_gap_co = np.min(np.abs(evals_Hco))
    mass_gap_L = lambda_co.min()
    print(f"\n  Mass gaps:")
    print(f"    L_coexact: {mass_gap_L:.6f}")
    print(f"    H_coexact: {mass_gap_co:.6f}")
    for name, val in [('7-4*phi', 7-4*PHI), ('phi', PHI), ('phi^2', PHI**2),
                      ('2*pi/a1', 2*np.pi/a1), ('1/phi', 1/PHI)]:
        print(f"    {name:>10s} = {val:.6f} (ratio = {mass_gap_co/val:.4f})")

t7 = time.time()
print(f"  Framework time: {t7-t6:.1f}s")


# ============================================================
# PART 8: SUMMARY
# ============================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

total_time = time.time() - t0

if n_neg_co == 0:
    stability = "POSITIVE DEFINITE (stable)"
else:
    stability = f"INDEFINITE ({n_neg_co} neg modes -- conformal factor problem, EXPECTED)"

print(f"""
4D REGGE HESSIAN ON S^4 (DOUBLE-CONE OVER 600-CELL)
====================================================

TRIANGULATION:
  S^4: V=122, E={N_E}, F={N_F}, T={N_T}, S4=1200
  Euler: {chi}, ALL congruent

PER-SIMPLEX HESSIAN:
  10x10 universal, {n_classes} symmetry classes, max spread {max_spread:.2e}

GLOBAL HESSIAN (H = 2*pi*H_area - H_f):
  ||H_edge|| = {np.linalg.norm(H_edge):.4f}
  FD validation: rel err = {max_rel_err:.2e}

FULL SPECTRUM:
  {n_pos} pos, {n_neg} neg, {n_zero_H} zero
  Range: [{evals_H[0]:.4f}, {evals_H[-1]:.4f}]
  Levels: {len(groups)}, Mults: {[g[1] for g in groups]}

COEXACT (PHYSICAL):
  {n_pos_co} pos, {n_neg_co} neg
  {stability}""")
if len(evals_Hco) > 0:
    print(f"  Range: [{evals_Hco.min():.6f}, {evals_Hco.max():.6f}]")
    print(f"  Levels: {len(groups_co)}, Mults: {[g[1] for g in groups_co]}")

print(f"""
EXACT SECTOR:
  {n_pos_ex} pos, {n_neg_ex} neg

COUPLING:
  ||H_co||/||L_co|| = {ratio_HL:.4f} -> {pt_status}

TOTAL TIME: {total_time:.1f}s
""")

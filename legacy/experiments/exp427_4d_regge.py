"""
exp427_4d_regge.py — 4D Regge Calculus on S^4 (double-cone over 600-cell)
=========================================================================
Emergent gravity: extract G_4D from 4D deficit angles at triangular hinges.

S^4 triangulation (double-cone):
  122 vertices: 120 equatorial (600-cell on S^3) + N pole + S pole
  1200 four-simplices: 600 north-cone + 600 south-cone
  ALL congruent: 4 edges sqrt(2), 6 edges 1/phi
  2640 triangular hinges: 1200 equatorial + 720 N-radial + 720 S-radial

4D Regge action: S = sum_triangles A_triangle * epsilon_triangle
  epsilon = 2*pi - sum(dihedral angles of surrounding 4-simplices)
"""

import numpy as np
import itertools
import time

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 7.2973525693e-3

print("=" * 72)
print("EXP-427: 4D REGGE CALCULUS ON S^4 (DOUBLE-CONE OVER 600-CELL)")
print("=" * 72)

# ============================================================
# PART 1: BUILD 600-CELL + S^4 TRIANGULATION
# ============================================================
print("\n" + "=" * 72)
print("PART 1: BUILD 600-CELL + S^4 TRIANGULATION")
print("=" * 72)
t0 = time.time()

# --- Build 600-cell vertices on unit S^3 in R^4 ---
verts_4d = []
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0,0,0,0]; v[i] = s
        verts_4d.append(v)
for signs in itertools.product([0.5, -0.5], repeat=4):
    verts_4d.append(list(signs))
vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]
even_perms = [p for p in itertools.permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]) % 2 == 0]
for perm in even_perms:
    base = [vals_base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for signs in itertools.product([1,-1], repeat=len(nz)):
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

# Build edges (chord distance 1/phi)
target = 1.0 / PHI
edges = []
adj = {i: set() for i in range(120)}
for i in range(120):
    for j in range(i+1, 120):
        if abs(np.linalg.norm(verts_4d[i] - verts_4d[j]) - target) < 1e-8:
            edges.append((i, j))
            adj[i].add(j); adj[j].add(i)
assert len(edges) == 720
edge_set = {(min(i,j), max(i,j)) for i,j in edges}
edge_idx = {}
for idx, (i,j) in enumerate(edges):
    edge_idx[(i,j)] = idx; edge_idx[(j,i)] = idx

# Build cells (tetrahedra = 4-cliques)
cells = set()
for (i,j) in edges:
    common = adj[i] & adj[j]
    for k in common:
        for l in common:
            if l > k and (min(k,l), max(k,l)) in edge_set:
                cells.add(tuple(sorted([i,j,k,l])))
cells = sorted(cells)
assert len(cells) == 600

# Build faces
faces = set()
for cell in cells:
    for tri in itertools.combinations(cell, 3):
        faces.add(tri)
faces = sorted(faces)
assert len(faces) == 1200

# Face-to-cells map (each triangle in exactly 2 tets)
face_idx = {f: i for i, f in enumerate(faces)}
face_to_cells = [[] for _ in range(1200)]
for ci, cell in enumerate(cells):
    for tri in itertools.combinations(cell, 3):
        face_to_cells[face_idx[tri]].append(ci)
assert all(len(fc) == 2 for fc in face_to_cells)

# Edge-to-cells map
edge_to_cells = [[] for _ in range(720)]
for ci, cell in enumerate(cells):
    for pair in itertools.combinations(cell, 2):
        ei = edge_idx.get((pair[0], pair[1]))
        if ei is not None:
            edge_to_cells[ei].append(ci)
assert all(len(ec) == 5 for ec in edge_to_cells)

print(f"  600-cell: V=120, E=720, F=1200, C=600")
print(f"  Edge length: 1/phi = {1/PHI:.6f}")

# --- S^4 double-cone: embed in R^5, add N and S poles ---
N_eq = 120
idx_N, idx_S = N_eq, N_eq + 1
verts_5d = np.zeros((122, 5))
verts_5d[:120, :4] = verts_4d
verts_5d[idx_N] = [0, 0, 0, 0, 1]   # North pole
verts_5d[idx_S] = [0, 0, 0, 0, -1]  # South pole

# Verify edge lengths
d_radial = np.linalg.norm(verts_5d[idx_N] - verts_5d[0])
print(f"  Radial edge (pole to equator): {d_radial:.6f} (expect sqrt(2)={np.sqrt(2):.6f})")
d_poles = np.linalg.norm(verts_5d[idx_N] - verts_5d[idx_S])
print(f"  Pole-pole distance: {d_poles:.6f} (NOT an edge)")

# Build 1200 four-simplices
four_simplices = []
for tet in cells:
    four_simplices.append((idx_N,) + tet)  # North cone
    four_simplices.append((idx_S,) + tet)  # South cone
four_simplices = np.array(four_simplices)
assert len(four_simplices) == 1200

print(f"\n  S^4 triangulation:")
print(f"    Vertices: 122 (120 equatorial + 2 poles)")
print(f"    Four-simplices: {len(four_simplices)}")

# Count edges, triangles, tetrahedra of the 4D complex
all_edges_5d = set()
all_tris_5d = set()
all_tets_5d = set()
for s4 in four_simplices:
    for e in itertools.combinations(s4, 2):
        all_edges_5d.add(tuple(sorted(e)))
    for t in itertools.combinations(s4, 3):
        all_tris_5d.add(tuple(sorted(t)))
    for tt in itertools.combinations(s4, 4):
        all_tets_5d.add(tuple(sorted(tt)))

n_edges_4d = len(all_edges_5d)
n_tris_4d = len(all_tris_5d)
n_tets_4d = len(all_tets_5d)
chi = 122 - n_edges_4d + n_tris_4d - n_tets_4d + 1200
print(f"    Edges: {n_edges_4d} (720 eq + 120 N-rad + 120 S-rad = {720+240})")
print(f"    Triangles: {n_tris_4d} (1200 eq + 720 N-rad + 720 S-rad = {1200+1440})")
print(f"    Tetrahedra: {n_tets_4d}")
print(f"    Euler char (V-E+F-T+S4): {chi} (S^4 expects 2)")

t1 = time.time()
print(f"  Build time: {t1-t0:.1f}s")


# ============================================================
# PART 2: UNIVERSAL 4-SIMPLEX DIHEDRAL ANGLES
# ============================================================
print("\n" + "=" * 72)
print("PART 2: UNIVERSAL 4-SIMPLEX DIHEDRAL ANGLES")
print("=" * 72)

def dihedral_angle_4d(p0, p1, p2, p3, p4):
    """
    Dihedral angle at triangle (p0,p1,p2) in 4-simplex (p0,p1,p2,p3,p4).
    Projects p3, p4 into the 2D orthogonal complement of the triangle
    within the simplex's 4D affine span.
    """
    u1 = p1 - p0
    u2 = p2 - p0
    u3 = p3 - p0
    u4 = p4 - p0

    # Gram-Schmidt: {e1,e2} span triangle, {e3,e4} span complement
    e1 = u1 / np.linalg.norm(u1)

    u2o = u2 - np.dot(u2, e1) * e1
    e2 = u2o / np.linalg.norm(u2o)

    u3o = u3 - np.dot(u3, e1)*e1 - np.dot(u3, e2)*e2
    e3 = u3o / np.linalg.norm(u3o)

    u4o = u4 - np.dot(u4, e1)*e1 - np.dot(u4, e2)*e2 - np.dot(u4, e3)*e3
    e4 = u4o / np.linalg.norm(u4o)

    # Project p3, p4 onto {e3, e4}
    d3 = p3 - p0
    d3c = np.array([np.dot(d3, e3), np.dot(d3, e4)])

    d4 = p4 - p0
    d4c = np.array([np.dot(d4, e3), np.dot(d4, e4)])

    cos_a = np.dot(d3c, d4c) / (np.linalg.norm(d3c) * np.linalg.norm(d4c))
    return np.arccos(np.clip(cos_a, -1.0, 1.0))


# Use a representative northern 4-simplex
# Pick first boundary tet -> 4-simplex (N, a, b, c, d)
rep_tet = cells[0]
p_N = verts_5d[idx_N]
pa, pb, pc, pd = [verts_5d[v] for v in rep_tet]

# Verify edge lengths
print(f"  Representative 4-simplex: (N, {rep_tet[0]}, {rep_tet[1]}, {rep_tet[2]}, {rep_tet[3]})")
radial_lens = [np.linalg.norm(p_N - p) for p in [pa, pb, pc, pd]]
eq_lens = []
for i, pi in enumerate([pa, pb, pc, pd]):
    for j, pj in enumerate([pa, pb, pc, pd]):
        if j > i:
            eq_lens.append(np.linalg.norm(pi - pj))
print(f"  Radial edges: {np.mean(radial_lens):.6f} +/- {np.std(radial_lens):.2e}")
print(f"  Equatorial edges: {np.mean(eq_lens):.6f} +/- {np.std(eq_lens):.2e}")

# Compute ALL 10 dihedral angles of this 4-simplex
simplex_verts = [p_N, pa, pb, pc, pd]
labels = ['N', 'a', 'b', 'c', 'd']

print(f"\n  Dihedral angles at each triangle:")
theta_A_list = []  # Radial triangles (contain N)
theta_B_list = []  # Equatorial triangles (no N)

for tri in itertools.combinations(range(5), 3):
    others = [k for k in range(5) if k not in tri]
    p = [simplex_verts[k] for k in tri]
    q = [simplex_verts[k] for k in others]
    theta = dihedral_angle_4d(p[0], p[1], p[2], q[0], q[1])
    tri_label = ''.join(labels[k] for k in tri)
    has_N = (0 in tri)
    if has_N:
        theta_A_list.append(theta)
    else:
        theta_B_list.append(theta)
    print(f"    ({tri_label}): theta = {np.degrees(theta):.4f} deg = {theta:.6f} rad"
          f"  [{'radial' if has_N else 'equatorial'}]")

theta_A = np.mean(theta_A_list)
theta_B = np.mean(theta_B_list)
print(f"\n  theta_A (radial, 6 per simplex):    {np.degrees(theta_A):.4f} deg = {theta_A:.6f} rad")
print(f"    spread: {np.std(theta_A_list):.2e} (should be ~0)")
print(f"  theta_B (equatorial, 4 per simplex): {np.degrees(theta_B):.4f} deg = {theta_B:.6f} rad")
print(f"    spread: {np.std(theta_B_list):.2e} (should be ~0)")

# Verify with a second 4-simplex
rep_tet2 = cells[100]
pa2, pb2, pc2, pd2 = [verts_5d[v] for v in rep_tet2]
sv2 = [p_N, pa2, pb2, pc2, pd2]
theta_check = []
for tri in itertools.combinations(range(5), 3):
    others = [k for k in range(5) if k not in tri]
    p = [sv2[k] for k in tri]
    q = [sv2[k] for k in others]
    theta_check.append(dihedral_angle_4d(p[0], p[1], p[2], q[0], q[1]))
print(f"\n  Cross-check (simplex #100): max deviation = {max(abs(a-b) for a,b in zip(sorted(theta_check), sorted(theta_A_list + theta_B_list))):.2e}")


# ============================================================
# PART 3: DEFICIT ANGLES AND REGGE ACTION
# ============================================================
print("\n" + "=" * 72)
print("PART 3: DEFICIT ANGLES AND 4D REGGE ACTION")
print("=" * 72)

# Triangle areas
A_eq = (np.sqrt(3)/4) * (1/PHI)**2  # Equilateral, side 1/phi
h_rad = np.sqrt(2 - (1/(2*PHI))**2)  # Height of isoceles triangle
A_rad = 0.5 * (1/PHI) * h_rad        # Isoceles: sides sqrt(2), sqrt(2), 1/phi

print(f"  A_eq  (equatorial triangle): {A_eq:.6f}")
print(f"  A_rad (radial triangle):     {A_rad:.6f}")

# Deficit angles
# Equatorial triangles: 4 surrounding 4-simplices (2N + 2S), each contributes theta_B
eps_eq = 2*np.pi - 4*theta_B
# Radial triangles: 5 surrounding 4-simplices (all N or all S), each contributes theta_A
eps_rad = 2*np.pi - 5*theta_A

print(f"\n  Equatorial deficit: eps_eq  = 2pi - 4*theta_B = {eps_eq:.6f} rad = {np.degrees(eps_eq):.4f} deg")
print(f"  Radial deficit:    eps_rad = 2pi - 5*theta_A = {eps_rad:.6f} rad = {np.degrees(eps_rad):.4f} deg")

# 4D Regge action
S_regge = 1200 * A_eq * eps_eq + 1440 * A_rad * eps_rad
print(f"\n  S_Regge = 1200*A_eq*eps_eq + 1440*A_rad*eps_rad")
print(f"         = {1200*A_eq*eps_eq:.4f} + {1440*A_rad*eps_rad:.4f}")
print(f"         = {S_regge:.4f}")

# Continuum comparison: for round S^4 of radius 1
# integral(R * sqrt(g) d^4x) = R * Vol(S^4) = 12 * 8*pi^2/3 = 32*pi^2
S_continuum = 32 * np.pi**2
print(f"\n  Continuum: integral(R*sqrt(g)) = 32*pi^2 = {S_continuum:.4f}")
print(f"  Ratio S_Regge/S_continuum = {S_regge/S_continuum:.6f}")

# ============================================================
# PART 4: 4D VOLUME AND COSMOLOGICAL CONSTANT
# ============================================================
print("\n" + "=" * 72)
print("PART 4: 4D VOLUME AND COSMOLOGICAL CONSTANT")
print("=" * 72)

# Compute 4D volume of each 4-simplex
# For vertices in R^5, the 4-volume uses: V = sqrt(det(G)) / 4!
# where G_ij = (p_i - p_0) . (p_j - p_0) is the 4x4 Gram matrix
rep_simplex = four_simplices[0]
ps = [verts_5d[v] for v in rep_simplex]
edge_vecs = np.array([ps[i] - ps[0] for i in range(1, 5)])  # 4x5
gram = edge_vecs @ edge_vecs.T  # 4x4 Gram matrix
V_simplex = np.sqrt(abs(np.linalg.det(gram))) / 24
V_total = 1200 * V_simplex
V_continuum = 8 * np.pi**2 / 3  # Vol(S^4, r=1)

print(f"  Volume per 4-simplex: {V_simplex:.6f}")
print(f"  Total 4D volume: 1200 * {V_simplex:.6f} = {V_total:.6f}")
print(f"  Continuum Vol(S^4): 8*pi^2/3 = {V_continuum:.6f}")
print(f"  Ratio V_Regge/V_cont = {V_total/V_continuum:.6f}")

# Cross-check with another simplex
ps2 = [verts_5d[v] for v in four_simplices[599]]
ev2 = np.array([ps2[i] - ps2[0] for i in range(1, 5)])
V_check = np.sqrt(abs(np.linalg.det(ev2 @ ev2.T))) / 24
print(f"  Cross-check (simplex #599): V = {V_check:.6f} (diff = {abs(V_check-V_simplex):.2e})")

# Average scalar curvature from Regge
# S_Regge = integral(R*sqrt(g)) ≈ R_avg * V_total
if V_total > 1e-10:
    R_avg = S_regge / V_total
    print(f"\n  Average scalar curvature: R_avg = S_Regge/V = {R_avg:.4f}")
    print(f"  Continuum R(S^4) = 12.0000")
    print(f"  Ratio: {R_avg/12:.4f}")

# Extract Lambda from EOM: dS/dl = 2*Lambda * dV/dl
# Compute dV/dl for equatorial and radial edges by FD
delta = 1e-6

def compute_total_vol(v5d):
    """Total 4D volume of the triangulation."""
    vol = 0.0
    for s4 in four_simplices:
        pp = [v5d[k] for k in s4]
        ev = np.array([pp[i] - pp[0] for i in range(1, 5)])
        gram = ev @ ev.T
        vol += np.sqrt(abs(np.linalg.det(gram))) / 24
    return vol

# dV/dl for equatorial edge
e_eq = edges[0]
vi, vj = e_eq
d_eq = verts_5d[vj] - verts_5d[vi]
d_eq = d_eq / np.linalg.norm(d_eq)
vp = verts_5d.copy(); vp[vj] = verts_5d[vj] + delta * d_eq
vm = verts_5d.copy(); vm[vj] = verts_5d[vj] - delta * d_eq
dVdl_eq = (compute_total_vol(vp) - compute_total_vol(vm)) / (2*delta)

# dV/dl for N-radial edge
vi_r = 0
d_Nr = verts_5d[vi_r] - verts_5d[idx_N]
d_Nr = d_Nr / np.linalg.norm(d_Nr)
vp = verts_5d.copy(); vp[vi_r] = verts_5d[vi_r] + delta * d_Nr
vm = verts_5d.copy(); vm[vi_r] = verts_5d[vi_r] - delta * d_Nr
dVdl_Nr = (compute_total_vol(vp) - compute_total_vol(vm)) / (2*delta)

print(f"\n  dV/dl (equatorial): {dVdl_eq:.6f}")
print(f"  dV/dl (N-radial):   {dVdl_Nr:.6f}")

# Lambda from EOM: dS/dl = 2*Lambda * dV/dl
# (This is the Regge form of Einstein eq R_ab = Lambda*g_ab for constant curvature)
if abs(dVdl_eq) > 1e-10:
    Lambda_eq = grad_eq / (2 * dVdl_eq) if 'grad_eq' in dir() else 0
    print(f"\n  (Will extract Lambda after gradient computation in Part 6)")

# ============================================================
# PART 5: LINEAR RESPONSE — EXTRACT G FROM POINT SOURCE
# ============================================================
print("\n" + "=" * 72)
print("PART 5: LINEAR RESPONSE (POINT MASS -> G_4D)")
print("=" * 72)

# Strategy: The Regge EOM with matter are
#   dS_Regge/dl_e = 16*pi*G * dS_matter/dl_e
# For a point mass M at vertex v: S_matter = M * (proper time at v)
# In Euclidean: S_matter = M * (4D volume associated to v)
#
# Linearized around the background:
#   H * delta_l = 16*pi*G * J
# where H_ee' = d^2 S_Regge / dl_e dl_e' and J_e = dS_matter/dl_e
#
# The response delta_l encodes how geometry deforms.
# From delta_l -> delta_epsilon -> Newtonian potential -> G

# First, compute the Hessian by finite differences
# Since ALL 4-simplices are congruent, we can use the tet-local trick:
# H[ea, eb] = sum over 4-simplices containing both ea, eb of (local Hessian contribution)

# For each 4-simplex: 10 edges, 10 triangles
# S_simplex = sum_triangles A(edges) * epsilon(edges, neighbors)
# But epsilon involves OTHER simplices too. So the Hessian is non-local.

# Simpler approach: direct finite difference of the total Regge action.
# Perturb edge e by delta, recompute all deficit angles, get S(l+delta) and S(l-delta).

# But we have 960 edges... and each perturbation requires recomputing dihedral angles
# for all 4-simplices containing that edge.

# Actually: by symmetry, there are only 3 types of edges:
#   1. Equatorial (720): between equatorial vertices, length 1/phi
#   2. N-radial (120): from N to equatorial vertex, length sqrt(2)
#   3. S-radial (120): from S to equatorial vertex, length sqrt(2)

# By the icosahedral symmetry of the 600-cell:
#   - All equatorial edges are equivalent
#   - All N-radial edges are equivalent
#   - All S-radial edges are equivalent (and same as N-radial by N<->S symmetry)

# So the Hessian has a 3x3 block structure (by symmetry classes).
# Actually within each class, all edges are equivalent, so:
#   dS/dl_e = same for all e in same class
#   d^2S/dl_e dl_e' depends on the symmetry relation between e and e'

# Let's just compute dS/dl_e first (should be 0 for equilibrium = EOM satisfied)

def compute_S_regge(verts_5d_mod, four_simplices_list, all_tris_dict):
    """
    Compute 4D Regge action for given vertex positions.
    all_tris_dict: {frozenset(tri): [list of simplex indices]}
    """
    S = 0.0
    for tri_key, simplex_indices in all_tris_dict.items():
        tri = list(tri_key)
        # Triangle area
        v0, v1, v2 = verts_5d_mod[tri[0]], verts_5d_mod[tri[1]], verts_5d_mod[tri[2]]
        cross = np.zeros(10)  # Cross product in 5D needs special handling
        # Use the standard formula: A = 0.5 * |u x v| where u=v1-v0, v=v2-v0
        u = v1 - v0
        v = v2 - v0
        # In any dimension: A = 0.5 * sqrt(|u|^2*|v|^2 - (u.v)^2)
        A = 0.5 * np.sqrt(np.dot(u,u)*np.dot(v,v) - np.dot(u,v)**2)

        # Deficit angle: 2*pi - sum of dihedral angles
        total_dihedral = 0.0
        for si in simplex_indices:
            s = four_simplices_list[si]
            others = [k for k in s if k not in tri_key]
            # Dihedral angle at this triangle in this 4-simplex
            p0, p1, p2 = verts_5d_mod[tri[0]], verts_5d_mod[tri[1]], verts_5d_mod[tri[2]]
            p3, p4 = verts_5d_mod[others[0]], verts_5d_mod[others[1]]
            theta = dihedral_angle_4d(p0, p1, p2, p3, p4)
            total_dihedral += theta

        eps = 2*np.pi - total_dihedral
        S += A * eps
    return S


# Build triangle-to-simplices dict
print("  Building triangle-to-simplices map...")
tri_to_simp = {}
for si, s4 in enumerate(four_simplices):
    for tri in itertools.combinations(s4, 3):
        key = frozenset(tri)
        if key not in tri_to_simp:
            tri_to_simp[key] = []
        tri_to_simp[key].append(si)

print(f"  Total triangles: {len(tri_to_simp)}")

# Classify triangles
n_eq_tri = 0
n_Nrad_tri = 0
n_Srad_tri = 0
for tri_key in tri_to_simp:
    tri = list(tri_key)
    has_N = idx_N in tri
    has_S = idx_S in tri
    if not has_N and not has_S:
        n_eq_tri += 1
    elif has_N and not has_S:
        n_Nrad_tri += 1
    elif has_S and not has_N:
        n_Srad_tri += 1
print(f"  Equatorial triangles: {n_eq_tri}")
print(f"  N-radial triangles: {n_Nrad_tri}")
print(f"  S-radial triangles: {n_Srad_tri}")

# Compute full Regge action (verify against analytical formula)
print("\n  Computing full Regge action (verify)...")
t2 = time.time()
S_full = compute_S_regge(verts_5d, four_simplices, tri_to_simp)
t3 = time.time()
print(f"  S_Regge (full computation): {S_full:.6f}")
print(f"  S_Regge (analytical):       {S_regge:.6f}")
print(f"  Difference: {abs(S_full - S_regge):.2e}")
print(f"  Computation time: {t3-t2:.1f}s")


# ============================================================
# PART 6: GRADIENT (REGGE EQUATIONS OF MOTION)
# ============================================================
print("\n" + "=" * 72)
print("PART 6: GRADIENT OF REGGE ACTION")
print("=" * 72)

# Check if the round S^4 configuration satisfies Regge EOM
# dS/dl_e should be proportional to dV/dl_e (cosmological constant term)
# For pure gravity without Lambda: dS/dl_e = 0 is Ricci-flat (not S^4)
# For S^4 (constant positive curvature): dS/dl_e = Lambda * dV/dl_e

# Compute gradient by FD: perturb one edge of each type
delta = 1e-6

# Equatorial edge gradient
e_eq = edges[0]  # Pick first equatorial edge
vi, vj = e_eq
direction = verts_5d[vj] - verts_5d[vi]
direction = direction / np.linalg.norm(direction)

verts_plus = verts_5d.copy()
verts_plus[vj] = verts_5d[vj] + delta * direction
S_plus = compute_S_regge(verts_plus, four_simplices, tri_to_simp)

verts_minus = verts_5d.copy()
verts_minus[vj] = verts_5d[vj] - delta * direction
S_minus = compute_S_regge(verts_minus, four_simplices, tri_to_simp)

grad_eq = (S_plus - S_minus) / (2*delta)
print(f"  dS/dl (equatorial edge): {grad_eq:.6f}")

# N-radial edge gradient
vi_rad = 0  # equatorial vertex
direction_N = verts_5d[vi_rad] - verts_5d[idx_N]
direction_N = direction_N / np.linalg.norm(direction_N)

verts_plus = verts_5d.copy()
verts_plus[vi_rad] = verts_5d[vi_rad] + delta * direction_N
S_plus = compute_S_regge(verts_plus, four_simplices, tri_to_simp)

verts_minus = verts_5d.copy()
verts_minus[vi_rad] = verts_5d[vi_rad] - delta * direction_N
S_minus = compute_S_regge(verts_minus, four_simplices, tri_to_simp)

grad_Nrad = (S_plus - S_minus) / (2*delta)
print(f"  dS/dl (N-radial edge):   {grad_Nrad:.6f}")

# Extract Lambda from EOM: dS_geom/dl = 2*Lambda * dV/dl
# (Regge form of R_ab = Lambda*g_ab on constant curvature space)
print(f"\n  dV/dl (equatorial): {dVdl_eq:.6f}")
print(f"  dV/dl (N-radial):   {dVdl_Nr:.6f}")

if abs(dVdl_eq) > 1e-10 and abs(dVdl_Nr) > 1e-10:
    Lambda_from_eq = grad_eq / (2 * dVdl_eq)
    Lambda_from_Nr = grad_Nrad / (2 * dVdl_Nr)
    print(f"\n  Lambda (from equatorial): {Lambda_from_eq:.4f}")
    print(f"  Lambda (from N-radial):   {Lambda_from_Nr:.4f}")
    print(f"  Continuum Lambda(S^4, r=1) = 3.0000")
    print(f"  Ratio eq:  {Lambda_from_eq/3:.4f}")
    print(f"  Ratio rad: {Lambda_from_Nr/3:.4f}")
    Lambda_avg = (Lambda_from_eq + Lambda_from_Nr) / 2
    print(f"  Average: {Lambda_avg:.4f}")

    # G from dimensional analysis on S^4:
    # Lambda = 3/r^2, and R = 4*Lambda = 12/r^2
    # With our Lambda_lattice and lattice spacing a = 1/phi:
    # Lambda_phys = Lambda_lattice / a^2 = Lambda_lattice * phi^2
    print(f"\n  Lambda_lattice * phi^2 = {Lambda_avg * PHI**2:.4f}")
    print(f"  (if a=1/phi is Planck length, this gives Lambda in Planck units)")

print(f"\n  Schlaefli identity: dS/dl_e = sum_tri (dA/dl_e) * epsilon_tri")
print(f"  S^4 is NOT Ricci-flat => gradient nonzero (needs Lambda)")


# ============================================================
# PART 7: FRAMEWORK CONNECTIONS
# ============================================================
print("\n" + "=" * 72)
print("PART 7: FRAMEWORK CONNECTIONS")
print("=" * 72)

print(f"\n  Key quantities:")
print(f"    theta_A (radial):     {theta_A:.6f} rad = {np.degrees(theta_A):.4f} deg")
print(f"    theta_B (equatorial): {theta_B:.6f} rad = {np.degrees(theta_B):.4f} deg")
print(f"    eps_eq:  {eps_eq:.6f} rad")
print(f"    eps_rad: {eps_rad:.6f} rad")
print(f"    S_Regge: {S_regge:.6f}")
print(f"    S_cont:  {S_continuum:.6f}")

# Test framework expressions
tests = {
    'pi/3': np.pi/3,
    'pi/4': np.pi/4,
    'pi/5': np.pi/5,
    'arccos(1/4)': np.arccos(0.25),
    'arccos(1/3)': np.arccos(1/3),
    'arccos(phi/3)': np.arccos(PHI/3),
    'arccos(1/phi^2)': np.arccos(1/PHI**2),
    'arccos(phi-1)': np.arccos(PHI-1),
    'pi/phi^2': np.pi/PHI**2,
    '2*pi/5': 2*np.pi/5,
}

print(f"\n  theta_A = {theta_A:.6f}. Close to:")
for name, val in sorted(tests.items(), key=lambda x: abs(x[1]-theta_A)):
    if abs(val - theta_A) < 0.1:
        print(f"    {name} = {val:.6f} (diff = {theta_A - val:.6f})")

print(f"\n  theta_B = {theta_B:.6f}. Close to:")
for name, val in sorted(tests.items(), key=lambda x: abs(x[1]-theta_B)):
    if abs(val - theta_B) < 0.1:
        print(f"    {name} = {val:.6f} (diff = {theta_B - val:.6f})")

# Check eps values against framework
print(f"\n  eps_eq = {eps_eq:.6f}. As fraction of 2*pi: {eps_eq/(2*np.pi):.6f}")
print(f"  eps_rad = {eps_rad:.6f}. As fraction of 2*pi: {eps_rad/(2*np.pi):.6f}")

# S_Regge in terms of framework quantities
print(f"\n  S_Regge / (2*pi)   = {S_regge/(2*np.pi):.4f}")
print(f"  S_Regge / (4*pi^2) = {S_regge/(4*np.pi**2):.4f}")
print(f"  S_Regge / (8*pi^2) = {S_regge/(8*np.pi**2):.4f}")
print(f"  S_Regge / N        = {S_regge/120:.4f}")
print(f"  S_Regge / phi      = {S_regge/PHI:.4f}")

# G_4D from various normalizations
if abs(S_regge) > 1e-10:
    G_from_Regge = 2*np.pi / abs(S_regge)
    print(f"\n  If S_Regge = 2*pi/G (on-shell S^4):")
    print(f"    G = 2*pi/S_Regge = {G_from_Regge:.6f}")
    print(f"    1/G = {1/G_from_Regge:.4f}")
    print(f"    Compare: 12 (degree), {PHI**2:.4f} (phi^2), {1/ALPHA:.4f} (1/alpha)")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"""
4D REGGE CALCULUS ON S^4 (DOUBLE-CONE OVER 600-CELL)
=====================================================

TRIANGULATION:
  S^4 = double-cone: 122 vertices, {n_edges_4d} edges, {n_tris_4d} triangles,
  {n_tets_4d} tetrahedra, 1200 four-simplices
  Euler characteristic: {chi} (S^4 = 2)
  ALL 4-simplices congruent: 4 edges sqrt(2), 6 edges 1/phi

4D VOLUME:
  V_simplex = {V_simplex:.6f}, V_total = {V_total:.6f}
  V_cont = 8*pi^2/3 = {V_continuum:.4f}, ratio = {V_total/V_continuum:.4f}

DIHEDRAL ANGLES (universal):
  theta_A (radial triangle):     {np.degrees(theta_A):.4f} deg = {theta_A:.6f} rad
  theta_B (equatorial triangle): {np.degrees(theta_B):.4f} deg = {theta_B:.6f} rad

DEFICIT ANGLES:
  eps_eq  = 2*pi - 4*theta_B = {np.degrees(eps_eq):.4f} deg = {eps_eq:.6f} rad
  eps_rad = 2*pi - 5*theta_A = {np.degrees(eps_rad):.4f} deg = {eps_rad:.6f} rad

REGGE ACTION:
  S_Regge = {S_regge:.4f}
  S_continuum (round S^4) = 32*pi^2 = {S_continuum:.4f}
  Ratio: {S_regge/S_continuum:.6f}

GRADIENT:
  dS/dl_eq  = {grad_eq:.6f}
  dS/dl_rad = {grad_Nrad:.6f}
""")

t_end = time.time()
print(f"TOTAL TIME: {t_end-t0:.1f}s")

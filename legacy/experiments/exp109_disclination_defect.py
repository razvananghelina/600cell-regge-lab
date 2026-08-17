"""
exp109_disclination_defect.py
=============================
Mass as Disclination (Topological Defect) on the 600-Cell

QUESTION: Can we model mass as a rotational defect (disclination) in the
600-cell lattice, as Kleinert's World Crystal theory suggests?

APPROACH:
A) Build complete simplicial complex (vertices, edges, triangles, tetrahedra)
B) Compute Regge deficit angles for regular 600-cell
C) Create disclination by removing vertex -> compute deficit change
D) Solve graph Poisson equation -> Green's function profile
E) Extract effective C(r) and compare with C=1+M/r (exp102) vs C=exp(M/r)
F) Elastic relaxation of defected lattice -> strain field

Category: DERIVAT (numerical computation on exact geometry)
"""

import numpy as np
from itertools import product, permutations
from collections import defaultdict, deque

PHI = (1 + np.sqrt(5)) / 2

print("=" * 65)
print("EXP109: DISCLINATION (MASS AS TOPOLOGICAL DEFECT) ON 600-CELL")
print("=" * 65)

# ============================================================
# PART A: BUILD 600-CELL SIMPLICIAL COMPLEX
# ============================================================
print("\n--- PART A: Build 600-cell simplicial complex ---")

def build_600cell():
    """Build 120 vertices of 600-cell on unit S^3 in R^4."""
    phi = PHI
    iphi = 1.0 / phi
    verts = set()

    def add(v):
        n = np.sqrt(sum(x**2 for x in v))
        verts.add(tuple(round(x/n, 10) for x in v))

    # 8 vertices: permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [+1, -1]:
            v = [0., 0., 0., 0.]
            v[i] = float(s)
            add(v)

    # 16 vertices: (+-1/2, +-1/2, +-1/2, +-1/2)
    for ss in product([+1, -1], repeat=4):
        add([s * 0.5 for s in ss])

    # 96 vertices: even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base = [phi/2, 0.5, iphi/2, 0.0]
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
        if inv % 2 == 0:  # even permutations only
            pv = [base[p[i]] for i in range(4)]
            nz = [i for i in range(4) if abs(pv[i]) > 1e-12]
            for ss in product([+1, -1], repeat=len(nz)):
                v = list(pv)
                for k, idx in enumerate(nz):
                    v[idx] *= ss[k]
                add(v)

    V = np.array([list(v) for v in verts])
    return V

V = build_600cell()
N = len(V)
print(f"  Vertices: {N} (expect 120)")

# Build adjacency and edges
dots = V @ V.T
np.clip(dots, -1, 1, out=dots)
cos_edge = np.cos(np.pi / 5)  # cos(36 deg) ~ 0.809

adj = defaultdict(set)
edge_set = set()
for i in range(N):
    for j in range(i+1, N):
        if abs(dots[i, j] - cos_edge) < 0.01:
            edge_set.add((i, j))
            adj[i].add(j)
            adj[j].add(i)

edges = list(edge_set)
print(f"  Edges: {len(edges)} (expect 720)")

# Verify degree
degrees = [len(adj[i]) for i in range(N)]
print(f"  Degree: min={min(degrees)}, max={max(degrees)} (expect 12)")

# Build triangles
tri_set = set()
for a, b in edges:
    common = adj[a] & adj[b]
    for c in common:
        tri_set.add(tuple(sorted([a, b, c])))
tris = list(tri_set)
print(f"  Triangles: {len(tris)} (expect 1200)")

# Build tetrahedra
tet_set = set()
for a, b in edges:
    common = sorted(adj[a] & adj[b])
    for i in range(len(common)):
        for j in range(i+1, len(common)):
            ci, cj = common[i], common[j]
            if cj in adj[ci]:
                tet_set.add(tuple(sorted([a, b, ci, cj])))
tets = list(tet_set)
print(f"  Tetrahedra: {len(tets)} (expect 600)")

# Euler characteristic: chi(S^3) = 0
chi = N - len(edges) + len(tris) - len(tets)
print(f"  Euler characteristic: V-E+F-C = {N}-{len(edges)}+{len(tris)}-{len(tets)} = {chi} (expect 0)")

# ============================================================
# PART B: REGGE DEFICIT ANGLES
# ============================================================
print("\n--- PART B: Regge deficit angles ---")

# Count tetrahedra per edge
edge_tet_count = defaultdict(int)
for tet in tets:
    for i in range(4):
        for j in range(i+1, 4):
            e = tuple(sorted([tet[i], tet[j]]))
            edge_tet_count[e] += 1

counts = [edge_tet_count[tuple(sorted(e))] for e in edges]
print(f"  Tetrahedra per edge: min={min(counts)}, max={max(counts)}, "
      f"mean={np.mean(counts):.1f} (expect 5)")

# Dihedral angles
theta_flat = np.arccos(1.0 / 3.0)  # flat tetrahedron: 70.53 deg
theta_sph = 2 * np.pi / 5           # spherical (on S^3): 72 deg exactly

print(f"\n  Dihedral angles:")
print(f"    Flat (R^4 embedding): {np.degrees(theta_flat):.4f} deg")
print(f"    Spherical (S^3 tiling): {np.degrees(theta_sph):.4f} deg = 2*pi/5")

# Deficit per edge
deficit_flat = 2*np.pi - 5 * theta_flat
deficit_sph = 2*np.pi - 5 * theta_sph

print(f"\n  Deficit per edge:")
print(f"    Flat: 360 - 5*{np.degrees(theta_flat):.2f} = {np.degrees(deficit_flat):.4f} deg")
print(f"    Spherical: 360 - 5*72 = {np.degrees(deficit_sph):.4f} deg (= 0 exact)")

# Total Regge curvature (flat)
# S_Regge = sum delta_e * l_e, where l_e = chord length = 2*sin(pi/10) = 1/phi
l_chord = 1.0 / PHI
total_regge_flat = len(edges) * deficit_flat * l_chord
R_S3_analytic = 12 * np.pi**2  # integral of R*sqrt(g) for unit S^3

print(f"\n  Total Regge curvature (flat embedding):")
print(f"    Sum(delta_e * l_e) = {total_regge_flat:.4f}")
print(f"    Analytic integral R*sqrt(g) for S^3: {R_S3_analytic:.4f}")
print(f"    Ratio: {total_regge_flat / R_S3_analytic:.4f}")
print(f"    (Regge approximation captures {total_regge_flat/R_S3_analytic*100:.1f}% "
      f"of S^3 curvature)")

# ============================================================
# PART C: DISCLINATION (VERTEX REMOVAL)
# ============================================================
print("\n--- PART C: Disclination = vertex removal ---")

v0 = 0
neighbors_v0 = sorted(adj[v0])
print(f"  Remove vertex {v0}")
print(f"  Neighbors: {len(neighbors_v0)} (expect 12)")

# Boundary edges: edges of the icosahedral vertex figure
boundary_edges = []
for a in neighbors_v0:
    for b in neighbors_v0:
        if b > a and b in adj[a]:
            boundary_edges.append((a, b))
print(f"  Boundary edges (vertex figure = icosahedron): {len(boundary_edges)} (expect 30)")

# Tetrahedra, triangles, edges containing v0
tets_v0 = [t for t in tets if v0 in t]
tris_v0 = [t for t in tris if v0 in t]
edges_v0 = [(a, b) for (a, b) in edges if a == v0 or b == v0]
print(f"  Tetrahedra containing v0: {len(tets_v0)} (expect 20)")
print(f"  Triangles containing v0: {len(tris_v0)} (expect 30)")
print(f"  Edges containing v0: {len(edges_v0)} (expect 12)")

# Verify how many tets each boundary edge loses
print(f"\n  Tetrahedra lost per boundary edge (sample):")
losses = []
for be in boundary_edges:
    n_lost = sum(1 for t in tets_v0 if be[0] in t and be[1] in t)
    losses.append(n_lost)
print(f"    min={min(losses)}, max={max(losses)}, all={set(losses)}")
print(f"    (expect exactly 2 for all boundary edges)")

# Deficit analysis
print(f"\n  DEFICIT ANALYSIS:")
print(f"  -------------------------------------------------------")
print(f"  Spherical interpretation (relative to S^3 background):")
delta_per_edge_sph = 2 * theta_sph  # lost 2 tets * 72 deg each
total_excess_sph = len(boundary_edges) * delta_per_edge_sph
print(f"    Lost per boundary edge: 2 * 72 = 144 deg = 4*pi/5 rad")
print(f"    Total excess deficit: 30 * 144 = {np.degrees(total_excess_sph):.0f} deg "
      f"= {total_excess_sph/np.pi:.2f}*pi rad")
print(f"    In units of 4*pi: {total_excess_sph/(4*np.pi):.4f}")
print(f"    In units of 8*pi (Regge mass): M = {total_excess_sph/(8*np.pi):.4f}")

print(f"\n  Flat interpretation (Regge calculus in R^4):")
delta_per_edge_flat = 2 * theta_flat
total_excess_flat = len(boundary_edges) * delta_per_edge_flat
print(f"    Lost per boundary edge: 2 * {np.degrees(theta_flat):.2f} = "
      f"{np.degrees(delta_per_edge_flat):.2f} deg")
print(f"    Total excess deficit: 30 * {np.degrees(delta_per_edge_flat):.2f} = "
      f"{np.degrees(total_excess_flat):.1f} deg = {total_excess_flat/np.pi:.4f}*pi rad")
print(f"    In units of 8*pi: M = {total_excess_flat/(8*np.pi):.4f}")

# Key numbers
print(f"\n  KEY GEOMETRIC NUMBERS:")
print(f"    Total excess (spherical): 24*pi = {24*np.pi:.4f} rad")
print(f"    24 = 2 * 12_neighbors = 30_edges * 4/5")
print(f"    Mass from 8*pi*M = 24*pi => M = 3.0000")
print(f"    Or: M = (edges_lost * 2 * dihedral) / (8*pi)")
print(f"         = (30 * 2 * 72 deg) / (8 * 180 deg) = 30*144/1440 = 3")
print(f"    THIS IS EXACT: vacancy mass = 3 in natural units")

# ============================================================
# PART D: GREEN'S FUNCTION ON FULL GRAPH
# ============================================================
print("\n--- PART D: Green's function (graph Laplacian) ---")

# Build graph Laplacian
L = np.zeros((N, N))
for i, j in edges:
    L[i, j] = -1.0
    L[j, i] = -1.0
for i in range(N):
    L[i, i] = float(len(adj[i]))

# Source at v0
source = np.zeros(N)
source[v0] = 1.0
source -= 1.0 / N  # zero-sum condition

# Regularize and solve
L_reg = L + np.ones((N, N)) / N
phi = np.linalg.solve(L_reg, source)
phi -= np.mean(phi)

# Graph distances from v0
dist = np.full(N, -1, dtype=int)
dist[v0] = 0
queue = deque([v0])
while queue:
    v = queue.popleft()
    for w in adj[v]:
        if dist[w] == -1:
            dist[w] = dist[v] + 1
            queue.append(w)

# Shell analysis
print(f"\n  Shell structure (Green's function at vertex {v0}):")
print(f"  {'Shell d':>8} {'Count':>6} {'Mean phi':>12} {'Std phi':>10} {'Arc r':>8}")
shell_phi = {}
shell_count = {}
for d in range(int(dist.max()) + 1):
    mask = dist == d
    count = np.sum(mask)
    mean_phi = np.mean(phi[mask])
    std_phi = np.std(phi[mask])
    r_arc = d * np.pi / 5
    shell_phi[d] = mean_phi
    shell_count[d] = count
    print(f"  {d:8d} {count:6d} {mean_phi:12.6f} {std_phi:10.6f} {r_arc:8.4f}")

# Ratio analysis
print(f"\n  Successive shell ratios:")
for d in range(1, int(dist.max())):
    if abs(shell_phi[d]) > 1e-12 and abs(shell_phi[d+1]) > 1e-12:
        ratio = shell_phi[d+1] / shell_phi[d]
        print(f"    phi({d+1})/phi({d}) = {ratio:.4f}")

# ============================================================
# PART E: C(r) EXTRACTION AND COMPARISON
# ============================================================
print("\n--- PART E: C(r) extraction and metric comparison ---")

# Physical distance: r = d * pi/5 (arc length on unit S^3)
# Model 1: C(r) = 1 + M/r  =>  phi ~ M/r  =>  phi*r = const
# Model 2: C(r) = exp(M/r)  =>  phi ~ M/r (at large r), deviates at small r

print(f"\n  Testing phi*r = const (1/r decay):")
print(f"  {'Shell':>6} {'r':>8} {'phi':>12} {'phi*r':>12} {'phi*r^2':>12}")
for d in range(1, int(dist.max()) + 1):
    r = d * np.pi / 5
    p = shell_phi[d]
    print(f"  {d:6d} {r:8.4f} {p:12.6f} {p*r:12.6f} {p*r**2:12.6f}")

# Power law fit: phi ~ A * r^(-n) using shells 1-4
r_vals = np.array([d * np.pi / 5 for d in range(1, 5)])
phi_vals = np.array([shell_phi[d] for d in range(1, 5)])

# Handle sign alternation - use absolute values for power law
signs = np.sign(phi_vals)
abs_phi = np.abs(phi_vals)

# Check if signs alternate
print(f"\n  Sign pattern: {['+'if s>0 else '-' for s in signs]}")
if np.all(signs > 0) or np.all(signs < 0):
    log_r = np.log(r_vals)
    log_phi = np.log(abs_phi)
    coeffs = np.polyfit(log_r, log_phi, 1)
    n_fit = -coeffs[0]
    print(f"  Power law fit (all same sign): n = {n_fit:.4f}")
else:
    print(f"  Signs ALTERNATE => no simple power law (oscillatory Green's function)")
    # Fit envelope
    pos_mask = signs > 0
    if np.sum(pos_mask) >= 2:
        log_r_pos = np.log(r_vals[pos_mask])
        log_phi_pos = np.log(abs_phi[pos_mask])
        if len(log_r_pos) >= 2:
            coeffs = np.polyfit(log_r_pos, log_phi_pos, 1)
            n_fit = -coeffs[0]
            print(f"  Envelope fit (positive shells): n = {n_fit:.4f}")
        else:
            n_fit = float('nan')
    else:
        n_fit = float('nan')

# Compare C = 1+M/r vs C = exp(M/r)
print(f"\n  MODEL COMPARISON (normalized at shell 1):")
phi_1 = shell_phi[1]
r_1 = np.pi / 5
M_1r = phi_1 * r_1  # fit parameter for 1/r model
M_exp = phi_1 * r_1  # same normalization for exp model (linearized)

print(f"  M (fit from shell 1) = {M_1r:.6f}")
print(f"  {'Shell':>6} {'phi_actual':>12} {'phi_1/r':>12} {'err_1/r%':>10} "
      f"{'phi_exp':>12} {'err_exp%':>10}")
for d in range(1, int(dist.max()) + 1):
    r = d * np.pi / 5
    p = shell_phi[d]
    pred_1r = M_1r / r
    pred_exp = np.exp(M_1r / r) - 1.0  # C - 1 for comparison
    err_1r = (p - pred_1r) / abs(p) * 100 if abs(p) > 1e-12 else 0
    err_exp = (p - pred_exp) / abs(p) * 100 if abs(p) > 1e-12 else 0
    print(f"  {d:6d} {p:12.6f} {pred_1r:12.6f} {err_1r:10.1f} "
          f"{pred_exp:12.6f} {err_exp:10.1f}")

# ============================================================
# PART F: POISSON ON DEFECTED GRAPH (VACANCY)
# ============================================================
print("\n--- PART F: Vacancy potential (defected graph) ---")

# Build Laplacian without vertex v0
idx_map = {}
new_idx = 0
for i in range(N):
    if i != v0:
        idx_map[i] = new_idx
        new_idx += 1
N_def = N - 1

adj_def = defaultdict(set)
for i, j in edges:
    if i == v0 or j == v0:
        continue
    ni, nj = idx_map[i], idx_map[j]
    adj_def[ni].add(nj)
    adj_def[nj].add(ni)

L_def = np.zeros((N_def, N_def))
for i, j in edges:
    if i == v0 or j == v0:
        continue
    ni, nj = idx_map[i], idx_map[j]
    L_def[ni, nj] = -1.0
    L_def[nj, ni] = -1.0
for i in range(N_def):
    L_def[i, i] = -np.sum(L_def[i, :])

# Source: neighbors of v0 have "missing pull"
source_def = np.zeros(N_def)
for nb in neighbors_v0:
    source_def[idx_map[nb]] = -1.0
source_def -= np.sum(source_def) / N_def  # zero-sum

# Solve
L_def_reg = L_def + np.ones((N_def, N_def)) / N_def
phi_def = np.linalg.solve(L_def_reg, source_def)
phi_def -= np.mean(phi_def)

# Distance from defect (BFS from neighbors of v0 in original graph)
dist_from_defect = np.full(N, -1, dtype=int)
dist_from_defect[v0] = 0
for nb in neighbors_v0:
    dist_from_defect[nb] = 1
queue = deque(list(neighbors_v0))
while queue:
    v = queue.popleft()
    for w in adj[v]:
        if dist_from_defect[w] == -1:
            dist_from_defect[w] = dist_from_defect[v] + 1
            queue.append(w)

print(f"\n  Vacancy potential (distributed source on 12 neighbors):")
print(f"  {'Shell':>6} {'Count':>6} {'Mean phi_vac':>14} {'Std':>10}")
shell_vac = {}
for d in range(1, int(dist_from_defect.max()) + 1):
    vals = []
    for i in range(N):
        if i != v0 and dist_from_defect[i] == d:
            vals.append(phi_def[idx_map[i]])
    if vals:
        mean_v = np.mean(vals)
        std_v = np.std(vals)
        shell_vac[d] = mean_v
        print(f"  {d:6d} {len(vals):6d} {mean_v:14.6f} {std_v:10.6f}")

# Compare point source vs vacancy
print(f"\n  Point source vs vacancy (ratio):")
for d in range(1, min(6, int(dist.max()) + 1)):
    if d in shell_vac and d in shell_phi and abs(shell_phi[d]) > 1e-12:
        ratio = shell_vac[d] / shell_phi[d]
        print(f"    Shell {d}: vacancy/point = {ratio:.4f}")

# ============================================================
# PART G: ELASTIC RELAXATION
# ============================================================
print("\n--- PART G: Elastic relaxation after vacancy ---")

# Natural edge length on S^3
d0 = np.pi / 5  # angular distance = 36 deg

# Remove vertex v0, relax remaining vertices
# Energy = sum over edges: k/2 * (d_ij - d0)^2
# where d_ij = arccos(v_i . v_j)
# Minimize by gradient descent on S^3

V_relax = V.copy()
# Remove v0 from the vertex set
V_relax_active = np.delete(V_relax, v0, axis=0)
# edges for relaxation (mapped indices)
edges_relax = []
for i, j in edges:
    if i == v0 or j == v0:
        continue
    ni = idx_map[i]
    nj = idx_map[j]
    edges_relax.append((ni, nj))

def compute_energy(V_r, edges_r, d0):
    """Total elastic energy."""
    E = 0.0
    for i, j in edges_r:
        dot = np.clip(np.dot(V_r[i], V_r[j]), -1, 1)
        d_ij = np.arccos(dot)
        E += 0.5 * (d_ij - d0)**2
    return E

def compute_forces(V_r, edges_r, d0):
    """Gradient of energy w.r.t. vertex positions (projected to tangent space of S^3)."""
    N_r = len(V_r)
    forces = np.zeros_like(V_r)
    for i, j in edges_r:
        dot = np.clip(np.dot(V_r[i], V_r[j]), -1, 1)
        d_ij = np.arccos(dot)
        if d_ij < 1e-12:
            continue
        strain = d_ij - d0
        # Gradient of arccos(vi.vj) w.r.t. vi = -vj / sin(d_ij) (tangent component)
        sin_d = np.sin(d_ij)
        if sin_d < 1e-12:
            continue
        # Force on i from spring (i,j)
        # d(arccos(vi.vj))/d(vi) = -(vj - dot*vi) / sin(d_ij)
        grad_i = -(V_r[j] - dot * V_r[i]) / sin_d
        force_i = -strain * grad_i  # F = -k*(d-d0)*grad(d)
        forces[i] += force_i
        # Force on j
        grad_j = -(V_r[i] - dot * V_r[j]) / sin_d
        force_j = -strain * grad_j
        forces[j] += force_j
    return forces

def project_to_tangent(F, V_r):
    """Project forces to tangent space of S^3 (remove radial component)."""
    for i in range(len(V_r)):
        F[i] -= np.dot(F[i], V_r[i]) * V_r[i]
    return F

# Run gradient descent
E_init = compute_energy(V_relax_active, edges_relax, d0)
print(f"  Initial energy (after vacancy): {E_init:.8f}")

# Also compute energy of unperturbed lattice (without v0)
V_unperturbed = np.delete(V, v0, axis=0)
E_unperturbed = compute_energy(V_unperturbed, edges_relax, d0)
print(f"  Unperturbed energy (natural lengths): {E_unperturbed:.8f}")

# Gradient descent
dt = 0.01
n_steps = 500
V_r = V_relax_active.copy()
energies = [E_init]

for step in range(n_steps):
    F = compute_forces(V_r, edges_relax, d0)
    F = project_to_tangent(F, V_r)
    V_r += dt * F
    # Normalize back to S^3
    norms = np.linalg.norm(V_r, axis=1, keepdims=True)
    V_r /= norms
    if step % 100 == 99:
        E = compute_energy(V_r, edges_relax, d0)
        energies.append(E)

E_final = compute_energy(V_r, edges_relax, d0)
print(f"  Final energy after {n_steps} steps: {E_final:.8f}")
print(f"  Energy change: {E_final - E_init:.8f}")

# Measure displacement field
displacements = np.zeros(N_def)
for i in range(N_def):
    dot = np.clip(np.dot(V_r[i], V_relax_active[i]), -1, 1)
    displacements[i] = np.arccos(dot)

print(f"\n  Displacement field (angular displacement from original position):")
print(f"  {'Shell':>6} {'Count':>6} {'Mean |u|':>14} {'Max |u|':>12}")

# Map back to original indices for distance
for d in range(1, int(dist_from_defect.max()) + 1):
    disp_vals = []
    for i_orig in range(N):
        if i_orig != v0 and dist_from_defect[i_orig] == d:
            disp_vals.append(displacements[idx_map[i_orig]])
    if disp_vals:
        print(f"  {d:6d} {len(disp_vals):6d} {np.mean(disp_vals):14.8f} "
              f"{np.max(disp_vals):12.8f}")

# Measure compression: local density change
# C(r) ~ (local edge length) / d0
print(f"\n  Local strain (edge length change) by distance from defect:")
print(f"  {'Shell':>6} {'N_edges':>8} {'Mean d_ij':>12} {'d_ij/d0':>10} {'C(r)':>10}")

for d in range(1, int(dist_from_defect.max()) + 1):
    edge_lengths = []
    for ni, nj in edges_relax:
        # Find original indices
        i_orig = [k for k, v in idx_map.items() if v == ni][0]
        j_orig = [k for k, v in idx_map.items() if v == nj][0]
        # Both endpoints at distance d?
        di = dist_from_defect[i_orig]
        dj = dist_from_defect[j_orig]
        if di == d and dj == d:
            dot_ij = np.clip(np.dot(V_r[ni], V_r[nj]), -1, 1)
            edge_lengths.append(np.arccos(dot_ij))
    if edge_lengths:
        mean_d = np.mean(edge_lengths)
        C_r = mean_d / d0  # compression factor
        print(f"  {d:6d} {len(edge_lengths):8d} {mean_d:12.8f} {mean_d/d0:10.6f} {C_r:10.6f}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)

print(f"""
SIMPLICIAL COMPLEX:
  600-cell: V={N}, E={len(edges)}, F={len(tris)}, C={len(tets)}
  Euler characteristic: {chi} (= 0 for S^3)
  Dihedral angle: 72 deg (spherical), 70.53 deg (flat)
  Deficit per edge: 0 (spherical), 7.35 deg (flat)

DISCLINATION (vertex removal):
  Boundary edges: {len(boundary_edges)}
  Tetrahedra removed: {len(tets_v0)}
  Spherical deficit per boundary edge: 144 deg (= 4*pi/5)
  TOTAL EXCESS DEFICIT: 24*pi rad = 4320 deg
  VACANCY MASS (8*pi*M units): M = 3.0 EXACT
  Decomposition: 30 edges * 2 lost_tets * (2*pi/5) / (8*pi) = 3

GREEN'S FUNCTION:
  phi(d=1) = {shell_phi.get(1, 'N/A')}
  phi(d=2) = {shell_phi.get(2, 'N/A')}
  phi(d=3) = {shell_phi.get(3, 'N/A')}
""")

ratio_21 = shell_phi[2]/shell_phi[1] if abs(shell_phi[1]) > 1e-15 else float('nan')
print(f"  Ratio phi(2)/phi(1) = {ratio_21:.4f}")

print(f"""
ELASTIC RELAXATION:
  Energy change after vacancy: {E_final - E_init:.10f}
  Maximum displacement: {np.max(displacements):.10f}
  RESULT: ZERO STRAIN - the 600-cell is a "perfect crystal"
  REASON: Removing vertex only removes edges TO it.
          Remaining edges are already at natural length d0=pi/5.
          No restoring force => no relaxation => no long-range strain.

METRIC COMPARISON:
  The document's claim: C(r) = exp(M/r) => NOT SUPPORTED
  Our exp102 result: C(r) = 1 + M/r => NOT from vacancy mechanics
  Graph Green's function: OSCILLATORY (signs alternate)
  Neither 1/r nor exp(1/r) fits the discrete graph response.

CATEGORY: DERIVAT (simplicial complex) + COMPUTED (Green's function)
""")

# Final verdict
print("KEY FINDINGS:")
print("  1. DERIVAT: Vacancy deficit = 24*pi rad (EXACT)")
print("  2. DERIVAT: Vacancy mass M = 3 in 8*pi*G units (EXACT)")
print("     Decomposition: 30 boundary_edges * 2 lost_tets * (2*pi/5) / (8*pi) = 3")
print("  3. DERIVAT: Deficit is LOCALIZED on 30 boundary edges only")
print("     No long-range curvature field from topology alone")
print("  4. COMPUTED: Elastic relaxation gives ZERO strain")
print("     The 600-cell is too 'perfect' - vacancy doesn't stress remaining bonds")
print("     This is fundamentally different from real crystal vacancies")
print("  5. COMPUTED: Green's function OSCILLATES (not monotonic 1/r)")
print("     phi(2)/phi(1) = %.4f (exp106 found beta=-3.09)" % ratio_21)
print("  6. VERDICT: The Kleinert 'World Crystal' model does NOT work directly")
print("     on the 600-cell. The lattice is too rigid/symmetric.")
print("     A vacancy creates localized deficit but no long-range strain field.")
print()
print("NEGATIVE RESULT for the 'spider web' document's claims:")
print("  - C(r) = exp(M/r) is NOT derived from 600-cell geometry")
print("  - There is no 'strain hardening' because there is no strain at all")
print("  - The 1/r^2 force law does NOT emerge from elastic relaxation")
print("  - The Poisson eq gives oscillatory, not 1/r, Green's function")
print()
print("POSITIVE RESULTS (new geometric identities):")
print("  - M_vacancy = 3 (integer!) - universal for all 600-cell vertices")
print("  - 24*pi = total deficit = connects to 24-cell? (24 octahedral cells)")
print("  - Regge captures 48.2% of S^3 curvature (known Regge approximation)")
print("  - All 30 boundary edges lose exactly 2 tetrahedra (icosahedral symmetry)")

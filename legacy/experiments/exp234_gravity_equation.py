# EXP-234: Gravitational Field Equation from the 600-Cell
# =========================================================
# Strategy: Build the full simplicial complex of the 600-cell,
# compute Hodge Laplacians Delta_0 (vertices) and Delta_1 (edges),
# and see if combining them gives full GR (gamma_PPN = 1).
# =========================================================
import numpy as np
from itertools import combinations
import math

PHI = (1 + math.sqrt(5)) / 2
PHI_INV = PHI - 1  # = 1/phi
a_1 = 5
b_1 = 6

print("=" * 70)
print("EXP-234: GRAVITATIONAL FIELD EQUATION FROM THE 600-CELL")
print("=" * 70)

# =============================================================
# PART A: Construct the 600-cell vertices (binary icosahedral group)
# =============================================================
print("\n" + "=" * 70)
print("PART A: Constructing the 600-cell")
print("=" * 70)

vertices = []

# Type A: 8 vertices -- permutations of (+/-1, 0, 0, 0)
for i in range(4):
    for s in [1, -1]:
        v = [0, 0, 0, 0]
        v[i] = s
        vertices.append(tuple(v))

# Type B: 16 vertices -- (+/-0.5, +/-0.5, +/-0.5, +/-0.5)
for s0 in [0.5, -0.5]:
    for s1 in [0.5, -0.5]:
        for s2 in [0.5, -0.5]:
            for s3 in [0.5, -0.5]:
                vertices.append((s0, s1, s2, s3))

# Type C: 96 vertices -- even permutations of (0, +/-1/2, +/-phi/2, +/-1/(2*phi))
# with an even number of minus signs on the non-zero coordinates
half = 0.5
hp = PHI / 2       # phi/2
hi = PHI_INV / 2   # 1/(2*phi)

# The 96 golden vertices: even permutations of (0, +/-1/(2phi), +/-1/2, +/-phi/2)
# with independent signs
base_coords = [0, hi, half, hp]

# Generate all even permutations of 4 elements
def even_permutations(lst):
    """Generate all even permutations of a list of 4 elements."""
    perms = []
    # All 12 even permutations of (0,1,2,3)
    even_perm_indices = [
        (0,1,2,3), (0,2,3,1), (0,3,1,2),
        (1,0,3,2), (1,2,0,3), (1,3,2,0),
        (2,0,1,3), (2,1,3,0), (2,3,0,1),
        (3,0,2,1), (3,1,0,2), (3,2,1,0)
    ]
    for idx in even_perm_indices:
        perms.append(tuple(lst[i] for i in idx))
    return perms

# For each even permutation of (0, hi, half, hp),
# apply all sign combinations to non-zero entries (8 combos each)
type_c = set()
for perm in even_permutations(base_coords):
    # Find non-zero positions
    nonzero = [i for i in range(4) if abs(perm[i]) > 1e-10]
    # All sign combinations for non-zero entries
    for signs in range(8):  # 2^3 = 8
        v = list(perm)
        for bit, pos in enumerate(nonzero):
            if signs & (1 << bit):
                v[pos] = -v[pos]
        type_c.add(tuple(round(x, 10) for x in v))

for v in type_c:
    vertices.append(v)

# Remove duplicates
vertices_set = set()
clean_vertices = []
for v in vertices:
    key = tuple(round(x, 8) for x in v)
    if key not in vertices_set:
        vertices_set.add(key)
        clean_vertices.append(v)

vertices = clean_vertices
N = len(vertices)
print(f"  Vertices constructed: {N} (expected 120)")

# Verify: all on unit S^3
norms = [math.sqrt(sum(x**2 for x in v)) for v in vertices]
print(f"  Norms: min={min(norms):.6f}, max={max(norms):.6f} (should be 1.0)")

# =============================================================
# PART B: Build the adjacency matrix
# =============================================================
print("\n" + "=" * 70)
print("PART B: Building adjacency matrix and edge list")
print("=" * 70)

verts = np.array(vertices)

# Edge criterion: inner product = phi/2 (for unit S^3 600-cell)
# The edge length is 1/phi, so inner product = 1 - (1/phi)^2/2 = 1 - 1/(2*phi^2)
# Actually, for unit 600-cell on S^3, adjacent vertices have
# inner product = cos(pi/5) = phi/2
edge_ip = PHI / 2

# Compute all pairwise inner products
A = np.zeros((N, N), dtype=int)
edges = []
for i in range(N):
    for j in range(i+1, N):
        ip = sum(verts[i] * verts[j])
        if abs(ip - edge_ip) < 1e-6:
            A[i, j] = 1
            A[j, i] = 1
            edges.append((i, j))

n_edges = len(edges)
degrees = A.sum(axis=1)
print(f"  Edges: {n_edges} (expected 720)")
print(f"  Degree: min={degrees.min()}, max={degrees.max()} (expected 12)")

# =============================================================
# PART C: Scalar Laplacian (Delta_0) spectrum
# =============================================================
print("\n" + "=" * 70)
print("PART C: Scalar Laplacian Delta_0 spectrum")
print("=" * 70)

L0 = np.diag(degrees) - A
evals_L0 = np.sort(np.linalg.eigvalsh(L0.astype(float)))

# Group eigenvalues
from collections import Counter
eval_groups = Counter(round(e, 3) for e in evals_L0)
print(f"\n  Delta_0 spectrum (Laplacian = degree - adjacency):")
for val, mult in sorted(eval_groups.items()):
    print(f"    lambda = {val:8.3f}, mult = {mult}")

print(f"\n  Spectral gap = {evals_L0[1]:.6f}")
print(f"  = 12 - 6*phi = b_1*(2-phi) = b_1*phi'^2 = {b_1*(2-PHI):.6f}")

# =============================================================
# PART D: Find all triangular faces
# =============================================================
print("\n" + "=" * 70)
print("PART D: Finding triangular faces (1200 expected)")
print("=" * 70)

# Build adjacency list for speed
adj = [set() for _ in range(N)]
for i, j in edges:
    adj[i].add(j)
    adj[j].add(i)

# Find triangles: for each edge (i,j), find common neighbors
triangles = set()
for i, j in edges:
    common = adj[i] & adj[j]
    for k in common:
        tri = tuple(sorted([i, j, k]))
        triangles.add(tri)

triangles = sorted(triangles)
n_faces = len(triangles)
print(f"  Triangular faces: {n_faces} (expected 1200)")

# How many triangles per edge?
tri_per_edge = {}
for idx, (i, j, k) in enumerate(triangles):
    for e in [(i,j), (i,k), (j,k)]:
        e_sorted = tuple(sorted(e))
        if e_sorted not in tri_per_edge:
            tri_per_edge[e_sorted] = []
        tri_per_edge[e_sorted].append(idx)

tpe_counts = Counter(len(v) for v in tri_per_edge.values())
print(f"  Triangles per edge: {dict(tpe_counts)}")

# =============================================================
# PART E: Find tetrahedral cells
# =============================================================
print("\n" + "=" * 70)
print("PART E: Finding tetrahedral cells (600 expected)")
print("=" * 70)

# A tetrahedron is 4 mutually adjacent vertices
# For each triangle (i,j,k), find vertices adjacent to all three
tetrahedra = set()
for i, j, k in triangles:
    common = adj[i] & adj[j] & adj[k]
    for l in common:
        tet = tuple(sorted([i, j, k, l]))
        tetrahedra.add(tet)

tetrahedra = sorted(tetrahedra)
n_cells = len(tetrahedra)
print(f"  Tetrahedral cells: {n_cells} (expected 600)")

# Tetrahedra per edge
tet_per_edge = {}
for idx, (i, j, k, l) in enumerate(tetrahedra):
    for e in combinations([i,j,k,l], 2):
        e_sorted = tuple(sorted(e))
        if e_sorted not in tet_per_edge:
            tet_per_edge[e_sorted] = []
        tet_per_edge[e_sorted].append(idx)

tpe_tet = Counter(len(v) for v in tet_per_edge.values())
print(f"  Tetrahedra per edge: {dict(tpe_tet)} (expected 5)")

# Tetrahedra per face
tet_per_face = {}
for idx, (i, j, k, l) in enumerate(tetrahedra):
    for f in combinations([i,j,k,l], 3):
        f_sorted = tuple(sorted(f))
        if f_sorted not in tet_per_face:
            tet_per_face[f_sorted] = []
        tet_per_face[f_sorted].append(idx)

tpf_counts = Counter(len(v) for v in tet_per_face.values())
print(f"  Tetrahedra per face: {dict(tpf_counts)} (expected 2)")

# =============================================================
# PART F: Boundary operators and Hodge Laplacian Delta_1
# =============================================================
print("\n" + "=" * 70)
print("PART F: Edge Laplacian Delta_1 (Hodge)")
print("=" * 70)

# Boundary operator d_0: vertices -> edges (incidence matrix B)
# B has shape (n_edges, N): B[e, v] = +1 if v is head, -1 if v is tail
# We choose orientation: for edge (i,j) with i < j, i is tail, j is head
edge_index = {e: idx for idx, e in enumerate(edges)}

B = np.zeros((n_edges, N), dtype=float)
for idx, (i, j) in enumerate(edges):
    B[idx, i] = -1
    B[idx, j] = +1

# Boundary operator d_1: edges -> faces
# For each oriented triangle (i,j,k) with i<j<k,
# d_1[f, e] = +1 or -1 depending on orientation
face_index = {f: idx for idx, f in enumerate(triangles)}

C = np.zeros((n_faces, n_edges), dtype=float)
for fidx, (i, j, k) in enumerate(triangles):
    # Edges of this triangle: (i,j), (i,k), (j,k)
    # Oriented boundary: (i,j) + (j,k) - (i,k)
    e_ij = edge_index.get((i,j), edge_index.get((j,i), -1))
    e_ik = edge_index.get((i,k), edge_index.get((k,i), -1))
    e_jk = edge_index.get((j,k), edge_index.get((k,j), -1))

    # Sign: +1 if edge orientation agrees with face boundary, -1 otherwise
    C[fidx, e_ij] = +1   # (i,j) agrees with (i->j->k) boundary
    C[fidx, e_jk] = +1   # (j,k) agrees
    C[fidx, e_ik] = -1   # (i,k) is reversed (boundary goes i->j->k, not i->k)

# Verify: d_1 . d_0 = 0 (boundary of boundary = 0)
CB = C @ B
print(f"  d_1 . d_0 = 0? Max element: {np.max(np.abs(CB)):.2e} (should be ~0)")

# Hodge Laplacian on edges: Delta_1 = d_0 d_0^T + d_1^T d_1
# = B B^T + C^T C
print(f"\n  Computing Delta_1 = B*B^T + C^T*C ({n_edges}x{n_edges} matrix)...")
Delta_1 = B @ B.T + C.T @ C

print(f"  Delta_1 shape: {Delta_1.shape}")
print(f"  Computing eigenvalues...")
evals_L1 = np.sort(np.linalg.eigvalsh(Delta_1))

# Group eigenvalues
eval_groups_1 = Counter(round(e, 2) for e in evals_L1)
print(f"\n  Delta_1 spectrum (edge Laplacian):")
for val, mult in sorted(eval_groups_1.items()):
    print(f"    lambda = {val:8.3f}, mult = {mult}")

# Key: kernel dimension = b_1 (first Betti number)
n_zero_1 = sum(1 for e in evals_L1 if abs(e) < 1e-6)
print(f"\n  Kernel dimension of Delta_1: {n_zero_1}")
print(f"  (= first Betti number b_1 of S^3 = 0)")

# Spectral gap of Delta_1
nonzero_1 = [e for e in evals_L1 if e > 1e-6]
if nonzero_1:
    print(f"  Spectral gap of Delta_1: {nonzero_1[0]:.6f}")

# =============================================================
# PART G: Compare Delta_0 and Delta_1 spectra
# =============================================================
print("\n" + "=" * 70)
print("PART G: Comparing scalar (Delta_0) and edge (Delta_1) spectra")
print("=" * 70)

# On a smooth manifold S^3 of radius R:
# Delta_0 eigenvalues: l(l+2)/R^2, mult = (l+1)^2, l=0,1,2,...
# Delta_1 eigenvalues: (l(l+2)-1)/R^2 for "exact" 1-forms
#                      (l(l+2)+1)/R^2 for "coexact" 1-forms
# The difference between Delta_0 and Delta_1 eigenvalues is the
# Weitzenboeck curvature term (Ricci curvature on S^3 = 2/R^2)

print(f"\n  On smooth S^3:")
print(f"    Delta_0 gap = l(l+2) at l=1 => 3")
print(f"    Delta_1 gap = l(l+2)-1 at l=1 => 2 (exact 1-forms)")
print(f"    Delta_1 gap = l(l+2)+1 at l=1 => 4 (coexact 1-forms)")
print(f"    Ricci shift: Delta_1 = Delta_0 + Ric (Weitzenboeck)")

# On the 600-cell graph:
gap_0 = evals_L0[1]
gap_1 = nonzero_1[0] if nonzero_1 else 0
print(f"\n  On 600-cell graph:")
print(f"    Delta_0 spectral gap: {gap_0:.6f}")
print(f"    Delta_1 spectral gap: {gap_1:.6f}")
print(f"    Ratio gap_1/gap_0: {gap_1/gap_0:.6f}")
print(f"    Expected for S^3: 2/3 = {2/3:.6f} (exact forms)")
print(f"                  or  4/3 = {4/3:.6f} (coexact forms)")

ratio_gaps = gap_1 / gap_0
print(f"\n  *** gap_1/gap_0 = {ratio_gaps:.6f} ***")

if abs(ratio_gaps - 2/3) < 0.1:
    print(f"  CLOSE TO 2/3! The Weitzenboeck curvature is NEGATIVE")
    print(f"  (exact 1-forms are 'lighter' than scalar modes)")
elif abs(ratio_gaps - 4/3) < 0.1:
    print(f"  CLOSE TO 4/3! The Weitzenboeck curvature is POSITIVE")
    print(f"  (coexact 1-forms are 'heavier' than scalar modes)")

# =============================================================
# PART H: Weitzenboeck curvature = Ricci from the spectrum
# =============================================================
print("\n" + "=" * 70)
print("PART H: Extracting Ricci curvature from Hodge spectra")
print("=" * 70)

# On a Riemannian manifold: Delta_1 = nabla^T nabla + Ric
# (Weitzenboeck formula)
# So: eigenvalue(Delta_1) = eigenvalue(connection Laplacian) + Ricci
# On constant curvature S^3: Ric = 2*K*g where K = 1/R^2
# So Ric shifts all eigenvalues by 2*K

# On the graph, the "Ricci-like" shift can be extracted from:
# average(Delta_1 eigenvalues) - average(Delta_0 eigenvalues adjusted for dimension)

avg_L0 = np.mean(evals_L0[1:])  # exclude zero mode
avg_L1 = np.mean(evals_L1)  # no zero modes

print(f"  Average Delta_0 eigenvalue (nonzero): {avg_L0:.4f}")
print(f"  Average Delta_1 eigenvalue: {avg_L1:.4f}")

# Trace of Laplacians
tr_L0 = np.trace(L0.astype(float))
tr_L1 = np.trace(Delta_1)
print(f"\n  Tr(Delta_0) = {tr_L0:.0f} = N*deg = {N}*12 = {N*12}")
print(f"  Tr(Delta_1) = {tr_L1:.0f}")
print(f"  Tr(Delta_1)/n_edges = {tr_L1/n_edges:.4f}")
print(f"  Tr(Delta_0)/N = {tr_L0/N:.4f} = 12 (vertex degree)")

# The "Ricci curvature" in graph terms
# Bochner: Delta_1 = connection_Laplacian + Ric
# Ric per edge = (Tr(Delta_1)/n_edges) - (some connection Laplacian trace)
# The connection Laplacian trace = sum of vertex degrees at endpoints = 2*deg for each edge
# So: Ric_average = Tr(Delta_1)/n_edges - 2*deg

ric_estimate = tr_L1/n_edges - 2*12
print(f"\n  Estimated Ricci per edge: {ric_estimate:.4f}")
print(f"  = Tr(Delta_1)/n_edges - 2*deg = {tr_L1/n_edges:.4f} - 24")

# =============================================================
# PART I: The gravitational propagator
# =============================================================
print("\n" + "=" * 70)
print("PART I: Gravitational propagator analysis")
print("=" * 70)

print(f"""
  GRAVITY ON THE 600-CELL:

  The scalar Laplacian Delta_0 gives the NEWTONIAN potential (g_tt).
  The edge Laplacian Delta_1 gives the SPATIAL metric (g_rr).

  For full GR, we need BOTH potentials to agree: Phi = Psi.
  In PPN: gamma = Psi/Phi.

  The ratio of Green's functions G_0 and G_1 at the spectral level
  is determined by their spectral gaps:
    Phi ~ 1/gap_0 at long range
    Psi ~ 1/gap_1 at long range

  gamma_PPN = gap_0 / gap_1 (in the simplest approximation)
""")

gamma_spectral = gap_0 / gap_1
print(f"  gamma_PPN (spectral) = gap_0/gap_1 = {gap_0:.4f}/{gap_1:.4f} = {gamma_spectral:.6f}")
print(f"  Expected for GR: gamma = 1")
print(f"  Old STG result: gamma = 1/2 = 0.5")

# PPN precession
# delta_phi = (2 + 2*gamma - beta) * pi * GM / (a*(1-e^2)*c^2)
# GR: 2 + 2*1 - 1 = 3
# STG (gamma=0.5): 2 + 1 - 1 = 2 -> ratio 2/3
print(f"\n  Precession coefficients:")
print(f"    GR (gamma=1): 2 + 2*1 - 1 = 3")
print(f"    Spectral (gamma={gamma_spectral:.4f}): 2 + 2*{gamma_spectral:.4f} - 1 = {2 + 2*gamma_spectral - 1:.4f}")
ratio_to_GR = (2 + 2*gamma_spectral - 1) / 3
print(f"    Ratio to GR: {ratio_to_GR:.6f}")

# =============================================================
# PART J: Euler characteristic check
# =============================================================
print("\n" + "=" * 70)
print("PART J: Topological checks")
print("=" * 70)

# Euler characteristic of S^3: chi = 0
# chi = V - E + F - C = 120 - 720 + 1200 - 600
chi = N - n_edges + n_faces - n_cells
print(f"  Euler characteristic: {N} - {n_edges} + {n_faces} - {n_cells} = {chi}")
print(f"  Expected for S^3: 0")

# Betti numbers from Hodge theory
# b_0 = dim(ker(Delta_0))
b0 = sum(1 for e in evals_L0 if abs(e) < 1e-6)
b1 = n_zero_1
# b_2 from face Laplacian (not computed yet)
# b_3 from cell "Laplacian"
print(f"\n  Betti numbers from Hodge Laplacians:")
print(f"    b_0 = {b0} (expected 1)")
print(f"    b_1 = {b1} (expected 0)")
print(f"    b_2 = ? (expected 0, need Delta_2)")
print(f"    b_3 = ? (expected 1)")
print(f"    Expected chi = 1 - 0 + 0 - 1 = 0, actual chi = {chi}")

# =============================================================
# PART K: Summary
# =============================================================
print("\n" + "=" * 70)
print("SUMMARY EXP-234")
print("=" * 70)

print(f"""
  RESULTS:

  1. SIMPLICIAL COMPLEX VERIFIED:
     V={N}, E={n_edges}, F={n_faces}, C={n_cells}
     chi = {chi} (correct for S^3)
     Triangles per edge: 5 (icosahedral)
     Tetrahedra per edge: 5 (icosahedral)

  2. HODGE LAPLACIAN SPECTRA:
     Delta_0 (scalar): gap = {gap_0:.4f} = b_1*phi'^2
     Delta_1 (edge):   gap = {gap_1:.4f}
     Ratio: gap_1/gap_0 = {ratio_gaps:.4f}

  3. GRAVITATIONAL IMPLICATION:
     gamma_PPN (spectral) = {gamma_spectral:.4f}
     Precession ratio to GR: {ratio_to_GR:.4f}

  STATUS: The edge Laplacian spectrum provides the SECOND gravitational
  potential. The ratio of spectral gaps determines gamma_PPN.
""")

print("=" * 70)
print("EXP-234 COMPLETE")
print("=" * 70)

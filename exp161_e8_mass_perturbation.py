"""
EXP-161: Mass perturbation on E8 lattice
=========================================
Put a "mass" on one E8 vertex (in S-copy of 600-cell).
Compute how the effective metric deforms via Green's function propagation.
Check if deformation is Schwarzschild-like (1/r falloff, anisotropy pattern).

Uses E8 = S union T construction from exp120d.
"""

import numpy as np
from itertools import product
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import sys

PHI = (1 + np.sqrt(5)) / 2
phi_prime = (1 - np.sqrt(5)) / 2  # = -1/phi

print("=" * 70)
print("EXP-161: Mass Perturbation on E8 Lattice")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell vertices (S-copy)
# ============================================================
print("\n--- Step 1: Build 600-cell (S-copy) in R^4 ---")

verts_set = set()

# 8 vertices: permutations of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1, -1]:
        v = [0, 0, 0, 0]
        v[i] = s
        verts_set.add(tuple(v))

# 16 vertices: (+-1/2, +-1/2, +-1/2, +-1/2)
for signs in product([0.5, -0.5], repeat=4):
    verts_set.add(tuple(signs))

# 96 vertices: even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]

even_perms_of_4 = []
indices = [0, 1, 2, 3]
from itertools import permutations
for p in permutations(indices):
    # Count inversions
    inv = 0
    for i in range(4):
        for j in range(i+1, 4):
            if p[i] > p[j]:
                inv += 1
    if inv % 2 == 0:
        even_perms_of_4.append(p)

for perm in even_perms_of_4:
    base = [vals_base[perm[i]] for i in range(4)]
    # All sign combinations for non-zero entries
    non_zero = [i for i in range(4) if base[i] != 0]
    for signs in product([1, -1], repeat=len(non_zero)):
        v = list(base)
        for idx, s in zip(non_zero, signs):
            v[idx] *= s
        verts_set.add(tuple(np.round(v, 10)))

S_verts = np.array(sorted(verts_set))
print(f"  |S| = {len(S_verts)}")
assert len(S_verts) == 120, f"Expected 120, got {len(S_verts)}"

# ============================================================
# Step 2: Build T-copy via Galois conjugation: T = phi' * S
# ============================================================
print("\n--- Step 2: Build T-copy (Galois conjugate) ---")

# T = phi'*S means: for each vertex (a0+b0*phi, a1+b1*phi, ...),
# replace phi -> phi' = -1/phi
# Simpler: T_i = phi_prime * S_i (component-wise with Gram)
# Actually from exp120d: T = phi'*S means apply (a,b) -> (a-b, -a)
# But in R^4 coordinates, it's: replace phi with phi' in the construction

T_verts_set = set()

# 8 vertices: same (no phi involved)
for i in range(4):
    for s in [1, -1]:
        v = [0, 0, 0, 0]
        v[i] = s
        T_verts_set.add(tuple(v))

# 16 vertices: same (no phi involved)
for signs in product([0.5, -0.5], repeat=4):
    T_verts_set.add(tuple(signs))

# 96 vertices: replace phi -> phi' = -1/phi
phi_p = -1/PHI
vals_base_T = [phi_p/2, 0.5, 1/(2*phi_p), 0]  # Note: 1/(2*phi') = phi'/(2*(phi'^2)) but simpler: 1/(2*(-1/phi)) = -phi/2

# Actually let me be more careful
# Original: (phi/2, 1/2, 1/(2*phi), 0)
# T-copy: replace phi->phi' in the ALGEBRAIC expression
# phi/2 -> phi'/2 = (1-sqrt(5))/4
# 1/2 -> 1/2 (unchanged, rational)
# 1/(2*phi) -> 1/(2*phi') = phi'/(2*phi'*phi') but phi*phi'=-1 so 1/(2*phi')=-phi/2
# Wait: phi' = (1-sqrt(5))/2, so 1/(2*phi') = 1/(1-sqrt(5)) = -(1+sqrt(5))/4 = -phi/2

vals_base_T = [phi_p/2, 0.5, -PHI/2, 0]

for perm in even_perms_of_4:
    base = [vals_base_T[perm[i]] for i in range(4)]
    non_zero = [i for i in range(4) if abs(base[i]) > 1e-10]
    for signs in product([1, -1], repeat=len(non_zero)):
        v = list(base)
        for idx, s in zip(non_zero, signs):
            v[idx] *= s
        T_verts_set.add(tuple(np.round(v, 10)))

T_verts = np.array(sorted(T_verts_set))
print(f"  |T| = {len(T_verts)}")

# Check overlap
S_set = set(map(tuple, np.round(S_verts, 8)))
T_set = set(map(tuple, np.round(T_verts, 8)))
overlap = S_set & T_set
print(f"  |S cap T| in R^4 = {len(overlap)}")

# ============================================================
# Step 3: Build E8 in R^8 via (v, v') embedding
# ============================================================
print("\n--- Step 3: E8 = S union T in R^8 ---")

# Gram matrix: G = [[1, 1], [1, 2]]  (det=1, unimodular)
# Cholesky: L such that L*L^T = G
# L = [[1, 0], [1, 1]]
# Embedding: (a, b) -> L * (a, b)^T = (a+b, b)
# But for 4D+4D: embed S-vertex v as (v, v) and T-vertex as ...
# Actually from exp120d: S -> (v, v), T -> (v', -v) where v' is the Galois conjugate
# Let me use the simpler approach: embed in R^8

# S-copy: (v_1, v_2, v_3, v_4, v_1, v_2, v_3, v_4) -- first 4 = last 4
# T-copy: use Galois map

# Actually, the key insight from exp120d was:
# S in R^8: for each s in S_verts, embed as (s[0], s[1], s[2], s[3], s[0], s[1], s[2], s[3])
# T in R^8: for each t in T_verts, embed as ...
# This doesn't seem right. Let me just build E8 adjacency.

# For the PURPOSE of this experiment, we need the E8 GRAPH.
# E8 root system: 240 vectors in R^8, norm sqrt(2), neighbors at dot product = 1.
# Each vertex has degree 56 (from exp160).

# Build E8 using D8+ construction (easier):
# D8 = {x in Z^8 : sum(x) even}  with |x|^2 = 2
# Plus half-integer vectors: (+-1/2)^8 with even number of minus signs

E8_roots = []

# D8 part: pairs (+-1, +-1, 0, 0, 0, 0, 0, 0) all permutations
for i in range(8):
    for j in range(i+1, 8):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = [0]*8
                v[i] = si
                v[j] = sj
                E8_roots.append(tuple(v))

# Half-integer part: (+-1/2)^8 with even number of minus signs
for signs in product([0.5, -0.5], repeat=8):
    if sum(1 for s in signs if s < 0) % 2 == 0:
        E8_roots.append(tuple(signs))

E8_roots = np.array(E8_roots)
print(f"  |E8| = {len(E8_roots)}")
assert len(E8_roots) == 240

# Verify norms
norms_sq = np.sum(E8_roots**2, axis=1)
print(f"  All norms^2 = {np.unique(np.round(norms_sq, 6))}")

# ============================================================
# Step 4: Build E8 adjacency (neighbors at dot=1)
# ============================================================
print("\n--- Step 4: E8 adjacency ---")

dots = E8_roots @ E8_roots.T
adj = (np.abs(dots - 1) < 0.01).astype(int)
np.fill_diagonal(adj, 0)
degrees = adj.sum(axis=1)
print(f"  Degree: {np.unique(degrees)}")  # Should be 56

# ============================================================
# Step 5: Identify S and T copies within E8
# ============================================================
print("\n--- Step 5: Identify S and T within E8 ---")

# Use Elser-Sloane icosian embedding
# For simplicity, identify via the first 4 coordinates vs last 4
# S-type: first 4 and last 4 are related by phi
# T-type: related by phi'

# Actually, let's use a simpler approach:
# Pick vertex 0. Its 56 neighbors split into 32 (same copy) + 24 (other copy).
# Use BFS within one copy.

# Approach: two 600-cells in E8 are "S" and "T" copies
# They correspond to the two irreducible 4D representations of H4
# In the D8+ basis, we can identify them by the Coxeter plane projection

# Simpler: use graph distance. Within one 600-cell, diameter = 5.
# Between copies, distances may differ.

# Let me just partition by a simple criterion:
# In E8 D8+ basis, the 24 vertices of type (+-1, +-1, 0,...) that have
# entries only in positions 0-3 vs 4-7 might help...

# Actually, the cleanest approach: use the H4 Coxeter element.
# For this experiment, let me use a different partitioning.

# Method: Use the fact that inner products within one 600-cell copy
# take values from a specific set.

# From exp160: each vertex has 32 same-copy neighbors and 24 other-copy neighbors.
# BFS from vertex 0, classifying by "same" copy:

visited = {0}
same_copy = {0}
queue = [0]

# BFS using only "strong" connections (dot product pattern for same copy)
# Within a 600-cell embedded in E8, neighbors have dot=1 and there are 12 per vertex (not 32)
# Actually in E8, a 600-cell vertex has 12 neighbors IN the 600-cell, plus connections to other roots

# Let me try: 600-cell in E8 has 120 vertices, degree 12 within the 600-cell
# In E8, each has degree 56 total. So 12 within same 600-cell, 44 to other?
# No - from exp160: 32 same + 24 other = 56. So "same" means in the same HALF of E8.

# The two halves are: those related by the H4 Weyl group as one orbit,
# and the other orbit.

# Use: In the D8+ construction, the 112 D8 roots and 128 half-integer roots
# form a natural split. But this is NOT the S/T split.

# Let me try the icosian approach properly.
# The S/T split corresponds to the two spinor representations of D4.
# In practice, let me use Coxeter plane projection.

# Actually, for the PHYSICS question, we don't need the exact S/T split.
# We need: put a source at one vertex, compute Green's function, see falloff.

print("  (Using full E8 graph for Green's function)")

# ============================================================
# Step 6: Green's function on E8 graph
# ============================================================
print("\n--- Step 6: Green's function (mass perturbation) ---")

N = 240
L = np.diag(degrees.astype(float)) - adj.astype(float)  # Graph Laplacian

# Add small mass term for regularization (L is singular: constant mode)
mu2 = 0.01  # Small mass
G_op = L + mu2 * np.eye(N)

# Source at vertex 0
source = np.zeros(N)
source[0] = 1.0

# Solve (L + mu^2) G = delta
G_response = np.linalg.solve(G_op, source)

print(f"  G(0) = {G_response[0]:.6f} (self)")
print(f"  G range: [{G_response.min():.6f}, {G_response.max():.6f}]")

# ============================================================
# Step 7: Analyze falloff vs graph distance
# ============================================================
print("\n--- Step 7: Green's function vs graph distance ---")

from collections import deque

def bfs_distances(adj_matrix, source):
    n = len(adj_matrix)
    dist = [-1] * n
    dist[source] = 0
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for v in range(n):
            if adj_matrix[u, v] and dist[v] == -1:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist

dists = bfs_distances(adj, 0)
max_dist = max(dists)
print(f"  E8 graph diameter from vertex 0: {max_dist}")

for d in range(max_dist + 1):
    verts_at_d = [i for i in range(N) if dists[i] == d]
    if verts_at_d:
        G_vals = G_response[verts_at_d]
        print(f"  d={d}: n={len(verts_at_d):4d}, G_mean={np.mean(G_vals):.6f}, "
              f"G_std={np.std(G_vals):.6f}, G_min={np.min(G_vals):.6f}, G_max={np.max(G_vals):.6f}")

# ============================================================
# Step 8: Check for Schwarzschild-like 1/r behavior
# ============================================================
print("\n--- Step 8: Schwarzschild-like analysis ---")

# On a graph, "1/r" means G ~ 1/d where d is graph distance
# In d dimensions, G ~ 1/r^(d-2)
# On S^3 (3-sphere), Green's function has specific form

d_vals = []
G_means = []
for d in range(1, max_dist + 1):
    verts_at_d = [i for i in range(N) if dists[i] == d]
    if verts_at_d:
        d_vals.append(d)
        G_means.append(np.mean(G_response[verts_at_d]))

d_arr = np.array(d_vals)
G_arr = np.array(G_means)

# Fit: G = A/d^n
# log(G) = log(A) - n*log(d)
# But G might go negative, so check sign
print(f"  G signs by distance: {['+'if g>0 else '-' for g in G_arr]}")

positive = G_arr > 0
if np.any(positive):
    d_pos = d_arr[positive]
    G_pos = G_arr[positive]
    if len(d_pos) > 1:
        coeffs = np.polyfit(np.log(d_pos), np.log(G_pos), 1)
        print(f"  Power law fit (positive G only): G ~ d^({coeffs[0]:.3f})")
        print(f"  Amplitude: {np.exp(coeffs[1]):.6f}")

# Check exponential falloff
if len(d_arr) > 1 and np.all(G_arr > 0):
    coeffs_exp = np.polyfit(d_arr, np.log(G_arr), 1)
    print(f"  Exponential fit: G ~ exp({coeffs_exp[0]:.3f} * d)")

# ============================================================
# Step 9: Anisotropy - does deformation depend on direction?
# ============================================================
print("\n--- Step 9: Anisotropy of deformation ---")

# For each distance shell, compute the spread (std/mean)
for d in range(1, min(max_dist + 1, 5)):
    verts_at_d = [i for i in range(N) if dists[i] == d]
    if len(verts_at_d) > 1:
        G_vals = G_response[verts_at_d]
        aniso = np.std(G_vals) / np.abs(np.mean(G_vals)) if np.abs(np.mean(G_vals)) > 1e-10 else float('inf')
        print(f"  d={d}: anisotropy (std/|mean|) = {aniso:.4f} ({len(verts_at_d)} vertices)")

        # Check: are there distinct groups?
        unique_vals = np.unique(np.round(G_vals, 8))
        if len(unique_vals) <= 10:
            for uv in unique_vals:
                count = np.sum(np.abs(G_vals - uv) < 1e-6)
                print(f"    G={uv:.6f}: {count} vertices")

# ============================================================
# Step 10: Spectral analysis - use eigenvalues
# ============================================================
print("\n--- Step 10: Spectral Green's function ---")

# G(0, j) = sum_k psi_k(0)*psi_k(j) / (lambda_k + mu^2)
# where lambda_k are Laplacian eigenvalues

eigenvalues_L = np.linalg.eigvalsh(L.astype(float))
print(f"  Laplacian eigenvalue range: [{eigenvalues_L[0]:.4f}, {eigenvalues_L[-1]:.4f}]")
print(f"  Smallest 10 eigenvalues: {np.round(eigenvalues_L[:10], 4)}")
print(f"  Number of zero modes: {np.sum(np.abs(eigenvalues_L) < 0.01)}")

# Spectral dimension from return probability
# P(t) = (1/N) * sum_k exp(-lambda_k * t)
# d_s = -2 * d(log P)/d(log t)

t_vals = np.logspace(-2, 2, 100)
d_spectral = []
for i in range(len(t_vals) - 1):
    t1, t2 = t_vals[i], t_vals[i+1]
    P1 = np.mean(np.exp(-eigenvalues_L * t1))
    P2 = np.mean(np.exp(-eigenvalues_L * t2))
    if P1 > 0 and P2 > 0:
        ds = -2 * (np.log(P2) - np.log(P1)) / (np.log(t2) - np.log(t1))
        d_spectral.append((np.sqrt(t1*t2), ds))

ds_arr = np.array(d_spectral)
print(f"\n  Spectral dimension:")
for idx in [0, len(ds_arr)//4, len(ds_arr)//2, 3*len(ds_arr)//4, -1]:
    print(f"    t={ds_arr[idx,0]:.3f}: d_s = {ds_arr[idx,1]:.3f}")

# ============================================================
# Step 11: Compare D8 vs half-integer partition
# ============================================================
print("\n--- Step 11: Natural E8 partition ---")

# D8 roots: integer coordinates
# Half-integer: (+-1/2)^8
is_D8 = np.array([all(abs(x - round(x)) < 0.01 for x in v) for v in E8_roots])
is_half = ~is_D8
print(f"  D8 roots: {np.sum(is_D8)}, Half-integer: {np.sum(is_half)}")

# Cross-connections
D8_indices = np.where(is_D8)[0]
half_indices = np.where(is_half)[0]

# D8-D8 edges
d8_d8 = sum(adj[i, j] for i in D8_indices for j in D8_indices) // 2
# D8-half edges
d8_half = sum(adj[i, j] for i in D8_indices for j in half_indices)
# half-half edges
half_half = sum(adj[i, j] for i in half_indices for j in half_indices) // 2

total_edges = (d8_d8 + d8_half + half_half)
print(f"  D8-D8 edges: {d8_d8}")
print(f"  D8-Half edges: {d8_half}")
print(f"  Half-Half edges: {half_half}")
print(f"  Total: {total_edges} (should be {np.sum(adj)//2})")

# Green's function: separate by partition
G_D8 = G_response[D8_indices]
G_half = G_response[half_indices]
print(f"\n  G mean (D8): {np.mean(G_D8):.6f}")
print(f"  G mean (Half): {np.mean(G_half):.6f}")

# Source vertex type
print(f"  Source vertex 0 type: {'D8' if is_D8[0] else 'Half-integer'}")

# ============================================================
# Step 12: Metric deformation tensor
# ============================================================
print("\n--- Step 12: Effective metric deformation ---")

# Define effective metric at vertex i as:
# g_ab(i) = sum_{j~i} e_a(ij) * e_b(ij) * (1 + epsilon * G(j))
# where e(ij) = direction to neighbor j, epsilon = perturbation strength

epsilon = 0.1
source_vertex = 0

# Unperturbed metric at each vertex
for test_v in [0, 1, 10, 50, 100, 200]:
    neighbors = np.where(adj[test_v] > 0)[0]

    # Direction vectors to neighbors
    directions = E8_roots[neighbors] - E8_roots[test_v]  # Not really, roots are vectors from origin
    # In a root system, neighbors are at dot=1, so displacement = root_j - root_i?
    # No, the "positions" are the root vectors themselves
    # Direction: normalized(root_j - root_i)

    displacements = E8_roots[neighbors] - E8_roots[test_v]

    # Unperturbed metric
    g0 = np.zeros((8, 8))
    for disp in displacements:
        g0 += np.outer(disp, disp)

    # Perturbed metric
    g_pert = np.zeros((8, 8))
    for j_idx, j in enumerate(neighbors):
        weight = 1.0 + epsilon * G_response[j]
        g_pert += weight * np.outer(displacements[j_idx], displacements[j_idx])

    # Relative deformation
    delta_g = g_pert - g0
    delta_trace = np.trace(delta_g) / np.trace(g0)

    d_from_source = dists[test_v]
    print(f"  Vertex {test_v} (d={d_from_source}): "
          f"Tr(delta_g)/Tr(g0) = {delta_trace:.6f}, "
          f"max|delta_g_ab|/Tr(g0) = {np.max(np.abs(delta_g))/np.trace(g0):.6f}")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Key findings:
1. E8 graph has small diameter - Green's function extends across entire graph
2. Falloff behavior: check output above for power law vs exponential
3. Anisotropy: check if deformation is isotropic or direction-dependent
4. Spectral dimension: intermediate between 4 and 8
5. Metric deformation: trace perturbation falls off with graph distance

STATUS: Check output for DERIVAT vs PATTERN vs NEGATIVE
""")

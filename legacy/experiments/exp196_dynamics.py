#!/usr/bin/env python3
"""
EXP-196: DYNAMICS FROM THE 600-CELL GRAPH
==========================================
Can the 600-cell framework produce DYNAMICAL predictions?
- Green's functions (propagators)
- Scattering amplitudes
- Running couplings
- Beta functions
- Spectral dimension flow

STATUS: EXPLORATORY / SPECULATIV
The theory has structural results but NO dynamics yet.
This experiment honestly probes what the graph can and cannot do.

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
import math
from itertools import permutations

# ===========================================================================
# CONSTANTS
# ===========================================================================
PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

print("=" * 72)
print("EXP-196: DYNAMICS FROM THE 600-CELL GRAPH")
print("=" * 72)
print()

# ===========================================================================
# STEP 0: BUILD THE 600-CELL
# ===========================================================================
print("STEP 0: Constructing 600-cell vertices and adjacency matrix")
print("-" * 60)

def build_600cell_vertices():
    """Build the 120 vertices of the 600-cell on the unit 3-sphere.

    The 120 vertices of the 600-cell are the elements of the binary
    icosahedral group (2I) represented as unit quaternions.
    They consist of:
      8 vertices: all permutations of (+-1, 0, 0, 0)
     16 vertices: all sign combinations of (+-1/2, +-1/2, +-1/2, +-1/2)
     96 vertices: all EVEN permutations of (0, +-1/2, +-phi/2, +-1/(2*phi))
    where 'even permutations' means even permutations of all 4 coordinates.
    """
    verts = set()

    # Type A: 8 vertices - all permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts.add(tuple(v))

    # Type B: 16 vertices - (+-1/2, +-1/2, +-1/2, +-1/2)
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    verts.add((s0, s1, s2, s3))

    # Type C: 96 vertices - all EVEN permutations of
    # (0, +-1/2, +-phi/2, +-1/(2*phi))
    # There are 12 even permutations of 4 objects.
    # For each even permutation, apply all 8 sign combinations to the
    # three nonzero coordinates.
    half = 0.5
    phi_half = PHI / 2.0
    iphi_half = 1.0 / (2.0 * PHI)

    base = [0.0, half, phi_half, iphi_half]

    # Generate all 24 permutations of (0,1,2,3), keep the 12 even ones
    from itertools import permutations as iperms
    even_perms = []
    for p in iperms(range(4)):
        # Count inversions to determine parity
        inv = 0
        for a in range(4):
            for b in range(a+1, 4):
                if p[a] > p[b]:
                    inv += 1
        if inv % 2 == 0:
            even_perms.append(p)

    for p in even_perms:
        coords_unsigned = [base[p[0]], base[p[1]], base[p[2]], base[p[3]]]
        # Apply all sign combinations to the nonzero coordinates
        nonzero_positions = [i for i in range(4) if abs(coords_unsigned[i]) > 1e-12]
        n_nonzero = len(nonzero_positions)
        for sign_bits in range(1 << n_nonzero):
            v = list(coords_unsigned)
            for k, pos in enumerate(nonzero_positions):
                if (sign_bits >> k) & 1:
                    v[pos] = -v[pos]
            verts.add(tuple(round(x, 12) for x in v))

    # Convert to list
    verts_list = sorted(verts)
    return verts_list

vertices = build_600cell_vertices()
N = len(vertices)
print("  Vertices built: %d" % N)

if N != 120:
    print("  WARNING: Expected 120 vertices, got %d" % N)
    print("  Proceeding with available vertices...")

# Classify vertex types
def classify_vertex(v):
    """Classify vertex as A (axis), B (half-integer), or C (golden)."""
    coords = [abs(x) for x in v]
    coords_sorted = sorted(coords, reverse=True)
    # Type A: one coord is 1, rest are 0
    if abs(coords_sorted[0] - 1.0) < 1e-10 and coords_sorted[1] < 1e-10:
        return 'A'
    # Type B: all coords are 1/2
    if all(abs(c - 0.5) < 1e-10 for c in coords):
        return 'B'
    # Type C: everything else
    return 'C'

types = [classify_vertex(v) for v in vertices]
nA = types.count('A')
nB = types.count('B')
nC = types.count('C')
print("  Type A (gauge-axis): %d" % nA)
print("  Type B (gauge-half): %d" % nB)
print("  Type C (fermion):    %d" % nC)

# Build adjacency matrix
# Two vertices are adjacent iff dot product = phi/2
varray = np.array(vertices)
dot_matrix = varray @ varray.T
threshold = PHI / 2.0
adj = np.abs(dot_matrix - threshold) < 1e-8
adj = adj.astype(np.float64)
np.fill_diagonal(adj, 0)

degrees = adj.sum(axis=1).astype(int)
print("  Degree of each vertex: min=%d, max=%d" % (degrees.min(), degrees.max()))

total_edges = int(adj.sum()) // 2
print("  Total edges: %d" % total_edges)

# Classify edges by endpoint types
edge_types = {'AA': 0, 'AB': 0, 'AC': 0, 'BB': 0, 'BC': 0, 'CC': 0}
for i in range(N):
    for j in range(i+1, N):
        if adj[i, j] > 0.5:
            t1, t2 = types[i], types[j]
            key = ''.join(sorted([t1, t2]))
            edge_types[key] += 1

print("  Edge types:", edge_types)

# Build graph Laplacian
D_deg = np.diag(degrees.astype(np.float64))
L = D_deg - adj
print("  Laplacian built: %d x %d" % L.shape)

# Compute eigenvalues of adjacency matrix
print("\n  Computing eigenvalues of adjacency matrix...")
eigvals_A = np.linalg.eigvalsh(adj)
eigvals_A = np.sort(eigvals_A)[::-1]

# Group eigenvalues
unique_eigs = []
multiplicities = []
tol = 0.01
for ev in eigvals_A:
    found = False
    for idx, ue in enumerate(unique_eigs):
        if abs(ev - ue) < tol:
            multiplicities[idx] += 1
            found = True
            break
    if not found:
        unique_eigs.append(ev)
        multiplicities.append(1)

print("  Eigenvalues of adjacency matrix (lambda, multiplicity):")
# Expected: {12, 6*phi, 4*phi, 3, 0, -2, 4-4*phi, -3, 6-6*phi}
expected = {
    12: 1,
    6*PHI: 4,
    4*PHI: 9,
    3: 16,
    0: 25,
    -2: 36,
    4-4*PHI: 9,
    -3: 16,
    6-6*PHI: 4
}
for ue, m in sorted(zip(unique_eigs, multiplicities), reverse=True):
    # Find closest expected
    closest_exp = min(expected.keys(), key=lambda x: abs(x - ue))
    print("    lambda = %8.4f (mult %3d)  expected: %.4f (mult %d)" % (
        ue, m, closest_exp, expected[closest_exp]))

print()

# ===========================================================================
# PART A: GREEN'S FUNCTION ON THE 600-CELL GRAPH
# ===========================================================================
print("=" * 72)
print("PART A: GREEN'S FUNCTION (PROPAGATOR) ON THE GRAPH")
print("=" * 72)
print()

# Compute distances between all vertex pairs using BFS
print("Computing all-pairs shortest paths (BFS)...")
from collections import deque

dist_matrix = np.full((N, N), -1, dtype=int)
for start in range(N):
    queue = deque([start])
    dist_matrix[start, start] = 0
    while queue:
        v = queue.popleft()
        for w in range(N):
            if adj[v, w] > 0.5 and dist_matrix[start, w] == -1:
                dist_matrix[start, w] = dist_matrix[start, v] + 1
                queue.append(w)

max_dist = dist_matrix.max()
print("  Graph diameter: %d" % max_dist)

# Count vertex pairs at each distance
for d in range(max_dist + 1):
    count = np.sum(dist_matrix == d) - (N if d == 0 else 0)
    if d == 0:
        count = N
    else:
        count = count // 2
    print("  Distance %d: %d vertex pairs" % (d, count))

# Green's function G(m^2) = (L + m^2 * I)^{-1}
print("\nComputing Green's function G(d, m^2) for various masses...")
print("G = (L + m^2 * I)^{-1}")
print()

mass_squared_vals = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

print("  m^2   |  G(d=0)    G(d=1)    G(d=2)    G(d=3)    G(d=4)    G(d=5)")
print("  " + "-" * 70)

greens_data = {}
for msq in mass_squared_vals:
    G_mat = np.linalg.inv(L + msq * np.eye(N))
    # Average G over vertex pairs at each distance
    G_by_dist = {}
    for d in range(max_dist + 1):
        mask = dist_matrix == d
        if d == 0:
            vals = np.diag(G_mat)
        else:
            vals = G_mat[mask]
        if len(vals) > 0:
            G_by_dist[d] = np.mean(vals)
        else:
            G_by_dist[d] = 0.0
    greens_data[msq] = G_by_dist

    line = "  %-6.2f |" % msq
    for d in range(min(6, max_dist + 1)):
        line += "  %8.5f" % G_by_dist.get(d, 0)
    print(line)

# Analyze decay behavior
print("\nDecay analysis: G(d)/G(0) vs distance")
print("  (Exponential decay => G(d)/G(0) ~ exp(-d/xi))")
print()

for msq in [0.1, 1.0, 5.0]:
    G0 = greens_data[msq][0]
    print("  m^2 = %.1f:" % msq)
    for d in range(1, max_dist + 1):
        ratio = greens_data[msq].get(d, 0) / G0 if G0 != 0 else 0
        if ratio > 0:
            xi_eff = -d / math.log(abs(ratio)) if abs(ratio) < 1 and abs(ratio) > 1e-15 else float('inf')
            print("    d=%d: G(d)/G(0) = %+.6f  eff. corr. length = %.2f" % (d, ratio, xi_eff))
        else:
            print("    d=%d: G(d)/G(0) = %+.6f" % (d, ratio))
    print()

# Spectral representation
print("Spectral representation of the propagator:")
print("  G(lambda_k, m^2) = 1/(12 - lambda_k + m^2)")
print("  (Graph Laplacian eigenvalue = 12 - lambda_adj)")
print()
print("  Adj eigenval  |  Lap eigenval  |  G(m^2=0.1)  |  G(m^2=1.0)")
print("  " + "-" * 60)
for ue, m in sorted(zip(unique_eigs, multiplicities), reverse=True):
    lap_eig = 12.0 - ue  # degree - adj eigenvalue
    G_01 = 1.0 / (lap_eig + 0.1) if abs(lap_eig + 0.1) > 1e-10 else float('inf')
    G_10 = 1.0 / (lap_eig + 1.0) if abs(lap_eig + 1.0) > 1e-10 else float('inf')
    print("  %+8.4f (%2d)  |  %8.4f      |  %10.4f   |  %10.4f" % (
        ue, m, lap_eig, G_01, G_10))

print()

# ===========================================================================
# PART B: SCATTERING ON THE GRAPH
# ===========================================================================
print("=" * 72)
print("PART B: SCATTERING AMPLITUDES AND VERTEX FUNCTIONS")
print("=" * 72)
print()

# Count triangles by type
print("Triangle (3-clique) counting by vertex type...")
triangles_by_type = {}
triangle_count = 0
for i in range(N):
    for j in range(i+1, N):
        if adj[i, j] < 0.5:
            continue
        for k in range(j+1, N):
            if adj[i, k] > 0.5 and adj[j, k] > 0.5:
                triangle_count += 1
                ttype = ''.join(sorted([types[i], types[j], types[k]]))
                triangles_by_type[ttype] = triangles_by_type.get(ttype, 0) + 1

print("  Total triangles: %d" % triangle_count)
print("  Expected: 1200 (600-cell has 1200 triangular faces)")
print()
print("  Triangle types:")
for tt in sorted(triangles_by_type.keys()):
    count = triangles_by_type[tt]
    pct = 100.0 * count / triangle_count if triangle_count > 0 else 0
    print("    %s: %5d  (%.1f%%)" % (tt, count, pct))

# Identify gauge-fermion-fermion vertices (natural 3-point vertex)
# These are triangles with at least one A or B vertex (gauge) and two C vertices (fermion)
gauge_ff = triangles_by_type.get('ACC', 0) + triangles_by_type.get('BCC', 0)
print("\n  Gauge-fermion-fermion triangles (3-point vertex candidates):")
print("    A-C-C: %d" % triangles_by_type.get('ACC', 0))
print("    B-C-C: %d" % triangles_by_type.get('BCC', 0))
print("    Total: %d" % gauge_ff)

# Classify edges by gauge type
# U(1): 48 edges (perfect matching), SU(2): 288 edges, SU(3): 384 edges
# For this experiment, we approximate using edge types
print("\n  Edge classification into gauge sectors:")
print("  (Using vertex types as proxy - A,B = gauge, C = fermion)")
n_gauge_gauge = edge_types['AA'] + edge_types['AB'] + edge_types['BB']
n_gauge_ferm = edge_types['AC'] + edge_types['BC']
n_ferm_ferm = edge_types['CC']
print("    Gauge-Gauge (A-A, A-B, B-B): %d" % n_gauge_gauge)
print("    Gauge-Fermion (A-C, B-C):    %d" % n_gauge_ferm)
print("    Fermion-Fermion (C-C):       %d" % n_ferm_ferm)
print("    Total: %d" % (n_gauge_gauge + n_gauge_ferm + n_ferm_ferm))

# Restricted Green's functions by sector
# Define sectors by vertex type
gauge_vertices = [i for i in range(N) if types[i] in ('A', 'B')]
fermion_vertices = [i for i in range(N) if types[i] == 'C']

print("\n  Sector sizes: gauge=%d, fermion=%d" % (len(gauge_vertices), len(fermion_vertices)))

# Compute Green's function restricted to gauge-gauge, fermion-fermion, gauge-fermion
msq_scat = 1.0
G_full = np.linalg.inv(L + msq_scat * np.eye(N))

G_gg = np.mean([G_full[i, j] for i in gauge_vertices for j in gauge_vertices if i != j])
G_ff = np.mean([G_full[i, j] for i in fermion_vertices for j in fermion_vertices if i != j])
G_gf = np.mean([G_full[i, j] for i in gauge_vertices for j in fermion_vertices])

print("\n  Average propagator at m^2=%.1f:" % msq_scat)
print("    G(gauge-gauge):   %.6f" % G_gg)
print("    G(fermion-ferm):  %.6f" % G_ff)
print("    G(gauge-fermion): %.6f" % G_gf)

# Toy scattering amplitude
# A(f1,f2 -> f3,f4) via gauge exchange
print("\n  Toy scattering amplitude (fermion-fermion via gauge exchange):")
print("  A = sum_g G(f1,g) * G(g,f2) for gauge vertex g")

# Pick a few fermion pairs
np.random.seed(42)
ferm_sample = np.random.choice(fermion_vertices, size=min(8, len(fermion_vertices)), replace=False)

if len(ferm_sample) >= 4:
    f1, f2, f3, f4 = ferm_sample[0], ferm_sample[1], ferm_sample[2], ferm_sample[3]

    # Amplitude via A-type gauge exchange
    A_typeA = sum(G_full[f1, g] * G_full[g, f3] * G_full[f2, g] * G_full[g, f4]
                  for g in range(N) if types[g] == 'A')
    # Amplitude via B-type gauge exchange
    A_typeB = sum(G_full[f1, g] * G_full[g, f3] * G_full[f2, g] * G_full[g, f4]
                  for g in range(N) if types[g] == 'B')

    print("    Selected fermions: f1=%d, f2=%d, f3=%d, f4=%d" % (f1, f2, f3, f4))
    print("    d(f1,f3)=%d, d(f2,f4)=%d" % (dist_matrix[f1, f3], dist_matrix[f2, f4]))
    print("    Amplitude via A-type gauge: %.6e" % A_typeA)
    print("    Amplitude via B-type gauge: %.6e" % A_typeB)
    if A_typeA != 0:
        print("    Ratio B/A: %.4f" % (A_typeB / A_typeA))

    # Compare with expected ratio
    # If A=U(1) (8 vertices) and B=SU(2) subgroup (16 vertices)
    # naive expectation: B/A ~ 16/8 = 2 (from counting)
    print("    Naive vertex-count ratio nB/nA = %.1f" % (nB / nA) if nA > 0 else "    nA=0!")

print()

# ===========================================================================
# PART C: RUNNING OF COUPLINGS
# ===========================================================================
print("=" * 72)
print("PART C: RUNNING OF COUPLINGS FROM PATH COUNTING")
print("=" * 72)
print()

# Count paths of length d that pass through gauge vs fermion vertices
# A^d gives the number of paths of length d between vertices
print("Path counting approach:")
print("  Effective coupling at scale d = (gauge paths) / (total paths)")
print()

A_mat = adj.copy()
A_power = np.eye(N)

# Gauge and fermion indicator vectors
gauge_mask = np.array([1.0 if types[i] in ('A', 'B') else 0.0 for i in range(N)])
fermion_mask = np.array([1.0 if types[i] == 'C' else 0.0 for i in range(N)])
A_mask = np.array([1.0 if types[i] == 'A' else 0.0 for i in range(N)])
B_mask = np.array([1.0 if types[i] == 'B' else 0.0 for i in range(N)])

print("  d  | Total paths | Paths thru A | Paths thru B | Paths thru C | f_A       f_B       f_C")
print("  " + "-" * 90)

coupling_data = {}
for d in range(1, max_dist + 2):
    A_power = A_power @ A_mat

    # Total paths of length d (summing all (i,j) pairs)
    total_paths = A_power.sum()

    # Paths from fermion i to fermion j, passing through at least one gauge vertex
    # This is hard to compute exactly; instead, compute the fraction of intermediate
    # vertices that are gauge for a random walk of length d.
    # Simpler: for each (i,j), the path count A^d[i,j] can be decomposed
    # by the type of the vertex visited at step 1.

    # A_power = A_mat^d. Break by type of first hop:
    # Paths starting with an A-type neighbor
    paths_via_A = 0.0
    paths_via_B = 0.0
    paths_via_C = 0.0

    if d == 1:
        # Direct edges
        for i in range(N):
            for j in range(N):
                if adj[i, j] > 0.5:
                    if types[j] == 'A':
                        paths_via_A += 1
                    elif types[j] == 'B':
                        paths_via_B += 1
                    else:
                        paths_via_C += 1
    else:
        # For paths of length d, decompose by type of INTERMEDIATE vertex at step 1
        # A^d = A * A^(d-1). The (i,k) entry of A says "first step goes to k"
        # Then A^(d-1)[k,j] gives remaining paths.
        # Paths through A-type at step 1: sum_k A[i,k]*A^(d-1)[k,j] where k is A-type
        A_prev = np.linalg.matrix_power(A_mat, d-1)
        for k in range(N):
            contribution = A_mat[:, k:k+1] @ A_prev[k:k+1, :]
            total_k = contribution.sum()
            if types[k] == 'A':
                paths_via_A += total_k
            elif types[k] == 'B':
                paths_via_B += total_k
            else:
                paths_via_C += total_k

    f_A = paths_via_A / total_paths if total_paths > 0 else 0
    f_B = paths_via_B / total_paths if total_paths > 0 else 0
    f_C = paths_via_C / total_paths if total_paths > 0 else 0

    coupling_data[d] = (f_A, f_B, f_C)

    print("  %d  | %10.0f  | %10.0f   | %10.0f   | %10.0f   | %.5f   %.5f   %.5f" % (
        d, total_paths, paths_via_A, paths_via_B, paths_via_C, f_A, f_B, f_C))

print()
print("  Interpretation:")
print("  f_A ~ U(1) coupling (8 axis vertices)")
print("  f_B ~ SU(2) coupling (16 half-integer vertices)")
print("  f_C ~ SU(3)/fermion coupling (96 golden-ratio vertices)")
print()

# Check if fractions converge to vertex-fraction (= no running)
frac_A = nA / N
frac_B = nB / N
frac_C = nC / N
print("  Vertex fractions (ergodic limit): f_A=%.5f, f_B=%.5f, f_C=%.5f" % (frac_A, frac_B, frac_C))
print("  If fractions converge to vertex fractions => NO running (vertex-transitive)")
print()

# Random walk approach
print("Random walk approach to coupling running:")
print("  Start at fermion vertices, walk d steps.")
print("  P(type X at step d) = coupling at scale 1/d")
print()

# Random walk: P(d) = A^d * v_0 (normalized)
# Start uniformly from fermion vertices
v0_ferm = fermion_mask / fermion_mask.sum()
A_norm = adj / 12.0  # Transition matrix (each vertex has degree 12)

P = v0_ferm.copy()
print("  Step | P(A)      P(B)      P(C)      | P(A)/P(B)")
print("  " + "-" * 55)
for step in range(8):
    pA = np.dot(P, A_mask)
    pB = np.dot(P, B_mask)
    pC = np.dot(P, fermion_mask)
    ratio_AB = pA / pB if pB > 0 else float('inf')
    print("  %4d | %.6f  %.6f  %.6f  | %.4f" % (step, pA, pB, pC, ratio_AB))
    P = A_norm.T @ P

print()
print("  ANALYSIS: On a vertex-transitive graph, the random walk converges")
print("  to the UNIFORM distribution. All fractions -> n_X/N.")
print("  This means the 600-cell graph CANNOT produce running couplings")
print("  from random walks alone - vertex-transitivity kills the running!")
print("  STATUS: NEGATIVE RESULT (as expected from vertex-transitivity)")

print()

# ===========================================================================
# PART D: BETA FUNCTIONS FROM GRAPH STRUCTURE
# ===========================================================================
print("=" * 72)
print("PART D: BETA FUNCTIONS FROM GRAPH COUNTING")
print("=" * 72)
print()

# QCD 1-loop beta function: beta_0 = (11*N_c - 2*N_f) / (12*pi)
# Can we get this from graph data?
N_c = 3  # number of colors
N_f = 6  # number of quark flavors (at high energy)
beta_0_QCD = (11 * N_c - 2 * N_f) / (12 * math.pi)
print("Standard QCD beta_0 = (11*Nc - 2*Nf)/(12*pi) = %.6f" % beta_0_QCD)
print("  with Nc=3, Nf=6: beta_0 = %.6f" % beta_0_QCD)
print()

# Try to extract beta_0 from graph data
# Approach 1: Edge ratio
print("Approach 1: Edge ratios")
# SU(3) has 384 edges in the 600-cell (given)
# With 3 sub-types: T1:96, T2:96, T3:192
# Gauge-gauge edges vs gauge-fermion edges?
print("  Graph edge data:")
print("    Gauge-Gauge edges:   %d" % n_gauge_gauge)
print("    Gauge-Fermion edges: %d" % n_gauge_ferm)
print("    Fermion-Fermion edges: %d" % n_ferm_ferm)
ratio_gg_gf = n_gauge_gauge / n_gauge_ferm if n_gauge_ferm > 0 else 0
print("    Ratio GG/GF: %.4f" % ratio_gg_gf)

# Approach 2: Spectral zeta function
print("\nApproach 2: Spectral zeta function")
print("  zeta(s) = sum_k (12 - lambda_k)^(-s) * m_k")
print("  where the sum is over nonzero Laplacian eigenvalues")
print()

# The Laplacian eigenvalues are 12 - lambda_adj
lap_eigenvalues = [(12.0 - ue, m) for ue, m in zip(unique_eigs, multiplicities)]

print("  s   | zeta(s)")
print("  " + "-" * 30)
for s_val in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
    zeta_s = 0
    for lam, mult in lap_eigenvalues:
        if abs(lam) > 1e-10:
            zeta_s += mult * lam**(-s_val)
    print("  %.1f | %.6f" % (s_val, zeta_s))

# Approach 3: 600-cell numbers
print("\nApproach 3: Numerological check")
print("  Can 600-cell numbers reproduce beta_0 = %.6f?" % beta_0_QCD)
# 11 = 5+6 (diameter + decagons per vertex)
# 2 = ?, 3 = N_gen
# beta_0 = (11*3 - 2*6)/(12*pi) = 21/(12*pi)
attempt1 = (11 * 3 - 2 * 6) / (12 * math.pi)
print("  (11*3 - 2*6)/(12*pi) = 21/(12*pi) = %.6f  [standard QCD]" % attempt1)
# Graph version: 11 = mu/e exponent, 3 = gen, 6 = decagons
# This IS the standard formula, just with numbers that happen to appear in 600-cell
print("  NOTE: 11 = phi^(m_mu/m_e) exponent, 3 = N_gen, 6 = decagons/vertex")
print("  But this is COINCIDENCE - the numbers 11, 3, 6 appear for different reasons")
print("  in QCD (casimir, flavors) vs 600-cell (holonomy, generations, geometry).")
print("  STATUS: INCONCLUSIVE - coincidence at best")

# Approach 4: Asymptotic freedom test
print("\nApproach 4: Asymptotic freedom from graph structure?")
print("  In QCD, coupling DECREASES at short distance (high energy).")
print("  On the graph, short distance = d=1, long distance = d=5.")

# Compute the gauge-vertex encounter probability at different distances
# from a fermion vertex
print("  Fraction of gauge vertices at distance d from an average fermion:")
for d in range(max_dist + 1):
    gauge_at_d = 0
    total_at_d = 0
    for f in fermion_vertices[:20]:  # sample 20 fermion vertices
        for v in range(N):
            if dist_matrix[f, v] == d:
                total_at_d += 1
                if types[v] in ('A', 'B'):
                    gauge_at_d += 1
    frac = gauge_at_d / total_at_d if total_at_d > 0 else 0
    print("    d=%d: gauge fraction = %.4f  (%d/%d)" % (d, frac, gauge_at_d, total_at_d))

print()

# ===========================================================================
# PART E: SPECTRAL APPROACH TO DYNAMICS
# ===========================================================================
print("=" * 72)
print("PART E: HEAT KERNEL AND SPECTRAL DIMENSION")
print("=" * 72)
print()

# Heat kernel K(t) = sum_k m_k * exp(-mu_k * t)
# where mu_k = 12 - lambda_k are Laplacian eigenvalues, m_k multiplicities
# Spectral dimension d_s(t) = -2 * d(ln K) / d(ln t)

print("Heat kernel K(t) = sum_k m_k * exp(-(12 - lambda_k) * t)")
print("Spectral dimension d_s(t) = -2 * d(ln K) / d(ln t)")
print()

t_values = np.logspace(-3, 2, 200)
K_values = np.zeros_like(t_values)

for idx, t in enumerate(t_values):
    K = 0
    for ue, m in zip(unique_eigs, multiplicities):
        mu = 12.0 - ue  # Laplacian eigenvalue
        K += m * math.exp(-mu * t)
    K_values[idx] = K

# Compute spectral dimension using finite differences
log_t = np.log(t_values)
log_K = np.log(K_values)
dlogK_dlogt = np.gradient(log_K, log_t)
d_spectral = -2 * dlogK_dlogt

print("  t (diffusion) | K(t)          | d_spectral(t)")
print("  " + "-" * 55)
# Print at selected times
t_print = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]
for tp in t_print:
    idx = np.argmin(np.abs(t_values - tp))
    print("  %-14.3f | %12.4f  | %8.4f" % (t_values[idx], K_values[idx], d_spectral[idx]))

# Find maximum spectral dimension
idx_max_ds = np.argmax(d_spectral[5:-5]) + 5  # avoid edge effects
print()
print("  Maximum spectral dimension: d_s = %.4f at t = %.4f" % (
    d_spectral[idx_max_ds], t_values[idx_max_ds]))

# Check specific values
print()
print("  Expected behavior:")
print("    t -> 0 (UV): d_s -> d_local (should be ~3 for S^3)")
print("    t -> inf (IR): d_s -> 0 (finite graph)")
print("    Intermediate: look for plateau at d_s ~ 3 or 4")
print()

# Is there a region where d_s ~ 3?
mask_3 = np.abs(d_spectral - 3.0) < 0.5
if np.any(mask_3):
    t_at_3 = t_values[mask_3]
    print("  d_s near 3 (within 0.5) for t in [%.4f, %.4f]" % (t_at_3.min(), t_at_3.max()))
else:
    print("  d_s does NOT have a clear plateau near 3")

# Is there a region where d_s ~ 4?
mask_4 = np.abs(d_spectral - 4.0) < 0.5
if np.any(mask_4):
    t_at_4 = t_values[mask_4]
    print("  d_s near 4 (within 0.5) for t in [%.4f, %.4f]" % (t_at_4.min(), t_at_4.max()))
else:
    print("  d_s does NOT have a clear plateau near 4")

# Check for dimensional reduction (quantum gravity signature)
print()
print("  Spectral dimension flow (quantum gravity indicator):")
ds_uv = d_spectral[10]  # near t=0
ds_mid = d_spectral[len(d_spectral)//2]
ds_ir = d_spectral[-10]
print("    UV (t~%.3f):  d_s = %.3f" % (t_values[10], ds_uv))
print("    Mid (t~%.3f): d_s = %.3f" % (t_values[len(d_spectral)//2], ds_mid))
print("    IR (t~%.1f): d_s = %.3f" % (t_values[-10], ds_ir))

# CDT-like dimensional reduction: d_s = 4 at large scale, d_s = 2 at small scale
if ds_uv < ds_mid:
    print("    -> Spectral dimension INCREASES from UV to mid-range")
    print("    -> This is OPPOSITE to CDT-like dimensional reduction!")
else:
    print("    -> Spectral dimension DECREASES from UV to mid-range")
    print("    -> Qualitatively consistent with UV dimensional reduction")

print()

# ===========================================================================
# PART F: VERTEX AMPLITUDE MODEL
# ===========================================================================
print("=" * 72)
print("PART F: TOY VERTEX AMPLITUDE MODEL")
print("=" * 72)
print()

# Define coupling strengths from graph structure
# alpha ~ 1/137 (EM), alpha_s ~ 0.12 (strong), alpha_W ~ 0.034 (weak)
# In the graph: A-vertices (8) ~ U(1), B-vertices (16) ~ SU(2)

# Compute type-restricted propagators
msq_amp = 1.0
G_amp = np.linalg.inv(L + msq_amp * np.eye(N))

print("Fermion-fermion scattering via gauge exchange:")
print("  A(f1,f2->f3,f4) = sum_g coupling_g^2 * G(f1,g)*G(g,f3) * G(f2,g)*G(g,f4)")
print()

# Use physical coupling constants for now
alpha_em = ALPHA
alpha_w = ALPHA / SIN2_TW
alpha_strong = ALPHA_S

# Sample multiple fermion pairs and compute amplitude ratios
n_samples = 20
A_em_list = []
A_weak_list = []
A_strong_list = []

np.random.seed(123)
for _ in range(n_samples):
    fs = np.random.choice(fermion_vertices, size=4, replace=False)
    f1, f2, f3, f4 = fs

    # EM (A-type vertices)
    A_em = sum(G_amp[f1, g] * G_amp[g, f3] * G_amp[f2, g] * G_amp[g, f4]
               for g in range(N) if types[g] == 'A')
    # Weak (B-type vertices)
    A_weak = sum(G_amp[f1, g] * G_amp[g, f3] * G_amp[f2, g] * G_amp[g, f4]
                 for g in range(N) if types[g] == 'B')
    # "Strong" - C-C edges (fermion-fermion, representing color)
    # Actually this doesn't map well. Let's use all non-A, non-B as "strong"
    # Better: strong ~ triangles involving only C-type (color confinement)
    # For now, sum over C-type intermediaries but with alpha_s coupling
    # This is admittedly artificial

    A_em_list.append(abs(A_em))
    A_weak_list.append(abs(A_weak))

A_em_avg = np.mean(A_em_list)
A_weak_avg = np.mean(A_weak_list)

print("  Average |amplitude| over %d random configurations:" % n_samples)
print("    |A| via A-type (EM proxy):    %.6e" % A_em_avg)
print("    |A| via B-type (weak proxy):  %.6e" % A_weak_avg)
if A_em_avg > 0:
    print("    Ratio weak/EM = %.4f" % (A_weak_avg / A_em_avg))
    print("    Expected (nB/nA = 16/8): %.4f" % (float(nB) / nA))
    print("    Expected (alpha_W/alpha_EM ~ 1/sin2tW): %.4f" % (1.0 / SIN2_TW))

print()

# Distance-dependent amplitude
print("Distance-dependent amplitude (confinement test):")
print("  Compute |A(d)| as function of fermion separation d")
print()

for d_target in range(1, max_dist + 1):
    # Find fermion pairs at distance d
    pairs = []
    for i in fermion_vertices[:30]:
        for j in fermion_vertices[:30]:
            if i < j and dist_matrix[i, j] == d_target:
                pairs.append((i, j))

    if len(pairs) < 2:
        continue

    # Compute average self-energy-like amplitude: G_gauge(f,f') at distance d
    amp_A_d = []
    amp_B_d = []
    for f1, f2 in pairs[:min(20, len(pairs))]:
        # Via A-type
        a_A = sum(G_amp[f1, g] * G_amp[g, f2] for g in range(N) if types[g] == 'A')
        a_B = sum(G_amp[f1, g] * G_amp[g, f2] for g in range(N) if types[g] == 'B')
        amp_A_d.append(abs(a_A))
        amp_B_d.append(abs(a_B))

    mean_A = np.mean(amp_A_d) if amp_A_d else 0
    mean_B = np.mean(amp_B_d) if amp_B_d else 0
    ratio = mean_B / mean_A if mean_A > 0 else 0

    print("  d=%d: |A_EM| = %.6e, |A_weak| = %.6e, ratio = %.4f  (pairs: %d)" % (
        d_target, mean_A, mean_B, ratio, len(pairs[:min(20, len(pairs))])))

print()
print("  ANALYSIS: If ratio changes with distance => different running for EM vs weak")
print("  If ratio ~ constant = nB/nA => no differential running (vertex-transitive)")

print()

# ===========================================================================
# PART G: WHAT WOULD BE NEEDED FOR TRUE DYNAMICS
# ===========================================================================
print("=" * 72)
print("PART G: ASSESSMENT - WHAT THE 600-CELL CAN AND CANNOT DO")
print("=" * 72)
print()

print("POSITIVE RESULTS:")
print("  1. Green's function (propagator) is well-defined on the graph")
print("  2. The propagator decays with distance (exponentially for m > 0)")
print("  3. Spectral representation gives natural momentum-space propagator")
print("  4. Triangle counting gives natural 3-point vertices")
print("  5. Heat kernel is computable and shows spectral dimension flow")
print("  6. The spectral dimension flows from ~0 (IR) to maximum ~3 (intermediate)")
print()

print("NEGATIVE RESULTS:")
print("  1. VERTEX-TRANSITIVITY kills coupling running from random walks")
print("     (All couplings converge to n_X/N = trivial fractions)")
print("  2. The graph is FINITE (120 vertices) - no continuum limit")
print("  3. No natural separation of gauge sectors from graph structure alone")
print("     (A,B,C classification is geometric, not dynamical)")
print("  4. Beta functions cannot be reproduced from pure graph counting")
print("  5. No scattering matrix - only propagators (2-point functions)")
print("  6. The amplitude ratios are determined by vertex COUNTING, not dynamics")
print()

print("FUNDAMENTAL OBSTACLES:")
print("  1. The 600-cell is a SINGLE polytope - dynamics requires a FIELD over it")
print("  2. Graph propagator G(i,j) is a STATIC object, not a scattering amplitude")
print("  3. Running couplings require RENORMALIZATION GROUP, which needs")
print("     a continuous energy parameter - the graph has only 6 distance scales")
print("  4. The vertex-transitive symmetry ensures all vertices are equivalent,")
print("     so gauge vs fermion distinction must come from SYMMETRY BREAKING,")
print("     not from the graph structure itself")
print("  5. To get dynamics, one would need:")
print("     a) A field theory ON the graph (like lattice gauge theory)")
print("     b) A path integral over field configurations")
print("     c) A mechanism for symmetry breaking")
print("     d) A continuum limit (lattice spacing -> 0)")
print()

print("SPECULATIVE DIRECTIONS (for future work):")
print("  1. Lattice gauge theory on the 600-cell graph")
print("     - Place SU(3) link variables on edges, fermion fields on vertices")
print("     - Wilson action: S = sum_plaquettes Re Tr(U_p)")
print("     - The 1200 triangles serve as 'plaquettes'")
print("     - This would give genuine dynamics but is essentially standard LGT")
print()
print("  2. Spectral action approach")
print("     - S = Tr f(D/Lambda) where D is the Dirac operator on the graph")
print("     - Requires defining a Dirac operator (not just Laplacian)")
print("     - Noncommutative geometry framework (Connes)")
print()
print("  3. Topological field theory")
print("     - Use the 600-cell as a triangulation of S^3")
print("     - Ponzano-Regge or Turaev-Viro state sum models")
print("     - Would give 3D quantum gravity partition function")
print()
print("  4. Soliton scattering")
print("     - The theory has solitons (E_flip=3, E_pair=4, E_bag=6)")
print("     - Scattering of solitons on the graph COULD give dynamics")
print("     - Requires simulating soliton propagation (cellular automaton-like)")
print()

# ===========================================================================
# NUMERICAL SUMMARY
# ===========================================================================
print("=" * 72)
print("NUMERICAL SUMMARY")
print("=" * 72)
print()

print("600-cell graph properties:")
print("  Vertices: %d  Edges: %d  Degree: %d  Diameter: %d" % (N, total_edges, 12, max_dist))
print("  Eigenvalues: 9 distinct, all of form a + b*phi")
print("  Triangles: %d" % triangle_count)
print()

print("Green's function (m^2=1.0):")
for d in range(max_dist + 1):
    G_d = greens_data[1.0].get(d, 0)
    print("  G(d=%d) = %.6f" % (d, G_d))
print()

print("Heat kernel spectral dimension:")
print("  Peak d_s = %.3f at t = %.4f" % (d_spectral[idx_max_ds], t_values[idx_max_ds]))
print("  UV limit: d_s -> %.3f" % d_spectral[5])
print("  IR limit: d_s -> %.3f" % d_spectral[-5])
print()

print("Coupling running (random walk fractions at step 7):")
last_step = list(coupling_data.keys())[-1]
fA, fB, fC = coupling_data[last_step]
print("  f_A = %.5f (ergodic: %.5f)" % (fA, frac_A))
print("  f_B = %.5f (ergodic: %.5f)" % (fB, frac_B))
print("  f_C = %.5f (ergodic: %.5f)" % (fC, frac_C))
print("  Convergence to ergodic => NO running (vertex-transitive graph)")
print()

print("=" * 72)
print("CONCLUSION")
print("=" * 72)
print()
print("The 600-cell graph provides a natural PROPAGATOR (Green's function)")
print("and a natural VERTEX FUNCTION (from triangles), but it CANNOT produce")
print("genuine dynamical predictions (running couplings, scattering amplitudes,")
print("cross-sections) from graph structure alone.")
print()
print("The fundamental reason is VERTEX-TRANSITIVITY: all 120 vertices are")
print("equivalent under the symmetry group, so any vertex-based distinction")
print("(gauge vs fermion) is necessarily broken by hand, not by dynamics.")
print()
print("To get dynamics, one needs additional structure ON TOP of the graph:")
print("field theory, path integral, or soliton dynamics.")
print()
print("The heat kernel analysis reveals an interesting spectral dimension")
print("flow, but this is a property of ANY regular graph, not specific")
print("to the 600-cell's role as a Theory of Everything.")
print()
print("STATUS: NEGATIVE / EXPLORATORY")
print("  - Green's functions: DERIVED (from graph Laplacian)")
print("  - Scattering amplitudes: NOT POSSIBLE (from graph alone)")
print("  - Running couplings: NEGATIVE (vertex-transitivity prevents running)")
print("  - Beta functions: NEGATIVE (cannot reproduce from graph counting)")
print("  - Spectral dimension: COMPUTED (but generic to regular graphs)")
print("  - Path forward: field theory on the graph (future work)")
print()
print("EXP-196 COMPLETE")

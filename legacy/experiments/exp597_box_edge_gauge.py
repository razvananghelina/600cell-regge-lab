"""
exp100: Box operator on the EDGE space (720-dimensional) of the 600-cell.

Goal: Derive gauge boson spectrum from geometry.
Prediction: ker(Box_edge) should contain 12 = 1+3+8 massless modes
           decomposing as U(1) + SU(2) + SU(3).

Approach: Build multiple definitions of Box_edge, compute spectra BLINDLY,
          then analyze structure.

Follows the step-by-step plan in astazi_bro.txt.
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh, norm
from collections import defaultdict, Counter
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
TOL = 1e-8

# =====================================================================
# STEP 1: BUILD 600-CELL + SIMPLICIAL COMPLEX
# =====================================================================
print("=" * 72)
print("STEP 1: BUILD 600-CELL + SIMPLICIAL COMPLEX")
print("=" * 72)

verts, adj, lap = build_600cell()
Nv = len(verts)
print(f"  Vertices: {Nv} (expect 120)")
assert Nv == 120, f"Expected 120 vertices, got {Nv}"

# Verify degree
degrees = np.sum(adj, axis=1)
print(f"  Degree range: [{degrees.min():.0f}, {degrees.max():.0f}] (expect 12)")
assert np.allclose(degrees, 12), "Not all vertices have degree 12"

# Build adjacency list and edge list
adj_list = defaultdict(set)
edges = []
edge_to_idx = {}
for i in range(Nv):
    for j in range(i+1, Nv):
        if adj[i, j] > 0.5:
            adj_list[i].add(j)
            adj_list[j].add(i)
            idx = len(edges)
            edges.append((i, j))
            edge_to_idx[(i, j)] = idx

Ne = len(edges)
print(f"  Edges: {Ne} (expect 720)")
assert Ne == 720, f"Expected 720 edges, got {Ne}"

# Build triangles (faces)
triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    triangles.append((i, j, k))

Nf = len(triangles)
print(f"  Faces (triangles): {Nf} (expect 1200)")
assert Nf == 1200, f"Expected 1200 faces, got {Nf}"

# Euler characteristic
Nc = 600  # cells (tetrahedra)
euler = Nv - Ne + Nf - Nc
print(f"  Euler characteristic: {euler} (expect 0 for S^3)")
assert euler == 0, f"Euler check failed: {euler}"

# Build d0: vertices -> edges (Ne x Nv)
d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = +1.0

# Build d1: edges -> faces (Nf x Ne)
d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

# Verify d1 @ d0 = 0
err_d1d0 = norm(d1 @ d0)
print(f"  ||d1 @ d0|| = {err_d1d0:.2e} (expect 0)")
assert err_d1d0 < 1e-10, f"d1 @ d0 != 0: {err_d1d0}"

# Verify d0^T d0 = graph Laplacian
Delta_0 = d0.T @ d0
err_lap = norm(Delta_0 - lap)
print(f"  ||d0^T d0 - L_graph|| = {err_lap:.2e} (expect 0)")
assert err_lap < 1e-10, f"Delta_0 != L_graph: {err_lap}"

print("\n  STEP 1 PASSED: All verifications OK.")


# =====================================================================
# STEP 2: HOPF FIBRATION + EDGE CLASSIFICATION
# =====================================================================
print("\n" + "=" * 72)
print("STEP 2: HOPF FIBRATION + EDGE CLASSIFICATION")
print("=" * 72)

def qmul(p, q):
    """Quaternion multiplication."""
    w = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]
    x = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2]
    y = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1]
    z = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0]
    return np.array([w, x, y, z])

def find_idx(v, verts, tol=1e-6):
    """Find vertex index closest to v."""
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

def find_fibration():
    """Find ordered Hopf fibration: 12 fibers of 10 vertices."""
    target_w = PHI / 2.0
    for i in range(Nv):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            p = g.copy()
            ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1, 0, 0, 0], atol=1e-6):
                    ok = False; break
                if k == 10 and not np.allclose(p, [1, 0, 0, 0], atol=1e-6):
                    ok = False
            if not ok:
                continue
            used = set()
            fibers = []
            subg = []
            pp = np.array([1.0, 0, 0, 0])
            for k in range(10):
                idx = find_idx(pp, verts)
                subg.append(idx)
                pp = qmul(pp, g)
            for s in range(Nv):
                if s in used:
                    continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si])
                    idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used:
                        fib.append(idx)
                        used.add(idx)
                if len(fib) == 10:
                    fibers.append(fib)
            if len(fibers) == 12:
                return fibers
    return None

fibers = find_fibration()
assert fibers is not None, "Could not find Hopf fibration"
assert len(fibers) == 12, f"Expected 12 fibers, got {len(fibers)}"

# Verify: each fiber has 10 vertices, all 120 vertices covered
all_verts_in_fibers = set()
for f in fibers:
    assert len(f) == 10, f"Fiber has {len(f)} vertices, expected 10"
    all_verts_in_fibers.update(f)
assert len(all_verts_in_fibers) == 120, "Fibers don't cover all vertices"
print(f"  Hopf fibration: {len(fibers)} fibers x {len(fibers[0])} vertices = {sum(len(f) for f in fibers)}")

# Build fiber lookup
vertex_to_fiber = {}
for fi, f in enumerate(fibers):
    for v in f:
        vertex_to_fiber[v] = fi

# Classify edges
is_fiber_edge = np.zeros(Ne, dtype=bool)
for e_idx, (i, j) in enumerate(edges):
    if vertex_to_fiber[i] == vertex_to_fiber[j]:
        is_fiber_edge[e_idx] = True

n_fiber = np.sum(is_fiber_edge)
n_cross = Ne - n_fiber
print(f"  Fiber edges: {n_fiber} (expect 120)")
print(f"  Cross edges: {n_cross} (expect 600)")
assert n_fiber == 120, f"Expected 120 fiber edges, got {n_fiber}"
assert n_cross == 600, f"Expected 600 cross edges, got {n_cross}"

# Build projectors
P_fiber = np.diag(is_fiber_edge.astype(float))
P_cross = np.eye(Ne) - P_fiber

# Verify fiber edges form 12 decagonal cycles
# Each fiber should have exactly 10 edges connecting consecutive vertices
fiber_edge_count = Counter()
for e_idx, (i, j) in enumerate(edges):
    if is_fiber_edge[e_idx]:
        fi = vertex_to_fiber[i]
        fiber_edge_count[fi] += 1

for fi in range(12):
    assert fiber_edge_count[fi] == 10, f"Fiber {fi} has {fiber_edge_count[fi]} edges, expected 10"
print(f"  Each fiber has exactly 10 edges (decagonal cycles verified)")

print("\n  STEP 2 PASSED: Hopf fibration and edge classification OK.")


# =====================================================================
# STEP 3: HODGE LAPLACIAN
# =====================================================================
print("\n" + "=" * 72)
print("STEP 3: HODGE LAPLACIAN")
print("=" * 72)

B = d0 @ d0.T       # Exact part (720 x 720)
C = d1.T @ d1       # Coexact part (720 x 720)
Delta_1 = B + C      # Full edge Laplacian

# Verify B*C = 0 (orthogonality of exact and coexact)
err_BC = norm(B @ C, 'fro')
print(f"  ||B @ C||_F = {err_BC:.2e} (expect 0)")
assert err_BC < 1e-8, f"B and C not orthogonal: {err_BC}"

# Compute dimensions
evals_B = eigvalsh(B)
evals_C = eigvalsh(C)
dim_exact = np.sum(evals_B > TOL)
dim_coexact = np.sum(evals_C > TOL)
dim_harmonic = Ne - dim_exact - dim_coexact

print(f"  Exact dim (rank B): {dim_exact} (expect 119 = Nv-1)")
print(f"  Coexact dim (rank C): {dim_coexact} (expect 601)")
print(f"  Harmonic dim (b1 of S^3): {dim_harmonic} (expect 0)")
assert dim_exact == 119, f"Expected 119, got {dim_exact}"
assert dim_coexact == 601, f"Expected 601, got {dim_coexact}"
assert dim_harmonic == 0, f"Expected 0 harmonic, got {dim_harmonic}"

# Verify: nonzero spectrum of B = nonzero spectrum of Delta_0
evals_Delta0 = eigvalsh(Delta_0)
nonzero_B = np.sort(evals_B[evals_B > TOL])
nonzero_D0 = np.sort(evals_Delta0[evals_Delta0 > TOL])
err_hodge = norm(nonzero_B - nonzero_D0)
print(f"  ||spec(B)_nz - spec(Delta_0)_nz|| = {err_hodge:.2e} (Hodge isomorphism)")
assert err_hodge < 1e-8, f"Hodge isomorphism failed: {err_hodge}"

# Commutativity check: does [C, P_fiber] = 0?
comm_CP = C @ P_fiber - P_fiber @ C
norm_CP = norm(comm_CP, 'fro')
comm_BP = B @ P_fiber - P_fiber @ B
norm_BP = norm(comm_BP, 'fro')
comm_D1P = Delta_1 @ P_fiber - P_fiber @ Delta_1
norm_D1P = norm(comm_D1P, 'fro')

print(f"\n  Commutativity with P_fiber:")
print(f"    ||[B, P_fiber]||_F = {norm_BP:.6f}")
print(f"    ||[C, P_fiber]||_F = {norm_CP:.6f}")
print(f"    ||[Delta_1, P_fiber]||_F = {norm_D1P:.6f}")
print(f"    B commutes: {norm_BP < TOL}")
print(f"    C commutes: {norm_CP < TOL}")
print(f"    Delta_1 commutes: {norm_D1P < TOL}")

# Classify triangles
n_tri_type = Counter()
for (i, j, k) in triangles:
    n_fib = sum([
        vertex_to_fiber[i] == vertex_to_fiber[j],
        vertex_to_fiber[j] == vertex_to_fiber[k],
        vertex_to_fiber[i] == vertex_to_fiber[k],
    ])
    n_tri_type[n_fib] += 1

print(f"\n  Triangle classification (fiber edges per triangle):")
for n_fib in sorted(n_tri_type):
    print(f"    {n_fib} fiber edges: {n_tri_type[n_fib]} triangles")

print("\n  STEP 3 COMPLETE.")


# =====================================================================
# STEP 4: BUILD LINE GRAPH L(G)
# =====================================================================
print("\n" + "=" * 72)
print("STEP 4: LINE GRAPH L(G)")
print("=" * 72)

# L(G) has 720 nodes (one per edge of G)
# Two nodes in L(G) are adjacent if the corresponding edges share a vertex

# Build adjacency matrix of line graph
print("  Building line graph adjacency (720 x 720)...")

# For each vertex, collect all edges incident to it
vertex_to_edges = defaultdict(list)
for e_idx, (i, j) in enumerate(edges):
    vertex_to_edges[i].append(e_idx)
    vertex_to_edges[j].append(e_idx)

# Build line graph adjacency
A_line = np.zeros((Ne, Ne), dtype=float)
for v in range(Nv):
    inc_edges = vertex_to_edges[v]
    for a in range(len(inc_edges)):
        for b in range(a+1, len(inc_edges)):
            A_line[inc_edges[a], inc_edges[b]] = 1.0
            A_line[inc_edges[b], inc_edges[a]] = 1.0

# Verify degree of line graph
line_degrees = np.sum(A_line, axis=1)
unique_degrees = Counter(line_degrees.astype(int))
print(f"  Line graph degrees: {dict(unique_degrees)}")
print(f"  Min degree: {line_degrees.min():.0f}, Max degree: {line_degrees.max():.0f}")

# For 600-cell (d=12 regular), each edge has degree d(i)-1 + d(j)-1 - 2*t_edge
# where t_edge = triangles containing the edge
# Expected: 11 + 11 - 2*t = 22 - 2*t (WAIT: line graph formula is 2*(d-1) - t)
# Actually: line graph degree = (d_i - 1) + (d_j - 1) - (edges between N(i)&N(j) that share both endpoints)
# For vertex-transitive graph: L(G) degree = 2*(d-1) - (common neighbors that form triangle)
# Let me just count triangles per edge empirically

triangles_per_edge = np.zeros(Ne, dtype=int)
for (i, j, k) in triangles:
    triangles_per_edge[edge_to_idx[(i, j)]] += 1
    triangles_per_edge[edge_to_idx[(j, k)]] += 1
    triangles_per_edge[edge_to_idx[(i, k)]] += 1

unique_tpe = Counter(triangles_per_edge)
print(f"  Triangles per edge: {dict(unique_tpe)}")

# Line graph degree should be 2*(d-1) - t_edge for each edge
expected_line_deg = 2 * (12 - 1) - triangles_per_edge
err_line_deg = norm(line_degrees - expected_line_deg)
print(f"  ||actual_deg - (22 - t_edge)|| = {err_line_deg:.2e}")
print(f"  Expected line graph degree: {22 - triangles_per_edge[0]} (if all edges have same t)")

# Build fiber adjacency in line graph
# A_fiber_line: two edges of G are fiber-adjacent if they share a FIBER vertex
# AND are in the same fiber
print("\n  Building fiber adjacency in line graph...")

# Method: edges in same fiber that share a vertex (consecutive in C10)
A_fiber_line = np.zeros((Ne, Ne), dtype=float)

for fi, fiber in enumerate(fibers):
    # Get all edges in this fiber
    fiber_edges_in_fiber = []
    for k in range(10):
        u, v = fiber[k], fiber[(k+1) % 10]
        i, j = (min(u, v), max(u, v))
        e_idx = edge_to_idx[(i, j)]
        fiber_edges_in_fiber.append(e_idx)

    # Connect consecutive fiber edges (they share a vertex)
    for k in range(10):
        e1 = fiber_edges_in_fiber[k]
        e2 = fiber_edges_in_fiber[(k+1) % 10]
        A_fiber_line[e1, e2] = 1.0
        A_fiber_line[e2, e1] = 1.0

fiber_line_degrees = np.sum(A_fiber_line, axis=1)
print(f"  Fiber line graph: fiber edges have degree {fiber_line_degrees[is_fiber_edge].mean():.1f}")
print(f"  Fiber line graph: cross edges have degree {fiber_line_degrees[~is_fiber_edge].mean():.1f}")
print(f"  (Expect: fiber edges = 2, cross edges = 0)")

# Cross adjacency in line graph
A_cross_line = A_line - A_fiber_line

# Verify A_cross_line is non-negative
min_cross = A_cross_line.min()
print(f"  Min of A_cross_line: {min_cross:.1f} (expect >= 0)")

# Verify it's symmetric
sym_err = norm(A_cross_line - A_cross_line.T)
print(f"  ||A_cross_line - A_cross_line^T|| = {sym_err:.2e}")

print("\n  STEP 4 COMPLETE.")


# =====================================================================
# STEP 5: DEFINE AND COMPUTE BOX_EDGE SPECTRA
# =====================================================================
print("\n" + "=" * 72)
print("STEP 5: BOX_EDGE DEFINITIONS AND SPECTRA")
print("=" * 72)

# Laplacian of line graph
D_line = np.diag(line_degrees)
L_line = D_line - A_line

# Fiber Laplacian of line graph
D_fiber_line = np.diag(np.sum(A_fiber_line, axis=1))
L_fiber_line = D_fiber_line - A_fiber_line

# Cross Laplacian of line graph
D_cross_line = np.diag(np.sum(A_cross_line, axis=1))
L_cross_line = D_cross_line - A_cross_line

# ------- Definition 1: Analog of vertex Box on line graph -------
# Box_vertex = b1 * A_fiber - A = b1 * A_fiber - A_fiber - A_cross = (b1-1)*A_fiber - A_cross
# Or equivalently: Box = L_cross - a1 * L_fiber (on vertices)
# On line graph: Box1_line = L_cross_line - a1 * L_fiber_line
print("\n--- Definition 1: Box1 = L_cross_line - a1 * L_fiber_line ---")
Box1 = L_cross_line - a1 * L_fiber_line
evals_Box1 = np.sort(eigvalsh(Box1))
ker1 = np.sum(np.abs(evals_Box1) < 1e-6)
print(f"  ker(Box1) = {ker1}")
print(f"  Tr(Box1) = {np.trace(Box1):.6f}")
print(f"  min(evals) = {evals_Box1[0]:.6f}")
print(f"  max(evals) = {evals_Box1[-1]:.6f}")

# Distinct eigenvalues
ev1_distinct = Counter(np.round(evals_Box1, 4))
print(f"  Distinct eigenvalues: {len(ev1_distinct)}")
print(f"  {'eigenvalue':>12s} {'mult':>5s}")
for val, mult in sorted(ev1_distinct.items()):
    print(f"  {val:12.4f} {mult:5d}")


# ------- Definition 2: Hodge version -------
# Project Delta_1 via fiber/cross
# Delta_1_fiber = P_fiber @ Delta_1 @ P_fiber (restricted to fiber subspace)
# But Delta_1 doesn't commute with P_fiber in general, so use block approach
# Box2 = P_cross @ Delta_1 @ P_cross - a1 * P_fiber @ Delta_1 @ P_fiber
print("\n--- Definition 2: Box2 = P_cross @ Delta_1 @ P_cross - a1 * P_fiber @ Delta_1 @ P_fiber ---")
Delta_1_ff = P_fiber @ Delta_1 @ P_fiber
Delta_1_cc = P_cross @ Delta_1 @ P_cross
Box2 = Delta_1_cc - a1 * Delta_1_ff
evals_Box2 = np.sort(eigvalsh(Box2))
ker2 = np.sum(np.abs(evals_Box2) < 1e-6)
print(f"  ker(Box2) = {ker2}")
print(f"  Tr(Box2) = {np.trace(Box2):.6f}")
print(f"  min(evals) = {evals_Box2[0]:.6f}")
print(f"  max(evals) = {evals_Box2[-1]:.6f}")

ev2_distinct = Counter(np.round(evals_Box2, 4))
print(f"  Distinct eigenvalues: {len(ev2_distinct)}")
print(f"  {'eigenvalue':>12s} {'mult':>5s}")
for val, mult in sorted(ev2_distinct.items()):
    print(f"  {val:12.4f} {mult:5d}")


# ------- Definition 3: b1 * A_fiber_line - A_line (direct analog of vertex Box) -------
print("\n--- Definition 3: Box3 = b1 * A_fiber_line - A_line ---")
Box3 = b1 * A_fiber_line - A_line
evals_Box3 = np.sort(eigvalsh(Box3))
ker3 = np.sum(np.abs(evals_Box3) < 1e-6)
print(f"  ker(Box3) = {ker3}")
print(f"  Tr(Box3) = {np.trace(Box3):.6f}")
print(f"  min(evals) = {evals_Box3[0]:.6f}")
print(f"  max(evals) = {evals_Box3[-1]:.6f}")

ev3_distinct = Counter(np.round(evals_Box3, 4))
print(f"  Distinct eigenvalues: {len(ev3_distinct)}")
print(f"  {'eigenvalue':>12s} {'mult':>5s}")
for val, mult in sorted(ev3_distinct.items()):
    print(f"  {val:12.4f} {mult:5d}")


# ------- Definition 4: Coexact Box (from exp564 insight) -------
# C = d1^T d1 (coexact Laplacian on edges)
# Box4 = C - a1 * B = coexact - a1 * exact (signed edge Laplacian)
print("\n--- Definition 4: Box4 = C - a1 * B (Hodge signed) ---")
Box4 = C - a1 * B
evals_Box4 = np.sort(eigvalsh(Box4))
ker4 = np.sum(np.abs(evals_Box4) < 1e-6)
print(f"  ker(Box4) = {ker4}")
print(f"  Tr(Box4) = {np.trace(Box4):.6f}")
print(f"  min(evals) = {evals_Box4[0]:.6f}")
print(f"  max(evals) = {evals_Box4[-1]:.6f}")

ev4_distinct = Counter(np.round(evals_Box4, 4))
print(f"  Distinct eigenvalues: {len(ev4_distinct)}")
print(f"  {'eigenvalue':>12s} {'mult':>5s}")
for val, mult in sorted(ev4_distinct.items()):
    print(f"  {val:12.4f} {mult:5d}")


# ------- Definition 5: Line graph with full fiber direction -------
# A_fiber_full: two edges of G are "fiber-connected" if they share a vertex
# AND that shared vertex's fiber-neighbor connects them in the fiber direction
# I.e., the shared vertex v has the two edges "going along the fiber through v"
print("\n--- Definition 5: Box5 = b1 * A_fiber_full - A_line ---")
# For each vertex v, look at all incident edges. Among them, exactly 2 are fiber edges.
# For cross edges incident to v, consider them fiber-connected if they share
# the same "base point" on S^2 projection. But that's too speculative.
# Instead: A_fiber_full connects edge e1 and e2 if they share vertex v,
# and BOTH e1 and e2 are fiber edges (in the same fiber through v).
# This is the same as A_fiber_line above.

# Alternative: connect ANY two edges that share a vertex IF that vertex
# contributes a fiber direction. I.e., for vertex v with fiber edges e_f1, e_f2:
# connect e_f1 and e_f2 to ALL other edges incident to v.
# This gives a "star" pattern per vertex in fiber direction.

A_fiber_star = np.zeros((Ne, Ne), dtype=float)
for v in range(Nv):
    inc = vertex_to_edges[v]
    # Find which of the incident edges are fiber edges
    fiber_inc = [e for e in inc if is_fiber_edge[e]]
    # Connect each fiber edge to all other incident edges at this vertex
    for ef in fiber_inc:
        for e_other in inc:
            if e_other != ef:
                A_fiber_star[ef, e_other] = 1.0
                A_fiber_star[e_other, ef] = 1.0

# Make sure diagonal is zero
np.fill_diagonal(A_fiber_star, 0)
# Clip to binary
A_fiber_star = np.minimum(A_fiber_star, 1.0)

Box5 = b1 * A_fiber_star - A_line
evals_Box5 = np.sort(eigvalsh(Box5))
ker5 = np.sum(np.abs(evals_Box5) < 1e-6)
print(f"  ker(Box5) = {ker5}")
print(f"  Tr(Box5) = {np.trace(Box5):.6f}")
print(f"  min(evals) = {evals_Box5[0]:.6f}")
print(f"  max(evals) = {evals_Box5[-1]:.6f}")

ev5_distinct = Counter(np.round(evals_Box5, 4))
print(f"  Distinct eigenvalues: {len(ev5_distinct)}")
print(f"  {'eigenvalue':>12s} {'mult':>5s}")
for val, mult in sorted(ev5_distinct.items()):
    print(f"  {val:12.4f} {mult:5d}")


# =====================================================================
# STEP 6: KERNEL ANALYSIS
# =====================================================================
print("\n" + "=" * 72)
print("STEP 6: KERNEL ANALYSIS")
print("=" * 72)

# For each definition with nontrivial kernel, analyze it
definitions = [
    ("Box1 (L_cross - a1*L_fiber on line graph)", Box1, evals_Box1),
    ("Box2 (Hodge projected)", Box2, evals_Box2),
    ("Box3 (b1*A_fiber - A on line graph)", Box3, evals_Box3),
    ("Box4 (C - a1*B Hodge signed)", Box4, evals_Box4),
    ("Box5 (b1*A_fiber_star - A)", Box5, evals_Box5),
]

for name, op, evals in definitions:
    ker_dim = np.sum(np.abs(evals) < 1e-6)
    print(f"\n  {name}:")
    print(f"    ker = {ker_dim}")

    if ker_dim > 0 and ker_dim <= 200:
        # Get kernel eigenvectors
        all_evals, all_evecs = eigh(op)
        ker_mask = np.abs(all_evals) < 1e-6
        ker_vecs = all_evecs[:, ker_mask]

        # Fiber/cross content of kernel
        fiber_content = np.sum(ker_vecs[is_fiber_edge, :]**2) / ker_dim
        cross_content = np.sum(ker_vecs[~is_fiber_edge, :]**2) / ker_dim
        print(f"    Average fiber content: {fiber_content:.4f}")
        print(f"    Average cross content: {cross_content:.4f}")
        print(f"    Fiber/cross ratio: {fiber_content/max(cross_content, 1e-15):.4f}")

        # Per-eigenvector fiber content
        print(f"    Per-eigenvector fiber fractions:")
        for k in range(ker_dim):
            v = ker_vecs[:, k]
            ff = np.sum(v[is_fiber_edge]**2)
            print(f"      mode {k}: fiber_frac = {ff:.4f}")

    # First few eigenvalues near zero (potential massive modes)
    near_zero = evals[np.abs(evals) < 5.0]
    if len(near_zero) > 0:
        ev_near = Counter(np.round(near_zero, 4))
        print(f"    Near-zero eigenvalues (|lambda| < 5):")
        for val, mult in sorted(ev_near.items()):
            print(f"      {val:10.4f} x {mult}")


# =====================================================================
# STEP 7: CROSS-CHECKS
# =====================================================================
print("\n" + "=" * 72)
print("STEP 7: CROSS-CHECKS")
print("=" * 72)

for name, op, evals in definitions:
    print(f"\n  {name}:")

    # Tracelessness
    tr = np.trace(op)
    print(f"    Tr = {tr:.6f}")

    # Trace moments
    tr2 = np.trace(op @ op)
    tr3 = np.trace(op @ op @ op)
    tr4 = np.trace(op @ op @ op @ op)
    print(f"    Tr(Box^2) = {tr2:.4f}")
    print(f"    Tr(Box^3) = {tr3:.4f}")
    print(f"    Tr(Box^4) = {tr4:.4f}")

    # Check if eigenvalues are in Z[phi]
    # An eigenvalue is in Z[phi] if it's of the form a + b*phi with a,b integer
    # Equivalent: lambda and (lambda - lambda_conjugate)/sqrt(5) are both in Z
    # Or: lambda^2 - (integer)*lambda + (integer) = 0
    # For each distinct eigenvalue, check if it satisfies a quadratic with integer coeffs
    n_zphi = 0
    n_total = 0
    distinct_evals = sorted(set(np.round(evals, 6)))
    for lam in distinct_evals:
        if abs(lam) < 1e-8:
            n_zphi += 1; n_total += 1
            continue
        # Check: is lam in Z[phi]?
        # lam = a + b*phi => lam satisfies x^2 - (2a+b)x + (a^2+ab-b^2) = 0
        # The trace (2a+b) and norm (a^2+ab-b^2) should be integers
        # Galois conjugate: lam' = a + b*(1-phi) = a+b - b*phi = (a+b) - b*phi
        # Sum: lam + lam' = 2a + b (integer)
        # Product: lam * lam' = a^2 + ab - b^2 (integer)
        # So check: lam + Galois_conj = integer, lam * Galois_conj = integer
        # But we don't know the conjugate. Instead: solve for a,b.
        # lam = a + b*phi => b = (lam - round(lam - round(lam)*PHI_frac))...
        # Simpler: b = round((lam - round(lam))/PHI) won't work.
        # Just check if lam is close to a + b*phi for small a,b
        found = False
        for b in range(-30, 31):
            a_test = lam - b * PHI
            if abs(a_test - round(a_test)) < 0.001:
                n_zphi += 1
                found = True
                break
        n_total += 1

    print(f"    Eigenvalues in Z[phi]: {n_zphi}/{n_total}")

    # Ratio with vertex Box
    # Vertex Tr(Box^2) = sum of eigenvalues^2 on vertex space
    Box_vertex = b1 * np.diag(np.sum(adj * np.array([[float(vertex_to_fiber[i] == vertex_to_fiber[j])
                     for j in range(Nv)] for i in range(Nv)]), axis=1)) \
                 - b1 * adj * np.array([[float(vertex_to_fiber[i] == vertex_to_fiber[j])
                     for j in range(Nv)] for i in range(Nv)])\
                 + np.diag(np.sum(adj * np.array([[float(vertex_to_fiber[i] != vertex_to_fiber[j])
                     for j in range(Nv)] for i in range(Nv)]), axis=1)) \
                 - adj * np.array([[float(vertex_to_fiber[i] != vertex_to_fiber[j])
                     for j in range(Nv)] for i in range(Nv)])
    # Wait, simpler: Box_vertex = L_cross - a1 * L_fiber
    A_fiber_v = adj * np.array([[float(vertex_to_fiber[i] == vertex_to_fiber[j])
                    for j in range(Nv)] for i in range(Nv)])
    A_cross_v = adj - A_fiber_v
    L_fiber_v = np.diag(A_fiber_v.sum(axis=1)) - A_fiber_v
    L_cross_v = np.diag(A_cross_v.sum(axis=1)) - A_cross_v
    Box_v = L_cross_v - a1 * L_fiber_v
    tr2_v = np.trace(Box_v @ Box_v)

    if abs(tr2) > 0.01:
        ratio = tr2 / tr2_v
        print(f"    Tr(Box_edge^2) / Tr(Box_vertex^2) = {ratio:.6f}")
        # Check if ratio is a framework constant
        for label, val in [("a1", 5), ("b1", 6), ("N/Nv", Ne/Nv), ("phi", PHI),
                           ("phi^2", PHI**2), ("a1*b1", 30), ("12", 12)]:
            if abs(ratio - val) < 0.01:
                print(f"      *** Matches {label} = {val} ***")


# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

print(f"\n  Definition                                  ker    Tr    min_eval  max_eval")
print(f"  {'-'*80}")
for name, op, evals in definitions:
    ker_dim = np.sum(np.abs(evals) < 1e-6)
    tr = np.trace(op)
    print(f"  {name:45s} {ker_dim:4d} {tr:8.2f} {evals[0]:9.4f} {evals[-1]:9.4f}")

print(f"\n  Target: ker = 12 = 1 + 3 + 8 (gauge bosons)")
print(f"  Target: ker = 9 = 1 + 8 (if SU(2) broken)")

# Check if any definition gives 9 or 12
for name, op, evals in definitions:
    ker_dim = np.sum(np.abs(evals) < 1e-6)
    if ker_dim in [9, 12]:
        print(f"\n  *** {name} has ker = {ker_dim} ***")

print("\n" + "=" * 72)
print("EXP100 COMPLETE")
print("=" * 72)

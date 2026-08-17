"""
exp564: Lorentzian operator Box on the EDGE space (720-dimensional).

The vertex Box acts on scalars (120 modes). W/Z are vector bosons -- they
live in the 1-form (edge) sector. The coexact Laplacian C = d1^T @ d1
acts on 720-dimensional edge space with 601 nonzero eigenvalues.

PLAN:
  1. Build full simplicial complex (verts, edges, faces, cells)
  2. Build boundary operators d0, d1
  3. Build Hodge decomposition: B = d0 @ d0^T, C = d1^T @ d1
  4. Classify edges as fiber (120) vs cross-fiber (600)
  5. Build fiber/cross projectors on edge space
  6. Check if [C, P_fiber] = 0 (do they commute?)
  7. If yes: decompose C into C_fiber + C_cross, build Box_edge
  8. Analyze spectrum for gauge boson structure
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh, norm
from collections import defaultdict
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
TOL = 1e-8

print("Building 600-cell...")
verts, adj, lap = build_600cell()
Nv = 120

# === Build Hopf fibration ===
def qmul(p, q):
    w = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]
    x = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2]
    y = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1]
    z = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0]
    return np.array([w, x, y, z])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

def find_fibration():
    target_w = PHI / 2.0
    for i in range(Nv):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            p = g.copy()
            ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6):
                    ok = False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6):
                    ok = False
            if not ok: continue
            used = set()
            fibers = []
            subg = []
            pp = np.array([1.0,0,0,0])
            for k in range(10):
                idx = find_idx(pp, verts)
                subg.append(idx); pp = qmul(pp, g)
            for s in range(Nv):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si])
                    idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used:
                        fib.append(idx); used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12: return fibers
    return None

fibers = find_fibration()
fiber_sets = [set(f) for f in fibers]
vertex_to_fiber = {}
for fi, f in enumerate(fibers):
    for v in f:
        vertex_to_fiber[v] = fi


# =====================================================================
# STEP 1: Build simplicial complex
# =====================================================================
print("\n" + "=" * 70)
print("STEP 1: SIMPLICIAL COMPLEX")
print("=" * 70)

# Edges
adj_list = defaultdict(set)
edges = []
edge_to_idx = {}
for i in range(Nv):
    for j in range(i+1, Nv):
        if adj[i,j] > 0.5:
            adj_list[i].add(j)
            adj_list[j].add(i)
            idx = len(edges)
            edges.append((i, j))
            edge_to_idx[(i,j)] = idx

Ne = len(edges)
print(f"  Vertices: {Nv}, Edges: {Ne}")

# Faces (triangles)
triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    triangles.append((i, j, k))

face_to_idx = {t: idx for idx, t in enumerate(triangles)}
Nf = len(triangles)
print(f"  Faces: {Nf}, Cells (expected): 600")
print(f"  Euler: {Nv - Ne + Nf - 600}")


# =====================================================================
# STEP 2: Boundary operators
# =====================================================================
print("\n" + "=" * 70)
print("STEP 2: BOUNDARY OPERATORS")
print("=" * 70)

# d0: vertices -> edges (Ne x Nv)
d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = +1.0

# d1: edges -> faces (Nf x Ne)
d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    # boundary of (i,j,k) = (j,k) - (i,k) + (i,j)
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

# Verify d1 @ d0 = 0
err_d1d0 = norm(d1 @ d0)
print(f"  ||d1 @ d0|| = {err_d1d0:.2e} (expect 0)")


# =====================================================================
# STEP 3: Hodge Laplacians
# =====================================================================
print("\n" + "=" * 70)
print("STEP 3: HODGE LAPLACIANS")
print("=" * 70)

B = d0 @ d0.T       # Exact part of edge Laplacian (720 x 720)
C = d1.T @ d1       # Coexact part of edge Laplacian (720 x 720)
Delta_1 = B + C      # Full edge Laplacian

# Verify Delta_0 matches known graph Laplacian
Delta_0 = d0.T @ d0
err_lap = norm(Delta_0 - lap)
print(f"  ||Delta_0 - graph Laplacian|| = {err_lap:.2e}")

# Edge Laplacian spectrum dimensions
evals_B = eigvalsh(B)
evals_C = eigvalsh(C)
dim_exact = np.sum(evals_B > TOL)
dim_coexact = np.sum(evals_C > TOL)
dim_harmonic = Ne - dim_exact - dim_coexact
print(f"  Exact dim (image d0): {dim_exact} (expect 119 = Nv-1)")
print(f"  Coexact dim (image d1^T): {dim_coexact} (expect 601)")
print(f"  Harmonic dim (b1): {dim_harmonic} (expect 0)")


# =====================================================================
# STEP 4: Classify edges as fiber vs cross
# =====================================================================
print("\n" + "=" * 70)
print("STEP 4: FIBER vs CROSS EDGE CLASSIFICATION")
print("=" * 70)

is_fiber_edge = np.zeros(Ne, dtype=bool)
for e_idx, (i, j) in enumerate(edges):
    if vertex_to_fiber[i] == vertex_to_fiber[j]:
        is_fiber_edge[e_idx] = True

n_fiber_edges = np.sum(is_fiber_edge)
n_cross_edges = Ne - n_fiber_edges
print(f"  Fiber edges: {n_fiber_edges} (expect 120)")
print(f"  Cross edges: {n_cross_edges} (expect 600)")

# Build projectors
P_fiber = np.diag(is_fiber_edge.astype(float))      # 720 x 720
P_cross = np.eye(Ne) - P_fiber                       # 720 x 720


# =====================================================================
# STEP 5: Commutativity tests
# =====================================================================
print("\n" + "=" * 70)
print("STEP 5: COMMUTATIVITY TESTS")
print("=" * 70)

# Does the coexact Laplacian C commute with the fiber projector?
comm_CP = C @ P_fiber - P_fiber @ C
norm_CP = norm(comm_CP, 'fro')
print(f"  ||[C, P_fiber]||_F = {norm_CP:.6f}")
print(f"  C commutes with P_fiber: {norm_CP < TOL}")

# Does B commute with P_fiber?
comm_BP = B @ P_fiber - P_fiber @ B
norm_BP = norm(comm_BP, 'fro')
print(f"  ||[B, P_fiber]||_F = {norm_BP:.6f}")
print(f"  B commutes with P_fiber: {norm_BP < TOL}")

# Does Delta_1 commute with P_fiber?
comm_D1P = Delta_1 @ P_fiber - P_fiber @ Delta_1
norm_D1P = norm(comm_D1P, 'fro')
print(f"  ||[Delta_1, P_fiber]||_F = {norm_D1P:.6f}")
print(f"  Delta_1 commutes with P_fiber: {norm_D1P < TOL}")

# Classify triangles by fiber/cross content
n_tri_type = {0: 0, 1: 0, 2: 0, 3: 0}
for f_idx, (i, j, k) in enumerate(triangles):
    n_fib = sum([
        vertex_to_fiber[i] == vertex_to_fiber[j],
        vertex_to_fiber[j] == vertex_to_fiber[k],
        vertex_to_fiber[i] == vertex_to_fiber[k],
    ])
    n_tri_type[n_fib] += 1

print(f"\n  Triangle classification (fiber edges per triangle):")
for n_fib, count in sorted(n_tri_type.items()):
    print(f"    {n_fib} fiber edges: {count} triangles")
print(f"  Total: {sum(n_tri_type.values())} (expect {Nf})")


# =====================================================================
# STEP 6: Restricted operators (even if they don't commute)
# =====================================================================
print("\n" + "=" * 70)
print("STEP 6: RESTRICTED OPERATORS ON EDGE SUBSPACES")
print("=" * 70)

# Project C onto fiber and cross subspaces
C_ff = P_fiber @ C @ P_fiber   # fiber-fiber block (120x120 effective)
C_cc = P_cross @ C @ P_cross   # cross-cross block (600x600 effective)
C_fc = P_fiber @ C @ P_cross   # fiber-cross mixing
C_cf = P_cross @ C @ P_fiber   # cross-fiber mixing

norm_mixing = norm(C_fc, 'fro')
print(f"  ||C_fc|| (fiber-cross mixing) = {norm_mixing:.6f}")
if norm_mixing < TOL:
    print(f"  C is block-diagonal in fiber/cross! Can define C_fiber and C_cross.")
else:
    print(f"  C has fiber-cross mixing (expected from triangles with mixed edges).")
    print(f"  Mixing relative to ||C||: {norm_mixing/norm(C, 'fro'):.4f}")

# Even with mixing, compute the restricted spectra
evals_Cff = eigvalsh(C_ff)
evals_Ccc = eigvalsh(C_cc)

n_Cff_nonzero = np.sum(np.abs(evals_Cff) > TOL)
n_Ccc_nonzero = np.sum(np.abs(evals_Ccc) > TOL)
print(f"\n  C_ff (fiber-fiber block): {n_Cff_nonzero} nonzero eigenvalues out of {Ne}")
print(f"  C_cc (cross-cross block): {n_Ccc_nonzero} nonzero eigenvalues out of {Ne}")

# Alternative: build vertex-space fiber/cross Laplacians and lift to edges
# The vertex Laplacian L_fiber has eigenvectors. The d0 operator maps these
# to the edge space. So we can build B_fiber = d0 @ L_fiber_vertex... no,
# that's not quite right.

# Actually: B = d0 @ d0.T. And d0.T @ d0 = Delta_0 = L_vertex.
# The fiber part of the vertex Laplacian is L_fiber_v = 2I - A_fiber_v.
# Can we build B_fiber = d0_fiber @ d0_fiber.T where d0_fiber only uses
# fiber edges?

# d0 restricted to fiber edges
d0_fiber = d0[is_fiber_edge, :]     # (120 x 120)
d0_cross = d0[~is_fiber_edge, :]    # (600 x 120)

# B_fiber and B_cross from restricted d0
B_fiber_small = d0_fiber @ d0_fiber.T   # (120 x 120)
B_cross_small = d0_cross @ d0_cross.T   # (600 x 600)

print(f"\n  d0_fiber shape: {d0_fiber.shape}")
print(f"  d0_cross shape: {d0_cross.shape}")

# Verify: d0_fiber.T @ d0_fiber = L_fiber_vertex
L_fiber_v = np.diag(np.sum(adj * (np.array([[vertex_to_fiber[i] == vertex_to_fiber[j]
                    for j in range(Nv)] for i in range(Nv)]).astype(float)), axis=1)) \
            - adj * np.array([[float(vertex_to_fiber[i] == vertex_to_fiber[j])
                    for j in range(Nv)] for i in range(Nv)])

err_Lf = norm(d0_fiber.T @ d0_fiber - L_fiber_v)
print(f"  ||d0_fiber^T @ d0_fiber - L_fiber_vertex|| = {err_Lf:.2e}")


# =====================================================================
# STEP 7: Coexact spectrum analysis
# =====================================================================
print("\n" + "=" * 70)
print("STEP 7: FULL COEXACT SPECTRUM")
print("=" * 70)

evals_C_full, evecs_C = eigh(C)
evals_C_sorted = np.sort(evals_C_full)
nonzero_C = evals_C_sorted[np.abs(evals_C_sorted) > TOL]

print(f"  Coexact spectrum: {len(nonzero_C)} nonzero eigenvalues")

# Distinct eigenvalues with multiplicities
from collections import Counter
C_distinct = Counter()
for e in evals_C_sorted:
    C_distinct[round(e, 4)] += 1

print(f"  Distinct eigenvalues: {len(C_distinct)}")
print(f"\n  {'lambda_C':>10s} {'mult':>5s}")
print(f"  {'-'*20}")
for val, mult in sorted(C_distinct.items()):
    if mult > 0:
        print(f"  {val:10.4f} {mult:5d}")

# For each coexact eigenvector, compute the fiber/cross content
print(f"\n  Fiber content of coexact eigenvectors:")
print(f"  {'lambda_C':>10s} {'fiber_frac':>12s} {'cross_frac':>12s} {'mult':>5s}")
print(f"  {'-'*45}")

for val, mult in sorted(C_distinct.items()):
    if abs(val) < TOL:
        continue
    # Find eigenvectors at this eigenvalue
    indices = [i for i in range(Ne) if abs(evals_C_full[i] - val) < 0.01]
    if indices:
        fiber_fracs = []
        for idx in indices:
            v = evecs_C[:, idx]
            fiber_content = np.sum(v[is_fiber_edge]**2)
            fiber_fracs.append(fiber_content)
        avg_fiber = np.mean(fiber_fracs)
        print(f"  {val:10.4f} {avg_fiber:12.4f} {1-avg_fiber:12.4f} {mult:5d}")


# =====================================================================
# STEP 8: KEY QUESTION -- MASS GAP
# =====================================================================
print("\n" + "=" * 70)
print("STEP 8: COEXACT MASS GAP AND EW COMPARISON")
print("=" * 70)

# The coexact spectral gap (smallest nonzero C eigenvalue)
gap_C = nonzero_C[0]
print(f"  Coexact spectral gap = {gap_C:.10f}")

# Known: coexact gap = 7 - 4*phi (from paper)
gap_expected = 7 - 4*PHI
print(f"  Expected 7 - 4*phi = {gap_expected:.10f}")
print(f"  Match: {abs(gap_C - gap_expected) < TOL}")

# The vertex Box gap is a1*(2-phi) = 1.9098 (the cross gap)
# The coexact gap is 7-4*phi = 0.5279
# Ratio
print(f"\n  Coexact gap / vertex gap = {gap_C / (2-PHI):.6f}")
print(f"  = (7-4*phi)/(2-phi)")

# EW comparison
# m_Z/m_W = sqrt(13/10) = 1.14018
# Can we build a "coexact Box" and find mass ratios?
# For now, just list the coexact eigenvalue ratios
print(f"\n  Coexact eigenvalue ratios (relative to gap):")
distinct_nonzero = sorted(set(round(e, 4) for e in nonzero_C))
for i, val in enumerate(distinct_nonzero[:8]):
    ratio = val / gap_C
    print(f"    C_{i} = {val:10.4f}, ratio = {ratio:10.4f}")


print("\n" + "=" * 70)
print("EXP-564 COMPLETE")
print("=" * 70)

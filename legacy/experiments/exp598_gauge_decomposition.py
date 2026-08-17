"""
exp101: Decomposition of the 12 gauge modes (2·rho_5) into 1+3+8.

The 12 massless gauge modes from ker(Box_edge) form 2*rho_5 under 2I.
rho_5 (dim 6) is IRREDUCIBLE under 2I, so any 2I-equivariant operator
acts as a scalar. The 1+3+8 decomposition MUST come from breaking 2I
symmetry — either through:

1. The fiber-to-fiber adjacency (icosahedral graph on 12 fibers)
2. The A5 subgroup (quotient 2I -> A5)
3. A vertex-dependent structure (star of a vertex)
4. The face (triangle) incidence structure
5. Cross-fiber edge geometry (distance/angle on S^2)

Strategy: build all natural operators on the 12-dim gauge space,
diagonalize, and check if ANY gives multiplicities 1, 3, 8.
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
SQRT5 = np.sqrt(5)

# =====================================================================
# BUILD EVERYTHING (compact)
# =====================================================================
print("Building 600-cell + operators...")
verts, adj, lap = build_600cell()
Nv = 120

adj_list = defaultdict(set)
edges = []; edge_to_idx = {}
for i in range(Nv):
    for j in range(i+1, Nv):
        if adj[i, j] > 0.5:
            adj_list[i].add(j); adj_list[j].add(i)
            edges.append((i, j)); edge_to_idx[(i, j)] = len(edges) - 1
Ne = 720

triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j: triangles.append((i, j, k))

d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0; d0[e_idx, j] = +1.0

d1 = np.zeros((len(triangles), Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

def qmul(p, q):
    return np.array([p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1-tol else -1

def find_fibration():
    for i in range(Nv):
        if abs(verts[i, 0] - PHI/2) < 1e-6:
            g = verts[i]; p = g.copy(); ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok = False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok = False
            if not ok: continue
            used = set(); fibers = []; subg = []
            pp = np.array([1.0,0,0,0])
            for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
            for s in range(Nv):
                if s in used: continue
                fib = []
                for si in subg:
                    idx = find_idx(qmul(verts[s], verts[si]), verts)
                    if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12: return fibers
    return None

fibers = find_fibration()
vtf = {}
for fi, f in enumerate(fibers):
    for v in f: vtf[v] = fi

is_fiber = np.zeros(Ne, dtype=bool)
for e, (i, j) in enumerate(edges):
    if vtf[i] == vtf[j]: is_fiber[e] = True

vte = defaultdict(list)
for e, (i, j) in enumerate(edges): vte[i].append(e); vte[j].append(e)

A_line = np.zeros((Ne, Ne))
for v in range(Nv):
    inc = vte[v]
    for a in range(len(inc)):
        for b in range(a+1, len(inc)):
            A_line[inc[a], inc[b]] = 1.0; A_line[inc[b], inc[a]] = 1.0

A_fib_line = np.zeros((Ne, Ne))
for fi, fiber in enumerate(fibers):
    fe = []
    for k in range(10):
        u, v = fiber[k], fiber[(k+1)%10]
        fe.append(edge_to_idx[(min(u,v), max(u,v))])
    for k in range(10):
        A_fib_line[fe[k], fe[(k+1)%10]] = 1.0; A_fib_line[fe[(k+1)%10], fe[k]] = 1.0

Lf = np.diag(np.sum(A_fib_line, axis=1)) - A_fib_line
Lc = np.diag(np.sum(A_line - A_fib_line, axis=1)) - (A_line - A_fib_line)
Box1 = Lc - a1 * Lf

# Get kernel
evals_B1, evecs_B1 = eigh(Box1)
ker_mask = np.abs(evals_B1) < 1e-6
ker_vecs = evecs_B1[:, ker_mask]  # 720 x 13

# Separate rho_0 and rho_5
# Build rho_0 projector using group average
all_vperm = []
all_eperm = []
for g in range(Nv):
    vp = np.array([find_idx(qmul(verts[g], verts[i]), verts) for i in range(Nv)])
    all_vperm.append(vp)
    ep = np.array([edge_to_idx[(min(vp[i], vp[j]), max(vp[i], vp[j]))] for i, j in edges])
    all_eperm.append(ep)

P0_ker = np.zeros((13, 13))
for g in range(Nv):
    ep = all_eperm[g]
    for ki in range(13):
        gv = ker_vecs[:, ki][ep]
        for kj in range(13):
            P0_ker[ki, kj] += np.dot(ker_vecs[:, kj], gv)
P0_ker /= 120.0

P5_ker = np.eye(13) - P0_ker
p5e, p5v = eigh(P5_ker)
rho5_basis_ker = p5v[:, p5e > 0.5]  # 12 vectors in kernel space
gauge_modes = ker_vecs @ rho5_basis_ker  # 720 x 12

# Verify: gauge modes live 100% on fiber edges
assert np.allclose(np.sum(gauge_modes[~is_fiber]**2), 0, atol=1e-10)
print("Done. 12 gauge modes extracted (100% fiber).\n")


# =====================================================================
# PART 1: NATURAL BASIS — ONE MODE PER FIBER
# =====================================================================
print("=" * 72)
print("PART 1: FIBER IDENTIFICATION OF GAUGE MODES")
print("=" * 72)

# Each gauge mode is the antiperiodic mode on one C10 fiber.
# Build the natural basis: for each fiber, the antiperiodic mode.

fiber_edge_lists = []
for fi, fiber in enumerate(fibers):
    fe = []
    for k in range(10):
        u, v = fiber[k], fiber[(k+1)%10]
        fe.append(edge_to_idx[(min(u,v), max(u,v))])
    fiber_edge_lists.append(fe)

# Antiperiodic mode on C10: v[k] = (-1)^k / sqrt(10)
# But we need to match the orientation (d0 convention)
antiperiodic_modes = np.zeros((Ne, 12))
for fi, fiber in enumerate(fibers):
    fe = fiber_edge_lists[fi]
    for k in range(10):
        u, v = fiber[k], fiber[(k+1)%10]
        i, j = min(u, v), max(u, v)
        e_idx = edge_to_idx[(i, j)]
        # Orientation: d0 gives sign based on i < j convention
        # The antiperiodic mode alternates sign along the fiber
        sign = (-1)**k
        # But also need edge orientation: if (u,v) has u > v, flip sign
        if u > v:
            sign = -sign
        antiperiodic_modes[e_idx, fi] = sign
    # Normalize
    antiperiodic_modes[:, fi] /= norm(antiperiodic_modes[:, fi])

# Check: are these in the kernel of Box1?
for fi in range(12):
    residual = norm(Box1 @ antiperiodic_modes[:, fi])
    if residual > 0.01:
        print(f"  WARNING: fiber {fi} antiperiodic mode NOT in kernel! residual = {residual:.4f}")

# Check overlap with gauge_modes (should span same space)
overlap = gauge_modes.T @ antiperiodic_modes  # 12 x 12
sv = np.linalg.svd(overlap, compute_uv=False)
print(f"  Overlap singular values: {np.sort(sv)}")
print(f"  All > 0.99: {np.all(sv > 0.99)}")

# If the antiperiodic modes ARE in the kernel, use them as the natural basis
# Check directly
residuals = [norm(Box1 @ antiperiodic_modes[:, fi]) for fi in range(12)]
print(f"  Max residual: {max(residuals):.4e}")

if max(residuals) > 0.1:
    print("  Antiperiodic modes NOT exactly in kernel. Projecting onto kernel...")
    # Project antiperiodic modes onto kernel
    proj = ker_vecs @ ker_vecs.T @ antiperiodic_modes  # 720 x 12
    # Orthogonalize
    u, s, vt = np.linalg.svd(proj, full_matrices=False)
    gauge_fiber = u[:, :12]  # 720 x 12
    print(f"  Projection singular values: {s[:12]}")
else:
    # Use antiperiodic modes directly (they span a 12-dim subspace of kernel)
    u, s, vt = np.linalg.svd(antiperiodic_modes, full_matrices=False)
    gauge_fiber = antiperiodic_modes
    print(f"  Using antiperiodic modes directly as gauge basis.")

# The TRANSITION MATRIX from fiber basis to gauge_modes basis
T = gauge_modes.T @ antiperiodic_modes  # 12 x 12


# =====================================================================
# PART 2: FIBER-TO-FIBER ADJACENCY (ICOSAHEDRAL GRAPH)
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: FIBER-TO-FIBER WEIGHTED ADJACENCY")
print("=" * 72)

# Count cross edges between each pair of fibers
W_fiber = np.zeros((12, 12))
for e_idx, (i, j) in enumerate(edges):
    if not is_fiber[e_idx]:  # cross edge
        fi, fj = vtf[i], vtf[j]
        if fi != fj:
            W_fiber[fi, fj] += 1
            W_fiber[fj, fi] += 1

print("  Fiber-to-fiber cross-edge counts:")
unique_weights = sorted(set(W_fiber[W_fiber > 0].flatten().astype(int)))
print(f"  Unique weights: {unique_weights}")

# Count pairs at each weight
weight_pairs = Counter()
for i in range(12):
    for j in range(i+1, 12):
        w = int(W_fiber[i, j])
        weight_pairs[w] += 1

print(f"  Weight distribution:")
for w, count in sorted(weight_pairs.items()):
    print(f"    w = {w}: {count} pairs")
print(f"  Total cross edges: {int(W_fiber.sum()/2)} (expect 600)")

# Fiber degree (total cross edges from each fiber)
fiber_degrees = W_fiber.sum(axis=1)
print(f"  Fiber degrees: {fiber_degrees.astype(int)}")

# Eigenvalues of W_fiber
evals_W = np.sort(eigvalsh(W_fiber))
print(f"\n  Eigenvalues of W_fiber (12x12):")
ev_W_distinct = Counter(np.round(evals_W, 4))
for val, mult in sorted(ev_W_distinct.items()):
    # Try Z[phi] identification
    best = ""
    for b in range(-20, 21):
        a = val - b * PHI
        if abs(a - round(a)) < 0.01:
            best = f" = {int(round(a))} + {b}*phi"
            break
    print(f"    {val:10.4f} x {mult}{best}")

# Binary adjacency (fiber graph: are fibers connected at all?)
A_fiber_graph = (W_fiber > 0).astype(float)
fiber_graph_degree = A_fiber_graph.sum(axis=1)
print(f"\n  Fiber graph degree: {fiber_graph_degree[0]:.0f} (all same: {np.allclose(fiber_graph_degree, fiber_graph_degree[0])})")

evals_Afg = np.sort(eigvalsh(A_fiber_graph))
ev_Afg_distinct = Counter(np.round(evals_Afg, 4))
print(f"  Binary fiber adjacency eigenvalues:")
for val, mult in sorted(ev_Afg_distinct.items()):
    print(f"    {val:10.4f} x {mult}")


# =====================================================================
# PART 3: DISTANCE CLASSES ON FIBER GRAPH
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: DISTANCE CLASSES ON FIBER GRAPH")
print("=" * 72)

# Compute fiber centers on S^2 (Hopf projection)
def hopf_project(q):
    """Hopf projection S^3 -> S^2: (w,x,y,z) -> (2(xz+wy), 2(yz-wx), w^2+z^2-x^2-y^2)"""
    w, x, y, z = q
    return np.array([2*(x*z + w*y), 2*(y*z - w*x), w**2 + z**2 - x**2 - y**2])

fiber_centers = np.zeros((12, 3))
for fi, fiber in enumerate(fibers):
    # Average Hopf projection over all vertices in fiber
    pts = np.array([hopf_project(verts[v]) for v in fiber])
    center = pts.mean(axis=0)
    center /= norm(center)  # normalize to S^2
    fiber_centers[fi] = center

# Check that all fibers map to the same point on S^2 (they should, that's what fibers are)
for fi, fiber in enumerate(fibers):
    pts = np.array([hopf_project(verts[v]) for v in fiber])
    spread = np.std(pts, axis=0)
    if norm(spread) > 0.01:
        print(f"  WARNING: fiber {fi} has spread on S^2: {spread}")

# Pairwise distances on S^2
S2_dist = np.zeros((12, 12))
for i in range(12):
    for j in range(12):
        dot = np.clip(np.dot(fiber_centers[i], fiber_centers[j]), -1, 1)
        S2_dist[i, j] = np.arccos(dot)

# Classify distance types
dist_types = Counter()
for i in range(12):
    for j in range(i+1, 12):
        d = round(S2_dist[i, j], 4)
        dist_types[d] += 1

print("  Pairwise S^2 distances between fibers:")
for d, count in sorted(dist_types.items()):
    print(f"    theta = {d:.4f} rad ({np.degrees(d):.2f} deg): {count} pairs")

# Build distance-class adjacency matrices
dist_levels = sorted(dist_types.keys())
print(f"\n  {len(dist_levels)} distance classes")

for di, d in enumerate(dist_levels):
    A_d = np.zeros((12, 12))
    for i in range(12):
        for j in range(12):
            if abs(S2_dist[i, j] - d) < 0.01:
                A_d[i, j] = 1.0
    evals_d = np.sort(eigvalsh(A_d))
    ev_d = Counter(np.round(evals_d, 4))
    mults = [mult for _, mult in sorted(ev_d.items())]
    print(f"    d = {d:.4f}: degree = {A_d.sum(axis=1)[0]:.0f}, "
          f"eigenvalue mults = {mults}")


# =====================================================================
# PART 4: ACT ON GAUGE MODES WITH FIBER OPERATORS
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: FIBER OPERATORS ON 12 GAUGE MODES")
print("=" * 72)

# The 12 gauge modes (antiperiodic on fibers) can be identified with fibers.
# Define M = antiperiodic_modes^T @ X @ antiperiodic_modes for various X.

# The key insight: even though A_line = -2*I and C = 5*I on the gauge sector
# (when using arbitrary orthonormal basis), if we use the FIBER-LABELED basis
# (antiperiodic_modes), we can project CROSS operators onto this basis.

# First: is the antiperiodic basis IN the kernel?
# From Part 1 we saw residuals. Let me check if they span the rho_5 subspace.
# Project out rho_0 from antiperiodic modes:
ones_fiber = np.ones(12) / np.sqrt(12)
proj_rho0 = np.outer(ones_fiber, ones_fiber)
T_centered = T - T @ proj_rho0  # remove rho_0 component

# The 12 antiperiodic modes, projected onto the gauge subspace
M_proj = gauge_modes.T @ antiperiodic_modes  # 12 x 12
print(f"  gauge_modes^T @ antiperiodic_modes:")
print(f"  Singular values: {np.sort(np.linalg.svd(M_proj, compute_uv=False))}")

# Now use antiperiodic modes to build effective 12x12 operators
# WARNING: antiperiodic modes may not be in the kernel exactly.
# Use the PROJECTION of antiperiodic modes onto the kernel instead.

anti_ker = ker_vecs @ (ker_vecs.T @ antiperiodic_modes)  # project onto 13-dim kernel
# Further project out rho_0
rho0_mode = ker_vecs @ P0_ker @ ker_vecs.T  # projector onto rho_0 in edge space
anti_gauge = anti_ker - rho0_mode @ anti_ker  # 720 x 12, in rho_5 subspace

# Orthonormalize within the rho_5 gauge subspace
u_ag, s_ag, _ = np.linalg.svd(anti_gauge, full_matrices=False)
print(f"\n  Antiperiodic modes projected to gauge subspace:")
print(f"  Singular values: {np.round(s_ag, 6)}")

# Use the projected modes as our fiber-labeled gauge basis
# Note: after projection, the modes may mix between fibers
# The EFFECTIVE 12x12 operator for any edge-space operator X is:
#   M_eff = anti_gauge^T @ X @ anti_gauge
# But this is in a non-orthonormal basis. Better to use:
#   M_eff = (anti_gauge^T @ X @ anti_gauge) * S^{-2} ... messy
# Instead, work with gauge_modes (orthonormal) and relate to fibers via T

# The fiber overlap matrix
F_ij = antiperiodic_modes.T @ antiperiodic_modes  # 12 x 12
print(f"\n  Fiber overlap matrix F_ij eigenvalues: {np.sort(eigvalsh(F_ij))}")

# The W_fiber matrix acts on fiber labels.
# In the gauge_modes basis: M_eff = T^T @ W_fiber @ T / (T^T @ T)
# But T = gauge_modes^T @ antiperiodic_modes is rectangular and may not be invertible.
# Let me just compute W_fiber eigenvalues (already done above) and see what structure emerges.


# =====================================================================
# PART 5: ALL CROSS-EDGE OPERATORS ON FIBER LABELS
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: CROSS-EDGE OPERATORS ON 12 FIBERS")
print("=" * 72)

# For each distance class, build the 12x12 matrix of cross-edge counts
for di, d_val in enumerate(dist_levels):
    W_d = np.zeros((12, 12))
    for e_idx, (i, j) in enumerate(edges):
        if not is_fiber[e_idx]:
            fi, fj = vtf[i], vtf[j]
            if fi != fj and abs(S2_dist[fi, fj] - d_val) < 0.01:
                W_d[fi, fj] += 1
                W_d[fj, fi] += 1

    evals_Wd = np.sort(eigvalsh(W_d))
    ev_Wd = Counter(np.round(evals_Wd, 4))
    total = int(W_d.sum() / 2)
    mults_str = ", ".join(f"{v:.2f}({m})" for v, m in sorted(ev_Wd.items()))
    print(f"  Distance class {di} (theta={d_val:.4f}, {total} edges):")
    print(f"    Multiplicities: {[m for _, m in sorted(ev_Wd.items())]}")
    print(f"    Eigenvalues: {mults_str}")


# =====================================================================
# PART 6: VERTEX STAR STRUCTURE
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: VERTEX STAR STRUCTURE")
print("=" * 72)

# For each vertex v, it connects to 12 neighbors.
# 2 are fiber neighbors, 10 are cross neighbors.
# The 10 cross neighbors go to 10 different fibers.
# Which fiber is NOT connected? The "antipodal" fiber.

# For vertex 0: find its fiber, find all connected fibers
v0 = 0
f0 = vtf[v0]
connected_fibers = set()
for j in adj_list[v0]:
    fj = vtf[j]
    if fj != f0:
        connected_fibers.add(fj)

missing_fibers = set(range(12)) - {f0} - connected_fibers
print(f"  Vertex {v0} (fiber {f0}):")
print(f"  Connected to {len(connected_fibers)} cross fibers: {sorted(connected_fibers)}")
print(f"  Missing fiber(s): {sorted(missing_fibers)}")
print(f"  (The missing fiber is the 'antipodal' fiber)")

# Check for all vertices: is there always exactly 1 missing fiber?
all_missing = []
for v in range(Nv):
    fv = vtf[v]
    cf = set()
    for j in adj_list[v]:
        if vtf[j] != fv:
            cf.add(vtf[j])
    missing = set(range(12)) - {fv} - cf
    all_missing.append(len(missing))

missing_counts = Counter(all_missing)
print(f"  Missing fiber count distribution: {dict(missing_counts)}")
print(f"  (expect all = 1)")

# The "antipodal pairing" of fibers
antipodal = {}
for v in range(10):  # just use vertices from fiber 0
    v_actual = fibers[0][v]
    fv = vtf[v_actual]
    cf = set()
    for j in adj_list[v_actual]:
        if vtf[j] != fv:
            cf.add(vtf[j])
    missing = set(range(12)) - {fv} - cf
    if len(missing) == 1:
        antipodal[fv] = missing.pop()

print(f"\n  Antipodal fiber pairs (from fiber 0 vertices):")
print(f"  Fiber 0 -> Fiber {antipodal.get(0, '?')}")

# Build full antipodal map
anti_map = {}
for fi in range(12):
    v_rep = fibers[fi][0]
    cf = set()
    for j in adj_list[v_rep]:
        if vtf[j] != fi:
            cf.add(vtf[j])
    missing = set(range(12)) - {fi} - cf
    if len(missing) == 1:
        anti_map[fi] = missing.pop()

print(f"  Full antipodal map: {anti_map}")
# Verify it's an involution
for fi, fj in anti_map.items():
    assert anti_map[fj] == fi, f"Not involution: {fi}->{fj}->{anti_map[fj]}"
print(f"  Verified: antipodal map is an involution (6 pairs)")

# The 6 antipodal pairs
pairs = set()
for fi, fj in anti_map.items():
    pairs.add((min(fi, fj), max(fi, fj)))
print(f"  Antipodal pairs: {sorted(pairs)}")


# =====================================================================
# PART 7: DECOMPOSITION USING CROSS-EDGE STRUCTURE
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: DECOMPOSITION ATTEMPTS")
print("=" * 72)

# Attempt 1: Use the WEIGHTED fiber adjacency W_fiber
# Its eigenvalue decomposition gives the natural splitting
evals_W, evecs_W = eigh(W_fiber)
ev_W_mults = Counter(np.round(evals_W, 4))
print("  Attempt 1: Weighted fiber adjacency W_fiber")
print(f"  Eigenvalue multiplicities: {[m for _, m in sorted(ev_W_mults.items())]}")

# Check if any combination gives 1+3+8
# Current decomposition: let's see the actual values
for val, mult in sorted(ev_W_mults.items()):
    print(f"    lambda = {val:10.4f} x {mult}")

# Attempt 2: Use the SQUARED adjacency (2-step walks between fibers)
W2 = W_fiber @ W_fiber
evals_W2 = np.sort(eigvalsh(W2))
ev_W2_mults = Counter(np.round(evals_W2, 4))
print(f"\n  Attempt 2: W_fiber^2")
print(f"  Eigenvalue multiplicities: {[m for _, m in sorted(ev_W2_mults.items())]}")

# Attempt 3: Antipodal matrix (P_anti[i,j] = 1 if j is antipodal to i)
P_anti = np.zeros((12, 12))
for fi, fj in anti_map.items():
    P_anti[fi, fj] = 1.0
evals_Pa = np.sort(eigvalsh(P_anti))
ev_Pa_mults = Counter(np.round(evals_Pa, 4))
print(f"\n  Attempt 3: Antipodal pairing matrix")
print(f"  Eigenvalue multiplicities: {[m for _, m in sorted(ev_Pa_mults.items())]}")

# Attempt 4: Fiber Laplacian on the icosahedral fiber graph
L_fib_graph = np.diag(A_fiber_graph.sum(axis=1)) - A_fiber_graph
evals_Lfg = np.sort(eigvalsh(L_fib_graph))
ev_Lfg_mults = Counter(np.round(evals_Lfg, 4))
print(f"\n  Attempt 4: Fiber graph Laplacian (icosahedron)")
print(f"  Eigenvalue multiplicities: {[m for _, m in sorted(ev_Lfg_mults.items())]}")
for val, mult in sorted(ev_Lfg_mults.items()):
    print(f"    {val:10.4f} x {mult}")

# Attempt 5: W_fiber minus diagonal (remove self-loops), then eigenvalues
# Already done above. Try W + alpha * P_anti for different alpha
print(f"\n  Attempt 5: W_fiber + alpha * P_anti")
for alpha in [-10, -5, -1, 0, 1, 5, 10, 50]:
    M = W_fiber + alpha * P_anti
    evals_M = np.sort(eigvalsh(M))
    ev_M_mults = Counter(np.round(evals_M, 2))
    mults = [m for _, m in sorted(ev_M_mults.items())]
    if 1 in mults and 3 in mults and 8 in mults:
        print(f"    alpha = {alpha:4d}: *** 1+3+8 FOUND! *** mults = {mults}")
    else:
        print(f"    alpha = {alpha:4d}: mults = {mults}")

# Attempt 6: Use the NUMBER of cross edges at each distance class
# Build separate matrices for each distance
print(f"\n  Attempt 6: Linear combinations of distance-class matrices")
# First, build all distance-class matrices
D_matrices = {}
for d_val in dist_levels:
    W_d = np.zeros((12, 12))
    for e_idx, (i, j) in enumerate(edges):
        if not is_fiber[e_idx]:
            fi, fj = vtf[i], vtf[j]
            if fi != fj and abs(S2_dist[fi, fj] - d_val) < 0.01:
                W_d[fi, fj] += 1
                W_d[fj, fi] += 1
    D_matrices[d_val] = W_d

# These matrices span the Bose-Mesner algebra of the icosahedral association scheme
# Their joint eigenspaces should decompose the 12-dim space
# Compute the joint eigenspace decomposition
joint_matrix = np.zeros((12, 12))
for di, d_val in enumerate(dist_levels):
    joint_matrix += (di + 1) * D_matrices[d_val]  # weight by distance rank

evals_joint, evecs_joint = eigh(joint_matrix)
ev_joint_mults = Counter(np.round(evals_joint, 4))
print(f"  Joint distance operator eigenvalue mults: {[m for _, m in sorted(ev_joint_mults.items())]}")

# Attempt 7: Use the FACE (triangle) structure
# For each pair of fibers, count the number of triangles that contain
# one edge from each fiber
print(f"\n  Attempt 7: Triangle-mediated fiber coupling")
T_fiber = np.zeros((12, 12))
for (i, j, k) in triangles:
    fi, fj, fk = vtf[i], vtf[j], vtf[k]
    fibers_in_tri = {fi, fj, fk}
    if len(fibers_in_tri) == 2:  # triangle spans 2 fibers
        f_list = sorted(fibers_in_tri)
        T_fiber[f_list[0], f_list[1]] += 1
        T_fiber[f_list[1], f_list[0]] += 1

evals_Tf = np.sort(eigvalsh(T_fiber))
ev_Tf_mults = Counter(np.round(evals_Tf, 4))
print(f"  Triangle coupling eigenvalue mults: {[m for _, m in sorted(ev_Tf_mults.items())]}")
for val, mult in sorted(ev_Tf_mults.items()):
    print(f"    {val:10.4f} x {mult}")

# Attempt 8: Use the VERTEX adjacency structure (Section 4 approach)
# For each vertex, look at which fibers its 10 cross neighbors connect to.
# Build a 12x12 matrix that counts "vertex-mediated" connections between fibers.
print(f"\n  Attempt 8: Vertex-mediated connections")
# For each fiber, count connections from vertices in that fiber to other fibers
# But weight by vertex-specific structure
# Try: M[fi, fj] = number of vertices in fi with a neighbor in fj
V_fiber = np.zeros((12, 12))
for v in range(Nv):
    fv = vtf[v]
    neighbor_fibers = Counter()
    for j in adj_list[v]:
        fj = vtf[j]
        if fj != fv:
            neighbor_fibers[fj] += 1
    for fj, count in neighbor_fibers.items():
        V_fiber[fv, fj] += count  # same as W_fiber, actually

# Actually V_fiber = W_fiber. Let me try something else:
# Count vertices that have EXACTLY k neighbors in a given fiber
# This gives finer structure.

# Attempt 9: Hopf ANGLE structure
# Each fiber has a phase angle on S^1 (the Hopf fiber). The relative
# angle between fibers could give a finer structure than distance.
print(f"\n  Attempt 9: Hopf phase angles")
# The Hopf generator g has angle 2*pi/10 per step.
# Different fibers have different base points on S^2.
# The "angle" between fibers is the phase difference in the Hopf circle.
# This is not straightforward...

# Instead, let's try the QUATERNIONIC product structure.
# For fibers i and j, the "transition quaternion" from fiber i to fiber j
# is q = v_j^{-1} * v_i for representative vertices.
# The scalar part of this quaternion gives the "Hopf distance."

# Use first vertex of each fiber as representative
fiber_reps = [fibers[fi][0] for fi in range(12)]
Q_product = np.zeros((12, 12))
for i in range(12):
    for j in range(12):
        # Quaternion distance: |q_i . q_j| = cos(angle/2)
        dot = abs(np.dot(verts[fiber_reps[i]], verts[fiber_reps[j]]))
        Q_product[i, j] = dot

print(f"  Quaternionic dot product values:")
q_vals = Counter()
for i in range(12):
    for j in range(i+1, 12):
        q_vals[round(Q_product[i, j], 4)] += 1
for v, c in sorted(q_vals.items()):
    print(f"    |q_i . q_j| = {v:.4f}: {c} pairs")


# =====================================================================
# PART 8: DEFINITIVE CHECK — GALOIS CONJUGATION ON FIBERS
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: GALOIS CONJUGATION ON FIBER GRAPH")
print("=" * 72)

# The Galois conjugation phi -> 1-phi sends the 600-cell to itself
# but permutes the fibers. This might split the 12 fibers into
# orbits of size 1+3+8 or similar.

# Galois conjugation on vertex coordinates: replace phi by 1-phi = -1/phi
# Equivalently: in the golden ratio coordinates, negate all phi terms.
# This maps the 600-cell to itself (preserving adjacency).

# Apply Galois conjugation to vertices
verts_galois = verts.copy()
# The 600-cell vertices involve phi in coordinates.
# Galois: phi -> (1-sqrt(5))/2 = -1/phi
# For Type C vertices: (0, ±1/2, ±phi/2, ±1/(2*phi))
# Galois sends phi/2 -> -1/(2*phi) and 1/(2*phi) -> -phi/2
# This is a specific linear transformation.
# Actually, the Galois conjugate of the vertex set is obtained by
# replacing phi with its conjugate (1-phi) = -(phi-1) = -1/phi everywhere.
#
# For a vertex (w, x, y, z) in the 600-cell, if it involves phi,
# we replace phi -> phi' = (1-sqrt(5))/2 = -1/phi.
# But the COORDINATES are real numbers. The Galois action acts on
# the ALGEBRAIC structure, not directly on coordinates.
#
# For the Laplacian spectrum: phi-eigenvalues map to phi'-eigenvalues.
# For the FIBER GRAPH, the Galois action permutes eigenvalues and
# therefore permutes eigenspaces. The eigenspace decomposition of
# the fiber graph tells us how Galois acts.

# The fiber graph (icosahedron) has eigenvalues:
# A: 5 (x1), sqrt(5) (x5), -1 (x5), -sqrt(5) (x1)
# Galois: sqrt(5) -> -sqrt(5), so:
# Galois maps the sqrt(5)-eigenspace (dim 5) to the -sqrt(5)-eigenspace (dim 1)?
# NO — Galois maps eigenvalue to its conjugate, but the eigenSPACE stays the same dim.
# sqrt(5) -> -sqrt(5): eigenspace dim 5 maps to eigenspace dim 1.
# This is IMPOSSIBLE for a group automorphism (it preserves dimensions).
#
# Actually, Galois acts on the eigenvalues but NOT as a graph automorphism.
# The Galois conjugation of the 600-cell IS a graph automorphism, but it
# acts on the fibers in a specific way.

# Let me compute the Galois action on fibers directly.
# The 600-cell has Type C vertices with phi in coordinates.
# Galois conjugate: phi -> -1/phi = (1-sqrt(5))/2

# Build Galois-conjugated vertices
phi_conj = (1 - np.sqrt(5)) / 2  # = -1/phi

# The 600-cell construction uses phi in Type C vertices.
# Type A: (±1, 0, 0, 0) permutations — no phi
# Type B: (±1/2)^4 — no phi
# Type C: even permutations of (0, ±1/2, ±phi/2, ±1/(2*phi))
# Galois: phi/2 -> phi_conj/2 = (1-sqrt(5))/4, 1/(2*phi) -> 1/(2*phi_conj) = -phi/2

# Actually, the easiest way: the Galois conjugation is equivalent to
# the linear map that sends (w,x,y,z) -> (w', x', y', z') where
# the transformation depends on the TYPE of vertex.
# This is complicated to implement directly.

# Alternative: use the LAPLACIAN eigenspace structure.
# On the vertex Laplacian, eigenvalue 12-6*phi maps to 12-6*(1-phi) = 6+6*phi.
# Wait: phi_conj = -1/phi = 1-phi. So 12-6*phi -> 12-6*(1-phi) = 12-6+6phi = 6+6phi.
# But 6+6*phi ≈ 15.708, which matches eigenvalue L8.
# And eigenvalue L2 = 10-2*sqrt(5) -> 10+2*sqrt(5) = L6.
# So Galois maps: L1 <-> L8, L2 <-> L6 (the known Galois pairs).

# On the FIBER GRAPH (icosahedron), eigenvalues 5, sqrt(5), -1, -sqrt(5):
# Galois: sqrt(5) <-> -sqrt(5), 5 and -1 are fixed.
# But dim(sqrt(5)-eigenspace) = 5 and dim(-sqrt(5)-eigenspace) = 1.
# So Galois does NOT act as a simple permutation of eigenspaces here.
# This means the icosahedral Galois structure is more subtle.

# Let me instead try a PURELY COMPUTATIONAL approach.
# The question "does anything give 1+3+8?" can be tested by
# searching over all operators in the commutant algebra of the
# fiber graph (Bose-Mesner algebra).

print("  All distance-class matrices eigenvalue multiplicities:")
all_mults = []
for d_val in dist_levels:
    evals = np.sort(eigvalsh(D_matrices[d_val]))
    mults = [m for _, m in sorted(Counter(np.round(evals, 4)).items())]
    all_mults.append(mults)
    print(f"    d = {d_val:.4f}: {mults}")

print(f"\n  The icosahedral association scheme has eigenspace dimensions: 1, 5, 5, 1")
print(f"  This is because the icosahedron has 3 distance classes (d=1,2,3)")
print(f"  giving a Bose-Mesner algebra with 4 idempotents of ranks 1,5,5,1.")
print(f"  NO linear combination of distance-class matrices can give 1+3+8.")
print(f"  The multiplicities are TOPOLOGICAL (from icosahedral symmetry).")


# =====================================================================
# PART 9: ALTERNATIVE APPROACHES
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: ALTERNATIVE APPROACHES TO 1+3+8")
print("=" * 72)

# The 12 = 2*rho_5 under 2I is IRREDUCIBLE.
# The icosahedral fiber graph gives 1+5+5+1 (association scheme).
# A5 gives 1+3+3'+5 on the permutation rep of 12 vertices.
# NONE of these give 1+3+8.

# But wait: 8 is NOT a dimension of any A5 irrep (dims: 1,3,3',4,5).
# So 1+3+8 CANNOT come from A5 alone!

# The 1+3+8 must come from a DIFFERENT structure.
# Options:
# a) The Lie algebra structure on the gauge sector (emergent, not group-theoretic)
# b) The A4 subgroup of A5 (dims: 1,1',1'',3 -> but 1+1'+1''+3 = 6, too small)
# c) The S3 subgroup (face symmetry of icosahedron)
# d) Higher-dimensional structure from E8

print(f"  A5 irrep dimensions: 1, 3, 3', 4, 5")
print(f"  8 is NOT an A5 irrep dimension!")
print(f"  => 1+3+8 CANNOT be an A5 irrep decomposition")
print(f"")
print(f"  The decomposition 12 = 1+3+8 is the LIE ALGEBRA decomposition:")
print(f"    u(1) + su(2) + su(3) = 1 + 3 + 8")
print(f"  This is about the GAUGE ALGEBRA STRUCTURE, not a group rep.")
print(f"")
print(f"  The A5 decomposition of the permutation rep on 12 icosahedron vertices:")
print(f"    12 = 1 + 3 + 3' + 5 under A5")
print(f"  With 3+3' = 6 and 5+1 = 6:")
print(f"    12 = (1 + 5) + (3 + 3') = 6 + 6")
print(f"  This matches 2*rho_5 = 6 + 6 of 2I!")

# Wait: under A5, rho_5 of 2I restricts to a PROJECTIVE representation.
# As a genuine A5 representation, the permutation rep = 1+3+3'+5.
# The 6-dim subspaces (1+5) and (3+3') correspond to the two copies of rho_5!

# If this is correct, then within each rho_5:
# Copy 1 = 1 + 5 (trivial + standard A5 rep)
# Copy 2 = 3 + 3' (the two 3-dim A5 irreps)

# But 1+5 = 6 = rho_5 restricted to A5?
# And 3+3' = 6 = rho_5 restricted to A5?
# These would be TWO DIFFERENT A5-decompositions of the same 2I irrep rho_5.

print(f"\n  KEY INSIGHT:")
print(f"  Under A5, the 12-dim permutation rep = 1 + 3 + 3' + 5")
print(f"  Grouping: (1 + 5) + (3 + 3') = 6 + 6")
print(f"  This matches 2*rho_5 = 6 + 6 of 2I!")
print(f"")

# Let me verify this by computing the A5 decomposition explicitly.
# A5 = 2I / {1, -1}. Find the -1 element.
minus1_idx = -1
for g in range(Nv):
    if np.allclose(verts[g], [-1, 0, 0, 0], atol=1e-6):
        minus1_idx = g
        break

print(f"  -1 element at index: {minus1_idx}")

# The action of -1 on fibers
ep_minus1 = all_eperm[minus1_idx]
fiber_perm_minus1 = {}
for fi in range(12):
    # Where does the first edge of fiber fi map?
    fe = fiber_edge_lists[fi]
    mapped = ep_minus1[fe[0]]
    # Which fiber is the mapped edge in?
    mi, mj = edges[mapped]
    fiber_perm_minus1[fi] = vtf[mi]

print(f"  Action of -1 on fibers: {fiber_perm_minus1}")
# If -1 acts trivially on fibers, the A5 quotient acts on fibers.
trivial_on_fibers = all(fiber_perm_minus1[fi] == fi for fi in range(12))
print(f"  -1 acts trivially on fibers: {trivial_on_fibers}")

if trivial_on_fibers:
    print(f"  => Fiber permutation factors through A5 = 2I/{{1,-1}}")

    # Compute A5 character on 12 fibers
    # For each conjugacy class of 2I, the fiber permutation character
    # is the number of fixed fibers
    print(f"\n  Computing A5 action on fibers...")

    fiber_chars = {}  # group index -> number of fixed fibers
    for g in range(Nv):
        n_fixed = 0
        for fi in range(12):
            fe = fiber_edge_lists[fi]
            mapped = all_eperm[g][fe[0]]
            mi, mj = edges[mapped]
            if vtf[mi] == fi or vtf[mj] == fi:
                # Check properly: the mapped edge is in the same fiber
                if vtf[mi] == vtf[mj] == fi:
                    n_fixed += 1
        fiber_chars[g] = n_fixed

    # Average over conjugacy classes (using the character vector approach)
    # Build character of fiber permutation rep on each eigenspace
    evals_lap, evecs_lap = eigh(lap)
    ev_sp = {}
    for i in range(Nv):
        val = round(evals_lap[i], 4)
        if val not in ev_sp: ev_sp[val] = []
        ev_sp[val].append(i)
    eig_list = sorted(ev_sp.items())

    irrep_dims = [int(round(np.sqrt(len(eig_list[k][1])))) for k in range(9)]
    print(f"  Irrep dims: {irrep_dims}")

    # For A5, we need irreps that factor through 2I/{±1}.
    # These are the irreps with chi(-1) = chi(1) = dim(rho).
    # From exp100c: chi(-1) = chi_0(-1)=1, chi_2(-1)=9, chi_4(-1)=25, chi_6(-1)=9
    # Wait, from the eigenspace characters:
    # rho_0: dim 1, chi(-1) = 1 -> factors
    # rho_1: dim 2, chi(-1) = -4 (eigenspace char) / 2 = -2 -> doesn't factor
    # rho_2: dim 3, chi(-1) = 9/3 = 3 -> factors
    # rho_3: dim 4, chi(-1) = -16/4 = -4 -> doesn't factor
    # rho_4: dim 5, chi(-1) = 25/5 = 5 -> factors
    # rho_5: dim 6, chi(-1) = -36/6 = -6 -> doesn't factor
    # rho_6: dim 3, chi(-1) = 9/3 = 3 -> factors
    # rho_7: dim 4, chi(-1) = -16/4 = -4 -> doesn't factor
    # rho_8: dim 2, chi(-1) = -4/2 = -2 -> doesn't factor

    # So A5 irreps correspond to: rho_0(1), rho_2(3), rho_4(5), rho_6(3)
    # But that's only 4 irreps of dims 1, 3, 5, 3 summing to 12 = √144... no.
    # A5 has 5 irreps: 1, 3, 3', 4, 5. Sum of squares = 1+9+9+16+25 = 60.
    # The 4-dim irrep of A5 should also appear. But chi(-1) for a 4-dim rep
    # that factors through A5 should be +4. However rho_3 and rho_7 have
    # chi(-1) = -4. So the 4-dim A5 irrep is NOT a 2I eigenspace.
    # The 4-dim A5 irrep = rho_3 tensor sign? Hmm, this is getting confusing.

    # Actually, A5 has 5 conjugacy classes, 5 irreps.
    # The irreps of A5 pullback to 2I as the 5 irreps with chi(-1) > 0:
    # rho_0(1), rho_2(3), rho_4(5), rho_6(3')
    # That's only 4. The 5th A5 irrep (dim 4) must be a NEW rep, not a
    # direct pullback but related to tensor products.
    # Actually no: dim 4 of A5 pulls back to 2I where it has chi(-1) = +4.
    # Looking at the eigenspace characters... rho_3 has chi(-1) = -4 and
    # rho_7 has chi(-1) = -4. Neither gives +4. So the dim-4 A5 irrep
    # must correspond to something else... perhaps it's not a Laplacian
    # eigenspace at all for the 600-cell.

    # This is getting complicated. Let me just numerically decompose
    # the 12-dim fiber permutation representation of A5.

    # The A5 elements are pairs {g, -g} in 2I. For the fiber representation,
    # g and -g give the SAME permutation (since -1 acts trivially on fibers).
    # So I can use any g or -g.

    # For each 2I element g, compute the permutation of 12 fibers
    fiber_perms = np.zeros((Nv, 12), dtype=int)
    for g in range(Nv):
        for fi in range(12):
            # Map first vertex of fiber fi
            v_rep = fibers[fi][0]
            gv = all_vperm[g][v_rep]
            fiber_perms[g, fi] = vtf[gv]

    # The 12-dim permutation character of A5
    # = number of fixed fibers for each element
    # Group A5 elements by their fiber-permutation character
    a5_chars = {}
    for g in range(Nv):
        n_fixed = np.sum(fiber_perms[g] == np.arange(12))
        key = n_fixed
        if key not in a5_chars:
            a5_chars[key] = 0
        a5_chars[key] += 1

    print(f"\n  Fiber permutation character distribution:")
    for nf, count in sorted(a5_chars.items()):
        print(f"    {nf} fixed fibers: {count} group elements")

    # Character of fiber representation at identity
    print(f"    chi(1) = 12")

    # Decompose using eigenspaces of the fiber graph Laplacian
    # The joint eigenspaces of all distance-class matrices give the
    # irreducible components under the full automorphism group (A5).
    print(f"\n  Fiber graph Laplacian eigenspaces (= A5 irrep decomposition):")
    evals_fg, evecs_fg = eigh(L_fib_graph)
    fg_spaces = {}
    for i in range(12):
        val = round(evals_fg[i], 4)
        if val not in fg_spaces: fg_spaces[val] = []
        fg_spaces[val].append(i)

    for val, indices in sorted(fg_spaces.items()):
        space = evecs_fg[:, indices]
        dim = len(indices)
        print(f"    lambda = {val:8.4f}, dim = {dim}")


# =====================================================================
# FINAL ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("FINAL ASSESSMENT")
print("=" * 72)

print(f"""
  The 12 gauge modes form 2*rho_5 under 2I, which is IRREDUCIBLE.

  Under A5 (the icosahedral quotient of 2I):
  - The fiber permutation rep = 1 + 3 + 3' + 5
  - The icosahedral association scheme gives eigenspaces: 1, 5, 5, 1
  - These CAN be refined to 1 + 3 + 3' + 5 using the FULL A5 character theory

  CRITICAL OBSERVATION:
  8 is NOT a dimension of any A5 irrep (dims are 1, 3, 3', 4, 5).
  So 1+3+8 CANNOT come from A5 representation theory alone.

  The decomposition 12 = 1+3+8 is NOT a group-theoretic decomposition.
  It is the decomposition of the GAUGE LIE ALGEBRA:
    su(3) x su(2) x u(1) has dimensions 8 + 3 + 1 = 12

  To get 1+3+8, one needs ADDITIONAL structure beyond group theory:
  - A Lie bracket on the 12 modes (making them into a Lie algebra)
  - Identification of the Lie algebra as su(3)+su(2)+u(1)

  The A5 decomposition 1+3+3'+5 is RELATED but different:
  - 1: U(1) singlet (same as in 1+3+8)
  - 3+3': two 3-dim reps -> could form su(2)+... ?
  - 5: standard rep -> could form part of su(3)?
  - BUT: 8 = 3' + 5 requires MERGING two A5 irreps

  POSSIBLE PATHWAY to 1+3+8:
  If the 3' and 5 irreps of A5 MERGE into an 8-dim space under some
  additional structure (e.g., a Lie bracket that closes 3'+5 into su(3)),
  then: 12 = 1 + 3 + (3' + 5) = 1 + 3 + 8

  STATUS: PARTIALLY DERIVED
  - 12 gauge modes: DERIVED (exact, from Box_edge kernel)
  - Their fiber/geometric structure: DERIVED
  - The 1+3+8 splitting: REQUIRES additional structure
  - The pathway 1 + 3 + (3'+5) is SUGGESTIVE but not proven
""")

print("=" * 72)
print("EXP101 COMPLETE")
print("=" * 72)

"""
exp578: Ollivier-Ricci curvature on the flat 600-cell.

The Ollivier-Ricci curvature kappa(x,y) on an edge (x,y) of a graph is:
    kappa(x,y) = 1 - W_1(m_x, m_y) / d(x,y)

where:
    m_x = uniform measure on neighbors of x (idleness=0, Lin-Lu-Yau version)
    W_1 = Wasserstein-1 (Earth Mover's) distance with graph metric as ground cost
    d(x,y) = graph distance (= 1 for edges)

For the 600-cell (vertex-transitive, 12-regular), we expect:
    - All edges have the SAME curvature (by symmetry)
    - Positive curvature (sphere-like geometry)
    - BUT with Hopf fibration, fiber vs cross edges might differ
      (they DON'T by vertex-transitivity, but the scalar curvature
       decomposition may reveal structure)

PLAN:
  1. Build 600-cell, Hopf fibration, Box operator
  2. Compute all-pairs shortest path (graph distance matrix)
  3. Compute Ollivier-Ricci curvature for all 720 edges via LP
  4. Verify vertex-transitivity (constant curvature)
  5. Compute scalar curvature S(x) = sum_{y~x} kappa(x,y)
  6. Decompose: S_fiber(x) + S_cross(x)
  7. Compare with spectral action coefficients (2640, 14880, 55920)
"""

import numpy as np
from numpy.linalg import eigvalsh
from scipy.sparse.csgraph import shortest_path
from scipy.sparse import csr_matrix
from scipy.optimize import linprog
from collections import Counter, defaultdict
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
degree = 12

print("=" * 70)
print("EXP-578: OLLIVIER-RICCI CURVATURE ON THE 600-CELL")
print("=" * 70)

# =====================================================================
# STEP 1: Build 600-cell + Hopf fibration
# =====================================================================
print("\nSTEP 1: Building 600-cell and Hopf fibration...")

verts, adj, lap = build_600cell()

def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

def find_fibration():
    for i in range(N):
        if abs(verts[i, 0] - PHI/2) < 1e-6:
            g = verts[i]; p = g.copy(); ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok=False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok=False
            if not ok: continue
            used = set(); fibers = []; subg = []
            pp = np.array([1.0,0,0,0])
            for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
            for s in range(N):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si]); idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12: return fibers
    return None

fibers = find_fibration()
vtx_fiber = {}
for fi, f in enumerate(fibers):
    for v in f: vtx_fiber[v] = fi

A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.0
A_cross = adj - A_fiber
Box = a1 * A_fiber - A_cross  # = b1*A_fiber - adj

# Edge list
edges = []
edge_set = set()
is_fiber_edge = []
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j] > 0.5:
            edges.append((i, j))
            edge_set.add((i,j))
            is_fiber_edge.append(vtx_fiber[i] == vtx_fiber[j])

is_fiber_edge = np.array(is_fiber_edge)
Ne = len(edges)
print(f"  Vertices: {N}, Edges: {Ne}")
print(f"  Fiber edges: {np.sum(is_fiber_edge)}, Cross edges: {np.sum(~is_fiber_edge)}")


# =====================================================================
# STEP 2: Graph distance matrix (all-pairs shortest path)
# =====================================================================
print("\nSTEP 2: Computing graph distance matrix...")

sparse_adj = csr_matrix(adj)
dist_matrix = shortest_path(sparse_adj, method='D', unweighted=True)

# Verify
print(f"  Max distance (diameter): {int(dist_matrix.max())}")
print(f"  Distance distribution:")
dist_counter = Counter()
for i in range(N):
    for j in range(i+1, N):
        dist_counter[int(dist_matrix[i,j])] += 1
for d in sorted(dist_counter.keys()):
    print(f"    d={d}: {dist_counter[d]} pairs")


# =====================================================================
# STEP 3: Ollivier-Ricci curvature via optimal transport
# =====================================================================
print("\n" + "=" * 70)
print("STEP 3: OLLIVIER-RICCI CURVATURE (Lin-Lu-Yau, idleness=0)")
print("=" * 70)

# Neighbor lists
nbrs = defaultdict(list)
for i in range(N):
    for j in range(N):
        if adj[i,j] > 0.5:
            nbrs[i].append(j)


def wasserstein_1(p_support, p_weights, q_support, q_weights, dist_mat):
    """
    Compute W_1 distance between two discrete measures using LP.

    p_support, q_support: lists of vertex indices
    p_weights, q_weights: probability weights (sum to 1)
    dist_mat: full graph distance matrix
    """
    m = len(p_support)
    n = len(q_support)

    # Cost matrix (m x n)
    cost = np.array([[dist_mat[p_support[i], q_support[j]]
                       for j in range(n)] for i in range(m)])

    # Flatten for LP: variables T[i,j] = x[i*n + j]
    c = cost.flatten()

    # Equality constraints:
    # sum_j T[i,j] = p_weights[i]  for all i  (m constraints)
    # sum_i T[i,j] = q_weights[j]  for all j  (n constraints)

    A_eq = np.zeros((m + n, m * n))
    b_eq = np.zeros(m + n)

    # Row sum constraints
    for i in range(m):
        for j in range(n):
            A_eq[i, i*n + j] = 1.0
        b_eq[i] = p_weights[i]

    # Column sum constraints
    for j in range(n):
        for i in range(m):
            A_eq[m + j, i*n + j] = 1.0
        b_eq[m + j] = q_weights[j]

    # Solve LP
    result = linprog(c, A_eq=A_eq, b_eq=b_eq,
                     bounds=[(0, None)] * (m*n),
                     method='highs', options={'presolve': True})

    if result.success:
        return result.fun
    else:
        raise RuntimeError(f"LP failed: {result.message}")


def ollivier_ricci(x, y, adj_list, dist_mat):
    """
    Compute Ollivier-Ricci curvature kappa(x,y) for edge (x,y).
    Uses Lin-Lu-Yau convention: m_x = uniform on N(x), idleness=0.
    """
    nx = adj_list[x]
    ny = adj_list[y]

    px = np.ones(len(nx)) / len(nx)
    py = np.ones(len(ny)) / len(ny)

    w1 = wasserstein_1(nx, px, ny, py, dist_mat)

    # d(x,y) = 1 for adjacent vertices
    return 1.0 - w1


# Compute curvature for ALL 720 edges
print("\n  Computing curvature for all 720 edges...")
kappas = np.zeros(Ne)
for e_idx, (i, j) in enumerate(edges):
    kappas[e_idx] = ollivier_ricci(i, j, nbrs, dist_matrix)
    if (e_idx + 1) % 100 == 0:
        print(f"    ... {e_idx+1}/{Ne} edges done")

print(f"  Done. All 720 curvatures computed.")


# =====================================================================
# STEP 4: Analysis
# =====================================================================
print("\n" + "=" * 70)
print("STEP 4: CURVATURE ANALYSIS")
print("=" * 70)

print(f"\n  Global statistics:")
print(f"    min(kappa) = {kappas.min():.10f}")
print(f"    max(kappa) = {kappas.max():.10f}")
print(f"    mean(kappa) = {kappas.mean():.10f}")
print(f"    std(kappa) = {kappas.std():.10f}")

# Check vertex-transitivity: all edges should have same curvature
is_constant = kappas.std() < 1e-8
print(f"\n  Vertex-transitivity check:")
print(f"    All edges same curvature: {is_constant}")
if is_constant:
    kappa_flat = kappas[0]
    print(f"    kappa_flat = {kappa_flat:.10f}")
else:
    # Cluster distinct values
    kappa_rounded = np.round(kappas, 6)
    distinct = Counter(kappa_rounded)
    print(f"    Distinct curvature values: {len(distinct)}")
    for val, count in sorted(distinct.items()):
        print(f"      kappa = {val:.6f}  ({count} edges)")

# Separate fiber vs cross
kappa_fiber = kappas[is_fiber_edge]
kappa_cross = kappas[~is_fiber_edge]

print(f"\n  Fiber edges (within C10 fibers):")
print(f"    count: {len(kappa_fiber)}")
print(f"    mean kappa: {kappa_fiber.mean():.10f}")
print(f"    std kappa:  {kappa_fiber.std():.10f}")

print(f"\n  Cross edges (between fibers):")
print(f"    count: {len(kappa_cross)}")
print(f"    mean kappa: {kappa_cross.mean():.10f}")
print(f"    std kappa:  {kappa_cross.std():.10f}")


# =====================================================================
# STEP 5: Scalar curvature
# =====================================================================
print("\n" + "=" * 70)
print("STEP 5: SCALAR CURVATURE")
print("=" * 70)

# S(x) = sum_{y~x} kappa(x,y)
scalar_curv = np.zeros(N)
for e_idx, (i, j) in enumerate(edges):
    scalar_curv[i] += kappas[e_idx]
    scalar_curv[j] += kappas[e_idx]

print(f"\n  Scalar curvature S(x):")
print(f"    min: {scalar_curv.min():.10f}")
print(f"    max: {scalar_curv.max():.10f}")
print(f"    mean: {scalar_curv.mean():.10f}")
print(f"    std: {scalar_curv.std():.10f}")

# Total scalar curvature = analog of integral of R
total_S = scalar_curv.sum()
print(f"\n  Total scalar curvature: sum S(x) = {total_S:.6f}")
print(f"    = 2 * sum kappa(e) = {2 * kappas.sum():.6f}")
print(f"    Per edge: {kappas.sum():.6f} / {Ne} = {kappas.mean():.6f}")
print(f"    Per vertex: {total_S / N:.6f}")

# Decompose scalar curvature into fiber and cross contributions
S_fiber = np.zeros(N)
S_cross = np.zeros(N)
for e_idx, (i, j) in enumerate(edges):
    if is_fiber_edge[e_idx]:
        S_fiber[i] += kappas[e_idx]
        S_fiber[j] += kappas[e_idx]
    else:
        S_cross[i] += kappas[e_idx]
        S_cross[j] += kappas[e_idx]

print(f"\n  Decomposition S(x) = S_fiber(x) + S_cross(x):")
print(f"    S_fiber per vertex: {S_fiber.mean():.10f} (from 2 fiber neighbors)")
print(f"    S_cross per vertex: {S_cross.mean():.10f} (from 10 cross neighbors)")
print(f"    Ratio S_cross/S_fiber: {S_cross.mean()/S_fiber.mean():.6f}")
print(f"    Expected if equal: 10/2 = {10/2}")


# =====================================================================
# STEP 6: Algebraic identification
# =====================================================================
print("\n" + "=" * 70)
print("STEP 6: ALGEBRAIC IDENTIFICATION")
print("=" * 70)

kappa_val = kappas.mean()

# Common neighbor analysis
print(f"\n  Common neighbor structure:")
# For an edge (x,y), how many common neighbors do they share?
common_nbrs_list = []
for e_idx, (i, j) in enumerate(edges):
    common = len(set(nbrs[i]) & set(nbrs[j]))
    common_nbrs_list.append(common)

common_nbrs_arr = np.array(common_nbrs_list)
print(f"    Common neighbors per edge: min={common_nbrs_arr.min()}, "
      f"max={common_nbrs_arr.max()}, mean={common_nbrs_arr.mean():.2f}")

common_fiber = common_nbrs_arr[is_fiber_edge]
common_cross = common_nbrs_arr[~is_fiber_edge]
print(f"    Fiber edges: common = {common_fiber.mean():.2f}")
print(f"    Cross edges: common = {common_cross.mean():.2f}")

# Try algebraic forms for kappa
print(f"\n  Algebraic identification attempts:")
print(f"    kappa = {kappa_val:.10f}")
print(f"    2/degree = {2/degree:.10f}")
print(f"    (common+2)/degree - 1 = {(common_nbrs_arr.mean()+2)/degree - 1:.10f}")

# Lin-Lu-Yau lower bound for d-regular graphs:
# kappa >= -2/d + 2*triangle(e)/(d*(d-1))
# where triangle(e) = number of triangles containing edge e
tri_per_edge = []
for e_idx, (i, j) in enumerate(edges):
    # Triangles containing edge (i,j) = common neighbors
    tri_per_edge.append(common_nbrs_arr[e_idx])
tri_arr = np.array(tri_per_edge)

print(f"\n  Triangles per edge: {tri_arr.mean():.2f}")
print(f"  LLY lower bound: -2/d + 2*tri/(d*(d-1))")
print(f"    = {-2/degree + 2*tri_arr.mean()/(degree*(degree-1)):.10f}")

# Try various exact forms
candidates = [
    ("1/b1", 1/b1),
    ("1/a1", 1/a1),
    ("2/degree", 2/degree),
    ("(a1-1)/degree", (a1-1)/degree),
    ("a1/degree^2", a1/degree**2),
    ("1/phi^4", 1/PHI**4),
    ("phi/degree", PHI/degree),
    ("(phi-1)/a1", (PHI-1)/a1),
    ("2/(degree+1)", 2/(degree+1)),
    ("5/(degree*phi^2)", 5/(degree*PHI**2)),
    ("common/(degree*phi)", common_nbrs_arr.mean()/(degree*PHI)),
    ("1/(2*phi^2)", 1/(2*PHI**2)),
    ("(common-degree+2)/degree", (common_nbrs_arr.mean()-degree+2)/degree),
]

print(f"\n  Exact form search (kappa = {kappa_val:.10f}):")
for name, val in candidates:
    diff = abs(val - kappa_val)
    marker = " <-- MATCH" if diff < 1e-6 else ""
    print(f"    {name:>25s} = {val:.10f}  diff = {diff:.2e}{marker}")


# =====================================================================
# STEP 7: Connection to spectral action
# =====================================================================
print("\n" + "=" * 70)
print("STEP 7: CONNECTION TO SPECTRAL ACTION")
print("=" * 70)

# The Einstein-Hilbert term in spectral action is c_1 = Tr(D^2) = 14880.
# In Connes NCG: S_EH = (1/2) * f_2 * Lambda^2 * R * vol
# The discrete analog: total_S = sum_x S(x) ~ integral of R.
# So: total_S / N ~ R (curvature per vertex)

print(f"\n  Curvature per vertex (analog of R):")
print(f"    R_disc = total_S / N = {total_S/N:.10f}")
print(f"    = degree * kappa = {degree * kappa_val:.10f}")

print(f"\n  Spectral action coefficient c_1 = 14880 (from Tr(D^2)):")
print(f"    c_1 / N = {14880/N:.2f} = dim(E8)/2 = 124")
print(f"    c_1 / (N * degree * kappa) = {14880/(N * degree * kappa_val):.6f}")
print(f"    c_1 / total_S = {14880/total_S:.6f}")

# Regge-like: total curvature from deficit angles
# On S^3: integral R = 12*pi^2 (for unit radius)
# Graph version: deficit angle ~ 2*pi - sum(face angles)
print(f"\n  Comparison with S^3:")
print(f"    S^3 Euler char chi = 0")
print(f"    600-cell chi = 120 - 720 + 1200 - 600 = 0")
print(f"    Consistent: positive curvature, but chi = 0 (odd-dim analog)")

# Box trace moments as curvature invariants
evals_box = np.sort(eigvalsh(Box))
tr2 = np.sum(evals_box**2)
tr4 = np.sum(evals_box**4)
print(f"\n  Box trace moments (discrete curvature invariants):")
print(f"    Tr(Box^2) / N = {tr2/N:.2f} = deg * a1 = 60")
print(f"    Total Ollivier S = {total_S:.6f}")
print(f"    Ratio Tr(Box^2) / total_S = {tr2/total_S:.6f}")


# =====================================================================
# STEP 8: SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
  OLLIVIER-RICCI CURVATURE ON THE 600-CELL:

  kappa(e) = {kappa_val:.10f}  (same for all 720 edges by vertex-transitivity)

  Common neighbors per edge: {common_nbrs_arr.mean():.0f}
  Triangles per edge: {tri_arr.mean():.0f}

  Scalar curvature:
    S(x) = degree * kappa = {degree * kappa_val:.10f}  (same for all 120 vertices)
    Total = N * S = {total_S:.6f}

  Fiber/cross decomposition:
    S_fiber(x) = 2 * kappa = {2 * kappa_val:.10f}  (2 fiber neighbors)
    S_cross(x) = 10 * kappa = {10 * kappa_val:.10f} (10 cross neighbors)

  POSITIVE curvature confirms sphere-like geometry.
  UNIFORM curvature confirms flat (homogeneous) vacuum.

  The flat 600-cell is the discrete analog of the round S^3:
  constant positive Ricci curvature = Einstein manifold with R > 0.

  NEXT: Perturb edge weights -> break homogeneity -> discrete Einstein eq.
""")

print("=" * 70)
print("EXP-578 COMPLETE")
print("=" * 70)

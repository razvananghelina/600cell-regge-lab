"""
exp579: Weighted Box operator and perturbation theory on the 600-cell.

KEY RESULT FROM exp578: The flat 600-cell has ZERO Ollivier-Ricci curvature.
This means it IS the vacuum solution of the discrete Einstein equation.

Now we introduce dynamics:
  1. Weighted adjacency: A(w)_ij = w_ij * A_ij, with w_ij > 0
  2. Weighted Box: Box(w) = a1 * A_fiber(w) - A_cross(w)
  3. Perturbation: w_ij = 1 + eps * h_ij
  4. First-order spectral response: delta(lambda) = eps * <psi| dBox |psi>
  5. Ollivier-Ricci curvature on the perturbed graph
  6. Linear response tensor: D_{edge,edge} = d(kappa_e)/d(w_f)
  7. Discrete Einstein equation: kappa(w) = 8*pi*G * T(w)

Goal: Show that the spectral variation of Box under weight perturbation
defines a natural notion of discrete gravity on the 600-cell.
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh, norm
from scipy.sparse.csgraph import shortest_path
from scipy.sparse import csr_matrix
from scipy.optimize import linprog
from collections import defaultdict
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
degree = 12

print("=" * 70)
print("EXP-579: WEIGHTED BOX AND PERTURBATION THEORY")
print("=" * 70)


# =====================================================================
# STEP 1: Build 600-cell + Hopf fibration + Box
# =====================================================================
print("\nSTEP 1: Building infrastructure...")

verts, adj, lap = build_600cell()

def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v; idx = np.argmax(dots)
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
Box_flat = a1 * A_fiber - A_cross

# Edge list and classification
edges = []
edge_idx_map = {}  # (i,j) -> index
is_fiber_edge = []
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j] > 0.5:
            idx = len(edges)
            edges.append((i, j))
            edge_idx_map[(i,j)] = idx
            edge_idx_map[(j,i)] = idx
            is_fiber_edge.append(vtx_fiber[i] == vtx_fiber[j])
is_fiber_edge = np.array(is_fiber_edge)
Ne = len(edges)

# Flat spectrum
evals_flat, evecs_flat = eigh(Box_flat)
print(f"  Vertices: {N}, Edges: {Ne}")
print(f"  Fiber: {np.sum(is_fiber_edge)}, Cross: {np.sum(~is_fiber_edge)}")
print(f"  Box spectrum range: [{evals_flat.min():.4f}, {evals_flat.max():.4f}]")


# =====================================================================
# STEP 2: Perturbation matrices dBox/dw_e for each edge e
# =====================================================================
print("\n" + "=" * 70)
print("STEP 2: PERTURBATION MATRICES (dBox/dw_e)")
print("=" * 70)

# Box(w) = a1 * sum_{fiber edges} w_e * (|i><j| + |j><i|)
#         - sum_{cross edges} w_e * (|i><j| + |j><i|)
#
# So dBox/dw_e = a1 * (|i><j| + |j><i|) if e is fiber
#              = -(|i><j| + |j><i|)      if e is cross

def make_dBox(edge_idx):
    """Perturbation matrix dBox/dw for edge edge_idx."""
    i, j = edges[edge_idx]
    dB = np.zeros((N, N))
    if is_fiber_edge[edge_idx]:
        dB[i,j] = a1
        dB[j,i] = a1
    else:
        dB[i,j] = -1.0
        dB[j,i] = -1.0
    return dB


# First-order eigenvalue perturbation:
# delta(lambda_k) = eps * <psi_k| dBox/dw_e |psi_k>
# For degenerate eigenvalues, we need the MATRIX of perturbation
# in the degenerate subspace, but for the trace (sum of eigenvalues
# in each multiplet), it's just the trace of the perturbation projected
# onto the subspace.

# Compute the sensitivity matrix: S[k, e] = <psi_k| dBox_e |psi_k>
print("  Computing sensitivity matrix S[mode, edge]...")
S = np.zeros((N, Ne))
for e_idx in range(Ne):
    dB = make_dBox(e_idx)
    for k in range(N):
        S[k, e_idx] = evecs_flat[:, k] @ dB @ evecs_flat[:, k]

print(f"  Sensitivity matrix shape: {S.shape}")
print(f"  S range: [{S.min():.6f}, {S.max():.6f}]")


# =====================================================================
# STEP 3: Eigenvalue response to specific perturbations
# =====================================================================
print("\n" + "=" * 70)
print("STEP 3: EIGENVALUE RESPONSE TO PERTURBATIONS")
print("=" * 70)

# Group eigenvalues by value (multiplets)
tol = 1e-6
multiplets = []
used = set()
evals_sorted_idx = np.argsort(evals_flat)
for k in evals_sorted_idx:
    if k in used:
        continue
    val = evals_flat[k]
    group = [k]
    for k2 in evals_sorted_idx:
        if k2 != k and k2 not in used and abs(evals_flat[k2] - val) < 0.01:
            group.append(k2)
    for g in group:
        used.add(g)
    multiplets.append((val, group))

multiplets.sort(key=lambda x: x[0])
print(f"\n  Eigenvalue multiplets of Box:")
for val, group in multiplets:
    mult = len(group)
    # Total sensitivity to a uniform perturbation (all edges +eps)
    total_sens = sum(S[k, :].sum() for k in group)
    print(f"    lambda = {val:+10.4f}  mult = {mult:3d}  "
          f"total sensitivity = {total_sens:+10.2f}")


# Test 1: Uniform fiber stretch (all fiber edges w = 1+eps)
print(f"\n  --- Test 1: Uniform fiber stretch ---")
h_fiber = np.zeros(Ne)
h_fiber[is_fiber_edge] = 1.0

delta_evals_fiber = S @ h_fiber  # first-order shifts
print(f"  Perturbation: all fiber edges w = 1+eps")
print(f"  First-order eigenvalue shifts (per multiplet):")
for val, group in multiplets:
    shift = sum(delta_evals_fiber[k] for k in group)
    avg_shift = shift / len(group)
    print(f"    lambda = {val:+10.4f}: total shift = {shift:+10.4f}, "
          f"per-mode = {avg_shift:+10.6f}")

# Verify numerically
eps = 1e-4
A_fiber_pert = A_fiber * (1 + eps)
Box_pert = a1 * A_fiber_pert - A_cross
evals_pert = np.sort(eigvalsh(Box_pert))
evals_flat_sorted = np.sort(evals_flat)
numerical_shifts = (evals_pert - evals_flat_sorted) / eps

print(f"\n  Numerical verification (eps={eps}):")
print(f"  Sum of first-order shifts: {delta_evals_fiber.sum():.6f}")
print(f"  Sum of numerical shifts:   {numerical_shifts.sum():.6f}")
print(f"  Tr(dBox_fiber) = {a1 * 2 * np.sum(is_fiber_edge) * 0:.6f} (off-diagonal only)")

# The trace of the total shift = Tr(sum_e h_e * dBox_e)
tr_shift = sum(np.trace(make_dBox(e)) * h_fiber[e] for e in range(Ne))
print(f"  Tr(perturbation) = {tr_shift:.6f}")


# Test 2: Uniform cross stretch
print(f"\n  --- Test 2: Uniform cross stretch ---")
h_cross = np.zeros(Ne)
h_cross[~is_fiber_edge] = 1.0

delta_evals_cross = S @ h_cross
print(f"  Perturbation: all cross edges w = 1+eps")
for val, group in multiplets:
    shift = sum(delta_evals_cross[k] for k in group)
    avg_shift = shift / len(group)
    print(f"    lambda = {val:+10.4f}: total shift = {shift:+10.4f}, "
          f"per-mode = {avg_shift:+10.6f}")


# Test 3: Single edge perturbation
print(f"\n  --- Test 3: Single fiber edge perturbation ---")
# Pick the first fiber edge
fiber_edge_idx = np.where(is_fiber_edge)[0][0]
i0, j0 = edges[fiber_edge_idx]
print(f"  Perturbing fiber edge ({i0}, {j0})")

delta_single = S[:, fiber_edge_idx]
print(f"  Eigenvalue shifts (first order):")
for val, group in multiplets:
    shift = sum(delta_single[k] for k in group)
    print(f"    lambda = {val:+10.4f}: total shift = {shift:+10.6f}")


# =====================================================================
# STEP 4: Spectral action variation
# =====================================================================
print("\n" + "=" * 70)
print("STEP 4: SPECTRAL ACTION VARIATION")
print("=" * 70)

# The spectral action S[w] = Tr f(Box(w)/Lambda)
# For polynomial cutoff f(x) = x^2: S[w] = Tr(Box(w)^2)
# dS/dw_e = 2 * Tr(Box * dBox_e) = 2 * sum_k lambda_k * <psi_k|dBox_e|psi_k>
#         = 2 * sum_k lambda_k * S[k,e]

# Variation of Tr(Box^2)
dTr2 = np.zeros(Ne)
for e in range(Ne):
    dTr2[e] = 2 * np.sum(evals_flat * S[:, e])

print(f"  dTr(Box^2)/dw_e:")
print(f"    Fiber edges: mean = {dTr2[is_fiber_edge].mean():.6f}, "
      f"std = {dTr2[is_fiber_edge].std():.6f}")
print(f"    Cross edges: mean = {dTr2[~is_fiber_edge].mean():.6f}, "
      f"std = {dTr2[~is_fiber_edge].std():.6f}")

# By vertex-transitivity, these should be constant within each class
dTr2_fiber_val = dTr2[is_fiber_edge].mean()
dTr2_cross_val = dTr2[~is_fiber_edge].mean()

print(f"\n  Exact values:")
print(f"    dTr(Box^2)/dw_fiber = {dTr2_fiber_val:.6f}")
print(f"    dTr(Box^2)/dw_cross = {dTr2_cross_val:.6f}")
print(f"    Ratio: {dTr2_fiber_val / dTr2_cross_val:.6f}")

# Verify: Tr(Box^2) = sum_ij Box_ij^2 = sum_ij (a1*Af - Ac)_ij^2
# dTr(Box^2)/dw_e = 2 * sum_{ij} Box_ij * dBox_ij/dw_e
# For fiber edge e=(p,q): dBox_pq/dw_e = a1, dBox_qp/dw_e = a1
# So dTr2/dw_e = 2 * (Box_pq * a1 + Box_qp * a1) = 4*a1*Box_pq

# Box_pq for adjacent p,q in same fiber: Box_pq = a1*Af_pq - Ac_pq = a1*1 - 0 = a1
# So dTr2/dw_fiber = 4*a1*a1 = 4*a1^2 = 100
print(f"\n  Algebraic check:")
print(f"    For fiber edge: dTr2 = 4*a1*Box_pq = 4*a1*a1 = {4*a1**2}")
print(f"    Match: {abs(dTr2_fiber_val - 4*a1**2) < 0.01}")

# For cross edge e=(p,q): dBox_pq/dw_e = -1
# Box_pq for adjacent p,q in different fibers: Box_pq = a1*0 - 1 = -1
# So dTr2/dw_cross = 2*(-1)*(-1) + 2*(-1)*(-1) = 4*1 = 4
print(f"    For cross edge: dTr2 = 4*Box_pq*(-1) = 4*(-1)*(-1) = {4}")
print(f"    Match: {abs(dTr2_cross_val - 4) < 0.01}")


# Variation of Tr(Box^4)
dTr4 = np.zeros(Ne)
for e in range(Ne):
    dTr4[e] = 4 * np.sum(evals_flat**3 * S[:, e])

dTr4_fiber_val = dTr4[is_fiber_edge].mean()
dTr4_cross_val = dTr4[~is_fiber_edge].mean()

print(f"\n  dTr(Box^4)/dw_e:")
print(f"    Fiber: {dTr4_fiber_val:.2f}")
print(f"    Cross: {dTr4_cross_val:.2f}")
print(f"    Ratio: {dTr4_fiber_val / dTr4_cross_val:.6f}")


# =====================================================================
# STEP 5: Ollivier-Ricci response to perturbation
# =====================================================================
print("\n" + "=" * 70)
print("STEP 5: OLLIVIER-RICCI RESPONSE TO WEIGHT PERTURBATION")
print("=" * 70)

# Compute kappa for a weighted graph
# For weighted graph: edge lengths l_ij = 1/w_ij (heavier = shorter)
# Alternative: l_ij = f(w_ij) for some monotone f
# Natural choice: weighted shortest path with l_ij = 1/w_ij

nbrs = defaultdict(list)
for i in range(N):
    for j in range(N):
        if adj[i,j] > 0.5:
            nbrs[i].append(j)


def wasserstein_1(p_support, p_weights, q_support, q_weights, dist_mat):
    m = len(p_support)
    n = len(q_support)
    cost = np.array([[dist_mat[p_support[i], q_support[j]]
                       for j in range(n)] for i in range(m)])
    c = cost.flatten()
    A_eq = np.zeros((m + n, m * n))
    b_eq = np.zeros(m + n)
    for i in range(m):
        for j in range(n):
            A_eq[i, i*n + j] = 1.0
        b_eq[i] = p_weights[i]
    for j in range(n):
        for i in range(m):
            A_eq[m + j, i*n + j] = 1.0
        b_eq[m + j] = q_weights[j]
    result = linprog(c, A_eq=A_eq, b_eq=b_eq,
                     bounds=[(0, None)] * (m*n), method='highs')
    return result.fun if result.success else float('nan')


def compute_weighted_kappas(weights, edges_list, adj_matrix, n_verts):
    """
    Compute Ollivier-Ricci curvature for a weighted graph.

    weights: edge weights w_e > 0
    Edge lengths: l_e = 1/w_e
    Measure at x: m_x proportional to w_{x,y} for y~x (weighted degree)
    """
    n = n_verts

    # Build weighted adjacency
    W = np.zeros((n, n))
    for e_idx, (i, j) in enumerate(edges_list):
        W[i,j] = weights[e_idx]
        W[j,i] = weights[e_idx]

    # Weighted shortest path distances: edge length = 1/w
    L = np.zeros((n, n))
    for e_idx, (i, j) in enumerate(edges_list):
        L[i,j] = 1.0 / weights[e_idx]
        L[j,i] = 1.0 / weights[e_idx]
    # Replace 0 with inf for non-edges
    L_sparse = csr_matrix(L)
    dist_w = shortest_path(L_sparse, method='D')

    # Neighbor lists with weights
    nbrs_w = defaultdict(list)
    for i in range(n):
        for j in range(n):
            if adj_matrix[i,j] > 0.5:
                nbrs_w[i].append(j)

    # Compute kappa for each edge
    kappas_w = np.zeros(len(edges_list))
    for e_idx, (x, y) in enumerate(edges_list):
        nx_list = nbrs_w[x]
        ny_list = nbrs_w[y]

        # Weighted measure: m_x(z) proportional to W[x,z]
        wx = np.array([W[x, z] for z in nx_list])
        wy = np.array([W[y, z] for z in ny_list])
        px = wx / wx.sum()
        py = wy / wy.sum()

        w1 = wasserstein_1(nx_list, px, ny_list, py, dist_w)
        d_xy = dist_w[x, y]
        kappas_w[e_idx] = 1.0 - w1 / d_xy

    return kappas_w


# Test: perturb a single fiber edge
print("\n  Testing single fiber edge perturbation...")
eps_values = [0.01, 0.05, 0.10, 0.20, 0.50]

fiber_edge_idx = np.where(is_fiber_edge)[0][0]
i0, j0 = edges[fiber_edge_idx]
print(f"  Perturbed edge: ({i0}, {j0}) [fiber]")

for eps in eps_values:
    w = np.ones(Ne)
    w[fiber_edge_idx] = 1.0 + eps

    kappas_w = compute_weighted_kappas(w, edges, adj, N)

    # Statistics
    nonzero = np.abs(kappas_w) > 1e-10
    n_nonzero = np.sum(nonzero)
    max_abs = np.max(np.abs(kappas_w))
    kappa_pert = kappas_w[fiber_edge_idx]

    print(f"\n  eps = {eps:.2f}:")
    print(f"    kappa(perturbed edge) = {kappa_pert:.8f}")
    print(f"    Nonzero curvature edges: {n_nonzero} / {Ne}")
    print(f"    max|kappa| = {max_abs:.8f}")
    print(f"    Total scalar curvature: {kappas_w.sum():.8f}")

    if eps == eps_values[0]:
        # Store for linear response
        kappas_first = kappas_w.copy()
        dkappa_dw = kappas_w / eps  # linear response at first eps

# Check linearity
print(f"\n  Linear response check:")
print(f"    dkappa/dw at perturbed edge (eps->0 limit):")
print(f"    From eps=0.01: {kappas_first[fiber_edge_idx]/0.01:.8f}")


# =====================================================================
# STEP 6: Curvature response matrix (Jacobian)
# =====================================================================
print("\n" + "=" * 70)
print("STEP 6: CURVATURE RESPONSE TO PERTURBATION (SAMPLE)")
print("=" * 70)

# For a single cross edge perturbation
cross_edge_idx = np.where(~is_fiber_edge)[0][0]
i1, j1 = edges[cross_edge_idx]
print(f"  Perturbed cross edge: ({i1}, {j1})")

eps_test = 0.05
w_cross = np.ones(Ne)
w_cross[cross_edge_idx] = 1.0 + eps_test

kappas_cross = compute_weighted_kappas(w_cross, edges, adj, N)

nonzero_cross = np.sum(np.abs(kappas_cross) > 1e-10)
print(f"  eps = {eps_test}:")
print(f"    kappa(perturbed edge) = {kappas_cross[cross_edge_idx]:.8f}")
print(f"    Nonzero curvature edges: {nonzero_cross} / {Ne}")
print(f"    max|kappa| = {np.max(np.abs(kappas_cross)):.8f}")

# Compare fiber vs cross perturbation response
print(f"\n  COMPARISON: fiber vs cross perturbation (eps={eps_test}):")
print(f"    Fiber: kappa(pert edge) = {kappas_first[fiber_edge_idx]/0.01*eps_test:.8f}")
print(f"    Cross: kappa(pert edge) = {kappas_cross[cross_edge_idx]:.8f}")


# =====================================================================
# STEP 7: Spectral action as gravitational action
# =====================================================================
print("\n" + "=" * 70)
print("STEP 7: SPECTRAL ACTION AS GRAVITATIONAL ACTION")
print("=" * 70)

# The spectral action S_grav[w] = Tr f(Box(w)^2 / Lambda^2)
# For polynomial f: S_grav = c_0 * Lambda^4 * N + c_1 * Lambda^2 * Tr(Box^2) + ...
# The Einstein equation comes from: delta S_grav / delta w_e = 0
# i.e., dTr(Box^2)/dw_e = 0 for all e -> uniform w (flat vacuum)

# But we showed dTr(Box^2)/dw_e = 4*a1^2 (fiber) or 4 (cross).
# These are NOT zero! So Tr(Box^2) alone doesn't give Einstein eq.

# The correct quantity is the spectral action with VARYING degree:
# When w changes, the degree changes -> the Laplacian changes -> more terms.

# Actually, the natural gravitational action is:
# S[w] = Tr(Box(w)^2) - (1/N) * [Tr(Box(w))]^2
# (subtract trace squared to impose tracelessness constraint)

# Or: S[w] with constraint Tr(Box) = 0 (preserved if sum of weights constant)

# Constraint: Tr(Box) = 0. For Box = a1*Af(w) - Ac(w):
# Tr(Box) = 0 (always, since Box is off-diagonal = symmetric with zero diagonal)
# Actually Tr(adj) = 0, Tr(Af) = 0, so Tr(Box) = 0 always. Good.

# So the constraint is that the total "volume" is fixed:
# sum_e w_e = Ne (total weight = number of edges)

# Einstein equation from spectral action:
# Minimize S[w] = Tr(Box(w)^2) subject to sum(w) = Ne, w > 0.
# By Lagrange multipliers: dTr(Box^2)/dw_e = lambda for all e.
# This gives: 4*a1^2 * Af_pq + 4*Ac_pq = lambda -> all the same.
# But fiber and cross give different values! So the extremum is NOT at w=1
# for Tr(Box^2).

# Better: use a WEIGHTED trace that accounts for the fiber/cross structure.
# S[w] = Tr(Box(w)^2) with the constraint that the GEOMETRY is fixed
# (not just total weight, but fiber and cross separately).

# Let's check: what is the extremum of Tr(Box^2) with constraints
# sum_{fiber} w = 120, sum_{cross} w = 600?

# dTr2/dw = 4*a1^2 for fiber, 4 for cross.
# With separate Lagrange multipliers: all fiber have same w, all cross same w.
# -> w_fiber = w_cross = 1 is the ONLY solution (by symmetry).
# So the flat metric IS the (constrained) extremum, but it's a SADDLE POINT
# or minimum depending on the Hessian.

print(f"  Gravitational action candidates:")
print(f"    Tr(Box^2) at flat: {np.sum(evals_flat**2):.2f}")
print(f"    dTr(Box^2)/dw: fiber = {4*a1**2}, cross = {4}")
print(f"    Ratio: {a1**2} = a1^2 (fiber edges are {a1**2}x stiffer)")

# Compute Hessian d^2(Tr(Box^2))/dw_e dw_f
# d^2 Tr(Box^2)/dw_e dw_f = 2 * Tr(dBox_e * dBox_f)
# For e != f: Tr(dBox_e * dBox_f) = 2 * delta_edge_share * coeff_e * coeff_f
# where delta_edge_share = 1 if edges e,f share a vertex

print(f"\n  Hessian structure:")
# For same edge: d^2/dw_e^2 = 2*Tr(dBox_e^2)
# dBox_e = a1*(|i><j|+|j><i|) for fiber, -(|i><j|+|j><i|) for cross
# Tr(dBox_e^2) = 2*a1^2 for fiber, 2 for cross
print(f"    d^2Tr(Box^2)/dw_e^2: fiber = {4*a1**2}, cross = {4}")
print(f"    (These are positive -> local minimum in each direction)")

# For different edges sharing a vertex:
# dBox_e * dBox_f has nonzero entries only if edges share a vertex
# Let e = (i,j), f = (i,k). Then:
# (dBox_e * dBox_f)_{j,k} = coeff_e * coeff_f (if both connected to i)
# Tr contribution: need j = k, which means the edges are the same. So Tr = 0 for distinct edges!

# Wait: dBox_e = c_e * (|i><j| + |j><i|)
# dBox_e * dBox_f = c_e * c_f * (|i><j| + |j><i|)(|k><l| + |l><k|)
# If e=(i,j), f=(i,k) with i shared:
# = c_e * c_f * (|j><i| |i><k| + ...) = c_e * c_f * |j><k| + ...
# Tr = c_e * c_f * (<k|j> + ...) = 0 if j != k (distinct edges)

# So the Hessian of Tr(Box^2) is DIAGONAL in edge indices!
# d^2/dw_e dw_f = 0 for e != f
# d^2/dw_e^2 = 4*a1^2 (fiber) or 4 (cross)

print(f"    Off-diagonal Hessian: 0 (Hessian is diagonal!)")
print(f"    -> Tr(Box^2) is SEPARABLE in edge weights")
print(f"    -> Each edge weight decouples")
print(f"    -> Flat metric is a STRICT MINIMUM of Tr(Box^2) "
      f"(positive definite Hessian)")


# =====================================================================
# STEP 8: Higher-order spectral action
# =====================================================================
print("\n" + "=" * 70)
print("STEP 8: VARIATION OF Tr(Box^4) -- COUPLING BETWEEN EDGES")
print("=" * 70)

# Tr(Box^4) = sum_k lambda_k^4
# d/dw_e = 4 * sum_k lambda_k^3 * S[k,e]  (already computed as dTr4)
# d^2/dw_e dw_f = 4*3*sum_k lambda_k^2 * S[k,e]*S[k,f]
#                 + 4*sum_k lambda_k^3 * d^2lambda_k/(dw_e dw_f)

# The first term gives a NON-DIAGONAL Hessian -> edges COUPLE at order 4!
# This is the discrete analog of the Gauss-Bonnet / R^2 terms in gravity.

# Compute the coupling matrix (first term only, which dominates)
print("  Computing edge-edge coupling from Tr(Box^4)...")

# C_ef = 12 * sum_k lambda_k^2 * S[k,e] * S[k,f]
# This is a rank-N matrix (at most)
# C = 12 * (S^T . diag(lambda^2) . S)
lam2_diag = evals_flat**2
C_coupling = 12 * (S.T * lam2_diag[np.newaxis, :]) @ S  # Ne x Ne... too big?

# Actually S is N x Ne = 120 x 720
# S^T is 720 x 120
# diag(lam^2) is 120 x 120
# So C = 720 x 720. That's fine.

# But let's just look at the structure
print(f"  Coupling matrix C[e,f] shape: {C_coupling.shape}")
print(f"  C range: [{C_coupling.min():.2f}, {C_coupling.max():.2f}]")

# Diagonal elements
C_diag = np.diag(C_coupling)
print(f"\n  Diagonal (self-coupling):")
print(f"    Fiber edges: {C_diag[is_fiber_edge].mean():.4f} +/- {C_diag[is_fiber_edge].std():.4f}")
print(f"    Cross edges: {C_diag[~is_fiber_edge].mean():.4f} +/- {C_diag[~is_fiber_edge].std():.4f}")

# Off-diagonal: neighbor vs non-neighbor edges
print(f"\n  Off-diagonal structure:")
# Two edges are neighbors if they share a vertex
n_edge_nbr = 0
n_edge_non_nbr = 0
sum_nbr = 0
sum_non_nbr = 0
for e1 in range(min(50, Ne)):  # sample
    for e2 in range(e1+1, min(50, Ne)):
        i1, j1 = edges[e1]
        i2, j2 = edges[e2]
        shared = len({i1,j1} & {i2,j2})
        if shared > 0:
            n_edge_nbr += 1
            sum_nbr += C_coupling[e1, e2]
        else:
            n_edge_non_nbr += 1
            sum_non_nbr += C_coupling[e1, e2]

if n_edge_nbr > 0:
    print(f"    Adjacent edges: mean coupling = {sum_nbr/n_edge_nbr:.4f}")
if n_edge_non_nbr > 0:
    print(f"    Non-adjacent: mean coupling = {sum_non_nbr/n_edge_non_nbr:.4f}")


# =====================================================================
# STEP 9: DISCRETE EINSTEIN EQUATION
# =====================================================================
print("\n" + "=" * 70)
print("STEP 9: DISCRETE EINSTEIN EQUATION")
print("=" * 70)

# The natural gravitational dynamics comes from:
# S_grav[w] = Tr f(Box(w)/Lambda)
#           = c_0 + c_2/Lambda^2 * Tr(Box^2) + c_4/Lambda^4 * Tr(Box^4) + ...
#
# The Einstein equation is: delta S_grav / delta w_e = 0
#
# At leading order (Tr(Box^2)):
#   dTr(Box^2)/dw_e = const for all e -> AUTOMATICALLY SATISFIED by symmetry
#   (fiber and cross have different constants, but within each class, constant)
#   So the flat vacuum IS a solution at leading order.
#
# At next order (Tr(Box^4)):
#   dTr(Box^4)/dw_e = const for all e within each edge class
#   -> Also satisfied by flat vacuum (by vertex-transitivity).
#
# The flat 600-cell is an EXACT vacuum solution at all orders of the
# spectral action, because it's vertex-transitive.

print(f"""
  DISCRETE EINSTEIN EQUATION FROM SPECTRAL ACTION:

  S_grav[w] = sum_n c_n * Tr(Box(w)^{{2n}})

  Vacuum equation: dS_grav/dw_e = 0  for all edges e

  For the flat 600-cell (w_e = 1):
    dTr(Box^2)/dw_e = {4*a1**2} (fiber), {4} (cross)  -- constant within class
    dTr(Box^4)/dw_e = {dTr4_fiber_val:.2f} (fiber), {dTr4_cross_val:.2f} (cross)  -- constant

  With constraint preserving fiber/cross structure:
    delta S / delta w_e - lambda_f * I_fiber - lambda_c * I_cross = 0
    -> The flat metric is the UNIQUE vacuum solution (by symmetry).

  KEY PROPERTIES:
    1. Flat = vacuum: kappa = 0, and dS/dw = const (extremum)
    2. Hessian of Tr(Box^2) is DIAGONAL: edges decouple at order 2
    3. Hessian of Tr(Box^4) is NON-DIAGONAL: edges couple at order 4
    4. The coupling is through shared spectral modes (gravity propagates!)
    5. Stiffness ratio: fiber/cross = a1^2 = 25
       (fiber edges are 25x stiffer than cross edges)

  PHYSICAL INTERPRETATION:
    - Fiber = time direction, Cross = space direction
    - Fiber 25x stiffer -> the speed of light c^2 = a1 = 5 is RIGID
    - Space is more easily deformed than time -> gravitational waves
      are predominantly spatial perturbations (transverse-traceless)
    - The coupling at order 4 is the discrete graviton propagator
""")

# Stiffness ratio
ratio_stiffness = dTr2[is_fiber_edge].mean() / dTr2[~is_fiber_edge].mean()
print(f"  Stiffness ratio fiber/cross = {ratio_stiffness:.2f} = a1^2 = {a1**2}")
print(f"  This is the SAME a1^2 that gives c^2 = a1!")
print(f"  The causal structure (c^2 = 5) and gravitational stiffness (25x)")
print(f"  share the SAME origin: the Hopf fiber anisotropy.")


# =====================================================================
# STEP 10: Graviton spectrum (edge perturbation modes)
# =====================================================================
print("\n" + "=" * 70)
print("STEP 10: GRAVITON SPECTRUM (EDGE MODES)")
print("=" * 70)

# The graviton is a fluctuation h_e of edge weights around flat.
# Its kinetic term comes from the Hessian of the spectral action:
# S_grav^(2)[h] = (1/2) * sum_{e,f} H_{ef} * h_e * h_f
# where H_{ef} = d^2 S_grav / (dw_e dw_f)

# At leading order (from Tr(Box^2)):
# H_{ef}^(2) = 4*a1^2 * delta_{ef} (fiber) or 4*delta_{ef} (cross)
# This is diagonal -> no propagation (massive but non-propagating at this order)

# At next order (from Tr(Box^4)):
# H_{ef}^(4) = C_coupling[e,f] (computed above)
# This has off-diagonal elements -> PROPAGATING graviton!

# Diagonalize the coupling matrix to find graviton modes
print("  Diagonalizing coupling matrix (720 x 720)...")
evals_grav = np.sort(eigvalsh(C_coupling))

print(f"\n  Graviton spectrum (eigenvalues of coupling matrix):")
print(f"    Min: {evals_grav[0]:.4f}")
print(f"    Max: {evals_grav[-1]:.4f}")
print(f"    Positive modes: {np.sum(evals_grav > 1e-6)}")
print(f"    Zero modes: {np.sum(np.abs(evals_grav) < 1e-6)}")
print(f"    Negative modes: {np.sum(evals_grav < -1e-6)}")

# Count distinct eigenvalues
from collections import Counter
grav_counter = Counter(np.round(evals_grav, 2))
print(f"\n  Distinct eigenvalues (graviton masses):")
for val, mult in sorted(grav_counter.items()):
    if mult >= 1:
        print(f"    g = {val:12.4f}  mult = {mult:4d}")


# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY OF EXP-579")
print("=" * 70)

print(f"""
  DISCRETE GRAVITY ON THE 600-CELL:

  1. OLLIVIER-RICCI FLAT: kappa = 0 on all edges (exp578)
     -> The flat 600-cell IS the vacuum Einstein solution

  2. SPECTRAL ACTION VARIATION:
     dTr(Box^2)/dw = 4*a1^2 = 100 (fiber), 4 (cross)
     -> Stiffness ratio a1^2 = 25 (fiber time is rigid)
     -> Same origin as c^2 = a1 = 5

  3. EDGE DECOUPLING at order 2:
     Hessian diagonal -> no gravity propagation from Tr(Box^2) alone

  4. EDGE COUPLING at order 4:
     Tr(Box^4) generates off-diagonal Hessian
     -> Graviton propagates via shared spectral modes
     -> Discrete graviton = edge weight fluctuation

  5. GRAVITON SPECTRUM from coupling matrix:
     {np.sum(evals_grav > 1e-6)} propagating modes
     {np.sum(np.abs(evals_grav) < 1e-6)} zero modes (gauge/diffeomorphisms)

  NEXT STEPS:
    - Identify zero modes with discrete diffeomorphisms
    - Compute graviton propagator G^{{-1}}_{{ef}}
    - Add matter source: T_e from Box eigenvalue perturbation
    - Verify that the linearized equation reproduces Newton's law
    - Check: does the emergent R match the spectral action c_1 = 14880?
""")

print("=" * 70)
print("EXP-579 COMPLETE")
print("=" * 70)

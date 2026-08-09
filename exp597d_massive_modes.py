"""
exp100d: Massive modes analysis and final cross-checks.

KEY RESULT FROM exp100c:
  ker(Box1_edge) = 13 = rho_0 + 2*rho_5 = 1 + 2*6
  rho_5 is the BRANCHING NODE of E8 Dynkin diagram (dim 6).

Now investigate:
1. What are the first massive eigenvalues? Do they relate to W/Z masses?
2. Irrep decomposition of each eigenspace
3. Do the eigenvalues relate to Z[phi]?
4. Is there a natural 1+12 = 1 + (8+3+1) gauge boson structure in 2*rho_5?
5. What happens with the 600 cross-edge dominated modes?
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
# REBUILD
# =====================================================================
print("Building...")
verts, adj, lap = build_600cell()
Nv = 120

adj_list = defaultdict(set)
edges = []
edge_to_idx = {}
for i in range(Nv):
    for j in range(i+1, Nv):
        if adj[i, j] > 0.5:
            adj_list[i].add(j)
            adj_list[j].add(i)
            edges.append((i, j))
            edge_to_idx[(i, j)] = len(edges) - 1
Ne = len(edges)

triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    triangles.append((i, j, k))

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
    target_w = PHI / 2.0
    for i in range(Nv):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            p = g.copy()
            ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok = False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok = False
            if not ok: continue
            used = set(); fibers = []; subg = []
            pp = np.array([1.0,0,0,0])
            for k in range(10):
                subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
            for s in range(Nv):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si])
                    idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12: return fibers
    return None

fibers = find_fibration()
vertex_to_fiber = {}
for fi, f in enumerate(fibers):
    for v in f: vertex_to_fiber[v] = fi

is_fiber_edge = np.zeros(Ne, dtype=bool)
for e_idx, (i, j) in enumerate(edges):
    if vertex_to_fiber[i] == vertex_to_fiber[j]:
        is_fiber_edge[e_idx] = True

vertex_to_edges = defaultdict(list)
for e_idx, (i, j) in enumerate(edges):
    vertex_to_edges[i].append(e_idx)
    vertex_to_edges[j].append(e_idx)

A_line = np.zeros((Ne, Ne))
for v in range(Nv):
    inc = vertex_to_edges[v]
    for a in range(len(inc)):
        for b in range(a+1, len(inc)):
            A_line[inc[a], inc[b]] = 1.0
            A_line[inc[b], inc[a]] = 1.0

A_fiber_line = np.zeros((Ne, Ne))
for fi, fiber in enumerate(fibers):
    fe = []
    for k in range(10):
        u, v = fiber[k], fiber[(k+1)%10]
        fe.append(edge_to_idx[(min(u,v), max(u,v))])
    for k in range(10):
        A_fiber_line[fe[k], fe[(k+1)%10]] = 1.0
        A_fiber_line[fe[(k+1)%10], fe[k]] = 1.0

A_cross_line = A_line - A_fiber_line
L_fiber_line = np.diag(np.sum(A_fiber_line, axis=1)) - A_fiber_line
L_cross_line = np.diag(np.sum(A_cross_line, axis=1)) - A_cross_line
Box1 = L_cross_line - a1 * L_fiber_line

# Vertex Box for reference
A_fiber_v = adj * np.array([[float(vertex_to_fiber[i] == vertex_to_fiber[j])
                for j in range(Nv)] for i in range(Nv)])
L_fiber_v = np.diag(A_fiber_v.sum(axis=1)) - A_fiber_v
L_cross_v = np.diag((adj - A_fiber_v).sum(axis=1)) - (adj - A_fiber_v)
Box_v = L_cross_v - a1 * L_fiber_v

print("Done.\n")


# =====================================================================
# PART 1: FULL SPECTRUM OF BOX1 — ALGEBRAIC STRUCTURE
# =====================================================================
print("=" * 72)
print("PART 1: FULL SPECTRUM — ALGEBRAIC STRUCTURE")
print("=" * 72)

evals_Box1, evecs_Box1 = eigh(Box1)
evals_v, evecs_v = eigh(Box_v)

# Vertex Box eigenvalues
ev_v = Counter(np.round(evals_v, 4))
print("\n  Vertex Box eigenvalues:")
for val, mult in sorted(ev_v.items()):
    # Try to identify in Z[phi]
    best_ab = None
    for b in range(-20, 21):
        a = val - b * PHI
        if abs(a - round(a)) < 0.005:
            best_ab = (int(round(a)), b)
            break
    label = f"  = {best_ab[0]} + {best_ab[1]}*phi" if best_ab else ""
    print(f"    {val:10.4f} x {mult:3d}{label}")

# Edge Box1 eigenvalues
ev_e = Counter(np.round(evals_Box1, 4))
print(f"\n  Edge Box1 eigenvalues ({len(ev_e)} distinct):")
total_modes = 0
for val, mult in sorted(ev_e.items()):
    total_modes += mult
    # Try to identify in Z[phi]
    best_ab = None
    for b in range(-30, 31):
        a = val - b * PHI
        if abs(a - round(a)) < 0.01:
            best_ab = (int(round(a)), b)
            break
    # Also try a + b*sqrt(5)
    best_abs5 = None
    for b2 in range(-30, 31):
        a2 = val - b2 * SQRT5
        if abs(a2 - round(a2)) < 0.01:
            best_abs5 = (int(round(a2)), b2)
            break
    label = ""
    if best_ab:
        label = f"  = {best_ab[0]} + {best_ab[1]}*phi"
    elif best_abs5:
        label = f"  = {best_abs5[0]} + {best_abs5[1]}*sqrt5"
    print(f"    {val:10.4f} x {mult:3d}  (cumul {total_modes:3d}){label}")


# =====================================================================
# PART 2: IRREP DECOMPOSITION OF EACH EIGENSPACE
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: IRREP DECOMPOSITION OF EIGENSPACES")
print("=" * 72)

# Build group action tools
all_vperm = []
for g in range(Nv):
    perm = np.zeros(Nv, dtype=int)
    for i in range(Nv):
        product = qmul(verts[g], verts[i])
        perm[i] = find_idx(product, verts)
    all_vperm.append(perm)

all_eperm = []
for g in range(Nv):
    eperm = np.zeros(Ne, dtype=int)
    for e_idx, (i, j) in enumerate(edges):
        gi, gj = all_vperm[g][i], all_vperm[g][j]
        key = (min(gi, gj), max(gi, gj))
        eperm[e_idx] = edge_to_idx[key]
    all_eperm.append(eperm)

# Build proper conjugacy classes using vertex Laplacian eigenspaces
evals_lap, evecs_lap = eigh(lap)
ev_spaces = {}
for i in range(Nv):
    val = round(evals_lap[i], 4)
    if val not in ev_spaces: ev_spaces[val] = []
    ev_spaces[val].append(i)
eigenspace_list = sorted(ev_spaces.items())

char_vectors = np.zeros((Nv, 9))
for g in range(Nv):
    for k, (val, indices) in enumerate(eigenspace_list):
        space = evecs_lap[:, indices]
        chi = sum(np.dot(space[:, col], space[:, col][all_vperm[g]]) for col in range(space.shape[1]))
        char_vectors[g, k] = chi

class_dict = defaultdict(list)
for g in range(Nv):
    key = tuple(np.round(char_vectors[g], 4))
    class_dict[key].append(g)

class_list = sorted(class_dict.items(), key=lambda x: (-len(x[1]), x[0]))
class_sizes = np.array([len(m) for _, m in class_list], dtype=float)
class_reps = [m[0] for _, m in class_list]

# True character table (normalized)
irrep_dims = [int(round(np.sqrt(len(eigenspace_list[k][1])))) for k in range(9)]
chi_true = np.zeros((9, len(class_list)))
for ri in range(9):
    for ci, (key, _) in enumerate(class_list):
        chi_true[ri, ci] = key[ri] / irrep_dims[ri]

print(f"  Irrep dims: {irrep_dims}")
print(f"  Irrep labels: rho_0(1), rho_1(2), rho_2(3), rho_3(4), rho_4(5), rho_5(6), rho_6(3'), rho_7(4'), rho_8(2')")

# For each eigenspace of Box1, compute irrep multiplicities
print(f"\n  {'eigenvalue':>10s} {'mult':>5s}  Decomposition")
print(f"  {'-'*70}")

for val, mult in sorted(ev_e.items()):
    # Get eigenvectors at this eigenvalue
    mask = np.abs(evals_Box1 - val) < 0.01
    space_vecs = evecs_Box1[:, mask]
    dim = space_vecs.shape[1]

    # Compute character on each class
    ev_chars = np.zeros(len(class_list))
    for ci, (_, members) in enumerate(class_list):
        g = members[0]
        eperm = all_eperm[g]
        tr = sum(np.dot(space_vecs[:, k], space_vecs[:, k][eperm]) for k in range(dim))
        ev_chars[ci] = tr

    # Decompose
    decomp = []
    for ri in range(9):
        n = np.sum(class_sizes * chi_true[ri] * ev_chars) / 120.0
        if abs(n) > 0.01:
            decomp.append((ri, irrep_dims[ri], round(n)))

    decomp_str = " + ".join(f"{n}*rho_{ri}({d})" for ri, d, n in decomp if n > 0)
    dim_check = sum(d * n for _, d, n in decomp)
    ok = "OK" if dim_check == mult else f"MISMATCH ({dim_check})"
    print(f"  {val:10.4f} {mult:5d}  {decomp_str}  [{ok}]")


# =====================================================================
# PART 3: GAUGE STRUCTURE IN 2*rho_5
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: GAUGE STRUCTURE IN THE 12-DIM (2*rho_5) SECTOR")
print("=" * 72)

# The kernel has rho_0(1) + 2*rho_5(6) = 13.
# Question: does the 12-dimensional 2*rho_5 subspace further decompose
# in a way that reveals 8+3+1 = SU(3)+SU(2)+U(1)?

# The A5 decomposition on vertices: 12 = 1 + 3 + 8 decomposes the
# vertex adjacency into gauge sectors.
# Here we ask: within the kernel of Box1_edge, restricted to the 2*rho_5
# subspace, is there a natural splitting?

# Get the 13-dim kernel
ker_mask = np.abs(evals_Box1) < 1e-6
ker_vecs = evecs_Box1[:, ker_mask]
ker_dim = 13

# Project out the rho_0 component (it's the unique trivial rep mode)
# The trivial rep character: chi_0(g) = 1 for all g
# Projector onto rho_0: P_0 = (1/|G|) sum_g R(g)
P_rho0 = np.zeros((ker_dim, ker_dim))
for g in range(Nv):
    eperm = all_eperm[g]
    for ki in range(ker_dim):
        gv = ker_vecs[:, ki][eperm]
        for kj in range(ker_dim):
            P_rho0[ki, kj] += np.dot(ker_vecs[:, kj], gv)
P_rho0 /= 120.0

# Verify P_rho0 is rank 1
p0_evals = eigvalsh(P_rho0)
print(f"  P_rho0 eigenvalues (sorted): {np.sort(p0_evals)[-5:]}")
print(f"  P_rho0 rank: {np.sum(np.abs(p0_evals) > 0.01)}")

# The rho_5 projector (complementary in kernel)
P_rho5 = np.eye(ker_dim) - P_rho0

# Analyze the 12-dim rho_5 subspace
# Get an orthonormal basis for this subspace
p5_evals, p5_evecs = eigh(P_rho5)
rho5_basis_ker = p5_evecs[:, p5_evals > 0.5]  # 12 vectors in 13-dim kernel space
print(f"  rho_5 subspace dimension: {rho5_basis_ker.shape[1]}")

# Lift to edge space
rho5_edge = ker_vecs @ rho5_basis_ker  # (720 x 12)

# Now analyze the rho5 sector using the fiber structure.
# Fiber content of each rho5 mode
print(f"\n  Fiber/cross content of rho_5 modes:")
for k in range(12):
    v = rho5_edge[:, k]
    ff = np.sum(v[is_fiber_edge]**2)
    print(f"    mode {k}: fiber = {ff:.4f}, cross = {1-ff:.4f}")

# Compute the LINE GRAPH adjacency restricted to the kernel
A_ker = rho5_edge.T @ A_line @ rho5_edge  # (12 x 12)
A_fiber_ker = rho5_edge.T @ A_fiber_line @ rho5_edge  # (12 x 12)
A_cross_ker = rho5_edge.T @ A_cross_line @ rho5_edge  # (12 x 12)

print(f"\n  A_line restricted to rho_5 kernel:")
print(f"    Eigenvalues: {np.sort(eigvalsh(A_ker))}")
print(f"\n  A_fiber restricted to rho_5 kernel:")
print(f"    Eigenvalues: {np.sort(eigvalsh(A_fiber_ker))}")


# =====================================================================
# PART 4: VERTEX BOX EIGENVALUE CORRESPONDENCE
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: VERTEX-EDGE EIGENVALUE CORRESPONDENCE")
print("=" * 72)

# Compare vertex and edge eigenvalues
# Vertex: 13 distinct eigenvalues
# Edge: 30 distinct eigenvalues
# Is there a map?

print("\n  Vertex Box spectrum:")
for val, mult in sorted(ev_v.items()):
    print(f"    {val:10.4f} x {mult}")

print(f"\n  Edge eigenvalue 1.5020:")
val = 1.5020
# Check: is this related to vertex eigenvalues?
for vval, vmult in sorted(ev_v.items()):
    if abs(vval) < 0.001:
        continue
    ratio = val / vval
    diff = val - vval
    if abs(diff) < 2:
        print(f"    ratio to vertex {vval:.4f}: {ratio:.6f}, diff: {diff:.6f}")

# Check: is 1.5020 = 3-phi? = 3 - 1.618 = 1.382? No.
# Is it 2-1/phi = 2 - 0.618 = 1.382? No.
# Is it phi^2 - phi = phi = 1.618? No.
# Let me try exact identification
test_val = sorted(ev_e.items())[1][0]  # first nonzero eigenvalue
print(f"\n  First nonzero eigenvalue: {test_val:.10f}")
# Test various algebraic forms
candidates = {
    "3-phi": 3 - PHI,
    "phi^2-phi": PHI**2 - PHI,
    "2-1/phi": 2 - 1/PHI,
    "3/2": 1.5,
    "3/phi": 3/PHI,
    "5-2phi": 5 - 2*PHI,
    "phi/phi": 1.0,
    "3-sqrt(5)/2": 3 - SQRT5/2,
    "7-4phi+2": 7 - 4*PHI + 2,
    "(12-6phi)/2": (12-6*PHI)/2,
    "2*(12-6phi)/(14-2)": 2*(12-6*PHI)/12,
    "12-6phi-phi^2+phi": 12-6*PHI-PHI**2+PHI,
    "(7-4phi)*phi+1": (7-4*PHI)*PHI+1,
}
for name, cv in candidates.items():
    if abs(test_val - cv) < 0.01:
        print(f"    MATCH: {name} = {cv:.10f} (diff = {test_val-cv:.2e})")


# =====================================================================
# PART 5: CHARACTERISTIC POLYNOMIAL CHECK (exact eigenvalue)
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: EXACT EIGENVALUE IDENTIFICATION")
print("=" * 72)

# Box1 is integer-valued. Its eigenvalues are algebraic integers.
# For the first nonzero eigenvalue cluster (multiplicity 10):
# The minimal polynomial should have integer coefficients.

# Eigenvalue ~1.5020: test if it satisfies various polynomial equations
lam = test_val
print(f"  Testing lambda = {lam:.12f}")

# Quadratic: lam^2 - p*lam + q = 0 (for p,q integer)
for p in range(-10, 30):
    q_test = lam**2 - p*lam
    if abs(q_test - round(q_test)) < 0.001:
        q = round(q_test)
        print(f"    Quadratic: lambda^2 - {p}*lambda + {q} = 0  (residual = {lam**2 - p*lam + q:.2e})")

# Also check quartic for non-rational algebraic
print(f"\n  lambda^2 = {lam**2:.10f}")
print(f"  lambda^3 = {lam**3:.10f}")
print(f"  lambda^4 = {lam**4:.10f}")

# Check: is lam a root of x^2 - 3x + q?
# x^2 - 3x + q = 0 => q = -(lam^2 - 3*lam)
q_3 = -(lam**2 - 3*lam)
print(f"  If x^2 - 3x + q = 0: q = {q_3:.10f}")
q_4 = -(lam**2 - 4*lam)
print(f"  If x^2 - 4x + q = 0: q = {q_4:.10f}")

# For each distinct eigenvalue, find its minimal polynomial (degree ≤ 4)
print(f"\n  Minimal polynomial identification for all eigenvalues:")
for val, mult in sorted(ev_e.items()):
    found = False
    # Integer?
    if abs(val - round(val)) < 0.001:
        print(f"    {val:10.4f} x {mult:3d}: INTEGER = {int(round(val))}")
        continue
    # Quadratic: x^2 - px + q = 0 with integer p,q
    for p in range(-50, 50):
        q = val**2 - p*val
        if abs(q - round(q)) < 0.001:
            q_int = int(round(q))
            disc = p**2 - 4*q_int
            if disc >= 0:
                root_str = f"({p} ± sqrt({disc}))/2"
            else:
                root_str = f"({p} ± i*sqrt({-disc}))/2"
            print(f"    {val:10.4f} x {mult:3d}: x^2 - {p}x + {q_int} = 0  => {root_str}")
            found = True
            break
    if not found:
        print(f"    {val:10.4f} x {mult:3d}: no simple quadratic found")


# =====================================================================
# PART 6: MASS RATIOS FROM EIGENVALUES
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: MASS RATIOS")
print("=" * 72)

# If eigenvalues represent mass^2 (in some units), the ratios matter.
# First nonzero eigenvalues: these are the "massive" modes.
# Expected: W±, Z masses => m_W/m_Z = cos(theta_W) = sqrt(20/26)

nonzero_evals = sorted(set(round(e, 4) for e in evals_Box1 if abs(e) > 0.01))
print(f"  First 10 nonzero eigenvalues:")
for i, val in enumerate(nonzero_evals[:10]):
    mult = ev_e[val]
    print(f"    lambda_{i} = {val:10.4f} x {mult}")

# Ratios between consecutive eigenvalues
if len(nonzero_evals) >= 3:
    print(f"\n  Ratios:")
    for i in range(min(5, len(nonzero_evals)-1)):
        r = nonzero_evals[i+1] / nonzero_evals[i] if nonzero_evals[i] != 0 else float('inf')
        print(f"    lambda_{i+1}/lambda_{i} = {nonzero_evals[i+1]:.4f}/{nonzero_evals[i]:.4f} = {r:.6f}")

# Expected mass ratios
cos2_tW = 20.0/26.0  # = 1 - sin^2(theta_W)
print(f"\n  Expected m_W^2/m_Z^2 = cos^2(theta_W) = 20/26 = {cos2_tW:.6f}")
print(f"  Expected m_Z^2/m_W^2 = 26/20 = {26/20:.6f}")

# Check all pairs of small eigenvalues for this ratio
for i in range(min(8, len(nonzero_evals))):
    for j in range(i+1, min(8, len(nonzero_evals))):
        r = nonzero_evals[j] / nonzero_evals[i]
        if abs(r - 26/20) < 0.05 or abs(r - 20/26) < 0.05:
            mi, mj = ev_e[nonzero_evals[i]], ev_e[nonzero_evals[j]]
            print(f"  *** {nonzero_evals[j]:.4f}/{nonzero_evals[i]:.4f} = {r:.6f} " +
                  f"(mults {mi}, {mj})")


# =====================================================================
# PART 7: THE NUMBER 13 AND THE FRAMEWORK
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: THE NUMBER 13 IN THE FRAMEWORK")
print("=" * 72)

print(f"""
  ker(Box_edge) = 13

  Framework connections:
  - 13 = (a1^2 + 1)/2 = 26/2
  - 26 = a1^2 + 1 = denominator of sin^2(theta_W) = 6/26
  - 13 = 1 + 12 = trivial + gauge generators
  - 12 = dim SU(3)xSU(2)xU(1) = 8 + 3 + 1
  - 12 = 2 * dim(rho_5) where rho_5 is branching node of E8
  - 12 = vertex degree of 600-cell
  - 13 is NOT a dimension of any 2I irrep

  Decomposition: rho_0 + 2*rho_5 = 1 + 12
  - rho_0: the constant mode (would be projected out in physics)
  - 2*rho_5: the 12 gauge boson modes
  - rho_5 is the node of E8 where the gauge symmetry breaks!

  Comparison with vertex Box:
  - ker(Box_vertex) = 9 = rho_0 + 2*rho_1 + 2*rho_8
  - ker(Box_edge) = 13 = rho_0 + 2*rho_5

  The vertex kernel uses Galois pairs (rho_1, rho_8): dim 2+2 = 4 -> 2*4=8+1=9
  The edge kernel uses the branching representation rho_5: dim 6 -> 2*6=12+1=13

  Trace ratios:
  - Tr(Box_edge^2) / Tr(Box_vertex^2) = 45 = 9 * a1 = 9 * 5
  - Note: 45 = dim(ker_vertex) * a1
""")


# =====================================================================
# PART 8: NEAR-ZERO STRUCTURE — ARE 3 MODES SEPARATED?
# =====================================================================
print("=" * 72)
print("PART 8: DETAILED NEAR-ZERO STRUCTURE (W/Z SEARCH)")
print("=" * 72)

# The plan predicted: 9 exactly massless (gluons+photon) + 3 near-zero (W±,Z)
# We found 13 exactly massless. But maybe the 13 -> 12+1 -> 9+3+1?
# Look at the structure within the kernel more carefully.

# Within the 13-dim kernel, compute the VERTEX Box operator projected
# Box_v lifts to edge space via d0: check Box_v_lift restricted to kernel
d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = +1.0

Box_v_lift = d0 @ Box_v @ d0.T  # 720 x 720
Box_ker = ker_vecs.T @ Box_v_lift @ ker_vecs  # 13 x 13

print(f"\n  Box_vertex_lift restricted to ker(Box_edge): 13x13 matrix")
bk_evals = np.sort(eigvalsh(Box_ker))
print(f"  Eigenvalues: {np.round(bk_evals, 4)}")

# Count distinct
bk_distinct = Counter(np.round(bk_evals, 4))
print(f"  Distinct:")
for val, mult in sorted(bk_distinct.items()):
    print(f"    {val:10.4f} x {mult}")

# Also: Hodge Laplacian restricted to kernel
# Need to rebuild d1
d1 = np.zeros((len(triangles), Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

B = d0 @ d0.T
C = d1.T @ d1
Delta_1 = B + C

B_ker = ker_vecs.T @ B @ ker_vecs
C_ker = ker_vecs.T @ C @ ker_vecs
D1_ker = ker_vecs.T @ Delta_1 @ ker_vecs

print(f"\n  Hodge B (exact) restricted to kernel:")
print(f"  Eigenvalues: {np.round(np.sort(eigvalsh(B_ker)), 4)}")

print(f"\n  Hodge C (coexact) restricted to kernel:")
print(f"  Eigenvalues: {np.round(np.sort(eigvalsh(C_ker)), 4)}")

print(f"\n  Hodge Delta_1 restricted to kernel:")
print(f"  Eigenvalues: {np.round(np.sort(eigvalsh(D1_ker)), 4)}")


print("\n" + "=" * 72)
print("EXP100d COMPLETE")
print("=" * 72)

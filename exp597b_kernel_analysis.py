"""
exp100b: Deep analysis of the 13-dimensional kernel of Box1 on line graph.

Key questions:
1. Is ker = 13 exact or numerical artifact (could it be 12+1)?
2. How does the kernel decompose under 2I symmetry?
3. What is the physical interpretation?
4. Does removing the trivial constant mode give 12 = 1+3+8?

Also: investigate Box3's eigenvalue 2.0 with mult 480, and
      the clean trace ratios found in exp100.
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
# REBUILD EVERYTHING (from exp100)
# =====================================================================
print("Building 600-cell and all operators...")
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
            idx = len(edges)
            edges.append((i, j))
            edge_to_idx[(i, j)] = idx
Ne = len(edges)

triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    triangles.append((i, j, k))
Nf = len(triangles)

d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = +1.0

d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

B = d0 @ d0.T
C = d1.T @ d1

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
                if k == 5 and not np.allclose(p, [-1, 0, 0, 0], atol=1e-6):
                    ok = False; break
                if k == 10 and not np.allclose(p, [1, 0, 0, 0], atol=1e-6):
                    ok = False
            if not ok: continue
            used = set()
            fibers = []
            subg = []
            pp = np.array([1.0, 0, 0, 0])
            for k in range(10):
                idx = find_idx(pp, verts)
                subg.append(idx)
                pp = qmul(pp, g)
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
vertex_to_fiber = {}
for fi, f in enumerate(fibers):
    for v in f:
        vertex_to_fiber[v] = fi

is_fiber_edge = np.zeros(Ne, dtype=bool)
for e_idx, (i, j) in enumerate(edges):
    if vertex_to_fiber[i] == vertex_to_fiber[j]:
        is_fiber_edge[e_idx] = True

P_fiber = np.diag(is_fiber_edge.astype(float))
P_cross = np.eye(Ne) - P_fiber

vertex_to_edges = defaultdict(list)
for e_idx, (i, j) in enumerate(edges):
    vertex_to_edges[i].append(e_idx)
    vertex_to_edges[j].append(e_idx)

A_line = np.zeros((Ne, Ne), dtype=float)
for v in range(Nv):
    inc_edges = vertex_to_edges[v]
    for a in range(len(inc_edges)):
        for bb in range(a+1, len(inc_edges)):
            A_line[inc_edges[a], inc_edges[bb]] = 1.0
            A_line[inc_edges[bb], inc_edges[a]] = 1.0

A_fiber_line = np.zeros((Ne, Ne), dtype=float)
for fi, fiber in enumerate(fibers):
    fiber_edges_in_fiber = []
    for k in range(10):
        u, v = fiber[k], fiber[(k+1) % 10]
        i, j = (min(u, v), max(u, v))
        e_idx = edge_to_idx[(i, j)]
        fiber_edges_in_fiber.append(e_idx)
    for k in range(10):
        e1 = fiber_edges_in_fiber[k]
        e2 = fiber_edges_in_fiber[(k+1) % 10]
        A_fiber_line[e1, e2] = 1.0
        A_fiber_line[e2, e1] = 1.0

A_cross_line = A_line - A_fiber_line
D_fiber_line = np.diag(np.sum(A_fiber_line, axis=1))
D_cross_line = np.diag(np.sum(A_cross_line, axis=1))
L_fiber_line = D_fiber_line - A_fiber_line
L_cross_line = D_cross_line - A_cross_line

Box1 = L_cross_line - a1 * L_fiber_line

print("All operators built.\n")


# =====================================================================
# ANALYSIS 1: PRECISE KERNEL EIGENVALUES OF BOX1
# =====================================================================
print("=" * 72)
print("ANALYSIS 1: PRECISE KERNEL OF BOX1")
print("=" * 72)

evals_Box1, evecs_Box1 = eigh(Box1)
sorted_idx = np.argsort(np.abs(evals_Box1))

# Show the 20 smallest |eigenvalues|
print("\n  20 smallest |eigenvalues| of Box1:")
for k in range(20):
    i = sorted_idx[k]
    print(f"    [{k:2d}] lambda = {evals_Box1[i]:16.12f}")

# How many are truly zero (threshold scan)?
print("\n  Kernel size vs threshold:")
for thr in [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2]:
    nk = np.sum(np.abs(evals_Box1) < thr)
    print(f"    |lambda| < {thr:.0e}: ker = {nk}")


# =====================================================================
# ANALYSIS 2: KERNEL VECTORS — DECOMPOSITION BY FIBER
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 2: KERNEL EIGENVECTORS OF BOX1")
print("=" * 72)

ker_mask = np.abs(evals_Box1) < 1e-6
ker_vecs = evecs_Box1[:, ker_mask]
ker_dim = ker_vecs.shape[1]
print(f"\n  Kernel dimension: {ker_dim}")

# For each kernel vector, show support on each of 12 fibers
# First, build fiber edge index
fiber_edge_indices = {}  # fiber_id -> list of edge indices
for fi, fiber in enumerate(fibers):
    eidxs = []
    for k in range(10):
        u, v = fiber[k], fiber[(k+1) % 10]
        i, j = (min(u, v), max(u, v))
        eidxs.append(edge_to_idx[(i, j)])
    fiber_edge_indices[fi] = eidxs

print("\n  Kernel vector support on each fiber (fiber edge weight):")
print(f"  {'mode':>4s}", end="")
for fi in range(12):
    print(f"  {'F'+str(fi):>6s}", end="")
print(f"  {'cross':>8s}")
print("  " + "-" * 90)

for k in range(ker_dim):
    v = ker_vecs[:, k]
    print(f"  {k:4d}", end="")
    for fi in range(12):
        weight = np.sum(v[fiber_edge_indices[fi]]**2)
        print(f"  {weight:6.3f}", end="")
    cross_weight = np.sum(v[~is_fiber_edge]**2)
    print(f"  {cross_weight:8.3f}")


# =====================================================================
# ANALYSIS 3: 2I SYMMETRY ACTION ON EDGE SPACE
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 3: 2I SYMMETRY ON EDGE SPACE")
print("=" * 72)

# The 2I group acts on vertices by left quaternion multiplication.
# This induces an action on edges: g(i,j) = (g*i, g*j).
# Build the 120 permutation matrices on the edge space.

# First, build the action on vertices
def vertex_perm(g_idx):
    """Permutation of vertices under left multiplication by verts[g_idx]."""
    g = verts[g_idx]
    perm = np.zeros(Nv, dtype=int)
    for i in range(Nv):
        product = qmul(g, verts[i])
        idx = find_idx(product, verts)
        assert idx >= 0, f"Vertex {i} mapped outside 600-cell by g={g_idx}"
        perm[i] = idx
    return perm

# Build the induced edge permutation
def edge_perm(vperm):
    """Given a vertex permutation, compute the induced edge permutation."""
    eperm = np.zeros(Ne, dtype=int)
    for e_idx, (i, j) in enumerate(edges):
        gi, gj = vperm[i], vperm[j]
        key = (min(gi, gj), max(gi, gj))
        assert key in edge_to_idx, f"Edge ({gi},{gj}) not found"
        eperm[e_idx] = edge_to_idx[key]
    return eperm

# Build the representation matrices on the kernel
# Use a few representative group elements to compute characters

# For character computation, we need Tr(g|_kernel) for each group element.
# This equals sum of kernel vector dot products under g permutation.

print("\n  Computing 2I action on kernel (this takes a moment)...")

# Group into conjugacy classes of 2I
# 2I has 9 conjugacy classes with sizes: 1, 1, 12, 12, 20, 20, 30, 12, 12
# Characters of the 9 irreps (dims 1,2,3,4,5,6,4,2,3):
# We'll compute the character of the kernel representation

# Compute character for all 120 group elements
characters_ker = np.zeros(Nv)
for g in range(Nv):
    vperm = vertex_perm(g)
    eperm = edge_perm(vperm)
    # Trace of g on kernel = sum_k <v_k, g*v_k>
    tr = 0.0
    for k in range(ker_dim):
        v = ker_vecs[:, k]
        gv = v[eperm]  # Permuted vector
        tr += np.dot(v, gv)
    characters_ker[g] = tr

# Find conjugacy classes by computing the character table signature
# Group elements with the same character value on ALL irreps are in the same class
# But we only have one representation. Let's use the trace of adj (vertex character)
# to classify conjugacy classes.

# Alternative: classify by the ORDER of the group element
# Using powers of the vertex permutation
def perm_order(perm):
    """Order of a permutation."""
    current = perm.copy()
    identity = np.arange(len(perm))
    for k in range(1, 121):
        if np.array_equal(current, identity):
            return k
        current = perm[current]
    return -1

print("  Computing conjugacy classes...")
orders = []
for g in range(Nv):
    vperm = vertex_perm(g)
    orders.append(perm_order(vperm))

# Group by (order, character on vertices)
vertex_chars = np.zeros(Nv)
for g in range(Nv):
    vperm = vertex_perm(g)
    vertex_chars[g] = np.sum(np.arange(Nv) == vperm)  # = number of fixed points

# Conjugacy classes: group by (order, fixed_points)
class_key = {}
for g in range(Nv):
    key = (orders[g], int(vertex_chars[g]))
    if key not in class_key:
        class_key[key] = []
    class_key[key].append(g)

print(f"\n  Conjugacy class summary (order, fixed_pts): size, chi_kernel")
print(f"  {'order':>5s} {'fix_pts':>7s} {'size':>5s} {'chi_kernel':>10s}")
print(f"  {'-'*35}")
for key in sorted(class_key.keys()):
    members = class_key[key]
    chi_vals = [characters_ker[g] for g in members]
    # All members should have same character (verify)
    chi_std = np.std(chi_vals)
    chi_avg = np.mean(chi_vals)
    print(f"  {key[0]:5d} {key[1]:7d} {len(members):5d} {chi_avg:10.4f}" +
          (f"  (std={chi_std:.6f})" if chi_std > 1e-6 else ""))

# Expected 2I character table (for reference):
# Irrep dims: 1, 2, 3_a, 3_b, 4_a, 4_b, 5, 6, (but actually 2I has 9 classes)
# The 9 irreps of 2I have dimensions: 1, 2, 3, 4, 5, 6, 4', 2', 3'
# These correspond to the nodes of the extended E8 Dynkin diagram

# Irrep dimensions of 2I (binary icosahedral group)
irrep_dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]  # rho_0 through rho_8
print(f"\n  2I irrep dimensions: {irrep_dims}")
print(f"  Sum of dim^2 = {sum(d**2 for d in irrep_dims)} (expect 120 = |2I|)")

# Character table of 2I (standard, 9 classes)
# Classes with sizes: 1, 1, 12, 12, 20, 20, 30, 12, 12
# Orders: 1, 2, 10, 10, 6, 6, 4, 10, 10
# Class labels: {1}, {-1}, {g5}, {g5^2}, {g3}, {-g3}, {g4}, {-g5}, {-g5^2}

# Using standard 2I character table (chi values)
# [class size, order, chi_1, chi_2, chi_3, chi_4, chi_5, chi_6, chi_4', chi_2', chi_3']
phi = PHI
ct = np.array([
    # 1,  2,    3,       4,       5,    6,    4',     2',      3'
    [1,   2,    3,       4,       5,    6,    4,      2,       3      ],  # class {1}, size 1
    [1,  -2,    3,      -4,       5,   -6,    4,     -2,       3      ],  # class {-1}, size 1
    [1,   phi,  1+phi,   1+2*phi, 2+2*phi, 2+3*phi, 2*phi-1, phi-1,  phi   ],  # class C5, size 12
    [1,   phi-1, phi,    2*phi-1, 2*phi, 3*phi-1,   1+phi,   phi,    1+phi ],  # class C5^2, size 12
    [1,   1,    0,       1,      -1,    0,   -1,     -1,       0      ],  # class C3, size 20
    [1,  -1,    0,      -1,      -1,    0,    1,      1,       0      ],  # class -C3, size 20
    [1,   0,   -1,       0,       1,    0,    0,      0,      -1      ],  # class C4, size 30
    [1,  1-phi, phi-2, -1-phi, 2+2*phi-4, 1-3*phi, 2*phi-1, 1-phi,  phi  ],  # class -C5, size 12
    [1,  -phi+1, -phi+1, 1-2*phi, 2-2*phi, -3*phi+2, 1+phi-2, phi-1, -phi+2  ],  # class -C5^2, size 12
])
# Hmm, this character table is getting complicated. Let me just use the numerical result.

# Instead of the analytic character table, let me decompose using:
# n_rho = (1/|G|) * sum_g chi_rho(g)* chi_kernel(g)
# I need the full character table. Let me compute it from the vertex representation.

# The VERTEX representation of 2I is the 120-dimensional regular representation.
# Actually no, 2I acts on the 120 vertices, but the vertex representation is NOT
# the regular representation (that would be 120-dimensional with every irrep
# appearing dim(rho) times).

# Actually, the 120-dimensional permutation representation on vertices IS related
# to the regular representation because 2I acts freely on itself (left multiplication).
# So the vertex permutation representation IS the regular representation!
# Therefore: chi_vertex = sum_rho dim(rho) * chi_rho

# Let me compute the character of the edge representation for each conjugacy class
# and then decompose.

print("\n\n  Computing FULL characters for edge representation...")
chars_edge = {}  # (order, fix_pts) -> chi_edge_full

# For the edge representation: count fixed edges under each g
for key, members in class_key.items():
    g = members[0]  # representative
    vperm = vertex_perm(g)
    eperm = edge_perm(vperm)
    n_fixed = np.sum(np.arange(Ne) == eperm)
    chars_edge[key] = n_fixed

print(f"\n  Edge representation characters (number of fixed edges):")
print(f"  {'order':>5s} {'fix_pts':>7s} {'size':>5s} {'chi_edge':>8s} {'chi_ker':>10s}")
print(f"  {'-'*45}")
for key in sorted(class_key.keys()):
    print(f"  {key[0]:5d} {key[1]:7d} {len(class_key[key]):5d} {chars_edge[key]:8d} {np.mean([characters_ker[g] for g in class_key[key]]):10.4f}")


# =====================================================================
# ANALYSIS 4: DECOMPOSE KERNEL USING PROJECTION FORMULA
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 4: IRREP DECOMPOSITION OF KERNEL")
print("=" * 72)

# Build all 120 representation matrices on the kernel space
# P_rho = (dim_rho / |G|) * sum_g chi_rho(g)* * R(g)
# where R(g) is the matrix of g acting on the kernel

print("  Building 120 kernel representation matrices...")
R_ker = np.zeros((Nv, ker_dim, ker_dim))
for g in range(Nv):
    vperm = vertex_perm(g)
    eperm = edge_perm(vperm)
    # R(g) on kernel: R_ker[g, i, j] = <ker_i | g | ker_j>
    for ki in range(ker_dim):
        gv = ker_vecs[:, ki][eperm]  # g acting on kernel vector ki
        for kj in range(ker_dim):
            R_ker[g, ki, kj] = np.dot(ker_vecs[:, kj], gv)

# Verify: R is a representation (check R(g1)*R(g2) = R(g1*g2) for a few)
g1, g2 = 1, 2
vperm1 = vertex_perm(g1)
vperm2 = vertex_perm(g2)
vperm12 = vertex_perm(find_idx(qmul(verts[g1], verts[g2]), verts))
prod_check = R_ker[g1] @ R_ker[g2]
g12 = find_idx(qmul(verts[g1], verts[g2]), verts)
err_rep = norm(prod_check - R_ker[g12])
print(f"  Representation check: ||R(g1)R(g2) - R(g1*g2)|| = {err_rep:.2e}")

# Compute multiplicities using character inner product
# First need proper conjugacy class structure
# n_rho = (1/120) * sum_g chi_rho(g)* * chi_kernel(g)
# = (1/120) * sum_classes |C| * chi_rho(rep)* * chi_kernel(rep)

# Since we have all 120 characters already:
# n_rho = (1/120) * sum_g Tr(R_rho(g))* * chi_kernel(g)
# But we need the individual irrep characters.

# Approach: use the fact that the vertex representation = regular representation
# Its character: chi_reg(1) = 120, chi_reg(g) = 0 for g != 1
# And: reg = sum_rho dim(rho) * rho

# We can get irrep characters from the Laplacian eigenspaces!
# The Laplacian eigenvalues correspond to irreps. Each eigenspace is
# an irrep (or sum of copies of an irrep).

# The 9 distinct Laplacian eigenvalues on vertices correspond to the 9 irreps.
# For each eigenspace, the character is the trace of g restricted to that eigenspace.

print("\n  Computing irrep characters from vertex Laplacian eigenspaces...")
evals_v, evecs_v = eigh(lap)

# Find distinct eigenvalues
ev_distinct = []
ev_to_space = {}
for i in range(Nv):
    val = round(evals_v[i], 6)
    if val not in ev_to_space:
        ev_to_space[val] = []
        ev_distinct.append(val)
    ev_to_space[val].append(i)

print(f"  {len(ev_distinct)} distinct eigenvalues, multiplicities:")
irrep_spaces = []  # list of (eigenvalue, multiplicity, eigenvectors)
for val in sorted(ev_distinct):
    indices = ev_to_space[val]
    mult = len(indices)
    space = evecs_v[:, indices]
    irrep_spaces.append((val, mult, space))
    print(f"    lambda = {val:10.4f}, mult = {mult:3d}")

# For each eigenspace (irrep), compute its character on each group element
# chi_rho(g) = Tr(P_rho * Pi_g) where Pi_g is the permutation matrix on vertices

print("\n  Computing irrep characters for each conjugacy class...")
irrep_chars = {}  # (irrep_index, class_key) -> chi value

# Build vertex permutation matrices (sparse: just use permutation)
# chi_rho(g) = sum_j |P_rho * e_{sigma(j)}|^2 ... no
# chi_rho(g) = Tr(g on eigenspace) = sum of <v_k, g * v_k> for basis {v_k}

for ri, (eigval, mult, space) in enumerate(irrep_spaces):
    for key, members in class_key.items():
        g = members[0]
        vperm = vertex_perm(g)
        chi = 0.0
        for k in range(mult):
            v = space[:, k]
            gv = v[vperm]
            chi += np.dot(v, gv)
        irrep_chars[(ri, key)] = chi

# Print character table
print(f"\n  Character table (rows = irreps, cols = conjugacy classes):")
keys_sorted = sorted(class_key.keys())
print(f"  {'irrep':>5s} {'dim':>4s}", end="")
for key in keys_sorted:
    print(f"  {f'({key[0]},{key[1]})':>8s}", end="")
print()
print(f"  {'-'*80}")

for ri, (eigval, mult, _) in enumerate(irrep_spaces):
    print(f"  rho_{ri:<2d} {mult:4d}", end="")
    for key in keys_sorted:
        chi = irrep_chars[(ri, key)]
        print(f"  {chi:8.3f}", end="")
    print()

# Now decompose kernel representation
print(f"\n  Kernel character:")
print(f"  {'':>11s}", end="")
for key in keys_sorted:
    chi = np.mean([characters_ker[g] for g in class_key[key]])
    print(f"  {chi:8.3f}", end="")
print()

# n_rho = (1/|G|) * sum_C |C| * chi_rho(C)* * chi_kernel(C)
print(f"\n  Irrep multiplicities in kernel:")
decomposition = []
for ri, (eigval, mult, _) in enumerate(irrep_spaces):
    n_rho = 0.0
    for key, members in class_key.items():
        size = len(members)
        chi_rho = irrep_chars[(ri, key)]
        chi_ker = np.mean([characters_ker[g] for g in members])
        n_rho += size * chi_rho * chi_ker  # real characters, no conjugate needed
    n_rho /= 120.0
    decomposition.append((ri, mult, eigval, n_rho))
    if abs(n_rho) > 0.001:
        print(f"    rho_{ri} (dim {mult}, lambda={eigval:.4f}): multiplicity = {n_rho:.4f}")

total_dim = sum(mult * round(n) for ri, mult, _, n in decomposition if abs(n) > 0.001)
print(f"\n  Total dimension check: {total_dim} (expect {ker_dim})")


# =====================================================================
# ANALYSIS 5: IS THERE A NATURAL 12+1 SPLITTING?
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 5: 12+1 SPLITTING INVESTIGATION")
print("=" * 72)

# Check if one of the 13 kernel modes is "special" (e.g., constant on fibers)
# A constant mode on the line graph would be a vector with all entries equal
const_vec = np.ones(Ne) / np.sqrt(Ne)
projections = []
for k in range(ker_dim):
    proj = abs(np.dot(ker_vecs[:, k], const_vec))
    projections.append(proj)
    if proj > 0.1:
        print(f"  Mode {k}: projection onto constant vector = {proj:.6f}")

print(f"  Max projection onto constant: {max(projections):.6f}")

# Check projection onto fiber-constant vector (constant on each fiber, zero on cross)
fiber_const = np.zeros(Ne)
fiber_const[is_fiber_edge] = 1.0 / np.sqrt(120)
for k in range(ker_dim):
    proj = abs(np.dot(ker_vecs[:, k], fiber_const))
    if proj > 0.1:
        print(f"  Mode {k}: projection onto fiber-constant = {proj:.6f}")

# Check each fiber's constant mode
for fi in range(12):
    fvec = np.zeros(Ne)
    for e_idx in fiber_edge_indices[fi]:
        fvec[e_idx] = 1.0 / np.sqrt(10)
    for k in range(ker_dim):
        proj = abs(np.dot(ker_vecs[:, k], fvec))
        if proj > 0.3:
            print(f"  Mode {k}: projection onto fiber {fi} constant = {proj:.4f}")


# =====================================================================
# ANALYSIS 6: BOX1 SHIFTED — SUBTRACT CONSTANT DIAGONAL
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 6: SHIFTED BOX1 (TRACELESS VERSION)")
print("=" * 72)

# Box1 has Tr = 14400 = N^2. The diagonal is not uniform:
# fiber edges: 10, cross edges: 22
# Try making it traceless: Box1_shifted = Box1 - (Tr/N_e)*I
shift = np.trace(Box1) / Ne
print(f"  Tr(Box1)/Ne = {shift:.4f}")
Box1_shifted = Box1 - shift * np.eye(Ne)
evals_shifted = np.sort(eigvalsh(Box1_shifted))
ker_shifted = np.sum(np.abs(evals_shifted) < 1e-6)
print(f"  ker(Box1 - {shift:.2f}*I) = {ker_shifted}")

# The vertex Box also isn't traceless. Let's check:
A_fiber_v = adj * np.array([[float(vertex_to_fiber[i] == vertex_to_fiber[j])
                for j in range(Nv)] for i in range(Nv)])
L_fiber_v = np.diag(A_fiber_v.sum(axis=1)) - A_fiber_v
L_cross_v = np.diag((adj - A_fiber_v).sum(axis=1)) - (adj - A_fiber_v)
Box_v = L_cross_v - a1 * L_fiber_v
tr_Bv = np.trace(Box_v)
print(f"  Tr(Box_vertex) = {tr_Bv:.4f}")
print(f"  Tr(Box_vertex)/Nv = {tr_Bv/Nv:.4f}")

# Try different shifts
for s_val in [10, 20, shift]:
    evals_s = np.sort(eigvalsh(Box1 - s_val * np.eye(Ne)))
    ker_s = np.sum(np.abs(evals_s) < 1e-6)
    if ker_s > 0:
        print(f"  Box1 - {s_val:.2f}*I: ker = {ker_s}")


# =====================================================================
# ANALYSIS 7: ALTERNATIVE — VERTEX BOX LIFTED TO EDGES VIA d0
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 7: VERTEX BOX LIFTED TO EDGE SPACE")
print("=" * 72)

# d0 Box_v d0^T maps vertex operator to edge space
# But this just gives eigenvalues of Box_v * Delta_0 (product of vertex operators)
# The LIFT: Box_edge_lift = d0 @ Box_v @ d0.T

Box_lift = d0 @ Box_v @ d0.T
evals_lift = np.sort(eigvalsh(Box_lift))
ker_lift = np.sum(np.abs(evals_lift) < 1e-6)
print(f"  Box_lift = d0 @ Box_v @ d0^T")
print(f"  ker(Box_lift) = {ker_lift}")
print(f"  Tr(Box_lift) = {np.trace(Box_lift):.4f}")

# Distinct eigenvalues
ev_lift_distinct = Counter(np.round(evals_lift, 4))
print(f"  Distinct eigenvalues: {len(ev_lift_distinct)}")
for val, mult in sorted(ev_lift_distinct.items()):
    if abs(val) < 50:
        print(f"    {val:12.4f} x {mult}")


# =====================================================================
# ANALYSIS 8: KEY RATIOS AND FRAMEWORK CONSTANTS
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 8: FRAMEWORK CONSTANTS")
print("=" * 72)

# Tr ratios with vertex Box
tr2_v = np.trace(Box_v @ Box_v)
tr3_v = np.trace(Box_v @ Box_v @ Box_v)
tr2_Box1 = np.trace(Box1 @ Box1)
tr3_Box1 = np.trace(Box1 @ Box1 @ Box1)

print(f"  Tr(Box_v^2) = {tr2_v:.4f}")
print(f"  Tr(Box1^2) = {tr2_Box1:.4f}")
print(f"  Ratio Tr(Box1^2)/Tr(Box_v^2) = {tr2_Box1/tr2_v:.6f}")
print(f"  = {tr2_Box1/tr2_v:.6f} = 45 = 9*a1 = 9*5")

print(f"\n  Tr(Box_v^3) = {tr3_v:.4f}")
print(f"  Tr(Box1^3) = {tr3_Box1:.4f}")
if abs(tr3_v) > 0.01:
    print(f"  Ratio Tr(Box1^3)/Tr(Box_v^3) = {tr3_Box1/tr3_v:.6f}")

# Number 13: is it related to any framework constant?
print(f"\n  ker(Box1) = 13")
print(f"  13 = 1 + 12 (trivial + gauge)?")
print(f"  13 = 4 + 9 (vertex ker + ?)?")
print(f"  13 = Nf/Ne * ... ?")
print(f"  a1^2 + 1 = 26, (a1^2+1)/2 = 13  ***")
print(f"  13 is the number of vertices of the icosahedron? No (12).")
print(f"  13 is the rank of SO(26)? Hmm.")
print(f"  sin^2(theta_W) = 6/26 = 6/(2*13). So 2*13 = a1^2+1 = 26.")

# The vertex Box eigenvalue structure:
evals_Bv = np.sort(eigvalsh(Box_v))
ev_Bv_distinct = Counter(np.round(evals_Bv, 4))
print(f"\n  Vertex Box eigenvalues:")
for val, mult in sorted(ev_Bv_distinct.items()):
    print(f"    {val:10.4f} x {mult}")


print("\n" + "=" * 72)
print("EXP100b COMPLETE")
print("=" * 72)

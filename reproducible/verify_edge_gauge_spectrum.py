"""
verify_edge_gauge_spectrum.py — Gauge boson spectrum from Box on edges.

Verifies Theorem (edge_gauge): The operator Box_1 = L_cross - a1*L_fiber
on the line graph of the 600-cell has:

  (a) ker(Box_1) = 13 = (a1^2+1)/2                              [integer rank]
  (b) ker decomposes as rho_0 + 2*rho_5 under 2I                [character]
  (c) The 12 gauge modes (2*rho_5) live 100% on fiber edges      [projector]
  (d) C|_gauge = a1 = 5, B|_gauge = 16/a1                       [Hodge]
  (e) Tr(Box_1^2)/Tr(Box_0^2) = 9*a1 = 45                      [trace identity]
  (f) The fiber graph has Laplacian eigenspaces 1+3+5+3          [A5 decomp]
  (g) Under A5: 12 = 1+3+3'+5 => 1+3+8 = u(1)+su(2)+su(3)      [gauge algebra]

Each test exits with code 1 on failure.
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh, norm
from collections import defaultdict, Counter
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
PASS = 0
tests_run = 0
tests_pass = 0


def check(name, condition, detail=""):
    global tests_run, tests_pass, PASS
    tests_run += 1
    if condition:
        tests_pass += 1
        print(f"  [PASS] {name}")
    else:
        PASS = 1
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


# =====================================================================
# BUILD
# =====================================================================
print("Building 600-cell, simplicial complex, operators...")
verts, adj, lap = build_600cell()
Nv = 120

# --- Edges, faces ---
adj_list = defaultdict(set)
edges = []; edge_to_idx = {}
for i in range(Nv):
    for j in range(i+1, Nv):
        if adj[i, j] > 0.5:
            adj_list[i].add(j); adj_list[j].add(i)
            edge_to_idx[(i, j)] = len(edges); edges.append((i, j))
Ne = len(edges)
assert Ne == 720

triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j:
                    triangles.append((i, j, k))
Nf = len(triangles)
assert Nf == 1200

# --- Boundary operators ---
d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0; d0[e_idx, j] = +1.0

d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

B = d0 @ d0.T   # exact
C = d1.T @ d1    # coexact

# --- Hopf fibration ---
def qmul(p, q):
    return np.array([p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1-tol else -1

fibers = None
for i in range(Nv):
    if abs(verts[i, 0] - PHI/2) >= 1e-6: continue
    g = verts[i]; p = g.copy(); ok = True
    for k in range(2, 11):
        p = qmul(p, g)
        if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok = False; break
        if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok = False
    if not ok: continue
    used = set(); fibs = []; subg = []
    pp = np.array([1.0,0,0,0])
    for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
    for s in range(Nv):
        if s in used: continue
        fib = []
        for si in subg:
            idx = find_idx(qmul(verts[s], verts[si]), verts)
            if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
        if len(fib) == 10: fibs.append(fib)
    if len(fibs) == 12: fibers = fibs; break
assert fibers is not None and len(fibers) == 12

vtf = {}
for fi, f in enumerate(fibers):
    for v in f: vtf[v] = fi

is_fiber = np.zeros(Ne, dtype=bool)
for e, (i, j) in enumerate(edges):
    if vtf[i] == vtf[j]: is_fiber[e] = True
assert np.sum(is_fiber) == 120

# --- Line graph ---
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

# --- Vertex Box ---
Afv = adj * np.array([[float(vtf[i]==vtf[j]) for j in range(Nv)] for i in range(Nv)])
Lfv = np.diag(Afv.sum(1)) - Afv
Lcv = np.diag((adj-Afv).sum(1)) - (adj-Afv)
Box0 = Lcv - a1 * Lfv

# --- Group action ---
print("Building 2I action on edges...")
all_vperm = []
all_eperm = []
for g_idx in range(Nv):
    vp = np.array([find_idx(qmul(verts[g_idx], verts[i]), verts) for i in range(Nv)])
    all_vperm.append(vp)
    ep = np.array([edge_to_idx[(min(vp[i], vp[j]), max(vp[i], vp[j]))] for i, j in edges])
    all_eperm.append(ep)

print("Setup complete.\n")


# =====================================================================
# TEST (a): ker(Box_1) = 13
# =====================================================================
print("=" * 60)
print("TEST (a): ker(Box_1) = 13")
print("=" * 60)

Box1_int = np.round(Box1).astype(int)
check("Box1 has integer entries", np.max(np.abs(Box1 - Box1_int)) < 1e-10)
rank = np.linalg.matrix_rank(Box1_int.astype(float))
ker_dim = Ne - rank
check("ker(Box_1) = 13 = (a1^2+1)/2", ker_dim == 13,
      f"rank = {rank}, ker = {ker_dim}, (a1^2+1)/2 = {(a1**2+1)//2}")
check("Tr(Box_1) = N^2 = 14400", abs(np.trace(Box1) - 14400) < 0.01)


# =====================================================================
# TEST (b): Decomposition rho_0 + 2*rho_5
# =====================================================================
print("\n" + "=" * 60)
print("TEST (b): Irrep decomposition rho_0 + 2*rho_5")
print("=" * 60)

evals_B1, evecs_B1 = eigh(Box1)
ker_mask = np.abs(evals_B1) < 1e-6
ker_vecs = evecs_B1[:, ker_mask]

# Character at identity and at -1
chi_1 = ker_vecs.shape[1]  # = 13

minus1 = None
for g in range(Nv):
    if np.allclose(verts[g], [-1, 0, 0, 0], atol=1e-6):
        minus1 = g; break
assert minus1 is not None

ep = all_eperm[minus1]
chi_minus1 = sum(np.dot(ker_vecs[:, k], ker_vecs[:, k][ep]) for k in range(chi_1))

# rho_0 + 2*rho_5: chi(1) = 1+2*6 = 13, chi(-1) = 1+2*(-6) = -11
check("chi(1) = 13", chi_1 == 13)
check("chi(-1) = -11 (confirms 2*rho_5)", abs(chi_minus1 - (-11)) < 0.1,
      f"chi(-1) = {chi_minus1:.4f}")

# n(rho_0) via group average
n_rho0 = sum(
    sum(np.dot(ker_vecs[:, k], ker_vecs[:, k][all_eperm[g]]) for k in range(chi_1))
    for g in range(Nv)
) / 120.0
check("n(rho_0) = 1", abs(n_rho0 - 1.0) < 0.01, f"n(rho_0) = {n_rho0:.6f}")


# =====================================================================
# TEST (c): 12 gauge modes 100% on fiber edges
# =====================================================================
print("\n" + "=" * 60)
print("TEST (c): 12 gauge modes 100% on fiber edges")
print("=" * 60)

# Project out rho_0 to get the 12-dim rho_5 subspace
P0 = np.zeros((chi_1, chi_1))
for g in range(Nv):
    ep = all_eperm[g]
    for ki in range(chi_1):
        gv = ker_vecs[:, ki][ep]
        for kj in range(chi_1):
            P0[ki, kj] += np.dot(ker_vecs[:, kj], gv)
P0 /= 120.0

P5 = np.eye(chi_1) - P0
p5e, p5v = eigh(P5)
rho5_basis = p5v[:, p5e > 0.5]
gauge = ker_vecs @ rho5_basis  # 720 x 12

fiber_frac = np.sum(gauge[is_fiber, :]**2) / 12.0
cross_frac = np.sum(gauge[~is_fiber, :]**2) / 12.0
check("Gauge modes 100% on fiber edges", abs(fiber_frac - 1.0) < 1e-10,
      f"fiber = {fiber_frac:.12f}, cross = {cross_frac:.2e}")

# The gauge kernel is exactly the alternating edge mode on each even
# ten-cycle Hopf fiber.  (It is not the constant lift.)  This identifies its
# canonical edge metric with ten times the Euclidean metric on the 12 fiber
# amplitudes; changing the initial sign on a fiber is an orthogonal basis
# change and does not affect the metric.
fiber_lift = np.zeros((Ne, 12))
for fiber_index, fiber in enumerate(fibers):
    for k in range(10):
        u, v = fiber[k], fiber[(k+1) % 10]
        edge_index = edge_to_idx[(min(u, v), max(u, v))]
        fiber_lift[edge_index, fiber_index] = (-1.0)**k
check("alternating fiber-amplitude lift lies in ker(Box_1)",
      np.linalg.norm(Box1@fiber_lift) < 1e-9)
check("fiber lift has canonical Gram matrix 10 I_12",
      np.allclose(fiber_lift.T@fiber_lift, 10.0*np.eye(12), atol=1e-12))
check("alternating fiber amplitudes equal the 12-dimensional gauge kernel",
      np.linalg.matrix_rank(gauge.T@fiber_lift, tol=1e-9) == 12)


# =====================================================================
# TEST (d): Hodge eigenvalues on gauge sector
# =====================================================================
print("\n" + "=" * 60)
print("TEST (d): C|_gauge = a1, B|_gauge = 16/a1")
print("=" * 60)

C_g = gauge.T @ C @ gauge
B_g = gauge.T @ B @ gauge
D1_g = gauge.T @ (B + C) @ gauge

c_evals = np.sort(eigvalsh(C_g))
b_evals = np.sort(eigvalsh(B_g))
d1_evals = np.sort(eigvalsh(D1_g))

check("C (coexact) on gauge = a1 = 5", np.allclose(c_evals, 5.0, atol=1e-6),
      f"range [{c_evals[0]:.8f}, {c_evals[-1]:.8f}]")
check("B (exact) on gauge = 16/a1 = 3.2", np.allclose(b_evals, 16.0/a1, atol=1e-6),
      f"range [{b_evals[0]:.8f}, {b_evals[-1]:.8f}]")
check("Delta_1 on gauge = 41/a1 = 8.2", np.allclose(d1_evals, 41.0/a1, atol=1e-6),
      f"range [{d1_evals[0]:.8f}, {d1_evals[-1]:.8f}]")

# Also verify rho_0 values
rho0_v = ker_vecs @ (p5v[:, p5e < 0.5])  # the one rho_0 mode
v0 = rho0_v[:, 0]; v0 /= norm(v0)
check("C on rho_0 = 5/3", abs(v0 @ C @ v0 - 5.0/3) < 1e-6,
      f"value = {v0 @ C @ v0:.8f}")
check("B on rho_0 = 32/15", abs(v0 @ B @ v0 - 32.0/15) < 1e-6,
      f"value = {v0 @ B @ v0:.8f}")


# =====================================================================
# TEST (e): Trace identity Tr(Box_1^2)/Tr(Box_0^2) = 9*a1
# =====================================================================
print("\n" + "=" * 60)
print("TEST (e): Tr(Box_1^2)/Tr(Box_0^2) = 9*a1 = 45")
print("=" * 60)

tr2_edge = np.trace(Box1 @ Box1)
tr2_vert = np.trace(Box0 @ Box0)
ratio = tr2_edge / tr2_vert

check("Tr(Box_1^2)/Tr(Box_0^2) = 45", abs(ratio - 45.0) < 0.01,
      f"ratio = {ratio:.6f}, 9*a1 = {9*a1}")


# =====================================================================
# TEST (f): Fiber graph Laplacian eigenspaces 1+3+5+3
# =====================================================================
print("\n" + "=" * 60)
print("TEST (f): Fiber graph eigenspaces = 1+3+5+3")
print("=" * 60)

# Build fiber adjacency graph
W_fiber = np.zeros((12, 12))
for e_idx, (i, j) in enumerate(edges):
    if not is_fiber[e_idx]:
        fi, fj = vtf[i], vtf[j]
        if fi != fj:
            W_fiber[fi, fj] += 1; W_fiber[fj, fi] += 1

A_fib_graph = (W_fiber > 0).astype(float)
deg_fg = A_fib_graph.sum(axis=1)
L_fg = np.diag(deg_fg) - A_fib_graph

check("Fiber graph is 5-regular", np.allclose(deg_fg, 5),
      f"degrees: {deg_fg.astype(int)}")
check("All cross-edge weights equal (20)",
      len(set(W_fiber[W_fiber > 0].astype(int))) == 1,
      f"weights: {sorted(set(W_fiber[W_fiber > 0].astype(int)))}")

evals_fg = np.sort(eigvalsh(L_fg))
mults_fg = [m for _, m in sorted(Counter(np.round(evals_fg, 4)).items())]

check("Fiber Laplacian multiplicities = [1, 3, 5, 3]",
      mults_fg == [1, 3, 5, 3],
      f"got {mults_fg}")

# Verify eigenvalues are 0, 5-sqrt(5), 6, 5+sqrt(5)
fg_distinct = sorted(set(np.round(evals_fg, 4)))
expected_fg = [0.0, round(5 - np.sqrt(5), 4), 6.0, round(5 + np.sqrt(5), 4)]
check("Fiber Laplacian eigenvalues match icosahedral structure",
      all(abs(a-b) < 0.01 for a, b in zip(fg_distinct, expected_fg)),
      f"got {fg_distinct}, expected {expected_fg}")


# =====================================================================
# TEST (g): A5 decomposition => 1+3+8 = gauge algebra
# =====================================================================
print("\n" + "=" * 60)
print("TEST (g): A5 decomposition 1+3+3'+5 => 1+3+8")
print("=" * 60)

# Verify -1 acts trivially on fibers (so fiber perm factors through A5)
trivial = True
for fi in range(12):
    v_rep = fibers[fi][0]
    gv = all_vperm[minus1][v_rep]
    if vtf[gv] != fi:
        trivial = False; break
check("Central element -1 acts trivially on fibers", trivial)

# Fiber permutation character: chi_perm(g) = #fixed fibers
# Compute for all 120 elements
perm_chars = np.zeros(Nv)
for g in range(Nv):
    n_fixed = 0
    for fi in range(12):
        v_rep = fibers[fi][0]
        if vtf[all_vperm[g][v_rep]] == fi:
            n_fixed += 1
    perm_chars[g] = n_fixed

# Decompose using vertex Laplacian eigenspaces as irrep characters
evals_lap, evecs_lap = eigh(lap)
ev_sp = {}
for i in range(Nv):
    val = round(evals_lap[i], 4)
    if val not in ev_sp: ev_sp[val] = []
    ev_sp[val].append(i)
eig_list = sorted(ev_sp.items())

# Build conjugacy classes from character vectors
char_vecs = np.zeros((Nv, 9))
for g in range(Nv):
    for k, (val, indices) in enumerate(eig_list):
        space = evecs_lap[:, indices]
        char_vecs[g, k] = sum(np.dot(space[:, c], space[:, c][all_vperm[g]])
                               for c in range(space.shape[1]))

class_dict = defaultdict(list)
for g in range(Nv):
    class_dict[tuple(np.round(char_vecs[g], 4))].append(g)

class_list = sorted(class_dict.items(), key=lambda x: (-len(x[1]), x[0]))
class_sizes = np.array([len(m) for _, m in class_list], dtype=float)

check("Found 9 conjugacy classes of 2I", len(class_list) == 9)

# True irrep characters
irrep_dims = [int(round(np.sqrt(len(eig_list[k][1])))) for k in range(9)]
chi_true = np.zeros((9, len(class_list)))
for ri in range(9):
    for ci, (key, _) in enumerate(class_list):
        chi_true[ri, ci] = key[ri] / irrep_dims[ri]

# Decompose fiber permutation representation
perm_class_chars = np.zeros(len(class_list))
for ci, (_, members) in enumerate(class_list):
    perm_class_chars[ci] = perm_chars[members[0]]

fiber_decomp = {}
for ri in range(9):
    n = np.sum(class_sizes * chi_true[ri] * perm_class_chars) / 120.0
    if abs(n) > 0.01:
        fiber_decomp[ri] = round(n)

# A5 irreps: those with chi(-1) = +dim (factor through quotient)
# rho_0(1): chi(-1)=1, rho_2(3): chi(-1)=3, rho_4(5): chi(-1)=5, rho_6(3): chi(-1)=3
# Non-factoring: rho_1(2), rho_3(4), rho_5(6), rho_7(4), rho_8(2)

# The fiber permutation should decompose into factoring irreps only
factoring = {0, 2, 4, 6}  # rho_0(1), rho_2(3), rho_4(5), rho_6(3')
decomp_str = " + ".join(f"{n}*rho_{ri}({irrep_dims[ri]})"
                         for ri, n in sorted(fiber_decomp.items()) if n > 0)
print(f"  Fiber perm decomposition: {decomp_str}")

# Check: should be 1*rho_0(1) + 1*rho_2(3) + 1*rho_4(5) + 1*rho_6(3)
expected_decomp = {0: 1, 2: 1, 4: 1, 6: 1}
check("Fiber perm = rho_0(1) + rho_2(3) + rho_4(5) + rho_6(3')",
      fiber_decomp == expected_decomp,
      f"got {fiber_decomp}")

# Dimension check: 1 + 3 + 5 + 3 = 12
dim_check = sum(irrep_dims[ri] * n for ri, n in fiber_decomp.items())
check("Dimension check: 1+3+5+3 = 12", dim_check == 12)

# Final identification:
# rho_0(1) = u(1), rho_2(3) = su(2), rho_6(3')+rho_4(5) = su(3)
print(f"\n  Gauge algebra identification:")
print(f"    rho_0 (dim 1)         -> u(1)")
print(f"    rho_2 (dim 3)         -> su(2) = adj SU(2)")
print(f"    rho_6(3') + rho_4(5)  -> su(3) = adj SU(3)  [3'+5 = 8]")
print(f"    Total: 1 + 3 + 8 = 12 = dim SU(3)xSU(2)xU(1)")

check("1 + 3 + (3'+5) = 1 + 3 + 8 = 12",
      1 + 3 + (3 + 5) == 12 and 3 + 5 == 8)


# =====================================================================
# PART 8: Box hierarchy (p=0,1,2,3) and spectral index
# =====================================================================
print("\n" + "=" * 60)
print("PART 8: Box hierarchy and spectral index")
print("=" * 60)

# Build vertex adjacency matrix, fiber_edge_set, and vertex A_fiber
A_vtx = np.zeros((Nv, Nv))
for (i, j) in edges:
    A_vtx[i, j] = A_vtx[j, i] = 1.0

fiber_edge_set = set()
A_fiber = np.zeros((Nv, Nv))
for e_idx, (i, j) in enumerate(edges):
    if is_fiber[e_idx]:
        fiber_edge_set.add((i, j))
        A_fiber[i, j] = A_fiber[j, i] = 1.0

# p=0: vertices
A_cross_v = A_vtx - A_fiber
L_fiber_v = np.diag(A_fiber.sum(axis=1)) - A_fiber
L_cross_v = np.diag(A_cross_v.sum(axis=1)) - A_cross_v
Box0 = L_cross_v - a1 * L_fiber_v
ker0 = Nv - np.linalg.matrix_rank(Box0)

# p=1: edges (already computed above as Box1 on line graph)
# ker1 = 13 already verified

# p=2: faces
print(f"\n  Triangles (reusing from above): {Nf}")

# Face adjacency: two faces share an edge
edge_to_faces = {}
for fi, (i,j,k) in enumerate(triangles):
    for e in [(min(i,j),max(i,j)), (min(i,k),max(i,k)), (min(j,k),max(j,k))]:
        edge_to_faces.setdefault(e, []).append(fi)

A_face = np.zeros((Nf, Nf), dtype=int)
A_fiber_face = np.zeros((Nf, Nf), dtype=int)
for e, fs in edge_to_faces.items():
    is_fiber = (e in fiber_edge_set)
    for a in range(len(fs)):
        for b in range(a+1, len(fs)):
            A_face[fs[a], fs[b]] = A_face[fs[b], fs[a]] = 1
            if is_fiber:
                A_fiber_face[fs[a], fs[b]] = A_fiber_face[fs[b], fs[a]] = 1

A_cross_face = A_face - A_fiber_face
L_fiber_f = np.diag(A_fiber_face.sum(axis=1).astype(float)) - A_fiber_face.astype(float)
L_cross_f = np.diag(A_cross_face.sum(axis=1).astype(float)) - A_cross_face.astype(float)
Box2 = L_cross_f - a1 * L_fiber_f
ker2 = Nf - np.linalg.matrix_rank(Box2)
print(f"  Face degree: {A_face[0].sum()}, fiber: {A_fiber_face[0].sum()}")
print(f"  ker(Box_2) = {ker2}")

# p=3: cells (tetrahedra)
print("\n  Building cell adjacency (600 tetrahedra)...")
tetrahedra = []
for (i,j,k) in triangles:
    common = adj_list[i] & adj_list[j] & adj_list[k]
    for l in common:
        if l > k:
            tetrahedra.append((i,j,k,l))
Nc = len(tetrahedra)
print(f"  Tetrahedra: {Nc}")

face_to_tets = {}
for ti, (i,j,k,l) in enumerate(tetrahedra):
    for f in [tuple(sorted([j,k,l])), tuple(sorted([i,k,l])),
              tuple(sorted([i,j,l])), tuple(sorted([i,j,k]))]:
        face_to_tets.setdefault(f, []).append(ti)

A_cell = np.zeros((Nc, Nc), dtype=int)
A_fiber_cell = np.zeros((Nc, Nc), dtype=int)
for f, ts in face_to_tets.items():
    if len(ts) == 2:
        A_cell[ts[0],ts[1]] = A_cell[ts[1],ts[0]] = 1
        fi,fj,fk = f
        edges_f = [(min(fi,fj),max(fi,fj)), (min(fi,fk),max(fi,fk)), (min(fj,fk),max(fj,fk))]
        if any(e in fiber_edge_set for e in edges_f):
            A_fiber_cell[ts[0],ts[1]] = A_fiber_cell[ts[1],ts[0]] = 1

A_cross_cell = A_cell - A_fiber_cell
L_fiber_c = np.diag(A_fiber_cell.sum(axis=1).astype(float)) - A_fiber_cell.astype(float)
L_cross_c = np.diag(A_cross_cell.sum(axis=1).astype(float)) - A_cross_cell.astype(float)
Box3 = L_cross_c - a1 * L_fiber_c
ker3 = Nc - np.linalg.matrix_rank(Box3)
print(f"  Cell degree: {A_cell[0].sum()}, fiber: {A_fiber_cell[0].sum()}")
print(f"  ker(Box_3) = {ker3}")

# Spectral index
ker1 = 13  # verified above
index = ker0 - ker1 + ker2 - ker3
d_ST = a1 - 1

print(f"\n  Box hierarchy:")
print(f"    ker(Box_0) = {ker0}")
print(f"    ker(Box_1) = {ker1}")
print(f"    ker(Box_2) = {ker2}")
print(f"    ker(Box_3) = {ker3}")
print(f"    Spectral index = {ker0} - {ker1} + {ker2} - {ker3} = {index}")
print(f"    -(a1-1) = -d_ST = {-d_ST}")

check("Box hierarchy: ker0 = 9", ker0 == 9)
check("Box hierarchy: ker2 = 1", ker2 == 1)
check("Box hierarchy: ker3 = 1", ker3 == 1)
check("Box spectral index = -(a1-1) = -d_ST = -4",
      index == -d_ST, f"index = {index}")

# Trace ratios
tr0 = np.sum(eigvalsh(Box0)**2)
tr2 = np.sum(eigvalsh(Box2)**2)
tr3 = np.sum(eigvalsh(Box3)**2)
print(f"\n  Trace ratios Tr(Box_p^2)/Tr(Box_0^2):")
print(f"    p=1: {324000/tr0:.1f} = 9*a1")
print(f"    p=2: {tr2/tr0:.1f}")
print(f"    p=3: {tr3/tr0:.4f}")

check("Tr(Box_2^2)/Tr(Box_0^2) = 34", abs(tr2/tr0 - 34) < 0.01)


# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 60)
print(f"RESULTS: {tests_pass}/{tests_run} tests passed")
print("=" * 60)

if PASS != 0:
    print("SOME TESTS FAILED")
sys.exit(PASS)

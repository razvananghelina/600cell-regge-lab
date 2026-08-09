"""
exp105: Two directions from Box.

PART A: Eigenvalue(Box_0) -> fermion mass? Direct algebraic formula?
  9 irreps, 9 eigenvalues, 9 fermion masses. Is there a map?

PART B: Box on FACES (1200-dimensional).
  Vertices -> scalars (c^2, photon)
  Edges -> vectors (gauge bosons)
  Faces -> rank-2 tensors (graviton?)
  Predict: ker(Box_2) encodes graviton content.
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
N = 120

# =====================================================================
# BUILD
# =====================================================================
print("Building 600-cell + full simplicial complex...")
verts, adj, lap = build_600cell()
Nv = 120

adj_list = defaultdict(set)
edges = []; edge_to_idx = {}
for i in range(Nv):
    for j in range(i+1, Nv):
        if adj[i, j] > 0.5:
            adj_list[i].add(j); adj_list[j].add(i)
            edge_to_idx[(i, j)] = len(edges); edges.append((i, j))
Ne = 720

triangles = []
tri_to_idx = {}
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j:
                    tri_to_idx[(i,j,k)] = len(triangles)
                    triangles.append((i, j, k))
Nf = len(triangles)

# Build tetrahedra (cells) for d2
print("  Finding 600 tetrahedra...")
tetrahedra = []
for i in range(Nv):
    for j in adj_list[i]:
        if j <= i: continue
        for k in adj_list[i] & adj_list[j]:
            if k <= j: continue
            for l in adj_list[i] & adj_list[j] & adj_list[k]:
                if l <= k:
                    continue
                tetrahedra.append((i, j, k, l))
Nc = len(tetrahedra)
print(f"  Tetrahedra: {Nc} (expect 600)")
assert Nc == 600

# Boundary operators
d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0; d0[e_idx, j] = +1.0

d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

# d2: faces -> tetrahedra (Nc x Nf)
d2 = np.zeros((Nc, Nf))
for c_idx, (i, j, k, l) in enumerate(tetrahedra):
    # Boundary of tetrahedron [i,j,k,l] = [j,k,l] - [i,k,l] + [i,j,l] - [i,j,k]
    faces_of_tet = [(j,k,l,+1), (i,k,l,-1), (i,j,l,+1), (i,j,k,-1)]
    for (a, b, c, sgn) in faces_of_tet:
        tri_key = tuple(sorted([a, b, c]))
        if tri_key in tri_to_idx:
            d2[c_idx, tri_to_idx[tri_key]] = sgn

# Verify d2 @ d1 = 0
err_d2d1 = norm(d2 @ d1)
print(f"  ||d2 @ d1|| = {err_d2d1:.2e} (expect 0)")
assert err_d2d1 < 1e-10

# Hopf fibration
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

vtf = {}
for fi, f in enumerate(fibers):
    for v in f: vtf[v] = fi

is_fiber_edge = np.zeros(Ne, dtype=bool)
for e, (i, j) in enumerate(edges):
    if vtf[i] == vtf[j]: is_fiber_edge[e] = True

print("Done.\n")


# =====================================================================
# PART A: EIGENVALUE(Box_0) -> FERMION MASS
# =====================================================================
print("=" * 72)
print("PART A: BOX EIGENVALUE vs FERMION MASS")
print("=" * 72)

# Vertex Box
Afv = adj * np.array([[float(vtf[i]==vtf[j]) for j in range(Nv)] for i in range(Nv)])
Lfv = np.diag(Afv.sum(1)) - Afv
Lcv = np.diag((adj-Afv).sum(1)) - (adj-Afv)
Box0 = Lcv - a1 * Lfv

evals0, evecs0 = eigh(Box0)

# Get eigenvalue per irrep (use vertex Laplacian eigenspaces)
evals_lap, evecs_lap = eigh(lap)
ev_sp = {}
for i in range(Nv):
    val = round(evals_lap[i], 4)
    if val not in ev_sp: ev_sp[val] = []
    ev_sp[val].append(i)
eig_list = sorted(ev_sp.items())
irrep_dims = [int(round(np.sqrt(len(eig_list[k][1])))) for k in range(9)]

# Box eigenvalue per irrep
box_per_irrep = []
for k, (lap_val, indices) in enumerate(eig_list):
    space = evecs_lap[:, indices]
    # Project Box onto this eigenspace
    box_restricted = space.T @ Box0 @ space
    box_eval = np.mean(np.diag(box_restricted))  # should be scalar on each irrep
    box_per_irrep.append(box_eval)

# McKay node labels and fermion assignments
# From the paper: irreps rho_0..rho_8 with dims 1,2,3,4,5,6,3,4,2
# Fermion mass exponents n_f (from paper Section 7):
# e(0), mu(11), tau(17), u(6), c(15), t(25), d(8), s(11), b(20)
# McKay assignment: each fermion to a node
# From memory: McKay mass-distance confirmed with r=0.93

# The 9 irreps correspond to McKay nodes 0-8
# Node assignment: rho_k -> fermion with mass exponent n_k
# From the paper (Thm uniqueness):
mckay_fermion = {
    0: ("e", 0),      # rho_0 dim 1
    1: ("u", 6),      # rho_1 dim 2
    2: ("d", 8),      # rho_2 dim 3
    3: ("s/mu", 11),  # rho_3 dim 4 (degenerate)
    4: ("c", 15),     # rho_4 dim 5
    5: ("tau", 17),   # rho_5 dim 6
    6: ("b", 20),     # rho_6 dim 3'
    7: ("t", 25),     # rho_7 dim 4'
    8: ("mu/s", 11),  # rho_8 dim 2' (degenerate with rho_3)
}

print(f"\n  {'irrep':>6s} {'dim':>4s} {'Box_eval':>10s} {'|Box|':>8s} {'fermion':>8s} {'n_f':>5s} {'m/m_e':>10s}")
print(f"  {'-'*60}")

box_abs = []
mass_exp = []
for k in range(9):
    dim = irrep_dims[k]
    bev = box_per_irrep[k]
    fermion, nf = mckay_fermion[k]
    mass_ratio = PHI**nf
    box_abs.append(abs(bev))
    mass_exp.append(nf)
    print(f"  rho_{k:d}  {dim:4d} {bev:10.4f} {abs(bev):8.4f} {fermion:>8s} {nf:5d} {mass_ratio:10.2f}")

# Correlation
from numpy import corrcoef
r_box_mass = corrcoef(box_abs, mass_exp)[0, 1]
print(f"\n  Correlation |Box eigenvalue| vs mass exponent: r = {r_box_mass:.4f}")

# Is it monotonic?
sorted_by_box = sorted(zip(box_abs, mass_exp))
monotonic = all(sorted_by_box[i][1] <= sorted_by_box[i+1][1] for i in range(len(sorted_by_box)-1))
print(f"  Monotonic (|Box| -> n_f)? {monotonic}")

# Check for algebraic formula: n_f = f(|Box|)?
# Try linear: n_f = a * |Box| + b
A_fit = np.column_stack([box_abs, np.ones(9)])
coeffs, residuals, _, _ = np.linalg.lstsq(A_fit, mass_exp, rcond=None)
print(f"\n  Linear fit: n_f = {coeffs[0]:.4f} * |Box| + {coeffs[1]:.4f}")
predicted = A_fit @ coeffs
for k in range(9):
    fermion, nf = mckay_fermion[k]
    print(f"    rho_{k}: predicted = {predicted[k]:.2f}, actual = {nf}, diff = {predicted[k]-nf:+.2f}")


# =====================================================================
# PART B: BOX ON FACES (1200-dimensional)
# =====================================================================
print("\n" + "=" * 72)
print("PART B: BOX ON FACES (1200-dimensional)")
print("=" * 72)

# Hodge Laplacian on faces: Delta_2 = d1 d1^T + d2^T d2
print("  Building face Hodge Laplacian (1200 x 1200)...")
B2 = d1 @ d1.T       # exact part of face Laplacian (from edges)
C2 = d2.T @ d2       # coexact part (from tetrahedra)
Delta_2 = B2 + C2

# Verify B2*C2 = 0
err_B2C2 = norm(B2 @ C2, 'fro')
print(f"  ||B2 @ C2||_F = {err_B2C2:.2e} (expect 0)")

# Dimensions
evals_B2 = eigvalsh(B2)
evals_C2 = eigvalsh(C2)
dim_exact2 = np.sum(evals_B2 > 1e-8)
dim_coexact2 = np.sum(evals_C2 > 1e-8)
dim_harmonic2 = Nf - dim_exact2 - dim_coexact2
print(f"  Exact dim (rank B2): {dim_exact2} (expect 601 = Ne - Nv + 1)")
print(f"  Coexact dim (rank C2): {dim_coexact2} (expect 599 = Nc - 1)")
print(f"  Harmonic dim (b2): {dim_harmonic2} (expect 0 for S^3)")

# Classify faces as fiber vs cross
# A face has 0, 1, 2, or 3 fiber edges
n_fiber_per_face = np.zeros(Nf, dtype=int)
for f_idx, (i, j, k) in enumerate(triangles):
    nf = int(vtf[i]==vtf[j]) + int(vtf[j]==vtf[k]) + int(vtf[i]==vtf[k])
    n_fiber_per_face[f_idx] = nf

face_types = Counter(n_fiber_per_face)
print(f"\n  Face classification:")
for nf in sorted(face_types):
    print(f"    {nf} fiber edges: {face_types[nf]} faces")

# Face "fiber" flag: a face is "fiber-type" if it has exactly 1 fiber edge
# (the only type that connects fiber to cross)
is_fiber_face = (n_fiber_per_face == 1)
is_cross_face = (n_fiber_per_face == 0)
print(f"  Fiber-type faces (1 fiber edge): {np.sum(is_fiber_face)}")
print(f"  Cross-type faces (0 fiber edges): {np.sum(is_cross_face)}")

# Build face adjacency graph: two faces are adjacent if they share an edge
print("\n  Building face adjacency (1200 x 1200)...")
# For each edge, find all faces containing it
edge_to_faces = defaultdict(list)
for f_idx, (i, j, k) in enumerate(triangles):
    for e_key in [(i,j), (j,k), (i,k)]:
        e_idx = edge_to_idx[e_key]
        edge_to_faces[e_idx].append(f_idx)

A_face = np.zeros((Nf, Nf))
for e_idx in range(Ne):
    faces = edge_to_faces[e_idx]
    for a in range(len(faces)):
        for b in range(a+1, len(faces)):
            A_face[faces[a], faces[b]] = 1.0
            A_face[faces[b], faces[a]] = 1.0

face_degrees = A_face.sum(axis=1)
print(f"  Face graph degrees: min={face_degrees.min():.0f}, max={face_degrees.max():.0f}")
print(f"  (Each face shares edges with other faces)")

# Fiber adjacency on faces: two faces share a FIBER edge
A_fiber_face = np.zeros((Nf, Nf))
for e_idx in range(Ne):
    if is_fiber_edge[e_idx]:
        faces = edge_to_faces[e_idx]
        for a in range(len(faces)):
            for b in range(a+1, len(faces)):
                A_fiber_face[faces[a], faces[b]] = 1.0
                A_fiber_face[faces[b], faces[a]] = 1.0

fiber_face_deg = A_fiber_face.sum(axis=1)
unique_ffd = Counter(fiber_face_deg.astype(int))
print(f"  Fiber face adjacency degrees: {dict(unique_ffd)}")

# Box on faces: analog of vertex/edge
# Box_2 = L_cross_face - a1 * L_fiber_face
A_cross_face = A_face - A_fiber_face
L_fiber_face = np.diag(A_fiber_face.sum(1)) - A_fiber_face
L_cross_face = np.diag(A_cross_face.sum(1)) - A_cross_face
Box2 = L_cross_face - a1 * L_fiber_face

print(f"\n  Box_2 = L_cross - a1*L_fiber on face graph")
print(f"  Box_2 is {Nf}x{Nf}, integer entries: {np.allclose(Box2, np.round(Box2))}")

evals2 = np.sort(eigvalsh(Box2))
ker2 = np.sum(np.abs(evals2) < 1e-4)
tr2 = np.trace(Box2)
tr2_sq = np.trace(Box2 @ Box2)

print(f"  ker(Box_2) = {ker2}")
print(f"  Tr(Box_2) = {tr2:.0f}")
print(f"  Tr(Box_2^2) = {tr2_sq:.0f}")
print(f"  min eigenvalue: {evals2[0]:.4f}")
print(f"  max eigenvalue: {evals2[-1]:.4f}")

# Distinct eigenvalues
ev2_distinct = Counter(np.round(evals2, 4))
print(f"  Distinct eigenvalues: {len(ev2_distinct)}")
print(f"\n  {'eigenvalue':>10s} {'mult':>5s}")
print(f"  {'-'*20}")
for val, mult in sorted(ev2_distinct.items()):
    print(f"  {val:10.4f} {mult:5d}")

# Trace ratios
tr2_v = np.trace(Box0 @ Box0)
tr2_e_sq = 324000  # from exp100
print(f"\n  Trace ratios:")
print(f"    Tr(Box_2^2)/Tr(Box_0^2) = {tr2_sq/tr2_v:.4f}")
print(f"    Tr(Box_2^2)/Tr(Box_1^2) = {tr2_sq/tr2_e_sq:.4f}")
print(f"    Tr(Box_2) = {tr2:.0f}")

# Check: is Tr(Box_2) = N^2 or something framework?
for name, val in [("N^2", N**2), ("Nf*a1", Nf*a1), ("Nf*b1", Nf*b1),
                   ("Ne*a1", Ne*a1), ("Ne*b1", Ne*b1), ("Nc*N", Nc*N)]:
    if abs(tr2 - val) < 1:
        print(f"    Tr(Box_2) = {name} = {val}")


# =====================================================================
# PART B2: KERNEL ANALYSIS OF Box_2
# =====================================================================
print("\n" + "=" * 72)
print("PART B2: KERNEL OF Box_2")
print("=" * 72)

if ker2 > 0:
    print(f"  ker(Box_2) = {ker2}")

    # Framework constant identification
    for name, val in [("1", 1), ("a1-1", a1-1), ("a1", a1), ("b1", b1),
                       ("N_eig", 9), ("13=(a1^2+1)/2", 13), ("degree", 12),
                       ("Nc/Ne", Nc/Ne), ("Nf/Ne", Nf/Ne)]:
        if abs(ker2 - val) < 0.5:
            print(f"  *** ker = {name} = {int(val)} ***")

    # Decompose under 2I
    print("\n  Computing 2I action on kernel...")
    evals_B2_full, evecs_B2_full = eigh(Box2)
    ker_mask = np.abs(evals_B2_full) < 1e-4
    ker_vecs2 = evecs_B2_full[:, ker_mask]

    # Build face permutations under 2I
    all_vperm = []
    for g_idx in range(Nv):
        vp = np.array([find_idx(qmul(verts[g_idx], verts[i]), verts) for i in range(Nv)])
        all_vperm.append(vp)

    # Face permutation: g maps triangle (i,j,k) to (g(i), g(j), g(k))
    all_fperm = []
    for g_idx in range(Nv):
        fp = np.zeros(Nf, dtype=int)
        vp = all_vperm[g_idx]
        for f_idx, (i, j, k) in enumerate(triangles):
            gi, gj, gk = vp[i], vp[j], vp[k]
            tri_key = tuple(sorted([gi, gj, gk]))
            fp[f_idx] = tri_to_idx[tri_key]
        all_fperm.append(fp)

    # Character of kernel on identity and -1
    chi_1 = ker2
    minus1 = None
    for g in range(Nv):
        if np.allclose(verts[g], [-1, 0, 0, 0], atol=1e-6):
            minus1 = g; break

    fp = all_fperm[minus1]
    chi_minus1 = sum(np.dot(ker_vecs2[:, k], ker_vecs2[:, k][fp]) for k in range(ker2))
    print(f"  chi(1) = {chi_1}, chi(-1) = {chi_minus1:.1f}")

    # Full decomposition using conjugacy classes
    evals_lap_v, evecs_lap_v = eigh(lap)
    ev_sp_v = {}
    for i in range(Nv):
        val = round(evals_lap_v[i], 4)
        if val not in ev_sp_v: ev_sp_v[val] = []
        ev_sp_v[val].append(i)
    eig_list_v = sorted(ev_sp_v.items())
    irrep_dims_v = [int(round(np.sqrt(len(eig_list_v[k][1])))) for k in range(9)]

    char_vecs_v = np.zeros((Nv, 9))
    for g in range(Nv):
        for k, (val, indices) in enumerate(eig_list_v):
            space = evecs_lap_v[:, indices]
            char_vecs_v[g, k] = sum(np.dot(space[:, c], space[:, c][all_vperm[g]])
                                     for c in range(space.shape[1]))

    class_dict_v = defaultdict(list)
    for g in range(Nv):
        class_dict_v[tuple(np.round(char_vecs_v[g], 4))].append(g)
    class_list_v = sorted(class_dict_v.items(), key=lambda x: (-len(x[1]), x[0]))
    class_sizes_v = np.array([len(m) for _, m in class_list_v], dtype=float)

    chi_true_v = np.zeros((9, len(class_list_v)))
    for ri in range(9):
        for ci, (key, _) in enumerate(class_list_v):
            chi_true_v[ri, ci] = key[ri] / irrep_dims_v[ri]

    # Kernel character on each class
    ker2_chars = np.zeros(len(class_list_v))
    for ci, (_, members) in enumerate(class_list_v):
        g = members[0]
        fp = all_fperm[g]
        ker2_chars[ci] = sum(np.dot(ker_vecs2[:, k], ker_vecs2[:, k][fp]) for k in range(ker2))

    print(f"\n  Irrep decomposition of ker(Box_2):")
    total_dim = 0
    for ri in range(9):
        n = np.sum(class_sizes_v * chi_true_v[ri] * ker2_chars) / 120.0
        if abs(n) > 0.01:
            print(f"    rho_{ri} (dim {irrep_dims_v[ri]}): n = {n:.4f}")
            total_dim += irrep_dims_v[ri] * round(n)
    print(f"  Total dim: {total_dim} (expect {ker2})")

else:
    print(f"  ker(Box_2) = 0. No massless face modes.")
    print(f"  Trying to find the structure in the near-zero modes...")
    near_zero = evals2[np.abs(evals2) < 2]
    print(f"  Eigenvalues |lambda| < 2: {len(near_zero)}")


# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

print(f"""
  PART A: Box eigenvalue vs fermion mass
  - Correlation |Box| vs n_f: r = {r_box_mass:.4f}
  - Monotonic: {monotonic}
  - (Already known from McKay mass-distance, r = 0.93)

  PART B: Box on faces (1200-dim)
  - ker(Box_2) = {ker2}
  - Tr(Box_2) = {tr2:.0f}
  - {len(ev2_distinct)} distinct eigenvalues

  THE BOX HIERARCHY:
    Box_0 (vertices, 120-dim): ker = 9 = rho_0 + 2*rho_1 + 2*rho_8
    Box_1 (edges, 720-dim):    ker = 13 = rho_0 + 2*rho_5
    Box_2 (faces, 1200-dim):   ker = {ker2}
""")

print("=" * 72)
print("EXP105 COMPLETE")
print("=" * 72)

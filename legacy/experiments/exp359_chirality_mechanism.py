"""
exp359_chirality_mechanism.py
==============================
OPEN-2: Can the chiral SU(2)_L coupling be derived?

The question: why does SU(2) couple ONLY to gamma = -1 states?

Strategy:
1. Analytical: Show Tr(gamma * D^{2k}) = 0 on product space
   (consequence of spectral balance + Euler(S^3) = 0)
2. Hopf T3 quadrupole: T3 = w^2 + z^2 - x^2 - y^2 on quaternions
   Check if T3 selects the right chiral structure
3. Spectral action on H+ vs H-: check if gauge terms differ
4. Check [T3, D_M], [T3, gamma_form], [T3, A] commutation relations

Dependencies: numpy, scipy
"""

import numpy as np
from scipy.linalg import eigh
from itertools import permutations, product as cartesian_product

a1 = 5
b1 = a1 + 1
phi = (1.0 + np.sqrt(a1)) / 2.0
N_gen = 3
TOL = 1e-8

print("=" * 72)
print("EXP-359: CHIRALITY MECHANISM (OPEN-2)")
print("=" * 72)

# ============================================================================
# BUILD 600-CELL
# ============================================================================
print("\n--- Building 600-cell ---")

def build_2I():
    verts = set()
    def add_vert(v):
        arr = np.array(v, dtype=float)
        n = np.linalg.norm(arr)
        if n > 1e-12:
            arr = arr / n
        verts.add(tuple(np.round(arr, 10)))
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            add_vert(v)
    for signs in cartesian_product([0.5, -0.5], repeat=4):
        add_vert(list(signs))
    base = [0.0, 0.5, phi / 2.0, 1.0 / (2.0 * phi)]
    even_perms = [p for p in permutations(range(4))
                  if sum(1 for i in range(4) for j in range(i+1, 4)
                         if p[i] > p[j]) % 2 == 0]
    for perm in even_perms:
        coords = [base[perm[i]] for i in range(4)]
        nz = [i for i in range(4) if abs(coords[i]) > 1e-12]
        for signs in cartesian_product([1, -1], repeat=len(nz)):
            v = list(coords)
            for idx, s in zip(nz, signs):
                v[idx] *= s
            add_vert(v)
    return np.array(sorted(verts))

verts = build_2I()
N = len(verts)
print(f"  Vertices: {N}")

# Adjacency
dot_threshold = phi / 2.0
dots = verts @ verts.T
np.clip(dots, -1.0, 1.0, out=dots)
A = np.zeros((N, N), dtype=float)
for i in range(N):
    for j in range(i + 1, N):
        if abs(dots[i, j] - dot_threshold) < 1e-6:
            A[i, j] = 1.0
            A[j, i] = 1.0

degree = 12

# ============================================================================
# PART 1: HOPF T3 QUADRUPOLE
# ============================================================================
print("\n" + "=" * 72)
print("PART 1: HOPF T3 QUADRUPOLE")
print("=" * 72)

# T3 = w^2 + z^2 - x^2 - y^2 for quaternion (w,x,y,z)
# This is one of the three Hopf coordinates on S3
# It gives a function T3: S3 -> [-1, 1]

T3_vals = np.array([v[0]**2 + v[3]**2 - v[1]**2 - v[2]**2 for v in verts])

# Check unique values and their multiplicities
T3_rounded = np.round(T3_vals, 8)
unique_T3, counts_T3 = np.unique(T3_rounded, return_counts=True)

print(f"\n  T3 values and multiplicities:")
for val, cnt in zip(unique_T3, counts_T3):
    print(f"    T3 = {val:+8.5f}, mult = {cnt}")

# Check the 48+48 split
n_pos = np.sum(T3_vals > TOL)
n_zero = np.sum(np.abs(T3_vals) <= TOL)
n_neg = np.sum(T3_vals < -TOL)
print(f"\n  T3 > 0: {n_pos} vertices")
print(f"  T3 = 0: {n_zero} vertices")
print(f"  T3 < 0: {n_neg} vertices")

# Does T3 give 48+48+24 or something like that?
# 16 per generation x 3 generations = 48 per sign = 96 + 24 neutral?

# T3 as a diagonal matrix on vertices
T3_diag = np.diag(T3_vals)

# Check [T3, A] - does T3 commute with the adjacency?
commutator_T3_A = T3_diag @ A - A @ T3_diag
comm_norm = np.linalg.norm(commutator_T3_A, 'fro')
print(f"\n  ||[T3, A]||_F = {comm_norm:.6f}")
if comm_norm < TOL:
    print("  T3 COMMUTES with A -> conserved quantity")
else:
    print("  T3 does NOT commute with A -> not conserved")
    # But check if it approximately commutes
    print(f"  ||[T3, A]||_F / ||A||_F = {comm_norm / np.linalg.norm(A, 'fro'):.6f}")

# ============================================================================
# PART 2: ALL THREE HOPF COORDINATES
# ============================================================================
print("\n" + "=" * 72)
print("PART 2: ALL THREE HOPF COORDINATES")
print("=" * 72)

# The three Hopf coordinates on S3 are:
# H1 = 2(wy + xz)     [maps to x on S2]
# H2 = 2(xy - wz)     [maps to y on S2]
# H3 = w^2+z^2-x^2-y^2 [maps to z on S2]
# These satisfy H1^2 + H2^2 + H3^2 = 1

H1 = np.array([2*(v[0]*v[2] + v[1]*v[3]) for v in verts])
H2 = np.array([2*(v[1]*v[2] - v[0]*v[3]) for v in verts])
H3 = T3_vals.copy()

# Verify S2 constraint
s2_check = H1**2 + H2**2 + H3**2
print(f"  H1^2 + H2^2 + H3^2: max deviation from 1 = {np.max(np.abs(s2_check - 1)):.2e}")

# Build all three Hopf diagonal operators
H1_diag = np.diag(H1)
H2_diag = np.diag(H2)
H3_diag = np.diag(H3)

# Check commutators with A
for name, H_diag in [("H1", H1_diag), ("H2", H2_diag), ("H3", H3_diag)]:
    comm = H_diag @ A - A @ H_diag
    cn = np.linalg.norm(comm, 'fro')
    print(f"  ||[{name}, A]||_F = {cn:.6f}")

# Check commutators between Hopf coordinates
for n1, H1d in [("H1", H1_diag), ("H2", H2_diag)]:
    for n2, H2d in [("H2", H2_diag), ("H3", H3_diag)]:
        if n1 >= n2:
            continue
        comm = H1d @ H2d - H2d @ H1d
        cn = np.linalg.norm(comm, 'fro')
        print(f"  ||[{n1}, {n2}]||_F = {cn:.6f}")

# ============================================================================
# PART 3: T3 SPECTRUM ON ADJACENCY EIGENSPACES
# ============================================================================
print("\n" + "=" * 72)
print("PART 3: T3 STRUCTURE ON ADJACENCY EIGENSPACES")
print("=" * 72)

# Get full A eigensystem
eigvals_A, eigvecs_A = eigh(A)

# Group by eigenvalue
tol_group = 0.1
eigenspaces = {}
for i in range(N):
    ev = round(eigvals_A[i], 4)
    found = False
    for key in eigenspaces:
        if abs(ev - key) < tol_group:
            eigenspaces[key].append(i)
            found = True
            break
    if not found:
        eigenspaces[ev] = [i]

# For each eigenspace, compute T3 restricted to that space
print(f"\n  T3 expectation values on adjacency eigenspaces:")
for ev in sorted(eigenspaces.keys(), reverse=True):
    indices = eigenspaces[ev]
    mult = len(indices)
    # Get eigenvectors for this eigenspace
    V_space = eigvecs_A[:, indices]  # (120, mult)
    # Project T3 onto this eigenspace
    T3_proj = V_space.T @ T3_diag @ V_space  # (mult, mult)
    # Eigenvalues of T3 in this subspace
    t3_evals = np.sort(np.linalg.eigvalsh(T3_proj))
    t3_evals_r = np.round(t3_evals, 6)
    unique_t3, cnt_t3 = np.unique(t3_evals_r, return_counts=True)
    t3_str = ", ".join(f"{v:+.4f}({c})" for v, c in zip(unique_t3, cnt_t3))
    print(f"    lambda_A = {ev:+8.4f} (mult {mult:2d}): T3 evals = {t3_str}")

# ============================================================================
# PART 4: ANALYTICAL RESULT - SPECTRAL SUPERTRACE VANISHES
# ============================================================================
print("\n" + "=" * 72)
print("PART 4: ANALYTICAL - Tr(gamma * D^{2k}) ON PRODUCT SPACE")
print("=" * 72)

print("""
  THEOREM (spectral supertrace vanishing):
  On the product space H = H_M(2640) x H_F(30) with
    D = D_M x 1 + gamma_form x D_F
    gamma = gamma_form x gamma_F

  Tr(gamma * D^{2k}) = 0 for all k >= 0.

  PROOF:
  Since {D_M, gamma_form} = 0 (simplicial Dirac anticommutes with
  form-degree parity), we have:
    D^2 = D_M^2 x 1 + 1 x D_F^2   (cross terms cancel)

  Therefore D^{2k} = sum_{j=0}^{k} C(k,j) D_M^{2j} x D_F^{2(k-j)}

  Tr(gamma * D^{2k}) = sum_j C(k,j) Tr(gamma_form * D_M^{2j}) * Tr(gamma_F * D_F^{2(k-j)})
                      + sum_j ... terms with Tr(gamma_form) or Tr(gamma_F * D_F^{...})

  Actually more carefully:
  Tr(gamma * D^{2k}) = Tr((gamma_form x gamma_F) * sum_j C(k,j) D_M^{2j} x D_F^{2(k-j)})
                      = sum_j C(k,j) Tr(gamma_form * D_M^{2j}) * Tr(gamma_F * D_F^{2(k-j)})

  Key facts:
  1. Tr(gamma_form * D_M^{2j}) = Str_form(D_M^{2j}) = 0 (spectral balance, Prop. 6)
  2. This holds for ALL j >= 0 (including j=0: Tr(gamma_form) = chi(S3) = 0)

  Therefore EVERY term in the sum vanishes. QED.

  CONSEQUENCE: The spectral action Tr(f(D/Lambda)) is IDENTICAL on H+ and H-.
  The chiral gauge coupling CANNOT come from the spectral action on S3.
  It requires either:
    (a) 4D continuum limit (where chi != 0)
    (b) A non-spectral-action mechanism
    (c) Inner fluctuations with algebra structure
""")

# Numerical verification: compute Str on simplicial D
print("  Numerical verification on D_M (2640x2640):")

# Build simplicial complex and D_M
edges = []
for i in range(N):
    for j in range(i + 1, N):
        if A[i, j] > 0.5:
            edges.append((i, j))
E = len(edges)
edge_idx = {e: k for k, e in enumerate(edges)}

faces = []
for k, (i, j) in enumerate(edges):
    for l in range(j + 1, N):
        if A[i, l] > 0.5 and A[j, l] > 0.5:
            faces.append((i, j, l))
F = len(faces)

cells = []
for fi, (i, j, k) in enumerate(faces):
    for l in range(k + 1, N):
        if A[i, l] > 0.5 and A[j, l] > 0.5 and A[k, l] > 0.5:
            cells.append((i, j, k, l))
C_cells = len(cells)

dim_H = N + E + F + C_cells

# Build gamma_form = (-1)^p as diagonal
gamma_form = np.zeros(dim_H)
gamma_form[:N] = 1.0          # 0-forms: +1
gamma_form[N:N+E] = -1.0      # 1-forms: -1
gamma_form[N+E:N+E+F] = 1.0   # 2-forms: +1
gamma_form[N+E+F:] = -1.0     # 3-forms: -1

# Tr(gamma_form) = chi(S3)
tr_gamma = np.sum(gamma_form)
print(f"    Tr(gamma_form) = {tr_gamma:.0f} = chi(S3)")

# Build D_M (need boundary operators)
d0 = np.zeros((E, N), dtype=float)
for k, (i, j) in enumerate(edges):
    d0[k, i] = -1.0
    d0[k, j] = 1.0

d1 = np.zeros((F, E), dtype=float)
for fi, (i, j, k) in enumerate(faces):
    e_ij = edge_idx.get((i, j), edge_idx.get((j, i), -1))
    e_ik = edge_idx.get((i, k), edge_idx.get((k, i), -1))
    e_jk = edge_idx.get((j, k), edge_idx.get((k, j), -1))
    d1[fi, e_ij] = 1.0
    d1[fi, e_ik] = -1.0
    d1[fi, e_jk] = 1.0

face_idx = {f: fi for fi, f in enumerate(faces)}
d2 = np.zeros((C_cells, F), dtype=float)
for ci, (i, j, k, l) in enumerate(cells):
    f_faces = [(j, k, l), (i, k, l), (i, j, l), (i, j, k)]
    f_signs = [1.0, -1.0, 1.0, -1.0]
    for face, sign in zip(f_faces, f_signs):
        fi = face_idx.get(face, -1)
        if fi >= 0:
            d2[ci, fi] = sign

D_M = np.zeros((dim_H, dim_H), dtype=float)
D_M[N:N+E, 0:N] = d0
D_M[0:N, N:N+E] = d0.T
D_M[N+E:N+E+F, N:N+E] = d1
D_M[N:N+E, N+E:N+E+F] = d1.T
D_M[N+E+F:dim_H, N+E:N+E+F] = d2
D_M[N+E:N+E+F, N+E+F:dim_H] = d2.T

# Compute Str(D_M^{2k}) for k=1,...,4
print("    Spectral supertraces Str(D_M^{2k}):")
D2 = D_M @ D_M
power = np.eye(dim_H)
for k in range(1, 5):
    power = power @ D2
    str_k = np.sum(gamma_form * np.diag(power))
    print(f"      Str(D_M^{2*k}) = {str_k:.6e}")

# ============================================================================
# PART 5: INNER CHIRALITY gamma_F
# ============================================================================
print("\n" + "=" * 72)
print("PART 5: INNER CHIRALITY gamma_F = (-1)^{2j}")
print("=" * 72)

# 2I irreps: rho_0(1), rho_1(2), rho_2(3), rho_3(4), rho_4(5),
#             rho_5(6), rho_6(3'), rho_7(4'), rho_8(2')
# Spin-k/2 values: k = 0, 1, 2, 3, 4, 5, 2, 3, 1
# (-1)^{2j} = (-1)^k = +1,-1,+1,-1,+1,-1,+1,-1,-1
# Wait - need to be careful. The McKay graph bipartition:

# From chirality_mckay.md / exp351:
# WHITE (gamma_F = +1): rho_0(1), rho_4(5), rho_5(6), rho_6(3'), rho_8(2')
#                        dims: 1+5+6+3+2 = 17... no wait
# Hmm, let me use the correct assignment from MEMORY:
# WHITE(dim 16): rho_1, rho_4, rho_5, rho_6, rho_8
# BLACK(dim 14): rho_0, rho_2, rho_3, rho_7, rho_9

# Actually from MEMORY: WHITE(dim 16): rho_1,4,5,6,8 = {e,d,b,tau,mu}
# But rho indices in the 2I are 0-8 (9 irreps).
# Let me use the correct bipartition from the McKay tree:

# The McKay graph for 2I (affine E8) has 9 nodes.
# Bipartition (tree is bipartite):
# Using the standard E8 Dynkin convention where rho_0 is the affine node:
# WHITE: rho_0(1), rho_2(3), rho_4(5), rho_6(3'), rho_8(2')
# BLACK: rho_1(2), rho_3(4), rho_5(6), rho_7(4'), rho_9... wait, there are only 9 irreps

# Let me be precise. The 9 irreps of 2I and their dimensions:
irrep_dims = [1, 2, 3, 4, 5, 6, 3, 4, 2]
# Indices: 0,1,2,3,4,5,6,7,8
# Total dim: 1+2+3+4+5+6+3+4+2 = 30

# McKay graph bipartition (from the affine E8 tree):
# The affine E8 graph is:
#   0 - 1 - 2 - 3 - 4 - 5 - 6
#                       |
#                       7 - 8
# Wait, that's not right either. Let me use the standard convention.
# The affine E8 has nodes 0-8 with the branching at node 4 (dim 5):
#   0 - 1 - 2 - 3 - 4 - 5
#                   |
#                   6 - 7 - 8
# No... The E8 Dynkin diagram is:
#   1 - 2 - 3 - 4 - 5 - 6 - 7
#                   |
#                   8
# And the affine E8 (E8~) adds node 0 connected to node 8.
# But for the McKay graph of 2I, the structure might differ.

# From MEMORY: McKay graph branch at rho_9(6), legs (1,2,5)
# This means the graph has a branching node with 3 legs of lengths 1, 2, 5
# That's the E8 shape: legs 1, 2, 5 (or equivalently 2, 3, 6 nodes)

# For the bipartition of a tree: alternate black/white from any root
# Let's use the structure: main chain rho_0 - rho_1 - ... - rho_8
# with branch at some node.

# Actually, for the key result we just need dim(WHITE) - dim(BLACK):
# From MEMORY: dim(WHITE) = 16, dim(BLACK) = 14
# WHITE has rho indices where gamma_F = +1
# BLACK has rho indices where gamma_F = -1

# The bipartition assignment (from exp351):
# WHITE: rho indices {0, 2, 4, 6, 8} -> dims 1+3+5+3+2 = 14
# BLACK: rho indices {1, 3, 5, 7} -> dims 2+4+6+4 = 16

# Wait, that gives WHITE=14 and BLACK=16, which is the opposite.
# From MEMORY: WHITE(dim 16), BLACK(dim 14)
# So WHITE = {1, 3, 5, 7} (dims 2+4+6+4=16)
#    BLACK = {0, 2, 4, 6, 8} (dims 1+3+5+3+2=14)

# But gamma_F = (-1)^{2j} and MEMORY says:
# "integer spin=WHITE, half-integer=BLACK"
# The spin values k/2 for each rho:
# rho_0: k=0 (spin 0) -> (-1)^0 = +1
# rho_1: k=1 (spin 1/2) -> (-1)^1 = -1
# rho_2: k=2 (spin 1) -> (-1)^2 = +1
# rho_3: k=3 (spin 3/2) -> (-1)^3 = -1
# rho_4: k=4 (spin 2) -> (-1)^4 = +1
# rho_5: k=5 (spin 5/2) -> (-1)^5 = -1
# rho_6: k=2 (spin 1) -> (-1)^2 = +1    [NOTE: rho_6 has k=2 from Sym^2(std x std')]
# rho_7: k=3 (spin 3/2) -> (-1)^3 = -1
# rho_8: k=1 (spin 1/2) -> (-1)^1 = -1  [NOTE: rho_8 has k=1 from Galois conjugate]

# Wait - this gives:
# gamma_F = +1 for k even: rho_0(1), rho_2(3), rho_4(5), rho_6(3) -> total dim = 12
# gamma_F = -1 for k odd: rho_1(2), rho_3(4), rho_5(6), rho_7(4), rho_8(2) -> total dim = 18
# That doesn't match 16/14 either!

# The issue is the k values for rho_6, rho_7, rho_8.
# From MEMORY lesson: "McKay spin degree: 2I irreps have k from construction (Sym^k)"
# rho_6 = std x std' has k=2 (not 3/2)
# So rho_6 has k=2, rho_7 has k=3, rho_8 has k=1? No...

# Let me use the correct k values from the McKay correspondence:
# The irreps of 2I as restrictions of SU(2) irreps:
# rho_0 = Sym^0(C^2) = trivial, k=0
# rho_1 = Sym^1(C^2) = std, k=1  (the fundamental 2D rep)
# rho_2 = Sym^2(C^2), k=2
# rho_3 = Sym^3(C^2), k=3
# rho_4 = Sym^4(C^2), k=4
# rho_5 = Sym^5(C^2), k=5
# For k >= 5, the restriction starts to decompose:
# Sym^6 = rho_6 + rho_4, so rho_6 appears with k=6 (but also built from k=2)
# Actually the correct statement is that the McKay graph encodes
# the tensor product with the fundamental:
# rho_1 x rho_k = rho_{k-1} + rho_{k+1} (for the chain part)

# This is getting complicated. Let me just use the bipartition directly.
# For a tree graph, vertices at even distance from root = one color,
# odd distance = other color.

# The McKay graph structure (affine E8 Dynkin):
# Linear chain: 0-1-2-3-4-5 with branch at 4: 4-6, 6-7, 7-8
# Actually the standard affine E8 is:
# 0-1-2-3-4-5-6-7 with 8 branching from 4
# But for 2I McKay, let me just check distances from rho_0.

# McKay adjacency (from tensor product with std=rho_1):
# rho_1 x rho_0 = rho_1 -> edge 0-1
# rho_1 x rho_1 = rho_0 + rho_2 -> edges 1-0, 1-2
# rho_1 x rho_2 = rho_1 + rho_3 -> edges 2-1, 2-3
# rho_1 x rho_3 = rho_2 + rho_4 -> edges 3-2, 3-4
# rho_1 x rho_4 = rho_3 + rho_5 -> edges 4-3, 4-5
# rho_1 x rho_5 = rho_4 + rho_6 + rho_7 -> edges 5-4, 5-6, 5-7 (BRANCH!)
# rho_1 x rho_6 = rho_5 -> edge 6-5
# rho_1 x rho_7 = rho_5 + rho_8 -> edges 7-5, 7-8
# rho_1 x rho_8 = rho_7 -> edge 8-7

# So the graph is:
# 0-1-2-3-4-5-6 (main chain of length 6, legs from 5)
#             |
#             7-8

# Wait, that gives legs of length 1 (to 6) and 2 (to 7-8) from the branch at 5.
# And the main chain from 0 to 5 has length 5.
# So legs: 5, 1, 2 from the branch node (rho_5, dim 6).

# Bipartition from rho_0:
# Distance 0 (even): rho_0
# Distance 1 (odd): rho_1
# Distance 2 (even): rho_2
# Distance 3 (odd): rho_3
# Distance 4 (even): rho_4
# Distance 5 (odd): rho_5
# Distance 6 (even): rho_6, rho_7 (both adjacent to rho_5)
# Wait - rho_7 is connected to rho_5, so distance = 6? No:
# rho_7 is connected to rho_5 (distance 5+1=6) AND to rho_8
# rho_8 is at distance 7 from rho_0 (through 0-1-2-3-4-5-7-8)
# rho_6 is at distance 6 from rho_0 (through 0-1-2-3-4-5-6)

# Bipartition:
# EVEN distance (WHITE): rho_0(1), rho_2(3), rho_4(5), rho_6(3), rho_7(4)
# dims: 1+3+5+3+4 = 16
# ODD distance (BLACK): rho_1(2), rho_3(4), rho_5(6), rho_8(2)
# dims: 2+4+6+2 = 14

# YES! This gives WHITE=16, BLACK=14. CORRECT!

white_irreps = [0, 2, 4, 6, 7]  # dims 1+3+5+3+4 = 16
black_irreps = [1, 3, 5, 8]     # dims 2+4+6+2 = 14

dim_white = sum(irrep_dims[i] for i in white_irreps)
dim_black = sum(irrep_dims[i] for i in black_irreps)

print(f"\n  McKay graph bipartition:")
print(f"    WHITE (gamma_F=+1): rho_{white_irreps}, dims = "
      f"{[irrep_dims[i] for i in white_irreps]}, total = {dim_white}")
print(f"    BLACK (gamma_F=-1): rho_{black_irreps}, dims = "
      f"{[irrep_dims[i] for i in black_irreps]}, total = {dim_black}")
print(f"    dim(WHITE) - dim(BLACK) = {dim_white - dim_black}")
print(f"    dim(WHITE) = 16 = (a1-1)^2 = {(a1-1)**2}")

# Tr(gamma_F) = dim(WHITE) - dim(BLACK)
tr_gamma_F = dim_white - dim_black
print(f"\n  Tr(gamma_F) = {tr_gamma_F}")
print(f"  This is NOT zero! (unlike Tr(gamma_form) = 0)")
print(f"  The imbalance comes from the non-symmetric tree structure.")

# ============================================================================
# PART 6: PRODUCT CHIRALITY ANALYSIS
# ============================================================================
print("\n" + "=" * 72)
print("PART 6: PRODUCT CHIRALITY ANALYSIS")
print("=" * 72)

dim_F = 30  # sum of all irrep dims
dim_product = dim_H * dim_F

print(f"  dim(H_M) = {dim_H}")
print(f"  dim(H_F) = {dim_F}")
print(f"  dim(H_total) = {dim_product}")

# H+ and H-
# gamma = gamma_form x gamma_F
# H+ = (gamma_form=+1, gamma_F=+1) + (gamma_form=-1, gamma_F=-1)
# H- = (gamma_form=+1, gamma_F=-1) + (gamma_form=-1, gamma_F=+1)

dim_even = N + F       # 120 + 1200 = 1320 (gamma_form = +1)
dim_odd = E + C_cells  # 720 + 600 = 1320 (gamma_form = -1)

dim_H_plus = dim_even * dim_white + dim_odd * dim_black
dim_H_minus = dim_even * dim_black + dim_odd * dim_white

print(f"\n  dim(H+) = {dim_even}*{dim_white} + {dim_odd}*{dim_black}"
      f" = {dim_H_plus}")
print(f"  dim(H-) = {dim_even}*{dim_black} + {dim_odd}*{dim_white}"
      f" = {dim_H_minus}")
print(f"  H+ = H- ? {dim_H_plus == dim_H_minus}")
print(f"  (This follows from dim_even = dim_odd = 1320"
      f" and dim_white + dim_black = 30)")

# ============================================================================
# PART 7: T3 AND GALOIS KERNEL INTERACTION
# ============================================================================
print("\n" + "=" * 72)
print("PART 7: T3 ON GALOIS KERNEL")
print("=" * 72)

# The Galois kernel has 9 modes (from Box = b1*A_fiber - A).
# How does T3 act on these modes?

# First, find a Hopf fibration
def quat_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])

def find_vertex_index(verts, q, tol=1e-6):
    dists = np.linalg.norm(verts - q, axis=1)
    idx = np.argmin(dists)
    return idx if dists[idx] < tol else -1

# Find order-10 element
gen_idx = -1
for i in range(N):
    if abs(verts[i, 0] - phi/2) < 1e-6:
        g = verts[i]
        power = g.copy()
        ok = True
        for k in range(2, 11):
            power = quat_mult(power, g)
            if k == 5 and not np.allclose(power, [-1,0,0,0], atol=1e-6):
                ok = False
                break
            if k == 10 and not np.allclose(power, [1,0,0,0], atol=1e-6):
                ok = False
        if ok:
            gen_idx = i
            break

if gen_idx >= 0:
    g = verts[gen_idx]
    subgroup = []
    power = np.array([1.0, 0.0, 0.0, 0.0])
    for k in range(10):
        idx = find_vertex_index(verts, power)
        subgroup.append(idx)
        power = quat_mult(power, g)

    # Build fibers as cosets
    assigned = np.full(N, -1, dtype=int)
    fibers = []
    fid = 0
    for i in range(N):
        if assigned[i] >= 0:
            continue
        coset = []
        for si in subgroup:
            q_prod = quat_mult(verts[i], verts[si])
            idx = find_vertex_index(verts, q_prod)
            if idx >= 0 and assigned[idx] < 0:
                coset.append(idx)
                assigned[idx] = fid
        fibers.append(coset)
        fid += 1

    # Build A_fiber and Box
    A_fiber = np.zeros((N, N), dtype=float)
    for fiber in fibers:
        for i in fiber:
            for j in fiber:
                if i != j and A[i, j] > 0.5:
                    A_fiber[i, j] = 1.0

    Box = b1 * A_fiber - A
    evals_box, evecs_box = eigh(Box)
    kernel_mask = np.abs(evals_box) < TOL
    kernel_vecs = evecs_box[:, kernel_mask]

    # Project T3 onto kernel
    T3_kernel = kernel_vecs.T @ T3_diag @ kernel_vecs
    t3_kernel_evals = np.sort(np.linalg.eigvalsh(T3_kernel))

    print(f"  T3 eigenvalues on Galois kernel (9 modes):")
    t3_kr = np.round(t3_kernel_evals, 6)
    unique_t3k, cnt_t3k = np.unique(t3_kr, return_counts=True)
    for v, c in zip(unique_t3k, cnt_t3k):
        print(f"    T3 = {v:+.6f}, mult = {c}")

    # Check if T3 splits kernel into +/- sectors
    n_pos_k = np.sum(t3_kernel_evals > TOL)
    n_zero_k = np.sum(np.abs(t3_kernel_evals) <= TOL)
    n_neg_k = np.sum(t3_kernel_evals < -TOL)
    print(f"\n  Kernel T3 split: {n_pos_k} positive, {n_zero_k} zero, {n_neg_k} negative")

    # Check Tr(T3) on kernel
    tr_T3_kernel = np.trace(T3_kernel)
    print(f"  Tr(T3|kernel) = {tr_T3_kernel:.6f}")

    # Project Hopf coordinates onto kernel
    for name, H_diag in [("H1", H1_diag), ("H2", H2_diag), ("H3", H3_diag)]:
        H_kernel = kernel_vecs.T @ H_diag @ kernel_vecs
        h_evals = np.sort(np.linalg.eigvalsh(H_kernel))
        h_tr = np.trace(H_kernel)
        print(f"  Tr({name}|kernel) = {h_tr:.6f}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

print("""
RESULTS:

1. SPECTRAL SUPERTRACE: Tr(gamma * D^{2k}) = 0 on product space.
   REASON: Str_form(D_M^{2k}) = 0 (spectral balance on S3).
   CONSEQUENCE: Spectral action is identical on H+ and H-.
   STATUS: NEGATIVE for deriving chirality from spectral action on S3.

2. Tr(gamma_F) = 2 (not zero): the internal chirality is unbalanced
   because the McKay tree is not symmetric. This is a KEY feature
   but it doesn't help with the spectral supertrace because
   chi(S3) = 0 kills the product.

3. H+ = H- = 39600 exactly (balanced by dim_even = dim_odd).

4. T3 Hopf quadrupole: does NOT commute with A. Provides a
   "charge-like" splitting of vertices but is not a conserved
   quantity on the graph.

5. The chirality mechanism likely requires either:
   (a) 4D continuum limit where chi(M4) != 0
   (b) Gauge algebra structure via inner fluctuations
   (c) The Hopf T3 as an external operator selecting the coupling

CONCLUSION: OPEN-2 is NOT RESOLVED by the spectral action on S3.
The 3D geometry cannot distinguish H+ from H- because S3 has
no chirality anomaly. This is a fundamental limitation of the
discrete 3D approach, not a computational gap.
""")

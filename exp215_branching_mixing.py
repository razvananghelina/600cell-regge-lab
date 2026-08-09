"""
EXP-215: Character Table of 2I and Branching Rules to 2T
=========================================================
PURPOSE: Compute the full character table of binary icosahedral group 2I
         and branching rules for restriction to binary tetrahedral 2T.
         Test whether 3 two-dim irreps of 2T = 3 generations.

KEY IDEA:
- 600-cell = Cayley graph of 2I (120 elements)
- 2I has 9 irreps -> McKay correspondence -> affine E8 Dynkin
- 2T (order 24) has 7 irreps: dims {1,1,1,2,2,2,3}
- The THREE 2-dim irreps of 2T could be 3 generations
- Branching 2I -> 2T gives decomposition structure
- Mixing angles from multiplicity space overlaps

METHOD:
1. Construct 2I as quaternions, find conjugacy classes
2. Build Cayley graph class matrices -> eigendecompose -> character table
3. Identify 2T subgroup, compute its character table
4. Compute branching rules
5. Analyze 3x3 sub-matrices for CKM/PMNS structure
"""

import numpy as np
from itertools import product as iprod

PHI = (1 + np.sqrt(5)) / 2
PSI = 1 - PHI  # = -1/phi
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179

print("=" * 70)
print("EXP-215: CHARACTER TABLE OF 2I AND BRANCHING TO 2T")
print("=" * 70)

# ======================================================================
# PART 1: Construct 2I (120 unit quaternions)
# ======================================================================
print("\n--- PART 1: Constructing 2I ---")

def qmul(q1, q2):
    """Quaternion multiplication: q = w + xi + yj + zk"""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return (
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    )

def qinv(q):
    """Inverse of unit quaternion = conjugate"""
    w, x, y, z = q
    return (w, -x, -y, -z)

# Generate 120 elements
elements = []

# Type A: 8 vertices, one coord +-1
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s
        elements.append(tuple(v))

# Type B: 16 vertices, all +-1/2
for signs in iprod([0.5, -0.5], repeat=4):
    elements.append(tuple(signs))

# Type C: 96 vertices, even permutations of (0, +-1/2, +-phi/2, +-1/(2phi))
even_perms = [
    (0,1,2,3), (0,2,3,1), (0,3,1,2),
    (1,0,3,2), (1,2,0,3), (1,3,2,0),
    (2,0,1,3), (2,1,3,0), (2,3,0,1),
    (3,0,2,1), (3,1,0,2), (3,2,1,0)
]
base_vals = [0.0, 0.5, PHI/2, 1/(2*PHI)]

for perm in even_perms:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                v = [0.0, 0.0, 0.0, 0.0]
                v[perm[0]] = 0.0
                v[perm[1]] = s1 * base_vals[1]
                v[perm[2]] = s2 * base_vals[2]
                v[perm[3]] = s3 * base_vals[3]
                elements.append(tuple(v))

# Remove duplicates (some type C might overlap)
unique = []
for v in elements:
    is_dup = False
    for u in unique:
        if sum((v[k]-u[k])**2 for k in range(4)) < 1e-10:
            is_dup = True
            break
    if not is_dup:
        unique.append(v)
elements = unique
print("Generated %d unique elements" % len(elements))

# Verify unit quaternions
for q in elements:
    norm2 = sum(c**2 for c in q)
    assert abs(norm2 - 1.0) < 1e-10

# Build multiplication table
print("Building multiplication table...")
N = len(elements)
mult_table = np.zeros((N, N), dtype=int)

# Precompute element array for fast lookup
elem_arr = np.array(elements)

for i in range(N):
    for j in range(N):
        prod = qmul(elements[i], elements[j])
        prod_arr = np.array(prod)
        dists = np.sum((elem_arr - prod_arr)**2, axis=1)
        best_k = np.argmin(dists)
        assert dists[best_k] < 1e-8, "Product not in group!"
        mult_table[i, j] = best_k

# Find identity
identity_idx = -1
for i in range(N):
    if all(mult_table[i, j] == j for j in range(N)):
        identity_idx = i
        break
print("Identity: index %d = %s" % (identity_idx, str(elements[identity_idx])))

# Build inverse table
inv_table = np.zeros(N, dtype=int)
for i in range(N):
    for k in range(N):
        if mult_table[i, k] == identity_idx:
            inv_table[i] = k
            break

# ======================================================================
# PART 2: Conjugacy Classes of 2I
# ======================================================================
print("\n--- PART 2: Conjugacy Classes of 2I ---")

visited = [False] * N
classes_2I = []

for g in range(N):
    if visited[g]:
        continue
    cc = set()
    for q in range(N):
        q_inv = inv_table[q]
        conj = mult_table[mult_table[q, g], q_inv]
        cc.add(conj)
    for c in cc:
        visited[c] = True
    classes_2I.append(sorted(cc))

# Sort by size then by w-coordinate of representative
classes_2I.sort(key=lambda cc: (len(cc), -elements[cc[0]][0]))

print("Found %d conjugacy classes:" % len(classes_2I))
class_sizes_2I = []
for i, cc in enumerate(classes_2I):
    rep = elements[cc[0]]
    g = cc[0]
    order = 1
    current = g
    while current != identity_idx and order < 200:
        current = mult_table[current, g]
        order += 1
    class_sizes_2I.append(len(cc))
    print("  C_%d: size %2d, order %2d, rep=(%.4f, %.4f, %.4f, %.4f)"
          % (i, len(cc), order, rep[0], rep[1], rep[2], rep[3]))

print("Total elements: %d (should be 120)" % sum(class_sizes_2I))

# ======================================================================
# PART 3: Class matrices and eigendecomposition -> Character table
# ======================================================================
print("\n--- PART 3: Character Table of 2I ---")

# Build class membership
class_of = np.zeros(N, dtype=int)
for k, cc in enumerate(classes_2I):
    for g in cc:
        class_of[g] = k

n_classes = len(classes_2I)

# Build association matrices R_k[i,j] = 1 iff i^{-1}*j in C_k
A_class = np.zeros((n_classes, N, N), dtype=float)
for i in range(N):
    i_inv = inv_table[i]
    for j in range(N):
        g = mult_table[i_inv, j]
        k = class_of[g]
        A_class[k, i, j] = 1.0

# All class matrices commute and are normal -> simultaneously diagonalizable
# Use a random real linear combination to split degeneracies
np.random.seed(12345)
coeffs = np.random.randn(n_classes)
A_combo = sum(coeffs[k] * A_class[k] for k in range(n_classes))

# Check symmetry: R_k is symmetric iff C_k is self-inverse
for k in range(n_classes):
    sym_err = np.max(np.abs(A_class[k] - A_class[k].T))
    if sym_err > 1e-10:
        print("  WARNING: A_%d not symmetric (err=%.2e)" % (k, sym_err))

# A_combo should be symmetric since all class matrices are symmetric for 2I
# (every element is conjugate to its inverse in 2I)
evals_combo, evecs = np.linalg.eigh(A_combo)

# Find distinct eigenvalues -> eigenspaces
distinct_combo = []
mults = []
tol = 1e-6
for ev in evals_combo:
    found = False
    for i_d, dev in enumerate(distinct_combo):
        if abs(ev - dev) < tol:
            mults[i_d] += 1
            found = True
            break
    if not found:
        distinct_combo.append(ev)
        mults.append(1)

dims = [int(round(np.sqrt(m))) for m in mults]
print("Irrep dimensions: %s" % dims)
print("Sum of dim^2 = %d (should be 120)" % sum(d**2 for d in dims))

# Build eigenspaces
eigenspaces = []
idx = 0
for i_d in range(len(distinct_combo)):
    m = mults[i_d]
    eigenspaces.append(evecs[:, idx:idx+m])
    idx += m

# P-matrix: eigenvalue of class matrix k on eigenspace i
P_matrix = np.zeros((len(distinct_combo), n_classes))
for i_d in range(len(distinct_combo)):
    V_i = eigenspaces[i_d]
    for k in range(n_classes):
        v = V_i[:, 0]
        Av = A_class[k] @ v
        P_matrix[i_d, k] = np.dot(v, Av) / np.dot(v, v)

# Character table: chi_i(C_k) = P[i,k] * dim_i / |C_k|
char_table = np.zeros((len(distinct_combo), n_classes))
for i_d in range(len(distinct_combo)):
    for k in range(n_classes):
        char_table[i_d, k] = P_matrix[i_d, k] * dims[i_d] / class_sizes_2I[k]

# Sort rows by dimension (ascending) to match McKay convention
sort_idx = np.argsort(dims)
# For same dimension, sort by eigenvalue on adjacency matrix
# Find which class matrix is the adjacency (row sum = 12)
adj_k = -1
for k in range(n_classes):
    if class_sizes_2I[k] == 12:
        # Pick the one with largest positive eigenvalue
        max_eval = max(P_matrix[i_d, k] for i_d in range(len(distinct_combo)))
        if max_eval > 10:  # Should be 12
            adj_k = k
            break

if adj_k >= 0:
    print("Adjacency matrix is class %d" % adj_k)
    # Sort: first by dim, then by adjacency eigenvalue (descending)
    sort_key = [(dims[i], -P_matrix[i, adj_k]) for i in range(len(distinct_combo))]
    sort_idx = sorted(range(len(distinct_combo)), key=lambda i: sort_key[i])

dims_sorted = [dims[i] for i in sort_idx]
char_table_sorted = char_table[sort_idx, :]
P_sorted = P_matrix[sort_idx, :]

print("\nCharacter Table of 2I (sorted by dimension):")
print("             ", end="")
for k in range(n_classes):
    print("C_%d(%d) " % (k, class_sizes_2I[k]), end="")
print()
for i_row, i_orig in enumerate(sort_idx):
    d = dims_sorted[i_row]
    print("rho_%d (d=%d): " % (i_row, d), end="")
    for k in range(n_classes):
        val = char_table_sorted[i_row, k]
        # Try a + b*phi decomposition
        best = None
        for a_try in range(-8, 9):
            for b_try in range(-8, 9):
                if abs(a_try + b_try*PHI - val) < 1e-4:
                    best = (a_try, b_try)
                    break
            if best:
                break
        if best:
            a_v, b_v = best
            if b_v == 0:
                print("%5d   " % a_v, end="")
            elif a_v == 0:
                print(" %d*phi " % b_v, end="")
            else:
                s = "+" if b_v > 0 else ""
                print("%d%s%dphi" % (a_v, s, b_v), end="")
        else:
            print("%7.3f " % val, end="")
    print()

# Verify orthogonality
print("\nOrthogonality check (<chi_i, chi_j> should be delta_ij):")
max_off_diag = 0
for i in range(len(dims_sorted)):
    for j in range(i, len(dims_sorted)):
        inner = sum(class_sizes_2I[k] * char_table_sorted[i, k] *
                     char_table_sorted[j, k] for k in range(n_classes)) / 120.0
        if i == j:
            if abs(inner - 1.0) > 0.01:
                print("  WARNING: <chi_%d, chi_%d> = %.4f (should be 1)" % (i, j, inner))
        else:
            if abs(inner) > max_off_diag:
                max_off_diag = abs(inner)
if max_off_diag < 0.01:
    print("  All off-diagonal < 0.01: PASS")
else:
    print("  Max off-diagonal = %.4f" % max_off_diag)

# Identify eigenvalues on adjacency matrix
print("\n600-cell eigenvalues (from adjacency class matrix):")
if adj_k >= 0:
    for i_row in range(len(dims_sorted)):
        ev = P_sorted[sort_idx.index(sort_idx[i_row]) if isinstance(sort_idx, list) else i_row, adj_k]
        # Actually use sorted P
        ev = char_table_sorted[i_row, adj_k] * class_sizes_2I[adj_k] / dims_sorted[i_row]
        # Try a+b*phi
        best = None
        for a_try in range(-10, 13):
            for b_try in range(-10, 11):
                if abs(a_try + b_try*PHI - ev) < 1e-3:
                    best = (a_try, b_try)
                    break
            if best:
                break
        label = "%d+%d*phi" % best if best else "%.4f" % ev
        print("  rho_%d (dim %d): lambda = %10.6f = %s, mult = %d"
              % (i_row, dims_sorted[i_row], ev, label, dims_sorted[i_row]**2))

# ======================================================================
# PART 4: Binary Tetrahedral Subgroup 2T
# ======================================================================
print("\n--- PART 4: Binary Tetrahedral Subgroup 2T ---")

# 2T = Type A (8) + Type B (16) = 24 elements
type_A = []
type_B = []
type_C = []
for i, q in enumerate(elements):
    n_zero = sum(1 for c in q if abs(c) < 1e-10)
    if n_zero == 3:
        type_A.append(i)
    elif n_zero == 0 and all(abs(abs(c) - 0.5) < 1e-10 for c in q):
        type_B.append(i)
    else:
        type_C.append(i)

subgroup_2T = sorted(type_A + type_B)
print("2T: %d elements (8 type-A + 16 type-B)" % len(subgroup_2T))

# Verify closure
closed = True
for i in subgroup_2T:
    for j in subgroup_2T:
        if mult_table[i, j] not in subgroup_2T:
            closed = False
            break
    if not closed:
        break
print("2T closure: %s" % ("PASS" if closed else "FAIL"))

# Conjugacy classes of 2T
visited_2T = set()
classes_2T = []
for g in subgroup_2T:
    if g in visited_2T:
        continue
    cc = set()
    for q in subgroup_2T:
        q_inv = inv_table[q]
        conj = mult_table[mult_table[q, g], q_inv]
        cc.add(conj)
    visited_2T |= cc
    classes_2T.append(sorted(cc))

classes_2T.sort(key=lambda cc: (len(cc), -elements[cc[0]][0]))

print("2T conjugacy classes:")
class_sizes_2T = []
for i, cc in enumerate(classes_2T):
    rep = elements[cc[0]]
    g = cc[0]
    order = 1
    current = g
    while current != identity_idx and order < 200:
        current = mult_table[current, g]
        order += 1
    class_sizes_2T.append(len(cc))
    print("  D_%d: size %d, order %d, rep=(%.4f, %.4f, %.4f, %.4f)"
          % (i, len(cc), order, rep[0], rep[1], rep[2], rep[3]))

# ======================================================================
# PART 5: Character Table of 2T
# ======================================================================
print("\n--- PART 5: Character Table of 2T ---")

n_2T = len(subgroup_2T)
n_classes_2T = len(classes_2T)

# Local index mapping
local_idx = {}
for loc, glob in enumerate(subgroup_2T):
    local_idx[glob] = loc

# 2T multiplication table
mult_2T = np.zeros((n_2T, n_2T), dtype=int)
for i_loc in range(n_2T):
    for j_loc in range(n_2T):
        i_glob = subgroup_2T[i_loc]
        j_glob = subgroup_2T[j_loc]
        prod_glob = mult_table[i_glob, j_glob]
        mult_2T[i_loc, j_loc] = local_idx[prod_glob]

# Find identity in local indexing
id_2T = local_idx[identity_idx]

# 2T inverse table
inv_2T = np.zeros(n_2T, dtype=int)
for i in range(n_2T):
    for k in range(n_2T):
        if mult_2T[i, k] == id_2T:
            inv_2T[i] = k
            break

# Class membership in 2T
class_of_2T = np.zeros(n_2T, dtype=int)
for k, cc in enumerate(classes_2T):
    for g in cc:
        class_of_2T[local_idx[g]] = k

# Build class matrices for 2T
A_class_2T = np.zeros((n_classes_2T, n_2T, n_2T), dtype=float)
for i in range(n_2T):
    i_inv = inv_2T[i]
    for j in range(n_2T):
        g = mult_2T[i_inv, j]
        k = class_of_2T[g]
        A_class_2T[k, i, j] = 1.0

# Check symmetry of 2T class matrices
print("2T class matrix symmetry:")
all_sym = True
for k in range(n_classes_2T):
    sym_err = np.max(np.abs(A_class_2T[k] - A_class_2T[k].T))
    if sym_err > 1e-10:
        print("  A_%d: NOT symmetric (err=%.2e)" % (k, sym_err))
        all_sym = False
if all_sym:
    print("  All symmetric: PASS")

# 2T has complex characters (omega = e^{2pi*i/3})
# Non-symmetric class matrices => use general eigendecomposition
# Or use the LEFT regular representation properly

# Alternative: use the known character table of 2T
# 2T has 7 irreps: dims 1,1,1,2,2,2,3
# Characters involve omega = e^{2pi*i/3}
# But since we want real branching rules, we can pair complex conjugate irreps

# Let's compute using simultaneous diagonalization
# Use random complex combination to break degeneracies
np.random.seed(54321)
coeffs_2T = np.random.randn(n_classes_2T)
A_combo_2T = sum(coeffs_2T[k] * A_class_2T[k] for k in range(n_classes_2T))

# Check if symmetric
is_sym = np.max(np.abs(A_combo_2T - A_combo_2T.T)) < 1e-10
if is_sym:
    evals_2T, evecs_2T = np.linalg.eigh(A_combo_2T)
else:
    print("  Using general eigensolver (non-symmetric class matrices)")
    evals_2T, evecs_2T = np.linalg.eig(A_combo_2T)
    # Sort by real part
    order_2T = np.argsort(evals_2T.real)
    evals_2T = evals_2T[order_2T]
    evecs_2T = evecs_2T[:, order_2T]

# Find distinct eigenvalues
distinct_2T = []
mults_2T = []
for ev in evals_2T:
    found = False
    for i_d, dev in enumerate(distinct_2T):
        if abs(ev - dev) < 1e-4:
            mults_2T[i_d] += 1
            found = True
            break
    if not found:
        distinct_2T.append(ev)
        mults_2T.append(1)

dims_2T = [int(round(np.sqrt(m))) for m in mults_2T]
print("\n2T irrep dimensions: %s" % dims_2T)
print("Sum of dim^2 = %d (should be 24)" % sum(d**2 for d in dims_2T))

# Build eigenspaces for 2T
eigenspaces_2T = []
idx = 0
for i_d in range(len(distinct_2T)):
    m = mults_2T[i_d]
    eigenspaces_2T.append(evecs_2T[:, idx:idx+m])
    idx += m

# Compute P-matrix for 2T
P_2T = np.zeros((len(distinct_2T), n_classes_2T), dtype=complex)
for i_d in range(len(distinct_2T)):
    V_i = eigenspaces_2T[i_d]
    for k in range(n_classes_2T):
        v = V_i[:, 0]
        Av = A_class_2T[k] @ v
        if np.abs(np.dot(np.conj(v), v)) > 1e-10:
            P_2T[i_d, k] = np.dot(np.conj(v), Av) / np.dot(np.conj(v), v)

# Character table of 2T
char_table_2T = np.zeros((len(distinct_2T), n_classes_2T), dtype=complex)
for i_d in range(len(distinct_2T)):
    for k in range(n_classes_2T):
        char_table_2T[i_d, k] = P_2T[i_d, k] * dims_2T[i_d] / class_sizes_2T[k]

# Sort by dimension
sort_idx_2T = sorted(range(len(distinct_2T)), key=lambda i: (dims_2T[i], distinct_2T[i].real))
dims_2T_sorted = [dims_2T[i] for i in sort_idx_2T]
char_table_2T_sorted = char_table_2T[sort_idx_2T, :]

print("\nCharacter Table of 2T:")
print("              ", end="")
for k in range(n_classes_2T):
    print("D_%d(%d)  " % (k, class_sizes_2T[k]), end="")
print()
for i_row, i_orig in enumerate(sort_idx_2T):
    d = dims_2T_sorted[i_row]
    print("psi_%d (d=%d): " % (i_row, d), end="")
    for k in range(n_classes_2T):
        val = char_table_2T_sorted[i_row, k]
        if abs(val.imag) < 1e-4:
            print("%7.3f " % val.real, end="")
        else:
            print("%5.2f%+5.2fi" % (val.real, val.imag), end="")
    print()

# Verify orthogonality
print("\n2T orthogonality check:")
max_off = 0
for i in range(len(dims_2T_sorted)):
    for j in range(i, len(dims_2T_sorted)):
        inner = sum(class_sizes_2T[k] * char_table_2T_sorted[i, k] *
                     np.conj(char_table_2T_sorted[j, k])
                     for k in range(n_classes_2T)) / 24.0
        if i == j and abs(inner - 1.0) > 0.05:
            print("  <psi_%d, psi_%d> = %.4f (should be 1)" % (i, j, abs(inner)))
        if i != j and abs(inner) > max_off:
            max_off = abs(inner)
if max_off < 0.05:
    print("  Max off-diagonal < 0.05: PASS")
else:
    print("  Max off-diagonal = %.4f" % max_off)

# ======================================================================
# PART 6: Branching Rules 2I -> 2T
# ======================================================================
print("\n--- PART 6: Branching Rules 2I -> 2T ---")

# Map 2T classes to 2I classes
print("Class embedding 2T -> 2I:")
class_map = []
for k_2T, cc_2T in enumerate(classes_2T):
    rep = cc_2T[0]
    k_2I = class_of[rep]
    class_map.append(k_2I)
    print("  D_%d (size %d) -> C_%d (size %d)"
          % (k_2T, len(cc_2T), k_2I, class_sizes_2I[k_2I]))

# Check if a 2I class splits into multiple 2T classes
print("\n2I class splitting under 2T:")
for k_2I in range(n_classes):
    sub_classes = [k_2T for k_2T in range(n_classes_2T) if class_map[k_2T] == k_2I]
    if len(sub_classes) > 0:
        sizes = [class_sizes_2T[k_2T] for k_2T in sub_classes]
        print("  C_%d (size %d) -> %d sub-classes: sizes %s"
              % (k_2I, class_sizes_2I[k_2I], len(sub_classes), sizes))

# Restrict 2I characters to 2T
restricted = np.zeros((len(dims_sorted), n_classes_2T))
for i_2I in range(len(dims_sorted)):
    for k_2T in range(n_classes_2T):
        restricted[i_2I, k_2T] = char_table_sorted[i_2I, class_map[k_2T]]

print("\nRestricted 2I characters on 2T classes:")
for i in range(len(dims_sorted)):
    print("  rho_%d|_{2T} (d=%d): [" % (i, dims_sorted[i]), end="")
    for k in range(n_classes_2T):
        print("%.3f " % restricted[i, k], end="")
    print("]")

# Decompose into 2T irreps using inner products
# m_{ij} = <chi_i|_{2T}, psi_j>_{2T} = (1/24) sum_k |D_k| chi_i(D_k) conj(psi_j(D_k))
branching = np.zeros((len(dims_sorted), len(dims_2T_sorted)), dtype=complex)
for i in range(len(dims_sorted)):
    for j in range(len(dims_2T_sorted)):
        inner = 0
        for k in range(n_classes_2T):
            inner += class_sizes_2T[k] * restricted[i, k] * np.conj(char_table_2T_sorted[j, k])
        branching[i, j] = inner / 24.0

print("\nBranching matrix (2I irreps -> 2T irreps):")
print("             ", end="")
for j in range(len(dims_2T_sorted)):
    print("psi_%d(d=%d) " % (j, dims_2T_sorted[j]), end="")
print()
for i in range(len(dims_sorted)):
    print("rho_%d (d=%d):" % (i, dims_sorted[i]), end="")
    for j in range(len(dims_2T_sorted)):
        val = branching[i, j]
        if abs(val.imag) < 0.01:
            print("  %5.2f  " % val.real, end="")
        else:
            print(" %4.1f%+4.1fi" % (val.real, val.imag), end="")
    print()

# Round to integers and verify
branching_int = np.round(branching.real).astype(int)
print("\nBranching (integer coefficients):")
print("             ", end="")
for j in range(len(dims_2T_sorted)):
    print("psi_%d(%d) " % (j, dims_2T_sorted[j]), end="")
print()
for i in range(len(dims_sorted)):
    dim_check = sum(branching_int[i, j] * dims_2T_sorted[j]
                    for j in range(len(dims_2T_sorted)))
    status = "OK" if dim_check == dims_sorted[i] else "FAIL(%d)" % dim_check
    print("rho_%d (d=%d): " % (i, dims_sorted[i]), end="")
    for j in range(len(dims_2T_sorted)):
        print("  %3d    " % branching_int[i, j], end="")
    print(" [dim: %s]" % status)

# ======================================================================
# PART 7: Three 2D Irreps and Generation Structure
# ======================================================================
print("\n--- PART 7: Generation Structure from 2D Irreps ---")

two_dim_indices = [j for j in range(len(dims_2T_sorted)) if dims_2T_sorted[j] == 2]
print("Three 2D irreps of 2T: indices %s" % two_dim_indices)

if len(two_dim_indices) == 3:
    print("\nBranching into THREE 2-dim irreps (= 3 generations?):")
    print("             ", end="")
    for j in two_dim_indices:
        print("  Gen_%d  " % (two_dim_indices.index(j) + 1), end="")
    print()
    for i in range(len(dims_sorted)):
        print("rho_%d (d=%d): " % (i, dims_sorted[i]), end="")
        for j in two_dim_indices:
            print("   %3d   " % branching_int[i, j], end="")
        print()

    # The 3x3 sub-matrix for specific E8 legs
    print("\n--- E8 Dynkin Leg Analysis ---")
    # McKay dims: 1,2,3,4,5,6,4,3,2 for nodes 0-8
    # Our sorted dims should match this
    mckay_labels = {1: 'node0', 2: 'node1/8', 3: 'node2/7',
                    4: 'node3/6', 5: 'node4', 6: 'node5'}

    # Group by dimension
    by_dim = {}
    for i in range(len(dims_sorted)):
        d = dims_sorted[i]
        if d not in by_dim:
            by_dim[d] = []
        by_dim[d].append(i)

    print("Irreps by dimension:")
    for d in sorted(by_dim.keys()):
        print("  dim %d: %d irreps (indices %s)" % (d, len(by_dim[d]), by_dim[d]))
        for i in by_dim[d]:
            gen_content = [branching_int[i, j] for j in two_dim_indices]
            print("    rho_%d: generation content = %s (total 2D = %d)"
                  % (i, gen_content, sum(gen_content)))

    # Leg B: dims 3, 2, 1 -> 3 nodes
    # Leg C: dims 4, 3, 2 -> 3 nodes
    # For paired dims (2,3,4), distinguish via eigenvalue
    print("\n--- Leg Assignments via Galois Symmetry ---")
    if adj_k >= 0:
        for d in [2, 3, 4]:
            if len(by_dim.get(d, [])) == 2:
                i1, i2 = by_dim[d]
                ev1 = P_sorted[sort_idx.index(i1) if isinstance(sort_idx, list) else i1, adj_k]
                ev2 = P_sorted[sort_idx.index(i2) if isinstance(sort_idx, list) else i2, adj_k]
                # Actually recompute from char table
                ev1 = char_table_sorted[i1, adj_k] * class_sizes_2I[adj_k] / dims_sorted[i1]
                ev2 = char_table_sorted[i2, adj_k] * class_sizes_2I[adj_k] / dims_sorted[i2]
                print("  dim %d: rho_%d (lambda=%.4f) vs rho_%d (lambda=%.4f)"
                      % (d, i1, ev1, i2, ev2))
                # Galois: phi-type (positive phi coeff) = Leg B, phi'-type = Leg C
                # 6phi ~ 9.71 vs 6-6phi ~ -3.71 -> larger = Leg B
                # 4phi ~ 6.47 vs 4-4phi ~ -2.47 -> larger = Leg B

    # Build leg assignments
    # Node 0 (dim 1) -> Leg B endpoint
    # Node 1 (dim 2, higher eval) -> Leg B
    # Node 2 (dim 3, higher eval) -> Leg B
    # Node 3 (dim 4, higher eval) -> Branch
    # Node 4 (dim 5) -> Leg A
    # Node 5 (dim 6) -> Leg A
    # Node 6 (dim 4, lower eval) -> Leg C start
    # Node 7 (dim 3, lower eval) -> Leg C
    # Node 8 (dim 2, lower eval) -> Leg C endpoint

    # Extract 3x3 sub-matrices for each leg
    print("\n--- 3x3 Generation Sub-Matrices ---")

    # Sort paired irreps by eigenvalue
    leg_B = []  # nodes 0, 1, 2 -> dims 1, 2, 3
    leg_C = []  # nodes 8, 7, 6 -> dims 2, 3, 4
    leg_A = []  # nodes 4, 5 -> dims 5, 6
    branch = [] # node 3 -> dim 4

    # dim 1: unique
    if 1 in by_dim:
        leg_B.append(by_dim[1][0])

    # dims 2, 3, 4: paired, sort by eigenvalue
    for d in [2, 3, 4]:
        if d in by_dim and len(by_dim[d]) == 2:
            i1, i2 = by_dim[d]
            ev1 = char_table_sorted[i1, adj_k] * class_sizes_2I[adj_k] / dims_sorted[i1]
            ev2 = char_table_sorted[i2, adj_k] * class_sizes_2I[adj_k] / dims_sorted[i2]
            if ev1 > ev2:
                phi_type, phip_type = i1, i2
            else:
                phi_type, phip_type = i2, i1

            if d == 2:
                leg_B.append(phi_type)   # node 1
                leg_C.append(phip_type)  # node 8
            elif d == 3:
                leg_B.append(phi_type)   # node 2
                leg_C.append(phip_type)  # node 7
            elif d == 4:
                branch.append(phi_type)  # node 3
                leg_C.append(phip_type)  # node 6

    # dim 5, 6: unique
    if 5 in by_dim:
        leg_A.extend(by_dim[5])
    if 6 in by_dim:
        leg_A.extend(by_dim[6])

    print("Leg B (nodes 0,1,2, dims 1,2,3): irreps %s" % leg_B)
    print("Leg C (nodes 8,7,6, dims 2,3,4): irreps %s" % leg_C)
    print("Leg A (nodes 4,5, dims 5,6): irreps %s" % leg_A)
    print("Branch (node 3, dim 4): irreps %s" % branch)

    # 3x3 matrix: Leg B irreps vs 3 generation (2D) irreps
    if len(leg_B) == 3 and len(two_dim_indices) == 3:
        M_B = np.array([[branching_int[i, j] for j in two_dim_indices] for i in leg_B], dtype=float)
        print("\nLeg B 3x3 matrix (dims 1,2,3 vs 3 generations):")
        for row in range(3):
            print("  rho_%d (d=%d): [%s]"
                  % (leg_B[row], dims_sorted[leg_B[row]],
                     "  ".join("%d" % int(M_B[row, c]) for c in range(3))))

    if len(leg_C) == 3 and len(two_dim_indices) == 3:
        M_C = np.array([[branching_int[i, j] for j in two_dim_indices] for i in leg_C], dtype=float)
        print("\nLeg C 3x3 matrix (dims 2,3,4 vs 3 generations):")
        for row in range(3):
            print("  rho_%d (d=%d): [%s]"
                  % (leg_C[row], dims_sorted[leg_C[row]],
                     "  ".join("%d" % int(M_C[row, c]) for c in range(3))))

    # Analysis of 3x3 matrices
    print("\n--- Matrix Analysis ---")
    for label, M, leg in [("Leg B", M_B if len(leg_B)==3 else None, leg_B),
                           ("Leg C", M_C if len(leg_C)==3 else None, leg_C)]:
        if M is None:
            continue
        print("\n%s matrix:" % label)
        print("  M = ")
        for r in range(3):
            print("    [%s]" % "  ".join("%.0f" % M[r,c] for c in range(3)))

        # Check if it has CKM-like structure
        # CKM is nearly diagonal with small off-diagonal elements
        # But branching coefficients are integers, so structure is different
        total_per_row = [sum(M[r, :]) for r in range(3)]
        total_per_col = [sum(M[:, c]) for c in range(3)]
        print("  Row sums: %s" % total_per_row)
        print("  Col sums: %s" % total_per_col)

        # Normalize rows to get probabilities
        M_norm = np.zeros_like(M)
        for r in range(3):
            s = sum(M[r, :])
            if s > 0:
                M_norm[r, :] = M[r, :] / s
        print("  Normalized (row-stochastic):")
        for r in range(3):
            print("    [%s]" % "  ".join("%.4f" % M_norm[r,c] for c in range(3)))

        # Compare to PMNS |U|^2 (row-stochastic)
        s12_sq = 0.307
        s23_sq = 0.546
        s13_sq = 0.0220
        c12_sq = 1 - s12_sq
        c23_sq = 1 - s23_sq
        c13_sq = 1 - s13_sq
        PMNS_sq = np.array([
            [c12_sq*c13_sq, s12_sq*c13_sq, s13_sq],
            [0.327, 0.353, 0.320],
            [0.174, 0.340, 0.487]
        ])
        # (simplified PMNS)
        print("  PMNS |U|^2 for comparison:")
        for r in range(3):
            print("    [%s]" % "  ".join("%.4f" % PMNS_sq[r,c] for c in range(3)))

# ======================================================================
# PART 8: Coset Structure and Generation Assignment
# ======================================================================
print("\n--- PART 8: 5 Cosets and Generation Content ---")

# The 5 cosets of 2I/2T
print("Computing 5 cosets of 2I/2T...")
coset_reps = []
assigned = set(subgroup_2T)
cosets = [sorted(subgroup_2T)]

while len(assigned) < N:
    for g in range(N):
        if g not in assigned:
            new_coset = set()
            for h in subgroup_2T:
                new_coset.add(mult_table[g, h])
            cosets.append(sorted(new_coset))
            assigned |= new_coset
            coset_reps.append(g)
            break

print("Found %d cosets of size %d each" % (len(cosets), len(cosets[0])))

# For each coset, compute which 2I conjugacy class its elements belong to
print("\nCoset content by conjugacy class:")
for ci, coset in enumerate(cosets):
    class_count = [0] * n_classes
    for g in coset:
        class_count[class_of[g]] += 1
    nonzero = [(k, class_count[k]) for k in range(n_classes) if class_count[k] > 0]
    print("  Coset %d: %s" % (ci, nonzero))

# ======================================================================
# PART 9: Key Structural Tests
# ======================================================================
print("\n--- PART 9: Key Tests ---")

# Test 1: Are the 3 two-dim irreps of 2T related by outer automorphisms?
print("\nTest 1: 2D irrep relationships")
if len(two_dim_indices) == 3:
    for i_idx, j_2T in enumerate(two_dim_indices):
        print("  psi_%d characters: [" % j_2T, end="")
        for k in range(n_classes_2T):
            val = char_table_2T_sorted[j_2T, k]
            if abs(val.imag) < 1e-4:
                print("%.3f " % val.real, end="")
            else:
                print("%.2f%+.2fi " % (val.real, val.imag), end="")
        print("]")

    # Check if any pair of 2D irreps are complex conjugates
    for i_a in range(3):
        for i_b in range(i_a+1, 3):
            j_a = two_dim_indices[i_a]
            j_b = two_dim_indices[i_b]
            diff_real = max(abs(char_table_2T_sorted[j_a, k].real - char_table_2T_sorted[j_b, k].real)
                           for k in range(n_classes_2T))
            diff_conj = max(abs(char_table_2T_sorted[j_a, k] - np.conj(char_table_2T_sorted[j_b, k]))
                           for k in range(n_classes_2T))
            print("  psi_%d vs psi_%d: max|diff_real|=%.4f, max|diff_conj|=%.4f"
                  % (j_a, j_b, diff_real, diff_conj))

# Test 2: Total generation content per E8 leg
print("\nTest 2: Total 2D content by E8 leg")
for label, leg in [("Leg A", leg_A), ("Leg B", leg_B), ("Leg C", leg_C), ("Branch", branch)]:
    total = [0, 0, 0]
    for i in leg:
        for idx_g, j in enumerate(two_dim_indices):
            total[idx_g] += branching_int[i, j] * dims_sorted[i]  # weighted by dim
    print("  %s: gen content (dim-weighted) = %s, total = %d" % (label, total, sum(total)))

# Test 3: Does the branching matrix have any phi structure?
print("\nTest 3: Phi structure in branching")
for i in range(len(dims_sorted)):
    gen_content = [branching_int[i, j] for j in two_dim_indices]
    if max(gen_content) > 0:
        # Try ratios
        nz = [g for g in gen_content if g > 0]
        if len(nz) >= 2:
            ratios = []
            for a in range(len(nz)):
                for b in range(a+1, len(nz)):
                    r = nz[a] / nz[b]
                    ratios.append(r)
            ratio_str = ", ".join("%.4f" % r for r in ratios)
        else:
            ratio_str = "N/A"
        print("  rho_%d (d=%d): gen = %s, ratios = %s"
              % (i, dims_sorted[i], gen_content, ratio_str))

# ======================================================================
# SUMMARY
# ======================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("\n1. 2I has %d conjugacy classes (should be 9)" % n_classes)
print("   Class sizes: %s" % class_sizes_2I)
print("   Irrep dims: %s" % dims_sorted)
print("\n2. 2T has %d conjugacy classes (should be 7)" % n_classes_2T)
print("   Class sizes: %s" % class_sizes_2T)
print("   Irrep dims: %s" % dims_2T_sorted)
print("\n3. Three 2D irreps of 2T: %s" % ("YES" if len(two_dim_indices) == 3 else "NO"))
print("   -> Potential 3-generation structure")
print("\n4. Branching rules computed: %dx%d matrix"
      % (len(dims_sorted), len(dims_2T_sorted)))
print("\n5. Key question: Do the 3x3 sub-matrices encode CKM/PMNS?")
print("   -> Branching coefficients are INTEGERS, not continuous")
print("   -> They constrain the MULTIPLICITY SPACE structure")
print("   -> Mixing angles arise from basis choice WITHIN multiplicity spaces")

print("\n" + "=" * 70)
print("EXP-215 COMPLETE")
print("=" * 70)

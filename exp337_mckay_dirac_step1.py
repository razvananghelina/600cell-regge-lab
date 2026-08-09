"""
exp337 Step 1: Build 2I group, character table, verify McKay graph
==================================================================
Binary icosahedral group 2I (order 120) embedded in SU(2).
The 120 elements are the vertices of the 600-cell as unit quaternions.
"""
import numpy as np
from itertools import product

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2  # = -1/PHI

# ================================================================
# Step 1: Generate all 120 elements of 2I as SU(2) matrices
# ================================================================
# Quaternion q = (q0, q1, q2, q3) -> SU(2) matrix:
# [[q0 + i*q1, -q2 + i*q3],
#  [q2 + i*q3,  q0 - i*q1]]

def quat_to_su2(q):
    """Convert unit quaternion (q0,q1,q2,q3) to SU(2) matrix."""
    a, b, c, d = q
    return np.array([[a + 1j*b, -c + 1j*d],
                      [c + 1j*d,  a - 1j*b]])

def generate_2I():
    """Generate all 120 elements of binary icosahedral group as SU(2) matrices."""
    elements = []
    quats = []

    # Part 1: 8 elements - permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [1, -1]:
            q = [0, 0, 0, 0]
            q[i] = s
            quats.append(tuple(q))

    # Part 2: 16 elements - (+-1/2, +-1/2, +-1/2, +-1/2)
    for signs in product([0.5, -0.5], repeat=4):
        quats.append(signs)

    # Part 3: 96 elements - even permutations of (0, +-1/(2*phi), +-1/2, +-phi/2)
    vals = [0, 1/(2*PHI), 0.5, PHI/2]

    # The 12 even permutations of (0,1,2,3)
    even_perms = [
        (0,1,2,3), (0,2,3,1), (0,3,1,2),
        (1,0,3,2), (1,2,0,3), (1,3,2,0),
        (2,0,1,3), (2,1,3,0), (2,3,0,1),
        (3,0,2,1), (3,1,0,2), (3,2,1,0),
    ]

    for perm in even_perms:
        base = [vals[perm[0]], vals[perm[1]], vals[perm[2]], vals[perm[3]]]
        # Find which positions are nonzero
        nonzero_pos = [i for i in range(4) if abs(base[i]) > 1e-10]
        # Apply all sign combinations to nonzero positions
        for signs in product([1, -1], repeat=len(nonzero_pos)):
            q = list(base)
            for idx, s in zip(nonzero_pos, signs):
                q[idx] = abs(q[idx]) * s
            quats.append(tuple(q))

    # Remove duplicates (within tolerance)
    unique_quats = []
    for q in quats:
        is_dup = False
        for uq in unique_quats:
            if np.allclose(q, uq, atol=1e-10):
                is_dup = True
                break
        if not is_dup:
            unique_quats.append(q)

    print(f"Generated {len(unique_quats)} unique quaternions (expect 120)")

    # Convert to SU(2) matrices
    matrices = [quat_to_su2(q) for q in unique_quats]

    # Verify: each should be unitary with det=1
    for i, M in enumerate(matrices):
        assert np.allclose(M @ M.conj().T, np.eye(2), atol=1e-10), f"Not unitary: {i}"
        assert np.allclose(np.linalg.det(M), 1.0, atol=1e-10), f"Det != 1: {i}"

    # Verify closure under multiplication
    print("Verifying group closure...")
    for i in range(min(20, len(matrices))):  # spot-check
        for j in range(min(20, len(matrices))):
            prod = matrices[i] @ matrices[j]
            found = False
            for k, Mk in enumerate(matrices):
                if np.allclose(prod, Mk, atol=1e-8):
                    found = True
                    break
            assert found, f"Product {i}*{j} not in group!"

    print("Group closure verified (spot-check).")
    return matrices, unique_quats


# ================================================================
# Step 2: Identify conjugacy classes
# ================================================================
def find_conjugacy_classes(matrices):
    """Identify conjugacy classes by trace (= 2*cos(alpha))."""
    traces = [np.trace(M).real for M in matrices]

    # Group by trace value
    class_traces = []
    class_indices = []
    for i, tr in enumerate(traces):
        found = False
        for j, ct in enumerate(class_traces):
            if abs(tr - ct) < 1e-8:
                class_indices[j].append(i)
                found = True
                break
        if not found:
            class_traces.append(tr)
            class_indices.append([i])

    # Sort by trace (descending)
    order = np.argsort(class_traces)[::-1]
    class_traces = [class_traces[i] for i in order]
    class_indices = [class_indices[i] for i in order]

    print(f"\nFound {len(class_traces)} conjugacy classes:")
    print(f"{'Class':>6} {'Size':>5} {'Trace':>10} {'cos(a)':>10} {'alpha':>10} {'Order':>6}")
    print("-"*55)
    for i, (tr, idx) in enumerate(zip(class_traces, class_indices)):
        cos_a = tr / 2
        if abs(cos_a) <= 1:
            alpha = np.arccos(np.clip(cos_a, -1, 1))
        else:
            alpha = 0 if cos_a > 0 else np.pi
        # Element order: smallest n such that M^n = +-I
        M = matrices[idx[0]]
        Mn = np.eye(2, dtype=complex)
        for n in range(1, 61):
            Mn = Mn @ M
            if np.allclose(Mn, np.eye(2), atol=1e-8):
                order = n
                break
            if np.allclose(Mn, -np.eye(2), atol=1e-8):
                order = 2*n
                break
        print(f"C_{i:>3} {len(idx):>5} {tr:>10.5f} {cos_a:>10.5f} {alpha:>10.5f} {order:>6}")

    return class_traces, class_indices


# ================================================================
# Step 3: Compute symmetric power representations
# ================================================================
def sym_power_matrix(M, k):
    """Compute the (k+1)x(k+1) matrix of Sym^k(M) for M in SU(2).
    Uses monomial basis: v_m = e1^(k-m) * e2^m, m=0,...,k
    With normalization: v_m_norm = v_m / sqrt(C(k,m))
    """
    if k == 0:
        return np.array([[1.0 + 0j]])
    if k == 1:
        return M.copy()

    a, b = M[0, 0], M[0, 1]
    c, d = M[1, 0], M[1, 1]

    from math import comb

    dim = k + 1
    R = np.zeros((dim, dim), dtype=complex)

    for m in range(dim):  # column index (input basis)
        for p in range(dim):  # row index (output basis)
            # R_pm = sum over i,j with i+j=p, 0<=i<=k-m, 0<=j<=m
            val = 0
            for i in range(max(0, p-m), min(k-m, p)+1):
                j = p - i
                if j < 0 or j > m:
                    continue
                val += (comb(k-m, i) * comb(m, j) *
                        a**(k-m-i) * c**i * b**(m-j) * d**j)
            # Normalize for orthonormal basis
            val *= np.sqrt(comb(k, p)) / np.sqrt(comb(k, m))
            R[p, m] = val

    return R


# ================================================================
# Step 4: Compute character table
# ================================================================
def compute_characters(matrices, class_traces, class_indices, max_k=11):
    """Compute characters of Sym^k for k=0,...,max_k."""
    n_classes = len(class_traces)
    chars = np.zeros((max_k+1, n_classes))

    for k in range(max_k+1):
        for ci in range(n_classes):
            rep_idx = class_indices[ci][0]  # representative
            M = matrices[rep_idx]
            Sk = sym_power_matrix(M, k)
            chars[k, ci] = np.trace(Sk).real

    return chars


# ================================================================
# MAIN
# ================================================================
print("="*72)
print("exp337 Step 1: Binary Icosahedral Group 2I")
print("="*72)

matrices, quats = generate_2I()
class_traces, class_indices = find_conjugacy_classes(matrices)
class_sizes = [len(idx) for idx in class_indices]

print("\n" + "="*72)
print("CHARACTER TABLE OF Sym^k REPRESENTATIONS")
print("="*72)

chars = compute_characters(matrices, class_traces, class_indices, max_k=11)

print(f"\n{'k':>3} {'dim':>4} |", end="")
for i in range(len(class_traces)):
    print(f" C{i:>1}({class_sizes[i]:>2})", end="")
print()
print("-"*80)
for k in range(12):
    print(f"{k:>3} {k+1:>4} |", end="")
    for ci in range(len(class_traces)):
        print(f" {chars[k, ci]:>6.3f}", end="")
    print()

# ================================================================
# Step 5: Decompose Sym^k into irreps using McKay recursion
# ================================================================
print("\n" + "="*72)
print("IRREP DECOMPOSITION via McKay recursion")
print("="*72)

# rho_0 through rho_5 = Sym^0 through Sym^5
# Sym^6 = rho_6 + rho_8 (from rho_1 x rho_5 = Sym^6 + Sym^4 = rho_4 + rho_6 + rho_8)
# Sym^7 = rho_5 + rho_7 (from rho_1 x Sym^6 = Sym^7 + Sym^5, etc.)

n_classes = len(class_traces)
irrep_chars = np.zeros((9, n_classes))

# rho_0 to rho_5 = Sym^0 to Sym^5
for k in range(6):
    irrep_chars[k] = chars[k]

# rho_6 + rho_8 = Sym^6
# rho_5 + rho_7 = Sym^7
# From McKay: rho_1 x rho_8 = rho_5, so chi_1 * chi_8 = chi_5
# This gives chi_8 = chi_5 / chi_1 where defined (careful with zeros)

# Use character inner products to decompose Sym^6
# First compute chi_8 from chi_1 * chi_8 = chi_5
chi_8 = np.zeros(n_classes)
for ci in range(n_classes):
    if abs(chars[1, ci]) > 1e-10:
        chi_8[ci] = chars[5, ci] / chars[1, ci]
    else:
        # At C4 class (trace=0), chi_1=0, chi_5=0, need another method
        # Use norm constraint
        pass

# For the class where chi_1=0: use <chi_8, chi_8> = 1
# Sum without the zero-trace class:
partial_sum = sum(chi_8[ci]**2 * class_sizes[ci]
                  for ci in range(n_classes) if abs(chars[1, ci]) > 1e-10)
# Find the index of C4 class (trace ~ 0)
c4_idx = None
for ci in range(n_classes):
    if abs(class_traces[ci]) < 1e-8 and class_sizes[ci] == 30:
        c4_idx = ci
        break

if c4_idx is not None:
    # <chi_8, chi_8> = 1 => partial_sum + chi_8[c4]^2 * 30 = 120
    remaining = 120 - partial_sum
    chi_8_c4_sq = remaining / 30
    print(f"chi_8 at C4: x^2 = {chi_8_c4_sq:.4f}, so x = +-{np.sqrt(abs(chi_8_c4_sq)):.4f}")
    # Use orthogonality with chi_2 to determine sign
    # <chi_8, chi_2> should be 0
    partial_orth = sum(chi_8[ci] * irrep_chars[2, ci] * class_sizes[ci]
                       for ci in range(n_classes) if ci != c4_idx)
    # partial_orth + chi_8[c4] * chi_2[c4] * 30 = 0
    chi_2_c4 = irrep_chars[2, c4_idx]
    if abs(chi_2_c4) > 1e-10:
        chi_8[c4_idx] = -partial_orth / (chi_2_c4 * 30)
    else:
        chi_8[c4_idx] = np.sqrt(abs(chi_8_c4_sq))

irrep_chars[8] = chi_8

# rho_6 = Sym^6 - rho_8
irrep_chars[6] = chars[6] - chi_8

# rho_7 = Sym^7 - rho_5
irrep_chars[7] = chars[7] - irrep_chars[5]

# Print character table
print("\nFull character table of 2I:")
names = ['rho_0', 'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5', 'rho_6', 'rho_7', 'rho_8']
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

print(f"{'':>8} {'dim':>4} |", end="")
for i in range(n_classes):
    print(f" C{i}({class_sizes[i]:>2})", end="")
print()
print("-"*90)
for i in range(9):
    print(f"{names[i]:>8} {dims[i]:>4} |", end="")
    for ci in range(n_classes):
        v = irrep_chars[i, ci]
        # Try to identify as phi, phi', integer, etc.
        print(f" {v:>7.3f}", end="")
    print()

# Verify orthogonality
print("\nOrthogonality check <chi_i, chi_j>:")
gram = np.zeros((9, 9))
for i in range(9):
    for j in range(9):
        gram[i, j] = sum(irrep_chars[i, ci] * irrep_chars[j, ci] * class_sizes[ci]
                         for ci in range(n_classes)) / 120
print("Should be identity matrix:")
print(np.round(gram, 4))

# ================================================================
# Step 6: Verify McKay graph
# ================================================================
print("\n" + "="*72)
print("VERIFY McKay GRAPH: rho_1 x rho_i = sum of neighbors")
print("="*72)

for i in range(9):
    # Compute chi_1 * chi_i
    product_char = irrep_chars[1] * irrep_chars[i]
    # Decompose into irreps
    print(f"\nrho_1 x rho_{i} (dim {dims[i]}) = ", end="")
    decomp = []
    for j in range(9):
        mult = sum(product_char[ci] * irrep_chars[j, ci] * class_sizes[ci]
                   for ci in range(n_classes)) / 120
        mult_int = int(round(mult))
        if mult_int > 0:
            decomp.append(f"{mult_int}*rho_{j}" if mult_int > 1 else f"rho_{j}")
    print(" + ".join(decomp))

# The McKay graph adjacency:
print("\nExpected McKay graph (affine E8):")
print("rho_0 - rho_1 - rho_2 - rho_3 - rho_4 - rho_5 - rho_8")
print("                                        |")
print("                                       rho_6 - rho_7")
print("\n(Branch at rho_5 (dim 6 = b1), legs (1, 2, 5))")

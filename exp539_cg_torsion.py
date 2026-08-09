"""
EXP-539: CG-Weighted Spectral Torsion on the McKay Graph
==========================================================

Build D_F with PHYSICAL Clebsch-Gordan coefficients of 2I,
not democratic Yukawa. Then recompute spectral torsion.

STEPS:
  1. Build all 9 irreps of 2I as explicit matrices (from quaternions)
  2. Compute CG intertwiners for each McKay edge: rho_1 x rho_k -> rho_l
  3. Build D_F with CG-weighted Yukawa
  4. Compute torsion, Ricci scalar, Tr(D^4)
  5. Compare with democratic Yukawa (exp537)
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons.cell600 import build_600cell
from commons import (PHI, a1, b1, N, rank_E8, IRREP_DIMS_2I, MCKAY_EDGES)

print("=" * 70)
print("EXP-539: CG-WEIGHTED SPECTRAL TORSION")
print("=" * 70)

# ============================================================
# SECTION 1: BUILD 2I ELEMENTS AS SU(2) MATRICES
# ============================================================
print(f"\n{'='*70}")
print("SECTION 1: 2I AS SU(2) MATRICES")
print(f"{'='*70}")

verts, _, _ = build_600cell()
n_group = len(verts)  # = 120

# Quaternion (a, b, c, d) -> SU(2) matrix:
# [[a + bi, c + di], [-c + di, a - bi]]
def quat_to_su2(q):
    a, b, c, d = q
    return np.array([[a + 1j*b, c + 1j*d],
                     [-c + 1j*d, a - 1j*b]])

group_su2 = np.array([quat_to_su2(v) for v in verts])  # (120, 2, 2)

# Verify: all should be unitary with det = 1
for i in range(n_group):
    g = group_su2[i]
    assert abs(np.linalg.det(g) - 1.0) < 1e-8, f"det(g_{i}) = {np.linalg.det(g)}"
    assert np.allclose(g @ g.conj().T, np.eye(2), atol=1e-8), f"g_{i} not unitary"

# Find identity element
id_idx = -1
for i in range(n_group):
    if np.allclose(group_su2[i], np.eye(2), atol=1e-8):
        id_idx = i
        break
print(f"  Identity at index: {id_idx}")

# Build multiplication table (for group algebra operations)
# g_i * g_j = g_{mult[i,j]}
print(f"  Building multiplication table...")
mult_table = np.zeros((n_group, n_group), dtype=int)
for i in range(n_group):
    for j in range(n_group):
        prod = group_su2[i] @ group_su2[j]
        # Find which element this is
        dists = np.array([np.linalg.norm(prod - group_su2[k]) for k in range(n_group)])
        mult_table[i, j] = np.argmin(dists)
        assert dists[mult_table[i, j]] < 1e-6, f"Product g_{i}*g_{j} not found in group"

print(f"  Multiplication table built: {n_group}x{n_group}")

# Find inverse table
inv_table = np.zeros(n_group, dtype=int)
for i in range(n_group):
    inv_table[i] = mult_table[i].tolist().index(id_idx) if id_idx in mult_table[i] else -1
    # Actually: g_i * g_j = e => j = inv(i)
    for j in range(n_group):
        if mult_table[i, j] == id_idx:
            inv_table[i] = j
            break

# ============================================================
# SECTION 2: BUILD ALL 9 IRREPS
# ============================================================
print(f"\n{'='*70}")
print("SECTION 2: ALL 9 IRREPS OF 2I")
print(f"{'='*70}")

dims = IRREP_DIMS_2I  # [1, 2, 3, 4, 5, 6, 4, 2, 3]

# rho_0 (dim 1): trivial representation
# rho_1 (dim 2): fundamental = the SU(2) matrices themselves
# rho_k (dim d): symmetric powers Sym^k(rho_1) for k=0,...,3 give dims 1,2,3,4
# But 2I has 9 irreps, and Sym^k for k>3 may decompose.

# Method: build Sym^k(rho_1) for k=0,...,8, then use character inner products
# to extract the 9 irreps.

# Symmetric power: Sym^n(V) has basis {e_{i1} ... e_{in}} with i1<=...<=in
# For V = C^2: Sym^n has dim n+1
# The representation on Sym^n is: rho_n(g) acts on the symmetric tensor.

def sym_power_rep(g_2x2, n):
    """Compute the (n+1)x(n+1) matrix of g acting on Sym^n(C^2).
    Basis: e1^n, e1^{n-1}*e2, ..., e2^n (with binomial normalization)."""
    dim = n + 1
    mat = np.zeros((dim, dim), dtype=complex)
    a, b = g_2x2[0, 0], g_2x2[0, 1]
    c, d = g_2x2[1, 0], g_2x2[1, 1]

    for i in range(dim):  # output: e1^{n-i} * e2^i
        for j in range(dim):  # input: e1^{n-j} * e2^j
            # Coefficient of e1^{n-i}*e2^i in g(e1^{n-j}*e2^j)
            # g(e1) = a*e1 + c*e2,  g(e2) = b*e1 + d*e2
            # g(e1^{n-j} * e2^j) = (a*e1+c*e2)^{n-j} * (b*e1+d*e2)^j
            # Coefficient of e1^{n-i}*e2^i = sum over partitions
            coeff = 0
            for p in range(max(0, n-j-i+j), min(n-j, n-i)+1):
                # p = number of e1 from first factor
                # q = n-i-p = number of e1 from second factor
                q = n - i - p
                if q < 0 or q > j:
                    continue
                from math import comb
                coeff += (comb(n-j, p) * a**p * c**(n-j-p) *
                         comb(j, q) * b**q * d**(j-q))
            # Normalize by sqrt of binomial coefficients
            from math import comb as comb_func
            norm_factor = np.sqrt(comb_func(n, i) / comb_func(n, j))
            mat[i, j] = coeff * norm_factor
    return mat


# Build character table using Sym^k representations
print(f"  Computing characters via Sym^k(rho_1)...")

# First: find conjugacy classes
# Two elements are conjugate if they have the same trace in rho_1
traces_rho1 = np.array([np.trace(group_su2[i]) for i in range(n_group)])
# Round to cluster
unique_traces = []
class_labels = np.zeros(n_group, dtype=int)
for i in range(n_group):
    found = False
    for ci, ut in enumerate(unique_traces):
        if abs(traces_rho1[i] - ut) < 1e-6:
            class_labels[i] = ci
            found = True
            break
    if not found:
        class_labels[i] = len(unique_traces)
        unique_traces.append(traces_rho1[i])

n_classes = len(unique_traces)
class_sizes = [np.sum(class_labels == c) for c in range(n_classes)]
print(f"  Conjugacy classes: {n_classes}")
print(f"  Class sizes: {class_sizes}")
print(f"  Sum: {sum(class_sizes)} (should be 120)")

# Representatives: first element of each class
class_reps = [np.where(class_labels == c)[0][0] for c in range(n_classes)]

# Compute characters of Sym^k for k=0,...,8
max_k = 8
sym_chars = np.zeros((max_k + 1, n_classes), dtype=complex)
for k in range(max_k + 1):
    for c in range(n_classes):
        g = group_su2[class_reps[c]]
        rho_k = sym_power_rep(g, k)
        sym_chars[k, c] = np.trace(rho_k)

# The irreps of 2I are NOT all Sym^k. For k >= 4, Sym^k decomposes.
# Use character orthogonality to extract irreducible characters.
# <chi_i, chi_j> = (1/|G|) * sum_g chi_i(g)* chi_j(g) = delta_ij

# Start with known irreps: Sym^0 (dim 1), Sym^1 (dim 2), Sym^2 (dim 3), Sym^3 (dim 4)
# These should be irreducible for 2I.
irrep_chars = []
irrep_names = []

for k in range(4):  # Sym^0, Sym^1, Sym^2, Sym^3 are irreducible
    chi = sym_chars[k].real  # characters are real for 2I
    # Verify irreducibility: <chi, chi> = 1
    inner = sum(class_sizes[c] * chi[c]**2 for c in range(n_classes)) / n_group
    irrep_chars.append(chi.copy())
    irrep_names.append(f"Sym^{k} (dim {k+1})")
    print(f"  Sym^{k} (dim {k+1}): <chi,chi> = {inner:.4f} {'IRR' if abs(inner-1)<0.1 else 'RED'}")

# Sym^4 (dim 5) should be irreducible (= rho_4, standard rep of A5)
chi_4 = sym_chars[4].real
inner_4 = sum(class_sizes[c] * chi_4[c]**2 for c in range(n_classes)) / n_group
print(f"  Sym^4 (dim 5): <chi,chi> = {inner_4:.4f} {'IRR' if abs(inner_4-1)<0.1 else 'RED'}")
irrep_chars.append(chi_4.copy())
irrep_names.append("Sym^4 (dim 5)")

# Sym^5 (dim 6): should be irreducible (= rho_5)
chi_5 = sym_chars[5].real
inner_5 = sum(class_sizes[c] * chi_5[c]**2 for c in range(n_classes)) / n_group
print(f"  Sym^5 (dim 6): <chi,chi> = {inner_5:.4f} {'IRR' if abs(inner_5-1)<0.1 else 'RED'}")

if abs(inner_5 - 1) < 0.1:
    irrep_chars.append(chi_5.copy())
    irrep_names.append("Sym^5 (dim 6)")
else:
    # Sym^5 is reducible. Decompose: Sym^5 = rho_5 + ...
    # Remove known components
    chi_5_reduced = chi_5.copy()
    for known_chi in irrep_chars:
        inner = sum(class_sizes[c] * chi_5_reduced[c] * known_chi[c] for c in range(n_classes)) / n_group
        if abs(inner) > 0.1:
            chi_5_reduced -= round(inner) * known_chi
    norm = sum(class_sizes[c] * chi_5_reduced[c]**2 for c in range(n_classes)) / n_group
    if abs(norm - 1) < 0.1:
        irrep_chars.append(chi_5_reduced)
        irrep_names.append(f"Sym^5 reduced (dim {int(chi_5_reduced[class_reps.index(class_reps[0])])})")
    # Also add Sym^5 itself if irreducible
    print(f"    Sym^5 reduced: <chi,chi> = {norm:.4f}")

# For remaining irreps: use Sym^6, Sym^7, Sym^8 and subtract known ones
for k in range(6, max_k + 1):
    chi_k = sym_chars[k].real
    # Subtract known irreps
    chi_new = chi_k.copy()
    for known_chi in irrep_chars:
        inner = sum(class_sizes[c] * chi_new[c] * known_chi[c] for c in range(n_classes)) / n_group
        mult = round(inner.real)
        if mult != 0:
            chi_new -= mult * known_chi

    norm = sum(class_sizes[c] * chi_new[c]**2 for c in range(n_classes)) / n_group
    dim_new = int(round(chi_new[0].real)) if abs(chi_new[0]) > 0.5 else 0

    if abs(norm - 1) < 0.1 and dim_new > 0:
        # Check if we already have this irrep
        is_new = True
        for known_chi in irrep_chars:
            overlap = sum(class_sizes[c] * chi_new[c] * known_chi[c] for c in range(n_classes)) / n_group
            if abs(overlap - 1) < 0.1:
                is_new = False
                break
        if is_new:
            irrep_chars.append(chi_new)
            irrep_names.append(f"from Sym^{k} (dim {dim_new})")
            print(f"  From Sym^{k}: new irrep dim {dim_new}, <chi,chi> = {norm:.4f}")
    elif norm > 1.1:
        # Multiple new irreps in here — try to split further
        print(f"  Sym^{k} residual: <chi,chi> = {norm:.4f} (reducible, dim {dim_new})")

print(f"\n  Found {len(irrep_chars)} irreps so far (need 9)")
print(f"  Dims: {[int(round(chi[0].real)) for chi in irrep_chars]}")

# ============================================================
# SECTION 3: CHARACTER TABLE AND CG DECOMPOSITION
# ============================================================
print(f"\n{'='*70}")
print("SECTION 3: CHARACTER TABLE")
print(f"{'='*70}")

print(f"\n  {'Irrep':>25s}", end="")
for c in range(n_classes):
    print(f" {'C'+str(c):>6s}", end="")
print()
print(f"  {'(class size)':>25s}", end="")
for c in range(n_classes):
    print(f" {class_sizes[c]:>6d}", end="")
print()

for idx, (chi, name) in enumerate(zip(irrep_chars, irrep_names)):
    print(f"  {name:>25s}", end="")
    for c in range(n_classes):
        print(f" {chi[c]:>6.2f}", end="")
    print()

# Verify orthogonality
print(f"\n  Orthogonality check:")
for i in range(min(len(irrep_chars), 6)):
    for j in range(i, min(len(irrep_chars), 6)):
        inner = sum(class_sizes[c] * irrep_chars[i][c] * irrep_chars[j][c]
                    for c in range(n_classes)) / n_group
        if abs(inner) > 0.01:
            print(f"    <chi_{i}, chi_{j}> = {inner:.4f}")

# ============================================================
# SECTION 4: CG COEFFICIENTS FOR McKAY EDGES
# ============================================================
print(f"\n{'='*70}")
print("SECTION 4: CG COEFFICIENTS (rho_1 x rho_k -> rho_l)")
print(f"{'='*70}")

# For each McKay edge (k,l): rho_l appears in rho_1 x rho_k
# The CG coefficient is the intertwiner I: V_1 x V_k -> V_l

# Character of rho_1 x rho_k:
# chi_{1 x k}(g) = chi_1(g) * chi_k(g)

if len(irrep_chars) >= 6:
    print(f"  Decomposition of rho_1 x rho_k:")
    for k in range(min(len(irrep_chars), 6)):
        # chi_{1 x k}
        chi_product = np.array([irrep_chars[1][c] * irrep_chars[k][c] for c in range(n_classes)])
        # Decompose into irreps
        decomp = []
        for l in range(len(irrep_chars)):
            mult = sum(class_sizes[c] * chi_product[c] * irrep_chars[l][c]
                       for c in range(n_classes)) / n_group
            if abs(mult) > 0.5:
                decomp.append((l, int(round(mult.real))))
        dim_k = int(round(irrep_chars[k][0].real))
        print(f"    rho_1 x rho_{k} (dim 2x{dim_k}={2*dim_k}) = "
              f"{' + '.join(f'{m}*rho_{l}' if m>1 else f'rho_{l}' for l,m in decomp)}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

n_found = len(irrep_chars)
dims_found = sorted([int(round(chi[0].real)) for chi in irrep_chars])
print(f"""
  Found {n_found}/9 irreps of 2I via symmetric powers.
  Dimensions found: {dims_found}
  Expected: [1, 2, 2, 3, 3, 4, 4, 5, 6]

  {'COMPLETE' if n_found == 9 else f'INCOMPLETE — need {9-n_found} more irreps'}

  Next steps:
  - If 9 irreps found: build explicit matrix reps, compute CG intertwiners
  - If incomplete: use alternative method (induced representations or
    direct character table lookup from literature)
""")

print("=" * 70)
print("EXP-539 COMPLETE")
print("=" * 70)

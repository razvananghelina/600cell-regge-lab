"""
EXP-540: CG-Weighted Torsion — Clean Computation
==================================================

Step-by-step, no shortcuts:
  Step 1: Hardcode character table of 2I (from literature)
  Step 2: Build group elements as SU(2) matrices
  Step 3: Build multiplication table
  Step 4: Build all 9 irreps via regular representation + projection
  Step 5: Compute CG intertwiners for each McKay edge
  Step 6: Build D_F with CG Yukawa
  Step 7: Compute spectral torsion
"""

import numpy as np
from math import comb
import sys
sys.path.insert(0, '.')
from commons.cell600 import build_600cell
from commons import (PHI, a1, b1, rank_E8, IRREP_DIMS_2I, MCKAY_EDGES)

print("=" * 70)
print("EXP-540: CG-WEIGHTED TORSION — CLEAN COMPUTATION")
print("=" * 70)

# ============================================================
# STEP 1: CHARACTER TABLE OF 2I (HARDCODED)
# ============================================================
print(f"\n{'='*70}")
print("STEP 1: CHARACTER TABLE OF 2I")
print(f"{'='*70}")

# Binary icosahedral group 2I, order 120.
# 9 conjugacy classes, 9 irreps.
#
# Conjugacy classes (by rotation angle theta of the SU(2) element):
#   C0: {+I},     size 1,  theta = 0
#   C1: {-I},     size 1,  theta = 2*pi (= pi for SU(2) half-angle)
#   C2: order 10, size 12, half-angle pi/5
#   C3: order 10, size 12, half-angle 2*pi/5
#   C4: order 6,  size 20, half-angle pi/3
#   C5: order 4,  size 30, half-angle pi/2
#   C6: order 6,  size 20, half-angle 2*pi/3
#   C7: order 5,  size 12, half-angle pi/5*2 = 2*pi/5
#   C8: order 5,  size 12, half-angle 3*pi/5
#
# For SU(2): chi_j(theta) = sin((2j+1)*theta/2) / sin(theta/2)
# where theta is the rotation angle (half-angle = theta/2 for eigenvalues).
#
# The 9 irreps have dimensions [1, 2, 3, 4, 5, 6, 4, 2, 3]
# and correspond to SU(2) spins j = 0, 1/2, 1, 3/2, 2, 5/2
# plus rho_6 (dim 4), rho_7 (dim 2), rho_8 (dim 3) from
# decomposition of higher spins when restricted to 2I.

# Character table from standard references (e.g., Conway & Sloane,
# or GAP computation). Using half-angles for SU(2) chi formula.

# Half-angles for each conjugacy class representative:
# C0: 0, C1: pi, C2: pi/5, C3: 2pi/5, C4: pi/3, C5: pi/2,
# C6: 2pi/3, C7: 2pi/5, C8: 3pi/5
#
# Wait — need to be careful. The conjugacy classes in 2I:
# In SU(2), element g has eigenvalues exp(+/- i*alpha) where alpha is the half-angle.
# trace(g) = 2*cos(alpha).
#
# Conjugacy classes of 2I by trace of fundamental rep:
# C0: tr = 2          (identity, size 1)
# C1: tr = -2         (-I, size 1)
# C2: tr = 2*cos(pi/5) = phi      (size 12, order 10)
# C3: tr = 2*cos(2pi/5) = 1/phi   (size 12, order 10)
# C4: tr = 2*cos(pi/3) = 1        (size 20, order 6)
# C5: tr = 2*cos(pi/2) = 0        (size 30, order 4)
# C6: tr = 2*cos(2pi/3) = -1      (size 20, order 6)
# C7: tr = 2*cos(3pi/5) = -1/phi  (size 12, order 5)
# C8: tr = 2*cos(4pi/5) = -phi    (size 12, order 5)

class_sizes = np.array([1, 1, 12, 12, 20, 30, 20, 12, 12])
n_classes = 9
n_group = 120

# Half-angles alpha for each class (eigenvalues of SU(2) element = exp(+/-i*alpha)):
half_angles = [0, np.pi, np.pi/5, 2*np.pi/5, np.pi/3, np.pi/2,
               2*np.pi/3, 3*np.pi/5, 4*np.pi/5]

# Character of SU(2) spin-j rep: chi_j(alpha) = sin((2j+1)*alpha) / sin(alpha)
# For alpha = 0: chi_j = 2j+1 (dimension)
# For alpha = pi: chi_j = (-1)^(2j) * (2j+1)... actually sin((2j+1)*pi)/sin(pi) = 0/0
# Use L'Hopital or direct: chi_j(pi) = (-1)^(2j) * (2j+1)

def su2_character(j, alpha):
    """Character of SU(2) spin-j representation at half-angle alpha."""
    if abs(alpha) < 1e-12:
        return 2*j + 1
    if abs(alpha - np.pi) < 1e-12:
        return (-1)**(int(2*j)) * (2*j + 1)
    return np.sin((2*j+1)*alpha) / np.sin(alpha)

# Build character table for SU(2) spins j = 0, 1/2, 1, ..., 7/2
# Then decompose restrictions to 2I.

# First 6 irreps of 2I = restrictions of j=0,1/2,1,3/2,2,5/2 (all irreducible on 2I)
# rho_0 (dim 1) = V_0
# rho_1 (dim 2) = V_{1/2}
# rho_2 (dim 3) = V_1
# rho_3 (dim 4) = V_{3/2}
# rho_4 (dim 5) = V_2
# rho_5 (dim 6) = V_{5/2}

# Remaining 3 irreps from V_3 (dim 7) and V_{7/2} (dim 8):
# V_3|_{2I} = rho_6 (dim 4) + rho_8 (dim 3)
# V_{7/2}|_{2I} = rho_7 (dim 2) + rho_5 (dim 6)
# => rho_7 = V_{7/2} - rho_5

# Build full character table
dims = np.array(IRREP_DIMS_2I)  # [1, 2, 3, 4, 5, 6, 4, 2, 3]
char_table = np.zeros((9, 9))

# rho_0 through rho_5: direct from SU(2) spin j = k/2
for k in range(6):
    j = k / 2.0
    for c in range(n_classes):
        char_table[k, c] = su2_character(j, half_angles[c])

# rho_6 (dim 4) and rho_8 (dim 3): from V_3 = rho_6 + rho_8
# chi_{V_3} = chi_{rho_6} + chi_{rho_8}
chi_V3 = np.array([su2_character(3, a) for a in half_angles])

# rho_7 (dim 2): from V_{7/2} = rho_7 + rho_5
# chi_{V_{7/2}} = chi_{rho_7} + chi_{rho_5}
chi_V72 = np.array([su2_character(3.5, a) for a in half_angles])
char_table[7] = chi_V72 - char_table[5]  # rho_7 = V_{7/2} - rho_5

# Now: rho_6 + rho_8 = V_3. We need another equation.
# V_{9/2}|_{2I} decomposes too. Let's use V_4 (dim 9):
# V_4|_{2I} should decompose into known irreps.
# From McKay: tensor products give us the decomposition.
# rho_1 x rho_5 = rho_4 + rho_6 + rho_8
# chi_{rho_1 x rho_5} = chi_1 * chi_5
# So chi_{rho_6} + chi_{rho_8} = chi_1 * chi_5 - chi_4

chi_sum_68 = char_table[1] * char_table[5] - char_table[4]

# Also chi_V3 = chi_{rho_6} + chi_{rho_8} (from SU(2) restriction)
# Check consistency:
print(f"  chi(rho_1*rho_5 - rho_4) = {[f'{x:.4f}' for x in chi_sum_68]}")
print(f"  chi(V_3)                 = {[f'{x:.4f}' for x in chi_V3]}")
print(f"  Match: {np.allclose(chi_sum_68, chi_V3, atol=1e-6)}")

# We need one more relation. Use rho_1 x rho_6 = rho_5 + rho_7
# chi_{rho_1 x rho_6} = chi_1 * chi_6
# chi_6 = chi_1 * chi_6... wait, we don't have chi_6 yet.

# Alternative: use V_4 (dim 9) = V_4|_{2I}
# V_4: chi_j=4 at each class
chi_V4 = np.array([su2_character(4, a) for a in half_angles])

# V_4 should decompose as rho_0 + rho_4 + rho_6 (? -- need to check via inner products)
# <chi_V4, chi_k> = (1/120) * sum class_size * chi_V4 * chi_k

# Actually let me just use McKay recursion to get rho_6 and rho_8.
# From the McKay graph: rho_1 x rho_6 = rho_5 + rho_7
# => chi_6 = (chi_1 * chi_6 - chi_5 - chi_7)... no, this is circular.

# Let me use: rho_6 + rho_8 = V_3
# AND: rho_1 x rho_7 = rho_6
# => chi_6 = chi_1 * chi_7

char_table[6] = char_table[1] * char_table[7]  # rho_6 = rho_1 x rho_7

# Then rho_8 = V_3 - rho_6
char_table[8] = chi_V3 - char_table[6]

print(f"\n  Character table of 2I:")
print(f"  {'':>8s}", end="")
for c in range(n_classes):
    print(f" C{c}({class_sizes[c]:>2d})", end="")
print()
for k in range(9):
    print(f"  rho_{k}({dims[k]:d})", end="")
    for c in range(n_classes):
        print(f" {char_table[k,c]:>6.3f}", end="")
    print()

# Verify orthogonality
print(f"\n  Orthogonality check (should be delta_ij):")
gram = np.zeros((9, 9))
for i in range(9):
    for j in range(9):
        gram[i, j] = np.sum(class_sizes * char_table[i] * char_table[j]) / n_group
print(f"  Max off-diagonal: {max(abs(gram[i,j]) for i in range(9) for j in range(9) if i != j):.6f}")
print(f"  Diagonal: {[f'{gram[i,i]:.4f}' for i in range(9)]}")

# Verify dimensions
print(f"  Dims from chi(identity): {[f'{char_table[k,0]:.1f}' for k in range(9)]}")
print(f"  Expected:                {list(dims)}")

# ============================================================
# STEP 2: BUILD GROUP ELEMENTS AND CLASSIFY
# ============================================================
print(f"\n{'='*70}")
print("STEP 2: GROUP ELEMENTS AND CONJUGACY CLASSES")
print(f"{'='*70}")

verts, _, _ = build_600cell()

def quat_to_su2(q):
    a, b, c, d = q
    return np.array([[a + 1j*b, c + 1j*d],
                     [-c + 1j*d, a - 1j*b]])

G = np.array([quat_to_su2(v) for v in verts])

# Classify elements by trace (= 2*cos(alpha))
traces = np.array([np.real(np.trace(G[i])) for i in range(n_group)])

# Map each element to its conjugacy class
class_of = np.zeros(n_group, dtype=int)
expected_traces = [2*np.cos(a) for a in half_angles]

for i in range(n_group):
    dists = [abs(traces[i] - et) for et in expected_traces]
    class_of[i] = np.argmin(dists)
    assert min(dists) < 1e-6, f"Element {i}: trace={traces[i]}, no class match"

# Verify class sizes
for c in range(n_classes):
    actual_size = np.sum(class_of == c)
    assert actual_size == class_sizes[c], f"Class {c}: expected {class_sizes[c]}, got {actual_size}"

print(f"  All 120 elements classified into 9 conjugacy classes. OK.")

# Find identity
id_idx = np.where(class_of == 0)[0][0]
print(f"  Identity at index: {id_idx}")

# ============================================================
# STEP 3: BUILD ALL 9 IRREPS VIA PROJECTION
# ============================================================
print(f"\n{'='*70}")
print("STEP 3: BUILD IRREPS VIA REGULAR REPRESENTATION")
print(f"{'='*70}")

# Regular representation: R(g) is the 120x120 permutation matrix
# R(g)_{ij} = delta(g*h_j, h_i) where h_j are group elements
# We need the multiplication table first.

print(f"  Building multiplication table...")
mult = np.zeros((n_group, n_group), dtype=int)
for i in range(n_group):
    for j in range(n_group):
        prod = G[i] @ G[j]
        dists = np.array([np.linalg.norm(prod - G[k]) for k in range(n_group)])
        mult[i, j] = np.argmin(dists)
        assert dists[mult[i, j]] < 1e-6

print(f"  OK. Now building irreps via character projection...")

# For each irrep rho_k of dimension d_k:
# P_k = (d_k / |G|) * sum_g chi_k(g)* R(g)
# P_k is the projection onto the d_k-dimensional irrep (appears d_k times in regular rep)
# To get one copy: take the first d_k columns of P_k's image.

# Instead of the 120x120 regular representation, we can work smarter:
# The irrep matrices rho_k(g) can be extracted from the PROJECTION.
# P_k = (d_k / |G|) sum_g chi_k(g)^* R(g) is a 120x120 matrix of rank d_k^2.
# The irrep rho_k(g) is the restriction of R(g) to the image of P_k,
# in a suitable d_k-dimensional subspace.

# Build P_k for each irrep
irreps = {}  # irreps[k] = dict with 'matrices': array of (120, dk, dk)

for k in range(9):
    dk = dims[k]

    # Build projection P_k = (dk/120) * sum_g chi_k(g)^* * R(g)
    P_k = np.zeros((n_group, n_group))
    for g_idx in range(n_group):
        c = class_of[g_idx]
        chi_val = char_table[k, c]
        # R(g) is the permutation: R(g)_{ij} = 1 if g * h_j = h_i
        # i.e., R(g)_{mult[g,j], j} = 1
        for j in range(n_group):
            P_k[mult[g_idx, j], j] += (dk / n_group) * chi_val

    # P_k should have rank dk^2
    rank_P = np.linalg.matrix_rank(P_k, tol=1e-6)

    # Extract dk orthonormal columns from the image of P_k
    U, S, Vt = np.linalg.svd(P_k)
    # Take first dk columns (one copy of the irrep)
    basis = U[:, :dk]  # 120 x dk

    # Build irrep matrices: rho_k(g) = basis^T @ R(g) @ basis
    rho_k_mats = np.zeros((n_group, dk, dk), dtype=complex)
    for g_idx in range(n_group):
        # R(g) @ basis: permute rows of basis
        Rg_basis = np.zeros((n_group, dk))
        for j in range(n_group):
            Rg_basis[mult[g_idx, j]] = basis[j]
        rho_k_mats[g_idx] = basis.T @ Rg_basis

    # Verify: trace should match character
    ok = True
    for c in range(n_classes):
        g_idx = np.where(class_of == c)[0][0]
        tr = np.real(np.trace(rho_k_mats[g_idx]))
        expected = char_table[k, c]
        if abs(tr - expected) > 0.1:
            ok = False

    irreps[k] = rho_k_mats
    print(f"  rho_{k} (dim {dk}): rank(P) = {rank_P}, traces OK: {ok}")

# ============================================================
# STEP 4: CG INTERTWINERS FOR MCKAY EDGES
# ============================================================
print(f"\n{'='*70}")
print("STEP 4: CG INTERTWINERS")
print(f"{'='*70}")

# For each McKay edge (k, l): rho_l appears in rho_1 x rho_k.
# The CG intertwiner Y_{kl}: C^{d_k} -> C^{d_1} x C^{d_l} ... no.
# More precisely: the Yukawa Y_{kl} is a d_k x d_l matrix such that
# rho_k(g) Y_{kl} = Y_{kl} rho_l(g) for all g  (if we think of it as
# the component of [D_F, a] between blocks k and l).
#
# Actually in the NCG setup: Y_{kl} is the intertwiner for the
# McKay decomposition rho_1 x rho_k = sum_l rho_l.
# It maps V_l -> V_1 tensor V_k, so it's d_1*d_k x d_l = 2*dk x dl.
#
# For D_F: the off-diagonal block is a dk x dl matrix (not 2dk x dl).
# This comes from contracting the rho_1 index with a fixed direction
# (the Higgs field direction in Connes' NCG).
#
# Simplest approach: project the tensor product rho_1 x rho_k onto rho_l.
# P_l = (d_l / |G|) sum_g chi_l(g)^* (rho_1(g) kron rho_k(g))
# This is a (2*dk) x (2*dk) matrix. Its image has dimension d_l.
# The intertwiner is the basis of this image.

cg_intertwiners = {}

for (k, l) in MCKAY_EDGES:
    dk, dl = dims[k], dims[l]
    d1 = dims[1]  # = 2

    # Projection onto rho_l in rho_1 x rho_k
    P_l = np.zeros((d1 * dk, d1 * dk), dtype=complex)
    for g_idx in range(n_group):
        c = class_of[g_idx]
        chi_l = char_table[l, c]
        tensor = np.kron(irreps[1][g_idx], irreps[k][g_idx])
        P_l += (dl / n_group) * np.conj(chi_l) * tensor

    # Extract the intertwiner: SVD of P_l
    U, S, Vt = np.linalg.svd(P_l)
    # Number of significant singular values should be dl
    n_sig = np.sum(S > 0.1)

    # The intertwiner is the first dl columns of U (or Vt^H)
    # Y: (2*dk) x dl matrix
    Y_raw = U[:, :dl].real  # take real part (characters are real for 2I)

    # Normalize to ||Y||_F = 1/sqrt(2)
    norm = np.linalg.norm(Y_raw, 'fro')
    Y = Y_raw * (1.0 / (np.sqrt(2) * norm)) if norm > 1e-10 else Y_raw

    cg_intertwiners[(k, l)] = Y
    print(f"  Edge ({k},{l}): dim {d1}*{dk} x {dl} = {d1*dk}x{dl}, "
          f"rank(P) = {n_sig}, ||Y||_F = {np.linalg.norm(Y, 'fro'):.6f}")

# ============================================================
# STEP 5: BUILD D_F WITH CG YUKAWA
# ============================================================
print(f"\n{'='*70}")
print("STEP 5: BUILD D_F WITH CG YUKAWA")
print(f"{'='*70}")

# D_F acts on H_F = direct_sum V_{rho_k}, dim = sum(dims) = 30.
# But the CG intertwiner Y_{kl} is a (2*dk) x dl matrix,
# not a dk x dl matrix!
#
# In Connes' NCG: D_F has blocks between V_k and V_l that are dk x dl.
# The CG intertwiner maps V_l -> V_1 tensor V_k, which is 2dk x dl.
# To get a dk x dl block, we need to CONTRACT the V_1 index.
#
# The Higgs field direction selects a vector v in V_1 = C^2.
# Then: D_F[k,l] = (v^T tensor I_dk) @ Y_{kl}, which is dk x dl.
#
# Natural choice for Higgs direction: v = (1, 0) or (1/sqrt(2), 1/sqrt(2)).
# In the framework: the Higgs direction should be determined by
# symmetry breaking. Let's try v = (1, 0) first.

v_higgs = np.array([1.0, 0.0])  # Higgs direction in V_1 = C^2

dim_H = sum(dims)  # = 30
offsets = [0]
for d in dims:
    offsets.append(offsets[-1] + d)

D_F = np.zeros((dim_H, dim_H))

for (k, l) in MCKAY_EDGES:
    dk, dl = dims[k], dims[l]
    Y_full = cg_intertwiners[(k, l)]  # (2*dk) x dl

    # Contract with Higgs direction: Y_contracted = (v^T kron I_dk) @ Y_full
    # v^T kron I_dk is a dk x (2*dk) matrix
    contraction = np.kron(v_higgs.reshape(1, -1), np.eye(dk))  # dk x 2dk
    Y_contracted = contraction @ Y_full  # dk x dl

    # Renormalize to ||Y||_F = 1/sqrt(2)
    norm = np.linalg.norm(Y_contracted, 'fro')
    if norm > 1e-10:
        Y_contracted *= 1.0 / (np.sqrt(2) * norm)

    D_F[offsets[k]:offsets[k]+dk, offsets[l]:offsets[l]+dl] = Y_contracted
    D_F[offsets[l]:offsets[l]+dl, offsets[k]:offsets[k]+dk] = Y_contracted.T

# Verify
print(f"  dim(H_F) = {dim_H}")
print(f"  ||D_F - D_F^T|| = {np.linalg.norm(D_F - D_F.T):.2e}")

D2 = D_F @ D_F
D4 = D2 @ D2
tr_D2 = np.trace(D2)
tr_D4 = np.trace(D4)
print(f"  Tr(D_F^2) = {tr_D2:.6f}")
print(f"  Tr(D_F^4) = {tr_D4:.6f}")

evals_D = np.linalg.eigvalsh(D_F)
nonzero = evals_D[np.abs(evals_D) > 1e-10]
print(f"  Nonzero eigenvalues: {len(nonzero)} / {dim_H}")
print(f"  Eigenvalues: {[f'{e:.4f}' for e in sorted(nonzero)]}")

# ============================================================
# STEP 6: SPECTRAL TORSION
# ============================================================
print(f"\n{'='*70}")
print("STEP 6: SPECTRAL TORSION WITH CG YUKAWA")
print(f"{'='*70}")

# Algebra elements
E = []
for k in range(9):
    e = np.zeros((dim_H, dim_H))
    e[offsets[k]:offsets[k]+dims[k], offsets[k]:offsets[k]+dims[k]] = np.eye(dims[k])
    E.append(e)

# 1-forms
one_forms = {}
for (k, l) in MCKAY_EDGES:
    dE_l = D_F @ E[l] - E[l] @ D_F
    one_forms[(k, l)] = E[k] @ dE_l
    dE_k = D_F @ E[k] - E[k] @ D_F
    one_forms[(l, k)] = E[l] @ dE_k

# Torsion
edge_keys = sorted(one_forms.keys())
torsion_vals = []
for e1 in edge_keys:
    for e2 in edge_keys:
        for e3 in edge_keys:
            T = np.trace(one_forms[e1] @ one_forms[e2] @ one_forms[e3] @ D_F)
            if abs(T) > 1e-12:
                torsion_vals.append((e1, e2, e3, T))

print(f"  Nonzero torsion triples: {len(torsion_vals)} / {len(edge_keys)**3}")

if torsion_vals:
    torsion_vals.sort(key=lambda x: -abs(x[3]))
    unique_T = sorted(set(round(abs(T), 8) for _, _, _, T in torsion_vals))

    print(f"  Distinct |T| values: {len(unique_T)}")
    print(f"  Max |T| = {max(abs(T) for _,_,_,T in torsion_vals):.8f}")
    print(f"  Min |T| = {min(abs(T) for _,_,_,T in torsion_vals):.8f}")

    print(f"\n  |T| distribution:")
    for t_val in unique_T[:10]:
        count = sum(1 for _,_,_,T in torsion_vals if abs(abs(T)-t_val) < 1e-7)
        print(f"    |T| = {t_val:.8f} (x{count})")

    total_torsion = sum(abs(T) for _,_,_,T in torsion_vals)
    print(f"\n  Total |torsion| = {total_torsion:.6f}")
    print(f"  Compare democratic: 12.000000")

    # Top entries
    print(f"\n  Top torsion values:")
    for e1, e2, e3, T in torsion_vals[:10]:
        print(f"    T({e1},{e2},{e3}) = {T:+.8f}")
else:
    print(f"  ALL TORSION VALUES ARE ZERO!")

# ============================================================
# STEP 7: COMPARE WITH DEMOCRATIC
# ============================================================
print(f"\n{'='*70}")
print("STEP 7: COMPARISON CG vs DEMOCRATIC")
print(f"{'='*70}")

print(f"  {'Quantity':>25s} {'Democratic':>12s} {'CG-weighted':>12s}")
print(f"  {'-'*25} {'-'*12} {'-'*12}")
print(f"  {'Tr(D^2)':>25s} {'8.000':>12s} {tr_D2:>12.4f}")
print(f"  {'Tr(D^4)':>25s} {'12.000':>12s} {tr_D4:>12.4f}")
print(f"  {'Nonzero evals':>25s} {'8':>12s} {len(nonzero):>12d}")
print(f"  {'Nonzero torsion':>25s} {'48':>12s} {len(torsion_vals):>12d}")
n_distinct = len(set(round(abs(T), 8) for _,_,_,T in torsion_vals)) if torsion_vals else 0
print(f"  {'Distinct |T|':>25s} {'1':>12s} {n_distinct:>12d}")

print(f"\n{'='*70}")
print("EXP-540 COMPLETE")
print(f"{'='*70}")

"""
EXP-131: CKM Matrix from H4 Representation Theory on 600-Cell Eigenspaces
===========================================================================

CONTEXT from exp129:
- The 600-cell adjacency matrix has 9 distinct eigenvalues with multiplicities
  {1^2, 2^2, 3^2, 4^2, 5^2, 6^2, 3^2, 4^2, 2^2} = H4 irrep dimensions
- The two 9-dim eigenspaces (eigenvalues 4*phi and 4-4*phi) are EXACTLY ORTHOGONAL
- exp129 conclusion: "The most promising path forward: understand the CKM as
  arising from the REPRESENTATION THEORY of H4 on the 9-dim eigenspaces"
- The CKM remains a PATTERN, not DERIVAT

QUESTION: Can the CKM mixing matrix be derived from H4 representation theory?

APPROACHES:
1. Build 600-cell, compute adjacency matrix eigenspaces
2. Focus on the 9-dim eigenspaces (mult=3^2): decompose as 3 x 3
3. Fermion wavefunctions via A^n applied to reference vertex
4. Up-Down overlap matrix as CKM candidate
5. Association scheme P-matrix spectral fingerprints
6. 3x3 block structure from H4 subgroup decomposition
7. Level-to-eigenspace mapping via spectral weights

Author: Claude (exp131)
Date: 2026-02-09
Category: DERIVAT / PATTERN / NEGATIVE (to be determined)
"""

import numpy as np
from itertools import product as iterproduct, permutations
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

print("=" * 80)
print("EXP-131: CKM FROM H4 REPRESENTATION THEORY ON 600-CELL EIGENSPACES")
print("=" * 80)

# ================================================================
# SETUP: PDG CKM values and quark quantum numbers
# ================================================================
print("\n" + "=" * 80)
print("SETUP: Experimental data and quantum numbers")
print("=" * 80)

# PDG 2024 CKM magnitudes
V_exp = np.array([
    [0.97373, 0.2243, 0.00382],
    [0.221,   0.975,  0.0408 ],
    [0.0086,  0.0415, 1.014  ]
])

V_exp_err = np.array([
    [0.00031, 0.0008, 0.00020],
    [0.004,   0.006,  0.0014 ],
    [0.0002,  0.0011, 0.029  ]
])

labels_up = ['u', 'c', 't']
labels_dn = ['d', 's', 'b']

# Quark quantum numbers (a,b) from exp099/100/112
a_up = np.array([3, 2, 4])    # u, c, t
b_up = np.array([-2, 1, 1])
a_dn = np.array([1, 1, -1])   # d, s, b
b_dn = np.array([0, 1, 4])

# Spectral levels n = 5a + 6b
n_up = 5 * a_up + 6 * b_up    # [3, 16, 26]
n_dn = 5 * a_dn + 6 * b_dn    # [5, 11, 19]

print(f"\nUp-type:   u(a={a_up[0]:+d},b={b_up[0]:+d},n={n_up[0]}), "
      f"c(a={a_up[1]:+d},b={b_up[1]:+d},n={n_up[1]}), "
      f"t(a={a_up[2]:+d},b={b_up[2]:+d},n={n_up[2]})")
print(f"Down-type: d(a={a_dn[0]:+d},b={b_dn[0]:+d},n={n_dn[0]}), "
      f"s(a={a_dn[1]:+d},b={b_dn[1]:+d},n={n_dn[1]}), "
      f"b(a={a_dn[2]:+d},b={b_dn[2]:+d},n={n_dn[2]})")

# Experimental phi-exponents
n_exp_12 = -np.log(V_exp[0, 1]) / np.log(PHI)
n_exp_23 = -np.log(V_exp[1, 2]) / np.log(PHI)
n_exp_13 = -np.log(V_exp[0, 2]) / np.log(PHI)

print(f"\nExperimental CKM as phi-exponents:")
print(f"  |V_us| = {V_exp[0,1]:.5f} = phi^(-{n_exp_12:.4f})")
print(f"  |V_cb| = {V_exp[1,2]:.5f} = phi^(-{n_exp_23:.4f})")
print(f"  |V_ub| = {V_exp[0,2]:.6f} = phi^(-{n_exp_13:.4f})")

# Utility functions
def build_ckm(s12, s23, s13, delta=0):
    """Build standard CKM parametrization."""
    c12 = np.sqrt(max(0, 1 - s12**2))
    c23 = np.sqrt(max(0, 1 - s23**2))
    c13 = np.sqrt(max(0, 1 - s13**2))
    V = np.array([
        [c12*c13, s12*c13, s13],
        [-s12*c23 - c12*s23*s13*np.exp(1j*delta),
         c12*c23 - s12*s23*s13*np.exp(1j*delta),
         s23*c13],
        [s12*s23 - c12*c23*s13*np.exp(1j*delta),
         -c12*s23 - s12*c23*s13*np.exp(1j*delta),
         c23*c13]
    ])
    return np.abs(V)

def compute_errors(V_pred, label="", show=True):
    """Compute errors against experimental CKM."""
    abs_err = np.abs(np.abs(V_pred) - V_exp)
    rel_err = abs_err / V_exp * 100
    sigma_err = abs_err / V_exp_err
    mean_rel = np.mean(rel_err)
    max_rel = np.max(rel_err)
    chi2 = np.sum(sigma_err**2)
    # Unitarity check
    prod = np.abs(V_pred) @ np.abs(V_pred).T
    unitarity_err = np.max(np.abs(prod - np.eye(3)))
    if show:
        print(f"\n--- {label} ---")
        print(f"{'Element':<10} {'Pred':>10} {'Exp':>10} {'Err%':>8} {'Sigma':>8}")
        print("-" * 50)
        for i in range(3):
            for j in range(3):
                name = f"V_{labels_up[i]}{labels_dn[j]}"
                print(f"  {name:<8} {np.abs(V_pred[i,j]):10.5f} {V_exp[i,j]:10.5f} "
                      f"{rel_err[i,j]:7.2f}% {sigma_err[i,j]:7.1f}s")
        print(f"\n  Mean relative error: {mean_rel:.2f}%")
        print(f"  Max  relative error: {max_rel:.2f}%")
        print(f"  Chi^2 (9 dof):       {chi2:.1f}")
        print(f"  Unitarity max error: {unitarity_err:.2e}")
    return mean_rel, max_rel, chi2


# ================================================================
# PART 1: BUILD 600-CELL AND EIGENDECOMPOSE
# ================================================================
print("\n" + "=" * 80)
print("PART 1: BUILD 600-CELL AND EIGENDECOMPOSE")
print("=" * 80)

def build_600cell():
    """Build 600-cell vertices on unit S^3."""
    phi = PHI
    iphi = 1.0 / phi
    verts = set()
    def add(v):
        n = np.sqrt(sum(x**2 for x in v))
        if n > 1e-10:
            verts.add(tuple(round(x/n, 10) for x in v))
    # 8 vertices: permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [+1, -1]:
            v = [0., 0., 0., 0.]; v[i] = float(s); add(v)
    # 16 vertices: (+-1/2)^4
    for ss in iterproduct([+1, -1], repeat=4):
        add([s * 0.5 for s in ss])
    # 96 vertices: even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base = [phi/2, 0.5, iphi/2, 0.0]
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            pv = [base[p[i]] for i in range(4)]
            nz = [i for i in range(4) if abs(pv[i]) > 1e-12]
            for ss in iterproduct([+1, -1], repeat=len(nz)):
                v = list(pv)
                for k, idx in enumerate(nz): v[idx] *= ss[k]
                add(v)
    return np.array(sorted([list(v) for v in verts]))

V600 = build_600cell()
N600 = len(V600)
print(f"  600-cell vertices: {N600}")
assert N600 == 120, f"Expected 120 vertices, got {N600}"

# Dot products and adjacency
dots = V600 @ V600.T
np.clip(dots, -1, 1, out=dots)

# Distance class 1: cos(theta) = phi/2 ~ 0.809
cos_nn = PHI / 2
A = (np.abs(dots - cos_nn) < 0.001).astype(float)
np.fill_diagonal(A, 0)
degrees = A.sum(axis=1)
print(f"  Adjacency: degree = {int(degrees[0])} (all same: {np.all(degrees == degrees[0])})")

# Full eigendecomposition
eigenvalues_full, eigenvectors_full = np.linalg.eigh(A)
idx_sort = np.argsort(-eigenvalues_full)
eigenvalues_full = eigenvalues_full[idx_sort]
eigenvectors_full = eigenvectors_full[:, idx_sort]

# Identify distinct eigenvalues
unique_evals = []
mults = []
tol = 0.01
i = 0
while i < len(eigenvalues_full):
    val = eigenvalues_full[i]
    count = 1
    while i + count < len(eigenvalues_full) and abs(eigenvalues_full[i + count] - val) < tol:
        count += 1
    unique_evals.append(val)
    mults.append(count)
    i += count

print(f"\n  9 distinct eigenvalues with multiplicities:")
eigenspace_names = ['12', '6phi', '4phi', '3', '0', '-2', '4-4phi', '-3', '6-6phi']
for k, (ev, m) in enumerate(zip(unique_evals, mults)):
    sqm = int(round(np.sqrt(m)))
    print(f"    lambda_{k} = {ev:10.4f} ({eigenspace_names[k]:>8s}), mult = {m:3d} = {sqm}^2")

# Build eigenspace index ranges
eigenspaces = []
idx_pos = 0
for k, (ev, m) in enumerate(zip(unique_evals, mults)):
    eigenspaces.append((ev, list(range(idx_pos, idx_pos + m)), m))
    idx_pos += m

# ================================================================
# PART 2: ASSOCIATION SCHEME P-MATRIX
# ================================================================
print("\n" + "=" * 80)
print("PART 2: ASSOCIATION SCHEME P-MATRIX")
print("=" * 80)

# Build distance class matrices
unique_cos_sorted = sorted(set(np.round(dots[0], 6)))
dist_matrices = []
shell_sizes = []
for cos_val in unique_cos_sorted:
    D_mat = (np.abs(dots - cos_val) < 0.001).astype(float)
    if abs(cos_val - 1.0) < 0.001:
        D_mat = np.eye(N600)
    dist_matrices.append(D_mat)
    shell_sizes.append(int(D_mat[0].sum()))

# Compute P-matrix: P[d][k] = eigenvalue of A_d in eigenspace k
P_matrix = np.zeros((len(dist_matrices), len(unique_evals)))
for d, Ad in enumerate(dist_matrices):
    for k, (ev, idxs, m) in enumerate(eigenspaces):
        v0 = eigenvectors_full[:, idxs[0]]
        P_matrix[d, k] = v0 @ Ad @ v0

# Sort distance classes by angle (reverse cos = increasing angle)
# Note: unique_cos_sorted goes from most negative (d=8, antipodal) to 1 (d=0, self)
# We want d=0 (self) first
dist_order = list(range(len(unique_cos_sorted)))[::-1]  # reverse order
P_sorted = P_matrix[dist_order]
shell_sorted = [shell_sizes[i] for i in dist_order]

print(f"\n  P-matrix (rows=distance class d, cols=eigenspace k):")
print(f"  {'d':>3s} {'size':>5s}", end="")
for k in range(len(unique_evals)):
    print(f"  {'k='+str(k):>8s}", end="")
print()
for d_idx, (row, size) in enumerate(zip(P_sorted, shell_sorted)):
    print(f"  {d_idx:3d} {size:5d}", end="")
    for val in row:
        print(f"  {val:8.3f}", end="")
    print()

# ================================================================
# PART 3: FERMION WAVEFUNCTIONS VIA A^n
# ================================================================
print("\n" + "=" * 80)
print("PART 3: FERMION WAVEFUNCTIONS VIA A^n")
print("=" * 80)

print("""
Strategy: For each fermion with level n, define a wavefunction psi_n on the
120 vertices by applying A^n to a reference vertex (vertex 0) and projecting
onto each eigenspace. Since A^n in eigenspace k gives eigenvalue theta_k^n,
the spectral weight of fermion n in eigenspace k is:

  w_k(n) = theta_k^n * mult_k / 120

The fermion's "spectral fingerprint" is the 9-vector (w_0, ..., w_8).
""")

# Reference vertex
ref = np.zeros(N600)
ref[0] = 1.0

# Compute A^n * ref for each level
# This is expensive for large n. Use eigendecomposition: A^n = V * diag(lam^n) * V^T
# psi_n = A^n * e_0 = sum_k lam_k^n * (v_k . e_0) * v_k (summed over all eigenvectors)

print("  Computing spectral fingerprints for each fermion level...")

all_levels = sorted(set(list(n_up) + list(n_dn) + [0, 17]))  # include e, tau
level_labels = {0: 'e', 3: 'u', 5: 'd', 11: 'mu/s', 16: 'c', 17: 'tau', 19: 'b', 26: 't'}

# For each eigenspace k, the projection of e_0:
# P_k * e_0 = sum_{i in k} (v_i . e_0) * v_i
# The squared norm of this projection = sum_{i in k} (v_i[0])^2 = mult_k / N (by vertex transitivity)

# The spectral weight at eigenspace k for level n:
# w_k(n) = theta_k^n * ||P_k * e_0||^2 = theta_k^n * mult_k / N

spectral_fingerprints = {}
for n in all_levels:
    fp = np.zeros(len(unique_evals))
    for k, (ev, idxs, m) in enumerate(eigenspaces):
        # ||P_k e_0||^2 = sum_i (v_i[0])^2 where i in eigenspace k
        proj_sq = sum(eigenvectors_full[0, idx]**2 for idx in idxs)
        fp[k] = ev**n * proj_sq
    # Normalize
    norm = np.linalg.norm(fp)
    if norm > 1e-30:
        fp_norm = fp / norm
    else:
        fp_norm = fp
    spectral_fingerprints[n] = (fp, fp_norm)

print(f"\n  Spectral fingerprints (raw, per eigenspace):")
print(f"  {'Level':>6s} {'Label':>6s}", end="")
for k in range(len(unique_evals)):
    print(f"  {'k='+str(k):>10s}", end="")
print()
for n in all_levels:
    fp_raw, fp_norm = spectral_fingerprints[n]
    lbl = level_labels.get(n, '?')
    print(f"  {n:6d} {lbl:>6s}", end="")
    for val in fp_raw:
        if abs(val) < 1e-6:
            print(f"  {'~0':>10s}", end="")
        else:
            print(f"  {val:10.3e}", end="")
    print()

# ================================================================
# PART 4: CKM FROM SPECTRAL FINGERPRINT OVERLAPS
# ================================================================
print("\n" + "=" * 80)
print("PART 4: CKM FROM SPECTRAL FINGERPRINT OVERLAPS")
print("=" * 80)

print("""
Approach: CKM_ij = <psi_up_i | psi_down_j> using normalized spectral fingerprints.
This tests whether the 9-eigenspace structure naturally gives CKM mixing.
""")

# 4a: Raw overlap of normalized fingerprints
print("--- 4a: Overlap of normalized 9-dim spectral fingerprints ---")
V_spec = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        _, fp_up = spectral_fingerprints[n_up[i]]
        _, fp_dn = spectral_fingerprints[n_dn[j]]
        V_spec[i, j] = np.dot(fp_up, fp_dn)

print(f"\n  Raw overlap matrix <psi_up | psi_dn>:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_spec[i,j]:8.5f}", end="")
    print()

print(f"\n  Analysis: This is NOT CKM-like (no clear diagonal dominance).")
print(f"  The issue: A^n grows exponentially with largest eigenvalue (12^n),")
print(f"  so ALL fingerprints are dominated by eigenspace k=0 (lambda=12).")

# 4b: Remove the trivial eigenspace (k=0, lambda=12)
print("\n--- 4b: Overlap excluding trivial eigenspace (k=0) ---")
V_spec2 = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        fp_up = spectral_fingerprints[n_up[i]][0].copy()
        fp_dn = spectral_fingerprints[n_dn[j]][0].copy()
        fp_up[0] = 0  # remove k=0
        fp_dn[0] = 0
        n1 = np.linalg.norm(fp_up)
        n2 = np.linalg.norm(fp_dn)
        if n1 > 1e-30 and n2 > 1e-30:
            V_spec2[i, j] = np.dot(fp_up, fp_dn) / (n1 * n2)

print(f"\n  Overlap without k=0:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_spec2[i,j]:8.5f}", end="")
    print()

# 4c: Use LOG of eigenvalue (logarithmic fingerprint)
print("\n--- 4c: Logarithmic spectral fingerprint ---")
print("  Instead of theta_k^n, use n*log|theta_k| to avoid exponential domination.")

spectral_log_fps = {}
for n in all_levels:
    fp = np.zeros(len(unique_evals))
    for k, (ev, idxs, m) in enumerate(eigenspaces):
        proj_sq = sum(eigenvectors_full[0, idx]**2 for idx in idxs)
        if abs(ev) > 1e-10:
            fp[k] = n * np.log(abs(ev)) * proj_sq * np.sign(ev)**n
        else:
            fp[k] = 0 if n > 0 else proj_sq
    norm = np.linalg.norm(fp)
    spectral_log_fps[n] = fp / norm if norm > 1e-30 else fp

V_log = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        V_log[i, j] = np.dot(spectral_log_fps[n_up[i]], spectral_log_fps[n_dn[j]])

print(f"\n  Log-spectral overlap:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_log[i,j]:8.5f}", end="")
    print()


# ================================================================
# PART 5: FOCUS ON 9-DIM EIGENSPACES (mult=3^2)
# ================================================================
print("\n" + "=" * 80)
print("PART 5: 9-DIM EIGENSPACES (mult=3^2) AND 3x3 STRUCTURE")
print("=" * 80)

print("""
The two eigenspaces with mult = 9 = 3^2 are our prime candidates for CKM.
  lambda_2 = 4*phi ~ 6.472  (mult=9)
  lambda_6 = 4-4*phi ~ -2.472  (mult=9)

If each decomposes as 3 x 3 under some subgroup, the CKM could be the
rotation between the two 3-dim factors. We explore this systematically.
""")

spaces_9 = []
for k, (ev, idxs, m) in enumerate(eigenspaces):
    if m == 9:
        vecs = eigenvectors_full[:, idxs]
        spaces_9.append((ev, vecs, k))
        print(f"  Eigenspace k={k}: lambda={ev:.4f}, 9 eigenvectors extracted")

# 5a: SVD between the two 9-dim spaces
print("\n--- 5a: Canonical angles between the two 9-dim eigenspaces ---")
ev1, vecs1, k1 = spaces_9[0]
ev2, vecs2, k2 = spaces_9[1]

overlap_9x9 = vecs1.T @ vecs2
U9, S9, V9t = np.linalg.svd(overlap_9x9)
print(f"  Singular values: {np.round(S9, 8)}")
print(f"  (These should all be 0 since eigenspaces of symmetric matrix are orthogonal)")
print(f"  Max |SV|: {np.max(np.abs(S9)):.2e}")
print(f"  --> CONFIRMED: The two 9-dim eigenspaces are EXACTLY ORTHOGONAL.")
print(f"  This means we CANNOT get CKM from overlap between them directly.")

# 5b: Decompose EACH 9-dim eigenspace into 3+3+3
print("\n--- 5b: Internal structure of each 9-dim eigenspace ---")

for space_idx, (ev_k, vecs_k, k_idx) in enumerate(spaces_9):
    print(f"\n  Eigenspace k={k_idx} (lambda={ev_k:.4f}):")

    # The 9 eigenvectors span a 9-dim subspace of R^120.
    # We need to find a natural 3x3 decomposition.

    # Method: Use the symmetry of the 600-cell.
    # The H4 Weyl group acts on this eigenspace.
    # Mult=9 suggests: 9-dim irrep, or 3-dim irrep with multiplicity 3.

    # Test: compute the Gram matrix of the 9 eigenvectors restricted to
    # specific vertex subsets. If the 600-cell has a 3-fold structure,
    # this should show up.

    # Method 1: Use adjacency in the eigenspace
    # A restricted to the eigenspace: A_k = vecs_k^T * A * vecs_k
    # This should be ev_k * I (since these are eigenvectors of A)
    A_restricted = vecs_k.T @ A @ vecs_k
    print(f"    A restricted to eigenspace: diag check = {np.allclose(A_restricted, ev_k * np.eye(9), atol=0.01)}")

    # Method 2: Use distance-2 matrix (A_2)
    # The distance-2 matrix has DIFFERENT eigenvalues in this space
    if len(dist_matrices) > 2:
        # dist_matrices are ordered by cos value (most negative = most distant first)
        # dist_order reverses: d=0 is self, d=1 is nearest, etc.
        A2 = dist_matrices[dist_order[2]]  # distance class 2
        A2_restricted = vecs_k.T @ A2 @ vecs_k
        eigs_A2 = np.linalg.eigvalsh(A2_restricted)
        eigs_A2_sorted = np.sort(eigs_A2)[::-1]
        print(f"    A_2 restricted eigenvalues: {np.round(eigs_A2_sorted, 4)}")

        # Check: does A_2 restricted break the 9-dim into 3+3+3?
        rounded = np.round(eigs_A2_sorted, 2)
        vals, counts = np.unique(rounded, return_counts=True)
        print(f"    A_2 eigenvalue multiplicities: {dict(zip(np.round(vals, 4), counts))}")

    # Method 3: Use ALL distance matrices to find which one splits 9 -> 3+3+3
    print(f"\n    Checking all distance matrices for 3+3+3 splitting:")
    for d_idx in range(len(dist_order)):
        Ad = dist_matrices[dist_order[d_idx]]
        Ad_res = vecs_k.T @ Ad @ vecs_k
        eigs_d = np.sort(np.linalg.eigvalsh(Ad_res))[::-1]
        rounded_d = np.round(eigs_d, 3)
        vals_d, counts_d = np.unique(rounded_d, return_counts=True)
        mults_str = "+".join(str(c) for c in sorted(counts_d, reverse=True))
        if 3 in counts_d:
            marker = " <-- HAS 3-FOLD!"
        else:
            marker = ""
        print(f"      d={d_idx}: mults = {mults_str}, vals = {dict(zip(np.round(vals_d,3), counts_d))}{marker}")


# ================================================================
# PART 6: H4 SUBGROUP DECOMPOSITION
# ================================================================
print("\n" + "=" * 80)
print("PART 6: H4 SUBGROUP DECOMPOSITION OF 9-DIM EIGENSPACE")
print("=" * 80)

print("""
The 600-cell has symmetry group H4 (the Coxeter group of type H4) with
|H4| = 14400. The 9-dim eigenspace transforms as some representation of H4.

Key question: does 9 = 3+3+3 under some physical subgroup?

H4 has subgroups including:
  - H3 (icosahedral symmetry, order 120)
  - D4, A4, etc.
  - Z5 (pentagonal rotations)

If we restrict H4 to Z3 (a cyclic subgroup), the 9-dim rep decomposes as
3 copies of 3 (one for each Z3 irrep). This is our 3x3 structure.

Approach: find a Z3 (or S3) symmetry of the 600-cell and use it to
decompose the 9-dim eigenspace.
""")

# 6a: Find a 3-fold rotational symmetry of the 600-cell
print("--- 6a: Finding 3-fold symmetry axes ---")

# A Z3 element permutes vertices in orbits of 3.
# For the 600-cell, Z3 rotations correspond to 120-degree rotations
# in 4D that preserve the polytope.

# Method: find pairs of vertices at 120-degree separation in some plane.
# The 600-cell has 200 triangular faces, giving many Z3 axes.

# Simpler approach: construct a specific Z3 from the geometry.
# The 600-cell has 1200 triangular faces. Each triangle defines a Z3.

# Find a triangle: 3 vertices mutually at edge distance
triangles = []
for i in range(N600):
    neighbors_i = set(np.where(A[i] > 0.5)[0])
    for j in neighbors_i:
        if j <= i:
            continue
        neighbors_j = set(np.where(A[j] > 0.5)[0])
        common = neighbors_i.intersection(neighbors_j)
        for k in common:
            if k > j:
                triangles.append((i, j, k))

print(f"  Found {len(triangles)} triangular faces (expect 1200)")

# Use the first triangle to define a Z3 permutation
if len(triangles) > 0:
    tri = triangles[0]
    v0, v1, v2 = V600[tri[0]], V600[tri[1]], V600[tri[2]]
    center = (v0 + v1 + v2) / 3
    center_norm = center / np.linalg.norm(center)
    print(f"  Triangle vertices: {tri}")
    print(f"  Triangle center direction: {np.round(center_norm, 4)}")

# 6b: Construct Z3 rotation matrix in R^4
print("\n--- 6b: Constructing Z3 rotation from triangle ---")

# The Z3 rotation that cycles v0->v1->v2->v0
# In R^4, a rotation that maps v0->v1, v1->v2, v2->v0
# We need the 4x4 rotation matrix R such that R*v0 = v1, R*v1 = v2, R*v2 = v0

# Build R from the constraint: R maps v0->v1, v1->v2
# We have 3 known image pairs, need the 4th direction.
# The 4th direction is the complement orthogonal to the triangle plane.

# Gram-Schmidt: build orthonormal basis in the triangle plane
e1 = v1 - v0
e1 = e1 / np.linalg.norm(e1)
e2 = v2 - v0
e2 = e2 - np.dot(e2, e1) * e1
e2 = e2 / np.linalg.norm(e2)
e3 = center_norm - np.dot(center_norm, e1) * e1 - np.dot(center_norm, e2) * e2
if np.linalg.norm(e3) > 1e-10:
    e3 = e3 / np.linalg.norm(e3)
else:
    # center_norm is in the span of e1, e2; find orthogonal direction differently
    e3 = np.zeros(4)
    # Use cross product analog in R^4: find vector orthogonal to e1, e2, center
    M3 = np.array([e1, e2, center_norm])
    # Null space of M3
    U3, S3, V3t = np.linalg.svd(M3)
    e3 = V3t[-1]

# Find e4 orthogonal to e1, e2, e3
M4 = np.array([e1, e2, e3])
U4, S4, V4t = np.linalg.svd(M4)
e4 = V4t[-1]

# In the basis {e1, e2, e3, e4}, the rotation should cycle the triangle
# Express v0, v1, v2 in this basis
basis = np.array([e1, e2, e3, e4])  # 4x4
coords0 = basis @ v0
coords1 = basis @ v1
coords2 = basis @ v2

print(f"  v0 in basis: {np.round(coords0, 6)}")
print(f"  v1 in basis: {np.round(coords1, 6)}")
print(f"  v2 in basis: {np.round(coords2, 6)}")

# The Z3 rotation in the e1-e2 plane is a 120-degree rotation
# R_{12}(120 deg): rotate by 2*pi/3 in the e1-e2 plane, leave e3, e4 fixed
# But we need to check this actually maps the triangle correctly.

# Build permutation matrix: for each vertex, find its image under R
# R * v_i = v_j means dot(R*v_i, v_j) ~ 1

# Actually, let's find the Z3 permutation by brute force:
# The rotation that maps v0->v1, v1->v2, v2->v0
# Solve for R: 4x4 matrix, with R*v_i = v_{i+1 mod 3}
# Plus R fixes the axis orthogonal to the triangle.

# Build: [v0 | v1 | v2 | n] * R_perm = [v1 | v2 | v0 | n]
# where n is the normal to the triangle (in the plane of center)
# Actually we need 4 points to define R.

# Use center as the 4th fixed point
c4 = center_norm
M_in = np.column_stack([v0, v1, v2, c4])
M_out = np.column_stack([v1, v2, v0, c4])

if abs(np.linalg.det(M_in)) > 1e-10:
    R_z3 = M_out @ np.linalg.inv(M_in)
    # Verify
    print(f"\n  R_z3 * v0 = v1? {np.allclose(R_z3 @ v0, v1, atol=1e-6)}")
    print(f"  R_z3 * v1 = v2? {np.allclose(R_z3 @ v1, v2, atol=1e-6)}")
    print(f"  R_z3 * v2 = v0? {np.allclose(R_z3 @ v2, v0, atol=1e-6)}")
    print(f"  R_z3^3 = I? {np.allclose(np.linalg.matrix_power(R_z3, 3), np.eye(4), atol=1e-6)}")

    # Check if R is a proper rotation (det=1, orthogonal)
    print(f"  det(R_z3) = {np.linalg.det(R_z3):.6f}")
    print(f"  R^T * R = I? {np.allclose(R_z3.T @ R_z3, np.eye(4), atol=1e-6)}")

    # Find the permutation of all 120 vertices under R_z3
    perm_z3 = np.zeros(N600, dtype=int)
    matched = 0
    for i in range(N600):
        Rv = R_z3 @ V600[i]
        dists = np.linalg.norm(V600 - Rv, axis=1)
        j = np.argmin(dists)
        if dists[j] < 1e-6:
            perm_z3[i] = j
            matched += 1
        else:
            perm_z3[i] = -1

    print(f"\n  Z3 permutation: {matched}/120 vertices matched")

    if matched == N600:
        # Orbit structure
        visited = set()
        orbits = {1: 0, 3: 0}
        for i in range(N600):
            if i not in visited:
                orbit = [i]
                visited.add(i)
                j = perm_z3[i]
                while j != i:
                    orbit.append(j)
                    visited.add(j)
                    j = perm_z3[j]
                orb_len = len(orbit)
                orbits[orb_len] = orbits.get(orb_len, 0) + 1

        print(f"  Orbit structure: {orbits}")
        print(f"  Fixed points: {orbits.get(1, 0)}, 3-cycles: {orbits.get(3, 0)}")
        n_fixed = orbits.get(1, 0)
        n_3cyc = orbits.get(3, 0)
        print(f"  Check: {n_fixed} + 3*{n_3cyc} = {n_fixed + 3*n_3cyc} (should be 120)")

        # 6c: Decompose 9-dim eigenspace under Z3
        print(f"\n--- 6c: Z3 decomposition of 9-dim eigenspaces ---")

        # The Z3 acts on eigenvectors by permuting vertex indices
        # P_z3 is the permutation matrix
        P_z3 = np.zeros((N600, N600))
        for i in range(N600):
            P_z3[i, perm_z3[i]] = 1.0

        # Verify P_z3 preserves adjacency
        print(f"  P_z3 * A * P_z3^T = A? {np.allclose(P_z3 @ A @ P_z3.T, A, atol=1e-6)}")

        omega = np.exp(2j * np.pi / 3)  # cube root of unity

        for space_idx, (ev_k, vecs_k, k_idx) in enumerate(spaces_9):
            print(f"\n  Eigenspace k={k_idx} (lambda={ev_k:.4f}):")

            # Z3 action on the 9-dim eigenspace
            # P_z3 acts on eigenvectors: for each eigvec v, P_z3 * v is in the same eigenspace
            # (since P_z3 commutes with A)
            P_z3_restricted = vecs_k.T @ P_z3 @ vecs_k  # 9x9 matrix

            # Eigenvalues of Z3 restricted: should be cube roots of unity
            z3_eigs = np.linalg.eigvals(P_z3_restricted)
            print(f"    Z3 eigenvalues: {np.sort(np.round(z3_eigs, 6))[::-1]}")

            # Count eigenvalues close to 1, omega, omega^2
            n_1 = sum(1 for e in z3_eigs if abs(e - 1) < 0.01)
            n_w = sum(1 for e in z3_eigs if abs(e - omega) < 0.01)
            n_w2 = sum(1 for e in z3_eigs if abs(e - omega**2) < 0.01)
            print(f"    Decomposition: {n_1} x trivial + {n_w} x omega + {n_w2} x omega^2")
            print(f"    Pattern: 9 = {n_1} + {n_w} + {n_w2}")

            # If we get 3+3+3, extract the three 3-dim subspaces
            if n_1 == 3 and n_w == 3 and n_w2 == 3:
                print(f"    ** PERFECT 3+3+3 DECOMPOSITION! **")

                # Diagonalize P_z3 restricted
                z3_vals, z3_vecs = np.linalg.eig(P_z3_restricted)

                # Group by eigenvalue
                idx_1 = [i for i, e in enumerate(z3_vals) if abs(e - 1) < 0.01]
                idx_w = [i for i, e in enumerate(z3_vals) if abs(e - omega) < 0.01]
                idx_w2 = [i for i, e in enumerate(z3_vals) if abs(e - omega**2) < 0.01]

                # These give three 3-dim subspaces of the 9-dim eigenspace
                # Lift back to R^120
                sub_1 = vecs_k @ z3_vecs[:, idx_1].real  # 120x3
                sub_w = vecs_k @ z3_vecs[:, idx_w]  # 120x3 (complex)
                sub_w2 = vecs_k @ z3_vecs[:, idx_w2]  # 120x3 (complex)

                print(f"    Sub_1 (trivial): shape {sub_1.shape}, orthonormal: {np.allclose(sub_1.T @ sub_1, np.eye(3), atol=0.01)}")

else:
    print(f"  WARNING: det(M_in) ~ 0, cannot build Z3 rotation this way.")
    print(f"  Using alternative approach...")


# ================================================================
# PART 7: PROJECTION OF FERMION LEVELS ONTO 9-DIM EIGENSPACES
# ================================================================
print("\n" + "=" * 80)
print("PART 7: FERMION LEVEL PROJECTION ONTO 9-DIM EIGENSPACES")
print("=" * 80)

print("""
Key idea: Each fermion at level n has a "wavefunction" psi_n on the 600-cell.
The most natural definition: psi_n is the indicator function of vertices at
graph distance n from a reference (for n <= diameter = 5), or more generally,
a function that encodes the spectral level.

Alternative: Use the P-matrix to define spectral wavefunctions.
The P-matrix row d gives the eigenvalue of the distance-d matrix in each eigenspace.
For level n = 5a + 6b, define the wavefunction amplitude in eigenspace k as:
  c_k(n) = P[d, k] where d is chosen to maximize discrimination.

But actually, the most physical approach: the fermion wavefunction IS the
projection onto the 9-dim eigenspace, weighted by the level n through the
Chebyshev polynomial relation.
""")

# 7a: Projection using P-matrix rows
print("--- 7a: P-matrix spectral profiles ---")

# The P-matrix maps distance classes to eigenspace amplitudes.
# For each eigenspace k, P[d,k] gives the eigenvalue of A_d in that eigenspace.
# The fermion at level n uses the n-th "spectral moment": theta_k^n.

# Within the 9-dim eigenspace (k=2, lambda=4*phi), the fermion at level n
# has amplitude proportional to (4*phi)^n. This is the same for all fermions
# in this eigenspace! So we need a DIFFERENT projection.

# The key insight: the fermion level n determines which LINEAR COMBINATION
# of eigenspaces to use. Within one eigenspace, all fermions look the same.
# The CKM comes from mixing between different eigenspaces.

# Let's construct the "spectral state" for each fermion across ALL 9 eigenspaces
print("  Spectral state: c_k(n) = theta_k^n (unnormalized)")
print(f"  Then the CKM would be: V_ij = <up_i | down_j> = sum_k c_k(n_up_i) * c_k(n_dn_j)")
print()

# But wait - this is just what we computed in Part 4. Let's be smarter.

# 7b: Relative weighting approach
print("--- 7b: Relative eigenspace weighting ---")
print("""
  Instead of theta_k^n (exponentially dominated), use the RELATIVE amplitude:
  r_k(n) = theta_k^n / sum_j theta_j^n
  This normalizes by the total "spectral weight" and prevents k=0 domination.
""")

for norm_mode in ['relative', 'signed_log', 'polynomial']:
    V_ckm = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            n_u = n_up[i]
            n_d = n_dn[j]

            if norm_mode == 'relative':
                # Relative: weight_k = theta_k^n / max_k(|theta_k|^n)
                w_u = np.array([ev**n_u for ev, _, _ in eigenspaces])
                w_d = np.array([ev**n_d for ev, _, _ in eigenspaces])
                w_u = w_u / np.max(np.abs(w_u)) if np.max(np.abs(w_u)) > 0 else w_u
                w_d = w_d / np.max(np.abs(w_d)) if np.max(np.abs(w_d)) > 0 else w_d

            elif norm_mode == 'signed_log':
                # Signed log: sgn(theta_k) * log(|theta_k|) * n
                w_u = np.array([np.sign(ev) * np.log(max(abs(ev), 1e-30)) * n_u for ev, _, _ in eigenspaces])
                w_d = np.array([np.sign(ev) * np.log(max(abs(ev), 1e-30)) * n_d for ev, _, _ in eigenspaces])

            elif norm_mode == 'polynomial':
                # Chebyshev-like: T_n(theta_k / 12) where 12 = max eigenvalue
                from numpy.polynomial.chebyshev import chebval
                w_u = np.array([chebval(ev/12.0, [0]*n_u + [1]) for ev, _, _ in eigenspaces])
                w_d = np.array([chebval(ev/12.0, [0]*n_d + [1]) for ev, _, _ in eigenspaces])

            # Weight by multiplicity
            mult_weights = np.array([m for _, _, m in eigenspaces], dtype=float)
            w_u *= np.sqrt(mult_weights)
            w_d *= np.sqrt(mult_weights)

            # Normalize
            nu = np.linalg.norm(w_u)
            nd = np.linalg.norm(w_d)
            if nu > 1e-30 and nd > 1e-30:
                V_ckm[i, j] = np.dot(w_u, w_d) / (nu * nd)
            else:
                V_ckm[i, j] = 0

    print(f"\n  Mode: {norm_mode}")
    print(f"         d         s         b")
    for ii in range(3):
        print(f"  {labels_up[ii]}  ", end="")
        for jj in range(3):
            print(f"  {V_ckm[ii,jj]:8.5f}", end="")
        print()

    # Check if it resembles CKM
    diag_avg = np.mean(np.abs(np.diag(V_ckm)))
    offdiag = [V_ckm[0,1], V_ckm[1,2], V_ckm[0,2]]
    print(f"  Diagonal average: {diag_avg:.4f}")
    print(f"  Key off-diag: V_us={V_ckm[0,1]:.5f}, V_cb={V_ckm[1,2]:.5f}, V_ub={V_ckm[0,2]:.5f}")


# ================================================================
# PART 8: VERTEX-BASED WAVEFUNCTIONS AND SPATIAL OVERLAP
# ================================================================
print("\n" + "=" * 80)
print("PART 8: VERTEX-BASED WAVEFUNCTIONS AND SPATIAL OVERLAP")
print("=" * 80)

print("""
Different approach: define fermion wavefunctions directly on vertices.

For fermion at level n, define:
  psi_n(v) = (A^n * delta_0)[v]  (propagation from vertex 0 in n steps)

This is the probability amplitude of reaching vertex v in n steps.
The CKM is the overlap:
  V_ij = <psi_up_i | psi_dn_j> / (||psi_up_i|| * ||psi_dn_j||)
""")

# Compute A^n * e_0 using eigendecomposition
# A^n = V * diag(lambda^n) * V^T
# (A^n * e_0)[v] = sum_k lambda_k^n * (V_{v,k} * V_{0,k})

def compute_wavefunction(n_level, ref_idx=0):
    """Compute psi_n = A^n * delta_{ref_idx} using eigendecomposition."""
    psi = np.zeros(N600)
    for k in range(len(eigenvalues_full)):
        lam_n = eigenvalues_full[k] ** n_level
        coeff = eigenvectors_full[ref_idx, k] * lam_n
        psi += coeff * eigenvectors_full[:, k]
    return psi

print("--- 8a: Raw A^n wavefunctions ---")
psi_up = []
psi_dn = []
for i in range(3):
    psi_u = compute_wavefunction(n_up[i])
    psi_d = compute_wavefunction(n_dn[i])
    psi_up.append(psi_u)
    psi_dn.append(psi_d)
    print(f"  psi_{labels_up[i]}(n={n_up[i]}): norm={np.linalg.norm(psi_u):.3e}, "
          f"max={np.max(np.abs(psi_u)):.3e}, #nonzero={np.sum(np.abs(psi_u) > 1e-10)}")
    print(f"  psi_{labels_dn[i]}(n={n_dn[i]}): norm={np.linalg.norm(psi_d):.3e}, "
          f"max={np.max(np.abs(psi_d)):.3e}, #nonzero={np.sum(np.abs(psi_d) > 1e-10)}")

# Overlap matrix
V_raw = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        V_raw[i, j] = np.dot(psi_up[i], psi_dn[j]) / (np.linalg.norm(psi_up[i]) * np.linalg.norm(psi_dn[j]))

print(f"\n  Raw overlap <psi_up | psi_dn> / norms:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_raw[i,j]:8.5f}", end="")
    print()

# 8b: Project out the dominant eigenspace (k=0)
print("\n--- 8b: Wavefunctions with k=0 projected out ---")

def compute_wavefunction_no_trivial(n_level, ref_idx=0):
    """Compute psi_n excluding the trivial eigenspace (largest eigenvalue)."""
    psi = np.zeros(N600)
    for k in range(len(eigenvalues_full)):
        if abs(eigenvalues_full[k] - unique_evals[0]) < 0.01:
            continue  # skip trivial eigenspace
        lam_n = eigenvalues_full[k] ** n_level
        coeff = eigenvectors_full[ref_idx, k] * lam_n
        psi += coeff * eigenvectors_full[:, k]
    return psi

psi_up_nt = [compute_wavefunction_no_trivial(n_up[i]) for i in range(3)]
psi_dn_nt = [compute_wavefunction_no_trivial(n_dn[i]) for i in range(3)]

V_nt = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        nu = np.linalg.norm(psi_up_nt[i])
        nd = np.linalg.norm(psi_dn_nt[j])
        V_nt[i, j] = np.dot(psi_up_nt[i], psi_dn_nt[j]) / (nu * nd)

print(f"\n  Non-trivial overlap:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_nt[i,j]:8.5f}", end="")
    print()

# 8c: Use ONLY the 9-dim eigenspace wavefunctions
print("\n--- 8c: Wavefunctions restricted to 9-dim eigenspaces ---")

for space_idx, (ev_k, vecs_k, k_idx) in enumerate(spaces_9):
    print(f"\n  Eigenspace k={k_idx} (lambda={ev_k:.4f}):")

    # Project psi_n onto this eigenspace
    # psi_n restricted = P_k * psi_n where P_k = vecs_k * vecs_k^T
    psi_up_9 = []
    psi_dn_9 = []

    for i in range(3):
        # Project onto eigenspace
        psi_u_full = compute_wavefunction(n_up[i])
        proj_u = vecs_k @ (vecs_k.T @ psi_u_full)  # 120-dim
        psi_up_9.append(proj_u)

        psi_d_full = compute_wavefunction(n_dn[i])
        proj_d = vecs_k @ (vecs_k.T @ psi_d_full)
        psi_dn_9.append(proj_d)

    # In the 9-dim eigenspace, psi_n is proportional to theta_k^n * P_k * e_0
    # So ALL wavefunctions in one eigenspace are PARALLEL (differ only by scalar)
    # This means the 3x3 overlap is rank-1! Not useful for CKM.

    V_9 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            nu = np.linalg.norm(psi_up_9[i])
            nd = np.linalg.norm(psi_dn_9[j])
            if nu > 1e-30 and nd > 1e-30:
                V_9[i, j] = np.dot(psi_up_9[i], psi_dn_9[j]) / (nu * nd)

    print(f"  Overlap matrix:")
    print(f"         d         s         b")
    for i in range(3):
        print(f"  {labels_up[i]}  ", end="")
        for j in range(3):
            print(f"  {V_9[i,j]:8.5f}", end="")
        print()

    # Check rank
    svd_vals = np.linalg.svd(V_9, compute_uv=False)
    print(f"  Singular values: {np.round(svd_vals, 6)}")
    print(f"  Rank: {np.sum(svd_vals > 1e-6)}")
    print(f"  --> Rank-1 as expected: all wavefunctions in one eigenspace are parallel.")


# ================================================================
# PART 9: MULTI-REFERENCE APPROACH
# ================================================================
print("\n" + "=" * 80)
print("PART 9: MULTI-REFERENCE APPROACH (BREAKING VERTEX-TRANSITIVITY)")
print("=" * 80)

print("""
The problem with single-reference wavefunctions: vertex-transitivity makes
all projections proportional. To get a non-trivial 3x3 matrix, we need
to BREAK the symmetry by using MULTIPLE reference vertices.

Idea: Use 3 reference vertices that correspond to the 3 generations.
The (a,b) quantum numbers define directions in the Z[phi] lattice.
Different reference vertices correspond to different "orientations"
of the fermion on the 600-cell.

Approach: For each generation g=1,2,3, pick a reference vertex at
graph distance related to the generation number, then compute wavefunctions.
""")

# 9a: Use vertices at specific graph distances as references
print("--- 9a: Graph-distance based references ---")

# Compute graph distances from vertex 0
# Using BFS
from collections import deque
def graph_distances(adj_matrix, source):
    N = adj_matrix.shape[0]
    dist = -np.ones(N, dtype=int)
    dist[source] = 0
    queue = deque([source])
    while queue:
        v = queue.popleft()
        for u in range(N):
            if adj_matrix[v, u] > 0.5 and dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist

gdist = graph_distances(A, 0)
print(f"  Graph distance distribution from vertex 0:")
for d in range(6):
    n_at_d = np.sum(gdist == d)
    print(f"    d={d}: {n_at_d} vertices")

# Pick 3 reference vertices at different distances
# Gen 1: vertex 0 (d=0)
# Gen 2: first vertex at d=1
# Gen 3: first vertex at d=2
ref_verts = [0]
for d in [1, 2]:
    candidates = np.where(gdist == d)[0]
    ref_verts.append(candidates[0])

print(f"\n  Reference vertices: {ref_verts}")
print(f"  At distances: {[gdist[v] for v in ref_verts]}")

# Compute wavefunctions from each reference
print("\n--- 9b: CKM from multi-reference wavefunctions ---")

# For up-type quark i: use reference ref_verts[i] with level n_up[i]
# For down-type quark j: use reference ref_verts[j] with level n_dn[j]
psi_up_mr = []
psi_dn_mr = []
for i in range(3):
    psi_u = compute_wavefunction(n_up[i], ref_idx=ref_verts[i])
    psi_d = compute_wavefunction(n_dn[i], ref_idx=ref_verts[i])
    psi_up_mr.append(psi_u / np.linalg.norm(psi_u))
    psi_dn_mr.append(psi_d / np.linalg.norm(psi_d))

V_mr = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        V_mr[i, j] = np.abs(np.dot(psi_up_mr[i], psi_dn_mr[j]))

print(f"\n  Multi-reference overlap |<psi_up_i | psi_dn_j>|:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_mr[i,j]:8.5f}", end="")
    print()

# 9c: Use (a,b) quantum numbers to define reference directions
print("\n--- 9c: (a,b)-directed wavefunctions ---")
print("""
  Each quark has quantum numbers (a,b). We use these to define a
  "direction" on the 600-cell via the Cholesky map from exp120d:
  (a,b) -> (a+b, b) in the E8 lattice.

  For each quark, we find the vertex closest to its (a,b) direction
  (projected to S^3) and use that as reference.
""")

# Map (a,b) to R^4 direction using the first 2 coordinates of each vertex
# This is ad hoc, but let's see what happens
# Actually, use the (a,b) as weights on two chosen basis vectors of R^4
# Basis v_a and v_b chosen from the 600-cell

# Use the first two basis vectors (1,0,0,0) and (0,1,0,0) as proxies
v_basis_a = np.array([1.0, 0.0, 0.0, 0.0])
v_basis_b = np.array([0.0, 1.0, 0.0, 0.0])

def find_nearest_vertex(direction):
    """Find vertex closest to given direction on S^3."""
    direction_norm = direction / np.linalg.norm(direction)
    dot_prods = V600 @ direction_norm
    return np.argmax(dot_prods)

# For each quark, define direction = a * v_a + b * v_b (on S^3)
quark_refs_up = []
quark_refs_dn = []
for i in range(3):
    dir_u = a_up[i] * v_basis_a + b_up[i] * v_basis_b
    dir_d = a_dn[i] * v_basis_a + b_dn[i] * v_basis_b
    if np.linalg.norm(dir_u) > 1e-10:
        ref_u = find_nearest_vertex(dir_u)
    else:
        ref_u = 0
    if np.linalg.norm(dir_d) > 1e-10:
        ref_d = find_nearest_vertex(dir_d)
    else:
        ref_d = 0
    quark_refs_up.append(ref_u)
    quark_refs_dn.append(ref_d)
    print(f"  {labels_up[i]}(a={a_up[i]:+d},b={b_up[i]:+d}): ref vertex {ref_u}, "
          f"v={np.round(V600[ref_u], 3)}")
    print(f"  {labels_dn[i]}(a={a_dn[i]:+d},b={b_dn[i]:+d}): ref vertex {ref_d}, "
          f"v={np.round(V600[ref_d], 3)}")

# Compute wavefunctions from (a,b)-directed references
psi_up_ab = []
psi_dn_ab = []
for i in range(3):
    psi_u = compute_wavefunction(n_up[i], ref_idx=quark_refs_up[i])
    psi_d = compute_wavefunction(n_dn[i], ref_idx=quark_refs_dn[i])
    psi_up_ab.append(psi_u / np.linalg.norm(psi_u))
    psi_dn_ab.append(psi_d / np.linalg.norm(psi_d))

V_ab = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        V_ab[i, j] = np.abs(np.dot(psi_up_ab[i], psi_dn_ab[j]))

print(f"\n  (a,b)-directed overlap |<psi_up_i | psi_dn_j>|:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_ab[i,j]:8.5f}", end="")
    print()

# Does it have CKM-like hierarchy?
diag_ab = [V_ab[i,i] for i in range(3)]
offdiag_12 = V_ab[0, 1]
print(f"\n  Diagonal: {[f'{x:.5f}' for x in diag_ab]}")
print(f"  V_us = {V_ab[0,1]:.5f}, V_cb = {V_ab[1,2]:.5f}, V_ub = {V_ab[0,2]:.5f}")


# ================================================================
# PART 10: P-MATRIX APPROACH TO CKM
# ================================================================
print("\n" + "=" * 80)
print("PART 10: P-MATRIX SPECTRAL DISTANCE APPROACH")
print("=" * 80)

print("""
The association scheme P-matrix encodes ALL spectral information.
P[d,k] = eigenvalue of the distance-d matrix in eigenspace k.

Idea: Define a "spectral distance" between two levels n1, n2 using the
P-matrix. The CKM element V_ij is related to exp(-spectral_distance).

For level n, the spectral profile is: f_n = (theta_0^n, ..., theta_8^n)
weighted by the P-matrix columns (which tell us how each eigenspace
contributes to each distance shell).

CKM_ij = exp(-D(up_i, dn_j)) / sum_j exp(-D(up_i, dn_j))
where D is the spectral distance.
""")

# 10a: Build spectral vectors using P-matrix
print("--- 10a: P-matrix spectral vectors ---")

# For level n, the amplitude at distance d is:
# A_d(n) = sum_k P[d,k] * theta_k^n * mult_k / N
# This is the expected number of paths of length n that end at distance d.

def spectral_vector_P(n_level):
    """Compute spectral vector for level n using P-matrix."""
    n_dist = len(dist_order)
    sv = np.zeros(n_dist)
    for d_idx in range(n_dist):
        for k in range(len(unique_evals)):
            ev = unique_evals[k]
            m = mults[k]
            sv[d_idx] += P_sorted[d_idx, k] * ev**n_level * m / N600
        # Normalize by shell size
        if shell_sorted[d_idx] > 0:
            sv[d_idx] /= shell_sorted[d_idx]
    return sv

print(f"  Spectral vectors (per distance shell):")
print(f"  {'Level':>6s} {'Label':>6s}", end="")
for d in range(len(dist_order)):
    print(f"  {'d='+str(d):>10s}", end="")
print()

sv_all = {}
for n in sorted(set(list(n_up) + list(n_dn))):
    sv = spectral_vector_P(n)
    sv_all[n] = sv
    lbl = level_labels.get(n, '?')
    print(f"  {n:6d} {lbl:>6s}", end="")
    for val in sv:
        print(f"  {val:10.3e}", end="")
    print()

# 10b: Spectral distance as CKM
print("\n--- 10b: CKM from spectral distance ---")

# Several distance measures
for dist_type in ['cosine', 'L2', 'KL']:
    V_pdist = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            sv_u = sv_all[n_up[i]]
            sv_d = sv_all[n_dn[j]]

            if dist_type == 'cosine':
                nu = np.linalg.norm(sv_u)
                nd = np.linalg.norm(sv_d)
                if nu > 1e-30 and nd > 1e-30:
                    V_pdist[i, j] = np.dot(sv_u, sv_d) / (nu * nd)
                else:
                    V_pdist[i, j] = 0

            elif dist_type == 'L2':
                d_L2 = np.linalg.norm(sv_u - sv_d) / max(np.linalg.norm(sv_u), np.linalg.norm(sv_d), 1e-30)
                V_pdist[i, j] = np.exp(-d_L2)

            elif dist_type == 'KL':
                # KL divergence on |sv| distributions
                p = np.abs(sv_u) + 1e-30
                q = np.abs(sv_d) + 1e-30
                p = p / np.sum(p)
                q = q / np.sum(q)
                kl = np.sum(p * np.log(p / q))
                V_pdist[i, j] = np.exp(-abs(kl))

    print(f"\n  Distance type: {dist_type}")
    print(f"         d         s         b")
    for ii in range(3):
        print(f"  {labels_up[ii]}  ", end="")
        for jj in range(3):
            print(f"  {V_pdist[ii,jj]:8.5f}", end="")
        print()


# ================================================================
# PART 11: DIRECT 3x3 EXTRACTION FROM 9-DIM EIGENSPACE
# ================================================================
print("\n" + "=" * 80)
print("PART 11: DIRECT 3x3 FROM 9-DIM EIGENSPACE VIA TENSOR DECOMPOSITION")
print("=" * 80)

print("""
The 9-dim eigenspace has mult = 3^2 = 9. This suggests a TENSOR PRODUCT
structure: 9 = 3 x 3, where the two factors correspond to "generation
space up" and "generation space down."

If we can find this tensor product decomposition, the CKM matrix is
the transformation between the two 3-dim factors.

Method: Reshape the 9 eigenvectors into a 3x3 matrix of R^120 vectors.
The "correct" reshaping is the one that makes the 3x3 structure
match the CKM hierarchy.
""")

# 11a: Try all splittings of 9-dim eigenspace into 3+3+3
print("--- 11a: Distance matrix splitting ---")

for space_idx, (ev_k, vecs_k, k_idx) in enumerate(spaces_9):
    print(f"\n  Eigenspace k={k_idx} (lambda={ev_k:.4f}):")

    best_split_D = None
    best_split_d = -1
    best_split_condition = 1e10

    for d_idx in range(1, len(dist_order)):  # skip d=0 (identity)
        Ad = dist_matrices[dist_order[d_idx]]
        Ad_res = vecs_k.T @ Ad @ vecs_k  # 9x9

        # Diagonalize
        evals_ad, evecs_ad = np.linalg.eigh(Ad_res)
        evals_sorted = np.sort(evals_ad)[::-1]
        rounded = np.round(evals_sorted, 3)
        vals, counts = np.unique(rounded, return_counts=True)

        # Check for 3+3+3 pattern
        if len(vals) == 3 and np.all(counts == 3):
            print(f"    d={d_idx}: PERFECT 3+3+3 split!")
            print(f"    Eigenvalues: {vals} (each with mult 3)")

            # Group eigenvectors by eigenvalue -> three 3-dim subspaces
            groups = {}
            for val in vals:
                mask = np.abs(evals_ad - val) < 0.05
                groups[round(val, 3)] = evecs_ad[:, mask]

            # For each pair of subspaces, compute the 3x3 overlap
            keys = sorted(groups.keys(), reverse=True)
            print(f"    Subspace eigenvalues (sorted): {keys}")

            for gi, ki in enumerate(keys):
                for gj, kj in enumerate(keys):
                    if gi >= gj:
                        continue
                    G_i = groups[ki]  # 9x3
                    G_j = groups[kj]  # 9x3
                    # Overlap in 9-dim space
                    O_ij = G_i.T @ G_j  # 3x3
                    print(f"\n    Overlap between subspace {ki} and {kj}:")
                    print(f"    {np.round(O_ij, 5)}")
                    sv_ij = np.linalg.svd(O_ij, compute_uv=False)
                    print(f"    SVD: {np.round(sv_ij, 5)}")

            best_split_D = groups
            best_split_d = d_idx
            break

        elif len(vals) <= 3 and 3 in counts:
            # Partial 3-fold structure
            cond = np.max(counts) - np.min(counts)
            if cond < best_split_condition:
                best_split_condition = cond
                best_split_d = d_idx

    if best_split_D is None:
        # No perfect 3+3+3 found. Try the best approximate
        print(f"    No perfect 3+3+3 split found.")
        print(f"    Best partial: d={best_split_d}")


# 11b: Tensor product decomposition via matrix reshaping
print("\n--- 11b: Tensor product decomposition ---")

for space_idx, (ev_k, vecs_k, k_idx) in enumerate(spaces_9):
    print(f"\n  Eigenspace k={k_idx} (lambda={ev_k:.4f}):")

    # The 9 eigenvectors form a 120x9 matrix.
    # Reshape as 120 x 3 x 3 tensor.
    # Then compute the "core tensor" to find the 3x3 structure.

    # Method: Higher-Order SVD (Tucker decomposition)
    # For a 3rd order tensor T of shape 120 x 3 x 3:
    # T ~ G x_1 U_1 x_2 U_2 x_3 U_3
    # where G is the core tensor.

    # Simpler: for each vertex v, we have a 3x3 matrix M_v = reshape(vecs_k[v, :], (3,3))
    # The SVD of the average |M_v| might reveal the CKM structure.

    # But the ordering of the 9 eigenvectors matters!
    # Let's try ALL possible orderings (9! is too many).
    # Instead, use the Gram matrix structure.

    # The Gram matrix of the 9 eigenvectors restricted to vertex subsets
    # might reveal the tensor structure.

    # Actually, let's use a more systematic approach:
    # The "entanglement" structure of the 120 x 9 matrix.
    # Reshape 120 x 9 -> (120*3) x 3 and compute its SVD.
    # This gives the best rank-3 approximation.

    for split_idx in range(3):
        # Try 3 different reshapings: [0:3|3:6|6:9], [0:3|6:9|3:6], etc
        perms_3 = [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]
        perm = perms_3[split_idx]
        cols = [list(range(3*p, 3*p+3)) for p in perm]

        # Block structure: M_ij = vecs_k[:, cols[i][j]]  -> three 120x3 blocks
        blocks = [vecs_k[:, cols[b]] for b in range(3)]

        # The 3x3 "transition matrix" between blocks i and j
        for bi in range(3):
            for bj in range(bi+1, bi+2):
                bj = bj % 3
                T_ij = blocks[bi].T @ blocks[bj]  # 3x3
                if split_idx == 0 and bi == 0:
                    print(f"    Perm {perm}, Block {bi}-{bj} overlap:")
                    print(f"    {np.round(T_ij, 5)}")
                    sv_t = np.linalg.svd(T_ij, compute_uv=False)
                    print(f"    SVD: {np.round(sv_t, 5)}")

    # The eigenvectors are orthonormal, so block overlaps are orthogonal
    # We need a different approach.

    # Method: Use the outer product structure.
    # If 9-dim = 3 x 3, then each eigenvector is a tensor product e_i x f_j.
    # Compute the 9x9 Gram matrix G = vecs_k^T @ vecs_k = I (orthonormal).
    # But in a SUBSPACE (e.g., restricted to 12 vertices), the Gram is NOT I.
    # The deviation from I reveals the tensor product structure.

    # Restrict to nearest neighbors of vertex 0 (12 vertices)
    nn_indices = np.where(A[0] > 0.5)[0]
    vecs_nn = vecs_k[nn_indices, :]  # 12 x 9

    # Gram of restricted vectors
    G_nn = vecs_nn.T @ vecs_nn  # 9x9
    print(f"\n    Gram matrix restricted to 12 nearest neighbors:")
    print(f"    Eigenvalues: {np.round(np.sort(np.linalg.eigvalsh(G_nn))[::-1], 4)}")

    # If tensor product, Gram has structure G = G_3 kron G_3
    # where G_3 is a 3x3 matrix.
    # Check: does G have eigenvalues with mult pattern 1+2+3+2+1?
    G_eigs = np.sort(np.linalg.eigvalsh(G_nn))[::-1]
    G_rounded = np.round(G_eigs, 3)
    G_vals, G_counts = np.unique(G_rounded, return_counts=True)
    print(f"    Gram eigenvalue multiplicities: {dict(zip(G_vals, G_counts))}")


# ================================================================
# PART 12: PHI-WEIGHTED WAVEFUNCTION APPROACH
# ================================================================
print("\n" + "=" * 80)
print("PART 12: PHI-WEIGHTED WAVEFUNCTIONS")
print("=" * 80)

print("""
New idea: Instead of A^n (which is dominated by largest eigenvalue),
use phi-WEIGHTED spectral decomposition:

  psi_n(v) = sum_k (theta_k / theta_max)^n * phi_k(v)

where phi_k(v) is the k-th eigenfunction at vertex v, normalized by theta_max.
This bounds all contributions to [-1, 1]^n.

Even better: use the P-matrix entries directly as weights:
  psi_n(v) = sum_k P[n mod 9, k] * phi_k(v)

where we use the level n modulo the diameter+1 (=6) as the distance class.
""")

# 12a: Normalized eigenvalue wavefunctions
print("--- 12a: Normalized eigenvalue wavefunctions ---")
theta_max = max(abs(ev) for ev, _, _ in eigenspaces)

def compute_wavefunction_normalized(n_level, ref_idx=0):
    """Compute wavefunction with eigenvalues normalized to [-1,1]."""
    psi = np.zeros(N600)
    for k_idx2 in range(len(eigenvalues_full)):
        lam_norm = (eigenvalues_full[k_idx2] / theta_max) ** n_level
        coeff = eigenvectors_full[ref_idx, k_idx2] * lam_norm
        psi += coeff * eigenvectors_full[:, k_idx2]
    return psi

psi_up_norm = [compute_wavefunction_normalized(n_up[i]) for i in range(3)]
psi_dn_norm = [compute_wavefunction_normalized(n_dn[i]) for i in range(3)]

V_norm_ckm = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        nu = np.linalg.norm(psi_up_norm[i])
        nd = np.linalg.norm(psi_dn_norm[j])
        V_norm_ckm[i, j] = np.dot(psi_up_norm[i], psi_dn_norm[j]) / (nu * nd)

print(f"\n  Normalized eigenvalue overlap:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_norm_ckm[i,j]:8.5f}", end="")
    print()

# 12b: Use A/12 (stochastic matrix) instead of A
print("\n--- 12b: Stochastic matrix (A/12)^n wavefunctions ---")

# (A/12)^n converges to uniform. The DEVIATIONS from uniform encode structure.
def compute_wavefunction_stochastic(n_level, ref_idx=0):
    """Compute wavefunction using stochastic matrix (A/12)^n."""
    psi = np.zeros(N600)
    for k_idx2 in range(len(eigenvalues_full)):
        lam_stoch = (eigenvalues_full[k_idx2] / 12.0) ** n_level
        coeff = eigenvectors_full[ref_idx, k_idx2] * lam_stoch
        psi += coeff * eigenvectors_full[:, k_idx2]
    # Subtract the uniform part (eigenspace k=0)
    psi -= 1.0 / N600  # uniform distribution
    return psi

psi_up_stoch = [compute_wavefunction_stochastic(n_up[i]) for i in range(3)]
psi_dn_stoch = [compute_wavefunction_stochastic(n_dn[i]) for i in range(3)]

V_stoch = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        nu = np.linalg.norm(psi_up_stoch[i])
        nd = np.linalg.norm(psi_dn_stoch[j])
        if nu > 1e-30 and nd > 1e-30:
            V_stoch[i, j] = np.dot(psi_up_stoch[i], psi_dn_stoch[j]) / (nu * nd)

print(f"\n  Stochastic deviation overlap:")
print(f"         d         s         b")
for i in range(3):
    print(f"  {labels_up[i]}  ", end="")
    for j in range(3):
        print(f"  {V_stoch[i,j]:8.5f}", end="")
    print()

# Does this show CKM hierarchy?
print(f"\n  Analysis:")
print(f"  Diagonal: {[f'{V_stoch[i,i]:.5f}' for i in range(3)]}")
print(f"  V_us={V_stoch[0,1]:.5f}, V_cb={V_stoch[1,2]:.5f}, V_ub={V_stoch[0,2]:.5f}")

# For CKM-like behavior, we need:
# diagonal ~ 1, V_us ~ 0.22, V_cb ~ 0.04, V_ub ~ 0.004
# Hierarchy: V_us >> V_cb >> V_ub


# ================================================================
# PART 13: GRAPH HEAT KERNEL APPROACH
# ================================================================
print("\n" + "=" * 80)
print("PART 13: GRAPH HEAT KERNEL APPROACH")
print("=" * 80)

print("""
The heat kernel K(t, i, j) = sum_k exp(-lambda_k * t) * phi_k(i) * phi_k(j)
where lambda_k are LAPLACIAN eigenvalues. At time t, K gives the probability
of diffusion from i to j.

For fermion at level n, use t = n * dt for some time scale dt.
The CKM is the overlap of heat-diffused wavefunctions.

This is more physical than A^n because the heat kernel decays (not grows)
with eigenvalue, giving natural suppression of high-energy modes.
""")

# Laplacian eigenvalues: lambda_L = degree - theta = 12 - theta
# So lambda_L ranges from 0 (theta=12) to 12-(-3.71) ~ 15.71

# Heat kernel from vertex 0 to vertex v at time t:
# K(t, 0, v) = sum_k exp(-(12 - theta_k) * t) * U_{0,k} * U_{v,k}

def heat_kernel_wavefunction(t_val, ref_idx=0):
    """Compute heat kernel diffusion from ref vertex at time t."""
    psi = np.zeros(N600)
    for k in range(len(eigenvalues_full)):
        lam_L = 12.0 - eigenvalues_full[k]  # Laplacian eigenvalue
        coeff = np.exp(-lam_L * t_val) * eigenvectors_full[ref_idx, k]
        psi += coeff * eigenvectors_full[:, k]
    return psi

# Try different time scales
for dt_name, dt_val in [("1/12", 1/12), ("1/5", 1/5), ("1", 1.0), ("phi", PHI), ("1/phi", 1/PHI)]:
    psi_up_heat = [heat_kernel_wavefunction(n_up[i] * dt_val) for i in range(3)]
    psi_dn_heat = [heat_kernel_wavefunction(n_dn[i] * dt_val) for i in range(3)]

    V_heat = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            nu = np.linalg.norm(psi_up_heat[i])
            nd = np.linalg.norm(psi_dn_heat[j])
            if nu > 1e-30 and nd > 1e-30:
                V_heat[i, j] = np.dot(psi_up_heat[i], psi_dn_heat[j]) / (nu * nd)

    # Check for CKM-like structure
    has_hierarchy = (V_heat[0,0] > V_heat[0,1] > V_heat[0,2] > 0 and
                     V_heat[1,1] > V_heat[1,0] and V_heat[1,1] > V_heat[1,2])

    marker = " <-- CKM-LIKE HIERARCHY!" if has_hierarchy else ""
    print(f"\n  dt={dt_name} (dt_val={dt_val:.4f}):{marker}")
    print(f"         d         s         b")
    for i in range(3):
        print(f"  {labels_up[i]}  ", end="")
        for j in range(3):
            print(f"  {V_heat[i,j]:8.5f}", end="")
        print()

    if has_hierarchy:
        # Compare with CKM more carefully
        # Normalize rows
        V_hn = V_heat.copy()
        for i in range(3):
            V_hn[i] /= np.linalg.norm(V_hn[i])
        print(f"  Row-normalized:")
        print(f"         d         s         b")
        for i in range(3):
            print(f"  {labels_up[i]}  ", end="")
            for j in range(3):
                print(f"  {V_hn[i,j]:8.5f}", end="")
            print()

# 13b: Heat kernel with different reference for each generation
print("\n--- 13b: Heat kernel with generation-specific references ---")

for dt_name, dt_val in [("1/12", 1/12), ("1/5", 1/5)]:
    psi_up_h2 = [heat_kernel_wavefunction(n_up[i] * dt_val, ref_idx=ref_verts[i]) for i in range(3)]
    psi_dn_h2 = [heat_kernel_wavefunction(n_dn[i] * dt_val, ref_idx=ref_verts[i]) for i in range(3)]

    V_h2 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            nu = np.linalg.norm(psi_up_h2[i])
            nd = np.linalg.norm(psi_dn_h2[j])
            if nu > 1e-30 and nd > 1e-30:
                V_h2[i, j] = np.dot(psi_up_h2[i], psi_dn_h2[j]) / (nu * nd)

    print(f"\n  dt={dt_name}, multi-reference:")
    print(f"         d         s         b")
    for i in range(3):
        print(f"  {labels_up[i]}  ", end="")
        for j in range(3):
            print(f"  {V_h2[i,j]:8.5f}", end="")
        print()


# ================================================================
# PART 14: SCAN FOR OPTIMAL CKM-GENERATING PARAMETERS
# ================================================================
print("\n" + "=" * 80)
print("PART 14: SYSTEMATIC SCAN FOR CKM-LIKE OVERLAP MATRICES")
print("=" * 80)

print("""
Brute-force scan over:
  - Time scale dt in heat kernel
  - Power law exponents
  - Different reference vertex configurations

Goal: Find if ANY combination of natural parameters produces a CKM-like matrix.
If yes, identify the parameters. If no, conclude that the CKM cannot be
derived from simple eigenspace overlaps.
""")

best_results = []

# 14a: Heat kernel dt scan
print("--- 14a: Heat kernel dt scan ---")
best_heat_err = 1e10
best_heat_dt = 0
best_heat_V = None

for dt_val in np.logspace(-2, 1, 200):
    psi_u = [heat_kernel_wavefunction(n_up[i] * dt_val) for i in range(3)]
    psi_d = [heat_kernel_wavefunction(n_dn[i] * dt_val) for i in range(3)]

    V_scan = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            nu = np.linalg.norm(psi_u[i])
            nd = np.linalg.norm(psi_d[j])
            if nu > 1e-30 and nd > 1e-30:
                V_scan[i, j] = np.abs(np.dot(psi_u[i], psi_d[j]) / (nu * nd))

    # Row-normalize
    for i in range(3):
        row_norm = np.linalg.norm(V_scan[i])
        if row_norm > 1e-30:
            V_scan[i] /= row_norm

    err = np.mean(np.abs(V_scan - V_exp) / V_exp * 100)
    if err < best_heat_err:
        best_heat_err = err
        best_heat_dt = dt_val
        best_heat_V = V_scan.copy()

print(f"  Best heat kernel dt = {best_heat_dt:.6f}")
print(f"  Best mean error = {best_heat_err:.2f}%")
if best_heat_V is not None:
    compute_errors(best_heat_V, f"14a: Best heat kernel (dt={best_heat_dt:.6f})")
    best_results.append(("Heat kernel", best_heat_err, best_heat_dt, "FITTING"))

    # Check if dt is a geometric quantity
    print(f"\n  dt = {best_heat_dt:.6f}")
    print(f"  dt * 12 = {best_heat_dt * 12:.6f}")
    print(f"  dt * phi = {best_heat_dt * PHI:.6f}")
    print(f"  1/dt = {1/best_heat_dt:.6f}")
    print(f"  dt ~ 1/phi^n: n = {-np.log(best_heat_dt)/np.log(PHI):.4f}")

# 14b: Power-law scan: V_ij = phi^(-alpha * |n_up - n_dn|) row-normalized
print("\n--- 14b: Power-law phi^(-alpha*|dn|) scan ---")
best_pow_err = 1e10
best_pow_alpha = 0

dn_matrix = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        dn_matrix[i, j] = abs(n_up[i] - n_dn[j])

for alpha in np.linspace(0.01, 3.0, 3000):
    V_pow = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            V_pow[i, j] = PHI**(-alpha * dn_matrix[i, j])
        V_pow[i] /= np.linalg.norm(V_pow[i])

    err = np.mean(np.abs(V_pow - V_exp) / V_exp * 100)
    if err < best_pow_err:
        best_pow_err = err
        best_pow_alpha = alpha
        best_pow_V = V_pow.copy()

print(f"  Best alpha = {best_pow_alpha:.6f}")
print(f"  Best mean error = {best_pow_err:.2f}%")
compute_errors(best_pow_V, f"14b: phi^(-{best_pow_alpha:.4f}*|dn|)")
best_results.append(("Power-law", best_pow_err, best_pow_alpha, "FITTING"))

# Check geometric significance of alpha
print(f"\n  alpha = {best_pow_alpha:.6f}")
for name, val in [("1/phi", 1/PHI), ("1/phi^2", 1/PHI**2), ("1/3", 1/3),
                   ("1/5", 1/5), ("1/6", 1/6), ("alpha_s", ALPHA_S),
                   ("sin^2(tW)/5", SIN2_TW/5), ("1/(2*phi)", 1/(2*PHI))]:
    err_pct = abs(best_pow_alpha - val) / best_pow_alpha * 100
    if err_pct < 20:
        print(f"    alpha ~ {name} = {val:.6f} ({err_pct:.1f}% off)")

# 14c: Two-parameter scan: V_ij = phi^(-alpha*|da| - beta*|db|)
print("\n--- 14c: Two-parameter phi^(-alpha*|da| - beta*|db|) scan ---")

da_matrix = np.zeros((3, 3))
db_matrix = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        da_matrix[i, j] = abs(a_up[i] - a_dn[j])
        db_matrix[i, j] = abs(b_up[i] - b_dn[j])

best_2p_err = 1e10
best_2p_alpha = 0
best_2p_beta = 0

for alpha in np.linspace(0.01, 3.0, 200):
    for beta in np.linspace(0.01, 3.0, 200):
        V_2p = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                V_2p[i, j] = PHI**(-(alpha * da_matrix[i, j] + beta * db_matrix[i, j]))
            V_2p[i] /= np.linalg.norm(V_2p[i])

        err = np.mean(np.abs(V_2p - V_exp) / V_exp * 100)
        if err < best_2p_err:
            best_2p_err = err
            best_2p_alpha = alpha
            best_2p_beta = beta
            best_2p_V = V_2p.copy()

print(f"  Best alpha = {best_2p_alpha:.4f}, beta = {best_2p_beta:.4f}")
print(f"  Best mean error = {best_2p_err:.2f}%")
compute_errors(best_2p_V, f"14c: phi^(-{best_2p_alpha:.3f}*|da| - {best_2p_beta:.3f}*|db|)")
best_results.append(("Two-param (da,db)", best_2p_err, (best_2p_alpha, best_2p_beta), "FITTING"))

# Geometric significance
print(f"\n  alpha/beta = {best_2p_alpha/best_2p_beta:.4f}")
print(f"  5*alpha = {5*best_2p_alpha:.4f}, 6*beta = {6*best_2p_beta:.4f}")
print(f"  alpha*6 + beta*5 = {best_2p_alpha*6 + best_2p_beta*5:.4f}")

# 14d: Compare with reference CKM formulas
print("\n--- 14d: Reference comparisons ---")

# Half-integer Coxeter
V_half = build_ckm(PHI**(-3), PHI**(-13/2), PHI**(-23/2))
m_half, _, _ = compute_errors(V_half, "Coxeter half-int: sin=phi^(-n), n=(3, 13/2, 23/2)")
best_results.append(("Coxeter half-int", m_half, (3, 6.5, 11.5), "PATTERN"))

# Integer Coxeter
V_int = build_ckm(PHI**(-3), PHI**(-7), PHI**(-12))
m_int, _, _ = compute_errors(V_int, "Coxeter integer: sin=phi^(-n), n=(3, 7, 12)")
best_results.append(("Coxeter integer", m_int, (3, 7, 12), "PATTERN"))


# ================================================================
# FINAL SUMMARY AND CONCLUSIONS
# ================================================================
print("\n" + "=" * 80)
print("FINAL SUMMARY: ALL APPROACHES")
print("=" * 80)

print(f"\n{'#':>3} {'Approach':<40s} {'Mean err%':>10s} {'Params':>20s} {'Status':<20s}")
print("-" * 95)
for idx, (name, err, params, status) in enumerate(sorted(best_results, key=lambda x: x[1])):
    print(f"  {idx+1:2d} {name:<40s} {err:10.2f} {str(params):>20s} {status}")

print(f"""

{'='*80}
CONCLUSIONS
{'='*80}

=== SYSTEMATIC RESULTS ===

1. SPECTRAL FINGERPRINT OVERLAP (Parts 4, 7):
   The spectral fingerprints c_k(n) = theta_k^n are dominated by the
   largest eigenvalue (12^n), making ALL overlap matrices nearly rank-1.
   Removing the trivial eigenspace or using logarithmic fingerprints
   does NOT produce CKM-like structure.
   STATUS: NEGATIVE RESULT

2. 9-DIM EIGENSPACE STRUCTURE (Parts 5, 6, 11):
   The two 9-dim eigenspaces (mult=3^2) are EXACTLY ORTHOGONAL (as expected
   for different eigenvalues of a symmetric matrix). Within each eigenspace,
   the distance matrices may or may not split as 3+3+3. When they do, the
   3x3 overlap between subspaces is orthogonal (trivially zero).
   The Z3 decomposition gives 3+3+3 but the subspaces encode symmetry,
   not generation structure.
   STATUS: NEGATIVE / INCONCLUSIVE

3. VERTEX-BASED WAVEFUNCTIONS (Parts 8, 9):
   A^n * delta_0 wavefunctions with single reference produce rank-1 overlaps
   in each eigenspace (all wavefunctions are parallel). Multi-reference
   approaches break this but introduce arbitrary choices.
   (a,b)-directed references give non-trivial overlaps but NOT CKM-like.
   STATUS: NEGATIVE RESULT

4. HEAT KERNEL (Part 13):
   Heat kernel wavefunctions at various time scales dt produce overlap
   matrices. Some dt values give hierarchical structure (V_us > V_cb > V_ub)
   but the NUMERICAL values do not match CKM without fitting dt.
   Best fit dt is a free parameter with no obvious geometric meaning.
   STATUS: FITTING (if dt is adjusted)

5. POWER-LAW FITS (Part 14):
   phi^(-alpha * |dn|) with 1 free parameter gives reasonable fits.
   phi^(-alpha*|da| - beta*|db|) with 2 free parameters gives better fits.
   These are FITTING, not derivation.

=== HONEST ASSESSMENT ===

The CKM matrix CANNOT be derived from simple eigenspace overlaps on the
600-cell. The fundamental obstacle is:

  ** In a vertex-transitive graph, single-reference wavefunctions
     in any eigenspace are ALL PARALLEL. There is no way to generate
     a non-trivial 3x3 overlap matrix within one eigenspace. **

The 9 = 3^2 structure of the eigenspace multiplicity is suggestive but
does NOT directly encode the CKM. The tensor product decomposition
3 x 3 relates to the H4 irrep structure, not to quark generations.

What DOES work (from exp129, confirmed here):
  - Cabibbo angle: sin(theta_12) = phi^(-3), 1.8% error. PATTERN.
  - Half-integer Coxeter: n = (3, 13/2, 23/2). PATTERN.
  - The hierarchy phi^(-|dn|) is natural but the coefficient requires fitting.

The CKM remains a PATTERN in the 600-cell framework.
""")

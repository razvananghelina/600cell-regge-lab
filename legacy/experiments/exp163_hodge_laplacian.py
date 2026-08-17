"""
EXP-163: Hodge Laplacian on 600-cell / 120-cell complex
=========================================================
Compute the Hodge Laplacian on k-forms (k=0,1,2) for the 600-cell simplicial complex.
The 1-form Laplacian encodes "vector field" dynamics - graviton-like modes.

Structure:
- 0-forms: functions on vertices (120)
- 1-forms: functions on edges (720)
- 2-forms: functions on faces (1200 triangles)
- 3-forms: functions on cells (600 tetrahedra)

Hodge Laplacian: Delta_k = d_{k-1} d_{k-1}^* + d_k^* d_k
where d_k: k-forms -> (k+1)-forms is the exterior derivative.

For gravitons: look at the spectrum of Delta_1 (1-form Laplacian).
Massless graviton = zero mode of Delta_1 on transverse-traceless sector.
"""

import numpy as np
from itertools import combinations, product, permutations
from collections import Counter, defaultdict
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
import sys

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-163: Hodge Laplacian on 600-Cell Complex")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell
# ============================================================
print("\n--- Step 1: Build 600-cell vertices ---")

verts_set = set()

for i in range(4):
    for s in [1, -1]:
        v = [0,0,0,0]; v[i] = s
        verts_set.add(tuple(v))

for signs in product([0.5, -0.5], repeat=4):
    verts_set.add(tuple(signs))

vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]
even_perms = []
for p in permutations(range(4)):
    inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
    if inv % 2 == 0:
        even_perms.append(p)

for perm in even_perms:
    base = [vals_base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for signs in product([1,-1], repeat=len(nz)):
        v = list(base)
        for idx, s in zip(nz, signs):
            v[idx] *= s
        verts_set.add(tuple(np.round(v, 10)))

verts = np.array(sorted(verts_set))
N_v = len(verts)
print(f"  Vertices: {N_v}")

# ============================================================
# Step 2: Build edges (720)
# ============================================================
print("\n--- Step 2: Build edges ---")

# Edge criterion: dot product = phi/2
dot_threshold = PHI / 2
dots = verts @ verts.T

edges = []
edge_set = set()
for i in range(N_v):
    for j in range(i+1, N_v):
        if abs(dots[i,j] - dot_threshold) < 0.001:
            edges.append((i, j))
            edge_set.add((i, j))

N_e = len(edges)
print(f"  Edges: {N_e}")

# Edge index lookup
edge_to_idx = {}
for idx, (i, j) in enumerate(edges):
    edge_to_idx[(i, j)] = idx
    edge_to_idx[(j, i)] = idx

# ============================================================
# Step 3: Build faces (triangles = 1200)
# ============================================================
print("\n--- Step 3: Build triangular faces ---")

# Build adjacency list
adj_list = defaultdict(set)
for i, j in edges:
    adj_list[i].add(j)
    adj_list[j].add(i)

# Triangle: three vertices mutually adjacent
triangles = []
for i in range(N_v):
    for j in adj_list[i]:
        if j > i:
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    triangles.append((i, j, k))

N_f = len(triangles)
print(f"  Faces (triangles): {N_f}")

# Face index lookup
face_to_idx = {}
for idx, (i, j, k) in enumerate(triangles):
    face_to_idx[(i, j, k)] = idx

# ============================================================
# Step 4: Build tetrahedra (cells = 600)
# ============================================================
print("\n--- Step 4: Build tetrahedral cells ---")

tetrahedra = []
for i in range(N_v):
    ni = adj_list[i]
    for j in ni:
        if j > i:
            common_ij = ni & adj_list[j]
            for k in common_ij:
                if k > j:
                    common_ijk = common_ij & adj_list[k]
                    for l in common_ijk:
                        if l > k:
                            tetrahedra.append((i, j, k, l))

N_c = len(tetrahedra)
print(f"  Cells (tetrahedra): {N_c}")
print(f"  Euler: V-E+F-C = {N_v}-{N_e}+{N_f}-{N_c} = {N_v - N_e + N_f - N_c}")

# ============================================================
# Step 5: Boundary operators
# ============================================================
print("\n--- Step 5: Build boundary operators ---")

# d0: 0-forms -> 1-forms (V x E)
# (d0 f)(e_{ij}) = f(j) - f(i) for oriented edge i->j (i<j)
print("  Building d0 (vertices -> edges)...")
d0 = lil_matrix((N_e, N_v), dtype=float)
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = 1.0
d0 = csr_matrix(d0)
print(f"    d0 shape: {d0.shape}")

# d1: 1-forms -> 2-forms (E x F)
# (d1 w)(f_{ijk}) = w(e_{ij}) + w(e_{jk}) - w(e_{ik})
# With orientation: triangle (i,j,k) with i<j<k
# boundary = (j,k) - (i,k) + (i,j)
print("  Building d1 (edges -> faces)...")
d1 = lil_matrix((N_f, N_e), dtype=float)
for f_idx, (i, j, k) in enumerate(triangles):
    # Boundary of (i,j,k): (j,k) - (i,k) + (i,j)
    e_ij = edge_to_idx.get((i, j))
    e_ik = edge_to_idx.get((i, k))
    e_jk = edge_to_idx.get((j, k))

    if e_ij is not None and e_ik is not None and e_jk is not None:
        # Orientation: (i,j) with i<j is +1
        d1[f_idx, e_ij] = 1.0   # (i,j)
        d1[f_idx, e_jk] = 1.0   # (j,k)
        d1[f_idx, e_ik] = -1.0  # -(i,k)

d1 = csr_matrix(d1)
print(f"    d1 shape: {d1.shape}")

# d2: 2-forms -> 3-forms (F x C)
print("  Building d2 (faces -> cells)...")
d2 = lil_matrix((N_c, N_f), dtype=float)

# For tetrahedron (i,j,k,l) with i<j<k<l
# boundary = (j,k,l) - (i,k,l) + (i,j,l) - (i,j,k)
for c_idx, (i, j, k, l) in enumerate(tetrahedra):
    faces_of_tet = [
        ((j, k, l), +1),
        ((i, k, l), -1),
        ((i, j, l), +1),
        ((i, j, k), -1),
    ]
    for face, sign in faces_of_tet:
        f_idx = face_to_idx.get(face)
        if f_idx is not None:
            d2[c_idx, f_idx] = sign

d2 = csr_matrix(d2)
print(f"    d2 shape: {d2.shape}")

# Verify: d1 * d0 = 0, d2 * d1 = 0
check1 = d1 @ d0
check2 = d2 @ d1
print(f"  |d1*d0| = {abs(check1).max():.2e} (should be 0)")
print(f"  |d2*d1| = {abs(check2).max():.2e} (should be 0)")

# ============================================================
# Step 6: Hodge Laplacians
# ============================================================
print("\n--- Step 6: Hodge Laplacians ---")

# Delta_0 = d0^T d0 (graph Laplacian on vertices)
# Delta_1 = d0 d0^T + d1^T d1 (Hodge-1 on edges)
# Delta_2 = d1 d1^T + d2^T d2 (Hodge-2 on faces)

print("  Computing Delta_0 = d0^T * d0...")
Delta_0 = d0.T @ d0
print(f"    Shape: {Delta_0.shape}")

# Eigenvalues of Delta_0
n_eigs_0 = min(30, N_v - 2)
evals_0 = eigsh(Delta_0.astype(float), k=n_eigs_0, which='SM', return_eigenvectors=False)
evals_0 = np.sort(evals_0)
print(f"    Smallest eigenvalues: {np.round(evals_0[:15], 4)}")
print(f"    Zero modes (< 0.01): {np.sum(np.abs(evals_0) < 0.01)}")

print("\n  Computing Delta_1 = d0*d0^T + d1^T*d1...")
Delta_1 = d0 @ d0.T + d1.T @ d1
print(f"    Shape: {Delta_1.shape}")

# Delta_1 is 720x720 - use sparse eigensolve
n_eigs_1 = min(50, N_e - 2)
evals_1 = eigsh(Delta_1.astype(float), k=n_eigs_1, which='SM', return_eigenvectors=False)
evals_1 = np.sort(evals_1)
print(f"    Smallest eigenvalues: {np.round(evals_1[:20], 4)}")
print(f"    Zero modes (< 0.01): {np.sum(np.abs(evals_1) < 0.01)}")
print(f"    First Betti number b1 = dim(ker Delta_1) = {np.sum(np.abs(evals_1) < 0.01)}")

print("\n  Computing Delta_2 = d1*d1^T + d2^T*d2...")
Delta_2 = d1 @ d1.T + d2.T @ d2
print(f"    Shape: {Delta_2.shape}")

n_eigs_2 = min(50, N_f - 2)
evals_2 = eigsh(Delta_2.astype(float), k=n_eigs_2, which='SM', return_eigenvectors=False)
evals_2 = np.sort(evals_2)
print(f"    Smallest eigenvalues: {np.round(evals_2[:20], 4)}")
print(f"    Zero modes (< 0.01): {np.sum(np.abs(evals_2) < 0.01)}")
print(f"    Second Betti number b2 = {np.sum(np.abs(evals_2) < 0.01)}")

# ============================================================
# Step 7: Betti numbers and topology
# ============================================================
print("\n--- Step 7: Betti numbers ---")
b0 = np.sum(np.abs(evals_0) < 0.01)
b1 = np.sum(np.abs(evals_1) < 0.01)
b2 = np.sum(np.abs(evals_2) < 0.01)
# b3 would need Delta_3, but we know chi = b0 - b1 + b2 - b3 = 0 (S^3)
b3_expected = b0 - b1 + b2  # from chi = 0

print(f"  b0 = {b0} (connected components)")
print(f"  b1 = {b1} (loops)")
print(f"  b2 = {b2} (voids)")
print(f"  b3 = {b3_expected} (from chi=0)")
print(f"  Expected for S^3: b0=1, b1=0, b2=0, b3=1")

# ============================================================
# Step 8: Spectral analysis of Delta_1 (graviton modes)
# ============================================================
print("\n--- Step 8: Delta_1 spectral analysis (graviton modes) ---")

# Get more eigenvalues with eigenvectors
n_eigs_detail = min(80, N_e - 2)
evals_1_full, evecs_1 = eigsh(Delta_1.astype(float), k=n_eigs_detail, which='SM')
sort_idx = np.argsort(evals_1_full)
evals_1_full = evals_1_full[sort_idx]
evecs_1 = evecs_1[:, sort_idx]

# Count multiplicities
evals_rounded = np.round(evals_1_full, 2)
unique_evals = []
current = evals_rounded[0]
count = 1
for i in range(1, len(evals_rounded)):
    if abs(evals_rounded[i] - current) < 0.05:
        count += 1
    else:
        unique_evals.append((current, count))
        current = evals_rounded[i]
        count = 1
unique_evals.append((current, count))

print(f"  Eigenvalue multiplicities (Delta_1):")
for val, mult in unique_evals[:15]:
    # Check if a + b*phi
    a_check = val
    b_check = 0
    for a in range(-20, 21):
        for b in range(-20, 21):
            if abs(a + b * PHI - val) < 0.05:
                a_check, b_check = a, b
                break
    phi_form = f"{a_check}+{b_check}*phi" if b_check != 0 else f"{a_check}"
    print(f"    lambda = {val:8.3f} (mult={mult:3d})  [{phi_form}]")

# ============================================================
# Step 9: Decompose Delta_1 spectrum into exact + co-exact
# ============================================================
print("\n--- Step 9: Exact / co-exact decomposition ---")

# Delta_1 = d0*d0^T + d1^T*d1
# Exact part: im(d0) -> eigenvalues of d0*d0^T
# Co-exact part: im(d1^T) -> eigenvalues of d1^T*d1
# Harmonic: ker(Delta_1) = ker(d0^T) cap ker(d1) [= b1 = 0 for S^3]

# Eigenvalues of d0^T*d0 (=Delta_0) give eigenvalues of d0*d0^T (same nonzero spectrum)
# Nonzero eigenvalues of d0*d0^T = nonzero eigenvalues of Delta_0
print(f"  Nonzero Delta_0 eigenvalues (= exact part of Delta_1):")
nz_0 = evals_0[np.abs(evals_0) > 0.01]
print(f"    Count: {len(nz_0)}")
print(f"    Values: {np.round(nz_0[:10], 3)}...")

# Eigenvalues of d1*d1^T give eigenvalues of d1^T*d1 (same nonzero spectrum)
# These contribute to the co-exact part
# Nonzero eigenvalues of d1^T*d1 = nonzero eigenvalues of d1*d1^T = nonzero eigenvalues of Delta_2 minus d2^T*d2 part
# Actually: nonzero evals of d1^T*d1 = nonzero evals of d1*d1^T

# Compute d1*d1^T spectrum
d1_d1T = d1 @ d1.T  # N_f x N_f = 1200 x 1200
n_eigs_coexact = min(50, N_f - 2)
evals_coexact = eigsh(d1_d1T.astype(float), k=n_eigs_coexact, which='SM', return_eigenvectors=False)
evals_coexact = np.sort(evals_coexact)
nz_coexact = evals_coexact[np.abs(evals_coexact) > 0.01]
print(f"\n  Co-exact part (from d1^T*d1):")
print(f"    Nonzero count in sample: {len(nz_coexact)}")
print(f"    Values: {np.round(nz_coexact[:10], 3)}...")

# The full Delta_1 spectrum = {exact eigenvalues} U {co-exact eigenvalues} U {0 x b1}
# Total: (N_v - b0) exact + (N_f - b2) co-exact + b1 harmonic = N_e
exact_count = N_v - b0
coexact_expected = N_f - b2  # assuming b2 from computation
harmonic_count = b1
print(f"\n  Decomposition count:")
print(f"    Exact: {exact_count} (from {N_v} vertices - {b0} zero modes)")
print(f"    Co-exact: N_f - b2 = {N_f} - {b2} = {N_f - b2}")
print(f"    Harmonic: {harmonic_count}")
print(f"    Total: {exact_count + N_f - b2 + harmonic_count} (should be {N_e})")

# ============================================================
# Step 10: Check for phi-structure in Delta_1 spectrum
# ============================================================
print("\n--- Step 10: Phi-structure in Delta_1 ---")

# Do the eigenvalues of Delta_1 live in Q(sqrt(5))?
# Check: are they of the form a + b*phi with a,b integer?

n_phi_form = 0
n_checked = 0
for val, mult in unique_evals:
    found = False
    for a in range(-50, 51):
        for b in range(-50, 51):
            if abs(a + b * PHI - val) < 0.03:
                found = True
                break
        if found:
            break
    if found:
        n_phi_form += 1
    n_checked += 1

print(f"  Of {n_checked} distinct eigenvalues, {n_phi_form} are of form a+b*phi")
print(f"  Spectral purity: {'PURE' if n_phi_form == n_checked else 'MIXED'}")

# ============================================================
# Step 11: Connection to 120-cell (dual complex)
# ============================================================
print("\n--- Step 11: Dual complex (120-cell) ---")

# The dual of the 600-cell triangulation of S^3 is related to the 120-cell
# Dual 0-cells = 600 tetrahedra (cells of 600-cell)
# Dual 1-cells = 1200 triangles
# Dual 2-cells = 720 edges
# Dual 3-cells = 120 vertices

# The Hodge star maps k-forms to (3-k)-forms on the dual
# So Delta_1 on primal <-> Delta_2 on dual (up to Hodge star)

# Check: spectrum of Delta_1 (primal, 720) vs Delta_2 (primal, 1200)
# They should share the same nonzero eigenvalues (Hodge duality on closed manifold)

print(f"  Primal Delta_1 nonzero eigenvalues (first 10): {np.round(evals_1_full[evals_1_full>0.01][:10], 3)}")
print(f"  Primal Delta_2 nonzero eigenvalues (first 10): {np.round(evals_2[evals_2>0.01][:10], 3)}")

# Check matching
nz1 = sorted(evals_1_full[evals_1_full > 0.01])
nz2 = sorted(evals_2[evals_2 > 0.01])
print(f"\n  Nonzero spectrum sizes: Delta_1 has {len(nz1)}, Delta_2 has {len(nz2)}")

# Hodge duality: nonzero spec(Delta_k) = nonzero spec(Delta_{n-k-1}) for n-manifold
# Here n=3, so spec(Delta_1) should match spec(Delta_2) (nonzero parts)
# This is because the Hodge star * commutes with the Laplacian on a Riemannian manifold

common = 0
for v1 in nz1[:30]:
    for v2 in nz2[:30]:
        if abs(v1 - v2) < 0.05:
            common += 1
            break

print(f"  Matching eigenvalues (in first 30): {common}/{min(30, len(nz1))}")

# ============================================================
# Step 12: Low-lying Delta_1 modes and their structure
# ============================================================
print("\n--- Step 12: Low-lying 1-form modes ---")

# The lowest nonzero eigenvalue of Delta_1 gives the "mass gap" for vector modes
# In gravity: massless graviton = zero mode (only if b1 > 0, which it isn't for S^3)
# The lowest mode gives the lightest massive vector excitation

lowest_nz = evals_1_full[evals_1_full > 0.01]
if len(lowest_nz) > 0:
    mass_gap = lowest_nz[0]
    print(f"  Mass gap (lowest nonzero Delta_1): {mass_gap:.4f}")

    # Check form
    for a in range(-20, 21):
        for b in range(-20, 21):
            if abs(a + b * PHI - mass_gap) < 0.05:
                print(f"    = {a} + {b}*phi (in Q(sqrt(5)))")

    # The mode structure
    mode_idx = np.argmin(np.abs(evals_1_full - mass_gap))
    mode = evecs_1[:, mode_idx]

    # How does this mode distribute over edge types?
    # Use edge classification (A-A, A-B, etc.)
    print(f"    Mode norm: {np.linalg.norm(mode):.4f}")
    print(f"    Mode max: {np.max(np.abs(mode)):.4f}")
    print(f"    Mode entries > 0.01: {np.sum(np.abs(mode) > 0.01)}/{N_e}")

# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
Simplicial complex: V={N_v}, E={N_e}, F={N_f}, C={N_c}
Euler characteristic: {N_v - N_e + N_f - N_c} (S^3 = 0)

Betti numbers: b0={b0}, b1={b1}, b2={b2}, b3={b3_expected}
Expected (S^3): 1, 0, 0, 1

Delta_0 (scalar): {n_eigs_0} eigenvalues computed
Delta_1 (vector/graviton): {n_eigs_detail} eigenvalues computed
Delta_2 (2-form): {n_eigs_2} eigenvalues computed

Key result: b1=0 means NO massless vector modes on S^3
(Consistent: S^3 is simply connected, no harmonic 1-forms)

For gravity: need symmetric 2-tensor modes, not just 1-forms.
Graviton = transverse-traceless part of metric perturbation =
  specific combination of 0-form and 2-form modes.

STATUS: DERIVED (topology), PATTERN/NEGATIVE (graviton modes)
""")

"""
exp470: Quantum Gravity from E8 Rotations -- EXPLORATORY
=========================================================
APPROACH: Work with standard E8 roots directly.
Find 600-cell subsets by analyzing the graph structure.

The icosian construction shows that E8 roots can be partitioned
into two sets of 120, each projecting to a 600-cell.
We want to count how many such partitions exist.

STATUS: EXPLORATORY -- not paper material
"""

import numpy as np
from itertools import product as iprod
from collections import Counter

phi = (1 + np.sqrt(5)) / 2
phi_c = (1 - np.sqrt(5)) / 2

print("=" * 70)
print("exp470: E8 ROOT SYSTEM AND 600-CELL DECOMPOSITIONS")
print("=" * 70)

# ============================================================
# PART A: Standard E8 root system (240 roots in R^8)
# ============================================================

E8_roots = []

# Type A: +-e_i +- e_j, i < j  (112 roots)
for i in range(8):
    for j in range(i+1, 8):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = np.zeros(8)
                v[i] = si; v[j] = sj
                E8_roots.append(v)

# Type B: (+-1/2)^8 with even number of minus signs (128 roots)
for signs in iprod([1, -1], repeat=8):
    if sum(1 for s in signs if s == -1) % 2 == 0:
        E8_roots.append(np.array(signs) / 2.0)

E8 = np.array(E8_roots)
print(f"E8 roots: {len(E8)} vectors in R^8")

# Inner product matrix
IP = E8 @ E8.T  # Euclidean inner product
ip_vals = sorted(set(np.round(IP[np.triu_indices(240, k=1)], 6)))
print(f"Inner products <r_i, r_j>: {ip_vals}")
print(f"Diagonal (norms^2): {sorted(set(np.round(np.diag(IP), 6)))}")

# ============================================================
# PART B: Graph structure of E8 root system
# ============================================================
print(f"\n{'='*70}")
print("PART B: E8 root graph (adjacency at inner product = 1)")
print(f"{'='*70}")

# In the root system, two roots are "adjacent" if <r_i, r_j> = 1
# (equivalently, the angle between them is 60 degrees)
# This gives the root graph.

adj_counts = Counter()
for i in range(240):
    n_adj = np.sum(np.abs(IP[i] - 1) < 0.01) - 0  # don't count self
    adj_counts[n_adj] += 1

print(f"Adjacency degree distribution (inner product = 1):")
for deg, count in sorted(adj_counts.items()):
    print(f"  degree {deg}: {count} roots")

# ============================================================
# PART C: 600-cell inner product structure
# ============================================================
print(f"\n{'='*70}")
print("PART C: What inner products should a 600-cell subset have?")
print(f"{'='*70}")

# Build the 600-cell in R^4 and compute its inner products
def construct_600cell():
    vertices = set()
    for i in range(4):
        for s in [1, -1]:
            v = [0.0]*4; v[i] = float(s)
            vertices.add(tuple(np.round(v, 10)))
    for signs in iprod([0.5, -0.5], repeat=4):
        vertices.add(tuple(np.round(signs, 10)))
    inv2phi = 1.0 / (2.0 * phi)
    base = [0.0, 0.5, phi/2.0, inv2phi]
    even_perms = [
        (0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
        (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0)
    ]
    for perm in even_perms:
        vals = [base[perm[i]] for i in range(4)]
        nz = [i for i in range(4) if abs(vals[i]) > 1e-10]
        for sc in iprod([1,-1], repeat=len(nz)):
            v = list(vals)
            for idx, pos in enumerate(nz):
                v[pos] *= sc[idx]
            vertices.add(tuple(np.round(v, 10)))
    return np.array(sorted(vertices))

V4 = construct_600cell()
print(f"600-cell: {len(V4)} vertices in R^4, all unit vectors")

# Inner products of 600-cell
IP_600 = V4 @ V4.T
ip_600_vals = sorted(set(np.round(IP_600[np.triu_indices(120, k=1)], 6)))
print(f"600-cell inner products (unit vectors in R^4):")
ip_600_count = Counter(np.round(IP_600[np.triu_indices(120, k=1)], 4))
for val, cnt in sorted(ip_600_count.items()):
    print(f"  <v_i, v_j> = {val:+.4f}: {cnt} pairs")

# ============================================================
# PART D: The icosian embedding -- using the CORRECT map
# ============================================================
print(f"\n{'='*70}")
print("PART D: Correct icosian -> E8 embedding")
print(f"{'='*70}")

print("""
The key insight: the 600-cell has 120 vertices with inner products
involving phi (golden ratio). These are NON-INTEGER, so the vertices
cannot be a subset of E8 roots (which have integer inner products)
in the SAME R^4 with the SAME Euclidean metric.

The icosian construction DOUBLES the dimension: R^4 -> R^8, mapping
each Z[phi]-coordinate to a Z^2 pair. The resulting 120 vectors in
R^8 have integer inner products under the TRACE bilinear form, but
NOT under the standard Euclidean form.

IMPORTANT: To get E8 roots (with standard Euclidean inner product),
we need to find the right change-of-basis matrix.
""")

# The correct embedding (following Wilson, "The Finite Simple Groups"):
# Map Z[phi] -> Z^2 using the "Minkowski embedding":
# z = a + b*phi  ->  (z, sigma(z)) = (a + b*phi, a + b*phi')
#
# Then the norm form is:
# z * sigma(z) = (a + b*phi)(a + b*phi') = a^2 + ab(phi+phi') + b^2*phi*phi'
# = a^2 + ab - b^2
# This is the NORM, not the TRACE.
#
# The TRACE is z + sigma(z) = 2a + b.
#
# The inner product <z_1, z_2> = z_1*sigma(z_2) + sigma(z_1)*z_2
# = 2*Re(z_1 * sigma(z_2))... hmm, this is getting complicated.
#
# Let me use Minkowski coordinates directly.
# Map quaternion q = (q_0, q_1, q_2, q_3) with q_m = a_m + b_m*phi to:
# R^8 vector: (q_0, sigma(q_0), q_1, sigma(q_1), q_2, sigma(q_2), q_3, sigma(q_3))
# = (a_0+b_0*phi, a_0+b_0*phi', a_1+b_1*phi, a_1+b_1*phi', ...)

def to_minkowski(v4):
    """Map R^4 vertex to R^8 via Minkowski embedding."""
    r8 = np.zeros(8)
    for i in range(4):
        x = v4[i]  # = a + b*phi
        # Find a, b
        for a2 in range(-4, 5):
            for b2 in range(-4, 5):
                a, b = a2/2.0, b2/2.0
                if abs(a + b*phi - x) < 1e-8:
                    r8[2*i] = a + b*phi      # z
                    r8[2*i+1] = a + b*phi_c   # sigma(z)
                    break
            else:
                continue
            break
    return r8

S_mink = np.array([to_minkowski(v) for v in V4])
print(f"Minkowski embedding: {S_mink.shape}")

# Check norms (standard Euclidean)
S_mink_norms = np.sqrt(np.sum(S_mink**2, axis=1))
print(f"Norms: {sorted(set(np.round(S_mink_norms, 6)))}")

# Check inner products
IP_mink = S_mink @ S_mink.T
ip_mink_vals = sorted(set(np.round(IP_mink[np.triu_indices(120, k=1)], 6)))
print(f"Inner products (Euclidean): {ip_mink_vals}")

# ============================================================
# PART E: Alternative -- find 600-cell GRAPH in E8 root graph
# ============================================================
print(f"\n{'='*70}")
print("PART E: Finding 600-cell subgraph in E8 root graph")
print(f"{'='*70}")

# Instead of finding an explicit embedding, let's look at the
# GRAPH structure.
#
# The 600-cell graph: 120 vertices, each with 12 neighbors.
# Two vertices are adjacent iff inner product = phi/2 in R^4.
#
# The E8 root graph: 240 vertices, each with degree...
# Let's compute it.

# E8 adjacency at inner product = 1 (angle 60 deg)
E8_adj_1 = (np.abs(IP - 1) < 0.01).astype(int)
np.fill_diagonal(E8_adj_1, 0)
print(f"E8 graph at <r,r'>=1: degree = {np.sum(E8_adj_1, axis=1)[0]}")

# E8 adjacency at inner product = -1 (angle 120 deg)
E8_adj_m1 = (np.abs(IP + 1) < 0.01).astype(int)
np.fill_diagonal(E8_adj_m1, 0)
print(f"E8 graph at <r,r'>=-1: degree = {np.sum(E8_adj_m1, axis=1)[0]}")

# Degree count for each inner product value
print(f"\nDegree at each inner product value:")
for val in sorted(set(np.round(IP[0], 2))):
    if abs(val - 2) > 0.01:  # skip self
        count = np.sum(np.abs(IP[0] - val) < 0.01)
        print(f"  <r_0, r_j> = {val:+.1f}: {count} roots")

# ============================================================
# PART F: The 120+120 decomposition via INVOLUTIONS
# ============================================================
print(f"\n{'='*70}")
print("PART F: Involutions of the E8 root system")
print(f"{'='*70}")

print("""
An involution sigma of E8 that maps r -> -sigma(r) for all roots
would partition roots into pairs {r, sigma(r)} with r != sigma(r).
But we want pairs {r, sigma(r)} where sigma preserves the root
system and has no fixed roots.

For the icosian decomposition, sigma is Galois conjugation:
it swaps phi <-> phi', preserving Z but changing the lattice embedding.

KEY QUESTION: Since we can't easily construct sigma in standard E8
coordinates, let's count involutions from the OPPOSITE direction.

An involution w in W(E8) with w^2 = 1 and no fixed roots partitions
240 roots into 120 pairs.

The conjugacy classes of involutions in W(E8) are known.
W(E8) has several classes of involutions, characterized by
their eigenvalue signatures on R^8 (number of +1 and -1 eigenvalues).
""")

# Let's construct some explicit involutions:

# INVOLUTION 1: negation of 4 coordinates
# w: (x1,...,x8) -> (-x1,...,-x4, x5,...,x8)
# This preserves E8 roots. Its fixed roots are those with x1=...=x4=0.
# But no Type B roots have x1=...=x4=0 (they're all +-1/2).
# Type A roots with both nonzero entries in {5,6,7,8}: that's C(4,2)*4 = 24.
# So 24 roots are fixed -> NOT a fixed-point-free involution.

# INVOLUTION 2: w(x) = -x (central involution)
# Fixed points: r such that -r = r, i.e. r = 0. Only trivial fixed point.
# But w maps each root to its negative, which IS also a root.
# This gives 120 pairs {r, -r}.
# Each pair gives a 1-dimensional root subsystem (A1).
# Choosing one from each pair gives a "half" of the root system.

# Let's see: if we pick all roots with (say) first nonzero coordinate > 0,
# we get a POSITIVE root system.

positive_roots = []
for r in E8:
    for c in r:
        if abs(c) > 1e-10:
            if c > 0:
                positive_roots.append(r)
            break

positive_roots = np.array(positive_roots)
print(f"\nPositive roots: {len(positive_roots)} (expect 120)")

# Check: are these a 600-cell? Check degree sequence.
IP_pos = positive_roots @ positive_roots.T
# Adjacency at <r, r'> = 1
adj_pos = (np.abs(IP_pos - 1) < 0.01).astype(int)
np.fill_diagonal(adj_pos, 0)
degrees_pos = np.sum(adj_pos, axis=1)
print(f"Positive root graph (adj at <r,r'>=1):")
print(f"  Degree distribution: {Counter(degrees_pos)}")
print(f"  (600-cell has degree 12 for all vertices)")

# So the positive root system does NOT form a regular graph.
# It's not the 600-cell graph.

# ============================================================
# PART G: The correct decomposition -- orbit of H4
# ============================================================
print(f"\n{'='*70}")
print("PART G: Understanding the icosian decomposition")
print(f"{'='*70}")

print("""
The icosian construction works as follows (Koca et al., 2014):

1. The 240 E8 roots can be decomposed under the H4 Coxeter group
   into ORBITS. Since |W(H4)| = 14400 and there are 240 roots,
   the orbit structure depends on how H4 embeds in W(E8).

2. The key mathematical fact: the E8 root system, when projected
   from R^8 to R^4 via a SPECIFIC projection, gives the union of
   TWO 600-cells (at different scales).

3. This projection maps: 240 roots -> 120 at norm sqrt(2) (on a
   600-cell) + 120 at norm sqrt(2*phi') (on another 600-cell).
   WRONG -- let me reconsider.

Actually, the correct statement (Moody & Patera, 1993):

   The E8 roots, under the projection R^8 -> R^4 that selects the
   "physical" subspace, map to:
   - 120 points forming a 600-cell (these ARE distinct roots)
   - 120 MORE points forming ANOTHER 600-cell at a different scale

   Both 600-cells are in the SAME R^4.
   The 240 roots project to 240 DISTINCT points in R^4.
""")

# Let's verify this by constructing the projection explicitly.
# We need the H4 embedding in R^8.

# The H4 root system has simple roots that can be expressed in R^4.
# In R^8, we embed H4 by identifying a 4D subspace.

# Actually, let's use a known result: the Coxeter element of E8
# has eigenvalues exp(2*pi*i*m/30) for m = 1, 7, 11, 13, 17, 19, 23, 29
# (the exponents of E8 are 1, 7, 11, 13, 17, 19, 23, 29).
# Under this Coxeter element, R^8 decomposes into 4 two-dimensional
# eigenspaces corresponding to pairs (m, 30-m).

# The pair (1, 29) and (11, 19) span a 4D subspace = "physical" H4
# The pair (7, 23) and (13, 17) span the complement = "internal" H4'

# Let's construct the Coxeter element.
# First, we need the simple roots of E8.

# E8 simple roots (standard Dynkin diagram ordering):
# alpha_1 = e_1 - e_2
# alpha_2 = e_1 + e_2
# alpha_3 = e_2 - e_3
# alpha_4 = e_3 - e_4
# alpha_5 = e_4 - e_5
# alpha_6 = e_5 - e_6
# alpha_7 = e_6 - e_7
# alpha_8 = -(1/2)(e_1+e_2+e_3+e_4+e_5+e_6+e_7-e_8)
# Wait, this might not be right. Let me use a standard reference.

# Standard E8 simple roots (Bourbaki):
# alpha_1 = (1/2)(e_1 - e_2 - e_3 - e_4 - e_5 - e_6 - e_7 + e_8)
# alpha_2 = e_1 + e_2
# alpha_i = e_{i-1} - e_{i-2} for i = 3,...,8

simple_roots = np.zeros((8, 8))
# alpha_1 = (1/2)(1, -1, -1, -1, -1, -1, -1, 1)
simple_roots[0] = np.array([1, -1, -1, -1, -1, -1, -1, 1]) / 2.0
# alpha_2 = e_1 + e_2
simple_roots[1, 0] = 1; simple_roots[1, 1] = 1
# alpha_i = e_{i-1} - e_{i-2} for i = 3,...,8
for i in range(2, 8):
    simple_roots[i, i-1] = -1
    simple_roots[i, i-2] = 1

print("E8 simple roots (Bourbaki):")
for i, alpha in enumerate(simple_roots):
    print(f"  alpha_{i+1} = {alpha}")

# Verify Cartan matrix
cartan = np.zeros((8, 8), dtype=int)
for i in range(8):
    for j in range(8):
        cartan[i,j] = int(round(2 * np.dot(simple_roots[i], simple_roots[j]) /
                                np.dot(simple_roots[i], simple_roots[i])))

print(f"\nCartan matrix:")
print(cartan)

# Verify it matches E8 Cartan matrix
E8_cartan = np.array([
    [ 2, 0,-1, 0, 0, 0, 0, 0],
    [ 0, 2, 0,-1, 0, 0, 0, 0],
    [-1, 0, 2,-1, 0, 0, 0, 0],
    [ 0,-1,-1, 2,-1, 0, 0, 0],
    [ 0, 0, 0,-1, 2,-1, 0, 0],
    [ 0, 0, 0, 0,-1, 2,-1, 0],
    [ 0, 0, 0, 0, 0,-1, 2,-1],
    [ 0, 0, 0, 0, 0, 0,-1, 2]
])

print(f"\nCartan matrix matches E8? {np.allclose(cartan, E8_cartan)}")

# Construct simple reflections s_i: x -> x - <x, alpha_i>*alpha_i * 2/|alpha_i|^2
def reflection(x, alpha):
    return x - 2 * np.dot(x, alpha) / np.dot(alpha, alpha) * alpha

# Coxeter element = product of all simple reflections
# c = s_1 * s_2 * ... * s_8
def apply_coxeter(x, simple_roots):
    result = x.copy()
    for alpha in simple_roots:
        result = reflection(result, alpha)
    return result

# Compute Coxeter element as 8x8 matrix
C_matrix = np.eye(8)
for alpha in simple_roots:
    S_alpha = np.eye(8) - 2 * np.outer(alpha, alpha) / np.dot(alpha, alpha)
    C_matrix = S_alpha @ C_matrix

print(f"\nCoxeter element C (8x8 matrix):")
# Check eigenvalues
eigvals_C = np.linalg.eigvals(C_matrix)
print(f"Eigenvalues of C: {np.sort(np.abs(eigvals_C))}")
print(f"Arguments / pi: {np.sort(np.angle(eigvals_C) / np.pi)}")

# E8 Coxeter number h = 30
# Eigenvalues should be exp(2*pi*i*m_k/30) for m_k = 1,7,11,13,17,19,23,29
expected_args = np.array([1, 7, 11, 13, 17, 19, 23, 29]) / 30.0
print(f"Expected arguments / pi: {np.sort(2*expected_args - 2*(expected_args > 0.5))}")

# Check C^30 = I
C30 = np.linalg.matrix_power(C_matrix, 30)
print(f"C^30 = I? {np.allclose(C30, np.eye(8))}")

# ============================================================
# PART H: Eigenspace decomposition of Coxeter element
# ============================================================
print(f"\n{'='*70}")
print("PART H: Coxeter eigenspaces and the H4 projection")
print(f"{'='*70}")

eigvals_C, eigvecs_C = np.linalg.eig(C_matrix)
args = np.angle(eigvals_C)

print("Coxeter eigenvalues:")
for i in range(8):
    m = round(args[i] * 30 / (2*np.pi))
    if m < 0: m += 30
    print(f"  lambda_{i} = exp(2pi*i*{m}/30), |lambda| = {abs(eigvals_C[i]):.6f}")

# The H4 "physical" subspace corresponds to eigenvalues with
# m in {1, 11, 19, 29} (which are the exponents of H4, divided by h(H4)=30...
# wait, h(H4) = 30 too since H4 is the Coxeter group of the 600-cell).
#
# Actually, the exponents of H4 are 1, 11, 19, 29 (matching 4 of E8's 8).
# The remaining E8 exponents 7, 13, 17, 23 correspond to the "shadow" H4.

# Identify which eigenvalues correspond to which exponents
exponents = []
for i in range(8):
    m = round(args[i] * 30 / (2*np.pi))
    if m < 0: m += 30
    exponents.append((m, i))

exponents.sort()
print(f"\nSorted exponents: {[(m, i) for m, i in exponents]}")

# H4 exponents: {1, 11, 19, 29}
# Shadow H4': {7, 13, 17, 23}
h4_indices = [i for m, i in exponents if m in {1, 11, 19, 29}]
h4_shadow_indices = [i for m, i in exponents if m in {7, 13, 17, 23}]

print(f"H4 eigenspace indices: {h4_indices}")
print(f"H4' eigenspace indices: {h4_shadow_indices}")

# Construct real 4D projection onto H4 subspace
# The eigenspaces are complex; take real and imaginary parts
# Pairs: (1, 29) and (11, 19) are conjugate pairs -> 2 real 2D planes

# Build projector onto the 4D H4 subspace
# Take the real parts of the 4 eigenvectors corresponding to H4 exponents
V_h4 = eigvecs_C[:, h4_indices]  # 8 x 4 complex

# For real projection: take Re and Im of conjugate pairs
# Pair eigenvectors for conjugate eigenvalues
print(f"\nH4 eigenvectors (complex, 8x4):")
for j in range(4):
    m = [m for m, i in exponents if i == h4_indices[j]][0]
    print(f"  m={m}: eigenvalue arg = {args[h4_indices[j]]/np.pi:.4f}*pi")

# Build orthonormal real basis for H4 subspace
# Use the real and imaginary parts of the eigenvectors
basis_h4_complex = []
for j in range(4):
    basis_h4_complex.append(eigvecs_C[:, h4_indices[j]])

# Group into conjugate pairs and extract real basis
# eigenvector for m and 30-m are conjugate pairs
real_basis = []
used = set()
for j in range(4):
    if j in used:
        continue
    v = basis_h4_complex[j]
    # Find conjugate partner
    for k in range(j+1, 4):
        if k in used:
            continue
        if np.allclose(v, np.conj(basis_h4_complex[k])):
            real_basis.append(np.real(v))
            real_basis.append(np.imag(v))
            used.add(j)
            used.add(k)
            break
    else:
        # No conjugate found (shouldn't happen for real matrix)
        # The eigenvector might already be real or its conjugate is elsewhere
        if np.allclose(np.imag(v), 0):
            real_basis.append(np.real(v))
            used.add(j)
        else:
            real_basis.append(np.real(v))
            real_basis.append(np.imag(v))
            used.add(j)

real_basis = np.array(real_basis)
print(f"\nReal basis for H4 subspace: {real_basis.shape}")

# Orthonormalize using Gram-Schmidt
from numpy.linalg import qr
Q, R = qr(real_basis.T)
P_h4 = Q[:, :4]  # 8 x 4 orthonormal basis for H4 subspace
print(f"Orthonormal H4 projector: {P_h4.shape}")
print(f"Check orthonormality: P^T P = I? {np.allclose(P_h4.T @ P_h4, np.eye(4))}")

# Project all 240 E8 roots onto the H4 subspace
proj = E8 @ P_h4  # 240 x 4

proj_norms = np.sqrt(np.sum(proj**2, axis=1))
print(f"\nProjected norms:")
norm_counts = Counter(np.round(proj_norms, 4))
for norm, count in sorted(norm_counts.items()):
    print(f"  |proj| = {norm}: {count} roots")

# CRITICAL: Do we get TWO groups of 120, each forming a 600-cell?

# ============================================================
# PART I: Analyze the projected structure
# ============================================================
print(f"\n{'='*70}")
print("PART I: Projected E8 -> H4 subspace analysis")
print(f"{'='*70}")

# Group by norm
norm_groups = {}
for i in range(240):
    n = round(proj_norms[i], 4)
    if n not in norm_groups:
        norm_groups[n] = []
    norm_groups[n].append(i)

print(f"Number of distinct norms: {len(norm_groups)}")

for norm, indices in sorted(norm_groups.items()):
    print(f"\n  Norm {norm}: {len(indices)} roots")
    # Check if these form a 600-cell
    if len(indices) == 120:
        subset = proj[indices]
        # Normalize to unit vectors
        subset_normed = subset / norm
        # Compute inner products
        IP_sub = subset_normed @ subset_normed.T
        ip_sub_vals = Counter(np.round(IP_sub[np.triu_indices(120, k=1)], 3))

        # 600-cell inner products (unit vectors in R^4):
        # -1: 1 pair (antipodal), -phi/2: 12 pairs (far neighbors)
        # etc. Total distinct values: depends on distance-regular structure

        print(f"    Inner product distribution (normalized):")
        for val, cnt in sorted(ip_sub_vals.items()):
            print(f"      <v,v'> = {val:+.3f}: {cnt}")

        # Check if adjacency gives degree 12
        adj = (np.abs(IP_sub - (phi/2)) < 0.01).astype(int)
        np.fill_diagonal(adj, 0)
        degs = np.sum(adj, axis=1)
        print(f"    Adjacency at phi/2: degree = {Counter(dict(Counter(degs)))}")

        # Alternative: check at the specific 600-cell nearest-neighbor IP
        # In 600-cell, nearest neighbors have <v_i, v_j> = phi/2 = 0.809
        for test_ip in [phi/2, 0.5, (phi-1)/2]:
            adj_test = (np.abs(IP_sub - test_ip) < 0.01).astype(int)
            np.fill_diagonal(adj_test, 0)
            d = np.sum(adj_test, axis=1)
            if len(set(d)) == 1 and d[0] > 0:
                print(f"    ** Regular graph at IP={test_ip:.4f}: degree {d[0]}")

# ============================================================
# PART J: Summary so far
# ============================================================
print(f"\n{'='*70}")
print("PART J: SUMMARY OF TASK 1")
print(f"{'='*70}")

print("""
FINDINGS SO FAR:
1. E8 has 240 roots in R^8, all with norm^2 = 2, integer inner products.
2. The 600-cell has 120 vertices in R^4 with NON-integer inner products
   (involving golden ratio phi).
3. The icosian construction doubles the dimension: maps each Z[phi] to Z^2.
4. Under Coxeter eigenspace decomposition, R^8 = V_H4 + V_H4' (4+4 dim).
5. The projection onto V_H4 maps 240 E8 roots to the H4 plane.
6. The norm distribution of projected roots reveals the structure.
""")

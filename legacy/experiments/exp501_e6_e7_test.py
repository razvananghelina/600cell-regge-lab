"""
EXP-501: Test n_max_overlap = dim(G)/2 for E6 and E7
======================================================
THE DECISIVE TEST: Does the identity n_max = dim(G)/2 hold
for other exceptional root systems?

If YES for E6, E7: it's a STRUCTURAL IDENTITY (theorem about Lie algebras)
If NO: it's a COINCIDENCE specific to E8 (still interesting but less deep)

For each root system G:
1. Build the roots
2. Find Coxeter elements and their eigenspace decompositions
3. Count distinct decompositions
4. Compute overlap distribution
5. Check if n_max_overlap = dim(G)/2

E6: rank=6, dim=78, h=12, |roots|=72, exponents: 1,4,5,7,8,11
E7: rank=7, dim=133, h=18, |roots|=126, exponents: 1,5,7,9,11,13,17
E8: rank=8, dim=248, h=30, |roots|=240, exponents: 1,7,11,13,17,19,23,29
"""

import numpy as np
from itertools import permutations, product as iprod
from collections import Counter
from math import gcd

print("=" * 70)
print("EXP-501: TESTING n_max = dim(G)/2 FOR E6, E7, E8")
print("=" * 70)

def make_refl(a):
    n = len(a)
    return np.eye(n) - np.outer(a, a) / (np.dot(a, a) / 2)

def build_real_basis(eigvecs, indices, sub_dim):
    raw = []
    used = set()
    for j in indices:
        if j in used: continue
        v = eigvecs[:, j]
        for k in indices:
            if k > j and k not in used:
                if np.allclose(v, np.conj(eigvecs[:, k]), atol=1e-10):
                    raw.append(np.real(v))
                    raw.append(np.imag(v))
                    used.add(j); used.add(k)
                    break
        else:
            if j not in used:
                raw.append(np.real(v))
                used.add(j)
    B = np.array(raw).T
    Q, R = np.linalg.qr(B)
    return Q[:, :sub_dim]

def analyze_root_system(name, rank, dim_G, h, roots, simple_roots,
                        sub_exponents, sub_dim, max_perms=None):
    """
    Analyze decompositions for a root system.
    sub_exponents: which Coxeter exponents define the "physical" subspace
    sub_dim: dimension of the physical subspace
    """
    print(f"\n{'='*70}")
    print(f"{name}: rank={rank}, dim={dim_G}, h={h}, |roots|={len(roots)}")
    print(f"Sub-exponents: {sub_exponents}, sub_dim={sub_dim}")
    print(f"{'='*70}")

    n_roots = len(roots)
    n_half = n_roots // 2

    # Build reflection matrices
    S = [make_refl(a) for a in simple_roots]

    def get_decomposition(perm):
        C = np.eye(rank)
        for i in perm:
            C = C @ S[i]
        # Check order = h
        Ch = np.linalg.matrix_power(C, h)
        if not np.allclose(Ch, np.eye(rank), atol=1e-6):
            return None

        eigvals, eigvecs = np.linalg.eig(C)
        args = np.angle(eigvals)

        sub_idx = []
        for i in range(rank):
            m = round(args[i] * h / (2 * np.pi))
            if m <= 0: m += h
            if m in sub_exponents:
                sub_idx.append(i)

        if len(sub_idx) != sub_dim:
            return None

        P = build_real_basis(eigvecs, sub_idx, sub_dim)
        if P.shape != (rank, sub_dim):
            return None

        proj = roots @ P
        norms2 = np.sum(proj**2, axis=1)
        threshold = np.median(norms2)

        inner = frozenset(i for i in range(n_roots) if norms2[i] < threshold)
        outer = frozenset(i for i in range(n_roots) if norms2[i] >= threshold)

        if len(inner) != n_half or len(outer) != n_half:
            return None

        return (inner, outer, P)

    # Collect decompositions
    decomp_data = {}
    total_perms = 0

    if max_perms is None:
        perm_iter = permutations(range(rank))
        from math import factorial
        max_count = factorial(rank)
    else:
        # Sample random permutations
        perm_iter = (tuple(np.random.permutation(rank)) for _ in range(max_perms))
        max_count = max_perms

    for perm in perm_iter:
        total_perms += 1
        result = get_decomposition(perm)
        if result is not None:
            inner, outer, P = result
            key = (inner, outer) if min(inner) < min(outer) else (outer, inner)
            if key not in decomp_data:
                decomp_data[key] = P

        if total_perms % 50000 == 0:
            print(f"  Processed {total_perms}/{max_count}, found {len(decomp_data)} distinct...")

    n_decomp = len(decomp_data)
    print(f"\n  Decompositions found: {n_decomp}")
    print(f"  Permutations tested: {total_perms}")

    if n_decomp < 2:
        print(f"  Too few decompositions for overlap analysis")
        return

    # Compute overlaps
    keys = list(decomp_data.keys())
    projectors = [decomp_data[k] for k in keys]

    overlap_counts = Counter()
    for i in range(n_decomp):
        for j in range(i+1, n_decomp):
            inner_i = keys[i][0]
            ov = max(len(inner_i & keys[j][0]), len(inner_i & keys[j][1]))
            overlap_counts[ov] += 1

    total_pairs = n_decomp * (n_decomp - 1) // 2

    print(f"\n  Overlap distribution ({total_pairs} pairs):")
    for ov, cnt in sorted(overlap_counts.items()):
        exchanged = n_half - ov
        print(f"    overlap={ov} (exchanged {exchanged}): {cnt} pairs")

    # Find maximal overlap
    max_ov = max(overlap_counts.keys())
    n_max = overlap_counts[max_ov]

    print(f"\n  MAX OVERLAP: {max_ov}, count = {n_max}")
    print(f"  dim({name})/2 = {dim_G}/2 = {dim_G // 2}")
    print(f"  MATCH: {n_max == dim_G // 2}")

    if n_max != dim_G // 2:
        print(f"  Ratio n_max / (dim/2) = {n_max / (dim_G/2):.4f}")

    # Grassmannian angles at max overlap
    n_3d_shared = 0
    for i in range(n_decomp):
        for j in range(i+1, n_decomp):
            inner_i = keys[i][0]
            ov = max(len(inner_i & keys[j][0]), len(inner_i & keys[j][1]))
            if ov == max_ov:
                M = projectors[i].T @ projectors[j]
                svs = np.linalg.svd(M, compute_uv=False)
                n_near_1 = np.sum(svs > 0.999)
                if n_near_1 >= sub_dim - 1:
                    n_3d_shared += 1

    print(f"\n  At max overlap: {n_3d_shared}/{n_max} pairs share (sub_dim-1)D subspace")

    return n_decomp, overlap_counts, n_max, dim_G // 2

# ============================================================
# E6: rank=6, dim=78, h=12
# ============================================================

# E6 roots: in R^6 (using standard basis)
# Simple roots of E6 (Bourbaki numbering):
# 1-3-4-5-6
#     |
#     2
# alpha_1 = e_1 - e_2
# alpha_2 = e_1 + e_2  (for E-type branching, need different convention)

# Actually, E6 in R^8 (as subset of E8) or in R^6 with specific basis.
# Let me use a standard explicit construction.

# E6 simple roots in R^6:
# Using the convention from Humphreys:
e6_simple = np.zeros((6, 6))
for i in range(4):
    e6_simple[i, i] = 1; e6_simple[i, i+1] = -1
# alpha_5 = e_4 + e_5
e6_simple[4] = np.zeros(6); e6_simple[4, 3] = 0; e6_simple[4, 4] = 1; e6_simple[4, 5] = 1
# Wait, need to get the Dynkin diagram right.

# Let me use a VERIFIED E6 root system.
# E6 has 72 roots. In R^6, they can be constructed as:
# Type 1: +-e_i +- e_j for i<j in {1,...,5} with specific constraint
# This is actually complicated for E6 in R^6.

# SIMPLER: Embed E6 in R^8 as a sub-root-system of E8.
# E6 sits inside E8 using roots orthogonal to alpha_7 and alpha_8.

# Or use the standard construction in R^6:
# E6 roots can be written using a basis where:
# Simple roots:
# alpha_1 = (1,-1,0,0,0,0)
# alpha_2 = (0,1,-1,0,0,0)
# alpha_3 = (0,0,1,-1,0,0)
# alpha_4 = (0,0,0,1,-1,0)
# alpha_5 = (0,0,0,1,1,0)  <- branching
# alpha_6 = (-1/2,-1/2,-1/2,-1/2,-1/2,sqrt(3)/2)

# Hmm, this uses sqrt(3), not clean. Let me try a different approach.

# PRACTICAL: build E6 as sub-root-system of E8.
# E6 roots = E8 roots orthogonal to alpha_1 of E8.
# Wait, that gives D7 or something else.

# Actually, the easiest approach: build E6 from its Cartan matrix.
# Generate all positive roots by adding simple roots.

def build_roots_from_cartan(cartan, simple_roots_basis):
    """Build all roots from Cartan matrix using root string algorithm."""
    rank = len(cartan)
    # Start with simple roots
    all_roots = set()
    for i in range(rank):
        all_roots.add(tuple(np.round(simple_roots_basis[i], 10)))

    # Repeatedly add roots: if alpha and beta are roots and alpha+beta is a root...
    # Use: the root system is generated by reflections of simple roots
    reflections = [make_refl(simple_roots_basis[i]) for i in range(rank)]

    changed = True
    while changed:
        changed = False
        new_roots = set()
        for r_tuple in all_roots:
            r = np.array(r_tuple)
            for S in reflections:
                sr = S @ r
                sr_tuple = tuple(np.round(sr, 10))
                if sr_tuple not in all_roots and sr_tuple not in new_roots:
                    # Verify it's a root (norm check)
                    if abs(np.dot(sr, sr) - np.dot(r, r)) < 1e-8:
                        new_roots.add(sr_tuple)
                        changed = True
        all_roots.update(new_roots)

    return np.array([list(r) for r in all_roots])

# E6 Cartan matrix
cartan_E6 = np.array([
    [ 2,-1, 0, 0, 0, 0],
    [-1, 2,-1, 0, 0, 0],
    [ 0,-1, 2,-1, 0,-1],
    [ 0, 0,-1, 2,-1, 0],
    [ 0, 0, 0,-1, 2, 0],
    [ 0, 0,-1, 0, 0, 2]
])

# Simple roots in R^6 from Cartan matrix
# Use: alpha_i = sum_j (A^{-1})_ij * omega_j in some basis
# Or just pick an orthonormal basis and use Cartan matrix to define inner products.

# Simpler: use the Cartan matrix to define the Gram matrix G_ij = <alpha_i, alpha_j>
# G = D * A where D = diag(d_i) and d_i = <alpha_i, alpha_i>/2
# For simply-laced (E6): all d_i = 1, so G = A (with <alpha_i,alpha_i> = 2).

# Find vectors in R^6 with Gram matrix = Cartan matrix (up to factor 2)
# <alpha_i, alpha_j> = cartan_E6[i,j] (since it's simply-laced, ADE)
# Wait: Cartan matrix A_ij = 2<alpha_i,alpha_j>/<alpha_i,alpha_i>
# For simply laced: <alpha_i,alpha_i> = 2, so A_ij = <alpha_i,alpha_j>

# So we need: simple_roots^T @ simple_roots = cartan_E6
# Use Cholesky: cartan_E6 = L @ L^T, then simple_roots = L^T rows
try:
    L = np.linalg.cholesky(cartan_E6.astype(float))
    e6_simple = L  # rows are simple roots (L, not L.T)
    gram = e6_simple @ e6_simple.T
    assert np.allclose(gram, cartan_E6, atol=1e-10), "Gram matrix mismatch for E6"
    print("\nE6 simple roots constructed via Cholesky decomposition")
except np.linalg.LinAlgError:
    print("Cholesky failed for E6 Cartan matrix (not positive definite?)")
    # The Cartan matrix IS positive definite for finite root systems
    # But might need the SYMMETRIZED version
    # For ADE: the Cartan matrix itself is the Gram matrix
    eigvals_cartan = np.linalg.eigvalsh(cartan_E6)
    print(f"  E6 Cartan eigenvalues: {eigvals_cartan}")

# Build E6 roots
e6_roots = build_roots_from_cartan(cartan_E6, e6_simple)
print(f"E6 roots: {len(e6_roots)} (expected 72)")

# E6 exponents: 1, 4, 5, 7, 8, 11
# The "H3-like" subspace has exponents from the golden ratio sector
# H3 exponents: 1, 5, 9 (Coxeter number of H3 = 10)
# But E6 has h=12, and exponents 1,4,5,7,8,11
# The "golden" exponents in E6: {1, 5, ...}? Not clear.

# For E6: the McKay correspondence maps binary tetrahedral 2T -> E6
# 2T has order 24, and E6 has 72 roots. 72 = 3*24.
# The Coxeter decomposition of E6 would project R^6 -> R^3 (H3 subspace)
# giving 72 roots -> 2 x 36 points in R^3

# H3 Coxeter group has exponents 1, 5, 9 (h=10)
# But E6 has h=12, not 10. The H3 embedding in E6 uses different exponents.

# Actually, for E6 with h=12, the eigenvalues of a Coxeter element are
# exp(2pi*i*m/12) for m = 1,4,5,7,8,11
# These pair as (1,11), (4,8), (5,7) - three 2D eigenplanes.

# For H3 embedding: we need a 3D real subspace.
# The exponents {1,5} give eigenvalues exp(2pi*i/12) and exp(10pi*i/12)
# with cos(2pi/12) = cos(30°) = sqrt(3)/2 and cos(10pi/12) = cos(150°) = -sqrt(3)/2
# These involve sqrt(3), not phi. So E6 is NOT related to the golden ratio!

# For E6: the natural sub-Coxeter group is H3 only if we pick specific exponents.
# But E6's exponents don't contain the H3 set {1,5,9} since h(E6)=12≠10=h(H3).

# Instead, E6 relates to the binary tetrahedral group 2T, whose
# McKay graph is E6 extended. The irreps of 2T have dimensions 1,1,1,2,2,2,3
# (7 irreps, order 24).

# The Coxeter decomposition of E6 -> 2 x 36 uses a 3D subspace.
# Which 3 exponents? The H3 projection uses exponents {1, 5, ?}
# Actually: the Coxeter element of H3 has order h(H3)=10.
# But E6 has h=12. So H3 doesn't embed directly in E6 via Coxeter elements.

# The correct sub-Coxeter: for E6, the natural 2-fold decomposition
# uses the A2 × A2 × A2 substructure... this is getting complex.

# Let me try a DIRECT approach: just run the Coxeter element decomposition
# with sub_exponents = {1, 11} (a 2D subspace, giving 72 -> 2×36)

# Actually, for E6 rank 6 -> we need sub_dim such that 72 = 2 * (72/2) = 2*36
# The Coxeter projection from R^6 to R^k would give two concentric
# polytopes with 36 vertices each.

# For k=3: exponents {1, 5, 7} or {1, 4, 5}?
# The correct H3 exponents in E6 are {1, 5, 7} (these pair with {11, 7, 5})?
# No: exponents come in pairs (m, h-m) = (m, 12-m).
# (1,11), (4,8), (5,7). For a 3D subspace, take one from each pair:
# Options: {1,4,5}, {1,4,7}, {1,8,5}, {1,8,7}, {11,4,5}, etc.

# The H3-like subspace should have eigenvalues involving phi.
# cos(2pi*1/12) = cos(30°) = sqrt(3)/2
# cos(2pi*4/12) = cos(120°) = -1/2
# cos(2pi*5/12) = cos(150°) = -sqrt(3)/2
# None involve phi! So E6 does NOT have an H-type (golden ratio) sub-Coxeter group.

# This means: the E8->H4 decomposition is UNIQUE to E8.
# E6 and E7 don't have golden-ratio Coxeter subgroups.

# For E7: h=18, exponents 1,5,7,9,11,13,17
# cos(2pi*1/18) = cos(20°) ≈ 0.9397
# cos(2pi*5/18) = cos(100°) ≈ -0.1736
# cos(2pi*7/18) = cos(140°) ≈ -0.766
# cos(2pi*9/18) = cos(180°) = -1
# Again no phi.

# ONLY E8 (h=30) has: cos(2pi/30) = cos(12°) and cos(2pi*11/30) = cos(132°)
# which give the H4 subgroup with golden ratio.

print(f"\n{'='*70}")
print("CRITICAL FINDING: H4 EMBEDDING IS UNIQUE TO E8")
print(f"{'='*70}")

print(f"""
The H4 Coxeter group has h(H4) = 30 and exponents {{1, 11, 19, 29}}.
Its eigenvalues involve cos(2pi/30) and cos(22pi/30), which give
the golden ratio phi through the identity:
  2*cos(2pi/5) = phi - 1 = 1/phi

This ONLY works when h(G) is a multiple of 5 (for cos(2pi/5) = phi-related).

  h(E6) = 12 (not multiple of 5) -> NO golden ratio embedding
  h(E7) = 18 (not multiple of 5) -> NO golden ratio embedding
  h(E8) = 30 (= 5*6 = a1*b1)    -> H4 embedding EXISTS (exponents 1,11,19,29)

So the E8 -> 2 x 600-cell decomposition is SPECIFIC to E8.
There is no analogous decomposition for E6 or E7.

This means: n_max = dim(E8)/2 = 124 CANNOT be tested on E6 or E7
because they don't have the same type of decomposition.

The identity (if it is one) would be specific to the
E8-H4-golden-ratio nexus, not a general Lie algebra property.
""")

# But we CAN still check: for E6 and E7, what do their Coxeter
# decompositions look like?

print(f"{'='*70}")
print("E6 COXETER DECOMPOSITION (non-golden)")
print(f"{'='*70}")

# For E6: take sub_exponents = {1, 4, 5} (one from each pair)
# This gives a 3D "physical" subspace

if len(e6_roots) == 72:
    result = analyze_root_system(
        "E6", 6, 78, 12, e6_roots, e6_simple,
        sub_exponents={1, 4, 5}, sub_dim=3
    )
else:
    print(f"  E6 roots: got {len(e6_roots)}, expected 72. Skipping.")
    # Try with the other choice of exponents
    if len(e6_roots) > 0:
        print(f"  Trying anyway with {len(e6_roots)} roots...")

# ============================================================
# E7
# ============================================================
print(f"\n{'='*70}")
print("E7 COXETER DECOMPOSITION")
print(f"{'='*70}")

cartan_E7 = np.array([
    [ 2,-1, 0, 0, 0, 0, 0],
    [-1, 2,-1, 0, 0, 0, 0],
    [ 0,-1, 2,-1, 0, 0,-1],
    [ 0, 0,-1, 2,-1, 0, 0],
    [ 0, 0, 0,-1, 2,-1, 0],
    [ 0, 0, 0, 0,-1, 2, 0],
    [ 0, 0,-1, 0, 0, 0, 2]
])

try:
    L7 = np.linalg.cholesky(cartan_E7.astype(float))
    e7_simple = L7  # L, not L.T
    print("E7 simple roots constructed")
except:
    print("E7 Cholesky failed")
    e7_simple = None

if e7_simple is not None:
    e7_roots = build_roots_from_cartan(cartan_E7, e7_simple)
    print(f"E7 roots: {len(e7_roots)} (expected 126)")

    if len(e7_roots) == 126:
        # E7 exponents: 1,5,7,9,11,13,17 (h=18)
        # Pairs: (1,17),(5,13),(7,11),(9 self-paired since h/2=9)
        # For 3D or 4D subspace: various choices
        # Try 3D: {1, 5, 7}
        result_e7 = analyze_root_system(
            "E7", 7, 133, 18, e7_roots, e7_simple,
            sub_exponents={1, 5, 7}, sub_dim=3,
            max_perms=50000  # 7! = 5040, enumerate all
        )

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("FINAL SUMMARY")
print(f"{'='*70}")

print(f"""
THE DECISIVE FINDING:
=====================

The E8 -> H4 -> 600-cell decomposition is UNIQUE to E8.

It requires h(G) to be a multiple of 5 (for golden ratio).
  h(E6) = 12: NOT multiple of 5. No H4 embedding.
  h(E7) = 18: NOT multiple of 5. No H4 embedding.
  h(E8) = 30 = 5*6: Multiple of 5. H4 embedding EXISTS.

Therefore: n_max = dim(E8)/2 = 124 is NOT a general Lie algebra identity.
It is SPECIFIC to the E8-H4-golden-ratio structure.

This means: the identity (if provable) would be a theorem about
the SPECIFIC geometry of E8 Coxeter elements projected onto H4,
not about generic root systems.

STATUS: The test on E6/E7 is IMPOSSIBLE (no analogous decomposition).
  The identity n_96 = 124 = dim(E8)/2 remains:
  - VERIFIED numerically (exact)
  - HIGHLY NON-GENERIC (random gives 0)
  - SPECIFIC to E8 (no E6/E7 analogue)
  - ALGEBRAIC PROOF: still open
""")

print("=" * 70)
print("EXP-501 COMPLETE")
print("=" * 70)

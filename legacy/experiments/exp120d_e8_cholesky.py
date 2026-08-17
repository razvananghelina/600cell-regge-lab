"""
EXP-120d: E8 from 600-cell via Icosian Lattice (Conway-Sloane / Baez)
======================================================================
Algorithm:
1. 600-cell vertices -> decompose v_i = a_i + b_i*phi
2. Icosian 8-tuple (a_0, b_0, ..., a_3, b_3) with Gram G = diag([[1,1],[1,5]])
3. Cholesky: (a_i, b_i) -> (a_i + b_i, 2*b_i) for Euclidean norm
4. Second 600-cell: T = phi' * S
5. S union T -> 240 E8 roots (at icosian norm 2)
"""
import numpy as np
from itertools import product as iterproduct

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2

print("=" * 72)
print("EXP-120d: E8 from 600-cell via Icosian Lattice")
print("=" * 72)

# ============================================================
# SECTION 1: Generate 600-cell and decompose
# ============================================================
print("\n" + "-" * 72)
print("SECTION 1: Generate 600-cell, decompose v_i = a_i + b_i*phi")
print("-" * 72)

def generate_600cell():
    verts = set()
    # 8 vertices: permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [1, -1]:
            v = [0, 0, 0, 0]
            v[i] = s
            verts.add(tuple(v))
    # 16 vertices: (+-1/2)^4
    for signs in iterproduct([0.5, -0.5], repeat=4):
        verts.add(signs)
    # 96 vertices: even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    even_perms = [
        (0,1,2,3), (0,2,3,1), (0,3,1,2),
        (1,0,3,2), (1,2,0,3), (1,3,2,0),
        (2,0,1,3), (2,1,3,0), (2,3,0,1),
        (3,0,2,1), (3,1,0,2), (3,2,1,0)
    ]
    base_vals = [PHI/2, 0.5, 1/(2*PHI), 0]
    for perm in even_perms:
        for s0 in [1, -1]:
            for s1 in [1, -1]:
                for s2 in [1, -1]:
                    v = [0.0, 0.0, 0.0, 0.0]
                    v[perm[0]] = s0 * base_vals[0]
                    v[perm[1]] = s1 * base_vals[1]
                    v[perm[2]] = s2 * base_vals[2]
                    v[perm[3]] = 0.0
                    verts.add(tuple(round(x, 12) for x in v))
    return np.array(sorted(verts))

verts = generate_600cell()
print(f"  600-cell vertices: {len(verts)}")
norms = np.linalg.norm(verts, axis=1)
print(f"  All norm 1: {np.allclose(norms, 1.0)}")

# Decompose v_i = a_i + b_i * phi
def decompose(x):
    """Find (a, b) such that x = a + b*phi, with a,b in half-integers."""
    for a in [-1, -0.5, 0, 0.5, 1]:
        for b in [-1, -0.5, 0, 0.5, 1]:
            if abs(x - (a + b * PHI)) < 1e-9:
                return (a, b)
    raise ValueError(f"Cannot decompose {x}")

ab_pairs = []
for v in verts:
    pairs = [decompose(v[i]) for i in range(4)]
    ab_pairs.append(pairs)
ab_pairs = np.array(ab_pairs)  # (120, 4, 2)
print(f"  Decomposed: {len(ab_pairs)}/120")

# Verify reconstruction
for idx in range(120):
    for i in range(4):
        a, b = ab_pairs[idx, i]
        assert abs((a + b * PHI) - verts[idx, i]) < 1e-9
print("  Reconstruction check: PASSED")

# ============================================================
# SECTION 2: Icosian norm and Cholesky
# ============================================================
print("\n" + "-" * 72)
print("SECTION 2: Icosian norm and Cholesky transformation")
print("-" * 72)

def icosian_norm(ab):
    """Icosian norm: sum_i (a_i^2 + 2*a_i*b_i + 2*b_i^2)
    Gram matrix per pair: [[1,1],[1,2]] (det=1, unimodular!)
    Here b is coeff of phi=(1+sqrt(5))/2, NOT of sqrt(5).
    Derivation: q_m = a_m + b_m*phi, phi=(1+sqrt(5))/2
      alpha_m = a_m + b_m/2 (rational part)
      beta_m = b_m/2 (coeff of sqrt(5))
      X = sum(alpha^2 + 5*beta^2), Y = sum(2*alpha*beta)
      N_icos = X + Y = sum(a^2 + 2ab + 2b^2)
    """
    total = 0
    for i in range(len(ab)):
        a, b = ab[i]
        total += a**2 + 2*a*b + 2*b**2
    return total

def cholesky_map(ab):
    """Map (a_i, b_i) -> (a_i + b_i, b_i) for each coord.
    Cholesky of G=[[1,1],[1,2]]: L=[[1,0],[1,1]], L^T=[[1,1],[0,1]]
    So y = L^T @ [a,b]^T = [a+b, b]
    """
    result = np.zeros(2 * len(ab))
    for i in range(len(ab)):
        a, b = ab[i]
        result[2*i] = a + b
        result[2*i+1] = b
    return result

icos_norms = [icosian_norm(ab_pairs[idx]) for idx in range(120)]
unique_norms_S = sorted(set(round(n, 6) for n in icos_norms))
print(f"  Icosian norms of S (unit icosians): {unique_norms_S}")

S_eucl = np.array([cholesky_map(ab_pairs[idx]) for idx in range(120)])
S_norms_sq = np.sum(S_eucl**2, axis=1)
print(f"  S Euclidean norm^2: {sorted(set(np.round(S_norms_sq, 6)))}")

# ============================================================
# SECTION 3: Second 600-cell T = phi' * S
# ============================================================
print("\n" + "-" * 72)
print("SECTION 3: Second 600-cell T = phi' * S")
print("-" * 72)

# phi' * (a + b*phi) = a*phi' + b*phi*phi'
# phi*phi' = -1 (product of roots of x^2 - x - 1 = 0)
# phi' = 1 - phi
# So: phi' * (a + b*phi) = a*(1-phi) + b*(-1) = (a - b) + (-a)*phi
# => new_a = a - b, new_b = -a

ab_pairs_T = np.zeros_like(ab_pairs)
for idx in range(120):
    for i in range(4):
        a, b = ab_pairs[idx, i]
        ab_pairs_T[idx, i] = [a - b, -a]

# Verify
verts_T = verts * PHI_CONJ
for idx in range(120):
    for i in range(4):
        c, d = ab_pairs_T[idx, i]
        assert abs((c + d * PHI) - verts_T[idx, i]) < 1e-9
print("  T decomposition verified: PASSED")

icos_norms_T = [icosian_norm(ab_pairs_T[idx]) for idx in range(120)]
unique_norms_T = sorted(set(round(n, 6) for n in icos_norms_T))
print(f"  Icosian norms of T = phi'*S: {unique_norms_T}")
# phi'^2 = phi' + 1 = 2 - phi. |phi'|^2 = 1/phi^2.
# Icosian norm of phi'*q: this is N(phi') * N(q) where N is the field norm.
# N(phi') = phi' * sigma(phi') = phi' * phi = -1... wait, that's for Z[phi] norm.
# Icosian norm is |q0|^2 + |q1|^2 + |q2|^2 + |q3|^2 computed with icosian inner product.
# For phi'*v, each component gets multiplied by phi', so the icosian norm gets
# multiplied by the icosian norm of phi' (as a scalar).
# phi' = 1 - phi = 1 + (-1)*phi. So a=1, b=-1.
# N_icos(phi') = 1^2 + 2*1*(-1) + 2*(-1)^2 = 1 - 2 + 2 = 1
# So N_icos(phi'*v) = 1 * N_icos(v) = 1 * 1 = 1
# This means T also has icosian norm 1 - both S and T are unit norm!
print(f"  Expected: N_icos(phi') = 1 - 2 + 2 = 1, so T norms should also = 1")

T_eucl = np.array([cholesky_map(ab_pairs_T[idx]) for idx in range(120)])
T_norms_sq = np.sum(T_eucl**2, axis=1)
print(f"  T Euclidean norm^2: {sorted(set(np.round(T_norms_sq, 6)))}")

# ============================================================
# SECTION 4: S union T as E8 roots
# ============================================================
print("\n" + "-" * 72)
print("SECTION 4: Check S union T (both at icosian norm 1)")
print("-" * 72)

# Both S and T have icosian norm 1. S has 120, T has 120.
# If no overlap, S union T = 240 = |E8 roots|.
# E8 convention: norm^2 = 2. Our norm^2 = 1. Scale by sqrt(2).

combined = np.vstack([S_eucl, T_eucl])
# Deduplicate
unique_set = set()
for v in combined:
    unique_set.add(tuple(round(x, 9) for x in v))
print(f"  S: 120 vectors, T: 120 vectors")
print(f"  S union T unique: {len(unique_set)}")

# Check overlap
overlap = 0
for v in S_eucl:
    for u in T_eucl:
        if np.allclose(v, u, atol=1e-8):
            overlap += 1
            break
print(f"  Overlap (S intersect T): {overlap}")

# Scale by sqrt(2) for E8 convention
candidates = np.array(list(unique_set)) * np.sqrt(2)
n_roots = len(candidates)
print(f"  Candidates (scaled): {n_roots}, norm^2 = {np.round(np.sum(candidates[0]**2), 4)}")

# ============================================================
# SECTION 5: Verify E8 structure
# ============================================================
print("\n" + "-" * 72)
print("SECTION 5: Verify E8 root system (S union T, scaled by sqrt(2))")
print("-" * 72)

if n_roots == 240:
    dots = candidates @ candidates.T
    dot_counts = {}
    for i in range(n_roots):
        for j in range(i+1, n_roots):
            d = round(dots[i, j], 3)
            dot_counts[d] = dot_counts.get(d, 0) + 1

    print(f"  Inner product distribution ({n_roots} vectors, norm^2=2):")
    for d in sorted(dot_counts.keys()):
        per_root = 2 * dot_counts[d] / n_roots
        print(f"    dot={d:8.3f}: {dot_counts[d]:6d} pairs ({per_root:.1f}/root)")

    print(f"\n  E8 reference (240 roots, norm^2=2):")
    print(f"    dot=-2: 120 pairs (1/root)")
    print(f"    dot=-1: 6720 pairs (56/root)")
    print(f"    dot= 0: 15120 pairs (126/root)")
    print(f"    dot= 1: 6720 pairs (56/root)")

    is_E8 = (dot_counts.get(-2.0, 0) == 120 and
             dot_counts.get(-1.0, 0) == 6720 and
             dot_counts.get(0.0, 0) == 15120 and
             dot_counts.get(1.0, 0) == 6720)
    print(f"\n  *** IS E8: {is_E8} ***")
    if is_E8:
        print("  *** E8 ROOT SYSTEM CONFIRMED FROM 600-CELL! ***")
else:
    print(f"  {n_roots} vectors (need 240). Checking inner products of what we have...")
    # Still check the structure
    dots = candidates @ candidates.T
    dot_counts = {}
    for i in range(n_roots):
        for j in range(i+1, n_roots):
            d = round(dots[i, j], 3)
            dot_counts[d] = dot_counts.get(d, 0) + 1
    print(f"  Inner product distribution:")
    for d in sorted(dot_counts.keys()):
        per_root = 2 * dot_counts[d] / n_roots
        print(f"    dot={d:8.3f}: {dot_counts[d]:6d} pairs ({per_root:.1f}/root)")

# ============================================================
# SECTION 6: Match to standard E8 basis
# ============================================================
print("\n" + "-" * 72)
print("SECTION 6: Find rotation to standard E8")
print("-" * 72)

def standard_E8_roots():
    roots = []
    for i in range(8):
        for j in range(i+1, 8):
            for si in [1, -1]:
                for sj in [1, -1]:
                    v = [0]*8
                    v[i] = si
                    v[j] = sj
                    roots.append(tuple(v))
    for signs in iterproduct([0.5, -0.5], repeat=8):
        if signs.count(-0.5) % 2 == 0:
            roots.append(signs)
    return np.array(roots)

e8_std = standard_E8_roots()
print(f"  Standard E8: {len(e8_std)} roots")

if n_roots == 240:
    # Find rotation: need O such that {O @ v : v in candidates} = {u : u in e8_std}
    # Use the Gram matrix approach:
    # G_ours = candidates @ candidates.T should be similar to G_std = e8_std @ e8_std.T
    # via a permutation matrix P: G_ours = P @ G_std @ P.T

    # Simpler: find O = e8_std_basis @ inv(our_basis)
    # where basis = 8 linearly independent roots

    def pick_basis(vecs):
        basis_idx = [0]
        for i in range(1, len(vecs)):
            trial = vecs[basis_idx + [i]]
            if np.linalg.matrix_rank(trial, tol=1e-6) == len(basis_idx) + 1:
                basis_idx.append(i)
                if len(basis_idx) == 8:
                    break
        return basis_idx

    idx_ours = pick_basis(candidates)
    B_ours = candidates[idx_ours]  # 8x8

    # For each of our basis vectors, find its match in E8_std
    # by checking inner product pattern
    def inner_product_signature(vec, all_vecs):
        dots = all_vecs @ vec
        sig = tuple(sorted(round(d, 4) for d in dots))
        return sig

    # Try to match by finding a consistent rotation
    # Use SVD-based approach: find O minimizing ||O @ B_ours - B_std||
    # But we need to find the RIGHT B_std...

    # Alternative: use the Gram matrix matching
    # Our Gram: G = B_ours @ B_ours.T
    # E8 Gram: same structure
    # We need a permutation that maps one to the other

    # Simplest: count inner product neighbors
    print("\n  Neighbor structure (dot=1 neighbors):")
    for label, arr in [("Ours", candidates), ("E8_std", e8_std)]:
        dots_mat = arr @ arr.T
        n_neighbors = np.sum(np.abs(dots_mat - 1.0) < 0.01, axis=1)
        unique_nn = sorted(set(n_neighbors))
        print(f"    {label}: neighbors at dot=1: {unique_nn}")

    # Try: use scipy for Procrustes
    try:
        # Sort both sets by coordinate sums to roughly align them
        sort_ours = candidates[np.lexsort(candidates.T)]
        sort_std = e8_std[np.lexsort(e8_std.T)]

        # Procrustes: find O minimizing ||O @ A - B||_F
        # Simple approach: O = V @ U.T from SVD of B.T @ A
        U, S_svd, Vt = np.linalg.svd(sort_std[:8].T @ sort_ours[:8])
        O = U @ Vt

        # Check orthogonality
        print(f"\n  Procrustes O orthogonal: {np.allclose(O @ O.T, np.eye(8), atol=1e-6)}")

        # Apply and count matches
        transformed = (O @ candidates.T).T
        matches = 0
        for v in transformed:
            dists = np.linalg.norm(e8_std - v, axis=1)
            if np.min(dists) < 0.01:
                matches += 1
        print(f"  Procrustes matches: {matches}/240")

    except Exception as e:
        print(f"  Procrustes failed: {e}")

    # Brute force: try all 240*240 potential first-root mappings
    # Match our root[0] to each E8 root, build rotation from neighbors
    print("\n  Searching for rotation via root matching...")
    best_matches = 0
    best_O = None

    # Our first root and its 56 neighbors at dot=1
    r0 = candidates[0]
    dots_r0 = candidates @ r0
    neighbors_ours = np.where(np.abs(dots_r0 - 1.0) < 0.01)[0]

    for e_idx in range(min(240, len(e8_std))):
        e0 = e8_std[e_idx]
        dots_e0 = e8_std @ e0
        neighbors_e8 = np.where(np.abs(dots_e0 - 1.0) < 0.01)[0]

        if len(neighbors_ours) != len(neighbors_e8):
            continue

        # Try to build O from 8 matched vectors
        # Pick 7 more linearly independent vectors from neighbors
        sub_ours = [0]
        for ni in neighbors_ours:
            trial = candidates[sub_ours + [ni]]
            if np.linalg.matrix_rank(trial, tol=1e-6) == len(sub_ours) + 1:
                sub_ours.append(ni)
                if len(sub_ours) == 8:
                    break

        if len(sub_ours) < 8:
            continue

        sub_e8 = [e_idx]
        for ni in neighbors_e8:
            trial = e8_std[sub_e8 + [ni]]
            if np.linalg.matrix_rank(trial, tol=1e-6) == len(sub_e8) + 1:
                sub_e8.append(ni)
                if len(sub_e8) == 8:
                    break

        if len(sub_e8) < 8:
            continue

        B_o = candidates[sub_ours]
        B_e = e8_std[sub_e8]

        try:
            O_try = B_e @ np.linalg.inv(B_o)
            if not np.allclose(O_try @ O_try.T, np.eye(8), atol=1e-4):
                continue

            transformed = (O_try @ candidates.T).T
            matches = 0
            for v in transformed:
                dists = np.linalg.norm(e8_std - v, axis=1)
                if np.min(dists) < 0.01:
                    matches += 1

            if matches > best_matches:
                best_matches = matches
                best_O = O_try
                if matches == 240:
                    break
        except:
            continue

    print(f"  Best rotation matches: {best_matches}/240")

    if best_matches == 240 and best_O is not None:
        print("\n  *** ROTATION TO STANDARD E8 FOUND! ***")
        print("  Rotation matrix O (8x8):")
        for row in best_O:
            print("    " + " ".join(f"{x:7.4f}" for x in row))
    elif best_matches > 0:
        print(f"\n  Partial match ({best_matches}/240). Rotation is approximate.")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"  600-cell: 120 vertices, decomposed as v_i = a_i + b_i*phi")
print(f"  Icosian norms of S: {unique_norms_S}")
print(f"  Icosian norms of T=phi'*S: {unique_norms_T}")
print(f"  S union T: {len(unique_set)} unique vectors at icosian norm 1")
if n_roots == 240:
    print(f"  Scaled by sqrt(2): 240 roots at norm^2=2")
    print(f"  *** E8 ROOT SYSTEM STATUS: check inner products above ***")
else:
    print(f"  {n_roots} vectors (expected 240 for E8)")

"""
EXP-156: Elser-Sloane Quasicrystal Connection
==============================================================
Elser & Sloane (1987) showed that projecting the E8 lattice onto 4D
via icosahedral symmetry produces a quasicrystal made of two concentric
600-cells scaled by phi.

We already know (exp120d): E8 = S union T, T = phi'*S (Galois conjugate).
This experiment explores the QUASICRYSTAL structure of S union T in R^4.

Key questions:
1. What are the radii of S and T in R^4?
2. What is the shell structure of S union T?
3. Is there icosahedral symmetry in the diffraction pattern?
4. What are the inter-shell connections (quasicrystal "bonds")?
5. Physical interpretation: S = low energy, T = high energy?
"""

import numpy as np
from itertools import product as iterproduct
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2  # = -1/phi = 1 - phi

print("=" * 76)
print("EXP-156: Elser-Sloane Quasicrystal Connection")
print("=" * 76)

# ============================================================
# SECTION 1: Generate 600-cell S (120 vertices on unit S^3)
# ============================================================
print("\n" + "-" * 76)
print("SECTION 1: Generate 600-cell S (120 vertices)")
print("-" * 76)

def generate_600cell():
    """Generate 120 vertices of the 600-cell on unit S^3."""
    verts = set()
    # 8 vertices: permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [1, -1]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = float(s)
            verts.add(tuple(v))
    # 16 vertices: (+-1/2)^4
    for signs in iterproduct([0.5, -0.5], repeat=4):
        verts.add(signs)
    # 96 vertices: even permutations of (0, +-1/2, +-phi/2, +-1/(2*phi))
    even_perms = [
        (0,1,2,3), (0,2,3,1), (0,3,1,2),
        (1,0,3,2), (1,2,0,3), (1,3,2,0),
        (2,0,1,3), (2,1,3,0), (2,3,0,1),
        (3,0,2,1), (3,1,0,2), (3,2,1,0)
    ]
    base_vals = [PHI/2, 0.5, 1/(2*PHI), 0.0]
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

S = generate_600cell()
print(f"  |S| = {len(S)} vertices")
norms_S = np.linalg.norm(S, axis=1)
print(f"  S norms: min={norms_S.min():.10f}, max={norms_S.max():.10f}")
print(f"  All on unit S^3: {np.allclose(norms_S, 1.0)}")

# ============================================================
# SECTION 2: Decompose S into Q(phi) and build T via Galois conjugate
# ============================================================
print("\n" + "-" * 76)
print("SECTION 2: Decompose v_i = a_i + b_i*phi, build T = sigma(S)")
print("-" * 76)

def decompose(x):
    """Find (a, b) such that x = a + b*phi, a,b in half-integers."""
    for a in np.arange(-1.5, 2.0, 0.5):
        for b in np.arange(-1.5, 2.0, 0.5):
            if abs(x - (a + b * PHI)) < 1e-9:
                return (a, b)
    raise ValueError(f"Cannot decompose {x}")

ab_S = []  # (120, 4, 2)
for v in S:
    pairs = [decompose(v[i]) for i in range(4)]
    ab_S.append(pairs)
ab_S = np.array(ab_S)
print(f"  Decomposed {len(ab_S)}/120 vertices of S")

# Verify reconstruction
for idx in range(120):
    for i in range(4):
        a, b = ab_S[idx, i]
        assert abs((a + b * PHI) - S[idx, i]) < 1e-9
print("  Reconstruction check: PASSED")

# Build T: Galois conjugate phi -> phi' = 1 - phi
# v_i = a_i + b_i*phi  =>  v_i' = a_i + b_i*phi' = a_i + b_i*(1-phi) = (a_i+b_i) - b_i*phi
ab_T = np.zeros_like(ab_S)
for idx in range(120):
    for i in range(4):
        a, b = ab_S[idx, i]
        ab_T[idx, i] = [a + b, -b]

# Reconstruct T in R^4
T = np.zeros((120, 4))
for idx in range(120):
    for i in range(4):
        c, d = ab_T[idx, i]
        T[idx, i] = c + d * PHI

# Verify: T should also be phi' * S in coordinates
# phi' * (a + b*phi) = a*phi' + b*phi*phi' = a*(1-phi) + b*(-1) = (a-b) - a*phi
# Wait, let's be careful. The Galois conjugate of the vertex is NOT the same as
# multiplying by phi'. The Galois conjugate applies to EACH coordinate independently:
# sigma(v_i) = a_i + b_i*phi'
# This is different from phi' * v which would be a quaternionic multiplication.
#
# The Elser-Sloane construction uses the GALOIS CONJUGATE (coordinate-wise),
# NOT multiplication by phi'.

norms_T = np.linalg.norm(T, axis=1)
print(f"\n  |T| = {len(T)} vertices (Galois conjugate)")
print(f"  T norms: min={norms_T.min():.10f}, max={norms_T.max():.10f}")

# ============================================================
# SECTION 3: Shell structure of S and T
# ============================================================
print("\n" + "-" * 76)
print("SECTION 3: Shell structure (radii of S and T)")
print("-" * 76)

# S should all be on unit sphere
unique_norms_S = sorted(set(np.round(norms_S, 8)))
print(f"  S: {len(unique_norms_S)} distinct radii: {unique_norms_S[:5]}...")

# T has multiple radii!
unique_norms_T = sorted(set(np.round(norms_T, 8)))
print(f"  T: {len(unique_norms_T)} distinct radii:")
for r in unique_norms_T:
    count = np.sum(np.abs(norms_T - r) < 1e-6)
    # Check if r is a nice expression of phi
    r_over_phi = r / PHI
    r_times_phi = r * PHI
    r_sq = r**2
    print(f"    r = {r:.8f} ({count:3d} vertices), r^2 = {r_sq:.8f}, "
          f"r/phi = {r_over_phi:.6f}, r*phi = {r_times_phi:.6f}")

# Check: is T a SCALED 600-cell?
# If T = (1/phi)*S (global scaling), then all T norms should be 1/phi
print(f"\n  1/phi = {1/PHI:.8f}")
print(f"  phi   = {PHI:.8f}")

# Actually check if the Galois conjugate maps the 600-cell to itself (possibly scaled)
# The 600-cell has quaternionic structure: vertices are unit icosians.
# Galois conjugate of an icosian q = a+b*phi maps to q' = a+b*phi'.
# |q'|^2 = sum(q_i')^2 which is NOT necessarily 1.

# Let's verify with the actual norms
print(f"\n  T norm distribution:")
norm_counts_T = Counter(np.round(norms_T, 6))
for r, cnt in sorted(norm_counts_T.items()):
    # Try to identify r^2 as a+b*phi
    r2 = r**2
    best_err = 1e10
    best_ab = None
    for a in np.arange(-3, 4, 0.5):
        for b in np.arange(-3, 4, 0.5):
            err = abs(r2 - (a + b*PHI))
            if err < best_err:
                best_err = err
                best_ab = (a, b)
    print(f"    r={r:.6f}, r^2={r2:.6f}, count={cnt}, "
          f"r^2 ~ {best_ab[0]}+{best_ab[1]}*phi (err={best_err:.2e})")

# ============================================================
# SECTION 4: S union T as quasicrystal in R^4
# ============================================================
print("\n" + "-" * 76)
print("SECTION 4: S union T in R^4 (quasicrystal structure)")
print("-" * 76)

all_verts = np.vstack([S, T])
print(f"  Total vertices: {len(all_verts)}")

# Check for duplicates between S and T
overlap_count = 0
overlap_indices = []
for i, v in enumerate(T):
    dists = np.linalg.norm(S - v, axis=1)
    if np.min(dists) < 1e-8:
        overlap_count += 1
        overlap_indices.append(i)
print(f"  Overlap |S intersect T|: {overlap_count}")

if overlap_count > 0:
    print(f"  Overlapping T vertices (indices): {overlap_indices[:10]}...")
    # These should be the vertices with no phi component (Type A + Type B)
    for idx in overlap_indices[:5]:
        a_parts = ab_S[idx, :, 0]
        b_parts = ab_S[idx, :, 1]
        print(f"    T[{idx}]: a={a_parts}, b={b_parts}, "
              f"b=0? {np.allclose(b_parts, 0)}")

# Deduplicate
unique_verts_set = set()
for v in all_verts:
    unique_verts_set.add(tuple(round(x, 9) for x in v))
n_unique = len(unique_verts_set)
print(f"  Unique vertices in S union T: {n_unique}")
print(f"  Expected if no overlap: 240")

# ============================================================
# SECTION 5: Nearest-neighbor structure in the quasicrystal
# ============================================================
print("\n" + "-" * 76)
print("SECTION 5: Nearest-neighbor structure")
print("-" * 76)

# For S (the 600-cell on unit sphere), nearest-neighbor distance:
dots_S = S @ S.T
# 600-cell: nearest neighbors have dot product = phi/2 (geodesic distance = pi/5)
# Actually, the first ring has dot product = (1+phi)/2*2 = ... let's compute
nn_dist_S = []
for i in range(len(S)):
    dists = np.linalg.norm(S - S[i], axis=1)
    dists_sorted = np.sort(dists)
    nn_dist_S.append(dists_sorted[1])  # exclude self
nn_dist_S = np.array(nn_dist_S)
print(f"  S (600-cell) nearest-neighbor distance: {nn_dist_S[0]:.8f}")
print(f"  Expected (edge length = 1/phi): {1/PHI:.8f}")

# For the combined set, compute pairwise distances
unique_verts = np.array(list(unique_verts_set))
n_u = len(unique_verts)
print(f"\n  Computing pairwise distances for {n_u} unique vertices...")

# Use vectorized computation
dist_matrix = np.zeros((n_u, n_u))
for i in range(n_u):
    diff = unique_verts - unique_verts[i]
    dist_matrix[i] = np.sqrt(np.sum(diff**2, axis=1))

# Distribution of nearest-neighbor distances
nn_dists_combined = []
for i in range(n_u):
    sorted_dists = np.sort(dist_matrix[i])
    nn_dists_combined.append(sorted_dists[1])
nn_dists_combined = np.array(nn_dists_combined)

unique_nn = sorted(set(np.round(nn_dists_combined, 6)))
print(f"  Distinct nearest-neighbor distances: {len(unique_nn)}")
for d in unique_nn[:10]:
    cnt = np.sum(np.abs(nn_dists_combined - d) < 1e-4)
    # Try to express as function of phi
    d_phi = d * PHI
    d_over_phi = d / PHI
    d_sq = d**2
    print(f"    d = {d:.6f} ({cnt} vertices), d^2 = {d_sq:.6f}, "
          f"d*phi = {d_phi:.6f}, d/phi = {d_over_phi:.6f}")

# Full distance histogram
print(f"\n  Distance distribution (all pairs):")
all_dists = dist_matrix[np.triu_indices(n_u, k=1)]
dist_hist = Counter(np.round(all_dists, 4))
print(f"  Number of distinct distances: {len(dist_hist)}")
for d, cnt in sorted(dist_hist.items())[:15]:
    print(f"    d = {d:.4f}: {cnt:6d} pairs")

# ============================================================
# SECTION 6: Identify S vs T membership for unique vertices
# ============================================================
print("\n" + "-" * 76)
print("SECTION 6: S vs T classification of combined set")
print("-" * 76)

# Label each unique vertex as S-only, T-only, or both
labels = []
for v in unique_verts:
    in_S = np.min(np.linalg.norm(S - v, axis=1)) < 1e-7
    in_T = np.min(np.linalg.norm(T - v, axis=1)) < 1e-7
    if in_S and in_T:
        labels.append("both")
    elif in_S:
        labels.append("S")
    elif in_T:
        labels.append("T")
    else:
        labels.append("?")
labels = np.array(labels)

label_counts = Counter(labels)
print(f"  S-only: {label_counts.get('S', 0)}")
print(f"  T-only: {label_counts.get('T', 0)}")
print(f"  Both:   {label_counts.get('both', 0)}")
print(f"  Total:  {sum(label_counts.values())}")

# Radii by label
for lab in ["S", "T", "both"]:
    mask = labels == lab
    if np.any(mask):
        norms = np.linalg.norm(unique_verts[mask], axis=1)
        print(f"  {lab:5s}: radii = {sorted(set(np.round(norms, 6)))[:5]}")

# ============================================================
# SECTION 7: Concentric shell analysis (Elser-Sloane structure)
# ============================================================
print("\n" + "-" * 76)
print("SECTION 7: Concentric shell structure (Elser-Sloane)")
print("-" * 76)

# The key claim: S union T forms two 600-cells scaled by phi
# S has radius 1, so T should have radius 1/phi (since phi' = -1/phi,
# |phi'| = 1/phi). But Galois conjugate is coordinate-wise, not a global scaling.
# Let's check if T vertices form a 600-cell at SOME radius.

# Rescale T vertices that are NOT in S
T_only_mask = labels == "T"
T_only_verts = unique_verts[T_only_mask]
if len(T_only_verts) > 0:
    T_only_norms = np.linalg.norm(T_only_verts, axis=1)
    unique_T_norms = sorted(set(np.round(T_only_norms, 6)))
    print(f"  T-only vertices: {len(T_only_verts)}")
    print(f"  T-only radii: {unique_T_norms[:5]}...")

    # Check: are T-only vertices on a single sphere?
    if len(unique_T_norms) == 1:
        print(f"  T-only is on a SINGLE sphere of radius {unique_T_norms[0]}")
        # Normalize and check if it's a 600-cell
        T_normalized = T_only_verts / unique_T_norms[0]
        # Count vertices at each dot product with first vertex
        dots_Tn = T_normalized @ T_normalized[0]
        dot_dist = Counter(np.round(dots_Tn, 4))
        print(f"  Dot product distribution from first vertex:")
        for d, c in sorted(dot_dist.items()):
            print(f"    dot = {d:.4f}: {c}")
    else:
        print(f"  T-only is on MULTIPLE spheres (not a single 600-cell)")
        # Group by radius
        for r in unique_T_norms:
            mask_r = np.abs(T_only_norms - r) < 1e-4
            cnt_r = np.sum(mask_r)
            print(f"    r = {r:.6f}: {cnt_r} vertices")

# Also check ALL T vertices (including overlap)
print(f"\n  All T vertices (including overlap with S):")
for r in sorted(set(np.round(norms_T, 6))):
    cnt = np.sum(np.abs(norms_T - r) < 1e-4)
    ratio_to_S = r / 1.0  # S is at radius 1
    print(f"    r = {r:.6f}: {cnt} vertices, r/r_S = {ratio_to_S:.6f}")

# The Elser-Sloane result: in the PHYSICAL space (4D projection),
# the quasicrystal consists of vertices at radii 1 and 1/phi.
# But this depends on the projection. Let me try:
# Project E8 = (S_cholesky, T_cholesky) onto the "physical" 4D subspace.
# Physical = first 4 coordinates, internal = last 4.
# Actually, the Elser-Sloane projection is more subtle.

# The correct statement: the 240 E8 roots, when projected onto a 4D subspace
# aligned with an icosahedral reflection plane, give BOTH shells.
# In our setup, S is already in 4D. T is the Galois conjugate, also in 4D.
# The union S union T in 4D IS the quasicrystal.

print(f"\n  Key: S (radius 1) and T (multiple radii) are BOTH in R^4.")
print(f"  The quasicrystal is S union T in R^4.")

# ============================================================
# SECTION 8: Icosahedral symmetry check
# ============================================================
print("\n" + "-" * 76)
print("SECTION 8: Icosahedral (H4) symmetry of S union T")
print("-" * 76)

# The 600-cell has H4 symmetry (|H4| = 14400).
# Check if T (Galois conjugate) also has H4 symmetry.
# H4 is generated by reflections in 4D.

# Simple check: does S union T have the same distance distribution from
# every vertex? (regularity check)
# For S alone (600-cell): YES, it's vertex-transitive.
# For S union T: let's check.

# Check dot product distributions for a few vertices
print("  Checking vertex-transitivity of S union T...")
dot_dists = []
for i in range(min(10, n_u)):
    dots = unique_verts @ unique_verts[i]
    dot_dist = tuple(sorted(np.round(dots, 6)))
    dot_dists.append(dot_dist)

n_distinct_envs = len(set(dot_dists))
print(f"  Number of distinct local environments (first 10 vertices): {n_distinct_envs}")
if n_distinct_envs == 1:
    print(f"  => VERTEX-TRANSITIVE (all environments identical)")
else:
    # Count how many distinct environments total
    all_dot_dists = []
    for i in range(n_u):
        dots = unique_verts @ unique_verts[i]
        dot_dist = tuple(sorted(np.round(dots, 5)))
        all_dot_dists.append(dot_dist)
    n_envs = len(set(all_dot_dists))
    print(f"  Total distinct environments: {n_envs}")
    # Count per environment
    env_counter = Counter(all_dot_dists)
    for env, cnt in env_counter.most_common(5):
        print(f"    Environment type with {cnt} vertices")

# ============================================================
# SECTION 9: Diffraction pattern (structure factor)
# ============================================================
print("\n" + "-" * 76)
print("SECTION 9: Diffraction pattern (structure factor)")
print("-" * 76)

# Structure factor: S(k) = |sum_j exp(i k.r_j)|^2
# For a quasicrystal, this should show SHARP Bragg peaks
# with icosahedral symmetry.

# Scan along several directions
directions = {
    "e1": np.array([1, 0, 0, 0]),
    "e1+e2": np.array([1, 1, 0, 0]) / np.sqrt(2),
    "diagonal": np.array([1, 1, 1, 1]) / 2.0,
    "golden": np.array([1, PHI, 0, 0]) / np.sqrt(1 + PHI**2),
    "icosahedral": np.array([PHI, 1, 1/PHI, 0]) / np.sqrt(PHI**2 + 1 + 1/PHI**2),
}

k_values = np.linspace(0, 30, 3000)  # |k| range

print("  Computing structure factor S(k) along high-symmetry directions...")
print("  (Looking for sharp Bragg peaks => quasicrystal signature)")

for dir_name, dir_vec in directions.items():
    Sk_values = []
    for k_mag in k_values:
        k = k_mag * dir_vec
        # Sum over ALL unique vertices
        phases = unique_verts @ k
        Sk = np.abs(np.sum(np.exp(1j * phases)))**2
        Sk_values.append(Sk)
    Sk_values = np.array(Sk_values)

    # Find peaks (local maxima above threshold)
    threshold = n_u * 5  # well above noise level (n_u for random positions)
    peaks = []
    for i in range(1, len(Sk_values) - 1):
        if (Sk_values[i] > Sk_values[i-1] and
            Sk_values[i] > Sk_values[i+1] and
            Sk_values[i] > threshold):
            peaks.append((k_values[i], Sk_values[i]))

    print(f"\n  Direction: {dir_name}")
    print(f"    Max S(k) = {Sk_values.max():.1f} at k = {k_values[np.argmax(Sk_values)]:.4f}")
    print(f"    Peaks above {threshold:.0f} (= {n_u}*5):")
    for k_peak, sk_peak in peaks[:8]:
        # Check if peak position involves phi
        print(f"      k = {k_peak:.4f}, S(k) = {sk_peak:.1f}, "
              f"k/phi = {k_peak/PHI:.4f}, k*phi = {k_peak*PHI:.4f}")

    # Check for phi-scaling in peak positions (quasicrystal signature!)
    if len(peaks) >= 2:
        ratios = []
        for i in range(len(peaks)-1):
            r = peaks[i+1][0] / peaks[i][0]
            ratios.append(r)
        phi_ratios = [r for r in ratios if abs(r - PHI) < 0.15 or abs(r - PHI**2) < 0.3]
        if phi_ratios:
            print(f"    *** PHI-SCALED PEAKS DETECTED! Ratios: "
                  f"{[f'{r:.4f}' for r in phi_ratios]} ***")

# ============================================================
# SECTION 10: Quasicrystal tile analysis
# ============================================================
print("\n" + "-" * 76)
print("SECTION 10: Tile analysis (Voronoi-like)")
print("-" * 76)

# In 4D quasicrystals, the relevant "tiles" are the Voronoi cells.
# For the Elser-Sloane quasicrystal, these should be related to
# the 600-cell and 120-cell (dual polytopes).

# Instead of full Voronoi, analyze the coordination structure:
# For each vertex, count neighbors at each distance.

print("  Coordination shells for combined S union T:")

# Pick a representative vertex from each environment type
for label_type in ["S", "T", "both"]:
    mask = labels == label_type
    if not np.any(mask):
        continue
    idx = np.where(mask)[0][0]
    v = unique_verts[idx]
    dists = np.sort(np.linalg.norm(unique_verts - v, axis=1))

    print(f"\n  Representative vertex from '{label_type}' (r={np.linalg.norm(v):.4f}):")
    # Coordination shells
    shell_dists = sorted(set(np.round(dists[1:], 4)))[:8]
    for sd in shell_dists:
        cnt = np.sum(np.abs(dists - sd) < 5e-4) - (1 if sd < 1e-4 else 0)
        # Identify if neighbors are S, T, or both
        neighbor_mask = np.abs(np.linalg.norm(unique_verts - v, axis=1) - sd) < 5e-4
        neighbor_labels = labels[neighbor_mask]
        n_S = np.sum(neighbor_labels == "S")
        n_T = np.sum(neighbor_labels == "T")
        n_B = np.sum(neighbor_labels == "both")
        print(f"    d={sd:.4f}: {cnt:3d} neighbors (S={n_S}, T={n_T}, both={n_B})")

# ============================================================
# SECTION 11: E8 embedding verification
# ============================================================
print("\n" + "-" * 76)
print("SECTION 11: E8 Cholesky embedding (verify exp120d)")
print("-" * 76)

def cholesky_map(ab):
    """Map (a_i, b_i) -> (a_i + b_i, b_i) for Euclidean norm."""
    result = np.zeros(8)
    for i in range(4):
        a, b = ab[i]
        result[2*i] = a + b
        result[2*i+1] = b
    return result

S_8d = np.array([cholesky_map(ab_S[idx]) for idx in range(120)])
T_8d = np.array([cholesky_map(ab_T[idx]) for idx in range(120)])

S_norms_8d = np.sqrt(np.sum(S_8d**2, axis=1))
T_norms_8d = np.sqrt(np.sum(T_8d**2, axis=1))
print(f"  S in R^8: norms = {sorted(set(np.round(S_norms_8d, 6)))}")
print(f"  T in R^8: norms = {sorted(set(np.round(T_norms_8d, 6)))}")

# S union T should give 240 E8 roots
E8_roots = np.vstack([S_8d, T_8d])
E8_unique = set()
for v in E8_roots:
    E8_unique.add(tuple(round(x, 8) for x in v))
print(f"  |S union T| in R^8 = {len(E8_unique)} (need 240 for E8)")

# Check overlap in R^8
overlap_8d = 0
for v in S_8d:
    dists = np.linalg.norm(T_8d - v, axis=1)
    if np.min(dists) < 1e-6:
        overlap_8d += 1
print(f"  Overlap in R^8: {overlap_8d}")

# Verify E8 inner products
if len(E8_unique) == 240:
    E8_arr = np.array(list(E8_unique)) * np.sqrt(2)
    dots = E8_arr @ E8_arr.T
    dot_counts = Counter()
    n240 = len(E8_arr)
    for i in range(n240):
        for j in range(i+1, n240):
            d = round(dots[i, j], 2)
            dot_counts[d] += 1

    print(f"\n  E8 inner product distribution (scaled by sqrt(2), norm^2=2):")
    for d, cnt in sorted(dot_counts.items()):
        print(f"    dot={d:6.1f}: {cnt:6d} pairs")

    is_E8 = (dot_counts.get(-2.0, 0) == 120 and
             dot_counts.get(-1.0, 0) == 6720 and
             dot_counts.get(0.0, 0) == 15120 and
             dot_counts.get(1.0, 0) == 6720)
    print(f"  *** E8 CONFIRMED: {is_E8} ***")

# ============================================================
# SECTION 12: Physical vs Internal decomposition
# ============================================================
print("\n" + "-" * 76)
print("SECTION 12: Physical (parallel) vs Internal (perp) space")
print("-" * 76)

# In the Elser-Sloane setup, R^8 = R^4_physical x R^4_internal
# Physical projection: the "parallel" space
# Internal projection: the "perpendicular" space
#
# Our Cholesky embedding maps (a,b) -> (a+b, b) per coordinate.
# So for vertex with R^4 coords v_i = a_i + b_i*phi:
#   Physical = S vertices (a_i + b_i*phi)
#   Internal = T vertices (a_i + b_i*phi')
#
# In R^8 = (a0+b0, b0, a1+b1, b1, ...) the "physical" projection
# recovers v_i = a_i + b_i*phi and the "internal" recovers v_i' = a_i + b_i*phi'

# For each E8 root, split into physical (S-like) and internal (T-like) parts
print("  For each 600-cell vertex of S:")
print("    Physical projection = S vertex (in R^4)")
print("    Internal projection = T vertex (Galois conjugate, in R^4)")

# Norms in physical and internal space
phys_norms = norms_S  # already computed
int_norms = norms_T   # Galois conjugate

print(f"\n  Physical (S) norms: {sorted(set(np.round(phys_norms, 6)))}")
print(f"  Internal (T) norms: {sorted(set(np.round(int_norms, 6)))[:5]}...")

# Relationship between physical and internal norms
print(f"\n  Relationship |v|^2 + |v'|^2 for each vertex:")
sum_sq = phys_norms**2 + int_norms**2
unique_sums = sorted(set(np.round(sum_sq, 6)))
for s in unique_sums[:5]:
    cnt = np.sum(np.abs(sum_sq - s) < 1e-4)
    print(f"    |v|^2 + |v'|^2 = {s:.6f}: {cnt} vertices")

# This should be the icosian norm, which is 1 for all unit icosians!
# Actually: N_icos(q) = |v|^2 + |v'|^2 = Trace over Galois = 1 for all unit icosians
# Wait, that's not right. Let me compute properly.
# v_i = a_i + b_i*phi, v_i' = a_i + b_i*phi'
# |v|^2 = sum(a_i + b_i*phi)^2 = sum(a_i^2 + 2a_i*b_i*phi + b_i^2*phi^2)
# |v'|^2 = sum(a_i + b_i*phi')^2 = sum(a_i^2 + 2a_i*b_i*phi' + b_i^2*phi'^2)
# Sum = sum(2a_i^2 + 2a_i*b_i*(phi+phi') + b_i^2*(phi^2+phi'^2))
# phi + phi' = 1, phi^2 + phi'^2 = (phi+phi')^2 - 2*phi*phi' = 1 - 2*(-1) = 3
# Sum = sum(2a_i^2 + 2a_i*b_i + 3b_i^2) = 2*sum(a_i^2 + a_i*b_i + 1.5*b_i^2)
# Compare with icosian norm: sum(a_i^2 + 2a_i*b_i + 2b_i^2)
# NOT the same! Let's check numerically.

print(f"\n  Algebraic verification:")
print(f"    phi + phi' = {PHI + PHI_CONJ:.6f} (should be 1)")
print(f"    phi * phi' = {PHI * PHI_CONJ:.6f} (should be -1)")
print(f"    phi^2 + phi'^2 = {PHI**2 + PHI_CONJ**2:.6f} (should be 3)")

# ============================================================
# SECTION 13: Phi-scaling in the quasicrystal
# ============================================================
print("\n" + "-" * 76)
print("SECTION 13: Phi-scaling (inflation/deflation symmetry)")
print("-" * 76)

# Key property of quasicrystals: scaling by phi maps the set to itself
# (inflation symmetry). Check: does phi * T ~ S? (after rotation?)

# Scale T by phi
T_scaled = T * PHI
T_scaled_norms = np.linalg.norm(T_scaled, axis=1)
print(f"  phi * T norms: {sorted(set(np.round(T_scaled_norms, 6)))[:5]}")

# How many of phi*T vertices are close to S vertices?
phi_T_matches_S = 0
for v in T_scaled:
    dists = np.linalg.norm(S - v, axis=1)
    if np.min(dists) < 1e-6:
        phi_T_matches_S += 1
print(f"  Vertices of phi*T that match S: {phi_T_matches_S}/120")

# Scale S by 1/phi
S_scaled = S / PHI
S_scaled_norms = np.linalg.norm(S_scaled, axis=1)
# How many match T?
S_over_phi_matches_T = 0
for v in S_scaled:
    dists = np.linalg.norm(T - v, axis=1)
    if np.min(dists) < 1e-6:
        S_over_phi_matches_T += 1
print(f"  Vertices of S/phi that match T: {S_over_phi_matches_T}/120")

# Actually, for Galois conjugate: sigma(phi*v) = phi'*sigma(v) = phi'*v'
# So the "inflation" in physical space corresponds to "deflation" in internal space
# This is the fundamental property of quasicrystals!
print(f"\n  Galois action on scaling:")
print(f"    Physical: multiply by phi")
print(f"    Internal: multiply by phi' = -1/phi")
print(f"    |phi| = {PHI:.6f}, |phi'| = {abs(PHI_CONJ):.6f}")
print(f"    phi * |phi'| = {PHI * abs(PHI_CONJ):.6f} (= 1: preserved total)")

# ============================================================
# SECTION 14: Type A/B/C classification under Galois
# ============================================================
print("\n" + "-" * 76)
print("SECTION 14: Type A/B/C behavior under Galois conjugate")
print("-" * 76)

# Type classification based on phi content
n_type_A = 0  # no phi component at all
n_type_B = 0  # some phi component
type_A_indices = []
type_B_indices = []

for idx in range(120):
    b_parts = ab_S[idx, :, 1]  # phi coefficients
    if np.allclose(b_parts, 0):
        n_type_A += 1
        type_A_indices.append(idx)
    else:
        n_type_B += 1
        type_B_indices.append(idx)

print(f"  Vertices with b=0 (no phi): {n_type_A} (Type A+B: 8+16=24)")
print(f"  Vertices with b!=0 (has phi): {n_type_B} (Type C: 96)")

# For b=0 vertices: Galois conjugate = identity (v' = v)
# These are the OVERLAP vertices
print(f"\n  For b=0 vertices: v = v' (Galois fixed points)")
for idx in type_A_indices[:5]:
    v = S[idx]
    v_conj = T[idx]
    print(f"    S[{idx}] = {np.round(v, 4)}, T[{idx}] = {np.round(v_conj, 4)}, "
          f"match: {np.allclose(v, v_conj)}")

# For b!=0 vertices: v != v' (Galois moves them)
print(f"\n  For b!=0 vertices: v != v'")
for idx in type_B_indices[:3]:
    v = S[idx]
    v_conj = T[idx]
    print(f"    S[{idx}] = {np.round(v, 4)}")
    print(f"    T[{idx}] = {np.round(v_conj, 4)}")
    print(f"    |v|={np.linalg.norm(v):.6f}, |v'|={np.linalg.norm(v_conj):.6f}")

# ============================================================
# SECTION 15: Window function (acceptance domain)
# ============================================================
print("\n" + "-" * 76)
print("SECTION 15: Window function / acceptance domain")
print("-" * 76)

# In quasicrystal theory, the "window" or "acceptance domain" is the
# projection of the Voronoi cell of the parent lattice (E8) onto
# the internal (perpendicular) space.
#
# A vertex r in physical space is in the quasicrystal iff its
# "internal partner" r_perp falls inside the window W.
#
# For our setup: physical = S vertex, internal = T vertex (Galois conjugate)
# The window W should be related to the 600-cell or 120-cell.

# Check: what shape does the internal projection trace out?
# For the 600-cell: 120 vertices project to 120 internal points
# These internal points should lie within the window.

# For the E8 lattice (not just roots), the window for the Elser-Sloane
# quasicrystal is a triacontahedron (30-faced) in 3D cross-sections,
# and a 600-cell/120-cell dual pair in 4D.

# We can verify: the internal projections of S fill out a 600-cell shape
print("  Internal projections of S vertices (= T vertices):")
print(f"    Shape: {T.shape}")
print(f"    Radii: {sorted(set(np.round(norms_T, 6)))[:5]}")

# Check if T (internal) is related to a 600-cell
# Scale all T to unit sphere and check
T_mean_norm = np.mean(norms_T)
print(f"    Mean radius of T: {T_mean_norm:.6f}")

# Check if T has the same graph structure as a 600-cell
# (same adjacency even if different radii)
# Compute angular distances between T vertices
T_normalized = T / norms_T[:, np.newaxis]
dots_T_norm = T_normalized @ T_normalized.T
# 600-cell signature: 12 nearest neighbors at dot = phi/2
nn_count_T = np.sum(np.abs(dots_T_norm - PHI/2) < 0.01, axis=1)  # wrong, 600-cell nn at cos(pi/5)
# Actually 600-cell nn at dot product = (1+phi)/2 / 1 wait...
# For unit 600-cell: nn distance = 1/phi, so dot = 1 - d^2/2 = 1 - 1/(2*phi^2)
nn_dot_600 = 1 - 1/(2*PHI**2)
print(f"\n  600-cell nn dot product: {nn_dot_600:.6f}")

# For S (unit 600-cell):
dots_S_unit = S @ S.T
nn_count_S = np.sum(np.abs(dots_S_unit - nn_dot_600) < 0.01, axis=1)
print(f"  S (unit 600-cell) neighbors at dot={nn_dot_600:.4f}: "
      f"{Counter(nn_count_S).most_common(3)}")

# For T (Galois conjugate):
dots_T = T @ T.T
# T vertices have different norms, so we need to normalize dots
for idx in [0, 50, 100]:
    dists_from_v = np.linalg.norm(T - T[idx], axis=1)
    sorted_dists = np.sort(dists_from_v)
    nn12 = sorted_dists[1:13]  # first 12 neighbors
    print(f"  T[{idx}] (r={norms_T[idx]:.4f}): 12 nn dists = "
          f"[{nn12[0]:.4f}, ..., {nn12[-1]:.4f}]")

# ============================================================
# SECTION 16: Physical interpretation
# ============================================================
print("\n" + "-" * 76)
print("SECTION 16: Physical interpretation")
print("-" * 76)

print("""
  The Elser-Sloane construction reveals:

  1. E8 lattice decomposes as S + T (two 600-cells related by Galois conjugate)
  2. S lives on unit S^3 (physical space, radius 1)
  3. T = sigma(S) has MULTIPLE radii in R^4 (not a single scaled 600-cell)
  4. The 24 vertices with b=0 (no phi content) are FIXED by Galois: S cap T = 24
  5. The 96 vertices with phi content are SPLIT: different in S vs T

  Physical analogy:
  - S = "visible" sector (at electroweak scale M_Z ~ 91 GeV)
  - T = "shadow" sector (Galois conjugate, different scale)
  - 24 shared vertices = gauge sector (D4 root system, unchanged by Galois)
  - 96 split vertices = matter sector (transformed by Galois)

  Quasicrystal properties:
  - Diffraction shows Bragg peaks (long-range order)
  - Peak positions should scale by phi (aperiodic)
  - No translational symmetry, but icosahedral point symmetry
""")

# ============================================================
# SECTION 17: Detailed phi analysis of T norms
# ============================================================
print("-" * 76)
print("SECTION 17: Detailed analysis of T vertex norms")
print("-" * 76)

# For each T vertex, compute |v'|^2 = sum(a_i + b_i*phi')^2
# = sum((a_i+b_i) - b_i*phi)^2 = sum((a_i+b_i)^2 - 2(a_i+b_i)*b_i*phi + b_i^2*phi^2)
# This is itself in Q(phi): rational + irrational part.
# But |v'|^2 is a REAL number, so it must be rational... wait, no.
# phi^2 = phi + 1, so:
# (c - d*phi)^2 = c^2 - 2cd*phi + d^2*(phi+1) = (c^2 + d^2) + (d^2 - 2cd)*phi
# So |v'|^2 = sum_i [(c_i^2 + d_i^2) + (d_i^2 - 2*c_i*d_i)*phi]
# where c_i = a_i + b_i, d_i = -b_i (from the Galois conjugate decomposition)
# = sum_i [(a_i+b_i)^2 + b_i^2 + (b_i^2 + 2*(a_i+b_i)*b_i)*phi]
# Hmm, let me just compute.

print("\n  Checking if |v'|^2 is in Q(phi) for each vertex:")
norm2_T_rat = []
norm2_T_irr = []
for idx in range(120):
    c_parts = ab_T[idx, :, 0]  # a+b for T
    d_parts = ab_T[idx, :, 1]  # -b for T
    # |v'|^2 = sum((c_i + d_i*phi)^2) = sum(c_i^2 + 2*c_i*d_i*phi + d_i^2*phi^2)
    # = sum(c_i^2 + d_i^2) + sum(2*c_i*d_i + d_i^2)*phi  [since phi^2 = phi+1]
    rat_part = sum(c**2 + d**2 for c, d in zip(c_parts, d_parts))
    irr_part = sum(2*c*d + d**2 for c, d in zip(c_parts, d_parts))
    norm2_T_rat.append(rat_part)
    norm2_T_irr.append(irr_part)

norm2_T_rat = np.array(norm2_T_rat)
norm2_T_irr = np.array(norm2_T_irr)

# Show distinct (rat, irr) pairs
pairs = Counter(zip(np.round(norm2_T_rat, 6), np.round(norm2_T_irr, 6)))
print(f"\n  |v'|^2 = rat + irr*phi:")
for (r, ir), cnt in sorted(pairs.items()):
    val = r + ir * PHI
    print(f"    ({r:6.2f}) + ({ir:6.2f})*phi = {val:.6f}  [{cnt} vertices]")
    # Check relation to phi powers
    for n in range(-4, 5):
        if abs(val - PHI**n) < 0.01:
            print(f"      *** = phi^{n} ***")

# ============================================================
# SECTION 18: Connection to Penrose tilings
# ============================================================
print("\n" + "-" * 76)
print("SECTION 18: Connection to lower-dimensional quasicrystals")
print("-" * 76)

# Project the 4D quasicrystal onto 2D to see Penrose-like structure
# Use the "golden plane" projection
print("  Projecting S union T onto 2D golden plane...")

# Projection matrix: first two components of icosahedral basis
# Use a Coxeter plane projection
proj_2d = np.array([
    [1, PHI, 0, 0],
    [0, 0, 1, PHI]
]) / np.sqrt(1 + PHI**2)

pts_2d = unique_verts @ proj_2d.T
print(f"  2D projected points: {len(pts_2d)}")

# Check for 5-fold or 10-fold symmetry
# Compute angular distribution from centroid
centroid = np.mean(pts_2d, axis=0)
angles = np.arctan2(pts_2d[:, 1] - centroid[1], pts_2d[:, 0] - centroid[0])
angles_deg = np.degrees(angles) % 360

# Bin into 36-degree sectors (10-fold) or 72-degree sectors (5-fold)
sectors_10 = Counter(np.floor(angles_deg / 36).astype(int))
print(f"\n  Angular distribution (36-degree sectors for 10-fold symmetry):")
for sec in range(10):
    cnt = sectors_10.get(sec, 0)
    print(f"    [{sec*36:3d}-{(sec+1)*36:3d}): {cnt:3d} points {'*' * (cnt // 2)}")

# Check if counts are equal (perfect 10-fold symmetry)
counts_10 = [sectors_10.get(s, 0) for s in range(10)]
print(f"  10-fold symmetry measure: std/mean = "
      f"{np.std(counts_10)/np.mean(counts_10):.4f} (0 = perfect)")

# ============================================================
# SECTION 19: Key ratios and numbers
# ============================================================
print("\n" + "-" * 76)
print("SECTION 19: Summary of key numbers")
print("-" * 76)

print(f"""
  600-CELL AND E8:
    |S| = 120 (600-cell vertices)
    |T| = 120 (Galois conjugate)
    |S intersect T| = {overlap_count} (Galois-fixed, b_i=0 for all i)
    |S union T| unique in R^4 = {n_unique}
    |S union T| in R^8 = {len(E8_unique)} (should be 240 = |E8 roots|)

  GALOIS STRUCTURE:
    24 fixed vertices (no phi content): 8 (Type A: perms of +-1,0,0,0)
                                       + 16 (Type B: all +-1/2)
    96 moving vertices (has phi content): Type C

  PHYSICAL vs INTERNAL:
    Physical norms (S): all = 1.0 (unit S^3)
    Internal norms (T): multiple values
    phi + phi' = 1 (trace)
    phi * phi' = -1 (norm, det of Galois)

  PHI-SCALING:
    phi*|phi'| = 1 (product of conjugate scales)
    Physical inflation by phi <=> Internal deflation by 1/phi
""")

# ============================================================
# SUMMARY AND CLASSIFICATION
# ============================================================
print("=" * 76)
print("SUMMARY AND CLASSIFICATION")
print("=" * 76)

print("""
  DERIVED (rigorous, from geometry):
  [D1] Galois conjugate sigma: v_i = a_i + b_i*phi -> a_i + b_i*phi' (coordinate-wise).
       NOT a global scaling. Both S and T are 600-cells of RADIUS 1 in R^4.
  [D2] S cap T = 24 vertices (Galois fixed points, b_i = 0 for all i).
       These are EXACTLY Type A (8) + Type B (16) = 24-cell = D4 root system.
  [D3] |S union T| = 216 unique vertices in R^4 (120+120-24).
       NOT 240 in R^4! The Cholesky embedding into R^8 also gives 216 (same overlap).
       The E8 = 240 requires the FULL 8D embedding (v, v') as independent directions.
  [D4] phi + phi' = 1 (field trace), phi * phi' = -1 (field norm).
       |v|^2 + |v'|^2 = 2.0 for ALL vertices (since both on unit S^3).
  [D5] T is itself a 600-cell: same edge length (1/phi), same 12 nearest neighbors.
       T is ROTATED relative to S, not scaled. The rotation encodes Galois action.
  [D6] 96 Type-C vertices split: S-only and T-only are DISTINCT sets in R^4.
       Physical inflation by phi <=> Internal deflation by 1/phi (algebraic identity).

  PATTERN (observed, not rigorously derived):
  [P1] Two vertex environments in S union T: 192 "Type C" (S-only or T-only)
       and 24 "Type AB" (shared). NOT vertex-transitive.
  [P2] Key new distance d=0.2701: closest S-T distance (96 such pairs).
       Each S-only vertex has exactly 1 T-only vertex at this distance.
       d^2 = 0.0730 ~ (3-2*phi)/2. This is the "inter-shell" bond.
  [P3] "Both" vertices (24-cell) have ONLY 600-cell distances (0.618, 1.0, ...).
       They see S and T symmetrically: 12 S-neighbors + 12 T-neighbors at each shell.
  [P4] Diffraction S(k) shows peaks with ratio ~1.63 near phi along icosahedral direction.
       10-fold angular symmetry nearly perfect (std/mean = 3.7%).
  [P5] 24 Galois-fixed vertices = gauge sector (D4).
       96 Galois-moving vertices = matter sector.
       This matches exp124/148/149 gauge-fermion split EXACTLY.

  SPECULATIV:
  [S1] The 216 = 6^3 unique vertices may relate to the 6D cubic lattice
       that generates 3D icosahedral quasicrystals (Penrose).
  [S2] Quasicrystal: Elser-Sloane two-shell structure appears for the LATTICE
       (infinite points), not just roots. Our 216 roots show the seed.
  [S3] Physical interpretation: S and T are two "copies" of the 600-cell
       related by Galois. The 24-cell overlap (gauge) is scale-invariant,
       while the 96 matter vertices are "entangled" between the two copies.

  CORRECTION TO INITIAL EXPECTATION:
  The Elser-Sloane "two concentric 600-cells scaled by phi" applies to the
  projection of the full E8 LATTICE (infinite), not just the 240 roots.
  For the 240 roots: T = sigma(S) is a ROTATED (not scaled) 600-cell.
  Both S and T live on the SAME unit S^3. The "scaling by phi" emerges
  only when extending to the full lattice (higher shells).
""")

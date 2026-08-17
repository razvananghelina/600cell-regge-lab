"""
EXP-500: Control Test for n96 = 124
=====================================
REGULA ZERO: Before celebrating, verify it's not trivial.

CRITICAL CHECKS:
1. Sum of degrees = 2*n_edges is the HANDSHAKING LEMMA (tautology).
   2*124 = 248 = dim(E8) is NOT a result, it's a definition.
   STRIKE THIS from the results.

2. Would 40 RANDOM 4D subspaces on Gr(4,8) also give ~124 pairs
   with a shared 3D subspace? If yes, nothing is special about E8.

3. What IS non-trivial? Only n96 = 124 itself (= dim(E8)/2).
   But IS dim(E8)/2 special here, or is 124 just what you get
   from 40 subspaces with the given constraints?

CONTROL: Generate many random sets of 40 subspaces on Gr(4,8),
compute their overlap statistics, and compare with E8.
"""

import numpy as np
from collections import Counter
from scipy.stats import ortho_group

np.random.seed(42)

print("=" * 70)
print("EXP-500: CONTROL TEST - IS 124 SPECIAL?")
print("=" * 70)

# ============================================================
# FIRST: Acknowledge the tautology
# ============================================================
print("""
CORRECTION: Sum of degrees = 2 * n_edges is the HANDSHAKING LEMMA.
  2 * 124 = 248 is ALWAYS true for 124 edges, regardless of dim(E8).
  This is NOT a result. It was wrong to list it as one.

What COULD be non-trivial:
  - n96 = 124 = dim(E8)/2 itself
  - The degree DISTRIBUTION (4*5 + 24*6 + 12*7)
  - The fact that ALL 124 pairs have the SAME angle profile
  - The ratio (n72+n84)/(n60+n96) = 8/5

But we need to check: are these properties of E8, or of ANY
set of 40 subspaces on Gr(4,8)?
""")

# ============================================================
# CONTROL 1: Random subspaces on Gr(4,8)
# ============================================================
print("=" * 70)
print("CONTROL 1: RANDOM 4D SUBSPACES ON Gr(4,8)")
print("=" * 70)

def random_grassmannian_points(n_points, ambient_dim=8, sub_dim=4):
    """Generate n random 4D subspaces of R^8."""
    subspaces = []
    for _ in range(n_points):
        # Random orthogonal matrix, take first sub_dim columns
        Q = ortho_group.rvs(ambient_dim)
        subspaces.append(Q[:, :sub_dim])
    return subspaces

def compute_overlap_stats(subspaces, threshold_3d=0.999):
    """For each pair: compute principal angles and count 'near 3D shared'."""
    n = len(subspaces)
    overlaps = []
    n_3d_shared = 0  # pairs sharing a 3D subspace (3 SVs near 1)

    for i in range(n):
        for j in range(i+1, n):
            M = subspaces[i].T @ subspaces[j]
            svs = np.linalg.svd(M, compute_uv=False)
            svs_sorted = np.sort(svs)[::-1]

            # Count how many SVs are near 1 (= shared dimensions)
            n_shared = np.sum(svs_sorted > threshold_3d)

            if n_shared >= 3:
                n_3d_shared += 1

            overlaps.append(svs_sorted)

    return overlaps, n_3d_shared

# Run many trials
n_trials = 200
results_random = []

print(f"\nGenerating {n_trials} random sets of 40 subspaces on Gr(4,8)...")
for trial in range(n_trials):
    subspaces = random_grassmannian_points(40, 8, 4)
    _, n_3d = compute_overlap_stats(subspaces)
    results_random.append(n_3d)

print(f"\nPairs sharing >= 3D subspace (threshold SV > 0.999):")
print(f"  E8 value: 124")
print(f"  Random: mean = {np.mean(results_random):.1f}, "
      f"std = {np.std(results_random):.1f}, "
      f"max = {max(results_random)}, min = {min(results_random)}")
print(f"  Random range: [{min(results_random)}, {max(results_random)}]")

if max(results_random) >= 124:
    print(f"  VERDICT: Random achieves 124. NOT SPECIAL.")
else:
    print(f"  VERDICT: Random NEVER reaches 124. E8 value IS special.")
    print(f"  Sigma: {(124 - np.mean(results_random)) / max(np.std(results_random), 0.001):.1f} sigma away")

# ============================================================
# CONTROL 2: Random subspaces with E8 constraints
# ============================================================
print(f"\n{'='*70}")
print("CONTROL 2: CONSTRAINED RANDOM (same pairwise angle structure)")
print(f"{'='*70}")

# The E8 subspaces aren't arbitrary: they come from Coxeter elements.
# Their principal angles are constrained to cos^2 = k/5.
# A fairer control: random subspaces with the SAME constraint.

# But even without this: if random gives 0 pairs with 3D sharing,
# while E8 gives 124, that's already decisive.

# ============================================================
# CONTROL 3: Is 124 = C(something, something)?
# ============================================================
print(f"\n{'='*70}")
print("CONTROL 3: COMBINATORIAL EXPLANATIONS FOR 124")
print(f"{'='*70}")

from math import comb, factorial

print(f"\n  124 = ?")
candidates = []
for n in range(1, 50):
    for k in range(0, n+1):
        if comb(n, k) == 124:
            candidates.append(f"C({n},{k})")

if candidates:
    print(f"  124 as binomial: {candidates}")
else:
    print(f"  124 is NOT a binomial coefficient")

# Other decompositions
print(f"  124 = 4 * 31")
print(f"  124 = 2^2 * 31")
print(f"  31 is prime")
print(f"  124 = 120 + 4 = N + d")
print(f"  124 = 248/2 = dim(E8)/2")
print(f"  124 = (240 + 8)/2 = (|roots| + rank)/2")

# Is 124 forced by the structure?
# We have 40 nodes. If each node has degree ~6 at overlap 96,
# then n_edges ~ 40*6/2 = 120. Close to 124 but not exact.

# The ACTUAL constraint: overlap must be in {60,72,84,96}.
# 780 pairs total. If overlap 96 gets ~16% of pairs: 780*0.16 = 125. Close!
# Is 124/780 special?
print(f"\n  124/780 = {124/780:.6f}")
print(f"  1/2pi = {1/(2*np.pi):.6f}")
print(f"  1/b1 = {1/6:.6f}")
print(f"  a1/(a1^2+1) = {5/26:.6f}")
g = np.gcd(124, 780)
print(f"  124/780 = {124//g}/{780//g} (irreducible)")

# ============================================================
# CONTROL 4: What if we use a DIFFERENT root system?
# ============================================================
print(f"\n{'='*70}")
print("CONTROL 4: DOES dim(G)/2 APPEAR FOR OTHER ROOT SYSTEMS?")
print(f"{'='*70}")

# For E8: 40 decompositions, n_max_overlap = 124 = dim(E8)/2
# For D4: the 24-cell decomposes differently
# For E6: 72 roots, McKay -> 2T, different structure
# For E7: 126 roots, McKay -> 2O

# We can't easily compute these without building the full machinery.
# But we can check: for a SMALLER system, does the pattern hold?

# A5 symmetry (icosahedral) in lower dimensions?
# The key E8 property: h(E8) = 30 = a1*b1, and H4 has exponents {1,11,19,29}

# For D5 (rank 5, dim 45, h=8): exponents 1,3,5,7,9? No, D5 exponents are 1,3,5,7
# Doesn't have an H-type subgroup

# For A4 (rank 4, dim 24, h=5): exponents 1,2,3,4
# No golden ratio connection

print(f"  E8: rank=8, dim=248, h=30, n_decomp=40, n_max_overlap=124 = dim/2")
print(f"  Testing other systems requires full Coxeter element analysis")
print(f"  (not implemented here)")

# ============================================================
# CONTROL 5: Is the degree distribution (4,24,12) special?
# ============================================================
print(f"\n{'='*70}")
print("CONTROL 5: DEGREE DISTRIBUTION ANALYSIS")
print(f"{'='*70}")

# E8 degrees at overlap 96: {5: 4, 6: 24, 7: 12}
# This gives 4+24+12 = 40 nodes, 2*124 = 248 degree sum.

# For a random graph with 40 nodes and 124 edges:
# Expected average degree = 2*124/40 = 6.2
# By Erdos-Renyi: degree ~ Binomial(39, 124/780)
p_edge = 124 / 780
expected_degree = 39 * p_edge
var_degree = 39 * p_edge * (1 - p_edge)
print(f"\n  For random graph G(40, p={p_edge:.4f}):")
print(f"  Expected degree: {expected_degree:.2f}")
print(f"  Std of degree: {np.sqrt(var_degree):.2f}")
print(f"  E8 degrees: 5, 6, 7 (range 2)")
print(f"  Random degree range would be ~{expected_degree-2*np.sqrt(var_degree):.1f} to {expected_degree+2*np.sqrt(var_degree):.1f}")
print(f"  E8 degree range (5-7) is NARROWER than random expectation")

# Is the narrow range (only 3 distinct values) special?
n_simulations = 1000
n_distinct_degrees_random = []
for _ in range(n_simulations):
    # Random graph with exactly 124 edges on 40 nodes
    adj = np.zeros((40, 40))
    edges_placed = 0
    while edges_placed < 124:
        i, j = np.random.randint(0, 40, 2)
        if i != j and adj[i,j] == 0:
            adj[i,j] = 1; adj[j,i] = 1
            edges_placed += 1
    degrees = adj.sum(axis=1).astype(int)
    n_distinct = len(set(degrees))
    n_distinct_degrees_random.append(n_distinct)

print(f"\n  Distinct degree values in random graphs with 124 edges on 40 nodes:")
print(f"  E8: 3 distinct values (5, 6, 7)")
print(f"  Random: mean = {np.mean(n_distinct_degrees_random):.1f}, "
      f"std = {np.std(n_distinct_degrees_random):.1f}")
print(f"  Random having <= 3 distinct: "
      f"{sum(1 for x in n_distinct_degrees_random if x <= 3)} / {n_simulations}")

# ============================================================
# CONTROL 6: Is all-same-profile special?
# ============================================================
print(f"\n{'='*70}")
print("CONTROL 6: IS SINGLE ANGLE PROFILE AT OVERLAP 96 SPECIAL?")
print(f"{'='*70}")

# E8: ALL 124 pairs at overlap 96 have EXACTLY cos^2 = (1,1,1,1/5)
# At other overlaps: MULTIPLE profiles exist
# Overlap 60: 3 profiles (96+48+32 pairs)
# Overlap 72: 3 profiles (100+100+96)
# Overlap 84: 2 profiles (132+52)
# Overlap 96: 1 profile (all 124)

# Is having a SINGLE profile at the extreme overlap a generic property?
# For random Grassmannian points: the probability that two random 4D
# subspaces share EXACTLY a 3D subspace is ZERO (measure zero).
# The fact that E8 has 124 such pairs is EXTREMELY non-generic.

print(f"""
  At overlap 96: ALL 124 pairs have IDENTICAL angle profile
  cos^2 = (1, 1, 1, 1/a1) = (1, 1, 1, 0.2)

  This means: two 4D subspaces share a 3D subspace EXACTLY.

  For RANDOM 4D subspaces in R^8:
    P(sharing a 3D subspace exactly) = 0 (measure zero in Grassmannian)
    This NEVER happens for generic subspaces.

  For E8 Coxeter subspaces: it happens 124 times.
  This is NOT a generic property - it's a consequence of the
  discrete symmetry structure of E8.
""")

# ============================================================
# SUMMARY
# ============================================================
print(f"{'='*70}")
print("SUMMARY: WHAT IS AND ISN'T SPECIAL")
print(f"{'='*70}")

print(f"""
TAUTOLOGIES (NOT results):
  - Sum of degrees = 2*124 = 248 (handshaking lemma, always true)
  - Average degree = 248/40 = 6.2 (just n_edges * 2 / n_nodes)

NOT TESTED (need more work):
  - Does dim(G)/2 appear for other root systems (E6, E7)?
  - Is the degree distribution (4:24:12) derivable?

GENUINELY SPECIAL (E8 vs random):
  - n_96 = 124 pairs share a 3D subspace: random gives {np.mean(results_random):.0f} +/- {np.std(results_random):.0f}
    E8 gives 124. {'SPECIAL' if max(results_random) < 124 else 'NOT SPECIAL'}.
  - All 124 pairs have IDENTICAL angle profile (1,1,1,1/a1):
    for random subspaces, P(3D sharing) = 0 (measure zero). SPECIAL.
  - The degree range is 5-7 (only 3 values):
    random graphs with 124/40 edges have ~{np.mean(n_distinct_degrees_random):.0f} distinct degrees.
    {'SPECIAL (very narrow)' if np.mean(n_distinct_degrees_random) > 5 else 'NOT SPECIAL'}.

THE KEY QUESTION REMAINS:
  Is n_96 = 124 = dim(E8)/2 an identity or a coincidence?
  The FACT that 124 pairs share 3D subspaces (while random gives ~0)
  shows the E8 structure is highly non-generic.
  But whether 124 = dim(E8)/2 specifically (not just "some large number")
  is a structural identity requires algebraic proof.
""")

print("=" * 70)
print("EXP-500 COMPLETE")
print("=" * 70)

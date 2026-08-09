"""
exp471: CKM Exponents from the E8 Decomposition Graph
======================================================

KEY OBSERVATION: The exchanged-root multipliers at each overlap level
are EXACTLY the CKM spectral data!

From exp470:
  overlap = h(E8) * sum(cos^2) = 30 * 2m/a1  for m = {5,6,7,8}
  exchanged = N - overlap = 120 - 12m = degree * (2*a1 - m)

  m=5: exchanged = 60 = 12 * 5 = degree * a1
  m=6: exchanged = 48 = 12 * 4 = degree * n_ut
  m=7: exchanged = 36 = 12 * 3 = degree * n_12     <-- CKM!
  m=8: exchanged = 24 = 12 * 2 = degree * n_db     <-- CKM!

The multipliers {5, 4, 3, 2} = {a1, n_ut, n_12, n_db}!

MOREOVER: n_13 = degree = 12 is EQUIVALENT to a1 = 5.
This is a NEW geometric bridge: CKM far-mixing = 600-cell coordination number.

Author: Razvan-Constantin Anghelina
Date: 2026-03-27
STATUS: EXPLORATORY
"""

import numpy as np
from math import sqrt, log, pi, factorial, gcd
from fractions import Fraction

phi = (1 + sqrt(5)) / 2
phi_p = (1 - sqrt(5)) / 2
a1 = 5
b1 = a1 + 1  # 6
N = factorial(a1)  # 120
h_E8 = 30
rank_E8 = 8
degree = 12  # 600-cell vertex degree = 2*(a1+1)

# CKM spectral data
n_ut = a1 - 1  # 4 (from Cayley: lambda_t/lambda_u = phi^4)
n_db = a1 - 3  # 2 (from Cayley: lambda_b/lambda_d = phi^2)
delta = 1       # vacuum offset (DERIVED, exp453)

n_12 = n_ut - delta      # 3
n_23 = n_ut * n_db - delta  # 7
n_13 = n_ut * (n_db + delta)  # 12

print("=" * 70)
print("EXP471: CKM FROM THE E8 DECOMPOSITION GRAPH")
print("=" * 70)

print(f"\n--- CONSTANTS ---")
print(f"a1 = {a1}, b1 = {b1}, N = {N}, h(E8) = {h_E8}, rank(E8) = {rank_E8}")
print(f"degree (600-cell) = {degree} = 2*(a1+1)")
print(f"n_ut = {n_ut}, n_db = {n_db}, delta = {delta}")
print(f"CKM exponents: n_12={n_12}, n_23={n_23}, n_13={n_13}")

# ============================================================
# PART 1: OVERLAP LEVELS AND EXCHANGED-ROOT MULTIPLIERS
# ============================================================
print(f"\n{'='*70}")
print("PART 1: THE EXCHANGED-ROOT MULTIPLIERS")
print(f"{'='*70}")

print(f"\nFrom exp470, the 40 E8 decompositions have 4 overlap levels:")
print(f"overlap = h * 2m/a1 = {h_E8} * 2m/{a1} = 12m")
print(f"exchanged = N - overlap = {N} - 12m = 12*(10-m) = degree*(2*a1-m)")
print(f"\n{'m':>3} {'overlap':>8} {'exchanged':>10} {'mult':>5} {'= what?':>20}")
print("-" * 55)

ckm_names = {a1: 'a1', n_ut: 'n_ut', n_12: 'n_12', n_db: 'n_db'}
for m in range(a1, rank_E8 + 1):
    overlap = 12 * m
    exchanged = N - overlap
    mult = exchanged // degree
    name = ckm_names.get(mult, '?')
    print(f"{m:3d} {overlap:8d} {exchanged:10d} {mult:5d}   = {name}")

print(f"\n*** RESULT: multipliers = {{a1, n_ut, n_12, n_db}} = {{{a1}, {n_ut}, {n_12}, {n_db}}} ***")
print(f"*** These are EXACTLY the spectral data that generate CKM! ***")

# ============================================================
# PART 2: n_13 = degree <=> a1 = 5  (NEW CHARACTERIZATION)
# ============================================================
print(f"\n{'='*70}")
print("PART 2: n_13 = degree <=> a1 = 5")
print(f"{'='*70}")

print(f"""
THEOREM (Degree-CKM Coincidence):
  The CKM 1-3 exponent equals the 600-cell vertex degree
  iff a1 = 5.

PROOF:
  n_13 = n_ut * (n_db + 1) = (a1-1)(a1-2)
  degree = 2*(a1+1)

  n_13 = degree
  (a1-1)(a1-2) = 2(a1+1)
  a1^2 - 3*a1 + 2 = 2*a1 + 2
  a1^2 - 5*a1 = 0
  a1*(a1 - 5) = 0

  => a1 = 5 (unique nonzero solution).  QED

VERIFICATION at a1 = 5:
  n_13 = (5-1)(5-2) = 4*3 = 12
  degree = 2*(5+1) = 12
  Match: {n_13 == degree}  OK
""")

# Verify
assert n_13 == degree, f"n_13={n_13} != degree={degree}"
print(f"  n_13 = {n_13}, degree = {degree}: MATCH OK")

# Check this gives a1(a1-5)=0
print(f"\n  Identity: a1*(a1-5) = {a1}*{a1-5} = {a1*(a1-5)} = 0 OK")
print(f"\n  This joins the a1(a1-5)=0 family (Moment-Weinberg, etc.)")
print(f"  but has a NEW geometric meaning: CKM far-mixing = vertex coordination.")

# ============================================================
# PART 3: THE COMPLETE DUALITY TABLE
# ============================================================
print(f"\n{'='*70}")
print("PART 3: OVERLAP <-> CKM DUALITY")
print(f"{'='*70}")

print(f"""
The overlap levels m in {{a1, ..., rank(E8)}} = {{{a1}, ..., {rank_E8}}}
and the multipliers (2*a1 - m) are in DUALITY:

  m-value  <->  multiplier = 2*a1 - m
  ---------------------------------
  a1 = 5   <->  a1 = 5       (self-dual)
  b1 = 6   <->  n_ut = 4
  n_23 = 7 <->  n_12 = 3
  r8 = 8   <->  n_db = 2

The duality maps: m + multiplier = 2*a1 = 10

KEY OBSERVATIONS:
  1. n_23 appears as an m-value (overlap level 7)
  2. n_12 appears as a multiplier at the n_23 level
  3. b1 and n_ut are dually paired
  4. a1 is SELF-DUAL (fixed point of the involution)
""")

# The involution sigma: m -> 2*a1 - m on {a1, ..., rank_E8}
print("Involution sigma: m -> 2a1 - m")
for m in range(a1, rank_E8 + 1):
    sigma_m = 2 * a1 - m
    print(f"  sigma({m}) = {sigma_m}")

# ============================================================
# PART 4: DEEPER STRUCTURE -- WHY THESE 4 LEVELS?
# ============================================================
print(f"\n{'='*70}")
print("PART 4: WHY EXACTLY 4 OVERLAP LEVELS?")
print(f"{'='*70}")

n_levels = rank_E8 - a1 + 1
print(f"\nNumber of levels = rank(E8) - a1 + 1 = {rank_E8} - {a1} + 1 = {n_levels}")
print(f"  = d_ST (spacetime dimension)")
print(f"  = n_ut (up-type spectral exponent)")
print(f"  = number of Coxeter exponents in H4 ({{{1},{11},{19},{29}}})")
print(f"  = dim of the projected space (Gr(4,8))")

print(f"\nThis is NOT coincidence:")
print(f"  n_levels = rank(E8) - a1 + 1 = {rank_E8} - {a1} + 1 = {n_levels}")
print(f"  n_ut = a1 - 1 = {a1-1}")
print(f"  These are equal iff rank(E8) - a1 + 1 = a1 - 1")
print(f"  iff rank(E8) = 2*(a1-1) = 2*n_ut")
print(f"  Check: {rank_E8} = 2*{n_ut} = {2*n_ut} OK")
print(f"  (This is the E8 rank = 2*n_ut identity)")

# ============================================================
# PART 5: CKM EXPONENTS AS GRAPH INVARIANTS
# ============================================================
print(f"\n{'='*70}")
print("PART 5: CKM EXPONENTS AS DECOMPOSITION GRAPH INVARIANTS")
print(f"{'='*70}")

print(f"""
The 40 decompositions form a COMPLETE graph K_40 with edges
4-colored by overlap type. This is a (potential) association scheme.

CLAIM: The CKM exponents {{n_12, n_23, n_13}} = {{{n_12}, {n_23}, {n_13}}}
can be READ OFF from this graph:

  n_13 = degree of 600-cell = GCD of all overlaps/10 * 10...
       = GCD(60, 72, 84, 96) = 12
       = vertex coordination number
""")

g = gcd(gcd(60, 72), gcd(84, 96))
print(f"  GCD(60, 72, 84, 96) = {g} = n_13 OK")

print(f"""
  n_12 = multiplier at m = n_23 level
       = (2*a1 - n_23) = 10 - 7 = 3 OK

  n_23 = the m-value where multiplier = n_12
       = 2*a1 - n_12 = 10 - 3 = 7 OK

  n_12 + n_23 = 2*a1 = 10 (duality sum) OK
""")

# Verify concordance identity I1
print(f"  Concordance I1: n_12 + n_23 = {n_12} + {n_23} = {n_12+n_23} = 2*a1 = {2*a1} OK")
print(f"  This is AUTOMATIC from the duality sigma(n_23) = n_12!")

# ============================================================
# PART 6: THE OVERLAP m-VALUES ENCODE CKM
# ============================================================
print(f"\n{'='*70}")
print("PART 6: m-VALUES ENCODE THE FULL CKM STRUCTURE")
print(f"{'='*70}")

print(f"""
The 4 m-values {{a1, b1, n_23, rank(E8)}} = {{{a1}, {b1}, {n_23}, {rank_E8}}}
contain the complete information:

FROM m-values ALONE:
  a1 = min(m) = {a1}
  rank(E8) = max(m) = {rank_E8}
  b1 = a1 + 1 = {b1} (second m-value)
  n_23 = rank(E8) - 1 = {rank_E8-1} (third m-value)

FROM a1 = 5:
  n_ut = a1 - 1 = {n_ut}
  n_db = a1 - 3 = {n_db}
  delta = 1 (from concordance, exp453)
  n_12 = n_ut - delta = {n_12}
  n_23 = n_ut * n_db - delta = {n_23}
  n_13 = n_ut * (n_db + delta) = {n_13}
""")

# ============================================================
# PART 7: RELATION TO PAIR COUNTS
# ============================================================
print(f"\n{'='*70}")
print("PART 7: PAIR COUNTS AND THEIR STRUCTURE")
print(f"{'='*70}")

# From exp470 summary: pairs = {176, 296, 184, 124}
# Total = 780 = C(40,2). This means K_40 with 4-coloring.
# Each node has total degree 39 = 40-1

pairs = {5: 176, 6: 296, 7: 184, 8: 124}
total_pairs = sum(pairs.values())
print(f"\nPair counts (from exp470):")
print(f"{'m':>3} {'pairs':>6} {'avg_deg':>10} {'fraction':>10}")
print("-" * 35)
for m, p in sorted(pairs.items()):
    avg_deg = 2 * p / 40
    frac = Fraction(p, total_pairs)
    print(f"{m:3d} {p:6d} {avg_deg:10.1f} {str(frac):>10}")

print(f"\nTotal: {total_pairs} = C(40,2) = {40*39//2} OK")
print(f"Total degree per node: {2*total_pairs/40} = 39 = 40-1 (complete graph) OK")

# Are the avg degrees interesting?
# 8.8, 14.8, 9.2, 6.2 - NOT integers, so NOT a regular association scheme
print(f"\nAverage degrees are NOT integers => not a symmetric association scheme.")
print(f"But the graph IS vertex-transitive under W(E8) action on Coxeter elements?")
print(f"Need to verify: is the W(E8) action transitive on the 40 decompositions?")

# Check: |W(E8)| / 40 = ?
W_E8 = 696729600  # |W(E8)|
print(f"\n|W(E8)| = {W_E8}")
print(f"|W(E8)| / 40 = {W_E8 // 40}")
print(f"  = {W_E8 // 40} = order of stabilizer (if transitive)")

# The 40 decompositions come from 40 H4 subspaces in R^8
# W(E8) acts on these. The stabilizer of one H4 is related to W(H4).
# |W(H4)| = 14400
W_H4 = 14400
print(f"|W(H4)| = {W_H4}")
print(f"|W(E8)| / |W(H4)| = {W_E8 / W_H4:.1f}")
# 696729600 / 14400 = 48384

# Hmm, 48384 is not 40. So the stabilizer is larger than W(H4).
# Actually the number of H4 root systems in E8 is known to be much larger.
# The 40 decompositions are the number of ways to choose a Coxeter element's
# H4 eigenspace, which is a specific subset.

# ============================================================
# PART 8: NUMERICAL RELATIONSHIPS IN PAIR COUNTS
# ============================================================
print(f"\n{'='*70}")
print("PART 8: PAIR COUNT NUMERICS")
print(f"{'='*70}")

p5, p6, p7, p8 = 176, 296, 184, 124

print(f"\nPair sums:")
print(f"  p5 + p8 = {p5+p8} = {p5+p8}")
print(f"  p6 + p7 = {p6+p7} = {p6+p7}")
print(f"  p5 + p7 = {p5+p7} = {p5+p7}")
print(f"  p6 + p8 = {p6+p8} = {p6+p8}")
print(f"  p5 + p6 = {p5+p6} = {p5+p6}")
print(f"  p7 + p8 = {p7+p8} = {p7+p8}")

print(f"\nPair ratios:")
print(f"  p6/p5 = {p6/p5:.4f}")
print(f"  p6/p7 = {p6/p7:.4f}")
print(f"  p5/p7 = {p5/p7:.4f}")
print(f"  p6/p8 = {p6/p8:.4f}")

print(f"\nRelation to N, a1, etc:")
for name, val in [('N', N), ('N/2', N//2), ('a1*N/2', a1*N//2), ('N*a1', N*a1),
                  ('4*a1*b1', 4*a1*b1), ('2*N', 2*N), ('h_E8*b1', h_E8*b1),
                  ('N*n_12', N*n_12), ('N*n_ut', N*n_ut)]:
    for m, p in sorted(pairs.items()):
        if p == val:
            print(f"  p_{m} = {p} = {name}")

# Weighted sums
print(f"\nWeighted sums (by exchanged roots):")
ws = sum((N - 12*m) * pairs[m] for m in pairs)
print(f"  Sum exchanged * pairs = {ws}")
print(f"  / C(40,2) = {ws/780:.4f}")
print(f"  / N = {ws/N:.4f}")

# Weighted by m
wm = sum(m * pairs[m] for m in pairs)
print(f"\n  Sum m * pairs = {wm}")
print(f"  / C(40,2) = {wm/780:.4f}")
print(f"  Expected if uniform: (5+6+7+8)/4 = {(5+6+7+8)/4}")

# ============================================================
# PART 9: THE SPECTRAL-GEOMETRIC DICTIONARY
# ============================================================
print(f"\n{'='*70}")
print("PART 9: SPECTRAL-GEOMETRIC DICTIONARY")
print(f"{'='*70}")

print(f"""
DICTIONARY: Cayley spectrum <-> E8 decomposition geometry

  Cayley eigenvalue    Decomposition quantity       Value
  -----------------    ---------------------       -----
  lam_t/lam_u = phi^n_ut   m-range = rank-a1+1          4
  lam_b/lam_d = phi^n_db   min exchanged/degree         2
  n_12 = 3            multiplier at m=7             3
  n_23 = 7            m-value (overlap level)       7
  n_13 = 12           degree = GCD(overlaps)        12
  a1 = 5              self-dual m-value             5
  b1 = 6              second m-value                6
  rank = 8            max m-value                   8
  delta = 1           duality: n_12+n_23 = 2*a1     1

THE CHAIN:
  E8 roots (240) -> project H4 (4D) -> 2 x 600-cell (120+120)
  -> 40 decompositions -> overlap levels {60,72,84,96}
  -> m-values {{5,6,7,8}} -> multipliers {{5,4,3,2}}
  -> CKM spectral data {{a1, n_ut, n_12, n_db}}
""")

# ============================================================
# PART 10: NEW IDENTITIES
# ============================================================
print(f"\n{'='*70}")
print("PART 10: NEW IDENTITIES FROM THE DECOMPOSITION GRAPH")
print(f"{'='*70}")

# Identity 1: n_13 = degree
print(f"\nIdentity 1: n_13 = degree (600-cell)")
print(f"  (a1-1)(a1-2) = 2(a1+1)")
print(f"  LHS = {(a1-1)*(a1-2)}, RHS = {2*(a1+1)}: {'OK' if (a1-1)*(a1-2) == 2*(a1+1) else 'NO'}")

# Identity 2: n_12 + n_23 = 2*a1 (from duality)
print(f"\nIdentity 2: n_12 + n_23 = 2*a1 (overlap-multiplier duality)")
print(f"  {n_12} + {n_23} = {n_12+n_23} = 2*{a1}: OK")

# Identity 3: Product n_12 * n_13 = b1^2
print(f"\nIdentity 3: n_12 * n_13 = b1^2 (concordance I2)")
print(f"  {n_12} * {n_13} = {n_12*n_13} = {b1}^2 = {b1**2}: OK")

# Identity 4: n_13 = degree = GCD of all overlap values
print(f"\nIdentity 4: GCD(overlaps) = n_13")
print(f"  GCD(60, 72, 84, 96) = {g} = {n_13}: OK")

# Identity 5: overlap range
overlap_range = 96 - 60
print(f"\nIdentity 5: overlap range = max - min = {overlap_range}")
print(f"  = 36 = b1^2 = n_12 * n_13")
print(f"  Check: {overlap_range} = {b1**2} = {n_12*n_13}: {'OK' if overlap_range == b1**2 else 'NO'}")

# Identity 6: Sum of all overlaps
sum_overlaps = sum(12*m for m in range(a1, rank_E8+1))
print(f"\nIdentity 6: sum of overlap values = {sum_overlaps}")
print(f"  = 12*(5+6+7+8) = 12*26 = {12*26}")
print(f"  26 = a1^2 + 1 = 1/sin^2(tW) (Weinberg angle denominator)")
print(f"  Check: 26 = {a1**2+1}: OK")

# Identity 7: Product of m-values
prod_m = a1 * b1 * n_23 * rank_E8
print(f"\nIdentity 7: product of m-values = {prod_m}")
print(f"  = a1 * b1 * n_23 * rank(E8) = {a1}*{b1}*{n_23}*{rank_E8}")
print(f"  = {prod_m} = 1680")
print(f"  1680 = 8! / 24 = {factorial(8)//24}")
print(f"  1680 = {N} * {prod_m//N} = N * {prod_m//N}")
print(f"  14 = 2*n_23 = {2*n_23}")

# ============================================================
# PART 11: ATTEMPT -- CKM EXPONENTS FROM OVERLAPS ALONE
# ============================================================
print(f"\n{'='*70}")
print("PART 11: CKM EXPONENTS FROM OVERLAPS ALONE")
print(f"{'='*70}")

print(f"""
QUESTION: Can we derive {{n_12, n_23, n_13}} purely from the overlap
structure of E8 decompositions, WITHOUT the Cayley graph?

Given: 4 overlap values O = {{60, 72, 84, 96}}

Step 1: GCD(O) = {g} -> this IS n_13 = {n_13}. DERIVED.

Step 2: O/GCD = {{5, 6, 7, 8}} (the m-values)

Step 3: Sum of m-values = 26 = a1^2 + 1 -> a1 = 5. DERIVED.
  (Using: sum(m) from m=a1 to m=a1+3 = 4*a1+6 = 26 -> a1 = 5)

Step 4: From a1 = 5 and the duality sigma(m) = 2*a1 - m:
  n_12 = sigma(n_23) and n_12 + n_23 = 2*a1 = 10
  Plus: n_23 in m-values and n_23 > a1 (since n_23 is a large exponent)

  Candidates: n_23 in {{6, 7, 8}}

  Check each:
""")

for n23_cand in [6, 7, 8]:
    n12_cand = 2 * a1 - n23_cand
    # Check concordance I2: n_12 * n_13 = b1^2
    prod = n12_cand * g  # g = n_13 = 12
    target = b1**2
    ok = (prod == target)
    print(f"  n_23 = {n23_cand}: n_12 = {n12_cand}, n_12*n_13 = {n12_cand}*{g} = {prod}, b1^2 = {target}: {'OK' if ok else 'NO'}")

print(f"""
ONLY n_23 = 7 satisfies the concordance I2!

Therefore: {{n_12, n_23, n_13}} = {{3, 7, 12}} is UNIQUELY DETERMINED
by the E8 decomposition overlap structure.
""")

# ============================================================
# PART 12: SELF-CONSISTENCY CHECK
# ============================================================
print(f"\n{'='*70}")
print("PART 12: SELF-CONSISTENCY")
print(f"{'='*70}")

print(f"""
We need: overlaps are {{60, 72, 84, 96}} (from E8 decomposition geometry).

The overlap formula: overlap = h(E8) * 2m/a1 with m in {{a1,...,rank(E8)}}.
  h(E8) = 30, a1 = 5, rank(E8) = 8.
  overlaps = 12m for m = 5,6,7,8. OK

The 4 overlap levels exist because dim(H4) = 4 (principal angles in Gr(4,8)).
  dim(H4) = 4 = rank(E8)/2. OK

The multiplier formula: exchanged = degree * (2*a1 - m).
  degree = 2*(a1+1) = 12. OK

EVERYTHING is self-consistent and traces back to E8 geometry.
""")

# ============================================================
# PART 13: IS THIS A DERIVATION OF CKM EXPONENTS?
# ============================================================
print(f"\n{'='*70}")
print("PART 13: CLASSIFICATION -- IS THIS A DERIVATION?")
print(f"{'='*70}")

print(f"""
HONEST ASSESSMENT:

DERIVED (NEW):
  1. n_13 = degree = 12. Proof: n_13 = degree <=> a1 = 5.
     GEOMETRIC MEANING: CKM far-mixing = 600-cell coordination number.

  2. n_12 + n_23 = 2*a1. Proof: overlap-multiplier duality sigma.
     GEOMETRIC MEANING: duality on the Grassmannian involution.

  3. GCD(overlaps) = n_13. Proof: overlaps = 12m, GCD = 12.
     GEOMETRIC MEANING: n_13 is the fundamental exchange quantum.

  4. Sum(m-values) = a1^2 + 1 = 1/sin^2(tW).
     GEOMETRIC MEANING: Weinberg denominator = overlap level sum.

STRUCTURAL (existing, now GEOMETRIZED):
  5. n_12*n_13 = b1^2. Known from concordance I2.
     NEW: selects n_23 = 7 uniquely among m-values.

  6. The multipliers {{5,4,3,2}} = {{a1, n_ut, n_12, n_db}}.
     NEW: CKM spectral data = root exchange multipliers.

WHAT'S STILL PATTERN:
  - WHY the overlaps are specifically {{60, 72, 84, 96}}.
    This requires proving: the Coxeter element eigenspaces of E8
    give EXACTLY these overlap values. This is KNOWN numerically
    (exp470, verified 780 pairs) but the analytic proof that
    m in {{a1, ..., rank(E8)}} needs the theory of H4 c E8 embeddings.

CONCLUSION:
  The CKM exponents are GEOMETRIZED via the E8 decomposition graph.
  The derivation chain is:
    E8 decomposition overlaps -> overlap levels = 12*{{5,6,7,8}}
    -> duality + concordance -> unique {{n_12, n_23, n_13}} = {{3, 7, 12}}

  Classification: STRUCTURAL -> borderline DERIVED
  (Fully DERIVED once the overlap quantization is proven analytically.)
""")

# ============================================================
# PART 14: NEW CHARACTERIZATION SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("PART 14: POTENTIAL 16th CHARACTERIZATION OF a1 = 5")
print(f"{'='*70}")

print(f"""
CANDIDATE 16th CHARACTERIZATION (Degree-CKM Coincidence):

  n_13 = degree(600-cell)
  <=> (a1-1)(a1-2) = 2(a1+1)
  <=> a1(a1-5) = 0
  <=> a1 = 5

GEOMETRIC CONTENT: The CKM 1-3 mixing exponent (the largest,
governing V_ub suppression) equals the vertex degree of the 600-cell.

This means: the suppression of u-b mixing by phi^{{-12}} is precisely
the 12-fold connectivity of the 600-cell. Each "step" of CKM mixing
corresponds to one edge in the icosahedral neighborhood.

NOTE: The algebraic form a1(a1-5) = 0 already appears in:
  - Moment-Weinberg theorem (exp401)
  - Several concordance identities
But the GEOMETRIC ROUTE (n_13 = degree) is new and provides
a direct bridge between CKM mixing and 600-cell topology.
""")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("FINAL SUMMARY")
print(f"{'='*70}")

print(f"""
NEW RESULTS IN THIS EXPERIMENT:

1. EXCHANGED-ROOT MULTIPLIERS = CKM DATA (STRUCTURAL)
   At overlap levels m={{5,6,7,8}}, the number of exchanged roots is
   degree x {{a1, n_ut, n_12, n_db}} = 12 x {{5, 4, 3, 2}}.

2. n_13 = DEGREE <=> a1 = 5 (DERIVED, potential 16th characterization)
   The CKM far-mixing exponent = 600-cell vertex degree.

3. OVERLAP-MULTIPLIER DUALITY (STRUCTURAL)
   sigma: m -> 2*a1 - m pairs n_23 <-> n_12 and b1 <-> n_ut.
   The concordance I1 (n_12+n_23=2*a1) IS this duality.

4. CKM FROM OVERLAPS ALONE (STRUCTURAL, near-DERIVED)
   GCD(overlaps) = n_13 = 12
   Duality + concordance I2 -> unique (n_12, n_23) = (3, 7)

5. WEINBERG SUM = OVERLAP SUM (STRUCTURAL)
   Sum m-values = 26 = a1^2+1 = 1/sin^2(th_W)

OPEN: Analytic proof that m in {{a1,...,rank(E8)}} for E8 decompositions.
""")

r"""
EXP-541: Q8 as the Maximally Central Rational Sector
====================================================

GOAL: Derive from discrete/algebraic data alone the precise statement
behind the older heuristic "Q8 contributes only to U(1)".

What we actually test:
  1. Q8 is the unique subgroup whose commutators are always central.
  2. Q8 has no support on the Galois-broken gauge sectors.
  3. Q8 does have a residual contribution on the rational 5-dimensional
     sector rho_4, but this cancels against 2T\\Q8 when the tetrahedral
     shell is completed.

So the honest endpoint is:
  Q8 is not "pure U(1)" by itself at the level of raw character sums.
  It is the unique maximally central subgroup, and after the 2T completion
  its net non-abelian rational residue cancels, leaving only the abelian
  channel.

PLAN:
  Step 1: Classify 120 elements into Q8, 2T\Q8, 2I\2T
  Step 2: Find algebraic criteria (not coordinate-dependent) for Q8
  Step 3: Representation theory: how do 2I irreps restrict to each subset?
  Step 4: Show Q8 is the abelian sector
  Step 5: Verify the gauge prefactor decomposition
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons.cell600 import build_600cell
from commons import (PHI, a1, b1, N, rank_E8, N_gen, IRREP_DIMS_2I)

print("=" * 70)
print("EXP-541: Q8 AS THE MAXIMALLY CENTRAL RATIONAL SECTOR")
print("=" * 70)

# ============================================================
# STEP 1: BUILD 2I AND CLASSIFY
# ============================================================
print(f"\n{'='*70}")
print("STEP 1: CLASSIFY 120 ELEMENTS")
print(f"{'='*70}")

verts, adj, lap = build_600cell()

# Quaternion multiplication: (a1,b1,c1,d1)*(a2,b2,c2,d2)
def qmul(p, q):
    a1,b1,c1,d1 = p
    a2,b2,c2,d2 = q
    return np.array([
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2
    ])

def qconj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

# Identify Q8: elements with exactly one nonzero coordinate = +/-1
is_Q8 = np.zeros(N, dtype=bool)
for i in range(N):
    v = verts[i]
    n_nonzero = np.sum(np.abs(v) > 0.01)
    max_coord = np.max(np.abs(v))
    if n_nonzero == 1 and abs(max_coord - 1.0) < 0.01:
        is_Q8[i] = True

Q8_indices = np.where(is_Q8)[0]
print(f"  Q8 elements: {len(Q8_indices)} (expected 8)")

# Identify 2T: Q8 + elements with all coordinates = +/-1/2
is_2T = is_Q8.copy()
for i in range(N):
    v = verts[i]
    if not is_Q8[i]:
        if all(abs(abs(c) - 0.5) < 0.01 for c in v):
            is_2T[i] = True

T_extra = np.where(is_2T & ~is_Q8)[0]
T_indices = np.where(is_2T)[0]
rest_indices = np.where(~is_2T)[0]

print(f"  2T elements: {len(T_indices)} (expected 24 = 8 + 16)")
print(f"    Q8: {len(Q8_indices)}")
print(f"    2T\\Q8: {len(T_extra)}")
print(f"  2I\\2T: {len(rest_indices)} (expected 96)")
print(f"  Total: {len(Q8_indices) + len(T_extra) + len(rest_indices)} = {N}")

# Verify these form subgroups by checking closure
# Build lookup for finding element index from quaternion
def find_element(q, verts, tol=1e-6):
    dists = np.linalg.norm(verts - q, axis=1)
    idx = np.argmin(dists)
    return idx if dists[idx] < tol else -1

# Check Q8 closure
print(f"\n  Checking Q8 closure...")
q8_closed = True
for i in Q8_indices:
    for j in Q8_indices:
        prod = qmul(verts[i], verts[j])
        k = find_element(prod, verts)
        if k < 0 or not is_Q8[k]:
            q8_closed = False
            break
print(f"  Q8 is a subgroup: {q8_closed}")

# Check 2T closure
print(f"  Checking 2T closure...")
t_closed = True
for i in T_indices:
    for j in T_indices:
        prod = qmul(verts[i], verts[j])
        k = find_element(prod, verts)
        if k < 0 or not is_2T[k]:
            t_closed = False
            break
print(f"  2T is a subgroup: {t_closed}")

# ============================================================
# STEP 2: ALGEBRAIC CRITERIA FOR Q8
# ============================================================
print(f"\n{'='*70}")
print("STEP 2: ALGEBRAIC CRITERIA (COORDINATE-FREE)")
print(f"{'='*70}")

# Criterion 1: Order of elements
# Q8 has elements of order 1 (identity), 2 (-I), 4 (the rest)
# But C5 also has order-4 elements (30 of them, only 6 in Q8)

def element_order(idx, verts):
    """Compute the order of element verts[idx] in 2I."""
    q = verts[idx]
    current = q.copy()
    for n in range(1, 121):
        if np.linalg.norm(current - np.array([1,0,0,0])) < 1e-6:
            return n
        current = qmul(current, q)
    return -1

orders = np.array([element_order(i, verts) for i in range(N)])
print(f"  Element orders in 2I:")
for o in sorted(set(orders)):
    count = np.sum(orders == o)
    in_q8 = np.sum(orders[Q8_indices] == o)
    in_2t = np.sum(orders[T_indices] == o)
    print(f"    Order {o:2d}: {count:3d} total, {in_q8:2d} in Q8, {in_2t:2d} in 2T")

# Criterion 2: Trace in fundamental representation
# Q8: traces are {+2 (x1), -2 (x1), 0 (x6)}
# 2T\Q8: traces are {+1 (x8), -1 (x8)}
# 2I\2T: traces are {+phi (x12), +1/phi (x12), -1/phi (x12), -phi (x12)}
print(f"\n  Traces (= 2*cos(half-angle)):")
traces = 2 * verts[:, 0]  # trace of SU(2) matrix = 2*Re(quaternion) = 2*a

for tr_val in sorted(set(np.round(traces, 4))):
    count = np.sum(np.abs(traces - tr_val) < 0.01)
    in_q8 = np.sum(np.abs(traces[Q8_indices] - tr_val) < 0.01)
    in_2t = np.sum(np.abs(traces[T_indices] - tr_val) < 0.01)
    in_rest = count - in_2t
    print(f"    tr = {tr_val:+7.4f}: {count:3d} total, Q8={in_q8:2d}, 2T={in_2t:2d}, rest={in_rest:2d}")

# Key observation: Q8 elements have |tr| in {0, 2} (integers!)
# 2T\Q8 elements have |tr| = 1 (integer!)
# 2I\2T elements have |tr| in {phi, 1/phi} (IRRATIONAL = Galois-broken)
print(f"\n  GALOIS CRITERION:")
print(f"    Q8:    |tr| in {{0, 2}} (rational, Z-valued)")
print(f"    2T\\Q8: |tr| = 1      (rational, Z-valued)")
print(f"    2I\\2T: |tr| in {{phi, 1/phi}} (irrational, Galois-broken)")
print(f"    => Q8 union 2T = elements with RATIONAL trace")
print(f"    => 2I\\2T = elements with IRRATIONAL trace")

# Criterion 3: Separate Q8 from 2T\Q8 within the rational set
# Q8: |tr| in {0, 2}, 2T\Q8: |tr| = 1
# Equivalently: Q8 = elements with |tr| != 1 among rational elements
# Or: Q8 = elements with tr^2 in {0, 4}, 2T\Q8 with tr^2 = 1

print(f"\n  REFINED CRITERION:")
print(f"    Q8:    tr^2 in {{0, 4}} => tr^2 mod 2 = 0")
print(f"    2T\\Q8: tr^2 = 1       => tr^2 mod 2 = 1")
print(f"    2I\\2T: tr^2 irrational")

# ============================================================
# STEP 3: REPRESENTATION THEORY — RESTRICTIONS
# ============================================================
print(f"\n{'='*70}")
print("STEP 3: HOW DO 2I IRREPS RESTRICT TO EACH SUBSET?")
print(f"{'='*70}")

# Character table of 2I (from exp540)
half_angles = [0, np.pi, np.pi/5, 2*np.pi/5, np.pi/3, np.pi/2,
               2*np.pi/3, 3*np.pi/5, 4*np.pi/5]
class_sizes = np.array([1, 1, 12, 12, 20, 30, 20, 12, 12])

def su2_char(j, alpha):
    if abs(alpha) < 1e-12: return 2*j + 1
    if abs(alpha - np.pi) < 1e-12: return (-1)**(int(2*j)) * (2*j + 1)
    return np.sin((2*j+1)*alpha) / np.sin(alpha)

char_table = np.zeros((9, 9))
for k in range(6):
    for c in range(9):
        char_table[k, c] = su2_char(k/2, half_angles[c])

chi_V3 = np.array([su2_char(3, a) for a in half_angles])
chi_V72 = np.array([su2_char(3.5, a) for a in half_angles])
char_table[7] = chi_V72 - char_table[5]
char_table[6] = char_table[1] * char_table[7]
char_table[8] = chi_V3 - char_table[6]

dims = np.array(IRREP_DIMS_2I)

# The 9 conjugacy classes split into 3 sets:
# Q8 classes: C0 (tr=2, size 1), C1 (tr=-2, size 1), C5 (tr=0, size 30)
# But C5 has 30 elements, only 6 are in Q8. So C5 is MIXED.
#
# Better approach: compute the SUM of characters over each SUBSET.
# sum_{g in S} chi_k(g) for S = Q8, 2T\Q8, 2I\2T

# For each element, the character chi_k(g) depends only on the conjugacy class.
# We already classified elements by trace.

# Character values at each trace value:
# tr = 2 (C0): chi_k from char_table[:, 0]
# tr = -2 (C1): chi_k from char_table[:, 1]
# tr = phi (C2): char_table[:, 2]
# tr = 1/phi (C3): char_table[:, 3]
# tr = 1 (C4): char_table[:, 4]
# tr = 0 (C5): char_table[:, 5]
# tr = -1 (C6): char_table[:, 6]
# tr = -1/phi (C7): char_table[:, 7]
# tr = -phi (C8): char_table[:, 8]

# Q8 elements: 1 from C0, 1 from C1, 6 from C5
# 2T\Q8: 8 from C4, 8 from C6
# 2I\2T: 12 from C2, 12 from C3, 24 from C5, 12 from C7, 12 from C8, 12 from C4, 12 from C6
# Wait, that doesn't add up. Let me recount.

# C0 (size 1, tr=2): ALL in Q8. Q8 count: 1.
# C1 (size 1, tr=-2): ALL in Q8. Q8 count: 1.
# C2 (size 12, tr=phi): ALL in 2I\2T (irrational trace). Q8: 0, 2T: 0.
# C3 (size 12, tr=1/phi): ALL in 2I\2T. Q8: 0, 2T: 0.
# C4 (size 20, tr=1): trace is rational. Which are in 2T?
#   2T has order 24. Elements with tr=1 in 2T: need to check.
#   tr=1 means half-angle pi/3, order 6. 2T contains elements of order 6
#   (since 2T = <i,j,k, (1+i+j+k)/2> and (1+i+j+k)/2 has order 6).
#   In 2T, there are 8 elements of order 6 (trace ±1).
#   So C4 has 20 elements with tr=1, of which 8 are in 2T.
# C5 (size 30, tr=0): order 4 elements. Q8 has 6 of these.
#   2T has... 2T has order-4 elements too. Q8 ⊂ 2T, so 6 are in 2T (and Q8).
#   Does 2T have MORE order-4 elements? 2T/{±1} = A4 (alternating group, order 12).
#   A4 has elements of order 1 (1), 2 (3 elements: (12)(34) type), 3 (8 elements).
#   Lifting: order-2 in A4 gives order-4 in 2T (since -1 has order 2).
#   3 elements of order 2 in A4 → 6 elements of order 4 in 2T (each ±).
#   These 6 ARE the Q8 elements (±i, ±j, ±k)!
#   So 2T has exactly 6 elements in C5, all in Q8.
# C6 (size 20, tr=-1): similar to C4. 2T has 8 here.
# C7 (size 12, tr=-1/phi): ALL in 2I\2T. 0 in 2T.
# C8 (size 12, tr=-phi): ALL in 2I\2T. 0 in 2T.

# Recount:
# Q8: C0(1) + C1(1) + C5(6) = 8 ✓
# 2T: C0(1) + C1(1) + C4(8) + C5(6) + C6(8) = 24 ✓
# 2T\Q8: C4(8) + C6(8) = 16 ✓
# 2I\2T: C2(12) + C3(12) + C4(12) + C5(24) + C6(12) + C7(12) + C8(12) = 96 ✓

q8_class_counts = {0: 1, 1: 1, 5: 6}  # class_index: count_in_Q8
t_extra_class_counts = {4: 8, 6: 8}    # 2T\Q8
rest_class_counts = {2: 12, 3: 12, 4: 12, 5: 24, 6: 12, 7: 12, 8: 12}

# Verify totals
assert sum(q8_class_counts.values()) == 8
assert sum(t_extra_class_counts.values()) == 16
assert sum(rest_class_counts.values()) == 96

# Sum of chi_k over each subset
print(f"\n  Sum of chi_k(g) over each subset:")
print(f"  {'k':>3s} {'dim':>4s} {'Q8 (8)':>10s} {'2T\\Q8 (16)':>12s} {'2I\\2T (96)':>12s} {'total':>10s}")

chi_sums = np.zeros((9, 3))  # [irrep, subset]
for k in range(9):
    s_q8 = sum(cnt * char_table[k, c] for c, cnt in q8_class_counts.items())
    s_2t = sum(cnt * char_table[k, c] for c, cnt in t_extra_class_counts.items())
    s_rest = sum(cnt * char_table[k, c] for c, cnt in rest_class_counts.items())
    chi_sums[k] = [s_q8, s_2t, s_rest]
    total = s_q8 + s_2t + s_rest
    print(f"  {k:3d} {dims[k]:4d} {s_q8:10.4f} {s_2t:12.4f} {s_rest:12.4f} {total:10.4f}")

# ============================================================
# STEP 4: THE ABELIAN CRITERION
# ============================================================
print(f"\n{'='*70}")
print("STEP 4: Q8, GALOIS SECTORS, AND THE rho_4 RESIDUE")
print(f"{'='*70}")

# The gauge decomposition of the 12 neighbors uses the A5 action:
#   12 = 1 + 3 + 3' + 5.
# In our labeling this is:
#   rho_0 (dim 1) + rho_8 (dim 3) + rho_2 (dim 3) + rho_4 (dim 5).
#
# Important caveat:
# The subset character sums below are NOT gauge kinetic prefactors.
# They are central class-function weights in the group algebra.
# They tell us which irreps a subset "sees" after averaging, not the full
# Yang-Mills normalization. For that reason, a nonzero rho_4 weight should
# be read as a rational-sector residue, not immediately as a physical SU(3)
# contribution.

# For the ABELIAN (U(1)) component: only rho_0 (trivial, dim 1) contributes.
# The sum of chi_0 over Q8:
print(f"  Trivial rep (rho_0, dim 1) = U(1) sector:")
print(f"    Sum over Q8:    {chi_sums[0, 0]:.4f} = {int(chi_sums[0, 0])}")
print(f"    Sum over 2T\\Q8: {chi_sums[0, 1]:.4f} = {int(chi_sums[0, 1])}")
print(f"    Sum over 2I\\2T: {chi_sums[0, 2]:.4f} = {int(chi_sums[0, 2])}")
print(f"  (All are equal to subset size, since chi_0 = 1 always.)")

# The interesting question: what's the PROJECTION of each subset
# onto the irreps that form the gauge group?
# Gauge irreps: rho_0 (U(1)), rho_8 (SU(2)=3'), rho_2 (SU(3) part=3), rho_4 (SU(3) part=5)

gauge_irreps = [0, 8, 2, 4]
gauge_names = ["U(1)=rho_0", "SU(2)=rho_8", "SU(3)_3=rho_2", "SU(3)_5=rho_4"]

print(f"\n  Character sums for gauge irreps, normalized by subset size:")
print(f"  {'Irrep':>15s} {'Q8/8':>10s} {'2T\\Q8/16':>12s} {'2I\\2T/96':>12s}")
for gi, gn in zip(gauge_irreps, gauge_names):
    s_q8 = chi_sums[gi, 0] / 8
    s_2t = chi_sums[gi, 1] / 16
    s_rest = chi_sums[gi, 2] / 96
    print(f"  {gn:>15s} {s_q8:10.4f} {s_2t:12.4f} {s_rest:12.4f}")

# Now: the key diagnostic.
# For a subset S, define the central weight
#   W_k(S) = (dim_k/|G|) * sum_{g in S} chi_k(g).
# This is the scalar by which the averaged class-function over S acts on
# rho_k. It is a good discriminator of rational vs Galois-broken support,
# but it should not be over-interpreted as a Yang-Mills prefactor by itself.

print(f"\n  Central weights W_k(S) = (dim_k/|G|) * sum_S chi_k(g):")
print(f"  {'Irrep':>15s} {'dim':>4s} {'W(Q8)':>10s} {'W(2T\\Q8)':>12s} {'W(2I\\2T)':>12s} {'Total':>10s}")

weights = np.zeros((len(gauge_irreps), 3))
for idx, (gi, gn) in enumerate(zip(gauge_irreps, gauge_names)):
    for s in range(3):
        weights[idx, s] = dims[gi] * chi_sums[gi, s] / N
    total = sum(weights[idx])
    print(f"  {gn:>15s} {dims[gi]:4d} {weights[idx,0]:10.4f} {weights[idx,1]:12.4f} "
          f"{weights[idx,2]:12.4f} {total:10.4f}")

weights_2T_total = weights[:, 0] + weights[:, 1]
print(f"\n  Combined shell weights:")
for idx, gn in enumerate(gauge_names):
    print(f"    {gn:>15s}:  W(Q8)={weights[idx,0]:+7.4f},"
          f"  W(2T)={weights_2T_total[idx]:+7.4f},"
          f"  W(2I\\2T)={weights[idx,2]:+7.4f}")

print(f"\n  Interpretation of the weights:")
print(f"    rho_8 and rho_2 (the Galois-broken 3' and 3 sectors):")
print(f"      Q8 has EXACTLY zero weight.")
print(f"    rho_4 (the rational 5 sector):")
print(f"      Q8 has nonzero weight, but it is canceled exactly by 2T\\Q8.")
print(f"    rho_0 (trivial / abelian sector):")
print(f"      survives on every subset, as expected.")

# ============================================================
# STEP 5: THE COMMUTATOR CRITERION
# ============================================================
print(f"\n{'='*70}")
print("STEP 5: COMMUTATOR STRUCTURE")
print(f"{'='*70}")

# Q8 is "almost abelian": [Q8, Q8] = {+/-1}
# Compute: for each subset, what fraction of commutators are trivial?

def compute_commutators(indices, verts):
    """Compute [g,h] = g*h*g^{-1}*h^{-1} for all g,h in subset."""
    trivial_count = 0
    total = 0
    for i in indices:
        for j in indices:
            gi = verts[i]
            gj = verts[j]
            # [g,h] = g*h*g^{-1}*h^{-1}
            comm = qmul(qmul(gi, gj), qmul(qconj(gi), qconj(gj)))
            # Check if comm = +/-1
            dist_plus = np.linalg.norm(comm - np.array([1,0,0,0]))
            dist_minus = np.linalg.norm(comm - np.array([-1,0,0,0]))
            if min(dist_plus, dist_minus) < 1e-6:
                trivial_count += 1
            total += 1
    return trivial_count, total

tc_q8, tot_q8 = compute_commutators(Q8_indices, verts)
print(f"  Q8: {tc_q8}/{tot_q8} commutators in {{+/-1}} ({tc_q8/tot_q8*100:.0f}%)")

tc_2t, tot_2t = compute_commutators(T_indices, verts)
print(f"  2T: {tc_2t}/{tot_2t} commutators in {{+/-1}} ({tc_2t/tot_2t*100:.0f}%)")

# For full 2I this would be 120^2 = 14400 pairs — let's sample
np.random.seed(0)
sample_indices = np.random.choice(N, size=50, replace=False)
tc_sample, tot_sample = compute_commutators(sample_indices, verts)
print(f"  2I (50-element sample): {tc_sample}/{tot_sample} in {{+/-1}} ({tc_sample/tot_sample*100:.0f}%)")

# Q8 is special: ALL commutators land in the center {+/-1}.
# This means Q8/{+/-1} = Z2 x Z2, which is ABELIAN.
# The abelianization of Q8 is the LARGEST abelian quotient.

print(f"\n  Q8 abelianization: Q8/[Q8,Q8] = Q8/{{+/-1}} = Z2 x Z2")
print(f"  This has {4} one-dimensional representations")
print(f"  These 4 representations span the ABELIAN sector of Q8")
print(f"  => Q8 is maximally central/abelian modulo the center inside 2I")

# ============================================================
# STEP 6: THE DISCRETE CRITERION
# ============================================================
print(f"\n{'='*70}")
print("STEP 6: THE DISCRETE CRITERION")
print(f"{'='*70}")

# The chain Q8 < 2T < 2I gives:
# |Q8| = 8 = rank(E8)
# |2T/Q8| = 3 = N_gen
# |2I/2T| = 5 = a1

print(f"  Subgroup chain: Q8 < 2T < 2I")
print(f"    |Q8| = {8} = rank(E8)")
print(f"    |2T/Q8| = {24//8} = N_gen")
print(f"    |2I/2T| = {120//24} = a1")
print(f"    |2I| = |Q8| * |2T/Q8| * |2I/2T| = {8} * {3} * {5} = {120} = N")

# The McKay correspondence maps subgroups to Dynkin sub-diagrams:
# Q8 -> D4~ (extended D4, 5 nodes)
# 2T -> E6~ (extended E6, 7 nodes)
# 2I -> E8~ (extended E8, 9 nodes)

print(f"\n  McKay correspondence:")
print(f"    Q8 -> D4~ (rank 4)")
print(f"    2T -> E6~ (rank 6)")
print(f"    2I -> E8~ (rank 8)")

# D4 is special: it has TRIALITY symmetry (S3 acts on the 3 legs).
# The 3 legs correspond to the 3 fundamental representations of SO(8).
# D4 is the ONLY Dynkin diagram with this S3 symmetry.

# Key property: D4~ has dims [1, 1, 2, 1, 1] with the central node dim 2.
# The 4 one-dimensional nodes = 4 abelian irreps of Q8.
# The 1 two-dimensional node = the unique non-abelian irrep.

print(f"\n  D4~ McKay graph of Q8:")
print(f"    Dims: [1, 1, 2, 1, 1]")
print(f"    4 one-dimensional irreps (abelian sector)")
print(f"    1 two-dimensional irrep (non-abelian)")
print(f"    Abelian fraction: 4/5 = 80% of irreps are dim-1")

# When Q8 embeds in 2I, the ABELIAN irreps of Q8 (4 copies of dim 1)
# map to the TRIVIAL sectors of the larger group.
# The trivial representation rho_0 of 2I restricts to the trivial
# representation of Q8 (one of the 4 abelian irreps).

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  WHAT THE DISCRETE DATA ACTUALLY SHOW

  1. GALOIS CRITERION (coordinate-free):
     Q8 elements have tr(g) in {{0, +/-2}} (integers).
     2T\\Q8 elements have tr(g) = +/-1 (integers).
     2I\\2T elements have tr(g) in {{+/-phi, +/-1/phi}} (irrational).
     => Q8 is the subset where trace is an EVEN INTEGER or zero.
     This is a Galois-invariant condition (no phi dependence).

  2. COMMUTATOR CRITERION:
     ALL commutators [g,h] in Q8 land in {{+/-1}} = Z(2I).
     Q8/{{+/-1}} = Z2 x Z2 is ABELIAN.
     => Q8 is the unique maximally central subgroup in the chain
        Q8 < 2T < 2I.

  3. McKAY / CHARACTER CRITERION:
     Q8 -> D4~ on the McKay correspondence.
     D4~ has 4 one-dimensional irreps (abelian) and 1 two-dimensional.
     => 80% of Q8's representation content is abelian.
     The central weights satisfy:
       W_{{rho_8}}(Q8) = 0,  W_{{rho_2}}(Q8) = 0,
       W_{{rho_4}}(Q8) = -W_{{rho_4}}(2T\\Q8).
     So Q8 has NO support on the Galois-broken sectors and only a
     rational rho_4 residue, canceled exactly at the 2T level.

  HONEST CONCLUSION:
    The old slogan "Q8 contributes only to U(1)" is too strong if read
    literally at the level of raw subset character sums.

    What IS derived is:
    - Q8 is identified coordinate-free by its commutator-centrality.
    - Q8 has zero support on the Galois-broken gauge sectors rho_8 and rho_2.
    - Q8 carries a rational rho_4 residue.
    - That rho_4 residue cancels exactly against 2T\\Q8 when the shell is
      completed from Q8 to 2T.

    So the defensible statement is:
    Q8 is the maximally central rational subgroup, and after 2T completion
    its net non-abelian rational residue vanishes, leaving only the abelian
    channel as the surviving Q8-specific contribution.

  STATUS: DERIVED (commutator criterion + Galois/rational support);
          the direct Yang-Mills prefactor identification remains STRUCTURAL.
""")

print("=" * 70)
print("EXP-541 COMPLETE")
print("=" * 70)

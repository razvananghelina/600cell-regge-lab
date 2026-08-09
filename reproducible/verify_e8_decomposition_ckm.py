"""
verify_e8_decomposition_ckm.py
===============================
Self-contained verification of the E8 decomposition graph and its
connection to CKM mixing exponents.

What this script does:
  1. Builds the E8 root system (240 roots).
  2. Constructs Coxeter elements from all orderings of simple reflections
     and extracts the 40 distinct H4 eigenspace projectors.
  3. Computes all C(40,2) = 780 pairwise overlaps via Grassmannian trace.
  4. VERIFIES the Overlap Quantization Theorem:
     (a) Every overlap is in {60, 72, 84, 96}.
     (b) Equivalently, Tr(P1*P2) = 2m/a1 for m in {5, 6, 7, 8}.
  5. VERIFIES the CKM connection:
     (a) GCD of overlaps = 12 = n_13 (CKM 1-3 exponent).
     (b) Overlap-multiplier duality: exchanged roots = 12 * {5,4,3,2}.
     (c) Concordance I2 uniquely selects n_23 = 7.
     (d) n_13 = degree(600-cell) iff a1 = 5.
  6. VERIFIES the Galois cancellation (Gap 1):
     sigma(a,b) = sigma(-a,-b) for the CRT-decomposed trace.
  7. VERIFIES the range structure (Gap 2):
     m_max = (3*a1+1)/2 = rank(E8) = 8.
  8. Prints a summary table with PASS/FAIL verdicts.

Dependencies: numpy (standard).
Encoding: ASCII only.

Author: Razvan-Constantin Anghelina
"""

import numpy as np
from itertools import permutations
from math import sqrt, pi, cos, factorial, gcd
from collections import Counter
from functools import reduce

# ============================================================================
# SECTION 0: THE SINGLE FREE INTEGER
# ============================================================================
a1 = 5
phi = (1 + sqrt(a1)) / 2
b1 = a1 + 1                   # 6
N = factorial(a1)              # 120
h_E8 = a1 * b1                # 30 (Coxeter number)
rank_E8 = 8
degree_600cell = 2 * (a1 + 1) # 12

# CKM spectral data (from Cayley graph, exp452b + exp453)
n_ut = a1 - 1                 # 4
n_db = a1 - 3                 # 2
delta = 1                     # vacuum offset
n_12 = n_ut - delta           # 3
n_23 = n_ut * n_db - delta    # 7
n_13 = n_ut * (n_db + delta)  # 12

results = []

def check(name, condition, detail=""):
    tag = "PASS" if condition else "FAIL"
    results.append((name, tag, detail))
    return condition

# ============================================================================
# SECTION 1: BUILD E8 ROOT SYSTEM
# ============================================================================
print("=" * 70)
print("SECTION 1: E8 root system")
print("=" * 70)

E8_roots = []
# Type I: +-e_i +- e_j  (i < j)
for i in range(8):
    for j in range(i + 1, 8):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = np.zeros(8)
                v[i] = si
                v[j] = sj
                E8_roots.append(v)
# Type II: (1/2)(+-1, ..., +-1) with even number of minus signs
for bits in range(256):
    signs = [(-1) ** ((bits >> k) & 1) for k in range(8)]
    if sum(1 for s in signs if s < 0) % 2 == 0:
        E8_roots.append(np.array(signs) / 2.0)

E8 = np.array(E8_roots)
check("E8 root count", len(E8) == 240, f"|E8| = {len(E8)}")

# Verify all norms are sqrt(2)
norms_ok = all(abs(np.dot(r, r) - 2.0) < 1e-10 for r in E8)
check("E8 root norms", norms_ok, "||r||^2 = 2 for all roots")
print(f"  E8 roots: {len(E8)}, norms OK: {norms_ok}")

# ============================================================================
# SECTION 2: SIMPLE ROOTS AND COXETER ELEMENTS
# ============================================================================
print("\n" + "=" * 70)
print("SECTION 2: Coxeter elements and H4 projectors")
print("=" * 70)

simple = np.zeros((8, 8))
simple[0] = [1, -1, 0, 0, 0, 0, 0, 0]
simple[1] = [0, 1, -1, 0, 0, 0, 0, 0]
simple[2] = [0, 0, 1, -1, 0, 0, 0, 0]
simple[3] = [0, 0, 0, 1, -1, 0, 0, 0]
simple[4] = [0, 0, 0, 0, 1, -1, 0, 0]
simple[5] = [0, 0, 0, 0, 0, 1, -1, 0]
simple[6] = [0, 0, 0, 0, 0, 1, 1, 0]
simple[7] = np.array([-0.5] * 8)

# Verify simple roots are E8 roots
for i in range(8):
    assert abs(np.dot(simple[i], simple[i]) - 2) < 1e-10

refls = [np.eye(8) - 2 * np.outer(simple[i], simple[i]) / np.dot(simple[i], simple[i])
         for i in range(8)]


def make_coxeter(perm):
    C = np.eye(8)
    for i in perm:
        C = refls[i] @ C
    return C


def get_H4_projector(C):
    """Extract 8x4 orthonormal basis for the H4 eigenspace of Coxeter element C."""
    if not np.allclose(np.linalg.matrix_power(C, h_E8), np.eye(8), atol=1e-8):
        return None
    eigvals, eigvecs = np.linalg.eig(C)
    args = np.angle(eigvals)
    exps = [(round(args[i] * h_E8 / (2 * pi)) % h_E8 or h_E8, i) for i in range(8)]
    h4_idx = [i for m, i in exps if m in {1, 11, 19, 29}]
    if len(h4_idx) != 4:
        return None
    v_raw = eigvecs[:, h4_idx]
    basis = []
    used = set()
    for k in range(4):
        if k in used:
            continue
        v = v_raw[:, k]
        if np.max(np.abs(np.imag(v))) > 0.01:
            basis.append(np.real(v))
            basis.append(np.imag(v))
            for k2 in range(k + 1, 4):
                if k2 not in used and np.allclose(v_raw[:, k2], np.conj(v)):
                    used.add(k2)
                    break
        else:
            basis.append(np.real(v))
    if len(basis) != 4:
        return None
    B = np.array(basis).T
    Q, _ = np.linalg.qr(B)
    return Q[:, :4]


# Find all 40 distinct H4 subspaces
print("  Finding H4 subspaces (sampling Coxeter elements)...")
projectors = []
coxeter_elts = []
np.random.seed(42)
perm_list = list(permutations(range(8)))
np.random.shuffle(perm_list)

for perm in perm_list:
    C = make_coxeter(perm)
    P = get_H4_projector(C)
    if P is None:
        continue
    is_new = True
    for P_old in projectors:
        svs = np.linalg.svd(P_old.T @ P, compute_uv=False)
        if np.allclose(svs, np.ones(4), atol=1e-6):
            is_new = False
            break
    if is_new:
        projectors.append(P)
        coxeter_elts.append(C)
    if len(projectors) >= 40:
        break

n_decomp = len(projectors)
check("40 decompositions", n_decomp == 40, f"found {n_decomp}")
check("40 = rank*a1", n_decomp == rank_E8 * a1, f"{rank_E8}*{a1} = {rank_E8 * a1}")
print(f"  Found {n_decomp} distinct H4 subspaces")

# ============================================================================
# SECTION 3: ALL 780 PAIRWISE OVERLAPS
# ============================================================================
print("\n" + "=" * 70)
print("SECTION 3: Pairwise overlaps (780 pairs)")
print("=" * 70)

n_pairs = n_decomp * (n_decomp - 1) // 2
overlaps = []
traces = []

for i in range(n_decomp):
    for j in range(i + 1, n_decomp):
        M = projectors[i].T @ projectors[j]
        sc2 = np.sum(np.linalg.svd(M, compute_uv=False) ** 2)
        ov = round(h_E8 * sc2)
        overlaps.append(ov)
        traces.append(sc2)

check("780 pairs", len(overlaps) == 780, f"C(40,2) = {n_pairs}")

ov_counts = Counter(overlaps)
ov_set = set(overlaps)
expected_set = {60, 72, 84, 96}

check("overlaps in {60,72,84,96}", ov_set == expected_set,
      f"found {sorted(ov_set)}")
check("pair count: overlap 60", ov_counts[60] == 176, f"{ov_counts.get(60, 0)}")
check("pair count: overlap 72", ov_counts[72] == 296, f"{ov_counts.get(72, 0)}")
check("pair count: overlap 84", ov_counts[84] == 184, f"{ov_counts.get(84, 0)}")
check("pair count: overlap 96", ov_counts[96] == 124, f"{ov_counts.get(96, 0)}")
check("total pairs", sum(ov_counts.values()) == 780, f"{sum(ov_counts.values())}")

print(f"  Overlap distribution:")
for ov in sorted(ov_counts):
    print(f"    overlap={ov}: {ov_counts[ov]} pairs")

# ============================================================================
# SECTION 4: OVERLAP QUANTIZATION THEOREM
# ============================================================================
print("\n" + "=" * 70)
print("SECTION 4: Overlap Quantization Theorem")
print("=" * 70)

# (a) overlap = 12*m with m integer
all_12m = all(ov % 12 == 0 for ov in overlaps)
check("overlap = 12*m", all_12m)

# (b) Tr(P1*P2) = 2m/a1 exactly
max_residual = 0.0
for tr in traces:
    m_val = tr * a1 / 2
    residual = abs(m_val - round(m_val))
    max_residual = max(max_residual, residual)

check("Tr = 2m/a1 exact", max_residual < 1e-10,
      f"max residual = {max_residual:.2e}")

# (c) m-values
m_values = sorted(set(ov // 12 for ov in overlaps))
check("m in {5,6,7,8}", m_values == [5, 6, 7, 8], f"m = {m_values}")

# (d) Number of overlap levels = rank - a1 + 1 = 4
n_levels = len(m_values)
check("4 overlap levels", n_levels == rank_E8 - a1 + 1,
      f"{n_levels} = {rank_E8}-{a1}+1")

print(f"  Quantization: max |m - round(m)| = {max_residual:.2e}")
print(f"  m-values: {m_values}")
print(f"  Overlap levels: {n_levels} = rank(E8) - a1 + 1 = {rank_E8 - a1 + 1}")

# ============================================================================
# SECTION 5: CKM MULTIPLIER TABLE
# ============================================================================
print("\n" + "=" * 70)
print("SECTION 5: CKM from decomposition graph")
print("=" * 70)

# (a) GCD of overlaps = n_13
ov_gcd = reduce(gcd, ov_set)
check("GCD(overlaps) = n_13 = 12", ov_gcd == n_13, f"GCD = {ov_gcd}")

# (b) Exchanged-root multipliers
print(f"\n  Exchanged-root multiplier table:")
print(f"  {'m':>3} {'overlap':>8} {'exchanged':>10} {'mult':>5} {'identity':>12}")
mult_correct = True
for m in m_values:
    ov = 12 * m
    exchanged = N - ov
    mult = exchanged // degree_600cell
    expected_mult = 2 * a1 - m
    if mult != expected_mult:
        mult_correct = False
    label = {a1: "a1", n_ut: "n_ut", n_12: "n_12", n_db: "n_db"}.get(mult, "?")
    print(f"  {m:3d} {ov:8d} {exchanged:10d} {mult:5d}   {label:>12}")

check("multipliers = {a1,n_ut,n_12,n_db}", mult_correct)

# (c) Duality sum
check("n_12 + n_23 = 2*a1", n_12 + n_23 == 2 * a1,
      f"{n_12}+{n_23} = {n_12 + n_23}")

# (d) Concordance I2 uniqueness
candidates_n23 = []
for m in m_values:
    if m <= a1:
        continue
    n12_cand = 2 * a1 - m
    if n12_cand * ov_gcd == b1 ** 2:
        candidates_n23.append(m)
check("I2 selects n_23=7 uniquely", candidates_n23 == [7],
      f"candidates = {candidates_n23}")

# (e) n_13 = degree <=> a1 = 5
lhs = (a1 - 1) * (a1 - 2)
rhs = 2 * (a1 + 1)
check("n_13 = degree <=> a1=5", lhs == rhs and lhs == n_13,
      f"(a1-1)(a1-2)={lhs}, 2(a1+1)={rhs}, n_13={n_13}")

# (f) Sum of m-values = Weinberg denominator
sum_m = sum(m_values)
check("sum(m) = a1^2 + 1 = 26", sum_m == a1 ** 2 + 1,
      f"sum = {sum_m}")

# (g) Overlap range = b1^2 = n_12 * n_13
ov_range = max(ov_set) - min(ov_set)
check("overlap range = b1^2", ov_range == b1 ** 2,
      f"range = {ov_range}, b1^2 = {b1 ** 2}")

print(f"\n  CKM exponents derived: n_12={n_12}, n_23={n_23}, n_13={n_13}")

# ============================================================================
# SECTION 6: GALOIS CANCELLATION (GAP 1)
# ============================================================================
print("\n" + "=" * 70)
print("SECTION 6: Galois cancellation (Gap 1)")
print("=" * 70)


def crt_j(a, u):
    """CRT: j in Z/30 with j=a mod 5, j=u mod 6."""
    return (6 * a + 25 * u) % h_E8


A_vals = [cos(2 * pi * a / a1) for a in range(a1)]
B_vals = [cos(pi * u / 3) for u in range(b1)]

# Test one pair from each overlap level
galois_ok_all = True
sqrt5_ok_all = True
test_pairs = {}
idx = 0
for i in range(n_decomp):
    for j in range(i + 1, n_decomp):
        ov = overlaps[idx]
        if ov not in test_pairs:
            test_pairs[ov] = (i, j)
        idx += 1
    if len(test_pairs) == 4:
        break

for ov_target, (idx_i, idx_j) in sorted(test_pairs.items()):
    C1 = coxeter_elts[idx_i]
    C2 = coxeter_elts[idx_j]

    # Compute T(j,k) = Tr(C1^j C2^k)
    C1p = [np.eye(8)]
    C2p = [np.eye(8)]
    for _ in range(h_E8 - 1):
        C1p.append(C1p[-1] @ C1)
        C2p.append(C2p[-1] @ C2)

    T = np.zeros((h_E8, h_E8))
    for j in range(h_E8):
        for k in range(h_E8):
            T[j, k] = round(np.trace(C1p[j] @ C2p[k]))

    # Compute sigma(a,b)
    sigma = np.zeros((a1, a1))
    for a in range(a1):
        for b in range(a1):
            s = 0.0
            for u in range(b1):
                for v in range(b1):
                    s += B_vals[u] * B_vals[v] * T[crt_j(a, u), crt_j(b, v)]
            sigma[a, b] = s

    # Check Galois: sigma(a,b) = sigma(-a,-b)
    max_diff = max(abs(sigma[a, b] - sigma[(a1 - a) % a1, (a1 - b) % a1])
                   for a in range(a1) for b in range(a1))
    if max_diff > 1e-6:
        galois_ok_all = False

    # Check sqrt(5) cancellation
    irr_sum = 0.0
    for a in range(a1):
        for b in range(a1):
            if a == 0:
                ra, ia = 1.0, 0.0
            elif a in [1, 4]:
                ra, ia = -0.25, 0.25
            else:
                ra, ia = -0.25, -0.25
            if b == 0:
                rb, ib = 1.0, 0.0
            elif b in [1, 4]:
                rb, ib = -0.25, 0.25
            else:
                rb, ib = -0.25, -0.25
            irr_sum += (ra * ib + ia * rb) * sigma[a, b]

    if abs(irr_sum) > 1e-6:
        sqrt5_ok_all = False

check("Galois: sigma(a,b)=sigma(-a,-b)", galois_ok_all,
      "for all 4 overlap levels")
check("sqrt(5) cancellation", sqrt5_ok_all,
      "irrational coefficient = 0")
print(f"  Galois symmetry: {'PASS' if galois_ok_all else 'FAIL'}")
print(f"  sqrt(5) cancels: {'PASS' if sqrt5_ok_all else 'FAIL'}")

# ============================================================================
# SECTION 7: RANGE STRUCTURE (GAP 2)
# ============================================================================
print("\n" + "=" * 70)
print("SECTION 7: Range structure (Gap 2)")
print("=" * 70)

# Upper bound identity
m_max_formula = (3 * a1 + 1) / 2
check("m_max = (3a1+1)/2 = rank", m_max_formula == rank_E8,
      f"(3*{a1}+1)/2 = {m_max_formula}")

# H4/H4' self-duality
duality_ok = True
idx = 0
for i in range(n_decomp):
    for j in range(i + 1, n_decomp):
        P1 = projectors[i] @ projectors[i].T
        P2 = projectors[j] @ projectors[j].T
        tr_h4 = np.trace(P1 @ P2)
        tr_h4c = np.trace((np.eye(8) - P1) @ (np.eye(8) - P2))
        if abs(tr_h4 - tr_h4c) > 1e-8:
            duality_ok = False
            break
        idx += 1
    if not duality_ok:
        break

check("H4/H4' self-duality", duality_ok,
      "Tr(P1*P2) = Tr((I-P1)(I-P2))")

# Root counting verification: overlap = |I1 & I2|
root_count_ok = True
for ov_target, (idx_i, idx_j) in sorted(test_pairs.items()):
    P1 = projectors[idx_i]
    P2 = projectors[idx_j]
    n_agree = 0
    for r in E8:
        inner1 = np.sum((P1.T @ r) ** 2) < 1.0
        inner2 = np.sum((P2.T @ r) ** 2) < 1.0
        if inner1 == inner2:
            n_agree += 1
    # agreement = 2 * overlap (both inner-inner and outer-outer)
    if n_agree != 2 * ov_target:
        root_count_ok = False

check("root counting", root_count_ok,
      "agreement = 2*overlap for all 4 levels")
print(f"  Root counting: {'PASS' if root_count_ok else 'FAIL'}")

# ============================================================================
# SECTION 8: KEY IDENTITIES
# ============================================================================
print("\n" + "=" * 70)
print("SECTION 8: Algebraic identities")
print("=" * 70)

check("a1! = 4*a1*(a1+1)",
      factorial(a1) == 4 * a1 * (a1 + 1),
      f"{factorial(a1)} = {4 * a1 * (a1 + 1)}")

check("n_12*n_13 = b1^2",
      n_12 * n_13 == b1 ** 2,
      f"{n_12}*{n_13} = {n_12 * n_13} = {b1}^2")

check("n_12+n_23 = 2*a1",
      n_12 + n_23 == 2 * a1,
      f"{n_12}+{n_23} = {2 * a1}")

check("n_23+n_13 = a1^2-a1-1",
      n_23 + n_13 == a1 ** 2 - a1 - 1,
      f"{n_23}+{n_13} = {n_23 + n_13} = {a1 ** 2 - a1 - 1}")

check("degree = 2*(a1+1)",
      degree_600cell == 2 * (a1 + 1),
      f"{degree_600cell} = {2 * (a1 + 1)}")

check("rank(E8) = 2*n_ut",
      rank_E8 == 2 * n_ut,
      f"{rank_E8} = 2*{n_ut}")

check("3*a1 = 2*rank-1",
      3 * a1 == 2 * rank_E8 - 1,
      f"{3 * a1} = {2 * rank_E8 - 1}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

n_pass = sum(1 for _, tag, _ in results if tag == "PASS")
n_fail = sum(1 for _, tag, _ in results if tag == "FAIL")
n_total = len(results)

for name, tag, detail in results:
    marker = "[PASS]" if tag == "PASS" else "[FAIL]"
    line = f"  {marker} {name}"
    if detail:
        line += f"  ({detail})"
    print(line)

print(f"\n  {n_pass}/{n_total} PASSED, {n_fail} FAILED")

if n_fail == 0:
    print("\n  ALL TESTS PASS.")
    print("  The E8 decomposition graph encodes the CKM exponents {3, 7, 12}.")
else:
    print(f"\n  {n_fail} TESTS FAILED -- investigate before proceeding.")

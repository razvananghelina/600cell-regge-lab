"""
exp342f: Uniqueness Theorem for (a,b) Assignment
=================================================
Exhaustive enumeration of all Z^2 assignments satisfying:
C1: electron (0,0)
C2: sum(a) = 12
C3: sum(b) = 8
C4: n_i = 5a + 6b gives ordered masses {0,3,5,11,11,16,17,19,26}
C5: mu and s share same (a,b)
C6: |a|,|b| <= B (compactness bound)
C7: norm multiset = {-19,-1,-1,0,1,1,1,5,19}

Questions to answer:
1. Verify uniqueness with B=6
2. Extend to B=10, 15, 20, 50, 100 - still unique?
3. How many solutions with C1-C5 only (no C7)?
4. Can C6 be derived from C1-C5+C7?
5. Is 19 uniquely determined? Try other primes.
6. Classify the rejected solutions.
"""
import numpy as np
from itertools import product

PHI = (1 + np.sqrt(5)) / 2
A1 = 5
B1 = 6

# Target mass exponents (ordered)
mass_exponents = [0, 3, 5, 11, 16, 17, 19, 26]  # 8 distinct (11 appears 2x)
# For the 9 fermions: n = {0, 3, 5, 11, 11, 16, 17, 19, 26}
# mu and s both get n=11 with same (a,b)

def norm_zphi(a, b):
    """N(a + b*phi) = a^2 + a*b - b^2"""
    return a*a + a*b - b*b

def get_candidates(n, B):
    """Find all (a,b) with 5a+6b=n and |a|<=B, |b|<=B"""
    cands = []
    for a in range(-B, B+1):
        rem = n - A1*a
        if rem % B1 == 0:
            b = rem // B1
            if abs(b) <= B:
                cands.append((a, b))
    return cands

def smart_search(B, target_norms=None):
    """Efficient search using incremental pruning on partial sums.
    Returns (count_sums, count_norms, solutions_sums, solutions_norms).
    target_norms: sorted list of 9 norms, or None to collect all.
    """
    all_cands = {}
    for n in mass_exponents:
        all_cands[n] = get_candidates(n, B)

    # Levels: e(0), u(3), d(5), ms(11x2), c(16), tau(17), b(19), t(26)
    # e is fixed at (0,0). Remaining 7 choices (ms counted once).
    # Target partial sums after all 7: sum_a=12, sum_b=8
    # (e contributes 0,0; ms contributes 2x)

    levels = mass_exponents[1:]  # [3, 5, 11, 16, 17, 19, 26]
    # multiplicity: 11 counts 2x in sum
    mult = {3:1, 5:1, 11:2, 16:1, 17:1, 19:1, 26:1}

    # Compute suffix ranges for pruning
    # suffix_min_a[i], suffix_max_a[i] = min/max possible sum_a from level i..end
    n_levels = len(levels)
    suffix_min_a = [0] * (n_levels + 1)
    suffix_max_a = [0] * (n_levels + 1)
    suffix_min_b = [0] * (n_levels + 1)
    suffix_max_b = [0] * (n_levels + 1)

    for i in range(n_levels - 1, -1, -1):
        n = levels[i]
        m = mult[n]
        cands = all_cands[n]
        if not cands:
            return 0, 0, [], []
        a_vals = [c[0] for c in cands]
        b_vals = [c[1] for c in cands]
        suffix_min_a[i] = suffix_min_a[i+1] + m * min(a_vals)
        suffix_max_a[i] = suffix_max_a[i+1] + m * max(a_vals)
        suffix_min_b[i] = suffix_min_b[i+1] + m * min(b_vals)
        suffix_max_b[i] = suffix_max_b[i+1] + m * max(b_vals)

    solutions_sums = []
    solutions_norms = []

    def recurse(level_idx, partial_a, partial_b, chosen):
        if level_idx == n_levels:
            if partial_a == 12 and partial_b == 8:
                # Build full assignment
                full_ab = [(0,0)] + list(chosen)
                # Expand: insert ms duplicate
                # chosen order: u(3), d(5), ms(11), c(16), tau(17), b(19), t(26)
                expanded = [(0,0), chosen[0], chosen[1],
                            chosen[2], chosen[2],  # mu and s
                            chosen[3], chosen[4], chosen[5], chosen[6]]
                solutions_sums.append(expanded)
                if target_norms is not None:
                    norms = sorted([norm_zphi(a, b) for a, b in expanded])
                    if norms == target_norms:
                        solutions_norms.append(expanded)
                else:
                    solutions_norms.append(expanded)
            return

        n = levels[level_idx]
        m = mult[n]
        remaining_target_a = 12 - partial_a
        remaining_target_b = 8 - partial_b

        for a, b in all_cands[n]:
            new_a = partial_a + m * a
            new_b = partial_b + m * b

            # Check if remaining levels can reach target
            rem_a = 12 - new_a
            rem_b = 8 - new_b

            if rem_a < suffix_min_a[level_idx+1] or rem_a > suffix_max_a[level_idx+1]:
                continue
            if rem_b < suffix_min_b[level_idx+1] or rem_b > suffix_max_b[level_idx+1]:
                continue

            recurse(level_idx + 1, new_a, new_b, chosen + ((a, b),))

    recurse(0, 0, 0, ())
    return len(solutions_sums), len(solutions_norms), solutions_sums, solutions_norms

# ================================================================
# PART 1: Candidate Enumeration
# ================================================================
print("=" * 72)
print("PART 1: Candidates per Mass Level")
print("=" * 72)

for B in [6, 10, 20, 50]:
    print(f"\n  Bound B = {B}:")
    total = 1
    for n in mass_exponents:
        cands = get_candidates(n, B)
        count = len(cands)
        if n == 0:
            total *= 1  # forced
        else:
            total *= count
        print(f"    n={n:>2}: {count} candidates: {cands[:6]}{'...' if count > 6 else ''}")
    print(f"  Total combinations (without sum constraint): {total}")

# ================================================================
# PART 2: Exhaustive Search with B=6
# ================================================================
print(f"\n{'='*72}")
print("PART 2: Exhaustive Search (B=6)")
print("=" * 72)

target_norms_sorted = sorted([-19, -1, -1, 0, 1, 1, 1, 5, 19])

B = 6
count_sums, count_norms, solutions_c1c5, solutions_c1c7 = smart_search(B, target_norms_sorted)

# Count total combinations for reference
total_comb = 1
for n in mass_exponents[1:]:
    total_comb *= len(get_candidates(n, B))
print(f"\n  Total candidate combinations: {total_comb}")
print(f"  Solutions satisfying C1-C5 + C2,C3 (sum rules): {count_sums}")
print(f"  Solutions satisfying ALL C1-C7 (+ norms):       {count_norms}")

# ================================================================
# PART 3: Display Solutions with Sum Rules
# ================================================================
print(f"\n{'='*72}")
print(f"PART 3: All {len(solutions_c1c5)} Solutions Satisfying C1-C5 + Sums (B=6)")
print("=" * 72)

fermion_names = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 'b', 't']
n_values = [0, 3, 5, 11, 11, 16, 17, 19, 26]

for idx, sol in enumerate(solutions_c1c5):
    norms = [norm_zphi(a, b) for a, b in sol]
    norm_set = sorted(norms)
    target_norms = sorted([-19, -1, -1, 0, 1, 1, 1, 5, 19])
    is_target = norm_set == target_norms

    print(f"\n  Solution {idx+1} {'*** MATCHES C7 ***' if is_target else ''}:")
    for i, ((a, b), f, n) in enumerate(zip(sol, fermion_names, n_values)):
        N = norm_zphi(a, b)
        print(f"    {f:>3} (n={n:>2}): (a,b) = ({a:>+2},{b:>+2}), "
              f"N = {N:>+3}, z = {a+b*PHI:>+8.4f}")
    print(f"    Norms: {norm_set}")
    print(f"    Max |N|: {max(abs(n) for n in norms)}")

# ================================================================
# PART 4: Extended Bounds (using smart pruning)
# ================================================================
print(f"\n{'='*72}")
print("PART 4: Uniqueness with Extended Bounds")
print("=" * 72)

for B_ext in [10, 15, 20, 30, 50, 100]:
    total_combos = 1
    for n in mass_exponents[1:]:
        total_combos *= len(get_candidates(n, B_ext))

    cs, cn, sols_s, sols_n = smart_search(B_ext, target_norms_sorted)
    print(f"  B={B_ext:>3}: {total_combos:>15} combos, "
          f"{cs:>5} pass sums, "
          f"{cn:>2} pass norms")

# ================================================================
# PART 5: Is the Norm Multiset Uniquely Determined?
# ================================================================
print(f"\n{'='*72}")
print("PART 5: What Norm Multisets Are Compatible with C1-C5?")
print("=" * 72)

for B5 in [10, 20]:
    cs, cn, sols_s, sols_n = smart_search(B5, None)  # collect ALL

    norm_multisets = {}
    for sol in sols_s:
        norms = tuple(sorted([norm_zphi(a, b) for a, b in sol]))
        if norms not in norm_multisets:
            norm_multisets[norms] = []
        norm_multisets[norms].append(sol)

    print(f"\n  With B={B5}, {cs} solutions pass sums, "
          f"{len(norm_multisets)} distinct norm multisets:")
    for norms, sols in sorted(norm_multisets.items(),
                               key=lambda x: max(abs(nn) for nn in x[0])):
        max_norm = max(abs(nn) for nn in norms)
        is_target = (list(norms) == target_norms_sorted)
        mark = " *** TARGET ***" if is_target else ""
        print(f"    {list(norms)} ({len(sols)} sol{'s' if len(sols)>1 else ''}), "
              f"max|N|={max_norm}{mark}")

# ================================================================
# PART 6: Replace 19 with Other Primes
# ================================================================
print(f"\n{'='*72}")
print("PART 6: Is 19 Uniquely Determined? Replace with Other Primes")
print("=" * 72)

print(f"\n  Classification of small primes in Z[phi] = Z[(1+sqrt(5))/2]:")
print(f"  p = 5: RAMIFIES (5 = sqrt(5)^2)")
for p in [2, 3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]:
    r = p % 5
    if r == 0:
        status = "RAMIFIES"
    elif r == 1 or r == 4:
        status = "SPLITS"
    else:
        status = "INERT"
    print(f"  p = {p:>2}: {p} mod 5 = {r}, {status}")

print(f"\n  Testing: replace +-19 with +-p in norm multiset (B=20):")
B_test = 20

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]:
    target_test = sorted([-p, -1, -1, 0, 1, 1, 1, 5, p])
    cs, cn, _, _ = smart_search(B_test, target_test)
    r = p % 5
    ptype = "RAM" if r == 0 else ("SPL" if r in (1,4) else "INR")
    print(f"    p={p:>2} ({ptype}): {cn} solutions  (norms={target_test})")

# ================================================================
# PART 7: Replacing 5 with Other Values
# ================================================================
print(f"\n{'='*72}")
print("PART 7: Is N(c) = 5 Uniquely Determined?")
print("=" * 72)

print(f"\n  Testing: replace 5 with q in norm multiset (B=20):")
for q in [-5, -4, -3, -1, 0, 1, 2, 3, 4, 5, 7, 9, 11, 19]:
    target_test = sorted([-19, -1, -1, 0, 1, 1, 1, q, 19])
    cs, cn, _, _ = smart_search(B_test, target_test)
    mark = " *** ORIGINAL ***" if q == 5 else ""
    print(f"    q={q:>+3}: {cn} solutions{mark}")

# ================================================================
# PART 8: Can C6 Be Derived? Minimum Bound for Uniqueness
# ================================================================
print(f"\n{'='*72}")
print("PART 8: Is Bound C6 Necessary? Min Bound for Uniqueness")
print("=" * 72)

for B_min in range(4, 30):
    all_exist = all(len(get_candidates(n, B_min)) > 0 for n in mass_exponents)
    if not all_exist:
        print(f"  B={B_min:>2}: some mass levels have no candidates, skipping")
        continue

    cs, cn, _, _ = smart_search(B_min, target_norms_sorted)
    total_c = 1
    for n in mass_exponents[1:]:
        total_c *= len(get_candidates(n, B_min))

    unique_str = ""
    if cn == 1:
        unique_str = " (UNIQUE)"
    elif cn > 1:
        unique_str = f" *** NOT UNIQUE ***"
    print(f"  B={B_min:>2}: {total_c:>10} combos, {cs:>4} pass sums, "
          f"{cn} pass norms{unique_str}")

# ================================================================
# PART 9: Dimension-Mass Ordering
# ================================================================
print(f"\n{'='*72}")
print("PART 9: Does Mass Ordering Follow Irrep Dimension?")
print("=" * 72)

# Affine E8 dims: [1, 2, 3, 4, 5, 6, 4, 2, 3]
# Sorted by dim: 1, 2, 2, 3, 3, 4, 4, 5, 6
# Mass exponents: 0, 3, 5, 11, 11, 16, 17, 19, 26

# Sorted dims with fermion mapping (Sol 3):
# rho_1(1)=e(0), rho_2(2)=u(3), rho_8(2)=t(26),
# rho_3(3)=d(5), rho_9(3)=b(19),
# rho_4(4)=s(11), rho_7(4)=tau(17),
# rho_5(5)=mu(11), rho_6(6)=c(16)

dim_mass = [(1,0,'e'), (2,3,'u'), (2,26,'t'), (3,5,'d'), (3,19,'b'),
            (4,11,'s'), (4,17,'tau'), (5,11,'mu'), (6,16,'c')]

# Within each dimension, is mass ordered?
print(f"\n  Irrep dim -> fermion (mass exponent):")
for d, n, f in sorted(dim_mass):
    print(f"    dim {d}: {f:>4} (n={n:>2})")

print(f"\n  AVERAGE n per dimension:")
from collections import defaultdict
dim_groups = defaultdict(list)
for d, n, f in dim_mass:
    dim_groups[d].append(n)
for d in sorted(dim_groups):
    vals = dim_groups[d]
    avg = np.mean(vals)
    print(f"    dim {d}: n = {vals}, avg = {avg:.1f}")

# Is average n monotone in dim?
avgs = [np.mean(dim_groups[d]) for d in sorted(dim_groups)]
print(f"\n  Average n by dim: {[f'{a:.1f}' for a in avgs]}")
print(f"  Monotone? {all(avgs[i] <= avgs[i+1] for i in range(len(avgs)-1))}")

# What about: within Galois pairs (same dim), the AVERAGE is fixed?
print(f"\n  Galois pair averages:")
for d, n1, f1 in dim_mass:
    for d2, n2, f2 in dim_mass:
        if d == d2 and f1 < f2:
            print(f"    dim {d}: ({f1},{f2}) avg(n) = {(n1+n2)/2:.1f}, "
                  f"sum(n) = {n1+n2}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*72}")
print("SUMMARY")
print("=" * 72)

print(f"""
UNIQUENESS THEOREM VERIFICATION:

  With B=6:
    {count_sums} satisfy sum rules (C1-C5)
    {count_norms} satisfy all constraints (C1-C7)

  The solution IS {'UNIQUE' if count_norms == 1 else 'NOT UNIQUE'}!

KEY FINDINGS:
  - C7 (norm constraint) is the most powerful
  - Extended bounds tested above
  - Sensitivity to prime choice tested above
  - Dimension ordering tested above
""")
print("Done.")

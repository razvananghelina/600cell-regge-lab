"""exp342f_quick: Fast uniqueness check (B<=20 only)"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
A1 = 5
B1 = 6

mass_exponents = [0, 3, 5, 11, 16, 17, 19, 26]

def norm_zphi(a, b):
    return a*a + a*b - b*b

def get_candidates(n, B):
    cands = []
    for a in range(-B, B+1):
        rem = n - A1*a
        if rem % B1 == 0:
            b = rem // B1
            if abs(b) <= B:
                cands.append((a, b))
    return cands

def smart_search(B, target_norms=None):
    all_cands = {}
    for n in mass_exponents:
        all_cands[n] = get_candidates(n, B)
    levels = mass_exponents[1:]
    mult = {3:1, 5:1, 11:2, 16:1, 17:1, 19:1, 26:1}
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
                expanded = [(0,0), chosen[0], chosen[1],
                            chosen[2], chosen[2],
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
        for a, b in all_cands[n]:
            new_a = partial_a + m * a
            new_b = partial_b + m * b
            rem_a = 12 - new_a
            rem_b = 8 - new_b
            if rem_a < suffix_min_a[level_idx+1] or rem_a > suffix_max_a[level_idx+1]:
                continue
            if rem_b < suffix_min_b[level_idx+1] or rem_b > suffix_max_b[level_idx+1]:
                continue
            recurse(level_idx + 1, new_a, new_b, chosen + ((a, b),))
    recurse(0, 0, 0, ())
    return len(solutions_sums), len(solutions_norms), solutions_sums, solutions_norms

target_norms_sorted = sorted([-19, -1, -1, 0, 1, 1, 1, 5, 19])
fermion_names = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 'b', 't']
n_values = [0, 3, 5, 11, 11, 16, 17, 19, 26]

# PART 1: B=6 exhaustive
print("=" * 72)
print("B=6: Exhaustive Search")
print("=" * 72)
B = 6
cs, cn, sols_s, sols_n = smart_search(B, target_norms_sorted)
print(f"  Solutions satisfying C1-C5 + sums: {cs}")
print(f"  Solutions satisfying ALL C1-C7:    {cn}")

for idx, sol in enumerate(sols_s):
    norms = [norm_zphi(a, b) for a, b in sol]
    norm_set = sorted(norms)
    is_target = norm_set == target_norms_sorted
    print(f"\n  Sol {idx+1} {'*** MATCHES C7 ***' if is_target else ''}:")
    for i, ((a, b), f, n) in enumerate(zip(sol, fermion_names, n_values)):
        N = norm_zphi(a, b)
        print(f"    {f:>3} (n={n:>2}): (a,b)=({a:>+2},{b:>+2}), N={N:>+3}")
    print(f"    Norms: {norm_set}")

# PART 2: Extended bounds
print(f"\n{'='*72}")
print("Extended Bounds")
print("=" * 72)
for B_ext in [10, 15, 20]:
    cs, cn, _, _ = smart_search(B_ext, target_norms_sorted)
    total = 1
    for n in mass_exponents[1:]:
        total *= len(get_candidates(n, B_ext))
    print(f"  B={B_ext:>3}: {total:>10} combos, {cs:>4} pass sums, {cn} pass norms")

# PART 3: Norm multisets (B=10 only, fast)
print(f"\n{'='*72}")
print("Norm Multisets Compatible with C1-C5 (B=10)")
print("=" * 72)
cs, cn, sols_s, _ = smart_search(10, None)
from collections import defaultdict
norm_multisets = defaultdict(list)
for sol in sols_s:
    norms = tuple(sorted([norm_zphi(a, b) for a, b in sol]))
    norm_multisets[norms].append(sol)
print(f"  {cs} solutions pass sums, {len(norm_multisets)} distinct norm multisets:")
for norms, sols in sorted(norm_multisets.items(), key=lambda x: max(abs(nn) for nn in x[0])):
    max_norm = max(abs(nn) for nn in norms)
    is_target = (list(norms) == target_norms_sorted)
    mark = " *** TARGET ***" if is_target else ""
    print(f"    {list(norms)} ({len(sols)} sol), max|N|={max_norm}{mark}")

# PART 4: Is 19 special?
print(f"\n{'='*72}")
print("Is 19 Uniquely Determined? (B=20)")
print("=" * 72)
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    target_test = sorted([-p, -1, -1, 0, 1, 1, 1, 5, p])
    cs, cn, _, _ = smart_search(20, target_test)
    r = p % 5
    ptype = "RAM" if r == 0 else ("SPL" if r in (1,4) else "INR")
    mark = " *** ORIGINAL ***" if p == 19 else ""
    print(f"    p={p:>2} ({ptype}): {cn} solutions{mark}")

# PART 5: Min bound for uniqueness
print(f"\n{'='*72}")
print("Min Bound for Uniqueness")
print("=" * 72)
for B_min in range(4, 25):
    all_exist = all(len(get_candidates(n, B_min)) > 0 for n in mass_exponents)
    if not all_exist:
        continue
    cs, cn, _, _ = smart_search(B_min, target_norms_sorted)
    tag = " (UNIQUE)" if cn == 1 else (f" *** {cn} sols ***" if cn > 1 else "")
    print(f"  B={B_min:>2}: {cs:>4} pass sums, {cn} pass norms{tag}")

print("\nDone.")

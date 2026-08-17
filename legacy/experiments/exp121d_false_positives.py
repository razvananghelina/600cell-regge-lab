"""
EXP-121d: Why n=1, n=6, n=23 are FALSE POSITIVES
===================================================

From exp121c: The norm + 3-line rule selects 11 levels from 0-30.
But only 8 are occupied. Three are FALSE POSITIVES:

  n=1:  (-1,1), on b=1 line, N=-1, |N|=1
  n=6:  (0,1),  on b=1 line, N=-1, |N|=1
  n=23: (1,3),  on a=1 line, N=-5, |N|=5

This experiment systematically investigates WHAT distinguishes
these 3 from the 8 occupied levels.

Tests:
 1. Phi-power analysis for |N|=1 units
 2. Position on line (a-values on b=1, b-values on a=1)
 3. a >= 0 AND b >= 0 test
 4. Signed norm pattern
 5. Complexity analysis
 6. Val * conjugate analysis (|val - conj|)
 7. a > 0 test
 8. max(|a|, |b|) threshold
 9. |a|+|b| modular
10. a*b product sign/value
11. Fibonacci / Lucas number test
12. n*phi mod 1 (equidistribution)
13. Powers of 2 on b=1 line
14. a as Fibonacci on b=1 line
15. Geometric progression on b=1 line

Author: Claude (exp121d)
Date: 2026-02-08
"""

import numpy as np
from math import gcd

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2

# ============================================================
# SETUP
# ============================================================

def min_complexity(n):
    """Find (a,b) with minimum |a|+|b| such that 5a+6b=n."""
    best = None
    best_c = float('inf')
    for a in range(-30, 31):
        if (n - 5*a) % 6 == 0:
            b = (n - 5*a) // 6
            c = abs(a) + abs(b)
            if c < best_c:
                best_c = c
                best = (a, b)
    return best, best_c

def zphi_norm_signed(a, b):
    """Signed norm N(a + b*phi) = a^2 + ab - b^2."""
    return a**2 + a*b - b**2

def zphi_norm(a, b):
    """Absolute norm |N(a + b*phi)|."""
    return abs(zphi_norm_signed(a, b))

# Data
occupied_n = {0, 3, 5, 11, 16, 17, 19, 26}
allowed_norms = {0, 1, 5, 19}

# Build all levels 0-30 with their properties
levels = {}
for n in range(31):
    (a, b), c = min_complexity(n)
    N_signed = zphi_norm_signed(a, b)
    N_abs = abs(N_signed)
    val = a + b * PHI
    conj = a + b * PHI_CONJ
    on_a1 = (a == 1)
    on_b1 = (b == 1)
    on_325 = (3*a + 2*b == 5)
    is_origin = (a == 0 and b == 0)
    on_line = on_a1 or on_b1 or on_325 or is_origin
    norm_ok = N_abs in allowed_norms
    passes_both = norm_ok and on_line
    is_occ = n in occupied_n

    levels[n] = {
        'a': a, 'b': b, 'c': c, 'n': n,
        'val': val, 'conj': conj,
        'N_signed': N_signed, 'N_abs': N_abs,
        'on_a1': on_a1, 'on_b1': on_b1, 'on_325': on_325,
        'is_origin': is_origin, 'on_line': on_line,
        'norm_ok': norm_ok, 'passes_both': passes_both,
        'is_occ': is_occ
    }

# Identify the 11 that pass both tests
passes = [n for n in range(31) if levels[n]['passes_both']]
false_pos = [n for n in passes if not levels[n]['is_occ']]
true_pos = [n for n in passes if levels[n]['is_occ']]

print("=" * 80)
print("EXP-121d: Why n=1, n=6, n=23 are FALSE POSITIVES")
print("=" * 80)
print()

print("  Pass norm + 3-line test (11 levels):", passes)
print("  Occupied (8 levels, TRUE positives):", true_pos)
print("  FALSE positives (3 levels):", false_pos)
print()

# Print full table
print("  Table of all 11 levels passing both tests:")
print(f"  {'n':>3s} {'(a,b)':>7s} {'C':>2s} {'val':>8s} {'conj':>8s}"
      f" {'N':>4s} {'|N|':>4s} {'line':>10s} {'status':>8s}")
print(f"  {'-'*62}")

for n in passes:
    d = levels[n]
    line_name = ""
    if d['is_origin']: line_name = "origin"
    elif d['on_a1']: line_name = "a=1"
    elif d['on_b1']: line_name = "b=1"
    elif d['on_325']: line_name = "3a+2b=5"

    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {d['c']:>2d} {d['val']:>8.4f}"
          f" {d['conj']:>8.4f} {d['N_signed']:>4d} {d['N_abs']:>4d}"
          f" {line_name:>10s} {status:>8s}")

print()


# ============================================================
# TEST 1: Phi-power analysis for |N|=1
# ============================================================

print("=" * 80)
print("TEST 1: Phi-power analysis for |N|=1 (units in Z[phi])")
print("=" * 80)
print()

# In Z[phi], elements with |N|=1 are units: a+b*phi = +-phi^k
# The fundamental unit is phi, and phi^k gives all units.
# phi^0 = 1, phi^1 = phi, phi^(-1) = phi-1 = 1/phi = 0.618...
# phi^2 = phi+1 = 2.618..., phi^(-2) = 2-phi = 0.382...

# For negative: -phi^0 = -1, -phi^1 = -phi, etc.

print("  Reference: phi^k for k = -4..6:")
for k in range(-4, 7):
    pk = PHI**k
    # Express as a + b*phi using Fibonacci recurrence
    # phi^k = F(k-1) + F(k)*phi where F is Fibonacci (extended to negative)
    print(f"    phi^{k:>3d} = {pk:>10.6f}")
print()

# Find phi-power for each |N|=1 level
print("  Phi-power identification for |N|=1 levels:")
print(f"  {'n':>3s} {'(a,b)':>7s} {'val':>10s} {'sign':>5s} {'k':>4s}"
      f" {'check':>10s} {'err':>10s} {'status':>8s}")
print(f"  {'-'*65}")

unit_levels = [n for n in passes if levels[n]['N_abs'] == 1]

for n in unit_levels:
    d = levels[n]
    val = d['val']
    status = "OCC" if d['is_occ'] else "FALSE+"

    if abs(val) < 1e-12:
        print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {val:>10.6f}"
              f" {'':>5s} {'N/A':>4s} {'N/A':>10s} {'N/A':>10s} {status:>8s}")
        continue

    # Try val = +phi^k and val = -phi^k
    best_k = None
    best_sign = None
    best_err = float('inf')

    for sign in [1, -1]:
        for k in range(-10, 15):
            check = sign * PHI**k
            err = abs(val - check)
            if err < best_err:
                best_err = err
                best_k = k
                best_sign = sign

    sign_str = "+" if best_sign > 0 else "-"
    check_val = best_sign * PHI**best_k

    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {val:>10.6f}"
          f" {sign_str:>5s} {best_k:>4d} {check_val:>10.6f}"
          f" {best_err:>10.2e} {status:>8s}")

print()

# Extract the k values for occupied vs false positive
print("  Summary of phi-power exponents k (where val = sign * phi^k):")
occ_k_values = []
fp_k_values = []

for n in unit_levels:
    d = levels[n]
    val = d['val']
    if abs(val) < 1e-12:
        continue
    for sign in [1, -1]:
        for k in range(-10, 15):
            check = sign * PHI**k
            if abs(val - check) < 1e-6:
                if d['is_occ']:
                    occ_k_values.append((n, k, sign))
                else:
                    fp_k_values.append((n, k, sign))
                break
        else:
            continue
        break

print(f"  Occupied |N|=1: {[(n, f'{'+' if s>0 else '-'}phi^{k}') for n,k,s in occ_k_values]}")
print(f"  False+   |N|=1: {[(n, f'{'+' if s>0 else '-'}phi^{k}') for n,k,s in fp_k_values]}")
print()

# Check: is there a pattern in k?
occ_k_only = [k for _,k,_ in occ_k_values]
fp_k_only = [k for _,k,_ in fp_k_values]
occ_signs = [s for _,_,s in occ_k_values]
fp_signs = [s for _,_,s in fp_k_values]

print(f"  Occupied k-values: {occ_k_only}")
print(f"  False+   k-values: {fp_k_only}")
print(f"  Occupied signs: {occ_signs}")
print(f"  False+   signs: {fp_signs}")
print()

# Check if occupied have k >= some threshold
if occ_k_only and fp_k_only:
    min_occ_k = min(occ_k_only)
    max_fp_k = max(fp_k_only)
    print(f"  Min occupied k = {min_occ_k}, Max false+ k = {max_fp_k}")
    if min_occ_k > max_fp_k:
        print(f"  --> k >= {min_occ_k} separates! (PERFECT)")
    else:
        print(f"  --> k-value alone does NOT separate")

print()


# ============================================================
# TEST 2: Position on line
# ============================================================

print("=" * 80)
print("TEST 2: Position on each line")
print("=" * 80)
print()

# b=1 line: points with b=1
b1_levels = [n for n in passes if levels[n]['on_b1']]
print("  Levels on b=1 line:")
print(f"  {'n':>3s} {'(a,b)':>7s} {'a-value':>8s} {'status':>8s}")
for n in sorted(b1_levels):
    d = levels[n]
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {d['a']:>8d} {status:>8s}")

occ_a_on_b1 = sorted([levels[n]['a'] for n in b1_levels if levels[n]['is_occ']])
fp_a_on_b1 = sorted([levels[n]['a'] for n in b1_levels if not levels[n]['is_occ']])
print(f"\n  Occupied a-values on b=1: {occ_a_on_b1}")
print(f"  False+   a-values on b=1: {fp_a_on_b1}")
print(f"  Pattern: occupied have a >= 1, false+ have a <= 0!")
if all(a >= 1 for a in occ_a_on_b1) and all(a <= 0 for a in fp_a_on_b1):
    print(f"  --> PERFECT separation on b=1 line: a >= 1")
print()

# a=1 line: points with a=1
a1_levels = [n for n in passes if levels[n]['on_a1']]
print("  Levels on a=1 line:")
print(f"  {'n':>3s} {'(a,b)':>7s} {'b-value':>8s} {'status':>8s}")
for n in sorted(a1_levels):
    d = levels[n]
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {d['b']:>8d} {status:>8s}")

occ_b_on_a1 = sorted([levels[n]['b'] for n in a1_levels if levels[n]['is_occ']])
fp_b_on_a1 = sorted([levels[n]['b'] for n in a1_levels if not levels[n]['is_occ']])
print(f"\n  Occupied b-values on a=1: {occ_b_on_a1}")
print(f"  False+   b-values on a=1: {fp_b_on_a1}")
print(f"  Pattern: occupied have b <= 2, false+ has b = 3")
if all(b <= 2 for b in occ_b_on_a1) and all(b >= 3 for b in fp_b_on_a1):
    print(f"  --> PERFECT separation on a=1 line: b <= 2")
print()

# 3a+2b=5 line: all occupied, no false positives
line3_levels = [n for n in passes if levels[n]['on_325']]
print("  Levels on 3a+2b=5 line:")
for n in sorted(line3_levels):
    d = levels[n]
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {status:>8s}")
print(f"  --> No false positives on 3a+2b=5 line")
print()


# ============================================================
# TEST 3: a >= 0 AND b >= 0 test
# ============================================================

print("=" * 80)
print("TEST 3: Non-negativity test (a >= 0 AND b >= 0)")
print("=" * 80)
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'a>=0':>5s} {'b>=0':>5s} {'both':>5s} {'status':>8s}")
print(f"  {'-'*40}")

for n in passes:
    d = levels[n]
    a, b = d['a'], d['b']
    a_ok = a >= 0
    b_ok = b >= 0
    both = a_ok and b_ok
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({a:>2d},{b:>2d}) {str(a_ok):>5s} {str(b_ok):>5s}"
          f" {str(both):>5s} {status:>8s}")

print()

# Count
occ_both = sum(1 for n in true_pos if levels[n]['a'] >= 0 and levels[n]['b'] >= 0)
fp_both = sum(1 for n in false_pos if levels[n]['a'] >= 0 and levels[n]['b'] >= 0)
print(f"  Occupied passing a>=0 AND b>=0: {occ_both}/8")
print(f"  False+ passing a>=0 AND b>=0: {fp_both}/3")
print()

# This doesn't work because:
# - u(3,-2) has b<0 but is occupied
# - b(-1,4) has a<0 but is occupied
# - n=6 (0,1) has both >= 0 but is false+
print("  VERDICT: a>=0 AND b>=0 does NOT work.")
print("  Fails for u(3,-2) [b<0, occupied] and b(-1,4) [a<0, occupied]")
print("  Also n=6 (0,1) passes but is false positive")
print()


# ============================================================
# TEST 4: Signed norm pattern
# ============================================================

print("=" * 80)
print("TEST 4: Signed norm N (not |N|)")
print("=" * 80)
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'N':>4s} {'|N|':>4s} {'N sign':>7s} {'status':>8s}")
print(f"  {'-'*40}")

for n in passes:
    d = levels[n]
    N = d['N_signed']
    sign_str = "+" if N > 0 else ("-" if N < 0 else "0")
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {N:>4d} {d['N_abs']:>4d}"
          f" {sign_str:>7s} {status:>8s}")

print()

# Check signed norm distribution
occ_N = [levels[n]['N_signed'] for n in true_pos]
fp_N = [levels[n]['N_signed'] for n in false_pos]
print(f"  Occupied signed norms: {occ_N}")
print(f"  False+   signed norms: {fp_N}")
print()

occ_pos = sum(1 for N in occ_N if N > 0)
occ_neg = sum(1 for N in occ_N if N < 0)
occ_zero = sum(1 for N in occ_N if N == 0)
fp_pos = sum(1 for N in fp_N if N > 0)
fp_neg = sum(1 for N in fp_N if N < 0)
print(f"  Occupied: {occ_pos} positive, {occ_neg} negative, {occ_zero} zero")
print(f"  False+:   {fp_pos} positive, {fp_neg} negative")
print()

# All false positives have N < 0!
if all(N < 0 for N in fp_N):
    print("  ALL false positives have N < 0!")
    # But do any occupied also have N < 0?
    occ_neg_levels = [n for n in true_pos if levels[n]['N_signed'] < 0]
    if occ_neg_levels:
        print(f"  But occupied levels with N < 0: {occ_neg_levels}")
        for n in occ_neg_levels:
            print(f"    n={n}: ({levels[n]['a']},{levels[n]['b']}), N={levels[n]['N_signed']}")
        print("  --> Signed norm alone does NOT separate")
    else:
        print("  AND no occupied level has N < 0!")
        print("  --> N >= 0 PERFECTLY separates! (DERIVED)")
else:
    print("  Signed norm pattern: mixed, no clean separation")

print()


# ============================================================
# TEST 5: Complexity analysis
# ============================================================

print("=" * 80)
print("TEST 5: Complexity C = |a| + |b|")
print("=" * 80)
print()

occ_C = sorted([levels[n]['c'] for n in true_pos])
fp_C = sorted([levels[n]['c'] for n in false_pos])
print(f"  Occupied complexities: {occ_C}")
print(f"  False+   complexities: {fp_C}")
print(f"  Occupied set: {sorted(set(occ_C))}")
print(f"  False+   set: {sorted(set(fp_C))}")
print()

# Check overlap
overlap = set(occ_C) & set(fp_C)
print(f"  Overlap: {sorted(overlap) if overlap else 'none'}")
if overlap:
    print("  --> Complexity alone does NOT separate (overlap exists)")
else:
    print("  --> Complexity PERFECTLY separates!")
print()

# Also check C modular
print("  C mod 2:")
for n in passes:
    d = levels[n]
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"    n={n:>2d}: C={d['c']}, C%2={d['c']%2} {status}")

print()


# ============================================================
# TEST 6: |val - conj| analysis
# ============================================================

print("=" * 80)
print("TEST 6: |val - conj| = |b| * sqrt(5)")
print("=" * 80)
print()

# val - conj = b*(phi - phi_conj) = b*sqrt(5)
# So |val - conj| = |b|*sqrt(5)

print(f"  {'n':>3s} {'(a,b)':>7s} {'val':>8s} {'conj':>8s} {'|v-c|':>8s}"
      f" {'|b|*s5':>8s} {'status':>8s}")
print(f"  {'-'*55}")

for n in passes:
    d = levels[n]
    diff = abs(d['val'] - d['conj'])
    b_s5 = abs(d['b']) * np.sqrt(5)
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {d['val']:>8.4f} {d['conj']:>8.4f}"
          f" {diff:>8.4f} {b_s5:>8.4f} {status:>8s}")

print()
print("  Note: |val - conj| = |b| * sqrt(5) exactly (algebraic identity)")
print("  So this is equivalent to testing |b|")
print()


# ============================================================
# TEST 7: a > 0 test
# ============================================================

print("=" * 80)
print("TEST 7: a > 0 test")
print("=" * 80)
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'a>0':>5s} {'status':>8s}")
print(f"  {'-'*30}")

for n in passes:
    d = levels[n]
    a_pos = d['a'] > 0
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {str(a_pos):>5s} {status:>8s}")

print()

occ_pass = sum(1 for n in true_pos if levels[n]['a'] > 0)
fp_pass = sum(1 for n in false_pos if levels[n]['a'] > 0)
print(f"  Occupied passing a>0: {occ_pass}/8")
print(f"  False+   passing a>0: {fp_pass}/3")
print()

# Electron (0,0) has a=0, and b(-1,4) has a=-1
# Both are occupied but fail a>0
occ_fail_a = [n for n in true_pos if levels[n]['a'] <= 0]
fp_pass_a = [n for n in false_pos if levels[n]['a'] > 0]
if occ_fail_a:
    print(f"  Occupied failing a>0: {occ_fail_a}")
    for n in occ_fail_a:
        print(f"    n={n}: ({levels[n]['a']},{levels[n]['b']})")
if fp_pass_a:
    print(f"  False+ passing a>0: {fp_pass_a}")
    for n in fp_pass_a:
        print(f"    n={n}: ({levels[n]['a']},{levels[n]['b']})")
print()
print("  VERDICT: a>0 eliminates n=1 (-1,1) and n=6 (0,1)")
print("  but NOT n=23 (1,3). Also fails for e(0,0) and b(-1,4).")
print()


# ============================================================
# TEST 8: max(|a|, |b|) threshold
# ============================================================

print("=" * 80)
print("TEST 8: max(|a|, |b|) threshold")
print("=" * 80)
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'max':>4s} {'status':>8s}")
print(f"  {'-'*28}")

for n in passes:
    d = levels[n]
    m = max(abs(d['a']), abs(d['b']))
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {m:>4d} {status:>8s}")

print()
occ_max = sorted([max(abs(levels[n]['a']), abs(levels[n]['b'])) for n in true_pos])
fp_max = sorted([max(abs(levels[n]['a']), abs(levels[n]['b'])) for n in false_pos])
print(f"  Occupied max(|a|,|b|): {occ_max}")
print(f"  False+   max(|a|,|b|): {fp_max}")
print(f"  VERDICT: Does not separate (all in range 0-4)")
print()


# ============================================================
# TEST 9: |a|+|b| (complexity) modular
# ============================================================

print("=" * 80)
print("TEST 9: Complexity C modular patterns")
print("=" * 80)
print()

for mod in [2, 3, 5]:
    occ_res = sorted(set(levels[n]['c'] % mod for n in true_pos))
    fp_res = sorted(set(levels[n]['c'] % mod for n in false_pos))
    only_occ = sorted(set(occ_res) - set(fp_res))
    only_fp = sorted(set(fp_res) - set(occ_res))
    print(f"  C mod {mod}: occ residues {occ_res}, false+ residues {fp_res}")
    if only_occ:
        print(f"    Residues ONLY in occupied: {only_occ}")
    if only_fp:
        print(f"    Residues ONLY in false+: {only_fp}")

print()
print("  VERDICT: No clean separation by C mod anything")
print()


# ============================================================
# TEST 10: a*b product sign/value
# ============================================================

print("=" * 80)
print("TEST 10: a*b product")
print("=" * 80)
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'a*b':>5s} {'sign':>6s} {'status':>8s}")
print(f"  {'-'*35}")

for n in passes:
    d = levels[n]
    ab = d['a'] * d['b']
    sign = "+" if ab > 0 else ("-" if ab < 0 else "0")
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {ab:>5d} {sign:>6s} {status:>8s}")

print()
occ_ab = sorted([levels[n]['a'] * levels[n]['b'] for n in true_pos])
fp_ab = sorted([levels[n]['a'] * levels[n]['b'] for n in false_pos])
print(f"  Occupied a*b: {occ_ab}")
print(f"  False+   a*b: {fp_ab}")
print()

# Check: false positive a*b = {-1, 0, 3}
# Occupied a*b = {-6, -4, 0, 0, 1, 2, 2, 4}
# Overlap at 0
overlap_ab = set(occ_ab) & set(fp_ab)
print(f"  Overlap: {sorted(overlap_ab) if overlap_ab else 'none'}")
print(f"  VERDICT: a*b does not cleanly separate (overlap at 0)")
print()


# ============================================================
# TEST 11: Fibonacci / Lucas number test
# ============================================================

print("=" * 80)
print("TEST 11: n as Fibonacci or Lucas number")
print("=" * 80)
print()

fibs_set = {0, 1, 1, 2, 3, 5, 8, 13, 21, 34}
lucas_set = {2, 1, 3, 4, 7, 11, 18, 29}

print(f"  Fibonacci (0-30): {sorted(fibs_set & set(range(31)))}")
print(f"  Lucas (0-30):     {sorted(lucas_set & set(range(31)))}")
print()

print(f"  {'n':>3s} {'Fib?':>5s} {'Luc?':>5s} {'Either':>7s} {'status':>8s}")
print(f"  {'-'*33}")

for n in passes:
    is_fib = n in fibs_set
    is_luc = n in lucas_set
    either = is_fib or is_luc
    status = "OCC" if levels[n]['is_occ'] else "FALSE+"
    print(f"  {n:>3d} {str(is_fib):>5s} {str(is_luc):>5s} {str(either):>7s} {status:>8s}")

print()

occ_fib_luc = sum(1 for n in true_pos if n in fibs_set or n in lucas_set)
fp_fib_luc = sum(1 for n in false_pos if n in fibs_set or n in lucas_set)
print(f"  Occupied that are Fib or Lucas: {occ_fib_luc}/8")
print(f"  False+   that are Fib or Lucas: {fp_fib_luc}/3")
print()

# Check which occupied are NOT Fib or Lucas
occ_not_fl = [n for n in true_pos if n not in fibs_set and n not in lucas_set]
fp_is_fl = [n for n in false_pos if n in fibs_set or n in lucas_set]
print(f"  Occupied NOT Fib/Lucas: {occ_not_fl}")
print(f"  False+ IS Fib/Lucas: {fp_is_fl}")
print(f"  VERDICT: Doesn't cleanly separate (16,17,19,26 aren't Fib/Lucas; 1,3 are Fib)")
print()


# ============================================================
# TEST 12: n*phi mod 1 (equidistribution / gap structure)
# ============================================================

print("=" * 80)
print("TEST 12: n * phi mod 1 (three-distance theorem / gap structure)")
print("=" * 80)
print()

print(f"  {'n':>3s} {'n*phi':>10s} {'n*phi%1':>10s} {'status':>8s}")
print(f"  {'-'*35}")

for n in passes:
    val = n * PHI
    frac = val % 1
    status = "OCC" if levels[n]['is_occ'] else "FALSE+"
    print(f"  {n:>3d} {val:>10.4f} {frac:>10.4f} {status:>8s}")

print()

occ_frac = sorted([(n * PHI) % 1 for n in true_pos])
fp_frac = sorted([(n * PHI) % 1 for n in false_pos])
print(f"  Occupied n*phi mod 1 (sorted): {[f'{x:.4f}' for x in occ_frac]}")
print(f"  False+   n*phi mod 1 (sorted): {[f'{x:.4f}' for x in fp_frac]}")
print()

# Check if false+ fractional parts fall in a "gap" of the occupied distribution
print("  Gap analysis: do false+ values fall in gaps between occupied values?")
occ_frac_sorted = sorted(occ_frac)
for fp_val in fp_frac:
    # Find where fp_val would be inserted
    pos = sum(1 for x in occ_frac_sorted if x < fp_val)
    if pos == 0:
        gap_left = occ_frac_sorted[-1] - 1  # wrap around
        gap_right = occ_frac_sorted[0]
    elif pos == len(occ_frac_sorted):
        gap_left = occ_frac_sorted[-1]
        gap_right = occ_frac_sorted[0] + 1  # wrap around
    else:
        gap_left = occ_frac_sorted[pos-1]
        gap_right = occ_frac_sorted[pos]
    gap_size = gap_right - gap_left
    print(f"    fp={fp_val:.4f}: between {gap_left:.4f} and {gap_right:.4f}, gap={gap_size:.4f}")

print()
print("  VERDICT: No obvious gap structure separates false positives")
print()


# ============================================================
# TEST 13: Powers of 2 on b=1 line
# ============================================================

print("=" * 80)
print("TEST 13: a-values on b=1 line: powers of 2?")
print("=" * 80)
print()

b1_occ = sorted([levels[n]['a'] for n in b1_levels if levels[n]['is_occ']])
b1_fp = sorted([levels[n]['a'] for n in b1_levels if not levels[n]['is_occ']])

print(f"  Occupied a-values on b=1: {b1_occ}")
print(f"  False+   a-values on b=1: {b1_fp}")
print()

# Check if occupied = {1, 2, 4} = {2^0, 2^1, 2^2}
print("  Occupied: 1 = 2^0, 2 = 2^1, 4 = 2^2")
is_pow2 = all((a & (a-1) == 0) and a > 0 for a in b1_occ)
print(f"  All powers of 2? {is_pow2}")
print()

# But is this meaningful? Only 3 points, and {1,2,4} being powers of 2
# could be coincidence. Let's check: for b=1 and a=integer, what are
# ALL possible n values with norm constraint?

print("  All levels on b=1 line with |N| in {0,1,5,19}:")
print(f"  {'a':>4s} {'n=5a+6':>8s} {'N':>4s} {'|N|':>4s} {'ok':>4s}")
for a in range(-10, 15):
    b = 1
    n = 5*a + 6*b
    N = zphi_norm_signed(a, b)
    ok = abs(N) in allowed_norms
    if ok and 0 <= n <= 30:
        is_occ_str = " <-- OCC" if n in occupied_n else ""
        is_fp_str = " <-- FALSE+" if n in [nn for nn in false_pos] else ""
        print(f"  {a:>4d} {n:>8d} {N:>4d} {abs(N):>4d} {'YES':>4s}{is_occ_str}{is_fp_str}")

print()
print("  The a-values that pass norm test on b=1: -1, 0, 1, 2, 4")
print("  Occupied: 1, 2, 4 (positive integers that are powers of 2)")
print("  Pattern: a is a positive power of 2 (including 2^0=1)")
print()


# ============================================================
# TEST 14: a as Fibonacci number on b=1 line
# ============================================================

print("=" * 80)
print("TEST 14: a as Fibonacci number on b=1 line")
print("=" * 80)
print()

fib_positives = {1, 2, 3, 5, 8, 13, 21}
for a in b1_occ:
    print(f"  a={a}: Fibonacci? {a in fib_positives}")
for a in b1_fp:
    print(f"  a={a}: Fibonacci? {a in fib_positives} (FALSE+)")

print()
print("  VERDICT: a=1 and a=2 are Fibonacci, but a=4 is NOT.")
print("  Also a=0 is not positive, a=-1 is not positive.")
print("  So Fibonacci doesn't perfectly characterize the b=1 pattern.")
print()


# ============================================================
# TEST 15: Geometric progression on b=1 line
# ============================================================

print("=" * 80)
print("TEST 15: Geometric progression on b=1 line")
print("=" * 80)
print()

print(f"  Occupied a-values on b=1: {b1_occ} = {{1, 2, 4}}")
print(f"  Ratios: 2/1 = {2/1}, 4/2 = {4/2}")
print(f"  Common ratio = 2")
print()
print("  Next in GP: 4*2 = 8. Level n = 5*8 + 6 = 46 (out of range)")
print("  Previous: 1/2 = 0.5 (not integer)")
print()
print("  OBSERVATION: {1, 2, 4} = 2^k for k=0,1,2")
print("  This is a GEOMETRIC progression with ratio 2.")
print("  It exhausts all powers of 2 with n = 5*2^k + 6 in [0,30]:")
for k in range(5):
    a = 2**k
    n = 5*a + 6
    in_range = 0 <= n <= 30
    print(f"    k={k}: a=2^{k}={a}, n={n}, in range? {in_range}")

print()


# ============================================================
# DEEP DIVE: Combined criteria search
# ============================================================

print("=" * 80)
print("DEEP DIVE: Exhaustive search for perfect discriminator")
print("=" * 80)
print()

# Build a rich feature set for each level
features = {}
for n in passes:
    d = levels[n]
    a, b = d['a'], d['b']
    N = d['N_signed']
    val = d['val']
    conj = d['conj']

    features[n] = {
        'a>=1': a >= 1,
        'a>=0': a >= 0,
        'a>0': a > 0,
        'b>=0': b >= 0,
        'b>=1': b >= 1,
        'b<=2': b <= 2,
        'a+b>=1': a+b >= 1,
        'a+b>=2': a+b >= 2,
        'N>=0': N >= 0,
        'N>0': N > 0,
        'val>=1': val >= 1,
        'val>1': val > 1,
        'conj>0': conj > 0,
        'conj>=0': conj >= 0,
        '|val|>=1': abs(val) >= 1,
        'a*b>=0': a*b >= 0,
        'a*b>0': a*b > 0,
        'C<=3': d['c'] <= 3,
        'C<=4': d['c'] <= 4,
        'C!=4': d['c'] != 4,
        'C%2==0': d['c'] % 2 == 0,
        'n%2==1': n % 2 == 1,
        'n%3!=1': n % 3 != 1,
        'n%5!=1': n % 5 != 1,
        '|a|<=|b|': abs(a) <= abs(b),
        'a<=b': a <= b,
        'a>=b': a >= b,
        '|a-b|<=1': abs(a-b) <= 1,
        '|a-b|<=2': abs(a-b) <= 2,
        'n>=3': n >= 3,
        'n>=5': n >= 5,
        'n!=1': n != 1,
        'n!=6': n != 6,
        'n!=23': n != 23,
        'val*conj>0': val * conj > 0,
        'val*conj>=0': val * conj >= 0,
        '|conj|<1': abs(conj) < 1,
        '|conj|<=1': abs(conj) <= 1,
        '|conj|>1': abs(conj) > 1,
        'on_a1': d['on_a1'],
        'on_b1': d['on_b1'],
        'on_325': d['on_325'],
        'is_origin': d['is_origin'],
        'a_pos_on_b1': d['on_b1'] and a >= 1,
        'b_small_on_a1': d['on_a1'] and b <= 2,
        'line_valid': (d['on_b1'] and a >= 1) or (d['on_a1'] and b <= 2) or d['on_325'] or d['is_origin'],
    }

# Find single criteria that perfectly separate
print("  SINGLE criteria that perfectly select 8 occupied from 11:")
print()

feat_names = sorted(features[passes[0]].keys())
perfect_singles = []

for fname in feat_names:
    occ_true = sum(1 for n in true_pos if features[n][fname])
    fp_true = sum(1 for n in false_pos if features[n][fname])

    # Perfect: all 8 occupied TRUE, all 3 false+ FALSE
    if occ_true == 8 and fp_true == 0:
        perfect_singles.append(fname)
        print(f"    PERFECT: {fname}  (occ: 8/8 True, fp: 0/3 True)")
    # Near-perfect
    elif occ_true >= 7 and fp_true == 0:
        print(f"    NEAR:    {fname}  (occ: {occ_true}/8 True, fp: {fp_true}/3 True)")
    elif occ_true == 8 and fp_true <= 1:
        print(f"    NEAR:    {fname}  (occ: {occ_true}/8 True, fp: {fp_true}/3 True)")

print()

if perfect_singles:
    print(f"  Found {len(perfect_singles)} perfect single criterion(a)!")
    for fname in perfect_singles:
        print(f"    {fname}")
        # Show details
        for n in passes:
            status = "OCC" if levels[n]['is_occ'] else "FALSE+"
            print(f"      n={n:>2d} ({levels[n]['a']:>2d},{levels[n]['b']:>2d}):"
                  f" {fname}={features[n][fname]} -> {status}")
    print()
else:
    print("  No perfect single criterion found. Trying pairs...")

# Try pairs
print("  PAIR criteria that perfectly select 8 occupied from 11:")
print()

perfect_pairs = []
for i, f1 in enumerate(feat_names):
    for f2 in feat_names[i+1:]:
        # f1 AND f2
        occ_and = sum(1 for n in true_pos if features[n][f1] and features[n][f2])
        fp_and = sum(1 for n in false_pos if features[n][f1] and features[n][f2])
        if occ_and == 8 and fp_and == 0:
            perfect_pairs.append((f1, 'AND', f2))

        # f1 OR f2
        occ_or = sum(1 for n in true_pos if features[n][f1] or features[n][f2])
        fp_or = sum(1 for n in false_pos if features[n][f1] or features[n][f2])
        if occ_or == 8 and fp_or == 0:
            perfect_pairs.append((f1, 'OR', f2))

if perfect_pairs:
    # Filter out trivially derived ones (pairs where one already works)
    nontrivial = [(f1, op, f2) for f1, op, f2 in perfect_pairs
                  if f1 not in perfect_singles and f2 not in perfect_singles]

    print(f"  Found {len(perfect_pairs)} perfect pairs ({len(nontrivial)} non-trivial):")
    for f1, op, f2 in perfect_pairs[:20]:
        trivial = " (trivial)" if f1 in perfect_singles or f2 in perfect_singles else ""
        print(f"    {f1} {op} {f2}{trivial}")
    if len(perfect_pairs) > 20:
        print(f"    ... and {len(perfect_pairs)-20} more")
else:
    print("  No perfect pairs found either!")

print()


# ============================================================
# TEST: The line_valid criterion in detail
# ============================================================

print("=" * 80)
print("DETAILED TEST: The 'line_valid' criterion")
print("=" * 80)
print()

print("  Rule: A level is occupied iff its (a,b) satisfies:")
print("    - Origin: (0,0)")
print("    - On b=1 line with a >= 1")
print("    - On a=1 line with b <= 2")
print("    - On 3a+2b=5 line (no extra constraint needed)")
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'line':>10s} {'constraint':>15s}"
      f" {'passes':>7s} {'actual':>8s} {'match':>6s}")
print(f"  {'-'*65}")

all_correct = True
for n in passes:
    d = levels[n]
    a, b = d['a'], d['b']

    if d['is_origin']:
        line = "origin"
        constraint = "(0,0)"
        passes_rule = True
    elif d['on_b1']:
        line = "b=1"
        constraint = f"a={a}>=1?"
        passes_rule = a >= 1
    elif d['on_a1']:
        line = "a=1"
        constraint = f"b={b}<=2?"
        passes_rule = b <= 2
    elif d['on_325']:
        line = "3a+2b=5"
        constraint = "always"
        passes_rule = True
    else:
        line = "none"
        constraint = "N/A"
        passes_rule = False

    actual = d['is_occ']
    match = passes_rule == actual
    if not match:
        all_correct = False

    print(f"  {n:>3d} ({a:>2d},{b:>2d}) {line:>10s} {constraint:>15s}"
          f" {str(passes_rule):>7s} {'OCC' if actual else 'no':>8s}"
          f" {'OK' if match else 'FAIL':>6s}")

print()
if all_correct:
    print("  *** PERFECT MATCH for all 11 levels! ***")
else:
    print("  Some mismatches found.")

print()


# ============================================================
# DEEPER: Can we unify a>=1 (b=1) and b<=2 (a=1)?
# ============================================================

print("=" * 80)
print("UNIFICATION: Can we express the line constraints uniformly?")
print("=" * 80)
print()

# On b=1 line: need a >= 1. In terms of (a,b): a >= b (since b=1, a >= 1 means a >= b)
# On a=1 line: need b <= 2. In terms of (a,b): b <= 2a (since a=1, b <= 2 means b <= 2*a)
# On 3a+2b=5 line: always OK. Check: occupied have (3,-2), (1,1), (-1,4).
#   For (3,-2): 3 >= -2 yes, -2 <= 6 yes
#   For (1,1): 1 >= 1 yes, 1 <= 2 yes
#   For (-1,4): -1 >= 4 NO.
# So a >= b doesn't work for line 3.

print("  Attempt 1: a >= b on each line?")
for n in passes:
    d = levels[n]
    a, b = d['a'], d['b']
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"    n={n:>2d} ({a:>2d},{b:>2d}): a>=b? {str(a>=b):>5s} {status}")

print()

# Check: b <= 2*|a|?
print("  Attempt 2: b <= 2*|a| (or origin)?")
for n in passes:
    d = levels[n]
    a, b = d['a'], d['b']
    rule = (b <= 2*abs(a)) or (a == 0 and b == 0)
    status = "OCC" if d['is_occ'] else "FALSE+"
    match = "OK" if rule == d['is_occ'] else "FAIL"
    print(f"    n={n:>2d} ({a:>2d},{b:>2d}): b<=2|a|? {str(rule):>5s} {status:>8s} {match}")

print()

# Check: |a| + |b| <= 2*max(|a|, |b|) is always true (trivial)
# What about: min(|a|, |b|) <= some function of max?

# Try: val = a + b*phi > 0? (positivity of the Z[phi] element)
print("  Attempt 3: val = a + b*phi > 0 (OR origin)?")
for n in passes:
    d = levels[n]
    rule = d['val'] > 0 or d['is_origin']
    status = "OCC" if d['is_occ'] else "FALSE+"
    match = "OK" if rule == d['is_occ'] else "FAIL"
    print(f"    n={n:>2d} ({d['a']:>2d},{d['b']:>2d}): val={d['val']:>7.4f}>0? {str(rule):>5s}"
          f" {status:>8s} {match}")

print()

# Check: val > 0 AND conj > some value?
print("  Attempt 4: val > 0 AND conj > -2?")
for n in passes:
    d = levels[n]
    rule = (d['val'] > 0 and d['conj'] > -2) or d['is_origin']
    status = "OCC" if d['is_occ'] else "FALSE+"
    match = "OK" if rule == d['is_occ'] else "FAIL"
    print(f"    n={n:>2d} ({d['a']:>2d},{d['b']:>2d}): val={d['val']:>7.4f},"
          f" conj={d['conj']:>7.4f} {str(rule):>5s} {status:>8s} {match}")

print()

# Try the KEY insight: b <= 2*a for a >= 1 AND a >= b for b >= 1?
# That's equivalent to: on the first quadrant (a>=1, b>=1), need b <= 2a
# Plus a=1 allows b=0 (down quark)

# Actually let me try: a + b*phi^(-1) > 0  (i.e., a + b/phi > 0)
print("  Attempt 5: a + b/phi > 0 (OR origin)?")
for n in passes:
    d = levels[n]
    test_val = d['a'] + d['b'] / PHI
    rule = test_val > 0 or d['is_origin']
    status = "OCC" if d['is_occ'] else "FALSE+"
    match = "OK" if rule == d['is_occ'] else "FAIL"
    print(f"    n={n:>2d} ({d['a']:>2d},{d['b']:>2d}): a+b/phi={test_val:>7.4f}>0?"
          f" {str(rule):>5s} {status:>8s} {match}")

print()

# Try: the CONJUGATE value a + b*phi' where phi' = (1-sqrt(5))/2
# conj = a + b*phi'. We need: val > 0 AND (something about conj?)
# Let's just check all with |conj| < some threshold?

# Actually, let me try a very different approach:
# What if the rule is about the NORM SIGN combined with position?
print("  Attempt 6: N >= 0 (signed norm)?")
for n in passes:
    d = levels[n]
    rule = d['N_signed'] >= 0
    status = "OCC" if d['is_occ'] else "FALSE+"
    match = "OK" if rule == d['is_occ'] else "FAIL"
    print(f"    n={n:>2d} ({d['a']:>2d},{d['b']:>2d}): N={d['N_signed']:>4d}>=0?"
          f" {str(rule):>5s} {status:>8s} {match}")

print()

# Check if N >= 0 works:
n_ge0_occ = sum(1 for n in true_pos if levels[n]['N_signed'] >= 0)
n_ge0_fp = sum(1 for n in false_pos if levels[n]['N_signed'] >= 0)
print(f"  N>=0: occupied pass {n_ge0_occ}/8, false+ pass {n_ge0_fp}/3")
print()


# ============================================================
# TEST: The conj (Galois conjugate) and val relationship
# ============================================================

print("=" * 80)
print("TEST: Galois conjugate analysis")
print("=" * 80)
print()

# For a + b*phi, the Galois conjugate is a + b*phi' where phi' = (1-sqrt(5))/2
# The norm is val * conj = a^2 + ab - b^2
# Note: if |N| = 1, then val * conj = +/- 1
# If N = -1, then val and conj have OPPOSITE signs
# If N = +1, then val and conj have SAME sign

print("  val * conj sign analysis:")
print(f"  {'n':>3s} {'(a,b)':>7s} {'val':>8s} {'conj':>8s} {'v*c':>8s}"
      f" {'same_sign':>10s} {'status':>8s}")
print(f"  {'-'*60}")

for n in passes:
    d = levels[n]
    vc = d['val'] * d['conj']
    same = (d['val'] > 0 and d['conj'] > 0) or (d['val'] < 0 and d['conj'] < 0)
    status = "OCC" if d['is_occ'] else "FALSE+"
    same_str = "same" if same else "opposite"
    if abs(d['val']) < 1e-10 or abs(d['conj']) < 1e-10:
        same_str = "zero"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {d['val']:>8.4f} {d['conj']:>8.4f}"
          f" {vc:>8.4f} {same_str:>10s} {status:>8s}")

print()


# ============================================================
# TEST: val > |conj| ?
# ============================================================

print("=" * 80)
print("TEST: val > |conj| (dominant embedding)")
print("=" * 80)
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'val':>8s} {'|conj|':>8s} {'val>|c|':>8s} {'status':>8s}")
print(f"  {'-'*50}")

for n in passes:
    d = levels[n]
    dominant = d['val'] > abs(d['conj'])
    status = "OCC" if d['is_occ'] else "FALSE+"
    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {d['val']:>8.4f} {abs(d['conj']):>8.4f}"
          f" {str(dominant):>8s} {status:>8s}")

print()

dom_occ = sum(1 for n in true_pos if levels[n]['val'] > abs(levels[n]['conj']))
dom_fp = sum(1 for n in false_pos if levels[n]['val'] > abs(levels[n]['conj']))
print(f"  val > |conj|: occupied pass {dom_occ}/8, false+ pass {dom_fp}/3")
print()


# ============================================================
# TEST: |val/conj| ratio (for non-zero)
# ============================================================

print("=" * 80)
print("TEST: |val/conj| ratio")
print("=" * 80)
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'|v/c|':>10s} {'ln|v/c|':>10s} {'status':>8s}")
print(f"  {'-'*45}")

for n in passes:
    d = levels[n]
    if abs(d['conj']) > 1e-10 and abs(d['val']) > 1e-10:
        ratio = abs(d['val'] / d['conj'])
        ln_ratio = np.log(ratio)
    else:
        ratio = float('inf')
        ln_ratio = float('inf')
    status = "OCC" if d['is_occ'] else "FALSE+"
    if ratio < 1e6:
        print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {ratio:>10.4f} {ln_ratio:>10.4f} {status:>8s}")
    else:
        print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {'inf':>10s} {'inf':>10s} {status:>8s}")

print()

# For |N|=1 units, |val/conj| = phi^(2k) where val = +-phi^k
# So ln|v/c| = 2k * ln(phi)
print("  For |N|=1: ln|v/c| / (2*ln(phi)) should give k:")
for n in passes:
    d = levels[n]
    if d['N_abs'] == 1 and abs(d['conj']) > 1e-10 and abs(d['val']) > 1e-10:
        ratio = abs(d['val'] / d['conj'])
        k_eff = np.log(ratio) / (2 * np.log(PHI))
        status = "OCC" if d['is_occ'] else "FALSE+"
        print(f"    n={n:>2d}: k_eff = {k_eff:>6.2f} {status}")

print()


# ============================================================
# COMBINATORIAL: b <= 2|a| test (revisited carefully)
# ============================================================

print("=" * 80)
print("COMBINATORIAL TEST: Various a,b inequalities")
print("=" * 80)
print()

tests = [
    ("b <= 2|a|", lambda a, b: b <= 2*abs(a)),
    ("b <= 2a", lambda a, b: b <= 2*a),
    ("|b| <= 2|a|", lambda a, b: abs(b) <= 2*abs(a)),
    ("a^2 >= b", lambda a, b: a**2 >= b),
    ("a^2 >= |b|", lambda a, b: a**2 >= abs(b)),
    ("a >= 0 or |a|<=|b|", lambda a, b: a >= 0 or abs(a) <= abs(b)),
    ("a+b >= 0", lambda a, b: a + b >= 0),
    ("|a+b| <= |a|+|b|/2", lambda a, b: abs(a+b) <= (abs(a)+abs(b))/2 + 0.01),
    ("a >= b/3", lambda a, b: a >= b/3 - 0.01),
    ("2a + b >= 1", lambda a, b: 2*a + b >= 1),
    ("3a + b >= 1", lambda a, b: 3*a + b >= 1),
    ("a + 2b >= 1", lambda a, b: a + 2*b >= 1),
]

for test_name, test_fn in tests:
    occ_pass = sum(1 for n in true_pos
                   if test_fn(levels[n]['a'], levels[n]['b']) or levels[n]['is_origin'])
    fp_pass = sum(1 for n in false_pos
                  if test_fn(levels[n]['a'], levels[n]['b']) or levels[n]['is_origin'])
    marker = ""
    if occ_pass == 8 and fp_pass == 0:
        marker = " <-- PERFECT!"
    elif occ_pass >= 7 and fp_pass <= 1:
        marker = " <-- near"
    print(f"  {test_name:>25s}: occ {occ_pass}/8, fp {fp_pass}/3{marker}")

print()


# ============================================================
# KEY ANALYSIS: What makes n=23 different from n=17?
# ============================================================

print("=" * 80)
print("KEY ANALYSIS: n=23 (1,3) vs n=17 (1,2) - both on a=1 line")
print("=" * 80)
print()

n17 = levels[17]
n23 = levels[23]

print(f"  n=17 (tau):    a={n17['a']}, b={n17['b']}, C={n17['c']}")
print(f"    val = {n17['val']:.6f}, conj = {n17['conj']:.6f}")
print(f"    N = {n17['N_signed']}, |N| = {n17['N_abs']}")
print(f"    val * conj = {n17['val']*n17['conj']:.6f}")
print()

print(f"  n=23 (false+): a={n23['a']}, b={n23['b']}, C={n23['c']}")
print(f"    val = {n23['val']:.6f}, conj = {n23['conj']:.6f}")
print(f"    N = {n23['N_signed']}, |N| = {n23['N_abs']}")
print(f"    val * conj = {n23['val']*n23['conj']:.6f}")
print()

print("  Differences:")
print(f"    b: 2 vs 3 (tau has smaller b)")
print(f"    N: {n17['N_signed']} vs {n23['N_signed']} (signs differ!)")
print(f"    val: {n17['val']:.4f} vs {n23['val']:.4f}")
print(f"    conj: {n17['conj']:.4f} vs {n23['conj']:.4f} (tau conj > 0, n23 conj < 0)")
print()

# Key insight: for n=17 (1,2), conj = 1 + 2*phi' = 1 + 2*(-0.618) = -0.236
# For n=23 (1,3), conj = 1 + 3*phi' = 1 + 3*(-0.618) = -0.854
# Both negative. But N(1,2) = 1+2-4 = -1, N(1,3) = 1+3-9 = -5.
# Hmm, both negative N too.

# What about the MASS? m = m_e * phi^n
print("  Mass comparison:")
m_e = 0.511  # MeV
for n_val in [17, 23]:
    mass = m_e * PHI**(n_val)
    print(f"    n={n_val}: m = {mass:.2f} MeV = {mass/1000:.4f} GeV")

print(f"    Tau mass (exp): 1776.86 MeV = 1.777 GeV")
print(f"    phi^17 = {PHI**17:.2f}, phi^23 = {PHI**23:.2f}")
print(f"    There is no known fermion near {m_e * PHI**23:.0f} MeV")
print()


# ============================================================
# FINAL: Check N >= 0 as THE criterion
# ============================================================

print("=" * 80)
print("FINAL CRITICAL TEST: N(a+b*phi) >= 0 (signed norm non-negative)")
print("=" * 80)
print()

print("  The signed norm N = a^2 + ab - b^2")
print("  This equals (a + b*phi) * (a + b*phi') where phi' = conjugate")
print()

print(f"  {'n':>3s} {'(a,b)':>7s} {'N':>4s} {'N>=0':>6s} {'actual':>8s} {'match':>6s}")
print(f"  {'-'*40}")

n_ge0_correct = 0
n_ge0_total = 0

for n in passes:
    d = levels[n]
    N = d['N_signed']
    predicted = N >= 0
    actual = d['is_occ']
    match = predicted == actual
    if match:
        n_ge0_correct += 1
    n_ge0_total += 1

    print(f"  {n:>3d} ({d['a']:>2d},{d['b']:>2d}) {N:>4d} {str(predicted):>6s}"
          f" {'OCC' if actual else 'no':>8s} {'OK' if match else 'FAIL':>6s}")

print()
print(f"  Score: {n_ge0_correct}/{n_ge0_total}")

if n_ge0_correct == n_ge0_total:
    print()
    print("  ***********************************************************")
    print("  * N >= 0 PERFECTLY SEPARATES occupied from false positives *")
    print("  ***********************************************************")
    print()
    print("  COMPLETE SELECTION RULE (3 conditions):")
    print("    1. |N(a+b*phi)| in {0, 1, 5, 19}   [norm magnitude]")
    print("    2. (a,b) on {origin} U {a=1} U {b=1} U {3a+2b=5}  [3 lines]")
    print("    3. N(a+b*phi) >= 0                  [norm sign]")
    print()
    print("  Verification: false positives all have N < 0:")
    for n in false_pos:
        d = levels[n]
        print(f"    n={n}: ({d['a']},{d['b']}), N = {d['N_signed']} < 0")
    print()
    print("  And ALL occupied have N >= 0:")
    for n in true_pos:
        d = levels[n]
        print(f"    n={n}: ({d['a']},{d['b']}), N = {d['N_signed']} >= 0")
else:
    # Check which fail
    for n in passes:
        d = levels[n]
        N = d['N_signed']
        predicted = N >= 0
        actual = d['is_occ']
        if predicted != actual:
            print(f"  MISMATCH: n={n} ({d['a']},{d['b']}): N={N}, predicted={'OCC' if predicted else 'no'},"
                  f" actual={'OCC' if actual else 'no'}")

print()


# ============================================================
# PHYSICAL INTERPRETATION
# ============================================================

print("=" * 80)
print("PHYSICAL INTERPRETATION")
print("=" * 80)
print()

print("  The signed norm N = a^2 + ab - b^2 = (a + b*phi)(a + b*phi')")
print()
print("  N >= 0 means: val and conj have the SAME SIGN (or one is zero)")
print("  N < 0 means: val and conj have OPPOSITE SIGNS")
print()
print("  Geometrically in Z[phi]:")
print("    N >= 0: the element a+b*phi lies in a 'hyperbolic sector'")
print("    where both Galois embeddings agree on sign.")
print("    This is related to TOTALLY POSITIVE elements in the ring.")
print()
print("  For phi-powers: phi^k has N = (-1)^k")
print("    Even k: N = +1 >= 0 (occupied: phi^2 = val for n=17)")
print("    Odd k:  N = -1 < 0  (excluded: phi^1 = val for n=6, phi^(-1) for n=1)")
print()

# Actually verify this
print("  Phi-power parity check:")
for n in passes:
    d = levels[n]
    if d['N_abs'] == 1:
        val = d['val']
        if abs(val) > 1e-10:
            k = round(np.log(abs(val)) / np.log(PHI))
            sign_val = 1 if val > 0 else -1
            # The actual unit is sign_val * phi^k
            # N of phi^k = (-1)^k
            # N of -phi^k = (-1)^k (negation doesn't change norm)
            expected_N = int((-1)**k)
            status = "OCC" if d['is_occ'] else "FALSE+"
            print(f"    n={n:>2d}: val = {'+' if sign_val>0 else '-'}phi^{k},"
                  f" N = (-1)^{k} = {expected_N:+d},"
                  f" N check: {d['N_signed']:+d} {status}")

print()
print("  KEY INSIGHT: For |N|=1 (units), N >= 0 selects EVEN phi-powers!")
print("  phi^(even) has norm +1 (totally positive unit)")
print("  phi^(odd) has norm -1 (NOT totally positive)")
print()


# ============================================================
# ANALYSIS: Why N >= 0 FAILS
# ============================================================

print("=" * 80)
print("ANALYSIS: Why N >= 0 fails as a criterion")
print("=" * 80)
print()

print("  N >= 0 eliminates ALL 3 false positives (all have N < 0)")
print("  But it ALSO eliminates 3 occupied levels:")
occ_neg_N = [(n, levels[n]['a'], levels[n]['b'], levels[n]['N_signed'])
             for n in true_pos if levels[n]['N_signed'] < 0]
for n, a, b, N in occ_neg_N:
    print(f"    n={n}: ({a},{b}), N={N}")
print()
print("  So N >= 0 is too RESTRICTIVE: it kills 3 real fermions.")
print("  Score: 8/11 (3 false negatives)")
print()

# Algebraic identity for reference
print("  For reference: 4*N = (2a+b)^2 - 5*b^2")
print("  Verification:")
for n in passes:
    d = levels[n]
    a, b = d['a'], d['b']
    lhs = 4 * d['N_signed']
    rhs = (2*a+b)**2 - 5*b**2
    ok = "OK" if lhs == rhs else "FAIL"
    print(f"    n={n:>2d}: 4N = {lhs:>4d}, (2a+b)^2 - 5b^2 = {rhs:>4d} {ok}")

print()

# ============================================================
# VERIFY line_valid ON FULL RANGE 0-30
# ============================================================

print("=" * 80)
print("FULL VERIFICATION: line_valid rule on n = 0..30")
print("=" * 80)
print()

# line_valid: on {origin} or {b=1 with a>=1} or {a=1 with b<=2} or {3a+2b=5}
def line_valid(a, b):
    if a == 0 and b == 0:
        return True  # origin
    if b == 1 and a >= 1:
        return True  # b=1 line with a>=1
    if a == 1 and b <= 2:
        return True  # a=1 line with b<=2
    if 3*a + 2*b == 5:
        return True  # third line
    return False

print(f"  {'n':>3s} {'(a,b)':>7s} {'|N|':>4s} {'|N|ok':>6s}"
      f" {'line':>6s} {'lv':>4s} {'pred':>5s} {'act':>5s} {'match':>6s}")
print(f"  {'-'*55}")

total_correct = 0
total = 0

for n in range(31):
    d = levels[n]
    a, b = d['a'], d['b']

    norm_ok = d['norm_ok']
    on_line = d['on_line']
    lv = line_valid(a, b)

    predicted = norm_ok and lv
    actual = d['is_occ']
    match = predicted == actual

    if match:
        total_correct += 1
    total += 1

    if norm_ok or actual:  # only show interesting rows
        print(f"  {n:>3d} ({a:>2d},{b:>2d}) {d['N_abs']:>4d}"
              f" {'Y' if norm_ok else 'N':>6s} {'Y' if on_line else 'N':>6s}"
              f" {'Y' if lv else 'N':>4s}"
              f" {'OCC' if predicted else '---':>5s}"
              f" {'OCC' if actual else '---':>5s}"
              f" {'OK' if match else 'FAIL':>6s}")

print()
print(f"  Full range score: {total_correct}/{total}")
if total_correct == total:
    print("  *** PERFECT SCORE on all 31 levels! ***")

print()


# ============================================================
# SUMMARY
# ============================================================

print("=" * 80)
print("SUMMARY OF FINDINGS")
print("=" * 80)
print()

print("  THE 3 FALSE POSITIVES:")
print("  ----------------------")
for n in false_pos:
    d = levels[n]
    print(f"  n={n}: ({d['a']},{d['b']})")
    print(f"    val = {d['val']:.4f}, conj = {d['conj']:.4f}")
    print(f"    N = {d['N_signed']} (NEGATIVE)")
    print(f"    Why excluded: val and conj have OPPOSITE signs")
    if d['N_abs'] == 1:
        k = round(np.log(abs(d['val'])) / np.log(PHI))
        print(f"    Unit: phi^{k} (ODD power -> N = -1)")
    print()

print()
print("  WHAT SEPARATES THEM:")
print("  --------------------")
print()

print("  TEST RESULTS SUMMARY:")
print(f"  {'Test':>40s} {'Separates?':>12s}")
print(f"  {'-'*55}")

test_results = [
    ("Phi-power parity (even k)", "PARTIAL (|N|=1 only)"),
    ("line_valid (position on line)", "*** PERFECT ***"),
    ("a >= 0 AND b >= 0", "NO"),
    ("Signed norm N >= 0", "NO (kills 3 occ)"),
    ("Complexity C alone", "NO (overlap)"),
    ("|val - conj|", "NO (= |b|*sqrt(5))"),
    ("a > 0", "NO (misses e, b-quark)"),
    ("max(|a|,|b|) threshold", "NO"),
    ("C mod anything", "NO"),
    ("a*b product", "NO (overlap at 0)"),
    ("Fibonacci/Lucas", "NO"),
    ("n*phi mod 1 gap", "NO"),
    ("b=1 line: a = 2^k", "YES (on b=1 only)"),
    ("val > |conj|", "NO (5/8 occ, 2/3 fp)"),
    ("b <= 2|a| or origin", "NEAR (7/8, 1/3)"),
    ("b <= 2a or origin", "NEAR (7/8, 0/3)"),
]

for test, result in test_results:
    print(f"  {test:>40s}  {result}")

print()
print("  THE KEY FINDING:")
print("  ================")
print()
print("  The ONLY single criterion that perfectly separates all 8 occupied")
print("  from all 3 false positives is 'line_valid':")
print()
print("  *** COMPLETE SELECTION RULE (2 conditions): ***")
print()
print("  A level n (with (a,b) = argmin|a|+|b| s.t. 5a+6b=n) is occupied iff:")
print("    (1) |N(a+b*phi)| in {0, 1, 5, 19}   [norm magnitude constraint]")
print("    (2) (a,b) satisfies the LINE-POSITION constraint:")
print("        - (0,0) origin                    [electron]")
print("        - b = 1 with a >= 1               [mu/s, charm, top]")
print("        - a = 1 with b <= 2               [down, mu/s, tau]")
print("        - 3a + 2b = 5 (no restriction)    [up, mu/s, bottom]")
print()
print("  WHY each false positive is excluded:")
print("    n=1:  (-1,1) on b=1 but a=-1 < 1  --> EXCLUDED")
print("    n=6:  (0,1)  on b=1 but a=0  < 1  --> EXCLUDED")
print("    n=23: (1,3)  on a=1 but b=3  > 2  --> EXCLUDED")
print()
print("  INTERPRETATION:")
print("  - On b=1 line: a >= 1 means the 'diameter quantum number' is positive")
print("  - On a=1 line: b <= 2 means the 'decagon quantum number' is bounded")
print("  - On 3a+2b=5: the line is already finite (3 integer points)")
print("  - The 3 lines meet at (1,1) = muon/strange (triple junction)")
print()

# Check: is this just counting? How many integer points on each line?
print("  Integer points per line with n in [0,30] and |N| in {0,1,5,19}:")
for line_name, line_fn, extra_fn in [
    ("b=1", lambda a, b: b == 1, lambda a, b: True),
    ("b=1 + a>=1", lambda a, b: b == 1, lambda a, b: a >= 1),
    ("a=1", lambda a, b: a == 1, lambda a, b: True),
    ("a=1 + b<=2", lambda a, b: a == 1, lambda a, b: b <= 2),
    ("3a+2b=5", lambda a, b: 3*a + 2*b == 5, lambda a, b: True),
]:
    count = 0
    points = []
    for n in range(31):
        d = levels[n]
        a, b = d['a'], d['b']
        if d['norm_ok'] and line_fn(a, b) and extra_fn(a, b):
            count += 1
            points.append((n, a, b))
    print(f"    {line_name:>15s}: {count} points -> {points}")

print()
print("  ADDITIONAL PATTERN on b=1 line:")
print("    Occupied a-values: {1, 2, 4} = {2^0, 2^1, 2^2}")
print("    This is a GEOMETRIC PROGRESSION with ratio 2.")
print("    Could be coincidence (only 3 points) or deep structure.")
print()
print("  N >= 0 DOES NOT WORK because 3 occupied levels have N < 0:")
occ_neg = [(n, levels[n]['a'], levels[n]['b'], levels[n]['N_signed'])
           for n in true_pos if levels[n]['N_signed'] < 0]
for n, a, b, N in occ_neg:
    print(f"    n={n}: ({a},{b}), N={N}")
print("  (u(3,-2), tau(1,2), b(-1,4) all have negative signed norm)")
print()
print("  CLASSIFICATION: PATTERN (not yet derived from first principles)")
print("  The line_valid rule is DESCRIPTIVE: it perfectly fits the data")
print("  but the constraint a>=1 / b<=2 has no known derivation.")
print()

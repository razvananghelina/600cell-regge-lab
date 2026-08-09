"""
EXP-121b: Z[phi] Norms and Average Level = 12
================================================

Two remarkable observations from exp121:
1. The AVERAGE occupied level = 12 = vertex degree of 600-cell (EXACT)
2. The Z[phi] norms of occupied fermions are ONLY {0, 1, 5, 19}

This experiment investigates:
A. Is average = 12 a CONSTRAINT or a COINCIDENCE?
B. Why only norms {0, 1, 5, 19} in Z[phi]?
C. Number theory of Z[phi] = Z[(1+sqrt(5))/2]
D. What makes 1, 5, 19 special primes/units in Z[phi]?
E. Can the norm constraint SELECT the occupied levels?

Author: Claude (exp121b)
Date: 2026-02-08
"""

import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
from math import gcd

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2

# Fermion data
fermions = {
    'e':   (0, 0, 0),    # (a, b, n)
    'mu':  (1, 1, 11),
    'tau': (1, 2, 17),
    'u':   (3, -2, 3),
    'c':   (2, 1, 16),
    't':   (4, 1, 26),
    'd':   (1, 0, 5),
    's':   (1, 1, 11),
    'b_q': (-1, 4, 19),
}

occupied_n = sorted(set(f[2] for f in fermions.values()))

def min_complexity(n):
    best = None
    best_c = float('inf')
    for a in range(-20, 21):
        if (n - 5*a) % 6 == 0:
            b = (n - 5*a) // 6
            c = abs(a) + abs(b)
            if c < best_c:
                best_c = c
                best = (a, b)
    return best, best_c

def zphi_norm(a, b):
    """Norm of a + b*phi in Z[phi]: N = a^2 + ab - b^2."""
    return a**2 + a*b - b**2

print("=" * 80)
print("EXP-121b: Z[phi] Norms and Average Level = 12")
print("=" * 80)
print()

# ============================================================
# SECTION 1: The average = 12 property
# ============================================================

print("-" * 80)
print("SECTION 1: Average occupied level = 12")
print("-" * 80)
print()

all_n = [f[2] for f in fermions.values()]  # with multiplicity
sum_n = sum(all_n)
avg_n = sum_n / len(all_n)

print(f"  9 fermion levels (with mu/s degeneracy): {sorted(all_n)}")
print(f"  Sum = {sum_n}")
print(f"  Average = {sum_n}/{len(all_n)} = {avg_n:.4f}")
print(f"  Vertex degree of 600-cell = 12")
print(f"  EXACT MATCH: {avg_n == 12.0}")
print()

# Is this forced by some structural constraint?
# Test: how many 9-element multisets from {0,...,26} have average = 12?

# We need 9 levels with sum = 108.
# But we also have the constraint:
# - 3 leptons, 3 up-type, 3 down-type
# - mu and s share n=11 (same (a,b)=(1,1))
# - electron is at n=0 (ground state)

# First, the unconstrained count: 9 numbers from {0,...,26} summing to 108
# This is a partition problem. Let's use Monte Carlo sampling.

print("  Monte Carlo: how common is sum=108 for 9 random levels from {0,...,26}?")
np.random.seed(42)
n_trials = 1000000
count_108 = 0
total_sums = []
for _ in range(n_trials):
    levels = np.random.choice(27, 9, replace=True)
    s = levels.sum()
    total_sums.append(s)
    if s == 108:
        count_108 += 1

total_sums = np.array(total_sums)
print(f"  Trials: {n_trials}")
print(f"  Sum = 108 count: {count_108} ({count_108/n_trials*100:.2f}%)")
print(f"  Expected sum (uniform): {13 * 9} = {13*9}")
print(f"  Sum 108 vs expected 117: our sum is LOWER")
print(f"  Mean sampled sum: {total_sums.mean():.1f}")
print(f"  Std sampled sum: {total_sums.std():.1f}")
print()

# More restrictive: 3 sectors of 3, each from {0,...,26}
print("  Monte Carlo: 3 sectors of 3, sum=108?")
count_108_sector = 0
for _ in range(n_trials):
    lep = sorted(np.random.choice(27, 3, replace=False))
    up = sorted(np.random.choice(27, 3, replace=False))
    down = sorted(np.random.choice(27, 3, replace=False))
    s = sum(lep) + sum(up) + sum(down)
    if s == 108:
        count_108_sector += 1

print(f"  Sum = 108 count: {count_108_sector} ({count_108_sector/n_trials*100:.2f}%)")
print()

# ============================================================
# SECTION 2: Z[phi] norm analysis
# ============================================================

print("-" * 80)
print("SECTION 2: Z[phi] norm analysis")
print("-" * 80)
print()

# The norm form: N(a + b*phi) = a^2 + ab - b^2
# This is the determinant of the 2x2 matrix representation.
# Properties:
# - N is multiplicative: N(xy) = N(x)N(y)
# - N(phi) = -1 (phi is a unit of norm -1)
# - N(1) = 1
# - Units of Z[phi]: +-phi^n, which have norm +-1

print("  Z[phi] = Z[(1+sqrt(5))/2] = ring of integers of Q(sqrt(5))")
print()
print("  Norm form: N(a + b*phi) = a^2 + ab - b^2")
print()

# Norms for ALL levels 0-26
print(f"  {'n':>4s} {'(a,b)':>8s} {'a+b*phi':>10s} {'N(a+b*phi)':>12s} {'|N|':>5s} {'occ':>4s}")
print(f"  {'-'*48}")

all_norms = {}
for n in range(27):
    (a, b), c = min_complexity(n)
    val = a + b * PHI
    N = zphi_norm(a, b)
    is_occ = '*' if n in occupied_n else ' '
    all_norms[n] = N
    print(f"  {n:>4d} ({a:>2d},{b:>2d}) {val:>10.4f} {N:>12d} {abs(N):>5d} {is_occ}")

print()

# Unique absolute norms
occ_abs_norms = sorted(set(abs(all_norms[n]) for n in occupied_n))
all_abs_norms = sorted(set(abs(all_norms[n]) for n in range(27)))
print(f"  Occupied |norms|: {occ_abs_norms}")
print(f"  All |norms| (0-26): {all_abs_norms}")
print(f"  Norms ONLY in unoccupied: {sorted(set(all_abs_norms) - set(occ_abs_norms))}")
print()

# ============================================================
# SECTION 3: Number theory of norms {0, 1, 5, 19}
# ============================================================

print("-" * 80)
print("SECTION 3: What's special about {0, 1, 5, 19}?")
print("-" * 80)
print()

# Prime factorization behavior in Z[phi]:
# Z[phi] is a PID (principal ideal domain), class number = 1
# A rational prime p behaves as follows in Z[phi]:
# - p = 5: RAMIFIES (5 = -sqrt(5)^2 * unit), since disc(Q(sqrt(5))) = 5
# - p = 2: INERT (since 5 mod 8 = 5, Legendre symbol (5/2) = ?)
# Actually, for Q(sqrt(5)):
# - p = 5: ramifies
# - p = ±1 mod 5: splits (Legendre symbol (5/p) = 1)
# - p = ±2 mod 5: inert (Legendre symbol (5/p) = -1)

print("  Prime behavior in Z[phi] = Z[(1+sqrt(5))/2]:")
print("  (Ring of integers of Q(sqrt(5)), discriminant = 5)")
print()

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    p_mod_5 = p % 5
    if p == 5:
        status = "RAMIFIES"
    elif p_mod_5 in [1, 4]:  # ±1 mod 5
        status = "SPLITS"
    else:  # ±2 mod 5
        status = "INERT"
    print(f"  p = {p:>3d}: p mod 5 = {p_mod_5}, {status}")

print()

# Our norms: 0, 1, 5, 19
print("  Analysis of occupied norms:")
print(f"  |N| = 0: trivial (electron at origin)")
print(f"  |N| = 1: UNIT in Z[phi] (element is +-phi^k)")
print(f"  |N| = 5: RAMIFIED prime (special in Q(sqrt(5)))")
print(f"  |N| = 19: 19 mod 5 = 4 = -1 mod 5 => SPLITS in Z[phi]")
print()

# Factorization of 19 in Z[phi]:
# 19 splits as 19 = pi * pi' where pi, pi' are conjugate primes
# pi has norm 19. Find pi:
print("  Finding pi such that N(pi) = 19:")
found_19 = []
for a in range(-10, 11):
    for b in range(-10, 11):
        if zphi_norm(a, b) == 19:
            found_19.append((a, b, a + b * PHI))
        elif zphi_norm(a, b) == -19:
            found_19.append((a, b, a + b * PHI))

print(f"  Elements with |N| = 19:")
for a, b, val in sorted(found_19, key=lambda x: abs(x[2])):
    N = zphi_norm(a, b)
    print(f"    {a} + {b}*phi = {val:.4f}, N = {N}")

print()

# Similarly for norm 5
print("  Elements with |N| = 5:")
found_5 = []
for a in range(-10, 11):
    for b in range(-10, 11):
        N = zphi_norm(a, b)
        if abs(N) == 5:
            found_5.append((a, b, a + b * PHI, N))

for a, b, val, N in sorted(found_5, key=lambda x: abs(x[2]))[:10]:
    print(f"    {a} + {b}*phi = {val:.4f}, N = {N}")

print()

# ============================================================
# SECTION 4: The norm constraint as selection rule
# ============================================================

print("-" * 80)
print("SECTION 4: Can norms SELECT the occupied levels?")
print("-" * 80)
print()

# Hypothesis: occupied levels are those where |N(a+b*phi)| is in {0, 1, 5, 19}
# Test this!

predicted_occ = [n for n in range(27) if abs(all_norms[n]) in {0, 1, 5, 19}]
actual_occ = occupied_n

print(f"  Levels with |N| in {{0, 1, 5, 19}}: {predicted_occ}")
print(f"  Actual occupied levels:              {actual_occ}")
print(f"  Match: {predicted_occ == actual_occ}")
print()

if predicted_occ != actual_occ:
    # What's the overlap?
    pred_set = set(predicted_occ)
    act_set = set(actual_occ)
    print(f"  Intersection: {sorted(pred_set & act_set)}")
    print(f"  In predicted but not occupied: {sorted(pred_set - act_set)}")
    print(f"  In occupied but not predicted: {sorted(act_set - pred_set)}")
    print()

# What if we try |N| in {0, 1, 5} only?
predicted_015 = [n for n in range(27) if abs(all_norms[n]) in {0, 1, 5}]
print(f"  Levels with |N| in {{0, 1, 5}}: {predicted_015}")
print()

# What if we require N itself (not |N|) to be in a specific set?
predicted_signed = [n for n in range(27) if all_norms[n] in {0, 1, 5, -1, 19, -19}]
print(f"  Levels with N in {{0, +-1, 5, +-19}}: {predicted_signed}")
print()

# ============================================================
# SECTION 5: Distribution of norms
# ============================================================

print("-" * 80)
print("SECTION 5: Norm distribution analysis")
print("-" * 80)
print()

# Group levels by their absolute norm
norm_groups = defaultdict(list)
for n in range(27):
    norm_groups[abs(all_norms[n])].append(n)

print(f"  Levels grouped by |N|:")
for norm_val in sorted(norm_groups.keys()):
    levels = norm_groups[norm_val]
    occ_in_group = [n for n in levels if n in occupied_n]
    print(f"  |N| = {norm_val:>3d}: levels = {levels}")
    print(f"         occupied: {occ_in_group} ({len(occ_in_group)}/{len(levels)})")
    print()

# ============================================================
# SECTION 6: Why {0, 1, 5, 19}? Algebraic structure
# ============================================================

print("-" * 80)
print("SECTION 6: Algebraic structure of {0, 1, 5, 19}")
print("-" * 80)
print()

# Are these related by phi operations?
# phi has norm -1. So multiplying by phi changes sign of norm.
# phi^2 = phi + 1 has norm 1 (since N(phi^2) = N(phi)^2 = 1)

# Key: N(a+b*phi) = a^2 + ab - b^2
# The discriminant of x^2 + x - 1 = 0 (min poly of phi) is 5.

# The numbers 0, 1, 5, 19:
# 0: trivial
# 1 = 1^2
# 5 = discriminant of Z[phi]
# 19 = ?

# Check: 19 = 4*5 - 1 = 4*disc - 1
print(f"  0 = trivial")
print(f"  1 = unit norm")
print(f"  5 = discriminant of Z[phi]")
print(f"  19 = 4 * 5 - 1 = 4 * disc - 1")
print()

# But also: 19 is a Lucas number? Let's check
lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47]
print(f"  Lucas numbers: {lucas}")
print(f"  19 is Lucas? {19 in lucas}")
print(f"  Fibonacci: 0,1,1,2,3,5,8,13,21...")
print(f"  5 is Fibonacci F_5: YES")
print()

# Relationship between 1, 5, 19:
# 1, 5, 19 differ by 4 and 14
# 5 = 1 + 4, 19 = 5 + 14
# Or: 5 = 5*1, 19 = 5*4 - 1

# Check: are these all representable as norms of SMALL elements?
# N = 1: a+b*phi with a,b in {-1,0,1}
# N = 5: a+b*phi with a,b in {-2,...,2}
# N = 19: a+b*phi with a,b in {-4,...,4}

# Better: what is the SMALLEST (a,b) (by |a|+|b|) achieving each norm?
for target in [0, 1, 5, 19]:
    min_elem = None
    min_c = float('inf')
    for a in range(-20, 21):
        for b in range(-20, 21):
            if abs(zphi_norm(a, b)) == target:
                c = abs(a) + abs(b)
                if c < min_c:
                    min_c = c
                    min_elem = (a, b)
    print(f"  |N| = {target}: smallest (a,b) = {min_elem}, complexity = {min_c}")

print()

# The occupied elements and their norms:
print("  Occupied fermions sorted by |norm|:")
for name, (a, b, n) in sorted(fermions.items(), key=lambda x: abs(zphi_norm(x[1][0], x[1][1]))):
    N = zphi_norm(a, b)
    print(f"  {name:>4s}: ({a:>2d},{b:>2d}), n={n:>2d}, N={N:>3d}, |N|={abs(N):>2d}")

print()

# ============================================================
# SECTION 7: The "unit fraction" pattern
# ============================================================

print("-" * 80)
print("SECTION 7: Units, associates, and prime factorization")
print("-" * 80)
print()

# In Z[phi], the group of units is {+-phi^k : k in Z}
# Two elements are ASSOCIATES if they differ by a unit.
# The factorization of our fermion elements:

print("  Fermion elements a + b*phi and their algebraic properties:")
print()

for name, (a, b, n) in sorted(fermions.items(), key=lambda x: x[1][2]):
    val = a + b * PHI
    val_conj = a + b * PHI_CONJ
    N = zphi_norm(a, b)

    # Is this a unit? (|N| = 1)
    is_unit = abs(N) == 1
    # If unit, find which power of phi
    if is_unit and abs(val) > 1e-10:
        k = np.log(abs(val)) / np.log(PHI)
        k_round = round(k)
        sign = 1 if val > 0 else -1
        check = sign * PHI**k_round
        if abs(check - val) < 1e-6:
            unit_str = f"{'+'if sign>0 else '-'}phi^{k_round}"
        else:
            unit_str = f"~phi^{k:.2f}"
    elif abs(val) < 1e-10:
        unit_str = "0"
    else:
        unit_str = f"norm={N}"

    print(f"  {name:>4s}: {a}+{b}*phi = {val:>8.4f}, conj = {val_conj:>8.4f}, "
          f"N = {N:>4d}, {unit_str}")

print()

# ============================================================
# SECTION 8: Connection between norm and mass hierarchy
# ============================================================

print("-" * 80)
print("SECTION 8: Norm vs mass hierarchy")
print("-" * 80)
print()

# The masses go as m_f = m_e * phi^n
# The norm goes as |a^2 + ab - b^2| where 5a+6b = n
# Is there a relationship between norm and generation?

print("  Fermion norms organized by generation:")
print()
print("  Gen 1 (lightest): e(N=0), u(N=-1), d(N=1)")
print("  Gen 2 (middle):   mu(N=1), c(N=5), s(N=1)")
print("  Gen 3 (heaviest): tau(N=-1), t(N=19), b(N=-19)")
print()

print("  PATTERN: |N| grows with generation!")
print("  Gen 1: max|N| = 1")
print("  Gen 2: max|N| = 5 = discriminant")
print("  Gen 3: max|N| = 19")
print()

# Gen 1 norms: {0, 1, 1} -> max = 1
# Gen 2 norms: {1, 5, 1} -> max = 5
# Gen 3 norms: {1, 19, 19} -> max = 19
# So: 1 -> 5 -> 19. Ratio: 5, 3.8. Hmm.

# Sum of |N| per generation:
gen1 = {'e': (0,0), 'u': (3,-2), 'd': (1,0)}
gen2 = {'mu': (1,1), 'c': (2,1), 's': (1,1)}
gen3 = {'tau': (1,2), 't': (4,1), 'b_q': (-1,4)}

for gen_name, gen_data in [('Gen 1', gen1), ('Gen 2', gen2), ('Gen 3', gen3)]:
    norms = [abs(zphi_norm(a, b)) for a, b in gen_data.values()]
    print(f"  {gen_name}: |norms| = {norms}, sum = {sum(norms)}, max = {max(norms)}")

print()

# ============================================================
# SECTION 9: Exhaustive search - can norms + average select levels?
# ============================================================

print("-" * 80)
print("SECTION 9: Combined constraint: avg=12 AND norm constraint")
print("-" * 80)
print()

# Find ALL 9-element selections from {0,...,26} (with repetition for mu/s)
# such that:
# 1. Sum = 108 (avg = 12)
# 2. Each element has |N| in {0, 1, 5, 19}

# First, which levels have |N| in {0, 1, 5, 19}?
allowed = [n for n in range(27) if abs(all_norms[n]) in {0, 1, 5, 19}]
print(f"  Levels with |N| in {{0, 1, 5, 19}}: {allowed}")
print(f"  Count: {len(allowed)}")
print()

# Find all 8-element subsets (allowing one duplicate for mu/s) of 'allowed'
# with specific sum constraints
# Actually: we need 9 values (not necessarily distinct) from 'allowed'
# that sum to 108, with the sector structure:
# - 3 leptons (one must be 0 = electron)
# - 3 up-type
# - 3 down-type
# - mu/s degeneracy means 1 lepton = 1 down-type

# Simplification: just find all 8-element subsets of 'allowed'
# that could give sum=108 when one element is counted twice.

# Or even simpler: find 9-tuples from 'allowed' (with possible repeats)
# summing to 108.

# But this is huge. Let's just count how many 3-sector configurations work.

print("  Searching for valid 3-sector configurations...")
print("  Constraint: 3 leptons (one=0) + 3 up + 3 down = 108")
print("  All levels from the 'allowed' set with |N| in {0,1,5,19}")
print()

# Leptons: e=0, then 2 more from allowed, ordered
lep_options = []
for l2 in allowed:
    for l3 in allowed:
        if l3 > l2 > 0:
            lep_options.append((0, l2, l3))

# Up-type: 3 from allowed, ordered
up_options = []
for u1 in allowed:
    for u2 in allowed:
        for u3 in allowed:
            if u3 > u2 > u1 > 0:
                up_options.append((u1, u2, u3))

# Down-type: 3 from allowed, ordered
down_options = []
for d1 in allowed:
    for d2 in allowed:
        for d3 in allowed:
            if d3 > d2 > d1 > 0:
                down_options.append((d1, d2, d3))

print(f"  Lepton options: {len(lep_options)}")
print(f"  Up-type options: {len(up_options)}")
print(f"  Down-type options: {len(down_options)}")
print()

# Find all combos summing to 108
valid_configs = []
for lep in lep_options:
    lep_sum = sum(lep)
    for up in up_options:
        up_sum = sum(up)
        remaining = 108 - lep_sum - up_sum
        # remaining must equal sum of a down triple from allowed
        for down in down_options:
            if sum(down) == remaining:
                valid_configs.append((lep, up, down))

print(f"  Valid configurations (sum=108, all from allowed): {len(valid_configs)}")
print()

if len(valid_configs) <= 50:
    for i, (lep, up, down) in enumerate(valid_configs[:20]):
        all_levels = sorted(lep + up + down)
        print(f"  Config {i+1}: L={lep}, U={up}, D={down}")
        print(f"            All levels: {all_levels}")

    # Is our actual configuration among them?
    our_config = ((0, 11, 17), (3, 16, 26), (5, 11, 19))
    if our_config in valid_configs:
        print(f"\n  *** OUR CONFIGURATION IS IN THE LIST! ***")
        print(f"  Position: {valid_configs.index(our_config) + 1}/{len(valid_configs)}")
    else:
        # Try permutations of sectors
        found = False
        from itertools import permutations
        for perm in permutations([our_config[0], our_config[1], our_config[2]]):
            if perm in valid_configs:
                found = True
                print(f"\n  *** OUR CONFIGURATION FOUND (permuted sectors)! ***")
                break
        if not found:
            print(f"\n  Our configuration NOT found directly.")
            print(f"  Our config: L=(0,11,17), U=(3,16,26), D=(5,11,19)")
            # Check if the levels are in allowed
            for n in [0, 3, 5, 11, 16, 17, 19, 26]:
                print(f"    n={n}: in allowed? {n in allowed}, |N|={abs(all_norms[n])}")
else:
    print(f"  (too many to list)")
    # Check if ours is there
    our_config = ((0, 11, 17), (3, 16, 26), (5, 11, 19))
    if our_config in valid_configs:
        print(f"\n  *** OUR CONFIGURATION IS CONFIG #{valid_configs.index(our_config)+1} ***")

print()

# ============================================================
# SECTION 10: Additional constraint - mu/s degeneracy
# ============================================================

print("-" * 80)
print("SECTION 10: Adding mu/s degeneracy constraint")
print("-" * 80)
print()

# Extra constraint: one level must be shared between leptons and down-type
# (the mu/s degeneracy at n=11)

degen_configs = []
for lep, up, down in valid_configs:
    shared = set(lep) & set(down)
    if shared:
        degen_configs.append((lep, up, down, shared))

print(f"  Configurations with lepton-down sharing: {len(degen_configs)}")
print()

for i, (lep, up, down, shared) in enumerate(degen_configs[:20]):
    all_levels = sorted(set(lep + up + down))
    print(f"  Config {i+1}: L={lep}, U={up}, D={down}, shared={shared}")

print()

# Even more restrictive: electron must be at n=0
# (already enforced in our search)

# ============================================================
# SECTION 11: The norm sequence and Fibonacci
# ============================================================

print("-" * 80)
print("SECTION 11: Norm sequence and algebraic patterns")
print("-" * 80)
print()

# Norms for all levels: compute and look for patterns
print("  Norm sequence N(n) for n = 0 to 26:")
norm_seq = [all_norms[n] for n in range(27)]
abs_norm_seq = [abs(all_norms[n]) for n in range(27)]
print(f"  N:     {norm_seq}")
print(f"  |N|:   {abs_norm_seq}")
print()

# The absolute norms: 0,1,4,1,1,1,1,5,11,5,4,1,4,11,11,9,5,1,9,19,16,11,4,5,16,25,19
# Unique absolute norms: 0,1,4,5,9,11,16,19,25
# These are: 0, 1, 4, 5, 9, 11, 16, 19, 25

unique_abs = sorted(set(abs_norm_seq))
print(f"  Unique |N| values: {unique_abs}")
print()

# Are these perfect squares?
# 0=0^2, 1=1^2, 4=2^2, 9=3^2, 16=4^2, 25=5^2 -> YES for {0,1,4,9,16,25}
# But 5, 11, 19 are NOT perfect squares.
# 5 = disc, 11 = Lucas_5, 19 = ?

# Split into squares and non-squares:
squares = [n for n in unique_abs if int(np.sqrt(n))**2 == n]
non_squares = [n for n in unique_abs if int(np.sqrt(n))**2 != n]
print(f"  Perfect squares: {squares}")
print(f"  Non-squares: {non_squares}")
print()

# The non-squares are {5, 11, 19}
# 5 mod 5 = 0 (ramified)
# 11 mod 5 = 1 (splits)
# 19 mod 5 = 4 = -1 (splits)
print(f"  Non-square norms: {{5, 11, 19}}")
print(f"  5 mod 5 = 0 (RAMIFIED)")
print(f"  11 mod 5 = 1 (SPLITS)")
print(f"  19 mod 5 = 4 (SPLITS)")
print()

# Interesting: 5, 11, 19 are all prime!
# And: 5 = L_3, 11 = L_5, 19 = ?
# Lucas: 2, 1, 3, 4, 7, 11, 18, 29
# 5 is NOT a Lucas number but IS L_3 = 4? No, L_3 = 4.
# Actually: L_0=2, L_1=1, L_2=3, L_3=4, L_4=7, L_5=11, L_6=18, L_7=29
# So 11 = L_5. 5 is F_5. 19 is neither Fibonacci nor Lucas.

# But: 5, 11, 19 differ by 6 and 8.
# 5 + 6 = 11, 11 + 8 = 19.
# 6 = b_1 (intersection number), 8 = ??

# Or: 5 = a_1 (intersection number!), 11 = a_1 + b_1, 19 = a_1 + 2*b_1 + 2
# 5, 5+6=11, 11+8=19. Not clean.

# Or: N for occupied levels: {0, 1, 5, 19}
# Sum: 0 + 1 + 5 + 19 = 25 = 5^2
print(f"  Sum of unique occupied |N|: {sum(occ_abs_norms)} = {sum(occ_abs_norms)}")
print(f"  5^2 = 25: {sum(occ_abs_norms) == 25}")
print()

# Product: 1 * 5 * 19 = 95
print(f"  Product (non-zero): {1*5*19} = 95")
print(f"  95 = 5 * 19")
print()

# ============================================================
# SECTION 12: The norm as "generation number"
# ============================================================

print("-" * 80)
print("SECTION 12: Norm as generation discriminant")
print("-" * 80)
print()

# Observation from Section 8: |N| grows with generation
# Gen 1: max|N| = 1 (all units)
# Gen 2: max|N| = 5 (charm has |N|=5)
# Gen 3: max|N| = 19 (top and bottom have |N|=19)

# This means the NORM in Z[phi] acts as a generation label!
# The generation structure IS the norm structure.

print("  Generation structure from Z[phi] norms:")
print()
print("  Generation 1 (|N| <= 1): e(0), u(1), d(1)")
print("  Generation 2 (|N| <= 5): mu(1), c(5), s(1)")
print("  Generation 3 (|N| <= 19): tau(1), t(19), b(19)")
print()
print("  The THRESHOLD norms per generation: 1, 5, 19")
print("  These are: unit, ramified prime, split prime")
print()

# Note: mu, s, tau all have |N|=1 (they're units in Z[phi])
# Only charm (gen 2) and top/bottom (gen 3) have non-unit norms.
# The norm discriminates QUARKS better than leptons.

# Summary: the algebraic number theory of Z[phi] provides
# a natural generation structure through the norm form.

# ============================================================
# SECTION 13: Summary
# ============================================================

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

print("  1. Average occupied level = 108/9 = 12 = vertex degree. EXACT.")
print(f"     Probability by random chance: ~{count_108/n_trials*100:.1f}%")
print()
print("  2. Z[phi] norms of occupied fermions: |N| in {{0, 1, 5, 19}}")
print(f"     Sum: 0+1+5+19 = 25 = 5^2")
print()
print("  3. Norm as generation discriminant:")
print("     Gen 1: all |N| <= 1 (units)")
print("     Gen 2: |N| up to 5 (ramified prime)")
print("     Gen 3: |N| up to 19 (split prime)")
print()
print(f"  4. Combined constraint (avg=12 + |N| in {{0,1,5,19}})")
print(f"     selects {len(valid_configs)} configurations from the full space.")
if degen_configs:
    print(f"     Adding mu/s degeneracy narrows to {len(degen_configs)}.")
print()
print("  5. The non-square norms {{5, 11, 19}} appearing in 0-26 are ALL PRIME.")
print("     5 = ramified, 11 = splits, 19 = splits in Z[phi].")
print()
print("  CLASSIFICATION:")
print("  - Average = 12: PATTERN (exact but may be coincidence)")
print("  - Norm structure {0,1,5,19}: PATTERN (algebraically natural)")
print("  - Generation from norms: PATTERN (interesting but not derivation)")
print("  - Combined selection power: needs further analysis")
print()

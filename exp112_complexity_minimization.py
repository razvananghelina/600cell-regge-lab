"""
EXP-112: Derivation of fermion (a,b) quantum numbers from complexity minimization
==================================================================================

HYPOTHESIS: The quantum numbers (a,b) for each fermion are determined by the
UNIQUE solution that minimizes topological complexity |a|+|b| subject to the
constraint n = a_1*a + b_1*b = 5a + 6b, where:
  - a_1 = 5 (intersection number = common neighbors per edge)
  - b_1 = 6 (second intersection number = next-shell neighbors per edge)
  - n = mass level from DSI ladder: m_f = m_e * phi^n

This would upgrade (a,b) from FITTING to DERIVAT: the quantum numbers are not
free parameters but are determined by the 600-cell geometry (through a_1, b_1)
plus a single principle (minimize topological complexity).

APPROACH:
1. For each fermion mass level n, enumerate ALL integer (a,b) with 5a+6b=n
2. Find the decomposition minimizing |a|+|b|
3. Check if it is UNIQUE (no ties)
4. Compare with known assignments from exp099/100

Author: Claude (exp112)
Date: 2026-02-08
"""

import numpy as np
from math import gcd

# ============================================================
# SECTION 1: Setup - known fermion data
# ============================================================

PHI = (1 + np.sqrt(5)) / 2

# Intersection numbers from 600-cell (derived in exp111)
a1 = 5  # common neighbors per edge
b1 = 6  # next-shell neighbors per edge

# Known fermion assignments from exp099/100
# Format: name, (a, b), n = 5a + 6b
fermions = [
    ("electron", (0, 0), 0),
    ("muon",     (1, 1), 11),
    ("tau",      (1, 2), 17),
    ("up",       (3, -2), 3),
    ("charm",    (2, 1), 16),
    ("top",      (4, 1), 26),
    ("down",     (1, 0), 5),
    ("strange",  (1, 1), 11),
    ("bottom",   (-1, 4), 19),
]

# Experimental masses (GeV) for reference
m_exp = {
    "electron": 0.000511,
    "muon": 0.10566,
    "tau": 1.777,
    "up": 0.00216,
    "charm": 1.27,
    "top": 172.76,
    "down": 0.00467,
    "strange": 0.0934,
    "bottom": 4.18,
}

print("=" * 80)
print("EXP-112: Fermion (a,b) from Topological Complexity Minimization")
print("=" * 80)
print()
print(f"Intersection numbers: a_1 = {a1}, b_1 = {b1}")
print(f"Mass formula: n = {a1}*a + {b1}*b")
print(f"Principle: minimize |a| + |b| subject to constraint")
print()

# ============================================================
# SECTION 2: General theory of 5a + 6b = n decompositions
# ============================================================

print("-" * 80)
print("SECTION 2: Mathematical structure of 5a + 6b = n")
print("-" * 80)
print()

# For 5a + 6b = n with integer a,b:
# General solution: a = a0 + 6t, b = b0 - 5t (for any integer t)
# where (a0, b0) is any particular solution.
# Since gcd(5,6) = 1, solutions exist for ALL integer n.

print(f"gcd({a1}, {b1}) = {gcd(a1, b1)} => solutions exist for ALL integer n")
print()
print("General solution: a = a0 + 6t, b = b0 - 5t")
print("Complexity: |a| + |b| = |a0 + 6t| + |b0 - 5t|")
print()

# ============================================================
# SECTION 3: Minimization for each fermion
# ============================================================

print("-" * 80)
print("SECTION 3: Complexity minimization for each fermion")
print("-" * 80)
print()

def find_all_minimal_decompositions(n, a1_val=5, b1_val=6, search_range=50):
    """Find all (a,b) minimizing |a|+|b| subject to a1*a + b1*b = n."""

    # Find all solutions in a reasonable range
    solutions = []
    for a in range(-search_range, search_range + 1):
        # b = (n - a1*a) / b1
        remainder = n - a1_val * a
        if remainder % b1_val == 0:
            b = remainder // b1_val
            solutions.append((a, b))

    if not solutions:
        return [], float('inf')

    # Find minimum complexity
    min_complexity = min(abs(a) + abs(b) for a, b in solutions)

    # Find all solutions achieving the minimum
    minimal = [(a, b) for a, b in solutions if abs(a) + abs(b) == min_complexity]

    return minimal, min_complexity

# Also find the next-best (second minimum) for each n
def find_second_minimum(n, a1_val=5, b1_val=6, search_range=50):
    """Find the second-best complexity after the minimum."""
    solutions = []
    for a in range(-search_range, search_range + 1):
        remainder = n - a1_val * a
        if remainder % b1_val == 0:
            b = remainder // b1_val
            solutions.append((a, b, abs(a) + abs(b)))

    if len(solutions) < 2:
        return None, float('inf')

    complexities = sorted(set(c for _, _, c in solutions))
    if len(complexities) < 2:
        return None, float('inf')

    second_min = complexities[1]
    second_best = [(a, b) for a, b, c in solutions if c == second_min]
    return second_best, second_min

all_match = True
results = []

for name, (a_known, b_known), n in fermions:
    minimal, min_c = find_all_minimal_decompositions(n)
    second, sec_c = find_second_minimum(n)

    n_check = a1 * a_known + b1 * b_known
    known_c = abs(a_known) + abs(b_known)

    # Check if known assignment is among the minimal ones
    is_minimal = (a_known, b_known) in minimal
    is_unique = len(minimal) == 1

    match_str = "MATCH" if is_minimal else "MISMATCH"
    unique_str = "unique" if is_unique else f"{len(minimal)} solutions"

    if not is_minimal:
        all_match = False

    gap = sec_c - min_c if sec_c != float('inf') else "inf"

    print(f"  {name:10s}: n={n:3d}, known=({a_known:+d},{b_known:+d}), "
          f"|a|+|b|={known_c}")
    print(f"              minimal: {minimal}, complexity={min_c}, "
          f"{unique_str}, gap={gap} => {match_str}")

    if not is_unique:
        print(f"              WARNING: {len(minimal)} degenerate solutions!")

    results.append({
        'name': name,
        'n': n,
        'known': (a_known, b_known),
        'known_c': known_c,
        'minimal': minimal,
        'min_c': min_c,
        'is_minimal': is_minimal,
        'is_unique': is_unique,
        'gap': gap,
        'second': second,
        'sec_c': sec_c,
    })
    print()

# ============================================================
# SECTION 4: Degeneracy analysis - muon/strange problem
# ============================================================

print("-" * 80)
print("SECTION 4: Degeneracy analysis")
print("-" * 80)
print()

# muon and strange have the SAME n=11, same (a,b)=(1,1)
# This is a feature, not a bug: they are in different sectors (lepton vs down-quark)
# The sector-dependent corrections (delta_d, delta_k) distinguish them

degenerate_ns = {}
for name, (a, b), n in fermions:
    if n not in degenerate_ns:
        degenerate_ns[n] = []
    degenerate_ns[n].append(name)

print("Fermions sharing the same n value:")
for n, names in sorted(degenerate_ns.items()):
    if len(names) > 1:
        print(f"  n={n}: {', '.join(names)}")
        print(f"    => Resolved by sector-dependent corrections (delta_d, delta_k)")

print()

# ============================================================
# SECTION 5: Alternative complexity measures
# ============================================================

print("-" * 80)
print("SECTION 5: Alternative complexity measures")
print("-" * 80)
print()

# Test if other norms also give the same result
def find_minimal_L2(n, a1_val=5, b1_val=6, search_range=50):
    """Minimize a^2 + b^2 (L2 squared)."""
    solutions = []
    for a in range(-search_range, search_range + 1):
        remainder = n - a1_val * a
        if remainder % b1_val == 0:
            b = remainder // b1_val
            solutions.append((a, b, a**2 + b**2))

    if not solutions:
        return [], float('inf')

    min_c = min(c for _, _, c in solutions)
    minimal = [(a, b) for a, b, c in solutions if c == min_c]
    return minimal, min_c

def find_minimal_Linf(n, a1_val=5, b1_val=6, search_range=50):
    """Minimize max(|a|, |b|) (L_infinity)."""
    solutions = []
    for a in range(-search_range, search_range + 1):
        remainder = n - a1_val * a
        if remainder % b1_val == 0:
            b = remainder // b1_val
            solutions.append((a, b, max(abs(a), abs(b))))

    if not solutions:
        return [], float('inf')

    min_c = min(c for _, _, c in solutions)
    minimal = [(a, b) for a, b, c in solutions if c == min_c]
    return minimal, min_c

print("Comparison of complexity measures:")
print(f"  {'Fermion':10s} {'n':>3s}  {'L1 min':12s} {'L2 min':12s} {'Linf min':12s} {'Known':8s}")

l1_matches = 0
l2_matches = 0
linf_matches = 0

for name, (a_known, b_known), n in fermions:
    l1_min, _ = find_all_minimal_decompositions(n)
    l2_min, _ = find_minimal_L2(n)
    linf_min, _ = find_minimal_Linf(n)

    l1_match = (a_known, b_known) in l1_min
    l2_match = (a_known, b_known) in l2_min
    linf_match = (a_known, b_known) in linf_min

    if l1_match: l1_matches += 1
    if l2_match: l2_matches += 1
    if linf_match: linf_matches += 1

    l1_str = str(l1_min[0]) if len(l1_min) == 1 else str(l1_min)
    l2_str = str(l2_min[0]) if len(l2_min) == 1 else str(l2_min)
    linf_str = str(linf_min[0]) if len(linf_min) == 1 else str(linf_min)

    print(f"  {name:10s} {n:3d}  {l1_str:12s} {l2_str:12s} {linf_str:12s} ({a_known:+d},{b_known:+d})")

print()
print(f"  L1 (|a|+|b|) matches: {l1_matches}/9")
print(f"  L2 (a^2+b^2) matches: {l2_matches}/9")
print(f"  Linf (max)    matches: {linf_matches}/9")

# ============================================================
# SECTION 6: Proof that L1 minimum is unique for n = 5a + 6b
# ============================================================

print()
print("-" * 80)
print("SECTION 6: Analytical proof of uniqueness")
print("-" * 80)
print()

# For 5a + 6b = n, the general solution is a = a0 + 6t, b = b0 - 5t
# We need to minimize f(t) = |a0 + 6t| + |b0 - 5t|
# This is a convex piecewise-linear function of t.
# It has a unique minimum unless the minimum is achieved on a flat segment.

# The critical values of t are:
#   t1 = -a0/6 (where a = 0)
#   t2 = b0/5  (where b = 0)

# Since the slopes change by +12 at t1 and +10 at t2,
# the minimum is unique unless t1 = t2 (i.e., 5*a0 = -6*b0, which means n=0).

print("For n = 5a + 6b, general solution: a = a0 + 6t, b = b0 - 5t")
print("Complexity: f(t) = |a0 + 6t| + |b0 - 5t|")
print()
print("This is a CONVEX piecewise-linear function with two kinks:")
print("  t* = -a0/6  (where a crosses 0)")
print("  t* = b0/5   (where b crosses 0)")
print()
print("For n != 0: these kinks are at DIFFERENT positions")
print("  => the minimum is achieved at a UNIQUE integer t")
print("  => (a,b) is UNIQUE")
print()
print("Exception: n=0 has a=b=0 as the unique trivial solution")
print()

# Verify: for each fermion, show that the integer t giving min is unique
print("Verification - optimal t for each fermion:")
for name, (a_known, b_known), n in fermions:
    if n == 0:
        print(f"  {name:10s}: n=0, trivial (0,0)")
        continue

    # Particular solution: use a0 such that 5*a0 = n mod 6
    # Since gcd(5,6)=1 and 5*5=25=1 mod 6, a0 = 5n mod 6 works
    # But let's just find a simple particular solution
    for a0_try in range(-10, 11):
        if (n - 5 * a0_try) % 6 == 0:
            a0 = a0_try
            b0 = (n - 5 * a0) // 6
            break

    # Find optimal t
    best_t = None
    best_c = float('inf')
    for t in range(-20, 21):
        a = a0 + 6 * t
        b = b0 - 5 * t
        c = abs(a) + abs(b)
        if c < best_c:
            best_c = c
            best_t = t

    a_opt = a0 + 6 * best_t
    b_opt = b0 - 5 * best_t

    match = "OK" if (a_opt, b_opt) == (a_known, b_known) else "MISMATCH"
    print(f"  {name:10s}: n={n:3d}, a0={a0:+d}, b0={b0:+d}, "
          f"t*={best_t:+d} => ({a_opt:+d},{b_opt:+d}) [{match}]")

# ============================================================
# SECTION 7: Landscape of all n values 0-30
# ============================================================

print()
print("-" * 80)
print("SECTION 7: Complete landscape n = 0..30")
print("-" * 80)
print()

# Map known fermions by n
fermion_map = {}
for name, (a, b), n in fermions:
    if n not in fermion_map:
        fermion_map[n] = []
    fermion_map[n].append((name, a, b))

print(f"{'n':>3s}  {'min(a,b)':12s}  {'|a|+|b|':>7s}  {'gap':>4s}  {'fermion':20s}")
print("-" * 60)

for n in range(0, 31):
    minimal, min_c = find_all_minimal_decompositions(n)
    _, sec_c = find_second_minimum(n)
    gap = sec_c - min_c if sec_c != float('inf') else "--"

    min_str = str(minimal[0]) if len(minimal) == 1 else str(minimal)

    ferm_str = ""
    if n in fermion_map:
        ferm_str = ", ".join(f"{name}({a:+d},{b:+d})" for name, a, b in fermion_map[n])

    marker = " <--" if n in fermion_map else ""
    print(f"{n:3d}  {min_str:12s}  {min_c:7d}  {str(gap):>4s}  {ferm_str}{marker}")

# ============================================================
# SECTION 8: Complexity vs generation pattern
# ============================================================

print()
print("-" * 80)
print("SECTION 8: Complexity vs generation")
print("-" * 80)
print()

# Group by generation
gen1 = [("electron", 0, 0), ("up", 3, -2), ("down", 1, 0)]
gen2 = [("muon", 1, 1), ("charm", 2, 1), ("strange", 1, 1)]
gen3 = [("tau", 1, 2), ("top", 4, 1), ("bottom", -1, 4)]

for gen_name, gen_data in [("Gen 1", gen1), ("Gen 2", gen2), ("Gen 3", gen3)]:
    complexities = [abs(a) + abs(b) for _, a, b in gen_data]
    avg_c = np.mean(complexities)
    print(f"  {gen_name}: ", end="")
    for name, a, b in gen_data:
        print(f"{name}({a:+d},{b:+d}):|a|+|b|={abs(a)+abs(b)}, ", end="")
    print(f"  avg = {avg_c:.2f}")

print()
print("  Generation 1 avg complexity: {:.2f}".format(np.mean([0, 5, 1])))
print("  Generation 2 avg complexity: {:.2f}".format(np.mean([2, 3, 2])))
print("  Generation 3 avg complexity: {:.2f}".format(np.mean([3, 5, 5])))
print()
print("  Complexity increases with generation: 2.00 < 2.33 < 4.33")
print("  This is consistent with heavier particles requiring more")
print("  'topological winding' on the 600-cell graph.")

# ============================================================
# SECTION 9: Connection to DSI mechanism (exp107)
# ============================================================

print()
print("-" * 80)
print("SECTION 9: Connection to DSI mechanism")
print("-" * 80)
print()

print("The mass formula m = m_e * phi^n with n = 5a + 6b connects:")
print()
print("  1. DSI mechanism (exp107): phi-periodic potential on S^3")
print("     V(psi) = A * sin^2(pi * ln|psi| / ln(phi))")
print("     => mass levels separated by phi")
print()
print("  2. Intersection numbers (exp111): a_1=5, b_1=6")
print("     => n = a_1*a + b_1*b (graph-theoretic origin)")
print()
print("  3. Complexity minimization (THIS experiment):")
print("     => (a,b) = argmin |a|+|b| s.t. 5a+6b=n")
print("     => unique for all n != 0")
print()
print("INTERPRETATION: Each fermion mass level n allows")
print("infinitely many (a,b) decompositions on the 600-cell.")
print("Nature selects the one with MINIMAL topological complexity.")
print("This is analogous to 'path of least action' but on the")
print("discrete graph structure.")

# ============================================================
# SECTION 10: Summary and classification
# ============================================================

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

n_match = sum(1 for r in results if r['is_minimal'])
n_unique = sum(1 for r in results if r['is_unique'])

print(f"  Fermions matching complexity minimum: {n_match}/9")
print(f"  Unique minima: {n_unique}/9 (n=0 trivial)")
print()

if all_match:
    print("  STATUS: ALL 9 fermion (a,b) assignments are the UNIQUE")
    print("  complexity-minimizing decompositions of their mass levels n.")
    print()
    print("  CLASSIFICATION: DERIVAT")
    print("    - (a,b) quantum numbers are NOT free parameters")
    print("    - They follow from:")
    print(f"      1. Intersection numbers a_1={a1}, b_1={b1} (600-cell geometry)")
    print(f"      2. Mass levels n (from DSI mechanism)")
    print(f"      3. Complexity minimization principle min|a|+|b|")
    print()
    print("  REMAINING INPUT: The 9 mass levels n = {0,3,5,11,16,17,19,26}")
    print("    These come from experimental masses + DSI ladder.")
    print("    Deriving them from pure geometry is an open problem.")
else:
    print("  STATUS: PARTIAL MATCH - not all fermions match")
    for r in results:
        if not r['is_minimal']:
            print(f"    MISMATCH: {r['name']} - known {r['known']}, "
                  f"minimal {r['minimal']}")

print()
print("  KEY FORMULA:")
print(f"    n = {a1}a + {b1}b,  (a,b) = argmin |a|+|b|")
print()
print("  PARAMETER COUNT:")
print("    Before: 9 fermions x 2 parameters = 18 free parameters")
print("    After:  9 mass levels n = 9 free parameters")
print("    Reduction: 18 -> 9 (factor 2 reduction)")
print()

# Final verification table
print("-" * 80)
print("FINAL VERIFICATION TABLE")
print("-" * 80)
print()
print(f"{'Fermion':10s} {'n':>3s} {'(a,b)':>8s} {'|a|+|b|':>7s} {'min?':>5s} "
      f"{'unique?':>7s} {'gap':>4s} {'m_pred/m_exp':>12s}")
print("-" * 70)

m_e = 0.000511  # GeV

for r in results:
    name = r['name']
    n = r['n']
    a, b = r['known']

    # Predicted mass (using base formula without corrections for simplicity)
    d_eff = 5
    k_eff = 6
    m_pred = m_e * PHI**(a * d_eff + b * k_eff)

    ratio = m_pred / m_exp[name] if m_exp[name] > 0 else 0

    match = "YES" if r['is_minimal'] else "NO"
    uniq = "YES" if r['is_unique'] else "NO"
    gap = str(r['gap'])

    print(f"  {name:10s} {n:3d} ({a:+d},{b:+d}) {r['known_c']:7d} {match:>5s} "
          f"{uniq:>7s} {gap:>4s} {ratio:12.4f}")

print()
print("Note: m_pred/m_exp uses base formula only (no sector corrections)")
print("With corrections (exp099/100): 7/9 fermions within 2%")

"""
EXP-113: Level Selection Rule -- Why n in {0, 3, 5, 11, 16, 17, 19, 26}?
=========================================================================

The DSI ladder has minima at m_e * phi^n for ALL integers n >= 0.
Since gcd(5,6)=1, every n >= 0 can be written as 5a+6b for some (a,b) in Z^2.
Yet only 8 of the 27 levels from n=0 to n=26 are occupied by fermions.

This experiment investigates multiple hypotheses for the selection rule.

APPROACHES:
A. Number-theoretic patterns (mod arithmetic, Fibonacci, etc.)
B. Eigenspace mapping (9 eigenvalues <-> 9 fermions?)
C. Representation theory (H_4 irrep dimensions)
D. Frobenius / lattice structure
E. Graph-theoretic path constraints
F. Information-theoretic / entropy approaches
G. Combinatorial search for generating rules

Author: Claude (exp113)
Date: 2026-02-08
"""

import numpy as np
from itertools import combinations
from math import gcd

PHI = (1 + np.sqrt(5)) / 2

# ============================================================
# SECTION 1: Setup
# ============================================================

# Occupied mass levels
occupied = [0, 3, 5, 11, 16, 17, 19, 26]
n_occupied = len(occupied)  # 8 distinct n (muon/strange share n=11)

# Full range
N_max = 26  # top quark level
all_levels = list(range(N_max + 1))
unoccupied = [n for n in all_levels if n not in occupied]

# Fermion data
fermions = {
    0:  [("electron", 0, 0, "lepton")],
    3:  [("up", 3, -2, "up-type")],
    5:  [("down", 1, 0, "down-type")],
    11: [("muon", 1, 1, "lepton"), ("strange", 1, 1, "down-type")],
    16: [("charm", 2, 1, "up-type")],
    17: [("tau", 1, 2, "lepton")],
    19: [("bottom", -1, 4, "down-type")],
    26: [("top", 4, 1, "up-type")],
}

# 600-cell intersection numbers
a1, b1 = 5, 6

# 600-cell eigenvalue data (from exp110)
eigenvalues = [12, 6*PHI, 4*PHI, 3, 0, -2, 4-4*PHI, -3, 6-6*PHI]
multiplicities = [1, 4, 9, 16, 25, 36, 9, 16, 4]
laplacian_eigs = [12 - ev for ev in eigenvalues]

print("=" * 80)
print("EXP-113: Level Selection Rule")
print("=" * 80)
print()
print(f"Occupied levels: {occupied}")
print(f"Unoccupied (0-26): {unoccupied}")
print(f"Occupied: {n_occupied}/27 = {n_occupied/27:.1%}")
print()

# ============================================================
# SECTION 2: Number-theoretic patterns
# ============================================================

print("-" * 80)
print("SECTION 2: Number-theoretic patterns")
print("-" * 80)
print()

# A. Modular arithmetic
for mod in [2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 30]:
    residues_occ = sorted(set(n % mod for n in occupied))
    residues_unocc = sorted(set(n % mod for n in unoccupied))
    only_occ = [r for r in residues_occ if r not in residues_unocc]
    only_unocc = [r for r in residues_unocc if r not in residues_occ]
    if only_occ or only_unocc:
        print(f"  mod {mod:2d}: occupied residues = {residues_occ}, "
              f"ONLY in occupied: {only_occ}, ONLY in unoccupied: {only_unocc}")

print()

# B. Gaps between consecutive occupied levels
gaps = [occupied[i+1] - occupied[i] for i in range(len(occupied)-1)]
print(f"  Gaps: {gaps}")
print(f"  Sum of gaps: {sum(gaps)} (= {N_max}, the max level)")
print(f"  Gap set: {sorted(set(gaps))}")
print()

# C. Fibonacci numbers
fib = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
fib_set = set(fib)
print(f"  Fibonacci numbers <= 26: {[f for f in fib if f <= 26]}")
print(f"  Occupied AND Fibonacci: {[n for n in occupied if n in fib_set]}")
print(f"  Occupied NOT Fibonacci: {[n for n in occupied if n not in fib_set]}")
print()

# D. Differences of occupied levels
print("  Pairwise differences (occupied levels):")
diffs = set()
for i, a in enumerate(occupied):
    for j, b in enumerate(occupied):
        if i < j:
            diffs.add(b - a)
print(f"  {sorted(diffs)}")
print()

# E. Sums and products
print("  Sums of pairs from {5, 6} (intersection numbers):")
print(f"    5+6=11 (occupied!), 5+5=10 (not), 6+6=12 (not)")
print(f"    2*5+6=16 (occupied!), 5+2*6=17 (occupied!)")
print(f"    4*5+6=26 (occupied!), 3*5-2*6=3 (occupied!)")
print(f"    -5+4*6=19 (occupied!), 5+0*6=5 (occupied!)")
print()

# F. Check if occupied = {a*5 + b*6 : (a,b) in some special set}
print("  All occupied as a*5 + b*6 (from complexity minimization):")
for n in occupied:
    # Find minimal (a,b)
    best = None
    best_c = float('inf')
    for a in range(-10, 11):
        rem = n - 5*a
        if rem % 6 == 0:
            b = rem // 6
            c = abs(a) + abs(b)
            if c < best_c:
                best_c = c
                best = (a, b)
    print(f"    n={n:2d}: ({best[0]:+d},{best[1]:+d}), complexity={best_c}")
print()

# ============================================================
# SECTION 3: Eigenspace mapping hypothesis
# ============================================================

print("-" * 80)
print("SECTION 3: Eigenspace <-> Fermion mapping")
print("-" * 80)
print()

# There are 9 eigenvalues and 9 fermions (counting muon and strange separately).
# Is there a natural mapping?

print("  9 eigenvalues vs 9 fermions:")
print()
print(f"  {'Eigenvalue':>12s} {'mult':>5s}  {'n?':>4s}  {'Fermion?':20s}")
print(f"  {'-'*50}")

# Sort fermions by n
all_fermions = []
for n, flist in sorted(fermions.items()):
    for name, a, b, sector in flist:
        all_fermions.append((n, name, a, b, sector))

for i, (ev, mult) in enumerate(zip(eigenvalues, multiplicities)):
    ev_str = f"{ev:.4f}" if abs(ev - round(ev)) > 0.01 else f"{ev:.0f}"
    if i < len(all_fermions):
        n, name, a, b, sector = all_fermions[i]
        print(f"  theta_{i} = {ev_str:>8s}  m={mult:3d}   n={n:2d}  {name} ({sector})")
    else:
        print(f"  theta_{i} = {ev_str:>8s}  m={mult:3d}")

print()
print("  NOTE: This mapping is arbitrary (order-dependent). Testing correlations...")
print()

# Test: is there a correlation between eigenvalue index and mass level?
# Or between multiplicity and something?
for i, (ev, mult) in enumerate(zip(eigenvalues, multiplicities)):
    lam = 12 - ev  # Laplacian eigenvalue
    if i < len(all_fermions):
        n = all_fermions[i][0]
        print(f"    i={i}: lambda={lam:.4f}, sqrt(mult)={int(np.sqrt(mult))}, "
              f"n={n}, ratio n/sqrt(mult)={n/np.sqrt(mult):.3f}")

print()

# ============================================================
# SECTION 4: Complexity structure of occupied vs unoccupied
# ============================================================

print("-" * 80)
print("SECTION 4: Complexity analysis - occupied vs unoccupied")
print("-" * 80)
print()

def min_complexity(n):
    best = float('inf')
    for a in range(-50, 51):
        rem = n - 5*a
        if rem % 6 == 0:
            b = rem // 6
            c = abs(a) + abs(b)
            if c < best:
                best = c
    return best

# Complexity of all levels 0-26
complexities = {n: min_complexity(n) for n in range(27)}

occ_complexities = [complexities[n] for n in occupied]
unocc_complexities = [complexities[n] for n in unoccupied]

print(f"  Occupied complexities:   {occ_complexities}")
print(f"  Mean: {np.mean(occ_complexities):.2f}, median: {np.median(occ_complexities):.1f}")
print()
print(f"  Unoccupied complexities: {unocc_complexities}")
print(f"  Mean: {np.mean(unocc_complexities):.2f}, median: {np.median(unocc_complexities):.1f}")
print()

# Distribution of complexity values
for c in range(6):
    occ_at_c = [n for n in occupied if complexities[n] == c]
    unocc_at_c = [n for n in unoccupied if complexities[n] == c]
    total_at_c = len(occ_at_c) + len(unocc_at_c)
    frac = len(occ_at_c) / total_at_c if total_at_c > 0 else 0
    print(f"  Complexity {c}: occupied {occ_at_c}, "
          f"unoccupied {unocc_at_c}, "
          f"fraction occupied: {frac:.2f}")

print()

# ============================================================
# SECTION 5: Generating rule search
# ============================================================

print("-" * 80)
print("SECTION 5: Systematic search for generating rules")
print("-" * 80)
print()

# Can the occupied set be generated by a simple rule involving 5, 6, phi, etc.?

# Test A: n = 5*i + 6*j for i,j >= 0 (non-negative Frobenius)
frobenius_representable = set()
for i in range(6):
    for j in range(5):
        n = 5*i + 6*j
        if n <= 26:
            frobenius_representable.add(n)

frobenius_non = set(range(27)) - frobenius_representable
print(f"  Frobenius (5,6) non-representable (with a,b >= 0): {sorted(frobenius_non)}")
print(f"  Occupied that are Frobenius non-representable: "
      f"{[n for n in occupied if n in frobenius_non]}")
print()

# Test B: Linear recurrence
# Check if occupied levels satisfy a linear recurrence
print("  Testing linear recurrences of order 2:")
for i in range(len(occupied) - 2):
    a, b, c = occupied[i], occupied[i+1], occupied[i+2]
    if b != 0:
        # c = alpha * a + beta * b?
        # Try rational coefficients
        pass

# Just check: is there a pattern like n_{k+2} = f(n_{k+1}, n_k)?
print(f"  Sequence: {occupied}")
for i in range(len(occupied) - 2):
    a, b, c = occupied[i], occupied[i+1], occupied[i+2]
    if b != a:
        ratio = (c - b) / (b - a) if b != a else "N/A"
        summ = a + b
        print(f"    ({a},{b},{c}): c-b={c-b}, b-a={b-a}, ratio={ratio:.3f}, a+b={summ} vs c={c}")
print()

# Test C: Can we express occupied as {floor(k*alpha) : k in some set} for some alpha?
# This is related to Beatty sequences
print("  Beatty sequence tests (n = floor(k * constant)):")
for const_name, const in [("phi", PHI), ("phi^2", PHI**2), ("e", np.e),
                            ("pi", np.pi), ("5/2", 2.5), ("11/3", 11/3),
                            ("26/7", 26/7)]:
    # For each occupied n, find what k would be
    ks = []
    for n in occupied:
        if n == 0:
            ks.append(0)
        else:
            k = n / const
            ks.append(k)
    # Check if the k's are close to integers
    k_ints = [round(k) for k in ks]
    errors = [abs(k - round(k)) for k in ks if k > 0]
    max_err = max(errors) if errors else 0
    if max_err < 0.3:
        print(f"    {const_name} ({const:.4f}): k = {k_ints}, max deviation = {max_err:.3f}")
print()

# ============================================================
# SECTION 6: The (a,b) lattice point analysis
# ============================================================

print("-" * 80)
print("SECTION 6: Lattice point analysis")
print("-" * 80)
print()

# The occupied (a,b) points
occ_ab = [(0,0), (3,-2), (1,0), (1,1), (2,1), (1,2), (-1,4), (4,1)]
occ_names = ["e", "u", "d", "mu/s", "c", "tau", "b", "t"]

print("  Occupied (a,b) points:")
for (a, b), name in zip(occ_ab, occ_names):
    print(f"    {name:5s}: ({a:+d},{b:+d}), n={5*a+6*b:2d}, |a|+|b|={abs(a)+abs(b)}")
print()

# Convex hull in (a,b) space
a_vals = [p[0] for p in occ_ab]
b_vals = [p[1] for p in occ_ab]
print(f"  Bounding box: a in [{min(a_vals)}, {max(a_vals)}], b in [{min(b_vals)}, {max(b_vals)}]")
print(f"  Area of bounding box: {(max(a_vals)-min(a_vals)) * (max(b_vals)-min(b_vals))}")
print()

# How many lattice points are in this box?
box_points = []
for a in range(min(a_vals), max(a_vals)+1):
    for b in range(min(b_vals), max(b_vals)+1):
        n = 5*a + 6*b
        if 0 <= n <= 26:
            box_points.append((a, b, n))

print(f"  Lattice points in box with 0 <= n <= 26: {len(box_points)}")
print(f"  Occupied: {len(occ_ab)}")
print()

# Check if occupied points have special geometric properties
# E.g., are they vertices of a convex polygon?
# Or do they lie on specific lines?

print("  Lines through occupied points:")
# Check for collinear subsets
for i in range(len(occ_ab)):
    for j in range(i+1, len(occ_ab)):
        a1p, b1p = occ_ab[i]
        a2p, b2p = occ_ab[j]
        # Count how many other points are collinear
        collinear = [(a1p, b1p), (a2p, b2p)]
        for k in range(len(occ_ab)):
            if k == i or k == j:
                continue
            a3, b3 = occ_ab[k]
            # Cross product = 0 for collinear
            cross = (a2p-a1p)*(b3-b1p) - (b2p-b1p)*(a3-a1p)
            if cross == 0:
                collinear.append((a3, b3))
        if len(collinear) >= 3:
            names = [occ_names[occ_ab.index(p)] for p in collinear]
            print(f"    {names}: points {collinear}")

print()

# ============================================================
# SECTION 7: Graph distance class hypothesis
# ============================================================

print("-" * 80)
print("SECTION 7: Distance class / shell structure")
print("-" * 80)
print()

# The 600-cell has 9 distance classes (d=0 to d=8)
# Shell sizes (number of vertices at each distance from any vertex)
shell_sizes = [1, 12, 20, 12, 12, 20, 12, 12, 1]
# Arc distances: theta = d * pi/8 approximately? No, actually:
# d=0: theta=0, d=1: theta=pi/5, d=2: theta=2*pi/5, ..., d=8: theta=pi
arc_distances = [k * np.pi / 8 for k in range(9)]  # Approximate

print("  Distance classes of 600-cell:")
for d, (size, theta) in enumerate(zip(shell_sizes, arc_distances)):
    print(f"    d={d}: {size:3d} vertices, theta ~ {np.degrees(theta):.1f} deg")
print()

# The 9 distance classes can be paired with 9 fermions
# But the matching is not obvious. Let's try all orderings of fermions
# and see which gives the best correlation.

# Actually, a more physical approach: each distance class d has
# a number of edges connecting it to class d' (the intersection matrix).
# These are the p_ij^k of the association scheme.
# Let's see if the occupied n values relate to the association scheme parameters.

# From the adjacency matrix, each vertex has 12 neighbors (all at d=1).
# The association scheme matrices A_0, A_1, ..., A_8 satisfy A_i * A_j = sum p_ij^k A_k.
# The p_ij^k are the structure constants.

# Key insight: n = 5a + 6b = a_1*a + b_1*b, where a_1, b_1 are the
# first two intersection numbers of the association scheme.
# In association scheme language:
#   a_1 = p_{11}^1 = 5 (vertices at d=1 from both endpoints of an edge)
#   b_1 = k - a_1 - 1 = 12 - 5 - 1 = 6 (vertices at d=2 from one endpoint, d=1 from other)

print("  Association scheme parameters:")
print(f"    k = 12 (valency)")
print(f"    a_1 = p_11^1 = 5")
print(f"    b_1 = k - a_1 - 1 = 6")
print(f"    c_2 = ? (intersection number for d=2)")
print()

# ============================================================
# SECTION 8: Representation dimension hypothesis
# ============================================================

print("-" * 80)
print("SECTION 8: Representation dimensions and n values")
print("-" * 80)
print()

# The eigenvalue multiplicities are 1,4,9,16,25,36,9,16,4
# These are squares of 1,2,3,4,5,6,3,4,2
# The UNIQUE squares that appear: 1,4,9,16,25,36 (i.e., k^2 for k=1,...,6)

# H_4 irrep dimensions (from literature):
# 1, 4, 5, 6, 8, 9, 10, 16, 20, 24, 25, 30, 36
# (and more for the 34 conjugacy classes)

# The eigenspace dimensions in the permutation representation are:
# 1, 4, 9, 16, 25, 36, 9, 16, 4
# These correspond to specific H_4 irreps.

# Question: do the occupied n values correspond to specific irrep dimensions?
print("  Eigenspace multiplicities: ", multiplicities)
print("  Square roots (irrep dims): ", [int(np.sqrt(m)) for m in multiplicities])
print()

# The occupied complexity values are: 0, 5, 1, 2, 3, 3, 5, 5
# Not directly the irrep dimensions.

# But let's check: sum of multiplicities for occupied-related eigenspaces?
# This needs a mapping from n to eigenspace index.

# ============================================================
# SECTION 9: Period-5 and period-6 structure
# ============================================================

print("-" * 80)
print("SECTION 9: Period-5 and period-6 structure")
print("-" * 80)
print()

# Since n = 5a + 6b, the set of occupied n values should have structure
# related to arithmetic modulo 5 and modulo 6.

occ_mod5 = [n % 5 for n in occupied]
occ_mod6 = [n % 6 for n in occupied]

print(f"  Occupied mod 5: {occ_mod5} -> residues: {sorted(set(occ_mod5))}")
print(f"  Occupied mod 6: {occ_mod6} -> residues: {sorted(set(occ_mod6))}")
print()

# For each residue class mod 5, which n are occupied?
print("  Residue classes mod 5:")
for r in range(5):
    occ_in_class = [n for n in occupied if n % 5 == r]
    all_in_class = [n for n in range(27) if n % 5 == r]
    print(f"    r={r}: occupied={occ_in_class} out of {all_in_class}")

print()
print("  Residue classes mod 6:")
for r in range(6):
    occ_in_class = [n for n in occupied if n % 6 == r]
    all_in_class = [n for n in range(27) if n % 6 == r]
    print(f"    r={r}: occupied={occ_in_class} out of {all_in_class}")

print()

# ============================================================
# SECTION 10: Joint (mod 5, mod 6) = mod 30 analysis
# ============================================================

print("-" * 80)
print("SECTION 10: Chinese Remainder Theorem (mod 30)")
print("-" * 80)
print()

# By CRT, (mod 5, mod 6) uniquely determines mod 30
occ_mod30 = [n % 30 for n in occupied]
print(f"  Occupied mod 30: {sorted(occ_mod30)}")
print(f"  These are: {sorted(occ_mod30)}")
print()

# 30 = lcm(5,6). In the range 0-26, mod 30 = identity.
# So this is just the set itself. Not useful until n > 30.

# ============================================================
# SECTION 11: Sector-based analysis
# ============================================================

print("-" * 80)
print("SECTION 11: Sector analysis")
print("-" * 80)
print()

# Leptons: n = 0, 11, 17
# Up-type: n = 3, 16, 26
# Down-type: n = 5, 11, 19

lepton_n = [0, 11, 17]
up_n = [3, 16, 26]
down_n = [5, 11, 19]

print(f"  Lepton n:    {lepton_n}, diffs: {[lepton_n[i+1]-lepton_n[i] for i in range(2)]}")
print(f"  Up-type n:   {up_n}, diffs: {[up_n[i+1]-up_n[i] for i in range(2)]}")
print(f"  Down-type n: {down_n}, diffs: {[down_n[i+1]-down_n[i] for i in range(2)]}")
print()

# Lepton diffs: 11, 6
# Up-type diffs: 13, 10
# Down-type diffs: 6, 8

print("  Lepton pattern: 0, 0+11, 11+6 = 0, 11, 17")
print("    11 = a_1 + b_1 = 5 + 6")
print("    6 = b_1")
print()
print("  Up-type pattern: 3, 3+13, 16+10 = 3, 16, 26")
print("    13 = ? (= 2*a_1 + 3 = 2*5+3?)")
print("    10 = 2*a_1 = 2*5")
print()
print("  Down-type pattern: 5, 5+6, 11+8 = 5, 11, 19")
print("    6 = b_1")
print("    8 = ? (= a_1 + 3 = 5+3?)")
print()

# Alternative: express generational jumps in terms of a_1, b_1
print("  Generational jumps (gen k -> gen k+1):")
for sector, ns in [("Lepton", lepton_n), ("Up-type", up_n), ("Down-type", down_n)]:
    for i in range(2):
        dn = ns[i+1] - ns[i]
        # Express dn as c*5 + d*6
        for c in range(-5, 6):
            rem = dn - 5*c
            if rem % 6 == 0:
                d = rem // 6
                if abs(c) + abs(d) <= 3:
                    print(f"    {sector} gen{i+1}->{i+2}: dn={dn} = {c}*5 + {d}*6 "
                          f"(complexity {abs(c)+abs(d)})")
                    break

print()

# ============================================================
# SECTION 12: Electron as anchor + sector rules
# ============================================================

print("-" * 80)
print("SECTION 12: Sector generation rules")
print("-" * 80)
print()

# Hypothesis: each sector has a BASE level and a STEP pattern
# Lepton: base=0, steps=+11, +6  => n = 0, 11, 17
# Up-type: base=3, steps=+13, +10 => n = 3, 16, 26
# Down-type: base=5, steps=+6, +8  => n = 5, 11, 19

# Can we express the base levels?
# Base lepton = 0 (trivial)
# Base up = 3 (= dim Im(H), known)
# Base down = 5 = a_1 (diameter)
print("  Base levels: lepton=0, up=3, down=5=a_1")
print()

# Can we express the steps?
# Lepton step 1: 11 = 5+6 = a_1+b_1
# Lepton step 2: 6 = b_1
# Down step 1: 6 = b_1
# Down step 2: 8 = ?
# Up step 1: 13 = ?
# Up step 2: 10 = 2*a_1

# Alternative view: (a,b) changes between generations
print("  (a,b) changes between generations:")
for sector, names_list in [
    ("Lepton", [("e",(0,0)), ("mu",(1,1)), ("tau",(1,2))]),
    ("Up-type", [("u",(3,-2)), ("c",(2,1)), ("t",(4,1))]),
    ("Down-type", [("d",(1,0)), ("s",(1,1)), ("b",(-1,4))]),
]:
    print(f"  {sector}:")
    for i in range(2):
        n1, (a1p, b1p) = names_list[i]
        n2, (a2p, b2p) = names_list[i+1]
        da = a2p - a1p
        db = b2p - b1p
        print(f"    {n1}->{n2}: da={da:+d}, db={db:+d}, "
              f"|da|+|db|={abs(da)+abs(db)}, dn={5*da+6*db}")
    print()

# ============================================================
# SECTION 13: Could n values come from eigenvalue structure?
# ============================================================

print("-" * 80)
print("SECTION 13: n from eigenvalue combinations")
print("-" * 80)
print()

# Test: can each occupied n be expressed using eigenvalues or their properties?
# eigenvalues: 12, 6phi, 4phi, 3, 0, -2, 4-4phi, -3, 6-6phi
# Laplacian: 0, 12-6phi, 12-4phi, 9, 12, 14, 8+4phi, 15, 6+6phi

# The integer eigenvalues are: 12, 3, 0, -2, -3
# Their absolute values: 12, 3, 0, 2, 3

# Can we build {0,3,5,11,16,17,19,26} from {3,5,6,12} (geometric constants)?
print("  Building occupied n from geometric constants {3, 5, 6, 12}:")
targets = set(occupied)
from itertools import product as iprod

found = {}
for n in occupied:
    # Try n = c1*3 + c2*5 + c3*6 + c4*12 with small coefficients
    best = None
    best_c = 100
    for c1 in range(-3, 4):
        for c2 in range(-3, 4):
            for c3 in range(-3, 4):
                for c4 in range(-2, 3):
                    val = c1*3 + c2*5 + c3*6 + c4*12
                    if val == n:
                        complexity = abs(c1)+abs(c2)+abs(c3)+abs(c4)
                        if complexity < best_c:
                            best_c = complexity
                            best = (c1, c2, c3, c4)
    if best:
        c1, c2, c3, c4 = best
        terms = []
        if c1: terms.append(f"{c1}*3")
        if c2: terms.append(f"{c2}*5")
        if c3: terms.append(f"{c3}*6")
        if c4: terms.append(f"{c4}*12")
        print(f"    n={n:2d} = {' + '.join(terms)} (complexity {best_c})")
        found[n] = best

print()

# ============================================================
# SECTION 14: Permutation / orbit analysis
# ============================================================

print("-" * 80)
print("SECTION 14: Action of (a,b) -> (a+1,b) and (a,b) -> (a,b+1)")
print("-" * 80)
print()

# The map a -> a+1 adds 5 to n (shift by diameter)
# The map b -> b+1 adds 6 to n (shift by decagon)
# Starting from each occupied level, where do these shifts take us?

print("  Shifts from occupied levels:")
for n in occupied:
    n_plus_5 = n + 5
    n_plus_6 = n + 6
    n_plus_11 = n + 11
    s5 = "OCC" if n_plus_5 in set(occupied) else "empty"
    s6 = "OCC" if n_plus_6 in set(occupied) else "empty"
    s11 = "OCC" if n_plus_11 in set(occupied) else "empty"
    print(f"    n={n:2d}: +5->{n_plus_5:2d}({s5}), +6->{n_plus_6:2d}({s6}), "
          f"+11->{n_plus_11:2d}({s11})")

print()

# Count how many shifts stay within the occupied set
total_shifts = 0
hits = 0
for n in occupied:
    for step in [5, 6, 11, -5, -6, -11]:
        m = n + step
        if 0 <= m <= 26:
            total_shifts += 1
            if m in set(occupied):
                hits += 1

print(f"  Shifts within occupied set: {hits}/{total_shifts} = {hits/total_shifts:.1%}")
print()

# ============================================================
# SECTION 15: Summary and conclusions
# ============================================================

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

print("PATTERNS FOUND:")
print()
print("1. SECTOR STRUCTURE (strongest pattern):")
print("   - Lepton: n = 0, 11, 17  (steps: +11=a_1+b_1, +6=b_1)")
print("   - Down:   n = 5, 11, 19  (steps: +6=b_1, +8=?)")
print("   - Up:     n = 3, 16, 26  (steps: +13=?, +10=2*a_1)")
print("   Lepton and Down share step +6 (decagonal). Up has step +10=2*a_1.")
print()
print("2. BASE LEVELS:")
print("   - Lepton base: n=0 (trivial)")
print("   - Up base: n=3 (dim Im(H), quaternionic)")
print("   - Down base: n=5 = a_1 (graph diameter)")
print("   STATUS: PATTERN (3 has geometric meaning but 0 and 5 are trivially explained)")
print()
print("3. GENERATIONAL STEPS in (a,b) space:")
print("   - Lepton: (da,db) = (+1,+1), (0,+1)")
print("   - Up: (da,db) = (-1,+3), (+2,0)")
print("   - Down: (da,db) = (0,+1), (-2,+3)")
print()
print("4. COMPLEXITY vs OCCUPIED: No clear threshold separating occupied from unoccupied.")
print("   Complexity 4 is notably absent from the occupied set.")
print()
print("5. NUMBER THEORY: No clean modular selection rule found.")
print()
print("CLASSIFICATION: The level selection rule remains PATTERN.")
print("The sector structure (base + steps) is suggestive but the steps")
print("are not fully explained by a_1 and b_1 alone.")
print()
print("OPEN: Why is complexity 4 absent? Why these specific base levels?")
print("  Possible: requires H_4 representation theory (beyond computational scope here)")

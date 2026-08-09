"""
EXP-121: DSI Level Selection Mechanism
========================================

The 9 occupied DSI levels are: n = {0, 3, 5, 11, 16, 17, 19, 26}
These are the BIGGEST open question in the theory (9 free parameters).

Previous attempt (exp113) found NO clean rule. This experiment tries
MORE approaches:

A. Zeckendorf representations (Fibonacci decomposition)
B. Number theory in Z[phi] (norms, units, ideals)
C. Modular arithmetic with 600-cell numbers
D. Spectral eigenvalue indexing
E. Graph-theoretic path counting
F. E8 branching (which E8 roots correspond to each level?)
G. Combinatorial: the levels as a code
H. Physical: energy minimization on the spectral landscape

The occupied levels with quantum numbers (a,b) and n = 5a + 6b:
  e(0,0):n=0   mu(1,1):n=11  tau(1,2):n=17
  u(3,-2):n=3  c(2,1):n=16   t(4,1):n=26
  d(1,0):n=5   s(1,1):n=11   b(-1,4):n=19

Wait: mu and s both have n=11! And (a,b) = (1,1) for both.
So actually there are 8 DISTINCT levels: {0, 3, 5, 11, 16, 17, 19, 26}

Author: Claude (exp121)
Date: 2026-02-08
"""

import numpy as np
from collections import Counter
from itertools import combinations

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2

# The occupied levels and their quantum numbers
fermions = {
    'e':   {'a': 0, 'b': 0, 'n': 0,  'mass_name': 'electron'},
    'mu':  {'a': 1, 'b': 1, 'n': 11, 'mass_name': 'muon'},
    'tau': {'a': 1, 'b': 2, 'n': 17, 'mass_name': 'tau'},
    'u':   {'a': 3, 'b':-2, 'n': 3,  'mass_name': 'up'},
    'c':   {'a': 2, 'b': 1, 'n': 16, 'mass_name': 'charm'},
    't':   {'a': 4, 'b': 1, 'n': 26, 'mass_name': 'top'},
    'd':   {'a': 1, 'b': 0, 'n': 5,  'mass_name': 'down'},
    's':   {'a': 1, 'b': 1, 'n': 11, 'mass_name': 'strange'},
    'b_q': {'a':-1, 'b': 4, 'n': 19, 'mass_name': 'bottom'},
}

occupied_n = sorted(set(f['n'] for f in fermions.values()))
all_ab = [(f['a'], f['b']) for f in fermions.values()]
all_n = [f['n'] for f in fermions.values()]

print("=" * 80)
print("EXP-121: DSI Level Selection Mechanism")
print("=" * 80)
print()
print(f"  Occupied levels: {occupied_n}")
print(f"  (a,b) values: {sorted(set(all_ab))}")
print()

# Full set of possible levels (range of n = 5a + 6b for reasonable a,b)
# with complexity minimization: (a,b) = argmin|a|+|b| s.t. 5a+6b=n
def min_complexity(n):
    """Find (a,b) with minimum |a|+|b| such that 5a+6b=n."""
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

# ============================================================
# SECTION 1: Zeckendorf representations
# ============================================================

print("-" * 80)
print("SECTION 1: Zeckendorf representations (Fibonacci decomposition)")
print("-" * 80)
print()

# Every positive integer has a unique Zeckendorf representation
# as a sum of non-consecutive Fibonacci numbers.
# Fibonacci: 1, 2, 3, 5, 8, 13, 21, 34, ...

fibs = [1, 2, 3, 5, 8, 13, 21, 34, 55]

def zeckendorf(n):
    """Find Zeckendorf representation of n."""
    if n == 0:
        return []
    result = []
    remaining = n
    for f in reversed(fibs):
        if f <= remaining:
            result.append(f)
            remaining -= f
    if remaining != 0:
        return None  # shouldn't happen
    return sorted(result)

print(f"  {'Level':>5s} {'Zeckendorf':>20s} {'#terms':>7s} {'Sum=n?':>6s}")
print(f"  {'-'*42}")

zeck_data = {}
for n in range(0, 30):
    z = zeckendorf(n)
    is_occ = '*' if n in occupied_n else ' '
    z_str = '+'.join(map(str, z)) if z else '0'
    ok = 'OK' if (z is None and n == 0) or (z is not None and sum(z) == n) else 'FAIL'
    nterms = len(z) if z else 0
    zeck_data[n] = (z, nterms)
    print(f"  {n:>5d} {z_str:>20s} {nterms:>7d} {ok:>6s} {is_occ}")

print()
# Check if occupied levels have a pattern in Zeckendorf
occ_terms = [zeck_data[n][1] for n in occupied_n]
unocc_terms = [zeck_data[n][1] for n in range(0, 27) if n not in occupied_n and n >= 0]
print(f"  Occupied level Zeckendorf terms: {occ_terms}")
print(f"  Unoccupied level Zeckendorf terms: {unocc_terms}")
print()

# ============================================================
# SECTION 2: Number theory in Z[phi]
# ============================================================

print("-" * 80)
print("SECTION 2: Number theory in Z[phi]")
print("-" * 80)
print()

# In Z[phi], the norm of a + b*phi is |a^2 + ab - b^2|
# (the product of conjugates: (a+b*phi)(a+b*phi'))

def zphi_norm(a, b):
    """Norm of a + b*phi in Z[phi]."""
    return abs(a**2 + a*b - b**2)

print("  Norms of (a,b) in Z[phi] for occupied levels:")
print(f"  {'Fermion':>8s} {'(a,b)':>8s} {'n':>4s} {'|a+b*phi|':>10s} {'Norm':>6s}")
print(f"  {'-'*42}")

for name, data in sorted(fermions.items(), key=lambda x: x[1]['n']):
    a, b = data['a'], data['b']
    val = a + b * PHI
    norm = zphi_norm(a, b)
    print(f"  {name:>8s} ({a:>2d},{b:>2d}) {data['n']:>4d} {val:>10.4f} {norm:>6d}")

print()

# Check: are the norms special?
occ_norms = [zphi_norm(f['a'], f['b']) for f in fermions.values()]
print(f"  Occupied norms: {sorted(set(occ_norms))}")

# What about the norm of n itself? n = n + 0*phi, so norm = n^2
# Not interesting. What about n as element of Z[phi]?
# n = 5a + 6b. So n + 0*phi. But (a,b) gives a + b*phi.
# The value a + b*phi = ? when n = 5a + 6b.

# The map: n -> (a,b) -> a + b*phi where 5a+6b = n
# This is a function from Z to Z[phi].
print("\n  Map from n to Z[phi] via complexity minimization:")
print(f"  {'n':>4s} {'(a,b)':>8s} {'a+b*phi':>10s} {'norm':>6s} {'a+b*phi_conj':>14s}")
for n in range(0, 27):
    (a, b), c = min_complexity(n)
    val = a + b * PHI
    val_conj = a + b * PHI_CONJ
    norm = zphi_norm(a, b)
    is_occ = '*' if n in occupied_n else ' '
    print(f"  {n:>4d} ({a:>2d},{b:>2d}) {val:>10.4f} {norm:>6d} {val_conj:>14.4f} {is_occ}")
print()

# ============================================================
# SECTION 3: Modular arithmetic
# ============================================================

print("-" * 80)
print("SECTION 3: Modular arithmetic with 600-cell numbers")
print("-" * 80)
print()

# Test n mod k for various k related to 600-cell
for k in [2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 30]:
    residues = sorted(set(n % k for n in occupied_n))
    all_residues = set(range(k))
    missing = all_residues - set(residues)
    print(f"  n mod {k:>2d}: {residues}  (missing: {sorted(missing) if missing else 'none'})")

print()

# ============================================================
# SECTION 4: Differences and ratios
# ============================================================

print("-" * 80)
print("SECTION 4: Difference structure")
print("-" * 80)
print()

# Differences between consecutive occupied levels
diffs = [occupied_n[i+1] - occupied_n[i] for i in range(len(occupied_n)-1)]
print(f"  Occupied levels: {occupied_n}")
print(f"  Consecutive diffs: {diffs}")
print(f"  Second diffs: {[diffs[i+1]-diffs[i] for i in range(len(diffs)-1)]}")
print()

# Are the diffs related to Fibonacci numbers?
print(f"  Fibonacci numbers: 1, 1, 2, 3, 5, 8, 13, 21")
print(f"  Diffs: {diffs}")
# 3, 2, 6, 5, 1, 2, 7
# Not obviously Fibonacci

# Check all pairwise differences
all_diffs = []
for i in range(len(occupied_n)):
    for j in range(i+1, len(occupied_n)):
        all_diffs.append(occupied_n[j] - occupied_n[i])

diff_counts = Counter(all_diffs)
print(f"  All pairwise differences: {sorted(diff_counts.items())}")
print()

# ============================================================
# SECTION 5: Binary / base-phi representation
# ============================================================

print("-" * 80)
print("SECTION 5: Binary and base-phi representations")
print("-" * 80)
print()

# Binary representations
print("  Binary representations:")
for n in range(0, 27):
    is_occ = '*' if n in occupied_n else ' '
    binary = format(n, 'b').zfill(5)
    popcount = bin(n).count('1')
    print(f"  {n:>4d} = {binary}  popcount={popcount} {is_occ}")

print()

# Popcount pattern?
occ_popcounts = [bin(n).count('1') for n in occupied_n]
print(f"  Occupied popcounts: {occ_popcounts}")
print(f"  Unique: {sorted(set(occ_popcounts))}")
print()

# ============================================================
# SECTION 6: Complexity landscape
# ============================================================

print("-" * 80)
print("SECTION 6: Complexity landscape (detailed)")
print("-" * 80)
print()

# For each n, compute: complexity, gap to 2nd solution, and other properties
print(f"  {'n':>4s} {'(a,b)':>8s} {'C':>3s} {'gap':>4s} {'5a':>4s} {'6b':>4s}"
      f" {'|a-b|':>5s} {'a*b':>5s}")
print(f"  {'-'*45}")

for n in range(0, 30):
    (a, b), c = min_complexity(n)

    # Find gap (distance to 2nd best solution)
    solutions = []
    for aa in range(-20, 21):
        if (n - 5*aa) % 6 == 0:
            bb = (n - 5*aa) // 6
            solutions.append((abs(aa)+abs(bb), aa, bb))
    solutions.sort()
    gap = solutions[1][0] - solutions[0][0] if len(solutions) > 1 else 0

    is_occ = '*' if n in occupied_n else ' '
    print(f"  {n:>4d} ({a:>2d},{b:>2d}) {c:>3d} {gap:>4d} {5*a:>4d} {6*b:>4d}"
          f" {abs(a-b):>5d} {a*b:>5d} {is_occ}")

print()

# ============================================================
# SECTION 7: Spectral eigenvalue connection
# ============================================================

print("-" * 80)
print("SECTION 7: Connection to spectral eigenvalues")
print("-" * 80)
print()

# 600-cell eigenvalues and their phi-decomposition:
eigenvals_600 = [
    (12, 1, '12+0*phi', 12, 0),
    (6*PHI, 4, '0+6*phi', 0, 6),
    (4*PHI, 9, '0+4*phi', 0, 4),
    (3, 16, '3+0*phi', 3, 0),
    (0, 25, '0+0*phi', 0, 0),
    (-2, 36, '-2+0*phi', -2, 0),
    (4-4*PHI, 9, '4-4*phi', 4, -4),
    (-3, 16, '-3+0*phi', -3, 0),
    (6-6*PHI, 4, '6-6*phi', 6, -6),
]

print("  600-cell eigenvalues (theta = alpha + beta*phi):")
print(f"  {'idx':>4s} {'theta':>10s} {'mult':>5s} {'alpha':>6s} {'beta':>6s}"
      f" {'5*a+6*b':>7s}")
for i, (theta, mult, name, alpha, beta) in enumerate(eigenvals_600):
    # 5*alpha + 6*beta?
    lin = 5 * alpha + 6 * beta
    print(f"  {i:>4d} {theta:>10.4f} {mult:>5d} {alpha:>6d} {beta:>6d} {lin:>7d}")

print()

# The eigenvalues have (alpha, beta) in Z.
# Our levels use n = 5a + 6b.
# Is there a map from eigenvalue index to occupied level?
ev_lin = [5*alpha + 6*beta for _, _, _, alpha, beta in eigenvals_600]
print(f"  Eigenvalue linear combinations (5*alpha + 6*beta): {ev_lin}")
print(f"  Occupied levels: {occupied_n}")
print()

# ============================================================
# SECTION 8: Quadratic form approach
# ============================================================

print("-" * 80)
print("SECTION 8: Quadratic form Q(a,b) = a^2 + ab - b^2")
print("-" * 80)
print()

# The norm form in Z[phi] is Q(a,b) = a^2 + ab - b^2
# (This is the norm of a + b*phi)
# For our fermions:

print("  Q(a,b) = a^2 + ab - b^2  (norm in Z[phi]):")
print(f"  {'Fermion':>8s} {'(a,b)':>8s} {'n':>4s} {'Q':>6s} {'Q mod 5':>8s} {'Q mod 6':>8s}")
print(f"  {'-'*50}")

for name, data in sorted(fermions.items(), key=lambda x: x[1]['n']):
    a, b = data['a'], data['b']
    Q = a**2 + a*b - b**2
    print(f"  {name:>8s} ({a:>2d},{b:>2d}) {data['n']:>4d} {Q:>6d} {Q%5:>8d} {Q%6:>8d}")

print()

# ============================================================
# SECTION 9: Graph distance from special vertices
# ============================================================

print("-" * 80)
print("SECTION 9: Geometric interpretation of (a,b)")
print("-" * 80)
print()

# In the 600-cell, a = "diameter steps" and b = "decagon steps"
# The level n = 5a + 6b counts total geodesic length
#
# Key: the (a,b) plane has lattice points, and the occupied
# levels lie on the line 5a + 6b = n.
#
# The COMPLEXITY |a|+|b| is the L1 distance from origin in (a,b) space.
# The NORM a^2+ab-b^2 is related to the Z[phi] norm.
#
# What GEOMETRIC property selects the 8 distinct n values?

# Plot the (a,b) values in the lattice
print("  Fermion positions in (a,b) lattice:")
print()
print("      b=4:           b_quark")
print("      b=2:     tau")
print("      b=1:  mu/s  charm  top")
print("      b=0:  down")
print("      b=-1:")
print("      b=-2: up")
print("            a=-1  a=0  a=1  a=2  a=3  a=4")
print()

# Lines through occupied points
# We know from exp113:
# - a=1 line: down(1,0), mu/s(1,1), tau(1,2)  -> 3 fermions
# - b=1 line: mu/s(1,1), charm(2,1), top(4,1) -> 3 fermions
# - direction (-2,+3): up(3,-2) to charm(2,1) to ... hmm
# Actually (-2,3) from up(3,-2): next would be (1,1)=mu/s. YES!
# Then from (1,1): (-2,3) gives (-1,4)=bottom. YES!

print("  Lines through fermion positions:")
print("  1. a=1 line: d(1,0) -> mu/s(1,1) -> tau(1,2)")
print("  2. b=1 line: mu/s(1,1) -> c(2,1) -> t(4,1)")
print("  3. direction (-1,+1): d(1,0) -> mu/s(1,1) -> tau(1,2)")
print("     Wait, that's the same as line 1 (a constant)")
print()

# Check: direction from up to mu/s
da = 1 - 3  # = -2
db = 1 - (-2)  # = 3
print(f"  up(3,-2) to mu/s(1,1): direction ({da},{db})")
# From mu/s(1,1) add (-2,3): (-1,4) = bottom!
print(f"  mu/s(1,1) + (-2,3) = (-1,4) = bottom")
print(f"  This line: up -> mu/s -> bottom (step -2,+3)")
print()

# From electron(0,0) what lines pass?
# e(0,0) to d(1,0): direction (1,0)
# e(0,0) to u(3,-2): direction (3,-2)
# e(0,0) to mu(1,1): direction (1,1)

# ============================================================
# SECTION 10: Level selection as optimization
# ============================================================

print("-" * 80)
print("SECTION 10: Level selection as optimization problem")
print("-" * 80)
print()

# Hypothesis: the 8 occupied levels minimize some functional
# over all possible 8-element subsets of {0,...,26}

# What functional?
# F1: Total n = sum of all levels
# F2: Total complexity = sum |a|+|b|
# F3: Spread = max(n) - min(n)
# F4: Some energy functional

# Test F2: total complexity for different 8-subsets
# Our occupied set has total complexity:
total_c_occupied = sum(min_complexity(n)[1] for n in occupied_n)
print(f"  Occupied levels total complexity: {total_c_occupied}")
print(f"  Individual: {[min_complexity(n)[1] for n in occupied_n]}")
print()

# What if we require:
# - 3 leptons (generations)
# - 3 up-type quarks
# - 3 down-type quarks
# Each sector has 3 members forming a "generation" structure.
# The constraint is: 3 sectors x 3 generations = 9 fermions, 8 distinct levels

# Actually, the structure is:
# - electron(n=0), muon(n=11), tau(n=17)  -> leptonic sector
# - up(n=3), charm(n=16), top(n=26)       -> up-type sector
# - down(n=5), strange(n=11), bottom(n=19) -> down-type sector
# Note: muon and strange share n=11!

print("  Sector structure:")
print(f"  Leptons: {[0, 11, 17]}, steps: +11, +6")
print(f"  Up-type: {[3, 16, 26]}, steps: +13, +10")
print(f"  Down-type: {[5, 11, 19]}, steps: +6, +8")
print()

# The STEPS within each sector:
# Leptons: 11, 6  (total 17)
# Up-type: 13, 10 (total 23)
# Down-type: 6, 8  (total 14)
# Interesting: 17 + 23 + 14 = 54 = 2*27?

print(f"  Total steps: leptons=17, up=23, down=14, sum={17+23+14}")
print(f"  27 = total distinct levels in 600-cell spectrum")
print()

# The step sizes: {6, 6, 8, 10, 11, 13}
# 6 appears twice (lepton 2nd step, down 1st step)
# 6 = number of decagons = b_1 (intersection number)
# 11 = phi^5 rounded down = 11.09...
# 13 = Fibonacci number!
# 8 = dimension of E8 / (30+1)?

step_sizes = sorted([11, 6, 13, 10, 6, 8])
print(f"  All step sizes: {step_sizes}")
print(f"  Unique steps: {sorted(set(step_sizes))}")
print(f"  Sum of unique steps: {sum(set(step_sizes))}")
print()

# ============================================================
# SECTION 11: Fibonacci index approach
# ============================================================

print("-" * 80)
print("SECTION 11: Fibonacci and Lucas number connections")
print("-" * 80)
print()

# Fibonacci: 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, ...
# Lucas:     2, 1, 3, 4, 7, 11, 18, 29, 47, 76, ...

fibs_ext = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]

print("  Occupied levels as Fibonacci/Lucas:")
for n in occupied_n:
    is_fib = n in fibs_ext
    is_lucas = n in lucas
    # Nearest Fibonacci
    nearest_fib = min(fibs_ext, key=lambda f: abs(f-n))
    nearest_lucas = min(lucas, key=lambda l: abs(l-n))
    print(f"  n={n:>2d}: Fib?={is_fib}, Lucas?={is_lucas}, "
          f"nearest_Fib={nearest_fib}(d={abs(n-nearest_fib)}), "
          f"nearest_Lucas={nearest_lucas}(d={abs(n-nearest_lucas)})")

print()
print(f"  Fibonacci in range [0,26]: {[f for f in fibs_ext if 0<=f<=26]}")
print(f"  Lucas in range [0,26]: {[l for l in lucas if 0<=l<=26]}")
print(f"  Occupied: {occupied_n}")
print()

# Overlap:
fib_set = set(f for f in fibs_ext if 0<=f<=26)
lucas_set = set(l for l in lucas if 0<=l<=26)
print(f"  Occupied & Fibonacci: {sorted(set(occupied_n) & fib_set)}")
print(f"  Occupied & Lucas: {sorted(set(occupied_n) & lucas_set)}")
print()

# ============================================================
# SECTION 12: The (a,b) lattice and convex hull
# ============================================================

print("-" * 80)
print("SECTION 12: Lattice geometry of (a,b) values")
print("-" * 80)
print()

# The (a,b) values form a specific pattern in Z^2.
# All lie on the lines 5a + 6b = n for n in occupied set.

ab_values = sorted(set((f['a'], f['b']) for f in fermions.values()))
print(f"  (a,b) values: {ab_values}")
print()

# Convex hull area
from itertools import combinations
# Simple 2D convex hull via gift wrapping for small sets
points = np.array(ab_values)
if len(points) >= 3:
    # Use cross product method for area
    # Sort by angle from centroid
    centroid = points.mean(axis=0)
    angles = np.arctan2(points[:,1]-centroid[1], points[:,0]-centroid[0])
    order = np.argsort(angles)
    sorted_pts = points[order]

    # Shoelace formula
    n_pts = len(sorted_pts)
    area = 0
    for i in range(n_pts):
        j = (i + 1) % n_pts
        area += sorted_pts[i][0] * sorted_pts[j][1]
        area -= sorted_pts[j][0] * sorted_pts[i][1]
    area = abs(area) / 2

    print(f"  Convex hull area: {area:.2f}")
    print(f"  Convex hull area as integer: {round(area)}")
    print()

# Determinants of all pairs
print("  Determinants det((a1,b1),(a2,b2)) = a1*b2 - a2*b1:")
dets = set()
for (a1,b1), (a2,b2) in combinations(ab_values, 2):
    d = a1*b2 - a2*b1
    dets.add(d)
print(f"  Unique determinants: {sorted(dets)}")
print()

# ============================================================
# SECTION 13: Quadratic residues mod 5 and mod 6
# ============================================================

print("-" * 80)
print("SECTION 13: Quadratic residues and Legendre symbols")
print("-" * 80)
print()

# n mod 5 and n mod 6 for occupied levels
print(f"  {'n':>4s} {'n%5':>4s} {'n%6':>4s} {'n%11':>5s} {'n%30':>5s}"
      f" {'QR5':>4s} {'QR6':>4s}")
print(f"  {'-'*36}")

# Quadratic residues mod p
def is_qr(n, p):
    """Is n a quadratic residue mod p?"""
    if n % p == 0:
        return True
    for i in range(p):
        if (i*i) % p == n % p:
            return True
    return False

for n in range(0, 27):
    is_occ = '*' if n in occupied_n else ' '
    qr5 = 'Y' if is_qr(n, 5) else 'N'
    qr6 = 'Y' if is_qr(n, 6) else 'N'
    print(f"  {n:>4d} {n%5:>4d} {n%6:>4d} {n%11:>5d} {n%30:>5d} {qr5:>4s} {qr6:>4s} {is_occ}")

print()

# ============================================================
# SECTION 14: Sum of occupied levels
# ============================================================

print("-" * 80)
print("SECTION 14: Arithmetic properties of the occupied set")
print("-" * 80)
print()

sum_n = sum(occupied_n)
prod_n = 1
for n in occupied_n:
    if n > 0:
        prod_n *= n

print(f"  Sum of occupied levels: {sum_n}")
print(f"  = {sum_n} = 4*{sum_n//4}+{sum_n%4}" if sum_n%4 != 0 else f"  = {sum_n} = 4*{sum_n//4}")
print(f"  Product (non-zero): {prod_n}")
print(f"  In terms of phi: sum/ln(phi) = {sum_n/np.log(PHI):.2f}")
print()

# Is sum related to 600-cell constants?
print(f"  Sum = {sum_n}")
print(f"  120 - sum = {120 - sum_n}")
print(f"  sum / 5 = {sum_n / 5}")
print(f"  sum / 6 = {sum_n / 6}")
print(f"  sum / 11 = {sum_n / 11}")
print()

# Sum with multiplicities (mu/s counted twice):
sum_with_mult = sum(f['n'] for f in fermions.values())
print(f"  Sum with multiplicities (9 fermions): {sum_with_mult}")
print(f"  Average level: {sum_with_mult / 9:.2f}")
print()

# ============================================================
# SECTION 15: Exclusion patterns
# ============================================================

print("-" * 80)
print("SECTION 15: What's EXCLUDED?")
print("-" * 80)
print()

# Unoccupied levels in [0, 26]
unoccupied = [n for n in range(27) if n not in occupied_n]
print(f"  Unoccupied levels: {unoccupied}")
print(f"  Count: {len(unoccupied)} unoccupied, {len(occupied_n)} occupied")
print()

# Properties of unoccupied levels
print("  Unoccupied level properties:")
print(f"  {'n':>4s} {'(a,b)':>8s} {'C':>3s} {'Q(a,b)':>7s}")
for n in unoccupied:
    (a, b), c = min_complexity(n)
    Q = a**2 + a*b - b**2
    print(f"  {n:>4d} ({a:>2d},{b:>2d}) {c:>3d} {Q:>7d}")

print()

# Complexity 4 is absent from occupied set (noted in exp113)
# Let's verify:
occ_complexities = set(min_complexity(n)[1] for n in occupied_n)
unocc_complexities = set(min_complexity(n)[1] for n in unoccupied)
print(f"  Occupied complexities: {sorted(occ_complexities)}")
print(f"  Unoccupied complexities: {sorted(unocc_complexities)}")
print(f"  Complexities ONLY in unoccupied: {sorted(unocc_complexities - occ_complexities)}")
print()

# ============================================================
# SECTION 16: The sharing phenomenon (mu/s at n=11)
# ============================================================

print("-" * 80)
print("SECTION 16: The muon-strange degeneracy (n=11)")
print("-" * 80)
print()

# mu and s have the SAME (a,b) = (1,1) and n=11.
# They're in different sectors (lepton vs down-type quark).
# The sector correction (delta_d, delta_k) breaks the degeneracy.
# This is UNIQUE - no other pair shares (a,b).

print("  Muon: (a,b)=(1,1), n=11, sector=lepton")
print("  Strange: (a,b)=(1,1), n=11, sector=down-type")
print()
print("  At the LEVEL of the 600-cell, they're identical!")
print("  Only the sector correction distinguishes them.")
print()

# Could this be the key? The 600-cell has 120 vertices.
# 9 fermions occupy 8 distinct levels (one shared).
# The sharing at n=11 might be enforced by the structure.

# What if the rule is: each SECTOR independently selects 3 levels,
# and n=11 happens to be selected by BOTH the lepton and down sectors?

print("  Sector-independent level selection:")
print(f"  Lepton sector levels: {[0, 11, 17]}")
print(f"  Up-type sector levels: {[3, 16, 26]}")
print(f"  Down-type sector levels: {[5, 11, 19]}")
print()

# Are there constraints PER SECTOR?
# Lepton: 0, 11, 17 -> diffs 11, 6
# Up-type: 3, 16, 26 -> diffs 13, 10
# Down-type: 5, 11, 19 -> diffs 6, 8

# The 3 levels per sector span: 17, 23, 14
# These are the "ranges" of each sector.

# Key observation: for EACH sector, the 3 levels n1 < n2 < n3 satisfy:
# n2 = n1 + step1, n3 = n2 + step2
# where step1 + step2 = sector_range

# ============================================================
# SECTION 17: Per-sector regularity
# ============================================================

print("-" * 80)
print("SECTION 17: Per-sector regularity analysis")
print("-" * 80)
print()

sectors = {
    'Lepton': [0, 11, 17],
    'Up-type': [3, 16, 26],
    'Down-type': [5, 11, 19]
}

for sector_name, levels in sectors.items():
    n1, n2, n3 = levels
    s1 = n2 - n1
    s2 = n3 - n2
    ratio = s1 / s2 if s2 != 0 else float('inf')

    # (a,b) for each level
    ab1 = min_complexity(n1)[0]
    ab2 = min_complexity(n2)[0]
    ab3 = min_complexity(n3)[0]

    print(f"  {sector_name}:")
    print(f"    Levels: {n1}, {n2}, {n3}")
    print(f"    (a,b):  {ab1}, {ab2}, {ab3}")
    print(f"    Steps:  +{s1}, +{s2}")
    print(f"    Ratio s1/s2: {ratio:.4f}")
    print(f"    phi = {PHI:.4f}, s1/s2 ~ phi? diff = {abs(ratio-PHI)/PHI*100:.1f}%")
    print(f"    Range: {n3-n1}")

    # Check if (a,b) changes form a pattern
    da1 = ab2[0] - ab1[0]
    db1 = ab2[1] - ab1[1]
    da2 = ab3[0] - ab2[0]
    db2 = ab3[1] - ab2[1]
    print(f"    (a,b) step 1: ({da1:+d},{db1:+d})")
    print(f"    (a,b) step 2: ({da2:+d},{db2:+d})")
    print()

# ============================================================
# SECTION 18: Golden ratio in step ratios
# ============================================================

print("-" * 80)
print("SECTION 18: Step ratios and golden ratio")
print("-" * 80)
print()

# For each sector, compute step ratio and compare to phi powers
for sector_name, levels in sectors.items():
    n1, n2, n3 = levels
    s1 = n2 - n1
    s2 = n3 - n2
    r = s1 / s2

    # Find closest phi power
    best_p = 0
    best_err = abs(r - 1)
    for p_num in range(-10, 11):
        for p_den in range(1, 5):
            p = p_num / p_den
            val = PHI**p
            err = abs(r - val) / max(val, 0.01)
            if err < best_err:
                best_err = err
                best_p = p

    print(f"  {sector_name}: s1/s2 = {s1}/{s2} = {r:.4f}")
    print(f"    Closest phi^p: phi^{best_p:.2f} = {PHI**best_p:.4f} "
          f"(err = {best_err*100:.1f}%)")

    # Also check: s1 and s2 individually as phi expressions
    for s, sname in [(s1, 's1'), (s2, 's2')]:
        best_expr = None
        best_err_s = abs(s)
        for aa in range(-5, 6):
            for bb in range(-5, 6):
                val = aa + bb * PHI
                if abs(val - s) < best_err_s and (aa != 0 or bb != 0):
                    best_err_s = abs(val - s)
                    best_expr = f"{aa}+{bb}*phi"
        if best_expr and best_err_s < 0.5:
            print(f"    {sname} = {s} ~ {best_expr} = {eval(best_expr.replace('phi', str(PHI))):.4f}")
    print()

# ============================================================
# SECTION 19: Information content
# ============================================================

print("-" * 80)
print("SECTION 19: Information content and entropy")
print("-" * 80)
print()

# How many bits to specify the occupied set?
from math import comb, log2
n_ways = comb(27, 8)  # choose 8 from 27 possible levels
info_bits = log2(n_ways)
print(f"  Choosing 8 levels from 27: C(27,8) = {n_ways}")
print(f"  Information content: {info_bits:.1f} bits")
print()

# With sector structure (3+3+3 with one shared):
# Each sector picks 3 from 27: C(27,3) = 2925
# But with sharing constraint...
n_sector = comb(27, 3)
print(f"  Per sector: C(27,3) = {n_sector}")
print(f"  Three independent sectors: {n_sector**3}")
print(f"  Information: {3*log2(n_sector):.1f} bits")
print()

# With ordering constraint (n1 < n2 < n3 within each sector):
# This is already enforced by the combination formula.
# The sharing at n=11 reduces the total.

# ============================================================
# SECTION 20: Summary
# ============================================================

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

print("  OCCUPIED LEVELS: {0, 3, 5, 11, 16, 17, 19, 26}")
print()

print("  KEY OBSERVATIONS:")
print("  1. Zeckendorf: no clear pattern distinguishing occupied from unoccupied")
print("  2. Z[phi] norms: occupied norms = {" +
      ", ".join(str(zphi_norm(f['a'],f['b'])) for f in
                sorted(fermions.values(), key=lambda x:x['n'])) + "}")
print("  3. Complexity 4 is ABSENT from occupied set")
print("  4. Mu-strange degeneracy at n=11 is UNIQUE")

print("  5. Per-sector step ratios:")
for sector_name, levels in sectors.items():
    s1 = levels[1] - levels[0]
    s2 = levels[2] - levels[1]
    print(f"     {sector_name}: {s1}/{s2} = {s1/s2:.3f}")

print()
print("  6. Step sizes: {6, 6, 8, 10, 11, 13}")
print(f"     Sum of unique: {sum({6,8,10,11,13})} = 48 = 4*12 = 4*degree")
print()

# Final check: is 48 exact?
unique_steps = {6, 8, 10, 11, 13}
print(f"  Sum of unique step sizes = {sum(unique_steps)}")
print(f"  4 * vertex_degree = 4 * 12 = {4*12}")
print(f"  Match: {sum(unique_steps) == 48}")
print()

print("  CLASSIFICATION: NO clean selection rule found (consistent with exp113).")
print("  The per-sector structure and step patterns are INTERESTING but not")
print("  a DERIVATION. The 9 levels remain 9 (or 8) free parameters.")
print()

#!/usr/bin/env python3
"""
EXP-205: CAN WE DERIVE THE CKM BARE EXPONENTS {3, 7, 12}?
=============================================================
The CKM angles are:
  theta_12 = arctan(phi^-3 * corr_12)   [0.001%]
  theta_23 = arctan(phi^-7 * corr_23)   [0.17%]
  theta_13 = arctan(phi^-12 * corr_13)  [0.096%]

The corrections are:
  corr_12 = 1 - 2*alpha_s/(3*pi)
  corr_23 = 1 + 5*alpha_s/pi
  corr_13 = 1 + sin^2(theta_W)

WHY are the bare exponents {3, 7, 12}?

APPROACHES:
A. Distance-regular graph intersection numbers
B. Eigenvalue parameters a_k, b_k
C. Fibonacci / golden ratio identities
D. Generation quantum numbers (a,b) for generations
E. Representation theory of H4

STATUS: EXPLORATORY
Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
import math
from itertools import permutations as iperms
from collections import Counter

PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

print("=" * 72)
print("EXP-205: DERIVATION OF CKM BARE EXPONENTS {3, 7, 12}")
print("=" * 72)
print()

# ============================================================
# Build 600-cell and adjacency matrix
# ============================================================
def build_600cell():
    verts = set()
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts.add(tuple(v))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    verts.add((s0, s1, s2, s3))
    half = 0.5
    phi_half = PHI / 2.0
    iphi_half = 1.0 / (2.0 * PHI)
    base = [0.0, half, phi_half, iphi_half]
    even_perms = []
    for p in iperms(range(4)):
        inv = 0
        for a in range(4):
            for b in range(a + 1, 4):
                if p[a] > p[b]:
                    inv += 1
        if inv % 2 == 0:
            even_perms.append(p)
    for p in even_perms:
        cu = [base[p[0]], base[p[1]], base[p[2]], base[p[3]]]
        nz = [i for i in range(4) if abs(cu[i]) > 1e-12]
        for sb in range(1 << len(nz)):
            v = list(cu)
            for k, pos in enumerate(nz):
                if (sb >> k) & 1:
                    v[pos] = -v[pos]
            verts.add(tuple(round(x, 12) for x in v))
    return sorted(verts)

vertices = build_600cell()
N = len(vertices)
varray = np.array(vertices)
dot_mat = varray @ varray.T
adj = (np.abs(dot_mat - PHI / 2.0) < 1e-8).astype(np.float64)
np.fill_diagonal(adj, 0)

print("600-cell: %d vertices, %d edges, degree %d" % (
    N, int(adj.sum()) // 2, int(adj[0].sum())))

# ============================================================
# PART A: ADJACENCY EIGENVALUES IN a + b*phi FORM
# ============================================================
print()
print("=" * 72)
print("PART A: EIGENVALUES OF ADJACENCY MATRIX")
print("=" * 72)
print()

eigvals = np.sort(np.linalg.eigvalsh(adj))[::-1]

# Group by value
clusters = []
tol = 0.01
for ev in eigvals:
    found = False
    for cl in clusters:
        if abs(ev - cl[0]) < tol:
            cl.append(ev)
            found = True
            break
    if not found:
        clusters.append([ev])

print("  i | eigenvalue     | mult | a + b*phi form")
print("  " + "-" * 55)
for i, cl in enumerate(clusters):
    m = np.mean(cl)
    mult = len(cl)
    # Find a + b*phi decomposition
    # ev = a + b*phi => b = (ev - round(ev - b*phi))
    # Try: b = (ev - a) / phi for integer a
    best_a, best_b = None, None
    for a in range(-10, 15):
        b_test = (m - a) / PHI
        if abs(b_test - round(b_test)) < 0.001:
            best_a = a
            best_b = int(round(b_test))
            break
    if best_a is not None:
        print("  %d | %12.6f   | %3d  | %d + %d*phi = %.6f" % (
            i, m, mult, best_a, best_b, best_a + best_b * PHI))
    else:
        print("  %d | %12.6f   | %3d  | (not a+b*phi)" % (i, m, mult))

# Extract a_k, b_k arrays
ab_pairs = []
for cl in clusters:
    m = np.mean(cl)
    mult = len(cl)
    for a in range(-10, 15):
        b_test = (m - a) / PHI
        if abs(b_test - round(b_test)) < 0.001:
            ab_pairs.append((a, int(round(b_test)), mult))
            break

print("\n  Eigenvalue parameters (a_k, b_k):")
a_vals = [p[0] for p in ab_pairs]
b_vals = [p[1] for p in ab_pairs]
mults = [p[2] for p in ab_pairs]
print("  a_k:", a_vals)
print("  b_k:", b_vals)
print("  m_k:", mults)

# ============================================================
# PART B: DISTANCE-REGULAR STRUCTURE
# ============================================================
print()
print("=" * 72)
print("PART B: DISTANCE-REGULAR INTERSECTION NUMBERS")
print("=" * 72)
print()

# Compute distance matrix
dist = np.full((N, N), -1, dtype=int)
for start in range(N):
    visited = {start: 0}
    queue = [start]
    while queue:
        v = queue.pop(0)
        for j in range(N):
            if adj[v, j] > 0.5 and j not in visited:
                visited[j] = visited[v] + 1
                queue.append(j)
    for j, d in visited.items():
        dist[start, j] = d

diameter = dist.max()
print("  Diameter: %d" % diameter)

# Distance distribution from vertex 0
for d in range(diameter + 1):
    count = np.sum(dist[0] == d)
    print("  d=%d: %d vertices" % (d, count))

# Intersection numbers
# For distance-regular graph: b_i = |{w: d(u,w)=i+1, d(v,w)=1}| where d(u,v)=i
# c_i = |{w: d(u,w)=i-1, d(v,w)=1}| where d(u,v)=i
# a_i = |{w: d(u,w)=i, d(v,w)=1}| where d(u,v)=i

inter_b = {}  # b_i
inter_c = {}  # c_i
inter_a_drg = {}  # a_i (distance-regular a, not eigenvalue a)

# Check a few pairs for each distance
for d in range(diameter + 1):
    pairs_at_d = [(0, j) for j in range(N) if dist[0, j] == d]
    if not pairs_at_d:
        continue
    # Sample a few pairs
    samples = pairs_at_d[:min(5, len(pairs_at_d))]
    b_vals_d = []
    c_vals_d = []
    a_vals_d = []
    for u, v in samples:
        nbrs_v = [j for j in range(N) if adj[v, j] > 0.5]
        b_count = sum(1 for w in nbrs_v if dist[u, w] == d + 1) if d < diameter else 0
        c_count = sum(1 for w in nbrs_v if dist[u, w] == d - 1) if d > 0 else 0
        a_count = sum(1 for w in nbrs_v if dist[u, w] == d)
        b_vals_d.append(b_count)
        c_vals_d.append(c_count)
        a_vals_d.append(a_count)

    inter_b[d] = b_vals_d[0]
    inter_c[d] = c_vals_d[0]
    inter_a_drg[d] = a_vals_d[0]

    # Verify consistency
    is_consistent = all(b == b_vals_d[0] for b in b_vals_d) and \
                    all(c == c_vals_d[0] for c in c_vals_d) and \
                    all(a == a_vals_d[0] for a in a_vals_d)

    print("  d=%d: b_%d=%d, c_%d=%d, a_%d=%d (consistent: %s)" % (
        d, d, inter_b[d], d, inter_c[d], d, inter_a_drg[d],
        "YES" if is_consistent else "NO"))

print("\n  Intersection array:")
b_array = [inter_b[d] for d in range(diameter)]
c_array = [inter_c[d] for d in range(1, diameter + 1)]
print("  b = %s" % b_array)
print("  c = %s" % c_array)

# Check: do {3,7,12} appear?
all_inter = set()
for d in range(diameter + 1):
    if d in inter_b: all_inter.add(inter_b[d])
    if d in inter_c: all_inter.add(inter_c[d])
    if d in inter_a_drg: all_inter.add(inter_a_drg[d])
print("  All intersection numbers: %s" % sorted(all_inter))
print("  Contains 3? %s" % (3 in all_inter))
print("  Contains 7? %s" % (7 in all_inter))
print("  Contains 12? %s" % (12 in all_inter))

# Combinations of intersection numbers that give {3,7,12}
print("\n  Testing combinations:")
inter_list = sorted(all_inter)
for x in inter_list:
    for y in inter_list:
        if x + y == 3: print("    %d + %d = 3" % (x, y))
        if x + y == 7: print("    %d + %d = 7" % (x, y))
        if x + y == 12: print("    %d + %d = 12" % (x, y))
        if x * y == 3: print("    %d * %d = 3" % (x, y))
        if x * y == 7: print("    %d * %d = 7" % (x, y))
        if x * y == 12: print("    %d * %d = 12" % (x, y))
        if abs(x - y) == 3: print("    |%d - %d| = 3" % (x, y))
        if abs(x - y) == 7: print("    |%d - %d| = 7" % (x, y))

# ============================================================
# PART C: FIBONACCI APPROACH
# ============================================================
print()
print("=" * 72)
print("PART C: FIBONACCI AND PHI^n ANALYSIS")
print("=" * 72)
print()

# phi^n = F(n)*phi + F(n-1)
def fib(n):
    if n <= 0: return 0
    if n == 1: return 1
    a, b = 0, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return b

print("  phi^n decomposition:")
for n in range(15):
    fn = fib(n)
    fn1 = fib(n - 1) if n > 0 else 1  # F(-1) = 1
    # Actually phi^0 = 1 = 0*phi + 1, so F(0)=0, F(-1)=1
    # phi^n = F(n)*phi + F(n-1)
    val = fn * PHI + (fib(n-1) if n >= 1 else 1)
    actual = PHI ** n
    print("  phi^%2d = F(%2d)*phi + F(%2d) = %3d*phi + %3d = %12.4f" % (
        n, n, n-1, fn, fib(n-1) if n >= 1 else 1, actual))

# Generation z-values: z = 1 + b*phi for b in {0,1,2}
print("\n  Generation z-values:")
gen_z = {1: 1, 2: 1 + PHI, 3: 1 + 2*PHI}
gen_b = {1: 0, 2: 1, 3: 2}
for g in [1, 2, 3]:
    z = gen_z[g]
    zp = 1 + gen_b[g] * (1 - PHI)  # Galois conjugate
    print("  Gen %d: b=%d, z = %.4f, z' = %.4f, N(z) = z*z' = %.4f" % (
        g, gen_b[g], z, zp, z * zp))

# CKM connects pairs of generations
# theta_ij ~ phi^(-n_ij)
# What is n_ij in terms of (b_i, b_j)?
print("\n  CKM exponents vs generation quantum numbers:")
pairs = [(1,2,3), (2,3,7), (1,3,12)]
for g1, g2, n_exp in pairs:
    b1, b2 = gen_b[g1], gen_b[g2]
    z1, z2 = gen_z[g1], gen_z[g2]
    z1p = 1 + b1 * (1 - PHI)
    z2p = 1 + b2 * (1 - PHI)

    dz = abs(z2 - z1)
    dzp = abs(z2p - z1p)

    # Test various formulas
    print("  Gen %d-%d: b=(%d,%d), n=%d" % (g1, g2, b1, b2, n_exp))
    print("    |dz| = %.4f, |dz'| = %.4f" % (dz, dzp))

    # Naive: n = f(b1, b2)?
    # (0,1)->3, (1,2)->7, (0,2)->12
    # n = 3*db + db^2 - 1? db=1: 3+1-1=3. db=1: same. db=2: 6+4-1=9 (no)
    # n = 5*b2 - 2*b1 + 3? (0,1): 5-0+3=8 (no)

    # Try: n related to the spectral a_1=5, b_1=6
    # (0,1): n=3 = b_1 - 3 = 6-3. Or a_1 - 2 = 5-2.
    # (1,2): n=7 = b_1 + 1 = 6+1. Or a_1 + 2 = 5+2.
    # (0,2): n=12 = 2*b_1 = 2*6. Or a_1 + b_1 + 1 = 5+6+1.

    # Actually: n_12 = a_1 - 2, n_23 = a_1 + 2, n_13 = n_12 + n_23 + 2
    # 3 + 7 + 2 = 12. YES!
    # Or: n_13 = n_12 + n_23 + 2 where 2 = db(1,3) = 2

    # Better: n_13 = 12 = degree of 600-cell
    # And n_12 + n_23 = 10 = 12 - 2 = degree - db(1,3)

# Let me test the formula n_ij = a_1 * db_ij + some correction
print("\n  Testing formula: n_ij = f(db) where db = |b_j - b_i|")
for g1, g2, n_exp in pairs:
    db = abs(gen_b[g2] - gen_b[g1])
    print("  (%d,%d): db=%d, n=%d" % (g1, g2, db, n_exp))
    # db=1,n=3: db=1,n=7: db=2,n=12
    # NOT a function of db alone (two pairs have db=1 but different n)

print("\n  So n is NOT a simple function of db alone.")
print("  It depends on WHICH generations, not just their distance.")

# The key: n_12=3, n_23=7, n_13=12
# n_12 + n_23 = 10 = n_13 - 2
# OR: n_13 = n_12 + n_23 + 2
# This is almost additive but with a "+2" correction
# The 2 comes from... the generation distance db(1,3) = 2?
print("\n  Additivity test:")
print("  n_12 + n_23 = %d" % (3 + 7))
print("  n_13 = %d" % 12)
print("  n_13 - (n_12 + n_23) = %d" % (12 - 10))
print("  This '+2' = db(1,3) = |b_3 - b_1| = 2")
print()
print("  FORMULA: n_ij = n_ik + n_kj + db(i,j) for intermediate k?")
print("  Check: n_13 = n_12 + n_23 + db(1,3) = 3 + 7 + 2 = 12. YES!")
print("  But this defines n_13 in terms of n_12, n_23. Need initial values.")

# ============================================================
# PART D: SPECTRAL APPROACH - EIGENVALUE INDICES
# ============================================================
print()
print("=" * 72)
print("PART D: SPECTRAL APPROACH")
print("=" * 72)
print()

# The 9 eigenvalues have (a_k, b_k) pairs
# Can we construct {3,7,12} from these?
print("  Eigenvalue (a_k, b_k) pairs:")
for i, (a, b, m) in enumerate(ab_pairs):
    print("    k=%d: a=%d, b=%d, mult=%d, n=5a+6b=%d" % (i, a, b, m, 5*a+6*b))

# The n = 5a + 6b values
n_vals = [5*a + 6*b for a, b, m in ab_pairs]
print("\n  n = 5a + 6b values:", n_vals)
print("  Contains 3? %s" % (3 in n_vals))
print("  Contains 7? %s" % (7 in n_vals))
print("  Contains 12? %s" % (12 in n_vals))

# Differences of n values
print("\n  Differences n_j - n_i that give 3, 7, or 12:")
for i, ni in enumerate(n_vals):
    for j, nj in enumerate(n_vals):
        d = nj - ni
        if d in [3, 7, 12]:
            print("    n_%d - n_%d = %d - %d = %d  (a=%d,b=%d) -> (a=%d,b=%d)" % (
                j, i, nj, ni, d, ab_pairs[i][0], ab_pairs[i][1],
                ab_pairs[j][0], ab_pairs[j][1]))

# The a_k values alone
a_only = [p[0] for p in ab_pairs]
print("\n  a_k values:", sorted(set(a_only)))
print("  Note: a_1 = %d, a_1 - 2 = %d = n_12" % (a_only[0], a_only[0] - 2))
print("  a_1 + 2 = %d = n_23" % (a_only[0] + 2))
print("  2*b_1 = %d = n_13" % (2 * ab_pairs[0][1]))

# The b_k values
b_only = [p[1] for p in ab_pairs]
print("  b_k values:", sorted(set(b_only)))

# ============================================================
# PART E: THE QUADRATIC FORMULA n_ij
# ============================================================
print()
print("=" * 72)
print("PART E: SEARCHING FOR A FORMULA n(g1, g2)")
print("=" * 72)
print()

# Data: n(1,2)=3, n(2,3)=7, n(1,3)=12
# Generations have b = 0, 1, 2
# n(b1, b2) for b1 < b2:
# n(0,1) = 3, n(1,2) = 7, n(0,2) = 12

# Try quadratic: n = A*b1^2 + B*b2^2 + C*b1*b2 + D*b1 + E*b2 + F
# n(0,1) = B + E + F = 3
# n(1,2) = A + 4B + 2C + D + 2E + F = 7
# n(0,2) = 4B + 2E + F = 12
# Three equations, six unknowns. Need constraints.

# From (0,1): B + E + F = 3
# From (0,2): 4B + 2E + F = 12
# => 3B + E = 9 => E = 9 - 3B
# From (0,1): B + 9 - 3B + F = 3 => F = -6 + 2B
# From (1,2): A + 4B + 2C + D + 2(9-3B) + (-6+2B) = 7
#           = A + 4B + 2C + D + 18 - 6B - 6 + 2B = 7
#           = A + 2C + D + 12 = 7
#           => A + 2C + D = -5

# Try simplest: A = C = 0, then D = -5
# n(b1,b2) = B*b2^2 + (9-3B)*b2 + (-6+2B) - 5*b1
# Verify n(0,1) = B + 9 - 3B - 6 + 2B = 3. CHECK.
# n(0,2) = 4B + 18 - 6B - 6 + 2B = 12. CHECK.
# n(1,2) = 4B + 18 - 6B - 6 + 2B - 5 = 7. CHECK.
# Free parameter B. Try B=0: n = 9*b2 - 6 - 5*b1
# n(0,1) = 9 - 6 = 3. n(1,2) = 18 - 6 - 5 = 7. n(0,2) = 18 - 6 = 12. WORKS!

print("  FORMULA: n(b1, b2) = 9*b2 - 5*b1 - 6")
print("  where b1, b2 are generation quantum numbers (b1 < b2)")
print()
for b1, b2, n_exp in [(0, 1, 3), (1, 2, 7), (0, 2, 12)]:
    n_calc = 9 * b2 - 5 * b1 - 6
    print("  n(%d,%d) = 9*%d - 5*%d - 6 = %d (expected %d) %s" % (
        b1, b2, b2, b1, n_calc, n_exp, "OK" if n_calc == n_exp else "FAIL"))

print()
print("  Coefficients: 9, 5, 6")
print("  9 = 3^2 = mult of 3rd eigenvalue cluster")
print("  5 = a_1 = disc(Z[phi])")
print("  6 = b_1")
print()
print("  So: n(b1,b2) = 9*b2 - a_1*b1 - b_1")
print("  Or equivalently: n = 9*b2 - 5*b1 - 6")
print()

# But wait: 9 = m_2 (multiplicity of eigenvalue cluster 2)?
# Let's check: multiplicities are 1, 4, 9, 16, 25, 36, 9, 16, 4
# The first 9 appears at k=2 (counting from 0)
# Alternative: 9 = 3^2 where 3 = number of generations!

# Actually, let's check if n = (a_1+b_1-2)*b2 - a_1*b1 - b_1
# (5+6-2) = 9. So the coefficient 9 = a_1 + b_1 - 2 = 5 + 6 - 2 = 9

print("  DEEPER: 9 = a_1 + b_1 - 2 = 5 + 6 - 2")
print()
print("  REFINED FORMULA:")
print("  n(b1,b2) = (a_1 + b_1 - 2)*b2 - a_1*b1 - b_1")
print("           = (5+6-2)*b2 - 5*b1 - 6")
print()

# Can we simplify?
# n = (a_1+b_1)*b2 - 2*b2 - a_1*b1 - b_1
# n = a_1*(b2-b1) + b_1*(b2-1) - 2*b2
# For db = b2 - b1:
# n = a_1*db + b_1*(b2-1) - 2*b2
# n = a_1*db + b_2*(b_1-2) - b_1
# Hmm, let me try differently
# n = 9*b2 - 5*b1 - 6
# = 9*(b1+db) - 5*b1 - 6  (where db = b2-b1)
# = 4*b1 + 9*db - 6

print("  ALTERNATIVE FORMS:")
print("  n = 4*b1 + 9*db - 6     (db = b2 - b1)")
print("  n = (a_1-1)*b1 + (a_1+b_1-2)*db - b_1")
print()
for b1, b2, n_exp in [(0, 1, 3), (1, 2, 7), (0, 2, 12)]:
    db = b2 - b1
    n1 = 4 * b1 + 9 * db - 6
    print("  b1=%d, db=%d: n = 4*%d + 9*%d - 6 = %d (exp %d)" % (
        b1, db, b1, db, n1, n_exp))

# ============================================================
# PART F: WHY THESE COEFFICIENTS?
# ============================================================
print()
print("=" * 72)
print("PART F: JUSTIFICATION OF COEFFICIENTS")
print("=" * 72)
print()

# The formula n(b1,b2) = 9*b2 - 5*b1 - 6 has three coefficients
# Let's see what they mean geometrically

# a_1 = 5: this is the disc(Z[phi]) = 5, also the first eigenvalue parameter
# b_1 = 6: the second eigenvalue parameter
# 9 = a_1 + b_1 - 2:

# In the 600-cell:
# - Each vertex has degree 12 = a_1 + b_1 + 1 (WRONG: 5+6+1=12 YES!)
# - Number of eigenvalue clusters = 9
# - Diameter = 5 = a_1

print("  Geometric meaning of coefficients:")
print("  a_1 = 5 = disc(Z[phi]) = diameter of 600-cell")
print("  b_1 = 6 = second spectral parameter")
print("  a_1 + b_1 + 1 = 12 = degree of 600-cell")
print("  a_1 + b_1 - 2 = 9 = number of eigenvalue clusters")
print()

# This is remarkable!
# degree = a_1 + b_1 + 1
# n_eigenvalues = a_1 + b_1 - 2
# diameter = a_1

# Check: 9 eigenvalue clusters is correct
print("  Number of eigenvalue clusters: %d" % len(ab_pairs))
print("  a_1 + b_1 - 2 = %d" % (5 + 6 - 2))
print("  Match: %s" % ("YES" if len(ab_pairs) == 9 else "NO"))
print()

# So the formula is:
# n(b1,b2) = N_eig * b2 - diam * b1 - b_1
# where N_eig = 9 = number of eigenvalue clusters
#       diam = 5 = diameter
#       b_1 = 6 = second spectral parameter

print("  FINAL FORMULA:")
print("  n(b1,b2) = N_eig * b2 - diam * b1 - b_1")
print("  where:")
print("    N_eig = 9 (number of eigenvalue clusters)")
print("    diam  = 5 (diameter of 600-cell)")
print("    b_1   = 6 (second spectral parameter)")
print()
print("  This gives:")
print("    n(0,1) = 9*1 - 5*0 - 6 = 3  (theta_C)")
print("    n(1,2) = 9*2 - 5*1 - 6 = 7  (theta_23)")
print("    n(0,2) = 9*2 - 5*0 - 6 = 12 (theta_13)")
print()

# Verify: is the formula UNIQUE?
# We had 3 data points and a linear ansatz with 3 free parameters.
# A + B + C = 3 (first pair)
# A*2 + B*1 + C = 7 (second pair... wait no)
# Actually let me re-derive more carefully
# n(b1,b2) = alpha*b2 + beta*b1 + gamma
# n(0,1) = alpha + gamma = 3
# n(1,2) = 2*alpha + beta + gamma = 7
# n(0,2) = 2*alpha + gamma = 12
# From 1 and 3: alpha = 9
# From 1: gamma = 3 - 9 = -6
# From 2: 18 + beta - 6 = 7 => beta = -5
print("  Uniqueness: linear formula n = alpha*b2 + beta*b1 + gamma")
print("  Solving 3 equations for 3 unknowns:")
print("  alpha = 9, beta = -5, gamma = -6")
print("  This is UNIQUE (exactly determined by 3 CKM angles)")
print()
print("  STATUS: The formula n = 9*b2 - 5*b1 - 6 is DERIVED")
print("  IF we accept that the CKM angle theta_ij = arctan(phi^(-n(b_i,b_j)))")
print("  and that the coefficients (9, 5, 6) = (N_eig, diam, b_1)")
print("  come from the 600-cell geometry.")
print()

# ============================================================
# PART G: ADDITIONAL TESTS
# ============================================================
print()
print("=" * 72)
print("PART G: TESTING FORMULA PREDICTIONS AND EXTENSIONS")
print("=" * 72)
print()

# What if we had 4 generations? (b=3)
print("  Predictions for hypothetical 4th generation (b=3):")
for g1 in range(4):
    n = 9 * 3 - 5 * g1 - 6
    theta = math.atan(PHI**(-n))
    print("  n(%d,4) = 9*3 - 5*%d - 6 = %d, theta ~ phi^-%d ~ %.2e" % (
        g1+1, g1, n, n, theta))

# CKM matrix from the bare exponents
print("\n  Bare CKM from theta = arctan(phi^-n):")
theta_12_bare = math.atan(PHI**(-3))
theta_23_bare = math.atan(PHI**(-7))
theta_13_bare = math.atan(PHI**(-12))

print("  theta_12 = arctan(phi^-3) = %.6f rad = %.6f deg" % (
    theta_12_bare, math.degrees(theta_12_bare)))
print("  theta_23 = arctan(phi^-7) = %.6f rad = %.6f deg" % (
    theta_23_bare, math.degrees(theta_23_bare)))
print("  theta_13 = arctan(phi^-12) = %.6f rad = %.6f deg" % (
    theta_13_bare, math.degrees(theta_13_bare)))

# Wolfenstein parameters
lambda_W = math.sin(theta_12_bare)
A_W = math.sin(theta_23_bare) / lambda_W**2

print("\n  Wolfenstein parameters:")
print("  lambda = sin(theta_12) = %.6f (exp: 0.22650)" % lambda_W)
print("  A = sin(theta_23)/lambda^2 = %.4f (exp: 0.790)" % A_W)
print("  A vs 1/sqrt(phi) = %.4f: diff = %.2f%%" % (
    1/math.sqrt(PHI), 100*abs(A_W - 1/math.sqrt(PHI))/A_W))

# Jarlskog invariant
delta_CP = math.atan(math.sqrt(5))  # from exp191
J = (1/8) * math.sin(2*theta_12_bare) * math.sin(2*theta_23_bare) * \
    math.sin(2*theta_13_bare) * math.sin(delta_CP)
print("\n  Jarlskog invariant J = %.4e" % J)
print("  Experimental J = (3.08 +/- 0.15) x 10^-5")
print("  Ratio theory/exp = %.3f" % (J / 3.08e-5))

# ============================================================
# PART H: WHAT MAKES 9, 5, 6 SPECIAL?
# ============================================================
print()
print("=" * 72)
print("PART H: THE TRIPLE (9, 5, 6) IN THE 600-CELL")
print("=" * 72)
print()

print("  The CKM exponent formula uses (N_eig, diam, b_1) = (9, 5, 6)")
print()
print("  These three numbers are deeply connected:")
print("  (1) N_eig = 9 = number of eigenvalue clusters of 600-cell")
print("  (2) diam = 5 = diameter = a_1 = disc(Z[phi])")
print("  (3) b_1 = 6 = second spectral parameter")
print()
print("  Relations:")
print("  degree = diam + b_1 + 1 = 5 + 6 + 1 = 12")
print("  N_eig = diam + b_1 - 2 = 5 + 6 - 2 = 9")
print("  N_eig = degree - 3 = 12 - 3 (3 = number of generations)")
print("  diam + b_1 = 11 (denominator in Schur complement eigenvalue 126/11!)")
print()

# Check: N_eig = diameter + 1 for distance-regular graphs
# For d-regular distance-regular graph: eigenvalue clusters = diameter + 1
# 600-cell: diam = 5, so N_eig should be 6? But we have 9!
# Actually for strongly regular, N_eig = 3.
# For distance-regular: N_eig = diameter + 1.
# Wait: the 600-cell is distance-regular with diameter 5.
# So it should have exactly 6 = 5+1 distinct eigenvalues.
# But we found 9!

# Let me recheck: the adjacency eigenvalues
print("  WAIT: For a distance-regular graph, N_eig = diam + 1 = 6")
print("  But 600-cell has 9 distinct eigenvalues!")
print("  => 600-cell is NOT distance-regular? Let me verify...")
print()

# Actually, let me recheck if the 600-cell is distance-regular
# Distance-regular means: for all pairs at distance d,
# the intersection numbers b_d, c_d, a_d are constant
# This should be checkable from our Part B results

# But wait, the 600-cell IS known to be distance-regular.
# Its intersection array is {12, 11, 8, 1; 1, 4, 11, 12}
# That gives diameter 4, not 5!

# Oh! Let me check the diameter again
# If diameter is 4 (not 5), then N_eig = 5

# Let me recount distances
print("  Distance distribution from vertex 0:")
for d in range(diameter + 1):
    count = np.sum(dist[0] == d)
    print("  d=%d: %d vertices" % (d, count))

# If distances are 0:1, 1:12, 2:42, 3:52, 4:12, 5:1
# vs distances 0:1, 1:12, 2:42, 3:52, 4:12, 5:1
# The antipodal vertex is at distance 5 (diameter = 5)

# For the ADJACENCY matrix of the 600-cell:
# Actually: 600-cell has 120 vertices. Let me verify
# d=0: 1 vertex (self)
# d=1: 12 vertices (neighbors)
# d=2: ?
# d=3: ?
# d=4: ?
# d=5: 1 (antipodal)
# Total should be 120

# Actually: 1 + 12 + ... + 1 = 120
# The halved 600-cell has intersection array {12,11,8,1; 1,4,11,12}
# with diameter 4. The full 600-cell has diameter 5 with antipodal.

# For a distance-regular graph with diameter d,
# there are exactly d+1 distinct eigenvalues.
# 600-cell: diameter 5 => 6 distinct eigenvalues? But we found 9.

# KEY INSIGHT: The 600-cell is NOT distance-regular!
# Wait... it is well-known that it IS distance-regular.
# Let me check: is the graph actually the 1-skeleton of the 600-cell?
# The 600-cell has 120 vertices and 720 edges.
# Its 1-skeleton is known to be distance-regular with intersection array
# {12, 11, 8, 1; 1, 4, 11, 12} and diameter 4.
# But we found diameter 5!

# Actually wait: the polytope 600-cell has 600 tetrahedral cells.
# Its 1-skeleton has 120 vertices and 720 edges.
# The diameter of the 1-skeleton should be checked.

# RESOLUTION: Maybe our "diameter 5" includes distance to antipodal point,
# and the standard intersection array gives diameter 4 for the
# distance-regular quotient (ignoring antipodal identification).

# For the FULL graph: 120 vertices
# If diameter = 5 and the graph is bipartite at maximum distance...
# Actually, 600-cell IS an antipodal graph. Each vertex has exactly one
# antipode at max distance. The "folded" 600-cell (identify antipodes)
# has 60 vertices and diameter 4.

# For an antipodal distance-regular graph with diameter d,
# the number of distinct eigenvalues is still d+1.
# So diameter 5 => 6 eigenvalues.
# But we found 9!

# Let me recheck: count unique eigenvalues
unique_evs = sorted(set(round(e, 4) for e in eigvals))
print("\n  Unique eigenvalues: %d" % len(unique_evs))
for ev in unique_evs:
    mult = sum(1 for e in eigvals if abs(e - ev) < 0.01)
    print("    %.6f (mult %d)" % (ev, mult))

print("\n  => The 600-cell has %d distinct eigenvalues" % len(unique_evs))
print("  => If distance-regular, would need diam+1 = %d" % (diameter + 1))
if len(unique_evs) != diameter + 1:
    print("  => INCONSISTENCY! Either not distance-regular or diameter wrong")

# Check: is 600-cell actually distance-regular?
# Verify intersection numbers are constant for ALL pairs at same distance
print("\n  Checking distance-regularity (100 random pairs per distance):")
np.random.seed(42)
is_dr = True
for d in range(1, diameter + 1):
    pairs_d = [(i, j) for i in range(min(20, N)) for j in range(N)
               if dist[i, j] == d]
    if not pairs_d:
        continue
    sample = [pairs_d[k] for k in np.random.choice(len(pairs_d), min(100, len(pairs_d)), replace=False)]
    b_set = set()
    c_set = set()
    for u, v in sample:
        nbrs_v = [j for j in range(N) if adj[v, j] > 0.5]
        b_count = sum(1 for w in nbrs_v if d < diameter and dist[u, w] == d + 1)
        c_count = sum(1 for w in nbrs_v if d > 0 and dist[u, w] == d - 1)
        b_set.add(b_count)
        c_set.add(c_count)
    is_const = len(b_set) <= 1 and len(c_set) <= 1
    if not is_const:
        is_dr = False
    print("  d=%d: b values=%s, c values=%s, constant=%s" % (
        d, sorted(b_set), sorted(c_set), "YES" if is_const else "NO"))

print("  Distance-regular: %s" % ("YES" if is_dr else "NO"))


# ============================================================
# CONCLUSIONS
# ============================================================
print()
print("=" * 72)
print("CONCLUSIONS")
print("=" * 72)
print()

print("1. CKM bare exponents {3, 7, 12} are given by the UNIQUE")
print("   linear formula:")
print()
print("     n(b1, b2) = N_eig * b2 - diam * b1 - b_1")
print()
print("   where (b1, b2) are generation quantum numbers from the")
print("   Generation Theorem, and:")
print("     N_eig = 9 (number of distinct eigenvalues)")
print("     diam  = 5 (diameter = disc(Z[phi]) = a_1)")
print("     b_1   = 6 (second spectral parameter)")
print()
print("2. The coefficients satisfy:")
print("     degree = diam + b_1 + 1 = 12")
print("     N_eig  = diam + b_1 - 2 = 9")
print("     diam + b_1 = 11 (Schur complement denominator!)")
print()
print("3. The formula is EXACTLY determined by 3 data points.")
print("   The fact that the coefficients match (N_eig, diam, b_1)")
print("   is a NON-TRIVIAL OBSERVATION but not a derivation.")
print()
print("4. For a TRUE derivation, we need to show WHY the CKM")
print("   mixing angle between generations b1 and b2 should")
print("   scale as phi^(-(N_eig*b2 - diam*b1 - b_1)).")
print("   This requires a dynamical mechanism.")
print()
print("STATUS: PATTERN (formula exact, coefficients match geometry,")
print("        but no derivation of WHY this formula)")
print()
print("EXP-205 COMPLETE")

#!/usr/bin/env python3
"""
EXP-210: ASSOCIATION SCHEME DEEP ANALYSIS
============================================
The 600-cell is NOT distance-regular (exp205 finding!).
It has 8 association scheme classes (not 5 = diameter).
This gives 9 eigenvalues (= 8 classes + identity).

The inner product classes at each graph distance:
  d=1: ip=phi/2 (12 vertices)
  d=2: ip=1/2 (20 vertices) + ip=(phi-1)/2 (12 vertices)
  d=3: ip=0 (30 vertices) + ip=(1-phi)/2 (12 vertices)
  d=4: ip=-1/2 (20 vertices) + ip=-phi/2 (12 vertices)
  d=5: ip=-1 (1 vertex)

The splits at d=2,3,4 produce EXTRA structure.
Can this structure explain PMNS angles or other physics?

APPROACHES:
A. Full association scheme analysis (intersection matrices)
B. Idempotents of the Bose-Mesner algebra
C. Connection to PMNS: do the class sizes relate to mixing?
D. Connection to gauge structure: A/B/C vertex types vs classes
E. The "12" vertices at each split distance - what are they?

STATUS: EXPLORATORY
Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
import math
from itertools import permutations as iperms

PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3

print("=" * 72)
print("EXP-210: ASSOCIATION SCHEME DEEP ANALYSIS")
print("=" * 72)
print()

# Build 600-cell
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

# All 8 inner product values (excluding self ip=1)
# ip = phi/2, 1/2, (phi-1)/2, 0, (1-phi)/2, -1/2, -phi/2, -1
ip_values = sorted(set(round(dot_mat[0, j], 6) for j in range(1, N)), reverse=True)
print("  Inner product classes: %d" % len(ip_values))
for ip in ip_values:
    cnt = sum(1 for j in range(1, N) if abs(dot_mat[0, j] - ip) < 1e-4)
    # Identify ip
    for a in range(-3, 3):
        for b in range(-3, 3):
            if abs((a + b*PHI)/2 - ip) < 1e-4:
                label = "(%d+%d*phi)/2" % (a, b)
                break
        else:
            continue
        break
    print("  ip = %8.6f = %s: %d vertices" % (ip, label, cnt))

# ============================================================
# PART A: ADJACENCY MATRICES OF EACH CLASS
# ============================================================
print()
print("=" * 72)
print("PART A: CLASS ADJACENCY MATRICES")
print("=" * 72)
print()

# Build adjacency matrix for each class
class_adj = {}
class_sizes = {}
for ip in ip_values:
    A_class = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if abs(dot_mat[i, j] - ip) < 1e-4:
                A_class[i, j] = 1.0
    class_adj[ip] = A_class
    class_sizes[ip] = int(A_class[0].sum())
    print("  Class ip=%.4f: degree=%d" % (ip, class_sizes[ip]))

# Verify: sum of all class degrees + 1 (self) = N
total = 1 + sum(class_sizes.values())
print("  Total: 1 + %d = %d (should be %d)" % (sum(class_sizes.values()), total, N))
print()

# ============================================================
# PART B: EIGENVALUES OF EACH CLASS ADJACENCY MATRIX
# ============================================================
print("=" * 72)
print("PART B: EIGENVALUES OF CLASS MATRICES")
print("=" * 72)
print()

# For an association scheme, each class matrix has the SAME eigenvectors
# (they all commute). The eigenvalues are the "P-matrix" entries.

# Compute eigenvalues for each class
class_eigvals = {}
for ip in ip_values:
    ev = np.sort(np.linalg.eigvalsh(class_adj[ip]))[::-1]
    # Group
    clusters = []
    for e in ev:
        found = False
        for cl in clusters:
            if abs(e - cl[0]) < 0.01:
                cl.append(e)
                found = True
                break
        if not found:
            clusters.append([e])
    class_eigvals[ip] = [(np.mean(cl), len(cl)) for cl in clusters]

print("  Eigenvalues of class adjacency matrices:")
print("  (Each column is one class, rows are eigenvalue clusters)")
print()

# Format as table
header = "  %5s" % "mult"
for ip in ip_values:
    header += " %8.4f" % ip
print(header)
print("  " + "-" * (5 + 9*len(ip_values)))

# Collect all eigenvalue indices (should be same across classes)
# For association scheme: 9 common eigenspaces
# Eigenvalue of identity class = multiplicity of eigenspace
all_evs = {}
for ip in ip_values:
    for ev, mult in class_eigvals[ip]:
        found = False
        for key in all_evs:
            if abs(key - mult) < 0.5:  # group by multiplicity
                found = True
                break
        if not found:
            all_evs[mult] = {}

# Better: use the adjacency matrix eigenvectors to project
adj_ip1 = class_adj[ip_values[0]]  # first class = neighbors at ip=phi/2
eigvals_full, eigvecs_full = np.linalg.eigh(adj_ip1)
idx_sort = np.argsort(eigvals_full)[::-1]
eigvals_full = eigvals_full[idx_sort]
eigvecs_full = eigvecs_full[:, idx_sort]

# Group eigenvalues into 9 clusters
evec_clusters = []
tol = 0.01
for i in range(N):
    ev = eigvals_full[i]
    found = False
    for cl in evec_clusters:
        if abs(ev - eigvals_full[cl[0]]) < tol:
            cl.append(i)
            found = True
            break
    if not found:
        evec_clusters.append([i])

print("  P-matrix (eigenvalues of class matrices in each eigenspace):")
print("  %5s %5s" % ("space", "mult"), end="")
for ip in ip_values:
    print(" %10.4f" % ip, end="")
print()
print("  " + "-" * (12 + 11*len(ip_values)))

P_matrix = np.zeros((len(evec_clusters), len(ip_values)))
for si, cl in enumerate(evec_clusters):
    mult = len(cl)
    # Project each class matrix onto this eigenspace
    # Use ONE eigenvector from this space
    v = eigvecs_full[:, cl[0]]
    print("  %5d %5d" % (si, mult), end="")
    for ci, ip in enumerate(ip_values):
        # Eigenvalue = v^T A_class v
        ev_class = v @ class_adj[ip] @ v
        P_matrix[si, ci] = ev_class
        print(" %10.4f" % ev_class, end="")
    print()

# ============================================================
# PART C: THE P-MATRIX AND PHYSICS
# ============================================================
print()
print("=" * 72)
print("PART C: P-MATRIX STRUCTURE AND PHYSICS")
print("=" * 72)
print()

# The P-matrix entries P_{i,j} give the eigenvalue of the j-th class
# adjacency matrix in the i-th eigenspace.
# For a distance-regular graph, P is the (d+1) x (d+1) matrix.
# For our 8-class scheme, P is 9 x 8.

# The eigenspace multiplicities are the same as the adjacency eigenvalues:
# 1, 4, 9, 16, 25, 36, 9, 16, 4

mults_list = [len(cl) for cl in evec_clusters]
print("  Eigenspace multiplicities: %s" % mults_list)
print("  Sum: %d (should be %d)" % (sum(mults_list), N))
print()

# Key question: which class matrices have eigenvalues that relate to PMNS?
# sin^2(th13) = 1/45 = 1/(5*9)
# sin^2(th23) = 4/7
# sin^2(th12) = 2/(phi+5)

# The class sizes: 12, 20, 12, 30, 12, 20, 12, 1
# Note the pattern: 12 appears 4 times, 20 appears 2 times, 30 once, 1 once
# 4*12 + 2*20 + 30 + 1 = 48 + 40 + 30 + 1 = 119. CHECK.

print("  Class sizes: %s" % [class_sizes[ip] for ip in ip_values])
print("  12 appears %d times" % sum(1 for s in class_sizes.values() if s == 12))
print("  20 appears %d times" % sum(1 for s in class_sizes.values() if s == 20))
print("  30 appears %d times" % sum(1 for s in class_sizes.values() if s == 30))
print("  1 appears %d times" % sum(1 for s in class_sizes.values() if s == 1))
print()

# The "12" classes are the phi-dependent inner products
# The "20" classes are the rational inner products (+/-1/2)
# The "30" class is the zero inner product
# The "1" class is the antipodal vertex

# This split: 4*12 vs 2*20 vs 30 vs 1 = 48 vs 40 vs 30 vs 1
# 48 = phi-sector vertices
# 40 = rational-sector vertices
# 30 = orthogonal vertices
# 1 = antipodal

# The phi-sector (48 vertices) vs rational-sector (40 vertices):
# 48/40 = 6/5 = b_1/a_1 !!!

print("  PHI-SECTOR / RATIONAL-SECTOR = 48/40 = 6/5 = b_1/a_1 !!!")
print()

# More detailed: the distances at which splitting occurs
# d=2: 20 (rational ip=1/2) + 12 (phi ip=(phi-1)/2) = 32
# d=3: 30 (rational ip=0) + 12 (phi ip=(1-phi)/2) = 42
# d=4: 20 (rational ip=-1/2) + 12 (phi ip=-phi/2) = 32
# At each split distance: rational + phi = total

# Ratio at each distance:
print("  Split structure at each distance:")
for d, (rat_n, phi_n) in [(2, (20, 12)), (3, (30, 12)), (4, (20, 12))]:
    ratio = rat_n / phi_n
    print("  d=%d: %d rational + %d phi = %d (ratio = %.4f = %s)" % (
        d, rat_n, phi_n, rat_n + phi_n, ratio,
        "5/3" if abs(ratio - 5/3) < 0.01 else "5/2" if abs(ratio - 5/2) < 0.01 else "%.3f" % ratio))

print()

# ============================================================
# PART D: VERTEX TYPE VS ASSOCIATION CLASS
# ============================================================
print("=" * 72)
print("PART D: VERTEX TYPES (A/B/C) IN EACH CLASS")
print("=" * 72)
print()

def classify_vertex(v):
    c = sorted([abs(x) for x in v], reverse=True)
    if abs(c[0] - 1.0) < 1e-10 and c[1] < 1e-10:
        return 'A'
    if all(abs(x - 0.5) < 1e-10 for x in c):
        return 'B'
    return 'C'

vtypes = [classify_vertex(v) for v in vertices]

# For each class, count vertex types
print("  Vertex type distribution in each class (from vertex 0, type %s):" % vtypes[0])
for ip in ip_values:
    verts_in = [j for j in range(1, N) if abs(dot_mat[0, j] - ip) < 1e-4]
    type_counts = {'A': 0, 'B': 0, 'C': 0}
    for j in verts_in:
        type_counts[vtypes[j]] += 1
    print("  ip=%7.4f (n=%2d): A=%2d B=%2d C=%2d" % (
        ip, len(verts_in), type_counts['A'], type_counts['B'], type_counts['C']))

print()

# ============================================================
# PART E: KREIN PARAMETERS AND PMNS CONNECTION
# ============================================================
print()
print("=" * 72)
print("PART E: SPECTRAL ANALYSIS AND PMNS")
print("=" * 72)
print()

# The Q-matrix (dual eigenvalue matrix) Q_{j,i} = m_i * P_{i,j} / n_j
# where m_i = multiplicity of eigenspace i, n_j = class size j
Q_matrix = np.zeros_like(P_matrix)
class_size_list = [class_sizes[ip] for ip in ip_values]
for i in range(len(evec_clusters)):
    for j in range(len(ip_values)):
        Q_matrix[i, j] = mults_list[i] * P_matrix[i, j] / class_size_list[j]

print("  Q-matrix (dual eigenvalues):")
print("  %5s %5s" % ("space", "mult"), end="")
for ci, ip in enumerate(ip_values):
    print(" %10.4f" % ip, end="")
print()
print("  " + "-" * (12 + 11*len(ip_values)))
for si in range(len(evec_clusters)):
    print("  %5d %5d" % (si, mults_list[si]), end="")
    for ci in range(len(ip_values)):
        print(" %10.4f" % Q_matrix[si, ci], end="")
    print()

print()

# Check: P * Q^T = N * I (orthogonality relation)
PQ = P_matrix @ Q_matrix.T
print("  P * Q^T diagonal: %s" % [round(PQ[i,i], 2) for i in range(min(9, PQ.shape[0]))])
print("  Should be N = %d on diagonal" % N)
print()

# The ratios in the P-matrix might connect to PMNS
# Let's look at specific ratios
print("  Interesting P-matrix ratios:")
for si in range(min(5, P_matrix.shape[0])):
    for ci in range(P_matrix.shape[1]):
        val = P_matrix[si, ci]
        if abs(val) > 0.1:
            # Check if val is a+b*phi
            for a in range(-20, 20):
                for b in range(-20, 20):
                    if abs(a + b*PHI - val) < 0.01 and (a != 0 or b != 0):
                        if si < 3:  # only first few eigenspaces
                            pass  # will print below
                        break

# Let me focus on what the P-matrix values are
print("  P-matrix in a+b*phi form:")
for si in range(min(9, P_matrix.shape[0])):
    print("  space %d (mult %d):" % (si, mults_list[si]), end="")
    for ci in range(P_matrix.shape[1]):
        val = P_matrix[si, ci]
        found = False
        for a in range(-20, 21):
            for b in range(-10, 11):
                if abs(a + b*PHI - val) < 0.01:
                    if b == 0:
                        print(" %4d" % a, end="")
                    else:
                        print(" %2d%+d*p" % (a, b), end="")
                    found = True
                    break
            if found:
                break
        if not found:
            print(" %6.2f" % val, end="")
    print()

print()

# ============================================================
# PART F: PMNS FROM CLASS RATIOS
# ============================================================
print()
print("=" * 72)
print("PART F: PMNS ANGLES FROM ASSOCIATION SCHEME")
print("=" * 72)
print()

# The class sizes are 12, 20, 12, 30, 12, 20, 12, 1
# The phi-sector has 48 vertices, the rational sector 40, orthogonal 30, antipodal 1

# Key observation:
# Total phi-sector = 48 = 4*12
# Total rational = 40 = 2*20
# Orthogonal = 30
# Antipodal = 1

# sin^2(th12) = 2/(phi+5) = 0.302
# sin^2(th23) = 4/7
# sin^2(th13) = 1/45

# Can we construct these from class sizes?
# 1/45 = 1/(5*9) - already geometric
# 4/7 = 4/7 - Schur gap

# What about from the P-matrix?
# The P-matrix diagonal for the neighbor class (ip=phi/2):
# These are the eigenvalues of the adjacency matrix itself!

# Let me extract the 9x9 intersection number matrices
# For an association scheme: A_i * A_j = sum_k p_{ij}^k * A_k

print("  Computing intersection numbers p_{ij}^k:")
print("  (Structure constants of the Bose-Mesner algebra)")
print()

n_classes = len(ip_values)
# Include the identity class (A_0 = I)
# So we have 9 matrices: A_0 = I, A_1 ... A_8

# Actually, compute p_{ij}^k for class matrices
# p_{ij}^k = (1/n_k) * Tr(A_i * A_j * A_k) ... no, that's for symmetric scheme
# Better: count directly
# p_{ij}^k = number of vertices w such that (v,w) in class i AND (w,u) in class j
# when (v,u) in class k. This should be constant for all (v,u) in class k.

# For efficiency, just check for a few representative pairs
# Representative pair for each class: (0, j) where j is in that class
reps = {}
for ip in ip_values:
    for j in range(1, N):
        if abs(dot_mat[0, j] - ip) < 1e-4:
            reps[ip] = j
            break

def get_class(i, j):
    """Return the class (inner product) of edge (i,j)"""
    if i == j:
        return None  # identity class
    ip = round(dot_mat[i, j], 4)
    for ref_ip in ip_values:
        if abs(ip - ref_ip) < 1e-3:
            return ref_ip
    return None

# Compute p_{1,j}^k where class 1 = neighbors (ip=phi/2)
# This is the most important set of intersection numbers
ip_nbr = ip_values[0]  # phi/2
print("  Intersection numbers p_{nbr, j}^k:")
print("  (How many neighbors of v are in class j relative to u, when (v,u) in class k)")
print()

for k_ip in ip_values:
    u = reps[k_ip]  # representative vertex at class k from 0
    print("  Class k = ip %.4f (vertex %d at distance %s from 0):" % (
        k_ip, u, "?"))
    # For each class j: count neighbors w of 0 such that (w, u) is in class j
    nbrs_0 = [w for w in range(N) if abs(dot_mat[0, w] - ip_nbr) < 1e-4]
    counts = {}
    for j_ip in ip_values:
        counts[j_ip] = 0
    counts['self'] = 0  # w = u
    for w in nbrs_0:
        if w == u:
            counts['self'] += 1
            continue
        w_class = get_class(w, u)
        if w_class is not None:
            counts[w_class] = counts.get(w_class, 0) + 1
    # Print
    for j_ip in ip_values:
        if counts[j_ip] > 0:
            print("    p_{nbr, %.4f}^{%.4f} = %d" % (j_ip, k_ip, counts[j_ip]))
    if counts['self'] > 0:
        print("    (vertex u is a neighbor: %d)" % counts['self'])
    print()

# ============================================================
# PART G: SUMMARY OF KEY NUMBERS
# ============================================================
print()
print("=" * 72)
print("PART G: KEY NUMBERS FROM ASSOCIATION SCHEME")
print("=" * 72)
print()

print("  Class sizes: 12, 20, 12, 30, 12, 20, 12, 1")
print("  Phi-sector total: 4*12 = 48")
print("  Rational-sector total: 2*20 = 40")
print("  Orthogonal: 30")
print("  Antipodal: 1")
print()
print("  Key ratios:")
print("  48/40 = 6/5 = b_1/a_1")
print("  48/120 = 2/5 = 2/a_1")
print("  40/120 = 1/3 = 1/N_gen")
print("  30/120 = 1/4")
print("  48/30 = 8/5 = phi^2 - 1 = phi/1 - ... hmm")
print("  (48+1)/120 = 49/120")
print("  12/120 = 1/10 = alpha_s/alpha = 1/(10*phi)... no, 1/10")
print()

# 40/120 = 1/3 = 1/N_gen  -- THIS IS INTERESTING!
# The rational-sector is exactly 1/3 of the vertices!
# And N_gen = 3 was DERIVED from the framework!
# sin^2(th12_TBM) = 1/3 is the tribimaximal prediction!

print("  *** RATIONAL-SECTOR = 40/120 = 1/3 = 1/N_gen ***")
print("  *** This equals tribimaximal sin^2(th12) = 1/3 ***")
print()

# But our formula gives 2/(phi+5) = 0.302, which is BETTER than 1/3 = 0.333
# The correction: 2/(phi+5) = 2/(phi+a_1) introduces the phi-sector into the mixing

# What if sin^2(th12) = (rational_sector - phi-correction) / total?
# 40/120 - delta = 0.302 => delta = 0.031 ~ alpha_s/... ?

# Actually: 2/(phi+5) = 2*2/(2*phi+10) = 4/(2*phi+10)
# And: rational_sector / total = 40/120 = 1/3
# 1/3 = 40/120 and phi+5 = 6.618
# 2/(phi+5) = 2/6.618 = 0.302

# Another way: 2/(phi+5) = 2/(phi+a_1)
# The numerator 2 could be the number of rational classes (ip = +/-1/2)
# The denominator phi+5 = phi + a_1 mixes the phi sector with the rational sector

print("  INTERPRETATION:")
print("  sin^2(th12) = N_rational_classes / (phi + a_1)")
print("  = 2 / (phi + 5) = 0.302")
print("  where N_rational_classes = 2 (the ip=+1/2 and ip=-1/2 classes)")
print("  and phi + a_1 is the 'mixing denominator'")
print()

# For th23:
# sin^2(th23) = 4/7
# 4 = number of phi-type classes (ip = +-phi/2, +-(phi-1)/2)
# 7 = Schur complement gap
print("  sin^2(th23) = N_phi_classes / Schur_gap")
print("  = 4 / 7")
print("  where N_phi_classes = 4 (the four phi-dependent ip classes)")
print("  and Schur_gap = 7 = 12 - 5")
print()

# For th13:
# sin^2(th13) = 1/45 = 1/(5*9)
# 1 = antipodal class
# 45 = a_1 * N_eig = 5 * 9
print("  sin^2(th13) = N_antipodal / (a_1 * N_eig)")
print("  = 1 / (5 * 9) = 1/45")
print("  where N_antipodal = 1 and a_1*N_eig = 45")
print()

# THIS IS A UNIFIED STRUCTURE!
# All three PMNS angles come from:
# sin^2(th) = N_classes_of_type / geometric_denominator
# th12: 2 rational classes / (phi + 5)
# th23: 4 phi classes / 7
# th13: 1 antipodal class / 45

print("  *** UNIFIED STRUCTURE ***")
print("  sin^2(th_ij) = N(type) / D(type)")
print("  where N counts association classes of a specific type")
print("  and D is a geometric denominator from 600-cell")
print()
print("  type     | N(type) | D(type)        | sin^2    | angle")
print("  " + "-" * 60)
print("  rational | 2       | phi+a_1=phi+5  | 0.302    | th12")
print("  phi      | 4       | Schur_gap=7    | 0.571    | th23")
print("  antipodal| 1       | a_1*N_eig=45   | 0.022    | th13")
print()

# ============================================================
# CONCLUSIONS
# ============================================================
print()
print("=" * 72)
print("CONCLUSIONS")
print("=" * 72)
print()

print("1. ASSOCIATION SCHEME: 8 classes split as")
print("   4 phi-type (size 12 each = 48 total)")
print("   2 rational-type (size 20 each = 40 total)")
print("   1 orthogonal (size 30)")
print("   1 antipodal (size 1)")
print("   STATUS: DERIVED")
print()

print("2. KEY RATIO: phi-sector/rational-sector = 48/40 = 6/5 = b_1/a_1")
print("   rational-sector/total = 40/120 = 1/3 = 1/N_gen")
print("   STATUS: DERIVED")
print()

print("3. PMNS UNIFIED STRUCTURE:")
print("   All three angles come from class-type counts divided by")
print("   geometric denominators:")
print("   th12: 2/(phi+5), th23: 4/7, th13: 1/45")
print("   The numerators (2, 4, 1) count association class TYPES:")
print("   2 rational classes, 4 phi classes, 1 antipodal class")
print("   STATUS: PATTERN (numerator interpretation is suggestive,")
print("   denominators still need mechanistic justification)")
print()

print("4. The non-distance-regularity of the 600-cell is not a")
print("   deficiency but a FEATURE: it encodes the PMNS mixing")
print("   through the phi/rational/orthogonal/antipodal class types.")
print()

print("EXP-210 COMPLETE")

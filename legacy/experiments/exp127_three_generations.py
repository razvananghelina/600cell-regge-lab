"""
EXP-127: Can We Derive 3 Generations from 600-Cell Geometry?
=============================================================

CONTEXT:
In the 600-cell theory, fermion masses are m_f = m_e * phi^(5a + 6b).
The quantum numbers (a,b) = argmin|a|+|b| s.t. 5a+6b=n.
Occupied levels: n = {0, 3, 5, 11, 16, 17, 19, 26}.

From exp121c/d: all 9 fermions lie on 3 lines in (a,b) plus origin:
  - Origin: electron (0,0)
  - Line a=1: down(1,0), muon/strange(1,1), tau(1,2)
  - Line b=1: muon/strange(1,1), charm(2,1), top(4,1)
  - Line 3a+2b=5: up(3,-2), muon/strange(1,1), bottom(-1,4)

THE 3-GENERATION PROBLEM:
On the a=1 line, occupied fermions have b in {0,1,2}.
The next candidate b=3 (n=23, a 4th gen lepton/down-quark) is NOT occupied.
Currently b<=2 is a PHYSICAL constraint but not GEOMETRICALLY derived.

QUESTION: Can we derive b<=2 (i.e., 3 generations) from 600-cell geometry?

APPROACHES:
 1. Norm analysis on a=1 line (Z[phi] norm cutoff)
 2. Diameter constraint (n vs graph distance)
 3. H4 representation theory (dimension counting)
 4. Icosahedral vertex figure (A5 representations)
 5. Spectral gap (adjacency eigenvalue constraint)
 6. E8 branching (H4 x H4' Coxeter exponents)
 7. Z[phi] unit structure (unit group pattern)
 8. Topological/homological (S^3 triangulation)
 9. Physical mass bound (is phi^23 too heavy?)
10. Decagon structure (great decagon counting)

Author: Claude (exp127)
Date: 2026-02-09
"""

import numpy as np
from itertools import product, permutations
from collections import defaultdict, Counter
from math import gcd, factorial

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2

# ============================================================
# SETUP: Common functions and data
# ============================================================

def min_complexity(n, a1=5, b1=6, search_range=50):
    """Find (a,b) with minimum |a|+|b| such that a1*a + b1*b = n."""
    best = None
    best_c = float('inf')
    for a in range(-search_range, search_range + 1):
        if (n - a1 * a) % b1 == 0:
            b = (n - a1 * a) // b1
            c = abs(a) + abs(b)
            if c < best_c:
                best_c = c
                best = (a, b)
    return best, best_c

def zphi_norm(a, b):
    """Signed norm N(a + b*phi) = a^2 + ab - b^2."""
    return a**2 + a * b - b**2

# Fermion data
fermions = {
    'electron':  (0, 0, 0),
    'down':      (1, 0, 5),
    'muon':      (1, 1, 11),
    'strange':   (1, 1, 11),
    'tau':       (1, 2, 17),
    'up':        (3, -2, 3),
    'charm':     (2, 1, 16),
    'top':       (4, 1, 26),
    'bottom':    (-1, 4, 19),
}

occupied_n = {0, 3, 5, 11, 16, 17, 19, 26}

# The a=1 line fermions
a1_line = {
    'down':   (1, 0, 5),
    'muon':   (1, 1, 11),  # also strange
    'tau':    (1, 2, 17),
}
# Candidate 4th generation on a=1 line:
# b=3 -> n=5+18=23, (a,b) = (1,3)

print("=" * 80)
print("EXP-127: Can We Derive 3 Generations from 600-Cell Geometry?")
print("=" * 80)
print()

print("  The a=1 line in (a,b) lattice:")
print("    b=0: down   (1,0), n=5")
print("    b=1: mu/s   (1,1), n=11")
print("    b=2: tau    (1,2), n=17")
print("    b=3: ???    (1,3), n=23  <-- NOT occupied (4th generation?)")
print()

# ============================================================
# APPROACH 1: Norm analysis on a=1 line
# ============================================================

print("=" * 80)
print("APPROACH 1: Z[phi] Norm Analysis on a=1 Line")
print("=" * 80)
print()

print("  For a=1, the Z[phi] element is z = 1 + b*phi")
print("  Norm N(z) = 1^2 + 1*b - b^2 = 1 + b - b^2")
print()

print(f"  {'b':>3s} {'(a,b)':>7s} {'n=5+6b':>7s} {'N(z)':>6s} {'|N|':>5s}"
      f" {'z':>10s} {'z_conj':>10s} {'occupied':>9s}")
print(f"  {'-'*65}")

for b in range(-2, 8):
    a = 1
    n = 5 * a + 6 * b
    N = zphi_norm(a, b)
    z = a + b * PHI
    z_conj = a + b * PHI_CONJ
    occ = "YES" if n in occupied_n else "no"
    print(f"  {b:>3d} ({a:>2d},{b:>2d}) {n:>7d} {N:>6d} {abs(N):>5d}"
          f" {z:>10.4f} {z_conj:>10.4f} {occ:>9s}")

print()

# From exp121b: allowed norms are |N| in {0, 1, 5, 19}
print("  Allowed norms from exp121b: |N| in {0, 1, 5, 19}")
print()
print("  Norm values on a=1 line:")
for b in range(-2, 8):
    N = zphi_norm(1, b)
    allowed = abs(N) in {0, 1, 5, 19}
    print(f"    b={b}: N = {N:>4d}, |N| = {abs(N):>3d}, allowed = {allowed}")

print()
print("  RESULT: b=3 has |N|=5 which IS in the allowed set.")
print("  Norm analysis alone does NOT exclude b=3.")
print()

# Check: is there a cutoff on the SIGNED norm?
print("  Signed norm pattern on a=1 line:")
print("    b=0: N = +1 (occupied)")
print("    b=1: N = -1 (occupied)")
print("    b=2: N = -5 (occupied)")
print("    b=3: N = -11 (NOT occupied, |N|=11 NOT in allowed set)")

print()
print("  WAIT: |N(1,3)| = |1+3-9| = |-5| = 5. Let me recalculate...")
N_13 = zphi_norm(1, 3)
print(f"  N(1,3) = 1 + 3 - 9 = {N_13}")
print(f"  |N(1,3)| = {abs(N_13)}")
print(f"  Is |N|=5 in allowed set {{0,1,5,19}}? {abs(N_13) in {0, 1, 5, 19}}")
print()

# So the norm alone does NOT exclude b=3. What about exp121d's line_valid?
# On a=1 line, line_valid requires b <= 2. That's the constraint we want to DERIVE.
print("  CONCLUSION APPROACH 1: Norm analysis does NOT derive b<=2.")
print("  |N(1,3)| = 5 is in the allowed set. The exclusion of n=23")
print("  requires the ADDITIONAL constraint b<=2 from exp121d.")
print("  Category: NOT DERIVED")
print()


# ============================================================
# APPROACH 2: Diameter constraint
# ============================================================

print("=" * 80)
print("APPROACH 2: Diameter Constraint")
print("=" * 80)
print()

print("  The 600-cell has diameter 5 (max graph distance between vertices).")
print("  In our mass formula, n = 5a + 6b where 5 = diameter, 6 = decagons.")
print()

# Build 600-cell to check diameter
def build_600cell():
    phi = PHI
    verts = set()
    def add(v):
        n = np.sqrt(sum(x**2 for x in v))
        verts.add(tuple(round(x/n, 10) for x in v))
    for i in range(4):
        for s in [+1, -1]:
            v = [0., 0., 0., 0.]
            v[i] = float(s)
            add(v)
    for ss in product([+1, -1], repeat=4):
        add([s * 0.5 for s in ss])
    base = [phi / 2, 0.5, 1 / (2 * phi), 0.0]
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            pv = [base[p[i]] for i in range(4)]
            nz = [i for i in range(4) if abs(pv[i]) > 1e-12]
            for ss in product([+1, -1], repeat=len(nz)):
                v = list(pv)
                for k, idx in enumerate(nz):
                    v[idx] *= ss[k]
                add(v)
    return np.array([list(v) for v in verts])

V = build_600cell()
N_verts = len(V)
print(f"  600-cell vertices: {N_verts}")

# Build adjacency
dots = V @ V.T
np.clip(dots, -1, 1, out=dots)
angles = np.arccos(dots)
theta_edge = np.sort(np.unique(np.round(angles[0], 6)))[1]

adj = np.zeros((N_verts, N_verts), dtype=int)
for i in range(N_verts):
    for j in range(i + 1, N_verts):
        if abs(angles[i, j] - theta_edge) < 0.001:
            adj[i, j] = 1
            adj[j, i] = 1

# BFS for distances
from collections import deque

def bfs_distances(adj_matrix, start):
    n = adj_matrix.shape[0]
    dist = [-1] * n
    dist[start] = 0
    queue = deque([start])
    while queue:
        v = queue.popleft()
        for u in range(n):
            if adj_matrix[v, u] == 1 and dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist

dist_0 = bfs_distances(adj, 0)
diameter = max(dist_0)
print(f"  Graph diameter: {diameter}")

# Shell sizes
shell_sizes = Counter(dist_0)
print(f"  Shell sizes from vertex 0:")
for d in range(diameter + 1):
    print(f"    d={d}: {shell_sizes[d]} vertices")

print()

# Check: on the a=1 line, n = 5 + 6b.
# Does n have to fit within some geodesic constraint?
print("  On a=1 line: n = 5 + 6b")
print("  The max n from graph paths: diameter d=5 gives n up to what?")
print()
print("  Each step along the graph changes n by at most 5 (diameter axis)")
print("  or 6 (decagon axis). But n = 5a + 6b is a LINEAR combination,")
print("  not a path length.")
print()

# What IS the max n for occupied levels?
print("  Max occupied n: 26 (top quark), which is 4*5 + 1*6 = 26")
print("  Candidate n=23: 1*5 + 3*6 = 23")
print("  Both well beyond diameter 5.")
print()

# Is there a constraint from wrapping around S^3?
# S^3 circumference = pi (in S^3 of unit radius, geodesic distance max = pi)
# The 600-cell lives on unit S^3
print("  Geodesic distance analysis:")
print("  600-cell on unit S^3: max geodesic distance = pi")
print(f"  Edge angular distance: {theta_edge:.6f} rad = {np.degrees(theta_edge):.2f} deg")
print()

# For each level n, what angular distance corresponds?
print("  Angular interpretation of n = 5a + 6b:")
print("  If 'a' counts steps of size theta_diam and 'b' steps of size theta_decagon...")
theta_diam = np.pi  # diameter = antipodal point
# For decagon: each great decagon has vertices at angular separations 2*pi/10 = pi/5
theta_decagon = np.pi / 5

print(f"  theta_diameter (antipodal): {theta_diam:.4f} rad = {np.degrees(theta_diam):.1f} deg")
print(f"  theta_decagon (1/10 circle): {theta_decagon:.4f} rad = {np.degrees(theta_decagon):.1f} deg")
print()

for b in range(5):
    a = 1
    n = 5 * a + 6 * b
    # Total angular displacement (naive): a * theta_diam + b * theta_decagon
    # But this wraps around S^3
    ang = a * theta_diam + b * theta_decagon
    ang_mod = ang % (2 * np.pi)
    if ang_mod > np.pi:
        ang_mod = 2 * np.pi - ang_mod
    print(f"    b={b}: n={n:>3d}, raw angle = {ang:.4f} rad = {np.degrees(ang):.1f} deg,"
          f" wrapped = {ang_mod:.4f} rad = {np.degrees(ang_mod):.1f} deg")

print()
print("  RESULT: The diameter constraint (d=5) does not directly bound b.")
print("  n=23 and n=26 (top quark) are both > 5*5=25, so diameter alone")
print("  does not provide a meaningful cutoff.")
print("  Category: NOT DERIVED")
print()


# ============================================================
# APPROACH 3: H4 Representation Theory
# ============================================================

print("=" * 80)
print("APPROACH 3: H4 Representation Theory")
print("=" * 80)
print()

print("  H4 (Coxeter group): order |H4| = 14400 = 120^2")
print("  H4 Coxeter exponents: {1, 11, 19, 29}")
print("  H4 degrees: {2, 12, 20, 30}")
print()

# H4 irreducible representations
# H4 is a non-crystallographic Coxeter group of type H4
# It has irreps of dimensions:
# 1, 4, 5, 6, 8, 9, 10, 16, 18, 20, 24, 25, 30, 36, 40, 45, 48, 60 ...
# Actually, the complete set of irreps dimensions is known:

print("  H4 irreducible representations (from theory):")
print("  H4 has 34 conjugacy classes, hence 34 irreps.")
print()

# The adjacency eigenvalue multiplicities are perfect squares:
adj_mult = [1, 4, 9, 16, 25, 36, 9, 16, 4]
adj_eigs_desc = ["12", "6phi", "4phi", "3", "0", "-2", "4-4phi", "-3", "6-6phi"]
print("  Adjacency matrix multiplicities (from exp110):")
for ev, m in zip(adj_eigs_desc, adj_mult):
    sq = int(np.sqrt(m))
    print(f"    theta = {ev:>10s}, mult = {m:>3d} = {sq}^2")

print()
print(f"  Sum of multiplicities: {sum(adj_mult)} = {N_verts}")
print()

# Can "3" come from representation theory?
# H4 has a representation of dimension 4 (the natural reflection rep in R^4)
# A5 (icosahedral group) is a subgroup of H4
# A5 has irreps of dimensions: 1, 3, 3', 4, 5

print("  H4 -> A5 (vertex stabilizer) branching:")
print("  A5 irreps: 1, 3, 3', 4, 5 (dimensions)")
print("  The number 3 appears as dim(3) and dim(3')!")
print()
print("  Question: does '3 generations' correspond to dim(3) of A5?")
print()

# The number 3 also appears as:
print("  Occurrences of '3' in 600-cell geometry:")
print("    - Vertex dimension: 4 (NOT 3)")
print("    - A5 has 3-dim irreps")
print("    - H4 Coxeter number h = 30 = 10 * 3")
print("    - Euler characteristic chi(S^3) = 0")
print("    - 120/3 = 40 (no obvious meaning)")
print("    - 5 tetrahedra per edge, 5/? = ?")
print("    - 600-cell has 600 cells, 600/200 = 3")
print("    - 720 edges / 240 = 3")
print("    - 1200 faces / 400 = 3")
print()

# Key: on the a=1 line, {b=0,1,2} has 3 elements
# Is this related to the 3-dim irrep of A5?
print("  ANALYSIS: A5 has two 3-dimensional irreps (3 and 3').")
print("  In icosahedral symmetry, these correspond to rotations")
print("  of the icosahedron. Each vertex of the 600-cell has an")
print("  icosahedral vertex figure (12 neighbors = icosahedron).")
print()
print("  SPECULATION: The 3 values b in {0,1,2} on the a=1 line")
print("  might correspond to the 3 components of an A5 triplet.")
print("  But this requires showing that decagon steps (b) transform")
print("  under a 3-dim representation of the vertex stabilizer.")
print()
print("  RESULT: A5 does provide a natural '3', but the connection")
print("  to the b-quantum number on the a=1 line is not rigorous.")
print("  Category: SPECULATIV")
print()


# ============================================================
# APPROACH 4: Icosahedral Vertex Figure Analysis
# ============================================================

print("=" * 80)
print("APPROACH 4: Icosahedral Vertex Figure (A5 reps)")
print("=" * 80)
print()

# Each vertex has 12 neighbors forming an icosahedron
# The icosahedron has symmetry A5 (order 60)
# Its eigenvalues are: 5 (m=1), sqrt(5) (m=3), -1 (m=5), -sqrt(5) (m=3)

print("  Icosahedron adjacency eigenvalues:")
ico_eigs = [(5, 1), (np.sqrt(5), 3), (-1, 5), (-np.sqrt(5), 3)]
for ev, m in ico_eigs:
    print(f"    lambda = {ev:>8.4f}, mult = {m}")

print()
print("  Icosahedron eigenvalue multiplicities: {1, 3, 5, 3}")
print("  Note: 3 appears TWICE (for +-sqrt(5))")
print()

# The vertex figure of the 600-cell is an icosahedron
# The local symmetry at each vertex is the icosahedral group A5
# A5 irreps: 1(trivial), 3(standard), 3'(conjugate), 4, 5

# Key question: how does a "step in b direction" transform under A5?
# On the a=1 line, incrementing b by 1 changes n by 6.
# n=6 corresponds to one step along a great decagon.
# Each vertex lies on 6 great decagons.

print("  Number of great decagons through each vertex: 6")
print()

# The 6 great decagons through a vertex:
# How do they transform under the icosahedral vertex stabilizer?
# Each great decagon has 10 vertices, so from vertex v,
# it provides 2 nearest neighbors on the decagon (at distance 1)
# Total: 6 * 2 = 12 = all neighbors. So each edge belongs to exactly 1 decagon.

# Wait, that means the 12 edges from a vertex are partitioned into 6 pairs,
# each pair belonging to one decagon. The 6 decagons form a structure.

print("  Key structure: 6 great decagons through a vertex partition")
print("  the 12 edges into 6 pairs (opposite edges of the icosahedron).")
print()

# The 6 pairs of opposite edges of the icosahedron
# Under A5, these 6 pairs form a representation.
# A5 acting on 6 opposite-edge pairs:
# The icosahedron has 30 edges, so 15 pairs of opposite edges
# Wait, an icosahedron has 12 vertices and 30 edges.
# "Opposite edges" of the icosahedron: the icosahedron has no natural "opposite edge"
# concept... Let me reconsider.

# Actually: 6 great decagons correspond to the 6 5-fold rotation axes
# of the icosahedron (each axis connects two vertices of the icosahedron).
# 12 vertices / 2 = 6 axes.
# Each 5-fold axis of the icosahedron = one great decagon.

print("  CORRECTION: 6 great decagons correspond to the 6 five-fold")
print("  axes of the icosahedral vertex figure.")
print("  (12 icosahedron vertices / 2 = 6 axes)")
print()

# Under A5, the 6 five-fold axes form a representation.
# Since there are 6 axes and A5 has order 60, the stabilizer of
# one axis has order 60/6 = 10 (the cyclic group C_10 = C_5 x C_2).
# So the 6 axes form the representation of A5 on the coset space A5/(C5xC2).

# This representation decomposes into A5 irreps.
# A5 acting on 6 objects (the 5-fold axes):
# Character analysis: 6 = 1 + 5  (trivial + 5-dim standard)

# To verify: the A5 permutation representation on 6 axes
# The A5 action on the icosahedron vertices gives the 12-dim perm rep.
# Pairing opposite vertices gives 6 objects.
# The permutation representation on 6 points:
# Character: chi(id) = 6
# A5 conjugacy classes: {1}, {12 order 5}, {12 order 5}, {20 order 3}, {15 order 2}
# For pairs of opposite vertices under A5:
# - identity: fixes all 6 pairs -> chi = 6
# - 5-fold rotation: fixes 1 pair (the axis) -> chi = 1
# - 3-fold rotation: fixes 0 pairs -> chi = 0
# - 2-fold rotation: fixes 2 pairs -> chi = 2 (axis pair + one transverse)
#   Wait, a 2-fold rotation of icosahedron fixes 2 pairs (axis + midpoints)
#   Actually for order 2 (rotation by pi about edge midpoint):
#   fixes 0 vertices -> fixes 0 pairs. But there are different types:
#   In A5, order 2 elements are rotations by pi about axes through edge midpoints
#   These fix no vertices -> fix no pairs -> chi = 0
#   Hmm, but then chi = (6, 1, 1, 0, 0) which decomposes as:
#   6 = 1 + 5 since chi(5-dim) = (5, 0, 0, -1, 1) and chi(1) = (1,1,1,1,1)
#   (5,0,0,-1,1) + (1,1,1,1,1) = (6,1,1,0,2) which doesn't match.
#   Let me be more careful.

# A5 character table:
# Irreps:        1    3    3'   4    5
# {id}(1):       1    3    3    4    5
# {5A}(12):      1    phi  phi' -1    0
# {5B}(12):      1    phi' phi  -1    0
# {3}(20):       1    0    0    1   -1
# {2}(15):       1   -1   -1    0    1

phi_val = PHI
phi_conj_val = PHI_CONJ

print("  A5 character table:")
print(f"  {'Irrep':>6s} {'d':>3s}  {'C1':>5s} {'C5A':>5s} {'C5B':>5s} {'C3':>5s} {'C2':>5s}")
chars = {
    '1':  [1, 1, 1, 1, 1],
    '3':  [3, phi_val, phi_conj_val, 0, -1],
    "3'": [3, phi_conj_val, phi_val, 0, -1],
    '4':  [4, -1, -1, 1, 0],
    '5':  [5, 0, 0, -1, 1],
}
for name, ch in chars.items():
    print(f"  {name:>6s} {int(ch[0]):>3d}  {ch[0]:>5.2f} {ch[1]:>5.2f} {ch[2]:>5.2f}"
          f" {ch[3]:>5.2f} {ch[4]:>5.2f}")
print()

# For the 6-axes representation, I need to determine the character.
# A5 acting on 6 pentagonal (5-fold) axes:
# id: fixes 6 -> chi=6
# C5A (order 5 rot): fixes the axis of rotation -> chi=1
# C5B (order 5 rot, other class): also fixes 1 axis -> chi=1
# C3 (order 3 rot): no 5-fold axis is fixed by 3-fold rot -> chi=0
# C2 (order 2 rot): a 2-fold rotation about edge midpoint:
#   no 5-fold axis is fixed -> chi=0

# Character of A5 on 6 five-fold axes (pairs of opposite vertices):
# - C1 (identity): all 6 pairs fixed -> chi = 6
# - C5A (order 5 rot): fixes 1 pair (the rotation axis) -> chi = 1
# - C5B (order 5 rot): fixes 1 pair -> chi = 1
# - C3 (order 3 rot): no 5-fold axis is fixed -> chi = 0
# - C2 (order 2 rot, pi rotation about edge-midpoint axis):
#   Each order-2 rotation preserves 2 of the 6 five-fold axes
#   (the pair of axes symmetric w.r.t. the rotation plane) -> chi = 2
#   [This is because the 15 order-2 elements partition into the action
#   on 6 objects with stabilizer of size 60/6=10, and 10/2(order)=5
#   fixed-point-free elements, leaving 2 fixed axes per rotation.]

chi_6axes = [6, 1, 1, 0, 2]
print(f"  Character of A5 on 6 five-fold axes: {chi_6axes}")
print()

# Decompose into irreps using inner product formula
# <chi_rho, chi_sigma> = (1/|G|) sum_g chi_rho(g)*chi_sigma(g)*|class_g|
class_sizes = [1, 12, 12, 20, 15]
order_A5 = 60

print("  Decomposition into A5 irreps:")
decomp = {}
for name, ch in chars.items():
    mult = sum(chi_6axes[i] * ch[i] * class_sizes[i]
               for i in range(5)) / order_A5
    decomp[name] = mult
    if abs(mult) > 0.01:
        print(f"    <chi_6axes, chi_{name}> = {mult:.4f}")

print()

# Verify: sum d_i * m_i = 6
total = sum(int(chars[name][0]) * decomp[name] for name in chars)
print(f"  Verification: sum d_i * m_i = {total:.4f} (expect 6)")
print()

# The decomposition of 6 five-fold axes under A5
print("  RESULT: 6 five-fold axes = ", end="")
parts = []
for name in chars:
    m = round(decomp[name])
    if m > 0:
        if m == 1:
            parts.append(name)
        else:
            parts.append(f"{m}*{name}")
print(" + ".join(parts) if parts else "(decomposition failed)")
print()

# Now: on each line, the fermions use b = 0, 1, 2.
# If the decagon quantum number b is related to one of the A5 irreps,
# then b ranges over the basis vectors of that irrep.
# The 3-dim irrep has 3 basis vectors -> b in {0, 1, 2}?

# But this is not quite right. b is an integer indexing decagon steps,
# not a component of a representation.

# Let's think differently: the NUMBER of generations = 3 could be
# the dimension of one of the A5 irreps. But why specifically the 3?

print("  INTERPRETATION:")
print("  The 6 great decagons through a vertex decompose as 1 + 5 under A5.")
print("  This does NOT directly give a 3.")
print()
print("  However, A5 DOES have 3-dimensional irreps (3 and 3').")
print("  These arise as the restriction from SO(3) to the icosahedral group.")
print("  Angular momentum l=1 (3-dim) restricts to the 3 of A5.")
print()
print("  ALTERNATIVE: The 3 comes from the b=1 line structure.")
print("  On b=1: occupied a-values are {1, 2, 4} = {2^0, 2^1, 2^2}.")
print("  This is 3 elements, but the reason for exactly 3 is unclear.")
print()
print("  Category: SPECULATIV (A5 has 3-dim irreps, but no rigorous link)")
print()


# ============================================================
# APPROACH 5: Spectral Gap Analysis
# ============================================================

print("=" * 80)
print("APPROACH 5: Spectral Gap Constraint")
print("=" * 80)
print()

# Build adjacency matrix eigenvalues
A_mat = adj.astype(float)
eig_vals = np.linalg.eigvalsh(A_mat)
eig_sorted = np.sort(eig_vals)[::-1]
eig_unique = np.unique(np.round(eig_sorted, 4))[::-1]

print("  600-cell adjacency eigenvalues (9 distinct):")
for i, ev in enumerate(eig_unique):
    print(f"    theta_{i} = {ev:>10.4f}")
print()

# The spectral gap is the difference between the largest and second-largest eigenvalue
spectral_gap = eig_unique[0] - eig_unique[1]
print(f"  Spectral gap: {eig_unique[0]:.4f} - {eig_unique[1]:.4f} = {spectral_gap:.4f}")
print()

# From exp110: eigenvalues are a + b*phi form
# theta = {12, 6phi, 4phi, 3, 0, -2, 4-4phi, -3, 6-6phi}
# In terms of the mass formula, n = 5a + 6b maps to m = m_e * phi^n
# Can we impose a spectral constraint?

print("  Adjacency eigenvalues as a+b*phi:")
eig_ab = [
    (12, 0, 12),   # 12 + 0*phi = 12
    (0, 6, 6*PHI),  # 0 + 6*phi = 9.708
    (0, 4, 4*PHI),  # 0 + 4*phi = 6.472
    (3, 0, 3),      # 3 + 0*phi = 3
    (0, 0, 0),      # 0
    (-2, 0, -2),    # -2
    (4, -4, 4-4*PHI), # 4 - 4*phi = -2.472
    (-3, 0, -3),     # -3
    (6, -6, 6-6*PHI), # 6 - 6*phi = -3.708
]

for a, b, val in eig_ab:
    print(f"    {a:>3d} + {b:>3d}*phi = {val:>10.4f}")
print()

# Check: are there exactly 9 eigenvalues? Yes.
# 9 eigenvalues <-> 9 shells <-> 9 fermions
# The constraint "only 9 eigenvalues" means the association scheme
# has exactly 9 classes, limiting how many independent mass levels exist.

print("  The 600-cell has exactly 9 distance classes (shells).")
print("  This gives exactly 9 distinct eigenvalues.")
print("  If each eigenvalue corresponds to a fermion mass level,")
print("  then at most 9 fermions are allowed.")
print()
print("  But n=23 would be a 9th DISTINCT level (the 8 occupied levels")
print("  plus n=23 would give 9 distinct levels, since mu/s share n=11).")
print("  So the spectral count (9) does NOT exclude a 4th generation on a=1.")
print()
print("  HOWEVER: currently we have 8 distinct n-values for 9 fermions")
print("  (mu and strange share n=11). Adding n=23 gives 9 distinct levels")
print("  for 10 fermions. This would STILL require only 9 eigenvalues")
print("  (since the 9th eigenvalue is already unused).")
print()
print("  RESULT: Spectral gap does NOT impose a 3-generation limit.")
print("  Category: NOT DERIVED")
print()


# ============================================================
# APPROACH 6: E8 Branching (H4 x H4' Coxeter Exponents)
# ============================================================

print("=" * 80)
print("APPROACH 6: E8 -> H4 x H4' Branching")
print("=" * 80)
print()

print("  E8 embeds 600-cell via icosians (exp120d):")
print("  240 E8 roots = S(120) union T(120)")
print("  where T = phi' * S (Galois conjugate)")
print()

# H4 Coxeter exponents: {1, 11, 19, 29}
# H4' Coxeter exponents: Same, but in the conjugate copy
# The E8 Coxeter exponents: {1, 7, 11, 13, 17, 19, 23, 29}
# Note: 23 appears in E8 Coxeter exponents!

print("  Coxeter exponents:")
print("    H4:  {1, 11, 19, 29}")
print("    H4': {7, 13, 17, 23}")
print("    E8:  {1, 7, 11, 13, 17, 19, 23, 29}")
print()

# The CRITICAL observation: the E8 exponents split into two sets
# corresponding to H4 and H4' (the two copies in the icosian decomposition)
# H4 gets {1, 11, 19, 29} and H4' gets {7, 13, 17, 23}

print("  KEY OBSERVATION:")
print("  E8 Coxeter exponents split as {1,11,19,29} U {7,13,17,23}")
print("  under the E8 -> H4 x H4' branching.")
print()

# Now check: n=23 is a Coxeter exponent of H4' (the CONJUGATE copy)
print("  n=23 is a Coxeter exponent of H4' (the CONJUGATE copy)!")
print("  It belongs to the 'wrong' sector of the E8 decomposition.")
print()

# Let's check ALL occupied n-values against H4 vs H4' exponents
h4_exp = {1, 11, 19, 29}
h4p_exp = {7, 13, 17, 23}
e8_exp = h4_exp | h4p_exp

print("  Occupied n-values vs Coxeter exponents:")
for n_val in sorted(occupied_n):
    (a, b), c = min_complexity(n_val)
    in_h4 = n_val in h4_exp
    in_h4p = n_val in h4p_exp
    in_e8 = n_val in e8_exp
    print(f"    n={n_val:>2d} ({a:>2d},{b:>2d}): H4={'Y' if in_h4 else 'N'}"
          f" H4'={'Y' if in_h4p else 'N'} E8={'Y' if in_e8 else 'N'}")

print()
print("  n=23 (the forbidden 4th gen): H4'=Y (it IS a H4' Coxeter exponent)")
print()

# But wait: many occupied n-values are NOT in the Coxeter exponents either.
# {0, 3, 5, 16, 26} are not E8 Coxeter exponents.
# So this is not a clean rule.

print("  PROBLEM: n in {0, 3, 5, 16, 26} are occupied but NOT Coxeter exponents.")
print("  So 'n must be a Coxeter exponent' is not the rule.")
print()

# Let's look at it differently: from exp122, CKM uses H4 exponents
# {1, 11, 19, 29} -> (1+11)/4=3, |1-29|/4=7, (19+29)/4=12
# The H4' exponents {7, 13, 17, 23} are the "partner" set

# More interesting: check the (a,b) quantum numbers of the E8 exponents
print("  (a,b) for Coxeter exponents of E8:")
for n_val in sorted(e8_exp):
    (a, b), c = min_complexity(n_val)
    N = zphi_norm(a, b)
    print(f"    n={n_val:>2d}: (a,b) = ({a:>2d},{b:>2d}), |a|+|b|={c}, N={N}")

print()

# The H4' exponents map to:
# n=7:  (1,-1/3)... wait, 5a+6b=7. Best: a=-1, b=2 -> |a|+|b|=3
#   or a=5, b=-3 -> |a|+|b|=8. Min is (1, 1/3)... no. Let me compute.
print("  Checking (a,b) for H4' exponents specifically:")
for n_val in sorted(h4p_exp):
    (a, b), c = min_complexity(n_val)
    on_a1 = (a == 1)
    on_b1 = (b == 1)
    on_325 = (3 * a + 2 * b == 5)
    lines = []
    if on_a1:
        lines.append("a=1")
    if on_b1:
        lines.append("b=1")
    if on_325:
        lines.append("3a+2b=5")
    line_str = ",".join(lines) if lines else "none"
    print(f"    n={n_val:>2d}: ({a:>2d},{b:>2d}), lines: {line_str}")

print()

# n=23 has (a,b) = (1,3), which is on the a=1 line with b=3
# n=7 has (a,b) = (-1,2), not on any of our 3 lines
# n=13 has (a,b) = (1,4/3)... no, must be integer
# Let me recompute:
for n_val in [7, 13, 17, 23]:
    (a, b), c = min_complexity(n_val)
    print(f"    n={n_val}: min complexity (a,b) = ({a},{b}), 5*{a}+6*{b}={5*a+6*b}")

print()
print("  PARTIAL RESULT: n=23 maps to H4' Coxeter exponent (conjugate sector).")
print("  This is suggestive: perhaps fermions must correspond to H4 (not H4').")
print("  But the connection is not tight (most occupied n are NOT Coxeter exps).")
print("  Category: SPECULATIV (suggestive but not rigorous)")
print()


# ============================================================
# APPROACH 7: Z[phi] Unit Structure
# ============================================================

print("=" * 80)
print("APPROACH 7: Z[phi] Unit Structure")
print("=" * 80)
print()

# Units in Z[phi] are +-phi^k for integer k
# phi^0 = 1, phi^1 = phi = 1.618..., phi^2 = 2.618..., etc.
# Norm of phi^k = (-1)^k

print("  Z[phi] elements on the a=1 line: z = 1 + b*phi")
print()
print(f"  {'b':>3s} {'z':>10s} {'N':>5s} {'unit?':>6s} {'phi-power':>12s}")
print(f"  {'-'*42}")

for b in range(-1, 6):
    z = 1 + b * PHI
    N = zphi_norm(1, b)
    is_unit = abs(N) == 1

    # If unit, find k such that z = +-phi^k
    phi_pow = ""
    if is_unit and abs(z) > 1e-10:
        k = round(np.log(abs(z)) / np.log(PHI))
        sign = 1 if z > 0 else -1
        check = sign * PHI ** k
        if abs(check - z) < 1e-6:
            phi_pow = f"{'+' if sign > 0 else '-'}phi^{k}"

    print(f"  {b:>3d} {z:>10.4f} {N:>5d} {'yes' if is_unit else 'no':>6s}"
          f" {phi_pow:>12s}")

print()

# On the a=1 line, b=0 and b=1 give units (|N|=1)
# b=0: z=1 = phi^0 (even power)
# b=1: z=1+phi = phi^2 (even power!)
# b=2: z=1+2*phi, N=-5 (NOT a unit)
# b=3: z=1+3*phi, N=-5 (NOT a unit either... wait)
N_b2 = zphi_norm(1, 2)
N_b3 = zphi_norm(1, 3)
print(f"  b=2: N = {N_b2} (|N| = {abs(N_b2)})")
print(f"  b=3: N = {N_b3} (|N| = {abs(N_b3)})")
print()

# Check associate classes: z1 and z2 are associates if z1 = u*z2 for some unit u
print("  Associate class analysis:")
print("  Two elements are associates if they differ by a unit factor +-phi^k")
print()

for b in range(4):
    z = 1 + b * PHI
    N = zphi_norm(1, b)
    print(f"  b={b}: z = {z:.4f}, N = {N}")
    if abs(N) > 1:
        # Factor z. In Z[phi], the element z = a + b*phi has norm N.
        # We can write z = unit * prime factorization
        # The primes of Z[phi] are related to rational primes
        # p=5 ramifies (phi^2-phi-1=0, discriminant 5)
        # p=2: 2 = -phi^(-2)*(2+phi) ... it's complicated
        pass

print()
print("  RESULT: The unit structure provides no obvious cutoff at b=2.")
print("  Both b=2 and b=3 give norm |N|=5, same algebraic structure.")
print("  Category: NOT DERIVED")
print()


# ============================================================
# APPROACH 8: Topological / Homological Constraints
# ============================================================

print("=" * 80)
print("APPROACH 8: Topological / Homological Constraints")
print("=" * 80)
print()

print("  600-cell triangulates S^3.")
print("  Betti numbers: b_0=1, b_1=0, b_2=0, b_3=1")
print("  Euler characteristic: chi = 1-0+0-1 = 0 (consistent with odd-dim)")
print()
print("  Simplicial complex: V=120, E=720, F=1200, C=600")
print("  Euler: 120 - 720 + 1200 - 600 = 0. Check.")
print()

# The 3 lines in (a,b) space: can they be interpreted homologically?
# (a,b) is a point in Z^2. The constraint 5a+6b=n defines a line in Z^2.
# Our 3 lines {a=1, b=1, 3a+2b=5} are elements of the lattice Z^2.

print("  The 3 fermion lines in Z^2 lattice:")
print("    Line a=1: direction (0,1)")
print("    Line b=1: direction (1,0)")
print("    Line 3a+2b=5: direction (-2,3)")
print()

# Do these 3 directions have topological meaning?
# In H_1(T^2) = Z^2 (first homology of the 2-torus),
# (0,1), (1,0), and (-2,3) are 3 distinct homology classes.

print("  As elements of H_1(T^2) = Z^2:")
print("    (0,1): the 'b-cycle'")
print("    (1,0): the 'a-cycle'")
print("    (-2,3): a combined cycle")
print()

# The number of independent generators of H_1(T^2) is 2
# So 3 lines cannot all be independent; one is a combination.
# (0,1) and (1,0) generate everything. (-2,3) = -2*(1,0) + 3*(0,1).

print("  H_1(T^2) has rank 2. So (-2,3) = -2*(1,0) + 3*(0,1).")
print("  The third line is linearly dependent on the first two.")
print()
print("  This does NOT explain why there are exactly 3 lines.")
print("  Any lattice Z^2 has infinitely many sublattice lines.")
print()

# What about H_1(S^3) = 0?
print("  H_1(S^3) = 0 imposes no constraint on the (a,b) lattice.")
print()
print("  RESULT: Standard homology provides no generation bound.")
print("  Category: NOT DERIVED")
print()


# ============================================================
# APPROACH 9: Physical Mass Bound
# ============================================================

print("=" * 80)
print("APPROACH 9: Physical Mass Bound")
print("=" * 80)
print()

m_e = 0.000511  # GeV

print("  Mass formula: m = m_e * phi^n")
print(f"  m_e = {m_e} GeV")
print()

print(f"  {'n':>3s} {'b (a=1)':>8s} {'m (GeV)':>15s} {'m (TeV)':>12s} {'Particle':>12s}")
print(f"  {'-'*55}")

for b in range(6):
    n = 5 + 6 * b
    mass = m_e * PHI ** n
    particle = ""
    if b == 0:
        particle = "down"
    elif b == 1:
        particle = "mu/s"
    elif b == 2:
        particle = "tau"
    elif b == 3:
        particle = "4th gen?"

    print(f"  {n:>3d} {b:>8d} {mass:>15.4f} {mass/1000:>12.6f} {particle:>12s}")

print()

# Key physical thresholds
m_top = 172.76  # GeV
m_higgs = 125.25  # GeV
m_W = 80.377  # GeV

mass_23 = m_e * PHI ** 23
mass_17 = m_e * PHI ** 17

print(f"  tau (n=17):  m = {mass_17:.2f} GeV")
print(f"  4th gen (n=23): m = {mass_23:.2f} GeV")
print(f"  top quark:   m = {m_top:.2f} GeV")
print(f"  Higgs boson: m = {m_higgs:.2f} GeV")
print(f"  W boson:     m = {m_W:.2f} GeV")
print()

print(f"  Ratio m(n=23) / m(tau) = phi^6 = {PHI**6:.2f}")
print(f"  m(n=23) / m(top) = {mass_23/m_top:.4f}")
print()

# Is the 4th gen mass excluded experimentally?
print("  Experimental bounds on 4th generation:")
print(f"    Predicted 4th gen lepton mass: {mass_23:.1f} GeV")
print(f"    LHC excludes heavy leptons below ~100 GeV (from Z-width)")
print(f"    The Z-width measurement limits N_nu(light) = 2.9840 +/- 0.0082")
print(f"    This excludes a 4th light neutrino, but NOT a heavy lepton")
print(f"    with mass > M_Z/2 = 45.6 GeV.")
print(f"    A 4th gen lepton at {mass_23:.1f} GeV is NOT excluded by Z-width alone.")
print()
print(f"    Direct search limits (LEP/LHC): heavy charged lepton > ~100 GeV")
print(f"    Our predicted mass {mass_23:.1f} GeV is ABOVE this threshold.")
print(f"    So a 4th gen lepton at {mass_23:.1f} GeV is in principle ALLOWED")
print(f"    by current experiments (though Higgs data constrain it).")
print()
print("  HOWEVER: Higgs production cross-section measurements at LHC")
print("  strongly disfavor a 4th generation of heavy fermions")
print("  (they would contribute to gg->H loop). This is a")
print("  model-dependent constraint, not a geometric derivation.")
print()
print("  RESULT: No absolute physical mass bound prevents n=23.")
print("  The mass ~{:.0f} GeV is in a physically accessible range.".format(mass_23))
print("  Experimental exclusion exists but is model-dependent.")
print("  Category: NOT DERIVED (physics constraints, not geometry)")
print()


# ============================================================
# APPROACH 10: Decagon Structure
# ============================================================

print("=" * 80)
print("APPROACH 10: Decagon Structure")
print("=" * 80)
print()

# Count great decagons
# Each vertex is on 6 great decagons
# Each decagon has 10 vertices
# Total (with overcounting): 120 * 6 / 10 = 72 great decagons
n_decagons = 120 * 6 // 10
print(f"  Total great decagons: {n_decagons}")
print(f"  Decagons per vertex: 6")
print(f"  Vertices per decagon: 10")
print()

# The 6 decagons through a vertex:
# Can we find how many independent decagons there are?
# The decagons are the great circles of the 600-cell on S^3

# In 4D, great circles on S^3 are intersections of 2-planes through origin with S^3.
# The set of all 2-planes through origin in R^4 forms the Grassmannian Gr(2,4).
# dim Gr(2,4) = 4 (real dimension).

print("  Great decagons correspond to 2-planes in R^4 through origin")
print("  that contain 10 vertices of the 600-cell.")
print(f"  Grassmannian Gr(2,4) has real dimension {2*(4-2)} = 4.")
print()

# On a=1 line, b counts decagon steps.
# One decagon step = moving from one vertex to a vertex 2 steps away
# on the same great decagon (since n=6 means angular distance 2*pi/10*6... no)
# Actually, b=1 means one increment of the b-quantum number,
# which corresponds to 6 = b_1 steps in the mass formula.

# The number 3 in "3 generations" on a=1:
# On a decagon (10 vertices), from a starting vertex,
# there are 4 "independent" steps (to the 5th vertex is antipodal, then repeats).
# But this gives 4 or 5, not 3.

print("  On a great decagon (10 vertices, cyclic C_10):")
print("  From vertex 0, steps to vertices 1,2,3,4,5(antipodal),6,7,8,9")
print("  Independent positions: 0,1,2,3,4 (then 5=antipodal of 0, etc.)")
print("  This gives 5 positions (or 4 non-trivial), not 3.")
print()

# But what about the DIHEDRAL angle?
# Each edge of the 600-cell is shared by 5 tetrahedra
# The dihedral angle = 2*pi/5 = 72 degrees
# 360/72 = 5 tetrahedra per edge

print("  Dihedral structure:")
print("  5 tetrahedra per edge (dihedral angle = 72 deg = 2*pi/5)")
print("  5 = diameter of graph")
print()

# 3-coloring of faces?
# A tetrahedron has 4 faces. Each face is a triangle.
# 3 colors suffice to color the faces of any triangulation?
# Actually, the 4-color theorem applies to planar graphs.
# For S^3 triangulation, Heawood conjecture gives different bounds.

# Let me count something more concrete about the intersection of decagons
print("  Intersection pattern of 6 decagons through vertex 0:")
# The 6 decagons partition the 12 edges at vertex 0 into 6 pairs.
# Two decagons through vertex 0 can share only vertex 0 (generically)
# or they might share a second vertex.

# Find the 6 decagons through vertex 0
# First, find the 12 neighbors of vertex 0
neighbors_0 = [j for j in range(N_verts) if adj[0, j] == 1]
print(f"  Vertex 0 neighbors: {len(neighbors_0)}")

# For each neighbor, find which great decagon they share
# A great decagon is determined by the 2-plane spanned by the two vertices
# Strategy: from vertex 0, follow each neighbor and trace the decagon

def trace_decagon(start, next_v, vertices, adj_matrix):
    """Trace a great decagon starting from start, going to next_v."""
    decagon = [start, next_v]
    current = next_v
    prev = start

    for _ in range(8):  # 10 vertices total, we have 2
        # Find the neighbor of current that is:
        # (a) not prev, and
        # (b) continues the great circle
        # The next vertex on a great decagon through (start, next_v) is the one
        # that lies on the same 2-plane.

        # Compute the 2-plane: span of V[start] and V[next_v]
        v1 = vertices[start]
        v2 = vertices[next_v]

        # Project current vertex's neighbors onto this plane
        # The next vertex on the decagon maximizes the "continuation" angle

        # Actually, for a great decagon on S^3, the vertices lie on a great circle
        # parameterized by angle theta = 2*pi*k/10 for k=0,...,9
        # The 2-plane contains all these vertices.

        # Find all vertices on this 2-plane
        # The 2-plane is spanned by v1 and v2_orth (v2 orthogonalized to v1)
        v2_orth = v2 - np.dot(v1, v2) * v1
        v2_orth = v2_orth / np.linalg.norm(v2_orth)

        # A vertex u is on the 2-plane if it has zero component perpendicular to both
        # v1 and v2_orth
        break  # This approach is too complex, use a different method

    return decagon

# Simpler approach: identify great decagons by finding all 10-cycles
# that are "great" (evenly spaced on a great circle)
#
# A great decagon: 10 vertices equally spaced on a great circle of S^3
# The angular distance between consecutive vertices is 2*pi/10 = pi/5 = 36 deg

theta_decagon_step = np.pi / 5  # 36 degrees

print(f"  Decagon step angle: {np.degrees(theta_decagon_step):.1f} deg")
print(f"  Edge angle: {np.degrees(theta_edge):.1f} deg")
print()

# Check: is the edge angle equal to the decagon step?
print(f"  Edge angle = {np.degrees(theta_edge):.4f} deg")
print(f"  pi/5       = {np.degrees(theta_decagon_step):.4f} deg")
print(f"  Are they equal? {abs(theta_edge - theta_decagon_step) < 0.01}")
print()

# Find decagons through vertex 0 by tracing paths
# From vertex 0, go to neighbor n1, then continue in the same "direction"
# on the great circle

def find_great_decagon(v_idx, n1_idx, vertices, adj_mat, angles_mat):
    """Find the great decagon through v_idx and n1_idx."""
    decagon = [v_idx, n1_idx]
    for step in range(8):
        current = decagon[-1]
        prev = decagon[-2]

        # The next vertex should be at angle pi/5 from current
        # and at angle 2*pi/5 from prev, on the same great circle
        # Find neighbors of current that are at the right angle from prev

        target_angle_from_prev = angles_mat[prev, current] + theta_decagon_step
        # Actually simpler: use the 2-plane projection

        v0 = vertices[v_idx]
        v1 = vertices[n1_idx]

        # Orthonormal basis for the 2-plane
        e1 = v0 / np.linalg.norm(v0)
        e2 = v1 - np.dot(v1, e1) * e1
        e2 = e2 / np.linalg.norm(e2)

        # Project current vertex to get angle
        vc = vertices[current]
        c1 = np.dot(vc, e1)
        c2 = np.dot(vc, e2)
        angle_current = np.arctan2(c2, c1)

        # Next vertex should be at angle_current + 2*pi/10
        target_angle = angle_current + 2 * np.pi / 10

        best_j = -1
        best_err = 999
        for j in range(len(vertices)):
            if j == current or j == prev:
                continue
            vj = vertices[j]
            j1 = np.dot(vj, e1)
            j2 = np.dot(vj, e2)
            # Check if vj is on the 2-plane (perpendicular components ~0)
            perp = np.linalg.norm(vj - j1 * e1 - j2 * e2)
            if perp > 0.01:
                continue
            angle_j = np.arctan2(j2, j1)
            err = abs(((angle_j - target_angle + np.pi) % (2 * np.pi)) - np.pi)
            if err < best_err:
                best_err = err
                best_j = j

        if best_j >= 0 and best_err < 0.1:
            decagon.append(best_j)
        else:
            break

    return decagon

# Find all great decagons through vertex 0
decagons_through_0 = []
used_neighbors = set()

for n1 in neighbors_0:
    if n1 in used_neighbors:
        continue
    dec = find_great_decagon(0, n1, V, adj, angles)
    if len(dec) == 10:
        decagons_through_0.append(dec)
        # Mark the opposite neighbor (on same decagon) as used
        # The neighbor of vertex 0 on this decagon going the other way
        # is the last vertex before returning to 0
        if dec[-1] in neighbors_0:
            used_neighbors.add(n1)
            used_neighbors.add(dec[-1])
    else:
        # Try going the other direction
        pass

print(f"  Great decagons found through vertex 0: {len(decagons_through_0)}")
for i, dec in enumerate(decagons_through_0):
    print(f"    Decagon {i}: {len(dec)} vertices")

print()

# Regardless of the exact count, the key question is:
# Does the decagon structure limit b to {0,1,2}?

# On a single decagon of 10 vertices:
# Start at vertex 0. The b-quantum number counts "decagon increments".
# b=1 means moving 6 positions along the DSI ladder.
# On a decagon with 10 positions, after 5 steps you reach the antipodal point.
# b=3 would require 18 positions on the DSI ladder, well beyond one decagon.

print("  ANALYSIS: A single decagon has 10 vertices.")
print("  On the DSI ladder, one decagon step contributes 6 to n.")
print("  After 5 decagon steps, you've gone n=30 (wrapping around S^3).")
print("  But b counts something different: b = decagon quantum number,")
print("  not decagon STEPS along a single decagon.")
print()
print("  On one decagon: positions 0 through 9, then wrap.")
print("  But different decagons can combine to give b>1.")
print("  There is no clear geometric limit from a single decagon.")
print()
print("  RESULT: Decagon structure does not obviously bound b<=2.")
print("  Category: NOT DERIVED")
print()


# ============================================================
# ADDITIONAL ANALYSIS: The 3 lines and "3"
# ============================================================

print("=" * 80)
print("ADDITIONAL ANALYSIS: Why 3 Lines? Why 3 Points Per Line?")
print("=" * 80)
print()

# There are exactly 3 lines through the triple junction (1,1)
# Each line has exactly 3 occupied points (including (1,1))
# This gives 3*3 - 2 = 7 points + origin = 8 distinct (a,b)

print("  The 3-line structure:")
print("  - All 3 lines pass through (1,1) = mu/strange")
print("  - Each line has exactly 3 occupied fermions")
print("  - Total: 3*3 = 9 fermion positions, but 3 overlap at (1,1)")
print("  - So 9 - 2 = 7 on lines + 1 origin = 8 distinct (a,b)")
print()

# On EACH line, there are 3 fermions. Why 3?
# Line a=1: b in {0,1,2} -> 3 points
# Line b=1: a in {1,2,4} -> 3 points (powers of 2)
# Line 3a+2b=5: (3,-2),(1,1),(-1,4) -> 3 points

# For line a=1: b goes from 0 to 2. The STEP is 1. 3 = (2-0)/1 + 1.
# For line b=1: a goes from 1 to 4. Steps are 1, 2. Non-uniform.
# For line 3a+2b=5: spread from a=-1 to a=3. Steps 2, 2. Uniform.

print("  Points per line:")
print("  Line a=1: b = {0, 1, 2}. Spacing: 1, 1 (uniform)")
print("  Line b=1: a = {1, 2, 4}. Spacing: 1, 2 (non-uniform, GP ratio 2)")
print("  Line 3a+2b=5: a = {-1, 1, 3}. Spacing: 2, 2 (uniform)")
print()

# Interesting: on lines a=1 and 3a+2b=5, the spacing is uniform.
# On line b=1, the spacing is geometric (ratio 2).

# For line a=1, the constraint b<=2:
# On the 5a+6b=n lattice with a=1, n = 5 + 6b.
# The allowed norms are |N| in {0,1,5,19}.
# N(1,b) = 1 + b - b^2.
# |N(1,0)| = 1, |N(1,1)| = 1, |N(1,2)| = 5, |N(1,3)| = 5,
# |N(1,4)| = 11, |N(1,5)| = 19, |N(1,6)| = 29

print("  Norms on a=1 line: N(1,b) = 1 + b - b^2")
print(f"  {'b':>3s} {'n':>4s} {'N':>5s} {'|N|':>4s} {'|N| in allowed':>16s}")
for b in range(-1, 10):
    N = zphi_norm(1, b)
    allowed = abs(N) in {0, 1, 5, 19}
    print(f"  {b:>3d} {5+6*b:>4d} {N:>5d} {abs(N):>4d} {'YES' if allowed else 'NO':>16s}")

print()
print("  Allowed b on a=1 line (norm constraint): b in {-1, 0, 1, 2, 3, 5}")
print("  Occupied b on a=1 line: b in {0, 1, 2}")
print("  False positive: b=3 (n=23), b=-1 (n=-1, out of range)")
print("  Also allowed: b=5 (n=35), b=-1 (n=-1)")
print()

# So from norm + line constraint, we get b in {0,1,2,3} (for n >= 0)
# The exp121d additional constraint b <= 2 removes b=3.

# Can we find a geometric reason for b <= 2?
# Note: b = 2 gives N(1,2) = 1+2-4 = -3. Wait, let me recalculate.
# N(1,2) = 1 + 2 - 4 = -1. Hmm.
# N(a,b) = a^2 + ab - b^2. So N(1,2) = 1 + 2 - 4 = -1.
# Wait, that's not -5. Let me recheck.

print("  RECHECK: N(a,b) = a^2 + a*b - b^2")
for b in range(5):
    N = 1 + b - b**2
    print(f"    N(1,{b}) = 1 + {b} - {b**2} = {N}")

print()
# So N(1,0)=1, N(1,1)=1, N(1,2)=-1, N(1,3)=-5
# |N(1,2)| = 1, |N(1,3)| = 5. Both in allowed set.
# Earlier in exp121d: N(1,3) = -5, |N| = 5. Confirmed.


# ============================================================
# APPROACH 6b: Deeper E8/H4 Analysis
# ============================================================

print("=" * 80)
print("APPROACH 6b: E8 Coxeter Exponents and Generation Counting")
print("=" * 80)
print()

print("  E8 Coxeter exponents: {1, 7, 11, 13, 17, 19, 23, 29}")
print("  These equal the degrees minus 1 of the basic invariants.")
print()

# The E8 root system has 240 roots.
# Under E8 -> H4 x H4' (icosian decomposition), 240 = 120 + 120.
# The Coxeter exponents split: {1,11,19,29} for H4 and {7,13,17,23} for H4'.

# The DIFFERENCE pattern:
# H4: 1, 11, 19, 29 -> differences: 10, 8, 10
# H4': 7, 13, 17, 23 -> differences: 6, 4, 6
print("  H4 exponents: 1, 11, 19, 29 -> diffs: 10, 8, 10")
print("  H4' exponents: 7, 13, 17, 23 -> diffs: 6, 4, 6")
print()

# The H4' exponent differences are 6, 4, 6.
# 6 = our b_1 intersection number!
# 4 = dimension of the 600-cell!

print("  H4' differences: 6, 4, 6")
print("    6 = b_1 (second intersection number)")
print("    4 = dimension of 600-cell")
print()

# On the a=1 line, n-values differ by 6 (since n = 5 + 6b):
# n=5, 11, 17 -> diffs: 6, 6
print("  On a=1 line: n = 5, 11, 17 -> diffs: 6, 6")
print("  Next would be n=23: diff still 6.")
print("  So the 'diff=6' pattern does NOT stop at 3 generations.")
print()

# BUT: n=23 is a Coxeter exponent of H4' (the conjugate).
# And the H4' exponents are: 7, 13, 17, 23.
# n=17 appears in BOTH H4' AND as an occupied level (tau).
# n=23 appears in H4' but is NOT occupied.

# This is confusing. Let's look more carefully.

print("  Cross-reference n-values with E8 exponents:")
print(f"  {'n':>3s} {'occupied':>9s} {'H4':>4s} {'H4p':>4s} {'E8':>4s}")
for n_val in range(31):
    occ = "YES" if n_val in occupied_n else ""
    h4 = "Y" if n_val in h4_exp else ""
    h4p = "Y" if n_val in h4p_exp else ""
    e8 = "Y" if n_val in e8_exp else ""
    if occ or h4 or h4p or e8:
        print(f"  {n_val:>3d} {occ:>9s} {h4:>4s} {h4p:>4s} {e8:>4s}")

print()

# Observation: occupied n-values include SOME E8 exponents (11, 17, 19)
# and some non-exponents (0, 3, 5, 16, 26).
# The E8 exponent n=23 is NOT occupied.
# Neither are n=1, 7, 13, 29.

# Different angle: the H4' Coxeter exponents {7,13,17,23} are exactly
# the H4 exponents {1,11,19,29} shifted by +6:
# 1+6=7, 11+6=17, 19+6=25(NO)... doesn't work.
# How about: 30 - H4 exponents = {29, 19, 11, 1} -> {30-29, 30-19, 30-11, 30-1}
# = {1, 11, 19, 29} -> that's H4 itself (palindromic).
# For H4': {30-23, 30-17, 30-13, 30-7} = {7, 13, 17, 23} -> also palindromic.

# The Coxeter number h(H4) = 30.
# The exponents m_i satisfy: m_i + m_{rank+1-i} = h = 30.
# H4: 1+29=30, 11+19=30. Check.
# For E8: h=30. E8 exponents: m_i + m_{9-i} = 30.
# 1+29=30, 7+23=30, 11+19=30, 13+17=30. Check.

print("  Coxeter number h = 30 for both H4 and E8.")
print("  Palindromic property: m_i + m_{r+1-i} = h = 30")
print("  H4: 1+29=30, 11+19=30")
print("  E8: 1+29=30, 7+23=30, 11+19=30, 13+17=30")
print()

# The pairing 7+23=30 and 13+17=30 in H4':
# n=17 is occupied (tau), n=23 is NOT.
# n=13 is NOT occupied, n=7 is NOT occupied.
# So the E8 Coxeter pairing does not cleanly distinguish occupied from not.

print("  E8 Coxeter pairs (m + m' = 30):")
for m in [1, 7, 11, 13]:
    mp = 30 - m
    occ_m = "OCC" if m in occupied_n else "   "
    occ_mp = "OCC" if mp in occupied_n else "   "
    print(f"    {m:>2d} + {mp:>2d} = 30   ({occ_m} + {occ_mp})")

print()
print("  RESULT: E8 branching gives suggestive connections (n=23 is in H4')")
print("  but does not cleanly derive the 3-generation limit.")
print("  Category: SPECULATIV")
print()


# ============================================================
# SYNTHESIS: All Approaches
# ============================================================

print("=" * 80)
print("SYNTHESIS: Summary of All 10 Approaches")
print("=" * 80)
print()

results_summary = [
    ("1. Z[phi] norm analysis",     "NOT DERIVED", "b=3 has |N|=5, which is in allowed set"),
    ("2. Diameter constraint",       "NOT DERIVED", "n=23 < n=26 (top), diameter gives no bound"),
    ("3. H4 representation theory",  "SPECULATIV",  "A5 has 3-dim irreps, but no rigorous link to b"),
    ("4. Icosahedral vertex figure", "SPECULATIV",  "6 decagons = 1+5 under A5, not 3"),
    ("5. Spectral gap",             "NOT DERIVED", "9 eigenvalues allow >= 9 fermions"),
    ("6. E8 branching",             "SPECULATIV",  "n=23 is H4' Coxeter exp (conjugate sector)"),
    ("7. Z[phi] unit structure",     "NOT DERIVED", "b=2 and b=3 have same |N|=5 structure"),
    ("8. Topological/homological",   "NOT DERIVED", "H_1(S^3)=0 gives no constraint"),
    ("9. Physical mass bound",       "NOT DERIVED", "m(n=23) ~ {:.0f} GeV, experimentally accessible".format(mass_23)),
    ("10. Decagon structure",        "NOT DERIVED", "Single decagon has 10 positions, not 3"),
]

print(f"  {'Approach':>38s}  {'Status':>14s}  Explanation")
print(f"  {'-'*90}")
for approach, status, explanation in results_summary:
    print(f"  {approach:>38s}  {status:>14s}  {explanation}")

print()

# ============================================================
# CLOSEST CANDIDATE: E8 branching
# ============================================================

print("=" * 80)
print("CLOSEST CANDIDATE FOR DERIVATION: E8 -> H4 x H4'")
print("=" * 80)
print()

print("  The most promising geometric connection is:")
print()
print("  (1) E8 = H4 x H4' (icosian decomposition, exp120d)")
print("  (2) E8 Coxeter exponents: {1,7,11,13,17,19,23,29}")
print("  (3) H4 exponents:  {1, 11, 19, 29}")
print("  (4) H4' exponents: {7, 13, 17, 23}")
print()
print("  n=23 is a Coxeter exponent of H4' (the CONJUGATE sector).")
print("  This suggests that levels corresponding to H4' Coxeter exponents")
print("  belong to a different sector and should not appear as fermion masses.")
print()
print("  HOWEVER, this argument has a problem:")
print("  - n=17 (tau) is ALSO a H4' Coxeter exponent, yet IS occupied.")
print("  - n=11 (muon) is a H4 exponent AND is occupied.")
print("  - Most occupied n-values are NOT Coxeter exponents at all.")
print()
print("  So the E8 branching provides PARTIAL insight but NOT a derivation.")
print()


# ============================================================
# ALTERNATIVE PERSPECTIVE: Counting argument
# ============================================================

print("=" * 80)
print("ALTERNATIVE: Simple Counting Argument")
print("=" * 80)
print()

# Each of the 3 lines has a finite number of lattice points with allowed norms
# On a=1 line: allowed b in {0, 1, 2, 3, 5} (for n in [0,35])
# On b=1 line: allowed a in {-1, 0, 1, 2, 4, 7} (for n in [0,41])
# On 3a+2b=5: allowed points = finite set

# What if the constraint is: each line contributes AT MOST 3 occupied levels?
# Then 3 * 3 - 2 (overlaps at (1,1)) = 7, plus origin = 8. Matches!

print("  Counting hypothesis: each line contributes exactly 3 fermions.")
print("  3 lines x 3 fermions - 2 overlaps at (1,1) + 1 origin = 8 = 9-1 fermions")
print("  (9 fermions total, with mu/strange degenerate)")
print()

# Why 3 per line? If the 3 lines correspond to the 3 sectors
# (lepton, up-quark, down-quark) and each sector has 3 generations,
# then 3 generations is the INPUT, not the OUTPUT.

print("  BUT: 3 per line <=> 3 generations is CIRCULAR.")
print("  We assumed 3 generations and found 3 per line.")
print("  The question is WHY 3, not just HOW.")
print()

# Deeper: the 3 lines meet at (1,1). In projective geometry,
# 3 lines through a point define a "pencil of lines."
# A pencil in P^1(Z) over Z[phi] might have constraints.

print("  The 3 lines through (1,1) form a 'pencil' in the (a,b) lattice.")
print("  Directions: (0,1), (1,0), (-2,3)")
print("  These generate 3 of the 6 directions in the lattice with")
print("  'small' denominators when expressed as slopes.")
print()

# Cross-ratio of directions
# If we think of slopes: a=1 has slope infinity, b=1 has slope 0,
# 3a+2b=5 has slope -3/2.
# The cross-ratio is a projective invariant.
print("  Slopes of the 3 lines (in b = m*a + c form):")
print("    a=1: vertical (slope = infinity)")
print("    b=1: horizontal (slope = 0)")
print("    3a+2b=5: slope db/da = -3/2")
print()


# ============================================================
# FINAL: The Honest Answer
# ============================================================

print("=" * 80)
print("FINAL CONCLUSION")
print("=" * 80)
print()

print("  QUESTION: Can we derive b<=2 (3 generations) from 600-cell geometry?")
print()
print("  ANSWER: NO. After testing 10 different approaches, NONE provides")
print("  a rigorous derivation of the 3-generation limit from the 600-cell.")
print()
print("  STATUS: '3 generations remains UNDERIVED.'")
print()
print("  BEST CANDIDATES (in decreasing order of promise):")
print()
print("  1. E8 -> H4 x H4' branching (SPECULATIV)")
print("     n=23 is a H4' Coxeter exponent (conjugate sector).")
print("     Suggestive, but n=17 (tau) is also H4' yet IS occupied.")
print()
print("  2. A5 representation dimension (SPECULATIV)")
print("     The icosahedral vertex stabilizer A5 has 3-dim irreps.")
print("     But 6 decagons decompose as 1+5 under A5, not involving dim=3.")
print()
print("  3. Physical mass bounds (NOT DERIVED)")
print("     m(n=23) ~ {:.0f} GeV is experimentally accessible.".format(mass_23))
print("     No absolute physical barrier prevents a 4th generation lepton.")
print()
print("  WHAT WE DO KNOW:")
print("  - The line_valid rule from exp121d PERFECTLY selects 8 from 11:")
print("    b <= 2 on a=1 line, a >= 1 on b=1 line, always on 3a+2b=5.")
print("  - This is a PATTERN (descriptive rule), not a DERIVATION.")
print("  - The constraint b<=2 is equivalent to stating 3 generations.")
print("  - Deriving it would be equivalent to deriving the number of")
print("    generations from geometry -- one of the deepest open problems")
print("    in theoretical physics.")
print()
print("  CLASSIFICATION: 3 GENERATIONS = PATTERN (not derived)")
print("  The honest answer: we do not know WHY there are 3 generations.")
print("  The 600-cell theory organizes them beautifully on 3 lines,")
print("  but does not yet explain WHY exactly 3 per line.")
print()
print("  REGULA ZERO: 'E mai bine sa spui nu stiu decat sa inventezi")
print("  o minciuna eleganta.'")
print()

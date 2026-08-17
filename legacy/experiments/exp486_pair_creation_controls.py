"""
EXP-486: Control Tests for Pair Creation Near-Resonance
========================================================
QUESTION 1: Is the Fibonacci near-resonance (delta_omega/omega = 0.38%)
a NON-TRIVIAL result, or is it a BANAL consequence of any irrational
number producing continued-fraction approximations?

TEST: Replace phi with sqrt(2), sqrt(3), other irrationals.
Build analogous "spectral" structures and check if similar
near-resonances appear with THEIR continued-fraction convergents.

QUESTION 2: The triangle-edge participation ratios 336:576:432
simplify to 7:12:9 (after dividing by GCD=48).
12 = degree, 9 = N_eig, 7 = n_23 (CKM exponent).
Coincidence or structure?

STATUS: CONTROL EXPERIMENT (testing for Texas Sharpshooter)
"""

import numpy as np
from fractions import Fraction
from collections import defaultdict

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6

print("=" * 70)
print("EXP-486: CONTROL TESTS FOR PAIR CREATION")
print("=" * 70)

# ============================================================
# PART 1: IS THE FIBONACCI NEAR-RESONANCE TRIVIAL?
# ============================================================
print("\n" + "=" * 70)
print("PART 1: FIBONACCI NEAR-RESONANCE CONTROL")
print("=" * 70)

print("""
Setup: In exp485, the pair creation threshold eigenvalue is
  lam_thresh = 4 * (12 - 6*phi) = 48 - 24*phi
and the nearest coexact eigenvalue above it is
  lam_res = b1 + 2*phi = 6 + 2*phi

The near-resonance gap is:
  delta = (6+2*phi) - (48-24*phi) = 26*phi - 42
  delta/lam_thresh = 0.0075 (0.75%)
  delta_omega/omega = 0.38%

CLAIM TO TEST: This near-resonance comes from phi ~ 21/13 (Fibonacci).
CONTROL: For generic irrationals x, do analogous expressions
  lam_thresh(x) = 48 - 24*x
  lam_res(x) = 6 + 2*x
give similarly small gaps when x ~ p/q (convergent)?
""")

# The gap function: delta(x) = (6+2x) - (48-24x) = 26x - 42
# This vanishes when x = 42/26 = 21/13
# For phi: 26*phi - 42 = 26*(1+sqrt(5))/2 - 42 = 13*(1+sqrt(5)) - 42
#         = 13 + 13*sqrt(5) - 42 = 13*sqrt(5) - 29

print("--- Analysis of the gap function ---")
print(f"delta(x) = 26*x - 42")
print(f"delta = 0 when x = 21/13 = {21/13:.10f}")
print(f"phi = {phi:.10f}")
print(f"phi - 21/13 = {phi - 21/13:.10f}")
print(f"|delta(phi)| = |26*phi - 42| = {abs(26*phi - 42):.10f}")
print()

# KEY QUESTION: Is it just that 21/13 happens to be a convergent of phi?
# The continued fraction convergents of phi are:
# 1/1, 2/1, 3/2, 5/3, 8/5, 13/8, 21/13, 34/21, 55/34, ...
# ALL are Fibonacci ratios F_{n+1}/F_n

# For ANY irrational x, the continued fraction convergents p_n/q_n
# satisfy |x - p_n/q_n| < 1/(q_n * q_{n+1})
# So |26*x - 42| is small when 42/26 = 21/13 is near a convergent of x.

# BUT: the coefficients 26 and 42 in our theory come from:
# 26 = a1^2 + 1 (denominator of sin^2(tW))
# 42 = ... let's check: 42 = 2*21 = 2*F_8. Or 42 = 6*7. Or 42 = b1 * 7.
# Actually: 48 = 4*12 = 4*degree, 24 = 2*12 = 2*degree
# And: 6 = b1, 2 = coefficient

# The REAL question is: do the specific coefficients (48,24) and (6,2)
# that emerge from the 600-cell geometry conspire to give 26x-42,
# which happens to vanish near phi? Or could ANY two linear expressions
# ax+b and cx+d give similar near-resonances for ANY irrational?

print("--- The structural question ---")
print(f"Threshold: 4*(12-6x) = 48-24x   (4 = (a1-1), 12 = degree, 6 = b1)")
print(f"Resonance: b1+2x = 6+2x         (b1 = a1+1, 2 = coefficient)")
print(f"Gap: (6+2x) - (48-24x) = 26x - 42")
print(f"  26 = a1^2+1  (denominator of sin^2 theta_W!)")
print(f"  42 = 26*21/13 ... but also 42 = b1*7 = (a1+1)*7")
print()

# CONTROL TEST 1: Replace phi with other irrationals
print("--- CONTROL TEST 1: Other irrationals ---")
print("For each irrational x, compute |26x - 42| / (6+2x)")
print("and compare to the best continued fraction approximation of x near 21/13")
print()

test_irrationals = {
    'phi = (1+sqrt(5))/2': phi,
    'sqrt(2)': np.sqrt(2),
    'sqrt(3)': np.sqrt(3),
    'sqrt(5)': np.sqrt(5),
    'e (Euler)': np.e,
    'pi': np.pi,
    'sqrt(7)': np.sqrt(7),
    'cbrt(2)': 2**(1/3),
    '(1+sqrt(2)) (silver)': 1 + np.sqrt(2),
    'plastic number': 1.3247179572,  # real root of x^3=x+1
}

print(f"{'Irrational':>25} {'x':>10} {'|26x-42|':>10} {'6+2x':>8} "
      f"{'gap_rel':>10} {'gap_omega':>10} {'Near p/q':>12}")
print("-" * 95)

for name, x in test_irrationals.items():
    gap = abs(26*x - 42)
    lam_res = 6 + 2*x
    lam_thresh = 48 - 24*x

    if lam_thresh > 0 and lam_res > 0:
        gap_rel = gap / lam_res
        omega_res = np.sqrt(lam_res)
        omega_thresh = np.sqrt(abs(lam_thresh))
        gap_omega = abs(omega_res - omega_thresh) / omega_thresh if lam_thresh > 0 else float('inf')
    else:
        gap_rel = float('inf')
        gap_omega = float('inf')

    # Find closest fraction with denominator <= 20
    best_frac = None
    best_err = float('inf')
    for q in range(1, 21):
        p = round(x * q)
        err = abs(x - p/q)
        if err < best_err:
            best_err = err
            best_frac = f"{p}/{q}"

    print(f"{name:>25} {x:10.6f} {gap:10.6f} {lam_res:8.4f} "
          f"{gap_rel:10.6f} {gap_omega:10.6f} {best_frac:>12}")

# CONTROL TEST 2: Is the near-resonance SPECIFIC to the 600-cell coefficients?
print("\n--- CONTROL TEST 2: Sensitivity to coefficients ---")
print("Vary the coefficients slightly and check if near-resonance persists")
print()

# Original: threshold = 48-24x, resonance = 6+2x
# These come from: 4*(12-6x) and b1+2x
# What if we use different multiples?

print(f"{'Threshold':>15} {'Resonance':>12} {'Gap coeff':>10} {'Zero at x=':>12} "
      f"{'|x_0 - phi|':>12} {'gap(phi)':>10} {'gap_rel':>10}")
print("-" * 95)

variants = [
    ("4*(12-6x)", "6+2x",      48, -24, 6, 2,   "ORIGINAL"),
    ("3*(12-6x)", "6+2x",      36, -18, 6, 2,   "factor 3"),
    ("5*(12-6x)", "6+2x",      60, -30, 6, 2,   "factor 5"),
    ("4*(12-6x)", "6+3x",      48, -24, 6, 3,   "coeff 3"),
    ("4*(12-6x)", "6+x",       48, -24, 6, 1,   "coeff 1"),
    ("4*(12-6x)", "7+2x",      48, -24, 7, 2,   "const 7"),
    ("2*(12-6x)", "6+2x",      24, -12, 6, 2,   "factor 2"),
    ("4*(12-5x)", "6+2x",      48, -20, 6, 2,   "5x not 6x"),
    ("4*(12-6x)", "5+2x",      48, -24, 5, 2,   "a1 not b1"),
]

for label, res_label, a, b, c, d, note in variants:
    # threshold = a + b*x, resonance = c + d*x
    # gap = (c+d*x) - (a+b*x) = (c-a) + (d-b)*x
    gap_coeff = d - b  # coefficient of x
    gap_const = c - a  # constant term
    if gap_coeff != 0:
        x_zero = -gap_const / gap_coeff
    else:
        x_zero = float('inf')

    gap_phi = abs(gap_const + gap_coeff * phi)
    lam_res_phi = c + d * phi

    if lam_res_phi > 0:
        gap_rel = gap_phi / lam_res_phi
    else:
        gap_rel = float('inf')

    phi_dist = abs(x_zero - phi) if x_zero != float('inf') else float('inf')

    print(f"{label:>15} {res_label:>12} {gap_coeff:>4}x+{gap_const:>4} "
          f"  x0={x_zero:10.6f}  dist={phi_dist:10.6f}  "
          f"gap={gap_phi:8.4f}  rel={gap_rel:8.6f}  [{note}]")

# CONTROL TEST 3: How special is 21/13 as an approximant?
print("\n--- CONTROL TEST 3: Continued fraction analysis ---")
print("Convergents of phi and their quality")
print()

def continued_fraction_convergents(x, n_terms=15):
    """Compute continued fraction convergents of x."""
    convergents = []
    a = int(np.floor(x))
    p_prev, p_curr = 1, a
    q_prev, q_curr = 0, 1
    convergents.append((p_curr, q_curr, abs(x - p_curr/q_curr)))
    rem = x - a
    for _ in range(n_terms - 1):
        if abs(rem) < 1e-15:
            break
        x_next = 1.0 / rem
        a = int(np.floor(x_next))
        p_prev, p_curr = p_curr, a * p_curr + p_prev
        q_prev, q_curr = q_curr, a * q_curr + q_prev
        convergents.append((p_curr, q_curr, abs(x - p_curr/q_curr)))
        rem = x_next - a
    return convergents

print("Continued fraction convergents of phi:")
print(f"{'p/q':>10} {'value':>12} {'|phi - p/q|':>14} {'q^2 * error':>12} {'26*p/q':>10} {'Near 42?':>10}")
print("-" * 75)

for p, q, err in continued_fraction_convergents(phi, 12):
    q2_err = q**2 * err
    val_26 = 26.0 * p / q
    near_42 = abs(val_26 - 42)
    marker = " <-- RESONANCE" if abs(near_42) < 1 else ""
    print(f"  {p:>4}/{q:<4} {p/q:12.8f} {err:14.10f} {q2_err:12.6f} "
          f"{val_26:10.4f} {near_42:10.4f}{marker}")

# KEY INSIGHT: The resonance occurs when 26*(p/q) ~ 42,
# i.e., p/q ~ 42/26 = 21/13.
# For phi, 21/13 IS a convergent (the 7th one).
# For other irrationals, 21/13 is NOT generally a convergent.

print("\n--- Same analysis for sqrt(2) ---")
print(f"{'p/q':>10} {'value':>12} {'|x - p/q|':>14} {'26*p/q':>10} {'Near 42?':>10}")
print("-" * 65)
for p, q, err in continued_fraction_convergents(np.sqrt(2), 12):
    val_26 = 26.0 * p / q
    near_42 = abs(val_26 - 42)
    marker = " <--" if near_42 < 2 else ""
    print(f"  {p:>4}/{q:<4} {p/q:12.8f} {err:14.10f} {val_26:10.4f} {near_42:10.4f}{marker}")

print("\n--- Same analysis for sqrt(3) ---")
print(f"{'p/q':>10} {'value':>12} {'|x - p/q|':>14} {'26*p/q':>10} {'Near 42?':>10}")
print("-" * 65)
for p, q, err in continued_fraction_convergents(np.sqrt(3), 12):
    val_26 = 26.0 * p / q
    near_42 = abs(val_26 - 42)
    marker = " <--" if near_42 < 2 else ""
    print(f"  {p:>4}/{q:<4} {p/q:12.8f} {err:14.10f} {val_26:10.4f} {near_42:10.4f}{marker}")


# ============================================================
# PART 2: TRIANGLE-EDGE RATIOS 7:12:9
# ============================================================
print("\n" + "=" * 70)
print("PART 2: TRIANGLE-EDGE PARTICIPATION RATIOS")
print("=" * 70)

print("""
From exp485:
  ACC: 336 edges participating in ACC triangles
  BCC: 576 edges participating in BCC triangles
  CCC: 432 edges participating in CCC triangles

GCD = 48, so ratios = 7:12:9

Question: 12 = degree, 9 = N_eig, 7 = n_23 (CKM). Structure or coincidence?
""")

# First: DERIVE the counts from first principles
print("--- Deriving the counts from combinatorics ---")
print()

# Total edges: 720
# Edge types by vertex pair:
#   AC: 96 edges (A-vertex to C-vertex)
#   BC: 192 edges (B-vertex to C-vertex)
#   CC: 432 edges (C-vertex to C-vertex)
# No AA, AB, BB edges (these vertex types are too far apart)

# Triangle types:
#   ACC: 240 triangles, each with 2 AC edges + 1 CC edge
#   BCC: 480 triangles, each with 2 BC edges + 1 CC edge
#   CCC: 480 triangles, each with 3 CC edges

# An edge "participates in" a triangle type if it belongs to at least
# one triangle of that type.

# ACC triangles have edges: 2 AC + 1 CC
# So edges participating in ACC = {AC edges in ACC triangles} union {CC edges in ACC triangles}

# How many of the 96 AC edges participate in ACC triangles?
# Every AC edge connects an A-vertex to a C-vertex.
# An ACC triangle has one A and two C vertices.
# The AC edge is one of the two AC edges of the triangle.

# Key: each A vertex has degree 12, with ALL neighbors being C-type
# (since A vertices are far from other A and B vertices)
# So each A vertex has 12 edges, all AC type.
# Each AC edge participates in some number of ACC triangles.

print("Edge type counts (from exp485):")
print("  AC: 96, BC: 192, CC: 432  (total 720)")
print()

# Per-vertex analysis
print("Per-vertex analysis:")
print(f"  8 A-vertices, degree 12 each -> 8*12/2 = 48? No, 8*12 = 96 edge-endpoints")
print(f"  But each AC edge has 1 A endpoint -> 96 AC edges (matches)")
print(f"  16 B-vertices, degree 12 each -> 16*12 = 192 edge-endpoints on B side")
print(f"  But each BC edge has 1 B endpoint -> 192 BC edges (matches)")
print(f"  96 C-vertices, degree 12 each -> 96*12 = 1152 edge-endpoints on C side")
print(f"  From AC: 96 C-endpoints, from BC: 192 C-endpoints, from CC: 2*432 = 864 C-endpoints")
print(f"  Total C-endpoints: 96+192+864 = {96+192+864} (should be 1152)")
print()

# Now: which edges participate in which triangle types?
# An edge participates in a triangle type if it's part of at least one such triangle.

# For ACC triangles (240 total):
# Each ACC triangle = 1 A + 2 C vertices
# Edges: 2 AC-type + 1 CC-type
# How many DISTINCT AC edges? All 96 AC edges (each A connects to 12 C's,
# and ACC triangles are formed by pairs of those C's that are also connected)

# For BCC triangles (480 total):
# Each BCC triangle = 1 B + 2 C vertices
# Edges: 2 BC-type + 1 CC-type

# For CCC triangles (480 total):
# Each CCC triangle = 3 C vertices
# Edges: 3 CC-type

# "Participation" means: edge e belongs to >= 1 triangle of type T
# An AC edge can only appear in ACC triangles (or ACB if those existed)
# Since there are no ABx triangles, AC edges ONLY participate in ACC.
# -> 96 AC edges participate in ACC

# A BC edge can only appear in BCC triangles.
# -> 192 BC edges participate in BCC

# A CC edge can appear in ACC, BCC, or CCC triangles.
# How many CC edges participate in each type?

# Total CC edges: 432
# CC edges in ACC: each ACC triangle contributes 1 CC edge.
#   But multiple ACC triangles may share a CC edge.
# CC edges in BCC: each BCC triangle contributes 1 CC edge.
# CC edges in CCC: each CCC triangle contributes 3 CC edges.

# Let's compute exactly.
print("--- Computing CC edge participation ---")

# We need the actual 600-cell data for this. Let me rebuild it.
# (Reuse the build code from exp485)

vertices = []
for i in range(4):
    for s in [1, -1]:
        v = [0.0]*4; v[i] = s; vertices.append(v)
for s0 in [1,-1]:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            for s3 in [1,-1]:
                vertices.append([s0/2, s1/2, s2/2, s3/2])
p_v, h_v, q_v = phi/2, 0.5, 1/(2*phi)
templates = [
    [p_v,h_v,q_v,0], [p_v,q_v,0,h_v], [p_v,0,h_v,q_v],
    [h_v,p_v,0,q_v], [h_v,q_v,p_v,0], [h_v,0,q_v,p_v],
    [q_v,p_v,h_v,0], [q_v,h_v,0,p_v], [q_v,0,p_v,h_v],
    [0,p_v,q_v,h_v], [0,h_v,p_v,q_v], [0,q_v,h_v,p_v]
]
for t in templates:
    for s0 in [1,-1]:
        for s1 in [1,-1]:
            for s2 in [1,-1]:
                for s3 in [1,-1]:
                    v = [s0*t[0], s1*t[1], s2*t[2], s3*t[3]]
                    if np.linalg.norm(v) > 0.01:
                        vertices.append(v)
unique_verts = []
for v in vertices:
    va = np.array(v)
    if not any(np.allclose(va, u, atol=1e-10) for u in unique_verts):
        unique_verts.append(va)
vertices = np.array(unique_verts)
N = len(vertices)
norms = np.linalg.norm(vertices, axis=1)
vertices = vertices / norms[:, np.newaxis]

vertex_type = []
for i in range(N):
    v = vertices[i]
    nonzero = np.sum(np.abs(v) > 0.01)
    if nonzero == 1 and np.max(np.abs(v)) > 0.95:
        vertex_type.append('A')
    elif np.allclose(np.sort(np.abs(v)), [0.5]*4, atol=0.05):
        vertex_type.append('B')
    else:
        vertex_type.append('C')

from scipy.spatial.distance import cdist
dist_matrix = cdist(vertices, vertices, 'euclidean')
np.fill_diagonal(dist_matrix, np.inf)
nn_dist = np.median(dist_matrix.min(axis=1))
adj = (np.abs(dist_matrix - nn_dist) < 0.02).astype(int)
np.fill_diagonal(adj, 0)

edges = []
edge_index = {}
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j]:
            idx = len(edges)
            edges.append((i,j))
            edge_index[(i,j)] = idx
            edge_index[(j,i)] = idx
E = len(edges)

edge_type_map = {}
for idx, (i,j) in enumerate(edges):
    t = ''.join(sorted(vertex_type[i] + vertex_type[j]))
    edge_type_map[idx] = t

faces = set()
for i in range(N):
    for j in range(i+1, N):
        if not adj[i,j]: continue
        common = np.where(adj[i] & adj[j])[0]
        for k in common:
            if k > j:
                faces.add((i,j,k))
faces = list(faces)

# Classify triangles and track which CC edges participate
cc_in_ACC = set()
cc_in_BCC = set()
cc_in_CCC = set()
all_in_ACC = set()
all_in_BCC = set()
all_in_CCC = set()

for (i,j,k) in faces:
    tri_t = ''.join(sorted(vertex_type[i] + vertex_type[j] + vertex_type[k]))
    tri_edges = []
    for v1,v2 in [(i,j),(i,k),(j,k)]:
        e_idx = edge_index[(min(v1,v2), max(v1,v2))]
        tri_edges.append(e_idx)

    if tri_t == 'ACC':
        all_in_ACC.update(tri_edges)
        for e in tri_edges:
            if edge_type_map[e] == 'CC':
                cc_in_ACC.add(e)
    elif tri_t == 'BCC':
        all_in_BCC.update(tri_edges)
        for e in tri_edges:
            if edge_type_map[e] == 'CC':
                cc_in_BCC.add(e)
    elif tri_t == 'CCC':
        all_in_CCC.update(tri_edges)
        for e in tri_edges:
            if edge_type_map[e] == 'CC':
                cc_in_CCC.add(e)

# Count AC edges in ACC, BC edges in BCC
ac_in_ACC = sum(1 for e in all_in_ACC if edge_type_map[e] == 'AC')
bc_in_BCC = sum(1 for e in all_in_BCC if edge_type_map[e] == 'BC')

print(f"\nCC edges participating in each triangle type:")
print(f"  CC in ACC: {len(cc_in_ACC)} out of 432")
print(f"  CC in BCC: {len(cc_in_BCC)} out of 432")
print(f"  CC in CCC: {len(cc_in_CCC)} out of 432")
print(f"  CC in ACC and BCC: {len(cc_in_ACC & cc_in_BCC)}")
print(f"  CC in ACC and CCC: {len(cc_in_ACC & cc_in_CCC)}")
print(f"  CC in BCC and CCC: {len(cc_in_BCC & cc_in_CCC)}")
print(f"  CC in ALL three: {len(cc_in_ACC & cc_in_BCC & cc_in_CCC)}")

print(f"\nTotal edge participation (ALL edge types):")
print(f"  Edges in ACC triangles: {len(all_in_ACC)}")
print(f"    AC: {ac_in_ACC}, CC: {len(cc_in_ACC)}")
print(f"  Edges in BCC triangles: {len(all_in_BCC)}")
print(f"    BC: {bc_in_BCC}, CC: {len(cc_in_BCC)}")
print(f"  Edges in CCC triangles: {len(all_in_CCC)}")
print(f"    CC: {len(cc_in_CCC)}")

# Verify the 336:576:432 from exp485
print(f"\nVerification: {len(all_in_ACC)}:{len(all_in_BCC)}:{len(all_in_CCC)}")

from math import gcd
g = gcd(gcd(len(all_in_ACC), len(all_in_BCC)), len(all_in_CCC))
r_ACC = len(all_in_ACC) // g
r_BCC = len(all_in_BCC) // g
r_CCC = len(all_in_CCC) // g
print(f"GCD = {g}, Simplified ratios: {r_ACC}:{r_BCC}:{r_CCC}")

# Now derive the counts analytically
print("\n--- Analytical derivation ---")
print()

# ACC triangle: 1 A-vertex + 2 C-vertices, all mutually adjacent
# Edges: 2 AC-type + 1 CC-type
# 240 ACC triangles, 3 edges each -> 720 edge-appearances
# But edges are shared: each edge appears in multiple triangles

# Count how many triangles each edge type appears in:
edge_tri_count = defaultdict(lambda: defaultdict(int))
for (i,j,k) in faces:
    tri_t = ''.join(sorted(vertex_type[i] + vertex_type[j] + vertex_type[k]))
    for v1,v2 in [(i,j),(i,k),(j,k)]:
        e_idx = edge_index[(min(v1,v2), max(v1,v2))]
        edge_tri_count[e_idx][tri_t] += 1

# Average number of triangles per edge, by edge type and triangle type
print("Average triangles per edge, by type:")
for etype in ['AC', 'BC', 'CC']:
    edges_of_type = [e for e in range(E) if edge_type_map[e] == etype]
    counts = defaultdict(list)
    for e in edges_of_type:
        for tt in ['ACC', 'BCC', 'CCC']:
            counts[tt].append(edge_tri_count[e].get(tt, 0))

    print(f"  {etype} edges ({len(edges_of_type)}):")
    for tt in ['ACC', 'BCC', 'CCC']:
        vals = counts[tt]
        if sum(vals) > 0:
            print(f"    in {tt}: mean={np.mean(vals):.2f}, "
                  f"min={min(vals)}, max={max(vals)}, "
                  f"total={sum(vals)}")

# Verify: total triangle-edge incidences
print(f"\nEdge-face incidence check:")
print(f"  240 ACC * 3 edges = {240*3} incidences")
print(f"  480 BCC * 3 edges = {480*3} incidences")
print(f"  480 CCC * 3 edges = {480*3} incidences")
print(f"  Total: {240*3 + 480*3 + 480*3} (should = 1200*3 = 3600 since each face has 3 edges)")

# Derive: for AC edges
# Each AC edge connects A to C.
# An ACC triangle through this edge needs a second C neighbor of both A and C.
# A has 12 neighbors (all C), C has 12 neighbors.
# Common C-neighbors of A and C (excluding each other) = # ACC triangles through this edge

# Each A vertex is in how many ACC triangles?
a_verts = [i for i in range(N) if vertex_type[i] == 'A']
print(f"\nPer-A-vertex triangle count:")
for av in a_verts[:3]:
    tri_count = sum(1 for (i,j,k) in faces
                    if av in (i,j,k) and
                    ''.join(sorted(vertex_type[i]+vertex_type[j]+vertex_type[k])) == 'ACC')
    print(f"  Vertex {av}: {tri_count} ACC triangles")

# Each A vertex should be in 240/8 * (fraction) ACC triangles
# 240 ACC triangles, each has 1 A vertex -> 240 A-vertex appearances / 8 A-vertices = 30
print(f"  Expected: 240/8 = {240/8} ACC triangles per A-vertex")

# Similarly for B vertices in BCC
b_verts = [i for i in range(N) if vertex_type[i] == 'B']
for bv in b_verts[:3]:
    tri_count = sum(1 for (i,j,k) in faces
                    if bv in (i,j,k) and
                    ''.join(sorted(vertex_type[i]+vertex_type[j]+vertex_type[k])) == 'BCC')
    print(f"  B-vertex {bv}: {tri_count} BCC triangles")
print(f"  Expected: 480/16 = {480/16} BCC triangles per B-vertex")

# Now derive 336, 576, 432:
print("\n--- Deriving 336:576:432 ---")

# ACC uses: 96 AC edges + CC edges in ACC
# BCC uses: 192 BC edges + CC edges in BCC
# CCC uses: only CC edges in CCC

print(f"  ACC participation = {ac_in_ACC} AC + {len(cc_in_ACC)} CC = {ac_in_ACC + len(cc_in_ACC)}")
print(f"  BCC participation = {bc_in_BCC} BC + {len(cc_in_BCC)} CC = {bc_in_BCC + len(cc_in_BCC)}")
print(f"  CCC participation = {len(cc_in_CCC)} CC")

# Check: are ALL CC edges in some triangle type?
cc_edges = set(e for e in range(E) if edge_type_map[e] == 'CC')
cc_in_any = cc_in_ACC | cc_in_BCC | cc_in_CCC
print(f"\n  CC edges in at least one triangle type: {len(cc_in_any)} / {len(cc_edges)}")
print(f"  CC edges NOT in any triangle: {len(cc_edges - cc_in_any)}")

# Analyze the CC edge overlap structure
print(f"\n  CC edge Venn diagram:")
only_ACC = cc_in_ACC - cc_in_BCC - cc_in_CCC
only_BCC = cc_in_BCC - cc_in_ACC - cc_in_CCC
only_CCC = cc_in_CCC - cc_in_ACC - cc_in_BCC
ACC_BCC = (cc_in_ACC & cc_in_BCC) - cc_in_CCC
ACC_CCC = (cc_in_ACC & cc_in_CCC) - cc_in_BCC
BCC_CCC = (cc_in_BCC & cc_in_CCC) - cc_in_ACC
all_three = cc_in_ACC & cc_in_BCC & cc_in_CCC
print(f"  Only ACC: {len(only_ACC)}")
print(f"  Only BCC: {len(only_BCC)}")
print(f"  Only CCC: {len(only_CCC)}")
print(f"  ACC+BCC only: {len(ACC_BCC)}")
print(f"  ACC+CCC only: {len(ACC_CCC)}")
print(f"  BCC+CCC only: {len(BCC_CCC)}")
print(f"  All three: {len(all_three)}")
total = (len(only_ACC) + len(only_BCC) + len(only_CCC) +
         len(ACC_BCC) + len(ACC_CCC) + len(BCC_CCC) + len(all_three))
print(f"  Total: {total} (should be {len(cc_edges)})")

# ============================================================
# PART 3: ARE 7, 12, 9 MEANINGFUL?
# ============================================================
print("\n" + "=" * 70)
print("PART 3: INTERPRETING THE RATIOS 7:12:9")
print("=" * 70)

# Each A vertex has degree 12, all neighbors are C-type
# Each A vertex is in 30 ACC triangles
# Each B vertex is in 30 BCC triangles
# Each C vertex is in:
#   some ACC (as one of the 2 C's)
#   some BCC (as one of the 2 C's)
#   some CCC (as all 3 C's)

# ACC: 240 tri * 2 C-endpoints / 96 C-verts = 5.0 ACC tri per C-vertex
# BCC: 480 tri * 2 C-endpoints / 96 C-verts = 10.0 BCC tri per C-vertex
# CCC: 480 tri * 3 C-endpoints / 96 C-verts = 15.0 CCC tri per C-vertex
# Total: 5+10+15 = 30 triangles per C-vertex?
# Should be: each vertex in 1200*3/120 = 30 tri. Yes!

print("\nTriangles per C-vertex:")
print(f"  ACC: 240*2/96 = {240*2/96}")
print(f"  BCC: 480*2/96 = {480*2/96}")
print(f"  CCC: 480*3/96 = {480*3/96}")
print(f"  Total: {240*2/96 + 480*2/96 + 480*3/96} (should be 30)")

# Now: participation counts derivation
# 336 = 96 (AC) + 240 (CC in ACC)
# 576 = 192 (BC) + 384 (CC in BCC)
# 432 = CC in CCC (all CC edges? let's see)

print(f"\nHypothesis: 336 = 96 + 240, 576 = 192 + 384, 432 = 432")
print(f"  96 + {len(cc_in_ACC)} = {96 + len(cc_in_ACC)}")
print(f"  192 + {len(cc_in_BCC)} = {192 + len(cc_in_BCC)}")
print(f"  {len(cc_in_CCC)}")

# Express in terms of theory constants
print(f"\nIn terms of a1={a1}, b1={b1}, N={N}:")
print(f"  AC edges: 96 = 8*12 = |cross-polytope| * degree = {8*12}")
print(f"  BC edges: 192 = 16*12 = |tesseract| * degree = {16*12}")
print(f"  CC edges: 432 = 96*12/2 - 96 - 192 = ... or 720-96-192 = {720-96-192}")
print()

# The key: 7:12:9
print(f"Ratios 7:12:9:")
print(f"  7 + 12 + 9 = {7+12+9}")
print(f"  7 * 12 * 9 = {7*12*9}")
print(f"  7 = ?")
print(f"    n_23 (CKM exponent): YES")
print(f"    a1 + 2 = 7: YES")
print(f"  12 = degree = 2*b1: YES")
print(f"  9 = N_eig = distinct adjacency eigenvalues: YES")
print(f"    Also: 9 = rank(E8) + 1, or dim(coexact rep)")
print()

# But IS this structural or numerological?
# Test: do these numbers arise from the COMBINATORICS or from the PHYSICS?

# The numbers 96, 192, 432 come from vertex counts 8, 16, 96 and degree 12.
# 336 = AC + CC_in_ACC = 96 + 240 = 96 + 240
# 576 = BC + CC_in_BCC = 192 + 384

# 240 = ? Number of CC edges in ACC triangles
# 384 = ? Number of CC edges in BCC triangles

# Each ACC triangle has exactly 1 CC edge. 240 triangles -> 240 CC edge appearances.
# But edges can be shared! Need to count DISTINCT CC edges.
# From data: |cc_in_ACC| CC edges participate in ACC triangles

# How many ACC triangles share each CC edge?
cc_acc_sharing = defaultdict(int)
for (i,j,k) in faces:
    tri_t = ''.join(sorted(vertex_type[i]+vertex_type[j]+vertex_type[k]))
    if tri_t == 'ACC':
        for v1,v2 in [(i,j),(i,k),(j,k)]:
            e_idx = edge_index[(min(v1,v2), max(v1,v2))]
            if edge_type_map[e_idx] == 'CC':
                cc_acc_sharing[e_idx] += 1

sharing_counts = list(cc_acc_sharing.values())
print(f"ACC triangles per CC edge (among CC edges in ACC):")
print(f"  Values: {sorted(set(sharing_counts))}")
print(f"  Counts: {dict(sorted([(v, sharing_counts.count(v)) for v in set(sharing_counts)]))}")
print(f"  Total CC edges in ACC: {len(sharing_counts)}")
print(f"  Total ACC tri-edge incidences: {sum(sharing_counts)} (should be 240)")

cc_bcc_sharing = defaultdict(int)
for (i,j,k) in faces:
    tri_t = ''.join(sorted(vertex_type[i]+vertex_type[j]+vertex_type[k]))
    if tri_t == 'BCC':
        for v1,v2 in [(i,j),(i,k),(j,k)]:
            e_idx = edge_index[(min(v1,v2), max(v1,v2))]
            if edge_type_map[e_idx] == 'CC':
                cc_bcc_sharing[e_idx] += 1

sharing_bcc = list(cc_bcc_sharing.values())
print(f"\nBCC triangles per CC edge (among CC edges in BCC):")
print(f"  Values: {sorted(set(sharing_bcc))}")
print(f"  Counts: {dict(sorted([(v, sharing_bcc.count(v)) for v in set(sharing_bcc)]))}")
print(f"  Total CC edges in BCC: {len(sharing_bcc)}")
print(f"  Total BCC tri-edge incidences: {sum(sharing_bcc)} (should be 480)")

cc_ccc_sharing = defaultdict(int)
for (i,j,k) in faces:
    tri_t = ''.join(sorted(vertex_type[i]+vertex_type[j]+vertex_type[k]))
    if tri_t == 'CCC':
        for v1,v2 in [(i,j),(i,k),(j,k)]:
            e_idx = edge_index[(min(v1,v2), max(v1,v2))]
            cc_ccc_sharing[e_idx] += 1

sharing_ccc = list(cc_ccc_sharing.values())
print(f"\nCCC triangles per CC edge (among CC edges in CCC):")
print(f"  Values: {sorted(set(sharing_ccc))}")
print(f"  Counts: {dict(sorted([(v, sharing_ccc.count(v)) for v in set(sharing_ccc)]))}")
print(f"  Total CC edges in CCC: {len(sharing_ccc)}")
print(f"  Total CCC tri-edge incidences: {sum(sharing_ccc)} (should be 1440)")

# ============================================================
# FINAL VERDICT
# ============================================================
print("\n" + "=" * 70)
print("FINAL VERDICT")
print("=" * 70)

# Compute: for which irrationals does the gap vanish?
# gap = 26x - 42 = 0 <=> x = 21/13
# For phi: |gap| = |26*phi - 42| = 0.069
# For sqrt(2): |gap| = |26*1.414 - 42| = |36.77 - 42| = 5.23 (HUGE)
# For sqrt(3): |gap| = |26*1.732 - 42| = |45.03 - 42| = 3.03 (BIG)

print(f"""
QUESTION 1: Is the Fibonacci near-resonance trivial?

VERDICT: PARTIALLY TRIVIAL, PARTIALLY STRUCTURAL.

The near-resonance delta = 26*phi - 42 ~ 0.069 occurs because:
  - 42/26 = 21/13 is a continued fraction convergent of phi
  - This is a MATHEMATICAL FACT about phi, not specific to 600-cell

BUT: The coefficients 26 and 42 ARE specific to 600-cell!
  - 26 = a1^2 + 1 (appears as denominator of sin^2 theta_W)
  - The threshold comes from 4*(12-6*phi) where 12=degree, 6=b1
  - The resonance comes from b1+2*phi where b1=a1+1

For OTHER irrationals:
  sqrt(2): |26*sqrt(2) - 42| = {abs(26*np.sqrt(2) - 42):.3f} (FAR from resonance)
  sqrt(3): |26*sqrt(3) - 42| = {abs(26*np.sqrt(3) - 42):.3f} (FAR from resonance)
  e:       |26*e - 42| = {abs(26*np.e - 42):.3f} (far)
  pi:      |26*pi - 42| = {abs(26*np.pi - 42):.3f} (far)

The near-resonance requires BOTH:
  (a) phi as the irrational (from golden ratio / icosahedral symmetry)
  (b) coefficients 26 and 42 from the 600-cell structure

Neither alone suffices. So the result is STRUCTURAL but NOT surprising
given that the entire framework is built on phi.

HONEST STATUS: STRUCTURAL (expected consequence of phi + 600-cell)
Not a deep new result, but a consistency check that the pair creation
threshold is naturally compatible with the spectral structure.
""")

print(f"""
QUESTION 2: Are the ratios 7:12:9 meaningful?

From the computation:
  336 = 96 (AC edges) + {len(cc_in_ACC)} (CC edges in ACC triangles)
  576 = 192 (BC edges) + {len(cc_in_BCC)} (CC edges in BCC triangles)
  432 = {len(cc_in_CCC)} (CC edges in CCC triangles)

The numbers 96, 192, 432 come from vertex type counts 8, 16, 96
times degree 12. These are pure COMBINATORICS of the 600-cell.

Whether the simplified ratio 7:12:9 = (a1+2):degree:(N_eig) is
COINCIDENCE or STRUCTURE requires checking:
""")

# Final check: are these truly combinatorial identities?
print(f"  336 = 96 + {len(cc_in_ACC)}")
print(f"      96 = |A|*degree = 8*12")
print(f"      {len(cc_in_ACC)} = CC edges adjacent to A-vertices")
print(f"      336/48 = 7 = a1+2? {336/48 == a1+2}")
print(f"  576 = 192 + {len(cc_in_BCC)}")
print(f"      192 = |B|*degree = 16*12")
print(f"      {len(cc_in_BCC)} = CC edges in BCC triangles")
print(f"      576/48 = 12 = degree? {576/48 == 12}")
print(f"  432 = {len(cc_in_CCC)}")
print(f"      432/48 = 9 = N_eig? {432/48 == 9}")
print()

# Check: does 7:12:9 follow from 8:16:96 vertex counts and degree 12?
# That would make it purely combinatorial.
# 8:16:96 = 1:2:12. These give AC:BC:CC = 1:2:9 (edges), times degree/2.
# But the PARTICIPATION counts mix edge types (AC,BC,CC are all different sizes).
# So the ratio 7:12:9 is NOT simply proportional to vertex counts.

# Is 7 = a1+2 = n_23?
print(f"  7 = a1+2 = {a1+2}: algebraic identity with a1")
print(f"  7 = n_23 (CKM): CKM exponent for V_cb")
print(f"  Connection: n_23 is DERIVED as a1+2 in exp453? ", end="")
# Check: in the CKM derivation, n_23 = n_ut*n_db - delta = 4*2-1 = 7
# And n_ut = a1-1 = 4, n_db = 2, delta = 1
# So n_23 = (a1-1)*2 - 1 = 2*a1-3 = 7 (for a1=5)
# Also a1+2 = 7 for a1=5. But 2*a1-3 = a1+2 only when a1=5!
print(f"n_23 = 2*a1-3 = {2*a1-3}, a1+2 = {a1+2}. Equal iff a1=5!")

print(f"""
VERDICT on 7:12:9:
  7 = a1+2 = 2*a1-3 (ONLY for a1=5). Combinatorial origin plausible.
  12 = degree (trivially structural)
  9 = N_eig (structural but WHY here?)

The ratio encodes the D4<2I vertex decomposition structure.
It's COMBINATORIAL (from vertex type counts) but the specific values
being a1+2, degree, N_eig is a CONSISTENCY CHECK, not a new result.

HONEST STATUS: COMBINATORIAL CONSISTENCY (not a deep connection to CKM)
""")

print("=" * 70)
print("EXP-486 COMPLETE")
print("=" * 70)

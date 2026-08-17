"""
EXP-160: E8 Dimensional Reduction and Gravity

E8 lives in R^8. The icosian construction splits R^8 = R^4_phys x R^4_int.
The 600-cell sits in R^4_int (fiber). The physical R^4 is the base.

Key idea: S^7 Hopf fibration S^3 -> S^7 -> S^4.
E8 roots live on S^7 (after rescaling). The fiber S^3 gives gauge (D4).
The base S^4 should give gravity!

Questions:
1. How does the E8 metric decompose under the 4+4 split?
2. What is the "gravitational" sector in the split?
3. Can we extract Schwarzschild-like behavior from the T copy?
4. What role does the Galois conjugate play?
"""

import numpy as np
from itertools import combinations
from collections import defaultdict

print("=" * 70)
print("EXP-160: E8 Dimensional Reduction and Gravity")
print("=" * 70)

# =====================================================================
# Generate 600-cell and E8
# =====================================================================
phi = (1 + np.sqrt(5)) / 2
phi_prime = (1 - np.sqrt(5)) / 2  # = 1 - phi = -1/phi

# Generate 600-cell vertices
vertices_600 = []
# Type A
for i in range(4):
    for s in [1, -1]:
        v = [0.0]*4; v[i] = float(s); vertices_600.append(v)
# Type B
for s0 in [1,-1]:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            for s3 in [1,-1]:
                vertices_600.append([s0/2, s1/2, s2/2, s3/2])
# Type C
vals = [phi/2, 0.5, 1/(2*phi), 0.0]
even_perms = [[0,1,2,3],[0,2,3,1],[0,3,1,2],[1,0,3,2],[1,2,0,3],[1,3,2,0],
              [2,0,1,3],[2,1,3,0],[2,3,0,1],[3,0,2,1],[3,1,0,2],[3,2,1,0]]
for perm in even_perms:
    base = [vals[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for s0 in [1,-1]:
        for s1 in [1,-1]:
            for s2 in [1,-1]:
                v = [0.0]*4
                for idx, ni in enumerate(nz):
                    v[ni] = base[ni] * [s0,s1,s2][idx]
                vertices_600.append(v)

# Deduplicate
verts = []
for v in vertices_600:
    v = np.array(v)
    if not any(np.linalg.norm(v-w) < 1e-10 for w in verts):
        verts.append(v)
S = np.array(verts)
print(f"|S| = {len(S)}")

# Decompose each vertex: v_i = a_i + b_i * phi
def decompose_phi(v):
    """Decompose v = a + b*phi, return (a, b)"""
    a = np.zeros(4)
    b = np.zeros(4)
    for i in range(4):
        x = v[i]
        # x = a_i + b_i * phi
        # Try all reasonable (a,b) pairs
        best_err = 1e10
        for bi in np.arange(-2, 2.5, 0.5):
            ai = x - bi * phi
            if abs(ai - round(2*ai)/2) < 1e-8:
                err = abs(x - round(2*ai)/2 - bi*phi)
                if err < best_err:
                    best_err = err
                    a[i] = round(2*ai)/2
                    b[i] = bi
    return a, b

# Decompose all vertices
A_coeffs = []
B_coeffs = []
for v in S:
    a, b = decompose_phi(v)
    A_coeffs.append(a)
    B_coeffs.append(b)
    assert abs(np.linalg.norm(a + b*phi) - 1.0) < 1e-8, f"Decomposition failed for {v}"

A_coeffs = np.array(A_coeffs)
B_coeffs = np.array(B_coeffs)

# Build T = phi' * S (Galois conjugate via icosian multiplication)
# T has decomposition (a', b') = (a-b, -a) from (a,b)
T_A = A_coeffs - B_coeffs  # a' = a - b
T_B = -A_coeffs              # b' = -a
T_4d = T_A + T_B * phi     # Reconstruct in R^4

print(f"|T| = {len(T_4d)}")

# Also build sigma(S) = Galois conjugate (replace phi -> phi')
sigma_S = A_coeffs + B_coeffs * phi_prime
print(f"|sigma(S)| = {len(sigma_S)}")

# Verify T and sigma(S) are both on unit sphere
T_norms = np.linalg.norm(T_4d, axis=1)
sigma_norms = np.linalg.norm(sigma_S, axis=1)
print(f"T norms: min={T_norms.min():.6f}, max={T_norms.max():.6f}")
print(f"sigma(S) norms: min={sigma_norms.min():.6f}, max={sigma_norms.max():.6f}")

# =====================================================================
# PART 1: E8 in R^8 via Cholesky
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: E8 Root System in R^8")
print("=" * 70)

# Cholesky embedding: (a_i, b_i) -> (a_i + b_i, b_i)
# Gram matrix G = [[1,1],[1,2]], Cholesky L = [[1,0],[1,1]]
def cholesky_embed(a_coeffs, b_coeffs):
    """Embed (a,b) decomposition into R^8 via Cholesky"""
    N = len(a_coeffs)
    result = np.zeros((N, 8))
    for i in range(4):
        result[:, 2*i] = a_coeffs[:, i] + b_coeffs[:, i]  # a+b
        result[:, 2*i+1] = b_coeffs[:, i]                   # b
    return result

S_8d = cholesky_embed(A_coeffs, B_coeffs)
T_8d = cholesky_embed(T_A, T_B)

E8_roots = np.vstack([S_8d, T_8d])
print(f"S in R^8: {S_8d.shape}")
print(f"T in R^8: {T_8d.shape}")
print(f"E8 total: {E8_roots.shape}")

# Check norms
S_8d_norms = np.linalg.norm(S_8d, axis=1)
T_8d_norms = np.linalg.norm(T_8d, axis=1)
print(f"S norms in R^8: {S_8d_norms.min():.6f} to {S_8d_norms.max():.6f}")
print(f"T norms in R^8: {T_8d_norms.min():.6f} to {T_8d_norms.max():.6f}")

# Verify no overlap in R^8
overlaps_8d = 0
for i in range(120):
    for j in range(120):
        if np.linalg.norm(S_8d[i] - T_8d[j]) < 1e-8:
            overlaps_8d += 1
print(f"S-T overlaps in R^8: {overlaps_8d}")

# =====================================================================
# PART 2: The 4+4 Split of R^8
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: The 4+4 Decomposition of R^8")
print("=" * 70)

# R^8 = R^4_physical x R^4_internal
# The Cholesky embedding interleaves: (a0+b0, b0, a1+b1, b1, a2+b2, b2, a3+b3, b3)
# Physical = even indices (a_i + b_i) = "total" coordinates
# Internal = odd indices (b_i) = "phi-component" coordinates

# Alternative split: separate a and b entirely
# Physical: (a0, a1, a2, a3) - rational part
# Internal: (b0, b1, b2, b3) - irrational part

# Physical projection (a-part):
phys_proj = A_coeffs  # (a0, a1, a2, a3)
int_proj = B_coeffs   # (b0, b1, b2, b3)

print("Split 1: a (rational) vs b (irrational) components")
print(f"  Physical (a): shape {phys_proj.shape}")
print(f"  Internal (b): shape {int_proj.shape}")

# For Type A+B (24 vertices with b=0):
mask_AB = np.all(np.abs(B_coeffs) < 1e-10, axis=1)
mask_C = ~mask_AB
print(f"\n  Type A+B (b=0): {mask_AB.sum()} vertices")
print(f"  Type C (b!=0): {mask_C.sum()} vertices")
print(f"  Type A+B have ZERO internal projection (Galois-fixed)")

# =====================================================================
# PART 3: S^7 Hopf Fibration
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: S^7 Hopf Fibration S^3 -> S^7 -> S^4")
print("=" * 70)

print("""
The Hopf fibration S^3 -> S^7 -> S^4 comes from the quaternionic structure.
View R^8 = H^2 (pairs of quaternions).
Each E8 root (q1, q2) in H^2 with |q1|^2 + |q2|^2 = 1 lives on S^7.

The Hopf map: (q1, q2) -> q1 * q2^{-1} (when q2 != 0)
maps S^7 to S^4 = HP^1 (quaternionic projective line).

The fiber over each point is S^3 = {(q1*g, q2*g) : g in unit quaternions}.
""")

# Interpret E8 roots as pairs of quaternions (q1, q2) in H^2
# S_8d[i] = (a0+b0, b0, a1+b1, b1, a2+b2, b2, a3+b3, b3)
# Rearrange as q1 = (x0, x1, x2, x3), q2 = (x4, x5, x6, x7)

# For the Cholesky embedding, the natural quaternion split is:
# q1 = (S_8d[0], S_8d[1], S_8d[2], S_8d[3])
# q2 = (S_8d[4], S_8d[5], S_8d[6], S_8d[7])

S_q1 = S_8d[:, :4]
S_q2 = S_8d[:, 4:]
T_q1 = T_8d[:, :4]
T_q2 = T_8d[:, 4:]

S_q1_norms = np.linalg.norm(S_q1, axis=1)
S_q2_norms = np.linalg.norm(S_q2, axis=1)

print(f"|q1| for S: min={S_q1_norms.min():.6f}, max={S_q1_norms.max():.6f}")
print(f"|q2| for S: min={S_q2_norms.min():.6f}, max={S_q2_norms.max():.6f}")
print(f"|q1|^2 + |q2|^2 for S: min={(S_q1_norms**2 + S_q2_norms**2).min():.6f}, "
      f"max={(S_q1_norms**2 + S_q2_norms**2).max():.6f}")

# Hopf map: (q1, q2) -> q1/q2 (quaternion division)
# But we need to handle q2 = 0 cases (north pole of S^4)

def quat_mult(a, b):
    """Quaternion multiplication a*b"""
    r = np.zeros(4)
    r[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
    r[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2]
    r[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1]
    r[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0]
    return r

def quat_conj(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

def quat_inv(a):
    return quat_conj(a) / np.dot(a, a)

# Compute Hopf projections for S roots
hopf_S = []
for i in range(120):
    q1, q2 = S_q1[i], S_q2[i]
    norm_q2 = np.linalg.norm(q2)
    if norm_q2 > 1e-10:
        # Hopf map: q1 * q2^{-1}
        h = quat_mult(q1, quat_inv(q2))
        hopf_S.append(h)
    else:
        hopf_S.append(np.array([np.inf, 0, 0, 0]))  # North pole

hopf_S = np.array(hopf_S)
finite_mask = np.all(np.isfinite(hopf_S), axis=1)
print(f"\nHopf projections (S):")
print(f"  Finite: {finite_mask.sum()}, Infinite (north pole): {(~finite_mask).sum()}")

if finite_mask.sum() > 0:
    hopf_norms = np.linalg.norm(hopf_S[finite_mask], axis=1)
    print(f"  |hopf(S)|: min={hopf_norms.min():.4f}, max={hopf_norms.max():.4f}")
    print(f"  Distinct |hopf| values: {len(np.unique(np.round(hopf_norms, 4)))}")

    # The Hopf image lives in R^4, representing points on S^4 via stereographic projection
    # Count distinct images
    hopf_finite = hopf_S[finite_mask]
    distinct = []
    for h in hopf_finite:
        if not any(np.linalg.norm(h - d) < 1e-4 for d in distinct):
            distinct.append(h)
    print(f"  Distinct Hopf images: {len(distinct)}")

# =====================================================================
# PART 4: Alternative 4+4 Split: Original v vs Galois v'
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: Physical vs Internal via Galois")
print("=" * 70)

print("""
The most natural 4+4 split for the 600-cell is:
  Physical = v (original 600-cell vertex in R^4)
  Internal = sigma(v) = Galois conjugate (replace phi -> phi')

Each vertex lives at (v, sigma(v)) in R^4 x R^4 = R^8.
The 24 gauge vertices have v = sigma(v) (on the "diagonal" of R^8).
The 96 matter vertices have v != sigma(v) (off-diagonal).

For the E8 roots:
  S: (v, sigma(v)) - pairs of Galois conjugates
  T: (phi'*v, phi'*sigma(v)) = (phi'*v, phi*v) - swapped!

Key insight: S lives on the "same side" (v, sigma(v))
             T lives on the "flipped side" (phi'*v, phi*v)
""")

# Construct (v, sigma(v)) pairs
V_phys_int = np.zeros((120, 8))
for i in range(120):
    V_phys_int[i, :4] = S[i]  # physical = v
    V_phys_int[i, 4:] = sigma_S[i]  # internal = sigma(v)

# Check: are these the same as the Cholesky embedding?
# Cholesky: (a+b, b, ...) for each coordinate pair
# Galois: (a+b*phi, a+b*phi') for each full vector

# The relationship:
# v_i = a_i + b_i*phi  (physical)
# sigma(v)_i = a_i + b_i*phi' = a_i + b_i*(1-phi) = (a_i+b_i) - b_i*phi (internal)

# So the pair is: ((a+b*phi)_i, (a+b*phi')_i) for i=0..3
# This is NOT the same as the Cholesky embedding ((a+b)_i, b_i)

# But it's related by a linear transformation:
# [v, sigma(v)] = [a+b*phi, a+b*phi'] = [[phi, phi'], [1,1]] . [b, a]
# Hmm, let me think...

# v = a + b*phi
# sigma(v) = a + b*phi'
# So [v, sigma(v)] = [1, phi; 1, phi'] . [a; b]
# = A . [a; b] where A = [[1, phi], [1, phi']]

# The Cholesky embedding is [a+b, b] = [[1,1],[0,1]] . [a, b]
# These are related by A * inv([[1,1],[0,1]]) = A * [[1,-1],[0,1]]

# In any case, both are valid R^8 embeddings of the same data.

# The physical vs internal split in the (v, sigma(v)) basis:
print("In the (v, sigma(v)) basis:")
print("  Physical = first 4 coordinates (v)")
print("  Internal = last 4 coordinates (sigma(v))")

# Norms in each sector
phys_norms = np.linalg.norm(V_phys_int[:, :4], axis=1)
int_norms = np.linalg.norm(V_phys_int[:, 4:], axis=1)
print(f"\n  Physical |v|: all = {phys_norms[0]:.6f}")
print(f"  Internal |sigma(v)|: all = {int_norms[0]:.6f}")
print(f"  |v|^2 + |sigma(v)|^2 = {(phys_norms**2 + int_norms**2)[0]:.6f}")

# =====================================================================
# PART 5: Metric from E8 Lattice
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: Effective Metric from E8 Structure")
print("=" * 70)

print("""
In Kaluza-Klein theory, the metric on M^4 x K (internal space K) is:
  ds^2 = g_mn dx^m dx^n + g_ab dy^a dy^b + 2 A_m^a dx^m dy^a

where:
  g_mn = 4D metric (gravity)
  g_ab = metric on K (determines gauge couplings)
  A_m^a = gauge fields

For our setup:
  K = S^3 (internal, 600-cell)
  g_ab determined by 600-cell geometry (gives gauge couplings)
  g_mn = 4D metric (what we want!)
  A_m^a = mixed terms (gauge-gravity mixing)

The E8 metric (Killing form) in R^8 is just the Euclidean metric
(since E8 is simply-laced). The 4+4 split gives:

ds^2_E8 = ds^2_phys + ds^2_int + cross_terms

The cross terms between physical and internal encode the gauge fields.
The physical part encodes gravity.
""")

# Compute the E8 neighbor graph
print("Building E8 neighbor graph...")
E8_all = np.vstack([S_8d, T_8d])
E8_dots = E8_all @ E8_all.T
# E8 neighbors: inner product = 1/2 (for norm-1 vectors)
# After rescaling to norm sqrt(2), neighbors have dot = 1
# For our norm-1 vectors, neighbors have dot = 1/2
E8_adj = defaultdict(set)
E8_edges = []
for i in range(240):
    for j in range(i+1, 240):
        if abs(E8_dots[i,j] - 0.5) < 1e-6:
            E8_adj[i].add(j)
            E8_adj[j].add(i)
            E8_edges.append((i,j))

E8_degrees = [len(E8_adj[i]) for i in range(240)]
print(f"E8 edges: {len(E8_edges)}")
print(f"E8 degree: min={min(E8_degrees)}, max={max(E8_degrees)}")

# Classify E8 edges by type: S-S, S-T, T-T
SS_edges = [(i,j) for i,j in E8_edges if i < 120 and j < 120]
TT_edges = [(i,j) for i,j in E8_edges if i >= 120 and j >= 120]
ST_edges = [(i,j) for i,j in E8_edges if (i < 120) != (j < 120)]

print(f"\nS-S edges: {len(SS_edges)} (gauge + matter within S)")
print(f"T-T edges: {len(TT_edges)} (gauge + matter within T)")
print(f"S-T edges: {len(ST_edges)} (cross-sector)")
print(f"Total: {len(SS_edges) + len(TT_edges) + len(ST_edges)}")

# The cross-sector (S-T) edges are the "gravitational" connections
# in the KK picture - they mix the physical and internal sectors.

# What is the degree distribution within each sector?
S_degree_in_S = [sum(1 for j in E8_adj[i] if j < 120) for i in range(120)]
S_degree_in_T = [sum(1 for j in E8_adj[i] if j >= 120) for i in range(120)]
T_degree_in_S = [sum(1 for j in E8_adj[i] if j < 120) for i in range(120, 240)]
T_degree_in_T = [sum(1 for j in E8_adj[i] if j >= 120) for i in range(120, 240)]

print(f"\nPer-vertex degree decomposition:")
print(f"  S vertex: {S_degree_in_S[0]} within S + {S_degree_in_T[0]} to T = {S_degree_in_S[0]+S_degree_in_T[0]} total")
print(f"  T vertex: {T_degree_in_S[0]} to S + {T_degree_in_T[0]} within T = {T_degree_in_S[0]+T_degree_in_T[0]} total")

# Check if all S vertices have the same decomposition
S_in_S_unique = set(S_degree_in_S)
S_in_T_unique = set(S_degree_in_T)
print(f"  S degree within S: {S_in_S_unique}")
print(f"  S degree to T: {S_in_T_unique}")

# =====================================================================
# PART 6: S-T Edge Analysis (Gravitational Sector)
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: S-T Edge Analysis (Cross-Sector)")
print("=" * 70)

# For each S-T edge, compute the physical and internal components
print(f"Total S-T edges: {len(ST_edges)}")

# Edge vectors in R^8
ST_vectors = []
for i, j in ST_edges:
    v = E8_all[j] - E8_all[i]
    ST_vectors.append(v)
ST_vectors = np.array(ST_vectors)

# Decompose into physical (first 4) and internal (last 4) components
ST_phys = ST_vectors[:, :4]  # Cholesky "physical"
ST_int = ST_vectors[:, 4:]   # Cholesky "internal"

ST_phys_norms = np.linalg.norm(ST_phys, axis=1)
ST_int_norms = np.linalg.norm(ST_int, axis=1)

print(f"S-T edge physical component |delta_phys|: {ST_phys_norms.min():.4f} to {ST_phys_norms.max():.4f}")
print(f"S-T edge internal component |delta_int|: {ST_int_norms.min():.4f} to {ST_int_norms.max():.4f}")

# Do the same for S-S edges (should be mostly internal)
SS_vectors = np.array([E8_all[j] - E8_all[i] for i,j in SS_edges])
SS_phys_norms = np.linalg.norm(SS_vectors[:, :4], axis=1)
SS_int_norms = np.linalg.norm(SS_vectors[:, 4:], axis=1)

print(f"\nS-S edge physical component: {SS_phys_norms.min():.4f} to {SS_phys_norms.max():.4f}")
print(f"S-S edge internal component: {SS_int_norms.min():.4f} to {SS_int_norms.max():.4f}")

# The physical component of S-T edges represents the "gravitational displacement"
# in the KK picture.

# =====================================================================
# PART 7: Type-Specific Analysis
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: S-T Connections by Vertex Type")
print("=" * 70)

# Classify S vertices
type_A = list(range(8))
type_B = list(range(8, 24))
type_C = list(range(24, 120))

# For each type, how many T-neighbors?
for label, indices in [("A", type_A), ("B", type_B), ("C", type_C)]:
    t_degrees = [S_degree_in_T[i] for i in indices]
    s_degrees = [S_degree_in_S[i] for i in indices]
    print(f"Type {label} ({len(indices)} vertices):")
    print(f"  Degree in S: {set(s_degrees)}")
    print(f"  Degree in T: {set(t_degrees)}")

# =====================================================================
# PART 8: The KK Metric Components
# =====================================================================
print("\n" + "=" * 70)
print("PART 8: Kaluza-Klein Metric Components")
print("=" * 70)

# In KK, the 8D metric decomposes as:
# G_MN = [[g_mn + A_m^a A_n^b g_ab, A_m^a g_ab],
#          [A_n^b g_ab, g_ab]]
# where M,N = 0..7, m,n = 0..3 (physical), a,b = 4..7 (internal)

# For our lattice, the "metric" is encoded in the inner products of roots.
# The E8 Cartan matrix (in our basis) is a specific 8x8 matrix.

# Let's compute the "effective metric" by looking at the inner products
# restricted to physical and internal sectors.

# The E8 inner product between roots r_i and r_j is:
# <r_i, r_j> = sum_M r_i^M r_j^M (Euclidean in Cholesky basis)

# Effective physical metric: G_phys = sum of outer products of physical components
# For nearest-neighbor roots:
G_phys_avg = np.zeros((4, 4))
G_int_avg = np.zeros((4, 4))
G_cross_avg = np.zeros((4, 4))

for i in range(240):
    for j in E8_adj[i]:
        dp = E8_all[j, :4] - E8_all[i, :4]  # physical displacement
        di = E8_all[j, 4:] - E8_all[i, 4:]  # internal displacement
        G_phys_avg += np.outer(dp, dp)
        G_int_avg += np.outer(di, di)
        G_cross_avg += np.outer(dp, di)

n_edges_total = 2 * len(E8_edges)  # each edge counted twice
G_phys_avg /= n_edges_total
G_int_avg /= n_edges_total
G_cross_avg /= n_edges_total

print("Effective physical metric (from E8 neighbor displacements):")
print(np.round(G_phys_avg, 6))

print("\nEffective internal metric:")
print(np.round(G_int_avg, 6))

print("\nCross (gauge field) components:")
print(np.round(G_cross_avg, 6))

# Check if physical metric is proportional to identity (isotropic)
print(f"\nPhysical metric diagonal: {np.diag(G_phys_avg)}")
print(f"Physical metric off-diag max: {np.max(np.abs(G_phys_avg - np.diag(np.diag(G_phys_avg)))):.8f}")

# =====================================================================
# PART 9: Regge-Like Action on E8
# =====================================================================
print("\n" + "=" * 70)
print("PART 9: Regge-Like Action from E8 Structure")
print("=" * 70)

# The E8 root system defines a crystallographic structure.
# The "curvature" can be measured by the deviation from flatness.
# For the 600-cell (S), we know the deficit angles are 0 (flat S^3).
# But the COMBINED S+T structure might have non-zero curvature
# when projected to the physical 4D.

# Compute the "angular deficit" in the physical projection:
# For each triple of S-T-S vertices forming a triangle,
# compute the deviation from flatness.

# Count S-T-S paths (S vertex -> T neighbor -> another S neighbor of T)
print("Analyzing S-T-S triangles (KK graviton paths)...")
n_STS_triangles = 0
physical_angles = []

for s1 in range(120):
    for t in E8_adj[s1]:
        if t >= 120:  # t is in T
            for s2 in E8_adj[t]:
                if s2 < 120 and s2 > s1:  # s2 in S, avoid double counting
                    if s2 in E8_adj[s1]:  # closed triangle
                        n_STS_triangles += 1
                        # Physical projection of this triangle
                        p1 = E8_all[s1, :4]
                        p2 = E8_all[t, :4]
                        p3 = E8_all[s2, :4]
                        # Angle at vertex t (in physical space)
                        v1 = p1 - p2
                        v3 = p3 - p2
                        n1, n3 = np.linalg.norm(v1), np.linalg.norm(v3)
                        if n1 > 1e-10 and n3 > 1e-10:
                            cos_a = np.dot(v1, v3) / (n1 * n3)
                            physical_angles.append(np.arccos(np.clip(cos_a, -1, 1)))
                    if n_STS_triangles >= 10000:
                        break
        if n_STS_triangles >= 10000:
            break
    if n_STS_triangles >= 10000:
        break

print(f"S-T-S closed triangles found: {n_STS_triangles}")
if physical_angles:
    physical_angles = np.array(physical_angles)
    print(f"Physical angle at T-vertex:")
    print(f"  Mean: {np.degrees(np.mean(physical_angles)):.2f} deg")
    print(f"  Std: {np.degrees(np.std(physical_angles)):.2f} deg")
    print(f"  Distinct: {len(np.unique(np.round(physical_angles, 4)))}")

# =====================================================================
# PART 10: Summary
# =====================================================================
print("\n" + "=" * 70)
print("PART 10: Summary")
print("=" * 70)

print(f"""
E8 STRUCTURE:
  |S| = 120, |T| = 120, |E8| = 240
  S-S edges: {len(SS_edges)}
  T-T edges: {len(TT_edges)}
  S-T edges: {len(ST_edges)} (cross-sector = "gravitational")

DEGREE DECOMPOSITION:
  Each S vertex: {S_degree_in_S[0]} within S + {S_degree_in_T[0]} to T
  Each T vertex: {T_degree_in_S[0]} to S + {T_degree_in_T[0]} within T

PHYSICAL METRIC (from E8 displacements):
  Diagonal: {np.diag(G_phys_avg).round(6)}
  Off-diagonal max: {np.max(np.abs(G_phys_avg - np.diag(np.diag(G_phys_avg)))):.8f}
  => {"ISOTROPIC" if np.max(np.abs(G_phys_avg - np.diag(np.diag(G_phys_avg)))) < 1e-4 else "ANISOTROPIC"}

CROSS-SECTOR (gauge field) COMPONENTS:
  Max: {np.max(np.abs(G_cross_avg)):.8f}
  => {"ZERO (no gauge-gravity mixing)" if np.max(np.abs(G_cross_avg)) < 1e-4 else "NON-ZERO (gauge-gravity mixing present)"}

KEY FINDINGS:
1. E8 naturally splits into S (120) + T (120) with cross-edges
2. Each vertex has BOTH S-type and T-type neighbors
3. The cross-sector (S-T) edges represent gauge-gravity mixing in KK
4. The physical metric from E8 is {"isotropic" if np.max(np.abs(G_phys_avg - np.diag(np.diag(G_phys_avg)))) < 1e-4 else "anisotropic"}
""")

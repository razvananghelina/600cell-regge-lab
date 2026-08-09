"""
EXP-263: ATTACKING THE FOUR OPEN PROBLEMS
==========================================
1. Derive phi - 8*alpha from spectral data
2. Build explicit D_F on McKay graph (off-diagonal Yukawa)
3. Find the 1+2+9 edge decomposition criterion
4. Cosmological constant sign and value

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from itertools import permutations, product
from collections import defaultdict, Counter
import time

PHI = (1 + np.sqrt(5)) / 2
phi = PHI
a_1 = 5
b_1 = 6
alpha_em = 1.0 / 137.035999177
alpha_s = 1.0 / (2 * phi**3)

print("=" * 72)
print("EXP-263: FOUR OPEN PROBLEMS")
print("=" * 72)

# =====================================================================
# Build 600-cell (compact, reused from exp262)
# =====================================================================
print("\n--- Building 600-cell ---")
t0 = time.time()
verts_set = set()
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0]*4; v[i] = s; verts_set.add(tuple(v))
for signs in product([0.5, -0.5], repeat=4):
    verts_set.add(signs)
base = [PHI/2, 0.5, 1/(2*PHI), 0.0]
even_perms = [p for p in permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]) % 2 == 0]
for perm in even_perms:
    coords = [base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if abs(coords[i]) > 1e-12]
    for signs in product([1, -1], repeat=len(nz)):
        v = list(coords)
        for idx, s in zip(nz, signs):
            v[idx] *= s
        verts_set.add(tuple(round(x, 10) for x in v))
verts = np.array(sorted(verts_set))
Nv = len(verts)

dots = verts @ verts.T
edges = []
adj_list = defaultdict(set)
for i in range(Nv):
    for j in range(i+1, Nv):
        if abs(dots[i, j] - PHI/2) < 0.001:
            edges.append((i, j))
            adj_list[i].add(j)
            adj_list[j].add(i)
Ne = len(edges)
edge_to_idx = {}
for idx, (i, j) in enumerate(edges):
    edge_to_idx[(i, j)] = idx
    edge_to_idx[(j, i)] = idx

triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j:
                    triangles.append((i, j, k))
Nf = len(triangles)

tetrahedra = []
for i in range(Nv):
    ni = adj_list[i]
    for j in ni:
        if j > i:
            cij = ni & adj_list[j]
            for k in cij:
                if k > j:
                    for l in cij & adj_list[k]:
                        if l > k:
                            tetrahedra.append((i, j, k, l))
Nc = len(tetrahedra)
N_total = Nv + Ne + Nf + Nc
print(f"  V={Nv}, E={Ne}, F={Nf}, C={Nc}, Total={N_total}")

# Classify vertex types: A (axis), B (half-integer), C (golden ratio)
def vertex_type(v, tol=1e-6):
    a, b, c, d = v
    if sum(abs(x) > tol for x in v) == 1 and any(abs(abs(x)-1) < tol for x in v):
        return 'A'
    if all(abs(abs(x) - 0.5) < tol for x in v):
        return 'B'
    return 'C'

vtypes = [vertex_type(verts[i]) for i in range(Nv)]
nA = vtypes.count('A')
nB = vtypes.count('B')
nC = vtypes.count('C')
print(f"  Vertex types: A={nA}, B={nB}, C={nC}")
print(f"  Time: {time.time()-t0:.1f}s")

# =====================================================================
# PROBLEM 3: The 1+2+9 edge decomposition
# =====================================================================
print("\n" + "=" * 72)
print("PROBLEM 3: THE 1+2+9 EDGE DECOMPOSITION")
print("=" * 72)

# Hypothesis: for each type C vertex, its 12 neighbors split as:
# 1 type A + 2 type B + 9 type C = 12

print("\n  Neighbor type distribution for each vertex type:")
for vt in ['A', 'B', 'C']:
    # Take the first vertex of this type
    v_indices = [i for i in range(Nv) if vtypes[i] == vt]
    # Count neighbor types (check all vertices of this type to confirm uniformity)
    type_counts = Counter()
    for vi in v_indices:
        nbr_types = Counter(vtypes[n] for n in adj_list[vi])
        for key in nbr_types:
            type_counts[key] = nbr_types[key]  # same for all (vertex-transitive within type)
    # Display for first vertex
    vi = v_indices[0]
    nbr_types = Counter(vtypes[n] for n in adj_list[vi])
    print(f"    Type {vt} vertex (count={len(v_indices)}): neighbors = {dict(sorted(nbr_types.items()))}")

# Check uniformity: do ALL vertices of the same type have the same neighbor distribution?
print("\n  Uniformity check:")
for vt in ['A', 'B', 'C']:
    v_indices = [i for i in range(Nv) if vtypes[i] == vt]
    distributions = set()
    for vi in v_indices:
        nbr_types = tuple(sorted(Counter(vtypes[n] for n in adj_list[vi]).items()))
        distributions.add(nbr_types)
    print(f"    Type {vt}: {len(distributions)} distinct neighbor distributions")

# Count edges by endpoint types
edge_type_counts = Counter()
for i, j in edges:
    t1, t2 = vtypes[i], vtypes[j]
    key = tuple(sorted([t1, t2]))
    edge_type_counts[key] += 1

print(f"\n  Edge counts by endpoint types:")
for key in sorted(edge_type_counts.keys()):
    print(f"    {key[0]}-{key[1]}: {edge_type_counts[key]} edges")

# The 1+3+8 GAUGE decomposition via CC* criterion
# For each type C vertex v, among its 9 type-C neighbors,
# find the UNIQUE one (CC*) that shares a type-A common neighbor
# with v but NOT a type-B common neighbor.

print(f"\n  CC* CRITERION:")
print(f"  For each type-C vertex v, among its 9 type-C neighbors n,")
print(f"  find the UNIQUE n such that:")
print(f"    - v and n share at least one type-A common neighbor")
print(f"    - v and n share NO type-B common neighbor")

cc_star_count = Counter()  # how many CC* neighbors does each C-vertex have?
cc_star_edges = set()  # set of CC* edges

for vi in range(Nv):
    if vtypes[vi] != 'C':
        continue
    c_neighbors = [n for n in adj_list[vi] if vtypes[n] == 'C']
    cc_star_for_vi = []
    for n in c_neighbors:
        # Common neighbors of vi and n
        common = adj_list[vi] & adj_list[n]
        common_A = [c for c in common if vtypes[c] == 'A']
        common_B = [c for c in common if vtypes[c] == 'B']
        if len(common_A) > 0 and len(common_B) == 0:
            cc_star_for_vi.append(n)
            edge_key = (min(vi, n), max(vi, n))
            cc_star_edges.add(edge_key)
    cc_star_count[len(cc_star_for_vi)] += 1

print(f"\n  CC* neighbor count distribution: {dict(cc_star_count)}")
print(f"  Total CC* edges: {len(cc_star_edges)}")

# So the decomposition of edges is:
# U(1): edges from C to A (shared A > 0 trivially, but A-C not C-C)
# SU(2): edges from C to B PLUS CC* edges
# SU(3): edges from C to C that are NOT CC*

n_AC = edge_type_counts.get(('A', 'C'), 0)
n_BC = edge_type_counts.get(('B', 'C'), 0)
n_CC = edge_type_counts.get(('C', 'C'), 0)
n_CC_star = len(cc_star_edges)
n_CC_regular = n_CC - n_CC_star

print(f"\n  GAUGE SECTOR DECOMPOSITION:")
print(f"    A-C edges (U(1) candidates):     {n_AC}")
print(f"    B-C edges (SU(2) candidates):     {n_BC}")
print(f"    CC* edges (SU(2) extension):      {n_CC_star}")
print(f"    CC regular edges (SU(3)):          {n_CC_regular}")
print(f"    Total: {n_AC} + {n_BC} + {n_CC_star} + {n_CC_regular} = {n_AC+n_BC+n_CC_star+n_CC_regular}")

# The 1+3+8 interpretation:
# dim 1 (U(1)): A-C edges -> 96 edges / 2 = 48 per direction?
# dim 3 (SU(2)): B-C + CC* edges -> 192 + CC* = ?
# dim 8 (SU(3)): regular CC edges -> n_CC_regular

# Try another approach: the gauge group dimensions match edge ORBITS
# Under the vertex stabilizer, edges from a given vertex decompose as 1+2+9
# The 1+2+9 -> 1+3+8 regrouping assigns:
#   1 A-neighbor -> 1 gauge boson (U(1))
#   2 B-neighbors + 1 special C-neighbor -> 3 gauge bosons (SU(2))
#   8 remaining C-neighbors -> 8 gauge bosons (SU(3))

# For each C-vertex, find which of its 9 C-neighbors is CC*
print(f"\n  Per-vertex check (type C vertex 24):")
vi = [i for i in range(Nv) if vtypes[i] == 'C'][0]
print(f"  Vertex {vi} (type C), neighbors:")
for n in sorted(adj_list[vi]):
    t = vtypes[n]
    is_cc_star = ""
    if t == 'C':
        common = adj_list[vi] & adj_list[n]
        cA = sum(1 for c in common if vtypes[c] == 'A')
        cB = sum(1 for c in common if vtypes[c] == 'B')
        cC = sum(1 for c in common if vtypes[c] == 'C')
        if cA > 0 and cB == 0:
            is_cc_star = " <-- CC*"
        print(f"    n={n:3d} type={t} common=(A:{cA}, B:{cB}, C:{cC-2}){is_cc_star}")
    else:
        print(f"    n={n:3d} type={t}")

# Local decomposition: 1(A) + 2(B) + 1(CC*) + 8(CC_regular) = 1 + 3 + 8 = 12
print(f"\n  LOCAL GAUGE DECOMPOSITION (per type-C vertex):")
print(f"    1 type-A neighbor  -> U(1)")
print(f"    2 type-B neighbors -> SU(2) (partial)")
print(f"    1 CC* neighbor     -> SU(2) (completing dim 3)")
print(f"    8 regular C-nbrs   -> SU(3)")
print(f"    Total: 1 + (2+1) + 8 = 1 + 3 + 8 = 12")

# Global edge count for gauge sectors
print(f"\n  GLOBAL EDGE COUNT:")
print(f"    U(1) sector:  {n_AC} edges (A-C)")
print(f"    SU(2) sector: {n_BC + n_CC_star} edges (B-C + CC*)")
print(f"    SU(3) sector: {n_CC_regular} edges (regular C-C)")
print(f"    Total: {n_AC + n_BC + n_CC_star + n_CC_regular}")

# =====================================================================
# PROBLEM 1: Derive phi - 8*alpha from spectral data
# =====================================================================
print("\n" + "=" * 72)
print("PROBLEM 1: DERIVING m_H/m_W = phi - 8*alpha")
print("=" * 72)

# Key insight: 8*alpha = 2/(a_1*phi^4) to leading order
# From alpha equation: 4*a_1*phi^4*alpha ~ 1
# So 8*alpha ~ 2/(a_1*phi^4)
# And m_H/m_W = phi - 2/(a_1*phi^4) to leading order

alpha_leading = 1.0 / (4 * a_1 * phi**4)
alpha_exact = alpha_em

print(f"\n  Alpha equation: 2*pi*alpha^2 - 4*a_1*phi^4*alpha + 1 = 0")
print(f"  Leading order: alpha ~ 1/(4*a_1*phi^4) = 1/{4*a_1*phi**4:.4f}")
print(f"  Leading alpha = {alpha_leading:.8f} = 1/{1/alpha_leading:.4f}")
print(f"  Exact alpha   = {alpha_exact:.8f} = 1/{1/alpha_exact:.4f}")
print(f"  Ratio: {alpha_leading/alpha_exact:.8f}")

print(f"\n  GEOMETRIC FORM of Higgs mass formula:")
val_geom = phi - 2.0 / (a_1 * phi**4)
val_alpha = phi - 8 * alpha_exact
val_alpha_theory = phi - 8 * alpha_leading

print(f"    phi - 2/(a_1*phi^4) = {val_geom:.8f}")
print(f"    phi - 8*alpha_exact = {val_alpha:.8f}")
print(f"    Difference: {abs(val_geom - val_alpha):.2e}")
print(f"    Relative diff: {abs(val_geom - val_alpha)/val_alpha*100:.4f}%")

# Now derive this GEOMETRICALLY:
# m_H/m_W = phi * (1 - 2/(a_1 * phi^5))
# phi^5 = a_1*phi + (a_1-2) [Fibonacci: F(5)=5=a_1, F(4)=3=a_1-2]
# So a_1*phi^5 = a_1*(a_1*phi + a_1-2) = a_1^2*phi + a_1*(a_1-2)
#              = 25*phi + 15

print(f"\n  DERIVATION:")
print(f"  phi^5 = F(5)*phi + F(4) = {a_1}*phi + {a_1-2} = {a_1*phi + a_1-2:.6f}")
print(f"  a_1*phi^5 = {a_1}^2*phi + {a_1}*{a_1-2} = {a_1**2}*phi + {a_1*(a_1-2)}")
print(f"           = {a_1**2*phi + a_1*(a_1-2):.4f}")
print(f"  2/(a_1*phi^5) = {2/(a_1*phi**5):.8f}")
print(f"  8*alpha       = {8*alpha_exact:.8f}")

# The O(alpha^2) correction:
# Exact alpha satisfies: 2*pi*alpha^2 - 4*a_1*phi^4*alpha + 1 = 0
# alpha = [4*a_1*phi^4 - sqrt((4*a_1*phi^4)^2 - 8*pi)] / (4*pi)
# Let x = 4*a_1*phi^4. Then alpha = [x - sqrt(x^2 - 8*pi)] / (4*pi)
# To next order: alpha = 1/x + 2*pi/x^3 + ...
# So 8*alpha = 8/x + 16*pi/x^3 + ...
# = 2/(a_1*phi^4) + 16*pi/(4*a_1*phi^4)^3 + ...
# = 2/(a_1*phi^4) + pi/(8*a_1^3*phi^12) + ...

x = 4 * a_1 * phi**4
correction_2 = np.pi / (8 * a_1**3 * phi**12)
alpha_order2 = 1/x + 2*np.pi/x**3

print(f"\n  Perturbative expansion:")
print(f"    8*alpha = 2/(a_1*phi^4) + pi/(8*a_1^3*phi^12) + ...")
print(f"    Leading: {2/(a_1*phi**4):.8f}")
print(f"    Next order: {8*correction_2:.2e}")
print(f"    Sum (order 2): {2/(a_1*phi**4) + 8*correction_2:.8f}")
print(f"    Exact 8*alpha: {8*alpha_exact:.8f}")

# PHYSICAL INTERPRETATION
print(f"\n  PHYSICAL INTERPRETATION:")
print(f"  m_H/m_W = phi - 2/(a_1*phi^4)")
print(f"")
print(f"  phi = golden ratio = base DSI scaling of the 600-cell")
print(f"  2/(a_1*phi^4) = quantum correction from finite geometry")
print(f"")
print(f"  In Z[phi]: phi^4 = 3*phi + 2 = (a_1-2)*phi + (a_1-3)")
print(f"  So a_1*phi^4 = a_1*(a_1-2)*phi + a_1*(a_1-3) = {a_1*(a_1-2)}*phi + {a_1*(a_1-3)}")
print(f"             = {a_1*(a_1-2)*phi + a_1*(a_1-3):.4f}")

# Connection to Connes framework
print(f"\n  CONNECTION TO CONNES FORMULA:")
print(f"  Standard Connes: m_H/m_W ~ sqrt(2) * m_t/m_W  (at GUT scale)")
print(f"  600-cell:        m_H/m_W = phi - 2/(a_1*phi^4)")
print(f"")
print(f"  Ratio: (phi - 2/(a_1*phi^4)) / (sqrt(2)*m_t_bare/m_W)")
m_t_bare = 0.51099895 * phi**26  # m_e * phi^26
m_W = 80379.0
ratio_connes = val_geom / (np.sqrt(2) * m_t_bare / m_W)
print(f"  = {val_geom:.6f} / (sqrt(2)*{m_t_bare:.0f}/{m_W:.0f})")
print(f"  = {val_geom:.6f} / {np.sqrt(2)*m_t_bare/m_W:.6f}")
print(f"  = {ratio_connes:.6f}")
print(f"")
print(f"  This ratio ~ m_W / (sqrt(2)*m_t_bare) = {m_W/(np.sqrt(2)*m_t_bare):.6f}")
print(f"  = 1/(sqrt(2)*phi^26/phi^25) * v/m_e  ... (scale dependent)")

# Check: is the ratio related to any Z[phi] quantity?
print(f"\n  Looking for Z[phi] expression of the ratio:")
print(f"  ratio = {ratio_connes:.6f}")
print(f"  phi^(-1) = {1/phi:.6f}")
print(f"  phi^(-2) = {1/phi**2:.6f}")
print(f"  1/phi^2 * sqrt(a_1) = {np.sqrt(a_1)/phi**2:.6f}")
# 0.6401 ~ 0.6395 ~ hmm
print(f"  (a_1-1)/(2*a_1+1) = {(a_1-1)/(2*a_1+1):.6f}")  # = 4/11
print(f"  b_1/(a_1+b_1-2) = {b_1/(a_1+b_1-2):.6f}")  # = 6/9 = 2/3

# STATUS for problem 1
print(f"\n  STATUS: PARTIALLY DERIVED")
print(f"  The formula m_H/m_W = phi - 8*alpha is equivalent to:")
print(f"    m_H/m_W = phi - 2/(a_1*phi^4)  [to O(alpha^2)]")
print(f"  This is PURELY GEOMETRIC (only a_1=5 and phi).")
print(f"  Physical meaning: phi (DSI base) minus quantum correction 2/(a_1*phi^4)")
print(f"  The correction 2/(a_1*phi^4) = 8*alpha comes from the alpha equation.")
print(f"  WHAT REMAINS: derive WHY m_H/m_W = phi - correction from spectral action")

# =====================================================================
# PROBLEM 2: Explicit D_F on McKay graph
# =====================================================================
print("\n" + "=" * 72)
print("PROBLEM 2: EXPLICIT D_F ON McKAY GRAPH")
print("=" * 72)

# The affine E8 adjacency
A_E8 = np.zeros((9, 9))
E8_edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (4,8)]
for i, j in E8_edges:
    A_E8[i,j] = 1
    A_E8[j,i] = 1

# Node assignment: 0(e), 1(d), 2(u), 3(s), 4(c), 5(b), 6(t), 7(tau), 8(mu)
node_names = ['e', 'd', 'u', 's', 'c', 'b', 't', 'tau', 'mu']
n_exp = np.array([0, 5, 3, 11, 16, 19, 26, 17, 11], dtype=float)  # mass exponents
coxeter = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3], dtype=float)

# Method 1: D_F = diagonal (mass basis, trivial)
D_F_diag = np.diag(n_exp)
print(f"  Method 1: D_F = diag(n_f) [mass basis]")
print(f"  Eigenvalues: {sorted(n_exp)}")
print(f"  This is trivial - no off-diagonal mixing information.")

# Method 2: D_F from adjacency + diagonal
# D_F = diag(n_f) + w * A_E8 where w is a coupling parameter
# This puts Yukawa mixing on the E8 edges
# The eigenvalues should be CLOSE to n_f for small w

print(f"\n  Method 2: D_F = diag(n_f) + w * A_E8")
print(f"  Scanning w to find eigenvalue structure:")

for w in [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]:
    D_F = np.diag(n_exp) + w * A_E8
    evals = np.sort(np.linalg.eigvalsh(D_F))
    max_shift = max(abs(e - n) for e, n in zip(sorted(evals), sorted(n_exp)))
    print(f"    w={w:4.1f}: evals={np.round(evals, 2)}, max_shift={max_shift:.2f}")

# Method 3: Coxeter-weighted adjacency
# D_F = diag(n_f) + w * C_adj where C_adj_ij = sqrt(cox_i * cox_j) * A_E8_ij
# This weights edges by the geometric mean of Coxeter labels

print(f"\n  Method 3: D_F = diag(n_f) + w * Coxeter-weighted adjacency")
C_adj = np.zeros((9, 9))
for i, j in E8_edges:
    C_adj[i, j] = np.sqrt(coxeter[i] * coxeter[j])
    C_adj[j, i] = C_adj[i, j]

for w in [0.0, 0.1, 0.5, 1.0]:
    D_F = np.diag(n_exp) + w * C_adj
    evals = np.sort(np.linalg.eigvalsh(D_F))
    print(f"    w={w:4.1f}: evals={np.round(evals, 2)}")

# Method 4: The NATURAL D_F from the mass formula
# The mass formula is m_f = m_e * phi^(5a + 6b)
# In the McKay graph basis, the (a,b) quantum numbers define D_F
# D_F should encode the generation structure: how generations mix

# The CKM matrix (leading order from 600-cell)
theta_12 = np.arctan(phi**(-3))
theta_23 = np.arctan(phi**(-7))
theta_13 = np.arctan(phi**(-12))
delta_CP = np.arctan(np.sqrt(5))

# Standard CKM parametrization
c12, s12 = np.cos(theta_12), np.sin(theta_12)
c23, s23 = np.cos(theta_23), np.sin(theta_23)
c13, s13 = np.cos(theta_13), np.sin(theta_13)
cd, sd = np.cos(delta_CP), np.sin(delta_CP)

V_CKM = np.array([
    [c12*c13,                    s12*c13,                    s13*np.exp(-1j*delta_CP)],
    [-s12*c23-c12*s23*s13*np.exp(1j*delta_CP), c12*c23-s12*s23*s13*np.exp(1j*delta_CP), s23*c13],
    [s12*s23-c12*c23*s13*np.exp(1j*delta_CP), -c12*s23-s12*c23*s13*np.exp(1j*delta_CP), c23*c13]
])

print(f"\n  Method 4: CKM-rotated mass matrix")
print(f"  |V_CKM| =")
for row in np.abs(V_CKM):
    print(f"    [{row[0]:.6f}  {row[1]:.6f}  {row[2]:.6f}]")

# The quark mass matrices in the interaction basis:
# M_up = V_CKM^dag * diag(m_u, m_c, m_t) * V_CKM (if down is diagonal)
# Or: M_down = V_CKM * diag(m_d, m_s, m_b) * V_CKM^dag (if up is diagonal)

# Using mass exponents:
n_up = np.array([3.0, 16.0, 26.0])    # u, c, t
n_down = np.array([5.0, 11.0, 19.0])  # d, s, b
n_lep = np.array([0.0, 11.0, 17.0])   # e, mu, tau

# Quark mass matrix in interaction basis (down-type diagonal convention)
M_up_interaction = V_CKM.conj().T @ np.diag(n_up) @ V_CKM
print(f"\n  Up-quark mass matrix in interaction basis (exponents):")
print(f"  M_up_int = V_CKM^dag * diag({n_up}) * V_CKM")
for row in M_up_interaction:
    print(f"    [{row[0].real:+8.4f} {row[1].real:+8.4f} {row[2].real:+8.4f}]")

# Build the full 9x9 D_F in the interaction basis
# Block structure: leptons (0,7,8) = (e,tau,mu), down quarks (1,3,5) = (d,s,b),
#                  up quarks (2,4,6) = (u,c,t)
#
# Reorder to sector blocks: [e,mu,tau, u,c,t, d,s,b]
# Node order: 0(e), 8(mu), 7(tau), 2(u), 4(c), 6(t), 1(d), 3(s), 5(b)
reorder = [0, 8, 7, 2, 4, 6, 1, 3, 5]  # node indices in sector order

# In the mass basis (sector-ordered), D_F is diagonal
n_sector = n_exp[reorder]
print(f"\n  Mass exponents in sector order: {n_sector}")
print(f"  Sectors: lep={n_sector[:3]}, up={n_sector[3:6]}, down={n_sector[6:]}")

# Full 9x9 D_F in interaction basis
D_F_full = np.zeros((9, 9), dtype=complex)
# Lepton block: diagonal (PMNS ignored for simplicity)
D_F_full[:3, :3] = np.diag(n_sector[:3])
# Up quark block: CKM-rotated
D_F_full[3:6, 3:6] = V_CKM.conj().T @ np.diag(n_sector[3:6]) @ V_CKM
# Down quark block: diagonal (convention: down = mass basis)
D_F_full[6:, 6:] = np.diag(n_sector[6:])

evals_full = np.sort(np.linalg.eigvalsh(D_F_full))
print(f"\n  D_F eigenvalues (CKM-rotated up quarks):")
print(f"    {np.round(evals_full, 4)}")
print(f"    Original: {sorted(n_exp)}")
print(f"    Match: {np.allclose(sorted(evals_full), sorted(n_exp), atol=0.01)}")

# The off-diagonal elements encode CKM mixing
print(f"\n  Off-diagonal elements (CKM mixing in up-quark block):")
up_block = D_F_full[3:6, 3:6]
for i in range(3):
    for j in range(3):
        if i != j:
            val = up_block[i, j]
            print(f"    ({['u','c','t'][i]},{['u','c','t'][j]}): {val.real:+.6f} + {val.imag:+.6f}i  |{abs(val):.6f}|")

# Project D_F back onto the McKay graph basis (original node ordering)
inv_reorder = [0]*9
for i, r in enumerate(reorder):
    inv_reorder[r] = i

D_F_mckay = np.zeros((9, 9), dtype=complex)
for i in range(9):
    for j in range(9):
        D_F_mckay[i, j] = D_F_full[inv_reorder[i], inv_reorder[j]]

# Check: which off-diagonal elements of D_F_mckay are non-zero?
# Compare with E8 adjacency
print(f"\n  D_F on McKay graph (off-diagonal |D_F_ij|):")
print(f"  {'':>6}", end="")
for name in node_names:
    print(f" {name:>6}", end="")
print()
for i in range(9):
    print(f"  {node_names[i]:>6}", end="")
    for j in range(9):
        val = abs(D_F_mckay[i, j])
        marker = "*" if A_E8[i, j] > 0 else " "
        if val > 0.001 and i != j:
            print(f" {val:5.3f}{marker}", end="")
        elif i == j:
            print(f" [{val:4.1f}]", end="")
        else:
            print(f"     . ", end="")
    print()
print(f"  (* = E8 adjacency edge)")

# How much mixing is ON the E8 edges vs OFF?
on_E8 = 0
off_E8 = 0
for i in range(9):
    for j in range(9):
        if i == j:
            continue
        val = abs(D_F_mckay[i, j])**2
        if A_E8[i, j] > 0:
            on_E8 += val
        else:
            off_E8 += val

total_offdiag = on_E8 + off_E8
print(f"\n  Mixing on E8 edges:  {on_E8:.4f} ({on_E8/total_offdiag*100:.1f}%)")
print(f"  Mixing off E8 edges: {off_E8:.4f} ({off_E8/total_offdiag*100:.1f}%)")

# =====================================================================
# PROBLEM 4: Cosmological constant
# =====================================================================
print("\n" + "=" * 72)
print("PROBLEM 4: COSMOLOGICAL CONSTANT")
print("=" * 72)

# Build Hodge Laplacians for eigenvalue computation
print("\n  Building Hodge Laplacians...")
d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = +1.0

face_to_idx = {tri: idx for idx, tri in enumerate(triangles)}
d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

d2 = np.zeros((Nc, Nf))
for c_idx, (i, j, k, l) in enumerate(tetrahedra):
    for face, sign in [((j,k,l),+1), ((i,k,l),-1), ((i,j,l),+1), ((i,j,k),-1)]:
        f_idx = face_to_idx.get(face)
        if f_idx is not None:
            d2[c_idx, f_idx] = sign

Delta_0 = d0.T @ d0
Delta_1 = d0 @ d0.T + d1.T @ d1
Delta_2 = d1 @ d1.T + d2.T @ d2
Delta_3 = d2 @ d2.T

print("  Computing eigenvalues...")
evals_0 = np.sort(np.linalg.eigvalsh(Delta_0))
evals_1 = np.sort(np.linalg.eigvalsh(Delta_1))
evals_2 = np.sort(np.linalg.eigvalsh(Delta_2))
evals_3 = np.sort(np.linalg.eigvalsh(Delta_3))
all_evals = np.sort(np.concatenate([evals_0, evals_1, evals_2, evals_3]))

# The cosmological constant in the spectral action
# S = Tr(f(D^2/Lambda^2))
# For the product geometry: D_total^2 = D_M^2 (x) 1 + 1 (x) D_F^2
# The cosmological constant term is the f(0) * dim(H) piece

# In the standard spectral action:
# Lambda_cc = f(0) * c_0 / (16*pi*G) (the vacuum energy)
# The SIGN depends on what you're computing.

# Method 1: from the spectral action expansion
# S = c_0 - c_1/Lambda^2 + c_2/Lambda^4 - ...
# The cosmological term is c_0 * f(0) WITH the appropriate sign convention.

# In GR: S_grav = integral (R - 2*Lambda_cc) * sqrt(g) / (16*pi*G)
# The spectral action gives: S ~ c_0 * f(0) (positive if f(0) > 0)
# This corresponds to a POSITIVE cosmological constant (de Sitter)
# if the spectral action is written as S = Tr(f(D^2/Lambda^2))

c_0 = N_total  # = 2640
c_1 = np.sum(all_evals)  # = 14880

print(f"\n  Spectral action coefficients (geometry only):")
print(f"    c_0 = {c_0} (cosmological term)")
print(f"    c_1 = {c_1:.0f} (Einstein-Hilbert)")

# Method 2: The EFFECTIVE cosmological constant from the spectral action
# The spectral action naturally gives:
# Lambda_cc_eff = (c_0 * f(0) * Lambda^3) / (c_1 * f'(0) * Lambda)
#               = c_0/c_1 * Lambda^2 * f(0)/f'(0)
# For Gaussian: f(0) = 1, f'(0) = -1
# Lambda_cc_eff = -c_0/c_1 * Lambda^2

ratio_c0_c1 = c_0 / c_1
print(f"\n  Ratio c_0/c_1 = {ratio_c0_c1:.6f}")
print(f"  Effective cosmological constant: Lambda_cc ~ -c_0/c_1 * Lambda^2")
print(f"  SIGN: NEGATIVE (for Gaussian cutoff f(x)=exp(-x))")

# Method 3: with the product geometry
c_0_prod = N_total * 9  # = 23760
c_1_prod = c_1 * 9 + N_total * 1858  # gravity + matter
n_exp_arr = np.array([0, 5, 3, 11, 16, 19, 26, 17, 11], dtype=float)
a_Y = np.sum(n_exp_arr**2)

print(f"\n  Product geometry:")
print(f"    c_0_prod = {c_0_prod}")
print(f"    c_1_prod = {c_1_prod:.0f}")
print(f"    c_0/c_1 ratio = {c_0_prod/c_1_prod:.6f}")

# Method 4: The physical cosmological constant
# In natural units of the 600-cell:
# Lambda_cc = (c_0/c_1) * (1/R^2) where R = radius of S^3
# At the Planck scale: R ~ l_Planck
# Lambda_cc ~ (c_0/c_1) / l_P^2 ~ (c_0/c_1) * M_P^2

# The observed Lambda_cc ~ 10^{-122} * M_P^2
# Our ratio c_0/c_1 ~ 0.177 which is O(1), NOT 10^{-122}
# This is the COSMOLOGICAL CONSTANT PROBLEM in any spectral action theory

print(f"\n  THE COSMOLOGICAL CONSTANT PROBLEM:")
print(f"    Spectral action gives: Lambda_cc ~ {ratio_c0_c1:.4f} * Lambda^2")
print(f"    This is O(1) in Planck units = way too large")
print(f"    Observed: Lambda_cc ~ 10^{{-122}} * M_P^2")
print(f"    PROBLEM: 122 orders of magnitude discrepancy")

# Method 5: eta invariant and cosmological constant
# The eta invariant eta(D) = regularized asymmetry of D spectrum
# For the simplicial D = d+d*, eta = 0 (symmetric spectrum)
# But in the PRODUCT geometry, D_total can have non-zero eta if D_F
# breaks the symmetry

# With D_F = diag(n_f), the D_total eigenvalues are +-sqrt(lambda_i + mu_j)
# The spectrum is STILL symmetric, so eta = 0

# The FERMIONIC contribution to the cosmological constant:
# Lambda_ferm = sum_f m_f^4 * log(m_f^2/Lambda^2) / (64*pi^2)
# This is the Coleman-Weinberg effective potential at the minimum

# With m_f = m_e * phi^n_f:
# sum(n_c * m_f^4) = m_e^4 * sum(n_c * phi^{4*n_f})
# Dominated by top: 3 * phi^{104} ~ 3 * 10^{21.7}

# The KEY question: does the geometry provide a CANCELLATION mechanism?
# In the 600-cell: the bosonic contribution c_0 is POSITIVE
# The fermionic contribution is NEGATIVE (fermion loop)
# Could they cancel?

# Fermion data for reference
fermion_data_local = {
    'e': {'n': 0, 'nc': 1}, 'mu': {'n': 11, 'nc': 1}, 'tau': {'n': 17, 'nc': 1},
    'u': {'n': 3, 'nc': 3}, 'c': {'n': 16, 'nc': 3}, 't': {'n': 26, 'nc': 3},
    'd': {'n': 5, 'nc': 3}, 's': {'n': 11, 'nc': 3}, 'b': {'n': 19, 'nc': 3},
}

b_Y_phys = sum(d['nc'] * d['n']**4 for d in fermion_data_local.values())

print(f"\n  BOSONIC vs FERMIONIC contributions:")
print(f"    Bosonic (spectral action): c_0 = +{c_0_prod} (positive)")
print(f"    Fermionic (one-loop): ~ -sum(n_c * n_f^4) * log-term")
print(f"    sum(n_c * n_f^4) = {b_Y_phys}")
print(f"    b_Y (with color) = {b_Y_phys}")

# Susy-like cancellation check
print(f"\n  Cancellation check:")
print(f"    Boson DOF: {c_0_prod} = {N_total}*{9}")
print(f"    Fermion DOF: 9 species * (color factors)")
n_ferm_dof = sum(d['nc'] * 2 for d in fermion_data_local.values())  # 2 for spin
print(f"    Effective fermion DOF (spin*color): {n_ferm_dof}")
print(f"    Boson/Fermion ratio: {c_0_prod/n_ferm_dof:.1f}")
print(f"    NO natural cancellation (not supersymmetric)")

# Connection to eta = -35
print(f"\n  Connection to spectral asymmetry:")
print(f"    eta(Delta_0) = -35 (from exp234)")
print(f"    h + eta = 30 + (-35) = -5 = -a_1")
print(f"    This gives the neutrino mass scale but NOT Lambda_cc")

print(f"\n  STATUS: OPEN")
print(f"  The cosmological constant problem persists in the 600-cell theory.")
print(f"  The spectral action gives Lambda_cc ~ O(1) * Lambda^2,")
print(f"  which is 122 orders of magnitude too large.")
print(f"  No natural cancellation mechanism found.")
print(f"  The sign is NEGATIVE (from the Gaussian spectral action).")

# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("SUMMARY: FOUR OPEN PROBLEMS")
print("=" * 72)

print(f"""
  PROBLEM 1: m_H/m_W = phi - 8*alpha
  ====================================
  RESULT: PARTIALLY RESOLVED
  - The formula is equivalent to m_H/m_W = phi - 2/(a_1*phi^4)
    to O(alpha^2) accuracy (diff = {abs(val_geom-val_alpha):.2e})
  - This is PURELY GEOMETRIC: only a_1=5 and phi appear
  - phi = base DSI scaling, 2/(a_1*phi^4) = quantum correction
  - STILL OPEN: derive this from the spectral action principle

  PROBLEM 2: Explicit D_F on McKay graph
  ========================================
  RESULT: CONSTRUCTED
  - D_F in mass basis: diagonal with n_f
  - D_F in interaction basis: CKM-rotated up-quark block
  - Off-diagonal mixing: {on_E8/(on_E8+off_E8)*100:.0f}% on E8 edges, {off_E8/(on_E8+off_E8)*100:.0f}% off
  - The CKM mixing is NOT confined to E8 edges
  - This suggests the node assignment needs refinement

  PROBLEM 3: 1+2+9 edge decomposition
  ======================================
  RESULT: FOUND!
  - For type-C vertices (96/120): neighbors split as
    1 type-A + 2 type-B + 9 type-C = 12
  - The CC* criterion: among 9 C-neighbors, the UNIQUE one sharing
    a type-A common neighbor but NO type-B common neighbor
  - This gives: 1(A) + 2(B) + 1(CC*) + 8(regular C) = 1 + 3 + 8
  - Global edge counts: A-C={n_AC}, B-C={n_BC}, CC*={n_CC_star}, CC_reg={n_CC_regular}
  - CC* neighbors per C-vertex: {dict(cc_star_count)}

  PROBLEM 4: Cosmological constant
  ==================================
  RESULT: OPEN (as expected)
  - Spectral action gives Lambda_cc ~ O(1) * Lambda^2 (too large)
  - Sign: NEGATIVE (from Gaussian cutoff ratio -c_0/c_1)
  - No natural boson-fermion cancellation
  - The 122-order-of-magnitude problem persists
  - This is a universal problem, not specific to the 600-cell theory
""")

print("=" * 72)
print("EXP-263 COMPLETE")
print("=" * 72)

"""
EXP-253: Dirac Operator from McKay Graph (Affine E8)
=====================================================
The McKay correspondence: 9 irreps of 2I <-> 9 nodes of affine E8.
Each node = one eigenvalue of 600-cell = one fermion species.

APPROACH: Build operators on the McKay graph (9x9 matrices)
and check if their eigenvalues give fermion mass exponents.

The affine E8 Dynkin diagram:

    [6]         (node labels = irrep dimensions = Coxeter labels)
     |
[1]-[2]-[3]-[4]-[5]-[6]-[4]-[2]

Nodes: 0=trivial(1), 1(2), 2(3), 3(4), 4(5), 5(6), 6(4'), 7(2'), branch=8(3')
The branch node connects to node 4(5).
Wait - let me use the standard labeling.

Standard affine E8 (extended Dynkin diagram):
Nodes labeled 0-8, with Coxeter labels (= irrep dims of 2I):

   3'(node 8)
    |
1 - 2 - 3 - 4 - 5 - 6 - 4' - 2'
(0) (1) (2) (3) (4) (5) (6)  (7)

Coxeter labels: [1, 2, 3, 4, 5, 6, 4, 2, 3]
These are the irrep dimensions of 2I.
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
m_e = 0.51099895  # MeV

# =====================================================================
# PART 1: Build the affine E8 Dynkin diagram
# =====================================================================
print("=" * 70)
print("EXP-253: McKAY GRAPH (AFFINE E8) DIRAC OPERATOR")
print("=" * 70)
print()

# Affine E8 adjacency matrix (9 nodes)
# Standard labeling: node 0 is the affine (extended) node
# Connections: 0-1, 1-2, 2-3, 3-4, 4-5, 5-6, 6-7, and branch 4-8
# Coxeter labels: a = [1, 2, 3, 4, 5, 6, 4, 2, 3]
# (These satisfy: sum(a_i) = 30 = h(E8) = Coxeter number)

# But wait - there are different conventions. Let me use the one where
# the McKay correspondence maps irreps to eigenvalues clearly.

# The irreps of 2I (binary icosahedral) have dimensions:
# 1, 2, 2', 3, 3', 4, 4', 5, 6
# Total: 1+2+2+3+3+4+4+5+6 = 30 = h(E8)

# The McKay graph: tensor product of each irrep with the
# 2-dimensional irrep (the "natural" rep of 2I in SU(2)):
# chi_natural * chi_i = sum_j A_ij chi_j

# For 2I, the natural rep is the 2-dim irrep.
# The McKay graph adjacency matrix A has A_ij = multiplicity of
# chi_j in chi_natural x chi_i.

# From literature, the affine E8 Dynkin diagram is:
# Nodes: 0(dim=1), 1(dim=2), 2(dim=3), 3(dim=4), 4(dim=5),
#         5(dim=6), 6(dim=4'), 7(dim=2'), 8(dim=3')
# Edges: 0-1, 1-2, 2-3, 3-4, 4-5, 5-6, 6-7, 4-8

n_nodes = 9
coxeter_labels = [1, 2, 3, 4, 5, 6, 4, 2, 3]  # irrep dimensions
node_labels = ['1', '2', '3', '4', '5', '6', "4'", "2'", "3'"]

# Adjacency matrix of affine E8
A_E8 = np.zeros((n_nodes, n_nodes))
edges_E8 = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (4,8)]
for i, j in edges_E8:
    A_E8[i,j] = 1
    A_E8[j,i] = 1

print("Affine E8 Dynkin diagram:")
print("  Nodes: %s" % node_labels)
print("  Coxeter labels: %s (sum = %d = h(E8))" % (coxeter_labels, sum(coxeter_labels)))
print("  Edges: %s" % edges_E8)
print()

# Eigenvalues of adjacency matrix of affine E8
eig_E8 = np.linalg.eigvalsh(A_E8)
eig_E8 = np.sort(eig_E8)[::-1]
print("  Eigenvalues of A(affine E8): %s" % np.round(eig_E8, 4))
print()

# Mapping: which eigenvalue of 600-cell corresponds to which node?
# From exp252:
# Node 0 (dim=1,  mult=1):   lambda = 12
# Node 1 (dim=2,  mult=4):   lambda = 6*phi  (Galois: 6-6*phi)
# Node 2 (dim=3,  mult=9):   lambda = 4*phi  (Galois: 4-4*phi)
# Node 3 (dim=4,  mult=16):  lambda = 3      (Galois: -3)
# Node 4 (dim=5,  mult=25):  lambda = 0
# Node 5 (dim=6,  mult=36):  lambda = -2
# Node 6 (dim=4', mult=16):  lambda = -3
# Node 7 (dim=2', mult=4):   lambda = 6-6*phi
# Node 8 (dim=3', mult=9):   lambda = 4-4*phi

# The Galois pairs: (node 1, node 7) and (node 2, node 8) and (node 3, node 6)
# These are exactly the primed/unprimed irrep pairs!

lambda_600 = [12, 6*PHI, 4*PHI, 3, 0, -2, -3, 6-6*PHI, 4-4*PHI]
lambda_ab = [(12,0), (0,6), (0,4), (3,0), (0,0), (-2,0), (-3,0), (6,-6), (4,-4)]

print("  McKay node -> 600-cell eigenvalue -> Z[phi] decomposition:")
for i in range(n_nodes):
    a, b = lambda_ab[i]
    n_mass = a1*a + b1*b
    print("    Node %d (dim=%d): lambda=%.4f = %d+%d*phi, n(mass)=%d"
          % (i, coxeter_labels[i], lambda_600[i], a, b, n_mass))

# =====================================================================
# PART 2: Operators on McKay graph
# =====================================================================
print()
print("=" * 70)
print("PART 2: Natural operators on the McKay graph")
print("=" * 70)
print()

# The Cartan matrix
C_E8 = 2 * np.eye(n_nodes) - A_E8
print("  Cartan matrix eigenvalues: %s" % np.round(np.sort(np.linalg.eigvalsh(C_E8)), 4))

# Coxeter-weighted adjacency
D_cox = np.diag(coxeter_labels)
print("  Coxeter label diagonal: %s" % coxeter_labels)

# Weighted Laplacian: D_cox^{-1/2} * A * D_cox^{1/2}
D_sqrt = np.diag(np.sqrt(coxeter_labels))
D_inv_sqrt = np.diag(1.0/np.sqrt(coxeter_labels))
A_weighted = D_inv_sqrt @ A_E8 @ D_sqrt
eig_Aw = np.sort(np.linalg.eigvalsh(A_weighted))[::-1]
print("  Weighted adjacency eigenvalues: %s" % np.round(eig_Aw, 4))

# =====================================================================
# PART 3: The KEY idea - eigenvalue-to-mass mapping via McKay
# =====================================================================
print()
print("=" * 70)
print("PART 3: Eigenvalue-to-fermion mapping")
print("=" * 70)
print()

# We need to map 9 eigenvalue nodes to 9 fermion species.
# The fermion species with their mass exponents:
fermions = [
    ('e',   0, 0, 0, 0),     # electron
    ('mu',  1, 1, 1, 11),    # muon
    ('tau', 2, 1, 2, 17),    # tau
    ('u',   0, 3,-2, 3),     # up
    ('c',   1, 2, 1, 16),    # charm
    ('t',   2, 4, 1, 26),    # top
    ('d',   0, 1, 0, 5),     # down
    ('s',   1, 1, 1, 11),    # strange
    ('b',   2,-1, 4, 19),    # bottom
]

print("  Fermions sorted by n: ", sorted([(f[0], f[4]) for f in fermions], key=lambda x: x[1]))
print()

# KEY INSIGHT: The E8 legs partition the nodes:
# Leg A (long leg): nodes 0-1-2-3-4-5-6-7 (length 8 = a_1+N_gen)
# Leg B (branch): node 8, connected to node 4
# In terms of physics:
# - The 3 legs of E8 have dimensions A=11, B=6, C=9 (from exp211)
# - CKM coefficients = leg dimensions

# From exp226b: generations localize on E8 legs
# psi_5 (gen 0) on Leg B, psi_3+psi_4 (gen 1,2) on Leg C
# all gens on Leg A via rho_8

# But let me try a different mapping based on MULTIPLICITIES.
# The adjacency eigenvalue multiplicities are dim^2:
# {1, 4, 9, 16, 25, 36, 16, 4, 9}
# The fermion mass exponents need multiplicities:
# n=0(1), n=3(1), n=5(1), n=11(2), n=16(1), n=17(1), n=19(1), n=26(1)

# n=11 has multiplicity 2 (muon + strange).
# Eigenvalue multiplicities: two nodes have same dim (1&7: dim 2, 2&8: dim 3, 3&6: dim 4)
# Galois pairs! So the n=11 degeneracy maps to a GALOIS PAIR.

# Which Galois pair? n=11 = 5*1 + 6*1, (a,b) = (1,1)
# As eigenvalue: which lambda has (a,b) close to (1,1)?
# None directly. But if we think of IRREPS not eigenvalues...

# The dim=3 and dim=3' irreps (nodes 2 and 8) form a Galois pair.
# dim=3 has lambda = 4*phi (mult 9)
# dim=3' has lambda = 4-4*phi (mult 9)
# Their n-values: n(0,4)=24 and n(4,-4)=-4. Not 11.

# What about dim=2 and dim=2' (nodes 1 and 7)?
# dim=2: lambda = 6*phi (mult 4), n = 36
# dim=2': lambda = 6-6*phi (mult 4), n = -6

# And dim=4 and dim=4' (nodes 3 and 6)?
# dim=4: lambda = 3 (mult 16), n = 15
# dim=4': lambda = -3 (mult 16), n = -15

# None give n=11. So the mapping is NOT through the mass functional.

# =====================================================================
# PART 4: What if the mass operator IS the Coxeter-weighted distance?
# =====================================================================
print()
print("=" * 70)
print("PART 4: Distance-based mass on McKay graph")
print("=" * 70)
print()

# On the affine E8, define the DISTANCE from a reference node.
# The natural reference: node 4 (dim=5 = a_1, the central node with highest label)
# or node 0 (dim=1, the trivial irrep = vacuum)

from collections import deque

def graph_distance(adj, source, n):
    """BFS distance from source."""
    dist = [-1] * n
    dist[source] = 0
    queue = deque([source])
    while queue:
        v = queue.popleft()
        for u in range(n):
            if adj[v, u] > 0 and dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist

for ref_node, ref_label in [(0, "trivial(1)"), (4, "center(5)"), (5, "max(6)")]:
    dists = graph_distance(A_E8, ref_node, n_nodes)
    print("  Distance from node %d (%s):" % (ref_node, ref_label))
    for i in range(n_nodes):
        print("    Node %d (dim=%d): dist=%d" % (i, coxeter_labels[i], dists[i]))
    print()

# =====================================================================
# PART 5: Coxeter element and its eigenvalues
# =====================================================================
print()
print("=" * 70)
print("PART 5: Coxeter element of E8")
print("=" * 70)
print()

# The Coxeter number h = 30 for E8.
# The exponents of E8 are: 1, 7, 11, 13, 17, 19, 23, 29
# These are the values m where exp(2*pi*i*m/h) are eigenvalues
# of the Coxeter element!

e8_exponents = [1, 7, 11, 13, 17, 19, 23, 29]
print("  E8 exponents: %s" % e8_exponents)
print("  Coxeter number h = 30")
print()

# REMARKABLE: The E8 exponents include 11, 17, 19 !!!
# These are EXACTLY the lepton/bottom mass exponents!
# n_mu = n_s = 11, n_tau = 17, n_b = 19

# What about the others?
# Fermion n-values: {0, 3, 5, 11, 16, 17, 19, 26}
# E8 exponents:     {1, 7, 11, 13, 17, 19, 23, 29}

# Differences:
print("  Fermion n-values:  {0, 3, 5, 11, 16, 17, 19, 26}")
print("  E8 exponents:      {1, 7, 11, 13, 17, 19, 23, 29}")
print()

# Check: are they related by a simple transformation?
n_fermion = sorted([0, 3, 5, 11, 16, 17, 19, 26])
e8_exp = sorted(e8_exponents)

print("  Difference (n_fermion - E8_exp): %s" % [n-e for n, e in zip(n_fermion, e8_exp)])
print("  Ratio: %s" % [round(n/e, 3) if e != 0 else 'inf' for n, e in zip(n_fermion, e8_exp)])
print()

# The E8 exponents satisfy: m + (h-m) = h = 30.
# So pairs: (1,29), (7,23), (11,19), (13,17)
# And our fermion pairs from Galois: n + n' should satisfy something...

print("  E8 exponent pairs (m, h-m): %s" %
      [(m, 30-m) for m in e8_exponents if m <= 15])
print("  Fermion exponent pairs (n, 30-n): %s" %
      [(n, 30-n) for n in n_fermion if n <= 15])
print()

# CHECK: 30 - n_fermion
complement = [30 - n for n in n_fermion]
print("  30 - n_fermion: %s" % sorted(complement))
print("  n_fermion:      %s" % n_fermion)
print("  Union:          %s" % sorted(set(n_fermion + complement)))
print()

# What about: n = E8_exponent - 1?
shifted = [e - 1 for e in e8_exp]
print("  E8_exp - 1:     %s" % shifted)
print("  n_fermion:      %s" % n_fermion)
print("  Difference:     %s" % [n - s for n, s in zip(n_fermion, shifted)])
print()

# What about using E8 exponents modulo something?
# Or: which 8 values from {0, 1, ..., 29} satisfy our constraints?
# E8 exponents are coprime to 30: gcd(m, 30) = 1
# Our n-values: gcd(0,30)=30, gcd(3,30)=3, gcd(5,30)=5, gcd(11,30)=1,
#               gcd(16,30)=2, gcd(17,30)=1, gcd(19,30)=1, gcd(26,30)=2

from math import gcd
print("  GCD with h=30:")
print("    E8 exponents: %s (ALL coprime to 30)" % [gcd(m, 30) for m in e8_exp])
print("    Fermion n:    %s (NOT all coprime!)" % [gcd(n, 30) for n in n_fermion])
print()

# =====================================================================
# PART 6: The DEEP connection - E8 exponents contain {11, 17, 19}
# =====================================================================
print()
print("=" * 70)
print("PART 6: E8 exponents and fermion masses - the overlap")
print("=" * 70)
print()

overlap = set(n_fermion) & set(e8_exponents)
print("  Overlap of fermion n-values and E8 exponents: %s" % sorted(overlap))
print("  These are: n=11 (muon/strange), n=17 (tau), n=19 (bottom)")
print()

# 3 out of 8 fermion mass exponents ARE E8 exponents!
# And these 3 are exactly the LEPTON exponents (11, 17) plus bottom (19)!
# The other E8 exponents are {1, 7, 13, 23, 29}.

# What if the QUARK exponents are obtained from E8 exponents by
# some transformation?
# Missing from fermions: {1, 7, 13, 23, 29}
# Missing from E8: {0, 3, 5, 16, 26}

# Check: are the quark exponents related to E8 exponents?
# n_u = 3 = ? , n_c = 16 = ?, n_t = 26 = ?, n_d = 5 = ?
#
# 3 = 30 - 27 ... no. 3 = 1 + 2?
# 16 = 30 - 14 ... 14 = h/2 - 1. Hmm.
# 26 = 30 - 4. And 4 = dim(irrep)...
# 5 = 30 - 25. And 25 = 5^2 = a_1^2...

# Or: 30 - n_quark = {27, 14, 4, 25}. Not illuminating.

# What about: n_up + n_down = {3+5, 16+11, 26+19} = {8, 27, 45}
# We already know n_u + n_d = 8 = rank(E8)!
# 27 = 3^3 or 27 = h - 3
# 45 = a_1 * N_eig

print("  Sum rules: n_up + n_down per generation:")
for g in range(3):
    ups = [f for f in fermions if f[0] in ['u','c','t'] and f[1] == g]
    downs = [f for f in fermions if f[0] in ['d','s','b'] and f[1] == g]
    if ups and downs:
        nu, nd = ups[0][4], downs[0][4]
        print("    Gen %d: n_up=%d + n_down=%d = %d" % (g, nu, nd, nu+nd))

print()

# =====================================================================
# PART 7: The Coxeter transformation on mass exponents
# =====================================================================
print()
print("=" * 70)
print("PART 7: Coxeter transformation: n -> h - n = 30 - n")
print("=" * 70)
print()

print("  If we apply the Coxeter involution (m -> h-m) to mass exponents:")
for name, g, a, b, n in fermions:
    n_cox = 30 - n
    is_e8 = n_cox in e8_exponents
    is_fermion = n_cox in [f[4] for f in fermions]
    print("    %s: n=%2d -> 30-n=%2d  %s%s"
          % (name, n, n_cox,
             "(E8 exponent!)" if is_e8 else "",
             "(= fermion n!)" if is_fermion else ""))

print()

# Remarkable findings:
# e(0) -> 30: not special
# mu/s(11) -> 19 = n_bottom! And 19 IS an E8 exponent!
# tau(17) -> 13 = E8 exponent (not a fermion mass)
# u(3) -> 27: not an E8 exponent
# c(16) -> 14: not an E8 exponent
# t(26) -> 4: not an E8 exponent
# d(5) -> 25: not an E8 exponent
# b(19) -> 11 = n_muon/strange! And 11 IS an E8 exponent!

print("  KEY: Coxeter involution swaps muon/strange (11) <-> bottom (19)!")
print("  Both 11 and 19 are E8 exponents AND fermion mass exponents!")
print("  They form a Coxeter pair: 11 + 19 = 30 = h(E8)")
print()
print("  Also: tau (17) maps to E8 exponent 13 (not a fermion mass)")
print("  And: 17 + 13 = 30 = h(E8)")
print()

# What are ALL the E8 exponent pairs?
# (1, 29), (7, 23), (11, 19), (13, 17)
# Of these, (11, 19) and (13, 17) contain fermion mass exponents.
# (1, 29) and (7, 23) do not.

# Fermion n-values: {0, 3, 5, 11, 16, 17, 19, 26}
# E8 exponents:     {1, 7, 11, 13, 17, 19, 23, 29}
# Intersection:     {11, 17, 19}

# The 3 fermions whose n IS an E8 exponent: muon(11), tau(17), bottom(19)
# Pattern: these are the HEAVY leptons and the HEAVY down-type quark!
# Or: these are the fermions on the a=1 line with b = {1, 2, 4}
# (muon: b=1, tau: b=2, bottom: b=4) -- b values are powers of 2!
# Wait no, bottom has b=4 but a=-1, not a=1.

# =====================================================================
# PART 8: Generate mass exponents from E8 exponents?
# =====================================================================
print()
print("=" * 70)
print("PART 8: Can fermion n-values be generated from E8 structure?")
print("=" * 70)
print()

# E8 exponents: {1, 7, 11, 13, 17, 19, 23, 29}
# They come in pairs summing to h=30.
# They are exactly the numbers coprime to 30 in {1,...,29}.

# Our fermion n-values: {0, 3, 5, 11, 16, 17, 19, 26}
# Observation: 0 = ground state (electron)
# The other 7: {3, 5, 11, 16, 17, 19, 26}

# Try: n = E8_exponent mod something?
# Or: n obtained from E8 exponents by some algebraic operation?

# What if: each fermion n is the CLOSEST E8 exponent + correction?
print("  Fermion n -> nearest E8 exponent -> difference:")
for name, g, a, b, n in fermions:
    nearest = min(e8_exponents, key=lambda e: abs(e - n))
    diff = n - nearest
    print("    %s: n=%2d, nearest E8 exp=%2d, diff=%+d" % (name, n, nearest, diff))

print()
# Results: e(0,nearest=1,diff=-1), mu(11,11,0), tau(17,17,0),
#          u(3,1,+2), c(16,17,-1), t(26,29,-3),
#          d(5,7,-2), s(11,11,0), b(19,19,0)

# Pattern: leptons/strange/bottom have diff=0 (ARE E8 exponents or ground state)
#          quarks have diff = {+2, -1, -3, -2} (small corrections)

# Can we express the quark corrections in terms of derived quantities?
# u: diff = +2 = F(3) = Fibonacci
# c: diff = -1
# t: diff = -3 = -N_gen
# d: diff = -2 = -F(3)
# These don't form a clean pattern.

# =====================================================================
# PART 9: The spectral zeta function approach
# =====================================================================
print()
print("=" * 70)
print("PART 9: Spectral zeta function and mass formula")
print("=" * 70)
print()

# In Connes' NCG, the spectral action is Tr(f(D/Lambda)).
# The mass information is encoded in the RESIDUES of the spectral
# zeta function: zeta_D(s) = Tr(|D|^{-s}).

# For our 600-cell, the zeta function of the adjacency matrix is:
# zeta_A(s) = sum_i mult_i * |lambda_i|^{-s}

# What if the MASS of each fermion species is related to the
# contribution of its eigenspace to the zeta function?

print("  Spectral zeta function of 600-cell adjacency matrix:")
print("  zeta_A(s) = sum_i mult_i * |lambda_i|^{-s}")
print()
print("  At s=1 (first moment):")

zeta_1 = 0
for i in range(n_nodes):
    lam = abs(lambda_600[i])
    mult = coxeter_labels[i]**2
    if lam > 0.01:
        contrib = mult / lam
        zeta_1 += contrib
        print("    Node %d: %d^2 / |%.4f| = %.4f" % (i, coxeter_labels[i], lambda_600[i], contrib))

print("  Total zeta_A(1) = %.4f" % zeta_1)
print()

# What if mass ~ multipliciy * eigenvalue processed through phi?
print("  Alternative: mass exponent = f(dim_irrep, eigenvalue)?")
print()
for i in range(n_nodes):
    dim_i = coxeter_labels[i]
    lam_i = lambda_600[i]
    a_i, b_i = lambda_ab[i]

    # Various derived combinations
    n1 = dim_i * (dim_i + 1) / 2  # triangular number
    n2 = a1 * dim_i - b1  # linear in dim
    n3 = dim_i**2 - 1  # dim^2 - 1

    print("    Node %d (dim=%d, lam=%.3f): T(dim)=%.0f, 5*dim-6=%.0f, dim^2-1=%.0f"
          % (i, dim_i, lam_i, n1, n2, n3))

# =====================================================================
# SUMMARY
# =====================================================================
print()
print("=" * 70)
print("SUMMARY OF DISCOVERIES")
print("=" * 70)
print()
print("1. E8 EXPONENTS contain 3 fermion mass exponents: {11, 17, 19}")
print("   These are muon/strange, tau, and bottom -- the HEAVY ones!")
print()
print("2. COXETER INVOLUTION (n -> 30-n) swaps:")
print("   muon/strange (11) <-> bottom (19) -- both E8 exponents!")
print("   tau (17) <-> 13 (E8 exponent, not a fermion)")
print()
print("3. The 4 Coxeter PAIRS of E8: (1,29), (7,23), (11,19), (13,17)")
print("   Pairs (11,19) and (13,17) contain ALL our E8-overlapping n-values")
print()
print("4. Quark n-values differ from nearest E8 exponent by SMALL amounts")
print("   u: +2, c: -1, t: -3, d: -2 -- could be generation corrections?")
print()
print("5. The McKay graph operators DON'T directly give mass eigenvalues,")
print("   but the E8 EXPONENTS partially overlap with fermion n-values.")
print()
print("OPEN DIRECTIONS:")
print("  - Why are lepton n-values E8 exponents but quark n-values are not?")
print("  - Can quark corrections (+2,-1,-3,-2) be derived from Casimir?")
print("  - The spectral triple (A, H, D) with D encoding masses")
print("  - Connection between Coxeter element and generation structure")

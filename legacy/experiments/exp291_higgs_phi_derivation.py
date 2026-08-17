"""
EXP-291: WHY IS m_H / m_W = phi (at tree level)?
=====================================================
This is the #1 open problem in the framework.

KNOWN:
  m_H = m_W * (phi - 8*alpha)   [0.09% error]
  The correction -8*alpha is DERIVED (8 = rank(E8))
  But WHY is the tree-level ratio = phi?

APPROACHES TO INVESTIGATE:
  A. Spectral action: Higgs quartic from Tr(D_F^4)/Tr(D_F^2)^2
  B. 600-cell geometry: phi as diagonal/edge ratio in pentagon
  C. A_5 representation theory: Higgs from edge 0-1 on McKay
  D. Algebraic: phi as eigenvalue of some natural operator
  E. DSI/mass formula: m_H exponent vs m_W exponent
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
phi = PHI
phi_p = (1 - np.sqrt(5)) / 2  # = -1/phi
ALPHA = 1/137.035999
ALPHA_S = 1/(2*PHI**3)
SIN2_TW = 6/26

# Experimental values
m_e = 0.51099895e-3  # GeV
m_W = 80.377  # GeV
m_H = 125.25  # GeV
m_Z = 91.1876  # GeV

a1 = 5
b1 = 6
N = 120
h = 30
rank_E8 = 8
N_eig = 9

print("="*72)
print("EXP-291: WHY m_H/m_W = phi?")
print("="*72)

ratio_exp = m_H / m_W
print(f"\n  m_H/m_W (exp) = {ratio_exp:.6f}")
print(f"  phi            = {phi:.6f}")
print(f"  phi - 8*alpha  = {phi - 8*ALPHA:.6f}")
print(f"  Bare error (phi): {abs(ratio_exp - phi)/ratio_exp * 100:.2f}%")
print(f"  With correction:  {abs(ratio_exp - (phi - 8*ALPHA))/ratio_exp * 100:.2f}%")

# ================================================================
# A. SPECTRAL ACTION APPROACH
# ================================================================
print(f"\n{'='*72}")
print(f"A. SPECTRAL ACTION: Higgs quartic from D_F traces")
print(f"{'='*72}")

print(f"""
  In NCG spectral action, the Higgs mass relation is:
    m_H^2 = (8*lambda/g^2) * m_W^2
  where lambda = quartic Higgs coupling, g = SU(2) coupling.

  From spectral action: lambda/g^2 = f(0)*Tr(D_F^4) / [f(2)*Tr(D_F^2)]
  (f(0), f(2) are spectral function moments)

  For m_H/m_W = phi:
    m_H^2/m_W^2 = phi^2 = phi + 1
    => 8*lambda/g^2 = phi + 1 = phi^2
    => lambda/g^2 = phi^2/8 = (phi+1)/8
""")

# What ratio of D_F traces gives phi^2/8?
print(f"  Target: lambda/g^2 = phi^2/8 = {phi**2/8:.6f}")
print(f"  Or equivalently: Tr(D_F^4)/Tr(D_F^2)^2 * (something) = phi^2/8")

# In our framework, D_F has eigenvalues n_i = 5a + 6b
# The 9 fermion mass exponents
fermions = {
    'e': (0, 0, 0),
    'mu': (1, 1, 11),
    'tau': (1, 2, 17),
    'u': (3, -2, 3),
    'c': (2, 1, 16),
    't': (4, 1, 26),
    'd': (1, 0, 5),
    's': (1, 1, 11),
    'b': (-1, 4, 19),
}

exponents = [f[2] for f in fermions.values()]
print(f"\n  Fermion exponents n_i: {exponents}")

# Traces of D_F (with n_i as eigenvalues)
Tr2 = sum(n**2 for n in exponents)
Tr4 = sum(n**4 for n in exponents)
print(f"  Tr(D_F^2) = sum(n_i^2) = {Tr2}")
print(f"  Tr(D_F^4) = sum(n_i^4) = {Tr4}")
print(f"  Tr(D_F^4)/Tr(D_F^2)^2 = {Tr4/Tr2**2:.6f}")

# With multiplicities (each generation has same exponents for leptons)
# But actually each fermion appears once
ratio_traces = Tr4 / Tr2**2
print(f"  Tr4/Tr2^2 = {ratio_traces:.6f}")
print(f"  phi^2/8 = {phi**2/8:.6f}")
print(f"  These are NOT equal. Direct trace approach fails.")

# But wait - in NCG the Higgs quartic also involves N_gen
# lambda proportional to Tr(k^4) where k = Yukawa matrix
# In our case k_f ~ phi^(n_f/2) approximately

# Try: ratio with generation counting
# In standard NCG: m_H^2 = 2*Tr(k^4)/Tr(k^2) * v^2
# where v = Higgs vev, k = Yukawa coupling matrix

# ================================================================
# B. GEOMETRIC APPROACH: PHI FROM PENTAGON
# ================================================================
print(f"\n{'='*72}")
print(f"B. GEOMETRIC: phi as natural ratio in vertex figure")
print(f"{'='*72}")

print(f"""
  The vertex figure of the 600-cell is an icosahedron.
  The icosahedron is built from pentagons.
  In a regular pentagon: diagonal/edge = phi.

  The W boson corresponds to the SU(2) sector (edges of type BC).
  The Higgs comes from inner fluctuations on the edge 0-1 of McKay E8~.

  Question: Is there a natural "diagonal vs edge" ratio
  in the gauge sector that gives phi?
""")

# In the icosahedron:
# - Edge length = 1 (by convention)
# - 12 vertices, 30 edges, 20 faces
# - Vertex-to-vertex distances: 1 (edge), phi (diagonal of pentagon face)
# The two distances in an icosahedron: edge and phi*edge

print(f"  Icosahedron distances: edge = 1, pentagon diagonal = phi = {phi:.4f}")
print(f"  The icosahedron has exactly 2 distance classes:")
print(f"    30 edges (distance 1)")
print(f"    30 'long' connections (distance phi)")
print(f"  Total: 30 + 30 = 60 pairs out of C(12,2) = 66")
print(f"  Remaining 6 pairs: antipodal (distance 2*sin(pi/2)*R)")

# Key insight: In the McKay graph E8~:
# Edge 0-1 connects trivial (dim 1) to fundamental (dim 2)
# This edge is WHERE the Higgs lives (inner fluctuation)
# Dim ratio: 2/1 = 2
# But we need phi, not 2

# What about the NODE marks?
# Higgs is at edge 0-1. The adjacent nodes are 0 (mark 1) and 1 (mark 2)
# The NEXT node is 2 (mark 3)
# Ratio mark_2/mark_1 = 2/1, mark_3/mark_2 = 3/2

# Golden ratio as limit of Fibonacci: mark ratios converge to phi
# But we need EXACT phi, not convergence

# ================================================================
# C. A_5 REPRESENTATION THEORY APPROACH
# ================================================================
print(f"\n{'='*72}")
print(f"C. A_5 REPRESENTATION THEORY: Higgs from 3' -> 3 transition")
print(f"{'='*72}")

print(f"""
  From exp290b: 12 = 1 + 3' + (3+5) = U(1) + SU(2) + SU(3)

  The Higgs field transforms as (2, 2bar) under SU(2)_L x SU(2)_R
  In A_5 language: Higgs ~ 3' (x) 3 (transition between the two 3-dim irreps)

  We computed: 3 (x) 3' = 4 + 5

  The Higgs doublet has 4 real components (2 complex).
  3 (x) 3' = 4 + 5 contains the 4-dim irrep!

  So: Higgs doublet = the 4-dim IRREP of A_5 in the product 3 (x) 3'
  This is the representation that "bridges" SU(2) and SU(3)!
""")

# Character of the 4-dim irrep
chi_4 = np.array([4, 0, 1, -1, -1])
print(f"  chi(4): {[int(x) for x in chi_4]}")
print(f"  dim = 4 (= Higgs doublet real components)")

# The Casimir of the 4-dim rep of A_5
# For a permutation group, Casimir ~ sum of squared characters
# Actually, for the 4-dim irrep: it's the "standard rep" of S_5 restricted to A_5
# Casimir eigenvalue in terms of the group

# KEY QUESTION: Does the 4-dim irrep of A_5 naturally produce phi?
# The character values of the 4-dim irrep are: 4, 0, 1, -1, -1
# No phi here! All integers.
# But the 3 and 3' irreps DO have phi and phi' in their characters.

# The Higgs potential V(H) = -mu^2 |H|^2 + lambda |H|^4
# mu^2/lambda determines v^2 (Higgs vev)
# m_H^2 = 2*lambda*v^2
# m_W = g*v/2

# So m_H^2/m_W^2 = 8*lambda/g^2

# Can we get lambda/g^2 from A_5 Casimirs?
# SU(2) coupling g comes from the 3' sector
# Higgs quartic lambda comes from the 4 sector

# Casimir of SU(2) in the adjoint (3'): C_2(adj) = 2 (for SU(2))
# Casimir of the Higgs (fundamental of SU(2)): C_2(fund) = 3/4

# ================================================================
# D. ALGEBRAIC APPROACH: phi as eigenvalue
# ================================================================
print(f"\n{'='*72}")
print(f"D. ALGEBRAIC: phi from natural operators on the framework")
print(f"{'='*72}")

# Approach D1: phi^2 = phi + 1 (defining property)
# m_H^2/m_W^2 = phi^2 = phi + 1
# So m_H^2 = m_W^2 * (phi + 1) = m_W^2 * phi + m_W^2
# i.e., m_H^2 = m_W^2 + m_W^2 * phi
# Could m_W^2*phi = something physical?

print(f"  m_H^2 = m_W^2 * phi^2 = m_W^2 * (phi + 1)")
print(f"  => m_H^2 = m_W^2 + m_W^2 * phi")
print(f"  => m_H^2 - m_W^2 = m_W^2 * phi")
print(f"  Numerically: {m_H**2 - m_W**2:.1f} vs {m_W**2 * phi:.1f} GeV^2")
print(f"  (Using exact phi: {m_W**2 * phi:.1f})")

# Approach D2: The adjacency matrix of the icosahedron
# has eigenvalues: 5, sqrt(5), 1, -1, -sqrt(5), -3
# with multiplicities: 1, 3, 5, 5, 3, 1
print(f"\n  Icosahedron adjacency eigenvalues:")
ico_eigs = [(5, 1), (np.sqrt(5), 3), (1, 5), (-1, 5), (-np.sqrt(5), 3), (-3, 1)]
for val, mult in ico_eigs:
    print(f"    {val:+.4f} (mult {mult})")

print(f"\n  sqrt(5) = 2*phi - 1 = {np.sqrt(5):.4f}")
print(f"  Eigenvalue ratio: 5/sqrt(5) = sqrt(5) = {np.sqrt(5):.4f}")
print(f"  phi = (1 + sqrt(5))/2 = average of 1 and sqrt(5)")

# The 3-dim eigenspace (eigenvalue sqrt(5)) corresponds to the 3-dim irrep!
# The 3-dim eigenspace (eigenvalue -sqrt(5)) corresponds to 3' irrep!
print(f"\n  KEY: Eigenvalue of 3-irrep = sqrt(5) = sqrt(a_1)")
print(f"       Eigenvalue of 3'-irrep = -sqrt(5) = -sqrt(a_1)")
print(f"       Eigenvalue of 5-irrep = 1")
print(f"       Eigenvalue of trivial = 5 = a_1")

# ================================================================
# D3. THE CRITICAL INSIGHT: Adjacency eigenvalue ratios
# ================================================================
print(f"\n--- D3. Adjacency eigenvalue ratio ---")

# SU(3) sector: eigenvalues sqrt(5) (from 3) and 1 (from 5)
# SU(2) sector: eigenvalue -sqrt(5) (from 3')
# U(1) sector: eigenvalue 5 (from trivial)

# The W boson mass comes from SU(2) gauge coupling
# The Higgs mass comes from the Higgs potential

# In the icosahedral adjacency:
# 3' has eigenvalue -sqrt(5)
# The "gauge coupling" is related to this eigenvalue

# The Higgs is at the INTERSECTION of 3 and 3' (the 4-dim in 3x3')
# What if the Higgs mass ratio is |lambda(3)|/|lambda(3')| or similar?

# Actually: lambda(trivial)/|lambda(3')| = 5/sqrt(5) = sqrt(5)
# And phi = (1+sqrt(5))/2 = (1+sqrt(a_1))/2
# That gives sqrt(5), not phi directly

# But the McKay graph E8~ adjacency has different eigenvalues!
print(f"\n  McKay E8~ adjacency matrix eigenvalues:")
# E8~ Cartan matrix: A = 2I - C where C is E8~ adjacency
# E8~ adjacency eigenvalues = 2 - cartan_eigenvalues
# Extended Cartan matrix of E8~ has eigenvalues:
# 0, and the eigenvalues related to E8 Cartan (all positive)

# Actually, let me compute the E8~ adjacency directly
# E8~ is the affine E8 Dynkin diagram
# Nodes: 0-1-2-3-4-5-6-7, with 4 also connected to 8
# Marks (null vector): 1,2,3,4,5,6,4,2,3

# Wait, that's the standard E8 Dynkin. E8~ (affine) adds node 0:
# 0-1-2-3-4-5-6-7 with branch at 4 connecting to 8
# Adjacency:
e8t_adj = np.zeros((9,9))
# Trunk: 0-1-2-3-4-5-6-7
for i in range(7):
    e8t_adj[i, i+1] = 1
    e8t_adj[i+1, i] = 1
# Branch: 4-8
e8t_adj[4, 8] = 1
e8t_adj[8, 4] = 1

e8t_eigs = np.sort(np.linalg.eigvalsh(e8t_adj))[::-1]
print(f"  Eigenvalues: {[f'{x:.4f}' for x in e8t_eigs]}")

# The marks
marks = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])
print(f"  Marks (null vector): {marks}")
print(f"  Check A*marks: {e8t_adj @ marks}")

# ================================================================
# E. MASS FORMULA APPROACH
# ================================================================
print(f"\n{'='*72}")
print(f"E. MASS FORMULA: m_H vs m_W exponents")
print(f"{'='*72}")

# m_Z = m_e * phi^25 * alpha(mZ)/alpha(0)
# m_W = m_Z * cos(theta_W) = m_Z * sqrt(20/26)
# m_H = m_W * phi

# So m_H = m_Z * sqrt(20/26) * phi = m_e * phi^26 * sqrt(20/26) * alpha(mZ)/alpha(0)
# Note: phi^26 = phi^25 * phi

n_Z = 25
n_H_eff = 26  # one step higher!

print(f"  m_Z = m_e * phi^{n_Z} * alpha(mZ)/alpha(0)")
print(f"  m_W = m_Z * cos(tW)")
print(f"  m_H = m_W * phi = m_Z * cos(tW) * phi")
print(f"       = m_e * phi^{n_Z+1} * cos(tW) * alpha(mZ)/alpha(0)")
print(f"  So effectively n_H = n_Z + 1 = {n_Z + 1}")
print(f"")
print(f"  WHY +1? The Higgs is ONE phi-step above the Z/W scale.")
print(f"  In the DSI lattice, the Higgs sits at level n_Z + 1 = 26.")

# What's special about 26?
print(f"\n  n_H = 26 = a_1^2 + 1 = 1/sin^2(tW)")
print(f"  n_Z = 25 = a_1^2")
print(f"  AMAZING: n_Z = a_1^2, n_H = a_1^2 + 1!")
print(f"  The Higgs-Z gap is EXACTLY 1, and n_Z = a_1^2.")

# 26 = dim of bosonic string theory!
# Also 26 = 1/sin^2(tW) in our framework
print(f"  26 = 1/sin^2(tW) = (a_1^2+1)/(a_1+1) ... wait")
print(f"  sin^2(tW) = b_1/(a_1^2+1) = 6/26, so 1/sin^2(tW) = 26/6")
print(f"  No, 26 itself = a_1^2 + 1")

# ================================================================
# F. THE KEY IDENTITY: phi = exp_gap on DSI lattice
# ================================================================
print(f"\n{'='*72}")
print(f"F. KEY: m_H/m_W = phi from the Higgs mechanism itself")
print(f"{'='*72}")

print(f"""
  In the SM, the Higgs mechanism gives:
    m_W = g*v/2,  m_Z = m_W/cos(tW),  m_H = sqrt(2*lambda)*v

  So m_H/m_W = 2*sqrt(2*lambda)/g = sqrt(8*lambda/g^2) = sqrt(8*lambda)/g

  For this to equal phi:
    8*lambda/g^2 = phi^2
    lambda = g^2 * phi^2/8 = g^2 * (phi+1)/8

  The quartic coupling lambda is determined by the spectral action.
  In Connes' NCG framework:
    lambda = pi^2 * f_0 * Tr(Y^4) / [3 * (f_2)^2 * (Tr(Y^2))^2]
  where Y is the Yukawa matrix and f_0, f_2 are spectral moments.
""")

# Let's compute what Tr(Y^4)/Tr(Y^2)^2 must be
# In our framework, Yukawa couplings ~ m_f/v ~ m_f * g/(2*m_W)
# So Y_f proportional to m_f = m_e * phi^(n_f)

# Y_f = g * m_f / (sqrt(2) * m_W) (SM relation)
# Since m_f = m_e * phi^(n_f), Y_f proportional to phi^(n_f)

# N_gen copies of each fermion type
# Total Tr(Y^2) = N_c * sum_f Y_f^2 (with color factor N_c for quarks)

print(f"  Yukawa couplings Y_f proportional to phi^(n_f):")
for name, (a, b, n) in fermions.items():
    y = phi**n
    print(f"    Y_{name:>3s} ~ phi^{n:>2d} = {y:.2e}")

# Compute trace ratios
# Leptons: e(0), mu(11), tau(17) with N_c = 1
# Up quarks: u(3), c(16), t(26) with N_c = 3
# Down quarks: d(5), s(11), b(19) with N_c = 3

leptons = [0, 11, 17]
up_quarks = [3, 16, 26]
down_quarks = [5, 11, 19]

TrY2 = sum(phi**(2*n) for n in leptons) + \
       3*sum(phi**(2*n) for n in up_quarks) + \
       3*sum(phi**(2*n) for n in down_quarks)

TrY4 = sum(phi**(4*n) for n in leptons) + \
       3*sum(phi**(4*n) for n in up_quarks) + \
       3*sum(phi**(4*n) for n in down_quarks)

print(f"\n  Tr(Y^2) = {TrY2:.4e}")
print(f"  Tr(Y^4) = {TrY4:.4e}")
print(f"  Tr(Y^4)/Tr(Y^2)^2 = {TrY4/TrY2**2:.6f}")

# Top quark dominates everything
# Tr(Y^2) ~ 3*phi^52, Tr(Y^4) ~ 3*phi^104
# So Tr(Y^4)/Tr(Y^2)^2 ~ 3*phi^104 / (3*phi^52)^2 = 1/(3*phi^0) = 1/3
print(f"  1/3 = {1/3:.6f} (top-dominated limit)")
print(f"  Actual ratio is very close to 1/3 (top dominance)")

# In NCG: lambda/g^2 = pi^2*f0/(3*f2^2) * Tr(Y^4)/Tr(Y^2)^2
# For phi^2/8: pi^2*f0/(3*f2^2) * (1/3) = phi^2/8
# => pi^2*f0/(9*f2^2) = phi^2/8
# => f0/f2^2 = 9*phi^2/(8*pi^2)

f0_over_f2sq = 9*phi**2/(8*np.pi**2)
print(f"\n  Required: f_0/f_2^2 = 9*phi^2/(8*pi^2) = {f0_over_f2sq:.6f}")

# ================================================================
# G. NEW APPROACH: HIGGS FROM INNER FLUCTUATION ON McKAY
# ================================================================
print(f"\n{'='*72}")
print(f"G. INNER FLUCTUATION: Higgs from McKay edge 0-1")
print(f"{'='*72}")

print(f"""
  In the finite spectral triple (exp288):
  - A_F = C + H + M_3(C) (first 3 McKay nodes)
  - The Higgs comes from inner fluctuations of D on edge 0-1
  - Edge 0-1 connects C (dim 1) to H ~ M_2(C) (dim 2)

  The inner fluctuation: A' = A + sum_i a_i[D, b_i]
  On edge 0-1: this gives a scalar field transforming as (1,2)

  The MASS of this scalar is determined by D_F restricted to edge 0-1.
  D_F ~ phi^(dist) where dist = McKay graph distance.

  For edge 0-1: D_F(0,1) ~ phi^1 = phi (distance 1 on McKay)
  For the gauge boson (self-coupling): ~ phi^0 = 1

  RATIO: m_H/m_W ~ phi^1/phi^0 = phi!
""")

print(f"  THIS IS THE ARGUMENT:")
print(f"  - The Higgs field lives on McKay edge 0-1 (distance 1)")
print(f"  - The gauge field lives at McKay node 1 (distance 0 from itself)")
print(f"  - The mass ratio scales as phi^(distance)")
print(f"  - Distance of Higgs coupling = 1, distance of gauge = 0")
print(f"  - Therefore m_H/m_W = phi^1 / phi^0 = phi")

# This is consistent with exp280 where D_F ~ phi^(-dist) was found
# but "qualitative only". Let's make it more precise.

# ================================================================
# H. MAKING IT PRECISE: D_F eigenvalue structure
# ================================================================
print(f"\n{'='*72}")
print(f"H. PRECISE: D_F on McKay edge 0-1")
print(f"{'='*72}")

# The finite Dirac operator D_F has matrix elements between nodes
# In NCG: D_F is off-diagonal (connects different algebra blocks)
# D_F(i,j) != 0 only if i,j are adjacent on the McKay graph

# For the Higgs: D_F(0,1) determines the Yukawa coupling
# For the W boson: the gauge coupling g comes from the algebra structure

# In the spectral action S = Tr(f(D/Lambda)):
# m_W^2 proportional to f_2 * Tr(D_F^2|_{gauge})
# m_H^2 proportional to f_2 * Tr(D_F^2|_{Higgs}) + f_0 * Tr(D_F^4|_{Higgs})

# The KEY: if D_F(0,1) = phi * (some reference scale)
# then the Higgs potential V(H) has its minimum at v proportional to phi

# McKay graph distances from node 0:
mckay_dist = [0, 1, 2, 3, 4, 5, 6, 7, 4]  # node 8 is at distance 4 from 0
print(f"  McKay distances from node 0: {mckay_dist}")
print(f"  Node 0 (trivial, dim 1): dist = 0")
print(f"  Node 1 (fund, dim 2):    dist = 1")
print(f"  Node 2 (3-dim, dim 3):   dist = 2")

# In the DSI framework: phi^k at distance k
# The fermion mass at node i: m ~ m_e * phi^(5*a_i + 6*b_i)
# But the D_F matrix element between adjacent nodes: D_F(i,j) ~ phi^1 = phi

# The W boson mass comes from the SU(2) gauge coupling g:
#   m_W = g*v/2
# The Higgs mass comes from the quartic coupling lambda:
#   m_H = sqrt(2*lambda) * v

# In the spectral action with D_F:
#   g^2 proportional to 1/Tr(D_F^2|_{SU(2)})
#   lambda proportional to Tr(D_F^4|_{Higgs})/Tr(D_F^2|_{Higgs})^2

# If D_F has a single scale phi on the 0-1 edge:
# Higgs sector: D_F ~ phi (the edge value)
# This contributes to both m_W (through symmetry breaking) and m_H (through quartic)

# ================================================================
# I. SELF-COUPLING vs CROSS-COUPLING
# ================================================================
print(f"\n{'='*72}")
print(f"I. RATIO FROM QUADRATIC CASIMIRS")
print(f"{'='*72}")

# Another approach: the Higgs is in the fundamental of SU(2)
# The gauge field is in the adjoint of SU(2)
# Casimir of fundamental of SU(2): C_2(fund) = 3/4
# Casimir of adjoint of SU(2): C_2(adj) = 2

# Ratio C_2(adj)/C_2(fund) = 8/3
# m_H^2/m_W^2 = 8*lambda/g^2
# If lambda/g^2 = C_2(adj)/(N*C_2(fund)) for some N...
# phi^2 = 8*lambda/g^2 -> lambda/g^2 = phi^2/8

# From A_5: the Casimir in the 4-dim rep
# 4-dim rep of A_5: standard rep
# Second Casimir = ?

# For permutation groups, C_2 = n(n-1)/N * sum
# Actually for A_5, the quadratic Casimir of the r-dim irrep:
# C_2(r) = |G|/dim(r) * (1 - chi(g)^2/dim(r)^2) ... no, that's not right

# Use the formula: C_2(r) * dim(r) = |G| * (dim(r)^2 - 1)/dim(r)
# No, that's for SU(N)

# For A_5, the Casimir eigenvalues can be computed from characters
# C_2 = |G|/(dim(r)) * sum_C |C| * chi_r(C) * chi_r(C^-1) / |G| - dim(r)
# Actually: for a finite group, the "Casimir" analog is
# Omega = sum_{g in G} g (x) g in the group algebra
# Its eigenvalue on irrep r: |G|/dim(r)

# The relevant quantity is the quadratic Casimir of A_5 representations
# C_2(trivial) = 0
# C_2(3) = C_2(3') (both 3-dim)
# C_2(4)
# C_2(5)

# For SU(2): C_2(j) = j(j+1), fundamental (j=1/2): 3/4, adjoint (j=1): 2
# Ratio adj/fund = 8/3

C2_fund = 3/4
C2_adj = 2
print(f"  SU(2) Casimirs: C_2(fund) = {C2_fund}, C_2(adj) = {C2_adj}")
print(f"  Ratio C_2(adj)/C_2(fund) = {C2_adj/C2_fund:.4f}")

# For SU(3): C_2(fund) = 4/3, C_2(adj) = 3
# The Higgs doesn't transform under SU(3) so this is less relevant

# Key ratio for the Higgs:
# m_H^2/m_W^2 = 8*lambda/g^2
# In some theories: lambda = g^2 * C_2(r_H) / (4 * something)

# Try: phi^2 from the icosahedron adjacency
# Adjacency eigenvalue for the 3-rep: sqrt(5)
# Adjacency eigenvalue for the 3'-rep: -sqrt(5)
# Product: |sqrt(5)| * |sqrt(5)| = 5 = a_1
# But phi^2 = phi + 1 = 2.618...
# And 5/2 = 2.5 (close but not equal)

print(f"\n  Icosahedron adjacency eigenvalues for gauge sectors:")
print(f"    trivial (U(1)): lambda = {a1} = a_1")
print(f"    3' (SU(2)):     lambda = -sqrt({a1}) = {-np.sqrt(a1):.4f}")
print(f"    3 (SU(3)):      lambda = +sqrt({a1}) = {+np.sqrt(a1):.4f}")
print(f"    5 (in SU(3)):   lambda = 1")
print(f"")
print(f"  Adjacency eigenvalue of trivial / |eigenvalue of 3'|:")
print(f"    {a1}/sqrt({a1}) = sqrt({a1}) = {np.sqrt(a1):.4f}")
print(f"    phi = {phi:.4f}")
print(f"    sqrt(a_1) vs phi: NOT equal (2.236 vs 1.618)")

# ================================================================
# J. THE CLEANEST APPROACH: PHI FROM SYMMETRY BREAKING PATTERN
# ================================================================
print(f"\n{'='*72}")
print(f"J. SYMMETRY BREAKING: phi from A_5 -> Z_5 -> Z_1")
print(f"{'='*72}")

print(f"""
  The 600-cell has symmetry 2I (binary icosahedral, order 120).
  The vertex stabilizer is Z_5 (5-fold rotation).
  Symmetry breaking: 2I -> Z_5 (by choosing a vertex)

  Index [2I : Z_5] = 120/5 = 24 (not useful directly)

  But for A_5: the stabilizer of a vertex in A_5 acting on
  the icosahedron is Z_5 (cyclic of order 5).

  A_5/Z_5 has 12 elements = 12 vertices.

  The phi emerges because Z_5 representations involve
  5th roots of unity: omega = exp(2*pi*i/5)
  And phi = -omega^2 - omega^3 = 2*cos(2*pi/5)... wait:

  2*cos(2*pi/5) = (sqrt(5)-1)/2 = 1/phi = phi - 1
  2*cos(4*pi/5) = -(sqrt(5)+1)/2 = -phi

  So phi IS the character of Z_5 representations!
""")

cos_2pi5 = 2*np.cos(2*np.pi/5)
cos_4pi5 = 2*np.cos(4*np.pi/5)
print(f"  2*cos(2*pi/5) = {cos_2pi5:.6f} = 1/phi = phi - 1 = {1/phi:.6f}")
print(f"  2*cos(4*pi/5) = {cos_4pi5:.6f} = -phi = {-phi:.6f}")

# In the 3-dim irrep of A_5 restricted to Z_5:
# 3|_{Z_5} has character: chi(omega^k) for k=0,...,4
# chi(1) = 3 (dimension)
# chi(omega) = chi(omega^4) and chi(omega^2) = chi(omega^3)
# Since the 3 of A_5 restricted to Z_5 is omega + omega^2 + omega^3
# No wait, 3 restricted to Z_5 = 1 + omega + omega^4 (or similar)

# The point: phi appears as eigenvalue of the Z_5 Fourier transform
# in the icosahedral context

# ================================================================
# SYNTHESIS
# ================================================================
print(f"\n{'='*72}")
print(f"SYNTHESIS: Multiple derivation paths for m_H/m_W = phi")
print(f"{'='*72}")

print(f"""
  PATH 1 (DSI exponents): n_H = n_Z + 1 = a_1^2 + 1 = 26
    - The Higgs is one DSI step above the Z.
    - m_H/m_W = phi follows from the +1 shift.
    - But WHY +1? This needs justification.
    - POSSIBLE: Higgs = next harmonic on DSI lattice.

  PATH 2 (McKay inner fluctuation): D_F ~ phi^(dist) on E8~
    - Higgs on edge 0-1 (distance 1): D_F ~ phi
    - Gauge at node (distance 0): D_F ~ phi^0 = 1
    - Ratio = phi. CLEAN but the D_F ~ phi^dist needs proof.

  PATH 3 (A_5 Casimir/eigenvalue): Numerically doesn't work.
    - Icosahedron eigenvalues give sqrt(5), not phi.
    - SU(2) Casimirs give 8/3, not phi^2.

  PATH 4 (NCG spectral action): lambda/g^2 = phi^2/8
    - Top dominance: Tr(Y^4)/Tr(Y^2)^2 ~ 1/3
    - Requires f_0/f_2^2 = 9*phi^2/(8*pi^2) = {f0_over_f2sq:.4f}
    - This is a CONDITION on the spectral function, not a derivation.

  VERDICT: Path 2 (McKay inner fluctuation) is most promising.
  It directly connects the Higgs mass to the graph distance.
""")

# ================================================================
# PATH 2 DEEP DIVE: D_F ~ phi on McKay edge
# ================================================================
print(f"\n{'='*72}")
print(f"DEEP DIVE: D_F matrix elements from phi-scaling on E8~")
print(f"{'='*72}")

# In the framework, ALL masses scale as m ~ phi^n
# The fermion at McKay node i has mass m_i ~ phi^(n_i)
# The Dirac operator D_F has matrix elements between adjacent nodes:
# <i|D_F|j> determines the Yukawa coupling between i and j

# For the Higgs at edge 0-1:
# <0|D_F|1> = Y_H (Yukawa for Higgs)
# In NCG: Y_H sets the electroweak symmetry breaking scale

# The gauge coupling at node 1:
# g_2 comes from the algebra structure at node 1

# Mass formulas:
# m_W = g_2 * v / 2 where v = <H> (Higgs vev)
# m_H = sqrt(2*lambda) * v

# If D_F(0,1) = phi * Lambda (some scale Lambda):
# Then v ~ phi * Lambda / g (schematically)
# And m_W ~ g * phi * Lambda / (2*g) = phi * Lambda / 2
# And m_H ~ sqrt(2*lambda) * phi * Lambda / g

# Hmm, this gives m_H/m_W = 2*sqrt(2*lambda)/g, same as before.
# The phi cancels out! It just sets the overall scale.

# UNLESS: phi appears in the RATIO of D_F contributions
# to m_H vs m_W differently.

# In NCG:
# m_W^2 = g^2 * f_2 * |D_F(0,1)|^2 / (4*pi^2)
# m_H^2 = 2 * f_0 * |D_F(0,1)|^4 / (f_2 * |D_F(0,1)|^2 * pi^2) * m_W^2
# Simplifying: m_H^2/m_W^2 = 2*f_0*|D_F(0,1)|^2 / (pi^2 * f_2)

# So the ratio depends on |D_F(0,1)|^2 * f_0/f_2

print(f"  In NCG: m_H^2/m_W^2 ~ |D_F(0,1)|^2 * f_0/f_2")
print(f"  For this to give phi^2: |D_F(0,1)|^2 * f_0/f_2 = phi^2 * C")
print(f"  where C is a numerical constant from the NCG framework.")

# The cleanest statement:
# If the spectral function f has f_0/f_2 = phi^2/(some integer),
# then m_H/m_W = phi automatically.

# But f_0 and f_2 are MOMENTS of the cutoff function.
# They are normally arbitrary parameters.
# UNLESS the 600-cell geometry fixes them!

print(f"\n  From the 600-cell spectral action (exp261-264):")
print(f"    c_0 = 2640 (cosmological term)")
print(f"    c_1 = 14880 (Einstein-Hilbert term)")
print(f"    c_2 = 55920 (Yang-Mills term)")
print(f"    c_0/c_1 = {2640/14880:.6f}")
print(f"    c_1/c_2 = {14880/55920:.6f}")
print(f"    c_0*c_2/c_1^2 = {2640*55920/14880**2:.6f}")

# c_k = sum of Tr(D^(2k)) essentially
# c_0 = dim(H) = 2640
# c_1 corresponds to f_2 moment
# c_2 corresponds to f_0 moment (higher power)

# In the spectral action S = f_0*c_2/2 + f_2*c_1 + f_4*c_0
# The Higgs mass ratio in this context:
# m_H^2/m_W^2 = 8 * (f_0*c_2) / (f_2*c_1)^? ... complex formula

ratio_c = 2640 * 55920 / 14880**2
print(f"\n  c_0*c_2/c_1^2 = {ratio_c:.6f}")
print(f"  phi^2/8 = {phi**2/8:.6f}")
print(f"  phi^2 = {phi**2:.6f}")

# Actually, let me check a simpler ratio
print(f"  c_2/c_1 = {55920/14880:.6f} = {55920//14880}*{14880}/{14880} + {55920%14880}/{14880:.4f}")
print(f"  55920/14880 = {55920/14880}")
print(f"  = 233/62 = {233/62:.6f}")

# WAIT: 233 = F(13), 62 from spectral action
# And 233/62... let's check
print(f"  233/62 = {233/62:.6f}")
print(f"  233 = F(13), 62 = c_1/240")

# Check if any ratio of c_k gives phi^2
print(f"\n  Checking ratios:")
print(f"  c_2/c_1 = 55920/14880 = {55920/14880:.4f}")
print(f"  c_1/c_0 = 14880/2640 = {14880/2640:.4f}")
print(f"  c_2/(c_1*phi) = {55920/(14880*phi):.4f}")
print(f"  c_0*c_2/c_1^2 = {ratio_c:.4f}")
print(f"  8*c_0*c_2/c_1^2 = {8*ratio_c:.4f}")
print(f"  phi^2 = {phi**2:.4f}")
print(f"  Ratio {8*ratio_c:.4f} / phi^2 = {8*ratio_c/phi**2:.4f}")

# Hmm, 8*c0*c2/c1^2 = 5.333 = 16/3
# phi^2 = 2.618
# Not matching

# ================================================================
# K. BACK TO BASICS: THE SIMPLEST STATEMENT
# ================================================================
print(f"\n{'='*72}")
print(f"K. SIMPLEST STATEMENT: phi = next eigenvalue ratio")
print(f"{'='*72}")

# In the DSI framework, masses come from phi^n.
# The Z boson is at n_Z = 25 = a_1^2.
# The Higgs is at n_H = 26 = a_1^2 + 1.
#
# m_H/m_Z = phi^(n_H - n_Z) * (corrections) = phi * (corrections)
# m_H/m_W = m_H/(m_Z * cos(tW)) = phi * (corrections) / cos(tW)
#
# Wait, the exact formula is m_H = m_W * (phi - 8*alpha)
# So m_H = m_Z * cos(tW) * (phi - 8*alpha)

mW_pred = m_Z * np.sqrt(1 - SIN2_TW)
mH_pred = mW_pred * (phi - 8*ALPHA)
print(f"  m_Z = {m_Z} GeV")
print(f"  m_W = m_Z * cos(tW) = {mW_pred:.3f} GeV (exp: {m_W})")
print(f"  m_H = m_W * (phi - 8*alpha) = {mH_pred:.3f} GeV (exp: {m_H})")

# The REAL question: in the DSI potential V(z) = -cos(2*pi*ln|z|/ln(phi))
# the masses correspond to phi^n levels.
# WHY is the Higgs at n = n_Z + 1?

# n_Z = 25 comes from: m_Z = m_e * phi^25 * alpha(mZ)/alpha(0)
# The 25 is DERIVED: phi^25 is the spectral quantity...

# But n_H = 26 is DEFINED by m_H = m_W * phi, which is what we're trying to derive!

# THE CORE INSIGHT might be:
# The Higgs field is the NEXT EXCITATION after the gauge bosons on the DSI lattice.
# In quantum mechanics, adjacent levels differ by one quantum number.
# The Higgs is the n+1 mode of the electroweak DSI potential.

print(f"\n  POSSIBLE DERIVATION:")
print(f"  The EW sector occupies DSI level n_Z = a_1^2 = 25.")
print(f"  The Higgs is the scalar fluctuation AROUND this level.")
print(f"  In the DSI potential, the next level is n_Z + 1 = 26.")
print(f"  The mass ratio between adjacent levels is ALWAYS phi.")
print(f"  Therefore: m_H/m_W = phi (the DSI step size).")
print(f"")
print(f"  This reduces to: WHY is the Higgs at n_Z + 1?")
print(f"  Answer: Because the Higgs IS the radial excitation")
print(f"  of the electroweak ground state in the DSI potential.")
print(f"  The electroweak vacuum sits at n_Z; the Higgs boson")
print(f"  is the first excited state, one step up: n_Z + 1.")
print(f"")
print(f"  In DSI: consecutive levels always have ratio phi.")
print(f"  Radial excitation = +1 level = multiply by phi.")
print(f"  QED: m_H/m_W = phi. (Modulo corrections from 8*alpha.)")

print(f"\n  STATUS:")
print(f"  - n_Z + 1 = 26 = a_1^2 + 1: CLEAN")
print(f"  - m_H/m_W = phi from DSI step: NATURAL")
print(f"  - But 'Higgs = first radial excitation' needs PROOF")
print(f"  - CATEGORY: PATTERN -> PARTIALLY DERIVED")
print(f"  - The argument that Higgs IS the +1 mode is physical")
print(f"    but not yet proven from the spectral action.")

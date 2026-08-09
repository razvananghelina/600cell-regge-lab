"""
exp362_higgs_8alpha.py
=======================
OPEN-8: Derive the Higgs mass correction -8*alpha.

KNOWN:
  m_H/m_W = phi - 8*alpha = 1.5597  (0.09% error with exp m_W)
  Tree-level: phi from L(3')/L(3) = phi^2 (exact, DERIVED)
  Correction: -8*alpha where 8 = rank(E8) = Tr(D_F^2)

APPROACHES:
  1. Spectral action on product space: does the 1-loop correction give -8*alpha?
  2. Tr(D_F^2) * alpha as perturbative correction
  3. Algebraic: phi - 8*alpha as a function of framework quantities
  4. NCG inner fluctuation: Higgs potential from D_F expansion
  5. Relation between 8*alpha and alpha_s/2 (observed ~1% agreement)

Dependencies: numpy, scipy
"""

import numpy as np
from scipy.linalg import eigh

a1 = 5
b1 = a1 + 1
phi = (1.0 + np.sqrt(a1)) / 2.0
phi_c = (1.0 - np.sqrt(a1)) / 2.0  # Galois conjugate
PI = np.pi

# Framework coupling constants
A_coef = 2 * PI
B_coef = -4 * a1 * phi**4
C_coef = 1.0
disc = B_coef**2 - 4 * A_coef * C_coef
alpha = (-B_coef - np.sqrt(disc)) / (2 * A_coef)
alpha_prime = (-B_coef + np.sqrt(disc)) / (2 * A_coef)  # Galois conjugate root
alpha_s = 1.0 / (2 * phi**3)
sin2tw = float(b1) / (a1**2 + 1)

# Experimental
m_W_exp = 80.377  # GeV
m_H_exp = 125.25  # GeV
ratio_exp = m_H_exp / m_W_exp

print("=" * 72)
print("EXP-362: HIGGS MASS CORRECTION -8*alpha MECHANISM")
print("=" * 72)

print(f"\n--- Framework constants ---")
print(f"  a1 = {a1}, b1 = {b1}")
print(f"  phi = {phi:.10f}")
print(f"  alpha = {alpha:.10f} (1/{1/alpha:.6f})")
print(f"  alpha' = {alpha_prime:.10f} (Galois conjugate)")
print(f"  alpha_s = {alpha_s:.10f}")
print(f"  rank(E8) = 8, Tr(D_F^2) = 8")

print(f"\n--- The Higgs formula ---")
print(f"  m_H/m_W = phi - 8*alpha = {phi - 8*alpha:.6f}")
print(f"  Experimental: {ratio_exp:.6f}")
print(f"  Error: {(phi - 8*alpha - ratio_exp)/ratio_exp*100:+.3f}%")
print(f"  Tree level (phi): {phi:.6f} ({(phi - ratio_exp)/ratio_exp*100:+.3f}%)")

# =====================================================================
# 1. ALGEBRAIC STRUCTURE OF THE CORRECTION
# =====================================================================
print(f"\n{'=' * 72}")
print("1. ALGEBRAIC STRUCTURE")
print(f"{'=' * 72}")

print(f"\n  8*alpha = {8*alpha:.10f}")
print(f"  alpha_s/2 = {alpha_s/2:.10f}")
print(f"  Ratio: 8*alpha/(alpha_s/2) = {8*alpha/(alpha_s/2):.6f}")
print(f"  Difference: {(8*alpha - alpha_s/2)/(alpha_s/2)*100:+.3f}%")

# Is there an exact relation?
# 8*alpha = 8 * (smaller root of 2*pi*x^2 - 20*phi^4*x + 1 = 0)
# alpha_s/2 = 1/(4*phi^3)
# Are these equal? 8*alpha = 1/(4*phi^3)?
# If so: alpha = 1/(32*phi^3) = alpha_s/16
# Check: alpha_s/16 = 0.007377, actual alpha = 0.007297. Not equal (1.1% off).

print(f"\n  alpha_s/16 = {alpha_s/16:.10f}")
print(f"  alpha     = {alpha:.10f}")
print(f"  NOT equal. Close (1.1%) but distinct.")

# Check: 8*alpha in terms of alpha equation
# From 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0:
# alpha = (20*phi^4 - sqrt(400*phi^8 - 8*pi)) / (4*pi)
# 8*alpha = 2*(20*phi^4 - sqrt(400*phi^8 - 8*pi)) / pi

val_8a = 2*(20*phi**4 - np.sqrt(400*phi**8 - 8*PI)) / PI
print(f"\n  8*alpha = 2*(20*phi^4 - sqrt(400*phi^8 - 8*pi))/pi = {val_8a:.10f}")
print(f"  Check: {8*alpha:.10f}")

# phi - 8*alpha = phi - 2*(20*phi^4 - sqrt(400*phi^8 - 8*pi))/pi
# = phi - 40*phi^4/pi + 2*sqrt(400*phi^8 - 8*pi)/pi
target = phi - 8*alpha
print(f"\n  phi - 8*alpha = {target:.10f}")
print(f"  (phi - 8*alpha)^2 = {target**2:.10f}")
print(f"  phi^2 = {phi**2:.10f}")
print(f"  phi^2 - (phi-8*alpha)^2 = {phi**2 - target**2:.10f}")
print(f"  = 16*alpha*(phi - 4*alpha) = {16*alpha*(phi - 4*alpha):.10f}")

# =====================================================================
# 2. SPECTRAL ACTION: Tr(D_F^2) AND THE CORRECTION
# =====================================================================
print(f"\n{'=' * 72}")
print("2. SPECTRAL ACTION CORRECTION FROM Tr(D_F^2)")
print(f"{'=' * 72}")

# In NCG, the Higgs mass formula involves the ratio Tr(Y^4)/Tr(Y^2)
# The "Y" is the Yukawa matrix on the finite space
# In our framework, D_F has Tr(D_F^2) = 8 = rank(E8)

# Hypothesis: the correction is a 1-loop effect
# delta(m_H^2/m_W^2) = -alpha * Tr(D_F^2) * (some geometric factor)

# For m_H/m_W = phi*(1 - delta/2) ~ phi - phi*delta/2:
# phi*delta/2 = 8*alpha
# delta = 16*alpha/phi

# Alternatively: m_H^2/m_W^2 = phi^2*(1 - delta)
# (phi-8*alpha)^2 ~ phi^2 - 16*phi*alpha
# delta = 16*phi*alpha/phi^2 = 16*alpha/phi

delta_needed = 16*alpha/phi
print(f"\n  For m_H^2/m_W^2 correction:")
print(f"  delta = 16*alpha/phi = {delta_needed:.6f}")
print(f"  = 2 * Tr(D_F^2) * alpha / phi")
print(f"  = 2 * rank(E8) * alpha / phi")

# In perturbation theory, a 1-loop correction to a mass ratio has the form:
# delta ~ (coupling)/(4*pi) * (group theory factor) * log(cutoff/mass)
# or on a discrete space: delta ~ coupling * (spectral factor)

# Check: does alpha/(4*pi) * something = 16*alpha/phi?
# something = 64*pi/phi = 124.2
# That's too large for a perturbative coefficient

# But on a FINITE space, there's no 1/(4*pi) suppression!
# The correction would be: delta = alpha * C_spectral
# where C_spectral is a pure number from the spectral geometry

C_spectral = 16.0/phi
print(f"\n  Required C_spectral = 16/phi = {C_spectral:.6f}")
print(f"  = 2*Tr(D_F^2)/phi")

# Possible origin: In the spectral action expansion
# S = Tr(f(D^2/Lambda^2)) = c_0*f_4 + c_1*f_2 + c_2*f_0 + ...
# The Higgs terms come from expanding D^2 = (D_M + D_F)^2
# On the product space: D^2 = D_M^2 x 1_F + 1_M x D_F^2 + cross terms

# The cross terms involve gamma_M x D_F and give the Yukawa couplings
# The D_F^2 term gives a mass-squared correction to the Higgs

# Specifically: Tr(f(D^2)) expanded around D_M^2:
# ~ Tr(f(D_M^2)) * dim(F) + f'(D_M^2) * Tr(D_F^2) + ...
# The Tr(D_F^2) = 8 correction shifts the effective potential

# =====================================================================
# 3. NCG HIGGS POTENTIAL ON THE PRODUCT SPACE
# =====================================================================
print(f"\n{'=' * 72}")
print("3. NCG HIGGS POTENTIAL ON PRODUCT SPACE")
print(f"{'=' * 72}")

# Build the McKay graph (affine E8) adjacency with CORRECT construction
# Using character inner products (not hardcoded!)

print(f"\n  Building McKay graph from 2I character table...")

# 2I character table (9 irreps)
# Conjugacy classes: 1, -1, 10A, 10B, 12A, 12B, 6A, 4A, 4B
# Order:             1,  2,   10,   10,   12,   12,   6,   4,   4
# Size:              1,  1,   12,   12,   20,   20,   12,  30,  12

# Characters of the 9 irreps (from 2I = <2,3,5>):
# rho_0 (dim 1): trivial
# rho_1 (dim 2): fundamental (z)
# rho_2 (dim 3): Sym^2(z)
# rho_3 (dim 4): Sym^3(z)
# rho_4 (dim 5): Sym^4(z)
# rho_5 (dim 6): Sym^5(z) / ...
# rho_6 (dim 4): z x z'
# rho_7 (dim 2): z'
# rho_8 (dim 3): Sym^2(z')

class_sizes = np.array([1, 1, 12, 12, 20, 20, 12, 30, 12])
order = 120

# Character values
w5 = np.exp(2j*PI/5)
w10 = np.exp(2j*PI/10)
w12 = np.exp(2j*PI/12)

# Fundamental (dim 2): chi = w^k + w^{-k} for class of order n
# Conjugacy classes by trace of SU(2) element:
# 1: tr=2, -1: tr=-2, 10A: tr=phi, 10B: tr=-1/phi,
# 12A: tr=sqrt(3), 12B: tr=-sqrt(3), 6A: tr=1, 4A: tr=0, 4B: tr=-1

traces = {
    '1': 2.0,
    '-1': -2.0,
    '10A': phi,
    '10B': -1.0/phi,
    '12A': np.sqrt(3),
    '12B': -np.sqrt(3),
    '6A': 1.0,
    '4A': 0.0,
    '4B': -1.0
}

# Build character table from Sym^k decomposition
# chi_Sym^k(g) is computed from tr(g) = t using:
# chi_1(t) = 1
# chi_2(t) = t  (= fundamental)
# chi_{k+1}(t) = t * chi_k(t) - chi_{k-1}(t)  (Chebyshev recurrence)

t_vals = np.array([2.0, -2.0, phi, -1.0/phi, np.sqrt(3), -np.sqrt(3), 1.0, 0.0, -1.0])
n_classes = len(t_vals)

# Compute Sym^k characters up to k=5 using Chebyshev recurrence
# S_0 = 1, S_1 = t, S_{k+1} = t*S_k - S_{k-1}
S = np.zeros((10, n_classes))
S[0] = 1.0  # dim 1 (trivial)
S[1] = t_vals  # dim 2 (fundamental z)

for k in range(1, 9):
    S[k+1] = t_vals * S[k] - S[k-1]

# The 9 irreps and their characters:
# rho_0 = S_0 (dim 1)
# rho_1 = S_1 (dim 2)  [= z, fundamental]
# rho_2 = S_2 (dim 3)  [= Sym^2(z)]
# rho_3 = S_3 (dim 4)  [= Sym^3(z)]
# rho_4 = S_4 (dim 5)  [= Sym^4(z)]
# rho_5 = S_5 (dim 6)  [this is reducible for 2I! Sym^5 decomposes]
# Actually, for 2I, Sym^k(z) for k <= 4 gives irreps of dim k+1
# For k >= 5, they start decomposing.
# The remaining irreps come from z' (the other fundamental) and its Sym^k

# Let me use the KNOWN character table of 2I directly:
# Classes: 1, -1, C10a, C10b, C12a, C12b, C6, C4, C4'
# Using traces as computed above

# From the character table of 2I (binary icosahedral group):
# dim:  1   2   3   4   5   6   4   2   3
chi_table = np.zeros((9, 9))

# rho_0 (trivial, dim 1)
chi_table[0] = [1, 1, 1, 1, 1, 1, 1, 1, 1]

# rho_1 (dim 2, fundamental z): chi = trace
chi_table[1] = t_vals

# rho_2 (dim 3, Sym^2 z)
chi_table[2] = S[2]

# rho_3 (dim 4, Sym^3 z)
chi_table[3] = S[3]

# rho_4 (dim 5, Sym^4 z)
chi_table[4] = S[4]

# For rho_5 through rho_8, we need to be more careful
# The key: z x z (tensor product of fundamental with itself)
# z x z = Sym^2(z) + Alt^2(z) = rho_2 + rho_0 (since Alt^2 of 2-dim = trivial)
# So tensor products decompose using characters

# For rho_7 (dim 2, z'): chi(g) = chi_1(g) * det(g)^{-1}
# In 2I, det of SU(2) elements = 1 always
# z' is obtained from z by the outer automorphism
# For 2I specifically: z' has chi(10A) = -1/phi, chi(10B) = phi
# (swapping the two 10-element classes)

# rho_7 (dim 2): z' (swaps 10A <-> 10B)
chi_table[7] = [2, -2, -1.0/phi, phi, np.sqrt(3), -np.sqrt(3), 1.0, 0.0, -1.0]

# rho_8 (dim 3): Sym^2(z')
t_prime = chi_table[7]
chi_table[8] = t_prime**2 - np.array([1,1,1,1,1,1,1,1,1])  # S_2(t') = t'^2 - 1

# rho_6 (dim 4): z x z'
chi_table[6] = chi_table[1] * chi_table[7]

# rho_5 (dim 6): need to find. Total sum of dim^2 = 1+4+9+16+25+36+16+4+9 = 120
# Sym^5(z) = dim 6 irrep... let's check S[5]
# S_5 for 2I: should give dim 6 character
# But Sym^5(z) may not be irreducible! For SU(2), it is irreducible.
# For 2I (finite subgroup), Sym^5(z) is irreducible if dim <= |2I|/2 = 60.
# dim(Sym^5) = 6, so yes it should be irreducible.

chi_table[5] = S[5]

# Verify: dimensions
dims = chi_table[:, 0].astype(int)
print(f"  Dimensions: {list(dims)}")
print(f"  Sum dim^2 = {sum(d**2 for d in dims)} (should be 120)")

# Verify orthogonality
print(f"\n  Orthogonality check (should be delta_ij * 120):")
for i in range(9):
    for j in range(i, min(i+3, 9)):
        inner = sum(class_sizes[k] * chi_table[i,k] * chi_table[j,k] for k in range(9))
        if abs(inner) > 0.01 or i == j:
            status = "OK" if (i == j and abs(inner - 120) < 0.01) or (i != j and abs(inner) < 0.01) else "FAIL"
            if i == j or abs(inner) > 0.01:
                print(f"    <rho_{i}|rho_{j}> = {inner:.2f} {status}")

# Build McKay graph: tensor product of each rho_k with rho_1 (fundamental)
print(f"\n  McKay graph edges (rho_k x rho_1 decomposition):")
mckay_adj = np.zeros((9, 9), dtype=int)
for k in range(9):
    product_chi = chi_table[k] * chi_table[1]  # chi(rho_k x rho_1)
    for j in range(9):
        mult = sum(class_sizes[c] * product_chi[c] * chi_table[j,c] for c in range(9)) / order
        mult_int = int(round(mult))
        if mult_int > 0:
            mckay_adj[k, j] = mult_int

print(f"  McKay adjacency matrix:")
for i in range(9):
    row = [str(mckay_adj[i,j]) for j in range(9)]
    print(f"    rho_{i} (dim {dims[i]}): [{', '.join(row)}]")

# Identify edges
edges = []
for i in range(9):
    for j in range(i+1, 9):
        if mckay_adj[i, j] > 0:
            edges.append((i, j, mckay_adj[i,j]))
print(f"\n  Edges: {[(f'rho_{i}-rho_{j}', w) for i,j,w in edges]}")

# Build D_F from McKay edge weights
# D_F has structure matching the McKay adjacency (with appropriate weights)
# For now, use the mass exponents as diagonal and McKay edges as off-diagonal

# McKay graph distances from rho_0
from collections import deque
def bfs_dist(adj, start):
    n = adj.shape[0]
    dist = [-1]*n
    dist[start] = 0
    q = deque([start])
    while q:
        u = q.popleft()
        for v in range(n):
            if adj[u,v] > 0 and dist[v] == -1:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist

dist_from_0 = bfs_dist(mckay_adj, 0)
print(f"\n  Distances from rho_0: {dist_from_0}")

# =====================================================================
# 4. THE D_F MATRIX AND HIGGS DIRECTION
# =====================================================================
print(f"\n{'=' * 72}")
print("4. D_F MATRIX AND HIGGS DIRECTION")
print(f"{'=' * 72}")

# In NCG, D_F is an off-diagonal operator on the fiber Hilbert space
# H_F = direct sum of rho_k (total dim = sum dims = 30)
# D_F connects representations that are adjacent in the McKay graph

# The Wilson line masses (from exp352):
# Edge weights w_e such that m_f = m_e * phi^{sum w_e along path}
# Solution 3 weights: {0, 1, 2, 3, 3, 5, 6, 9}
# The mass of each fermion = product of phi^{w_e} along path from rho_0

# Wilson line edge weights (Solution 3 from exp352)
# These need to be assigned to the 8 edges of the McKay tree
# The assignment depends on which edges correspond to which weights

# From the mass formula n_f = 5a + 6b:
# The exponents are: e(0), u(3), d(5), s(11), c(16), b(19), t(26), tau(17), mu(11)
# Path from rho_0: the sum of weights along the unique path gives n_f

# McKay assignment (from chirality work):
# rho_0=trivial, rho_1=fund(2), rho_2=Sym^2(3), rho_3=Sym^3(4),
# rho_4=Sym^4(5), rho_5=Sym^5(6), rho_6=z*z'(4), rho_7=z'(2), rho_8=Sym^2(z')(3)

# Fermion assignment (from exp351):
# WHITE: rho_0(e), rho_4(d), rho_5(b), rho_6(tau), rho_8(mu)
# BLACK: rho_1(u), rho_3(t), rho_7(s), rho_2(c)  ... wait, need to check

# Actually, let me just use the GRAPH structure and edge weights
# McKay tree structure: 0-1-2-3-4-5-6, branch: 6-7, 6-8
# Wait, need to check. Let me read the edge list.

# From the computed adjacency, the edges with weight 1 are:
actual_edges = [(i,j) for i,j,w in edges if w == 1]
print(f"  McKay tree edges: {actual_edges}")

# E8~ affine Dynkin:
# Standard: 0-1-2-3-4-5-6-7, branch at 4 connecting to 8
# But OUR convention may differ. Let me check the graph structure.

# Check degrees
degrees_mckay = mckay_adj.sum(axis=1)
print(f"  Degrees: {[(f'rho_{i}(dim{dims[i]})', int(degrees_mckay[i])) for i in range(9)]}")

# The branch point has degree 3
branch_nodes = [i for i in range(9) if degrees_mckay[i] == 3]
print(f"  Branch point: rho_{branch_nodes[0]} (dim {dims[branch_nodes[0]]})")

# Legs from branch point
branch = branch_nodes[0]
neighbors_of_branch = [j for j in range(9) if mckay_adj[branch, j] > 0]
print(f"  Neighbors of branch: {[f'rho_{j}(dim{dims[j]})' for j in neighbors_of_branch]}")

# =====================================================================
# 5. HIGGS AS INNER FLUCTUATION AND THE CORRECTION
# =====================================================================
print(f"\n{'=' * 72}")
print("5. HIGGS AS INNER FLUCTUATION")
print(f"{'=' * 72}")

# In NCG, the Higgs arises from inner fluctuations of D along the fiber
# The key edge is rho_0 -- rho_1 (trivial to fundamental)
# This edge connects C (dim 1) to M_2(C) (dim 2) = SU(2) sector

# The Higgs doublet has 4 real DOF = dim(rho_3) or dim(rho_6)
# Actually in NCG: Higgs = inner fluctuation along the edge connecting
# the "lepton" algebra to the "quark" algebra

# For the m_H/m_W ratio:
# Tree level: phi from L(3')/L(3) = phi^2 on the icosahedron (DERIVED)
# This is a property of the BASE MANIFOLD (600-cell), not the fiber

# The correction must come from the FIBER (McKay graph)
# Specifically: the spectral action on the product space gives
# corrections proportional to Tr(D_F^{2k}) / Tr(D_M^{2k})

# On the 600-cell: Tr(D_M^2) = c_1 = 14880
# On McKay: Tr(D_F^2) = 8 = rank(E8)
# Ratio: Tr(D_F^2)/Tr(D_M^2) = 8/14880 = 1/1860

# The correction to the Higgs mass squared:
# delta(m_H^2/m_W^2) = -f(Tr(D_F^2), alpha)

# Hypothesis: delta = Tr(D_F^2) * alpha / phi_factor
# For delta = 16*alpha/phi (needed):
# 16*alpha/phi = Tr(D_F^2) * (2*alpha/phi)
# So the coupling is 2*alpha/phi per unit of Tr(D_F^2)

print(f"  Tr(D_F^2) = 8 (= rank E8)")
print(f"  Need: correction = 16*alpha/phi = {16*alpha/phi:.6f}")
print(f"  Per unit Tr: 2*alpha/phi = {2*alpha/phi:.6f}")

# =====================================================================
# 6. KEY IDEA: SELF-ENERGY ON THE McKAY GRAPH
# =====================================================================
print(f"\n{'=' * 72}")
print("6. HIGGS SELF-ENERGY ON McKAY GRAPH")
print(f"{'=' * 72}")

# On the McKay graph, the Higgs "lives" at the edge rho_0 -- rho_1.
# Its self-energy receives contributions from ALL other edges (loops on the tree).
#
# In the spectral action framework:
# S_Higgs = <H, (D_F^2 + gauge corrections) H>
# The gauge correction to D_F^2 involves the gauge coupling alpha.
#
# On a tree graph, a "1-loop" correction to an edge involves:
# - The edge itself (mass renormalization)
# - Paths of length 3 through the edge (triangle = 0 on tree!)
# - Paths of length 5 through the edge (pentagon)
#
# On the McKay tree with 8 edges:
# The correction to the Higgs edge (0-1) from gauge fluctuations is:
# delta_H = alpha * sum_over_paths * path_weight

# For a tree: no triangles, so the leading correction comes from
# paths of length 2 through the Higgs vertex (vertex corrections)

# Vertex corrections at rho_0: only one edge (0-1), so no vertex correction
# Vertex corrections at rho_1: edges (0-1) and (1-2), path 0-1-2

# Self-energy on edge 0-1:
# Sigma = alpha * [vertex correction at 0] + alpha * [vertex correction at 1]
# = alpha * 0 + alpha * (contribution from edge 1-2)

# The total correction involves the FULL propagator on the McKay graph
# which sums over all paths. On a tree, this is the Green's function G(i,j).

# McKay graph Laplacian
L_mckay = np.diag(degrees_mckay.astype(float)) - mckay_adj.astype(float)
eig_L, vec_L = eigh(L_mckay)
print(f"  McKay Laplacian eigenvalues: {[f'{e:.4f}' for e in eig_L]}")

# Green's function (regularized, excluding zero mode)
# G = sum_{k>0} |v_k><v_k| / lambda_k
G = np.zeros((9, 9))
for k in range(1, 9):  # skip zero mode
    G += np.outer(vec_L[:, k], vec_L[:, k]) / eig_L[k]

print(f"\n  Green's function G(0,0) = {G[0,0]:.6f}")
print(f"  Green's function G(1,1) = {G[1,1]:.6f}")
print(f"  Green's function G(0,1) = {G[0,1]:.6f}")
print(f"  G(0,0) - 2*G(0,1) + G(1,1) = {G[0,0] - 2*G[0,1] + G[1,1]:.6f}")
print(f"  (This is the effective resistance of edge 0-1)")

# Effective resistance of all edges
print(f"\n  Effective resistances of McKay edges:")
for i, j, w in edges:
    R_eff = G[i,i] - 2*G[i,j] + G[j,j]
    print(f"    rho_{i}-rho_{j}: R_eff = {R_eff:.6f}")

# Total effective resistance (Kirchhoff index)
K = sum(G[i,i] - 2*G[i,j] + G[j,j] for i,j,_ in edges)
print(f"  Kirchhoff index (sum R_eff) = {K:.6f}")

# =====================================================================
# 7. SPECTRAL ACTION COEFFICIENTS AND HIGGS MASS
# =====================================================================
print(f"\n{'=' * 72}")
print("7. SPECTRAL ACTION COEFFICIENTS")
print(f"{'=' * 72}")

# The 600-cell spectral action coefficients
c0 = 2640.0  # dim(H_M) = (720+1200+600+120) = 2640
c1 = 14880.0
c2 = 55920.0

# The NCG Higgs mass formula (Chamseddine-Connes):
# m_H^2 = (2*f_0/f_2) * b/a * v^2
# where a = Tr(Y^*Y), b = Tr((Y^*Y)^2)
# f_0 and f_2 are spectral moments

# On the 600-cell: f_2 relates to c_1, f_0 relates to c_2
# The exact relation depends on normalization

# Key ratio:
print(f"  c_0 = {c0:.0f}")
print(f"  c_1 = {c1:.0f}")
print(f"  c_2 = {c2:.0f}")
print(f"  c_2/c_1 = {c2/c1:.6f}")
print(f"  c_1/c_0 = {c1/c0:.6f}")
print(f"  c_0*c_2/c_1^2 = {c0*c2/c1**2:.6f}")

# The Higgs mass in terms of spectral coefficients:
# m_H^2/m_W^2 = C * (c_2/c_1) * (Tr(D_F^4)/Tr(D_F^2))

TrDF2 = 8.0  # = rank(E8), from ||Y||_F = 1/sqrt(2) and 9 fermions
# Tr(D_F^4) depends on the specific D_F matrix
# For D_F with eigenvalues proportional to mass exponents n_f:
# This is dominated by the top quark

# But wait - in our framework, Tr(D_F^2) = 8 is EXACT (from McKay)
# regardless of the specific eigenvalues!

# Let me check: what ratio of c_k gives phi^2?
print(f"\n  Looking for phi^2 = {phi**2:.6f} in spectral ratios:")
print(f"  c_2/(c_1 * sqrt(c_0/c_2)) = {c2/(c1*np.sqrt(c0/c2)):.6f}")
print(f"  sqrt(c_2/c_0) = {np.sqrt(c2/c0):.6f}")
print(f"  (c_2/c_1)/(c_1/c_0) = {(c2/c1)/(c1/c0):.6f}")
print(f"  c_0*c_2/c_1^2 * 4 = {4*c0*c2/c1**2:.6f}")
print(f"  c_2/(a1*c_1) = {c2/(a1*c1):.6f}")
print(f"  (c_2/c_1)/a1 = {(c2/c1)/a1:.6f}")

# =====================================================================
# 8. DIRECT ALGEBRAIC APPROACH
# =====================================================================
print(f"\n{'=' * 72}")
print("8. DIRECT ALGEBRAIC IDENTITIES")
print(f"{'=' * 72}")

# Check if phi - 8*alpha can be expressed as a closed form
# involving only a1, phi, pi

# From the alpha equation: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
# So: alpha = (1 + 2*pi*alpha^2)/(4*a1*phi^4)  [from rearranging]
# And: 8*alpha = 2*(1 + 2*pi*alpha^2)/(a1*phi^4)
# = 2/(a1*phi^4) + 4*pi*alpha^2/(a1*phi^4)
# The second term is O(alpha^2) ~ 10^{-4}, negligible

# Leading order: 8*alpha ~ 2/(a1*phi^4) = 2/(5*phi^4)
val_approx = 2.0/(a1*phi**4)
print(f"  8*alpha (exact)  = {8*alpha:.10f}")
print(f"  2/(a1*phi^4)     = {val_approx:.10f}")
print(f"  Difference: {(8*alpha - val_approx)/8/alpha*100:.4f}%")

# So to leading order:
# phi - 8*alpha ~ phi - 2/(a1*phi^4)
# = phi - 2/(a1*(phi^2)^2)
# = phi - 2/(a1*(phi+1)^2)

# Since phi^4 = 3*phi + 2:
# phi - 2/(a1*(3*phi+2))
# = phi - 2/(15*phi + 10)
# = (phi*(15*phi+10) - 2)/(15*phi+10)
# = (15*phi^2 + 10*phi - 2)/(15*phi+10)
# = (15*(phi+1) + 10*phi - 2)/(15*phi+10)
# = (25*phi + 13)/(15*phi+10)

val_frac = (25*phi + 13)/(15*phi + 10)
print(f"\n  Leading order: phi - 2/(a1*phi^4)")
print(f"  = (25*phi+13)/(15*phi+10) = {val_frac:.10f}")
print(f"  Exact phi-8*alpha        = {phi - 8*alpha:.10f}")
print(f"  Difference: {abs(val_frac - (phi-8*alpha)):.2e}")

# Even more: 25 = a1^2, 13 = (a1^2+1)/2, 15 = 3*a1, 10 = 2*a1
# So: (a1^2*phi + (a1^2+1)/2) / (3*a1*phi + 2*a1)
# = (a1^2*phi + (a1^2+1)/2) / (a1*(3*phi + 2))
# = (a1^2*phi + (a1^2+1)/2) / (a1*phi^4)  [since 3*phi+2 = phi^4]

# Numerator: a1^2*phi + (a1^2+1)/2 = 25*1.618 + 13 = 40.45 + 13 = 53.45
# Denominator: a1*phi^4 = 5*6.854 = 34.27

print(f"\n  Expressing in framework terms:")
print(f"  Numerator: a1^2*phi + (a1^2+1)/2 = {a1**2*phi + (a1**2+1)/2:.6f}")
print(f"  Denominator: a1*phi^4 = {a1*phi**4:.6f}")

# This can be further simplified:
# = (2*a1^2*phi + a1^2 + 1) / (2*a1*phi^4)
# = (a1^2*(2*phi+1) + 1) / (2*a1*phi^4)
# Note: 2*phi + 1 = 1 + sqrt(5) + 1 = 2 + sqrt(5) = sqrt(5) + 2

# Wait, 2*phi + 1 = 2*(1+sqrt(5))/2 + 1 = 1 + sqrt(5) + 1 = 2 + sqrt(5) = phi + phi^2 = phi^3
# YES! 2*phi + 1 = phi^3 (since phi^3 = phi^2 + phi = (phi+1) + phi = 2*phi+1)

print(f"\n  KEY: 2*phi + 1 = phi^3 = {2*phi+1:.6f} vs {phi**3:.6f}")
print(f"  So numerator = a1^2*phi^3 + 1")
print(f"  Leading order: phi - 8*alpha ~ (a1^2*phi^3 + 1) / (2*a1*phi^4)")
print(f"  = a1/(2*phi) + 1/(2*a1*phi^4)")
print(f"  = a1/(2*phi) + alpha/4  (since alpha ~ 1/(2*a1*phi^4) leading)")

# Check: a1/(2*phi) = 5/(2*1.618) = 5/3.236 = 1.5451
# 1/(2*a1*phi^4) = 1/34.27 = 0.02918
# alpha/4 = 0.001824
# Hmm, these don't add up right. Let me recalculate.

# Actually: (a1^2*phi^3 + 1)/(2*a1*phi^4)
# = a1*phi^3/(2*phi^4) + 1/(2*a1*phi^4)
# = a1/(2*phi) + 1/(2*a1*phi^4)
val_check = a1/(2*phi) + 1/(2*a1*phi**4)
print(f"  a1/(2*phi) = {a1/(2*phi):.6f}")
print(f"  1/(2*a1*phi^4) = {1/(2*a1*phi**4):.6f}")
print(f"  Sum = {val_check:.6f}")
print(f"  phi - 8*alpha = {phi - 8*alpha:.6f}")
print(f"  Difference: {abs(val_check - (phi - 8*alpha)):.2e}")

# =====================================================================
# 9. THE VIETA APPROACH
# =====================================================================
print(f"\n{'=' * 72}")
print("9. VIETA'S FORMULAS AND THE HIGGS CORRECTION")
print(f"{'=' * 72}")

# From the alpha equation: alpha + alpha' = 4*a1*phi^4/(2*pi) = 2*a1*phi^4/pi
# And: alpha * alpha' = 1/(2*pi)

# So: 8*alpha = 8*alpha
# And: 8*alpha' = 8*alpha'

# phi - 8*alpha = phi - 8*alpha
# phi - 8*alpha' = phi - 8*alpha' (Galois conjugate version)

# Product: (phi - 8*alpha)(phi - 8*alpha')
# = phi^2 - 8*phi*(alpha+alpha') + 64*alpha*alpha'
# = phi^2 - 8*phi*2*a1*phi^4/pi + 64/(2*pi)
# = phi^2 - 16*a1*phi^5/pi + 32/pi

product_vieta = (phi - 8*alpha) * (phi - 8*alpha_prime)
product_formula = phi**2 - 16*a1*phi**5/PI + 32/PI

print(f"  alpha + alpha' = {alpha + alpha_prime:.10f}")
print(f"  2*a1*phi^4/pi = {2*a1*phi**4/PI:.10f}")
print(f"  alpha * alpha' = {alpha * alpha_prime:.10f}")
print(f"  1/(2*pi) = {1/(2*PI):.10f}")
print(f"")
print(f"  (phi-8*alpha)(phi-8*alpha') = {product_vieta:.10f}")
print(f"  phi^2 - 16*a1*phi^5/pi + 32/pi = {product_formula:.10f}")

# Simplify: phi^5 = 5*phi + 3 (from phi^5 = phi^4*phi = (3*phi+2)*phi = 3*phi^2+2*phi = 3(phi+1)+2*phi = 5*phi+3)
# So: -16*a1*phi^5/pi = -16*5*(5*phi+3)/pi = -80*(5*phi+3)/pi
# = (-400*phi - 240)/pi

print(f"\n  phi^5 = 5*phi + 3 = {5*phi+3:.6f} vs {phi**5:.6f}")

# Full product:
# = phi^2 + (-400*phi - 240 + 32)/pi
# = phi^2 + (-400*phi - 208)/pi
# = phi^2 - (400*phi + 208)/pi

val_product = phi**2 - (400*phi + 208)/PI
print(f"  Product = phi^2 - (400*phi + 208)/pi = {val_product:.10f}")
print(f"  Direct: {product_vieta:.10f}")

# Hmm, this doesn't simplify nicely. The transcendental pi prevents closure.

# =====================================================================
# 10. THE CORRECT FRAMEWORK FORMULA
# =====================================================================
print(f"\n{'=' * 72}")
print("10. SEARCHING FOR THE CORRECT FORMULA")
print(f"{'=' * 72}")

# Let me check systematically: is there an expression involving only
# framework quantities that gives phi - 8*alpha to very high precision?

target = phi - 8*alpha
print(f"  Target: phi - 8*alpha = {target:.10f}")

candidates = [
    ("a1/(2*phi) + alpha/2", a1/(2*phi) + alpha/2),
    ("phi*(1 - Tr(D_F^2)*alpha/phi)", phi*(1 - TrDF2*alpha/phi)),
    ("phi - rank(E8)*alpha", phi - 8*alpha),  # tautology, for reference
    ("phi*(1 - 2/(a1*phi^3))", phi*(1 - 2/(a1*phi**3))),
    ("phi*(1 - alpha_s/phi)", phi*(1 - alpha_s/phi)),
    ("phi - alpha_s/2", phi - alpha_s/2),
    ("phi - 1/(4*phi^3)", phi - 1/(4*phi**3)),
    ("phi^2 / (phi + 8*alpha)", phi**2 / (phi + 8*alpha)),  # rearrangement
    ("(phi^2 - 8*alpha*phi)", phi**2 - 8*alpha*phi),
    ("sqrt(phi^2 - 16*alpha*phi)", np.sqrt(phi**2 - 16*alpha*phi)),
    ("phi*cos(8*alpha)", phi*np.cos(8*alpha)),
    ("phi*exp(-8*alpha/phi)", phi*np.exp(-8*alpha/phi)),
    ("(a1^2*phi^3+1)/(2*a1*phi^4)", (a1**2*phi**3+1)/(2*a1*phi**4)),
    ("phi*(phi^4-2)/(phi^4)", phi*(phi**4-2)/phi**4),
    ("phi - 2*phi^{-4}/a1", phi - 2*phi**(-4)/a1),
]

print(f"\n  {'Expression':45s} {'Value':14s} {'Err from target':14s}")
for name, val in candidates:
    err = val - target
    print(f"  {name:45s} {val:14.10f} {err:+14.2e}")

# =====================================================================
# 11. THE alpha_s CONNECTION
# =====================================================================
print(f"\n{'=' * 72}")
print("11. CONNECTION TO alpha_s")
print(f"{'=' * 72}")

# We observed: 8*alpha ~ alpha_s/2 (1.1% accuracy)
# Let's check more carefully

print(f"  8*alpha   = {8*alpha:.10f}")
print(f"  alpha_s/2 = {alpha_s/2:.10f}")
print(f"  Ratio: {8*alpha/(alpha_s/2):.8f}")
print(f"  Difference: {(8*alpha - alpha_s/2):.6e}")

# From alpha equation: alpha ~ 1/(2*pi * 2*a1*phi^4/pi) = 1/(4*a1*phi^4)
# Exact leading order: alpha ~ 1/(4*a1*phi^4) = 1/(20*phi^4)
# alpha_s = 1/(2*phi^3)
# 8*alpha ~ 8/(20*phi^4) = 2/(5*phi^4) = 2/(5*phi*phi^3) = 2*alpha_s*phi/(5*phi)
# Wait: = 2/(5*phi^4) and alpha_s/2 = 1/(4*phi^3)
# Ratio: (2/(5*phi^4)) / (1/(4*phi^3)) = 8/(5*phi) = 8/8.09 = 0.989
# So 8*alpha / (alpha_s/2) ~ 8/(5*phi) = 0.989

print(f"\n  Leading order: 8*alpha ~ 2/(a1*phi^4)")
print(f"  alpha_s/2 = 1/(4*phi^3)")
print(f"  Ratio (leading) = 2/(a1*phi^4) / (1/(4*phi^3))")
print(f"                   = 8/(a1*phi) = 8/{a1*phi:.4f} = {8/(a1*phi):.6f}")
print(f"  Exact ratio: {8*alpha/(alpha_s/2):.6f}")
print(f"  Leading order ratio: {8/(a1*phi):.6f}")

# So 8*alpha = (alpha_s/2) * 8/(a1*phi) * (1 + O(alpha))
# = 4*alpha_s/(a1*phi) * (1 + small correction)

# =====================================================================
# 12. THE GRAND SYNTHESIS
# =====================================================================
print(f"\n{'=' * 72}")
print("12. SYNTHESIS: WHAT MECHANISM GIVES -8*alpha?")
print(f"{'=' * 72}")

print(f"""
  FINDING 1: Tree-level m_H/m_W = phi (DERIVED from eigenvalue ratio)

  FINDING 2: 8*alpha ~ 2/(a1*phi^4) to leading order in alpha
    This equals 2*alpha/phi * phi^{-3}/a1... hmm, messy.

  FINDING 3: 8*alpha ~ alpha_s/2 * 8/(a1*phi)
    The ratio 8/(a1*phi) = {8/(a1*phi):.6f} ~ 1 but NOT exactly 1.

  FINDING 4: The exact correction is
    phi - 8*alpha = phi - 8*(root of 2*pi*x^2 - 4*a1*phi^4*x + 1)
    This involves pi (transcendental), so cannot be purely algebraic.

  FINDING 5: Leading-order formula
    phi - 8*alpha ~ (a1^2*phi^3 + 1) / (2*a1*phi^4)
    = (25*phi^3 + 1) / (10*phi^4)
    Accuracy: {abs((a1**2*phi**3+1)/(2*a1*phi**4) - (phi-8*alpha)):.2e} (excellent)

  MECHANISM CANDIDATES:

  (A) SPECTRAL: The correction -8*alpha = -Tr(D_F^2)*alpha arises from the
      fiber contribution to the spectral action. On the product space
      M x F, the spectral action receives a correction:
      delta S ~ alpha * Tr(D_F^2) * (Higgs sector terms)
      This would give m_H/m_W = phi - Tr(D_F^2)*alpha = phi - 8*alpha.
      STATUS: PLAUSIBLE but needs explicit computation.

  (B) RADIATIVE: The correction is a 1-loop electroweak correction
      on the 600-cell lattice. The factor 8 counts the number of
      independent gauge directions (rank E8).
      STATUS: POSSIBLE but 1-loop on lattice is non-standard.

  (C) ALGEBRAIC: The correction is a consequence of the alpha equation.
      The leading order phi - 2/(a1*phi^4) has a purely framework form,
      but the exact expression involves pi (from the alpha equation).
      STATUS: The formula IS framework-derived (alpha is framework),
      the question is whether the MECHANISM can be identified.
""")

# =====================================================================
# 13. NUMERICAL CHECK: m_H PREDICTIONS
# =====================================================================
print(f"{'=' * 72}")
print("13. NUMERICAL PREDICTIONS")
print(f"{'=' * 72}")

m_W_fw = 91.1876 * np.sqrt(20.0/26)  # m_Z * cos(tW) with exp m_Z
m_W_fw2 = 91.44 * np.sqrt(20.0/26)  # with m_Z from 16/15 formula

m_H_1 = m_W_exp * phi  # tree level with exp m_W
m_H_2 = m_W_exp * (phi - 8*alpha)  # with correction, exp m_W
m_H_3 = m_W_fw * (phi - 8*alpha)  # with correction, fw m_W
m_H_4 = m_W_fw2 * (phi - 8*alpha)  # with correction, fw m_W from 16/15

print(f"\n  m_H predictions:")
print(f"  Tree (phi, exp m_W):          {m_H_1:.3f} GeV ({(m_H_1/m_H_exp-1)*100:+.2f}%)")
print(f"  phi-8*alpha, exp m_W:         {m_H_2:.3f} GeV ({(m_H_2/m_H_exp-1)*100:+.2f}%)")
print(f"  phi-8*alpha, fw m_W (exp mZ): {m_H_3:.3f} GeV ({(m_H_3/m_H_exp-1)*100:+.2f}%)")
print(f"  phi-8*alpha, fw m_W (16/15):  {m_H_4:.3f} GeV ({(m_H_4/m_H_exp-1)*100:+.2f}%)")
print(f"  Experimental:                 {m_H_exp:.3f} GeV")

# =====================================================================
# SUMMARY
# =====================================================================
print(f"\n{'=' * 72}")
print("SUMMARY")
print(f"{'=' * 72}")

print(f"""
  STATUS: OPEN-8 PARTIALLY RESOLVED

  1. The correction -8*alpha is CONSISTENT with the framework:
     - 8 = Tr(D_F^2) = rank(E8) (DERIVED from McKay graph)
     - alpha is DERIVED from the alpha equation
     - So phi - 8*alpha uses only framework quantities

  2. The MECHANISM is not fully derived:
     - Most promising: spectral action on product space gives
       a correction proportional to alpha * Tr(D_F^2) to the
       Higgs mass ratio
     - The factor Tr(D_F^2) = 8 naturally emerges from the fiber
     - The coupling alpha naturally enters from the gauge sector

  3. Leading-order algebraic form:
     phi - 8*alpha ~ (a1^2*phi^3 + 1) / (2*a1*phi^4)
     This is purely algebraic in a1 and phi (no pi).

  4. Near-relation: 8*alpha ~ alpha_s/2 (within 1.1%)
     But this is NOT exact: ratio = 8/(a1*phi) = 0.989 != 1

  CATEGORY: IDENTIFIED -> PARTIALLY DERIVED
  The formula uses framework quantities.
  The mechanism (Tr(D_F^2)*alpha from spectral action) is PLAUSIBLE
  but the explicit spectral action computation has not been done.

  HONEST ASSESSMENT: This is between PATTERN and DERIVATION.
  We can identify 8 = Tr(D_F^2) and alpha = framework constant,
  but cannot yet prove that the spectral action produces exactly
  this combination as a correction to the tree-level phi.
""")

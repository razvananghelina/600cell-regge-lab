"""
exp475: Can the gauge prefactors (8/15, 1/3, 2/15) be derived from the
finite Dirac operator D_F on the McKay graph?

The gauge-channel skeleton 1+3+8 = 12 IS derived.
The question: can we get the SPECIFIC fractions from Tr over gauge sectors?
"""
import numpy as np
from itertools import product as iproduct

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120

print("="*70)
print("  exp475: GAUGE PREFACTORS FROM THE FINITE DIRAC OPERATOR")
print("="*70)

# ============================================================
# McKay graph of 2I (affine E8)
# ============================================================
# Node labeling (by spin j = 0, 1/2, 1, 3/2, 2, 5/2, 3, 7/2, branch):
# Standard physics labeling:
#   0: j=0   (dim 1, trivial)     WHITE
#   1: j=1/2 (dim 2, std)         BLACK
#   2: j=1   (dim 3, Sym^2)       WHITE   <-- 3 of A5, part of SU(3)
#   3: j=3/2 (dim 4)              BLACK
#   4: j=2   (dim 5)              WHITE   <-- 5 of A5, part of SU(3)
#   5: j=5/2 (dim 6)              BLACK
#   6: j=3   (dim 4)              WHITE   (note: NOT dim 4, this is the branch arm)
#   7: j=7/2 (dim 2)              BLACK
#   8: branch (dim 3)             WHITE   <-- 3' of A5 = ad(SU(2))

# Wait, I need to be precise. The 2I irreps and their dimensions:
# The binary icosahedral group 2I has 9 irreps with dimensions:
# 1, 2, 3, 4, 5, 6, 4, 2, 3
# These sum to 1+2+3+4+5+6+4+2+3 = 30

dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
n_irreps = len(dims)
print(f"\n  Irrep dimensions: {dims}")
print(f"  Total: {sum(dims)} (should be 30)")
print(f"  Sum of squares: {sum(d**2 for d in dims)} (should be |2I| = 120)")

# McKay graph adjacency (affine E8):
# Main chain: 0-1-2-3-4-5-6-7
# Branch: 4-8
adj = {
    0: [1],
    1: [0, 2],
    2: [1, 3],
    3: [2, 4],
    4: [3, 5, 8],  # branch point
    5: [4, 6],
    6: [5, 7],
    7: [6],
    8: [4],         # branch
}

print("\n  McKay graph (affine E8):")
for i in range(n_irreps):
    parity = "WHITE" if i % 2 == 0 or i == 8 else "BLACK"
    # Actually parity is by spin: integer spin = WHITE, half-integer = BLACK
    # j = 0, 1/2, 1, 3/2, 2, 5/2, 3, 7/2, branch(?)
    # For the branch node 8 (dim 3), it's at the end of a length-1 leg from node 4
    # In 2I: the branch rep has j that makes it WHITE
    # Actually: bipartite coloring of E8 tree: 0(W)-1(B)-2(W)-3(B)-4(W)-5(B)-6(W)-7(B), 8(B)
    # Wait: 4 connects to 8, and 4 is WHITE, so 8 must be BLACK.
    # But the paper says 3' (the branch node) is WHITE...
    # Let me re-check: the bipartite coloring is based on distance from root.
    # 0: distance 0 -> WHITE
    # 1: distance 1 -> BLACK
    # 2: distance 2 -> WHITE
    # 3: distance 3 -> BLACK
    # 4: distance 4 -> WHITE
    # 5: distance 5 -> BLACK
    # 6: distance 6 -> WHITE
    # 7: distance 7 -> BLACK
    # 8: distance 1 from 4 -> BLACK (since 4 is WHITE)

    if i in [0, 2, 4, 6]:  # even distance from root
        parity = "WHITE"
    else:
        parity = "BLACK"

    print(f"    rho_{i} (dim {dims[i]:d}, {parity}): neighbors = {adj[i]}")

# Chirality assignment (from paper):
# WHITE (integer spin, gamma_F = +1): rho_0, rho_2, rho_4, rho_6
# BLACK (half-integer, gamma_F = -1): rho_1, rho_3, rho_5, rho_7, rho_8
WHITE = [0, 2, 4, 6]
BLACK = [1, 3, 5, 7, 8]

print(f"\n  WHITE nodes: {WHITE} (dims {[dims[i] for i in WHITE]}, total {sum(dims[i] for i in WHITE)})")
print(f"  BLACK nodes: {BLACK} (dims {[dims[i] for i in BLACK]}, total {sum(dims[i] for i in BLACK)})")

# ============================================================
# Gauge group decomposition
# ============================================================
# The 12 nearest neighbors decompose under A5:
# 12 = 1 + 3 + 3' + 5
# where:
#   1 (trivial) -> U(1)
#   3' (branch dim-3 rep, rho_8) -> ad(SU(2))
#   3 + 5 (rho_2 + rho_4) -> ad(SU(3)) [by Burnside]

# So the gauge sectors in terms of McKay nodes:
U1_nodes = [0]          # rho_0, dim 1
SU2_nodes = [8]         # rho_8, dim 3 (the branch 3')
SU3_nodes = [2, 4]      # rho_2 (dim 3) + rho_4 (dim 5) = 8

print("\n  Gauge sector decomposition:")
print(f"    U(1):  nodes {U1_nodes}, dims {[dims[i] for i in U1_nodes]}, total = {sum(dims[i] for i in U1_nodes)}")
print(f"    SU(2): nodes {SU2_nodes}, dims {[dims[i] for i in SU2_nodes]}, total = {sum(dims[i] for i in SU2_nodes)}")
print(f"    SU(3): nodes {SU3_nodes}, dims {[dims[i] for i in SU3_nodes]}, total = {sum(dims[i] for i in SU3_nodes)}")
print(f"    Total: 1 + 3 + 8 = {1+3+8}")

# ============================================================
# Build the finite Dirac operator D_F
# ============================================================
print("\n" + "="*70)
print("  BUILDING D_F ON McKAY HILBERT SPACE (dim 30)")
print("="*70)

# H_F = direct sum of V_k (dim d_k for each irrep k)
# D_F has off-diagonal blocks Y_{kl} for each edge (k,l) of McKay graph
# |Y_{kl}|_F = 1/sqrt(2) universally

# Build block structure
offsets = [0]
for d in dims:
    offsets.append(offsets[-1] + d)
total_dim = offsets[-1]  # should be 30

print(f"  Hilbert space dim = {total_dim}")
print(f"  Block offsets: {offsets}")

# D_F is self-adjoint with Y_{kl} in the (k,l) block
# For simplicity, use Y_{kl} = (1/sqrt(2)) * (1/sqrt(min(dk,dl))) * I_{min x min}
# padded with zeros, adjusted so that ||Y||_F = 1/sqrt(2)

# Actually, the simplest choice preserving ||Y||_F = 1/sqrt(2):
# Y_{kl} is a dk x dl matrix with ||Y_{kl}||_F = 1/sqrt(2)
# The simplest: Y_{kl} has a single nonzero entry = 1/sqrt(2) if dk=dl=1
# For general dims: distribute uniformly

# For the Frobenius norm: ||Y||_F^2 = sum |y_ij|^2 = 1/2
# Simplest: Y_{kl} = (1/sqrt(2*dk*dl)) * J_{dk x dl} where J is all-ones
# Then ||Y||_F^2 = dk*dl * 1/(2*dk*dl) = 1/2 OK

# But this is a specific choice. The paper says the trace Tr(D_F^2) = 8
# is INDEPENDENT of the CG coefficients (only depends on ||Y||_F).
# So we can use ANY Y with the right norm.

D_F = np.zeros((total_dim, total_dim))

edges = []
for k in range(n_irreps):
    for l in adj[k]:
        if l > k:  # each edge once
            edges.append((k, l))

print(f"  Edges: {edges} ({len(edges)} edges)")

for k, l in edges:
    dk, dl = dims[k], dims[l]
    # Y_{kl}: dk x dl matrix with ||Y||_F = 1/sqrt(2)
    # Use uniform: each entry = 1/sqrt(2*dk*dl)
    Y = np.ones((dk, dl)) / np.sqrt(2 * dk * dl)

    # Place Y in D_F at (k,l) block and Y^T at (l,k) block
    ik, il = offsets[k], offsets[l]
    D_F[ik:ik+dk, il:il+dl] = Y
    D_F[il:il+dl, ik:ik+dk] = Y.T

# Verify basic properties
print(f"\n  D_F is symmetric: {np.allclose(D_F, D_F.T)}")
print(f"  Tr(D_F^2) = {np.trace(D_F @ D_F):.6f} (should be 8)")
print(f"  Tr(D_F^4) = {np.trace(np.linalg.matrix_power(D_F, 4)):.6f}")

TrD2 = np.trace(D_F @ D_F)
TrD4 = np.trace(np.linalg.matrix_power(D_F, 4))

# ============================================================
# Decompose traces by gauge sector
# ============================================================
print("\n" + "="*70)
print("  TRACE DECOMPOSITION BY GAUGE SECTOR")
print("="*70)

# Method 1: Block-diagonal restriction
# Tr_sector(D_F^2) = sum over edges (k,l) where k or l is in the sector
# Actually, Tr(D_F^2) = sum over edges 2*||Y_{kl}||_F^2 = 2 * |edges| * 1/2 = |edges| = 8

# For each edge, which gauge sector does it connect?
# An edge (k,l) connects two McKay nodes.
# The gauge "charge" flows through the edge.

# In NCG, the gauge field A_mu is an inner fluctuation:
# A = sum_i a_i [D_F, b_i]
# The gauge field lives in the TANGENT space of the gauge group,
# which is determined by the commutator structure of D_F.

# The trace of |F|^2 in the spectral action decomposes by gauge sector.
# For each sector G_i, the contribution is proportional to:
# Tr(pi_i(D_F^2)) where pi_i is the projection onto the G_i subspace.

# But wait - the "gauge sector" is about the VERTEX FIGURE (icosahedron),
# not about the McKay graph edges.

# Let me think differently.
# In the standard NCG model (Chamseddine-Connes-Marcolli):
# The gauge couplings are:
# 1/g_3^2 = f_0 * Tr(Q_3^2)   for SU(3)
# 1/g_2^2 = f_0 * Tr(Q_2^2)   for SU(2)
# 1/g_1^2 = f_0 * Tr(Q_1^2)   for U(1)

# where Q_i are the gauge generators acting on H_F,
# and the trace is over the FULL Hilbert space H_F.

# The RATIO of couplings at unification:
# g_3^2 / g_2^2 = Tr(Q_2^2) / Tr(Q_3^2)

# For the SM, this gives the GUT normalization condition:
# g_1^2 : g_2^2 : g_3^2 = 5/3 : 1 : 1 (at GUT scale)

# In our framework, the gauge generators are the generators of the
# isotropy group of the gauge decomposition 12 = 1 + 3 + 3' + 5.

# Method 2: Dimension-weighted contributions
# The fraction for each sector might simply be dim^2 / sum(dim^2)
# for the relevant irreps.

print("\n  Method 1: dim^2 / sum(dim^2) for gauge irreps")
dim_U1 = 1
dim_SU2 = 3   # ad(SU(2)) = 3' irrep
dim_SU3 = 8   # ad(SU(3)) = 3 + 5

total_dim_gauge = dim_U1**2 + dim_SU2**2 + dim_SU3**2
# But wait, dim_SU3 = 8, so dim_SU3^2 = 64. Total = 1+9+64 = 74. Not useful.

# Actually for the gauge KINETIC term, the coefficient is:
# proportional to dim(G_i) / dim(G_total) where G_total = the full gauge group

# dim(U(1)) = 1, dim(SU(2)) = 3, dim(SU(3)) = 8, total = 12
frac_U1 = 1/12
frac_SU2 = 3/12
frac_SU3 = 8/12
print(f"    U(1):  {frac_U1:.6f} = 1/12")
print(f"    SU(2): {frac_SU2:.6f} = 3/12 = 1/4")
print(f"    SU(3): {frac_SU3:.6f} = 8/12 = 2/3")
print(f"    Sum:   {frac_U1+frac_SU2+frac_SU3:.6f}")
print(f"    Target: 2/15, 1/3, 8/15")
print(f"    MISMATCH: 1/12 != 2/15, 1/4 != 1/3, 2/3 != 8/15")

# Method 2: Standard NCG trace over fermion representations
# In NCG, each fermion rep contributes to Tr(Q_i^2).
# For SU(3): quarks contribute C_2(3) = 4/3 per quark
# For SU(2): doublets contribute C_2(2) = 3/4 per doublet
# For U(1): each fermion contributes Y^2

# In our framework, the "fermion reps" are the 9 McKay irreps.
# But which irreps are "quarks" and which are "leptons"?

print("\n  Method 2: NCG-style Tr(Q_i^2) over McKay irreps")
# The 9 irreps of 2I, under the gauge decomposition:
# Each irrep rho_k, when tensored with the standard rep rho_1,
# gives a sum over neighbors: rho_1 x rho_k = sum_{l adj k} rho_l
# This IS the McKay correspondence.

# The SU(3) generators act on the {rho_2, rho_4} part of the spectrum.
# The SU(2) generators act on the {rho_8} part.
# The U(1) generator acts on {rho_0}.

# For each irrep rho_k, the "gauge charge" under sector G_i is
# determined by the decomposition of the tensor product.

# Method 3: Edge-based decomposition
# Each McKay edge connects two irreps. The "gauge content" of an edge
# is determined by which gauge sector it belongs to.
#
# In the 12 = 1+3+3'+5 decomposition of the vertex figure:
# The vertex has 12 neighbors, split as 1+3+3'+5
# Each irrep corresponds to a "direction" in gauge space

# The 8 McKay edges can be classified:
# Main chain edges: 0-1, 1-2, 2-3, 3-4, 4-5, 5-6, 6-7
# Branch edge: 4-8

print("\n  Method 3: Edge classification by gauge content")
# Each edge of the McKay graph corresponds to the tensor product
# rho_1 x rho_k -> rho_{neighbors of k}
# The EDGE (k,l) carries the representation rho_1 restricted to
# the (rho_k, rho_l) matrix element.

# In the Chamseddine-Connes framework:
# The inner fluctuation A = sum a_i [D_F, b_i] generates gauge fields
# The commutator [D_F, .] acts on the off-diagonal blocks (the edges)
# Each edge contributes to the gauge field strength

# For each edge (k,l), define its "multiplied dimension" dk*dl
for k, l in edges:
    print(f"    Edge ({k},{l}): dims ({dims[k]},{dims[l]}), product = {dims[k]*dims[l]}")

# Method 4: The NCG gauge coupling formula
# In the product geometry M x F:
# S_gauge = f_0 * sum_i (1/g_i^2) * int |F_i|^2
# where 1/g_i^2 = f_0 * sum_{fermions f} Tr_f(Q_i^2) * dim_f

# The trace over the finite space decomposes as:
# Tr(D_F^4) = sum over paths of length 2 on McKay graph
# Each path k-l-m contributes Y_{kl}^T Y_{kl} Y_{lm}^T Y_{lm}

print("\n" + "="*70)
print("  QUADRATIC CASIMIR APPROACH")
print("="*70)

# The gauge prefactors in NCG come from:
# c_i = Tr_F(Q_i^2) / sum_j Tr_F(Q_j^2)
# where Q_i are generators of gauge group G_i acting on H_F

# For our McKay-graph H_F (dim 30):
# - The SU(3) generators act on the 3+5 = 8-dim subspace (rho_2 + rho_4)
#   but also the quarks (other reps that transform under SU(3))
# - The SU(2) generators act on the 3'-dim subspace (rho_8)
#   but also the doublets

# Actually, in NCG the gauge generators act on ALL of H_F,
# not just the adjoint subspace. Each fermion has charges under all gauge groups.

# Let me compute the Casimir per EDGE:
# For each edge (k,l), the contribution to the gauge kinetic term of sector G_i is:
# proportional to the projection of that edge onto the G_i sector.

# The key identity from the paper:
# The vertex figure decomposes as 12 = 1(U1) + 3(SU2) + 8(SU3)
# So each of the 12 neighbors carries a specific gauge charge.
# In the Cayley graph, the 12 edges FROM each vertex are:
#   - 2 fiber edges (part of the Hopf S^1) -> U(1)
#   - 10 cross-fiber edges -> SU(2) x SU(3)

# But the McKay graph is a DIFFERENT object from the Cayley graph!
# The McKay graph has 9 nodes and 8 edges.
# The Cayley graph has 120 nodes and 720 edges.

# The connection: each McKay edge represents a TENSOR PRODUCT decomposition.
# Edge (k,l) means rho_1 x rho_k contains rho_l (and vice versa).

# Method 5: Direct dimension counting
print("\n  Method 5: Dimension-based fractions")
# In the 12 = 1+3+3'+5 decomposition:
# Each direction in the 12-dim space has a gauge assignment
# The spectral action traces over this 12-dim space
# Contribution of sector i = dim(sector_i) / 12

# But the PREFACTORS in the Lagrangian also involve the GROUP THEORY:
# L = (8/15)|G|^2 + (1/3)|W|^2 + (2/15)|B|^2
# These fractions arise from: N(rep) * C_2(rep) / normalization

# For the SM with N_gen = 3 generations:
# SU(3): 3 gens * (2 quarks * C_2(3) + ...)
# SU(2): 3 gens * (doublets * C_2(2) + ...)
# U(1): 3 gens * (sum Y^2 over all fermions)

# In our framework with N_gen = 3 and the McKay irreps:
# The fermions are the 9 irreps of 2I
# Each fermion has quantum numbers (color, isospin, hypercharge)
# determined by the gauge decomposition

# Standard NCG result for the SM:
# Tr(Q_3^2) = 2 * N_gen  (from quarks in fundamental of SU(3))
# Tr(Q_2^2) = N_gen * (N_c + 1) = 3*(3+1) = 12... varies by convention

# Let me just compute the RATIOS from first principles
# using the dimension data.

# The N_gen = 3 generations each have:
# - quarks: transform under SU(3) fundamental (dim 3)
# - leptons: SU(3) singlets
# - left doublets: SU(2) fundamental (dim 2)
# - right singlets: SU(2) singlets

# Quadratic Casimir per fermion:
# SU(3) fund: C_2 = 4/3
# SU(2) fund: C_2 = 3/4
# U(1): Y^2 (varies)

# The trace ratio that gives the prefactors is:
# In the spectral action: the gauge kinetic term coefficient for G_i is
# proportional to sum over fermions of dim(rep under G_i) * C_2(rep)

# For one generation:
# SU(3): 2 quarks (u_L, d_L doublet) + 2 quarks (u_R, d_R singlets) = 4 Weyl fermions
#   Each in fund of SU(3): Tr_f(T_a^2) = 1/2 per Weyl fermion
#   Total per gen: 4 * 1/2 = 2
# SU(2): 1 quark doublet (3 colors) + 1 lepton doublet (1 color) = 4 Weyl doublets
#   Each: Tr(T_i^2) = 1/2 per doublet
#   Total per gen: (3+1) * 1/2 = 2  (hmm, varies by normalization)
# U(1): sum Y^2 over all fermions

# Actually in the standard NCG normalization:
# The gauge kinetic terms are:
# S = f_0/(2*pi^2) * [a * int F_3^2 + b * int F_2^2 + c * int F_1^2]
# where:
# a = Tr(rho_3(1)^2) = sum over irreps of dim(color rep)^2 contributions
# ... this is getting complicated.

# Let me just check: do the TARGET fractions 8/15, 1/3, 2/15 arise from
# simple dimension counting?

print("\n  TARGET fractions: 8/15, 5/15, 2/15")
print(f"  Ratio: 8 : 5 : 2 = {8}:{5}:{2}")
print(f"  Sum: 15 = 3*a1 = a1*N_gen")
print()

# What gives 8:5:2?
# 8 = rank(E8) = dim(SU(3))
# 5 = a1 = dim of irrep 5 of A5
# 2 = first nontrivial Fibonacci
# Sum = 15 = 3*5

# Alternative: 8:5:2 where 5 = 3+2?
# Or: 8 = 2*4 (quarks * colors?), 5 = ?, 2 = leptons?

# What if the prefactor for sector i is:
# dim(G_i) * (dim(G_i) + something) / normalization?

# dim(SU3)*something_3 : dim(SU2)*something_2 : dim(U1)*something_1
# = 8*x : 3*y : 1*z with 8x+3y+z = 15

# If x=y=z=1: 8+3+1 = 12 != 15
# If x=1, y=5/3, z=2: 8+5+2 = 15 OK
#   So SU(3): 8*1=8, SU(2): 3*(5/3)=5, U(1): 1*2=2

# The SU(2) factor 5/3 = a1/N_gen. Interesting!
# The U(1) factor 2.

# What if: prefactor_i = dim(G_i) * n_i where n_i is the number of
# "charged" sectors?
# SU(3): 8 * 1 = 8 (one adjoint)
# SU(2): 3 * (5/3) = 5
# U(1): 1 * 2 = 2

# The 5/3 for SU(2) could be: a1/N_gen, or (dim of 5-irrep)/N_gen

# OR: could it be that:
# prefactor_i = Tr(C_2(G_i)) over all 2I irreps?
# Let's compute!

print("  Casimir traces over 2I irreps:")

# For each 2I irrep rho_k (dim d_k), its decomposition under the gauge group
# determines its charges. But the gauge group SU(3)xSU(2)xU(1) comes from
# the A_5 decomposition of the vertex figure, not from 2I directly.

# The relationship is through the McKay correspondence:
# rho_1 (std, dim 2) tensoring with each rho_k gives the neighbors on McKay graph.
# The gauge decomposition 12 = 1+3+3'+5 is about A_5, not 2I.

# Let me try a DIFFERENT approach entirely.
# The prefactors might come from the ADJACENCY SPECTRUM of the McKay graph.

# The adjacency matrix of affine E8 has eigenvalues:
# 2*cos(k*pi/30) for specific k values... actually it's more subtle.

# For a tree graph, Tr(A^2) = 2*|edges| = 16
# For individual nodes: (A^2)_{kk} = degree(k)

A = np.zeros((9, 9))
for k in range(9):
    for l in adj[k]:
        A[k][l] = 1

A2 = A @ A
A4 = A2 @ A2

print(f"\n  McKay adjacency matrix A (9x9):")
print(f"  Tr(A^2) = {np.trace(A2):.0f} (= 2*|edges| = 16)")
print(f"  Tr(A^4) = {np.trace(A4):.0f}")
print(f"  Eigenvalues of A: {sorted(np.linalg.eigvalsh(A), reverse=True)}")

# Per-node (A^2)_{kk} = degree:
print(f"\n  Per-node degree (A^2)_kk:")
for k in range(9):
    print(f"    rho_{k} (dim {dims[k]}): degree = {A2[k,k]:.0f}")

# Weighted by dim^2:
print(f"\n  dim^2-weighted contributions:")
for k in range(9):
    print(f"    rho_{k}: dim^2 = {dims[k]**2}, degree = {A2[k,k]:.0f}, weighted = {dims[k]**2 * A2[k,k]:.0f}")

# ============================================================
# THE KEY COMPUTATION: dim-weighted edge contributions
# ============================================================
print("\n" + "="*70)
print("  EDGE CONTRIBUTIONS WEIGHTED BY DIMENSIONS")
print("="*70)

# For each edge (k,l), the contribution to the trace is:
# In Tr(D_F^2): each edge contributes 2*||Y||_F^2 = 1 (since ||Y||=1/sqrt(2))
# But the GAUGE decomposition depends on which nodes the edge connects.

# The gauge content of an edge (k,l):
# In the tensor product rho_1 x rho_k -> sum of rho_l's,
# the "coupling" to the gauge group is determined by
# the A_5 content of the irreps.

# Under A_5 (icosahedral rotation group):
# rho_0 (dim 1) -> trivial
# rho_1 (dim 2) -> standard 2-dim rep (does not exist in A_5, but in 2I)
# rho_2 (dim 3) -> the 3 of A_5 (part of SU(3))
# rho_3 (dim 4) -> not in A_5 (half-integer spin)
# rho_4 (dim 5) -> the 5 of A_5 (part of SU(3))
# rho_5 (dim 6) -> not in A_5
# rho_6 (dim 4) -> (integer spin, but paired)
# rho_7 (dim 2) -> not in A_5
# rho_8 (dim 3) -> the 3' of A_5 (= ad(SU(2)))

# The gauge-RELEVANT nodes are: 0 (U(1)), 2+4 (SU(3)), 8 (SU(2))
# Total gauge dimension: 1 + 8 + 3 = 12

# The "gauge fraction" for each edge could be:
# f_edge = dim(gauge-relevant node) / 12

# For each edge, identify the gauge content:
print("  Edge gauge assignments:")
total_gauge_weight = 0
gauge_weights = {'U1': 0, 'SU2': 0, 'SU3': 0, 'other': 0}

for k, l in edges:
    # Which gauge sectors do k and l belong to?
    sectors_k = []
    sectors_l = []
    for node, sectors in [(k, sectors_k), (l, sectors_l)]:
        if node == 0:
            sectors.append('U1')
        elif node == 8:
            sectors.append('SU2')
        elif node in [2, 4]:
            sectors.append('SU3')
        else:
            sectors.append('matter')

    print(f"    ({k},{l}): {sectors_k[0]:>6s} -- {sectors_l[0]:>6s}  (dims {dims[k]}x{dims[l]})")

# ============================================================
# Method 6: The Connes-Lott formula
# ============================================================
print("\n" + "="*70)
print("  CONNES-LOTT: Tr over fermion reps")
print("="*70)

# In NCG, the gauge field strength comes from:
# F = pi(dA + A^2) where A is the inner fluctuation
# The inner fluctuation for gauge group G_i has:
# A_i = sum_a gamma^mu A_mu^a otimes Q_i^a
# where Q_i^a are generators of G_i acting on H_F

# The spectral action gives:
# S_i = f_0/(2pi^2) * Tr_F(Q_i^a Q_i^a) * int |F_i|^2

# So the coupling ratios are:
# (1/g_i^2) / (1/g_j^2) = Tr_F(Q_i^2) / Tr_F(Q_j^2)

# For the 600-cell framework:
# H_F = sum_{k=0}^{8} V_k with dim V_k = d_k
# The gauge generators Q_i act on each V_k with some representation

# For SU(3): Q_3^a acts on V_2 (dim 3) as the FUNDAMENTAL and on V_4 (dim 5)
# as part of the adjoint. On other V_k, it acts trivially or via tensor products.

# For SU(2): Q_2^i acts on V_8 (dim 3') as the ADJOINT.
# On V_1 (dim 2) as the FUNDAMENTAL, on V_3 (dim 4)... etc.

# For U(1): Q_1 is a diagonal matrix with hypercharges.

# This is getting into representation theory of 2I under A_5.
# Let me compute the EMBEDDING of A_5 irreps into 2I irreps.

# Under A_5 (the quotient 2I/Z_2):
# A_5 has 5 irreps: 1, 3, 3', 4, 5
# 2I has 9 irreps: 1, 2, 3, 4, 5, 6, 4', 2', 3'
# The A_5 irreps embed as:
# 1_A5 = rho_0
# 3_A5 = rho_2
# 3'_A5 = rho_8  (or rho_4? need to check)
# 4_A5 = rho_6  (dim 4, integer spin)
# 5_A5 = rho_4

# The half-integer spin reps (rho_1, rho_3, rho_5, rho_7) are SPINORS
# and don't appear in A_5.

# So under A_5:
# Gauge group: SU(3) x SU(2) x U(1) from 12 = 1+3+3'+5
# SU(3) adjoint = 3+5 (from rho_2 + rho_4 of 2I)
# SU(2) adjoint = 3' (from rho_8 of 2I)
# U(1) = 1 (from rho_0 of 2I)

# Now, each 2I irrep rho_k decomposes under SU(3)xSU(2)xU(1):
# But this depends on how the gauge group embeds.

# The KEY is: what are the CHARGES of each McKay node under the gauge group?

# rho_0 (dim 1): gauge singlet (it IS the U(1) direction)
# rho_2 (dim 3): transforms as 3 under SU(3) (it IS part of SU(3))
# rho_4 (dim 5): transforms as 5 under A_5 -> ad(SU(3)) complement
# rho_6 (dim 4): transforms as 4 under A_5
# rho_8 (dim 3'): transforms as ad(SU(2))

# For the GAUGE COUPLING computation, we need Tr_F(Q_i^2) where
# the trace is over the FERMIONIC Hilbert space.
# In our framework, the fermionic content is the 9 McKay irreps.

# For Tr(Q_3^2) (SU(3)):
# Only reps that transform nontrivially under SU(3) contribute.
# rho_2 (dim 3): fundamental of SU(3), C_2(3) = 4/3, contribution = 3*(4/3) = 4
# rho_4 (dim 5): the 5 of A_5, which under SU(3) is... ?
# In the decomposition 8 = 3+5, the 5-dim irrep of A_5 is the SECOND
# part of the SU(3) adjoint. So rho_4 transforms as 5_A5 under A_5.
# Under SU(3), this 5_A5 is NOT a standard SU(3) rep. It's part of the adjoint.
# Actually ad(SU(3)) = 8_SU3 which under A_5 restricts to 3_A5 + 5_A5.
# So the 5_A5 is the complement of 3_A5 in 8_SU3.

# This is subtle. Let me just compute the Casimir ratios from the
# standard NCG result and see if they match 8:5:2.

print("  Standard NCG for SM (per generation):")
# For N_gen = 3 generations:
# Each gen has: (Q_L, u_R, d_R, L_L, e_R, nu_R)
# where Q_L = (u_L, d_L) in (3, 2, 1/6), etc.

# Tr(T_3^a T_3^a) per generation (SU(3)):
# Q_L: 3 x 2 Weyl, each in fund(3): 2 * (1/2) = 1 per color doublet
#   Total: 2 * 1/2 = 1 (per quark doublet in 3)
# u_R, d_R: each in fund(3): 1/2 each
# Total SU(3) per gen: 2*(1/2) + 1/2 + 1/2 = 2

# Tr(T_2^i T_2^i) per generation (SU(2)):
# Q_L: (3 colors) x fund(2): 3 * (1/2) = 3/2
# L_L: fund(2): 1/2
# Total SU(2) per gen: 3/2 + 1/2 = 2

# Tr(Y^2) per generation (U(1), GUT normalized with sqrt(5/3)):
# Q_L: 6 * (1/6)^2 = 6/36 = 1/6
# u_R: 3 * (2/3)^2 = 3*4/9 = 4/3
# d_R: 3 * (1/3)^2 = 3/9 = 1/3
# L_L: 2 * (1/2)^2 = 1/2
# e_R: 1 * 1^2 = 1
# nu_R: 0
# Total Y^2 per gen: 1/6+4/3+1/3+1/2+1 = (1+8+2+3+6)/6 = 20/6 = 10/3
# With GUT normalization (5/3): Tr = (5/3)*(10/3) = 50/9... hmm

# Standard result: at GUT unification g_1 = g_2 = g_3 with GUT normalization
# This gives equal coupling coefficients, NOT 8:5:2.

# The 8:5:2 ratio is NOT the GUT normalization.
# It's something specific to the 600-cell framework.

# Let me check: 8/15 + 5/15 + 2/15 = 1
# The total sums to 1 (normalized).
# In the spectral action: S_YM = c_2 * f_0 * [a*|G|^2 + b*|W|^2 + c*|B|^2]
# with a+b+c = 1 (if normalized to c_2*f_0*Tr(F^2))

# What determines a:b:c?
# In the Chamseddine-Connes model:
# a = Tr(rho^2(G_3)) / Tr_total
# etc.

# For our model, let me try the SIMPLEST possibility:
# The fraction for sector i is proportional to:
# sum over fermions f of: dim(rep of f under G_i)^2 / dim(G_i)

# Or simply: fraction_i = dim(ad(G_i)) * (something per generation)

# ============================================================
# Method 7: Direct from A_5 representation theory
# ============================================================
print("\n" + "="*70)
print("  A_5 REPRESENTATION THEORY")
print("="*70)

# The vertex figure is the icosahedron = Cayley graph of A_5
# Its 12 vertices decompose under A_5 as 12 = 1 + 3 + 3' + 5
# The EDGE decomposition: each vertex has 5 neighbors, total 30 directed edges
# These edges decompose under A_5 as well.

# The natural inner product on the 12-dim rep:
# <rho_i | rho_j> = delta_ij * dim(rho_i)
# So Tr(P_i) = dim(rho_i) for projection P_i onto sector i

# The gauge kinetic term involves Tr(F^2) which has TWO powers of the connection.
# The coefficient for sector i is: Tr(P_i^2) / Tr(I) = dim(rho_i) / 12

# But that gives 1:3:3:5 (for 1, 3, 3', 5), not 2:5:8.

# Wait - 3' and 3 are in DIFFERENT sectors:
# 3' -> SU(2), 3+5 -> SU(3), 1 -> U(1)
# So: U(1):SU(2):SU(3) = 1:3:8 for the dimensions.
# The fractions 1/12:3/12:8/12 = 1/12:1/4:2/3

# But the TARGET is 2/15:1/3:8/15 = 2:5:8 (times 1/15)

# The difference: 2:5:8 vs 1:3:8
# The SU(3) part is the SAME (8)!
# The U(1) and SU(2) parts differ: 1->2 and 3->5

# Ratio 2:5 vs 1:3 -> factor 2/1=2 for U(1), 5/3 for SU(2)
# Why would U(1) get doubled and SU(2) get multiplied by 5/3?

# In the standard NCG model, the GUT normalization gives:
# g_1^2/g_3^2 = (5/3) at the GUT scale
# g_2^2/g_3^2 = 1

# With this normalization: 8:(5/3*3):((5/3)*1) = 8:5:(5/3)
# Hmm, not quite 8:5:2 either.

# Let me check: if g_1 has GUT normalization sqrt(5/3):
# Then the coefficient for U(1) in the Lagrangian is:
# (g_1^2/g_GUT^2) * dim(U(1)) = (5/3) * 1 = 5/3
# No...

# Actually the standard NCG normalization gives:
# 1/g_1^2 : 1/g_2^2 : 1/g_3^2 = sum Y^2 : N_gen+1 : N_gen
# But this depends on the specific fermion content.

# Let me just CHECK whether 8:5:2 = 8 : a1 : 2 = 8 : 5 : 2
# and whether 15 = a1 * N_gen
print(f"  8:5:2 = rank(E8) : a1 : 2")
print(f"  Sum = 15 = a1 * N_gen = {a1*3}")
print(f"  8/15 = rank(E8)/(a1*N_gen) = {8/(a1*3):.10f}")
print(f"  5/15 = a1/(a1*N_gen) = 1/N_gen = {1/3:.10f}")
print(f"  2/15 = 2/(a1*N_gen) = {2/(a1*3):.10f}")

print(f"\n  Alternative: 8 = dim(SU(3)), 5 = dim(A5 irrep), 2 = dim(std)")
print(f"  Or: 8 = 2*rank, 5 = a1, 2 = half of (a1-1)")

# ============================================================
# FINAL: Check if 8:5:2 comes from dim^2 normalization
# ============================================================
print("\n" + "="*70)
print("  SQUARED-DIMENSION NORMALIZATION?")
print("="*70)

# What if: fraction_i = dim(sector_i)^2 / sum dim^2 but with different sector assignment?
# Standard: SU(3)=8, SU(2)=3, U(1)=1
# dim^2: 64:9:1, sum=74. Fractions: 64/74, 9/74, 1/74. NOT 8:5:2.

# What if we use the A_5 irrep dims instead?
# A_5 has: 1, 3, 3', 4, 5
# Gauge sectors: U(1)=1, SU(2)=3', SU(3)=3+5
# dim^2: 1^2=1, 3^2=9, 3^2+5^2=34
# Fractions: 1/44, 9/44, 34/44. NOT useful.

# What about dim * (dim+1) / 2?
# SU(3): 8*9/2 = 36
# SU(2): 3*4/2 = 6
# U(1): 1*2/2 = 1
# Total: 43. NOT useful.

# Or: use C_2 * dim?
# SU(N): C_2(adj) = N, dim(adj) = N^2-1
# SU(3): C_2=3, dim=8 -> 24
# SU(2): C_2=2, dim=3 -> 6
# U(1): C_2=0, dim=1 -> 0 (fails for U(1))

# What about N_gen-dependent?
N_gen = 3
# For each sector: n_fermions * C_2(fund)
# SU(3): N_gen * 2 quarks * C_2(3) = 3*2*(4/3) = 8 !!
# SU(2): N_gen * (N_c+1) * C_2(2) = 3*(3+1)*(3/4) = 3*4*3/4 = 9... no

# Wait: SU(3): N_gen * 2 * C_2(fund of SU(3)) = 3 * 2 * (4/3) = 8 OKOKOK
# SU(2): N_gen * (3+1) * C_2(fund of SU(2)) = 3 * 4 * (3/4) = 9... not 5

# Hmm. But 8 = N_gen * 2 * C_2(3) is interesting!
print("  N_gen * n_quarks * C_2(fund_SU3) = 3*2*(4/3) = 8  OK")

# For SU(2), what gives 5?
# N_gen * n_doublets * C_2(fund_SU2) = 3 * n * (3/4) = 5
# => n = 20/9... not integer

# Or: a1 = 5 directly?
# In the framework: the SU(2) coefficient is just a1.
# Is there a reason for this?

# The DUAL Coxeter number of SU(2) is 2.
# dim(SU(2)) * (dual Coxeter - 1) = 3 * 1 = 3... not 5
# dim(SU(2)) + dual Coxeter = 3 + 2 = 5 = a1!

print(f"\n  SU(2) coefficient: dim(SU(2)) + h_dual(SU(2)) = 3 + 2 = {3+2} = a1  ?!")
print(f"  SU(3) coefficient: dim(SU(3)) + h_dual(SU(3))... = 8 + 3 = 11 != 8")
print(f"  That doesn't work for SU(3).")

# ============================================================
print("\n" + "="*70)
print("  ASSESSMENT")
print("="*70)

print("""
  RESULT: The gauge prefactors 8:5:2 (summing to 15 = a1*N_gen) have
  the following partial decomposition:

  - SU(3): 8 = rank(E8) = N_gen * 2 * C_2(fund_3) = 3 * 2 * (4/3)
    This can be DERIVED from the fermion content (3 generations,
    2 quark chiralities, SU(3) Casimir 4/3).

  - SU(2): 5 = a1. The origin is UNCLEAR.
    Candidate: dim(SU(2)) + h_dual(SU(2)) = 3+2 = 5, but this
    doesn't generalize to SU(3).

  - U(1): 2 = dim(std rep of 2I). The origin is UNCLEAR.
    Candidate: related to the 2 Hopf fiber neighbors per vertex.

  - Total: 15 = a1 * N_gen = 5 * 3. STRUCTURAL.

  STATUS: The prefactors are NOT derivable from D_F traces alone.
  The D_F traces are CG-independent (depend only on ||Y||_F), so they
  give no gauge-sector resolution. The decomposition requires the
  A_5 representation theory of the vertex figure, which determines
  the gauge-channel skeleton (1+3+8) but NOT the specific fractions.

  The fractions 8:5:2 are a PATTERN with interesting algebraic structure
  (8=rank(E8), 5=a1, 2=dim(std), sum=15=a1*N_gen) but no derivation
  from the spectral action or D_F has been found.

  HONEST CONCLUSION: NEGATIVE. Prefactors remain PATTERN.
""")

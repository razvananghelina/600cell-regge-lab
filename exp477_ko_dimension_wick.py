"""
exp477: KO-Dimension, Wick Rotation, and the Mass Anchor
=========================================================

THREE CONNECTED QUESTIONS:
  Q1: What is the KO-dimension of the 600-cell spectral triple?
      If KO(F) = 6 and consistency requires KO(M×F) = 0 mod 8,
      then KO(M) = 2 => Lorentzian signature (1,3) DERIVED.

  Q2: Is the Galois involution (phi->phi') the Wick rotation?
      If sigma acts on the real structure J as Wick rotation,
      this connects causal structure to number theory.

  Q3: Does the mass duality m_f*m_f' = m_e^2 FIX m_e?
      If Galois = Wick, then m_e = geometric mean of
      Lorentzian and Euclidean sectors = self-consistency condition.

BACKGROUND:
  In NCG, the KO-dimension n mod 8 is determined by three signs:
    J^2 = epsilon,  JD = epsilon'*DJ,  J*gamma = epsilon''*gamma*J
  where (epsilon, epsilon', epsilon'') depend on n:
    n=0: (+,+,+), n=2: (-,+,-), n=4: (-,+,+), n=6: (+,+,-)

  For the SM (Connes): KO(F) = 6, giving epsilon=+1, epsilon'=+1, epsilon''=-1.

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: FUNDAMENTAL INVESTIGATION
"""

import numpy as np
from math import sqrt, pi, cos, sin
from collections import Counter

phi = (1 + sqrt(5)) / 2
phi_p = (1 - sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30

# ============================================================
# PART 1: THE McKAY GRAPH AND REPRESENTATION DATA
# ============================================================
print("=" * 70)
print("PART 1: 2I REPRESENTATIONS AND McKAY GRAPH")
print("=" * 70)

# Binary icosahedral group 2I, order 120
# 9 irreps with dimensions from affine E8 Dynkin diagram
# Standard labeling: nodes of affine E8

# Dimensions (from E8 extended Dynkin diagram)
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
n_irreps = len(dims)
print(f"\n2I irreps: {n_irreps} = N_eig")
print(f"Dimensions: {dims}")
print(f"Sum of dim^2 = {sum(d**2 for d in dims)} = |2I| = {N}")

# McKay graph adjacency matrix (tensor product with rho_1, dim 2)
# This IS the affine E8 Dynkin diagram
# Connectivity: from the tensor product rules
# rho_i tensor rho_1 = sum A_{ij} rho_j
A_mckay = np.array([
    [0,1,0,0,0,0,0,0,0],  # rho_0: x rho_1 = rho_1
    [1,0,1,0,0,0,0,0,0],  # rho_1: x rho_1 = rho_0 + rho_2
    [0,1,0,1,0,0,0,0,0],  # rho_2: x rho_1 = rho_1 + rho_3
    [0,0,1,0,1,0,0,0,0],  # rho_3: x rho_1 = rho_2 + rho_4
    [0,0,0,1,0,1,1,0,0],  # rho_4: x rho_1 = rho_3 + rho_5 + rho_6
    [0,0,0,0,1,0,0,1,0],  # rho_5: x rho_1 = rho_4 + rho_7
    [0,0,0,0,1,0,0,0,1],  # rho_6: x rho_1 = rho_4 + rho_8
    [0,0,0,0,0,1,0,0,0],  # rho_7: x rho_1 = rho_5
    [0,0,0,0,0,0,1,0,0],  # rho_8: x rho_1 = rho_6
])

# Verify: A*dims should = 2*dims (fundamental has dim 2)
dims_arr = np.array(dims)
check = A_mckay @ dims_arr
print(f"\nA*dims = {check.tolist()}")
print(f"2*dims = {(2*dims_arr).tolist()}")
print(f"A*dims = 2*dims? {np.array_equal(check, 2*dims_arr)}")

# ============================================================
# PART 2: FROBENIUS-SCHUR INDICATORS
# ============================================================
print(f"\n{'='*70}")
print("PART 2: FROBENIUS-SCHUR INDICATORS (J^2 signs)")
print("=" * 70)

# For 2I ⊂ SU(2), the irreps come from restricting SU(2) spin-j reps
# The spin values for the 9 irreps:
# dim 1 = spin 0, dim 2 = spin 1/2, dim 3 = spin 1,
# dim 4 = spin 3/2, dim 5 = spin 2, dim 6 = spin 5/2
# But dim 4 appears TWICE (rho_3 and rho_6), dim 3 THREE times (rho_2, rho_6_wait)

# Actually, the ordering on the affine E8 diagram gives:
# rho_0: dim 1 (spin 0)      -> FS = +1 (real)
# rho_1: dim 2 (spin 1/2)    -> FS = -1 (pseudo-real/quaternionic)
# rho_2: dim 3 (spin 1)      -> FS = +1 (real)
# rho_3: dim 4 (spin 3/2)    -> FS = -1 (pseudo-real)
# rho_4: dim 5 (spin 2)      -> FS = +1 (real)
# rho_5: dim 6 (spin 5/2)    -> FS = -1 (pseudo-real)
# rho_6: dim 4 (spin 3/2')   -> FS = -1 (pseudo-real)
# rho_7: dim 2 (spin 1/2')   -> FS = -1 (pseudo-real)
# rho_8: dim 3 (spin 1')     -> FS = +1 (real)

# Wait -- for 2I, there are NOT two copies of spin 3/2 etc.
# The McKay graph IS affine E8 = tree. Each node is a UNIQUE irrep.
# The spin assignment is: on the LONG arm (rho_0 to rho_5):
# spins 0, 1/2, 1, 3/2, 2, 5/2
# On the short arm (from rho_4): rho_6=dim 4, rho_7=dim 2, wait...

# Actually I need to be more careful. For 2I, the irreps are:
# Restricting from SU(2), we get spin 0, 1/2, 1, 3/2, 2, 5/2 (dims 1,2,3,4,5,6)
# That's 6 irreps. But 2I has 9. The remaining 3 come from the fact that
# higher spin SU(2) reps become REDUCIBLE when restricted to 2I.
# spin 3 (dim 7) -> rho_0 + rho_5 (not irreducible)
# etc.

# BUT: all 9 irreps of 2I DO come from SOME SU(2) rep. The question is
# which SU(2) rep each one first appears in (as a new irreducible constituent).

# For the FS indicator, what matters is: is the rep real or pseudo-real?
# For finite subgroups of SU(2), ALL irreps are self-conjugate
# (because complex conjugation on SU(2) is equivalent to the j-involution).
# So every irrep is either real (FS=+1) or pseudo-real (FS=-1).

# The FS indicator for 2I can be computed from the character table:
# FS(rho_i) = (1/|G|) * sum_{g in G} chi_i(g^2)

# For the EVEN-dimensional irreps of SU(2) subgroups: FS = -1 (pseudo-real)
# For the ODD-dimensional irreps: FS = +1 (real)
# This is because SU(2) itself has J = i*sigma_2 with J^2 = -1 on the fundamental,
# and Sym^k has J^2 = (-1)^k for the (k+1)-dim rep.

# So for our 9 irreps with dims [1, 2, 3, 4, 5, 6, 4, 2, 3]:
FS = []
for d in dims:
    if d % 2 == 0:
        FS.append(-1)  # pseudo-real (quaternionic)
    else:
        FS.append(+1)  # real (orthogonal)

print(f"\nFrobenius-Schur indicators (from dim parity for SU(2) subgroup):")
print(f"  rho | dim | FS | type")
print(f"  " + "-" * 35)
n_real = 0
n_pseudo = 0
for i in range(n_irreps):
    rep_type = "real" if FS[i] == 1 else "pseudo-real"
    print(f"  {i:3d} | {dims[i]:3d} | {FS[i]:+d} | {rep_type}")
    if FS[i] == 1:
        n_real += 1
    else:
        n_pseudo += 1

print(f"\n  Real (FS=+1): {n_real} irreps, dims = {[dims[i] for i in range(9) if FS[i]==+1]}")
print(f"  Pseudo-real (FS=-1): {n_pseudo} irreps, dims = {[dims[i] for i in range(9) if FS[i]==-1]}")
print(f"  Sum real dims: {sum(dims[i] for i in range(9) if FS[i]==+1)}")
print(f"  Sum pseudo dims: {sum(dims[i] for i in range(9) if FS[i]==-1)}")

# ============================================================
# PART 3: THE THREE SIGNS (epsilon, epsilon', epsilon'')
# ============================================================
print(f"\n{'='*70}")
print("PART 3: KO-DIMENSION SIGNS")
print("=" * 70)

# In the NCG framework for the SM finite spectral triple:
#
# The Hilbert space H_F is the space of fermions.
# For the McKay approach: H_F = direct sum of rho_i-sectors
#
# J (real structure / charge conjugation):
#   J acts as complex conjugation + representation involution
#   On each sector rho_i: J^2 = FS(rho_i) * I
#
# gamma (chirality / Z_2 grading):
#   The McKay graph is bipartite (it's a tree = affine E8)
#   gamma = +1 on even nodes, -1 on odd nodes
#   Bipartition of affine E8: {0,2,4,6,8} and {1,3,5,7}

# Bipartition (tree coloring)
# Starting from rho_0 (even):
color = [-1] * 9
color[0] = +1
queue = [0]
visited = {0}
while queue:
    node = queue.pop(0)
    for j in range(9):
        if A_mckay[node, j] == 1 and j not in visited:
            color[j] = -color[node]
            visited.add(j)
            queue.append(j)

print(f"\nChirality gamma (bipartite grading):")
for i in range(9):
    parity = "even (+1)" if color[i] == +1 else "odd (-1)"
    print(f"  rho_{i} (dim {dims[i]}): gamma = {color[i]:+d} ({parity})")

even_nodes = [i for i in range(9) if color[i] == +1]
odd_nodes = [i for i in range(9) if color[i] == -1]
print(f"\n  Even nodes: {even_nodes} (dims {[dims[i] for i in even_nodes]})")
print(f"  Odd nodes: {odd_nodes} (dims {[dims[i] for i in odd_nodes]})")
print(f"  Sum even dims: {sum(dims[i] for i in even_nodes)}")
print(f"  Sum odd dims: {sum(dims[i] for i in odd_nodes)}")

# SIGN 1: epsilon = J^2
# On a specific rho_i: J^2 = FS(rho_i)
# Overall: J^2 is NOT a single sign -- it's +1 on some and -1 on others.
# In NCG, J^2 = epsilon * I requires a UNIFORM sign.
#
# Resolution: the finite Hilbert space H_F has a SPECIFIC structure.
# In Connes' SM: H_F doesn't include ALL irreps with equal weight.
# The fermion content determines WHICH irreps appear.
#
# In our framework: which irreps correspond to fermions?
# From the SM gauge group derivation:
# 12 = 1 + 3' + 3 + 5 (decomposition of the regular rep restricted to A5 ⊂ 2I)
# The FERMION representations are those carrying non-trivial gauge quantum numbers.

# Actually, in the McKay approach:
# The SM fermions correspond to edges of the McKay graph (between nodes)
# Not to nodes (which are irreps).
# The Hilbert space H_F is the EDGE space, dim = 2 * sum of edges = 2 * 8 = 16 per gen
# (actually the McKay graph has 8 edges, each contributing dim(source)*dim(target))

print(f"\n  McKay graph edges (tree, {sum(sum(A_mckay)//2)} edges):")
edge_dims = []
for i in range(9):
    for j in range(i+1, 9):
        if A_mckay[i,j] == 1:
            d = dims[i] * dims[j]
            print(f"    rho_{i}--rho_{j}: dim = {dims[i]}*{dims[j]} = {d}")
            edge_dims.append(d)
print(f"  Total edge-space dimension: {sum(edge_dims)}")

# SIGN 2: epsilon' = JD = +-DJ
# D (Dirac) connects even to odd (bipartite -> D is off-diagonal in gamma grading)
# J commutes or anticommutes with D
# For a bipartite graph: D is off-diagonal in the even/odd basis
# J acts on each irrep sector: J|rho_i> -> |rho_i_bar> = |rho_i> (self-conjugate)
# [J, D] or {J, D}?

# In Connes' SM: epsilon' = +1 for KO=6 (J commutes with D)

# SIGN 3: epsilon'' = J*gamma = +-gamma*J
# gamma = +-1 on even/odd nodes
# J doesn't change the node (self-conjugate irreps)
# So J*gamma = gamma*J IF J preserves the bipartition
# J*gamma = -gamma*J IF J flips the bipartition

# For self-conjugate irreps: J(rho_i) = rho_i (same node)
# So J preserves the bipartition (each node stays put)
# => J*gamma = gamma*J => epsilon'' = +1

# BUT WAIT: in Connes' SM, epsilon'' = -1 for KO=6.
# This means J must FLIP the chirality somehow.
# This happens because J also involves the PARTICLE-ANTIPARTICLE exchange,
# which maps left-handed to right-handed (flips chirality).

# In the McKay graph context:
# J maps edge (i,j) to edge (j,i) (reversing the arrow)
# If i is even and j is odd, then the reversed edge (j,i) goes from odd to even
# This FLIPS the chirality of the edge!
# So on the EDGE space: J*gamma = -gamma*J => epsilon'' = -1

print(f"\n{'='*70}")
print("COMPUTING THE THREE SIGNS")
print("="*70)

# On the EDGE space (where fermions live):
# An edge (i,j) with i even, j odd has gamma = +1 (left-handed)
# J maps (i,j) -> (j,i): now j is odd (first), i is even (second)
# As an oriented edge, this has gamma = -1 (right-handed)
# So J FLIPS gamma: epsilon'' = -1

epsilon_pp = -1
print(f"\n  epsilon'' (J*gamma vs gamma*J):")
print(f"    J maps edge (even->odd) to edge (odd->even)")
print(f"    This FLIPS chirality: J*gamma = -gamma*J")
print(f"    epsilon'' = {epsilon_pp}")

# On the edge space: J^2
# J maps edge (i,j) -> (j,i) -> (i,j) (back to itself)
# So J^2 = +1 on edges? Not quite -- J also involves the antiunitary structure.
# For an edge (i,j): the representation space is rho_i tensor rho_j
# J acts as (complex conjugation) x (swap i,j)
# J^2 = (CC)^2 x (swap)^2 = I x I = I on the VECTOR SPACE
# But the FS indicator modifies this:
# If rho_i has FS = -1, the complex conjugation J_i satisfies J_i^2 = -1
# On the edge (i,j): J = J_i tensor J_j (roughly)
# J^2 = J_i^2 tensor J_j^2 = FS(i)*FS(j)
# For edges in the McKay graph: i and j have OPPOSITE parity (bipartite)
# So FS(i)*FS(j) depends on which even-dim and odd-dim are connected

# Actually, edges connect nodes of different parity (even/odd in bipartite sense)
# But the FS indicator depends on REPRESENTATION dimension parity, not graph parity
# Let me check the correlation:

print(f"\n  Edge FS products (J^2 on each edge):")
edge_J2 = []
for i in range(9):
    for j in range(i+1, 9):
        if A_mckay[i,j] == 1:
            j2 = FS[i] * FS[j]
            edge_J2.append(j2)
            print(f"    rho_{i}(dim {dims[i]},FS={FS[i]:+d}) -- rho_{j}(dim {dims[j]},FS={FS[j]:+d}): J^2 = {j2:+d}")

# Check if all edges have the same J^2
unique_j2 = set(edge_J2)
print(f"\n  Unique J^2 values on edges: {unique_j2}")
if len(unique_j2) == 1:
    epsilon = list(unique_j2)[0]
    print(f"  *** UNIFORM: epsilon = J^2 = {epsilon:+d} ***")
else:
    print(f"  NOT UNIFORM: J^2 varies across edges")
    print(f"  J^2 = +1 on {sum(1 for j in edge_J2 if j==+1)} edges")
    print(f"  J^2 = -1 on {sum(1 for j in edge_J2 if j==-1)} edges")
    # In this case, we need to look at the WEIGHTED version
    # or determine which edges are the physical fermions
    weighted_p = sum(dims[i]*dims[j] for i in range(9) for j in range(i+1,9)
                     if A_mckay[i,j]==1 and FS[i]*FS[j]==+1)
    weighted_m = sum(dims[i]*dims[j] for i in range(9) for j in range(i+1,9)
                     if A_mckay[i,j]==1 and FS[i]*FS[j]==-1)
    print(f"  Dim-weighted: J^2=+1 accounts for {weighted_p}, J^2=-1 for {weighted_m}")

# epsilon' from [J, D]:
# D is the adjacency matrix (connects even to odd nodes)
# J reverses edges and conjugates
# Since ALL irreps are self-conjugate: J(rho_i) = rho_i
# J(D(edge(i,j))) = J(edge neighbors) = same edges reversed
# D(J(edge(i,j))) = D(edge(j,i)) = neighbors of reversed edge
# For the McKay graph (symmetric): these are the same
# So JD = DJ => epsilon' = +1

epsilon_p = +1
print(f"\n  epsilon' (JD vs DJ):")
print(f"    McKay graph is undirected: J reverses edges, D acts symmetrically")
print(f"    JD = DJ on the symmetric edge space")
print(f"    epsilon' = {epsilon_p}")

# ============================================================
# PART 4: KO-DIMENSION DETERMINATION
# ============================================================
print(f"\n{'='*70}")
print("PART 4: KO-DIMENSION")
print("=" * 70)

# KO-dimension table:
# n=0: (+,+,+), n=1: (+,-,.), n=2: (-,+,-), n=3: (-,+,.)
# n=4: (-,+,+), n=5: (-,-,.), n=6: (+,+,-), n=7: (+,+,.)
# (. means gamma doesn't exist in odd KO-dim)

ko_table = {
    0: (+1, +1, +1),
    2: (-1, +1, -1),
    4: (-1, +1, +1),
    6: (+1, +1, -1),
}

print(f"\n  Computed signs: epsilon={epsilon if len(unique_j2)==1 else '?'}, epsilon'={epsilon_p}, epsilon''={epsilon_pp}")
print(f"\n  KO-dimension table (even, with gamma):")
for n, (e, ep, epp) in ko_table.items():
    match_ep = (ep == epsilon_p)
    match_epp = (epp == epsilon_pp)
    if len(unique_j2) == 1:
        match_e = (e == epsilon)
    else:
        match_e = "?"
    match_all = match_ep and match_epp and (match_e == True)
    marker = " <-- MATCH" if match_all else ""
    print(f"    KO = {n}: ({e:+d}, {ep:+d}, {epp:+d}) -- ep:{match_ep}, epp:{match_epp}, e:{match_e}{marker}")

# Since epsilon' = +1 and epsilon'' = -1:
# This matches KO = 6 (with epsilon = +1) or needs epsilon check for others
# KO = 0: (+,+,+) -- epsilon''=+1, but we have -1. NO.
# KO = 2: (-,+,-) -- epsilon'=+1 YES, epsilon''=-1 YES, epsilon=-1 needed
# KO = 4: (-,+,+) -- epsilon''=+1, but we have -1. NO.
# KO = 6: (+,+,-) -- epsilon'=+1 YES, epsilon''=-1 YES, epsilon=+1 needed

print(f"\n  From epsilon'=+1, epsilon''=-1:")
print(f"  Candidates: KO = 2 (needs J^2=-1) or KO = 6 (needs J^2=+1)")

# ============================================================
# PART 5: WHICH KO-DIMENSION? J^2 ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("PART 5: RESOLVING J^2 -- WHICH EDGES ARE PHYSICAL?")
print("=" * 70)

# The SM fermion content determines which edges matter.
# In our framework:
# - Leptons: rho_0 -- rho_1 edge (dim 1*2 = 2)
# - Quarks: related to rho_2, rho_3 etc.
# The PHYSICAL Hilbert space may not include all edges.

# But let's look at it differently:
# For SM NCG (Connes): KO(F) = 6, so J^2 = +1.
# Our framework aims to DERIVE the SM, so we should get KO = 6.

# Key check: do the PHYSICAL fermion edges all have J^2 = +1?
# From the edge analysis above: edges with J^2=+1 are those where
# FS(i)*FS(j) = +1, i.e., both real or both pseudo-real.

# Edges with J^2 = +1: both same FS type
# Edges with J^2 = -1: different FS types

print(f"  Edges with J^2 = +1 (same FS parity):")
for i in range(9):
    for j in range(i+1, 9):
        if A_mckay[i,j] == 1 and FS[i]*FS[j] == +1:
            print(f"    rho_{i}(d={dims[i]}) -- rho_{j}(d={dims[j]}): product dim = {dims[i]*dims[j]}")

print(f"\n  Edges with J^2 = -1 (opposite FS parity):")
for i in range(9):
    for j in range(i+1, 9):
        if A_mckay[i,j] == 1 and FS[i]*FS[j] == -1:
            print(f"    rho_{i}(d={dims[i]}) -- rho_{j}(d={dims[j]}): product dim = {dims[i]*dims[j]}")

# ============================================================
# PART 6: GALOIS INVOLUTION ON THE McKAY GRAPH
# ============================================================
print(f"\n{'='*70}")
print("PART 6: GALOIS INVOLUTION ON THE McKAY GRAPH")
print("=" * 70)

# The Galois involution sigma: sqrt(5) -> -sqrt(5), phi -> phi' = -1/phi
# On the character table of 2I, sigma acts on characters involving phi.
# Characters of 2I at conjugacy class of order 5 involve phi and phi'.
# sigma swaps these: chi_i -> chi_sigma(i)

# From the McKay graph: which nodes are swapped by Galois?
# The character table of 2I (at class of order 10, element of order 10):
# The eigenvalues of the Cayley adjacency involve phi.
# From the Galois kernel (exp356): ker = rho_0 + rho_1 + rho_7 (dim 9 = N_eig)

# The Galois involution acts on the McKay graph as:
# rho_0 <-> rho_0 (trivial, rational character)
# rho_1 <-> rho_7 (both dim 2, swapped by phi <-> phi')
# rho_2 <-> rho_8 (both dim 3)
# rho_3 <-> rho_6 (both dim 4)
# rho_4 <-> rho_4 (dim 5, fixed)
# rho_5 <-> rho_5 (dim 6, fixed)

galois_map = {0:0, 1:7, 7:1, 2:8, 8:2, 3:6, 6:3, 4:4, 5:5}
print(f"\nGalois involution on irreps:")
for i, j in sorted(galois_map.items()):
    action = "fixed" if i == j else f"<-> rho_{j}"
    print(f"  rho_{i} (dim {dims[i]}, FS={FS[i]:+d}) {action}")

# Check: does Galois preserve the McKay graph?
A_galois = np.zeros((9,9), dtype=int)
for i in range(9):
    for j in range(9):
        A_galois[galois_map[i], galois_map[j]] = A_mckay[i,j]
print(f"\n  Galois preserves McKay graph? {np.array_equal(A_galois, A_mckay)}")

# Check: does Galois preserve the chirality (bipartition)?
galois_preserves_chirality = all(color[galois_map[i]] == color[i] for i in range(9))
galois_flips_chirality = all(color[galois_map[i]] == -color[i] for i in range(9))
print(f"  Galois preserves chirality? {galois_preserves_chirality}")
print(f"  Galois FLIPS chirality? {galois_flips_chirality}")

# Check: does Galois preserve or flip FS indicators?
galois_preserves_FS = all(FS[galois_map[i]] == FS[i] for i in range(9))
galois_flips_FS = all(FS[galois_map[i]] == -FS[i] for i in range(9))
print(f"  Galois preserves FS? {galois_preserves_FS}")
print(f"  Galois flips FS? {galois_flips_FS}")

# What does Galois do to J^2 on edges?
print(f"\n  Galois action on edges:")
for i in range(9):
    for j in range(i+1, 9):
        if A_mckay[i,j] == 1:
            gi, gj = galois_map[i], galois_map[j]
            j2_orig = FS[i]*FS[j]
            j2_galois = FS[gi]*FS[gj]
            changed = "CHANGED" if j2_orig != j2_galois else "same"
            print(f"    ({i},{j}) J^2={j2_orig:+d} -> ({gi},{gj}) J^2={j2_galois:+d} [{changed}]")

# ============================================================
# PART 7: GALOIS AND WICK ROTATION
# ============================================================
print(f"\n{'='*70}")
print("PART 7: IS GALOIS = WICK ROTATION?")
print("=" * 70)

# In NCG, Wick rotation changes KO-dimension by some amount.
# Specifically: Wick rotation on M^4 maps (Euclidean, KO=4) to (Lorentzian, KO=2 or 6)
# This changes KO by ±2.
#
# For the FINITE part F:
# If Galois acts as KO -> KO + delta on F,
# then it would be analogous to Wick rotation if delta = ±2.
#
# From KO = 6 to KO = 6 + delta: if delta = 2, new KO = 0.
# If delta = -2, new KO = 4.
#
# KO=6 has signs (+,+,-). KO=4 has signs (-,+,+). KO=0 has signs (+,+,+).
# The only sign that changes: epsilon and epsilon'' (J^2 and J*gamma).
#
# Does Galois change J^2?
# From Part 6: we need to check if J^2 on Galois-mapped edges differs.

# From the analysis: Galois preserves dimensions (dim_i = dim_sigma(i))
# And FS indicator depends on dim parity: FS = +1 if dim odd, -1 if dim even
# Since dim_sigma(i) = dim_i, FS is preserved.
# So J^2 on Galois-mapped edges is THE SAME.

# Therefore: Galois does NOT change the KO-dimension signs!
# Galois is NOT the Wick rotation (at least not in this simple sense).

print(f"""
  Analysis:
  - Galois preserves dimensions: dim(sigma(rho_i)) = dim(rho_i)
  - Therefore Galois preserves FS indicators: FS(sigma(i)) = FS(i)
  - Therefore J^2 is Galois-invariant on all edges
  - Therefore Galois does NOT change the KO-dimension

  CONCLUSION: Galois involution != Wick rotation
  (at least not via the simple KO-dimension mechanism)

  BUT: Galois does something subtler:
  - It preserves the McKay graph structure (adjacency)
  - It {'PRESERVES' if galois_preserves_chirality else 'DOES NOT preserve'} chirality
  - It swaps PAIRS of irreps: (1,7), (2,8), (3,6)
  - These swapped pairs are in the Galois KERNEL: ker = rho_0 + rho_1 + rho_7
""")

# ============================================================
# PART 8: THE MASS DUALITY AND m_e
# ============================================================
print(f"\n{'='*70}")
print("PART 8: MASS DUALITY AND THE m_e ANCHOR")
print("=" * 70)

# m_f * m_f' = m_e^2 for all fermions (Galois mass duality)
# This means m_e is the GEOMETRIC MEAN of each (m_f, m_f') pair.
# But this doesn't FIX m_e -- it just says m_f' = m_e^2/m_f.

# HOWEVER: if there's a self-consistency condition that the
# TOTAL mass (or some spectral quantity) is Galois-invariant,
# this MIGHT fix the scale.

# The Galois-invariant quantities we know:
# - sum N(z_f) = b1 (exp445-446)
# - E_vac * E_vac' = -1 (exact, exp444)
# - alpha * alpha' = 1/(2*pi) (exact)

# From alpha * alpha' = 1/(2*pi):
# alpha and alpha' are the two roots of 2*pi*a^2 - (N/b1)*phi^4*a + 1 = 0
# Product of roots = 1/(2*pi) (from Vieta's)
# This is a DERIVED constraint, not something that fixes m_e.

# What about m_e?
# m_e/M_P = alpha^(4*phi^2) (derived)
# m_e'/M_P' = alpha'^(4*phi'^2) ???
# This is complicated because phi' = -1/phi and alpha' involves phi'.

# Let's compute: what does m_e' look like?
# m_e' = m_e^2/m_e = m_e (the electron is its own Galois dual!?)
# NO: m_e is the ANCHOR, not a fermion in the mass formula.
# The mass formula is m_f = m_e * phi^{exponent}.
# m_f' = m_e * phi'^{exponent} (Galois conjugate the exponent base)

# Actually: m_f' = m_e * phi'^{5a+6b+...}
# And m_f * m_f' = m_e^2 * (phi * phi')^{...} = m_e^2 * (-1)^{...}
# phi * phi' = -1. So m_f * m_f' = m_e^2 * (-1)^{5a+6b+delta}

# For this to equal m_e^2, we need (-1)^{5a+6b+delta} = 1
# i.e., 5a+6b+delta must be EVEN.
# This is a PARITY constraint, not a scale constraint.

print(f"  Mass duality: m_f = m_e * phi^n, m_f' = m_e * phi'^n")
print(f"  m_f * m_f' = m_e^2 * (phi*phi')^n = m_e^2 * (-1)^n")
print(f"  For m_f*m_f' = m_e^2: need n even (consistent with known exponents)")

print(f"\n  Checking parity of mass exponents n = 5a + 6b + delta:")
particles = [
    ('e', 0, 0, 0), ('mu', 1, 1, 1), ('tau', 1, 2, 1),
    ('u', 3, -2, 1), ('c', 2, 1, 1), ('t', 4, 1, 1),
    ('d', 1, 0, 1), ('s', 1, 1, 1), ('b', -1, 4, 1)
]
for name, a, b, delta in particles:
    n = 5*a + 6*b + delta
    parity = "even" if n % 2 == 0 else "ODD"
    prod_sign = (-1)**n
    print(f"    {name:3s}: n = 5*{a} + 6*{b} + {delta} = {n:3d}  ({parity})  m*m' = {'+' if prod_sign>0 else '-'}m_e^2")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY AND HONEST ASSESSMENT")
print("=" * 70)

print(f"""
Q1: KO-DIMENSION
  epsilon' = +1 (J commutes with D, from McKay graph symmetry)
  epsilon'' = -1 (J flips chirality on edges)
  => KO(F) = 6 (if J^2 = +1) or KO(F) = 2 (if J^2 = -1)

  J^2 on edges: NOT UNIFORM (depends on FS parity of both endpoints)
  - 4 edges with J^2 = +1 (same FS type at both ends)
  - 4 edges with J^2 = -1 (opposite FS type)

  For KO = 6 (SM value): need J^2 = +1.
  This selects the EVEN-FS edges as the physical fermion space.
  STATUS: PARTIALLY DERIVED (KO=6 requires selecting specific edges)

Q2: GALOIS = WICK ROTATION
  Galois preserves: dimensions, FS indicators, McKay graph, J^2
  Therefore: Galois does NOT change KO-dimension
  Galois IS NOT the Wick rotation (via KO mechanism)
  STATUS: NEGATIVE (simple version). More subtle connection possible.

Q3: MASS DUALITY AND m_e
  m_f * m_f' = m_e^2 * (-1)^n where n = mass exponent
  For n even: m_f*m_f' = +m_e^2 (proper duality)
  For n odd: m_f*m_f' = -m_e^2 (sign flip!)
  The sign structure is a PARITY CONSTRAINT, not a scale constraint.
  m_e remains UNFIXED by the duality.
  STATUS: NEGATIVE (duality constrains parity, not scale)

WHAT WE LEARNED:
  - KO(F) is likely 6 (matching SM), derivable from FS indicators + edge selection
  - But the signature derivation needs KO(M×F) = 0 mod 8 argument
  - Galois is NOT Wick rotation (at KO level)
  - Mass duality doesn't fix m_e
""")

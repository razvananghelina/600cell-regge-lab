"""
exp351: Chirality from the McKay Graph Bipartite Structure
==========================================================
KEY INSIGHT: The McKay graph of 2I (= affine E8 Dynkin diagram) is a TREE.
Trees are BIPARTITE. This gives a natural Z/2 grading gamma_F on the internal
Hilbert space H_F with {gamma_F, D_F} = 0 EXACT.

APPROACH: Use ANALYTICAL character table (known SU(2) angles per class)
to avoid numerical issues with Galois conjugation.
"""

import numpy as np
from collections import deque

PHI = (1 + np.sqrt(5)) / 2
PHI_PRIME = (1 - np.sqrt(5)) / 2
a1 = 5
b1 = 6

# ============================================================
# PART 1: Character Table of 2I (Analytical)
# ============================================================
print("=" * 70)
print("PART 1: Analytical Character Table of 2I")
print("=" * 70)

# 9 conjugacy classes of 2I
# SU(2) half-angle alpha: eigenvalues e^{+-i*alpha}
# Class name, size, alpha
class_data = [
    ("1A",   1, 0),           # identity
    ("2A",   1, np.pi),       # central element -I
    ("4A",  30, np.pi/2),     # order 4 (pi-rotations)
    ("6A",  20, np.pi/3),     # order 6
    ("3A",  20, 2*np.pi/3),   # order 3
    ("10A", 12, np.pi/5),     # order 10 (2pi/5 rotations)
    ("5A",  12, 4*np.pi/5),   # order 5 (= -1 * 10A)
    ("5B",  12, 2*np.pi/5),   # order 5 (4pi/5 rotations)
    ("10B", 12, 3*np.pi/5),   # order 10 (= -1 * 5B)
]
class_names = [c[0] for c in class_data]
class_sizes = [c[1] for c in class_data]
alphas = [c[2] for c in class_data]
N_classes = 9
N_group = 120

print(f"\n9 conjugacy classes of 2I (|2I| = {N_group}):")
for name, size, alpha in class_data:
    chi2 = 2*np.cos(alpha)
    print(f"  {name:4s}: size={size:2d}, alpha={alpha/np.pi:.4f}*pi, "
          f"chi_std={chi2:+.4f}")

# 9 irreps of 2I
# Sym^k(std) for k=0,1,...,5 gives dims 1,2,3,4,5,6
# Sym^k(std') for k=1,2 gives another dim-2 and dim-3
# std tensor std' gives the second dim-4
# Character formulas:
#   Sym^k(std): chi = sin((k+1)*alpha) / sin(alpha)
#   Sym^k(std'): chi = sin((k+1)*alpha') / sin(alpha')
#     where alpha' = Galois conjugate angle
#   std tensor std': chi = chi_std * chi_std'

# Galois action on the angle: at class with chi_std = 2*cos(alpha),
# chi_std' = sigma(chi_std) where sigma: sqrt(5) -> -sqrt(5)

print("\nComputing characters analytically...")

def chi_spin(k, alpha):
    """Character of Sym^k(std) = spin-k/2 at SU(2) half-angle alpha."""
    j = k / 2.0  # spin
    dim = k + 1
    if abs(np.sin(alpha)) < 1e-10:
        # At identity (alpha=0) or central (alpha=pi)
        if abs(alpha) < 1e-10:
            return float(dim)
        else:  # alpha = pi
            return float(dim * (-1)**k)
    return np.sin((k+1)*alpha) / np.sin(alpha)

# For Galois conjugate: need sigma(2*cos(alpha))
# sigma maps phi -> phi', which maps 2cos(pi/5)=phi to -1/phi=2cos(3pi/5)
# So sigma swaps alpha and alpha' at 5-fold classes
# At rational classes (1A,2A,4A,6A,3A): sigma is identity

def galois_angle(alpha):
    """Galois conjugate of the SU(2) angle."""
    # At 5-fold and 10-fold classes, Galois swaps the two
    # pairs: (pi/5) <-> (3pi/5) and (2pi/5) <-> (4pi/5)
    # At other classes: alpha is unchanged
    eps = 1e-8
    if abs(alpha - np.pi/5) < eps:
        return 3*np.pi/5
    elif abs(alpha - 3*np.pi/5) < eps:
        return np.pi/5
    elif abs(alpha - 2*np.pi/5) < eps:
        return 4*np.pi/5
    elif abs(alpha - 4*np.pi/5) < eps:
        return 2*np.pi/5
    else:
        return alpha  # rational class, Galois-fixed

# Build the 9x9 character table
# Irreps in order:
# rho_1 = Sym^0(std) = trivial (dim 1)
# rho_2 = Sym^1(std) = std (dim 2)
# rho_3 = Sym^1(std') (dim 2)
# rho_4 = Sym^2(std) (dim 3)
# rho_5 = Sym^2(std') (dim 3)
# rho_6 = std tensor std' (dim 4)
# rho_7 = Sym^3(std) (dim 4)
# rho_8 = Sym^4(std) (dim 5)
# rho_9 = Sym^5(std) (dim 6)

irrep_names = [
    "rho_1: Sym^0(std), trivial",
    "rho_2: Sym^1(std) = std",
    "rho_3: Sym^1(std') = std'",
    "rho_4: Sym^2(std)",
    "rho_5: Sym^2(std')",
    "rho_6: std x std'",
    "rho_7: Sym^3(std)",
    "rho_8: Sym^4(std)",
    "rho_9: Sym^5(std)",
]
irrep_dims = [1, 2, 2, 3, 3, 4, 4, 5, 6]

chars = np.zeros((9, N_classes))
for c in range(N_classes):
    a = alphas[c]
    a_gal = galois_angle(a)

    chars[0, c] = chi_spin(0, a)         # rho_1 = trivial
    chars[1, c] = chi_spin(1, a)         # rho_2 = std
    chars[2, c] = chi_spin(1, a_gal)     # rho_3 = std'
    chars[3, c] = chi_spin(2, a)         # rho_4 = Sym^2(std)
    chars[4, c] = chi_spin(2, a_gal)     # rho_5 = Sym^2(std')
    chars[5, c] = chars[1,c] * chars[2,c]  # rho_6 = std x std'
    chars[6, c] = chi_spin(3, a)         # rho_7 = Sym^3(std)
    chars[7, c] = chi_spin(4, a)         # rho_8 = Sym^4(std)
    chars[8, c] = chi_spin(5, a)         # rho_9 = Sym^5(std)

print("\nCharacter table of 2I:")
header = "      dim  " + "  ".join(f"{n:>6s}" for n in class_names)
print(header)
sizes_row = "     sizes " + "  ".join(f"{s:>6d}" for s in class_sizes)
print(sizes_row)
for k in range(9):
    vals = "  ".join(f"{chars[k,c]:+6.3f}" for c in range(N_classes))
    print(f"  rho_{k+1}  {irrep_dims[k]}  {vals}")

# Verify orthogonality
print("\nOrthogonality verification:")
max_err = 0
for i in range(9):
    for j in range(9):
        ip = sum(class_sizes[c] * chars[i,c] * chars[j,c] for c in range(N_classes))
        expected = N_group if i == j else 0
        err = abs(ip - expected)
        max_err = max(max_err, err)
        if err > 0.1:
            print(f"  ERROR: <rho_{i+1}|rho_{j+1}> = {ip:.2f} (expected {expected})")
print(f"  Max orthogonality error: {max_err:.2e}")
if max_err < 0.1:
    print("  ALL orthogonality relations VERIFIED")

# Verify: sum d_i^2 = |G|
sum_d2 = sum(d**2 for d in irrep_dims)
print(f"\n  sum(d_i^2) = {sum_d2} (expected {N_group})")


# ============================================================
# PART 2: McKay Graph from Character Inner Products
# ============================================================
print("\n" + "=" * 70)
print("PART 2: McKay Graph = Affine E8")
print("=" * 70)

# McKay graph: edge from rho_i to rho_j iff rho_j appears in rho_2 x rho_i
# Multiplicity = (1/|G|) sum_g chi_j(g)* chi_2(g) chi_i(g)
fund_idx = 1  # rho_2 is the fundamental

mckay_adj = np.zeros((9, 9), dtype=int)
for i in range(9):
    product_chars = chars[fund_idx] * chars[i]
    for j in range(9):
        mult = sum(class_sizes[c] * chars[j,c] * product_chars[c]
                   for c in range(N_classes)) / N_group
        mult_int = int(round(mult))
        if mult_int > 0:
            mckay_adj[i, j] = mult_int

print("\nMcKay adjacency matrix:")
for i in range(9):
    nbrs = [f"rho_{j+1}(d={irrep_dims[j]})" for j in range(9) if mckay_adj[i,j] > 0]
    print(f"  rho_{i+1}(d={irrep_dims[i]}): {', '.join(nbrs)}")

# McKay property check: 2*d_i = sum of neighbor dims
print("\nMcKay property check (2*d_i = sum neighbor dims):")
all_ok = True
for i in range(9):
    nbr_sum = sum(mckay_adj[i,j]*irrep_dims[j] for j in range(9))
    ok = (nbr_sum == 2*irrep_dims[i])
    if not ok:
        all_ok = False
    print(f"  rho_{i+1}: 2*{irrep_dims[i]} = {2*irrep_dims[i]}, sum = {nbr_sum} {'OK' if ok else 'FAIL'}")

# Count edges and check tree
n_edges = sum(1 for i in range(9) for j in range(i+1,9) if mckay_adj[i,j] > 0)
print(f"\nTotal edges: {n_edges} (expected 8 for tree on 9 nodes)")
is_tree = (n_edges == 8)
print(f"Is a TREE: {is_tree}")

# Find branch node (degree 3)
branch = -1
for i in range(9):
    degree = sum(1 for j in range(9) if mckay_adj[i,j] > 0)
    if degree == 3:
        branch = i
        print(f"\nBranch node: rho_{i+1} (dim {irrep_dims[i]}, degree {degree})")

# Find legs from branch
def find_legs(adj, br):
    legs = []
    nbrs = [j for j in range(9) if adj[br,j] > 0]
    for start in nbrs:
        leg = [start]
        current, prev = start, br
        while True:
            nexts = [j for j in range(9) if adj[current,j] > 0 and j != prev]
            if not nexts:
                break
            prev, current = current, nexts[0]
            leg.append(current)
        legs.append(leg)
    return legs

if branch >= 0:
    legs = find_legs(mckay_adj, branch)
    legs.sort(key=len)
    print(f"\nLegs from branch (sorted by length):")
    for idx, leg in enumerate(legs):
        chain = " -- ".join(f"rho_{n+1}(d={irrep_dims[n]})" for n in leg)
        print(f"  Leg {idx} (len {len(leg)}): {chain}")
    print(f"  Leg lengths: {[len(l) for l in legs]} (expected [1, 2, 5] for E8)")


# ============================================================
# PART 3: Bipartite Structure
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Bipartite Partition and Chirality Grading")
print("=" * 70)

if not is_tree:
    print("WARNING: Graph is NOT a tree! Checking bipartiteness anyway...")

# BFS coloring from rho_1 (index 0)
color = [-1] * 9
color[0] = 0  # WHITE
queue = deque([0])
is_bipartite = True
while queue:
    node = queue.popleft()
    for j in range(9):
        if mckay_adj[node, j] > 0:
            if color[j] < 0:
                color[j] = 1 - color[node]
                queue.append(j)
            elif color[j] == color[node]:
                is_bipartite = False

white = [i for i in range(9) if color[i] == 0]
black = [i for i in range(9) if color[i] == 1]
dim_white = sum(irrep_dims[i] for i in white)
dim_black = sum(irrep_dims[i] for i in black)

print(f"\nBipartite: {is_bipartite}")
print(f"\nWHITE (gamma_F = +1):")
for i in white:
    print(f"  rho_{i+1} (dim {irrep_dims[i]}): {irrep_names[i]}")
print(f"  Total dim = {dim_white}")

print(f"\nBLACK (gamma_F = -1):")
for i in black:
    print(f"  rho_{i+1} (dim {irrep_dims[i]}): {irrep_names[i]}")
print(f"  Total dim = {dim_black}")

print(f"\n  dim(WHITE) = {dim_white}")
print(f"  dim(BLACK) = {dim_black}")
print(f"  Sum = {dim_white + dim_black}")

# Verify all edges cross partition
print("\nEdge crossing check:")
for i in range(9):
    for j in range(i+1, 9):
        if mckay_adj[i,j] > 0:
            crosses = (color[i] != color[j])
            ci = 'W' if color[i]==0 else 'B'
            cj = 'W' if color[j]==0 else 'B'
            print(f"  rho_{i+1}--rho_{j+1}: {ci}--{cj} {'CROSSES' if crosses else 'SAME (NOT BIPARTITE!)'}")


# ============================================================
# PART 4: THEOREM - {gamma_F, D_F} = 0
# ============================================================
print("\n" + "=" * 70)
print("PART 4: THEOREM - {gamma_F, D_F} = 0")
print("=" * 70)

if is_bipartite:
    print("""
THEOREM: The McKay graph of 2I is bipartite (it is a tree = affine E8).
The bipartite partition defines gamma_F on H_F (dim """ + str(dim_white+dim_black) + """):
  gamma_F|_{V_rho} = +1 if rho in WHITE, -1 if rho in BLACK

For ANY Hermitian D_F with blocks only between adjacent McKay nodes:
  {gamma_F, D_F} = 0   (EXACT)

Proof: For adjacent nodes i,j: color[i] != color[j] (bipartite).
  So gamma_F[i] = -gamma_F[j].
  Block (i,j) of anticommutator: gamma_F[i]*D[i,j] + D[i,j]*gamma_F[j]
    = chi_i*D[i,j] + D[i,j]*(-chi_i) = 0.
  Non-adjacent blocks: D[i,j] = 0. Diagonal: D[i,i] = 0.  QED.
""")

    # Numerical verification
    print("Numerical verification with random D_F...")
    total_dim = sum(irrep_dims)
    gamma_F = np.zeros((total_dim, total_dim))
    block_starts = {}
    offset = 0
    for k in range(9):
        block_starts[k] = offset
        sign = 1.0 if color[k] == 0 else -1.0
        for i in range(irrep_dims[k]):
            gamma_F[offset+i, offset+i] = sign
        offset += irrep_dims[k]

    np.random.seed(42)
    D_F = np.zeros((total_dim, total_dim), dtype=complex)
    for i in range(9):
        for j in range(i+1, 9):
            if mckay_adj[i,j] > 0:
                di, dj = irrep_dims[i], irrep_dims[j]
                block = np.random.randn(di, dj) + 1j*np.random.randn(di, dj)
                si, sj = block_starts[i], block_starts[j]
                D_F[si:si+di, sj:sj+dj] = block
                D_F[sj:sj+dj, si:si+di] = block.conj().T

    anticomm = gamma_F @ D_F + D_F @ gamma_F
    print(f"  ||{{gamma_F, D_F}}|| = {np.max(np.abs(anticomm)):.2e} (should be ~0)")
    print(f"  gamma_F^2 = I: {np.allclose(gamma_F@gamma_F, np.eye(total_dim))}")
    comm_df2 = gamma_F @ D_F @ D_F - D_F @ D_F @ gamma_F
    print(f"  ||[gamma_F, D_F^2]|| = {np.max(np.abs(comm_df2)):.2e} (should be ~0)")

    # Eigenvalue pairing
    evals = np.sort(np.linalg.eigvalsh(D_F))
    print(f"\n  D_F eigenvalue pairing (lambda_k + lambda_{{n+1-k}}):")
    for k in range(total_dim//2):
        s = evals[k] + evals[total_dim-1-k]
        print(f"    pair {k+1}: {evals[k]:+8.4f} + {evals[total_dim-1-k]:+8.4f} = {s:+.6f}")
else:
    print("  Graph is NOT bipartite - theorem does not apply in simple form.")
    print("  Need to investigate modified grading.")


# ============================================================
# PART 5: Fermion Content (Solution 3)
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Fermion Content of Each Partition")
print("=" * 70)

if branch >= 0 and len(legs) == 3:
    # Solution 3 from exp323:
    # Main chain from endpoint of long leg to branch:
    # e -> u -> d -> s -> mu -> c(branch)
    # Mid leg: tau -> t
    # Short leg: b
    long_leg = [l for l in legs if len(l) == 5][0]
    mid_leg = [l for l in legs if len(l) == 2][0]
    short_leg = [l for l in legs if len(l) == 1][0]

    # Long leg goes from branch outward. Reverse to get endpoint first.
    chain_from_e = list(reversed(long_leg))

    sol3 = {}
    fermion_list = ['e', 'u', 'd', 's', 'mu']
    for i, f in enumerate(fermion_list):
        sol3[chain_from_e[i]] = f
    sol3[branch] = 'c'
    sol3[mid_leg[0]] = 'tau'
    sol3[mid_leg[1]] = 't'
    sol3[short_leg[0]] = 'b'

    print("\nSolution 3 fermion-node assignment:")
    for node in range(9):
        f = sol3.get(node, '?')
        col = 'WHITE' if color[node]==0 else 'BLACK'
        print(f"  rho_{node+1} (dim {irrep_dims[node]}) = {f:5s}  [{col}]  {irrep_names[node]}")

    print("\n  CHIRAL PARTITION:")
    wf = [sol3[i] for i in white]
    bf = [sol3[i] for i in black]
    print(f"    WHITE (gamma_F=+1, dim {dim_white}): {', '.join(wf)}")
    print(f"    BLACK (gamma_F=-1, dim {dim_black}): {', '.join(bf)}")

    # Classify
    print("\n  By fermion type:")
    for label, group in [("charged leptons", ['e','mu','tau']),
                         ("up-type quarks", ['u','c','t']),
                         ("down-type quarks", ['d','s','b'])]:
        in_w = [f for f in wf if f in group]
        in_b = [f for f in bf if f in group]
        print(f"    {label}: WHITE={in_w}, BLACK={in_b}")


# ============================================================
# PART 6: Galois Action on McKay Nodes
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Galois Action on Irreps")
print("=" * 70)

# Apply sigma to character table: sigma(chi(g)) at each class
# sigma swaps 5-fold classes: 10A<->10B, 5A<->5B
# (because sigma maps phi to phi', swapping cos(pi/5)<->cos(3pi/5) etc.)

galois_partner = list(range(9))
for k in range(9):
    sigma_chars = np.zeros(N_classes)
    for c in range(N_classes):
        val = chars[k, c]
        # Decompose: val = a + b*phi, sigma(val) = (a+b) - b*phi
        best_err = 1e10
        best_sigma = val
        for b_try in np.arange(-6, 7, 1):
            b = b_try / 2.0
            a = val - b * PHI
            if abs(a - round(2*a)/2) < 1e-6:
                sigma_val = round(2*a)/2 + b - b * PHI
                if abs(sigma_val - round(sigma_val)) < 1e-6:
                    sigma_val = round(sigma_val)
                best_sigma = sigma_val
                best_err = 0
                break
        sigma_chars[c] = best_sigma

    # Find matching irrep
    for l in range(9):
        if all(abs(sigma_chars[c] - chars[l,c]) < 0.01 for c in range(N_classes)):
            galois_partner[k] = l
            break

print("\nGalois action sigma on 2I irreps:")
for k in range(9):
    p = galois_partner[k]
    status = "FIXED" if p == k else f"-> rho_{p+1}"
    f = sol3.get(k, '?') if branch >= 0 else '?'
    print(f"  rho_{k+1}(d={irrep_dims[k]}, {f}): sigma {status}")

# Summary
fixed = [k for k in range(9) if galois_partner[k] == k]
pairs = []
seen = set()
for k in range(9):
    if galois_partner[k] != k and k not in seen:
        pairs.append((k, galois_partner[k]))
        seen.update([k, galois_partner[k]])

print(f"\nGalois-FIXED ({len(fixed)} irreps, dim {sum(irrep_dims[k] for k in fixed)}):")
for k in fixed:
    print(f"  rho_{k+1}(d={irrep_dims[k]})")
print(f"Galois-SWAPPED ({len(pairs)} pairs):")
for k, l in pairs:
    ck = 'W' if color[k]==0 else 'B'
    cl = 'W' if color[l]==0 else 'B'
    print(f"  rho_{k+1}(d={irrep_dims[k]},{ck}) <-> rho_{l+1}(d={irrep_dims[l]},{cl})")

# Check if swapped pairs have same or different bipartite color
print("\nBipartite color of Galois-swapped pairs:")
for k, l in pairs:
    ck = 'W' if color[k]==0 else 'B'
    cl = 'W' if color[l]==0 else 'B'
    same = "SAME" if color[k]==color[l] else "DIFFERENT"
    print(f"  rho_{k+1}({ck}) <-> rho_{l+1}({cl}): {same} color")


# ============================================================
# PART 7: Summary and Framework Connections
# ============================================================
print("\n" + "=" * 70)
print("PART 7: SUMMARY")
print("=" * 70)

print(f"""
MAIN RESULT: FINITE CHIRALITY FROM MCKAY GRAPH
================================================

The McKay graph of 2I = affine E8 Dynkin diagram.
It is a TREE on 9 nodes with 8 edges.
Trees are bipartite => natural Z/2 grading gamma_F.

  WHITE (gamma_F = +1): dims {[irrep_dims[i] for i in white]}, total = {dim_white}
  BLACK (gamma_F = -1): dims {[irrep_dims[i] for i in black]}, total = {dim_black}

THEOREM: {{gamma_F, D_F}} = 0 for ANY D_F respecting McKay topology.
  This is the chirality condition for NCG: D = D_M x 1 + gamma_5 x D_F

THE {dim_white}/{dim_black} SPLIT:
  dim(WHITE) = {dim_white}
  dim(BLACK) = {dim_black}
  Sum = {dim_white+dim_black} = h(E8) = a1 * b1 = Coxeter number
  {dim_white} = (a1-1)^2 = dim(ST)^2 = N(z_Planck) = fermions/SM generation
  {dim_white}/{dim_black} = {dim_white/dim_black if dim_black > 0 else 'inf'}

GALOIS vs BIPARTITE:
  Galois swaps {len(pairs)} pair(s), fixes {len(fixed)} irrep(s).
  Bipartite and Galois gradings: INDEPENDENT Z/2 structures.

ASSESSMENT:
  - gamma_F existence: DERIVED (tree structure of affine E8)
  - {{gamma_F, D_F}} = 0: PROVED (theorem, not numerical)
  - gamma_F^2 = I: trivial (diagonal +/-1)
  - [gamma_F, D_F^2] = 0: corollary
  - dim(WHITE) = {dim_white} = fermions/generation: {'MATCHES' if dim_white == 16 else 'DOES NOT MATCH'}
  - Physical L/R identification: OPEN (which partition is L?)

CATEGORY: DERIVED (chirality operator exists with correct properties)
          OPEN (physical identification of L vs R sectors)
""")

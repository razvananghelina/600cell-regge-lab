"""
exp337 Step 3: Deep spectral analysis of D_F
=============================================
Now that D_F is built (from v2), analyze its structure in depth:
1. Eigenvector projections onto each node
2. Spectral functions: D^2, |D|, sign(D)*|D|^p
3. Laplacian-weighted D_F (Cayley eigenvalues as diagonal mass)
4. T3-graded version (chirality grading from Hopf)
5. Green's function / resolvent analysis
6. Connection between 13 positive eigenvalues and 9 fermion masses
"""
import numpy as np
from itertools import product, permutations
from math import comb

PHI = (1 + np.sqrt(5)) / 2

# ================================================================
# Reproduce the group and D_F quickly (from v2 verified code)
# ================================================================
def quat_to_su2(q):
    a, b, c, d = q
    return np.array([[a + 1j*b, -c + 1j*d],
                      [c + 1j*d,  a - 1j*b]])

def generate_2I():
    quats = []
    for i in range(4):
        for s in [1, -1]:
            q = [0, 0, 0, 0]
            q[i] = s
            quats.append(tuple(q))
    for signs in product([0.5, -0.5], repeat=4):
        quats.append(signs)
    vals = [0, 1/(2*PHI), 0.5, PHI/2]
    even_perms = [
        (0,1,2,3), (0,2,3,1), (0,3,1,2),
        (1,0,3,2), (1,2,0,3), (1,3,2,0),
        (2,0,1,3), (2,1,3,0), (2,3,0,1),
        (3,0,2,1), (3,1,0,2), (3,2,1,0),
    ]
    for perm in even_perms:
        base = [vals[perm[0]], vals[perm[1]], vals[perm[2]], vals[perm[3]]]
        nonzero_pos = [i for i in range(4) if abs(base[i]) > 1e-10]
        for signs in product([1, -1], repeat=len(nonzero_pos)):
            q = list(base)
            for idx, s in zip(nonzero_pos, signs):
                q[idx] = abs(q[idx]) * s
            quats.append(tuple(q))
    unique = []
    for q in quats:
        if not any(np.allclose(q, u, atol=1e-10) for u in unique):
            unique.append(q)
    return [quat_to_su2(q) for q in unique]

def extract_irrep_from_tensor(group, reps, i, j, chi_j_vals, d_j, N):
    d_i = reps[i][0].shape[0]
    d_prod = 2 * d_i
    P = np.zeros((d_prod, d_prod), dtype=complex)
    for g in range(N):
        tensor = np.kron(reps[1][g], reps[i][g])
        P += np.conj(chi_j_vals[g]) * tensor
    P *= d_j / N
    eigvals, eigvecs = np.linalg.eigh(P)
    basis_indices = np.where(eigvals > 0.5)[0]
    V = eigvecs[:, basis_indices]
    if V.shape[1] != d_j:
        print(f"  WARNING: Expected {d_j} basis vectors, got {V.shape[1]}")
        U, S, _ = np.linalg.svd(P)
        V = U[:, :d_j]
    rep_j = []
    for g in range(N):
        tensor = np.kron(reps[1][g], reps[i][g])
        R = V.conj().T @ tensor @ V
        rep_j.append(R)
    return rep_j, V

print("="*72)
print("exp337 Step 3: Deep spectral analysis")
print("="*72)

group = generate_2I()
N = len(group)
traces = np.array([np.trace(g).real for g in group])
trace_vals = [2.0, PHI, 1.0, 1/PHI, 0.0, -1/PHI, -1.0, -PHI, -2.0]
char_table = np.array([
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [2, PHI, 1, 1/PHI, 0, -1/PHI, -1, -PHI, -2],
    [3, PHI, 0, -1/PHI, -1, -1/PHI, 0, PHI, 3],
    [4, 1, -1, -1, 0, 1, 1, -1, -4],
    [5, 0, -1, 0, 1, 0, -1, 0, 5],
    [6, -1, 0, 1, 0, -1, 0, 1, -6],
    [4, -1, 1, -1, 0, -1, 1, -1, 4],
    [2, -1/PHI, 1, -PHI, 0, PHI, -1, 1/PHI, -2],
    [3, -1/PHI, 0, PHI, -1, PHI, 0, -1/PHI, 3],
])
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

def element_to_class(tr):
    for ci, tv in enumerate(trace_vals):
        if abs(tr - tv) < 1e-6: return ci
    return -1

element_class = np.array([element_to_class(t) for t in traces])
chi = np.zeros((9, N))
for j in range(9):
    for g in range(N):
        chi[j, g] = char_table[j, element_class[g]]

# Build all reps
print("Building representations...")
reps = {}
cg = {}

reps[0] = [np.array([[1.0+0j]]) for _ in range(N)]
reps[1] = [g.copy() for g in group]
reps[2], cg[(1,2)] = extract_irrep_from_tensor(group, reps, 1, 2, chi[2], 3, N)
_, cg[(1,0)] = extract_irrep_from_tensor(group, reps, 1, 0, chi[0], 1, N)
reps[3], cg[(2,3)] = extract_irrep_from_tensor(group, reps, 2, 3, chi[3], 4, N)
_, cg[(2,1)] = extract_irrep_from_tensor(group, reps, 2, 1, chi[1], 2, N)
reps[4], cg[(3,4)] = extract_irrep_from_tensor(group, reps, 3, 4, chi[4], 5, N)
_, cg[(3,2)] = extract_irrep_from_tensor(group, reps, 3, 2, chi[2], 3, N)
reps[5], cg[(4,5)] = extract_irrep_from_tensor(group, reps, 4, 5, chi[5], 6, N)
_, cg[(4,3)] = extract_irrep_from_tensor(group, reps, 4, 3, chi[3], 4, N)
reps[6], cg[(5,6)] = extract_irrep_from_tensor(group, reps, 5, 6, chi[6], 4, N)
reps[8], cg[(5,8)] = extract_irrep_from_tensor(group, reps, 5, 8, chi[8], 3, N)
_, cg[(5,4)] = extract_irrep_from_tensor(group, reps, 5, 4, chi[4], 5, N)
reps[7], cg[(6,7)] = extract_irrep_from_tensor(group, reps, 6, 7, chi[7], 2, N)
_, cg[(6,5)] = extract_irrep_from_tensor(group, reps, 6, 5, chi[5], 6, N)
_, cg[(0,1)] = extract_irrep_from_tensor(group, reps, 0, 1, chi[1], 2, N)

offsets = [0]
for d in dims:
    offsets.append(offsets[-1] + d)

edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (5,8), (6,7)]

def build_dirac(v_higgs, normalize=None):
    D = np.zeros((30, 30), dtype=complex)
    for i_src, j_tgt in edges:
        V = cg[(i_src, j_tgt)]
        d_i = dims[i_src]
        d_j = dims[j_tgt]
        V3 = V.reshape(2, d_i, d_j)
        block = np.conj(v_higgs[0]) * V3[0] + np.conj(v_higgs[1]) * V3[1]
        if normalize == 'unit_frob':
            nrm = np.linalg.norm(block, 'fro')
            if nrm > 1e-10: block /= nrm
        elif normalize == 'inv_sqrt2':
            block /= np.sqrt(2)
        oi, oj = offsets[i_src], offsets[j_tgt]
        D[oi:oi+d_i, oj:oj+d_j] = block
        D[oj:oj+d_j, oi:oi+d_i] = block.conj().T
    return D

D = build_dirac(np.array([1.0, 0.0]))
print("D_F built (30x30)")

# ================================================================
# 1. EIGENVECTOR ANALYSIS
# ================================================================
print("\n" + "="*72)
print("1. EIGENVECTOR ANALYSIS")
print("="*72)

evals, evecs = np.linalg.eigh(D)
print(f"\nEigenvalues: {np.round(evals, 5)}")
print(f"\nPositive eigenvalues: {np.round(evals[evals > 1e-8], 5)}")

# For each eigenvector, compute its projection weight on each node
# Weight_k(lambda) = ||psi_lambda restricted to rho_k||^2 / d_k
print(f"\n{'lambda':>10} |", end="")
for k in range(9):
    print(f" rho_{k}({dims[k]})", end="")
print()
print("-"*110)

proj_matrix = np.zeros((30, 9))  # (eigenvalue index, node index) = weight
for ev_idx in range(30):
    psi = evecs[:, ev_idx]
    weights = []
    for k in range(9):
        ok = offsets[k]
        w = np.sum(np.abs(psi[ok:ok+dims[k]])**2)
        weights.append(w)
        proj_matrix[ev_idx, k] = w

    # Only print positive eigenvalues and zeros
    if evals[ev_idx] >= -1e-8:
        print(f"{evals[ev_idx]:>10.5f} |", end="")
        for w in weights:
            if w > 0.01:
                print(f"  {w:>6.3f} ", end="")
            else:
                print(f"    .    ", end="")
        print()

# ================================================================
# 2. WHICH EIGENVALUES "BELONG" TO WHICH NODES?
# ================================================================
print("\n" + "="*72)
print("2. NODE ASSIGNMENT OF EIGENVALUES")
print("="*72)

# For each node k, the "spectral weight" is the set of eigenvalues
# weighted by projection
print("\nPrimary eigenvalue(s) per node (largest projection):")
for k in range(9):
    # Which positive eigenvalue has the largest projection on node k?
    pos_mask = evals > 1e-8
    pos_idx = np.where(pos_mask)[0]
    weights_k = proj_matrix[pos_idx, k]
    best_idx = pos_idx[np.argmax(weights_k)]
    print(f"  rho_{k} (dim {dims[k]}): lambda={evals[best_idx]:>8.5f} (weight={proj_matrix[best_idx, k]:.3f})")

# Effective eigenvalue per node (weighted mean)
print("\nEffective eigenvalue per node (weighted by projection):")
pos_idx = np.where(evals > 1e-8)[0]
eff_evals = []
for k in range(9):
    weights_k = proj_matrix[pos_idx, k]
    if np.sum(weights_k) > 1e-10:
        eff = np.sum(evals[pos_idx] * weights_k) / np.sum(weights_k)
    else:
        eff = 0.0
    eff_evals.append(eff)
    print(f"  rho_{k}: <lambda>_k = {eff:.5f}")


# ================================================================
# 3. SPECTRAL FUNCTIONS
# ================================================================
print("\n" + "="*72)
print("3. SPECTRAL FUNCTIONS f(D_F)")
print("="*72)

# Compute various functions of D_F and extract per-node values
D2 = D @ D
D_abs = evecs @ np.diag(np.abs(evals)) @ evecs.conj().T

# Per-node mean of |lambda|
abs_per_node = []
for k in range(9):
    ok = offsets[k]
    tr = np.trace(D_abs[ok:ok+dims[k], ok:ok+dims[k]]).real
    abs_per_node.append(tr / dims[k])

# Per-node mean of lambda^2
sq_per_node = []
for k in range(9):
    ok = offsets[k]
    tr = np.trace(D2[ok:ok+dims[k], ok:ok+dims[k]]).real
    sq_per_node.append(tr / dims[k])

# Per-node mean of lambda^4
D4 = D2 @ D2
q4_per_node = []
for k in range(9):
    ok = offsets[k]
    tr = np.trace(D4[ok:ok+dims[k], ok:ok+dims[k]]).real
    q4_per_node.append(tr / dims[k])

print(f"\n{'Node':>6} {'dim':>4} {'Tr|D|/d':>10} {'Tr(D^2)/d':>10} {'Tr(D^4)/d':>10} {'<lambda>':>10}")
print("-"*55)
for k in range(9):
    print(f"rho_{k:>2} {dims[k]:>4} {abs_per_node[k]:>10.5f} {sq_per_node[k]:>10.5f} {q4_per_node[k]:>10.5f} {eff_evals[k]:>10.5f}")


# ================================================================
# 4. LAPLACIAN-WEIGHTED D_F
# ================================================================
print("\n" + "="*72)
print("4. LAPLACIAN-WEIGHTED D_F")
print("="*72)

# Cayley eigenvalues
chi_C1 = char_table[:, 1]  # class with trace=phi
cayley_eigs = 12 * (1 - chi_C1 / np.array(dims))

# Build D_F with diagonal Cayley mass terms
# D_tot = D_F + diag(cayley_i * I_{d_i})
for scale in [0.01, 0.05, 0.1, 0.5, 1.0]:
    D_weighted = D.copy()
    for k in range(9):
        ok = offsets[k]
        for i in range(dims[k]):
            D_weighted[ok+i, ok+i] = scale * cayley_eigs[k]

    evals_w = np.sort(np.linalg.eigvalsh(D_weighted))
    pos_w = evals_w[evals_w > 1e-8]
    print(f"\n  scale={scale}: {len(pos_w)} positive eigs")
    if len(pos_w) <= 15:
        print(f"    {np.round(pos_w, 4)}")

# ================================================================
# 5. D_F WITH DIMENSION GRADING
# ================================================================
print("\n" + "="*72)
print("5. DIMENSION-GRADED D_F")
print("="*72)

# Weight each edge by some function of the node dimensions
# Attempt: the "Yukawa" as Y_{ij} = V[0]/sqrt(d_i * d_j)
for weight_fn in ['none', 'inv_geom_mean', 'inv_d_src', 'inv_d_tgt', 'sqrt_d_ratio']:
    D_w = np.zeros((30, 30), dtype=complex)
    v = np.array([1.0, 0.0])
    for i_src, j_tgt in edges:
        V = cg[(i_src, j_tgt)]
        d_i = dims[i_src]
        d_j = dims[j_tgt]
        V3 = V.reshape(2, d_i, d_j)
        block = np.conj(v[0]) * V3[0] + np.conj(v[1]) * V3[1]

        if weight_fn == 'inv_geom_mean':
            block /= np.sqrt(d_i * d_j)
        elif weight_fn == 'inv_d_src':
            block /= d_i
        elif weight_fn == 'inv_d_tgt':
            block /= d_j
        elif weight_fn == 'sqrt_d_ratio':
            block *= np.sqrt(d_j / d_i) if d_i > 0 else 1.0

        oi, oj = offsets[i_src], offsets[j_tgt]
        D_w[oi:oi+d_i, oj:oj+d_j] = block
        D_w[oj:oj+d_j, oi:oi+d_i] = block.conj().T

    evals_w = np.sort(np.linalg.eigvalsh(D_w))

    # Per-node Tr(D^2)/d
    D2_w = D_w @ D_w
    node_vals = []
    for k in range(9):
        ok = offsets[k]
        tr = np.trace(D2_w[ok:ok+dims[k], ok:ok+dims[k]]).real
        node_vals.append(tr / dims[k])

    print(f"\n  weight={weight_fn}:")
    print(f"    Tr(D^2)/d: {[f'{v:.4f}' for v in node_vals]}")


# ================================================================
# 6. TRY TO MATCH 13 -> 9 EIGENVALUES
# ================================================================
print("\n" + "="*72)
print("6. MATCHING 13 POSITIVE EIGENVALUES TO 9 FERMIONS")
print("="*72)

pos_evals = np.sort(evals[evals > 1e-8])
print(f"\n13 positive eigenvalues:")
for i, ev in enumerate(pos_evals):
    print(f"  lambda_{i:>2} = {ev:.6f}")

# The 4 zero modes suggest dim(ker) = 4
# What node(s) contribute to the kernel?
zero_mask = np.abs(evals) < 1e-8
zero_evecs = evecs[:, zero_mask]
print(f"\n4 zero modes - node projections:")
for z in range(zero_evecs.shape[1]):
    psi = zero_evecs[:, z]
    print(f"  zero mode {z}: ", end="")
    for k in range(9):
        ok = offsets[k]
        w = np.sum(np.abs(psi[ok:ok+dims[k]])**2)
        if w > 0.01:
            print(f"rho_{k}({w:.3f}) ", end="")
    print()

# Mass data
delta_vals = np.array([0.000, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278])
fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
mass_exponents = np.array([0, 11, 17, 3, 16, 26, 5, 11, 19])

# Strategy: Use eigenvector projections to assign each node its
# "primary" positive eigenvalue. With 13 evals and 9 nodes,
# some nodes claim >1 eigenvalue. The node's effective mass contribution
# is some function of its eigenvalues.

# Method A: Each node takes the eigenvalue with highest projection
print("\nMethod A: Max-projection assignment")
assigned = {}
for k in range(9):
    best_eval = -1
    best_weight = 0
    for i, ev_idx in enumerate(np.where(evals > 1e-8)[0]):
        w = proj_matrix[ev_idx, k]
        if w > best_weight:
            best_weight = w
            best_eval = evals[ev_idx]
    assigned[k] = best_eval
    print(f"  rho_{k} -> lambda={best_eval:.5f} (weight={best_weight:.3f})")

# Method B: Weighted mean of ALL positive eigenvalues per node
print("\nMethod B: Weighted mean eigenvalue per node")
for k in range(9):
    pos_idx = np.where(evals > 1e-8)[0]
    ws = proj_matrix[pos_idx, k]
    mean_ev = np.sum(evals[pos_idx] * ws) / np.sum(ws) if np.sum(ws) > 0 else 0
    print(f"  rho_{k} -> <lambda>={mean_ev:.5f}")


# ================================================================
# 7. THE KEY: Express per-node values in terms of phi
# ================================================================
print("\n" + "="*72)
print("7. EXACT VALUES IN TERMS OF PHI")
print("="*72)

# Per-node Tr(D^2)/d is topological. Let's compute EXACT rational form.
# Formula: Tr(D^2 at k)/d_k = (1/(2*d_k)) * [sum over edges touching k of d_target]
# where d_target is the dimension of the EDGE TARGET in the CG.

print("\nPer-node Tr(D^2)/d (exact):")
for k in range(9):
    contrib = 0
    for i_src, j_tgt in edges:
        if k == i_src:
            # k is source, target dim = d_j_tgt
            contrib += dims[j_tgt]
        elif k == j_tgt:
            # k is target, contributes d_k
            contrib += dims[k]
    exact = contrib / (2 * dims[k])
    print(f"  rho_{k}: sum_neighbor_d = {contrib}, Tr(D^2)/d = {contrib}/(2*{dims[k]}) = {exact:.6f}")
    # Express as fraction
    from fractions import Fraction
    frac = Fraction(contrib, 2*dims[k])
    print(f"           = {frac}")


# ================================================================
# 8. COMBINED OPERATOR: D_F + alpha * L_Cayley
# ================================================================
print("\n" + "="*72)
print("8. COMBINED OPERATORS")
print("="*72)

# Build diagonal Cayley Laplacian in the 30-dim space
L_cayley_30 = np.zeros((30, 30))
for k in range(9):
    ok = offsets[k]
    for i in range(dims[k]):
        L_cayley_30[ok+i, ok+i] = cayley_eigs[k]

# D_tot = D_F + alpha * L_Cayley
# This breaks chiral symmetry (adds mass-like diagonal terms)
for alpha_val in [0.01, 0.02, 0.05, 0.1]:
    D_tot = D + alpha_val * L_cayley_30
    evals_tot = np.sort(np.linalg.eigvalsh(D_tot))

    # Now eigenvalues are NOT in +/- pairs
    # Take the 9 largest positive as candidates
    pos_tot = evals_tot[evals_tot > 0]

    # Per-node "effective mass": diagonal of D_tot
    diag_vals = []
    for k in range(9):
        ok = offsets[k]
        diag_vals.append(D_tot[ok, ok].real)

    print(f"\nalpha={alpha_val}:")
    print(f"  9 largest pos eigs: {np.round(pos_tot[-9:], 4)}")


# ================================================================
# 9. THE DEFINITIVE TEST: eigenvalues of (D_F^2)_restricted
# ================================================================
print("\n" + "="*72)
print("9. D_F^2 RESTRICTED TO EACH NODE")
print("="*72)

# D_F^2 restricted to node k: the d_k x d_k block of D^2
# Its eigenvalues tell us the mass-squared spectrum at node k
D2 = D @ D
print("\nEigenvalues of D_F^2 restricted to each node:")
for k in range(9):
    ok = offsets[k]
    block = D2[ok:ok+dims[k], ok:ok+dims[k]]
    block_eigs = np.sort(np.linalg.eigvalsh(block))
    print(f"  rho_{k} (dim {dims[k]}): {np.round(block_eigs, 5)}")
    # Are they all equal? If so, Schur's lemma applies
    if dims[k] > 1:
        spread = np.max(block_eigs) - np.min(block_eigs)
        mean = np.mean(block_eigs)
        print(f"    spread = {spread:.6f}, mean = {mean:.6f}, {'DEGENERATE' if spread < 1e-4 else 'NON-DEGENERATE'}")


# ================================================================
# 10. COMMUTATOR STRUCTURE
# ================================================================
print("\n" + "="*72)
print("10. COMMUTATOR AND ANTICOMMUTATOR")
print("="*72)

# Gamma_5 grading: assign +1 to even-dim nodes, -1 to odd-dim nodes?
# Or: chirality from the graph bipartition
# McKay graph (E8 affine) IS bipartite!
# Bipartition: {rho_0, rho_2, rho_4, rho_6} vs {rho_1, rho_3, rho_5, rho_7, rho_8}
# Check: dims 1+3+5+4=13 vs 2+4+6+2+3=17. Total 30.
# This is NOT equal, so gamma_5 is not a standard Z/2 grading.

# Actually for affine E8: the graph IS bipartite (it's a tree).
# Partition A = {0, 2, 4, 6, 8} (sum dims = 1+3+5+4+3=16)
# Partition B = {1, 3, 5, 7} (sum dims = 2+4+6+2=14)
# Hmm 16+14=30. NOT equal.

# Build chirality operator
gamma = np.zeros(30)
partition_A = [0, 2, 4, 6, 8]  # Even distance from rho_0
partition_B = [1, 3, 5, 7]     # Odd distance from rho_0

# Wait - but rho_8 connects to rho_5 (B), so rho_8 should be in A.
# And rho_6 connects to rho_5 (B), so rho_6 should be in A.
# rho_7 connects to rho_6 (A), so rho_7 is in B. Correct.
# So A = {0, 2, 4, 6, 8}, B = {1, 3, 5, 7}

for k in partition_A:
    ok = offsets[k]
    gamma[ok:ok+dims[k]] = 1.0
for k in partition_B:
    ok = offsets[k]
    gamma[ok:ok+dims[k]] = -1.0

Gamma = np.diag(gamma)

# Check: D anticommutes with Gamma if the graph is bipartite
anticomm = D @ Gamma + Gamma @ D
print(f"\n||{{D, Gamma}}|| = {np.linalg.norm(anticomm):.2e} {'(=0: anticommutes!)' if np.linalg.norm(anticomm) < 1e-10 else '(NOT zero)'}")

# This means D has a chiral structure!
# The chirality splits H = H+ (dim 16) + H- (dim 14)
print(f"H+ (A-nodes, dim {sum(dims[k] for k in partition_A)}): rho_0, rho_2, rho_4, rho_6, rho_8")
print(f"H- (B-nodes, dim {sum(dims[k] for k in partition_B)}): rho_1, rho_3, rho_5, rho_7")

# D maps H+ -> H- and H- -> H+
# The "Dirac mass" is the off-diagonal part
# D = [[0, M], [M^dag, 0]] where M is 16x14 (or 14x16)

# Extract M
idx_A = []
for k in partition_A:
    idx_A.extend(range(offsets[k], offsets[k]+dims[k]))
idx_B = []
for k in partition_B:
    idx_B.extend(range(offsets[k], offsets[k]+dims[k]))

M = D[np.ix_(idx_A, idx_B)]
print(f"\nM matrix (H+ -> H-): shape {M.shape}")
print(f"Singular values of M: {np.round(np.linalg.svd(M, compute_uv=False), 5)}")

# M M^dag eigenvalues = positive eigenvalues of D^2 restricted to H+
MM = M @ M.conj().T
MM_eigs = np.sort(np.linalg.eigvalsh(MM))
print(f"\nEigenvalues of M*M^dag (dim {MM.shape[0]}): {np.round(MM_eigs, 5)}")

# M^dag M eigenvalues = positive eigenvalues of D^2 restricted to H-
MdM = M.conj().T @ M
MdM_eigs = np.sort(np.linalg.eigvalsh(MdM))
print(f"Eigenvalues of M^dag*M (dim {MdM.shape[0]}): {np.round(MdM_eigs, 5)}")


# ================================================================
# 11. FINAL COMPARISON TABLE
# ================================================================
print("\n" + "="*72)
print("11. FINAL COMPARISON")
print("="*72)

# Compute various per-node quantities
print(f"\n{'Node':>6} {'dim':>4} {'part':>5} {'D2/d':>8} {'|D|/d':>8} {'<lam>':>8} {'delta':>8} {'n':>4}")
print("-"*65)
for k in range(9):
    ok = offsets[k]
    d2_d = np.trace(D2[ok:ok+dims[k], ok:ok+dims[k]]).real / dims[k]
    abs_d = np.trace(D_abs[ok:ok+dims[k], ok:ok+dims[k]]).real / dims[k]
    part = 'A' if k in partition_A else 'B'
    print(f"rho_{k:>2} {dims[k]:>4} {part:>5} {d2_d:>8.4f} {abs_d:>8.4f} {eff_evals[k]:>8.4f} {delta_vals[k]:>+8.4f} {mass_exponents[k]:>4}")

# Key insight: the branch nodes (A-partition rho_6, rho_8)
# vs main chain nodes - does the partition match T3 sectors?
print("\nPartition A (rho_0,2,4,6,8): includes BOTH main chain AND branch")
print("Partition B (rho_1,3,5,7): includes main chain AND branch endpoint")
print("=> Bipartition does NOT match lepton/quark or T3 splitting")

# Check if any linear combination of spectral quantities matches deltas
print("\n--- Linear fit: delta ~ a * D2/d + b * Cayley + c ---")
# Assemble feature matrix
X = np.column_stack([
    sq_per_node,
    cayley_eigs,
    np.ones(9)
])
# Least squares
coeffs, residuals, rank, sv = np.linalg.lstsq(X, delta_vals, rcond=None)
fitted = X @ coeffs
residual_rms = np.sqrt(np.mean((delta_vals - fitted)**2))
print(f"Coefficients: D2/d: {coeffs[0]:.4f}, Cayley: {coeffs[1]:.4f}, const: {coeffs[2]:.4f}")
print(f"Residual RMS: {residual_rms:.4f}")
print(f"Fitted vs actual:")
for k in range(9):
    print(f"  rho_{k}: fitted={fitted[k]:+.4f}, actual={delta_vals[k]:+.4f}, err={fitted[k]-delta_vals[k]:+.4f}")

print("\n--- Linear fit: delta ~ a * <lambda> + b * Cayley + c ---")
X2 = np.column_stack([
    eff_evals,
    cayley_eigs,
    np.ones(9)
])
coeffs2, _, _, _ = np.linalg.lstsq(X2, delta_vals, rcond=None)
fitted2 = X2 @ coeffs2
rms2 = np.sqrt(np.mean((delta_vals - fitted2)**2))
print(f"Coefficients: <lam>: {coeffs2[0]:.4f}, Cayley: {coeffs2[1]:.4f}, const: {coeffs2[2]:.4f}")
print(f"Residual RMS: {rms2:.4f}")

# Also try dim alone
print("\n--- Linear fit: delta ~ a * dim + b ---")
X3 = np.column_stack([np.array(dims, dtype=float), np.ones(9)])
coeffs3, _, _, _ = np.linalg.lstsq(X3, delta_vals, rcond=None)
fitted3 = X3 @ coeffs3
rms3 = np.sqrt(np.mean((delta_vals - fitted3)**2))
print(f"Coefficients: dim: {coeffs3[0]:.4f}, const: {coeffs3[1]:.4f}")
print(f"Residual RMS: {rms3:.4f}")

# Try all pairs of spectral quantities
print("\n--- Best 2-parameter fits ---")
features = {
    'D2/d': np.array(sq_per_node),
    '|D|/d': np.array(abs_per_node),
    '<lam>': np.array(eff_evals),
    'Cayley': cayley_eigs,
    'dim': np.array(dims, dtype=float),
    'dist': np.array([5,4,3,2,1,0,1,2,1], dtype=float),
    'D4/d': np.array(q4_per_node),
}

best_rms = 999
best_combo = None
for n1, f1 in features.items():
    for n2, f2 in features.items():
        if n1 >= n2:
            continue
        X = np.column_stack([f1, f2, np.ones(9)])
        c, _, _, _ = np.linalg.lstsq(X, delta_vals, rcond=None)
        fit = X @ c
        rms = np.sqrt(np.mean((delta_vals - fit)**2))
        if rms < best_rms:
            best_rms = rms
            best_combo = (n1, n2, c, fit)

n1, n2, c, fit = best_combo
print(f"Best: {n1} + {n2}")
print(f"  Coeffs: {c[0]:.4f}*{n1} + {c[1]:.4f}*{n2} + {c[2]:.4f}")
print(f"  RMS = {best_rms:.4f}")
for k in range(9):
    print(f"  rho_{k}: fitted={fit[k]:+.4f}, actual={delta_vals[k]:+.4f}")

print("\n\nDone.")

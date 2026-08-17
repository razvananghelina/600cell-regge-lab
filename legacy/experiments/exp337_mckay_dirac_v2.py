"""
exp337 v2: D_F on McKay graph - CORRECT recursive construction
================================================================
Build all 9 irreps RECURSIVELY from rho_1 using CG projection.
This guarantees unitary representations at every step.

McKay graph (affine E8):
rho_0(1) - rho_1(2) - rho_2(3) - rho_3(4) - rho_4(5) - rho_5(6) - rho_8(3)
                                                          |
                                                         rho_6(4) - rho_7(2)
"""
import numpy as np
from itertools import product
from math import comb

PHI = (1 + np.sqrt(5)) / 2

# ================================================================
# Generate 2I group (same as before, verified)
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

# ================================================================
# Character computation (for projectors)
# ================================================================
def character_at_element(reps, k, g_idx):
    """Trace of reps[k][g_idx]."""
    return np.trace(reps[k][g_idx]).real

# ================================================================
# CG projection: extract irrep j from rho_1 x rho_i
# ================================================================
def extract_irrep_from_tensor(group, reps, i, j, chi_j_vals, d_j, N):
    """Extract irrep rho_j from rho_1 x rho_i using character projection.

    Returns:
      - rep_j: list of d_j x d_j unitary matrices (the irrep)
      - V: (2*d_i) x d_j isometry (the CG intertwiner)
    """
    d_i = reps[i][0].shape[0]
    d_prod = 2 * d_i

    # Build projector P = (d_j/|G|) sum_g conj(chi_j(g)) * (rho_1(g) kron rho_i(g))
    P = np.zeros((d_prod, d_prod), dtype=complex)
    for g in range(N):
        tensor = np.kron(reps[1][g], reps[i][g])
        P += np.conj(chi_j_vals[g]) * tensor
    P *= d_j / N

    # P is Hermitian idempotent, rank = d_j
    # Use eigendecomposition to find the range
    eigvals, eigvecs = np.linalg.eigh(P)

    # Select eigenvectors with eigenvalue close to 1
    basis_indices = np.where(eigvals > 0.5)[0]
    V = eigvecs[:, basis_indices]  # (2*d_i) x d_j isometry

    if V.shape[1] != d_j:
        print(f"  WARNING: Expected {d_j} basis vectors, got {V.shape[1]}")
        # Use SVD as fallback
        U, S, _ = np.linalg.svd(P)
        V = U[:, :d_j]

    # Compute the representation in this basis
    # rho_j(g) = V^dag * (rho_1(g) kron rho_i(g)) * V
    rep_j = []
    for g in range(N):
        tensor = np.kron(reps[1][g], reps[i][g])
        R = V.conj().T @ tensor @ V
        rep_j.append(R)

    return rep_j, V

# ================================================================
# MAIN COMPUTATION
# ================================================================
print("="*72)
print("exp337 v2: McKay D_F with recursive CG construction")
print("="*72)

group = generate_2I()
N = len(group)
print(f"Group order: {N}")
assert N == 120

# Known character table values (verified in Step 1)
# Classes ordered by trace: 2, 1.618, 1, 0.618, 0, -0.618, -1, -1.618, -2
# We need chi_j(g) for each GROUP ELEMENT g (not just class representatives)

# First, compute traces of rho_1 for all elements (= class identifier)
traces = np.array([np.trace(g).real for g in group])

# Character table of 2I (9 irreps x 9 classes)
# Ordered by trace: 2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2
trace_vals = [2.0, PHI, 1.0, 1/PHI, 0.0, -1/PHI, -1.0, -PHI, -2.0]

char_table = np.array([
    [1, 1, 1, 1, 1, 1, 1, 1, 1],              # rho_0
    [2, PHI, 1, 1/PHI, 0, -1/PHI, -1, -PHI, -2],  # rho_1
    [3, PHI, 0, -1/PHI, -1, -1/PHI, 0, PHI, 3],   # rho_2
    [4, 1, -1, -1, 0, 1, 1, -1, -4],            # rho_3
    [5, 0, -1, 0, 1, 0, -1, 0, 5],              # rho_4
    [6, -1, 0, 1, 0, -1, 0, 1, -6],             # rho_5
    [4, -1, 1, -1, 0, -1, 1, -1, 4],            # rho_6
    [2, -1/PHI, 1, -PHI, 0, PHI, -1, 1/PHI, -2], # rho_7
    [3, -1/PHI, 0, PHI, -1, PHI, 0, -1/PHI, 3],  # rho_8
])

dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

# Map each group element to its class index (by trace)
def element_to_class(tr):
    for ci, tv in enumerate(trace_vals):
        if abs(tr - tv) < 1e-6:
            return ci
    return -1

element_class = np.array([element_to_class(t) for t in traces])
assert np.all(element_class >= 0), "Some elements not classified!"

# chi_j(g) for each irrep j and element g
chi = np.zeros((9, N))
for j in range(9):
    for g in range(N):
        chi[j, g] = char_table[j, element_class[g]]

# ================================================================
# Build representations recursively
# ================================================================
print("\nBuilding representations recursively...")

reps = {}
cg_intertwiners = {}  # (i, j) -> V, the CG isometry

# rho_0: trivial
reps[0] = [np.array([[1.0+0j]]) for _ in range(N)]
print("  rho_0 (dim 1): trivial")

# rho_1: fundamental SU(2)
reps[1] = [g.copy() for g in group]
print("  rho_1 (dim 2): fundamental SU(2)")

# rho_2 from rho_1 x rho_1 (dim 4 = 3+1, extract 3-dim piece)
reps[2], cg_intertwiners[(1, 2)] = extract_irrep_from_tensor(
    group, reps, 1, 2, chi[2], 3, N)
print("  rho_2 (dim 3): extracted from rho_1 x rho_1")

# rho_0 is also in rho_1 x rho_1, get its CG intertwiner
_, cg_intertwiners[(1, 0)] = extract_irrep_from_tensor(
    group, reps, 1, 0, chi[0], 1, N)

# rho_3 from rho_1 x rho_2 (dim 6 = 4+2, extract 4-dim piece)
reps[3], cg_intertwiners[(2, 3)] = extract_irrep_from_tensor(
    group, reps, 2, 3, chi[3], 4, N)
print("  rho_3 (dim 4): extracted from rho_1 x rho_2")

# Also get CG for (2, 1)
_, cg_intertwiners[(2, 1)] = extract_irrep_from_tensor(
    group, reps, 2, 1, chi[1], 2, N)

# rho_4 from rho_1 x rho_3 (dim 8 = 5+3, extract 5-dim piece)
reps[4], cg_intertwiners[(3, 4)] = extract_irrep_from_tensor(
    group, reps, 3, 4, chi[4], 5, N)
print("  rho_4 (dim 5): extracted from rho_1 x rho_3")

_, cg_intertwiners[(3, 2)] = extract_irrep_from_tensor(
    group, reps, 3, 2, chi[2], 3, N)

# rho_5 from rho_1 x rho_4 (dim 10 = 6+4, extract 6-dim piece)
reps[5], cg_intertwiners[(4, 5)] = extract_irrep_from_tensor(
    group, reps, 4, 5, chi[5], 6, N)
print("  rho_5 (dim 6): extracted from rho_1 x rho_4")

_, cg_intertwiners[(4, 3)] = extract_irrep_from_tensor(
    group, reps, 4, 3, chi[3], 4, N)

# rho_6 and rho_8 from rho_1 x rho_5 (dim 12 = 5+4+3)
reps[6], cg_intertwiners[(5, 6)] = extract_irrep_from_tensor(
    group, reps, 5, 6, chi[6], 4, N)
print("  rho_6 (dim 4): extracted from rho_1 x rho_5")

reps[8], cg_intertwiners[(5, 8)] = extract_irrep_from_tensor(
    group, reps, 5, 8, chi[8], 3, N)
print("  rho_8 (dim 3): extracted from rho_1 x rho_5")

_, cg_intertwiners[(5, 4)] = extract_irrep_from_tensor(
    group, reps, 5, 4, chi[4], 5, N)

# rho_7 from rho_1 x rho_6 (dim 8 = 6+2, extract 2-dim piece)
reps[7], cg_intertwiners[(6, 7)] = extract_irrep_from_tensor(
    group, reps, 6, 7, chi[7], 2, N)
print("  rho_7 (dim 2): extracted from rho_1 x rho_6")

_, cg_intertwiners[(6, 5)] = extract_irrep_from_tensor(
    group, reps, 6, 5, chi[5], 6, N)

# Also need: CG for (7, 6) and (8, 5) - the "reverse" edges
_, cg_intertwiners[(7, 6)] = extract_irrep_from_tensor(
    group, reps, 7, 6, chi[6], 4, N)

_, cg_intertwiners[(8, 5)] = extract_irrep_from_tensor(
    group, reps, 8, 5, chi[5], 6, N)

# Also CG for edge (0, 1)
_, cg_intertwiners[(0, 1)] = extract_irrep_from_tensor(
    group, reps, 0, 1, chi[1], 2, N)


# ================================================================
# Verify all representations
# ================================================================
print("\nVerification:")
for k in range(9):
    # Unitarity
    unit_ok = all(np.allclose(reps[k][g] @ reps[k][g].conj().T,
                               np.eye(dims[k]), atol=1e-9) for g in range(N))
    # Character
    char_ok = all(abs(np.trace(reps[k][g]).real - chi[k, g]) < 1e-6
                  for g in range(N))
    # Rep property (spot check)
    rep_ok = True
    np.random.seed(42)
    for _ in range(30):
        i, j = np.random.randint(N, size=2)
        prod_mat = reps[k][i] @ reps[k][j]
        g_prod = group[i] @ group[j]
        found = False
        for m in range(N):
            if np.allclose(group[m], g_prod, atol=1e-8):
                if not np.allclose(prod_mat, reps[k][m], atol=1e-6):
                    rep_ok = False
                found = True
                break
        if not found:
            rep_ok = False

    status = "OK" if (unit_ok and char_ok and rep_ok) else "FAIL"
    details = f"unit={'OK' if unit_ok else 'FAIL'} char={'OK' if char_ok else 'FAIL'} rep={'OK' if rep_ok else 'FAIL'}"
    print(f"  rho_{k} (dim {dims[k]}): {status} [{details}]")


# ================================================================
# Verify CG intertwiners
# ================================================================
print("\n" + "="*72)
print("CG INTERTWINER VERIFICATION")
print("="*72)

edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (5,8), (6,7)]

for i_src, j_tgt in edges:
    # The CG intertwiner is stored as (i_src, j_tgt)
    V = cg_intertwiners[(i_src, j_tgt)]  # (2*d_i) x d_j
    d_i = dims[i_src]
    d_j = dims[j_tgt]

    # Verify: V^dag (rho_1 kron rho_i)(g) V = rho_j(g) for all g
    ok = True
    max_err = 0
    for g in range(N):
        tensor = np.kron(reps[1][g], reps[i_src][g])
        lhs = V.conj().T @ tensor @ V
        rhs = reps[j_tgt][g]
        err = np.max(np.abs(lhs - rhs))
        max_err = max(max_err, err)
        if err > 1e-6:
            ok = False

    frob = np.linalg.norm(V, 'fro')
    print(f"  Edge ({i_src},{j_tgt}): V is ({2*d_i}x{d_j}), "
          f"||V||_F={frob:.4f} (sqrt(d_j)={np.sqrt(d_j):.4f}), "
          f"intertwining: {'OK' if ok else 'FAIL'} (max_err={max_err:.2e})")


# ================================================================
# BUILD D_F: The 30x30 finite Dirac operator
# ================================================================
print("\n" + "="*72)
print("BUILDING D_F")
print("="*72)

offsets = [0]
for d in dims:
    offsets.append(offsets[-1] + d)
print(f"Block offsets: {offsets}")

# D_F construction: for each edge (i,j), the Yukawa block is obtained
# by contracting the CG intertwiner V: C^{d_j} -> C^2 x C^{d_i}
# with the Higgs vev v in C^2.
#
# D_F[i-block, j-block] = v^dag_tensor V (reshaped)
# where v^dag_tensor means: contract the C^2 index of V with v.
#
# V has shape (2*d_i, d_j). Reshape as (2, d_i, d_j).
# Then D_block = sum_a v_a^* * V[a, :, :] = d_i x d_j matrix.

def build_dirac(cg_intertwiners, edges, dims, offsets, v_higgs, normalize=None):
    """Build 30x30 D_F by contracting CG maps with Higgs vev."""
    D = np.zeros((30, 30), dtype=complex)

    for i_src, j_tgt in edges:
        V = cg_intertwiners[(i_src, j_tgt)]
        d_i = dims[i_src]
        d_j = dims[j_tgt]
        V3 = V.reshape(2, d_i, d_j)

        # Contract with Higgs vev
        block = np.conj(v_higgs[0]) * V3[0] + np.conj(v_higgs[1]) * V3[1]

        if normalize == 'unit_frob':
            nrm = np.linalg.norm(block, 'fro')
            if nrm > 1e-10:
                block /= nrm
        elif normalize == 'inv_sqrt2':
            block /= np.sqrt(2)

        # D_F[rho_i block, rho_j block] = block
        oi, oj = offsets[i_src], offsets[j_tgt]
        D[oi:oi+d_i, oj:oj+d_j] = block
        D[oj:oj+d_j, oi:oi+d_i] = block.conj().T

    return D


# Test multiple Higgs vev directions and normalizations
vevs = {
    '(1,0)': np.array([1.0, 0.0]),
    '(0,1)': np.array([0.0, 1.0]),
    '(1,1)/sqrt2': np.array([1.0, 1.0]) / np.sqrt(2),
    '(1,i)/sqrt2': np.array([1.0, 1j]) / np.sqrt(2),
}

norms = ['raw', 'unit_frob', 'inv_sqrt2']

print("\n--- Eigenvalue scan ---")
for v_label, v in vevs.items():
    for norm in norms:
        D = build_dirac(cg_intertwiners, edges, dims, offsets, v, norm)

        # Check Hermiticity
        herm = np.allclose(D, D.conj().T, atol=1e-10)

        eigs = np.sort(np.linalg.eigvalsh(D))
        n_pos = np.sum(eigs > 1e-8)
        n_neg = np.sum(eigs < -1e-8)
        n_zero = 30 - n_pos - n_neg

        trD2 = np.trace(D @ D).real

        pos_eigs = eigs[eigs > 1e-8]

        print(f"\n  vev={v_label}, norm={norm}, herm={herm}")
        print(f"    Eigs: {n_neg}(-) + {n_zero}(0) + {n_pos}(+), Tr(D^2)={trD2:.4f}")
        if len(pos_eigs) > 0:
            print(f"    Pos eigs: {np.round(pos_eigs, 5)}")


# ================================================================
# DETAILED ANALYSIS with best D_F
# ================================================================
print("\n" + "="*72)
print("DETAILED ANALYSIS")
print("="*72)

# Use vev=(1,0), raw normalization as the natural choice
D_main = build_dirac(cg_intertwiners, edges, dims, offsets,
                      np.array([1.0, 0.0]), 'raw')
eigs_main = np.sort(np.linalg.eigvalsh(D_main))
print(f"\nD_F eigenvalues (vev=(1,0), raw):")
print(f"  {np.round(eigs_main, 6)}")

# Per-node analysis: Tr(D_F^2) contribution from each node
D2 = D_main @ D_main
print(f"\nPer-node Tr(D^2)/d_i:")
node_mass2 = []
for k in range(9):
    ok = offsets[k]
    block = D2[ok:ok+dims[k], ok:ok+dims[k]]
    tr = np.trace(block).real
    node_mass2.append(tr / dims[k])
    print(f"  rho_{k} (dim {dims[k]}): Tr(D^2)/d = {tr/dims[k]:.6f}")

# Mass deltas to compare
delta_vals = [0.000, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278]
fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
mass_exponents = [0, 11, 17, 3, 16, 26, 5, 11, 19]

# The nodes of the McKay graph:
# rho_0(dim 1), rho_1(dim 2), ..., rho_8(dim 3)
# Fermion assignment is not fixed. Try all permutations? Too many (9! = 362880).
# Instead, try the natural ordering and a few motivated assignments.

print("\n--- Correlation analysis ---")
print(f"Node values (Tr(D^2)/d): {[f'{v:.4f}' for v in node_mass2]}")
print(f"Mass deltas:              {[f'{v:+.4f}' for v in delta_vals]}")
print(f"Mass exponents:           {mass_exponents}")

# Check correlation with natural ordering
arr_n = np.array(node_mass2)
arr_d = np.array(delta_vals)
arr_m = np.array(mass_exponents)

if np.std(arr_n) > 1e-10:
    corr_d = np.corrcoef(arr_n, arr_d)[0, 1]
    corr_m = np.corrcoef(arr_n, arr_m)[0, 1]
    print(f"\nCorrelation (Tr(D^2)/d vs deltas): {corr_d:.4f}")
    print(f"Correlation (Tr(D^2)/d vs mass exp): {corr_m:.4f}")

# Try: use eigenvalues of the 9x9 "contracted" Dirac operator
# where the (i,j) entry is ||D_F block||_F between nodes i and j
print("\n--- 9x9 contracted D_F ---")
D9 = np.zeros((9, 9))
for i_src, j_tgt in edges:
    oi, oj = offsets[i_src], offsets[j_tgt]
    block = D_main[oi:oi+dims[i_src], oj:oj+dims[j_tgt]]
    frob = np.linalg.norm(block, 'fro')
    D9[i_src, j_tgt] = frob
    D9[j_tgt, i_src] = frob

eigs9 = np.sort(np.linalg.eigvalsh(D9))
print(f"9x9 contracted D_F eigenvalues: {np.round(eigs9, 6)}")

# Also try: weighted adjacency with dim factors
D9w = np.zeros((9, 9))
for i_src, j_tgt in edges:
    D9w[i_src, j_tgt] = 1.0 / np.sqrt(dims[i_src] * dims[j_tgt])
    D9w[j_tgt, i_src] = D9w[i_src, j_tgt]

eigs9w = np.sort(np.linalg.eigvalsh(D9w))
print(f"9x9 dim-weighted adj eigs: {np.round(eigs9w, 6)}")

# Cayley Laplacian eigenvalues (from character theory)
# For the Cayley graph of 2I with generating set = 12 order-5 elements:
# lambda_i = 12 * (1 - chi_i(5-cycle) / d_i)
# Character of rho_i on 5-cycle class (class with trace = PHI, which is C1 in our ordering)
# Wait, class C1 has trace PHI = 1.618, size 12. This is the order-10 elements.
# The 5-cycle class: in the icosahedral group, 5-cycles have order 5.
# In 2I, the lift has order 10. The class with trace PHI has size 12, order 10.
# There are also the elements with trace 1/PHI (size 12, order 5).

# Actually, the Cayley graph uses the generating set. For the McKay graph,
# the "generator" is the 2-dim fundamental rho_1. The adjacency is:
# A_ij = multiplicity of rho_j in rho_1 x rho_i.

# The CAYLEY graph Laplacian eigenvalues of 2I on the Cayley graph
# with generating set S (sum over 12 nearest neighbors) are:
# lambda_i = |S| - chi_i(sum_S) / d_i = 12 - (sum_{s in S} chi_i(s)) / d_i

print("\n--- Cayley graph Laplacian (12-neighbor) ---")
# The 12 nearest neighbors in the 600-cell are the order-10 elements
# (the ones with trace = phi, class C1 of size 12)
# chi_i(C1) from our character table:
chi_C1 = char_table[:, 1]  # column 1 is the class with trace=phi
print(f"chi_i at C1 (trace=phi): {np.round(chi_C1, 4)}")

cayley_eigs = 12 * (1 - chi_C1 / np.array(dims))
print(f"Cayley Laplacian eigs (lambda_i = 12*(1-chi_i(5c)/d_i)):")
for k in range(9):
    print(f"  rho_{k}: lambda = {cayley_eigs[k]:.6f}")

print(f"\nCayley eigs sorted: {np.round(np.sort(cayley_eigs), 4)}")
print(f"Mass deltas sorted: {np.round(np.sort(delta_vals), 4)}")
print(f"Mass exponents:     {sorted(mass_exponents)}")

# Check: do Cayley eigenvalues correlate with mass exponents?
corr_cm = np.corrcoef(cayley_eigs, arr_m)[0, 1]
corr_cd = np.corrcoef(cayley_eigs, arr_d)[0, 1]
print(f"\nCorrelation (Cayley eig vs mass exp): {corr_cm:.4f}")
print(f"Correlation (Cayley eig vs delta):    {corr_cd:.4f}")

# Express in terms of phi
print("\n--- Cayley eigenvalues in terms of phi ---")
for k in range(9):
    lam = cayley_eigs[k]
    # Try to express as a + b*phi
    # lam = 12 - 12*chi/d
    # For rho_0: 12*(1-1/1) = 0
    # For rho_1: 12*(1-phi/2) = 12 - 6*phi = 12 - 6*1.618 = 2.292
    # Let's compute exactly
    chi_val = chi_C1[k]
    d = dims[k]
    # lambda = 12 - 12*chi/d
    # chi values are: 1, phi, phi, 1, 0, -1, -1, -1/phi, -1/phi
    # So lambda = 12 - 12*chi/d
    print(f"  rho_{k}: chi={chi_val:>7.4f}, d={d}, lambda = 12-12*{chi_val:.4f}/{d} = {lam:.6f}")


# ================================================================
# FINAL: Exhaustive comparison of D_F spectrum with mass data
# ================================================================
print("\n" + "="*72)
print("EXHAUSTIVE COMPARISON")
print("="*72)

# Approach 1: D_F eigenvalues (30 values) -> look for 9 that match deltas
# Approach 2: Per-node D_F^2 values (9 values) -> compare with deltas
# Approach 3: Cayley eigenvalues (9 values) -> compare with deltas
# Approach 4: Graph distance from branch point -> compare with mass exponents

# Graph distances from rho_5 (the branch point, dim 6 = b1):
distances_from_5 = {0: 5, 1: 4, 2: 3, 3: 2, 4: 1, 5: 0, 6: 1, 7: 2, 8: 1}
dist_arr = np.array([distances_from_5[k] for k in range(9)])

print(f"\nGraph distance from rho_5: {[distances_from_5[k] for k in range(9)]}")
print(f"Mass exponents:            {mass_exponents}")
print(f"Mass deltas:               {[f'{d:+.3f}' for d in delta_vals]}")
print(f"Cayley eigenvalues:        {[f'{v:.3f}' for v in cayley_eigs]}")
print(f"Per-node Tr(D^2)/d:        {[f'{v:.4f}' for v in node_mass2]}")

# Summary table
print(f"\n{'Node':>6} {'dim':>4} {'dist':>5} {'Cayley':>8} {'D^2/d':>8} {'delta':>8} {'n_mass':>7}")
print("-"*55)
for k in range(9):
    print(f"rho_{k:>2} {dims[k]:>4} {dist_arr[k]:>5} {cayley_eigs[k]:>8.4f} "
          f"{node_mass2[k]:>8.4f} {delta_vals[k]:>+8.4f} {mass_exponents[k]:>7}")

print("\n\nDone. See analysis above.")

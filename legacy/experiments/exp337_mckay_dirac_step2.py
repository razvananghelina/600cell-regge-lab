"""
exp337 Step 2: Build all 9 irreps, CG intertwiners, and D_F
=============================================================
Uses the 2I group from Step 1. Constructs:
1. Explicit representation matrices for all 9 irreps
2. CG intertwiners for all 8 edges of the McKay graph
3. The 30x30 finite Dirac operator D_F
4. Eigenvalues and comparison with mass deltas
"""
import numpy as np
from math import comb
from itertools import product
import sys

PHI = (1 + np.sqrt(5)) / 2

# ================================================================
# Reuse group generation from Step 1
# ================================================================
def quat_to_su2(q):
    a, b, c, d = q
    return np.array([[a + 1j*b, -c + 1j*d],
                      [c + 1j*d,  a - 1j*b]])

def generate_2I():
    elements = []
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

    unique_quats = []
    for q in quats:
        is_dup = False
        for uq in unique_quats:
            if np.allclose(q, uq, atol=1e-10):
                is_dup = True
                break
        if not is_dup:
            unique_quats.append(q)

    matrices = [quat_to_su2(q) for q in unique_quats]
    return matrices

def sym_power_matrix(M, k):
    if k == 0:
        return np.array([[1.0 + 0j]])
    if k == 1:
        return M.copy()

    a, b = M[0, 0], M[0, 1]
    c, d = M[1, 0], M[1, 1]

    dim = k + 1
    R = np.zeros((dim, dim), dtype=complex)

    for m_col in range(dim):
        for p_row in range(dim):
            val = 0
            for i in range(max(0, p_row - m_col), min(k - m_col, p_row) + 1):
                j = p_row - i
                if j < 0 or j > m_col:
                    continue
                val += (comb(k - m_col, i) * comb(m_col, j) *
                        a**(k - m_col - i) * c**i * b**(m_col - j) * d**j)
            val *= np.sqrt(comb(k, p_row)) / np.sqrt(comb(k, m_col))
            R[p_row, m_col] = val

    return R


# ================================================================
# Build all 9 representation matrices
# ================================================================
print("="*72)
print("exp337 Step 2: Building irreps, CG, and D_F")
print("="*72)

group = generate_2I()
N = len(group)
print(f"Group order: {N}")

dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
total_dim = sum(dims)  # = 30

print(f"\nComputing Sym^k representations for k=0,...,7...")

# Compute Sym^k matrices for all group elements
sym_reps = {}
for k in range(8):
    sym_reps[k] = [sym_power_matrix(g, k) for g in group]
    # Verify unitarity
    for i, R in enumerate(sym_reps[k]):
        if not np.allclose(R @ R.conj().T, np.eye(k+1), atol=1e-8):
            print(f"  WARNING: Sym^{k} not unitary for element {i}")
            break
    else:
        print(f"  Sym^{k} (dim {k+1}): OK")

# rho_0 through rho_5 are just Sym^0 through Sym^5
reps = {}
for k in range(6):
    reps[k] = sym_reps[k]

# ================================================================
# Extract rho_6, rho_7, rho_8 from Sym^6 and Sym^7
# ================================================================
print("\nExtracting exceptional irreps...")

# rho_6 + rho_8 = Sym^6 (dim 4+3=7)
# Use projection: P_j = (d_j/|G|) * sum_g conj(chi_j(g)) * rho(g)

# First compute character values for each group element
def compute_element_characters(group, sym_reps, dims_list):
    """Compute chi_k(g) for each irrep k and group element g."""
    n_irreps = 9
    n_elts = len(group)

    # Character table (from Step 1 verification)
    # We need chi(g) for each element, not just class representatives
    # Since chi is a class function, we need to identify each element's class
    traces_1 = np.array([np.trace(g).real for g in group])

    # Map elements to classes by trace
    class_trace_vals = sorted(set(np.round(traces_1, 8)))
    class_trace_vals = sorted(set(np.round(traces_1, 6)))

    # For Sym^k, chi_k(g) = Tr(Sym^k(g))
    chi = np.zeros((12, n_elts))  # up to Sym^11
    for k in range(8):
        for i in range(n_elts):
            chi[k, i] = np.trace(sym_reps[k][i]).real

    return chi

chi_sym = compute_element_characters(group, sym_reps, dims)

# Now compute characters of the 9 irreps for each element
chi_irrep = np.zeros((9, N))
for k in range(6):
    chi_irrep[k] = chi_sym[k]

# chi_8 = chi_5 / chi_1 (where chi_1 != 0), and -1 where chi_1 = 0
# From McKay: rho_1 x rho_8 = rho_5, so chi_1 * chi_8 = chi_5
for i in range(N):
    if abs(chi_sym[1, i]) > 1e-8:
        chi_irrep[8, i] = chi_sym[5, i] / chi_sym[1, i]
    else:
        chi_irrep[8, i] = -1.0  # determined from orthogonality in step 1

# chi_6 = chi_Sym6 - chi_8
chi_irrep[6] = chi_sym[6] - chi_irrep[8]

# chi_7 = chi_Sym7 - chi_5
chi_irrep[7] = chi_sym[7] - chi_irrep[5]

# Verify
for i in range(9):
    norm = sum(chi_irrep[i, g]**2 for g in range(N)) / N
    print(f"  rho_{i}: dim={dims[i]}, ||chi||^2 = {norm:.6f} (expect 1.000)")


# ================================================================
# Build explicit representation matrices using projection
# ================================================================
print("\nBuilding explicit irrep matrices via projection...")

def build_irrep_from_projection(group, source_rep, chi_target, d_target, N):
    """Extract irrep from a reducible representation using character projection.

    P = (d/|G|) * sum_g conj(chi(g)) * rho(g)

    This projects onto the isotypic component. Then orthonormalize the image.
    Returns: list of d_target x d_target matrices for each group element.
    """
    d_source = source_rep[0].shape[0]

    # Compute projector
    P = np.zeros((d_source, d_source), dtype=complex)
    for g_idx in range(N):
        P += np.conj(chi_target[g_idx]) * source_rep[g_idx]
    P *= d_target / N

    # P should be a projector: P^2 = P, rank = d_target
    rank = int(round(np.trace(P).real))

    # Find orthonormal basis of range(P) via SVD
    U, S, Vh = np.linalg.svd(P)
    # Keep vectors with singular value ~ 1
    basis_vecs = []
    for j in range(len(S)):
        if S[j] > 0.5:
            basis_vecs.append(U[:, j])

    if len(basis_vecs) != d_target:
        print(f"  WARNING: Expected rank {d_target}, got {len(basis_vecs)} (rank from trace: {rank})")
        # Try eigendecomposition instead
        eigvals, eigvecs = np.linalg.eigh(P)
        basis_vecs = []
        for j in range(len(eigvals)):
            if eigvals[j] > 0.5:
                basis_vecs.append(eigvecs[:, j])

    if len(basis_vecs) != d_target:
        print(f"  FATAL: Cannot extract {d_target}-dim subspace, got {len(basis_vecs)}")
        return None

    # Form the change-of-basis matrix V: columns are the basis vectors
    V = np.column_stack(basis_vecs)  # d_source x d_target

    # The representation in this basis: rho_target(g) = V^dag * rho_source(g) * V
    target_rep = []
    for g_idx in range(N):
        R = V.conj().T @ source_rep[g_idx] @ V
        target_rep.append(R)

    # Verify unitarity
    ok = True
    for g_idx in range(min(10, N)):
        R = target_rep[g_idx]
        if not np.allclose(R @ R.conj().T, np.eye(d_target), atol=1e-7):
            ok = False
            break
    if not ok:
        print(f"  WARNING: Unitarity check failed")

    return target_rep, V


# Extract rho_6 (dim 4) from Sym^6 (dim 7)
print("  Extracting rho_6 from Sym^6...")
reps[6], V6 = build_irrep_from_projection(group, sym_reps[6], chi_irrep[6], 4, N)

# Extract rho_8 (dim 3) from Sym^6 (dim 7)
print("  Extracting rho_8 from Sym^6...")
reps[8], V8 = build_irrep_from_projection(group, sym_reps[6], chi_irrep[8], 3, N)

# Extract rho_7 (dim 2) from Sym^7 (dim 8)
# First verify Sym^7 = rho_5 + rho_7
print("  Extracting rho_7 from Sym^7...")
reps[7], V7 = build_irrep_from_projection(group, sym_reps[7], chi_irrep[7], 2, N)

# Verify all irreps
print("\nVerification of all 9 irreps:")
for k in range(9):
    # Check character matches
    char_ok = True
    for g_idx in range(N):
        tr = np.trace(reps[k][g_idx]).real
        expected = chi_irrep[k, g_idx]
        if abs(tr - expected) > 1e-6:
            char_ok = False
            break
    # Check unitarity
    unit_ok = all(np.allclose(reps[k][g] @ reps[k][g].conj().T,
                               np.eye(dims[k]), atol=1e-7) for g in range(N))
    # Check representation property: rho(g1*g2) = rho(g1)*rho(g2) (spot check)
    rep_ok = True
    for _ in range(50):
        i, j = np.random.randint(N, size=2)
        prod_mat = reps[k][i] @ reps[k][j]
        # Find g_i * g_j in the group
        g_prod = group[i] @ group[j]
        for m in range(N):
            if np.allclose(group[m], g_prod, atol=1e-8):
                if not np.allclose(prod_mat, reps[k][m], atol=1e-6):
                    rep_ok = False
                break
    print(f"  rho_{k} (dim {dims[k]}): char={'OK' if char_ok else 'FAIL'}, "
          f"unitary={'OK' if unit_ok else 'FAIL'}, rep={'OK' if rep_ok else 'FAIL'}")


# ================================================================
# Compute CG intertwiners for all 8 edges
# ================================================================
print("\n" + "="*72)
print("CLEBSCH-GORDAN INTERTWINERS")
print("="*72)

# McKay graph edges (i, j) where rho_j appears in rho_1 x rho_i:
# rho_1 x rho_0 = rho_1           => edge (0, 1)
# rho_1 x rho_1 = rho_0 + rho_2   => edges (1, 0), (1, 2)
# rho_1 x rho_2 = rho_1 + rho_3   => edges (2, 1), (2, 3)
# rho_1 x rho_3 = rho_2 + rho_4   => edges (3, 2), (3, 4)
# rho_1 x rho_4 = rho_3 + rho_5   => edges (4, 3), (4, 5)
# rho_1 x rho_5 = rho_4 + rho_6 + rho_8 => edges (5, 4), (5, 6), (5, 8)
# rho_1 x rho_6 = rho_5 + rho_7   => edges (6, 5), (6, 7)
# rho_1 x rho_7 = rho_6           => edge (7, 6)
# rho_1 x rho_8 = rho_5           => edge (8, 5)

# Unique undirected edges of the McKay graph:
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (5,8), (6,7)]

def compute_cg_intertwiner(group, reps, i, j, dims, N):
    """Compute CG intertwiner Y: C^{d_j} -> C^2 tensor C^{d_i}
    such that Y * rho_j(g) = (rho_1(g) tensor rho_i(g)) * Y for all g.

    Returns Y as a (2*d_i) x d_j matrix.
    """
    d_i = dims[i]
    d_j = dims[j]
    d_prod = 2 * d_i

    # Compute the isotypic projector for rho_j in rho_1 tensor rho_i
    # P = (d_j / |G|) * sum_g conj(chi_j(g)) * (rho_1(g) tensor rho_i(g))
    P = np.zeros((d_prod, d_prod), dtype=complex)
    for g_idx in range(N):
        R1 = reps[1][g_idx]  # 2x2
        Ri = reps[i][g_idx]  # d_i x d_i
        tensor = np.kron(R1, Ri)  # (2*d_i) x (2*d_i)
        chi_j_g = np.trace(reps[j][g_idx]).real
        P += np.conj(chi_j_g) * tensor
    P *= d_j / N

    # P should be idempotent with rank d_j
    rank = int(round(np.trace(P).real))

    # Extract orthonormal basis of range(P)
    eigvals, eigvecs = np.linalg.eigh(P)
    basis = []
    for k_idx in range(len(eigvals)):
        if eigvals[k_idx] > 0.5:
            basis.append(eigvecs[:, k_idx])

    if len(basis) != d_j:
        print(f"  Edge ({i},{j}): Expected rank {d_j}, got {len(basis)} (trace rank: {rank})")
        return None

    # Y: columns are the basis vectors
    Y = np.column_stack(basis)  # (2*d_i) x d_j

    # Verify intertwining: Y * rho_j(g) = (rho_1(g) tensor rho_i(g)) * Y
    ok = True
    for g_idx in range(min(20, N)):
        lhs = Y @ reps[j][g_idx]
        R1 = reps[1][g_idx]
        Ri = reps[i][g_idx]
        rhs = np.kron(R1, Ri) @ Y
        # lhs should equal rhs up to a global phase/unitary on the d_j space
        # Since we chose an arbitrary basis for the projected subspace,
        # rho_j might differ from our chosen reps[j] by a unitary.
        # Check: lhs = rhs * U for some fixed U
        if g_idx == 0:
            # First element: compute U = Y^dag * tensor * Y vs reps[j]
            U_check = Y.conj().T @ np.kron(reps[1][g_idx], reps[i][g_idx]) @ Y
            rj_check = reps[j][g_idx]
            # U_check should be similar to rj_check
        if not np.allclose(np.abs(np.abs(lhs) - np.abs(rhs)), 0, atol=1e-5):
            ok = False

    return Y, rank


print("\nComputing CG intertwiners for all 8 edges...")
cg_maps = {}
for i, j in edges:
    Y, rank = compute_cg_intertwiner(group, reps, i, j, dims, N)
    cg_maps[(i, j)] = Y

    # Compute Frobenius norm
    frob = np.linalg.norm(Y, 'fro')
    print(f"  Edge ({i},{j}): rho_{j}(dim {dims[j]}) in rho_1 x rho_{i}(dim {dims[i]}), "
          f"rank={rank}, ||Y||_F = {frob:.6f} (expect sqrt(d_j)={np.sqrt(dims[j]):.4f})")

    # Verify intertwining property more carefully
    # The action of rho_1(g) x rho_i(g) restricted to the image of Y
    # should be equivalent to rho_j(g).
    # Check: Y^dag * (rho_1(g) x rho_i(g)) * Y should be a d_j x d_j unitary
    test_g = 1  # arbitrary nonidentity element
    R_sub = Y.conj().T @ np.kron(reps[1][test_g], reps[i][test_g]) @ Y
    is_unit = np.allclose(R_sub @ R_sub.conj().T, np.eye(dims[j]), atol=1e-7)
    print(f"    Subspace action unitary: {is_unit}")


# ================================================================
# Normalize CG maps: ||Y||_F = 1/sqrt(2) (framework convention)
# ================================================================
# Actually, let's try different normalizations and see which gives
# the best eigenvalue spectrum. First, the "natural" normalization
# where ||Y||_F = 1 for each edge.

print("\n" + "="*72)
print("BUILDING D_F (30x30 Dirac operator)")
print("="*72)

# D_F block structure: rows/cols indexed by irreps
# Block offsets:
offsets = [0]
for d in dims:
    offsets.append(offsets[-1] + d)
print(f"Block offsets: {offsets}")
print(f"Total dimension: {offsets[-1]} (expect 30)")

def build_D_F(cg_maps, edges, dims, normalization='unit_frob'):
    """Build the 30x30 Dirac operator.

    For each edge (i,j), the block [D_F]_{j,i} is proportional to Y^dag
    and [D_F]_{i,j} = Y (with appropriate normalization).

    Actually, D_F has off-diagonal blocks between the two endpoints of each edge.
    The block from rho_i to rho_j is Y_{ij}^dag (a d_i x d_j matrix? No...)

    Let me think about this more carefully.
    Y_{ij}: C^{d_j} -> C^2 x C^{d_i} is a (2*d_i) x d_j matrix.
    This is NOT directly a block of D_F.

    In the spectral triple framework:
    D_F acts on H_F = direct_sum C^{d_i}.
    The Yukawa coupling Y_{ij} connects irrep i to irrep j.

    The standard prescription: for each edge (i,j) in the McKay graph,
    the off-diagonal block of D_F connecting the d_i and d_j subspaces
    is given by a d_i x d_j matrix derived from the CG intertwiner.

    Specifically, if Y: C^{d_j} -> C^2 x C^{d_i}, we can extract a
    d_i x d_j matrix by "tracing out" the C^2 factor:
    [D_F]_{ij} = sum_a <e_a| Y |.> where e_a are basis vectors of C^2.

    Actually, Y is a (2*d_i) x d_j matrix. We can reshape it as
    Y_reshaped[a, m, n] where a=0,1 (C^2 index), m=0,...,d_i-1, n=0,...,d_j-1.
    Then the "Yukawa matrix" connecting rho_i and rho_j is:
    D_ij = sum_a Y_reshaped[a, :, :] * weight_a

    For the simplest D_F, we trace over the C^2:
    D_ij[m, n] = sum_a Y[a*d_i + m, n]

    Or we can use just one component.
    """
    D = np.zeros((30, 30), dtype=complex)

    for i, j in edges:
        Y = cg_maps[(i, j)]  # (2*d_i) x d_j
        d_i = dims[i]
        d_j = dims[j]

        # Reshape Y as (2, d_i, d_j)
        Y_3d = Y.reshape(2, d_i, d_j)

        # Method 1: Trace over C^2 (sum over the 2-dim index)
        # This gives a d_i x d_j matrix
        D_block = Y_3d[0, :, :] + Y_3d[1, :, :]  # simple sum

        # Normalize
        if normalization == 'unit_frob':
            norm = np.linalg.norm(D_block, 'fro')
            if norm > 1e-10:
                D_block /= norm
        elif normalization == 'sqrt2':
            D_block /= np.sqrt(2)
        elif normalization == 'raw':
            pass

        # Place in D_F
        oi, oj = offsets[i], offsets[j]
        D[oi:oi+d_i, oj:oj+d_j] = D_block
        D[oj:oj+d_j, oi:oi+d_i] = D_block.conj().T

    return D


# Also try: use individual components Y_3d[0] and Y_3d[1] with Pauli structure
def build_D_F_pauli(cg_maps, edges, dims):
    """Build D_F using Pauli structure: D = sigma_x part + sigma_z part, etc."""
    D = np.zeros((30, 30), dtype=complex)

    for i, j in edges:
        Y = cg_maps[(i, j)]
        d_i = dims[i]
        d_j = dims[j]
        Y_3d = Y.reshape(2, d_i, d_j)

        # The two components correspond to the two spinor indices
        # Use them as real and imaginary parts, or as separate blocks
        # Method: D_block = Y_3d[0] (just the "up" component)
        D_block = Y_3d[0, :, :]

        norm = np.linalg.norm(D_block, 'fro')
        if norm > 1e-10:
            D_block /= norm

        oi, oj = offsets[i], offsets[j]
        D[oi:oi+d_i, oj:oj+d_j] = D_block
        D[oj:oj+d_j, oi:oi+d_i] = D_block.conj().T

    return D


def build_D_F_null_space(group, reps, edges, dims, N):
    """Build D_F using the null-space method for CG maps.

    For each edge (i,j): find the d_i x d_j matrix M such that
    M * rho_j(g) = rho_i(g) * M for all g in 2I (intertwining rho_j -> rho_i).

    If rho_i and rho_j are inequivalent, Schur's lemma says M=0.
    But we want the CG map: the intertwiner of rho_j into the
    restriction of rho_1 x rho_i to the rho_j-isotypic component.

    Alternative approach: directly solve for D_F blocks.
    For each edge (i,j), find a d_i x d_j matrix B such that
    there exists a way to embed it as a Yukawa coupling.

    Actually, the simplest correct approach:
    The Yukawa coupling Y_ij is a d_i x d_j matrix satisfying
    the equivariance: rho_i(g) * Y_ij = Y_ij * rho_j(g) for all g.
    This is the space of 2I-equivariant maps from rho_j to rho_i.
    By Schur's lemma, this is 0-dimensional if rho_i != rho_j.

    But that's not right for the Yukawa coupling! The Yukawa couples
    through the Higgs (rho_1). The coupling is:
    Y: rho_j -> rho_1 x rho_i, not rho_j -> rho_i directly.

    For D_F, we need: a map from the rho_j-component of H_F to
    the rho_i-component, mediated by the Higgs vev <H> in rho_1.

    After the Higgs gets a vev v = (v, 0)^T in C^2 (or general direction),
    the Yukawa becomes:
    [D_F]_{ij} = <v| Y_{ij} = v^dag * Y_{ij}

    where v is a 2-component vector and Y_{ij} is (2*d_i) x d_j.
    This gives a d_i x d_j matrix.

    With v = (1, 0)^T: [D_F]_{ij} = Y_3d[0, :, :]
    With v = (0, 1)^T: [D_F]_{ij} = Y_3d[1, :, :]
    With v = (1, 1)^T/sqrt(2): [D_F]_{ij} = (Y_3d[0] + Y_3d[1])/sqrt(2)
    """
    # The Higgs vev direction in C^2
    # In the standard model, v = (0, v)^T or (v, 0)^T depending on convention
    # In SU(2), the vev breaks to U(1). The natural choice is v = (1, 0)^T.

    results = {}
    for v_label, v_higgs in [("(1,0)", np.array([1, 0])),
                               ("(0,1)", np.array([0, 1])),
                               ("(1,1)/sqrt2", np.array([1, 1])/np.sqrt(2))]:
        D = np.zeros((30, 30), dtype=complex)

        for i, j in edges:
            Y = cg_maps[(i, j)]  # (2*d_i) x d_j
            d_i = dims[i]
            d_j = dims[j]
            Y_3d = Y.reshape(2, d_i, d_j)

            # Contract with Higgs vev
            D_block = v_higgs[0] * Y_3d[0] + v_higgs[1] * Y_3d[1]

            oi, oj = offsets[i], offsets[j]
            D[oi:oi+d_i, oj:oj+d_j] = D_block
            D[oj:oj+d_j, oi:oi+d_i] = D_block.conj().T

        results[v_label] = D

    return results


# ================================================================
# Build and analyze D_F
# ================================================================

# Method 1: Trace over C^2
D1 = build_D_F(cg_maps, edges, dims, normalization='unit_frob')
eigs1 = np.sort(np.linalg.eigvalsh(D1))
print("\nMethod 1: Trace over C^2, unit Frobenius norm per edge")
print(f"Eigenvalues: {np.round(eigs1, 6)}")
print(f"Nonzero eigenvalues: {eigs1[np.abs(eigs1) > 1e-8]}")

# Method 2: Raw normalization
D2 = build_D_F(cg_maps, edges, dims, normalization='raw')
eigs2 = np.sort(np.linalg.eigvalsh(D2))
print("\nMethod 2: Raw (no normalization)")
print(f"Eigenvalues: {np.round(eigs2, 6)}")

# Method 3: Higgs vev contraction
D_higgs = build_D_F_null_space(group, reps, edges, dims, N)
for label, D in D_higgs.items():
    eigs = np.sort(np.linalg.eigvalsh(D))
    nz = eigs[np.abs(eigs) > 1e-8]
    print(f"\nMethod 3 (Higgs vev={label}):")
    print(f"  Eigenvalues: {np.round(eigs, 6)}")
    print(f"  Nonzero: {len(nz)}, Zeros: {30-len(nz)}")
    print(f"  Positive eigs: {np.round(nz[nz>0], 6)}")
    if len(nz) > 0:
        print(f"  Tr(D^2) = {np.trace(D @ D).real:.6f}")
        print(f"  Tr(D^4) = {np.trace(D @ D @ D @ D).real:.6f}")


# ================================================================
# Method 4: Direct intertwiner (null-space method)
# ================================================================
print("\n" + "="*72)
print("Method 4: Direct null-space intertwiners rho_j -> rho_i")
print("="*72)

def find_intertwiner_direct(group, reps, i, j, dims, N):
    """Find d_i x d_j matrix M such that rho_i(g) M = M rho_j(g) for all g.
    This uses the averaging (Reynolds) operator:
    M = (1/|G|) sum_g rho_i(g) M_0 rho_j(g)^{-1}
    for random M_0. If <rho_i, rho_j> = 0, result is 0.
    If <rho_i, rho_j> = 1, result is the unique intertwiner (up to scalar).
    """
    d_i = dims[i]
    d_j = dims[j]

    # For edges in McKay graph, rho_j is NOT in rho_i (they're different irreps).
    # The intertwiner rho_j -> rho_i is 0 by Schur.
    # What we want is the CG coefficient: the component of the product
    # representation rho_1 x rho_i restricted to rho_j.

    # Instead, use the "matrix element" approach:
    # For the CG map, we want matrices M_a (a=0,1) such that
    # sum_a rho_1(g)_{ba} * M_a * rho_j(g) = rho_i(g) * M_b
    # This is the tensor-product intertwiner equation.

    # Equivalently: M_a = (d_j/|G|) sum_g [rho_1(g)^{-1}]_{a,b0} * rho_i(g) * M * rho_j(g)^{-1}

    # Let's just use the projection method we already have
    # and extract the d_i x d_j blocks from the Higgs contraction.
    pass


# ================================================================
# Comparison with mass deltas
# ================================================================
print("\n" + "="*72)
print("COMPARISON WITH MASS DELTAS")
print("="*72)

# Mass deltas from exp336b
deltas = {
    'e': 0.000, 'mu': 0.080, 'tau': -0.055,
    'u': -0.004, 'c': 0.247, 't': 0.456,
    'd': -0.402, 's': -0.177, 'b': -0.278
}

# Fermion-to-node assignment (from the framework)
# The 9 nodes of the McKay graph correspond to the 9 fermions
# Assignment: rho_i -> fermion_i
# Standard: rho_0(dim1)->e, rho_1(dim2)->mu?, etc.
# But the assignment is not unique. Let's try the one from exp323.

# From memory: McKay Mass Correspondence (exp323):
# Main chain FORCED: [N_gen, F(3), b1, 0, a1] = [3, 2, 6, 0, 5]
# w=0 fuses rho_7+rho_8
# Solution 3 best: weights = {0, 1, 2, 3, 3, 5, 6, 9}

# The mass exponents are: e(0), mu(11), tau(17), u(3), c(16), t(26), d(5), s(11), b(19)
# Let's try multiple assignments and see which one correlates best with eigenvalues.

print("\nUsing Higgs vev (1,0) D_F eigenvalues:")
D_best = D_higgs["(1,0)"]
eigs_best = np.sort(np.linalg.eigvalsh(D_best))
pos_eigs = np.sort(eigs_best[eigs_best > 1e-8])
neg_eigs = np.sort(eigs_best[eigs_best < -1e-8])[::-1]
zero_count = np.sum(np.abs(eigs_best) < 1e-8)

print(f"  Positive eigenvalues ({len(pos_eigs)}): {np.round(pos_eigs, 6)}")
print(f"  Negative eigenvalues ({len(neg_eigs)}): {np.round(neg_eigs, 6)}")
print(f"  Zero eigenvalues: {zero_count}")

# The eigenvalues of D_F should give some structure.
# Compare with deltas:
delta_vals = [0, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278]
fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']

print(f"\nMass deltas: {[f'{d:+.3f}' for d in delta_vals]}")
print(f"Sorted |deltas|: {sorted(np.abs(delta_vals))}")

# Try: eigenvalues of the ADJACENCY matrix of McKay graph (simpler)
print("\n" + "="*72)
print("ADJACENCY and LAPLACIAN of McKay graph")
print("="*72)

A = np.zeros((9, 9))
for i, j in edges:
    A[i, j] = 1
    A[j, i] = 1

degrees = A.sum(axis=1)
L = np.diag(degrees) - A  # Laplacian
D_adj = np.diag(1/np.sqrt(degrees)) @ A @ np.diag(1/np.sqrt(degrees))  # normalized adj

eigs_A = np.sort(np.linalg.eigvalsh(A))
eigs_L = np.sort(np.linalg.eigvalsh(L))
eigs_D = np.sort(np.linalg.eigvalsh(D_adj))

print(f"Adjacency eigenvalues: {np.round(eigs_A, 6)}")
print(f"Laplacian eigenvalues: {np.round(eigs_L, 6)}")
print(f"Normalized adj eigs:   {np.round(eigs_D, 6)}")

# Weighted Laplacian using dims
D_w = np.diag(np.array(dims, dtype=float))
L_w = D_w @ L  # or various forms
eigs_Lw = np.sort(np.linalg.eigvalsh(L_w))
print(f"Dim-weighted Laplacian: {np.round(eigs_Lw, 6)}")

# Cayley graph Laplacian eigenvalues (from character theory)
# lambda_i = degree - chi_i(generator_sum) / d_i
# For McKay: degree = 2 (each node has 2 connections for the simple chain...)
# Actually, degree varies: rho_5 has degree 3, rho_7 has degree 1, etc.
print(f"\nNode degrees: {dict(zip(range(9), degrees.astype(int)))}")

# Try: D_F^2 eigenvalues
D_sq = D_best @ D_best
eigs_sq = np.sort(np.linalg.eigvalsh(D_sq))
print(f"\nD_F^2 eigenvalues: {np.round(eigs_sq, 6)}")
nz_sq = eigs_sq[eigs_sq > 1e-8]
print(f"Nonzero D_F^2 eigs: {np.round(nz_sq, 6)}")

# Per-node "mass" from D_F^2: sum of D_F^2 entries in the diagonal block
print("\nPer-node Tr(D_F^2) contribution:")
for k in range(9):
    ok = offsets[k]
    block = D_sq[ok:ok+dims[k], ok:ok+dims[k]]
    tr = np.trace(block).real
    print(f"  rho_{k} (dim {dims[k]}): Tr = {tr:.6f}, Tr/dim = {tr/dims[k]:.6f}")

print("\n" + "="*72)
print("FULL D_F EIGENVALUE ANALYSIS")
print("="*72)

# Test ALL three Higgs vev directions
for label, D in D_higgs.items():
    eigs = np.sort(np.linalg.eigvalsh(D))
    pos = np.sort(eigs[eigs > 1e-8])
    print(f"\nHiggs vev = {label}:")
    print(f"  All eigs: {np.round(eigs, 5)}")
    print(f"  Tr(D^2) = {np.trace(D @ D).real:.6f}")
    print(f"  Tr(D^4) = {(np.trace(D @ D @ D @ D)).real:.6f}")

    # Per-node analysis
    Dsq = D @ D
    node_vals = []
    for k in range(9):
        ok = offsets[k]
        block = Dsq[ok:ok+dims[k], ok:ok+dims[k]]
        tr = np.trace(block).real / dims[k]
        node_vals.append(tr)
    print(f"  Per-node Tr(D^2)/d: {[f'{v:.4f}' for v in node_vals]}")

    # Check if per-node values correlate with deltas
    delta_arr = np.array(delta_vals)
    node_arr = np.array(node_vals)
    if np.std(node_arr) > 1e-10:
        corr = np.corrcoef(delta_arr, node_arr)[0, 1]
        print(f"  Correlation with deltas: {corr:.4f}")

print("\n\nDone.")

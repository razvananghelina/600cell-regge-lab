"""
exp419_higgs_coeff_exact.py
============================
GOAL: Derive the coefficient 1 in m_H/m_W = phi - 1*8*alpha.

Strategy:
  1. Build exact CG intertwining maps for all 8 McKay edges of 2I
  2. Normalize each Y_{ij} to ||Y||_F = 1/sqrt(2) (framework Yukawa)
  3. Compute Tr(D_F^2) and Tr(D_F^4) EXACTLY
  4. Use Connes formula: m_H^2/m_W^2 = 8*Tr(D_F^4)/Tr(D_F^2)^2
  5. Check if this gives phi^2 (=> coefficient 1 is DERIVED)

BLIND PROCEDURE: compute first, compare after.

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from scipy.special import comb as binomial
from itertools import product as iprod
from math import sqrt, pi

# Framework constants
a1 = 5
b1 = a1 + 1  # = 6
PHI = (1 + sqrt(a1)) / 2
PHI_PRIME = (1 - sqrt(a1)) / 2
N_GROUP = 120  # = a1!
ALPHA = (4*a1*PHI**4 - sqrt((4*a1*PHI**4)**2 - 8*pi)) / (4*pi)

print("=" * 72)
print("exp419: EXACT Tr(D_F^4) FOR HIGGS COEFFICIENT DERIVATION")
print("=" * 72)
print(f"  a1 = {a1}, phi = {PHI:.10f}, alpha = {ALPHA:.10f}")
print(f"  Target: Tr(D_F^2) = 8, Tr(D_F^4) = 8*phi^2 = {8*PHI**2:.6f}")

# =====================================================================
# SECTION 1: Build 2I (120 elements as unit quaternions)
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 1: Building 2I")
print("=" * 72)


def qmul(q1, q2):
    """Multiply two quaternions."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return (
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    )


def build_2I():
    """Construct all 120 elements of the binary icosahedral group."""
    elements = []

    # Type A: 8 cross-polytope vertices (+-1, 0, 0, 0)
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            elements.append(tuple(v))

    # Type B: 16 tesseract vertices (+-1/2, +-1/2, +-1/2, +-1/2)
    for signs in iprod([0.5, -0.5], repeat=4):
        elements.append(tuple(signs))

    # Type C: 96 golden vertices - even permutations of (0, 1/2, phi/2, 1/(2phi))
    even_perms = [
        (0,1,2,3), (0,2,3,1), (0,3,1,2),
        (1,0,3,2), (1,2,0,3), (1,3,2,0),
        (2,0,1,3), (2,1,3,0), (2,3,0,1),
        (3,0,2,1), (3,1,0,2), (3,2,1,0)
    ]
    bv = [0.0, 0.5, PHI/2, 1/(2*PHI)]

    for perm in even_perms:
        base = [bv[perm[k]] for k in range(4)]
        nonzero = [k for k in range(4) if abs(base[k]) > 1e-10]
        for signs in iprod([1, -1], repeat=len(nonzero)):
            v = list(base)
            for idx, k in enumerate(nonzero):
                v[k] *= signs[idx]
            elements.append(tuple(v))

    # Deduplicate
    unique = []
    for v in elements:
        is_dup = False
        for u in unique:
            if sum((v[k] - u[k])**2 for k in range(4)) < 1e-10:
                is_dup = True
                break
        if not is_dup:
            unique.append(v)

    return unique


elements_2I = build_2I()
print(f"|2I| = {len(elements_2I)} (expected: {N_GROUP})")
assert len(elements_2I) == N_GROUP, f"Got {len(elements_2I)} elements!"


def find_element(q):
    """Find the index of quaternion q in the group."""
    for k in range(N_GROUP):
        if sum((q[m] - elements_2I[k][m])**2 for m in range(4)) < 1e-8:
            return k
    return -1


# =====================================================================
# SECTION 2: Build all 9 irreducible representations
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 2: Representations")
print("=" * 72)


def quat_to_su2(q):
    """Convert quaternion (w,x,y,z) to SU(2) matrix."""
    w, x, y, z = q
    a = complex(w, x)
    b = complex(y, z)
    return np.array([[a, b], [-np.conj(b), np.conj(a)]])


def galois_real(x):
    """Apply Galois conjugation phi -> phi' to a real number."""
    for b_num in range(-6, 7):
        b = b_num / 2.0
        a_val = x - b * PHI
        a2 = round(2 * a_val)
        if abs(a_val - a2/2) < 1e-8:
            return (a2/2 + b) - b * PHI
    return x  # fallback


def galois_matrix(M):
    """Apply Galois conjugation phi -> phi' to all entries of a matrix."""
    result = np.empty_like(M)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            z = M[i, j]
            result[i, j] = galois_real(z.real) + 1j * galois_real(z.imag)
    return result


def spin_j_matrix(U, j):
    """Compute spin-j representation matrix of SU(2) element U."""
    two_j = int(round(2 * j))
    dim = two_j + 1
    a, b = U[0, 0], U[0, 1]
    c, d = U[1, 0], U[1, 1]

    D = np.zeros((dim, dim), dtype=complex)

    for m_idx in range(dim):
        m = j - m_idx
        p = int(round(j + m))
        q_val = int(round(j - m))

        for mp_idx in range(dim):
            mp = j - mp_idx
            pp = int(round(j + mp))

            val = 0.0 + 0.0j
            for s in range(p + 1):
                t = pp - s
                if 0 <= t <= q_val:
                    coeff = (binomial(p, s, exact=True) *
                             binomial(q_val, t, exact=True) *
                             a**s * c**(p - s) * b**t * d**(q_val - t))
                    val += coeff
            D[mp_idx, m_idx] = val

    return D


# dims in our 0-indexed node ordering matching McKay graph
# Node: 0   1   2   3   4   5   6   7   8
# dim:  1   2   3   4   5   6   4   2   3
irrep_dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

# Build rep matrices
# Node 0 = trivial (dim 1)
# Node 1 = std (dim 2, spin-1/2)
# Node 2 = Sym^2(std) (dim 3, spin-1)
# Node 3 = Sym^3(std) (dim 4, spin-3/2)
# Node 4 = Sym^4(std) (dim 5, spin-2)
# Node 5 = Sym^5(std) (dim 6, spin-5/2)
# Node 6 = std tensor galois (dim 4)
# Node 7 = galois (dim 2, spin-1/2 Galois)
# Node 8 = Sym^2(galois) (dim 3, spin-1 Galois)

# WAIT - need to match the correct McKay graph ordering.
# From MEMORY: edges (0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(5,8)
# With dims [1,2,3,4,5,6,4,2,3]
# Node 5 has degree 3 (branch node), connects to 4, 6, 8
# Long leg: 0-1-2-3-4-5 (dims 1,2,3,4,5,6)
# Medium arm: 5-6-7 (dims 6,4,2)
# Short arm: 5-8 (dims 6,3)

# The representations:
# spin-j representation for j = 0, 1/2, 1, 3/2, 2, 5/2 gives dims 1,2,3,4,5,6
# That maps to nodes 0,1,2,3,4,5
# Node 6 (dim 4): This is std x galois_std (spin-1/2 x galois spin-1/2)
# Node 7 (dim 2): galois_std (Galois conjugate of spin-1/2)
# Node 8 (dim 3): Sym^2(galois_std) (Galois conjugate of spin-1)

print("Building all 120 representation matrices...")

reps = [[] for _ in range(9)]  # 0-indexed nodes

for g_idx, q in enumerate(elements_2I):
    U = quat_to_su2(q)
    U_gal = galois_matrix(U)

    # Node 0: trivial (dim 1)
    reps[0].append(np.array([[1.0 + 0j]]))

    # Node 1: std = spin-1/2 (dim 2)
    reps[1].append(U.copy())

    # Node 2: Sym^2(std) = spin-1 (dim 3)
    reps[2].append(spin_j_matrix(U, 1))

    # Node 3: Sym^3(std) = spin-3/2 (dim 4)
    reps[3].append(spin_j_matrix(U, 1.5))

    # Node 4: Sym^4(std) = spin-2 (dim 5)
    reps[4].append(spin_j_matrix(U, 2))

    # Node 5: Sym^5(std) = spin-5/2 (dim 6)
    reps[5].append(spin_j_matrix(U, 2.5))

    # Node 7: galois_std (dim 2)
    reps[7].append(U_gal.copy())

    # Node 8: Sym^2(galois_std) = galois spin-1 (dim 3)
    reps[8].append(spin_j_matrix(U_gal, 1))

    # Node 6: std tensor galois_std (dim 4)
    reps[6].append(np.kron(U, U_gal))

# Verify characters at identity
print("  Character at identity (should be dims):")
for k in range(9):
    chi = np.trace(reps[k][0]).real
    print(f"    Node {k}: chi = {chi:.1f} (dim {irrep_dims[k]})")

# Compute characters for all elements
characters = []
for k in range(9):
    chars = np.array([np.trace(reps[k][g_idx]) for g_idx in range(N_GROUP)])
    characters.append(chars)

# Verify orthogonality
print("\n  Character orthogonality (sample):")
for i in range(9):
    for j in range(i, min(i+2, 9)):
        inner = np.sum(np.conj(characters[i]) * characters[j]).real / N_GROUP
        expected = 1.0 if i == j else 0.0
        print(f"    <chi_{i}, chi_{j}> = {inner:.6f} (expected {expected:.0f})")

# =====================================================================
# SECTION 3: Verify McKay graph
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 3: McKay Graph Verification")
print("=" * 72)

# McKay graph: decompose rho_1 tensor rho_k for each k
# rho_1 = node 1 (std, dim 2)
print("Computing McKay adjacency from tensor products with rho_1...")

mckay_adj = np.zeros((9, 9), dtype=int)
for k in range(9):
    product_chars = characters[1] * characters[k]
    for j in range(9):
        mult = np.sum(np.conj(characters[j]) * product_chars).real / N_GROUP
        mult_int = int(round(mult))
        if mult_int > 0:
            mckay_adj[k, j] = mult_int

expected_edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
print("  Expected edges: ", expected_edges)
print("  Computed adjacency:")
actual_edges = []
for i in range(9):
    for j in range(i+1, 9):
        if mckay_adj[i, j] > 0:
            actual_edges.append((i, j))
print("  Actual edges: ", actual_edges)

# Check McKay null condition: 2*d_i = sum of neighbor dims
print("\n  McKay null condition (2*d_i = sum neighbor dims):")
for i in range(9):
    nbr_sum = sum(mckay_adj[i, j] * irrep_dims[j] for j in range(9))
    ok = (nbr_sum == 2 * irrep_dims[i])
    print(f"    Node {i} (d={irrep_dims[i]}): 2*{irrep_dims[i]}={2*irrep_dims[i]}, "
          f"sum={nbr_sum} {'OK' if ok else 'FAIL'}")

# Use actual edges from here on
mckay_edges = actual_edges
print(f"\n  Total edges: {len(mckay_edges)} (expected: 8)")

# =====================================================================
# SECTION 4: Compute CG intertwining maps for all edges
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 4: CG Intertwining Maps")
print("=" * 72)

# For each edge (i,j), compute CG map C: V_j -> V_1 tensor V_i
# satisfying: (rho_1(g) tensor rho_i(g)) C = C rho_j(g) for all g
# C has shape (2*d_i, d_j)

cg_maps = {}  # (i, j) -> C matrix of shape (2*d_i, d_j)

for node_i, node_j in mckay_edges:
    # For each direction of the edge
    for source, target in [(node_i, node_j), (node_j, node_i)]:
        d_s = irrep_dims[source]
        d_t = irrep_dims[target]

        # Check multiplicity of target in rho_1 tensor source
        product_chars = characters[1] * characters[source]
        mult = np.sum(np.conj(characters[target]) * product_chars).real / N_GROUP
        mult_int = int(round(mult))

        if mult_int < 1:
            continue

        dim_product = 2 * d_s  # dim(V_1 tensor V_source) = 2 * d_source
        n_vars = dim_product * d_t

        # Build intertwining equation: (rho_1(g) x rho_s(g)) C = C rho_t(g)
        # Vectorized: [kron(rho_1(g) x rho_s(g), I_dt) - kron(I_2ds, rho_t(g)^T)] vec(C) = 0
        equations = []
        # Use a subset of group elements (sufficient for determining null space)
        for g_idx in range(N_GROUP):
            kron_g = np.kron(reps[1][g_idx], reps[source][g_idx])
            rho_t_g = reps[target][g_idx]
            A_g = (np.kron(kron_g, np.eye(d_t, dtype=complex)) -
                   np.kron(np.eye(dim_product, dtype=complex), rho_t_g.T))
            equations.append(A_g)

        A = np.vstack(equations)
        U_svd, S_svd, Vh_svd = np.linalg.svd(A, full_matrices=True)

        # Find null space
        if S_svd[0] > 1e-15:
            null_dim = int(np.sum(S_svd / S_svd[0] < 1e-8))
        else:
            null_dim = n_vars

        if null_dim >= 1:
            C_vec = Vh_svd[-1].conj()
            C_map = C_vec.reshape(dim_product, d_t)

            # Phase fix: largest entry real and positive
            max_idx = np.argmax(np.abs(C_map))
            max_val = C_map.flat[max_idx]
            if abs(max_val) > 1e-15:
                C_map *= np.abs(max_val) / max_val

            # Verify intertwining
            max_err = 0
            for g_idx in range(0, N_GROUP, 10):  # sample
                kron_g = np.kron(reps[1][g_idx], reps[source][g_idx])
                lhs = kron_g @ C_map
                rhs = C_map @ reps[target][g_idx]
                err = np.max(np.abs(lhs - rhs))
                max_err = max(max_err, err)

            cg_maps[(source, target)] = C_map
            print(f"  CG({source}->{target}): shape {C_map.shape}, "
                  f"null_dim={null_dim}, err={max_err:.2e}")

print(f"\n  Total CG maps: {len(cg_maps)}")

# =====================================================================
# SECTION 5: Build D_F with ||Y||_F = 1/sqrt(2) normalization
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 5: Build D_F with Universal Yukawa Norm")
print("=" * 72)

total_dim = sum(irrep_dims)  # = 30
print(f"  H_F dimension: {total_dim}")

# Block offsets
offsets = [0]
for d in irrep_dims:
    offsets.append(offsets[-1] + d)

# For each McKay edge (i,j), we need a d_i x d_j Yukawa matrix Y_{ij}
# Extract from CG map: the CG map C: V_j -> V_1 tensor V_i has shape (2*d_i, d_j)
# The Yukawa is the "Higgs vev slice": first d_i rows of C (vev direction (1,0))
# Then normalize to ||Y||_F = 1/sqrt(2)

D_F = np.zeros((total_dim, total_dim), dtype=complex)

print("\n  Edge Yukawa matrices (normalized to ||Y||_F = 1/sqrt(2)):")
for node_i, node_j in mckay_edges:
    di = irrep_dims[node_i]
    dj = irrep_dims[node_j]

    # Get CG map for this edge: C: V_j -> V_1 tensor V_i
    if (node_i, node_j) in cg_maps:
        C = cg_maps[(node_i, node_j)]
        # Y_{ij} = first d_i rows of C (Higgs vev = (1,0) in V_1)
        Y_raw = C[:di, :]  # shape (d_i, d_j)
    elif (node_j, node_i) in cg_maps:
        C = cg_maps[(node_j, node_i)]
        # Y_{ji} = first d_j rows of C, shape (d_j, d_i)
        # Y_{ij} = Y_{ji}^dagger, shape (d_i, d_j)
        Y_raw = C[:dj, :].conj().T
    else:
        print(f"    WARNING: No CG map for edge ({node_i},{node_j})!")
        continue

    # Normalize to ||Y||_F = 1/sqrt(2)
    norm = np.linalg.norm(Y_raw, 'fro')
    if norm > 1e-15:
        Y = Y_raw * (1.0 / (sqrt(2) * norm))
    else:
        print(f"    WARNING: Zero Y for edge ({node_i},{node_j})!")
        continue

    Y_norm = np.linalg.norm(Y, 'fro')

    # Singular values
    svs = np.linalg.svd(Y, compute_uv=False)
    svs_nonzero = svs[svs > 1e-10]
    rank = len(svs_nonzero)

    print(f"    Edge ({node_i},{node_j}): shape ({di},{dj}), "
          f"||Y||_F = {Y_norm:.6f}, rank = {rank}, "
          f"sv = [{', '.join(f'{s:.4f}' for s in svs_nonzero)}]")

    # Place in D_F
    ri, rj = offsets[node_i], offsets[node_j]
    D_F[ri:ri+di, rj:rj+dj] = Y
    D_F[rj:rj+dj, ri:ri+di] = Y.conj().T

# Make exactly Hermitian
D_F = (D_F + D_F.conj().T) / 2

print(f"\n  D_F Hermitian: {np.allclose(D_F, D_F.conj().T)}")

# =====================================================================
# SECTION 6: Compute Tr(D_F^2) and Tr(D_F^4)
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 6: KEY RESULTS - Tr(D_F^2) and Tr(D_F^4)")
print("=" * 72)

D2 = D_F @ D_F
D4 = D2 @ D2

tr_DF2 = np.trace(D2).real
tr_DF4 = np.trace(D4).real

print(f"\n  Tr(D_F^2) = {tr_DF2:.10f}")
print(f"    Expected: 8 = rank(E8)")
print(f"    Match: {abs(tr_DF2 - 8.0) < 1e-6}")

print(f"\n  Tr(D_F^4) = {tr_DF4:.10f}")
print(f"    Target for Connes (b = 8*phi^2): {8*PHI**2:.10f}")
print(f"    Ratio Tr(D_F^4)/(8*phi^2) = {tr_DF4/(8*PHI**2):.10f}")

# Connes formula: m_H^2/m_W^2 = 8*b/a^2 = 8*Tr(D_F^4)/Tr(D_F^2)^2
connes_ratio = 8 * tr_DF4 / tr_DF2**2
connes_mH_over_mW = sqrt(connes_ratio)

print(f"\n  Connes: m_H^2/m_W^2 = 8*Tr(D_F^4)/Tr(D_F^2)^2 = {connes_ratio:.10f}")
print(f"    = {connes_mH_over_mW:.10f}")
print(f"    phi^2 = {PHI**2:.10f}")
print(f"    phi   = {PHI:.10f}")
print(f"    (phi-8*alpha)^2 = {(PHI-8*ALPHA)**2:.10f}")

# Also compute other trace ratios
print(f"\n  Other trace ratios:")
print(f"    Tr(D_F^4)/Tr(D_F^2) = {tr_DF4/tr_DF2:.10f}")
print(f"    2*Tr(D_F^4)/Tr(D_F^2) = {2*tr_DF4/tr_DF2:.10f}")
print(f"    Tr(D_F^4)/Tr(D_F^2)^2 = {tr_DF4/tr_DF2**2:.10f}")
print(f"    b1 = 6: {b1}")
print(f"    phi^2/8 = {PHI**2/8:.10f}")

# D_F eigenvalues
eigs = np.sort(np.linalg.eigvalsh(D_F))[::-1]
print(f"\n  D_F eigenvalues:")
for i, ev in enumerate(eigs):
    print(f"    lambda_{i+1:2d} = {ev:12.8f}")

# =====================================================================
# SECTION 7: Alternative - Use BOTH Higgs directions
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 7: Full CG D_F (both Higgs directions)")
print("=" * 72)

# The CG map C has shape (2*d_i, d_j).
# Direction 1: rows 0..d_i-1 (Higgs vev (1,0))
# Direction 2: rows d_i..2*d_i-1 (Higgs vev (0,1))
# The FULL D_F should include BOTH, or we should understand
# why only one direction matters.

# Actually, in the Connes framework, D_F includes the FULL CG structure.
# The Higgs vev selects a direction, but the TRACES of D_F^2, D_F^4
# involve the full operator.

# Let me build D_F using the FULL CG map (not just one vev slice).
# D_F as 30x30 with blocks being the FULL CG projectors.

# Actually, the issue is that the CG map C: V_j -> V_1 tensor V_i
# is a (2*d_i) x d_j matrix. The Yukawa Y_{ij} that goes into D_F
# is a d_i x d_j matrix, obtained by contracting with one slot of V_1.

# In the Connes framework, D_F includes ALL Yukawa couplings.
# For the mass formula, the relevant traces are:
# a = Tr(Y^dag Y) where Y is the full Yukawa MATRIX on H_F
# b = Tr((Y^dag Y)^2)

# Y^dag Y = D_F^2 (diagonal blocks), so a = Tr(D_F^2)
# (Y^dag Y)^2 involves paths through the graph.

# The key question: should we use ONE Higgs direction or BOTH?
# In the SM, there's ONE Higgs doublet H = (H+, H0).
# The vev is <H> = (0, v/sqrt(2)). So one direction.
# But before EWSB, both directions contribute to the traces.

# Let's compute with the FULL CG projector approach.
# For each edge (i,j), the FULL Yukawa is the CG map C itself,
# viewed as an operator on V_1 tensor V_i -> V_j.
# In the Connes framework: Y_{ij}^dag Y_{ij} = C^dag C is a d_j x d_j matrix.
# And Y_{ij} Y_{ij}^dag = C C^dag is a (2*d_i) x (2*d_i) matrix.

# For traces: Tr(Y^dag Y) = ||C||_F^2
# With our normalization: ||Y (one slice)||_F = 1/sqrt(2)
# The FULL CG has ||C||_F^2 = ||slice1||_F^2 + ||slice2||_F^2

# Let me compute ||C||_F for each edge to understand the structure.
print("  Full CG map norms:")
for node_i, node_j in mckay_edges:
    di = irrep_dims[node_i]
    dj = irrep_dims[node_j]

    if (node_i, node_j) in cg_maps:
        C = cg_maps[(node_i, node_j)]
        norm_full = np.linalg.norm(C, 'fro')
        slice1 = C[:di, :]
        slice2 = C[di:2*di, :]
        norm1 = np.linalg.norm(slice1, 'fro')
        norm2 = np.linalg.norm(slice2, 'fro')
        print(f"    Edge ({node_i},{node_j}): ||C||_F = {norm_full:.6f}, "
              f"||slice1|| = {norm1:.6f}, ||slice2|| = {norm2:.6f}, "
              f"ratio = {norm1/norm2:.6f}")

# =====================================================================
# SECTION 8: Analytical computation of Tr(D_F^4)
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 8: Analytical Tr(D_F^4) Structure")
print("=" * 72)

# On the McKay TREE, Tr(D_F^4) = sum_i Tr((D_F^2)_{ii}^2) + off-diagonal
# (D_F^2)_{ii} = sum_{j~i} Y_{ij} Y_{ij}^dag
# For normalized Y with ||Y||_F = 1/sqrt(2):
# Tr(Y Y^dag) = ||Y||_F^2 = 1/2

# The DIAGONAL part of Tr(D_F^4):
# Tr((sum_j Y_{ij} Y_{ij}^dag)^2) = sum_j Tr((Y_{ij} Y_{ij}^dag)^2) + cross terms
# Cross terms: sum_{j!=k} Tr(Y_{ij} Y_{ij}^dag Y_{ik} Y_{ik}^dag)

# For CG-derived Y maps: are the cross terms zero?
# Y_{ij} Y_{ij}^dag is the projection onto the j-component in V_1 tensor V_i
# For different j, these projections are ORTHOGONAL (direct sum decomposition)
# So Y_{ij} Y_{ij}^dag * Y_{ik} Y_{ik}^dag = 0 for j != k

# Verify this numerically:
print("  Cross-term orthogonality check (CG projectors):")
neighbors = {i: [] for i in range(9)}
for i, j in mckay_edges:
    neighbors[i].append(j)
    neighbors[j].append(i)

max_cross = 0
for i in range(9):
    nbrs = neighbors[i]
    if len(nbrs) < 2:
        continue
    for a_idx in range(len(nbrs)):
        for b_idx in range(a_idx+1, len(nbrs)):
            j, k = nbrs[a_idx], nbrs[b_idx]
            di = irrep_dims[i]

            # Get Y_{ij} from D_F
            Y_ij = D_F[offsets[i]:offsets[i]+di, offsets[j]:offsets[j]+irrep_dims[j]]
            Y_ik = D_F[offsets[i]:offsets[i]+di, offsets[k]:offsets[k]+irrep_dims[k]]

            cross = np.trace(Y_ij @ Y_ij.conj().T @ Y_ik @ Y_ik.conj().T)
            max_cross = max(max_cross, abs(cross))
            if abs(cross) > 1e-8:
                print(f"    Node {i}, nbrs ({j},{k}): cross = {cross:.6e} NON-ZERO!")

print(f"  Max cross term: {max_cross:.6e}")
if max_cross < 1e-8:
    print("  ALL cross terms vanish! (CG orthogonality confirmed)")
else:
    print("  WARNING: Cross terms are NOT zero!")

# With orthogonality, Tr(D_F^4) = sum_i sum_{j~i} Tr((Y_{ij} Y_{ij}^dag)^2)
print("\n  Per-edge contributions to Tr(D_F^4):")
tr4_edge_total = 0
for i in range(9):
    di = irrep_dims[i]
    for j in neighbors[i]:
        dj = irrep_dims[j]
        Y_ij = D_F[offsets[i]:offsets[i]+di, offsets[j]:offsets[j]+dj]
        # Tr((Y Y^dag)^2) where Y is d_i x d_j
        P = Y_ij @ Y_ij.conj().T  # d_i x d_i
        contrib = np.trace(P @ P).real
        tr4_edge_total += contrib

        # What is the rank and singular value structure?
        svs = np.linalg.svd(Y_ij, compute_uv=False)
        svs_nz = svs[svs > 1e-10]
        rank = len(svs_nz)
        # Tr((YY^dag)^2) = sum sigma^4
        sum_s4 = sum(s**4 for s in svs_nz)

        print(f"    ({i}->{j}): rank={rank}, sv={[f'{s:.6f}' for s in svs_nz]}, "
              f"Tr((YY^dag)^2) = {contrib:.8f}, 1/(4*rank) = {1/(4*rank):.8f}")

print(f"\n  Sum of edge contributions: {tr4_edge_total:.10f}")
print(f"  Direct Tr(D_F^4): {tr_DF4:.10f}")
print(f"  Match: {abs(tr4_edge_total - tr_DF4) < 1e-6}")

# Off-diagonal check
print(f"\n  Off-diagonal contribution: {tr_DF4 - tr4_edge_total:.10f}")

# =====================================================================
# SECTION 9: Compare with Connes prediction
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 9: COMPARISON WITH CONNES PREDICTION")
print("=" * 72)

# Connes: m_H^2/m_W^2 = 8*b/a^2
a_val = tr_DF2
b_val = tr_DF4

connes_mH2_mW2 = 8 * b_val / a_val**2
connes_mH_mW = sqrt(connes_mH2_mW2) if connes_mH2_mW2 > 0 else 0

print(f"  a = Tr(D_F^2) = {a_val:.10f}")
print(f"  b = Tr(D_F^4) = {b_val:.10f}")
print(f"  Connes: m_H^2/m_W^2 = 8*b/a^2 = {connes_mH2_mW2:.10f}")
print(f"  Connes: m_H/m_W = {connes_mH_mW:.10f}")
print(f"  phi = {PHI:.10f}")
print(f"  phi - 8*alpha = {PHI - 8*ALPHA:.10f}")
print(f"  Connes/phi = {connes_mH_mW/PHI:.10f}")

# What would b need to be for Connes to give phi?
b_for_phi = a_val**2 * PHI**2 / 8
print(f"\n  For Connes = phi: b = {b_for_phi:.10f} = 8*phi^2 = {8*PHI**2:.10f}")
print(f"  Actual b = {b_val:.10f}")
print(f"  Ratio b/b_target = {b_val/b_for_phi:.10f}")

# What does the ACTUAL Connes ratio give in GeV?
m_W_exp = 80.3692  # GeV
m_H_exp = 125.25   # GeV
m_H_connes = m_W_exp * connes_mH_mW
print(f"\n  m_H(Connes) = {m_H_connes:.2f} GeV")
print(f"  m_H(phi-8a) = {m_W_exp * (PHI - 8*ALPHA):.2f} GeV")
print(f"  m_H(exp) = {m_H_exp:.2f} GeV")

# =====================================================================
# SECTION 10: Alternative Connes normalizations
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 10: Alternative Trace Ratios")
print("=" * 72)

# Maybe the Connes formula is not 8*b/a^2 but some other combination
# Test various ratios to see which gives phi^2 or (phi-8alpha)^2
ratios_to_test = [
    ("8*b/a^2", 8*b_val/a_val**2),
    ("b/a", b_val/a_val),
    ("2*b/a", 2*b_val/a_val),
    ("b/a^2", b_val/a_val**2),
    ("4*b/a^2", 4*b_val/a_val**2),
    ("16*b/a^2", 16*b_val/a_val**2),
    ("b/(a/2)", b_val/(a_val/2)),
    ("2*b/a^2", 2*b_val/a_val**2),
]

print(f"  {'Ratio':.<30} = {'Value':>12}  {'sqrt':>12}  {'phi err':>10}  {'(phi-8a) err':>12}")
for name, val in ratios_to_test:
    sq = sqrt(val) if val > 0 else 0
    err_phi = abs(sq - PHI)/PHI*100
    err_corr = abs(sq - (PHI-8*ALPHA))/(PHI-8*ALPHA)*100
    print(f"  {name:.<30} = {val:12.8f}  {sq:12.8f}  {err_phi:8.4f}%  {err_corr:10.4f}%")

# =====================================================================
# SECTION 11: Honest Assessment
# =====================================================================
print("\n" + "=" * 72)
print("SECTION 11: HONEST ASSESSMENT")
print("=" * 72)

print(f"""
KEY RESULTS:
  Tr(D_F^2) = {tr_DF2:.6f} (expected 8: {'MATCH' if abs(tr_DF2-8)<1e-6 else 'MISMATCH'})
  Tr(D_F^4) = {tr_DF4:.6f}
  Target 8*phi^2 = {8*PHI**2:.6f}

  Connes m_H/m_W = sqrt(8*b/a^2) = {connes_mH_mW:.6f}
  Framework m_H/m_W = phi - 8*alpha = {PHI-8*ALPHA:.6f}
  Experiment m_H/m_W = {m_H_exp/m_W_exp:.6f}

  The Connes formula {'DOES' if abs(connes_mH_mW - PHI) < 0.01 else 'does NOT'} reproduce phi at tree level.
""")

print("=" * 72)
print("EXPERIMENT COMPLETE")
print("=" * 72)

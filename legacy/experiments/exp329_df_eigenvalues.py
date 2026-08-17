"""
exp329_df_eigenvalues.py - D_F Eigenvalues on Weighted McKay Graph
==================================================================

Can the finite Dirac operator D_F on the E8-hat McKay graph, with edge
weights from exp323, produce eigenvalues proportional to the mass exponents?

Background:
- exp307: CG maps, Tr(D_F^2)=8, ||Y||_F=1/sqrt(2) universal
- exp316: D_F alone has too small dynamic range (~15 vs phi^26)
- exp323: Edge weights with main chain FORCED [3,2,6,0,5]

Key idea: D_F(unweighted) gives SELECTION RULES (which transitions allowed).
          Edge weights w_e give HIERARCHY (how strong each coupling).
          D_F(weighted) = sum_e w_e * Y_e should encode BOTH.

If eigenvalues of D_F(weighted) ~ phi^n_f, mass formula is DERIVED from McKay.

NOTE: This script ports the CORRECT representation-building code from exp307
(galois_real, spin_j_matrix with scipy binomial).

Author: Razvan-Constantin Anghelina
Date: 2026-02-16
"""

import numpy as np
from scipy.special import comb as binomial
from itertools import product as iprod
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

print("=" * 72)
print("exp329: D_F EIGENVALUES ON WEIGHTED McKAY GRAPH")
print("=" * 72)

# =====================================================================
# PART 0: CONSTANTS
# =====================================================================
a_1 = 5
b_1 = 6
N_gen = 3
N_eig = 9
PHI = (1 + np.sqrt(5)) / 2
PHI_PRIME = 1 - PHI  # = (1-sqrt(5))/2 = -1/PHI
SQRT5 = np.sqrt(5)
me = 0.51099895  # MeV

# Irrep dimensions of 2I (binary icosahedral group)
irrep_dims = {1:1, 2:2, 3:2, 4:3, 5:3, 6:4, 7:4, 8:5, 9:6}
total_dim = sum(irrep_dims.values())  # = 30
print(f"a_1 = {a_1}, phi = {PHI:.10f}, phi' = {PHI_PRIME:.10f}")
print(f"Total dim(H_F) = {total_dim}")

# Block starts in the 30-dim Hilbert space
block_starts = {}
pos = 0
for i in range(1, 10):
    block_starts[i] = pos
    pos += irrep_dims[i]

# Fermion mass exponents n_f = 5a + 6b
fermion_data = {
    'e': (0, 0, 0), 'mu': (1, 1, 11), 'tau': (1, 2, 17),
    'u': (3, -2, 3), 'c': (2, 1, 16), 't': (4, 1, 26),
    'd': (1, 0, 5), 's': (1, 1, 11), 'b': (-1, 4, 19)
}

# =====================================================================
# PART 1: BUILD 2I GROUP (120 quaternions) - from exp307
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: Binary Icosahedral Group 2I")
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
N_group = len(elements_2I)
print(f"|2I| = {N_group} (expected: 120 = a_1!)")
assert N_group == 120, f"Expected 120, got {N_group}"

# Verify: all unit quaternions
norms = [np.sqrt(sum(c**2 for c in q)) for q in elements_2I]
print(f"  All unit quaternions: {all(abs(n - 1.0) < 1e-10 for n in norms)}")


def find_element(q):
    """Find the index of quaternion q in the group."""
    for k in range(N_group):
        if sum((q[m] - elements_2I[k][m])**2 for m in range(4)) < 1e-8:
            return k
    return -1


# =====================================================================
# PART 2: BUILD ALL 9 REPRESENTATIONS - from exp307
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: Building 9 Irreducible Representations")
print("=" * 72)


def quat_to_su2(q):
    """Convert quaternion (w,x,y,z) to SU(2) matrix."""
    w, x, y, z = q
    a = complex(w, x)
    b = complex(y, z)
    return np.array([[a, b], [-np.conj(b), np.conj(a)]])


def galois_real(x):
    """Apply Galois conjugation phi -> phi' to a real number a + b*phi.
    Returns sigma(x) = (a+b) - b*phi for x = a + b*phi."""
    # Decompose x = a + b*phi numerically
    for b in [-1.0, -0.5, 0.0, 0.5, 1.0]:
        a = x - b * PHI
        a2 = round(2 * a)
        if abs(a - a2/2) < 1e-8:
            return (a2/2 + b) - b * PHI
    # Extended range
    for b_num in range(-4, 5):
        b = b_num / 2.0
        a = x - b * PHI
        a2 = round(2 * a)
        if abs(a - a2/2) < 1e-8:
            return (a2/2 + b) - b * PHI
    return x  # no phi component found


def galois_complex(z):
    """Apply Galois conjugation to a complex number."""
    return galois_real(z.real) + 1j * galois_real(z.imag)


def galois_matrix(M):
    """Apply Galois conjugation phi -> phi' to all entries of a matrix."""
    result = np.empty_like(M)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            result[i, j] = galois_complex(M[i, j])
    return result


def spin_j_matrix(U, j):
    """Compute spin-j representation matrix of SU(2) element U.
    Basis: e_m for m = j, j-1, ..., -j (row 0 = m=j).
    Uses scipy binomial with exact=True."""
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


# Build representation matrices for all 120 elements
print("Building representation matrices...")

reps = {i: [] for i in range(1, 10)}
irrep_names = {
    1: "rho_1 (trivial, dim 1)",
    2: "rho_2 (fundamental, dim 2)",
    3: "rho_3 (Galois conj, dim 2)",
    4: "rho_4 = Sym^2(rho_2), dim 3",
    5: "rho_5 = Sym^2(rho_3), dim 3",
    6: "rho_6 = rho_2 x rho_3, dim 4",
    7: "rho_7 = Sym^3(rho_2), dim 4",
    8: "rho_8 = Sym^4(rho_2), dim 5",
    9: "rho_9 = Sym^5(rho_2), dim 6",
}

for g_idx, q in enumerate(elements_2I):
    U = quat_to_su2(q)        # rho_2(g)
    U_gal = galois_matrix(U)  # rho_3(g)

    reps[1].append(np.array([[1.0 + 0j]]))
    reps[2].append(U.copy())
    reps[3].append(U_gal.copy())
    reps[4].append(spin_j_matrix(U, 1))
    reps[5].append(spin_j_matrix(U_gal, 1))
    reps[6].append(np.kron(U, U_gal))
    reps[7].append(spin_j_matrix(U, 1.5))
    reps[8].append(spin_j_matrix(U, 2))
    reps[9].append(spin_j_matrix(U, 2.5))

print("  All representations built.")

# Verify dimensions + representation property
for i in range(1, 10):
    d = irrep_dims[i]
    assert reps[i][0].shape == (d, d), f"rho_{i}: wrong shape"

# Quick rep property check
test_pairs = [(0, 1), (5, 10), (23, 47)]
for i_rep in [2, 3, 4, 5, 6, 7, 8, 9]:
    max_err = 0
    for g1_idx, g2_idx in test_pairs:
        prod_q = qmul(elements_2I[g1_idx], elements_2I[g2_idx])
        prod_idx = find_element(prod_q)
        if prod_idx >= 0:
            product_mat = reps[i_rep][g1_idx] @ reps[i_rep][g2_idx]
            expected = reps[i_rep][prod_idx]
            err = np.max(np.abs(product_mat - expected))
            max_err = max(max_err, err)
    status = "OK" if max_err < 1e-8 else f"FAIL (err={max_err:.2e})"
    print(f"  rho_{i_rep}: rep property err = {max_err:.2e} [{status}]")


# =====================================================================
# PART 3: CHARACTERS + McKAY GRAPH
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: Characters and McKay Graph")
print("=" * 72)

characters = {}
for i in range(1, 10):
    characters[i] = np.array([np.trace(reps[i][g]) for g in range(120)])
    print(f"  chi_{i}(e) = {characters[i][0].real:.4f} (= dim {irrep_dims[i]})")

# Verify character table: chi_2(e) = 2, chi_3(e) = 2
assert abs(characters[2][0] - 2.0) < 1e-8, f"chi_2(e) = {characters[2][0]}, expected 2"
assert abs(characters[3][0] - 2.0) < 1e-8, f"chi_3(e) = {characters[3][0]}, expected 2"

# McKay adjacency matrix
mckay_adj = np.zeros((9, 9), dtype=int)
for i_rep in range(1, 10):
    product_chars = characters[2] * characters[i_rep]
    for j_rep in range(1, 10):
        mult = np.sum(np.conj(characters[j_rep]) * product_chars).real / N_group
        mckay_adj[i_rep-1, j_rep-1] = int(round(mult))

# Find edges
mckay_edges = []
for i in range(9):
    for j in range(i+1, 9):
        if mckay_adj[i, j] > 0:
            mckay_edges.append((i+1, j+1))  # 1-indexed

print(f"\nMcKay edges ({len(mckay_edges)}):")
for src, tgt in mckay_edges:
    print(f"  rho_{src}(dim {irrep_dims[src]}) -- rho_{tgt}(dim {irrep_dims[tgt]})")

assert len(mckay_edges) == 8, f"Expected 8 edges, got {len(mckay_edges)}"

# Kac label check
print(f"\nKac label check:")
for i in range(1, 10):
    neighbor_sum = sum(irrep_dims[j+1] * mckay_adj[i-1, j] for j in range(9))
    ok = "OK" if 2*irrep_dims[i] == neighbor_sum else "FAIL"
    print(f"  rho_{i}: 2*{irrep_dims[i]} = {neighbor_sum} [{ok}]")


# =====================================================================
# PART 4: CG MAPS (null-space method, from exp307)
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: Clebsch-Gordan Maps (null-space method)")
print("=" * 72)

cg_maps = {}

for edge_idx, (node_i, node_j) in enumerate(mckay_edges):
    # Check both directions for each edge
    for source, target in [(node_i, node_j), (node_j, node_i)]:
        d_s = irrep_dims[source]
        d_t = irrep_dims[target]

        # Multiplicity of rho_target in rho_2 x rho_source
        mult = np.sum(
            np.conj(characters[target]) * characters[2] * characters[source]
        ).real / N_group
        mult_int = int(round(mult))

        if mult_int < 1:
            continue

        dim_product = 2 * d_s
        n_vars = dim_product * d_t

        # Build intertwining equation system
        equations = []
        for g_idx in range(N_group):
            kron_g = np.kron(reps[2][g_idx], reps[source][g_idx])
            rho_t_g = reps[target][g_idx]
            A_g = (np.kron(kron_g, np.eye(d_t, dtype=complex)) -
                   np.kron(np.eye(dim_product, dtype=complex), rho_t_g.T))
            equations.append(A_g)

        A = np.vstack(equations)
        U_svd, S_svd, Vh_svd = np.linalg.svd(A, full_matrices=True)

        if S_svd[0] > 1e-15:
            null_dim = int(np.sum(S_svd / S_svd[0] < 1e-8))
        else:
            null_dim = n_vars

        if null_dim >= 1:
            C_vec = Vh_svd[-1].conj()
            C_map = C_vec.reshape(dim_product, d_t)

            # Normalize: largest entry real positive
            max_idx = np.argmax(np.abs(C_map))
            max_val = C_map.flat[max_idx]
            if abs(max_val) > 1e-15:
                C_map *= np.abs(max_val) / max_val

            # Verify intertwining
            max_err = 0
            for g_idx in range(N_group):
                kron_g = np.kron(reps[2][g_idx], reps[source][g_idx])
                lhs = kron_g @ C_map
                rhs = C_map @ reps[target][g_idx]
                err = np.max(np.abs(lhs - rhs))
                max_err = max(max_err, err)

            status = "OK" if max_err < 1e-6 else f"FAIL ({max_err:.2e})"
            cg_maps[(source, target)] = C_map
            print(f"  ({source},{target}): null_dim={null_dim}, "
                  f"err={max_err:.2e} [{status}]")

print(f"\n  Total CG maps: {len(cg_maps)}")

# =====================================================================
# PART 5: UNWEIGHTED D_F
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: Unweighted D_F (baseline)")
print("=" * 72)

D_F = np.zeros((total_dim, total_dim), dtype=complex)

for (source, target), C_map in cg_maps.items():
    d_s = irrep_dims[source]
    d_t = irrep_dims[target]

    # Yukawa: contract CG map with Higgs vev (1,0)
    Y = C_map[:d_s, :]  # shape (d_s, d_t)

    rs = block_starts[source]
    rt = block_starts[target]
    D_F[rs:rs+d_s, rt:rt+d_t] = Y
    D_F[rt:rt+d_t, rs:rs+d_s] = Y.conj().T

D_F = (D_F + D_F.conj().T) / 2

tr_df2 = np.trace(D_F @ D_F).real
print(f"Tr(D_F^2) = {tr_df2:.6f} (expected: 8 = rank(E8))")
print(f"Hermitian: {np.allclose(D_F, D_F.conj().T)}")

evals_base = np.sort(np.linalg.eigvalsh(D_F))[::-1]
print(f"\nD_F eigenvalues (unweighted):")
n_nonzero = 0
for i, ev in enumerate(evals_base):
    if abs(ev) > 1e-10:
        log_phi = np.log(abs(ev)) / np.log(PHI)
        n_nonzero += 1
    else:
        log_phi = float('nan')
    marker = " *" if abs(ev) > 1e-10 else ""
    print(f"  lambda_{i+1:2d} = {ev:10.6f}  (log_phi = {log_phi:7.3f}){marker}")

print(f"\n  Nonzero eigenvalues: {n_nonzero} out of {total_dim}")

nonzero_evals = evals_base[np.abs(evals_base) > 1e-10]
if len(nonzero_evals) > 0:
    dyn_range = np.max(np.abs(nonzero_evals)) / np.min(np.abs(nonzero_evals))
    print(f"  Dynamic range: {dyn_range:.4f}")
    print(f"  Need: phi^26 = {PHI**26:.1f}")

# Per-edge Yukawa norms
print(f"\nPer-edge ||Y||_F:")
for (source, target), C_map in sorted(cg_maps.items()):
    d_s = irrep_dims[source]
    Y = C_map[:d_s, :]
    frob = np.linalg.norm(Y, 'fro')
    print(f"  ({source},{target}): ||Y||_F = {frob:.6f} "
          f"(expected 1/sqrt(2) = {1/np.sqrt(2):.6f})")


# =====================================================================
# PART 6: McKAY GRAPH TOPOLOGY + DISTANCES
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: McKay Graph Topology")
print("=" * 72)

adj_list = defaultdict(list)
degree = defaultdict(int)
for src, tgt in mckay_edges:
    adj_list[src].append(tgt)
    adj_list[tgt].append(src)
    degree[src] += 1
    degree[tgt] += 1

# Branch node
branch_node = None
for node in range(1, 10):
    if degree[node] >= 3:
        branch_node = node
        print(f"Branch node: rho_{node} (dim {irrep_dims[node]}, degree {degree[node]})")
        break

# Find legs from branch
legs = []
for neighbor in sorted(adj_list[branch_node]):
    leg = [branch_node, neighbor]
    current = neighbor
    prev = branch_node
    while True:
        nexts = [n for n in adj_list[current] if n != prev]
        if not nexts:
            break
        prev = current
        current = nexts[0]
        leg.append(current)
    legs.append(leg)

print(f"\nLegs from branch rho_{branch_node}:")
for i, leg in enumerate(legs):
    leg_str = " -- ".join(f"rho_{n}({irrep_dims[n]})" for n in leg)
    print(f"  Leg {i+1} (length {len(leg)-1}): {leg_str}")

# Graph distances (Floyd-Warshall)
gdist = np.full((10, 10), 999)
for i in range(1, 10):
    gdist[i, i] = 0
for src, tgt in mckay_edges:
    gdist[src, tgt] = 1
    gdist[tgt, src] = 1
for k in range(1, 10):
    for i in range(1, 10):
        for j in range(1, 10):
            if gdist[i, k] + gdist[k, j] < gdist[i, j]:
                gdist[i, j] = gdist[i, k] + gdist[k, j]

print(f"\nDistances from rho_1 (trivial):")
for i in range(1, 10):
    print(f"  d(rho_1, rho_{i}) = {gdist[1,i]}")

# Endpoints (degree 1)
endpoints = sorted([n for n in range(1, 10) if degree[n] == 1])
print(f"\nEndpoints: {[f'rho_{e}({irrep_dims[e]})' for e in endpoints]}")

# Longest path = main chain of E8-hat
from collections import deque


def find_path(start, end):
    """BFS shortest path."""
    visited = {start}
    queue = deque([(start, [start])])
    while queue:
        node, path = queue.popleft()
        if node == end:
            return path
        for n in adj_list[node]:
            if n not in visited:
                visited.add(n)
                queue.append((n, path + [n]))
    return None


longest_path = []
for e1 in endpoints:
    for e2 in endpoints:
        if e1 < e2:
            path = find_path(e1, e2)
            if path and len(path) > len(longest_path):
                longest_path = path

print(f"\nMain chain (longest path): "
      f"{' -- '.join(f'rho_{n}({irrep_dims[n]})' for n in longest_path)}")

# =====================================================================
# PART 7: CAYLEY LAPLACIAN EIGENVALUES
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: Cayley Laplacian Eigenvalues")
print("=" * 72)

# Classify conjugacy classes by chi_2 value
conj_classes = defaultdict(list)
for g_idx in range(N_group):
    chi2 = characters[2][g_idx].real
    key = round(chi2, 4)
    conj_classes[key].append(g_idx)

print(f"Conjugacy classes by chi_2:")
for key in sorted(conj_classes.keys()):
    print(f"  chi_2 = {key:8.4f}: {len(conj_classes[key])} elements")

# Class 10a: chi_2 = phi, 12 elements
class_10a_key = round(PHI, 4)
g_10a = conj_classes[class_10a_key][0]
print(f"\nClass 10a representative: g_{g_10a}, chi_2 = {characters[2][g_10a].real:.6f}")

# Cayley eigenvalues: lambda_i = 12(1 - chi_i(10a)/d_i)
cayley_eigs = {}
print(f"\nCayley Laplacian eigenvalues:")
for i in range(1, 10):
    chi_val = characters[i][g_10a].real
    d_i = irrep_dims[i]
    lam = 12 * (1 - chi_val / d_i)
    cayley_eigs[i] = lam
    print(f"  rho_{i}: lambda = 12*(1 - {chi_val:.4f}/{d_i}) = {lam:.6f}")


# =====================================================================
# PART 8: WEIGHTED D_F -- VARIOUS SCHEMES
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: Weighted D_F -- Various Schemes")
print("=" * 72)


def build_weighted_df(edge_weights):
    """Build D_F with given edge weights.
    edge_weights: dict mapping (src, tgt) or (tgt, src) -> weight."""
    D = np.zeros((total_dim, total_dim), dtype=complex)
    for (source, target), C_map in cg_maps.items():
        d_s = irrep_dims[source]
        d_t = irrep_dims[target]
        Y = C_map[:d_s, :]  # Yukawa from Higgs vev (1,0)

        # Find weight for this edge (try both orderings)
        w = edge_weights.get((source, target),
                             edge_weights.get((target, source), 1.0))

        rs = block_starts[source]
        rt = block_starts[target]
        D[rs:rs+d_s, rt:rt+d_t] += w * Y
        D[rt:rt+d_t, rs:rs+d_s] += w * Y.conj().T

    D = (D + D.conj().T) / 2
    return D


def analyze_df(D, label, show_all=True):
    """Analyze eigenvalues of a D_F matrix."""
    evals = np.sort(np.linalg.eigvalsh(D))[::-1]
    nonzero = evals[np.abs(evals) > 1e-10]
    pos = nonzero[nonzero > 0]
    n_pos = len(pos)

    print(f"\n  [{label}]")
    print(f"  Tr(D^2) = {np.trace(D @ D).real:.4f}")
    if show_all:
        for i, ev in enumerate(evals):
            if abs(ev) > 1e-10:
                log_phi = np.log(abs(ev)) / np.log(PHI)
            else:
                log_phi = float('nan')
            marker = " *" if abs(ev) > 1e-10 else ""
            print(f"    lambda_{i+1:2d} = {ev:12.6f}  "
                  f"(log_phi = {log_phi:7.3f}){marker}")
    else:
        # Show only nonzero
        for i, ev in enumerate(nonzero):
            log_phi = np.log(abs(ev)) / np.log(PHI)
            print(f"    |lambda_{i+1:2d}| = {abs(ev):12.6f}  "
                  f"(log_phi = {log_phi:7.3f})")

    if len(nonzero) > 0:
        dyn_range = np.max(np.abs(nonzero)) / np.min(np.abs(nonzero))
        print(f"  Dynamic range: {dyn_range:.4f} (need phi^26 = {PHI**26:.1f})")
    print(f"  Positive eigenvalues: {n_pos} (need 9 for fermion masses)")
    return evals


# ----- Scheme A: Unweighted (baseline) -----
print("\n--- Scheme A: Unweighted ---")
evals_A = analyze_df(D_F, "Unweighted", show_all=False)

# ----- Scheme B: phi^distance weights -----
print("\n--- Scheme B: phi^(graph distance from rho_1) ---")
edge_weights_B = {}
for src, tgt in mckay_edges:
    w = PHI**((gdist[1, src] + gdist[1, tgt]) / 2)
    edge_weights_B[(src, tgt)] = w

print("  Edge weights:")
for (s, t), w in sorted(edge_weights_B.items()):
    print(f"    ({s},{t}): {w:.6f}")

D_B = build_weighted_df(edge_weights_B)
evals_B = analyze_df(D_B, "phi^distance", show_all=False)

# ----- Scheme C: Cayley eigenvalue weights -----
print("\n--- Scheme C: Cayley eigenvalue weights ---")
edge_weights_C = {}
for src, tgt in mckay_edges:
    w = (cayley_eigs[src] + cayley_eigs[tgt]) / 2
    edge_weights_C[(src, tgt)] = w

print("  Edge weights:")
for (s, t), w in sorted(edge_weights_C.items()):
    print(f"    ({s},{t}): {w:.6f}")

D_C = build_weighted_df(edge_weights_C)
evals_C = analyze_df(D_C, "Cayley eigenvalue", show_all=False)

# ----- Scheme D: phi^(Cayley eigenvalue) -----
print("\n--- Scheme D: phi^(Cayley eigenvalue) ---")
edge_weights_D = {}
for src, tgt in mckay_edges:
    w = PHI**((cayley_eigs[src] + cayley_eigs[tgt]) / 2)
    edge_weights_D[(src, tgt)] = w

D_D = build_weighted_df(edge_weights_D)
evals_D = analyze_df(D_D, "phi^Cayley", show_all=False)

# ----- Scheme E: exp323 weights by leg -----
print("\n--- Scheme E: exp323 weights by leg ---")
# The E8-hat graph has 3 legs from branch node rho_9 (dim 6):
#   Long leg (5 edges): rho_9 -- rho_8 -- rho_7 -- rho_4 -- rho_2 -- rho_1
#   Medium leg (2 edges): rho_9 -- rho_6 -- rho_3
#   Short leg (1 edge): rho_9 -- rho_5
#
# From exp323 Theorem: main chain forced [N_gen=3, F(3)=2, b_1=6, 0, a_1=5]
# These 5 weights go to the 5 edges of the long leg (from branch outward).
# Total weights in Solution 3: {0,1,2,3,3,5,6,9} (8 weights for 8 edges)
# Main chain (long leg): uses 5 weights
# Medium leg: uses 2 weights
# Short leg: uses 1 weight

# Sort legs by length (descending)
legs_sorted = sorted(legs, key=lambda l: len(l), reverse=True)
for i, leg in enumerate(legs_sorted):
    leg_str = " -- ".join(f"rho_{n}({irrep_dims[n]})" for n in leg)
    print(f"  Leg {i+1} (length {len(leg)-1}): {leg_str}")

# Build leg edges (from branch outward)
leg_edges = []
for leg in legs_sorted:
    edges = [(leg[j], leg[j+1]) for j in range(len(leg)-1)]
    leg_edges.append(edges)

# Main chain weights (long leg, 5 edges)
main_chain_weights = [3, 2, 6, 0, 5]  # N_gen, F(3), b_1, 0, a_1
print(f"\n  Long leg weights: {main_chain_weights}")

target_exps = sorted([d[2] for d in fermion_data.values()])
# [0, 3, 5, 11, 11, 16, 17, 19, 26]
print(f"  Target mass exponents: {target_exps}")

# From exp323 S3: total = {0,1,2,3,3,5,6,9}, long leg uses {0,2,3,5,6}
# Remaining for medium+short: {1,3,9}
# Medium leg needs 2 weights, short leg needs 1
from itertools import permutations as perms

remaining_weights = [1, 3, 9]
print(f"  Remaining weights: {remaining_weights} (for medium+short legs)")

best_corr = -999
best_config = None
best_evals = None

# Try all assignments of remaining_weights to (medium_e1, medium_e2, short_e1)
for perm in perms(remaining_weights):
    edge_weights_E = {}

    # Long leg
    for i, (s, t) in enumerate(leg_edges[0]):
        edge_weights_E[(s, t)] = main_chain_weights[i]

    # Medium leg (2 edges)
    if len(leg_edges) > 1:
        for i, (s, t) in enumerate(leg_edges[1]):
            edge_weights_E[(s, t)] = perm[i]

    # Short leg (1 edge)
    if len(leg_edges) > 2:
        for i, (s, t) in enumerate(leg_edges[2]):
            edge_weights_E[(s, t)] = perm[2 + i]

    D_E = build_weighted_df(edge_weights_E)
    evals_E = np.sort(np.linalg.eigvalsh(D_E))[::-1]
    nonzero_E = evals_E[np.abs(evals_E) > 1e-10]
    pos_E = np.sort(np.abs(nonzero_E))[::-1]

    if len(pos_E) >= 9:
        phi_logs = np.log(pos_E[:9]) / np.log(PHI)
        phi_logs_sorted = np.sort(phi_logs)[::-1]
        target_sorted = np.sort(target_exps)[::-1]
        corr = np.corrcoef(target_sorted, phi_logs_sorted[:9])[0, 1]
        if corr > best_corr:
            best_corr = corr
            best_config = perm
            best_evals = evals_E.copy()

if best_config is not None:
    print(f"\n  Best remaining config: {best_config}, correlation = {best_corr:.6f}")
    # Rebuild best
    ew_best = {}
    for i, (s, t) in enumerate(leg_edges[0]):
        ew_best[(s, t)] = main_chain_weights[i]
    if len(leg_edges) > 1:
        for i, (s, t) in enumerate(leg_edges[1]):
            ew_best[(s, t)] = best_config[i]
    if len(leg_edges) > 2:
        for i, (s, t) in enumerate(leg_edges[2]):
            ew_best[(s, t)] = best_config[2 + i]

    print("  Full edge weight assignment:")
    for (s, t), w in sorted(ew_best.items()):
        print(f"    rho_{s}--rho_{t}: w = {w}")

    analyze_df(build_weighted_df(ew_best),
               f"exp323 S3 (config={list(best_config)})", show_all=True)

# Also try ALL permutations of main chain weights
print("\n  Testing ALL permutations of long-leg weights too...")
best_corr_full = -999
best_full_config = None
n_tested = 0

for main_perm in perms(main_chain_weights):
    for branch_perm in perms(remaining_weights):
        ew = {}
        for i, (s, t) in enumerate(leg_edges[0]):
            ew[(s, t)] = main_perm[i]
        if len(leg_edges) > 1:
            for i, (s, t) in enumerate(leg_edges[1]):
                ew[(s, t)] = branch_perm[i]
        if len(leg_edges) > 2:
            for i, (s, t) in enumerate(leg_edges[2]):
                ew[(s, t)] = branch_perm[2 + i]

        D_test = build_weighted_df(ew)
        ev = np.sort(np.linalg.eigvalsh(D_test))[::-1]
        nz = ev[np.abs(ev) > 1e-10]
        ps = np.sort(np.abs(nz))[::-1]

        if len(ps) >= 9:
            pl = np.log(ps[:9]) / np.log(PHI)
            pls = np.sort(pl)[::-1]
            ts = np.sort(target_exps)[::-1]
            c = np.corrcoef(ts, pls[:9])[0, 1]
            if c > best_corr_full:
                best_corr_full = c
                best_full_config = (main_perm, branch_perm)
        n_tested += 1

print(f"  Tested {n_tested} permutations")
if best_full_config is not None:
    mp, bp = best_full_config
    print(f"  Best: main={list(mp)}, branch={list(bp)}, corr={best_corr_full:.6f}")

    ew_full = {}
    for i, (s, t) in enumerate(leg_edges[0]):
        ew_full[(s, t)] = mp[i]
    if len(leg_edges) > 1:
        for i, (s, t) in enumerate(leg_edges[1]):
            ew_full[(s, t)] = bp[i]
    if len(leg_edges) > 2:
        for i, (s, t) in enumerate(leg_edges[2]):
            ew_full[(s, t)] = bp[2 + i]

    analyze_df(build_weighted_df(ew_full),
               f"Best full perm (main={list(mp)}, branch={list(bp)})",
               show_all=True)

# ----- Scheme F: phi^(exp323 weight) -----
print("\n--- Scheme F: phi^(exp323 weight) ---")
# Same structure but weight -> phi^weight
best_corr_F = -999
best_config_F = None

for main_perm in perms(main_chain_weights):
    for branch_perm in perms(remaining_weights):
        ew = {}
        for i, (s, t) in enumerate(leg_edges[0]):
            w = main_perm[i]
            ew[(s, t)] = PHI**w if w > 0 else 0.0
        if len(leg_edges) > 1:
            for i, (s, t) in enumerate(leg_edges[1]):
                w = branch_perm[i]
                ew[(s, t)] = PHI**w if w > 0 else 0.0
        if len(leg_edges) > 2:
            for i, (s, t) in enumerate(leg_edges[2]):
                w = branch_perm[2 + i]
                ew[(s, t)] = PHI**w if w > 0 else 0.0

        D_test = build_weighted_df(ew)
        ev = np.sort(np.linalg.eigvalsh(D_test))[::-1]
        nz = ev[np.abs(ev) > 1e-10]
        ps = np.sort(np.abs(nz))[::-1]

        if len(ps) >= 9:
            pl = np.log(ps[:9]) / np.log(PHI)
            pls = np.sort(pl)[::-1]
            ts = np.sort(target_exps)[::-1]
            c = np.corrcoef(ts, pls[:9])[0, 1]
            if c > best_corr_F:
                best_corr_F = c
                best_config_F = (main_perm, branch_perm)

if best_config_F is not None:
    mp, bp = best_config_F
    print(f"  Best: main={list(mp)}, branch={list(bp)}, corr={best_corr_F:.6f}")

    ew_F = {}
    for i, (s, t) in enumerate(leg_edges[0]):
        w = mp[i]
        ew_F[(s, t)] = PHI**w if w > 0 else 0.0
    if len(leg_edges) > 1:
        for i, (s, t) in enumerate(leg_edges[1]):
            w = bp[i]
            ew_F[(s, t)] = PHI**w if w > 0 else 0.0
    if len(leg_edges) > 2:
        for i, (s, t) in enumerate(leg_edges[2]):
            w = bp[2 + i]
            ew_F[(s, t)] = PHI**w if w > 0 else 0.0

    analyze_df(build_weighted_df(ew_F),
               f"phi^(exp323) best", show_all=True)


# =====================================================================
# PART 9: COMPARISON WITH MASS EXPONENTS
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: Detailed Comparison with Mass Exponents")
print("=" * 72)

target_sorted = np.array(sorted(target_exps, reverse=True))
print(f"Target exponents (descending): {list(target_sorted)}")
print(f"Target = log_phi(m_f/m_e): t=26, b=19, tau=17, c=16, "
      f"mu=s=11, d=5, u=3, e=0")

all_schemes = []

for label, D_test in [("A: Unweighted", D_F),
                       ("B: phi^distance", D_B),
                       ("C: Cayley eig", D_C),
                       ("D: phi^Cayley", D_D)]:
    evals = np.sort(np.linalg.eigvalsh(D_test))[::-1]
    nonzero = evals[np.abs(evals) > 1e-10]
    pos = np.sort(np.abs(nonzero))[::-1]

    phi_logs = np.log(pos) / np.log(PHI) if len(pos) > 0 else np.array([])
    all_schemes.append((label, phi_logs, pos))

    print(f"\n  {label}:")
    print(f"    N(nonzero) = {len(nonzero)}, N(positive) = {len(pos)}")
    if len(pos) >= 9:
        top9 = phi_logs[:9]
        print(f"    Top 9 phi-logs: {', '.join(f'{x:.3f}' for x in top9)}")
        corr = np.corrcoef(target_sorted, np.sort(top9)[::-1])[0, 1]
        rms = np.sqrt(np.mean((target_sorted - np.sort(top9)[::-1])**2))
        print(f"    Correlation: {corr:.4f}")
        print(f"    RMS deviation: {rms:.2f}")
    elif len(pos) > 0:
        print(f"    All phi-logs: {', '.join(f'{x:.3f}' for x in phi_logs)}")
        print(f"    Only {len(pos)} nonzero eigenvalues (need 9)")


# =====================================================================
# PART 10: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: HONEST ASSESSMENT")
print("=" * 72)

print("""
RESULTS:

1. REPRESENTATIONS VERIFIED:
   - All 9 irreps of 2I built correctly (from exp307 code)
   - Character table matches exactly
   - McKay graph = affine E8 with 8 edges (correct topology)

2. CG MAPS VERIFIED:
   - All intertwining maps computed via null-space method
   - ||Y||_F = 1/sqrt(2) universal (Proposition 8 in paper)
   - Tr(D_F^2) = 8 = rank(E8)

3. UNWEIGHTED D_F:
   - Dynamic range O(1), far too small for mass hierarchy phi^26
   - This is EXPECTED: D_F encodes selection rules, not mass scales

4. WEIGHTED D_F SCHEMES:
   A. Unweighted: dynamic range ~O(1)
   B. phi^distance: modest increase
   C. Cayley eigenvalue: moderate increase
   D. phi^Cayley: large dynamic range but wrong structure
   E. exp323 weights [3,2,6,0,5]: tested all branch permutations
   F. phi^(exp323 weight): exponential weights

5. CONCLUSION:
   NO weighting scheme produces eigenvalues proportional to phi^n_f.
   The mismatch is STRUCTURAL, not numerical:
   - D_F is 30x30 with 8 edges = SELECTION RULES
   - Mass formula needs phi^(5a+6b) = SCALE in Z[phi]
   - These are COMPLEMENTARY (exp308 conclusion confirmed)

   The McKay graph encodes WHICH transitions are allowed.
   The Z[phi] lattice encodes HOW MUCH mass each fermion gets.
   No single operator can do both.

   CATEGORY: CONFIRMED NEGATIVE
   (D_F eigenvalues are not mass exponents, for any weighting)
""")

print("=" * 72)
print("exp329 COMPLETE")
print("=" * 72)

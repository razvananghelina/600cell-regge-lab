"""
EXP-315: McKay D_F Sub-Graph Eigenvalues vs Gauge Corrections
=============================================================

From exp313-314: Berry phase on 600-cell gives UNIVERSAL baseline 1/phi^4
for ALL gauge sectors. Sector-dependent corrections must come from elsewhere.

HYPOTHESIS: The finite Dirac operator D_F on the McKay graph (affine E8)
produces sector-dependent corrections through its sub-graph spectra.

The 3 legs of affine E8 (from branch at rho_9, dim 6):
  - Short (1 edge): rho_5(3')  -- linked to SU(2) adjoint
  - Medium (2 edges): rho_3(2)--rho_6(4)
  - Long (5 edges): rho_1(1)--rho_2(2)--rho_4(3)--rho_7(4)--rho_8(5)

Target corrections (from holonomy-corrected masses):
  Lepton:    delta_k = -1/phi^4 = -0.14590,  delta_d = sin^2(tW) = 6/26
  Up quark:  delta_k = +alpha_s = +1/(2*phi^3),  delta_d = alpha_s*2/pi
  Down quark: delta_k = -alpha_s,  delta_d = -1/a1

Tests:
1. Build D_F (9x9 adjacency and 30x30 CG-based) with multiple weight schemes
2. Restrict to sub-graphs for each leg (with and without branch node)
3. Compute sub-graph eigenvalues and spectral invariants
4. Compare with alpha_s, sin^2(tW), 1/a1, and 1/phi^4

Author: Computational exploration for R.-C. Anghelina's framework
"""

import numpy as np
from scipy.special import comb as binomial
from itertools import product as iprod

# ============================================================
# Constants from a1 = 5
# ============================================================
a1 = 5
b1 = a1 + 1
PHI = (1 + np.sqrt(5)) / 2
PHI_PRIME = 1 - PHI
N = 120  # = a1!
h_E8 = a1 * b1  # = 30

# Framework coupling constants
alpha_s_fw = 1.0 / (2.0 * PHI**3)
sin2tW_fw = float(b1) / (a1**2 + 1)  # = 6/26
inv_a1 = 1.0 / a1
inv_phi4 = 1.0 / PHI**4

print("=" * 70)
print("EXP-315: McKAY D_F SUB-GRAPH EIGENVALUES vs GAUGE CORRECTIONS")
print("=" * 70)
print(f"a1 = {a1}, phi = {PHI:.10f}, phi' = {PHI_PRIME:.10f}")
print(f"\nTarget corrections:")
print(f"  alpha_s     = 1/(2*phi^3) = {alpha_s_fw:.6f}")
print(f"  sin^2(tW)   = {b1}/{a1**2+1} = {sin2tW_fw:.6f}")
print(f"  1/a1        = 1/{a1} = {inv_a1:.6f}")
print(f"  1/phi^4     = {inv_phi4:.6f}  (Berry phase baseline)")


# ============================================================
# SECTION 1: Build 2I group and representations (compact)
# ============================================================
print("\n" + "=" * 70)
print("SECTION 1: Building 2I group and McKay graph")
print("=" * 70)


def qmul(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
            w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)


def build_2I():
    elements = []
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0]*4; v[i] = s; elements.append(tuple(v))
    for signs in iprod([0.5, -0.5], repeat=4):
        elements.append(tuple(signs))
    even_perms = [(0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
                  (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0)]
    bv = [0.0, 0.5, PHI/2, 1/(2*PHI)]
    for perm in even_perms:
        base = [bv[perm[k]] for k in range(4)]
        nonzero = [k for k in range(4) if abs(base[k]) > 1e-10]
        for signs in iprod([1, -1], repeat=len(nonzero)):
            v = list(base)
            for idx, k in enumerate(nonzero):
                v[k] *= signs[idx]
            elements.append(tuple(v))
    unique = []
    for v in elements:
        if not any(sum((v[k]-u[k])**2 for k in range(4)) < 1e-10 for u in unique):
            unique.append(v)
    return unique


def quat_to_su2(q):
    w, x, y, z = q
    return np.array([[complex(w,x), complex(y,z)], [-complex(y,-z), complex(w,-x)]])


def galois_complex(z):
    def gal_real(x):
        for b_num in range(-8, 9):
            b = b_num / 2.0
            a = x - b * PHI
            a2 = round(2 * a)
            if abs(a - a2/2) < 1e-8:
                return (a2/2 + b) - b * PHI
        return x
    return gal_real(z.real) + 1j * gal_real(z.imag)


def galois_matrix(M):
    result = np.empty_like(M)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            result[i,j] = galois_complex(M[i,j])
    return result


def spin_j_matrix(U, j):
    two_j = int(round(2*j)); dim = two_j + 1
    a, b = U[0,0], U[0,1]; c, d = U[1,0], U[1,1]
    D = np.zeros((dim, dim), dtype=complex)
    for m_idx in range(dim):
        m = j - m_idx; p = int(round(j+m)); q = int(round(j-m))
        for mp_idx in range(dim):
            mp = j - mp_idx; pp = int(round(j+mp))
            val = 0j
            for s in range(p+1):
                t = pp - s
                if 0 <= t <= q:
                    val += (binomial(p,s,exact=True) * binomial(q,t,exact=True) *
                            a**s * c**(p-s) * b**t * d**(q-t))
            D[mp_idx, m_idx] = val
    return D


# Build group
elements_2I = build_2I()
assert len(elements_2I) == 120, f"Expected 120, got {len(elements_2I)}"
print(f"  |2I| = {len(elements_2I)}")

# Build representations
irrep_dims = {1:1, 2:2, 3:2, 4:3, 5:3, 6:4, 7:4, 8:5, 9:6}
dims = [irrep_dims[i] for i in range(1, 10)]
reps = {i: [] for i in range(1, 10)}

for q in elements_2I:
    U = quat_to_su2(q)
    U_gal = galois_matrix(U)
    reps[1].append(np.array([[1.0+0j]]))
    reps[2].append(U.copy())
    reps[3].append(U_gal.copy())
    reps[4].append(spin_j_matrix(U, 1))
    reps[5].append(spin_j_matrix(U_gal, 1))
    reps[6].append(np.kron(U, U_gal))
    reps[7].append(spin_j_matrix(U, 1.5))
    reps[8].append(spin_j_matrix(U, 2))
    reps[9].append(spin_j_matrix(U, 2.5))

# Characters
characters = {}
for i in range(1, 10):
    characters[i] = np.array([np.trace(reps[i][g]) for g in range(120)])

# McKay adjacency: rho_j in rho_2 x rho_i
mckay_adj = np.zeros((9, 9), dtype=int)
for i in range(1, 10):
    prod_chars = characters[2] * characters[i]
    for j in range(1, 10):
        mult = int(round(np.sum(np.conj(characters[j]) * prod_chars).real / 120))
        if mult > 0:
            mckay_adj[i-1, j-1] = mult

# Print adjacency
print("\n  McKay adjacency matrix (affine E8):")
for i in range(9):
    row = " ".join(str(mckay_adj[i,j]) for j in range(9))
    print(f"    rho_{i+1}(d={dims[i]}): [{row}]")

# Edge list
mckay_edges = [(i, j) for i in range(9) for j in range(i+1, 9) if mckay_adj[i,j] > 0]
print(f"\n  Edges ({len(mckay_edges)}):")
for i, j in mckay_edges:
    print(f"    rho_{i+1}(d={dims[i]}) -- rho_{j+1}(d={dims[j]})")

# Find branch and legs
from collections import deque
degrees = [sum(1 for j in range(9) if mckay_adj[i,j] > 0) for i in range(9)]
branch_node = next(i for i in range(9) if degrees[i] == 3)
print(f"\n  Branch node: rho_{branch_node+1} (dim {dims[branch_node]}, degree {degrees[branch_node]})")

# BFS from branch
def bfs_dist(adj, src):
    d = [-1]*9; d[src] = 0; q = deque([src])
    while q:
        u = q.popleft()
        for v in range(9):
            if adj[u,v] > 0 and d[v] == -1:
                d[v] = d[u]+1; q.append(v)
    return d

dist_from_branch = bfs_dist(mckay_adj, branch_node)

# Find 3 legs (connected components when branch removed)
def find_legs(adj, branch):
    legs = []
    neighbors = [j for j in range(9) if adj[branch,j] > 0]
    for start in neighbors:
        leg = [start]; current, prev = start, branch
        while True:
            nexts = [j for j in range(9) if adj[current,j] > 0 and j != prev]
            if not nexts:
                break
            prev, current = current, nexts[0]
            leg.append(current)
        legs.append(leg)
    return legs

legs = find_legs(mckay_adj, branch_node)
legs.sort(key=len)

print(f"\n  Three legs from branch (sorted by length):")
leg_names = ['Short', 'Medium', 'Long']
for li, (leg, name) in enumerate(zip(legs, leg_names)):
    leg_str = " -- ".join(f"rho_{n+1}(d={dims[n]})" for n in leg)
    total_dim = sum(dims[n] for n in leg)
    print(f"    {name} ({len(leg)} nodes, dim {total_dim}): {leg_str}")

# Verify dimensions
print(f"\n  Dimension check:")
print(f"    sum(dims) = {sum(dims)} (expected h(E8) = {h_E8})")
print(f"    sum(dims^2) = {sum(d**2 for d in dims)} (expected |2I| = {N})")
print(f"    Perron-Frobenius check: A*d = 2*d?")
d_vec = np.array(dims, dtype=float)
Ad = mckay_adj @ d_vec
pf_check = np.allclose(Ad, 2 * d_vec)
print(f"      {pf_check} (Ad = {Ad.astype(int).tolist()}, 2d = {(2*d_vec).astype(int).tolist()})")


# ============================================================
# SECTION 2: Build D_F matrices
# ============================================================
print("\n" + "=" * 70)
print("SECTION 2: Building D_F (9x9 and 30x30)")
print("=" * 70)

# --- 2A: Simple 9x9 adjacency-based D_F ---
print("\n--- 2A: 9x9 adjacency matrix (unit weights) ---")
A9 = mckay_adj.astype(float)
eigs_A9 = np.sort(np.linalg.eigvalsh(A9))[::-1]
print(f"  Eigenvalues of A(9x9): {np.round(eigs_A9, 6).tolist()}")
print(f"  Largest eigenvalue: {eigs_A9[0]:.6f} (should be 2.0)")

# 9x9 dimension-weighted: D[i,j] = sqrt(d_i*d_j) if connected
print("\n--- 2B: 9x9 dimension-weighted D_F ---")
D9_dimwt = np.zeros((9, 9))
for i, j in mckay_edges:
    w = np.sqrt(dims[i] * dims[j])
    D9_dimwt[i,j] = w
    D9_dimwt[j,i] = w
eigs_D9 = np.sort(np.linalg.eigvalsh(D9_dimwt))[::-1]
tr2_D9 = np.trace(D9_dimwt @ D9_dimwt)
print(f"  Tr(D^2) = {tr2_D9:.4f}")
print(f"  Eigenvalues: {np.round(eigs_D9, 6).tolist()}")

# --- 2C: 30x30 CG-based D_F (from exp308) ---
print("\n--- 2C: 30x30 CG-based D_F ---")
print("  Computing CG maps...")

# CG maps via null-space method
cg_maps = {}
mckay_edges_1idx = [(i+1, j+1) for i, j in mckay_edges]

for (ni, nj) in mckay_edges_1idx:
    for source, target in [(ni, nj), (nj, ni)]:
        d_s, d_t = irrep_dims[source], irrep_dims[target]
        mult = np.sum(np.conj(characters[target]) * characters[2] * characters[source]).real / 120
        if round(mult) < 1:
            continue
        dim_p = 2 * d_s
        equations = []
        for g in range(120):
            kron_g = np.kron(reps[2][g], reps[source][g])
            rho_t_g = reps[target][g]
            A_g = (np.kron(kron_g, np.eye(d_t, dtype=complex)) -
                   np.kron(np.eye(dim_p, dtype=complex), rho_t_g.T))
            equations.append(A_g)
        A_sys = np.vstack(equations)
        _, S, Vh = np.linalg.svd(A_sys, full_matrices=True)
        null_dim = int(np.sum(S / S[0] < 1e-8)) if S[0] > 1e-15 else dim_p * d_t
        if null_dim >= 1:
            C_vec = Vh[-1].conj()
            C_map = C_vec.reshape(dim_p, d_t)
            max_val = C_map.flat[np.argmax(np.abs(C_map))]
            if abs(max_val) > 1e-15:
                C_map *= np.abs(max_val) / max_val
            cg_maps[(source, target)] = C_map

print(f"  {len(cg_maps)} CG maps computed.")

# Build 30x30 D_F
total_dim = sum(dims)  # = 30
block_starts = {}
offset = 0
for i in range(1, 10):
    block_starts[i] = offset
    offset += irrep_dims[i]

# Per-edge D_F contributions
edge_D = {}
edge_list_1idx = sorted(set(tuple(sorted([s,t])) for (s,t) in cg_maps.keys()))

for (n1, n2) in edge_list_1idx:
    D_e = np.zeros((total_dim, total_dim), dtype=complex)
    for (src, tgt) in [(n1,n2), (n2,n1)]:
        if (src, tgt) in cg_maps:
            C = cg_maps[(src, tgt)]
            d_s = irrep_dims[src]
            Y = C[:d_s, :]
            rs = block_starts[src]
            cs = block_starts[tgt]
            D_e[rs:rs+d_s, cs:cs+irrep_dims[tgt]] = Y
            D_e[cs:cs+irrep_dims[tgt], rs:rs+d_s] = Y.conj().T
    D_e = (D_e + D_e.conj().T) / 2
    edge_D[(n1, n2)] = D_e

# D_F with unit weights
D_F30 = sum(edge_D.values())
D_F30 = (D_F30 + D_F30.conj().T) / 2

eigs_30, evecs_30 = np.linalg.eigh(D_F30)
idx = np.argsort(np.abs(eigs_30))[::-1]
eigs_30 = eigs_30[idx]
evecs_30 = evecs_30[:, idx]

tr2_30 = np.real(np.trace(D_F30 @ D_F30))
print(f"  D_F(30x30): Tr(D^2) = {tr2_30:.4f} (should be 8 = rank(E8))")
print(f"  Eigenvalues (sorted by |lambda|):")
for k in range(len(eigs_30)):
    if abs(eigs_30[k]) > 1e-10:
        print(f"    lambda_{k+1:2d} = {eigs_30[k]:+.6f}")

# Count zero eigenvalues
n_zero = sum(1 for e in eigs_30 if abs(e) < 1e-10)
print(f"  Zero eigenvalues: {n_zero} (of 30)")


# ============================================================
# SECTION 3: Sub-graph decomposition
# ============================================================
print("\n" + "=" * 70)
print("SECTION 3: Sub-Graph Eigenvalue Analysis")
print("=" * 70)

# Physical identification from exp290b:
# 12 = 1 + 3' + (3+5) = 1 + 3 + 8 = U(1) + SU(2) + SU(3)
# rho_1(1) -> U(1), rho_5(3') -> ad(SU(2)), rho_4(3)+rho_8(5) -> ad(SU(3))

print("""
  Gauge sector identification (from exp290b):
    U(1):  rho_1 (dim 1)               -- trivial rep
    SU(2): rho_5 (dim 3')              -- Galois adjoint
    SU(3): rho_4 (dim 3) + rho_8 (dim 5)  -- adjoint decomp

  Graph decomposition (remove branch rho_{0}):
    Short leg:  [{1}]
    Medium leg: [{2}]
    Long leg:   [{3}]
""".format(branch_node+1,
           ", ".join(f"rho_{n+1}(d={dims[n]})" for n in legs[0]),
           ", ".join(f"rho_{n+1}(d={dims[n]})" for n in legs[1]),
           ", ".join(f"rho_{n+1}(d={dims[n]})" for n in legs[2])))

# --- 3A: 9x9 sub-graph analysis ---
print("--- 3A: 9x9 adjacency sub-graph eigenvalues ---")

for variant_name, include_branch in [("without branch", False), ("with branch", True)]:
    print(f"\n  Sub-graphs {variant_name}:")
    for li, (leg, name) in enumerate(zip(legs, leg_names)):
        nodes = list(leg)
        if include_branch:
            nodes = [branch_node] + nodes
        nodes.sort()
        n_nodes = len(nodes)
        total_leg_dim = sum(dims[n] for n in nodes)

        # Extract sub-matrix of A9
        sub_A = A9[np.ix_(nodes, nodes)]
        eigs_sub = np.sort(np.linalg.eigvalsh(sub_A))[::-1]
        n_edges_sub = int(np.trace(sub_A @ sub_A) / 2)

        # Also try dimension-weighted sub-matrix
        sub_Dw = D9_dimwt[np.ix_(nodes, nodes)]
        eigs_sub_dw = np.sort(np.linalg.eigvalsh(sub_Dw))[::-1]
        tr2_sub_dw = np.trace(sub_Dw @ sub_Dw)

        print(f"\n    {name} leg ({n_nodes} nodes, {n_edges_sub} edges, dim_sum={total_leg_dim}):")
        print(f"      Nodes: {['rho_'+str(n+1)+'('+str(dims[n])+')' for n in nodes]}")
        print(f"      Unit-wt eigenvalues: {np.round(eigs_sub, 6).tolist()}")
        if n_nodes > 1:
            spectral_gap = eigs_sub[0] - eigs_sub[1] if n_nodes > 1 else 0
            print(f"      Spectral gap: {spectral_gap:.6f}")

        print(f"      Dim-wt eigenvalues: {np.round(eigs_sub_dw, 6).tolist()}")
        if n_nodes > 1:
            print(f"      Dim-wt Tr(D^2) = {tr2_sub_dw:.6f}")
            print(f"      <lambda^2> = Tr(D^2)/n_nodes = {tr2_sub_dw/n_nodes:.6f}")
            print(f"      <lambda^2>/dim_sum = {tr2_sub_dw/total_leg_dim:.6f}")


# --- 3B: 30x30 sub-graph analysis ---
print("\n\n--- 3B: 30x30 CG-based D_F sub-graph eigenvalues ---")

for variant_name, include_branch in [("without branch", False), ("with branch", True)]:
    print(f"\n  Sub-graphs {variant_name}:")
    for li, (leg, name) in enumerate(zip(legs, leg_names)):
        nodes_0idx = list(leg)
        if include_branch:
            nodes_0idx = [branch_node] + nodes_0idx

        # Build row/column index set for 30x30
        idx_set = []
        for n in sorted(nodes_0idx):
            s = block_starts[n+1]
            e = s + irrep_dims[n+1]
            idx_set.extend(range(s, e))
        idx_set = sorted(idx_set)
        sub_dim = len(idx_set)

        # Extract sub-matrix
        sub_D = D_F30[np.ix_(idx_set, idx_set)]
        sub_D = np.real((sub_D + sub_D.T) / 2)  # ensure hermitian real
        eigs_sub = np.sort(np.linalg.eigvalsh(sub_D))[::-1]
        tr2_sub = np.trace(sub_D @ sub_D)

        # Count edges in this sub-graph
        n_edges_sub = 0
        for (n1, n2) in edge_list_1idx:
            if (n1-1) in nodes_0idx and (n2-1) in nodes_0idx:
                n_edges_sub += 1

        total_leg_dim = sub_dim
        print(f"\n    {name} leg ({len(nodes_0idx)} nodes, {n_edges_sub} edges, "
              f"sub_dim={sub_dim}):")
        print(f"      Nodes: {['rho_'+str(n+1)+'('+str(dims[n])+')' for n in sorted(nodes_0idx)]}")
        print(f"      Tr(D^2) = {tr2_sub:.6f}")
        print(f"      <lambda^2> = Tr(D^2)/sub_dim = {tr2_sub/sub_dim:.6f}")

        nonzero = [e for e in eigs_sub if abs(e) > 1e-10]
        n_zero_sub = sub_dim - len(nonzero)
        print(f"      Nonzero eigenvalues ({len(nonzero)} of {sub_dim}):")
        for e in sorted(nonzero, key=lambda x: -abs(x)):
            print(f"        {e:+.6f}")
        print(f"      Zero eigenvalues: {n_zero_sub}")

        if len(nonzero) > 0:
            max_eig = max(abs(e) for e in nonzero)
            min_eig = min(abs(e) for e in nonzero)
            ratio = max_eig / min_eig if min_eig > 1e-15 else float('inf')
            log_phi_max = np.log(max_eig) / np.log(PHI) if max_eig > 0 else 0
            log_phi_min = np.log(min_eig) / np.log(PHI) if min_eig > 0 else 0
            print(f"      |lambda|_max = {max_eig:.6f} (phi^{log_phi_max:.4f})")
            print(f"      |lambda|_min = {min_eig:.6f} (phi^{log_phi_min:.4f})")
            print(f"      Spectral ratio = {ratio:.6f}")


# ============================================================
# SECTION 4: Comparison with correction parameters
# ============================================================
print("\n" + "=" * 70)
print("SECTION 4: Spectral Quantities vs Correction Parameters")
print("=" * 70)

print(f"""
  Targets:
    alpha_s   = 1/(2*phi^3) = {alpha_s_fw:.6f}
    sin^2(tW) = 6/26       = {sin2tW_fw:.6f}
    1/a1      = 1/5        = {inv_a1:.6f}
    1/phi^4                = {inv_phi4:.6f}
""")

# For each sub-graph, compute spectral quantities and compare
print("  Computing spectral invariants for each leg (with branch):")

leg_data = {}
for li, (leg, name) in enumerate(zip(legs, leg_names)):
    nodes_0idx = [branch_node] + list(leg)

    # 30x30 sub-graph
    idx_set = []
    for n in sorted(nodes_0idx):
        s = block_starts[n+1]
        e = s + irrep_dims[n+1]
        idx_set.extend(range(s, e))
    idx_set = sorted(idx_set)
    sub_dim = len(idx_set)

    sub_D = D_F30[np.ix_(idx_set, idx_set)]
    sub_D = np.real((sub_D + sub_D.T) / 2)
    eigs_sub = np.linalg.eigvalsh(sub_D)

    # Count edges
    n_edges = 0
    for (n1, n2) in edge_list_1idx:
        if (n1-1) in nodes_0idx and (n2-1) in nodes_0idx:
            n_edges += 1

    nonzero = [e for e in eigs_sub if abs(e) > 1e-10]
    tr2 = np.sum(np.array(eigs_sub)**2)

    leg_data[name] = {
        'nodes': nodes_0idx,
        'n_edges': n_edges,
        'sub_dim': sub_dim,
        'tr2': tr2,
        'eigenvalues': sorted(eigs_sub, key=lambda x: -abs(x)),
        'nonzero': sorted(nonzero, key=lambda x: -abs(x)),
    }

# Spectral invariant 1: Tr(D^2)/dim
print("\n  --- Invariant 1: <lambda^2> = Tr(D^2)/sub_dim ---")
full_avg = tr2_30 / total_dim  # = 8/30
print(f"    Full graph: {tr2_30:.4f}/{total_dim} = {full_avg:.6f}")
for name in leg_names:
    d = leg_data[name]
    avg = d['tr2'] / d['sub_dim']
    ratio_to_full = avg / full_avg
    print(f"    {name:7s}: {d['tr2']:.4f}/{d['sub_dim']} = {avg:.6f}  "
          f"(ratio to full: {ratio_to_full:.6f})")

# Spectral invariant 2: n_edges / sub_dim
print("\n  --- Invariant 2: n_edges/sub_dim ---")
for name in leg_names:
    d = leg_data[name]
    val = d['n_edges'] / d['sub_dim']
    print(f"    {name:7s}: {d['n_edges']}/{d['sub_dim']} = {val:.6f}")

# Spectral invariant 3: n_edges / dim_leg (total irrep dim of leg nodes)
print("\n  --- Invariant 3: n_edges/dim_leg (leg without branch) ---")
for li, (leg, name) in enumerate(zip(legs, leg_names)):
    dim_leg = sum(dims[n] for n in leg)
    n_edges = leg_data[name]['n_edges']
    val = n_edges / dim_leg if dim_leg > 0 else 0
    print(f"    {name:7s}: {n_edges}/{dim_leg} = {val:.6f}")

# Spectral invariant 4: spectral zeta at s=2
print("\n  --- Invariant 4: Spectral zeta zeta_leg(2) = sum(1/lambda^2) ---")
for name in leg_names:
    d = leg_data[name]
    nonzero = d['nonzero']
    if len(nonzero) > 0:
        zeta2 = sum(1.0/e**2 for e in nonzero)
        zeta2_norm = zeta2 / len(nonzero)
        print(f"    {name:7s}: zeta(2) = {zeta2:.6f}, "
              f"zeta(2)/n_nonzero = {zeta2_norm:.6f}")

# Spectral invariant 5: ratios of largest eigenvalues
print("\n  --- Invariant 5: Largest |eigenvalue| ratios ---")
full_max = max(abs(e) for e in eigs_30 if abs(e) > 1e-10)
for name in leg_names:
    d = leg_data[name]
    if d['nonzero']:
        leg_max = max(abs(e) for e in d['nonzero'])
        ratio = leg_max / full_max
        print(f"    {name:7s}: lambda_max = {leg_max:.6f}, "
              f"ratio to full = {ratio:.6f}")

# Spectral invariant 6: Tr(D^2) fraction
print("\n  --- Invariant 6: Tr(D^2)_leg / Tr(D^2)_full ---")
for name in leg_names:
    d = leg_data[name]
    frac = d['tr2'] / tr2_30
    print(f"    {name:7s}: {d['tr2']:.4f}/{tr2_30:.4f} = {frac:.6f}")

# Spectral invariant 7: dim fraction
print("\n  --- Invariant 7: dim_leg / dim_total ---")
for name in leg_names:
    d = leg_data[name]
    frac = d['sub_dim'] / total_dim
    print(f"    {name:7s}: {d['sub_dim']}/{total_dim} = {frac:.6f}")


# ============================================================
# SECTION 5: Direct comparison matrix
# ============================================================
print("\n" + "=" * 70)
print("SECTION 5: Direct Comparison - All Invariants vs Targets")
print("=" * 70)

targets = {
    'alpha_s': alpha_s_fw,
    'sin2tW': sin2tW_fw,
    '1/a1': inv_a1,
    '1/phi^4': inv_phi4,
    'alpha_s*2/pi': alpha_s_fw * 2 / np.pi,
}

print(f"\n  {'Invariant':>30} | {'Short':>10} | {'Medium':>10} | {'Long':>10} |")
print(f"  {'-'*30}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+")

# Compute many invariants
invariants = {}
for name in leg_names:
    d = leg_data[name]
    inv = {}

    # Basic
    inv['n_edges/sub_dim'] = d['n_edges'] / d['sub_dim']
    inv['Tr(D^2)/sub_dim'] = d['tr2'] / d['sub_dim']
    inv['Tr(D^2)/Tr_full'] = d['tr2'] / tr2_30

    # Eigenvalue-based (if nonzero exist)
    if d['nonzero']:
        inv['lambda_max'] = max(abs(e) for e in d['nonzero'])
        inv['lambda_min'] = min(abs(e) for e in d['nonzero'])
        inv['1/lambda_max'] = 1.0 / inv['lambda_max']
        inv['1/lambda_max^2'] = 1.0 / inv['lambda_max']**2
        inv['gap/lambda_max'] = (inv['lambda_max'] - inv['lambda_min']) / inv['lambda_max']

        # Spectral zeta
        inv['zeta(2)/n'] = sum(1.0/e**2 for e in d['nonzero']) / len(d['nonzero'])

        # Sum of positive eigenvalues
        pos_eigs = [e for e in d['nonzero'] if e > 0]
        neg_eigs = [e for e in d['nonzero'] if e < 0]
        if pos_eigs:
            inv['sum_pos/n_pos'] = sum(pos_eigs) / len(pos_eigs)
        if neg_eigs:
            inv['sum_neg/n_neg'] = sum(neg_eigs) / len(neg_eigs)

    # Leg-specific (without branch)
    leg_only_dim = sum(dims[n] for n in legs[leg_names.index(name)])
    inv['leg_dim/total_dim'] = leg_only_dim / total_dim

    invariants[name] = inv

for inv_name in sorted(set().union(*[set(v.keys()) for v in invariants.values()])):
    vals = []
    for name in leg_names:
        v = invariants[name].get(inv_name, float('nan'))
        vals.append(v)
    print(f"  {inv_name:>30} | {vals[0]:10.6f} | {vals[1]:10.6f} | {vals[2]:10.6f} |")

# Now check which invariant best matches which target
print(f"\n\n  MATCHES (looking for target values in invariant differences):")
print(f"  {'Target':>15} = {'Value':>10}")
for t_name, t_val in sorted(targets.items()):
    print(f"  {t_name:>15} = {t_val:10.6f}")

print(f"\n  Cross-leg DIFFERENCES and RATIOS:")
for inv_name in sorted(set().union(*[set(v.keys()) for v in invariants.values()])):
    vals = {}
    for name in leg_names:
        vals[name] = invariants[name].get(inv_name, float('nan'))
    if any(np.isnan(v) for v in vals.values()):
        continue

    # Differences
    diffs = {
        'L-S': vals['Long'] - vals['Short'],
        'L-M': vals['Long'] - vals['Medium'],
        'M-S': vals['Medium'] - vals['Short'],
    }
    # Ratios
    ratios = {}
    for k1 in leg_names:
        for k2 in leg_names:
            if k1 != k2 and abs(vals[k2]) > 1e-15:
                ratios[f'{k1[0]}/{k2[0]}'] = vals[k1] / vals[k2]

    # Check if any diff or ratio matches a target
    for label, val in list(diffs.items()) + list(ratios.items()):
        for t_name, t_val in targets.items():
            if abs(t_val) > 1e-15 and abs(val - t_val) / abs(t_val) < 0.05:
                pct = abs(val - t_val) / abs(t_val) * 100
                print(f"    ** {inv_name}[{label}] = {val:.6f} ~ {t_name} = {t_val:.6f} "
                      f"({pct:.2f}%)")
            # Also check negative
            if abs(t_val) > 1e-15 and abs(val + t_val) / abs(t_val) < 0.05:
                pct = abs(val + t_val) / abs(t_val) * 100
                print(f"    ** {inv_name}[{label}] = {val:.6f} ~ -{t_name} = {-t_val:.6f} "
                      f"({pct:.2f}%)")


# ============================================================
# SECTION 6: Normalized sub-graph spectra relative to 1/phi^4
# ============================================================
print("\n" + "=" * 70)
print("SECTION 6: Sub-Graph Spectra Relative to 1/phi^4 Baseline")
print("=" * 70)

print(f"\n  Baseline: 1/phi^4 = {inv_phi4:.6f}")
print(f"\n  For each leg, compute: delta = (spectral_quantity) - 1/phi^4")

for name in leg_names:
    d = leg_data[name]
    print(f"\n  {name} leg:")

    # Several spectral quantities to subtract baseline from
    quantities = {}
    quantities['Tr(D^2)/sub_dim'] = d['tr2'] / d['sub_dim']
    quantities['n_edges/sub_dim'] = d['n_edges'] / d['sub_dim']

    if d['nonzero']:
        quantities['1/lambda_max'] = 1.0 / max(abs(e) for e in d['nonzero'])
        quantities['1/lambda_max^2'] = 1.0 / max(abs(e) for e in d['nonzero'])**2
        pos = [e for e in d['nonzero'] if e > 0]
        if pos:
            quantities['1/sum_pos'] = 1.0 / sum(pos)
            quantities['n_pos/sum_pos^2'] = len(pos) / sum(pos)**2

    for q_name, q_val in quantities.items():
        delta = q_val - inv_phi4
        print(f"    {q_name:>25} = {q_val:.6f}, delta = {delta:+.6f}")
        # Check if delta matches any target
        for t_name, t_val in targets.items():
            if abs(t_val) > 1e-15 and abs(delta) > 1e-15:
                if abs(delta - t_val) / abs(t_val) < 0.10:
                    pct = abs(delta - t_val) / abs(t_val) * 100
                    print(f"      ** delta ~ {t_name} = {t_val:.6f} ({pct:.1f}%)")
                if abs(delta + t_val) / abs(t_val) < 0.10:
                    pct = abs(delta + t_val) / abs(t_val) * 100
                    print(f"      ** delta ~ -{t_name} = {-t_val:.6f} ({pct:.1f}%)")


# ============================================================
# SECTION 7: Dimension-weighted Laplacian approach
# ============================================================
print("\n" + "=" * 70)
print("SECTION 7: Dimension-Weighted Graph Laplacian")
print("=" * 70)

print("""
  The graph Laplacian L = D_deg - A captures diffusion on the McKay graph.
  With dimension weights: L[i,j] = -sqrt(d_i*d_j) if connected,
                          L[i,i] = sum of adjacent weights.
  Eigenvalues of L give the "harmonic frequencies" of each sub-graph.
""")

# Build 9x9 Laplacian (unit weights)
L9_unit = -A9.copy()
for i in range(9):
    L9_unit[i,i] = sum(A9[i,:])
eigs_L9 = np.sort(np.linalg.eigvalsh(L9_unit))
print(f"  Unit-weight Laplacian eigenvalues: {np.round(eigs_L9, 6).tolist()}")
print(f"  Smallest nonzero (Fiedler): {min(e for e in eigs_L9 if e > 1e-10):.6f}")

# Build 9x9 Laplacian (dim weights)
L9_dim = -D9_dimwt.copy()
for i in range(9):
    L9_dim[i,i] = sum(D9_dimwt[i,:])
eigs_L9d = np.sort(np.linalg.eigvalsh(L9_dim))
print(f"\n  Dim-weight Laplacian eigenvalues: {np.round(eigs_L9d, 4).tolist()}")

# Sub-graph Laplacians
print("\n  Sub-graph Laplacians (with branch):")
for li, (leg, name) in enumerate(zip(legs, leg_names)):
    nodes = sorted([branch_node] + list(leg))
    sub_L_unit = L9_unit[np.ix_(nodes, nodes)]
    sub_L_dim = L9_dim[np.ix_(nodes, nodes)]

    eigs_sub_u = np.sort(np.linalg.eigvalsh(sub_L_unit))
    eigs_sub_d = np.sort(np.linalg.eigvalsh(sub_L_dim))

    print(f"\n    {name} leg ({len(nodes)} nodes):")
    print(f"      Unit Laplacian eigs: {np.round(eigs_sub_u, 6).tolist()}")
    print(f"      Dim Laplacian eigs:  {np.round(eigs_sub_d, 4).tolist()}")

    # Fiedler value (algebraic connectivity)
    nonzero_u = [e for e in eigs_sub_u if e > 1e-10]
    if nonzero_u:
        fiedler_u = min(nonzero_u)
        print(f"      Fiedler (unit): {fiedler_u:.6f}")
        # Normalized Fiedler
        print(f"      Fiedler/n_nodes: {fiedler_u/len(nodes):.6f}")


# ============================================================
# SECTION 8: Kac label analysis
# ============================================================
print("\n" + "=" * 70)
print("SECTION 8: Kac Label / Dimension Ratios per Leg")
print("=" * 70)

print("\n  The Kac labels ARE the irrep dimensions.")
print(f"  sum(d_i) = {sum(dims)} = h(E8) = {h_E8}")
print(f"  sum(d_i^2) = {sum(d**2 for d in dims)} = |2I| = {N}")

for li, (leg, name) in enumerate(zip(legs, leg_names)):
    leg_dims = [dims[n] for n in leg]
    s = sum(leg_dims)
    s2 = sum(d**2 for d in leg_dims)
    print(f"\n  {name} leg: dims = {leg_dims}")
    print(f"    sum = {s}, sum^2 = {s2}")
    print(f"    sum/h(E8) = {s}/{h_E8} = {s/h_E8:.6f}")
    print(f"    sum^2/|2I| = {s2}/{N} = {s2/N:.6f}")
    print(f"    sum/(sum+b1) = {s}/{s+b1} = {s/(s+b1):.6f}")
    print(f"    b1/(sum+a1^2+1) = {b1}/{s+a1**2+1} = {b1/(s+a1**2+1):.6f}")

# Kac label sums per leg (without branch)
leg_sums = [sum(dims[n] for n in leg) for leg in legs]
print(f"\n  Leg dimension sums (without branch): {leg_sums}")
print(f"  Leg dim sums + branch: {[s + dims[branch_node] for s in leg_sums]}")
print(f"  Total (3 legs + branch): {sum(leg_sums) + dims[branch_node]} "
      f"(= {sum(dims)})")

# Key ratios
print(f"\n  Key ratios:")
for li, (s, name) in enumerate(zip(leg_sums, leg_names)):
    for lj, (s2_val, name2) in enumerate(zip(leg_sums, leg_names)):
        if li < lj:
            print(f"    {name}/{name2} = {s}/{s2_val} = {s/s2_val:.6f}")

# Ratios involving a1, b1
print(f"\n  Ratios with framework constants:")
for s, name in zip(leg_sums, leg_names):
    print(f"    {name}: dim_sum = {s}")
    print(f"      s/a1^2 = {s/a1**2:.6f}")
    print(f"      s/N = {s/N:.6f}")
    print(f"      b1/s = {b1/s:.6f}" if s > 0 else "      (empty)")
    print(f"      1/(s+1) = {1/(s+1):.6f}")
    print(f"      s/(a1*b1) = {s/(a1*b1):.6f}")


# ============================================================
# SECTION 9: Heat kernel at t = phi^4
# ============================================================
print("\n" + "=" * 70)
print("SECTION 9: Heat Kernel K(t) = Tr(exp(-t*D^2)) at Special Values")
print("=" * 70)

# Full D_F heat kernel
eigs2_full = np.array([e**2 for e in eigs_30])

for t_name, t_val in [('1', 1.0), ('phi', PHI), ('phi^2', PHI**2),
                       ('phi^4', PHI**4), ('1/phi^4', inv_phi4),
                       ('2*phi^3', 2*PHI**3), ('a1', float(a1))]:
    K = sum(np.exp(-t_val * eigs2_full))
    print(f"  K({t_name:>8} = {t_val:8.4f}) = {K:.6f} / {total_dim}")

# Sub-graph heat kernels
print(f"\n  Sub-graph heat kernels at t = phi^4 = {PHI**4:.4f}:")
t_phi4 = PHI**4
for name in leg_names:
    d = leg_data[name]
    eigs2_leg = np.array([e**2 for e in d['eigenvalues']])
    K_leg = sum(np.exp(-t_phi4 * eigs2_leg))
    K_norm = K_leg / d['sub_dim']
    print(f"    {name:7s}: K = {K_leg:.6f}, K/dim = {K_norm:.6f}")

print(f"\n  Sub-graph heat kernels at t = 2*phi^3 = {2*PHI**3:.4f}:")
t_2phi3 = 2 * PHI**3
for name in leg_names:
    d = leg_data[name]
    eigs2_leg = np.array([e**2 for e in d['eigenvalues']])
    K_leg = sum(np.exp(-t_2phi3 * eigs2_leg))
    K_norm = K_leg / d['sub_dim']
    print(f"    {name:7s}: K = {K_leg:.6f}, K/dim = {K_norm:.6f}")


# ============================================================
# SECTION 10: Summary
# ============================================================
print("\n" + "=" * 70)
print("SECTION 10: Summary and Assessment")
print("=" * 70)

print("""
  QUESTION: Do sub-graph eigenvalues of D_F on affine E8 produce
  the sector-dependent corrections alpha_s, sin^2(tW), 1/a1?

  METHODOLOGY:
  - Built D_F in two versions: 9x9 adjacency, 30x30 CG-based
  - Decomposed into 3 legs (short/medium/long) with and without branch
  - Computed ~15 spectral invariants per leg
  - Searched for matches with alpha_s, sin^2(tW), 1/a1, 1/phi^4

  KEY FINDINGS:
""")

# Automated summary of best matches
all_matches = []
for inv_name in sorted(set().union(*[set(v.keys()) for v in invariants.values()])):
    vals = {}
    for name in leg_names:
        vals[name] = invariants[name].get(inv_name, float('nan'))
    if any(np.isnan(v) for v in vals.values()):
        continue
    for label_pair in [('Long', 'Short'), ('Long', 'Medium'), ('Medium', 'Short')]:
        diff = vals[label_pair[0]] - vals[label_pair[1]]
        ratio = vals[label_pair[0]] / vals[label_pair[1]] if abs(vals[label_pair[1]]) > 1e-15 else float('inf')
        for t_name, t_val in targets.items():
            if abs(t_val) > 1e-15:
                for val, vtype in [(diff, 'diff'), (ratio, 'ratio')]:
                    for sign in [1, -1]:
                        if abs(sign*val - t_val) / abs(t_val) < 0.05:
                            pct = abs(sign*val - t_val) / abs(t_val) * 100
                            s = '+' if sign > 0 else '-'
                            all_matches.append((pct, f"{s}{inv_name}[{label_pair[0][0]}-{label_pair[1][0]}] "
                                               f"({vtype}) = {sign*val:.6f} ~ {t_name} = {t_val:.6f} "
                                               f"({pct:.2f}%)"))

if all_matches:
    all_matches.sort()
    print("  Best matches (<5% error):")
    for pct, desc in all_matches[:20]:
        print(f"    {desc}")
else:
    print("  No matches found within 5% tolerance.")

print(f"""
  CONCLUSION:
  The McKay graph D_F eigenvalues provide a rich spectral structure.
  Whether they produce the EXACT sector-dependent corrections
  alpha_s, sin^2(tW), 1/a1 is assessed above.

  CATEGORY: [determined by results]
""")

print("\n" + "=" * 70)
print("EXP-315 COMPLETE")
print("=" * 70)

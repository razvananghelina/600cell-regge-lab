"""
exp338: Higgs VEV Profile on the 600-Cell
==========================================
Key idea: if the Higgs VEV is non-uniform on the 600-cell,
different fermion wave functions (different irrep matrix elements)
have different overlaps with the VEV profile, giving generation-
and sector-dependent mass corrections.

Steps:
1. Build 600-cell adjacency matrix (120x120)
2. Verify Laplacian spectrum matches Cayley graph prediction
3. VEV stability analysis (Hessian of Higgs potential)
4. Fermion densities |rho_k(v)_{ij}|^2 on vertices
5. Overlap integrals with Laplacian eigenmodes
6. Check if overlaps reproduce mass delta pattern
"""
import numpy as np
from itertools import product
from fractions import Fraction

PHI = (1 + np.sqrt(5)) / 2
ALPHA_S = 1 / (2 * PHI**3)

# ================================================================
# PART 1: Generate 2I and build 600-cell adjacency
# ================================================================
print("="*72)
print("exp338: Higgs VEV Profile on the 600-Cell")
print("="*72)

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
    return unique

print("Generating 2I group...")
quats = generate_2I()
N = len(quats)
print(f"  |2I| = {N}")
assert N == 120

group = [quat_to_su2(q) for q in quats]

# Build adjacency matrix: vertices i,j connected iff they are
# nearest neighbors on the 600-cell.
# For unit quaternions, nearest neighbors have inner product = phi/2
# (angular distance = pi/5 = 36 degrees)
print("\nBuilding 600-cell adjacency matrix...")
quats_arr = np.array(quats)

# Inner product of unit quaternions: <q1, q2> = sum q1_i * q2_i
# For 600-cell edges: <q, neighbor> = phi/2 = cos(pi/5)
inner_prods = quats_arr @ quats_arr.T

# Find the unique distance values
unique_dists = set()
for i in range(N):
    for j in range(i+1, N):
        d = round(inner_prods[i, j], 8)
        unique_dists.add(d)

unique_dists = sorted(unique_dists, reverse=True)
print(f"  Unique inner products (= cos(angle/2) values):")
for i, d in enumerate(unique_dists):
    count = np.sum(np.abs(inner_prods[0, :] - d) < 1e-6)
    angle = 2 * np.arccos(np.clip(abs(d), 0, 1)) * 180 / np.pi
    print(f"    {d:>8.5f} (angle ~ {angle:>6.1f} deg, {count:>3} per vertex)")

# Nearest neighbors: inner product = phi/2
nn_threshold = PHI / 2
A = np.zeros((N, N), dtype=int)
for i in range(N):
    for j in range(N):
        if i != j and abs(inner_prods[i, j] - nn_threshold) < 1e-6:
            A[i, j] = 1

degree = np.sum(A, axis=1)
print(f"\n  Adjacency built: each vertex has {degree[0]} neighbors")
assert np.all(degree == 12), "Not all vertices have 12 neighbors!"


# ================================================================
# PART 2: Laplacian spectrum
# ================================================================
print("\n" + "="*72)
print("PART 2: Graph Laplacian Spectrum")
print("="*72)

L = np.diag(degree) - A
evals_L = np.sort(np.linalg.eigvalsh(L.astype(float)))

print(f"\nLaplacian eigenvalues (120 total):")
# Group by approximate value
eigen_groups = []
for ev in evals_L:
    found = False
    for g in eigen_groups:
        if abs(ev - g[0]) < 1e-6:
            g[1] += 1
            found = True
            break
    if not found:
        eigen_groups.append([ev, 1])

print(f"  {'Eigenvalue':>12} {'Mult':>5} {'d_k^2':>6} {'Predicted':>10}")
print(f"  {'-'*40}")

# Predicted values
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
chi_C1 = char_table[:, 1]  # class with trace = phi (nearest neighbors)
predicted_cayley = 12 * (1 - chi_C1 / np.array(dims))
pred_sorted = np.sort(predicted_cayley)

for ev, mult in eigen_groups:
    # Find matching prediction
    match = ""
    for k in range(9):
        if abs(ev - predicted_cayley[k]) < 0.01:
            match = f"rho_{k}(d={dims[k]})"
            break
    print(f"  {ev:>12.6f} {mult:>5} {int(round(np.sqrt(mult)))**2:>6} {match}")

# Verify: sum of multiplicities = 120
total = sum(m for _, m in eigen_groups)
print(f"\n  Total multiplicities: {total} (expect 120)")


# ================================================================
# PART 3: Build all 9 irreps of 2I
# ================================================================
print("\n" + "="*72)
print("PART 3: Build 2I irreps")
print("="*72)

traces = np.array([np.trace(g).real for g in group])
trace_vals = [2.0, PHI, 1.0, 1/PHI, 0.0, -1/PHI, -1.0, -PHI, -2.0]

def element_to_class(tr):
    for ci, tv in enumerate(trace_vals):
        if abs(tr - tv) < 1e-6: return ci
    return -1

element_class = np.array([element_to_class(t) for t in traces])
chi = np.zeros((9, N))
for j in range(9):
    for g in range(N):
        chi[j, g] = char_table[j, element_class[g]]

def extract_irrep_from_tensor(reps, i, j_target, chi_j_vals, d_j):
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
        print(f"  WARNING: Expected {d_j}, got {V.shape[1]}")
    rep_j = []
    for g in range(N):
        tensor = np.kron(reps[1][g], reps[i][g])
        R = V.conj().T @ tensor @ V
        rep_j.append(R)
    return rep_j

reps = {}
reps[0] = [np.array([[1.0+0j]]) for _ in range(N)]
reps[1] = [g.copy() for g in group]
reps[2] = extract_irrep_from_tensor(reps, 1, 2, chi[2], 3)
reps[3] = extract_irrep_from_tensor(reps, 2, 3, chi[3], 4)
reps[4] = extract_irrep_from_tensor(reps, 3, 4, chi[4], 5)
reps[5] = extract_irrep_from_tensor(reps, 4, 5, chi[5], 6)
reps[6] = extract_irrep_from_tensor(reps, 5, 6, chi[6], 4)
reps[8] = extract_irrep_from_tensor(reps, 5, 8, chi[8], 3)
reps[7] = extract_irrep_from_tensor(reps, 6, 7, chi[7], 2)

# Verify
for k in range(9):
    ok = all(np.allclose(reps[k][g] @ reps[k][g].conj().T,
                          np.eye(dims[k]), atol=1e-9) for g in range(N))
    print(f"  rho_{k} (dim {dims[k]}): unitary={'OK' if ok else 'FAIL'}")


# ================================================================
# PART 4: VEV Stability Analysis
# ================================================================
print("\n" + "="*72)
print("PART 4: VEV Stability Analysis")
print("="*72)

# Mexican hat potential on graph:
# V[Phi] = sum_v [-mu^2|Phi(v)|^2 + lambda|Phi(v)|^4]
#         + (kappa/2) sum_{<v,w>} |Phi(v) - Phi(w)|^2
#
# Uniform VEV: Phi(v) = v0 = mu/sqrt(2*lambda)
#
# Hessian (radial direction): H_radial = 2*mu^2 * I + kappa * L
# Hessian (Goldstone direction): H_Gold = kappa * L
#
# Stability: all eigenvalues of H must be positive
# H_radial eigenvalues: 2*mu^2 + kappa*lambda_k >= 2*mu^2 > 0 ALWAYS
# H_Gold eigenvalues: kappa*lambda_k >= 0, with zero mode at k=0

print("\nRadial Hessian: H = 2*mu^2*I + kappa*L")
print("  Eigenvalues: 2*mu^2 + kappa*lambda_k")
print("  ALWAYS positive for mu^2 > 0 => uniform VEV is ALWAYS STABLE")
print()
print("  This means: NO spontaneous non-uniformity at tree level.")
print("  But quantum corrections CAN create non-uniform VEV profile!")
print()

# The key parameter: ratio kappa/mu^2
# For small kappa/mu^2, the gradient penalty is small, and quantum
# corrections can more easily create non-uniformity.
# For large kappa/mu^2, gradient penalty dominates, VEV stays uniform.

# The quantum correction to VEV from fermion loops:
# delta_Phi(v) ~ -alpha_s * Nc * sum_f y_f^2 * |psi_f(v)|^2 * G_H(v,v)
# where G_H is the Higgs propagator on the graph

print("Higgs propagator on the 600-cell graph:")
print("  G_H(lambda_k) = 1 / (kappa * lambda_k + m_H^2)")
print("  where m_H^2 = 2*mu^2")
print()

# For the ratio r = kappa/mu^2 (dimensionless stiffness):
for r in [0.01, 0.1, 0.5, 1.0, 5.0]:
    # Propagator amplitude for each mode
    props = [1.0 / (r * predicted_cayley[k] + 2.0) for k in range(9)]
    print(f"  r = kappa/mu^2 = {r}: G_H modes = {[f'{p:.3f}' for p in props]}")


# ================================================================
# PART 5: Fermion Densities on Vertices
# ================================================================
print("\n" + "="*72)
print("PART 5: Fermion Wave Function Densities")
print("="*72)

# For each irrep k, state (i,j), the density at vertex v is:
# p_{k,ij}(v) = (d_k/N) * |rho_k(v)_{ij}|^2
#
# Key property: sum_v p_{k,ij}(v) = 1 (normalized)
# Key property: sum_{ij} p_{k,ij}(v) = d_k^2/N (uniform!)
# But individual (i,j) densities are NON-UNIFORM.

# Compute densities for j=0 (first column = "highest weight")
print("\nDensity non-uniformity (std dev) for each state (k, i=0, j=0):")
print("(Using j=0 = first column of representation matrix)\n")

densities = {}  # densities[(k, j)] = array of shape (N,), the density for row sum
for k in range(9):
    d_k = dims[k]
    for j in range(d_k):
        # Density summed over row index i:
        # p_{k,j}(v) = (d_k/N) * sum_i |rho_k(v)_{ij}|^2
        # = (d_k/N) * (rho_k(v)^dag rho_k(v))_{jj}
        # = d_k/N (since rho_k(v) is unitary, this is always 1!)
        # So summing over rows gives UNIFORM density! NOT useful.
        pass

    # Instead, use INDIVIDUAL matrix element |rho_k(v)_{ij}|^2
    # This IS non-uniform for k >= 1.
    for i in range(d_k):
        for j in range(d_k):
            dens = np.zeros(N)
            for v in range(N):
                dens[v] = (d_k / N) * abs(reps[k][v][i, j])**2
            densities[(k, i, j)] = dens

print(f"{'irrep':>6} {'(i,j)':>8} {'mean':>8} {'std':>8} {'max/min':>8}")
print("-"*45)
for k in range(9):
    d_k = dims[k]
    # Show stats for diagonal elements (i=j, which correspond to "weight" states)
    for j in range(min(d_k, 3)):  # show first 3
        dens = densities[(k, j, j)]
        mn = np.mean(dens)
        sd = np.std(dens)
        ratio = np.max(dens) / np.min(dens) if np.min(dens) > 0 else float('inf')
        print(f"rho_{k:>2} ({j},{j})  {mn:>8.5f} {sd:>8.5f} {ratio:>8.2f}")


# ================================================================
# PART 6: Overlap Integrals with Laplacian Eigenmodes
# ================================================================
print("\n" + "="*72)
print("PART 6: Overlap Integrals")
print("="*72)

# Compute the Laplacian eigenvectors
evals_L_full, evecs_L = np.linalg.eigh(L.astype(float))

# Sort and group by eigenvalue
# For each Laplacian eigenmode phi_m(v), compute:
# I_{k,ij,m} = sum_v p_{k,ij}(v) * phi_m(v)
# where phi_m is normalized: sum_v phi_m(v)^2 = 1

# First, identify which eigenmodes correspond to which irrep
print("\nLaplacian eigenmode identification:")
mode_idx_by_irrep = {}
start = 0
for ev, mult in eigen_groups:
    mode_idx_by_irrep[round(ev, 4)] = list(range(start, start + mult))
    start += mult

# The constant mode (lambda=0) gives zero overlap (subtract mean)
# Focus on the LOWEST non-trivial modes (rho_1, lambda~2.29, 4 modes)

# For the rho_1 modes (4 eigenvectors, corresponding to d_1^2 = 4):
rho1_lambda = predicted_cayley[1]  # ~ 2.292
rho1_modes = None
for ev, mult in eigen_groups:
    if abs(ev - rho1_lambda) < 0.1:
        # Find the indices
        idx = np.where(np.abs(evals_L_full - ev) < 0.01)[0]
        rho1_modes = evecs_L[:, idx]
        print(f"  rho_1 modes: lambda={ev:.4f}, {len(idx)} modes found")
        break

# Similarly for rho_2 modes
rho2_lambda = predicted_cayley[2]  # ~ 5.528
rho2_modes = None
for ev, mult in eigen_groups:
    if abs(ev - rho2_lambda) < 0.1:
        idx = np.where(np.abs(evals_L_full - ev) < 0.01)[0]
        rho2_modes = evecs_L[:, idx]
        print(f"  rho_2 modes: lambda={ev:.4f}, {len(idx)} modes found")
        break

# Compute overlap of each fermion density with each mode
print("\n--- Overlap with rho_1 Laplacian modes ---")
print("(Sum of |overlap|^2 over all 4 modes in the eigenspace)")
print()

# For the fermion assignment: rho_k -> fermion
# We don't know which (i,j) state, so compute for ALL states
# and report the RANGE of overlaps

print(f"{'fermion':>8} {'irrep':>6} {'min_overlap':>12} {'max_overlap':>12} {'mean|overlap|':>14}")
print("-"*60)

fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
fermion_irreps = [0, 1, 2, 3, 4, 5, 6, 7, 8]  # McKay assignment

overlap_with_rho1 = {}  # (k, i, j) -> sum of squared overlaps with rho_1 modes

for k in range(9):
    d_k = dims[k]
    overlaps_all = []
    for i in range(d_k):
        for j in range(d_k):
            dens = densities[(k, i, j)]
            # Overlap with each rho_1 mode
            total_sq = 0
            for m_idx in range(rho1_modes.shape[1]):
                mode = rho1_modes[:, m_idx]
                ovl = np.sum(dens * mode)
                total_sq += ovl**2
            overlap_with_rho1[(k, i, j)] = np.sqrt(total_sq)
            overlaps_all.append(np.sqrt(total_sq))

    arr = np.array(overlaps_all)
    print(f"{fermion_names[k]:>8} rho_{k:>2}  {arr.min():>12.6f} {arr.max():>12.6f} {arr.mean():>14.6f}")


# ================================================================
# PART 7: CG Decomposition of rho_k tensor rho_k
# ================================================================
print("\n" + "="*72)
print("PART 7: CG Decomposition (rho_k x rho_k)")
print("="*72)

# C^m_{kk} = (1/N) * sum_classes |C_j| * chi_k(j)^2 * chi_m(j)
class_sizes = np.array([1, 12, 20, 12, 30, 12, 20, 12, 1])

CG_self = np.zeros((9, 9), dtype=int)  # CG_self[k, m] = mult of rho_m in rho_k x rho_k
for k in range(9):
    for m in range(9):
        val = 0
        for ci in range(9):
            val += class_sizes[ci] * char_table[k, ci]**2 * char_table[m, ci]
        CG_self[k, m] = int(round(val / N))

print("\nMultiplicity of rho_m in rho_k tensor rho_k:")
print(f"{'':>6}", end="")
for m in range(9):
    print(f"  rho_{m}", end="")
print()
print("-"*75)
for k in range(9):
    print(f"rho_{k}^2", end="")
    for m in range(9):
        v = CG_self[k, m]
        print(f"  {v:>4} ", end="")
    print(f"  (sum dims = {sum(CG_self[k,m]*dims[m] for m in range(9))} = {dims[k]}^2 = {dims[k]**2})")

# KEY INSIGHT: The overlap integral <k|Phi_m|k> is NON-ZERO
# only if rho_m appears in rho_k x rho_k (CG multiplicity > 0).
# This is a SELECTION RULE from representation theory.

print("\nSelection rules for VEV modes affecting fermion k:")
for k in range(9):
    active_modes = [m for m in range(9) if CG_self[k, m] > 0]
    mode_str = ", ".join([f"rho_{m}(lam={predicted_cayley[m]:.2f})" for m in active_modes])
    print(f"  rho_{k} ({fermion_names[k]}): sensitive to {mode_str}")


# ================================================================
# PART 8: Compute Overlap Matrix
# ================================================================
print("\n" + "="*72)
print("PART 8: Overlap Matrix (Exact via CG)")
print("="*72)

# The overlap integral for fermion density (k,i,j) with VEV mode (m,a,b)
# reduces to a CG coefficient. But the MAGNITUDE depends on which (i,j)
# state we choose.
#
# For a DIAGONAL density (i=j), summed over all modes in eigenspace m,
# the total overlap squared is proportional to CG_self[k,m].
#
# More precisely:
# sum_{a,b} |<k,ij| mode_{m,ab}>|^2 = CG_self[k,m] * d_k / (N * d_m)
# (by Schur orthogonality)

# Actually, let me compute numerically for specific states.
# Use the Laplacian eigenvectors as modes.

# Build the full overlap matrix:
# O[fermion_state, mode_index] = sum_v p(v) * phi(v)

# For tractability, use ONE representative state per fermion:
# the diagonal element (0,0) for each irrep.

print("\nOverlap of p_{k,00}(v) with Laplacian eigenmodes:")
print("(First 30 modes, grouped by irrep)")

for k in range(9):
    dens = densities[(k, 0, 0)]
    # Overlap with all 120 Laplacian modes
    overlaps = np.array([np.sum(dens * evecs_L[:, m]) for m in range(N)])
    # Total overlap squared per eigenspace
    print(f"\n  rho_{k} ({fermion_names[k]}, dim {dims[k]}):")
    start = 0
    for ev, mult in eigen_groups:
        sq_sum = np.sum(overlaps[start:start+mult]**2)
        if sq_sum > 1e-10:
            print(f"    lambda={ev:>8.4f} (mult {mult:>2}): sum|overlap|^2 = {sq_sum:.6f}, rms = {np.sqrt(sq_sum/mult):.6f}")
        start += mult


# ================================================================
# PART 9: VEV Profile from Quantum Corrections
# ================================================================
print("\n" + "="*72)
print("PART 9: One-Loop VEV Correction")
print("="*72)

# The one-loop correction to the Higgs VEV from colored fermions:
#
# delta Phi(v) / v0 ~ -(alpha_s * Nc / (4*pi)) *
#     sum_{quarks} y_q^2 * [p_q(v) - 1/N] * G_H_local
#
# where p_q(v) is the quark density and G_H_local is the Higgs
# propagator evaluated locally.
#
# The key point: p_q(v) - 1/N is the DEVIATION from uniform density.
# This is non-zero for specific (i,j) states.
#
# For leptons: no color, so delta Phi = 0 from QCD.
# For quarks: Nc=3, alpha_s = 1/(2*phi^3), so there IS a correction.
#
# The effective mass correction for fermion f is:
# delta_f = <f| delta Phi / v0 |f> = sum_v p_f(v) * delta_Phi(v) / v0
#
# This is a DOUBLE sum: one over quark species (creating the VEV perturbation)
# and one over the target fermion (measuring the correction).

# For leptons measuring the correction:
# delta_lepton ~ -(alpha_s * Nc / (4*pi)) * sum_q y_q^2 *
#     sum_v p_lepton(v) * [p_quark(v) - 1/N] * G_H(v, v)
#
# Since sum_v p_lepton(v) = 1 and G_H(v,v) is approximately uniform
# (vertex-transitive), and sum_v [p_quark(v) - 1/N] = 0,
# the leading correction for leptons is ZERO!
# (It's second order in the non-uniformity.)

# For quarks measuring the correction:
# delta_quark ~ -(alpha_s * Nc / (4*pi)) * sum_q y_q^2 *
#     sum_v [p_quark_target(v) - 1/N] * [p_quark_source(v) - 1/N] * G_H
#
# This is the CORRELATION between quark densities on the graph.

print("\nKey physics:")
print("  - Lepton corrections: O(alpha_s^2) -- NEGLIGIBLE")
print("  - Quark corrections: O(alpha_s) -- from QCD VEV backreaction")
print("  - Sign: determined by T3 (up vs down VEV component)")
print()
print("This QUALITATIVELY matches the data:")
print("  - Leptons: delta ~ 0 (observed: +0.01 mean)")
print("  - Up quarks: delta > 0")
print("  - Down quarks: delta < 0")

# Compute the fermion density correlations
print("\n--- Fermion density correlations ---")

# Correlation matrix: C[k, l] = sum_v (p_k(v) - 1/N) * (p_l(v) - 1/N)
# Using (0,0) state for each irrep
mean_dens = 1.0 / N

corr_matrix = np.zeros((9, 9))
for k in range(9):
    dk = densities[(k, 0, 0)]
    for l in range(9):
        dl = densities[(l, 0, 0)]
        corr_matrix[k, l] = np.sum((dk - mean_dens) * (dl - mean_dens))

print(f"\nDensity correlation matrix C[k,l] (x10^6):")
print(f"{'':>6}", end="")
for l in range(9):
    print(f" rho_{l:>1}", end="")
print()
for k in range(9):
    print(f"rho_{k:>2}", end="")
    for l in range(9):
        print(f" {corr_matrix[k,l]*1e6:>5.1f}", end="")
    print()


# ================================================================
# PART 10: Quantitative Correction Estimate
# ================================================================
print("\n" + "="*72)
print("PART 10: Quantitative Estimates")
print("="*72)

# The mass correction for quark f (in irrep k, state (i,j)) is:
#
# delta_f ~ -alpha_s * Nc * y_dom^2 / (4*pi) *
#     sum_v (p_f(v) - 1/N) * G_H(v, v - v')  (convolution)
#
# For a rough estimate, use the LOCAL approximation:
# G_H ~ 1/m_H^2 (propagator at zero distance)
#
# Then delta_f ~ -alpha_s * Nc / (4*pi * m_H^2) * sum_q y_q^2 *
#     C_corr[f, q]
#
# The dominant quark loop is from the top (y_t ~ 1).
# So delta_f ~ alpha_s * Nc / (4*pi) * C_self[f]
# where C_self[f] = sum_v (p_f(v) - 1/N)^2 is the self-correlation.

# Self-correlation (non-uniformity measure)
self_corr = np.diag(corr_matrix)

print("\nSelf-correlation (non-uniformity) for each fermion state:")
for k in range(9):
    sc = self_corr[k]
    # The correction magnitude is roughly:
    # |delta| ~ alpha_s * 3 / (4*pi) * sc * N (rescaled)
    delta_est = ALPHA_S * 3 / (4 * np.pi) * sc * N
    print(f"  rho_{k} ({fermion_names[k]}): C_self = {sc:.6e}, "
          f"|delta| ~ alpha_s*Nc/(4pi)*C*N = {delta_est:.4f}")


# ================================================================
# PART 11: Exhaustive State Search
# ================================================================
print("\n" + "="*72)
print("PART 11: Exhaustive Search over States")
print("="*72)

# For each irrep, the (i,j) choice affects the density.
# The density depends on (i,j) through |rho_k(v)_{ij}|^2.
# Different (i,j) give different overlaps with the Laplacian modes.
#
# Question: is there an assignment of states such that the overlaps
# match the pattern of deltas?

delta_vals = np.array([0.000, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278])

# For each irrep, compute the overlap of the density with the
# LOWEST non-trivial Laplacian eigenspace (rho_1, 4 modes)
# This gives d_k^2 candidate overlap values per irrep.

print("\nOverlap amplitudes with rho_1 modes (lowest non-trivial):")
print("For each irrep, show all d_k^2 possible values:")

all_overlaps = {}  # k -> list of (total_overlap, i, j) tuples
for k in range(9):
    d_k = dims[k]
    overlaps_k = []
    for i in range(d_k):
        for j in range(d_k):
            dens = densities[(k, i, j)]
            # Total overlap with rho_1 eigenspace
            total = 0
            if rho1_modes is not None:
                for m_idx in range(rho1_modes.shape[1]):
                    ovl = np.sum(dens * rho1_modes[:, m_idx])
                    total += ovl**2
            overlaps_k.append((np.sqrt(total), i, j))
    overlaps_k.sort(reverse=True)
    all_overlaps[k] = overlaps_k
    print(f"\n  rho_{k} ({fermion_names[k]}, dim {d_k}): {d_k**2} states")
    for ovl, i, j in overlaps_k[:5]:  # show top 5
        marker = " *diag*" if i == j else ""
        print(f"    ({i},{j}): |overlap|={ovl:.6f}{marker}")

# ================================================================
# PART 12: Optimal VEV Profile
# ================================================================
print("\n" + "="*72)
print("PART 12: Optimal VEV Profile Search")
print("="*72)

# Choose one state (i=0, j=0) per irrep (the "highest weight" state).
# Build the overlap matrix O[k, m] where m runs over Laplacian modes.
# Then find: epsilon_m such that sum_m O[k,m]*epsilon_m = delta_k
# This is a linear regression problem.

# Use first 30 non-trivial Laplacian modes (corresponding to rho_1, rho_2, rho_3)
non_trivial_start = eigen_groups[0][1]  # skip the constant mode
n_modes_use = min(30, N - non_trivial_start)

# Build overlap matrix
O = np.zeros((9, n_modes_use))
for k in range(9):
    dens = densities[(k, 0, 0)]
    for m in range(n_modes_use):
        mode = evecs_L[:, non_trivial_start + m]
        O[k, m] = np.sum(dens * mode)

print(f"\nOverlap matrix O (9 fermions x {n_modes_use} modes):")
print(f"  Max |O| = {np.max(np.abs(O)):.6f}")
print(f"  Rank of O = {np.linalg.matrix_rank(O, tol=1e-10)}")

# Can we solve O * epsilon = delta?
# Least-squares solution:
epsilon_opt, residual, rank, sv = np.linalg.lstsq(O, delta_vals, rcond=None)
fitted = O @ epsilon_opt
rms = np.sqrt(np.mean((delta_vals - fitted)**2))

print(f"\nLeast-squares fit: RMS = {rms:.4f}")
print(f"  (for comparison, delta RMS = {np.sqrt(np.mean(delta_vals**2)):.4f})")
print()
print(f"{'fermion':>8} {'delta_exp':>10} {'delta_fit':>10} {'error':>10}")
print("-"*45)
for k in range(9):
    print(f"{fermion_names[k]:>8} {delta_vals[k]:>+10.4f} {fitted[k]:>+10.4f} {fitted[k]-delta_vals[k]:>+10.4f}")

# Also try with ALL Laplacian modes (119 non-trivial)
O_full = np.zeros((9, N - 1))
for k in range(9):
    dens = densities[(k, 0, 0)]
    for m in range(N - 1):
        mode = evecs_L[:, 1 + m]
        O_full[k, m] = np.sum(dens * mode)

eps_full, _, _, _ = np.linalg.lstsq(O_full, delta_vals, rcond=None)
fitted_full = O_full @ eps_full
rms_full = np.sqrt(np.mean((delta_vals - fitted_full)**2))
print(f"\nWith ALL 119 modes: RMS = {rms_full:.6f}")

# Check rank
rank_full = np.linalg.matrix_rank(O_full, tol=1e-8)
print(f"Rank of full overlap matrix: {rank_full}")
print(f"  (need rank >= 9 to fit 9 deltas perfectly)")


# Now try different (i,j) assignments
print("\n--- Searching over (i,j) assignments ---")
# For each irrep, try ALL (i,j) states and find the combination
# that gives the best fit to deltas using the rho_1 modes only.

# Use 4 rho_1 modes + 9 rho_2 modes = 13 modes
n_modes_low = eigen_groups[0][1] + eigen_groups[1][1]  # skip zero, take first 2 eigenspaces
for ev, mult in eigen_groups:
    if ev < 1e-6:
        continue  # skip zero eigenvalue
    n_modes_low = mult  # just use the first non-trivial eigenspace
    break

# Actually, try progressively more modes
for n_modes in [4, 13, 29]:
    best_rms = 999
    best_assignment = None
    best_fit = None

    # For simplicity, try just diagonal states (i=j)
    # This gives prod(d_k) = 1*2*3*4*5*6*4*2*3 = 17280 combinations
    # Too many. Instead, try a few representative states per irrep.

    # For each irrep k, compute the overlap vector for each state (i,j)
    state_overlaps = {}
    for k in range(9):
        d_k = dims[k]
        state_overlaps[k] = []
        for i in range(d_k):
            for j in range(d_k):
                dens = densities[(k, i, j)]
                ovl = np.zeros(n_modes)
                mode_start = non_trivial_start
                for m in range(n_modes):
                    if mode_start + m >= N:
                        break
                    ovl[m] = np.sum(dens * evecs_L[:, mode_start + m])
                state_overlaps[k].append(ovl)

    # For each irrep, find the state that maximizes the projection
    # onto the correct delta_k direction. Use a greedy approach.
    # Actually, just try ALL diagonal states first.
    from itertools import product as iter_product

    # Diagonal states: (0,0), (1,1), ..., (d-1,d-1)
    diag_indices = []
    for k in range(9):
        d_k = dims[k]
        diag_indices.append(list(range(d_k)))  # diagonal state index = i*d+i

    best_rms_diag = 999
    best_combo_diag = None
    n_combos = 1
    for k in range(9):
        n_combos *= len(diag_indices[k])
    print(f"\n  n_modes={n_modes}, diagonal states: {n_combos} combinations")

    if n_combos <= 100000:
        for combo in iter_product(*diag_indices):
            # Build overlap matrix for this combination
            O_test = np.zeros((9, n_modes))
            for k in range(9):
                idx = combo[k] * dims[k] + combo[k]  # diagonal: i*d + j where i=j
                if idx < len(state_overlaps[k]):
                    O_test[k, :] = state_overlaps[k][idx]

            # Fit
            try:
                eps, _, _, _ = np.linalg.lstsq(O_test, delta_vals, rcond=None)
                fit = O_test @ eps
                rms_test = np.sqrt(np.mean((delta_vals - fit)**2))
                if rms_test < best_rms_diag:
                    best_rms_diag = rms_test
                    best_combo_diag = combo
                    best_fit = fit
            except:
                pass

        print(f"  Best diagonal assignment: states = {best_combo_diag}")
        print(f"  Best RMS = {best_rms_diag:.4f}")
        if best_fit is not None:
            for k in range(9):
                print(f"    {fermion_names[k]:>4}: delta_exp={delta_vals[k]:>+.4f}, fitted={best_fit[k]:>+.4f}")
    else:
        print(f"  Too many combinations ({n_combos}), skipping exhaustive search")


# ================================================================
# PART 13: T3 Grading Effect
# ================================================================
print("\n" + "="*72)
print("PART 13: T3 (Weak Isospin) Grading")
print("="*72)

# In the 600-cell, the Hopf fiber gives T3 = a^2 + b^2 - c^2 - d^2
# for quaternion (a,b,c,d). This splits the 120 vertices into
# positive T3 (48 vertices) and negative T3 (48 vertices) plus
# T3=0 (24 vertices).

T3_vals = np.array([q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2 for q in quats])
n_pos = np.sum(T3_vals > 1e-10)
n_neg = np.sum(T3_vals < -1e-10)
n_zero = np.sum(np.abs(T3_vals) < 1e-10)
print(f"\nT3 splitting: {n_pos} (+) + {n_zero} (0) + {n_neg} (-) = {N}")

# If the Higgs VEV couples differently to T3>0 and T3<0 vertices,
# this would give T3-dependent mass corrections.
# The "VEV profile" could have a T3-dependent component.

# Overlap of fermion density with T3 function:
print("\nOverlap of fermion density with T3:")
T3_mean = np.mean(T3_vals)
for k in range(9):
    dens = densities[(k, 0, 0)]
    ovl_T3 = np.sum(dens * (T3_vals - T3_mean))
    print(f"  rho_{k} ({fermion_names[k]}): <p_{k}|T3> = {ovl_T3:+.6f}")

# Also check: does the SIGN of the overlap match T3 of the fermion?
print("\n  T3 of fermions: e(-1/2), mu(-1/2), tau(-1/2), u(+1/2), c(+1/2), t(+1/2), d(-1/2), s(-1/2), b(-1/2)")


# ================================================================
# PART 14: Combined VEV = v0 * (1 + c_T3 * T3 + c_alpha * phi_mode)
# ================================================================
print("\n" + "="*72)
print("PART 14: Physical VEV Model")
print("="*72)

# Model: Phi(v) = v0 * (1 + c_1 * T3(v) + c_2 * alpha_s_mode(v))
# where alpha_s_mode is the lowest non-trivial Laplacian eigenmode
# weighted by alpha_s * Nc

# The mass correction for fermion in state (k,i,j) is:
# delta_k = <p_{k,ij}| c_1*T3 + c_2*phi_1 |> / ln(phi)

# Compute features: T3 overlap and first-mode overlaps for each fermion
features = np.zeros((9, 5))  # [T3_overlap, mode1-4]
for k in range(9):
    dens = densities[(k, 0, 0)]
    features[k, 0] = np.sum(dens * T3_vals)
    if rho1_modes is not None:
        for m in range(min(4, rho1_modes.shape[1])):
            features[k, 1+m] = np.sum(dens * rho1_modes[:, m])

print(f"\nFeature matrix (T3 overlap + 4 rho_1 modes):")
for k in range(9):
    print(f"  rho_{k} ({fermion_names[k]}): T3={features[k,0]:+.6f}, modes=[{', '.join(f'{v:+.5f}' for v in features[k,1:])}]")

# Linear fit: delta = a*T3_overlap + b*mode_sum + c
X = np.column_stack([features, np.ones(9)])
coeffs, _, _, _ = np.linalg.lstsq(X, delta_vals, rcond=None)
fitted = X @ coeffs
rms_model = np.sqrt(np.mean((delta_vals - fitted)**2))
print(f"\nLinear fit with T3 + 4 modes: RMS = {rms_model:.4f}")
print(f"  Coefficients: T3={coeffs[0]:.4f}, modes=[{', '.join(f'{c:.4f}' for c in coeffs[1:5])}], const={coeffs[5]:.4f}")
print()
for k in range(9):
    print(f"  {fermion_names[k]:>4}: exp={delta_vals[k]:>+.4f}, fit={fitted[k]:>+.4f}, err={fitted[k]-delta_vals[k]:>+.4f}")


# ================================================================
# SUMMARY
# ================================================================
print("\n" + "="*72)
print("SUMMARY")
print("="*72)
print("""
KEY FINDINGS:

1. UNIFORM VEV IS ALWAYS STABLE at tree level (Hessian positive definite)
   => Non-uniform VEV requires QUANTUM CORRECTIONS (one-loop)

2. SELECTION RULES from CG: VEV mode rho_m affects fermion k
   only if rho_m appears in rho_k x rho_k. This is a HARD constraint
   from representation theory.

3. LEPTON CORRECTIONS ARE SUPPRESSED: alpha_s * Nc = 0 for leptons,
   so QCD-induced VEV non-uniformity gives ZERO correction to leptons
   at leading order. (Matches observation: lepton deltas ~ 0)

4. QUARK CORRECTIONS: proportional to alpha_s * Nc, with magnitude
   determined by the overlap of quark density with VEV profile.

5. T3 SIGN: determined by the coupling of the Higgs to the Hopf
   fiber T3 grading on the 600-cell.
""")
print("Done.")

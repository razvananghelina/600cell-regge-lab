"""
exp338e: Corrections using McKay NEIGHBOR information
=====================================================
Key insight: rho_3 and rho_6 have identical CG_self patterns,
but DIFFERENT neighbors on the McKay graph. So neighbor-based
corrections can break this degeneracy.

Uses the COMPLETE 2I character table for all tensor products.
"""
import numpy as np
from itertools import permutations

PHI = (1 + np.sqrt(5)) / 2
PHIp = -1/PHI  # Galois conjugate
ALPHA_S = 1 / (2 * PHI**3)

# ================================================================
# Complete 2I character table (9 irreps x 9 conjugacy classes)
# ================================================================
# Classes: C1(1), C2(1), C10a(12), C10b(12), C5a(12), C5b(12), C6(20), C3(20), C4(30)
class_sizes = np.array([1, 1, 12, 12, 12, 12, 20, 20, 30])
N = 120  # |2I|

# Character table: chartab[k][c] = chi_{rho_k}(class c)
chartab = np.array([
    [1,  1,    1,     1,     1,     1,    1,   1,   1],    # rho_0 (dim 1)
    [2, -2,    PHI,  PHIp,  1/PHI, -PHI,  1,  -1,   0],   # rho_1 (dim 2)
    [3,  3,    PHI,  PHIp, PHIp,   PHI,   0,   0,  -1],   # rho_2 (dim 3)
    [4, -4,    1,     1,    -1,    -1,   -1,   1,   0],    # rho_3 (dim 4)
    [5,  5,    0,     0,     0,     0,   -1,  -1,   1],    # rho_4 (dim 5)
    [6, -6,   -1,    -1,     1,     1,    0,   0,   0],    # rho_5 (dim 6)
    [4,  4,   -1,    -1,    -1,    -1,    1,   1,   0],    # rho_6 (dim 4)
    [2, -2,   PHIp,  PHI,  -PHI,  1/PHI, 1,  -1,   0],    # rho_7 (dim 2)
    [3,  3,   PHIp,  PHI,   PHI,  PHIp,  0,   0,  -1],    # rho_8 (dim 3)
])

dims = chartab[:, 0].astype(int)

print("="*72)
print("exp338e: McKay Neighbor Corrections")
print("="*72)

# Verify character table: orthogonality
print("\nVerifying character table orthogonality...")
gram = np.zeros((9, 9))
for i in range(9):
    for j in range(9):
        gram[i,j] = np.sum(class_sizes * chartab[i] * chartab[j]) / N
max_off = max(abs(gram[i,j] - (1 if i==j else 0)) for i in range(9) for j in range(9))
print(f"  Max orthogonality error: {max_off:.2e}")

# Verify: sum of d_k^2 = N
print(f"  Sum d_k^2 = {sum(d**2 for d in dims)} (should be {N})")

# ================================================================
# Compute ALL tensor product multiplicities
# ================================================================
def tensor_mult(k, j, m):
    """Multiplicity of rho_m in rho_k tensor rho_j."""
    return int(round(np.sum(class_sizes * chartab[k] * chartab[j] * chartab[m]) / N))

print("\nComputing all tensor product decompositions...")
# Full tensor product table: tp[k][j][m] = mult of rho_m in rho_k x rho_j
tp = np.zeros((9, 9, 9), dtype=int)
for k in range(9):
    for j in range(9):
        for m in range(9):
            tp[k, j, m] = tensor_mult(k, j, m)

# Verify McKay: rho_1 x rho_k should give neighbors of k
print("\nMcKay graph verification (rho_1 x rho_k):")
for k in range(9):
    neighbors = [m for m in range(9) if tp[1, k, m] > 0]
    print(f"  rho_1 x rho_{k} = " + " + ".join(f"rho_{m}" for m in neighbors))

# CG_self verification
print("\nCG_self verification (rho_k x rho_k):")
CG_self = np.zeros((9, 9), dtype=int)
for k in range(9):
    for m in range(9):
        CG_self[k, m] = tp[k, k, m]
    even = [CG_self[k, m] for m in [0, 2, 4, 6, 8]]
    print(f"  rho_{k}xrho_{k} even modes: {even}")

# ================================================================
# McKay graph structure
# ================================================================
mckay_adj = np.zeros((9, 9), dtype=int)
for k in range(9):
    for m in range(9):
        if tp[1, k, m] > 0:
            mckay_adj[k, m] = 1

# McKay graph distance
mckay_dist = np.full((9, 9), 99)
for i in range(9):
    mckay_dist[i, i] = 0
for i in range(9):
    for j in range(9):
        if mckay_adj[i, j]:
            mckay_dist[i, j] = 1
for k in range(9):
    for i in range(9):
        for j in range(9):
            if mckay_dist[i, k] + mckay_dist[k, j] < mckay_dist[i, j]:
                mckay_dist[i, j] = mckay_dist[i, k] + mckay_dist[k, j]

# ================================================================
# Cayley eigenvalues
# ================================================================
cayley_eigs = np.array([12 * (1 - chartab[k, 2] / dims[k]) for k in range(9)])
even_modes = [2, 4, 6, 8]
lam = np.array([cayley_eigs[m] for m in even_modes])
d_m = np.array([dims[m] for m in even_modes])

# ================================================================
# Fermion data
# ================================================================
fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
fermion_T3 = np.array([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5])
fermion_Nc = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1])
delta_exp = np.array([0.000, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278])

# McKay assignment (Solution 3, Cayley numbering):
# e->0, u->1, d->2, mu->3or4, s->4or3, c->5, tau->6, t->7, b->8
# Test all mu/s swaps
mckay_assigns = {
    'McKay_A': [0, 3, 6, 1, 5, 7, 2, 4, 8],  # mu->3, s->4
    'McKay_B': [0, 4, 6, 1, 5, 7, 2, 3, 8],  # mu->4, s->3
}

# Also the best from exp338d
mckay_assigns['McKay_best'] = [0, 4, 6, 1, 5, 8, 2, 3, 7]  # tau->6,t->8,b->7

# VEV optimal
assign_optimal = [0, 3, 5, 4, 2, 8, 7, 1, 6]

print(f"\n{'='*72}")
print("CROSS-TERM CG: rho_k x rho_j for McKay neighbors")
print("="*72)

# For each node k, compute the CG with its neighbors
print("\nNeighbor CG vectors (even modes only):")
print(f"{'Node':>6} {'Neighbors':>20} {'Self CG':>20} {'Avg Nbr CG':>20}")
for k in range(9):
    nbrs = [j for j in range(9) if mckay_adj[k, j]]
    self_cg = [CG_self[k, m] for m in even_modes]
    # Average neighbor CG
    nbr_cg = np.zeros(4)
    for j in nbrs:
        for idx, m in enumerate(even_modes):
            nbr_cg[idx] += CG_self[j, m]
    nbr_cg /= len(nbrs)
    print(f"rho_{k:>2} {str(nbrs):>20} {str(self_cg):>20} {str(np.round(nbr_cg,2).tolist()):>20}")

# ================================================================
# MECHANISM 1: Diffused CG (self + neighbor average)
# ================================================================
print(f"\n{'='*72}")
print("MECHANISM 1: Diffused CG = alpha*self + beta*avg_neighbor")
print("="*72)

def build_design_matrix(assignment, cg_matrix):
    """Build M[f, idx] = T3*Nc * cg_matrix[k_f, even_modes[idx]] * d_m."""
    M = np.zeros((9, 4))
    for f in range(9):
        k = assignment[f]
        for idx, m_idx in enumerate(even_modes):
            M[f, idx] = fermion_T3[f] * fermion_Nc[f] * cg_matrix[k, m_idx] * dims[m_idx]
    return M

# Build diffused CG matrices for different alpha/beta
def diffused_cg(alpha_val, beta_val):
    """CG_eff[k, m] = alpha*CG_self[k,m] + beta*avg_j~k(CG_self[j,m])"""
    cg_eff = np.zeros((9, 9))
    for k in range(9):
        nbrs = [j for j in range(9) if mckay_adj[k, j]]
        for m in range(9):
            nbr_avg = sum(CG_self[j, m] for j in nbrs) / len(nbrs)
            cg_eff[k, m] = alpha_val * CG_self[k, m] + beta_val * nbr_avg
    return cg_eff

# Check: does diffusion break rho_3/rho_6 degeneracy?
print("\nDiffused CG at alpha=1, beta=1 (even modes):")
cg_diff = diffused_cg(1.0, 1.0)
for k in range(9):
    vals = [cg_diff[k, m] for m in even_modes]
    print(f"  rho_{k}: {[round(v,2) for v in vals]}")

print(f"\n  rho_3 vs rho_6 at alpha=1, beta=1:")
print(f"  rho_3: {[round(cg_diff[3, m],2) for m in even_modes]}")
print(f"  rho_6: {[round(cg_diff[6, m],2) for m in even_modes]}")
degeneracy_broken = any(abs(cg_diff[3, m] - cg_diff[6, m]) > 0.01 for m in even_modes)
print(f"  Degeneracy broken: {degeneracy_broken}")

# Scan alpha, beta for each McKay assignment
for name, assign in mckay_assigns.items():
    print(f"\n--- {name}: {assign} ---")
    best = (999, 0, 0, None)
    for alpha_v in np.linspace(0, 2, 50):
        for beta_v in np.linspace(-2, 2, 50):
            cg_eff = diffused_cg(alpha_v, beta_v)
            M = build_design_matrix(assign, cg_eff)
            result = np.linalg.lstsq(M, delta_exp, rcond=None)
            c_opt = result[0]
            pred = M @ c_opt
            rms = np.sqrt(np.mean((delta_exp - pred)**2))
            if rms < best[0]:
                best = (rms, alpha_v, beta_v, c_opt)
    rms, alph, bet, c_opt = best
    cg_eff = diffused_cg(alph, bet)
    M = build_design_matrix(assign, cg_eff)
    pred = M @ c_opt
    print(f"  Best: alpha={alph:.3f}, beta={bet:.3f}, RMS={rms:.4f}")
    print(f"  c_m = {np.round(c_opt, 4)}")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}, err={pred[f]-delta_exp[f]:>+.4f}")

# ================================================================
# MECHANISM 2: Cross-CG (tensor product with neighbors)
# ================================================================
print(f"\n{'='*72}")
print("MECHANISM 2: Cross-CG tensor products rho_k x rho_j (j=neighbor)")
print("="*72)

# For fermion at node k, the correction involves CG of rho_k with its
# McKay neighbors. Define:
# CG_cross[k, m] = sum_{j~k} tp[k, j, m]  (total cross-CG with neighbors)
CG_cross = np.zeros((9, 9), dtype=int)
for k in range(9):
    for m in range(9):
        for j in range(9):
            if mckay_adj[k, j]:
                CG_cross[k, m] += tp[k, j, m]

print("\nCross-CG (sum over neighbors) even modes:")
print(f"{'Node':>6} {'cross(2,4,6,8)':>25}")
for k in range(9):
    vals = [CG_cross[k, m] for m in even_modes]
    print(f"rho_{k:>2} {str(vals):>25}")

# Check rho_3 vs rho_6
print(f"\nrho_3 cross: {[CG_cross[3, m] for m in even_modes]}")
print(f"rho_6 cross: {[CG_cross[6, m] for m in even_modes]}")
print(f"DIFFERENT: {any(CG_cross[3,m] != CG_cross[6,m] for m in even_modes)}")

# Use cross-CG as design matrix
for name, assign in mckay_assigns.items():
    M = np.zeros((9, 4))
    for f in range(9):
        k = assign[f]
        for idx, m_idx in enumerate(even_modes):
            M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_cross[k, m_idx] * dims[m_idx]
    result = np.linalg.lstsq(M, delta_exp, rcond=None)
    c_opt = result[0]
    pred = M @ c_opt
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    print(f"\n{name} with cross-CG: RMS={rms:.4f}")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# ================================================================
# MECHANISM 3: Mixed self + cross
# ================================================================
print(f"\n{'='*72}")
print("MECHANISM 3: alpha*self_CG + beta*cross_CG (8 columns total)")
print("="*72)

for name, assign in mckay_assigns.items():
    # Build 8-column design matrix: 4 from self, 4 from cross
    M = np.zeros((9, 8))
    for f in range(9):
        k = assign[f]
        for idx, m_idx in enumerate(even_modes):
            M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
            M[f, idx+4] = fermion_T3[f] * fermion_Nc[f] * CG_cross[k, m_idx] * dims[m_idx]
    result = np.linalg.lstsq(M, delta_exp, rcond=None)
    c_opt = result[0]
    pred = M @ c_opt
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    print(f"\n{name}: RMS={rms:.4f} (8 params)")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# ================================================================
# MECHANISM 4: McKay Laplacian eigenmodes
# ================================================================
print(f"\n{'='*72}")
print("MECHANISM 4: McKay graph Laplacian eigenmodes")
print("="*72)

# McKay 9x9 Laplacian
mckay_degree = np.sum(mckay_adj, axis=1)
L_mckay = np.diag(mckay_degree) - mckay_adj
evals_mck, evecs_mck = np.linalg.eigh(L_mckay.astype(float))

print("\nMcKay Laplacian eigenvalues:")
for i in range(9):
    print(f"  lambda_{i} = {evals_mck[i]:.4f}")

print(f"\nMcKay Laplacian eigenvectors (columns):")
print(f"{'node':>6}", end="")
for i in range(9):
    print(f"  v_{i:>2}  ", end="")
print()
for k in range(9):
    print(f"rho_{k:>2}", end="")
    for i in range(9):
        print(f" {evecs_mck[k,i]:>+6.3f}", end="")
    print()

# For each McKay assignment, use McKay eigenmodes as features
for name, assign in mckay_assigns.items():
    # Feature: f_i = T3*Nc * v_mode(k_f)
    # Use non-trivial McKay modes (skip mode 0 which is constant)
    n_modes = 8  # skip the zero mode
    M = np.zeros((9, n_modes))
    for f in range(9):
        k = assign[f]
        for i in range(n_modes):
            M[f, i] = fermion_T3[f] * fermion_Nc[f] * evecs_mck[k, i+1]
    result = np.linalg.lstsq(M, delta_exp, rcond=None)
    c_opt = result[0]
    pred = M @ c_opt
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    print(f"\n{name} with McKay eigenmodes: RMS={rms:.4f} ({n_modes} modes)")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# ================================================================
# MECHANISM 5: Heat kernel on McKay graph
# ================================================================
print(f"\n{'='*72}")
print("MECHANISM 5: Heat kernel diffusion from specific source nodes")
print("="*72)

# Heat kernel: K(k, source; t) = sum_i exp(-t*lambda_i) * v_i(k) * v_i(source)
# The correction for fermion at node k could be K(k, source; t)

def heat_kernel(source, t):
    """Heat kernel from source node at time t."""
    K = np.zeros(9)
    for k in range(9):
        for i in range(9):
            K[k] += np.exp(-t * evals_mck[i]) * evecs_mck[k, i] * evecs_mck[source, i]
    return K

# Try different source nodes and diffusion times
print("\nBest heat kernel source + time for each McKay assignment:")
for name, assign in mckay_assigns.items():
    best = (999, 0, 0)
    for source in range(9):
        for t in np.logspace(-2, 2, 200):
            K = heat_kernel(source, t)
            # Correction: delta_f = A * T3*Nc * K(k_f, source; t)
            pred = np.zeros(9)
            for f in range(9):
                pred[f] = fermion_T3[f] * fermion_Nc[f] * K[assign[f]]
            denom = np.dot(pred, pred)
            if denom > 1e-20:
                A = np.dot(delta_exp, pred) / denom
                pred_scaled = A * pred
                rms = np.sqrt(np.mean((delta_exp - pred_scaled)**2))
                if rms < best[0]:
                    best = (rms, source, t)
    rms, source, t = best
    K = heat_kernel(source, t)
    pred = np.zeros(9)
    for f in range(9):
        pred[f] = fermion_T3[f] * fermion_Nc[f] * K[assign[f]]
    A = np.dot(delta_exp, pred) / np.dot(pred, pred)
    pred = A * pred
    print(f"\n{name}: source=rho_{source}, t={t:.4f}, A={A:.4f}, RMS={rms:.4f}")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# ================================================================
# MECHANISM 6: Neighbor mass feedback
# ================================================================
print(f"\n{'='*72}")
print("MECHANISM 6: Neighbor mass / exponent feedback")
print("="*72)

# The correction at node k depends on the masses of its McKay neighbors:
# delta_k ~ T3 * Nc * f(sum of neighbor exponents / masses)
mass_exp = np.array([0, 11, 17, 3, 16, 26, 5, 11, 19])

for name, assign in mckay_assigns.items():
    # Build node -> fermion exponent map
    node_exp = np.zeros(9)
    for f in range(9):
        node_exp[assign[f]] = mass_exp[f]

    # Feature 1: sum of neighbor exponents
    feat1 = np.zeros(9)
    for k in range(9):
        nbrs = [j for j in range(9) if mckay_adj[k, j]]
        feat1[k] = sum(node_exp[j] for j in nbrs) / len(nbrs)

    # Feature 2: max neighbor exponent
    feat2 = np.zeros(9)
    for k in range(9):
        nbrs = [j for j in range(9) if mckay_adj[k, j]]
        feat2[k] = max(node_exp[j] for j in nbrs)

    # Feature 3: difference from neighbor average
    feat3 = np.zeros(9)
    for k in range(9):
        nbrs = [j for j in range(9) if mckay_adj[k, j]]
        feat3[k] = node_exp[k] - sum(node_exp[j] for j in nbrs) / len(nbrs)

    # Use 3 features
    M = np.zeros((9, 3))
    for f in range(9):
        k = assign[f]
        M[f, 0] = fermion_T3[f] * fermion_Nc[f] * feat1[k]
        M[f, 1] = fermion_T3[f] * fermion_Nc[f] * feat2[k]
        M[f, 2] = fermion_T3[f] * fermion_Nc[f] * feat3[k]
    result = np.linalg.lstsq(M, delta_exp, rcond=None)
    c_opt = result[0]
    pred = M @ c_opt
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    print(f"\n{name}: RMS={rms:.4f} (3 neighbor-mass features)")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# ================================================================
# MECHANISM 7: D_F spectral (per-node eigenvalue spectrum)
# ================================================================
print(f"\n{'='*72}")
print("MECHANISM 7: Cayley eigenvalue weighting by CG")
print("="*72)

# Instead of mode amplitudes, weight by Cayley eigenvalue functions
# CG_weighted[k] = sum_m CG_self[k,m] * f(lambda_m)
# This is the same as exp338c but now WITH neighbor corrections

# Define effective CG that includes neighbor info
for name, assign in mckay_assigns.items():
    # Scan: delta = T3*Nc * [a*sum_m CG_self*f(lam) + b*sum_m CG_cross*g(lam)]
    # Where f and g are simple functions
    best = (999, None)

    # Try: f = 1/lam, g = 1/lam (same function, different CG)
    for a in np.linspace(-2, 2, 100):
        for b in np.linspace(-2, 2, 100):
            pred = np.zeros(9)
            for f_idx in range(9):
                k = assign[f_idx]
                for idx, m in enumerate(even_modes):
                    val = a * CG_self[k, m] + b * CG_cross[k, m]
                    pred[f_idx] += fermion_T3[f_idx] * fermion_Nc[f_idx] * val * dims[m] / cayley_eigs[m]
                    # weight by 1/lambda
            rms = np.sqrt(np.mean((delta_exp - pred)**2))
            if rms < best[0]:
                best = (rms, (a, b))

    rms, (a, b) = best
    pred = np.zeros(9)
    for f_idx in range(9):
        k = assign[f_idx]
        for idx, m in enumerate(even_modes):
            val = a * CG_self[k, m] + b * CG_cross[k, m]
            pred[f_idx] += fermion_T3[f_idx] * fermion_Nc[f_idx] * val * dims[m] / cayley_eigs[m]
    print(f"\n{name}: a*self+b*cross, 1/lam weight")
    print(f"  a={a:.4f}, b={b:.4f}, RMS={rms:.4f}")
    for f_idx in range(9):
        print(f"  {fermion_names[f_idx]:>4}: pred={pred[f_idx]:>+.4f}, exp={delta_exp[f_idx]:>+.4f}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*72}")
print("SUMMARY: Best RMS per mechanism and assignment")
print("="*72)
print("\nNote: Leptons always get delta=0 (Nc=0), contributing 0.032 to RMS.")
print("Bare formula RMS = 0.247. VEV optimal (free assign) RMS = 0.033.")
print("Best achievable with McKay: see above.")

print("\n\nDone.")

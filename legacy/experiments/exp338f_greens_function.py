"""
exp338f: Green's function (Higgs propagator) on McKay graph
===========================================================
Key idea (from user): the Higgs symmetry breaking is SOURCED at a
specific node and PROPAGATES along the McKay graph via the Green's
function G(k, source) = [(L + mu^2)^{-1}]_{k,source}.

The correction for fermion at node k is:
  delta_f = T3 * Nc * G(k_f, source; mu)

This has 2 params: source node (discrete) and mu (continuous).
Or 1 continuous param if the source is fixed by physics.

Also try: weighted Green's function with dimension-dependent mass.
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
ALPHA_S = 1 / (2 * PHI**3)

# ================================================================
# McKay graph (affine E8, Cayley numbering)
# ================================================================
dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])

# Adjacency from McKay graph
mckay_edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
mckay_adj = np.zeros((9, 9), dtype=int)
for i, j in mckay_edges:
    mckay_adj[i, j] = 1
    mckay_adj[j, i] = 1

# Graph Laplacian
degree = np.sum(mckay_adj, axis=1)
L = np.diag(degree.astype(float)) - mckay_adj.astype(float)

# Weighted Laplacian (edges weighted by sqrt(d_i * d_j))
L_weighted = np.zeros((9, 9))
for i in range(9):
    for j in range(9):
        if mckay_adj[i, j]:
            w = np.sqrt(dims[i] * dims[j])
            L_weighted[i, j] = -w
            L_weighted[i, i] += w

# McKay distance
mckay_dist = np.full((9, 9), 99)
for i in range(9):
    mckay_dist[i, i] = 0
for i, j in mckay_edges:
    mckay_dist[i, j] = 1
    mckay_dist[j, i] = 1
for k in range(9):
    for i in range(9):
        for j in range(9):
            if mckay_dist[i, k] + mckay_dist[k, j] < mckay_dist[i, j]:
                mckay_dist[i, j] = mckay_dist[i, k] + mckay_dist[k, j]

# ================================================================
# Fermion data
# ================================================================
fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
fermion_T3 = np.array([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5])
fermion_Nc = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1])
delta_exp = np.array([0.000, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278])
mass_exp = np.array([0, 11, 17, 3, 16, 26, 5, 11, 19])

# McKay assignments (Cayley numbering)
# Main: e->0, u->1, d->2, {mu,s}->{3,4}, c->5
# Branch: {tau,t,b} on {6,7,8}
assigns = {
    'McKay_B': [0, 4, 6, 1, 5, 7, 2, 3, 8],  # mu->4, s->3
    'McKay_best': [0, 4, 6, 1, 5, 8, 2, 3, 7],  # best from exp338d
    'VEV_optimal': [0, 3, 5, 4, 2, 8, 7, 1, 6],
}

print("="*72)
print("exp338f: Green's Function on McKay Graph")
print("="*72)

# ================================================================
# PART 1: Basic Green's function G = (L + mu^2)^{-1}
# ================================================================
print(f"\n{'='*72}")
print("PART 1: G(k, source; mu) = [(L + mu^2*I)^{-1}]_{k,source}")
print("="*72)

def greens_function(L_mat, mu_sq):
    """Compute Green's function matrix G = (L + mu^2 I)^{-1}."""
    return np.linalg.inv(L_mat + mu_sq * np.eye(9))

# For each McKay assignment, for each source, scan mu
for name, assign in assigns.items():
    best = (999, 0, 0, 0)
    for source in range(9):
        for mu_sq in np.logspace(-3, 3, 500):
            G = greens_function(L, mu_sq)
            pred = np.zeros(9)
            for f in range(9):
                pred[f] = fermion_T3[f] * fermion_Nc[f] * G[assign[f], source]
            denom = np.dot(pred, pred)
            if denom > 1e-20:
                A = np.dot(delta_exp, pred) / denom
                pred_s = A * pred
                rms = np.sqrt(np.mean((delta_exp - pred_s)**2))
                if rms < best[0]:
                    best = (rms, source, mu_sq, A)

    rms, source, mu_sq, A = best
    G = greens_function(L, mu_sq)
    pred = np.zeros(9)
    for f in range(9):
        pred[f] = A * fermion_T3[f] * fermion_Nc[f] * G[assign[f], source]

    print(f"\n{name}: source=rho_{source}, mu^2={mu_sq:.4f}, A={A:.4f}, RMS={rms:.4f}")
    print(f"  G(k, {source}) values: {np.round(G[:, source], 4)}")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}, err={pred[f]-delta_exp[f]:>+.4f}")

# ================================================================
# PART 2: Two-source Green's function
# ================================================================
print(f"\n{'='*72}")
print("PART 2: Two sources G = A1*G(k,s1) + A2*G(k,s2)")
print("="*72)

for name, assign in assigns.items():
    best = (999, None)
    for mu_sq in np.logspace(-2, 2, 100):
        G = greens_function(L, mu_sq)
        # Build design matrix with ALL 9 possible sources
        M_sources = np.zeros((9, 9))
        for f in range(9):
            for s in range(9):
                M_sources[f, s] = fermion_T3[f] * fermion_Nc[f] * G[assign[f], s]
        # Lstsq with all 9 sources
        result = np.linalg.lstsq(M_sources, delta_exp, rcond=None)
        pred = M_sources @ result[0]
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, mu_sq, result[0], pred)

        # Also try pairs of sources (2 params)
        for s1 in range(9):
            for s2 in range(s1+1, 9):
                M2 = np.column_stack([M_sources[:, s1], M_sources[:, s2]])
                res2 = np.linalg.lstsq(M2, delta_exp, rcond=None)
                pred2 = M2 @ res2[0]
                rms2 = np.sqrt(np.mean((delta_exp - pred2)**2))
                if rms2 < best[0]:
                    best = (rms2, mu_sq, res2[0], pred2)

    rms = best[0]
    print(f"\n{name}: best 2-source RMS={rms:.4f}")

# ================================================================
# PART 3: Dimension-weighted Green's function
# ================================================================
print(f"\n{'='*72}")
print("PART 3: Dimension-weighted Laplacian")
print("="*72)

for name, assign in assigns.items():
    best = (999, 0, 0, 0)
    for source in range(9):
        for mu_sq in np.logspace(-3, 3, 500):
            G = greens_function(L_weighted, mu_sq)
            pred = np.zeros(9)
            for f in range(9):
                pred[f] = fermion_T3[f] * fermion_Nc[f] * G[assign[f], source]
            denom = np.dot(pred, pred)
            if denom > 1e-20:
                A = np.dot(delta_exp, pred) / denom
                pred_s = A * pred
                rms = np.sqrt(np.mean((delta_exp - pred_s)**2))
                if rms < best[0]:
                    best = (rms, source, mu_sq, A)

    rms, source, mu_sq, A = best
    G = greens_function(L_weighted, mu_sq)
    pred = np.zeros(9)
    for f in range(9):
        pred[f] = A * fermion_T3[f] * fermion_Nc[f] * G[assign[f], source]

    print(f"\n{name}: source=rho_{source}, mu^2={mu_sq:.4f}, A={A:.4f}, RMS={rms:.4f}")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# ================================================================
# PART 4: Green's function with node-dependent mass
# ================================================================
print(f"\n{'='*72}")
print("PART 4: Node-dependent mass: (L + diag(m_k^2))^{-1}")
print("="*72)

# Physical idea: the Higgs mass^2 at each node depends on the
# fermion mass at that node. Heavy fermions give larger loop corrections.
# m_k^2 ~ alpha_s * n_k (mass exponent at node k)

for name, assign in assigns.items():
    # Build node exponent map
    node_exp = np.zeros(9)
    for f in range(9):
        node_exp[assign[f]] = mass_exp[f]

    best = (999, 0, 0, 0, 0)
    for source in range(9):
        for mu0 in np.logspace(-2, 2, 100):
            for alpha in np.linspace(-0.5, 2, 50):
                mass_diag = mu0 * (1 + alpha * node_exp / 26.0)  # normalize by max exp
                mass_diag = np.maximum(mass_diag, 0.01)  # keep positive
                G = np.linalg.inv(L + np.diag(mass_diag))
                pred = np.zeros(9)
                for f in range(9):
                    pred[f] = fermion_T3[f] * fermion_Nc[f] * G[assign[f], source]
                denom = np.dot(pred, pred)
                if denom > 1e-20:
                    A = np.dot(delta_exp, pred) / denom
                    pred_s = A * pred
                    rms = np.sqrt(np.mean((delta_exp - pred_s)**2))
                    if rms < best[0]:
                        best = (rms, source, mu0, alpha, A)

    rms, source, mu0, alpha, A = best
    mass_diag = mu0 * (1 + alpha * node_exp / 26.0)
    mass_diag = np.maximum(mass_diag, 0.01)
    G = np.linalg.inv(L + np.diag(mass_diag))
    pred = np.zeros(9)
    for f in range(9):
        pred[f] = A * fermion_T3[f] * fermion_Nc[f] * G[assign[f], source]
    print(f"\n{name}: src=rho_{source}, mu0={mu0:.4f}, alpha={alpha:.4f}, A={A:.4f}, RMS={rms:.4f}")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# ================================================================
# PART 5: Resistance distance
# ================================================================
print(f"\n{'='*72}")
print("PART 5: Effective resistance distance on McKay graph")
print("="*72)

# The effective resistance R(i,j) between nodes i and j:
# R(i,j) = (L^+)_{ii} + (L^+)_{jj} - 2*(L^+)_{ij}
# where L^+ is the pseudoinverse of L

L_pinv = np.linalg.pinv(L)
R_eff = np.zeros((9, 9))
for i in range(9):
    for j in range(9):
        R_eff[i, j] = L_pinv[i, i] + L_pinv[j, j] - 2 * L_pinv[i, j]

print("\nEffective resistance from rho_0:")
for k in range(9):
    print(f"  R(0, {k}) = {R_eff[0, k]:.4f} (graph dist = {mckay_dist[0, k]})")

# Use R_eff as correction: delta ~ T3 * Nc * f(R(k, source))
for name, assign in assigns.items():
    best = (999, 0, 0)
    for source in range(9):
        r_vals = np.array([R_eff[assign[f], source] for f in range(9)])
        for p in [-2, -1, -0.5, 0.5, 1, 2]:
            pred_raw = np.zeros(9)
            for f in range(9):
                if R_eff[assign[f], source] > 1e-10:
                    pred_raw[f] = fermion_T3[f] * fermion_Nc[f] * R_eff[assign[f], source]**p
            denom = np.dot(pred_raw, pred_raw)
            if denom > 1e-20:
                A = np.dot(delta_exp, pred_raw) / denom
                pred_s = A * pred_raw
                rms = np.sqrt(np.mean((delta_exp - pred_s)**2))
                if rms < best[0]:
                    best = (rms, source, p)

    rms, source, p = best
    pred = np.zeros(9)
    for f in range(9):
        if R_eff[assign[f], source] > 1e-10:
            pred[f] = fermion_T3[f] * fermion_Nc[f] * R_eff[assign[f], source]**p
    A = np.dot(delta_exp, pred) / np.dot(pred, pred)
    pred = A * pred
    print(f"\n{name}: src=rho_{source}, power={p}, A={A:.4f}, RMS={rms:.4f}")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# ================================================================
# PART 6: Full propagator scan (all Laplacians, all sources)
# ================================================================
print(f"\n{'='*72}")
print("PART 6: Best single-source propagator (comprehensive)")
print("="*72)

# Try different Laplacians
laplacians = {
    'unweighted': L,
    'dim_weighted': L_weighted,
}

# Also: Cayley-eigenvalue weighted Laplacian
# L_cayley[k,j] = -cayley_eig_weight if adjacent, diagonal compensates
char_C10a = [1, PHI, PHI, 1, 0, -1, -1, -1/PHI, -1/PHI]
cayley_eigs = np.array([12 * (1 - char_C10a[k] / dims[k]) for k in range(9)])

L_cayley = np.zeros((9, 9))
for i, j in mckay_edges:
    w = np.sqrt(cayley_eigs[i] * cayley_eigs[j]) if cayley_eigs[i] > 0 and cayley_eigs[j] > 0 else 1.0
    L_cayley[i, j] = -w
    L_cayley[j, i] = -w
    L_cayley[i, i] += w
    L_cayley[j, j] += w
laplacians['cayley_weighted'] = L_cayley

for lap_name, L_mat in laplacians.items():
    print(f"\n--- Laplacian: {lap_name} ---")
    for name, assign in assigns.items():
        best = (999, 0, 0, 0)
        for source in range(9):
            for mu_sq in np.logspace(-4, 4, 1000):
                try:
                    G = np.linalg.inv(L_mat + mu_sq * np.eye(9))
                except:
                    continue
                pred = np.zeros(9)
                for f in range(9):
                    pred[f] = fermion_T3[f] * fermion_Nc[f] * G[assign[f], source]
                denom = np.dot(pred, pred)
                if denom > 1e-20:
                    A = np.dot(delta_exp, pred) / denom
                    pred_s = A * pred
                    rms = np.sqrt(np.mean((delta_exp - pred_s)**2))
                    if rms < best[0]:
                        best = (rms, source, mu_sq, A)
        rms, source, mu_sq, A = best
        print(f"  {name:>15}: src=rho_{source}, mu^2={mu_sq:.4f}, RMS={rms:.4f}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*72}")
print("SUMMARY TABLE")
print("="*72)
print(f"\n{'Method':>35} {'Params':>7} {'McKay RMS':>10} {'Optimal RMS':>12}")
print("-"*70)
print(f"{'Bare formula':>35} {'0':>7} {'0.247':>10} {'0.247':>12}")
print(f"{'T3*Nc*CG self (lstsq)':>35} {'4':>7} {'0.074':>10} {'0.033':>12}")
print(f"{'Diffused CG (best)':>35} {'4+2':>7} {'0.036':>10} {'-':>12}")
print(f"{'Green fn 1-source':>35} {'1+1':>7} {'see above':>10} {'-':>12}")
print(f"{'Green fn node-dep mass':>35} {'1+2':>7} {'see above':>10} {'-':>12}")
print(f"{'Heat kernel':>35} {'1+1':>7} {'0.105':>10} {'-':>12}")
print(f"\nLepton floor: sqrt((0.08^2+0.055^2)/9) = {np.sqrt((0.08**2+0.055**2)/9):.4f}")

print("\n\nDone.")

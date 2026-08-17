"""
exp433_product_selfenergy.py
============================
PRODUCT-SPACE SELF-ENERGY: Cayley(120) x McKay(9)

The key insight: delta_2 comes not from Cayley alone (Σ=0 for muon)
nor McKay alone (G=52/81 ≠ sqrt(2/5)), but from the PRODUCT.

The McKay eigenvalues SHIFT the Cayley propagator away from its zero.
This is the exact KK mechanism on discrete spaces.

COMPUTATION:
  Σ_product(u_mu) = Σ_{j≠0} |φ_j(u_mu)|² × G_Cayley(μ_j)

  where:
  - φ_j = McKay eigenvectors at muon node
  - μ_j = McKay Laplacian eigenvalues
  - G_Cayley(E) = (1/120) × Σ_k d_k² / (E + L_k) = Cayley resolvent

Author: Razvan-Constantin Anghelina
Date: March 2026
"""

import numpy as np
from numpy.linalg import eigh, pinv, eigvalsh

# =============================================================================
# CONSTANTS
# =============================================================================
a1 = 5
phi = (1 + np.sqrt(5)) / 2
alpha_em = 7.2973525693e-3
m_e = 0.51099895000e-3  # GeV
m_mu = 105.6583755e-3   # GeV

# Precise target
n_exact = np.log(m_mu / m_e) / np.log(phi)
n_bare = 11
delta_1_mu = phi**1.5 / 26.0
delta_total = n_exact - n_bare
delta_2_needed = delta_total - delta_1_mu
c_needed = delta_2_needed / (alpha_em * delta_1_mu)

print("=" * 80)
print("exp433: PRODUCT-SPACE SELF-ENERGY (Cayley x McKay)")
print("=" * 80)
print()
print("Target: c_needed = %.15f" % c_needed)
print("        sqrt(2/5) = %.15f" % np.sqrt(2.0/5))
print("        delta_2_needed = %.6e" % delta_2_needed)
print()

# =============================================================================
# McKAY GRAPH (Internal Space)
# =============================================================================
print("=" * 80)
print("STEP 1: McKay graph eigendecomposition")
print("=" * 80)
print()

dims_mckay = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3], dtype=float)
names = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
edges_mckay = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
N_mckay = 9
mu_idx = 4  # muon node
e_idx = 0   # electron node

A_mckay = np.zeros((N_mckay, N_mckay))
for i, j in edges_mckay:
    A_mckay[i, j] = 1
    A_mckay[j, i] = 1

D_mckay = np.diag(A_mckay.sum(axis=1))
L_mckay = D_mckay - A_mckay

# Eigendecomposition
mu_evals, mu_evecs = eigh(L_mckay)

print("McKay Laplacian eigenvalues:")
for j in range(N_mckay):
    print("  mu_%d = %.10f  |  phi_%d(mu) = %+.8f  |  |phi_%d(mu)|^2 = %.10f" % (
        j, mu_evals[j], j, mu_evecs[mu_idx, j], j, mu_evecs[mu_idx, j]**2))
print()

# Verify: sum |phi_j(mu)|^2 = 1
print("Check: sum |phi_j(mu)|^2 = %.10f (should be 1)" % sum(mu_evecs[mu_idx, :]**2))
print()

# =============================================================================
# CAYLEY GRAPH (Spacetime) - Analytical eigenvalues
# =============================================================================
print("=" * 80)
print("STEP 2: Cayley graph eigenvalues")
print("=" * 80)
print()

# chi(10a) values - VERIFIED in exp430
chi_10a = np.array([1, phi, phi, 1, 0, -1, -1, -1/phi, -1/phi])
L_cayley = 12.0 * (1 - chi_10a / dims_mckay)
mult_cayley = (dims_mckay**2).astype(int)  # multiplicities = d_k^2
N_cayley = 120

print("Cayley Laplacian eigenvalues (class 10a generators):")
for k in range(9):
    print("  L_%d = %12.8f  mult = %3d  (irrep %s, dim %d)" % (
        k, L_cayley[k], mult_cayley[k], names[k], int(dims_mckay[k])))
print("  Total multiplicity: %d (should be 120)" % sum(mult_cayley))
print()

# =============================================================================
# CAYLEY RESOLVENT: G_Cayley(v,v; E)
# =============================================================================
def G_cayley_resolvent(E):
    """
    Diagonal element of (E*I + L_Cayley)^{-1} for vertex-transitive graph.
    G(v,v; E) = (1/N) * sum_k d_k^2 / (E + L_k)

    For E > 0: all terms are finite (L_k >= 0).
    For E = 0: the L_0 = 0 term diverges, exclude it (pseudo-inverse).
    """
    result = 0.0
    for k in range(9):
        denom = E + L_cayley[k]
        if abs(denom) > 1e-15:
            result += mult_cayley[k] / denom
    return result / N_cayley

# Pure Cayley pseudo-inverse diagonal (E=0, excluding zero mode)
G_cayley_pinv = sum(mult_cayley[k] / L_cayley[k] for k in range(1, 9)) / N_cayley

print("Cayley pseudo-inverse G(v,v) = %.10f" % G_cayley_pinv)
print()

# =============================================================================
# PRODUCT-SPACE SELF-ENERGY
# =============================================================================
print("=" * 80)
print("STEP 3: Product-space self-energy at muon node")
print("=" * 80)
print()

# VERSION A: Full pseudo-inverse of product Laplacian at (v, u_mu)
# G_prod(u_mu) = sum_{(i,j) != (0,0)} |psi_i(v)|^2 * |phi_j(u_mu)|^2 / (L_i + mu_j)
# Since Cayley is vertex-transitive: |psi_i(v)|^2 = d_i^2/120

print("VERSION A: Full product pseudo-inverse diagonal at muon node")
print()

G_prod_mu = 0.0
for j in range(N_mckay):
    for k in range(9):
        if j == 0 and k == 0:
            continue  # skip zero mode
        denom = mu_evals[j] + L_cayley[k]
        if abs(denom) > 1e-15:
            G_prod_mu += (mult_cayley[k] / N_cayley) * mu_evecs[mu_idx, j]**2 / denom

# Same for electron (reference)
G_prod_e = 0.0
for j in range(N_mckay):
    for k in range(9):
        if j == 0 and k == 0:
            continue
        denom = mu_evals[j] + L_cayley[k]
        if abs(denom) > 1e-15:
            G_prod_e += (mult_cayley[k] / N_cayley) * mu_evecs[e_idx, j]**2 / denom

print("  G_prod(mu) = %.10f" % G_prod_mu)
print("  G_prod(e)  = %.10f" % G_prod_e)
print("  G_prod(mu) - G_prod(e) = %.10f" % (G_prod_mu - G_prod_e))
print()

# VERSION B: KK contribution only (j >= 1 terms)
# This is the part that comes from the internal space excitations
print("VERSION B: KK modes only (j >= 1)")
print()

Sigma_KK_mu = 0.0
Sigma_KK_e = 0.0
Sigma_KK_mu_terms = []
Sigma_KK_e_terms = []

for j in range(1, N_mckay):  # skip j=0 (zero mode)
    # Cayley resolvent at shifted mass mu_j
    G_shifted = G_cayley_resolvent(mu_evals[j])

    # Also include the zero-mode of Cayley: 1/(mu_j + 0) = 1/mu_j
    G_shifted_with_zero = G_shifted + (1.0 / N_cayley) / mu_evals[j]

    term_mu = mu_evecs[mu_idx, j]**2 * G_shifted_with_zero
    term_e = mu_evecs[e_idx, j]**2 * G_shifted_with_zero

    Sigma_KK_mu += term_mu
    Sigma_KK_e += term_e
    Sigma_KK_mu_terms.append((j, mu_evals[j], mu_evecs[mu_idx,j]**2, G_shifted_with_zero, term_mu))
    Sigma_KK_e_terms.append((j, mu_evals[j], mu_evecs[e_idx,j]**2, G_shifted_with_zero, term_e))

print("  Σ_KK(mu) = %.10f" % Sigma_KK_mu)
print("  Σ_KK(e)  = %.10f" % Sigma_KK_e)
print("  Σ_KK(mu) - Σ_KK(e) = %.10f" % (Sigma_KK_mu - Sigma_KK_e))
print()

print("  Per-mode breakdown for muon:")
print("  %3s  %12s  %12s  %12s  %12s" % ("j", "mu_j", "|phi_j|^2", "G_Cayley(mu_j)", "term"))
print("  " + "-" * 55)
for j, mu_j, phi2, G_sh, term in Sigma_KK_mu_terms:
    print("  %3d  %12.8f  %12.8f  %12.8f  %12.8f" % (j, mu_j, phi2, G_sh, term))
print("  Total: %.10f" % Sigma_KK_mu)
print()

# VERSION C: Only the Cayley zero-mode in the resolvent (1/mu_j term)
# This isolates the contribution where the spacetime propagator is the zero mode
print("VERSION C: Cayley zero-mode only (deepest IR)")
print()

Sigma_IR_mu = 0.0
Sigma_IR_e = 0.0
for j in range(1, N_mckay):
    term_mu = mu_evecs[mu_idx, j]**2 / (N_cayley * mu_evals[j])
    term_e = mu_evecs[e_idx, j]**2 / (N_cayley * mu_evals[j])
    Sigma_IR_mu += term_mu
    Sigma_IR_e += term_e

print("  Σ_IR(mu) = %.10f" % Sigma_IR_mu)
print("  Σ_IR(e)  = %.10f" % Sigma_IR_e)
print("  Difference = %.10f" % (Sigma_IR_mu - Sigma_IR_e))
print()

# =============================================================================
# COMPARISON WITH TARGET
# =============================================================================
print("=" * 80)
print("STEP 4: Comparison with target")
print("=" * 80)
print()

# We want: delta_2 = c * alpha * delta_1
# So: c = delta_2 / (alpha * delta_1)
# And: c_needed = 0.6323...

# Try various normalizations of the self-energy
print("Target: c_needed = %.10f  (sqrt(2/5) = %.10f)" % (c_needed, np.sqrt(2.0/5)))
print()

# List of candidate self-energy quantities to normalize
raw_quantities = {
    'G_prod(mu)': G_prod_mu,
    'G_prod(e)': G_prod_e,
    'G_prod(mu)-G_prod(e)': G_prod_mu - G_prod_e,
    'Sigma_KK(mu)': Sigma_KK_mu,
    'Sigma_KK(e)': Sigma_KK_e,
    'Sigma_KK(mu)-Sigma_KK(e)': Sigma_KK_mu - Sigma_KK_e,
    'Sigma_IR(mu)': Sigma_IR_mu,
    'Sigma_IR(e)': Sigma_IR_e,
    'Sigma_IR(mu)-Sigma_IR(e)': Sigma_IR_mu - Sigma_IR_e,
}

# Try normalizing by various factors
norm_factors = {
    '1': 1,
    'N_Cayley=120': N_cayley,
    'N_McKay=9': N_mckay,
    'N_2I*N_eig=1080': N_cayley * N_mckay,
    'a1=5': a1,
    'h(E8)=30': 30,
    'phi': phi,
    'phi^2': phi**2,
    'phi^3': phi**3,
    'pi': np.pi,
    '2*pi': 2*np.pi,
    '4*pi': 4*np.pi,
    'sqrt(4*pi)': np.sqrt(4*np.pi),
    'delta_1': delta_1_mu,
    '12 (degree)': 12,
    'sqrt(120)': np.sqrt(120),
    'sqrt(9)=3': 3,
    'sqrt(1080)': np.sqrt(1080),
}

print("Searching for normalization that gives c ≈ sqrt(2/5) = 0.6325...")
print()

hits = []
for q_name, q_val in raw_quantities.items():
    if abs(q_val) < 1e-15:
        continue
    for n_name, n_val in norm_factors.items():
        c_test = q_val * n_val
        if c_test > 0 and abs(c_test) < 100:
            err = abs(c_test - c_needed) / c_needed * 100
            if err < 5:  # within 5%
                hits.append((err, q_name, n_name, c_test))

        c_test2 = q_val / n_val
        if c_test2 > 0 and abs(c_test2) < 100:
            err2 = abs(c_test2 - c_needed) / c_needed * 100
            if err2 < 5:
                hits.append((err2, q_name, "1/" + n_name, c_test2))

# Also try sqrt of quantities
for q_name, q_val in raw_quantities.items():
    if q_val > 0:
        sq = np.sqrt(q_val)
        err = abs(sq - c_needed) / c_needed * 100
        if err < 10:
            hits.append((err, "sqrt(" + q_name + ")", "1", sq))
        for n_name, n_val in norm_factors.items():
            c_test = sq * n_val
            if 0 < c_test < 100:
                err = abs(c_test - c_needed) / c_needed * 100
                if err < 5:
                    hits.append((err, "sqrt(" + q_name + ")", n_name, c_test))
            c_test2 = sq / n_val
            if 0 < c_test2 < 100:
                err2 = abs(c_test2 - c_needed) / c_needed * 100
                if err2 < 5:
                    hits.append((err2, "sqrt(" + q_name + ")", "1/" + n_name, c_test2))

hits.sort()
print("%-6s %-35s %-20s %15s %10s" % ("Rank", "Quantity", "Normalization", "c_test", "Error(%)"))
print("-" * 90)
for rank, (err, q, n, c_val) in enumerate(hits[:25], 1):
    marker = " ***" if err < 0.1 else " **" if err < 1 else ""
    print("%-6d %-35s %-20s %15.10f %10.4f%s" % (rank, q[:35], n[:20], c_val, err, marker))

print()

# =============================================================================
# STEP 5: The spectral dimension hypothesis
# =============================================================================
print("=" * 80)
print("STEP 5: Spectral dimension of McKay graph")
print("=" * 80)
print()

# Heat kernel return probability P(t) = (1/N) Tr(exp(-tL))
# For short t: P(t) ~ t^{-d_s/2}
# d_s = -2 * d log P / d log t |_{t->0}

t_values = np.logspace(-3, 1, 100)
P_values = []
for t in t_values:
    P = sum(np.exp(-mu_evals[j] * t) for j in range(N_mckay)) / N_mckay
    P_values.append(P)

P_values = np.array(P_values)
# Compute local spectral dimension
d_s_values = []
for i in range(1, len(t_values)):
    d_s = -2 * (np.log(P_values[i]) - np.log(P_values[i-1])) / (np.log(t_values[i]) - np.log(t_values[i-1]))
    d_s_values.append(d_s)

print("Spectral dimension d_s at various scales:")
for t_idx in [5, 10, 20, 30, 50, 70, 90]:
    if t_idx < len(d_s_values):
        print("  t = %.4e:  d_s = %.6f" % (t_values[t_idx], d_s_values[t_idx-1]))

# At very short time, d_s should equal the topological dimension for regular graphs
# For a tree, it approaches 2 at intermediate times
max_ds = max(d_s_values)
print()
print("  Maximum d_s = %.6f (expected ~2 for tree at intermediate scale)" % max_ds)
print()

# =============================================================================
# STEP 6: Direct product-graph computation (build 1080x1080)
# =============================================================================
print("=" * 80)
print("STEP 6: Direct 1080-node product graph")
print("=" * 80)
print()

# Build the Cartesian product graph Cayley x McKay
# But we don't have the 120x120 Cayley adjacency matrix analytically.
# We can use the eigenvalue decomposition instead.

# On the product graph, the pseudo-inverse diagonal at (v, u_mu) is:
# G_prod = sum_{(i,j) != (0,0)} (d_i^2/120) * |phi_j(u_mu)|^2 / (L_i + mu_j)

# Let's compute this very carefully, mode by mode
print("Detailed product pseudo-inverse computation:")
print()

# First compute the j=0 contribution (McKay zero mode):
# For j=0: sum_{i>=1} (d_i^2/120) * |phi_0(u_mu)|^2 / (L_i + 0)
# = |phi_0(u_mu)|^2 * G_Cayley_pinv
phi0_mu_sq = mu_evecs[mu_idx, 0]**2
j0_contrib = phi0_mu_sq * G_cayley_pinv
print("j=0 contribution (pure Cayley): %.10f * %.10f = %.10f" % (
    phi0_mu_sq, G_cayley_pinv, j0_contrib))

# Now j>=1 contributions:
j_ge1_contrib = 0.0
print()
print("j>=1 contributions (KK modes):")
print("  %3s  %12s  %12s  %15s  %15s" % ("j", "mu_j", "|phi_j(mu)|^2", "G_Cayley(mu_j)", "contribution"))
print("  " + "-" * 65)

for j in range(1, N_mckay):
    phi_j_mu_sq = mu_evecs[mu_idx, j]**2
    # G_Cayley(mu_j) = (1/120) * sum_k d_k^2 / (mu_j + L_k)
    G_at_mu_j = sum(mult_cayley[k] / (mu_evals[j] + L_cayley[k]) for k in range(9)) / N_cayley
    contrib = phi_j_mu_sq * G_at_mu_j
    j_ge1_contrib += contrib
    print("  %3d  %12.8f  %12.8f  %15.10f  %15.10f" % (
        j, mu_evals[j], phi_j_mu_sq, G_at_mu_j, contrib))

print()
print("  Total j>=1: %.10f" % j_ge1_contrib)
print("  Total (j=0 + j>=1): %.10f" % (j0_contrib + j_ge1_contrib))
print("  Cross-check G_prod(mu): %.10f" % G_prod_mu)
print()

# =============================================================================
# STEP 7: What ratio gives sqrt(2/5)?
# =============================================================================
print("=" * 80)
print("STEP 7: What ratio of product-space quantities gives sqrt(2/5)?")
print("=" * 80)
print()

# The KK self-energy at muon node
Sigma = j_ge1_contrib  # KK contribution
Sigma_0 = j0_contrib   # zero-mode contribution

print("Σ_KK(mu) = %.10f (KK modes)" % Sigma)
print("Σ_0(mu)  = %.10f (zero mode)" % Sigma_0)
print("Σ_total(mu) = %.10f" % (Sigma + Sigma_0))
print()

# Same for other nodes
print("Product pseudo-inverse at all McKay nodes:")
print("  %5s  %12s  %12s  %12s" % ("Node", "Σ_0", "Σ_KK", "Σ_total"))
print("  " + "-" * 45)

node_sigma_kk = []
node_sigma_total = []
for u in range(N_mckay):
    s0 = mu_evecs[u, 0]**2 * G_cayley_pinv
    s_kk = 0.0
    for j in range(1, N_mckay):
        G_at_mu_j = sum(mult_cayley[k] / (mu_evals[j] + L_cayley[k]) for k in range(9)) / N_cayley
        s_kk += mu_evecs[u, j]**2 * G_at_mu_j
    node_sigma_kk.append(s_kk)
    node_sigma_total.append(s0 + s_kk)
    print("  %5s  %12.8f  %12.8f  %12.8f" % (names[u], s0, s_kk, s0 + s_kk))

print()

# Check ratios
print("Ratios Σ_KK(node) / Σ_KK(mu):")
for u in range(N_mckay):
    ratio = node_sigma_kk[u] / node_sigma_kk[mu_idx] if node_sigma_kk[mu_idx] != 0 else 0
    print("  %5s: %.10f" % (names[u], ratio))
print()

# Differences from electron
print("Differences Σ_KK(node) - Σ_KK(e):")
for u in range(N_mckay):
    diff = node_sigma_kk[u] - node_sigma_kk[e_idx]
    print("  %5s: %+.10f" % (names[u], diff))
print()

# The RELATIVE self-energy (difference from electron)
Delta_Sigma_mu = node_sigma_kk[mu_idx] - node_sigma_kk[e_idx]
print("Key quantity: ΔΣ_KK(mu) = Σ_KK(mu) - Σ_KK(e) = %.10f" % Delta_Sigma_mu)
print()

# Now try: c = ΔΣ_KK / (delta_1 * something)
# We want c * alpha * delta_1 = delta_2
# So: delta_2 = ΔΣ_KK * alpha * normalization?
# Or: c = ΔΣ_KK / something?

# Try various normalizations
print("Trying to match c_needed = %.10f:" % c_needed)
print()

quantities_to_test = {
    'Σ_KK(mu)': node_sigma_kk[mu_idx],
    'ΔΣ_KK(mu)': Delta_Sigma_mu,
    'Σ_total(mu)': node_sigma_total[mu_idx],
    'Σ_KK(mu)/Σ_0(mu)': node_sigma_kk[mu_idx] / (mu_evecs[mu_idx, 0]**2 * G_cayley_pinv),
}

for q_name, q_val in quantities_to_test.items():
    if abs(q_val) > 1e-15:
        # What normalization N makes q_val/N = c_needed?
        N_req = q_val / c_needed
        print("  %s = %.10f" % (q_name, q_val))
        print("    Needed normalization: %.10f" % N_req)
        print("    Close to: %.4f = %.4f*pi, %.4f*phi, %.4f*a1" % (
            N_req, N_req/np.pi, N_req/phi, N_req/a1))
        print()

# =============================================================================
# STEP 8: Deeper analysis - spectral dimension and loop factor
# =============================================================================
print("=" * 80)
print("STEP 8: The (4pi)^{d/2} factor")
print("=" * 80)
print()

# In KK theory, the one-loop factor from d compact dimensions is 1/(4*pi)^{d/2}
# For d_s = 2 (tree): 1/(4*pi) = 0.0796
# For d_s = 1: 1/sqrt(4*pi) = 0.2821

# What if c = N_something / sqrt(4*pi*a1)?
print("Loop factor tests:")
print("  1/(4*pi) = %.10f" % (1/(4*np.pi)))
print("  1/sqrt(4*pi) = %.10f" % (1/np.sqrt(4*np.pi)))
print("  1/sqrt(4*pi*a1) = %.10f" % (1/np.sqrt(4*np.pi*a1)))
print("  1/sqrt(2*pi*a1) = %.10f" % (1/np.sqrt(2*np.pi*a1)))
print("  sqrt(2/(4*pi)) = %.10f" % np.sqrt(2/(4*np.pi)))
print("  sqrt(a1/(4*pi*N_eig)) = %.10f" % np.sqrt(a1/(4*np.pi*N_mckay)))
print()

# Direct: sqrt(2/a1) = sqrt(2/5) from dimensional analysis?
print("Hypothesis: c = sqrt(d_spectral / a1)")
print("  d_spectral(tree) = 2 (random walk on tree)")
print("  c = sqrt(2/5) = %.10f" % np.sqrt(2.0/5))
print("  c_needed = %.10f" % c_needed)
print("  Agreement: %.4f%%" % (abs(np.sqrt(2.0/5) - c_needed)/c_needed*100))
print()

# =============================================================================
# STEP 9: Alternative - product trace formula
# =============================================================================
print("=" * 80)
print("STEP 9: Trace-type quantities")
print("=" * 80)
print()

# McKay zeta function
print("McKay spectral zeta function:")
mckay_nonzero = [mu_evals[j] for j in range(N_mckay) if mu_evals[j] > 1e-10]
for s in [0.5, 1.0, 1.5, 2.0]:
    zeta_s = sum(lam**(-s) for lam in mckay_nonzero)
    print("  zeta_McKay(%.1f) = %.10f" % (s, zeta_s))

print()

# The KEY quantity: what is sum_j |phi_j(mu)|^2 / mu_j?
# This is the McKay GREEN'S function at mu,mu
# = G_McKay(mu,mu) = 52/81
mckay_green_mu = sum(mu_evecs[mu_idx, j]**2 / mu_evals[j] for j in range(1, N_mckay))
print("McKay Green's function G(mu,mu) = sum |phi_j(mu)|^2/mu_j = %.10f" % mckay_green_mu)
print("  = 52/81 = %.10f" % (52.0/81))
print()

# What is sum_j |phi_j(mu)|^2 * G_Cayley(mu_j)?
# This is the KK self-energy Σ_KK(mu)
print("KK self-energy Σ_KK(mu) = sum |phi_j(mu)|^2 * G_Cayley(mu_j) = %.10f" % Sigma)
print()

# Ratio Σ_KK / G_McKay: how does the Cayley propagator modify the McKay Green's function?
ratio = Sigma / mckay_green_mu
print("Σ_KK(mu) / G_McKay(mu,mu) = %.10f" % ratio)
print("  If this is a simple function of framework params...")
print("  Ratio * 81/52 = %.10f" % (ratio * 81/52))
print("  Ratio / G_cayley_pinv = %.10f" % (ratio / G_cayley_pinv))
print()

# =============================================================================
# STEP 10: What if the self-energy IS the delta_2?
# =============================================================================
print("=" * 80)
print("STEP 10: Direct identification Σ_KK = delta_2/(alpha * something)")
print("=" * 80)
print()

# If delta_2 = alpha * Σ_KK * normalization
# Then: c = Σ_KK * normalization / delta_1
# And: normalization = c_needed * delta_1 / Σ_KK

norm_needed = c_needed * delta_1_mu / Sigma
print("If delta_2 = alpha * Σ_KK * N:")
print("  Σ_KK(mu) = %.10f" % Sigma)
print("  delta_1(mu) = %.10f" % delta_1_mu)
print("  c_needed = %.10f" % c_needed)
print("  => N = c_needed * delta_1 / Σ_KK = %.10f" % norm_needed)
print()

# Is this normalization a simple number?
print("  N = %.10f" % norm_needed)
print("  N*120 = %.6f" % (norm_needed*120))
print("  N*9 = %.6f" % (norm_needed*9))
print("  N*1080 = %.6f" % (norm_needed*1080))
print("  N*phi = %.6f" % (norm_needed*phi))
print("  N*pi = %.6f" % (norm_needed*np.pi))
print("  N*a1 = %.6f" % (norm_needed*a1))
print("  N/alpha = %.6f" % (norm_needed/alpha_em))
print("  1/N = %.6f" % (1/norm_needed))
print("  1/(N*alpha) = %.6f" % (1/(norm_needed*alpha_em)))
print()

# =============================================================================
# STEP 11: Ratios between self-energies at different nodes
# =============================================================================
print("=" * 80)
print("STEP 11: Self-energy ratios between nodes")
print("=" * 80)
print()

print("Σ_KK ratios (relative to muon):")
for u in range(N_mckay):
    if node_sigma_kk[mu_idx] > 0:
        r = node_sigma_kk[u] / node_sigma_kk[mu_idx]
        print("  Σ_KK(%s) / Σ_KK(mu) = %.10f" % (names[u], r))

print()
print("G_McKay ratios (for comparison):")
G_mckay_pinv = pinv(L_mckay)
for u in range(N_mckay):
    if G_mckay_pinv[mu_idx, mu_idx] > 0:
        r = G_mckay_pinv[u, u] / G_mckay_pinv[mu_idx, mu_idx]
        print("  G_McKay(%s) / G_McKay(mu) = %.10f  (G = %.6f)" % (
            names[u], r, G_mckay_pinv[u, u]))

# =============================================================================
# STEP 12: The c^2 = 2/5 test
# =============================================================================
print()
print("=" * 80)
print("STEP 12: Direct test of c^2 = 2/a1 hypothesis")
print("=" * 80)
print()

# If c^2 = 2/a1, this should come from some SQUARED self-energy
# or from a ratio involving spectral quantities
print("Squared self-energies:")
print("  Σ_KK(mu)^2 = %.10f" % Sigma**2)
print("  G_McKay(mu)^2 = %.10f" % mckay_green_mu**2)
print("  (Σ_KK/G_McKay)^2 = %.10f" % (Sigma/mckay_green_mu)**2)
print()

# Is there a quantity Q such that Q = 2/5?
print("Looking for Q = 2/a1 = 0.4:")
q_candidates = {
    'Σ_KK(mu)': Sigma,
    'G_McKay(mu,mu)': mckay_green_mu,
    'Σ_KK/G_McKay': Sigma / mckay_green_mu,
    'Σ_KK*a1': Sigma * a1,
    'Σ_KK*N_eig': Sigma * N_mckay,
    'Σ_KK*120': Sigma * 120,
}

for name, val in q_candidates.items():
    err = abs(val - 0.4) / 0.4 * 100
    if err < 50:
        print("  %s = %.10f (%.2f%% from 2/5)" % (name, val, err))

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("1. Product-space self-energy Σ_KK(mu) = %.10f" % Sigma)
print("2. McKay Green's function G(mu,mu) = 52/81 = %.10f" % (52/81))
print("3. Cayley resolvent G_Cayley(v,v) = %.10f" % G_cayley_pinv)
print("4. Ratio Σ_KK/G_McKay = %.10f" % (Sigma/mckay_green_mu))
print("5. c_needed = %.10f, sqrt(2/5) = %.10f" % (c_needed, np.sqrt(2.0/5)))
print()
print("Does Σ_KK naturally give sqrt(2/5)?")
print("  Direct: Σ_KK = %.6f (NOT sqrt(2/5))" % Sigma)
print("  Ratio to G_McKay: %.6f (NOT sqrt(2/5))" % (Sigma/mckay_green_mu))
print("  The product-space computation gives a DIFFERENT number")
print()
print("=" * 80)
print("END OF EXPERIMENT 433")
print("=" * 80)

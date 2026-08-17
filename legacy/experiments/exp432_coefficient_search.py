"""
exp432_coefficient_search.py
============================
SYSTEMATIC SEARCH for the true one-loop coefficient c in:
    delta_2 = c * alpha * delta_1

The experimentally needed c_needed ≈ 0.6323 must come from computable
framework quantities. sqrt(2/5) = 0.6325 is close but might not be
the true answer. We search BROADLY among:

1. Simple algebraic expressions in a1=5, phi, pi
2. McKay graph spectral quantities (Laplacian, Green's, heat kernel)
3. Cayley graph / group theory quantities
4. Z[phi] arithmetic
5. Trigonometric/hyperbolic combinations
6. Eigenvalue ratios and character values

TARGET: c_needed = 0.632254... (computed precisely below)

Author: Razvan-Constantin Anghelina
Date: March 2026
"""

import numpy as np
from numpy.linalg import pinv, eigvalsh, eigh
from fractions import Fraction

# =============================================================================
# CONSTANTS
# =============================================================================
a1 = 5
phi = (1 + np.sqrt(5)) / 2
alpha_em = 7.2973525693e-3  # CODATA 2018
m_e = 0.51099895000e-3  # GeV
m_mu = 105.6583755e-3   # GeV (PDG 2024)
m_mu_err = 0.0000023e-3 # GeV

# =============================================================================
# COMPUTE c_needed PRECISELY
# =============================================================================
n_exact = np.log(m_mu / m_e) / np.log(phi)
n_bare = 11  # 5*1 + 6*1 for mu

# Exact delta_1(mu) = phi^(3/2) / 26
delta_1_mu = phi**(1.5) / 26.0

delta_total = n_exact - n_bare
delta_2_needed = delta_total - delta_1_mu
c_needed = delta_2_needed / (alpha_em * delta_1_mu)

print("=" * 80)
print("exp432: SYSTEMATIC COEFFICIENT SEARCH")
print("=" * 80)
print()
print("PRECISE TARGET:")
print("  n_exact(mu) = %.15f" % n_exact)
print("  delta_1(mu) = phi^(3/2)/26 = %.15f" % delta_1_mu)
print("  delta_total = %.15f" % delta_total)
print("  delta_2_needed = %.15e" % delta_2_needed)
print("  c_needed = %.15f" % c_needed)
print()

# =============================================================================
# McKAY GRAPH SETUP
# =============================================================================
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
names = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
edges_mckay = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
N_mckay = 9

A_mckay = np.zeros((N_mckay, N_mckay))
for i, j in edges_mckay:
    A_mckay[i, j] = 1
    A_mckay[j, i] = 1

D_mckay = np.diag(A_mckay.sum(axis=1))
L_mckay = D_mckay - A_mckay
G_mckay = pinv(L_mckay)

# Eigendecomposition
evals_mckay, evecs_mckay = eigh(L_mckay)

# Graph distances
dist_mckay = np.zeros((N_mckay, N_mckay), dtype=int)
for start in range(N_mckay):
    visited = {start}
    queue = [(start, 0)]
    while queue:
        node, d = queue.pop(0)
        dist_mckay[start, node] = d
        for j in range(N_mckay):
            if A_mckay[node, j] > 0 and j not in visited:
                visited.add(j)
                queue.append((j, d + 1))

sigma_mckay = dist_mckay.sum(axis=1)
W_mckay = dist_mckay.sum() // 2  # Wiener index

# Normalized Laplacian
D_inv_sqrt = np.diag(1.0 / np.sqrt(A_mckay.sum(axis=1)))
L_norm = D_inv_sqrt @ L_mckay @ D_inv_sqrt
G_norm = pinv(L_norm)

# Random walk Laplacian
D_inv = np.diag(1.0 / A_mckay.sum(axis=1))
L_rw = D_inv @ L_mckay
G_rw = pinv(L_rw)

# Dimension-weighted Laplacian
D_dim = np.diag(np.array(dims, dtype=float))
D_dim_sqrt = np.diag(np.sqrt(dims))
D_dim_inv_sqrt = np.diag(1.0 / np.sqrt(dims))
L_dim = D_dim_sqrt @ L_mckay @ D_dim_sqrt  # weight by sqrt(d_i)*sqrt(d_j)
G_dim = pinv(L_dim)

# McKay adjacency spectrum
evals_adj_mckay = eigvalsh(A_mckay)

# =============================================================================
# HEAT KERNEL on McKay graph
# =============================================================================
def heat_kernel_diag(L_mat, t, node):
    """K(t, node, node) = sum_k |v_k(node)|^2 * exp(-lambda_k * t)"""
    evals, evecs = eigh(L_mat)
    return sum(evecs[node, k]**2 * np.exp(-evals[k] * t) for k in range(len(evals)))

# =============================================================================
# CAYLEY GRAPH (600-cell) eigenvalues
# =============================================================================
# chi(10a) values (VERIFIED in exp430)
chi_10a = [1, phi, phi, 1, 0, -1, -1, -1/phi, -1/phi]
L_cayley = [12 * (1 - chi_10a[k] / dims[k]) for k in range(9)]

# chi(5a) values (Galois conjugate)
chi_5a = [1, -1/phi, -1/phi, 1, 0, -1, -1, phi, phi]
L_cayley_5a = [12 * (1 - chi_5a[k] / dims[k]) for k in range(9)]

# Galois anti-symmetric
G_galois = [(L_cayley[k] - L_cayley_5a[k]) / 2 for k in range(9)]

# =============================================================================
# CANDIDATE EXPRESSIONS
# =============================================================================
candidates = []

def add(name, value, category="misc"):
    candidates.append((name, value, category))

# --- CATEGORY 1: Simple algebraic ---
add("sqrt(2/a1) = sqrt(2/5)", np.sqrt(2.0/5), "algebraic")
add("sqrt(2)/phi", np.sqrt(2)/phi, "algebraic")
add("2/pi", 2/np.pi, "algebraic")
add("(a1-1)/(2*pi)", (a1-1)/(2*np.pi), "algebraic")
add("1/phi", 1/phi, "algebraic")
add("phi/phi^2 = 1/phi", phi/phi**2, "algebraic")
add("phi^2 - 2", phi**2 - 2, "algebraic")
add("2 - phi", 2 - phi, "algebraic")
add("phi - 1 = 1/phi", phi - 1, "algebraic")
add("(phi+1)/(2*phi+1)", (phi+1)/(2*phi+1), "algebraic")
add("phi/(1+phi) = 1/phi", phi/(1+phi), "algebraic")
add("sqrt(phi) - 1", np.sqrt(phi) - 1, "algebraic")
add("a1/(2*a1-1) = 5/9", a1/(2*a1-1), "algebraic")
add("(a1-1)/a1^(3/2) = 4/sqrt(125)", (a1-1)/a1**1.5, "algebraic")
add("sqrt((a1-1)/a1^2) = 2/5", np.sqrt((a1-1)/a1**2), "algebraic")
add("1/phi^(2/3)", phi**(-2.0/3), "algebraic")
add("phi^(-0.95)", phi**(-0.95), "algebraic")
add("phi^(-0.953)", phi**(-0.953), "algebraic")

# Rational approximations near 0.6323
add("12/19", 12.0/19, "rational")
add("19/30 = 19/h(E8)", 19.0/30, "rational")
add("31/49", 31.0/49, "rational")
add("43/68", 43.0/68, "rational")
add("55/87", 55.0/87, "rational")
add("74/117", 74.0/117, "rational")

# Golden ratio related
add("(phi^3 - 3)/phi^2", (phi**3 - 3)/phi**2, "golden")
add("phi^2 / (phi^2 + phi + 1)", phi**2/(phi**2 + phi + 1), "golden")
add("(2*phi - 1)/(phi^2 + 1)", (2*phi-1)/(phi**2+1), "golden")
add("phi/(2 + phi)", phi/(2 + phi), "golden")
add("phi^2/(a1 + phi)", phi**2/(a1 + phi), "golden")
add("(phi-1)*(phi+1)/phi^2", (phi-1)*(phi+1)/phi**2, "golden")
add("phi^(3/2)/pi", phi**1.5/np.pi, "golden")
add("2*phi/(phi+a1)", 2*phi/(phi+a1), "golden")

# Fibonacci related
add("F(8)/F(10) = 21/55", 21.0/55, "fibonacci")
add("L(5)/L(7) = 11/29", 11.0/29, "fibonacci")

# --- CATEGORY 2: McKay graph spectral ---
mu_idx = 4

# Green's function values
add("G_McKay(mu,mu) = 52/81", G_mckay[mu_idx, mu_idx], "mckay")
add("G_McKay(e,e) = 196/81", G_mckay[0, 0], "mckay")
add("G_McKay(mu,mu)/G_McKay(e,e)*N/a1",
    G_mckay[mu_idx,mu_idx]/G_mckay[0,0]*N_mckay/a1, "mckay")

# Normalized Laplacian Green's function
add("G_norm(mu,mu)", G_norm[mu_idx, mu_idx], "mckay_norm")

# Random walk Green's function
add("G_rw(mu,mu)", G_rw[mu_idx, mu_idx], "mckay_rw")

# Dimension-weighted Green's function
add("G_dim(mu,mu)", G_dim[mu_idx, mu_idx], "mckay_dim")

# Off-diagonal Green's function
add("G_McKay(e,mu)", G_mckay[0, mu_idx], "mckay")
add("|G_McKay(e,mu)|", abs(G_mckay[0, mu_idx]), "mckay")

# Ratios of Green's function values
for j in range(N_mckay):
    if j != mu_idx and abs(G_mckay[mu_idx, mu_idx]) > 0:
        add("G(mu,mu)/G(%s,%s)" % (names[j],names[j]),
            G_mckay[mu_idx,mu_idx]/G_mckay[j,j] if G_mckay[j,j] != 0 else 999, "mckay_ratio")

# sigma / combinations
add("sigma_mu/(N*a1) = 18/45", sigma_mckay[mu_idx]/(N_mckay*a1), "mckay")
add("sigma_mu/W = 18/110", sigma_mckay[mu_idx]/W_mckay, "mckay")
add("W/(N^2*a1-1) = 110/365", W_mckay/(N_mckay**2*(a1-1)), "mckay")
add("(sigma_mu - N)/W = 9/110", (sigma_mckay[mu_idx]-N_mckay)/W_mckay, "mckay")
add("sqrt(sigma_mu/(2*W)) = sqrt(18/220)", np.sqrt(sigma_mckay[mu_idx]/(2*W_mckay)), "mckay")
add("sigma_mu^2/(N*W)", sigma_mckay[mu_idx]**2/(N_mckay*W_mckay), "mckay")

# McKay adjacency eigenvalues
evals_adj_sorted = sorted(evals_adj_mckay)
for k in range(N_mckay):
    add("McKay_adj_eval[%d] = %.4f" % (k, evals_adj_sorted[k]), evals_adj_sorted[k], "mckay_eval")

# Eigenvalue combinations
add("(eval_max - eval_min)/(2*eval_max)",
    (evals_adj_sorted[-1] - evals_adj_sorted[0])/(2*evals_adj_sorted[-1]), "mckay_eval")
add("sum(|evals_adj|)/N", sum(abs(e) for e in evals_adj_sorted)/N_mckay, "mckay_eval")

# Heat kernel at specific times
for t_val in [0.1, 0.2, 0.5, 1.0, alpha_em, 1/phi, phi, 2.0, np.pi, a1]:
    val = heat_kernel_diag(L_mckay, t_val, mu_idx)
    add("K_McKay(t=%.4f, mu, mu)" % t_val, val, "mckay_heat")

# Resistance distance combinations
R_e_mu = dist_mckay[0, mu_idx]  # = 4 on tree
add("R(e,mu)/phi^2 = 4/phi^2", R_e_mu/phi**2, "resistance")
add("R(e,mu)/N = 4/9", R_e_mu/N_mckay, "resistance")
add("R(e,mu)/(2*pi)", R_e_mu/(2*np.pi), "resistance")
add("sqrt(R(e,mu))/pi = 2/pi", np.sqrt(R_e_mu)/np.pi, "resistance")

# --- CATEGORY 3: Cayley graph / group theory ---
N_2I = 120
h_E8 = 30
rank_E8 = 8

add("h(E8)/N_2I * phi = 30*phi/120", h_E8*phi/N_2I, "cayley")
add("rank(E8)/N_2I * phi^3", rank_E8*phi**3/N_2I, "cayley")
add("a1/rank(E8) = 5/8", a1/rank_E8, "cayley")
add("sqrt(a1/N_2I) = sqrt(1/24)", np.sqrt(a1/N_2I), "cayley")
add("h(E8)/(N_2I/phi^2)", h_E8/(N_2I/phi**2), "cayley")

# Cayley eigenvalue combinations
L_mu_cayley = L_cayley[mu_idx]  # = 12
add("L_mu/N_2I * a1 = 12*5/120", L_mu_cayley*a1/N_2I, "cayley")
add("1/sqrt(L_mu/a1) = 1/sqrt(12/5)", 1/np.sqrt(L_mu_cayley/a1), "cayley_eval")

# Character value combinations
add("chi_mu(10a)/d_mu = 0/5 = 0", chi_10a[mu_idx]/dims[mu_idx], "character")
add("chi_d(10a)/d_d = phi/3", chi_10a[2]/dims[2], "character")

# Spectral zeta function of McKay Laplacian
mckay_Levals = sorted(eigvalsh(L_mckay))
mckay_nonzero = [e for e in mckay_Levals if e > 1e-10]
add("zeta_McKay(1) = sum(1/lambda)", sum(1/e for e in mckay_nonzero), "mckay_zeta")
add("zeta_McKay(1)/N = sum(1/lambda)/9", sum(1/e for e in mckay_nonzero)/N_mckay, "mckay_zeta")
add("1/zeta_McKay(1)", 1/sum(1/e for e in mckay_nonzero), "mckay_zeta")
add("zeta_McKay(2) = sum(1/lambda^2)", sum(1/e**2 for e in mckay_nonzero), "mckay_zeta")
add("sqrt(1/zeta_McKay(2))", 1/np.sqrt(sum(1/e**2 for e in mckay_nonzero)), "mckay_zeta")

# Kirchhoff's theorem: number of spanning trees
det_reduced = np.prod(mckay_nonzero) / N_mckay
add("det'(L_McKay)/N^2 = tau/N^2", det_reduced/N_mckay**2, "mckay")
add("tau^(1/8)/N", det_reduced**(1.0/8)/N_mckay, "mckay")

# --- CATEGORY 4: Z[phi] arithmetic ---
# z_mu = 1 + phi, N(z_mu) = (1+phi)(1+phi') = (1+phi)(2-phi) = 2+phi-phi^2 = 2+phi-(phi+1) = 1
# Wait, N(z) = z*z' = (1+phi)(1+(1-phi)) = (1+phi)(2-phi)
# = 2 - phi + 2*phi - phi^2 = 2 + phi - (phi+1) = 1
z_mu = 1 + phi
z_mu_conj = 1 + (1-phi)  # = 2 - phi = 1/phi^2
N_z_mu = z_mu * z_mu_conj  # should be 1

add("N(z_mu) = (1+phi)(2-phi)", N_z_mu, "zphi")
add("sqrt(|N(z_mu)|/a1)", np.sqrt(abs(N_z_mu)/a1), "zphi")

# z values for all fermions
z_vals = {
    'e': (0, 0),     # z = 0
    'u': (3, -2),    # z = 3 - 2*phi
    'd': (1, 0),     # z = 1
    's': (1, 1),     # z = 1 + phi
    'mu': (1, 1),    # z = 1 + phi
    'c': (2, 1),     # z = 2 + phi
    'tau': (1, 2),   # z = 1 + 2*phi
    't': (4, 1),     # z = 4 + phi
    'b': (-1, 4),    # z = -1 + 4*phi
}

# Z[phi] norms
zphi_norms = {}
for name, (a, b) in z_vals.items():
    z = a + b*phi
    zp = a + b*(1-phi)  # Galois conjugate
    N_z = z * zp
    zphi_norms[name] = N_z

# Norm combinations
add("|N(z_tau)|/|N(z_c)|", abs(zphi_norms['tau'])/abs(zphi_norms['c']), "zphi")
add("sqrt(|N(z_mu)|)", np.sqrt(abs(zphi_norms['mu'])), "zphi")

# --- CATEGORY 5: Trigonometric ---
add("sin(pi/a1) = sin(pi/5)", np.sin(np.pi/a1), "trig")
add("cos(pi/a1) = cos(pi/5)", np.cos(np.pi/a1), "trig")
add("sin(2*pi/a1) = sin(2pi/5)", np.sin(2*np.pi/a1), "trig")
add("cos(2*pi/a1) = cos(2pi/5)", np.cos(2*np.pi/a1), "trig")
add("sin(pi/a1)/sqrt(phi)", np.sin(np.pi/a1)/np.sqrt(phi), "trig")
add("cos(pi/a1)/phi", np.cos(np.pi/a1)/phi, "trig")
add("sin(pi/(2*a1))", np.sin(np.pi/(2*a1)), "trig")
add("cos(pi/(a1-1)) = cos(pi/4)", np.cos(np.pi/(a1-1)), "trig")
add("tan(pi/a1)", np.tan(np.pi/a1), "trig")
add("sin(pi/a1)*cos(pi/a1) = sin(2pi/a1)/2", np.sin(np.pi/a1)*np.cos(np.pi/a1), "trig")

# --- CATEGORY 6: Composite / dimensional ---
add("sqrt(2*sin^2(tW)) = sqrt(6/13)", np.sqrt(2*6.0/26), "composite")
add("sin^2(tW)*phi = (6/26)*phi", (6.0/26)*phi, "composite")
add("sqrt(alpha_s) = 1/sqrt(2*phi^3)", np.sqrt(1/(2*phi**3)), "composite")
add("alpha_s*phi^2", 1/(2*phi**3)*phi**2, "composite")
add("C*phi^2 = 2*phi^2/13", (2.0/13)*phi**2, "composite")
add("C*phi^(3/2)", (2.0/13)*phi**1.5, "composite")
add("sqrt(C) = sqrt(2/13)", np.sqrt(2.0/13), "composite")
add("sqrt(C)*phi^(1/4)", np.sqrt(2.0/13)*phi**0.25, "composite")
add("C*a1 = 10/13", (2.0/13)*a1, "composite")
add("c_ell*a1 = phi^3*a1/26", phi**3*a1/26, "composite")
add("sqrt(c_ell*a1)", np.sqrt(phi**3*a1/26), "composite")
add("c_ell*phi^(3/2) = phi^3/(26/phi^1.5)", phi**3*phi**1.5/26, "composite")

# (1 + pi*alpha)/phi and variants
add("(1 + pi*alpha)/phi", (1 + np.pi*alpha_em)/phi, "composite_alpha")
add("(1 + pi*alpha/2)/phi", (1 + np.pi*alpha_em/2)/phi, "composite_alpha")
add("(1 + a1*alpha)/phi", (1 + a1*alpha_em)/phi, "composite_alpha")
add("(1 + 2*pi*alpha)/phi", (1 + 2*np.pi*alpha_em)/phi, "composite_alpha")
add("1/phi + alpha/(2*phi)", 1/phi + alpha_em/(2*phi), "composite_alpha")
add("1/phi + alpha", 1/phi + alpha_em, "composite_alpha")
add("1/phi * (1 + alpha*phi^2)", (1 + alpha_em*phi**2)/phi, "composite_alpha")

# --- CATEGORY 7: Number-theoretic ---
add("(a1^2-1)/(2*a1^2) = 24/50", (a1**2-1)/(2*a1**2), "number")
add("(a1^2+1)/(2*a1^2+a1) = 26/55", (a1**2+1)/(2*a1**2+a1), "number")
add("(2*a1-1)/(2*a1+a1) = 9/15", (2*a1-1)/(3*a1), "number")
add("a1/(2*a1+phi) = 5/(10+phi)", a1/(2*a1+phi), "number")
add("sqrt(a1)/(2+sqrt(a1))", np.sqrt(a1)/(2+np.sqrt(a1)), "number")

# E8 related
add("rank(E8)/h(E8)*phi^2 = 8*phi^2/30", rank_E8*phi**2/h_E8, "E8")
add("(rank(E8)-a1)/(a1-2) = 3/3 = 1", (rank_E8-a1)/(a1-2), "E8")
add("sqrt(rank(E8)/N_2I) = sqrt(1/15)", np.sqrt(rank_E8/N_2I), "E8")
add("N_eig/(N_eig + a1) = 9/14", N_mckay/(N_mckay + a1), "E8")
add("N_eig*phi/h(E8) = 9*phi/30", N_mckay*phi/h_E8, "E8")

# --- CATEGORY 8: Sophisticated McKay combinations ---
# Spectral gap
lambda_1_mckay = min(mckay_nonzero)
add("1/lambda_1(McKay)", 1/lambda_1_mckay, "mckay_spectral")
add("lambda_1(McKay)/phi", lambda_1_mckay/phi, "mckay_spectral")
add("sqrt(lambda_1(McKay))", np.sqrt(lambda_1_mckay), "mckay_spectral")

# Cheeger constant approximation
deg_max = int(A_mckay.sum(axis=1).max())
deg_min = int(A_mckay.sum(axis=1).min())
add("deg_min/(deg_min+deg_max) = %d/%d" % (deg_min, deg_min+deg_max),
    deg_min/(deg_min+deg_max), "mckay_spectral")

# Effective resistance between specific nodes
# On tree: R_eff(i,j) = d(i,j)
add("sqrt(R(e,mu)/N_eig) = sqrt(4/9) = 2/3", np.sqrt(dist_mckay[0,4]/N_mckay), "resistance")
add("2/3", 2.0/3, "resistance")
add("R(e,mu)*phi/(N_eig*phi) = 4/9", dist_mckay[0,4]*phi/(N_mckay*phi), "resistance")

# Kirchhoff index
K_index = sum(1/e for e in mckay_nonzero)
add("N/K_index = 9/sum(1/lambda)", N_mckay/K_index, "mckay_spectral")
add("K_index/N^2 = sum(1/lambda)/81", K_index/N_mckay**2, "mckay_spectral")

# --- CATEGORY 9: McKay graph + Cayley cross-terms ---
add("G_mckay(mu,mu) * L_cayley(mu)/N_2I", G_mckay[mu_idx,mu_idx]*L_cayley[mu_idx]/N_2I, "cross")
add("sqrt(G_mckay(mu,mu) * L_cayley(mu)/N_2I)",
    np.sqrt(G_mckay[mu_idx,mu_idx]*L_cayley[mu_idx]/N_2I), "cross")
add("sigma_mu * a1 / N^3", sigma_mckay[mu_idx]*a1/N_mckay**3, "cross")

# --- CATEGORY 10: Morrey-Sobolev related ---
d_ST = 4
alpha_morrey = 3.0/4  # = (d_ST-1)/d_ST
add("alpha_morrey / sqrt(phi) = (3/4)/sqrt(phi)", alpha_morrey/np.sqrt(phi), "morrey")
add("alpha_morrey * sqrt(phi/a1)", alpha_morrey*np.sqrt(phi/a1), "morrey")
add("(d_ST-1)*phi / (d_ST*phi+1)", (d_ST-1)*phi/(d_ST*phi+1), "morrey")
add("alpha_morrey / phi^(1/4)", alpha_morrey/phi**0.25, "morrey")
add("(d_ST-1)/(d_ST+1/phi) = 3/(4+1/phi)", (d_ST-1)/(d_ST+1/phi), "morrey")
add("alpha_morrey * (1 - 1/(a1*d_ST))", alpha_morrey*(1 - 1/(a1*d_ST)), "morrey")

# --- CATEGORY 11: Interesting near-misses ---
# These involve 2-parameter combinations but could reveal structure
add("(phi^2 - 1)/phi^2 = 1/phi^2 = same as 2-phi", 1/phi**2, "structural")
add("1 - alpha_morrey*1/phi^2", 1 - alpha_morrey/phi**2, "structural")
add("phi*(1-1/phi^3)", phi*(1-1/phi**3), "structural")
add("(a1-1)/((a1-1)*phi + a1-3) = 4/(4*phi+2)", (a1-1)/((a1-1)*phi + a1-3), "structural")
add("phi^2/(phi^2+phi+1/phi) = phi^2/(phi^2+phi+1/phi)", phi**2/(phi**2+phi+1/phi), "structural")

# =============================================================================
# SORT AND DISPLAY RESULTS
# =============================================================================
print("=" * 80)
print("ALL CANDIDATES (sorted by |c - c_needed|)")
print("=" * 80)
print()

# Compute relative errors
results = []
for name, value, cat in candidates:
    if isinstance(value, (int, float, np.floating)) and np.isfinite(value) and value > 0:
        rel_err = abs(value - c_needed) / c_needed * 100
        results.append((rel_err, name, value, cat))

results.sort()

# Show top 40
print("%-6s %-55s %15s %12s  %s" % ("Rank", "Expression", "Value", "Error(%)", "Category"))
print("-" * 100)
for rank, (err, name, value, cat) in enumerate(results[:40], 1):
    marker = ""
    if err < 0.01:
        marker = " <-- EXCELLENT"
    elif err < 0.1:
        marker = " <-- GOOD"
    elif err < 0.5:
        marker = " <-- OK"
    print("%-6d %-55s %15.10f %12.6f  %-15s%s" % (rank, name[:55], value, err, cat, marker))

print()

# =============================================================================
# DETAILED ANALYSIS OF TOP CANDIDATES
# =============================================================================
print("=" * 80)
print("DETAILED ANALYSIS OF TOP CANDIDATES")
print("=" * 80)
print()

# Filter top 10
top_N = 15
print("c_needed = %.15f" % c_needed)
print()

for rank, (err, name, value, cat) in enumerate(results[:top_N], 1):
    # Compute mass prediction
    delta_2_test = value * alpha_em * delta_1_mu
    n_total = n_bare + delta_1_mu + delta_2_test
    m_pred = m_e * phi**n_total * 1e3  # back to GeV
    m_diff = (m_pred - m_mu*1e3) * 1e6  # in eV
    pull = m_diff / (m_mu_err * 1e6)

    print("%2d. %s" % (rank, name))
    print("    value = %.15f, error = %.6f%%" % (value, err))
    print("    m_pred = %.10f GeV, diff = %.1f eV, pull = %.1f sigma" % (m_pred, m_diff, pull))
    print()

# =============================================================================
# INVERSE SYMBOLIC CALCULATOR: what IS c_needed?
# =============================================================================
print("=" * 80)
print("INVERSE SYMBOLIC: What algebraic number IS c_needed?")
print("=" * 80)
print()

c = c_needed
print("c_needed = %.15f" % c)
print()

# Test: is c = phi^k for rational k?
k_phi = np.log(c) / np.log(phi)
print("phi^k test: k = %.10f" % k_phi)
print("  k close to -1? diff = %.6f" % abs(k_phi - (-1)))
print("  k close to -19/20? diff = %.6f" % abs(k_phi - (-19/20)))

# Test: is c = p/q for small p, q?
print()
print("Best rational approximations:")
best_rats = []
for q in range(1, 200):
    p = round(c * q)
    if p > 0:
        err_rat = abs(p/q - c)
        best_rats.append((err_rat, p, q))
best_rats.sort()
for err_rat, p, q in best_rats[:10]:
    print("  %d/%d = %.10f (error = %.2e, %.6f%%)" % (p, q, p/q, err_rat, err_rat/c*100))

# Test: is c^2 = p/q?
print()
print("Best rational approximations for c^2 = %.15f:" % c**2)
c2 = c**2
best_rats2 = []
for q in range(1, 200):
    p = round(c2 * q)
    if p > 0:
        err_rat = abs(p/q - c2)
        best_rats2.append((err_rat, p, q))
best_rats2.sort()
for err_rat, p, q in best_rats2[:5]:
    print("  %d/%d = %.10f (error = %.2e, %.6f%%)" % (p, q, p/q, err_rat, err_rat/c2*100))

# Test: is c related to pi?
print()
print("Pi-related tests:")
for expr_name, expr_val in [
    ("c*pi", c*np.pi), ("c*pi/2", c*np.pi/2), ("c*pi^2", c*np.pi**2),
    ("c/pi", c/np.pi), ("c*sqrt(pi)", c*np.sqrt(np.pi)),
    ("c^2*pi", c**2*np.pi), ("c^2*2*pi", c**2*2*np.pi),
]:
    # Check if expr_val is close to a simple rational or sqrt
    frac = Fraction(expr_val).limit_denominator(100)
    print("  %s = %.10f ~ %s (err %.6f%%)" % (
        expr_name, expr_val, frac, abs(float(frac) - expr_val)/expr_val*100))

# Test: is c related to phi?
print()
print("Phi-related tests:")
for expr_name, expr_val in [
    ("c*phi", c*phi), ("c*phi^2", c*phi**2), ("c*phi^3", c*phi**3),
    ("c/phi", c/phi), ("c*sqrt(phi)", c*np.sqrt(phi)),
    ("c^2*phi", c**2*phi), ("c^2*phi^2", c**2*phi**2),
    ("c*phi + c", c*phi + c), ("1/(c*phi)", 1/(c*phi)),
]:
    frac = Fraction(expr_val).limit_denominator(100)
    print("  %s = %.10f ~ %s (err %.6f%%)" % (
        expr_name, expr_val, frac, abs(float(frac) - expr_val)/expr_val*100))

# =============================================================================
# SPECIFIC TEST: c*phi ~ 1?
# =============================================================================
print()
print("=" * 80)
print("INTERESTING RELATION: c * phi")
print("=" * 80)
print()
print("c_needed * phi = %.15f" % (c * phi))
print("  deviation from 1: %.6e (%.4f%%)" % (c*phi - 1, (c*phi-1)*100))
print()
print("If c = 1/phi, then delta_2 = (alpha/phi) * delta_1")
print("  c_1_over_phi = 1/phi = %.15f" % (1/phi))
print("  c_needed =            %.15f" % c)
print("  difference: %.6e (%.4f%%)" % (c - 1/phi, (c - 1/phi)/c * 100))
print()
print("c_needed/c_sqrt25 = %.15f" % (c / np.sqrt(2.0/5)))
print("  => c_needed is %.6f%% BELOW sqrt(2/5)" % ((1 - c/np.sqrt(2.0/5))*100))
print()

# =============================================================================
# SPECIAL: Renormalization group correction to Morrey exponent?
# =============================================================================
print("=" * 80)
print("ALTERNATIVE INTERPRETATION: Is delta_2 even the right decomposition?")
print("=" * 80)
print()
print("What if the total correction is NOT delta_1 + delta_2")
print("but delta_1 * (1 + c*alpha)?")
print()

# n_total = n_bare + delta_1 * (1 + c*alpha)
# Same algebra, same c. But the PHYSICS is different:
# c*alpha is a MULTIPLICATIVE renormalization of delta_1.
print("Multiplicative form: delta = delta_1 * (1 + c*alpha)")
print("  gives same c_needed = %.10f" % c_needed)
print()
print("But what if the exponent itself is renormalized?")
print("  alpha_eff = 3/4 * (1 + c'*alpha)")
print("  Then delta = c_ell * |z'|^{alpha_eff}")
print()

# Solve for alpha_eff
# delta_total = c_ell * |z'|^{alpha_eff}
# = (phi^3/26) * (1/phi^2)^{alpha_eff}
# = phi^{3 - 2*alpha_eff} / 26
# So: delta_total * 26 = phi^{3 - 2*alpha_eff}
# 3 - 2*alpha_eff = ln(delta_total * 26) / ln(phi)
val_rhs = delta_total * 26
alpha_eff = (3 - np.log(val_rhs) / np.log(phi)) / 2

print("  alpha_eff = %.15f" % alpha_eff)
print("  3/4 = %.15f" % (3.0/4))
print("  deviation = %.6e (%.4f%%)" % (alpha_eff - 3.0/4, (alpha_eff/0.75-1)*100))
print()

c_prime = (alpha_eff / 0.75 - 1) / alpha_em
print("  c' in alpha_eff = (3/4)*(1+c'*alpha) = %.10f" % c_prime)
print("  This c' is DIFFERENT from c = %.10f" % c_needed)
print()

# What if alpha_eff = (d_ST - 1)/(d_ST + alpha/something)?
# alpha_eff = 3/(4 + x) where x is small
# 3/(4+x) = alpha_eff => x = 3/alpha_eff - 4
x_val = 3.0/alpha_eff - 4
print("  If alpha_eff = 3/(4+x): x = %.10f" % x_val)
print("  x/alpha = %.10f" % (x_val/alpha_em))
print("  Is x = alpha * something? x/alpha = %.6f" % (x_val/alpha_em))
print()

# =============================================================================
# DEEP SEARCH: expressions of form f(a1, phi) with small integers
# =============================================================================
print("=" * 80)
print("DEEP SEARCH: f(a1=5, phi, pi) = c_needed?")
print("=" * 80)
print()

deep_hits = []

# sqrt(p/q) for q up to 50
for p in range(1, 50):
    for q in range(1, 50):
        val = np.sqrt(p/q)
        err = abs(val - c_needed) / c_needed
        if err < 0.005:  # within 0.5%
            deep_hits.append((err, "sqrt(%d/%d)" % (p, q), val))

# (p + q*phi)/r for small integers
for p in range(-5, 10):
    for q in range(-5, 5):
        for r in range(1, 20):
            val = (p + q*phi) / r
            if val > 0:
                err = abs(val - c_needed) / c_needed
                if err < 0.001:  # within 0.1%
                    deep_hits.append((err, "(%d + %d*phi)/%d" % (p, q, r), val))

# p*phi^(q/r) for small rationals
for p_num in range(1, 5):
    for p_den in range(1, 5):
        for q_num in range(-5, 5):
            for q_den in range(1, 10):
                val = (p_num/p_den) * phi**(q_num/q_den)
                if val > 0:
                    err = abs(val - c_needed) / c_needed
                    if err < 0.001:
                        name = "(%d/%d)*phi^(%d/%d)" % (p_num, p_den, q_num, q_den)
                        deep_hits.append((err, name, val))

# p/(q*pi) and p*pi/q
for p in range(1, 20):
    for q in range(1, 20):
        for func in [(p/(q*np.pi), "%d/(%d*pi)" % (p, q)),
                      (p*np.pi/q, "%d*pi/%d" % (p, q))]:
            val, name = func
            if 0 < val < 10:
                err = abs(val - c_needed) / c_needed
                if err < 0.002:
                    deep_hits.append((err, name, val))

# sqrt((p + q*sqrt(5))/r) for small integers
for p in range(-10, 20):
    for q in range(-5, 5):
        for r in range(1, 20):
            inner = (p + q*np.sqrt(5)) / r
            if inner > 0:
                val = np.sqrt(inner)
                err = abs(val - c_needed) / c_needed
                if err < 0.001:
                    deep_hits.append((err, "sqrt((%d+%d*sqrt5)/%d)" % (p, q, r), val))

deep_hits.sort()
print("Top hits (within 0.5%% of c_needed = %.10f):" % c_needed)
print()
# Remove duplicates
seen_vals = set()
unique_hits = []
for err, name, val in deep_hits:
    rounded = round(val, 8)
    if rounded not in seen_vals:
        seen_vals.add(rounded)
        unique_hits.append((err, name, val))

for err, name, val in unique_hits[:30]:
    print("  %-40s = %.10f  (%.6f%%)" % (name, val, err*100))

print()
print("=" * 80)
print("END OF EXPERIMENT 432")
print("=" * 80)

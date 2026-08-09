"""
exp362c_connes_exact.py
========================
Find the EXACT D_F weighting on E8~ McKay graph that gives
the Connes ratio m_H/m_W = phi.

Strategy:
1. Parameterize weights by decay w_e = phi^{-beta*k}
2. Find beta that gives Tr(D^4)/Tr(D^2) = phi^2/2 exactly
3. Try physically motivated weightings (CG, dim-based)
4. Check if the solution has algebraic interpretation
"""

import numpy as np
from scipy.optimize import brentq

a1 = 5
b1 = 6
phi = (1.0 + np.sqrt(a1)) / 2.0
alpha = 0.0072973425  # framework
target_ratio = phi**2 / 2  # = 1.309017... (gives m_H/m_W = phi)
target_corrected = (phi - 8*alpha)**2 / 2  # gives m_H/m_W = phi - 8*alpha

print("=" * 72)
print("EXP-362c: EXACT CONNES RATIO ON McKAY GRAPH")
print("=" * 72)
print(f"  phi = {phi:.10f}")
print(f"  phi^2/2 = {target_ratio:.10f}  (tree-level target)")
print(f"  (phi-8a)^2/2 = {target_corrected:.10f}  (corrected target)")
print(f"  8*alpha = {8*alpha:.6f}")

# E8~ McKay graph
dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
n_edges = len(edges)


def connes_ratio(weights):
    """Compute Tr(D^4)/Tr(D^2) for given edge weights on E8~ tree."""
    D = np.zeros((9, 9))
    for idx, (i, j) in enumerate(edges):
        D[i, j] = weights[idx]
        D[j, i] = weights[idx]
    evals = np.linalg.eigvalsh(D)
    t2 = np.sum(evals**2)
    t4 = np.sum(evals**4)
    if t2 < 1e-15:
        return float('inf')
    return t4 / t2


# =====================================================================
# 1. EXPONENTIAL DECAY: w_e = phi^{-beta*k}
# =====================================================================
print("\n" + "=" * 72)
print("1. EXPONENTIAL DECAY: w_e = phi^{-beta*k}")
print("=" * 72)

def ratio_for_beta(beta):
    weights = [phi**(-beta * k) for k in range(n_edges)]
    return connes_ratio(weights)

# Scan beta
print(f"\n  {'beta':>8s}  {'ratio':>12s}  {'2*ratio':>12s}  {'sqrt(2R)':>12s}")
for beta in [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0]:
    r = ratio_for_beta(beta)
    print(f"  {beta:8.2f}  {r:12.6f}  {2*r:12.6f}  {np.sqrt(2*r):12.6f}")

# Find exact beta for phi^2/2
print(f"\n  Target ratio = phi^2/2 = {target_ratio:.10f}")

# Need to find beta where ratio(beta) = phi^2/2
try:
    beta_exact = brentq(lambda b: ratio_for_beta(b) - target_ratio, 0.5, 2.0)
    print(f"  beta_exact = {beta_exact:.10f}")
    print(f"  1/beta = {1/beta_exact:.10f}")
    print(f"  beta*phi = {beta_exact*phi:.10f}")
    print(f"  beta*a1 = {beta_exact*a1:.10f}")
    print(f"  beta/phi = {beta_exact/phi:.10f}")
    print(f"  phi^beta = {phi**beta_exact:.10f}")

    r_check = ratio_for_beta(beta_exact)
    print(f"  Verification: ratio = {r_check:.10f} (target {target_ratio:.10f})")
    print(f"  sqrt(2R) = {np.sqrt(2*r_check):.10f} (target phi = {phi:.10f})")
except Exception as e:
    print(f"  Failed: {e}")
    beta_exact = None

# Also for corrected target
print(f"\n  Corrected target = (phi-8a)^2/2 = {target_corrected:.10f}")
try:
    beta_corr = brentq(lambda b: ratio_for_beta(b) - target_corrected, 0.5, 2.0)
    print(f"  beta_corrected = {beta_corr:.10f}")
    print(f"  Delta beta = {beta_corr - beta_exact:.10f}")
except Exception as e:
    print(f"  Failed: {e}")
    beta_corr = None


# =====================================================================
# 2. DEPTH-BASED DECAY: w_e = phi^{-beta*depth(e)}
# =====================================================================
print("\n" + "=" * 72)
print("2. DEPTH-BASED DECAY: w_e = phi^{-beta*depth(e)}")
print("=" * 72)

# Depth of each edge = min distance from node 0
edge_depths = []
for i, j in edges:
    # distances from node 0 in E8~ tree
    # 0-1-2-3-4-5-6-7, 5-8
    dist = {0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7, 8:6}
    d = min(dist[i], dist[j])
    edge_depths.append(d)

print(f"  Edge depths: {edge_depths}")

def ratio_for_beta_depth(beta):
    weights = [phi**(-beta * d) for d in edge_depths]
    return connes_ratio(weights)

print(f"\n  {'beta':>8s}  {'ratio':>12s}  {'sqrt(2R)':>12s}")
for beta in [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]:
    r = ratio_for_beta_depth(beta)
    print(f"  {beta:8.2f}  {r:12.6f}  {np.sqrt(2*r):12.6f}")

try:
    beta_depth = brentq(lambda b: ratio_for_beta_depth(b) - target_ratio, 0.2, 3.0)
    print(f"\n  beta_depth_exact = {beta_depth:.10f}")
    print(f"  Check: sqrt(2R) = {np.sqrt(2*ratio_for_beta_depth(beta_depth)):.10f}")
except Exception as e:
    print(f"  Failed: {e}")


# =====================================================================
# 3. DIMENSION-BASED WEIGHTINGS
# =====================================================================
print("\n" + "=" * 72)
print("3. DIMENSION-BASED WEIGHTINGS")
print("=" * 72)

def try_scheme(name, weights):
    r = connes_ratio(weights)
    print(f"  {name:40s}  ratio={r:10.6f}  sqrt(2R)={np.sqrt(2*r):10.6f}")
    return r

# Various physically motivated schemes
schemes = {}

# CG-coefficient-like: 1/sqrt(product of dims)
schemes['1/sqrt(d_i*d_j)'] = [1.0/np.sqrt(dims[i]*dims[j]) for i,j in edges]

# Normalized: 1/(d_i*d_j)
schemes['1/(d_i*d_j)'] = [1.0/(dims[i]*dims[j]) for i,j in edges]

# Harmonic mean of dims
schemes['2/(d_i+d_j)'] = [2.0/(dims[i]+dims[j]) for i,j in edges]

# Geometric mean
schemes['sqrt(d_i*d_j)'] = [np.sqrt(dims[i]*dims[j]) for i,j in edges]

# Max dim
schemes['max(d_i,d_j)'] = [float(max(dims[i],dims[j])) for i,j in edges]

# Coxeter labels difference
schemes['|d_i-d_j|+1'] = [float(abs(dims[i]-dims[j])+1) for i,j in edges]

# Product / sum
schemes['d_i*d_j/(d_i+d_j)'] = [dims[i]*dims[j]/(dims[i]+dims[j]) for i,j in edges]

# phi^{-Casimir} where C_2(k) = k(k+2)/4 for SU(2) irrep dim k+1
casimirs = [(d**2 - 1)/4.0 for d in dims]
schemes['phi^{-C_2(min)}'] = [phi**(-min(casimirs[i], casimirs[j])) for i,j in edges]
schemes['phi^{-C_2(max)}'] = [phi**(-max(casimirs[i], casimirs[j])) for i,j in edges]
schemes['phi^{-(C_i+C_j)/2}'] = [phi**(-(casimirs[i]+casimirs[j])/2) for i,j in edges]

# Spin-based: j_k = (d_k - 1)/2 = SU(2) spin
spins = [(d-1)/2.0 for d in dims]
schemes['phi^{-j_max}'] = [phi**(-max(spins[i], spins[j])) for i,j in edges]
schemes['phi^{-j_min}'] = [phi**(-min(spins[i], spins[j])) for i,j in edges]
schemes['phi^{-(j_i+j_j)}'] = [phi**(-(spins[i]+spins[j])) for i,j in edges]

# Perron-Frobenius related
pf = dims / np.sum(dims)  # normalized PF eigenvector
schemes['pf_i*pf_j'] = [pf[i]*pf[j] for i,j in edges]
schemes['sqrt(pf_i*pf_j)'] = [np.sqrt(pf[i]*pf[j]) for i,j in edges]
schemes['1/sqrt(pf_i*pf_j)'] = [1.0/np.sqrt(pf[i]*pf[j]) for i,j in edges]

# Using the COXETER NUMBER h=30
h = 30
schemes['phi^{-k/h}'] = [phi**(-k/h) for k in range(8)]
schemes['phi^{-k/a1}'] = [phi**(-k/a1) for k in range(8)]

for name, w in schemes.items():
    try_scheme(name, w)


# =====================================================================
# 4. CLEBSCH-GORDAN INSPIRED
# =====================================================================
print("\n" + "=" * 72)
print("4. CLEBSCH-GORDAN / WIGNER 6j INSPIRED")
print("=" * 72)

# For SU(2) CG: <j1 m1 j2 m2 | j3 m3>
# The McKay graph edge (k, l) represents rho_std x rho_k -> rho_l
# where rho_std has j = 1/2 (dim 2)
# Wigner 3j squared gives probability
# For j1=1/2, j2=j, j3=j+1/2 or j-1/2:
# |CG|^2 ~ (j+m+1)/(2j+1) or (j-m)/(2j+1)
# Summed over m: overall factor ~ dim_output / dim_input

# Simplified: CG coefficient for std x rho_k -> rho_{k+1}
# has norm ~ sqrt(d_{k+1}/(d_std * d_k)) = sqrt(d_{k+1}/(2*d_k))
# and for std x rho_k -> rho_{k-1}: ~ sqrt(d_{k-1}/(2*d_k))

# Let's try the "branching fraction"
print("\n  Branching fraction model: w_{k->l} = sqrt(d_l / (2*d_k))")
for i, j in edges:
    # edge i-j means rho_std x rho_i has rho_j (and vice versa)
    w_ij = np.sqrt(dims[j] / (2.0 * dims[i]))
    w_ji = np.sqrt(dims[i] / (2.0 * dims[j]))
    print(f"    ({i},{j}): w_ij = {w_ij:.4f}, w_ji = {w_ji:.4f}")

# For symmetric D_F, use geometric mean:
schemes_cg = {}
schemes_cg['CG geom: (d_j*d_i)^{1/4}/sqrt(2*d_i*d_j)^{1/2}'] = [
    (dims[j]*dims[i])**0.25 / np.sqrt(2.0*(dims[i]*dims[j])**0.5) for i,j in edges
]

# Actually simpler: branching ratio for each edge
# In rho_std x rho_k = sum rho_l, the coefficient of rho_l is 0 or 1 (for A-type McKay)
# The CG matrix element has norm 1/sqrt(2*d_k) * sqrt(d_l)
# Symmetric average: sqrt(d_l/(2*d_k) * d_k/(2*d_l)) = 1/2
# This is trivial. The actual CG values depend on m quantum numbers.

# Better approach: for the ICOSAHEDRAL group, use the character theory directly
# The "Frobenius-Schur" type norm: ||CG_{kl}||^2 = delta(McKay) / (d_k * d_std)
# Normalized CG: c_{kl} = 1/sqrt(d_std) = 1/sqrt(2)
# This gives all edges weight 1/sqrt(2) -> same as unit up to scaling!

# Let's try a 1-parameter FAMILY of dim-weighted adjacencies
print("\n  1-parameter family: w_e = (d_i * d_j)^{-gamma}")
print(f"  {'gamma':>8s}  {'ratio':>12s}  {'sqrt(2R)':>12s}")

for gamma in np.arange(-1.0, 2.01, 0.1):
    weights = [(dims[i]*dims[j])**(-gamma) for i,j in edges]
    r = connes_ratio(weights)
    sr = np.sqrt(2*r)
    if abs(sr - phi) < 0.1:
        print(f"  {gamma:8.2f}  {r:12.6f}  {sr:12.6f}  <-- CLOSE to phi!")
    elif abs(gamma % 0.5) < 0.01:
        print(f"  {gamma:8.2f}  {r:12.6f}  {sr:12.6f}")

# Find exact gamma
try:
    gamma_exact = brentq(
        lambda g: connes_ratio([(dims[i]*dims[j])**(-g) for i,j in edges]) - target_ratio,
        0.0, 2.0
    )
    w_exact = [(dims[i]*dims[j])**(-gamma_exact) for i,j in edges]
    r_check = connes_ratio(w_exact)
    print(f"\n  gamma_exact = {gamma_exact:.10f}")
    print(f"  Check: ratio = {r_check:.10f}, sqrt(2R) = {np.sqrt(2*r_check):.10f}")
    print(f"  gamma*2 = {2*gamma_exact:.10f}")
    print(f"  gamma*a1 = {a1*gamma_exact:.10f}")
    print(f"  gamma*phi = {phi*gamma_exact:.10f}")
    print(f"  1/gamma = {1/gamma_exact:.10f}")
except Exception as e:
    print(f"  gamma search failed: {e}")


# =====================================================================
# 5. 1-PARAM FAMILY: w_e = (d_min/d_max)^gamma
# =====================================================================
print("\n" + "=" * 72)
print("5. ASYMMETRY-BASED: w_e = (d_min/d_max)^gamma")
print("=" * 72)

def ratio_asym(gamma):
    weights = [(min(dims[i],dims[j])/max(dims[i],dims[j]))**gamma for i,j in edges]
    return connes_ratio(weights)

print(f"  {'gamma':>8s}  {'ratio':>12s}  {'sqrt(2R)':>12s}")
for gamma in np.arange(0, 3.01, 0.25):
    r = ratio_asym(gamma)
    print(f"  {gamma:8.2f}  {r:12.6f}  {np.sqrt(2*r):12.6f}")

try:
    gamma_asym = brentq(lambda g: ratio_asym(g) - target_ratio, 0, 5)
    print(f"\n  gamma_asym = {gamma_asym:.10f}")
    r_check = ratio_asym(gamma_asym)
    print(f"  Check: sqrt(2R) = {np.sqrt(2*r_check):.10f}")
except Exception as e:
    print(f"  Failed: {e}")


# =====================================================================
# 6. EIGENVALUE ANALYSIS
# =====================================================================
print("\n" + "=" * 72)
print("6. EIGENVALUE STRUCTURE")
print("=" * 72)

# Unit adjacency eigenvalues
A = np.zeros((9, 9))
for i, j in edges:
    A[i, j] = 1
    A[j, i] = 1
evals_A = np.sort(np.linalg.eigvalsh(A))[::-1]
print(f"  Unit adjacency eigenvalues: {[f'{x:.6f}' for x in evals_A]}")

# What eigenvalues would give Tr(lambda^4)/Tr(lambda^2) = phi^2/2?
# For symmetric bipartite: eigenvalues come in +/- pairs plus 0
# Let eigenvalues be {a, b, c, d, 0, -d, -c, -b, -a} with a>b>c>d>0
# Tr(lambda^2) = 2(a^2 + b^2 + c^2 + d^2)
# Tr(lambda^4) = 2(a^4 + b^4 + c^4 + d^4)
# ratio = (a^4+b^4+c^4+d^4)/(a^2+b^2+c^2+d^2) = phi^2/2

# For unit: a=2, b=phi, c=1, d=1/phi
# ratio = (16+phi^4+1+1/phi^4)/(4+phi^2+1+1/phi^2) = (16+7+1)/(4+3+1) = 24/8 = 3
# phi^2/2 = 1.309

# What if we want the "natural" eigenvalues?
# If all |lambda_k| = phi^{p_k} for some powers p_k:
# ratio = sum phi^{4p_k} / sum phi^{2p_k}

# For the STANDARD adjacency, the eigenvalues are:
# 2 ~ phi^{1.44}, phi, 1 = phi^0, 1/phi = phi^{-1}, 0
# So "phi-powers" are roughly {1.44, 1, 0, -1, -inf}

print(f"\n  Eigenvalues in phi-units:")
for ev in evals_A:
    if abs(ev) > 1e-10:
        p = np.log(abs(ev)) / np.log(phi)
        print(f"    {ev:+.6f} = phi^{p:.6f}")
    else:
        print(f"    {ev:+.6f} = 0")

# Now try: what if eigenvalues are exactly {phi^2, phi, 1, 1/phi, 0, ...}?
print(f"\n  If eigenvalues = {{phi^2, phi, 1, 1/phi, 0}} (bipartite):")
evals_test = [phi**2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -phi**2]
t2 = sum(x**2 for x in evals_test)
t4 = sum(x**4 for x in evals_test)
print(f"  Tr(D^2) = {t2:.6f}")
print(f"  Tr(D^4) = {t4:.6f}")
print(f"  Tr(D^4)/Tr(D^2) = {t4/t2:.6f}")
print(f"  sqrt(2R) = {np.sqrt(2*t4/t2):.6f}")

# What about {phi^(3/2), phi, phi^(1/2), phi^0, 0} type?
print(f"\n  If eigenvalues = {{phi^{3/2}, phi, phi^{1/2}, 1, 0}} (bipartite):")
evals_test2 = [phi**1.5, phi, phi**0.5, 1, 0, -1, -phi**0.5, -phi, -phi**1.5]
t2 = sum(x**2 for x in evals_test2)
t4 = sum(x**4 for x in evals_test2)
print(f"  ratio = {t4/t2:.6f}, sqrt(2R) = {np.sqrt(2*t4/t2):.6f}")


# =====================================================================
# 7. KEY TEST: NORMALIZED ADJACENCY D = A / spectral_radius
# =====================================================================
print("\n" + "=" * 72)
print("7. SPECTRAL RADIUS NORMALIZATION")
print("=" * 72)

# Normalize so max eigenvalue = 1 (divide by 2)
D_norm = A / 2.0
evals_norm = np.sort(np.linalg.eigvalsh(D_norm))[::-1]
print(f"  Eigenvalues: {[f'{x:.6f}' for x in evals_norm]}")
t2 = np.sum(evals_norm**2)
t4 = np.sum(evals_norm**4)
print(f"  ratio = {t4/t2:.6f}, sqrt(2R) = {np.sqrt(2*t4/t2):.6f}")

# Normalize so Tr(D^2) = 8 (rank E8)
scale = np.sqrt(8.0 / np.sum(np.linalg.eigvalsh(A)**2))
D_8 = A * scale
evals_8 = np.sort(np.linalg.eigvalsh(D_8))[::-1]
print(f"\n  Tr(D^2) = 8 normalization:")
print(f"  scale = {scale:.6f}")
print(f"  Eigenvalues: {[f'{x:.6f}' for x in evals_8]}")
t2 = np.sum(evals_8**2)
t4 = np.sum(evals_8**4)
print(f"  Tr(D^2) = {t2:.6f}, Tr(D^4) = {t4:.6f}")
print(f"  ratio = {t4/t2:.6f}, sqrt(2R) = {np.sqrt(2*t4/t2):.6f}")
print(f"  (scaling doesn't change the ratio - as expected)")


# =====================================================================
# 8. GOLDEN RATIO IN TRACES: ALGEBRAIC ANALYSIS
# =====================================================================
print("\n" + "=" * 72)
print("8. ALGEBRAIC ANALYSIS: WHEN DOES RATIO = phi^2/2?")
print("=" * 72)

# For E8~ adjacency, eigenvalues {2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2}
# phi^2 + 1/phi^2 = 3, phi^4 + 1/phi^4 = 7

# Tr(A^{2k}) = 2*(2^{2k} + phi^{2k} + 1 + phi^{-2k})
# Let P_k = phi^{2k} + phi^{-2k} = Lucas number L_{2k}
# P_0 = 2, P_1 = 3, P_2 = 7, P_3 = 18, P_4 = 47, P_5 = 123

# Tr(A^2) = 2*(4 + P_1 + 1) = 2*8 = 16
# Tr(A^4) = 2*(16 + P_2 + 1) = 2*24 = 48
# ratio = (16+P_2+1)/(4+P_1+1) = 24/8 = 3

# For weighted adjacency with eigenvalues {a*2, a*phi, a*1, a/phi, 0}:
# ratio = a^2 * 3 (same, scaling invariant)

# For non-trivially weighted adjacency:
# eigenvalues = {mu_0, mu_1, mu_2, mu_3, 0, -mu_3, -mu_2, -mu_1, -mu_0}
# ratio = (mu_0^4+mu_1^4+mu_2^4+mu_3^4)/(mu_0^2+mu_1^2+mu_2^2+mu_3^2) = phi^2/2

# CONSTRAINT: this is Cauchy-Schwarz type.
# By Cauchy-Schwarz: sum(x_i^4)/sum(x_i^2) >= (sum x_i^2)^2 / (n * sum x_i^2) = sum(x_i^2)/n
# where n=4 (nonzero pairs). So ratio >= sum(x_i^2)/4.
# Also: ratio <= max(x_i^2) (since x^4/x^2 = x^2 for dominant term).

# For unit adjacency: min ratio = 8/4 = 2 (if all equal), actual = 3.
# We need phi^2/2 = 1.309... < 2, which is BELOW the all-equal case!
# This means we need the eigenvalues to be MORE concentrated (less spread).

print("  For BIPARTITE eigenvalues {mu_0, mu_1, mu_2, mu_3, 0, -mu_3, -mu_2, -mu_1, -mu_0}:")
print(f"  Unit adjacency: mu = (2, phi, 1, 1/phi) -> ratio = 3")
print(f"  Target: ratio = phi^2/2 = {target_ratio:.6f}")
print(f"")
print(f"  NOTE: Cauchy-Schwarz gives ratio >= sum(mu^2)/4")
print(f"  For unit: sum(mu^2)/4 = 8/4 = 2 > phi^2/2 = {target_ratio:.6f}")
print(f"")
print(f"  THIS MEANS: for 4 nonzero eigenvalue pairs, ratio < 2 is IMPOSSIBLE")
print(f"  unless some eigenvalue pairs have ZERO contribution (i.e., eigenvalues = 0)!")
print(f"")

# Can we achieve phi^2/2 with FEWER nonzero eigenvalues?
# If only 1 pair nonzero: ratio = mu_0^2 (arbitrary)
# If 2 pairs: ratio = (a^4 + b^4)/(a^2 + b^2)
# Let r = b/a: ratio = a^2*(1+r^4)/(1+r^2)

# For phi^2/2 with 2 pairs: a^2*(1+r^4)/(1+r^2) = phi^2/2
# If a = 1: (1+r^4)/(1+r^2) = phi^2/2
# This means: 2(1+r^4) = phi^2*(1+r^2) = (phi^2+phi^2*r^2)
# 2 + 2r^4 = phi^2 + phi^2*r^2
# 2r^4 - phi^2*r^2 + (2-phi^2) = 0
# Let u = r^2: 2u^2 - phi^2*u + (2-phi^2) = 0
# u = (phi^2 +/- sqrt(phi^4 - 8(2-phi^2))) / 4
# = (phi^2 +/- sqrt(phi^4 - 16 + 8*phi^2)) / 4
# phi^4 = 3+2phi, phi^2 = 1+phi
# phi^4 + 8*phi^2 - 16 = 3+2phi + 8+8phi - 16 = -5 + 10phi = 5(2phi-1) = 5*sqrt(5)
# sqrt(5*sqrt(5))... not clean

# Actually let me reconsider. phi^2/2 = 1.309...
# For 4 pairs with ratio = phi^2/2 < 2: IMPOSSIBLE by C-S (as shown above)

# WAIT - the C-S lower bound assumes all 4 are nonzero.
# sum(mu_i^4)/sum(mu_i^2) >= (sum mu_i^2)^2 / (4 * sum mu_i^2) = sum(mu_i^2)/4
# This is wrong! C-S says sum(a_i^2)/sum(b_i^2) in different form.
# Actually: By power mean, sum(x^4)/sum(x^2) = <x^2>_weighted >= min(x^2) when all positive.
# But also sum(x^4)/sum(x^2) = sum(x^2 * x^2)/sum(x^2) which is a weighted average of x^2.
# So phi^2/2 must be in the range [min(mu_i^2), max(mu_i^2)].
# For unit adj: range = [1/phi^2, 4] = [0.382, 4]. phi^2/2 = 1.309 is IN this range.
# So it IS achievable with appropriate weighting!

# The issue is not C-S but rather: CAN we weight the E8~ edges to shift eigenvalues
# such that the ratio becomes phi^2/2?

print(f"  CORRECTION: The ratio is a WEIGHTED average of eigenvalue squares.")
print(f"  Range for unit adj: [{1/phi**2:.4f}, {4:.4f}]")
print(f"  phi^2/2 = {target_ratio:.4f} is WITHIN this range -> achievable!")

# From Section 1: beta_exact gives the exponential decay that achieves this.
if beta_exact is not None:
    print(f"\n  Exponential decay phi^{{-beta*k}} with beta = {beta_exact:.10f} achieves phi exactly.")


# =====================================================================
# 9. PHYSICAL INTERPRETATION: WHAT IS beta?
# =====================================================================
print("\n" + "=" * 72)
print("9. PHYSICAL INTERPRETATION")
print("=" * 72)

if beta_exact is not None:
    print(f"  beta_exact = {beta_exact:.10f}")
    print(f"")
    print(f"  Testing algebraic relations:")

    # Test many algebraic combinations
    tests = {
        '1/phi': 1/phi,
        'phi-1': phi-1,
        '2/phi^2': 2/phi**2,
        '1/phi^2': 1/phi**2,
        'sqrt(2)/phi': np.sqrt(2)/phi,
        '(a1-1)/a1': (a1-1.0)/a1,
        '2/(a1+1)': 2.0/(a1+1),
        'log(2)/log(phi)': np.log(2)/np.log(phi),
        'a1/(a1^2-1)': a1/(a1**2 - 1.0),
        '1/sqrt(a1)': 1/np.sqrt(a1),
        'phi/a1': phi/a1,
        '2*phi/a1': 2*phi/a1,
        'sqrt(phi)/a1': np.sqrt(phi)/a1,
        '1/(phi+1)': 1/(phi+1),
        '1/phi^3': 1/phi**3,
        'alpha_s': 1/(2*phi**3),
        '2*alpha_s': 2/(2*phi**3),
        'sqrt(2/a1)': np.sqrt(2/a1),
        '1/(a1-1)': 1/(a1-1.0),
        '2/a1^2': 2/a1**2,
        'phi/(2*a1)': phi/(2*a1),
        '(phi-1)/phi': (phi-1)/phi,
        'ln(phi)': np.log(phi),
        '2*ln(phi)': 2*np.log(phi),
        'ln(phi)/ln(2)': np.log(phi)/np.log(2),
        'pi/(2*a1)': np.pi/(2*a1),
        'sin(pi/a1)': np.sin(np.pi/a1),
        '1/sqrt(phi*a1)': 1/np.sqrt(phi*a1),
        'phi^2/a1': phi**2/a1,
        'sqrt(1/3)': np.sqrt(1/3),
        '(sqrt(5)-1)/3': (np.sqrt(5)-1)/3,
    }

    print(f"\n  {'Expression':>25s}  {'Value':>12s}  {'Error':>12s}")
    sorted_tests = sorted(tests.items(), key=lambda x: abs(x[1] - beta_exact))
    for name, val in sorted_tests[:15]:
        err = abs(val - beta_exact) / beta_exact
        marker = " <-- MATCH!" if err < 0.01 else (" <--" if err < 0.05 else "")
        print(f"  {name:>25s}  {val:12.8f}  {err:12.6f}{marker}")


# =====================================================================
# 10. EXACT SYMBOLIC CHECK: beta = specific values
# =====================================================================
print("\n" + "=" * 72)
print("10. TESTING SPECIFIC beta VALUES")
print("=" * 72)

# Test beta = values near the exact solution
special_betas = {
    '1': 1.0,
    '1/phi': 1/phi,
    'ln(phi)': np.log(phi),
    'ln(2)/ln(phi)': np.log(2)/np.log(phi),
    '2/pi': 2/np.pi,
    'pi/a1': np.pi/a1,
    '(sqrt(5)-1)/3': (np.sqrt(a1)-1)/3,
    'phi/(phi+2)': phi/(phi+2),
    '3/(2*a1)': 3/(2.0*a1),
}

if beta_exact is not None:
    print(f"  beta_exact = {beta_exact:.10f}")

for name, beta in sorted(special_betas.items(), key=lambda x: abs(x[1] - (beta_exact or 1))):
    r = ratio_for_beta(beta)
    sr = np.sqrt(2 * r)
    err_phi = abs(sr - phi) / phi * 100
    print(f"  beta={name:>20s} = {beta:.8f}:  sqrt(2R) = {sr:.8f}, err = {err_phi:.3f}%")


# =====================================================================
# 11. THE TWO KEY RESULTS
# =====================================================================
print("\n" + "=" * 72)
print("11. SUMMARY OF KEY RESULTS")
print("=" * 72)

print(f"\n  CONNES FORMULA: m_H^2/m_W^2 = 2*Tr(D_F^4)/Tr(D_F^2)")
print(f"  TARGET: m_H/m_W = phi = {phi:.6f} (tree-level)")
print(f"          m_H/m_W = phi - 8*alpha = {phi - 8*alpha:.6f} (corrected)")
print(f"")
print(f"  RESULT 1: Unit McKay adjacency gives sqrt(6) = {np.sqrt(6):.6f} (too high)")
print(f"  RESULT 2: Standard diagonal Yukawa gives ~ 3 (too high)")

if beta_exact is not None:
    print(f"  RESULT 3: Exponential decay phi^{{-{beta_exact:.4f}*k}} gives EXACT phi")
    print(f"")
    print(f"  The decay parameter beta = {beta_exact:.10f}")
    print(f"  Physical meaning: couplings DECAY along McKay chain")
    print(f"  Edge k has weight phi^{{{-beta_exact:.4f}*k}}")
    print(f"  This is the HIERARCHY of Yukawa couplings encoded in graph distance!")

print(f"\n  FULL CHAIN (with 16/15 running):")
m_e = 0.51099895  # MeV
m_Z_fw = m_e * phi**25 * 16/15 / 1000  # GeV
m_W_fw = m_Z_fw * np.sqrt(1 - 6.0/26)
m_H_phi = m_W_fw * phi
m_H_corr = m_W_fw * (phi - 8*alpha)
print(f"  m_Z = m_e*phi^25*16/15 = {m_Z_fw:.2f} GeV (exp: 91.1876)")
print(f"  m_W = m_Z*cos(tW) = {m_W_fw:.2f} GeV (exp: 80.3692)")
print(f"  m_H = m_W*phi = {m_H_phi:.2f} GeV")
print(f"  m_H = m_W*(phi-8a) = {m_H_corr:.2f} GeV (exp: 125.25)")
print(f"  Error: {(m_H_corr - 125.25)/125.25*100:.2f}%")

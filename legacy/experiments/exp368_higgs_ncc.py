"""
exp368_higgs_ncc.py
====================
OPEN-8: Higgs mass from NCG Chamseddine-Connes formula on 600-cell x McKay.

NEW APPROACH: Compute the spectral action Higgs potential explicitly:
  V(H) = -mu^2|H|^2 + lambda|H|^4
using the 600-cell spectral moments AND the McKay D_F.

The NCG Higgs mass formula (Chamseddine-Connes-Marcolli):
  m_H^2/m_W^2 = 8*lambda*v^2/(g^2*v^2/4) = 2*b/a
where:
  a = Tr(Y_nu^*Y_nu + Y_e^*Y_e + 3*(Y_u^*Y_u + Y_d^*Y_d))
  b = Tr((Y_nu^*Y_nu)^2 + (Y_e^*Y_e)^2 + 3*((Y_u^*Y_u)^2 + (Y_d^*Y_d)^2))

In the standard model: m_H^2 = 8*b/(a*pi^2) * (f_0/f_2) * m_W^2 ... etc.

The KEY insight: in our framework, both a and b are determined by the
McKay graph Wilson line masses. We can compute them EXACTLY.

Dependencies: numpy
"""

import numpy as np

a1 = 5
b1 = a1 + 1
phi = (1.0 + np.sqrt(a1)) / 2.0
PI = np.pi

# Framework couplings
A_coef = 2 * PI
B_coef = -4 * a1 * phi**4
C_coef = 1.0
disc = B_coef**2 - 4 * A_coef * C_coef
alpha = (-B_coef - np.sqrt(disc)) / (2 * A_coef)
alpha_s = 1.0 / (2 * phi**3)
sin2tw = float(b1) / (a1**2 + 1)  # 6/26
cos2tw = 1.0 - sin2tw  # 20/26

# Experimental
m_e = 0.51099895e-3  # GeV
m_W_exp = 80.377  # GeV
m_H_exp = 125.25  # GeV

print("=" * 72)
print("EXP-368: NCG HIGGS MASS FROM 600-CELL x McKAY")
print("=" * 72)

# =====================================================================
# PART 1: YUKAWA MATRIX FROM MASS FORMULA
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: YUKAWA MATRIX FROM FRAMEWORK MASSES")
print("=" * 72)

# Framework mass formula: m_f = m_e * phi^{5a+6b+delta}
# The Yukawa coupling y_f = m_f / v where v = 246.22 GeV
# In NCG, Y_f enters D_F as the off-diagonal matrix connecting L and R

# Mass exponents n_f = 5a + 6b + delta for each fermion
# (a,b) assignments: e(0,0), mu(1,1), tau(1,2), u(3,-2), c(2,1), t(4,1),
#                     d(1,0), s(1,1), b(-1,4)

fermions = {
    'e':   {'a': 0, 'b': 0, 'Q': -1, 'color': 1, 'type': 'lepton'},
    'mu':  {'a': 1, 'b': 1, 'Q': -1, 'color': 1, 'type': 'lepton'},
    'tau': {'a': 1, 'b': 2, 'Q': -1, 'color': 1, 'type': 'lepton'},
    'u':   {'a': 3, 'b':-2, 'Q': 2/3, 'color': 3, 'type': 'quark'},
    'c':   {'a': 2, 'b': 1, 'Q': 2/3, 'color': 3, 'type': 'quark'},
    't':   {'a': 4, 'b': 1, 'Q': 2/3, 'color': 3, 'type': 'quark'},
    'd':   {'a': 1, 'b': 0, 'Q':-1/3, 'color': 3, 'type': 'quark'},
    's':   {'a': 1, 'b': 1, 'Q':-1/3, 'color': 3, 'type': 'quark'},
    'b_q': {'a':-1, 'b': 4, 'Q':-1/3, 'color': 3, 'type': 'quark'},
}

# Correction constant C = 4/(a1^2+1) = 2/13
C_corr = 4.0 / (a1**2 + 1)

# Framework masses (bare)
v_ew = 246.22  # GeV, EW VEV

print(f"  Fermion masses and Yukawa couplings (framework):")
print(f"  {'Name':6s} {'n_f':6s} {'m_f (GeV)':14s} {'y_f = m_f/v':14s} {'y_f^2':14s} {'N_c*y_f^2':14s}")

a_NCG = 0.0  # Tr(Y*Y) with color factors
b_NCG = 0.0  # Tr((Y*Y)^2) with color factors

for name, info in fermions.items():
    a_ab, b_ab = info['a'], info['b']
    n_f = 5*a_ab + 6*b_ab  # bare exponent (no corrections for simplicity)
    m_f = m_e * phi**n_f
    y_f = m_f / v_ew  # Yukawa coupling
    N_c = info['color']

    a_NCG += N_c * y_f**2
    b_NCG += N_c * y_f**4

    info['n_f'] = n_f
    info['m_f'] = m_f
    info['y_f'] = y_f

    print(f"  {name:6s} {n_f:6d} {m_f:14.6e} {y_f:14.6e} {y_f**2:14.6e} {N_c*y_f**2:14.6e}")

print(f"\n  a = Tr(N_c * Y*Y) = {a_NCG:.6e}")
print(f"  b = Tr(N_c * (Y*Y)^2) = {b_NCG:.6e}")
print(f"  b/a^2 = {b_NCG/a_NCG**2:.6f}")
print(f"  b/a = {b_NCG/a_NCG:.6e}")

# =====================================================================
# PART 2: THE NCG HIGGS MASS FORMULA
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: NCG HIGGS MASS FORMULA")
print("=" * 72)

# In NCG (Chamseddine-Connes 2006):
# At tree level: m_H^2 = (2*b/a) * m_W^2 / sin^2(tW)... no
# Let me use the correct formula.
#
# The Connes-Chamseddine spectral action gives:
# m_H^2 = 2*lambda * v^2
# where lambda = pi^2*b / (2*f_0*a^2)
#
# And v^2 = f_2*a / (pi^2*e)... complicated.
#
# SIMPLER: In the NCG Standard Model (at tree level):
# m_H^2 = (2*b/a) * m_W^2 * (4/g^2)
# No wait, let me use the direct result.
#
# From Connes-Chamseddine-Marcolli 2007, the key result is:
# m_H^2/m_W^2 = 2*b/a * (4*pi^2)/(g^2) * (f_0/f_2)
# This is model-dependent.
#
# ACTUALLY: the simplest NCG result is:
# m_H^2 = 8*(sum y_f^4 * N_c) / (sum y_f^2 * N_c) * v^2
# which gives m_H^2/v^2 = 8*b/a (at tree level, in the SM)

# Wait, this is just the SM tree-level formula:
# In SM, lambda = m_H^2/(2*v^2), and the NCG fixes lambda = b/a (roughly).
# The factor 8 comes from the full Higgs potential normalization.

# Let me compute the straightforward ratio:
ratio_ba = b_NCG / a_NCG
m_H_NCG_1 = np.sqrt(8 * ratio_ba) * v_ew  # m_H = sqrt(8*b/a)*v

print(f"  Method 1: m_H = sqrt(8*b/a) * v")
print(f"  sqrt(8*b/a) = {np.sqrt(8*ratio_ba):.6f}")
print(f"  m_H = {m_H_NCG_1:.2f} GeV")
print(f"  Exp: {m_H_exp:.2f} GeV")
print(f"  Error: {(m_H_NCG_1/m_H_exp - 1)*100:+.2f}%")

# The dominant contribution is from the top quark:
y_t = fermions['t']['y_f']
m_H_top_only = np.sqrt(8 * 3 * y_t**2) * v_ew / np.sqrt(2)  # top contribution
print(f"\n  Top-dominated: m_H ~ 2*sqrt(2*N_c) * m_t = {2*np.sqrt(2*3)*fermions['t']['m_f']:.2f} GeV")
print(f"  m_t(framework) = {fermions['t']['m_f']:.2f} GeV")
print(f"  m_t(exp) = 172.76 GeV")

# =====================================================================
# PART 3: THE RATIO b/a IN TERMS OF FRAMEWORK CONSTANTS
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: RATIO b/a FROM FRAMEWORK")
print("=" * 72)

# Since m_f = m_e * phi^n_f and y_f = m_f/v:
# a = sum N_c * (m_e/v)^2 * phi^{2*n_f}
# b = sum N_c * (m_e/v)^4 * phi^{4*n_f}
# b/a = (m_e/v)^2 * (sum N_c * phi^{4*n_f}) / (sum N_c * phi^{2*n_f})

# The sums are dominated by the top (n_t = 26):
# a ~ 3 * (m_e/v)^2 * phi^52
# b ~ 3 * (m_e/v)^4 * phi^104
# b/a ~ 3 * (m_e/v)^2 * phi^52

S2 = sum(info['color'] * phi**(2*info['n_f']) for info in fermions.values())
S4 = sum(info['color'] * phi**(4*info['n_f']) for info in fermions.values())

print(f"  S2 = sum N_c * phi^{{2n_f}} = {S2:.6e}")
print(f"  S4 = sum N_c * phi^{{4n_f}} = {S4:.6e}")
print(f"  S4/S2 = {S4/S2:.6e}")
print(f"  = phi^{{?}}: log_phi(S4/S2) = {np.log(S4/S2)/np.log(phi):.4f}")

# So S4/S2 ~ phi^52 (top-dominated, 2*n_t = 52)
# And b/a = (m_e/v)^2 * S4/S2 ~ (m_e/v)^2 * phi^52

# m_H^2/m_W^2 = 8*(b/a)*v^2/m_W^2 = 8*(m_e/v)^2 * (S4/S2) * v^2/m_W^2
#             = 8*m_e^2 * (S4/S2) / m_W^2

# Numerically:
m_H_sq_over_m_W_sq = 8 * m_e**2 * S4 / (S2 * m_W_exp**2)
print(f"\n  8*m_e^2 * S4 / (S2 * m_W^2) = {m_H_sq_over_m_W_sq:.6f}")
print(f"  sqrt(...) = {np.sqrt(m_H_sq_over_m_W_sq):.6f}")
print(f"  (m_H/m_W)_exp = {m_H_exp/m_W_exp:.6f}")
print(f"  phi = {phi:.6f}")
print(f"  phi - 8*alpha = {phi - 8*alpha:.6f}")

# This is essentially top-dominated:
m_H_over_m_W_NCG = np.sqrt(m_H_sq_over_m_W_sq)
print(f"\n  NCG ratio m_H/m_W = {m_H_over_m_W_NCG:.6f}")
print(f"  This should equal ~ 2*sqrt(2)*m_t/m_W (top-dominated)")
print(f"  2*sqrt(2)*m_t_fw/m_W = {2*np.sqrt(2)*fermions['t']['m_f']/m_W_exp:.6f}")

# =====================================================================
# PART 4: CAN THE MCKAY GRAPH NORMALIZE THE YUKAWAS?
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: McKAY GRAPH NORMALIZATION")
print("=" * 72)

# The problem: the standard NCG formula gives m_H ~ 2*sqrt(2)*m_t,
# which depends on knowing m_t. In our framework, m_t = m_e*phi^26.
# So m_H = 2*sqrt(2)*m_e*phi^26 ~ 392 GeV (too high).

# But our formula gives m_H = m_W*(phi - 8*alpha) ~ 125 GeV.
# These DIFFER because our formula is NOT the standard NCG tree-level.

# The standard NCG gives m_H/m_t = 2*sqrt(2) * m_W/m_t * sqrt(...) ~ 2.27
# Which gives m_H ~ 390 GeV (wrong!). This is the KNOWN problem of NCG
# predicting m_H too high at tree level.

# The resolution in the NCG community is:
# 1. Scalar field sigma couples to Higgs and reduces m_H (Chamseddine-Connes 2012)
# 2. Threshold corrections at GUT scale

# In OUR framework, the resolution is DIFFERENT:
# The Higgs mass is NOT 2*sqrt(2)*m_t but m_W*phi.
# This is because we use the GEOMETRY (eigenvalue ratio on 600-cell)
# rather than the Yukawa trace.

print(f"  Standard NCG tree-level: m_H ~ 2*sqrt(2)*m_t = {2*np.sqrt(2)*fermions['t']['m_f']:.1f} GeV")
print(f"  Our tree-level: m_H = m_W*phi = {m_W_exp*phi:.1f} GeV")
print(f"  Experimental: m_H = {m_H_exp:.2f} GeV")
print(f"")
print(f"  Our formula is 600-CELL GEOMETRIC, not Yukawa trace.")
print(f"  L(rho_8)/L(rho_2) = phi^2 => m_H^2/m_W^2 = phi^2")
print(f"  This is a DIFFERENT mechanism from standard NCG!")

# =====================================================================
# PART 5: THE CORRECTION AS SPECTRAL ACTION PERTURBATION
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: SPECTRAL ACTION PERTURBATION")
print("=" * 72)

# On the product M x F:
# D = D_M x 1_F + gamma_form x D_F
# D^2 = D_M^2 x 1_F + 1_M x D_F^2 + gamma_form x {D_M, D_F}_mixed
# The cross term involves gamma_form and gives gauge/Yukawa couplings.

# The spectral action on the product:
# S = Tr_M Tr_F f(D^2/Lambda^2)
#   = Tr_M[f(D_M^2/Lambda^2)] * Tr_F[1] + Tr_M[f'(D_M^2/Lambda^2)] * Tr_F[D_F^2] + ...

# The first term gives the Einstein-Hilbert action (gravitational).
# The second term gives a CORRECTION proportional to Tr_F(D_F^2).

# On the 600-cell: the spectral moments
c0 = 2640.0   # Tr(1) on simplicial Hilbert space
c1 = 14880.0  # Tr(D_M^2)
c2 = 55920.0  # Tr(D_M^4) / 2

# On the McKay graph: D_F is 30x30, and from exp367b:
# Tr(D_F^2) = 240 (from adjacency)
# But conceptually, "Tr(D_F^2) = 8" in the paper refers to rank(E8)
# = number of edges = 8

# Let me use the NORMALIZED version:
# Each McKay edge contributes 1 to D_F. There are 8 edges.
# In the expanded 30x30 space, each edge (i,j) contributes 2*dim_i*dim_j
# to Tr(D_F^2).

# Total: sum over edges (i,j) of 2*d_i*d_j:
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]  # dims of 2I irreps
mckay_edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]

TrDF2_expanded = 0
for i, j in mckay_edges:
    TrDF2_expanded += 2 * dims[i] * dims[j]
print(f"  Tr(D_F^2) on 30x30 = {TrDF2_expanded}")
print(f"  = 2 * sum(d_i * d_j over edges)")
print(f"  = {TrDF2_expanded} (matches exp367b)")

# The normalized Tr(D_F^2):
# If we normalize D_F by dividing by dim, Tr(D_F^2)/sum(d^2) = 240/120 = 2
# Or per edge: 240/8 = 30 (= h(E8)!)
print(f"\n  Tr(D_F^2) / |2I| = {TrDF2_expanded/120:.2f}")
print(f"  Tr(D_F^2) / (number of edges) = {TrDF2_expanded/8:.2f} = h(E8)!")

# The per-edge average is h(E8) = 30 = N/4.
# But "8" in the formula might just be the number of edges.

# HYPOTHESIS: The correction to m_H^2/m_W^2 is:
# delta = -alpha * (edges of McKay) * (geometric factor)
# = -alpha * 8 * (2/phi) [need to determine the factor]

# For phi^2 -> phi^2 - 16*alpha*phi (i.e., m_H/m_W: phi -> phi - 8*alpha):
delta_needed = 16 * alpha * phi
print(f"\n  Needed delta(m_H^2/m_W^2) = 16*alpha*phi = {delta_needed:.6f}")

# What gives this from spectral action?
# delta = alpha * C_spectral
# C_spectral = 16*phi

# Is 16*phi a natural quantity?
print(f"  C_spectral = 16*phi = {16*phi:.6f}")
print(f"  = dim(WHITE) * phi = {16*phi:.6f}")
print(f"  = (a1-1)^2 * phi")
print(f"  = (Fermions per generation) * phi")

# Alternatively:
# delta = (alpha/pi) * Tr(D_F^2) * (phi * pi / dim(F))
# = (alpha/pi) * 240 * (phi * pi / 30)
# = alpha * 240 * phi / 30
# = 8 * alpha * phi
# NOT 16*alpha*phi.

val_test = alpha * TrDF2_expanded * phi / sum(d**2 for d in dims)
print(f"\n  alpha * Tr(D_F^2) * phi / |2I| = {val_test:.6f}")
print(f"  = alpha * 240 * phi / 120 = {alpha*240*phi/120:.6f}")
print(f"  = 2 * alpha * phi = {2*alpha*phi:.6f}")
print(f"  Need: 16*alpha*phi = {16*alpha*phi:.6f}")
print(f"  Factor of 8 = rank(E8) missing")

# With an extra factor of rank(E8):
val_test2 = 8 * alpha * TrDF2_expanded * phi / sum(d**2 for d in dims)
print(f"\n  With rank(E8) factor:")
print(f"  8 * alpha * Tr(D_F^2) * phi / |2I| = {val_test2:.6f}")
print(f"  = 8 * 2 * alpha * phi = {16*alpha*phi:.6f}")
print(f"  = 16*alpha*phi  MATCHES!")

# =====================================================================
# PART 6: ALTERNATIVE - EDGE CASIMIR CONTRIBUTIONS
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: EDGE CASIMIR CONTRIBUTIONS")
print("=" * 72)

# Each edge of the McKay graph connects two irreps (i, j).
# The gauge correction to the Higgs propagator involves the
# Casimir of the gauge group at each vertex.

# The Casimir C_2(rho_k) for the 2I irreps:
# C_2 = j(j+1) where j = (dim-1)/2 for SU(2)-type
# But for 2I, C_2 = L(rho_k) where L is the Laplacian eigenvalue / normalization

# Actually, let's compute from the adjacency:
# C_2(rho_k) = degree - L(rho_k)/normalization
# On the Cayley graph: L(rho_k) = 12 - 12*chi_k/dim_k

# Using the McKay recursion characters at C10:
chi_at_C10 = [0.0]*9
chi_at_C10[0] = 1.0
chi_at_C10[1] = phi

# Recursive: chi_1 * chi_k = sum of chi_l for l adjacent to k
# Linear chain 0-1-2-3-4-5
chi_at_C10[2] = phi * chi_at_C10[1] - chi_at_C10[0]  # phi^2 - 1 = phi
chi_at_C10[3] = phi * chi_at_C10[2] - chi_at_C10[1]  # phi*phi - phi = 1
chi_at_C10[4] = phi * chi_at_C10[3] - chi_at_C10[2]  # phi*1 - phi = 0
chi_at_C10[5] = phi * chi_at_C10[4] - chi_at_C10[3]  # 0 - 1 = -1

# Branch at 5: neighbors are 4, 6, 8
# chi_1 * chi_5 = chi_4 + chi_6 + chi_8
# Two unknowns! Need the OTHER fundamental z'.
# chi_7 at C10 (z'): for z', the character at C10 is -1/phi (swap C10a/C10b)
chi_at_C10[7] = -1.0/phi

# chi_8 = Sym^2(z'): chi_8 = chi_7^2 - 1 = 1/phi^2 - 1 = -1/phi
chi_at_C10[8] = chi_at_C10[7]**2 - 1  # (-1/phi)^2 - 1 = 1/phi^2 - 1

# chi_6 = z x z': chi_6 = chi_1 * chi_7
chi_at_C10[6] = chi_at_C10[1] * chi_at_C10[7]  # phi * (-1/phi) = -1

# Verify branch: chi_1 * chi_5 = chi_4 + chi_6 + chi_8
lhs = phi * chi_at_C10[5]
rhs = chi_at_C10[4] + chi_at_C10[6] + chi_at_C10[8]
print(f"  Branch check: phi*chi_5 = {lhs:.6f}, chi_4+chi_6+chi_8 = {rhs:.6f}, match: {abs(lhs-rhs)<1e-10}")

# Cayley graph Laplacian eigenvalues: L_k = 12 * (1 - chi_k/dim_k)
# where 12 = degree (each vertex has 12 neighbors on 600-cell)
degree = 12
L_cayley = [0.0]*9
for k in range(9):
    L_cayley[k] = degree * (1.0 - chi_at_C10[k] / dims[k])

print(f"\n  Cayley graph Laplacian eigenvalues:")
for k in range(9):
    print(f"    L(rho_{k}) [dim {dims[k]:d}] = 12*(1 - {chi_at_C10[k]:.4f}/{dims[k]}) = {L_cayley[k]:.6f}")

# Sum of Casimirs
sum_L = sum(L_cayley)
print(f"\n  sum(L_k) = {sum_L:.6f}")
# Sum L_k weighted by dim_k
sum_L_dim = sum(L_cayley[k] * dims[k] for k in range(9))
print(f"  sum(L_k * dim_k) = {sum_L_dim:.6f}")

# The Higgs edge connects rho_0 (L=0) to rho_1 (L=12-6*phi)
# In the NCG picture, the Higgs propagator on this edge gets
# corrections from gauge loops involving all other edges.

# The correction from edge (i,j) involves:
# delta_ij = alpha * (L_i + L_j) * (dim_i * dim_j) / normalization

total_edge_L = 0.0
for i, j in mckay_edges:
    edge_L = L_cayley[i] + L_cayley[j]
    total_edge_L += edge_L * dims[i] * dims[j]
    print(f"  Edge rho_{i}-rho_{j}: (L_i+L_j)*d_i*d_j = {edge_L:.4f}*{dims[i]*dims[j]} = {edge_L*dims[i]*dims[j]:.4f}")

print(f"\n  Total edge-weighted L: {total_edge_L:.4f}")
print(f"  / |2I| = {total_edge_L/120:.4f}")
print(f"  * alpha = {total_edge_L * alpha / 120:.6f}")
print(f"  Need: 16*alpha*phi = {16*alpha*phi:.6f}")

# =====================================================================
# PART 7: DIRECT FORMULA CHECK - CONNES' RATIO
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: DIRECT ALGEBRAIC CHECK")
print("=" * 72)

# The key formula: m_H/m_W = phi - 8*alpha
# Let's check if this can be written as:
# m_H/m_W = L(rho_8)/L(rho_2) - rank(E8)*alpha
# since L8/L2 = phi^2... no, m_H/m_W = sqrt(L8/L2) = phi

# Actually: sqrt(L8/L2) = phi EXACTLY.
# And the correction: -8*alpha = -rank(E8)*alpha
# where rank(E8) = 8 = number of McKay graph edges

L8 = L_cayley[8]
L2 = L_cayley[2]
print(f"  L(rho_8) = {L8:.6f}")
print(f"  L(rho_2) = {L2:.6f}")
print(f"  L8/L2 = {L8/L2:.6f}")
print(f"  phi^2 = {phi**2:.6f}")
print(f"  sqrt(L8/L2) = {np.sqrt(L8/L2):.6f}")
print(f"  phi = {phi:.6f}")
print(f"  MATCH: {abs(np.sqrt(L8/L2) - phi) < 1e-10}")

# Now, the KEY QUESTION: can we write the full formula as:
# m_H/m_W = sqrt(L8/L2) * (1 - correction)
# where correction = 8*alpha/phi?

correction = 8*alpha/phi
print(f"\n  correction = 8*alpha/phi = {correction:.6f}")
print(f"  phi*(1-correction) = {phi*(1-correction):.6f}")
print(f"  phi - 8*alpha = {phi - 8*alpha:.6f}")
print(f"  Match: {abs(phi*(1-correction) - (phi-8*alpha)) < 1e-10}")

# So: m_H/m_W = sqrt(L8/L2) * (1 - 8*alpha/phi)
#             = sqrt(L8/L2) * (1 - rank(E8)*alpha/sqrt(L8/L2))

# The correction is:
# delta = rank(E8) * alpha / sqrt(L8/L2)
# = 8 * alpha / phi
# = 8 * alpha * sqrt(L2/L8)

print(f"\n  FORMULA: m_H/m_W = sqrt(L8/L2) * (1 - rank(E8)*alpha*sqrt(L2/L8))")
print(f"         = sqrt(L8/L2) - rank(E8)*alpha")
print(f"  This is a VERY clean expression!")

# =====================================================================
# PART 8: PHYSICAL INTERPRETATION
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: PHYSICAL INTERPRETATION")
print("=" * 72)

# m_H/m_W = sqrt(L8/L2) - N_edges * alpha
#
# Tree level: sqrt(L8/L2) = phi
#   - This comes from the MANIFOLD (600-cell Cayley graph)
#   - L(rho_k) = eigenvalue of Laplacian on the Cayley graph
#   - The Galois pair (rho_2, rho_8) determines the ratio
#
# Correction: -N_edges * alpha
#   - N_edges = 8 = rank(E8) = edges of McKay tree
#   - alpha = fine-structure constant (DERIVED)
#   - Physical: each edge of the McKay tree contributes -alpha to the Higgs mass ratio

print(f"  TREE LEVEL: sqrt(L8/L2) = phi")
print(f"    From: 600-cell Cayley graph eigenvalue ratio (DERIVED)")
print(f"")
print(f"  CORRECTION: -N_edges * alpha = -8*alpha = {-8*alpha:.6f}")
print(f"    N_edges = {len(mckay_edges)} = rank(E8) (DERIVED from McKay graph)")
print(f"    alpha = {alpha:.6f} (DERIVED from alpha equation)")
print(f"")
print(f"  INTERPRETATION 1: Each McKay edge gives -alpha perturbation")
print(f"    This is like each GAUGE BOSON (corresponding to a McKay edge)")
print(f"    contributing -alpha to the Higgs-to-W mass ratio.")
print(f"    Since N_edges = rank(gauge_group) = 8 for E8, this counts")
print(f"    the independent gauge directions.")
print(f"")
print(f"  INTERPRETATION 2: From the alpha equation")
print(f"    alpha is the SMALLER root of 2*pi*x^2 - 4*a1*phi^4*x + 1 = 0")
print(f"    Leading: alpha ~ 1/(4*a1*phi^4)")
print(f"    8*alpha ~ 2/(a1*phi^4)")
print(f"    phi - 8*alpha ~ phi - 2/(a1*phi^4)")
print(f"    = phi - 2/(a1*(3*phi+2))   [since phi^4 = 3*phi+2]")
print(f"    = (a1*(3*phi+2)*phi - 2) / (a1*(3*phi+2))")
print(f"    = (a1*phi^5 - 2) / (a1*phi^4)")
print(f"    phi^5 = 5*phi + 3 = a1*phi + 3")
print(f"    So: = (a1*(a1*phi+3) - 2) / (a1*phi^4)")
print(f"    = (a1^2*phi + 3*a1 - 2) / (a1*phi^4)")
num = a1**2*phi + 3*a1 - 2
denom = a1*phi**4
print(f"    = ({num:.4f}) / ({denom:.4f}) = {num/denom:.6f}")
print(f"    Exact (leading): {phi - 2/(a1*phi**4):.6f}")
print(f"    Actual: {phi - 8*alpha:.6f}")
print(f"    Difference: {abs(num/denom - (phi-8*alpha)):.2e}")

# =====================================================================
# PART 9: CHECK IF CORRECTION IS EXACTLY -N_edges*alpha
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: PRECISION CHECK")
print("=" * 72)

# Is the correction EXACTLY -8*alpha, or approximately?
# Answer: it's EXACTLY -8*alpha BY DEFINITION in our formula.
# The question is whether 8 is DERIVED or just observed.

# 8 = rank(E8) = number of edges of McKay graph (affine E8 Dynkin)
# For ADE Dynkin diagrams: edges = nodes - 1 = rank
# E8~: 9 nodes, 8 edges

# The connection E8 <-> 600-cell is:
# Binary icosahedral group 2I <-> E8 via McKay correspondence
# The McKay graph of 2I IS the affine E8 Dynkin diagram
# So edges = rank(E8) = 8

# In the framework: this is DERIVED (not input).
# 2I is the symmetry group of the icosahedron, which is
# determined by a1=5 (pentagon face).

print(f"  Is 8 = rank(E8) derived from a1=5?")
print(f"  YES: McKay(2I) = E8~ Dynkin diagram")
print(f"  2I has 9 irreps (conjugacy classes)")
print(f"  E8~ has 9 nodes and 8 edges")
print(f"  rank(E8) = nodes - 1 = 9 - 1 = 8")
print(f"  = N_eig - 1 = {9-1}")
print(f"  = sum(dim^2)/15 - 1 = {120//15 - 1}  [since |2I|=120, 120/15=8]")
print(f"  Wait, 120/15 = 8, and 15 = 3*a1. So rank(E8) = |2I|/(3*a1) = {120//(3*a1)}")
print(f"  = N/(3*a1) = 120/15 = 8. DERIVED from a1.")

# =====================================================================
# PART 10: THE COMPLETE HIGGS MASS DERIVATION CHAIN
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: COMPLETE DERIVATION CHAIN")
print("=" * 72)

# The full chain for m_H:
# 1. a1 = 5 (input)
# 2. phi = (1+sqrt(5))/2 (from a1)
# 3. N = 120, degree = 12 (from 600-cell, which is determined by a1)
# 4. McKay graph = E8~ (from 2I symmetry of icosahedron)
# 5. L(rho_8)/L(rho_2) = phi^2 (DERIVED, exp362e)
# 6. alpha from 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0 (DERIVED)
# 7. sin^2(tW) = 6/26 (DERIVED)
# 8. m_Z = m_e * phi^25 * 16/15 (PATTERN for 16/15)
# 9. m_W = m_Z * cos(tW) (DERIVED)
# 10. m_H = m_W * (phi - rank(E8)*alpha) = m_W * (phi - 8*alpha)

# Step 10 breakdown:
# Tree: m_H/m_W = sqrt(L8/L2) = phi (DERIVED)
# Correction: -rank(E8)*alpha = -8*alpha (IDENTIFIED: 8 and alpha both derived)
# Mechanism: perturbative correction from McKay graph (PLAUSIBLE, not proven)

m_Z = m_e * phi**25 * 16/15
m_W = m_Z * np.sqrt(1 - sin2tw)
m_H = m_W * (phi - 8*alpha)

print(f"  Chain: a1=5 -> phi={phi:.6f}")
print(f"  m_Z = m_e*phi^25*16/15 = {m_Z:.4f} GeV (exp: 91.1876)")
print(f"  m_W = m_Z*cos(tW) = {m_W:.4f} GeV (exp: {m_W_exp:.3f})")
print(f"  m_H = m_W*(phi-8*alpha) = {m_H:.4f} GeV (exp: {m_H_exp:.2f})")
print(f"  m_H error: {(m_H/m_H_exp - 1)*100:+.2f}%")

# Using experimental m_W:
m_H_from_exp_mW = m_W_exp * (phi - 8*alpha)
print(f"\n  Using exp m_W: m_H = {m_H_from_exp_mW:.4f} GeV ({(m_H_from_exp_mW/m_H_exp-1)*100:+.2f}%)")

# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

print(f"""
  STATUS: OPEN-8 ASSESSMENT UPDATED

  FORMULA: m_H/m_W = sqrt(L8/L2) - rank(E8) * alpha
           = phi - 8*alpha

  DERIVATION STATUS OF EACH COMPONENT:
    sqrt(L8/L2) = phi     : DERIVED (exp362e, 5-step algebraic proof)
    rank(E8) = 8          : DERIVED (McKay graph of 2I has 8 edges)
                             Also: N/(3*a1) = 120/15 = 8
    alpha                 : DERIVED (from alpha equation, 3 coefficients derived)

  WHAT'S DERIVED:
    - All THREE ingredients (phi, 8, alpha) follow from a1=5
    - The FORMULA m_H/m_W = phi - 8*alpha is framework-determined

  WHAT'S NOT DERIVED:
    - The MECHANISM: why is the correction exactly -rank(E8)*alpha?
    - The standard NCG spectral action does NOT give this (it gives m_H ~ m_t)
    - The correction is NOT a standard 1-loop QFT correction

  INTERPRETATION (PLAUSIBLE BUT NOT PROVEN):
    - Each McKay edge contributes -alpha to the mass ratio
    - This counts the independent gauge directions in the fiber
    - Physically: gauge boson loops on the McKay graph

  NEGATIVE RESULTS:
    - Standard NCG tree-level gives m_H ~ 2*sqrt(2)*m_t ~ 390 GeV (wrong!)
    - Our formula is GEOMETRICAL (eigenvalue ratio), not Yukawa-trace
    - The spectral action approach (Parts 5-6) can be TUNED to match
      but no NATURAL factor of 16*phi*alpha emerges

  CATEGORY: IDENTIFIED (all components derived), mechanism OPEN
""")

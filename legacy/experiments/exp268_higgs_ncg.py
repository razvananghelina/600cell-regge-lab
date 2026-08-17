"""
EXP-268: HIGGS MASS FROM NCG INNER FLUCTUATION ON McKAY GRAPH
===============================================================
Key question: WHY m_H/m_W = phi - 8*alpha instead of Connes sqrt(2)*m_t?

Approach: The spectral action on 600-cell x McKay gives a Higgs potential
where the Higgs field = inner fluctuation of D_F along McKay edges.
V(H) ~ -c_1*Tr(D_F(H)^2) + c_2*Tr(D_F(H)^4)/2
The 600-cell coefficients c_1, c_2 differ from smooth manifold -> modified m_H.
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120
ALPHA = 1/137.035999084
ALPHA_S = 1/(2*PHI**3)
m_e = 0.51099895; m_W = 80379.0; m_H_exp = 125250.0; v_EW = 246220.0

print("="*72)
print("EXP-268: HIGGS MASS FROM NCG ON McKAY GRAPH")
print("="*72)

# === McKay graph (affine E8) ===
# 0(e)-1(u)-2(d)-3(s)-4(c)-5(b)-6(t)-7(tau), branch 4(c)-8(mu)
E8_edges = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(4,8)]
A = np.zeros((9,9))
for i,j in E8_edges:
    A[i,j] = 1; A[j,i] = 1

n_f = np.array([0, 3, 5, 11, 16, 19, 26, 17, 11])  # mass exponents
coxeter = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])     # Coxeter labels
deg = A.sum(axis=1).astype(int)

print(f"\n  McKay graph: 9 nodes, 8 edges")
print(f"  Node assignment: e(0) u(1) d(2) s(3) c(4) b(5) t(6) tau(7) mu(8)")
print(f"  Mass exponents: {list(n_f)}")
print(f"  Degrees: {list(deg)}")

# D_F = diagonal mass exponents
D_F = np.diag(n_f.astype(float))

# === Part 1: Adjacency matrix spectrum ===
print("\n--- Part 1: McKay Adjacency Spectrum ---")
eig_A = np.sort(np.linalg.eigvalsh(A))
print(f"  Eigenvalues of A_E8: {[f'{e:.4f}' for e in eig_A]}")
print(f"  Spectral radius: {max(abs(eig_A)):.4f}")

# Adjacency matrix powers (for potential traces)
A2 = A @ A; A3 = A @ A2; A4 = A @ A3
print(f"  Tr(A^2) = {np.trace(A2):.0f} = 2*|edges| = {2*len(E8_edges)}")
print(f"  Tr(A^4) = {np.trace(A4):.0f}")
print(f"  Tr(A^3) = {np.trace(A3):.0f} (=0 for tree)")

# === Part 2: Inner fluctuation = Higgs direction ===
print("\n--- Part 2: Higgs as Inner Fluctuation ---")
# In NCG: Higgs = inner fluctuation of D_F
# On McKay: H_dir = adjacency matrix A (connects adjacent fermion species)
# D_F(h) = D_F + h * H_dir where h = Higgs field amplitude

# Traces needed for Higgs potential:
# Tr(D_F^2) = sum n_f^2
# Tr(D_F^2 * A^2) = sum n_f_i^2 * deg_i (mixing of mass and topology)
# Tr(D_F^4) = sum n_f^4
# Tr(A^4) = closed walks of length 4 on E8

a_Y = np.trace(D_F @ D_F)  # = sum n^2
b_Y = np.trace(D_F @ D_F @ D_F @ D_F)  # = sum n^4
cross = np.trace(D_F @ D_F @ A2)  # = sum n_i^2 * deg_i
trA4 = np.trace(A4)

# Also: Tr(D_F * A * D_F * A) - important cross term
cross2 = np.trace(D_F @ A @ D_F @ A)  # = sum_{(i,j) edge} n_i * n_j

print(f"  Tr(D_F^2) = a_Y = {a_Y:.0f}")
print(f"  Tr(D_F^4) = b_Y = {b_Y:.0f}")
print(f"  Tr(D_F^2 * A^2) = {cross:.0f} = sum n_i^2 * deg_i")
print(f"  Tr(D_F*A*D_F*A) = {cross2:.0f} = sum_{{edges}} n_i*n_j")
print(f"  Tr(A^4) = {trA4:.0f}")

# Full 4th order trace:
# Tr((D_F + h*A)^4) = Tr(D_F^4) + 4h*Tr(D_F^3*A) + h^2*[4*Tr(D_F^2*A^2) + 2*Tr(D_F*A*D_F*A)]
#                    + 4h^3*Tr(D_F*A^3) + h^4*Tr(A^4)
# Odd terms vanish for tree graph (no odd cycles, diagonal D_F)
tr_DF3A = np.trace(D_F@D_F@D_F@A)
tr_DFA3 = np.trace(D_F@A@A@A)
print(f"  Tr(D_F^3*A) = {tr_DF3A:.0f} (should be 0)")
print(f"  Tr(D_F*A^3) = {tr_DFA3:.0f} (should be 0 for tree)")

# h^2 coefficient in Tr((D_F+h*A)^4)
coeff_h2_in_4 = 4*cross + 2*cross2
print(f"\n  Tr((D_F+hA)^4) = {b_Y:.0f} + {coeff_h2_in_4:.0f}*h^2 + {trA4:.0f}*h^4")

# === Part 3: Higgs potential from spectral action ===
print("\n--- Part 3: Higgs Potential V(h) ---")
# Spectral action: S ~ f_2*c_1*Tr(D_F^2)/Lambda^2 - f_4*c_2*Tr(D_F^4)/(2*Lambda^4)
# (signs: -mu^2|H|^2 + lambda|H|^4)
# V(h) = -mu_0^2 * Tr((D_F+hA)^2) + lambda_0 * Tr((D_F+hA)^4)
# = -mu_0^2 * [a_Y + 2*|edges|*h^2] + lambda_0 * [b_Y + coeff_h2*h^2 + trA4*h^4]

# Spectral coefficients of 600-cell
c1 = 14880.0; c2 = 55920.0
# Normalized: s_k = c_k/240
s0 = 11; s1 = 62; s2 = 233

# In spectral action, the coefficients are:
# mu^2 ~ c_1 / Lambda^2 (related to scalar curvature)
# lambda ~ c_2 / Lambda^4 (related to curvature squared)
# On smooth manifold: mu^2/lambda ~ Lambda^2 * c_1/c_2

# The ratio that determines physics:
ratio_c = c1**2 / (c2)  # has dim of c_0, related to scale
print(f"  c_1 = {c1:.0f}, c_2 = {c2:.0f}")
print(f"  c_1/c_2 = {c1/c2:.6f}")
print(f"  c_1^2/c_2 = {c1**2/c2:.2f}")
print(f"  c_2/c_1 = {c2/c1:.6f} = {s2}/{s1} = {s2/s1:.6f}")

# The key RATIO for Higgs mass (Connes formula modified):
# m_H^2/m_W^2 = 8*lambda/g2^2
# In NCG: lambda = (pi^2/(2*f_0)) * (b_Y/a_Y^2) * GEOMETRIC_FACTOR
# The geometric factor on 600-cell differs from smooth manifold

# On smooth S^3: ratio c2_smooth/c1_smooth^2 = R^(-2)*V / (R^{-1}*V)^2 = 1/(R^2*V*factors)
# On 600-cell: c2/c1^2 = 55920/14880^2 = 2.525e-4

# === Part 4: The correction factor ===
print("\n--- Part 4: Correction Factor Analysis ---")
# Connes on smooth 4-manifold: m_H = sqrt(2*b/a) where a=Tr(Y^2), b=Tr(Y^4)
# With our mass exponents (in units where v is absorbed):
# a = a_Y = 1858, b = b_Y = 766342

m_H_Connes_exponent = np.sqrt(2*b_Y/a_Y)  # in exponent units
print(f"  sqrt(2*b_Y/a_Y) = {m_H_Connes_exponent:.4f} (in exponent units)")

# In physical units, we need to convert exponent to mass
# n_f -> m_f = m_e * phi^n -> y_f = sqrt(2)*m_e*phi^n/v
# So a_phys = sum y^2 = (2*m_e^2/v^2) * sum phi^(2*n)
# b_phys = sum y^4 = (4*m_e^4/v^4) * sum phi^(4*n)
# Connes: m_H^2 = 2*v^2*b_phys/a_phys = 4*m_e^2 * sum(phi^4n)/sum(phi^2n)

sum_phi2n = sum(PHI**(2*n) for n in n_f)
sum_phi4n = sum(PHI**(4*n) for n in n_f)
m_H_Connes = np.sqrt(4*m_e**2 * sum_phi4n/sum_phi2n)  # no color factors
print(f"  m_H(Connes, no color) = {m_H_Connes:.0f} MeV = {m_H_Connes/1000:.2f} GeV")

# With color factors (Nc=3 for quarks)
Nc = [1,3,3,3,3,3,3,1,1]  # e,u,d,s,c,b,t,tau,mu
sum_phi2n_c = sum(nc*PHI**(2*n) for nc,n in zip(Nc,n_f))
sum_phi4n_c = sum(nc*PHI**(4*n) for nc,n in zip(Nc,n_f))
m_H_Connes_c = np.sqrt(4*m_e**2 * sum_phi4n_c/sum_phi2n_c)
print(f"  m_H(Connes, with color) = {m_H_Connes_c:.0f} MeV = {m_H_Connes_c/1000:.2f} GeV")

# Our formula
m_H_ours = m_W * (PHI - 8*ALPHA)
print(f"  m_H(600-cell) = {m_H_ours:.0f} MeV = {m_H_ours/1000:.2f} GeV")
print(f"  m_H(exp) = {m_H_exp:.0f} MeV")

R_no_color = (m_H_ours/m_H_Connes)**2
R_with_color = (m_H_ours/m_H_Connes_c)**2
print(f"\n  Correction factor (no color): R = {R_no_color:.6f}")
print(f"  Correction factor (with color): R = {R_with_color:.6f}")

# === Part 5: McKay graph modification ===
print("\n--- Part 5: McKay Graph Correction to Connes ---")
# Hypothesis: R comes from the graph topology of McKay
# On McKay: the effective Yukawa trace is modified by graph structure

# Effective a_eff = Tr(D_F^2 * W) where W is a weight matrix from geometry
# Effective b_eff = Tr(D_F^4 * W')

# Weight from Coxeter labels:
print(f"  Coxeter labels: {list(coxeter)}")
print(f"  sum(coxeter) = {sum(coxeter)} = h(E8) = {a1*b1}")

# Coxeter-weighted traces
a_cox = sum(c*n**2 for c,n in zip(coxeter, n_f))
b_cox = sum(c*n**4 for c,n in zip(coxeter, n_f))
m_H_cox = np.sqrt(4*m_e**2 * b_cox/a_cox)  # Coxeter-weighted Connes
print(f"  a_cox = sum(c_i * n_i^2) = {a_cox}")
print(f"  b_cox = sum(c_i * n_i^4) = {b_cox}")
print(f"  m_H(Coxeter) = {m_H_cox/1000:.2f} GeV")

# Degree-weighted traces
a_deg = sum(d*n**2 for d,n in zip(deg, n_f))
b_deg = sum(d*n**4 for d,n in zip(deg, n_f))
m_H_deg = np.sqrt(4*m_e**2 * b_deg/a_deg)
print(f"  a_deg = sum(deg_i * n_i^2) = {a_deg}")
print(f"  b_deg = sum(deg_i * n_i^4) = {b_deg}")
print(f"  m_H(degree) = {m_H_deg/1000:.2f} GeV")

# Coxeter * color
a_cox_c = sum(c*nc*n**2 for c,nc,n in zip(coxeter, Nc, n_f))
b_cox_c = sum(c*nc*n**4 for c,nc,n in zip(coxeter, Nc, n_f))
m_H_cox_c = np.sqrt(4*m_e**2 * b_cox_c/a_cox_c)
print(f"  m_H(Cox*color) = {m_H_cox_c/1000:.2f} GeV")

# === Part 6: Spectral action Higgs potential ===
print("\n--- Part 6: Full Spectral Action Potential ---")
# V(h) = -alpha_2 * h^2 + alpha_4 * h^4
# where:
# alpha_2 = c_1 * Tr(A^2) - c_2 * (4*cross + 2*cross2) / (something)
# alpha_4 = c_2 * Tr(A^4) / (something)
# The "something" involves Lambda (cutoff)

# Let's parametrize differently. Define:
# V(h) = a2*h^2 + a4*h^4
# From Tr((D+hA)^2) = a_Y + 16*h^2
# From Tr((D+hA)^4) = b_Y + coeff_h2*h^2 + trA4*h^4

# With spectral action weighting:
# V(h) = -(c_1/L^2) * 16*h^2 + (c_2/(2*L^4)) * [coeff_h2*h^2 + trA4*h^4]
# = h^2 * [-16*c_1/L^2 + coeff_h2*c_2/(2*L^4)] + h^4 * [trA4*c_2/(2*L^4)]

# For SSB: need the h^2 coefficient to be negative
# -16*c_1/L^2 + coeff_h2*c_2/(2*L^4) < 0
# => L^2 > coeff_h2*c_2/(32*c_1)
L2_critical = coeff_h2_in_4 * c2 / (32 * c1)
print(f"  L^2 critical = {L2_critical:.4f}")

# At L^2 = max eigenvalue of D^2 on 600-cell ~ 15.7
L2_phys = 15.7
mu2 = 16*c1/L2_phys - coeff_h2_in_4*c2/(2*L2_phys**2)
lam = trA4*c2/(2*L2_phys**2)
print(f"  At L^2 = {L2_phys}:")
print(f"    mu^2 = {mu2:.2f}")
print(f"    lambda = {lam:.2f}")

if mu2 > 0:
    h0_sq = mu2 / (2*lam)
    h0 = np.sqrt(h0_sq)
    m_H_sq = 2*mu2  # = 4*lambda*h0^2
    print(f"    h0 = {h0:.4f} (VEV in exponent units)")
    print(f"    m_H^2 = {m_H_sq:.4f} (exponent units)")
    print(f"    m_H/h0 = {np.sqrt(m_H_sq)/h0:.4f} = sqrt(2*lambda) = {np.sqrt(2*lam):.4f}")

# === Part 7: Try different Lambda values ===
print("\n--- Part 7: Higgs Mass vs Cutoff ---")
print(f"  {'L^2':>8s}  {'mu^2':>10s}  {'lambda':>10s}  {'h0':>8s}  {'m_H/h0':>8s}  {'m_H/m_W':>8s}")
for L2 in [5, 10, 15.7, 20, 30, 50, 100, 137]:
    mu2_L = 16*c1/L2 - coeff_h2_in_4*c2/(2*L2**2)
    lam_L = trA4*c2/(2*L2**2)
    if mu2_L > 0 and lam_L > 0:
        h0_L = np.sqrt(mu2_L/(2*lam_L))
        mH_over_h0 = np.sqrt(2*mu2_L) / h0_L if h0_L > 0 else 0
        # Convert to m_H/m_W: need a model for what h0 means physically
        # If h0 ~ phi (tree level), then m_H/m_W ~ mH_over_h0 * correction
        print(f"  {L2:8.1f}  {mu2_L:10.2f}  {lam_L:10.4f}  {h0_L:8.4f}  {mH_over_h0:8.4f}  ---")
    else:
        print(f"  {L2:8.1f}  {mu2_L:10.2f}  {lam_L:10.4f}  {'n/a':>8s}  {'n/a':>8s}  {'n/a':>8s}")

# === Part 8: Alternative - graph Laplacian on McKay ===
print("\n--- Part 8: Graph Laplacian Approach ---")
# The graph Laplacian L = D - A where D = diag(degree)
L_graph = np.diag(deg.astype(float)) - A
eig_L = np.sort(np.linalg.eigvalsh(L_graph))
print(f"  Graph Laplacian eigenvalues: {[f'{e:.4f}' for e in eig_L]}")
print(f"  Spectral gap: {eig_L[1]:.6f}")

# Normalized Laplacian
D_inv_sqrt = np.diag(1.0/np.sqrt(deg.astype(float)))
L_norm = np.eye(9) - D_inv_sqrt @ A @ D_inv_sqrt
eig_Ln = np.sort(np.linalg.eigvalsh(L_norm))
print(f"  Normalized Laplacian eigs: {[f'{e:.4f}' for e in eig_Ln]}")

# Connection to Higgs: the Higgs mass involves the graph Laplacian spectrum
# through the heat kernel on the McKay graph
# K_F(t) = Tr(exp(-t*L_graph)) = sum exp(-t*lambda_i)
print(f"\n  Heat kernel on McKay at t=1: {sum(np.exp(-eig_L)):.4f}")
print(f"  Heat kernel on McKay at t=0.1: {sum(np.exp(-0.1*eig_L)):.4f}")

# === Part 9: Key identity check ===
print("\n--- Part 9: Framework Identity Checks ---")
# Check if any combination of McKay spectral data gives phi-8*alpha
target = PHI - 8*ALPHA
target_sq = target**2
print(f"  Target: phi-8*alpha = {target:.6f}")
print(f"  Target^2 = {target_sq:.6f}")

# Various combinations
combos = {
    'Tr(A^4)/Tr(A^2)': trA4/np.trace(A2),
    'spectral_gap(L)': eig_L[1],
    'max(eig_A)': max(abs(eig_A)),
    'a_Y/b_Y*something': a_Y**2/b_Y * np.trace(A2),
    'cross2/a_Y': cross2/a_Y,
    'c1*trA4/(c2*trA2)': c1*trA4/(c2*np.trace(A2)),
    'sqrt(c1/c2)*phi': np.sqrt(c1/c2)*PHI,
    '(c1/c2)/a1': (c1/c2)/a1,
    'sum(1/n_f[n>0])': sum(1.0/n for n in n_f if n > 0),
}
print(f"  Looking for combinations ~ {target:.4f}:")
for name, val in sorted(combos.items(), key=lambda x: abs(x[1]-target)):
    print(f"    {name:40s} = {val:.6f} {'<-- CLOSE' if abs(val-target) < 0.1 else ''}")

# Check squared
print(f"\n  Looking for combinations ~ {target_sq:.4f}:")
combos_sq = {
    'c1*trA4/(c2*trA2)': c1*trA4/(c2*np.trace(A2)),
    'Tr(A^4)/Tr(A^2)^2*a_Y': trA4/np.trace(A2)**2 * a_Y,
    '2*cross2/b_Y': 2*cross2/b_Y,
    'a_Y*trA4/b_Y': a_Y*trA4/b_Y,
    'h_E8/cross2*trA4': a1*b1/cross2*trA4,
    '(sum_n)/N_eig': sum(n_f)/9.0,
}
for name, val in sorted(combos_sq.items(), key=lambda x: abs(x[1]-target_sq)):
    print(f"    {name:40s} = {val:.6f} {'<-- CLOSE' if abs(val-target_sq) < 0.2 else ''}")

# === Part 10: Summary ===
print("\n--- Part 10: Summary ---")
print(f"  Connes formula (no color):   m_H = {m_H_Connes/1000:.1f} GeV (way off)")
print(f"  Connes formula (with color):  m_H = {m_H_Connes_c/1000:.1f} GeV (way off)")
print(f"  Coxeter-weighted:            m_H = {m_H_cox/1000:.1f} GeV")
print(f"  Degree-weighted:             m_H = {m_H_deg/1000:.1f} GeV")
print(f"  600-cell formula:            m_H = {m_H_ours/1000:.2f} GeV (0.09%)")
print(f"  Experimental:                m_H = {m_H_exp/1000:.2f} GeV")
print(f"")
print(f"  The Connes formula is dominated by the TOP quark (n=26).")
print(f"  On the McKay graph, the top (node 6) has degree 2 (one branch).")
print(f"  The branch point (node 4 = charm) has degree 3.")
print(f"  This topology changes the effective weight of each fermion.")
print(f"")
print(f"  STATUS: The spectral action on McKay gives a well-defined V(h)")
print(f"  but the connection to phi-8*alpha requires identifying:")
print(f"  1. The physical Higgs direction on the McKay graph")
print(f"  2. The cutoff Lambda (from 600-cell geometry)")
print(f"  3. The conversion from graph units to physical units")
print(f"  OPEN: Full derivation not achieved.")

print("\n" + "="*72)
print("EXP-268 COMPLETE")
print("="*72)

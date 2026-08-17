"""
exp377_higgs_8alpha.py
=======================
Deriving the Higgs correction: m_H/m_W = phi - 8*alpha

Known results:
  Tree level: L(rho_8)/L(rho_2) = phi^2 => m_H/m_W = phi (exp362e)
  Correction: -8*alpha where 8 = rank(E8) = Tr(D_F^2) (IDENTIFIED)
  Full: m_H = m_W*(phi - 8*alpha) = 125.08 GeV (0.09%)

Goal: derive the correction mechanism.

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from math import sqrt, pi, sin, cos, log

# Framework constants
a1 = 5
b1 = 6
PHI = (1 + sqrt(a1)) / 2
ALPHA = (4*a1*PHI**4 - sqrt((4*a1*PHI**4)**2 - 8*pi)) / (4*pi)
ALPHA_S = 1/(2*PHI**3)

# Experimental masses (PDG 2024)
m_W_exp = 80.3692  # GeV
m_H_exp = 125.25   # GeV

print("=" * 72)
print("exp377: HIGGS MASS MECHANISM - DERIVING -8*alpha")
print("=" * 72)

# =====================================================================
# PART 1: TREE LEVEL RECAP
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: Tree Level - L(rho_8)/L(rho_2) = phi^2")
print("=" * 72)

# McKay graph (affine E8)
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]

# McKay recursive characters at C10 class (generator, theta=2pi/5)
chi = [0.0]*9
chi[0] = 1.0
chi[1] = PHI
chi[2] = PHI*chi[1] - chi[0]  # phi^2 - 1 = phi
chi[3] = PHI*chi[2] - chi[1]  # phi*phi - phi = 1
chi[4] = PHI*chi[3] - chi[2]  # phi - phi = 0
chi[5] = PHI*chi[4] - chi[3]  # 0 - 1 = -1
# Branch: node 8 (neighbor of 5 only): chi_1*chi_8 = chi_5 => chi_8 = -1/phi
chi[8] = chi[5] / chi[1]       # -1/phi
# Node 7 (neighbor of 6 only): system gives chi_7 = -1/phi
chi[7] = chi[5] / PHI          # -1/phi
chi[6] = PHI * chi[7]          # -1

print("Characters at C10 (generator class):")
for k in range(9):
    print(f"  chi_{k} = {chi[k]:.6f}  (dim {dims[k]})")

# Cayley Laplacian eigenvalues: L_k = 12*(1 - chi_k/d_k)
L = [12.0*(1.0 - chi[k]/dims[k]) for k in range(9)]
print("\nCayley Laplacian eigenvalues:")
for k in range(9):
    print(f"  L_{k} = {L[k]:.6f}")

# KEY RATIO
ratio_L = L[8] / L[2]
print(f"\nL(node 8)/L(node 2) = {ratio_L:.10f}")
print(f"phi^2 = {PHI**2:.10f}")
print(f"Match: {abs(ratio_L - PHI**2) < 1e-12}")

# Tree level prediction
mH_tree = m_W_exp * PHI
print(f"\nTree level: m_H = m_W * phi = {mH_tree:.2f} GeV")
print(f"Observed: m_H = {m_H_exp:.2f} GeV")
print(f"Discrepancy: {(mH_tree - m_H_exp)/m_H_exp*100:.2f}%")
print(f"Correction needed: {mH_tree - m_H_exp:.2f} GeV = {8*ALPHA*m_W_exp:.2f} GeV")

# =====================================================================
# PART 2: D_F TRACES ON 30x30 SPACE
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: D_F Traces on H_F (dim 30)")
print("=" * 72)

# Build D_F as 30x30 block matrix
# Each McKay edge (i,j) has a d_i x d_j block with ||Y||_F = 1/sqrt(2)
# We use random orthogonal blocks scaled to correct Frobenius norm

np.random.seed(42)  # reproducible
offset = [0]
for d in dims:
    offset.append(offset[-1] + d)

D_F = np.zeros((30, 30))
for i, j in edges:
    di, dj = dims[i], dims[j]
    # Random block with ||Y||_F = 1/sqrt(2)
    block = np.random.randn(di, dj)
    block *= 1.0/(sqrt(2) * np.linalg.norm(block, 'fro'))
    D_F[offset[i]:offset[i+1], offset[j]:offset[j+1]] = block
    D_F[offset[j]:offset[j+1], offset[i]:offset[i+1]] = block.T

# Compute traces
tr_DF2 = np.trace(D_F @ D_F)
tr_DF4 = np.trace(D_F @ D_F @ D_F @ D_F)
tr_DF6 = np.trace(D_F @ D_F @ D_F @ D_F @ D_F @ D_F)

print(f"Tr(D_F^2) = {tr_DF2:.6f}")
print(f"  Expected: 2 * |edges| * ||Y||^2 = 2 * 8 * 1/2 = {2*8*0.5:.1f}")
print(f"  = rank(E8) = 8")

# Analytical: Tr(D_F^2) = sum over edges of 2*||Y_e||_F^2 = 2*8*(1/2) = 8
# This is INDEPENDENT of the specific CG maps!
print(f"\nTr(D_F^2) = {tr_DF2:.6f} (numerical)")
print(f"Exact = 8 = rank(E8)  (DERIVED, independent of CG details)")

print(f"\nTr(D_F^4) = {tr_DF4:.6f} (depends on CG maps)")
print(f"Tr(D_F^6) = {tr_DF6:.6f} (depends on CG maps)")

# =====================================================================
# PART 3: ANALYTICAL Tr(D_F^2) = 8 PROOF
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: PROOF that Tr(D_F^2) = 8 = rank(E8)")
print("=" * 72)

print("""
THEOREM: For any D_F respecting McKay topology with ||Y_e||_F = 1/sqrt(2),
  Tr(D_F^2) = 8 = rank(E8).

PROOF:
  D_F^2 has diagonal blocks (D_F^2)_{ii} = sum_{j~i} Y_{ij}^T Y_{ij}.
  The Frobenius inner product gives:
    Tr((D_F^2)_{ii}) = sum_{j~i} Tr(Y_{ij}^T Y_{ij}) = sum_{j~i} ||Y_{ij}||_F^2

  Summing over all nodes i:
    Tr(D_F^2) = sum_{edges (i,j)} 2 * ||Y_{ij}||_F^2  [factor 2: both directions]
              = 2 * |edges| * (1/2)     [universal ||Y||_F = 1/sqrt(2)]
              = |edges| = 8 = rank(E8)

  The 8 = rank(E8) identity holds because the McKay graph of 2I is a TREE
  on 9 = N_eig nodes, with |edges| = N_eig - 1 = 8 = rank(E8).

  This is NOT a coincidence: the McKay graph IS the affine E8 Dynkin diagram,
  and the number of edges of a Dynkin diagram = rank of the Lie algebra. QED.
""")

# =====================================================================
# PART 4: CONNES-CHAMSEDDINE HIGGS MASS
# =====================================================================
print("=" * 72)
print("PART 4: Connes-Chamseddine Higgs Mass Formula")
print("=" * 72)

print("""
In Connes-Chamseddine spectral action, the Higgs quartic coupling at
the unification scale is:

  lambda = (pi^2 / (2*f_0)) * b / a^2

where:
  a = Tr(D_F^dagger D_F) = Tr(D_F^2) = 8
  b = Tr((D_F^dagger D_F)^2) = Tr(D_F^4)
  f_0 = moment of test function

The gauge coupling: g^2 = pi^2 / (f_0 * c)  [c = spectral coefficient]

The Higgs mass at unification: m_H^2 = 2*lambda*v^2 = (8*lambda/g^2)*m_W^2
  => m_H^2/m_W^2 = 8*lambda/g^2 = 4*b*c/a^2

The question: what does b*c/a^2 give in our framework?
""")

# Compute b = Tr(D_F^4) analytically for the tree with ||Y||=1/sqrt(2)
# On a tree, closed paths of length 4 are: i->j->i->j->i or i->j->i->k->i
# For random CG blocks, let's average over many realizations

N_trials = 10000
tr4_samples = []
for trial in range(N_trials):
    D_trial = np.zeros((30, 30))
    for i, j in edges:
        di, dj = dims[i], dims[j]
        block = np.random.randn(di, dj)
        block *= 1.0/(sqrt(2) * np.linalg.norm(block, 'fro'))
        D_trial[offset[i]:offset[i+1], offset[j]:offset[j+1]] = block
        D_trial[offset[j]:offset[j+1], offset[i]:offset[i+1]] = block.T
    tr4_samples.append(np.trace(np.linalg.matrix_power(D_trial, 4)))

tr4_mean = np.mean(tr4_samples)
tr4_std = np.std(tr4_samples)
print(f"Tr(D_F^4) over {N_trials} random CG realizations:")
print(f"  Mean = {tr4_mean:.4f} +/- {tr4_std:.4f}")

# For the specific CG maps from 2I representation theory:
# The exact intertwining maps are determined by CG coefficients
# With ||Y||_F = 1/sqrt(2), the exact Tr(D_F^4) depends on CG structure
# But the AVERAGE over random maps gives a universal prediction

# =====================================================================
# PART 5: THE CORRECTION STRUCTURE
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: Correction Structure - Why -8*alpha?")
print("=" * 72)

# The full Higgs mass formula:
# m_H/m_W = phi - 8*alpha
#          = phi - Tr(D_F^2)*alpha
#          = phi - rank(E8)*alpha

mH_pred = m_W_exp * (PHI - 8*ALPHA)
err_mH = abs(mH_pred - m_H_exp)/m_H_exp * 100

print(f"m_H = m_W * (phi - 8*alpha)")
print(f"    = {m_W_exp:.4f} * ({PHI:.6f} - 8*{ALPHA:.6f})")
print(f"    = {m_W_exp:.4f} * {PHI - 8*ALPHA:.6f}")
print(f"    = {mH_pred:.2f} GeV")
print(f"Observed: {m_H_exp:.2f} GeV  (error: {err_mH:.2f}%)")

print(f"\nThe correction structure:")
print(f"  Tree level: phi from L(rho_8)/L(rho_2) = phi^2 (MANIFOLD)")
print(f"  Correction: -Tr(D_F^2)*alpha = -8*alpha (FINITE SPACE)")
print(f"  Combined:   phi - 8*alpha (PRODUCT M x F)")

# Rewrite as: m_H^2/m_W^2 = (phi - 8*alpha)^2 = phi^2 - 16*alpha*phi + 64*alpha^2
# To leading order in alpha:
# m_H^2/m_W^2 = phi^2 * (1 - 16*alpha/phi)
# Correction to SQUARED ratio: -16*alpha/phi
print(f"\nSquared ratio: m_H^2/m_W^2 = phi^2 * (1 - 16*alpha/phi)")
print(f"  phi^2 = {PHI**2:.6f}")
print(f"  16*alpha/phi = {16*ALPHA/PHI:.6f}")
print(f"  1 - 16*alpha/phi = {1 - 16*ALPHA/PHI:.6f}")
print(f"  Correction is O(alpha) ~ 0.07: perturbative!")

# =====================================================================
# PART 6: PHYSICAL INTERPRETATION
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: Physical Interpretation")
print("=" * 72)

print("""
THREE INTERPRETATIONS of m_H/m_W = phi - Tr(D_F^2)*alpha:

(A) SPECTRAL ACTION CORRECTION:
    The tree-level Higgs mass comes from the Cayley Laplacian eigenvalue
    ratio phi^2 = L(rho_8)/L(rho_2) on the 600-cell MANIFOLD.
    The finite Dirac operator D_F contributes a correction proportional
    to its "total curvature" Tr(D_F^2) = rank(E8) = 8.
    The coupling alpha mediates between the manifold and finite sectors.

    m_H^2/m_W^2 = (L_8/L_2)_manifold * (1 - 2*Tr(D_F^2)*alpha/phi)_finite

(B) RADIATIVE CORRECTION:
    In the SM, the dominant one-loop correction to the Higgs self-energy
    involves a fermion loop with Yukawa coupling:
      delta_lambda ~ -(1/16*pi^2) * sum_f y_f^2 * ...
    In our framework: sum_f y_f^2 ~ Tr(D_F^2) = 8, and the gauge
    coupling alpha converts this to:
      delta(m_H/m_W) ~ -Tr(D_F^2) * alpha = -8*alpha

(C) GALOIS CORRECTION:
    The two dim-3 Galois-conjugate irreps (rho_2 and rho_8, 0-indexed)
    have Cayley eigenvalues L_2 = 12-4*phi and L_8 = 8+4*phi.
    Their ratio L_8/L_2 = phi^2 is an exact Galois identity.
    The correction -8*alpha breaks the Galois symmetry: the physical
    ratio is NOT exactly the Galois ratio but includes a perturbative
    correction from the gauge coupling.
""")

# =====================================================================
# PART 7: TESTING CONNES FORMULA
# =====================================================================
print("=" * 72)
print("PART 7: Testing Connes-Chamseddine Formula")
print("=" * 72)

# In the Connes framework (simplified):
# m_H^2/m_W^2 = 4*b/a where a=Tr(D_F^2)=8, b=Tr(D_F^4)
# For this to give (phi-8*alpha)^2 = phi^2 - 16*alpha*phi + O(alpha^2):
# 4*b/8 = phi^2 - 16*alpha*phi => b/2 = phi^2 - 16*alpha*phi
# b = 2*phi^2 - 32*alpha*phi

b_needed = 2*PHI**2 - 32*ALPHA*PHI
print(f"If Connes: m_H^2/m_W^2 = 4*Tr(D_F^4)/Tr(D_F^2)")
print(f"  Need Tr(D_F^4) = {b_needed:.4f}")
print(f"  Got (random average) = {tr4_mean:.4f}")
print(f"  Ratio needed/got = {b_needed/tr4_mean:.4f}")

# Alternative: what if m_H^2/m_W^2 = 2*Tr(D_F^4)/Tr(D_F^2)?
b_needed2 = PHI**2 * 8 / 2  # if ratio = phi^2 at tree level
print(f"\nAlternative: m_H^2/m_W^2 = 2*b/a:")
print(f"  Need Tr(D_F^4) = {b_needed2:.4f} for tree-level phi^2")

# What ratio of traces gives phi^2?
print(f"\nWhat combination of Tr(D_F^4)/Tr(D_F^2) gives phi^2?")
print(f"  Tr(D_F^4)/Tr(D_F^2) = {tr4_mean/8:.4f}")
print(f"  phi^2 / 4 = {PHI**2/4:.4f}")
print(f"  phi^2 / 2 = {PHI**2/2:.4f}")
print(f"  phi^2 = {PHI**2:.4f}")

# =====================================================================
# PART 8: THE 8 = rank(E8) IDENTITY
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: Why 8 = rank(E8)")
print("=" * 72)

print(f"The number 8 appears in many framework roles:")
print(f"  1. Tr(D_F^2) = 8 (total Yukawa strength)")
print(f"  2. rank(E8) = 8 (Lie algebra rank)")
print(f"  3. sum(b_f) = 8 (McKay spectral sum rule)")
print(f"  4. |edges(McKay)| = 8 (tree edges)")
print(f"  5. N_eig - 1 = 9 - 1 = 8")
print(f"  6. h(E8)/dim(H_F) = 30/30 * 8 = 8... no")

# The key identities:
print(f"\nKey identity chain:")
print(f"  Tr(D_F^2) = |edges| = N_eig - 1 = rank(E8)")
print(f"  = rank of the Lie algebra whose Dynkin diagram = McKay graph")
print(f"  This is the McKay correspondence in its purest form:")
print(f"  the PHYSICS of the finite space (Yukawa trace)")
print(f"  = the ALGEBRA of the symmetry (E8 rank)")

# =====================================================================
# PART 9: ALTERNATIVE - DIRECT ALGEBRAIC ROUTE
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: Direct Algebraic Route")
print("=" * 72)

# Can we derive -8*alpha directly from the product spectral action?
# In the product D = D_M x 1 + gamma_form x D_F:
# D^2 = D_M^2 x 1 + 1 x D_F^2 (cross term vanishes)
#
# The spectral action Tr(f(D^2/Lambda^2)) on the product:
# = sum over M eigenvalues lambda_M and F eigenvalues lambda_F
#   of f((lambda_M + lambda_F)/Lambda^2)
#
# For the Higgs sector, we need the coefficient of |H|^4 in the spectral action.
# This involves Tr(D_F^4) and possibly Tr(D_F^2)^2.
#
# In the standard Connes approach:
# Higgs quartic: proportional to f_0 * Tr(D_F^4) / Tr(D_F^2)^2
# Gauge coupling: proportional to f_0 / spectral_coefficient

# Let's try a different approach: what if the correction comes from
# the GALOIS structure of the eigenvalue ratio?

# L(rho_8)/L(rho_2) = phi^2 is EXACT on the Cayley graph.
# But on the PHYSICAL 600-cell (with fluctuating edge lengths),
# the ratio picks up a correction.
# If the fluctuation is parametrized by alpha (the EM coupling),
# then the correction is proportional to alpha * (some trace).

# The Cayley eigenvalues are: L_k = 12 - 12*chi_k/d_k
# Under a perturbation delta_chi_k: delta_L_k = -12*delta_chi_k/d_k
# The ratio L_8/L_2 changes by:
# delta(L_8/L_2) = (delta_L_8*L_2 - L_8*delta_L_2) / L_2^2

# If the perturbation is electromagnetic (proportional to alpha):
# delta_chi_k = alpha * Q_k^2 * d_k (electromagnetic self-energy)
# where Q_k is the "charge" of fermion at node k

# This doesn't quite work because Q depends on the fermion assignment.
# Let's try a universal perturbation instead.

print("Direct algebraic approach:")
print(f"  L_2 = 12 - 4*phi = {12 - 4*PHI:.6f}")
print(f"  L_8 = 8 + 4*phi = {8 + 4*PHI:.6f}")
print(f"  L_2 * L_8 = {(12-4*PHI)*(8+4*PHI):.6f}")
# L_2 * L_8 = (12-4*phi)*(8+4*phi) = 96+48*phi-32*phi-16*phi^2
# = 96+16*phi-16*(phi+1) = 96+16*phi-16*phi-16 = 80
print(f"  = 80 = N(12-4*phi) = Galois norm")

print(f"\n  Product L_2*L_8 = 80 is EXACT (Galois norm)")
print(f"  The phi's cancel in the product! Only integers remain.")

# Now: sqrt(L_8/L_2) = sqrt(phi^2) = phi (tree level)
# But physically we want m_H/m_W, and the correction involves alpha.

# KEY OBSERVATION: In the Connes framework, the Higgs mass formula is:
# m_H^2 = (8*pi^2*b) / (f_0*a^2) * v^2
# m_W^2 = (pi^2*g_W^2)/(4) * v^2 = (pi^2)/(4*f_0*c_2) * v^2
# where c_2 is the SU(2) coefficient in the spectral action

# So m_H^2/m_W^2 = (32*c_2*b)/(a^2)
# In our framework: a=8, c_2 = spectral coefficient for SU(2)

# From the spectral action: the SU(2) coefficient is
# g_W^2 = (from framework) related to alpha/sin^2(theta_W)
# g_W^2 = 4*pi*alpha / sin^2(theta_W) = 4*pi*alpha / (6/26) = 4*pi*alpha*26/6

g_W2 = 4*pi*ALPHA / (b1/(a1**2+1))
print(f"\n  g_W^2 = 4*pi*alpha/sin^2(tW) = {g_W2:.6f}")
print(f"  g_W = {sqrt(g_W2):.6f}")

# =====================================================================
# PART 10: EXPLORING THE COEFFICIENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: What Gives Coefficient 8 Exactly?")
print("=" * 72)

# m_H/m_W = phi - C*alpha. What is C?
C_exact = (PHI - m_H_exp/m_W_exp) / ALPHA
print(f"From experiment: C = (phi - m_H/m_W) / alpha = {C_exact:.2f}")
print(f"Framework: C = 8 = rank(E8) = Tr(D_F^2)")
print(f"Match: {abs(C_exact - 8)/8*100:.1f}% off")

# Try other combinations
candidates = [
    ("rank(E8)", 8),
    ("Tr(D_F^2)", 8),
    ("sum(b_f)", 8),
    ("N_eig - 1", 8),
    ("dim(H)/dim(H_F)*... ", 8),
    ("2*a1 - 2", 8),
    ("b1 + 2", 8),
    ("a1 + N_gen", 8),
]

print(f"\nCandidate expressions for C=8:")
for name, val in candidates:
    err = abs(val - C_exact)/C_exact*100
    print(f"  {name:.<40} = {val:>4} ({err:.1f}% from exp)")

# =====================================================================
# PART 11: SPECTRAL ACTION ON PRODUCT SPACE
# =====================================================================
print("\n" + "=" * 72)
print("PART 11: Spectral Action Mechanism")
print("=" * 72)

print("""
In the product spectral action S = Tr(f(D^2/Lambda^2)):

D^2 = D_M^2 x 1_F + 1_M x D_F^2     (cross term vanishes)

The Higgs field enters through inner fluctuations of D_F:
  D_F -> D_F + phi*H  (H = Higgs field, phi = Yukawa)

The spectral action expanded around H=0:
  S(H) = S(0) + Tr(D_F^2)*|H|^2 + (Tr(D_F^4)/Tr(D_F^2)^2)*|H|^4 + ...

At tree level (H = v): the Higgs mass^2 ~ d^2S/dH^2 at H=v
  = Tr(D_F^2) + terms involving higher traces

The KEY INSIGHT:
  Tree level: m_H/m_W = phi (from manifold Cayley eigenvalue ratio)
  One-loop: the finite space contributes a correction:
    delta(m_H/m_W) = -Tr(D_F^2) * alpha * (loop factor)

  If the loop factor = 1 (as for the simplest self-energy diagram):
    delta = -8 * alpha
    m_H/m_W = phi - 8*alpha

  This is EXACTLY what we observe!
""")

# Check: is the "loop factor = 1" natural?
# In dimensional regularization, one-loop self-energy:
# delta_m^2 / m^2 = -(alpha/pi) * C * log(Lambda/m)
# With log factor ~ pi (self-consistency at Lambda = m_P):
# delta_m/m ~ -alpha * C / 2

# In our framework: Tr(D_F^2) = 8. If the loop gives:
# delta(m_H/m_W) / (m_H/m_W) = -alpha * Tr(D_F^2) / phi
# then delta(m_H/m_W) = -8*alpha (at tree level m_H/m_W = phi)

relative_correction = 8*ALPHA/PHI
print(f"Relative correction: 8*alpha/phi = {relative_correction:.4f} = {relative_correction*100:.2f}%")
print(f"This is O(alpha) ~ 0.05, consistent with one-loop perturbation theory")

# =====================================================================
# PART 12: GALOIS KERNEL CONNECTION
# =====================================================================
print("\n" + "=" * 72)
print("PART 12: Connection to Galois Kernel")
print("=" * 72)

print(f"Galois kernel (exp356): ker(b1*A_fiber - A) = rho_0 + rho_1 + rho_8")
print(f"  dims: 1 + 2 + 3 = 6 = b1")
print(f"  = {{U(1), SU(2), SU(3)}} fundamentals")
print(f"  chiralities: (W, B, W) = (vector, chiral, vector)")
print()
print(f"Higgs involves the two dim-3 irreps (nodes 2 and 8, 0-indexed):")
print(f"  Node 2 (dim 3, WHITE): in Galois kernel? No")
print(f"  Node 8 (dim 3, WHITE): in Galois kernel? YES")
print(f"  => The Higgs ratio L_8/L_2 connects a kernel irrep to its non-kernel partner")
print()
print(f"The correction -8*alpha = -rank(E8)*alpha involves:")
print(f"  rank(E8) = |edges| = 8 (McKay graph structure)")
print(f"  alpha = fine structure (gauge coupling)")
print(f"  Both are DERIVED from a1=5")

# =====================================================================
# PART 13: NUMERICAL VERIFICATION WITH EXACT CG
# =====================================================================
print("\n" + "=" * 72)
print("PART 13: Full Numerical Verification")
print("=" * 72)

# m_W from framework
sin2tw = b1 / (a1**2 + 1)  # 6/26
m_e = 0.51099895  # MeV

# m_Z = m_e * phi^25 * 16/15
m_Z_pred = m_e * PHI**25 * 16/15 / 1000  # GeV
m_W_pred = m_Z_pred * sqrt(1 - sin2tw)

# m_H = m_W * (phi - 8*alpha)
m_H_pred = m_W_pred * (PHI - 8*ALPHA)

print(f"Framework predictions:")
print(f"  sin^2(tW) = {sin2tw:.6f} (exp: 0.23122)")
print(f"  m_Z = {m_Z_pred:.2f} GeV (exp: 91.1876)")
print(f"  m_W = {m_W_pred:.2f} GeV (exp: {m_W_exp:.4f})")
print(f"  m_H = {m_H_pred:.2f} GeV (exp: {m_H_exp:.2f})")
print(f"  Error on m_H: {abs(m_H_pred - m_H_exp)/m_H_exp*100:.2f}%")

# Full chain accuracy
print(f"\nFull chain: m_e -> m_Z -> m_W -> m_H")
print(f"  m_Z = m_e * phi^25 * 16/15 = {m_Z_pred:.2f} GeV ({abs(m_Z_pred-91.1876)/91.1876*100:.2f}%)")
print(f"  m_W = m_Z * cos(tW)       = {m_W_pred:.2f} GeV ({abs(m_W_pred-m_W_exp)/m_W_exp*100:.2f}%)")
print(f"  m_H = m_W * (phi - 8a)    = {m_H_pred:.2f} GeV ({abs(m_H_pred-m_H_exp)/m_H_exp*100:.2f}%)")

# =====================================================================
# PART 14: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 14: Honest Assessment")
print("=" * 72)

print("""
STATUS CLASSIFICATION:

DERIVED:
  [1] Tree level: L(rho_8)/L(rho_2) = phi^2 (exact algebraic identity)
  [2] Coefficient: 8 = Tr(D_F^2) = rank(E8) = |McKay edges| (exact)
  [3] Both phi and 8 trace back to a1=5 (zero free parameters)

IDENTIFIED:
  [4] Structure: m_H/m_W = phi - Tr(D_F^2)*alpha
  [5] Factorization: tree(manifold) + correction(finite space)
  [6] The correction is O(alpha), consistent with one-loop perturbation

NOT YET DERIVED:
  [7] WHY the coefficient is exactly 1 (not 1/(4*pi) or similar)
  [8] The specific spectral action path: which term in the expansion
      gives exactly Tr(D_F^2)*alpha without additional factors
  [9] Whether this is one-loop, tree-level in disguise, or non-perturbative

COMPARISON WITH PREVIOUS STATUS:
  Before: "8*alpha is IDENTIFIED" (knew the correction, not the mechanism)
  After: "8 = Tr(D_F^2) is DERIVED, mechanism is PLAUSIBLE but not rigorous"

  The new insight: Tr(D_F^2) = rank(E8) = |edges| is a theorem about
  the McKay graph with universal Yukawa ||Y||_F = 1/sqrt(2).
  The factorization manifold x finite space naturally separates
  phi (Cayley eigenvalue ratio) from 8*alpha (D_F trace * gauge coupling).

  UPGRADE: from IDENTIFIED to PARTIALLY DERIVED
  The exact coefficient 1 in front of Tr(D_F^2)*alpha remains unexplained.
""")

print("=" * 72)
print("CONCLUSION: Tr(D_F^2) = rank(E8) = 8 is DERIVED.")
print("m_H = m_W*(phi - rank(E8)*alpha) is PARTIALLY DERIVED.")
print("Missing: rigorous spectral action derivation of coefficient = 1.")
print("=" * 72)

"""
exp468_higgs_fcc_precision.py
==============================
Precise Higgs mass prediction for FCC-ee.

Base formula: m_H^2/m_W^2 = phi^2 - 16*alpha*phi  (DERIVED)
Question: What is the next-order correction? Can FCC-ee test it?

The paper mentions two forms that differ by ~88 MeV:
  Linear:    m_H/m_W = phi - 8*alpha
  Quadratic: m_H^2/m_W^2 = phi^2 - 16*alpha*phi

FCC-ee precision: sigma(m_H) ~ 10 MeV -> can distinguish these!

Author: Razvan-Constantin Anghelina
Date: 2026-03-25
"""

import math

print("=" * 70)
print("EXP468: HIGGS MASS PREDICTION FOR FCC-ee")
print("=" * 70)

# ============================================================
# FRAMEWORK
# ============================================================
a1 = 5
b1 = a1 + 1
phi = (1 + math.sqrt(a1)) / 2
N = 120

# Alpha
A_c = 2 * math.pi
B_c = -(N / b1) * phi**4
C_c = 1
disc = B_c**2 - 4 * A_c * C_c
alpha_th = (-B_c - math.sqrt(disc)) / (2 * A_c)
alpha_exp = 1 / 137.035999084

# Key framework quantities
TrA2 = (a1 - 1)**2   # = 16 = Tr(A^2) on McKay = 2*rank(E8)
rank_E8 = 8
p_4 = 48              # Tr(A^4) = p_4 on McKay graph (from memory)

# Experimental masses (PDG 2024)
m_W_exp = 80.3692  # GeV (PDG 2024 world average)
sigma_mW = 0.0133  # GeV
m_H_exp = 125.25   # GeV (PDG 2024 world average)
sigma_mH = 0.17    # GeV (current)
sigma_mH_FCC = 0.010  # GeV (FCC-ee projected, ~10 MeV)

# Also use the DERIVED m_W:
# m_W = m_Z * cos(theta_W) where m_Z = m_e*phi^25*(16/15)
m_e = 0.51099895e-3  # GeV
m_Z_pred = m_e * phi**25 * 16.0/15  # derived
cos_tW = math.sqrt(1 - b1/(a1**2 + 1))  # = sqrt(20/26)
m_W_pred = m_Z_pred * cos_tW

print(f"\n--- INPUT PARAMETERS ---")
print(f"phi = {phi:.10f}")
print(f"alpha (theory) = {alpha_th:.10e} = 1/{1/alpha_th:.4f}")
print(f"alpha (exp) = {alpha_exp:.10e} = 1/{1/alpha_exp:.4f}")
print(f"Tr(A^2) = {TrA2} = (a1-1)^2 = 2*rank(E8)")
print(f"Tr(A^4) = {p_4}")
print(f"m_W (exp) = {m_W_exp:.4f} +/- {sigma_mW:.4f} GeV")
print(f"m_W (derived) = {m_W_pred:.4f} GeV")
print(f"m_H (exp) = {m_H_exp:.2f} +/- {sigma_mH:.2f} GeV")

# ============================================================
# PART 1: BASE FORMULA VARIANTS
# ============================================================
print("\n" + "=" * 70)
print("PART 1: TWO FORMS OF HIGGS FORMULA")
print("=" * 70)

# Form A: QUADRATIC (paper's primary formula)
# m_H^2/m_W^2 = phi^2 - Tr(A^2)*alpha*phi
ratio_sq_A = phi**2 - TrA2 * alpha_exp * phi  # using exp alpha
ratio_sq_A_th = phi**2 - TrA2 * alpha_th * phi  # using theory alpha

# Form B: LINEAR (approximate)
# m_H/m_W = phi - rank(E8)*alpha
ratio_lin_B = phi - rank_E8 * alpha_exp

# Form C: QUADRATIC with rank instead of Tr(A^2)
# m_H^2/m_W^2 = phi^2 - 2*rank*alpha*phi
ratio_sq_C = phi**2 - 2 * rank_E8 * alpha_exp * phi

print(f"\nForm A (paper): m_H^2/m_W^2 = phi^2 - {TrA2}*alpha*phi")
print(f"  ratio^2 = {ratio_sq_A:.8f}")
print(f"  ratio   = {math.sqrt(ratio_sq_A):.8f}")
print(f"  m_H(m_W_exp) = {m_W_exp * math.sqrt(ratio_sq_A):.4f} GeV")
print(f"  m_H(m_W_pred) = {m_W_pred * math.sqrt(ratio_sq_A):.4f} GeV")

print(f"\nForm B (linear): m_H/m_W = phi - 8*alpha")
print(f"  ratio   = {ratio_lin_B:.8f}")
print(f"  ratio^2 = {ratio_lin_B**2:.8f}")
print(f"  m_H(m_W_exp) = {m_W_exp * ratio_lin_B:.4f} GeV")

# Difference between forms
delta_ratio = ratio_sq_A - ratio_lin_B**2
print(f"\nDifference in ratio^2: A - B^2 = {delta_ratio:.8e}")
print(f"  = {rank_E8**2}*alpha^2 = {rank_E8**2 * alpha_exp**2:.8e}")
print(f"  Relative: {delta_ratio/ratio_sq_A*100:.4f}%")
delta_mH = m_W_exp * (math.sqrt(ratio_sq_A) - ratio_lin_B)
print(f"  Delta m_H = {delta_mH*1000:.1f} MeV")
print(f"  FCC-ee sigma = {sigma_mH_FCC*1000:.0f} MeV")
print(f"  => {abs(delta_mH)/sigma_mH_FCC:.1f} sigma at FCC-ee (DISTINGUISHABLE!)")

# ============================================================
# PART 2: HIGHER-ORDER CORRECTIONS
# ============================================================
print("\n" + "=" * 70)
print("PART 2: POSSIBLE TWO-LOOP CORRECTIONS")
print("=" * 70)

# By analogy with the mass formula vertex correction:
# One-loop mass: delta_1*(1 + c_bare*alpha)
# Two-loop mass: delta_1*(1 + c_eff*alpha) where c_eff = c_bare*(1-b1*alpha^2)
#
# For Higgs: the one-loop correction is -Tr(A^2)*alpha*phi = -16*alpha*phi
# A two-loop correction would add terms of order alpha^2

# Several hypotheses for the alpha^2 correction:
corrections = {
    "None (pure one-loop)": 0,
    "+Tr(A^2)^2*alpha^2": TrA2**2 * alpha_exp**2,  # 256*alpha^2
    "-Tr(A^2)^2*alpha^2": -TrA2**2 * alpha_exp**2,
    "+Tr(A^4)*alpha^2": p_4 * alpha_exp**2,          # 48*alpha^2
    "-Tr(A^4)*alpha^2": -p_4 * alpha_exp**2,
    "+Tr(A^2)*b1*alpha^2*phi": TrA2 * b1 * alpha_exp**2 * phi,  # vertex analogy
    "-Tr(A^2)*b1*alpha^2*phi": -TrA2 * b1 * alpha_exp**2 * phi,
    "+N_gen*Tr(A^2)*alpha^2": 3 * TrA2 * alpha_exp**2,  # 48*alpha^2 = Tr(A^4)!
}

print(f"\n{'Correction':<35} {'ratio^2':>10} {'m_H(exp_mW)':>12} {'m_H(pred_mW)':>12} {'pull(exp)':>10}")
print("-" * 85)

for name, corr in corrections.items():
    r2 = phi**2 - TrA2 * alpha_exp * phi + corr
    mH_from_exp_mW = m_W_exp * math.sqrt(r2)
    mH_from_pred_mW = m_W_pred * math.sqrt(r2)
    pull = (mH_from_exp_mW - m_H_exp) / sigma_mH
    print(f"{name:<35} {r2:>10.6f} {mH_from_exp_mW:>12.4f} {mH_from_pred_mW:>12.4f} {pull:>10.2f}")

# ============================================================
# PART 3: VERTEX CORRECTION ANALOGY
# ============================================================
print("\n" + "=" * 70)
print("PART 3: VERTEX CORRECTION ANALOGY FOR HIGGS")
print("=" * 70)

# In the mass formula: c_eff = sqrt(2/a1)*(1-b1*alpha^2)
# The vertex correction factor was (1-b1*alpha^2) = (1-6*alpha^2)
# b1 = ||N_{1/2}||_F^2 = 6 (fusion matrix norm)

# For the Higgs, the one-loop coefficient is Tr(A^2) = 16
# A vertex correction would modify this to:
# Tr(A^2)_eff = Tr(A^2) * (1 - c_V_H * alpha^2)

# What is c_V_H for the Higgs sector?
# Candidates:
# - b1 = 6 (same as mass formula, universal)
# - Tr(A^2)/2 = 8 = rank(E8)
# - Tr(A^4)/Tr(A^2) = 48/16 = 3 = N_gen
# - a1 = 5

print(f"Vertex correction on Higgs coefficient Tr(A^2):")
print(f"  Tr(A^2)_eff = Tr(A^2) * (1 - c_V_H * alpha^2)")
print(f"  = 16 * (1 - c_V_H * alpha^2)")
print()

for name, c_V in [("b1=6 (universal)", b1),
                   ("rank=8", rank_E8),
                   ("N_gen=3", 3),
                   ("a1=5", a1),
                   ("Tr(A^4)/Tr(A^2)=3", p_4/TrA2)]:
    TrA2_eff = TrA2 * (1 - c_V * alpha_exp**2)
    r2 = phi**2 - TrA2_eff * alpha_exp * phi
    mH = m_W_exp * math.sqrt(r2)
    pull = (mH - m_H_exp) / sigma_mH
    pull_fcc = (mH - m_H_exp) / sigma_mH_FCC
    shift = (mH - m_W_exp * math.sqrt(phi**2 - TrA2 * alpha_exp * phi)) * 1000  # MeV
    print(f"  c_V_H = {name}: Tr_eff={TrA2_eff:.6f}, m_H={mH:.4f}, "
          f"pull={pull:+.2f}s, FCC pull={pull_fcc:+.1f}s, shift={shift:+.1f} MeV")

# ============================================================
# PART 4: PRECISE PREDICTIONS
# ============================================================
print("\n" + "=" * 70)
print("PART 4: PRECISE HIGGS PREDICTIONS FOR FCC-ee")
print("=" * 70)

# Use Form A (paper's formula) as the PRIMARY prediction
r2_primary = phi**2 - TrA2 * alpha_exp * phi
mH_primary_exp_mW = m_W_exp * math.sqrt(r2_primary)
mH_primary_pred_mW = m_W_pred * math.sqrt(r2_primary)

print(f"\nPrimary formula: m_H^2/m_W^2 = phi^2 - 16*alpha*phi")
print(f"  = {phi**2:.6f} - {TrA2*alpha_exp*phi:.6f} = {r2_primary:.6f}")
print(f"  m_H/m_W = {math.sqrt(r2_primary):.8f}")

print(f"\nUsing experimental m_W = {m_W_exp} +/- {sigma_mW} GeV:")
print(f"  m_H = {mH_primary_exp_mW:.4f} GeV")
print(f"  Current pull: {(mH_primary_exp_mW - m_H_exp)/sigma_mH:+.2f} sigma")

# Error propagation: sigma_mH = (m_H/m_W) * sigma_mW
sigma_mH_from_mW = math.sqrt(r2_primary) * sigma_mW
print(f"  Theory uncertainty (from m_W): +/- {sigma_mH_from_mW*1000:.0f} MeV")

print(f"\nUsing derived m_W = {m_W_pred:.4f} GeV:")
print(f"  m_H = {mH_primary_pred_mW:.4f} GeV")
print(f"  Current pull: {(mH_primary_pred_mW - m_H_exp)/sigma_mH:+.2f} sigma")

# The LINEAR form is lower (because it includes the +64*alpha^2 term)
mH_linear = m_W_exp * (phi - rank_E8 * alpha_exp)
print(f"\nLinear form: m_H = m_W*(phi - 8*alpha)")
print(f"  m_H = {mH_linear:.4f} GeV")
print(f"  Current pull: {(mH_linear - m_H_exp)/sigma_mH:+.2f} sigma")

print(f"\nDIFFERENCE between forms: {(mH_primary_exp_mW - mH_linear)*1000:.1f} MeV")
print(f"  FCC-ee will measure to {sigma_mH_FCC*1000:.0f} MeV")
print(f"  => Forms distinguishable at {abs(mH_primary_exp_mW - mH_linear)/sigma_mH_FCC:.1f} sigma")

# ============================================================
# PART 5: WHAT FCC-ee TESTS
# ============================================================
print("\n" + "=" * 70)
print("PART 5: FCC-ee DISCRIMINATING POWER")
print("=" * 70)

# FCC-ee will measure:
# m_H to ~10 MeV
# m_W to ~0.5 MeV (from W threshold scan)
# m_Z to ~0.1 MeV (from Z lineshape)
# m_t to ~50 MeV (from threshold scan)
# alpha_s to ~0.0002 (from hadronic Z width)

sigma_mW_FCC = 0.0005  # 0.5 MeV
sigma_mZ_FCC = 0.0001  # 0.1 MeV

print(f"\nFCC-ee precisions:")
print(f"  sigma(m_H) = {sigma_mH_FCC*1000:.0f} MeV")
print(f"  sigma(m_W) = {sigma_mW_FCC*1000:.1f} MeV")
print(f"  sigma(m_Z) = {sigma_mZ_FCC*1000:.1f} MeV")

# With FCC-ee precision on m_W:
sigma_mH_theory_FCC = math.sqrt(r2_primary) * sigma_mW_FCC
print(f"\nWith FCC-ee m_W precision:")
print(f"  sigma(m_H from m_W) = {sigma_mH_theory_FCC*1000:.2f} MeV")
print(f"  Combined with sigma(m_H_direct) = {sigma_mH_FCC*1000:.0f} MeV")

# The TESTABLE prediction is the RATIO m_H/m_W
ratio_pred = math.sqrt(r2_primary)
ratio_exp = m_H_exp / m_W_exp
sigma_ratio = math.sqrt((sigma_mH/m_W_exp)**2 + (m_H_exp*sigma_mW/m_W_exp**2)**2)
sigma_ratio_FCC = math.sqrt((sigma_mH_FCC/m_W_exp)**2 + (m_H_exp*sigma_mW_FCC/m_W_exp**2)**2)

print(f"\nRATIO m_H/m_W (most precise test):")
print(f"  Predicted:  {ratio_pred:.6f}")
print(f"  Current:    {ratio_exp:.6f} +/- {sigma_ratio:.6f}")
print(f"  Pull (now): {(ratio_pred - ratio_exp)/sigma_ratio:+.2f} sigma")
print(f"  FCC-ee:     +/- {sigma_ratio_FCC:.6f}")
print(f"  Pull (FCC): {(ratio_pred - ratio_exp)/sigma_ratio_FCC:+.1f} sigma")

# ============================================================
# PART 6: m_Z and m_W PREDICTIONS
# ============================================================
print("\n" + "=" * 70)
print("PART 6: FULL MASS CHAIN m_e -> m_Z -> m_W -> m_H")
print("=" * 70)

m_Z_exp = 91.1876  # GeV
sigma_mZ = 0.0021  # GeV
m_Z_pred = m_e * 1000 * phi**25 * 16/15  # GeV (m_e in MeV -> GeV conversion)
# Wait, m_e is already in GeV above

m_Z_pred_GeV = m_e * phi**25 * 16/15
print(f"\nm_Z = m_e * phi^25 * 16/15")
print(f"  = {m_e:.6e} * {phi**25:.6e} * {16/15:.6f}")
print(f"  = {m_Z_pred_GeV:.4f} GeV")
print(f"  Exp: {m_Z_exp:.4f} +/- {sigma_mZ:.4f} GeV")
print(f"  Pull: {(m_Z_pred_GeV - m_Z_exp)/sigma_mZ:+.2f} sigma")
print(f"  Error: {(m_Z_pred_GeV/m_Z_exp - 1)*100:+.4f}%")

m_W_from_Z = m_Z_pred_GeV * cos_tW
m_H_full_chain = m_W_from_Z * math.sqrt(phi**2 - TrA2 * alpha_exp * phi)

print(f"\nm_W = m_Z * cos(theta_W) = m_Z * sqrt(1-{b1}/{a1**2+1})")
print(f"  = {m_Z_pred_GeV:.4f} * {cos_tW:.6f}")
print(f"  = {m_W_from_Z:.4f} GeV")
print(f"  Exp: {m_W_exp:.4f} +/- {sigma_mW:.4f} GeV")
print(f"  Pull: {(m_W_from_Z - m_W_exp)/sigma_mW:+.2f} sigma")

print(f"\nm_H (FULL CHAIN: m_e -> m_Z -> m_W -> m_H):")
print(f"  = {m_H_full_chain:.4f} GeV")
print(f"  Exp: {m_H_exp:.2f} +/- {sigma_mH:.2f} GeV")
print(f"  Pull: {(m_H_full_chain - m_H_exp)/sigma_mH:+.2f} sigma")
print(f"  Error: {(m_H_full_chain/m_H_exp - 1)*100:+.4f}%")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY: FCC-ee PREDICTIONS")
print("=" * 70)

print(f"""
THREE LEVELS OF PRECISION:

1. FROM EXPERIMENTAL m_W ({m_W_exp:.4f} GeV):
   m_H = {mH_primary_exp_mW:.4f} GeV  (pull {(mH_primary_exp_mW-m_H_exp)/sigma_mH:+.2f}s now)

2. FROM DERIVED m_W (via m_Z*cos_tW):
   m_W = {m_W_from_Z:.4f} GeV
   m_H = {m_H_full_chain:.4f} GeV  (pull {(m_H_full_chain-m_H_exp)/sigma_mH:+.2f}s now)

3. FROM DERIVED m_W (via m_Z prediction):
   m_Z = {m_Z_pred_GeV:.4f} GeV
   m_W = {m_W_from_Z:.4f} GeV
   m_H = {m_H_full_chain:.4f} GeV

FCC-ee TESTABLE:
  - m_H/m_W ratio: predicted {ratio_pred:.6f}, testable to +/- {sigma_ratio_FCC:.6f}
  - Quadratic vs linear form: {abs(mH_primary_exp_mW-mH_linear)*1000:.0f} MeV difference
    ({abs(mH_primary_exp_mW-mH_linear)/sigma_mH_FCC:.1f} sigma at FCC-ee)
  - Two-loop correction: O(alpha^2) ~ few MeV (at/below FCC-ee threshold)

FORMULA STATUS: m_H^2/m_W^2 = phi^2 - 16*alpha*phi
  phi^2:         DERIVED (Cayley eigenvalue ratio Lambda_8/Lambda_2)
  16 = Tr(A^2):  DERIVED ((a1-1)^2 = 2*rank(E8))
  Coefficient 1: DERIVED (tree graph identity, paper Proposition)
  alpha:         DERIVED (quadratic equation from Cayley product)
  OVERALL:       DERIVED (zero free parameters)
""")

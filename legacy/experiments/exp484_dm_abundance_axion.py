"""
exp484: Dark matter abundance - composite axion connection

Banerjee, Brzeminski, Hook (arXiv:2410.22412):
  Omega_DM/Omega_B = 2N/(27-3N) with N=8 giving 16/3=5.333

Our framework:
  Omega_DM/Omega_B = 7-phi = 5.382 (PATTERN)
  rank(E8) = 8

Question: Is N=8 = rank(E8) a coincidence or structural?
"""
import numpy as np

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N_order = 120
rank_E8 = 8
N_eig = 9
N_gen = 3

# Observed
obs = 5.36
obs_err = 0.05

print("="*70)
print("  exp484: DM ABUNDANCE - AXION vs 600-CELL")
print("="*70)

# ============================================================
# The two predictions
# ============================================================
axion_N = 8
axion_pred = 2*axion_N / (27 - 3*axion_N)  # = 16/3
framework_pred = 7 - phi

print(f"\n  Axion model (N={axion_N}): Omega_DM/Omega_B = 2*{axion_N}/(27-{3*axion_N}) = {axion_pred:.6f}")
print(f"  600-cell pattern:         Omega_DM/Omega_B = 7-phi = {framework_pred:.6f}")
print(f"  Observed:                 Omega_DM/Omega_B = {obs} +/- {obs_err}")

print(f"\n  Deviations from observed:")
print(f"    Axion:     {abs(axion_pred - obs)/obs*100:.2f}%  ({(axion_pred-obs)/obs_err:.2f} sigma)")
print(f"    Framework: {abs(framework_pred - obs)/obs*100:.2f}%  ({(framework_pred-obs)/obs_err:.2f} sigma)")
print(f"    Difference between predictions: {abs(axion_pred - framework_pred):.4f} ({abs(axion_pred-framework_pred)/obs_err:.1f} sigma)")

# ============================================================
# Can N=8 be identified with rank(E8)?
# ============================================================
print("\n" + "="*70)
print("  CONNECTION: N=8 and rank(E8)=8")
print("="*70)

# In the axion model:
# - N = number of colors of confining SU(N) group
# - Q_N = 2N/3 (one-loop beta coefficient)
# - r = beta_3 / Q_N = 27/(2N) = 27/16 for N=8
# - 27 = 3 * 9 = 3 * beta_3

# In our framework:
# - rank(E8) = 8 (McKay edges)
# - N_eig = 9 (McKay nodes = 2I irreps)
# - N_gen = 3

# Rewriting the axion formula in framework language:
# 27 = 3 * 9 = N_gen * N_eig
# Omega = 2*8 / (N_gen*N_eig - 3*8) = 2*rank / (N_gen*N_eig - N_gen*rank)
#       = 2*rank / (N_gen * (N_eig - rank))
#       = 2*rank / (N_gen * 1)     [since N_eig - rank = 9-8 = 1]
#       = 16/3

print(f"  Axion formula: Omega = 2N/(27-3N)")
print(f"  Rewrite 27 = N_gen * N_eig = {N_gen} * {N_eig} = {N_gen * N_eig}")
print(f"  With N = rank(E8) = {rank_E8}:")
print(f"    Omega = 2*rank / (N_gen*N_eig - N_gen*rank)")
print(f"          = 2*{rank_E8} / ({N_gen}*{N_eig} - {N_gen}*{rank_E8})")
print(f"          = {2*rank_E8} / ({N_gen*N_eig} - {N_gen*rank_E8})")
print(f"          = {2*rank_E8} / {N_gen*N_eig - N_gen*rank_E8}")
print(f"          = {2*rank_E8/(N_gen*N_eig - N_gen*rank_E8):.6f}")
print(f"          = 2*rank(E8) / N_gen  [since N_eig - rank = 1 for a TREE]")
print(f"          = {2*rank_E8/N_gen:.6f}")

# KEY INSIGHT: N_eig - rank = 9 - 8 = 1 because the McKay graph is a TREE
# (a tree on n nodes has n-1 edges)
# So the denominator simplifies to N_gen * 1 = N_gen = 3

print(f"\n  KEY: N_eig - rank(E8) = {N_eig} - {rank_E8} = {N_eig - rank_E8}")
print(f"  This = 1 because the McKay graph (affine E8) is a TREE")
print(f"  A tree on V vertices has V-1 edges: N_eig - 1 = rank(E8)")
print(f"  So: Omega = 2*(N_eig-1)/N_gen = 2*rank/N_gen = 16/3")

# ============================================================
# The number 27 in context
# ============================================================
print("\n" + "="*70)
print("  THE NUMBER 27")
print("="*70)

# In the axion model: 27 = 3 * beta_3 where beta_3 = 9 is QCD beta coefficient
# In our framework: 27 = N_gen * N_eig = 3 * 9

# But in QCD: beta_3 = (11*N_c - 2*n_f)/3 = (11*3 - 2*3)/3 = 27/3 = 9
# With N_c = 3 = N_gen and n_f = 3 (light flavors) = N_gen
# So: beta_3 = (11*N_gen - 2*N_gen)/3 = 9*N_gen/3 = 3*N_gen = 9

# Hmm, that gives beta_3 = 3*N_gen = 9, and 27 = 3*9 = 9*N_gen = N_gen * 3*N_gen?
# No: 27 = 3*beta_3 = 3*9 = 27. But in the axion formula it's just "27" as a constant.

# In our framework: N_eig = 9 = the number of irreps of 2I
# And beta_3 = (11*3-2*3)/3 = 9 = N_eig (!)

print(f"  QCD beta_3 = (11*N_c - 2*n_f)/3 = (33-6)/3 = 9")
print(f"  N_eig (2I irreps) = 9")
print(f"  beta_3 = N_eig?  9 = 9  YES")
print(f"\n  This is suggestive: the QCD one-loop beta coefficient equals")
print(f"  the number of irreps of the binary icosahedral group.")
print(f"  Both involve the number 9 = N_eig = a1^2 - (a1^2-a1-1+1) = ...")
print(f"  Actually: 9 = N_eig = sum of Galois orbit sizes 1+2+2 + 1+3 = ...")
print(f"  N_eig = (a1^2+3)/a1+1... no, N_eig = 9 for 2I specifically.")

# N_eig for 2I: the binary icosahedral group has |2I|/sum(d_i^2) irreps
# sum(d_i^2) = |G| for any finite group (Burnside), so sum(d_i^2) = 120
# Number of irreps = number of conjugacy classes = 9

# For SU(3)_QCD: beta_3 = (11*3 - 2*3)/3 = 9
# With N_c = 3 (from the framework's N_gen)
# and n_f = 3 (the number of light quark flavors = N_gen)

# So beta_3 = (11*N_gen - 2*N_gen)/3 = 3*N_gen = 9 = N_eig
# This gives: 3*N_gen = N_eig, i.e., N_eig = 3*N_gen = 9

print(f"\n  Connection chain:")
print(f"    beta_3 = (11*N_c - 2*n_f)/3")
print(f"    With N_c = N_gen = 3 and n_f = N_gen = 3:")
print(f"    beta_3 = (11*3 - 2*3)/3 = 27/3 = 9 = N_eig")
print(f"    So: 27 = 3*beta_3 = 3*N_eig = N_gen*N_eig = {N_gen*N_eig}")

# ============================================================
# Framework formula for Omega
# ============================================================
print("\n" + "="*70)
print("  FRAMEWORK FORMULA")
print("="*70)

# Combining:
# Omega = 2*rank(E8) / (N_gen*N_eig - N_gen*rank(E8))
#       = 2*(N_eig-1) / (N_gen*(N_eig - N_eig + 1))  [rank = N_eig-1 for tree]
#       = 2*(N_eig-1) / N_gen
# OR equivalently:
# Omega = 2*|edges(McKay)| / N_gen

Omega_mckay = 2 * rank_E8 / N_gen
print(f"  Omega = 2*|edges(McKay)| / N_gen = 2*{rank_E8}/{N_gen} = {Omega_mckay:.6f}")
print(f"  = 2*(N_eig-1)/N_gen = 2*{N_eig-1}/{N_gen}")
print(f"  = 16/3 = {16/3:.6f}")

# But our PATTERN is 7-phi = 5.382, not 16/3 = 5.333
# Which is right?
print(f"\n  COMPARISON:")
print(f"    McKay formula: 16/3 = {16/3:.4f} ({abs(16/3-obs)/obs_err:.2f}s from obs)")
print(f"    Pattern 7-phi: {7-phi:.4f} ({abs(7-phi-obs)/obs_err:.2f}s from obs)")
print(f"    Observed: {obs} +/- {obs_err}")
print(f"\n  Both within 1 sigma. Future precision will discriminate.")

# ============================================================
# Could the CORRECT formula be 2*rank/N_gen?
# ============================================================
print("\n" + "="*70)
print("  WHICH PREDICTION IS RIGHT?")
print("="*70)

# 7-phi = b1+1-phi = 6+1-phi = 7-phi
# 16/3 = 2*rank/N_gen = 2*8/3

# Are these related?
# 7-phi = 5.382, 16/3 = 5.333
# Difference = 7-phi-16/3 = 7-phi-16/3 = (21-3phi-16)/3 = (5-3phi)/3
delta = 7 - phi - 16/3
print(f"  7-phi - 16/3 = {delta:.6f} = (5-3*phi)/3 = {(5-3*phi)/3:.6f}")
print(f"  = (a1-3*phi)/3 = (a1-3*phi)/N_gen")
print(f"  = (a1 - 3*(1+sqrt(5))/2) / N_gen")

# Not obviously meaningful. They're just two different numbers.

# The AXION formula has a derivation chain (beta functions, relaxation).
# Our 7-phi is a pattern (no mechanism).
# The McKay rewriting 2*rank/N_gen = 16/3 is a TRANSLATION of the axion formula.

print(f"\n  Honest assessment:")
print(f"  - The axion formula Omega = 2N/(27-3N) has a physical derivation.")
print(f"  - Setting N = rank(E8) and 27 = N_gen*N_eig is a TRANSLATION")
print(f"    into framework language, not an independent derivation.")
print(f"  - The simplification 2*rank/N_gen works because N_eig-rank=1 (tree).")
print(f"  - This is STRUCTURAL (algebraic identity) not DERIVED (no dynamics).")
print(f"  - Our pattern 7-phi = 5.382 is CLOSER to observation but has NO mechanism.")
print(f"  - The axion's 16/3 = 5.333 is slightly further but HAS a mechanism.")

# ============================================================
# FINAL VERDICT
# ============================================================
print("\n" + "="*70)
print("  VERDICT")
print("="*70)
print(f"""
  The composite axion model gives Omega_DM/Omega_B = 2N/(27-3N) with N=8.

  In framework language:
    N = rank(E8) = 8 (number of McKay edges)
    27 = N_gen * N_eig = 3*9 (SM beta coefficient = number of 2I irreps)
    Omega = 2*rank(E8)/N_gen = 16/3 = 5.333

  The simplification to 2*rank/N_gen uses N_eig - rank = 1 (McKay tree).

  This provides a possible MECHANISM for the DM/baryon ratio:
  if the confining gauge group has N_colors = rank(E8) = 8,
  the axion relaxation mechanism gives the right ratio.

  HOWEVER: This is a translation, not a derivation from the 600-cell.
  The 600-cell framework does not predict the existence of a confining
  SU(8) sector or the composite axion mechanism.

  The framework's own prediction 7-phi = 5.382 remains a PATTERN.
  The axion's 16/3 = 5.333 has a physical mechanism but N=8 is ad hoc
  (in the original paper). In our framework, N=8 WOULD be derived
  (as rank(E8)), but the axion mechanism itself is not part of the framework.

  STATUS: SUGGESTIVE CONNECTION. Not a derivation. Worth citing.
  The coincidence 27 = N_gen * N_eig is the most interesting structural link.
""")

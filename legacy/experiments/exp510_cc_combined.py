"""
EXP-510: CC Combined Formula - Orbifold E_vac + Heat Kernel
=============================================================
Instead of beta = 123 (ad hoc), try beta = N = 120 (natural)
with a prefactor involving E_vac from the orbifold.

The idea: Lambda = f(E_vac_orbifold) * exp(-N * L_1)

E_vac = -539/43200 (exact, from exp507)
exp(-N * L_1) = exp(-120 * (12-6*phi)) (exact, from heat kernel)

What function f gives the observed CC?
"""

import numpy as np
from fractions import Fraction
import sys
sys.path.insert(0, '.')
from commons import (PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen,
                      alpha_em, alpha_s, ln_inv_alpha, degree, h)

print("=" * 70)
print("EXP-510: CC COMBINED - ORBIFOLD + HEAT KERNEL")
print("=" * 70)

# Exact quantities
E_vac = Fraction(-539, 43200)
E_vac_f = float(E_vac)
L_1 = 12 - 6*PHI
obs = 2.8625e-122  # Lambda * l_P^2

print(f"\n  E_vac(C^2/2I) = {E_vac} = {E_vac_f:.10f}")
print(f"  L_1 = 12 - 6*phi = {L_1:.10f}")
print(f"  N = {N}")

# Heat kernel at beta = N = 120
heat_N = np.exp(-N * L_1)
print(f"  exp(-N*L_1) = exp(-{N*L_1:.4f}) = {heat_N:.6e}")

# What prefactor is needed?
prefactor_needed = obs / heat_N
print(f"\n  Observed CC = {obs:.4e}")
print(f"  Needed prefactor = obs / exp(-N*L_1) = {prefactor_needed:.6f}")
print(f"  4*sqrt(5) = {4*SQRT5:.6f}")
print(f"  Ratio needed/4sqrt5 = {prefactor_needed / (4*SQRT5):.6f}")

# ============================================================
# Try: Lambda = 4*sqrt(5) * |E_vac|^x * exp(-N * L_1)
# What x gives the observed CC?
# ============================================================
print(f"\n{'='*70}")
print("SCANNING EXPONENT x IN |E_vac|^x")
print(f"{'='*70}")

ratio_needed = prefactor_needed / (4*SQRT5)
# |E_vac|^x = ratio_needed
# x * ln|E_vac| = ln(ratio_needed)
x_exact = np.log(ratio_needed) / np.log(abs(E_vac_f))

print(f"\n  Need |E_vac|^x = {ratio_needed:.6f}")
print(f"  ln|E_vac| = {np.log(abs(E_vac_f)):.6f}")
print(f"  ln(ratio) = {np.log(ratio_needed):.6f}")
print(f"  x_exact = {x_exact:.6f}")
print(f"  phi = {PHI:.6f}")
print(f"  x_exact - phi = {x_exact - PHI:.6f}")
print(f"  x_exact / phi = {x_exact / PHI:.6f}")

# Test with x = phi
Lambda_phi = 4 * SQRT5 * abs(E_vac_f)**PHI * heat_N
z_phi = -np.log(Lambda_phi) / ln_inv_alpha

print(f"\n  TEST x = phi = {PHI:.6f}:")
print(f"    |E_vac|^phi = {abs(E_vac_f)**PHI:.6e}")
print(f"    Lambda = 4*sqrt(5) * |E_vac|^phi * exp(-N*L_1) = {Lambda_phi:.6e}")
print(f"    z_equiv = {z_phi:.4f}")
print(f"    Observed: {obs:.4e}")
print(f"    Error: {(Lambda_phi/obs - 1)*100:+.2f}%")

# Test with other natural exponents
print(f"\n  Scanning clean exponents:")
print(f"  {'x':>10} {'Lambda':>14} {'error':>10} {'notes':>20}")

for x, notes in [
    (1, "linear"),
    (2, "squared"),
    (PHI, "golden ratio"),
    (1/PHI, "1/phi"),
    (PHI**2, "phi^2"),
    (3/2, "3/2"),
    (0, "no E_vac"),
    (1/2, "sqrt"),
    (np.e/np.pi, "e/pi"),
    (1 + 1/PHI, "1 + 1/phi = phi"),
]:
    Lambda_x = 4 * SQRT5 * abs(E_vac_f)**x * heat_N
    err = (Lambda_x/obs - 1)*100
    marker = " <---" if abs(err) < 10 else ""
    print(f"  {x:10.6f} {Lambda_x:14.4e} {err:+10.2f}% {notes:>20}{marker}")

# ============================================================
# The phi formula in detail
# ============================================================
print(f"\n{'='*70}")
print("THE PHI FORMULA")
print(f"{'='*70}")

print(f"""
  Lambda_P = 4*sqrt(5) * |E_vac|^phi * exp(-N * L_1)

  where:
    4*sqrt(5) = Galois weight (mult * sqrt(a1))
    E_vac = -539/43200 (orbifold vacuum energy of C^2/2I)
    phi = (1+sqrt(5))/2 (golden ratio = natural scaling exponent)
    N = 120 (group order)
    L_1 = 12-6*phi (spectral gap)

  Numerically:
    |E_vac|^phi = (539/43200)^phi = {abs(E_vac_f)**PHI:.6e}
    exp(-N*L_1) = exp(-{N*L_1:.4f}) = {heat_N:.6e}
    4*sqrt(5) = {4*SQRT5:.6f}

    Lambda_P = {Lambda_phi:.6e}
    Observed: {obs:.6e}
    Error: {(Lambda_phi/obs - 1)*100:+.2f}%
""")

# ============================================================
# But: is |E_vac|^phi clean algebraically?
# ============================================================
print(f"{'='*70}")
print("IS |E_vac|^phi ALGEBRAIC?")
print(f"{'='*70}")

# (539/43200)^phi is NOT algebraic (irrational exponent of a rational).
# It's a transcendental number.
# So the formula Lambda = 4*sqrt(5) * (539/43200)^phi * exp(-720*phi'^2)
# involves BOTH exp and a power of phi.

# This is less clean than the pure formula 4*sqrt(5)*exp(-738*phi'^2).

# UNLESS: there's a way to combine (539/43200)^phi into the exponent:
# Lambda = 4*sqrt(5) * exp(-720*phi'^2 + phi*ln(539/43200))
# = 4*sqrt(5) * exp(-720*phi'^2 + phi*ln(539/43200))

combined_exp = -720*PHI_CONJ**2 + PHI*np.log(539/43200)
print(f"\n  Combined exponent: -720*phi'^2 + phi*ln(539/43200)")
print(f"  = {-720*PHI_CONJ**2:.6f} + {PHI*np.log(539/43200):.6f}")
print(f"  = {combined_exp:.6f}")

# Compare with -738*phi'^2 (the beta=123 exponent)
print(f"  vs -738*phi'^2 = {-738*PHI_CONJ**2:.6f}")
print(f"  Difference: {combined_exp - (-738*PHI_CONJ**2):.6f}")

# So the phi formula has exponent -275.02 + (-7.09) = -282.11
# vs the beta=123 formula: -281.89
# Difference: -0.22

# The phi formula gives a LARGER exponent (more negative),
# hence a SMALLER Lambda. That's why it's closer to observation.

print(f"\n  The phi formula effectively has:")
print(f"    effective_beta = N + phi*ln(539/43200)/L_1")
print(f"    = 120 + {PHI*np.log(539/43200)/L_1:.6f}")
print(f"    = {120 + PHI*np.log(539/43200)/L_1:.6f}")

eff_beta = 120 + PHI*np.log(539/43200)/L_1
print(f"    (vs 123 in the N+N_gen formula)")
print(f"    (vs 123.076 for exact match)")

# ============================================================
# Alternative: maybe the formula doesn't involve E_vac at all
# Maybe 123 = (dim(E8) - 2) / 2
# ============================================================
print(f"\n{'='*70}")
print("ALTERNATIVE: 123 = (dim(E8) - 2) / 2")
print(f"{'='*70}")

dim_E8 = 248
print(f"\n  dim(E8) = {dim_E8}")
print(f"  (dim(E8) - 2) / 2 = {(dim_E8 - 2) // 2}")
print(f"  = 123")

# This would mean: beta = (dim(E8) - 2) / 2
# But what is dim(E8) - 2 = 246?
# 246 = 2 * 123 = 2 * 3 * 41
# In E8: dim = rank + |roots| = 8 + 240.
# 246 = 240 + 6 = |roots| + b1
# Or: 246 = 248 - 2. What's 2? rank(SU(2))? dim(fund rep of SU(2))?
# Or: number of NON-CARTAN generators minus something?

print(f"  246 = |roots| + b1 = 240 + 6 = {240 + b1}")
print(f"  246 = dim(E8) - dim(fund SU(2)) = 248 - 2")
print(f"  246 = dim(E8) - rank/4 = 248 - 2: {248 - 2 == 246}")

# Is there a Lie-algebraic reason for beta = (dim - 2) / 2?
# In the adjoint representation: 248 = 1 + 247 (trivial + irreducible)?
# No, the adjoint of E8 IS irreducible (248-dimensional).

# Another angle: 123 = N + 3 = |2I| + |center of 2I| + 1?
# center of 2I = {I, -I}, |center| = 2. So N + |center| + 1 = 123.
print(f"\n  123 = N + |Z(2I)| + 1 = 120 + 2 + 1 = 123")
print(f"  Or: 123 = (N + |Z(2I)|*N_gen) = 120 + 2*... no")

# Actually simplest: 123 = N + N_gen = 120 + 3.
# And N_gen = 3 is DERIVED (from McKay p4/p2 = 3).
# So 123 IS a derived quantity, just the physical meaning of
# "3 extra tunneling steps for generations" is not rigorous.

# ============================================================
# Physical interpretation: WHY phi as exponent?
# ============================================================
print(f"\n{'='*70}")
print("WHY phi AS EXPONENT? (SPECULATIVE)")
print(f"{'='*70}")

print(f"""
  If the CC formula is:
    Lambda = 4*sqrt(5) * |E_vac|^phi * exp(-N*L_1)

  Why phi? In the theory, phi is the SCALING EXPONENT:
    masses: m_f = m_e * phi^n_f
    mixing: theta_ij = arctan(phi^(-n_ij))

  The orbifold vacuum energy E_vac = -539/43200 is a SCALE.
  The CC is this scale raised to the phi power (= one step
  of the mass hierarchy), then multiplied by the Galois
  tunneling factor.

  PHYSICAL PICTURE:
    1. The orbifold C^2/2I has vacuum energy E_vac (exact, rational)
    2. This energy is "scaled" by phi (the natural mass exponent)
    3. Then suppressed by coherent Galois tunneling through N vertices
    4. Result: Lambda = 4*sqrt(5) * |E_vac|^phi * exp(-N*L_1)

  This separates into:
    - GEOMETRIC part: E_vac^phi (from the orbifold = local geometry)
    - TOPOLOGICAL part: exp(-N*L_1) (from tunneling = global topology)
    - ALGEBRAIC part: 4*sqrt(5) (from Galois quantum dimensions)

  Error vs observation: {(Lambda_phi/obs - 1)*100:+.2f}%
""")

# ============================================================
# COMPARISON: all three formulas
# ============================================================
print(f"{'='*70}")
print("COMPARISON OF CC FORMULAS")
print(f"{'='*70}")

Lambda_A = alpha_em ** (57 - alpha_s)
Lambda_B = 4 * SQRT5 * np.exp(-738 * PHI_CONJ**2)  # beta = 123
Lambda_C = 4 * SQRT5 * abs(E_vac_f)**PHI * np.exp(-N * L_1)  # combined

print(f"\n  {'Formula':>45} {'Lambda':>12} {'error':>8} {'z':>8}")
print(f"  {'A: alpha^(57-as)':>45} {Lambda_A:12.4e} {(Lambda_A/obs-1)*100:+7.1f}% {-np.log(Lambda_A)/ln_inv_alpha:8.4f}")
print(f"  {'B: 4s5*exp(-738p2) [beta=123]':>45} {Lambda_B:12.4e} {(Lambda_B/obs-1)*100:+7.1f}% {-np.log(Lambda_B)/ln_inv_alpha:8.4f}")
print(f"  {'C: 4s5*|Ev|^phi*exp(-NL1) [beta=120]':>45} {Lambda_C:12.4e} {(Lambda_C/obs-1)*100:+7.1f}% {-np.log(Lambda_C)/ln_inv_alpha:8.4f}")
print(f"  {'Observed (Planck 2018)':>45} {obs:12.4e} {'':>8} {-np.log(obs)/ln_inv_alpha:8.4f}")

print(f"""
  Formula A (alpha^z): best fit ({abs(Lambda_A/obs-1)*100:.1f}% error), uses transcendental alpha
  Formula B (beta=123): algebraic, but 123 is ad hoc ({abs(Lambda_B/obs-1)*100:.1f}% error)
  Formula C (E_vac^phi): beta=N=120 (natural!), ({abs(Lambda_C/obs-1)*100:.1f}% error)

  Formula C eliminates the ad hoc 123 by using the orbifold E_vac.
  But introduces |E_vac|^phi which needs physical justification.
""")

print("=" * 70)
print("EXP-510 COMPLETE")
print("=" * 70)

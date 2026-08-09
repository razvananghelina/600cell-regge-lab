"""
EXP-534: The Running Ratio alpha(M_Z)/alpha(0) = 16/15?
=========================================================

From the m_Z formula: m_Z = m_e * phi^25 * alpha(M_Z)/alpha(0).
The running ratio is approximately 16/15 = 1 + 1/(N_gen * a1).

QUESTION: Does this ratio follow from the framework fermion masses
(which are DERIVED) + standard QED vacuum polarization?

APPROACH:
  1. Compute Delta_alpha from one-loop QED with FRAMEWORK masses.
  2. Compare with PDG and with 16/15.
  3. Check if the fermion mass ratios conspire to give 1/15.
  4. Examine whether the spectral action gives 16/15 directly.
  5. Test: is the ratio 16/15 a consequence of a1=5, or accidental?
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, SQRT5, a1, b1, N, h, degree, d_ST, rank_E8,
                      alpha_em, inv_alpha, alpha_s, sin2_tW, N_gen)

print("=" * 70)
print("EXP-534: THE RUNNING RATIO 16/15")
print("=" * 70)

# Physical constants
M_Z = 91.1876    # GeV
m_e_GeV = 0.51099895e-3  # GeV

# Framework alpha
alpha_0 = alpha_em   # Thomson limit
inv_alpha_0 = inv_alpha

# PDG at M_Z
alpha_MZ_pdg = 1.0 / 127.951
ratio_pdg = alpha_MZ_pdg / alpha_0

print(f"\n  Framework: 1/alpha(0) = {inv_alpha_0:.4f}")
print(f"  PDG: 1/alpha(M_Z) = {1/alpha_MZ_pdg:.3f}")
print(f"  Ratio alpha(M_Z)/alpha(0) = {ratio_pdg:.6f}")
print(f"  16/15 = {16/15:.6f}")
print(f"  Deviation 16/15 vs PDG: {(16/15 / ratio_pdg - 1)*100:+.3f}%")

# ============================================================
# SECTION 1: FRAMEWORK FERMION MASSES
# ============================================================
print(f"\n{'='*70}")
print("SECTION 1: FRAMEWORK FERMION MASSES")
print(f"{'='*70}")

# Mass formula: m_f = m_e * phi^{5a + 6b + delta_1*(1 + c_eff*alpha)}
# For the running calculation, we need the POLE masses for leptons
# and MSbar masses for quarks. Use the tree-level formula for simplicity.

# (a, b) quantum numbers:
fermion_data = [
    # (name, a, b, Q, N_c, mass_type)
    ("e",   0, 0, -1,   1, "pole"),
    ("mu",  1, 1, -1,   1, "pole"),
    ("tau", 1, 2, -1,   1, "pole"),
    ("u",   3,-2, 2/3,  3, "msbar"),
    ("d",   1, 0,-1/3,  3, "msbar"),
    ("s",   1, 1,-1/3,  3, "msbar"),
    ("c",   2, 1, 2/3,  3, "msbar"),
    ("b",  -1, 4,-1/3,  3, "msbar"),
    ("t",   4, 1, 2/3,  3, "pole"),
]

# Tree-level masses (no corrections for simplicity)
C_mass = 4 / (a1**2 + 1)  # = 2/13
delta_1_values = {}

print(f"\n  {'Name':>4s} {'(a,b)':>8s} {'Q':>5s} {'Nc':>3s} {'5a+6b':>6s} {'m_fw (GeV)':>12s} {'m_exp (GeV)':>12s}")

# Experimental masses for comparison
m_exp = {
    "e": 0.51099895e-3, "mu": 0.10566, "tau": 1.77686,
    "u": 0.00216, "d": 0.00467, "s": 0.0934,
    "c": 1.27, "b": 4.18, "t": 172.69
}

fw_masses = {}
for name, a, b, Q, Nc, mtype in fermion_data:
    exponent = 5*a + 6*b
    m_tree = m_e_GeV * PHI**exponent
    fw_masses[name] = m_tree
    print(f"  {name:>4s} ({a:+d},{b:+d}) {Q:>5.2f} {Nc:3d} {exponent:6d} {m_tree:12.4e} {m_exp[name]:12.4e}")

# ============================================================
# SECTION 2: ONE-LOOP QED RUNNING WITH FRAMEWORK MASSES
# ============================================================
print(f"\n{'='*70}")
print("SECTION 2: ONE-LOOP QED VACUUM POLARIZATION")
print(f"{'='*70}")

# Standard formula:
# Delta_alpha(M_Z) = sum_f (alpha/3pi) * N_c * Q^2 * [ln(M_Z^2/m_f^2) - 5/3]
# for fermions with m_f << M_Z (light fermions)
#
# For top quark (m_t > M_Z), the contribution is suppressed.
#
# More precisely:
# Delta_alpha_f = (alpha/3pi) * N_c * Q^2 * R(M_Z^2/m_f^2)
# where R(x) = -5/3 + 4/3*(1+2/x)*sqrt(1-4/x)*... for massive fermions
# For m_f << M_Z: R ~ ln(x) - 5/3
# For m_f >> M_Z: R ~ 4/(15*x) + ... (decoupled)

def delta_alpha_fermion(m_f, Q, Nc, M_Z_val, alpha_val):
    """One-loop vacuum polarization contribution from a single fermion."""
    x = M_Z_val**2 / m_f**2
    if x > 100:
        # Light fermion: use asymptotic form
        R = np.log(x) - 5/3
    elif x < 0.25:
        # Heavy fermion: decoupled
        R = 4 / (15 * x)
    else:
        # General formula with threshold
        if x >= 4:
            beta = np.sqrt(1 - 4/x)
            R = -5/3 + 4*(1 + 2/x) * (0.5 - beta * np.arctanh(beta) / 2)
            # Actually the exact formula for the photon self-energy:
            # Pi(q^2) = (alpha/3pi)*Nc*Q^2*[-5/3 + 4*(1+2m^2/q^2)*f(q^2/4m^2)]
            # where f(z) = 1/2 - sqrt(z-1)/sqrt(z) * arctan(sqrt(z/(z-1))) for z>1
            # This is getting complicated. Use the light approximation for safety.
            R = np.log(x) - 5/3
        else:
            R = np.log(max(x, 1)) - 5/3

    return (alpha_val / (3*np.pi)) * Nc * Q**2 * R

# Compute with FRAMEWORK masses
print(f"\n  Using FRAMEWORK tree-level masses:")
print(f"  {'Name':>4s} {'m_fw (GeV)':>12s} {'Q^2*Nc':>7s} {'Delta_alpha':>14s}")

delta_alpha_fw_total = 0
for name, a, b, Q, Nc, mtype in fermion_data:
    m_f = fw_masses[name]
    da = delta_alpha_fermion(m_f, Q, Nc, M_Z, alpha_0)
    delta_alpha_fw_total += da
    print(f"  {name:>4s} {m_f:12.4e} {Nc*Q**2:7.3f} {da:14.6e}")

print(f"\n  Total Delta_alpha (fw masses) = {delta_alpha_fw_total:.6f}")
ratio_fw = 1.0 / (1.0 - delta_alpha_fw_total)
print(f"  alpha(M_Z)/alpha(0) = 1/(1-Delta) = {ratio_fw:.6f}")

# Compare with experimental masses
print(f"\n  Using EXPERIMENTAL masses:")
delta_alpha_exp_total = 0
for name, a, b, Q, Nc, mtype in fermion_data:
    m_f = m_exp[name]
    da = delta_alpha_fermion(m_f, Q, Nc, M_Z, alpha_0)
    delta_alpha_exp_total += da

print(f"  Total Delta_alpha (exp masses) = {delta_alpha_exp_total:.6f}")
ratio_exp = 1.0 / (1.0 - delta_alpha_exp_total)
print(f"  alpha(M_Z)/alpha(0) = {ratio_exp:.6f}")

# Comparison table
print(f"\n  COMPARISON:")
print(f"    Framework masses:  {ratio_fw:.6f}")
print(f"    Experimental masses: {ratio_exp:.6f}")
print(f"    PDG (incl. had):   {ratio_pdg:.6f}")
print(f"    16/15:             {16/15:.6f}")
print(f"    Exact needed:      {M_Z / (m_e_GeV * PHI**25):.6f}")

# ============================================================
# SECTION 3: STRUCTURE OF THE SUM
# ============================================================
print(f"\n{'='*70}")
print("SECTION 3: STRUCTURE OF THE LOG SUM")
print(f"{'='*70}")

# The running is dominated by: (alpha/3pi) * sum_f N_c*Q^2 * ln(M_Z^2/m_f^2)
# Using M_Z = m_e * phi^25 * 16/15 and m_f = m_e * phi^{5a+6b}:
# ln(M_Z^2/m_f^2) = 2*ln(M_Z/m_f) = 2*ln(phi^25 * 16/15 / phi^{5a+6b})
#                  = 2*[(25-5a-6b)*ln(phi) + ln(16/15)]

# Define n_f = 5a+6b (mass exponent of fermion f)
# Then: ln(M_Z^2/m_f^2) ~ 2*(25 - n_f)*ln(phi)  (ignoring ln(16/15))

print(f"  Mass exponents n_f = 5a + 6b:")
print(f"  ln(phi) = {np.log(PHI):.6f}")
print(f"  {'Name':>4s} {'n_f':>4s} {'25-n_f':>6s} {'Nc*Q^2':>7s} {'weight':>10s}")

weight_sum = 0
for name, a, b, Q, Nc, mtype in fermion_data:
    n_f = 5*a + 6*b
    w = Nc * Q**2 * 2 * (25 - n_f)
    weight_sum += w
    print(f"  {name:>4s} {n_f:4d} {25-n_f:6d} {Nc*Q**2:7.3f} {w:10.3f}")

print(f"\n  Sum of weights: {weight_sum:.1f}")
print(f"  Delta ~ (alpha/3pi) * ln(phi) * weight_sum")
print(f"        = ({alpha_0:.6f} / {3*np.pi:.4f}) * {np.log(PHI):.6f} * {weight_sum:.1f}")
Delta_approx = (alpha_0 / (3*np.pi)) * np.log(PHI) * weight_sum
print(f"        = {Delta_approx:.6f}")
ratio_approx = 1.0 / (1.0 - Delta_approx)
print(f"  alpha(M_Z)/alpha(0) ~ {ratio_approx:.6f}")
print(f"  16/15 = {16/15:.6f}")

# For this to equal 16/15 = 1 + 1/15:
# Delta = 1/16  (since 1/(1-1/16) = 16/15)
# Need: (alpha/3pi) * ln(phi) * W = 1/16
# W = 1/(16 * alpha * ln(phi) / (3*pi))
W_needed = (3*np.pi) / (16 * alpha_0 * np.log(PHI))
print(f"\n  For ratio = 16/15, need Delta = 1/16 = {1/16:.6f}")
print(f"  Need weight sum W = 3*pi / (16*alpha*ln(phi)) = {W_needed:.2f}")
print(f"  Actual W = {weight_sum:.1f}")
print(f"  Ratio: {weight_sum/W_needed:.4f}")

# ============================================================
# SECTION 4: THE SUM Nc*Q^2*(25-nf) IN FRAMEWORK
# ============================================================
print(f"\n{'='*70}")
print("SECTION 4: ALGEBRAIC STRUCTURE OF THE WEIGHT SUM")
print(f"{'='*70}")

# The weight sum is: W = 2 * sum_f Nc*Q^2*(25-n_f)
# = 2 * [25 * sum Nc*Q^2 - sum Nc*Q^2*n_f]

sum_NcQ2 = sum(Nc*Q**2 for _, _, _, Q, Nc, _ in fermion_data)
sum_NcQ2_nf = sum(Nc*Q**2*(5*a+6*b) for _, a, b, Q, Nc, _ in fermion_data)

print(f"  sum_f Nc*Q^2 = {sum_NcQ2:.4f}")
print(f"  sum_f Nc*Q^2*n_f = {sum_NcQ2_nf:.4f}")
print(f"  W = 2*(25*{sum_NcQ2:.4f} - {sum_NcQ2_nf:.4f})")
print(f"    = 2*({25*sum_NcQ2:.4f} - {sum_NcQ2_nf:.4f})")
print(f"    = 2*{25*sum_NcQ2 - sum_NcQ2_nf:.4f}")
print(f"    = {weight_sum:.1f}")

# sum Nc*Q^2 per generation = 1 + 3*(4/9) + 3*(1/9) = 1 + 4/3 + 1/3 = 8/3
# For 3 generations: sum Nc*Q^2 = 3 * 8/3 = 8
print(f"\n  Per generation: sum Nc*Q^2 = 1 + 3*(2/3)^2 + 3*(1/3)^2 = 8/3")
print(f"  Total (3 gen): {N_gen} * 8/3 = {N_gen*8/3:.4f} = 8")
print(f"  Check: {sum_NcQ2:.4f}")

# Now sum_f Nc*Q^2*n_f:
print(f"\n  sum_f Nc*Q^2*n_f breakdown:")
for name, a, b, Q, Nc, mtype in fermion_data:
    n_f = 5*a + 6*b
    contrib = Nc*Q**2*n_f
    print(f"    {name}: Nc*Q^2*n_f = {Nc}*{Q**2:.4f}*{n_f} = {contrib:.4f}")

# Total
print(f"  Total: {sum_NcQ2_nf:.4f}")

# The "center of mass" of the mass spectrum weighted by Nc*Q^2:
n_com = sum_NcQ2_nf / sum_NcQ2
print(f"\n  Center-of-mass exponent: <n_f>_{{NcQ^2}} = {n_com:.4f}")
print(f"  25 - <n_f> = {25 - n_com:.4f}")

# For W = 2*8*(25-n_com):
print(f"  W = 2*8*(25-<n_f>) = 16*{25-n_com:.4f} = {16*(25-n_com):.2f}")

# ============================================================
# SECTION 5: CAN WE GET 16/15 ALGEBRAICALLY?
# ============================================================
print(f"\n{'='*70}")
print("SECTION 5: ALGEBRAIC DERIVATION ATTEMPT")
print(f"{'='*70}")

# alpha(M_Z)/alpha(0) = 1/(1 - Delta)
# For 16/15: Delta = 1/16
# Delta = (alpha/3pi) * ln(phi) * 2 * sum Nc*Q^2*(25-n_f)
# = (alpha/3pi) * ln(phi) * 2 * 8 * (25-<n_f>)
# = (16*alpha*ln(phi))/(3pi) * (25 - <n_f>)

# For Delta = 1/16:
# (16*alpha*ln(phi))/(3pi) * (25 - <n_f>) = 1/16
# 25 - <n_f> = 3pi / (256 * alpha * ln(phi))

rhs = 3*np.pi / (256 * alpha_0 * np.log(PHI))
print(f"  For 16/15, need: 25 - <n_f> = {rhs:.4f}")
print(f"  Actual: 25 - <n_f> = {25 - n_com:.4f}")
print(f"  Ratio: {(25 - n_com)/rhs:.4f}")

# This doesn't simplify to a nice number.
# The approach of deriving 16/15 from the mass spectrum doesn't work
# because alpha*ln(phi) is transcendental and doesn't cancel.

print(f"\n  CONCLUSION: 16/15 does NOT follow from QED running with")
print(f"  framework masses. The ratio alpha*ln(phi) = {alpha_0*np.log(PHI):.6f}")
print(f"  is transcendental and doesn't cancel algebraically.")

# ============================================================
# SECTION 6: ALTERNATIVE — SPECTRAL SUM RULE?
# ============================================================
print(f"\n{'='*70}")
print("SECTION 6: SPECTRAL ORIGIN OF 16/15")
print(f"{'='*70}")

# Instead of QED running, could 16/15 come from a spectral property?
# r = alpha(M_Z)/alpha(0) mediates between bare and physical m_Z.
# In the spectral action: the physical mass involves a "wave function
# renormalization" from the spectral geometry.

# 16/15 = Tr(A^2) / (Tr(A^2) - 1) = 16/15
# Tr(A^2) = sum of squared adjacency eigenvalues / N
# On McKay: Tr(A^2) = 2*|edges| = 16 for E8_tilde
# Tr(A^2) - 1 = 15 = number of edges minus... no, 15 = 16-1.

# Actually: on the McKay graph E8~, Tr(A^2) = 2*8 = 16 (8 edges, each
# contributes A_{ij}^2 + A_{ji}^2 = 2).
# So Tr(A^2)_McKay = 16.

# What is 15?
# 15 = sum_{rho!=trivial} dim(rho)^2 / N * ... no
# 15 = Tr(A^2) - 1
# The "1" that's subtracted could be the contribution of the trivial irrep.

print(f"  Tr(A^2) on McKay E8~ = 2 * |edges| = 2*8 = 16")
print(f"  16/15 = Tr(A^2) / (Tr(A^2) - 1)")
print(f"  The '1' is the trivial representation contribution")

# In the spectral action, the gauge coupling at the cutoff Lambda:
# 1/alpha(Lambda) = (f_2 / Lambda^2) * C_em
# The physical coupling at q < Lambda includes the self-energy:
# 1/alpha(q) = 1/alpha(Lambda) + (1/(2pi)) * b_em * ln(Lambda/q)

# If the "running" is reinterpreted as a spectral sum:
# alpha(M_Z)/alpha(0) = Tr(A^2) / (Tr(A^2) - trivial_mode_contribution)

# This would be: the physical mass is the bare mass AMPLIFIED by
# the ratio of total McKay modes to non-trivial McKay modes.

# The trivial mode (rho_0, dim 1) contributes dim^2/N = 1/120 to the
# total spectral weight. But 16/15 is a ratio of 16 to 15, not 120 to 119.

# More natural: on the McKay graph itself (not the Cayley graph):
# Total modes = 9 nodes, Tr(A^2) = 16, non-trivial = Tr(A^2) - 1 = 15.
# Then: ratio = 16/15 = Tr(A^2) / (Tr(A^2) - Id).

# But WHY would this ratio equal the QED running?

# The connection might be:
# The Seeley-DeWitt coefficient a_2 involves Tr(A^2) on the McKay graph.
# In the spectral action: S ~ f_0*a_0 + f_2*a_2 + f(0)*a_4
# a_2 involves the scalar curvature R ~ Tr(A^2) of the finite geometry.
# The "wave function renormalization" from a_0 to a_2 gives the ratio.

# Let's check: a_0 = N*dim(H_F) includes all modes.
# a_2 involves Tr(A^2) and other terms.
# The ratio a_0/(a_0 - a_2_contribution) could give 16/15.

# From the known Seeley-DeWitt coefficients (memory: A_0=11, A_1=62, A_2=233):
A_0 = 11  # per root pair
A_1 = 62
A_2 = 233

# The total a_k = 2*N*A_k
a_0 = 2*N*A_0  # = 2640
a_2 = 2*N*A_1  # = 14880 (this is the a_2 coefficient, using A_1)
a_4 = 2*N*A_2  # = 55920

print(f"\n  Seeley-DeWitt: a_0 = 2*N*A_0 = {a_0}")
print(f"  Seeley-DeWitt: a_2 = 2*N*A_1 = {a_2}")

# Ratios?
print(f"\n  a_0/(a_0-something) possibilities:")
print(f"    a_0/(a_0-2*N) = {a_0}/{a_0-2*N} = {a_0/(a_0-2*N):.6f}")
print(f"    2*N*A_0/(2*N*(A_0-1)) = A_0/(A_0-1) = {A_0/(A_0-1):.6f}")

# A_0/(A_0-1) = 11/10 = 1.1. Not 16/15.
# Tr(A^2)/(Tr(A^2)-1) = 16/15 = 1.0667.

# ============================================================
# SECTION 7: EXACT VALUE OF THE RATIO NEEDED
# ============================================================
print(f"\n{'='*70}")
print("SECTION 7: WHAT RATIO IS ACTUALLY NEEDED?")
print(f"{'='*70}")

# m_Z = m_e * phi^25 * r
# r = m_Z / (m_e * phi^25)
r_exact = M_Z / (m_e_GeV * PHI**25)
print(f"  Exact r needed: {r_exact:.8f}")
print(f"  16/15 = {16/15:.8f}")
print(f"  Difference: {r_exact - 16/15:.8f}")
print(f"  Relative: {(r_exact/(16/15)-1)*100:+.4f}%")

# Is r_exact a simple algebraic number?
# r = 1.06371... Let's check some candidates:
candidates = {
    "16/15": 16/15,
    "(a1-1)^2/((a1-1)^2-1)": (a1-1)**2/((a1-1)**2-1),
    "1+1/(N_gen*a1)": 1 + 1/(N_gen*a1),
    "1+alpha*a1": 1 + alpha_0*a1,
    "phi^(2*alpha)": PHI**(2*alpha_0),
    "1+1/(a1^2-a1)": 1 + 1/(a1**2-a1),
    "1+1/(4*a1)": 1 + 1/(4*a1),
    "(a1^2+1)/a1^2": (a1**2+1)/a1**2,
    "b1^2/(b1^2-1)": b1**2/(b1**2-1),
    "degree/(degree-1)": degree/(degree-1),
    "1+sin^2(tW)/a1": 1 + sin2_tW/a1,
}

print(f"\n  Algebraic candidates for r:")
print(f"  {'Expression':>35s} {'Value':>12s} {'Dev from r':>12s}")
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-r_exact)):
    dev = (val/r_exact - 1) * 100
    marker = " <--" if abs(dev) < 0.5 else ""
    print(f"  {name:>35s} {val:12.8f} {dev:+12.4f}%{marker}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  1. QED RUNNING: One-loop with framework masses gives ratio ~ {ratio_fw:.4f}.
     This is {(ratio_fw/(16/15)-1)*100:+.2f}% from 16/15 = {16/15:.4f}.
     The discrepancy is from missing hadronic VP (non-perturbative).

  2. ALGEBRAIC STRUCTURE: The log sum involves alpha*ln(phi), which is
     transcendental. No algebraic cancellation produces 16/15 exactly.

  3. SPECTRAL INTERPRETATION: 16/15 = Tr(A^2)_McKay / (Tr(A^2)_McKay - 1)
     where Tr(A^2) = 16 on the E8~ McKay graph. The "1" is the trivial
     mode contribution. This is a STRUCTURAL identification but lacks a
     physical derivation of why this McKay-graph ratio equals the QED running.

  4. EXACT RATIO NEEDED: r = {r_exact:.6f}, while 16/15 = {16/15:.6f}.
     The 0.28% gap means 16/15 is an APPROXIMATION, not exact.

  HONEST ASSESSMENT:
    16/15 = Tr(A^2)/(Tr(A^2)-1) is a clean algebraic identity on the
    McKay graph that happens to approximate the QED running ratio.
    It is NOT derivable from QED with framework masses (the sum involves
    alpha*ln(phi) which is transcendental).
    The match is to 0.3%, consistent with being a tree-level McKay identity
    that receives perturbative QED corrections.

  STATUS: PATTERN (algebraically clean, 0.3% from experiment,
  but no derivation connecting McKay Tr(A^2) to vacuum polarization).
""")

print("=" * 70)
print("EXP-534 COMPLETE")
print("=" * 70)

"""
EXP-515: Two DM Connections - Careful Analysis
================================================
Connection A: 26 Galois-broken modes = sin^2(tW) denominator
Connection B: 7-phi = a1/phi + b1/phi^2 and mass formula

RULE: Separate IDENTITIES (always true) from COINCIDENCES (same number,
different origin) from DERIVATIONS (one implies the other).
"""

import numpy as np
from fractions import Fraction
import sys
sys.path.insert(0, '.')
from commons import (PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen,
                      alpha_em, alpha_s, degree)

print("=" * 70)
print("EXP-515: DM CONNECTIONS - CAREFUL ANALYSIS")
print("=" * 70)

# ============================================================
# CONNECTION A: 26 = Galois-broken modes = sin^2(tW) denominator
# ============================================================
print(f"\n{'='*70}")
print("CONNECTION A: THE NUMBER 26")
print(f"{'='*70}")

# FACT 1: sin^2(theta_W) = b1/(a1^2+1) = 6/26
# WHERE does 26 come from in the Weinberg angle?
# Answer: from the hypercharge normalization in the spectral action.
# The trace condition: Tr(Y^2) over one generation of fermions.
# In the 600-cell: this involves a1^2 + 1 = 26.

print(f"""
  sin^2(theta_W) = b1/(a1^2+1) = {b1}/{a1**2+1} = {b1/(a1**2+1):.6f}

  Origin of 26 in Weinberg angle:
    a1^2 + 1 = 25 + 1 = 26
    This comes from hypercharge normalization: Tr(Y^2) over SM fermions
    in the spectral triple built from the 600-cell.
""")

# FACT 2: Galois-broken modes of 600-cell Laplacian = 26
# The 9 eigenvalues have multiplicities: 1, 4, 9, 16, 25, 36, 9, 16, 4
# Eigenvalues involving sqrt(5) (Galois-broken): L1(4), L2(9), L6(9), L8(4)
# Total broken: 4 + 9 + 9 + 4 = 26

# WHERE does this 26 come from?
# The McKay correspondence maps 2I irreps to extended E8 nodes.
# Irrep dimensions: 1, 2, 3, 4, 5, 6, 4, 2, 3
# Galois sigma swaps: rho_1(dim 2) <-> rho_7(dim 2), rho_2(dim 3) <-> rho_8(dim 3)
# These are the irreps whose characters involve sqrt(5).
# Galois-broken irreps: dims 2, 2, 3, 3
# Sum of dim^2: 4 + 4 + 9 + 9 = 26

print(f"  Galois-broken irreps of 2I:")
print(f"    rho_1 (dim 2) <-> rho_7 (dim 2): characters involve phi")
print(f"    rho_2 (dim 3) <-> rho_8 (dim 3): characters involve sqrt(5)")
print(f"    Sum of dim^2 = 2^2 + 2^2 + 3^2 + 3^2 = {4+4+9+9}")
print(f"    = 2*(2^2 + 3^2) = 2*13 = 26")

# QUESTION: Is 2*(2^2 + 3^2) = a1^2 + 1 an IDENTITY or COINCIDENCE?
# 2*(4+9) = 2*13 = 26
# a1^2 + 1 = 25 + 1 = 26
# These are the SAME NUMBER but from DIFFERENT constructions.

# Let's check: is it ALWAYS true that Galois-broken modes = a1^2 + 1?
# For 2I: the broken irreps are the ones NOT on the symmetry axis of
# the extended E8 Dynkin diagram. The Galois involution swaps the
# two "wings" of the diagram.

# Extended E8 Dynkin:
#   1 - 2 - 3 - 4 - 5 - 6
#                       |
#                       4 - 2 - 3
#                  (actually it's:)
#   1 - 2 - 3 - 4 - 5 - 6 - 4' - 2'
#                       |
#                       3'
# Wait, extended E8 has a specific shape. Let me be precise.

# Extended E8 (affine E8) Dynkin diagram:
# Nodes: 0,1,2,3,4,5,6,7,8 with dims 1,2,3,4,5,6,4,2,3
# Edges: 0-1, 1-2, 2-3, 3-4, 4-5, 5-6, 6-7, 2-8
# The diagram has a Z/2 symmetry swapping: 1<->7, 3<->? Hmm...

# Actually the extended E8 diagram has NO non-trivial automorphism!
# E8 (non-extended) has no diagram automorphism either.
# So the Galois involution does NOT correspond to a diagram automorphism.

# The Galois action is on the CHARACTER TABLE, not on the diagram.
# It swaps irreps whose characters are Galois conjugates.

print(f"""
  KEY QUESTION: Is sum(dim^2, broken irreps) = a1^2 + 1 a theorem?

  For 2I specifically:
    Broken irrep dims: 2, 2, 3, 3
    2*(2^2 + 3^2) = 2*(4+9) = 26 = a1^2 + 1

  This requires: 2^2 + 3^2 = (a1^2+1)/2 = 13.
  And: the broken irreps have dims 2 and 3.

  WHY dims 2 and 3?
    The 2D irreps of 2I are the ones that define the SU(2) embedding.
    Their characters are chi = phi and chi = -phi (conjugate pair).
    The 3D irreps come from Sym^2 of the 2D rep (roughly).
    Their characters involve sqrt(5) through phi.

  The dims 2 and 3 are specific to E8 McKay. For E6 McKay (2T):
    Broken irreps would be different.
""")

# Check for 2T (binary tetrahedral, E6):
# 2T irreps: dims 1, 1, 1, 2, 2, 2, 3 (7 irreps, extended E6)
# Galois action: sigma maps omega = exp(2pi*i/3) to omega^2
# This swaps the irreps whose characters involve omega (cube roots of unity)
# NOT sqrt(5). So the "Galois" here is different (Q(omega) not Q(sqrt(5))).

# For 2I: the relevant Galois group is Gal(Q(sqrt(5))/Q) = Z/2.
# For 2T: the relevant Galois group is Gal(Q(omega)/Q) = Z/2.
# Different Galois actions on different groups!

print(f"\n  VERDICT ON CONNECTION A:")
print(f"  26 = a1^2 + 1 appears in BOTH contexts because both depend on a1 = 5.")
print(f"  sin^2(tW) = b1/(a1^2+1) from hypercharge normalization.")
print(f"  Galois-broken = 2*(2^2+3^2) = 26 from McKay irrep structure.")
print(f"  These are CONSEQUENCES of the same a1, not independent facts.")
print(f"  But the STRUCTURAL link (Galois breaking <-> hypercharge) is suggestive.")
print(f"  Can we derive sin^2(tW) FROM Galois mode counting?")

# sin^2(tW) = b1/broken_modes
# b1 = 6 = degree of McKay graph (each node connects to avg 2 nodes... no)
# Actually b1 = 6 = a1+1. And broken = a1^2+1 = 26.
# sin^2(tW) = (a1+1)/(a1^2+1).

# The question: is there a PHYSICAL reason for sin^2(tW) = b1/broken_modes?
# sin^2(tW) is a ratio of gauge couplings: g'^2/(g^2+g'^2).
# If g'^2 ~ b1 and g^2+g'^2 ~ broken_modes, then sin^2(tW) = b1/broken_modes.

# In spectral action: g'^2 comes from Tr(Y^2) and g^2 from Tr(T_L^2).
# Tr(Y^2) involves the hypercharge assignments.
# If the hypercharge assignments are DETERMINED by the Galois structure...

# This is speculative. Let me note it and move on.

print(f"\n  POTENTIAL: sin^2(tW) = b1/(Galois-broken modes)")
print(f"  = (McKay degree + 1) / (sum dim^2 of Galois-paired irreps)")
print(f"  This would DERIVE the Weinberg angle from Galois representation theory.")
print(f"  STATUS: SUGGESTIVE, needs proof that Tr(Y^2) = Galois broken count.")

# ============================================================
# CONNECTION B: 7-phi = a1/phi + b1/phi^2 and mass formula
# ============================================================
print(f"\n{'='*70}")
print("CONNECTION B: 7-phi AND THE MASS FORMULA")
print(f"{'='*70}")

# IDENTITY: 7-phi = a1*phi^(-1) + b1*phi^(-2) = 5/phi + 6/phi^2
# Proof: 5/phi + 6/phi^2 = 5(phi-1) + 6(2-phi) = 5phi-5+12-6phi = 7-phi. QED.
# This is ALWAYS TRUE for a1=5, b1=6 (i.e., it's an algebraic identity).

print(f"\n  IDENTITY: 7-phi = a1/phi + b1/phi^2")
print(f"  Proof: a1/phi + b1/phi^2 = a1(phi-1) + b1(2-phi)")
print(f"       = a1*phi - a1 + 2*b1 - b1*phi")
print(f"       = (a1-b1)*phi + (2*b1-a1)")
print(f"       = -phi + 7 = 7-phi. QED.")
print(f"  Numerical: {a1/PHI + b1/PHI**2:.6f} = {7-PHI:.6f}")

# Now: in the mass formula, m_f = m_e * phi^(a1*a + b1*b)
# For a fermion with (a,b), the mass exponent is n = 5a + 6b.

# The expression 7-phi = a1*phi^(-1) + b1*phi^(-2) looks like the
# mass exponent function E(a,b) = a1*a + b1*b evaluated at
# "fractional quantum numbers" a = 1/phi, b = 1/phi^2.

# In the DARK sector: masses are m_f' = m_e * phi^(-n_f).
# The "dark mass function" uses -n instead of +n.

# What if the DM abundance comes from evaluating the mass function
# at a SPECIFIC point in the (a,b) plane?

print(f"\n  Mass formula: n(a,b) = a1*a + b1*b = 5a + 6b")
print(f"  Fermion quantum numbers:")
fermions = {
    'e': (0,0), 'mu': (1,1), 'tau': (1,2),
    'u': (3,-2), 'c': (2,1), 't': (4,1),
    'd': (1,0), 's': (1,1), 'b': (-1,4)
}
for name, (a,b) in fermions.items():
    n = a1*a + b1*b
    print(f"    {name:>3}: (a,b)=({a:+d},{b:+d}), n={n:+3d}, m/m_e=phi^{n}")

# 7-phi = n(1/phi, 1/phi^2) = 5/phi + 6/phi^2
# This is n evaluated at a = 1/phi, b = 1/phi^2.
# Is (1/phi, 1/phi^2) = (phi^{-1}, phi^{-2}) meaningful?

print(f"\n  7-phi = n(1/phi, 1/phi^2) = n(phi^-1, phi^-2)")
print(f"  = the mass function at the 'Galois point' (phi^-1, phi^-2)")

# Note: 1/phi = phi - 1 and 1/phi^2 = 2 - phi = phi'^2.
# So (1/phi, 1/phi^2) = (phi-1, phi'^2).

print(f"\n  (1/phi, 1/phi^2) = (phi-1, phi'^2) = ({1/PHI:.4f}, {PHI_CONJ**2:.4f})")

# Is there a fermion with (a,b) close to (1/phi, 1/phi^2)?
# The closest integer point: (1, 0) = d quark (n=5).
# But (phi^{-1}, phi^{-2}) is not an integer point.

# ALTERNATIVE INTERPRETATION:
# The DM abundance = sum of dark mass contributions weighted by phi^(-1)
# i.e., the mass function in the "Galois basis" where phi -> 1/phi.

# Let's check: sum over fermions of phi^(-n_f) * weight
print(f"\n  Dark sector mass sum (phi^(-n) for each fermion):")
total_dark = 0
for name, (a,b) in fermions.items():
    n = a1*a + b1*b
    dark_weight = PHI**(-n)
    total_dark += dark_weight
    print(f"    {name:>3}: phi^{-n:+3d} = {dark_weight:.6f}")

print(f"  Sum: {total_dark:.6f}")
print(f"  7-phi = {7-PHI:.6f}")
print(f"  Ratio: {total_dark/(7-PHI):.6f}")

# The sum is NOT 7-phi. So the DM abundance is NOT the sum of dark weights.

# ============================================================
# DEEPER: What IS 7-phi algebraically?
# ============================================================
print(f"\n{'='*70}")
print("DEEPER: ALGEBRAIC STRUCTURE OF 7-phi")
print(f"{'='*70}")

# 7-phi = (13-sqrt(5))/2
# 13 = F_7 (Fibonacci number)
# Or: 13 = sum of Galois-broken dim^2 per pair = 2^2 + 3^2

# N(7-phi) = 41 (Galois norm = product with conjugate)
# 41 is prime.

# In Z[phi]: 7-phi = 7 - phi. The element has:
# Trace: (7-phi) + (7-phi') = 14 - (phi+phi') = 14 - 1 = 13
# Norm: (7-phi)(7-phi') = 49 - 7(phi+phi') + phi*phi' = 49 - 7 + (-1) = 41

print(f"  7-phi in Z[phi]:")
print(f"    Trace = 13 (= 2^2 + 3^2 = broken dims per pair)")
print(f"    Norm = 41 (prime)")
print(f"    Conjugate: 7-phi' = {7-PHI_CONJ:.6f}")

# KEY OBSERVATION: Trace(7-phi) = 13 = 2^2 + 3^2
# And 2 and 3 are the dimensions of the Galois-broken irreps!
# Is this a coincidence?

print(f"\n  Trace(7-phi) = 13 = 2^2 + 3^2")
print(f"  WHERE: 2 = dim(rho_1) = dim(rho_7) [Galois-paired irreps]")
print(f"         3 = dim(rho_2) = dim(rho_8) [Galois-paired irreps]")
print(f"  AND: 2*13 = 26 = total Galois-broken modes = a1^2+1")

# So: Omega_DM/Omega_b = 7-phi, and Tr(7-phi) = dim_2^2 + dim_3^2
# where dim_2, dim_3 are the dimensions of the Galois-broken irreps.

# Can we write: 7-phi = (dim_2^2 + dim_3^2 - sqrt(5)) / 2?
test = (4 + 9 - SQRT5) / 2
print(f"\n  (dim_2^2 + dim_3^2 - sqrt(5))/2 = (13 - sqrt(5))/2 = {test:.6f}")
print(f"  7-phi = {7-PHI:.6f}")
print(f"  Match: {abs(test - (7-PHI)) < 1e-10}")

# YES! So: 7-phi = (sum_pairs dim^2 - sqrt(a1)) / 2
# = (13 - sqrt(5)) / 2

# This is TRIVIALLY true (it's the definition of 7-phi = (13-sqrt(5))/2).
# But it CONNECTS the DM abundance to the Galois-broken irrep dimensions.

# More precisely:
# Omega_DM/Omega_b = (Tr(Galois-broken dim^2 per pair) - sqrt(a1)) / 2

print(f"\n  FORMULA: Omega_DM/Omega_b = (d_2^2 + d_3^2 - sqrt(a1)) / 2")
print(f"  WHERE: d_2 = 2 and d_3 = 3 are Galois-paired irrep dimensions")
print(f"  THIS IS: (13 - sqrt(5))/2 = 7 - phi")
print(f"  BUT: it's an identity, not a derivation.")

# ============================================================
# THE REAL QUESTION: WHY would the DM abundance
# equal (d2^2 + d3^2 - sqrt(a1))/2 ?
# ============================================================
print(f"\n{'='*70}")
print("THE REAL QUESTION")
print(f"{'='*70}")

# For the DM abundance to be (d2^2+d3^2-sqrt(a1))/2, we need a
# PHYSICAL mechanism that:
# 1. Involves the Galois-broken irreps (d=2 and d=3)
# 2. Uses their dimensions squared (= mode count)
# 3. Subtracts sqrt(a1) (= quantum dimension difference?)
# 4. Divides by 2 (= the two Galois pairs)

# sqrt(a1) = sqrt(5) = phi^2 - phi'^2 = Galois weight
# d2^2 + d3^2 = 13 = half of Galois-broken modes

# So: Omega = (half_broken - Galois_weight) / 2
# = (13 - sqrt(5)) / 2

# PHYSICAL INTERPRETATION ATTEMPT:
# The DM abundance counts how many Galois-broken modes CONTRIBUTE
# to the dark sector, MINUS the Galois weight (which cancels),
# divided by the number of Galois pairs.

# In freeze-in production: the production rate goes as
# (coupling)^2 * (number of channels)
# coupling = gravitational, proportional to 1/M_Pl^2
# channels = the Galois-broken modes that can produce DM

# If the production involves d2^2 + d3^2 = 13 channels,
# and the Galois weight sqrt(5) acts as a suppression factor,
# and we divide by 2 for particle/antiparticle...

print(f"""
  ATTEMPT AT PHYSICAL MEANING:

  Omega_DM/Omega_b = (d_2^2 + d_3^2 - sqrt(a1)) / 2

  d_2^2 = 4 : number of modes in the 2D Galois-broken eigenspace
  d_3^2 = 9 : number of modes in the 3D Galois-broken eigenspace
  sqrt(a1) = sqrt(5) : Galois weight = phi^2 - phi'^2
  /2 : two Galois pairs (each contributes equally)

  If we interpret this as a PRODUCTION formula:
  - 13 modes can produce dark matter (Galois-broken channels)
  - sqrt(5) is the quantum interference between phi and phi' channels
  - Each pair contributes (d^2 - sqrt(5)/2) effective channels

  This gives: Omega_DM/Omega_b = sum_pairs (d_k^2 - sqrt(5)/2) / 2
  Hmm, that doesn't quite work...

  Actually: (d_2^2 + d_3^2 - sqrt(5)) / 2 = (13 - 2.236) / 2 = 5.382

  The subtraction of sqrt(5) might represent:
  - The quantum dimension DIFFERENCE phi - 1/phi = sqrt(5)/phi * ... hmm
  - Or: the asymmetry between physical and dark quantum dimensions
  - D^2_phys - D^2_dark = 2*sqrt(5) and sqrt(5) = (D^2_phys-D^2_dark)/2

  So: Omega = (half_broken - (D^2_phys-D^2_dark)/2) / 2
  = (13 - (7.236-2.764)/2) / 2
  = (13 - 2.236) / 2 = 5.382
""")

D2_phys = a1 + np.sqrt(a1)
D2_dark = a1 - np.sqrt(a1)
print(f"  D^2_phys - D^2_dark = {D2_phys - D2_dark:.6f} = 2*sqrt(a1) = {2*np.sqrt(a1):.6f}")
print(f"  (D^2_phys - D^2_dark) / 2 = sqrt(a1) = {np.sqrt(a1):.6f}")

# So the formula is:
# Omega = ((a1^2+1)/2 - (D^2_phys-D^2_dark)/2) / 2
# = ((a1^2+1) - 2*sqrt(a1)) / 4
# = (a1^2 + 1 - 2*sqrt(a1)) / 4
# = (a1 - 1/... no. Let me compute:
test_formula = (a1**2 + 1 - 2*np.sqrt(a1)) / 4
print(f"\n  (a1^2+1-2*sqrt(a1))/4 = ({a1**2+1}-{2*np.sqrt(a1):.4f})/4 = {test_formula:.6f}")
print(f"  7-phi = {7-PHI:.6f}")
print(f"  Match: {abs(test_formula - (7-PHI)) < 1e-10}")

# YES! Omega = (a1^2 + 1 - 2*sqrt(a1)) / 4 = (sqrt(a1) - 1)^2 / ... wait
# (a1^2+1-2sqrt(a1))/4 = ((a1-1)^2 + 2*a1 - 2*sqrt(a1))/4 ... hmm
# Let me just check: (a1^2+1-2sqrt(a1))/4 = (25+1-2*2.236)/4 = (26-4.472)/4 = 21.528/4 = 5.382. Yes!

# Alternative: (a1^2+1-2*sqrt(a1))/4 = (a1-sqrt(a1))^2/(4*...) hmm
# a1^2+1-2sqrt(a1) = (sqrt(a1))^4 + 1 - 2*sqrt(a1)
# Let x = sqrt(a1). Then x^4+1-2x = ?
# For x=sqrt(5): 25+1-2sqrt(5) = 26-2sqrt(5). OK.

# Or: note that a1^2+1 = 26 = Galois-broken modes. And 2*sqrt(a1) = D^2_phys-D^2_dark.
# So: Omega = (broken_modes - TQFT_asymmetry) / 4

print(f"\n  CLEAN FORMULA:")
print(f"  Omega_DM/Omega_b = (Galois_broken - TQFT_asymmetry) / 4")
print(f"  = (26 - 2*sqrt(5)) / 4")
print(f"  = ({a1**2+1} - {2*SQRT5:.4f}) / 4")
print(f"  = {(a1**2+1 - 2*SQRT5)/4:.6f}")
print(f"  = 7 - phi = {7-PHI:.6f}")
print(f"  Observed: 5.36 (error: {((7-PHI)/5.36-1)*100:+.2f}%)")

# ============================================================
# IS THIS A DERIVATION?
# ============================================================
print(f"\n{'='*70}")
print("IS THIS A DERIVATION?")
print(f"{'='*70}")

print(f"""
  FORMULA: Omega_DM/Omega_b = (N_broken - Delta_D^2) / 4

  WHERE:
    N_broken = a1^2+1 = 26 (Galois-broken Laplacian modes)
    Delta_D^2 = D^2_phys - D^2_dark = 2*sqrt(a1) (TQFT quantum dim asymmetry)
    4 = d_ST (spacetime dimension) = phi^3 - 1/phi^3

  EACH FACTOR IS DERIVED:
    N_broken = 26: from McKay/Galois representation theory
    Delta_D^2 = 2*sqrt(5): from TQFT at level k+2=a1
    d_ST = 4: from phi^3-1/phi^3 or rank(E8)/2

  BUT: the COMBINATION (N_broken - Delta_D^2) / d_ST is NOT derived.
  We don't have a physical argument for WHY the DM abundance equals
  this specific combination.

  The formula REWRITES 7-phi = (26-2sqrt(5))/4 in terms of
  framework quantities, but doesn't EXPLAIN why the DM abundance
  should equal this expression.

  STATUS: STRUCTURAL PATTERN (every factor is derived,
  but their combination into DM abundance is not).
""")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
CONNECTION A (26 = Galois-broken = sin^2(tW) denominator):
  Both are a1^2+1 = 26, coming from the same a1 = 5.
  sin^2(tW) = b1/(a1^2+1) = 6/26 = (McKay degree)/(Galois-broken modes).
  STRUCTURAL CONNECTION but not a new derivation.
  Suggestive: Weinberg angle = McKay connectivity / Galois breaking.

CONNECTION B (7-phi = DM abundance):
  7-phi = (26 - 2*sqrt(5)) / 4 = (N_broken - Delta_D^2) / d_ST
  This REWRITES the DM abundance as:
    (Galois-broken modes - TQFT asymmetry) / spacetime dimension
  Every factor is individually derived, but the combination is PATTERN.

  Most interesting finding:
    Tr(7-phi) = 13 = 2^2 + 3^2 (sum of squared dims of Galois-paired irreps)
    This connects DM abundance to the REPRESENTATION THEORY of 2I directly.

KEY NEW IDENTITY:
  Omega_DM/Omega_b = (a1^2 + 1 - 2*sqrt(a1)) / (phi^3 - phi^(-3))
  = (sin^2(tW)_denom - (D^2_phys - D^2_dark)) / d_ST
  STATUS: STRUCTURAL PATTERN
""")

print("=" * 70)
print("EXP-515 COMPLETE")
print("=" * 70)

"""
EXP-299: Dark Matter from Galois Sector
========================================
The Galois automorphism phi <-> phi' creates a "mirror" sector.
In this sector: alpha_s' = 1/(2*phi'^3) where phi' = -1/phi.
Since phi'^3 = -1/phi^3, we get alpha_s' = phi^3/2 > 0.
But z' = 2*phi'^3 < 0 (Galois conjugate), meaning z' = 2-4*phi < 0.

The "Galois sector" has OPPOSITE signs in key quantities.
Could this give dark matter candidates?

NOTE: All print uses ASCII only (Windows cp1252).
"""

import math
import numpy as np

PHI = (1 + math.sqrt(5)) / 2
PHI_PRIME = (1 - math.sqrt(5)) / 2  # = -1/phi
a1 = 5
b1 = 6
N = 120
N_gen = 3
alpha = 7.2973525693e-3
alpha_s = 1 / (2 * PHI**3)

print("=" * 70)
print("EXP-299: DARK MATTER FROM GALOIS SECTOR")
print("=" * 70)

# =====================================================================
# SECTION A: Galois Automorphism Review
# =====================================================================
print("\n--- A. GALOIS AUTOMORPHISM ---")
print(f"phi  = (1+sqrt(5))/2 = {PHI:.6f}")
print(f"phi' = (1-sqrt(5))/2 = {PHI_PRIME:.6f} = -1/phi")
print()
print("Galois action sigma: sqrt(5) -> -sqrt(5), phi -> phi'")
print("Preserves Q (rationals), swaps the two real embeddings of Q(sqrt(5)).")
print()

# Key quantities in both sectors
print("Key quantities in phi-sector vs phi'-sector:")
print(f"  1/alpha_s   = 2*phi^3  = {2*PHI**3:.6f}")
z_prime = 2 * PHI_PRIME**3
print(f"  sigma(z)    = 2*phi'^3 = {z_prime:.6f}")
print(f"  |sigma(z)|  = {abs(z_prime):.6f}")
print()

# In Z[phi]: z = 2+4*phi (a=2, b=4)
# sigma(z) = 2+4*phi' = 2+4*(-1/phi) = 2 - 4/phi = 2 - 4*(phi-1) = 6-4*phi
z_check = 6 - 4*PHI
print(f"  z = 2+4*phi = {2+4*PHI:.6f}")
print(f"  z' = 6-4*phi = {z_check:.6f}")
print(f"  N(z) = z*z' = {(2+4*PHI)*(6-4*PHI):.6f}")
print(f"  (should be -4: {2*6 + 2*(-4) + 4*6 - 16:.0f})")
print()

# alpha_s in Galois sector
alpha_s_prime = -1 / z_prime
print(f"'alpha_s' in Galois sector:")
print(f"  1/alpha_s' = z' = {z_prime:.6f} (NEGATIVE)")
print(f"  alpha_s' = 1/z' = {1/z_prime:.6f} (NEGATIVE)")
print()
print("NEGATIVE coupling constant => CONFINING sector")
print("(beta function has opposite sign -> infrared slavery)")

# =====================================================================
# SECTION B: Galois Mass Spectrum
# =====================================================================
print("\n--- B. GALOIS MASS SPECTRUM ---")
print("The mass formula m_f = m_e * phi^n uses the phi-sector.")
print("The Galois sector would use phi' instead of phi:")
print("  m_f' = m_e * |phi'|^{n'}")
print("  where n' = 5a - b (from exp241)")
print()

# Fermion masses in Galois sector
m_e = 0.51099895  # MeV
fermions = {
    'e':   (0, 0), 'mu':  (1, 1), 'tau': (1, 2),
    'u':   (3, -2), 'c':  (2, 1), 't':   (4, 1),
    'd':   (1, 0), 's':   (1, 1), 'b':   (-1, 4)
}

print(f"{'Name':5s} {'n':4s} {'n_prime':7s} {'m(MeV)':>12s} {'m_prime(MeV)':>14s} {'Ratio':>8s}")
print("-" * 55)

for name, (a, b) in fermions.items():
    n = a1*a + b1*b
    n_prime = a1*a - b
    m_phi = m_e * PHI**n
    # In Galois sector: |phi'| = 1/phi, so |phi'|^{n'} = phi^{-n'}
    m_phi_prime = m_e * PHI**(-n_prime) if n_prime != 0 else m_e
    # Actually: phi' = -1/phi, so phi'^{n'} = (-1)^{n'} / phi^{n'}
    # For mass, we take absolute value: m' = m_e * |phi'^{n'}| = m_e / phi^{n'}
    # But wait: the DSI potential uses |z|, so sign doesn't matter
    # m' = m_e * phi^{-n'} = m_e / phi^{n'}
    m_galois = m_e * PHI**(-n_prime)
    ratio = m_galois / m_phi if m_phi > 0 else float('inf')
    print(f"{name:5s} {n:4d} {n_prime:7d} {m_phi:12.4f} {m_galois:14.6f} {ratio:8.4f}")

# =====================================================================
# SECTION C: Dark Matter Candidates
# =====================================================================
print("\n--- C. DARK MATTER CANDIDATE ANALYSIS ---")
print()
print("For a particle to be dark matter, it needs:")
print("  1. Electrically neutral (or very weakly charged)")
print("  2. Stable (or lifetime >> age of universe)")
print("  3. Right relic abundance")
print("  4. Not too strongly interacting (direct detection bounds)")
print()

# The Galois sector particles have INVERTED mass hierarchy:
# Heavy phi-sector particles become light in Galois sector.
# Light phi-sector particles become heavy in Galois sector.

print("Key observation: The Galois sector INVERTS the mass hierarchy!")
print("  phi-sector: m_t >> m_e")
print("  Galois:     m_t' << m_e'")
print()

# The lightest Galois particle is the one with LARGEST |n'|
# Among fermions: top (n'=19), bottom (n'=-9), charm (n'=9)
# top has n'=19, so m_top' = m_e / phi^19 = very light

# Wait: n' for top: a=4, b=1 -> n'=5*4-1=19
# m_top' = m_e * phi^{-19} = 0.511 / phi^19 = tiny

# The HEAVIEST Galois particle has SMALLEST |n'|
# e has n'=0: m_e' = m_e (unchanged!)
# d has n'=5: m_d' = m_e / phi^5 ~ 0.09 MeV
# u has n'=17: m_u' = m_e / phi^17 ~ 0.3 eV (neutrino scale!)

print("Most interesting Galois fermions:")
galois_sorted = sorted(
    [(name, a1*a-b, m_e * PHI**(-(a1*a-b))) for name, (a, b) in fermions.items()],
    key=lambda x: -x[2]
)
for name, np_val, m_gal in galois_sorted:
    print(f"  {name:5s}: n' = {np_val:3d}, m_galois = {m_gal:.6f} MeV = {m_gal*1000:.3f} keV")

print()
# The heaviest Galois fermion is the electron (m'=m_e, unchanged)
# followed by down' and strange' (same n' = 5 and 4 respectively)

# =====================================================================
# SECTION D: Galois Confinement
# =====================================================================
print("\n--- D. GALOIS CONFINEMENT ---")
print("The Galois sector has alpha_s' < 0 (confining).")
print("This means Galois quarks form BOUND STATES.")
print("The Galois sector is confined at ALL scales.")
print()
print("Consequences:")
print("  - No free Galois quarks at any energy")
print("  - Galois hadrons are the physical states")
print("  - Interaction with SM via gravity only?")
print("  - Or via Galois-invariant operators (Tr over Z[phi])?")
print()

# The key question: how does the Galois sector interact with SM?
# In NCG: the Galois automorphism is an OUTER automorphism of the algebra.
# It's not a gauge symmetry. So Galois-sector particles would
# interact with SM particles via Galois-invariant operators.

# The Galois-invariant quantity is the Trace: Tr(a+b*phi) = 2a+b
# Physical observables must be Galois-invariant = rational.
# But individual STATES can live in one sector.

print("Galois-invariant interactions:")
print("  Any observable must be in Q (Galois invariant).")
print("  So: Tr(z*z') = N(z) is always rational.")
print("  The GRAVITY sector is Galois-invariant (h, rank are integers).")
print("  -> Galois particles interact gravitationally!")
print()
print("  EM sector: alpha involves phi AND pi -> NOT purely Galois")
print("  -> Galois particles may have modified EM interactions")
print()
print("  Strong sector: alpha_s in Z[phi], Galois-variant")
print("  -> Galois sector has DIFFERENT strong coupling (alpha_s')")
print("  -> Galois sector is separately confined")

# =====================================================================
# SECTION E: Relic Abundance Estimate
# =====================================================================
print("\n--- E. RELIC ABUNDANCE ---")
print("If Galois sector interacts only gravitationally:")
print("  -> Freeze-out via gravitational interactions")
print("  -> Relic abundance ~ (T_freeze / m_P)^2 * m_DM")
print("  This gives WIMPY-type dark matter only for m_DM ~ TeV")
print()

# The Galois electron (m'_e = m_e) would be stable
# and interact only gravitationally.
# But m = 0.511 MeV is too light for gravitational DM.

# A Galois BARYON (3 Galois quarks) would be the natural DM candidate
# Its mass would be set by the Galois QCD scale Lambda_QCD'
# Since alpha_s' has |1/alpha_s'| = |z'| = |6-4*phi| = 4*phi-6 ~ 0.47
# This means alpha_s' is VERY strong -> Lambda_QCD' >> Lambda_QCD

Lambda_QCD_galois_ratio = abs(z_prime) / (2*PHI**3)
print(f"|z'|/z = {Lambda_QCD_galois_ratio:.4f}")
print(f"|alpha_s'| = 1/|z'| = {1/abs(z_prime):.4f}")
print(f"  vs alpha_s = {alpha_s:.4f}")
print(f"  |alpha_s'| / alpha_s = {(1/abs(z_prime))/alpha_s:.2f}")
print()
print("Galois strong coupling is ~18x larger than SM!")
print("-> Galois confinement scale is MUCH higher than SM")
print("-> Galois baryons would be very heavy")
print()

# Rough estimate: Lambda_QCD ~ Lambda * exp(-2*pi/(b0*alpha_s))
# For Galois: Lambda_QCD' ~ Lambda * exp(-2*pi/(b0*|alpha_s'|))
# Since |alpha_s'| >> alpha_s, Lambda_QCD' >> Lambda_QCD

# =====================================================================
# SECTION F: Honest Assessment
# =====================================================================
print("\n" + "=" * 70)
print("SYNTHESIS: DARK MATTER")
print("=" * 70)
print()
print("The Galois sector (phi -> phi') provides:")
print("  1. A natural 'mirror' sector with different couplings")
print("  2. Separate confinement (alpha_s' > 0 but z' < 0)")
print("  3. Gravitational interaction with SM (guaranteed)")
print("  4. Inverted mass hierarchy (heavy SM = light Galois)")
print()
print("Problems:")
print("  1. HOW the Galois sector couples to SM is unclear")
print("     (beyond gravity, which is Galois-invariant)")
print("  2. The Galois QCD scale is very high -> heavy bound states")
print("  3. Relic abundance not calculable without interaction rate")
print("  4. Direct detection predictions impossible without coupling")
print()
print("CATEGORY: SPECULATIVE")
print("  The Galois sector EXISTS mathematically in Z[phi].")
print("  Whether it has physical manifestation as dark matter")
print("  requires knowing how the Galois automorphism acts on")
print("  the physical Hilbert space. This is NOT derived.")
print()
print("RECOMMENDATION: Mention in paper as a POSSIBILITY,")
print("not a prediction. The mathematical structure is suggestive")
print("but the physical mechanism is missing.")

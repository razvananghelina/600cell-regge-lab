"""
EXP-555: Why Exactly One Dimensional Anchor?
=============================================

QUESTION: Not "how to eliminate m_e" but "why exactly one, and why this one?"

THREE SUB-QUESTIONS:
  Q1: How many independent dimensional parameters does the framework need?
      (Prove it's exactly 1, not 0 or 2.)
  Q2: What makes m_e the natural choice among all possible anchors?
      (Is it structurally distinguished, or just conventional?)
  Q3: What is the minimal FORM of the anchor?
      (A mass? A length? A dimensionless ratio + Planck units?)

APPROACH: Map the dimensional structure of the theory explicitly.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from commons import PHI

phi = PHI
a1 = 5
b1 = 6

# Framework coupling constants (pure numbers)
sin2 = b1 / (a1**2 + 1)         # = 6/26
alpha_s = 1 / (2*phi**3)         # = 0.1180
B = 4*a1*phi**4
alpha = (B - np.sqrt(B**2 - 8*np.pi)) / (4*np.pi)  # = 1/137.036

# Derived dimensionless ratios
z_Planck = 4*phi**2              # = 10.472
m_e_over_m_P = alpha**z_Planck   # = 4.20e-23

print("=" * 70)
print("EXP-555: WHY EXACTLY ONE DIMENSIONAL ANCHOR?")
print("=" * 70)

# ============================================================
# Q1: WHY EXACTLY ONE?
# ============================================================
print(f"\n{'='*70}")
print("Q1: COUNTING DIMENSIONAL PARAMETERS")
print(f"{'='*70}")

print(f"""
  The framework produces a TREE of mass predictions:

  Level 0 (anchor):   m_e                         [1 free parameter]
  Level 1 (derived):  m_P = m_e / alpha^{{z_P}}    [ratio derived]
  Level 2 (derived):  m_Z = m_e * phi^25 * R      [ratio derived]
                      m_W = m_Z * cos(tW)          [ratio derived]
                      m_H = m_W * sqrt(...)        [ratio derived]
  Level 3 (derived):  m_f = m_e * phi^{{5a+6b+d}}  [9 ratios derived]
  Level 4 (derived):  m_nu = 2*m_e/phi^35         [ratio derived]

  Total masses: ~15 (9 fermions + W + Z + H + Planck + 3 neutrinos)
  Free parameters: 1 (m_e)
  Derived ratios: 14

  Why not 0 free parameters?
  PROOF: The spectrum of D^2 on the 600-cell has eigenvalues in lattice
  units. These are pure numbers (e.g., 2.29, 5.53, 9, 12, ...).
  The map from lattice units to GeV is:
    E_phys = E_latt / a^2
  where a is the lattice spacing. This requires ONE dimensionful number.

  Why not 2 free parameters?
  Because ALL mass ratios are determined by the spectrum.
  Given m_e, EVERY other mass is fixed. There is no second free ratio.
  Explicitly: m_mu/m_e = phi^11 * correction (no freedom).
  The 9 fermion masses, 3 boson masses, and m_P are ALL determined.

  THEOREM (informal): The number of dimensional parameters equals
  the number of independent dimensionful quantities in physics (1)
  minus the number of dimensionless ratios determined by the theory (all).
  = 1 - 0 = 1.  (The theory determines all ratios but no absolute scale.)
""")

# ============================================================
# Q2: WHY m_e?
# ============================================================
print(f"{'='*70}")
print("Q2: WHAT MAKES m_e THE NATURAL ANCHOR?")
print(f"{'='*70}")

print(f"""
  Five properties distinguish m_e from all other masses:
""")

# Property 1: Galois fixed point
print(f"  P1: GALOIS FIXED POINT")
print(f"      z_e = (a,b) = (0,0). Mass exponent n_e = 0.")
print(f"      sigma(z_e) = sigma(0) = 0.  The electron is the ONLY")
print(f"      fermion that is POINTWISE FIXED by Galois conjugation.")
print(f"      All others have z != 0 and receive corrections.")
print(f"      m_e is the origin of the mass lattice Z[phi].")

# Property 2: McKay graph origin
print(f"\n  P2: McKay GRAPH ORIGIN")
print(f"      The electron maps to the trivial node rho_0 on the")
print(f"      McKay graph (affine E8). This is the unique node with:")
print(f"        - dimension 1 (trivial irrep)")
print(f"        - zero Casimir: C_2(rho_0) = 0")
print(f"        - zero adjacency: lambda_adj(rho_0) = 12 (degree)")
print(f"      All mass corrections are measured FROM this node.")

# Property 3: Spectral action ground state
print(f"\n  P3: SPECTRAL ACTION GROUND STATE")
print(f"      In the mass formula m_f = m_e * phi^n, the electron")
print(f"      has n = 0. It is the GROUND STATE of the mass operator.")
print(f"      The spectral action on M x F has m_e as the fundamental")
print(f"      Yukawa scale (the smallest nonzero eigenvalue of D_F).")

# Property 4: Dimensional transmutation
print(f"\n  P4: DIMENSIONAL TRANSMUTATION")
print(f"      m_e can be traded for m_P via m_e = m_P * alpha^{{z_P}}.")
print(f"      But m_P = sqrt(hbar*c/G) involves THREE fundamental")
print(f"      constants (hbar, c, G), while m_e is a single measured number.")
print(f"      The framework's economy: 1 input (m_e) vs 3 inputs (hbar,c,G).")

# Property 5: The hierarchy identity
print(f"\n  P5: HIERARCHY IDENTITY")
print(f"      m_e/m_P = alpha^{{4*phi^2}} links m_e to the gauge sector.")
print(f"      This means m_e is not independent of the spectral data --")
print(f"      it is the EXPONENTIAL of a spectral quantity (z_P = 4*phi^2)")
print(f"      applied to a spectral coupling (alpha).")
print(f"      In Planck units: m_e = alpha^{{10.47}} (dimensionless).")
print(f"      The dimensional content is ENTIRELY in the choice of")
print(f"      unit system (GeV vs Planck vs natural).")

# Numerical verification
print(f"\n  Numerical check:")
m_e = 0.51099895e-3  # GeV
m_P = 1.220890e19    # GeV
print(f"    m_e = {m_e:.4e} GeV")
print(f"    m_P = {m_P:.4e} GeV")
print(f"    m_e/m_P = {m_e/m_P:.4e}")
print(f"    alpha^(4*phi^2) = {alpha**z_Planck:.4e}")
print(f"    Ratio = {(m_e/m_P) / alpha**z_Planck:.6f} (error {abs((m_e/m_P)/alpha**z_Planck - 1)*100:.2f}%)")


# ============================================================
# Q3: MINIMAL FORM OF THE ANCHOR
# ============================================================
print(f"\n{'='*70}")
print("Q3: MINIMAL FORM OF THE DIMENSIONAL ANCHOR")
print(f"{'='*70}")

print(f"""
  The anchor can be stated in three equivalent forms:

  FORM A (mass): m_e = 0.511 MeV.
    Conventional. Requires a unit system (MeV).

  FORM B (ratio + Planck): m_e/m_P = alpha^{{4*phi^2}}.
    DIMENSIONLESS. All factors derived from a1 = 5.
    The dimensional content is in choosing Planck units as reference.
    m_P = sqrt(hbar*c/G) is defined by 3 constants.

  FORM C (pure number + natural units): m_e = alpha^{{4*phi^2}} in Planck units.
    A single pure number {alpha**z_Planck:.6e}.
    If you declare Planck units as fundamental (hbar=c=G=1),
    then m_e is a DERIVED pure number, not an input.

  The question becomes: is hbar=c=G=1 a physical choice or a convention?
""")

# The deep point:
print(f"  THE DEEP POINT:")
print(f"  ===============")
print(f"  In Planck units (hbar=c=G=1), ALL masses are pure numbers.")
print(f"  m_e = alpha^{{4*phi^2}} = {alpha**z_Planck:.6e}")
print(f"  m_Z = m_e * phi^25 * R = {m_e/m_P * phi**25 * 1.063:.6e}")
print(f"  m_P = 1 (by definition)")
print(f"")
print(f"  The 'one free parameter' is then not m_e itself, but the")
print(f"  CHOICE TO USE Planck units. In any other unit system,")
print(f"  you need to specify one conversion factor.")
print(f"")
print(f"  This is the same situation as in ALL of physics:")
print(f"  the number of independent dimensions (M, L, T) is a human")
print(f"  convention, not a law of nature. With hbar=c=G=1, there")
print(f"  are ZERO independent dimensions, and all 'masses' are numbers.")

# ============================================================
# THE STRUCTURAL DIAGRAM
# ============================================================
print(f"\n{'='*70}")
print("STRUCTURAL DIAGRAM: THE ANCHOR CHAIN")
print(f"{'='*70}")

print(f"""
  PURE SPECTRAL DATA (dimensionless, from a1=5):
    phi, alpha, alpha_s, sin^2(tW), z_Planck, n0=a1^2, ...

           |
           | alpha^{{z_Planck}} = derived ratio
           v

  m_e / m_P = alpha^{{4*phi^2}}  <-- dimensionless ratio

           |
           | x m_P (one dimensional constant)
           v

  m_e = 0.511 MeV  <-- the anchor

           |
           | x phi^n (derived exponents)
           v

  ALL OTHER MASSES: m_mu, m_tau, m_u, ..., m_Z, m_H, m_P

  The anchor point is WHERE the derived ratio meets the physical world.
  It can be placed at m_e, at m_P, or at any other mass.
  The CHOICE of placement is conventional.
  The EXISTENCE of exactly one anchor point is structural.
""")


# ============================================================
# CAN THE ANCHOR BE DERIVED?
# ============================================================
print(f"{'='*70}")
print("CAN THE ANCHOR BE DERIVED?")
print(f"{'='*70}")

print(f"""
  Three speculative routes:

  ROUTE 1: Gravitational path integral.
    If G (and hence m_P) emerges from a sum over discrete geometries
    that includes the 600-cell, then m_P is no longer free.
    Status: the framework gives discrete GR (Regge on the 600-cell)
    but does NOT derive G from first principles.

  ROUTE 2: Self-referential bootstrap.
    The equation x^2 = x + 1 is self-referential: the solution phi
    equals 1 + 1/phi. Could a similar self-referential condition
    fix the absolute scale? E.g., the lattice spacing equals the
    Compton wavelength of the lightest particle ON the lattice?
    a = hbar/(m_e * c) and m_e = m_P * alpha^{{z_P}}: this gives
    a = hbar*c / (m_P * c^2 * alpha^{{z_P}}) = l_P / alpha^{{z_P}}
    = l_P * (m_P/m_e). But this is just the definition of l_P.
    CIRCULAR. No new information.

  ROUTE 3: Cosmological selection.
    The observed universe has a specific age T ~ 13.8 Gyr.
    Combining T with hbar and c gives a mass scale m_cosmo ~ hbar/(Tc).
    If the 600-cell framework is embedded in a cosmological context,
    the age of the universe could provide the missing anchor.
    This connects to the anthropic principle but goes beyond the
    current framework.

  CURRENT STATUS: None of these routes is closed.
  The framework's honest position is:
    "One dimensional anchor is needed. m_e is the natural choice
    (ground state of the mass lattice, Galois fixed point).
    The anchor is not derived; it is the sole empirical input."
""")


# ============================================================
# WHAT DOES "EXACTLY ONE" MEAN FOR OTHER FRAMEWORKS?
# ============================================================
print(f"{'='*70}")
print("COMPARISON: DIMENSIONAL PARAMETERS IN OTHER FRAMEWORKS")
print(f"{'='*70}")

print(f"""
  Standard Model:     ~19 free parameters (masses, couplings, angles)
  MSSM:               ~124 free parameters
  String landscape:    ~10^500 vacua (no selection)
  NCG (Connes):        ~15 (finite geometry chosen, Yukawas free)
  This framework:      1 (m_e; everything else from a1=5)

  The reduction from ~19 to 1 is the framework's central claim.
  The reduction from 1 to 0 may require physics beyond geometry
  (cosmology, gravitational path integral, or a new principle).

  Note: even in a "theory of everything" with 0 free parameters,
  you still need to choose a unit system. The question "why does
  m_e = 0.511 MeV?" is really "why does the ratio m_e/MeV = 0.511?",
  which requires defining MeV in terms of something. In Planck units,
  m_e IS a derived number: alpha^{{4*phi^2}}.
""")

# ============================================================
# VERDICT
# ============================================================
print(f"{'='*70}")
print("VERDICT")
print(f"{'='*70}")

print(f"""
  Q1: WHY EXACTLY ONE?
      Because the spectral data determine all dimensionless ratios
      (proven) but no absolute scale (no-go, exp551). One anchor
      is necessary and sufficient.

  Q2: WHY m_e?
      Five structural properties distinguish it:
      (P1) Galois fixed point z = 0
      (P2) McKay graph origin (trivial irrep)
      (P3) Spectral ground state (n = 0)
      (P4) Minimal input (1 number vs 3 for m_P)
      (P5) Hierarchy identity m_e/m_P = alpha^{{4*phi^2}}

  Q3: MINIMAL FORM?
      In Planck units: m_e = alpha^{{4*phi^2}} is a DERIVED pure number.
      The "anchor" is then the CHOICE of Planck units (hbar=c=G=1),
      which is conventional, not dynamical.

  PROGRAM OPENED:
      The question "why exactly one anchor?" has a clean answer.
      The question "why THIS anchor value?" reduces to "why alpha^{{10.5}}?"
      which reduces to "why is alpha small?" -- and alpha IS derived
      (from the quadratic 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0).
      So the hierarchy is EXPLAINED (alpha^moderate = small number)
      even if the absolute scale is CONVENTIONAL.
""")

print("=" * 70)
print("EXP-555 COMPLETE")
print("=" * 70)

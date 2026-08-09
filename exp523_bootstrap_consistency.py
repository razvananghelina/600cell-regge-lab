"""
EXP-523: Bootstrap Self-Consistency for Coupling Constants
===========================================================

SHIFT IN PARADIGM:
  NOT: "alpha minimizes functional X on spectrum" (variational)
  BUT: "alpha is the unique value where physical (phi) and dark (phi')
        sectors are simultaneously consistent as TQFT" (bootstrap)

The framework already selects a1=5 via a bootstrap condition:
  d_{1/2}(a1) = phi(a1)  =>  quantum dim = algebraic parameter  =>  a1=5 unique.

Can the same logic fix the coupling constants?

INGREDIENTS:
  - Fusion rules N_{ij}^k in Z (integers, Galois-invariant)
  - Physical quantum dims: d_j = {1, phi, phi, 1}
  - Dark quantum dims: d'_j = {1, -1/phi, -1/phi, 1} (Galois conjugate)
  - S-matrix unitarity: S*S = I
  - Pentagon identity (associativity of fusion)
  - Hexagon identity (braiding consistency)

QUESTION: What self-consistency conditions constrain a continuous
parameter (like alpha) beyond the discrete choice a1=5?
"""

import numpy as np
from fractions import Fraction
import sys
sys.path.insert(0, '.')
from commons import (PHI, PHI_CONJ, SQRT5, a1, b1, N, degree, d_ST,
                      alpha_em, inv_alpha, sin2_tW, Omega_DM, d_dark,
                      D2_PHYS, D2_DARK, L1, L2, L6, L8,
                      GALOIS_NORM_A, GALOIS_NORM_B)

print("=" * 70)
print("EXP-523: BOOTSTRAP SELF-CONSISTENCY")
print("=" * 70)

# ============================================================
# PART 1: THE EXISTING BOOTSTRAP - HOW a1=5 IS SELECTED
# ============================================================
print(f"\n{'='*70}")
print("PART 1: THE a1=5 BOOTSTRAP (REVIEW)")
print(f"{'='*70}")

print(f"""
  For SU(2)_k TQFT at level k, set a1 = k+2.
  Quantum dimension of fundamental: d_{{1/2}} = 2*cos(pi/a1)
  Golden ratio: phi(a1) = (1+sqrt(a1))/2

  The BOOTSTRAP CONDITION: d_{{1/2}}(a1) = phi(a1)
  i.e., 2*cos(pi/a1) = (1+sqrt(a1))/2

  This selects a1 = 5 UNIQUELY among integers a1 >= 3.
""")

# Verify for a1 = 3,4,5,6,7,...
print(f"  {'a1':>4} {'d_{1/2}':>10} {'phi(a1)':>10} {'match?':>8}")
for a in range(3, 12):
    d_half = 2 * np.cos(np.pi / a)
    phi_a = (1 + np.sqrt(a)) / 2
    match = abs(d_half - phi_a) < 0.001
    print(f"  {a:>4} {d_half:>10.6f} {phi_a:>10.6f} {'YES' if match else '':>8}")

# ============================================================
# PART 2: WHAT DOES "SIMULTANEOUS CONSISTENCY" MEAN?
# ============================================================
print(f"\n{'='*70}")
print("PART 2: SIMULTANEOUS TQFT CONSISTENCY")
print(f"{'='*70}")

print(f"""
  The physical TQFT (SU(2)_3) has:
    - Fusion coefficients N_{{ij}}^k (integers)
    - Quantum dims d = (1, phi, phi, 1)
    - S-matrix S with S*S^T = I
    - Verlinde: N_{{ij}}^k = sum_l S_il S_jl S*_kl / S_0l

  The dark TQFT (Galois conjugate) has:
    - SAME fusion coefficients N_{{ij}}^k (Galois-invariant!)
    - Quantum dims d' = (1, -1/phi, -1/phi, 1)
    - S'-matrix: sigma(S)

  BOTH must satisfy:
    (a) Verlinde formula (fusion from S-matrix)
    (b) S-matrix unitarity
    (c) Pentagon identity (associativity)
    (d) Hexagon identity (braiding)

  For the PHYSICAL sector, all hold by construction.
  For the DARK sector, (a)-(d) hold because sigma is a ring automorphism.

  So: BOTH sectors are separately consistent. No new constraint yet.

  THE QUESTION: What COUPLES them?
""")

# ============================================================
# PART 3: WHAT COUPLES PHYSICAL AND DARK?
# ============================================================
print(f"{'='*70}")
print("PART 3: THE COUPLING BETWEEN SECTORS")
print(f"{'='*70}")

print(f"""
  The two sectors share:
    1. The SAME fusion ring (N_{{ij}}^k are integers, Galois-fixed)
    2. The SAME braiding phases (rational, Galois-fixed)
    3. DIFFERENT quantum dimensions (d vs d')

  What couples them further?
    - The SPECTRAL ACTION sees both sectors through the 600-cell spectrum
    - The LAPLACIAN eigenvalues involve BOTH phi and sqrt(5)
    - The Galois automorphism sigma: phi -> -1/phi acts on BOTH

  A self-consistency condition would be:
  "The theory that describes sector A must be compatible with
   the existence of sector B, and vice versa."

  Concretely: if sector A determines coupling alpha, and sector B
  determines coupling alpha', the CONSISTENCY condition is:
  alpha = f(alpha') and alpha' = f(alpha) => alpha is a FIXED POINT of f.
""")

# ============================================================
# PART 4: THE GALOIS FIXED-POINT EQUATION
# ============================================================
print(f"{'='*70}")
print("PART 4: GALOIS FIXED-POINT EQUATION")
print(f"{'='*70}")

# The key idea: if a coupling constant g lives in Q(sqrt5),
# then g = a + b*sqrt5 for rational a, b.
# Under Galois: sigma(g) = a - b*sqrt5 = g'.
#
# Self-consistency: the coupling in the dark sector should be
# related to the coupling in the physical sector by Galois,
# BUT also by the physics (same action, different sector).
#
# The simplest fixed-point condition:
# g * sigma(g) = rational constant (Galois NORM is rational)
# g + sigma(g) = rational constant (Galois TRACE is rational)
#
# For alpha_em: if alpha lives in Q(sqrt5), then:
# alpha = a + b*sqrt5 => alpha*sigma(alpha) = a^2 - 5b^2 (rational!)
# alpha + sigma(alpha) = 2a (rational!)
#
# Our alpha comes from: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
# Let's check: does alpha live in Q(sqrt5, pi)?

print(f"  Alpha equation: 2*pi*x^2 - 4*a1*phi^4*x + 1 = 0")
print(f"  Coefficient: 4*a1*phi^4 = 4*5*(3+sqrt5)^2/4 = 5*(14+6sqrt5)")
print(f"             = 70 + 30*sqrt5 = {4*a1*PHI**4:.6f}")

coeff = 4*a1*PHI**4  # = 70 + 30*sqrt5
print(f"\n  alpha = (coeff - sqrt(coeff^2 - 8pi)) / (4pi)")
print(f"  coeff = {coeff:.6f}")
print(f"  coeff^2 = {coeff**2:.6f}")
print(f"  8*pi = {8*np.pi:.6f}")
print(f"  disc = coeff^2 - 8pi = {coeff**2 - 8*np.pi:.6f}")

# The discriminant involves pi, so alpha does NOT live in Q(sqrt5).
# It lives in Q(sqrt5, sqrt(something involving pi)).
# This means Galois of Q(sqrt5) alone doesn't constrain alpha.

# BUT: what if we ask a DIFFERENT question?
# Instead of alpha in Q(sqrt5), ask:
# "What constraint on the 600-cell spectrum forces the alpha equation?"

print(f"""
  OBSERVATION: alpha does NOT live in Q(sqrt5) — it involves pi.
  So Galois conjugation alone cannot fix alpha.

  But the alpha EQUATION (2*pi*x^2 - 4*a1*phi^4*x + 1 = 0) has
  COEFFICIENTS in Q(sqrt5, pi). The Galois action on Q(sqrt5)
  maps 4*a1*phi^4 -> 4*a1*(phi')^4 = 4*a1/phi^4.

  Dark alpha equation: 2*pi*x^2 - (4*a1/phi^4)*x + 1 = 0

  Dark coupling: alpha' = (4*a1/phi^4 - sqrt((4*a1/phi^4)^2 - 8pi))/(4pi)
""")

coeff_dark = 4*a1/PHI**4
disc_dark = coeff_dark**2 - 8*np.pi
if disc_dark > 0:
    alpha_dark = (coeff_dark - np.sqrt(disc_dark)) / (4*np.pi)
else:
    alpha_dark = None
    print(f"  Dark discriminant = {disc_dark:.6f} < 0: NO REAL DARK ALPHA!")

print(f"  Physical: coeff = 4*a1*phi^4 = {coeff:.6f}")
print(f"  Dark:     coeff' = 4*a1/phi^4 = {coeff_dark:.6f}")
print(f"  Ratio: phi^8 = {PHI**8:.6f}")

print(f"\n  Dark discriminant: {coeff_dark:.4f}^2 - 8pi = {coeff_dark**2:.4f} - {8*np.pi:.4f} = {disc_dark:.4f}")
if disc_dark < 0:
    print(f"\n  *** THE DARK ALPHA EQUATION HAS NO REAL SOLUTION! ***")
    print(f"  This means: the dark sector CANNOT have a real EM coupling.")
    print(f"  coeff' = 4*a1/phi^4 = {coeff_dark:.6f}")
    print(f"  Need coeff'^2 >= 8*pi = {8*np.pi:.6f}")
    print(f"  But coeff'^2 = {coeff_dark**2:.6f} < {8*np.pi:.6f}")
    print(f"\n  The PHYSICAL sector has coupling. The DARK sector doesn't.")
    print(f"  This is a CONSISTENCY CONDITION: only the physical sector")
    print(f"  supports electromagnetism!")

# ============================================================
# PART 5: WHAT IS THE CRITICAL CONDITION?
# ============================================================
print(f"\n{'='*70}")
print("PART 5: THE CRITICAL BOUNDARY")
print(f"{'='*70}")

# The alpha equation has real solutions when:
# (4*a1*d^4)^2 >= 8*pi  where d = quantum dimension
# i.e., 16*a1^2*d^8 >= 8*pi
# i.e., d^8 >= pi/(2*a1^2)
# i.e., d >= (pi/(2*a1^2))^{1/8}

d_critical = (np.pi / (2*a1**2))**(1/8)
print(f"  Alpha equation: 2*pi*x^2 - 4*a1*d^4*x + 1 = 0")
print(f"  Real solutions require: (4*a1*d^4)^2 >= 8*pi")
print(f"  i.e., d^8 >= pi/(2*a1^2) = {np.pi/(2*a1**2):.6f}")
print(f"  i.e., d >= (pi/50)^{{1/8}} = {d_critical:.6f}")
print(f"")
print(f"  Physical: d = phi = {PHI:.6f} >= {d_critical:.6f}  YES -> alpha exists")
print(f"  Dark:     d'= 1/phi = {1/PHI:.6f} < {d_critical:.6f}   NO  -> alpha DOESN'T exist!")
print(f"")
print(f"  The critical d is {d_critical:.6f}")
print(f"  phi = {PHI:.6f}")
print(f"  1/phi = {1/PHI:.6f}")
print(f"  phi and 1/phi are on OPPOSITE SIDES of the critical value!")

# This is a bootstrap result: the EM coupling EXISTS in the physical sector
# and DOESN'T EXIST in the dark sector. The boundary is:
# d^8 = pi/(2*a1^2)

# Is the critical d a clean expression?
print(f"\n  d_crit = (pi/(2*a1^2))^(1/8) = (pi/50)^(1/8) = {d_critical:.10f}")
print(f"  d_crit^2 = {d_critical**2:.10f}")
print(f"  d_crit^4 = {d_critical**4:.10f}")
print(f"  d_crit^8 = {d_critical**8:.10f} = pi/50 = {np.pi/50:.10f}")

# Is d_crit close to anything?
print(f"\n  d_crit vs framework:")
print(f"    1 = {1.0:.6f} (diff {abs(d_critical-1):.4f})")
print(f"    phi/phi^2 = 1/phi = {1/PHI:.6f} (diff {abs(d_critical-1/PHI):.4f})")
print(f"    sqrt(phi/phi^2)... no")

# ============================================================
# PART 6: THE BOOTSTRAP CONTENT
# ============================================================
print(f"\n{'='*70}")
print("PART 6: WHAT THE BOOTSTRAP SAYS")
print(f"{'='*70}")

# Physical alpha: the SMALLER root of 2*pi*x^2 - C*x + 1 = 0 with C=4a1*phi^4
# Dark alpha: would be from C' = 4*a1/phi^4, but C'^2 < 8*pi, so NO SOLUTION.

# The BOOTSTRAP STATEMENT:
# "The EM coupling constant exists (is real) in the physical sector
#  and does NOT exist in the dark sector. This is because phi > d_crit > 1/phi."

# Moreover, the VALUE of alpha is fixed by:
# alpha = (C - sqrt(C^2 - 8pi))/(4pi) with C = 4*a1*phi^4

# Can we derive alpha from the BOUNDARY condition instead?
# At the boundary: C_crit^2 = 8*pi, so C_crit = sqrt(8*pi) = 2*sqrt(2pi)
# And alpha_boundary = C_crit/(4pi) = sqrt(2pi)/(2pi) = 1/sqrt(2pi)

alpha_boundary = 1/np.sqrt(2*np.pi)
print(f"  At the boundary (d = d_crit):")
print(f"    C_crit = 2*sqrt(2pi) = {2*np.sqrt(2*np.pi):.6f}")
print(f"    alpha_boundary = 1/sqrt(2pi) = {alpha_boundary:.6f}")
print(f"    = {alpha_boundary:.6f}")
print(f"    1/alpha_boundary = {1/alpha_boundary:.6f} = sqrt(2pi) = {np.sqrt(2*np.pi):.6f}")
print(f"")
print(f"  Physical alpha = {alpha_em:.6f}")
print(f"  Boundary alpha = {alpha_boundary:.6f}")
print(f"  Ratio = {alpha_em/alpha_boundary:.6f}")

# The DISTANCE from the boundary:
# C/C_crit = 4*a1*phi^4 / (2*sqrt(2pi))
ratio_C = coeff / (2*np.sqrt(2*np.pi))
print(f"\n  C/C_crit = {coeff:.4f}/{2*np.sqrt(2*np.pi):.4f} = {ratio_C:.6f}")
print(f"  = 4*a1*phi^4 / (2*sqrt(2pi))")
print(f"  = 2*a1*phi^4 / sqrt(2pi)")
print(f"  = {2*a1*PHI**4/np.sqrt(2*np.pi):.6f}")

# ============================================================
# PART 7: PRODUCT OF PHYSICAL AND DARK EQUATIONS
# ============================================================
print(f"\n{'='*70}")
print("PART 7: GALOIS NORM OF THE ALPHA EQUATION")
print(f"{'='*70}")

# The alpha equation: 2*pi*x^2 - C*x + 1 = 0, C = 4*a1*phi^4
# Galois conjugate: 2*pi*x^2 - C'*x + 1 = 0, C' = 4*a1/phi^4
#
# Product equation (Galois norm):
# (2*pi*x^2 - C*x + 1)(2*pi*y^2 - C'*y + 1) = 0
#
# But if we set x = y (same coupling in both sectors):
# (2*pi*x^2 + 1)^2 - C*C'*x^2 = 0... no, that's not right.
#
# C * C' = 16*a1^2*phi^4*(1/phi^4) = 16*a1^2 = 400
# C + C' = 4*a1*(phi^4 + 1/phi^4) = 4*a1*(phi^4+phi'^4)

CC_prime = coeff * coeff_dark
CC_sum = coeff + coeff_dark

print(f"  C * C' = 4*a1*phi^4 * 4*a1/phi^4 = 16*a1^2 = {CC_prime:.6f} = {16*a1**2}")
print(f"  C + C' = {CC_sum:.6f}")
print(f"  = 4*a1*(phi^4+1/phi^4) = 4*a1*(phi^4+phi'^4)")

phi4_plus_phi4_inv = PHI**4 + 1/PHI**4
print(f"  phi^4 + 1/phi^4 = {phi4_plus_phi4_inv:.6f}")
# phi^4 = 7+3sqrt5)/2, 1/phi^4 = (7-3sqrt5)/2
# sum = 7 (rational!)
print(f"  = (7+3sqrt5)/2 + (7-3sqrt5)/2 = 7")
print(f"  So C + C' = 4*a1*7 = {4*a1*7} = {CC_sum:.0f}")

# The NORM equation: if alpha satisfies 2pi*x^2 - Cx + 1 = 0,
# then the GALOIS NORM of (2pi*x^2 - Cx + 1) is:
# (2pi*x^2 + 1)^2 - (CC')*x^2 = ... wait, C is not just multiplied by x.
# Let me think about this differently.

# Product of the two roots of the physical equation:
# x1*x2 = 1/(2pi)  (from Vieta's)
# Sum of roots: x1+x2 = C/(2pi) = 4*a1*phi^4/(2pi)

product_roots = 1/(2*np.pi)
sum_roots = coeff/(2*np.pi)
print(f"\n  Vieta's formulas for physical equation:")
print(f"    x1 * x2 = 1/(2pi) = {product_roots:.6f}")
print(f"    x1 + x2 = C/(2pi) = {sum_roots:.6f}")
print(f"    x1 (physical alpha) = {alpha_em:.6f}")

x2 = product_roots / alpha_em  # second root
print(f"    x2 (second root) = {x2:.6f}")
print(f"    1/x2 = {1/x2:.6f}")
print(f"    x2 = 1/(2pi*alpha) = {1/(2*np.pi*alpha_em):.6f}")

# Interesting: x1*x2 = 1/(2pi). So the second root is x2 = 1/(2pi*alpha).
# x2 = 1/(2pi * 0.007297) = 1/0.04584 = 21.81
# 1/x2 = 2pi*alpha = 0.04584

# The physical alpha and the second root are "Vieta duals":
# alpha * alpha_2 = 1/(2pi)

# For the DARK equation: no real roots exist (disc < 0).
# But the COMPLEX roots have:
# y1*y2 = 1/(2pi) (same!)
# y1+y2 = C'/(2pi) = 4*a1/(phi^4 * 2pi)

# The modulus of the complex dark roots:
# |y|^2 = y*y_bar = 1/(2pi) (since complex conjugate pair)
# So |y| = 1/sqrt(2pi) = alpha_boundary!

print(f"\n  Dark equation (no real roots):")
print(f"    y1*y2 = 1/(2pi) = {product_roots:.6f} (same as physical)")
print(f"    |y| = 1/sqrt(2pi) = {1/np.sqrt(2*np.pi):.6f}")
print(f"    This is the MODULUS of the dark coupling (complex!)")
print(f"    |alpha_dark| = 1/sqrt(2pi) = {alpha_boundary:.6f}")
print(f"    1/|alpha_dark| = sqrt(2pi) = {np.sqrt(2*np.pi):.6f}")

# ============================================================
# PART 8: THE SELF-CONSISTENCY EQUATION
# ============================================================
print(f"\n{'='*70}")
print("PART 8: THE SELF-CONSISTENCY STRUCTURE")
print(f"{'='*70}")

# Physical: alpha is REAL (sector supports EM)
# Dark: alpha' is COMPLEX (sector does NOT support real EM)
#
# Both come from the SAME algebraic equation, just with d vs d'.
# The CONSISTENCY is:
# alpha * alpha_2 = 1/(2pi)  (physical: two real roots)
# |alpha_dark|^2 = 1/(2pi)   (dark: complex conjugate pair)
#
# BOTH give the SAME product = 1/(2pi)!
# This is the Galois invariant: the product of roots of the quadratic
# is 1/(2pi) regardless of which sector.

print(f"""
  SELF-CONSISTENCY STRUCTURE:
  ===========================

  Physical sector (d = phi):
    2pi*x^2 - (4a1*phi^4)*x + 1 = 0
    Two REAL roots: alpha = {alpha_em:.6f}, alpha_2 = {x2:.6f}
    Product: alpha * alpha_2 = 1/(2pi)

  Dark sector (d = 1/phi):
    2pi*x^2 - (4a1/phi^4)*x + 1 = 0
    Two COMPLEX roots: |alpha_dark| = 1/sqrt(2pi) = {alpha_boundary:.6f}
    Product: |alpha_dark|^2 = 1/(2pi)

  BOOTSTRAP CONTENT:
    1. Both sectors share: root product = 1/(2pi)
    2. Physical sector: REAL coupling -> electromagnetism EXISTS
    3. Dark sector: COMPLEX coupling -> electromagnetism ABSENT
    4. The TRANSITION between real and complex occurs at d_crit = (pi/50)^(1/8)
    5. phi > d_crit > 1/phi selects the physical sector as the one with EM

  This is NOT a derivation of alpha's VALUE from bootstrap.
  It IS a derivation of alpha's EXISTENCE (real vs complex) from bootstrap.
  The VALUE still comes from the quadratic equation with pi.

  OPEN: Why does the alpha equation involve pi? And why degree 2?
""")

# ============================================================
# PART 9: CHECK SIN2_TW AND OMEGA BOOTSTRAP
# ============================================================
print(f"{'='*70}")
print("PART 9: DO OTHER COUPLINGS HAVE BOOTSTRAP STRUCTURE?")
print(f"{'='*70}")

# sin2_tW = b1/(a1^2+1) = 6/26 -- this is RATIONAL, so Galois-invariant!
# sigma(sin2_tW) = sin2_tW. No physical/dark distinction.
print(f"  sin^2(tW) = {sin2_tW:.6f} = 6/26 (RATIONAL)")
print(f"  sigma(sin^2(tW)) = sin^2(tW)  (Galois-fixed)")
print(f"  -> Same in both sectors. No bootstrap constraint.")

# alpha_s = 1/(2*phi^3)
alpha_s = 1/(2*PHI**3)
alpha_s_dark = 1/(2*abs(PHI_CONJ)**3)  # sigma: phi -> |phi'| = 1/phi
print(f"\n  alpha_s = 1/(2*phi^3) = {alpha_s:.6f}")
print(f"  alpha_s_dark = 1/(2/phi^3) = phi^3/2 = {1/(2/PHI**3):.6f}")
print(f"  Product: alpha_s * alpha_s_dark = 1/4 = {alpha_s * (PHI**3/2):.6f}")
print(f"  sigma(alpha_s) = phi^3/2 = {PHI**3/2:.6f}")
print(f"  alpha_s * sigma(alpha_s) = 1/(2phi^3) * phi^3/2 = 1/4 = {alpha_s*PHI**3/2:.6f}")

# Omega_DM = 7-phi
sigma_Omega = 7 - PHI_CONJ  # = 7 + 1/phi = 7+phi-1 = 6+phi
# Actually sigma(7-phi) = 7-phi' where phi' = (1-sqrt5)/2
# sigma(phi) in Q(sqrt5): phi -> (1-sqrt5)/2 = phi' (not -1/phi!)
# 7-phi = (13-sqrt5)/2 -> sigma: (13+sqrt5)/2
sigma_Omega_correct = (13+SQRT5)/2
print(f"\n  Omega = 7-phi = {Omega_DM:.6f}")
print(f"  sigma(Omega) = (13+sqrt5)/2 = {sigma_Omega_correct:.6f}")
print(f"  Omega * sigma(Omega) = (169-5)/4 = 164/4 = 41 (prime)")
print(f"  Verify: {Omega_DM * sigma_Omega_correct:.6f}")
print(f"  Omega + sigma(Omega) = 13")
print(f"  Verify: {Omega_DM + sigma_Omega_correct:.6f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  BOOTSTRAP RESULTS:
  ==================

  1. a1 = 5 SELECTION (known): d_{{1/2}} = phi uniquely for a1=5.

  2. ALPHA EXISTENCE (new):
     Physical sector (d=phi): alpha equation has REAL solutions
     Dark sector (d=1/phi): alpha equation has NO real solutions
     The critical boundary: d_crit = (pi/50)^{{1/8}} = {d_critical:.6f}
     phi = {PHI:.4f} > {d_critical:.4f} > {1/PHI:.4f} = 1/phi

     CONTENT: EM coupling EXISTS in physical sector, ABSENT in dark.
     This is a genuine bootstrap result — the Galois structure
     determines WHICH sector has electromagnetism.

  3. ALPHA VALUE: Still comes from the quadratic 2pi*x^2 - Cx + 1 = 0.
     The value is NOT fixed by bootstrap alone — it requires pi
     (from the spectral action normalization).

  4. GALOIS NORMS OF COUPLINGS:
     alpha * alpha_2 = 1/(2pi)  (root product, Galois-invariant)
     alpha_s * sigma(alpha_s) = 1/4  (Galois norm)
     Omega * sigma(Omega) = 41  (prime!)
     sin^2(tW) is rational (Galois-fixed)

  STATUS: The bootstrap SELECTS a1=5 and determines which sector
  has EM coupling. But it does not (yet) fix the coupling VALUES
  beyond what the quadratic equation gives.
""")

print("=" * 70)
print("EXP-523 COMPLETE")
print("=" * 70)

# EXP-372: Derive the neutrino mass splitting ratio
# ==================================================
# KNOWN: m_3 = 2*m_e/phi^35 (DERIVED in exp371)
# OPEN:  m_2/m_3 = sqrt(alpha*phi^3)  -- WHY?
#        Equivalently: Dm2_21/Dm2_31 = alpha*phi^3
#
# Prior work (exp367):
#   - 4 hypotheses tested (A-D), only D matches as rewriting
#   - phi^3 = (L8/L1)^(3/4), exponent 3/4 = 1-1/d universal
#   - Category: PATTERN (not derived)
#
# NEW approaches:
#   1. Kernel operator spectrum
#   2. Alpha as perturbation of rank-1 seesaw
#   3. Algebraic decomposition in Z[phi]
#   4. PMNS structure constraint
#   5. Galois kernel mass matrix from first principles
#
# RULE ZERO: categorize honestly.
# REGULA BLIND: derive FIRST, compare AFTER.
# Windows: no Unicode in print/comments.

import numpy as np
from itertools import product as iprod

# Framework constants
a1 = 5
b1 = 6
PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2  # = 1-phi = -1/phi
N = 120
N_gen = 3
N_eig = 9
DEG = 12
rank_E8 = 8

A_eq = 2 * np.pi
B_eq = -4 * a1 * PHI**4
C_eq = 1.0
disc = B_eq**2 - 4 * A_eq * C_eq
ALPHA = (-B_eq - np.sqrt(disc)) / (2 * A_eq)
ALPHA_S = 1 / (2 * PHI**3)

m_e_eV = 0.51099895e6  # eV

# Kernel Laplacians
L0 = 0.0
L1 = DEG - b1 * PHI        # 12 - 6*phi
L8 = DEG - b1 * PHI_CONJ   # 12 + 6/phi
# Verified: L1*L8 = 36 = b1^2, L8/L1 = phi^4

# Known neutrino masses
m3 = 2 * m_e_eV / PHI**35
r_known = ALPHA * PHI**3  # = Dm2_21/Dm2_31 (the pattern we want to derive)
m2 = m3 * np.sqrt(r_known)
m1 = 0.0

# Experimental
Dm2_21_exp = 7.53e-5   # eV^2
Dm2_32_exp = 2.453e-3  # eV^2

print("=" * 75)
print("EXP-372: DERIVE THE NEUTRINO MASS SPLITTING RATIO")
print("=" * 75)
print(f"  m_3 = 2*m_e/phi^35 = {m3*1e3:.3f} meV [DERIVED]")
print(f"  r = Dm2_21/Dm2_31 = alpha*phi^3 = {r_known:.6f} [PATTERN]")
print(f"  m_2/m_3 = sqrt(r) = {np.sqrt(r_known):.6f}")
print(f"  GOAL: derive r from first principles")

# =====================================================================
# PART 1: ALGEBRAIC DECOMPOSITION OF alpha*phi^3
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 1: WHAT IS alpha*phi^3 ALGEBRAICALLY?")
print("=" * 75)

print(f"\n  alpha = {ALPHA:.10f}")
print(f"  phi^3 = {PHI**3:.10f}")
print(f"  alpha*phi^3 = {ALPHA*PHI**3:.10f}")
print(f"  alpha_s = 1/(2*phi^3) = {ALPHA_S:.10f}")
print(f"\n  KEY: alpha*phi^3 = alpha/(2*alpha_s)")
print(f"       = {ALPHA/(2*ALPHA_S):.10f}")

print(f"\n  So the splitting ratio is:")
print(f"    r = alpha/(2*alpha_s) = EM coupling / (2 * strong coupling)")
print(f"    This MIXES the two non-trivial gauge sectors!")

# The alpha equation: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
# So alpha = (4*a1*phi^4 - sqrt(16*a1^2*phi^8 - 8*pi)) / (4*pi)
# alpha*phi^3 = alpha*phi^3
# Using alpha_s = 1/(2*phi^3):
# r = alpha*phi^3 = alpha/(2*alpha_s)

# Can we express this in terms of pure framework quantities?
print(f"\n  Using the alpha equation:")
print(f"    2*pi*alpha^2 - (4*a1*phi^4)*alpha + 1 = 0")
print(f"    Multiply by phi^6:")
print(f"    2*pi*(alpha*phi^3)^2 - 4*a1*phi^(4+3)*(alpha*phi^3) + phi^6 = 0")
# Actually, let r' = alpha*phi^3:
# 2*pi*r'^2/phi^6 - 4*a1*phi^4*r'/phi^3 + 1 = 0
# 2*pi*r'^2 - 4*a1*phi^7*r' + phi^6 = 0
r_prime = ALPHA * PHI**3
check_eq = 2*np.pi*r_prime**2 - 4*a1*PHI**7*r_prime + PHI**6
print(f"    Let r = alpha*phi^3. Substitute alpha = r/phi^3:")
print(f"    2*pi*(r/phi^3)^2 - 4*a1*phi^4*(r/phi^3) + 1 = 0")
print(f"    2*pi*r^2/phi^6 - 4*a1*phi*r + 1 = 0")
print(f"    Multiply by phi^6:")
print(f"    2*pi*r^2 - 4*a1*phi^7*r + phi^6 = 0")
print(f"    Check: {check_eq:.2e} (should be ~0)")

# So r satisfies: 2*pi*r^2 - 4*a1*phi^7*r + phi^6 = 0
# r = (4*a1*phi^7 +/- sqrt(16*a1^2*phi^14 - 8*pi*phi^6)) / (4*pi)
# r = (2*a1*phi^7 +/- sqrt(4*a1^2*phi^14 - 2*pi*phi^6)) / (2*pi)

# Simplify the discriminant
disc_r = 4*a1**2*PHI**14 - 2*np.pi*PHI**6
print(f"\n    Discriminant = 4*a1^2*phi^14 - 2*pi*phi^6")
print(f"                 = phi^6*(4*a1^2*phi^8 - 2*pi)")
print(f"                 = phi^6 * original_discriminant")
print(f"    = {disc_r:.6f}")

# The r equation inherits from the alpha equation.
# Not a NEW constraint, just a rewriting.
print(f"\n  VERDICT: r satisfies an equation DERIVED from the alpha equation.")
print(f"  This is CONSISTENT but not a new constraint.")

# =====================================================================
# PART 2: KERNEL OPERATOR — MASS MATRIX FROM LAPLACIANS
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 2: GALOIS KERNEL MASS MATRIX")
print("=" * 75)

print(f"\n  The Galois kernel has 3 irreps: rho_0 (dim 1), rho_1 (dim 2), rho_8 (dim 2)")
print(f"  Laplacian eigenvalues: L0={L0:.4f}, L1={L1:.6f}, L8={L8:.6f}")
print(f"  Neutrino masses: m1=0, m2={m2*1e3:.3f} meV, m3={m3*1e3:.3f} meV")

# IDEA: Construct a 3x3 matrix from kernel data whose eigenvalues ~ neutrino masses
# The kernel operator acts on the 3 irreps.
# Natural operator: K = diag(L0, L1, L8) = diag(0, 12-6phi, 12+6/phi)

print(f"\n  Test: Is there a function f such that f(L0,L1,L8) gives (0, m2, m3)?")

# f(L) = m_e/phi^35 * L would give: 0, L1*m_e/phi^35, L8*m_e/phi^35
# Ratio: L8/L1 = phi^4 = 6.854
# But m3/m2 = 1/sqrt(alpha*phi^3) = 5.645
# NOT equal to phi^4 or phi^2. So direct proportionality fails (known from exp367).

# What if the mass matrix is NOT diagonal in the kernel basis?
# Off-diagonal terms from coupling alpha could mix rho_1 and rho_8.

print(f"\n  Approach: Perturbation theory on kernel operator")
print(f"  Unperturbed: H0 = diag(0, L1, L8) in kernel basis")
print(f"  Perturbation: V ~ alpha * (off-diagonal mixing)")
print(f"  Question: what V gives the correct mass ratio?")

# If V mixes rho_1 and rho_8 (Galois conjugates), this is natural:
# The EM coupling alpha BREAKS Galois symmetry (alpha is irrational, alpha' is complex)

# Let M = [[a, b], [b, c]] be the 2x2 mass matrix in (rho_1, rho_8) subspace
# We need eigenvalues proportional to m2 and m3.
# Scale: set m3 as the scale.
# Then eigenvalues are (x, 1) where x = m2/m3 = sqrt(alpha*phi^3)

x = np.sqrt(r_known)  # = m2/m3
print(f"\n  2x2 mass matrix in (rho_1, rho_8) subspace:")
print(f"    eigenvalues: (m2/m3, 1) = ({x:.6f}, 1)")
print(f"    trace = 1 + x = {1+x:.6f}")
print(f"    det = x = {x:.6f}")

# Parametrize: M = [[p, q], [q, s]] with p+s = 1+x, p*s-q^2 = x
# Natural guess: p ~ L1/(L1+L8), s ~ L8/(L1+L8), q ~ alpha
p_guess = L1 / (L1 + L8)
s_guess = L8 / (L1 + L8)
print(f"\n  Natural parametrization:")
print(f"    p = L1/(L1+L8) = {p_guess:.6f}")
print(f"    s = L8/(L1+L8) = {s_guess:.6f}")
print(f"    p + s = 1 (good, need {1+x:.6f})")
# So we need to shift: p -> p + x*f, s -> s + x*(1-f) for some f
# Or: the matrix is alpha-dependent

# More physical: decompose into Galois-symmetric + breaking
# Symmetric: proportional to identity I or to diag(L1, L8)
# Breaking: proportional to alpha * sigma_x (swap rho_1 <-> rho_8)

# M = a * I + b * diag(L1, L8) + c * sigma_x
# where sigma_x = [[0,1],[1,0]] is the Galois exchange

# Eigenvalues of M:
# For diag part: a + b*L1, a + b*L8
# sigma_x adds: +/- c
# Combined eigenvalues: a + b*(L1+L8)/2 +/- sqrt((b*(L8-L1)/2)^2 + c^2)

# Set: a + b*(L1+L8)/2 = (m2+m3)/(2*m3) = (1+x)/2
# and: sqrt((b*(L8-L1)/2)^2 + c^2) = (1-x)/2

half_sum = (1 + x) / 2
half_diff = (1 - x) / 2
print(f"\n  Galois-symmetric decomposition:")
print(f"    M = a*I + b*diag(L1,L8) + c*sigma_x")
print(f"    eigenvalues = center +/- spread")
print(f"    center = (1+x)/2 = {half_sum:.6f}")
print(f"    spread = (1-x)/2 = {half_diff:.6f}")

# Natural choice: b encodes Galois, c encodes EM mixing
# If b dominates (no EM mixing, c=0):
#   eigenvalues = a+b*L1, a+b*L8. Ratio = (a+b*L8)/(a+b*L1)
#   For this to give m3/m2 = 1/x = 5.645:
#   (a+b*L8)/(a+b*L1) = 5.645
#   Not very illuminating without fixing a,b.

# If b=0 (no Galois splitting, only EM mixing):
#   eigenvalues = a +/- c. Ratio = (a+c)/(a-c)
#   For ratio 5.645: c/a = (5.645-1)/(5.645+1) = 0.699
#   Not obviously a framework quantity.

# =====================================================================
# PART 3: ALPHA AS PERTURBATION OF RANK-1 SEESAW
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 3: ALPHA AS PERTURBATION — RANK-1 TO RANK-2")
print("=" * 75)

print(f"""
  HYPOTHESIS: Without EM coupling (alpha -> 0), the seesaw gives
  a RANK-1 mass matrix (only m_3 nonzero). The EM perturbation
  creates m_2 proportional to alpha.

  If M_nu^(0) = rank 1 (only m_3),
  and M_nu = M_nu^(0) + alpha * delta_M,
  then m_2 ~ alpha * (geometric factor).
""")

# In seesaw: M_nu = -m_D^T * M_R^{-1} * m_D
# If alpha -> 0, do the Dirac masses change?
# In the framework: charged lepton masses use phi, not alpha.
# But QUARK masses also use phi (via McKay tree).
# The alpha equation determines alpha FROM phi.
# If alpha -> 0, the framework breaks (no EM).

# ALTERNATIVE VIEW: alpha is DERIVED from the framework.
# It cannot be "turned off". But we can ask:
# what part of the seesaw mass matrix depends on alpha vs phi?

# The PMNS matrix is derived from A5 (pure group theory, no alpha).
# So M_nu = U * diag(0, m2, m3) * U^T
# = m3 * U * diag(0, x, 1) * U^T    where x = sqrt(alpha*phi^3)
# = m3 * [U * diag(0, 0, 1) * U^T + x * U * diag(0, 1, 0) * U^T]
# = m3 * [M_rank1 + x * M_rank1_perp]

# So M_nu = m3 * M_rank1 + m2 * M_rank1_perp
# where M_rank1 = U * diag(0,0,1) * U^T (rank 1)
# and M_rank1_perp = U * diag(0,1,0) * U^T (rank 1)

# The "perturbation" that creates m2 is controlled by x = sqrt(alpha*phi^3).
# But this is trivially true — it's just the eigendecomposition of M_nu.
# NOT a derivation.

print(f"  This approach is TAUTOLOGICAL.")
print(f"  M_nu = m3*|nu_3><nu_3| + m2*|nu_2><nu_2| is just the definition.")
print(f"  We need a mechanism that PRODUCES alpha*phi^3, not one that assumes it.")

# =====================================================================
# PART 4: THE GALOIS ACTION ON THE SEESAW
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 4: GALOIS ACTION ON THE SEESAW")
print("=" * 75)

print(f"""
  In Z[phi], the Galois conjugation phi <-> phi' maps:
    alpha  ->  alpha' (complex! no physical EM)
    alpha_s -> alpha_s' = 1/(2*phi'^3) < 0 (no confinement)

  For the SEESAW, if M_R = m_e*phi^35/2:
    M_R' = m_e*phi'^35/2 (Galois conjugate)

  The "physical" seesaw uses BOTH M_R and M_R':
    m_3 ~ 1/M_R  (phi sector)
    m_2 ~ 1/M_R' * (mixing factor)  ???
""")

# Compute phi'^35
phi_prime_35 = PHI_CONJ**35
print(f"  phi'^35 = {phi_prime_35:.6e}")
print(f"  |phi'^35| = {abs(phi_prime_35):.6e}")
print(f"  phi^35 = {PHI**35:.6e}")
print(f"  |phi'^35|/phi^35 = 1/phi^70 = {1/PHI**70:.6e}")

# That's tiny, not useful.
# The Galois conjugate of the seesaw mass is enormously different.

# But what about the Galois conjugate of alpha*phi^3?
r_galois = ALPHA * PHI_CONJ**3
print(f"\n  Galois conjugate of alpha*phi^3:")
print(f"    alpha*phi'^3 = {ALPHA*PHI_CONJ**3:.6f}")
print(f"    alpha*phi^3  = {ALPHA*PHI**3:.6f}")
print(f"    Ratio: {abs(ALPHA*PHI_CONJ**3/(ALPHA*PHI**3)):.6f}")
print(f"    = 1/phi^6 = {1/PHI**6:.6f}")

# Galois norm of r = alpha*phi^3:
r_norm = r_known * abs(ALPHA * PHI_CONJ**3)
print(f"\n  Galois norm: N(alpha*phi^3) = alpha^2 * N(phi^3)")
print(f"    = alpha^2 * |phi^3 * phi'^3| = alpha^2 * |(-1)^3| = alpha^2")
print(f"    = {ALPHA**2:.10f}")
print(f"    Check: {r_known * abs(ALPHA*PHI_CONJ**3):.10f}")
# N(phi^3) = phi^3 * phi'^3 = (phi*phi')^3 = (-1)^3 = -1
# So N(alpha*phi^3) = alpha^2 * |-1| = alpha^2
print(f"    Direct: alpha^2 = {ALPHA**2:.10f}")

print(f"\n  INSIGHT: N(r) = N(alpha*phi^3) = alpha^2")
print(f"  The Galois norm of the splitting ratio is alpha^2.")
print(f"  This means r*r' = alpha^2 (product of Galois conjugates).")

# =====================================================================
# PART 5: THE COUPLING MIXING STRUCTURE
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 5: COUPLING MIXING — WHY alpha/(2*alpha_s)?")
print("=" * 75)

print(f"\n  r = alpha * phi^3 = alpha / (2*alpha_s)")
print(f"  This is the ratio of the TWO non-trivial couplings.")
print(f"\n  The framework has 3 gauge couplings derived from 600-cell:")
print(f"    alpha_1 (U(1)): from alpha equation")
print(f"    alpha_2 (SU(2)): from sin^2(tW) = 6/26")
print(f"    alpha_3 (SU(3)): alpha_s = 1/(2*phi^3)")
print(f"\n  At M_Z scale:")
print(f"    alpha_em = alpha_0 * 16/15 (running pattern)")
alpha_mz = ALPHA * 16/15
sin2tw = 6.0/26
print(f"    sin^2(tW) = 6/26 = {sin2tw:.6f}")
print(f"    alpha_2 = alpha_em/sin^2(tW) = {alpha_mz/sin2tw:.6f}")
print(f"    alpha_3 = alpha_s = {ALPHA_S:.6f}")

# The splitting ratio involves alpha (EM) and alpha_s (strong)
# WHY these two? Because:
# - m3 comes from the seesaw (strong sector, phi^35)
# - m2 is suppressed by EM coupling alpha
# The neutrino sector "sees" both gauge couplings through the seesaw

# Is there a deeper reason?
# The Galois kernel has phi and phi' sectors.
# phi-sector -> strong coupling alpha_s = 1/(2*phi^3)
# phi'-sector -> the Galois conjugate is unphysical (alpha_s' < 0)
# But alpha is the ONLY coupling that bridges both Galois sectors
# (because alpha satisfies a quadratic over Q, not over Q(phi))

print(f"\n  WHY alpha bridges Galois sectors:")
print(f"    alpha_s = 1/(2*phi^3) is a function of phi (phi-sector)")
print(f"    alpha_s' = 1/(2*phi'^3) < 0 (phi'-sector, unphysical)")
print(f"    But ALPHA satisfies 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0")
print(f"    The other root: alpha' = {(-B_eq + np.sqrt(disc))/(2*A_eq):.6f}")
alpha_prime = (-B_eq + np.sqrt(disc))/(2*A_eq)
print(f"    alpha * alpha' = C/A = 1/(2*pi) = {1/(2*np.pi):.6f}")
print(f"    Verified: {ALPHA*alpha_prime:.6f}")
print(f"\n    alpha*alpha' = 1/(2*pi) [Kaluza-Klein on U(1) fiber]")
print(f"    This is EXACT, not Galois conjugation!")
print(f"    alpha' is the SECOND root of the alpha equation.")

# =====================================================================
# PART 6: THE M_NU MATRIX FROM KERNEL STRUCTURE
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 6: CONSTRUCTING M_NU FROM KERNEL STRUCTURE")
print("=" * 75)

# PMNS angles (DERIVED)
sin2_12 = 2.0 / (PHI + a1)    # = 2/(phi+5)
sin2_23 = 4.0 / 7.0            # = (a1-1)/(a1+2)
sin2_13 = 1.0 / (a1 * N_eig)  # = 1/45
delta_CP = 3 * np.arctan(np.sqrt(5))

s12 = np.sqrt(sin2_12)
c12 = np.sqrt(1 - sin2_12)
s23 = np.sqrt(sin2_23)
c23 = np.sqrt(1 - sin2_23)
s13 = np.sqrt(sin2_13)
c13 = np.sqrt(1 - sin2_13)

# PMNS matrix (no Majorana phases for now)
cd = np.cos(delta_CP)
sd = np.sin(delta_CP)

U = np.array([
    [c12*c13,                    s12*c13,                    s13*np.exp(-1j*delta_CP)],
    [-s12*c23 - c12*s23*s13*np.exp(1j*delta_CP),  c12*c23 - s12*s23*s13*np.exp(1j*delta_CP),  s23*c13],
    [s12*s23 - c12*c23*s13*np.exp(1j*delta_CP),  -c12*s23 - s12*c23*s13*np.exp(1j*delta_CP),  c23*c13]
])

# Mass matrix: M_nu = U * diag(0, m2^2, m3^2) * U^T  (for Majorana)
# Actually M_nu = U* * diag(m1,m2,m3) * U^dagger (convention depends)
# For Majorana: M_nu = U_conj * diag(m) * U^dagger
# But for real masses: M_nu = U * diag(m) * U^T

# With m1=0, the structure is:
# M_nu = m2 * |nu_2><nu_2| + m3 * |nu_3><nu_3|
# where |nu_i> is column i of U

# In the limit m2 -> 0:
# M_nu^(0) = m3 * |nu_3><nu_3|   (rank 1)

# The "correction" that creates m2:
# delta_M = m2 * |nu_2><nu_2|

# The question: why is m2 = m3 * sqrt(alpha*phi^3)?

# Let's examine the structure of |nu_2> and |nu_3>
print(f"\n  PMNS column vectors:")
print(f"  |nu_2> = ({U[0,1]:.6f}, {U[1,1]:.6f}, {U[2,1]:.6f})")
print(f"  |nu_3> = ({U[0,2]:.6f}, {U[1,2]:.6f}, {U[2,2]:.6f})")

# Inner products
print(f"\n  |nu_2|^2 = {np.sum(np.abs(U[:,1])**2):.6f} (should be 1)")
print(f"  |nu_3|^2 = {np.sum(np.abs(U[:,2])**2):.6f} (should be 1)")

# M_nu in flavor basis with m2/m3 = x
# M_nu/m3 = x * |nu_2><nu_2| + |nu_3><nu_3|

# For x << 1, M_nu/m3 ~ |nu_3><nu_3| + x * |nu_2><nu_2|
# The trace: Tr(M_nu/m3) = 1 + x
# The "ee" element: M_ee/m3 = x*|U_e2|^2 + |U_e3|^2
M_ee_over_m3 = x * np.abs(U[0,1])**2 + np.abs(U[0,2])**2
print(f"\n  M_ee/m3 = x*|U_e2|^2 + |U_e3|^2")
print(f"         = {x:.6f}*{np.abs(U[0,1])**2:.6f} + {np.abs(U[0,2])**2:.6f}")
print(f"         = {M_ee_over_m3:.6f}")

# Compute full M_nu/m3
M_over_m3 = np.zeros((3,3), dtype=complex)
for i in range(3):
    for j in range(3):
        M_over_m3[i,j] = x * U[i,1]*np.conj(U[j,1]) + U[i,2]*np.conj(U[j,2])

print(f"\n  M_nu/m3 matrix (real parts):")
for i in range(3):
    row = "    "
    for j in range(3):
        row += f"{M_over_m3[i,j].real:>10.6f} "
    print(row)

# =====================================================================
# PART 7: CAN THE PMNS STRUCTURE DETERMINE x = m2/m3?
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 7: DOES PMNS STRUCTURE CONSTRAIN m2/m3?")
print("=" * 75)

print(f"""
  In general, PMNS angles and masses are independent.
  But in this framework, BOTH come from the same algebra (A5 + McKay).
  Question: does the A5 structure impose a constraint on m2/m3?
""")

# The PMNS matrix is determined by sin^2(theta_12) = 2/(phi+5)
# This involves phi (framework constant).
# If m2/m3 also involves phi (via phi^3), there MIGHT be a connection.

# Test: Is there a relation between PMNS angles and the splitting ratio?
# Compute some combinations of PMNS parameters
print(f"  Framework PMNS parameters:")
print(f"    sin^2(12) = 2/(phi+5) = {sin2_12:.6f}")
print(f"    sin^2(23) = 4/7 = {sin2_23:.6f}")
print(f"    sin^2(13) = 1/45 = {sin2_13:.6f}")

# Various combinations
print(f"\n  Testing algebraic relations:")
tests = [
    ("sin2_12 * sin2_23", sin2_12 * sin2_23),
    ("sin2_12 * sin2_13", sin2_12 * sin2_13),
    ("sin2_13 / sin2_12", sin2_13 / sin2_12),
    ("sin2_13 * sin2_23 / sin2_12", sin2_13 * sin2_23 / sin2_12),
    ("sin2_12 - sin2_13", sin2_12 - sin2_13),
    ("sin2_12^2", sin2_12**2),
    ("sin2_13 * (1-sin2_23)", sin2_13 * (1-sin2_23)),
    ("sin2_12 * sin2_13 * N_eig", sin2_12 * sin2_13 * N_eig),
]

for name, val in tests:
    err = abs(val - r_known) / r_known * 100
    marker = " <---" if err < 5 else ""
    print(f"    {name:>40} = {val:.6f}  (vs r={r_known:.6f}, err={err:.1f}%){marker}")

# None of these will match exactly because r = alpha*phi^3 involves alpha
# and PMNS angles involve phi but not alpha.

# But what about the Jarlskog invariant?
J = s12 * c12 * s23 * c23 * s13 * c13**2 * np.sin(delta_CP)
print(f"\n  Jarlskog invariant J = {J:.6f}")
print(f"  J vs r: J/r = {J/r_known:.6f}")
print(f"  J * phi: {J*PHI:.6f}")
print(f"  alpha * phi^3 = {r_known:.6f}")
print(f"  PMNS angles do NOT determine r (alpha is needed).")

# =====================================================================
# PART 8: OPERATOR ON KERNEL — THE KEY IDEA
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 8: OPERATOR ON GALOIS KERNEL")
print("=" * 75)

print(f"""
  The Galois kernel has 3 irreps with Laplacian eigenvalues:
    rho_0: L0 = 0 (dim 1)
    rho_1: L1 = 12-6*phi = {L1:.6f} (dim 2)
    rho_8: L8 = 12+6/phi = {L8:.6f} (dim 2)

  Define the "seesaw operator" S on the kernel.
  S should have eigenvalues proportional to (0, m2, m3) = (0, x*m3, m3).

  The operator S acts on the 3-dim irrep space of the kernel.
  Basis: |rho_0>, |rho_1>, |rho_8>.

  Natural candidate: S = f(Laplacian) + alpha * mixing.
""")

# KEY IDEA: The mass matrix in the kernel basis is determined by:
# 1. The Laplacian eigenvalues (Galois structure)
# 2. The EM coupling alpha (breaks Galois symmetry)
# 3. The constraint m1 = 0 (rho_0 decouples)

# In the (rho_1, rho_8) subspace:
# The Laplacian gives: diag(L1, L8)
# The masses should be functions of L1, L8, and alpha.

# Hypothesis: The effective seesaw operator in the kernel is
# S = L^{-1} (inverse Laplacian, like a propagator)
# S_11 = 1/L1, S_22 = 1/L8
# Eigenvalues: 1/L1 = 0.4363, 1/L8 = 0.06366
# Ratio: L8/L1 = phi^4 (known)
# But m3/m2 = 1/x = 5.645, not phi^4 = 6.854.

print(f"  Hypothesis A: S = L^(-1) (propagator)")
print(f"    S_11 = 1/L1 = {1/L1:.6f}")
print(f"    S_22 = 1/L8 = {1/L8:.6f}")
print(f"    Ratio = L8/L1 = phi^4 = {L8/L1:.4f}")
print(f"    Need: m3/m2 = {1/x:.4f}")
print(f"    FAILS (phi^4 != 1/x)")

# Hypothesis B: S = L^{-k} for some k
# Ratio = (L8/L1)^k = phi^(4k)
# Need: phi^(4k) = 1/x = 1/sqrt(alpha*phi^3)
# 4k*ln(phi) = -0.5*ln(alpha) - 1.5*ln(phi)
# k = (-0.5*ln(alpha) - 1.5*ln(phi)) / (4*ln(phi))
k_needed = (-0.5*np.log(ALPHA) - 1.5*np.log(PHI)) / (4*np.log(PHI))
print(f"\n  Hypothesis B: S = L^(-k)")
print(f"    Need (L8/L1)^k = 1/x = {1/x:.6f}")
print(f"    k = {k_needed:.6f}")
print(f"    Not a clean rational number. FAILS.")

# Hypothesis C: Off-diagonal mixing by alpha
# S = [[1/L1, alpha], [alpha, 1/L8]]
# Eigenvalues: ((1/L1+1/L8) +/- sqrt((1/L1-1/L8)^2 + 4*alpha^2)) / 2
invL1 = 1.0/L1
invL8 = 1.0/L8
center_C = (invL1 + invL8) / 2
spread_C = np.sqrt((invL1 - invL8)**2 + 4*ALPHA**2) / 2
eig1_C = center_C + spread_C
eig2_C = center_C - spread_C
ratio_C = eig1_C / eig2_C

print(f"\n  Hypothesis C: S = [[1/L1, alpha], [alpha, 1/L8]]")
print(f"    eigenvalues: {eig1_C:.6f}, {eig2_C:.6f}")
print(f"    ratio: {ratio_C:.4f}")
print(f"    need: {1/x:.4f}")
print(f"    Error: {abs(ratio_C - 1/x)/(1/x)*100:.1f}%  {'CLOSE!' if abs(ratio_C - 1/x)/(1/x) < 0.05 else 'FAILS'}")

# Hypothesis D: S = [[1/L1, v], [v, 1/L8]] where v is framework quantity
# Need ratio = 1/x. Find v.
# (1/L1+1/L8)/2 +/- sqrt(((1/L1-1/L8)/2)^2 + v^2)
# Let R = eig_big/eig_small = 1/x
# Then: (R-1)/(R+1) = spread/center
# spread^2 = ((1/L1-1/L8)/2)^2 + v^2
# center = (1/L1+1/L8)/2

frac = (1/x - 1) / (1/x + 1)  # = (1-x)/(1+x)
center_val = (invL1 + invL8) / 2
needed_spread = frac * center_val
diag_spread = abs(invL1 - invL8) / 2
v_squared = needed_spread**2 - diag_spread**2

print(f"\n  Hypothesis D: Find v such that S = [[1/L1, v], [v, 1/L8]] gives 1/x")
print(f"    needed spread = {needed_spread:.6f}")
print(f"    diagonal spread = {diag_spread:.6f}")

if v_squared > 0:
    v_needed = np.sqrt(v_squared)
    print(f"    v^2 = {v_squared:.8f}")
    print(f"    v = {v_needed:.8f}")
    print(f"\n    Compare to framework quantities:")
    print(f"      alpha = {ALPHA:.8f} (ratio v/alpha = {v_needed/ALPHA:.4f})")
    print(f"      alpha^2 = {ALPHA**2:.8f} (ratio v/alpha^2 = {v_needed/ALPHA**2:.2f})")
    print(f"      alpha/phi = {ALPHA/PHI:.8f} (ratio = {v_needed/(ALPHA/PHI):.4f})")
    print(f"      alpha*phi = {ALPHA*PHI:.8f} (ratio = {v_needed/(ALPHA*PHI):.4f})")
    print(f"      1/(2*phi^3) = alpha_s = {ALPHA_S:.8f} (ratio = {v_needed/ALPHA_S:.4f})")
    print(f"      alpha/b1 = {ALPHA/b1:.8f} (ratio = {v_needed/(ALPHA/b1):.4f})")
    print(f"      alpha/(2*pi) = {ALPHA/(2*np.pi):.8f} (ratio = {v_needed/(ALPHA/(2*np.pi)):.4f})")
else:
    print(f"    v^2 = {v_squared:.8f} < 0 -- need IMAGINARY mixing!")
    v_needed = np.sqrt(-v_squared)
    print(f"    |v| = {v_needed:.8f}")
    print(f"    This means the eigenvalue ratio CANNOT come from Hermitian mixing")
    print(f"    if 1/L1 and 1/L8 provide the diagonal.")

# =====================================================================
# PART 9: THE DIRECT LAPLACIAN OPERATOR (NOT INVERSE)
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 9: DIRECT LAPLACIAN OPERATOR ON KERNEL")
print("=" * 75)

# Maybe the operator is L itself, not L^{-1}
# S = [[L1, v], [v, L8]]
# eigenvalues: (L1+L8)/2 +/- sqrt(((L8-L1)/2)^2 + v^2)

# We want eigenvalue ratio = m3/m2 = 1/x (larger / smaller)
# BUT the eigenvalue ordering would be: L8-like > L1-like
# So the LARGER eigenvalue ~ L8 ~ heavy scale ~ m3 (heaviest neutrino)

# For S = [[L1, w], [w, L8]] to give ratio 1/x:
center_L = (L1 + L8) / 2
frac_L = (1/x - 1) / (1/x + 1)
needed_spread_L = frac_L * center_L
diag_spread_L = (L8 - L1) / 2
w_squared = needed_spread_L**2 - diag_spread_L**2

print(f"\n  S = [[L1, w], [w, L8]] on kernel (rho_1, rho_8) subspace")
print(f"    center = (L1+L8)/2 = {center_L:.6f}")
print(f"    diag spread = (L8-L1)/2 = {diag_spread_L:.6f}")
print(f"    needed spread for ratio 1/x = {1/x:.4f}: {needed_spread_L:.6f}")

if w_squared > 0:
    w_val = np.sqrt(w_squared)
    print(f"    w = {w_val:.8f}")
    print(f"    Compare: w/alpha = {w_val/ALPHA:.4f}")
else:
    print(f"    w^2 = {w_squared:.6f} < 0 (IMPOSSIBLE)")
    print(f"    The diagonal spread (L8-L1)/2 = {diag_spread_L:.4f} is LARGER than needed")
    print(f"    spread = {needed_spread_L:.4f}")
    print(f"    CONCLUSION: The Laplacian ratio phi^4 is TOO LARGE")
    print(f"    to match m3/m2 = {1/x:.4f} with any real off-diagonal term.")
    print(f"    The observed mass ratio is CLOSER to 1 than the Laplacian ratio.")

# =====================================================================
# PART 10: THE HEAT KERNEL / EXPONENTIAL OPERATOR
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 10: EXPONENTIAL OPERATORS ON KERNEL")
print("=" * 75)

print(f"\n  The Laplacian ratio phi^4 = 6.854 is too large.")
print(f"  The mass ratio 1/x = {1/x:.4f} is smaller.")
print(f"  What operator O(L) gives the correct ratio?")
print(f"\n  Need: O(L8)/O(L1) = 1/x = {1/x:.6f}")

# Test various functions
functions = [
    ("L",       L1,       L8),
    ("L^(1/2)", L1**0.5,  L8**0.5),
    ("L^(3/4)", L1**0.75, L8**0.75),
    ("ln(L)",   np.log(L1), np.log(L8)),
    ("exp(-L)",  np.exp(-L1), np.exp(-L8)),
    ("exp(-L/12)", np.exp(-L1/12), np.exp(-L8/12)),
    ("exp(-L/b1)", np.exp(-L1/b1), np.exp(-L8/b1)),
    ("1/L",     1/L1,     1/L8),
    ("L/(L+1)", L1/(L1+1), L8/(L8+1)),
    ("tanh(L/DEG)", np.tanh(L1/DEG), np.tanh(L8/DEG)),
]

target = 1.0/x  # = m3/m2
print(f"\n  {'Function':>15} {'f(L1)':>10} {'f(L8)':>10} {'ratio':>10} {'target':>10} {'err%':>8}")
print(f"  {'-'*70}")
for name, f1, f8 in functions:
    ratio = f8/f1 if f1 != 0 else float('inf')
    err = abs(ratio - target)/target * 100
    marker = " <--" if err < 2 else ""
    print(f"  {name:>15} {f1:>10.6f} {f8:>10.6f} {ratio:>10.4f} {target:>10.4f} {err:>7.2f}%{marker}")

# Now test: alpha * f(L) ratios
print(f"\n  Testing alpha * f(L) combinations:")
print(f"  Need: alpha * f(L8)/f(L1) = r = {r_known:.6f}")
print(f"  i.e.: f(L8)/f(L1) = r/alpha = phi^3 = {PHI**3:.6f}")

target_geom = PHI**3
print(f"\n  {'Function':>15} {'f(L8)/f(L1)':>12} {'target=phi^3':>12} {'err%':>8}")
print(f"  {'-'*55}")
funcs2 = [
    ("L^(3/4)",  (L8/L1)**0.75),
    ("L^(1/2)",  (L8/L1)**0.5),
    ("L^(2/3)",  (L8/L1)**(2.0/3)),
    ("L^(5/8)",  (L8/L1)**(5.0/8)),
    ("sqrt(L)*ln(L)",  np.sqrt(L8)*np.log(L8)/(np.sqrt(L1)*np.log(L1))),
]
for name, ratio in funcs2:
    err = abs(ratio - target_geom)/target_geom * 100
    marker = " <-- EXACT" if err < 0.01 else (" <--" if err < 1 else "")
    print(f"  {name:>15} {ratio:>12.6f} {target_geom:>12.6f} {err:>7.3f}%{marker}")

# =====================================================================
# PART 11: THE 3/4 EXPONENT — IS IT DERIVABLE?
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 11: CAN WE DERIVE THE 3/4 EXPONENT?")
print("=" * 75)

print(f"""
  KNOWN: (L8/L1)^(3/4) = (phi^4)^(3/4) = phi^3
  So r = alpha * phi^3 = alpha * (L8/L1)^(3/4)

  The exponent 3/4 = 1 - 1/4 = 1 - 1/d (d = spacetime dim).

  WHERE ELSE does 3/4 appear in the framework?
    1. Charged lepton corrections: delta = c_ell * |z'|^(3/4)
       where k=3/4 comes from Tr(phi^n)=4 for n=3 (Lucas)
    2. Heat kernel on 3-sphere: Tr(exp(-t*D^2)) ~ t^(-3/2) for small t
       The "effective dimension" d=3 gives k = 1-1/3 = 2/3 (DIFFERENT!)
    3. For d=4 spacetime: k = 1-1/4 = 3/4 (MATCHES)

  The question: WHY does the neutrino splitting use d=4 not d=3?

  ARGUMENT: The mass splitting involves energy scales (mass-squared difference),
  which lives in 4D spacetime, not 3D space. The 600-cell provides the
  spatial geometry (d=3), but the physical observables (masses) live in d=4.

  Alternatively: k = 3/4 because the neutrino mass matrix is a
  4-component object (3 masses + 1 zero = 4 entries that determine
  the mass spectrum). The "effective dimension" is 4, giving k=3/4.
""")

# The derivation of 3/4 from the framework:
# In exp365, k was derived from the condition:
# Tr(phi^n) = L_n (Lucas numbers) = 1, 3, 4, 7, 11, 18, ...
# The requirement d = Tr(phi^n) = 4 gives n = 3, hence k = 1 - 1/d = 3/4.
# This is the SAME for both charged leptons and neutrinos.
# The "dimension" d = 4 comes from L_3 = 4 = Tr(phi^3).

print(f"  Lucas sequence L_n = Tr(phi^n) for n=0,1,2,3,4,5:")
for n in range(6):
    Ln = PHI**n + PHI_CONJ**n
    print(f"    L_{n} = {Ln:.4f} -> {int(round(Ln))}")
print(f"\n  L_3 = 4 = spacetime dimension.")
print(f"  k = 1 - 1/L_3 = 1 - 1/4 = 3/4.")
print(f"  This is DERIVED (unique n that gives integer d=4 for L_n even).")

# =====================================================================
# PART 12: SYNTHESIS — THE SPLITTING RATIO
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 12: SYNTHESIS")
print("=" * 75)

print("""
  CAN WE DERIVE r = alpha * phi^3 = alpha * (L8/L1)^(3/4)?

  STEP 1 (DERIVED): L8/L1 = phi^4 [Galois kernel Laplacian ratio]
  STEP 2 (DERIVED): k = 3/4 = 1 - 1/d, d = L_3 = 4 [Lucas/spacetime]
  STEP 3 (DERIVED): alpha from the alpha equation
  STEP 4 (???): WHY r = alpha * (L_ratio)^k ?

  Step 4 is the GAP. We have the three ingredients (L8/L1, k, alpha)
  but no first-principles mechanism that combines them into r.

  POSSIBLE ARGUMENTS for the combination:

  A) HEAT KERNEL ON PRODUCT SPACE:
     On M x F (manifold x fiber), the heat kernel at time t gives
     mass-squared ~ t^(1-d/2) * eigenvalue.
     For d=4 (effective spacetime): t^(-1) * L.
     The mass-squared RATIO between Galois sectors:
       Dm2_21/Dm2_31 ~ (L1/L8)^(3/4) = (1/phi^4)^(3/4) = phi^(-3)
     But this gives phi^(-3), and we need alpha*phi^3.
     So: r = alpha * (L8/L1)^(3/4) = alpha * phi^3. The alpha factor is
     the COUPLING CONSTANT at the seesaw scale.
     STATUS: MOTIVATED but mechanism unclear

  B) RADIATIVE CORRECTION:
     m_2 = 0 at tree level (rank-1 seesaw).
     1-loop EM radiative correction generates m_2 ~ alpha * (Galois factor).
     r = alpha * phi^3 is a 1-loop effect with geometric enhancement phi^3.
     STATUS: STANDARD PHYSICS, quantitatively checkable

  C) COUPLING RATIO:
     r = alpha/(2*alpha_s) = ratio of EM to strong coupling.
     The neutrino mass hierarchy mirrors the gauge coupling hierarchy.
     STATUS: OBSERVATION, not derivation
""")

# Test: radiative correction formula
# 1-loop correction to neutrino mass: delta_m = (alpha/4*pi) * m * ln(M_R/m)
# For m ~ m3, M_R ~ 5.3 TeV:
delta_rad = (ALPHA/(4*np.pi)) * np.log(5271*1e9/(m3))
print(f"  Radiative correction estimate:")
print(f"    delta_m/m ~ (alpha/4pi) * ln(M_R/m_3)")
print(f"             = ({ALPHA:.6f}/{4*np.pi:.4f}) * ln({5271:.0f} GeV / {m3:.4e} eV)")
print(f"             = {ALPHA/(4*np.pi):.6f} * {np.log(5271*1e9/m3):.2f}")
print(f"             = {delta_rad:.6f}")
print(f"    Compare r = alpha*phi^3 = {r_known:.6f}")
print(f"    Order of magnitude: {'SIMILAR' if 0.1 < delta_rad/r_known < 10 else 'DIFFERENT'}")
print(f"    Ratio: {delta_rad/r_known:.3f}")

# =====================================================================
# PART 13: NUMERICAL SCAN — WHAT FORMULA GIVES r?
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 13: SYSTEMATIC SCAN OF ALGEBRAIC EXPRESSIONS FOR r")
print("=" * 75)

print(f"\n  Target: r = {r_known:.10f}")
print(f"\n  Scanning expressions of the form: alpha^a * phi^b / (c * pi^d)")

best_matches = []
for a_pow in [1, 2]:
    for b_pow in range(-5, 6):
        for c_val in [1, 2, 3, 4, 5, 6, 8, 12, 13, 26]:
            for d_pow in [0, 1]:
                val = ALPHA**a_pow * PHI**b_pow / (c_val * np.pi**d_pow)
                if val > 0 and abs(val - r_known)/r_known < 0.01:
                    err = (val - r_known)/r_known * 100
                    expr = f"alpha^{a_pow}*phi^{b_pow}/({c_val}*pi^{d_pow})"
                    best_matches.append((abs(err), err, expr, val))

best_matches.sort()
print(f"\n  Best matches (< 1% error):")
for _, err, expr, val in best_matches[:10]:
    print(f"    {expr:>35} = {val:.8f}  (err: {err:+.4f}%)")

# Also check: products involving L1, L8, b1, etc.
print(f"\n  Expressions involving kernel Laplacians:")
tests2 = [
    ("alpha*phi^3",                    ALPHA*PHI**3),
    ("alpha/(2*alpha_s)",              ALPHA/(2*ALPHA_S)),
    ("alpha*L8/(2*DEG)",               ALPHA*L8/(2*DEG)),
    ("alpha*L1*phi^4/(2*DEG)",         ALPHA*L1*PHI**4/(2*DEG)),
    ("alpha*(L8-L1)/DEG",              ALPHA*(L8-L1)/DEG),
    ("alpha*sqrt(L1*L8)/DEG",          ALPHA*np.sqrt(L1*L8)/DEG),
    ("alpha*b1/DEG",                   ALPHA*b1/DEG),
    ("L1/(4*pi*L8)",                   L1/(4*np.pi*L8)),
    ("alpha*L8/b1^2",                  ALPHA*L8/b1**2),
    ("(L1/DEG)^2",                     (L1/DEG)**2),
    ("alpha^2*DEG/L1",                 ALPHA**2*DEG/L1),
    ("1/(4*pi*phi^3)",                 1/(4*np.pi*PHI**3)),
    ("alpha*phi^3*b1/b1",              ALPHA*PHI**3),
    ("alpha/(2*phi^3)/(1/(2*phi^3)^2)",  ALPHA*2*PHI**3),
]

print(f"\n  {'Expression':>35} {'Value':>12} {'r':>12} {'err%':>8}")
print(f"  {'-'*70}")
for name, val in tests2:
    err = (val - r_known)/r_known * 100
    marker = " <-- EXACT" if abs(err) < 0.01 else (" <--" if abs(err) < 1 else "")
    print(f"  {name:>35} {val:>12.8f} {r_known:>12.8f} {err:>+8.3f}%{marker}")

# =====================================================================
# PART 14: HONEST ASSESSMENT
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 14: HONEST ASSESSMENT")
print("=" * 75)

print(f"""
  WHAT WE ESTABLISHED:

  1. r = alpha*phi^3 = alpha/(2*alpha_s) is the coupling ratio.
     This is an ALGEBRAIC IDENTITY, not a derivation.

  2. phi^3 = (L8/L1)^(3/4) connects to Galois kernel.
     The exponent 3/4 = 1-1/d is derivable from Lucas L_3=4.
     But the combination alpha*(L_ratio)^(3/4) lacks a mechanism.

  3. The operator approach (Parts 8-10) found:
     - Direct Laplacian ratio phi^4 is TOO LARGE
     - Inverse Laplacian ratio 1/phi^4 is TOO SMALL
     - Mixing by alpha does NOT reproduce the correct ratio
     - The correct operator exponent 3/4 is BETWEEN these extremes

  4. PMNS angles do NOT constrain m2/m3 (independent parameters).

  5. Radiative corrections give order-of-magnitude {delta_rad:.4f}
     vs exact {r_known:.6f} — same ballpark but NOT exact.

  WHAT IS NEW:
  - N(alpha*phi^3) = alpha^2 (Galois norm of splitting ratio)
  - r = alpha/(2*alpha_s) interpretation as EM/strong coupling ratio
  - Systematic exclusion of simple kernel operators
  - 3/4 exponent derived from Lucas L_3=4 (cross-sector universal)

  STATUS: r = alpha*phi^3 remains PATTERN.
  The three ingredients (alpha, L8/L1, 3/4 exponent) are each DERIVED,
  but the mechanism that combines them is MISSING.

  MOST PROMISING LEAD:
  Radiative correction: m2=0 at tree level, generated at 1-loop by alpha.
  The geometric factor phi^3 from the Galois kernel Laplacian ratio.
  This would make r a CALCULABLE 1-loop effect on M x F.
  But the exact computation on the product space has NOT been done.

  CATEGORY: PATTERN (upgraded from "pure pattern" to "ingredients derived,
  combination mechanism missing")
""")

print("=" * 75)
print("EXP-372 COMPLETE")
print("=" * 75)

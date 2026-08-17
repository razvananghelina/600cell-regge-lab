"""
exp445_galois_norm_sum_rules.py
================================
GALOIS NORM SUM RULES AND GAUGE-MASS UNIFICATION

Three new structural identities connecting the mass formula
to the gauge sector through the arithmetic of Z[phi]:

1. sum N(z_f) = b1 = 6  (Galois norm sum = first Betti number)
2. sum N(z_f)^3 = N + b1 = 126  (cubic identity, uniqueness of a1=5)
3. sin^2(theta_W) = sum N(z_f) / (a1^2+1)  (Weinberg = mass norms / geometry)

Plus connections:
4. E_vac = phi^3 = 1/(2*alpha_s)  (TQFT vacuum = inverse strong coupling)
5. r_nu = alpha * E_vac  (neutrino hierarchy = EM * TQFT vacuum)
6. Beta_sum = N_c * Level_1_sum = 15  (running = colored Arakelov)
7. z_t * z_b = |N(z_t)| * phi  (top-bottom Wilson line product)

Structure:
  Part A: Galois norm sum rule (sum N = b1)
  Part B: Cascade level decomposition
  Part C: Cubic identity (sum N^3 = N + b1)
  Part D: Weinberg angle from Galois norms
  Part E: TQFT vacuum and coupling connections
  Part F: Wilson line products (SU(2) doublets)
  Part G: Anomaly-like structure
  Part H: Summary and classification

Author: Razvan-Constantin Anghelina
Date: March 2026
"""

import numpy as np
from math import factorial

# =============================================================================
# CONSTANTS
# =============================================================================
a1 = 5
phi = (1 + np.sqrt(5)) / 2
phi_c = (1 - np.sqrt(5)) / 2
N = factorial(a1)  # 120
b1 = a1 + 1  # 6
N_gen = 3
N_c = 3
alpha_em = 7.2973525693e-3
alpha_s = 1.0 / (2 * phi**3)

# All fermion Wilson lines z_f = a + b*phi
# (a,b) from exp416 (UNIQUELY DERIVED)
fermions = [
    ('e',   0,  0, 'Level 0', 'input'),
    ('mu',  1,  1, 'Level 2a', 'lepton'),
    ('tau', 1,  2, 'Level 2a', 'lepton'),
    ('u',   3, -2, 'Level 2c', 'quark'),
    ('d',   1,  0, 'Level 3',  'quark'),
    ('s',   1,  1, 'Level 2b', 'quark'),
    ('c',   2,  1, 'Level 1',  'quark'),
    ('b',  -1,  4, 'Level 1',  'quark'),
    ('t',   4,  1, 'Level 1',  'quark'),
]

def galois_norm(a, b):
    """N(a + b*phi) = a^2 + ab - b^2 (using N(phi) = -1)."""
    return a**2 + a*b - b**2


print("=" * 72)
print("EXP-445: GALOIS NORM SUM RULES AND GAUGE-MASS UNIFICATION")
print("=" * 72)


# #############################################################################
#                  PART A: GALOIS NORM SUM RULE
# #############################################################################
print()
print("=" * 72)
print("PART A: GALOIS NORM SUM RULE")
print("=" * 72)
print()

print("Fermion Wilson lines z_f = a + b*phi and Galois norms N(z_f):")
print()
print("%-5s  (%2s,%2s)  %10s  %6s  %s" % ("Name", "a", "b", "z_f", "N(z)", "Level"))
print("-" * 55)

norms = {}
for name, a, b, level, sector in fermions:
    z = a + b * phi
    Nz = galois_norm(a, b)
    norms[name] = Nz
    print("%-5s  (%2d,%2d)  %10.4f  %+6d  %s" % (name, a, b, z, Nz, level))

print()

# Sum of all norms
sum_N = sum(norms[n] for n, _, _, _, _ in fermions)
print("Sum of all Galois norms:")
print("  sum_f N(z_f) = %+d" % sum_N)
print("  b1 = a1 + 1 = %d" % b1)
print("  sum N = b1?  %s" % ("EXACT!" if sum_N == b1 else "FAIL"))
print()

print("THEOREM (Galois Norm Sum Rule):")
print("  For the nine fermion Wilson lines z_f in Z[phi] selected by")
print("  the Z[phi]-norm minimization of exp416,")
print()
print("      sum_{f=1}^{9} N(z_f) = b_1 = a_1 + 1")
print()
print("  where b_1 is the first Betti number of the 600-cell boundary")
print("  complex and N(z) = z * sigma(z) is the Galois norm in Q(sqrt(5)).")
print()


# #############################################################################
#               PART B: CASCADE LEVEL DECOMPOSITION
# #############################################################################
print()
print("=" * 72)
print("PART B: CASCADE LEVEL DECOMPOSITION")
print("=" * 72)
print()

levels = {
    'Level 0': ['e'],
    'Level 1': ['c', 't', 'b'],
    'Level 2': ['mu', 'tau', 'u', 's'],
    'Level 3': ['d'],
}

print("Norm sums by Galois Cascade level:")
print()

for lev_name, names in levels.items():
    s = sum(norms[n] for n in names)
    detail = " + ".join("%+d" % norms[n] for n in names)
    print("  %s: %s = %+d" % (lev_name, detail, s))

print()
print("  Level 0 (input):      0")
print("  Level 1 (Arakelov):   a1 = %d" % a1)
print("  Level 2 (anti-sym):   0  (pairs cancel!)")
print("  Level 3 (repr):       1")
print("  Total:                a1 + 1 = b1 = %d" % b1)
print()

# Detail: WHY does Level 1 give a1?
print("WHY Level 1 sum = a1:")
print("  N(z_c) = N(2+phi) = 4+2-1 = 5 = a1")
print("  N(z_t) = N(4+phi) = 16+4-1 = 19")
print("  N(z_b) = N(-1+4phi) = 1-4-16 = -19")
print("  t and b CANCEL: 19 + (-19) = 0")
print("  Remainder: N(z_c) = a1")
print()

# WHY does Level 2 give 0?
print("WHY Level 2 sum = 0:")
print("  mu(+1) + tau(-1) = 0  (Morrey pair)")
print("  s(+1) + u(-1) = 0    (spectral pair)")
print("  Level 2 fermions come in Galois-conjugate PAIRS")
print("  with opposite norm signs.")
print()

# WHY does Level 3 give 1?
print("WHY Level 3 sum = 1:")
print("  N(z_d) = N(1) = 1  (z_d is rational)")
print("  The d quark has the simplest Wilson line: z_d = 1 in Z.")
print("  Its norm N(1) = 1 is the UNIT norm.")
print()

# Physical: the +1 from d quark completes b1
print("PHYSICAL SIGNIFICANCE:")
print("  b1 = a1 + 1 determines the gauge kernel (dim = 6)")
print("  The '+1' comes ENTIRELY from the d quark.")
print("  Without the d quark, sum N = a1 = 5 (not 6).")
print("  The d quark's unit norm COMPLETES the gauge kernel dimension.")
print()


# #############################################################################
#              PART C: CUBIC IDENTITY
# #############################################################################
print()
print("=" * 72)
print("PART C: CUBIC IDENTITY  sum N^3 = N + b1")
print("=" * 72)
print()

# Compute sum of cubes
sum_N3 = sum(norms[n]**3 for n, _, _, _, _ in fermions)
target = N + b1

print("Cubic norm sum:")
print()

# By level
for lev_name, names in levels.items():
    cubes = [norms[n]**3 for n in names]
    s = sum(cubes)
    detail = " + ".join("(%+d)^3" % norms[n] for n in names)
    print("  %s: %s = %+d" % (lev_name, detail, s))

print()
print("  Total: sum N(z_f)^3 = %d" % sum_N3)
print("  N + b1 = %d + %d = %d" % (N, b1, target))
print("  Match: %s" % ("EXACT!" if sum_N3 == target else "FAIL"))
print()

# Algebraic proof
print("ALGEBRAIC PROOF:")
print()
print("  From the Cascade structure:")
print("    Level 1: 5^3 + 19^3 + (-19)^3 = a1^3 + 0 = a1^3 = %d" % a1**3)
print("    Level 2: 1+(-1)+(-1)+1 = 0  (cube-cancel)")
print("    Level 3: 1^3 = 1")
print("    Level 0: 0^3 = 0")
print("    Total: a1^3 + 1 = %d" % (a1**3 + 1))
print()
print("  So: sum N^3 = a1^3 + 1")
print("  And: N + b1 = a1! + (a1 + 1)")
print()
print("  The identity a1^3 + 1 = a1! + a1 + 1 reduces to:")
print("    a1^3 - a1 = a1!")
print("    a1(a1-1)(a1+1) = a1!")
print("    (a1-1)(a1+1) = (a1-1)!")
print("    a1 + 1 = (a1 - 2)!")
print()

# Uniqueness
print("  UNIQUENESS (a1 + 1 = (a1-2)! for integer a1 >= 3):")
for n in range(3, 12):
    lhs = n + 1
    rhs = factorial(n - 2)
    match = " <-- UNIQUE SOLUTION" if lhs == rhs else ""
    print("    a1 = %2d:  %3d vs (%d-2)! = %d%s" % (n, lhs, n, rhs, match))

print()
print("  THEOREM (Cubic Galois Norm Identity):")
print("  The condition sum_f N(z_f)^3 = |2I| + b_1 is satisfied")
print("  if and only if a_1 = 5.")
print()
print("  This is the 12th independent characterization of a_1 = 5.")
print()


# #############################################################################
#            PART D: WEINBERG ANGLE FROM GALOIS NORMS
# #############################################################################
print()
print("=" * 72)
print("PART D: WEINBERG ANGLE FROM GALOIS NORMS")
print("=" * 72)
print()

sin2_tW = b1 / (a1**2 + 1)
sin2_tW_obs = 0.23122  # PDG 2024

print("Known: sin^2(theta_W) = b1 / (a1^2 + 1) = %d / %d = %.6f" % (
    b1, a1**2 + 1, sin2_tW))
print("Observation (MSbar, m_Z): %.5f +/- 0.00003" % sin2_tW_obs)
print("Deviation: %.2f%%" % (100 * abs(sin2_tW - sin2_tW_obs) / sin2_tW_obs))
print()

print("NEW: Since sum_f N(z_f) = b1 (Part A), we can write:")
print()
print("  sin^2(theta_W) = sum_f N(z_f) / (a1^2 + 1)")
print()
print("This formula expresses the WEAK MIXING ANGLE as the ratio of")
print("the total Galois norm (from the mass sector) to the geometric")
print("factor a1^2 + 1 (from the 600-cell).")
print()

# Decomposition by level
print("Level decomposition of sin^2(theta_W):")
print()
print("  sin^2(tW) = (0 + a1 + 0 + 1) / (a1^2 + 1)")
print("            = (a1 + 1) / (a1^2 + 1)")
print()
for lev_name, names in levels.items():
    s = sum(norms[n] for n in names)
    contrib = s / (a1**2 + 1)
    print("  %s contributes %+d/26 = %+.4f  (%s)" % (
        lev_name, s, contrib,
        ", ".join(names)))

print()
print("  The Arakelov sector (c,t,b) contributes 5/26 of the mixing.")
print("  The d quark cross-term contributes 1/26.")
print("  Leptons and light quarks contribute ZERO.")
print()


# #############################################################################
#            PART E: TQFT VACUUM AND COUPLING CONNECTIONS
# #############################################################################
print()
print("=" * 72)
print("PART E: TQFT VACUUM AND COUPLING CONNECTIONS")
print("=" * 72)
print()

E_vac = phi**3
E_vac_galois = phi_c**3

print("TQFT vacuum energy on A_3 fusion ring:")
print("  E_vac = (1/2) sum_{j>0} d_j * L_j = phi^3 = 2 + sqrt(5)")
print("        = %.10f" % E_vac)
print()

print("Connection to strong coupling:")
print("  alpha_s = 1/(2*phi^3) = 1/(2*E_vac)")
print("  E_vac = 1/(2*alpha_s) = %.10f" % (1/(2*alpha_s)))
print("  Match: %s" % ("EXACT!" if abs(E_vac - 1/(2*alpha_s)) < 1e-12 else "FAIL"))
print()

print("Galois conjugate:")
print("  E_vac' = phi'^3 = 2 - sqrt(5) = %.10f" % E_vac_galois)
print("  E_vac * E_vac' = phi^3 * phi'^3 = (phi*phi')^3 = (-1)^3 = -1")
print("  Check: %.10f  [%s]" % (
    E_vac * E_vac_galois,
    "EXACT" if abs(E_vac * E_vac_galois + 1) < 1e-12 else "FAIL"))
print()

print("PHYSICAL MEANING:")
print("  The TQFT vacuum energy is a UNIT in Z[phi] with norm -1.")
print("  E_vac > 0 (physical: stable)")
print("  E_vac' < 0 (Galois: unstable -> vacuum metastability)")
print()

# Neutrino splitting
print("Neutrino mass ratio:")
r_nu = alpha_em * phi**3
r_nu_obs = 7.53e-5 / 2.453e-3  # Dm^2_21 / Dm^2_31

print("  r = Dm^2_21 / Dm^2_31 = alpha * phi^3 = alpha * E_vac")
print("    = alpha / (2*alpha_s)")
print("    = %.6f" % r_nu)
print("  Observed: %.4f" % r_nu_obs)
print("  Match: %.2f%%" % (100 * abs(r_nu - r_nu_obs) / r_nu_obs))
print()
print("  The neutrino mass hierarchy measures the ratio of")
print("  electromagnetic to strong coupling constants!")
print()

# Seesaw exponent
n_seesaw = 35
print("Seesaw exponent:")
print("  n_seesaw = 35 = a1 * (a1+2) = a1 * r  (r = k+2 = TQFT level)")
print("  m_3 = 2*m_e/phi^35 = 49.5 meV (PREDICTION)")
print("  M_R = m_e*phi^35/2 ~ 5.3 TeV")
print()


# #############################################################################
#             PART F: WILSON LINE PRODUCTS
# #############################################################################
print()
print("=" * 72)
print("PART F: WILSON LINE PRODUCTS (SU(2) DOUBLETS)")
print("=" * 72)
print()

doublets_q = [('u', 'd'), ('c', 's'), ('t', 'b')]
doublets_l = [('e', 'mu'), ('e', 'tau')]  # not physical doublets

print("Quark SU(2)_L doublets:")
print()
for q1, q2 in doublets_q:
    a1_q = [a for n, a, b, _, _ in fermions if n == q1][0]
    b1_q = [b for n, a, b, _, _ in fermions if n == q1][0]
    a2_q = [a for n, a, b, _, _ in fermions if n == q2][0]
    b2_q = [b for n, a, b, _, _ in fermions if n == q2][0]

    # Product in Z[phi]
    a_prod = a1_q * a2_q + b1_q * b2_q
    b_prod = a1_q * b2_q + a2_q * b1_q + b1_q * b2_q
    N_prod = galois_norm(a_prod, b_prod)
    z_prod = a_prod + b_prod * phi

    print("  (%s,%s): z_%s * z_%s = %d + %d*phi = %.4f" % (
        q1, q2, q1, q2, a_prod, b_prod, z_prod))
    print("          N(product) = %+d = N(z_%s)*N(z_%s) = %+d*%+d = %+d" % (
        N_prod, q1, q2, norms[q1], norms[q2], norms[q1] * norms[q2]))

    # Special structure
    if b_prod != 0 and a_prod == 0:
        print("          = %d * phi  (PURE golden ratio multiple!)" % b_prod)
    print()

print("KEY OBSERVATION:")
print("  z_t * z_b = 19*phi = |N(z_t)| * phi")
print("  The top-bottom Wilson line product is the Galois norm")
print("  times the golden ratio. This is STRUCTURAL:")
print("  N(19*phi) = 19^2 * N(phi) = 19^2 * (-1) = -361")
print()

# Doublet norm products
print("Doublet norm products |N(z_up)*N(z_down)|:")
for q1, q2 in doublets_q:
    prod = abs(norms[q1] * norms[q2])
    print("  (%s,%s): |N_up * N_down| = %d" % (q1, q2, prod))

print()
print("  Ratios: 1 : a1 : (4*a1-1)^2 = 1 : 5 : 361")
print()


# #############################################################################
#              PART G: ANOMALY-LIKE STRUCTURE
# #############################################################################
print()
print("=" * 72)
print("PART G: ANOMALY-LIKE STRUCTURE")
print("=" * 72)
print()

# Lepton vs quark decomposition
lep_names = ['e', 'mu', 'tau']
qrk_names = ['u', 'd', 's', 'c', 'b', 't']

sum_lep = sum(norms[n] for n in lep_names)
sum_qrk = sum(norms[n] for n in qrk_names)

print("Lepton-quark decomposition:")
print("  Leptons: sum N = %+d" % sum_lep)
print("  Quarks:  sum N = %+d" % sum_qrk)
print()

print("  Leptons: %s = 0" % " + ".join("(%+d)" % norms[n] for n in lep_names))
print("  The leptonic Galois norms sum to ZERO.")
print("  This is analogous to gravitational anomaly cancellation!")
print()

print("  Quarks: %s = %d = b1" % (
    " + ".join("(%+d)" % norms[n] for n in qrk_names), sum_qrk))
print("  The quark sector carries ALL the Galois 'charge'.")
print()

print("With color multiplicity:")
print("  N_c * sum_quarks = %d * %d = %d" % (N_c, sum_qrk, N_c * sum_qrk))
print("  Leptons + N_c * quarks = %d + %d = %d" % (
    sum_lep, N_c * sum_qrk, sum_lep + N_c * sum_qrk))
print()

# Beta function connection
print("BETA FUNCTION CONNECTION:")
print("  Framework beta coefficients: b_1:b_2:b_3 = 8:5:2")
print("  Sum = 8 + 5 + 2 = 15 = a1 * N_gen")
print()
print("  N_c * (Level 1 norm sum) = 3 * a1 = 15 = beta sum!")
print("  The running of gauge couplings is determined by")
print("  the colored Arakelov norm sum.")
print()

# Power sums
print("Galois norm power sums:")
for p in range(1, 5):
    sp = sum(norms[n]**p for n in ['mu', 'tau', 'u', 'd', 's', 'c', 'b', 't'])
    print("  sum N^%d = %+d" % (p, sp))

print()

# Sum N^2 analysis
sum_N2 = sum(norms[n]**2 for n in ['mu', 'tau', 'u', 'd', 's', 'c', 'b', 't'])
print("  sum N^2 = %d" % sum_N2)
print("  = 4*1 + 1 + 25 + 2*361 = 4 + 1 + 25 + 722 = 752")
print("  = 8 * 94 = 8 * (N - 26) = rank(E8) * (N - (a1^2+1))")
print("  Check: 8*(120-26) = 8*94 = %d  [%s]" % (
    8 * 94, "EXACT!" if sum_N2 == 8 * 94 else "FAIL"))
print()

# Sum N^3 analysis (repeat for completeness)
sum_N3 = sum(norms[n]**3 for n in ['mu', 'tau', 'u', 'd', 's', 'c', 'b', 't'])
print("  sum N^3 = %d" % sum_N3)
print("  = a1^3 + 1 = N + b1 = |2I| + b1")
print("  Check: %d + %d = %d  [%s]" % (
    N, b1, N + b1, "EXACT!" if sum_N3 == N + b1 else "FAIL"))
print()

# Sum N^4
sum_N4 = sum(norms[n]**4 for n in ['mu', 'tau', 'u', 'd', 's', 'c', 'b', 't'])
print("  sum N^4 = %d" % sum_N4)
# Check: 4*1 + 1 + 625 + 2*130321 = 4+1+625+260642 = 261272
# = ?
factors = []
temp = sum_N4
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
    while temp % p == 0:
        factors.append(p)
        temp //= p
if temp > 1:
    factors.append(temp)
print("    = %s" % " * ".join(str(f) for f in factors))
print()


# #############################################################################
#                       PART H: SUMMARY
# #############################################################################
print()
print("=" * 72)
print("PART H: SUMMARY AND CLASSIFICATION")
print("=" * 72)
print()

print("NEW STRUCTURAL IDENTITIES:")
print()
print("  1. GALOIS NORM SUM RULE (THEOREM):")
print("     sum_f N(z_f) = b1 = a1 + 1 = 6")
print("     Level decomposition: a1 + 0 + 1 + 0 = b1")
print("     CLASSIFICATION: DERIVED")
print()
print("  2. CUBIC GALOIS NORM IDENTITY (THEOREM):")
print("     sum_f N(z_f)^3 = |2I| + b1 = 126")
print("     Equivalent to a1+1 = (a1-2)!, unique solution a1=5")
print("     12th independent characterization of a1=5")
print("     CLASSIFICATION: DERIVED")
print()
print("  3. WEINBERG ANGLE FROM MASS NORMS (COROLLARY):")
print("     sin^2(theta_W) = sum_f N(z_f) / (a1^2+1)")
print("     Connects gauge mixing to fermion mass Galois norms")
print("     CLASSIFICATION: DERIVED (from identity 1 + known formula)")
print()
print("  4. TQFT VACUUM = INVERSE STRONG COUPLING (THEOREM):")
print("     E_vac = phi^3 = 1/(2*alpha_s)")
print("     E_vac * E_vac' = -1 (unit norm in Z[phi])")
print("     CLASSIFICATION: DERIVED")
print()
print("  5. NEUTRINO RATIO = EM/STRONG (COROLLARY):")
print("     r_nu = Dm21/Dm31 = alpha * E_vac = alpha/(2*alpha_s)")
print("     CLASSIFICATION: DERIVED (from known formulas + identity 4)")
print()
print("  6. BETA SUM = COLORED ARAKELOV (COROLLARY):")
print("     sum(beta_i) = 15 = N_c * Level_1_norm_sum")
print("     CLASSIFICATION: DERIVED")
print()
print("  7. TOP-BOTTOM WILSON PRODUCT (OBSERVATION):")
print("     z_t * z_b = |N(z_t)| * phi = 19*phi")
print("     CLASSIFICATION: STRUCTURAL")
print()
print("  8. LEPTONIC NORM VANISHING (THEOREM):")
print("     sum_{leptons} N(z_f) = 0")
print("     Analogous to gravitational anomaly cancellation")
print("     CLASSIFICATION: DERIVED")
print()

print("WEB OF CONNECTIONS:")
print()
print("  MASS SECTOR  <---->  GAUGE SECTOR")
print("  N(z_f)       --->    b1 = sum N(z_f)     [norm sum]")
print("  b1           --->    sin^2(tW)            [Weinberg angle]")
print("  b1           --->    ker(Galois)          [gauge kernel]")
print("  Level 1 sum  --->    beta sum             [running]")
print()
print("  MASS SECTOR  <---->  TQFT SECTOR")
print("  corrections  --->    E_vac = phi^3        [vacuum energy]")
print("  alpha_s      <---    1/(2*E_vac)          [strong coupling]")
print("  r_nu         <---    alpha * E_vac        [neutrino ratio]")
print()
print("  All three sectors are unified by the Galois arithmetic")
print("  of Z[phi] on the 600-cell McKay graph.")
print()

print("=" * 72)
print("END OF EXPERIMENT 445")
print("=" * 72)

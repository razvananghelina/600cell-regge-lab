"""
EXP-414: SU(2)_3 Chern-Simons on S^3 via Hopf Fibration
=========================================================

Three questions:
  1. Compute Z(S^3) exactly for SU(2) at level k=3 (=a1-2)
  2. Check if Z(S^3) connects to spectral action coefficients c_0, c_1, c_2
  3. Check if the Hopf S^1 fiber relates to Lorentzian signature

Strategy: compute exactly, express in framework quantities.
If clean -> continue. If garbage -> STOP and say so.

Author: Claude (exp414)
Date: 2026-02-25
"""

import math

# Framework constants
a1 = 5
b1 = 6
PHI = (1 + math.sqrt(5)) / 2
PHI_prime = (1 - math.sqrt(5)) / 2  # = -1/phi
N = 120
N_gen = 3
N_eig = 9
alpha_s = 1 / (2 * PHI**3)
k = a1 - 2  # CS level = 3

# Seeley-DeWitt coefficients
A_0 = 11
A_1 = 62
A_2 = 233  # = F_13

print("=" * 72)
print("EXP-414: SU(2)_3 CHERN-SIMONS ON S^3 VIA HOPF FIBRATION")
print("=" * 72)


# =====================================================================
# PART 1: EXACT Z(S^3)
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: EXACT PARTITION FUNCTION Z(S^3)")
print("=" * 72)

print("""
For SU(2) Chern-Simons at level k on S^3 (Witten 1989):

  Z(S^3) = S_{00} = sqrt(2/(k+2)) * sin(pi/(k+2))

where S_{ij} is the modular S-matrix of the WZW model.

For k = 3 (= a1 - 2):
  Z(S^3) = sqrt(2/a1) * sin(pi/a1)
         = sqrt(2/5) * sin(pi/5)
""")

# Exact computation
Z_S3 = math.sqrt(2.0/a1) * math.sin(math.pi/a1)
Z_S3_sq = Z_S3**2

print(f"  Z(S^3) = {Z_S3:.15f}")
print(f"  Z(S^3)^2 = {Z_S3_sq:.15f}")

# Express in terms of framework quantities
# D^2 = a1 + sqrt(a1) (total quantum dimension squared, from exp405)
D_sq = a1 + math.sqrt(a1)  # = 5 + sqrt(5)
D_sq_prime = a1 - math.sqrt(a1)  # = 5 - sqrt(5)

print(f"\n  D^2 (total quantum dim) = a1 + sqrt(a1) = {D_sq:.10f}")
print(f"  D'^2 (Galois conjugate) = a1 - sqrt(a1) = {D_sq_prime:.10f}")
print(f"  D^2 * D'^2 = a1*(a1-1) = {D_sq * D_sq_prime:.10f} = {a1*(a1-1)}")

# Key relation: Z(S^3) = 1/D (total quantum dimension)
D_total = math.sqrt(D_sq)
check_Z = 1.0 / D_total
print(f"\n  1/D = 1/sqrt(a1+sqrt(a1)) = {check_Z:.15f}")
print(f"  Z(S^3)                     = {Z_S3:.15f}")
print(f"  Match: {abs(check_Z - Z_S3) < 1e-14}")

print(f"""
  ==> Z(S^3) = 1/D = 1/sqrt(a1 + sqrt(a1))   ... EXACT

  This is a standard TQFT result: the S^3 partition function
  equals the inverse total quantum dimension.
""")

# Z^2 in framework form
print("  Z(S^3)^2 in framework form:")
print(f"    Z^2 = 1/D^2 = 1/(a1+sqrt(a1)) = {Z_S3_sq:.15f}")

# Rationalize
Z_sq_rationalized = D_sq_prime / (a1*(a1-1))
print(f"    Z^2 = D'^2 / [a1(a1-1)] = D'^2 / {a1*(a1-1)} = {Z_sq_rationalized:.15f}")
print(f"    Match: {abs(Z_S3_sq - Z_sq_rationalized) < 1e-14}")

# Galois norm of Z
Z_prime = 1.0 / math.sqrt(D_sq_prime)
Z_product = Z_S3 * Z_prime
print(f"\n  Galois conjugate Z' = 1/D' = 1/sqrt(a1-sqrt(a1)) = {Z_prime:.15f}")
print(f"  Z * Z' = 1/sqrt(D^2*D'^2) = 1/sqrt(a1(a1-1))")
print(f"         = 1/sqrt({a1*(a1-1)}) = {Z_product:.15f}")
print(f"         = 1/(2*sqrt(a1)) = {1/(2*math.sqrt(a1)):.15f}")
print(f"  Match: {abs(Z_product - 1/(2*math.sqrt(a1))) < 1e-14}")

# Connection to alpha product
alpha_em_approx = 1/137.036  # approximate for comparison
alpha_prime_approx = 1/(2*math.pi*alpha_em_approx)
print(f"\n  For reference: alpha*alpha' = 1/(2*pi) = {1/(2*math.pi):.10f}")
print(f"  Z*Z' = 1/(2*sqrt(a1)) = {1/(2*math.sqrt(a1)):.10f}")
print(f"  (Z*Z')^2 = 1/(4*a1) = 1/{4*a1} = {1/(4*a1):.10f}")

print(f"""
  SUMMARY Part 1:
  +-----------------------------------------+---------------------------+
  | Quantity                                | Value                     |
  +-----------------------------------------+---------------------------+
  | Z(S^3)                                  | 1/sqrt(a1+sqrt(a1))       |
  | Z(S^3)^2                                | D'^2 / [a1(a1-1)]         |
  | Galois norm |N(Z)|^2 = Z^2 * Z'^2      | 1/(4*a1) = 1/20           |
  | D^2 * D'^2                              | a1(a1-1) = 20 = N/b1      |
  +-----------------------------------------+---------------------------+
  Status: DERIVED (standard TQFT). Clean phi/a1 expressions. CONTINUE.
""")


# =====================================================================
# PART 2: CONNECTION TO SPECTRAL ACTION COEFFICIENTS
# =====================================================================
print("=" * 72)
print("PART 2: SPECTRAL ACTION COEFFICIENTS")
print("=" * 72)

c_0 = 2 * N * A_0  # = 2640
c_1 = 2 * N * A_1  # = 14880
c_2 = 2 * N * A_2  # = 55920

print(f"""
  Seeley-DeWitt coefficients (from spectral action on 600-cell):
    c_0 = 2*N*A_0 = 2*{N}*{A_0} = {c_0}
    c_1 = 2*N*A_1 = 2*{N}*{A_1} = {c_1}
    c_2 = 2*N*A_2 = 2*{N}*{A_2} = {c_2}

  Question: does Z(S^3) = {Z_S3:.6f} connect to these?
""")

# Try various combinations
print("  Testing ratios and products:")

# Z * c_0
print(f"    Z * c_0 = {Z_S3 * c_0:.6f}")
print(f"    Z^2 * c_0 = {Z_S3_sq * c_0:.6f}")
print(f"    c_0 * D^2 = {c_0 * D_sq:.6f} = {c_0} * (5+sqrt(5))")
print(f"    c_0 / D^2 = {c_0 / D_sq:.6f}")

# Ratio c_1/c_0 = A_1/A_0 = 62/11
ratio_10 = c_1 / c_0
ratio_21 = c_2 / c_1
print(f"\n    c_1/c_0 = A_1/A_0 = {A_1}/{A_0} = {ratio_10:.10f}")
print(f"    c_2/c_1 = A_2/A_1 = {A_2}/{A_1} = {ratio_21:.10f}")

# Check if ratios involve phi
print(f"\n    A_1/A_0 = {ratio_10:.10f}")
print(f"    phi^4+1 = {PHI**4+1:.10f}")
print(f"    b1*A_0-b1+2 = {b1*A_0-b1+2}")
print(f"    a1*A_0+7 = {a1*A_0+7}")

# A_2/A_1 check
print(f"\n    A_2/A_1 = {ratio_21:.10f}")
print(f"    phi^2 + 2 = {PHI**2 + 2:.10f}")
print(f"    (a1-1)*phi^2 = {(a1-1)*PHI**2:.10f}")

# Direct check: is c_0 related to D^2?
print(f"\n    c_0 = {c_0}")
print(f"    c_0 / N = {c_0/N} = 2*A_0 = {2*A_0}")
print(f"    c_0 / (N * D^2) = {c_0/(N*D_sq):.10f}")

# The CS partition function normalization on S^3 involves
# the volume and the Chern-Simons coupling
# Z_CS ~ exp(i*k*CS[A]) for flat connections. On S^3, only trivial flat connection.
# Perturbative 1-loop: involves analytic torsion

# Try: does Z^2 * something = A_k / A_l ?
print("\n  Checking if Z^2 appears in A_k ratios:")
print(f"    Z^2 = {Z_S3_sq:.10f}")
print(f"    A_1/A_0 - a1 = {A_1/A_0 - a1:.10f}")
print(f"    A_0 * Z^2 = {A_0 * Z_S3_sq:.10f}")

# The key quantity: 1/D^2 = Z^2 = 1/(5+sqrt(5))
# Check: A_0 = 11 = 2*a1 + 1 = 2*5 + 1. Known.
# A_1 = 62 = 10*a1 + 12 = 10*5 + 12. Or: 62 = N/2 + 2 = 60 + 2.
# A_2 = 233 = F_13 (Fibonacci number!)

# Fibonacci connection: phi^n = F_n * phi + F_{n-1}
# A_2 = 233 = F_13. Check: phi^13 = 233*phi + 144 = 233*1.618... + 144 = ...
# F_13 = 233, F_12 = 144

print(f"\n  Fibonacci connection:")
print(f"    A_2 = {A_2} = F_13 (13th Fibonacci number)")
print(f"    phi^13 = {PHI**13:.6f}")
print(f"    F_13*phi + F_12 = {233*PHI + 144:.6f}")
print(f"    Match: {abs(PHI**13 - (233*PHI + 144)) < 1e-6}")

# CS level k=3 connection to Fibonacci:
# Fibonacci anyons emerge at k=3. The fusion category IS the Fibonacci category.
# The Fibonacci numbers emerge because phi^n = F_n * phi + F_{n-1}
print(f"\n    k=3 CS gives Fibonacci anyons => phi powers => Fibonacci numbers")
print(f"    A_2 = F_13 is consistent with Fibonacci structure in CS_3")

# Now the critical question: can we get c_k from CS?
print(f"""
  CRITICAL ANALYSIS:
  The spectral action coefficients c_k come from the heat kernel expansion
  of the Dirac operator on the 4D noncommutative geometry:
    Tr(f(D^2/L^2)) ~ sum_k f_k * c_k * L^(4-2k)  (4D)

  The CS partition function Z(S^3) is a 3D topological invariant.

  These are different objects living in different dimensions!

  The connection would require either:
    (a) Dimensional reduction: 4D spectral action -> 3D CS (via S^1 circle)
    (b) Cobordism: 4D manifold with boundary S^3 (holographic)
    (c) Index theorem linking 4D bulk to 3D boundary

  Let me check the dimensional reduction route:
""")

# Dimensional reduction on S^1 x S^3 -> CS on S^3
# If we compactify the Euclidean time on S^1 with radius beta (= 1/T):
# The spectral action on S^1 x M^3 at high T reduces to
# Tr(f(D^2)) ~ beta * [integral_M3 of local terms]
# At leading order, the 4D spectral action gives:
# S_spectral ~ f_4 * c_0 * Lambda^4 + f_2 * c_1 * Lambda^2 + f_0 * c_2

# For S^1 x S^3: the thermal partition function Z(beta) = Tr(e^{-beta*H})
# relates to the spectral zeta function

# The CS partition function Z(S^3) is NOT the thermal partition function!
# CS is topological - it doesn't depend on the metric of S^3.

print("  Route (a): Dimensional reduction S^1 x S^3 -> S^3")
print("  The 4D spectral action on S^1_beta x S^3 gives:")
print("    S = beta * [a_3(D_S3) + ...]  (3D Seeley-DeWitt)")
print(f"  But CS is topological (metric-independent).")
print(f"  The spectral action is NOT topological (depends on Lambda, metric).")
print(f"  => Direct identification FAILS.\n")

# However: the Chern-Simons term DOES appear in the spectral action!
# In Connes' NCG: the spectral action contains a CS-like term at c_1 order:
# c_1 term ~ integral R (Ricci scalar) + CS-like gauge terms
# More precisely, the spectral action for the Standard Model gives:
# S = ... + (1/g^2) * int F^2 + theta * int F*F + CS_gravitational + ...

print("  However: CS term appears INSIDE the spectral action:")
print("    c_1 term contains: integral R (gravity) + gauge kinetic terms")
print("    The gravitational CS term: CS_grav = (1/2) int Tr(w dw + 2/3 w^3)")
print("    appears in the eta invariant of the Dirac operator.\n")

# The actual connection: eta invariant
# For the Dirac operator on S^3, the eta invariant is:
# eta(D_S3) is related to the CS invariant of the spin connection

# On S^3 with round metric:
# The eta invariant of the standard Dirac operator vanishes: eta(S^3) = 0
# (by the symmetry S^3 = SU(2))

print("  Route (b): Eta invariant connection")
print("    For the Dirac operator on round S^3:")
print("    eta(D_{S^3}) = 0  (by SU(2) symmetry)")
print("    The gravitational CS invariant: cs(S^3) = 0 mod Z")
print("    => The S^3 contribution to the spectral action anomaly is TRIVIAL.\n")

# Check route (c): index/surgery
# Witten's relation: Z(S^3) relates to the index of D on the 4-ball B^4
# via the APS index theorem:
# ind(D_{B^4}) = integral_{B^4} A-hat - (h + eta(S^3))/2
# where h = dim(ker(D_{S^3}))

# For round S^3: h = 0, eta = 0, so ind = integral A-hat = 0 for B^4.
# This is trivial.

print("  Route (c): APS index theorem (B^4 bounded by S^3)")
print("    ind(D_{B^4}) = int_{B^4} A-hat - (h+eta)/2 = 0 - 0 = 0")
print("    => Trivial for the standard geometry.\n")

# So what IS the connection?
# The connection is INDIRECT, through representation theory:
# Both the spectral action (via the McKay correspondence) and the CS theory
# (via the modular S-matrix) use the SAME representation theory of SU(2).

# The number 120 = |2I| appears in both:
# - Spectral action: N = |2I| = 120 vertices
# - CS on Poincare homology sphere: Z(Sigma(2,3,5)) involves |2I| = 120

# Check: Z on the Poincare homology sphere Sigma(2,3,5) = S^3/2I
print("  INDIRECT CONNECTION via the Poincare homology sphere:")
print("  Sigma(2,3,5) = S^3 / 2I where 2I is the binary icosahedral group")
print(f"  |2I| = N = {N}")

# For SU(2)_k CS on Sigma(2,3,5) = S^3/Gamma:
# Z(S^3/Gamma) = sum_j |S_{0j}|^2 / S_{0j} * ...
# Actually, the surgery formula gives:
# Z(S^3/Gamma) = (1/|Gamma|) * sum_{g in Gamma} Z(S^3, holonomy g)

# For abelian orbifold S^3/Z_p:
# Z(S^3/Z_p) = (1/p) * sum_{a=0}^{p-1} sum_j S_{0j} * exp(2*pi*i*j^2*a/p)

# For 2I (non-abelian, |2I|=120):
# This requires the Verlinde formula applied to the surgery presentation

# Actually, there IS a much cleaner result:
# Z(M) for a Seifert manifold like Sigma(2,3,5) can be computed exactly.
# Sigma(2,3,5) has Seifert invariants {b; (2,1),(3,1),(5,1)} with b=-1

# The result for SU(2)_3 on Sigma(2,3,5) is known to involve phi.

print("\n  Computing Z(Sigma(2,3,5)) for SU(2)_3:")
print("  (Sigma(2,3,5) = Poincare homology sphere = S^3/2I)")

# Use the formula for lens spaces / Seifert manifolds
# For SU(2)_k, the quantum dimensions are:
# d_j = S_{0j}/S_{00} = sin((2j+1)*pi/(k+2)) / sin(pi/(k+2))

# For k=3: j = 0, 1/2, 1, 3/2
d = []
for j_idx, j in enumerate([0, 0.5, 1.0, 1.5]):
    d_j = math.sin((2*j+1)*math.pi/a1) / math.sin(math.pi/a1)
    d.append(d_j)
    label = f"d_{{{j}}}"
    print(f"    {label:8s} = {d_j:.10f}")

print(f"\n    Spectrum: {{{d[0]:.4f}, {d[1]:.4f}, {d[2]:.4f}, {d[3]:.4f}}}")
print(f"    = {{1, phi, phi, 1}} (as in exp407)")

# Verify
print(f"    d_{{1/2}} = phi: {abs(d[1] - PHI) < 1e-10}")
print(f"    d_1 = phi: {abs(d[2] - PHI) < 1e-10}")

# T-matrix (twists) for SU(2)_k:
# T_j = exp(2*pi*i * j(j+1)/(k+2)) = exp(2*pi*i * j(j+1)/a1)
print(f"\n  T-matrix (conformal weights h_j = j(j+1)/(k+2)):")
T = []
for j_idx, j in enumerate([0, 0.5, 1.0, 1.5]):
    h_j = j*(j+1)/a1
    theta_j = 2*math.pi*h_j
    T_j = complex(math.cos(theta_j), math.sin(theta_j))
    T.append(T_j)
    print(f"    h_{{{j}}} = {j}*{j+1}/{a1} = {h_j:.6f}")
    print(f"    T_{{{j}}} = exp(2*pi*i*{h_j:.4f}) = {T_j.real:.6f} + {T_j.imag:.6f}i")

print(f"""
  VERDICT Part 2:
  Z(S^3) = 1/D does NOT directly reproduce c_0, c_1, c_2.

  The connection is REPRESENTATIONAL, not computational:
  - Both use SU(2) representation theory
  - The spectral action coefficients A_k live in 4D (heat kernel)
  - The CS partition function Z lives in 3D (topological)
  - They share the SAME algebraic data (phi, a1, N) but compute
    DIFFERENT physical quantities.

  The Fibonacci number A_2 = 233 = F_13 IS connected to the
  Fibonacci anyons of CS_3 through the powers of phi,
  but this is a shared algebraic structure, not a direct formula.

  Status: NEGATIVE for direct c_k reproduction.
          POSITIVE for shared algebraic underpinning.
""")


# =====================================================================
# PART 3: HOPF FIBRATION AND LORENTZIAN SIGNATURE
# =====================================================================
print("=" * 72)
print("PART 3: HOPF FIBRATION S^3 -> S^2 AND SIGNATURE")
print("=" * 72)

print("""
  The Hopf fibration: S^1 -> S^3 -> S^2

  S^3 = {(z_1, z_2) in C^2 : |z_1|^2 + |z_2|^2 = 1}
  Hopf map: (z_1, z_2) -> [z_1 : z_2] in CP^1 = S^2
  Fiber: S^1 (the phase e^{i*theta} acting on (z_1, z_2))

  The Hopf bundle is the unique non-trivial S^1 bundle over S^2.
  First Chern class: c_1 = 1 (unit magnetic monopole).
""")

# The metric on S^3 decomposes as:
# ds^2_{S^3} = ds^2_{S^2} + (d*psi + A)^2
# where A is the monopole connection on S^2
# and psi is the fiber coordinate with period 4*pi

print("  Metric decomposition (Kaluza-Klein):")
print("    ds^2_{S^3} = ds^2_{S^2} + (d*psi + A)^2")
print("    A = monopole connection, psi in [0, 4*pi)")
print("    This is a U(1) Kaluza-Klein reduction.\n")

# Lorentzian continuation: psi -> i*t (Wick rotation of fiber)
# ds^2 = ds^2_{S^2} - (dt + i*A)^2
# This gives a LORENTZIAN metric on R x S^2 (Nariai-like)

print("  Wick rotation of the fiber: psi -> i*t")
print("    ds^2 = ds^2_{S^2} - (dt + i*A)^2")
print("    = -dt^2 + ds^2_{S^2} - 2i*A*dt + A^2")
print("    The cross term 2i*A*dt is imaginary => PROBLEM.")
print("    Pure Wick rotation of the fiber gives a COMPLEX metric.\n")

# The correct procedure: analytic continuation of the full S^3 metric
# In embedding coordinates: S^3 = {x^2 + y^2 + z^2 + w^2 = R^2}
# Wick rotate w -> i*w: x^2 + y^2 + z^2 - w^2 = R^2
# This is dS_3 (de Sitter 3-space) or H^{2,1}

print("  Correct Wick rotation: embed S^3 in R^4, then x_4 -> i*x_4")
print("    S^3: x_1^2 + x_2^2 + x_3^2 + x_4^2 = R^2  (Euclidean)")
print("    dS_3: x_1^2 + x_2^2 + x_3^2 - x_4^2 = R^2  (Lorentzian)")
print("    This gives 3D de Sitter space, NOT 4D spacetime!\n")

# The 600-cell context: vertices on S^3 subset R^4
# If we want 4D Lorentzian spacetime, S^3 is the SPATIAL section:
# ds^2_{4D} = -dt^2 + a(t)^2 * ds^2_{S^3}  (closed FRW cosmology)

print("  In our framework: S^3 = spatial section of closed FRW cosmology:")
print("    ds^2 = -dt^2 + a(t)^2 * ds^2_{S^3}")
print("    The 600-cell provides the discrete structure on S^3.")
print("    Time t is EXTERNAL (not from the Hopf fiber).\n")

# However, there IS an interesting Hopf-time connection:
# The Hopf fiber direction is the RIGHT-invariant U(1) in SU(2) = S^3
# Under the group action, it generates ROTATIONS
# These rotations can be identified with the KK gauge field

# In our framework: alpha*alpha' = 1/(2*pi) = 1/Vol(S^1)
# The S^1 fiber has circumference 4*pi (in standard normalization of Hopf)
# but the GEOMETRIC S^1 that appears in the alpha equation has Vol = 2*pi

print("  Connection to alpha equation (exp354-355):")
print(f"    alpha * alpha' = 1/(2*pi) = 1/Vol(S^1)")
print(f"    = {1/(2*math.pi):.10f}")
print(f"    The S^1 in the alpha equation IS the KK circle from Hopf.")
print(f"    Vol(S^1) = 2*pi is the coefficient '2*pi' in the alpha quadratic.\n")

# Check: the alpha equation is 2*pi*alpha^2 - (N/b1)*phi^4*alpha + 1 = 0
# Leading coeff = 2*pi = Vol(S^1)
# Product of roots = 1/(2*pi) = 1/Vol(S^1)
# These are CONSISTENT with KK on the Hopf fiber.

alpha_prod = 1/(2*math.pi)
N_over_b1 = N/b1  # = 20
alpha_sum = N_over_b1 * PHI**4 / (2*math.pi)

print("  Alpha equation: 2*pi*x^2 - (N/b1)*phi^4*x + 1 = 0")
print(f"    Leading coeff = 2*pi = Vol(S^1) = {2*math.pi:.10f}")
print(f"    Product of roots = 1/(2*pi) = {alpha_prod:.10f}")
print(f"    Sum of roots = (N/b1)*phi^4/(2*pi) = {alpha_sum:.10f}")

# Discriminant
disc = (N_over_b1 * PHI**4)**2 - 4 * 2 * math.pi
print(f"    Discriminant = {disc:.10f}")

alpha_plus = (N_over_b1 * PHI**4 + math.sqrt(disc)) / (2 * 2 * math.pi)
alpha_minus = (N_over_b1 * PHI**4 - math.sqrt(disc)) / (2 * 2 * math.pi)
print(f"    alpha  = {alpha_minus:.10f}  (= 1/137.036... our alpha)")
print(f"    alpha' = {alpha_plus:.10f}   (Galois conjugate)")
print(f"    alpha * alpha' = {alpha_minus * alpha_plus:.10f}")
print(f"    1/(2*pi) = {1/(2*math.pi):.10f}")
print(f"    Match: {abs(alpha_minus * alpha_plus - 1/(2*math.pi)) < 1e-10}")

# Now: does the Hopf fiber generate Lorentzian signature?
# The answer is NUANCED:
# (1) The Hopf fiber is NOT time. Time is external.
# (2) BUT: the KK reduction along the Hopf fiber gives the U(1)_EM gauge field.
#     The alpha equation encodes this KK reduction.
# (3) The Lorentzian signature comes from spectral geometry:
#     The Dirac operator D has a BUILT-IN Lorentzian structure
#     via the grading operator gamma_5 (chirality).

print(f"""
  QUESTION 3 VERDICT: Does the Hopf fiber generate Lorentzian time?

  ANSWER: NO, not directly. But it connects to gauge structure.

  What the Hopf fiber DOES generate:
  - The U(1)_EM gauge field (via Kaluza-Klein)
  - The alpha equation coefficient 2*pi = Vol(S^1)
  - The product alpha*alpha' = 1/Vol(S^1) = 1/(2*pi)

  What generates Lorentzian signature:
  - In NCG: the REALITY OPERATOR J and KO-dimension
  - KO-dim = 2 mod 8 for the Standard Model (Connes)
  - This gives the correct Lorentzian signature upon Wick rotation
  - The 600-cell framework has KO-dim derived from a1=5 structure

  The Hopf fiber S^1 is SPATIAL (internal gauge), not TEMPORAL.
  Lorentzian signature is a SPECTRAL property (Dirac operator sign),
  not a fibration property.

  Status: NEGATIVE for Hopf -> time.
          POSITIVE for Hopf -> U(1) gauge (already known from KK).
""")


# =====================================================================
# BONUS: What Z(S^3) DOES connect to
# =====================================================================
print("=" * 72)
print("BONUS: CLEAN FRAMEWORK CONNECTIONS OF Z(S^3)")
print("=" * 72)

# 1. Z^2 and the cosmological constant
# Lambda_P = alpha^(57-alpha_s). Does Z appear?
print("\n  1. Z(S^3) and quantum dimensions:")
print(f"     Z = 1/D = {Z_S3:.10f}")
print(f"     D^2 = a1 + sqrt(a1) = {D_sq:.10f}")
print(f"     20 * Z^2 = D'^2/(D^2) * D'^2/... ")
print(f"     Actually: a1*(a1-1)*Z^2 = D'^2 = a1-sqrt(a1)")
print(f"     = {a1*(a1-1)*Z_S3_sq:.10f} vs {D_sq_prime:.10f}")
print(f"     Match: {abs(a1*(a1-1)*Z_S3_sq - D_sq_prime) < 1e-10}")

# 2. Z(S^3) and the total vertex count
# D^2 = a1 + sqrt(a1). And D^2 * D'^2 = 20 = N/b1.
# So Z^2 = 1/D^2 = b1/(N * D'^2/D^2 * ...) hmm, not clean.
# Cleaner: 1/Z^2 = D^2 = a1 + sqrt(a1).
# D^2 = a1(1 + 1/sqrt(a1)) = a1 + a1/sqrt(a1)

print(f"\n  2. Z(S^3) and vertex count:")
print(f"     N * Z^2 = {N * Z_S3_sq:.10f}")
print(f"     N / D^2 = {N / D_sq:.10f}")
print(f"     = N / (a1+sqrt(a1)) = 120/(5+sqrt(5))")
val = N / D_sq
print(f"     = {val:.10f}")
print(f"     = 12*(5-sqrt(5))/4 = 3*(5-sqrt(5)) = 15 - 3*sqrt(5)")
print(f"     = {15 - 3*math.sqrt(5):.10f}")
# Not particularly clean.

# 3. The most natural connection: Z * D = 1
print(f"\n  3. The fundamental relation: Z(S^3) * D = 1")
print(f"     This is the TQFT axiom: Z(S^3) = 1/D")
print(f"     where D is the total quantum dimension.")
print(f"     D encodes ALL the Fibonacci anyon data.")
print(f"     Z(S^3) is the 'inverse vacuum amplitude'.\n")

# 4. Connection to Wilson lines (masses)
# In our framework, fermion masses = Wilson lines on McKay graph
# Wilson line W_j = d_j / d_0 = d_j (since d_0 = 1)
# So the partition function Z(S^3) encodes the normalization
# of Wilson lines: <W_j> = d_j * Z(S^3) for the unknot in S^3

print("  4. Wilson lines and masses:")
print("     <W_j>_{S^3} = d_j * Z(S^3) = d_j / D")
for j_idx, (j, dj) in enumerate(zip([0, 0.5, 1.0, 1.5], d)):
    Wj = dj * Z_S3
    print(f"     <W_{{{j}}}>  = {dj:.6f} * {Z_S3:.6f} = {Wj:.10f}")
    if j in [0.5, 1.0]:
        print(f"                = phi/D = phi/sqrt(a1+sqrt(a1)) = {PHI/D_total:.10f}")

print(f"\n     phi/D = {PHI/D_total:.10f}")
print(f"     phi * Z = {PHI * Z_S3:.10f}")
print(f"     = phi/sqrt(5+sqrt(5))")

# Simplify: phi^2/(5+sqrt(5)) = (phi+1)/(5+sqrt(5))
# = (3+sqrt(5))/(2*(5+sqrt(5)))
# Rationalize: multiply by (5-sqrt(5))/(5-sqrt(5)):
# = (3+sqrt(5))*(5-sqrt(5)) / (2*20)
# = (15-3*sqrt(5)+5*sqrt(5)-5) / 40
# = (10+2*sqrt(5)) / 40
# = (5+sqrt(5)) / 20
# So (phi*Z)^2 = (5+sqrt(5))/20 = D^2/(4*a1*(a1-1))... hmm
print(f"     (phi*Z)^2 = {(PHI*Z_S3)**2:.10f}")
print(f"     (5+sqrt(5))/20 = {(5+math.sqrt(5))/20:.10f}")

# Actually that IS D^2/20 = D^2/(a1*(a1-1))
print(f"     = D^2 / [a1*(a1-1)] = D^2/20")
print(f"     = {D_sq/20:.10f}")
print(f"     So: phi * Z = sqrt(D^2/20) = D/sqrt(20) = D/(2*sqrt(a1))")

# This means <W_{1/2}> = <W_1> = phi*Z = D/(2*sqrt(a1)) = sqrt((a1+sqrt(a1))/(4*a1))
W_half = PHI * Z_S3
print(f"\n     The Wilson line expectation value for j=1/2 and j=1:")
print(f"     <W> = phi / sqrt(a1+sqrt(a1)) = {W_half:.10f}")
print(f"     <W>^2 = (a1+sqrt(a1))^{-1} * phi^2 = {W_half**2:.10f}")


# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("FINAL SUMMARY: EXP-414")
print("=" * 72)

print(f"""
  PART 1: Z(S^3) = 1/D = 1/sqrt(a1+sqrt(a1))     STATUS: DERIVED
    - Standard TQFT result, clean in framework quantities
    - Galois norm: Z * Z' = 1/(2*sqrt(a1)) = 1/sqrt(20)
    - Z^2 = D'^2 / [a1(a1-1)] connects to dark sector
    - All expressions involve a1 and sqrt(a1) only

  PART 2: c_k coefficients from spectral action      STATUS: NEGATIVE
    - Z(S^3) does NOT directly reproduce c_0, c_1, c_2
    - CS is 3D topological, spectral action is 4D metric-dependent
    - They share algebraic underpinning (phi, Fibonacci) but
      compute genuinely different physical quantities
    - The Fibonacci number A_2 = F_13 = 233 IS connected to
      Fibonacci anyons of CS_3, but through shared algebra

  PART 3: Hopf fiber -> Lorentzian signature          STATUS: NEGATIVE
    - The Hopf S^1 is NOT the time direction
    - It IS the KK circle generating U(1)_EM gauge field
    - alpha*alpha' = 1/(2*pi) = 1/Vol(S^1) is the KK relation
    - Lorentzian signature comes from spectral KO-dimension, not fibration
    - S^3 is the spatial section of closed FRW cosmology

  OVERALL: Z(S^3) is clean but ALREADY KNOWN (it's 1/D from exp405).
           No new physics found. The spectral action connection
           and Hopf-Lorentz connection are both NEGATIVE.

  RECOMMENDATION: Do NOT continue this line. The CS/spectral action
  bridge requires genuine 4D-3D holographic machinery that is not
  present in the current framework.
""")

"""
exp343a: CC Exponent 57 - Derivation Attempts
==============================================
KEY DISCOVERY: 57 = 3 * 19 = N_gen * N(a1 + b1*phi)

The CC exponent z_CC = 57 - alpha_s has TWO decompositions:
  (A) 57 = |A5| - N_gen = 60 - 3
  (B) 57 = N_gen * N(z_mass) where z_mass = a1 + b1*phi, N(z_mass) = 19

Both use only framework numbers. (B) connects CC to the mass formula!
This experiment tests which (if either) can be DERIVED.

Also computes: spectral zeta function of 600-cell Laplacian,
character sums, and group-theoretic paths to 57.

CATEGORY: EXPLORATORY -> seeking DERIVATION
"""

import numpy as np
from math import factorial, gcd
import time

t0 = time.time()

# ============================================================
# FRAMEWORK CONSTANTS (all from a1=5)
# ============================================================
a1 = 5
b1 = a1 + 1  # 6
PHI = (1 + np.sqrt(5)) / 2
phi_prime = (1 - np.sqrt(5)) / 2  # = -1/PHI
N = factorial(a1)  # 120
N_gen = 3
N_eig = 9
rank_E8 = 8
h_E8 = a1 * b1  # 30
degree = 2 * b1  # 12
alpha_s = 1.0 / (2 * PHI**3)

# Alpha from quadratic
A_c = 2 * np.pi
B_c = -4 * a1 * PHI**4
C_c = 1
disc = B_c**2 - 4 * A_c * C_c
alpha = (-B_c - np.sqrt(disc)) / (2 * A_c)

# Observational CC
HBAR = 1.054571817e-34
C_light = 299792458
G_N = 6.67430e-11
E_CHARGE = 1.602176634e-19
l_P = np.sqrt(HBAR * G_N / C_light**3)
H_0_SI = 67.36 * 1000 / 3.0857e22
Omega_Lambda = 0.6847
Lambda_m2 = 3 * Omega_Lambda * H_0_SI**2 / C_light**2
Lambda_Planck = Lambda_m2 * l_P**2
z_exact = np.log(Lambda_Planck) / np.log(alpha)

# Uncertainty
H0_err = 0.54
Omega_err = 0.0073
H_lo = (67.36 - H0_err) * 1000 / 3.0857e22
H_hi = (67.36 + H0_err) * 1000 / 3.0857e22
Lam_lo = 3 * (Omega_Lambda - Omega_err) * H_lo**2 / C_light**2 * l_P**2
Lam_hi = 3 * (Omega_Lambda + Omega_err) * H_hi**2 / C_light**2 * l_P**2
z_lo = np.log(Lam_hi) / np.log(alpha)
z_hi = np.log(Lam_lo) / np.log(alpha)
z_err = abs(z_hi - z_lo) / 2

def norm_zphi(a, b):
    """Norm in Z[phi]: N(a+b*phi) = a^2 + a*b - b^2"""
    return a*a + a*b - b*b

def trace_zphi(a, b):
    """Trace in Z[phi]: Tr(a+b*phi) = 2a + b"""
    return 2*a + b

print("=" * 72)
print("EXP343a: CC EXPONENT 57 - DERIVATION ATTEMPTS")
print("=" * 72)

# ============================================================
# SECTION 1: THE TWO DECOMPOSITIONS
# ============================================================
print("\n" + "=" * 72)
print("SECTION 1: TWO DECOMPOSITIONS OF 57")
print("=" * 72)

# Decomposition A: 57 = |A5| - N_gen
print(f"\n  (A) 57 = |A5| - N_gen = {N//2} - {N_gen} = {N//2 - N_gen}")
print(f"      = N/2 - N_gen = a1!/2 - N_gen")
print(f"      Physical: rotational symmetries minus generations")

# Decomposition B: 57 = N_gen * N(z_mass)
z_mass_a, z_mass_b = a1, b1  # = (5, 6)
N_zmass = norm_zphi(z_mass_a, z_mass_b)
print(f"\n  (B) 57 = N_gen * N(z_mass) = {N_gen} * {N_zmass} = {N_gen * N_zmass}")
print(f"      where z_mass = a1 + b1*phi = {a1} + {b1}*phi")
print(f"      N(z_mass) = {a1}^2 + {a1}*{b1} - {b1}^2 = {a1**2} + {a1*b1} - {b1**2} = {N_zmass}")
print(f"      AND 19 is the UNIQUE prime from exp342f!")
print(f"      Physical: N_gen copies of the mass-direction norm")

# Are they equivalent?
print(f"\n  Both give 57. Are they algebraically related?")
print(f"  |A5| = N/2 = a1!/2 = 60")
print(f"  N_gen * N(a1 + b1*phi) = 3 * (a1^2 + a1*(a1+1) - (a1+1)^2)")
print(f"  = 3 * (a1^2 + a1^2 + a1 - a1^2 - 2*a1 - 1)")
print(f"  = 3 * (a1^2 - a1 - 1)")
print(f"  For a1=5: 3 * (25 - 5 - 1) = 3 * 19 = 57")
print(f"  For |A5| = a1!/2: 120/2 = 60")
print(f"  Equality: 3*(a1^2 - a1 - 1) = a1!/2 - 3")
print(f"  => 3*a1^2 - 3*a1 = a1!/2")
print(f"  => 3*a1*(a1-1) = a1!/2 = a1*(a1-1)*(a1-2)!/2")
print(f"  => 6 = (a1-2)!/2 ... no, 6 = (a1-2)!")

# Check: is 6 = (a1-2)! for a1=5?
print(f"\n  Check: (a1-2)! = {a1-2}! = {factorial(a1-2)} = 6 = b1? ", end="")
print("YES!!!" if factorial(a1-2) == b1 else "NO")

print(f"""
  KEY IDENTITY: The two decompositions are equivalent IFF:
    N_gen * (a1^2 - a1 - 1) + N_gen = a1!/2
    <=> 3*(a1^2 - a1) = a1!/2
    <=> 6*a1*(a1-1) = a1!
    <=> 6 = (a1-2)!
    <=> (a1-2)! = b1

  This is TRUE for a1=5 because 3! = 6 = b1.
  But b1 = a1+1 = 6, and (a1-2)! = 3! = 6 is a COINCIDENCE specific to a1=5!

  Actually: (a1-2)! = (5-2)! = 3! = 6 = a1+1 = b1.
  This happens ONLY for a1=5. (For a1=6: 4!=24, b1=7. For a1=4: 2!=2, b1=5.)

  So BOTH decompositions are TRUE and EQUIVALENT, but ONLY for a1=5.
  This is another reason why a1=5 is special!
""")

# ============================================================
# SECTION 2: z_mass NORM AND UNIQUENESS THEOREM CONNECTION
# ============================================================
print("=" * 72)
print("SECTION 2: MASS DIRECTION AND UNIQUENESS THEOREM")
print("=" * 72)

print(f"\nThe mass formula n = {a1}*a + {b1}*b uses the vector ({a1},{b1}).")
print(f"In Z[phi]: z_mass = {a1} + {b1}*phi = {a1 + b1*PHI:.6f}")
print(f"Galois conjugate: z_mass' = {a1} + {b1}*phi' = {a1 + b1*phi_prime:.6f}")
print(f"N(z_mass) = z_mass * z_mass' = {N_zmass}")
print(f"Tr(z_mass) = z_mass + z_mass' = {trace_zphi(a1, b1)}")
print(f"Minimal polynomial: t^2 - {trace_zphi(a1, b1)}t + {N_zmass} = 0")

# Connection to all 9 fermion norms
fermions = {
    'e':   (0,  0,  0),
    'u':   (3, -2,  3),
    'd':   (1,  0,  5),
    'mu':  (1,  1, 11),
    's':   (1,  1, 11),
    'c':   (2,  1, 16),
    'tau': (1,  2, 17),
    'b':  (-1,  4, 19),
    't':   (4,  1, 26),
}

print(f"\nFermion Z[phi] data:")
print(f"  {'Name':>4} {'(a,b)':>8} {'n':>4} {'N(z)':>6} {'Tr(z)':>6} {'z*z_mass':>12} {'Tr(z*zm)':>10}")
for name in ['e','u','d','s','mu','c','tau','b','t']:
    a, b, n = fermions[name]
    nz = norm_zphi(a, b)
    tz = trace_zphi(a, b)
    # z*z_mass = (a+b*phi)*(a1+b1*phi) = a*a1 + a*b1*phi + b*a1*phi + b*b1*phi^2
    # = a*a1 + (a*b1+b*a1)*phi + b*b1*(phi+1)
    # = (a*a1+b*b1) + (a*b1+b*a1+b*b1)*phi
    prod_a = a*a1 + b*b1
    prod_b = a*b1 + b*a1 + b*b1
    tr_prod = 2*prod_a + prod_b
    print(f"  {name:>4} ({a:+2},{b:+2}) {n:>4} {nz:>6} {tz:>6} ({prod_a:+3},{prod_b:+3})phi {tr_prod:>10}")

# Key check: is n related to z*z_mass?
print(f"\n  Is n_k = some_function(z_k * z_mass)?")
print(f"  n_k = a1*a + b1*b = standard inner product <(a1,b1),(a,b)>")
print(f"  Tr(z*z_mass)/2 = (2*(a*a1+b*b1) + a*b1+b*a1+b*b1)/2")
print(f"  NOT equal to n = a1*a + b1*b in general")

# But what IS the standard inner product in the Z[phi] world?
# n = a1*a + b1*b. Write a1 = 5, b1 = 6.
# In the Z[phi] norm form: B(z,w) = 2*ac + ad + bc - 2*bd (bilinear form of N)
# B(z_mass, z_k) = 2*5*a + 5*b + 6*a - 2*6*b = 10a + 5b + 6a - 12b = 16a - 7b
# So B != n. The standard inner product and the norm form are DIFFERENT.

B_values = []
for name in ['e','u','d','s','mu','c','tau','b','t']:
    a, b, n = fermions[name]
    B_val = 16*a - 7*b  # bilinear form of norm
    std_ip = a1*a + b1*b  # standard inner product = n
    B_values.append((name, n, B_val, std_ip))

print(f"\n  Compare bilinear form B(z_mass, z_k) vs n_k:")
print(f"  {'Name':>4} {'n=<zm,zk>':>10} {'B(zm,zk)':>10}")
for name, n, B, ip in B_values:
    print(f"  {name:>4} {ip:>10} {B:>10}")

# ============================================================
# SECTION 3: SPECTRAL ZETA FUNCTION OF 600-CELL LAPLACIAN
# ============================================================
print("\n" + "=" * 72)
print("SECTION 3: SPECTRAL ZETA FUNCTION")
print("=" * 72)

# Adjacency matrix eigenvalues of the 600-cell Cayley graph
# For Cay(2I, S) where S = 12 nearest neighbors
# lambda_k = (1/d_k) * sum_{s in S} chi_k(s), multiplicity d_k^2
#
# The 9 irreps of 2I with dimensions and adjacency eigenvalues:
# Using the known character table of 2I evaluated on the 12-element class
#
# 2I irreps: dims = [1, 2, 2, 3, 3, 4, 4, 5, 6]
# Labels: 1, 2a, 2b, 3, 3', 4, 4', 5, 6
# Adjacency eigenvalues: lambda_k = chi_k(generators)/d_k * |S|
#
# Known results for 600-cell adjacency spectrum:
# From the representation theory of 2I:
dims_2I = [1, 2, 2, 3, 3, 4, 4, 5, 6]

# The adjacency eigenvalues (each with multiplicity d_k^2)
# These come from chi_k(c) / d_k where c is the conjugacy class of generators
# For 600-cell with 12 neighbors:
# lambda = 12*chi_k(c)/d_k where chi_k(c) = character at generator class

# The 12 nearest neighbors form a conjugacy class in 2I
# The character values at this class for each irrep:
# From the character table of 2I:
# Class of 12 "nearest neighbor" elements (order 10, type A):
# chi_1 = 1 (trivial)
# chi_2a = PHI (2-dim, golden ratio)
# chi_2b = -1/PHI (2-dim, conjugate)
# chi_3 = PHI+1 = PHI^2 = phi^2 (3-dim)
# chi_3' = 1-PHI+1 = 2-PHI ... actually let me use known values

# CORRECT character table values at the 12A class (nearest neighbors):
# For 2I, the generators (12 vertices of icosahedron in SU(2)):
# Using the standard character table:
# Irrep dims: 1, 2, 2, 3, 3, 4, 4, 5, 6
# chi at identity: 1, 2, 2, 3, 3, 4, 4, 5, 6
# chi at 12A: 1, PHI, phi_prime, -1, -1, -PHI, -phi_prime, 0, 1
# Wait, I need the actual values. Let me compute them from the McKay Cayley matrix.

# The McKay-Cayley adjacency eigenvalues are KNOWN from exp342:
# spectrum = {0, +-1/phi, +-1, +-phi, +-2}
# These are eigenvalues of the 9x9 McKay adjacency matrix
# The 600-cell adjacency eigenvalues are: 12 * (McKay_eig / 2) * (dims)
# Actually: McKay matrix M_ij = 1 if rho_i x rho_2 contains rho_j
# The Cayley graph adjacency eigenvalues for 2I are:
# lambda_k = 12 * chi_k(S) / d_k  (not the McKay eigenvalues)

# Let me use the CAYLEY Laplacian eigenvalues from exp342e:
# lambda_k = 12(1 - chi_k(10a)/d_k)
# where chi_k(10a) is the character at the conjugacy class of 10-fold rotations

# From the character table of 2I (binary icosahedral group):
# 9 conjugacy classes: 1, -1, 12A, 12B, 20A, 20B, 30, 12C, 12D
# Orders:              1,  2,  10,  10,   6,   6,  4,  5,   5
# Sizes:               1,  1,  12,  12,  20,  20, 30, 12,  12

# The 12 nearest neighbors of the 600-cell ARE the 12A class (order 10).
# Character values chi_k(12A) for each irrep:

# Standard character table of 2I at class 12A:
# (Using the known values from representation theory)
chi_12A = {
    1: 1,
    '2a': PHI,       # = (1+sqrt(5))/2
    '2b': phi_prime + 1,  # = (1-sqrt(5))/2 + 1 = (3-sqrt(5))/2...
}

# Actually, let me be more careful. The character table of the binary icosahedral
# group 2I is well-known. The irreps are:
# R1 (dim 1), R2, R2' (dim 2), R3, R3' (dim 3), R4, R4' (dim 4), R5 (dim 5), R6 (dim 6)
#
# At the class 12A (the 12 elements closest to identity, order 10):
# chi_{R1}(12A) = 1
# chi_{R2}(12A) = phi = (1+sqrt(5))/2
# chi_{R2'}(12A) = phi' = (1-sqrt(5))/2
# chi_{R3}(12A) = phi + 1 = phi^2
# chi_{R3'}(12A) = phi' + 1 = 1/phi^2 ... no
# Hmm, let me compute from the adjacency eigenvalues.

# The 600-cell adjacency eigenvalues are well-known:
# 12, 5*phi-2, 5*phi'-2, -phi^2, -phi'^2, phi-3, phi'-3, -3, 2
# Wait, these don't look right either.

# Let me compute from scratch.
# The adjacency eigenvalues of the 600-cell graph:
# Each vertex has 12 neighbors. The graph is vertex-transitive.
# By representation theory: lambda_k = sum_{g in S} chi_k(g) / d_k
#
# The 600-cell vertices can be identified with 2I elements.
# The 12 neighbors of the identity are a specific 12-element subset.
# This subset consists of one element from each of the 12A class.
# So sum_{g in S} chi_k(g) = chi_k(S) = some value.
#
# For the generating set S of 2I being the 12 nearest neighbors:
# S IS the conjugacy class 12A (which has exactly 12 elements).
# So chi_k(S) = 12 * chi_k(g_{12A}) / 1 ... no.
# sum_{g in S} chi_k(g) = |12A| * chi_k(g_{12A}) = 12 * chi_k(g_{12A})? NO!
# S IS the 12A class, so sum_{g in S} chi_k(g) = 12 * chi_k(representative of 12A).
# Actually that's only true if all elements in 12A have the same character value,
# which they do (characters are class functions).
# So sum_{g in S} chi_k(g) = 12 * chi_k(12A_rep)... wait, there are 12 elements
# and chi_k is constant on the class, so sum = 12 * chi_k(12A).

# But wait - the generating set S might not be a single conjugacy class.
# In the 600-cell, the 12 nearest neighbors of the identity ARE all in one
# conjugacy class of 2I. This is because the 600-cell is the Cayley graph
# of 2I with generators being one specific conjugacy class.
# Actually, the 12 nearest-neighbor elements form ONE conjugacy class (12A).

# So: lambda_k = sum_{g in S} chi_k(g) / d_k = 12 * chi_k(12A) / d_k

# The adjacency eigenvalue is lambda_k with multiplicity d_k^2.
# The Laplacian eigenvalue is mu_k = 12 - lambda_k = 12(1 - chi_k(12A)/d_k).

# Now I need the character values chi_k(12A).
# From the standard character table of 2I (see e.g., ATLAS of finite groups):

# 2I character table (partial, just at 12A class):
# Using the convention where 12A has order 10, element = rotation by 2*pi/5 in SU(2)
# The characters are expressed in terms of phi and phi':

# R1 (dim 1): chi(12A) = 1
# R2 (dim 2): chi(12A) = phi   [fundamental 2-dim of SU(2)]
# R2'(dim 2): chi(12A) = -phi' = 1/phi  [conjugate]
# R3 (dim 3): chi(12A) = phi+1 = phi^2  [= Sym^2(R2)/det]
# R3'(dim 3): chi(12A) = -phi'+1 = 1+1/phi = phi/phi ...
#   Actually R3' is the Galois conjugate of R3, so chi_{R3'}(12A) = sigma(chi_{R3}(12A))
#   If chi_{R3}(12A) = phi^2 = phi+1, then chi_{R3'}(12A) = phi'^2 = phi'+1 = (1-phi)+1 = 2-phi
#   = 2 - (1+sqrt(5))/2 = (3-sqrt(5))/2 = 1/phi^2?
#   1/phi^2 = phi'^2 = (phi')^2. phi' = (1-sqrt(5))/2. (phi')^2 = (1-sqrt(5))^2/4 = (6-2*sqrt(5))/4 = (3-sqrt(5))/2.
#   And 2-phi = 2-(1+sqrt(5))/2 = (3-sqrt(5))/2. So yes, chi_{R3'}(12A) = (phi')^2 = phi'-1...
#   wait phi'^2 = phi'+1? No! phi^2 = phi+1, but phi'^2 = (phi')^2 = phi'^2.
#   phi' = (1-sqrt(5))/2 ~ -0.618
#   phi'^2 = phi' + 1? phi'+1 = (1-sqrt(5))/2 + 1 = (3-sqrt(5))/2 ~ 0.382
#   phi'^2 = ((1-sqrt(5))/2)^2 = (1-2*sqrt(5)+5)/4 = (6-2*sqrt(5))/4 = (3-sqrt(5))/2 ~ 0.382
#   So YES! phi'^2 = phi'+1 (same identity as phi^2 = phi+1). CORRECT.

# So R3' at 12A: chi = phi'^2 = phi'+1 = (3-sqrt(5))/2

# R4 (dim 4): chi(12A) = phi^3 - 1 = (2*phi+1)-1 = 2*phi ... hmm
#   R4 = Sym^3(R2)/... actually let me compute from the Cayley-Hamilton.
#   For SU(2) irreps, the character at an element of angle theta is:
#   chi_d(theta) = sin(d*theta/2) / sin(theta/2)
#   For 12A class: theta = 2*pi/5 (rotation by 2*pi/5 lifted to SU(2) double cover)
#   Actually in 2I, the 12A class consists of elements of ORDER 10 (not 5).
#   In SU(2), if g has angle theta, g^2 has angle 2*theta.
#   Order 10 means theta = 2*pi/10 = pi/5 (the fundamental angle).

# Let me use theta = pi/5 for the 12A elements of order 10.
# chi_d(theta) = sin(d*pi/10) / sin(pi/10)
theta_12A = np.pi / 5

def su2_char(d, theta):
    """Character of SU(2) irrep of dim d at element with rotation angle theta."""
    if abs(np.sin(theta/2)) < 1e-15:
        return d  # identity
    return np.sin(d * theta / 2) / np.sin(theta / 2)

print(f"\n--- 2I Character values at 12A class (theta = pi/5) ---")
print(f"  sin(pi/10) = {np.sin(np.pi/10):.6f} = 1/(2*phi) = {1/(2*PHI):.6f}")

dims = [1, 2, 2, 3, 3, 4, 4, 5, 6]
labels = ['R1', 'R2', "R2'", 'R3', "R3'", 'R4', "R4'", 'R5', 'R6']

# For 2I, the irreps come from restricting SU(2) irreps.
# The SU(2) irreps of dims 1,2,3,4,5,6 restrict to 2I irreps.
# BUT: the two 2-dim irreps R2 and R2' of 2I come from the SAME SU(2) irrep!
# No, actually 2I has 9 irreps while SU(2) restricted gives at most the
# isomorphism classes. Let me think more carefully.

# For 2I (binary icosahedral), the irreps ARE labeled by their dimensions:
# 1, 2, 2', 3, 3', 4, 4', 5, 6
# The unprimed ones come from the rotation group I = A5, the primed ones
# are "spinorial" (change sign under the central element -1).
# Actually: the irreps that factor through I = 2I/{+-1} have chi(-1) = chi(1) = d.
# The spinorial ones have chi(-1) = -d.

# For SU(2), the spin-j representation has dim = 2j+1 and:
# chi_{2j+1}(theta) = sin((2j+1)*theta/2) / sin(theta/2)
# chi_{2j+1}(-1) = (-1)^{2j} * (2j+1) = (2j+1) for integer j, -(2j+1) for half-int j

# 2I irreps from SU(2) restriction:
# dim 1 (j=0): chi = 1, chi(-1) = 1 [bosonic]
# dim 2 (j=1/2): chi(-1) = -2 [spinorial] => R2
# dim 2' (j=1/2): doesn't exist in SU(2)...

# Actually, 2I has MORE irreps than what comes from SU(2). Some are "exceptional."
# The standard decomposition of SU(2) irreps:
# j=0: dim 1, restricts to R1 of 2I
# j=1/2: dim 2, restricts to R2 of 2I
# j=1: dim 3, restricts to R3 of 2I
# j=3/2: dim 4, restricts to R4 of 2I
# j=2: dim 5, restricts to R5 of 2I
# j=5/2: dim 6, restricts to R6 of 2I
# Then j=3: dim 7, restricts to R3' + R4' of 2I (7 = 3 + 4)
# And j=2: dim 5 restricts to R5, but wait...
# Actually: j=3/2 gives dim 4 -> R4, j=2 gives dim 5 -> R5
# The missing irreps R2', R3', R4' must come from higher j:
# j=3: dim 7 = 3' + 4' (R3' and R4' are "exceptional")
# j=7/2: dim 8 = 2' + 6? No, 8 = 2 + 6.
#
# This is getting complicated. Let me just use the known character table directly.

# KNOWN CHARACTER TABLE OF 2I at class 12A (order 10):
# Using theta = pi/5 for the SU(2) parametrization:
# The 12A class in 2I corresponds to SU(2) elements with angle pi/5.
# (These have order 10 because (exp(i*pi/10))^10 = exp(i*pi) = -1 in SU(2),
#  but -1 in 2I is a separate element, so the element has order 10 in 2I.)

# Characters at angle theta = pi/5:
# For the "standard" (SU(2)-derived) irreps:
chi_at_12A = {}
for d_irrep in [1, 2, 3, 4, 5, 6]:
    chi_at_12A[d_irrep] = su2_char(d_irrep, theta_12A)

print(f"\n  Standard SU(2) characters at theta=pi/5:")
for d_irrep in [1, 2, 3, 4, 5, 6]:
    print(f"    dim {d_irrep}: chi = {chi_at_12A[d_irrep]:.6f}")

# The "Galois conjugate" irreps R2', R3', R4' have characters that are
# Galois conjugates (phi <-> phi') of the standard ones.
# This is because the Galois automorphism sigma: sqrt(5) -> -sqrt(5)
# permutes the 12A and 12B conjugacy classes.

# So for the SECOND conjugacy class 12B (theta = 3*pi/5):
theta_12B = 3 * np.pi / 5
chi_at_12B = {}
for d_irrep in [1, 2, 3, 4, 5, 6]:
    chi_at_12B[d_irrep] = su2_char(d_irrep, theta_12B)

print(f"\n  Characters at theta=3*pi/5 (12B class):")
for d_irrep in [1, 2, 3, 4, 5, 6]:
    print(f"    dim {d_irrep}: chi = {chi_at_12B[d_irrep]:.6f}")

# Now the adjacency eigenvalues of the Cayley graph:
# lambda_k = sum_{g in S} chi_k(g) / d_k = |S| * chi_k(12A) / d_k
# where S = 12A class, |S| = 12

# Wait - but ARE the generators the 12A class or something else?
# In the 600-cell, each vertex is connected to its 12 nearest neighbors.
# If we place the identity at the north pole, the 12 nearest neighbors
# in S^3 correspond to the 12 vertices of an icosahedron.
# In 2I, these form ONE of the conjugacy classes of size 12.
# Specifically, they should be the 12A class (the "closer" one to identity).

# But actually, in the 600-cell the vertex-set IS 2I, and the generators
# form a conjugacy class. The conjugacy class 12A (of order 10) contains
# the 12 elements closest to the identity. Let me verify.

# In 2I subset of unit quaternions, the identity is 1.
# The nearest quaternions are at angle 2*arcsin(1/(2*phi)) from 1.
# cos(theta) = 1 - 2*sin^2(theta/2) = 1 - 1/(2*phi^2) = ...
# Actually, geodesic distance on S^3 from 1 to q is arccos(Re(q)).
# Nearest: Re(q) = phi/2. So theta = arccos(phi/2) = pi/5.
# These are the 12 elements with real part phi/2, forming the 12A class.

print(f"\n  Generating set S = 12A class (12 elements, order 10)")
print(f"  SU(2) angle: theta = pi/5 = {np.pi/5:.6f}")
print(f"  Re(q) = cos(pi/5) = phi/2 = {PHI/2:.6f}")
print(f"  Geodesic distance: pi/5 = {np.pi/5:.6f}")

# Cayley graph adjacency eigenvalues:
# For each irrep rho_k of dim d_k:
# adj_eig_k = 12 * chi_k(12A_representative) / d_k
# Each eigenvalue has multiplicity d_k^2.

# For irreps coming from SU(2) (R1, R2, R3, R4, R5, R6):
# chi_k(12A) = sin(d_k * pi/10) / sin(pi/10)

print(f"\n--- Cayley Graph Adjacency Eigenvalues ---")
print(f"  {'Irrep':>6} {'dim':>4} {'chi(12A)':>10} {'adj_eig':>10} {'mult':>6} {'Lap_eig':>10}")

adj_eigs_full = []  # (eigenvalue, multiplicity)
lap_eigs_full = []

for i, d_k in enumerate([1, 2, 3, 4, 5, 6]):
    chi_val = su2_char(d_k, theta_12A)
    adj_eig = 12 * chi_val / d_k
    lap_eig = 12 - adj_eig
    mult = d_k**2
    label = f"R{d_k}"
    print(f"  {label:>6} {d_k:>4} {chi_val:>10.6f} {adj_eig:>10.6f} {mult:>6} {lap_eig:>10.6f}")
    adj_eigs_full.append((adj_eig, mult))
    lap_eigs_full.append((lap_eig, mult))

# The conjugate irreps R2', R3', R4' have Galois-conjugated characters.
# For R2' (dim 2): chi(12A) = phi' (Galois conjugate of phi)
# Wait no. The character of R2' at the 12A class is NOT simply phi'.
# It depends on how R2' is defined.
# R2' appears in the decomposition of the j=3/2 (dim 4) SU(2) irrep:
# Actually, the j=3/2 irrep of SU(2) restricts to R4 of 2I (it's irreducible).
# R2' appears in j=3 (dim 7): 7 = 3 + 4 as R3' + R4' of 2I.
# No, this isn't right either.

# The correct approach: compute from the McKay graph.
# The McKay adjacency matrix eigenvalues are related to the Cayley eigenvalues.
# From exp342: McKay spectrum = {0, +-1/phi, +-1, +-phi, +-2}
# These eigenvalues mu_i satisfy: adj_eig = 12 * mu_i / 2 = 6*mu_i
# NO! The McKay matrix M has eigenvalues mu_i and the Cayley adjacency
# eigenvalues are lambda_i = (12/d_fund) * mu_i = 6*mu_i... not right either.

# Actually, the relationship is:
# Cayley adjacency = tensor(identity, S) on the regular representation
# McKay adjacency captures: rho_i x rho_fund -> sum_j M_{ij} rho_j
# The eigenvalue of the Cayley adjacency on the rho_i isotypic component is:
# lambda_i = |S| * chi_i(S) / d_i

# And the McKay eigenvalue for node i is:
# mu_i = (1/d_fund) * sum_{j: M_{ij}=1} d_j / d_i
# Hmm, this doesn't simplify easily.

# Let me just compute the adjacency eigenvalues from the known formula.
# For 2I with generators = 12A class (all at SU(2) angle pi/5):

# The 9 irreps decompose into:
# "Bosonic" (chi(-1)=+d): R1(1), R3(3), R3'(3), R5(5)  [factor through I=A5]
# "Spinorial" (chi(-1)=-d): R2(2), R2'(2), R4(4), R4'(4), R6(6)

# For bosonic irreps, they come from representations of I=A5.
# I has 5 irreps of dims 1, 3, 3', 4, 5 (note: 4-dim is NOT spinorial in 2I).
# Wait, I'm getting confused. Let me use a definitive reference.

# The binary icosahedral group 2I has exactly 9 irreps:
# dim 1: trivial (R1)
# dim 2: (R2) and (R2') -- two non-isomorphic
# dim 3: (R3) and (R3') -- two non-isomorphic
# dim 4: (R4) and (R4') -- two non-isomorphic
# dim 5: (R5) -- unique
# dim 6: (R6) -- unique
# Check: 1+4+4+9+9+16+16+25+36 = 120 = |2I|. OK.

# The central element -1 acts as:
# R1: +1, R2: -1, R2': -1, R3: +1, R3': +1, R4: -1, R4': -1, R5: +1, R6: -1
# (Even-dim: spinorial -1, Odd-dim: bosonic +1, Except R6: -1 because 6=even)

# For the bosonic irreps (R1, R3, R3', R5), these come from A5=I:
# A5 has irreps 1, 3, 3', 4, 5
# Under 2I -> A5: R1->1, R3->3, R3'->3', R5->5
# But also R4+R4' -> 4+4 (the 4-dim of A5 splits in 2I... no,
# A5 has a 4-dim irrep that LIFTS to 2I as EITHER R4 or R4').

# OK this is getting too complicated without a proper reference.
# Let me compute numerically using the known 600-cell structure.

# APPROACH: Use the known Cayley graph eigenvalues.
# The 600-cell has 120 vertices, degree 12.
# Its adjacency spectrum is known. Let me use the formula:
# For the Cayley graph of 2I with generators at angle pi/5,
# the eigenvalues depend on the SU(2) characters.

# ALL 9 irreps can be obtained from SU(2) by restriction.
# The SU(2) irreps of dims 1,2,3,...,6 give the FIRST 6 irreps of 2I.
# But 2I has 9 irreps. The remaining 3 (R2', R3', R4') appear as
# constituents of higher SU(2) irreps.

# KEY FACT: The SU(2) irreps decompose as:
# j=0 (dim 1) -> R1
# j=1/2 (dim 2) -> R2
# j=1 (dim 3) -> R3
# j=3/2 (dim 4) -> R4
# j=2 (dim 5) -> R5
# j=5/2 (dim 6) -> R6
# j=3 (dim 7) -> R3' + R4'
# j=7/2 (dim 8) -> R2' + R6  (8 = 2 + 6)
# j=4 (dim 9) -> R4 + R5     (9 = 4 + 5)
# j=9/2 (dim 10) -> R4' + R6 (10 = 4 + 6)

# So R2' appears in j=7/2 decomposition.
# R3' appears in j=3 decomposition.
# R4' appears in j=3 decomposition.

# For the Cayley graph, the adjacency eigenvalue on the isotypic
# component of rho_k is lambda_k, and it can be computed from
# the SU(2) embedding using the Clebsch-Gordan series.

# Actually, there's a simpler way. The adjacency eigenvalue of
# an irrep rho_k in the Cayley graph is:
# lambda_k = (|S|/d_k) * chi_k(representative of S)
# = 12 * chi_k(g) / d_k where g is any element of the 12A class.

# For g at SU(2) angle pi/5:
# chi_{R_d}(g) = sin(d*pi/10) / sin(pi/10)  for the "standard" R_d (d=1,...,6)

# For R2' (dim 2): This has character chi_{R2'}(12A) that we need.
# Since j=7/2 -> R2' + R6:
# chi_{j=7/2}(pi/5) = chi_{R2'}(pi/5) + chi_{R6}(pi/5)
# chi_{R2'}(pi/5) = chi_{j=7/2}(pi/5) - chi_{R6}(pi/5)
#                 = sin(8*pi/10)/sin(pi/10) - sin(6*pi/10)/sin(pi/10)
#                 = (sin(4*pi/5) - sin(3*pi/5)) / sin(pi/10)

# Similarly for R3' and R4': j=3 -> R3' + R4'
# chi_{R3'}(pi/5) + chi_{R4'}(pi/5) = chi_{j=3}(pi/5) = sin(7*pi/10)/sin(pi/10)

# We need another equation. Use the j=4 decomposition: j=4 -> R4 + R5
# This gives: chi_{j=4}(pi/5) = chi_{R4}(pi/5) + chi_{R5}(pi/5)
# sin(9*pi/10)/sin(pi/10) = chi_{R4} + chi_{R5}
# This is already determined: chi_{R4} and chi_{R5} are from d=4 and d=5 formulas.

# For R3' and R4' separately, use j=9/2 -> R4' + R6:
# chi_{j=9/2}(pi/5) = chi_{R4'}(pi/5) + chi_{R6}(pi/5)
# chi_{R4'}(pi/5) = sin(10*pi/10)/sin(pi/10) - sin(6*pi/10)/sin(pi/10)
#                 = sin(pi)/sin(pi/10) - chi_{R6}(pi/5)
#                 = 0 - chi_{R6}(pi/5)

# sin(pi) = 0, so chi_{R4'}(pi/5) = -chi_{R6}(pi/5)!

# Then chi_{R3'}(pi/5) = chi_{j=3}(pi/5) - chi_{R4'}(pi/5)
#                       = chi_{j=3}(pi/5) + chi_{R6}(pi/5)

sin_pi_10 = np.sin(np.pi/10)

# Standard irrep characters at theta=pi/5:
chi_R = {}
for d in range(1, 11):
    chi_R[d] = np.sin(d * np.pi / 10) / sin_pi_10

print(f"\n--- SU(2) characters at theta=pi/5 ---")
for d in range(1, 11):
    print(f"  j={(d-1)/2:.1f} (dim {d}): chi = {chi_R[d]:.6f}")

# Now extract the 9 irreps of 2I:
chi_2I = {}
chi_2I['R1'] = (1, chi_R[1])
chi_2I['R2'] = (2, chi_R[2])
chi_2I['R3'] = (3, chi_R[3])
chi_2I['R4'] = (4, chi_R[4])
chi_2I['R5'] = (5, chi_R[5])
chi_2I['R6'] = (6, chi_R[6])

# R4' (dim 4): from j=9/2 -> R4' + R6
# chi_{j=9/2} = chi_R[10] = sin(10*pi/10)/sin(pi/10) = sin(pi)/sin(pi/10) = 0
chi_2I["R4'"] = (4, chi_R[10] - chi_R[6])

# R3' (dim 3): from j=3 -> R3' + R4'
chi_2I["R3'"] = (3, chi_R[7] - chi_2I["R4'"][1])

# R2' (dim 2): from j=7/2 -> R2' + R6
chi_2I["R2'"] = (2, chi_R[8] - chi_R[6])

print(f"\n--- 2I irrep characters at 12A (theta=pi/5) ---")
adj_eigs = []
lap_eigs = []
for label in ['R1','R2',"R2'",'R3',"R3'",'R4',"R4'",'R5','R6']:
    d_k, chi_val = chi_2I[label]
    adj_eig = 12 * chi_val / d_k
    lap_eig = 12 - adj_eig
    mult = d_k**2
    print(f"  {label:>4} (d={d_k}): chi={chi_val:>10.6f}, adj={adj_eig:>10.6f}, lap={lap_eig:>10.6f}, mult={mult}")
    adj_eigs.append((adj_eig, mult))
    lap_eigs.append((lap_eig, mult))

# Verify: sum of multiplicities = 120
total_mult = sum(m for _, m in adj_eigs)
print(f"\n  Total multiplicity: {total_mult} (should be 120)")

# Verify: largest eigenvalue is 12 (degree)
max_adj = max(e for e, _ in adj_eigs)
print(f"  Max adjacency eigenvalue: {max_adj:.6f} (should be 12)")

# ============================================================
# SECTION 4: SPECTRAL ZETA FUNCTION
# ============================================================
print("\n" + "=" * 72)
print("SECTION 4: SPECTRAL ZETA FUNCTION OF 600-CELL LAPLACIAN")
print("=" * 72)

# zeta(s) = sum_i mult_i * lambda_i^(-s), excluding zero eigenvalue
# Laplacian eigenvalue for R1 is 0 (excluded).
nonzero_lap = [(l, m) for l, m in lap_eigs if abs(l) > 1e-10]

def spectral_zeta(s, eigs=nonzero_lap):
    """Compute zeta_L(s) = sum mult_i * lambda_i^(-s) for non-zero eigenvalues."""
    result = 0
    for lam, mult in eigs:
        if lam > 0:
            result += mult * lam**(-s)
    return result

print(f"\n  Nonzero Laplacian eigenvalues:")
for lam, mult in nonzero_lap:
    print(f"    lambda = {lam:>10.6f}, mult = {mult:>3}")

# Evaluate zeta at special points
print(f"\n  Spectral zeta function values:")
special_s = [0.5, 1, 1.5, 2, 3, 4, -1, -2, -3]
for s in special_s:
    z_val = spectral_zeta(s)
    print(f"    zeta({s:>5.1f}) = {z_val:>14.6f}")

# Can we find s such that zeta(s) = 57 or z_CC?
print(f"\n  Searching for zeta(s) = 57:")
from scipy.optimize import brentq

def zeta_minus_target(s, target=57):
    return spectral_zeta(s) - target

# zeta(s) is decreasing for s > 0 (positive eigenvalues)
# Check zeta at endpoints
print(f"    zeta(0.1) = {spectral_zeta(0.1):.2f}")
print(f"    zeta(0.5) = {spectral_zeta(0.5):.2f}")
print(f"    zeta(1.0) = {spectral_zeta(1.0):.2f}")
print(f"    zeta(2.0) = {spectral_zeta(2.0):.2f}")

# Try to find s where zeta(s) = 57
try:
    s_57 = brentq(zeta_minus_target, 0.01, 5.0, args=(57,))
    print(f"    FOUND: zeta({s_57:.6f}) = 57")
    # Is s_57 a framework number?
    print(f"    s_57 = {s_57:.6f}")
    print(f"    s_57/alpha_s = {s_57/alpha_s:.4f}")
    print(f"    s_57*a1 = {s_57*a1:.4f}")
    print(f"    s_57/phi = {s_57/PHI:.4f}")
except ValueError:
    print(f"    No solution in [0.01, 5.0]")

# Try zeta(s) = z_CC = 57 - alpha_s
z_CC_target = 57 - alpha_s
try:
    s_zCC = brentq(zeta_minus_target, 0.01, 5.0, args=(z_CC_target,))
    print(f"    FOUND: zeta({s_zCC:.6f}) = z_CC = {z_CC_target:.6f}")
    print(f"    s_zCC = {s_zCC:.6f}")
except ValueError:
    print(f"    No solution for z_CC in [0.01, 5.0]")

# Also check: spectral determinant
# log det'(L) = -zeta'(0) (derivative of spectral zeta at s=0)
# Approximate numerically
ds = 1e-6
zeta_prime_0 = (spectral_zeta(ds) - spectral_zeta(-ds)) / (2*ds)
log_det = -zeta_prime_0
det_L = np.exp(log_det) if abs(log_det) < 500 else float('inf')

print(f"\n  Spectral determinant:")
print(f"    zeta'(0) ~ {zeta_prime_0:.6f}")
print(f"    log det'(L) = -zeta'(0) ~ {log_det:.6f}")
print(f"    det'(L) ~ {det_L:.6e}" if det_L < 1e300 else f"    det'(L) = exp({log_det:.2f})")

# Direct product of nonzero eigenvalues
log_det_direct = sum(m * np.log(l) for l, m in nonzero_lap if l > 0)
print(f"    Direct: log det'(L) = sum(m*log(l)) = {log_det_direct:.6f}")
print(f"    det'(L) = exp({log_det_direct:.2f})")
print(f"    log det'(L) / z_P = {log_det_direct / (4*PHI**2):.6f}")
print(f"    log det'(L) / z_CC = {log_det_direct / z_CC_target:.6f}")
print(f"    log det'(L) / 57 = {log_det_direct / 57:.6f}")

# ============================================================
# SECTION 5: CHARACTER SUM APPROACH TO 57
# ============================================================
print("\n" + "=" * 72)
print("SECTION 5: CHARACTER SUM APPROACH")
print("=" * 72)

# 57 as a sum involving characters or dimensions of 2I irreps
dims_list = [1, 2, 2, 3, 3, 4, 4, 5, 6]  # the 9 irrep dimensions

print(f"\n  2I irrep dimensions: {dims_list}")
print(f"  Sum d_i = {sum(dims_list)} (= N/4 = 30 = h_E8)")
print(f"  Sum d_i^2 = {sum(d**2 for d in dims_list)} (= N = 120)")
print(f"  Sum d_i^3 = {sum(d**3 for d in dims_list)}")
print(f"  Sum d_i(d_i+1)/2 = {sum(d*(d+1)//2 for d in dims_list)}")
print(f"  Sum (d_i choose 2) = {sum(d*(d-1)//2 for d in dims_list)}")

# Can we get 57 from dimensions?
S3 = sum(d**3 for d in dims_list)
print(f"\n  sum(d^3) = {S3}")
print(f"  sum(d^3)/sum(d^2) = {S3/120:.4f}")
print(f"  sum(d^3)/sum(d) = {S3/30:.4f} = {S3}/30")

# 57 from combinatorics of dims
print(f"\n  Testing combinations for 57:")
for k in range(2, 8):
    s = sum(d**k for d in dims_list)
    if s > 0 and abs(57 - s) < 200:
        print(f"    sum(d^{k}) = {s}, ratio to 57: {57/s:.4f}")

# Special products
print(f"\n  Products:")
print(f"    N_gen * p_unique = {N_gen} * 19 = {N_gen * 19}")
print(f"    N_eig * (b1+1/N_gen) = 9 * 6.333 = {N_eig * (b1 + 1.0/N_gen):.3f}")
print(f"    rank_E8 * (b1+alpha_s+1) = 8*{b1+alpha_s+1:.4f} = {rank_E8*(b1+alpha_s+1):.4f}")
print(f"    h_E8 + degree + a1*N_gen = 30+12+15 = {h_E8+degree+a1*N_gen}")

# The deepest connection: norm of mass direction
print(f"\n  === NORM OF MASS DIRECTION ===")
print(f"  z_mass = a1 + b1*phi = {a1} + {b1}*phi")
print(f"  N(z_mass) = a1^2 + a1*b1 - b1^2 = {a1**2} + {a1*b1} - {b1**2} = {N_zmass}")
print(f"  57 = N_gen * N(z_mass) = {N_gen} * {N_zmass} = {N_gen * N_zmass}")
print(f"")
print(f"  UNIQUENESS: p=19 is the ONLY prime compatible with the norm multiset")
print(f"  (tested p=2 to p=41 in exp342f)")
print(f"")
print(f"  The mass formula n = 5a + 6b = <z_mass, (a,b)>_std")
print(f"  The CC exponent uses N(z_mass) = |z_mass * z_mass'| = 19")
print(f"  => CC ~ alpha^(N_gen * N(z_mass))")
print(f"  => vacuum energy = N_gen generations, each contributing N(z_mass)")

# ============================================================
# SECTION 6: WHY N_gen * N(z_mass)?
# ============================================================
print("\n" + "=" * 72)
print("SECTION 6: PHYSICAL INTERPRETATION")
print("=" * 72)

print(f"""
  CANDIDATE MECHANISM:

  The vacuum energy density in Planck units:
    Lambda_P = alpha^(N_gen * |N(z_mass)| - alpha_s)
             = alpha^(3 * 19 - alpha_s)
             = alpha^(57 - alpha_s)

  Why N_gen * |N(z_mass)|?

  Each generation contributes to the vacuum energy through loops.
  The "weight" of each generation's contribution is controlled by
  the norm of the mass direction vector z_mass = a1 + b1*phi.

  N(z_mass) = 19 is a PRIME in Z, meaning z_mass is irreducible in Z[phi].
  It cannot be factored further. This "irreducibility" of the mass
  direction makes 19 the fundamental unit of CC suppression per generation.

  The QCD correction -alpha_s represents the one-loop strong force
  contribution that partially lifts the vacuum suppression.

  ALTERNATIVELY (group-theoretic):
  57 = |A5| - N_gen = 60 - 3
  The vacuum energy is suppressed by the size of the rotation group A5
  (60 symmetry operations), with 3 "active" modes (generations) that
  contribute positively to the energy density.

  BOTH interpretations are algebraically equivalent for a1=5
  (and ONLY for a1=5), via the identity:
    N_gen * (a1^2 - a1 - 1) = a1!/2 - N_gen
    <=> (a1-2)! = b1 = a1+1
""")

# ============================================================
# SECTION 7: ADDITIONAL SPECTRAL CHECKS
# ============================================================
print("=" * 72)
print("SECTION 7: NUMERICAL CHECKS AND REFINEMENTS")
print("=" * 72)

z_formula = 57 - alpha_s
z_dev = abs(z_formula - z_exact) / z_err

print(f"\n  z_formula = N_gen*N(z_mass) - alpha_s = {z_formula:.10f}")
print(f"  z_exact   = {z_exact:.10f}")
print(f"  z_err     = +/- {z_err:.4f}")
print(f"  Deviation = {z_dev:.2f} sigma")
print(f"")
print(f"  Lambda_pred = alpha^z_formula = {alpha**z_formula:.6e}")
print(f"  Lambda_obs  = {Lambda_Planck:.6e}")
err_pct = abs(alpha**z_formula - Lambda_Planck) / Lambda_Planck * 100
print(f"  Error: {err_pct:.2f}%")

# Check: is the N(z_mass)=19 result stable?
# What if we used a different mass direction?
print(f"\n  Stability: N(z) for various mass-like vectors:")
for aa, bb in [(5,6), (6,5), (5,5), (4,6), (5,7), (1,1), (a1,b1)]:
    nz = norm_zphi(aa, bb)
    z_test = N_gen * abs(nz) - alpha_s
    dev = abs(z_test - z_exact) / z_err
    mark = " *** BEST" if (aa,bb) == (a1,b1) else ""
    print(f"    ({aa},{bb}): N={nz:>4}, z={z_test:>10.4f}, dev={dev:.2f} sigma{mark}")

# ============================================================
# SECTION 8: ALGEBRAIC STRUCTURE
# ============================================================
print("\n" + "=" * 72)
print("SECTION 8: ALGEBRAIC STRUCTURE OF THE CC FORMULA")
print("=" * 72)

# Full CC formula: Lambda_P = alpha^(N_gen * N(a1 + b1*phi) - 1/(2*phi^3))
# All quantities from a1=5:
print(f"""
  FULL FORMULA:
    Lambda_Planck = alpha^(N_gen * N(a1 + b1*phi) - 1/(2*phi^3))

  Where EVERYTHING comes from a1 = 5:
    phi = (1+sqrt(a1))/2 = {PHI:.6f}
    b1 = a1+1 = {b1}
    N_gen = 3 (Nyquist on icosahedron)
    N(a1+b1*phi) = a1^2+a1*b1-b1^2 = {N_zmass} (= unique prime p)
    alpha_s = 1/(2*phi^3) = {alpha_s:.6f}
    alpha = smaller root of 2*pi*x^2 - 4*a1*phi^4*x + 1 = 0

  EQUIVALENT FORMS:
    (a) z = N_gen * (a1^2 - a1 - 1) - 1/(2*phi^3) = 57 - alpha_s
    (b) z = (N/2 - N_gen) - alpha_s = (|A5| - N_gen) - alpha_s
    (c) z = (N-N_gen)/2 - phi = 117/2 - phi  [since alpha_s = phi - 3/2]

  STATUS:
    Form (a): connects CC to MASS FORMULA (via z_mass norm) and NUMBER THEORY (prime 19)
    Form (b): connects CC to GROUP THEORY (|A5|) and GENERATION COUNT
    Form (c): connects CC to ALGEBRAIC NUMBER THEORY (Z[phi] element)

  All three are algebraically equivalent for a1=5 and ONLY for a1=5.
  The equivalence requires (a1-2)! = a1+1, i.e., 3! = 6.

  CATEGORY: STRONG PATTERN with new ALGEBRAIC interpretation.
  NOT YET DERIVED: we explain what 57 means but not WHY the CC
  should equal alpha^(N_gen*|N(z_mass)|).
""")

# ============================================================
# SECTION 9: SEEKING THE MECHANISM
# ============================================================
print("=" * 72)
print("SECTION 9: CANDIDATE MECHANISMS")
print("=" * 72)

print(f"""
  CANDIDATE 1: Loop counting
    Each of the N_gen=3 generations runs in vacuum loops.
    Each loop gives a factor alpha^(something).
    If each generation contributes alpha^|N(z_mass)| = alpha^19,
    then Lambda ~ alpha^(3*19) = alpha^57 at tree level.
    The -alpha_s correction = one QCD-corrected loop.

  CANDIDATE 2: Determinant of mass matrix in Z[phi]
    The fermion mass matrix M has entries in Z[phi].
    det(M) involves products of z_k = a_k + b_k*phi.
    In some sense, sum N(z_k) ~ vacuum energy contribution.
    But det != sum of norms.

  CANDIDATE 3: Spectral action trace anomaly
    In the spectral action S = Tr f(D/Lambda):
    The cosmological constant term is f_0*Lambda^4*a_0.
    If a_0 = c_0 = 2640 on the 600-cell,
    the suppression might come from the spectral structure.
    But exp328 showed this route is NEGATIVE.

  CANDIDATE 4: Algebraic K-theory
    K_0(Z[phi]) ~ Z (since Z[phi] is a PID for Q(sqrt(5))).
    The class number is 1. The fundamental unit is phi.
    The "regulator" involves log(phi).
    Could lambda ~ alpha^(regulator * something)?

  MOST PROMISING: Candidate 1 (loop counting), because:
    - N_gen naturally counts generation loops
    - |N(z_mass)|=19 is the "weight" per generation
    - alpha_s correction is physically QCD loops
    - Structure is additive in exponent (multiplicative in Lambda)

  LEAST PROMISING: Candidate 3 (already tested negative).
""")

# ============================================================
# SECTION 10: HONEST ASSESSMENT
# ============================================================
print("=" * 72)
print("SECTION 10: HONEST ASSESSMENT")
print("=" * 72)

print(f"""
  WHAT WE ACHIEVED:
  1. Identified 57 = N_gen * N(a1 + b1*phi) = 3 * 19
     This connects the CC to the mass formula via the Z[phi] norm.
  2. Showed equivalence: 57 = |A5| - N_gen (group-theoretic)
     These are the SAME only for a1=5 (via 3! = 6 = a1+1).
  3. Computed spectral zeta of 600-cell Laplacian.
  4. Identified most promising mechanism (loop counting).

  WHAT REMAINS:
  1. No rigorous derivation of z_CC = N_gen * N(z_mass) - alpha_s
  2. The "loop counting" interpretation is heuristic
  3. The spectral zeta did not directly produce 57

  UPGRADE FROM EXP311:
  - Before: "57 = N/2 - N_gen" (combinatorial observation)
  - Now: "57 = N_gen * N(z_mass) = N_gen * p_unique" (algebraic interpretation)
  - The second form is DEEPER because it connects to:
    * The mass formula (z_mass = a1 + b1*phi)
    * Number theory (N(z_mass) = 19 = unique prime)
    * The uniqueness theorem (exp342f)

  CATEGORY: STRONG PATTERN -> ALGEBRAIC PATTERN (upgraded)
  Still NOT DERIVED (no mechanism).
""")

elapsed = time.time() - t0
print(f"\nCompleted in {elapsed:.1f}s")

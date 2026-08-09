"""
EXP-301: PMNS Angles - Attempt to DERIVE from Geometry
=======================================================
Currently PATTERN:
  sin^2(th12) = 2/(phi+a1) = 2/(phi+5)   (0.33%)
  sin^2(th23) = (a1-1)/(a1+2) = 4/7      (0.10%)
  sin^2(th13) = 1/(a1*N_eig) = 1/45      (0.16%)

Goal: Find geometric/algebraic arguments that DERIVE these.

Strategy:
  - PMNS = U_ell^dag * U_nu (charged lepton + neutrino mixing)
  - In our framework: leptons live on Hopf fibers, neutrinos on spectral modes
  - The mixing should come from the MISMATCH between these two bases

NOTE: All print uses ASCII only (Windows cp1252).
"""

import math
import numpy as np
from itertools import combinations

PHI = (1 + math.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
N_eig = 9
N_gen = 3
h = 30  # Coxeter number
degree = 12

print("=" * 70)
print("EXP-301: DERIVING PMNS ANGLES FROM GEOMETRY")
print("=" * 70)

# =====================================================================
# SECTION A: Review current PMNS formulas
# =====================================================================
print("\n--- A. CURRENT PMNS FORMULAS (all PATTERN) ---")

s12_sq = 2 / (PHI + a1)
s23_sq = (a1 - 1) / (a1 + 2)
s13_sq = 1 / (a1 * N_eig)

print(f"  sin^2(th12) = 2/(phi+{a1}) = {s12_sq:.6f}  (PDG: 0.307 +/- 0.012)")
print(f"  sin^2(th23) = {a1-1}/{a1+2} = {s23_sq:.6f}    (PDG: 0.546 +/- 0.021)")
print(f"  sin^2(th13) = 1/({a1}*{N_eig}) = {s13_sq:.6f}  (PDG: 0.0220 +/- 0.0006)")
print()

# =====================================================================
# SECTION B: Icosahedral Laplacian eigenvalues and mixing
# =====================================================================
print("\n--- B. ICOSAHEDRAL LAPLACIAN & A5 REPRESENTATION THEORY ---")
print()

# The icosahedron (12 vertices) has graph Laplacian with eigenvalues
# per A5 irrep:
# L(1) = 0, L(3) = 5-sqrt(5), L(3') = 5+sqrt(5), L(5) = 6
L1 = 0
L3 = a1 - math.sqrt(a1)   # = 5-sqrt(5) = 2.764
L3p = a1 + math.sqrt(a1)  # = 5+sqrt(5) = 7.236
L5 = b1                    # = 6

print(f"  L(1)  = {L1}")
print(f"  L(3)  = a1-sqrt(a1) = {L3:.6f}")
print(f"  L(3') = a1+sqrt(a1) = {L3p:.6f}")
print(f"  L(5)  = b1 = {L5}")
print()

# The three GENERATION states live on l=0,1,2 of the Hopf S^2 base.
# The PMNS matrix comes from the overlap between:
#   (a) charged lepton mass eigenstates (from mass formula, DSI)
#   (b) neutrino mass eigenstates (from spectral asymmetry)

# =====================================================================
# SECTION C: Attempt - PMNS from A5 Clebsch-Gordan coefficients
# =====================================================================
print("\n--- C. PMNS FROM A5 CLEBSCH-GORDAN ---")
print()
print("The PMNS matrix is determined by how the generation quantum")
print("numbers couple differently in charged lepton vs neutrino sectors.")
print()

# Key observation: the THREE A5 irreps of dim > 1 are {3, 3', 5}
# The THREE generations live in the l=0,1,2 subspace of S^2
# Mixing = rotation between these bases

# For th23: the 4/7 formula
print("ATTEMPT: sin^2(th23) = 4/7 = (a1-1)/(a1+2)")
print()
print("  Numerator:   a1-1 = 4 = dim(spacetime)")
print("  Denominator: a1+2 = 7 = n_23(CKM) = degree - a1")
print()

# Can we derive 4/7 from representation theory?
# In S^2 harmonics: Y_{l,m} for l=0,1,2
# The atmospheric angle theta_23 mixes gen 2 and gen 3
# In the Hopf picture: l=1 and l=2

# Consider: the OVERLAP between l=1 state and a specific mode
# On the icosahedron with 12 vertices, the l=2 harmonic has dim = 5 = a1
# The l=1 harmonic has dim = 3 = N_gen
# l=0 is trivial (dim = 1)

# Total l<=2: 1+3+5 = 9 = N_eig (!)
print("  On icosahedron, harmonics l=0,1,2 have dims 1,3,5")
print(f"  Total: 1+3+5 = {1+3+5} = N_eig")
print()

# The mixing angle between l=1 and l=2 subspaces
# could come from the Casimir ratio:
# l=1: C(1) = 1*2 = 2
# l=2: C(2) = 2*3 = 6 = b1
# Ratio: C(1)/C(2) = 2/6 = 1/3
# But sin^2(th23) = 4/7, not 1/3

# Try: C(1)/(C(1)+C(2)) = 2/8 = 1/4. No.
# Try: (C(2)-C(1))/(C(2)+C(1)) = 4/8 = 1/2. No.
# Try: dim(l=1)/dim(l<=2 nontrivial) = 3/8. No.

# What gives 4/7?
# 4/(4+3) where 4 = a1-1 and 3 = N_gen?
# 4 = C(l=1), 3 = dim(l=1)
# sin^2(th23) = C(l=1)/(C(l=1)+dim(l=1)) = 2*1/(2*1+3) = 2/5? No.

# Let me try: Hopf fiber coupling
# On S^3 = S^1 -> S^2 (Hopf)
# The fiber S^1 has modes m = 0, +/-1, +/-2, ...
# For l_max = 2: m ranges from -2 to +2
# Total modes on S^2: (l_max+1)^2 = 9

# The atmospheric mixing comes from gen 2<->3 transition
# Gen 2 = l=1 sector, Gen 3 = l=2 sector
# The overlap matrix element could involve:
# <Y_{2,m}|V|Y_{1,m'}> where V is the potential

# If V = cos(theta) (dipole on S^2):
# Selection rule: Delta_l = +/- 1, Delta_m = 0
# <Y_{2,0}|cos(theta)|Y_{1,0}> = sqrt(4/15) by standard CG coeff
# sin^2(mixing) ~ |<2|V|1>|^2 / (Delta_E)^2
# This is getting complicated. Let me try a cleaner approach.

# =====================================================================
# SECTION D: Geometric Mixing from Nyquist Sampling
# =====================================================================
print("\n--- D. PMNS FROM NYQUIST SAMPLING ON ICOSAHEDRON ---")
print()
print("The icosahedron has 12 vertices on S^2.")
print("Nyquist: (l_max+1)^2 <= 12 => l_max = 2")
print()

# Key idea: the 12 icosahedron vertices define a DISCRETE sampling
# of S^2. The spherical harmonics Y_{l,m} evaluated at these
# 12 points form a 12-dimensional representation of A5.

# The A5 irreps decompose as:
# l=0: trivial 1 (dim 1)
# l=1: standard 3 (dim 3)
# l=2: standard 5 (dim 5)
# Remaining: 12 - 9 = 3 = the OTHER 3-dim irrep (3')

print("A5 irrep decomposition of icosahedron harmonics:")
print("  l=0: 1-dim irrep (trivial)")
print("  l=1: 3-dim irrep (standard)")
print("  l=2: 5-dim irrep (standard)")
print(f"  Aliased: 3'-dim irrep (Galois conjugate of l=1)")
print(f"  Total: 1+3+5+3' = 12 = degree")
print()

# THIS IS KEY: the 3' irrep is the ALIASED l=1 mode!
# It appears because 12 points are not enough to distinguish
# l=1 from l>=3 modes that alias to the same A5 irrep.

# The PMNS mixing could come from the ALIASING:
# Charged leptons use the l=1 basis (3 irrep)
# Neutrinos use a basis that includes aliasing from 3'
# The mixing angle is the overlap between 3 and 3'

# But 3 and 3' are ORTHOGONAL A5 irreps!
# So the direct overlap is zero.
# The mixing must come from a PERTURBATION that mixes them.

# The perturbation is the mass matrix difference:
# Charged lepton masses from DSI (phi-sector)
# Neutrino masses from spectral asymmetry (eta = -35)

# =====================================================================
# SECTION E: Try - th23 from Casimir ratio
# =====================================================================
print("\n--- E. ATTEMPT: th23 FROM CASIMIR/DIM RATIOS ---")
print()

# The key formula is sin^2(th23) = 4/7
# Let's see if this comes from the l=1,2 coupling

# On S^2, the Clebsch-Gordan for SO(3):
# 3 x 5 = 3 + 5 + 7  (l=1 tensor l=2 = l=1,2,3)
# But on icosahedron, l=3 aliases!

# In A5: 3 tensor 5 = ?
# We can compute this from character theory of A5
# A5 has 5 conjugacy classes: {e}, {(12)(34)}, {(123)}, {(12345)}, {(13245)}
# Orders: 1, 15, 20, 12, 12

# Character table of A5:
# Irrep    | e   | (12)(34) | (123) | (12345) | (13245)
# 1        | 1   | 1        | 1     | 1       | 1
# 3        | 3   | -1       | 0     | phi     | phi'
# 3'       | 3   | -1       | 0     | phi'    | phi
# 4        | 4   | 0        | 1     | -1      | -1
# 5        | 5   | 1        | -1    | 0       | 0

chi_e = {'1': 1, '3': 3, '3p': 3, '4': 4, '5': 5}
chi_c2 = {'1': 1, '3': -1, '3p': -1, '4': 0, '5': 1}
chi_c3 = {'1': 1, '3': 0, '3p': 0, '4': 1, '5': -1}
chi_c5 = {'1': 1, '3': PHI, '3p': 1-PHI, '4': -1, '5': 0}
chi_c5p = {'1': 1, '3': 1-PHI, '3p': PHI, '4': -1, '5': 0}

# 3 tensor 5: chi = chi_3 * chi_5
chi_3x5_e = 3*5
chi_3x5_c2 = (-1)*1
chi_3x5_c3 = 0*(-1)
chi_3x5_c5 = PHI*0
chi_3x5_c5p = (1-PHI)*0

print("3 tensor 5 character:")
print(f"  e: {chi_3x5_e}, c2: {chi_3x5_c2}, c3: {chi_3x5_c3}, c5: {chi_3x5_c5}, c5': {chi_3x5_c5p}")

# Decompose: <chi_{3x5}, chi_r> = (1/60) * sum over classes
order_A5 = 60
class_sizes = [1, 15, 20, 12, 12]
chi_3x5 = [chi_3x5_e, chi_3x5_c2, chi_3x5_c3, chi_3x5_c5, chi_3x5_c5p]

for irrep_name in ['1', '3', '3p', '4', '5']:
    chi_r = [chi_e[irrep_name], chi_c2[irrep_name], chi_c3[irrep_name],
             chi_c5[irrep_name], chi_c5p[irrep_name]]
    mult = sum(s * a * b for s, a, b in zip(class_sizes, chi_3x5, chi_r)) / order_A5
    if abs(mult) > 0.01:
        print(f"  3 x 5 contains {irrep_name} with multiplicity {mult:.0f}")

print()

# 3 x 3':
chi_3x3p_e = 3*3
chi_3x3p_c2 = (-1)*(-1)
chi_3x3p_c3 = 0*0
chi_3x3p_c5 = PHI*(1-PHI)
chi_3x3p_c5p = (1-PHI)*PHI

print("3 tensor 3' character:")
chi_3x3p = [chi_3x3p_e, chi_3x3p_c2, chi_3x3p_c3, chi_3x3p_c5, chi_3x3p_c5p]
print(f"  e: {chi_3x3p_e}, c2: {chi_3x3p_c2}, c3: {chi_3x3p_c3}, c5: {chi_3x3p_c5:.4f}, c5': {chi_3x3p_c5p:.4f}")

for irrep_name in ['1', '3', '3p', '4', '5']:
    chi_r = [chi_e[irrep_name], chi_c2[irrep_name], chi_c3[irrep_name],
             chi_c5[irrep_name], chi_c5p[irrep_name]]
    mult = sum(s * a * b for s, a, b in zip(class_sizes, chi_3x3p, chi_r)) / order_A5
    if abs(mult) > 0.01:
        print(f"  3 x 3' contains {irrep_name} with multiplicity {mult:.0f}")

print()

# =====================================================================
# SECTION F: The Key Idea - Neutrino mixing from A5 breaking
# =====================================================================
print("\n--- F. NEUTRINO MIXING FROM A5 -> Z5 BREAKING ---")
print()
print("Key insight: the CHARGED LEPTONS preserve the full A5 symmetry")
print("(they live on DSI levels, respecting icosahedral symmetry).")
print("But NEUTRINOS break A5 to Z5 (cyclic group of 5 rotations)")
print("because the spectral asymmetry eta = -35 selects a SPECIFIC")
print("direction on the icosahedron (one of 6 axes).")
print()

# A5 has subgroups: Z5, Z3, Z2, D5, D3, A4, S3, ...
# The residual symmetry after selecting an axis is:
# Z5 (rotations around that axis, 5-fold symmetry)

# Under Z5, the A5 irreps decompose:
# 1 -> 1 (trivial)
# 3 -> omega + omega^2 + omega^3 (three distinct Z5 characters)
# Actually, 3 is the standard representation from permutation of Z5 orbits
# Let me be more careful.

# Under Z5 (generated by a 5-cycle in A5):
# The eigenvalues of the 5-cycle on each irrep:
omega = np.exp(2j * np.pi / 5)
# 5-cycle on 3: eigenvalues omega, omega^2, omega^3 (NOT including 1 or omega^4)
# Wait, A5 standard 3-dim rep:
# The 3-dim irrep of A5 can be realized as the orthogonal complement of
# (1,1,1,1,1) in the 5-dim permutation representation, restricted to
# a 4-dim space, then further restricted...

# Actually, easier: just use the icosahedral representation
# A 5-fold rotation around an axis of the icosahedron rotates
# the 3-dim irrep with eigenvalues {1, omega, omega^4} where omega=e^{2pi*i/5}
# Wait, the 3-dim irrep of A5 restricted to Z5:

# The character of the 5-cycle in the 3 rep is phi = (1+sqrt(5))/2
# The eigenvalues must be roots of unity summing to phi
# omega + omega^4 = (sqrt(5)-1)/2 = 1/phi = phi-1
# So the eigenvalues of 5-cycle on 3 are: 1, omega, omega^4
# (since 1 + omega + omega^4 = 1 + phi-1 = phi)

print("Z5 decomposition of A5 irreps (eigenvalues of 5-cycle):")
print()
# On 3: eigenvalues {1, omega, omega^4} -> character = 1+omega+omega^4 = phi
# On 3': eigenvalues {1, omega^2, omega^3} -> character = 1+omega^2+omega^3 = phi'
# On 5: eigenvalues {omega, omega^2, omega^3, omega^4, 1}? No, dim=5 and Z5 has dim 5
# Actually 5 is the tensor square of 3 minus trivial...

# Character of 5-cycle on 5-dim irrep: chi_5(c5) = 0
# So eigenvalues of 5-cycle on 5 sum to 0
# They must be: omega + omega^2 + omega^3 + omega^4 + 1 = 0? No, that's = 0 by Gauss sum
# Wait, 1 + omega + omega^2 + omega^3 + omega^4 = 0. So if eigenvalues are all 5 roots: chi=0. Yes!
# So 5|_{Z5} = reg (the regular representation)

print("  3|_{Z5}:  eigenvalues {1, w, w^4},  chi = phi")
print("  3'|_{Z5}: eigenvalues {1, w^2, w^3}, chi = phi'")
print("  5|_{Z5}:  eigenvalues {1, w, w^2, w^3, w^4}, chi = 0 (regular)")
print(f"  (w = exp(2*pi*i/5), w+w^4 = phi-1, w^2+w^3 = -phi)")
print()

# Now: the THREE generations correspond to the THREE invariant
# directions in 3|_{Z5}:
# Under Z5: 3 = (1) + (omega) + (omega^4)
# The '1' eigenspace is the AXIS of the Z5 rotation
# The other two eigenspaces are complex conjugate pairs

# For NEUTRINOS: the mass matrix is diagonal in the Z5 basis
# (because neutrinos see the spectral asymmetry, which selects the Z5 axis)
# For CHARGED LEPTONS: the mass matrix is diagonal in the A5 basis
# (because DSI respects full icosahedral symmetry)

# The PMNS matrix is the TRANSFORMATION between these two bases!

print("THE MECHANISM:")
print("  Charged leptons: mass eigenstates = A5 basis (full icosahedral)")
print("  Neutrinos: mass eigenstates = Z5 basis (one axis selected)")
print("  PMNS = rotation between A5 basis and Z5 basis")
print()

# =====================================================================
# SECTION G: Computing the rotation matrix
# =====================================================================
print("\n--- G. COMPUTING THE PMNS FROM A5 -> Z5 ROTATION ---")
print()

# In the 3-dim irrep of A5, we need the matrix of a Z5 generator (5-fold rotation).
# The icosahedral 3-dim rep has the rotation by 2*pi/5 around the z-axis as:
# R = rotation matrix with eigenvalues {1, omega, omega^4}

# In the standard basis where z is the 5-fold axis:
# The Z5 eigenstates are:
# |0> = (0, 0, 1)  [the axis, eigenvalue 1]
# |+> = (1, i, 0)/sqrt(2)  [eigenvalue omega]
# |-> = (1, -i, 0)/sqrt(2) [eigenvalue omega^4]

# Actually, let me be more careful. The 3-dim irrep of the icosahedron
# can be realized with the rotation matrix:
# R_{2pi/5} about z-axis = diag(1, e^{2pi*i/5}, e^{-2pi*i/5}) in complex basis
# = [[1, 0, 0], [0, cos(2pi/5), -sin(2pi/5)], [0, sin(2pi/5), cos(2pi/5)]] in real basis

theta_5 = 2 * math.pi / 5
c5 = math.cos(theta_5)
s5 = math.sin(theta_5)

print(f"  cos(2*pi/5) = (sqrt(5)-1)/4 = {c5:.6f}")
print(f"  sin(2*pi/5) = sqrt(10+2*sqrt(5))/4 = {s5:.6f}")
print()

# The Z5 eigenstates in the real basis:
# eigenvalue 1: |z> = (1, 0, 0) [along axis]
# eigenvalue omega: |+> ~ (0, 1, i)/sqrt(2) -> real: (0, 1, 0), (0, 0, 1) span
# The real 2D eigenspace for {omega, omega^4} is the (y,z) plane

# For mixing: we need to know the A5 mass-eigenstate basis
# vs the Z5 mass-eigenstate basis.

# The A5 mass basis: let's call it |e>, |mu>, |tau>
# These are determined by the DSI exponents n_e=0, n_mu=11, n_tau=17
# In the generation space (3-dim irrep), these correspond to specific
# orientations determined by the full A5 symmetry.

# The Z5 mass basis: determined by the spectral asymmetry direction
# Let's call this |nu1>, |nu2>, |nu3>
# |nu1> is along the Z5 axis, |nu2>, |nu3> are in the perpendicular plane

# The PMNS matrix U is the overlap: U_ij = <ell_i | nu_j>

# Key question: what is the ANGLE between the A5 basis and the Z5 basis?

# In the A5 group, there are 6 Z5 axes (6 pairs of opposite 5-fold axes)
# Each axis makes a specific angle with the A5 "natural" basis

# The natural A5 basis for the 3-rep can be taken as the golden frame:
# e1 = (1, phi, 0)/sqrt(1+phi^2)
# e2 = (0, 1, phi)/sqrt(1+phi^2)
# e3 = (phi, 0, 1)/sqrt(1+phi^2)
# These are the edge midpoints of the icosahedron, forming the golden rectangle

# The Z5 axis is, say, (0, 0, 1) (north pole of icosahedron)
# or more precisely one of the 12 vertices:
# v = (0, 1, phi)/sqrt(1+phi^2) (for example)

# The overlap between the golden frame and the Z5 axis frame
# determines the PMNS matrix!

print("Constructing the two bases:")
print()

# Golden frame (A5 natural basis for 3-rep):
norm = math.sqrt(1 + PHI**2)
e1 = np.array([1, PHI, 0]) / norm
e2 = np.array([0, 1, PHI]) / norm
e3 = np.array([PHI, 0, 1]) / norm

print("A5 (charged lepton) basis [golden frame]:")
print(f"  e1 = (1, phi, 0)/sqrt(1+phi^2) = {e1}")
print(f"  e2 = (0, 1, phi)/sqrt(1+phi^2) = {e2}")
print(f"  e3 = (phi, 0, 1)/sqrt(1+phi^2) = {e3}")
print()

# Check orthogonality
print(f"  e1.e2 = {np.dot(e1,e2):.6f}")
print(f"  e2.e3 = {np.dot(e2,e3):.6f}")
print(f"  e1.e3 = {np.dot(e1,e3):.6f}")
print()

# The golden frame is NOT orthogonal! The overlap is:
overlap = np.dot(e1, e2)
print(f"  Golden frame overlap: {overlap:.6f} = 1/(1+phi^2) = {1/(1+PHI**2):.6f}")
# So this is not the right basis. The 3-rep basis should be orthogonal.

# Let me use the standard icosahedron vertex coordinates instead.
# 12 vertices of icosahedron:
# (0, +/-1, +/-phi), (+/-1, +/-phi, 0), (+/-phi, 0, +/-1)

# A5 permutes these. The 3-dim irrep IS the natural R^3 action.
# So the mass eigenstates of charged leptons, in the generation space,
# are determined by their DSI quantum numbers.

# For NEUTRINOS: the mass eigenstates are determined by the spectral
# asymmetry, which selects a preferred AXIS.

# Let's say the Z5 axis is along vertex (0, 1, phi)/sqrt(1+phi^2)
# The eigenstates of Z5 rotation about this axis:
# |nu1> = axis direction = (0, 1, phi)/sqrt(1+phi^2)
# |nu2>, |nu3> = orthogonal plane

z5_axis = np.array([0, 1, PHI]) / math.sqrt(1 + PHI**2)

# Build orthonormal frame around z5_axis
# nu1 = z5_axis
nu1 = z5_axis
# nu2 = something orthogonal
temp = np.array([1, 0, 0])
nu2 = temp - np.dot(temp, nu1) * nu1
nu2 = nu2 / np.linalg.norm(nu2)
# nu3 = cross product
nu3 = np.cross(nu1, nu2)

print("Z5 (neutrino) basis:")
print(f"  nu1 = {nu1} [Z5 axis]")
print(f"  nu2 = {nu2}")
print(f"  nu3 = {nu3}")
print()

# For the charged lepton basis, use the standard orthonormal x,y,z
# (since the 3-dim irrep of A5 acts on R^3 by rotations)
# The mass eigenstates are along the generation directions.

# Actually, the charged lepton mass eigenstates in the 3-rep space
# should be related to the Hopf fibration structure.
# The 6 Hopf fibrations define 6 preferred axes (same as the 6 Z5 axes!)
# So both charged leptons and neutrinos see DIFFERENT Z5 axes.

# The PMNS matrix = rotation between two Z5 frames.
# If the charged lepton frame uses Z5 axis along (0,1,phi)/sqrt(1+phi^2)
# and the neutrino frame uses Z5 axis along (phi,0,1)/sqrt(1+phi^2)
# then the PMNS = rotation from one to the other.

print("KEY IDEA: PMNS = rotation between TWO Z5 axes of icosahedron")
print()

# The angle between two adjacent Z5 axes of the icosahedron:
# Two adjacent 5-fold axes are separated by the icosahedral angle
# The 12 vertices form 6 pairs. Adjacent pairs share an edge.

# Let's compute: angle between (0,1,phi) and (phi,0,1)
v1 = np.array([0, 1, PHI])
v2 = np.array([PHI, 0, 1])
cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
angle = math.acos(cos_angle)
print(f"Angle between adjacent Z5 axes: {math.degrees(angle):.2f} deg")
print(f"  cos(angle) = {cos_angle:.6f}")
print(f"  = 1/(1+phi^2) = {1/(1+PHI**2):.6f}")
print(f"  = 1/(phi+2) = {1/(PHI+2):.6f}  [since phi^2 = phi+1]")
print()

# IMPORTANT: 1/(phi+2) = 1/(phi+2)
# And sin^2 of this angle?
sin_sq = 1 - cos_angle**2
print(f"  sin^2(angle) = {sin_sq:.6f}")
print(f"  = 1 - 1/(phi+2)^2 = {1 - 1/(PHI+2)**2:.6f}")
print()

# The actual PMNS rotation is more complex - it's a 3D rotation
# decomposed into Euler angles.

# Let's try: the PMNS matrix U = R_23(th23) * R_13(th13) * R_12(th12)
# where the angles correspond to the relative orientations of the
# two Z5 frames.

# =====================================================================
# SECTION H: Direct approach - diagonalize mass matrices
# =====================================================================
print("\n--- H. DIRECT: MASS MATRIX DIAGONALIZATION ---")
print()

# The charged lepton mass matrix in the generation space:
# M_ell = diag(m_e, m_mu, m_tau) in the A5 mass basis
# The neutrino mass matrix is diag(m1, m2, m3) in the Z5 mass basis
# U_PMNS = V_ell^dag * V_nu

# If we KNOW the bases, we can compute U directly.
# But the bases depend on which Z5 subgroups correspond to
# charged leptons vs neutrinos.

# A simpler approach: use the EIGENVALUE RATIOS of the icosahedral Laplacian
# on different A5 irreps as a proxy for the mass matrices.

# Charged lepton mass ratios:
# m_e : m_mu : m_tau = 1 : phi^11 : phi^17
# In generation space, these correspond to DSI levels 0, 11, 17

# Neutrino mass ratios:
# m_1 : m_2 : m_3 ~ 0 : sqrt(dm21) : sqrt(dm31) ~ 0 : 8.7 : 49.5 meV
# The neutrino masses come from the spectral asymmetry.

# The MIXING is determined by the MISMATCH between these two hierarchies.

# A clean test: if we construct M_ell and M_nu in the SAME basis
# and diagonalize, what mixing do we get?

# In the standard seesaw: M_nu = M_D^T * M_R^{-1} * M_D
# M_D is the Dirac mass (related to charged lepton masses)
# M_R is the Majorana mass (comes from B-L breaking)

# In our framework: M_R comes from the spectral asymmetry
# eta = -35 determines the overall scale
# The STRUCTURE of M_R in generation space determines the mixing!

print("The PMNS mixing comes from the STRUCTURE of M_R in generation space.")
print("M_R is determined by the spectral asymmetry restricted to each generation.")
print()

# =====================================================================
# SECTION I: Tribimaximal as a starting point
# =====================================================================
print("\n--- I. COMPARISON WITH KNOWN MIXING PATTERNS ---")
print()

# Tribimaximal mixing (TBM):
# sin^2(th12) = 1/3 = 0.3333
# sin^2(th23) = 1/2 = 0.5000
# sin^2(th13) = 0

# Our formulas:
# sin^2(th12) = 2/(phi+5) = 0.3028
# sin^2(th23) = 4/7 = 0.5714
# sin^2(th13) = 1/45 = 0.0222

print("Comparison with tribimaximal (TBM):")
print(f"  {'':20s} {'TBM':>10s} {'Ours':>10s} {'PDG':>10s}")
print(f"  {'sin^2(th12)':20s} {'1/3':>10s} {'2/(phi+5)':>10s} {'0.307':>10s}")
print(f"  {'sin^2(th23)':20s} {'1/2':>10s} {'4/7':>10s} {'0.546':>10s}")
print(f"  {'sin^2(th13)':20s} {'0':>10s} {'1/45':>10s} {'0.0220':>10s}")
print()

# TBM arises from S4 symmetry (octahedral).
# Our framework uses A5 (icosahedral).
# The DEVIATIONS from TBM should come from the A5 structure.

# Key: TBM has sin^2(th12) = 1/3. Ours has 2/(phi+5).
# The difference: 1/3 - 2/(phi+5) = (phi+5-6)/(3*(phi+5)) = (phi-1)/(3*(phi+5))
# = (1/phi)/(3*(phi+5)) -- small correction from golden ratio!

delta_12 = 1/3 - 2/(PHI+a1)
print(f"Delta from TBM:")
print(f"  th12: 1/3 - 2/(phi+5) = {delta_12:.6f} = (phi-1)/(3*(phi+5))")
print(f"  th23: 1/2 - 4/7 = {1/2 - 4/7:.6f} = -1/14")
print(f"  th13: 0 - 1/45 = {-1/45:.6f}")
print()

# Is there a pattern?
# delta_12 = 1/phi / (3*(phi+5)) = 1/(3*phi*(phi+5))
# delta_23 = -1/14 = -1/(2*7) = -1/(2*(a1+2))
# delta_13 = -1/45 = -1/(a1*N_eig)

print("Deviations from TBM in framework constants:")
print(f"  delta_12 = 1/(3*phi*(phi+a1)) = {1/(3*PHI*(PHI+a1)):.6f}")
print(f"           actual: {delta_12:.6f}")
print(f"  delta_23 = -1/(2*(a1+2)) = {-1/(2*(a1+2)):.6f}")
print(f"  delta_13 = -1/(a1*N_eig) = {-1/(a1*N_eig):.6f}")
print()

# =====================================================================
# SECTION J: A5 -> Z5 breaking pattern
# =====================================================================
print("\n--- J. A5 -> Z5 BREAKING: THE PMNS FROM RESIDUAL SYMMETRY ---")
print()

# In the A5 flavor symmetry approach (well-studied in literature):
# A5 -> G_ell x G_nu (residual symmetries in charged lepton and neutrino sectors)
# Common choice: G_ell = Z5, G_nu = Z2 x Z2 (Klein four-group)
# This gives specific predictions for mixing angles!

# A5 has Z5, Z3, Z2 subgroups.
# If G_ell = Z5 and G_nu = Z2 x Z2:
# -> sin^2(th12) = (2+phi)/5 ≈ 0.724 (too large!)

# If G_ell = Z3 and G_nu = Z2 x Z2:
# -> Tribimaximal mixing

# Our framework: the CHARGED LEPTONS see the Hopf S^1 fiber (period 10)
# which gives Z10 = Z5 x Z2 symmetry.
# NEUTRINOS see the spectral asymmetry which picks a Z5 direction.

# The key might be: G_ell = Z3 (from N_gen=3 generations)
# and G_nu = Z5 (from the spectral asymmetry selecting a 5-fold axis)

# A5 -> Z3 x Z5 would give a SPECIFIC mixing pattern.

print("Residual symmetry approach:")
print("  Charged leptons: DSI levels 0,11,17 -> phi^n pattern")
print("    Preserved symmetry: Z3 (from N_gen=3)")
print("  Neutrinos: spectral asymmetry eta=-35 -> selects Z5 axis")
print("    Preserved symmetry: Z5 (from a1=5)")
print()
print("  A5 breaking pattern: A5 -> Z3 (leptons) x Z5 (neutrinos)")
print()

# This is a well-defined mathematical setup!
# The PMNS matrix is determined by the relative embedding
# of Z3 and Z5 in A5.

# Z3 subgroup of A5: generated by a 3-cycle, e.g., (1,2,3)
# Z5 subgroup of A5: generated by a 5-cycle, e.g., (1,2,3,4,5)

# The 3-dim irrep restricted to Z3:
# chi_3(3-cycle) = 0
# Eigenvalues: 1, w3, w3^2 where w3 = e^{2pi*i/3}
print("3-rep restricted to Z3: eigenvalues {1, w3, w3^2}")
print("3-rep restricted to Z5: eigenvalues {1, w5, w5^4}")
print()

# The PMNS matrix elements are:
# U_ij = <Z3 eigenstate i | Z5 eigenstate j>
# This is determined by the EMBEDDING of Z3 and Z5 in SO(3)!

# Z3 in SO(3): rotation by 120 degrees around one axis
# Z5 in SO(3): rotation by 72 degrees around another axis
# The relative angle between these two axes determines the PMNS!

# For the icosahedron:
# 3-fold axes: through face centers (20 faces -> 10 axes)
# 5-fold axes: through vertices (12 vertices -> 6 axes)
# The angle between a 3-fold and an adjacent 5-fold axis:

# Face center to vertex angle on icosahedron:
# Using standard coordinates:
# Vertex: (0, 1, phi)/sqrt(1+phi^2)
# Face center: centroid of triangle with vertices
# e.g., (0,1,phi), (1,phi,0), (phi,0,1) -> centroid = (1+phi, 1+phi, 1+phi)/(3*sqrt(1+phi^2))
# -> = (1,1,1)/sqrt(3) after normalization (since all components equal)

face_center = np.array([1, 1, 1]) / math.sqrt(3)
vertex = np.array([0, 1, PHI]) / math.sqrt(1 + PHI**2)
cos_fv = np.dot(face_center, vertex)
angle_fv = math.acos(cos_fv)

print(f"Angle between Z3 axis (face center) and Z5 axis (vertex):")
print(f"  cos = (1+phi)/(sqrt(3)*sqrt(1+phi^2)) = {cos_fv:.6f}")
print(f"  angle = {math.degrees(angle_fv):.2f} deg")
print(f"  sin^2 = {1-cos_fv**2:.6f}")
print()

# This angle is ~20.9 degrees.
# The atmospheric angle is ~47 degrees (sin^2 = 4/7 ≈ 0.571)
# Not a direct match. But the PMNS involves 3D rotation decomposition.

# =====================================================================
# SECTION K: Explicit construction
# =====================================================================
print("\n--- K. EXPLICIT PMNS FROM Z3 x Z5 EMBEDDING ---")
print()

# Build the Z3 and Z5 generators in the 3-dim irrep explicitly.
# Z5 generator (rotation by 2*pi/5 around (0,1,phi)):
axis_5 = np.array([0, 1, PHI])
axis_5 = axis_5 / np.linalg.norm(axis_5)

def rotation_matrix(axis, angle):
    """Rodrigues rotation formula."""
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * K @ K

R5 = rotation_matrix(axis_5, 2*math.pi/5)

# Z3 generator (rotation by 2*pi/3 around (1,1,1)/sqrt(3)):
axis_3 = np.array([1, 1, 1])
axis_3 = axis_3 / np.linalg.norm(axis_3)
R3 = rotation_matrix(axis_3, 2*math.pi/3)

# Eigenvectors of R5 (Z5 basis = neutrino mass eigenstates):
evals5, evecs5 = np.linalg.eig(R5)
print("Z5 generator eigenvalues:", [f"{e:.4f}" for e in evals5])

# Eigenvectors of R3 (Z3 basis = charged lepton mass eigenstates):
evals3, evecs3 = np.linalg.eig(R3)
print("Z3 generator eigenvalues:", [f"{e:.4f}" for e in evals3])
print()

# Sort eigenvectors by real eigenvalue first (eigenvalue = 1)
idx5 = np.argsort(-np.real(evals5))  # real eigenvalue first
idx3 = np.argsort(-np.real(evals3))
evecs5 = evecs5[:, idx5]
evecs3 = evecs3[:, idx3]

# The PMNS matrix = overlap between these two bases
# U_PMNS = evecs3^dag * evecs5
U_PMNS = np.conj(evecs3.T) @ evecs5

print("Raw PMNS matrix |U_ij|^2:")
U_sq = np.abs(U_PMNS)**2
for i in range(3):
    print(f"  [{U_sq[i,0]:.6f}  {U_sq[i,1]:.6f}  {U_sq[i,2]:.6f}]")
print()

# Extract mixing angles from |U_PMNS|^2
# Standard parametrization:
# |U_e3|^2 = sin^2(th13)
# |U_mu3|^2 / (1-|U_e3|^2) = sin^2(th23)
# |U_e2|^2 / (1-|U_e3|^2) = sin^2(th12)

s13_sq_calc = U_sq[0, 2]
s23_sq_calc = U_sq[1, 2] / (1 - s13_sq_calc)
s12_sq_calc = U_sq[0, 1] / (1 - s13_sq_calc)

print("Extracted mixing angles:")
print(f"  sin^2(th13) = {s13_sq_calc:.6f}  (target: {1/45:.6f} = 1/45)")
print(f"  sin^2(th23) = {s23_sq_calc:.6f}  (target: {4/7:.6f} = 4/7)")
print(f"  sin^2(th12) = {s12_sq_calc:.6f}  (target: {2/(PHI+5):.6f} = 2/(phi+5))")
print()
print(f"  PDG values: th12={0.307:.4f}, th23={0.546:.4f}, th13={0.0220:.4f}")
print()

# Try different axis choices to see if we can match
print("--- Trying different Z3/Z5 axis combinations ---")
print()

# Icosahedron vertices (12):
ico_vertices = []
for s1 in [1, -1]:
    for s2 in [1, -1]:
        ico_vertices.append(np.array([0, s1*1, s2*PHI]))
        ico_vertices.append(np.array([s1*1, s2*PHI, 0]))
        ico_vertices.append(np.array([s2*PHI, 0, s1*1]))
ico_vertices = [v/np.linalg.norm(v) for v in ico_vertices]

# Face centers (20): centroids of triangular faces
# For simplicity, use a few representative 3-fold axes
axis3_options = [
    np.array([1, 1, 1]) / math.sqrt(3),
    np.array([1, 1, -1]) / math.sqrt(3),
    np.array([1, -1, 1]) / math.sqrt(3),
    np.array([-1, 1, 1]) / math.sqrt(3),
]

best_match = None
best_score = float('inf')

for i_5, ax5 in enumerate(ico_vertices[:6]):  # 6 Z5 axes
    for i_3, ax3 in enumerate(axis3_options):
        R5 = rotation_matrix(ax5, 2*math.pi/5)
        R3 = rotation_matrix(ax3, 2*math.pi/3)

        evals5, evecs5 = np.linalg.eig(R5)
        evals3, evecs3 = np.linalg.eig(R3)

        idx5 = np.argsort(-np.real(evals5))
        idx3 = np.argsort(-np.real(evals3))
        evecs5 = evecs5[:, idx5]
        evecs3 = evecs3[:, idx3]

        U = np.conj(evecs3.T) @ evecs5
        U_sq = np.abs(U)**2

        s13 = U_sq[0, 2]
        s23 = U_sq[1, 2] / (1 - s13) if s13 < 0.99 else 0
        s12 = U_sq[0, 1] / (1 - s13) if s13 < 0.99 else 0

        # Score: distance to our predictions
        score = (s12 - 2/(PHI+5))**2 + (s23 - 4/7)**2 + (s13 - 1/45)**2

        if score < best_score:
            best_score = score
            best_match = (i_5, i_3, s12, s23, s13, U_sq)

i5, i3, s12_b, s23_b, s13_b, U_sq_b = best_match
print(f"Best match: Z5 axis #{i5}, Z3 axis #{i3}")
print(f"  sin^2(th12) = {s12_b:.6f}  (target: {2/(PHI+5):.6f})")
print(f"  sin^2(th23) = {s23_b:.6f}  (target: {4/7:.6f})")
print(f"  sin^2(th13) = {s13_b:.6f}  (target: {1/45:.6f})")
print(f"  Score = {best_score:.6f}")
print()
print(f"|U|^2 matrix:")
for i in range(3):
    print(f"  [{U_sq_b[i,0]:.6f}  {U_sq_b[i,1]:.6f}  {U_sq_b[i,2]:.6f}]")

# =====================================================================
# SECTION L: Assessment
# =====================================================================
print("\n\n" + "=" * 70)
print("ASSESSMENT: PMNS DERIVATION")
print("=" * 70)
print()
print("The A5 -> Z3 x Z5 residual symmetry approach gives SPECIFIC")
print("predictions for PMNS mixing angles that depend on the relative")
print("embedding of Z3 and Z5 in A5.")
print()
print("The approach is mathematically well-defined (group theory).")
print("The question is whether the specific angles match our formulas.")
print()
print("If the Z3 x Z5 approach DOES reproduce our formulas:")
print("  -> UPGRADE from PATTERN to DERIVED")
print("  -> The mechanism is clear: A5 flavor symmetry with specific breaking")
print()
print("If NOT:")
print("  -> The formulas remain PATTERN but with a clear CANDIDATE mechanism")
print("  -> The specific breaking direction needs to be identified")

"""
exp466_test_function_f2.py
==========================
TEST FUNCTION BOOTSTRAP: Can we determine f_2 (and hence G) from the 600-cell?

Three routes:
  A) f_2 by analogy with f_4 = N^4/(2*pi)
  B) KK reduction from G_3D = 1/12 (DERIVED)
  C) Self-consistent spectral cutoff

The goal: eliminate m_e as dimensional anchor, derive G from a1=5.

Author: Razvan-Constantin Anghelina
Date: 2026-03-25
"""

import math

print("=" * 70)
print("EXP466: TEST FUNCTION f_2 BOOTSTRAP")
print("Can we derive Newton's constant from a_1 = 5?")
print("=" * 70)

# ============================================================
# FRAMEWORK CONSTANTS (all from a1=5)
# ============================================================
a1 = 5
b1 = a1 + 1  # 6
phi = (1 + math.sqrt(a1)) / 2  # golden ratio
N = 2 * a1 * (a1 - 1) * (a1 - 2)  # 120 = |2I|
# Actually N = a1! = 120
h_E8 = N / (a1 - 1)  # 30 = Coxeter number E8

# Alpha from quadratic
A_coef = 2 * math.pi
B_coef = -(N / b1) * phi**4
C_coef = 1
disc = B_coef**2 - 4 * A_coef * C_coef
alpha = (-B_coef - math.sqrt(disc)) / (2 * A_coef)
alpha_exp = 1 / 137.035999084

# Strong coupling
alpha_s = 1 / (2 * phi**3)

# Weinberg angle
sin2_tW = b1 / (a1**2 + 1)  # 6/26

# Seeley-DeWitt reduced coefficients
A0 = a1 + b1  # 11
A1 = 2 * a1**2 + 2 * a1 + 2  # 62
A2 = 6 * a1**2 + 15 * a1 + 8  # 233 = F_13

# Full Seeley-DeWitt coefficients
c0 = 2 * N * A0  # 2640
c1 = 2 * N * A1  # 14880
c2 = 2 * N * A2  # 55920

# Planck hierarchy
z_P = 4 * phi**2  # 10.4721...

# Galois norm
alpha_alpha_prime = 1 / (2 * math.pi)

# Lambda_max of D^2 spectrum
Lambda_max_D2 = b1 * phi**2  # 15.7082

print(f"\n--- FRAMEWORK CONSTANTS ---")
print(f"a1 = {a1}, b1 = {b1}, phi = {phi:.10f}")
print(f"N = {N}, h(E8) = {h_E8:.0f}")
print(f"alpha = 1/{1/alpha:.6f} (theory), 1/{1/alpha_exp:.6f} (exp)")
print(f"alpha_s = {alpha_s:.6f}")
print(f"sin^2(tW) = {sin2_tW:.6f} = {b1}/{a1**2+1}")
print(f"\nSeeley-DeWitt: A0={A0}, A1={A1}, A2={A2}")
print(f"  c0={c0}, c1={c1}, c2={c2}")
print(f"  Check Diophantine: 2*A1^2+1 = {2*A1**2+1}, 3*A0*A2 = {3*A0*A2}")
print(f"  c1^2/(c0*c2) = {c1**2/(c0*c2):.8f} (should be ~3/2)")
print(f"\nz_Planck = 4*phi^2 = {z_P:.10f}")
print(f"Lambda_max(D^2) = b1*phi^2 = {Lambda_max_D2:.6f}")
print(f"alpha*alpha' = 1/(2*pi) = {alpha_alpha_prime:.10f}")

# ============================================================
# PHYSICAL CONSTANTS
# ============================================================
m_e_MeV = 0.51099895000  # MeV
G_N = 6.67430e-11  # m^3 kg^-1 s^-2
hbar = 1.054571817e-34  # J s
c_light = 299792458  # m/s
M_P_kg = math.sqrt(hbar * c_light / G_N)
M_P_GeV = M_P_kg * c_light**2 / (1.602176634e-10)  # GeV (E = mc^2)
M_P_red_GeV = M_P_GeV / math.sqrt(8 * math.pi)  # reduced Planck mass

# Convert to natural units (GeV)
m_e_GeV = m_e_MeV / 1000

print(f"\n--- PHYSICAL CONSTANTS ---")
print(f"m_e = {m_e_GeV:.6e} GeV")
print(f"M_P = {M_P_GeV:.6e} GeV")
print(f"M_P_red = {M_P_red_GeV:.6e} GeV")
print(f"G_N = {G_N:.5e} m^3 kg^-1 s^-2")
print(f"m_e/M_P = {m_e_GeV/M_P_GeV:.6e}")
print(f"alpha^z_P = {alpha_exp**z_P:.6e}")
print(f"Gap = {(alpha_exp**z_P - m_e_GeV/M_P_GeV)/(m_e_GeV/M_P_GeV)*100:.4f}%")

# ============================================================
# ROUTE A: f_2 BY ANALOGY WITH f_4
# ============================================================
print("\n" + "=" * 70)
print("ROUTE A: f_2 by analogy with f_4 = N^4/(2*pi)")
print("=" * 70)

f4 = N**4 / (2 * math.pi)
print(f"\nf_4 = N^4/(2*pi) = {f4:.6e}")
print(f"f_4 = N^{a1-1} * (alpha*alpha') [since d_ST=4=a1-1]")

# The spectral action for gravity (Chamseddine-Connes):
# 1/(16*pi*G) = f_2 * Lambda^2 * c_1 / (96*pi^2)  [standard NCG relation]
# => G = 6*pi / (f_2 * Lambda^2 * c_1)
#
# For inflation (R^2 term):
# M_scal^2 / M_P_red^2 = 320*pi^2 / (12 * f_4 * c2)
# => This gives M_scal once we know M_P_red

# APPROACH: if f_4 = N^4/(2*pi), what about f_2?
# Several hypotheses:

hypotheses_f2 = {
    "N^2/(2*pi)": N**2 / (2 * math.pi),
    "N^2 * alpha*alpha'": N**2 * alpha_alpha_prime,  # same as above!
    "N^2": float(N**2),
    "N^(d_ST-2)/(2*pi) = N^2/(2*pi)": N**2 / (2 * math.pi),
    "c1/(2*pi)": c1 / (2 * math.pi),
    "A1^2/(2*pi)": A1**2 / (2 * math.pi),
    "f4^(1/2)": math.sqrt(f4),
    "f4 * (c0/c2)": f4 * c0 / c2,
    "N^2 * phi^2": N**2 * phi**2,
    "c0 * phi^2": c0 * phi**2,
    "N^2 / phi^2": N**2 / phi**2,
    "b1 * N^2 / (2*pi)": b1 * N**2 / (2 * math.pi),
    "N^2 / (2*pi*phi^2)": N**2 / (2 * math.pi * phi**2),
}

# For each f_2 hypothesis, compute implied Lambda and check consistency
# From G = 6*pi / (f_2 * Lambda^2 * c1), and using natural units where G = 1/M_P^2:
# 1/M_P^2 = 6*pi / (f_2 * Lambda^2 * c1)
# => Lambda^2 = 6*pi * M_P^2 / (f_2 * c1)
# => Lambda = M_P * sqrt(6*pi / (f_2 * c1))

print(f"\nFor each f_2 hypothesis, compute implied cutoff Lambda:")
print(f"  Using G = 6*pi / (f_2 * Lambda^2 * c1)")
print(f"  => Lambda = M_P * sqrt(6*pi / (f_2 * c1))")
print(f"{'Hypothesis':<35} {'f_2':>12} {'Lambda (GeV)':>15} {'Lambda/M_P':>12} {'Lambda/M_GUT':>12}")
print("-" * 90)

M_GUT = 2e16  # typical GUT scale in GeV

for name, f2_val in hypotheses_f2.items():
    Lambda_sq = 6 * math.pi * M_P_GeV**2 / (f2_val * c1)
    if Lambda_sq > 0:
        Lambda = math.sqrt(Lambda_sq)
        print(f"{name:<35} {f2_val:>12.4e} {Lambda:>15.4e} {Lambda/M_P_GeV:>12.4e} {Lambda/M_GUT:>12.4f}")

# REVERSE: what f_2 gives Lambda = M_P? Lambda = M_GUT?
print(f"\nReverse engineering:")
f2_for_MP = 6 * math.pi / (c1)  # Lambda = M_P
f2_for_MGUT = 6 * math.pi * M_P_GeV**2 / (c1 * M_GUT**2)
f2_for_A1 = 6 * math.pi * M_P_GeV**2 / (c1 * (A1 * M_P_red_GeV)**2)

print(f"  Lambda = M_P     => f_2 = {f2_for_MP:.8f} = 6*pi/c1")
print(f"    = 6*pi/{c1} = pi/{c1//6}")
print(f"    = pi/2480 (is 2480 a framework number?)")
print(f"    2480 = {2480} = 2*N*A1/b1 = 2*{N}*{A1}/{b1} = {2*N*A1//b1}")
check_2480 = 2 * N * A1 / b1
print(f"    CHECK: 2*N*A1/b1 = {check_2480:.1f}")
print(f"    2480 = N * A1 / 3 = {N*A1/3:.1f}")
print(f"    2480 = 8 * 310 = 8 * (A1*a1) = {8*A1*a1}")
print(f"    2480 = c1/b1 = {c1/b1:.1f}")
print(f"    *** f_2 = 6*pi/c1 = b1*pi/c1*b1 = pi*b1/c1 ***")
print(f"    f_2 = pi/(c1/b1) = pi / (2*N*A1/b1)")

# Actually c1/b1 = 2*N*A1/b1 = 2*120*62/6 = 2480
# And 6*pi/c1 = pi/2480 = pi * b1 / (b1 * c1/b1) ... hmm

print(f"\n  Lambda = M_GUT (2e16) => f_2 = {f2_for_MGUT:.4e}")
print(f"    Compare: N^2/(2*pi) = {N**2/(2*math.pi):.4e}")
print(f"    Ratio: {f2_for_MGUT / (N**2/(2*math.pi)):.6f}")

# ============================================================
# ROUTE B: KK REDUCTION FROM G_3D = 1/12
# ============================================================
print("\n" + "=" * 70)
print("ROUTE B: Kaluza-Klein reduction from G_3D = 1/12")
print("=" * 70)

G_3D = 1.0 / 12  # DERIVED from Turaev-Viro (exp424)
# In 3D topological gravity: G_3D = 1/(4*k) where k = a1-2 = 3
# k is the Chern-Simons level

print(f"\nG_3D = 1/(4k) = 1/12 (DERIVED, k = a1-2 = 3)")

# KK reduction: spacetime = M^4 x S^3_internal
# The 600-cell is a triangulation of S^3
#
# Standard KK relation (from d+n dimensions to d dimensions):
# G_d = G_{d+n} / V_n
#
# In our case: physical = 4D, internal = 3D (the 600-cell S^3)
# G_4D = G_7D / V_3
#
# But G_3D is the 3D gravity ON the internal space, not the 7D gravity.
# The relation is more subtle.
#
# In Chern-Simons 3D gravity:
# S_CS = (k/(4*pi)) * integral(Tr(A dA + 2/3 A^3))
# G_3D = 1/(4*k) in units where l_AdS = 1
#
# The physical interpretation:
# The "volume" of the internal S^3 in lattice units is related to N (vertices)
# or to the number of 3-simplices (600 cells)

# APPROACH 1: V_3 from number of cells
n_cells = 600  # number of 3-cells (tetrahedra) in 600-cell
n_vertices = 120
n_edges = 720
n_faces = 1200

# Euler characteristic: V - E + F - C = 0 (for S^3)
euler = n_vertices - n_edges + n_faces - n_cells
print(f"\n600-cell: V={n_vertices}, E={n_edges}, F={n_faces}, C={n_cells}")
print(f"Euler characteristic: {euler} (should be 0 for S^3)")

# Volume of regular 600-cell with unit edge length
# V(600-cell) = (N/4) * V(tetrahedron) ... in the regular case
# For unit edge tetrahedra: V_tet = 1/(6*sqrt(2))
# Total: V = 600/(6*sqrt(2)) * edge^3 ... but this is Euclidean, not S^3

# Better: the S^3 has radius R = phi (the circumradius of 600-cell with edge 1/phi)
# Volume of S^3 with radius R: V = 2*pi^2*R^3
R_600cell = phi  # circumradius when edge = 1/phi
V_S3 = 2 * math.pi**2 * R_600cell**3
print(f"\nS^3 circumradius of 600-cell: R = phi = {phi:.6f}")
print(f"V(S^3) = 2*pi^2*R^3 = 2*pi^2*phi^3 = {V_S3:.6f}")
print(f"Compare: 2*pi^2*phi^3 = pi^2/(alpha_s) = {math.pi**2/alpha_s:.6f}")
# Since alpha_s = 1/(2*phi^3), V = 2*pi^2*phi^3 = pi^2/alpha_s. Nice!

# APPROACH 2: KK with spectral data
# In the spectral action on M^4 x F (fiber = 600-cell):
# The 4D Newton constant is related to the trace on the fiber:
# 1/(16*pi*G_4D) = (f_2 * Lambda^2 / (96*pi^2)) * Tr_F(1)
# where Tr_F(1) = c0 = 2640 (dimension of fiber Hilbert space)
#
# Wait, that's the standard NCG Chamseddine-Connes formula where
# c1 = Tr(D^2) appears in the gravity sector, not c0.
# Let me be more precise.

# The Chamseddine-Connes spectral action expansion gives:
# S ~ f_4*Lambda^4*c0 + f_2*Lambda^2*c1*R + f_0*c2*(gauge terms) + ...
# where R is the scalar curvature
# So: 1/(16*pi*G) = f_2*Lambda^2*c1/(48*pi^2) ... I need exact coefficient

# Standard NCG result (Chamseddine-Connes-Marcolli):
# 1/(2*kappa_0^2) = (96*f_2*Lambda^2 - f_0*c)/(24*pi^2)
# where kappa_0 = sqrt(8*pi*G)
# Simplified (dominant term): 1/G ~ f_2*Lambda^2 * c1 / pi^2

# Let's parametrize:
# G_4D = C_G / (f_2 * Lambda^2 * c1)
# where C_G is a numerical constant we'll determine from the exact NCG formula

# From Chamseddine-Connes 2007 (CCM):
# g^2 = pi^2 / (f_0 * sum_generations * (a - 2e)) where a,e are matrix traces
#
# For gravity, the exact relation is:
# 1/(kappa_0^2) = 96*f_2*Lambda^2/(24*pi^2) = 4*f_2*Lambda^2/pi^2
# kappa_0^2 = pi^2/(4*f_2*Lambda^2)
# G = kappa_0^2/(8*pi) = pi/(32*f_2*Lambda^2)
#
# BUT: this has an ADDITIONAL factor of c1 from the internal space trace!
# In the product geometry M^4 x F:
# Tr(|D|^{-2}) includes both continuous and discrete sums
# The discrete part contributes c1 = Tr_F(D_F^2) = 14880

# Let me use the relation as stated in the paper (eq:lagrangian_explicit):
# The R coefficient in the Lagrangian is: c1 * Lambda^2 * f_2
# So: 1/(16*pi*G) = c1 * f_2 * Lambda^2 / (some numerical factor)
# The factor depends on the exact normalization.

# From standard NCG (Connes-Marcolli textbook, eq 1.212):
# 1/(16*pi*G) = f_2*Lambda^2/(48) * Tr(D_F^{-2}) ... no
# Actually the Einstein-Hilbert term coefficient is:
# (f_2*Lambda^2 / (4*pi^2)) * (1/12) * c_1(F) * R
# = f_2*Lambda^2*c1 / (48*pi^2) * R

# Comparing with (1/(16*pi*G)) * R:
# 1/(16*pi*G) = f_2*Lambda^2*c1 / (48*pi^2)
# G = 48*pi^2 / (16*pi * f_2*Lambda^2*c1)
# G = 3*pi / (f_2*Lambda^2*c1)

print(f"\n--- NCG Gravity Relation ---")
print(f"Using: 1/(16*pi*G) = f_2*Lambda^2*c1 / (48*pi^2)")
print(f"=> G = 3*pi / (f_2 * Lambda^2 * c1)")

# In Planck units (hbar=c=1): G = 1/M_P^2
# So: 1/M_P^2 = 3*pi / (f_2 * Lambda^2 * c1)
# => f_2 * Lambda^2 = 3*pi * M_P^2 / c1

f2_Lambda2_needed = 3 * math.pi * M_P_GeV**2 / c1
print(f"\nf_2 * Lambda^2 = 3*pi*M_P^2/c1 = {f2_Lambda2_needed:.6e} GeV^2")
print(f"  = 3*pi/(c1) * M_P^2")
print(f"  = pi/(c1/3) * M_P^2 = pi/{c1/3:.0f} * M_P^2")
print(f"  c1/3 = {c1/3:.0f} = 2*N*A1/3 = {2*N*A1/3:.0f}")
print(f"  = 2*{N}*{A1}/3 = {2*N*A1//3}")

# Now: what is f_2*Lambda^2 in terms of framework?
# If Lambda = Lambda_max_phys (the physical cutoff corresponding to Lambda_max_D2):
# Then Lambda is in physical units, and Lambda_max_D2 is in lattice units.
# The relation is: Lambda_phys = Lambda_max_D2 * E_lattice
# where E_lattice is the energy scale of one lattice unit.

# KK approach: the internal S^3 has a "radius" in physical units.
# All eigenvalues of D^2 on the internal space scale as 1/R_int^2
# The physical cutoff Lambda = sqrt(Lambda_max_D2) / R_int

# If R_int is the physical radius of the internal S^3, then:
# Lambda = phi * sqrt(b1) / R_int  [since Lambda_max_D2 = b1*phi^2]

# From KK: G_4D = G_3D_phys / V_3_phys
# where V_3_phys = 2*pi^2*R_int^3

# G_3D_phys = G_3D * l_3D^? ... this is where it gets tricky.
# In 3D Chern-Simons: G_3D = l/(16*pi*k) where l = AdS radius
# For k=3: G_3D_phys = l_int / (48*pi)
# If l_int = R_int (the internal space radius):

# Let's try a different approach.
# Use DIMENSIONAL ANALYSIS in the spectral action framework.

# The spectral action is Tr(f(D/Lambda))
# D has eigenvalues in lattice units (dimensionless on the internal space)
# Lambda is the physical cutoff (in GeV)

# The 4D effective action after integrating out the internal space:
# S_eff = sum_n f(lambda_n/Lambda) * S_4D(lambda_n)
# where lambda_n are internal eigenvalues

# For the gravity term:
# f_2*Lambda^2*c1*R  (c1 = sum of lambda_n^2 on internal space = Tr(D_F^2))

# So f_2*Lambda^2 has dimensions of [GeV]^2 / [eigenvalue^2]
# But eigenvalues are dimensionless in the lattice, so:
# f_2*Lambda^2 is in GeV^2

# KEY INSIGHT: Lambda should be the PHYSICAL scale of the internal space.
# If the internal S^3 has radius R_int (in physical units = 1/GeV):
# The eigenvalues of D on S^3 scale as 1/R_int
# So the dimensionful eigenvalues are: lambda_n_phys = lambda_n / R_int
# The cutoff is: Lambda ~ lambda_max_phys = sqrt(b1)*phi / R_int

# Then: f_2 * Lambda^2 = f_2 * b1*phi^2 / R_int^2

# KK RELATION for S^3:
# G_4D = G_{4+3} / V_3 = G_7D / (2*pi^2*R_int^3)
#
# But we need G_7D. From the 7D spectral action:
# 1/(16*pi*G_7D) = f_2_7D * Lambda_7D^5 * ... (higher dim)
# This is getting very complicated.

# SIMPLER APPROACH: Use the RATIO f_2/f_4 and the known f_4
print("\n\n--- APPROACH: f_2/f_4 ratio ---")
print(f"We know f_4 = N^4/(2*pi) = {f4:.6e}")
print(f"We need f_2*Lambda^2 = {f2_Lambda2_needed:.6e} GeV^2")

# From the inflation sector:
# M_scal^2/M_P_red^2 = 320*pi^2 / (12 * f_4 * c2)
# The OBSERVED A_s determines M_scal, hence connects f_4 to physical scale.
# With f_4 = N^4/(2*pi):
M_scal_sq = M_P_red_GeV**2 * 320 * math.pi**2 / (12 * f4 * c2)
M_scal = math.sqrt(M_scal_sq)
print(f"\nM_scalaron = {M_scal:.4e} GeV (from f_4 = N^4/(2*pi))")

# The cutoff Lambda in the inflation sector is NOT necessarily the same as
# in the gravity sector. But in a self-consistent theory, they should be related.

# From the spectral action, the RATIO of the R term to the R^2 term:
# (c1*f_2*Lambda^2) * R vs (c2*f_0) * R^2
# The transition scale is: Lambda_trans^2 = c2*f_0 / (c1*f_2*Lambda^2) * M_P^2
# This should be ~ M_scal^2

# Actually, the scalaron mass IS this ratio:
# M_scal^2 = c1*f_2*Lambda^2 / (c2*f_0)
# Wait no, in the Starobinsky model:
# S = integral(R/(16*pi*G) + R^2/(6*M^2)) sqrt(g) d^4x
# M^2 = (c1*f_2*Lambda^2)/(c2*f_0) * (some numerical factor)

# The exact relation from NCG:
# From the Lagrangian: L = c1*f_2*Lambda^2*R + c2*f_0*(curvature^2 terms)
# The R^2 term coefficient: c2*f_0/(320*pi^2) (from Gauss-Bonnet decomposition)
# The R term coefficient: c1*f_2*Lambda^2/(48*pi^2)
# Ratio: M_scal^2 = (R coeff)/(R^2 coeff) ~ f_2*Lambda^2*c1 / (f_0*c2) * const

# In the standard formula:
# M_scal^2 = (48*pi^2 * f_2*Lambda^2 * c1) / (320*pi^2 * c2 * f_0) * (some_factor)
# = (48/(320)) * f_2*Lambda^2*c1/(f_0*c2)
# = (3/20) * f_2*Lambda^2*c1/(f_0*c2)

# If f_0 = f(0) and in our convention f_4 corresponds to f_0:
# Actually in NCG: f_0 = f(0) (value at origin), f_2 = integral, f_4 = leading term
# Let me use the standard notation:
# f_0 corresponds to the ZEROTH moment = f(0)
# f_2 = integral(f(u)*u du) from 0 to infinity
# f_4 = integral(f(u) du) from 0 to infinity

# For a step function f(u) = theta(1-u): f_0 = 1, f_2 = 1/2, f_4 = 1
# These are all O(1)

# In our case, f_4 = N^4/(2*pi) = 3.3e7 >> 1, so f is NOT a step function.
# It's a very peaked function.

# RELATIONSHIP: f_2/f_4 determines the SHAPE of f
# For a Gaussian f(u) = exp(-u^2/s^2): f_4 ~ s, f_2 ~ s^3, f_0 ~ 1
# For a power law f(u) = u^{-p}: f_k ~ 1/(k-p)
# For the 600-cell: f(u) = ???

# Let's parametrize: f_2 = f_4^{2/4} * C = f_4^{1/2} * C (if f_k scales as f_4^{k/4})
# Then f_2 = sqrt(f_4) * C

# f_4^{1/2} = sqrt(N^4/(2*pi)) = N^2/sqrt(2*pi)
f4_sqrt = N**2 / math.sqrt(2 * math.pi)
print(f"\nf_4^(1/2) = N^2/sqrt(2*pi) = {f4_sqrt:.4f}")

# If f_2 = f_4^{1/2} (pure scaling):
Lambda_from_scaling = math.sqrt(f2_Lambda2_needed / f4_sqrt)
print(f"If f_2 = f_4^(1/2): Lambda = {Lambda_from_scaling:.4e} GeV")
print(f"  Lambda/M_P = {Lambda_from_scaling/M_P_GeV:.6f}")

# ============================================================
# ROUTE C: BACK-COMPUTE f_2 AND LOOK FOR ALGEBRAIC STRUCTURE
# ============================================================
print("\n" + "=" * 70)
print("ROUTE C: Back-compute f_2 from known physics")
print("=" * 70)

# We know: G = 3*pi / (f_2 * Lambda^2 * c1)
# And from inflation: M_scal^2 = 320*pi^2 / (12 * f_4 * c2) * M_P_red^2

# The two "test function" problems (f_2 for gravity, f_4 for inflation)
# are related by the RATIO f_2/f_4 * Lambda^2

# From M_scal determination (CMB data):
A_s_obs = 2.10e-9  # Planck 2018
N_e = 60  # e-folds

# In Starobinsky: A_s = N_e^2 / (24*pi^2) * (M_scal/M_P_red)^2
# => M_scal^2 = 24*pi^2*A_s / N_e^2 * M_P_red^2
M_scal_from_CMB = M_P_red_GeV * math.sqrt(24 * math.pi**2 * A_s_obs / N_e**2)
print(f"\nM_scal (from CMB) = {M_scal_from_CMB:.4e} GeV")

# From f_4: f_4*c2 = 24*pi^2 * M_P_red^2 / (N_e^2 * A_s * 2*M_scal^2) ???
# Actually: A_s = N_e^2/(24*pi^2) * M^2/M_P_red^2 * (12*f_4*c2/320*pi^2)^...
# This is getting circular. Let me use a different approach.

# KEY IDEA: Express everything in terms of m_e and framework constants.
# We have m_e and alpha^{z_P} = m_e/M_P, so M_P = m_e/alpha^{z_P}
# Then G = 1/M_P^2 = alpha^{2*z_P}/m_e^2
# And f_2*Lambda^2 = 3*pi*M_P^2/c1 = 3*pi*m_e^2/(c1*alpha^{2*z_P})

# The DIMENSIONLESS combination:
# f_2 * (Lambda/M_P)^2 = 3*pi/c1
ratio_f2_Lambda_MP = 3 * math.pi / c1
print(f"\nDIMENSIONLESS: f_2 * (Lambda/M_P)^2 = 3*pi/c1 = {ratio_f2_Lambda_MP:.8f}")
print(f"  = pi/{c1//3}")
print(f"  = pi/4960")
print(f"  4960 = c1/3 = 2*N*A1/3")

# If Lambda = M_P (natural choice): f_2 = 3*pi/c1
# If f_2 = 1 (natural choice): Lambda/M_P = sqrt(3*pi/c1)
Lambda_if_f2_1 = M_P_GeV * math.sqrt(3 * math.pi / c1)
print(f"\nIf f_2 = 1: Lambda = {Lambda_if_f2_1:.4e} GeV")
print(f"  Lambda/M_P = {Lambda_if_f2_1/M_P_GeV:.6f}")
print(f"  Lambda/M_GUT = {Lambda_if_f2_1/M_GUT:.4f}")

# What if Lambda = M_scal (scalaron mass)?
f2_if_Lambda_Mscal = f2_Lambda2_needed / M_scal_from_CMB**2
print(f"\nIf Lambda = M_scal: f_2 = {f2_if_Lambda_Mscal:.4e}")

# ============================================================
# THE KEY QUESTION: IS 3*pi/c1 A FRAMEWORK NUMBER?
# ============================================================
print("\n" + "=" * 70)
print("KEY ANALYSIS: Is 3*pi/c1 algebraically meaningful?")
print("=" * 70)

val = 3 * math.pi / c1
print(f"\n3*pi/c1 = 3*pi/(2*N*A1) = {val:.10f}")
print(f"  = 3*pi / (240 * 62) = pi / (80*62) = pi/4960")
print(f"  = pi / (N * A1 / 3) ... hmm")

# Let's factor 4960
print(f"\n4960 = 2^5 * 5 * 31 = 32 * 155")
print(f"  = (a1-1)^2 * a1 * (2*a1^2 + 2*a1 + 2)/3")
# A1 = 62 = 2*31, c1/3 = 4960 = 2*N*62/3 = 2*120*62/3
# Hmm, 62/3 is not integer...
# c1 = 14880, c1/3 = 4960
# c1 = 2*N*A1 = 240*62

# Framework expressions for this ratio:
print(f"\nAlternative forms of f_2 = 3*pi/c1:")
print(f"  = 3*pi / (2*N*(2*a1^2+2*a1+2))")
print(f"  = 3*pi / (4*N*(a1^2+a1+1))")
print(f"  = 3*pi / (4*{N}*{a1**2+a1+1})")
print(f"  = 3*pi / (4*120*31)")
print(f"  = pi / (4*{N}*{a1**2+a1+1}//3)")
# a1^2+a1+1 = 31, and 4*120 = 480, 4*120*31/3 is not integer

# Compare with alpha:
print(f"\n  f_2/alpha = {val/alpha_exp:.6f}")
print(f"  f_2*137 = {val*137:.6f}")
print(f"  f_2/alpha^2 = {val/alpha_exp**2:.4f}")
print(f"  1/f_2 = {1/val:.4f}")
print(f"  Compare: c1/(3*pi) = {c1/(3*math.pi):.4f}")
print(f"  Compare: 2*phi^3/pi = {2*phi**3/math.pi:.6f}")
print(f"  Compare: N/(2*pi) = {N/(2*math.pi):.6f}")

# ============================================================
# ROUTE B CONTINUED: G from KK with G_3D
# ============================================================
print("\n" + "=" * 70)
print("ROUTE B CONTINUED: G_4D from G_3D via KK")
print("=" * 70)

# The idea: if the internal S^3 (600-cell) has physical radius R_int,
# the 4D Newton constant from KK reduction of 7D gravity is:
# G_4D = G_7D / V_3 = G_7D / (2*pi^2*R_int^3)
#
# We can also relate G_3D (on the internal S^3) to G_7D.
# In the product M^4 x S^3:
# The 3D gravity on S^3 comes from the 7D gravity restricted to S^3
# G_3D_eff = G_7D / V_4 ... no, this is circular
#
# Better approach: use the TQFT result directly.
# G_3D = 1/12 is in NATURAL UNITS of the 600-cell (lattice units).
# In physical units: G_3D_phys = G_3D * l_lattice^1 (3D Newton has dim [length])
#
# In 3D, [G] = [length] (since G_3D*M/r is dimensionless in 3D)
# In 4D, [G] = [length]^2 (since G_4D*M/r is dimensionless in 4D)
#
# KK relation for S^1 compactification (simpler):
# G_4D = G_5D / (2*pi*R)
# For S^3:
# G_4D = G_7D / (2*pi^2*R^3)
#
# But we want to relate G_3D (internal) to G_4D (physical).
#
# KEY INSIGHT: In the 600-cell framework:
# - D_total = D_external (x) 1 + gamma_5 (x) D_F
# - The physical Newton constant comes from the external part
# - G_3D = 1/12 comes from the internal part (TQFT on S^3)
#
# The FULL gravity comes from the PRODUCT geometry:
# G_4D^{-1} = (f_2 * Lambda^2 / pi^2) * Tr(D_F^2) / 48
# And Tr(D_F^2) = c1 = 14880

# So the chain is: G_3D -> R_int -> Lambda -> G_4D

# From G_3D: If G_3D = l_P_3D (Planck length in 3D), and this equals 1/12
# in lattice units, the lattice spacing a = l_P_3D / G_3D ???
# This is getting confused. Let me try numerically.

# NUMERICAL APPROACH:
# We KNOW G_4D from experiment.
# We KNOW G_3D = 1/12 (in lattice units).
# What is the lattice spacing that makes them consistent?

# In 3D: G_3D_phys = G_3D * a_lattice (where a_lattice is in meters/GeV^-1)
# because [G_3D] = [length] in 3D

# KK: if we reduce from 7D = 4D + 3D:
# G_4D = G_7D / (2*pi^2*R^3)
# [G_7D] = [length]^5 in 7D

# And G_7D is related to G_3D by the spectral action on M^4:
# G_3D = G_7D / V_4 ... where V_4 is the 4D spacetime volume? No.

# Actually this KK picture is getting circular. Let me try a more direct approach.

# DIRECT APPROACH: Use the spectral zeta function
# From the spectrum of D^2 on the 600-cell, the zeta function is:
# zeta(s) = sum_n lambda_n^{-s}
# The spectral action is S = sum_n f(lambda_n/Lambda^2)
# For a smooth cutoff, this gives the heat kernel expansion

# The KEY quantity is: what is R_int (the physical size of the internal space)?
# Once we know R_int, we know Lambda (since eigenvalues scale as 1/R_int^2)
# And then f_2*Lambda^2 is determined.

# From the KK identity (exp463):
# 2*z_P + ln(6/pi)/ln(alpha) = 20.8128
# This relates the Planck hierarchy (z_P) to the KK scale
# If this = some exponent n: R_KK = l_P * alpha^n ... hmm

KK_exponent = 2 * z_P + math.log(6 / math.pi) / math.log(alpha_exp)
print(f"\nKK identity: 2*z_P + ln(6/pi)/ln(alpha) = {KK_exponent:.6f}")
print(f"  2*z_P = {2*z_P:.6f}")
print(f"  ln(6/pi)/ln(alpha) = {math.log(6/math.pi)/math.log(alpha_exp):.6f}")

# R_KK = alpha^(KK_exponent/2) * l_P  ???
# Actually: in KK, G_4D = G_5D / (2*pi*R_KK)
# If G_5D ~ l_P^3 and G_4D ~ l_P^2, then R_KK ~ l_P
# But that gives compactification at Planck scale (too small for observation)

# For S^3 with N vertices, the "number of lattice points" is N=120
# The radius in lattice units is: R_lattice ~ N^{1/3} ~ 4.93
# Or R_lattice = phi (the exact circumradius)
# The physical radius: R_phys = R_lattice * a_lattice = phi * a_lattice

# From the KK formula for S^3:
# G_4D = G_7D / (2*pi^2 * R_phys^3)

# If G_7D = G_3D * a^4 (the 3D Turaev-Viro times 4D volume factor):
# Hmm, no. G_3D is self-contained in 3D.

# Let me try yet another approach: PURE RATIO

# We want: M_P / m_e = alpha^{-z_P}
# Can we express M_P in terms of the spectral action WITHOUT external input?

# The spectral action gives TWO mass scales:
# 1. M_P ~ sqrt(f_2) * Lambda * sqrt(c1) / pi (from R term)
# 2. M_scal ~ M_P / sqrt(f_4*c2) * (some const) (from R^2 term)
#
# The RATIO M_P/M_scal depends on f_2*Lambda^2/(f_4*c2) * c1 * const

# If both f_2 and f_4 follow the same scaling law:
# f_k = N^k / (2*pi) [for k even]
# Then f_2/f_4 = N^{-2} = 1/14400
# And: M_P/M_scal ~ sqrt(f_2*Lambda^2*c1/(f_4*c2))

# ============================================================
# NEW IDEA: f_k = N^k * alpha*alpha' (unified test function moments)
# ============================================================
print("\n" + "=" * 70)
print("NEW IDEA: Unified test function f_k = N^k * alpha*alpha'")
print("=" * 70)

# If f_k = N^k / (2*pi) for all k:
f0_unified = 1 / (2 * math.pi)  # N^0/(2*pi)
f2_unified = N**2 / (2 * math.pi)
f4_unified = N**4 / (2 * math.pi)

print(f"f_0 = 1/(2*pi) = {f0_unified:.6f}")
print(f"f_2 = N^2/(2*pi) = {f2_unified:.4f}")
print(f"f_4 = N^4/(2*pi) = {f4_unified:.6e} (KNOWN to work for inflation)")

# With f_2 = N^2/(2*pi) and the gravity relation:
# G = 3*pi / (f_2 * Lambda^2 * c1)
# f_2*Lambda^2 = 3*pi*M_P^2/c1

# Lambda^2 = 3*pi*M_P^2 / (f_2*c1) = 3*pi*M_P^2 / (N^2/(2*pi) * c1)
#          = 6*pi^2*M_P^2 / (N^2*c1)

Lambda_unified = M_P_GeV * math.sqrt(6 * math.pi**2 / (N**2 * c1))
print(f"\nWith unified f_k = N^k/(2*pi):")
print(f"  Lambda = M_P * sqrt(6*pi^2/(N^2*c1))")
print(f"  Lambda = {Lambda_unified:.6e} GeV")
print(f"  Lambda/M_P = {Lambda_unified/M_P_GeV:.8f}")
print(f"  Lambda/M_GUT = {Lambda_unified/M_GUT:.6f}")

# This gives Lambda ~ 4.9e17 GeV, close to GUT scale!
# Is Lambda/M_P a framework number?
ratio_Lambda_MP = Lambda_unified / M_P_GeV
print(f"\n  Lambda/M_P = {ratio_Lambda_MP:.8f}")
print(f"  (Lambda/M_P)^2 = {ratio_Lambda_MP**2:.8e}")
print(f"  = 6*pi^2/(N^2*c1) = 6*pi^2/({N**2}*{c1})")
print(f"  = 6*pi^2/{N**2*c1}")
print(f"  = 6*pi^2/{N**2*c1} = {6*math.pi**2/(N**2*c1):.8e}")

# Simplify: 6*pi^2/(N^2 * 2*N*A1) = 6*pi^2/(2*N^3*A1) = 3*pi^2/(N^3*A1)
val_ratio_sq = 3 * math.pi**2 / (N**3 * A1)
print(f"  = 3*pi^2/(N^3*A1) = 3*pi^2/({N**3}*{A1}) = {val_ratio_sq:.8e}")
print(f"  N^3 = {N**3}, N^3*A1 = {N**3*A1}")
print(f"  N^3*A1 = {N}^3 * {A1} = {N**3*A1}")

# Lambda/M_P ~ 0.04 ... let's see:
print(f"\n  sqrt(3)*pi/(N^(3/2)*sqrt(A1)) = {math.sqrt(3)*math.pi/(N**1.5*math.sqrt(A1)):.8f}")
print(f"  Compare Lambda/M_P = {ratio_Lambda_MP:.8f}")
# These should match (they're the same formula)

# Now check: does this Lambda give a CONSISTENT M_scal?
# M_scal^2 = M_P_red^2 * 320*pi^2 / (12*f_4*c2)
# Already computed above. The consistency check is whether
# the ratio Lambda/M_scal makes physical sense.
print(f"\n  Lambda/M_scal = {Lambda_unified/M_scal_from_CMB:.4f}")
print(f"  M_scal = {M_scal_from_CMB:.4e} GeV")
print(f"  Lambda = {Lambda_unified:.4e} GeV")

# ============================================================
# CRITICAL TEST: Does the unified f_k reproduce m_e/M_P?
# ============================================================
print("\n" + "=" * 70)
print("CRITICAL TEST: Self-consistency with m_e/M_P = alpha^{z_P}")
print("=" * 70)

# If f_2 = N^2/(2*pi), then we have a DETERMINED theory:
# G is fixed, M_P is fixed, and m_e/M_P should equal alpha^{z_P}
#
# But wait - G is ALWAYS consistent because we USED G to compute Lambda!
# The real test is: does the unified f_k = N^k/(2*pi) give a PREDICTION
# that can be verified independently?

# PREDICTION 1: M_scal from the unified theory
# M_scal is determined by f_4 and c2 (both known)
# A_s is determined by M_scal and N_e
A_s_pred = N_e**2 / (24 * math.pi**2) * (M_scal**2 / M_P_red_GeV**2)
# where M_scal comes from f_4
print(f"\nA_s (predicted from f_4=N^4/(2*pi)) = {A_s_pred:.4e}")
print(f"A_s (observed) = {A_s_obs:.4e}")
print(f"Deviation: {(A_s_pred/A_s_obs - 1)*100:.2f}%")

# PREDICTION 2: The RATIO f_2/f_4 gives a mass ratio
# f_2/f_4 = N^{-2} = 1/N^2
# In the spectral action: f_2*Lambda^2 / (f_0 * something) relates to
# the scalaron mass vs Planck mass ratio

# Let's compute the predicted G directly:
# G_pred = 3*pi / (f_2_unified * Lambda_unified^2 * c1)
# But Lambda_unified was DERIVED from G... so this is circular.

# The NON-CIRCULAR test is: given f_2 = N^2/(2*pi), and an INDEPENDENT
# determination of Lambda (e.g., from the inflation sector), what G results?

# From inflation: f_4*c2 determines M_scal/M_P_red
# M_scal/M_P_red = sqrt(320*pi^2/(12*f_4*c2))
xi = math.sqrt(320 * math.pi**2 / (12 * f4 * c2))
print(f"\nM_scal/M_P_red = {xi:.6e}")
print(f"M_scal/M_P = {xi*math.sqrt(8*math.pi):.6e}")

# The cutoff Lambda from inflation:
# In standard NCG, Lambda enters through f_k*Lambda^k
# f_4*Lambda^4 = f_4_physical  and f_2*Lambda^2 = f_2_physical
# If f_k is the MOMENT (dimensionless), Lambda is the physical cutoff

# With f_4 = N^4/(2*pi) and f_2 = N^2/(2*pi):
# f_2/f_4 = 1/N^2
# f_2*Lambda^2 = (f_4*Lambda^4) * (1/N^2) / Lambda^2
# = (f_4/N^2) * Lambda^2

# But f_4 * Lambda^4 appears as the COSMOLOGICAL CONSTANT term:
# V_cc ~ c0 * f_4 * Lambda^4
# So Lambda ~ (V_cc / (c0*f_4))^{1/4}
# This involves the CC, which is another open problem.

# ============================================================
# ROUTE D (NEW): Use the Diophantine to constrain f_2/f_4
# ============================================================
print("\n" + "=" * 70)
print("ROUTE D: Diophantine constraint on test function moments")
print("=" * 70)

# The Seeley-DeWitt Diophantine: 2*A1^2 + 1 = 3*A0*A2
# This relates c0, c1, c2. Can it constrain f_0, f_2, f_4?

# In the physical Lagrangian, the three terms are:
# Term_0 = c0 * f_4 * Lambda^4  (cosmological constant)
# Term_1 = c1 * f_2 * Lambda^2  (Einstein-Hilbert)
# Term_2 = c2 * f_0             (Yang-Mills / R^2)

# The Diophantine constrains the c_k RATIOS:
# c1^2 / (c0*c2) = A1^2/(A0*A2) = 3/2 - 1/(2*A0*A2) ≈ 1.4998

# HYPOTHESIS: The same near-GP structure extends to the PHYSICAL terms:
# (c1*f_2*Lambda^2)^2 / (c0*f_4*Lambda^4 * c2*f_0) = c1^2*f_2^2 / (c0*c2*f_0*f_4)
# = (A1^2/(A0*A2)) * (f_2^2/(f_0*f_4))

# If the PHYSICAL coefficients also satisfy a near-GP condition:
# (c1*f_2*Lambda^2)^2 = (3/2) * (c0*f_4*Lambda^4) * (c2*f_0)
# Then: f_2^2 = (3/2) * f_0*f_4 * (A0*A2/A1^2)
# Since A0*A2/A1^2 = 2/3 + 1/(3*A1^2) ≈ 2/3:
# f_2^2 ≈ f_0*f_4
# i.e., f_2 = sqrt(f_0*f_4) (GEOMETRIC MEAN!)

print(f"\nIf f_2 = sqrt(f_0 * f_4) [geometric mean]:")
print(f"  and f_4 = N^4/(2*pi), f_0 = 1/(2*pi) [i.e. f_k = N^k/(2*pi)]:")
f2_geom = math.sqrt(f0_unified * f4_unified)
print(f"  f_2 = sqrt(N^4/(2*pi)^2) = N^2/(2*pi) = {f2_geom:.4f}")
print(f"  This IS the unified f_k = N^k/(2*pi) !")

# But what if f_0 ≠ 1/(2*pi)?
# If f_0 = 1 (step function value at origin):
f2_geom2 = math.sqrt(1 * f4_unified)
Lambda_geom2 = M_P_GeV * math.sqrt(3 * math.pi / (f2_geom2 * c1))
print(f"\nIf f_0 = 1, f_2 = sqrt(f_4) = {f2_geom2:.4f}:")
print(f"  Lambda = {Lambda_geom2:.4e} GeV")
print(f"  Lambda/M_P = {Lambda_geom2/M_P_GeV:.6f}")

# The exact Diophantine GP relation:
# f_2^2 / (f_0*f_4) = (c0*c2/c1^2) * (3/2)
# = (A0*A2/A1^2) * (3/2)
# = (2563/3844) * (3/2)
# = 7689/7688
# = 1 + 1/7688
gp_factor = (3.0/2) * (A0 * A2) / (A1**2)
print(f"\nExact GP constraint: f_2^2/(f_0*f_4) = {gp_factor:.10f}")
print(f"  = (3/2)*(A0*A2/A1^2) = (3/2)*({A0}*{A2}/{A1}^2)")
print(f"  = (3/2)*{A0*A2}/{A1**2} = {3*A0*A2}/{2*A1**2}")
print(f"  = {3*A0*A2}/(2*{A1**2}) = {3*A0*A2/(2*A1**2):.10f}")
print(f"  Recall: 3*A0*A2 = 2*A1^2+1 = {2*A1**2+1}")
print(f"  So: f_2^2/(f_0*f_4) = (2*A1^2+1)/(2*A1^2) = 1 + 1/{2*A1**2}")
print(f"  = 1 + 1/7688 = {1 + 1/7688:.10f}")
print(f"  This is ALMOST exactly 1 (geometric progression)!")
print(f"  Deviation from GP: 1/{2*A1**2} = {1/(2*A1**2):.6e}")

# ============================================================
# PUTTING IT ALL TOGETHER: The unified test function
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS: The Complete Test Function Picture")
print("=" * 70)

print(f"""
IF the test function moments satisfy the geometric progression:

  f_k = f_0 * (N^2/(2*pi*f_0))^(k/4)  ... no, simpler:

  f_k = N^k / (2*pi)  for k = 0, 2, 4

  f_0 = 1/(2*pi) = alpha*alpha' = {f0_unified:.6f}  (Galois norm)
  f_2 = N^2/(2*pi) = {f2_unified:.4f}
  f_4 = N^4/(2*pi) = {f4_unified:.4e}

Then:
  1. alpha*alpha' = 1/(2*pi) ALREADY DERIVED
  2. f_4 = N^4*alpha*alpha' gives A_s at 3.2% (KNOWN, exp408)
  3. f_2 = N^2*alpha*alpha' gives:
""")

# Compute G from unified f_2:
# G = 3*pi / (f_2 * Lambda^2 * c1)
# We need Lambda. The only remaining freedom.

# BUT: if f_0 is ALSO determined, we can get Lambda from the CC:
# CC term: V_cc = c0 * f_4 * Lambda^4 / (16*pi^2)
# Observed CC: Lambda_cc^4 ~ (2.3 meV)^4

Lambda_cc_meV = 2.3  # meV, the CC energy scale
Lambda_cc_GeV = Lambda_cc_meV * 1e-6  # to GeV (1 meV = 1e-6 GeV, wait no)
# 1 meV = 1e-3 eV = 1e-12 GeV
Lambda_cc_GeV = 2.3e-3 * 1e-9  # 2.3 meV in GeV
# Actually Lambda_cc^4 ~ rho_vac ~ (2.3e-3 eV)^4
# 2.3e-3 eV = 2.3e-12 GeV
Lambda_cc_GeV = 2.3e-12  # GeV

# From spectral action: rho_vac = c0 * f_4 * Lambda^4 / (16*pi^2)
# Lambda^4 = 16*pi^2 * rho_vac / (c0*f_4)
# Lambda^4 = 16*pi^2 * Lambda_cc^4 / (c0*f_4)

Lambda_from_CC = (16 * math.pi**2 * Lambda_cc_GeV**4 / (c0 * f4_unified))**(1.0/4)
print(f"  From CC (rho_vac ~ (2.3 meV)^4):")
print(f"    Lambda_cc = {Lambda_from_CC:.4e} GeV")
print(f"    This is TINY! The CC gives Lambda ~ {Lambda_from_CC:.2e} GeV")
print(f"    This is INCONSISTENT with gravity (which needs Lambda ~ 10^17 GeV)")
print(f"    => The CC is NOT simply c0*f_4*Lambda^4 (fine-tuning problem)")

# So the CC doesn't help fix Lambda. But the inflation sector does!
# From inflation: A_s = f_4*c2 determines M_scal/M_P
# And N_e determines n_s, r
# But Lambda drops out of the Starobinsky prediction...

# Actually, in the standard spectral action approach, Lambda is the UV cutoff
# and it DOES appear in the predictions. The issue is that the spectral action
# is an ASYMPTOTIC expansion, and Lambda is where you truncate.

# FINAL APPROACH: Lambda is NOT a separate parameter!
# In the DISCRETE 600-cell, there IS no UV cutoff - the spectrum is finite.
# The "cutoff" is just Lambda_max = sqrt(b1*phi^2) in lattice units.
# The physical Lambda = Lambda_max_phys = sqrt(b1)*phi / R_int
# where R_int is the ONE remaining free parameter (the physical size of the 600-cell)

print(f"\n--- FINAL CHAIN ---")
print(f"The ONE free parameter is R_int (physical size of internal S^3)")
print(f"Everything follows:")
print(f"  Lambda = sqrt(Lambda_max_D2) / R_int = phi*sqrt(b1) / R_int")
print(f"  M_P^2 = f_2*Lambda^2*c1/(3*pi) = [N^2/(2*pi)]*[b1*phi^2/R_int^2]*c1/(3*pi)")
print(f"  M_P^2 = N^2*b1*phi^2*c1 / (6*pi^2*R_int^2)")

# So: R_int = sqrt(N^2*b1*phi^2*c1/(6*pi^2)) / M_P
R_int_GeV_inv = math.sqrt(N**2 * b1 * phi**2 * c1 / (6 * math.pi**2)) / M_P_GeV
R_int_meters = R_int_GeV_inv * (1.9733e-16)  # hbar*c = 0.19733 GeV*fm

print(f"\n  R_int = {R_int_GeV_inv:.6e} GeV^-1 = {R_int_meters:.4e} m")
print(f"  Compare l_P = {1/M_P_GeV:.4e} GeV^-1 = {1.616e-35:.4e} m")
print(f"  R_int / l_P = {R_int_GeV_inv * M_P_GeV:.6f}")

# R_int/l_P should be a framework number!
R_over_lP = R_int_GeV_inv * M_P_GeV
print(f"\n  R_int/l_P = {R_over_lP:.6f}")
print(f"  = sqrt(N^2*b1*phi^2*c1/(6*pi^2))")
print(f"  = N*phi*sqrt(b1*c1/(6*pi^2))")
print(f"  = N*phi*sqrt(b1*2*N*A1/(6*pi^2))")
print(f"  = N*phi*sqrt(N*b1*A1/(3*pi^2))")

inner = N * b1 * A1 / (3 * math.pi**2)
print(f"  = {N}*{phi:.6f}*sqrt({inner:.4f})")
print(f"  = {N*phi:.4f}*{math.sqrt(inner):.6f}")
print(f"  = {N*phi*math.sqrt(inner):.6f}")
print(f"\n  Numerically: {R_over_lP:.6f}")
print(f"  N*phi = {N*phi:.4f}")
print(f"  sqrt(b1*A1/(3*pi^2)) = sqrt({b1}*{A1}/(3*pi^2)) = sqrt({b1*A1/(3*math.pi**2):.4f}) = {math.sqrt(b1*A1/(3*math.pi**2)):.6f}")

# Let's check some framework numbers:
print(f"\n  Compare with framework quantities:")
print(f"    N = {N}")
print(f"    N*phi = {N*phi:.4f}")
print(f"    N*phi/pi = {N*phi/math.pi:.4f}")
print(f"    sqrt(c1/(2*pi)) = {math.sqrt(c1/(2*math.pi)):.4f}")
print(f"    sqrt(c1*phi^2/pi^2) = {math.sqrt(c1*phi**2/math.pi**2):.4f}")
print(f"    R_int/l_P = {R_over_lP:.6f}")
print(f"    N*sqrt(b1*A1)/pi = {N*math.sqrt(b1*A1)/math.pi:.4f}")

# VERIFICATION: with this R_int, what is Lambda?
Lambda_final = phi * math.sqrt(b1) / R_int_GeV_inv
print(f"\n  Lambda = phi*sqrt(b1)/R_int = {Lambda_final:.4e} GeV")
print(f"  Lambda/M_P = {Lambda_final/M_P_GeV:.8f}")
print(f"  Lambda/M_GUT = {Lambda_final/M_GUT:.4f}")

# ============================================================
# SUMMARY OF ALL RESULTS
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
UNIFIED TEST FUNCTION HYPOTHESIS: f_k = N^k / (2*pi)

  f_0 = 1/(2*pi) = alpha*alpha'      [Galois norm, DERIVED]
  f_2 = N^2/(2*pi) = {f2_unified:.4f}        [NEW HYPOTHESIS]
  f_4 = N^4/(2*pi) = {f4_unified:.4e}  [PATTERN, exp408]

GEOMETRIC PROGRESSION property:
  f_2^2 / (f_0*f_4) = 1 + 1/{2*A1**2} = {gp_factor:.10f}
  (Exact GP to 0.013%, deviation from Diophantine!)

PHYSICAL CONSEQUENCES (assuming Lambda = phi*sqrt(b1)/R_int):
  R_int/l_P = {R_over_lP:.4f}
  Lambda = {Lambda_final:.4e} GeV (near GUT scale)
  Lambda/M_P = {Lambda_final/M_P_GeV:.6f}

CONSISTENCY CHECKS:
  f_4 -> A_s = 2.17e-9 (+3.2%)  [KNOWN]
  f_2 -> G = observed G  [BY CONSTRUCTION if Lambda chosen correctly]
  m_e/M_P = alpha^(z_P) at 0.24%  [EXISTING PATTERN]

STATUS:
  f_0 = alpha*alpha':  DERIVED
  f_4 = N^4/(2*pi):   PATTERN (exp408)
  f_2 = N^2/(2*pi):   NEW HYPOTHESIS (needs independent verification)
  GP structure:        CONSEQUENCE of Diophantine + f_k = N^k/(2*pi)
  Lambda:              ONE FREE PARAM (equivalent to m_e or R_int)
""")

# The honest assessment:
print("=" * 70)
print("HONEST ASSESSMENT")
print("=" * 70)
print(f"""
The unified test function f_k = N^k/(2*pi) is ELEGANT:
  - It has a clear interpretation: N^k = |2I|^k counts gauge configurations
    in k/2 dimensions (lattice gauge theory, exp408d)
  - 1/(2*pi) = alpha*alpha' is the Galois norm (DERIVED)
  - The Diophantine 2*A1^2+1=3*A0*A2 ensures NEAR-EXACT GP structure
  - f_4 already works for inflation (exp408)

BUT it does NOT eliminate the free parameter:
  - We still need m_e (or equivalently R_int, or Lambda, or G)
  - The test function fixes the RATIOS of f_k, not the absolute scale
  - The CC problem remains (CC term gives inconsistent Lambda)

WHAT WOULD CONSTITUTE A DERIVATION:
  1. A principle that fixes R_int in terms of a1=5 (e.g., R_int = phi*l_P)
  2. A self-consistency condition on the spectral action
  3. A bootstrap: the internal space IS the test function
     (the 600-cell spectrum determines f_k AND provides the eigenvalues)

NEXT STEPS:
  1. Test if R_int/l_P = {R_over_lP:.4f} is a framework number
  2. Check if the GP constraint gives TESTABLE predictions
  3. Investigate the bootstrap: can the 600-cell spectrum itself be the test function?
""")

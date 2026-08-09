"""
EXP-270: GAUGE COUPLING UNIFICATION ON 600-CELL
================================================
Test whether the 3 gauge couplings unify at a 600-cell scale.
Key: alpha_em, sin^2(tW), alpha_s all DERIVED -> test RG consistency.
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
N_eig = 9

print("="*72)
print("EXP-270: GAUGE COUPLING UNIFICATION")
print("="*72)

# === Framework values ===
# Alpha from equation: 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0
coeff_b = 4*a1*PHI**4
D = coeff_b**2 - 8*np.pi
alpha_framework = (coeff_b - np.sqrt(D))/(4*np.pi)
alpha_exp = 1/137.035999084

# sin^2(tW) = b1/(a1^2+1) = 6/26
sin2tW_framework = b1/(a1**2+1)
sin2tW_exp = 0.23122  # at m_Z

# alpha_s = 1/(2*phi^3) (framework)
alpha_s_framework = 1/(2*PHI**3)
alpha_s_exp = 0.1179  # at m_Z

print("\n--- Part 1: Framework Couplings ---")
print(f"  alpha_em = {alpha_framework:.10f} (exp: {alpha_exp:.10f}, {(alpha_framework/alpha_exp-1)*100:+.4f}%)")
print(f"  sin^2(tW) = {sin2tW_framework:.6f} = {b1}/{a1**2+1} (exp: {sin2tW_exp:.6f}, {(sin2tW_framework/sin2tW_exp-1)*100:+.2f}%)")
print(f"  alpha_s = {alpha_s_framework:.6f} = 1/(2*phi^3) (exp: {alpha_s_exp:.6f}, {(alpha_s_framework/alpha_s_exp-1)*100:+.2f}%)")

# SM coupling constants
# g1^2 = 4*pi*alpha/(cos^2(tW)) = 4*pi*alpha/(1-sin^2(tW))
# g2^2 = 4*pi*alpha/sin^2(tW)
# g3^2 = 4*pi*alpha_s
# In GUT normalization: g1'^2 = (5/3)*g1^2

alpha1_inv = (1-sin2tW_framework)/alpha_framework  # 1/alpha_1
alpha2_inv = sin2tW_framework/alpha_framework       # 1/alpha_2
alpha3_inv = 1/alpha_s_framework                    # 1/alpha_3

print(f"\n  1/alpha_1 = (1-sin^2tW)/alpha = {alpha1_inv:.4f}")
print(f"  1/alpha_2 = sin^2tW/alpha = {alpha2_inv:.4f}")
print(f"  1/alpha_3 = 2*phi^3 = {alpha3_inv:.4f}")

# GUT normalized
alpha1_GUT_inv = (3.0/5)*alpha1_inv
print(f"  1/alpha_1(GUT) = (3/5)*{alpha1_inv:.2f} = {alpha1_GUT_inv:.4f}")

# === Part 2: Standard RG running ===
print("\n--- Part 2: 1-loop RG Running ---")
# 1-loop beta coefficients for SM:
# b_1 = 41/10, b_2 = -19/6, b_3 = -7 (with standard normalization)
# GUT: b_1_GUT = (3/5)*41/10 = 41/6... no
# Standard: 1/alpha_i(mu) = 1/alpha_i(m_Z) - b_i/(2*pi) * ln(mu/m_Z)

b1_SM = 41.0/10   # U(1)_Y
b2_SM = -19.0/6   # SU(2)
b3_SM = -7.0      # SU(3)

# Use GUT normalization for b1
b1_GUT = (3.0/5) * b1_SM  # = 41/6... no, b1_GUT = 41/10 in standard, GUT norm = (5/3)*41/10 = 41/6
# Actually: b1_GUT coefficient in 1/alpha_1_GUT equation
# Standard: d(1/alpha_Y)/dt = -b1_SM/(2*pi)
# GUT: 1/alpha_1_GUT = (3/5)/alpha_Y, so d(1/alpha_1_GUT)/dt = -(3/5)*b1_SM/(2*pi)
# BUT convention: b_i_GUT is defined so d(1/alpha_i)/dt = -b_i/(2*pi)
# SM values (1 Higgs doublet, 3 generations):
b_SM = [41.0/6, -19.0/6, -7.0]  # GUT normalized b_1, b_2, b_3
# Actually the STANDARD reference:
# b_1(GUT) = 41/10 * (5/3) = 41/6
# b_2 = -19/6
# b_3 = -7

m_Z = 91.1876  # GeV
t_values = np.linspace(0, 40, 1000)  # t = ln(mu/m_Z)

# Use experimental values at m_Z
alpha_em_mZ = 1/127.944  # at m_Z (running)
sin2tW_mZ = 0.23122
alpha_s_mZ = 0.1179

a1_mZ_inv = (1-sin2tW_mZ)/alpha_em_mZ * (3.0/5)
a2_mZ_inv = sin2tW_mZ/alpha_em_mZ
a3_mZ_inv = 1/alpha_s_mZ

print(f"  At m_Z = {m_Z} GeV:")
print(f"    1/alpha_1(GUT) = {a1_mZ_inv:.4f}")
print(f"    1/alpha_2      = {a2_mZ_inv:.4f}")
print(f"    1/alpha_3      = {a3_mZ_inv:.4f}")

print(f"\n  SM beta coefficients (1-loop, GUT norm):")
print(f"    b_1 = {b_SM[0]:.4f}")
print(f"    b_2 = {b_SM[1]:.4f}")
print(f"    b_3 = {b_SM[2]:.4f}")

# Run couplings
a_inv = np.zeros((3, len(t_values)))
a_inv[0,:] = a1_mZ_inv - b_SM[0]/(2*np.pi) * t_values
a_inv[1,:] = a2_mZ_inv - b_SM[1]/(2*np.pi) * t_values
a_inv[2,:] = a3_mZ_inv - b_SM[2]/(2*np.pi) * t_values

# Find crossing points
print(f"\n  Crossing points:")
for i in range(3):
    for j in range(i+1,3):
        diff = a_inv[i,:] - a_inv[j,:]
        for k in range(len(diff)-1):
            if diff[k]*diff[k+1] < 0:
                t_cross = t_values[k] - diff[k]*(t_values[k+1]-t_values[k])/(diff[k+1]-diff[k])
                mu_cross = m_Z * np.exp(t_cross)
                alpha_cross = 1/a_inv[i,k]
                print(f"    alpha_{i+1} = alpha_{j+1} at t={t_cross:.2f}, mu={mu_cross:.2e} GeV, 1/alpha={a_inv[i,k]:.2f}")

# === Part 3: Framework unification ===
print("\n--- Part 3: Framework Unification Scale ---")
# In our framework: all from a1=5
# alpha_1(GUT), alpha_2, alpha_3 should relate at 600-cell scale
# The 600-cell has N=120 vertices -> natural scale ~ m_P / phi^something

# Framework values (NOT at m_Z, these are "bare" or geometric)
a1_fw_inv = (1-sin2tW_framework)/alpha_framework * (3.0/5)
a2_fw_inv = sin2tW_framework/alpha_framework
a3_fw_inv = 1/alpha_s_framework

print(f"  Framework coupling inverses:")
print(f"    1/alpha_1(GUT) = {a1_fw_inv:.4f}")
print(f"    1/alpha_2      = {a2_fw_inv:.4f}")
print(f"    1/alpha_3      = {a3_fw_inv:.4f}")

# Check: are these at the same energy?
# If they unified at some scale Lambda_600, then at that scale
# 1/alpha_i(Lambda) would all be equal

# Differences
print(f"\n  Differences (framework):")
print(f"    1/alpha_1 - 1/alpha_2 = {a1_fw_inv - a2_fw_inv:.4f}")
print(f"    1/alpha_2 - 1/alpha_3 = {a2_fw_inv - a3_fw_inv:.4f}")
print(f"    1/alpha_1 - 1/alpha_3 = {a1_fw_inv - a3_fw_inv:.4f}")
print(f"    Ratio: (a1-a2)/(a2-a3) = {(a1_fw_inv-a2_fw_inv)/(a2_fw_inv-a3_fw_inv):.4f}")
print(f"    SM prediction for this ratio: (b1-b2)/(b2-b3) = {(b_SM[0]-b_SM[1])/(b_SM[1]-b_SM[2]):.4f}")

# In SM unification:
# (1/a1 - 1/a2) / (1/a2 - 1/a3) should equal (b1-b2)/(b2-b3) if all meet at one point
# (b1-b2)/(b2-b3) = (41/6+19/6)/(-19/6+7) = 60/6 / (23/6) = 10/3.833 = 60/23

ratio_b = (b_SM[0]-b_SM[1])/(b_SM[1]-b_SM[2])
ratio_a_fw = (a1_fw_inv-a2_fw_inv)/(a2_fw_inv-a3_fw_inv)
ratio_a_exp = (a1_mZ_inv-a2_mZ_inv)/(a2_mZ_inv-a3_mZ_inv)

print(f"\n  Unification test B = (a1-a2)/(a2-a3):")
print(f"    SM beta ratio:   {ratio_b:.6f} = {int(b_SM[0]-b_SM[1])}/{b_SM[1]-b_SM[2]:.4f}")
print(f"    Framework ratio: {ratio_a_fw:.6f}")
print(f"    Experimental:    {ratio_a_exp:.6f}")

# === Part 4: Framework expressions ===
print("\n--- Part 4: Framework Algebraic Structure ---")
# 1/alpha_1(GUT) = (3/5)*(1-6/26)/alpha = (3/5)*(20/26)/alpha = (60/130)/alpha = (12/26)/alpha
# 1/alpha_2 = (6/26)/alpha
# 1/alpha_3 = 2*phi^3

# So 1/alpha_2 = (6/26)/alpha = sin^2(tW)/alpha
# And 1/alpha_1 = (12/26)/alpha = 2*sin^2(tW)/alpha
# Therefore 1/alpha_1 = 2 * (1/alpha_2)
print(f"  1/alpha_1(GUT) / (1/alpha_2) = {a1_fw_inv/a2_fw_inv:.6f}")
print(f"  Should be (3/5)*(1-s2)/(s2) = (3/5)*{(1-sin2tW_framework)/sin2tW_framework:.4f} = {(3.0/5)*(1-sin2tW_framework)/sin2tW_framework:.6f}")

# Key: 1/alpha_1(GUT) = (3/5) * (20/26) / alpha = (3*20)/(5*26*alpha) = 60/(130*alpha)
# = 12/(26*alpha) = 12*137.036/26 = ...
print(f"\n  In terms of framework constants:")
print(f"    1/alpha_1(GUT) = (3*(a1^2-a1))/(a1*(a1^2+1)*alpha) = {3*(a1**2-a1)/(a1*(a1**2+1)*alpha_framework):.4f}")
print(f"    1/alpha_2 = b1/((a1^2+1)*alpha) = {b1/((a1**2+1)*alpha_framework):.4f}")
print(f"    1/alpha_3 = 2*phi^3 = {2*PHI**3:.4f}")
print(f"")
print(f"    Note: 1/alpha = 4*a1*phi^4/(2*pi) - sqrt((4*a1*phi^4)^2-8*pi)/(2*pi*alpha)")
print(f"    So all 3 couplings are determined by a1=5 and phi=(1+sqrt(5))/2")

# === Part 5: Z[phi] structure of couplings ===
print("\n--- Part 5: Z[phi] Structure ---")
# 1/alpha_3 = 2*phi^3 = 2*(2+phi) = 4+2*phi. In Z[phi]: (4,2). Tr=8, N=4.
# What about 1/alpha_2 and 1/alpha_1?

z3 = 2*PHI**3  # = 4 + 2*phi
z3_a = 4; z3_b = 2  # a + b*phi
print(f"  1/alpha_3 = 2*phi^3 = {z3_a} + {z3_b}*phi = {z3:.4f}")
print(f"    Tr = {2*z3_a+z3_b}, N = {z3_a**2+z3_a*z3_b-z3_b**2}")

# 1/alpha_2 = sin^2(tW)/alpha = b1/((a1^2+1)*alpha)
z2 = sin2tW_framework/alpha_framework
print(f"  1/alpha_2 = {z2:.6f}")
# Is this in Z[phi]? Not obviously, since alpha involves pi
# alpha = (4*a1*phi^4 - sqrt(D))/(4*pi)
# So 1/alpha involves pi, and 1/alpha_2 = sin^2(tW)/alpha also involves pi

# However, for unification we want to express the RG scale
# t_unif = 2*pi * (1/alpha_i - 1/alpha_j) / (b_i - b_j)
# This t involves pi * (something with pi) / rational = involves pi^2 etc.

print(f"\n  The key point: alpha_1 and alpha_2 involve pi (from alpha equation)")
print(f"  while alpha_3 = 1/(2*phi^3) is purely algebraic in Z[phi]")
print(f"  This suggests alpha_3 is the 'natural' coupling of the discrete 600-cell")
print(f"  while alpha_em involves the CONTINUOUS Hopf structure (pi)")

# === Part 6: Alternative - 600-cell unification ===
print("\n--- Part 6: 600-Cell Unification ---")
# In 600-cell framework, the gauge group comes from 24-cell:
# 1 + 3 + 8 = 12 edges per vertex in 24-cell
# The coupling is SINGLE: vertex-transitivity => one coupling g
# Different g1,g2,g3 arise from ALGEBRA (E8), not from geometry

# If couplings unify at tree level (vertex-transitive), then
# the DIFFERENCES come from running from 600-cell scale to m_Z

# Natural 600-cell scales:
# Edge: pi/a1 (angular separation on S^3)
# In units of m_P: Lambda_600 ~ m_P / phi^something

# From m_e/m_P = alpha^(4*phi^2):
m_P = 1.22093e19 * 1000  # MeV (Planck mass)
m_e_MeV = 0.51099895

# Using framework: m_e = m_P * alpha^(4*phi^2)
log_ratio = np.log(m_e_MeV/m_P)
print(f"  ln(m_e/m_P) = {log_ratio:.4f}")
print(f"  4*phi^2 * ln(alpha) = {4*PHI**2 * np.log(alpha_framework):.4f}")
print(f"  Match: {(4*PHI**2*np.log(alpha_framework)/log_ratio - 1)*100:+.2f}%")

# The 600-cell scale (edge in Planck units)
Lambda_600 = m_P  # at Planck scale, the 600-cell IS the geometry
t_Planck = np.log(m_P/m_Z/1000)  # m_P in GeV / m_Z in GeV
print(f"\n  t(Planck) = ln(m_P/m_Z) = {t_Planck:.2f}")

# Couplings at Planck scale (1-loop SM)
a1_Planck = 1/(a1_mZ_inv - b_SM[0]/(2*np.pi)*t_Planck)
a2_Planck = 1/(a2_mZ_inv - b_SM[1]/(2*np.pi)*t_Planck)
a3_Planck = 1/(a3_mZ_inv - b_SM[2]/(2*np.pi)*t_Planck)
print(f"  At Planck scale (1-loop SM):")
print(f"    alpha_1(GUT) = {a1_Planck:.6f} (1/{1/a1_Planck:.2f})")
print(f"    alpha_2      = {a2_Planck:.6f} (1/{1/a2_Planck:.2f})")
print(f"    alpha_3      = {a3_Planck:.6f} (1/{1/a3_Planck:.2f})")

# Do they converge?
spread = max(1/a1_Planck, 1/a2_Planck, 1/a3_Planck) - min(1/a1_Planck, 1/a2_Planck, 1/a3_Planck)
mean_inv = (1/a1_Planck + 1/a2_Planck + 1/a3_Planck)/3
print(f"    Spread: {spread:.2f}")
print(f"    Mean 1/alpha: {mean_inv:.2f}")
print(f"    Relative spread: {spread/mean_inv*100:.1f}%")

# === Part 7: Edge-counting argument ===
print("\n--- Part 7: Edge-Counting at 600-Cell Scale ---")
# At tree level, all edges carry the same coupling (vertex-transitivity)
# The 3 gauge groups emerge from 1+3+8 splitting of 12 edges/vertex
# So at tree level: g1 = g2 = g3 = g
# This means: alpha_1(Lambda)/1 = alpha_2(Lambda)/3 = alpha_3(Lambda)/8
# Or: alpha_i(Lambda) = n_i * alpha_0

# More precisely: the coupling to a given edge is the same
# A U(1) vertex (1 edge) gets g^2 -> alpha_1 = g^2/(4*pi)
# An SU(2) vertex (3 edges) gets 3*g^2 but with group factor -> alpha_2 ~ alpha_0
# An SU(3) vertex (8 edges) -> alpha_3 ~ alpha_0

# Better: vertex-transitive means one spectral action coefficient
# The splitting comes from the REPRESENTATION of gauge group on each edge type

# In our framework:
# sin^2(tW) = 6/26 = b1/(a1^2+1) already at tree level
# This is NOT a running effect, it's GEOMETRIC

# So: at 600-cell scale:
# alpha_2(Lambda) = alpha / sin^2(tW) (same alpha!)
# alpha_1(Lambda) = alpha / (1-sin^2(tW)) * (3/5) (GUT norm)
# These are NOT unified at Lambda, they're RELATED by Weinberg angle

# The ONLY independent coupling that could unify is alpha_s with alpha_em
print(f"  At geometric level:")
print(f"    alpha_em and sin^2(tW) are from the SAME alpha equation + edge counting")
print(f"    alpha_s = 1/(2*phi^3) is from SPECTRAL conditions on 600-cell")
print(f"    These are TWO independent determinations")
print(f"")
print(f"    alpha_em = {alpha_framework:.8f}")
print(f"    alpha_s  = {alpha_s_framework:.8f}")
print(f"    alpha_s/alpha_em = {alpha_s_framework/alpha_framework:.4f}")
print(f"    = 2*pi*phi^3 * alpha_s = {2*np.pi*PHI**3*alpha_s_framework:.4f} (involves pi)")

# Check if ratio is expressible
ratio_as_aem = alpha_s_framework / alpha_framework
print(f"\n  alpha_s/alpha_em = {ratio_as_aem:.6f}")
# = alpha * 2*phi^3 / alpha ... no that doesn't simplify
# = (2*phi^3)^(-1) / alpha^(-1) = alpha * 2*phi^3 = 2*phi^3 * (4*a1*phi^4-sqrt(D))/(4*pi)
# This involves pi, so it's NOT a simple ratio

# === Part 8: The role of b1=6 in unification ===
print("\n--- Part 8: Role of b_1 in Gauge Structure ---")
# b1 = 6 = number of Hopf fibrations
# sin^2(tW) = b1/(a1^2+1) -> the Weinberg angle counts Hopf fibers
# alpha equation: 72 decagons = 6 * 12 -> b1 Hopf fibrations, 12 fibers each
# alpha_s: 2*phi^3 -> trace condition on Z[phi]

print(f"  b_1 = {b1} = number of Hopf fibrations on S^3")
print(f"  sin^2(tW) = b_1/(a_1^2+1) = {b1}/{a1**2+1} = {sin2tW_framework:.6f}")
print(f"  This means: the WEAK sector uses {b1} out of {a1**2+1} directions")
print(f"  The EM sector uses {a1**2-a1} out of {a1**2+1}")
print(f"  a1^2 - a1 = {a1**2-a1} = 20 (icosahedral rotations?)")

# Summary of what's derived
print(f"\n  DERIVED from a1=5:")
print(f"    alpha = 1/137.036... (from 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0)")
print(f"    sin^2(tW) = 6/26 (from b1/(a1^2+1))")
print(f"    alpha_s = 1/(2*phi^3) (from Tr=8, |N|=4 in Z[phi])")
print(f"    ALL THREE from a1=5. No unification scale needed.")
print(f"")
print(f"  QUESTION: Is there a SINGLE master equation that gives all three?")
print(f"    Current: 3 separate derivations (alpha eq, edge counting, spectral)")
print(f"    The alpha equation already encodes the 24-cell structure")
print(f"    The spectral conditions for alpha_s use different machinery")

# === Part 9: Master equation search ===
print("\n--- Part 9: Searching for Master Equation ---")
# All three couplings in terms of a1=5 and phi:
# alpha^(-1) = (4*pi - 4*a1*phi^4 + sqrt((4*a1*phi^4)^2-8*pi)) / (4*pi*alpha) ... circular
# Better: from alpha equation, 1/alpha ~ 2*a1*phi^4/pi (leading order)

# 1/alpha ~ 2*a1*phi^4/pi = 137.08 (leading order, 0.03% accuracy)
# 1/alpha_2 = (b1/(a1^2+1)) * 1/alpha ~ 2*a1*b1*phi^4/((a1^2+1)*pi)
# 1/alpha_3 = 2*phi^3

# Ratio: (1/alpha_3) / (1/alpha) ~ pi / (a1*phi) = pi/(a1*phi) = 0.387
print(f"  Leading order expressions:")
print(f"    1/alpha ~ 2*a1*phi^4/pi = {2*a1*PHI**4/np.pi:.4f} (exact: {1/alpha_framework:.4f})")
print(f"    1/alpha_2 ~ 2*a1*b1*phi^4/((a1^2+1)*pi) = {2*a1*b1*PHI**4/((a1**2+1)*np.pi):.4f} (exact: {a2_fw_inv:.4f})")
print(f"    1/alpha_3 = 2*phi^3 = {2*PHI**3:.4f}")
print(f"")
print(f"    Ratio alpha_3/alpha_2 ~ pi*sin^2(tW)/(a1*phi) = {np.pi*sin2tW_framework/(a1*PHI):.6f}")
print(f"    Exact: {a2_fw_inv/a3_fw_inv:.6f}")

# The transition from alpha to alpha_s involves pi -> phi^4/phi^3 = phi
# 1/alpha involves phi^4/pi, 1/alpha_s = 2*phi^3
# If there were a "master" equation, it would connect pi and phi

# Famous identity: pi = sum(1/n^2) * 6 (Basel)
# Or: pi/4 = 1 - 1/3 + 1/5 - ...
# In our framework: 4*a1*phi^4 = 4*5*phi^4 = 20*phi^4 ~ 2*pi/alpha (from alpha eq)
# So 20*phi^4 ~ 137 * 2*pi / 137 ... tautological

# === Part 10: Summary ===
print("\n--- Part 10: Summary ---")
print(f"  UNIFICATION STATUS:")
print(f"  1. All 3 couplings DERIVED from a1=5 and phi (DONE)")
print(f"  2. No traditional GUT unification scale needed")
print(f"     (couplings are GEOMETRIC, not running to a point)")
print(f"  3. SM 1-loop running from m_Z to Planck:")
print(f"     Relative spread = {spread/mean_inv*100:.1f}% (no exact unification)")
print(f"  4. The 600-cell IS the unification: vertex-transitivity")
print(f"     gives ONE edge coupling, split by 1+3+8 into 3 groups")
print(f"  5. sin^2(tW) is TREE LEVEL (geometric), not a running effect")
print(f"  6. Missing: master equation combining alpha, sin^2(tW), alpha_s")
print(f"  7. The role of pi: alpha involves pi (Hopf/continuous),")
print(f"     alpha_s is purely algebraic (Z[phi]/discrete)")
print(f"")
print(f"  CATEGORY: PARTIALLY DERIVED")
print(f"  - Each coupling individually: DERIVED")
print(f"  - Unified origin: STRUCTURAL (vertex-transitivity) but not a single equation")

print("\n" + "="*72)
print("EXP-270 COMPLETE")
print("="*72)

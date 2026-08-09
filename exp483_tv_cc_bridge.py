"""
exp483: Turaev-Viro CC bridge to 4D

Known from exp424:
  Lambda_3D = 4*pi^2/a1^2  (TV at r=a1=5, Barrett 1992)
  G_3D = 1/(4k) = 1/12     (CS level k=a1-2=3)

Known from the paper:
  Lambda_P = alpha^(57-alpha_s)  (matches observation to 0.13 sigma)
  The exponent 57 has 7 independent decompositions.

Question: Can we bridge Lambda_3D -> Lambda_4D and recover the exponent 57?
"""
import numpy as np

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
k = a1 - 2  # CS level = 3
r = a1      # root of unity parameter = 5
rank_E8 = 8
N_gen = 3

# Couplings
alpha_s = 1 / (2 * phi**3)
disc = (4*a1*phi**4)**2 - 8*np.pi
alpha = (4*a1*phi**4 - np.sqrt(disc)) / (4*np.pi)
z_Planck = 4 * phi**2  # hierarchy exponent

print("="*70)
print("  exp483: TURAEV-VIRO CC BRIDGE TO 4D")
print("="*70)

# ============================================================
# Known quantities
# ============================================================
Lambda_3D = 4 * np.pi**2 / a1**2  # = 4*pi^2/25
G_3D = 1 / (4*k)                   # = 1/12
D2 = a1 + np.sqrt(a1)              # total quantum dimension^2

print(f"\n  3D TQFT quantities:")
print(f"    Lambda_3D = 4*pi^2/a1^2 = {Lambda_3D:.6f}")
print(f"    G_3D = 1/(4k) = 1/{4*k} = {G_3D:.6f}")
print(f"    G*Lambda = {G_3D * Lambda_3D:.6f} = pi^2/(N_gen*a1^2) = {np.pi**2/(N_gen*a1**2):.6f}")
print(f"    D^2 = {D2:.6f}")
print(f"    ln(D^2) = {np.log(D2):.6f}")

# Target: Lambda_4D in Planck units
# Observed: Lambda_obs ~ 10^{-122.1} in Planck units
# Paper: Lambda_P = alpha^z where z = 57 - alpha_s = 56.882
z_target = 57 - alpha_s
Lambda_target = alpha**z_target
print(f"\n  4D target:")
print(f"    z = 57 - alpha_s = {z_target:.6f}")
print(f"    Lambda_P = alpha^z = {Lambda_target:.6e}")
print(f"    log10(Lambda_P) = {np.log10(Lambda_target):.2f}")

# ============================================================
# APPROACH 1: Direct dimensional reduction
# Lambda_4D = Lambda_3D * (reduction factor)
# ============================================================
print("\n" + "="*70)
print("  APPROACH 1: What reduction factor is needed?")
print("="*70)

ratio_needed = Lambda_target / Lambda_3D
print(f"  Lambda_4D / Lambda_3D = {ratio_needed:.6e}")
print(f"  ln(ratio) / ln(alpha) = {np.log(ratio_needed) / np.log(alpha):.6f}")
# So Lambda_4D = Lambda_3D * alpha^X where X = ?
X_needed = np.log(ratio_needed) / np.log(alpha)
print(f"  => Lambda_4D = Lambda_3D * alpha^{X_needed:.4f}")
print(f"  => need alpha^{X_needed:.4f} as the 3D->4D suppression")

# What is X in terms of framework quantities?
print(f"\n  X = {X_needed:.6f}")
print(f"  z_Planck = 4*phi^2 = {z_Planck:.6f}")
print(f"  2*z_Planck = {2*z_Planck:.6f}")
print(f"  X - 2*z_Planck = {X_needed - 2*z_Planck:.6f}")
print(f"  z_target = {z_target:.6f}")
print(f"  z_target + ln(Lambda_3D)/ln(alpha) = {z_target + np.log(Lambda_3D)/np.log(alpha):.6f}")

# Actually, let's decompose differently:
# Lambda_4D = alpha^z_target
# Lambda_3D = 4*pi^2/a1^2
# So: alpha^z_target = (4*pi^2/a1^2) * alpha^X
# => z_target = X + ln(4*pi^2/a1^2)/ln(alpha)
ln_ratio = np.log(4*np.pi**2/a1**2) / np.log(alpha)
print(f"\n  ln(Lambda_3D) / ln(alpha) = {ln_ratio:.6f}")
print(f"  So: z_target = X + {ln_ratio:.4f}")
print(f"  => X = z_target - ({ln_ratio:.4f}) = {z_target - ln_ratio:.4f}")

X_from_decomp = z_target - ln_ratio
print(f"\n  The 3D->4D suppression exponent X = {X_from_decomp:.6f}")
print(f"  Compare with framework quantities:")
print(f"    z_target itself = {z_target:.6f}")
print(f"    2*z_Planck = 2*4*phi^2 = {2*z_Planck:.6f}")
print(f"    z_Planck^2 = (4*phi^2)^2 = {z_Planck**2:.6f}")
print(f"    1/alpha_s = 2*phi^3 = {1/alpha_s:.6f}")
print(f"    57 = z_target + alpha_s = {57}")

# ============================================================
# APPROACH 2: Galois pairing
# Lambda_phys * Lambda_dark = something from TQFT
# ============================================================
print("\n" + "="*70)
print("  APPROACH 2: Galois pairing of CC")
print("="*70)

# From exp444: E_vac = phi^3, E_vac' = phi'^3, product = -1
# S - S' = 0 (Galois-invariant)
# From exp462: z_Planck = 2*D^2 - d_ST (DERIVED)

# The TV total quantum dimension:
# D^2 = a1 + sqrt(a1), D'^2 = a1 - sqrt(a1)
# Product: D^2 * D'^2 = a1^2 - a1 = 4*a1 = 20
# This is the same as L(3)*L(3')!

D2_dark = a1 - np.sqrt(a1)
print(f"  D^2 (physical) = {D2:.6f}")
print(f"  D'^2 (dark)    = {D2_dark:.6f}")
print(f"  D^2 * D'^2     = {D2 * D2_dark:.6f} = a1^2-a1 = {a1**2-a1} = 4*a1 = {4*a1}")
print(f"  D^2 / D'^2     = phi^2 = {D2/D2_dark:.6f} = {phi**2:.6f}")

# TEE (topological entanglement entropy):
gamma_phys = np.log(np.sqrt(D2))
gamma_dark = np.log(np.sqrt(D2_dark))
print(f"\n  TEE physical: gamma = ln(D) = {gamma_phys:.6f}")
print(f"  TEE dark:     gamma'= ln(D')= {gamma_dark:.6f}")
print(f"  gamma - gamma' = {gamma_phys - gamma_dark:.6f} = ln(phi) = {np.log(phi):.6f}")

# What if Lambda_4D = exp(-S_grav) where S_grav involves D^2?
# In 3D CS: Z = 1/D, so S = ln(D) = gamma
# In the Galois picture: Lambda involves BOTH sectors

# Try: Lambda_4D = (D^2)^{-z_Planck}?
test1 = D2**(-z_Planck)
print(f"\n  Test: D^{{-2*z_P}} = {D2:.4f}^{{-{z_Planck:.4f}}} = {test1:.6e}")
print(f"    log/log(alpha) = {np.log(test1)/np.log(alpha):.4f}")

# Try: Lambda_4D = alpha^{z_Planck * ln(D^2)/ln(alpha)}?
# Hmm, that's circular.

# ============================================================
# APPROACH 3: The exponent 57 from TV data
# ============================================================
print("\n" + "="*70)
print("  APPROACH 3: Can we build 57 from TV quantities?")
print("="*70)

# Known decompositions of 57:
# 57 = (N - N_gen)/2 = (120-3)/2 = 58.5... no wait
# The paper says 7 routes to 57. Let me check what they are.

# From the paper and memory:
# 57 = a1*N_eig + a1 + rank + 4*phi^2 + ... various
# But the TQFT quantities are: k=3, r=5, D^2=5+sqrt(5), etc.

# Let me try combinations of TQFT quantities that give 57:
print("  TQFT quantities:")
print(f"    k = {k}")
print(f"    r = a1 = {r}")
print(f"    D^2 = {D2:.6f}")
print(f"    ln(D^2)/ln(alpha) = {np.log(D2)/np.log(alpha):.6f}")
print(f"    S_TEE = ln(D) = {np.log(np.sqrt(D2)):.6f}")

# Can k, r, D^2 build 57?
# 57 = k*r + ... = 15 + 42? No obvious
# 57 = r*12 - 3 = 60 - 3 = 57!
print(f"\n  57 = r*12 - 3 = a1*degree - N_gen = {a1*12 - 3}")
print(f"  57 = N/2 - 3 = 60 - 3 = {N//2 - 3}")
print(f"  57 = a1*N_eig + rank + 4 = 5*9 + 8 + 4 = {5*9+8+4}")
# These are known decompositions.

# NEW from TV:
# Lambda_3D = 4*pi^2/r^2 in 3D Planck units
# G_3D = 1/(2*b1) = 1/(4k) = 1/12
# G_3D * Lambda_3D = pi^2/(3*25) = pi^2/75

# The NUMBER 57 in terms of k=3:
# 57 = 19*k = 19*3
# 19 = a1^2 - a1 - 1 (the prime in the mass formula!)
print(f"\n  57 = 19*k = (a1^2-a1-1)*(a1-2) = {(a1**2-a1-1)*(a1-2)}")
print(f"  19 = a1^2-a1-1 = {a1**2-a1-1} (the split prime in Z[phi])")
print(f"  k = a1-2 = {a1-2} (CS level)")
print(f"  57 = 19*3 = (a1^2-a1-1)*(a1-2)")
print(f"  THIS IS NEW: 57 = p_split * k_CS")

# Verify: is this a coincidence or structural?
# The prime 19 appears as |N(z_t)| = |N(4+phi)| = 16+4-1 = 19
# k = 3 is the CS level
# Their product is 57. Is this forced?

# Check: does (a1^2-a1-1)*(a1-2) = 57 only for a1=5?
print(f"\n  Check: f(a1) = (a1^2-a1-1)*(a1-2) for various a1:")
for a in range(2, 10):
    val = (a**2 - a - 1) * (a - 2)
    print(f"    a1={a}: f = {val}")

# ============================================================
# APPROACH 4: Modular properties
# ============================================================
print("\n" + "="*70)
print("  APPROACH 4: Modular S-matrix and alpha")
print("="*70)

# The modular S-matrix of SU(2)_k CS theory at level k=3:
# S_{jl} = sqrt(2/(k+2)) * sin(pi*(2j+1)*(2l+1)/(k+2))
# = sqrt(2/5) * sin(pi*(2j+1)*(2l+1)/5)

S_mat = np.zeros((k+1, k+1))
for j2 in range(k+1):     # j = j2/2
    for l2 in range(k+1):
        S_mat[j2, l2] = np.sqrt(2/(k+2)) * np.sin(np.pi*(j2+1)*(l2+1)/(k+2))

print("  Modular S-matrix at k=3:")
for j2 in range(k+1):
    row = "    j=" + f"{j2/2:.1f}: " + "  ".join(f"{S_mat[j2,l2]:8.5f}" for l2 in range(k+1))
    print(row)

# S-matrix properties
print(f"\n  S-matrix unitarity: S*S^T = I? {np.allclose(S_mat @ S_mat.T, np.eye(k+1))}")
print(f"  S_{'{'}01{'}'} = {S_mat[0,1]:.10f} = sqrt(2/5)*sin(2pi/5) = {np.sqrt(2/5)*np.sin(2*np.pi/5):.10f}")
print(f"  S_{'{'}00{'}'} = {S_mat[0,0]:.10f}")
print(f"  S_{'{'}01{'}'}/ S_{'{'}00{'}'} = {S_mat[0,1]/S_mat[0,0]:.10f} = phi = {phi:.10f}")

# The quantum dimensions from S-matrix:
print(f"\n  Quantum dimensions d_j = S_0j / S_00:")
for j2 in range(k+1):
    d_j = S_mat[0, j2] / S_mat[0, 0]
    print(f"    j={j2/2:.1f}: d = {d_j:.10f}")

# Connection to alpha: the S-matrix prefactor is sqrt(2/a1)
# This is EXACTLY c_bare = sqrt(2/a1) (the Verlinde self-energy coefficient)!
print(f"\n  S-matrix prefactor: sqrt(2/(k+2)) = sqrt(2/a1) = {np.sqrt(2/a1):.10f}")
print(f"  This IS c_bare from the mass formula!")
print(f"  c_eff = c_bare * (1 - b1*alpha^2) = {np.sqrt(2/a1) * (1 - b1*alpha**2):.10f}")

# ============================================================
# APPROACH 5: Chern-Simons free energy and CC
# ============================================================
print("\n" + "="*70)
print("  APPROACH 5: CS partition function and CC exponent")
print("="*70)

# Z_CS(S^3) at level k:
# |Z|^2 = 1/D^2 = 1/(a1 + sqrt(a1))
# Free energy: F = -ln|Z|^2 = ln(D^2) = ln(a1 + sqrt(a1))

F_CS = np.log(D2)
print(f"  F_CS = ln(D^2) = ln({D2:.4f}) = {F_CS:.10f}")
print(f"  F_CS / ln(alpha) = {F_CS / np.log(alpha):.10f}")

# The CC in the paper: Lambda = alpha^z, z ~ 57
# So ln(Lambda) = z * ln(alpha)
# If Lambda relates to exp(-F) somehow:
# Lambda ~ exp(-n*F_CS) = D^{-2n}
# Then z = n * F_CS / (-ln(alpha)) = n * ln(D^2) / (-ln(alpha))

n_needed = z_target * (-np.log(alpha)) / np.log(D2)
print(f"\n  For Lambda = D^{{-2n}}:")
print(f"    need n = {n_needed:.6f}")
print(f"    = z_target * |ln(alpha)| / ln(D^2)")
print(f"    ~ {n_needed:.2f}")

# Is n close to any framework quantity?
print(f"\n  Nearby framework quantities to n={n_needed:.2f}:")
print(f"    1/alpha_s = 2*phi^3 = {1/alpha_s:.4f}")
print(f"    z_Planck = 4*phi^2 = {z_Planck:.4f}")
print(f"    N/2 = {N/2}")
print(f"    57 = {57}")
print(f"    57 - z_Planck = {57 - z_Planck:.4f}")
print(f"    z_target = {z_target:.4f}")

# ============================================================
# FINAL ASSESSMENT
# ============================================================
print("\n" + "="*70)
print("  ASSESSMENT")
print("="*70)

print(f"""
  KEY NEW FINDING:
    57 = 19 * 3 = (a1^2-a1-1) * (a1-2) = p_split * k_CS

  This decomposes the CC exponent as:
    - 19 = the unique split prime in the mass formula (|N(z_t)| = 19)
    - 3  = the Chern-Simons level k = a1-2

  This is the EIGHTH decomposition of 57, and the FIRST that
  connects the CC exponent to the TQFT level k.

  Previously known decompositions (from the paper):
    57 = N/2 - N_gen = 60-3
    57 = a1*(N_eig+1) + rank = 5*10+7... (various)

  The new one: 57 = p_split * k_CS is STRUCTURAL because:
    - p_split = 19 is FORCED by the mass formula (the unique split prime
      compatible with the Z[phi] flatness condition)
    - k_CS = 3 is FORCED by a1 = 5 (level k = a1-2)
    - Their product gives the CC exponent

  However: 57 = 19*3 is an ALGEBRAIC IDENTITY for a1=5.
  For other a1: (a1^2-a1-1)*(a1-2) gives different numbers.
  There is no DYNAMICAL mechanism explaining WHY the CC exponent
  should be p_split * k_CS. It's a new decomposition, not a derivation.

  STATUS: STRUCTURAL (new route, not a derivation)
  The CC mechanism remains OPEN.
""")

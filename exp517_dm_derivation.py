"""
EXP-517: Deriving Omega_DM/Omega_b = 7 - phi from TQFT + Spectral Action
==========================================================================

KNOWN FORMULAS (equivalent):
  (A) 7 - phi = b1 - 1/phi = b1 - |sigma(phi)|           = 5.381966...
  (B) (N_broken - Delta_D^2) / d_ST = (26 - 2*sqrt(5))/4 = 5.381966...

KEY IDENTITY DISCOVERED:
  A_phys = sum N_{ij}^k d_i d_j d_k = D^4  (general MTC identity)
  A_dark_signed = D'^4  (Galois conjugate)
  A_dark_unsigned = 2*(D'^2 + 2)  (using |d'| instead of d')

STRATEGY:
  1. Verify all identities algebraically and numerically
  2. Show A_phys = D^4 is a THEOREM (unitarity of S-matrix)
  3. Explore the "quantum anomaly" interpretation:
     Omega = (classical_modes - quantum_correction) / d_ST
  4. Compute Seeley-DeWitt coefficients for C^2/2I orbifold
  5. Find which spectral action order produces the formula
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen, degree

print("=" * 70)
print("EXP-517: DERIVING DM ABUNDANCE FROM TQFT + SPECTRAL ACTION")
print("=" * 70)

# ============================================================
# PART 1: IDENTITIES
# ============================================================
print(f"\n{'='*70}")
print("PART 1: FUNDAMENTAL IDENTITIES")
print(f"{'='*70}")

k = 3
n_reps = k + 1  # 4 reps: j=0, 1/2, 1, 3/2

# Quantum dimensions
d_q = np.array([1.0, PHI, PHI, 1.0])
d_dark_unsigned = np.array([1.0, 1/PHI, 1/PHI, 1.0])
d_dark_signed = np.array([1.0, -1/PHI, -1/PHI, 1.0])

D2 = np.sum(d_q**2)
D2_dark = np.sum(d_dark_unsigned**2)

# S-matrix
S = np.zeros((n_reps, n_reps))
for i in range(n_reps):
    for j in range(n_reps):
        S[i,j] = np.sqrt(2.0/(k+2)) * np.sin((2*(i/2.0)+1)*(2*(j/2.0)+1)*np.pi/(k+2))

# Fusion coefficients
fusion = np.zeros((n_reps, n_reps, n_reps))
for i in range(n_reps):
    for j in range(n_reps):
        for kk in range(n_reps):
            val = sum(S[i,l]*S[j,l]*S[kk,l]/S[0,l] for l in range(n_reps))
            fusion[i,j,kk] = round(val)

# Compute amplitudes
A_phys = sum(fusion[i,j,kk]*d_q[i]*d_q[j]*d_q[kk]
             for i in range(n_reps) for j in range(n_reps) for kk in range(n_reps)
             if fusion[i,j,kk] > 0)

A_dark_u = sum(fusion[i,j,kk]*d_dark_unsigned[i]*d_dark_unsigned[j]*d_dark_unsigned[kk]
               for i in range(n_reps) for j in range(n_reps) for kk in range(n_reps)
               if fusion[i,j,kk] > 0)

A_dark_s = sum(fusion[i,j,kk]*d_dark_signed[i]*d_dark_signed[j]*d_dark_signed[kk]
               for i in range(n_reps) for j in range(n_reps) for kk in range(n_reps)
               if fusion[i,j,kk] > 0)

D4 = D2**2
D4_dark = D2_dark**2

print(f"""
  IDENTITY 1: sum N_{{ij}}^k d_i d_j d_k = D^4
    A_phys            = {A_phys:.6f}
    D^4 = (5+sqrt5)^2 = {D4:.6f}
    Match: {abs(A_phys - D4) < 1e-10}

  PROOF: By unitarity of S-matrix,
    sum_i S_{{il}} S_{{0i}} = delta_{{l0}}  (orthogonality)
    sum_{{ijk}} N_{{ij}}^k d_i d_j d_k
    = sum_{{ijk}} [sum_l S_il S_jl S_kl / S_0l] * [S_0i/S_00] * [S_0j/S_00] * [S_0k/S_00]
    = 1/S_00^3 * sum_l (1/S_0l) * [sum_i S_il S_0i]^2 * [sum_k S_kl S_0k]
    = 1/S_00^3 * (1/S_00) * 1 * 1 * 1  (only l=0 survives)
    = 1/S_00^4 = D^4.  QED.

  IDENTITY 2: sum N_{{ij}}^k sigma(d_i) sigma(d_j) sigma(d_k) = sigma(D^4) = D'^4
    (Since N_{{ij}}^k are integers, Galois-invariant)
    A_dark_signed     = {A_dark_s:.6f}
    D'^4 = (5-sqrt5)^2 = {D4_dark:.6f}
    Match: {abs(A_dark_s - D4_dark) < 1e-10}

  UNSIGNED dark amplitude (using |d'| instead of d'):
    A_dark_unsigned   = {A_dark_u:.6f}
    2*(D'^2 + 2)      = {2*(D2_dark+2):.6f}
    Match: {abs(A_dark_u - 2*(D2_dark+2)) < 1e-10}
""")

# ============================================================
# PART 2: THREE RATIOS
# ============================================================
print(f"{'='*70}")
print("PART 2: THE THREE RATIOS")
print(f"{'='*70}")

target = 7 - PHI

r_signed = A_phys / A_dark_s
r_unsigned = A_phys / A_dark_u
r_mixed = (A_phys - A_dark_s) / (A_phys + A_dark_s)

print(f"""
  Target: 7-phi = {target:.6f}

  Ratio 1 (signed):   D^4/D'^4 = phi^4 = {r_signed:.6f}  (err: {(r_signed/target-1)*100:+.2f}%)
  Ratio 2 (unsigned): D^4/(2(D'^2+2)) = {r_unsigned:.6f}  (err: {(r_unsigned/target-1)*100:+.2f}%)
  Ratio 3 (asymmetry): (D^4-D'^4)/(D^4+D'^4) = {r_mixed:.6f}

  phi^4 = {PHI**4:.6f}    (27.4% too high)
  5.496                     (2.1% too high)

  NEITHER ratio gives 7-phi exactly.
  The formula 7-phi requires DIFFERENT ingredients.
""")

# ============================================================
# PART 3: THE QUANTUM ANOMALY FORMULA
# ============================================================
print(f"{'='*70}")
print("PART 3: THE QUANTUM ANOMALY FORMULA")
print(f"{'='*70}")

# The 9 irreps of 2I (binary icosahedral group)
irrep_dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]  # McKay: extended E8 nodes

# Eigenvalues of the 600-cell Laplacian (one per irrep)
L = [0,                  # rho_0, dim 1
     12 - 6*PHI,         # rho_1, dim 2, GALOIS-BROKEN
     12 - 2*SQRT5,       # rho_2, dim 3, GALOIS-BROKEN
     9,                  # rho_3, dim 4
     12,                 # rho_4, dim 5
     14,                 # rho_5, dim 6
     12 + 2*SQRT5,       # rho_6, dim 3, GALOIS-BROKEN (partner of rho_2)
     15,                 # rho_7, dim 4 -- WAIT
     6 + 6*PHI]          # rho_8, dim 2, GALOIS-BROKEN (partner of rho_1)

# Correction: rho_7 has dim 2 and rho_8... let me use the standard ordering
# McKay quiver of 2I -> extended E8:
# Node:  0  1  2  3  4  5  6  7  8 (where 8 is the extending node)
# dim:   1  2  3  4  5  6  4  2  3
# Wait - for the binary icosahedral group 2I of order 120:
# The 9 irreps have dimensions: 1, 2, 3, 4, 5, 6, 4, 2, 3
# The conjugate pairs under Galois (sqrt(5) -> -sqrt(5)):
# rho_1 (dim 2) <-> rho_7 (dim 2)
# rho_2 (dim 3) <-> rho_8 (dim 3)

# Multiplicities = dim^2
mults = [d**2 for d in irrep_dims]

# Galois classification
galois_broken_irreps = [1, 2, 7, 8]  # indices
galois_fixed_irreps = [0, 3, 4, 5, 6]

N_broken = sum(mults[i] for i in galois_broken_irreps)
N_fixed = sum(mults[i] for i in galois_fixed_irreps)

Delta_D2 = D2 - D2_dark  # = 2*sqrt(5)
d_ST = 4

Omega_formula_B = (N_broken - Delta_D2) / d_ST
Omega_formula_A = b1 - 1/PHI

print(f"""
  2I GROUP STRUCTURE:
    9 irreps, dims: {irrep_dims}
    Multiplicities (dim^2): {mults}
    Total: {sum(mults)} = N = {N}

  GALOIS CLASSIFICATION:
    Broken irreps: rho_1(dim 2), rho_2(dim 3), rho_7(dim 2), rho_8(dim 3)
    Broken modes: 4 + 9 + 4 + 9 = {N_broken} = a1^2 + 1 = {a1**2+1}
    Fixed modes: {N_fixed}

  TQFT DATA:
    D^2 = {D2:.6f} = 5 + sqrt(5)
    D'^2 = {D2_dark:.6f} = 5 - sqrt(5)
    Delta_D^2 = D^2 - D'^2 = {Delta_D2:.6f} = 2*sqrt(5) = {2*SQRT5:.6f}

  THE FORMULA:
    Omega = (N_broken - Delta_D^2) / d_ST
          = ({N_broken} - 2*sqrt(5)) / {d_ST}
          = (26 - {2*SQRT5:.4f}) / 4
          = {Omega_formula_B:.6f}
          = 7 - phi = {target:.6f}
    Match: {abs(Omega_formula_B - target) < 1e-10}

  EQUIVALENT FORM:
    Omega = b1 - 1/phi = {b1} - {1/PHI:.6f} = {Omega_formula_A:.6f}
    Match: {abs(Omega_formula_A - target) < 1e-10}
""")

# Verify algebraic equivalence
print(f"  ALGEBRAIC EQUIVALENCE:")
print(f"    (a1^2+1 - 2*sqrt(a1))/4")
print(f"    = ({a1**2+1} - 2*{SQRT5:.4f}) / 4")
print(f"    = (26 - 2sqrt5) / 4")
print(f"    = 13/2 - sqrt5/2")
print(f"    = 13/2 - (2phi-1)/2")
print(f"    = (13 - 2phi + 1)/2")
print(f"    = (14 - 2phi)/2")
print(f"    = 7 - phi  QED")
print(f"")
print(f"    b1 - 1/phi")
print(f"    = 6 - (phi-1)")
print(f"    = 7 - phi  QED")
print(f"")
print(f"    Both representations are EXACT IDENTITIES.")

# ============================================================
# PART 4: PHYSICAL INTERPRETATION - QUANTUM ANOMALY
# ============================================================
print(f"\n{'='*70}")
print("PART 4: QUANTUM ANOMALY INTERPRETATION")
print(f"{'='*70}")

print(f"""
  THE STRUCTURE: Omega = (classical - quantum) / normalization

  CLASSICAL TERM: N_broken = 26
    = Number of Galois-broken spectral modes of the 600-cell Laplacian
    = Seeley-DeWitt a_0 coefficient for the Galois-broken sector
    = Same 26 that appears in sin^2(tW) = 6/26

  QUANTUM TERM: Delta_D^2 = 2*sqrt(5) = {2*SQRT5:.4f}
    = D^2(physical) - D'^2(dark) = total quantum dimension excess
    = TQFT correction to mode counting (irrational!)
    = Measures the "quantum dimension imbalance" between sectors

  NORMALIZATION: d_ST = 4
    = Number of spacetime dimensions
    = Per-dimension averaging of effective degrees of freedom
""")

# ============================================================
# PART 5: SEELEY-DEWITT COEFFICIENTS ON C^2/2I
# ============================================================
print(f"{'='*70}")
print("PART 5: SEELEY-DEWITT ANALYSIS OF C^2/2I")
print(f"{'='*70}")

# Heat kernel on C^2/2I:
# K(t) = (1/|2I|) * sum_{g in 2I} K_0(g, t)
# where K_0(g,t) is the flat-space heat kernel twisted by g.
#
# This decomposes over irreps:
# K(t) = sum_rho dim(rho)^2 * K_rho(t)
#
# The Seeley-DeWitt expansion:
# K(t) ~ (4*pi*t)^{-d/2} * sum_n a_n * t^n
#
# For our orbifold (d=4, real dimension of C^2):
# K(t) ~ (4*pi*t)^{-2} * (a_0 + a_1*t + a_2*t^2 + ...)

# a_0 = total number of modes = 120
# a_0(broken) = 26
# a_0(fixed) = 94

# The Galois asymmetry at each order:
# At order a_0: N_broken = 26 (mode counting)
# At order a_1: involves first moment of eigenvalues
# At order a_2: involves second moment (curvature term)

print(f"\n  HEAT KERNEL MOMENTS OF GALOIS-BROKEN SECTOR:")
print(f"  {'Order':>8} {'Description':>20} {'Value':>12} {'Contribution':>15}")

# For each Galois pair, compute moments
# Pair 1: (rho_1, rho_7): dims (2,2), eigenvalues (L1, L8)
L1 = 12 - 6*PHI    # = 2.292
L8 = 6 + 6*PHI     # = 15.708... wait, let me recalculate

# Actually: eigenvalue of rho_7 (dim 2). Let me recheck.
# The 9 eigenvalues of the 600-cell graph Laplacian are:
# lambda_j = degree - sum of characters, where degree = 12
# For rho_1 (dim 2): characters involve phi
# lambda_1 = 12 - 2*phi - 2*phi - 2*1 = 12 - 4phi - 2? No...

# The eigenvalues are:
# lambda = 12 - chi(C_5)/dim(rho) ... actually it's more complex.
# From the data: 0, 12-6phi, 12-2sqrt5, 9, 12, 14, 12+2sqrt5, 15, 6+6phi
# These are in order of the irreps.

# Galois pairs:
# (L_1 = 12-6phi, L_8 = 6+6phi) with mult 4 each
# (L_2 = 12-2sqrt5, L_6 = 12+2sqrt5) with mult 9 each

L_1 = 12 - 6*PHI       # rho_1, mult 4
L_8 = 6 + 6*PHI        # rho_8 (partner), mult 4  -- actually rho_7 partner
L_2 = 12 - 2*SQRT5     # rho_2, mult 9
L_6 = 12 + 2*SQRT5     # rho_6 (partner), mult 9  -- actually rho_8 partner

print(f"\n  Galois pairs:")
print(f"    Pair A: L_1 = 12-6phi = {L_1:.4f}, L_8 = 6+6phi = {L_8:.4f}, mult = 4")
print(f"    Pair B: L_2 = 12-2s5 = {L_2:.4f}, L_6 = 12+2s5 = {L_6:.4f}, mult = 9")

# Moments: M_n = sum_broken mult * lambda^n
# Galois DIFFERENCE: delta_M_n = sum_pairs mult*(L_phys^n - L_dark^n)
# where "phys" is the smaller eigenvalue and "dark" is the larger

# Order 0: delta_M_0 (mode counting asymmetry)
# Both partners have same multiplicity, so delta_M_0 = 0
# The TOTAL broken modes = 26, but the ASYMMETRY at order 0 is 0.
delta_M0 = 4*(1-1) + 9*(1-1)
print(f"\n  delta_M_0 (mode asymmetry): {delta_M0}")

# Order 1: delta_M_1 = sum mult*(L_low - L_high)
delta_M1 = 4*(L_1 - L_8) + 9*(L_2 - L_6)
print(f"  delta_M_1 (eigenvalue asymmetry): {delta_M1:.4f}")
# L_1 - L_8 = (12-6phi) - (6+6phi) = 6-12phi = -6*(2phi-1) = -6*sqrt5
# L_2 - L_6 = (12-2s5) - (12+2s5) = -4*sqrt5
# delta_M1 = 4*(-6sqrt5) + 9*(-4sqrt5) = -24sqrt5 - 36sqrt5 = -60sqrt5
print(f"  = -60*sqrt(5) = {-60*SQRT5:.4f}")
print(f"  = -N/2 * Delta_D^2 * sqrt(a1) = {-N/2*Delta_D2*SQRT5:.4f}?")
print(f"  Actually: -60*sqrt5. And N = 120, so -N*sqrt5/2 = {-N*SQRT5/2:.4f}")
print(f"  Match: {abs(delta_M1 - (-N*SQRT5/2)) < 1e-10}")

# Order 2: delta_M_2 = sum mult*(L_low^2 - L_high^2)
delta_M2 = 4*(L_1**2 - L_8**2) + 9*(L_2**2 - L_6**2)
print(f"  delta_M_2 (curvature asymmetry): {delta_M2:.4f}")
# L_1^2 - L_8^2 = (L_1-L_8)(L_1+L_8) = (-6sqrt5)(18-6phi+6+6phi)... wait
# L_1+L_8 = 12-6phi+6+6phi = 18
# L_1^2-L_8^2 = (L_1-L_8)(L_1+L_8) = -6sqrt5 * 18 = -108sqrt5
# L_2+L_6 = (12-2s5)+(12+2s5) = 24
# L_2^2-L_6^2 = (L_2-L_6)(L_2+L_6) = -4sqrt5 * 24 = -96sqrt5
# delta_M2 = 4*(-108sqrt5) + 9*(-96sqrt5) = -432sqrt5 - 864sqrt5 = -1296sqrt5
print(f"  = 4*(-108sqrt5) + 9*(-96sqrt5) = {delta_M2:.4f}")
print(f"  = -1296*sqrt5 = {-1296*SQRT5:.4f}")
print(f"  = -(6^4)*sqrt5 = {-(6**4)*SQRT5:.4f}")
print(f"  Match: {abs(delta_M2 - (-(6**4)*SQRT5)) < 1e-10}")

# ============================================================
# PART 6: WHAT COMBINATION GIVES 7-phi?
# ============================================================
print(f"\n{'='*70}")
print("PART 6: SEARCHING FOR 7-phi IN SPECTRAL DATA")
print(f"{'='*70}")

# Key insight: 7-phi involves BOTH 26 (total broken modes) AND 2sqrt5 (TQFT).
# The heat kernel moments give only sqrt5 multiples (asymmetries within pairs).
# The NUMBER 26 comes from the TOTAL count, not the asymmetry.

# So the formula Omega = (26 - 2sqrt5)/4 combines:
# - TOTAL Galois-broken modes (not just the asymmetry)
# - TQFT quantum dimension difference
# - Spacetime dimensions

# These are THREE DIFFERENT structures:
# 26 : 2I representation theory (McKay/E8)
# 2sqrt5 : SU(2)_3 TQFT quantum dimensions
# 4 : spacetime

print(f"""
  THE THREE INGREDIENTS:
    26 = sum dim(rho_broken)^2       [2I representation theory]
    2sqrt5 = D^2 - D'^2             [SU(2)_3 TQFT]
    4 = d_ST                        [spacetime]

  These come from THREE levels of the theory:
    26 from the GROUP (2I, binary icosahedral, order 120)
    2sqrt5 from the TQFT (SU(2)_3, level k=3, k+2=a1=5)
    4 from SPACETIME (orbifold C^2/2I has real dim 4)

  All three are controlled by a1 = 5:
    26 = a1^2 + 1 = 5^2 + 1
    2sqrt5 = 2*sqrt(a1)
    4 = a1 - 1
""")

# Check: d_ST = a1 - 1?
print(f"  d_ST = {d_ST} = a1-1 = {a1-1}: {d_ST == a1-1}")

# So: Omega = (a1^2+1 - 2*sqrt(a1)) / (a1-1)
Omega_from_a1 = (a1**2 + 1 - 2*np.sqrt(a1)) / (a1 - 1)
print(f"\n  Omega(a1) = (a1^2+1-2sqrt(a1))/(a1-1)")
print(f"           = ({a1**2+1} - {2*np.sqrt(a1):.4f}) / {a1-1}")
print(f"           = {Omega_from_a1:.6f}")
print(f"  7 - phi  = {target:.6f}")
print(f"  Match: {abs(Omega_from_a1 - target) < 1e-10}")

# ============================================================
# PART 7: CAN WE DERIVE THIS FROM f(a1)?
# ============================================================
print(f"\n{'='*70}")
print("PART 7: UNIQUENESS FROM a1 = 5")
print(f"{'='*70}")

# The formula: Omega(a1) = (a1^2 + 1 - 2*sqrt(a1)) / (a1-1)
# Simplify: numerator = (sqrt(a1) - 1)^2 + a1^2 - a1
#           = (sqrt(a1)-1)^2 + a1*(a1-1)
# So: Omega = (sqrt(a1)-1)^2/(a1-1) + a1

# Or: Omega = a1 + (sqrt(a1)-1)^2/(a1-1)
check = a1 + (np.sqrt(a1)-1)**2/(a1-1)
print(f"  Omega = a1 + (sqrt(a1)-1)^2/(a1-1) = {check:.6f}")

# Alternative factorization?
# (a1^2+1-2sqrt(a1))/(a1-1) = ?
# Let x = sqrt(a1). Then: (x^4+1-2x)/(x^2-1)
# = (x^4-2x+1)/(x^2-1)
# = (x^4-1-2x+2)/(x^2-1)
# = (x^4-1)/(x^2-1) + (-2x+2)/(x^2-1)
# = (x^2+1) + (-2(x-1))/((x-1)(x+1))
# = x^2+1 - 2/(x+1)
# = a1+1 - 2/(sqrt(a1)+1)
# = b1 - 2/(sqrt(a1)+1)

check2 = b1 - 2/(np.sqrt(a1)+1)
print(f"  Omega = b1 - 2/(sqrt(a1)+1) = {b1} - {2/(np.sqrt(a1)+1):.6f} = {check2:.6f}")
print(f"  Match: {abs(check2 - target) < 1e-10}")

# And 2/(sqrt(5)+1) = 2/(sqrt(5)+1) * (sqrt(5)-1)/(sqrt(5)-1) = 2(sqrt(5)-1)/4 = (sqrt(5)-1)/2 = 1/phi
print(f"\n  KEY: 2/(sqrt(a1)+1) = 2/(sqrt(5)+1) = (sqrt(5)-1)/2 = 1/phi = {1/PHI:.6f}")
print(f"  So: Omega = b1 - 1/phi")
print(f"  This is the CLEANEST form of the formula!")

# ============================================================
# PART 8: WHY b1 - 1/phi?
# ============================================================
print(f"\n{'='*70}")
print("PART 8: PHYSICAL MEANING OF b1 - 1/phi")
print(f"{'='*70}")

# b1 = 6 = a1+1 = k+3 = ||N_{1/2}||^2_F (Frobenius norm of fusion matrix)
# 1/phi = |sigma(phi)| = dark quantum dimension of fundamental rep
# 1/phi = 2/(sqrt(a1)+1) (rationalized form)

# In the Verlinde algebra:
# b1 = ||N_{1/2}||^2 counts the total fusion channels of the fundamental rep
# 1/phi = |sigma(d_{1/2})| is the dark sector quantum dimension

# So: Omega = (fusion channels of fundamental) - (dark quantum dimension)
#           = (topological integer) - (Galois-conjugate quantum number)

print(f"""
  b1 = {b1} = ||N_{{1/2}}||^2_F = Frobenius norm of fundamental fusion matrix
     = Number of nonzero fusion channels in N_{{1/2}}
     = Total fusion "capacity" of the fundamental representation

  1/phi = {1/PHI:.6f} = |sigma(phi)| = dark quantum dimension
     = Galois conjugate of the physical quantum dimension phi
     = "Effective dark particle weight"

  Omega = b1 - 1/phi = {target:.6f}
     = "Total fusion capacity minus dark weight"
     = "Net effective Galois-dark degrees of freedom"

  INTERPRETATION:
    Each fusion channel contributes 1 unit of "interaction capacity".
    The dark sector's quantum dimension (1/phi) acts as a "quantum tax".
    The DM abundance = fusion_channels - quantum_tax.
""")

# ============================================================
# PART 9: THE FORMULA AS AN ALGEBRAIC INTEGER
# ============================================================
print(f"{'='*70}")
print("PART 9: NUMBER-THEORETIC PROPERTIES")
print(f"{'='*70}")

# 7-phi lives in Q(sqrt(5))
# 7-phi = (13-sqrt(5))/2
# Its Galois conjugate: 7-phi' = (13+sqrt(5))/2
# Norm: N(7-phi) = (169-5)/4 = 164/4 = 41 (PRIME)
# Trace: Tr(7-phi) = 13

sigma_Omega = 7 - PHI_CONJ  # = 7+(1/phi) = 7+phi-1 = 6+phi
Omega_norm = (7-PHI)*(7-PHI_CONJ)  # wait, PHI_CONJ = -1/PHI
# 7-phi' = 7 - (1-sqrt5)/2 = (14-1+sqrt5)/2 = (13+sqrt5)/2

sigma_target = (13+SQRT5)/2
norm_Omega = target * sigma_target
trace_Omega = target + sigma_target

print(f"""
  Omega = 7-phi = (13-sqrt5)/2 in Q(sqrt(5))
  sigma(Omega) = (13+sqrt5)/2 = {sigma_target:.6f}

  Algebraic properties:
    Norm = Omega * sigma(Omega) = {norm_Omega:.6f} = 164/4 = 41
    Trace = Omega + sigma(Omega) = {trace_Omega:.6f} = 13
    Discriminant = Trace^2 - 4*Norm = {trace_Omega**2 - 4*norm_Omega:.0f} = 5

  41 is PRIME.
  13 = F_7 (Fibonacci number!)
  5 = a1 (fundamental integer)

  Minimal polynomial: x^2 - 13x + 41 = 0
  Verification: {target**2 - 13*target + 41:.10f} = 0
""")

# ============================================================
# PART 10: NORMALIZED VERLINDE OBSERVABLES - SYSTEMATIC SEARCH
# ============================================================
print(f"{'='*70}")
print("PART 10: SYSTEMATIC TQFT OBSERVABLE SEARCH")
print(f"{'='*70}")

# Compute all useful TQFT quantities
# Traces of fusion matrices
Tr_N = [np.trace(fusion[j].astype(int)) for j in range(n_reps)]

# Frobenius norms
Frob_N = [np.sum(fusion[j]**2) for j in range(n_reps)]

# d-weighted traces
d_Tr_phys = sum(d_q[j] * Tr_N[j] for j in range(n_reps))
d_Tr_dark = sum(d_dark_unsigned[j] * Tr_N[j] for j in range(n_reps))

# d^2-weighted Frobenius
d2_Frob_phys = sum(d_q[j]**2 * Frob_N[j] for j in range(n_reps))
d2_Frob_dark = sum(d_dark_unsigned[j]**2 * Frob_N[j] for j in range(n_reps))

# Moments of quantum dimensions
sum_d = sum(d_q)
sum_d_dark = sum(d_dark_unsigned)
sum_d3 = sum(d_q**3)
sum_d3_dark = sum(d_dark_unsigned**3)
sum_d4 = sum(d_q**4)
sum_d4_dark = sum(d_dark_unsigned**4)

# Von Neumann entropy
p_phys = d_q**2 / D2
S_vN_phys = -sum(p_phys * np.log(p_phys))
p_dark = d_dark_unsigned**2 / D2_dark
S_vN_dark = -sum(p_dark * np.log(p_dark))

# Topological entanglement entropy
S_topo_phys = np.log(np.sqrt(D2))
S_topo_dark = np.log(np.sqrt(D2_dark))

# Rényi entropy (order 2)
S_R2_phys = -np.log(sum(p_phys**2))
S_R2_dark = -np.log(sum(p_dark**2))

# Conformal weights
h = np.array([0, 3/20, 2/5, 3/4])

# Thermal partition function at various beta
def Z_thermal(beta, d_arr, D2_val):
    return sum(d_arr[j]**2 * np.exp(-2*np.pi*beta*h[j]) for j in range(n_reps))

# Collect ALL expressions
expressions = {}

# Basic ratios
expressions['D^2/D\'^2 = phi^2'] = D2/D2_dark
expressions['D^4/D\'^4 = phi^4'] = D4/D4_dark
expressions['A_phys/A_dark_unsigned'] = A_phys/A_dark_u
expressions['A_phys/A_dark_signed'] = A_phys/A_dark_s

# Sum/diff of quantum dims
expressions['sum d'] = sum_d
expressions['sum d\''] = sum_d_dark
expressions['sum d - sum d\''] = sum_d - sum_d_dark
expressions['D^2 - D\'^2'] = D2 - D2_dark
expressions['D^2 + D\'^2'] = D2 + D2_dark
expressions['D^2 * D\'^2'] = D2 * D2_dark
expressions['sum d / sum d\''] = sum_d / sum_d_dark
expressions['sum d^3 / sum d\'^3'] = sum_d3 / sum_d3_dark

# Verlinde data
expressions['sum Frob(N)'] = sum(Frob_N)
expressions['d_Tr_phys'] = d_Tr_phys
expressions['d_Tr_dark'] = d_Tr_dark
expressions['d_Tr_phys/d_Tr_dark'] = d_Tr_phys/d_Tr_dark

# Entropy differences
expressions['S_topo_phys/S_topo_dark'] = S_topo_phys/S_topo_dark
expressions['exp(S_topo_phys-S_topo_dark)'] = np.exp(S_topo_phys-S_topo_dark)
expressions['S_vN_phys/S_vN_dark'] = S_vN_phys/S_vN_dark
expressions['exp(S_vN_phys/S_vN_dark)'] = np.exp(S_vN_phys/S_vN_dark)

# The target formula components
expressions['b1 - 1/phi'] = b1 - 1/PHI
expressions['(26-2sqrt5)/4'] = (26-2*SQRT5)/4
expressions['(a1^2+1-2sqrt(a1))/(a1-1)'] = (a1**2+1-2*np.sqrt(a1))/(a1-1)

# Mixed TQFT-group expressions
expressions['N_broken/d_ST'] = N_broken/d_ST
expressions['N_broken/d_ST - sqrt(a1)/2'] = N_broken/d_ST - SQRT5/2
expressions['N_broken/(D^2-D\'^2)'] = N_broken/(D2-D2_dark)
expressions['(N_broken-1)/a1'] = (N_broken-1)/a1
expressions['N_broken*1/phi/d_ST'] = N_broken/(PHI*d_ST)
expressions['a1 + phi^(-2)'] = a1 + 1/PHI**2

# Fusion-weighted quantities
expressions['sum_j d_j Tr(N_j)'] = d_Tr_phys
expressions['sum_j d\'_j Tr(N_j)'] = d_Tr_dark
expressions['(D^2+2phi)/d_ST'] = (D2+2*PHI)/d_ST

# Thermal partition function ratios
for beta_name, beta_val in [('1', 1.0), ('1/2pi', 1/(2*np.pi)), ('1/phi', 1/PHI),
                              ('phi', PHI), ('1/a1', 1/a1), ('1/12', 1/12.0)]:
    Z_p = Z_thermal(beta_val, d_q, D2)
    Z_d = Z_thermal(beta_val, d_dark_unsigned, D2_dark)
    if Z_d > 1e-15:
        expressions[f'Z_phys/Z_dark(beta={beta_name})'] = Z_p/Z_d

# b1, N_gen combinations
expressions['b1 - phi + phi^2'] = b1 - PHI + PHI**2
expressions['b1 * (1 - 1/(b1*phi))'] = b1*(1-1/(b1*PHI))
expressions['N_gen * phi + 1/(N_gen*phi)'] = N_gen*PHI+1/(N_gen*PHI)
expressions['2*phi^2'] = 2*PHI**2
expressions['phi^3 + 1'] = PHI**3 + 1
expressions['4*phi - 1'] = 4*PHI - 1
expressions['a1/phi + b1/phi^2'] = a1/PHI + b1/PHI**2

print(f"\n  {'Expression':>42} {'Value':>10} {'Error vs 7-phi':>12} {'Match?':>7}")
print(f"  {'-'*42} {'-'*10} {'-'*12} {'-'*7}")

for name, val in sorted(expressions.items(), key=lambda x: abs(x[1] - target)):
    err = (val/target - 1) * 100
    match = "EXACT" if abs(err) < 0.001 else ("<5%" if abs(err) < 5 else "")
    print(f"  {name:>42} {val:10.4f} {err:+12.4f}% {match:>7}")
    if abs(err) > 100:
        break  # stop printing wildly off values

# ============================================================
# PART 11: THE DERIVATION PATH
# ============================================================
print(f"\n{'='*70}")
print("PART 11: TOWARD A DERIVATION")
print(f"{'='*70}")

print(f"""
  THE CHALLENGE:
    We have the FORMULA: Omega = b1 - 1/phi = 7-phi = (26-2sqrt5)/4
    We need to show WHY the DM abundance equals this.

  WHAT WE KNOW:
    1. sin^2(tW) = b1/(a1^2+1) = 6/26 is DERIVED from:
       - b1 = ||N_{{1/2}}||^2 (TQFT Frobenius norm)
       - a1^2+1 = N_broken (Galois-broken spectral modes)
       Both are structural facts of the 600-cell geometry.

    2. The TQFT quantum dimensions satisfy:
       D^2 = a1 + sqrt(a1) (physical)
       D'^2 = a1 - sqrt(a1) (dark, Galois conjugate)
       This is a THEOREM of SU(2)_k TQFT at k+2 = a1.

    3. d_ST = 4 = a1-1 is the dimension of the orbifold C^2/2I.

  DERIVATION ATTEMPT:
    In the spectral action on C^2/2I:

    The effective number of dark DoF = N_broken (classical modes)
    corrected by the quantum dimension mismatch Delta_D^2 = 2sqrt(a1).

    The correction arises because the dark sector's TQFT
    has smaller quantum dimensions (1/phi < phi), so each
    dark mode carries LESS weight than its physical counterpart.

    The deficit per mode = Delta_D^2/N_broken = 2sqrt5/26.
    Total deficit = Delta_D^2 = 2sqrt5.
    Effective dark DoF = N_broken - Delta_D^2 = 26-2sqrt5.
    Per spacetime dimension = (26-2sqrt5)/4 = 7-phi.

  STATUS: The formula is STRUCTURAL (depends only on a1=5).
    It combines group theory (26), TQFT (2sqrt5), and geometry (4).
    All three ingredients are independently derived.
    The COMBINATION Omega = (26-2sqrt5)/4 needs a physical principle
    linking dark matter abundance to these spectral quantities.
""")

# ============================================================
# PART 12: VERLINDE AMPLITUDE CORRECTION
# ============================================================
print(f"{'='*70}")
print("PART 12: WHY A_phys/A_dark != 7-phi (AND WHAT DOES)")
print(f"{'='*70}")

# A_phys/A_dark_unsigned = D^4/(2(D'^2+2)) = 5.496
# This is NOT 7-phi = 5.382
# The difference comes from the unsigned quantum dimensions

# The unsigned A_dark uses |d'| = 1/phi, while the signed uses d' = -1/phi.
# The difference is in the n_phi=3 channels (odd number of phi-type indices).

# With SIGNED dims: ratio = phi^4 = 6.854
# With UNSIGNED dims: ratio = 5.496
# Target: 7-phi = 5.382

# What is the geometric mean?
geo_mean = np.sqrt(PHI**4 * (A_phys/A_dark_u))
print(f"  Geometric mean of signed and unsigned ratios:")
print(f"    sqrt(phi^4 * 5.496) = sqrt({PHI**4 * (A_phys/A_dark_u):.4f}) = {geo_mean:.4f}")

# What is the harmonic mean?
harm_mean = 2 / (1/PHI**4 + A_dark_u/A_phys)
print(f"  Harmonic mean: {harm_mean:.4f}")

# What linear combination?
# alpha * phi^4 + (1-alpha) * 5.496 = 7-phi
# alpha * 6.854 + (1-alpha) * 5.496 = 5.382
# alpha * (6.854-5.496) = 5.382-5.496
# alpha * 1.358 = -0.114
# alpha = -0.084
alpha_mix = (target - A_phys/A_dark_u) / (PHI**4 - A_phys/A_dark_u)
print(f"  Linear mixing parameter: alpha = {alpha_mix:.4f}")
print(f"  (Negative alpha means neither is the right observable)")

# The TRUE quantity that gives 7-phi is NOT a Verlinde amplitude ratio.
# It's the spectral formula (26-2sqrt5)/4.

# Can we RELATE the Verlinde amplitude to the spectral formula?
# A_phys = D^4 = (5+sqrt5)^2 = 30+10sqrt5
# (26-2sqrt5)/4 = 7-phi

# D^4 / what = 7-phi?
X_needed = D4 / target
print(f"\n  D^4 / X = 7-phi requires X = D^4/(7-phi) = {X_needed:.4f}")
print(f"    X = (30+10sqrt5)/(7-phi)")
print(f"    = (30+10sqrt5)*2/(13-sqrt5)")
print(f"    = (60+20sqrt5)(13+sqrt5)/164")
print(f"    = (780+60sqrt5+260sqrt5+100)/164")
print(f"    = (880+320sqrt5)/164")
print(f"    = (220+80sqrt5)/41")
print(f"    = {(220+80*SQRT5)/41:.4f}")
print(f"    = {X_needed:.4f}  (not a clean expression)")

# ============================================================
# PART 13: THE CLEANEST DERIVATION PATH
# ============================================================
print(f"\n{'='*70}")
print("PART 13: THE a1/phi + b1/phi^2 DECOMPOSITION")
print(f"{'='*70}")

# 7-phi = a1/phi + b1/phi^2 = 5/phi + 6/phi^2
# = 5*(phi-1) + 6*(2-phi)
# = 5phi - 5 + 12 - 6phi
# = 7 - phi  (checks out)

val_decomp = a1/PHI + b1/PHI**2
print(f"  7-phi = a1/phi + b1/phi^2 = {a1}/phi + {b1}/phi^2 = {val_decomp:.6f}")

# This decomposition suggests:
# Omega = a1 * d_dark + b1 * d_dark^2
# where d_dark = 1/phi = dark quantum dimension

# This is a POLYNOMIAL in d_dark evaluated at d_dark = 1/phi:
# P(x) = a1*x + b1*x^2 = 5x + 6x^2

# What is P(x) = 5x + 6x^2?
# = x*(5 + 6x)
# = x*(a1 + b1*x)

# At x=phi (physical): P(phi) = phi*(5+6phi) = phi*(5+6phi)
# = 1.618*(5+9.708) = 1.618*14.708 = 23.798
# At x=1/phi (dark): P(1/phi) = (1/phi)*(5+6/phi) = (1/phi)*(5+6(phi-1))
# = (1/phi)*(5+6phi-6) = (1/phi)*(6phi-1) = (6phi-1)/phi = 6-1/phi = b1-1/phi = 7-phi

print(f"\n  P(x) = a1*x + b1*x^2 = {a1}x + {b1}x^2")
print(f"  P(phi)  = {a1*PHI + b1*PHI**2:.4f}  (physical sector)")
print(f"  P(1/phi) = {a1/PHI + b1/PHI**2:.4f} = 7-phi  (dark sector)")
print(f"  P(phi)/P(1/phi) = {(a1*PHI+b1*PHI**2)/target:.4f}")

# The polynomial P(x) = x*(a1+b1*x) has a beautiful interpretation:
# a1 = k+2 = number of primary fields (including 0)... well, n_reps = k+1 = 4
# b1 = ||N_{1/2}||^2 = fusion norm
# x = quantum dimension of the fundamental rep

# P(x) = d_{fund} * (level_param + fusion_norm * d_{fund})
# In the physical sector: P(phi) = phi*(5+6phi) = interaction strength
# In the dark sector: P(1/phi) = (1/phi)*(5+6/phi) = dark interaction strength

# The DM abundance IS the dark sector interaction polynomial!

print(f"""
  INTERPRETATION:
    P(d) = d * (a1 + b1*d) = "interaction polynomial"

    Physical sector: P(phi) = phi*(5+6phi) = {a1*PHI+b1*PHI**2:.4f}
    Dark sector:     P(1/phi) = (1/phi)*(5+6/phi) = {target:.4f} = 7-phi

    The DM abundance equals the dark sector's interaction polynomial!

    P(d_dark) combines:
    - d_dark = 1/phi: the dark quantum dimension (propagator weight)
    - a1 = 5: the TQFT level parameter (available channels)
    - b1 = 6: the fusion norm (coupling strength)

    Omega_DM/Omega_b = P(d_dark) = d_dark * (a1 + b1*d_dark)
""")

# ============================================================
# PART 14: CHECKING AGAINST OBSERVATION
# ============================================================
print(f"{'='*70}")
print("PART 14: COMPARISON WITH PLANCK DATA")
print(f"{'='*70}")

Omega_DM_obs = 5.36  # Planck 2018
Omega_DM_err = 0.12  # approximate 1-sigma

prediction = 7 - PHI
pull = (prediction - Omega_DM_obs) / Omega_DM_err

print(f"""
  PREDICTION:  Omega_DM/Omega_b = 7-phi = {prediction:.4f}
  OBSERVED:    Omega_DM/Omega_b = {Omega_DM_obs} +/- {Omega_DM_err} (Planck 2018)
  DEVIATION:   {(prediction-Omega_DM_obs):.4f}
  PULL:        {pull:.2f} sigma
  RELATIVE:    {(prediction/Omega_DM_obs-1)*100:+.2f}%

  The prediction is within ~0.2 sigma of the observed value.
""")

# More precise Planck 2018 values
Omega_b_h2 = 0.02237  # +/- 0.00015
Omega_DM_h2 = 0.1200  # +/- 0.0012

ratio_obs = Omega_DM_h2 / Omega_b_h2
ratio_err = ratio_obs * np.sqrt((0.0012/0.1200)**2 + (0.00015/0.02237)**2)

print(f"  More precise (Planck 2018 TT,TE,EE+lowE+lensing):")
print(f"    Omega_b h^2  = {Omega_b_h2} +/- 0.00015")
print(f"    Omega_DM h^2 = {Omega_DM_h2} +/- 0.0012")
print(f"    Ratio = {ratio_obs:.4f} +/- {ratio_err:.4f}")
print(f"    7-phi = {prediction:.4f}")
print(f"    Pull  = {(prediction-ratio_obs)/ratio_err:.2f} sigma")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY OF EXP-517")
print(f"{'='*70}")

print(f"""
  ESTABLISHED:
    1. A_phys = D^4 (THEOREM: general MTC identity from S-matrix unitarity)
    2. A_dark_signed = D'^4 (Galois invariance of fusion coefficients)
    3. A_phys/A_dark != 7-phi for ANY normalization
       (signed: phi^4, unsigned: 5.496, target: 5.382)

  THE FORMULA: Omega = 7-phi = b1 - 1/phi
    = (26 - 2sqrt5)/4
    = (a1^2+1 - 2sqrt(a1))/(a1-1)
    = P(d_dark) where P(x) = a1*x + b1*x^2

  THREE EQUIVALENT REPRESENTATIONS:
    (A) b1 - 1/phi : fusion norm minus dark quantum dimension
    (B) (N_broken-Delta_D^2)/d_ST : quantum-corrected mode counting
    (C) d_dark*(a1+b1*d_dark) : dark interaction polynomial

  PHYSICAL INTERPRETATION:
    The DM abundance is the "interaction polynomial" evaluated
    at the dark quantum dimension d_dark = 1/phi = |sigma(phi)|.

    This polynomial P(x) = x*(a1+b1*x) encodes:
    - x = d_dark: propagator weight in dark TQFT
    - a1 = 5: available channels (TQFT level parameter)
    - b1 = 6: coupling strength (fusion norm ||N_{{1/2}}||^2)

  STATUS: STRUCTURAL (all ingredients derived from a1=5)
    The formula emerges from THREE independent structures:
    - Group theory of 2I (mode counting: 26)
    - SU(2)_3 TQFT (quantum dimensions: 2sqrt5)
    - C^2/2I geometry (spacetime dim: 4)
    All controlled by a1 = 5.

    OPEN: Physical principle linking P(d_dark) to Omega_DM/Omega_b.
    CANDIDATE: Spectral action vacuum energy partition.
""")

print("=" * 70)
print("EXP-517 COMPLETE")
print("=" * 70)

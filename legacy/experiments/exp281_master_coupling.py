"""
EXP-281: MASTER EQUATION FOR ALL THREE COUPLINGS
==================================================
Open problem: Is there ONE principle that gives alpha, alpha_s, sin^2(tW)?

Currently:
- alpha: 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0 (spectral + Hopf)
- alpha_s: Tr(z)=8, |N(z)|=4 -> z=2*phi^3, unique (spectral conditions)
- sin^2(tW) = b1/(a1^2+1) = 6/26 (direct algebraic)

Can we unify these into a single spectral/algebraic principle?
"""
import numpy as np
from itertools import product as iprod

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
N_eig = 9; N_gen = 3
ALPHA = 1/137.035999084
ALPHA_S_fw = 1/(2*PHI**3)  # framework value
sin2tW = b1/(a1**2 + 1)    # = 6/26

print("="*72)
print("EXP-281: MASTER EQUATION FOR ALL THREE COUPLINGS")
print("="*72)

# === Part 1: Current status ===
print("\n--- Part 1: Three Coupling Derivations ---")
print(f"  alpha^(-1) = {1/ALPHA:.6f}")
print(f"  alpha_s^(-1) = 2*phi^3 = {2*PHI**3:.6f} (framework)")
print(f"  sin^2(tW) = {b1}/{a1**2+1} = {sin2tW:.6f}")

# Derived quantities
alpha_1 = ALPHA / (1 - sin2tW)  # U(1)_Y coupling (GUT normalization)
alpha_2 = ALPHA / sin2tW         # SU(2)_L coupling
alpha_3 = ALPHA_S_fw             # SU(3)_C coupling

print(f"\n  alpha_1 = alpha/(1-sin^2tW) = {alpha_1:.6f}, 1/alpha_1 = {1/alpha_1:.4f}")
print(f"  alpha_2 = alpha/sin^2tW = {alpha_2:.6f}, 1/alpha_2 = {1/alpha_2:.4f}")
print(f"  alpha_3 = alpha_s = {alpha_3:.6f}, 1/alpha_3 = {1/alpha_3:.4f}")

# === Part 2: Z[phi] representations ===
print("\n--- Part 2: Couplings in Z[phi] ---")

# 1/alpha_s = 2*phi^3 = 2*(2*phi+1) = 2+4*phi: (a,b)=(2,4)
# Tr = 2*2+4 = 8 = rank(E8)
# N = 4+8-16 = -4

z_s = (2, 4)  # 1/alpha_s in Z[phi]
print(f"  1/alpha_s = {z_s[0]}+{z_s[1]}*phi = {z_s[0]+z_s[1]*PHI:.6f}")
print(f"    Tr = {2*z_s[0]+z_s[1]}, N = {z_s[0]**2+z_s[0]*z_s[1]-z_s[1]**2}")

# 1/alpha = 137.036... NOT in Z[phi] (involves pi)
# But the alpha equation: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
# gives alpha in terms of pi and phi

# sin^2(tW) = 6/26 = 3/13
# In Z[phi]: this is a RATIONAL number, not a Z[phi] element per se
# But 6 = b1 and 26 = a1^2+1
# 1/sin^2(tW) = 26/6 = 13/3

# Weinberg angle in Z[phi]:
# cos(2*tW) = 1 - 2*sin^2(tW) = 1 - 12/26 = 14/26 = 7/13
print(f"\n  sin^2(tW) = {b1}/{a1**2+1} = {sin2tW:.6f}")
print(f"  1/sin^2(tW) = {(a1**2+1)/b1:.6f}")
print(f"  cos(2*tW) = {1-2*sin2tW:.6f} = {a1**2-a1}/{a1**2+1}")

# === Part 3: Spectral conditions approach ===
print("\n--- Part 3: Unified Spectral Conditions ---")

# For alpha_s: find z in Z[phi] with Tr(z) = rank(SU(3))_adj, |N(z)| = dim constraint
# Can we do the same for alpha_1 and alpha_2?

# For SU(2): rank_adj = 3, dim = 3
# For U(1): rank_adj = 1, dim = 1

# Attempt: 1/g_i^2 = z_i where Tr(z_i) = C_2(G_i) and N(z_i) = some constraint

# Casimir C_2(fundamental):
# SU(3): C_2(3) = 4/3, C_2(adj) = 3
# SU(2): C_2(2) = 3/4, C_2(adj) = 2
# U(1): C_2 = Y^2

print(f"  Looking for z in Z[phi] with specific (Tr, N) for each gauge group:")
print(f"  SU(3): Tr=8, |N|=4 -> z=2+4*phi, 1/alpha_s = {2+4*PHI:.4f}")

# For SU(2): try Tr = 3 (dim of SU(2))
print(f"\n  SU(2): searching Tr = 3, various N:")
for a in range(-20, 21):
    b = 3 - 2*a  # Tr = 2a+b = 3
    z_val = a + b*PHI
    n_val = a**2 + a*b - b**2
    if z_val > 0 and abs(z_val) < 200:
        # Check if 1/z_val corresponds to alpha_2
        coupling = 1/z_val
        err = abs(coupling - alpha_2)/alpha_2 * 100
        if err < 50:
            print(f"    a={a:3d}, b={b:3d}: z={z_val:.4f}, N={n_val:5d}, 1/z={coupling:.6f}, err={err:.1f}%")

# For U(1): try Tr = 1 (dim of U(1))
print(f"\n  U(1): searching Tr = 1, various N:")
for a in range(-20, 21):
    b = 1 - 2*a  # Tr = 2a+b = 1
    z_val = a + b*PHI
    n_val = a**2 + a*b - b**2
    if z_val > 0 and abs(z_val) < 200:
        coupling = 1/z_val
        err = abs(coupling - alpha_1)/alpha_1 * 100
        if err < 50:
            print(f"    a={a:3d}, b={b:3d}: z={z_val:.4f}, N={n_val:5d}, 1/z={coupling:.6f}, err={err:.1f}%")

# === Part 4: Try other trace values ===
print("\n--- Part 4: Systematic Search ---")

# What if the trace is the DIMENSION of the representation, not the group?
# SU(3): fund = 3, adj = 8
# SU(2): fund = 2, adj = 3
# U(1): charge = 1

# Also: what about the Dynkin index?
# SU(3): T(fund) = 1/2, T(adj) = 3
# SU(2): T(fund) = 1/2, T(adj) = 2

# Let's try: 1/coupling = z with Tr(z) = specific values
print(f"  Actual inverse couplings:")
print(f"    1/alpha_1 = {1/alpha_1:.4f}")
print(f"    1/alpha_2 = {1/alpha_2:.4f}")
print(f"    1/alpha_3 = {1/alpha_3:.4f}")

# Find CLOSEST Z[phi] element to each
print(f"\n  Closest Z[phi] to 1/alpha_2:")
best_err = 999
for a in range(-5, 50):
    for b in range(-20, 50):
        z = a + b*PHI
        if z > 0:
            err = abs(z - 1/alpha_2)
            if err < best_err:
                best_err = err
                best = (a, b, z, 2*a+b, a**2+a*b-b**2)

a, b, z, tr, norm = best
print(f"    Best: {a}+{b}*phi = {z:.6f}, Tr={tr}, N={norm}")
print(f"    Actual 1/alpha_2 = {1/alpha_2:.6f}")
print(f"    Error: {abs(z-1/alpha_2):.6f}")

# 1/alpha_2 involves pi (from alpha equation)
# alpha_2 = alpha/sin^2(tW) = alpha * (a1^2+1)/b1
# So 1/alpha_2 = b1/(a1^2+1) * 1/alpha
# = sin^2(tW) / alpha
# = sin^2(tW) * (solution of quadratic involving pi)

# === Part 5: The three equations ===
print("\n--- Part 5: Three Equations from One ---")

# All three come from the alpha equation + sin^2(tW):
# Equation: 2*pi*x^2 - 4*a1*phi^4*x + 1 = 0
# Solution: alpha
# Then: alpha_2 = alpha/sin^2(tW) = alpha*(a1^2+1)/b1
# And: alpha_1 = alpha/(1-sin^2(tW)) = alpha*(a1^2+1)/(a1^2+1-b1)

# So alpha_2 and alpha_1 are NOT independent - they come from alpha + tW

# The independent equations are:
# 1. alpha equation (quadratic, involves pi)
# 2. sin^2(tW) = b1/(a1^2+1) (algebraic, exact)
# 3. alpha_s = 1/(2*phi^3) (Z[phi], exact)

# Can we get all three from spectral data?
print(f"  Equation 1: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0")
print(f"    Coefficients: (2*pi, 4*a1*phi^4, 1)")
print(f"    = (2*pi, {4*a1*PHI**4:.4f}, 1)")
print(f"    Physical: Hopf holonomy (2*pi), spectral (4*a1*phi^4), self-consistency (1)")

print(f"\n  Equation 2: sin^2(tW) = b1/(a1^2+1)")
print(f"    Physical: |Leg_B|/|total| = 6/26")
print(f"    = dim(Leg_B) / (a1^2+1)")

print(f"\n  Equation 3: 1/alpha_s = 2*phi^3")
print(f"    Physical: Tr=rank(adj SU(3))=8, |N|=N/h=4")

# === Part 6: Common structure ===
print("\n--- Part 6: Common Algebraic Structure ---")

# All three involve products of framework constants:
# sin^2(tW) = b1/(a1^2+1) - pure INTEGERS from a1=5
# 1/alpha_s = 2*phi^3 - pure Z[phi] from a1=5
# 1/alpha = root of quadratic with Z[phi] + pi coefficients

# The ONLY transcendental is pi, and it appears ONLY in alpha
# alpha_s and sin^2(tW) are ALGEBRAIC (in Z[phi] and Q respectively)

print(f"  sin^2(tW): RATIONAL (Q), from E8 legs")
print(f"  alpha_s: ALGEBRAIC (Z[phi]), from spectral conditions")
print(f"  alpha: TRANSCENDENTAL (involves pi), from Hopf holonomy")
print(f"")
print(f"  Hierarchy: Q < Z[phi] < Z[phi,pi]")
print(f"  Each coupling lives in a different number field!")

# === Part 7: The 5/3 GUT normalization ===
print("\n--- Part 7: GUT Normalization and Relations ---")

# With GUT normalization: alpha_1_GUT = (5/3)*alpha_1
alpha_1_GUT = (5./3) * alpha_1
print(f"  alpha_1 (SM) = {alpha_1:.6f}, 1/alpha_1 = {1/alpha_1:.4f}")
print(f"  alpha_1 (GUT) = {alpha_1_GUT:.6f}, 1/alpha_1_GUT = {1/alpha_1_GUT:.4f}")
print(f"  alpha_2 = {alpha_2:.6f}, 1/alpha_2 = {1/alpha_2:.4f}")
print(f"  alpha_3 = {alpha_3:.6f}, 1/alpha_3 = {1/alpha_3:.4f}")

# Check relation: 1/alpha_1_GUT = ? * 1/alpha_2
ratio_12 = (1/alpha_1_GUT) / (1/alpha_2)
print(f"\n  (1/alpha_1_GUT) / (1/alpha_2) = {ratio_12:.6f}")
# From paper: 1/alpha_1(GUT) = 2*(1/alpha_2) exactly??
# Let me check:
# 1/alpha_1_GUT = 3/(5*alpha_1) = 3*(1-sin^2tW)/(5*alpha)
# 1/alpha_2 = sin^2tW/alpha
# ratio = 3*(1-sin^2tW)/(5*sin^2tW) = 3*(1-6/26)/(5*6/26) = 3*(20/26)/(30/26) = 3*20/30 = 2

print(f"  = 3*(1-sin^2tW)/(5*sin^2tW)")
print(f"  = 3*(1-{b1}/{a1**2+1})/(5*{b1}/{a1**2+1})")
print(f"  = 3*{a1**2+1-b1}/(5*{b1})")
print(f"  = 3*{a1**2+1-b1}/{5*b1}")
print(f"  = {3*(a1**2+1-b1)}/{5*b1} = {3*(a1**2+1-b1)/(5*b1):.6f}")

# With a1=5, b1=6: 3*20/30 = 60/30 = 2
print(f"  = 2 EXACTLY (for a1=5, b1=6)")
print(f"  This is the GUT relation: 1/alpha_1(GUT) = 2/alpha_2")

# === Part 8: Sum rules ===
print("\n--- Part 8: Sum Rules for Couplings ---")

# 1/alpha_1 + 1/alpha_2 = 1/alpha * (1/(1-sin^2tW) + 1/sin^2tW)
# = (1/alpha) * (1/(1-sin^2tW) + 1/sin^2tW) = (1/alpha) / (sin^2tW*(1-sin^2tW))
# Actually: 1/alpha_1 + 1/alpha_2 = (1/alpha) * ((1-s^2+s^2)/(s^2*(1-s^2)))
# Nope, simpler: 1/alpha = 1/alpha_2 * sin^2tW = 1/alpha_1 * (1-sin^2tW)
# So: 1/alpha_1 = 1/alpha * 1/(1-s^2), 1/alpha_2 = 1/alpha * 1/s^2

sum_12 = 1/alpha_1 + 1/alpha_2
print(f"  1/alpha_1 + 1/alpha_2 = {sum_12:.4f}")
print(f"  1/alpha / (sin^2tW * (1-sin^2tW)) = {1/ALPHA / (sin2tW*(1-sin2tW)):.4f}")
# Hmm that's not right either

# Actually:
# 1/alpha_2 - 1/alpha_1 = 1/alpha * (1/sin^2tW - 1/(1-sin^2tW))
# = 1/alpha * (1-2*sin^2tW)/(sin^2tW*(1-sin^2tW))
diff_12 = 1/alpha_2 - 1/alpha_1
print(f"  1/alpha_2 - 1/alpha_1 = {diff_12:.4f}")

# More interesting: what about combinations with alpha_3?
print(f"\n  1/alpha_3 - 1/alpha_2 = {1/alpha_3 - 1/alpha_2:.4f}")
print(f"  1/alpha_2 - 1/alpha_1 = {1/alpha_2 - 1/alpha_1:.4f}")
print(f"  1/alpha_3 - 1/alpha_1 = {1/alpha_3 - 1/alpha_1:.4f}")

# Ratios
print(f"\n  1/alpha_3 / 1/alpha_2 = {(1/alpha_3)/(1/alpha_2):.6f}")
print(f"  1/alpha_2 / 1/alpha_1 = {(1/alpha_2)/(1/alpha_1):.6f}")

# Check: (1/alpha_3)/(1/alpha_2) = alpha_2/alpha_3 = alpha_s * sin^2tW / alpha
ratio_32 = ALPHA_S_fw * sin2tW / ALPHA
print(f"  alpha_s * sin^2tW / alpha = {ratio_32:.6f}")

# Check weighted sums
# In SU(5) GUT: (3/5)/alpha_1 + 1/alpha_2 = ... nah
# Check: b_i * 1/alpha_i
b_coeffs = [41./10, -19./6, -7.]  # SM beta function coefficients
print(f"\n  SM beta coefficients: b = {b_coeffs}")
for i, (name, inv_a) in enumerate(zip(['alpha_1','alpha_2','alpha_3'], [1/alpha_1, 1/alpha_2, 1/alpha_3])):
    print(f"    b_{name} * 1/{name} = {b_coeffs[i] * inv_a:.4f}")

# === Part 9: Key insight ===
print("\n--- Part 9: Key Insight ---")
print(f"  The three couplings come from THREE DIFFERENT mathematical structures:")
print(f"")
print(f"  sin^2(tW) = b1/(a1^2+1)")
print(f"    -> GRAPH THEORY: E8 leg dimensions. Rational number.")
print(f"")
print(f"  alpha_s = 1/(2*phi^3)")
print(f"    -> ALGEBRAIC NUMBER THEORY: Z[phi] spectral conditions.")
print(f"    -> Tr = rank(SU(3)_adj), |N| = |H4|/h(E8) = Coxeter orbits.")
print(f"")
print(f"  alpha from 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0")
print(f"    -> TOPOLOGY + ALGEBRA: Hopf holonomy (2*pi) + E8 spectrum.")
print(f"    -> The only coupling involving pi (continuous topology).")
print(f"")
print(f"  UNIFICATION: NOT via a single equation, but via a single INTEGER a1=5")
print(f"  that generates all three number fields: Q, Z[phi], Z[phi,pi].")
print(f"")
print(f"  The GUT relation 1/alpha_1(GUT) = 2/alpha_2 is AUTOMATIC for a1=5.")
print(f"  No running needed---vertex-transitivity IS the unification.")

# === Summary ===
print("\n" + "="*72)
print("SUMMARY")
print("="*72)
print(f"  1. Three couplings, three number fields: Q, Z[phi], Z[phi,pi]")
print(f"  2. sin^2(tW) = b1/(a1^2+1): graph combinatorics (Leg_B/total)")
print(f"  3. alpha_s = 1/(2*phi^3): Z[phi] spectral (Tr=8, |N|=4)")
print(f"  4. alpha: quadratic with pi (Hopf holonomy) + phi^4 (spectrum)")
print(f"  5. GUT relation 1/alpha_1_GUT = 2/alpha_2 is EXACT for a1=5")
print(f"  6. No 'master equation' - the master is a1=5 itself")
print(f"  STATUS: CLARIFIED (not one equation, but one integer)")
print("="*72)

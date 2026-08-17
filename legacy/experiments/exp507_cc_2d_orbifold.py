"""
EXP-507: CC from 2D Orbifold C^2/2I
=====================================
IDEA: Everything is fundamentally 2D (complex).

The 600-cell = binary icosahedral group 2I acting on C^2.
C^2/2I = Du Val singularity of type E8.
Resolution -> 8 exceptional curves in E8 pattern = McKay correspondence.

In 2D CFT on the orbifold C^2/Gamma:
  - Each element g of Gamma defines a "twisted sector"
  - Ground state energy E_0(g) depends on the eigenvalues of g
  - Total vacuum energy = (1/|Gamma|) * sum_g E_0(g)

For g in SU(2) with eigenvalues exp(+/-2*pi*i*nu):
  E_0(g) = (nu - 1/2)^2 - 1/12  (per complex boson pair on C^2)

We have the 600-cell vertices = the 120 elements of 2I as unit quaternions.
The eigenvalue nu of each element is determined by the real part w:
  cos(2*pi*nu) = w (real part of quaternion)

COMPUTE: E_vac = (1/120) sum_{g in 2I} E_0(g)
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (build_600cell, PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen,
                      alpha_em, alpha_s, ln_inv_alpha, degree, dim_E8, rank_E8, h)

print("=" * 70)
print("EXP-507: CC FROM 2D ORBIFOLD C^2/2I")
print("=" * 70)

# ============================================================
# Get 2I elements as unit quaternions
# ============================================================
verts, adj, lap = build_600cell()
print(f"\n2I group elements: {len(verts)} unit quaternions")

# Each quaternion q = (w, x, y, z) corresponds to SU(2) matrix
# with eigenvalues exp(+/-i*theta) where cos(theta) = w
# So nu = theta / (2*pi) = arccos(w) / (2*pi)

# ============================================================
# Compute nu for each group element
# ============================================================
print(f"\n{'='*70}")
print("CONJUGACY CLASSES AND TWIST PARAMETERS")
print(f"{'='*70}")

nus = []
for q in verts:
    w = q[0]  # real part of quaternion
    w_clipped = np.clip(w, -1.0, 1.0)
    theta = np.arccos(w_clipped)  # in [0, pi]
    nu = theta / (2 * np.pi)     # in [0, 0.5]
    nus.append(nu)

nus = np.array(nus)

# Group by conjugacy class (same nu up to rounding)
nu_classes = {}
for i, nu in enumerate(nus):
    nu_r = round(nu, 6)
    if nu_r not in nu_classes:
        nu_classes[nu_r] = []
    nu_classes[nu_r].append(i)

print(f"\nConjugacy classes (by twist parameter nu):")
print(f"{'nu':>10} {'|class|':>8} {'order':>6} {'2*cos(2pi*nu)':>14} {'E_0(nu)':>12}")

total_E = 0
for nu_val in sorted(nu_classes.keys()):
    size = len(nu_classes[nu_val])
    # Order of element: smallest n such that n*nu is integer (or n*nu = half-integer for -I)
    if nu_val < 1e-8:
        order = 1
    else:
        # Find order: eigenvalue = exp(2*pi*i*nu), order = smallest n with n*nu in Z
        for n in range(1, 31):
            if abs(n * nu_val - round(n * nu_val)) < 1e-6:
                order = n
                break

    char_val = 2 * np.cos(2 * np.pi * nu_val)

    # Ground state energy for C^2 twisted by g:
    # E_0(nu) = (nu - 1/2)^2 - 1/12
    E_0 = (nu_val - 0.5)**2 - 1/12

    total_E += size * E_0

    print(f"{nu_val:10.6f} {size:8d} {order:6d} {char_val:14.6f} {E_0:12.6f}")

# Total vacuum energy
E_vac = total_E / N
print(f"\nTotal weighted E: sum |C_i| * E_0(nu_i) = {total_E:.6f}")
print(f"E_vac = total / |2I| = {total_E:.6f} / {N} = {E_vac:.10f}")

# ============================================================
# What IS this number?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS OF E_vac")
print(f"{'='*70}")

print(f"\n  E_vac = {E_vac:.10f}")
print(f"  1/E_vac = {1/E_vac:.6f}")

# Check against known quantities
tests = {
    '1/12': 1/12,
    '-1/12': -1/12,
    '1/24': 1/24,
    '-1/24': -1/24,
    '1/N': 1/N,
    '1/(2N)': 1/(2*N),
    'a1/N': a1/N,
    '1/(a1*N)': 1/(a1*N),
    'b1/N': b1/N,
    '1/(b1*N)': 1/(b1*N),
    '1/h': 1/h,
    '1/(2h)': 1/(2*h),
    '(a1-1)/(12*N)': (a1-1)/(12*N),
    'rank/(12*N)': rank_E8/(12*N),
    '1/(12*degree)': 1/(12*degree),
    '-c_eff/24': -(a1-1)/(24),  # effective central charge?
}

print(f"\n  Checking against framework quantities:")
for name, val in sorted(tests.items(), key=lambda x: abs(x[1] - E_vac)):
    if abs(val) > 1e-15:
        ratio = E_vac / val
        if abs(ratio - round(ratio)) < 0.01 or abs(1/ratio - round(1/ratio)) < 0.01:
            print(f"    {name:>20} = {val:.10f}  ratio = {ratio:.6f} ***")
        else:
            print(f"    {name:>20} = {val:.10f}  ratio = {ratio:.6f}")

# Try to express E_vac as a simple fraction
# E_vac should be a rational number (sum of rationals)
# Let's compute exactly using fractions
from fractions import Fraction

# Recompute with exact arithmetic
# nu values are arccos(w)/(2*pi), which are NOT rational in general
# for quaternion vertices involving phi.
# But the sum might simplify.

# Actually, the nu values for 2I elements are determined by the
# conjugacy classes. Let me identify them exactly.

print(f"\n{'='*70}")
print("EXACT CONJUGACY CLASS DATA")
print(f"{'='*70}")

# For 2I, the character of the 2D fundamental rep gives:
# chi(g) = 2*cos(2*pi*nu_g)
# The character values are: 2, -2, phi, -phi, 1, -1, 1/phi, -1/phi, 0
# (from McKay eigenvalues of the extended E8)

# McKay eigenvalues = {2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2}
# with multiplicities = {1, 4, 9, 16, 25, 16, 9, 4, 1}
# Wait no - those are Laplacian eigenvalue multiplicities, not class sizes.

# Character values of the 2D rep of 2I:
# 2I has 9 conjugacy classes. The 2D rep has characters:
# {2, -2, phi, -phi, 1, -1, -1/phi, 1/phi, 0}
# But which sizes? Let me use the known class sizes.

# 2I conjugacy classes (standard):
# Class: Identity  -I   C5   C5^2  C5^3  C5^4  C3   C3^2  C2
# |C|:      1       1   12    12    12    12    20    20    30
# chi_2:    2      -2   phi  1/phi -1/phi -phi   1    -1    0
# Order:    1       2   10    10    10    10    6     6     4

# But 1+1+12+12+12+12+20+20+30 = 120. Good.

# Let me verify: chi_2 = 2*cos(2*pi*nu)
# chi = 2 -> nu = 0 (identity)
# chi = -2 -> nu = 1/2 (-I)
# chi = phi = (1+sqrt(5))/2 -> cos(2*pi*nu) = phi/2 -> 2*pi*nu = arccos(phi/2) = pi/5 -> nu = 1/10
# chi = 1/phi = (sqrt(5)-1)/2 -> cos(2*pi*nu) = 1/(2*phi) = (sqrt(5)-1)/4 -> nu = ?
#   cos(2*pi*nu) = (sqrt(5)-1)/4 = cos(2*pi/5) ... wait
#   Actually 1/(2*phi) = 1/(1+sqrt(5)) = (sqrt(5)-1)/4? Let me check:
#   1/phi = (sqrt(5)-1)/2, so 1/(2*phi) = (sqrt(5)-1)/4
#   cos(72 deg) = cos(2*pi/5) = (sqrt(5)-1)/4. Yes!
#   So nu = 1/5 (= 2/10)
# chi = -1/phi -> cos(2*pi*nu) = -1/(2*phi) = -(sqrt(5)-1)/4 = -cos(72 deg) = cos(108 deg)
#   108 deg = 3*pi/5 -> nu = 3/10
# chi = -phi -> cos(2*pi*nu) = -phi/2 -> arccos(-phi/2) = pi - pi/5 = 4*pi/5 -> nu = 2/5
# chi = 1 -> cos(2*pi*nu) = 1/2 -> nu = 1/6
# chi = -1 -> cos(2*pi*nu) = -1/2 -> nu = 1/3
# chi = 0 -> cos(2*pi*nu) = 0 -> nu = 1/4

conj_classes = [
    # (chi_2, size, nu, label)
    (2,    1,  Fraction(0, 1),   "Identity"),
    (-2,   1,  Fraction(1, 2),   "-I"),
    (PHI, 12,  Fraction(1, 10),  "C5"),
    (1/PHI, 12, Fraction(1, 5),  "C5^2"),
    (-1/PHI, 12, Fraction(3, 10), "C5^3"),
    (-PHI, 12, Fraction(2, 5),   "C5^4"),
    (1,   20,  Fraction(1, 6),   "C3"),
    (-1,  20,  Fraction(1, 3),   "C3^2"),
    (0,   30,  Fraction(1, 4),   "C2"),
]

# Verify sizes sum to 120
total_size = sum(size for _, size, _, _ in conj_classes)
assert total_size == 120, f"Sizes sum to {total_size}"

print(f"\nExact conjugacy class data:")
print(f"{'Label':>8} {'chi_2':>8} {'|C|':>5} {'nu':>6} {'E_0(nu)':>15} {'|C|*E_0':>15}")

total_E_exact = Fraction(0)
for chi, size, nu, label in conj_classes:
    # E_0(nu) = (nu - 1/2)^2 - 1/12
    E_0_exact = (nu - Fraction(1, 2))**2 - Fraction(1, 12)
    contribution = size * E_0_exact
    total_E_exact += contribution

    print(f"{label:>8} {chi:8.4f} {size:5d} {str(nu):>6} {str(E_0_exact):>15} {str(contribution):>15}")

E_vac_exact = total_E_exact / 120
print(f"\nExact total: sum |C|*E_0 = {total_E_exact} = {float(total_E_exact):.10f}")
print(f"E_vac = total/120 = {E_vac_exact} = {float(E_vac_exact):.10f}")

# ============================================================
# Interpret the result
# ============================================================
print(f"\n{'='*70}")
print("INTERPRETATION")
print(f"{'='*70}")

E_val = float(E_vac_exact)
print(f"\n  E_vac(C^2/2I) = {E_vac_exact} = {E_val:.10f}")
print(f"  = {E_vac_exact.numerator}/{E_vac_exact.denominator}")

# Simplify
num = E_vac_exact.numerator
den = E_vac_exact.denominator
print(f"  = {num}/{den}")

# Factor
import math
g = math.gcd(abs(num), den)
print(f"  = {num//g}/{den//g} (simplified)")

# Check against theory constants
print(f"\n  Framework interpretations:")
print(f"    1/(12*N) = 1/{12*N} = {1/(12*N):.10f}")
print(f"    1/h = 1/{h} = {1/h:.10f}")
print(f"    a1/(12*N) = {a1}/{12*N} = {a1/(12*N):.10f}")
print(f"    E_vac / (1/N) = {E_val * N:.6f}")
print(f"    E_vac * N = {E_val * N:.6f}")
print(f"    E_vac * 12 = {E_val * 12:.6f}")
print(f"    E_vac * 24 = {E_val * 24:.6f}")
print(f"    E_vac * N * 12 = {E_val * N * 12:.6f}")
print(f"    E_vac * h = {E_val * h:.6f}")

# Central charge interpretation: E_vac = -c_eff / 24
c_eff = -24 * E_val
print(f"\n  If E_vac = -c_eff/24: c_eff = {c_eff:.6f}")
print(f"    c_eff = {c_eff:.6f}")
print(f"    c_eff * N = {c_eff * N:.6f}")
print(f"    c_eff / a1 = {c_eff / a1:.6f}")

# ============================================================
# Connection to CC
# ============================================================
print(f"\n{'='*70}")
print("CONNECTION TO COSMOLOGICAL CONSTANT")
print(f"{'='*70}")

# In 2D CFT, the partition function on a torus of modular parameter tau:
# Z(tau) = exp(-2*pi*tau_2 * E_vac) * (1 + ...)
# The vacuum energy gives the exponential suppression.

# If we identify: the "size" of the orbifold torus with 1/alpha somehow,
# then the CC would be:
# Lambda ~ exp(-2*pi * (1/alpha) * E_vac) or similar

# Let's check: what exponent does E_vac * something give?
# We need: the exponent to be ~ 282 (= z_CC * ln(1/alpha) ~ 57 * 4.92)

print(f"\n  E_vac = {E_val:.10f}")
print(f"  E_vac * 2*pi/alpha = {E_val * 2*np.pi/alpha_em:.4f}")
print(f"  E_vac * N * ln(1/alpha) = {E_val * N * ln_inv_alpha:.4f}")
print(f"  E_vac * N^2 = {E_val * N**2:.4f}")
print(f"  E_vac * dim(E8) = {E_val * dim_E8:.6f}")

# The z-equivalent:
# Lambda = exp(-S) where S = E_vac * (something)
# If Lambda = alpha^z, then S = z * ln(1/alpha)
# So z = S / ln(1/alpha) = E_vac * (something) / ln(1/alpha)

print(f"\n  For z = 57: need E_vac * X = 57 * ln(1/alpha) = {57*ln_inv_alpha:.4f}")
print(f"  X = {57*ln_inv_alpha/E_val:.4f}")
print(f"  X/N = {57*ln_inv_alpha/(E_val*N):.6f}")
print(f"  X/N^2 = {57*ln_inv_alpha/(E_val*N**2):.6f}")

# ============================================================
# Alternative: Casimir energy as CC directly
# ============================================================
print(f"\n{'='*70}")
print("CASIMIR ENERGY AS CC")
print(f"{'='*70}")

# In Kaluza-Klein: if extra dimensions are compactified on C^2/2I,
# the Casimir energy of the compact space IS the CC.
# Lambda = (E_vac / V_compact) in appropriate units.

# For a compact space of "radius" R:
# E_Casimir ~ E_vac / R^4 (for 4 real compact dimensions = 2 complex)
# In Planck units: Lambda_P = E_vac * (l_P/R)^4

# If R = l_P / alpha^z for some z:
# Lambda_P = E_vac * alpha^(4z)

# Or: if R = l_P * phi^n:
# Lambda_P = E_vac / phi^(4n)

# For Lambda_P ~ 10^{-122}:
# E_vac / phi^(4n) ~ alpha^57
# phi^(4n) ~ E_vac / alpha^57 ~ E_val / 10^{-122}

print(f"\n  E_vac = {E_val:.6e}")
print(f"  alpha^57 = {alpha_em**57:.6e}")
print(f"  E_vac / alpha^57 = {E_val / alpha_em**57:.6e}")
print(f"  ln(E_vac/alpha^57) / ln(phi) = {np.log(E_val/alpha_em**57)/np.log(PHI):.4f}")
print(f"  So phi^{np.log(E_val/alpha_em**57)/np.log(PHI):.0f} = E_vac/alpha^57")

# ============================================================
# The orbifold Euler characteristic
# ============================================================
print(f"\n{'='*70}")
print("ORBIFOLD EULER CHARACTERISTIC")
print(f"{'='*70}")

# For the orbifold C^2/Gamma:
# chi_orb = (1/|Gamma|) * sum_g chi(Fix(g))
# For C^2: chi(C^2) = 1 (contractible)
# Fix(g) for g != I: the origin only (for elements acting freely away from 0)
# So chi_orb = chi(C^2)/|Gamma| + sum_{g!=I} chi(origin)/|Gamma|
# = 1/120 + 119/120 = 1? No, this isn't right.

# Actually, for the RESOLUTION of C^2/2I:
# chi = 1 + 8 = 9 (one for the base + 8 exceptional curves)
# = 1 + rank(E8) = 9 = N_eig (number of eigenvalues!)

chi_resolution = 1 + rank_E8
print(f"  chi(resolution of C^2/2I) = 1 + rank(E8) = {chi_resolution}")
print(f"  = N_eig = 9 (number of distinct Laplacian eigenvalues!)")

# String-theoretic Euler characteristic of the orbifold:
# chi_string = (1/|Gamma|) * sum_{g,h: [g,h]=1} 1
# = (1/|Gamma|) * sum_g |Z(g)|
# = number of conjugacy classes
chi_string = len(conj_classes)
print(f"  chi_string(C^2/2I) = # conjugacy classes = {chi_string}")
print(f"  = N_eig = 9 (same!)")

# ============================================================
# KEY: E_vac in terms of the resolution
# ============================================================
print(f"\n{'='*70}")
print("E_vac AND THE RESOLUTION")
print(f"{'='*70}")

# On the resolved surface, the Euler characteristic is 9.
# The Hirzebruch signature is 8 (= rank E8).
# The arithmetic genus is chi_a = (1 + K^2)/12 where K is canonical class.

# For the minimal resolution of C^2/2I (= ALE space of type E8):
# K = 0 (the resolution is Calabi-Yau, i.e., c_1 = 0)
# chi_a = 1/12 * (chi_top + sigma) = (9 + (-8))/12 = 1/12

# So: chi_a = 1/12 for the resolved E8 surface!
chi_a = Fraction(1, 12)
print(f"  Arithmetic genus chi_a = {chi_a} = {float(chi_a):.6f}")
print(f"  E_vac = {E_vac_exact}")
print(f"  E_vac / chi_a = {E_vac_exact / chi_a} = {float(E_vac_exact/chi_a):.6f}")
print(f"  E_vac * 12 = {float(E_vac_exact * 12):.6f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
2D ORBIFOLD C^2/2I = E8 Du Val singularity

COMPUTED:
  E_vac = (1/120) * sum_{{g in 2I}} [(nu_g - 1/2)^2 - 1/12]
        = {E_vac_exact} = {float(E_vac_exact):.10f}

  Effective central charge: c_eff = -24*E_vac = {c_eff:.6f}

  Orbifold Euler characteristic = {chi_string} = N_eig
  Resolution Euler characteristic = {chi_resolution} = 1 + rank(E8)
  Arithmetic genus = {chi_a}

WHAT THIS TELLS US:
  The 2D orbifold perspective gives a DEFINITE vacuum energy {E_vac_exact}.
  This is a pure NUMBER from the group theory of 2I.
  All inputs are derived from a1 = 5.

  The CC connection requires identifying the "modular parameter"
  of the orbifold with a physical scale. This is where
  alpha and the exponential form would enter.

STATUS: The 2D perspective is mathematically clean and gives
  a definite E_vac. The connection to the observed CC
  requires additional physics (KK compactification scale,
  modular identification, etc.).
""")

print("=" * 70)
print("EXP-507 COMPLETE")
print("=" * 70)

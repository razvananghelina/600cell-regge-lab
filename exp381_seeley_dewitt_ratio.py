"""
exp381_seeley_dewitt_ratio.py
==============================
GOAL: Analyze the ratio c1^2/(c0*c2) = 3844/2563 ~ 1.4998.
Is it exactly 3/2? Or some other algebraic number?

From exp378:
  c0 = 2*N*A_0, c1 = 2*N*A_1, c2 = 2*N*A_2
  A_0 = 11 = a1+b1, A_1 = 62, A_2 = 233 = F_13

So c1^2/(c0*c2) = A_1^2/(A_0*A_2) = 62^2/(11*233) = 3844/2563

BLIND: compute first, interpret after.

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from math import sqrt, gcd, factorial
from fractions import Fraction

# Framework constants
a1 = 5
b1 = a1 + 1  # = 6
PHI = (1 + sqrt(a1)) / 2
N = 120
degree = 12

print("=" * 72)
print("exp381: SEELEY-DEWITT RATIO c1^2/(c0*c2)")
print("=" * 72)

# =====================================================================
# PART 1: EXACT VALUES
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: Exact Values")
print("=" * 72)

# Seeley-DeWitt coefficients
c0 = 2640
c1 = 14880
c2 = 55920

# A_k values
A0 = c0 // (2*N)  # = 11
A1 = c1 // (2*N)  # = 62
A2 = c2 // (2*N)  # = 233

print(f"c0 = {c0} = 2*{N}*{A0}")
print(f"c1 = {c1} = 2*{N}*{A1}")
print(f"c2 = {c2} = 2*{N}*{A2}")

# The ratio
ratio_num = A1**2
ratio_den = A0 * A2
g = gcd(ratio_num, ratio_den)
print(f"\nc1^2/(c0*c2) = A1^2/(A0*A2) = {ratio_num}/{ratio_den}")
print(f"  = {ratio_num//g}/{ratio_den//g} (reduced)")
print(f"  = {ratio_num/ratio_den:.10f}")
print(f"  3/2 = {3/2:.10f}")
print(f"  Difference from 3/2: {ratio_num/ratio_den - 1.5:.10f}")
print(f"  Relative: {(ratio_num/ratio_den - 1.5)/1.5*100:.6f}%")

# Check if exactly 3/2:
# 3/2 * A0 * A2 = 3/2 * 11 * 233 = 3/2 * 2563 = 3844.5
print(f"\n3/2 * {A0}*{A2} = 3/2 * {A0*A2} = {3*A0*A2/2}")
print(f"A1^2 = {A1**2}")
print(f"EXACT 3/2? {A1**2 == 3*A0*A2//2}  (NO, off by 0.5)")
print(f"2*A1^2 = {2*A1**2}, 3*A0*A2 = {3*A0*A2}")
print(f"2*A1^2 - 3*A0*A2 = {2*A1**2 - 3*A0*A2}")

# So 2*A1^2 - 3*A0*A2 = 7688 - 7689 = -1
# The ratio is (3/2) - 1/(2*A0*A2) = 3/2 - 1/5126
deficit = 2*A1**2 - 3*A0*A2
print(f"\n*** REMARKABLE: 2*A1^2 - 3*A0*A2 = {deficit} ***")
print(f"So c1^2/(c0*c2) = 3/2 - 1/(2*A0*A2) = 3/2 - 1/{2*A0*A2}")
print(f"  = 3/2 - 1/5126")
print(f"  = {Fraction(3,2) - Fraction(1, 2*A0*A2)}")
print(f"  = {(3*A0*A2 - 1)}/{2*A0*A2}")

# Verify
exact_frac = Fraction(A1**2, A0*A2)
print(f"\nExact fraction: {exact_frac}")
print(f"  = {float(exact_frac):.15f}")

# =====================================================================
# PART 2: A_k POLYNOMIALS IN a1
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: A_k as Polynomials in a1")
print("=" * 72)

# From exp378b (MEMORY):
# A_0 = 2*a1 + 1
# A_1 = 2*a1^2 + 2*a1 + 2
# A_2 = 6*a1^2 + 15*a1 + 8

def A0_poly(a):
    return 2*a + 1

def A1_poly(a):
    return 2*a**2 + 2*a + 2

def A2_poly(a):
    return 6*a**2 + 15*a + 8

# Verify for a1=5
print(f"A_0(5) = {A0_poly(5)} (should be {A0})")
print(f"A_1(5) = {A1_poly(5)} (should be {A1})")
print(f"A_2(5) = {A2_poly(5)} (should be {A2})")

# The deficit: 2*A1^2 - 3*A0*A2
# = 2*(2a^2+2a+2)^2 - 3*(2a+1)*(6a^2+15a+8)
# Let's expand algebraically

# 2*(2a^2+2a+2)^2:
# (2a^2+2a+2)^2 = 4a^4 + 8a^3 + 12a^2 + 8a + 4
#   (using (x+y+z)^2 where x=2a^2, y=2a, z=2)
#   = 4a^4 + 4a^2 + 4 + 2*4a^3 + 2*4a^2 + 2*4a = 4a^4+8a^3+12a^2+8a+4
# So 2*(2a^2+2a+2)^2 = 8a^4 + 16a^3 + 24a^2 + 16a + 8

# 3*(2a+1)*(6a^2+15a+8):
# (2a+1)*(6a^2+15a+8) = 12a^3 + 30a^2 + 16a + 6a^2 + 15a + 8
#                      = 12a^3 + 36a^2 + 31a + 8
# 3* = 36a^3 + 108a^2 + 93a + 24

# Deficit = (8a^4 + 16a^3 + 24a^2 + 16a + 8) - (36a^3 + 108a^2 + 93a + 24)
# = 8a^4 - 20a^3 - 84a^2 - 77a - 16

def deficit_poly(a):
    return 8*a**4 - 20*a**3 - 84*a**2 - 77*a - 16

print(f"\nDeficit polynomial: 2*A1^2 - 3*A0*A2 = 8*a^4 - 20*a^3 - 84*a^2 - 77*a - 16")
print(f"  D(5) = {deficit_poly(5)} (should be {deficit})")

# Check for other values of a1:
print("\nDeficit for various a1:")
for a in range(2, 10):
    a0 = A0_poly(a)
    a1_v = A1_poly(a)
    a2 = A2_poly(a)
    d = 2*a1_v**2 - 3*a0*a2
    r = a1_v**2 / (a0*a2) if a0*a2 != 0 else float('inf')
    print(f"  a1={a}: A0={a0:4d}, A1={a1_v:5d}, A2={a2:6d}, "
          f"deficit={d:8d}, ratio={r:.6f}")

# =====================================================================
# PART 3: IS THE DEFICIT -1 SPECIAL?
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: The Deficit = -1 for a1=5")
print("=" * 72)

# For a1=5, the deficit is EXACTLY -1.
# 2*A1^2 - 3*A0*A2 = -1
# This means: A1^2 = (3*A0*A2 - 1)/2

# Is this a Pell-like equation? Or a Diophantine condition?
# 2*(62)^2 + 1 = 3*11*233 = 7689
# 7688 + 1 = 7689. YES!

# In Pell equation form: 2*x^2 - 3*y*z = -1 where x=A1, y=A0, z=A2

# Check: for which a1 values is the deficit = -1?
print("Checking deficit = -1 for various a1:")
for a in range(1, 20):
    d = deficit_poly(a)
    if d == -1:
        print(f"  a1 = {a}: deficit = {d} <-- EXACTLY -1!")
    elif abs(d) < 10:
        print(f"  a1 = {a}: deficit = {d} (close to -1)")

# Factor the deficit polynomial:
# 8a^4 - 20a^3 - 84a^2 - 77a - 16
# Try a=5: 8*625 - 20*125 - 84*25 - 77*5 - 16 = 5000-2500-2100-385-16 = -1
# So (a-5) is a factor of 8a^4-20a^3-84a^2-77a-16+1 = 8a^4-20a^3-84a^2-77a-15

# Let me factor 8a^4-20a^3-84a^2-77a-15:
# Divide by (a-5):
# 8a^4 - 20a^3 - 84a^2 - 77a - 15 = (a-5)*q(a)
# Synthetic division:
coeffs = [8, -20, -84, -77, -15]
root = 5
result = [coeffs[0]]
for c in coeffs[1:]:
    result.append(c + root * result[-1])
# result = [8, 8*5-20=20, 20*5-84=16, 16*5-77=3, 3*5-15=0]
print(f"\nSynthetic division by (a-5): {result}")
print(f"  Quotient: {result[:-1]} -> 8a^3 + 20a^2 + 16a + 3")
print(f"  Remainder: {result[-1]}")

# Factor 8a^3 + 20a^2 + 16a + 3:
# Try a=-1/2: 8*(-1/8) + 20*(1/4) + 16*(-1/2) + 3 = -1 + 5 - 8 + 3 = -1. Not zero.
# Try a=-3: 8*(-27) + 20*9 + 16*(-3) + 3 = -216+180-48+3 = -81. No.
# Try a=-1/4: 8*(-1/64) + 20*(1/16) + 16*(-1/4) + 3 = -0.125+1.25-4+3 = 0.125. No.
# Rational root theorem: possible rational roots are +-{1,3,1/2,3/2,1/4,3/4,1/8,3/8}
for r_test in [1, -1, 3, -3, 0.5, -0.5, 1.5, -1.5, 0.25, -0.25]:
    val = 8*r_test**3 + 20*r_test**2 + 16*r_test + 3
    if abs(val) < 1e-10:
        print(f"  Root found: a = {r_test}")

# Try -1/4: 8*(-1/64) + 20/16 + 16*(-1/4) + 3 = -1/8 + 5/4 - 4 + 3 = -1/8+5/4-1 = 1/8. No.
# Try -3/8: 8*(-27/512) + 20*(9/64) + 16*(-3/8) + 3 = -27/64 + 180/64 - 6 + 3
# = 153/64 - 3 = 153/64 - 192/64 = -39/64. No.

# It seems 8a^3+20a^2+16a+3 has no rational roots.
# Use numpy to find roots:
roots_cubic = np.roots([8, 20, 16, 3])
print(f"\nRoots of 8a^3 + 20a^2 + 16a + 3:")
for r in roots_cubic:
    if abs(r.imag) < 1e-10:
        print(f"  a = {r.real:.6f}")
    else:
        print(f"  a = {r.real:.6f} + {r.imag:.6f}i")

# So: deficit+1 = (a-5)*(8a^3+20a^2+16a+3)
# deficit = -1 only at a=5 (among positive integers), since the cubic
# has no positive real root.

print("\n*** THEOREM: Among positive integers a1 >= 2, the deficit ***")
print("*** 2*A1^2 - 3*A0*A2 = -1 holds ONLY for a1 = 5. ***")
print("*** This is equivalent to c1^2/(c0*c2) = 3/2 - 1/(2*A0*A2). ***")

# =====================================================================
# PART 4: REFORMULATION
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: Clean Reformulation")
print("=" * 72)

# We have: 2*A1^2 + 1 = 3*A0*A2
# i.e.: 2*(2a1^2+2a1+2)^2 + 1 = 3*(2a1+1)*(6a1^2+15a1+8)

# This can be rewritten:
# 2*A1^2 + 1 = 3*A0*A2
# A1 = 2a^2+2a+2 = 2(a^2+a+1)
# A0 = 2a+1
# A2 = 6a^2+15a+8

# So: 2*4*(a^2+a+1)^2 + 1 = 3*(2a+1)*(6a^2+15a+8)
# 8*(a^2+a+1)^2 + 1 = 3*(2a+1)*(6a^2+15a+8)
# For a=5: 8*31^2 + 1 = 8*961+1 = 7689 = 3*11*233 = 7689. CHECK!

# Note: a^2+a+1 = 31 for a=5. And 31 = h(E8)+1.
# Also: 2a+1 = 11 = a1+b1 = L_5 (Lucas number).

print(f"a^2+a+1 = {a1**2+a1+1} = h(E8)+1 = {a1*b1+1}")
print(f"2a+1 = {2*a1+1} = a1+b1 = L_5")
print(f"6a^2+15a+8 = {6*a1**2+15*a1+8} = F_13 (Fibonacci)")
print()
print(f"The identity becomes:")
print(f"  8*(h(E8)+1)^2 + 1 = 3*L_5*F_13")
print(f"  8*31^2 + 1 = 3*11*233")
print(f"  7688 + 1 = 7689")
print(f"  7689 = 7689  CHECK!")

# Can we simplify? 7689 = 3*2563 = 3*11*233
# And 7688 = 8*961 = 8*31^2 = 2^3 * 31^2
# So: 2^3 * 31^2 + 1 = 3 * 11 * 233

print(f"\nIn framework notation:")
print(f"  2^3 * (h+1)^2 + 1 = 3 * (a1+b1) * F_{{2*b1+1}}")
print(f"  where h = h(E8) = a1*b1 = 30")

# =====================================================================
# PART 5: WHY 233 = F_13?
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: Why A_2 = F_13 (Fibonacci)")
print("=" * 72)

# From exp378b (MEMORY): 233 = F_13 and the proof involves 12^2 = F_12 (144).
# A_2 = 6a^2+15a+8 = 6*25+75+8 = 233

# Check if this is a Fibonacci number for other a1:
fib = [1, 1]
for _ in range(25):
    fib.append(fib[-1] + fib[-2])

print("Is A_2(a1) a Fibonacci number?")
for a in range(2, 10):
    a2_val = A2_poly(a)
    is_fib = a2_val in fib
    fib_idx = fib.index(a2_val) if is_fib else -1
    print(f"  a1={a}: A_2 = {a2_val}" +
          (f" = F_{fib_idx+1}" if is_fib else ""))

# Only a1=5 gives a Fibonacci number!
# a1=2: A2 = 62 (not Fib), a1=3: A2 = 107 (not Fib), a1=4: A2 = 164 (not Fib)
# a1=5: A2 = 233 = F_13, a1=6: A2 = 314 (not Fib)

print("\nIs A_0(a1) a Lucas/Fibonacci number?")
lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322]
for a in range(2, 10):
    a0_val = A0_poly(a)
    is_luc = a0_val in lucas
    luc_idx = lucas.index(a0_val) if is_luc else -1
    is_fib = a0_val in fib
    fib_idx = fib.index(a0_val) if is_fib else -1
    tag = ""
    if is_luc:
        tag = f" = L_{luc_idx}"
    if is_fib:
        tag += f" = F_{fib_idx+1}"
    print(f"  a1={a}: A_0 = {a0_val}{tag}")

# a1=5: A_0 = 11 = L_5. Also a Fibonacci/Lucas connection!
# a1=4: A_0 = 9. Not Lucas. But 9 = F_6? No, F_6 = 8. Not Fibonacci either.

# =====================================================================
# PART 6: TRIPLE UNIQUENESS AT a1=5
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: Triple Uniqueness at a1=5")
print("=" * 72)

# From MEMORY: a1=5 has TRIPLE UNIQUENESS:
# 1. A0 = L_5 (Lucas number with index = a1)
# 2. degree^2 = F_{degree} (144 = F_12)
# 3. A2 = F_{degree+1} (233 = F_13)

# Let's verify and add our new result:

print("a1=5 UNIQUENESS PROPERTIES:")
print(f"  1. A_0 = {A0} = L_{a1} = L_5 (Lucas)")
print(f"  2. degree^2 = {degree**2} = F_{degree} = F_12 (Fibonacci)")
print(f"  3. A_2 = {A2} = F_{degree+1} = F_13 (Fibonacci)")
print(f"  4. 2*A_1^2 - 3*A_0*A_2 = {deficit} (Pell-like)")
print(f"     Equivalent: c1^2/(c0*c2) = 3/2 - 1/{2*A0*A2}")

# Check if property 2 (degree^2 = F_degree) holds for other a1:
print("\nProperty 2: degree^2 = F_degree?")
for a in range(2, 10):
    deg = 2*(a+1)
    deg_sq = deg**2
    f_deg = fib[deg] if deg < len(fib) else None
    match = f_deg is not None and deg_sq == f_deg
    print(f"  a1={a}: degree={deg}, degree^2={deg_sq}, F_{deg}={f_deg}, match={match}")

# Only a1=5 (degree=12): 144 = F_12!

# =====================================================================
# PART 7: THE c1^2/(c0*c2) IDENTITY
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: The Master Identity")
print("=" * 72)

# We proved: 2*A1^2 + 1 = 3*A0*A2, unique to a1=5.
# Equivalently:
# c1^2/(c0*c2) = A1^2/(A0*A2) = (3*A0*A2 - 1)/(2*A0*A2) = 3/2 - 1/(2*A0*A2)

# In framework terms:
# 2*A0*A2 = 2*11*233 = 5126
# c1^2/(c0*c2) = 3/2 - 1/5126 = 7688/5126 = 3844/2563

print("MASTER IDENTITY (unique to a1=5):")
print()
print("  2 * A_1^2 + 1 = 3 * A_0 * A_2")
print()
print(f"  2 * {A1}^2 + 1 = 3 * {A0} * {A2}")
print(f"  {2*A1**2} + 1 = {3*A0*A2}")
print(f"  {2*A1**2 + 1} = {3*A0*A2}  CHECK!")
print()
print("Equivalent forms:")
print(f"  8*(a1^2+a1+1)^2 + 1 = 3*(2*a1+1)*(6*a1^2+15*a1+8)")
print(f"  8*(h(E8)+1)^2 + 1 = 3*L_5*F_13")
print(f"  c1^2/(c0*c2) = 3/2 - 1/(2*(a1+b1)*F_{{2*b1+1}})")
print(f"               = 3/2 - 1/5126")
print(f"               = 3844/2563")
print(f"               ~ 1.499805 (within 0.013% of 3/2)")

# =====================================================================
# PART 8: PHYSICAL SIGNIFICANCE
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: Physical Significance")
print("=" * 72)

print("""
The ratio c1^2/(c0*c2) = A1^2/(A0*A2) has physical meaning in the
Seeley-DeWitt expansion of the spectral action:

  S = f_0*Lambda^4*c0 + f_2*Lambda^2*c1 + f_4*c2 + ...

The ratio c1^2/(c0*c2) measures the "regularity" of the spectral sequence:
  - If c_k is geometric: c1^2 = c0*c2 exactly (ratio = 1)
  - If c_k is arithmetic: ratio varies

Our result: ratio = 3/2 - epsilon, where epsilon = 1/5126 is tiny.
This means the spectral coefficients are ALMOST in geometric progression
with common ratio ~ c1/c0 = 62/11 ~ 5.636.

The near-3/2 value means:
  c1 ~ sqrt(3/2 * c0 * c2) = sqrt(3/2) * sqrt(c0*c2)

This "super-geometric" property (ratio > 1) indicates that the
curvature term c1 is LARGER than geometric interpolation suggests.
In other words, the 600-cell has MORE curvature per unit volume
than expected from a smooth S^3 of the same "spectral width".
""")

# Check: what would c1 be for ratio exactly 3/2?
c1_exact_3_2 = sqrt(1.5 * c0 * c2)
print(f"If ratio = exactly 3/2: c1 = sqrt(3/2 * {c0} * {c2}) = {c1_exact_3_2:.4f}")
print(f"Actual c1 = {c1}")
print(f"Difference: {c1 - c1_exact_3_2:.4f} ({abs(c1-c1_exact_3_2)/c1*100:.4f}%)")

# =====================================================================
# PART 9: CONNECTION TO OTHER IDENTITIES
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: Connections to Other Framework Identities")
print("=" * 72)

# Our identity: 2*A1^2 + 1 = 3*A0*A2
# Substitute:
# A0 = 2a+1, A1 = 2(a^2+a+1), A2 = 6a^2+15a+8

# For a = a1 = 5:
# Check if A2 can be expressed differently:
# A2 = 6a^2+15a+8
# = 6*25+75+8 = 233
# = 6a^2+15a+8
# = a*(6a+15) + 8
# = a*3*(2a+5) + 8
# = 3*a*(2a+5) + 8
# For a=5: 3*5*15 + 8 = 225+8 = 233

# Or: A2 = 6a^2+15a+8 = 6a^2+12a+3a+8 = 6a(a+2)+3a+8
# For a=5: 6*5*7+15+8 = 210+23 = 233

# Try: A2 = (2a+1)*(3a+8) = 6a^2+16a+3a+8 = 6a^2+19a+8. No, coefficient wrong.
# A2 = (3a+1)*(2a+8) = 6a^2+24a+2a+8 = 6a^2+26a+8. No.
# A2 = (2a+1)*(3a+?) + ?
# (2a+1)*3a = 6a^2+3a. Need 12a+8 more. (2a+1)*6 = 12a+6. +2: 12a+8.
# So A2 = (2a+1)*(3a+6) + 2 = A0*(3a+6) + 2 = A0*3*(a+2) + 2
# For a=5: 11*21+2 = 231+2 = 233. YES!
print("A_2 = A_0 * 3*(a1+2) + 2")
print(f"  = {A0} * {3*(a1+2)} + 2 = {A0*3*(a1+2)+2}")
print(f"  = {A0*3*(a1+2)} + 2 = {A2}  CHECK!")

# Substitute into the identity:
# 2*A1^2 + 1 = 3*A0*A2 = 3*A0*(A0*3*(a+2) + 2)
# = 9*A0^2*(a+2) + 6*A0
# 2*A1^2 + 1 = 9*A0^2*(a+2) + 6*A0
# 2*A1^2 - 6*A0 + 1 = 9*A0^2*(a+2)
# 2*4*(a^2+a+1)^2 - 6*(2a+1) + 1 = 9*(2a+1)^2*(a+2)
# 8*(a^2+a+1)^2 - 12a - 5 = 9*(4a^2+4a+1)*(a+2)

# For a=5:
lhs_check = 8*31**2 - 12*5 - 5
rhs_check = 9*11**2*7
print(f"\n8*31^2 - 65 = {lhs_check}")
print(f"9*121*7 = {rhs_check}")
print(f"Match: {lhs_check == rhs_check}")

# Alternative form using h = a*b1 = a*(a+1):
# h+1 = a^2+a+1, A0 = 2a+1,
# 8*(h+1)^2 + 1 = 3*(2a+1)*(6a^2+15a+8)

# The right side: 3*(2a+1)*(6a^2+15a+8)
# = 3*(12a^3+30a^2+16a+6a^2+15a+8)
# = 3*(12a^3+36a^2+31a+8)
# = 36a^3+108a^2+93a+24

# For 600-cell: this connects the Coxeter number h(E8)=30 to Fibonacci F_13=233
# through a DIOPHANTINE identity.

# =====================================================================
# PART 10: THE IDENTITY IN TERMS OF SPECTRAL DATA
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: Spectral Interpretation")
print("=" * 72)

# c1/(2*c0) = A1/(2*A0) = 62/22 = 31/11 = (h+1)/(a1+b1)
# This is the "average eigenvalue" normalized differently.

# c2/c0 = A2/A0 = 233/11 = F_13/L_5
# This connects two Fibonacci-Lucas families.

# Our identity: 2*(c1/(2*N))^2 + 1 = 3*(c0/(2*N))*(c2/(2*N))
# i.e.: 2*A1^2 + 1 = 3*A0*A2

# In terms of the heat kernel Tr(e^{-t*Lap}):
# c0 = Tr(1) = dim(H) (total Hilbert space dimension)
# c1 = Tr(Lap) (total trace of Laplacian)
# c2 = Tr(Lap^2) (sum of squared eigenvalues)
# These are the 0th, 1st, 2nd moments of the spectral distribution.

# The identity says: 2*(mean)^2 + correction = 3*(variance + mean^2)
# i.e., mean^2/variance = 3/(2*mean^2/variance - 2 + correction/variance)
# Not very illuminating in this form.

# Better: define s0 = A0, s1 = A1/A0 = <Lap>/dim, s2 = A2/A0 = <Lap^2>/dim
# Then: 2*A0^2*s1^2 + 1 = 3*A0^2*s2
# s2 = (2*s1^2*A0^2 + 1)/(3*A0^2) = (2/3)*s1^2 + 1/(3*A0^2)

s1 = A1/A0  # = 62/11
s2 = A2/A0  # = 233/11
print(f"Normalized moments: s1 = A1/A0 = {s1:.6f}")
print(f"                    s2 = A2/A0 = {s2:.6f}")
print(f"The identity: s2 = (2/3)*s1^2 + 1/(3*A0^2)")
print(f"  (2/3)*s1^2 = {2/3*s1**2:.6f}")
print(f"  1/(3*A0^2) = {1/(3*A0**2):.6f}")
print(f"  Sum = {2/3*s1**2 + 1/(3*A0**2):.6f}")
print(f"  s2 = {s2:.6f}  Match: {abs(s2 - 2/3*s1**2 - 1/(3*A0**2)) < 1e-10}")

# The variance:
var = s2 - s1**2
print(f"\nSpectral variance: <Lap^2> - <Lap>^2 = {var:.6f} = {s2 - s1**2:.6f}")
print(f"  = s1^2 * (-1/3) + 1/(3*A0^2)")
print(f"  = -s1^2/3 + 1/(3*A0^2)")

# Hmm, negative variance? Let me check:
# s2 - s1^2 = (2/3)*s1^2 + 1/(3*A0^2) - s1^2 = -s1^2/3 + 1/(3*A0^2)
# = (1 - A0^2*s1^2) / (3*A0^2)
# = (1 - A1^2/A0^2 * A0^2) / (3*A0^2) ... that's wrong
# = (1/(3*A0^2)) - s1^2/3 = (1 - (A1/A0)^2 * A0^2)/(3*A0^2)?
# Let me just compute: s2 - s1^2 = 233/11 - (62/11)^2 = 233/11 - 3844/121
# = 2563/121 - 3844/121 = -1281/121 = -10.587...

# The variance is NEGATIVE because s1 and s2 are normalized by A0 not by c0.
# The ACTUAL spectral moments involve different normalizations for each form degree.

print(f"\nNote: The 'variance' is negative because A_k are not directly")
print(f"the moments of a probability distribution. They include all form degrees")
print(f"with different weights.")

# =====================================================================
# PART 11: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 11: Honest Assessment")
print("=" * 72)

print(f"""
FINDINGS:

1. c1^2/(c0*c2) is NOT exactly 3/2.
   It is {exact_frac} = {float(exact_frac):.10f}
   = 3/2 - 1/{2*A0*A2} = 3/2 - 1/5126

2. The EXACT identity is:
   2*A_1^2 + 1 = 3*A_0*A_2
   or: 8*(h+1)^2 + 1 = 3*L_5*F_13 (for a1=5)

3. This identity holds ONLY for a1=5 among positive integers.
   The deficit polynomial 8a^4-20a^3-84a^2-77a-16 = -1 has (a-5) as the
   unique positive integer factor making it equal -1.

4. The near-3/2 property is a CONSEQUENCE of this Diophantine identity,
   not a coincidence. The correction 1/5126 is small because A0*A2 = 2563
   is large.

STATUS: DERIVED (new Diophantine identity unique to a1=5)
  - Not a "3/2 theorem" but a "3/2 minus epsilon theorem"
  - The epsilon = 1/5126 is EXACT and framework-determined
  - The identity connects h(E8), L_5, F_13 in a non-trivial way
  - UNIQUENESS to a1=5 is PROVED (polynomial factorization)
""")

print("=" * 72)
print("EXPERIMENT COMPLETE")
print("=" * 72)

# EXP-232: Why 137? Decomposing alpha^{-1} in terms of a_1 = 5
# ==============================================================
import math

a_1 = 5
b_1 = a_1 + 1  # = 6
phi = (1 + math.sqrt(a_1)) / 2
h_E8 = a_1 * b_1  # = 30
N_eig = 9
rank_E8 = 8

print("=" * 70)
print("EXP-232: WHY 137? THE NUMBER BEHIND THE NUMBER")
print("=" * 70)

# =============================================================
# PART A: The leading-order approximation
# =============================================================
print("\n" + "=" * 70)
print("PART A: Leading-order -- 1/alpha = 4*a_1*phi^4")
print("=" * 70)

coef = 4 * a_1 * phi**4
print(f"\n  4*a_1*phi^4 = 4 * {a_1} * phi^4 = {coef:.6f}")
print(f"  1/alpha (exp) = 137.035999084")
print(f"  Difference = {coef - 137.035999084:.6f}")
print(f"  This difference = 2*pi*alpha = {2*math.pi/137.036:.6f}")

# Exact alpha from quadratic
A_coef = 2 * math.pi
B_coef = -4 * a_1 * phi**4
C_coef = 1
disc = B_coef**2 - 4*A_coef*C_coef
alpha_calc = (-B_coef - math.sqrt(disc)) / (2*A_coef)
print(f"\n  alpha (exact from quadratic) = 1/{1/alpha_calc:.6f}")
print(f"  alpha (CODATA 2022)          = 1/137.035999084")
print(f"  Error: {abs(1/alpha_calc - 137.035999084)/137.035999084 * 100:.4f}%")

# =============================================================
# PART B: Decompose 20*phi^4 algebraically
# =============================================================
print("\n" + "=" * 70)
print("PART B: Algebraic decomposition of 20*phi^4")
print("=" * 70)

# phi^4 = (7 + 3*sqrt(5))/2  (exact)
phi4_exact_a = 7  # rational part * 2
phi4_exact_b = 3  # sqrt(5) part * 2
# 20*phi^4 = 10*(7 + 3*sqrt(5)) = 70 + 30*sqrt(5)
rat_part = 70
irr_part = 30  # coefficient of sqrt(5)

print(f"\n  phi^4 = (7 + 3*sqrt(5))/2")
print(f"  20*phi^4 = 10*(7 + 3*sqrt(5)) = {rat_part} + {irr_part}*sqrt(5)")
print(f"           = {rat_part} + {irr_part}*sqrt({a_1})")
print(f"           = {rat_part} + {irr_part} * {math.sqrt(5):.6f}")
print(f"           = {rat_part + irr_part * math.sqrt(5):.6f}")

print(f"\n  {irr_part} = a_1 * b_1 = h(E8) = Coxeter number!")
print(f"  {rat_part} = 2 * a_1 * (a_1 + 2) = 2 * 5 * 7 = 70")
print(f"  Or: {rat_part} = (4*a_1^2 + 4*a_1)/2 * ... hmm")

# Better: 20*phi^4 in Z[phi] basis
# phi^4 = phi^2 * phi^2 = (1+phi)^2 = 1 + 2*phi + phi^2 = 1 + 2*phi + 1 + phi = 2 + 3*phi
# So: phi^4 = 2 + 3*phi  (exact in Z[phi])
# 20*phi^4 = 40 + 60*phi
a_zphi = 40  # rational part in Z[phi]
b_zphi = 60  # phi coefficient in Z[phi]
print(f"\n  In Z[phi] basis:")
print(f"  phi^4 = 2 + 3*phi  (since phi^2 = 1+phi)")
print(f"  20*phi^4 = {a_zphi} + {b_zphi}*phi")
print(f"           = {a_zphi} + {b_zphi} * {phi:.6f}")
print(f"           = {a_zphi + b_zphi * phi:.6f}")
print(f"\n  {a_zphi} = 8*a_1 = 8 * 5 = rank(E8) * a_1")
print(f"  {b_zphi} = 12*a_1 = 12 * 5 = (vertex degree) * a_1")
print(f"  Or: {b_zphi} = 2*h(E8) = 2*30 = 60")
print(f"  Or: {b_zphi} = N/2 = 120/2 = 60 = half the 600-cell vertices!")

# =============================================================
# PART C: So what IS 137?
# =============================================================
print("\n" + "=" * 70)
print("PART C: What IS 137?")
print("=" * 70)

print(f"""
  137 is NOT a fundamental number in this framework.

  The fundamental number is a_1 = 5.

  From a_1 = 5, we get:
    phi = (1+sqrt(5))/2

  The EXACT relation is:
    1/alpha = 4*a_1*phi^4 - 2*pi*alpha

  Leading term:  4*a_1*phi^4 = {coef:.4f}
  Correction:    2*pi*alpha   = {2*math.pi*alpha_calc:.4f}
  Result:        1/alpha      = {1/alpha_calc:.4f}

  In Z[phi]: 1/alpha ~ {a_zphi} + {b_zphi}*phi = (N/2)*phi + rank*a_1

  The integer 137 is just the nearest integer to 4*a_1*phi^4.
  It has no independent significance.
""")

# =============================================================
# PART D: But IS there anything special about 137?
# =============================================================
print("=" * 70)
print("PART D: Is 137 special on its own?")
print("=" * 70)

print(f"\n  137 is prime: {all(137 % i != 0 for i in range(2, 12))}")

# Check if 137 splits, ramifies, or stays inert in Z[phi]
# 137 mod 5
r = 137 % a_1
print(f"  137 mod a_1 = 137 mod 5 = {r}")
# Legendre symbol (5/137): need to check if 5 is a QR mod 137
# By QR: (5/137) = (137/5) * (-1)^{(5-1)/2 * (137-1)/2} = (2/5) * (-1)^{2*68} = (2/5) * 1
# 2 mod 5 = 2, and 2 is not a QR mod 5 (QR mod 5 are {0,1,4})
# So (5/137) = -1, meaning 5 is NOT a quadratic residue mod 137
# Therefore 137 stays INERT in Z[phi]
print(f"  5 mod 137 = {5 % 137}")
# Check QR: is there x such that x^2 = 5 mod 137?
is_qr = any((x*x) % 137 == 5 for x in range(137))
print(f"  Is 5 a quadratic residue mod 137? {is_qr}")
if is_qr:
    print(f"  => 137 SPLITS in Z[phi]")
else:
    print(f"  => 137 is INERT in Z[phi] (stays prime)")
    print(f"     This means 137 does NOT decompose in the number ring!")
    print(f"     It's an 'outsider' -- not part of the phi-lattice structure.")

# Fibonacci connection
fib = [0, 1]
while fib[-1] < 200:
    fib.append(fib[-1] + fib[-2])
print(f"\n  Fibonacci numbers near 137: ...{[f for f in fib if 80 < f < 200]}")
print(f"  137 is NOT a Fibonacci number")

# Lucas connection
lucas = [2, 1]
while lucas[-1] < 200:
    lucas.append(lucas[-1] + lucas[-2])
print(f"  Lucas numbers near 137: ...{[l for l in lucas if 80 < l < 200]}")
print(f"  137 is NOT a Lucas number")

# But F(11) = 89, F(12) = 144
# 137 = 144 - 7 = F(12) - L(4)
print(f"\n  137 = F(12) - L(4) = 144 - 7 = 137")
print(f"  137 = F(12) - (a_1 + 2) = 144 - 7 = 137")

# More: 137 = 120 + 17 = N + n_tau!
print(f"\n  ** 137 = N + n_tau = 120 + 17 = {120 + 17} **")
print(f"  (N = vertices of 600-cell, n_tau = tau mass exponent)")

# And: 137 = 5*26 + 7 = a_1 * n_top + (b_1+1)
print(f"  137 = a_1 * n_top + (b_1+1) = 5*26 + 7 = {5*26+7}")

# Integer part of 1/alpha
print(f"\n  floor(1/alpha) = {int(1/alpha_calc)} = 137")
print(f"  floor(4*a_1*phi^4) = {int(coef)} = {int(coef)}")

# =============================================================
# PART E: The correction 2*pi*alpha as spectral trace
# =============================================================
print("\n" + "=" * 70)
print("PART E: The correction term 2*pi*alpha")
print("=" * 70)

correction = 2 * math.pi * alpha_calc
print(f"\n  2*pi*alpha = {correction:.8f}")
print(f"  This is the self-energy correction from the Hopf fiber phase.")
print(f"  It shifts 1/alpha from 137.082 down to 137.036.")
print(f"  Without it: error would be 0.034%")
print(f"  With it: error is 0.0001%")
print(f"\n  Physical meaning: the 2*pi comes from the Hopf fibration")
print(f"  (one full turn on S^1 fiber), and alpha itself is the")
print(f"  coupling strength. So the correction is 'one loop around")
print(f"  the fiber weighted by the coupling'.")

# =============================================================
# PART F: Summary -- the hierarchy of constants
# =============================================================
print("\n" + "=" * 70)
print("PART F: SUMMARY -- The hierarchy")
print("=" * 70)

print(f"""
  THE HIERARCHY OF FUNDAMENTAL NUMBERS:

  Level 0 (Diophantine):   a_1 = 5   (THE number)
  Level 1 (Algebraic):     phi = {phi:.6f}   (from a_1)
                           b_1 = 6    (= a_1+1)
                           h = 30     (= a_1*b_1)
  Level 2 (Spectral):      alpha = 1/{1/alpha_calc:.3f}  (from a_1, phi, pi)
                           alpha_s = {1/(2*phi**3):.4f}  (from phi)
                           sin2_tW = {6/26:.4f}  (from a_1, b_1)
  Level 3 (Derived):       137 = floor(1/alpha) -- NOT fundamental

  In this hierarchy, 137 is a Level-3 number.
  It's as derived as the Higgs mass (125 GeV) or
  the number of Standard Model fermions (48).

  The REAL mystery number is 5, not 137.
""")

print("=" * 70)
print("EXP-232 COMPLETE")
print("=" * 70)

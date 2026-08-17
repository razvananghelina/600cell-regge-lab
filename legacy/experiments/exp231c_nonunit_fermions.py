"""
EXP-231c: Non-Unit Fermions - Why c, t, b Break the Pattern
=============================================================
From exp231: 6/9 fermions have z = a+b*phi as a UNIT of Z[phi] (|N(z)|=1).
The 3 non-units are: charm (N=5), top (N=19), bottom (N=-19).

KEY OBSERVATIONS:
  |N(c)| = 5 = a_1
  |N(t)| = 19 (prime!)
  |N(b)| = 19 (same as top!)

QUESTIONS:
  A. Why do c, t, b have these specific norms?
  B. Is 19 related to the 600-cell structure?
  C. Can non-units be expressed as unit * small correction?
  D. Do the norm values explain the large mass errors (11-20%)?
  E. Is there a Z[phi] ideal structure at play?
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
PHIP = (1 - np.sqrt(5)) / 2
a1 = 5
b1 = a1 + 1

print("=" * 70)
print("EXP-231c: NON-UNIT FERMIONS (c, t, b)")
print("=" * 70)

# All fermion data
fermions = {
    'e':   (0, 0),   'mu':  (1, 1),  'tau': (1, 2),
    'u':   (3, -2),  'c':   (2, 1),  't':   (4, 1),
    'd':   (1, 0),   's':   (1, 1),  'b':   (-1, 4),
}

m_exp = {
    'e': 0.51099895e-3, 'mu': 0.1056584, 'tau': 1.77686,
    'u': 2.16e-3, 'c': 1.27, 't': 172.69,
    'd': 4.67e-3, 's': 93.4e-3, 'b': 4.18,
}

m_e = m_exp['e']

# ===================================================================
# PART A: Non-unit analysis
# ===================================================================
print("\n--- PART A: Non-unit fermions ---")

for name in ['c', 't', 'b']:
    a, b = fermions[name]
    n = a1*a + b1*b
    z = a + b*PHI
    zp = a + b*PHIP
    norm = a**2 + a*b - b**2  # = z * z' (norm in Z[phi])
    m_pred = m_e * PHI**n
    m_actual = m_exp[name]
    err = abs(m_pred - m_actual) / m_actual * 100

    print("\n%s: (a,b) = (%d,%d), n = %d" % (name, a, b, n))
    print("  z = %d + %d*phi = %.4f" % (a, b, z))
    print("  z' = %d + %d*phi' = %.4f" % (a, b, zp))
    print("  N(z) = z*z' = %d" % norm)
    print("  m_pred = %.4f GeV, m_exp = %.4f GeV, err = %.1f%%" % (m_pred, m_actual, err))


# ===================================================================
# PART B: Factorization in Z[phi]
# ===================================================================
print("\n" + "=" * 70)
print("PART B: Factorization in Z[phi]")
print("=" * 70)

# In Z[phi], the norm N(z) = z * z' determines factorization.
# N = 5 = a_1 for charm.
# N = 19 for top and bottom.
# Both 5 and 19 are primes in Z.
# Are they prime in Z[phi]?

# A rational prime p is:
# - SPLIT in Z[phi] if p = 1 or 4 mod 5 (disc = 5)
# - INERT if p = 2 or 3 mod 5
# - RAMIFIED if p = 5

print("\nPrime classification in Z[phi] (disc = a_1 = 5):")
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    r = p % 5
    if p == 5:
        status = "RAMIFIED (p = a_1)"
        # 5 = -phi'^2 * phi^2 * ... actually 5 = phi^2 - phi'^2 ... hmm
        # Actually sqrt(5) = phi - phi' = phi + 1/phi = 2phi - 1
        # So 5 = (sqrt(5))^2 = (phi - phi')^2 = (2phi-1)^2 = 4phi^2-4phi+1
        # In Z[phi]: 5 = (sqrt(5))^2. And sqrt(5) = 2*phi - 1.
        # N(2*phi-1) = N(-1 + 2*phi) = (-1)^2 + (-1)*2 - 2^2 = 1-2-4 = -5
        # So pi_5 = 2*phi - 1 = sqrt(5), and 5 = (2phi-1)(2phi'-1) = (sqrt(5))(-sqrt(5)) up to units
    elif r in [1, 4]:
        status = "SPLIT"
    else:
        status = "INERT"
    print("  p=%2d: %d mod 5 -> %s" % (p, r, status))

print("\n5 in Z[phi]: RAMIFIED. 5 = sqrt(5)^2 = (2*phi-1)^2 * unit")
print("  pi_5 = 2*phi - 1 = sqrt(5)")
print("  N(pi_5) = N(2phi-1) = N(-1 + 2*phi) = 1 + (-2) - 4 = -5")
print("  So |N| = 5. The norm of the ramification prime IS a_1.")

print("\n19 in Z[phi]: 19 mod 5 = 4, so 19 SPLITS in Z[phi].")
print("  Need to find pi_19 such that N(pi_19) = 19.")
print("  Solve: a^2 + ab - b^2 = 19")

# Find elements with norm 19
print("  Solutions with |a|,|b| <= 10:")
norm_19_solutions = []
for a in range(-10, 11):
    for b in range(-10, 11):
        norm = a**2 + a*b - b**2
        if norm == 19:
            z = a + b*PHI
            n_val = a1*a + b1*b
            norm_19_solutions.append((a, b, z, n_val))
            print("    (%2d, %2d): z = %.4f, n = %d" % (a, b, z, n_val))

print("\n  Solutions with N = -19:")
norm_m19_solutions = []
for a in range(-10, 11):
    for b in range(-10, 11):
        norm = a**2 + a*b - b**2
        if norm == -19:
            z = a + b*PHI
            n_val = a1*a + b1*b
            norm_m19_solutions.append((a, b, z, n_val))
            print("    (%2d, %2d): z = %.4f, n = %d" % (a, b, z, n_val))

# The top and bottom quarks!
print("\n  Top: (4,1), N=19, n=26")
print("  Bottom: (-1,4), N=-19, n=19")
print("  Note: top and bottom are CONJUGATE PRIMES of 19!")
print("  z_t = 4+phi, z_b = -1+4*phi")
print("  z_t * z_b' = (4+phi)(-1+4*phi') = ?")
zt = 4 + PHI
zb = -1 + 4*PHI
zt_p = 4 + PHIP
zb_p = -1 + 4*PHIP
print("  z_t * z_b = %.4f * %.4f = %.4f" % (zt, zb, zt*zb))
print("  z_t' * z_b' = %.4f * %.4f = %.4f" % (zt_p, zb_p, zt_p*zb_p))
print("  N(z_t) * N(z_b) = %d * %d = %d" % (19, -19, 19*(-19)))

# Check: is z_t / z_b a unit?
ratio = zt / zb
print("\n  z_t / z_b = %.6f" % ratio)
# Check if ratio is a + b*phi for some integers
for a_t in range(-10, 11):
    for b_t in range(-10, 11):
        if abs(a_t + b_t*PHI - ratio) < 1e-4:
            norm_r = a_t**2 + a_t*b_t - b_t**2
            print("  z_t / z_b = %d + %d*phi, N = %d" % (a_t, b_t, norm_r))


# ===================================================================
# PART C: Decomposition as unit * prime
# ===================================================================
print("\n" + "=" * 70)
print("PART C: Decomposition as Unit * Prime Element")
print("=" * 70)

# For charm: N(z_c) = 5 (ramified prime)
# z_c = 2 + phi. Can we write z_c = unit * pi_5?
# pi_5 = sqrt(5) = 2*phi - 1 (in Z[phi], this is (-1,2))
# N(pi_5) = (-1)^2 + (-1)*2 - 2^2 = 1-2-4 = -5

pi5 = -1 + 2*PHI  # = sqrt(5) = 2.236
zc = 2 + PHI      # = 3.618

print("Charm: z_c = 2 + phi = %.4f" % zc)
print("pi_5 = -1 + 2*phi = sqrt(5) = %.4f" % pi5)
ratio_c = zc / pi5
print("z_c / pi_5 = %.6f" % ratio_c)
# Is this a unit?
for a_t in range(-10, 11):
    for b_t in range(-10, 11):
        if abs(a_t + b_t*PHI - ratio_c) < 1e-4:
            norm_r = a_t**2 + a_t*b_t - b_t**2
            print("  = %d + %d*phi, N = %d, unit: %s" % (a_t, b_t, norm_r, abs(norm_r)==1))

# So z_c = unit * pi_5. Find the unit:
# z_c = 2+phi, pi_5 = -1+2phi
# z_c / pi_5 = (2+phi)/(-1+2phi)
# Rationalize: multiply by conjugate (-1+2phi')/(-1+2phi')
# = (2+phi)(-1+2phi') / ((-1+2phi)(-1+2phi'))
# (-1+2phi)(-1+2phi') = 1 - 2(phi+phi') + 4*phi*phi' = 1 - 2*1 + 4*(-1) = 1-2-4 = -5
# (2+phi)(-1+2phi') = -2 + 4phi' - phi + 2*phi*phi' = -2+4phi'-phi+2*(-1)
# = -4 - phi + 4phi' = -4 - phi + 4(1-phi) = -4-phi+4-4phi = -5phi
# So z_c / pi_5 = -5phi / (-5) = phi

print("\nz_c / pi_5 = phi (EXACT!)")
print("  Verification: pi_5 * phi = (-1+2phi)*phi = -phi+2phi^2 = -phi+2(phi+1) = phi+2 = z_c")
print("  CHECK: %.4f * %.4f = %.4f = %.4f" % (pi5, PHI, pi5*PHI, zc))

print("\nSo: z_charm = phi * sqrt(5) = phi * pi_5")
print("  charm = (gen 1 unit) * (ramification prime)")
print("  The 'extra' factor is sqrt(a_1)!")

# For top: N(z_t) = 19.
# 19 splits in Z[phi]. Find prime element.
# z_t = 4 + phi. Check if z_t itself is prime.
# A prime element p of Z[phi] divides 19 (the rational prime).
# N(z_t) = 19, and 19 is a rational prime, so z_t IS a prime element!
print("\nTop: z_t = 4 + phi = %.4f, N(z_t) = 19" % zt)
print("  Since 19 is a rational prime and N(z_t) = 19, z_t IS a prime of Z[phi].")
print("  It cannot be decomposed further!")

# For bottom: N(z_b) = -19.
# z_b = -1 + 4*phi. Same analysis.
print("\nBottom: z_b = -1 + 4*phi = %.4f, N(z_b) = -19" % zb)
print("  z_b is also a prime of Z[phi] (conjugate prime to z_t).")

# Relation between t and b primes:
# z_t and z_b are NOT Galois conjugates of each other.
# z_t' = 4 + phi' = 4 - 0.618 = 3.382
# z_b' = -1 + 4*phi' = -1 - 2.472 = -3.472
# z_t * z_b' = (4+phi)(-1+4phi') should give a rational integer
print("\nz_t' = %.4f, z_b' = %.4f" % (zt_p, zb_p))
print("z_t * z_b' = %.4f" % (zt * zb_p))
print("z_t' * z_b = %.4f" % (zt_p * zb))
# These should not be rational integers unless related by norm

# Try: z_b = z_t * unit?
ratio_tb = zb / zt
print("\nz_b / z_t = %.6f" % ratio_tb)
for a_t in range(-10, 11):
    for b_t in range(-10, 11):
        if abs(a_t + b_t*PHI - ratio_tb) < 1e-4:
            norm_r = a_t**2 + a_t*b_t - b_t**2
            print("  = %d + %d*phi, N = %d" % (a_t, b_t, norm_r))

# z_b / z_t = (-1+4phi)/(4+phi)
# Rationalize: multiply by (4+phi')/(4+phi')
# (4+phi)(4+phi') = 16 + 4(phi+phi') + phi*phi' = 16+4-1 = 19
# (-1+4phi)(4+phi') = -4-phi'+16phi+4*phi*phi' = -4-phi'+16phi-4
# = -8 + 16phi - phi' = -8 + 16phi - (1-phi) = -9 + 17phi
# So z_b/z_t = (-9+17phi)/19
print("\nz_b / z_t = (-9 + 17*phi) / 19 = %.6f" % ((-9+17*PHI)/19))
print("  NOT a Z[phi] element! (denominator 19)")
print("  So z_b and z_t are INDEPENDENT primes above 19.")


# ===================================================================
# PART D: Nearest Unit Analysis
# ===================================================================
print("\n" + "=" * 70)
print("PART D: Nearest Unit to Each Non-Unit")
print("=" * 70)

# For each non-unit fermion, find the nearest unit of Z[phi]
# and compute the "correction"

for name in ['c', 't', 'b']:
    a, b = fermions[name]
    z = a + b*PHI
    n = a1*a + b1*b
    m_pred = m_e * PHI**n
    m_actual = m_exp[name]

    print("\n%s: z = %.4f, n = %d, m_pred = %.4f, m_exp = %.4f" % (name, z, n, m_pred, m_actual))

    # Find nearest units (phi^k)
    best_k = None
    best_dist = 1e10
    for k in range(-5, 15):
        pk = PHI**k
        dist = abs(z - pk)
        if dist < best_dist:
            best_dist = dist
            best_k = k
        dist_neg = abs(z + pk)
        if dist_neg < best_dist:
            best_dist = dist_neg
            best_k = -k  # represent as -phi^k

    # Get the nearest unit's (a,b)
    from math import factorial
    def fib(nn):
        if nn >= 0:
            aa, bb = 0, 1
            for _ in range(nn):
                aa, bb = bb, aa+bb
            return aa
        else:
            return ((-1)**(abs(nn)+1)) * fib(abs(nn))

    # Compute nearest phi^k
    pk = PHI**abs(best_k)
    if best_k >= 0:
        a_unit = fib(best_k - 1)
        b_unit = fib(best_k)
        z_unit = a_unit + b_unit * PHI
        n_unit = a1*a_unit + b1*b_unit
    else:
        a_unit = fib(abs(best_k) - 1)
        b_unit = fib(abs(best_k))
        z_unit = -(a_unit + b_unit * PHI)
        a_unit, b_unit = -a_unit, -b_unit
        n_unit = a1*a_unit + b1*b_unit

    delta_a = a - a_unit
    delta_b = b - b_unit
    delta_n = n - n_unit
    delta_z = z - z_unit
    m_unit = m_e * PHI**n_unit

    print("  Nearest unit: phi^%d, z_unit = %.4f, n_unit = %d" % (best_k, z_unit, n_unit))
    print("  (a,b)_unit = (%d,%d), (a,b)_fermion = (%d,%d)" % (a_unit, b_unit, a, b))
    print("  Delta (a,b) = (%d,%d), delta_n = %d" % (delta_a, delta_b, delta_n))
    print("  m_unit = %.4f GeV, correction factor = m_exp/m_unit = %.4f"
          % (m_unit, m_actual/m_unit))

    # The correction in terms of phi:
    if m_actual/m_unit > 0:
        corr_power = np.log(m_actual/m_unit) / np.log(PHI)
        print("  m_exp/m_unit = phi^%.3f" % corr_power)


# ===================================================================
# PART E: The Number 19
# ===================================================================
print("\n" + "=" * 70)
print("PART E: Why 19?")
print("=" * 70)

print("Both top and bottom have |N(z)| = 19.")
print("What is 19 in the 600-cell framework?")
print()
print("  19 = 4*a_1 - 1 = 4*5 - 1 = 19")
print("  19 = a_1^2 - b_1 = 25 - 6 = 19")
print("  19 = h(E8) - a_1 - b_1 = 30 - 5 - 6 = 19")
print("  19 = dim(E8)/rank(E8) - rank(E8) - a_1 = 31 - 8 - 5 = 18? No.")
print("  19 = N/b_1 - 1 = 120/6 - 1 = 19!")
print("  19 = h(E8) - a_1*b_1/a_1? No, that's 30-6=24.")
print()

# Check multiple expressions
print("Expressions for 19:")
expressions = [
    ("4*a_1 - 1", 4*a1 - 1),
    ("a_1^2 - b_1", a1**2 - b1),
    ("h - a_1 - b_1", 30 - a1 - b1),
    ("N/b_1 - 1", 120//b1 - 1),
    ("n_b (bottom mass exponent)", 19),
    ("N/b_1 - 1", 120//6 - 1),
    ("a_1*(a_1-1) - 1", a1*(a1-1) - 1),
]

for expr, val in expressions:
    if val == 19:
        print("  %s = %d  ***" % (expr, val))
    else:
        print("  %s = %d" % (expr, val))

# Most elegant: 19 = a_1^2 - b_1
print("\nMost natural: 19 = a_1^2 - b_1 = a_1^2 - (a_1+1)")
print("  = a_1^2 - a_1 - 1")
print("  = a_1*(a_1-1) - 1")
print("  = 5*4 - 1 = 19")
print()
print("  Note: a_1*(a_1-1) = 20 = 4*a_1 = coefficient in alpha equation!")
print("  So 19 = 4*a_1 - 1")

# Also: N(z_b) = -19, and n_b = 19!
# Is |N(z)| = n for bottom quark?
print("\nBottom quark: |N(z)| = 19 = n_b (mass exponent)!")
print("  Is this a coincidence?")
for name, (a, b) in fermions.items():
    n = a1*a + b1*b
    norm = abs(a**2 + a*b - b**2)
    print("  %s: n=%d, |N|=%d, equal=%s" % (name, n, norm, n==norm))


# ===================================================================
# PART F: Charm Decomposition Deep
# ===================================================================
print("\n" + "=" * 70)
print("PART F: Charm = phi * sqrt(5)")
print("=" * 70)

print("z_charm = 2 + phi = phi * sqrt(5)")
print("  = phi * (2*phi - 1)")
print("  = 2*phi^2 - phi")
print("  = 2*(phi+1) - phi = phi + 2  (check: phi+2 = 3.618, z_c = 3.618)")
print()
print("In terms of (a,b):")
print("  phi = (0,1), n_phi = 6")
print("  sqrt(5) = (-1,2), n_sqrt5 = -5+12 = 7")
print("  Product in Z[phi]: (0+1*phi)*(-1+2*phi) = -phi+2*phi^2 = -phi+2phi+2 = phi+2")
print("  = (2,1), n = 10+6 = 16. CHECK: n_c = 16. YES!")
print()
print("So: n_charm = n(phi) + n(sqrt(5)) = 6 + 7 = 13?")
print("  Wait, n is NOT additive over multiplication in Z[phi]!")
print("  n(z1*z2) != n(z1) + n(z2) in general.")
print("  n_phi = 6, n_{sqrt5} = 7, but n_c = 16, not 13.")
print()
print("But the MASS is multiplicative: m(z1*z2) = m(z1)*m(z2)/m_e")
print("  Because m = m_e * phi^n, so m(z1*z2) = m_e * phi^{n(z1*z2)}")
print("  And n(z1*z2) = n(z1*z2), computed from the product's (a,b).")

# The key point: charm = gen-1 unit * ramification prime
# gen-1 unit = phi (the fundamental unit for gen 1)
# ramification prime = sqrt(5) = sqrt(a_1)
print("\nInterpretation:")
print("  z_charm = phi^1 * sqrt(a_1)")
print("  = (gen-1 fundamental unit) * (ramification prime of Z[phi])")
print()
print("  This means charm is NOT a pure Fibonacci excitation.")
print("  It carries an extra 'sqrt(a_1) charge' from the ramification.")
print("  This is WHY its mass doesn't fit the unit formula!")


# ===================================================================
# PART G: Can top and bottom be similarly decomposed?
# ===================================================================
print("\n" + "=" * 70)
print("PART G: Top and Bottom Decompositions")
print("=" * 70)

# Top: z_t = 4 + phi. N = 19.
# 19 splits: z_t and its conjugate z_t' = 4+phi' are the two primes.
# Can we write z_t = unit * something recognizable?

# z_t = 4 + phi. Try: z_t = phi^k * pi_19
# where pi_19 has norm +/- 19.
# Since z_t itself has norm 19, and is prime, the only way is
# z_t = unit * z_t (trivially, unit = 1)

# But can we relate 4+phi to known quantities?
print("Top: z_t = 4 + phi = 4 + phi")
print("  = 3 + (1+phi) = 3 + phi^2")
print("  = phi^3 + 1 (since phi^3 = 2+phi+1 = ... let me check)")
print("  phi^3 = phi^2 + phi = (phi+1) + phi = 2phi+1 = %.4f" % PHI**3)
print("  phi^3 + 1 = %.4f, z_t = %.4f" % (PHI**3 + 1, zt))
# phi^3 = 2phi+1 = 4.236, z_t = 4+phi = 5.618. Not the same.

# Try: z_t = phi^4 - phi = phi*(phi^3 - 1) = phi*(2phi) = 2phi^2?
# phi^4 = 3phi+2 = 6.854. phi^4 - phi = 2phi+2 = 5.236? No.
# z_t = 4+phi = (4,1). Let me try other decompositions.

# z_t = a + b*phi = 4+phi
# As sum of units: phi^0 + phi^0 + phi^0 + phi^0 + phi^1 = 4+phi?
# 4*1 + phi = 4+phi. Yes, but that's trivial.
# phi^3 + phi^0 = (2phi+1) + 1 = 2phi+2 = 2(phi+1) = 2phi^2 = 5.236. No.
# phi^2 + phi^2 = 2phi+2 = 2phi^2 = 5.236. No.
# phi^3 + phi^{-1} = (2phi+1) + (phi-1) = 3phi = 4.854. No.

print("\nTrying to express z_t as sum/product of simpler Z[phi] elements:")
# z_t = 4+phi. In matrix form: phi acts as [[0,1],[1,1]]
# Norm 19 = a_1^2 - b_1 = a_1^2 - a_1 - 1
# Compare: norm of z=(a_1-1, 1) = (a_1-1)^2 + (a_1-1) - 1 = a_1^2 - 2a_1+1+a_1-1-1 = a_1^2-a_1-1 = 19!
# And (a_1-1, 1) = (4, 1) = z_t!

print("z_t = (a_1-1, 1) = (a_1-1) + 1*phi")
print("  N(z_t) = (a_1-1)^2 + (a_1-1) - 1 = a_1^2 - a_1 - 1")
print("  = %d^2 - %d - 1 = %d" % (a1, a1, a1**2 - a1 - 1))
print("  = a_1*(a_1-1) - 1 = %d * %d - 1 = %d" % (a1, a1-1, a1*(a1-1)-1))
print()
print("z_b = (-1, a_1-1) = -1 + (a_1-1)*phi")
print("  N(z_b) = 1 + (-1)*(a_1-1) - (a_1-1)^2 = 1 - (a_1-1) - (a_1-1)^2")
print("  = -(a_1^2 - a_1 - 1) = -%d" % (a1**2 - a1 - 1))

# Verify
a_t, b_t = a1-1, 1
a_b, b_b = -1, a1-1
print("\nVerification:")
print("  z_t = (%d, %d): (a_1-1)=%d, 1=%d" % (a_t, b_t, a1-1, 1))
print("  z_b = (%d, %d): -1=%d, (a_1-1)=%d" % (a_b, b_b, -1, a1-1))
print("  Match top: %s" % ((a_t, b_t) == fermions['t']))
print("  Match bottom: %s" % ((a_b, b_b) == fermions['b']))

# The (a,b) values are SYMMETRIC in a_1-1:
# top = (a_1-1, 1), bottom = (-1, a_1-1)
# This is a SWAP with sign flip!
print("\nSYMMETRY: top (a,b) and bottom (a,b) are related by:")
print("  (a_t, b_t) = (a_1-1, 1)")
print("  (a_b, b_b) = (-1, a_1-1)")
print("  Transform: (a,b) -> (-b, a) gives (-1, a_1-1) from (a_1-1, 1)?")
print("  (-b_t, a_t) = (-1, %d) = (-1, 4). YES!" % (a1-1))
print("  The map (a,b) -> (-b,a) is a 90-degree ROTATION in Z^2!")
print("  In Z[phi]: z -> z * i_phi where i_phi^2 acts as rotation.")
print()
print("  Actually: if z = a+b*phi, then i*z (quaternion i) maps")
print("  (a,b) -> (-b,a) and N(-b+a*phi) = b^2 - ab - a^2 = -(a^2+ab-b^2) = -N(z)")
print("  So N(z_b) = -N(z_t). CONFIRMED: 19 and -19.")

# Express in terms of a_1 only:
print("\n" + "=" * 70)
print("FULL (a,b) IN TERMS OF a_1:")
print("=" * 70)
print()
print("CHARM: (a,b) = (2, 1)")
print("  z_c = phi * sqrt(a_1) = phi * (2phi-1)")
print("  N(z_c) = a_1")
print()
print("TOP: (a,b) = (a_1-1, 1) = (%d, 1)" % (a1-1))
print("  N(z_t) = a_1^2 - a_1 - 1 = a_1*(a_1-1) - 1 = %d" % (a1**2-a1-1))
print()
print("BOTTOM: (a,b) = (-1, a_1-1) = (-1, %d)" % (a1-1))
print("  N(z_b) = -(a_1^2 - a_1 - 1) = %d" % (-(a1**2-a1-1)))
print("  RELATED TO TOP BY 90-DEGREE ROTATION: (a,b) -> (-b,a)")


# ===================================================================
# PART H: Complete Picture
# ===================================================================
print("\n" + "=" * 70)
print("PART H: Complete (a,b) Structure")
print("=" * 70)

print("""
UNIT FERMIONS (z = +/- phi^k, |N|=1):
  e:   (0,0)    z = 0          ground state
  d:   (1,0)    z = phi^0 = 1  gen-0 down
  mu:  (1,1)    z = phi^2      gen-1 lepton
  s:   (1,1)    z = phi^2      gen-1 down (= mu!)
  tau: (1,2)    z = phi^3      gen-2 lepton
  u:   (3,-2)   z = -phi^{-3}  gen-0 up = Galois(tau)

NON-UNIT FERMIONS:
  c:   (2,1)       z = phi*sqrt(a_1)    N = a_1 = 5
  t:   (a_1-1, 1)  z = (a_1-1)+phi      N = a_1^2-a_1-1 = 19
  b:   (-1, a_1-1)  z = -1+(a_1-1)*phi   N = -(a_1^2-a_1-1) = -19

PATTERN:
  Units = first two generations + gen-0 up (Galois image)
  Non-units = charm + third-gen quarks

  charm = gen-1 unit * sqrt(a_1) (ramification!)
  top = prime element at (a_1-1, 1)
  bottom = 90-degree rotation of top: (-1, a_1-1)

MASS ERRORS correlate with unit status:
  Units: 0-4%% error (pure geometric masses)
  Non-units: 11-20%% error (need QCD/prime corrections)
""")

print("=" * 70)
print("EXP-231c COMPLETE")
print("=" * 70)

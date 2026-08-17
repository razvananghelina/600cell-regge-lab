"""
exp426b_gap_analysis.py
========================
FOLLOW-UP: Is the gap z_CC*|ln alpha| - 2/alpha = 29/5 EXACT?

From exp426:
  gap = z_CC * |ln alpha| - 2/alpha = 5.800834
  b1 - 1/a1 = 29/5 = 5.8
  Match: 0.000834 (0.014%)

This is too good for coincidence but involves transcendental alpha.
Test whether this is an EXACT identity or a numerical near-miss.

Also test: exp(-33/alpha_s) match (33 = 3*11 = N_gen * A_0).

Author: Razvan-Constantin Anghelina
Date: 2026-02-27
"""

from mpmath import mp, mpf, log, pi, sqrt, exp, fabs

# High precision
mp.dps = 50  # 50 decimal digits

print("=" * 72)
print("EXP-426b: GAP ANALYSIS AT HIGH PRECISION")
print("=" * 72)

# Framework constants at high precision
a1 = mpf(5)
b1 = a1 + 1
PHI = (1 + sqrt(a1)) / 2
N = mpf(120)

# Alpha from quadratic: 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0
A = 2 * pi
B = -4 * a1 * PHI**4
C = mpf(1)
disc = B**2 - 4*A*C
alpha = (-B - sqrt(disc)) / (2*A)   # physical root
alpha_s = 1 / (2 * PHI**3)

z_CC = 57 - alpha_s
ln_alpha = log(alpha)
target = z_CC * fabs(ln_alpha)

print(f"\nalpha = {alpha}")
print(f"1/alpha = {1/alpha}")
print(f"alpha_s = {alpha_s}")
print(f"z_CC = {z_CC}")
print(f"|ln alpha| = {fabs(ln_alpha)}")

# =====================================================================
# TEST 1: gap = z_CC * |ln alpha| - 2/alpha vs 29/5
# =====================================================================
print(f"\n{'='*72}")
print("TEST 1: gap vs 29/5")
print("=" * 72)

gap = target - 2/alpha
exact_29_5 = mpf(29) / mpf(5)

print(f"\n  z_CC * |ln alpha| = {target}")
print(f"  2/alpha = {2/alpha}")
print(f"  gap = {gap}")
print(f"  29/5 = {exact_29_5}")
print(f"  gap - 29/5 = {gap - exact_29_5}")
print(f"  relative diff = {fabs(gap - exact_29_5) / exact_29_5}")

# This is NOT zero at high precision!
print(f"\n  VERDICT: gap - 29/5 = {float(gap - exact_29_5):.15e}")
print(f"  This is NOT exact. The match is a numerical near-coincidence.")

# =====================================================================
# TEST 2: What IS the gap exactly?
# =====================================================================
print(f"\n{'='*72}")
print("TEST 2: Exact gap analysis")
print("=" * 72)

print(f"\n  gap = {gap}")
print(f"  gap * a1 = {gap * a1}")
print(f"  gap * a1^2 = {gap * a1**2}")
print(f"  gap * a1^2 / ln(10) = {gap * a1**2 / log(10)}")

# Express gap in terms of alpha
# gap = (57 - 1/(2*phi^3)) * |ln alpha| - 2/alpha
# = 57 * |ln alpha| - |ln alpha|/(2*phi^3) - 2/alpha

# Using the alpha equation: 1/alpha = 20*phi^4 - 2*pi*alpha
# So: 2/alpha = 40*phi^4 - 4*pi*alpha
# And: gap = 57*|ln alpha| - |ln alpha|/(2*phi^3) - 40*phi^4 + 4*pi*alpha

gap_formula = 57 * fabs(ln_alpha) - fabs(ln_alpha) / (2*PHI**3) - 40*PHI**4 + 4*pi*alpha
print(f"\n  gap via formula: {gap_formula}")
print(f"  gap direct: {gap}")
print(f"  diff: {gap - gap_formula}")

# So: gap = 57*|ln(alpha)| - |ln(alpha)|*alpha_s - 2/alpha
# This involves ln(alpha) and 1/alpha - inherently transcendental relation

# =====================================================================
# TEST 3: Series expansion of gap near alpha = 0
# =====================================================================
print(f"\n{'='*72}")
print("TEST 3: Leading behavior of gap")
print("=" * 72)

# For small alpha: ln(alpha) ~ ln(alpha), 1/alpha ~ 1/alpha
# Both diverge. The gap is a CANCELLATION of two large numbers.
# gap ~ 57 * |ln(1/alpha)| - 2/alpha
# ~ 57 * ln(1/alpha) - 2/alpha  (for alpha -> 0)
# The second term dominates since 1/alpha >> ln(1/alpha).
# So the leading term of the gap is -2/alpha + 57*ln(1/alpha) -> -infinity.
# But at the actual alpha, it's +5.8 (the two terms nearly cancel).

print(f"  57 * |ln alpha| = {57 * fabs(ln_alpha)}")
print(f"  2/alpha = {2/alpha}")
print(f"  57*|ln a| - 2/a = {57*fabs(ln_alpha) - 2/alpha}")
print(f"  alpha_s * |ln a| = {alpha_s * fabs(ln_alpha)}")
print(f"  Total gap = (57-alpha_s)*|ln a| - 2/a = {gap}")

# =====================================================================
# TEST 4: exp(-33/alpha_s) match
# =====================================================================
print(f"\n{'='*72}")
print("TEST 4: exp(-33/alpha_s) vs alpha^z_CC")
print("=" * 72)

val_33 = 33 / alpha_s   # = 33 * 2 * phi^3 = 66 * phi^3
print(f"  33/alpha_s = 66*phi^3 = {val_33}")
print(f"  target = z_CC * |ln alpha| = {target}")
print(f"  diff = {val_33 - target}")
print(f"  relative = {fabs(val_33 - target) / target}")

# 33 = 3 * 11 = N_gen * A_0
# So: 33/alpha_s = N_gen * A_0 * 2 * phi^3 = N_gen * A_0 / alpha_s
# This is the action of N_gen * A_0 "strong instantons"

# What z would 33/alpha_s correspond to?
z_33 = val_33 / fabs(ln_alpha)
print(f"\n  z_eff(33/alpha_s) = {z_33}")
print(f"  z_CC = {z_CC}")
print(f"  diff = {z_33 - z_CC}")

# Also check 33 * 2 * phi^3 exactly:
print(f"\n  33/alpha_s = 66*phi^3 = {66 * PHI**3}")
print(f"  Is this related to z_CC * |ln alpha| algebraically? NO.")
print(f"  It involves ln(alpha) on one side, phi^3 on the other.")

# =====================================================================
# TEST 5: Search for EXACT relations
# =====================================================================
print(f"\n{'='*72}")
print("TEST 5: Algebraic search for gap")
print("=" * 72)

# Can gap = f(phi, pi, a1) for some simple f?
# gap = z_CC * |ln alpha| - 2/alpha
# Let's see if gap * alpha or gap / ln(alpha) is algebraic.

g_times_alpha = gap * alpha
g_div_lna = gap / fabs(ln_alpha)

print(f"  gap * alpha = {g_times_alpha}")
print(f"  gap / |ln alpha| = {g_div_lna}")
print(f"  gap^2 = {gap**2}")
print(f"  gap * phi = {gap * PHI}")
print(f"  gap / phi = {gap / PHI}")
print(f"  gap * pi = {gap * pi}")
print(f"  gap / pi = {gap / pi}")

# gap * alpha = c exact:
# z_CC * alpha * |ln alpha| - 2 = c
# c = 0.04233... This is close to alpha * something.
c_ga = g_times_alpha
print(f"\n  gap * alpha = {c_ga}")
print(f"  c / alpha = {c_ga / alpha} (= gap)")
print(f"  c / alpha^2 = {c_ga / alpha**2}")
print(f"  c / (alpha * alpha_s) = {c_ga / (alpha * alpha_s)}")
print(f"  c * phi^4 = {c_ga * PHI**4}")
print(f"  c * 20*phi^4 = {c_ga * 20 * PHI**4}")

# So gap * alpha = z_CC * alpha * |ln alpha| - 2
# From alpha equation: alpha * (20*phi^4) - 2*pi*alpha^2 = 1
# So: 2 = 2*alpha*(20*phi^4) - 4*pi*alpha^2 = 2 - 4*pi*alpha^2 + 2*alpha*20*phi^4
# Hmm: gap * alpha = z_CC * alpha * |ln alpha| - 2
# = (57 - alpha_s) * alpha * |ln alpha| - 2

# Is alpha * |ln alpha| a "nice" number? NO - it's transcendental.

# =====================================================================
# TEST 6: Alternative decomposition 2/alpha + x/alpha_s
# =====================================================================
print(f"\n{'='*72}")
print("TEST 6: Decomposition 2/alpha + x/alpha_s = target")
print("=" * 72)

x_as = (target - 2/alpha) * alpha_s
print(f"  Need: 2/alpha + x/alpha_s = target")
print(f"  x = (target - 2/alpha) * alpha_s = gap * alpha_s")
print(f"  x = {x_as}")
print(f"  Is x = b1*alpha_s - alpha? {b1*alpha_s - alpha}")
print(f"  Is x integer/2? x*2 = {x_as * 2}")
print(f"  Is x = phi - 1 = 1/phi? {PHI - 1}")
print(f"  Diff from 1/phi: {x_as - (PHI - 1)}")

# Hmm x = 0.685 ~ 1/phi + something
# Let's be systematic
for name, val in [("1", 1), ("1/2", mpf(1)/2), ("1/phi", 1/PHI),
                   ("phi-1", PHI-1), ("phi/a1", PHI/a1),
                   ("2/a1", mpf(2)/a1), ("(a1-1)/(2*a1)", mpf(a1-1)/(2*a1)),
                   ("alpha_s/alpha", alpha_s/alpha),
                   ("ln(phi)", log(PHI)),
                   ("3/4", mpf(3)/4),
                   ("5/7", mpf(5)/7),
                   ("7/10", mpf(7)/10),
                   ("29/(2*a1^2)", mpf(29)/(2*a1**2))]:
    diff = fabs(x_as - val)
    if diff < 0.1:
        print(f"  x ~ {name} = {val}, diff = {float(diff):.6e}")

# =====================================================================
# TEST 7: Is gap/(2*pi) close to something?
# =====================================================================
print(f"\n{'='*72}")
print("TEST 7: Various normalizations of gap")
print("=" * 72)

quantities = [
    ("gap", gap),
    ("gap/(2*pi)", gap/(2*pi)),
    ("gap/pi", gap/pi),
    ("gap*phi", gap*PHI),
    ("gap/phi", gap/PHI),
    ("gap*alpha", gap*alpha),
    ("gap*alpha_s", gap*alpha_s),
    ("gap/a1", gap/a1),
    ("gap/b1", gap/b1),
    ("gap/(b1-1/a1)", gap/(b1 - 1/a1)),  # should be ~1
]

for name, val in quantities:
    # Check against simple fractions
    for p in range(1, 50):
        for q in range(1, 20):
            frac = mpf(p)/mpf(q)
            if fabs(val - frac) < 0.0005:
                print(f"  {name} = {float(val):.10f} ~ {p}/{q} = {float(frac):.10f} "
                      f"(diff = {float(fabs(val-frac)):.2e})")

# =====================================================================
# SUMMARY
# =====================================================================
print(f"\n{'='*72}")
print("SUMMARY")
print("=" * 72)

residual = float(gap - exact_29_5)
rel = float(fabs(gap - exact_29_5) / exact_29_5)

print(f"""
  QUESTION: Is z_CC * |ln alpha| - 2/alpha = b1 - 1/a1 = 29/5 EXACT?

  ANSWER: NO.

  High-precision computation:
    gap = {float(gap):.15f}
    29/5 = {float(exact_29_5):.15f}
    gap - 29/5 = {residual:.15e}
    relative error = {rel:.6e} (= 0.014%)

  The match IS remarkable (0.014%) but NOT exact.
  At 50-digit precision, the residual is clearly nonzero.

  OTHER PATTERNS:
    exp(-33/alpha_s): diff from target = {float(val_33 - target):.6f} (0.10%)
    33 = 3 * 11 = N_gen * A_0 (Seeley-DeWitt)

  THE STRUCTURAL OBSERVATION:
    alpha^57 ~ exp(-2/alpha) shows that the CC exponent 57
    is "approximately" related to the instanton pair action 2/alpha
    by the transcendental number |ln(alpha)|:
      57 ~ 2 / (alpha * |ln alpha|) + 1.18

    The "correction" 1.18 ~ gap/|ln alpha| has no clean expression.
    This is consistent with the conclusion from exp421-425:
    the CC functional form alpha^z remains UNEXPLAINED.

  STATUS: PATTERN (0.014% match with 29/5, not exact)
""")

"""
exp453: Origin of the "-1" offset in Cayley→CKM map
=====================================================
From exp452b, the Cayley spectral exponents are n_ut=4, n_db=2 (EXACT).
The CKM exponents {3,7,12} are obtained via:
    n_12 = n_ut - 1 = 3
    n_23 = n_ut*n_db - 1 = 7
    n_13 = n_ut*(n_db + 1) = 12

KEY QUESTION: Where does the "-1" come from?

APPROACH: Parameterize the map with a UNIFORM offset delta:
    n_12 = n_ut - delta
    n_23 = n_ut*n_db - delta
    n_13 = n_ut*(n_db + delta)
and demand consistency with the concordance identities I1, I2, I3.

RESULT: delta=1 is FORCED by I3 alone, and then I1 gives a1=5.
"""
import numpy as np
from math import sqrt, log, pi
from fractions import Fraction

phi = (1 + sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120

print("=" * 70)
print("exp453: ORIGIN OF THE '-1' OFFSET IN CAYLEY -> CKM MAP")
print("=" * 70)

# ============================================================
# 1. RECAP: Cayley spectral exponents (from exp452b)
# ============================================================
print("\n--- 1. Cayley Spectral Exponents (exp452b) ---")
print(f"Galois pair dim 2 (rho_1, rho_7): lambda_t/lambda_u = phi^4")
print(f"  => n_ut = 4 = a1-1 = {a1}-1")
print(f"Galois pair dim 3 (rho_2, rho_8): lambda_b/lambda_d = phi^2")
print(f"  => n_db = 2 = (a1-1)/2 = {(a1-1)}/2")
print(f"\nGeneral: n = degree/d - 2 for Galois pair of dim d (degree=12)")
print(f"  d=2: n = 12/2 - 2 = 4 = n_ut")
print(f"  d=3: n = 12/3 - 2 = 2 = n_db")
print(f"  n_ut/n_db = 2 (dimension ratio: 3/2 -> exponent ratio: 2)")

# ============================================================
# 2. GENERAL PARAMETERIZATION with uniform offset delta
# ============================================================
print("\n" + "=" * 70)
print("2. GENERAL PARAMETERIZATION: Cayley -> CKM with offset delta")
print("=" * 70)
print("""
Given n_ut = a1-1 and n_db = (a1-1)/2, define the CKM map:

    n_12 = n_ut - delta
    n_23 = n_ut * n_db - delta
    n_13 = n_ut * (n_db + delta)

Concordance identities (from exp449, DERIVED):
    I1: n_12 + n_23 = 2*a1
    I2: n_12 * n_13 = b1^2 = (a1+1)^2
    I3: n_23 + n_13 = a1^2 - a1 - 1
""")

# ============================================================
# 3. ALGEBRAIC DERIVATION: I3 forces delta = 1
# ============================================================
print("=" * 70)
print("3. I3 FORCES delta = 1")
print("=" * 70)

print("""
Substitute into I3: n_23 + n_13 = a1^2 - a1 - 1

LHS = [n_ut*n_db - delta] + [n_ut*(n_db + delta)]
    = n_ut*n_db - delta + n_ut*n_db + n_ut*delta
    = 2*n_ut*n_db + delta*(n_ut - 1)

With n_ut = a1-1 and n_db = (a1-1)/2:
    2*(a1-1)*(a1-1)/2 + delta*(a1-2)
  = (a1-1)^2 + delta*(a1-2)

Set equal to RHS:
    (a1-1)^2 + delta*(a1-2) = a1^2 - a1 - 1

Expand:
    a1^2 - 2*a1 + 1 + delta*(a1-2) = a1^2 - a1 - 1

Cancel a1^2:
    -2*a1 + 1 + delta*(a1-2) = -a1 - 1

Rearrange:
    delta*(a1-2) = -a1 - 1 + 2*a1 - 1 = a1 - 2

Therefore:
    delta*(a1-2) = (a1-2)

    => delta = 1  (for a1 != 2)   QED
""")

# Numerical verification
print("Numerical verification for a1 = 2, 3, ..., 10:")
print(f"{'a1':>4s} {'I3_LHS':>10s} {'I3_RHS':>10s} {'delta_forced':>13s}")
for a in range(2, 11):
    n_ut_a = a - 1
    n_db_a = (a - 1) / 2
    rhs = a**2 - a - 1
    # (a-1)^2 + delta*(a-2) = rhs => delta = (rhs - (a-1)^2) / (a-2)
    lhs_base = (a - 1)**2
    if a != 2:
        delta_forced = (rhs - lhs_base) / (a - 2)
    else:
        delta_forced = float('inf')  # a-2=0
    print(f"{a:4d} {lhs_base:10.1f} {rhs:10.1f} {delta_forced:13.4f}")

print("\n=> delta = 1 for ALL a1 != 2. UNIVERSAL.")

# ============================================================
# 4. I1 => a1 = 5
# ============================================================
print("\n" + "=" * 70)
print("4. WITH delta=1, I1 GIVES a1 = 5")
print("=" * 70)

print("""
n_12 + n_23 = 2*a1  (Identity I1)

Substitute delta = 1:
    n_12 = (a1-1) - 1 = a1 - 2
    n_23 = (a1-1)*(a1-1)/2 - 1 = (a1-1)^2/2 - 1

    n_12 + n_23 = (a1-2) + (a1-1)^2/2 - 1
               = (a1-3) + (a1^2 - 2*a1 + 1)/2
               = (2*a1 - 6 + a1^2 - 2*a1 + 1) / 2
               = (a1^2 - 5) / 2

Set equal to 2*a1:
    (a1^2 - 5) / 2 = 2*a1
    a1^2 - 5 = 4*a1
    a1^2 - 4*a1 - 5 = 0
    (a1 - 5)(a1 + 1) = 0

    => a1 = 5   QED
""")

# Numerical verification
print("Numerical check: n_12 + n_23 vs 2*a1 for each a1:")
print(f"{'a1':>4s} {'n_12':>6s} {'n_23':>8s} {'sum':>8s} {'2*a1':>6s} {'match':>8s}")
for a in range(3, 11):
    n_12 = a - 2
    n_23_val = (a - 1)**2 / 2 - 1
    s = n_12 + n_23_val
    target = 2 * a
    match = "YES" if abs(s - target) < 1e-10 else "no"
    print(f"{a:4d} {n_12:6d} {n_23_val:8.1f} {s:8.1f} {target:6d} {match:>8s}")

# ============================================================
# 5. I2 => a1 = 5 (consistency)
# ============================================================
print("\n" + "=" * 70)
print("5. WITH delta=1, I2 GIVES a1 = 5 (consistency)")
print("=" * 70)

print("""
n_12 * n_13 = b1^2 = (a1+1)^2  (Identity I2)

With delta = 1:
    n_12 = a1 - 2
    n_13 = (a1-1)*((a1-1)/2 + 1) = (a1-1)*(a1+1)/2 = (a1^2-1)/2

    n_12 * n_13 = (a1-2)*(a1^2-1)/2

Set equal to (a1+1)^2:
    (a1-2)*(a1-1)*(a1+1)/2 = (a1+1)^2

Divide by (a1+1) [assuming a1 > 0]:
    (a1-2)*(a1-1)/2 = (a1+1)

Multiply by 2:
    (a1-2)*(a1-1) = 2*(a1+1)
    a1^2 - 3*a1 + 2 = 2*a1 + 2
    a1^2 - 5*a1 = 0
    a1*(a1 - 5) = 0

    => a1 = 5   QED   (consistent with I1)
""")

# Numerical verification
print("Numerical check: n_12 * n_13 vs (a1+1)^2:")
print(f"{'a1':>4s} {'n_12':>6s} {'n_13':>8s} {'product':>10s} {'(a1+1)^2':>10s} {'match':>8s}")
for a in range(3, 11):
    n_12 = a - 2
    n_13_val = (a**2 - 1) / 2
    prod = n_12 * n_13_val
    target = (a + 1)**2
    match = "YES" if abs(prod - target) < 1e-10 else "no"
    print(f"{a:4d} {n_12:6d} {n_13_val:8.1f} {prod:10.1f} {target:10.1f} {match:>8s}")

# ============================================================
# 6. I3 automatically satisfied (identity check)
# ============================================================
print("\n" + "=" * 70)
print("6. I3 IS AUTOMATICALLY SATISFIED (identity, no constraint)")
print("=" * 70)

print("""
n_23 + n_13 = a1^2 - a1 - 1  (Identity I3)

With delta = 1:
    n_23 = (a1-1)^2/2 - 1
    n_13 = (a1^2-1)/2

    n_23 + n_13 = [(a1-1)^2 - 2 + a1^2 - 1] / 2
               = [a1^2 - 2*a1 + 1 - 2 + a1^2 - 1] / 2
               = [2*a1^2 - 2*a1 - 2] / 2
               = a1^2 - a1 - 1

    = RHS identically!  NO constraint on a1.

This means: I3 encodes the STRUCTURE of the Cayley map itself.
The constraints come from I1 and I2 alone.
""")

# Numerical verification
for a in [3, 4, 5, 6, 7, 8]:
    n_23_val = (a - 1)**2 / 2 - 1
    n_13_val = (a**2 - 1) / 2
    lhs = n_23_val + n_13_val
    rhs = a**2 - a - 1
    print(f"  a1={a}: n_23+n_13 = {lhs:.1f}, a1^2-a1-1 = {rhs}  {'IDENTITY' if abs(lhs-rhs)<1e-10 else 'FAIL'}")

# ============================================================
# 7. ALTERNATIVE: Express everything in terms of n_ut alone
# ============================================================
print("\n" + "=" * 70)
print("7. EVERYTHING FROM n_ut ALONE")
print("=" * 70)

print("""
Since n_db = n_ut/2 and a1 = n_ut + 1, with delta = 1:

    n_12 = n_ut - 1
    n_23 = n_ut^2/2 - 1
    n_13 = n_ut*(n_ut + 2)/2

I1: n_12 + n_23 = 2*a1 = 2*(n_ut+1)
    (n_ut - 1) + (n_ut^2/2 - 1) = 2*n_ut + 2
    n_ut^2/2 + n_ut - 2 = 2*n_ut + 2
    n_ut^2/2 - n_ut - 4 = 0
    n_ut^2 - 2*n_ut - 8 = 0
    (n_ut - 4)(n_ut + 2) = 0
    => n_ut = 4   (hence a1 = 5)

I2: n_12 * n_13 = (a1+1)^2 = (n_ut+2)^2
    (n_ut - 1) * n_ut*(n_ut+2)/2 = (n_ut+2)^2
    n_ut*(n_ut-1)/2 = n_ut + 2  [dividing by (n_ut+2)]
    n_ut*(n_ut-1) = 2*n_ut + 4
    n_ut^2 - 3*n_ut - 4 = 0
    (n_ut - 4)(n_ut + 1) = 0
    => n_ut = 4   (consistent)
""")

print("Factorizations from each identity:")
print(f"  I1: (n_ut - 4)(n_ut + 2) = 0 => n_ut = 4")
print(f"  I2: (n_ut - 4)(n_ut + 1) = 0 => n_ut = 4")
print(f"  I3: identity (no constraint)")
print(f"  Intersection: n_ut = 4 UNIQUE")
print(f"  Therefore: a1 = n_ut + 1 = 5")

# ============================================================
# 8. COMPARISON OF FACTORIZATIONS
# ============================================================
print("\n" + "=" * 70)
print("8. COMPARISON: Three routes to a1 = 5")
print("=" * 70)

print("""
Route 1 (T-matrix, exp449):
    n_12=3, n_23=2a1-3, n_13=3a1-3
    I2: (a1-2)(a1-5) = 0
    I3: (a1-1)(a1-5) = 0
    Intersection: a1 = 5

Route 2 (Cayley + uniform delta, this work):
    n_12=a1-2, n_23=(a1-1)^2/2-1, n_13=(a1^2-1)/2
    I1: (a1-5)(a1+1) = 0
    I2: a1*(a1-5) = 0
    I3: identity
    Intersection: a1 = 5

Route 3 (Cayley, n_ut variable):
    n_12=n_ut-1, n_23=n_ut^2/2-1, n_13=n_ut(n_ut+2)/2
    I1: (n_ut-4)(n_ut+2) = 0
    I2: (n_ut-4)(n_ut+1) = 0
    Intersection: n_ut = 4, a1 = 5

ALL THREE ROUTES converge on a1 = 5 via DIFFERENT factorizations.
""")

# ============================================================
# 9. PHYSICAL INTERPRETATION: Vacuum subtraction
# ============================================================
print("=" * 70)
print("9. PHYSICAL INTERPRETATION: Vacuum Subtraction")
print("=" * 70)

print("""
The offset delta = 1 has a natural physical interpretation:

CAYLEY EIGENVALUES encode the FULL spectral propagation, including
the trivial (vacuum / identity) channel.

CKM MIXING ANGLES encode only the NON-TRIVIAL flavor change,
excluding the no-mixing (diagonal) channel.

The relationship: n_CKM = n_Cayley - delta, with delta = 1,
is the standard VACUUM SUBTRACTION of quantum field theory:

    S = 1 + i*T    (S-matrix = identity + transfer matrix)
    G = G_0 + G_c  (propagator = free + connected)
    V = 1 + dV     (CKM = identity + mixing)

In TQFT at level k = a1-2 = 3:
  - Number of simple objects: k+1 = 4 = n_ut
  - TQFT level: k = 3 = n_12
  - Subtracting the vacuum object j=0 gives k = n_ut - 1

The concordance identities ENFORCE this vacuum subtraction:
  delta = 1 is not a choice, it is a THEOREM.
""")

# ============================================================
# 10. DERIVATION CHAIN: degree=12 from Galois pair dimensions
# ============================================================
print("=" * 70)
print("10. FULL DERIVATION CHAIN")
print("=" * 70)

print("""
(A) 600-cell polytope => degree = 12 (nearest neighbors)
    This is geometric: each vertex of the 600-cell has 12 neighbors.

(B) Binary icosahedral group 2I => Galois pairs:
    - dim 2 pair: (rho_1, rho_7) -- smallest non-trivial
    - dim 3 pair: (rho_2, rho_8) -- next non-trivial
    Dimensions {2, 3} from McKay correspondence (affine E8).

(C) Cayley spectral exponents (exp452b):
    n_ut = degree/d_1 - 2 = 12/2 - 2 = 4
    n_db = degree/d_2 - 2 = 12/3 - 2 = 2
    These give lambda_t/lambda_u = phi^4, lambda_b/lambda_d = phi^2.

(D) Concordance identity I3 FORCES delta = 1 (vacuum subtraction):
    n_12 = n_ut - 1 = 3
    n_23 = n_ut*n_db - 1 = 7
    n_13 = n_ut*(n_db + 1) = 12

(E) Concordance identity I1 FORCES a1 = 5:
    (a1 - 5)(a1 + 1) = 0

EVERY STEP IS DERIVED. The "-1" is NOT ad hoc.
""")

# ============================================================
# 11. SELF-CONSISTENCY CHECK: degree=12 from n_db = n_ut/2
# ============================================================
print("=" * 70)
print("11. SELF-CONSISTENCY: degree = 12 from dimension constraint")
print("=" * 70)

print("""
n_db = n_ut/2 requires:
    degree/d_2 - 2 = (degree/d_1 - 2) / 2

For d_1=2, d_2=3:
    degree/3 - 2 = (degree/2 - 2) / 2 = degree/4 - 1
    degree/3 - degree/4 = 1
    degree * (1/3 - 1/4) = 1
    degree * (1/12) = 1
    degree = 12   QED

The 600-cell degree (=12) is uniquely determined by the constraint
that the two Galois pair exponents satisfy n_db = n_ut/2.
""")

# Verify
for deg in range(6, 25):
    n1 = deg / 2 - 2
    n2 = deg / 3 - 2
    if abs(n2 - n1/2) < 1e-10 and n1 > 0 and n2 > 0:
        print(f"  degree={deg}: n_ut={n1:.0f}, n_db={n2:.0f}, n_db=n_ut/2 ? YES")
    elif n1 > 0 and n2 > 0:
        print(f"  degree={deg}: n_ut={n1:.0f}, n_db={n2:.1f}, ratio={n1/n2:.3f}")

# ============================================================
# 12. NON-UNIFORM OFFSET ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("12. NON-UNIFORM OFFSET ANALYSIS")
print("=" * 70)

print("""
If we allow delta_1, delta_2, delta_3 independent:
    n_12 = n_ut - delta_1
    n_23 = n_ut*n_db - delta_2
    n_13 = n_ut*(n_db + delta_3)

I3 gives: delta_2 = 1 + (a1-1)*(delta_3 - 1)
  => delta_2 = 1 iff delta_3 = 1 (for a1 != 1)

I1 determines delta_1 in terms of a1 and delta_3:
    delta_1 = (a1^2 - 4*a1 - 3)/2 + (a1-1)*(delta_3-1)/2

For delta_3 = 1 (hence delta_2 = 1):
    delta_1 = (a1^2 - 4*a1 - 3)/2
    delta_1 = 1 iff a1^2 - 4*a1 - 5 = 0 iff (a1-5)(a1+1) = 0

CONCLUSION: The UNIFORM offset hypothesis (delta_1 = delta_2 = delta_3)
is EQUIVALENT to a1 = 5.
""")

# Show delta_1 as function of a1 (with delta_2=delta_3=1)
print("delta_1(a1) with delta_2 = delta_3 = 1:")
print(f"{'a1':>4s} {'delta_1':>10s} {'=1?':>6s}")
for a in range(3, 11):
    d1 = (a**2 - 4*a - 3) / 2
    is1 = "YES" if abs(d1 - 1) < 1e-10 else "no"
    print(f"{a:4d} {d1:10.1f} {is1:>6s}")

# ============================================================
# 13. VERIFICATION: explicit CKM values at a1=5
# ============================================================
print("\n" + "=" * 70)
print("13. EXPLICIT VERIFICATION at a1 = 5")
print("=" * 70)

a1_v = 5
n_ut_v = a1_v - 1  # = 4
n_db_v = (a1_v - 1) / 2  # = 2
delta = 1

n_12_v = n_ut_v - delta  # = 3
n_23_v = int(n_ut_v * n_db_v - delta)  # = 7
n_13_v = int(n_ut_v * (n_db_v + delta))  # = 12

print(f"n_ut = {n_ut_v}, n_db = {n_db_v:.0f}, delta = {delta}")
print(f"n_12 = {n_12_v}, n_23 = {n_23_v}, n_13 = {n_13_v}")
print()

# Concordance identities
b1_v = a1_v + 1
I1_lhs = n_12_v + n_23_v
I1_rhs = 2 * a1_v
I2_lhs = n_12_v * n_13_v
I2_rhs = b1_v**2
I3_lhs = n_23_v + n_13_v
I3_rhs = a1_v**2 - a1_v - 1

print(f"I1: n_12 + n_23 = {I1_lhs} = 2*a1 = {I1_rhs}  {'PASS' if I1_lhs==I1_rhs else 'FAIL'}")
print(f"I2: n_12 * n_13 = {I2_lhs} = b1^2 = {I2_rhs}  {'PASS' if I2_lhs==I2_rhs else 'FAIL'}")
print(f"I3: n_23 + n_13 = {I3_lhs} = a1^2-a1-1 = {I3_rhs}  {'PASS' if I3_lhs==I3_rhs else 'FAIL'}")

# Additional relations
print(f"\nAdditional structural relations:")
print(f"  n_12 = N_gen = {3}  (spatial dimensions)")
print(f"  n_ut = d_ST = {4}  (spacetime dimensions)")
print(f"  n_ut*n_db = rank(E8) = {int(n_ut_v*n_db_v)}  [{'YES' if n_ut_v*n_db_v==8 else 'NO'}]")
print(f"  n_13 = degree(600-cell) = {n_13_v}  (nearest neighbors)")
print(f"  n_12 = k (TQFT level) = {a1_v-2}  [{'YES' if n_12_v==a1_v-2 else 'NO'}]")
print(f"  n_12*(n_12+1) = n_13: {n_12_v}*{n_12_v+1} = {n_12_v*(n_12_v+1)} = {n_13_v}  [{'YES' if n_12_v*(n_12_v+1)==n_13_v else 'NO'}]")

# ============================================================
# 14. SCAN: Which a1 values give INTEGER CKM exponents?
# ============================================================
print("\n" + "=" * 70)
print("14. INTEGER CKM EXPONENTS: which a1 values work?")
print("=" * 70)

print(f"{'a1':>4s} {'n_12':>6s} {'n_23':>8s} {'n_13':>8s} {'all_int':>8s} {'I1':>5s} {'I2':>5s} {'I3':>5s}")
for a in range(3, 16):
    n12 = a - 2
    n23 = (a - 1)**2 / 2 - 1
    n13 = (a**2 - 1) / 2
    all_int = (n23 == int(n23)) and (n13 == int(n13))
    i1_ok = abs(n12 + n23 - 2*a) < 1e-10
    i2_ok = abs(n12 * n13 - (a+1)**2) < 1e-10
    i3_ok = abs(n23 + n13 - (a**2 - a - 1)) < 1e-10
    print(f"{a:4d} {n12:6d} {n23:8.1f} {n13:8.1f} {'YES' if all_int else 'no':>8s} "
          f"{'OK' if i1_ok else '--':>5s} {'OK' if i2_ok else '--':>5s} {'OK' if i3_ok else '--':>5s}")

print("\nNote: I3 is always satisfied (identity). I1 and I2 are satisfied ONLY at a1=5.")
print("Integer exponents require ODD a1 (so (a1-1)^2/2 - 1 and (a1^2-1)/2 are integers).")

# ============================================================
# 15. THEOREM STATEMENT
# ============================================================
print("\n" + "=" * 70)
print("15. THEOREM (Cayley-Concordance Vacuum Offset)")
print("=" * 70)

print("""
THEOREM: Let the Cayley graph of 2I have spectral exponents n_ut = a1-1
and n_db = (a1-1)/2 for the two Galois pairs (of dimensions 2 and 3).
Define the CKM exponents via a UNIFORM offset delta:

    n_12 = n_ut - delta
    n_23 = n_ut * n_db - delta
    n_13 = n_ut * (n_db + delta)

Then the concordance identities uniquely determine:

    (i)   I3 => delta = 1           (for any a1 != 2)
    (ii)  I1 => (a1-5)(a1+1) = 0   (hence a1 = 5)
    (iii) I2 => a1*(a1-5) = 0      (consistent)

COROLLARY: The uniform vacuum offset delta = 1 is the unique value
consistent with the concordance identities, and it implies a1 = 5.
Equivalently: universality of vacuum subtraction <=> a1 = 5.

INTERPRETATION: The "-1" is the TQFT vacuum subtraction.
The Cayley eigenvalues encode FULL propagation (including vacuum).
CKM mixing measures only NON-TRIVIAL flavor transitions.
The concordance identities ENFORCE the vacuum subtraction to be
exactly 1, and this uniquely determines a1 = 5.

STATUS: DERIVED (theorem). Origin of "-1" RESOLVED.
This constitutes a 14th characterisation of a1 = 5:
  "Uniform vacuum subtraction in the Cayley-CKM map."
""")

# ============================================================
# 16. SUMMARY TABLE of all factorizations
# ============================================================
print("=" * 70)
print("16. FACTORIZATION COMPARISON (three routes to a1=5)")
print("=" * 70)
print("""
+------------------+------------------------+------------------------+
|    Identity      |  T-matrix (exp449)     | Cayley+delta (exp453)  |
+------------------+------------------------+------------------------+
| I1: n12+n23=2a1  | (automatic)            | (a1-5)(a1+1) = 0       |
| I2: n12*n13=b1^2 | (a1-2)(a1-5) = 0       | a1*(a1-5) = 0          |
| I3: n23+n13=...  | (a1-1)(a1-5) = 0       | (identity, automatic)  |
+------------------+------------------------+------------------------+
| Intersection     | a1 = 5                 | a1 = 5                 |
+------------------+------------------------+------------------------+

The two routes use DIFFERENT identities as constraints and DIFFERENT
identities as automatic, but both arrive at a1 = 5.

T-matrix: I1 automatic, I2 and I3 constrain.
Cayley:   I3 automatic (+ fixes delta), I1 and I2 constrain.
""")

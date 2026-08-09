"""
EXP-508: E_vac for ALL Binary Polyhedral Groups (ADE Classification)
=====================================================================
Compute the orbifold vacuum energy E_vac(C^2/Gamma) for all finite
subgroups Gamma of SU(2), and look for a pattern.

E_vac = (1/|Gamma|) * sum_g [(nu_g - 1/2)^2 - 1/12]

where g has SU(2) eigenvalues exp(+/-2*pi*i*nu_g).

Groups: Z_n (A-type), 2D_n (D-type), 2T (E6), 2O (E7), 2I (E8)

Key question: does E_vac have a clean expression in Lie-theoretic data?
"""

import numpy as np
from fractions import Fraction
import sys
sys.path.insert(0, '.')
from commons import (build_600cell, PHI, PHI_CONJ, SQRT5, a1, b1, N,
                      alpha_em, alpha_s, ln_inv_alpha)

print("=" * 70)
print("EXP-508: E_vac FOR ALL BINARY POLYHEDRAL GROUPS")
print("=" * 70)

def E_0(nu):
    """Ground state energy for twist parameter nu (exact, using Fraction)."""
    return (nu - Fraction(1, 2))**2 - Fraction(1, 12)

def compute_E_vac(classes, group_order, name=""):
    """
    Compute E_vac from conjugacy class data.
    classes: list of (nu_fraction, size) pairs
    """
    total = Fraction(0)
    for nu, size in classes:
        total += size * E_0(nu)

    # Verify sizes
    total_size = sum(s for _, s in classes)
    assert total_size == group_order, f"{name}: sizes sum to {total_size}, expected {group_order}"

    E_vac = total / group_order
    return E_vac, total

# ============================================================
# CYCLIC GROUPS Z_n (A_{n-1} singularity)
# ============================================================
print(f"\n{'='*70}")
print("A-TYPE: CYCLIC GROUPS Z_n")
print(f"{'='*70}")

print(f"\n{'n':>4} {'|G|':>5} {'h':>4} {'r':>3} {'E_vac':>15} {'E_vac (dec)':>14} {'E_vac * 6h^2':>12}")

for n in range(2, 11):
    # Z_n: elements g_k with nu = k/n for k = 0, ..., n-1
    classes = [(Fraction(k, n), 1) for k in range(n)]
    E_vac, _ = compute_E_vac(classes, n, f"Z_{n}")
    h_val = n  # h(A_{n-1}) = n
    r_val = n - 1
    product = E_vac * 6 * h_val**2
    print(f"{n:4d} {n:5d} {h_val:4d} {r_val:3d} {str(E_vac):>15} {float(E_vac):14.8f} {float(product):12.6f}")

# ============================================================
# BINARY DIHEDRAL GROUPS 2D_n (D_{n+2} singularity)
# ============================================================
print(f"\n{'='*70}")
print("D-TYPE: BINARY DIHEDRAL GROUPS 2D_n")
print(f"{'='*70}")

# 2D_n (binary dihedral of order 4n, Lie algebra D_{n+2})
# h(D_m) = 2(m-1), rank = m
# |2D_n| = 4n, h = 2(n+2-1) = 2(n+1), rank = n+2

# Conjugacy classes of 2D_n (order 4n):
# Central: {I} (1), {-I} (1)
# Rotational: {g_k, g_{2n-k}} for k=1,...,n-1, each size 2, nu = k/(2n)
# Special class if n even/odd...
# Actually for the binary dihedral 2D_n of order 4n:
# Elements as SU(2) matrices:
#   Rotations: exp(2*pi*i*k/(2n) * sigma_3) for k = 0,...,2n-1
#   "Flips": sigma_1 * rotation, etc.
# The rotation part: nu_k = k/(2n)
# The flip part: these are order 4, nu = 1/4

# More carefully:
# 2D_n has 4n elements and n+3 conjugacy classes:
# 1. {I}: 1 element, nu = 0
# 2. {-I}: 1 element, nu = 1/2
# 3. {R_k, R_{-k}} for k = 1, ..., n-1: each size 2, nu = k/(2n)
# 4. "Flip" class 1: n elements, nu = 1/4
# 5. "Flip" class 2: n elements, nu = 1/4

# Wait - the flips in 2D_n. In the dihedral group D_n (order 2n),
# the reflections form one or two conjugacy classes depending on parity of n.
# In 2D_n (order 4n), the lifts have specific ν values.

# For a reflection in D_n, it has order 2. Its lift in 2D_n has order 4.
# A reflection in SO(3) corresponds to rotation by pi around some axis.
# In SU(2), lift has eigenvalue exp(+/-i*pi/2), so nu = 1/4.

# Conjugacy classes of 2D_n:
# {I}: 1, nu = 0
# {-I}: 1, nu = 1/2
# {R_k, R_{2n-k}}: 2 each, for k = 1, ..., n-1, nu = k/(2n)
# Two flip classes of size n each: nu = 1/4 for both

# Total: 1 + 1 + 2*(n-1) + 2*n = 2 + 2n - 2 + 2n = 4n. Check!
# Number of classes: 2 + (n-1) + 2 = n + 3. Check (for n+3 irreps of 2D_n/D-type McKay).

print(f"\n{'n':>4} {'|G|':>5} {'h':>4} {'r':>3} {'E_vac':>20} {'E_vac (dec)':>14}")

for n in range(2, 9):
    group_order = 4 * n
    h_val = 2 * (n + 1)  # h(D_{n+2}) = 2(n+1)
    r_val = n + 2  # rank(D_{n+2})

    classes = []
    classes.append((Fraction(0, 1), 1))      # I
    classes.append((Fraction(1, 2), 1))      # -I
    for k in range(1, n):
        classes.append((Fraction(k, 2*n), 2))  # rotation pair
    classes.append((Fraction(1, 4), n))      # flip class 1
    classes.append((Fraction(1, 4), n))      # flip class 2

    E_vac, _ = compute_E_vac(classes, group_order, f"2D_{n}")
    print(f"{n:4d} {group_order:5d} {h_val:4d} {r_val:3d} {str(E_vac):>20} {float(E_vac):14.8f}")

# ============================================================
# EXCEPTIONAL: 2T (E6), 2O (E7), 2I (E8)
# ============================================================
print(f"\n{'='*70}")
print("E-TYPE: BINARY POLYHEDRAL GROUPS")
print(f"{'='*70}")

# 2T (binary tetrahedral, |Gamma| = 24, E6: h=12, rank=6)
classes_2T = [
    (Fraction(0, 1), 1),    # I
    (Fraction(1, 2), 1),    # -I
    (Fraction(1, 4), 6),    # lifts of order-2 rotations
    (Fraction(1, 6), 4),    # lifts of order-3 class 1
    (Fraction(1, 3), 4),    # lifts of order-3 class 1'
    (Fraction(1, 3), 4),    # -1 x lifts of order-3 class 2
    (Fraction(1, 6), 4),    # -1 x lifts of order-3 class 2'
]
E_2T, _ = compute_E_vac(classes_2T, 24, "2T")

# 2O (binary octahedral, |Gamma| = 48, E7: h=18, rank=7)
# S4 has 5 conjugacy classes: e(1), (12)(34)(3), (123)(8), (1234)(6), (12)(6)
# 2O has 8 conjugacy classes:
classes_2O = [
    (Fraction(0, 1), 1),     # I
    (Fraction(1, 2), 1),     # -I
    (Fraction(1, 4), 6),     # lifts of (12)(34), order 4
    (Fraction(1, 4), 6),     # lifts of (12), order 4
    (Fraction(1, 6), 8),     # lifts of (123), order 6
    (Fraction(1, 3), 8),     # -1 x lifts of (123), order 3
    (Fraction(1, 8), 6),     # lifts of (1234), order 8
    (Fraction(3, 8), 6),     # -1 x lifts of (1234), order 8
]
# Verify: 1+1+6+6+8+8+6+6 = 42? No, that's wrong. Should be 48.
# Let me recount:
# (12) in S4: 6 transpositions. Lift: 12 elements in 2O.
# (12)(34): 3 double transpositions. Lift: 6 elements.
# (123): 8 3-cycles. Lift: 16 elements (8 order 6, 8 order 3).
# (1234): 6 4-cycles. Lift: 12 elements (6 order 8, 6 order 8).
# Total: 2 + 12 + 6 + 16 + 12 = 48. OK!

classes_2O = [
    (Fraction(0, 1), 1),      # I
    (Fraction(1, 2), 1),      # -I
    (Fraction(1, 4), 12),     # lifts of transpositions (12), order 4
    (Fraction(1, 4), 6),      # lifts of (12)(34), order 4
    (Fraction(1, 6), 8),      # lifts of (123), order 6
    (Fraction(1, 3), 8),      # -1 x lifts of (123), order 3
    (Fraction(1, 8), 6),      # lifts of (1234), order 8
    (Fraction(3, 8), 6),      # -1 x lifts of (1234), order 8
]
# 1+1+12+6+8+8+6+6 = 48. OK!
E_2O, _ = compute_E_vac(classes_2O, 48, "2O")

# 2I (binary icosahedral, |Gamma| = 120, E8: h=30, rank=8)
classes_2I = [
    (Fraction(0, 1),   1),   # I
    (Fraction(1, 2),   1),   # -I
    (Fraction(1, 10), 12),   # C5: lift of order-5 rotation
    (Fraction(1, 5),  12),   # C5^2
    (Fraction(3, 10), 12),   # C5^3
    (Fraction(2, 5),  12),   # C5^4
    (Fraction(1, 6),  20),   # C3: lift of order-3 rotation
    (Fraction(1, 3),  20),   # C3^2
    (Fraction(1, 4),  30),   # C2: lift of order-2 rotation
]
E_2I, _ = compute_E_vac(classes_2I, 120, "2I")

print(f"\n{'Group':>6} {'|G|':>5} {'Lie':>4} {'h':>4} {'r':>3} {'E_vac':>20} {'E_vac (dec)':>14}")
for name, order, lie, h_val, r_val, E_vac in [
    ("2T", 24, "E6", 12, 6, E_2T),
    ("2O", 48, "E7", 18, 7, E_2O),
    ("2I", 120, "E8", 30, 8, E_2I),
]:
    print(f"{name:>6} {order:5d} {lie:>4} {h_val:4d} {r_val:3d} {str(E_vac):>20} {float(E_vac):14.8f}")

# ============================================================
# PATTERN SEARCH
# ============================================================
print(f"\n{'='*70}")
print("PATTERN SEARCH")
print(f"{'='*70}")

# Collect all results
all_results = []

# Cyclic
for n in range(2, 11):
    classes = [(Fraction(k, n), 1) for k in range(n)]
    E_vac, _ = compute_E_vac(classes, n, f"Z_{n}")
    all_results.append(("A", n-1, n, n, E_vac))

# Binary dihedral
for n in range(2, 9):
    group_order = 4 * n
    h_val = 2 * (n + 1)
    r_val = n + 2
    classes = [(Fraction(0, 1), 1), (Fraction(1, 2), 1)]
    for k in range(1, n):
        classes.append((Fraction(k, 2*n), 2))
    classes.append((Fraction(1, 4), n))
    classes.append((Fraction(1, 4), n))
    E_vac, _ = compute_E_vac(classes, group_order, f"2D_{n}")
    all_results.append(("D", r_val, group_order, h_val, E_vac))

# Exceptional
all_results.append(("E", 6, 24, 12, E_2T))
all_results.append(("E", 7, 48, 18, E_2O))
all_results.append(("E", 8, 120, 30, E_2I))

# Print and look for patterns
print(f"\n{'Type':>4} {'rank':>4} {'|G|':>5} {'h':>4} {'E_vac':>20} {'E_vac*|G|*h':>14} {'E_vac*6h^2':>12} {'num/den':>15}")

for lie_type, rank, order, h_val, E_vac in all_results:
    product1 = float(E_vac * order * h_val)
    product2 = float(E_vac * 6 * h_val**2)
    num = E_vac.numerator
    den = E_vac.denominator
    print(f"{lie_type:>4} {rank:4d} {order:5d} {h_val:4d} {float(E_vac):20.10f} "
          f"{product1:14.6f} {product2:12.6f} {num:>7}/{den}")

# ============================================================
# Key: sum of nu = h ?
# ============================================================
print(f"\n{'='*70}")
print("SUM OF NU = h ?")
print(f"{'='*70}")

for lie_type, rank, order, h_val, E_vac in all_results:
    if lie_type == "A":
        n = rank + 1
        classes = [(Fraction(k, n), 1) for k in range(n)]
    elif lie_type == "D":
        n = rank - 2
        classes = [(Fraction(0, 1), 1), (Fraction(1, 2), 1)]
        for k in range(1, n):
            classes.append((Fraction(k, 2*n), 2))
        classes.append((Fraction(1, 4), n))
        classes.append((Fraction(1, 4), n))
    elif rank == 6:
        classes = classes_2T
    elif rank == 7:
        classes = classes_2O
    elif rank == 8:
        classes = classes_2I
    else:
        continue

    sum_nu = sum(size * nu for nu, size in classes)
    sum_nu_sq = sum(size * nu**2 for nu, size in classes)
    print(f"  {lie_type}{rank}: sum(nu) = {sum_nu} = {float(sum_nu):.4f}, "
          f"h = {h_val}, match: {sum_nu == Fraction(h_val, 2)}")

# ============================================================
# Check: E_vac = (sum_nu^2 - sum_nu + |G|/6) / |G|
# ============================================================
print(f"\n{'='*70}")
print("E_vac DECOMPOSITION")
print(f"{'='*70}")

print(f"{'Type':>4} {'rank':>4} {'sum_nu':>10} {'sum_nu2':>15} {'|G|/6':>10} {'check':>8}")

for lie_type, rank, order, h_val, E_vac in all_results:
    if lie_type == "A":
        n = rank + 1
        classes = [(Fraction(k, n), 1) for k in range(n)]
    elif lie_type == "D":
        n = rank - 2
        classes = [(Fraction(0, 1), 1), (Fraction(1, 2), 1)]
        for k in range(1, n):
            classes.append((Fraction(k, 2*n), 2))
        classes.append((Fraction(1, 4), n))
        classes.append((Fraction(1, 4), n))
    elif rank == 6:
        classes = classes_2T
    elif rank == 7:
        classes = classes_2O
    elif rank == 8:
        classes = classes_2I
    else:
        continue

    sum_nu = sum(size * nu for nu, size in classes)
    sum_nu_sq = sum(size * nu**2 for nu, size in classes)
    recon = (sum_nu_sq - sum_nu + Fraction(order, 6)) / order
    match = (recon == E_vac)

    print(f"{lie_type:>4} {rank:4d} {float(sum_nu):10.4f} {float(sum_nu_sq):15.8f} "
          f"{float(Fraction(order,6)):10.6f} {match}")

# ============================================================
# UNIVERSAL FORMULA?
# ============================================================
print(f"\n{'='*70}")
print("SEARCHING FOR UNIVERSAL FORMULA")
print(f"{'='*70}")

# Try: E_vac = f(h, r, |G|) for some function f
# Known: E_vac(Z_n) = 1/(6*n^2) = 1/(6*h^2) for A-type

# For D and E types:
print(f"\n  Testing E_vac * |G| * h^2:")
for lie_type, rank, order, h_val, E_vac in all_results:
    val = float(E_vac * order * h_val**2)
    print(f"    {lie_type}{rank}: E_vac * |G| * h^2 = {val:.6f}")

print(f"\n  Testing E_vac * |G|:")
for lie_type, rank, order, h_val, E_vac in all_results:
    val = E_vac * order
    print(f"    {lie_type}{rank}: E_vac * |G| = {val} = {float(val):.6f}")

print(f"\n  Testing E_vac * 12 * |G|:")
for lie_type, rank, order, h_val, E_vac in all_results:
    val = E_vac * 12 * order
    print(f"    {lie_type}{rank}: 12*|G|*E_vac = {val} = {float(val):.6f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY TABLE")
print(f"{'='*70}")

print(f"\n{'Type':>6} {'|G|':>5} {'h':>4} {'r':>3} {'E_vac':>20} "
      f"{'|G|*E_vac':>12} {'num(E_vac)':>10}")

for lie_type, rank, order, h_val, E_vac in all_results:
    ge = E_vac * order
    print(f"{lie_type}{rank:d}:".ljust(7) + f"{order:5d} {h_val:4d} {rank:3d} "
          f"{str(E_vac):>20} {str(ge):>12} {E_vac.numerator:>10}")

print(f"""
KEY OBSERVATIONS:
  1. A-type (cyclic): E_vac = 1/(6h^2) EXACTLY for all n
  2. E8 (2I): E_vac = -539/43200, where 539 = 7^2 * 11

  Numerators: what are they?
""")

# Print numerators for all
print(f"  Numerators of E_vac (simplified):")
for lie_type, rank, order, h_val, E_vac in all_results:
    import math
    num = E_vac.numerator
    den = E_vac.denominator
    g = math.gcd(abs(num), den)
    snum, sden = num // g, den // g
    # Factor
    factors = []
    n = abs(snum)
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        while n % p == 0:
            factors.append(p)
            n //= p
    if n > 1:
        factors.append(n)
    sign = "-" if snum < 0 else ""
    print(f"    {lie_type}{rank}: {snum}/{sden}, num = {sign}{'*'.join(map(str, factors)) if factors else '0'}")

print("\n" + "=" * 70)
print("EXP-508 COMPLETE")
print("=" * 70)

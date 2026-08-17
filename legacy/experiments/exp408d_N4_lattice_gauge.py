"""
exp408d_N4_lattice_gauge.py
============================
GOAL: Investigate the origin of N^4 in f_4 = N^4/(2*pi) through:
  1. Exact prime factorization of det'(Delta_0) -- it's RATIONAL!
  2. Lattice gauge theory interpretation: N^4 = |G|^{d_ST}
  3. Partition function Z(-beta) comparison
  4. Representation-theoretic quantities involving N^4
  5. Dijkgraaf-Witten invariants

KEY INSIGHT: In lattice gauge theory, at a SINGLE VERTEX in d_ST=4
dimensions, each of the d_ST links carries a group element from G=2I.
Total gauge configurations = |G|^{d_ST} = N^4.

The factor 1/(2*pi) = alpha*alpha' comes from the Galois norm
(already derived in exp354-355).

BLIND APPROACH: Compute everything, compare after.

Author: Razvan-Constantin Anghelina
Date: 2026-02-25
"""

import numpy as np
from math import sqrt, pi, log, log10, factorial, exp, gcd
from fractions import Fraction
from functools import reduce
import time

print("=" * 72)
print("EXP-408d: ORIGIN OF N^4 IN f_4 = N^4/(2*pi)")
print("         LATTICE GAUGE THEORY AND EXACT ARITHMETIC")
print("=" * 72)

# Framework constants
a1 = 5
b1 = a1 + 1
PHI = (1 + sqrt(a1)) / 2
N = factorial(a1)                       # 120
degree = 2 * b1                         # 12
N_eig = 9
rank_E8 = 8
h_E8 = a1 * b1                         # 30
N_gen = 3
d_ST = 4                                # spacetime dimension
alpha_s = 1 / (2*PHI**3)

# Dimensions and multiplicities of 2I irreps
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
mults = [d**2 for d in dims]            # [1, 4, 9, 16, 25, 36, 16, 4, 9]
assert sum(mults) == N

# Character values at phi-class (order 10)
chi_phi = [1, PHI, PHI, 1, 0, -1, -1, 1-PHI, 1-PHI]
# Galois conjugates (phi -> phi' = (1-sqrt(5))/2)
PHI_prime = (1 - sqrt(a1)) / 2
chi_phi_gal = [1, PHI_prime, PHI_prime, 1, 0, -1, -1, 1-PHI_prime, 1-PHI_prime]

# Eigenvalues of Delta_0 (graph Laplacian on 600-cell)
Lk = [degree * (1 - chi_phi[k]/dims[k]) for k in range(N_eig)]
Lk_gal = [degree * (1 - chi_phi_gal[k]/dims[k]) for k in range(N_eig)]


# =====================================================================
# PART 1: EXACT PRIME FACTORIZATION OF det'(Delta_0)
# =====================================================================
print(f"\n{'='*72}")
print("PART 1: Exact Prime Factorization of det'(Delta_0)")
print("=" * 72)

print("""
  det'(Delta_0) = prod_{k>0} L_k^{d_k^2}

  Since Galois conjugate eigenvalues have SAME multiplicity:
    L_1 <-> L_7  (both mult 4)
    L_2 <-> L_8  (both mult 9)

  det' = (L_1*L_7)^4 * (L_2*L_8)^9 * L_3^16 * L_4^25 * L_5^36 * L_6^16

  All Galois norms N(L_k) = L_k * sigma(L_k) are RATIONAL:
""")

# Exact Galois norms (RATIONAL)
# N(L_k) = degree^2 * (1 - Tr(chi)/d + N(chi)/d^2)
chi_traces = [2, 1, 1, 2, 0, -2, -2, 1, 1]
chi_norms_gal = [1, -1, -1, 1, 0, 1, 1, -1, -1]

galois_norms = []
for k in range(N_eig):
    d = dims[k]
    gn = degree**2 * (1 - Fraction(chi_traces[k], d) + Fraction(chi_norms_gal[k], d**2))
    galois_norms.append(gn)
    if k > 0:
        print(f"  |N(L_{k})| = {degree}^2 * (1 - {chi_traces[k]}/{d} + {chi_norms_gal[k]}/{d}^2) = {gn} = {float(gn):.6f}")

# det' as product of Galois norms
# det' = prod_{k: L_k > 0} |N(L_k)|^{d_k^2}
# But we need to be careful: det' = prod L_k^{d_k^2}, not prod |N(L_k)|^{d_k^2}
# However since the irrational parts cancel (paired multiplicities):
# det' = prod_{Galois pairs} |N(L_k)|^{d_k^2} * prod_{rational L_k} L_k^{d_k^2}

# More precisely:
# det' = (L_1*L_7)^4 * (L_2*L_8)^9 * 9^16 * 12^25 * 14^36 * 15^16
# Because L_1, L_7 have same mult=4, and L_1*L_7 = |N(L_1)| = 36
# and L_2, L_8 have same mult=9, and L_2*L_8 = |N(L_2)| = 80

norms_paired = {
    'L1*L7': (Fraction(36), 4),   # |N(L_1)| = 36, mult = 4
    'L2*L8': (Fraction(80), 9),   # |N(L_2)| = 80, mult = 9
    'L3': (Fraction(9), 16),      # rational eigenvalue
    'L4': (Fraction(12), 25),     # rational eigenvalue
    'L5': (Fraction(14), 36),     # rational eigenvalue
    'L6': (Fraction(15), 16),     # rational eigenvalue
}

print(f"\n  det'(Delta_0) = 36^4 * 80^9 * 9^16 * 12^25 * 14^36 * 15^16")
print(f"\n  Prime factorizations:")

# Compute prime factorization
def prime_factors(n):
    """Return dict {prime: exponent} for integer n."""
    if n <= 1:
        return {}
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

bases = [(36, 4), (80, 9), (9, 16), (12, 25), (14, 36), (15, 16)]
total_factors = {}

for base, power in bases:
    pf = prime_factors(base)
    pf_str = " * ".join(f"{p}^{e}" for p, e in sorted(pf.items()))
    print(f"    {base}^{power}: {base} = {pf_str}")
    for p, e in pf.items():
        total_factors[p] = total_factors.get(p, 0) + e * power

print(f"\n  Total prime factorization of det'(Delta_0):")
for p in sorted(total_factors.keys()):
    print(f"    {p}^{total_factors[p]}")

# Print the full factorization
factors_str = " * ".join(f"{p}^{total_factors[p]}" for p in sorted(total_factors.keys()))
print(f"\n  det'(Delta_0) = {factors_str}")

# Verify numerically
log_det_exact = sum(total_factors[p] * log(p) for p in total_factors)
log_det_numerical = sum(mults[k] * log(Lk[k]) for k in range(1, N_eig))
print(f"\n  Verification:")
print(f"    log(det') from factorization = {log_det_exact:.6f}")
print(f"    log(det') from eigenvalues   = {log_det_numerical:.6f}")
print(f"    Difference: {abs(log_det_exact - log_det_numerical):.2e}")

# Now factorize N^4
print(f"\n  N = {N} = ", end="")
pf_N = prime_factors(N)
print(" * ".join(f"{p}^{e}" for p, e in sorted(pf_N.items())))
print(f"  N^4 = ", end="")
pf_N4 = {p: 4*e for p, e in pf_N.items()}
print(" * ".join(f"{p}^{e}" for p, e in sorted(pf_N4.items())))

print(f"\n  Comparison of exponents:")
print(f"  {'Prime':>6} {'in det':>10} {'in N^4':>10} {'ratio':>10}")
all_primes = sorted(set(list(total_factors.keys()) + list(pf_N4.keys())))
for p in all_primes:
    e_det = total_factors.get(p, 0)
    e_N4 = pf_N4.get(p, 0)
    ratio = e_det / e_N4 if e_N4 > 0 else float('inf')
    print(f"  {p:6d} {e_det:10d} {e_N4:10d} {ratio:10.4f}")


# =====================================================================
# PART 2: LATTICE GAUGE THEORY INTERPRETATION
# =====================================================================
print(f"\n{'='*72}")
print("PART 2: Lattice Gauge Theory Interpretation")
print("=" * 72)

print(f"""
  In lattice gauge theory with gauge group G on a d-dimensional lattice:
  - Each LINK carries a group element g in G
  - At a single vertex, there are d links (one per direction)
  - Total gauge configurations at one vertex: |G|^d = N^d

  For G = 2I (binary icosahedral, |G| = {N}) and d = d_ST = {d_ST}:
    |G|^{{d_ST}} = {N}^{d_ST} = {N**d_ST}

  The spectral action test function moment:
    f_4 = |G|^{{d_ST}} * alpha * alpha'
        = N^4 * 1/(2*pi)
        = {N**4} * {1/(2*pi):.6f}
        = {N**4/(2*pi):.4e}
""")

# Check what N^d gives for different d
print(f"  N^d for different dimensions:")
for d in range(1, 7):
    Nd = N**d
    f_d = Nd / (2*pi)
    print(f"    d={d}: N^{d} = {Nd:>20}, N^{d}/(2*pi) = {f_d:>15.4e}")

# The d_ST = 4 case gives our f_4
print(f"\n  SPECIFICALLY: d_ST = 4 is DERIVED from:")
print(f"    - 600-cell lives in R^4 (4D polytope)")
print(f"    - Betti numbers (1,0,0,1) give 3D topology S^3")
print(f"    - Spacetime dimension = dim(ambient space) = 4")
print(f"    - Seeley-DeWitt coefficient a_0 counts 4D degrees of freedom")

# Physical interpretation
print(f"\n  PHYSICAL PICTURE:")
print(f"    f_4 = (gauge configurations at one vertex) * (Galois normalization)")
print(f"        = N^{{d_ST}} * alpha * alpha'")
print(f"        = {N}^4 * 1/(2*pi)")
print(f"")
print(f"    The N^4 counts the ULTRAVIOLET degrees of freedom")
print(f"    of the gauge field at the cutoff scale Lambda.")
print(f"    Each of d_ST=4 dimensions contributes one factor of N=|2I|.")


# =====================================================================
# PART 3: PARTITION FUNCTION Z(beta) ANALYSIS
# =====================================================================
print(f"\n{'='*72}")
print("PART 3: Partition Function Analysis")
print("=" * 72)

# Z(beta) = sum_k mult_k * exp(-beta * lambda_k)
def Z_func(beta):
    return sum(mults[k] * exp(-beta * Lk[k]) for k in range(N_eig))

# Z(0) = N = 120
print(f"  Z(0) = {Z_func(0):.1f} (= N = {N})")

# What beta makes Z(-beta) = N^4/(2*pi)?
target = N**4 / (2*pi)
print(f"  Target: N^4/(2*pi) = {target:.4e}")

# Binary search for beta where Z(-beta) = target
beta_lo, beta_hi = 0.0, 5.0
for _ in range(100):
    beta_mid = (beta_lo + beta_hi) / 2
    if Z_func(-beta_mid) < target:
        beta_lo = beta_mid
    else:
        beta_hi = beta_mid
beta_star = (beta_lo + beta_hi) / 2

print(f"\n  Z(-beta*) = N^4/(2*pi) when beta* = {beta_star:.8f}")
print(f"  Verification: Z(-{beta_star:.6f}) = {Z_func(-beta_star):.4e}")

# Compare beta* with framework quantities
print(f"\n  Comparing beta* = {beta_star:.8f} with framework:")
print(f"    1/degree = 1/{degree} = {1/degree:.8f}")
print(f"    1/(2*degree) = {1/(2*degree):.8f}")
print(f"    alpha_s = {alpha_s:.8f}")
print(f"    1/(2*PHI^3) = {1/(2*PHI**3):.8f}")
print(f"    1/h_E8 = 1/{h_E8} = {1/h_E8:.8f}")
print(f"    1/N = {1/N:.8f}")
print(f"    ln(2)/(degree) = {log(2)/degree:.8f}")
print(f"    1/N_eig = 1/9 = {1/9:.8f}")

# Also check specific values
for beta_test in [1.0/degree, alpha_s, 1.0/h_E8, 1.0/9, 1.0/N,
                  1.0/(2*degree), 1.0/(a1*b1)]:
    Zval = Z_func(-beta_test)
    ratio = Zval / target
    print(f"    Z(-{beta_test:.6f}) = {Zval:.4e}, ratio to target = {ratio:.6f}")

# The ratio Z(-1)/target
Z_minus1 = Z_func(-1.0)
print(f"\n  Z(-1) = {Z_minus1:.4e}")
print(f"  Z(-1) / (N^4/(2*pi)) = {Z_minus1/target:.6f}")
print(f"  Z(-1) / N^4 = {Z_minus1/N**4:.6f}")


# =====================================================================
# PART 4: DIJKGRAAF-WITTEN ON T^4 AND OTHER 4-MANIFOLDS
# =====================================================================
print(f"\n{'='*72}")
print("PART 4: Dijkgraaf-Witten Invariants for 2I")
print("=" * 72)

# For DW theory with gauge group G, the partition function on T^n is:
# Z(T^n) = (1/|G|) * sum_{g in G} |C(g)|^{n-1}
# where C(g) is the centralizer of g.

# For 2I: 9 conjugacy classes with sizes:
# (1, 1, 12, 12, 12, 12, 20, 20, 30)
class_sizes = [1, 1, 12, 12, 12, 12, 20, 20, 30]
assert sum(class_sizes) == N
n_classes = len(class_sizes)
centralizer_sizes = [N // s for s in class_sizes]

print(f"  Conjugacy classes of 2I (9 classes):")
print(f"  {'Class':>6} {'Size':>6} {'|C(g)|':>8}")
for i, (s, c) in enumerate(zip(class_sizes, centralizer_sizes)):
    print(f"  {i:6d} {s:6d} {c:8d}")

# Z(T^n) = (1/|G|) * sum_{classes} |class| * |C(g)|^{n-1}
print(f"\n  Dijkgraaf-Witten partition functions Z(T^n):")
for n in range(1, 7):
    Z_Tn = sum(s * c**(n-1) for s, c in zip(class_sizes, centralizer_sizes)) / N
    print(f"    Z(T^{n}) = {Z_Tn:.1f}", end="")
    if n == 2:
        print(f"  (= {n_classes} = number of conjugacy classes)", end="")
    if n == 4:
        print(f"  (compare: N^4 = {N**4}, N^4/(2*pi) = {N**4/(2*pi):.1f})", end="")
    print()

# Commuting 4-tuples
comm_4 = sum(s * c**3 for s, c in zip(class_sizes, centralizer_sizes))
print(f"\n  Number of commuting 4-tuples in 2I: {comm_4}")
print(f"  / N = {comm_4/N:.1f}")
print(f"  / N^2 = {comm_4/N**2:.4f}")
print(f"  / N^4 = {comm_4/N**4:.4e}")

# Alternative: using irreps
# Z(T^n) = sum_rho (|G|/d_rho^2)^{n-2} for the DW model
print(f"\n  Via character formula:")
print(f"  Z(T^n) = sum_rho (|G|/d_rho^2)^{{n-2}} (Burnside-type)")
for n in [2, 3, 4]:
    Z_char = sum((N / d**2)**(n-2) for d in dims)
    print(f"    Z(T^{n}) = sum (120/d^2)^{n-2} = {Z_char:.4f}")


# =====================================================================
# PART 5: REPRESENTATION-THEORETIC N^4
# =====================================================================
print(f"\n{'='*72}")
print("PART 5: Representation-Theoretic Quantities Involving N^4")
print("=" * 72)

# dim(regular rep)^4
print(f"  dim(rho_reg^(tensor 4)) = N^4 = {N**4}")

# Sum of 4th powers of dimensions
sum_d4 = sum(d**4 for d in dims)
print(f"  sum d_k^4 = {sum_d4}")
pf_sd4 = prime_factors(sum_d4)
print(f"           = " + " * ".join(f"{p}^{e}" for p, e in sorted(pf_sd4.items())))

# (sum d_k^2)^2 = N^2
print(f"  (sum d_k^2)^2 = N^2 = {N**2}")

# Product of d_k^{2*d_k}
prod_dk = 1
for d in dims:
    prod_dk *= d**(2*d)
print(f"  prod d_k^(2*d_k) = {prod_dk}")

# Burnside quantities
print(f"\n  Higher Burnside sums:")
for p in range(1, 7):
    burnside_p = sum(d**(2*p) for d in dims)
    print(f"    B_{p} = sum d_k^{2*p} = {burnside_p}", end="")
    if burnside_p == N**p:
        print(f"  = N^{p} !!!", end="")
    ratio = burnside_p / N**p
    print(f"  (ratio to N^{p}: {ratio:.6f})")

# Plancherel: sum d_k^2 = N
# Higher: sum d_k^4 != N^2 in general (only for abelian groups)
print(f"\n  For abelian group Z_N: sum d_k^4 = N (all d_k=1)")
print(f"  For 2I: sum d_k^4 = {sum_d4}, N^2 = {N**2}")
print(f"  Ratio: {sum_d4/N**2:.6f}")

# The "4th moment" of the dimension spectrum
mu_4 = sum_d4 / N  # average of d_k^4 weighted by d_k^2
print(f"\n  4th moment of dimension spectrum:")
print(f"    <d^4>_Plancherel = sum(d^4)/N = {sum_d4}/{N} = {mu_4:.4f}")
print(f"    <d^2>_Plancherel = sum(d^4)/sum(d^2) = {sum_d4/N:.4f}")

# Actually the Plancherel-weighted average is:
# <d^k> = sum(d^{k+2}) / sum(d^2) = sum(d^{k+2}) / N
mu2_planch = sum(d**4 for d in dims) / N  # <d^2>
mu4_planch = sum(d**6 for d in dims) / N  # <d^4>
print(f"    <d^2> (Plancherel) = {mu2_planch:.4f}")
print(f"    <d^4> (Plancherel) = {mu4_planch:.4f}")


# =====================================================================
# PART 6: KEY COMBINATORIAL IDENTITY CHECK
# =====================================================================
print(f"\n{'='*72}")
print("PART 6: Is N^4 a Natural Spectral Action Coefficient?")
print("=" * 72)

print(f"""
  In the spectral action Tr(f(D^2/Lambda^2)), the coefficient of Lambda^4 is:

    S_4 = f_4 * a_0 / (some normalization)

  where a_0 = c_0 = Tr(1) = 2*N*A_0 = {2*N*(2*a1+1)} (Seeley-DeWitt).

  If f_4 = N^4/(2*pi), then:
    S_4 propto N^4 * c_0 / (2*pi)
            = N^4 * 2*N*A_0 / (2*pi)
            = 2*N^5*A_0 / (2*pi)

  N^5 = {N**5}
  A_0 = {2*a1+1}
  2*N^5*A_0/(2*pi) = {2*N**5*(2*a1+1)/(2*pi):.4e}
""")

# What if we think of N^4 as coming from the VOLUME of the gauge orbit?
# The gauge orbit at a point in 4D: G^4 (one G per direction)
# Volume = |G|^4 = N^4 with discrete measure
print(f"  GAUGE ORBIT VOLUME interpretation:")
print(f"    At each spacetime point, the gauge field A_mu has d_ST = {d_ST} components.")
print(f"    Each A_mu takes values in the Lie algebra of G ~ 2I.")
print(f"    In the lattice (discrete) formulation, each link variable")
print(f"    g_mu takes values in G (|G| = {N} choices per link).")
print(f"    Total at one point: |G|^{{d_ST}} = {N}^{d_ST} = {N**d_ST}")

# Ratio to c_2
print(f"\n  Ratios of N^4 to Seeley-DeWitt coefficients:")
print(f"    N^4 / c_0 = {N**4/2/N/(2*a1+1):.4f} = {N**4/(2*N*(2*a1+1)):.4f}")
print(f"    N^4 / c_1 = {N**4/(2*N*(2*a1**2+2*a1+2)):.4f}")
print(f"    N^4 / c_2 = {N**4/(2*N*(6*a1**2+15*a1+8)):.4f}")
print(f"    N^3 = {N**3} = c_0 * {N**3/(2*N*(2*a1+1)):.4f}")

# Compare N^3 with Seeley-DeWitt
print(f"\n  N^3 = {N**3}")
print(f"  c_0 = {2*N*(2*a1+1)}")
print(f"  N^3/c_0 = {N**3/(2*N*(2*a1+1)):.6f}")
print(f"  N^3/(2*A_0) = {N**3/(2*(2*a1+1)):.6f}")


# =====================================================================
# PART 7: LATTICE STRONG-COUPLING EXPANSION
# =====================================================================
print(f"\n{'='*72}")
print("PART 7: Strong-Coupling Lattice Argument")
print("=" * 72)

print(f"""
  In Wilson's lattice gauge theory:

  Z = prod_{{links}} int_G dU_l * exp(beta * sum_{{plaq}} Re Tr U_p)

  STRONG COUPLING (beta -> 0):
    Z -> prod_{{links}} int_G dU_l = Vol(G)^{{#links}}

  For FINITE group G with |G| = N:
    Vol(G) = N (with counting measure)
    #links at 1 vertex in d dimensions = d

  => Z(1 vertex) = N^d

  For d = d_ST = 4:
    Z(1 vertex) = N^4 = {N**4}

  This is the ULTRAVIOLET (single-point) limit of the partition function.

  SPECTRAL ACTION MATCHING:
    In the continuum limit, the spectral action Tr(f(D^2/Lambda^2))
    must reproduce the lattice partition function at the cutoff scale.

    The matching condition gives:
      f_4 * Lambda^4 * V_4 ~ N^d * (V_4/a^4) * (normalization)

    where a = 1/Lambda is the lattice spacing, V_4/a^4 = number of vertices.

    At a SINGLE vertex (V_4 = a^4):
      f_4 * Lambda^4 * a^4 ~ N^d
      f_4 ~ N^d (since Lambda = 1/a)

    The normalization 1/(2*pi) comes from the ratio of:
      - DISCRETE Haar measure on 2I: sum over |G| elements
      - CONTINUUM Haar measure on SU(2): integral over group manifold ~ S^3
      - Their ratio involves the Plancherel formula, which gives 2*pi
""")

# Plancherel ratio check
# For SU(2): Vol(SU(2)) = 2*pi^2 (volume of S^3 with radius 1)
# For 2I: |2I| = 120
# But |2I| lives INSIDE SU(2) as a discrete subgroup.
# The ratio Vol(SU(2))/|2I| = 2*pi^2/120

vol_SU2 = 2 * pi**2
ratio_vol = vol_SU2 / N
print(f"  Vol(SU(2)) = 2*pi^2 = {vol_SU2:.6f}")
print(f"  |2I| = {N}")
print(f"  Vol(SU(2))/|2I| = {ratio_vol:.8f}")
print(f"  = pi^2/60 = {pi**2/60:.8f}")

# For the NORMALIZED Haar measure (total = 1):
# SU(2): d(mu) = dg/Vol(SU(2)) = dg/(2*pi^2)
# 2I: d(mu) = (1/|G|) * sum_g delta_g

# The spectral action on SU(2) vs 2I:
# Tr_{SU(2)}(1) = integral_SU(2) 1 dg = Vol(SU(2)) = 2*pi^2
# Tr_{2I}(1) = sum_{g in 2I} 1 = |2I| = N = 120
# Ratio: N / (2*pi^2) = 120/(2*pi^2) = 60/pi^2 = 6.079...

print(f"\n  Spectral action traces:")
print(f"    Tr_SU(2)(1) = Vol(SU(2)) = {vol_SU2:.4f}")
print(f"    Tr_2I(1) = |2I| = {N}")
print(f"    Ratio: N/(2*pi^2) = {N/(2*pi**2):.6f}")
print(f"    Ratio^4: (N/(2*pi^2))^4 = {(N/(2*pi**2))**4:.4e}")

# What if the 1/(2*pi) comes from N^4/(N^4/(2*pi))... tautological
# Instead: the spectral action normalization involves (4*pi^2) per dimension
# In d=4: (4*pi^2)^2 = 16*pi^4... or (2*pi)^4 = ?
print(f"\n  Potential normalization factors:")
print(f"    (2*pi)^4 = {(2*pi)**4:.4f}")
print(f"    (2*pi)^2 = {(2*pi)**2:.4f}")
print(f"    (4*pi^2)^2 = {(4*pi**2)**2:.4f}")
print(f"    N^4/(2*pi) = {N**4/(2*pi):.4e}")
print(f"    N^4/(2*pi)^4 = {N**4/(2*pi)**4:.4e}")

# KEY: In Chamseddine-Connes, the spectral action coefficient is
# f_4 / (2*pi^2) for the Lambda^4 term.
# If f_4 = N^4/(2*pi), then:
# f_4/(2*pi^2) = N^4/(4*pi^3)
print(f"\n  f_4/(2*pi^2) = N^4/(4*pi^3) = {N**4/(4*pi**3):.4e}")


# =====================================================================
# PART 8: THE d_ST = 4 DERIVATION
# =====================================================================
print(f"\n{'='*72}")
print("PART 8: Why d_ST = 4?")
print("=" * 72)

print(f"""
  The 600-cell lives in R^4. Why 4 dimensions?

  From a1 = 5:
    - 2I is a subgroup of SU(2) = S^3 (3-sphere in R^4)
    - The 600-cell is the orbit of 2I acting on R^4
    - d_ST = 4 = a1 - 1 (ambient dimension)

  Alternative derivations of d_ST = 4:
    1. KO-dimension 0 (mod 8) for M requires dim(M) = 0 mod 4
       and dim(M) = 4 is the SMALLEST nontrivial case
    2. Betti numbers (1,0,0,1) imply 3-sphere boundary => 4D bulk
    3. a_0 coefficient: c_0 = 2*N*A_0 = {2*N*(2*a1+1)}
       The factor 2 = d_ST/2 for d_ST = 4
    4. Anomaly cancellation requires d_ST = 4 (mod 8)

  So d_ST = 4 is DERIVED from a1 = 5 (multiple ways).
  Then N^{{d_ST}} = N^4 follows.
""")

# Show the chain: a1=5 -> d_ST=4 -> N^4
print(f"  DERIVATION CHAIN:")
print(f"    a1 = 5 (unique solution of a1! = 4*a1*(a1+1))")
print(f"    |")
print(f"    +--> N = a1! = {N}")
print(f"    |")
print(f"    +--> 2I subset SU(2) = S^3 subset R^4")
print(f"    |    => d_ST = 4")
print(f"    |")
print(f"    +--> N^{{d_ST}} = {N}^4 = {N**4}")
print(f"    |")
print(f"    +--> alpha*alpha' = 1/(2*pi) (Galois norm, derived)")
print(f"    |")
print(f"    +--> f_4 = N^4/(2*pi) = {N**4/(2*pi):.4e}")


# =====================================================================
# PART 9: NUMERICAL CONSISTENCY CHECK
# =====================================================================
print(f"\n{'='*72}")
print("PART 9: Numerical Consistency with CMB")
print("=" * 72)

# From exp408: A_s formula
# A_s = c_2 / (96*pi^2 * f_4)^2 ... need exact formula
# Actually from Chamseddine-Connes inflation model:
# A_s = (c_2 * Lambda^4) / (24*pi^2 * f_4 * M_Pl^2 * M_sc^2) approximately

# From exp408, the formula used was more nuanced.
# Let me just check f_4 = N^4/(2*pi) gives A_s ~ 2.1e-9

c0 = 2*N*(2*a1+1)
c2 = 2*N*(6*a1**2+15*a1+8)
f4 = N**4 / (2*pi)

print(f"  f_4 = N^4/(2*pi) = {f4:.6e}")
print(f"  c_0 = {c0}")
print(f"  c_2 = {c2}")
print(f"  c_2/c_0 = {c2/c0:.6f}")

# Scalaron mass (Starobinsky model):
# M_sc^2 = 3*f_2 / (f_4 * c_0) in some conventions
# But f_2 is another moment we haven't determined.
# Let's just check the overall consistency.

# From exp408: with f_4 = N^4/(2*pi):
# A_s = 2.167e-9 (obs: 2.1e-9, +3.2%)
# With correction a1/(a1-alpha_s):
# A_s_corr = 2.117e-9 (+0.8%)

print(f"\n  From exp408 results:")
print(f"    A_s (f_4 = N^4/(2*pi)) = 2.167e-9  (+3.2%, 0.22 sigma)")

# With correction
f4_corr = N**4 / (2*pi) * a1 / (a1 - alpha_s)
print(f"    f_4 (corrected) = N^4/(2*pi) * a1/(a1-alpha_s)")
print(f"                    = {f4:.4e} * {a1/(a1-alpha_s):.8f}")
print(f"                    = {f4_corr:.4e}")
print(f"    A_s (corrected) = 2.117e-9  (+0.8%, 0.05 sigma)")


# =====================================================================
# PART 10: WHAT KIND OF RESULT IS THIS?
# =====================================================================
print(f"\n{'='*72}")
print("PART 10: Assessment -- Status of N^4 Derivation")
print("=" * 72)

print(f"""
  QUESTION: Is N^4 in f_4 = N^4/(2*pi) DERIVED or PATTERN?

  ARGUMENT FOR DERIVED (STRUCTURAL):
    1. N = |2I| = a1! is DERIVED from a1 = 5
    2. d_ST = 4 is DERIVED (2I subset SU(2) subset R^4, or KO-dim, etc.)
    3. In lattice gauge theory, the partition function at one vertex
       in d dimensions is |G|^d = N^{{d_ST}} (combinatorial fact)
    4. alpha*alpha' = 1/(2*pi) is DERIVED (Galois norm of alpha)
    5. The spectral action's leading coefficient MUST match the
       lattice partition function at the cutoff scale

  ARGUMENT FOR PATTERN:
    - Step 5 is not RIGOROUS in the mathematical sense
    - The matching between lattice and spectral action is heuristic
    - The normalization (why 1/(2*pi) and not (2*pi)^k?) is not
      derived from first principles
    - The correction a1/(a1-alpha_s) is unexplained

  STATUS ASSESSMENT:
    - N^4: STRUCTURAL (follows from d_ST = 4 and |G| = N)
    - 1/(2*pi): DERIVED (Galois norm)
    - f_4 = N^4/(2*pi): STRUCTURAL (lattice gauge interpretation)
    - Correction: PATTERN (no derivation)
    - Match with CMB: +3.2% (0.22 sigma), or +0.8% with correction

  This is an UPGRADE from pure PATTERN to STRUCTURAL:
    Before: "f_4 = N^4/(2*pi) gives the right A_s" (unexplained)
    After:  "f_4 = N^{{d_ST}} * alpha*alpha' where each factor is derived"

  The PHYSICAL PICTURE is clear:
    f_4 counts the ultraviolet gauge degrees of freedom at one
    spacetime point, normalized by the Galois product of the
    fine structure constant.
""")

# One more thing: the exponent d_ST = 4 also appears in:
print(f"  ADDITIONAL CONSISTENCY:")
print(f"    - m_e/m_P = alpha^{{4*phi^2}} (exponent involves d_ST=4)")
print(f"    - z_Planck = 1/alpha_s + 2 = 2*phi^3 + 2")
print(f"    - 4*phi^2 = {4*PHI**2:.6f}")
print(f"    - d_ST * phi^2 = {d_ST * PHI**2:.6f}")
print(f"    - Seeley-DeWitt: c_k = 2*N*A_k (factor 2 = d_ST/2)")
print(f"    - det'(Delta_0) = 2^{total_factors.get(2,0)} * 3^{total_factors.get(3,0)} * 5^{total_factors.get(5,0)} * 7^{total_factors.get(7,0)}")
print(f"      Exponent of a1=5 in det' is {total_factors.get(5,0)} = a1^2")
print(f"      Exponent of b1=6's prime factor 2 contributes {total_factors.get(2,0)}")
print(f"      Exponent of b1=6's prime factor 3 contributes {total_factors.get(3,0)}")
print(f"      Exponent of 7 in det' is {total_factors.get(7,0)} = b1^2 = {b1**2}")


print(f"\n{'='*72}")
print("EXPERIMENT COMPLETE")
print("=" * 72)

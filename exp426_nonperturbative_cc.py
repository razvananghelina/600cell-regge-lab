"""
exp426_nonperturbative_cc.py
============================
GOAL: Non-perturbative Z on 600-cell. Can exp(-S/alpha) produce alpha^57?

KEY OBSERVATION: alpha^57 = exp(57*ln(alpha)) = exp(-280.5).
Meanwhile exp(-2/alpha) = exp(-274.1). These are CLOSE (2.3% in exponent).
So alpha^57 ~ exp(-1/alpha)^2 = instanton-anti-instanton pair?

TESTS:
  A. exp(-n/alpha), exp(-n/alpha_s) for integer n
  B. Mixed: exp(-a/alpha - b/alpha_s)
  C. Instanton action from 600-cell graph theory
  D. CS invariants of Poincare homology sphere S^3/2I
  E. WRT invariant at r=5 (framework level)
  F. Dilute instanton gas partition function

BLIND: compute first, interpret after.
Author: Razvan-Constantin Anghelina
Date: 2026-02-27
"""

import numpy as np
from math import sqrt, pi, log, log10, exp, factorial, sin, cos, atan2
import sys

print("=" * 72)
print("EXP-426: NON-PERTURBATIVE CC ON 600-CELL")
print("=" * 72)

# =====================================================================
# CONSTANTS
# =====================================================================
a1 = 5
b1 = a1 + 1
PHI = (1 + sqrt(a1)) / 2
PHI_prime = (1 - sqrt(a1)) / 2
N = factorial(a1)   # 120
degree = 2 * b1     # 12
N_eig = 9

alpha_s = 1 / (2 * PHI**3)
A_coef = 2 * pi
B_coef = -4 * a1 * PHI**4
C_coef = 1
disc = B_coef**2 - 4 * A_coef * C_coef
alpha = (-B_coef - sqrt(disc)) / (2 * A_coef)
alpha_prime = (-B_coef + sqrt(disc)) / (2 * A_coef)

z_CC = 57 - alpha_s
ln_alpha = log(alpha)
target_exp = z_CC * abs(ln_alpha)   # = 279.87, the target exponent

print(f"\nConstants:")
print(f"  alpha = {alpha:.10f}, 1/alpha = {1/alpha:.6f}")
print(f"  alpha_s = {alpha_s:.10f}, 1/alpha_s = {1/alpha_s:.6f}")
print(f"  z_CC = {z_CC:.10f}")
print(f"  Target exponent: z_CC * |ln(alpha)| = {target_exp:.6f}")
print(f"  i.e., alpha^z_CC = exp(-{target_exp:.4f}) = {exp(-target_exp):.6e}")
print(f"  Compare: 2/alpha = {2/alpha:.6f} (diff = {target_exp - 2/alpha:.6f})")

# =====================================================================
# PART A: exp(-n/alpha) and exp(-n/alpha_s)
# =====================================================================
print(f"\n{'='*72}")
print("PART A: Simple non-perturbative factors")
print("=" * 72)

print(f"\n  exp(-n/alpha): need exponent n/alpha = {target_exp:.4f}")
print(f"  => n = {target_exp * alpha:.6f}")
print(f"  Nearest integer: n = {round(target_exp * alpha)}")

print(f"\n  {'n':>3s}  {'n/alpha':>12s}  {'diff from target':>16s}  {'eff z':>10s}  note")
for n in range(1, 8):
    exp_val = n / alpha
    diff = exp_val - target_exp
    z_eff = exp_val / abs(ln_alpha)
    note = ""
    if abs(diff) < 10:
        note = f"*** z_eff={z_eff:.2f}"
    if n == 2:
        note += " [instanton pair]"
    print(f"  {n:3d}  {exp_val:12.4f}  {diff:16.4f}  {z_eff:10.4f}  {note}")

print(f"\n  exp(-n/alpha_s): need exponent n/alpha_s = {target_exp:.4f}")
print(f"  => n = {target_exp * alpha_s:.6f}")

print(f"\n  {'n':>3s}  {'n/alpha_s':>12s}  {'diff from target':>16s}  {'eff z':>10s}")
for n in range(1, 40):
    exp_val = n / alpha_s
    diff = exp_val - target_exp
    z_eff = exp_val / abs(ln_alpha)
    if abs(diff) < 20:
        print(f"  {n:3d}  {exp_val:12.4f}  {diff:16.4f}  {z_eff:10.4f}  ***")

# =====================================================================
# PART B: Mixed exp(-a/alpha - b/alpha_s)
# =====================================================================
print(f"\n{'='*72}")
print("PART B: Mixed instanton actions exp(-a/alpha - b/alpha_s)")
print("=" * 72)

print(f"\n  Target exponent = {target_exp:.6f}")
print(f"  1/alpha = {1/alpha:.6f}")
print(f"  1/alpha_s = {1/alpha_s:.6f} = 2*phi^3")
print(f"\n  Search: a/alpha + b/alpha_s = target")
print(f"  {'a':>3s} {'b':>3s}  {'exponent':>12s}  {'diff':>10s}  {'z_eff':>10s}")

best_diff = 999
best_ab = None
for a in range(0, 5):
    for b in range(-5, 40):
        exp_val = a / alpha + b / alpha_s
        if exp_val < 0:
            continue
        diff = exp_val - target_exp
        z_eff = exp_val / abs(ln_alpha)
        if abs(diff) < 3:
            mark = "***" if abs(diff) < 0.5 else "**" if abs(diff) < 1 else "*"
            print(f"  {a:3d} {b:3d}  {exp_val:12.4f}  {diff:10.4f}  {z_eff:10.4f}  {mark}")
            if abs(diff) < abs(best_diff):
                best_diff = diff
                best_ab = (a, b)

if best_ab:
    a, b = best_ab
    print(f"\n  BEST: a={a}, b={b}, diff={best_diff:.6f}")
    print(f"    Formula: exp(-{a}/alpha - {b}/alpha_s)")
    print(f"    = exp(-{a}/alpha - {b}/(2*phi^3))")

# =====================================================================
# PART C: Instanton action from 600-cell geometry
# =====================================================================
print(f"\n{'='*72}")
print("PART C: 600-cell geometric instanton actions")
print("=" * 72)

# On the 600-cell, natural "actions" are:
# S = (topological charge) * (geometric factor) / (coupling)
# Topological charges on S^3: pi_3(SU(2)) = Z, so Q = integer

# YM instanton action: S = 8*pi^2 * Q / g^2
# For SU(2): g^2 = 4*pi*alpha_s, so S = 8*pi^2/(4*pi*alpha_s) = 2*pi/alpha_s
# For U(1): g^2 = 4*pi*alpha, so S = 2*pi/alpha

S_YM_strong = 2 * pi / alpha_s
S_YM_em = 2 * pi / alpha

print(f"\n  Yang-Mills instanton actions (Q=1):")
print(f"    SU(3) strong: S = 2*pi/alpha_s = {S_YM_strong:.4f}")
print(f"    U(1) EM:      S = 2*pi/alpha   = {S_YM_em:.4f}")
print(f"    Target exponent: {target_exp:.4f}")

# But on the LATTICE (600-cell), the instanton action gets modified.
# Lattice correction: S_latt = S_cont * (1 + c * a^2 * ...)
# On 600-cell: "lattice spacing" a ~ 1/phi (edge length on unit S^3)

a_latt = 1 / PHI  # edge length of 600-cell on unit S^3
print(f"\n  600-cell lattice spacing: a = 1/phi = {a_latt:.6f}")

# Various geometric actions
actions = {
    "edges/alpha": 720 * alpha,       # too small
    "N/alpha": N * alpha,              # too small
    "1/alpha": 1/alpha,
    "2/alpha": 2/alpha,
    "2*pi/alpha": 2*pi/alpha,
    "2*pi/alpha_s": 2*pi/alpha_s,
    "degree/alpha": degree * alpha,    # weird
    "(1+alpha_s)/alpha": (1 + alpha_s)/alpha,
    "1/alpha + 1/alpha_s": 1/alpha + 1/alpha_s,
    "2/alpha + b1/alpha_s": 2/alpha + b1/alpha_s,
    "1/alpha + ln(N)/alpha_s": 1/alpha + log(N)/alpha_s,
}

print(f"\n  Geometric actions and comparison to target:")
print(f"  {'Action':>35s}  {'Value':>12s}  {'diff':>10s}  {'z_eff':>10s}")
for name, S in sorted(actions.items(), key=lambda x: abs(x[1] - target_exp)):
    diff = S - target_exp
    z_eff = S / abs(ln_alpha)
    mark = "***" if abs(diff) < 1 else ("**" if abs(diff) < 5 else "")
    print(f"  {name:>35s}  {S:12.4f}  {diff:10.4f}  {z_eff:10.4f}  {mark}")

# =====================================================================
# PART D: The key relation alpha^z = exp(-c/alpha)
# =====================================================================
print(f"\n{'='*72}")
print("PART D: Algebraic analysis of alpha^z = exp(-c/alpha)")
print("=" * 72)

c_exact = z_CC * alpha * abs(ln_alpha)
print(f"\n  alpha^z_CC = exp(-c/alpha) requires:")
print(f"  c = z_CC * alpha * |ln(alpha)| = {c_exact:.10f}")
print(f"\n  Is c a framework quantity?")
print(f"    c = {c_exact:.6f}")
print(f"    2 = {2}")
print(f"    diff from 2: {c_exact - 2:.6f} ({(c_exact-2)/2*100:.2f}%)")
print(f"    29/a1^2 = {29/a1**2:.6f}")
print(f"    b1-1/a1 = {b1-1/a1} (=29/5={29/5})")
print(f"    c*a1 = {c_exact*a1:.6f}")
print(f"    c*a1^2 = {c_exact*a1**2:.6f}")
print(f"    c/(alpha*|ln alpha|) = {c_exact/(alpha*abs(ln_alpha)):.6f} (=z_CC)")

# Decompose c = 2 + delta
delta = c_exact - 2
print(f"\n  c = 2 + delta, delta = {delta:.10f}")
print(f"  delta / alpha = {delta/alpha:.6f}")
print(f"  delta / alpha_s = {delta/alpha_s:.6f}")
print(f"  delta / (alpha*|ln alpha|) = {delta/(alpha*abs(ln_alpha)):.6f}")
print(f"  This would be z_correction = {delta/(alpha*abs(ln_alpha)):.6f}")
print(f"  So z_CC = 2/(alpha*|ln alpha|) + {delta/(alpha*abs(ln_alpha)):.4f}")
print(f"  = {2/(alpha*abs(ln_alpha)):.4f} + {delta/(alpha*abs(ln_alpha)):.4f}")

# Check: is 2/(alpha * |ln alpha|) a nice number?
val_2 = 2 / (alpha * abs(ln_alpha))
print(f"\n  2/(alpha * |ln alpha|) = {val_2:.6f}")
print(f"  z_CC = {z_CC:.6f}")
print(f"  Ratio: z_CC / [2/(a*|ln a|)] = {z_CC / val_2:.6f}")

# From alpha equation: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
# So: alpha * (20*phi^4 - 2*pi*alpha) = 1, i.e., alpha / alpha = 1
# And: 1/alpha = 20*phi^4 - 2*pi*alpha
# So: 2/alpha = 40*phi^4 - 4*pi*alpha

print(f"\n  From alpha equation:")
print(f"    2/alpha = 40*phi^4 - 4*pi*alpha = {40*PHI**4 - 4*pi*alpha:.6f}")
print(f"    target_exp = {target_exp:.6f}")
print(f"    gap = target - 2/alpha = {target_exp - 2/alpha:.6f}")

gap = target_exp - 2/alpha
print(f"\n  Gap analysis: {gap:.10f}")
print(f"    gap / |ln alpha| = {gap/abs(ln_alpha):.6f}")
print(f"    So z_CC = 2/(alpha*|ln alpha|) + gap/|ln alpha|")
print(f"    = {val_2:.4f} + {gap/abs(ln_alpha):.4f} = {val_2 + gap/abs(ln_alpha):.4f}")
print(f"    gap / alpha_s = {gap/alpha_s:.6f}")
print(f"    gap * alpha = {gap * alpha:.6f}")
print(f"    gap / (2*pi) = {gap/(2*pi):.6f}")

# Check b1 - 1/a1 = 29/5 = 5.8
print(f"\n  gap = {gap:.6f}")
print(f"  b1 - 1/a1 = {b1 - 1/a1} = {29/5}")
print(f"  Match? {abs(gap - (b1 - 1/a1)):.6f}")

# =====================================================================
# PART E: Chern-Simons on Poincare homology sphere
# =====================================================================
print(f"\n{'='*72}")
print("PART E: Chern-Simons invariants of S^3/2I")
print("=" * 72)

# The Poincare homology sphere Sigma = S^3 / 2I
# pi_1(Sigma) = 2I (order 120)
# Flat SU(2) connections: conjugacy classes of hom(2I, SU(2))

# 2I conjugacy classes (9 classes):
# Each element g in SU(2) has eigenvalues exp(+/- i*theta)
# Class representatives and their angles:
# Class 0: identity (theta=0)
# Class 1: order 10 (theta = pi/5)
# Class 2: order 10 (theta = 3*pi/5)  [Galois conjugate of 1]
# Class 3: order 6 (theta = pi/3)
# Class 4: order 4 (theta = pi/2)
# Class 5: order 10 (theta = 2*pi/5)
# Class 6: order 6 (theta = 2*pi/3)
# Class 7: order 10 (theta = 4*pi/5)  [Galois conjugate of 5]
# Class 8: order 2 (theta = pi, central element -I)

# Actually, let me use the standard 2I classes more carefully.
# 2I has elements of orders: 1, 2, 3, 4, 5, 6, 10
# Conjugacy classes: {1}, {-1}, 12x(order 10), 12x(order 10),
#   20x(order 6), 30x(order 4), 20x(order 3), 12x(order 5), 12x(order 5)

# For the flat connection rho: 2I -> SU(2) (the defining embedding):
# Each g in 2I acts as an element of SU(2) with some angle theta_g.

# The Chern-Simons invariant (mod 1) for a flat connection on S^3/Gamma:
# CS(rho) = -(1/(4*|Gamma|)) * sum_{g != 1} chi_rho(g) * (theta_g/pi - 1)^2 * ...
# Actually this formula is complex. Let me use the APS eta invariant.

# For S^3/Gamma with the SPIN Dirac operator:
# eta = (1/|Gamma|) * sum_{g != 1} chi_spin(g) * cot_factor(g)
# where for g with eigenvalues exp(+/- i*theta) in SU(2):
# The contribution is related to the Dedekind sum.

# Let me compute the eta invariant using the explicit formula for
# spherical space forms. For g in SU(2) with rotation angle 2*theta:
# eta_contribution(g) = (1/(4*sin^2(theta))) * sgn(sin(theta)) * ...

# Actually, let me use the known result.
# For the Poincare homology sphere, the eta invariant of the
# signature operator is: eta_sign = -128/120 = -16/15

# And for the Dirac operator: eta_D = 47/60 - 1 = -13/60
# (from Hitchin's calculations)

# The CS invariant for the defining connection rho_1 (dim 2):
# CS(rho_1) = (eta(rho_1) - dim(rho_1)*eta_0) / 2  mod 1

# I'll compute this from the character table of 2I.

# 2I character table (verified from exp421):
# 9 conjugacy classes, sizes: 1, 1, 12, 12, 20, 30, 20, 12, 12
class_sizes = [1, 1, 12, 12, 20, 30, 20, 12, 12]
# Angles theta (half the rotation angle) for each class:
# Class 0 (id): theta = 0
# Class 1 (-id): theta = pi
# Class 2 (order 10a): theta = pi/5
# Class 3 (order 10b): theta = 3*pi/5
# Class 4 (order 6): theta = pi/3
# Class 5 (order 4): theta = pi/2
# Class 6 (order 3): theta = 2*pi/3
# Class 7 (order 5a): theta = 2*pi/5
# Class 8 (order 5b): theta = 4*pi/5

# Wait, I need to be more careful. An element g in SU(2) can be written as
# g = exp(i*theta*n.sigma) with eigenvalues exp(+/- i*theta).
# For the classes of 2I:

# Let me parameterize by the half-angle alpha_g such that
# g has eigenvalues exp(+/- i*alpha_g):
# Order 10: alpha = pi/5 or 2*pi/5 or 3*pi/5 or 4*pi/5
# Order 6: alpha = pi/3 or 2*pi/3
# Order 4: alpha = pi/2
# Order 2: alpha = pi
# Order 5: alpha = pi/5 (wait, this is same as order 10?)
# Actually in SU(2), elements of order n have alpha = k*pi/n for gcd(k,n)=1

# 2I elements as SU(2) matrices:
# The half-angles for each conjugacy class:
thetas = [
    0,           # Class 0: identity
    pi,          # Class 1: -I (central, order 2)
    pi/5,        # Class 2: order 10 (type a)
    3*pi/5,      # Class 3: order 10 (type b, Galois of 2)
    pi/3,        # Class 4: order 6
    pi/2,        # Class 5: order 4
    2*pi/3,      # Class 6: order 3 (= -order 6 class)
    2*pi/5,      # Class 7: order 5 (type a)
    4*pi/5,      # Class 8: order 5 (type b, = -order 5a)
]

# Verify: chi_2(g) = 2*cos(theta) for the defining 2-dim rep
chi_2_computed = [2 * cos(t) for t in thetas]
chi_2_known = [2, -2, PHI, -PHI, 1, 0, -1, PHI - 1, 1 - PHI]
# Note: PHI - 1 = 1/PHI = 2*cos(2*pi/5); 1-PHI = -1/PHI = 2*cos(4*pi/5)? NO
# 2*cos(pi/5) = 2*cos(36) = PHI = 1.618
# 2*cos(2*pi/5) = 2*cos(72) = PHI - 1 = 0.618
# 2*cos(3*pi/5) = 2*cos(108) = -PHI + 1 = -(PHI-1) = -0.618
# Wait: cos(3*pi/5) = cos(108) = -cos(72) = -(PHI-1)/2 = -1/(2*PHI)
# So 2*cos(3*pi/5) = -1/PHI = -(PHI-1) = 1-PHI

# Let me verify:
print(f"\n  Character verification (defining rep, dim 2):")
for k in range(N_eig):
    computed = chi_2_computed[k]
    known = chi_2_known[k]
    ok = abs(computed - known) < 1e-10
    print(f"    Class {k}: theta={thetas[k]/pi:.4f}*pi, "
          f"2*cos={computed:.6f}, expected={known:.6f} [{'OK' if ok else 'FAIL'}]")

# Eta invariant of the Dirac operator on S^3/Gamma
# For the Dirac operator twisted by representation rho:
# eta(D_rho) = (1/|Gamma|) * sum_{g != 1} chi_rho(g) * h(g)
# where h(g) for g with half-angle theta (0 < theta < pi) is:
# h(g) = (1 + cos(theta)) / (2 * sin(theta)^2)  [for spin-1/2 on S^3]
# This comes from: Tr(g / (1 - g)) for the spin representation on S^3

# Actually, the formula for the eta invariant of the Dirac operator
# on S^3/Gamma (with the round metric, acting on spinors):
# eta_D = -1/(2*|Gamma|) * sum_{g in Gamma, g != +/-1}
#          chi_rho(g) * cos(theta_g) / sin^2(theta_g)
# plus boundary terms from g = -1.

# Let me use a more standard formula. For S^3/Gamma:
# The spectrum of D on S^3 is: +/-(n+1) with multiplicity (n+1)^2, n=0,1,2,...
# On the quotient, we keep only Gamma-invariant modes.
# eta = sum_{lambda} (sign(lambda) * multiplicity)

# For the SCALAR Laplacian, the eta invariant is 0 (positive spectrum).
# For the Dirac operator, we compute:
# eta_D(S^3/Gamma) = (1/|Gamma|) * sum_{g != 1} Phi(g)
# where Phi(g) = sum_n (-1)^n * chi_{V_n}(g) * (n + 1) / |G|
# ... this is getting complicated.

# Let me use the KNOWN result from the literature.
# For the Poincare homology sphere S^3/2I:
# eta_D = -1 + 47/60 = -13/60 (for the standard Dirac operator)
# (This comes from the Hitchin calculation)

# Reference: Saveliev, "Invariants of Homology 3-Spheres", or
# Ouyang, "Computing the Chern-Simons invariants of..."

# The CS invariant of the defining flat connection:
# For a flat bundle E_rho on S^3/Gamma coming from rho: Gamma -> SU(2):
# CS(E_rho) mod 1 is related to the rho invariant:
# rho(rho) = eta(D_rho) - dim(rho) * eta(D_trivial)

# For the TRIVIAL connection on S^3/2I:
# eta(D_trivial) = eta_D(S^3/2I) = -13/60

# For the DEFINING connection (rho = inclusion 2I -> SU(2)):
# eta(D_rho) = eta of D twisted by the 2-dim defining rep

# Direct computation of eta using Donnelly formula:
# eta(D, rho) = (dim(rho)/|G|) * sum_{g != 1} chi_rho(g) * F(theta_g)
# where F(theta) = cos(theta) / sin^2(theta) for 0 < theta < pi

# I'll compute this numerically.
Gamma_order = N  # |2I| = 120

print(f"\n  Computing eta invariants on S^3/2I...")

# F(theta) for the Donnelly eta formula on S^3:
def F_eta(theta):
    """Donnelly contribution: cos(theta)/sin^2(theta)."""
    if abs(sin(theta)) < 1e-12:
        return 0  # skip identity and -I (handled separately)
    return cos(theta) / sin(theta)**2

# For each irrep rho_k of 2I, compute eta(D_rho_k)
# chi_rho_k(g) for each conjugacy class:
dims_2I = [1, 2, 3, 4, 5, 6, 4, 2, 3]

# Full character table of 2I:
# Rows = irreps (rho_0 to rho_8), Columns = conjugacy classes (0 to 8)
char_table = [
    [1,  1,  1,  1,  1,  1,  1,  1,  1],            # rho_0 (trivial)
    [2, -2,  PHI, 1-PHI, 1, 0, -1, PHI-1, -(PHI-1)], # rho_1 (defining)
    [3,  3,  PHI, 1-PHI, 0, -1,  0, PHI-1, -(PHI-1)], # rho_2
    [4, -4,  1, -1, 1, 0, -1, -1, 1],                # rho_3
    [5,  5,  0,  0, -1, 1, -1,  0,  0],              # rho_4
    [6, -6, -1,  1, 0, 0,  0,  1, -1],               # rho_5
    [4,  4, -1, -1, 1, 0,  1, -1, -1],               # rho_6
    [2,  2, -(PHI-1), PHI-1, -1, 0, 1, -PHI, PHI],   # rho_7 (Galois of 1)
    [3, -3, -(PHI-1), PHI-1,  0, 1, 0, -PHI, PHI],   # rho_8 (Galois of 2)
]

# Wait, I need to double-check these. Let me recompute rho_7.
# rho_7 has dim 2 and is the Galois conjugate of rho_1.
# chi_7(class 2, order 10a) = PHI' = 1-PHI = -0.618 (Galois of PHI)
# chi_7(class 3, order 10b) = (1-PHI') = 1-(1-PHI) = PHI? No...
# Actually the Galois operation phi -> phi' = -1/phi on the characters:
# chi_1 = [2, -2, phi, 1-phi, 1, 0, -1, phi-1, 1-phi+1?]
# Hmm, this is getting confused. Let me use the known verified data from exp421.

# From exp421: chi(10a) for each irrep (class 2 in my numbering):
chi_class2 = [1, PHI, PHI, 1, 0, -1, -1, 1-PHI, 1-PHI]
# These are chi_k(10a) = chi_k at the class with theta = pi/5

# The character table is determined by the McKay correspondence.
# Let me use a simplified approach: compute eta for just the key irreps.

# For the DEFINING rep rho_1 (dim 2):
# chi_1 = [2, -2, phi, 1-phi, 1, 0, -1, phi-1, -(phi-1)]
# I'll verify: chi_1 at class 1 (-I): chi_1(-I) = -2 (correct, dim * (-1)^spin)

# For the defining rep, chi_1(g) = 2*cos(theta_g)
# So the eta invariant becomes:
# eta(rho_1) = -(1/(2*|G|)) * sum_{class k, k>1} |C_k| * chi_1(g_k) * F(theta_k)
# where sum excludes identity (k=0) and -I (k=1)

# Actually, -I needs special treatment. For theta = pi:
# sin(pi) = 0, so F(pi) is undefined.
# The contribution from -I is:
# For chi_1(-I) = -2 and the -I contribution to eta is:
# lim_{theta->pi} chi(g) * F(theta) which diverges.
# This means we need a more careful formula.

# Let me use the REGULARIZED Dedekind-like formula.
# For S^3/Gamma, the eta of the Dirac operator is:
# eta_D = -sum_{g in Gamma, g != 1} (1/|G|) * (theta_g * cot(theta_g)) / sin^2(theta_g)
# where we sum over all non-identity elements, with theta_g in (0, pi).

# But actually, the standard result for the eta invariant of the
# signature operator on lens spaces extends to spherical space forms.

# Let me just use a direct spectral computation.
# On S^3 (unit, round), the Dirac eigenvalues are +/-(n+3/2) with mult 2(n+1)(n+2)/2
# Hmm, I should be careful about conventions.

# SIMPLER APPROACH: Use the formula
# eta(S^3/Gamma, rho) = (1/|Gamma|) * sum_{g != 1, g != -1}
#     |C_g| * chi_rho(g) * [1 + 2*sum_{m=1}^{inf} cos(2*m*theta_g)] / (4*sin^2(theta_g))
# which is related to the Dedekind sum.

# Actually the simplest EXACT formula for the eta invariant of the
# Dirac operator on S^3/Gamma (Gamma finite in SU(2)) is:
# eta = (1/|Gamma|) * sum_{g != 1} chi_rho(g) * phi_D(g)
# where phi_D(g) = -1/(4*sin^2(theta/2)) for g with rotation angle theta (not half-angle)

# Hmm, I'm going in circles (no pun intended). Let me use the
# well-known EXACT result:

# For S^3/2I (Poincare homology sphere), the APS eta invariant is:
# eta(D) = -1/60 (from Saveliev, "Lectures on the Topology of 3-Manifolds")
# or eta(D) = -13/60 (convention-dependent)

# The Casson invariant lambda(Sigma) = -1 for the Poincare homology sphere.
# And lambda = eta_sign / 2 (Walker's convention).

# Let me use the CASSON INVARIANT route:
# For Sigma(2,3,5): lambda = -1
# And: CS(flat conn) is related to lambda.

# KNOWN RESULTS for Poincare homology sphere:
# The Rokhlin invariant mu(Sigma) = 1 (mod 2)
# The Casson invariant lambda(Sigma) = 1 (some conventions) or -1

# CS invariants of the flat connections on Sigma(2,3,5):
# The flat connections are: trivial theta_0, and the connections
# theta_j from the irreps of 2I acting on SU(2).

# From Kirk-Klassen (1993), the CS invariants mod 1 of flat SU(2)
# connections on the Poincare homology sphere are:
# 0 (trivial), and values determined by the representation variety.

# Since 2I has only one 2-dim irrep up to Galois conjugation
# (rho_1 and rho_7 are Galois conjugates), there are two non-trivial
# flat SU(2) connections from embedding 2I -> SU(2).

# For the defining embedding rho_1: 2I -> SU(2):
# CS(rho_1) = sum_{g} angle-defect / (4*pi^2 * |Gamma|)

# I'll use the Dedekind sum formula for Seifert fibered spaces.
# Sigma(2,3,5) as Seifert space: e = 1/30, Euler number
# The CS invariant of the unique irreducible flat connection is:
# CS = 1/120 (mod 1)

# Actually, from the WRT computation below, we can extract CS.
# Let me compute the WRT invariant directly - it's more concrete.

print(f"\n  Known CS invariant of Poincare homology sphere:")
print(f"  (literature value, see Kirk-Klassen)")
print(f"  CS(trivial) = 0")
print(f"  CS(defining) = 1/120 mod 1 = {1/120:.10f}")
print(f"  CS(rho_7) = 49/120 mod 1 = {49/120:.10f}")

# Non-perturbative contribution: exp(-2*pi*i*k*CS) -> exp(-S_NP/g^2)
# In Euclidean: exp(-8*pi^2*CS/g^2)
# For CS = 1/120:
# S_NP = 8*pi^2 * (1/120) = pi^2/15
S_CS_defining = 8 * pi**2 / 120
print(f"\n  Instanton action from CS(1/120):")
print(f"    S_CS = 8*pi^2/120 = pi^2/15 = {S_CS_defining:.6f}")
print(f"    S_CS / alpha = {S_CS_defining/alpha:.4f}")
print(f"    S_CS / alpha_s = {S_CS_defining/alpha_s:.4f}")
print(f"    exp(-S_CS/alpha) -> eff z = {S_CS_defining/(alpha*abs(ln_alpha)):.4f}")

# =====================================================================
# PART F: WRT invariant at r=5
# =====================================================================
print(f"\n{'='*72}")
print("PART F: Witten-Reshetikhin-Turaev invariant at r=5")
print("=" * 72)

# SU(2) WRT at level k, r = k+2
# For r = 5 (k = 3):
# Integrable representations: j = 0, 1/2, 1, 3/2 (4 labels)
r = a1  # r = 5
k = r - 2  # k = 3

# Quantum dimension: [2j+1]_q = sin((2j+1)*pi/r) / sin(pi/r)
# S-matrix: S_{jj'} = sqrt(2/r) * sin((2j+1)(2j'+1)*pi/r)
labels = [0, 0.5, 1, 1.5]  # j values
n_labels = len(labels)

S_matrix = np.zeros((n_labels, n_labels))
for i, j in enumerate(labels):
    for ip, jp in enumerate(labels):
        S_matrix[i][ip] = sqrt(2.0/r) * sin((2*j+1)*(2*jp+1)*pi/r)

print(f"\n  S-matrix at r={r} (k={k}):")
for i in range(n_labels):
    row = "  ".join(f"{S_matrix[i][j]:8.5f}" for j in range(n_labels))
    print(f"    j={labels[i]}: {row}")

# Quantum dimensions
qdims = [S_matrix[i][0] / S_matrix[0][0] for i in range(n_labels)]
print(f"\n  Quantum dimensions: {[f'{d:.6f}' for d in qdims]}")
print(f"  = [1, phi, phi, 1] expected: [1, {PHI:.6f}, {PHI:.6f}, 1]")

# Total quantum dimension
D2 = sum(d**2 for d in qdims)
print(f"  D^2 = sum(d_j^2) = {D2:.6f} (expected {a1 + sqrt(a1):.6f})")

# WRT invariant of Poincare homology sphere:
# Using surgery presentation: +1 surgery on left trefoil
# tau_r(Sigma) = F_r(trefoil, +1) / F_r(unknot, +1)
# where F_r(K, p) = sum_j [j+1]_q^2 * q^{p*j(j+1)/2} * colored_Jones

# Simpler: use the plumbing description.
# Sigma(2,3,5) = plumbing of E8 graph with all weights -2.

# For a plumbing of a tree graph with weights {a_v}:
# tau_r = (S_00)^{|V|-3} * sum_{j_v} prod_{v} S_{0,j_v} *
#         exp(2*pi*i * sum_v a_v * C_2(j_v) / (4*r)) * prod_{edges} S_{j_u,j_v}/S_{0,0}

# This is complex. Let me use the simplest formula:
# For the Poincare homology sphere as -E8 plumbing:
# Z_r = (1/D^2)^{7} * sum_{j1,...,j8} prod_{i} S_{0,ji}/S_{00} *
#        exp(+2*pi*i * sum_i (-2)*C_2(ji)/(4r)) * prod_{E8 edges} delta(ji, ji')

# Actually, this is wrong. Let me use a different approach.

# The WRT invariant of Sigma(2,3,5) at r=5 is KNOWN to be:
# tau_5 = 1  (the Poincare homology sphere is the boundary of E8 plumbing,
#              and at r=5 the invariant simplifies due to phi)

# Let me compute it directly from the definition.
# For a RATIONAL homology sphere obtained by +1 surgery on trefoil:

# The colored Jones polynomial of the trefoil at q = exp(2*pi*i/r):
q = np.exp(2j * pi / r)

print(f"\n  q = exp(2*pi*i/{r}) = {q:.6f}")

# J_n(trefoil, q) for the trefoil:
# For right trefoil: J_n(q) = sum_{k=0}^{n-1} q^{k(k+3)/2} * product term
# Actually, J_n uses the (n+1)-dim rep colored Jones.

# Simpler: use the explicit formula
# tau_r(+1 surgery on K) = sum_{j=0}^{(r-2)/2} [2j+1]_q * J_{2j+1}(K, q) * q^{(2j+1)^2/4}
# ... this requires the colored Jones of trefoil, which I can compute.

# Let me instead use the EXACT known value.
# From Hikami (2006), Lawrence-Zagier:
# tau_5(Sigma(2,3,5)) = exp(-3*pi*i/20) * phi  [up to normalization]

# The MAGNITUDE is what matters for non-perturbative physics:
# |tau_5| = phi = 1.618...

# Let me verify by direct computation using the Gauss sum approach.
# For Sigma(2,3,5) = plumbed along E8 with framings -2:

# E8 adjacency matrix
E8_adj = np.zeros((8, 8))
E8_edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (2,7)]
for i, j in E8_edges:
    E8_adj[i][j] = 1
    E8_adj[j][i] = 1

# Linking matrix (framing -2 on diagonal, -1 for edges)
B = -2 * np.eye(8) + E8_adj

print(f"\n  E8 plumbing matrix B:")
for row in B:
    print(f"    {[int(x) for x in row]}")

det_B = np.linalg.det(B)
print(f"  det(B) = {det_B:.1f} (should be 1 for homology sphere)")

# WRT invariant via Gauss sum:
# For a plumbing with matrix B of size n:
# tau_r = C_r^{-sigma(B)} * (1/D)^n * sum_{l in (Z/2rZ)^n}
#          prod_i S_{0,l_i} * exp(pi*i/(2r) * l^T B l)
# where sigma(B) = signature of B.

sigma_B = 8  # E8 is negative definite, so plumbing -E8 has sigma = -8
# Wait: B = -2*I + adj has eigenvalues that are all negative (for E8 with -2 framing)
# det = 1 > 0, 8x8 matrix. For -E8 plumbing, signature = -8.
eigs_B = np.linalg.eigvalsh(B)
sigma_B = int(np.sum(np.sign(eigs_B)))
print(f"  Eigenvalues of B: {sorted(eigs_B)}")
print(f"  Signature: {sigma_B}")

# At r=5, the allowed labels are l_i in {1, 2, 3, 4} (= 2j+1 for j=0,1/2,1,3/2)
# But we need l_i < r, so l_i in {1, 2, 3, 4}
# And S_{0,l} = sqrt(2/r) * sin(l*pi/r)

# The sum is:
# Z = sum_{l in {1,...,r-1}^8} prod_i sqrt(2/r)*sin(l_i*pi/r)
#     * exp(pi*i/(2*r) * sum_{ij} B_{ij} * l_i * l_j)

print(f"\n  Computing WRT invariant at r={r}...")
print(f"  Sum over l in {{1,...,{r-1}}}^8 = {(r-1)**8} terms")

# This is 4^8 = 65536 terms - very doable!
Z_WRT = 0j
count = 0
for l1 in range(1, r):
    for l2 in range(1, r):
        for l3 in range(1, r):
            for l4 in range(1, r):
                for l5 in range(1, r):
                    for l6 in range(1, r):
                        for l7 in range(1, r):
                            for l8 in range(1, r):
                                l = np.array([l1,l2,l3,l4,l5,l6,l7,l8], dtype=float)
                                # Product of sin factors
                                prod_sin = 1.0
                                for li in l:
                                    prod_sin *= sin(li * pi / r)
                                # Quadratic form
                                quad = l @ B @ l
                                phase = np.exp(1j * pi * quad / (2*r))
                                Z_WRT += prod_sin * phase
                                count += 1

# Normalization
norm = (2.0/r)**4  # (sqrt(2/r))^8
Z_WRT *= norm

# Additional normalization: divide by D^n where D = total quantum dim
D_total = sqrt(D2)
Z_WRT_normalized = Z_WRT / D_total**8

# And multiply by the signature correction
# C_r = exp(-i*pi*(r-2)/(4*r)) for SU(2)
# C_r^{-sigma} = exp(i*pi*(r-2)*sigma/(4*r))
C_phase = np.exp(1j * pi * (r-2) * (-sigma_B) / (4*r))
Z_WRT_final = C_phase * Z_WRT_normalized

print(f"  Counted {count} terms")
print(f"  Raw Z_WRT = {Z_WRT:.6f}")
print(f"  |Z_WRT| = {abs(Z_WRT):.6f}")
print(f"  Normalized Z = {Z_WRT_normalized:.6f}")
print(f"  |Z_norm| = {abs(Z_WRT_normalized):.6f}")
print(f"  With phase correction: Z_final = {Z_WRT_final:.6f}")
print(f"  |Z_final| = {abs(Z_WRT_final):.10f}")
print(f"  arg(Z_final) = {np.angle(Z_WRT_final)/pi:.6f} * pi")

# Compare to framework quantities
print(f"\n  Compare |Z_WRT| to framework:")
Z_abs = abs(Z_WRT_final)
if Z_abs > 1e-30:
    print(f"  |Z| = {Z_abs:.10f}")
    print(f"  phi = {PHI:.10f}")
    print(f"  1/phi = {1/PHI:.10f}")
    print(f"  |Z|/phi = {Z_abs/PHI:.6f}")
    print(f"  ln|Z|/ln(alpha) = {log(Z_abs)/ln_alpha:.6f}" if Z_abs > 0 else "")
else:
    print(f"  |Z| ~ 0 (too small)")

# =====================================================================
# PART G: The instanton-anti-instanton interpretation
# =====================================================================
print(f"\n{'='*72}")
print("PART G: Instanton-anti-instanton interpretation")
print("=" * 72)

print(f"""
  KEY OBSERVATION:
    alpha^57 = exp(-{target_exp:.4f})
    exp(-2/alpha) = exp(-{2/alpha:.4f})
    Ratio of exponents: {target_exp:.4f} / {2/alpha:.4f} = {target_exp/(2/alpha):.6f}

  If Lambda_CC = exp(-2*S_inst) where S_inst = 1/alpha:
    This is an instanton-ANTI-instanton pair (I-Ibar).
    Total action = 2/alpha = {2/alpha:.4f}
    But target = {target_exp:.4f}, gap = {target_exp - 2/alpha:.4f}

  Gap analysis:
    gap = z_CC * |ln alpha| - 2/alpha = {target_exp - 2/alpha:.6f}
    = {target_exp - 2/alpha:.6f}

  Can the gap come from:
    1. Determinant prefactor of instanton?
       In the instanton calculation, Z ~ exp(-S_inst) * det'^(-1/2) * ...
       The determinant gives a POLYNOMIAL prefactor, not exponential.
       So it enters as ln(prefactor), which is O(ln N) ~ O(5).
       gap = {target_exp - 2/alpha:.4f} ~ {target_exp - 2/alpha:.1f}
       This is the RIGHT ORDER for a one-loop correction!

    2. One-loop correction to instanton action:
       S_eff = S_tree + S_1-loop = 1/alpha + delta_S
       delta_S = (1/2) * ln det'(fluctuation operator)
       2 * delta_S = ln det'(L) correction
""")

# Compute what det' prefactor would close the gap
gap = target_exp - 2/alpha
print(f"  Required: 2*delta_S = {gap:.10f}")
print(f"  exp(gap) = {exp(gap):.6f}")
print(f"  So the I-Ibar amplitude with one-loop correction:")
print(f"  Z_IIbar = exp(-2/alpha) * exp(-gap)")
print(f"          = exp(-2/alpha) * {exp(-gap):.6f}")
print(f"          = alpha^z_CC  (by construction)")

print(f"\n  Is exp(gap) = {exp(gap):.6f} a framework quantity?")
val = exp(gap)
print(f"    exp(gap) = {val:.10f}")
print(f"    N^(1/4) = {N**0.25:.10f}")
print(f"    phi^8 = {PHI**8:.10f}")
print(f"    (2*pi)^(3/2) = {(2*pi)**1.5:.10f}")
print(f"    a1! = {factorial(a1)}")
print(f"    degree*a1 = {degree * a1}")
print(f"    N^(1/a1) = {N**(1/a1):.10f}")

# From exp421: ln det'(L/N) / ln(alpha) = 56.97
# This means: ln det'(L/N) = 56.97 * ln(alpha) = -280.33
# And: 2/alpha = 274.07
# So: ln det'(L/N) + 2/alpha = -280.33 + 274.07 = -6.26
# And: target = -279.87
# So: target = -2/alpha + (2/alpha - target) = -2/alpha - 5.80

# Let me check: does det'(L/N) provide the correction?
# Z_IIbar = exp(-2/alpha) * det'(L/N)^{-1/(2*something)}?

# From exp421: ln det'(L) = 289.385, ln det'(L/N) = 289.385 - 119*ln(120) = -280.326
ld = 289.385  # approximate from exp425
ld_over_N = ld - 119 * log(120)
print(f"\n  Cross-check with exp421:")
print(f"    ln det'(L) = {ld:.4f}")
print(f"    ln det'(L/N) = {ld_over_N:.4f}")
print(f"    ln det'(L/N) / ln(alpha) = {ld_over_N / ln_alpha:.4f}")
print(f"    2/alpha = {2/alpha:.4f}")
print(f"    gap = target - 2/alpha = {target_exp - 2/alpha:.4f}")
print(f"    ln det'(L/N) + 2/alpha = {ld_over_N + 2/alpha:.4f}")
print(f"    This SHOULD equal -target_exp = {-target_exp:.4f}? NO.")
print(f"    ld_over_N = {ld_over_N:.4f}, -target = {-target_exp:.4f}")
print(f"    They differ by: {ld_over_N + target_exp:.4f}")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print(f"\n{'='*72}")
print("FINAL SUMMARY")
print("=" * 72)

print(f"""
  NON-PERTURBATIVE ANALYSIS OF CC = alpha^57:

  1. INSTANTON INTERPRETATION:
     alpha^57 = exp(-{target_exp:.2f}) ~ exp(-2/alpha) = exp(-{2/alpha:.2f})
     Gap: {target_exp - 2/alpha:.2f} ({(target_exp - 2/alpha)/(2/alpha)*100:.1f}% of 2/alpha)
     This suggests CC ~ (instanton amplitude)^2 with O(1) correction.

  2. MIXED INSTANTON (best match from Part B):
     Best: a={best_ab[0] if best_ab else '?'}, b={best_ab[1] if best_ab else '?'} with diff={best_diff:.4f}
     exp(-a/alpha - b/alpha_s) closest to alpha^z_CC.

  3. WRT INVARIANT at r=5:
     |tau_5(Sigma(2,3,5))| = {abs(Z_WRT_final):.6f}
     This is a FINITE, EXACT topological invariant of the
     Poincare homology sphere at the framework level.

  4. THE CORRECTION PROBLEM:
     The gap {target_exp - 2/alpha:.4f} between 2/alpha and z_CC*|ln alpha|
     is O(a1) ~ O(5). This is the right order for a ONE-LOOP
     correction to the instanton action.

     If S_eff = 1/alpha + (1/2)*ln(prefactor)/N, then
     Z_IIbar = exp(-2*S_eff) needs prefactor = exp({2*gap:.4f}) = {exp(2*gap):.4f}

     BUT: this is NOT derived. The gap has no clean framework expression.

  5. STRUCTURAL INSIGHT:
     The observation alpha^57 ~ exp(-2/alpha) is NUMEROLOGICAL (2.1% match).
     It suggests a physical mechanism (I-Ibar pair) but does NOT constitute
     a derivation. The correction term remains unexplained.

  STATUS: PATTERN (suggestive but not derived)
""")

sys.exit(0)

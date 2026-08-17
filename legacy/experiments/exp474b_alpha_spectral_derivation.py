"""
exp474b_alpha_spectral_derivation.py
=====================================
SPECTRAL ACTION ROUTE TO THE LINEAR COEFFICIENT 4*a1*phi^4

GOAL: Find whether the linear coefficient in the alpha equation
      2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
can be DERIVED from the spectral action on the 600-cell, or
on the product geometry (600-cell x McKay graph).

CURRENT STATUS:
  - 2*pi: Hopf holonomy (DERIVED, exp293)
  - 1: Vieta/normalization (DERIVED)
  - 4*a1*phi^4 = 137.082...: "STRUCTURAL" (from icosahedral Laplacian
    eigenvalues L(3)*L(3')*phi^4, but phi^4 factor not from spectral action)

APPROACH:
  Systematically enumerate ALL spectral quantities and test combinations
  that could yield 4*a1*phi^4 = 20*phi^4 = 137.082...

Author: Razvan-Constantin Anghelina
Date: March 2026
"""

import numpy as np
from math import sqrt, pi, factorial, log, log10, exp, gcd
from fractions import Fraction
from itertools import combinations, product as iterproduct
from collections import defaultdict

# =============================================================================
# FRAMEWORK CONSTANTS (all from a1=5)
# =============================================================================
a1 = 5
b1 = a1 + 1                    # 6
phi = (1 + sqrt(a1)) / 2       # golden ratio = 1.618034...
phi_c = (1 - sqrt(a1)) / 2     # Galois conjugate = -0.618034...
N = factorial(a1)               # 120 = |2I|
degree = 2 * b1                 # 12 (vertex degree)
h_E8 = a1 * b1                 # 30 (Coxeter number E8)
rank_E8 = 8
N_gen = 3
N_eig = 9

# Couplings
sin2tW = b1 / (a1**2 + 1)      # 6/26
alpha_s = 1 / (2 * phi**3)
# Alpha from quadratic equation
disc_alpha = (4*a1*phi**4)**2 - 4*2*pi*1
alpha_em = (4*a1*phi**4 - sqrt(disc_alpha)) / (2*2*pi)
alpha_exp = 1 / 137.035999084

# TARGET
TARGET = 4 * a1 * phi**4       # 137.082039...

print("=" * 76)
print("EXP-474b: SPECTRAL ACTION ROUTE TO 4*a1*phi^4")
print("=" * 76)
print()
print(f"TARGET: 4*a1*phi^4 = {TARGET:.10f}")
print(f"  = 20*phi^4 = 20*(3phi+2) = 60*phi + 40")
print(f"  In Q(sqrt(5)): 70 + 30*sqrt(5)")
print(f"  1/alpha(theory) = {1/alpha_em:.6f}")
print(f"  1/alpha(CODATA) = 137.035999084")
print()

# =============================================================================
# SECTION 1: SPECTRAL DATA - SIMPLICIAL HODGE-DIRAC
# =============================================================================
print("=" * 76)
print("SECTION 1: SIMPLICIAL SPECTRAL DATA (Hodge-Dirac on 600-cell)")
print("=" * 76)
print()

# Seeley-DeWitt reduced coefficients
A0 = 2*a1 + 1                  # 11
A1 = 2*a1**2 + 2*a1 + 2        # 62
A2 = 6*a1**2 + 15*a1 + 8       # 233 = F_13

# Full Seeley-DeWitt coefficients
c0 = 2 * N * A0                # 2640
c1 = 2 * N * A1                # 14880
c2 = 2 * N * A2                # 55920

# f-vector of 600-cell
f0 = N                          # 120 vertices
f1 = N * b1                     # 720 edges
f2 = 2 * N * (a1 - 1) + N      # 1200 faces (triangles)
f3 = N * a1                     # 600 cells (tetrahedra)
# Actually: f0=120, f1=720, f2=1200, f3=600
f0, f1, f2, f3 = 120, 720, 1200, 600
hilbert_dim = f0 + f1 + f2 + f3  # 2640

print(f"f-vector: ({f0}, {f1}, {f2}, {f3})")
print(f"Hilbert space dim = {hilbert_dim} = 2*N*A0")
print()
print(f"Seeley-DeWitt reduced: A0={A0}, A1={A1}, A2={A2}")
print(f"  Full: c0={c0}, c1={c1}, c2={c2}")
print(f"  Identity: 2*A1^2+1 = {2*A1**2+1}, 3*A0*A2 = {3*A0*A2}")
print(f"  CHECK: {2*A1**2+1 == 3*A0*A2}")
print()

# Ratios between Seeley-DeWitt coefficients
r01 = c1 / c0                  # 14880/2640 = 56/10 = 28/5
r12 = c2 / c1                  # 55920/14880 = 233/62 * 120 ... no
r02 = c2 / c0                  # 55920/2640
r_sq = c1**2 / (c0 * c2)       # Should be ~A1^2/(A0*A2)

print("Spectral ratios:")
print(f"  c1/c0 = {r01:.6f} = {Fraction(c1,c0)}")
print(f"  c2/c1 = {r12:.6f} = {Fraction(c2,c1)}")
print(f"  c2/c0 = {r02:.6f} = {Fraction(c2,c0)}")
print(f"  c1^2/(c0*c2) = {r_sq:.10f} = {Fraction(c1**2, c0*c2)}")
print(f"  A1/A0 = {A1/A0:.6f} = {Fraction(A1,A0)}")
print(f"  A2/A1 = {A2/A1:.6f} = {Fraction(A2,A1)}")
print(f"  A2/A0 = {A2/A0:.6f} = {Fraction(A2,A0)}")
print(f"  A1^2/(A0*A2) = {A1**2/(A0*A2):.10f} = {Fraction(A1**2, A0*A2)}")
print()

# =============================================================================
# SECTION 2: ICOSAHEDRAL (VERTEX FIGURE) LAPLACIAN
# =============================================================================
print("=" * 76)
print("SECTION 2: ICOSAHEDRAL LAPLACIAN EIGENVALUES")
print("=" * 76)
print()

# Vertex figure = icosahedron (12 vertices, 30 edges, 20 faces)
# Adjacency eigenvalues: 5, sqrt(5), 1, -1, -sqrt(5), -5 (from icosahedron)
# Laplacian = degree - A, degree = 5 for icosahedron
# L_ico eigenvalues: L(1)=0, L(3)=5-sqrt(5), L(5)=6, L(3')=5+sqrt(5)
# with multiplicities 1, 3, 5, 3 (matching irreps of A5)

L1 = 0
L3 = a1 - sqrt(a1)             # 5 - sqrt(5) = 2.7639...
L5 = b1                        # 6
L3p = a1 + sqrt(a1)            # 5 + sqrt(5) = 7.2360...

print("Icosahedral Laplacian eigenvalues (12 vertices, degree 5):")
print(f"  L(1)  = {L1} (mult 1)")
print(f"  L(3)  = a1 - sqrt(a1) = {L3:.6f} (mult 3)")
print(f"  L(5)  = b1 = {L5} (mult 5)")
print(f"  L(3') = a1 + sqrt(a1) = {L3p:.6f} (mult 3)")
print()

# Key products and ratios
print("KEY IDENTITIES:")
print(f"  L(3)*L(3') = a1^2 - a1 = {a1*(a1-1)} = 4*a1 = {4*a1}")
print(f"    (uses a1*(a1-1) = 4*a1 iff a1=5)")
print(f"  L(3')/L(3) = {L3p/L3:.10f} = phi^2 = {phi**2:.10f}")
print(f"    CHECK: {abs(L3p/L3 - phi**2) < 1e-12}")
print(f"  L(3)+L(3') = 2*a1 = {2*a1}")
print(f"  L(3')^2/L(3) = {L3p**2/L3:.6f} (Galois-asymmetric)")
print(f"  L(3)*L(3')*phi^4 = 4*a1*phi^4 = {L3*L3p*phi**4:.10f} = TARGET")
print()

# =============================================================================
# SECTION 3: CAYLEY GRAPH (600-CELL 1-SKELETON)
# =============================================================================
print("=" * 76)
print("SECTION 3: CAYLEY GRAPH ADJACENCY EIGENVALUES")
print("=" * 76)
print()

# 600-cell 1-skeleton: 120 vertices, degree 12
# 9 distinct adjacency eigenvalues in Z[phi]
# Sorted descending:
cayley_eigs_zphi = [
    (12, 0),    # 12
    (0, 6),     # 6*phi
    (0, 4),     # 4*phi
    (3, 0),     # 3
    (0, 0),     # 0
    (-2, 0),    # -2
    (-3, 0),    # -3
    (4, -4),    # 4 - 4*phi = 4*(1-phi) = -4*phi_c... wait
    (6, -6),    # 6 - 6*phi
]
# Correction: let me use the proper list from the context
cayley_eigs_vals = [12, 6*phi, 4*phi, 3, 0, -2, -3, 4-4*phi, 6-6*phi]
cayley_mults = [1, 4, 9, 16, 25, 25, 16, 9, 4]
# Wait, the problem statement says: {12, 6phi, 4phi, 3, 0, -2, -3, 3-4phi, 6-6phi}
# with mults {1,4,9,16,25,25,16,9,4}
# Let me recheck: 4-4phi or 3-4phi?
# From exp117b: 12, 6phi, 4phi, 3, 0, -2, 4-4phi, -3, 6-6phi
# But context says: {12, 6phi, 4phi, 3, 0, -2, -3, 3-4phi, 6-6phi}
# There might be a sorting discrepancy. Let me compute both and verify sum.

# Using the set from the context:
cayley_eigs = {
    12:         1,
    6*phi:      4,
    4*phi:      9,
    3:          16,
    0:          25,
    -2:         25,
    -3:         16,
    3-4*phi:    9,   # = 3 - 4*1.618 = 3 - 6.472 = -3.472
    6-6*phi:    4,   # = 6 - 9.708 = -3.708
}

print("Cayley graph adjacency eigenvalues:")
total_mult = 0
total_trace = 0
for eig, mult in sorted(cayley_eigs.items(), reverse=True):
    print(f"  lambda = {eig:+10.6f}  (mult {mult:2d})")
    total_mult += mult
    total_trace += eig * mult
print(f"  Total multiplicity: {total_mult} (should be {N})")
print(f"  Tr(A) = {total_trace:.6f} (should be 0)")
print()

# Laplacian eigenvalues: mu_i = degree - lambda_i
cayley_lap_eigs = {}
for eig, mult in cayley_eigs.items():
    mu = degree - eig
    cayley_lap_eigs[mu] = mult

print("Cayley Laplacian eigenvalues (mu = 12 - lambda):")
for mu, mult in sorted(cayley_lap_eigs.items()):
    print(f"  mu = {mu:+10.6f}  (mult {mult:2d})")
print()

# Spectral gap of Cayley Laplacian
cayley_lap_sorted = sorted(cayley_lap_eigs.keys())
cayley_gap = cayley_lap_sorted[1] if cayley_lap_sorted[0] == 0 else cayley_lap_sorted[0]
print(f"Cayley Laplacian spectral gap: {cayley_gap:.6f}")
print(f"  = 12 - 6*phi = {12 - 6*phi:.6f}")
print(f"  = 6*(2-phi) = 6/phi^2 = {6/phi**2:.6f}")
print()

# =============================================================================
# SECTION 4: COEXACT LAPLACIAN
# =============================================================================
print("=" * 76)
print("SECTION 4: COEXACT LAPLACIAN SPECTRAL GAP")
print("=" * 76)
print()

# The coexact part of Delta_1 has spectral gap
lam1_coexact = 7 - 4*phi       # = 0.5278...
N_lam1 = 7**2 - 7*4 - 4**2     # wait, Galois norm of 7-4phi
# N(a+b*phi) = a^2 + a*b - b^2 (for Z[phi] element a+b*phi = a+b*(1+sqrt5)/2)
# Actually for element a - b*phi where phi=(1+sqrt5)/2:
# 7 - 4*phi = 7 - 4*(1+sqrt5)/2 = 7 - 2 - 2*sqrt5 = 5 - 2*sqrt5
# Norm = (5-2*sqrt5)(5+2*sqrt5) = 25-20 = 5 = a1
# Yes: N(lam1_coexact) = a1

print(f"Coexact spectral gap: lambda_1 = 7 - 4*phi = {lam1_coexact:.10f}")
print(f"  Galois norm: N(7-4*phi) = (7-4*phi)*(7-4*phi') = {lam1_coexact*(7-4*phi_c):.6f} = a1")
print(f"  Galois conjugate: 7-4*phi' = {7-4*phi_c:.6f}")
print()

# Hopf fiber spectral gap
hopf_gap = 1 / phi**2
print(f"Hopf fiber (C_10) spectral gap: 1/phi^2 = {hopf_gap:.10f}")
print(f"  = 2 - 2*cos(2*pi/10) = 2*(1-cos(pi/5))")
print()

# Lambda_max
Lambda_max_D2 = b1 * phi**2    # 15.7082...
print(f"Lambda_max(D^2) = b1*phi^2 = {Lambda_max_D2:.6f}")
print()

# =============================================================================
# SECTION 5: McKAY GRAPH (FINITE DIRAC)
# =============================================================================
print("=" * 76)
print("SECTION 5: McKAY GRAPH / FINITE DIRAC")
print("=" * 76)
print()

# McKay graph of 2I has 9 nodes (irreps of 2I)
# D_F eigenvalues give Tr(D_F^2) = 8 = rank(E8)
Tr_DF2 = rank_E8                # 8
dim_F = N_eig                   # 9

print(f"McKay graph: {dim_F} nodes (irreps of 2I)")
print(f"  Tr(D_F^2) = {Tr_DF2} = rank(E8)")
print(f"  dim(F) = {dim_F} = N_eig")
print()

# =============================================================================
# SECTION 6: SYSTEMATIC SEARCH FOR EXACT MATCHES
# =============================================================================
print("=" * 76)
print("SECTION 6: SYSTEMATIC SEARCH FOR 4*a1*phi^4")
print("=" * 76)
print()

def check_match(name, value, category="?"):
    """Check if value matches TARGET = 4*a1*phi^4"""
    ratio = value / TARGET if TARGET != 0 else float('inf')
    err_pct = abs(ratio - 1) * 100
    match = err_pct < 1e-8
    tag = "EXACT" if match else f"off by {err_pct:.4f}%"
    symbol = "***" if match else "   "
    print(f"  {symbol} {name:60s} = {value:14.6f}  (ratio={ratio:.10f})  [{tag}]  [{category}]")
    return match

found_exact = []

print("--- 6A: Simple ratios of Seeley-DeWitt coefficients ---")
print()

matches = []

# Direct ratios
m = check_match("c1/c0 = A1/A0 * 2N/2N", c1/c0, "TRIVIAL")
m = check_match("c2/c1", c2/c1, "TRIVIAL")
m = check_match("c2/c0", c2/c0, "TRIVIAL")

# Scaled by small integers and powers of phi
for p in range(-4, 5):
    for q_num in range(1, 30):
        for q_den in range(1, 20):
            val = (c1/c0) * phi**p * q_num / q_den
            if abs(val/TARGET - 1) < 1e-10:
                found_exact.append(f"(c1/c0)*phi^{p}*{q_num}/{q_den}")
                check_match(f"(c1/c0)*phi^{p}*{q_num}/{q_den}", val, "SEARCH")
            val = (c2/c1) * phi**p * q_num / q_den
            if abs(val/TARGET - 1) < 1e-10:
                found_exact.append(f"(c2/c1)*phi^{p}*{q_num}/{q_den}")
                check_match(f"(c2/c1)*phi^{p}*{q_num}/{q_den}", val, "SEARCH")
            val = (c2/c0) * phi**p * q_num / q_den
            if abs(val/TARGET - 1) < 1e-10:
                found_exact.append(f"(c2/c0)*phi^{p}*{q_num}/{q_den}")
                check_match(f"(c2/c0)*phi^{p}*{q_num}/{q_den}", val, "SEARCH")

if not found_exact:
    print("  No exact matches from simple Seeley-DeWitt ratios * phi^p * (small rational)")
print()

# Now try with reduced coefficients
print("--- 6B: Reduced coefficients A0, A1, A2 with eigenvalue expressions ---")
print()

found_6b = []

# Build a list of candidate "eigenvalue expressions"
eig_exprs = {
    "L(3)":         L3,
    "L(3')":        L3p,
    "L(5)=b1":      L5,
    "4*a1":         4*a1,
    "phi":          phi,
    "phi^2":        phi**2,
    "phi^3":        phi**3,
    "phi^4":        phi**4,
    "phi^5":        phi**5,
    "1/phi":        1/phi,
    "1/phi^2":      1/phi**2,
    "sqrt(a1)":     sqrt(a1),
    "a1":           float(a1),
    "b1":           float(b1),
    "N":            float(N),
    "degree":       float(degree),
    "h_E8":         float(h_E8),
    "rank_E8":      float(rank_E8),
    "N_gen":        float(N_gen),
    "N_eig":        float(N_eig),
    "Tr_DF2":       float(Tr_DF2),
    "lam1_coexact": lam1_coexact,
    "hopf_gap":     hopf_gap,
    "Lambda_max":   Lambda_max_D2,
    "2*pi":         2*pi,
    "pi":           pi,
    "pi/a1":        pi/a1,
}

reduced_coeffs = {
    "A0": A0,
    "A1": A1,
    "A2": A2,
    "A1/A0": A1/A0,
    "A2/A1": A2/A1,
    "A2/A0": A2/A0,
}

for rc_name, rc_val in reduced_coeffs.items():
    for ee_name, ee_val in eig_exprs.items():
        if ee_val == 0:
            continue
        # rc * ee
        val = rc_val * ee_val
        if abs(val/TARGET - 1) < 1e-8:
            name = f"{rc_name} * {ee_name}"
            check_match(name, val, "CHECK")
            found_6b.append(name)
        # rc / ee
        val = rc_val / ee_val
        if abs(val/TARGET - 1) < 1e-8:
            name = f"{rc_name} / {ee_name}"
            check_match(name, val, "CHECK")
            found_6b.append(name)

# Two eigenvalue expressions
for ee1_name, ee1_val in eig_exprs.items():
    for ee2_name, ee2_val in eig_exprs.items():
        if ee2_val == 0:
            continue
        for rc_name, rc_val in reduced_coeffs.items():
            val = rc_val * ee1_val * ee2_val
            if abs(val/TARGET - 1) < 1e-8:
                name = f"{rc_name} * {ee1_name} * {ee2_name}"
                check_match(name, val, "SEARCH")
                found_6b.append(name)
            val = rc_val * ee1_val / ee2_val
            if abs(val/TARGET - 1) < 1e-8:
                name = f"{rc_name} * {ee1_name} / {ee2_name}"
                check_match(name, val, "SEARCH")
                found_6b.append(name)

if not found_6b:
    print("  No exact matches from A_k * eigenvalue expressions")
print()

# =============================================================================
# SECTION 7: PRODUCT GEOMETRY M x F
# =============================================================================
print("=" * 76)
print("SECTION 7: PRODUCT GEOMETRY (600-cell x McKay)")
print("=" * 76)
print()

# In product geometry: D_total = D_M tensor 1_F + gamma tensor D_F
# D_total^2 = D_M^2 tensor 1_F + 1_M tensor D_F^2  (cross terms vanish)
#
# Seeley-DeWitt decomposition:
#   c0_total = c0_M * dim_F
#   c1_total = c1_M * dim_F + dim_M * Tr(D_F^2)
#   c2_total = c2_M * dim_F + c1_M * Tr(D_F^2) + dim_M * Tr(D_F^4) / 2

c0_M = c0
c1_M = c1
c2_M = c2  # This is (1/2)*Tr(D_M^4)

# McKay finite Dirac
dim_McKay = dim_F               # 9
Tr_DF2_val = Tr_DF2             # 8
# Tr(D_F^4) - need to estimate. From the McKay graph structure:
# D_F has eigenvalues that give Tr(D_F^2)=8. For binary icosahedral:
# The eigenvalues of D_F are typically the "mass exponents" n_f
# From the paper: Tr(D_F^2) = 8 = rank(E8)
# We need Tr(D_F^4). The McKay graph is the affine E8 Dynkin diagram.
# Eigenvalues of the adjacency matrix of affine E8 (10 nodes):
# For the extended Dynkin diagram (9+1 nodes), the adjacency gives specific values.
# But D_F is specifically the McKay finite Dirac, which is built differently.
# Let me use what we know: D_F eigenvalues sum of squares = 8.
# A common set for binary icosahedral McKay: dimensions of irreps give the McKay graph
# The irreps have dimensions: 1, 2, 3, 4, 5, 6, 4, 2, 3
# (the affine E8 Dynkin diagram has 9 nodes for 2I)
# D_F is the finite Dirac operator whose square gives mass eigenvalues.
# In Connes' framework, D_F encodes Yukawa couplings.
# From the memory: ||Y||_F = 1/sqrt(2), Tr(D_F^2) = 8
# We'll treat Tr(D_F^4) as unknown for now and see what constraint alpha gives.

print("Product geometry M x F:")
print(f"  M = 600-cell simplicial complex, dim_M = {c0_M}")
print(f"  F = McKay graph, dim_F = {dim_McKay}")
print(f"  Tr(D_F^2) = {Tr_DF2_val} = rank(E8)")
print()

c0_total = c0_M * dim_McKay
c1_total = c1_M * dim_McKay + c0_M * Tr_DF2_val
print(f"Product Seeley-DeWitt:")
print(f"  c0_total = c0*dim_F = {c0_M}*{dim_McKay} = {c0_total}")
print(f"  c1_total = c1*dim_F + c0*Tr(D_F^2) = {c1_M}*{dim_McKay} + {c0_M}*{Tr_DF2_val}")
print(f"           = {c1_M*dim_McKay} + {c0_M*Tr_DF2_val} = {c1_total}")
print()

# In Connes' spectral action, the gauge coupling comes from:
# 1/g^2 ~ f_0 * a_4 / (2*Lambda^4)
# where a_4 involves Tr(F^2) terms
# The KEY relation is: at unification, all couplings are equal to g^2/(2*f_0)
# In our framework: alpha = e^2/(4*pi) relates to the spectral action through
# specific combinations of the coefficients.

# Let's see what ratios of product coefficients give TARGET
print("Testing ratios involving product geometry:")
print()

prod_quantities = {
    "c0_total":     c0_total,
    "c1_total":     c1_total,
    "c0_M":         c0_M,
    "c1_M":         c1_M,
    "c2_M":         c2_M,
    "dim_F":        float(dim_McKay),
    "Tr_DF2":       float(Tr_DF2_val),
    "c1_M/c0_M":    c1_M / c0_M,
    "c2_M/c1_M":    c2_M / c1_M,
    "c1_total/c0_total": c1_total / c0_total,
    "gravity_part": float(c1_M * dim_McKay),    # gravity piece of c1_total
    "matter_part":  float(c0_M * Tr_DF2_val),   # matter piece of c1_total
}

found_7 = []
for pq_name, pq_val in prod_quantities.items():
    for ee_name, ee_val in eig_exprs.items():
        if ee_val == 0:
            continue
        val = pq_val * ee_val
        if abs(val/TARGET - 1) < 1e-8:
            name = f"{pq_name} * {ee_name}"
            check_match(name, val, "PRODUCT")
            found_7.append(name)
        val = pq_val / ee_val
        if abs(val/TARGET - 1) < 1e-8:
            name = f"{pq_name} / {ee_name}"
            check_match(name, val, "PRODUCT")
            found_7.append(name)

if not found_7:
    print("  No exact matches from product geometry single-step.")
print()

# =============================================================================
# SECTION 8: THE KEY ROUTE - EIGENVALUE RATIOS
# =============================================================================
print("=" * 76)
print("SECTION 8: EIGENVALUE-BASED ROUTES")
print("=" * 76)
print()

# Route 1: From the Cayley graph perspective
# lambda_4/lambda_1 = 12 / (6/phi^2) = 2*phi^2  (exp073)
# 5 * (2*phi^2)^2 = 5 * 4*phi^4 = 20*phi^4 = TARGET

val_route1 = a1 * (degree / (6/phi**2))**2
print(f"ROUTE 1: a1 * (degree / gap_Cayley)^2")
print(f"  = {a1} * ({degree} / {6/phi**2:.6f})^2")
print(f"  = {a1} * ({degree * phi**2/6:.6f})^2")
print(f"  = {a1} * (2*phi^2)^2 = {a1} * 4*phi^4")
print(f"  = {val_route1:.10f}")
check_match("a1 * (degree/gap)^2", val_route1, "STRUCTURAL")
print()

# Route 2: From icosahedral Laplacian
# L(3') * L(3) * phi^4 = (a1+sqrt(a1))*(a1-sqrt(a1)) * phi^4 = 4a1*phi^4
val_route2 = L3 * L3p * phi**4
print(f"ROUTE 2: L(3)*L(3')*phi^4 = {val_route2:.10f}")
check_match("L(3)*L(3')*phi^4", val_route2, "STRUCTURAL")
print()

# Route 3: L(3')^3 / L(3)
# L(3') = a1+sqrt(a1), L(3) = a1-sqrt(a1)
# L(3')^3 = (a1+sqrt(a1))^3
# L(3')^3/L(3) = (a1+sqrt(a1))^3/(a1-sqrt(a1))
#              = (a1+sqrt(a1))^3*(a1+sqrt(a1)) / [(a1-sqrt(a1))*(a1+sqrt(a1))]
#              = (a1+sqrt(a1))^4 / (a1^2-a1)
#              = (a1+sqrt(a1))^4 / (4*a1)
# Hmm, that doesn't simplify nicely. Let me just compute:
val_route3 = L3p**3 / L3
print(f"ROUTE 3: L(3')^3/L(3)")
print(f"  = {L3p:.6f}^3 / {L3:.6f}")
print(f"  = {L3p**3:.6f} / {L3:.6f}")
print(f"  = {val_route3:.10f}")
check_match("L(3')^3/L(3)", val_route3, "CHECK")
# Verify: L(3')^3/L(3) = L(3')*L(3')*(L(3')/L(3)) = L(3')*L(3')*phi^2
#        = L(3')^2*phi^2 = (a1+sqrt(a1))^2 * phi^2
# But TARGET = L(3)*L(3')*phi^4 = 4*a1*phi^4
# L(3')^3/L(3) = L(3')^2 * L(3')/L(3) = L(3')^2 * phi^2
# L(3')^2 = (a1+sqrt(a1))^2 = a1^2+2*a1*sqrt(a1)+a1 = a1^2+a1+2*a1*sqrt(a1)
#         = a1*(a1+1)+2*a1*sqrt(a1) = a1*b1 + 2*a1*sqrt(a1) = 30 + 10*sqrt(5)
# L(3')^2 * phi^2 = (30+10*sqrt(5)) * (3+sqrt(5))/2 = ...
# Let me compute numerically:
print(f"  L(3')^2 = {L3p**2:.6f}")
print(f"  L(3')^2 * phi^2 = {L3p**2 * phi**2:.10f}")
print(f"  TARGET = {TARGET:.10f}")
# Check algebraically: L(3')^2*phi^2 vs 4*a1*phi^4
# 4*a1*phi^4 = 4*5*(3+sqrt5)^2/4 = 5*(3+sqrt5)^2 = 5*(9+6sqrt5+5) = 5*14+30*sqrt5 = 70+30*sqrt5
# L(3')^2*phi^2 = (5+sqrt5)^2*(3+sqrt5)/2 = (30+10sqrt5)*(3+sqrt5)/2
#   = (90+30sqrt5+30sqrt5+50)/2 = (140+60sqrt5)/2 = 70+30sqrt5
# YES! Matches exactly.
print(f"  Algebraic check: both = 70 + 30*sqrt(5) = {70+30*sqrt(5):.10f}")
print(f"  MATCH: {abs(val_route3 - TARGET) < 1e-8}")
print()

# Route 4: From context hint: L(3)*L(3')*(L(3')/L(3))^2
val_route4 = L3 * L3p * (L3p/L3)**2
print(f"ROUTE 4: L(3)*L(3')*(L(3')/L(3))^2")
print(f"  = {L3*L3p:.6f} * {(L3p/L3)**2:.6f}")
print(f"  = 4*a1 * phi^4")
print(f"  = {val_route4:.10f}")
check_match("L(3)*L(3')*(L(3')/L(3))^2", val_route4, "STRUCTURAL")
print(f"  This is just the identity 4*a1*phi^4 rewritten using L eigenvalues.")
print()

# =============================================================================
# SECTION 9: THE SPECTRAL ACTION PERSPECTIVE
# =============================================================================
print("=" * 76)
print("SECTION 9: SPECTRAL ACTION -- CAN IT DETERMINE 4*a1*phi^4?")
print("=" * 76)
print()

# In Connes-Chamseddine spectral action:
# S = Tr(f(D/Lambda)) ~ f_4*Lambda^4*a_0 + f_2*Lambda^2*a_2 + f_0*a_4 + ...
#
# The gauge coupling comes from the a_4 term:
# a_4 = (1/16*pi^2) * int [... alpha_s * Tr(G^2) + alpha * Tr(F^2) + ...]
#
# For the PRODUCT geometry M x F:
# a_4 gets contributions from: R^2 terms (gravity), F^2 terms (gauge), Yukawa terms
#
# The coupling unification condition at the spectral cutoff Lambda is:
#   g_1^2 = g_2^2 = g_3^2 = pi^2 / (2*f_0*N_gen)  [Connes normalization]
# or equivalently: alpha = g^2/(4*pi) = pi/(8*f_0*N_gen)
#
# But this gives alpha AT the cutoff, not at low energy.
# The RG running from Lambda down to m_Z determines the physical alpha.
#
# KEY QUESTION: Does the spectral action ON THE 600-CELL directly give
# 4*a1*phi^4 without needing RG running?

print("In Connes' framework:")
print("  alpha_GUT = pi / (8 * f_0 * N_gen)")
print(f"  If f_0 = N^k/(2*pi) (our framework): at k=0, f_0 = 1/(2*pi)")
print(f"    => alpha_GUT = pi / (8 * (1/(2*pi)) * 3) = pi^2/(12) = {pi**2/12:.6f}")
print(f"    => 1/alpha_GUT = {12/pi**2:.6f}")
print(f"  This is NOT 137.08. (Connes gives unification coupling, not alpha_EM)")
print()

# Instead, try the DISCRETE spectral action on the 600-cell
# The 600-cell is NOT a smooth manifold -- it's a simplicial complex.
# The spectral action coefficients c_k are EXACT integers (up to factors of 240).

print("Discrete spectral action coefficients:")
print(f"  c0 = 2*N*A0 = {c0} = dim(H)")
print(f"  c1 = 2*N*A1 = {c1} = Tr(D^2)")
print(f"  c2 = 2*N*A2 = {c2} = (1/2)*Tr(D^4)")
print()

# Can we extract 4*a1*phi^4 from the RATIO of spectral coefficients?
print("CRITICAL TEST: Does c2/c1 * (something purely from a1) = 4*a1*phi^4?")
print()

# c2/c1 = A2/A1 = 233/62
ratio_21 = Fraction(A2, A1)
print(f"  c2/c1 = A2/A1 = {ratio_21} = {float(ratio_21):.10f}")
print(f"  TARGET / (A2/A1) = {TARGET / float(ratio_21):.10f}")
# TARGET / (233/62) = 137.082... * 62/233 = 36.4607...
needed = TARGET * A1 / A2
print(f"    Need to multiply by: {needed:.10f}")

# Check if this number has a nice form
for p in range(-6, 7):
    val = phi**p
    ratio_check = needed / val
    if abs(ratio_check - round(ratio_check)) < 0.001:
        k = round(ratio_check)
        print(f"    ~ {k} * phi^{p} = {k*val:.6f} (err = {abs(needed-k*val)/needed*100:.4f}%)")

# Try: needed ~ a1*phi^p
for p in range(-6, 7):
    ratio_check = needed / (a1 * phi**p)
    if abs(ratio_check - round(ratio_check)) < 0.01:
        k = round(ratio_check)
        result_val = float(ratio_21) * k * a1 * phi**p
        err = abs(result_val/TARGET - 1) * 100
        print(f"    ~ {k}*a1*phi^{p} = {k*a1*phi**p:.6f} -> c2/c1*{k}*a1*phi^{p} = {result_val:.6f} ({err:.6f}%)")

print()

# Alternative: c1/c0 * something
ratio_10 = Fraction(c1, c0)  # 14880/2640 = 28/5
print(f"  c1/c0 = {ratio_10} = {float(ratio_10):.10f}")
needed_10 = TARGET / float(ratio_10)
print(f"  TARGET / (c1/c0) = {needed_10:.10f}")
# 137.082 / 5.6 = 24.4789...
for p in range(-6, 7):
    val = phi**p
    ratio_check = needed_10 / val
    if abs(ratio_check - round(ratio_check)) < 0.01:
        k = round(ratio_check)
        result_val = float(ratio_10) * k * phi**p
        err = abs(result_val/TARGET - 1) * 100
        if err < 1:
            print(f"    ~ {k}*phi^{p}: c1/c0*{k}*phi^{p} = {result_val:.6f} ({err:.6f}%)")

print()

# =============================================================================
# SECTION 10: DEEP ANALYSIS -- THE ROLE OF phi^4
# =============================================================================
print("=" * 76)
print("SECTION 10: WHERE DOES phi^4 COME FROM?")
print("=" * 76)
print()

# The problem: 4*a1 = L(3)*L(3') is DERIVED (from icosahedral Laplacian)
# But phi^4 is the "spectral amplification" factor. Where does it arise?
#
# Observation 1: The Cayley adjacency eigenvalue ratio
# lambda_max / gap = 12 / (12-6*phi) = 12 / (6/phi^2) = 2*phi^2
# Then (lambda_max/gap)^2 = 4*phi^4
# So: TARGET = a1 * (lambda_max/gap)^2

val_10a = a1 * (degree / cayley_gap)**2
print(f"ROUTE 5: a1 * (degree/Cayley_gap)^2")
print(f"  = a1 * (12 / (12-6*phi))^2")
print(f"  = a1 * (2*phi^2)^2")
print(f"  = a1 * 4*phi^4 = {val_10a:.10f}")
check_match("a1 * (degree/Cayley_gap)^2", val_10a, "STRUCTURAL")
print()

# Observation 2: Cheeger ratio
# In spectral geometry, the Cheeger constant h relates to the spectral gap:
# lambda_1 >= h^2/2  (Cheeger inequality)
# The spectral gap of the Cayley Laplacian is 12 - 6*phi = 6/phi^2
# The expansion ratio: degree / gap = 12*phi^2/6 = 2*phi^2
# This measures "how quickly the heat spreads on the 600-cell graph"

print("OBSERVATION: The spectral expansion ratio of the Cayley graph")
print(f"  R = degree / gap = {degree} / {cayley_gap:.6f} = {degree/cayley_gap:.10f}")
print(f"  = 2*phi^2 = {2*phi**2:.10f}")
print(f"  R^2 = 4*phi^4 = {4*phi**4:.10f}")
print(f"  TARGET = a1 * R^2 = {a1 * (degree/cayley_gap)**2:.10f}")
print()
print("  Physical interpretation: R = degree/gap measures the spectral")
print("  spreading efficiency of the Cayley graph. The coupling is determined")
print("  by a1 times the SQUARE of this spreading ratio.")
print("  CATEGORY: STRUCTURAL (R from spectral data, but why a1*R^2?)")
print()

# Observation 3: From the coexact Laplacian perspective
# The coexact spectral gap is lambda_1_coex = 7 - 4*phi
# N(lambda_1_coex) = a1 = 5
# What about: (degree / lambda_1_coex)^2 ?
val_10c = (degree / lam1_coexact)**2
print(f"  (degree / lam1_coexact)^2 = ({degree}/{lam1_coexact:.6f})^2 = {val_10c:.6f}")
print(f"  TARGET = {TARGET:.6f}")
print(f"  Ratio: {val_10c/TARGET:.6f} (not a match)")
print()

# =============================================================================
# SECTION 11: KALUZA-KLEIN PERSPECTIVE
# =============================================================================
print("=" * 76)
print("SECTION 11: KALUZA-KLEIN ON HOPF FIBRATION")
print("=" * 76)
print()

# S^3 -> S^2 with U(1) fiber
# In KK theory: 1/g^2 ~ Vol(fiber) / G_higher
# The Hopf fiber in the 600-cell is C_10 (decagonal cycle, 10 vertices)
# Volume of fiber (as fraction of total) = 10/120 = 1/12
# But the coupling relates to the fiber's eigenvalue structure

hopf_length = 2 * a1             # 10 vertices in C_10
vol_fiber = hopf_length / N       # 10/120 = 1/12

print(f"Hopf fiber: C_10 cycle")
print(f"  Length = 2*a1 = {hopf_length} vertices")
print(f"  Vol(fiber)/Vol(total) = {hopf_length}/{N} = {vol_fiber:.6f} = 1/{N//hopf_length}")
print()

# KK coupling: in standard KK on S^1 of radius R,
# 1/g^2 = 2*pi*R / G_5D
# The coupling "runs" with energy scale. The 4D coupling at the KK scale is:
# alpha_KK = g^2/(4*pi) = G_5D / (8*pi^2*R)
# In our discrete setting: R -> 1/gap(C_10)

# Spectral gap of C_10: 2*(1-cos(2*pi/10)) = 2*(1-cos(pi/5)) = 2*(1-(1+sqrt5)/4)
# Actually: eigenvalues of C_n Laplacian are 2*(1-cos(2*pi*k/n)) for k=0,...,n-1
# For C_10: gap = 2*(1-cos(2*pi/10)) = 2*(1-cos(pi/5))
cos_pi5 = (1+sqrt(5))/4  # cos(pi/5) = phi/2
gap_C10 = 2*(1-cos_pi5)
print(f"Spectral gap of C_10:")
print(f"  gap = 2*(1-cos(pi/5)) = 2*(1-phi/2) = 2-phi = {gap_C10:.10f}")
print(f"  = 1/phi^2 = {1/phi**2:.10f}")
print(f"  CHECK: {abs(gap_C10 - 1/phi**2) < 1e-12}")
print()

# Number of fibers: N/hopf_length = 120/10 = 12
n_fibers = N // hopf_length     # 12
print(f"Number of parallel fibers: {n_fibers} = degree = 2*b1")
print()

# KK-type combination:
# 1/alpha ~ (N_fibers) * (1/gap)^2 * a1
# = 12 * phi^4 * 5 = 60*phi^4 ...hmm that's 3*TARGET
val_kk1 = n_fibers * (1/gap_C10)**2 * a1
print(f"  N_fibers * (1/gap)^2 * a1 = {n_fibers} * phi^4 * {a1} = {val_kk1:.6f}")
print(f"  Ratio to TARGET: {val_kk1/TARGET:.6f}")

# Try: (N/hopf_length) * (1/gap)^2 * (L(3)*L(3')/degree)
val_kk2 = n_fibers * (1/gap_C10)**2 * L3 * L3p / degree
print(f"  N_fibers * (1/gap)^2 * L(3)*L(3')/degree = {val_kk2:.6f}")
print(f"  Ratio to TARGET: {val_kk2/TARGET:.6f}")

# The simplest KK relation:
# 1/alpha_0 = (fiber_length/2*pi) * (1/gap) * degree * ???
# Let me think more carefully.

# In standard KK: gauge coupling ~ 1/(R*Lambda_KK) where R=fiber radius, Lambda_KK=compactification scale
# For discrete: R_eff ~ 1/gap(fiber) = phi^2
# Lambda_KK ~ degree (or something spectral)
# So 1/alpha ~ R_eff^2 * (something)

# The exact KK identity from our framework:
# alpha * alpha' = 1/(2*pi)
# where alpha' is the Galois conjugate
# This is the KK PRODUCT relation.
# Combined with the quadratic equation, this FIXES alpha.
# But the quadratic equation already contains 4*a1*phi^4.

print()
print("KK product relation: alpha * alpha' = 1/(2*pi)")
print(f"  alpha' = {1/(2*pi*alpha_em):.10f}")
print(f"  Product: {alpha_em * (1/(2*pi*alpha_em)):.10f} = 1/(2*pi) = {1/(2*pi):.10f}")
print()

# Vieta's formulas:
# alpha + alpha' = 4*a1*phi^4 / (2*pi) = TARGET / (2*pi)
# alpha * alpha' = 1 / (2*pi)
# So: TARGET = 2*pi * (alpha + alpha') ... but this is circular.
# Unless we can determine (alpha + alpha') independently.

print("Vieta's formulas for the alpha equation:")
print(f"  alpha + alpha' = 4*a1*phi^4 / (2*pi) = {TARGET/(2*pi):.10f}")
print(f"  alpha * alpha' = 1/(2*pi) = {1/(2*pi):.10f}")
print(f"  This means: 4*a1*phi^4 = 2*pi*(alpha+alpha')")
print(f"  CIRCULAR unless we determine the sum independently.")
print()

# =============================================================================
# SECTION 12: HEAT KERNEL AND SPECTRAL ZETA APPROACH
# =============================================================================
print("=" * 76)
print("SECTION 12: HEAT KERNEL / SPECTRAL ZETA")
print("=" * 76)
print()

# On a manifold M, the spectral action Tr(f(D^2/Lambda^2)) can be computed
# via the heat kernel: K(t) = Tr(exp(-t*D^2)) = sum_n exp(-t*lambda_n)
# The small-t expansion: K(t) ~ a_0*t^{-d/2} + a_1*t^{-(d-2)/2} + a_2 + ...
# For the 600-cell (d=3, simplicial):
# K(t) ~ c_0*t^0 - c_1*t + c_2*t^2/2 - ... (discrete version)

# The spectral zeta function:
# zeta_D(s) = Tr(D^{-2s}) = sum_{lambda>0} lambda^{-s}
# zeta_D(0) = number of nonzero eigenvalues - dim(ker)

# What we can compute from c0, c1, c2:
# The "effective dimension" from the spectral action:
# d_eff = 2 * c1^2 / (c0*c2 - c1^2/3)
# Actually, for the heat kernel K(t) = c0 - c1*t + c2*t^2 + ...
# The "spectral dimension" at scale t is d_s(t) = -2*t*K'(t)/K(t)
# At small t: d_s ~ 2*c1*t/c0 (to leading order)

# Let's use a different approach: the RATIO c1^2/(c0*c2) encodes the spectral geometry
ratio_sq = c1**2 / (c0 * c2)
print(f"Spectral ratio: c1^2/(c0*c2) = {ratio_sq:.10f}")
print(f"  = A1^2/(A0*A2) = {A1**2}/{A0*A2} = {Fraction(A1**2, A0*A2)}")
print(f"  = {A1**2/(A0*A2):.10f}")
print()

# From the Diophantine identity: 2*A1^2 + 1 = 3*A0*A2
# => A1^2/(A0*A2) = (3*A0*A2 - 1)/(2*A0*A2) = 3/2 - 1/(2*A0*A2)
print(f"From Diophantine: A1^2/(A0*A2) = 3/2 - 1/(2*A0*A2)")
print(f"  = 3/2 - 1/{2*A0*A2} = 3/2 - {Fraction(1, 2*A0*A2)}")
print(f"  = {1.5 - 1/(2*A0*A2):.10f}")
print()

# =============================================================================
# SECTION 13: CONNES-TYPE UNIFICATION APPROACH
# =============================================================================
print("=" * 76)
print("SECTION 13: CONNES-TYPE COUPLING FROM SPECTRAL ACTION")
print("=" * 76)
print()

# In Connes' NCG, for the product M x F with F = Standard Model:
# The spectral action at cutoff Lambda gives:
#   S = 48*f_4*Lambda^4/(pi^2) * int R*sqrt(g)*d^4x
#     + f_0/(2*pi^2) * int [a*|R|^2 + b*(F_mu_nu)^2 + ...] * sqrt(g)*d^4x
#     + ...
#
# where a, b involve traces over the finite Hilbert space H_F.
#
# Specifically: 1/g^2 = 2*f_0*Tr(T_a^2) where T_a are generators.
# For U(1): Tr(Y^2) = 10*N_gen/3 = 10  (with SM hypercharge normalization)
# For SU(2): Tr(T^2) = 2*N_gen = 6
# For SU(3): Tr(T^2) = 2*N_gen = 6
#
# Unification condition: g_1^2 = g_2^2 = g_3^2 at Lambda
# gives sin^2(theta_W) = 3/8 (GUT value, not our b1/(a1^2+1) = 6/26)

# OUR framework differs: sin^2(tW) = 6/26 (not 3/8)
# This means the finite geometry is NOT the standard SM F.
# The McKay graph with Tr(D_F^2) = 8 gives a DIFFERENT structure.

print("Connes SM: sin^2(tW) = 3/8 at unification")
print(f"Our value: sin^2(tW) = b1/(a1^2+1) = {sin2tW:.6f} = 6/26")
print(f"Ratio: (6/26)/(3/8) = {sin2tW/(3/8):.6f}")
print()

# In Connes, the coupling at cutoff Lambda is:
# alpha_GUT = g^2/(4*pi) relates to spectral action as:
# 1/alpha_GUT ~ f_0 * Tr(generators^2)
# For f_0 = 1/(2*pi): 1/alpha_GUT = Tr(T^2)/(2*pi)

# In OUR framework, the spectral action on the 600-cell gives:
# Tr(D^2) = c1 = 14880
# Tr(D^4) = 2*c2 = 111840
# The coupling should relate to some combination of these.

# KEY ATTEMPT: Spectral action coupling extraction
# On a discrete graph, the "gauge field" is the connection on edges.
# The gauge kinetic term comes from: Tr(F^2) where F is the curvature 2-form.
# In the discrete setting: F_{triangle} = holonomy around triangle
# The spectral action gives Tr(F^2) ~ sum over triangles.
#
# For the 600-cell: 1200 triangles (faces).
# Each edge is shared by a1 = 5 triangles.
# The gauge coupling enters as: S_gauge = (1/(4*g^2)) * Tr(F^2)
# Comparing with spectral action: 1/g^2 ~ c2 / (something * Lambda^2)

# Try: 1/alpha_0 = c2 / (c1 * something)
# c2/c1 = A2/A1 = 233/62 = 3.758...
# TARGET/A2*A1 = 137.08*62/233 = 36.46
# Hmm, 36.46 ~ 6*b1 = 36. Close but not exact.

# What if: 1/alpha_0 = A2 * phi^4 / A1?
val_13a = A2 * phi**4 / A1
print(f"A2*phi^4/A1 = {A2}*{phi**4:.6f}/{A1} = {val_13a:.10f}")
print(f"  Ratio to TARGET: {val_13a/TARGET:.10f}")
print(f"  = (233/62)*phi^4 = {233/62*phi**4:.10f}")
print()

# What about: TARGET = A1 * phi^4 * (something)?
# A1*phi^4 = 62*6.854 = 424.95 = 3.1*TARGET ... 3.1 ~ pi
val_13b = A1 * phi**4 / pi
print(f"A1*phi^4/pi = {A1}*{phi**4:.6f}/{pi:.6f} = {val_13b:.6f}")
print(f"  Ratio to TARGET: {val_13b/TARGET:.10f}")
print()

# Try A0*phi^4:
val_13c = A0 * phi**4
print(f"A0*phi^4 = {A0}*{phi**4:.6f} = {val_13c:.6f}")
print(f"  Ratio to TARGET: {val_13c/TARGET:.10f}")
# A0*phi^4 / TARGET = 11*phi^4/(20*phi^4) = 11/20. Not interesting.
print(f"  = A0/(4*a1) = {A0/(4*a1):.6f}")
print()

# =============================================================================
# SECTION 14: COMPREHENSIVE EXHAUSTIVE SEARCH
# =============================================================================
print("=" * 76)
print("SECTION 14: EXHAUSTIVE COMBINATORIAL SEARCH")
print("=" * 76)
print()

# Build a dictionary of ALL framework quantities
all_quantities = {
    # Simplicial spectral data
    "A0": A0,
    "A1": A1,
    "A2": A2,
    "c0": c0,
    "c1": c1,
    "c2": c2,
    "2N": 2*N,
    "N": N,
    # Icosahedral eigenvalues
    "L3": L3,
    "L3p": L3p,
    "L5": L5,
    "L3*L3p": L3*L3p,
    # Framework constants
    "a1": float(a1),
    "b1": float(b1),
    "degree": float(degree),
    "h_E8": float(h_E8),
    "rank8": float(rank_E8),
    "N_gen": float(N_gen),
    "N_eig": float(N_eig),
    # Spectral gaps
    "gap_Cayley": cayley_gap,
    "gap_coexact": lam1_coexact,
    "gap_Hopf": gap_C10,
    "Lambda_max": Lambda_max_D2,
    # McKay
    "Tr_DF2": float(Tr_DF2),
    "dim_F": float(dim_F),
    # Transcendentals
    "pi": pi,
    "2pi": 2*pi,
    # Powers of phi
    "phi": phi,
    "phi2": phi**2,
    "phi3": phi**3,
    "phi4": phi**4,
    "phi5": phi**5,
    "1/phi": 1/phi,
    "1/phi2": 1/phi**2,
    # sqrt(5)
    "sqrt5": sqrt(5),
}

# Now search for all pairs (q1 * q2 = TARGET) or (q1 / q2 = TARGET)
print("Searching for q1 * q2 = TARGET or q1 * q2 * q3 = TARGET ...")
print()

found_pairs = []
for (n1, v1), (n2, v2) in combinations(all_quantities.items(), 2):
    if v2 == 0 or v1 == 0:
        continue
    # Product
    val = v1 * v2
    if abs(val/TARGET - 1) < 1e-8:
        found_pairs.append((f"{n1} * {n2}", val))
    # Ratio v1/v2
    val = v1 / v2
    if abs(val/TARGET - 1) < 1e-8:
        found_pairs.append((f"{n1} / {n2}", val))
    # Ratio v2/v1
    val = v2 / v1
    if abs(val/TARGET - 1) < 1e-8:
        found_pairs.append((f"{n2} / {n1}", val))

if found_pairs:
    print(f"Found {len(found_pairs)} two-quantity matches:")
    for name, val in found_pairs:
        check_match(name, val, "PAIR")
else:
    print("  No two-quantity exact matches found.")
print()

# Three-quantity search (more selective to avoid combinatorial explosion)
print("Searching for q1 * q2 * q3 = TARGET (selective) ...")
print()

# Use a smaller set for the triple search
key_quantities = {
    "A0": A0, "A1": A1, "A2": A2,
    "a1": float(a1), "b1": float(b1), "N": float(N), "degree": float(degree),
    "L3": L3, "L3p": L3p, "L5": L5,
    "phi": phi, "phi2": phi**2, "phi3": phi**3, "phi4": phi**4,
    "1/phi": 1/phi, "1/phi2": 1/phi**2,
    "sqrt5": sqrt(5),
    "gap_C": cayley_gap, "gap_cx": lam1_coexact, "gap_H": gap_C10,
    "pi": pi, "2pi": 2*pi,
    "Tr_DF2": float(Tr_DF2), "dim_F": float(dim_F),
    "rank8": float(rank_E8), "N_gen": float(N_gen),
    "h_E8": float(h_E8),
}

found_triples = []
keys = list(key_quantities.keys())
vals = list(key_quantities.values())

for i in range(len(keys)):
    for j in range(i, len(keys)):
        for k in range(j, len(keys)):
            v = vals[i] * vals[j] * vals[k]
            if v != 0 and abs(v/TARGET - 1) < 1e-8:
                found_triples.append((f"{keys[i]} * {keys[j]} * {keys[k]}", v))

if found_triples:
    print(f"Found {len(found_triples)} three-quantity product matches:")
    for name, val in found_triples:
        check_match(name, val, "TRIPLE")
else:
    print("  No three-quantity product matches found.")
print()

# Also search q1 * q2 / q3
print("Searching for q1 * q2 / q3 = TARGET ...")
print()

found_triple_ratio = []
for i in range(len(keys)):
    for j in range(i, len(keys)):
        for k in range(len(keys)):
            if vals[k] == 0:
                continue
            v = vals[i] * vals[j] / vals[k]
            if abs(v/TARGET - 1) < 1e-8:
                name = f"{keys[i]} * {keys[j]} / {keys[k]}"
                # Skip trivial reformulations
                if "phi4" in name and ("a1" in name or "b1" in name) and len(set([keys[i],keys[j],keys[k]])) == 3:
                    found_triple_ratio.append((name, v))
                elif "phi4" not in name:
                    found_triple_ratio.append((name, v))

# Deduplicate
seen = set()
unique_triples = []
for name, val in found_triple_ratio:
    if name not in seen:
        seen.add(name)
        unique_triples.append((name, val))

if unique_triples:
    print(f"Found {len(unique_triples)} q1*q2/q3 matches (showing first 30):")
    for name, val in unique_triples[:30]:
        check_match(name, val, "RATIO3")
else:
    print("  No q1*q2/q3 matches found.")
print()

# =============================================================================
# SECTION 15: THE KEY INSIGHT -- SPECTRAL ACTION ON PRODUCT GEOMETRY
# =============================================================================
print("=" * 76)
print("SECTION 15: PRODUCT GEOMETRY SPECTRAL ACTION")
print("=" * 76)
print()

# For the product D = D_M tensor 1_F + gamma tensor D_F:
# Tr(D^2) = Tr(D_M^2)*dim_F + dim_M*Tr(D_F^2)
# = c1*9 + c0*8 = 14880*9 + 2640*8 = 133920 + 21120 = 155040
#
# Tr(D^4) = Tr(D_M^4)*dim_F + 2*Tr(D_M^2)*Tr(D_F^2) + dim_M*Tr(D_F^4)
# = 2*c2*9 + 2*c1*8 + c0*Tr(D_F^4)

c1_prod = c1 * dim_F + c0 * Tr_DF2
print(f"Product c1 = c1_M*dim_F + c0_M*Tr(D_F^2)")
print(f"  = {c1}*{dim_F} + {c0}*{Tr_DF2}")
print(f"  = {c1*dim_F} + {c0*Tr_DF2}")
print(f"  = {c1_prod}")
print()

# Ratio of product coefficients
c0_prod = c0 * dim_F
ratio_prod = c1_prod / c0_prod
print(f"Product c0 = c0_M*dim_F = {c0}*{dim_F} = {c0_prod}")
print(f"Product c1/c0 = {ratio_prod:.10f}")
print()

# CRITICAL: Does product c1/c0 relate to 1/alpha?
print(f"  c1_prod/c0_prod = {ratio_prod:.10f}")
print(f"  = (c1*dim_F + c0*Tr_DF2) / (c0*dim_F)")
print(f"  = c1/c0 + Tr_DF2/dim_F")
print(f"  = {c1/c0:.6f} + {Tr_DF2/dim_F:.6f}")
print(f"  = {c1/c0 + Tr_DF2/dim_F:.10f}")
print()

# c1/c0 = 28/5, Tr_DF2/dim_F = 8/9
# Sum = 28/5 + 8/9 = 252/45 + 40/45 = 292/45
ratio_exact = Fraction(28*9 + 8*5, 5*9)
print(f"  = {ratio_exact} = {float(ratio_exact):.10f}")
print(f"  TARGET / this = {TARGET/float(ratio_exact):.10f}")
print()

# That's 6.489 and TARGET = 137.08, so ratio = 21.12...
# 21.12 ~ 8*phi^2 = 20.94... not quite
# 21.12 ~ 4*a1*(phi^4-1)/... nah

# Let me try a DIFFERENT product: the ratio of gravity to matter parts
grav_part_c1 = c1 * dim_F      # 133920
matt_part_c1 = c0 * Tr_DF2     # 21120
ratio_gm = grav_part_c1 / matt_part_c1
print(f"Gravity/Matter in c1_product = {grav_part_c1}/{matt_part_c1} = {ratio_gm:.10f}")
print(f"  = (c1*dim_F)/(c0*Tr_DF2) = (A1*dim_F)/(A0*Tr_DF2)")
print(f"  = ({A1}*{dim_F})/({A0}*{Tr_DF2}) = {A1*dim_F}/{A0*Tr_DF2}")
print(f"  = {Fraction(A1*dim_F, A0*Tr_DF2)} = {A1*dim_F/(A0*Tr_DF2):.10f}")
print()

# 62*9/(11*8) = 558/88 = 279/44 = 6.34090909...
# TARGET / 6.34 = 21.62 ~ not clean

# =============================================================================
# SECTION 16: THE CONNES COUPLING EXTRACTION FORMULA
# =============================================================================
print("=" * 76)
print("SECTION 16: CONNES COUPLING FORMULA ON 600-CELL")
print("=" * 76)
print()

# In Connes' spectral action, the gauge coupling is extracted from:
# g^{-2} = f_0 * Tr_{H_F}(T_R^2)
# where T_R is the generator in representation R and f_0 = moments of cutoff.
#
# For the 600-cell framework:
# f_0 = N^0/(2*pi) = 1/(2*pi)  (from exp466)
# f_2 = N^2/(2*pi) = 14400/(2*pi)
# f_4 = N^4/(2*pi) = 2.0736e8/(2*pi)
#
# Then: 1/alpha_0 = (f_0/pi) * Tr(hypercharge^2 on 600-cell spectrum)
# The question is: what is Tr(Y^2) on the 600-cell?

f0 = 1 / (2*pi)
f2 = N**2 / (2*pi)
f4 = N**4 / (2*pi)

print("Test function moments (from exp466):")
print(f"  f_0 = 1/(2*pi) = {f0:.10f}")
print(f"  f_2 = N^2/(2*pi) = {f2:.6f}")
print(f"  f_4 = N^4/(2*pi) = {f4:.6f}")
print()

# If 1/alpha = f_0 * Tr(Y^2) / pi:
# TARGET = 4*a1*phi^4
# => Tr(Y^2) = TARGET * pi / f_0 = TARGET * pi * 2*pi = TARGET * 2*pi^2
Tr_Y2_needed = TARGET * pi / f0
print(f"If 1/alpha_0 = f_0*Tr(Y^2)/pi:")
print(f"  Tr(Y^2) = TARGET*pi/f_0 = {Tr_Y2_needed:.6f}")
print(f"  = TARGET*2*pi^2 = {TARGET*2*pi**2:.6f}")
print(f"  = {TARGET*2*pi**2:.2f}")
print(f"  This is ~2704. No obvious meaning.")
print()

# Alternative: 1/alpha = (2*f_0) * Tr(Y^2) (Connes normalization)
Tr_Y2_alt = TARGET / (2*f0)
print(f"If 1/alpha_0 = 2*f_0*Tr(Y^2):")
print(f"  Tr(Y^2) = TARGET/(2*f_0) = {Tr_Y2_alt:.6f}")
print(f"  = TARGET*pi = {TARGET*pi:.6f}")
print(f"  = {TARGET*pi:.2f}")
print(f"  ~ 430.5. Not obvious.")
print()

# What if the test function is NOT f_0=1/(2*pi)?
# What if f_0 is determined by the spectral action itself?
# E.g., f_0 = 1/c0 or f_0 = something from the spectrum?

# Try: 1/alpha_0 = c1 * Tr(D_F^2) / c0
val_16 = c1 * Tr_DF2 / c0
print(f"c1*Tr(D_F^2)/c0 = {c1}*{Tr_DF2}/{c0} = {val_16:.10f}")
print(f"  = A1*Tr(D_F^2)/A0 = {A1}*{Tr_DF2}/{A0} = {A1*Tr_DF2/A0:.10f}")
print(f"  TARGET = {TARGET:.10f}")
print(f"  Ratio: {val_16/TARGET:.10f}")
print()

# Try: 1/alpha_0 = c2 * something / (c1 * something_else)
# c2/c1 = 233/62. We need TARGET * 62/233 = 36.461...
# Hmm, 36.461 ~ 6*b1 + 0.46... not clean

# =============================================================================
# SECTION 17: THE DEFINITIVE ROUTE (if it exists)
# =============================================================================
print("=" * 76)
print("SECTION 17: DEFINITIVE SPECTRAL ROUTE ANALYSIS")
print("=" * 76)
print()

# Let me carefully analyze what we have:
#
# 4*a1*phi^4 = L(3)*L(3') * phi^4  [L(3)*L(3') from icosahedral Laplacian]
#            = 4*a1 * phi^4          [since L(3)*L(3')=a1(a1-1)=4a1 for a1=5]
#            = a1*(degree/gap)^2     [from Cayley spectral expansion ratio]
#            = L(3')^3/L(3)          [from icosahedral eigenvalue asymmetry]
#
# All these are EQUIVALENT algebraic identities. The question is:
# which one connects to the spectral ACTION (the operator trace)?
#
# The spectral action coefficients c0,c1,c2 are INDEPENDENT of the
# icosahedral Laplacian eigenvalues L(3),L(3'),L(5). These are:
# - c0,c1,c2: from the HODGE-DIRAC operator on the whole 600-cell
# - L(k): from the LAPLACIAN on the VERTEX FIGURE (icosahedron)
#
# The connection must come through the LOCAL geometry:
# The vertex figure determines the LOCAL structure, while the
# Hodge-Dirac gives the GLOBAL spectral data.
#
# In spectral geometry, the Seeley-DeWitt coefficients encode
# local curvature invariants integrated over the manifold.
# For a homogeneous space (like S^3), the local and global are related.

# KEY: c1 = Tr(D^2) = sum of all Hodge Laplacian eigenvalues
# On a graph, Tr(Delta_0) = sum of degrees = N*degree = 120*12 = 1440
# But c1 = 14880 includes Delta_0 + Delta_1 + Delta_2 + Delta_3

# The vertex figure eigenvalues appear in the EDGE Laplacian decomposition.
# Delta_1 splits into exact + coexact parts.
# The coexact part relates to the link (vertex figure) Laplacian.

# ROUTE: Can we factor c1 as something involving L(3)*L(3')?
# c1/N = 14880/120 = 124
# 124 = 2*A1 = 2*(2*a1^2+2*a1+2) = 2*62
# Also: 124 = 4*31 = 4*(a1^2+a1+1)
# Also: 124 = a1^2 + a1^2 - 1 + N ... no
# 124 = 4*a1*(a1+...hmm

# Direct: does c1/(N*L(3)*L(3')) give phi^4?
val_17a = c1 / (N * L3 * L3p)
print(f"c1/(N*L(3)*L(3')) = {c1}/({N}*{L3*L3p:.0f}) = {c1}/{N*int(L3*L3p)} = {val_17a:.10f}")
print(f"  phi^4 = {phi**4:.10f}")
print(f"  Ratio: {val_17a/phi**4:.10f}")
# c1/(N*20) = 14880/2400 = 6.2
# phi^4 = 6.854
# Not a match.
print()

# Try c1/(2*N*L(3)*L(3')) = 14880/4800 = 3.1
val_17b = c1 / (2*N*L3*L3p)
print(f"c1/(2*N*L(3)*L(3')) = {val_17b:.10f}")
print(f"  phi^3 = {phi**3:.10f}")
print(f"  pi = {pi:.10f}")
print(f"  A1/L(3)*L(3') = {A1/(L3*L3p):.10f} = {A1}/20 = {Fraction(A1,20)}")
# 62/20 = 31/10. Not interesting.
print()

# Alternative factoring of TARGET:
# 4*a1*phi^4 = 20*(3*phi+2) = 60*phi + 40
# In terms of N: = N/2 * phi + N/3 ... no, 60=N/2, 40=N/3
# 60 = N/2, 40 = N/3
# So: TARGET = N*phi/2 + N/3 = N*(phi/2 + 1/3) = N*(3*phi+2)/6
val_17c = N * (3*phi+2) / 6
print(f"N*(3*phi+2)/6 = {N}*(3*{phi:.4f}+2)/6 = {val_17c:.10f}")
print(f"  = N*phi^4/b1 = {N*phi**4/b1:.10f}")
print(f"  TARGET = {TARGET:.10f}")
print(f"  MATCH: {abs(val_17c-TARGET)<1e-8}")
print()

# So: TARGET = N*phi^4/b1 = |2I| * phi^4 / (a1+1)
# This IS the standard identity: 4*a1*phi^4 = N*phi^4/b1
# since N = 120 = 4*a1*(a1+1)/1 = 4*a1*b1... wait
# N = a1! = 120. 4*a1 = 20. phi^4/1 * 20 = 20*phi^4.
# N/b1 = 120/6 = 20 = 4*a1. So this is the same identity.
# Already known: the linear coefficient of alpha equation is (N/b1)*phi^4.

print("The alpha equation in standard form:")
print(f"  2*pi*alpha^2 - (N/b1)*phi^4*alpha + 1 = 0")
print(f"  where N/b1 = {N}/{b1} = {N//b1} = 4*a1")
print(f"  This is the established form. N/b1 from spectral action needs:")
print(f"    N = |2I| = 120 (size of the group)")
print(f"    b1 = a1+1 = 6 (homological/Betti)")
print(f"    phi^4 = golden ratio to the 4th power")
print()

# =============================================================================
# SECTION 18: SPECTRAL DETERMINANT ROUTE
# =============================================================================
print("=" * 76)
print("SECTION 18: SPECTRAL DETERMINANT AND ZETA")
print("=" * 76)
print()

# On a compact Riemannian manifold, the effective coupling at low energy
# relates to the spectral zeta function:
# alpha(mu) = alpha_0 * Z(mu/Lambda)
# where Z involves det'(D^2/Lambda^2)
#
# For the DISCRETE case, we can compute everything:
# log det' = sum_{lambda>0} log(lambda)
# zeta(s) = sum_{lambda>0} lambda^{-s}
#
# From the Seeley-DeWitt expansion:
# zeta(s) ~ c0/s + c1 + ... (schematically)
#
# The ANALYTIC TORSION T relates to log det':
# log T = (1/2) * sum_p (-1)^{p+1} * p * log det'(Delta_p)

# We know from exp463:
# Lambda_max(D^2) = b1*phi^2 = 15.7082...
# The spectral gap and maximum determine the "spectral window"

# Spectral action at a specific "scale" t:
# S(t) = Tr(exp(-t*D^2)) = c0 - c1*t + c2*t^2 - ...
# At t = 1/Lambda_max: S(t) ~ c0*(1 - c1/(c0*Lambda_max) + ...)
# The coupling emerges at the scale where S(t) = something

# Critical: at what scale t does the spectral action see alpha?
# If t* = 1/(4*a1*phi^4)^2 = 1/TARGET^2 (i.e., Lambda = TARGET):
# Then S(t*) involves c0 - c1*t* + c2*t*^2
t_star = 1/TARGET**2
S_tstar = c0 - c1*t_star + c2*t_star**2
print(f"Heat kernel at t* = 1/TARGET^2 = {t_star:.10e}:")
print(f"  S(t*) = {c0} - {c1}*{t_star:.6e} + {c2}*{t_star**2:.6e}")
print(f"  = {c0} - {c1*t_star:.6f} + {c2*t_star**2:.10f}")
print(f"  = {S_tstar:.6f}")
print()

# =============================================================================
# SECTION 19: THE MOST PROMISING ROUTE -- SPECTRAL ACTION ON LOCAL GEOMETRY
# =============================================================================
print("=" * 76)
print("SECTION 19: LOCAL-GLOBAL CONNECTION")
print("=" * 76)
print()

# The VERTEX FIGURE (icosahedron, 12 vertices) is the LOCAL structure.
# The 600-CELL (120 vertices) is the GLOBAL structure.
# Connection: N = degree * (N_vertices/degree) but this is trivial.
#
# More subtle: the Seeley-DeWitt coefficients of the 600-cell
# FACTORIZE through the vertex figure spectrum.
#
# For a homogeneous simplicial complex (like the 600-cell), the local
# Laplacian (on the link/vertex figure) determines the global spectrum
# through a "transfer matrix" or "representation-theoretic" decomposition.
#
# The key relation is the GAUSS-BONNET type:
# c1 = Tr(D^2) = sum over all p-simplices of their "local curvature contribution"
# For a vertex-transitive complex:
# c1 = N * sum_p f_p/N * <trace of link Laplacian on p-simplices>
#
# Specifically for the Hodge Laplacian on the 600-cell:
# c1 = N * [degree + edge_contribution + face_contribution + cell_contribution]
# where each contribution involves the vertex figure eigenvalues.

# The vertex figure (icosahedron) has:
# Adjacency eigenvalues: a1, sqrt(a1), 1, -1, -sqrt(a1), -a1
# (with multiplicities 1, 3, 5, 5, 3, 1)
# Wait, the icosahedron has 12 vertices, so its adjacency has 4 distinct eigs
# (by the irrep decomposition of A5):
# lambda_1=5, lambda_3=sqrt(5), lambda_5=1, lambda_3'=-sqrt(5)...
# Actually for icosahedron adjacency (5-regular):
# Eigenvalues: 5 (x1), sqrt(5) (x3), 1 (x5), -1 (x5), -sqrt(5) (x3), -5...
# No wait. Icosahedron has 12 vertices, degree 5.
# Its adjacency eigenvalues are: 5, sqrt(5), sqrt(5), 1, 1, -1, -1, -sqrt(5), -sqrt(5), -5,...
# Actually there are 12 eigenvalues. By A5 symmetry:
# 5(x1), sqrt(5)(x3), -1(x5), -sqrt(5)(x3) ... that's 1+3+5+3=12. Yes!

ico_adj_eigs = {5: 1, sqrt(5): 3, -1: 5, -sqrt(5): 3}
print("Icosahedron adjacency eigenvalues:")
total_ico = 0
for eig, mult in sorted(ico_adj_eigs.items(), reverse=True):
    print(f"  {eig:+8.4f}  (mult {mult})")
    total_ico += mult
print(f"  Total: {total_ico}")
print()

# Icosahedron Laplacian eigenvalues (degree - adjacency):
# 0(x1), 5-sqrt(5)(x3), 6(x5), 5+sqrt(5)(x3)
# Confirming L(3) = 5-sqrt(5), L(5)=6, L(3')=5+sqrt(5)

# KEY: The spectral zeta of the icosahedron Laplacian at s=-2:
# zeta_ico(-2) = sum lambda_i^2 (excluding zero mode)
zeta_ico_neg2 = 3*(a1-sqrt(a1))**2 + 5*b1**2 + 3*(a1+sqrt(a1))**2
print(f"Icosahedron Laplacian: zeta(-2) = sum lambda^2")
print(f"  = 3*L(3)^2 + 5*L(5)^2 + 3*L(3')^2")
print(f"  = 3*{(a1-sqrt(a1))**2:.4f} + 5*{b1**2} + 3*{(a1+sqrt(a1))**2:.4f}")
print(f"  = {3*(a1-sqrt(a1))**2:.4f} + {5*b1**2} + {3*(a1+sqrt(a1))**2:.4f}")
print(f"  = {zeta_ico_neg2:.6f}")
# 3*(25-10sqrt5+5) + 5*36 + 3*(25+10sqrt5+5)
# = 3*30 - 30sqrt5 + 180 + 3*30 + 30sqrt5 = 90 + 180 + 90 = 360
print(f"  = 3*(a1^2+a1) + 5*b1^2 + 3*(a1^2+a1)")
print(f"  = 6*(a1^2+a1) + 5*b1^2 = 6*30 + 5*36 = {6*30 + 5*36}")
print(f"  = {int(zeta_ico_neg2)}")
print()

# zeta(-2) = 360 = N*N_gen = 120*3
print(f"  zeta_ico(-2) = 360 = N*N_gen = {N}*{N_gen}")
print()

# zeta(-1):
zeta_ico_neg1 = 3*(a1-sqrt(a1)) + 5*b1 + 3*(a1+sqrt(a1))
print(f"Icosahedron Laplacian: zeta(-1) = Tr(L_ico\\0)")
print(f"  = 3*L(3) + 5*L(5) + 3*L(3')")
print(f"  = 3*{a1-sqrt(a1):.4f} + 5*{b1} + 3*{a1+sqrt(a1):.4f}")
print(f"  = {3*(a1-sqrt(a1)):.4f} + {5*b1} + {3*(a1+sqrt(a1)):.4f}")
# = 3*5-3sqrt5 + 30 + 3*5+3sqrt5 = 15+30+15 = 60
print(f"  = {zeta_ico_neg1:.0f}")
print(f"  = N/2 = {N//2}")
print()

# Now: does zeta_ico(-1) * phi^4 / b1 give something?
val_19a = zeta_ico_neg1 * phi**4 / b1
print(f"zeta_ico(-1)*phi^4/b1 = {zeta_ico_neg1:.0f}*{phi**4:.6f}/{b1}")
print(f"  = (N/2)*phi^4/b1 = {val_19a:.10f}")
print(f"  TARGET/2 = {TARGET/2:.10f}")
print(f"  MATCH: {abs(val_19a - TARGET/2) < 1e-8}")
print()
# N/2 * phi^4/b1 = 60*phi^4/6 = 10*phi^4 = TARGET/2. Not full match.

# OK so TARGET = 2 * zeta_ico(-1) * phi^4 / b1
# The factor 2 could come from: double cover (2I vs I), or Tr/N relation

# But: TARGET = N*phi^4/b1 = 2*(N/2)*phi^4/b1 = 2*zeta_ico(-1)*phi^4/b1
# Since zeta_ico(-1) = N/2, this is just 20*phi^4 again. CIRCULAR.

# =============================================================================
# SECTION 20: THE HONEST ASSESSMENT
# =============================================================================
print("=" * 76)
print("SECTION 20: HONEST ASSESSMENT")
print("=" * 76)
print()

print("SUMMARY OF ALL ROUTES TO 4*a1*phi^4:")
print()
print("=" * 76)
print()

routes = [
    ("a1*(degree/gap_Cayley)^2 = a1*(2*phi^2)^2",
     "STRUCTURAL",
     "Uses a1 as a free parameter and the Cayley graph spectral ratio.\n"
     "     The spectral ratio degree/gap = 2*phi^2 IS from spectral data.\n"
     "     But the factor a1 is the input parameter, not derived from Tr(D^k)."),

    ("L(3)*L(3')*phi^4 = 4*a1*phi^4",
     "STRUCTURAL",
     "L(3)*L(3') comes from the icosahedral Laplacian (vertex figure).\n"
     "     The phi^4 is (L(3')/L(3))^2, also from vertex figure.\n"
     "     This IS a spectral identity, but on the VERTEX FIGURE, not the full D."),

    ("L(3')^3/L(3)",
     "STRUCTURAL",
     "Same as above: purely from vertex figure eigenvalues.\n"
     "     Algebraically equivalent to L(3)*L(3')*(L(3')/L(3))^2."),

    ("N*phi^4/b1 = (N/b1)*phi^4",
     "STRUCTURAL",
     "From the alpha equation standard form. N and b1 are framework constants.\n"
     "     phi^4 comes from the eigenvalue structure. But N/b1=20 is just 4*a1."),

    ("c1/(2*N*L(3)*L(3')) * TARGET != phi^4",
     "NEGATIVE",
     "Attempted to extract phi^4 from ratio of c1 (global) to L(3)*L(3') (local).\n"
     "     Gives 3.1, not phi^4=6.854. The Seeley-DeWitt coefficients\n"
     "     do NOT factorize through the vertex figure eigenvalues."),

    ("Product geometry c1_prod/c0_prod",
     "NEGATIVE",
     "The product M x F spectral ratio gives 292/45 = 6.489,\n"
     "     which does NOT equal TARGET/something-simple."),

    ("Connes coupling: f_0 * Tr(Y^2)",
     "NEGATIVE",
     "Requires Tr(Y^2) ~ 2704 or 431, neither of which has clean meaning\n"
     "     in the 600-cell framework."),
]

for i, (route, cat, detail) in enumerate(routes, 1):
    print(f"  ROUTE {i}: {route}")
    print(f"  Category: [{cat}]")
    print(f"     {detail}")
    print()

print("=" * 76)
print("CRITICAL FINDING:")
print("=" * 76)
print()

print("The linear coefficient 4*a1*phi^4 CANNOT be extracted from the")
print("Seeley-DeWitt coefficients (c0, c1, c2) alone. The spectral action")
print("on the full Hodge-Dirac operator D = d + d* on the 600-cell gives")
print("integer coefficients (2640, 14880, 55920) that encode TOPOLOGICAL")
print("data (Euler characteristic, Betti numbers) but NOT the golden ratio")
print("structure needed for the coupling constant.")
print()
print("The phi^4 factor comes from the VERTEX FIGURE (icosahedral Laplacian),")
print("which is a LOCAL geometric quantity -- the LINK structure at each")
print("vertex. The Seeley-DeWitt coefficients are GLOBAL sums.")
print()
print("On a smooth manifold, the local-global connection is provided by")
print("the asymptotic expansion of the heat kernel: a_k(x) are LOCAL")
print("curvature invariants, and c_k = integral of a_k(x). For a homogeneous")
print("space, c_k = Vol(M) * a_k (constant on M). But for the DISCRETE case,")
print("the link Laplacian eigenvalues enter through the COMBINATORIAL")
print("Laplacians, not through the spectral action coefficients directly.")
print()

# =============================================================================
# SECTION 21: THE RIGHT QUESTION
# =============================================================================
print("=" * 76)
print("SECTION 21: WHAT WOULD A SPECTRAL DERIVATION LOOK LIKE?")
print("=" * 76)
print()

print("For a true spectral-action derivation of 4*a1*phi^4, we would need:")
print()
print("1. A SPECTRAL ACTION FUNCTIONAL S[D] = Tr(chi(D^2/Lambda^2))")
print("   on the product geometry D = D_600cell tensor 1_McKay + gamma tensor D_McKay")
print()
print("2. EXTRACT the gauge coupling from the Seeley-DeWitt expansion:")
print("   1/g^2 = f_0 * (c_2_gauge / Lambda^4)")
print("   where c_2_gauge is the coefficient of the F^2 term")
print()
print("3. The F^2 coefficient comes from the FOURTH moment of D_total^2:")
print("   c_2_total = c2_M*dim_F + c1_M*Tr(D_F^2) + c0_M*Tr(D_F^4)/2")
print()
print("4. The coupling would then be:")
print("   1/alpha = (some function of c2_total, c1_total, c0_total, Lambda)")
print()

# Let's check if the MIXED term c1_M * Tr(D_F^2) gives TARGET
mixed_term = c1 * Tr_DF2
print(f"Mixed term c1_M * Tr(D_F^2) = {c1} * {Tr_DF2} = {mixed_term}")
print(f"  / c0 = {mixed_term/c0:.10f}")
print(f"  / N = {mixed_term/N:.10f}")
print(f"  / (2*N) = {mixed_term/(2*N):.10f}")
print(f"  TARGET = {TARGET:.10f}")
print()

# c1*8/c0 = 14880*8/2640 = 119040/2640 = 45.09
# c1*8/N = 14880*8/120 = 992
# c1*8/(2N) = 496

# What about: mixed_term / (c0 * phi^4)?
val_21 = mixed_term / (c0 * phi**4)
print(f"  c1*Tr_DF2 / (c0*phi^4) = {val_21:.10f}")
# = 119040 / (2640*6.854) = 119040/18094.6 = 6.578
# Not TARGET.

# What about: mixed_term * phi^4 / (N * something)?
# 119040 * 6.854 / 120 = 6797. Too big.
# 119040 * phi^4 / (c1) = 8*phi^4 = 54.83
# 119040 * phi^4 / (c2) = 14.59

# What if we use Tr(D_F^4)?
# For the McKay graph of 2I:
# If D_F eigenvalues satisfy: Tr(D_F^2)=8, and the eigenvalues are related to
# the Dynkin diagram structure, then Tr(D_F^4) depends on the specific D_F.
# From exp262: b_Y = Tr(D_F^4) was computed there. Let me check what value is natural.

# For the affine E8 Dynkin diagram (McKay graph), the adjacency matrix has
# specific eigenvalues. The D_F that gives Tr(D_F^2)=8 could be:
# D_F eigenvalues: {0, +/-1, +/-1, +/-sqrt(2), +/-sqrt(2)}
# -> Tr(D_F^2) = 0+2+2+4+4 = 12. No.
# Or: some mass exponents. In our framework, the mass formula gives
# n_f values that square-sum to something.
# Actually from the paper: the n_f values (mass quantum numbers) satisfy
# ||Y||_F^2 = Tr(D_F^2) = 8.

# If D_F^2 has eigenvalues proportional to n_f^2, and sum(n_f^2) = 8,
# then a natural set is: the squared norms of the (a,b) mass quantum numbers
# For (a,b) pairs: e(0,0), mu(1,1), tau(1,2), u(3,-2), c(2,1), t(4,1), d(1,0), s(1,1), b(-1,4)
# n_f^2 = a^2 + a*b + b^2 (possibly)... or 5*a + 6*b + delta...
# The actual mass exponents in the formula are: m_f = m_e * phi^(5a+6b+delta)
# The mass "eigenvalue" is the exponent n = 5a+6b+delta.
# From the (a,b) pairs:
ab_pairs = [(0,0), (1,1), (1,2), (3,-2), (2,1), (4,1), (1,0), (1,1), (-1,4)]
labels = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
exponents = [5*a + 6*b for a, b in ab_pairs]

print("Mass exponents (tree-level, no delta):")
for label, (a,b), n in zip(labels, ab_pairs, exponents):
    print(f"  {label:3s}: (a,b)=({a:+d},{b:+d})  n=5a+6b={n:+d}")
sum_n2 = sum(n**2 for n in exponents)
print(f"Sum of n^2 = {sum_n2}")
print(f"  Tr(D_F^2) (if eigenvalues are n_f) = {sum_n2}")
print(f"  But the KNOWN value is Tr(D_F^2) = 8 = rank(E8)")
print(f"  So the mass exponents are NOT the D_F eigenvalues directly.")
print(f"  The D_F that gives Tr(D_F^2)=8 encodes the STRUCTURE of the McKay graph,")
print(f"  not the mass exponents (which come from Wilson lines on D_F).")
print()

# =============================================================================
# SECTION 22: FINAL COMPREHENSIVE SCAN
# =============================================================================
print("=" * 76)
print("SECTION 22: RATIO SCAN - c_k WITH EIGENVALUE FUNCTIONS")
print("=" * 76)
print()

# The most general ansatz:
# TARGET = c_a^p * c_b^q / c_c^r  * phi^s * (small rational)
# where a,b,c in {0,1,2} and p,q,r,s are small integers

# Since TARGET ~ 137 and c0=2640, c1=14880, c2=55920:
# Let's systematically check c_a/c_b * phi^s for small s

print("Checking c_a/c_b * phi^s * (integer/integer):")
print()

for (na, ca), (nb, cb) in [
    (("c0",c0), ("1",1)),
    (("c1",c1), ("1",1)),
    (("c2",c2), ("1",1)),
    (("c1",c1), ("c0",c0)),
    (("c2",c2), ("c0",c0)),
    (("c2",c2), ("c1",c1)),
    (("c0",c0), ("c1",c1)),
    (("c0",c0), ("c2",c2)),
    (("c1",c1), ("c2",c2)),
    (("c1^2",c1**2), ("c0*c2",c0*c2)),
]:
    base = ca / cb if cb != 0 else 0
    if base == 0:
        continue
    for s in range(-6, 7):
        val = base * phi**s
        if val == 0:
            continue
        ratio = TARGET / val
        if 0.01 < abs(ratio) < 1000:
            # Check if ratio is close to a simple fraction p/q with small p,q
            for denom in range(1, 50):
                numer = round(ratio * denom)
                if abs(numer) > 200:
                    continue
                if numer == 0:
                    continue
                test_val = base * phi**s * numer / denom
                err = abs(test_val/TARGET - 1)
                if err < 1e-10:
                    desc = f"({na}/{nb})*phi^{s}*{numer}/{denom}"
                    if len(desc) < 60:
                        desc = desc.ljust(60)
                    print(f"  *** {desc} = {test_val:.10f} (err={err:.2e})")

print()

# =============================================================================
# SECTION 23: DEEP SPECTRAL FACTORIZATION
# =============================================================================
print("=" * 76)
print("SECTION 23: DEEP ANALYSIS - A_k FACTORIZATION")
print("=" * 76)
print()

# A0 = 11 = 2*a1+1
# A1 = 62 = 2*a1^2+2*a1+2 = 2*(a1^2+a1+1)
# A2 = 233 = F_13 (Fibonacci!)
# A2 = 6*a1^2+15*a1+8

# The fact that A2=233=F_13 is suggestive. F_n/F_{n-1} -> phi.
# phi^4 ~ F_5/F_4 = 5/3... no. phi^4 = 6.854 != 5/3.
# Actually phi^n = F_n*phi + F_{n-1} (the Fibonacci representation)
# phi^4 = F_4*phi + F_3 = 3*phi + 2 = 6.854...

# Does A2 relate to TARGET? A2*phi^4/A1 already computed above.
# A2/A1 = 233/62 = 3.7580... TARGET needs 137.08.
# A2*phi^4/A1 = 233*6.854/62 = 25.76. Not TARGET.
# A2*phi^4/A0 = 233*6.854/11 = 145.16. Close to TARGET? Off by 6%.
val_23a = A2 * phi**4 / A0
print(f"A2*phi^4/A0 = {A2}*{phi**4:.6f}/{A0} = {val_23a:.6f}")
print(f"  TARGET = {TARGET:.6f}")
print(f"  Ratio: {val_23a/TARGET:.6f} ({(val_23a/TARGET-1)*100:.2f}%)")
print()

# Interesting: 145.16/137.08 = 1.059 ~ 6/A0*phi? Nah.
# Let me try: A2*phi^k for various k
print("A2*phi^k for various k:")
for k in range(-5, 6):
    val = A2 * phi**k
    print(f"  A2*phi^{k:+d} = {val:12.4f}  (ratio to TARGET: {val/TARGET:.4f})")
print()

# A2*phi^{-1} = 233/1.618 = 143.95. Still 5% off.
# A2*phi^{-2} = 233/2.618 = 89.00 = F_11! (Because A2=F_13 and phi^{-2}=F_11/F_13... approximately)
# Actually F_11=89, F_13=233. phi^2 = 233/89.0... no. phi^2 = 2.618. 233/2.618 = 89.0 exactly!
# Because F_{n}/F_{n-2} = phi^2 (asymptotically), and for Fibonacci: F_13/phi^2 = F_11 exactly? NO.
# F_13 = 233, F_11 = 89. 233/89 = 2.617978... while phi^2 = 2.618034...
# Close but NOT exact. The Fibonacci ratio converges to phi.

print(f"A2/phi^2 = {A2/phi**2:.10f}")
print(f"F_11 = 89")
print(f"Difference: {A2/phi**2 - 89:.10f} (NOT zero)")
print()

# =============================================================================
# SECTION 24: THE VERTEX FIGURE SPECTRAL ACTION
# =============================================================================
print("=" * 76)
print("SECTION 24: SPECTRAL ACTION ON THE VERTEX FIGURE (ICOSAHEDRON)")
print("=" * 76)
print()

# What if we compute the spectral action directly on the icosahedron?
# The icosahedron has 12 vertices, 30 edges, 20 faces.
# Hilbert space: 12 + 30 + 20 = 62 = A1!

print(f"Icosahedron f-vector: (12, 30, 20)")
print(f"  dim(H_ico) = 12+30+20 = 62 = A1!")
print(f"  This is a KNOWN identity: A1 = dim(H_vertex_figure)")
print()

# The link between A1 and the vertex figure is STRUCTURAL:
# A1 = dim(H_icosahedron) = number of simplices in the vertex figure
# This means: c1/(2*N) = A1 = dim(H_ico)

# Now: the spectral action on the icosahedron would give its own coefficients.
# For the icosahedron (12 vertices, 30 edges, 20 faces, S^2 topology):
# Euler char = 12 - 30 + 20 = 2 (S^2)
# Betti: b0=1, b1=0, b2=1

# Icosahedron Seeley-DeWitt (from its Hodge Laplacian):
# c0_ico = 62 (dim of Hilbert space)
# c1_ico = Tr(D_ico^2)
# For the icosahedral Laplacian on 0-forms (12x12):
# Tr(Delta_0_ico) = 12*5 = 60 (all vertices have degree 5)

# On 1-forms (edges): 30 edges, each in some triangles
# The edge Laplacian of the icosahedron...
# Let me just compute: Tr(Delta_p) for the icosahedron.

# Delta_0: 12x12 graph Laplacian, Tr = 60
# Delta_2: 20x20, dual to Delta_0 on sphere, Tr(Delta_2) = sum of dual degrees
# Each face (triangle) shares 3 edges, and each edge is in 2 faces.
# For 2-forms on S^2: Delta_2 eigenvalues same as Delta_0 (by Poincare duality)
# So Tr(Delta_2) = 60 (same as Tr(Delta_0))
# Wait: Poincare duality gives isospectral Delta_p ~ Delta_{d-p} for d=2
# So Delta_0 ~ Delta_2 (isospectrally). But dim(0-forms) = 12, dim(2-forms) = 20.
# That can't be fully isospectral. Poincare duality says:
# nonzero eigenvalues of Delta_p = nonzero eigenvalues of Delta_{d-p}
# So the nonzero spectra match, but the kernel dimensions differ.

# Actually for closed surface: b0=b2=1, b1=0
# Delta_0 has 12 eigenvalues: 1 zero + 11 nonzero
# Delta_2 has 20 eigenvalues: 1 zero + 19 nonzero
# Nonzero eigenvalues of Delta_0 = 11 values
# Nonzero eigenvalues of Delta_2 should CONTAIN those 11 values plus 8 more.
# No wait: Hodge decomposition for edges (1-forms):
# 30 = exact_1 + coexact_1 + harmonic_1
# exact_1 = rank(d_0) = 12-1 = 11
# harmonic_1 = b_1 = 0
# coexact_1 = 30 - 11 - 0 = 19
# Exact eigenvalues of Delta_1 = nonzero eigenvalues of Delta_0 (11 values)
# Coexact eigenvalues of Delta_1 = nonzero eigenvalues of Delta_2 (19 values)
# So Delta_2 has 20 eigenvalues: 1 zero + 19 nonzero.
# The 19 nonzero eigenvalues of Delta_2 = coexact eigenvalues of Delta_1.
# The 11 nonzero eigenvalues of Delta_0 = exact eigenvalues of Delta_1.

# Tr(Delta_0) = 60 (all eigenvalues, including zero)
# Tr(nonzero Delta_0) = 60 (since the one zero eigenvalue contributes nothing)
# Tr(Delta_1) = Tr(exact) + Tr(coexact) = Tr(nonzero Delta_0) + Tr(nonzero Delta_2)
#             = 60 + Tr(nonzero Delta_2)
# Tr(Delta_2) = Tr(nonzero Delta_2) + 0 = Tr(nonzero Delta_2)

# For the icosahedron, what is Tr(Delta_2)?
# Each face-face pair shares an edge if they are adjacent.
# The dual graph of icosahedron = dodecahedron (20 vertices, degree 3).
# Delta_2 on 2-forms = Laplacian of the dual graph = 20x20 with degree 3.
# Tr(Delta_2) = 20*3 = 60.

# So: c1_ico = Tr(D_ico^2) = Tr(Delta_0) + Tr(Delta_1) + Tr(Delta_2)
# = 60 + (60 + 60) + 60 = 240
# Wait: Tr(Delta_1) = Tr(exact) + Tr(coexact)
# = Tr(nonzero Delta_0) + Tr(nonzero Delta_2)
# = 60 + 60 = 120
# c1_ico = 60 + 120 + 60 = 240

# Actually we need to be careful: D^2 = bigoplus Delta_p.
# c1_ico = Tr(D_ico^2) = Tr(Delta_0) + Tr(Delta_1) + Tr(Delta_2)

# For Delta_1 on edges of icosahedron:
# Delta_1 = d_0^T*d_0 + d_1*d_1^T
# Tr(Delta_1) = Tr(d_0^T*d_0) + Tr(d_1*d_1^T)
# Tr(d_0^T*d_0) = Tr(d_0*d_0^T) = Tr(Delta_0) = 60
# Tr(d_1*d_1^T) = Tr(d_1^T*d_1) = Tr(Delta_2) = 60
# So Tr(Delta_1) = 60 + 60 = 120. Yes.

c1_ico = 60 + 120 + 60  # = 240
print(f"Icosahedron spectral data:")
print(f"  c0_ico = dim(H_ico) = {A1}")
print(f"  c1_ico = Tr(D_ico^2) = Tr(Delta_0)+Tr(Delta_1)+Tr(Delta_2)")
print(f"         = 60 + 120 + 60 = {c1_ico}")
print(f"  c1_ico / c0_ico = {c1_ico}/{A1} = {Fraction(c1_ico, A1)} = {c1_ico/A1:.10f}")
print()

# So the spectral density of the vertex figure: c1_ico/c0_ico = 240/62 = 120/31
# TARGET / (c1_ico/c0_ico) = 137.08 / (240/62) = 137.08 * 62/240 = 35.42
# Hmm, not clean.

# But wait: c1_ico = 240 = 2*N. And c0_ico = 62 = A1.
# c1_ico = 2*N is interesting: the number of E8 roots!
print(f"REMARKABLE: c1_ico = Tr(D_ico^2) = 240 = 2*N = |E8 roots|!")
print()

# This is a STRUCTURAL connection between the vertex figure (icosahedron)
# and the E8 root system. The trace of D^2 on the icosahedron equals
# the number of E8 roots = 2*|2I| = 240.

# Can we use this? TARGET = 4*a1*phi^4
# = (c1_ico/degree) * phi^4 = (240/12)*phi^4 = 20*phi^4. Same as before.
# Or: TARGET = c1_ico * phi^4 / degree
val_24 = c1_ico * phi**4 / degree
print(f"c1_ico * phi^4 / degree = {c1_ico}*{phi**4:.6f}/{degree} = {val_24:.10f}")
print(f"TARGET = {TARGET:.10f}")
print(f"MATCH: {abs(val_24 - TARGET) < 1e-8}")
print()

print("This rewrites as:")
print(f"  4*a1*phi^4 = Tr(D_ico^2) * phi^4 / degree")
print(f"             = (|E8 roots|) * phi^4 / (vertex degree)")
print(f"  where Tr(D_ico^2) = 240 = |E8 roots| = 2*N")
print()
print("  INTERPRETATION: The tree-level coupling 1/alpha_0 = 4*a1*phi^4")
print("  is the spectral trace of the VERTEX FIGURE Hodge-Dirac, amplified")
print("  by the golden ratio factor phi^4 and normalized by the degree.")
print()
print("  But: 240/12 = 20 = 4*a1, so this is STILL algebraically equivalent")
print("  to the known form. The phi^4 factor remains unexplained by spectral traces.")
print()

# =============================================================================
# SECTION 25: phi^4 FROM SPECTRAL GAP
# =============================================================================
print("=" * 76)
print("SECTION 25: phi^4 FROM SPECTRAL GAP RATIO")
print("=" * 76)
print()

# The ONLY spectral source of phi in the framework is:
# 1. Cayley adjacency eigenvalues (many contain phi)
# 2. Spectral gap of Cayley Laplacian: 12-6*phi = 6/phi^2
# 3. Spectral gap of Hopf fiber: 1/phi^2
# 4. Coexact gap: 7-4*phi
# 5. Lambda_max(D^2) = b1*phi^2

# The spectral gap of the Cayley Laplacian is: gap = 6*(2-phi) = 6/phi^2
# => phi^2 = 6/gap
# => phi^4 = 36/gap^2

# So: TARGET = 4*a1*phi^4 = 4*a1*36/gap^2 = 144*a1/gap^2
val_25a = 144 * a1 / cayley_gap**2
print(f"TARGET = 144*a1/gap_Cayley^2")
print(f"  = 144*{a1}/{cayley_gap:.6f}^2")
print(f"  = {val_25a:.10f}")
check_match("144*a1/gap^2", val_25a, "STRUCTURAL")
print()

# Or: 4*a1 * (degree/gap)^2 = TARGET  [since degree/gap = 12/(6/phi^2) = 2*phi^2]
# (degree/gap)^2 = 4*phi^4
# TARGET = 4*a1 * 4*phi^4 ... no, that gives 16*a1*phi^4.
# Wait: degree/gap = 12/(6/phi^2) = 12*phi^2/6 = 2*phi^2.
# (degree/gap)^2 = 4*phi^4. And TARGET = a1*(degree/gap)^2 = a1*4*phi^4.
# Correct! So:

# TARGET = a1 * (degree/gap_Cayley)^2

# The spectral expansion ratio R = degree/gap is a well-defined spectral
# invariant of the Cayley graph. R^2 = (expansion ratio)^2.
# And the coupling is a1 * R^2.

# But where does a1 come from in this product?
# a1 = 5 is the GALOIS NORM of the coexact spectral gap (N(7-4*phi) = 5).
# It's also the BASIC INPUT of the framework.

# Can we write a1 in terms of spectral data?
# a1 = N(gap_coexact) = (7-4*phi)*(7-4*phi') = 49-28*phi-28*phi'+16*phi*phi'
# = 49 - 28*(phi+phi') + 16*phi*phi' = 49 - 28*1 + 16*(-1) = 49-28-16 = 5
# Using phi+phi' = 1, phi*phi' = -1.

print("a1 from spectral data:")
print(f"  a1 = N(gap_coexact) = (7-4*phi)*(7-4*phi')")
print(f"     = 49 - 28*(phi+phi') + 16*phi*phi'")
print(f"     = 49 - 28*1 + 16*(-1) = {49-28-16}")
print()

print("So: TARGET = N(gap_coexact) * (degree/gap_Cayley)^2")
print()
val_25b = (lam1_coexact*(7-4*phi_c)) * (degree/cayley_gap)**2
print(f"  = {lam1_coexact*(7-4*phi_c):.6f} * {(degree/cayley_gap)**2:.6f}")
print(f"  = {val_25b:.10f}")
check_match("N(gap_coexact) * (degree/gap)^2", val_25b, "SPECTRAL-STRUCTURAL")
print()

print("ASSESSMENT: This formula expresses 4*a1*phi^4 ENTIRELY from spectral data:")
print(f"  - N(gap_coexact) = Galois norm of coexact spectral gap = {a1}")
print(f"  - degree = 12 (vertex degree of Cayley graph)")
print(f"  - gap_Cayley = 12-6*phi = smallest nonzero Cayley Laplacian eigenvalue")
print(f"  ALL three are spectral invariants of the 600-cell complex.")
print()
print("CATEGORY: [STRUCTURAL]")
print("  The formula IS purely spectral, but the combination")
print("  'Galois norm * (expansion ratio)^2' lacks a clear physical derivation")
print("  from a spectral action principle Tr(f(D/Lambda)).")
print("  It is an IDENTITY relating spectral quantities, not a DERIVATION")
print("  from a variational principle.")
print()

# =============================================================================
# SECTION 26: ATTEMPT - RATIO OF HEAT KERNEL DETERMINANTS
# =============================================================================
print("=" * 76)
print("SECTION 26: LOG-DETERMINANT ROUTE")
print("=" * 76)
print()

# In quantum field theory, the effective coupling at scale mu is:
# 1/alpha(mu) = 1/alpha_0 + (b_0/2*pi)*log(Lambda/mu)
# where b_0 is the one-loop beta function coefficient.
#
# For the DISCRETE spectral action, the one-loop correction involves:
# log det'(D^2) = sum_{lambda>0} log(lambda)
#
# We can try: 1/alpha_0 = log det'(D_M^2) / (something involving D_F)
#
# From the known spectral data (exp463):
# log det'(D_0^2-D_3^2) = log det'(D_1^2-D_2^2) (Poincare duality)
# Even = Odd = 2013.68

# Without computing the full spectrum, let's use the MOMENTS to estimate.
# For a spectrum with c_0 = number of eigenvalues, c_1 = sum, c_2 = sum of squares/2:
# A rough estimate: if eigenvalues were all equal to c_1/c_0:
# log det' ~ c_0 * log(c_1/c_0) [very rough]
log_det_rough = c0 * log(c1/c0)
print(f"Rough: log det' ~ c0*log(c1/c0) = {c0}*log({c1/c0:.4f}) = {log_det_rough:.4f}")
print(f"  Compare with exp463: Even=Odd=2013.68")
print()

# The "spectral conductance":
# sigma = c1/c0 - c2/c1 = mean(D^2) - mean(D^4)/mean(D^2)
sigma = c1/c0 - c2/c1
print(f"Spectral variance proxy: c1/c0 - c2/c1 = {c1/c0:.4f} - {c2/c1:.4f} = {sigma:.6f}")
print(f"  Compare TARGET = {TARGET:.6f}")
print(f"  Ratio: {sigma/TARGET:.6f}")
print()

# =============================================================================
# SECTION 27: FINAL VERDICT
# =============================================================================
print()
print("=" * 76)
print("=" * 76)
print("FINAL VERDICT")
print("=" * 76)
print("=" * 76)
print()

print("QUESTION: Can 4*a1*phi^4 be DERIVED from the spectral action?")
print()
print("ANSWER: PARTIALLY.")
print()
print("We found THREE spectral formulas for 4*a1*phi^4:")
print()
print("  (1) 4*a1*phi^4 = N(gap_coexact) * (degree/gap_Cayley)^2  [SPECTRAL]")
print(f"      = {a1} * ({degree}/{cayley_gap:.4f})^2 = {a1*(degree/cayley_gap)**2:.6f}")
print()
print("  (2) 4*a1*phi^4 = L(3)*L(3')*(L(3')/L(3))^2              [VERTEX FIGURE]")
print(f"      = {L3*L3p:.4f} * {(L3p/L3)**2:.4f} = {L3*L3p*(L3p/L3)**2:.6f}")
print()
print("  (3) 4*a1*phi^4 = Tr(D_ico^2) * phi^4 / degree            [MIXED]")
print(f"      = {c1_ico} * {phi**4:.4f} / {degree} = {c1_ico*phi**4/degree:.6f}")
print()
print("CLASSIFICATION:")
print()
print("  (1) [STRUCTURAL] All ingredients are spectral invariants of the 600-cell.")
print("      But the specific combination 'norm * ratio^2' is not derived from")
print("      a Lagrangian or variational principle. It's a spectral IDENTITY.")
print()
print("  (2) [DERIVED - vertex figure level] Everything from icosahedral spectrum.")
print("      The product L(3)*L(3') = a1(a1-1) = 4*a1 uses a1=5 specifically.")
print("      The ratio L(3')/L(3) = phi^2 is an EXACT algebraic identity.")
print("      This is the CLEANEST form: the alpha equation coefficients are")
print("      entirely determined by the Galois-conjugate pair of vertex-figure")
print("      eigenvalues (L(3), L(3')).")
print()
print("  (3) [MIXED] Uses both vertex figure spectral data (via c1_ico=240) and")
print("      the golden ratio phi^4, which enters as a spectral amplification.")
print("      The phi^4 factor is NOT from a spectral trace.")
print()
print("KEY NEGATIVE RESULT:")
print("  The Seeley-DeWitt coefficients (c0=2640, c1=14880, c2=55920) of the")
print("  Hodge-Dirac operator on the full 600-cell are INTEGER multiples of 240.")
print("  They do NOT contain phi (or sqrt(5)) in any ratio or combination.")
print("  Therefore, 4*a1*phi^4 CANNOT be extracted from c0, c1, c2 alone.")
print()
print("  The irrational part (phi^4) enters ONLY through the eigenvalue structure")
print("  of the Cayley graph or vertex figure, which is more refined than the")
print("  spectral action moments (which are sums over ALL eigenvalues).")
print()
print("WHAT REMAINS OPEN:")
print("  - A variational derivation of WHY the combination")
print("    'N(gap)*R^2' determines the coupling.")
print("  - Whether the full eigenvalue spectrum (not just moments) of the")
print("    Hodge-Dirac operator, combined with the McKay finite Dirac,")
print("    gives a spectral-action derivation of the LINEAR coefficient.")
print("  - The spectral zeta function route: zeta_D(s) for specific s values")
print("    might encode phi through analytic continuation.")
print()
print("HONEST CONCLUSION per REGULA ZERO:")
print("  4*a1*phi^4 is STRUCTURAL -- expressed entirely in terms of spectral")
print("  data of the 600-cell -- but NOT DERIVED from a spectral action")
print("  principle (i.e., not from Tr(f(D/Lambda)) for some cutoff function f).")
print("  The gap between STRUCTURAL and DERIVED remains open.")
print()

# =============================================================================
# NEW FINDING: c1_ico = 240 = |E8 roots|
# =============================================================================
print("=" * 76)
print("NEW FINDING: Tr(D_ico^2) = 240 = |E8 roots|")
print("=" * 76)
print()
print("The trace of D^2 on the icosahedral simplicial complex equals 240,")
print("which is the number of E8 roots. This is:")
print("  - 12 vertices * degree 5 = 60 for 0-forms")
print("  - 60 + 60 = 120 for 1-forms (exact + coexact)")
print("  - 20 faces * degree 3 = 60 for 2-forms")
print("  Total: 60 + 120 + 60 = 240 = 2*N")
print()
print("  Combined with: dim(H_ico) = 62 = A1 (Seeley-DeWitt reduced)")
print("  we get: c1_ico/c0_ico = 240/62 = 120/31")
print()
print("  CATEGORY: [STRUCTURAL] - Interesting numerical coincidence that")
print("  connects vertex figure spectral data to E8 root count, but the")
print("  identity 2*N = 2*a1! = 240 = Tr(D_ico^2) needs deeper understanding.")
print("  It follows from: Tr(Delta_p) = (number of p-simplices)*(degree of p-simplex)")
print("  For the icosahedron: 12*5 + (12*5+20*3) + 20*3 = 60+120+60 = 240.")
print("  And 12*5 = 60 = N/2 (not obviously 'E8').")
print()

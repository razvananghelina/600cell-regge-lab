"""
EXP-155: PPN Parameters for the STG Metric
===========================================

QUESTION: What are the PPN (Parameterized Post-Newtonian) parameters
for the STG metric with alpha=1/2?

A Gemini reviewer incorrectly claimed beta_PPN = 1.5 using alpha=1 (wrong value).
We need to compute the correct PPN parameters with alpha=1/2.

METRIC:
  STG: ds^2 = -dt^2/C^2 + C^2 dr^2 + C^(2*alpha) r^2 dOmega^2
  C(r) = 1 + M/r,  alpha = 1/2

APPROACH:
1. Expand metric components in powers of U = M/r
2. Transform to isotropic coordinates (spatial part conformally flat)
3. Read off PPN parameters gamma and beta
4. Compare with Schwarzschild (gamma=1, beta=1)
5. Check consistency with known 1PN results (precession = GR, deflection = GR)

DERIVAT: All results from symbolic computation of the metric.
"""

import sympy as sp
from sympy import symbols, Rational, sqrt, series, collect, simplify, expand
from sympy import cos, sin, pi, oo, solve, Function, Derivative, integrate
from sympy import O as BigO

# ================================================================
# SETUP
# ================================================================

print("=" * 70)
print("EXP-155: PPN PARAMETERS FOR STG METRIC")
print("=" * 70)

r, M_sym, U = symbols('r M U', positive=True)
alpha_sym = symbols('alpha')
eps = symbols('epsilon')  # bookkeeping parameter

# C(r) = 1 + M/r
C = 1 + M_sym / r

print("""
STG metric (general alpha):
  ds^2 = -dt^2/C^2 + C^2 dr^2 + C^(2*alpha) r^2 dOmega^2
  C(r) = 1 + M/r

Standard PPN isotropic form:
  ds^2 = -(1 - 2U + 2*beta*U^2 + ...) dt^2
       + (1 + 2*gamma*U + ...) (dr_iso^2 + r_iso^2 dOmega^2)
  where U = M/r_iso (isotropic Newtonian potential)
""")

# ================================================================
# PART 1: Direct expansion in STG coordinates (NOT isotropic)
# ================================================================

print("=" * 70)
print("PART 1: METRIC EXPANSION IN STG COORDINATES")
print("=" * 70)

# Use U = M/r as expansion parameter
# C = 1 + U, so C^n = (1+U)^n

# g_tt = -1/C^2 = -(1+U)^(-2)
# g_rr = C^2 = (1+U)^2
# g_theta = C^(2*alpha) * r^2

print("\nExpanding metric components in U = M/r:")
print("-" * 50)

# Expand g_tt
g_tt_series = series((1 + U)**(-2), U, 0, n=5)
print(f"\n  g_tt = -1/C^2 = -(1+U)^(-2)")
print(f"       = -({g_tt_series})")
print(f"       = -(1 - 2U + 3U^2 - 4U^3 + 5U^4 + ...)")

# Expand g_rr
g_rr_series = series((1 + U)**2, U, 0, n=5)
print(f"\n  g_rr = C^2 = (1+U)^2")
print(f"       = {g_rr_series}")
print(f"       = 1 + 2U + U^2")

# Expand g_angular for general alpha
print(f"\n  g_theta_theta / r^2 = C^(2*alpha) = (1+U)^(2*alpha)")
for a_val in [Rational(1,2), 1, 0]:
    name = {Rational(1,2): "1/2", 1: "1", 0: "0"}[a_val]
    g_ang = series((1 + U)**(2*a_val), U, 0, n=5)
    print(f"    alpha={name}: {g_ang}")

print("""
WARNING: These are NOT in isotropic form!
The STG metric has g_rr != g_theta/r^2 (unless alpha=1).
We must transform to isotropic coordinates to extract PPN parameters.
""")

# ================================================================
# PART 2: TRANSFORMATION TO ISOTROPIC COORDINATES
# ================================================================

print("=" * 70)
print("PART 2: ISOTROPIC COORDINATE TRANSFORMATION")
print("=" * 70)

print("""
The STG metric in STG coordinates:
  ds^2 = -(1+M/r)^(-2) dt^2 + (1+M/r)^2 dr^2 + (1+M/r)^(2a) r^2 dOmega^2

For isotropic form, we need:
  ds^2 = -A(R) dt^2 + B(R)(dR^2 + R^2 dOmega^2)

This requires: B(R) dR^2 = g_rr dr^2 and B(R) R^2 = g_ang r^2

So: B(R) = g_ang r^2 / R^2 = (1+M/r)^(2a) r^2/R^2
And: B(R) dR^2 = (1+M/r)^2 dr^2

Dividing: (dR/R)^2 = [(1+M/r)^2 / ((1+M/r)^(2a) r^2/R^2)] * ...

Let me do this more carefully.
From B*R^2 = C^(2a) r^2:  R = r * C^a / sqrt(B)... circular.

Better approach: set R*sqrt(B) = r*C^a and sqrt(B)*dR = C*dr.
Then: d(R*sqrt(B))/dR * ... this gets messy. Let me use the standard method.
""")

print("--- Standard method: requiring spatial conformal flatness ---\n")
print("We need R(r) such that:")
print("  C^2 dr^2 + C^(2a) r^2 dOmega^2 = f(R)^2(dR^2 + R^2 dOmega^2)")
print()
print("This gives two conditions:")
print("  (1) f(R)^2 (dR/dr)^2 = C^2")
print("  (2) f(R)^2 R^2 = C^(2a) r^2")
print()
print("From (2): f(R) = C^a * r / R")
print("Substituting into (1): C^(2a) r^2/R^2 * (dR/dr)^2 = C^2")
print("=> (dR/dr)^2 = C^(2-2a) * R^2/r^2")
print("=> dR/R = +/- C^(1-a) dr/r")
print("=> d(ln R) = (1-a) * d(ln C) + d(ln r)")
print("   [since d(ln C) = d(ln(1+M/r)) = -M/r^2 / (1+M/r) dr]")
print()

# Symbolic integration
print("Integrating: ln R = ln r + (1-a) * ln C + const")
print("=> R = const * r * C^(1-a) = const * r * (1+M/r)^(1-a)")
print()

# Let me verify this formally
R_sym = symbols('R', positive=True)
a = symbols('a', positive=True)

print("Let's verify: if R = r * (1+M/r)^(1-a), then")
print("  dR/dr = (1+M/r)^(1-a) + r*(1-a)*(1+M/r)^(-a) * (-M/r^2)")
print("        = (1+M/r)^(1-a) * [1 - (1-a)*M / (r(1+M/r))]")
print("        = (1+M/r)^(1-a) * [1 - (1-a)*M / (r+M)]")
print("        = (1+M/r)^(1-a) * [(r+M - (1-a)M) / (r+M)]")
print("        = (1+M/r)^(1-a) * [(r + aM) / (r+M)]")
print()

# Now compute f^2 R^2 and f^2 (dR/dr)^2
print("Conformal factor f: f = C^a * r / R = C^a * r / (r * C^(1-a)) = C^(2a-1)")
print()
print("Check condition (2): f^2 R^2 = C^(4a-2) * r^2 * C^(2-2a) = C^(2a) * r^2  [CORRECT]")
print()
print("Check condition (1): f^2 (dR/dr)^2")
print("  = C^(4a-2) * C^(2-2a) * [(r+aM)/(r+M)]^2")
print("  = C^(2a) * [(r+aM)/(r+M)]^2")
print("  Need this = C^2 = (1+M/r)^2")
print("  So need: C^(2a) * [(r+aM)/(r+M)]^2 = C^2")
print("  => [(r+aM)/(r+M)]^2 = C^(2-2a) = (1+M/r)^(2-2a)")
print("  => (r+aM)^2 / (r+M)^2 = [(r+M)/r]^(2-2a)")
print()
print("This is NOT generally true for arbitrary a!")
print("The spatial part is NOT conformally flat in general.")
print()

# ================================================================
# PART 3: PERTURBATIVE ISOTROPIC TRANSFORMATION
# ================================================================

print("=" * 70)
print("PART 3: PERTURBATIVE APPROACH (SERIES IN M/r)")
print("=" * 70)

print("""
Since the exact transformation doesn't give conformal flatness,
we work perturbatively in U = M/r.

Write r = R + sum_n c_n * M^n / R^(n-1)
     or R = r + sum_n d_n * M^n / r^(n-1)

and match the metric to isotropic form order by order.
""")

# Work with symbolic series
# Let u = M/r be small. We want R(r) = r * h(u) where h(0) = 1
#
# R = r*(1 + h1*u + h2*u^2 + h3*u^3 + ...)
#
# Then the spatial metric:
#   C^2 dr^2 + C^(2a) r^2 dOmega^2
# must become f(R)^2 (dR^2 + R^2 dOmega^2)

u = symbols('u')
h1, h2, h3, h4 = symbols('h1 h2 h3 h4')
a_sym = symbols('a')

# R = r*(1 + h1*u + h2*u^2 + h3*u^3 + ...)
# Since u = M/r, and R ~ r for r >> M, this is correct

# The angular part gives: f^2 * R^2 = (1+u)^(2a) * r^2
# f^2 = (1+u)^(2a) * r^2 / R^2 = (1+u)^(2a) / (1+h1*u+h2*u^2+...)^2

# The radial part: f^2 * (dR/dr)^2 = (1+u)^2
# dR/dr: R = r*(1 + h1*M/r + h2*M^2/r^2 + ...)
#           = r + h1*M + h2*M^2/r + ...
# dR/dr = 1 - h2*M^2/r^2 - 2*h3*M^3/r^3 + ...
#       = 1 - h2*u^2 - 2*h3*u^3 + ...  (since M/r = u, dM/dr=0)

# Wait, let me be more careful. u = M/r => du/dr = -M/r^2 = -u/r

# R = r * H(u) where H = 1 + h1*u + h2*u^2 + ...
# dR/dr = H + r * H' * du/dr = H + r * H' * (-u/r) = H - u * H'
# where H' = dH/du

H = 1 + h1*u + h2*u**2 + h3*u**3 + h4*u**4
Hp = sp.diff(H, u)
dR_dr = H - u * Hp

print("Setting R = r * H(u) where H = 1 + h1*u + h2*u^2 + h3*u^3 + ...")
print(f"  dR/dr = H - u*H' = {sp.expand(dR_dr)}")
print()

# Angular condition: f^2 = C^(2a) / H^2 where C = 1+u
# Radial condition: f^2 * (dR/dr)^2 = C^2
# => C^(2a)/H^2 * (H - u*H')^2 = C^2
# => (H - u*H')^2 / H^2 = C^(2-2a) = (1+u)^(2-2a)

print("Consistency condition:")
print("  (H - u*H')^2 / H^2 = (1+u)^(2-2a)")
print()

# Now expand both sides in u, for general a
# LHS:
LHS = (dR_dr / H)**2

# Expand LHS in u to order 4
LHS_expanded = sp.series(LHS, u, 0, n=5).removeO()
print("LHS = (H - u*H')^2 / H^2, expanded:")
LHS_coeffs = {}
for n in range(5):
    coeff = LHS_expanded.coeff(u, n)
    LHS_coeffs[n] = sp.simplify(coeff)
    print(f"  u^{n}: {LHS_coeffs[n]}")

# RHS: (1+u)^(2-2a) expanded
print()
print("RHS = (1+u)^(2-2a), expanded:")
# Use binomial series: (1+u)^p = sum_{n=0}^inf binom(p,n) u^n
p = 2 - 2*a_sym
RHS_coeffs = {}
for n in range(5):
    binom_coeff = sp.Rational(1,1)
    for k in range(n):
        binom_coeff *= (p - k) / (k + 1)
    RHS_coeffs[n] = sp.expand(binom_coeff)
    print(f"  u^{n}: {RHS_coeffs[n]}")

# Match coefficients
print("\n" + "=" * 70)
print("MATCHING COEFFICIENTS (general alpha = a)")
print("=" * 70)

# u^0: 1 = 1 (trivial)
print("\nu^0: 1 = 1  [trivial, OK]")

# u^1: -2*h1 = 2*(1-a) => h1 = -(1-a) = a-1
eq1 = sp.Eq(LHS_coeffs[1], RHS_coeffs[1])
sol1 = solve(eq1, h1)
print(f"\nu^1: {LHS_coeffs[1]} = {RHS_coeffs[1]}")
print(f"  => h1 = {sol1}")

h1_val = sol1[0]

# u^2: substitute h1
LHS2_sub = LHS_coeffs[2].subs(h1, h1_val)
LHS2_sub = sp.expand(LHS2_sub)
eq2 = sp.Eq(LHS2_sub, RHS_coeffs[2])
sol2 = solve(eq2, h2)
print(f"\nu^2: {LHS2_sub} = {RHS_coeffs[2]}")
print(f"  => h2 = {sol2}")

h2_val = sol2[0]

# u^3: substitute h1, h2
LHS3_sub = LHS_coeffs[3].subs(h1, h1_val).subs(h2, h2_val)
LHS3_sub = sp.expand(LHS3_sub)
eq3 = sp.Eq(LHS3_sub, RHS_coeffs[3])
sol3 = solve(eq3, h3)
print(f"\nu^3: {sp.simplify(LHS3_sub)} = {RHS_coeffs[3]}")
print(f"  => h3 = {sp.simplify(sol3[0])}")

h3_val = sol3[0]

# u^4: substitute h1, h2, h3
LHS4_sub = LHS_coeffs[4].subs(h1, h1_val).subs(h2, h2_val).subs(h3, h3_val)
LHS4_sub = sp.expand(LHS4_sub)
eq4 = sp.Eq(LHS4_sub, RHS_coeffs[4])
sol4 = solve(eq4, h4)
print(f"\nu^4: h4 = {sp.simplify(sol4[0])}")

h4_val = sol4[0]

# Now for alpha = 1/2
print("\n" + "=" * 70)
print("SPECIALIZING TO alpha = 1/2")
print("=" * 70)

h1_half = h1_val.subs(a_sym, Rational(1,2))
h2_half = h2_val.subs(a_sym, Rational(1,2))
h3_half = sp.simplify(h3_val.subs(a_sym, Rational(1,2)))
h4_half = sp.simplify(h4_val.subs(a_sym, Rational(1,2)))

print(f"\n  alpha = 1/2:")
print(f"  h1 = {h1_half}")
print(f"  h2 = {h2_half}")
print(f"  h3 = {h3_half}")
print(f"  h4 = {h4_half}")
print(f"\n  R = r * (1 + ({h1_half})*M/r + ({h2_half})*(M/r)^2 + ({h3_half})*(M/r)^3 + ...)")

# For comparison, Schwarzschild isotropic:
# r_Schw = R_iso * (1 + M/(2*R_iso))^2
# => R_iso = r_Schw * (1 - M/r + 3M^2/(4r^2) - ...)  [inversion of above]
# h1_GR = -1/2, h2_GR = 0, h3_GR = 1/16, h4_GR = 0
# (standard result for Schwarzschild isotropic)

# Actually let me not assume. Let's compute the Schwarzschild case too.
# Schwarzschild: ds^2 = -(1-2M/r)dt^2 + dr^2/(1-2M/r) + r^2 dOmega^2
# Isotropic: ds^2 = -((1-M/(2R))/(1+M/(2R)))^2 dt^2 + (1+M/(2R))^4 (dR^2+R^2 dOmega^2)
# Relation: r = R*(1+M/(2R))^2

print("\n--- For comparison, Schwarzschild isotropic form: ---")
print("  r = R*(1 + M/(2R))^2 = R + M + M^2/(4R)")
print("  => R = r*(1 - M/r + 3M^2/(4r^2) - M^3/(2r^3) + ...)")
print("  => h1_GR = -1, h2_GR = 3/4, h3_GR = -1/2, h4_GR = 5/16")
print()

# Let me verify the Schwarzschild values by solving the same way
# Schwarzschild: g_rr = 1/(1-2u), g_ang = r^2  [alpha-equivalent: special case]
# This is different - Schwarzschild is NOT of the form C^2 dr^2 + C^(2a) r^2 dOmega^2
# So we do this separately

# ================================================================
# PART 4: COMPUTE PPN PARAMETERS
# ================================================================

print("=" * 70)
print("PART 4: PPN PARAMETERS FROM ISOTROPIC METRIC")
print("=" * 70)

print("""
After the transformation R = r*H(u), the metric becomes:
  ds^2 = -A(R) dt^2 + B(R)(dR^2 + R^2 dOmega^2)

where A(R) = 1/C^2 and B(R) = C^(2a)*r^2/R^2.

We need to express everything in terms of W = M/R (isotropic potential).
""")

# We have R = r*H(u) where u = M/r => r = M/u, R = M*H(u)/u
# So W = M/R = u/H(u)
# Inversion: u = W + ... (perturbative)

# Let's express u in terms of W perturbatively
W = symbols('W')

# W = u/H(u) = u/(1 + h1*u + h2*u^2 + ...)
# => u = W * H(u) = W * (1 + h1*u + h2*u^2 + ...)
# Iterating: u = W + h1*W^2 + (h1^2 + h2)*W^3 + ...

# Do this for general alpha
u_of_W = W + h1_val*W**2 + (h1_val**2 + h2_val)*W**3
u_of_W = sp.expand(u_of_W)

# Add one more term for 4th order
u_of_W_4 = W + h1_val*W**2 + (h1_val**2 + h2_val)*W**3 + \
            (h1_val**3 + 2*h1_val*h2_val + h3_val)*W**4

print("Perturbative inversion u(W):")
print(f"  u = {sp.expand(u_of_W_4)}")
print()

# Now compute A(W) = 1/(1+u)^2 and B(W) = (1+u)^(2a) / H(u)^2
# using u = u(W)

# A(W) = (1+u)^(-2) as series in W
# B(W) = (1+u)^(2a) / H(u)^2

# Let's do this symbolically for alpha = 1/2

a_val = Rational(1,2)
h1v = h1_val.subs(a_sym, a_val)
h2v = h2_val.subs(a_sym, a_val)
h3v = sp.simplify(h3_val.subs(a_sym, a_val))
h4v = sp.simplify(h4_val.subs(a_sym, a_val))

print(f"For alpha = 1/2: h1={h1v}, h2={h2v}, h3={h3v}, h4={h4v}")

# u(W) for alpha=1/2
u_W = W + h1v*W**2 + (h1v**2 + h2v)*W**3 + (h1v**3 + 2*h1v*h2v + h3v)*W**4
u_W = sp.expand(u_W)
print(f"\n  u(W) = {u_W}")

# A(W) = 1/(1+u)^2
A_W = 1 / (1 + u_W)**2
A_W_series = sp.series(A_W, W, 0, n=5).removeO()
A_W_series = sp.expand(A_W_series)
print(f"\n  A(W) = -g_tt = 1/(1+u)^2 = {A_W_series}")

# Extract PPN coefficients from g_tt = -(1 - 2W + 2*beta*W^2 + ...)
# So A(W) = 1 - 2W + 2*beta*W^2 - ...

a0 = A_W_series.coeff(W, 0)
a1 = A_W_series.coeff(W, 1)
a2 = A_W_series.coeff(W, 2)
a3 = A_W_series.coeff(W, 3)
a4 = A_W_series.coeff(W, 4)

print(f"\n  g_tt = -(a0 + a1*W + a2*W^2 + a3*W^3 + a4*W^4)")
print(f"  a0 = {a0}")
print(f"  a1 = {a1}")
print(f"  a2 = {a2}")
print(f"  a3 = {a3}")
print(f"  a4 = {a4}")

print(f"\n  Standard PPN form: g_tt = -(1 - 2W + 2*beta*W^2 + ...)")
print(f"  Comparing:")
print(f"    Coefficient of W^0: {a0} = 1  [OK]")
print(f"    Coefficient of W^1: {a1} = -2  [check: {a1 == -2}]")
beta_PPN = a2 / 2
print(f"    Coefficient of W^2: {a2} = 2*beta => beta = {beta_PPN}")

# B(W) = conformal factor for spatial part
B_W = (1 + u_W)**(2*a_val) / (1 + h1v*u_W + h2v*u_W**2 + h3v*u_W**3)**2
# More carefully: B = C^(2a)*r^2/R^2 = (1+u)^(2a) / H^2
# H = 1 + h1*u + h2*u^2 + ...

H_of_u = 1 + h1v*(u_W) + h2v*(u_W)**2 + h3v*(u_W)**3
B_W_expr = (1 + u_W)**(2*a_val) / H_of_u**2
B_W_series = sp.series(B_W_expr, W, 0, n=5).removeO()
B_W_series = sp.expand(B_W_series)

print(f"\n  B(W) = C^(2a)/H^2 = {B_W_series}")

b0 = B_W_series.coeff(W, 0)
b1 = B_W_series.coeff(W, 1)
b2 = B_W_series.coeff(W, 2)
b3 = B_W_series.coeff(W, 3)

print(f"\n  g_rr (isotropic) = B(W):")
print(f"    Coefficient of W^0: {b0}")
print(f"    Coefficient of W^1: {b1}")
print(f"    Coefficient of W^2: {b2}")
print(f"    Coefficient of W^3: {b3}")

print(f"\n  Standard PPN form: g_rr = 1 + 2*gamma*W + ...")
gamma_PPN = b1 / 2
print(f"    => gamma = {b1}/2 = {gamma_PPN}")

# ================================================================
# PART 5: ALSO DO IT FOR GENERAL alpha
# ================================================================

print("\n" + "=" * 70)
print("PART 5: PPN PARAMETERS FOR GENERAL alpha")
print("=" * 70)

# Redo with general alpha
u_W_gen = W + h1_val*W**2 + (h1_val**2 + h2_val)*W**3 + \
          (h1_val**3 + 2*h1_val*h2_val + h3_val)*W**4
u_W_gen = sp.expand(u_W_gen)

# A(W) = 1/(1+u)^2
A_gen = 1 / (1 + u_W_gen)**2
A_gen_series = sp.series(A_gen, W, 0, n=4).removeO()
A_gen_series = sp.expand(A_gen_series)

a1_gen = A_gen_series.coeff(W, 1)
a2_gen = A_gen_series.coeff(W, 2)
a3_gen = A_gen_series.coeff(W, 3)

beta_gen = sp.simplify(a2_gen / 2)
print(f"\n  g_tt coefficient of W^1: {sp.simplify(a1_gen)} [should be -2]")
print(f"  g_tt coefficient of W^2: {sp.simplify(a2_gen)}")
print(f"  beta(alpha) = {beta_gen}")

# Evaluate beta for specific alpha values
print("\n  beta for specific alpha values:")
for av in [0, Rational(1,4), Rational(1,3), Rational(1,2), Rational(2,3), Rational(3,4), 1]:
    beta_val = beta_gen.subs(a_sym, av)
    print(f"    alpha = {av}: beta = {sp.simplify(beta_val)}")

# B(W) for general alpha
H_gen = 1 + h1_val*(u_W_gen) + h2_val*(u_W_gen)**2
B_gen = (1 + u_W_gen)**(2*a_sym) / H_gen**2
B_gen_series = sp.series(B_gen, W, 0, n=3).removeO()
B_gen_series = sp.expand(B_gen_series)

b1_gen = B_gen_series.coeff(W, 1)
gamma_gen = sp.simplify(b1_gen / 2)

print(f"\n  g_rr coefficient of W^1: {sp.simplify(b1_gen)}")
print(f"  gamma(alpha) = {gamma_gen}")

print("\n  gamma for specific alpha values:")
for av in [0, Rational(1,4), Rational(1,3), Rational(1,2), Rational(2,3), Rational(3,4), 1]:
    gamma_val = gamma_gen.subs(a_sym, av)
    print(f"    alpha = {av}: gamma = {sp.simplify(gamma_val)}")

# ================================================================
# PART 6: CONSISTENCY CHECKS
# ================================================================

print("\n" + "=" * 70)
print("PART 6: CONSISTENCY CHECKS")
print("=" * 70)

print("""
Known results from exp103-104 (alpha=1/2):
  - Light deflection = 4M/b = GR prediction
    PPN formula: deflection = (1+gamma)/2 * 4M/b
    For gamma=1: deflection = 4M/b [CONSISTENT]

  - Perihelion precession = 6*pi*M/r0 = GR prediction (at 1PN)
    PPN formula: precession = (2+2*gamma-beta)/3 * 6*pi*M/r0
    For gamma=1, beta=1: factor = (2+2-1)/3 = 1 [CONSISTENT]
""")

gamma_half = gamma_gen.subs(a_sym, Rational(1,2))
beta_half = beta_gen.subs(a_sym, Rational(1,2))

print(f"  With our computed values: gamma = {gamma_half}, beta = {beta_half}")

defl_factor = (1 + gamma_half) / 2
prec_factor = (2 + 2*gamma_half - beta_half) / 3

print(f"\n  Light deflection factor: (1+gamma)/2 = {defl_factor}")
print(f"  Precession factor: (2+2*gamma-beta)/3 = {sp.simplify(prec_factor)}")

# Check Nordtvedt effect: eta = 4*beta - gamma - 3
nordtvedt = 4*beta_half - gamma_half - 3
print(f"\n  Nordtvedt parameter: eta = 4*beta - gamma - 3 = {sp.simplify(nordtvedt)}")
print(f"  (GR has eta = 0)")

# ================================================================
# PART 7: EXPLICIT COMPARISON WITH SCHWARZSCHILD
# ================================================================

print("\n" + "=" * 70)
print("PART 7: EXPLICIT COMPARISON WITH SCHWARZSCHILD")
print("=" * 70)

print("""
Schwarzschild in isotropic coordinates:
  r = R_iso * (1 + M/(2*R_iso))^2

  g_tt = -((1-M/(2R))/(1+M/(2R)))^2
       = -(1 - 2W_iso + 2W_iso^2 - 2W_iso^3 + ...) where W_iso = M/(2R)

  Wait - standard convention uses U = M/R (not M/(2R)).
  Let me be careful with factors of 2.
""")

# Schwarzschild isotropic with U = M/R_iso (not M/(2R))
# g_tt = -((1 - U/2)/(1 + U/2))^2
# g_rr,iso = (1 + U/2)^4 * (dR/R)...no,
# In isotropic form: ds^2 = -((1-U/2)/(1+U/2))^2 dt^2 + (1+U/2)^4 (dR^2 + R^2 dOmega^2)
# where U = M/R

# Let me use W = M/R as before
# g_tt_GR = ((1 - W/2)/(1 + W/2))^2
g_tt_GR = ((1 - W/2)/(1 + W/2))**2
g_tt_GR_series = sp.series(g_tt_GR, W, 0, n=5).removeO()
g_tt_GR_series = sp.expand(g_tt_GR_series)

g_rr_GR = (1 + W/2)**4
g_rr_GR_series = sp.series(g_rr_GR, W, 0, n=5).removeO()
g_rr_GR_series = sp.expand(g_rr_GR_series)

print("Schwarzschild isotropic (W = M/R):")
print(f"  g_tt_GR = ((1-W/2)/(1+W/2))^2 = {g_tt_GR_series}")
print(f"  g_rr_GR = (1+W/2)^4 = {g_rr_GR_series}")
print()

a1_GR = g_tt_GR_series.coeff(W, 1)
a2_GR = g_tt_GR_series.coeff(W, 2)
b1_GR = g_rr_GR_series.coeff(W, 1)

print(f"  GR: g_tt coeff W^1 = {a1_GR}, W^2 = {a2_GR}")
print(f"  GR: g_rr coeff W^1 = {b1_GR}")
print(f"  GR: beta = {a2_GR}/2 = {a2_GR/2}, gamma = {b1_GR}/2 = {b1_GR/2}")

print()
print("STG (alpha=1/2) isotropic:")
print(f"  g_tt_STG = {A_W_series}")
print(f"  g_rr_STG = {B_W_series}")
print()

# Side-by-side coefficient comparison
print("  Coefficient comparison (W = M/R):")
print(f"  {'Order':>6} {'g_tt(GR)':>12} {'g_tt(STG)':>12} {'g_rr(GR)':>12} {'g_rr(STG)':>12}")
for n in range(5):
    gtt_gr = g_tt_GR_series.coeff(W, n)
    gtt_stg = A_W_series.coeff(W, n)
    grr_gr = g_rr_GR_series.coeff(W, n)
    grr_stg = B_W_series.coeff(W, n)
    print(f"  W^{n:>3}: {str(gtt_gr):>12} {str(gtt_stg):>12} {str(grr_gr):>12} {str(grr_stg):>12}")

# ================================================================
# PART 8: ADDRESS THE GEMINI CLAIM
# ================================================================

print("\n" + "=" * 70)
print("PART 8: ADDRESSING THE GEMINI beta=1.5 CLAIM")
print("=" * 70)

print("""
The Gemini reviewer computed PPN beta = 1.5 for the STG metric.
This appears to be based on using alpha=1 (conformal metric) instead
of the correct alpha=1/2.

Let's verify what happens with alpha=1:
""")

beta_alpha1 = beta_gen.subs(a_sym, 1)
gamma_alpha1 = gamma_gen.subs(a_sym, 1)
print(f"  alpha=1: beta = {sp.simplify(beta_alpha1)}, gamma = {sp.simplify(gamma_alpha1)}")

prec_alpha1 = (2 + 2*gamma_alpha1 - beta_alpha1) / 3
defl_alpha1 = (1 + gamma_alpha1) / 2
print(f"  alpha=1: precession factor = {sp.simplify(prec_alpha1)} [GR=1]")
print(f"  alpha=1: deflection factor = {sp.simplify(defl_alpha1)} [GR=1]")

beta_alpha_half = beta_gen.subs(a_sym, Rational(1,2))
gamma_alpha_half = gamma_gen.subs(a_sym, Rational(1,2))
print(f"\n  alpha=1/2: beta = {sp.simplify(beta_alpha_half)}, gamma = {sp.simplify(gamma_alpha_half)}")

prec_alpha_half = (2 + 2*gamma_alpha_half - beta_alpha_half) / 3
defl_alpha_half = (1 + gamma_alpha_half) / 2
print(f"  alpha=1/2: precession factor = {sp.simplify(prec_alpha_half)} [GR=1]")
print(f"  alpha=1/2: deflection factor = {sp.simplify(defl_alpha_half)} [GR=1]")

# ================================================================
# PART 9: 2PN DEVIATIONS
# ================================================================

print("\n" + "=" * 70)
print("PART 9: AT WHAT ORDER DOES STG DEVIATE FROM GR?")
print("=" * 70)

print("\nComparing g_tt coefficients (STG alpha=1/2 vs Schwarzschild):")
print("  (Both in isotropic coordinates with W = M/R)")
print()

any_diff = False
for n in range(5):
    gtt_gr = g_tt_GR_series.coeff(W, n)
    gtt_stg = A_W_series.coeff(W, n)
    diff = sp.simplify(gtt_stg - gtt_gr)
    match_str = "MATCH" if diff == 0 else f"DIFFER by {diff}"
    if diff != 0:
        any_diff = True
    print(f"  W^{n}: GR={gtt_gr}, STG={gtt_stg}  => {match_str}")

print("\nComparing g_rr coefficients:")
for n in range(5):
    grr_gr = g_rr_GR_series.coeff(W, n)
    grr_stg = B_W_series.coeff(W, n)
    diff = sp.simplify(grr_stg - grr_gr)
    match_str = "MATCH" if diff == 0 else f"DIFFER by {diff}"
    if diff != 0:
        any_diff = True
    print(f"  W^{n}: GR={grr_gr}, STG={grr_stg}  => {match_str}")

# ================================================================
# PART 10: ALTERNATIVE DIRECT COMPUTATION
# ================================================================

print("\n" + "=" * 70)
print("PART 10: DIRECT PPN FROM GEODESIC EQUATION")
print("=" * 70)

print("""
As a cross-check, we compute PPN parameters directly from the
equations of motion, without coordinate transformation.

For a general static spherically symmetric metric:
  ds^2 = -e^(2*Phi) dt^2 + e^(2*Lambda) dr^2 + e^(2*Psi) r^2 dOmega^2

the geodesic equation gives:
  d^2r/dt^2 = -dPhi/dr * (1 - v^2*(1 + e^(2*Lambda-2*Phi))) + ...

For PPN, what matters is the effective potential and the orbital dynamics.

Actually, the cleanest check is:
1. gamma determines the spatial curvature (Shapiro delay, light deflection)
2. beta determines the nonlinearity of g_tt

For the STG metric with C = 1+M/r and alpha=1/2:
  g_tt = -1/C^2 = -1/(1+M/r)^2

The KEY question is whether this g_tt, when expressed in isotropic
coordinates, gives beta=1.

This is exactly what we computed above. Let us verify with a
completely independent method.
""")

# Alternative: use the PPN definition directly
# In Will's PPN formalism, for a static spherically symmetric spacetime:
# g_tt = -(1 - 2*M_eff/r_iso + 2*beta*(M_eff/r_iso)^2 + ...)
# g_rr = 1 + 2*gamma*M_eff/r_iso + ...
# where r_iso is the isotropic radial coordinate

# Let's also check by computing the effective potential
# Vpot = E^2/(2) where the orbit equation is (dr/dtau)^2 + Veff = 0

print("Computing effective potential for circular orbits (alpha=1/2):")
print()

# From the metric: ds^2 = -dt^2/C^2 + C^2 dr^2 + C r^2 dOmega^2  (alpha=1/2)
# Killing vectors: E = dt/dtau * 1/C^2, L = C^(2*1/2) * r^2 * dphi/dtau = C*r^2*dphi/dtau
# Normalization: -1 = g_mu_nu dx^mu/dtau dx^nu/dtau
# -1 = -E^2*C^2 + C^2*(dr/dtau)^2 + L^2/(C*r^2) * C*r^2
# Wait, let me be more careful.

# g_tt = -1/C^2, g_rr = C^2, g_phi = C^(2*1/2)*r^2 = C*r^2 (in equatorial plane)
# E_conserved = |g_tt| * dt/dtau = (1/C^2) * dt/dtau
# L_conserved = g_phi * dphi/dtau = C*r^2 * dphi/dtau

# Normalization:
# -1 = g_tt (dt/dtau)^2 + g_rr (dr/dtau)^2 + g_phi (dphi/dtau)^2
# -1 = -(1/C^2)(E*C^2)^2 + C^2*(dr/dtau)^2 + C*r^2*(L/(C*r^2))^2
# -1 = -E^2*C^2 + C^2*(dr/dtau)^2 + L^2/(C*r^2)

# So: C^2*(dr/dtau)^2 = E^2*C^2 - L^2/(C*r^2) - 1
# (dr/dtau)^2 = E^2 - L^2/(C^3*r^2) - 1/C^2

r_s = symbols('r', positive=True)
M_s = symbols('M', positive=True)
E_s, L_s = symbols('E L', positive=True)

C_s = 1 + M_s/r_s

Veff = L_s**2 / (C_s**3 * r_s**2) + 1/C_s**2
print(f"  (dr/dtau)^2 = E^2 - V_eff(r)")
print(f"  V_eff = L^2/(C^3*r^2) + 1/C^2")

# Circular orbit: dVeff/dr = 0
dVeff_dr = sp.diff(Veff, r_s)
dVeff_dr_simplified = sp.simplify(dVeff_dr)

# For circular orbit at r=r0:
# Also: dr/dtau = 0 => E^2 = Veff(r0)
# And: dVeff/dr|_{r0} = 0

r0 = symbols('r0', positive=True)
C0 = 1 + M_s/r0

# dVeff/dr = 0 gives L^2 in terms of r0
# From the explicit form:
# V = L^2 * r^(-2) * (1+M/r)^(-3) + (1+M/r)^(-2)
# dV/dr = L^2 * [-2r^(-3)(1+M/r)^(-3) + 3M*r^(-4)*(1+M/r)^(-4)] + 2M*r^(-2)*(1+M/r)^(-3)
# = L^2 * r^(-4)*(1+M/r)^(-4) * [-2r(1+M/r) + 3M] + 2M*r^(-2)*(1+M/r)^(-3)
# = L^2 * r^(-4)*(1+M/r)^(-4) * [-2r - 2M + 3M] + 2M*r^(-2)*(1+M/r)^(-3)
# = L^2 * r^(-4)*(1+M/r)^(-4) * [M - 2r] + 2M*r^(-2)*(1+M/r)^(-3) = 0

# So: L^2 * (M-2r) / (r^4 * C^4) + 2M / (r^2 * C^3) = 0
# L^2 = 2M * r^2 * C / (2r - M)  [for r > M/2]

L2_circ = 2*M_s * r0**2 * C0 / (2*r0 - M_s)
print(f"\n  Circular orbit L^2 = {L2_circ}")

# Expand in M/r0
u0 = M_s/r0
L2_series = sp.series(L2_circ.rewrite(u0), u0, 0, n=4)
# Actually sympy doesn't know u0 is M/r0. Let me substitute.
L2_sub = L2_circ.subs(r0, M_s/u0)
L2_sub = sp.simplify(L2_sub)
# This gives L^2 in terms of u0 = M/r0

# Precession: omega_r^2 / omega_phi^2, where
# omega_r^2 = d^2Veff/dr^2 / (2*g_rr) and omega_phi^2 related to L

# The precession per orbit is:
# delta_phi = 2*pi * (omega_phi/omega_r - 1) ~ 2*pi * (1 - omega_r^2/omega_phi^2) / 2

# For PPN: delta_phi = (2+2*gamma-beta)/3 * 6*pi*M/a(1-e^2)
# For circular orbit (e=0): delta_phi = (2+2*gamma-beta)/3 * 6*pi*M/r0

print("\n  PPN precession formula: delta_phi = (2+2*gamma-beta)/3 * 6*pi*M/r0")
print(f"  With gamma=1, beta=1: delta_phi = 6*pi*M/r0 = GR prediction")
print(f"  With gamma={gamma_half}, beta={beta_half}:")
print(f"  Factor = (2+2*{gamma_half}-{beta_half})/3 = {sp.simplify(prec_alpha_half)}")
print()

# ================================================================
# PART 11: SUMMARY
# ================================================================

print("=" * 70)
print("SUMMARY OF RESULTS")
print("=" * 70)

print(f"""
STG METRIC: ds^2 = -dt^2/C^2 + C^2 dr^2 + C^(2*alpha) r^2 dOmega^2
  C(r) = 1 + M/r

ISOTROPIC TRANSFORMATION: R = r * (1 + h1*M/r + h2*(M/r)^2 + ...)
  where h_n are determined order-by-order to make spatial part
  conformally flat.

PPN PARAMETERS (alpha = 1/2):
  gamma = {gamma_half}
  beta  = {beta_half}

PPN PARAMETERS (alpha = 1, used by Gemini):
  gamma = {sp.simplify(gamma_alpha1)}
  beta  = {sp.simplify(beta_alpha1)}

CONSISTENCY CHECKS:
  Light deflection factor (1+gamma)/2 = {defl_alpha_half} [GR = 1]
  Precession factor (2+2*gamma-beta)/3 = {sp.simplify(prec_alpha_half)} [GR = 1]
  Nordtvedt parameter 4*beta-gamma-3 = {sp.simplify(nordtvedt)} [GR = 0]

CONCLUSION ON GEMINI'S CLAIM:
  The Gemini reviewer's beta=1.5 was computed with alpha=1 (conformal metric).
  With the CORRECT value alpha=1/2 (derived geometrically from 600-cell):
    beta = {beta_half}
""")

# Final explicit demonstration
print("EXPLICIT CALCULATION (alpha=1/2):")
print("-" * 50)
print(f"  Step 1: h1 = {h1_half} (from matching u^1)")
print(f"  Step 2: h2 = {h2_half} (from matching u^2)")
print(f"  Step 3: R = r*(1 + ({h1_half})*M/r + ({h2_half})*(M/r)^2 + ...)")
print(f"  Step 4: W = M/R => u = M/r as series in W")
print(f"  Step 5: g_tt = -1/(1+u)^2 expanded in W:")
print(f"         = -(1 - 2W + {a2}*W^2 + {a3}*W^3 + ...)")
print(f"  Step 6: beta = {a2}/2 = {beta_PPN}")
print(f"  Step 7: g_rr_iso = {B_W_series}")
print(f"         gamma = {b1}/2 = {gamma_PPN}")

status = "DERIVAT" if (beta_half == 1 and gamma_half == 1) else "DERIVAT (non-GR)"
print(f"\n  Status: {status}")

if beta_half == 1 and gamma_half == 1:
    print("\n  STG with alpha=1/2 has IDENTICAL 1PN parameters to GR.")
    print("  The theories only diverge at 2PN order or higher.")
else:
    print(f"\n  STG with alpha=1/2 has beta = {beta_half} != 1 in isotropic gauge.")
    print(f"  gamma = {gamma_half} (GR: 1)")
    print(f"  beta = {beta_half} (GR: 1)")

# ================================================================
# PART 11: PPN IN AREAL (STANDARD PPN) GAUGE
# ================================================================

print("\n" + "=" * 70)
print("PART 11: RECONCILING WITH EXP104 (PRECESSION = GR AT 1PN)")
print("=" * 70)

print("""
CRITICAL ISSUE: exp104 showed that STG with alpha=1/2 gives
  precession = 6*pi*M/r0 = GR prediction.

But the standard PPN formula gives:
  precession = (2+2*gamma-beta)/3 * 6*pi*M/r0 = (2+2-2)/3 * 6*pi*M/r0
             = (2/3) * 6*pi*M/r0 = 4*pi*M/r0

This is a CONTRADICTION. Let's resolve it.

The key insight is that the standard PPN precession formula
  delta_phi = (2+2*gamma-beta)/3 * 6*pi*M/p
is derived for metrics where g_theta = r^2 in the ORIGINAL coordinates
(i.e., the areal radius). In the STG metric, the angular part is
C^(2*alpha)*r^2, so r is NOT the areal radius.

The PPN formalism assumes a specific gauge (standard PPN gauge) where
the isotropic coordinate R_iso has a SPECIFIC physical meaning.
When we transform STG to isotropic form, the mapping between the
coordinates introduces corrections that modify the effective orbit equation.

Let me verify by computing precession DIRECTLY from the geodesic
equation in the ORIGINAL STG coordinates, then comparing.
""")

# Direct precession from the STG metric (alpha=1/2)
# ds^2 = -dt^2/C^2 + C^2 dr^2 + C r^2 dOmega^2
# C = 1 + M/r

# Geodesic equation using Binet substitution:
# u = 1/r, L = C * r^2 * dphi/dtau

# From exp104, the orbit equation for "Metric B" (which is our metric with alpha=1/2
# when rewritten with the general alpha form) was:
# d^2u/dphi^2 + u = ... (complicated function)

# For the angular piece: C^(2*alpha)*r^2 with alpha=1/2 gives C*r^2
# Conservation: L = g_phi_phi * dphi/dtau = C * r^2 * dphi/dtau

# From the normalization:
# -1 = -E^2*C^2 + C^2*(dr/dtau)^2 + L^2/(C*r^2)

# Let u = 1/r. Then r' = dr/dphi = (dr/dtau)/(dphi/dtau)
# dphi/dtau = L/(C*r^2) = L*u^2/C
# dr/dtau = r'*dphi/dtau = r'*L*u^2/C

# Also du/dphi = -u^2 * dr/dphi = -u^2 * r' = -u^2 * (dr/dtau)/(dphi/dtau)

# Let me denote ' = d/dphi
# du/dphi = -1/r^2 * dr/dphi

# From normalization:
# C^2*(dr/dtau)^2 = E^2*C^2 - 1 - L^2/(C*r^2)
# (dr/dtau)^2 = E^2 - 1/C^2 - L^2/(C^3*r^2)
# Also (dr/dtau) = (dr/dphi)*(dphi/dtau) = (dr/dphi) * L/(C*r^2)
# (dr/dphi)^2 * L^2/(C^2*r^4) = E^2 - 1/C^2 - L^2/(C^3*r^2)

# u'' = -(d/dphi)(u^2 * r') where r' = dr/dphi...
# This gets messy. Let me use sympy for the orbit equation.

print("Computing precession directly from STG geodesics...")
print()

# Work with u = M/r as the variable (not 1/r, to keep things dimensionless)
# Then r = M/u
u_var = symbols('u', positive=True)
phi = symbols('phi')
L_var = symbols('L', positive=True)
E_var = symbols('E', positive=True)
M_var = symbols('M', positive=True)

C_u = 1 + u_var  # C = 1 + M/r = 1 + u where u = M/r

# The effective potential approach:
# (dr/dtau)^2 = E^2 - L^2/(C^3*r^2) - 1/C^2
# r = M/u, so dr = -M/u^2 * du
# (M^2/u^4)*(du/dtau)^2 = E^2 - L^2*u^2/(M^2*C^3) - 1/C^2

# Also du/dphi = (du/dtau)/(dphi/dtau)
# dphi/dtau = L/(C*r^2) = L*u^2/(M^2*C)

# So du/dtau = (du/dphi) * L*u^2/(M^2*C)
# (M^2/u^4) * (du/dphi)^2 * L^2*u^4/(M^4*C^2) = E^2 - L^2*u^2/(M^2*C^3) - 1/C^2
# L^2/(M^2*C^2) * (du/dphi)^2 = E^2 - L^2*u^2/(M^2*C^3) - 1/C^2

# Let h = L^2/M^2 (dimensionless)
# (du/dphi)^2 = C^2*(E^2 - 1/C^2) - h*u^2/C = E^2*C^2 - 1 - h*u^2/C

# Take derivative d/dphi:
# 2*u''*u' = 2*E^2*C*C'*u' - (-h*2*u*u'/C + h*u^2*C'*u'/C^2)
# where C' = dC/du = 1 (since C = 1+u)
# Divide by 2*u':
# u'' = E^2*C - (-h*u/C + h*u^2/(2*C^2))
# Wait this isn't right. Let me be more careful.

# (du/dphi)^2 = E^2*(1+u)^2 - 1 - h*u^2/(1+u)
# Differentiate both sides w.r.t. phi:
# 2*(du/dphi)*(d^2u/dphi^2) = [2*E^2*(1+u) - h*(2*u*(1+u) - u^2)/(1+u)^2] * du/dphi
# Cancel du/dphi:
# 2*(d^2u/dphi^2) = 2*E^2*(1+u) - h*(2*u + 2*u^2 - u^2)/(1+u)^2
# 2*u'' = 2*E^2*(1+u) - h*(2*u + u^2)/(1+u)^2
# u'' = E^2*(1+u) - h*u*(2+u)/(2*(1+u)^2)

# For circular orbit at u = u0: u' = 0, u'' = 0
# From u'' = 0: E^2*(1+u0) = h*u0*(2+u0)/(2*(1+u0)^2)
# E^2 = h*u0*(2+u0)/(2*(1+u0)^3)

# From (du/dphi)^2 = 0: E^2*(1+u0)^2 - 1 = h*u0^2/(1+u0)
# Substitute E^2:
# h*u0*(2+u0)*(1+u0)^2/(2*(1+u0)^3) - 1 = h*u0^2/(1+u0)
# h*u0*(2+u0)/(2*(1+u0)) - 1 = h*u0^2/(1+u0)
# h*u0/(2*(1+u0)) * [(2+u0) - 2*u0] = 1
# h*u0*(2-u0)/(2*(1+u0)) = 1
# h = 2*(1+u0)/(u0*(2-u0))

h_circ = 2*(1+u_var)/(u_var*(2-u_var))

print(f"  Circular orbit: h = L^2/M^2 = {h_circ}")
print(f"  (Note: u = M/r, so u << 1 for weak field)")
print(f"  h ~ 1/u for u << 1 => L^2 ~ M*r [Newtonian!]")

# E^2 for circular:
E2_circ = h_circ*u_var*(2+u_var)/(2*(1+u_var)**3)
E2_circ_simp = sp.simplify(E2_circ)
print(f"  E^2_circ = {E2_circ_simp}")

# Radial oscillation frequency:
# Perturb: u = u0 + delta_u, du/dphi = delta_u'
# delta_u'' + omega_r^2 * delta_u = 0
# omega_r^2 = -d(RHS)/du at u=u0 where RHS = u'' expression

# u'' = E^2*(1+u) - h*u*(2+u)/(2*(1+u)^2)
# d(u'')/du = E^2 - h*[d/du(u*(2+u)/(1+u)^2)] / 2

# d/du[u*(2+u)/(1+u)^2] = [(2+2u)(1+u)^2 - u*(2+u)*2*(1+u)] / (1+u)^4
# = [(2+2u)(1+u) - 2u(2+u)] / (1+u)^3
# = [2+4u+2u^2 - 4u-2u^2] / (1+u)^3
# = 2/(1+u)^3

deriv_bracket = 2/(1+u_var)**3

# omega_r^2 = -(E^2 - h/(1+u)^3) evaluated at circular orbit
# omega_r^2 = h/(1+u0)^3 - E^2
# Substitute E^2 = h*u0*(2+u0)/(2*(1+u0)^3)
# omega_r^2 = h/(1+u0)^3 * [1 - u0*(2+u0)/2]
# = h/(1+u0)^3 * [(2 - 2u0 - u0^2)/2]

omega_r_sq = h_circ/(1+u_var)**3 * (2 - 2*u_var - u_var**2)/2
omega_r_sq_simp = sp.simplify(omega_r_sq)
print(f"\n  omega_r^2 = {omega_r_sq_simp}")

# Angular frequency: omega_phi = 1 (by definition, phi is the angle)
# But wait - in Binet approach, the "angle" is coordinate phi.
# The precession is: delta_phi_per_orbit = 2*pi*(1/omega_r - 1)
# ~ 2*pi*(1 - omega_r^2/2) for omega_r close to 1

# Actually: period of radial oscillation T_r = 2*pi/omega_r
# Number of angular turns per radial period = T_r/(2*pi) * omega_phi = 1/omega_r
# Wait no. The orbit equation gives u as function of phi.
# If u'' + omega_r^2 * u ~ 0, then the radial oscillation
# has period 2*pi/omega_r in phi.
# So the apsidal angle is pi/omega_r (half period).
# Precession per orbit = 2*(pi/omega_r) - 2*pi = 2*pi*(1/omega_r - 1)

# Expand omega_r^2 in powers of u0:
omega_r_sq_series = sp.series(omega_r_sq_simp, u_var, 0, n=3)
print(f"  omega_r^2 expanded: {omega_r_sq_series}")

# omega_r = sqrt(omega_r^2) ~ 1 - (1-omega_r^2)/2 for omega_r^2 close to 1
# 1/omega_r ~ 1 + (1-omega_r^2)/2 + ...
# precession = 2*pi*(1/omega_r - 1) ~ pi*(1 - omega_r^2) + ...

omega_r_sq_coeff0 = omega_r_sq_series.coeff(u_var, 0)
omega_r_sq_coeff1 = omega_r_sq_series.coeff(u_var, 1)
omega_r_sq_coeff2 = omega_r_sq_series.coeff(u_var, 2)

print(f"\n  omega_r^2 = {omega_r_sq_coeff0} + {omega_r_sq_coeff1}*u + {omega_r_sq_coeff2}*u^2 + ...")
print(f"  1 - omega_r^2 = {1-omega_r_sq_coeff0} + {-omega_r_sq_coeff1}*u + {-omega_r_sq_coeff2}*u^2 + ...")

prec_direct = sp.pi * (1 - omega_r_sq_coeff0) + sp.pi * (-omega_r_sq_coeff1) * u_var
print(f"\n  Precession (1PN) = {sp.simplify(prec_direct)} (where u = M/r0)")
print(f"  = {sp.simplify(prec_direct).coeff(u_var)} * pi * M/r0")

# Multiply: coefficient * pi * M/r0
prec_coeff = sp.simplify(-omega_r_sq_coeff1)
print(f"\n  Precession = {prec_coeff} * pi * M/r0 per orbit")
print(f"  GR prediction = 6*pi*M/r0")

# Now let's also compute the FULL omega_r^2 without approximation
print(f"\n  Full omega_r^2 = {omega_r_sq_simp}")
# Let's factor it
omega_r_sq_factored = sp.factor(omega_r_sq_simp)
print(f"  Factored: omega_r^2 = {omega_r_sq_factored}")

# ================================================================
# PART 12: THE AREAL RADIUS ISSUE
# ================================================================

print("\n" + "=" * 70)
print("PART 12: AREAL RADIUS AND PPN GAUGE")
print("=" * 70)

print("""
The PPN formalism in its standard form (Will 1993, 2014) uses the
"standard PPN gauge" where the radial coordinate R satisfies:
  g_theta_theta = R^2

In the STG metric, g_theta_theta = C^(2*alpha)*r^2 = C*r^2 (for alpha=1/2).
So the areal radius is:
  R_areal = sqrt(g_theta_theta) = r * sqrt(C) = r * sqrt(1+M/r)

This is DIFFERENT from both r (STG coordinate) and R_iso (isotropic).

The PPN parameters should be computed in terms of R_areal.
Let's redo the computation using the areal radius.
""")

# In the areal radius: R_a = r*sqrt(1+M/r) for alpha=1/2
# C = 1 + M/r => M/r = C-1, r = M/(C-1)
# R_a = r*C^(1/2) = M*C^(1/2)/(C-1)

# Express metric in terms of R_a:
# First find r(R_a) perturbatively.
# R_a = r*sqrt(1+M/r) = r*(1 + M/(2r) - M^2/(8r^2) + ...) = r + M/2 - M^2/(8r) + ...
# So r = R_a - M/2 + M^2/(8*R_a) + ...  (perturbative)

# Let V = M/R_a (weak field parameter in areal coords)
V = symbols('V')

# r = R_a - M/2 + ... => M/r = M/(R_a - M/2 + ...) = V/(1 - V/2 + ...)
# = V*(1 + V/2 + 3V^2/8 + ...) = V + V^2/2 + 3*V^3/8 + ...
# Actually let me invert R_a = r + M/2 - M^2/(8r) + M^3/(16r^2) + ...
# perturbatively.

# R_a^2 = r^2 * (1+M/r) = r^2 + M*r = r^2 + M*r
# So R_a^2 = r^2 + M*r => r = (R_a^2)/(r+M)... circular.
# Perturbatively: r = R_a - M/2 + d2*M^2/R_a + d3*M^3/R_a^2 + ...
# R_a = r + M/2 - M^2/(8r) + M^3/(16r^2) - 5*M^4/(128r^3) + ...

# Let me just do this symbolically
# u_r = M/r, V = M/R_a
# R_a = r*(1+u_r)^(1/2) => V = u_r/(1+u_r)^(1/2)
# u_r = V*(1+u_r)^(1/2) => u_r^2 = V^2*(1+u_r) => u_r^2 - V^2*u_r - V^2 = 0
# u_r = (V^2 + sqrt(V^4 + 4V^2))/2 = (V^2 + V*sqrt(V^2+4))/2
# For small V: sqrt(V^2+4) ~ 2*(1+V^2/8-V^4/128+...)
# u_r ~ (V^2 + 2V + V^3/4 + ...)/2 = V + V^2/2 + V^3/8 + ...

# Let's expand to sufficient order
u_of_V = V + V**2/2 + V**3/8 - V**5/128  # to 5th order in V... let me compute properly

# Actually let me solve u_r^2 - V^2*u_r - V^2 = 0 perturbatively
# u = V + a2*V^2 + a3*V^3 + a4*V^4 + ...
a2v, a3v, a4v, a5v = symbols('a2v a3v a4v a5v')
u_trial = V + a2v*V**2 + a3v*V**3 + a4v*V**4 + a5v*V**5
# u^2 - V^2*u - V^2 = 0
constraint = sp.expand(u_trial**2 - V**2*u_trial - V**2)
# Match each power of V to zero
for n in range(2, 7):
    cn = constraint.coeff(V, n)
    if n == 2:
        # V^2: 2*a2v - 1 = 0 (wait, let me compute)
        pass

# Easier: just use sympy series
u_exact = (V**2 + V*sp.sqrt(V**2 + 4))/2  # exact solution of u^2 - V^2*u - V^2 = 0
u_areal_series = sp.series(u_exact, V, 0, n=6).removeO()
u_areal_series = sp.expand(u_areal_series)
print(f"  u = M/r as function of V = M/R_areal:")
print(f"  u(V) = {u_areal_series}")

# Now compute g_tt in areal coords
# g_tt = -1/(1+u)^2
g_tt_areal = 1/(1 + u_areal_series)**2
g_tt_areal_series = sp.series(g_tt_areal, V, 0, n=5).removeO()
g_tt_areal_series = sp.expand(g_tt_areal_series)
print(f"\n  g_tt (areal) = -1/(1+u)^2 = -{g_tt_areal_series}")

a2_areal = g_tt_areal_series.coeff(V, 2)
a1_areal = g_tt_areal_series.coeff(V, 1)
beta_areal = a2_areal / 2
print(f"\n  Coeff of V^1: {a1_areal}")
print(f"  Coeff of V^2: {a2_areal}")
print(f"  => beta_areal = {a2_areal}/2 = {beta_areal}")

# For the spatial part: in areal coordinates, g_theta = R_a^2 by definition.
# g_rr_areal * dR_a^2 = g_rr_STG * dr^2
# g_rr_areal = g_rr_STG * (dr/dR_a)^2 = C^2 * (dr/dR_a)^2

# dR_a/dr = d/dr[r*C^(1/2)] = C^(1/2) + r*(1/2)*C^(-1/2)*(-M/r^2)
#          = C^(1/2)*(1 - M/(2r*C)) = C^(1/2)*(1 - u/(2*(1+u)))
#          = C^(1/2)*((2+u)/(2*(1+u)))

# g_rr_areal = C^2 / (dR_a/dr)^2 = C^2 / (C * (2+u)^2/(4*(1+u)^2))
#            = 4*C*(1+u)^2 / (2+u)^2

g_rr_areal_func = 4*(1+u_areal_series)*(1+u_areal_series)**2 / (2+u_areal_series)**2
g_rr_areal_series = sp.series(g_rr_areal_func, V, 0, n=5).removeO()
g_rr_areal_series = sp.expand(g_rr_areal_series)
print(f"\n  g_rr (areal) = 4*C*(1+u)^2/(2+u)^2 = {g_rr_areal_series}")

b1_areal = g_rr_areal_series.coeff(V, 1)
gamma_areal = b1_areal / 2
print(f"  Coeff of V^1: {b1_areal}")
print(f"  => gamma_areal = {b1_areal}/2 = {gamma_areal}")

# PPN precession with areal-coordinate parameters
prec_factor_areal = (2 + 2*gamma_areal - beta_areal) / 3
print(f"\n  Precession factor (areal gauge): (2+2*gamma-beta)/3 = {sp.simplify(prec_factor_areal)}")

# Compare with Schwarzschild in areal gauge
# Schwarzschild: g_tt = -(1-2V), g_rr = 1/(1-2V)
# g_tt = -(1 - 2V)  => coeff V^1=-2, V^2=0 => beta=0??
# Wait, that's Schwarzschild in Schwarzschild coordinates where r IS the areal radius.
# g_tt = -(1 - 2M/r) = -(1 - 2V)
# But PPN convention: g_tt = -(1 - 2V + 2*beta*V^2)
# For Schwarzschild: g_tt = -(1-2V), so the V^2 coefficient is 0, giving beta=0?
# That can't be right. The issue is that in Schwarzschild coords, the spatial
# metric is NOT isotropic: g_rr = 1/(1-2V) != g_theta/r^2 = 1.

print("""
IMPORTANT NOTE ON PPN GAUGE:
The standard PPN formalism by Will uses ISOTROPIC coordinates,
not areal radius coordinates. The PPN parameters are defined via:
  g_tt = -(1 - 2U + 2*beta*U^2 + ...)
  g_{ij} = (1 + 2*gamma*U)*delta_{ij} + ...
where U = M/R_iso is the potential in isotropic coordinates.

The reason is that PPN deals with general matter configurations,
not just spherical symmetry, and uses a specific gauge.

For Schwarzschild in isotropic coords: U = M/R_iso,
  r_Schw = R_iso*(1+U/2)^2
  g_tt = -((1-U/2)/(1+U/2))^2 = -(1 - 2U + 2U^2 - ...) => beta=1
  g_rr_iso = (1+U/2)^4 = 1 + 2U + ... => gamma=1

For STG in isotropic coords (as computed in Part 4):
  g_tt = -(1 - 2W + 4W^2 - ...) => beta=2
  g_rr_iso = 1 + 2W + ... => gamma=1

So the PPN result IS beta=2 for STG with alpha=1/2.
""")

# ================================================================
# PART 13: RESOLVING THE PRECESSION PARADOX
# ================================================================

print("=" * 70)
print("PART 13: RESOLVING THE PRECESSION PARADOX")
print("=" * 70)

print("""
The paradox: exp104 found precession = 6*pi*M/r0 matching GR,
but PPN gives (2+2-2)/3 * 6*pi*M/r0 = 4*pi*M/r0.

RESOLUTION: The standard PPN precession formula
  delta_phi = (2+2*gamma-beta)/3 * 6*pi*M/p
uses the semi-latus rectum p IN ISOTROPIC COORDINATES.

In exp104, precession was computed using the STG coordinate r0.
The semi-latus rectum differs between these coordinate systems.

For a circular orbit at STG coordinate r0:
  R_iso = r0 * H(M/r0) = r0 * (1 - M/(2*r0) + ...)
  p_iso = R_iso for circular orbit

So M/p_iso = M/R_iso = M/(r0*(1-M/(2*r0)+...)) = (M/r0)*(1 + M/(2*r0) + ...)

The precession in terms of M/r0:
  delta_phi = (2+2*gamma-beta)/3 * 6*pi * M/p_iso
            = (2/3) * 6*pi * (M/r0)*(1 + M/(2*r0) + ...)
            = 4*pi*M/r0 + higher order

But exp104 found 6*pi*M/r0. Let me re-examine.
""")

# Actually, the issue is more subtle. Let me compute omega_r^2 properly.
# From Part 11:
print("From direct geodesic computation (Part 11):")
print(f"  omega_r^2 = {omega_r_sq_simp}")
omega_r_sq_series2 = sp.series(omega_r_sq_simp, u_var, 0, n=4)
print(f"  Expanded: omega_r^2 = {omega_r_sq_series2}")

# Extract 1PN coefficient
c0 = omega_r_sq_series2.coeff(u_var, 0)
c1 = omega_r_sq_series2.coeff(u_var, 1)
c2 = omega_r_sq_series2.coeff(u_var, 2)

print(f"\n  omega_r^2 = {c0} + ({c1})*u + ({c2})*u^2 + ...")
print(f"  (where u = M/r in STG coordinates)")

# precession per orbit = 2*pi*(1/omega_r - 1) = 2*pi*(1-omega_r^2)/(2*omega_r^2) ~ pi*(1-omega_r^2)
prec_1pn = sp.pi * (-c1)
print(f"\n  Precession (1PN) = {prec_1pn} * u = {prec_1pn} * M/r0")

if prec_1pn == 6*sp.pi:
    print("  = 6*pi*M/r0  [MATCHES GR in STG coordinates!]")
elif prec_1pn == 4*sp.pi:
    print("  = 4*pi*M/r0  [DIFFERENT from GR's 6*pi*M/r0!]")
else:
    print(f"  = {prec_1pn/sp.pi}*pi*M/r0")

# If NOT 6*pi*M/r0, then exp104 was WRONG about matching GR at 1PN!
# Let's be very explicit about this.

print(f"""
DIRECT COMPUTATION RESULT:
  omega_r^2 = 1 + ({c1})*u + O(u^2)  where u = M/r_STG
  precession = {-c1}*pi*M/r_STG per orbit

COMPARISON WITH GR:
  GR: omega_r^2 = 1 - 6*u + O(u^2)  where u = M/r_Schw
  GR precession = 6*pi*M/r_Schw per orbit

NOTE: r_STG and r_Schw have DIFFERENT physical meanings!
  r_Schw = areal radius (area of sphere = 4*pi*r^2)
  r_STG = coordinate such that g_theta = C^(2*alpha)*r^2
""")

# Let's also verify with GR using same method
# GR: V_eff = L^2/(2r^2) * (1-2M/r) + (1-2M/r)/2 ... no, let me use
# the standard GR orbit equation
# (du/dphi)^2 + u^2 = E^2 - (1-2u)(1 + h_GR*u^2) ... actually the standard
# GR result is well-known: omega_r^2 = 1 - 6u => precession = 6*pi*u per orbit

print("VERIFYING: For GR (Schwarzschild), omega_r^2 = 1 - 6u is standard.")
print(f"For STG (alpha=1/2), omega_r^2 = 1 + ({c1})u = 1 - {-c1}u")

# ================================================================
# PART 14: NUMERICAL VERIFICATION
# ================================================================

print("\n" + "=" * 70)
print("PART 14: NUMERICAL VERIFICATION OF PRECESSION")
print("=" * 70)

import numpy as np
from scipy.integrate import solve_ivp

M_num = 1.0

def compute_stg_precession(r0, e=0.01, alpha_metric=0.5, n_orbits=20):
    """Compute precession by integrating the orbit ODE for STG metric."""
    u0 = M_num / r0

    # For the general metric: ds^2 = -dt^2/C^2 + C^2 dr^2 + C^(2*alpha)*r^2 dOmega^2
    # with C = 1 + M/r = 1 + u (where u = M/r)
    #
    # Conservation: L = C^(2*alpha)*r^2*dphi/dtau, E = (1/C^2)*dt/dtau
    # Normalization: -1 = -(EC^2)^2/C^2 + C^2*(dr/dtau)^2 + L^2/(C^(2*alpha)*r^2)
    #              = -E^2*C^2 + C^2*(dr/dtau)^2 + L^2/(C^(2*alpha)*r^2)
    #
    # Using u = M/r, phi as independent variable:
    # L = C^(2*alpha)*r^2*dphi/dtau => dphi/dtau = L*u^2/(M^2*C^(2*alpha))
    #
    # (du/dphi)^2 = [E^2*C^2 - 1 - L^2*u^2/(M^2*C^(2*alpha))] * C^(4*alpha) / (L^2*C^2/M^2)
    # Hmm this is getting messy for general alpha. Let me use dr/dphi directly.

    # Use r as variable, phi as parameter
    # (dr/dphi)^2 = (dr/dtau)^2 / (dphi/dtau)^2

    # From normalization:
    # C^2*(dr/dtau)^2 = E^2*C^2 - 1 - L^2/(C^(2*alpha)*r^2)
    # (dr/dtau)^2 = E^2 - 1/C^2 - L^2/(C^(2+2*alpha)*r^2)

    # dphi/dtau = L/(C^(2*alpha)*r^2)
    # (dr/dphi)^2 = (E^2 - 1/C^2 - L^2/(C^(2+2*alpha)*r^2)) * C^(4*alpha)*r^4/L^2

    # For circular orbit at r0:
    # E^2 = 1/C0^2 + L^2/(C0^(2+2*alpha)*r0^2)  [from dr/dtau = 0]
    # dVeff/dr = 0 gives another equation for L

    # Let me just integrate numerically with alpha=1/2
    C0 = 1 + M_num/r0

    # For alpha=1/2: Veff = L^2/(C^3*r^2) + 1/C^2
    # dVeff/dr = 0: L^2 = 2*M*r^2*C/(2r-M) (from Part 10)
    L2 = 2*M_num*r0**2*C0/(2*r0 - M_num)
    E2 = L2/(C0**3*r0**2) + 1/C0**2

    # Orbit equation: u = M/r, u'' = E^2*(1+u) - h*u*(2+u)/(2*(1+u)^2)
    # where h = L^2/M^2
    h = L2/M_num**2

    def rhs(u):
        return E2*(1+u) - h*u*(2+u)/(2*(1+u)**2)

    u_peri = u0 * (1 + e)

    def ode(phi, y):
        u, up = y
        return [up, rhs(u) - 0]  # u'' = rhs(u), NOT u'' + u = something

    # Wait, I need to re-derive. The orbit eq is:
    # u'' = E^2*(1+u) - h*u*(2+u)/(2*(1+u)^2)
    # This is NOT in the form u'' + u = f(u).

    # For small oscillation around u0:
    # u'' = rhs(u0) + rhs'(u0)*(u-u0) + ...
    # At circular orbit, rhs(u0) = 0 (by construction)
    # u'' + (-rhs'(u0))*(u-u0) = 0
    # omega_r^2 = -rhs'(u0)

    # For numerical integration, event detection at u' = 0 (perihelion)
    def event(phi, y):
        return y[1]  # u' = 0
    event.direction = -1  # u goes from increasing to decreasing

    phi_max = n_orbits * 2 * np.pi * 1.5
    sol = solve_ivp(ode, [0, phi_max], [u_peri, 0.0],
                    events=event, rtol=1e-13, atol=1e-15, max_step=0.01)

    peri_phis = sol.t_events[0]
    if len(peri_phis) >= 4:
        d_phis = np.diff(peri_phis)
        # Average, skip first
        mean_period = np.mean(d_phis[1:])
        precession = mean_period - 2*np.pi
        return precession
    return None

print("Numerical precession for STG (alpha=1/2):")
print(f"  {'r0/M':>8} {'prec_num':>14} {'6*pi*M/r0':>14} {'prec_predict':>14} {'ratio':>10}")

# Also compute omega_r^2 analytically for comparison
for r0 in [50, 100, 200, 500, 1000]:
    u0 = M_num/r0
    prec_num = compute_stg_precession(r0, e=0.001)
    prec_gr = 6*np.pi*M_num/r0

    # omega_r^2 from our formula
    omega2 = float(omega_r_sq_simp.subs(u_var, u0))
    omega_r = np.sqrt(omega2)
    prec_analytic = 2*np.pi*(1/omega_r - 1)

    if prec_num is not None:
        ratio = prec_num / prec_gr
        print(f"  {r0:>8.0f} {prec_num:>14.8f} {prec_gr:>14.8f} {prec_analytic:>14.8f} {ratio:>10.6f}")

# ================================================================
# PART 15: FINAL RESOLUTION
# ================================================================

print("\n" + "=" * 70)
print("PART 15: FINAL RESOLUTION AND CONCLUSIONS")
print("=" * 70)

print(f"""
RESULTS SUMMARY:

1. PPN PARAMETERS (STG, alpha=1/2, isotropic gauge):
   gamma = {gamma_half}
   beta = {beta_half}
   General formula: beta(alpha) = 5/2 - alpha

2. GEMINI'S CLAIM (beta=1.5):
   Correct computation for alpha=1 (conformal metric).
   But STG uses alpha=1/2, so beta=2 (not 1.5).
   Gemini used the WRONG alpha value.

3. PRECESSION:
   Standard PPN formula with beta=2, gamma=1 gives:
   factor = (2+2-2)/3 = 2/3, so precession = 4*pi*M/r0

   But this uses isotropic radial coordinate.
   Direct geodesic computation (omega_r^2) gives the result
   in STG coordinates, which may differ.

4. KEY FINDING: beta(alpha) = 5/2 - alpha for all alpha.
   For beta=1 (GR value) we need alpha = 3/2.
   NO value of alpha gives both beta=1 AND gamma=1 simultaneously
   (gamma=1 for all alpha, but beta depends on alpha).

5. STATUS: DERIVAT
   The computation is rigorous symbolic algebra.
   The result beta=2 for alpha=1/2 is mathematically correct.

6. IMPLICATIONS:
   STG is NOT equivalent to GR at 1PN level.
   beta=2 is strongly excluded by Solar System tests
   (Cassini: |gamma-1| < 2.3e-5, LLR: |beta-1| < 1.1e-4).

   HOWEVER: the exp104 precession result (matching GR) was computed
   in STG coordinates. The coordinate transformation to isotropic
   form introduces the beta!=1 effect. The physical observables
   (precession angle as measured by an observer) depend on the
   full coordinate transformation, not just the leading PPN terms.

   This needs further investigation: compute the PHYSICAL precession
   angle (coordinate-invariant) to determine if STG truly deviates
   from GR at 1PN.
""")

print("=" * 70)
print("END OF EXP-155")
print("=" * 70)

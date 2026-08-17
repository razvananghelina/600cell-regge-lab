"""
exp427b_analytical_angles.py — Analytical dihedral angles for S^4 double-cone
==============================================================================
Derives exact formulas for theta_A, theta_B using Gram matrix + Schur complement.

4-simplex: (N, a, b, c, d) with edges r=sqrt(2) (radial), a=1/phi (equatorial).
All 1200 four-simplices are congruent.
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
sqrt5 = np.sqrt(5)

print("=" * 72)
print("ANALYTICAL DIHEDRAL ANGLES FOR S^4 DOUBLE-CONE")
print("=" * 72)

# ============================================================
# GRAM MATRIX from pole N as reference vertex
# ============================================================
# Edge vectors: e_i = v_i - N, for i = a,b,c,d (equatorial vertices)
# G_ii = |e_i|^2 = r^2 = 2
# G_ij = (r^2 + r^2 - a^2)/2 = 2 - a^2/2  where a^2 = 1/phi^2 = 2-phi
#
# So: G_ij = 2 - (2-phi)/2 = 1 + phi/2   for i != j
#     G_ii = 2

r2 = 2          # r^2
a2 = 2 - PHI    # a^2 = 1/phi^2 = 2 - phi
c = r2 - a2/2   # off-diagonal Gram entry = 1 + phi/2

print(f"\n  r^2 = {r2}")
print(f"  a^2 = 1/phi^2 = 2-phi = {a2:.6f}")
print(f"  c = 1 + phi/2 = {c:.6f}")
print(f"  Verify: c = (5+sqrt(5))/4 = {(5+sqrt5)/4:.6f}")

# ============================================================
# THETA_A: dihedral at RADIAL triangle (N, v_i, v_j)
# ============================================================
# Triangle spanned by {e_1, e_2}. Others: {e_3, e_4}.
# Schur complement G' = G_NN - G_NT * G_TT^{-1} * G_TN
# By symmetry (G_TT = [[2,c],[c,2]], G_TN = c*J):
#
#   cos(theta_A) = c / (r^2 + 2c)
#
# With c = 1+phi/2, r^2 = 2:
#   r^2 + 2c = 2 + 2 + phi = 4 + phi
#   cos(theta_A) = (1+phi/2) / (4+phi) = (2+phi) / (2*(4+phi))
#
# CLEAN FORM:  cos(theta_A) = 1/2 - 1/(4+phi)

print(f"\n{'='*72}")
print("THETA_A (radial triangle)")
print(f"{'='*72}")

cos_A = c / (r2 + 2*c)
theta_A = np.arccos(cos_A)

print(f"\n  DERIVATION: Schur complement of Gram matrix")
print(f"  cos(theta_A) = c/(r^2 + 2c)")
print(f"               = (1+phi/2)/(4+phi)")
print(f"               = (2+phi)/(2*(4+phi))")
print(f"               = 1/2 - 1/(4+phi)")

# Verify all forms
form1 = (2 + PHI) / (2*(4 + PHI))
form2 = 0.5 - 1/(4 + PHI)
form3 = (10 + sqrt5) / 38  # rationalized

print(f"\n  Numerical verification:")
print(f"    (2+phi)/(2*(4+phi))  = {form1:.10f}")
print(f"    1/2 - 1/(4+phi)     = {form2:.10f}")
print(f"    (10+sqrt(5))/38     = {form3:.10f}")
print(f"    From Gram matrix    = {cos_A:.10f}")
print(f"    All agree: {np.allclose([form1, form2, form3], cos_A)}")

print(f"\n  theta_A = arccos((2+phi)/(2*(4+phi)))")
print(f"         = {np.degrees(theta_A):.6f} deg")
print(f"         = {theta_A:.10f} rad")

# ============================================================
# THETA_B: dihedral at EQUATORIAL triangle (v_i, v_j, v_k)
# ============================================================
# Reference vertex: v_i (equatorial). Vectors to: {N, v_j-origin, v_k, v_l}
# Triangle spanned by {g_c, g_d} (equatorial dirs). Others: {g_a, g_N}.
# Gram matrix K has K_TT = [[2-phi, q],[q, 2-phi]] with q = (2-phi)/2
# Schur complement gives:
#
#   cos^2(theta_B) = (2-phi) / (8*(4+phi))
#
# Since 2-phi = 1/phi^2:
#   cos^2(theta_B) = 1 / (8*phi^2*(phi^2+3))
#
# And phi^2*(phi^2+3) = (phi+1)*(phi+4) = 6*phi+5 = 8+3*sqrt(5)
#   cos^2(theta_B) = 1 / (8*(6*phi+5))
#
# Rationalized: cos^2(theta_B) = (8-3*sqrt(5)) / (8*19)

print(f"\n{'='*72}")
print("THETA_B (equatorial triangle)")
print(f"{'='*72}")

cos2_B = (2 - PHI) / (8*(4 + PHI))
cos_B = np.sqrt(cos2_B)
theta_B = np.arccos(cos_B)

print(f"\n  DERIVATION: Schur complement with equatorial reference")
print(f"  cos^2(theta_B) = (2-phi) / (8*(4+phi))")
print(f"                 = 1 / (8*phi^2*(phi^2+3))")
print(f"                 = 1 / (8*(6*phi+5))")
print(f"                 = (8-3*sqrt(5)) / (8*19)")

# Verify all forms
f1 = 1 / (8*PHI**2*(PHI**2+3))
f2 = 1 / (8*(6*PHI+5))
f3 = (8 - 3*sqrt5) / (8*19)

print(f"\n  Numerical verification:")
print(f"    (2-phi)/(8*(4+phi))        = {cos2_B:.10f}")
print(f"    1/(8*phi^2*(phi^2+3))      = {f1:.10f}")
print(f"    1/(8*(6*phi+5))            = {f2:.10f}")
print(f"    (8-3*sqrt(5))/(8*19)       = {f3:.10f}")
print(f"    All agree: {np.allclose([f1, f2, f3], cos2_B)}")

print(f"\n  theta_B = arccos(sqrt((2-phi)/(8*(4+phi))))")
print(f"         = {np.degrees(theta_B):.6f} deg")
print(f"         = {theta_B:.10f} rad")

# ============================================================
# DEFICIT ANGLES (exact)
# ============================================================
print(f"\n{'='*72}")
print("DEFICIT ANGLES")
print(f"{'='*72}")

eps_eq = 2*np.pi - 4*theta_B
eps_rad = 2*np.pi - 5*theta_A

print(f"  eps_eq  = 2*pi - 4*arccos(sqrt(1/(8*(6*phi+5))))")
print(f"         = {eps_eq:.10f} rad = {np.degrees(eps_eq):.6f} deg")
print(f"  eps_rad = 2*pi - 5*arccos((2+phi)/(2*(4+phi)))")
print(f"         = {eps_rad:.10f} rad = {np.degrees(eps_rad):.6f} deg")

# ============================================================
# THE NUMBER 19
# ============================================================
print(f"\n{'='*72}")
print("THE NUMBER 19 = 4*a1 - 1")
print(f"{'='*72}")
a1 = 5
print(f"\n  In theta_A (rationalized):")
print(f"    cos(theta_A) = (10+sqrt(5))/38 = (2*a1+sqrt(a1))/(2*(4*a1-1))")
print(f"    Denominator: 2*19 = 2*(4*a1-1)")
print(f"\n  In theta_B (rationalized):")
print(f"    cos^2(theta_B) = (8-3*sqrt(5))/(8*19)")
print(f"    Denominator: 8*19 = 8*(4*a1-1)")
print(f"\n  Same 19 as in Z[phi] norm formula (exp416):")
print(f"    N(z) = -(n^2 + 9nt + 19t^2)")
print(f"    prod(N) = -a1 * 19^2 = -{a1*19**2}")
print(f"\n  Origin: 19 = (9+sqrt(5))(9-sqrt(5))/4 = (2*a1-1)^2-a1)/(a1-1)")
print(f"        = (4*a1^2 - 5*a1 + 1)/(a1-1) = 4*a1 - 1")

# ============================================================
# GENERAL a1 FORMULAS
# ============================================================
print(f"\n{'='*72}")
print("GENERAL a1 FORMULAS")
print(f"{'='*72}")
print(f"""
  For general a1 (phi = (1+sqrt(a1))/2, a^2 = 2-phi, r^2 = 2):

  cos(theta_A) = (2+phi)/(2*(4+phi))
               = (2*a1 + sqrt(a1)) / (2*(4*a1-1))

  cos^2(theta_B) = (2-phi)/(8*(4+phi))
                 = (a1+3-3*sqrt(a1)) / (8*(4*a1-1))

  Both denominators contain 4*a1-1.
  For a1=5: 4*5-1 = 19.
""")

# ============================================================
# REGGE ACTION (analytical)
# ============================================================
print(f"{'='*72}")
print("REGGE ACTION (analytical)")
print(f"{'='*72}")

A_eq = np.sqrt(3)/4 * a2  # sqrt(3)/4 * (1/phi)^2
A_rad = 0.5 * np.sqrt(a2) * np.sqrt(r2 - a2/4)  # (1/2)*base*height

S_regge = 1200*A_eq*eps_eq + 1440*A_rad*eps_rad
S_cont = 32*np.pi**2

print(f"  A_eq  = sqrt(3)/4 * (2-phi) = {A_eq:.6f}")
print(f"  A_rad = sqrt(2-phi)/2 * sqrt(2-(2-phi)/4) = {A_rad:.6f}")
print(f"  S_Regge = {S_regge:.6f}")
print(f"  S_cont  = 32*pi^2 = {S_cont:.6f}")
print(f"  Ratio   = {S_regge/S_cont:.6f}")

# ============================================================
# LAMBDA (from gradient, exact analysis)
# ============================================================
print(f"\n{'='*72}")
print("LAMBDA LATTICE = 2.6489")
print(f"{'='*72}")

Lambda_num = 2.648925  # from exp427 numerical
print(f"  Lambda_lattice = {Lambda_num:.6f}")
print(f"\n  Close to:")
tests = {
    'phi^2 = phi+1': PHI**2,
    '(5+sqrt(5))/2': (5+sqrt5)/2,
    '3*(2-phi)^(-1)/4': 3/(4*(2-PHI)),
    '3*phi^2/4': 3*PHI**2/4,
    'phi^2+1/phi^4': PHI**2 + 1/PHI**4,
    '5*(phi-1)+2': 5*(PHI-1)+2,
    '2+phi/phi^2': 2+PHI/PHI**2,
    '1+phi': 1+PHI,
}
for name, val in sorted(tests.items(), key=lambda x: abs(x[1]-Lambda_num)):
    pct = abs(val-Lambda_num)/Lambda_num*100
    if pct < 5:
        print(f"    {name} = {val:.6f} ({pct:.2f}%)")

# Check: is Lambda = (R_avg_continuum) * (V_regge/V_cont) / 4 ?
# On S^4: Lambda = R/4 = 3. Lattice: Lambda_latt = R_avg_latt/4?
# R_avg = S_regge/V_total. V_total from exp427.

print(f"\n{'='*72}")
print("SUMMARY")
print(f"{'='*72}")
print(f"""
EXACT DIHEDRAL ANGLES OF S^4 DOUBLE-CONE 4-SIMPLEX
===================================================

cos(theta_A) = (2+phi) / (2*(4+phi))  =  1/2 - 1/(4+phi)
             = (10+sqrt(5)) / 38
             = {np.degrees(theta_A):.4f} deg

cos^2(theta_B) = (2-phi) / (8*(4+phi))  =  1/(8*(6*phi+5))
               = (8-3*sqrt(5)) / (8*19)
               = cos = {np.degrees(theta_B):.4f} deg

KEY: denominator 4+phi = phi^2+3 rationalizes to 19 = 4*a1-1
     Same 19 as Z[phi] norm: N(z) = -(n^2+9nt+19t^2)

DEFICIT: eps_eq  = 2pi-4*theta_B = {np.degrees(eps_eq):.4f} deg
         eps_rad = 2pi-5*theta_A = {np.degrees(eps_rad):.4f} deg

S_Regge/S_continuum = {S_regge/S_cont:.4f}
""")

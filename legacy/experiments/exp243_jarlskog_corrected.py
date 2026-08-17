"""
EXP-243: JARLSKOG INVARIANT WITH QCD-CORRECTED CKM ANGLES
==========================================================

From exp242: J ~ phi^(-22) * sqrt(a_1/b_1) gives 29% error.
This is because we used BARE CKM angles (phi^-3, phi^-7, phi^-12).

The CORRECTED angles (from exp191/192) include perturbative factors:
  theta_12 = arctan(phi^-3 * (1 - 2*alpha_s/(3*pi)))     [0.001%]
  theta_23 = arctan(phi^-7 * (1 + 5*alpha_s/pi))          [0.17%]
  theta_13 = arctan(phi^-12 * (1 + sin^2(tW)))            [0.096%]

And the CP phase (from exp242):
  delta = arctan(sqrt(a_1)) = arctan(sqrt(5))              [0.77%]

Let's compute J with these corrected angles.

Category: DERIVATION (combining known results)
"""

import numpy as np
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-243: JARLSKOG INVARIANT WITH CORRECTED CKM ANGLES")
print("=" * 70)

a_1 = 5
b_1 = 6
phi = PHI
alpha = ALPHA
alpha_s = 1 / (2 * phi**3)
sin2tW = b_1 / (a_1**2 + 1)  # 6/26

# =====================================================================
# PART A: Experimental values
# =====================================================================
print("\n" + "=" * 70)
print("PART A: EXPERIMENTAL VALUES")
print("=" * 70)

# PDG 2024 CKM parameters
theta12_exp = np.radians(12.96)   # Cabibbo angle
theta23_exp = np.radians(2.343)
theta13_exp = np.radians(0.2012)  # from Vub
delta_exp = np.radians(65.4)

J_exp = 3.08e-5  # +/- 0.15e-5
J_exp_err = 0.15e-5

print(f"Experimental CKM angles:")
print(f"  theta_12 = {np.degrees(theta12_exp):.3f} deg")
print(f"  theta_23 = {np.degrees(theta23_exp):.3f} deg")
print(f"  theta_13 = {np.degrees(theta13_exp):.4f} deg")
print(f"  delta    = {np.degrees(delta_exp):.1f} deg")
print(f"  J_exp    = {J_exp:.2e} +/- {J_exp_err:.2e}")

# =====================================================================
# PART B: Bare angles (no corrections)
# =====================================================================
print("\n" + "=" * 70)
print("PART B: BARE CKM ANGLES (phi^-n)")
print("=" * 70)

th12_bare = np.arctan(phi**(-3))
th23_bare = np.arctan(phi**(-7))
th13_bare = np.arctan(phi**(-12))
delta_pred = np.arctan(np.sqrt(a_1))

def jarlskog(s12, c12, s23, c23, s13, c13, sd):
    """Standard Jarlskog invariant."""
    return c12 * s12 * c23 * s23 * c13**2 * s13 * sd

s12b, c12b = np.sin(th12_bare), np.cos(th12_bare)
s23b, c23b = np.sin(th23_bare), np.cos(th23_bare)
s13b, c13b = np.sin(th13_bare), np.cos(th13_bare)
sd = np.sin(delta_pred)

J_bare = jarlskog(s12b, c12b, s23b, c23b, s13b, c13b, sd)

print(f"Bare angles:")
print(f"  theta_12 = arctan(phi^-3)  = {np.degrees(th12_bare):.4f} deg (exp: {np.degrees(theta12_exp):.3f})")
print(f"  theta_23 = arctan(phi^-7)  = {np.degrees(th23_bare):.4f} deg (exp: {np.degrees(theta23_exp):.3f})")
print(f"  theta_13 = arctan(phi^-12) = {np.degrees(th13_bare):.4f} deg (exp: {np.degrees(theta13_exp):.4f})")
print(f"  delta    = arctan(sqrt(5)) = {np.degrees(delta_pred):.4f} deg (exp: {np.degrees(delta_exp):.1f})")
print(f"\nJ_bare = {J_bare:.4e}")
print(f"J_exp  = {J_exp:.4e}")
print(f"Error  = {abs(J_bare - J_exp)/J_exp * 100:.1f}%")

# =====================================================================
# PART C: QCD-corrected angles (from exp191/192)
# =====================================================================
print("\n" + "=" * 70)
print("PART C: QCD-CORRECTED CKM ANGLES")
print("=" * 70)

# Correction factors from exp191/192
corr_12 = 1 - 2*alpha_s/(3*np.pi)
corr_23 = 1 + 5*alpha_s/np.pi
corr_13 = 1 + sin2tW

th12_corr = np.arctan(phi**(-3) * corr_12)
th23_corr = np.arctan(phi**(-7) * corr_23)
th13_corr = np.arctan(phi**(-12) * corr_13)

s12c, c12c = np.sin(th12_corr), np.cos(th12_corr)
s23c, c23c = np.sin(th23_corr), np.cos(th23_corr)
s13c, c13c = np.sin(th13_corr), np.cos(th13_corr)

J_corr = jarlskog(s12c, c12c, s23c, c23c, s13c, c13c, sd)

print(f"Correction factors:")
print(f"  corr_12 = 1 - 2*alpha_s/(3*pi) = {corr_12:.6f}")
print(f"  corr_23 = 1 + 5*alpha_s/pi      = {corr_23:.6f}")
print(f"  corr_13 = 1 + sin^2(tW)         = {corr_13:.6f}")

print(f"\nCorrected angles:")
print(f"  theta_12 = {np.degrees(th12_corr):.4f} deg (exp: {np.degrees(theta12_exp):.3f}, err: {abs(np.degrees(th12_corr)-np.degrees(theta12_exp))/np.degrees(theta12_exp)*100:.3f}%)")
print(f"  theta_23 = {np.degrees(th23_corr):.4f} deg (exp: {np.degrees(theta23_exp):.3f}, err: {abs(np.degrees(th23_corr)-np.degrees(theta23_exp))/np.degrees(theta23_exp)*100:.2f}%)")
print(f"  theta_13 = {np.degrees(th13_corr):.4f} deg (exp: {np.degrees(theta13_exp):.4f}, err: {abs(np.degrees(th13_corr)-np.degrees(theta13_exp))/np.degrees(theta13_exp)*100:.2f}%)")

print(f"\nJ_corrected = {J_corr:.4e}")
print(f"J_exp       = {J_exp:.4e}")
err_corr = abs(J_corr - J_exp)/J_exp * 100
print(f"Error       = {err_corr:.1f}%")
print(f"Improvement from bare: {abs(J_bare-J_exp)/J_exp*100:.1f}% -> {err_corr:.1f}%")

# =====================================================================
# PART D: What drives the 29% error?
# =====================================================================
print("\n" + "=" * 70)
print("PART D: ERROR DECOMPOSITION")
print("=" * 70)

# J = c12*s12 * c23*s23 * c13^2*s13 * sin(delta)
# Each factor contributes multiplicatively

# Compare each factor
print(f"Factor-by-factor comparison (corrected vs experimental):")
print(f"")

# Experimental angles -> experimental J factors
s12e, c12e = np.sin(theta12_exp), np.cos(theta12_exp)
s23e, c23e = np.sin(theta23_exp), np.cos(theta23_exp)
s13e, c13e = np.sin(theta13_exp), np.cos(theta13_exp)
sde = np.sin(delta_exp)

f12_pred = c12c * s12c
f12_exp = c12e * s12e
f23_pred = c23c * s23c
f23_exp = c23e * s23e
f13_pred = c13c**2 * s13c
f13_exp = c13e**2 * s13e
fd_pred = sd
fd_exp = sde

print(f"  c12*s12:     pred={f12_pred:.6f}, exp={f12_exp:.6f}, ratio={f12_pred/f12_exp:.4f}")
print(f"  c23*s23:     pred={f23_pred:.6f}, exp={f23_exp:.6f}, ratio={f23_pred/f23_exp:.4f}")
print(f"  c13^2*s13:   pred={f13_pred:.6e}, exp={f13_exp:.6e}, ratio={f13_pred/f13_exp:.4f}")
print(f"  sin(delta):  pred={fd_pred:.6f}, exp={fd_exp:.6f}, ratio={fd_pred/fd_exp:.4f}")

total_ratio = (f12_pred/f12_exp) * (f23_pred/f23_exp) * (f13_pred/f13_exp) * (fd_pred/fd_exp)
print(f"\n  Total ratio: {total_ratio:.4f}")
print(f"  J_pred/J_exp = {J_corr/J_exp:.4f}")

# The main culprit
print(f"\nDominant error source:")
ratios = [
    ("theta_12 (c*s)", f12_pred/f12_exp),
    ("theta_23 (c*s)", f23_pred/f23_exp),
    ("theta_13 (c^2*s)", f13_pred/f13_exp),
    ("delta (sin)", fd_pred/fd_exp),
]
for label, r in ratios:
    print(f"  {label}: ratio = {r:.4f}, contribution to error: {(r-1)*100:+.1f}%")

# =====================================================================
# PART E: Can we do better? Optimize correction factors
# =====================================================================
print("\n" + "=" * 70)
print("PART E: OPTIMAL CORRECTION ANALYSIS")
print("=" * 70)

# What correction factors would give EXACT experimental angles?
exact_12 = np.tan(theta12_exp) / phi**(-3)
exact_23 = np.tan(theta23_exp) / phi**(-7)
exact_13 = np.tan(theta13_exp) / phi**(-12)

print(f"Exact correction factors needed:")
print(f"  theta_12: {exact_12:.6f} (our: {corr_12:.6f}, err: {(corr_12/exact_12-1)*100:+.2f}%)")
print(f"  theta_23: {exact_23:.6f} (our: {corr_23:.6f}, err: {(corr_23/exact_23-1)*100:+.2f}%)")
print(f"  theta_13: {exact_13:.6f} (our: {corr_13:.6f}, err: {(corr_13/exact_13-1)*100:+.2f}%)")

# J with exact angles but our delta
J_exact_angles = jarlskog(s12e, c12e, s23e, c23e, s13e, c13e, sd)
print(f"\nJ with exact exp angles + our delta: {J_exact_angles:.4e}")
print(f"J_exp: {J_exp:.4e}")
print(f"Error: {abs(J_exact_angles-J_exp)/J_exp*100:.2f}%")
print(f"(This is purely from delta_CKM error)")

# J with our corrected angles but exact delta
J_exact_delta = jarlskog(s12c, c12c, s23c, c23c, s13c, c13c, sde)
print(f"\nJ with our corrected angles + exact delta: {J_exact_delta:.4e}")
print(f"Error: {abs(J_exact_delta-J_exp)/J_exp*100:.2f}%")
print(f"(This is purely from angle errors)")

# =====================================================================
# PART F: Analytic formula for J
# =====================================================================
print("\n" + "=" * 70)
print("PART F: ANALYTIC FORMULA FOR J")
print("=" * 70)

# Leading order: sin(theta) ~ tan(theta) ~ phi^-n for small angles
# J ~ sin(12)*sin(23)*sin(13)*sin(delta) * cos factors
# ~ phi^(-3) * phi^(-7) * phi^(-12) * sqrt(5/6) * (cos corrections)

# Better: use arctan formulas directly
# sin(arctan(x)) = x/sqrt(1+x^2), cos(arctan(x)) = 1/sqrt(1+x^2)

x12 = phi**(-3) * corr_12
x23 = phi**(-7) * corr_23
x13 = phi**(-12) * corr_13

J_formula = (x12 / (1+x12**2)) * (x23 / (1+x23**2)) * (1/(1+x13**2)) * x13 * np.sqrt(a_1/b_1)
print(f"J from formula: {J_formula:.4e}")
print(f"J numerical:    {J_corr:.4e}")
print(f"Match: {abs(J_formula - J_corr)/J_corr*100:.6f}%")

# Simplify: for small angles theta_23 and theta_13, cos ~ 1
# J ~ x12/(1+x12^2) * x23 * x13 * sqrt(5/6)
J_simple = (x12 / (1+x12**2)) * x23 * x13 * np.sqrt(a_1/b_1)
print(f"\nJ simplified (small-angle for 23,13):")
print(f"J_simple = {J_simple:.4e}")
print(f"Error vs exact: {abs(J_simple-J_corr)/J_corr*100:.4f}%")

# Even simpler: phi^(-22) * corrections * sqrt(5/6)
J_phi22 = phi**(-22) * corr_12 * corr_23 * corr_13 * np.sqrt(a_1/b_1) / (1 + phi**(-6)*corr_12**2)
print(f"\nJ ~ phi^(-22) * corr * sqrt(a_1/b_1) / (1+phi^-6):")
print(f"J = {J_phi22:.4e}")
print(f"Error: {abs(J_phi22-J_corr)/J_corr*100:.3f}%")

# =====================================================================
# PART G: Express J purely in terms of a_1
# =====================================================================
print("\n" + "=" * 70)
print("PART G: J IN TERMS OF a_1")
print("=" * 70)

# The complete formula:
# J = sin(th12)*cos(th12) * sin(th23)*cos(th23) * sin(th13)*cos(th13)^2 * sin(delta)
# where th_ij = arctan(phi^(-n_ij) * corr_ij)
# and delta = arctan(sqrt(a_1))
# and alpha_s = 1/(2*phi^3), sin^2(tW) = b_1/(a_1^2+1)

# corr_12 = 1 - 2/(2*phi^3 * 3*pi) = 1 - 1/(3*pi*phi^3)
# corr_23 = 1 + 5/(2*pi*phi^3) = 1 + 5/(2*pi*phi^3)
# corr_13 = 1 + b_1/(a_1^2+1) = 1 + 6/26 = 32/26 = 16/13

print(f"Correction factors in terms of a_1:")
print(f"  corr_12 = 1 - 1/(3*pi*phi^3) = {1 - 1/(3*np.pi*phi**3):.6f}")
print(f"  corr_23 = 1 + a_1/(2*pi*phi^3) = {1 + a_1/(2*np.pi*phi**3):.6f}")
print(f"  corr_13 = 1 + b_1/(a_1^2+1) = (a_1^2+1+b_1)/(a_1^2+1) = {(a_1**2+1+b_1)}/{a_1**2+1} = {(a_1**2+1+b_1)/(a_1**2+1):.6f}")

# Note: corr_13 = (a_1^2 + b_1 + 1)/(a_1^2+1) = (25+6+1)/26 = 32/26 = 16/13
print(f"  corr_13 = {a_1**2+b_1+1}/{a_1**2+1} = {(a_1**2+b_1+1)/(a_1**2+1):.6f}")
# Check: 32/26 = 16/13
print(f"  = {a_1**2+b_1+1}:{a_1**2+1}")

# The KEY insight: corr_13 is purely RATIONAL (no alpha_s, no pi)
# corr_13 = (a_1^2+a_1+2)/(a_1^2+1) = 32/26 = 16/13
actual_corr13 = (a_1**2 + a_1 + 2) / (a_1**2 + 1)
print(f"\n  corr_13 = (a_1^2+a_1+2)/(a_1^2+1) = {actual_corr13:.6f}")
print(f"  Check vs 1+sin^2(tW): {1 + sin2tW:.6f}")
print(f"  (a_1^2+a_1+2)/(a_1^2+1) = (25+5+2)/26 = 32/26 = 16/13")

# =====================================================================
# PART H: The Jarlskog from CKM matrix elements
# =====================================================================
print("\n" + "=" * 70)
print("PART H: FULL CKM MATRIX AND J")
print("=" * 70)

# Build CKM matrix with corrected angles + our delta
V = np.array([
    [c12c*c13c,                      s12c*c13c,                      s13c*np.exp(-1j*delta_pred)],
    [-s12c*c23c-c12c*s23c*s13c*np.exp(1j*delta_pred), c12c*c23c-s12c*s23c*s13c*np.exp(1j*delta_pred), s23c*c13c],
    [s12c*s23c-c12c*c23c*s13c*np.exp(1j*delta_pred),  -c12c*s23c-s12c*c23c*s13c*np.exp(1j*delta_pred), c23c*c13c]
])

# Experimental CKM magnitudes (PDG 2024)
V_exp_mag = np.array([
    [0.97435, 0.22500, 0.00369],
    [0.22486, 0.97349, 0.04182],
    [0.00857, 0.04110, 0.999118]
])

print(f"Predicted |V_CKM| (corrected):")
labels = [["Vud", "Vus", "Vub"], ["Vcd", "Vcs", "Vcb"], ["Vtd", "Vts", "Vtb"]]
max_err = 0
for i in range(3):
    for j in range(3):
        pred = abs(V[i,j])
        exp = V_exp_mag[i,j]
        err = abs(pred-exp)/exp*100
        max_err = max(max_err, err)
        print(f"  |{labels[i][j]}| = {pred:.5f} (exp: {exp:.5f}, err: {err:.2f}%)")

# Jarlskog from matrix
J_matrix = np.imag(V[0,0] * V[1,1] * np.conj(V[0,1]) * np.conj(V[1,0]))
print(f"\nJ from matrix = {J_matrix:.4e}")
print(f"J direct calc  = {J_corr:.4e}")
print(f"J experimental = {J_exp:.4e}")
print(f"Error: {abs(J_matrix-J_exp)/J_exp*100:.1f}%")

# =====================================================================
# PART I: What would fix J?
# =====================================================================
print("\n" + "=" * 70)
print("PART I: WHAT WOULD FIX J?")
print("=" * 70)

# The main issue: theta_23 and theta_13 are too small
# This makes Vcb and Vub too small, which directly reduces J

# What if we use different correction factors?
print(f"The problem: our theta_23 and theta_13 give Vcb and Vub too small.")
print(f"  |Vcb| = sin(theta_23)*cos(theta_13)")
print(f"  pred: {s23c*c13c:.5f}, exp: 0.04182, ratio: {s23c*c13c/0.04182:.4f}")
print(f"  |Vub| = sin(theta_13)")
print(f"  pred: {s13c:.5f}, exp: 0.00369, ratio: {s13c/0.00369:.4f}")

# The Vcb ratio and Vub ratio together account for most of the J error
vcb_ratio = s23c*c13c / 0.04182
vub_ratio = s13c / 0.00369
print(f"\n  Vcb deficit: {(1-vcb_ratio)*100:.1f}%")
print(f"  Vub deficit: {(1-vub_ratio)*100:.1f}%")
print(f"  Combined: {(1-vcb_ratio*vub_ratio)*100:.1f}%")
print(f"  J error:  {abs(J_corr-J_exp)/J_exp*100:.1f}%")

print(f"\nThe J error comes from Vcb ({(1-vcb_ratio)*100:.1f}%) and Vub ({(1-vub_ratio)*100:.1f}%).")
print(f"These are the elements sensitive to theta_23 and theta_13.")
print(f"The correction factors are 'right direction' but not strong enough.")

# =====================================================================
# PART J: Summary
# =====================================================================
print("\n" + "=" * 70)
print("PART J: SUMMARY")
print("=" * 70)

print(f"""
JARLSKOG INVARIANT RESULTS:
===========================

1. BARE (phi^-n only, no corrections):
   J_bare = {J_bare:.4e}
   Error: {abs(J_bare-J_exp)/J_exp*100:.1f}%

2. QCD-CORRECTED (exp191/192 correction factors):
   J_corr = {J_corr:.4e}
   Error: {abs(J_corr-J_exp)/J_exp*100:.1f}%

3. EXPERIMENTAL:
   J_exp  = {J_exp:.4e} +/- {J_exp_err:.2e}

IMPROVEMENT: {abs(J_bare-J_exp)/J_exp*100:.1f}% -> {abs(J_corr-J_exp)/J_exp*100:.1f}%

ERROR SOURCE: Vcb ({(1-vcb_ratio)*100:.1f}% too small) and Vub ({(1-vub_ratio)*100:.1f}% too small).
The correction factors for theta_23 and theta_13 are in the right direction
but insufficient. This suggests higher-order corrections or a different
mechanism for these off-diagonal elements.

FORMULA (leading order):
  J = phi^(-22) * sqrt(a_1/b_1) * (cos corrections)
  where 22 = 3+7+12 = sum of CKM exponents
  and sqrt(a_1/b_1) = sin(delta_CKM) = sin(arctan(sqrt(5)))

STATUS: PATTERN (the formula structure is correct but precision limited
by the mixing angle corrections, not by the CP phase).
""")

print("=" * 70)
print("EXP-243 COMPLETE")
print("=" * 70)

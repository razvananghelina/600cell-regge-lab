"""
Electron anomalous magnetic moment (a_e) computed from:
  1. Framework algebraic alpha (from 600-cell theory)
  2. CODATA 2022 alpha
Compared with experimental a_e.

QED perturbation theory through 5 loops + hadronic + electroweak.
"""

import math

print("=" * 72)
print("  ELECTRON g-2 FROM FRAMEWORK ALPHA vs CODATA ALPHA")
print("=" * 72)

# =====================================================================
# 1. Framework alpha from: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
# =====================================================================
a1 = 5
phi = (1 + math.sqrt(5)) / 2

A_coeff = 2 * math.pi
B_coeff = -4 * a1 * phi**4
C_coeff = 1

disc = B_coeff**2 - 4 * A_coeff * C_coeff
sqrt_disc = math.sqrt(disc)

root1 = (-B_coeff + sqrt_disc) / (2 * A_coeff)
root2 = (-B_coeff - sqrt_disc) / (2 * A_coeff)

# Pick root closer to 1/137
target = 1.0 / 137.0
if abs(root1 - target) < abs(root2 - target):
    alpha_fw = root1
    alpha_other = root2
else:
    alpha_fw = root2
    alpha_other = root1

print("\n--- STEP 1: Framework Alpha ---")
print("  Equation: 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0")
print("  a1 = %d, phi = (1+sqrt(5))/2 = %.15f" % (a1, phi))
print("  Coefficients: A = 2*pi = %.15f" % A_coeff)
print("                B = -4*a1*phi^4 = %.15f" % B_coeff)
print("                C = 1")
print("  Discriminant = %.15f" % disc)
print("  Root 1 = %.15e" % root1)
print("  Root 2 = %.15e" % root2)
print("  Selected (closer to 1/137): alpha_fw = %.15e" % alpha_fw)
print("  1/alpha_fw = %.10f" % (1/alpha_fw))

# =====================================================================
# 2. CODATA 2022 alpha
# =====================================================================
alpha_codata = 1.0 / 137.035999177
alpha_codata_inv = 137.035999177

print("\n--- STEP 2: CODATA 2022 Alpha ---")
print("  alpha_CODATA = 1/137.035999177 = %.15e" % alpha_codata)
print("  1/alpha_CODATA = %.10f" % alpha_codata_inv)

# Comparison
diff_alpha = alpha_fw - alpha_codata
diff_ppm = diff_alpha / alpha_codata * 1e6
print("\n  alpha_fw - alpha_CODATA = %.6e" % diff_alpha)
print("  Relative difference    = %.4f ppm" % diff_ppm)
print("  1/alpha_fw - 1/alpha_CODATA = %.10f" % (1/alpha_fw - alpha_codata_inv))

# =====================================================================
# 3. QED perturbation theory for a_e
# =====================================================================
print("\n--- STEP 3: QED Coefficients ---")

C1 = 0.5
C2 = -0.328478965579193
C3 = 1.181241456587
C4 = -1.9113
C5 = 6.737

print("  C1 = %s" % C1)
print("  C2 = %s" % C2)
print("  C3 = %s" % C3)
print("  C4 = %s (uncertainty: 0.0018)" % C4)
print("  C5 = %s (uncertainty: 0.159)" % C5)

a_e_had = 1.87e-12
a_e_ew  = 0.030e-12

print("\n  a_e(hadronic)     = %.3e" % a_e_had)
print("  a_e(electroweak)  = %.3e" % a_e_ew)

def compute_ae(alpha):
    x = alpha / math.pi
    ae_qed = 0.0
    terms = []
    coeffs = [C1, C2, C3, C4, C5]
    for n in range(1, 6):
        term = coeffs[n-1] * x**n
        ae_qed += term
        terms.append((n, coeffs[n-1], term))
    ae_total = ae_qed + a_e_had + a_e_ew
    return ae_qed, ae_total, terms

ae_qed_fw, ae_tot_fw, terms_fw = compute_ae(alpha_fw)
ae_qed_co, ae_tot_co, terms_co = compute_ae(alpha_codata)

print("\n--- QED Series with Framework Alpha ---")
print("  alpha/pi = %.15e" % (alpha_fw/math.pi))
for n, Cn, term in terms_fw:
    print("  %d-loop: C%d*(alpha/pi)^%d = %+.15e" % (n, n, n, term))
print("  a_e(QED)   = %.15e" % ae_qed_fw)
print("  a_e(had)   = %.15e" % a_e_had)
print("  a_e(EW)    = %.15e" % a_e_ew)
print("  a_e(total) = %.15e" % ae_tot_fw)

print("\n--- QED Series with CODATA Alpha ---")
print("  alpha/pi = %.15e" % (alpha_codata/math.pi))
for n, Cn, term in terms_co:
    print("  %d-loop: C%d*(alpha/pi)^%d = %+.15e" % (n, n, n, term))
print("  a_e(QED)   = %.15e" % ae_qed_co)
print("  a_e(had)   = %.15e" % a_e_had)
print("  a_e(EW)    = %.15e" % a_e_ew)
print("  a_e(total) = %.15e" % ae_tot_co)

# =====================================================================
# 4. Experimental value
# =====================================================================
ae_exp = 0.00115965218059

print("\n--- STEP 4: Experimental a_e ---")
print("  a_e(exp) = %.14e" % ae_exp)
print("  Uncertainty: 1.3 x 10^-13")

# =====================================================================
# 5. Comparison
# =====================================================================
print("\n" + "=" * 72)
print("  FINAL COMPARISON")
print("=" * 72)

diff_fw = ae_tot_fw - ae_exp
diff_co = ae_tot_co - ae_exp

ppt_fw = diff_fw * 1e12
ppt_co = diff_co * 1e12

rel_fw = diff_fw / ae_exp * 1e12
rel_co = diff_co / ae_exp * 1e12

print("\n  1/alpha_fw     = %.10f" % (1/alpha_fw))
print("  1/alpha_CODATA = %.10f" % alpha_codata_inv)
print("  Difference     = %+.10f" % (1/alpha_fw - alpha_codata_inv))
print("  Relative       = %+.4f ppm" % diff_ppm)

print("\n  a_e(framework)  = %.15e" % ae_tot_fw)
print("  a_e(CODATA)     = %.15e" % ae_tot_co)
print("  a_e(experiment) = %.15e" % ae_exp)

print("\n  --- Absolute differences (x 10^-12) ---")
print("  a_e(fw)     - a_e(exp) = %+.4f x 10^-12" % ppt_fw)
print("  a_e(CODATA) - a_e(exp) = %+.4f x 10^-12" % ppt_co)

print("\n  --- Relative differences (parts per 10^12 of a_e) ---")
print("  [a_e(fw) - a_e(exp)] / a_e(exp)     = %+.4f ppt" % rel_fw)
print("  [a_e(CODATA) - a_e(exp)] / a_e(exp) = %+.4f ppt" % rel_co)

print("\n  --- In units of experimental uncertainty (1.3 x 10^-13) ---")
exp_unc = 1.3e-13
sigma_fw = abs(diff_fw) / exp_unc
sigma_co = abs(diff_co) / exp_unc
print("  |a_e(fw) - a_e(exp)| / sigma_exp     = %.2f sigma" % sigma_fw)
print("  |a_e(CODATA) - a_e(exp)| / sigma_exp = %.2f sigma" % sigma_co)

print("\n  *** VERDICT ***")
if abs(diff_fw) < abs(diff_co):
    print("  Framework alpha gives BETTER a_e prediction!")
    print("  Improvement factor: %.2fx closer to experiment" % (abs(diff_co)/abs(diff_fw)))
elif abs(diff_fw) > abs(diff_co):
    print("  CODATA alpha gives BETTER a_e prediction.")
    print("  CODATA is %.2fx closer to experiment" % (abs(diff_fw)/abs(diff_co)))
else:
    print("  Both give identical predictions.")

# Additional context
print("\n  --- Theory Uncertainty Budget ---")
c4_unc_ae = 0.0018 * (alpha_fw/math.pi)**4
c5_unc_ae = 0.159 * (alpha_fw/math.pi)**5
had_unc = 0.02e-12
print("  C4 uncertainty contribution to a_e: %.3e" % c4_unc_ae)
print("  C5 uncertainty contribution to a_e: %.3e" % c5_unc_ae)
print("  Hadronic uncertainty:                %.3e" % had_unc)
total_thy_unc = math.sqrt(c4_unc_ae**2 + c5_unc_ae**2 + had_unc**2)
print("  Total theory uncertainty (quadrature): %.3e" % total_thy_unc)

print("\n" + "=" * 72)
print("  SENSITIVITY ANALYSIS: da_e/dalpha")
print("=" * 72)

x_fw = alpha_fw / math.pi
dae_dalpha = (C1/math.pi + 2*C2*x_fw/math.pi + 3*C3*x_fw**2/math.pi
              + 4*C4*x_fw**3/math.pi + 5*C5*x_fw**4/math.pi)
print("\n  da_e/dalpha = %.10f" % dae_dalpha)
print("  Dominated by Schwinger: 1/(2*pi) = %.10f" % (1/(2*math.pi)))

# What alpha would perfectly match experiment?
alpha_star = alpha_codata
coeffs = [C1, C2, C3, C4, C5]
for _ in range(50):
    x = alpha_star / math.pi
    ae_val = sum(coeffs[n-1] * x**n for n in range(1, 6)) + a_e_had + a_e_ew
    dae = sum(n * coeffs[n-1] * x**(n-1) / math.pi for n in range(1, 6))
    alpha_star -= (ae_val - ae_exp) / dae

print("\n  Alpha that perfectly reproduces a_e(exp):")
print("  1/alpha* = %.10f" % (1/alpha_star))
print("  vs 1/alpha_fw     = %.10f  (diff: %+.10f)" % (1/alpha_fw, 1/alpha_fw - 1/alpha_star))
print("  vs 1/alpha_CODATA = %.10f  (diff: %+.10f)" % (alpha_codata_inv, alpha_codata_inv - 1/alpha_star))

print("\n" + "=" * 72)

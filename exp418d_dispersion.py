"""
EXP-418d: Dispersion Relation - Massless vs Massive Graviton
==============================================================
Test: Are coexact eigenvalues proportional to l(l+2)?

If lambda_coexact(l) = c * l(l+2) for constant c => MASSLESS graviton
If lambda_coexact(l) = c * l(l+2) + m^2  => MASSIVE graviton
If neither => non-trivial dispersion (lattice effects)

We already have the level assignment from exp418c:
  Scalar:  l=0: 0, l=1: 2.292, l=2: 5.528, l=3: 9.000, l=4: 12.000, l=5: 14.000
  Coexact: l=1: 0.528, l=2: 1.146, l=3: 1.933, l=4: 2.807, l=5: (3.556+3.697)/2

Also compare with scalar eigenvalues: if both scale the same way,
the tensor/scalar ratio of eigenvalues should be constant.
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-418d: Dispersion Relation Test")
print("=" * 70)

# From exp418c - exact eigenvalues (verified numerically)
# Scalar Delta_0 eigenvalues by level
scalar = {
    0: 0.0,
    1: 12 - 6*PHI,      # 2.29180
    2: 12 - 4*PHI,      # 5.52786
    3: 9.0,
    4: 12.0,
    5: 14.0,
}

# Coexact eigenvalues by level (single eigenvalue per level for l=1-4)
# l=5 splits into two: 30 modes at 3.5563 + 40 modes at 3.6972
coexact = {
    1: 7 - 4*PHI,       # 0.52786
    2: 6 - 3*PHI,       # 1.14590
    3: 44 - 26*PHI,     # 1.93319  (or -15+11*phi for the a+b*phi form)
    4: -15 + 11*PHI,    # 2.80742
    5: None,            # split: will use weighted average
}

# l=5 weighted average: (30*3.5563 + 40*3.6972) / 70
l5_avg = (30*(11*PHI - 11) + 40*(49 - 28*PHI)) / 70
# Actually let me just use the numerical values
l5_a = -11 + 9*PHI   # 3.55634
l5_b = 49 - 28*PHI   # 3.69722
l5_weighted = (30*l5_a + 40*l5_b) / 70
coexact[5] = l5_weighted

print("\n--- Raw data ---")
print(f"{'l':>3s} {'l(l+2)':>8s} {'lam_scalar':>12s} {'lam_coexact':>12s}")
for l in range(0, 6):
    ll2 = l*(l+2)
    ls = scalar[l]
    lc = coexact.get(l, None)
    lc_str = f"{lc:12.5f}" if lc is not None else "        ---"
    print(f"{l:3d} {ll2:8d} {ls:12.5f} {lc_str}")

# ============================================================
# Test 1: lambda_coexact = c * l(l+2) ?
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: lambda_coexact(l) / l(l+2) = const? (massless)")
print("=" * 70)

print(f"\n{'l':>3s} {'l(l+2)':>8s} {'lam_co':>10s} {'lam/l(l+2)':>12s}")
ratios_co = []
for l in range(1, 6):
    ll2 = l*(l+2)
    lc = coexact[l]
    r = lc / ll2
    ratios_co.append(r)
    print(f"{l:3d} {ll2:8d} {lc:10.5f} {r:12.6f}")

mean_r = np.mean(ratios_co)
std_r = np.std(ratios_co)
print(f"\n  Mean ratio: {mean_r:.6f}")
print(f"  Std:        {std_r:.6f}")
print(f"  Variation:  {std_r/mean_r*100:.1f}%")

if std_r/mean_r < 0.05:
    print("  => MASSLESS (constant ratio to < 5%)")
elif std_r/mean_r < 0.20:
    print("  => APPROXIMATELY MASSLESS (mild lattice distortion)")
else:
    print("  => NOT proportional to l(l+2)")

# ============================================================
# Test 2: lambda_scalar = c * l(l+2) ?
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: lambda_scalar(l) / l(l+2) = const? (scalar comparison)")
print("=" * 70)

print(f"\n{'l':>3s} {'l(l+2)':>8s} {'lam_sc':>10s} {'lam/l(l+2)':>12s}")
ratios_sc = []
for l in range(1, 6):
    ll2 = l*(l+2)
    ls = scalar[l]
    r = ls / ll2
    ratios_sc.append(r)
    print(f"{l:3d} {ll2:8d} {ls:10.5f} {r:12.6f}")

mean_sc = np.mean(ratios_sc)
std_sc = np.std(ratios_sc)
print(f"\n  Mean ratio: {mean_sc:.6f}")
print(f"  Std:        {std_sc:.6f}")
print(f"  Variation:  {std_sc/mean_sc*100:.1f}%")

# ============================================================
# Test 3: lambda_coexact / lambda_scalar = const?
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: lambda_coexact / lambda_scalar = const?")
print("(Same dispersion relation => same effective speed)")
print("=" * 70)

print(f"\n{'l':>3s} {'lam_sc':>10s} {'lam_co':>10s} {'ratio co/sc':>12s}")
ratios_cs = []
for l in range(1, 6):
    ls = scalar[l]
    lc = coexact[l]
    r = lc / ls
    ratios_cs.append(r)
    print(f"{l:3d} {ls:10.5f} {lc:10.5f} {r:12.6f}")

mean_cs = np.mean(ratios_cs)
std_cs = np.std(ratios_cs)
print(f"\n  Mean ratio: {mean_cs:.6f}")
print(f"  Std:        {std_cs:.6f}")
print(f"  Variation:  {std_cs/mean_cs*100:.1f}%")

# ============================================================
# Test 4: Fit lambda = a * l(l+2) + b (massive dispersion)
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Fit lambda = a * l(l+2) + b")
print("(b=0 => massless, b>0 => massive with m^2=b)")
print("=" * 70)

# Coexact
x_co = np.array([l*(l+2) for l in range(1,6)], dtype=float)
y_co = np.array([coexact[l] for l in range(1,6)])
A_co = np.vstack([x_co, np.ones(len(x_co))]).T
result_co = np.linalg.lstsq(A_co, y_co, rcond=None)
a_co, b_co = result_co[0]
residuals_co = y_co - (a_co * x_co + b_co)
rms_co = np.sqrt(np.mean(residuals_co**2))

print(f"\n  COEXACT: lambda = {a_co:.6f} * l(l+2) + ({b_co:.6f})")
print(f"  Effective speed^2: a = {a_co:.6f}")
print(f"  Mass^2 term:       b = {b_co:.6f}")
print(f"  RMS residual:      {rms_co:.6f}")
print(f"  Residuals: {np.round(residuals_co, 5)}")

# Scalar
y_sc = np.array([scalar[l] for l in range(1,6)])
result_sc = np.linalg.lstsq(A_co, y_sc, rcond=None)
a_sc, b_sc = result_sc[0]
residuals_sc = y_sc - (a_sc * x_co + b_sc)
rms_sc = np.sqrt(np.mean(residuals_sc**2))

print(f"\n  SCALAR:  lambda = {a_sc:.6f} * l(l+2) + ({b_sc:.6f})")
print(f"  Effective speed^2: a = {a_sc:.6f}")
print(f"  Mass^2 term:       b = {b_sc:.6f}")
print(f"  RMS residual:      {rms_sc:.6f}")

# ============================================================
# Test 5: Fit lambda = a * l(l+2) only (force massless)
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: Fit lambda = a * l(l+2) (forced massless)")
print("=" * 70)

a_co_0 = np.dot(x_co, y_co) / np.dot(x_co, x_co)
res_co_0 = y_co - a_co_0 * x_co
rms_co_0 = np.sqrt(np.mean(res_co_0**2))

a_sc_0 = np.dot(x_co, y_sc) / np.dot(x_co, x_co)
res_sc_0 = y_sc - a_sc_0 * x_co
rms_sc_0 = np.sqrt(np.mean(res_sc_0**2))

print(f"\n  COEXACT: lambda = {a_co_0:.6f} * l(l+2)")
print(f"  RMS residual:      {rms_co_0:.6f}")
print(f"  Residuals: {np.round(res_co_0, 5)}")

print(f"\n  SCALAR:  lambda = {a_sc_0:.6f} * l(l+2)")
print(f"  RMS residual:      {rms_sc_0:.6f}")
print(f"  Residuals: {np.round(res_sc_0, 5)}")

print(f"\n  Speed ratio (coexact/scalar): {a_co_0/a_sc_0:.6f}")

# ============================================================
# Test 6: Non-linear fit lambda = a * [l(l+2)]^p
# ============================================================
print("\n" + "=" * 70)
print("TEST 6: Power-law lambda = a * [l(l+2)]^p")
print("=" * 70)

log_x = np.log(x_co)
log_y_co = np.log(y_co)
log_y_sc = np.log(y_sc)

p_co, log_a_co = np.polyfit(log_x, log_y_co, 1)
p_sc, log_a_sc = np.polyfit(log_x, log_y_sc, 1)

print(f"\n  COEXACT: lambda = {np.exp(log_a_co):.6f} * l(l+2)^{p_co:.4f}")
print(f"  Power: {p_co:.4f} (1.0 = massless linear dispersion)")

print(f"\n  SCALAR:  lambda = {np.exp(log_a_sc):.6f} * l(l+2)^{p_sc:.4f}")
print(f"  Power: {p_sc:.4f}")

# ============================================================
# Test 7: Is there a universal effective radius R?
# ============================================================
print("\n" + "=" * 70)
print("TEST 7: Effective S^3 radius")
print("=" * 70)

# On S^3 of radius R: lambda = l(l+2)/R^2
# Our discrete eigenvalues should give R^2 = l(l+2)/lambda

print(f"\n  {'l':>3s} {'l(l+2)':>8s} {'R^2_scalar':>12s} {'R^2_coexact':>12s}")
for l in range(1, 6):
    ll2 = l*(l+2)
    R2_sc = ll2 / scalar[l]
    R2_co = ll2 / coexact[l]
    print(f"  {l:3d} {ll2:8d} {R2_sc:12.5f} {R2_co:12.5f}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
DISPERSION RELATION ANALYSIS:

1. Massless test (lambda/l(l+2) = const):
   Coexact: variation {std_r/mean_r*100:.1f}% around mean {mean_r:.4f}
   Scalar:  variation {std_sc/mean_sc*100:.1f}% around mean {mean_sc:.4f}

2. Linear fit lambda = a*l(l+2) + b:
   Coexact: a={a_co:.5f}, b={b_co:.5f}, RMS={rms_co:.5f}
   Scalar:  a={a_sc:.5f}, b={b_sc:.5f}, RMS={rms_sc:.5f}

3. Forced massless lambda = a*l(l+2):
   Coexact: a={a_co_0:.5f}, RMS={rms_co_0:.5f}
   Scalar:  a={a_sc_0:.5f}, RMS={rms_sc_0:.5f}
   Speed ratio: {a_co_0/a_sc_0:.5f}

4. Power-law lambda ~ l(l+2)^p:
   Coexact: p={p_co:.4f}
   Scalar:  p={p_sc:.4f}
""")

# Verdict
if abs(b_co) < 0.05 and rms_co < 0.05:
    print("VERDICT: MASSLESS graviton (b ~ 0, good linear fit)")
elif abs(b_co) < abs(a_co) and rms_co < 0.1:
    print("VERDICT: APPROXIMATELY MASSLESS (small mass term, dominated by linear)")
    print(f"  Effective mass^2 = {b_co:.4f} (in lattice units)")
else:
    print("VERDICT: NEITHER massless nor simply massive")
    print("  Non-trivial dispersion (lattice effects dominant)")
    if abs(p_co - 1) < 0.15:
        print(f"  But power-law exponent p={p_co:.3f} is close to 1 (within {abs(p_co-1)*100:.0f}%)")
        print(f"  This suggests approximate linear dispersion with lattice corrections")
        print(f"  In the continuum limit (l << l_max), the graviton is effectively massless")

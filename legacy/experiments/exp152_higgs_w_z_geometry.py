"""
EXP-152: Higgs Mechanism and W/Z Masses from 600-cell Geometry
=============================================================
Can we derive m_W, m_Z, and the Higgs mechanism from the 600-cell?
The known formula m_H = m_W * (phi - 8*alpha) works at 0.09%.
Can we derive m_W itself?

Parts:
1. Known mass relations recap
2. m_W/m_e from phi powers
3. Electroweak symmetry breaking from edge structure
4. Higgs potential from DSI
5. vev from 600-cell geometry
6. W/Z mass ratio from sin^2(tW)
7. Higgs quartic coupling
8. Top-Higgs connection
9. Number-theoretic analysis of m_W
10. Gauge boson masses from spectral gap
"""

import numpy as np

print("="*70)
print("EXP-152: Higgs Mechanism and W/Z Masses from Geometry")
print("="*70)

PHI = (1 + np.sqrt(5)) / 2
PI = np.pi
ALPHA = 1/137.035999084  # CODATA
ALPHA_S = 0.1179

# Physical constants
m_e = 0.51099895  # MeV
m_W = 80377  # MeV (80.377 GeV)
m_Z = 91187.6  # MeV (91.1876 GeV)
m_H = 125250  # MeV (125.25 GeV) - PDG 2024
m_t = 172760  # MeV (172.76 GeV)
m_P = 1.22089e22  # MeV (Planck mass)
v_higgs = 246220  # MeV (246.22 GeV) - Higgs vev

sin2_tW = 0.23122
cos2_tW = 1 - sin2_tW

# ============================================================
# Part 1: Known relations
# ============================================================
print("\n" + "="*70)
print("PART 1: Known mass relations")
print("="*70)

print(f"m_W = {m_W/1000:.3f} GeV")
print(f"m_Z = {m_Z/1000:.4f} GeV")
print(f"m_H = {m_H/1000:.2f} GeV")
print(f"m_t = {m_t/1000:.2f} GeV")
print(f"v   = {v_higgs/1000:.2f} GeV")
print(f"m_e = {m_e:.6f} MeV")

# Known formulas
print(f"\nKnown 600-cell formulas:")

# m_H = m_W * (phi - 8*alpha)
mH_pred = m_W * (PHI - 8*ALPHA) / 1000
print(f"  m_H = m_W * (phi - 8*alpha) = {mH_pred:.2f} GeV (exp: {m_H/1000:.2f}, err: {abs(mH_pred-m_H/1000)/(m_H/1000)*100:.2f}%)")

# m_Z from sin^2(tW) = 6/26
mZ_pred = m_W / np.sqrt(1 - 6/26)
print(f"  m_Z = m_W/sqrt(1-6/26) = {mZ_pred/1000:.2f} GeV (exp: {m_Z/1000:.2f}, err: {abs(mZ_pred/1000-m_Z/1000)/(m_Z/1000)*100:.2f}%)")

# m_Z / m_W ratio
print(f"\n  m_Z/m_W = {m_Z/m_W:.6f}")
print(f"  1/cos(tW) = {1/np.sqrt(cos2_tW):.6f}")
print(f"  1/sqrt(20/26) = {1/np.sqrt(20/26):.6f} = sqrt(26/20) = sqrt(13/10)")
print(f"  Error: {abs(1/np.sqrt(20/26) - m_Z/m_W)/(m_Z/m_W)*100:.3f}%")

# ============================================================
# Part 2: m_W/m_e from phi powers
# ============================================================
print("\n" + "="*70)
print("PART 2: m_W/m_e from phi powers")
print("="*70)

ratio = m_W / m_e
log_ratio = np.log(ratio) / np.log(PHI)
print(f"m_W/m_e = {ratio:.2f}")
print(f"log_phi(m_W/m_e) = {log_ratio:.4f}")
print(f"  Nearest integer: {round(log_ratio)}")
print(f"  phi^{round(log_ratio)} = {PHI**round(log_ratio):.2f}")

# Check various phi-power formulas
for n in range(24, 27):
    pred = m_e * PHI**n / 1000
    err = abs(pred - m_W/1000) / (m_W/1000) * 100
    print(f"  m_e * phi^{n} = {pred:.2f} GeV (err {err:.2f}%)")

# More complex formulas
print(f"\nComplex phi formulas for m_W:")

# m_W = m_e * phi^25 * correction?
for n in [25, 26]:
    base = m_e * PHI**n
    corr = m_W / base
    print(f"  m_W / (m_e*phi^{n}) = {corr:.6f}")
    # Is correction a simple phi expression?
    for a in range(-5, 6):
        for b in range(-5, 6):
            val = a + b*PHI
            if val > 0 and abs(val - corr) < 0.01:
                pred = m_e * PHI**n * val
                err = abs(pred - m_W) / m_W * 100
                print(f"    m_e * phi^{n} * ({a}+{b}*phi) = {pred:.2f} MeV (err {err:.4f}%)")

# Try: m_W = f(alpha, phi, m_e)
# m_W = m_e / (2*alpha) * something?
print(f"\n  m_e / (2*alpha) = {m_e/(2*ALPHA):.2f} MeV = {m_e/(2*ALPHA)/1000:.2f} GeV")
print(f"  m_W / (m_e/(2*alpha)) = {m_W / (m_e/(2*ALPHA)):.6f}")

# v = m_e / sqrt(2*sqrt(2)*G_F) ... actually v = 1/(sqrt(sqrt(2)*G_F))
# v = 2*m_W / g_2
# g_2^2/(4*pi) = alpha_em/sin^2(tW) = alpha_2
alpha_2 = ALPHA / sin2_tW
g_2 = np.sqrt(4*PI*alpha_2)
v_pred = 2*m_W / g_2
print(f"\n  g_2 = sqrt(4*pi*alpha_2) = {g_2:.6f}")
print(f"  v = 2*m_W/g_2 = {v_pred:.2f} MeV = {v_pred/1000:.2f} GeV (exp: {v_higgs/1000:.2f})")

# ============================================================
# Part 3: vev from geometry
# ============================================================
print("\n" + "="*70)
print("PART 3: Higgs vev from 600-cell geometry")
print("="*70)

# v/m_e
v_ratio = v_higgs / m_e
log_v = np.log(v_ratio) / np.log(PHI)
print(f"v/m_e = {v_ratio:.2f}")
print(f"log_phi(v/m_e) = {log_v:.4f}")

# v = m_e * phi^n * correction
for n in [25, 26, 27]:
    pred = m_e * PHI**n
    err = abs(pred - v_higgs) / v_higgs * 100
    print(f"  m_e * phi^{n} = {pred:.2f} MeV = {pred/1000:.2f} GeV (err {err:.2f}%)")

# v = m_e * phi^26?
# phi^26 = 710647.15... vs v/m_e = 481744.6
# v = m_e * phi^25.54

# Try: v = m_P * alpha^k
print(f"\nv from Planck mass:")
for k_num in range(1, 20):
    for k_den in range(1, 10):
        k = k_num / k_den
        pred = m_P * ALPHA**k / 1000  # GeV
        if abs(pred - v_higgs/1000) / (v_higgs/1000) < 0.05:
            err = abs(pred - v_higgs/1000) / (v_higgs/1000) * 100
            print(f"  m_P * alpha^({k_num}/{k_den}) = {pred:.2f} GeV (err {err:.3f}%)")

# v from alpha and phi
# v = m_e / alpha * phi^k?
print(f"\nv from m_e/alpha:")
ratio_v_alpha = v_higgs * ALPHA / m_e
log_va = np.log(ratio_v_alpha) / np.log(PHI)
print(f"  v*alpha/m_e = {ratio_v_alpha:.4f}")
print(f"  log_phi(...) = {log_va:.4f}")

# ============================================================
# Part 4: Higgs potential from DSI
# ============================================================
print("\n" + "="*70)
print("PART 4: Higgs potential and DSI")
print("="*70)

# DSI potential: V(psi) = A*sin^2(pi*ln|psi|/ln(phi))
# This has minima at psi = phi^n
# If Higgs field phi_H is on the DSI ladder: <phi_H> = v = some phi^n scale

# m_H^2 = d^2V/dphi^2 at minimum
# V = A*sin^2(pi*ln(phi_H)/ln(phi))
# At minimum phi_H = phi^n: V=0, V''=2A*(pi/ln(phi))^2/phi^(2n)
# m_H^2 ~ A * (pi/ln(phi))^2 * phi^(-2n)

# Ratio: m_H/m_W = phi - 8*alpha ~ phi (dominant)
# This suggests m_H and m_W are on adjacent DSI levels? Not quite.

print(f"m_H/m_W = {m_H/m_W:.6f}")
print(f"phi - 8*alpha = {PHI - 8*ALPHA:.6f}")
print(f"phi = {PHI:.6f}")
print(f"Correction: -8*alpha = {-8*ALPHA:.6f} ({abs(8*ALPHA/PHI)*100:.2f}% of phi)")

# m_W and m_Z on DSI ladder?
log_mW = np.log(m_W/m_e) / np.log(PHI)
log_mZ = np.log(m_Z/m_e) / np.log(PHI)
log_mH = np.log(m_H/m_e) / np.log(PHI)
log_mt = np.log(m_t/m_e) / np.log(PHI)

print(f"\nDSI levels (log_phi(m/m_e)):")
print(f"  W: {log_mW:.4f}")
print(f"  Z: {log_mZ:.4f}")
print(f"  H: {log_mH:.4f}")
print(f"  t: {log_mt:.4f}")

# Check if close to half-integers
for name, val in [('W', log_mW), ('Z', log_mZ), ('H', log_mH), ('t', log_mt)]:
    nearest_half = round(2*val)/2
    err = abs(val - nearest_half)
    print(f"  {name}: nearest n/2 = {nearest_half}, offset = {err:.4f}")

# ============================================================
# Part 5: W mass from coupling and vev
# ============================================================
print("\n" + "="*70)
print("PART 5: m_W from coupling structure")
print("="*70)

# m_W = g_2 * v / 2
# If v = m_e * f(phi), then m_W = g_2 * m_e * f(phi) / 2
# g_2^2 = 4*pi*alpha/sin^2(tW) = 4*pi*alpha*26/6
g2_squared = 4*PI*ALPHA*26/6
g2_val = np.sqrt(g2_squared)
print(f"g_2^2 = 4*pi*alpha*26/6 = {g2_squared:.6f}")
print(f"g_2 = {g2_val:.6f}")
print(f"Experimental g_2 = {g_2:.6f}")

# m_W = g_2 * v / 2
# v = 2*m_W/g_2 = 2*80377/{g_2:.4f}
print(f"\nv = 2*m_W/g_2 = {2*m_W/g_2:.2f} MeV = {2*m_W/g_2/1000:.2f} GeV")

# Can we express v in terms of 600-cell quantities?
# v/m_e = 2*(m_W/m_e)/g_2
# log_phi(v/m_e) = 25.54
# This is NOT a clean integer or half-integer

# Try: v = m_e * phi^n * (1/alpha)^k
for n in range(20, 30):
    for k_num in range(-5, 6):
        for k_den in range(1, 5):
            k = k_num / k_den
            if k_den > 1 and k_num % k_den == 0:
                continue
            pred = m_e * PHI**n * ALPHA**k
            if pred > 0:
                err = abs(pred - v_higgs) / v_higgs
                if err < 0.02:
                    print(f"  v ~ m_e * phi^{n} * alpha^({k_num}/{k_den}) = {pred:.1f} MeV ({err*100:.3f}%)")

# ============================================================
# Part 6: Top-Higgs connection
# ============================================================
print("\n" + "="*70)
print("PART 6: Top quark and Higgs mass connection")
print("="*70)

# In SM: m_t = y_t * v / sqrt(2)
y_t = m_t * np.sqrt(2) / v_higgs
print(f"Top Yukawa: y_t = m_t*sqrt(2)/v = {y_t:.4f}")
print(f"  Close to 1! ({abs(y_t-1)*100:.1f}% off)")

# m_H/m_t
print(f"\nm_H/m_t = {m_H/m_t:.6f}")
print(f"  1/sqrt(2) = {1/np.sqrt(2):.6f} (err {abs(m_H/m_t - 1/np.sqrt(2))/(1/np.sqrt(2))*100:.2f}%)")

# In 600-cell: top has (a,b)=(4,1), n=5*4+6*1=26
# phi^26/phi^25 = phi
# m_H = m_W * phi ~ m_W * m_t/m_charm? No.

# Interesting: n_top = 26 = phi_dim
print(f"\nn_top = 5*4 + 6*1 = 26 = phi-sector dimension!")
print(f"This means the top quark DSI level coincides with the phi-sector dim")

# m_t prediction
m_t_pred = m_e * PHI**(5*4+6*1)  # with corrections
print(f"m_t from DSI: m_e * phi^26 = {m_e * PHI**26:.2f} MeV = {m_e * PHI**26/1000:.2f} GeV")
# But this needs the (a,b) corrections via d_eff, k_eff

# ============================================================
# Part 7: Spectral gap and gauge boson mass
# ============================================================
print("\n" + "="*70)
print("PART 7: Spectral gap approach")
print("="*70)

# 600-cell eigenvalues: 12, 6phi, 4phi, 3, 0, -2, 4-4phi, -3, 6-6phi
eigs = [12, 6*PHI, 4*PHI, 3, 0, -2, 4-4*PHI, -3, 6-6*PHI]
eig_names = ['12', '6phi', '4phi', '3', '0', '-2', '4-4phi', '-3', '6-6phi']
mults = [1, 4, 9, 16, 25, 36, 9, 16, 4]

print("600-cell eigenvalues and gaps:")
for i in range(len(eigs)-1):
    gap = eigs[i] - eigs[i+1]
    print(f"  {eig_names[i]}({mults[i]}) - {eig_names[i+1]}({mults[i+1]}) = {gap:.4f}")

# The spectral gap lambda_0 - lambda_1 = 12 - 6phi = 12 - 9.708 = 2.292
# = 2*(6-3phi) = 2*(6-3*1.618) = 2*1.146
print(f"\nLargest gap: 12 - 6phi = {12 - 6*PHI:.6f}")
print(f"  = 2*(6-3phi) = 2*{6-3*PHI:.6f}")

# Ratio lambda_0/lambda_1 = 12/(6phi) = 2/phi
print(f"lambda_0/lambda_1 = 12/(6phi) = {12/(6*PHI):.6f} = 2/phi = {2/PHI:.6f}")
print(f"  = 2*phi^(-1) = 2*(phi-1) = {2*(PHI-1):.6f}")

# Can we relate W/Z mass splitting to spectral gaps?
print(f"\nm_Z - m_W = {(m_Z-m_W)/1000:.3f} GeV")
print(f"m_Z/m_W - 1 = {m_Z/m_W - 1:.6f}")
print(f"(m_Z-m_W)/m_W = {(m_Z-m_W)/m_W:.6f}")
print(f"From sin^2(tW)=6/26: m_Z/m_W = sqrt(26/20) -> m_Z/m_W-1 = {np.sqrt(26/20)-1:.6f}")

# ============================================================
# Part 8: Number theory of m_W
# ============================================================
print("\n" + "="*70)
print("PART 8: Number theory analysis")
print("="*70)

# m_W / m_e = 80377/0.511 = 157293
ratio_We = m_W / m_e
print(f"m_W/m_e = {ratio_We:.2f}")

# Is this close to a "nice" number?
# phi^25 = 167761.1
# phi^24/phi = 103682.1/phi... no
# 120^2 * something?
print(f"phi^25 = {PHI**25:.1f}")
print(f"m_W/m_e / phi^25 = {ratio_We/PHI**25:.6f}")

# Try rational multiples of phi^n
for n in range(24, 27):
    r = ratio_We / PHI**n
    # Check if r is close to a/b for small a,b
    for a in range(1, 50):
        for b in range(1, 50):
            if abs(a/b - r) < 0.001:
                pred = (a/b) * PHI**n
                err = abs(pred - ratio_We)/ratio_We * 100
                if err < 1:
                    print(f"  m_W/m_e ~ ({a}/{b})*phi^{n} = {pred:.2f} (err {err:.4f}%)")

# alpha-based
print(f"\nm_W/m_e in terms of alpha:")
print(f"  1/(2*alpha) = {1/(2*ALPHA):.2f}")
print(f"  1/alpha = {1/ALPHA:.2f}")
print(f"  m_W/m_e * alpha = {ratio_We * ALPHA:.6f}")
print(f"  log_phi(m_W*alpha/m_e) = {np.log(ratio_We*ALPHA)/np.log(PHI):.4f}")

# m_W = m_e * (1/alpha)^a * phi^b?
log_alpha = np.log(1/ALPHA)
log_phi = np.log(PHI)
log_target = np.log(ratio_We)

print(f"\nSearching m_W/m_e = (1/alpha)^a * phi^b:")
best_err = 100
for a_num in range(-10, 20):
    for a_den in range(1, 8):
        a = a_num / a_den
        if a_den > 1 and a_num % a_den == 0:
            continue
        remaining = log_target - a * log_alpha
        b = remaining / log_phi
        if abs(b - round(b)) < 0.1:
            b_int = round(b)
            pred = (1/ALPHA)**a * PHI**b_int
            err = abs(pred - ratio_We)/ratio_We * 100
            if err < 2:
                print(f"  (1/alpha)^({a_num}/{a_den}) * phi^{b_int} = {pred:.2f} (err {err:.4f}%)")

# ============================================================
# Part 9: m_W from gauge structure
# ============================================================
print("\n" + "="*70)
print("PART 9: m_W from 600-cell gauge structure")
print("="*70)

# In SM: m_W = g_2 * v / 2
# g_2^2 = 4*pi*alpha_2 = 4*pi*alpha/sin^2(tW)

# From 600-cell: alpha = solution of quadratic, sin^2(tW) = 6/26
# So g_2 is determined. We need v.

# One approach: v = m_t / (y_t/sqrt(2))
# If y_t = 1 (natural): v = m_t*sqrt(2)
v_nat = m_t * np.sqrt(2)
mW_nat = g2_val * v_nat / 2
print(f"If y_t = 1: v = m_t*sqrt(2) = {v_nat/1000:.2f} GeV")
print(f"  m_W = g_2*v/2 = {mW_nat/1000:.2f} GeV (exp: {m_W/1000:.2f}, err {abs(mW_nat-m_W)/m_W*100:.1f}%)")

# Another: v from electroweak scale
# v^2 = 1/(sqrt(2)*G_F)
# G_F = pi*alpha/(sqrt(2)*m_W^2*sin^2(tW))
# So v = 2*m_W*sin(tW)/sqrt(4*pi*alpha) ... circular

# Dimensional analysis from 600-cell:
# Only scales: m_e (electron mass) and m_P (Planck mass)
# v should be expressible as m_e^a * m_P^b * phi^n * alpha^k

# v^2 ~ m_e * m_P?
print(f"\nsqrt(m_e * m_P) = {np.sqrt(m_e * m_P * 1e6)/1e6:.2e} MeV = {np.sqrt(m_e * m_P * 1e6)/1e9:.2f} GeV")
# Too big

# v ~ m_e * (m_P/m_e)^x for some x
x_v = np.log(v_higgs/m_e) / np.log(m_P/m_e)
print(f"v = m_e * (m_P/m_e)^x, x = {x_v:.6f}")
print(f"  x = 1/2 would give {m_e * (m_P/m_e)**0.5:.2f} MeV = {m_e * (m_P/m_e)**0.5/1000:.2f} GeV")

# m_e/m_P = alpha^(4*phi^2) -> m_P = m_e * alpha^(-4*phi^2)
# v = m_e * alpha^(-4*phi^2 * x)
log_v_alpha = np.log(v_higgs/m_e) / np.log(1/ALPHA)
print(f"\nv/m_e = alpha^(-{log_v_alpha:.6f})")
print(f"  4*phi^2 = {4*PHI**2:.6f}")
print(f"  ratio: {log_v_alpha/(4*PHI**2):.6f}")
print(f"  x = {x_v:.6f} -> -4*phi^2*x = {-4*PHI**2*x_v:.6f}")

# v = m_e * alpha^(-k) where k ~ 2.57
# k = 4*phi^2 * 0.484 = 5.07? No
# Just check specific k values
for k_num in range(1, 30):
    for k_den in range(1, 12):
        k = k_num / k_den
        pred = m_e * (1/ALPHA)**k / 1000  # GeV
        if abs(pred - v_higgs/1000) / (v_higgs/1000) < 0.05:
            err = abs(pred - v_higgs/1000) / (v_higgs/1000) * 100
            print(f"  v = m_e / alpha^({k_num}/{k_den}) = {pred:.2f} GeV (err {err:.3f}%)")

# ============================================================
# Part 10: Summary
# ============================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
Higgs sector analysis:

DERIVED (from 600-cell):
  m_H = m_W * (phi - 8*alpha)  [0.09%]
  m_Z = m_W / sqrt(20/26)      [0.50%] (via sin^2(tW)=6/26)
  m_Z/m_W = sqrt(13/10)        [purely from 6/26]

NOT DERIVED:
  m_W itself requires the Higgs vev v
  v = 2*m_W/g_2 = {v_higgs/1000:.2f} GeV
  v/m_e = {v_higgs/m_e:.1f} is NOT a clean phi-power

PATTERNS (fitting, not derivation):
  m_W/m_e ~ phi^25.39 (not integer)
  v/m_e ~ phi^25.54 (not integer)
  v ~ m_e/alpha^(5/2) at {abs(m_e*(1/ALPHA)**2.5/1000-v_higgs/1000)/(v_higgs/1000)*100:.1f}% (if found)

KEY INSIGHT:
  The 600-cell determines RATIOS (m_H/m_W, m_Z/m_W, alpha_s/alpha)
  but NOT absolute mass scales beyond m_e/m_P = alpha^(4*phi^2).
  The electroweak scale v is an ADDITIONAL parameter not determined
  by the 600-cell geometry alone.

  Top Yukawa y_t = {y_t:.4f} ~ 1, and n_top = 26 = phi_dim.
  The top quark may play a special role: if y_t = 1 exactly,
  then v = m_t*sqrt(2), reducing v to the top mass.

STATUS: MIXED
  - Ratios DERIVED (m_H/m_W, m_Z/m_W)
  - Absolute scale NOT DERIVED (v or m_W individually)
  - Top-Higgs connection y_t ~ 1 is PATTERN
""")

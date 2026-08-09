"""
EXP-189: Correction Patterns - Deep Dive
=========================================
Following exp188 findings:
- Leptons: delta_d/delta_k ~ -phi (Galois direction)
- Up quarks: delta_k ~ alpha_s (almost exact)
Key question: is there a DERIVABLE pattern?
"""

import math

PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
SIN2_TW = 0.23122
ALPHA_S = 0.1179

# PDG masses in MeV
M = {'e': 0.51099895, 'mu': 105.6583755, 'tau': 1776.86,
     'u': 2.16, 'c': 1270.0, 't': 172690.0,
     'd': 4.67, 's': 93.4, 'b': 4180.0}

def lp(x):
    return math.log(abs(x)) / math.log(PHI)

print("=" * 70)
print("EXP-189: CORRECTION PATTERNS")
print("=" * 70)

# Exact sector d_eff, k_eff (from exp188)
# Leptons: (mu-e): Da=(1,1), (tau-e): Da=(1,2) => d=5.214, k=5.865
exp_mu = lp(M['mu']/M['e'])
exp_tau = lp(M['tau']/M['e'])
d_lep = exp_mu - (exp_tau - exp_mu)  # d = exp_mu - k
k_lep = exp_tau - exp_mu
dd_lep = d_lep - 5
dk_lep = k_lep - 6

# Up quarks: u(3,-2)->c(2,1): Da=(-1,3), u(3,-2)->t(4,1): Da=(1,3)
exp_c = lp(M['c']/M['u'])
exp_t = lp(M['t']/M['u'])
det_up = (-1)*3 - 1*3
d_up = (exp_c*3 - exp_t*3) / det_up
k_up = ((-1)*exp_t - 1*exp_c) / det_up
dd_up = d_up - 5
dk_up = k_up - 6

# Down quarks: d(1,0)->s(1,1): Da=(0,1), d(1,0)->b(-1,4): Da=(-2,4)
exp_s = lp(M['s']/M['d'])
exp_b = lp(M['b']/M['d'])
det_dn = 0*4 - (-2)*1
d_dn = (exp_s*4 - exp_b*1) / det_dn
k_dn = (0*exp_b - (-2)*exp_s) / det_dn
dd_dn = d_dn - 5
dk_dn = k_dn - 6

print("\n=== SECTOR CORRECTIONS ===")
print(f"{'Sector':12} {'dd':>10} {'dk':>10} {'dd/dk':>10} {'dd+dk':>10}")
print("-" * 55)
print(f"{'Leptons':12} {dd_lep:10.6f} {dk_lep:10.6f} {dd_lep/dk_lep:10.4f} {dd_lep+dk_lep:10.6f}")
print(f"{'Up quarks':12} {dd_up:10.6f} {dk_up:10.6f} {dd_up/dk_up:10.4f} {dd_up+dk_up:10.6f}")
print(f"{'Down quarks':12} {dd_dn:10.6f} {dk_dn:10.6f} {dd_dn/dk_dn:10.4f} {dd_dn+dk_dn:10.6f}")

# =====================================================================
# TEST 1: Up quarks and alpha_s
# =====================================================================
print("\n" + "=" * 70)
print("TEST 1: Up quarks and alpha_s")
print("=" * 70)

print(f"\n  delta_k(up) = {dk_up:.8f}")
print(f"  alpha_s     = {ALPHA_S:.8f}")
print(f"  Ratio: {dk_up/ALPHA_S:.8f}")
print(f"  => delta_k(up) = alpha_s * {dk_up/ALPHA_S:.6f}")

print(f"\n  delta_d(up) = {dd_up:.8f}")
print(f"  alpha_s*(1-alpha_s) = {ALPHA_S*(1-ALPHA_S):.8f}")
print(f"  Ratio: {dd_up/(ALPHA_S*(1-ALPHA_S)):.8f}")

# More precise tests
tests_up = {
    'alpha_s':                      (ALPHA_S, dk_up),
    'alpha_s^2':                    (ALPHA_S**2, dd_up),
    'alpha_s*(1-alpha_s)':          (ALPHA_S*(1-ALPHA_S), dd_up),
    'alpha_s/phi':                  (ALPHA_S/PHI, dd_up),
    'alpha_s*sin^2(tW)':           (ALPHA_S*SIN2_TW, dd_up),
    'alpha_s*(1-sin^2(tW))':       (ALPHA_S*(1-SIN2_TW), dd_up),
    'alpha_s*(phi-1)':             (ALPHA_S*(PHI-1), dd_up),
    'alpha_s/phi^2':               (ALPHA_S/PHI**2, dd_up),
}

print(f"\n  Expressions for up quark corrections:")
print(f"  {'Expression':30} {'Value':>12} {'Target':>12} {'Err%':>8}")
print(f"  {'-'*65}")
for name, (val, target) in tests_up.items():
    err = abs(val - target) / abs(target) * 100
    marker = " <--" if err < 3 else ""
    print(f"  {name:30} {val:12.8f} {target:12.8f} {err:8.2f}%{marker}")

# =====================================================================
# TEST 2: Lepton corrections and Galois direction
# =====================================================================
print("\n" + "=" * 70)
print("TEST 2: Lepton corrections - Galois direction analysis")
print("=" * 70)

# If correction is along (phi, -1): (delta_d, delta_k) = epsilon * (phi, -1)
# Then epsilon = delta_k / (-1) = -delta_k
epsilon_from_k = -dk_lep
epsilon_from_d = dd_lep / PHI

print(f"\n  If (dd, dk) = epsilon * (phi, -1):")
print(f"    epsilon from dk: {epsilon_from_k:.8f}")
print(f"    epsilon from dd: {epsilon_from_d:.8f}")
print(f"    Average: {(epsilon_from_k + epsilon_from_d)/2:.8f}")
print(f"    Discrepancy: {abs(epsilon_from_k - epsilon_from_d):.8f}")
print(f"    Relative: {abs(epsilon_from_k - epsilon_from_d)/epsilon_from_k*100:.2f}%")

# What is epsilon in terms of known constants?
eps = (epsilon_from_k + epsilon_from_d) / 2
print(f"\n  epsilon = {eps:.8f}")
print(f"  Known quantities:")
print(f"    alpha = {ALPHA:.8f}, eps/alpha = {eps/ALPHA:.4f}")
print(f"    sin^2(tW) = {SIN2_TW:.8f}, eps/sin^2(tW) = {eps/SIN2_TW:.4f}")
print(f"    1/phi^4 = {1/PHI**4:.8f}, eps/(1/phi^4) = {eps*PHI**4:.4f}")
print(f"    1/phi^3 = {1/PHI**3:.8f}, eps/(1/phi^3) = {eps*PHI**3:.4f}")
print(f"    alpha_s = {ALPHA_S:.8f}, eps/alpha_s = {eps/ALPHA_S:.4f}")

# eps ~ 0.1325. Test specific:
print(f"\n  Testing expressions for epsilon:")
tests_eps = {
    '1/phi^4':                PHI**(-4),
    'alpha_s':                ALPHA_S,
    'sin^2(tW)/phi':         SIN2_TW/PHI,
    'alpha_s*phi^0.3':       ALPHA_S*PHI**0.3,
    '3/(4*phi^3)':           3/(4*PHI**3),
    '1/(2*phi^3)':           1/(2*PHI**3),
    'alpha_s/phi^0.5':       ALPHA_S/PHI**0.5,
    'phi^(-4)*phi^(0.1)':    PHI**(-3.9),
    '1/(phi^3+1)':           1/(PHI**3+1),
}
for name, val in tests_eps.items():
    err = abs(val - eps) / eps * 100
    marker = " <--" if err < 5 else ""
    print(f"    {name:25} = {val:.8f}, err = {err:.2f}%{marker}")

# =====================================================================
# TEST 3: Correction in Z[phi] coordinates
# =====================================================================
print("\n" + "=" * 70)
print("TEST 3: Corrections in algebraic coordinates z = d + k*phi")
print("=" * 70)

# For each sector, compute the correction to z and z'
for name, dd, dk in [('Leptons', dd_lep, dk_lep),
                      ('Up', dd_up, dk_up),
                      ('Down', dd_dn, dk_dn)]:
    dz = dd + dk * PHI
    dz_prime = dd + dk * (1 - PHI)  # phi' = 1-phi = -1/phi
    print(f"\n  {name}:")
    print(f"    delta_z  = dd + dk*phi  = {dd:.6f} + {dk:.6f}*{PHI:.4f} = {dz:.6f}")
    print(f"    delta_z' = dd + dk*phi' = {dd:.6f} + {dk:.6f}*{1-PHI:.4f} = {dz_prime:.6f}")
    print(f"    |delta_z|  = {abs(dz):.6f}")
    print(f"    |delta_z'| = {abs(dz_prime):.6f}")
    print(f"    Ratio |dz'|/|dz| = {abs(dz_prime)/abs(dz):.4f}" if abs(dz) > 1e-10 else "    dz ~ 0!")

# =====================================================================
# TEST 4: Quadratic form analysis
# =====================================================================
print("\n" + "=" * 70)
print("TEST 4: Norm of correction vector")
print("=" * 70)

# The correction (dd, dk) in Z[phi] has norm N(dz) = dz * dz'
# = (dd + dk*phi)(dd + dk*phi') = dd^2 + dd*dk - dk^2
for name, dd, dk in [('Leptons', dd_lep, dk_lep),
                      ('Up', dd_up, dk_up),
                      ('Down', dd_dn, dk_dn)]:
    norm = dd**2 + dd*dk - dk**2
    print(f"  {name:12}: N(dz) = dd^2 + dd*dk - dk^2 = {norm:.8f}")

# =====================================================================
# TEST 5: Can inter-sector ratios be predicted?
# =====================================================================
print("\n" + "=" * 70)
print("TEST 5: Inter-sector correction ratios")
print("=" * 70)

# Ratios of corrections
print(f"\n  dd(down)/dd(up) = {dd_dn/dd_up:.4f}")
print(f"  dk(down)/dk(up) = {dk_dn/dk_up:.4f}")
print(f"  dd(down)/dd(lep) = {dd_dn/dd_lep:.4f}")
print(f"  dk(down)/dk(lep) = {dk_dn/dk_lep:.4f}")

# Electric charges: Q(up) = +2/3, Q(down) = -1/3, Q(e) = -1
# Color factor: quarks have Casimir C_F = 4/3
# Maybe corrections scale with Q^2 or C_F?
print(f"\n  Charge-squared ratios:")
print(f"  Q(e)^2 / Q(up)^2 = 1 / (4/9) = {9/4:.4f}")
print(f"  Q(e)^2 / Q(down)^2 = 1 / (1/9) = {9:.4f}")
print(f"  Q(up)^2 / Q(down)^2 = (4/9) / (1/9) = {4:.4f}")

print(f"\n  dd ratios / charge-squared ratios:")
print(f"  dd(up)/dd(lep) = {dd_up/dd_lep:.4f} vs Q(up)^2/Q(e)^2 = {4/9:.4f}")
print(f"  dd(dn)/dd(lep) = {dd_dn/dd_lep:.4f} vs Q(down)^2/Q(e)^2 = {1/9:.4f}")

# That doesn't work. Try total gauge correction factor
# Leptons: only alpha. Quarks: alpha + alpha_s * C_F
C_F = 4.0/3
print(f"\n  Gauge correction factors:")
print(f"  Leptons: f = alpha = {ALPHA:.6f}")
print(f"  Up quarks: f = alpha*(2/3)^2 + alpha_s*C_F = {ALPHA*4/9 + ALPHA_S*C_F:.6f}")
print(f"  Down quarks: f = alpha*(1/3)^2 + alpha_s*C_F = {ALPHA/9 + ALPHA_S*C_F:.6f}")

f_lep = ALPHA
f_up = ALPHA*4/9 + ALPHA_S*C_F
f_dn = ALPHA/9 + ALPHA_S*C_F

print(f"\n  dd scaled by gauge factor:")
print(f"  dd(lep)/f(lep) = {dd_lep/f_lep:.4f}")
print(f"  dd(up)/f(up) = {dd_up/f_up:.4f}")
print(f"  dd(dn)/f(dn) = {dd_dn/f_dn:.4f}")
print(f"\n  dk scaled by gauge factor:")
print(f"  dk(lep)/f(lep) = {dk_lep/f_lep:.4f}")
print(f"  dk(up)/f(up) = {dk_up/f_up:.4f}")
print(f"  dk(dn)/f(dn) = {dk_dn/f_dn:.4f}")

# =====================================================================
# TEST 6: Mass anomalous dimensions
# =====================================================================
print("\n" + "=" * 70)
print("TEST 6: Anomalous dimensions interpretation")
print("=" * 70)

print("""
In QFT, masses run with energy scale:
  m(mu_R) = m(mu_0) * (alpha_s(mu_R)/alpha_s(mu_0))^(gamma_m / beta_0)
where gamma_m is the mass anomalous dimension.

For QCD at 1-loop: gamma_m = 6*C_F/(16*pi^2) = 8/(16*pi^2)
And beta_0 = (11*N_c - 2*N_f)/(16*pi^2)

The CORRECTIONS to integer exponents could be mass running effects.
If the 600-cell prediction is for a high-energy scale (Planck?),
then the corrections are the running to low energy.
""")

# 1-loop running: m(low)/m(high) = (alpha_s(low)/alpha_s(high))^(gamma_m/2*beta_0)
# gamma_m = 8/(16*pi^2) for QCD
# beta_0 = (33-2*N_f)/(12*pi) for SU(3)
N_f = 6
beta_0 = (33 - 2*N_f) / (12*math.pi)
gamma_m = 8 / (16*math.pi**2)
running_exp = gamma_m / (2 * beta_0)

print(f"  QCD mass anomalous dimension exponent:")
print(f"  gamma_m / (2*beta_0) = {running_exp:.6f}")
print(f"  = {gamma_m:.6f} / (2 * {beta_0:.6f})")

# alpha_s runs from ~0.12 (M_Z) to ~0.3 (1 GeV)
# (alpha_s(1GeV)/alpha_s(M_Z))^running_exp
alpha_s_low = 0.35  # ~ at 1 GeV
ratio_alphas = alpha_s_low / ALPHA_S
correction_from_running = running_exp * math.log(ratio_alphas)

print(f"\n  Running from M_Z to 1 GeV:")
print(f"  alpha_s(1GeV) ~ {alpha_s_low}")
print(f"  Running correction to log(m): {correction_from_running:.6f}")
print(f"  Compare with quark corrections dk(up) = {dk_up:.6f}, dk(dn) = {dk_dn:.6f}")

# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
FINDINGS:

1. LEPTON GALOIS DIRECTION:
   delta_d / delta_k = {dd_lep/dk_lep:.4f} ~ -phi ({-PHI:.4f})
   => Correction is along Galois conjugate direction (phi, -1)
   => delta_z = dd + dk*phi ~ 0 (correction vanishes in phi-sector!)
   => delta_z' = dd + dk*phi' = {dd_lep + dk_lep*(1-PHI):.6f} (all correction in phi'-sector)
   STATUS: PATTERN (suggestive but not derived)

2. UP QUARK AND alpha_s:
   delta_k(up) = {dk_up:.6f} ~ alpha_s = {ALPHA_S:.6f} (ratio {dk_up/ALPHA_S:.4f})
   delta_d(up) = {dd_up:.6f} ~ alpha_s*(1-alpha_s) = {ALPHA_S*(1-ALPHA_S):.6f} (ratio {dd_up/(ALPHA_S*(1-ALPHA_S)):.4f})
   STATUS: INTERESTING but NOT proven

3. DOWN QUARK CORRECTIONS LARGER:
   dd(dn) = {dd_dn:.4f}, dk(dn) = {dk_dn:.4f} (both larger than up quarks)
   dd(dn)/dd(up) = {dd_dn/dd_up:.2f}, dk(dn)/dk(up) = {dk_dn/dk_up:.2f}
   Not a simple charge-squared ratio. Possibly mass running effect.

4. GAUGE SCALING FAILS:
   Corrections do NOT simply scale with the gauge coupling strength.
   dd/f is NOT constant across sectors.
   => Corrections are NOT purely perturbative gauge effects.

5. MASS RUNNING:
   QCD anomalous dimension gives correction ~ 0.04 (from M_Z to 1 GeV)
   This is of the right ORDER but not the right VALUE for dk(up) = 0.12

BOTTOM LINE:
- Integer levels n = 5a + 6b: DERIVED
- Base phi: DERIVED
- Corrections to d and k: PATTERN at best
- Specific delta_d and delta_k values: NOT DERIVED
- Sector dependence suggests gauge corrections, but no clean formula found
""")

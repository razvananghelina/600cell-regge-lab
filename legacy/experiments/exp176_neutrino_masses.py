"""
EXP-176: Neutrino Masses from 600-Cell Theory
==============================================
Open problem: neutrino masses are tiny (sub-eV) while our mass formula
gives MeV-GeV scale masses. How can neutrinos fit?

Previous: exp116 was NEGATIVE (direct application of mass formula fails).

Approaches:
A. Seesaw mechanism with 600-cell parameters
B. Conjugate sector (phi' instead of phi)
C. Different lattice line or negative exponents
D. Mass-squared differences from phi-tower
E. Connection to non-unit fermions and generation structure
F. DSI (discrete scale invariance) in neutrino sector

Experimental data:
- Delta m^2_21 = 7.53e-5 eV^2 (solar)
- Delta m^2_32 = 2.453e-3 eV^2 (atmospheric, normal ordering)
- Sum m_nu < 0.12 eV (cosmological bound, Planck 2018)
- m_nu_e < 0.8 eV (KATRIN)

Category: RESEARCH (open problem, honest about status)
"""

import numpy as np
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2
ALPHA = 1/137.036
M_E = 0.000511  # GeV (electron mass)

print("=" * 70)
print("EXP-176: Neutrino Masses from 600-Cell Theory")
print("=" * 70)

# Experimental data
dm2_21 = 7.53e-5  # eV^2
dm2_32 = 2.453e-3  # eV^2 (normal ordering)
sum_nu_max = 0.12  # eV (Planck upper bound)

# Estimated individual masses (normal ordering, minimal)
# m1 ~ 0, m2 ~ sqrt(dm2_21) ~ 0.0087 eV, m3 ~ sqrt(dm2_32) ~ 0.0495 eV
m1_est = 0.0  # lightest (could be ~0)
m2_est = np.sqrt(dm2_21)
m3_est = np.sqrt(dm2_32)

print(f"\n  Experimental constraints:")
print(f"  Delta m^2_21 = {dm2_21:.2e} eV^2")
print(f"  Delta m^2_32 = {dm2_32:.2e} eV^2")
print(f"  Sum m_nu < {sum_nu_max} eV")
print(f"")
print(f"  Estimated masses (normal ordering, minimal):")
print(f"  m1 ~ {m1_est:.4f} eV (could be 0)")
print(f"  m2 ~ {m2_est:.4f} eV = sqrt(dm2_21)")
print(f"  m3 ~ {m3_est:.4f} eV = sqrt(dm2_32)")

# ============================================================
# APPROACH A: Seesaw Mechanism
# ============================================================
print("\n" + "=" * 70)
print("APPROACH A: Seesaw Mechanism with 600-Cell Parameters")
print("=" * 70)

print(f"""
  Type-I Seesaw: m_nu ~ m_D^2 / M_R

  where m_D = Dirac mass (from phi-tower) and M_R = heavy right-handed mass.

  In our theory, the natural scale for M_R could be:
  1. The "maximum mass" from the phi-tower: m_e * phi^n_max
  2. The E8 scale or Planck scale
  3. A scale derived from the 600-cell geometry

  Testing with m_D = charged lepton masses (same generation):
""")

# Seesaw with different M_R scales
for M_R_label, M_R_GeV in [
    ("phi^30 * m_e", M_E * PHI**30),
    ("phi^31 * m_e", M_E * PHI**31),
    ("GUT scale (2e16 GeV)", 2e16),
    ("Planck (1.22e19 GeV)", 1.22e19),
    ("120 * m_top (600-cell)", 120 * 172.76),
]:
    print(f"\n  M_R = {M_R_label} = {M_R_GeV:.3e} GeV:")
    m_D = [M_E, 0.10566, 1.777]  # e, mu, tau masses in GeV
    names = ['nu_e', 'nu_mu', 'nu_tau']
    for i, (name, mD) in enumerate(zip(names, m_D)):
        m_nu = mD**2 / M_R_GeV  # GeV
        m_nu_eV = m_nu * 1e9  # convert to eV
        print(f"    {name}: m_D={mD:.4e} GeV -> m_nu = {m_nu_eV:.4e} eV", end="")
        # Compare to experimental
        if i == 1:
            print(f"  (exp ~ {m2_est:.4e} eV, ratio = {m_nu_eV/m2_est:.2f})")
        elif i == 2:
            print(f"  (exp ~ {m3_est:.4f} eV, ratio = {m_nu_eV/m3_est:.2f})")
        else:
            print()

# Find the M_R that gives correct neutrino masses
print(f"\n  Finding M_R for correct neutrino masses:")
for name, mD, m_exp in [('nu_mu', 0.10566, m2_est*1e-9), ('nu_tau', 1.777, m3_est*1e-9)]:
    if m_exp > 0:
        M_R_needed = mD**2 / m_exp
        log_phi = np.log(M_R_needed/M_E) / np.log(PHI)
        print(f"  {name}: M_R = {M_R_needed:.3e} GeV = m_e * phi^{log_phi:.1f}")

# ============================================================
# APPROACH B: Conjugate Sector (phi' = -1/phi)
# ============================================================
print("\n" + "=" * 70)
print("APPROACH B: Conjugate Sector Masses")
print("=" * 70)

print(f"""
  Instead of m = m_e * phi^n, use the GALOIS CONJUGATE:
  m' = m_e * |phi'|^n = m_e * (1/phi)^n = m_e * phi^(-n)

  For charged leptons: n_e=0, n_mu=11, n_tau=17
  Conjugate masses: m_e * phi^(-n)
""")

print(f"  Conjugate sector masses (m_e * phi^(-n)):")
for name, n, m_charged in [('e', 0, M_E), ('mu', 11, 0.10566), ('tau', 17, 1.777)]:
    m_conj = M_E * PHI**(-n)  # GeV
    m_conj_eV = m_conj * 1e9  # eV
    print(f"  {name}: n={n:>2}, m_charged = {m_charged:.4e} GeV, m_conj = {m_conj_eV:.4e} eV")

# Compare to experimental neutrino masses
print(f"\n  Comparison with experimental neutrino masses:")
print(f"  m_e * phi^(-11) = {M_E * PHI**(-11) * 1e9:.4e} eV (vs m_nu2 ~ {m2_est:.4e} eV)")
print(f"  m_e * phi^(-17) = {M_E * PHI**(-17) * 1e9:.4e} eV (vs m_nu3 ~ {m3_est:.4e} eV)")

# These are way too small. Try different exponents
print(f"\n  Scanning for n that gives correct neutrino mass scale:")
for n in range(25, 45):
    m_nu = M_E * PHI**(-n) * 1e9  # eV
    if 1e-4 < m_nu < 1.0:
        print(f"  n={n}: m = {m_nu:.4e} eV")

# ============================================================
# APPROACH C: Phi-Tower with Negative Exponents
# ============================================================
print("\n" + "=" * 70)
print("APPROACH C: Direct Phi-Tower for Neutrinos")
print("=" * 70)

print(f"""
  If neutrinos follow the same formula m = m_e * phi^(5a + 6b),
  what (a,b) pairs give sub-eV masses?

  Need: m_e * phi^n < 1 eV => phi^n < 1/(m_e in eV) = 1/511000
  => n < ln(1/511000)/ln(phi) = {np.log(1/511000)/np.log(PHI):.1f}

  So n must be negative (large negative).
""")

# Find (a,b) pairs for neutrino-scale masses
target_masses_eV = [0.001, 0.0087, 0.05]  # rough estimates for nu1, nu2, nu3
target_names = ['nu_1', 'nu_2', 'nu_3']

print(f"  Searching for (a,b) with correct neutrino masses:")
for i, (target, name) in enumerate(zip(target_masses_eV, target_names)):
    target_GeV = target * 1e-9
    n_needed = np.log(target_GeV/M_E) / np.log(PHI)
    print(f"\n  {name}: target = {target:.4e} eV -> n = {n_needed:.2f}")

    # Find (a,b) with 5a+6b close to n_needed
    best = None
    best_err = 1e10
    for a in range(-10, 10):
        b_float = (n_needed - 5*a) / 6
        b = round(b_float)
        n = 5*a + 6*b
        m_pred = M_E * PHI**n * 1e9  # eV
        err = abs(n - n_needed)
        if err < best_err:
            best_err = err
            best = (a, b, n, m_pred)

    a, b, n, m_pred = best
    norm = abs(a**2 + a*b - b**2)
    print(f"    Best: (a,b) = ({a},{b}), n = {n}, m_pred = {m_pred:.4e} eV, |N| = {norm}")

# ============================================================
# APPROACH D: Mass-Squared Differences from Phi
# ============================================================
print("\n" + "=" * 70)
print("APPROACH D: Mass-Squared Differences")
print("=" * 70)

print(f"""
  The observable quantities are mass-SQUARED differences:
  R = dm2_32 / dm2_21 = {dm2_32/dm2_21:.4f}

  Test: is R related to phi?
""")

R = dm2_32 / dm2_21
print(f"  R = {R:.4f}")
print(f"\n  Phi-based candidates for R:")

candidates = [
    ("phi^5", PHI**5),
    ("phi^6", PHI**6),
    ("phi^7", PHI**7),
    ("5*phi^2", 5*PHI**2),
    ("10*phi", 10*PHI),
    ("phi^4 * phi", PHI**4 * PHI),
    ("12*phi^(1/2)", 12*PHI**0.5),
    ("20", 20),
    ("5^2", 25),
    ("30", 30),
    ("phi^7 / phi", PHI**7 / PHI),
    ("4*phi^3", 4*PHI**3),
    ("8*phi^2", 8*PHI**2),
    ("26", 26),  # phi-sector dimension
    ("2*phi^4", 2*PHI**4),
    ("6*phi^(3/2)", 6*PHI**1.5),
]

results = []
for label, val in candidates:
    err = abs(val - R)/R * 100
    results.append((err, val, label))
results.sort()

print(f"  {'err%':>6} {'value':>10} {'formula':<20}")
print("  " + "-" * 40)
for err, val, label in results[:10]:
    marker = " ***" if err < 5 else ""
    print(f"  {err:>5.1f}% {val:>10.4f} {label:<20}{marker}")

# Best match
print(f"\n  Best: {results[0][2]} = {results[0][1]:.4f} (err = {results[0][0]:.1f}%)")

# ============================================================
# APPROACH E: Neutrinos as "Conjugate Generations"
# ============================================================
print("\n" + "=" * 70)
print("APPROACH E: Neutrinos from Seesaw with Phi-Tower")
print("=" * 70)

print(f"""
  HYBRID APPROACH: Combine seesaw with phi-tower structure.

  If the heavy scale M_R also follows the phi-tower:
    M_R = m_e * phi^N_R (for some level N_R)

  Then: m_nu ~ m_D^2 / M_R = m_e * phi^(2*n_D) / (m_e * phi^N_R)
        m_nu = m_e * phi^(2*n_D - N_R)

  For the seesaw exponent to be correct:
    2*n_D - N_R = n_nu (a very negative number)

  Example: nu_tau with m_D = tau mass (n_D = 17):
    m_nu3 ~ 0.05 eV = m_e * phi^n_nu3
    n_nu3 = ln(0.05e-9/m_e) / ln(phi) = {np.log(0.05e-9/M_E)/np.log(PHI):.1f}
    N_R = 2*17 - ({np.log(0.05e-9/M_E)/np.log(PHI):.1f}) = {2*17 - np.log(0.05e-9/M_E)/np.log(PHI):.1f}
""")

# Find N_R for each neutrino
print(f"  Seesaw N_R determination:")
for name, n_D, m_nu_eV in [('nu_e', 0, 0.001), ('nu_mu', 11, 0.0087), ('nu_tau', 17, 0.05)]:
    m_nu_GeV = m_nu_eV * 1e-9
    n_nu = np.log(m_nu_GeV/M_E) / np.log(PHI)
    N_R = 2*n_D - n_nu
    M_R = M_E * PHI**N_R
    print(f"  {name}: n_D={n_D:>2}, n_nu={n_nu:>7.1f}, N_R={N_R:>6.1f}, M_R={M_R:.2e} GeV")

# Check if N_R values are consistent
N_R_values = []
for name, n_D, m_nu_eV in [('nu_mu', 11, 0.0087), ('nu_tau', 17, 0.05)]:
    m_nu_GeV = m_nu_eV * 1e-9
    n_nu = np.log(m_nu_GeV/M_E) / np.log(PHI)
    N_R = 2*n_D - n_nu
    N_R_values.append(N_R)

print(f"\n  N_R values: {[f'{x:.1f}' for x in N_R_values]}")
print(f"  Difference: {abs(N_R_values[1]-N_R_values[0]):.1f}")
print(f"  For a UNIVERSAL seesaw scale, N_R should be the same for all.")
print(f"  Difference of {abs(N_R_values[1]-N_R_values[0]):.1f} levels suggests")
print(f"  N_R is NOT universal -> seesaw with single scale FAILS.")

# What if N_R depends on generation?
print(f"\n  Generation-dependent seesaw:")
print(f"  If M_R_gen = m_e * phi^(N_R_base + delta_gen):")
for name, n_D, m_nu_eV in [('nu_mu', 11, 0.0087), ('nu_tau', 17, 0.05)]:
    m_nu_GeV = m_nu_eV * 1e-9
    n_nu = np.log(m_nu_GeV/M_E) / np.log(PHI)
    N_R = 2*n_D - n_nu
    print(f"    {name}: N_R = {N_R:.1f}")

# ============================================================
# APPROACH F: Neutrino Mass Ratio from Phi
# ============================================================
print("\n" + "=" * 70)
print("APPROACH F: Neutrino Mass Ratios")
print("=" * 70)

# The most model-independent test: mass RATIO between nu2 and nu3
r_nu = m3_est / m2_est
print(f"  m3/m2 = {m3_est:.4f}/{m2_est:.4f} = {r_nu:.4f}")
print(f"")
print(f"  Phi-based candidates for m3/m2:")

candidates_r = [
    ("phi", PHI),
    ("phi^2", PHI**2),
    ("phi^3", PHI**3),
    ("phi^4", PHI**4),
    ("2phi", 2*PHI),
    ("phi + 1 = phi^2", PHI+1),
    ("sqrt(5)", np.sqrt(5)),
    ("sqrt(5)*phi", np.sqrt(5)*PHI),
    ("phi^(5/2)", PHI**2.5),
    ("2*phi^2", 2*PHI**2),
    ("3*phi", 3*PHI),
    ("5", 5.0),
    ("6", 6.0),
    ("phi^(7/2)", PHI**3.5),
    ("sqrt(dm2_32/dm2_21)", np.sqrt(dm2_32/dm2_21)),
]

res_r = []
for label, val in candidates_r:
    err = abs(val - r_nu)/r_nu * 100
    res_r.append((err, val, label))
res_r.sort()

print(f"  {'err%':>6} {'value':>10} {'formula':<30}")
print("  " + "-" * 50)
for err, val, label in res_r[:10]:
    marker = " ***" if err < 5 else ""
    print(f"  {err:>5.1f}% {val:>10.4f} {label:<30}{marker}")

# Also test: is dm2_32/dm2_21 a nice phi expression?
ratio_dm2 = dm2_32/dm2_21
print(f"\n  dm2_32/dm2_21 = {ratio_dm2:.4f}")
print(f"  sqrt(dm2_32/dm2_21) = {np.sqrt(ratio_dm2):.4f} = m3/m2 (if m1<<m2)")

# ============================================================
# APPROACH G: PMNS Angles + Neutrino Masses
# ============================================================
print("\n" + "=" * 70)
print("APPROACH G: Connecting PMNS Angles to Mass Ratios")
print("=" * 70)

print(f"""
  PMNS angles: theta_12=33.41, theta_23=49.0, theta_13=8.54 deg
  These are well-predicted by our phi-tower.

  Can the MASS RATIO be derived from the MIXING ANGLES?

  In many neutrino models, there's a connection:
  tan(theta_12) ~ sqrt(m2/m3) or similar.

  Test:
""")

theta_12 = 33.41  # deg
theta_23 = 49.0
theta_13 = 8.54

print(f"  tan(theta_12) = {np.tan(np.radians(theta_12)):.4f}")
print(f"  m2/m3 = {m2_est/m3_est:.4f}")
print(f"  sqrt(m2/m3) = {np.sqrt(m2_est/m3_est):.4f}")
print(f"  tan(theta_12) vs sqrt(m2/m3): {abs(np.tan(np.radians(theta_12)) - np.sqrt(m2_est/m3_est))/np.sqrt(m2_est/m3_est)*100:.1f}% off")

# Phi-tower prediction for tan(theta_12)
# theta_12 = arctan(phi^(-1)) -> tan(theta_12) = phi^(-1) = 1/phi
tan_12_pred = 1/PHI
print(f"\n  Phi-tower: tan(theta_12) = 1/phi = {tan_12_pred:.4f}")
print(f"  If tan(theta_12)^2 = m2/m3:")
print(f"    m2/m3 = 1/phi^2 = {1/PHI**2:.4f}")
print(f"    Actual m2/m3 = {m2_est/m3_est:.4f}")
print(f"    Error: {abs(1/PHI**2 - m2_est/m3_est)/(m2_est/m3_est)*100:.1f}%")

print(f"\n  If m3/m2 = phi^2:")
m3_pred = m2_est * PHI**2
print(f"    m3 = m2 * phi^2 = {m2_est:.4f} * {PHI**2:.4f} = {m3_pred:.4f} eV")
print(f"    Actual m3 = {m3_est:.4f} eV")
print(f"    Error: {abs(m3_pred - m3_est)/m3_est*100:.1f}%")

print(f"\n  If m3/m2 = phi^3:")
m3_pred2 = m2_est * PHI**3
print(f"    m3 = m2 * phi^3 = {m2_est:.4f} * {PHI**3:.4f} = {m3_pred2:.4f} eV")
print(f"    Actual m3 = {m3_est:.4f} eV")
print(f"    Error: {abs(m3_pred2 - m3_est)/m3_est*100:.1f}%")

# ============================================================
# APPROACH H: DSI in Neutrino Sector
# ============================================================
print("\n" + "=" * 70)
print("APPROACH H: DSI (Discrete Scale Invariance)")
print("=" * 70)

print(f"""
  The DSI mechanism gives mass ratios as powers of phi.
  For charged leptons: m_tau/m_mu ~ phi^6, m_mu/m_e ~ phi^11/phi^0 = phi^11.

  For neutrinos, if the SAME DSI operates but in the conjugate sector:
  m_nu3/m_nu2 ~ phi^k for some k.

  Experimental: m3/m2 ~ {r_nu:.2f}
""")

# Find best phi^k
best_k = None
best_err_k = 1e10
for k_num in range(1, 40):
    for k_den in range(1, 5):
        k = k_num/k_den
        val = PHI**k
        err = abs(val - r_nu)/r_nu * 100
        if err < best_err_k:
            best_err_k = err
            best_k = k

print(f"  Best phi^k for m3/m2: k = {best_k:.2f}, phi^k = {PHI**best_k:.4f}, err = {best_err_k:.1f}%")

# Also check mass-squared ratio
r_dm2 = dm2_32/dm2_21
for k_num in range(1, 40):
    for k_den in range(1, 5):
        k = k_num/k_den
        val = PHI**k
        err = abs(val - r_dm2)/r_dm2 * 100
        if err < 2:
            print(f"  dm2_32/dm2_21 ~ phi^{k:.2f} = {val:.4f} (err = {err:.1f}%)")

# ============================================================
# CONCLUSIONS
# ============================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)

print(f"""
NEUTRINO MASS RESULTS:

1. DIRECT PHI-TOWER: Needs n ~ -27 to -30. No natural (a,b) pairs.
   Status: NEGATIVE (as expected from exp116)

2. SEESAW: Works if M_R ~ m_e * phi^{N_R_values[0]:.0f} to phi^{N_R_values[1]:.0f}.
   But N_R is NOT the same for all generations (differs by ~{abs(N_R_values[1]-N_R_values[0]):.0f}).
   Single-scale seesaw FAILS. Generation-dependent seesaw possible but ugly.
   Status: PARTIAL (needs generation-dependent heavy scale)

3. CONJUGATE SECTOR: phi^(-n) with n=11,17 gives masses too small.
   Status: NEGATIVE

4. MASS RATIOS:
   m3/m2 ~ {r_nu:.2f}
   Best phi match: phi^{best_k:.1f} at {best_err_k:.1f}%
   Status: WEAK (many candidates, not distinctive)

5. PMNS-MASS CONNECTION:
   m3/m2 = phi^2 gives {abs(m3_pred - m3_est)/m3_est*100:.1f}% error.
   tan(theta_12)^2 = 1/phi^2 predicts m2/m3, but {abs(1/PHI**2 - m2_est/m3_est)/(m2_est/m3_est)*100:.1f}% off.
   Status: SUGGESTIVE but not precise

6. DSI: Mass-squared ratio dm2_32/dm2_21 might match phi^k.
   Status: SEE RESULTS ABOVE

OVERALL: Neutrino masses remain the HARDEST open problem.
The theory's phi-tower structure doesn't naturally accommodate sub-eV masses.
A seesaw mechanism with phi-tower heavy scale is the best option,
but requires generation-dependent M_R (not elegant).

HONEST ASSESSMENT: Neutrinos likely require ADDITIONAL PHYSICS
beyond the 600-cell mass formula. Possible directions:
- Majorana mechanism with a NEW scale
- Radiative mass generation (loop-suppressed)
- Connection to cosmological constant (both unexplained)
""")

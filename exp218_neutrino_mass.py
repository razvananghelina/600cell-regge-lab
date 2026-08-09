# EXP-218: Neutrino Mass Mechanism from 600-Cell Framework
# ========================================================
# STATUS: EXPLORATORY
#
# The framework predicts PMNS angles (all 3 within 1-sigma) but NOT
# neutrino masses. Previous attempts (exp134): NEGATIVE.
#
# This experiment explores multiple approaches:
# 1. Galois seesaw revisited (shadow sector suppression)
# 2. Real vs complex irrep distinction (psi_5 vs psi_3/psi_4)
# 3. Radiative/loop suppression on the graph
# 4. DSI negative levels
# 5. Spectral combinations
#
# Key experimental data (PDG 2024, normal hierarchy):
#   Dm2_21 = 7.53e-5 eV^2
#   Dm2_31 = 2.455e-3 eV^2
#   m3 ~ 0.0496 eV, m2 ~ 0.00868 eV (if m1~0)
#   m2/m3 ~ 0.175
#
# RULE ZERO: no invention. Report what works and what fails.

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
PSI = 1 - PHI  # = -1/phi
ALPHA = 1/137.035999084
ALPHA_S = 0.1179
SIN2TW = 6.0/26.0

# Masses in eV
m_e = 0.51099895e6
m_mu = 105.6583755e6
m_tau = 1776.86e6
m_u = 2.16e6
m_d = 4.67e6
m_s = 93.4e6
m_c = 1.27e9
m_b = 4.18e9
m_t = 172.69e9

# Neutrino oscillation data (PDG 2024)
Dm2_21 = 7.53e-5  # eV^2
Dm2_31 = 2.455e-3  # eV^2  (normal hierarchy)
Dm2_ratio_exp = Dm2_21 / Dm2_31  # = 0.03067

# Derived neutrino masses (normal hierarchy, m1~0)
m3 = np.sqrt(Dm2_31)   # ~ 0.04955 eV
m2 = np.sqrt(Dm2_21)   # ~ 0.008678 eV
m1_approx = 0.0  # assumption

# 600-cell constants
N_eig = 9       # distinct eigenvalues
diam = 5        # graph diameter
a1 = 5          # first intersection number
b1 = 6          # second intersection number
N_vert = 120
N_edge = 720
deg = 12

print("="*70)
print("EXP-218: NEUTRINO MASS MECHANISM")
print("="*70)

print("\n--- Experimental data ---")
print(f"Dm2_21 = {Dm2_21:.4e} eV^2")
print(f"Dm2_31 = {Dm2_31:.4e} eV^2")
print(f"Dm2 ratio = {Dm2_ratio_exp:.5f}")
print(f"4*alpha  = {4*ALPHA:.5f} (our prediction, error {abs(4*ALPHA - Dm2_ratio_exp)/Dm2_ratio_exp*100:.2f}%)")
print(f"m3 (NH, m1~0) = {m3:.5f} eV")
print(f"m2 (NH, m1~0) = {m2:.6f} eV")
print(f"m2/m3 = {m2/m3:.5f}")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 1: POWER-OF-PHI SCAN FOR NEUTRINO MASSES")
print("="*70)
# What power of phi gives the ratio m_nu/m_e?
# m_f = m_e * phi^n => n = log(m_f/m_e) / log(phi)

for label, m_nu in [("m3", m3), ("m2", m2)]:
    ratio = m_nu / m_e
    n_phi = np.log(ratio) / np.log(PHI)
    print(f"{label}/m_e = {ratio:.4e}, n_phi = {n_phi:.3f}")
    # Check if n = 5a + 6b for small integers
    print(f"  Closest integer n = {round(n_phi)}")
    # Try n = 5a + 6b decomposition
    n_round = round(n_phi)
    for a in range(-10, 5):
        for b in range(-10, 5):
            if 5*a + 6*b == n_round:
                z = a + b*PHI
                zp = a + b*PSI
                N_norm = abs(a*a + a*b - b*b)
                print(f"  n={n_round}: (a,b)=({a},{b}), z={z:.3f}, z'={zp:.3f}, |N|={N_norm}")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 2: SEESAW-LIKE FORMULAS")
print("="*70)
# Standard seesaw: m_nu ~ m_D^2 / M_R
# What M_R gives correct masses?

print("\n--- Type-I seesaw with m_D = m_lepton ---")
for label, m_D, m_nu in [("e->nu1", m_e, m2), ("mu->nu2", m_mu, m2),
                           ("tau->nu3", m_tau, m3)]:
    M_R = m_D**2 / m_nu if m_nu > 0 else float('inf')
    print(f"  {label}: M_R = m_D^2/m_nu = {M_R:.3e} eV = {M_R/1e9:.3e} GeV")
    # Check if M_R is a nice power of phi * m_e or m_P
    m_P = 1.22e28  # Planck mass in eV
    if M_R > 0:
        n_MR = np.log(M_R/m_e) / np.log(PHI)
        print(f"    M_R/m_e = phi^{n_MR:.2f}")
        n_MP = np.log(M_R/m_P) / np.log(PHI)
        print(f"    M_R/m_P = phi^{n_MP:.2f}")

# Geometric seesaw: m_nu = m_lepton^2 / (m_e * phi^N)
print("\n--- Geometric seesaw: m_nu = m_lepton^2 / (m_e * phi^N) ---")
for label, m_lep, m_nu_exp in [("tau->nu3", m_tau, m3), ("mu->nu2", m_mu, m2)]:
    # phi^N = m_lep^2 / (m_e * m_nu)
    phiN = m_lep**2 / (m_e * m_nu_exp)
    N = np.log(phiN) / np.log(PHI)
    print(f"  {label}: phi^N = {phiN:.4e}, N = {N:.3f}")
    N_round = round(N)
    m_nu_pred = m_lep**2 / (m_e * PHI**N_round)
    err = abs(m_nu_pred - m_nu_exp)/m_nu_exp * 100
    print(f"    N={N_round}: m_nu_pred = {m_nu_pred:.5f} eV (exp: {m_nu_exp:.5f}, error: {err:.1f}%)")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 3: ALPHA-POWER SUPPRESSION")
print("="*70)
# m_nu = m_e * alpha^p * phi^q

print("\n--- m_nu = m_e * alpha^p * phi^q scan ---")
best_results = []
for p in range(1, 6):
    for q in range(-10, 10):
        m_pred = m_e * ALPHA**p * PHI**q
        for label, m_exp in [("m3", m3), ("m2", m2)]:
            if m_exp > 0:
                err = abs(m_pred - m_exp)/m_exp * 100
                if err < 15:
                    best_results.append((err, label, p, q, m_pred, m_exp))

best_results.sort()
print("Top matches (error < 15%):")
for err, label, p, q, m_pred, m_exp in best_results[:15]:
    n_eff = np.log(m_pred/m_e)/np.log(PHI)
    print(f"  {label} = m_e * alpha^{p} * phi^{q:+d} = {m_pred:.5f} eV "
          f"(exp: {m_exp:.5f}, error: {err:.1f}%, n_eff={n_eff:.2f})")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 4: REAL vs COMPLEX IRREP MECHANISM")
print("="*70)
# In 2T branching: psi_5 = REAL 2D, psi_3/psi_4 = COMPLEX pair
# Charged leptons: Dirac mass from complex irreps
# Neutrinos: Majorana mass from real irrep
# Suppression factor from representation theory

print("\n--- Representation-theoretic suppression ---")
# The Schur projector eigenvalues for the unitary rep:
# H eigenvalues were {8/15, 8/15, 16/15, 16/15, 80/15, 80/15}
# Ratios 1:2:10
H_eigs = np.array([8/15, 8/15, 16/15, 16/15, 80/15, 80/15])
print(f"H eigenvalues: {H_eigs}")
print(f"H eigenvalue ratios: 1 : {16/8:.0f} : {80/8:.0f}")
print(f"Product of H eigs: {np.prod(H_eigs):.6f}")
print(f"  = (8/15)^2 * (16/15)^2 * (80/15)^2 = {(8*16*80)**2/15**6:.6f}")

# The three Schur projectors have different "weights"
# P_3: projects onto psi_3 (complex), H-weight ~ related to 80/15
# P_4: projects onto psi_4 (complex conjugate), same weight
# P_5: projects onto psi_5 (real), H-weight ~ related to 8/15
#
# Ratio of real to complex weight: 8/80 = 1/10
print(f"\nWeight ratio (real/complex): 8/80 = {8/80}")
print(f"  = 1/{80//8} = 1/10")
print(f"  = 1/(2*a1) = 1/(2*{a1})")

# If neutrino mass ~ (real weight / complex weight) * charged lepton mass:
suppression = 1.0/10.0
for label, m_lep in [("e", m_e), ("mu", m_mu), ("tau", m_tau)]:
    m_nu_pred = m_lep * suppression
    print(f"  m_nu({label}) = m_{label}/10 = {m_nu_pred:.3e} eV")
# This gives way too large masses. Need more suppression.

# What about (real/complex)^k for some power k?
print("\n--- Power of suppression factor ---")
for k in range(1, 8):
    sup = (1.0/10.0)**k
    m_nu_tau = m_tau * sup
    print(f"  k={k}: m_tau * (1/10)^{k} = {m_nu_tau:.4e} eV", end="")
    if m_nu_tau > 0:
        err3 = abs(m_nu_tau - m3)/m3 * 100
        print(f"  (vs m3={m3:.4e}, error: {err3:.0f}%)", end="")
    print()

# ===================================================================
print("\n" + "="*70)
print("APPROACH 5: DSI POTENTIAL IN GALOIS SHADOW")
print("="*70)
# V(z) = sin^2(pi * ln|z| / ln(phi))
# For charged leptons: V(z) = 0, z = phi^k (unit)
# For neutrinos: what if mass comes from V(z') != 0?

print("\n--- Galois shadow potential for lepton units ---")
for gen, b_val, k_val in [(1, 0, 0), (2, 1, 2), (3, 2, 3)]:
    z = 1 + b_val * PHI  # on a=1 line
    zp = 1 + b_val * PSI  # Galois conjugate
    V_z = np.sin(np.pi * np.log(abs(z)) / np.log(PHI))**2
    V_zp = np.sin(np.pi * np.log(abs(zp)) / np.log(PHI))**2
    print(f"  Gen {gen}: b={b_val}, z=phi^{k_val}={z:.4f}, z'={zp:.4f}")
    print(f"    V(z) = {V_z:.6e}, V(z') = {V_zp:.6e}")
    print(f"    |z'|/|z| = {abs(zp)/abs(z):.6f}")
    # For units, V(z)=0 exactly. V(z') may not be 0.
    # What if neutrino mass ~ V(z') * m_lepton?
    if V_zp > 0:
        m_nu_pred = V_zp * m_e * PHI**(5*1 + 6*b_val)  # mass level of lepton
        print(f"    V(z') * m_lepton = {m_nu_pred:.4e} eV")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 6: COMBINATORIAL FORMULAS")
print("="*70)
# Try simple combinations of framework constants

print("\n--- Testing m3 = m_e * f(framework constants) ---")
candidates = []

# m3 ~ m_e * alpha^3 / phi^3
val = m_e * ALPHA**3 / PHI**3
err = abs(val - m3)/m3 * 100
candidates.append(("m_e * alpha^3 / phi^3", val, err))

# m3 ~ m_e * alpha^3 * phi^(-3)  (same as above)
# m3 ~ m_e * (alpha/phi)^3
val = m_e * (ALPHA/PHI)**3
candidates.append(("m_e * (alpha/phi)^3", val, err))

# m3 ~ m_e / (phi^(5*6))  = m_e / phi^30
val = m_e / PHI**30
err = abs(val - m3)/m3 * 100
candidates.append(("m_e / phi^30", val, err))

# m3 ~ m_e / (phi^(a1*b1))
val = m_e / PHI**(a1*b1)  # phi^30
candidates.append(("m_e / phi^(a1*b1)", val, err))

# m3 ~ m_e * alpha^(a1/2)
val = m_e * ALPHA**(a1/2.0)
err = abs(val - m3)/m3 * 100
candidates.append(("m_e * alpha^(5/2)", val, err))

# m3 ~ m_tau * alpha^2
val = m_tau * ALPHA**2
err = abs(val - m3)/m3 * 100
candidates.append(("m_tau * alpha^2", val, err))

# m3 ~ m_tau * alpha^2 / phi
val = m_tau * ALPHA**2 / PHI
err = abs(val - m3)/m3 * 100
candidates.append(("m_tau * alpha^2 / phi", val, err))

# m3 ~ m_tau * alpha^2 * phi
val = m_tau * ALPHA**2 * PHI
err = abs(val - m3)/m3 * 100
candidates.append(("m_tau * alpha^2 * phi", val, err))

# m3 ~ m_mu^2 / (m_tau * phi^n) -- inter-generation seesaw
for n in range(-3, 10):
    val = m_mu**2 / (m_tau * PHI**n)
    err = abs(val - m3)/m3 * 100
    if err < 20:
        candidates.append((f"m_mu^2 / (m_tau * phi^{n})", val, err))

# m3 ~ m_e^2 / (m_mu * phi^n)
for n in range(-5, 15):
    val = m_e**2 / (m_mu * PHI**n)
    err = abs(val - m3)/m3 * 100
    if err < 20:
        candidates.append((f"m_e^2 / (m_mu * phi^{n})", val, err))

# m3 ~ sqrt(m_e * m_mu * alpha^n)
for n in range(1, 6):
    val = np.sqrt(m_e * m_mu * ALPHA**n)
    err = abs(val - m3)/m3 * 100
    if err < 30:
        candidates.append((f"sqrt(m_e * m_mu * alpha^{n})", val, err))

# m3 ~ m_e * alpha * phi^(-n) -- alpha suppression
for n in range(0, 20):
    val = m_e * ALPHA * PHI**(-n)
    err = abs(val - m3)/m3 * 100
    if err < 10:
        candidates.append((f"m_e * alpha / phi^{n}", val, err))

# m3 ~ m_e * alpha^2 * phi^n -- alpha^2 suppression
for n in range(-10, 20):
    val = m_e * ALPHA**2 * PHI**n
    err = abs(val - m3)/m3 * 100
    if err < 10:
        candidates.append((f"m_e * alpha^2 * phi^{n}", val, err))

# m3 ~ m_e * sin^2(tW) * alpha * phi^n
for n in range(-10, 10):
    val = m_e * SIN2TW * ALPHA * PHI**n
    err = abs(val - m3)/m3 * 100
    if err < 10:
        candidates.append((f"m_e * sin2tW * alpha * phi^{n}", val, err))

# m3 ~ m_e / (phi^n * N_eig * diam)
for n in range(20, 35):
    val = m_e / (PHI**n)
    err = abs(val - m3)/m3 * 100
    if err < 10:
        candidates.append((f"m_e / phi^{n}", val, err))

candidates.sort(key=lambda x: x[2])
print("\nBest formulas for m3 (sorted by error):")
for formula, val, err in candidates[:20]:
    n_eff = np.log(val/m_e)/np.log(PHI) if val > 0 else 0
    print(f"  {formula:40s} = {val:.5f} eV  (error: {err:.2f}%, n_eff={n_eff:.2f})")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 7: MASS RATIO PATTERNS")
print("="*70)
# Instead of absolute masses, look at ratios

print("\n--- Neutrino mass ratios ---")
r23 = m2/m3
print(f"m2/m3 = {r23:.5f}")

# Check various expressions
tests = [
    ("sqrt(4*alpha)", np.sqrt(4*ALPHA)),
    ("sqrt(Dm2_ratio)", np.sqrt(Dm2_ratio_exp)),
    ("2*sqrt(alpha)", 2*np.sqrt(ALPHA)),
    ("phi^(-2) * sqrt(2)", PHI**(-2) * np.sqrt(2)),
    ("1/(phi^2 * sqrt(5))", 1/(PHI**2 * np.sqrt(5))),
    ("sin(theta_C)", np.sin(np.arctan(PHI**(-3)))),
    ("sin^2(theta_C)", np.sin(np.arctan(PHI**(-3)))**2),
    ("alpha_s / phi", ALPHA_S / PHI),
    ("1/phi^(a1-1)", 1/PHI**(a1-1)),
    ("1/phi^4", 1/PHI**4),
    ("alpha_s * phi^(-1)", ALPHA_S * PHI**(-1)),
    ("1/(2*phi^3)", 1/(2*PHI**3)),
    ("alpha_s / (2*phi)", ALPHA_S / (2*PHI)),
    ("sin^2(tW) / phi", SIN2TW / PHI),
    ("alpha^(1/2) * phi^(-1)", ALPHA**(0.5) / PHI),
    ("alpha^(1/3)", ALPHA**(1.0/3)),
    ("phi^(-3.64)", PHI**(-3.64)),
]

print(f"\nm2/m3 = {r23:.5f}")
for label, val in tests:
    err = abs(val - r23)/r23 * 100
    print(f"  {label:30s} = {val:.5f}  (error: {err:.1f}%)")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 8: GALOIS SEESAW WITH NORM")
print("="*70)
# For charged leptons: z = phi^k (unit), |N(z)| = 1
# For neutrinos: maybe mass ~ |N(z)| / something?
# Or: m_nu / m_lepton ~ |z'| / |z| evaluated at some specific point

print("\n--- z and z' for lepton (a,b) assignments ---")
leptons = [
    ("e",   0, 0, m_e),
    ("mu",  1, 1, m_mu),
    ("tau", 1, 2, m_tau),
]

for name, a, b, m_lep in leptons:
    z = a + b*PHI
    zp = a + b*PSI
    N_z = z * zp  # = a^2 + ab - b^2
    print(f"  {name}: (a,b)=({a},{b}), z={z:.4f}, z'={zp:.4f}, N(z)={N_z:.4f}")
    if abs(z) > 0:
        # Shadow suppression
        ratio = abs(zp) / abs(z)
        m_nu_shadow = m_lep * ratio
        print(f"    |z'|/|z| = {ratio:.6f}")
        print(f"    m_lep * |z'|/|z| = {m_nu_shadow:.4e} eV")

        # Deeper suppression: (|z'|/|z|)^k
        for k in range(2, 6):
            m_pred = m_lep * ratio**k
            print(f"    m_lep * (|z'|/|z|)^{k} = {m_pred:.4e} eV")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 9: NORM-BASED MASS FORMULA")
print("="*70)
# For charged leptons: m = m_e * phi^(5a+6b), z = unit (|N|=1)
# For neutrinos: what if m_nu = m_lepton * |N(z_nu)|^(-p) for some z_nu?
# Or: m_nu = m_e * phi^n_nu where n_nu is on a DIFFERENT line (not a=1)?

print("\n--- Exploring non-unit elements for neutrinos ---")
print("If neutrinos are on a=0 line (pure phi-sector):")
for b in range(1, 5):
    z = b * PHI
    zp = b * PSI
    N_z = z * zp  # = -b^2
    n = 6*b
    m_pred = m_e * PHI**n
    print(f"  b={b}: z={z:.3f}, z'={zp:.3f}, N(z)={N_z:.1f}, "
          f"n=6b={n}, m=m_e*phi^{n}={m_pred:.3e} eV")

print("\nIf neutrinos are on a=2 line:")
for b in range(-2, 5):
    z = 2 + b * PHI
    zp = 2 + b * PSI
    N_z = z * zp  # = 4 + 2b - b^2
    n = 10 + 6*b
    if n > -40 and n < 0:  # only negative n for small masses
        m_pred = m_e * PHI**n
        print(f"  b={b}: z={z:.3f}, z'={zp:.3f}, N(z)={N_z:.1f}, "
              f"n=10+6b={n}, m=m_e*phi^{n}={m_pred:.4e} eV")

print("\nIf neutrinos have NEGATIVE n on a=1:")
for b in range(-6, 0):
    z = 1 + b * PHI
    zp = 1 + b * PSI
    N_z = z * zp  # = 1 + b - b^2
    n = 5 + 6*b
    m_pred = m_e * PHI**n if n != 0 else m_e
    print(f"  b={b}: z={z:.4f}, z'={zp:.4f}, N(z)={N_z:.0f}, "
          f"n=5+6b={n}, m=m_e*phi^{n}={m_pred:.4e} eV")

# ===================================================================
print("\n" + "="*70)
print("APPROACH 10: DEMOCRATIC MIXING EIGENVALUE")
print("="*70)
# The democratic bare mixing matrix has eigenvalues {1, -1/4, -1/4}
# The ratio |-1/4| / 1 = 1/4
# What if neutrino mass involves this?

M_bare = np.array([[1/6, 5/12, 5/12],
                    [5/12, 1/6, 5/12],
                    [5/12, 5/12, 1/6]])
eigs = np.linalg.eigvalsh(M_bare)
print(f"Bare mixing eigenvalues: {sorted(eigs, reverse=True)}")
print(f"Ratio |lam2/lam1| = {abs(eigs[0]/eigs[-1]):.4f}")

# Democratic matrix structure:
# M = (1-p)*I_3/3 + p*J/3... no.
# M = pI + (1-p)/(n-1) * (J-I) where n=3 (generations)
# Hmm, actually M_ii = 1/6, M_ij = 5/12
# det(M) = ...
det_M = np.linalg.det(M_bare)
print(f"det(M) = {det_M:.6f}")
print(f"  = -1/64 = {-1/64:.6f}")
print(f"  Ratio: {det_M / (-1/64):.4f}")

# Trace
print(f"Tr(M) = {np.trace(M_bare):.6f} = 1/2 = {0.5:.6f}")

# ===================================================================
print("\n" + "="*70)
print("SUMMARY OF BEST CANDIDATES")
print("="*70)

print("""
Key findings:
- n_phi(m3) = log(m3/m_e)/log(phi) shows what DSI level neutrinos would be at
- Need to check if any formula has geometric meaning

Best formulas will be listed above from the scan.

OPEN QUESTIONS:
1. Why are neutrinos massive at all? (framework gives massless naturally)
2. Is the mass from Galois shadow (z' sector)?
3. Is the mass from real irrep (psi_5) vs complex (psi_3/psi_4)?
4. Does the democratic mixing eigenvalue -1/4 play a role?
""")

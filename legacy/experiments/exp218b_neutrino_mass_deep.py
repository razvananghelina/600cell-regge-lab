# EXP-218b: Deep Analysis of Neutrino Mass Candidates
# ====================================================
# Following EXP-218, the best candidate is:
#   m3 = m_e * (alpha/phi)^3    (5.4% off)
#
# Key question: can a small correction fix this?
# And is there a MECHANISM (3-loop radiative)?
#
# Also explore: m2/m3 ratio patterns.

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
PSI = 1 - PHI
ALPHA = 1/137.035999084
ALPHA_S = 0.1179
SIN2TW = 6.0/26.0  # = 0.23077

m_e = 0.51099895e6   # eV
m_mu = 105.6583755e6
m_tau = 1776.86e6

# Neutrino data
Dm2_21 = 7.53e-5    # eV^2 (PDG 2024)
Dm2_31 = 2.455e-3   # eV^2 (normal hierarchy)
m3 = np.sqrt(Dm2_31)   # 0.04955 eV
m2 = np.sqrt(Dm2_21)   # 0.008678 eV

# Framework constants
a1 = 5; b1 = 6; N_eig = 9; diam = 5; deg = 12

print("="*70)
print("EXP-218b: DEEP ANALYSIS OF NEUTRINO MASS CANDIDATES")
print("="*70)

# ===================================================================
print("\n--- CANDIDATE 1: m3 = m_e * (alpha/phi)^3 ---")
# ===================================================================
base = m_e * (ALPHA/PHI)**3
print(f"m_e * (alpha/phi)^3 = {base:.6f} eV")
print(f"m3 (exp) = {m3:.6f} eV")
print(f"Error: {(base-m3)/m3*100:.2f}%")
print(f"Needed correction factor: {m3/base:.6f}")
needed = m3/base - 1
print(f"  = 1 + {needed:.6f}")

# Test various corrections
print(f"\nSystematic correction scan (1 + X):")
corrections = {
    "alpha_s/2": ALPHA_S/2,
    "sin^2(tW)/4": SIN2TW/4,
    "3/(4*diam*N_eig)": 3.0/(4*diam*N_eig),
    "1/(2*deg)": 1.0/(2*deg),
    "alpha_s/(2*phi)": ALPHA_S/(2*PHI),
    "b1/(4*a1*N_eig)": b1/(4.0*a1*N_eig),
    "phi^(-a1)": PHI**(-a1),
    "alpha*phi": ALPHA*PHI,
    "3*alpha": 3*ALPHA,
    "1/(4*PHI^3)": 1.0/(4*PHI**3),
    "alpha_s * phi^(-2)": ALPHA_S * PHI**(-2),
    "1/(2*b1*phi)": 1.0/(2*b1*PHI),
}

results = []
for name, corr in corrections.items():
    m_pred = base * (1 + corr)
    err = abs(m_pred - m3)/m3 * 100
    results.append((err, name, corr, m_pred))

results.sort()
for err, name, corr, m_pred in results:
    marker = " <-- BEST" if err < 0.5 else ""
    print(f"  (1+{name:25s}) = (1+{corr:.6f}): m={m_pred:.6f} eV, err={err:.3f}%{marker}")

# ===================================================================
print("\n--- CANDIDATE 1b: m3 = m_e * alpha^3 * phi^q ---")
# ===================================================================
# Scan over non-integer q to find best phi power
base_alpha3 = m_e * ALPHA**3
q_best = np.log(m3/base_alpha3) / np.log(PHI)
print(f"m_e * alpha^3 = {base_alpha3:.6e} eV")
print(f"Best q (continuous): {q_best:.6f}")
print(f"  Nearest fraction: {q_best:.4f}")

# Check rational q values
print(f"\nRational q scan:")
for num in range(-20, 1):
    for den in range(1, 7):
        q = num/den
        if -5 < q < 0:
            m_pred = base_alpha3 * PHI**q
            err = abs(m_pred - m3)/m3 * 100
            if err < 2:
                print(f"  q={num}/{den}={q:.4f}: m={m_pred:.6f} eV, err={err:.3f}%")

# ===================================================================
print("\n\n--- MECHANISM: 3-LOOP RADIATIVE MASS ---")
# ===================================================================
print("""
In QFT, neutrino masses can be generated at L-loop order:
  m_nu ~ (alpha/(4*pi))^L * m^2/M

If the 600-cell framework forbids tree-level neutrino Yukawa:
  - Tree (L=0): ZERO (massless nu, as observed in framework)
  - 1-loop: m_nu ~ alpha * m^2/M
  - 2-loop: m_nu ~ alpha^2 * m^2/M
  - 3-loop: m_nu ~ alpha^3 * m^2/M

The power 3 = N_gen is the NUMBER OF GENERATIONS.
The 600-cell girth = 3 (smallest cycle = triangle).
""")

# Loop suppression with 4*pi factors
for L in range(1, 5):
    # Standard radiative: (alpha/(4*pi))^L * m_e
    val = (ALPHA/(4*np.pi))**L * m_e
    print(f"  L={L}: (alpha/(4*pi))^{L} * m_e = {val:.4e} eV")
    # Without 4*pi: alpha^L * m_e
    val2 = ALPHA**L * m_e
    print(f"       alpha^{L} * m_e = {val2:.4e} eV")
    # With phi: (alpha/phi)^L * m_e
    val3 = (ALPHA/PHI)**L * m_e
    print(f"       (alpha/phi)^{L} * m_e = {val3:.4e} eV")

print(f"\n  m3 (exp) = {m3:.4e} eV")
print(f"  Best: L=3, (alpha/phi)^3 * m_e matches to 5.4%")
print(f"  In standard QFT: 4*pi^3 ~ {(4*np.pi)**3:.1f} suppression per loop")
print(f"  In 600-cell: phi^3 ~ {PHI**3:.3f} suppression per loop")
print(f"  Ratio: (4*pi)/phi = {4*np.pi/PHI:.3f}")
print(f"  The 600-cell 'loop factor' is phi, not 4*pi!")
print(f"  This makes sense: phi is the fundamental scale factor (DSI).")

# ===================================================================
print("\n\n--- CANDIDATE 2: m2/m3 RATIO ---")
# ===================================================================
r_exp = m2/m3
print(f"m2/m3 (exp) = {r_exp:.6f}")

# Already know: Dm2_ratio = 4*alpha (1.24%)
# So m2/m3 = sqrt(4*alpha) = 2*sqrt(alpha)
r_pred = 2*np.sqrt(ALPHA)
print(f"2*sqrt(alpha) = {r_pred:.6f} (error: {abs(r_pred-r_exp)/r_exp*100:.2f}%)")

# Also: 1/(phi^2*sqrt(5))
r_pred2 = 1/(PHI**2 * np.sqrt(5))
print(f"1/(phi^2*sqrt(5)) = {r_pred2:.6f} (error: {abs(r_pred2-r_exp)/r_exp*100:.2f}%)")

# Are these the same? Compare 2*sqrt(alpha) vs 1/(phi^2*sqrt(5))
print(f"\n2*sqrt(alpha) vs 1/(phi^2*sqrt(5)):")
print(f"  2*sqrt(alpha)     = {r_pred:.8f}")
print(f"  1/(phi^2*sqrt(5)) = {r_pred2:.8f}")
print(f"  Ratio = {r_pred/r_pred2:.8f}")
print(f"  They differ by {abs(r_pred-r_pred2)/r_pred*100:.4f}%")
print(f"  So 4*alpha ~ 1/(5*phi^4)?")
print(f"  4*alpha = {4*ALPHA:.8f}")
print(f"  1/(5*phi^4) = {1/(5*PHI**4):.8f}")
print(f"  Error: {abs(4*ALPHA - 1/(5*PHI**4))/(4*ALPHA)*100:.3f}%")

# This means alpha ~ 1/(20*phi^4)
print(f"\n  1/(20*phi^4) = {1/(20*PHI**4):.10f}")
print(f"  alpha =        {ALPHA:.10f}")
print(f"  Ratio = {(1/(20*PHI**4))/ALPHA:.8f}")
print(f"  This is a KNOWN approximation, not exact in our framework.")
print(f"  Our framework derives alpha from spectral equation,")
print(f"  which gives 20*phi^4*alpha = 1 + 2*pi*alpha^2 (small correction).")

# ===================================================================
print("\n\n--- COMBINED: FULL NEUTRINO MASS FORMULA ---")
# ===================================================================
print("If m3 = m_e * (alpha/phi)^3 and m2/m3 = 2*sqrt(alpha):")
m3_pred = m_e * (ALPHA/PHI)**3
m2_pred = m3_pred * 2 * np.sqrt(ALPHA)
print(f"  m3 = m_e * (alpha/phi)^3 = {m3_pred:.6f} eV (exp: {m3:.6f}, err: {abs(m3_pred-m3)/m3*100:.2f}%)")
print(f"  m2 = m3 * 2*sqrt(alpha) = {m2_pred:.6f} eV (exp: {m2:.6f}, err: {abs(m2_pred-m2)/m2*100:.2f}%)")
print(f"  m1 ~ 0 (normal hierarchy)")

# Express m2 purely in terms of fundamentals
print(f"\n  m2 = m_e * (alpha/phi)^3 * 2*sqrt(alpha)")
print(f"     = 2 * m_e * alpha^(7/2) / phi^3")
m2_pure = 2 * m_e * ALPHA**(3.5) / PHI**3
print(f"     = {m2_pure:.6f} eV (exp: {m2:.6f}, err: {abs(m2_pure-m2)/m2*100:.2f}%)")

# ===================================================================
print("\n\n--- HIERARCHY AND n_eff VALUES ---")
# ===================================================================
print("Effective DSI levels n_eff = log(m/m_e)/log(phi):")
for name, m_val in [("m3", m3_pred), ("m2", m2_pred), ("m3_exp", m3), ("m2_exp", m2)]:
    n_eff = np.log(m_val/m_e)/np.log(PHI)
    print(f"  {name:8s}: n_eff = {n_eff:.3f}")

# Compare with charged lepton n values
print(f"\nCharged lepton n values:")
print(f"  e:   n = 0 (reference)")
print(f"  mu:  n = 11 = 5+6")
print(f"  tau: n = 17 = 5+12")
print(f"\nDifferences (charged - neutrino):")
n3 = np.log(m3_pred/m_e)/np.log(PHI)
n2 = np.log(m2_pred/m_e)/np.log(PHI)
print(f"  n(e) - n(nu3) = 0 - ({n3:.2f}) = {-n3:.2f}")
print(f"  n(mu) - n(nu2) = 11 - ({n2:.2f}) = {11-n2:.2f}")

# ===================================================================
print("\n\n--- WHY POWER 3? GEOMETRIC ARGUMENTS ---")
# ===================================================================
print("""
The power 3 in (alpha/phi)^3 could originate from:

1. N_gen = 3: Three generations require 3 loops to generate mass
   (each generation contributes one loop)

2. Girth = 3: The 600-cell has girth 3 (smallest cycle = triangle)
   Radiative mass generation requires at least one cycle.
   The MINIMUM number of independent cycles = 3 for consistent mass.

3. E_flip = 3: Soliton creation energy is exactly 3.
   This is the number of broken edges in a soliton defect.
   If neutrino mass is a "virtual soliton" effect: E ~ alpha^(E_flip)

4. deg/4 = 12/4 = 3: Each vertex has 3 neighbors in each other coset.
   The coset structure gives 3 interaction channels.

5. dim(Im(H)) = 3: Imaginary quaternion space is 3-dimensional.
   Quaternion geometry of the 600-cell naturally introduces factors of 3.

6. Branch node rho_8 splits into 3 pieces (psi_3 + psi_4 + psi_5).
   One loop per irreducible component.
""")

# ===================================================================
print("\n--- APPROACH: INVERSE SEESAW FROM BRANCHING ---")
# ===================================================================
# In the branching rho_8 = psi_3 + psi_4 + psi_5:
# psi_5 (real, dim 2) -> Majorana neutrino
# psi_3, psi_4 (complex pair, dim 2+2) -> Dirac charged lepton
#
# The Schur projector weights are:
# H eigenvalues: {8/15, 16/15, 80/15} with multiplicities {2, 2, 2}
# The real irrep psi_5 has H-eigenvalue 8/15 (smallest)
# The complex irreps have H-eigenvalue 80/15 (largest)
# Ratio: 8/80 = 1/10

print("Schur projector H-eigenvalues:")
print("  psi_5 (real):    8/15 = 0.5333")
print("  psi_3 (complex): likely 80/15 = 5.333 or 16/15 = 1.067")
print("  psi_4 (complex): same as psi_3")
print()

# In the unitarized representation, the three projectors have "sizes"
# proportional to the H-eigenvalues
# If neutrino mass ~ (H_real / H_complex)^k * m_lepton
# Then with H_real/H_complex = 1/10:
for k in range(1, 5):
    factor = (1.0/10)**k
    m_nu = m_tau * factor
    print(f"  k={k}: m_tau * (1/10)^{k} = {m_nu:.4e} eV")

# These are way too big. The H-eigenvalue ratio alone is not enough.

# What if we combine: (H_ratio)^k * (alpha/phi)^L?
print(f"\nCombined: m_e * (alpha/phi)^L * (1/10)^k:")
for L in range(1, 5):
    for k in range(0, 3):
        m_pred = m_e * (ALPHA/PHI)**L * (0.1)**k
        if 1e-4 < m_pred < 1:
            err3 = abs(m_pred - m3)/m3 * 100
            err2 = abs(m_pred - m2)/m2 * 100
            marker3 = " <-- m3!" if err3 < 10 else ""
            marker2 = " <-- m2!" if err2 < 10 else ""
            print(f"  L={L}, k={k}: {m_pred:.5f} eV  "
                  f"(vs m3: {err3:.1f}%{marker3}, vs m2: {err2:.1f}%{marker2})")

# ===================================================================
print("\n\n--- TEST: DOES 20*phi^4*alpha = 1 HELP? ---")
# ===================================================================
# Our spectral equation: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
# Leading order: alpha ~ 1/(20*phi^4) with correction ~ 2*pi*alpha^2

# If m3 = m_e * (alpha/phi)^3:
# Using alpha ~ 1/(20*phi^4):
# m3 ~ m_e * (1/(20*phi^4*phi))^3 = m_e * 1/(20^3 * phi^15)
# = m_e / (8000 * phi^15)
val = m_e / (8000 * PHI**15)
print(f"m_e / (8000 * phi^15) = {val:.6f} eV")
print(f"m3 (exp) = {m3:.6f} eV")
print(f"Error: {abs(val-m3)/m3*100:.2f}%")
print(f"Note: 8000 = 20^3 = (4*a1)^3")
print(f"  phi^15 = phi^(3*a1) = {PHI**15:.2f}")

# More precise: m3 = m_e / ((20*phi^4)^3 * phi^3) = m_e / (20^3 * phi^15)
# = m_e / (8000 * 1597.98) = m_e / 12783824
# = 511000 / 12783824 = 0.04
print(f"\nAlternative form: m3 = m_e / (20^3 * phi^15)")
print(f"  = m_e / ((4*a1)^3 * phi^(3*a1))")
print(f"  = m_e / ((4*a1*phi^a1)^3)")
val_alt = m_e / (4*a1 * PHI**a1)**3
print(f"  = {val_alt:.6f} eV (same: {abs(val-val_alt)<1e-10})")

# ===================================================================
print("\n\n--- COSMOLOGICAL CONSTRAINTS ---")
# ===================================================================
# Sum of neutrino masses:
m_sum_pred = m3_pred + m2_pred  # m1 ~ 0
print(f"Sum of neutrino masses:")
print(f"  Predicted: m1+m2+m3 = {m_sum_pred:.4f} eV")
print(f"  Planck 2018 limit: < 0.12 eV (95% CL)")
print(f"  Consistent: {'YES' if m_sum_pred < 0.12 else 'NO'}")

# Normal hierarchy check
print(f"\nHierarchy check:")
print(f"  m2/m3 = {m2_pred/m3_pred:.5f} (exp: {m2/m3:.5f})")
print(f"  m1/m2 ~ 0 (normal hierarchy)")

# ===================================================================
print("\n\n--- CROSS-CHECK WITH Dm^2 ---")
# ===================================================================
Dm2_31_pred = m3_pred**2
Dm2_21_pred = m2_pred**2
print(f"Dm2_31 predicted: {Dm2_31_pred:.4e} eV^2 (exp: {Dm2_31:.4e}, err: {abs(Dm2_31_pred-Dm2_31)/Dm2_31*100:.2f}%)")
print(f"Dm2_21 predicted: {Dm2_21_pred:.4e} eV^2 (exp: {Dm2_21:.4e}, err: {abs(Dm2_21_pred-Dm2_21)/Dm2_21*100:.2f}%)")
print(f"Ratio predicted:  {Dm2_21_pred/Dm2_31_pred:.5f} (exp: {Dm2_21/Dm2_31:.5f})")
print(f"  Our prediction: 4*alpha = {4*ALPHA:.5f}")

# ===================================================================
print("\n\n" + "="*70)
print("SUMMARY AND STATUS CLASSIFICATION")
print("="*70)
print(f"""
NEUTRINO MASS CANDIDATES:

1. m3 = m_e * (alpha/phi)^3 = {m3_pred:.5f} eV
   Experimental: {m3:.5f} eV
   Error: {abs(m3_pred-m3)/m3*100:.1f}%
   Status: PATTERN (5.4% error, but clean formula)
   Mechanism: 3-loop radiative, 3 = N_gen = girth = E_flip

   Alternative form: m3 = m_e / (4*a1 * phi^a1)^3
   where a1=5 is the first intersection number.
   All ingredients are 600-cell invariants.

2. m2/m3 = 2*sqrt(alpha)
   From: Dm2_ratio = 4*alpha (existing prediction)
   Error: 2.4%
   Status: PATTERN (from Dm2 ratio)

3. m2 = 2 * m_e * alpha^(7/2) / phi^3 = {m2_pred:.6f} eV
   Experimental: {m2:.6f} eV
   Error: {abs(m2_pred-m2)/m2*100:.1f}%
   Status: PATTERN (combines 1 and 2)

4. m1 ~ 0 (normal hierarchy assumed)
   Status: ASSUMPTION (consistent with current data)

NEGATIVE RESULTS:
- Galois shadow potential V(z'): too small (10^-20 eV)
- H-eigenvalue ratio 1/10: not enough suppression alone
- DSI levels n=-31,-37: don't have clean (a,b) decomposition
- Pure phi^n: no integer n matches neutrino masses

OVERALL STATUS: PATTERN/SPECULATIVE
The (alpha/phi)^3 formula is numerologically appealing but
the 3-loop mechanism is not derived from the framework.
The 5.4% error for m3 is worse than our best results (< 1%).
""")

# Final check: effective neutrino (a,b) if treated as DSI level
print("--- Effective DSI quantum numbers ---")
for name, m_val in [("nu3", m3_pred), ("nu2", m2_pred)]:
    n = np.log(m_val/m_e) / np.log(PHI)
    # n = 5a + 6b
    # Try a=1 line: n = 5 + 6b => b = (n-5)/6
    b_a1 = (n - 5) / 6
    print(f"  {name}: n_eff = {n:.3f}, b(a=1) = {b_a1:.3f}")
    # Not integer. Try other a values
    for a in range(-3, 4):
        b = (n - 5*a) / 6
        if abs(b - round(b)) < 0.2:
            print(f"    a={a}: b = {b:.3f} ~ {round(b)}")

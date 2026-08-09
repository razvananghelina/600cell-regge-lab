"""
EXP-514: Dark Matter Abundance from Orbifold Vacuum Energy
============================================================
IDEA: The DM/baryon ratio comes from the Galois structure of C^2/2I.

Physical sector: eigenvalues involve phi
Dark sector: eigenvalues involve phi' = -1/phi (Galois conjugate)

The orbifold vacuum energy E_vac = -539/43200 was computed for the
PHYSICAL sector. What is E_vac for the DARK sector?

The ratio E_vac(dark)/E_vac(physical) might give Omega_DM/Omega_b.
"""

import numpy as np
from fractions import Fraction
import sys
sys.path.insert(0, '.')
from commons import (PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen,
                      alpha_em, alpha_s, ln_inv_alpha, degree)

print("=" * 70)
print("EXP-514: DM ABUNDANCE FROM ORBIFOLD GALOIS STRUCTURE")
print("=" * 70)

# Observed
Omega_DM_over_b = 5.36  # Planck 2018: Omega_DM / Omega_b

# ============================================================
# STEP 1: Physical vs Dark eigenvalues
# ============================================================
print(f"\n{'='*70}")
print("STEP 1: PHYSICAL VS DARK EIGENVALUES")
print(f"{'='*70}")

# 600-cell Laplacian eigenvalues (physical sector):
# Galois pairs:
#   L_1 = 12 - 6*phi <-> L_8 = 6 + 6*phi = 12 - 6*phi'  (mult 4)
#   L_2 = 12 - 2*sqrt(5) <-> L_6 = 12 + 2*sqrt(5)       (mult 9)
# Fixed: 0(1), 9(16), 12(25), 14(36), 15(16)

# Under Galois sigma: phi -> phi', sqrt(5) -> -sqrt(5)
# Physical eigenvalues: L_1 = 12-6*phi, L_2 = 12-2*sqrt(5)
# Dark eigenvalues:     L_1' = 12-6*phi' = 6+6*phi = L_8
#                       L_2' = 12+2*sqrt(5) = L_6
# Fixed eigenvalues: UNCHANGED (rational)

# So sigma just SWAPS the Galois pairs!
# The SPECTRUM is the SAME (just reordered).
# Therefore: E_vac(physical) = E_vac(dark) EXACTLY.

print(f"""
  Under Galois sigma: phi -> phi', sqrt(5) -> -sqrt(5)

  Physical spectrum: {{0, 12-6phi, 12-2s5, 9, 12, 14, 12+2s5, 15, 6+6phi}}
  Dark spectrum:     {{0, 6+6phi, 12+2s5, 9, 12, 14, 12-2s5, 15, 12-6phi}}

  These are the SAME SET of eigenvalues, just reordered!
  The spectrum is Galois-INVARIANT (as a set).

  Therefore: E_vac(physical) = E_vac(dark) = -539/43200

  Ratio: E_vac(dark) / E_vac(physical) = 1

  This does NOT give the DM abundance.
""")

# ============================================================
# STEP 2: Weighted vacuum energy (quantum dimensions)
# ============================================================
print(f"{'='*70}")
print("STEP 2: QUANTUM-DIMENSION WEIGHTED VACUUM ENERGY")
print(f"{'='*70}")

# The orbifold E_vac used UNIFORM weights (each group element counts equally).
# But in TQFT, each irrep has a quantum dimension d_q.
# Under Galois: d_q(phi) -> d_q(phi')

# McKay correspondence: 2I irreps have dimensions 1,2,3,4,5,6,4,2,3
# The quantum dimensions at level k=3 (k+2=a1=5):
# d_q(j) = sin((2j+1)*pi/a1) / sin(pi/a1)

# j = 0: d_q = sin(pi/5)/sin(pi/5) = 1
# j = 1/2: d_q = sin(2pi/5)/sin(pi/5) = phi
# j = 1: d_q = sin(3pi/5)/sin(pi/5) = phi (same!)
# j = 3/2: d_q = sin(4pi/5)/sin(pi/5) = 1

print(f"  Quantum dimensions at k+2 = a1 = 5:")
for j2 in range(4):  # j = 0, 1/2, 1, 3/2
    j = j2 / 2
    dq = np.sin((2*j+1)*np.pi/a1) / np.sin(np.pi/a1)
    dq_galois = np.sin((2*j+1)*np.pi/a1) / np.sin(np.pi/a1)
    # Under Galois: sin(pi/5) involves sqrt(5+sqrt(5))/...
    # Actually: sin(pi/5) = sqrt(10-2*sqrt(5))/4
    # sigma: sqrt(5) -> -sqrt(5), so sin(pi/5) -> sin(pi/5)? No...
    # The quantum dimensions at k+2=5 are: 1, phi, phi, 1
    # Under sigma: phi -> phi' = -1/phi (magnitude 1/phi)
    # So |d_q'| = 1, 1/phi, 1/phi, 1

    dq_dark = abs(PHI_CONJ) if abs(dq - PHI) < 0.01 else dq
    print(f"    j={j:.1f}: d_q = {dq:.6f}, d_q(dark) = {dq_dark:.6f}")

# The WEIGHTED vacuum energy would use d_q^2 as weights:
# E_vac_weighted(physical) = (1/D^2) * sum_rho d_q(rho)^2 * E_0(rho)
# E_vac_weighted(dark)     = (1/D'^2) * sum_rho d_q'(rho)^2 * E_0(rho)

# Where D^2 = sum d_q^2 = a1 + sqrt(a1) = 5 + sqrt(5) (physical)
# And D'^2 = sum d_q'^2 = a1 - sqrt(a1) = 5 - sqrt(5) (dark)

D2_phys = a1 + np.sqrt(a1)
D2_dark = a1 - np.sqrt(a1)

print(f"\n  D^2(physical) = a1 + sqrt(a1) = {D2_phys:.6f}")
print(f"  D^2(dark) = a1 - sqrt(a1) = {D2_dark:.6f}")
print(f"  Ratio D^2(phys)/D^2(dark) = {D2_phys/D2_dark:.6f}")
print(f"  = (a1+sqrt(a1))/(a1-sqrt(a1)) = {(a1+np.sqrt(a1))/(a1-np.sqrt(a1)):.6f}")

# Rationalize: (5+s5)/(5-s5) = (5+s5)^2 / (25-5) = (30+10s5)/20 = (3+s5)/2 = phi^2
ratio_D2 = D2_phys / D2_dark
print(f"  = phi^2 = {PHI**2:.6f}")
print(f"  Match: {abs(ratio_D2 - PHI**2) < 0.001}")

print(f"\n  D^2(physical) / D^2(dark) = phi^2 = {PHI**2:.6f}")
print(f"  This is the TQFT quantum dimension ratio!")

# ============================================================
# STEP 3: What ratios give the DM abundance?
# ============================================================
print(f"\n{'='*70}")
print("STEP 3: GALOIS RATIOS VS DM ABUNDANCE")
print(f"{'='*70}")

print(f"\n  Observed: Omega_DM/Omega_b = {Omega_DM_over_b}")
print(f"\n  Framework ratios:")

ratios = {
    'phi': PHI,
    'phi^2': PHI**2,
    'phi^3': PHI**3,
    'phi^4': PHI**4,
    'D^2_phys/D^2_dark = phi^2': PHI**2,
    'a1': float(a1),
    'a1 + 1/N_gen = 16/3': 16/3,
    '7 - phi': 7 - PHI,
    'b1 - 1/phi': b1 - 1/PHI,
    'sqrt(a1) * phi': np.sqrt(a1) * PHI,
    '2*phi + 2': 2*PHI + 2,
    'a1 * phi - 3': a1 * PHI - 3,
    'E_vac_phys/E_vac_dark (unweighted)': 1.0,
    'N_gen * phi': N_gen * PHI,
    'D^4_phys/D^4_dark = phi^4': PHI**4,
    'sum_dq^2(phys) / sum_dq^2(dark)': D2_phys / D2_dark,
}

print(f"\n  {'Ratio':>35} {'Value':>10} {'Error':>10} {'Match?':>8}")
for name, val in sorted(ratios.items(), key=lambda x: abs(x[1] - Omega_DM_over_b)):
    err = (val/Omega_DM_over_b - 1) * 100
    match = abs(err) < 2
    print(f"  {name:>35} {val:10.4f} {err:+10.2f}% {'<---' if match else ''}")

# ============================================================
# STEP 4: The 7-phi formula
# ============================================================
print(f"\n{'='*70}")
print("STEP 4: WHY 7-phi?")
print(f"{'='*70}")

val_7phi = 7 - PHI
print(f"\n  7 - phi = {val_7phi:.6f}")
print(f"  Observed: {Omega_DM_over_b}")
print(f"  Error: {(val_7phi/Omega_DM_over_b - 1)*100:+.2f}%")

# 7-phi = 7 - (1+sqrt(5))/2 = (13-sqrt(5))/2
# Galois conjugate: 7-phi' = 7-(1-sqrt(5))/2 = (13+sqrt(5))/2
print(f"\n  7-phi = (13-sqrt(5))/2")
print(f"  sigma(7-phi) = 7-phi' = (13+sqrt(5))/2 = {(13+SQRT5)/2:.6f}")
print(f"  Product: (7-phi)(7-phi') = (169-5)/4 = {164/4} = 41")
print(f"  Sum: (7-phi)+(7-phi') = 13")
print(f"  N(7-phi) = 41 (prime)")

# 7 = ? in the framework
print(f"\n  7 in the framework:")
print(f"    7 = a1 + 2 = 5 + 2")
print(f"    7 = E8 exponent m_2")
print(f"    7 = CKM exponent n_23")
print(f"    7 = degree - a1 = 12 - 5")
print(f"    7 = number of DERIVED dark properties")

# Can 7-phi come from TQFT quantum dimensions?
# D^2 = a1 + sqrt(a1) = 5 + sqrt(5)
# 7-phi = (13-sqrt(5))/2

# Check: is 7-phi = f(D^2, D'^2)?
# D^2 = 5+sqrt(5), D'^2 = 5-sqrt(5)
# D^2 + D'^2 = 10, D^2 * D'^2 = 20, D^2 - D'^2 = 2*sqrt(5)
# 7-phi = (13-sqrt(5))/2

# (D^2 + D'^2 + N_gen) / 2 = (10+3)/2 = 13/2 = 6.5. Not 7-phi.
# (D^2 + D'^2 + N_gen - sqrt(5)) / 2 = (10+3-2.236)/2 = 5.382. YES!

test1 = (D2_phys + D2_dark + N_gen - SQRT5) / 2
print(f"\n  (D^2+D'^2+N_gen-sqrt(5))/2 = ({D2_phys:.4f}+{D2_dark:.4f}+{N_gen}-{SQRT5:.4f})/2 = {test1:.6f}")
print(f"  7-phi = {val_7phi:.6f}")
print(f"  Match: {abs(test1 - val_7phi) < 0.001}")

# Simplify: D^2+D'^2 = 2*a1 = 10
# So: (2*a1 + N_gen - sqrt(5))/2 = (10+3-sqrt(5))/2 = (13-sqrt(5))/2 = 7-phi
# This is TRIVIALLY true since sqrt(5) = 2*phi-1 and 7-phi = (13-(2*phi-1))/2 = (14-2*phi)/2 = 7-phi.

# The decomposition (2*a1 + N_gen - sqrt(a1))/2 is just a REWRITING.
# Not a derivation.

# ============================================================
# STEP 5: Mass budget approach
# ============================================================
print(f"\n{'='*70}")
print("STEP 5: MASS BUDGET")
print(f"{'='*70}")

# Dark fermion masses: m_f' = m_e^2 / m_f = m_e * phi^(-n_f)
# Physical fermion masses: m_f = m_e * phi^(n_f)

# The DM abundance depends on:
# Omega_DM / Omega_b = (sum dark masses) / (proton mass)  * (n_DM / n_b)
# For freeze-in: n_DM/n_b depends on production rate

# But: the MASS RATIO of dark to physical sectors:
# sum m_f' / sum m_f = ?

m_e = 0.511  # MeV
fermion_n = {
    'e': 0, 'mu': 11, 'tau': 17,
    'u': 3, 'c': 16, 't': 26,
    'd': 5, 's': 11, 'b': 19
}

sum_phys = sum(m_e * PHI**n for n in fermion_n.values())
sum_dark = sum(m_e * PHI**(-n) for n in fermion_n.values())

print(f"\n  Sum of physical masses: {sum_phys:.2f} MeV")
print(f"  Sum of dark masses: {sum_dark:.6f} MeV")
print(f"  Ratio physical/dark: {sum_phys/sum_dark:.4f}")
print(f"  Ratio dark/physical: {sum_dark/sum_phys:.6e}")

# The dark sector is MUCH lighter. The ratio is dominated by top:
print(f"\n  Dominant physical: top = {m_e * PHI**26:.0f} MeV")
print(f"  Dominant dark: e' = {m_e:.3f} MeV")
print(f"  Lightest dark: t' = {m_e * PHI**(-26)*1e6:.4f} eV")

# The NUMBER density ratio determines the abundance.
# If n_DM = n_b (equal numbers), then:
# Omega_DM/Omega_b = m_DM_eff / m_proton
# For m_DM_eff ~ m_e (dominated by e'): Omega_DM/Omega_b ~ 0.511/938 ~ 0.0005
# Way too small!

# So we need n_DM >> n_b. The production mechanism matters.

# ============================================================
# STEP 6: Galois vacuum energy RATIO
# ============================================================
print(f"\n{'='*70}")
print("STEP 6: VACUUM ENERGIES FROM EXP444")
print(f"{'='*70}")

# From exp444 (in memory):
# E_vac(physical) = phi^3 = 4.236
# E_vac(dark) = phi'^3 = (1-sqrt(5))^3/8 = ...
# phi'^3 = -1/phi^3 (since phi*phi' = -1, so (phi*phi')^3 = -1)
# Product: phi^3 * phi'^3 = -1 EXACT

E_phys = PHI**3
E_dark = PHI_CONJ**3  # = -1/phi^3

print(f"  E_vac(physical) = phi^3 = {E_phys:.6f}")
print(f"  E_vac(dark) = phi'^3 = {E_dark:.6f}")
print(f"  Product = {E_phys * E_dark:.6f} (should be -1)")
print(f"  |E_dark/E_phys| = {abs(E_dark/E_phys):.6f} = 1/phi^6 = {1/PHI**6:.6f}")
print(f"  phi^6 = {PHI**6:.6f}")

# The energy ASYMMETRY:
# (E_phys - |E_dark|) / (E_phys + |E_dark|) = (phi^3 - 1/phi^3)/(phi^3 + 1/phi^3)
asym = (E_phys - abs(E_dark)) / (E_phys + abs(E_dark))
print(f"\n  Asymmetry: (E_p - |E_d|) / (E_p + |E_d|) = {asym:.6f}")

# The ratio E_phys / |E_dark| = phi^6
# phi^6 = 8*phi + 5 = 17.944
print(f"  E_phys / |E_dark| = phi^6 = {PHI**6:.4f}")

# What combinations give 5.36?
print(f"\n  phi^3 - 1 = {PHI**3 - 1:.4f}")
print(f"  phi^3 + 1/phi^3 = {PHI**3 + 1/PHI**3:.4f} = 2*(phi^3+phi'^3)/2... = {2*(E_phys+abs(E_dark))/2:.4f}")
print(f"  Actually phi^3 + 1/phi^3 = sqrt(a1)*(a1+2)/a1... let me just check numerically")

# Direct check: phi^3 = 2+sqrt(5), 1/phi^3 = sqrt(5)-2
# phi^3 + 1/phi^3 = 2*sqrt(5) = 2*sqrt(a1)
sum_E = PHI**3 + 1/PHI**3
print(f"  phi^3 + 1/phi^3 = {sum_E:.6f} = 2*sqrt(5) = {2*SQRT5:.6f}")

# phi^3 - 1/phi^3 = (2+sqrt(5)) - (sqrt(5)-2) = 4
diff_E = PHI**3 - 1/PHI**3
print(f"  phi^3 - 1/phi^3 = {diff_E:.6f} = 4 = d_ST")

# Interesting! phi^3 - 1/phi^3 = 4 EXACT.
# And phi^3 + 1/phi^3 = 2*sqrt(5) EXACT.

# So: E_phys - |E_dark| = 4 = d_ST (spacetime dimension!)
# And: E_phys + |E_dark| = 2*sqrt(a1)

# Can we get 7-phi from these?
# (E_phys - |E_dark|) + ... ?

print(f"\n  Combinations of E_vac with DM abundance:")
print(f"    phi^3 = {PHI**3:.6f}")
print(f"    phi^3 - 1/phi = {PHI**3 - 1/PHI:.6f}")
print(f"    phi^3 + 1/phi = {PHI**3 + 1/PHI:.6f}")
print(f"    phi^3 / phi = phi^2 = {PHI**2:.6f}")
print(f"    phi^3 * (1 - 1/phi^6) = {PHI**3 * (1-1/PHI**6):.6f}")
print(f"    (phi^6-1)/phi^3 = phi^3 - 1/phi^3 = {diff_E:.6f}")
print(f"    phi^3 - 1 = {PHI**3-1:.6f} = 1+sqrt(5) = {1+SQRT5:.6f}")
print(f"    2*phi + 1 = {2*PHI+1:.6f} = sqrt(5)+2 = phi^3")
print(f"    sqrt(5) + 2 = {SQRT5+2:.6f}")

# None of these give 5.36 directly.

# ============================================================
# STEP 7: Spectral approach - number of modes
# ============================================================
print(f"\n{'='*70}")
print("STEP 7: MODE COUNTING")
print(f"{'='*70}")

# In the 600-cell spectrum:
# Galois-broken modes: 4 + 9 + 9 + 4 = 26 = a1^2 + 1 (!)
# Galois-fixed modes: 1 + 16 + 25 + 36 + 16 = 94

galois_broken = 4 + 9 + 9 + 4  # = 26
galois_fixed = 1 + 16 + 25 + 36 + 16  # = 94

print(f"  Galois-broken modes: {galois_broken} = a1^2 + 1 = {a1**2+1}")
print(f"  Galois-fixed modes: {galois_fixed}")
print(f"  Total: {galois_broken + galois_fixed} = N = {N}")

print(f"\n  Ratio fixed/broken = {galois_fixed}/{galois_broken} = {galois_fixed/galois_broken:.6f}")
print(f"  = 47/13 = {47/13:.6f}")

# The DM abundance might be related to the FRACTION of Galois-broken modes
fraction_broken = galois_broken / N
print(f"  Fraction broken = {galois_broken}/{N} = {fraction_broken:.6f}")
print(f"  = 13/60 = {13/60:.6f}")

# Or the ratio of broken to total minus broken:
print(f"  broken / (total - broken) = {galois_broken}/{galois_fixed} = {galois_broken/galois_fixed:.6f}")

# What about: (N - galois_broken) / galois_broken = 94/26
print(f"  (N - broken)/broken = {galois_fixed/galois_broken:.6f}")

# None of these give 5.36 directly.

# But: 26 = a1^2 + 1 = sin^2(tW) denominator!
# And the DM abundance involves the DARK sector coupling structure.

# Try: N / broken = 120/26 = 60/13 = 4.615
print(f"\n  N / broken = {N/galois_broken:.6f}")
print(f"  broken * phi / degree = {galois_broken * PHI / degree:.6f}")
print(f"  a1 + 1/(phi-1) = a1 + phi = {a1 + PHI:.6f}") # = 6.618
print(f"  a1 + phi^(-1) = {a1 + 1/PHI:.6f}") # = 5.618
print(f"  a1 * phi^(-1) + b1*phi^(-2) = {a1/PHI + b1/PHI**2:.6f}") # checking
print(f"  broken / a1 = {galois_broken/a1:.6f} = {26/5}")

# ============================================================
# STEP 8: The Galois heat kernel at t=1
# ============================================================
print(f"\n{'='*70}")
print("STEP 8: GALOIS HEAT KERNEL AT NATURAL TIMES")
print(f"{'='*70}")

# From exp513: the Galois heat kernel Delta_S(t) at various t.
# Maybe at some NATURAL time (t=1, t=1/degree, t=phi, etc.)
# the ratio Delta_S / S_total gives the DM abundance?

L_1 = 12 - 6*PHI
L_8 = 6 + 6*PHI
L_2 = 12 - 2*SQRT5
L_6 = 12 + 2*SQRT5

# All eigenvalues with multiplicities
all_evals = [(0, 1), (L_1, 4), (L_2, 9), (9, 16), (12, 25),
             (14, 36), (L_6, 9), (15, 16), (L_8, 4)]

print(f"\n  {'t':>10} {'S_total':>12} {'Delta_S':>12} {'|Delta/S|':>12} {'vs 5.36':>10}")

for t in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 1/PHI, PHI_CONJ**2, 1/degree,
          1/a1, 1/b1, 1/(2*np.pi)]:
    S_total = sum(m * np.exp(-ev*t) for ev, m in all_evals)
    delta_S = 4*(np.exp(-L_1*t) - np.exp(-L_8*t)) + 9*(np.exp(-L_2*t) - np.exp(-L_6*t))

    if abs(delta_S) > 1e-15 and S_total > 1e-15:
        ratio = abs(delta_S / S_total)
        inv_ratio = S_total / abs(delta_S) if abs(delta_S) > 0 else float('inf')
        err = (inv_ratio / Omega_DM_over_b - 1) * 100
        print(f"  {t:10.4f} {S_total:12.4f} {delta_S:12.4f} {ratio:12.6f} {err:+10.2f}%")
    else:
        print(f"  {t:10.4f} {S_total:12.4e} {delta_S:12.4e} {'---':>12} {'---':>10}")

# ============================================================
# STEP 9: Direct DM formula search
# ============================================================
print(f"\n{'='*70}")
print("STEP 9: DIRECT FORMULA SEARCH")
print(f"{'='*70}")

# Target: 5.36 = Omega_DM / Omega_b
# What expressions in the framework give ~5.36?

E_vac_exact = Fraction(-539, 43200)

print(f"\n  Searching for expressions giving {Omega_DM_over_b}:")
print(f"\n  {'Expression':>40} {'Value':>10} {'Error':>8}")

expressions = {
    '7 - phi': 7 - PHI,
    'a1 + 1/N_gen': a1 + Fraction(1, N_gen),
    'b1 - 1/phi': b1 - 1/PHI,
    '2*phi^2': 2*PHI**2,
    'a1/phi + 2/a1': a1/PHI + 2/a1,
    'phi^3 + 1': PHI**3 + 1,
    '(a1^2+1)/a1': (a1**2+1)/a1,
    'sqrt(29)': np.sqrt(29),
    '4*phi - 1': 4*PHI - 1,
    '3*phi + 1/phi': 3*PHI + 1/PHI,
    'N_gen*phi + 1/(N_gen*phi)': N_gen*PHI + 1/(N_gen*PHI),
    'galois_broken/a1': galois_broken/a1,
    '(galois_broken-1)/a1': (galois_broken-1)/a1,
    'b1*phi^(-1) + a1*phi^(-2)': b1/PHI + a1/PHI**2,
    'phi^3 + phi^(-3)': PHI**3 + PHI**(-3),
    '2*sqrt(a1)': 2*SQRT5,
    'a1 - 1/phi^3': a1 - 1/PHI**3,
    'b1/phi + phi': b1/PHI + PHI,
    '1/|E_vac| * phi/degree': (1/abs(float(E_vac_exact))) * PHI/degree,
    '|E_vac|^(-1/phi)': abs(float(E_vac_exact))**(-1/PHI),
    '(h-1)/a1 - phi^(-2)': (30-1)/a1 - 1/PHI**2,
}

for name, val in sorted(expressions.items(), key=lambda x: abs(float(x[1]) - Omega_DM_over_b)):
    val_f = float(val)
    err = (val_f/Omega_DM_over_b - 1) * 100
    marker = " <---" if abs(err) < 1.5 else ""
    print(f"  {name:>40} {val_f:10.4f} {err:+8.2f}%{marker}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
RESULTS:
========

1. E_vac(physical) = E_vac(dark) = -539/43200 (spectrum is Galois-invariant).
   The unweighted orbifold vacuum energy does NOT distinguish sectors.

2. TQFT quantum dimension ratio: D^2(phys)/D^2(dark) = phi^2 = {PHI**2:.4f}
   Close to target {Omega_DM_over_b} but NOT matching ({(PHI**2/Omega_DM_over_b-1)*100:+.1f}%).

3. Galois-broken modes: {galois_broken} = a1^2 + 1 = {a1**2+1} out of {N} total.
   The 26 broken modes are the same 26 that appear in sin^2(tW) = 6/26.

4. Best formulas for DM abundance remain:
   - 7-phi = {7-PHI:.4f} ({(val_7phi/Omega_DM_over_b-1)*100:+.2f}%)
   - 16/3 = {16/3:.4f} ({(16/3/Omega_DM_over_b-1)*100:+.2f}%)
   - 2*sqrt(5) = {2*SQRT5:.4f} ({(2*SQRT5/Omega_DM_over_b-1)*100:+.2f}%)

   BUT none are DERIVED. They are all PATTERNS.

5. KEY FINDING: phi^3 - 1/phi^3 = 4 = d_ST EXACT.
   phi^3 + 1/phi^3 = 2*sqrt(5) EXACT.
   Product E_phys * E_dark = -1 EXACT.
   These relate vacuum energies to spacetime dimension.

STATUS: Omega_DM/Omega_b remains PATTERN.
  The orbifold/Galois structure does not directly give the abundance.
  The production mechanism (gravitational freeze-in) is the bottleneck.
""")

print("=" * 70)
print("EXP-514 COMPLETE")
print("=" * 70)

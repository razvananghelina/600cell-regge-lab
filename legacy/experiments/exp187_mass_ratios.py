"""
EXP-187: Mass Ratios Between Generations
=========================================
Can we DERIVE the inter-generation mass ratios from 600-cell geometry?

CONTEXT:
- Generation Theorem: b={0,1,2} on a=1 line, unit powers phi^{0,2,3}
- Mass formula: m_f = m_e * phi^(a*d_eff + b*k_eff)
- d_eff = 5 + sin^2(theta_W), k_eff = 6 - 1/phi^4 (FITTING)
- Question: Can d_eff and k_eff be DERIVED from spectral data?

PARTS:
A) Pure phi^n vs experimental ratios (residual analysis)
B) Inter-generation ratios: pattern search
C) Attempt to derive d_eff, k_eff from graph invariants
D) Alternative: use unit powers directly (phi^0, phi^2, phi^3)
E) Koide-like relations from 600-cell
"""

import math
import sys

PHI = (1 + math.sqrt(5)) / 2
PHI_CONJ = (1 - math.sqrt(5)) / 2  # phi' = -1/phi
ALPHA = 7.2973525693e-3
ALPHA_INV = 1.0 / ALPHA
SIN2_TW = 0.23122

# PDG 2022 masses in MeV
M_E = 0.51099895
M_MU = 105.6583755
M_TAU = 1776.86
M_U = 2.16       # MS-bar at 2 GeV, large uncertainty
M_D = 4.67       # MS-bar at 2 GeV
M_S = 93.4       # MS-bar at 2 GeV
M_C = 1270.0     # MS-bar
M_B = 4180.0     # MS-bar
M_T = 172690.0   # pole mass (172.69 GeV)

# Assignment from the unified formula (a,b) pairs
FERMIONS = {
    'e':   {'mass': M_E,   'a': 0, 'b': 0, 'n': 0,  'family': 'lepton'},
    'mu':  {'mass': M_MU,  'a': 1, 'b': 1, 'n': 11, 'family': 'lepton'},
    'tau': {'mass': M_TAU, 'a': 1, 'b': 2, 'n': 17, 'family': 'lepton'},
    'u':   {'mass': M_U,   'a': 3, 'b': -2,'n': 3,  'family': 'up'},
    'c':   {'mass': M_C,   'a': 2, 'b': 1, 'n': 16, 'family': 'up'},
    't':   {'mass': M_T,   'a': 4, 'b': 1, 'n': 26, 'family': 'up'},
    'd':   {'mass': M_D,   'a': 1, 'b': 0, 'n': 5,  'family': 'down'},
    's':   {'mass': M_S,   'a': 1, 'b': 1, 'n': 11, 'family': 'down'},
    'b':   {'mass': M_B,   'a': -1,'b': 4, 'n': 19, 'family': 'down'},
}


def log_phi(x):
    """Logarithm base phi."""
    return math.log(abs(x)) / math.log(PHI)


print("=" * 70)
print("EXP-187: MASS RATIOS BETWEEN GENERATIONS")
print("=" * 70)

# =====================================================================
# PART A: Pure phi^n ratios vs experimental
# =====================================================================
print("\n" + "=" * 70)
print("PART A: Pure phi^n prediction vs experiment")
print("=" * 70)

print("\nFor each fermion: effective exponent = log_phi(m_f/m_e)")
print("Compare with integer level n = 5a + 6b")
print()
print(f"{'Fermion':8} {'m (MeV)':>12} {'log_phi':>10} {'n=5a+6b':>8} {'delta':>8} {'err%':>8}")
print("-" * 60)

for name, info in FERMIONS.items():
    m = info['mass']
    n = info['n']
    if m > M_E:
        eff_exp = log_phi(m / M_E)
    else:
        eff_exp = 0.0
    delta = eff_exp - n
    if n > 0:
        err = 100 * delta / n
    else:
        err = 0.0
    print(f"{name:8} {m:12.4f} {eff_exp:10.4f} {n:8d} {delta:8.4f} {err:8.2f}%")

# =====================================================================
# PART B: Inter-generation ratios
# =====================================================================
print("\n" + "=" * 70)
print("PART B: Inter-generation mass ratios")
print("=" * 70)

families = {
    'Leptons (e -> mu -> tau)': [('e', 'mu'), ('mu', 'tau'), ('e', 'tau')],
    'Up quarks (u -> c -> t)': [('u', 'c'), ('c', 't'), ('u', 't')],
    'Down quarks (d -> s -> b)': [('d', 's'), ('s', 'b'), ('d', 'b')],
}

for family_name, pairs in families.items():
    print(f"\n--- {family_name} ---")
    print(f"{'Ratio':10} {'m2/m1':>12} {'log_phi':>10} {'dn':>6} {'delta':>8}")
    print("-" * 50)
    for f1, f2 in pairs:
        m1 = FERMIONS[f1]['mass']
        m2 = FERMIONS[f2]['mass']
        n1 = FERMIONS[f1]['n']
        n2 = FERMIONS[f2]['n']
        ratio = m2 / m1
        eff = log_phi(ratio)
        dn = n2 - n1
        delta = eff - dn
        print(f"{f2}/{f1:3}    {ratio:12.4f} {eff:10.4f} {dn:6d} {delta:8.4f}")

# =====================================================================
# PART C: Derive d_eff and k_eff from graph invariants
# =====================================================================
print("\n" + "=" * 70)
print("PART C: Deriving d_eff and k_eff")
print("=" * 70)

# Current fitting values
d_eff_fit = 5 + SIN2_TW  # = 5.23122
k_eff_fit = 6 - 1.0/PHI**4  # = 6 - 0.14590... = 5.85410

print(f"\nCurrent FITTING values:")
print(f"  d_eff = 5 + sin^2(theta_W) = 5 + {SIN2_TW:.5f} = {d_eff_fit:.5f}")
print(f"  k_eff = 6 - 1/phi^4 = 6 - {1/PHI**4:.5f} = {k_eff_fit:.5f}")

# What d_eff and k_eff would need to be for exact mass reproduction
# Using leptons (cleanest, no QCD complications):
# m_mu = m_e * phi^(d_eff + k_eff) => d_eff + k_eff = log_phi(m_mu/m_e)
# m_tau = m_e * phi^(d_eff + 2*k_eff) => d_eff + 2*k_eff = log_phi(m_tau/m_e)

exp_mu = log_phi(M_MU / M_E)
exp_tau = log_phi(M_TAU / M_E)

# Solve: d + k = exp_mu, d + 2k = exp_tau
k_exact = exp_tau - exp_mu
d_exact = exp_mu - k_exact

print(f"\nExact values from lepton masses:")
print(f"  log_phi(m_mu/m_e) = {exp_mu:.6f} = d_eff + k_eff")
print(f"  log_phi(m_tau/m_e) = {exp_tau:.6f} = d_eff + 2*k_eff")
print(f"  => k_eff = {k_exact:.6f}")
print(f"  => d_eff = {d_exact:.6f}")

# How close are these to graph invariants?
print(f"\nComparison with graph invariants:")
print(f"  d_eff = {d_exact:.6f} vs a_1 = 5, difference = {d_exact - 5:.6f}")
print(f"  k_eff = {k_exact:.6f} vs b_1 = 6, difference = {k_exact - 6:.6f}")

# Test various geometric expressions for the corrections
print(f"\nCandidate expressions for corrections:")
delta_d = d_exact - 5
delta_k = k_exact - 6

candidates_d = {
    'sin^2(tW) = 6/26':            6.0/26,
    'alpha':                         ALPHA,
    '3*alpha':                       3*ALPHA,
    'b1/(a1^2+1)':                  6.0/26,  # same as 6/26
    'phi - 1 = 1/phi':             1.0/PHI,
    '1/(2*phi^2)':                  0.5/PHI**2,
    '2*alpha*phi':                  2*ALPHA*PHI,
    'alpha*pi':                     ALPHA*math.pi,
    '1/phi^3':                      1.0/PHI**3,
    'alpha*phi^2':                  ALPHA*PHI**2,
}

candidates_k = {
    '-1/phi^4':                     -1.0/PHI**4,
    '-alpha':                       -ALPHA,
    '-phi+1 = -1/phi':             -1.0/PHI,
    '-1/(2*phi)':                   -0.5/PHI,
    '-2*alpha':                     -2*ALPHA,
    '-alpha*phi':                   -ALPHA*PHI,
    '-phi^(-2)':                    -1.0/PHI**2,
    '-1/(phi^2+1)':                 -1.0/(PHI**2+1),
    '-2/(phi^4+phi^2)':            -2.0/(PHI**4+PHI**2),
    '-(phi-1)^2':                   -(PHI-1)**2,
}

print(f"\n  delta_d (exact) = {delta_d:.6f}")
print(f"  {'Expression':30} {'Value':>10} {'Error':>10}")
print(f"  {'-'*55}")
for name, val in candidates_d.items():
    err = abs(val - delta_d)
    print(f"  {name:30} {val:10.6f} {err:10.6f}")

print(f"\n  delta_k (exact) = {delta_k:.6f}")
print(f"  {'Expression':30} {'Value':>10} {'Error':>10}")
print(f"  {'-'*55}")
for name, val in candidates_k.items():
    err = abs(val - delta_k)
    print(f"  {name:30} {val:10.6f} {err:10.6f}")

# =====================================================================
# PART D: Use unit powers directly
# =====================================================================
print("\n" + "=" * 70)
print("PART D: Mass ratios from unit powers phi^k")
print("=" * 70)

# Generations have phi^0, phi^2, phi^3 as units on a=1 line
# The MASS ratio is NOT phi^0, phi^2, phi^3 directly
# But the LEVEL ratio is n = 5, 11, 17 (= 5 + 6*b)

# Key question: is there a function f such that m_gen/m_e = phi^f(k)?
# For k=0: m_e/m_e = 1 = phi^0 => f(0) = 0 (trivial)
# For k=2: m_mu/m_e = phi^11.06 => f(2) = 11.06
# For k=3: m_tau/m_e = phi^16.92 => f(3) = 16.92

# Ratios of phi-exponents:
# f(2)/f(0) = undefined (0/0)
# f(3)/f(2) = 16.92/11.06 = 1.530

f2 = exp_mu
f3 = exp_tau
print(f"\nUnit power exponents:")
print(f"  phi^0 -> level 0 (electron)")
print(f"  phi^2 -> effective exponent f(2) = {f2:.6f}")
print(f"  phi^3 -> effective exponent f(3) = {f3:.6f}")
print(f"  Ratio f(3)/f(2) = {f3/f2:.6f}")

# Test: f(k) = k * something?
print(f"\n  f(2)/2 = {f2/2:.6f}")
print(f"  f(3)/3 = {f3/3:.6f}")
# Not constant, so f is not linear in k

# Test: f(k) = a1*k + b1 (affine)?
# f(2) = 2a + b, f(3) = 3a + b
# => a = f(3) - f(2) = k_eff_exact, b = f(2) - 2*a = d_exact - a
a_affine = f3 - f2
b_affine = f2 - 2*a_affine
print(f"\n  Affine f(k) = {a_affine:.6f}*k + {b_affine:.6f}")
print(f"  Check: f(2) = {a_affine*2 + b_affine:.6f} (should be {f2:.6f})")
print(f"  Check: f(3) = {a_affine*3 + b_affine:.6f} (should be {f3:.6f})")
print(f"  Prediction f(0) = {b_affine:.6f} (should be 0)")
print(f"    => electron mass from affine: m_e * phi^{b_affine:.4f} = {M_E * PHI**b_affine:.4f} MeV")
print(f"       actual m_e = {M_E:.4f} MeV")

# This is just restating d_eff + k_eff decomposition.
# The REAL question: why these specific values?

# =====================================================================
# PART E: Alternative geometric derivations
# =====================================================================
print("\n" + "=" * 70)
print("PART E: Searching for geometric mass ratio formulas")
print("=" * 70)

# Approach 1: Eigenvalue ratios
# 600-cell eigenvalues: 12, 6phi, 4phi, 3, 0, -2, 4-4phi, -3, 6-6phi
eigenvalues = [12, 6*PHI, 4*PHI, 3, 0, -2, 4-4*PHI, -3, 6-6*PHI]
multiplicities = [1, 4, 9, 16, 25, 36, 9, 16, 4]

print("\n--- Approach 1: Eigenvalue ratios as mass ratios ---")
print(f"\nTarget: m_mu/m_e = {M_MU/M_E:.4f}")
print(f"Target: m_tau/m_mu = {M_TAU/M_MU:.4f}")
print(f"Target: m_tau/m_e = {M_TAU/M_E:.4f}")

# Test all ratios of positive eigenvalues
pos_eigs = [(i, e) for i, e in enumerate(eigenvalues) if e > 0]
print(f"\nPositive eigenvalues: {[(i, round(e, 4)) for i, e in pos_eigs]}")

print(f"\nRatios of eigenvalue products/powers:")
# Test phi^(eigenvalue) ratios
for i, ei in enumerate(eigenvalues):
    for j, ej in enumerate(eigenvalues):
        if i != j and ei > 0 and ej > 0:
            ratio = ei / ej
            if 10 < ratio < 300:
                # Close to m_mu/m_e?
                err_mu = abs(ratio - M_MU/M_E) / (M_MU/M_E) * 100
                if err_mu < 5:
                    print(f"  lambda_{i}/lambda_{j} = {ei:.4f}/{ej:.4f} = {ratio:.4f} ~ m_mu/m_e ({err_mu:.2f}%)")
            if 10 < ratio < 25:
                err_tau = abs(ratio - M_TAU/M_MU) / (M_TAU/M_MU) * 100
                if err_tau < 5:
                    print(f"  lambda_{i}/lambda_{j} = {ei:.4f}/{ej:.4f} = {ratio:.4f} ~ m_tau/m_mu ({err_tau:.2f}%)")

# Approach 2: Multiplicity-weighted expressions
print(f"\n--- Approach 2: Spectral invariants ---")
# Spectral zeta: sum of |lambda_i|^s * m_i
for s in [1, 2, -1, 0.5]:
    zeta = sum(abs(e)**s * m if e != 0 else 0
               for e, m in zip(eigenvalues, multiplicities))
    print(f"  zeta({s:4.1f}) = sum |lambda_i|^s * mult_i = {zeta:.4f}")

# Approach 3: Graph distance and mass
print(f"\n--- Approach 3: Level spacing from graph parameters ---")
# The levels on a=1 are: n = 5, 11, 17 (spacing = 6)
# But effective exponents are: 0, 11.06, 16.92 (spacing = 11.06, 5.86)
# The ratio 11.06/5.86 = 1.888... is close to 2
# Because (1,1) has a=1 contributing d_eff and b=1 contributing k_eff
# While (1,2)-(1,1) has only Delta_b = 1, contributing k_eff only

print(f"\n  Spacing analysis:")
print(f"  e -> mu:  Delta(a,b) = (1,1), delta_exp = {exp_mu:.4f}")
print(f"  mu -> tau: Delta(a,b) = (0,1), delta_exp = {exp_tau - exp_mu:.4f}")
print(f"  Ratio of spacings = {exp_mu / (exp_tau - exp_mu):.4f}")
print(f"  Compare with d_eff/k_eff = {d_exact/k_exact:.4f}")

# Approach 4: Norm and mass
print(f"\n--- Approach 4: Galois norm and mass scale ---")
for name in ['e', 'mu', 'tau', 'd', 's', 'b', 'u', 'c', 't']:
    info = FERMIONS[name]
    a, b = info['a'], info['b']
    z = a + b * PHI
    z_conj = a + b * PHI_CONJ
    norm = z * z_conj  # = a^2 + ab - b^2
    norm_int = a**2 + a*b - b**2
    print(f"  {name:5}: z = {a}+{b}*phi = {z:8.4f}, "
          f"N(z) = {norm_int:4d}, |N(z)| = {abs(norm_int):4d}, "
          f"m = {info['mass']:.2f} MeV")

# Approach 5: phi-adic valuation
print(f"\n--- Approach 5: Key number theory observation ---")
print(f"\nThe unit powers on a=1 are phi^{{0,2,3}}.")
print(f"These give b = {{0,1,2}}.")
print(f"The MASS exponent is NOT k, but rather n = 5a + 6b = 5+6b.")
print(f"So the question reduces to: why does m_f scale as phi^(5+6b)?")
print(f"Answer: because d_eff ~ 5 and k_eff ~ 6 are the intersection")
print(f"numbers a_1 and b_1 of the 600-cell distance-regular graph.")
print(f"")
print(f"The RESIDUALS (delta_d, delta_k) encode radiative corrections.")
print(f"For leptons:")
print(f"  delta_d = {delta_d:.6f}")
print(f"  delta_k = {delta_k:.6f}")

# Test: are residuals related to running of alpha?
# alpha(m_mu) ~ alpha/(1 - alpha*ln(m_mu/m_e)/(3*pi))
alpha_mu = ALPHA / (1 - ALPHA * math.log(M_MU/M_E) / (3*math.pi))
alpha_tau = ALPHA / (1 - ALPHA * math.log(M_TAU/M_E) / (3*math.pi))

print(f"\n--- Approach 6: Running coupling corrections ---")
print(f"  alpha(m_e)   = {ALPHA:.8f}")
print(f"  alpha(m_mu)  = {alpha_mu:.8f} (1-loop running)")
print(f"  alpha(m_tau) = {alpha_tau:.8f} (1-loop running)")

# The correction to the exponent from running would be:
# Instead of phi^n, we get phi^(n + correction)
# correction ~ alpha * something

# =====================================================================
# PART F: Can we derive delta_d and delta_k from first principles?
# =====================================================================
print("\n" + "=" * 70)
print("PART F: First-principles derivation attempt")
print("=" * 70)

print(f"\nThe exact lepton corrections are:")
print(f"  delta_d = {delta_d:.8f}")
print(f"  delta_k = {delta_k:.8f}")
print(f"  Sum d_eff + k_eff = {d_exact + k_exact:.8f}")

# Key observation: d_eff + k_eff = log_phi(m_mu/m_e)
# This is the MUON-ELECTRON ratio, one of the most precisely known

# Test: d_eff + k_eff = 11 + alpha*f(phi)?
residual_sum = (d_exact + k_exact) - 11
print(f"\n  d_eff + k_eff - 11 = {residual_sum:.8f}")
print(f"  = {residual_sum/ALPHA:.4f} * alpha")
print(f"  = {residual_sum/ALPHA/PHI:.4f} * alpha * phi")
print(f"  = {residual_sum/(1.0/PHI):.4f} * 1/phi")
print(f"  = {residual_sum/SIN2_TW:.4f} * sin^2(tW)")

# Test specific combinations
tests = {
    'alpha*(12+phi)':              ALPHA*(12+PHI),
    '2*alpha*pi*phi':              2*ALPHA*math.pi*PHI,
    'sin^2(tW) - 1/phi^4':        SIN2_TW - 1.0/PHI**4,
    '6/26 - 1/phi^4':             6.0/26 - 1.0/PHI**4,
    'alpha*(phi+8)':               ALPHA*(PHI+8),
    '1/phi^3 - 1/phi^4':          1.0/PHI**3 - 1.0/PHI**4,
    '1/phi^5':                     1.0/PHI**5,
    'phi^(-2) - phi^(-4)':        PHI**(-2) - PHI**(-4),
    '2/phi^3':                     2.0/PHI**3,
    'alpha_s (= 1/(2*phi^3))':    0.5/PHI**3,
    '(phi-1)^3':                   (PHI-1)**3,
}

print(f"\n  Target: delta_d + delta_k = {delta_d + delta_k:.8f}")
print(f"  {'Expression':30} {'Value':>12} {'Error':>12}")
print(f"  {'-'*58}")
for name, val in tests.items():
    err = abs(val - (delta_d + delta_k))
    marker = " <--" if err < 0.005 else ""
    print(f"  {name:30} {val:12.8f} {err:12.8f}{marker}")

# =====================================================================
# PART G: Koide-like analysis
# =====================================================================
print("\n" + "=" * 70)
print("PART G: Koide-like relations from 600-cell")
print("=" * 70)

# Koide formula: (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3
koide_num = M_E + M_MU + M_TAU
koide_den = (math.sqrt(M_E) + math.sqrt(M_MU) + math.sqrt(M_TAU))**2
koide = koide_num / koide_den
print(f"\nKoide ratio (leptons): {koide:.8f}")
print(f"  Exact 2/3:           {2/3:.8f}")
print(f"  Error:               {abs(koide - 2/3)/koide*100:.4f}%")

# Can we get 2/3 from 600-cell?
# 2/3 = 2*a_1 / (a_1^2 - a_1 + 6) ? No... let's check
# 2/3 could come from: number of cosets minus 2 / number of cosets = 3/5? No.
# Or: 8/(8+16) = 8/24 = 1/3... 2/3 = 16/24

# Koide-like for quarks
koide_u = (M_U + M_C + M_T) / (math.sqrt(M_U) + math.sqrt(M_C) + math.sqrt(M_T))**2
koide_d = (M_D + M_S + M_B) / (math.sqrt(M_D) + math.sqrt(M_S) + math.sqrt(M_B))**2
print(f"\nKoide ratio (up quarks):   {koide_u:.8f}")
print(f"Koide ratio (down quarks): {koide_d:.8f}")

# 600-cell version: use phi-weighted Koide
# Replace sqrt with phi-power: phi^(n/2)
# For leptons: n = 0, 11, 17
phi_koide_num = PHI**0 + PHI**11 + PHI**17
phi_koide_den = (PHI**0 + PHI**(11.0/2) + PHI**(17.0/2))**2
phi_koide = phi_koide_num / phi_koide_den
print(f"\nPhi-Koide (integer levels): {phi_koide:.8f}")

# Alternatively, use b-values (0,1,2)
# sqrt(m) ~ phi^(n/2), with n proportional to b
# Koide with b: (b=0: 1, b=1: phi^(k/2), b=2: phi^k)
# where k = k_eff
# This gives Koide = (1 + phi^k + phi^(2k)) / (1 + phi^(k/2) + phi^k)^2
# = (phi^(2k) + phi^k + 1) / (phi^k + phi^(k/2) + 1)^2

k = k_exact
geom_koide_num = PHI**(2*k) + PHI**k + 1
geom_koide_den = (PHI**k + PHI**(k/2) + 1)**2
geom_koide = geom_koide_num / geom_koide_den
print(f"Geometric Koide (k_eff spacing): {geom_koide:.8f}")

# What value of k gives EXACT 2/3?
# We need: 3*(x^2 + x + 1) = 2*(x + sqrt(x) + 1)^2
# where x = phi^k
# This is transcendental, solve numerically
print(f"\nSolving for k such that Koide = 2/3 exactly:")
from math import log

# Binary search
k_lo, k_hi = 4.0, 8.0
for _ in range(100):
    k_mid = (k_lo + k_hi) / 2
    x = PHI**k_mid
    sq = math.sqrt(x)
    koide_test = (x**2 + x + 1) / (x + sq + 1)**2
    if koide_test > 2.0/3:
        k_lo = k_mid
    else:
        k_hi = k_mid
k_koide = (k_lo + k_hi) / 2
print(f"  k_Koide = {k_koide:.6f}")
print(f"  k_eff   = {k_exact:.6f}")
print(f"  Difference: {abs(k_koide - k_exact):.6f}")

# =====================================================================
# PART H: Summary and status
# =====================================================================
print("\n" + "=" * 70)
print("PART H: SUMMARY")
print("=" * 70)

print(f"""
RESULTS:

1. PURE LEVEL STRUCTURE (n = 5a + 6b):
   - Integer levels 5, 11, 17 give phi^5, phi^11, phi^17
   - phi^11 = {PHI**11:.2f} vs m_mu/m_e = {M_MU/M_E:.2f} (error {abs(PHI**11 - M_MU/M_E)/(M_MU/M_E)*100:.1f}%)
   - phi^6 = {PHI**6:.2f} vs m_tau/m_mu = {M_TAU/M_MU:.2f} (error {abs(PHI**6 - M_TAU/M_MU)/(M_TAU/M_MU)*100:.1f}%)
   STATUS: DERIVED that n=5,11,17 from graph + Generation Theorem
   STATUS: But phi^n gives ~55% error (need corrections d_eff, k_eff)

2. EFFECTIVE EXPONENTS (from lepton masses):
   d_eff = {d_exact:.6f} (close to a_1 = 5)
   k_eff = {k_exact:.6f} (close to b_1 = 6)
   delta_d = {delta_d:.6f}
   delta_k = {delta_k:.6f}

3. BEST MATCHES FOR CORRECTIONS:
   delta_d ~ sin^2(theta_W) = 6/26 = {6.0/26:.6f} (err {abs(delta_d - 6.0/26):.6f})
   delta_k ~ -1/phi^4 = {-1.0/PHI**4:.6f} (err {abs(delta_k + 1.0/PHI**4):.6f})
   STATUS: FITTING (these match but are not derived)

4. KOIDE RELATION:
   Koide ratio = {koide:.8f} vs 2/3 = {2/3:.8f}
   600-cell Koide with k_eff: {geom_koide:.8f}
   k that gives exact 2/3: {k_koide:.6f} vs k_eff = {k_exact:.6f}
   STATUS: INTERESTING but not a derivation

5. HIERARCHY:
   - Generation Theorem: DERIVED (which levels are occupied)
   - Level spacing = 6 = b_1: DERIVED (graph parameter)
   - Corrections (delta_d, delta_k): FITTING
   - Mass ratios: PATTERN (d_eff = 5+sin^2(tW), k_eff = 6-1/phi^4)
""")

# What CAN we say is derived?
print("WHAT IS DERIVED vs WHAT IS NOT:")
print("-" * 50)
print("DERIVED:")
print("  - 3 generations (b=0,1,2)")
print("  - Level spacing proportional to b_1=6 (graph)")
print("  - Base of hierarchy = phi (DSI)")
print("  - Integer structure n = 5a + 6b")
print("")
print("NOT DERIVED (PATTERN/FITTING):")
print("  - Why d_eff = 5 + sin^2(theta_W)")
print("  - Why k_eff = 6 - 1/phi^4")
print("  - Absolute mass scale (m_e as input)")
print("  - Quark (a,b) assignments")

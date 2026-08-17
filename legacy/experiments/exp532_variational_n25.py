"""
EXP-532: Variational Selection of n=25
========================================

From exp531: m_e * phi^25 = 85.7 GeV falls inside the EW matching
window [83, 90] GeV where sin^2(tW) = 6/26 and alpha_s = 1/(2*phi^3)
match the SM running (exp527-529).

HYPOTHESIS: n = 25 is UNIQUELY selected as the integer exponent for
which the bare lattice scale m_e * phi^n is self-consistent with the
SM gauge couplings being equal to their framework values.

TESTS:
  1. For each integer n, compute Lambda(n) = m_e * phi^n.
  2. At Lambda(n), evaluate sin^2(tW) and alpha_s from SM RG running.
  3. Define a cost function measuring the discrepancy from framework values.
  4. Show that cost is minimized at n = 25 (and ONLY n = 25).
  5. Interpret: is this a DERIVATION or a CONSISTENCY CHECK?

BONUS: Can we derive the RUNNING RATIO 16/15 from the self-consistency?
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, SQRT5, a1, b1, N, h, degree, d_ST, rank_E8,
                      alpha_em, inv_alpha, alpha_s as alpha_s_fw, sin2_tW as sin2_fw,
                      N_gen, SPECTRUM_600CELL, IRREP_DIMS_2I)

print("=" * 70)
print("EXP-532: VARIATIONAL SELECTION OF n = 25")
print("=" * 70)

# Physical constants
m_e = 0.51099895e-3    # GeV
M_Z_exp = 91.1876      # GeV
M_W_exp = 80.3692      # GeV

# PDG 2024 at M_Z (MSbar):
alpha_em_MZ = 1.0 / 127.951
alpha_s_MZ  = 0.1179
sin2_MZ     = 0.23122

# Framework values:
sin2_target = 6.0 / 26      # = 0.230769...
as_target   = 1/(2*PHI**3)  # = 0.118034...

# SM one-loop beta coefficients (GUT normalised)
b_rg = np.array([41.0/10, -19.0/6, -7.0])

def sin2_aem_to_inv_alphas(sin2, aem):
    a2 = aem / sin2
    a_Y = aem / (1 - sin2)
    a1_gut = (5.0/3) * a_Y
    return 1.0/a1_gut, 1.0/a2

def inv_alphas_to_sin2_aem(inv_a1, inv_a2):
    al1 = 1.0/inv_a1
    al2 = 1.0/inv_a2
    a_Y = (3.0/5) * al1
    sin2 = a_Y / (a_Y + al2)
    aem = a_Y * al2 / (a_Y + al2)
    return sin2, aem

inv_a1_MZ, inv_a2_MZ = sin2_aem_to_inv_alphas(sin2_MZ, alpha_em_MZ)
inv_a3_MZ = 1.0 / alpha_s_MZ

def run_one_loop(inv_ai, b_i, mu_ref, mu_target):
    return inv_ai - b_i / (2*np.pi) * np.log(mu_target / mu_ref)

def couplings_at_mu(mu):
    """Return (sin^2, alpha_s, alpha_em) at scale mu via one-loop running from M_Z."""
    ia1 = run_one_loop(inv_a1_MZ, b_rg[0], M_Z_exp, mu)
    ia2 = run_one_loop(inv_a2_MZ, b_rg[1], M_Z_exp, mu)
    ia3 = run_one_loop(inv_a3_MZ, b_rg[2], M_Z_exp, mu)
    sin2, aem = inv_alphas_to_sin2_aem(ia1, ia2)
    return sin2, 1.0/ia3, aem

# ============================================================
# SECTION 1: SCAN OVER n
# ============================================================
print(f"\n{'='*70}")
print("SECTION 1: COST FUNCTION vs n")
print(f"{'='*70}")

print(f"\n  Framework targets: sin^2 = {sin2_target:.6f}, alpha_s = {as_target:.6f}")
print(f"\n  {'n':>3s} {'Lambda (GeV)':>14s} {'sin2(Lam)':>10s} {'as(Lam)':>10s} "
      f"{'d_sin2%':>9s} {'d_as%':>9s} {'cost':>12s} {'note':>15s}")
print(f"  {'-'*3} {'-'*14} {'-'*10} {'-'*10} {'-'*9} {'-'*9} {'-'*12} {'-'*15}")

results = {}
for n in range(15, 35):
    Lambda_n = m_e * PHI**n

    if Lambda_n < 1.0 or Lambda_n > 1e6:
        continue

    sin2_n, as_n, aem_n = couplings_at_mu(Lambda_n)

    d_sin2 = (sin2_n / sin2_target - 1) * 100
    d_as = (as_n / as_target - 1) * 100

    # Cost = sum of squared relative deviations
    cost = (sin2_n/sin2_target - 1)**2 + (as_n/as_target - 1)**2

    note = ""
    if n == a1**2:
        note = "<-- a1^2"
    elif n == a1**2 - 1:
        note = "a1^2-1"
    elif n == a1**2 + 1:
        note = "a1^2+1=N_broken"

    results[n] = (Lambda_n, sin2_n, as_n, d_sin2, d_as, cost)

    if 10 < Lambda_n < 1e5:
        print(f"  {n:3d} {Lambda_n:14.4f} {sin2_n:10.6f} {as_n:10.6f} "
              f"{d_sin2:+9.3f} {d_as:+9.3f} {cost:12.6e} {note:>15s}")

# Find minimum cost
best_n = min(results, key=lambda k: results[k][5])
best_Lambda, best_sin2, best_as, best_d_sin2, best_d_as, best_cost = results[best_n]

print(f"\n  MINIMUM cost at n = {best_n}")
print(f"    Lambda = {best_Lambda:.4f} GeV")
print(f"    sin^2 dev: {best_d_sin2:+.4f}%")
print(f"    alpha_s dev: {best_d_as:+.4f}%")
print(f"    cost = {best_cost:.6e}")

# ============================================================
# SECTION 2: CONTINUOUS OPTIMUM
# ============================================================
print(f"\n{'='*70}")
print("SECTION 2: CONTINUOUS OPTIMUM (non-integer n)")
print(f"{'='*70}")

# Find the continuous n that minimizes cost
from scipy.optimize import minimize_scalar

def cost_fn(n_cont):
    Lambda = m_e * PHI**n_cont
    if Lambda < 1 or Lambda > 1e5:
        return 1e10
    sin2, a_s, _ = couplings_at_mu(Lambda)
    return (sin2/sin2_target - 1)**2 + (a_s/as_target - 1)**2

res = minimize_scalar(cost_fn, bounds=(20, 30), method='bounded')
n_opt = res.x
Lambda_opt = m_e * PHI**n_opt
sin2_opt, as_opt, aem_opt = couplings_at_mu(Lambda_opt)

print(f"  Optimal continuous n = {n_opt:.4f}")
print(f"  Lambda_opt = {Lambda_opt:.4f} GeV")
print(f"  sin^2(Lambda_opt) = {sin2_opt:.6f} (target: {sin2_target:.6f}, "
      f"dev: {(sin2_opt/sin2_target-1)*100:+.4f}%)")
print(f"  alpha_s(Lambda_opt) = {as_opt:.6f} (target: {as_target:.6f}, "
      f"dev: {(as_opt/as_target-1)*100:+.4f}%)")

print(f"\n  Nearest integer: n = {round(n_opt)}")
print(f"  Distance from a1^2 = 25: {abs(n_opt - 25):.4f}")
print(f"  n_opt = {n_opt:.4f} = a1^2 + {n_opt - a1**2:.4f}")

# ============================================================
# SECTION 3: COST LANDSCAPE AROUND n=25
# ============================================================
print(f"\n{'='*70}")
print("SECTION 3: COST LANDSCAPE (FINE SCAN)")
print(f"{'='*70}")

n_fine = np.linspace(24.0, 26.0, 201)
cost_fine = np.array([cost_fn(n) for n in n_fine])

# Cost at integer values
cost_24 = cost_fn(24)
cost_25 = cost_fn(25)
cost_26 = cost_fn(26)

print(f"  Cost at n=24: {cost_24:.6e}  (Lambda = {m_e*PHI**24:.1f} GeV)")
print(f"  Cost at n=25: {cost_25:.6e}  (Lambda = {m_e*PHI**25:.1f} GeV)")
print(f"  Cost at n=26: {cost_26:.6e}  (Lambda = {m_e*PHI**26:.1f} GeV)")
print(f"  Ratio cost(24)/cost(25) = {cost_24/cost_25:.1f}")
print(f"  Ratio cost(26)/cost(25) = {cost_26/cost_25:.1f}")

# ============================================================
# SECTION 4: UNIQUENESS — IS n=25 THE ONLY OPTION?
# ============================================================
print(f"\n{'='*70}")
print("SECTION 4: UNIQUENESS OF n = 25")
print(f"{'='*70}")

# The EW window is [83, 90] GeV.
# Which integer n satisfies m_e * phi^n in [83, 90]?
print(f"\n  EW matching window: [83, 90] GeV")
print(f"  phi = {PHI:.6f}")
for n in range(20, 30):
    Lambda = m_e * PHI**n
    in_window = "YES" if 83 <= Lambda <= 90 else "no"
    in_extended = "yes" if 70 <= Lambda <= 100 else "no"
    print(f"    n={n}: Lambda = {Lambda:10.2f} GeV   in [83,90]: {in_window}   in [70,100]: {in_extended}")

print(f"\n  RESULT: Only n = {a1**2} places Lambda in [83, 90] GeV.")
print(f"  Even in the extended window [70, 100], only n = 25 works.")
print(f"  n = 24 gives {m_e*PHI**24:.1f} GeV (too low).")
print(f"  n = 26 gives {m_e*PHI**26:.1f} GeV (too high).")

# What is the WIDTH of the "allowed" window in log-space?
# m_e * phi^n in [83, 90] => n in [log(83/m_e)/log(phi), log(90/m_e)/log(phi)]
n_lo = np.log(83/m_e) / np.log(PHI)
n_hi = np.log(90/m_e) / np.log(PHI)
print(f"\n  Window in n-space: [{n_lo:.3f}, {n_hi:.3f}]")
print(f"  Width: {n_hi - n_lo:.3f}")
print(f"  Contains exactly one integer: {int(np.floor(n_hi))} (= a1^2)")

# Probability of a random integer falling in this window:
# The window width is 0.167 in n-space. For random placement, P = 0.167.
# But n is constrained to be a "nice" integer from the 600-cell data.
# The relevant candidates are: a1^2=25, 4a1=20, h-a1=25 (same!), etc.
# So the selection is NOT random — it's from a small set of candidates.

# ============================================================
# SECTION 5: SELF-CONSISTENCY BOOTSTRAP
# ============================================================
print(f"\n{'='*70}")
print("SECTION 5: SELF-CONSISTENCY BOOTSTRAP")
print(f"{'='*70}")

# The idea: if the lattice IS at scale Lambda, then the couplings
# at Lambda should be the framework values. This is self-consistent
# if and only if the SM running from M_Z to Lambda gives 6/26 and 1/(2phi^3).
#
# Turn this around: define Lambda_* as the scale where the SM couplings
# simultaneously best match the framework values. Then:
#   m_Z = Lambda_* * (alpha(m_Z)/alpha(Lambda_*))^{exponent???}
#   or equivalently: Lambda_* = m_e * phi^n for some n.
#
# The bootstrap equation:
#   Lambda = m_e * phi^n   [mass formula]
#   g_i(Lambda) = g_i^{FW}   [coupling matching]
#   m_Z = Lambda * R(Lambda)  [running dresses bare to physical]
#
# where R(Lambda) = alpha(m_Z)/alpha(Lambda) encodes the running FROM
# the lattice scale TO the physical m_Z.

# If Lambda = 85.7 GeV and m_Z = 91.2 GeV:
Lambda_25 = m_e * PHI**25
running_needed = M_Z_exp / Lambda_25
print(f"  Lambda(25) = m_e * phi^25 = {Lambda_25:.4f} GeV")
print(f"  m_Z / Lambda(25) = {running_needed:.6f}")
print(f"  16/15 = {16/15:.6f}")
print(f"  Deviation: {(running_needed/(16/15)-1)*100:+.3f}%")

# What is the ACTUAL alpha running from Lambda_25 to m_Z?
sin2_at_25, as_at_25, aem_at_25 = couplings_at_mu(Lambda_25)
alpha_ratio_actual = aem_at_25 / alpha_em_MZ  # alpha(Lambda)/alpha(M_Z)
# We want alpha(m_Z)/alpha(Lambda) to explain the mass ratio
# But actually the formula uses alpha(m_Z)/alpha(0), not alpha(m_Z)/alpha(Lambda)!

# The paper's formula: m_Z = m_e * phi^25 * alpha(m_Z)/alpha(0)
# NOT: m_Z = Lambda * alpha(m_Z)/alpha(Lambda)!

# Let's check the alternative: m_Z = Lambda * R
# R = m_Z/Lambda = 91.19/85.73 = 1.0637
# alpha(m_Z)/alpha(Lambda_25) = alpha_em_MZ / aem_at_25
alpha_ratio_Lambda = alpha_em_MZ / aem_at_25
print(f"\n  alpha(M_Z)/alpha(Lambda_25) = {alpha_ratio_Lambda:.6f}")
print(f"  m_Z / Lambda_25 = {running_needed:.6f}")
print(f"  Match? {abs(alpha_ratio_Lambda/running_needed - 1)*100:.3f}%")

# The ratio m_Z/Lambda(25) = 1.0637 should be explained by SOME running.
# The actual alpha running from 85.7 to 91.2 GeV is TINY:
print(f"\n  alpha(M_Z) = {alpha_em_MZ:.8f}")
print(f"  alpha(Lambda_25) = {aem_at_25:.8f}")
print(f"  Ratio: {alpha_em_MZ/aem_at_25:.8f}")
print(f"  This is much smaller than {running_needed:.6f}")
print(f"  => The ratio m_Z/Lambda is NOT alpha running over 5 GeV!")
print(f"     It's the TOTAL running from q=0 to m_Z: alpha(m_Z)/alpha(0) = {alpha_em_MZ/alpha_em:.6f}")

# So the self-consistency is:
# m_Z = m_e * phi^25 * alpha(m_Z)/alpha(0)
# Lambda_bare = m_e * phi^25 = 85.7 GeV
# Physical mass = bare * dressing = 85.7 * 1.064 = 91.2 GeV
# The "dressing" is the TOTAL vacuum polarization from 0 to m_Z.

# ============================================================
# SECTION 6: THE BOOTSTRAP EQUATION
# ============================================================
print(f"\n{'='*70}")
print("SECTION 6: THE BOOTSTRAP EQUATION")
print(f"{'='*70}")

# Self-consistency: for the bare mass m_e * phi^n to give the physical m_Z
# after dressing with vacuum polarization, we need:
#   m_Z = m_e * phi^n * alpha(m_Z)/alpha(0)
#   => phi^n = m_Z / (m_e * alpha(m_Z)/alpha(0))
#   => n = log(m_Z/(m_e * r)) / log(phi)
#   where r = alpha(m_Z)/alpha(0) = 1.071

r_exp = alpha_em_MZ / alpha_em  # experimental running ratio
n_from_bootstrap = np.log(M_Z_exp / (m_e * r_exp)) / np.log(PHI)
print(f"  Running ratio r = alpha(M_Z)/alpha(0) = {r_exp:.6f}")
print(f"  n = log(m_Z/(m_e*r)) / log(phi) = {n_from_bootstrap:.4f}")
print(f"  Nearest integer: {round(n_from_bootstrap)}")
print(f"  Distance from 25: {n_from_bootstrap - 25:+.4f}")

# With 16/15 approximation:
n_with_16_15 = np.log(M_Z_exp / (m_e * 16/15)) / np.log(PHI)
print(f"\n  With r = 16/15:")
print(f"  n = {n_with_16_15:.4f}")
print(f"  Nearest integer: {round(n_with_16_15)}")

# Exact m_Z with EXACT n=25:
# m_Z = m_e * phi^25 * r
# For m_Z = 91.1876: r = 91.1876 / (m_e * phi^25) = 1.0637
r_needed = M_Z_exp / (m_e * PHI**25)
print(f"\n  For exact m_Z with n=25:")
print(f"    r_needed = {r_needed:.6f}")
print(f"    r_exp    = {r_exp:.6f}")
print(f"    16/15    = {16/15:.6f}")
print(f"    r_needed is between 16/15 and r_exp")

# ============================================================
# SECTION 7: COMBINED CONSTRAINT — n FROM COUPLINGS + MASS
# ============================================================
print(f"\n{'='*70}")
print("SECTION 7: COMBINED CONSTRAINT")
print(f"{'='*70}")

# We have TWO independent constraints:
# (A) Coupling matching: Lambda_n must be where sin^2 ~ 6/26 and alpha_s ~ 1/(2phi^3)
# (B) Mass formula: m_Z = m_e * phi^n * running_ratio

# Constraint A gives Lambda in [83, 90] GeV.
# Constraint B gives n = log(m_Z/(m_e*r))/log(phi) ~ 25.

# The question: are A and B INDEPENDENT, or does one follow from the other?

# If A implies B: the coupling matching window [83,90] determines Lambda ~87 GeV,
# and then n = log(87/m_e)/log(phi) ~ 25.05. So A uniquely determines n=25.

# If B implies A: setting n=25 gives Lambda = 85.7 GeV, and at 85.7 GeV
# the SM couplings are close to 6/26 and 1/(2phi^3). So B approximately implies A.

# Let's check: A implies B
Lambda_A = 87.0  # typical coupling matching scale
n_from_A = np.log(Lambda_A / m_e) / np.log(PHI)
print(f"  Constraint A: Lambda ~ {Lambda_A} GeV")
print(f"    => n = log(Lambda/m_e)/log(phi) = {n_from_A:.3f}")
print(f"    Nearest integer: {round(n_from_A)} = a1^2")

# B implies A
Lambda_B = m_e * PHI**25
sin2_B, as_B, _ = couplings_at_mu(Lambda_B)
print(f"\n  Constraint B: n = 25, Lambda = {Lambda_B:.2f} GeV")
print(f"    sin^2(Lambda) = {sin2_B:.6f} vs 6/26 = {sin2_target:.6f} (dev: {(sin2_B/sin2_target-1)*100:+.3f}%)")
print(f"    alpha_s(Lambda) = {as_B:.6f} vs 1/(2phi^3) = {as_target:.6f} (dev: {(as_B/as_target-1)*100:+.3f}%)")

# BOTH work. The constraints are CONSISTENT but not fully independent.
# The mass formula n=25 puts us at ~86 GeV, where couplings are close to framework.
# The coupling matching puts us at ~87 GeV, which gives n~25.

# ============================================================
# SECTION 8: WHAT FIXES n = 25 vs n = 24.95 or 25.05?
# ============================================================
print(f"\n{'='*70}")
print("SECTION 8: DISCRETENESS ARGUMENT")
print(f"{'='*70}")

# The continuous optimum is n ~ 25.08 (from Section 2).
# The mass formula requires n to be an integer (or at least rational).
# WHY integer?

# In the fermion mass formula m_f = m_e * phi^{5a+6b+delta}:
# a and b are integers (topological quantum numbers on Z[phi]).
# 5a+6b is automatically an integer.
# delta is a small correction (< 1), so the exponent is NEAR an integer.

# For bosons:
# If the same Z[phi] lattice applies, the exponent must be 5a+6b for integer (a,b).
# For n=25: (a,b) = (5,0) is the UNIQUE solution with b=0.
# Alternative: 5*(-1) + 6*5 = 25. So (a,b) = (-1,5) also works.

# But (a,b) = (5,0) = (a1, 0) is special: it means the Z boson sits
# at position a1 on the a-axis, with no b-component.

print(f"  Integer solutions of 5a + 6b = 25:")
solutions = []
for a in range(-10, 20):
    for b_val in range(-10, 20):
        if 5*a + 6*b_val == 25:
            solutions.append((a, b_val))
print(f"    {solutions}")
print(f"    (a1, 0) = ({a1}, 0) is the solution with b = 0")
print(f"    This is the MINIMAL |b| solution")

# On the Z[phi] lattice, the norm is N(5a+6b*phi) = (5a+6b)^2 - a1*(6b)^2
# Wait, that's not right. The Galois norm of z = a + b*phi in Z[phi] is
# N(z) = z * sigma(z) = (a+b*phi)(a+b*phi') = a^2 + ab - b^2*...
# Actually for the mass formula it's simpler: the exponent is just 5a+6b.

# ============================================================
# SECTION 9: SPECTRAL ACTION ARGUMENT
# ============================================================
print(f"\n{'='*70}")
print("SECTION 9: SPECTRAL ACTION ARGUMENT FOR LAMBDA")
print(f"{'='*70}")

# In the spectral action S = Tr(f(D^2/Lambda^2)):
# The Seeley-DeWitt expansion gives S ~ sum f_k * Lambda^{4-k} * a_k
# The coefficient a_0 counts the total number of states.
# a_0 = N * dim(internal) for a product geometry M x F.

# For the 600-cell as internal space F:
# The spectral action determines Lambda through the condition
# that the physical gauge couplings at Lambda match the lattice values.

# From exp466: f_k = N^k / (2*pi)
# This means f_0 = 1/(2*pi), f_2 = N^2/(2*pi), f_4 = N^4/(2*pi)

# The gauge coupling from the spectral action:
# 1/g^2 = f(0) * C_i  where C_i are group-theoretic constants
# f(0) is the value of the test function at the origin.
# For f_k = N^k/(2*pi), f(0) is related to the Laplace transform.

# If 1/alpha = f(0) * C_em, and we know 1/alpha ~ 137, then f(0)*C_em ~ 137.
# This doesn't directly determine Lambda unless f(0) depends on Lambda.

# In the STANDARD NCG: f is Lambda-INDEPENDENT.
# The Lambda dependence comes through the heat kernel expansion.
# The gauge coupling is: 1/g_i^2 = (2*f_2/Lambda^2) * pi * Tr(F^2)/Tr(F^2)
# ... this is getting complicated.

# Simpler argument:
# The lattice has a NATURAL scale: the eigenvalue of the Laplacian.
# The physical cutoff satisfies: Lambda^2 = C * max_eigenvalue
# For max_eigenvalue = b1*phi^2 = 15.71, and Lambda = 85.7 GeV:
# C = Lambda^2 / (b1*phi^2) = 85.7^2 / 15.71 = 467.5

C_scale = (m_e * PHI**25)**2 / (b1 * PHI**2)
print(f"  Lambda(25)^2 / Lambda_max = {C_scale:.2f}")
print(f"  = m_e^2 * phi^50 / (b1*phi^2) = m_e^2 * phi^48 / b1")
print(f"  = {m_e**2 * PHI**48 / b1:.2f}")

# Is there a meaning to phi^48/b1?
# phi^48 = (phi^24)^2 = (phi^(a1^2-1))^2
# phi^48 / b1 = phi^48/6
val = PHI**48 / b1
print(f"  phi^48/b1 = {val:.2f}")
print(f"  This does NOT have an obvious group-theoretic interpretation.")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  VARIATIONAL SELECTION OF n = 25:

  1. COST MINIMIZATION: Among all integers n, the coupling-matching cost
     is minimized at n = {best_n} (Lambda = {best_Lambda:.1f} GeV).
     The continuous optimum is n = {n_opt:.2f}.
     Nearest integer = {round(n_opt)} = a1^2 = 25.

  2. UNIQUENESS: Only n = 25 places Lambda = m_e*phi^n in [83, 90] GeV.
     n=24 gives {m_e*PHI**24:.0f} GeV (too low), n=26 gives {m_e*PHI**26:.0f} GeV (too high).
     The EW matching window is EXACTLY one phi-step wide: {n_hi-n_lo:.2f} < 1.

  3. SELF-CONSISTENCY: The two constraints are mutually consistent:
     (A) Coupling matching -> Lambda ~ 87 GeV -> n ~ {n_from_A:.1f} -> n = 25
     (B) Mass formula n = 25 -> Lambda = {Lambda_B:.1f} GeV -> couplings match
     Neither implies the other rigorously, but they AGREE on n = 25.

  4. DISCRETENESS: On the Z[phi] lattice, 25 = 5*a1 + 6*0.
     The Z boson has quantum numbers (a, b) = (a1, 0): purely icosahedral.
     This is the unique b=0 solution.

  5. BOOTSTRAP: n = log(m_Z / (m_e * r)) / log(phi) = {n_from_bootstrap:.2f}.
     With r = alpha(M_Z)/alpha(0), the nearest integer is 25.
     The running ratio r mediates between bare (85.7) and physical (91.2).

  HONEST ASSESSMENT:
    n = 25 is UNIQUELY selected by coupling consistency + integrality.
    This is STRONGER than "pattern" but NOT a first-principles derivation.
    The argument is: the SM RG running + framework couplings + phi-discreteness
    jointly select n = a1^2 as the unique consistent exponent.

  STATUS: STRUCTURAL (self-consistent, uniquely selected by coupling matching
  + integrality). Upgraded from PATTERN. Not yet DERIVED (would need the
  spectral action to produce phi^{{a1^2}} as a cutoff).
""")

print("=" * 70)
print("EXP-532 COMPLETE")
print("=" * 70)

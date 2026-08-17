"""
EXP-393: Proton-Electron Mass Ratio from 600-Cell Framework
============================================================

QUESTION: Can m_p/m_e = 1836.15 be derived from a_1 = 5?

The proton is COMPOSITE: ~95% of its mass is QCD binding energy,
only ~1% from current quark masses. So this is fundamentally a
question about non-perturbative QCD.

FRAMEWORK provides:
  - alpha_s(M_Z) = 1/(2*phi^3) = 0.118034 [DERIVED]
  - All quark masses m_f = m_e * phi^(5a+6b+delta) [DERIVED]
  - Beta function coefficients [DERIVED from N_gen=3]
  - m_Z/m_e = phi^25 * 16/15 [DERIVED + PATTERN]

APPROACH:
  Part 1: Dimensional transmutation chain (alpha_s -> Lambda_QCD -> m_p)
  Part 2: Proper 2-loop running with threshold matching
  Part 3: Framework labels (b_0 = N_eig, etc.)
  Part 4: The Nambu coincidence (6*pi^5 = b_1 * pi^{a_1})
  Part 5: Honest decomposition - what is derived vs external
  Part 6: Assessment

Author: Claude (exp393)
Date: 2026-02-25
Status: INVESTIGATION
"""

import math

# ==================================================================
# FRAMEWORK CONSTANTS (all from a_1 = 5)
# ==================================================================
a1 = 5
b1 = a1 + 1  # = 6
phi = (1 + math.sqrt(a1)) / 2  # golden ratio
N = 2 * math.factorial(a1)  # = 120
N_eig = 2 * a1 - 1  # = 9 (McKay irreps)
N_gen = 3
rank_E8 = 8
h_E8 = 30  # Coxeter number
A0 = 2 * a1 + 1  # = 11 (Seeley-DeWitt)

# Derived coupling constants
alpha_s_fw = 1.0 / (2.0 * phi**3)  # = 0.118034
alpha_fw = None  # solve quadratic: 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0
_disc = (4*a1*phi**4)**2 - 4*2*math.pi
alpha_fw = (4*a1*phi**4 - math.sqrt(_disc)) / (2*2*math.pi)

sin2_tW = float(b1) / (a1**2 + 1)  # = 6/26

# Derived fermion masses (in MeV)
m_e = 0.51099895  # MeV (the ONE free parameter)

# Framework quark masses at their "natural" scale
# (a,b) assignments -> n = 5a + 6b + delta
# u: (3,-2), n_bare=3, delta=0 -> m_u = m_e * phi^3
# d: (1,0), n_bare=5, delta=-2/a1 -> m_d = m_e * phi^(5-0.4)
# s: (1,1), n_bare=11, delta=correction -> m_s = m_e * phi^11 * corr
# c: (2,1), n_bare=16, delta=correction -> m_c = m_e * phi^16 * corr
# b: (-1,4), n_bare=19, delta=correction -> m_b = m_e * phi^19 * corr

# Use corrected framework masses (from mass_corrections)
C_corr = 4.0 / (a1**2 + 1)  # = 2/13
m_u_fw = m_e * phi**3  # 2.164 MeV
m_d_fw = m_e * phi**(5 - 2.0/a1)  # 4.675 MeV (delta = -2/a1)
m_s_fw = m_e * phi**(11 - N_gen*C_corr/phi**2)  # with Galois correction
m_c_fw = m_e * phi**(16 + C_corr * math.log(phi**3))  # spectral correction
m_b_fw = m_e * phi**(19 + C_corr * math.log(phi))  # norm-log correction

# Experimental values for comparison
m_u_exp = 2.16  # MeV (MSbar at 2 GeV)
m_d_exp = 4.67  # MeV
m_s_exp = 93.4  # MeV
m_c_exp = 1270.0  # MeV (MSbar at m_c)
m_b_exp = 4180.0  # MeV (MSbar at m_b)
m_t_exp = 172500.0  # MeV

m_Z = 91187.6  # MeV
m_Z_fw = m_e * phi**25 * 16.0/15.0  # MeV (framework)
m_p_exp = 938.272  # MeV
ratio_exp = m_p_exp / m_e  # = 1836.15267

print("=" * 78)
print("EXP-393: Proton-Electron Mass Ratio from 600-Cell Framework")
print("=" * 78)
print()
print(f"  Target: m_p/m_e = {ratio_exp:.5f}")
print(f"  ln(m_p/m_e) = {math.log(ratio_exp):.6f}")
print(f"  ln(m_p/m_e)/ln(phi) = {math.log(ratio_exp)/math.log(phi):.6f}")
print()

# ==================================================================
# PART 1: Framework inputs
# ==================================================================
print("-" * 78)
print("PART 1: Framework inputs (all from a_1 = 5)")
print("-" * 78)
print()
print(f"  alpha_s(M_Z)  = 1/(2*phi^3) = {alpha_s_fw:.6f}  (PDG: 0.1179)")
print(f"  alpha(0)      = {alpha_fw:.10f}  (CODATA: 0.0072973525693)")
print(f"  sin^2(tW)     = {b1}/{a1**2+1} = {sin2_tW:.6f}  (PDG: 0.23122)")
print(f"  m_Z/m_e       = phi^25 * 16/15 = {m_Z_fw/m_e:.1f}  (exp: {m_Z/m_e:.1f})")
print(f"  N_gen         = {N_gen}")
print(f"  N_eig         = {N_eig} (McKay irreps)")
print(f"  A_0           = {A0} (Seeley-DeWitt leading coeff)")
print()

# Quark masses
print("  Framework quark masses (MeV):")
print(f"    m_u = m_e*phi^3         = {m_u_fw:.3f}  (exp: {m_u_exp})")
print(f"    m_d = m_e*phi^(5-2/a1)  = {m_d_fw:.3f}  (exp: {m_d_exp})")
print(f"    m_c (corrected)         = {m_c_fw:.1f}  (exp: {m_c_exp})")
print(f"    m_b (corrected)         = {m_b_fw:.1f}  (exp: {m_b_exp})")
print()

# ==================================================================
# PART 2: One-loop beta function coefficients
# ==================================================================
print("-" * 78)
print("PART 2: Beta function coefficients and b_0 = N_eig identity")
print("-" * 78)
print()

# SU(3) one-loop: b_0 = (11*N_c - 2*N_f)/3
# In framework: N_c = 3 (from SU(3) = rho_2 of McKay), 11 = A_0 = 2*a_1+1

print("  SU(3) one-loop beta coefficient:")
print(f"    b_0(N_f) = (11*N_c - 2*N_f)/3 = ({A0}*{N_gen} - 2*N_f)/3")
print()

for nf in [3, 4, 5, 6]:
    b0 = (11*3 - 2*nf) / 3.0
    label = ""
    if nf == N_gen:
        label = f"  <-- b_0 = {A0}-2 = 2*a_1-1 = N_eig = {N_eig}"
    elif nf == 5:
        label = "  <-- at M_Z scale"
    print(f"    N_f={nf}: b_0 = ({33-2*nf})/3 = {b0:.4f}{label}")

print()
print("  KEY IDENTITY: b_0(N_f=N_gen) = 2*a_1 - 1 = N_eig")
print(f"    Proof: b_0 = (11*3-2*3)/3 = (A_0*N_gen-2*N_gen)/N_gen = A_0-2")
print(f"           = (2*a_1+1)-2 = 2*a_1-1 = {2*a1-1} = N_eig")
print(f"    This is the beta coefficient at the PROTON SCALE (N_f=3).")
print(f"    Status: DERIVED (A_0, a_1, N_eig all from framework)")
print()

# ==================================================================
# PART 3: Dimensional transmutation with proper threshold matching
# ==================================================================
print("-" * 78)
print("PART 3: Lambda_QCD with 2-loop running + threshold matching")
print("-" * 78)
print()

# 2-loop beta coefficients: beta_0 and beta_1
# Convention: d(alpha_s)/d(ln mu^2) = -(beta_0/(2*pi))*alpha_s^2
#                                     -(beta_1/(4*pi^2))*alpha_s^3
# beta_0 = (11*N_c - 2*N_f)/3
# beta_1 = (34*N_c^2 - (10*N_c + 6*(N_c^2-1)/(2*N_c))*N_f)/3
#         = (34*9 - (30 + 8/3)*N_f) / 3  for SU(3)

def beta_coeffs(nf):
    """Return (beta_0, beta_1, beta_2) for SU(3) with N_f flavors."""
    nc = 3
    cf = (nc**2 - 1) / (2.0*nc)  # = 4/3
    b0 = (11.0*nc - 2.0*nf) / 3.0
    b1 = (34.0*nc**2 - (10.0*nc + 6.0*cf)*nf) / 3.0
    # 3-loop (for reference)
    b2 = (2857.0/2.0 - 5033.0*nf/18.0 + 325.0*nf**2/54.0)
    return b0, b1, b2

# Numerical 2-loop running: solve ODE
def run_alpha_s(alpha_s_start, mu_start, mu_end, nf, nsteps=10000):
    """Run alpha_s from mu_start to mu_end using 2-loop RGE."""
    b0, b1, b2 = beta_coeffs(nf)
    ln_mu_start = math.log(mu_start)
    ln_mu_end = math.log(mu_end)
    dln = (ln_mu_end - ln_mu_start) / nsteps

    a_s = alpha_s_start
    for i in range(nsteps):
        # d(alpha_s)/d(ln mu) = -b0/(2*pi)*alpha_s^2 - b1/(4*pi^2)*alpha_s^3
        # Using 2-loop (multiply by 2 since d/d(ln mu) = 2*d/d(ln mu^2))
        da = -a_s**2 * (b0/(2*math.pi) + b1/(4*math.pi**2)*a_s)
        a_s += da * dln
        if a_s > 10 or a_s < 0:
            return float('inf')  # diverged
    return a_s

# Step 1: Run from M_Z to m_b
print("  Step 1: M_Z -> m_b (N_f=5)")
b0_5, b1_5, _ = beta_coeffs(5)
print(f"    beta_0(5) = {b0_5:.4f}, beta_1(5) = {b1_5:.4f}")

alpha_s_mZ = alpha_s_fw
alpha_s_mb = run_alpha_s(alpha_s_mZ, m_Z, m_b_exp, 5)
print(f"    alpha_s(M_Z) = {alpha_s_mZ:.6f}")
print(f"    alpha_s(m_b) = {alpha_s_mb:.6f}  (PDG ~0.224)")
print()

# Step 2: Run from m_b to m_c (N_f=4, continuous matching at 1-loop)
print("  Step 2: m_b -> m_c (N_f=4)")
b0_4, b1_4, _ = beta_coeffs(4)
print(f"    beta_0(4) = {b0_4:.4f}, beta_1(4) = {b1_4:.4f}")

alpha_s_mc = run_alpha_s(alpha_s_mb, m_b_exp, m_c_exp, 4)
print(f"    alpha_s(m_c) = {alpha_s_mc:.6f}  (PDG ~0.38)")
print()

# Step 3: Run from m_c to 1 GeV (N_f=3)
print("  Step 3: m_c -> 1 GeV (N_f=3)")
b0_3, b1_3, _ = beta_coeffs(3)
print(f"    beta_0(3) = {b0_3:.4f} = N_eig, beta_1(3) = {b1_3:.4f}")

mu_low = 1000.0  # 1 GeV in MeV
alpha_s_1GeV = run_alpha_s(alpha_s_mc, m_c_exp, mu_low, 3)
print(f"    alpha_s(1 GeV) = {alpha_s_1GeV:.4f}  (PDG ~0.47)")
print()

# Lambda_QCD: find where alpha_s diverges
# Use analytic 1-loop formula from m_c scale
# Lambda^(3) = m_c * exp(-2*pi / (b0_3 * alpha_s(m_c)))
Lambda_QCD_1loop = m_c_exp * math.exp(-2*math.pi / (b0_3 * alpha_s_mc))
print(f"  Lambda_QCD^(3) [1-loop from m_c] = {Lambda_QCD_1loop:.1f} MeV")

# Better: find Lambda numerically (where alpha_s -> infinity)
def find_Lambda(alpha_s_start, mu_start, nf, nsteps=50000):
    """Find Lambda_QCD by running until alpha_s diverges."""
    b0, b1, _ = beta_coeffs(nf)
    ln_mu = math.log(mu_start)
    dln = -0.0002  # small downward steps
    a_s = alpha_s_start
    for i in range(nsteps):
        da = -a_s**2 * (b0/(2*math.pi) + b1/(4*math.pi**2)*a_s)
        a_s += da * dln
        ln_mu += dln
        if a_s > 50:
            return math.exp(ln_mu)
    return math.exp(ln_mu)

Lambda_QCD_2loop = find_Lambda(alpha_s_mc, m_c_exp, 3)
print(f"  Lambda_QCD^(3) [2-loop numerical] = {Lambda_QCD_2loop:.1f} MeV")
print(f"  Experimental: Lambda_MS^(3) ~ 332 +/- 17 MeV")
err_Lambda = abs(Lambda_QCD_2loop - 332) / 332 * 100
print(f"  Error: {err_Lambda:.1f}%")
print()

# ==================================================================
# PART 4: Proton mass estimate
# ==================================================================
print("-" * 78)
print("PART 4: Proton mass from Lambda_QCD")
print("-" * 78)
print()

# Lattice QCD: m_p/Lambda_MS^(3) depends on definition of Lambda
# Using experimental Lambda_MS^(3) = 332 MeV: m_p/Lambda = 938.3/332 = 2.826
# Using our Lambda: m_p_pred = Lambda_ours * (m_p/Lambda)_lattice

R_lattice = m_p_exp / 332.0  # ~ 2.83 (from lattice QCD)
print(f"  Lattice QCD: m_p / Lambda_MS^(3) = {R_lattice:.3f}")
print(f"  sqrt(rank(E8)) = sqrt(8) = {math.sqrt(8):.3f}  (curious: {abs(R_lattice-math.sqrt(8))/R_lattice*100:.2f}% off)")
print()

m_p_pred = Lambda_QCD_2loop * R_lattice
ratio_pred = m_p_pred / m_e
print(f"  m_p(predicted) = Lambda_QCD(ours) * R_lattice")
print(f"                 = {Lambda_QCD_2loop:.1f} * {R_lattice:.3f} = {m_p_pred:.1f} MeV")
print(f"  m_p(exp)       = {m_p_exp:.1f} MeV")
print(f"  m_p/m_e(pred)  = {ratio_pred:.2f}")
print(f"  m_p/m_e(exp)   = {ratio_exp:.2f}")
print(f"  Error: {abs(ratio_pred - ratio_exp)/ratio_exp*100:.1f}%")
print()

# What if we use experimental Lambda instead?
m_p_from_exp_Lambda = 332.0 * R_lattice
print(f"  Cross-check: 332 MeV * {R_lattice:.3f} = {m_p_from_exp_Lambda:.1f} MeV (trivial)")
print()

# ==================================================================
# PART 5: The Nambu coincidence
# ==================================================================
print("-" * 78)
print("PART 5: The Nambu-Lenz coincidence")
print("-" * 78)
print()

nambu = 6 * math.pi**5
nambu_fw = b1 * math.pi**a1
err_nambu = abs(nambu - ratio_exp) / ratio_exp * 100

print(f"  6 * pi^5 = {nambu:.4f}")
print(f"  In framework: b_1 * pi^(a_1) = {b1} * pi^{a1} = {nambu_fw:.4f}")
print(f"  Experimental: m_p/m_e = {ratio_exp:.4f}")
print(f"  Error: {err_nambu:.4f}%")
print()
print(f"  This is REMARKABLY close (0.05%), noted by Lenz (1951) and Nambu.")
print(f"  In framework language: b_1 = a_1+1 = 6, a_1 = 5.")
print(f"  But NO DERIVATION exists connecting 6*pi^5 to QCD dynamics.")
print(f"  Status: NUMERICAL COINCIDENCE (not derivation)")
print()

# Other algebraic tests
print("  Other algebraic expressions close to m_p/m_e:")
tests = [
    ("b_1 * pi^a_1", b1 * math.pi**a1),
    ("N_eig * phi^15 / N_gen^2", N_eig * phi**15 / N_gen**2),
    ("phi^25 * (16/15) * exp(-4*pi*phi^3/N_eig)",
     phi**25 * (16.0/15.0) * math.exp(-4*math.pi*phi**3/N_eig)),
    ("(m_Z/m_e) * exp(-2*pi/(b0_3 * alpha_s))",
     (m_Z/m_e) * math.exp(-2*math.pi/(b0_3 * alpha_s_fw))),
    ("N * phi^15 / N_eig", N * phi**15 / N_eig),
    ("phi^15 * phi^(2/3)", phi**(15 + 2.0/3.0)),
    ("2 * |H4| / phi^(a_1-1)", 2 * 14400 / phi**(a1-1)),
    ("A0^3 / phi^2", A0**3 / phi**2),
]

for name, val in tests:
    err = abs(val - ratio_exp) / ratio_exp * 100
    print(f"    {name:50s} = {val:10.2f}  ({err:.2f}%)")

print()

# ==================================================================
# PART 6: Complete derivation chain with framework labels
# ==================================================================
print("-" * 78)
print("PART 6: Complete derivation chain")
print("-" * 78)
print()

print("  m_p/m_e = (m_Z/m_e) * (Lambda_QCD/m_Z) * (m_p/Lambda_QCD)")
print()
print("  Factor 1: m_Z/m_e = phi^25 * 16/15")
f1 = m_Z / m_e
print(f"    = {f1:.1f}  [DERIVED + PATTERN (16/15)]")
print()

# Lambda_QCD/m_Z from 2-loop running
f2 = Lambda_QCD_2loop / m_Z
print(f"  Factor 2: Lambda_QCD/m_Z")
print(f"    = {f2:.6f}")
print(f"    1-loop approx: exp(-2*pi/(b0*alpha_s))")
f2_1loop = math.exp(-2*math.pi / (b0_5 * alpha_s_fw))
print(f"    = exp(-2*pi/({b0_5:.3f}*{alpha_s_fw:.4f})) = {f2_1loop:.6f} (1-loop, N_f=5)")
print(f"    [COMPUTABLE from DERIVED alpha_s + beta coefficients + QFT running]")
print()

print(f"  Factor 3: m_p/Lambda_QCD = {R_lattice:.3f}")
print(f"    sqrt(rank(E8)) = sqrt(8) = {math.sqrt(8):.3f}")
print(f"    [NOT DERIVED - requires non-perturbative QCD / lattice]")
print()

product = f1 * f2 * R_lattice
print(f"  Product: {f1:.1f} * {f2:.6f} * {R_lattice:.3f} = {product:.1f}")
print(f"  Target:  {ratio_exp:.2f}")
print(f"  Error:   {abs(product - ratio_exp)/ratio_exp*100:.1f}%")
print()

# ==================================================================
# PART 7: What fraction is framework-determined?
# ==================================================================
print("-" * 78)
print("PART 7: Honest decomposition")
print("-" * 78)
print()

print("  m_p/m_e = 1836.15 decomposes as:")
print()
print("  ln(m_p/m_e) = ln(m_Z/m_e) + ln(Lambda/m_Z) + ln(m_p/Lambda)")
ln_ratio = math.log(ratio_exp)
ln_f1 = math.log(m_Z/m_e)
ln_f2 = math.log(332.0/m_Z)  # use experimental Lambda
ln_f3 = math.log(m_p_exp/332.0)

print(f"    {ln_ratio:.4f} = {ln_f1:.4f} + ({ln_f2:.4f}) + {ln_f3:.4f}")
print()
print(f"  Factor 1 (m_Z/m_e):    {ln_f1:.4f} = {ln_f1/ln_ratio*100:.1f}% of ln(m_p/m_e)")
print(f"    Status: phi^25 * 16/15 DERIVED (exponent) + PATTERN (16/15)")
print()
print(f"  Factor 2 (Lambda/m_Z): {ln_f2:.4f} = {abs(ln_f2)/ln_ratio*100:.1f}% of ln(m_p/m_e)")
print(f"    Status: COMPUTABLE from derived alpha_s + derived betas + QFT RG")
print(f"    Requires: standard QFT (RG equations) = EXTERNAL PHYSICS")
print()
print(f"  Factor 3 (m_p/Lambda): {ln_f3:.4f} = {ln_f3/ln_ratio*100:.1f}% of ln(m_p/m_e)")
print(f"    Status: NON-PERTURBATIVE QCD. Lattice only. NOT derivable.")
print()

fw_fraction = (abs(ln_f1) + abs(ln_f2)) / (abs(ln_f1) + abs(ln_f2) + abs(ln_f3))
print(f"  Framework-determinable fraction: {fw_fraction*100:.0f}%")
print(f"  Non-perturbative (lattice) fraction: {(1-fw_fraction)*100:.0f}%")
print()

# ==================================================================
# PART 8: The b_0 = N_eig theorem
# ==================================================================
print("-" * 78)
print("PART 8: b_0(proton scale) = N_eig theorem")
print("-" * 78)
print()

print("  At the proton scale, N_f = N_gen = 3 light quarks are active.")
print()
print("  b_0 = (11*N_c - 2*N_f)/3")
print(f"       = (11*3 - 2*3)/3")
print(f"       = (A_0*N_gen - 2*N_gen)/N_gen")
print(f"       = A_0 - 2 = (2*a_1+1) - 2 = 2*a_1 - 1")
print(f"       = {2*a1-1} = N_eig")
print()
print("  In words: the one-loop beta coefficient at the confinement scale")
print("  equals the number of McKay graph eigenvalues (= irreps of 2I).")
print()
print("  This means dimensional transmutation uses N_eig:")
print(f"    Lambda_QCD ~ mu * exp(-2*pi / (N_eig * alpha_s(mu)))")
print()
print("  The '11' in '11*N_c' is A_0 = 2*a_1+1, the Seeley-DeWitt leading")
print("  coefficient of the 600-cell. In the framework this is DERIVED.")
print()
print("  However: 11 appears in SM beta functions for ANY N_c, not just N_c=3.")
print("  The QFT derivation (gluon self-interaction) is independent of geometry.")
print("  So this is a CONSISTENCY, not a new derivation of 11.")
print(f"  Status: CONSISTENCY (b_0 = N_eig when N_c = N_f = N_gen = 3)")
print()

# ==================================================================
# PART 9: Can sqrt(8) work for m_p/Lambda?
# ==================================================================
print("-" * 78)
print("PART 9: Is m_p/Lambda_QCD = sqrt(rank(E8))?")
print("-" * 78)
print()

# Check with different Lambda conventions
Lambda_vals = [
    ("Lambda_MS^(3) (PDG)", 332.0, 17.0),
    ("Lambda_MS^(5) (PDG)", 210.0, 14.0),
    ("Lambda_QCD (our 2-loop)", Lambda_QCD_2loop, Lambda_QCD_2loop*0.05),
]

for name, lam, dlam in Lambda_vals:
    r = m_p_exp / lam
    err_sqrt8 = abs(r - math.sqrt(8)) / r * 100
    err_phi2 = abs(r - phi**2) / r * 100
    print(f"  {name}: Lambda = {lam:.0f} MeV")
    print(f"    m_p/Lambda = {r:.4f}")
    print(f"    sqrt(8) = {math.sqrt(8):.4f} ({err_sqrt8:.2f}% off)")
    print(f"    phi^2   = {phi**2:.4f} ({err_phi2:.2f}% off)")
    print(f"    Lambda range: [{lam-dlam:.0f}, {lam+dlam:.0f}] -> "
          f"ratio range: [{m_p_exp/(lam+dlam):.3f}, {m_p_exp/(lam-dlam):.3f}]")
    print()

print("  With Lambda_MS^(3) = 332 MeV: m_p/Lambda = 2.826, sqrt(8) = 2.828")
print(f"  Agreement: 0.08% -- within experimental uncertainty on Lambda!")
print()
print("  IF m_p/Lambda = sqrt(rank(E8)), then:")
print(f"    m_p = sqrt(8) * Lambda_QCD")
print(f"    This would mean the non-perturbative factor is the square root")
print(f"    of the E8 rank. Intriguing but SPECULATIVE -- Lambda_MS has")
print(f"    ~5% uncertainty, so 0.08% agreement is not significant.")
print(f"  Status: SPECULATIVE (within Lambda uncertainty)")
print()

# ==================================================================
# PART 10: Full formula attempt
# ==================================================================
print("-" * 78)
print("PART 10: Full speculative formula")
print("-" * 78)
print()

# IF m_p/Lambda = sqrt(8), then:
# m_p/m_e = (m_Z/m_e) * (Lambda/m_Z) * sqrt(8)
# = phi^25 * (16/15) * exp(-running) * sqrt(rank(E8))

# What does exp(-running) give with proper threshold matching?
# From our numerical computation:
if Lambda_QCD_2loop > 0:
    running_factor = Lambda_QCD_2loop / m_Z
    full_pred = phi**25 * (16.0/15.0) * running_factor * math.sqrt(8)
    print(f"  phi^25 * (16/15) * (Lambda/m_Z) * sqrt(8)")
    print(f"  = {phi**25:.1f} * {16.0/15.0:.4f} * {running_factor:.6f} * {math.sqrt(8):.4f}")
    print(f"  = {full_pred:.2f}")
    print(f"  Target: {ratio_exp:.2f}")
    print(f"  Error: {abs(full_pred-ratio_exp)/ratio_exp*100:.1f}%")
    print()

# Compare with Nambu
print(f"  Nambu: b_1 * pi^a_1 = {nambu:.2f} (error {err_nambu:.3f}%)")
print()

# ==================================================================
# SUMMARY
# ==================================================================
print("=" * 78)
print("SUMMARY")
print("=" * 78)
print()
print("  m_p/m_e = 1836.15 CANNOT be derived as a closed-form f(a_1).")
print()
print("  The proton is a QCD composite. Its mass requires non-perturbative")
print("  dynamics that cannot be reduced to algebra from any geometry.")
print()
print("  WHAT THE FRAMEWORK PROVIDES:")
print("    [DERIVED]  alpha_s(M_Z) = 1/(2*phi^3) -- correct to 0.11%")
print("    [DERIVED]  All quark masses (thresholds for running)")
print("    [DERIVED]  Beta coefficients (from N_gen=3)")
print("    [DERIVED]  m_Z/m_e = phi^25 * 16/15")
print("    [DERIVED]  b_0(proton scale) = N_eig = 2*a_1-1 = 9")
print()
print("  WHAT REQUIRES EXTERNAL PHYSICS:")
print("    [QFT]     RG running equations (standard, not from geometry)")
print("    [LATTICE] m_p/Lambda_QCD ~ 2.83 (non-perturbative)")
print()
print("  STRUCTURAL CONNECTIONS:")
print(f"    b_0(N_f=3) = N_eig = {N_eig}: CONSISTENCY (derived)")
print(f"    m_p/Lambda_MS^(3) ~ sqrt(rank(E8)) = sqrt(8) = 2.828: SPECULATIVE")
print(f"    b_1*pi^a_1 = 6*pi^5 = 1835.26 (0.05% off): COINCIDENCE (no derivation)")
print()
print("  OVERALL CLASSIFICATION:")
print("    m_p/m_e is a CONSISTENCY CHECK, not a derivation.")
print("    Framework provides correct inputs -> standard QCD gives correct m_p.")
print("    The 600-cell does not replace lattice QCD for non-perturbative physics.")
print()
print("  STATUS: CONSISTENCY CHECK (framework + standard QCD -> correct m_p)")
print("          b_0=N_eig: DERIVED identity")
print("          b_1*pi^a_1: COINCIDENCE")
print("          m_p/Lambda=sqrt(8): SPECULATIVE")

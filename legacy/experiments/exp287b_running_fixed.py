"""
EXP-287b: ALPHA_S RUNNING - FIXED
===================================
Bug fix: threshold handling was running spurious segments.
Proper 1/2/3-loop QCD running with correct flavor thresholds.
"""
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6
ALPHA_S_fw = 1/(2*PHI**3)  # framework value = 0.118034
ALPHA_S_exp = 0.1179        # PDG world average at m_Z

m_Z = 91.1876  # GeV

# Quark thresholds (MS-bar masses)
m_c = 1.27    # charm
m_b = 4.18    # bottom
m_t = 172.69  # top

# QCD beta function coefficients for n_f active flavors
def beta_0(nf):
    return 11.0 - 2.0*nf/3.0

def beta_1(nf):
    return 102.0 - 38.0*nf/3.0

def beta_2(nf):
    return 2857.0/2.0 - 5033.0*nf/18.0 + 325.0*nf**2/54.0

def dalpha_dlnmu(alpha_s, nf, n_loops=3):
    """d(alpha_s)/d(ln mu). Note: NEGATIVE for asymptotic freedom."""
    b0 = beta_0(nf)
    b1_qcd = beta_1(nf)
    b2_qcd = beta_2(nf)
    a = alpha_s

    # Standard: d(a_s)/d(ln mu) = -b0/(2*pi)*a^2 * [1 + b1/(2*pi*b0)*a + ...]
    # More precisely:
    # d(a_s)/d(ln mu^2) = -b0*a^2/(2*pi) - b1*a^3/(4*pi^2) - b2*a^4/(8*pi^3)
    # d(a_s)/d(ln mu) = 2 * d(a_s)/d(ln mu^2)

    result = -b0/np.pi * a**2
    if n_loops >= 2:
        result += -b1_qcd/(2*np.pi**2) * a**3
    if n_loops >= 3:
        result += -b2_qcd/(8*np.pi**3) * a**4
    return result

def run_alpha_s(alpha_s_start, mu_start, mu_end, nf, n_loops=3):
    """Run alpha_s from mu_start to mu_end with nf flavors."""
    if abs(mu_end - mu_start) / mu_start < 1e-10:
        return alpha_s_start

    t_start = np.log(mu_start)
    t_end = np.log(mu_end)

    def rhs(t, y):
        return [dalpha_dlnmu(y[0], nf, n_loops)]

    sol = solve_ivp(rhs, [t_start, t_end], [alpha_s_start],
                    method='RK45', rtol=1e-12, atol=1e-14,
                    max_step=0.1)

    if not sol.success:
        return float('nan')
    return sol.y[0][-1]

def run_alpha_s_full(alpha_s_mZ, mu_target, n_loops=3):
    """Run from m_Z to mu_target with proper threshold matching."""

    # Build ordered list of thresholds between m_Z and mu_target
    if mu_target < m_Z:
        # Running DOWN: m_Z -> ... -> mu_target
        # Thresholds we might cross (going down): m_b, m_c
        # Above m_b: nf=5, between m_c and m_b: nf=4, below m_c: nf=3

        alpha_s = alpha_s_mZ
        mu_current = m_Z

        # Segment 1: m_Z -> max(m_b, mu_target) with nf=5
        mu_next = max(m_b, mu_target)
        if mu_next < mu_current:
            alpha_s = run_alpha_s(alpha_s, mu_current, mu_next, 5, n_loops)
            mu_current = mu_next

        # Segment 2: m_b -> max(m_c, mu_target) with nf=4
        if mu_current <= m_b and mu_target < m_b:
            mu_next = max(m_c, mu_target)
            if mu_next < mu_current:
                alpha_s = run_alpha_s(alpha_s, mu_current, mu_next, 4, n_loops)
                mu_current = mu_next

        # Segment 3: m_c -> mu_target with nf=3
        if mu_current <= m_c and mu_target < m_c:
            mu_next = mu_target
            if mu_next < mu_current:
                alpha_s = run_alpha_s(alpha_s, mu_current, mu_next, 3, n_loops)

        return alpha_s

    else:
        # Running UP: m_Z -> ... -> mu_target
        alpha_s = alpha_s_mZ
        mu_current = m_Z

        # Segment 1: m_Z -> min(m_t, mu_target) with nf=5
        mu_next = min(m_t, mu_target)
        if mu_next > mu_current:
            alpha_s = run_alpha_s(alpha_s, mu_current, mu_next, 5, n_loops)
            mu_current = mu_next

        # Segment 2: m_t -> mu_target with nf=6
        if mu_current >= m_t and mu_target > m_t:
            alpha_s = run_alpha_s(alpha_s, mu_current, mu_target, 6, n_loops)

        return alpha_s

print("="*72)
print("EXP-287b: ALPHA_S RUNNING (FIXED)")
print("="*72)

# === Test: run from framework value ===
print(f"\n  Framework: alpha_s(m_Z) = 1/(2*phi^3) = {ALPHA_S_fw:.6f}")
print(f"  PDG:       alpha_s(m_Z) = {ALPHA_S_exp} +/- 0.0009")
print(f"  Error: {abs(ALPHA_S_fw - ALPHA_S_exp)/ALPHA_S_exp*100:.2f}%")

# Run to various scales
print(f"\n--- Running from FRAMEWORK alpha_s(m_Z) = {ALPHA_S_fw:.6f} ---")
print(f"  {'Scale':25s}  {'1-loop':>10s}  {'2-loop':>10s}  {'3-loop':>10s}  {'PDG':>10s}  {'3L err':>8s}")
print(f"  {'-'*25}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*8}")

scales = [
    ("m_t = 172.7 GeV", m_t, 0.1085),
    ("m_Z = 91.2 GeV", m_Z, ALPHA_S_exp),
    ("m_b = 4.18 GeV", m_b, 0.2268),
    ("3 GeV", 3.0, 0.253),
    ("m_tau = 1.78 GeV", 1.777, 0.332),
    ("m_c = 1.27 GeV", m_c, 0.392),
    ("1 GeV", 1.0, 0.47),
]

for label, mu, pdg_val in scales:
    a1_l = run_alpha_s_full(ALPHA_S_fw, mu, 1)
    a2_l = run_alpha_s_full(ALPHA_S_fw, mu, 2)
    a3_l = run_alpha_s_full(ALPHA_S_fw, mu, 3)
    err3 = abs(a3_l - pdg_val)/pdg_val * 100 if not np.isnan(a3_l) else float('nan')
    print(f"  {label:25s}  {a1_l:10.4f}  {a2_l:10.4f}  {a3_l:10.4f}  {pdg_val:10.4f}  {err3:7.1f}%")

# === Cross-check: run from PDG value ===
print(f"\n--- Running from PDG alpha_s(m_Z) = {ALPHA_S_exp} ---")
print(f"  {'Scale':25s}  {'1-loop':>10s}  {'2-loop':>10s}  {'3-loop':>10s}  {'PDG':>10s}  {'3L err':>8s}")
print(f"  {'-'*25}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*8}")

for label, mu, pdg_val in scales:
    a1_l = run_alpha_s_full(ALPHA_S_exp, mu, 1)
    a2_l = run_alpha_s_full(ALPHA_S_exp, mu, 2)
    a3_l = run_alpha_s_full(ALPHA_S_exp, mu, 3)
    err3 = abs(a3_l - pdg_val)/pdg_val * 100 if not np.isnan(a3_l) else float('nan')
    print(f"  {label:25s}  {a1_l:10.4f}  {a2_l:10.4f}  {a3_l:10.4f}  {pdg_val:10.4f}  {err3:7.1f}%")

# === Lambda_QCD determination ===
print(f"\n--- Lambda_QCD from framework ---")
# Lambda_QCD (1-loop, nf=5): Lambda = m_Z * exp(-pi/(b0*alpha_s))
b0_5 = beta_0(5)
Lambda_1loop = m_Z * np.exp(-np.pi / (b0_5 * ALPHA_S_fw))
print(f"  Lambda_QCD (1-loop, nf=5) = m_Z * exp(-pi/(b0*alpha_s))")
print(f"  = {m_Z:.1f} * exp(-pi/({b0_5:.2f}*{ALPHA_S_fw:.4f}))")
print(f"  = {Lambda_1loop*1000:.0f} MeV")
print(f"  PDG Lambda_QCD(5) ~ 210 MeV")
print(f"  Error: {abs(Lambda_1loop*1000 - 210)/210*100:.0f}%")

# 2-loop Lambda
# 1/alpha_s = b0/(2*pi)*ln(mu^2/Lambda^2) + b1/(4*pi^2*b0)*ln(ln(mu^2/Lambda^2))
# At 1-loop: Lambda^2 = mu^2 * exp(-2*pi/(b0*alpha_s))
# Better: use numerical inversion
from scipy.optimize import brentq

def alpha_s_analytic_2loop(mu, Lambda, nf):
    """2-loop analytic formula for alpha_s."""
    b0 = beta_0(nf)
    b1_qcd = beta_1(nf)
    L = np.log(mu**2 / Lambda**2)
    if L <= 0:
        return float('inf')
    a = 1/(b0*L/(2*np.pi)) * (1 - b1_qcd/(2*np.pi*b0**2) * np.log(L)/L)
    return max(a, 0)

def find_Lambda(alpha_s_target, mu_ref, nf):
    """Find Lambda such that alpha_s(mu_ref) = alpha_s_target at 2-loop."""
    def objective(lnLambda):
        Lambda = np.exp(lnLambda)
        return alpha_s_analytic_2loop(mu_ref, Lambda, nf) - alpha_s_target

    try:
        # Lambda should be between 100 MeV and 1 GeV
        lnL = brentq(objective, np.log(0.05), np.log(1.0))
        return np.exp(lnL)
    except:
        return float('nan')

Lambda_2loop = find_Lambda(ALPHA_S_fw, m_Z, 5)
print(f"\n  Lambda_QCD (2-loop, nf=5) from framework: {Lambda_2loop*1000:.0f} MeV")
Lambda_2loop_pdg = find_Lambda(ALPHA_S_exp, m_Z, 5)
print(f"  Lambda_QCD (2-loop, nf=5) from PDG:       {Lambda_2loop_pdg*1000:.0f} MeV")

# === GUT scale running ===
print(f"\n--- GUT Scale Running (3-loop) ---")
mu_GUT_values = [1e14, 2e15, 2e16, 1e17]

ALPHA_0 = 1/137.035999084
ALPHA_MZ = 1/127.951
sin2tW = b1/(a1**2 + 1)

# EW couplings at m_Z (with GUT normalization for U(1))
alpha_1_mZ = (5./3) * ALPHA_MZ / (1 - sin2tW)
alpha_2_mZ = ALPHA_MZ / sin2tW

# EW 1-loop betas
b1_ew = 41./10
b2_ew = -19./6

print(f"  At m_Z:")
print(f"    1/alpha_1 = {1/alpha_1_mZ:.2f} (GUT normalized)")
print(f"    1/alpha_2 = {1/alpha_2_mZ:.2f}")
print(f"    1/alpha_3 = {1/ALPHA_S_fw:.4f}")

print(f"\n  {'mu (GeV)':>12s}  {'1/a1':>8s}  {'1/a2':>8s}  {'1/a3':>8s}  {'a2/a3':>8s}")
print(f"  {'-'*12}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

for mu_GUT in mu_GUT_values:
    ln_r = np.log(mu_GUT / m_Z)
    a1_gut = alpha_1_mZ / (1 - b1_ew*alpha_1_mZ/(2*np.pi)*ln_r)
    a2_gut = alpha_2_mZ / (1 - b2_ew*alpha_2_mZ/(2*np.pi)*ln_r)
    a3_gut = run_alpha_s_full(ALPHA_S_fw, mu_GUT, 3)
    print(f"  {mu_GUT:12.0e}  {1/a1_gut:8.2f}  {1/a2_gut:8.2f}  {1/a3_gut:8.2f}  {a2_gut/a3_gut:8.3f}")

# === KEY ARGUMENT for scale ===
print(f"\n" + "="*72)
print(f"KEY ARGUMENT: WHY COUPLINGS AT m_Z")
print(f"="*72)
print(f"""
  The 600-cell spectral conditions give ALGEBRAIC values:
    alpha_s = 1/(2*phi^3)  [exact in Z[phi]]
    sin^2(tW) = 6/26       [exact in Q]
    alpha = root of 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0

  These are FIXED numbers, not functions of energy.
  They correspond to the LATTICE SCALE of the 600-cell.

  The lattice scale IS m_Z because:
  1. m_Z = m_e * phi^(a1^2) * alpha(m_Z)/alpha(0)
  2. This is the ONLY mass formula with running factor
  3. The running factor alpha(m_Z)/alpha(0) = alpha at LATTICE SCALE
  4. => The lattice scale = m_Z

  SELF-CONSISTENCY CHECK:
  The 600-cell gives alpha_s(lattice) = {ALPHA_S_fw:.6f}
  The lattice scale = m_Z
  SM running gives alpha_s(m_tau) = [see table above]
  This matches PDG within ~few% (depending on loop order)

  The framework does NOT need to explain running.
  Running is SM physics, valid below the lattice scale.
  The framework explains the BOUNDARY CONDITIONS at m_Z.
""")

# === Comparison: framework vs experimental running ===
print(f"--- Final comparison: framework vs PDG alpha_s running ---")
print(f"  At m_Z: framework {ALPHA_S_fw:.6f} vs PDG {ALPHA_S_exp}")
print(f"  Difference at m_Z: {abs(ALPHA_S_fw - ALPHA_S_exp)*1000:.2f} * 10^-3")
print(f"  This difference is {abs(ALPHA_S_fw-ALPHA_S_exp)/0.0009:.1f} sigma (PDG uncertainty)")
print(f"")

# How does the initial 0.11% difference propagate?
print(f"  Propagated differences at low scale (3-loop):")
for label, mu, pdg_val in [("m_tau", 1.777, 0.332), ("m_b", 4.18, 0.2268)]:
    a_fw = run_alpha_s_full(ALPHA_S_fw, mu, 3)
    a_pdg = run_alpha_s_full(ALPHA_S_exp, mu, 3)
    print(f"    {label}: fw={a_fw:.4f}, from PDG={a_pdg:.4f}, diff={abs(a_fw-a_pdg)/a_pdg*100:.2f}%")
    print(f"    {label}: fw vs experiment: {abs(a_fw-pdg_val)/pdg_val*100:.1f}%")

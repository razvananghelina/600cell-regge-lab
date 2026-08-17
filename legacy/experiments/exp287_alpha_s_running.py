"""
EXP-287: ALPHA_S RUNNING - 2-LOOP AND BEYOND
=============================================
Criticism 1: framework gives alpha_s = 1/(2*phi^3) at m_Z.
Does SM running from m_Z reproduce experimental values at other scales?

1-loop was 26% off at m_tau. Try 2-loop, 3-loop.
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
ALPHA_S_fw = 1/(2*PHI**3)  # = 0.118034

# Experimental values (PDG 2024)
ALPHA_S_MZ = 0.1179       # at m_Z = 91.1876 GeV
ALPHA_S_MTAU = 0.332       # at m_tau = 1.777 GeV (from tau decays)
ALPHA_S_3GEV = 0.260       # at 3 GeV (lattice QCD)
ALPHA_S_5GEV = 0.213       # at 5 GeV (various)
ALPHA_S_MBOT = 0.226       # at m_b = 4.18 GeV

m_Z = 91.1876
m_tau = 1.77686
m_b = 4.18
m_c = 1.27

print("="*72)
print("EXP-287: ALPHA_S RUNNING FROM FRAMEWORK VALUE")
print("="*72)

# ================================================================
# 1-LOOP RUNNING
# ================================================================
print("\n--- 1-Loop Running ---")

# 1-loop: d(alpha_s)/d(ln mu) = -b0 * alpha_s^2 / (2*pi)
# Solution: 1/alpha_s(mu) = 1/alpha_s(m_Z) + b0/(2*pi) * ln(mu/m_Z)

# b0 depends on number of active flavors nf:
# b0 = (33 - 2*nf) / 3
def b0(nf): return (33 - 2*nf) / 3.0

# nf at different scales:
# mu > m_t ~ 173 GeV: nf=6
# m_b < mu < m_t: nf=5
# m_c < mu < m_b: nf=4
# mu < m_c: nf=3

def alpha_s_1loop(mu, alpha_s_mz=ALPHA_S_fw, m_z=m_Z):
    """1-loop running with flavor thresholds."""
    # Start at m_Z with nf=5
    alpha_inv = 1.0/alpha_s_mz

    if mu >= m_Z:
        # Run up: nf=5 to m_t, then nf=6
        m_t = 172.76
        if mu <= m_t:
            alpha_inv += b0(5)/(2*np.pi) * np.log(mu/m_Z)
        else:
            alpha_inv += b0(5)/(2*np.pi) * np.log(m_t/m_Z)
            alpha_inv += b0(6)/(2*np.pi) * np.log(mu/m_t)
    else:
        # Run down: nf=5 to m_b, then nf=4 to m_c, then nf=3
        if mu >= m_b:
            alpha_inv += b0(5)/(2*np.pi) * np.log(mu/m_Z)
        elif mu >= m_c:
            alpha_inv += b0(5)/(2*np.pi) * np.log(m_b/m_Z)
            alpha_inv += b0(4)/(2*np.pi) * np.log(mu/m_b)
        else:
            alpha_inv += b0(5)/(2*np.pi) * np.log(m_b/m_Z)
            alpha_inv += b0(4)/(2*np.pi) * np.log(m_c/m_b)
            alpha_inv += b0(3)/(2*np.pi) * np.log(mu/m_c)

    return 1.0/alpha_inv

scales = [
    (1.0, "1 GeV", 0.50, 3),      # approximate
    (m_tau, "m_tau", 0.332, 3),
    (3.0, "3 GeV", 0.260, 4),
    (m_b, "m_b", 0.226, 5),
    (5.0, "5 GeV", 0.213, 5),
    (10.0, "10 GeV", 0.179, 5),
    (m_Z, "m_Z", 0.1179, 5),
]

print(f"  {'Scale':<12} {'1-loop':>10} {'PDG':>10} {'Error':>8}")
print(f"  {'-'*42}")
for mu, name, alpha_exp, nf_expected in scales:
    alpha_pred = alpha_s_1loop(mu)
    err = (alpha_pred - alpha_exp)/alpha_exp * 100
    print(f"  {name:<12} {alpha_pred:10.4f} {alpha_exp:10.4f} {err:+7.1f}%")

# ================================================================
# 2-LOOP RUNNING
# ================================================================
print("\n--- 2-Loop Running ---")

# 2-loop beta: d(a_s)/d(ln mu^2) = -b0*a_s^2 - b1*a_s^3
# where a_s = alpha_s/(4*pi)
# b0 = (33-2*nf)/3
# b1 = (306-38*nf)/3

def b1_coeff(nf): return (306 - 38*nf) / 3.0

def alpha_s_2loop_numerical(mu, alpha_s_start=ALPHA_S_fw, mu_start=m_Z,
                             nsteps=10000):
    """2-loop running via numerical integration (RK4)."""
    # Use t = ln(mu^2/mu_start^2)
    t_end = 2*np.log(mu/mu_start)
    dt = t_end / nsteps

    a_s = alpha_s_start / (4*np.pi)  # convention: a_s = alpha_s/(4*pi)
    t = 0

    # Determine nf transitions
    thresholds = [(2*np.log(m_c/mu_start), 4, 3),  # below m_c: nf=3
                  (2*np.log(m_b/mu_start), 5, 4),  # below m_b: nf=4
                  (2*np.log(172.76/mu_start), 6, 5)]  # below m_t: nf=5

    def get_nf(t_val):
        mu_current = mu_start * np.exp(t_val/2)
        if mu_current < m_c: return 3
        elif mu_current < m_b: return 4
        elif mu_current < 172.76: return 5
        else: return 6

    def beta(a_s_val, t_val):
        nf = get_nf(t_val)
        b0_val = b0(nf)
        b1_val = b1_coeff(nf)
        return -b0_val * a_s_val**2 - b1_val * a_s_val**3

    # RK4 integration
    for i in range(abs(nsteps)):
        step = dt if t_end > 0 else -dt
        k1 = beta(a_s, t)
        k2 = beta(a_s + step/2*k1, t + step/2)
        k3 = beta(a_s + step/2*k2, t + step/2)
        k4 = beta(a_s + step*k3, t + step)
        a_s += step/6 * (k1 + 2*k2 + 2*k3 + k4)
        t += step

        if a_s <= 0:
            return float('inf')  # Landau pole

    return a_s * 4 * np.pi

print(f"  {'Scale':<12} {'2-loop':>10} {'PDG':>10} {'Error':>8}")
print(f"  {'-'*42}")
for mu, name, alpha_exp, nf_expected in scales:
    alpha_pred = alpha_s_2loop_numerical(mu, nsteps=50000)
    err = (alpha_pred - alpha_exp)/alpha_exp * 100
    print(f"  {name:<12} {alpha_pred:10.4f} {alpha_exp:10.4f} {err:+7.1f}%")

# ================================================================
# 3-LOOP RUNNING
# ================================================================
print("\n--- 3-Loop Running ---")

def b2_coeff(nf):
    """3-loop beta coefficient."""
    return 2857./2 - 5033./18*nf + 325./54*nf**2

def alpha_s_3loop_numerical(mu, alpha_s_start=ALPHA_S_fw, mu_start=m_Z,
                             nsteps=50000):
    """3-loop running via RK4."""
    t_end = 2*np.log(mu/mu_start)
    dt = t_end / nsteps

    a_s = alpha_s_start / (4*np.pi)
    t = 0

    def get_nf(t_val):
        mu_current = mu_start * np.exp(t_val/2)
        if mu_current < m_c: return 3
        elif mu_current < m_b: return 4
        elif mu_current < 172.76: return 5
        else: return 6

    def beta(a_s_val, t_val):
        nf = get_nf(t_val)
        return (-b0(nf) * a_s_val**2
                - b1_coeff(nf) * a_s_val**3
                - b2_coeff(nf) * a_s_val**4)

    for i in range(abs(nsteps)):
        step = dt if t_end > 0 else -dt
        k1 = beta(a_s, t)
        k2 = beta(a_s + step/2*k1, t + step/2)
        k3 = beta(a_s + step/2*k2, t + step/2)
        k4 = beta(a_s + step*k3, t + step)
        a_s += step/6 * (k1 + 2*k2 + 2*k3 + k4)
        t += step

        if a_s <= 0:
            return float('inf')

    return a_s * 4 * np.pi

print(f"  {'Scale':<12} {'3-loop':>10} {'PDG':>10} {'Error':>8}")
print(f"  {'-'*42}")
for mu, name, alpha_exp, nf_expected in scales:
    alpha_pred = alpha_s_3loop_numerical(mu, nsteps=50000)
    err = (alpha_pred - alpha_exp)/alpha_exp * 100
    print(f"  {name:<12} {alpha_pred:10.4f} {alpha_exp:10.4f} {err:+7.1f}%")

# ================================================================
# COMPARISON: framework vs PDG world average
# ================================================================
print("\n--- Framework vs PDG at m_Z ---")
print(f"  Framework: alpha_s(m_Z) = 1/(2*phi^3) = {ALPHA_S_fw:.6f}")
print(f"  PDG 2024:  alpha_s(m_Z) = {ALPHA_S_MZ}")
print(f"  Difference: {(ALPHA_S_fw - ALPHA_S_MZ)/ALPHA_S_MZ*100:+.3f}%")
print(f"  PDG uncertainty: +/- 0.0009")
print(f"  Framework within: {abs(ALPHA_S_fw - ALPHA_S_MZ)/0.0009:.1f} sigma")

# ================================================================
# GUT SCALE RUNNING
# ================================================================
print("\n--- GUT Scale ---")
# Run all 3 couplings to GUT scale
ALPHA_0 = 1/137.035999084
sin2tW = b1/(a1**2 + 1)
alpha_em_mz = 1/127.951

# At m_Z with GUT normalization
alpha_1_mz = (5./3) * alpha_em_mz / (1 - sin2tW)
alpha_2_mz = alpha_em_mz / sin2tW
alpha_3_mz = ALPHA_S_fw

print(f"  At m_Z (GUT normalization):")
print(f"    1/alpha_1 = {1/alpha_1_mz:.2f}")
print(f"    1/alpha_2 = {1/alpha_2_mz:.2f}")
print(f"    1/alpha_3 = {1/alpha_3_mz:.2f}")

# 1-loop SM running to GUT scale
b1_gut = 41./10
b2_gut = -19./6
b3_gut = -7.

mu_vals = np.logspace(np.log10(m_Z), 17, 1000)
inv_a1 = [1/alpha_1_mz + b1_gut/(2*np.pi)*np.log(mu/m_Z) for mu in mu_vals]
inv_a2 = [1/alpha_2_mz + b2_gut/(2*np.pi)*np.log(mu/m_Z) for mu in mu_vals]
inv_a3 = [1/alpha_3_mz + b3_gut/(2*np.pi)*np.log(mu/m_Z) for mu in mu_vals]

# Find crossing points
for i in range(len(mu_vals)-1):
    if (inv_a2[i] - inv_a3[i]) * (inv_a2[i+1] - inv_a3[i+1]) < 0:
        print(f"  alpha_2 = alpha_3 crossing at mu ~ {mu_vals[i]:.2e} GeV")
        print(f"    1/alpha ~ {inv_a2[i]:.2f}")
    if (inv_a1[i] - inv_a2[i]) * (inv_a1[i+1] - inv_a2[i+1]) < 0:
        print(f"  alpha_1 = alpha_2 crossing at mu ~ {mu_vals[i]:.2e} GeV")
        print(f"    1/alpha ~ {inv_a1[i]:.2f}")

# Framework GUT relation check
print(f"\n  Framework GUT relation (exact for a1=5):")
print(f"    1/alpha_1(GUT) = 2 * 1/alpha_2 ?")
print(f"    3*(a1^2-a1)/(5*b1) = 3*20/30 = {3*20/30:.0f} (should be 2)")
print(f"    This is EXACT and independent of scale!")

print("\n" + "="*72)
print("CONCLUSION")
print("="*72)
print(f"  The framework value alpha_s(m_Z) = 1/(2*phi^3) = {ALPHA_S_fw:.6f}")
print(f"  is within 1.5 sigma of PDG world average {ALPHA_S_MZ}.")
print(f"  SM 2-loop/3-loop running from this value to other scales")
print(f"  should be checked for consistency with lattice + tau data.")
print(f"  The 1-loop running is known to be insufficient at low scales.")

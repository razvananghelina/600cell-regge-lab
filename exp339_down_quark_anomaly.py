"""
exp339: The Down Quark Anomaly -- Why is b=0 Special?
=====================================================
The down quark is the ONLY fermion with b=0 (besides e at a=b=0).
It's Galois-invariant (z=z'=1) and the biggest outlier in every model.

Systematic investigation of:
Q1: delta_d vs sin^2(theta_W)
Q2: CKM mixing effect
Q4: Antipodal anomaly on 600-cell
Q6/Q7: delta_d = -1/phi^2 ?
Q5: Systematic framework expression search
Q8: Combined model
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
PHIp = (1 - np.sqrt(5)) / 2  # = -1/PHI
ALPHA_S = 1 / (2 * PHI**3)
ALPHA = 1/137.036
SIN2TW = 6/26
A1 = 5
B1 = 6
N_GEN = 3
N_EIG = 9

# PDG masses (MeV)
m_e = 0.51100
m_mu = 105.658
m_tau = 1776.86
m_u = 2.16
m_c = 1270.0
m_t = 172760.0
m_d = 4.67
m_s = 93.4
m_b = 4180.0

pdg_masses = np.array([m_e, m_mu, m_tau, m_u, m_c, m_t, m_d, m_s, m_b])
fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
a_vals = np.array([0, 1, 1, 3, 2, 4, 1, 1, -1])
b_vals = np.array([0, 1, 2, -2, 1, 1, 0, 1, 4])
n_bare = 5 * a_vals + 6 * b_vals  # [0, 11, 17, 3, 16, 26, 5, 11, 19]
z_vals = a_vals + b_vals * PHI

# Effective exponents
n_eff = np.log(pdg_masses / m_e) / np.log(PHI)
delta = n_eff - n_bare

# Galois norms N(z) = a^2 + ab - b^2
norms = a_vals**2 + a_vals*b_vals - b_vals**2

print("="*72)
print("exp339: The Down Quark Anomaly")
print("="*72)

print("\nFermion data:")
print(f"{'Name':>5} {'(a,b)':>8} {'n':>4} {'n_eff':>8} {'delta':>8} {'N(z)':>6} {'z':>10} {'b=0':>5}")
print("-"*65)
for i in range(9):
    z_str = f"{a_vals[i]}+{b_vals[i]}*phi" if b_vals[i] != 0 else f"{a_vals[i]}"
    b0 = "YES" if b_vals[i] == 0 else ""
    print(f"{fermion_names[i]:>5} ({a_vals[i]:>2},{b_vals[i]:>2}) {n_bare[i]:>4} {n_eff[i]:>8.3f} {delta[i]:>+8.4f} {norms[i]:>+6d} {z_str:>10} {b0:>5}")

# ================================================================
# Q1: delta_d vs sin^2(theta_W)
# ================================================================
print(f"\n{'='*72}")
print("Q1: delta_d vs electroweak parameters")
print("="*72)

d_delta = delta[6]
print(f"\ndelta_d = {d_delta:.4f}")
print(f"\nComparisons:")
candidates_q1 = {
    'sin2tW': SIN2TW,
    '2*sin2tW': 2*SIN2TW,
    'alpha_s': ALPHA_S,
    '2*alpha_s': 2*ALPHA_S,
    '3*alpha_s': 3*ALPHA_S,
    '1/phi^2': 1/PHI**2,
    'phi-1': PHI-1,
    '2-phi': 2-PHI,
    '1/phi': 1/PHI,
    'alpha_s*pi': ALPHA_S*np.pi,
    '1/a1': 1/A1,
    'b1/(a1^2+1)': B1/(A1**2+1),
    'alpha_s+sin2tW': ALPHA_S+SIN2TW,
    '1/(2*phi)': 1/(2*PHI),
    'sin2tW+alpha_s': SIN2TW+ALPHA_S,
    '(phi-1)^2': (PHI-1)**2,
}

print(f"{'Expression':>25} {'Value':>10} {'|delta_d|':>10} {'Diff':>10} {'%err':>8}")
print("-"*70)
for name, val in sorted(candidates_q1.items(), key=lambda x: abs(abs(d_delta)-x[1])):
    diff = abs(d_delta) - val
    pct = 100 * diff / abs(d_delta)
    print(f"{name:>25} {val:>10.4f} {abs(d_delta):>10.4f} {diff:>+10.4f} {pct:>+8.1f}%")

# ================================================================
# Q7: Is delta_d = -1/phi^2 exactly?
# ================================================================
print(f"\n{'='*72}")
print("Q7: delta_d = -1/phi^2 = -(phi-1) = -(2-phi)?")
print("="*72)

print(f"\ndelta_d = {d_delta:.6f}")
print(f"-1/phi^2 = {-1/PHI**2:.6f}")
print(f"Difference: {d_delta - (-1/PHI**2):.6f}")
print(f"Relative error: {100*(d_delta + 1/PHI**2)/abs(d_delta):.2f}%")

# If delta_d = -1/phi^2, then:
n_eff_pred = 5 - 1/PHI**2
m_pred = m_e * PHI**n_eff_pred
print(f"\nIf n_eff(d) = 5 - 1/phi^2 = {n_eff_pred:.4f}:")
print(f"  m_d(pred) = {m_pred:.3f} MeV (PDG: {m_d:.2f} MeV)")
print(f"  Error: {100*(m_pred-m_d)/m_d:.1f}%")

# Other phi expressions for delta_d
print(f"\nOther Z[phi] candidates for delta_d:")
phi_candidates = {
    '-1/phi^2 = -(phi-1)': -1/PHI**2,
    '-phi/phi^3 = -1/phi^2': -PHI/PHI**3,
    '-(phi^2-phi-1)/phi': -(PHI**2-PHI-1)/PHI,  # = 0/phi = 0 (phi^2=phi+1)
    '-alpha_s*phi^2': -ALPHA_S*PHI**2,
    '-alpha_s*(a1-1)': -ALPHA_S*(A1-1),
    '-2*alpha_s*phi': -2*ALPHA_S*PHI,
    '-1/(phi*a1)': -1/(PHI*A1),
    '-alpha_s*N_gen': -ALPHA_S*N_GEN,
    '-sin2tW*phi': -SIN2TW*PHI,
    '-1/(2*phi^2)*(1+phi)': -1/(2*PHI**2)*(1+PHI),  # = -phi/(2*phi^2) = -1/(2*phi)
    '-phi/(a1+1)': -PHI/(A1+1),
    '-phi^2/(a1^2+1)': -PHI**2/(A1**2+1),
    '-(a1-1)/(a1^2+1)': -(A1-1)/(A1**2+1),
    '-4/(a1^2+1)': -4/(A1**2+1),
    '-2/(a1+phi^2)': -2/(A1+PHI**2),
    '-N_gen/phi^4': -N_GEN/PHI**4,
    '-1/(phi*phi^2)': -1/(PHI*PHI**2),  # = -1/phi^3 = -alpha_s*2
    '-phi/(2*a1)': -PHI/(2*A1),
    '-b1/(a1^2*phi)': -B1/(A1**2*PHI),
}

print(f"{'Expression':>30} {'Value':>10} {'delta_d':>10} {'diff':>10} {'%err':>8}")
print("-"*75)
for name, val in sorted(phi_candidates.items(), key=lambda x: abs(d_delta-x[1])):
    diff = d_delta - val
    pct = 100 * diff / abs(d_delta) if abs(d_delta) > 0 else 0
    print(f"{name:>30} {val:>10.6f} {d_delta:>10.6f} {diff:>+10.6f} {pct:>+8.2f}%")

# ================================================================
# SYSTEMATIC SEARCH: f(a1, b1, phi, alpha_s, sin2tW) for delta_d
# ================================================================
print(f"\n{'='*72}")
print("SYSTEMATIC SEARCH: framework expressions for |delta_d|")
print("="*72)

target = abs(d_delta)  # 0.4021

# Build dictionary of "atoms" and try all ratios/products of 2 atoms
atoms = {
    '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6,
    'phi': PHI, '1/phi': 1/PHI, 'phi^2': PHI**2, '1/phi^2': 1/PHI**2,
    'phi^3': PHI**3, '1/phi^3': 1/PHI**3, 'phi^4': PHI**4,
    'alpha_s': ALPHA_S, 'sin2tW': SIN2TW, 'pi': np.pi,
    'sqrt5': np.sqrt(5), 'a1': A1, 'b1': B1, 'N_gen': N_GEN,
    'a1^2': A1**2, 'a1+1': A1+1, 'a1-1': A1-1,
    'a1^2+1': A1**2+1,
}

# Single ratios: a/b
best_single = []
for n1, v1 in atoms.items():
    for n2, v2 in atoms.items():
        if abs(v2) > 1e-10 and n1 != n2:
            val = v1 / v2
            if 0.01 < abs(val) < 10:
                err = abs(val - target) / target
                if err < 0.03:  # within 3%
                    best_single.append((err, f"{n1}/{n2}", val))

best_single.sort()
print(f"\nBest single ratios (within 3% of |delta_d|={target:.4f}):")
print(f"{'Expression':>25} {'Value':>10} {'%err':>8}")
print("-"*50)
for err, name, val in best_single[:20]:
    print(f"{name:>25} {val:>10.6f} {100*err:>+8.3f}%")

# Products of 2: a*b
best_prod = []
for n1, v1 in atoms.items():
    for n2, v2 in atoms.items():
        val = v1 * v2
        if 0.01 < abs(val) < 10:
            err = abs(val - target) / target
            if err < 0.03:
                name = f"{n1}*{n2}"
                best_prod.append((err, name, val))

best_prod.sort()
print(f"\nBest products (within 3%):")
print(f"{'Expression':>30} {'Value':>10} {'%err':>8}")
print("-"*55)
seen = set()
for err, name, val in best_prod[:30]:
    # Avoid duplicates like a*b and b*a
    key = tuple(sorted(name.split('*')))
    if key not in seen:
        seen.add(key)
        print(f"{name:>30} {val:>10.6f} {100*err:>+8.3f}%")

# ================================================================
# Q2: CKM mixing effect on down quark mass
# ================================================================
print(f"\n{'='*72}")
print("Q2: CKM mixing correction")
print("="*72)

# Framework CKM parameters
# V_us = phi^(-n_12) with n_12 = N_gen = 3
# V_us ~ phi^(-3) ~ 0.236 (PDG: 0.2243)
# V_cb ~ phi^(-7) ~ 0.034 (PDG: 0.0422)
# V_ub ~ phi^(-12) ~ 0.0035 (PDG: 0.00394)

V_us_fw = PHI**(-N_GEN)  # phi^-3
V_cb_fw = PHI**(-7)
V_ub_fw = PHI**(-12)
V_ud_fw = np.sqrt(1 - V_us_fw**2)
V_cs_fw = np.sqrt(1 - V_us_fw**2 - V_cb_fw**2)
V_td_fw = V_cb_fw  # approx

print(f"\nFramework CKM elements:")
print(f"  |V_ud| = {V_ud_fw:.4f} (PDG: 0.9743)")
print(f"  |V_us| = {V_us_fw:.4f} (PDG: 0.2243)")
print(f"  |V_ub| = {V_ub_fw:.6f} (PDG: 0.00394)")

# In the mass matrix, the physical d mass is shifted by off-diagonal elements
# The leading correction to m_d from CKM mixing is:
# delta_m_d / m_d ~ |V_us|^2 * (m_s/m_d - 1) + |V_ub|^2 * (m_b/m_d - 1)
# This is a perturbative estimate; the actual effect comes from diagonalization

# Bare masses (from framework)
m_d_bare = m_e * PHI**5
m_s_bare = m_e * PHI**11
m_b_bare = m_e * PHI**19

print(f"\nBare masses: m_d={m_d_bare:.2f}, m_s={m_s_bare:.1f}, m_b={m_b_bare:.0f} MeV")

# Simple CKM correction estimate:
# The mass matrix eigenvalue shift for the lightest eigenvalue:
# delta(m_d) ~ -|V_us|^2 * m_s - |V_ub|^2 * m_b (Gell-Mann-Oakes-Renner-like)
# But this is for mass-SQUARED matrices in Higgs sector

# Actually, in the Yukawa sector:
# M_down = V_CKM^dag * diag(y_d, y_s, y_b) * v
# The "geometric" masses are the diagonal elements. The physical masses
# are the eigenvalues of M_down^dag * M_down.
# For nearly diagonal case: m_d(phys) ~ y_d*v * (1 - |V_us|^2/2 * ...)

# Let's compute the actual effect numerically
# If the bare mass matrix is M = diag(m_d_bare, m_s_bare, m_b_bare)
# and CKM mixing rotates it: M_phys = V * M * V^dag (simplified)
# Eigenvalues of M_phys are still the same as M (CKM doesn't change eigenvalues!)

# Wait - this is wrong. CKM connects UP-type to DOWN-type, it doesn't mix
# the down-type masses among themselves. The mass matrix diagonalization
# IS what defines the CKM matrix.

# The correct statement: if the framework gives "geometric" masses, these
# ARE the physical mass eigenvalues. CKM is derived FROM the misalignment
# of up and down sector mass matrices. So CKM doesn't "correct" individual masses.

print("\nCKM analysis:")
print("  CKM mixing does NOT shift mass eigenvalues.")
print("  CKM is DERIVED from the mass matrix misalignment.")
print("  The down quark anomaly is NOT a CKM effect.")
print("  (CKM angles are already derived separately in the framework.)")

# ================================================================
# Q4: Antipodal point on 600-cell
# ================================================================
print(f"\n{'='*72}")
print("Q4: Antipodal structure on the 600-cell")
print("="*72)

# Distance classes on 600-cell:
# The 120 vertices of the 600-cell, viewed as quaternions in S^3,
# have angular distances cos(angle) = inner product.
# From any vertex, the other 119 vertices fall into distance classes:
# d=0: 1 vertex (itself)
# d=1: 12 vertices (nearest neighbors, angle=pi/5)
# d=2: 20 vertices (angle=2*pi/5? or pi/3?)
# ...
# d=5: 1 vertex (antipodal, angle=pi)

# The 600-cell has 9 distance classes corresponding to 9 2I irreps.
# The antipodal point is at maximum distance.

# From the character table, the contribution of each irrep at the antipodal point:
# chi_k(-1)/d_k = (-1)^{2j}/1 for spin-j reps
# For the main chain (j=0,...,5/2): chi_k(-1) = (-1)^{2j} * d_k
# = d_k for even 2j, -d_k for odd 2j

# Class C_2 (the central element -1, i.e., antipodal):
chartab_C2 = np.array([1, -2, 3, -4, 5, -6, 4, -2, 3])
dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])

print("\nCharacters at antipodal point (C_2 = {-1}):")
for k in range(9):
    ratio = chartab_C2[k] / dims[k]
    print(f"  rho_{k} (dim {dims[k]}): chi(-1) = {chartab_C2[k]:>3}, chi(-1)/d = {ratio:>+.3f}")

# The Green's function at the antipodal point:
# G(-1) = (1/N) * sum_k d_k^2 * chi_k(-1)/(d_k * (lambda_k + mu^2))
# = (1/N) * sum_k d_k * chi_k(-1)/(lambda_k + mu^2)

char_C10a = [1, PHI, PHI, 1, 0, -1, -1, -1/PHI, -1/PHI]
cayley_eigs = np.array([12 * (1 - char_C10a[k] / dims[k]) for k in range(9)])

print(f"\nGreen's function at antipodal point vs self:")
for mu_sq in [0.1, 1.0, 5.0, 10.0]:
    G_self = sum(dims[k]**2 / (cayley_eigs[k] + mu_sq) for k in range(9)) / 120
    G_anti = sum(dims[k] * chartab_C2[k] / (cayley_eigs[k] + mu_sq) for k in range(9)) / 120
    print(f"  mu^2={mu_sq:>5.1f}: G(self)={G_self:.4f}, G(anti)={G_anti:.4f}, ratio={G_anti/G_self:.4f}")

# Key: at the antipodal point, odd-dim irreps contribute POSITIVELY
# and even-dim contribute NEGATIVELY (or vice versa).
# The ALTERNATING SIGN structure at the antipodal point could give
# a qualitatively different correction.

print(f"\n  The antipodal Green's function has ALTERNATING signs")
print(f"  from chi_k(-1) = (-1)^(2j) * d_k")
print(f"  This is NOT the same as a simple distance correction.")

# ================================================================
# Q6/Q7: The key test - framework expressions for delta_d
# ================================================================
print(f"\n{'='*72}")
print("Q6/Q7: BEST framework expressions for delta_d")
print("="*72)

target_d = d_delta  # -0.4021

# Extended search with 3-atom expressions: a*b/c, a/(b*c), (a+b)/c, etc.
best_3 = []

atoms_small = {
    '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6,
    'phi': PHI, '1/phi': 1/PHI, 'phi^2': PHI**2, '1/phi^2': 1/PHI**2,
    'alpha_s': ALPHA_S, 'sin2tW': SIN2TW,
    'a1': 5, 'b1': 6, 'N_gen': 3, 'N_eig': 9,
    'a1+1': 6, 'a1^2+1': 26,
}

# Check sums/differences of 2 atoms
for n1, v1 in atoms_small.items():
    for n2, v2 in atoms_small.items():
        if n1 <= n2:
            continue
        for sign in [+1, -1]:
            val = -(v1 + sign * v2)  # negative because delta_d < 0
            if 0.3 < abs(val) < 0.5:
                err = abs(val - target_d) / abs(target_d)
                if err < 0.05:
                    op = '+' if sign > 0 else '-'
                    best_3.append((err, f"-({n1}{op}{n2})", val))

# Check ratios with sign
for n1, v1 in atoms_small.items():
    for n2, v2 in atoms_small.items():
        if abs(v2) < 0.01:
            continue
        val = -v1 / v2
        if 0.3 < abs(val) < 0.5:
            err = abs(val - target_d) / abs(target_d)
            if err < 0.05:
                best_3.append((err, f"-{n1}/{n2}", val))

best_3.sort()
print(f"\nBest expressions for delta_d = {target_d:.6f} (within 5%):")
print(f"{'Expression':>30} {'Value':>12} {'%err':>8}")
print("-"*55)
seen = set()
for err, name, val in best_3[:30]:
    if name not in seen:
        seen.add(name)
        print(f"{name:>30} {val:>12.6f} {100*err:>+8.3f}%")

# ================================================================
# Q8: Can we fix ALL deltas with separate mechanisms?
# ================================================================
print(f"\n{'='*72}")
print("Q8: Combined model - separate mechanisms per sector")
print("="*72)

# Observations:
# - Leptons: deltas are small (0, +0.08, -0.055). Consistent with zero.
# - Up quarks: deltas grow with generation (0, +0.25, +0.46)
# - Down quarks: large negative (-0.40, -0.18, -0.28)

# Pattern in up quarks: delta_up ~ alpha_s * g*(g+1) ?
print("\nUp quark deltas vs generation Casimir g*(g+1):")
up_deltas = [delta[3], delta[4], delta[5]]  # u, c, t
for g in range(3):
    cas = g * (g + 1)
    pred = ALPHA_S * cas * (A1-1)  # scale by something
    print(f"  g={g}: delta={up_deltas[g]:>+.4f}, g*(g+1)={cas}, alpha_s*cas*(a1-1)={pred:.4f}")

# Down quarks: delta_down ~ -alpha_s * (something + generation)
print(f"\nDown quark deltas:")
down_deltas = [delta[6], delta[7], delta[8]]  # d, s, b
for g in range(3):
    print(f"  g={g}: delta={down_deltas[g]:>+.4f}")

# Try: delta_up(g) = alpha_s * g*(g+1) * C_up
# delta_down(g) = -1/phi^2 * (1 - g*alpha_s)? or something similar

# Fit: delta_quark = T3 * [A + B*g + C*g^2]
gens = np.array([0, 1, 2, 0, 1, 2])  # generation for u,c,t,d,s,b
T3_q = np.array([0.5, 0.5, 0.5, -0.5, -0.5, -0.5])
delta_q = np.array([delta[3], delta[4], delta[5], delta[6], delta[7], delta[8]])

# Design matrix: T3*(1, g, g^2)
M = np.column_stack([T3_q, T3_q*gens, T3_q*gens**2])
result = np.linalg.lstsq(M, delta_q, rcond=None)
params = result[0]
pred_q = M @ params
rms_q = np.sqrt(np.mean((delta_q - pred_q)**2))

print(f"\nFit: delta = T3 * (A + B*g + C*g^2)")
print(f"  A = {params[0]:.4f}, B = {params[1]:.4f}, C = {params[2]:.4f}")
print(f"  RMS (quarks only) = {rms_q:.4f}")
for i, name in enumerate(['u', 'c', 't', 'd', 's', 'b']):
    print(f"  {name:>4}: pred={pred_q[i]:>+.4f}, exp={delta_q[i]:>+.4f}")

# Check if A, B, C are framework numbers
print(f"\n  A = {params[0]:.4f} ~= ? (compare: phi/2={PHI/2:.4f}, 2*alpha_s={2*ALPHA_S:.4f})")
print(f"  B = {params[1]:.4f} ~= ? (compare: alpha_s={ALPHA_S:.4f})")
print(f"  C = {params[2]:.4f} ~= ? (compare: alpha_s/2={ALPHA_S/2:.4f})")

# ================================================================
# CRITICAL TEST: delta_d = -(phi-1) corrected by running
# ================================================================
print(f"\n{'='*72}")
print("CRITICAL: delta_d and the Galois boundary")
print("="*72)

# The key insight: for b=0, the Galois conjugate mass is:
# m_d' = m_e * phi'^5 = m_e * (-1/phi)^5 = -m_e/phi^5
# |m_d'| = m_e/phi^5 = m_e * phi^(-5)
# So m_d * |m_d'| = m_e^2 (the mass duality!)

# The "Galois correction" for other fermions comes from the
# difference between phi^n and the Galois-invariant norm.
# For b=0: the correction is ZERO from Galois (z'=z).

# But there IS a correction from the EXPONENT being in the
# irrational part: phi^5 is irrational even though z=1 is rational.

# The continued fraction of phi^5:
phi5 = PHI**5
print(f"\nphi^5 = {phi5:.6f}")
print(f"Nearest integers: 11 (diff={phi5-11:.4f})")
print(f"phi^5 = 5*phi + 3 = {5*PHI+3:.6f}")
print(f"phi^5 = F(7) + F(5)*phi = 13 + 5*phi... no")
print(f"Actually: phi^5 = 8*phi + 5 = 5 + 8*phi = 5 + 8*{PHI:.4f} = {5+8*PHI:.4f}... no")
# phi^5 = F(5)*phi + F(4) = 5*phi + 3 = 5*1.618+3 = 11.09
print(f"phi^5 = F(5)*phi + F(4) = 5*phi + 3 = {5*PHI+3:.6f}")

# The Galois conjugate: phi'^5 = F(5)*phi' + F(4) = 5*phi' + 3
# phi'^5 = 5*(-0.618) + 3 = -3.09 + 3 = -0.09
phi_p_5 = PHIp**5
print(f"phi'^5 = {phi_p_5:.6f}")
print(f"|phi'^5| = {abs(phi_p_5):.6f}")

# The "Galois splitting" of phi^5:
# phi^5 - phi'^5 = (5*phi+3) - (5*phi'+3) = 5*(phi-phi') = 5*sqrt(5) = 11.180
# phi^5 + phi'^5 = 2*3 + 5*(phi+phi') = 6 + 5*1 = 11 = L(5) (Lucas number!)
galois_diff = PHI**5 - PHIp**5
galois_sum = PHI**5 + PHIp**5
print(f"\nphi^5 - phi'^5 = {galois_diff:.4f} = 5*sqrt(5) = {5*np.sqrt(5):.4f}")
print(f"phi^5 + phi'^5 = {galois_sum:.4f} = L(5) = 11")

# Ratio: phi'^5/phi^5 = (phi'/phi)^5 = (-1/phi^2)^5 = -1/phi^10
ratio = phi_p_5 / phi5
print(f"\nphi'^5/phi^5 = {ratio:.6f} = -1/phi^10 = {-1/PHI**10:.6f}")

# The correction delta_d could be related to this ratio:
# delta_d = -|phi'^5/phi^5| = -1/phi^10?
print(f"\n-1/phi^10 = {-1/PHI**10:.6f} vs delta_d = {d_delta:.6f}")
print(f"  Error: {100*(d_delta+1/PHI**10)/abs(d_delta):.1f}% (way off)")

# Or: delta_d = -phi'^5/phi^5... but that's = +1/phi^10 (positive, wrong sign)

# Better: phi^5 = L(5)/2 + F(5)*sqrt(5)/2 = 11/2 + 5*sqrt(5)/2
# The "irrational part" is 5*sqrt(5)/2 = 5.59
# The "rational part" is 11/2 = 5.5
# Ratio: irrational/rational = 5*sqrt(5)/11 = sqrt(5)*5/11

# For the down quark, the "correction" could be the ratio of
# irrational to total: F(5)*sqrt(5) / (2*phi^5)
irr_frac = 5*np.sqrt(5) / (2*phi5)
print(f"\nIrrational fraction of phi^5: {irr_frac:.6f}")
print(f"  vs |delta_d| = {abs(d_delta):.6f}")

# ================================================================
# THE COMBINED MODEL
# ================================================================
print(f"\n{'='*72}")
print("COMBINED MODEL: Best fit with derived structure")
print("="*72)

# Model: delta_f = T3 * Nc * alpha_s * [A + B*g(g+1)] + C * delta_{b=0}
# where delta_{b=0} = 1 for d quark only
# This has 3 params: A, B, C (or derive C = -1/phi^2)

# For leptons: T3*Nc = 0 -> delta = 0 (correct to ~3%)
# For quarks: T3*Nc*alpha_s*[A + B*g(g+1)]
# + C for d only

# First: fit with C = -1/phi^2 (fixed, "derived")
delta_b0 = np.zeros(9)
delta_b0[6] = 1  # d quark only

delta_q_all = delta.copy()

# Fit: delta = a * T3*Nc + b * T3*Nc*g*(g+1) + c * delta_{b=0}
gens_all = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
# T3: 0 for leptons, +1/2 for up quarks, -1/2 for down quarks
fermion_T3 = np.array([0, 0, 0, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5])
# Nc: 0 for leptons, 1 for quarks (factor)
T3Nc = fermion_T3 * np.array([0, 0, 0, 1, 1, 1, 1, 1, 1])
casimir = gens_all * (gens_all + 1)

# With C fixed
C_fixed = -1/PHI**2
delta_corrected = delta_q_all - C_fixed * delta_b0

M3 = np.column_stack([T3Nc, T3Nc * casimir])
result3 = np.linalg.lstsq(M3, delta_corrected, rcond=None)
A, B = result3[0]
pred3 = M3 @ result3[0] + C_fixed * delta_b0
rms3 = np.sqrt(np.mean((delta_q_all - pred3)**2))

print(f"\nModel: delta = A*T3*Nc + B*T3*Nc*g(g+1) + (-1/phi^2)*delta_{{b=0}}")
print(f"  A = {A:.4f}, B = {B:.4f}, C = -1/phi^2 = {C_fixed:.4f} (FIXED)")
print(f"  RMS = {rms3:.4f}")
print(f"\n  {'Name':>5} {'pred':>8} {'exp':>8} {'err':>8}")
for f in range(9):
    print(f"  {fermion_names[f]:>5} {pred3[f]:>+8.4f} {delta_q_all[f]:>+8.4f} {pred3[f]-delta_q_all[f]:>+8.4f}")

# Check: are A and B framework numbers?
print(f"\n  A = {A:.4f}")
print(f"    A/alpha_s = {A/ALPHA_S:.4f}")
print(f"    Compare: phi = {PHI:.4f}, 2 = 2, a1-1 = 4")
print(f"  B = {B:.4f}")
print(f"    B/alpha_s = {B/ALPHA_S:.4f}")
print(f"    Compare: 1 = 1, phi-1 = {PHI-1:.4f}")

# Now fit with C free too
M4 = np.column_stack([T3Nc, T3Nc * casimir, delta_b0])
result4 = np.linalg.lstsq(M4, delta_q_all, rcond=None)
A4, B4, C4 = result4[0]
pred4 = M4 @ result4[0]
rms4 = np.sqrt(np.mean((delta_q_all - pred4)**2))

print(f"\nFree C: A={A4:.4f}, B={B4:.4f}, C={C4:.4f}, RMS={rms4:.4f}")
print(f"  C = {C4:.4f} vs -1/phi^2 = {-1/PHI**2:.4f}")
print(f"  Difference in C: {C4+1/PHI**2:.4f}")

# Best possible with T3*Nc structure + free down correction
# delta = a*T3*Nc + b*T3*Nc*g + c*T3*Nc*g^2 + d*delta_{b=0}
M5 = np.column_stack([T3Nc, T3Nc*gens_all, T3Nc*gens_all**2, delta_b0])
result5 = np.linalg.lstsq(M5, delta_q_all, rcond=None)
pred5 = M5 @ result5[0]
rms5 = np.sqrt(np.mean((delta_q_all - pred5)**2))

print(f"\n4-param model (T3*Nc*[1,g,g^2] + delta_b0): RMS={rms5:.4f}")
print(f"  Params: {np.round(result5[0], 4)}")
for f in range(9):
    print(f"  {fermion_names[f]:>5} {pred5[f]:>+8.4f} {delta_q_all[f]:>+8.4f} {pred5[f]-delta_q_all[f]:>+8.4f}")

# ================================================================
# MASS PREDICTIONS
# ================================================================
print(f"\n{'='*72}")
print("MASS PREDICTIONS: Combined model")
print("="*72)

# Use the 3-param model with C = -1/phi^2 (1 derived + 2 fit)
print(f"\nModel: delta = {A:.4f}*T3*Nc + {B:.4f}*T3*Nc*g(g+1) + (-1/phi^2)*1_{{b=0}}")
print(f"\n{'Name':>5} {'n_bare':>6} {'n_eff':>8} {'m_pred(MeV)':>12} {'m_PDG(MeV)':>12} {'%err':>8}")
print("-"*55)
for f in range(9):
    n_pred = n_bare[f] + pred3[f]
    m_pred = m_e * PHI**n_pred
    pct = 100 * (m_pred - pdg_masses[f]) / pdg_masses[f]
    print(f"{fermion_names[f]:>5} {n_bare[f]:>6} {n_pred:>8.3f} {m_pred:>12.3f} {pdg_masses[f]:>12.3f} {pct:>+8.1f}%")

rms_mass_pct = np.sqrt(np.mean([(100*(m_e*PHI**(n_bare[f]+pred3[f]) - pdg_masses[f])/pdg_masses[f])**2 for f in range(9)]))
print(f"\nRMS mass error: {rms_mass_pct:.1f}%")

print("\n\nDone.")

# ================================================================
# ZERO-PARAMETER MODEL TEST
# ================================================================
print(f"\n{'='*72}")
print("ZERO-PARAMETER MODEL: delta = (g/2)*T3 - (2/a1)*1_{b=0}")
print("="*72)

print("\nAll coefficients are framework numbers:")
print(f"  g/2 * T3: generation number / 2, times isospin")
print(f"  2/a1 = 2/{A1} = {2/A1:.4f} for the b=0 rational sector")
print()

pred_zero = np.zeros(9)
for f in range(9):
    g = gens_all[f]
    T3 = fermion_T3[f]
    is_b0 = 1 if (f == 6) else 0  # only d quark
    pred_zero[f] = (g / 2) * T3 - (2 / A1) * is_b0

print(f"{'Name':>5} {'g':>3} {'T3':>5} {'b=0':>4} {'pred':>8} {'actual':>8} {'diff':>8} {'%m_err':>8}")
print("-"*60)
for f in range(9):
    n_pred = n_bare[f] + pred_zero[f]
    m_pred = m_e * PHI**n_pred
    pct_m = 100 * (m_pred - pdg_masses[f]) / pdg_masses[f]
    print(f"{fermion_names[f]:>5} {gens_all[f]:>3} {fermion_T3[f]:>+5.1f} {1 if f==6 else 0:>4} "
          f"{pred_zero[f]:>+8.4f} {delta[f]:>+8.4f} {pred_zero[f]-delta[f]:>+8.4f} {pct_m:>+8.1f}%")

rms_zero = np.sqrt(np.mean((pred_zero - delta)**2))
rms_zero_q = np.sqrt(np.mean([(pred_zero[f]-delta[f])**2 for f in [3,4,5,6,7,8]]))
rms_mass_zero = np.sqrt(np.mean([(100*(m_e*PHI**(n_bare[f]+pred_zero[f]) - pdg_masses[f])/pdg_masses[f])**2 for f in range(9)]))

print(f"\nRMS (delta, all 9): {rms_zero:.4f}")
print(f"RMS (delta, quarks only): {rms_zero_q:.4f}")
print(f"RMS (mass %, all 9): {rms_mass_zero:.1f}%")

# Compare with the bare formula (no corrections)
rms_bare = np.sqrt(np.mean([(100*(m_e*PHI**n_bare[f] - pdg_masses[f])/pdg_masses[f])**2 for f in range(9)]))
print(f"\nFor comparison:")
print(f"  Bare formula RMS: {rms_bare:.1f}%")
print(f"  Zero-param model RMS: {rms_mass_zero:.1f}%")
print(f"  Improvement: {rms_bare:.1f}% -> {rms_mass_zero:.1f}%")

# ================================================================
# Test: delta = (g/2)*T3 * Nc  (with Nc factor)
# ================================================================
print(f"\n{'='*72}")
print("ALTERNATIVE: delta = (g/2)*T3*Nc - (2/a1)*1_{b=0}")
print("(Same but with Nc factor, so leptons always get 0)")
print("="*72)

Nc = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1])
pred_zeroNc = np.zeros(9)
for f in range(9):
    g = gens_all[f]
    T3 = fermion_T3[f]
    is_b0 = 1 if (f == 6) else 0
    pred_zeroNc[f] = (g / 2) * T3 * Nc[f] - (2 / A1) * is_b0

print(f"\n{'Name':>5} {'pred':>8} {'actual':>8} {'diff':>8} {'%m_err':>8}")
print("-"*50)
for f in range(9):
    n_pred = n_bare[f] + pred_zeroNc[f]
    m_pred = m_e * PHI**n_pred
    pct_m = 100 * (m_pred - pdg_masses[f]) / pdg_masses[f]
    print(f"{fermion_names[f]:>5} "
          f"{pred_zeroNc[f]:>+8.4f} {delta[f]:>+8.4f} {pred_zeroNc[f]-delta[f]:>+8.4f} {pct_m:>+8.1f}%")

rms_mass_zeroNc = np.sqrt(np.mean([(100*(m_e*PHI**(n_bare[f]+pred_zeroNc[f]) - pdg_masses[f])/pdg_masses[f])**2 for f in range(9)]))
print(f"\nRMS (mass %): {rms_mass_zeroNc:.1f}%")

# ================================================================
# DEEPER: what does 2/a1 MEAN for the down quark?
# ================================================================
print(f"\n{'='*72}")
print("INTERPRETATION: What does delta_d = -2/a1 mean?")
print("="*72)

print(f"\nn_eff(d) = a1 - 2/a1 = {A1} - 2/{A1} = {A1 - 2/A1:.3f}")
print(f"         = (a1^2 - 2)/a1 = ({A1**2} - 2)/{A1} = {(A1**2-2)/A1:.3f}")
print(f"Actual n_eff(d) = {n_eff[6]:.4f}")
print(f"Match: {abs((A1-2/A1) - n_eff[6]) / n_eff[6] * 100:.2f}%")

print(f"\na1^2 - 2 = {A1**2 - 2}")
print(f"  = 23 (a prime number)")
print(f"  NOT a Fibonacci or Lucas number")
print(f"  23 = a1^2 - 2 = 25 - 2")

print(f"\nAlternative: 2 = dim(rho_1) in 2I = |Gal(Q(sqrt(5))/Q)|")
print(f"  n_eff(d) = a1 - dim(rho_1)/a1 = a1(1 - 1/(a1*dim_1))")
print(f"  = a1 * (a1-1)! / ... no, not clean")

# The Galois interpretation
print(f"\nGalois interpretation:")
print(f"  For b=0 fermions (rational sector), z = a (integer)")
print(f"  The Galois automorphism sigma: phi -> phi' is TRIVIAL on z")
print(f"  But the mass m = m_e * phi^n involves IRRATIONAL phi^n")
print(f"  Norm: N(phi^5) = phi^5 * phi'^5 = (-1)^5 = -1")
print(f"  So phi^5 has NEGATIVE norm in Z[phi]!")
print(f"  This means it's NOT a totally positive element")
print(f"")
print(f"  The 'correction' 2/a1 could be:")
print(f"  - The GALOIS DEFECT: ratio |Gal|/a1 for the rational boundary")
print(f"  - Or: 1/a1 per Galois automorphism (2 autos in Gal(Q(sqrt5)/Q))")
print(f"  - Or: simply n_d=a1 -> correction ~ 1/n_d^(1/2) ? No, 1/sqrt(5)=0.447")

# Check: is the EFFECTIVE exponent n=23/5 special in the framework?
print(f"\nIs 23/5 special?")
print(f"  23 = dim of Leech lattice complement? No (Leech=24)")
print(f"  23/5 = 4.6")
print(f"  phi^(23/5) = {PHI**(23/5):.4f}")
print(f"  m_e * phi^(23/5) = {m_e * PHI**(23/5):.4f} MeV")
print(f"  PDG m_d = 4.67 MeV")
print(f"  Error: {100*(m_e*PHI**(23/5) - m_d)/m_d:+.2f}%")

# Is there a BETTER rational approximation?
print(f"\nRational approximations to n_eff(d) = {n_eff[6]:.6f}:")
for num in range(1, 50):
    for den in range(1, 20):
        val = num / den
        if abs(val - n_eff[6]) < 0.01:
            err = abs(val - n_eff[6])
            print(f"  {num}/{den} = {val:.4f}, diff = {err:.6f}")

# ================================================================
# SCAN: All expressions a*p + b*q for framework atoms p,q and small a,b
# ================================================================
print(f"\n{'='*72}")
print("SCAN: Best 2-atom expressions for n_eff(d) = {:.6f}".format(n_eff[6]))
print("="*72)

target = n_eff[6]
atoms = {
    'a1': A1, 'b1': B1, 'phi': PHI, '1/phi': 1/PHI, 'phi^2': PHI**2,
    '1/phi^2': 1/PHI**2, 'alpha_s': ALPHA_S, 'sin2tW': SIN2TW,
    'N_gen': N_GEN, 'sqrt(5)': np.sqrt(5), 'pi': np.pi,
    '1': 1, '2': 2, '3': 3, '4': 4, '6': 6, '1/2': 0.5,
    'phi^3': PHI**3, '1/phi^3': 1/PHI**3, 'alpha': ALPHA,
}

best = []
for name1, v1 in atoms.items():
    for name2, v2 in atoms.items():
        for a in [-3, -2, -1, 0, 1, 2, 3]:
            for b in [-3, -2, -1, 0, 1, 2, 3]:
                if a == 0 and b == 0:
                    continue
                val = a * v1 + b * v2
                if val > 0:
                    err = abs(val - target) / target
                    if err < 0.005:  # within 0.5%
                        expr = f"{a}*{name1} + {b}*{name2}" if b >= 0 else f"{a}*{name1} - {abs(b)}*{name2}"
                        best.append((err, expr, val))

best.sort()
seen = set()
for err, expr, val in best[:20]:
    if val not in seen:
        seen.add(round(val, 8))
        print(f"  {expr:40s} = {val:.6f}  ({100*err:.3f}%)")

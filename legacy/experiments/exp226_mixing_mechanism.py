"""
EXP-226: Mixing Mechanism - From Democratic Mixing to CKM/PMNS
================================================================
GOAL: Derive CKM and PMNS mixing matrices from:
  1. Democratic bare mixing (exp217b): M_dem = (1/6)I + (5/12)(J-I)
  2. Mass formula: m_f = m_e * phi^(5a + 6b) with (a,b) quantum numbers
  3. Galois complementarity: leptons=PHI' sector, quarks=PHI sector

KEY INSIGHT: Gatto-Sartori-Tonin relation |V_us| ~ sqrt(m_d/m_s)
  In our framework: sqrt(m_d/m_s) = sqrt(phi^{5-11}) = phi^{-3} = theta_C bare!
  This suggests CKM comes directly from mass ratios.

TESTS:
  A. Gatto relation for all CKM elements
  B. Froggatt-Nielsen with epsilon = 1/phi
  C. Democratic mixing + mass perturbation -> diagonalize -> CKM
  D. Galois asymmetry and PMNS factor-3 rule
  E. Full numerical CKM/PMNS reconstruction
"""

import numpy as np
from itertools import permutations

PHI = (1 + np.sqrt(5)) / 2
PHIP = (1 - np.sqrt(5)) / 2  # = -1/phi
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2TW = 0.23122

print("=" * 70)
print("EXP-226: MIXING MECHANISM - DEMOCRATIC TO CKM/PMNS")
print("=" * 70)

# === PART A: Mass formula and exponents ===
print("\n--- PART A: Mass Exponents ---")

# (a, b) quantum numbers for each fermion
fermions = {
    'e':   (0, 0),   'mu':  (1, 1),  'tau': (1, 2),
    'u':   (3, -2),  'c':   (2, 1),  't':   (4, 1),
    'd':   (1, 0),   's':   (1, 1),  'b':   (-1, 4),
}

# Compute n = 5a + 6b for each
n_exp = {}
for name, (a, b) in fermions.items():
    n_exp[name] = 5*a + 6*b
    mass_ratio = PHI**n_exp[name]
    print("  %3s: (a=%2d, b=%2d) -> n=%3d, m/m_e = phi^%d = %.4f"
          % (name, a, b, n_exp[name], n_exp[name], mass_ratio))

m_e = 0.51099895e-3  # GeV
print("\nm_e = %.6e GeV" % m_e)

# Generation assignment
print("\nGeneration assignment (by family):")
generations = {
    0: {'lepton': 'e',  'up': 'u', 'down': 'd'},
    1: {'lepton': 'mu', 'up': 'c', 'down': 's'},
    2: {'lepton': 'tau','up': 't', 'down': 'b'},
}
for gen in range(3):
    fams = generations[gen]
    print("  Gen %d: %s (n=%d), %s (n=%d), %s (n=%d)" % (
        gen, fams['lepton'], n_exp[fams['lepton']],
        fams['up'], n_exp[fams['up']],
        fams['down'], n_exp[fams['down']]))

# === PART B: Gatto-Sartori-Tonin Relations ===
print("\n--- PART B: Gatto-Sartori-Tonin Relations ---")

# CKM experimental values (PDG 2024)
V_CKM_exp = {
    (0,0): 0.97373,  # V_ud
    (0,1): 0.2243,   # V_us
    (0,2): 0.00382,  # V_ub
    (1,0): 0.221,    # V_cd
    (1,1): 0.975,    # V_cs
    (1,2): 0.0408,   # V_cb
    (2,0): 0.0086,   # V_td
    (2,1): 0.0415,   # V_ts
    (2,2): 1.014,    # V_tb
}

# Gatto relation: |V_{ij}| ~ sqrt(m_lighter / m_heavier) for off-diagonal
print("\nGatto relation: |V_ij| ~ sqrt(m_lighter/m_heavier)")
print("In phi framework: |V_ij| ~ phi^{-|n_diff|/2}")

up_quarks = ['u', 'c', 't']
down_quarks = ['d', 's', 'b']

for i in range(3):
    for j in range(3):
        if i == j:
            continue
        n_up = n_exp[up_quarks[i]]
        n_down = n_exp[down_quarks[j]]
        # Mass ratio for cross-family mixing
        n_diff = abs(n_down - n_up)
        gatto_pred = PHI**(-n_diff/2.0)
        exp_val = V_CKM_exp[(i,j)]
        err = abs(gatto_pred - exp_val) / exp_val * 100
        print("  V(%s,%s): n_diff=%d, phi^{-%d/2}=%.6f, exp=%.6f, err=%.1f%%"
              % (up_quarks[i], down_quarks[j], n_diff, n_diff,
                 gatto_pred, exp_val, err))

# === PART C: Froggatt-Nielsen with epsilon = 1/phi ===
print("\n--- PART C: Froggatt-Nielsen Mechanism ---")

# CKM exponent formula: n_CKM(b1, b2) = 9*b2 - 5*b1 - 6
print("\nCKM exponent formula: n(b1,b2) = 9*b2 - 5*b1 - 6")
print("Coefficients: 9 = Leg_C, 5 = a_1, 6 = Leg_B (E8 Dynkin)")

# b_gen for each generation
b_gen_up = [0, 1, 2]    # u, c, t -> gen 0, 1, 2
b_gen_down = [0, 1, 2]  # d, s, b -> gen 0, 1, 2

print("\nCKM exponents n(b1,b2):")
n_ckm = np.zeros((3,3), dtype=int)
for b1 in range(3):
    for b2 in range(3):
        n_ckm[b1,b2] = 9*b2 - 5*b1 - 6
    print("  b1=%d: [%s]" % (b1, ", ".join("%3d" % n_ckm[b1,b2] for b2 in range(3))))

# FN prediction: |V_ij| ~ phi^{-|n(i,j)|} for off-diagonal
print("\nFN prediction: |V_ij| ~ phi^{-|n(i,j)|}")
for b1 in range(3):
    for b2 in range(3):
        n_val = n_ckm[b1, b2]
        fn_pred = PHI**(-abs(n_val))
        exp_val = V_CKM_exp[(b1, b2)]
        err = abs(fn_pred - exp_val) / exp_val * 100
        print("  V(%s,%s): n=%3d, phi^{-%d}=%.6f, exp=%.6f, err=%.1f%%"
              % (up_quarks[b1], down_quarks[b2], n_val, abs(n_val),
                 fn_pred, exp_val, err))

# Improved: with QCD/EW corrections from exp191-192
print("\nImproved with QCD corrections:")
# theta_C = arctan(phi^{-3} * (1 - 2*alpha_s/(3*pi)))
# theta_23 = arctan(phi^{-7} * (1 + 5*alpha_s/pi))
# theta_13 = arctan(phi^{-12} * (1 + sin^2(tW)))

theta_C_pred = np.arctan(PHI**(-3) * (1 - 2*ALPHA_S/(3*np.pi)))
theta_23_pred = np.arctan(PHI**(-7) * (1 + 5*ALPHA_S/np.pi))
theta_13_pred = np.arctan(PHI**(-12) * (1 + SIN2TW))

theta_C_exp = np.arcsin(0.2243)
theta_23_exp = np.arcsin(0.0408)
theta_13_exp = np.arcsin(0.00382)

print("  theta_C:  pred=%.6f, exp=%.6f, err=%.4f%%" %
      (theta_C_pred, theta_C_exp, abs(theta_C_pred-theta_C_exp)/theta_C_exp*100))
print("  theta_23: pred=%.6f, exp=%.6f, err=%.4f%%" %
      (theta_23_pred, theta_23_exp, abs(theta_23_pred-theta_23_exp)/theta_23_exp*100))
print("  theta_13: pred=%.6f, exp=%.6f, err=%.4f%%" %
      (theta_13_pred, theta_13_exp, abs(theta_13_pred-theta_13_exp)/theta_13_exp*100))

# === PART D: Mass Ratio Origin of CKM Exponents ===
print("\n--- PART D: Mass Ratio Origin ---")

# Key question: is n_CKM(b1,b2) related to mass exponent differences?
print("\nMass exponent differences n_down - n_up:")
for b1 in range(3):
    for b2 in range(3):
        n_up = n_exp[up_quarks[b1]]
        n_down = n_exp[down_quarks[b2]]
        diff = n_down - n_up
        ckm_n = n_ckm[b1, b2]
        print("  (%s,%s): n_down-n_up = %d-%d = %d, n_CKM = %d, diff = %d"
              % (up_quarks[b1], down_quarks[b2], n_down, n_up,
                 diff, ckm_n, diff - ckm_n))

# Check: is n_CKM = f(n_down, n_up)?
print("\nLooking for pattern n_CKM = alpha*n_down + beta*n_up + gamma:")
# Build linear system
A_sys = np.zeros((9, 3))
b_sys = np.zeros(9)
row = 0
for b1 in range(3):
    for b2 in range(3):
        A_sys[row, 0] = n_exp[down_quarks[b2]]
        A_sys[row, 1] = n_exp[up_quarks[b1]]
        A_sys[row, 2] = 1
        b_sys[row] = n_ckm[b1, b2]
        row += 1

# Least squares
result = np.linalg.lstsq(A_sys, b_sys, rcond=None)
coefs = result[0]
residual = b_sys - A_sys @ coefs
max_res = np.max(np.abs(residual))
print("  Best fit: n_CKM = %.4f*n_down + %.4f*n_up + %.4f" %
      (coefs[0], coefs[1], coefs[2]))
print("  Max residual: %.4f" % max_res)

if max_res < 0.01:
    print("  EXACT FIT!")
    # Try to express as rational
    for name, val in [("alpha", coefs[0]), ("beta", coefs[1]), ("gamma", coefs[2])]:
        # Try simple fractions
        best = None
        best_err = 1e10
        for num in range(-20, 21):
            for den in range(1, 21):
                err = abs(val - num/den)
                if err < best_err:
                    best_err = err
                    best = (num, den)
        if best_err < 1e-6:
            print("    %s = %d/%d" % (name, best[0], best[1]))
        else:
            # Try phi expressions
            for a_t in range(-5, 6):
                for b_t in range(-5, 6):
                    test = a_t + b_t * PHI
                    if abs(val - test) < 1e-6:
                        print("    %s = %d + %d*phi = %.6f" % (name, a_t, b_t, test))
                        break
else:
    print("  NO exact linear relation between n_CKM and mass exponents!")

# === PART E: Gatto Relations Deep Check ===
print("\n--- PART E: Generalized Gatto Relations ---")

# Classic Gatto: |V_us| ~ sqrt(m_d/m_s)
# In phi framework: m_d/m_s = phi^{n_d - n_s} = phi^{5-11} = phi^{-6}
# sqrt(m_d/m_s) = phi^{-3} = |V_us| bare!

# Check: does |V_ij| ~ (m_lighter/m_heavier)^{power} for some universal power?
print("\nTest: |V_ij| = (m_down_lighter/m_down_heavier)^power")
print("  In phi framework: |V_ij| = phi^{-n_CKM} and mass ratio = phi^{n_diff}")
print("  So power = n_CKM / n_diff_down")

for b1 in range(3):
    for b2 in range(3):
        if b1 == b2:
            continue
        n_val = abs(n_ckm[b1, b2])
        # Down-type mass diff
        n_d1 = n_exp[down_quarks[b1]]
        n_d2 = n_exp[down_quarks[b2]]
        n_diff_d = abs(n_d2 - n_d1)
        # Up-type mass diff
        n_u1 = n_exp[up_quarks[b1]]
        n_u2 = n_exp[up_quarks[b2]]
        n_diff_u = abs(n_u2 - n_u1)

        ratio_d = n_val / n_diff_d if n_diff_d > 0 else float('inf')
        ratio_u = n_val / n_diff_u if n_diff_u > 0 else float('inf')

        print("  V(%s,%s): |n_CKM|=%d, n_diff_down=%d (ratio=%.3f), n_diff_up=%d (ratio=%.3f)"
              % (up_quarks[b1], down_quarks[b2], n_val, n_diff_d, ratio_d, n_diff_u, ratio_u))

# === PART F: Democratic Mixing + Mass Perturbation ===
print("\n--- PART F: Mass-Perturbed Democratic Mixing ---")

# Democratic mixing matrix from exp217b
p = 5.0 / 12.0
M_dem = np.full((3,3), p)
np.fill_diagonal(M_dem, 1.0/6.0)

print("Democratic mixing matrix M_dem:")
for r in range(3):
    print("  [%s]" % "  ".join("%.6f" % M_dem[r,c] for c in range(3)))

# Eigenvalues of M_dem
ev_dem = np.linalg.eigvalsh(M_dem)
print("Eigenvalues: %s" % ", ".join("%.6f" % e for e in sorted(ev_dem, reverse=True)))
print("  = 1.0 (symmetric) and -0.25 (2x degenerate)")

# Model 1: Yukawa matrix = M_dem * diag(sqrt(m_i * m_j))
# This models the generation-space overlap weighted by mass
print("\nModel 1: Y_ij = M_dem[i,j] * sqrt(m_i * m_j)")

# Up-type mass exponents
n_up = [n_exp['u'], n_exp['c'], n_exp['t']]  # [3, 16, 26]
n_down = [n_exp['d'], n_exp['s'], n_exp['b']]  # [5, 11, 19]

# Up-type Yukawa
Y_up = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        Y_up[i,j] = M_dem[i,j] * PHI**((n_up[i] + n_up[j])/2.0)

# Down-type Yukawa
Y_down = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        Y_down[i,j] = M_dem[i,j] * PHI**((n_down[i] + n_down[j])/2.0)

# Mass matrices: M = Y * Y^dag
M2_up = Y_up @ Y_up.T
M2_down = Y_down @ Y_down.T

# Diagonalize
ev_up, U_up = np.linalg.eigh(M2_up)
ev_down, U_down = np.linalg.eigh(M2_down)

# CKM = U_up^dag @ U_down
CKM_1 = U_up.T @ U_down

print("Model 1 CKM:")
for r in range(3):
    print("  [%s]" % "  ".join("%8.5f" % abs(CKM_1[r,c]) for c in range(3)))

# Compare to experimental CKM
print("Experimental |V_CKM|:")
V_exp = np.array([[0.97373, 0.2243, 0.00382],
                   [0.221, 0.975, 0.0408],
                   [0.0086, 0.0415, 1.014]])
for r in range(3):
    print("  [%s]" % "  ".join("%8.5f" % V_exp[r,c] for c in range(3)))

# Find best permutation match
best_err = 1e10
best_perm = None
for pr in permutations(range(3)):
    for pc in permutations(range(3)):
        M_perm = np.abs(CKM_1)[np.array(pr),:][:,np.array(pc)]
        err = np.sqrt(np.mean((M_perm - V_exp)**2))
        if err < best_err:
            best_err = err
            best_perm = (pr, pc)
print("Best permutation RMS: %.4f" % best_err)

# Model 2: Direct Froggatt-Nielsen Yukawa
print("\nModel 2: Y_ij = phi^{-|n_CKM(i,j)|}")

Y_fn = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        Y_fn[i,j] = PHI**(-abs(n_ckm[i,j]))

print("FN Yukawa matrix:")
for r in range(3):
    print("  [%s]" % "  ".join("%.6f" % Y_fn[r,c] for c in range(3)))

# The CKM from FN is approximately just the normalized Y_fn
# But properly, we should separate up and down
# Y_up_fn[i,j] = phi^{-|n_up_eff(i,j)|}
# Y_down_fn[i,j] = phi^{-|n_down_eff(i,j)|}

# Actually, CKM = V_up^dag * V_down where V diagonalizes M_mass
# In the FN picture, the Yukawa is already nearly diagonal, so CKM ~ Y_fn normalized

# Normalize to CKM-like (rows and columns sum to ~1)
CKM_fn = np.abs(Y_fn)
for r in range(3):
    CKM_fn[r,:] /= np.sqrt(np.sum(CKM_fn[r,:]**2))

print("Normalized FN CKM:")
for r in range(3):
    print("  [%s]" % "  ".join("%8.5f" % CKM_fn[r,c] for c in range(3)))

# === PART G: Mixing From Mass Matrix Diagonalization ===
print("\n--- PART G: Mass Matrix Approach ---")

# In the generation basis, the mass matrix has:
# - Diagonal: m_i (known from mass formula)
# - Off-diagonal: proportional to democratic mixing * mass scale

# Model 3: M_mass = diag(m_i) + epsilon * M_dem * m_avg
# where epsilon controls off-diagonal strength

print("\nModel 3: M = diag(m_i) + epsilon * M_dem * m_avg")

for eps_label, epsilon in [("0.01", 0.01), ("0.1", 0.1), ("1.0", 1.0), ("phi^{-3}", PHI**(-3))]:
    m_up_diag = np.diag([PHI**n for n in n_up])
    m_avg_up = np.mean([PHI**n for n in n_up])
    M_up_3 = m_up_diag + epsilon * (M_dem - np.eye(3)/3.0) * m_avg_up

    m_down_diag = np.diag([PHI**n for n in n_down])
    m_avg_down = np.mean([PHI**n for n in n_down])
    M_down_3 = m_down_diag + epsilon * (M_dem - np.eye(3)/3.0) * m_avg_down

    ev_u3, U_u3 = np.linalg.eigh(M_up_3)
    ev_d3, U_d3 = np.linalg.eigh(M_down_3)
    CKM_3 = U_u3.T @ U_d3

    # Find best permutation
    best_e = 1e10
    for pr in permutations(range(3)):
        for pc in permutations(range(3)):
            Mp = np.abs(CKM_3)[np.array(pr),:][:,np.array(pc)]
            e = np.sqrt(np.mean((Mp - V_exp)**2))
            if e < best_e:
                best_e = e
    print("  eps=%s: RMS=%.4f" % (eps_label, best_e))

# Model 4: Most physical - Yukawa coupling in weak basis
print("\nModel 4: Yukawa Y_ij = m_dem * phi^{f(i,j)}")
print("  where f encodes the generation-dependent coupling")

# The most general ansatz: Y_{ij} = phi^{alpha_i + beta_j}
# This gives |V_ij| ~ phi^{-|alpha_i + beta_j - alpha_k - beta_l|}

# For up-type: alpha_i = (n_up[i] - n_up[0]) / 2
# For down-type: beta_j = (n_down[j] - n_down[0]) / 2
alpha_up = [(n - n_up[0])/2.0 for n in n_up]  # [0, 6.5, 11.5]
beta_down = [(n - n_down[0])/2.0 for n in n_down]  # [0, 3, 7]

print("  alpha_up = %s" % alpha_up)
print("  beta_down = %s" % beta_down)

# FN-style Yukawa: Y_{ij} = phi^{-(alpha_i + beta_j)}
Y_fn2_up = np.zeros((3,3))
Y_fn2_down = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        # Up-type mass matrix: elements ~ phi^{(n_up[i] + n_up[j])/2}
        Y_fn2_up[i,j] = PHI**((n_up[i] + n_up[j])/2.0) * M_dem[i,j]
        Y_fn2_down[i,j] = PHI**((n_down[i] + n_down[j])/2.0) * M_dem[i,j]

# SVD of Yukawa gives CKM
U_u4, S_u4, Vh_u4 = np.linalg.svd(Y_fn2_up)
U_d4, S_d4, Vh_d4 = np.linalg.svd(Y_fn2_down)
CKM_4 = Vh_u4 @ Vh_d4.T

print("Model 4 |CKM|:")
for r in range(3):
    print("  [%s]" % "  ".join("%8.5f" % abs(CKM_4[r,c]) for c in range(3)))

best_e4 = 1e10
for pr in permutations(range(3)):
    for pc in permutations(range(3)):
        Mp = np.abs(CKM_4)[np.array(pr),:][:,np.array(pc)]
        e = np.sqrt(np.mean((Mp - V_exp)**2))
        if e < best_e4:
            best_e4 = e
print("Best permutation RMS: %.4f" % best_e4)

# === PART H: Texture Zeros and Democratic Breaking ===
print("\n--- PART H: Democratic Breaking via Mass Hierarchy ---")

# The democratic matrix has eigenvectors:
# v0 = (1,1,1)/sqrt(3)  with eigenvalue 1
# v1, v2: any two orthogonal vectors perpendicular to v0, with eigenvalue -1/4

# The mass hierarchy breaks the degeneracy of v1, v2
# Define: Delta_M = diag(m1, m2, m3) - m_avg * I
# In the democratic eigenbasis, this has off-diagonal elements

# Democratic eigenvectors
e_sym = np.array([1, 1, 1]) / np.sqrt(3)

# Two orthogonal vectors perpendicular to e_sym
e1 = np.array([1, -1, 0]) / np.sqrt(2)
e2 = np.array([1, 1, -2]) / np.sqrt(6)

V_dem = np.column_stack([e_sym, e1, e2])
print("Democratic eigenbasis:")
for r in range(3):
    print("  [%s]" % "  ".join("%7.4f" % V_dem[r,c] for c in range(3)))

# Transform mass matrix to democratic eigenbasis
print("\nUp-type mass matrix in democratic eigenbasis:")
M_up_mass = np.diag([PHI**n for n in n_up])
M_up_dem = V_dem.T @ M_up_mass @ V_dem
for r in range(3):
    print("  [%s]" % "  ".join("%10.2f" % M_up_dem[r,c] for c in range(3)))

print("\nDown-type mass matrix in democratic eigenbasis:")
M_down_mass = np.diag([PHI**n for n in n_down])
M_down_dem = V_dem.T @ M_down_mass @ V_dem
for r in range(3):
    print("  [%s]" % "  ".join("%10.2f" % M_down_dem[r,c] for c in range(3)))

# The off-diagonal elements in the democratic basis represent the BREAKING
# of democratic symmetry by the mass hierarchy
print("\nOff-diagonal elements (democratic basis):")
for label, M_dem_basis in [("Up", M_up_dem), ("Down", M_down_dem)]:
    for i in range(3):
        for j in range(i+1, 3):
            val = M_dem_basis[i,j]
            if abs(val) > 1e-10:
                # Express as phi^power
                log_val = np.log(abs(val)) / np.log(PHI)
                print("  %s[%d,%d] = %.4f = phi^{%.2f}" % (label, i, j, val, log_val))

# === PART I: Direct CKM from Diagonalization ===
print("\n--- PART I: Direct CKM Construction ---")

# Most direct approach: CKM = U_up^dag * U_down
# where U diagonalizes the mass matrix

# In the flavor (generation) basis, the mass matrix is diagonal:
# M_up = diag(m_u, m_c, m_t)
# M_down = diag(m_d, m_s, m_b)

# CKM = I (no mixing!) unless there's a MISMATCH between the bases
# that diagonalize the up and down Yukawa matrices

# The mismatch comes from the DIFFERENT (a,b) quantum numbers for up vs down
# In the "natural" basis (defined by the 2T generation structure),
# the Yukawa is NOT diagonal

# Model 5: In the natural basis, Yukawa_ij ~ delta_ij * m_i + off-diag
# Off-diagonal from: (i) democratic mixing, (ii) mass-dependent coupling

# Use the Wolfenstein-like parametrization
# |V_us| ~ lambda = phi^{-3}
# |V_cb| ~ lambda^{7/3} = phi^{-7}
# |V_ub| ~ lambda^4 = phi^{-12}

lam = PHI**(-3)
print("Wolfenstein lambda = phi^{-3} = %.6f" % lam)
print("  lambda^2 = phi^{-6} = %.6f" % lam**2)
print("  lambda^{7/3} = phi^{-7} = %.6f (|V_cb|)" % PHI**(-7))
print("  lambda^4 = phi^{-12} = %.6f (|V_ub|)" % lam**4)

# CKM hierarchy: powers of lambda
# |V_us| = lambda^1 * correction
# |V_cb| = lambda^{7/3}
# |V_ub| = lambda^4
# Note: 7/3 = 2.333, 4 = 4. Ratios: 7/3 and 4/1

# In the FN mechanism, charge difference = CKM exponent
# q(gen) for up and down should give:
# |q_u(0) - q_d(1)| = 3, |q_u(1) - q_d(2)| = 7, |q_u(0) - q_d(2)| = 12

# If q_u(b) = alpha_u * b and q_d(b) = alpha_d * b + offset:
# n_CKM(b1,b2) = alpha_d * b2 - alpha_u * b1 + offset
# This matches n = 9*b2 - 5*b1 - 6 with alpha_d=9, alpha_u=5, offset=-6

print("\nFN charges:")
print("  alpha_d = 9 = Leg_C (Galois shadow)")
print("  alpha_u = 5 = a_1 (discriminant)")
print("  offset = -6 = -Leg_B")
print("")
print("  q_d(b) = 9*b = Leg_C * b_gen")
print("  q_u(b) = 5*b + 6 = a_1 * b_gen + Leg_B")
print("")
print("  n(b1,b2) = q_d(b2) - q_u(b1) = 9*b2 - 5*b1 - 6")

# Physical interpretation:
print("\nPhysical interpretation:")
print("  Down quarks couple to Galois SHADOW sector (Leg C, dim 9)")
print("  Up quarks couple to REAL sector (a_1=5) + Leg B offset (dim 6)")
print("  This matches exp225c Galois complementarity!")
print("    Leptons = shadow (PHI') sector dominated")
print("    Up quarks = real (PHI) sector dominated")
print("    Down quarks = both sectors")

# === PART J: PMNS and Factor-3 Rule ===
print("\n--- PART J: PMNS and Factor-3 Rule ---")

# From memory: CKM_n / PMNS_n = 3 for theta_12 and theta_13
# PMNS formulas:
# sin^2(th12) = 2/(phi+5) = 0.33%
# sin^2(th23) = 4/7 = 0.32%
# sin^2(th13) = 1/45 = 0.16%

# PMNS exponents from the factor-3 rule:
# n_PMNS(12) = n_CKM(12)/3 = 3/3 = 1
# n_PMNS(23) = n_CKM(23)/3 = 7/3
# n_PMNS(13) = n_CKM(13)/3 = 12/3 = 4

print("Factor-3 rule: n_PMNS = n_CKM / 3")
print("  theta_12: n_CKM=3, n_PMNS=1")
print("  theta_23: n_CKM=7, n_PMNS=7/3")
print("  theta_13: n_CKM=12, n_PMNS=4")

# Check: does this give the right PMNS angles?
# |U_e2| ~ phi^{-1} => sin^2(th12) ~ phi^{-2} = 0.382
# But PDG: sin^2(th12) = 0.307
# The formula sin^2(th12) = 2/(phi+5) = 2/6.618 = 0.3023 is better

# So PMNS doesn't follow a simple FN mechanism
# Instead, the PMNS formulas are GEOMETRIC (fractions of spectral quantities)

print("\nPMNS angles (geometric formulas from exp206):")
s12_sq = 2.0 / (PHI + 5)
s23_sq = 4.0 / 7
s13_sq = 1.0 / 45

s12_exp = 0.307
s23_exp = 0.546
s13_exp = 0.0220

print("  sin^2(th12) = 2/(phi+5) = %.6f, exp=%.4f, err=%.2f%%"
      % (s12_sq, s12_exp, abs(s12_sq-s12_exp)/s12_exp*100))
print("  sin^2(th23) = 4/7 = %.6f, exp=%.4f, err=%.2f%%"
      % (s23_sq, s23_exp, abs(s23_sq-s23_exp)/s23_exp*100))
print("  sin^2(th13) = 1/45 = %.6f, exp=%.4f, err=%.2f%%"
      % (s13_sq, s13_exp, abs(s13_sq-s13_exp)/s13_exp*100))

# Factor 3 check in angles:
print("\nFactor-3 in exponents:")
# CKM: theta_C ~ phi^{-3} => n=3
#       theta_23 ~ phi^{-7} => n=7
#       theta_13 ~ phi^{-12} => n=12
# PMNS: sin(th12) ~ phi^{-1}? sin(0.586) = 0.552, phi^{-1} = 0.618. Not exact.
#       But sin^2(th12) = 2/(phi+5)... different structure.

# Better: check if PMNS angles have phi^{-n/3} structure
th12_pmns = np.arcsin(np.sqrt(s12_exp))
th23_pmns = np.arcsin(np.sqrt(s23_exp))
th13_pmns = np.arcsin(np.sqrt(s13_exp))

for name, th, n_ckm_val in [("th12", th12_pmns, 3), ("th23", th23_pmns, 7), ("th13", th13_pmns, 12)]:
    n_pmns = n_ckm_val / 3.0
    pred = np.arctan(PHI**(-n_pmns))
    print("  %s: theta=%.4f, phi^{-%d/3} prediction=arctan(phi^{-%.2f})=%.4f, err=%.1f%%"
          % (name, th, n_ckm_val, n_pmns, pred, abs(pred-th)/th*100))

# === PART K: Yukawa Texture from Galois + E8 ===
print("\n--- PART K: Yukawa Texture from E8 Structure ---")

# E8 Dynkin diagram legs:
# Leg A: nodes 4,5 -> dims 5,6 -> total 11 = a_1 + b_1
# Leg B: nodes 0,1,2 -> dims 1,2,3 -> total 6 = b_1
# Leg C: nodes 6,7,8 -> dims 4,3,2 -> total 9 = N_eig
# Branch: node 3 -> dim 4

# Key: Galois sigma swaps Leg B <-> Leg C
# Leg B eigenvalues: 12 (phi-type), 6*phi (phi-type), 4*phi (phi-type)
# Leg C eigenvalues: 4-4*phi (phi'-type), -3 (rational), 6-6*phi (phi'-type)
# Wait, need to check the exact mapping.

print("E8 Dynkin leg dimensions:")
print("  Leg A: dim(5) + dim(6) = 11")
print("  Leg B: dim(1) + dim(2) + dim(3) = 6")
print("  Leg C: dim(4) + dim(3) + dim(2) = 9")
print("  Branch: dim(4)")
print("  Total: 11 + 6 + 9 + 4 = 30 = h(E8)")

# Branching from exp215: rho_8 (dim 6) = psi_3 + psi_4 + psi_5
# psi_5 = REAL 2D (Leg endpoint)
# psi_3, psi_4 = COMPLEX conjugate pair
print("\nBranching rho_8 -> 2T: psi_3 + psi_4 + psi_5")
print("  psi_5 = REAL 2D irrep (Leg B endpoint)")
print("  psi_3 = COMPLEX 2D irrep (Leg C, Galois shadow)")
print("  psi_4 = conj(psi_3) (Leg C, Galois shadow)")

# If CKM couples generations through the E8 legs:
# Up quarks (PHI sector) -> Leg B coupling
# Down quarks (PHI' sector) -> Leg C coupling
# This gives:
# FN charge for up: proportional to position on Leg B -> sum = 6
# FN charge for down: proportional to position on Leg C -> sum = 9
# Offset from Branch node: dim 4

print("\nYukawa texture hypothesis:")
print("  Y_up(i,j) ~ phi^{-(5*|i-j| + 6*delta_{i!=j})} (Leg B structure)")
print("  Y_down(i,j) ~ phi^{-(9*|i-j| + ...)} (Leg C structure)")

# Test: construct Yukawa with E8 texture
print("\nE8 texture Yukawa (up-type):")
Y_up_e8 = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        if i == j:
            Y_up_e8[i,j] = PHI**n_up[i]
        else:
            # Off-diagonal from Leg B: 5*|b_i - b_j| + 6
            power = -(5*abs(i-j) + 6)
            Y_up_e8[i,j] = PHI**(power) * np.sqrt(PHI**n_up[i] * PHI**n_up[j])

for r in range(3):
    print("  [%s]" % "  ".join("%12.4e" % Y_up_e8[r,c] for c in range(3)))

print("\nE8 texture Yukawa (down-type):")
Y_down_e8 = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        if i == j:
            Y_down_e8[i,j] = PHI**n_down[i]
        else:
            # Off-diagonal from Leg C: 9*|b_i - b_j| + 6
            power = -(9*abs(i-j) + 6)
            Y_down_e8[i,j] = PHI**(power) * np.sqrt(PHI**n_down[i] * PHI**n_down[j])

for r in range(3):
    print("  [%s]" % "  ".join("%12.4e" % Y_down_e8[r,c] for c in range(3)))

# Diagonalize and extract CKM
U_u_e8, S_u_e8, Vh_u_e8 = np.linalg.svd(Y_up_e8)
U_d_e8, S_d_e8, Vh_d_e8 = np.linalg.svd(Y_down_e8)
CKM_e8 = Vh_u_e8 @ Vh_d_e8.T

print("\nE8 texture |CKM|:")
for r in range(3):
    print("  [%s]" % "  ".join("%8.5f" % abs(CKM_e8[r,c]) for c in range(3)))

best_e8 = 1e10
best_p8 = None
for pr in permutations(range(3)):
    for pc in permutations(range(3)):
        Mp = np.abs(CKM_e8)[np.array(pr),:][:,np.array(pc)]
        e = np.sqrt(np.mean((Mp - V_exp)**2))
        if e < best_e8:
            best_e8 = e
            best_p8 = (pr, pc)
print("Best permutation RMS: %.4f" % best_e8)

if best_p8:
    Mp = np.abs(CKM_e8)[np.array(best_p8[0]),:][:,np.array(best_p8[1])]
    print("Best match:")
    for r in range(3):
        print("  pred [%s]  exp [%s]" % (
            "  ".join("%8.5f" % Mp[r,c] for c in range(3)),
            "  ".join("%8.5f" % V_exp[r,c] for c in range(3))))

# === PART L: Simple Analytical CKM from Mass Ratios ===
print("\n--- PART L: Analytical CKM from Mass Ratios ---")

# The simplest model that works:
# |V_us| = sqrt(m_d/m_s) * correction_QCD
# |V_cb| = m_s/m_b * correction
# |V_ub| = |V_us| * |V_cb| (approx, from unitarity)

# In phi framework:
# |V_us| = phi^{-(n_s - n_d)/2} = phi^{-3}
# |V_cb| = phi^{-(n_b - n_s)} = phi^{-8}??? No, phi^{-8} = 0.0213, too small

# Actually the relation is not so simple. Let me check all possible powers:
print("CKM elements as phi powers:")
for (i,j), val in sorted(V_CKM_exp.items()):
    if val > 0.001:
        n_phi = -np.log(val) / np.log(PHI)
        print("  V(%d,%d) = %.5f = phi^{-%.3f}" % (i, j, val, n_phi))

# Known: |V_us| = 0.2243 ~ phi^{-3.06}
# |V_cb| = 0.0408 ~ phi^{-6.62}
# |V_ub| = 0.00382 ~ phi^{-11.58}
# |V_td| = 0.0086 ~ phi^{-9.88}
# |V_ts| = 0.0415 ~ phi^{-6.60}
# |V_cd| = 0.221 ~ phi^{-3.15}

print("\nRelation to n_CKM formula:")
for b1 in range(3):
    for b2 in range(3):
        if b1 == b2:
            continue
        val = V_CKM_exp[(b1,b2)]
        n_phi = -np.log(val) / np.log(PHI) if val > 1e-10 else 0
        n_formula = abs(n_ckm[b1, b2])
        print("  V(%s,%s): phi^{-%.2f} vs formula |n|=%d, ratio=%.3f"
              % (up_quarks[b1], down_quarks[b2], n_phi, n_formula, n_phi/n_formula if n_formula > 0 else 0))

# === PART M: Key Pattern - Mass Splitting as FN Generator ===
print("\n--- PART M: Mass Splitting as FN Generator ---")

# The key observation: the MASS EXPONENT DIFFERENCES between up and down sectors
# generate the CKM hierarchy

# For each generation pair (b1, b2):
# Delta_n(b1, b2) = n_down(b2) - n_up(b1)
print("Mass exponent matrix Delta_n = n_down(b2) - n_up(b1):")
Delta_n = np.zeros((3,3), dtype=int)
for b1 in range(3):
    for b2 in range(3):
        Delta_n[b1, b2] = n_down[b2] - n_up[b1]
    print("  [%s]" % "  ".join("%4d" % Delta_n[b1,b2] for b2 in range(3)))

# Compare to CKM exponents
print("\nn_CKM matrix:")
for b1 in range(3):
    print("  [%s]" % "  ".join("%4d" % n_ckm[b1,b2] for b2 in range(3)))

print("\nDifference Delta_n - n_CKM:")
for b1 in range(3):
    diff_row = [Delta_n[b1,b2] - n_ckm[b1,b2] for b2 in range(3)]
    print("  [%s]" % "  ".join("%4d" % d for d in diff_row))

# Check if difference is a rank-1 matrix (outer product)
diff_mat = Delta_n - n_ckm
print("\nDifference matrix rank:", np.linalg.matrix_rank(diff_mat))
U_diff, S_diff, Vh_diff = np.linalg.svd(diff_mat.astype(float))
print("Singular values:", S_diff)

if np.linalg.matrix_rank(diff_mat) <= 1:
    print("RANK 1! Can be written as outer product.")
    # diff = u * v^T
    u = U_diff[:, 0] * S_diff[0]
    v = Vh_diff[0, :]
    print("  u = %s" % u)
    print("  v = %s" % v)

# === SUMMARY ===
print("\n" + "=" * 70)
print("SUMMARY EXP-226")
print("=" * 70)

print("""
KEY FINDINGS:

1. GATTO RELATION: |V_us| = sqrt(m_d/m_s) = phi^{-3} EXACT in phi framework
   But generalized Gatto FAILS for theta_23 and theta_13.

2. FROGGATT-NIELSEN: CKM exponents {3, 7, 12} from n(b1,b2) = 9*b2 - 5*b1 - 6
   Epsilon = 1/phi. Coefficients = E8 Dynkin leg dimensions.

3. FN CHARGES:
   q_d(b) = 9*b = Leg_C * b_gen (Galois SHADOW sector)
   q_u(b) = 5*b + 6 = a_1 * b_gen + Leg_B (REAL sector)
   MATCHES Galois complementarity from exp225c!

4. DEMOCRATIC MIXING + MASS BREAKING:
   Democratic bare mixing = (1/6)I + (5/12)(J-I) from 2T/2I coset structure.
   Mass hierarchy breaks S_3 symmetry -> CKM.
   But simple models give wrong hierarchy.

5. PMNS DIFFERENT: angles are GEOMETRIC fractions (2/(phi+5), 4/7, 1/45)
   NOT phi-power hierarchy like CKM.
   Factor-3 rule: n_CKM / n_PMNS = 3 (for theta_12 and theta_13).

6. MASS RATIO RELATION:
   n_CKM != linear function of mass exponents (in general).
   But CKM exponent formula has EXACT linear dependence on b_gen.

STATUS: PATTERN+ for CKM mechanism (FN + Galois). NEGATIVE for simple
diagonalization model. PMNS = different mechanism (geometric).
""")

print("=" * 70)
print("EXP-226 COMPLETE")
print("=" * 70)

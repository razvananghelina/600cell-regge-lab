"""
exp333_cohn_kumar.py - Cohn-Kumar Energy and Physical Scale
============================================================
Direction F from alte_ddddd.txt: Self-consistency between gravitational
Riesz energy of 600-cell on S^3(R) and Casimir energy fixes R, hence m_e.

Also tests Direction D revisited: spectral action self-consistency -> Lambda.
Also tests Direction E: lattice theta function / norm analysis.

Author: Razvan-Constantin Anghelina
Date: 2026-02-16
"""

import numpy as np
from fractions import Fraction
from itertools import product as iprod
from math import gcd
import time

t0 = time.time()

print("=" * 72)
print("exp333: COHN-KUMAR ENERGY AND PHYSICAL SCALE")
print("Direction F + D + E from alte_ddddd.txt")
print("=" * 72)

# Framework constants
a1 = 5
phi = (1 + np.sqrt(5)) / 2
phi_sq = phi + 1  # phi^2
b1 = a1 + 1  # 6
N = 4 * a1 * b1  # 120 = |2I| = a_1!
degree = 2 * b1  # 12
N_eig = a1 + 2 * (a1 - 1)  # 9 (distinct eigenvalues)

# Physical constants
alpha_inv = 137.035999177
alpha = 1.0 / alpha_inv
alpha_s_fw = 1.0 / (2 * phi**3)
m_e_exp = 0.51099895000e-3  # GeV
m_P = 1.220890e19  # GeV
hbar_c = 0.19732698  # GeV*fm

z_Planck = 4 * phi**2  # m_e/m_P = alpha^z_Planck
me_mp_exp = m_e_exp / m_P
me_mp_fw = alpha**z_Planck

print(f"\nFramework: a_1={a1}, phi={phi:.10f}, b_1={b1}, N={N}, degree={degree}")
print(f"alpha = 1/{alpha_inv:.6f}, alpha_s(fw) = {alpha_s_fw:.8f}")
print(f"m_e/m_P: experimental = {me_mp_exp:.6e}, framework = {me_mp_fw:.6e}")
print(f"z_Planck = 4*phi^2 = {z_Planck:.6f}")

# ============================================================
# PART 1: Build 600-cell vertices
# ============================================================
print("\n" + "=" * 72)
print("PART 1: 600-CELL VERTICES")
print("=" * 72)

vertices = []

# 8 vertices: permutations of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s
        vertices.append(tuple(v))

# 16 vertices: (+-1/2, +-1/2, +-1/2, +-1/2)
for signs in iprod([0.5, -0.5], repeat=4):
    vertices.append(signs)

# 96 golden ratio vertices: even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
vals = [phi / 2.0, 0.5, 1.0 / (2.0 * phi), 0.0]

even_perms = [
    (0,1,2,3), (0,2,3,1), (0,3,1,2),
    (1,0,3,2), (1,2,0,3), (1,3,2,0),
    (2,0,1,3), (2,1,3,0), (2,3,0,1),
    (3,0,2,1), (3,1,0,2), (3,2,1,0)
]

for perm in even_perms:
    base = [vals[perm[0]], vals[perm[1]], vals[perm[2]], vals[perm[3]]]
    nonzero_idx = [i for i in range(4) if abs(base[i]) > 1e-12]
    for signs in iprod([1, -1], repeat=len(nonzero_idx)):
        v = list(base)
        for k, idx in enumerate(nonzero_idx):
            v[idx] *= signs[k]
        vertices.append(tuple(v))

# Remove duplicates
unique = []
for v in vertices:
    found = False
    for u in unique:
        if sum((a - b)**2 for a, b in zip(v, u)) < 1e-10:
            found = True
            break
    if not found:
        unique.append(v)

vertices = np.array(unique)
print(f"Vertices: {len(vertices)} (expected {N})")
assert len(vertices) == N, f"Got {len(vertices)}, expected {N}"

norms = np.sqrt(np.sum(vertices**2, axis=1))
assert np.allclose(norms, 1.0), "Not all on unit S^3!"
print("All on unit S^3: VERIFIED")

# Verify degree = 12
inner_prods = vertices @ vertices.T
# Count neighbors within chord distance < 2(1-phi/2) + epsilon
d2_nearest = 2 * (1 - phi / 2)  # = 2 - phi = 1/phi^2
for i in range(3):
    d2_all = 2 * (1 - inner_prods[i, :])
    nn = np.sum((d2_all > 1e-8) & (d2_all < d2_nearest + 0.01))
    assert nn == degree, f"Vertex {i}: degree={nn}, expected {degree}"
print(f"Degree = {degree}: VERIFIED")

# ============================================================
# PART 2: Distance spectrum and Riesz energies
# ============================================================
print("\n" + "=" * 72)
print("PART 2: DISTANCE SPECTRUM AND RIESZ ENERGIES")
print("=" * 72)

# Inner product shells from any vertex
# cos(theta) values and their multiplicities
ip0 = inner_prods[0, 1:]
ip_rounded = np.round(ip0, 8)
unique_ip = np.sort(np.unique(ip_rounded))[::-1]

# Analytical shells: (cos_theta, multiplicity, d^2, 1/d^2, label)
shells = [
    (phi/2,      12, 2 - phi,  phi**2,               "2-phi=1/phi^2"),
    (0.5,        20, 1.0,      1.0,                   "1"),
    (1/(2*phi),  12, 3 - phi,  (5 + np.sqrt(5))/10,  "3-phi"),
    (0.0,        30, 2.0,      0.5,                   "2"),
    (-1/(2*phi), 12, phi_sq,   1/phi**2,              "phi^2"),
    (-0.5,       20, 3.0,      1.0/3,                 "3"),
    (-phi/2,     12, 2 + phi,  (5 - np.sqrt(5))/10,  "2+phi"),
    (-1,          1, 4.0,      0.25,                  "4"),
]

print("\nDistance shells from any vertex:")
print(f"{'d^2':>16} {'mult':>5} {'1/d^2':>14} {'m/d^2':>14}")
print("-" * 55)
total_mult = 0
for _, m, d2, inv_d2, lab in shells:
    total_mult += m
    print(f"{lab:>16} {m:5d} {inv_d2:14.8f} {m*inv_d2:14.8f}")
print(f"{'TOTAL':>16} {total_mult:5d}")
assert total_mult == N - 1

# Compute S_p for p = 1, 2, 3, 4
print("\n--- Riesz energies S_p = sum_{i<j} |x_i - x_j|^{-p} ---")

for p in [1, 2, 3, 4]:
    sp_pv = 0  # per vertex
    for _, m, d2, _, _ in shells:
        sp_pv += m * d2**(-p/2.0)
    S_p = sp_pv * N / 2.0
    print(f"\nS_{p} per vertex = {sp_pv:.10f}")
    print(f"S_{p} total     = {S_p:.6f}")

    if p == 2:
        S2_pv = sp_pv
        S2 = S_p

# ============================================================
# PART 3: Analytical S_2
# ============================================================
print("\n" + "=" * 72)
print("PART 3: ANALYTICAL S_2")
print("=" * 72)

# S_2_per_vertex decomposition:
# Galois pairs (phi + phi' cancel):
#   shells 1+5: 12*phi^2 + 12/phi^2 = 12*(phi^2 + 2-phi) = 12*3 = 36
#   shells 3+7: 12*(5+sqrt5)/10 + 12*(5-sqrt5)/10 = 12
# Rational:
#   shells 2,4,6,8: 20 + 15 + 20/3 + 1/4

pair1 = 12 * (phi**2 + (2 - phi))
pair2 = 12 * ((5 + np.sqrt(5))/10 + (5 - np.sqrt(5))/10)
rational = 20 + 15 + 20.0/3 + 0.25

print(f"Galois pair 1 (nearest + 5th shell): 12*(phi^2 + 1/phi^2)")
print(f"  phi^2 + 1/phi^2 = phi^2 + (2-phi) = (phi+1) + (2-phi) = 3 EXACT")
print(f"  Contribution: 12 * 3 = 36")
print(f"  Numerical: {pair1:.10f}")

print(f"\nGalois pair 2 (3rd + 7th shell):")
print(f"  (5+sqrt5)/10 + (5-sqrt5)/10 = 10/10 = 1")
print(f"  Contribution: 12 * 1 = 12")
print(f"  Numerical: {pair2:.10f}")

print(f"\nRational terms: 20 + 15 + 20/3 + 1/4")
rat_exact = Fraction(20) + Fraction(15) + Fraction(20, 3) + Fraction(1, 4)
print(f"  = {rat_exact} = {float(rat_exact):.10f}")
print(f"  Numerical: {rational:.10f}")

S2_pv_exact = Fraction(36) + Fraction(12) + rat_exact
S2_exact_int = S2_pv_exact * N // 2  # integer since N/2 = 60

print(f"\nS_2 per vertex = 36 + 12 + {rat_exact} = {S2_pv_exact}")
print(f"  = {S2_pv_exact} = {float(S2_pv_exact):.10f}")
print(f"  Numerical: {S2_pv:.10f}")
print(f"  Match: {abs(float(S2_pv_exact) - S2_pv) < 1e-8}")

S2_total_exact = S2_pv_exact * Fraction(N, 2)
print(f"\nS_2 total = {S2_pv_exact} * {N}/2 = {S2_total_exact}")
print(f"  Numerical: {S2:.6f}")
print(f"  Match: {abs(float(S2_total_exact) - S2) < 0.01}")

# Key result
print(f"\n** S_2 = {S2_total_exact} (RATIONAL -- all phi dependence cancels!) **")
print(f"** S_2 is Galois-invariant: sigma(S_2) = S_2 **")

# ============================================================
# PART 4: Number analysis
# ============================================================
print("\n" + "=" * 72)
print("PART 4: FACTORIZATION AND FRAMEWORK ANALYSIS")
print("=" * 72)

S2_int = int(S2_total_exact)
print(f"\nS_2 = {S2_int}")

# Factor
def factorize(n):
    factors = {}
    d = 2
    x = abs(n)
    while d * d <= x:
        while x % d == 0:
            factors[d] = factors.get(d, 0) + 1
            x //= d
        d += 1
    if x > 1:
        factors[x] = factors.get(x, 0) + 1
    return factors

f = factorize(S2_int)
print(f"  = {' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f.items()))}")

# S_2_per_vertex = 1079/12
pv_num = S2_pv_exact.numerator
pv_den = S2_pv_exact.denominator
print(f"\nS_2 per vertex = {pv_num}/{pv_den}")
f_num = factorize(pv_num)
f_den = factorize(pv_den)
print(f"  numerator {pv_num} = {' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f_num.items()))}")
print(f"  denominator {pv_den} = {' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f_den.items()))}")

# Framework number checks
print(f"\nFramework number analysis:")
print(f"  {S2_int} = a_1 * (degree+1) * 83 = {a1} * {degree+1} * 83")
print(f"  degree + 1 = {degree+1} = a_1 + rank(E8) = 5 + 8 = 13")
print(f"  83 is PRIME, 23rd prime")

# Check if 83 appears in any framework context
print(f"\n  Is 83 a framework number?")
print(f"    N - 37 = {N - 37} = 83? 37 = prime, not obvious")
print(f"    Fibonacci: ...55, 89... 83 NOT Fibonacci")
print(f"    Lucas: ...76, 123... 83 NOT Lucas")
print(f"    a1^2 + b1^2 + a1*b1 = {a1**2 + b1**2 + a1*b1} (not 83)")
print(f"    N - a1*b1 - 7 = {N - a1*b1 - 7} (=83!)")
print(f"    Actually: 83 = N - a1*b1 - 7 = 120 - 30 - 7")

# But 7 is not a clean framework number in this context
# Let's try other decompositions of S_2
print(f"\n  Alternative decompositions of {S2_int}:")
for a in range(1, 100):
    if S2_int % a == 0:
        b = S2_int // a
        if b <= 200:
            print(f"    {a} * {b}", end="")
            # Check if both are framework-ish
            notes = []
            for x in [a, b]:
                if x == a1: notes.append(f"{x}=a1")
                elif x == b1: notes.append(f"{x}=b1")
                elif x == N: notes.append(f"{x}=N")
                elif x == degree: notes.append(f"{x}=deg")
                elif x == degree + 1: notes.append(f"{x}=deg+1")
                elif x == 30: notes.append(f"{x}=h(E8)")
            if notes:
                print(f"  [{', '.join(notes)}]", end="")
            print()

# Spectral coefficient ratios
c0 = 2640
c1 = 14880
c2 = 55920
print(f"\nRatios with spectral coefficients:")
print(f"  S_2/c_0 = {Fraction(S2_int, c0)} = {S2_int/c0:.6f}")
print(f"  S_2/c_1 = {Fraction(S2_int, c1)} = {S2_int/c1:.8f}")
print(f"  c_1/S_2 = {Fraction(c1, S2_int)} = {c1/S2_int:.6f}")
print(f"  c_2/S_2 = {Fraction(c2, S2_int)} = {c2/S2_int:.6f}")
print(f"  c_0*c_2/c_1^2 = {Fraction(c0*c2, c1**2)} (from exp328)")

# ============================================================
# PART 5: Self-consistency equations for R
# ============================================================
print("\n" + "=" * 72)
print("PART 5: SELF-CONSISTENCY EQUATIONS")
print("=" * 72)

# Laplacian eigenvalues of 600-cell
eig_data = [
    (0,           1,  "0"),
    (12 - 6*phi,  4,  "12-6phi"),
    (12 - 4*phi,  9,  "12-4phi"),
    (9,          16,  "9"),
    (12,         25,  "12"),
    (14,         36,  "14"),
    (8 + 4*phi,   9,  "8+4phi"),
    (15,         16,  "15"),
    (6 + 6*phi,   4,  "6+6phi"),
]

tr_Delta = sum(m * l for l, m, _ in eig_data)
sqrt_sum = sum(m * np.sqrt(l) for l, m, _ in eig_data if l > 0)
log_sum = sum(m * np.log(l) for l, m, _ in eig_data if l > 0)

print(f"\nSpectral data of Delta_0:")
print(f"  Tr(Delta_0) = {tr_Delta:.4f} = c_1/degree = {c1/degree:.4f}")
print(f"  sum m_k sqrt(lam_k) = {sqrt_sum:.6f}")
print(f"  sum m_k ln(lam_k) = {log_sum:.6f} (= -zeta'(0))")

# --- Condition A: E_grav/vertex = E_Casimir/vertex ---
# E_grav per vertex ~ S_2_pv / R^2 (gravitational, 4D: V~1/r^2)
# E_Casimir per vertex ~ (1/R) * Tr(Delta_0)/N (zero-point energy)
# Balance: S_2_pv/R^2 = Tr(Delta_0)/(N*R)
# => R = N * S_2_pv / Tr(Delta_0)
# Hmm, that's one version. The text says S_2/(N*R^2) = c_1/(2R).

print(f"\n--- Condition A (from text): S_2/(N*R^2) = c_1/(2*R) ---")
R_A = 2.0 * S2_int / (N * c1)
R_A_exact = Fraction(2 * S2_int, N * c1)
print(f"  R = 2*S_2/(N*c_1) = {R_A_exact} = {R_A:.10f}")
# Simplify
g_A = gcd(R_A_exact.numerator, R_A_exact.denominator)
print(f"  = {R_A_exact.numerator//g_A}/{R_A_exact.denominator//g_A}")

print(f"\n--- Condition B (minimize E = -S_2/R^2 + c_1/(2R)) ---")
print(f"  dE/dR = 2*S_2/R^3 - c_1/(2*R^2) = 0")
R_B = 4.0 * S2_int / c1
R_B_exact = Fraction(4 * S2_int, c1)
print(f"  R = 4*S_2/c_1 = {R_B_exact} = {R_B:.10f}")

print(f"\n--- Condition C (use sqrt spectrum): S_2/(N*R^2) = sqrt_sum/(2*R) ---")
R_C = 2.0 * S2_int / (N * sqrt_sum)
print(f"  R = 2*S_2/(N*sqrt_sum) = {R_C:.10f}")

print(f"\n--- Condition D (3D gravity V~1/r): S_1/(N*R) = c_1/(2*R) ---")
# For V(r)=1/r: E_grav ~ S_1/R, E_Casimir ~ 1/R
# Balance: S_1/(N*R) = c_1/(2R) => S_1/N = c_1/2 (R cancels!)
S1_pv = sum(m * d2**(-0.5) for _, m, d2, _, _ in shells)
S1 = S1_pv * N / 2.0
print(f"  S_1/N = {S1/N:.6f}")
print(f"  c_1/2 = {c1/2}")
print(f"  Ratio: {2*S1/(N*c1):.6f}")
print(f"  ** R cancels: no constraint on R (just number check) **")

# ============================================================
# PART 6: Physical scale from R
# ============================================================
print("\n" + "=" * 72)
print("PART 6: PHYSICAL SCALE ANALYSIS")
print("=" * 72)

print(f"\nm_e/m_P needed: {me_mp_exp:.6e}")
print(f"= alpha^(4*phi^2) = alpha^{z_Planck:.4f} = {me_mp_fw:.6e}")
print(f"=> need R ~ m_P/m_e ~ {m_P/m_e_exp:.2e} (in Planck units)")

for label, R_val in [("A", R_A), ("B", R_B), ("C", R_C)]:
    # Interpret: m_e/m_P = 1/R if R is in Planck units
    if R_val > 0:
        z = np.log(1.0/R_val) / np.log(alpha)
        print(f"\nCondition {label}: R = {R_val:.6e}")
        print(f"  m_e/m_P = 1/R = {1.0/R_val:.4f}")
        print(f"  log_alpha(1/R) = {z:.4f} (need {z_Planck:.4f})")
        print(f"  ** R ~ O(1) => m_e/m_P ~ O(1): FAILS **")

print(f"\n** STRUCTURAL NO-GO: **")
print(f"  Balancing power-law energies E ~ 1/R^p gives R = (pure number).")
print(f"  All pure numbers from 600-cell spectral data are O(1)-O(10).")
print(f"  Need R ~ 10^22 for m_e. Requires EXPONENTIAL mechanism.")
print(f"  The hierarchy alpha^(4*phi^2) ~ 10^-23 is exponential,")
print(f"  not polynomial. Energy balance CANNOT generate it.")

# ============================================================
# PART 7: Direction D revisited - spectral action
# ============================================================
print("\n" + "=" * 72)
print("PART 7: DIRECTION D - SPECTRAL ACTION SELF-CONSISTENCY")
print("=" * 72)

# S = c_0 f_4 Lambda^4 - c_1 f_2 Lambda^2 + c_2 f_0/2 + ...
# Critical point: Lambda^2 = c_1 f_2 / (2 c_0 f_4)

lambda_max = 6 + 6 * phi

for name, f0, f2, f4 in [("Heat kernel (f_k=1)", 1, 1, 1),
                          ("Sharp cutoff (f_4=2)", 1, 1, 2),
                          ("Chamseddine (f_4=1/2)", 1, 1, 0.5)]:
    L2 = c1 * f2 / (2.0 * c0 * f4)
    L2_exact = Fraction(c1 * f2, int(2 * c0 * f4)) if f4 == int(f4) else None

    print(f"\n--- {name} ---")
    if L2_exact:
        print(f"  Lambda^2 = c_1*f_2/(2*c_0*f_4) = {L2_exact} = {float(L2_exact):.8f}")
    else:
        print(f"  Lambda^2 = c_1*f_2/(2*c_0*f_4) = {L2:.8f}")
    print(f"  Lambda = {np.sqrt(L2):.8f}")
    print(f"  Lambda^2/lambda_max = {L2/lambda_max:.6f}")

    # f_2 * Lambda^2 (the Higgs-related parameter)
    f2L2 = f2 * L2
    print(f"  f_2*Lambda^2 = {f2L2:.6f}")

    # Check if it's a framework number
    for n in range(-5, 15):
        if abs(L2 - phi**n) / max(L2, abs(phi**n)) < 0.05:
            print(f"  ** Lambda^2 ~ phi^{n} = {phi**n:.6f} ({(L2/phi**n-1)*100:.2f}%) **")

# The key ratio c_1/(2*c_0)
ratio_key = Fraction(c1, 2 * c0)
print(f"\nc_1/(2*c_0) = {ratio_key} = {float(ratio_key):.8f}")
f_rk = factorize(ratio_key.numerator)
print(f"  {ratio_key.numerator} = {' * '.join(str(p) for p in sorted(f_rk.keys()))}")
f_rd = factorize(ratio_key.denominator)
print(f"  {ratio_key.denominator} = {' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f_rd.items()))}")

# Does Lambda^2 = 31/11 have framework meaning?
# 31 = a1^2 + b1 = 25 + 6
# 11 = a1 + b1 = L(5)
print(f"\n  31 = a_1^2 + b_1 = {a1**2} + {b1} = {a1**2 + b1}")
print(f"  11 = a_1 + b_1 = L(5) (5th Lucas number)")
print(f"  So Lambda^2 = (a_1^2 + b_1) / (a_1 + b_1) = (a_1^2 + b_1) / L(a_1)")
print(f"  ** THIS IS A CLEAN FRAMEWORK EXPRESSION! **")

# Check: Lambda^2 = 31/11 in terms of phi
# (a1^2 + b1)/(a1 + b1) = (25 + 6)/(5 + 6) = 31/11
# Not immediately a phi power but a clean RATIONAL framework number

# ============================================================
# PART 8: Direction E - Lattice norms
# ============================================================
print("\n" + "=" * 72)
print("PART 8: DIRECTION E - Z[phi] LATTICE ANALYSIS")
print("=" * 72)

# Fermion quantum numbers and derived quantities
fermions = [
    ("e",   0, 0),
    ("mu",  1, 1),
    ("tau", 1, 2),
    ("u",   3,-2),
    ("c",   2, 1),
    ("t",   4, 1),
    ("d",   1, 0),
    ("s",   1, 1),
    ("b",  -1, 4),
]

print(f"\nFermion lattice data: n_f = {a1}a + {b1}b, N(a+b*phi) = a^2+ab-b^2")
print(f"{'name':>5} {'(a,b)':>8} {'n_f':>5} {'N(a+bphi)':>10} {'|N|':>5} {'phi^n_f':>12}")
print("-" * 50)

norms = []
for name, a, b in fermions:
    nf = a1 * a + b1 * b
    norm = a**2 + a * b - b**2  # N(a + b*phi) = (a+b*phi)(a+b*phi')
    norms.append((name, a, b, nf, norm))
    print(f"{name:>5} ({a:2d},{b:2d}) {nf:5d} {norm:10d} {abs(norm):5d} {phi**nf:12.2f}")

# Check: is norm monotonically related to mass?
nf_list = [x[3] for x in norms]
norm_list = [x[4] for x in norms]
abs_norm_list = [abs(x[4]) for x in norms]

# Spearman correlation between |N| and n_f
from scipy.stats import spearmanr
rho_nf_norm, p_val = spearmanr(nf_list, abs_norm_list)
print(f"\nSpearman correlation |N| vs n_f: rho = {rho_nf_norm:.4f}, p = {p_val:.4f}")

# Key observation
print(f"\nKey observations:")
print(f"  N(c) = a_1 = 5 (charm has norm a_1!)")
print(f"  |N(t)| = |N(b)| = 19 (top and bottom have equal |norm|)")
print(f"  N(mu) = N(s) = N(d) = 1 (three fermions with unit norm)")
print(f"  N(tau) = N(u) = -1 (two with norm -1)")
print(f"  N(e) = 0 (electron has norm 0 -- the ring ZERO element)")

# The discriminant of Q(a,b) = a^2 + ab - b^2 is Delta = 1 + 4 = 5 = a_1
print(f"\nQuadratic form Q(a,b) = a^2 + ab - b^2:")
print(f"  Discriminant = 1 + 4*1 = 5 = a_1")
print(f"  Class number h(5) = 1 (unique class)")

# Theta function of this form
# Theta_Q(tau) = sum_{a,b in Z} q^{a^2+ab-b^2} where q = exp(2*pi*i*tau)
# This is INDEFINITE (signature (1,1)), not positive definite!
# Indefinite forms don't give convergent theta series in the usual sense.
print(f"\n** PROBLEM: Q has signature (1,1) -- INDEFINITE **")
print(f"  Indefinite theta series diverge (no absolute convergence)")
print(f"  Need regularization (Hecke, Zwegers, mock modular forms)")
print(f"  This is a NON-TRIVIAL obstruction to Direction E")

# Alternative: use the linear form L(a,b) = 5a + 6b
# Count lattice points with given exponent
print(f"\n--- Counting lattice points by exponent ---")
exponent_counts = {}
for a_val in range(-10, 11):
    for b_val in range(-10, 11):
        nf = a1 * a_val + b1 * b_val
        nf_abs = abs(nf)
        exponent_counts[nf_abs] = exponent_counts.get(nf_abs, 0) + 1

print(f"{'n':>4} {'count':>6} {'fermion?':>10}")
for n in sorted(exponent_counts.keys())[:30]:
    ferm = ""
    for name, a, b, nf, norm in norms:
        if nf == n:
            ferm += name + " "
    print(f"{n:4d} {exponent_counts[n]:6d} {ferm:>10}")

# The generating function sum q^|5a+6b| is NOT modular
# because 5a+6b is LINEAR, not quadratic
print(f"\n** Linear form 5a+6b: generating function NOT modular **")
print(f"   (Modularity requires quadratic forms)")

# ============================================================
# PART 9: Are there ANY framework numbers in S_2?
# ============================================================
print("\n" + "=" * 72)
print("PART 9: DEEPER S_2 STRUCTURE")
print("=" * 72)

# S_2 = 5395. Let me check ALL possible framework expressions
# that evaluate to 5395

# S_2_per_vertex = 1079/12
# Per-vertex, per-neighbor: 1079/(12*119) = 1079/1428

# Connection to c_k:
# S_2 * degree = 5395 * 12 = 64740 = c_1 * 31/69 ... nah
# S_2 * 12 = 64740. 64740/c_1 = 64740/14880 = 4.3508... not clean

# But: c_1 = 14880, S_2 = 5395
# c_1 - 2*S_2 = 14880 - 10790 = 4090 = 2 * 5 * 409
# c_1 + S_2 = 20275 = 5^2 * 811
# 3*S_2 = 16185, 3*S_2 - c_1 = 1305 = 5*261 = 5*9*29

# Perhaps S_2 relates to S_2 on OTHER polytopes?
# For comparison, compute S_2 for simpler polytopes

# 5-cell (simplex in R^4): 5 vertices, all equidistant
# On unit S^3: inner product between vertices = -1/4
# d^2 = 2(1+1/4) = 5/2
# S_2 = C(5,2) * 2/5 = 10 * 2/5 = 4

# 8-cell (hypercube): 16 vertices
# 24-cell: 24 vertices

# 600-cell: S_2 = 5395
# Normalize per pair: S_2/C(120,2) = 5395/7140
s2_per_pair = Fraction(S2_int, N * (N-1) // 2)
print(f"S_2 per pair = {s2_per_pair} = {float(s2_per_pair):.8f}")
f_pp = factorize(s2_per_pair.numerator)
print(f"  numerator {s2_per_pair.numerator} = {' * '.join(str(p) for p in sorted(f_pp.keys()))}")
f_pd = factorize(s2_per_pair.denominator)
print(f"  denominator {s2_per_pair.denominator} = {' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f_pd.items()))}")

# S_2 * 12 / N = 5395 * 12 / 120 = 539.5 -- not integer even
print(f"\nS_2 * degree / N = {S2_int * degree / N}")

# Check: is S_2 related to Tr(Delta_0^{-1})?
# Tr(Delta_0^{-1}) = sum m_k / lambda_k (over nonzero eigenvalues)
tr_inv = sum(m / l for l, m, _ in eig_data if l > 0)
print(f"\nTr(Delta_0^{{-1}}) = {tr_inv:.10f}")
print(f"S_2 / Tr(Delta_0^{{-1}}) = {S2_int / tr_inv:.6f}")
print(f"N * Tr(Delta_0^{{-1}}) = {N * tr_inv:.6f}")

# Interesting: Green's function at coincident point
# G(v,v) = (1/N) * sum_{k>0} m_k/lambda_k (averaged)
G_vv = tr_inv / N
print(f"Green's function G(v,v) = Tr(Delta^{{-1}})/N = {G_vv:.10f}")
print(f"S_2 * G(v,v) = {S2_int * G_vv:.6f}")

# ============================================================
# PART 10: Riesz energy for OTHER regular 4-polytopes
# ============================================================
print("\n" + "=" * 72)
print("PART 10: COMPARISON WITH OTHER POLYTOPES")
print("=" * 72)

# 24-cell: 24 vertices on unit S^3
# Inner products: 0 (8 vertices at pi/2), +-1/2 (each 6), +-1 (each 1)
# Actually: from any vertex, the 23 others are:
# 8 at cos=1/2 (d^2=1), 6 at cos=0 (d^2=2), 8 at cos=-1/2 (d^2=3), 1 at cos=-1 (d^2=4)
# S_2_24_pv = 8/1 + 6/2 + 8/3 + 1/4 = 8 + 3 + 8/3 + 1/4 = 11 + 8/3 + 1/4 = (132+32+3)/12 = 167/12

s2_24_pv = 8.0/1 + 6.0/2 + 8.0/3 + 1.0/4
S2_24 = s2_24_pv * 24 / 2
s2_24_exact = Fraction(167, 12)
S2_24_exact = s2_24_exact * Fraction(24, 2)

print(f"24-cell: N=24, S_2_pv = {s2_24_exact} = {float(s2_24_exact):.6f}")
print(f"  S_2 = {S2_24_exact} = {float(S2_24_exact):.6f}")

# Normalized: S_2 / C(N,2)
s2_24_norm = S2_24_exact / Fraction(24 * 23, 2)
s2_120_norm = Fraction(S2_int, N * (N - 1) // 2)
print(f"  S_2/C(N,2): 24-cell = {float(s2_24_norm):.6f}, 600-cell = {float(s2_120_norm):.6f}")
print(f"  Ratio: {float(s2_120_norm/s2_24_norm):.6f}")

# 120-cell: 600 vertices -- too complex to build here, just note
print(f"\n120-cell: N=600 (dual of 600-cell) -- not computed")

# ============================================================
# PART 11: Lambda^2 = 31/11 deeper analysis
# ============================================================
print("\n" + "=" * 72)
print("PART 11: Lambda^2 = (a_1^2 + b_1)/(a_1 + b_1) ANALYSIS")
print("=" * 72)

L2 = Fraction(31, 11)
print(f"\nLambda^2 = {L2} = {float(L2):.10f}")
print(f"  = (a_1^2 + b_1)/(a_1 + b_1)")
print(f"  = ({a1}^2 + {b1})/({a1} + {b1})")
print(f"  = (25 + 6)/(5 + 6)")

# Lambda
L_val = np.sqrt(float(L2))
print(f"\nLambda = sqrt(31/11) = {L_val:.10f}")

# Check against phi powers
for n_num in range(-20, 21):
    for n_den in range(1, 13):
        trial = phi**(n_num / n_den)
        if abs(L_val - trial) / L_val < 0.01:
            print(f"  Lambda ~ phi^({n_num}/{n_den}) = {trial:.6f} ({(L_val/trial-1)*100:.3f}%)")

# f_2 * Lambda^2 in the spectral action determines the Higgs VEV
# In Connes' NCG: m_H^2 = 2*lambda * f_2*Lambda^2 / (pi^2 * a)
# where lambda is the quartic coupling and a is a normalization

# The key insight: Lambda^2 = 31/11 is a RATIONAL framework number
# Compare with lambda_max = 6+6*phi = 15.708...
print(f"\nLambda^2 / lambda_max = {float(L2)/(6+6*phi):.6f}")
print(f"  = (31/11)/(6+6*phi) = 31/(11*6*(1+phi)) = 31/(66*phi^2)")
print(f"  = {31/(66*phi**2):.6f}")

# How many eigenvalues below Lambda^2?
count_below = sum(m for l, m, _ in eig_data if l <= float(L2))
count_above = sum(m for l, m, _ in eig_data if l > float(L2))
print(f"\nEigenvalues of Delta_0 below Lambda^2={float(L2):.4f}:")
for l, m, lab in eig_data:
    if l <= float(L2) + 0.01:
        marker = " <-- cutoff" if abs(l - float(L2)) < 0.5 else ""
        print(f"  {lab:>12}: lambda={l:.4f}, mult={m}{marker}")
print(f"  Total below: {count_below}, above: {count_above}")
print(f"  Fraction below: {count_below}/{N} = {count_below/N:.4f}")

# ============================================================
# PART 12: Honest assessment
# ============================================================
print("\n" + "=" * 72)
print("PART 12: HONEST ASSESSMENT")
print("=" * 72)

print("""
DIRECTION F (Cohn-Kumar energy -> m_e): NEGATIVE (structural no-go)

  1. S_2 = 5395 is EXACTLY RATIONAL (all phi-dependence cancels).
     Galois-invariant: sigma(S_2) = S_2. This is topological.

  2. S_2 = 5 * 13 * 83. Factor 83 is NOT a framework number.
     5 = a_1, 13 = degree + 1 = a_1 + rank(E8).

  3. ** STRUCTURAL NO-GO **: Energy balance between power-law terms
     (E_grav ~ 1/R^2, E_Casimir ~ 1/R) gives R = pure number ~ O(1).
     The mass hierarchy m_e/m_P ~ 10^{-23} requires EXPONENTIAL
     suppression (alpha^{4*phi^2}), not polynomial ratios.
     No energy balance on finite polytopes can produce this.

  4. The no-go is GENERAL: applies to ANY Riesz energy S_p on
     ANY regular polytope balanced against ANY spectral quantity.

DIRECTION E (theta function on Z[phi]): NEGATIVE (two obstructions)

  1. The mass exponent 5a+6b is LINEAR, not quadratic.
     Linear forms don't have modular properties.

  2. The actual norm Q(a,b) = a^2+ab-b^2 is INDEFINITE (signature 1,1).
     Indefinite theta series diverge without regularization.
     Even regularized (mock modular), the connection to masses is absent.

  3. Norm values {0,1,-1,-1,5,19,1,1,-19} show NO correlation with
     mass hierarchy. |N(t)| = |N(b)| = 19 but masses differ by phi^7.

DIRECTION D (spectral action self-consistency): PARTIALLY POSITIVE

  1. Lambda^2 = c_1/(2*c_0) = 31/11 (for f_4=1 heat kernel)
     = (a_1^2 + b_1)/(a_1 + b_1) -- CLEAN FRAMEWORK EXPRESSION!
     This determines the spectral action cutoff in natural units.

  2. But this STILL doesn't fix m_e. The dimensionful cutoff
     Lambda_phys = Lambda * (energy scale) requires knowing m_e first.
     The self-referential bootstrap FAILS (same as exp327).

POSITIVE DELIVERABLES:
  - S_2 = 5395 (new topological invariant, Galois-invariant)
  - Lambda^2 = (a_1^2+b_1)/(a_1+b_1) = 31/11 (clean expression)
  - Structural no-go for polynomial m_e derivation
  - Z[phi] norm analysis for all fermions
""")

elapsed = time.time() - t0
print(f"Completed in {elapsed:.1f}s")

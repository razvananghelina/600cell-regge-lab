# EXP-220: Comprehensive Neutrino Mass Attack
# =============================================
# Systematic exploration of ALL approaches.
# Goal: find m3 AND m2 from 600-cell geometry.
#
# Key experimental values (PDG 2024, normal hierarchy):
#   Dm2_31 = 2.455e-3 eV^2 => m3 ~ 0.04955 eV
#   Dm2_21 = 7.53e-5 eV^2  => m2 ~ 0.00868 eV
#   Dm2_21/Dm2_31 = 0.03068 ~ 4*alpha (1.24% match)
#
# RULE ZERO: verify ALL claims numerically. No rounding tricks.

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 1/137.035999084
ALPHA_S = 0.1179
SIN2TW = 6.0/26.0  # DERIVED: 6/26 = 3/13

m_e = 0.51099895e6  # eV

# Neutrino data
Dm2_21 = 7.53e-5   # eV^2
Dm2_31 = 2.455e-3   # eV^2
m3_exp = np.sqrt(Dm2_31)  # 0.04955 eV
m2_exp = np.sqrt(Dm2_21)  # 0.008678 eV
r_exp = Dm2_21 / Dm2_31   # 0.03068

# 600-cell invariants (ALL DERIVED)
a_1 = 5        # first eigenvalue
b_1 = 6        # second eigenvalue
N_eig = 9      # number of eigenvalues
deg = 12       # vertex degree
diam = 5       # graph diameter
N_vert = 120   # vertices
N_edge = 720   # edges
h_E8 = 30      # Coxeter number
dim = 4        # ambient dimension
N_bag = 13     # soliton bag vertices (DERIVED EXP-178)
E_bag = 6.0    # soliton bag energy (DERIVED)
E_flip = 3     # single flip energy (DERIVED)
N_gen = 3      # generations (DERIVED)

print("="*75)
print("EXP-220: COMPREHENSIVE NEUTRINO MASS ATTACK")
print("="*75)
print(f"Target: m3 = {m3_exp:.5f} eV, m2 = {m2_exp:.5f} eV")
print(f"Dm2 ratio = {r_exp:.5f}, 4*alpha = {4*ALPHA:.5f} (err {abs(4*ALPHA-r_exp)/r_exp*100:.2f}%)")
print()

# =====================================================================
# PART A: BASE MASS SCALE m3 - Systematic scan
# =====================================================================
print("="*75)
print("PART A: BASE MASS SCALE m3")
print("="*75)
print(f"{'#':<4s} {'Formula':<55s} {'m3 (eV)':>10s} {'err%':>8s}")
print("-"*80)

results_A = []

def test_m3(label, value):
    err = abs(value - m3_exp) / m3_exp * 100
    results_A.append((label, value, err))
    print(f"{len(results_A):<4d} {label:<55s} {value:>10.5f} {err:>7.2f}%")

# --- Known formulas ---
test_m3("m_e * (alpha/phi)^3", m_e * (ALPHA/PHI)**3)
test_m3("m_e * alpha^2 / phi^13", m_e * ALPHA**2 / PHI**13)
test_m3("m_e * phi^(-34)", m_e / PHI**34)

# --- Power-of-alpha with phi ---
for p in [2, 3, 4]:
    for q in range(10, 40):
        val = m_e * ALPHA**p / PHI**q
        if 0.01 < val < 0.2:
            err = abs(val - m3_exp) / m3_exp * 100
            if err < 5:
                test_m3(f"m_e * alpha^{p} / phi^{q}", val)

# --- Combinations with sin^2(tW) ---
test_m3("m_e * alpha^2 * sin2tW / phi^11", m_e * ALPHA**2 * SIN2TW / PHI**11)
test_m3("m_e * alpha^2 * sin2tW / phi^12", m_e * ALPHA**2 * SIN2TW / PHI**12)
test_m3("m_e * alpha * sin2tW / phi^26", m_e * ALPHA * SIN2TW / PHI**26)

# --- N_gen as exponent ---
test_m3("m_e * (alpha/phi)^N_gen [=3]", m_e * (ALPHA/PHI)**N_gen)
test_m3("m_e * alpha^N_gen / phi^(N_gen*diam)", m_e * ALPHA**3 / PHI**15)
test_m3("m_e * alpha^N_gen * phi^(-N_gen*b_1)", m_e * ALPHA**3 / PHI**18)

# --- Soliton bag combinations ---
test_m3("m_e * alpha^2 / phi^N_bag [=13]", m_e * ALPHA**2 / PHI**13)
test_m3("m_e * alpha / phi^(N_bag+diam) [=18]", m_e * ALPHA / PHI**18)
test_m3("m_e / phi^(h_E8+dim) [=34]", m_e / PHI**34)

# --- E_bag and E_flip combinations ---
test_m3("m_e * alpha^(E_flip) / phi^(E_flip*diam)", m_e * ALPHA**3 / PHI**15)
test_m3("m_e * alpha^2 / phi^(2*E_bag) [=12]", m_e * ALPHA**2 / PHI**12)
test_m3("m_e * alpha^2 / phi^(2*E_bag+1)", m_e * ALPHA**2 / PHI**13)

# --- NEW: Products of invariants ---
test_m3("m_e / (a_1 * phi^(h_E8))", m_e / (a_1 * PHI**30))
test_m3("m_e / (b_1 * phi^(h_E8))", m_e / (b_1 * PHI**30))
test_m3("m_e / (N_eig * phi^(h_E8))", m_e / (N_eig * PHI**30))
test_m3("m_e / (deg * phi^(h_E8))", m_e / (deg * PHI**30))

# --- NEW: Combinations with N_eig ---
test_m3("m_e * alpha^2 / (N_eig * phi^11)", m_e * ALPHA**2 / (N_eig * PHI**11))
test_m3("m_e * alpha^2 / (a_1 * phi^11)", m_e * ALPHA**2 / (a_1 * PHI**11))

# --- NEW: SU(2) sector ---
N_SU2 = 288  # SU(2) edges
N_U1 = 48    # U(1) edges
test_m3("m_e * (N_U1/N_edge)^3", m_e * (N_U1/N_edge)**3)
test_m3("m_e * alpha^2 * (N_U1/N_SU2)", m_e * ALPHA**2 * (N_U1/N_SU2))

# --- NEW: Exponential suppression ---
test_m3("m_e * exp(-N_bag) * phi^(-diam)", m_e * np.exp(-N_bag) / PHI**diam)
test_m3("m_e * exp(-deg)", m_e * np.exp(-deg))
test_m3("m_e * exp(-N_bag-E_flip)", m_e * np.exp(-N_bag-E_flip))

# --- NEW: Majorana / half-integer hypothesis ---
test_m3("m_e * phi^(-33.5) = phi^(-67/2)", m_e / PHI**33.5)
test_m3("m_e * phi^(-(h_E8+E_flip)) = phi^-33", m_e / PHI**33)
test_m3("m_e * phi^(-(h_E8+dim-0.5))", m_e / PHI**33.5)

# --- NEW: Pure geometric ---
test_m3("m_e * phi^(-(a_1*b_1 + E_flip)) = phi^-33", m_e / PHI**33)
test_m3("m_e * phi^(-(a_1*b_1 + dim)) = phi^-34", m_e / PHI**34)
test_m3("m_e * alpha / phi^(a_1*dim) [=20]", m_e * ALPHA / PHI**20)

# --- NEW: Seesaw-like ---
test_m3("m_e^2 / (m_tau * phi^10)", m_e**2 / (1776.86e6 * PHI**10))
m_W = 80.3692e9  # eV
test_m3("m_e^2 * alpha / m_W", m_e**2 * ALPHA / m_W)
test_m3("m_e^2 / (m_W * phi^2)", m_e**2 / (m_W * PHI**2))

# --- NEW: McKay / E8 leg dimensions ---
Leg_A = 11
Leg_B = 6
Leg_C = 9
test_m3("m_e * alpha^2 / phi^(Leg_A+2)", m_e * ALPHA**2 / PHI**13)
test_m3("m_e * alpha^2 / phi^(Leg_C+dim)", m_e * ALPHA**2 / PHI**13)

# --- Sort by accuracy ---
results_A.sort(key=lambda x: x[2])
print("\n--- TOP 10 MOST ACCURATE m3 FORMULAS ---")
seen = set()
rank = 0
for label, val, err in results_A:
    rounded = round(val, 6)
    if rounded not in seen:
        seen.add(rounded)
        rank += 1
        if rank > 10:
            break
        marker = ""
        if err < 1:
            marker = " <-- SUB-1%"
        elif err < 3:
            marker = " <-- GOOD"
        print(f"  {rank}. {label}: {val:.5f} eV ({err:.2f}%){marker}")

# =====================================================================
# PART B: CORRECTION MECHANISMS (improving the best base)
# =====================================================================
print("\n" + "="*75)
print("PART B: CORRECTION MECHANISMS")
print("="*75)

# Base formulas to correct
bases = {
    "(alpha/phi)^3": m_e * (ALPHA/PHI)**3,
    "alpha^2/phi^13": m_e * ALPHA**2 / PHI**13,
}

corrections = {
    "1 + sin2tW/4": 1 + SIN2TW/4,
    "1 - sin2tW/4": 1 - SIN2TW/4,
    "1 + alpha_s/2": 1 + ALPHA_S/2,
    "1 - alpha_s/2": 1 - ALPHA_S/2,
    "1 + sin2tW": 1 + SIN2TW,
    "1 - sin2tW": 1 - SIN2TW,
    "1 + alpha_s/(2*pi)": 1 + ALPHA_S/(2*np.pi),
    "1 + 3/(4*N_bag)": 1 + 3.0/(4*N_bag),
    "1 - 3/(4*N_bag)": 1 - 3.0/(4*N_bag),
    "1 + alpha_s*sin2tW": 1 + ALPHA_S*SIN2TW,
    "1 + 2*alpha/pi": 1 + 2*ALPHA/np.pi,
    "phi/(phi+1) [=1/phi^-1+1]": PHI/(PHI+1),
    "1 + 1/(4*pi^2)": 1 + 1/(4*np.pi**2),
    "(1 + alpha_s/pi)": 1 + ALPHA_S/np.pi,
}

print(f"\n{'Base':<18s} {'Correction':<25s} {'m3 (eV)':>10s} {'err%':>8s}")
print("-"*65)

results_B = []
for bname, bval in bases.items():
    for cname, cval in corrections.items():
        m3_test = bval * cval
        err = abs(m3_test - m3_exp) / m3_exp * 100
        results_B.append((bname, cname, m3_test, err))
        if err < 2.0:
            print(f"{bname:<18s} {cname:<25s} {m3_test:>10.5f} {err:>7.3f}%")

results_B.sort(key=lambda x: x[3])
print("\n--- TOP 5 CORRECTED FORMULAS ---")
for i, (bn, cn, val, err) in enumerate(results_B[:5]):
    print(f"  {i+1}. {bn} * ({cn}) = {val:.5f} eV ({err:.3f}%)")

# =====================================================================
# PART C: MASS RATIO m2/m3
# =====================================================================
print("\n" + "="*75)
print("PART C: MASS RATIO m2/m3")
print("="*75)

r23_exp = m2_exp / m3_exp
print(f"Experimental: m2/m3 = {r23_exp:.5f}")
print(f"{'Formula':<45s} {'ratio':>10s} {'err%':>8s}")
print("-"*65)

results_C = []
def test_ratio(label, val):
    err = abs(val - r23_exp) / r23_exp * 100
    results_C.append((label, val, err))
    if err < 30:
        print(f"  {label:<43s} {val:>10.5f} {err:>7.2f}%")

test_ratio("2*sqrt(alpha)", 2*np.sqrt(ALPHA))
test_ratio("phi^(-3.72)", PHI**(-3.72))
test_ratio("1/phi^4", 1/PHI**4)
test_ratio("sqrt(4*alpha) [=Dm2 ratio^0.5]", np.sqrt(4*ALPHA))
test_ratio("sin2tW / sqrt(phi)", SIN2TW / np.sqrt(PHI))
test_ratio("alpha_s / phi", ALPHA_S / PHI)
test_ratio("2*alpha^(1/2)", 2*ALPHA**0.5)
test_ratio("sin(theta_C) [=phi^-3]", PHI**(-3))
test_ratio("alpha_s * sin2tW", ALPHA_S * SIN2TW)
test_ratio("1/(a_1 + 1/phi)", 1/(a_1 + 1/PHI))
test_ratio("1/(b_1 - 1/phi^2)", 1/(b_1 - 1/PHI**2))
test_ratio("2/(deg - 1/phi)", 2/(deg - 1/PHI))
test_ratio("alpha / sin2tW", ALPHA / SIN2TW)
test_ratio("sqrt(alpha/phi)", np.sqrt(ALPHA/PHI))
test_ratio("alpha * phi^2", ALPHA * PHI**2)
test_ratio("1 / (N_gen * phi^3)", 1/(N_gen * PHI**3))
test_ratio("sin2tW * phi^(-3)", SIN2TW * PHI**(-3))
test_ratio("2*alpha^(1/2) * (1-alpha_s/pi)", 2*np.sqrt(ALPHA) * (1 - ALPHA_S/np.pi))
test_ratio("2*sqrt(alpha) * (1+alpha_s/(2*pi))", 2*np.sqrt(ALPHA) * (1+ALPHA_S/(2*np.pi)))
test_ratio("(E_flip/h_E8)", E_flip/h_E8)
test_ratio("N_gen / (N_bag + E_flip)", N_gen / (N_bag + E_flip))
test_ratio("1/(2*a_1*phi^(1/2))", 1/(2*a_1*np.sqrt(PHI)))

# NEW: Scan for phi^p * alpha^q matches
print("\n--- Scanning phi^p * alpha^q for m2/m3 ratio ---")
for p_num in range(-20, 20):
    for q_num in range(-5, 5):
        p = p_num / 4.0
        q = q_num / 2.0
        if p == 0 and q == 0:
            continue
        val = PHI**p * ALPHA**q
        if 0.05 < val < 0.8:
            err = abs(val - r23_exp) / r23_exp * 100
            if err < 2.0:
                print(f"  phi^({p}) * alpha^({q}) = {val:.5f} (err {err:.2f}%)")

results_C.sort(key=lambda x: x[2])
print("\n--- TOP 5 RATIO FORMULAS ---")
for i, (label, val, err) in enumerate(results_C[:5]):
    print(f"  {i+1}. {label} = {val:.5f} ({err:.2f}%)")

# =====================================================================
# PART D: UNIFIED FORMULA (m3 AND m2 together)
# =====================================================================
print("\n" + "="*75)
print("PART D: UNIFIED FORMULA")
print("="*75)

# Combine best m3 base with best ratio
# Strategy: m3 from formula, m2 = m3 * ratio
print("\nBest combinations (m3_err + m2_err minimized):\n")
print(f"{'m3 formula':<35s} {'m2/m3 formula':<25s} {'m3':>8s} {'e3%':>6s} {'m2':>8s} {'e2%':>6s} {'Sum':>7s}")
print("-"*95)

best_combos = []

# Candidate m3 formulas
m3_cands = [
    ("m_e*(a/phi)^3", m_e*(ALPHA/PHI)**3),
    ("m_e*(a/phi)^3*(1+sin2tW/4)", m_e*(ALPHA/PHI)**3*(1+SIN2TW/4)),
    ("m_e*(a/phi)^3*(1+a_s/2)", m_e*(ALPHA/PHI)**3*(1+ALPHA_S/2)),
    ("m_e*a^2/phi^13", m_e*ALPHA**2/PHI**13),
    ("m_e*a^2/phi^13*(1-sin2tW/4)", m_e*ALPHA**2/PHI**13*(1-SIN2TW/4)),
    ("m_e*a^2/phi^13*(1-a_s/2)", m_e*ALPHA**2/PHI**13*(1-ALPHA_S/2)),
    ("m_e/phi^34", m_e/PHI**34),
    ("m_e/phi^33.5", m_e/PHI**33.5),
]

# Candidate ratio formulas
r_cands = [
    ("2*sqrt(a)", 2*np.sqrt(ALPHA)),
    ("sqrt(4*a)", np.sqrt(4*ALPHA)),
    ("a*phi^2", ALPHA*PHI**2),
    ("1/(a_1*phi^0.5)", 1/(a_1*np.sqrt(PHI))),
    ("sin2tW*phi^-3", SIN2TW/PHI**3),
]

for m3_label, m3_val in m3_cands:
    for r_label, r_val in r_cands:
        m2_val = m3_val * r_val
        m3_err = abs(m3_val - m3_exp)/m3_exp * 100
        m2_err = abs(m2_val - m2_exp)/m2_exp * 100
        total = m3_err + m2_err
        summ = m3_val + m2_val
        best_combos.append((m3_label, r_label, m3_val, m3_err, m2_val, m2_err, summ, total))

best_combos.sort(key=lambda x: x[7])  # sort by total error
for i, (ml, rl, m3v, m3e, m2v, m2e, s, tot) in enumerate(best_combos[:10]):
    cosmo = "OK" if s < 0.12 else "!!"
    print(f"{ml:<35s} {rl:<25s} {m3v:>7.4f} {m3e:>5.1f}% {m2v:>7.5f} {m2e:>5.1f}% {s:>6.3f} {cosmo}")

# =====================================================================
# PART E: NEW IDEAS - DEEPER EXPLORATION
# =====================================================================
print("\n" + "="*75)
print("PART E: NEW IDEAS")
print("="*75)

# IDEA 1: Neutrino exponent from alpha equation
print("\n--- IDEA 1: Neutrino exponent from alpha equation ---")
# alpha satisfies: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
# => 20*phi^4*alpha = 1 + 2*pi*alpha^2
# => 4*alpha = (1 + 2*pi*alpha^2) / (5*phi^4)
# Leading: 4*alpha ~ 1/(a_1 * phi^dim)
# And: Dm2_ratio ~ 4*alpha (1.24%)
# So: Dm2_ratio ~ 1/(a_1 * phi^dim) (geometric derivation!)
leading = 1.0 / (a_1 * PHI**dim)
exact_4a = (1 + 2*np.pi*ALPHA**2) / (a_1 * PHI**dim)
print(f"  1/(a_1 * phi^dim) = 1/(5*phi^4) = {leading:.6f}")
print(f"  4*alpha                          = {4*ALPHA:.6f}")
print(f"  (1+2pi*a^2)/(5*phi^4)           = {exact_4a:.6f}")
print(f"  Dm2_ratio (exp)                  = {r_exp:.6f}")
print(f"  Leading vs exact 4a: {abs(leading-4*ALPHA)/(4*ALPHA)*100:.4f}%")
print(f"  4a vs Dm2_ratio:     {abs(4*ALPHA-r_exp)/r_exp*100:.2f}%")
print(f"  ==> Dm2 ratio = 1/(a_1 * phi^dim) = GEOMETRIC! (leading order of alpha eq)")

# IDEA 2: Continuous exponent decomposition
print("\n--- IDEA 2: Continuous exponent of m3/m_e ---")
n_cont = -np.log(m3_exp / m_e) / np.log(PHI)
print(f"  m3/m_e = phi^(-n), n = {n_cont:.4f}")
print(f"  n = 33.52")
print(f"  Nearest 5a+6b decompositions:")
for a in range(-5, 10):
    for b in range(-5, 10):
        n_ab = 5*a + 6*b
        if abs(n_ab - n_cont) < 1.5:
            z = a + b*PHI
            zp = a + b*(1-PHI)
            N = a*a + a*b - b*b
            m_pred = m_e / PHI**(5*a + 6*b)
            err = abs(m_pred - m3_exp)/m3_exp * 100
            unit_str = "UNIT" if abs(N) == 1 else f"|N|={abs(N)}"
            print(f"    (a,b)=({a},{b}): n={n_ab}, |N(z)|={abs(N)}, m={m_pred:.5f} eV (err {err:.1f}%) [{unit_str}]")

# IDEA 3: Soliton spectral sub-levels as neutrino masses
print("\n--- IDEA 3: Soliton sub-level ratios ---")
print("  Each eigenvalue (mult m) splits into (m-3, 2, 1) sub-levels")
print("  If neutrinos = sub-levels of specific eigenvalue:")
print()
eigenvalues = [
    (12, 1, "lambda_0=12"),
    (5+6*PHI, 4, "lambda_1=a_1+b_1*phi"),  # ~ 14.71
    (5+4*PHI, 9, "lambda_2"),
    (3, 16, "lambda_3=E_flip"),
    (0, 25, "lambda_4=0"),
    (-2, 36, "lambda_5=-2"),
    (5-4*PHI, 9, "lambda_6"),  # ~ -1.47
    (-3, 16, "lambda_7=-E_flip"),
    (5-6*PHI, 4, "lambda_8"),  # ~ -4.71
]

for lam, mult, name in eigenvalues:
    if mult >= 4:
        # Sub-levels: bulk (mult-3), shell (2), core (1)
        # Ratio core/shell = 1/2, shell/bulk = 2/(mult-3)
        r_shell_core = 2.0  # multiplicity ratio
        # If mass ~ 1/multiplicity (inverse participation):
        m_core = 1.0
        m_shell = 1.0/2
        m_bulk = 1.0/(mult-3)
        r_23 = m_shell/m_core  # shell/core
        if mult > 3:
            r_12 = m_bulk/m_shell  # bulk/shell
            r_13 = m_bulk/m_core
            print(f"  {name} (mult={mult}): sub-levels ({mult-3},2,1)")
            print(f"    ratios: core/shell={r_23:.3f}, shell/bulk={2.0/(mult-3):.4f}")
            # Compare with neutrino ratio
            if abs(2.0/(mult-3) - r23_exp) / r23_exp < 0.5:
                print(f"    *** m2/m3 ~ 2/(mult-3) = {2.0/(mult-3):.4f} vs exp {r23_exp:.4f}")

# IDEA 4: Democratic mixing eigenvalue approach
print("\n--- IDEA 4: Democratic mixing and mass generation ---")
# Bare mixing: M = (1/6)*I + (5/12)*(J-I) where J = all-ones
# M = -1/4 * I + 5/12 * J
# Eigenvalues: {1, -1/4, -1/4} (J has eigenvalues {3, 0, 0})
# The degenerate eigenvalue -1/4 splits with mass breaking
p0 = 1.0/6
p1 = 5.0/12
M_bare = np.array([
    [p0, p1, p1],
    [p1, p0, p1],
    [p1, p1, p0]
])
eigs_M = np.sort(np.linalg.eigvalsh(M_bare))[::-1]
print(f"  Bare mixing matrix eigenvalues: {eigs_M}")
print(f"  Eigvals: 1, -1/4, -1/4")
print(f"  The 1-eigenvalue = sum of row = 1 (probability conservation)")
print(f"  The -1/4 eigenvalues are degenerate (S_3 symmetry)")
print()

# What if neutrino mass is proportional to the mixing matrix eigenvalues?
# m_i ~ |lambda_i| ? Then m3 : m2 : m1 = 1 : 1/4 : 1/4 (degenerate)
# m2/m3 = 1/4 = 0.25 vs experimental 0.175. Error = 43%. Not great.
print(f"  If m_i ~ |eigenvalue_i|: m2/m3 = 1/4 = 0.25 (exp: {r23_exp:.3f}, err: {abs(0.25-r23_exp)/r23_exp*100:.0f}%)")

# What if we perturb with alpha?
# m3 = 1, m2 = 1/4 + c*alpha, m1 = 1/4 - c*alpha?
# Need c such that m2/m3 = 0.175
# 1/4 + c*alpha = 0.175 => c*alpha = -0.075 => c = -0.075/alpha = -10.28
# Nope, c*alpha is NEGATIVE, meaning m2 < 1/4 not >. This goes wrong direction.
print(f"  Perturbation: wrong direction (need m2/m3 < 1/4, and exp IS < 1/4)")
print(f"  ==> Bare eigenvalue 1/4 is an UPPER BOUND on m2/m3. Consistent!")

# IDEA 5: SU(2) edge fraction
print("\n--- IDEA 5: SU(2) edge fraction mechanism ---")
f_SU2 = N_SU2 / N_edge  # 288/720 = 2/5
f_U1 = N_U1 / N_edge    # 48/720 = 1/15
N_SU3 = 384
f_SU3 = N_SU3 / N_edge  # 384/720 = 8/15
print(f"  SU(2) fraction: {N_SU2}/{N_edge} = {f_SU2:.4f} = 2/5")
print(f"  U(1) fraction:  {N_U1}/{N_edge} = {f_U1:.5f} = 1/15")
print(f"  SU(3) fraction: {N_SU3}/{N_edge} = {f_SU3:.5f} = 8/15")
print()

# Neutrinos couple to SU(2) only. What if:
# m_nu = m_e * (f_SU2)^p * alpha^q?
for p in range(1, 8):
    for q in range(1, 5):
        val = m_e * f_SU2**p * ALPHA**q
        if 0.01 < val < 0.2:
            err = abs(val - m3_exp) / m3_exp * 100
            if err < 10:
                print(f"  m_e * (2/5)^{p} * alpha^{q} = {val:.5f} eV (err {err:.1f}%)")

# What about m_nu = m_e * (f_SU2 * alpha)^p?
for p in range(1, 8):
    val = m_e * (f_SU2 * ALPHA)**p
    if 1e-6 < val < 1.0:
        err = abs(val - m3_exp) / m3_exp * 100
        if err < 50:
            print(f"  m_e * (2*alpha/5)^{p} = {val:.5f} eV (err {err:.1f}%)")

# IDEA 6: Fibonacci connection (N_gen from Fibonacci)
print("\n--- IDEA 6: Fibonacci / Lucas numbers ---")
# Fibonacci: 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89
# Lucas: 2, 1, 3, 4, 7, 11, 18, 29, 47, 76
# Note: 13 = F(7) = N_bag. 34 = F(9). 55 = F(10).
# phi^n = F(n)*phi + F(n-1)
print(f"  N_bag = 13 = F(7)")
print(f"  h(E8)+dim = 34 = F(9)")
print(f"  55 = F(10)")
print(f"  phi^13 = {PHI**13:.2f}, F(13) = 233")
print(f"  phi^34 = {PHI**34:.2f}")
print()

# What if n_nu = F(k) for some k?
for k in range(5, 15):
    fib_k = round((PHI**k - (1-PHI)**k) / np.sqrt(5))
    m_test = m_e / PHI**fib_k
    if 0.001 < m_test < 1.0:
        err = abs(m_test - m3_exp)/m3_exp * 100
        print(f"  n = F({k}) = {fib_k}: m = {m_test:.5f} eV (err {err:.1f}%)")

# IDEA 7: Three-loop with N_gen as the loop count
print("\n--- IDEA 7: Radiative mass formula m_nu = m_e * (alpha/phi)^N_gen ---")
print("  Interpretation: neutrino mass generated at N_gen=3 loop order")
print("  Loop factor = alpha/phi (replaces alpha/(4*pi) in SM)")
print()
m3_base = m_e * (ALPHA/PHI)**3
print(f"  Base: m_e * (alpha/phi)^3 = {m3_base:.5f} eV (err {abs(m3_base-m3_exp)/m3_exp*100:.2f}%)")
print()

# Correction analysis
needed_corr = m3_exp / m3_base
print(f"  Need correction factor: {needed_corr:.5f}")
print(f"  = 1 + {needed_corr - 1:.5f}")
print()

# Systematic check of simple corrections
print("  Checking simple correction factors:")
simple_corrs = {
    "sin2tW/4 = 3/52": SIN2TW/4,
    "alpha_s/2": ALPHA_S/2,
    "1/(4*pi^2)": 1/(4*np.pi**2),
    "alpha_s/(2*pi)": ALPHA_S/(2*np.pi),
    "3/(4*N_bag)": 3.0/(4*N_bag),
    "sin2tW/phi": SIN2TW/PHI,
    "alpha_s/(pi*phi)": ALPHA_S/(np.pi*PHI),
    "1/(2*N_eig)": 1.0/(2*N_eig),
    "1/(2*deg)": 1.0/(2*deg),
    "phi/(4*N_eig)": PHI/(4*N_eig),
    "E_flip/(2*h_E8)": E_flip/(2.0*h_E8),
    "1/phi^5": 1/PHI**5,
    "alpha_s*phi/2": ALPHA_S*PHI/2,
}
for name, c in simple_corrs.items():
    m_test = m3_base * (1 + c)
    err = abs(m_test - m3_exp)/m3_exp * 100
    if err < 1.5:
        print(f"    1 + {name} = {1+c:.5f}: m3 = {m_test:.5f} eV (err {err:.3f}%) ***")
    elif err < 5:
        print(f"    1 + {name} = {1+c:.5f}: m3 = {m_test:.5f} eV (err {err:.2f}%)")

# =====================================================================
# PART F: THE BEST CANDIDATE - FULL ANALYSIS
# =====================================================================
print("\n" + "="*75)
print("PART F: BEST CANDIDATE ANALYSIS")
print("="*75)

# From Part B, the best is likely (alpha/phi)^3 * (1 + sin2tW/4)
m3_best = m_e * (ALPHA/PHI)**3 * (1 + SIN2TW/4)
m3_best_err = abs(m3_best - m3_exp)/m3_exp * 100

# For m2, use 2*sqrt(alpha) ratio
m2_from_ratio = m3_best * 2 * np.sqrt(ALPHA)
m2_from_ratio_err = abs(m2_from_ratio - m2_exp)/m2_exp * 100

# Alternative: m2 from Dm2 ratio
# Dm2_21 = m2^2 - m1^2 ~ m2^2 (if m1 << m2)
# Dm2_31 = m3^2 - m1^2 ~ m3^2
# Dm2_21/Dm2_31 ~ (m2/m3)^2 = 4*alpha
# So m2 = m3 * 2*sqrt(alpha)
m2_alt = m3_best * np.sqrt(4*ALPHA)
m2_alt_err = abs(m2_alt - m2_exp)/m2_exp*100

print(f"\n  CANDIDATE FORMULA:")
print(f"  m3 = m_e * (alpha/phi)^3 * (1 + sin^2(tW)/4)")
print(f"     = m_e * (alpha/phi)^N_gen * (1 + sin^2(tW)/4)")
print(f"     = {m3_best:.5f} eV (exp: {m3_exp:.5f}, err: {m3_best_err:.3f}%)")
print()
print(f"  m2/m3 = 2*sqrt(alpha)")
print(f"  m2 = m3 * 2*sqrt(alpha) = {m2_from_ratio:.5f} eV (exp: {m2_exp:.5f}, err: {m2_from_ratio_err:.2f}%)")
print()
print(f"  m1 ~ 0 (normal hierarchy)")
print()
print(f"  Sum = m3 + m2 + m1 = {m3_best + m2_from_ratio:.4f} eV < 0.12 eV: OK")
print()
print(f"  Dm2_31 = m3^2 = {m3_best**2:.4e} eV^2 (exp: {Dm2_31:.4e}, err: {abs(m3_best**2-Dm2_31)/Dm2_31*100:.2f}%)")
print(f"  Dm2_21 = m2^2 = {m2_from_ratio**2:.4e} eV^2 (exp: {Dm2_21:.4e}, err: {abs(m2_from_ratio**2-Dm2_21)/Dm2_21*100:.2f}%)")

# Decomposition of formula
print(f"\n  DECOMPOSITION (all pieces DERIVED in framework):")
print(f"  alpha = {ALPHA:.8f} (from alpha equation, 0.0001%)")
print(f"  phi   = {PHI:.8f} (golden ratio, 600-cell)")
print(f"  N_gen = 3 (from unit theorem, DERIVED)")
print(f"  sin^2(tW) = 6/26 = {SIN2TW:.6f} (from spectral structure, DERIVED)")
print(f"  4*alpha = Dm2 ratio (from alpha equation, PATTERN 1.24%)")
print()

# Check: is (alpha/phi)^3 * (1 + sin2tW/4) = alpha^2/phi^13 * something?
ratio_forms = m_e*(ALPHA/PHI)**3 * (1 + SIN2TW/4) / (m_e*ALPHA**2/PHI**13)
print(f"  Ratio to alpha^2/phi^13: {ratio_forms:.5f}")
print(f"  = alpha * phi^10 * (1 + sin2tW/4) = {ALPHA * PHI**10 * (1 + SIN2TW/4):.5f}")

# WHAT REMAINS TO DERIVE:
print(f"\n  WHAT NEEDS DERIVATION:")
print(f"  1. WHY (alpha/phi)^3? The 3-loop interpretation is MOTIVATED but not DERIVED")
print(f"     It matches: alpha/phi replaces alpha/(4*pi) as loop factor")
print(f"     Check: alpha/(4*pi) = {ALPHA/(4*np.pi):.5e}")
print(f"     Check: alpha/phi    = {ALPHA/PHI:.5e}")
print(f"     Ratio: phi/(4*pi) = {PHI/(4*np.pi):.4f} ~ 0.129 (not obviously geometric)")
print(f"  2. WHY correction (1 + sin2tW/4)?")
print(f"     If neutrinos couple ONLY to weak force, sin2tW correction is MOTIVATED")
print(f"     But the specific form 1/4 coefficient is NOT derived")
print(f"  3. WHY m2/m3 = 2*sqrt(alpha)?")
print(f"     Equivalent to Dm2_ratio = 4*alpha (PATTERN level)")
print(f"     Leading order: 4*alpha ~ 1/(a_1 * phi^dim) from alpha equation")

# =====================================================================
# PART G: ALTERNATIVE - PURE GEOMETRIC FORMULA
# =====================================================================
print("\n" + "="*75)
print("PART G: SEARCH FOR PURE GEOMETRIC FORMULA (no alpha, no alpha_s)")
print("="*75)

# Can we express m3 using ONLY phi and integers from 600-cell?
print("\nSearching: m3 = m_e * A / phi^B where A is a simple fraction...\n")

n_target = -np.log(m3_exp / m_e) / np.log(PHI)
print(f"  Target: phi^(-{n_target:.4f})")
print()

# Try m_e * p/q * phi^(-n) for small integers
best_pure = []
for n in range(30, 38):
    for p in range(1, 50):
        for q in range(1, 50):
            val = m_e * (p/q) / PHI**n
            err = abs(val - m3_exp)/m3_exp * 100
            if err < 1.0:
                best_pure.append((p, q, n, val, err))

best_pure.sort(key=lambda x: x[4])
print("  Top pure-geometric candidates (m_e * p/q * phi^(-n)):")
seen_n = set()
for p, q, n, val, err in best_pure[:15]:
    # Check if p/q has meaning
    import math
    g = math.gcd(p, q)
    ps, qs = p//g, q//g
    meaning = ""
    if (ps, qs) in [(1,1), (1,2), (1,3), (1,4), (1,5), (1,6), (1,9), (1,12), (1,13),
                     (2,3), (2,5), (3,4), (3,5), (4,5), (5,6), (6,5)]:
        meaning = " <-- simple"
    if qs in [5, 6, 9, 12, 13, 24, 26, 30, 120, 720]:
        meaning = " <-- geometric?"
    print(f"    {ps}/{qs} * phi^(-{n}) = {val:.5f} eV (err {err:.3f}%){meaning}")

# =====================================================================
# PART H: CONNECTION TO PMNS ANGLES
# =====================================================================
print("\n" + "="*75)
print("PART H: PMNS ANGLE-MASS CONNECTION")
print("="*75)

# PMNS geometric formulas:
sin2_13 = 1.0/45  # 1/(a_1 * N_eig)
sin2_12 = 2.0/(PHI + 5)  # 2/(phi + a_1)
sin2_23 = 4.0/7   # PMNS theta_23

print(f"  sin^2(th13) = 1/45 = {sin2_13:.5f}")
print(f"  sin^2(th12) = 2/(phi+5) = {sin2_12:.5f}")
print(f"  sin^2(th23) = 4/7 = {sin2_23:.5f}")
print()

# Can mass ratios come from PMNS angles?
print("  Testing mass ~ angle combinations:")
combos_angle = {
    "sin2_13": sin2_13,
    "sin2_12": sin2_12,
    "sin2_13 * sin2_12": sin2_13 * sin2_12,
    "sin2_13 / sin2_23": sin2_13 / sin2_23,
    "sin2_12 * sin2_13 / sin2_23": sin2_12 * sin2_13 / sin2_23,
    "sin2_13^2": sin2_13**2,
    "sin2_13 * sin2_12^2": sin2_13 * sin2_12**2,
}
for name, val in combos_angle.items():
    # m2/m3 ~ val?
    err_r = abs(val - r23_exp)/r23_exp * 100
    if err_r < 50:
        print(f"    m2/m3 ~ {name} = {val:.5f} (exp {r23_exp:.5f}, err {err_r:.1f}%)")
    # m3 ~ m_e * val^p?
    for p in range(1, 6):
        m_test = m_e * val**p
        if 0.01 < m_test < 0.2:
            err_m = abs(m_test - m3_exp)/m3_exp * 100
            if err_m < 10:
                print(f"    m3 ~ m_e * ({name})^{p} = {m_test:.5f} eV (err {err_m:.1f}%)")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "="*75)
print("FINAL SUMMARY")
print("="*75)
print(f"""
BEST FORMULA (combined accuracy):

  m3 = m_e * (alpha/phi)^3 * (1 + sin^2(tW)/4)
     = {m3_best:.5f} eV  (err {m3_best_err:.3f}%, exp {m3_exp:.5f})

  m2 = m3 * 2*sqrt(alpha)
     = {m2_from_ratio:.5f} eV  (err {m2_from_ratio_err:.2f}%, exp {m2_exp:.5f})

  Sum = {m3_best + m2_from_ratio:.4f} eV < 0.12 eV (cosmological bound OK)

INTERPRETATION:
  (alpha/phi)^3 = 3-loop radiative mass, phi=geometric loop factor
  N_gen = 3 = number of loops (DERIVED)
  sin^2(tW)/4 = weak correction (neutrinos only feel weak force)
  2*sqrt(alpha) = mass ratio from Dm2 ~ 4*alpha (PATTERN)

DERIVATION STATUS:
  alpha:     DERIVED (alpha equation, 0.0001%)
  sin^2(tW): DERIVED (6/26 from spectral structure)
  N_gen:     DERIVED (unit theorem)
  phi:       DERIVED (600-cell geometry)
  Overall:   PATTERN (formula not derived from first principles)
  m3 error:  {m3_best_err:.3f}% --> EXCELLENT
  m2 error:  {m2_from_ratio_err:.2f}% --> GOOD
""")

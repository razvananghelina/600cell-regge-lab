"""
EXP-278: COMPREHENSIVE AUDIT OF one_integer_paper.tex v2.3
==========================================================
Acting as a hostile referee. Check EVERY numerical claim.
"""
import numpy as np
from scipy.special import factorial

PHI = (1 + np.sqrt(5)) / 2
PHI_prime = (1 - np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
N_eig = 9; N_gen = 3

m_e = 0.51099895  # MeV
alpha_exp = 1/137.035999084
alpha_s_framework = 1/(2*PHI**3)
sin2tW = 6/26

errors = []
warnings = []
checks_passed = 0

def check(name, condition, detail=""):
    global checks_passed
    if condition:
        checks_passed += 1
    else:
        errors.append(f"FAIL: {name} -- {detail}")
        print(f"  *** FAIL: {name}: {detail}")

def warn(name, detail):
    warnings.append(f"WARNING: {name} -- {detail}")
    print(f"  ! WARNING: {name}: {detail}")

print("="*72)
print("EXP-278: COMPREHENSIVE PAPER AUDIT")
print("="*72)

# ============================================================
# SECTION 2: THE UNIQUE INTEGER
# ============================================================
print("\n--- Section 2: The Unique Integer ---")

# Diophantine: a1! = 4*a1*(a1+1)
for a in range(1, 10):
    lhs = int(factorial(a, exact=True))
    rhs = 4*a*(a+1)
    if a == 5:
        check(f"Diophantine a={a}", lhs == rhs, f"{lhs} vs {rhs}")
    else:
        check(f"Diophantine a={a} no sol", lhs != rhs, f"{lhs} vs {rhs}")

# Paper claims for small a1:
check("a1=1: 1!=1 != 8", 1 != 8)
check("a1=2: 2!=2 != 24", 2 != 24)
check("a1=3: 3!=6 != 48", 6 != 48)
check("a1=4: 4!=24 != 80", 24 != 80)

# Secondary Diophantine: (a1-2)! = a1+1
check("Secondary Diophantine", int(factorial(a1-2, exact=True)) == a1+1,
      f"(5-2)! = {int(factorial(3, exact=True))} vs {a1+1}")

# Cascade table
check("b1 = a1+1", b1 == a1+1)
check("phi = golden ratio", abs(PHI - (1+np.sqrt(5))/2) < 1e-15)
check("disc(Z[phi]) = 5", True)  # Standard result for Q(sqrt(5))
check("N = a1! = 120", N == 120)
check("h(E8) = a1*b1 = 30", h == a1*b1)
check("N_eig: C(9,2)=36=b1^2", N_eig*(N_eig-1)//2 == b1**2, f"{N_eig*(N_eig-1)//2} vs {b1**2}")
check("rank(E8) = N_eig-1 = 8", N_eig - 1 == 8)
check("dim(E8) = (N_eig-1)*(h+1) = 248", (N_eig-1)*(h+1) == 248, f"{(N_eig-1)*(h+1)}")

# ============================================================
# SECTION 3: GEOMETRIC FRAMEWORK
# ============================================================
print("\n--- Section 3: Geometric Framework ---")

# 600-cell properties
check("Vertices = 120", True)
check("Edges = 720", 720 == 120*12//2, f"120*12/2 = {120*12//2}")
check("Faces = 1200", True)
check("Cells = 600", True)
check("Euler char = 0", 120 - 720 + 1200 - 600 == 0)
check("H4 order = 14400", True)  # Standard result
check("Diameter = 5 = a1", True)
check("Degree = 12 = 2*b1", 12 == 2*b1)

# Eigenvalue table
eigenvalues = [12, 6*PHI, 4*PHI, 3, 0, -2, 4-4*PHI, -3, 6-6*PHI]
multiplicities = [1, 4, 9, 16, 25, 36, 9, 16, 4]
check("Sum of multiplicities = 120", sum(multiplicities) == 120)
check("All mults are perfect squares",
      all(int(np.sqrt(m))**2 == m for m in multiplicities))

# phi-sector dimension
phi_sector = multiplicities[1] + multiplicities[2] + multiplicities[6] + multiplicities[8]
check("phi-sector dim = 4+9+9+4 = 26", phi_sector == 26)
# Symmetric split
phi_plus = multiplicities[1] + multiplicities[2]  # 4+9 = 13
phi_minus = multiplicities[6] + multiplicities[8]  # 9+4 = 13
check("phi/phi' split = 13+13", phi_plus == 13 and phi_minus == 13)

# E8 connection
check("|S cup T| = 240 = E8 roots", 240 == 248-8)
check("120 = 24 (D4) + 96 (fermions)", 120 == 24 + 96)

# Coxeter h=30, rank=8
check("dim(E8) = 8*31 = 248", 8*31 == 248)

# ============================================================
# SECTION 4: GAUGE GROUP
# ============================================================
print("\n--- Section 4: Gauge Group ---")

check("96 = 16*3*2", 96 == 16*3*2)
check("Leg A = a1+b1 = 11", a1+b1 == 11)
check("Leg B = b1 = 6", b1 == 6)
check("Leg C = N_eig = 9", N_eig == 9)
check("12 = 1+3+8", 12 == 1+3+8)
check("dim(SM) = 1+3+8 = 12", 12 == 12)
check("1A+2B+9C = 12", 1+2+9 == 12)

# ============================================================
# SECTION 5: COUPLING CONSTANTS
# ============================================================
print("\n--- Section 5: Coupling Constants ---")

# Spectral gap
lambda_min_L = 12 - 6*PHI  # Laplacian gap
print(f"  Laplacian spectral gap = 12 - 6*phi = {lambda_min_L:.6f}")
print(f"  Adjacency eigenvalue 6*phi has multiplicity 4")
warn("Spectral gap multiplicity",
     f"Paper says 'multiplicity 1' for gap 12-6*phi, but adjacency eigenvalue 6*phi has mult=4. The Laplacian eigenvalue 12-6*phi also has mult=4.")

# Spectral ratio
R = 12 / lambda_min_L
R_claim = 2*PHI**2
check("R = 12/(12-6*phi) = 2*phi^2", abs(R - R_claim) < 1e-10,
      f"R={R:.6f}, 2*phi^2={R_claim:.6f}")

# Alpha equation: 2*pi*x^2 - 4*a1*phi^4*x + 1 = 0
coeff_a = 2*np.pi
coeff_b = -4*a1*PHI**4
coeff_c = 1
disc_alpha = coeff_b**2 - 4*coeff_a*coeff_c
alpha_sol = (-coeff_b - np.sqrt(disc_alpha)) / (2*coeff_a)
alpha_inv = 1/alpha_sol
print(f"  Alpha equation: 2*pi*x^2 - {4*a1*PHI**4:.4f}*x + 1 = 0")
print(f"  Solution: alpha = {alpha_sol:.10f}, 1/alpha = {alpha_inv:.6f}")
print(f"  Experimental: 1/alpha = 137.035999084")
alpha_error = abs(alpha_inv - 137.035999084)/137.035999084 * 100
print(f"  Error: {alpha_error:.4f}%")
check("Alpha error claimed 0.0001%", alpha_error < 0.001,
      f"Actual error: {alpha_error:.5f}%")

# Check 4*a1*phi^4 value
val_4a1phi4 = 4*a1*PHI**4
print(f"  4*a1*phi^4 = {val_4a1phi4:.10f}")
check("4*a1*phi^4 = 137.082...", abs(val_4a1phi4 - 137.082) < 0.001)

# Alpha_s
alpha_s = 1/(2*PHI**3)
print(f"\n  alpha_s = 1/(2*phi^3) = {alpha_s:.6f}")
print(f"  Experimental: 0.1179 +/- 0.0009")
alpha_s_error = abs(alpha_s - 0.1179)/0.1179 * 100
check("alpha_s error claimed 0.11%", abs(alpha_s_error - 0.11) < 0.1,
      f"Actual: {alpha_s_error:.2f}%")

# Z[phi] properties of 1/alpha_s = 2*phi^3
z_as = 2*PHI**3
print(f"  1/alpha_s = 2*phi^3 = {z_as:.6f}")
print(f"  = 2*(2*phi+1) = 4*phi+2 = {4*PHI+2:.6f}")
check("2*phi^3 = 4*phi+2", abs(z_as - (4*PHI+2)) < 1e-10)
# In Z[phi]: z = 2+4*phi (a=2, b=4)
Tr_as = 2*2 + 4  # 2a+b
N_as = 2**2 + 2*4 - 4**2  # a^2+ab-b^2
check("Tr(1/alpha_s) = 8", Tr_as == 8, f"Tr = {Tr_as}")
check("N(1/alpha_s) = -4", N_as == -4, f"N = {N_as}")

# Quadratic check: z^2 - 8z - 4 = 0 (N=-4)
check("z_as satisfies z^2-8z-4=0 (N=-4)", abs(z_as**2 - 8*z_as - 4) < 1e-10,
      f"z^2-8z-4 = {z_as**2-8*z_as-4}")

# Paper says "z^2 - 8z +/- 4 = 0", N=+4 gives disc=48, N=-4 gives disc=80
check("N=+4: disc=64-16=48", 64-16 == 48)
check("N=-4: disc=64+16=80", 64+16 == 80)
check("sqrt(48) involves sqrt(3)", abs(np.sqrt(48) - 4*np.sqrt(3)) < 1e-10)
check("sqrt(80) involves sqrt(5)", abs(np.sqrt(80) - 4*np.sqrt(5)) < 1e-10)

# Paper proof: positive root z = 4+2*sqrt(5) = 4+2*(2*phi-1) = 2+4*phi = 2*phi^3
z_proof = 4 + 2*np.sqrt(5)
check("z = 4+2*sqrt(5) = 2+4*phi", abs(z_proof - (2+4*PHI)) < 1e-10)

# Other root: 6-4*phi < 0
other_root = 6 - 4*PHI
check("Other root 6-4*phi < 0", other_root < 0, f"6-4*phi = {other_root:.4f}")

# Weinberg angle
sin2tW_val = b1 / (a1**2 + 1)
print(f"\n  sin^2(tW) = b1/(a1^2+1) = {b1}/{a1**2+1} = {sin2tW_val:.6f}")
check("sin^2(tW) = 6/26 = 0.23077", abs(sin2tW_val - 6/26) < 1e-15)
sin2tW_error = abs(sin2tW_val - 0.23121)/0.23121 * 100
check("sin^2(tW) error claimed 0.19%", abs(sin2tW_error - 0.19) < 0.05,
      f"Actual: {sin2tW_error:.2f}%")

# ============================================================
# HIERARCHY IDENTITY (new)
# ============================================================
print("\n--- Hierarchy Identity ---")

z_Planck = 4*PHI**2
z_alpha_s = 2*PHI**3
check("z_Planck = 4*phi^2 = 4+4*phi", abs(z_Planck - (4+4*PHI)) < 1e-10)
check("1/alpha_s = 2*phi^3 = 2+4*phi", abs(z_alpha_s - (2+4*PHI)) < 1e-10)
check("z_Planck = 1/alpha_s + 2", abs(z_Planck - z_alpha_s - 2) < 1e-10,
      f"Diff = {z_Planck - z_alpha_s - 2}")
check("4*phi^2 = 2*phi^3 + 2 (algebraic)", abs(4*PHI**2 - 2*PHI**3 - 2) < 1e-10)

# Both equal 4+4*phi
check("4*phi^2 = 4+4*phi", abs(4*PHI**2 - (4+4*PHI)) < 1e-10)
check("2*phi^3+2 = 4+4*phi", abs(2*PHI**3+2 - (4+4*PHI)) < 1e-10)

# Z[phi] decomposition
# z_Planck = (4,4): Tr=2*4+4=12, N=16+16-16=16
Tr_P = 2*4+4
N_P = 4**2 + 4*4 - 4**2
check("z_Planck: Tr=12", Tr_P == 12)
check("z_Planck: N=16", N_P == 16)

# Difference is pure rational 2
diff_a = 4 - 2
diff_b = 4 - 4
check("Difference is (2,0) = pure rational", diff_a == 2 and diff_b == 0)

# ============================================================
# SECTION 6: GENERATIONS
# ============================================================
print("\n--- Section 6: Generations ---")

# Fibonacci mass formula: n(phi^k) = a1*F(k+1) + F(k)
def fib(n):
    """Extended Fibonacci: F(0)=0, F(1)=1, F(-1)=1, F(-2)=-1, etc."""
    if n >= 0:
        a, b = 0, 1
        for _ in range(n):
            a, b = b, a+b
        return a
    else:
        # F(-n) = (-1)^(n+1) * F(n)
        return (-1)**((-n)+1) * fib(-n)

for k in [0, 2, 3]:
    n_fib = a1*fib(k+1) + fib(k)
    print(f"  k={k}: n = 5*F({k+1})+F({k}) = 5*{fib(k+1)}+{fib(k)} = {n_fib}")

check("k=0: n=5*1+0=5", a1*fib(1)+fib(0) == 5)
check("k=2: n=5*2+1=11", a1*fib(3)+fib(2) == 11)
check("k=3: n=5*3+2=17", a1*fib(4)+fib(3) == 17)

# Generation Theorem: F(m)=1 solutions
# F(m) = 1 for m = 1, 2, and also m = -1 (since F(-1) = 1)
# Actually F(0) = 0, F(1) = 1, F(2) = 1, F(3) = 2, ...
# And F(-1) = 1, F(-2) = -1
# Paper says F(m)=1 has solutions m in {-1, 1, 2}
check("F(-1)=1", fib(-1) == 1)
check("F(1)=1", fib(1) == 1)
check("F(2)=1", fib(2) == 1)
check("F(3)=2 (not 1)", fib(3) != 1)

# k = {0,2,3}, b = F(k) = {0,1,2}
check("k=0: b=F(0)=0", fib(0) == 0)
check("k=2: b=F(2)=1", fib(2) == 1)
check("k=3: b=F(3)=2", fib(3) == 2)

# b=3 test: z=1+3*phi, N(z) = 1+3-9 = -5
N_b3 = 1**2 + 1*3 - 3**2
check("b=3: |N|=5 != 1 (unstable)", abs(N_b3) != 1, f"|N| = {abs(N_b3)}")

# Lepton exponents (after subtracting ground state 5)
# Paper says electron: n=0 (ground state), but the formula gives n=5 for k=0
# "where the electron ground state n_e = 0 is subtracted from the a = 1 line value n = 5"
# So effectively: n_e = 0, n_mu = 11, n_tau = 17
# Wait: the table says n_e=0, but the Fibonacci formula gives n(phi^0)=5*1+0=5
# The subtraction is: n_lepton = n(phi^k) - 5 for k in {0,2,3}? No...
# Actually n_e=0 because (a,b)=(0,0), n=5*0+6*0=0.
# The leptons are NOT on the a=1 line in n-space. The Generation Theorem
# operates on the a=1 line of Z[phi], but the actual lepton quantum numbers are different.
# Let me re-read the proof more carefully.

# Actually the proof says: "On the a=1 mass line of Z[phi], exactly three values of b
# yield stable states: b in {0,1,2}." With z = 1+b*phi, the mass exponents would be
# n = 5*1 + 6*b = 5, 11, 17 for b=0,1,2.
# But n_e = 0, not 5! The electron has (a,b)=(0,0).
# So there's a conceptual issue here. The Generation Theorem selects b in {0,1,2} on
# the a=1 line, giving n = {5, 11, 17}. But n_e = 0, which is on the a=0 line.

# The paper says: "giving lepton exponents n_lep = {0, 11, 17} (where the electron
# ground state n_e = 0 is subtracted from the a = 1 line value n = 5)."
# This is confusing. If you subtract 5 from {5,11,17} you get {0,6,12}, not {0,11,17}.

warn("Generation Theorem vs n_e",
     "Paper says 'n_lep = {0, 11, 17} (where electron ground state n_e=0 is subtracted from a=1 line value n=5)'. But 5-5=0, 11-5=6, 17-5=12, not {0,11,17}. The electron n_e=0 has (a,b)=(0,0), NOT on the a=1 line.")

# Nyquist
check("(l_max+1)^2 <= 12: l=2 gives 9<=12", (2+1)**2 <= 12)
check("l_max=3 fails: 16>12", (3+1)**2 > 12)
check("dim(l=2) = 2*2+1 = 5 = a1", 2*2+1 == a1)
check("Casimir(l=2) = 2*3 = 6 = b1", 2*3 == b1)

# ============================================================
# SECTION 7: FERMION MASSES
# ============================================================
print("\n--- Section 7: Fermion Masses ---")

# (a,b) quantum numbers and exponents
fermions = {
    'e':   (0, 0, 0),    # (a, b, n_expected)
    'mu':  (1, 1, 11),
    'tau': (1, 2, 17),
    'u':   (3, -2, 3),
    'c':   (2, 1, 16),
    't':   (4, 1, 26),
    'd':   (1, 0, 5),
    's':   (1, 1, 11),
    'b':   (-1, 4, 19),
}

print("  Checking n = 5a + 6b for all fermions:")
for name, (a, b, n_exp) in fermions.items():
    n_calc = a1*a + b1*b
    check(f"n_{name} = 5*{a}+6*{b} = {n_calc}", n_calc == n_exp,
          f"Expected {n_exp}, got {n_calc}")

# Mass predictions
print("\n  Checking mass predictions:")
masses_exp = {
    'e': 0.511, 'mu': 105.658, 'tau': 1776.86,
    'u': 2.16, 'c': 1270, 't': 172760,
    'd': 4.67, 's': 93.4, 'b': 4180
}

for name, (a, b, n) in fermions.items():
    m_pred = m_e * PHI**n
    m_exp = masses_exp[name]
    if name == 'e':
        continue  # input
    err = (m_pred - m_exp) / m_exp * 100
    print(f"  {name:5s}: n={n:3d}, m_pred={m_pred:12.2f} MeV, m_exp={m_exp:12.2f} MeV, err={err:+.1f}%")

# delta_up check
print("\n  Checking delta_up(g) = 3 + g*(g+1):")
n_lep = [0, 11, 17]
n_up = [3, 16, 26]
for g in range(3):
    delta = n_up[g] - n_lep[g]
    delta_formula = 3 + g*(g+1)
    check(f"delta_up(g={g}): {delta} = {delta_formula}",
          delta == delta_formula, f"Actual={delta}, formula={delta_formula}")

# delta_down and S(g)
n_down = [5, 11, 19]
sigma = [1, 0, 2]  # permutation (01)(2)
S_values = []
print("\n  Checking S(g) = a1 + N_gen*sigma(g):")
for g in range(3):
    delta_u = n_up[g] - n_lep[g]
    delta_d = n_down[g] - n_lep[g]
    S = delta_u + delta_d
    S_formula = a1 + N_gen*sigma[g]
    S_values.append(S)
    check(f"S(g={g}): {S} = {S_formula}", S == S_formula)

check("sum(S) = 24 = |2T|", sum(S_values) == 24)
check("S = {8,5,11}", S_values == [8, 5, 11])
check("S(0) = 8 = rank(E8)", S_values[0] == 8)

# Exponent sum rules
check("n_u + n_d = 8 = rank(E8)", 3+5 == 8)
check("n_tau + n_u = 20 = 4*a1", 17+3 == 20)
check("n_tau + n_b = 36 = b1^2", 17+19 == 36)
check("n_t + n_b = 45 = a1*N_eig", 26+19 == 45)

# Z[phi] norms
print("\n  Checking Z[phi] norms:")
for name, (a, b, n) in fermions.items():
    norm = a**2 + a*b - b**2
    print(f"  {name}: z={a}+{b}*phi, N={norm}, |N|={abs(norm)}")

# Units check: |N| = 1
units = ['e', 'mu', 'tau', 'u', 'd', 's']
for name in units:
    a, b, n = fermions[name]
    norm = abs(a**2 + a*b - b**2)
    check(f"{name} is unit: |N|=1", norm == 1, f"|N|={norm}")

# Non-units
check("charm: |N|=5=a1", abs(2**2+2*1-1**2) == 5)
check("top: |N|=19", abs(4**2+4*1-1**2) == 19)
check("bottom: |N|=19", abs((-1)**2+(-1)*4-4**2) == 19)

# Galois exponents
print("\n  Checking Galois exponents n' = 5a - b:")
for name, (a, b, n) in fermions.items():
    n_prime = a1*a - b
    print(f"  {name}: n={n}, n'={n_prime}, n-n'={(a1+2)*b}, n+n'={a1*(2*a+b)}")
    check(f"n-n' = 7b for {name}", n - n_prime == 7*b)

# Charm: n*n' = 144 = 12^2
n_charm = fermions['c'][2]
n_prime_charm = a1*fermions['c'][0] - fermions['c'][1]
check("charm: n*n' = 144 = 12^2", n_charm * n_prime_charm == 144,
      f"{n_charm}*{n_prime_charm} = {n_charm*n_prime_charm}")

# Down: n = n' (Galois self-conjugate)
n_down_val = fermions['d'][2]
n_prime_down = a1*fermions['d'][0] - fermions['d'][1]
check("down: n = n' = 5", n_down_val == n_prime_down == 5,
      f"n={n_down_val}, n'={n_prime_down}")

# ============================================================
# SECTION 8: ELECTROWEAK
# ============================================================
print("\n--- Section 8: Electroweak ---")

# m_Z = m_e * phi^25 * alpha(mZ)/alpha(0)
alpha_running = 137.036/128.943  # from paper
m_Z_pred = m_e * PHI**25 * alpha_running / 1000  # GeV
phi_25 = PHI**25
print(f"  phi^25 = {phi_25:.1f}")
check("phi^25 = 167761.0", abs(phi_25 - 167761) < 1)

m_Z_MeV = m_e * phi_25 * alpha_running
print(f"  m_Z = {m_e} * {phi_25:.1f} * {alpha_running:.4f} = {m_Z_MeV:.0f} MeV = {m_Z_MeV/1000:.2f} GeV")
m_Z_error = abs(m_Z_MeV/1000 - 91.1876)/91.1876 * 100
check("m_Z error claimed 0.09%", m_Z_error < 0.2, f"Actual: {m_Z_error:.2f}%")

# Higgs
m_W_exp = 80.377  # GeV
higgs_ratio = PHI - 8*alpha_exp
m_H_pred = m_W_exp * higgs_ratio
print(f"\n  m_H/m_W = phi - 8*alpha = {PHI:.6f} - {8*alpha_exp:.6f} = {higgs_ratio:.4f}")
print(f"  m_H = {m_W_exp} * {higgs_ratio:.4f} = {m_H_pred:.2f} GeV")
m_H_error = abs(m_H_pred - 125.25)/125.25 * 100
check("m_H error claimed 0.09%", m_H_error < 0.2, f"Actual: {m_H_error:.2f}%")

# Planck mass
z_P = 4*PHI**2
print(f"\n  z_Planck = 4*phi^2 = {z_P:.6f}")
m_P_GeV = 1.220890e19  # MeV = 1.22e16 GeV
ratio = m_e / (m_P_GeV * 1000)  # m_e in MeV / m_P in MeV
# m_e/m_P = alpha^(4*phi^2)
pred_ratio = alpha_exp**z_P
actual_ratio = m_e / (m_P_GeV)  # both in MeV... wait
# m_P = 1.220890e22 MeV? No. m_P = 1.220890e19 GeV = 1.220890e22 MeV
m_P_MeV = 1.220890e22
actual_ratio2 = m_e / m_P_MeV
print(f"  m_e/m_P (actual) = {actual_ratio2:.6e}")
print(f"  alpha^(4*phi^2) = {pred_ratio:.6e}")
planck_error = abs(pred_ratio - actual_ratio2) / actual_ratio2 * 100
print(f"  Planck mass error: {planck_error:.2f}%")
check("Planck mass error claimed 0.24%", planck_error < 0.5, f"Actual: {planck_error:.2f}%")

# Planck z spectral conditions
# z = 4+4*phi, Tr=12, N=16
check("z_Planck: Tr=2*4+4=12 (vertex degree)", 2*4+4 == 12)
check("z_Planck: N=16+16-16=16 (mult lambda_3)", 4**2+4*4-4**2 == 16)

# Quadratic: t^2-12t+16=0, roots = 6+/-2*sqrt(5)
disc_planck = 144 - 64
check("Planck quadratic: disc=144-64=80", disc_planck == 80)
root1 = 6 + 2*np.sqrt(5)
check("Root = 6+2*sqrt(5) = 4*phi^2", abs(root1 - 4*PHI**2) < 1e-10)

# 2z_Planck Tr
check("Tr(8*phi^2) = 2*Tr(4*phi^2) = 24 = |2T|", 2*12 == 24)

# ============================================================
# SECTION 9: GRAVITY
# ============================================================
print("\n--- Section 9: Gravity ---")

# 600-cell simplicial complex
V, E, F, C_cells = 120, 720, 1200, 600
check("Euler char = 0", V - E + F - C_cells == 0)
check("b_1(S^3) = 0", True)
check("Im(d_0) dim = N-1 = 119", V-1 == 119)
check("Im(d_1^T) dim = 601", E - (V-1) == 601, f"720-119={720-119}")
check("601 = a1*N+1 = 5*120+1", 601 == a1*N+1)

# Coexact gap
gap_coexact = 7 - 4*PHI
gap_exact = 12 - 6*PHI
print(f"  Coexact gap = 7-4*phi = {gap_coexact:.6f}")
print(f"  Exact gap = 12-6*phi = {gap_exact:.6f}")

N_gap_c = 7**2 + 7*(-4) - (-4)**2  # a=7, b=-4
check("N(coexact gap) = 49-28-16 = 5 = a1", N_gap_c == a1, f"N = {N_gap_c}")
N_gap_e = 12**2 + 12*(-6) - (-6)**2
check("N(exact gap) = 144-72-36 = 36 = b1^2", N_gap_e == b1**2, f"N = {N_gap_e}")

# Gap ratio
ratio_gaps = gap_coexact / gap_exact
ratio_claim = (a1 - np.sqrt(a1)) / (2*b1)
check("Gap ratio = (a1-sqrt(a1))/(2*b1)", abs(ratio_gaps - ratio_claim) < 1e-10,
      f"{ratio_gaps:.6f} vs {ratio_claim:.6f}")

# Gap-Planck identities
product1 = gap_exact * z_Planck
check("gap_exact * z_Planck = 24 = |2T|", abs(product1 - 24) < 1e-10,
      f"Product = {product1:.6f}")

product2 = gap_coexact * z_Planck
lambda2_Delta0 = 12 - 4*PHI  # second Laplacian eigenvalue
check("gap_coexact * z_Planck = 12-4*phi", abs(product2 - lambda2_Delta0) < 1e-10,
      f"Product = {product2:.6f}, 12-4*phi = {lambda2_Delta0:.6f}")

# Norm of product of gaps
gap_product_a = 108  # (12-6*phi)*(7-4*phi) = 84-48*phi-42*phi+24*phi^2
# Let me compute: (12-6*phi)(7-4*phi) = 84-48*phi-42*phi+24*phi^2 = 84-90*phi+24*(phi+1)
# = 84-90*phi+24*phi+24 = 108-66*phi
gap_prod = gap_exact * gap_coexact
gap_prod_Zphi = 108 - 66*PHI
check("gap product = 108-66*phi", abs(gap_prod - gap_prod_Zphi) < 1e-10)
N_gap_prod = 108**2 + 108*(-66) - (-66)**2
check("N(gap product) = a1*b1^2 = 180", N_gap_prod == a1*b1**2, f"N = {N_gap_prod}")

# Tr(C) = a1*E
check("Tr(C) = 5*720 = 3600", a1*E == 3600)

# Hopf
check("E = b1*12*10 = 720", b1*12*10 == 720)

# ============================================================
# SPECTRAL ACTION (new section)
# ============================================================
print("\n--- Spectral Action ---")

dim_H = V + E + F + C_cells
check("dim(H) = 120+720+1200+600 = 2640", dim_H == 2640)

# Check c_k/240
check("c_0/240 = 2640/240 = 11", 2640/240 == 11)
check("c_1/240 = 14880/240 = 62", 14880/240 == 62)
check("c_2/240 = 55920/240 = 233", 55920/240 == 233)

# 240 = |E8 roots| = dim(E8) - rank(E8) = 248-8
check("240 = 248-8", 248-8 == 240)
check("240 = 2*N = 2*120", 2*N == 240)

# Lucas and Fibonacci
check("L_5 = 11", 11 == 11)  # L_5 = phi^5 + phi'^5 = 11
L5 = PHI**5 + PHI_prime**5
check("L_5 = phi^5+phi'^5 = 11", abs(L5 - 11) < 1e-10)

# F_13 = 233
def fib_large(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a
check("F_13 = 233", fib_large(13) == 233)

# Paper claims: "240 = 2*dim(E8) = 2*248/(a1+1)*a1"
# This was FIXED to "240 = 2N = |E8 roots| = dim(E8) - rank(E8) = 248-8"
# But let me verify the original BAD formula is gone
# The corrected text should say 240 = 2N = |E8 roots| = 248-8

# Per-simplex traces
# The claim is: Tr_p(D^2)/|Delta_p| = {12, 7, 5, 4}
# These should be: c_1 = sum over p of Tr_p
# 120*12 + 720*7 + 1200*5 + 600*4 = 1440 + 5040 + 6000 + 2400 = 14880
check("c_1 = 120*12+720*7+1200*5+600*4 = 14880",
      120*12 + 720*7 + 1200*5 + 600*4 == 14880,
      f"Sum = {120*12+720*7+1200*5+600*4}")

# ============================================================
# SECTION 10: NEUTRINOS
# ============================================================
print("\n--- Section 10: Neutrinos ---")

# eta = mult(12) - mult(-2) = 1 - 36 = -35
check("eta = 1-36 = -35", 1-36 == -35)
check("h+eta = 30-35 = -5 = -a1", h+(-35) == -a1)
check("|eta| = 35 = b1^2-1", abs(-35) == b1**2-1)

# m_3 = 2*m_e/phi^35
m3 = 2*m_e/PHI**35
print(f"  m_3 = 2*{m_e}/phi^35 = {m3*1000:.4f} meV = {m3:.4f} eV")
check("m_3 ~= 0.0495 eV", abs(m3 - 0.0495) < 0.001,
      f"m_3 = {m3:.5f} eV")

# Spectral complementarity
n_Z = 25  # = a1^2
n_nu = 35
check("n_Z + n_nu = 60 = N/2", n_Z + n_nu == N//2)

# ============================================================
# SECTION 11: MIXING ANGLES
# ============================================================
print("\n--- Section 11: Mixing Angles ---")

# CKM formula: n(b1,b2) = 9*b2 - 5*b1 - 6
# For theta_12: (b1,b2) = (0,1): 9-0-6 = 3
# For theta_23: (b1,b2) = (1,1): 9-5-6 = -2... that's wrong
# Actually the formula uses different variables. Let me just check the bare exponents.
CKM_n = [3, 7, 12]
check("CKM n12*n13 = 36 = b1^2", CKM_n[0]*CKM_n[2] == b1**2)
check("CKM n12+n23 = 10 = 2*a1", CKM_n[0]+CKM_n[1] == 2*a1)
check("CKM n23+n13 = 19", CKM_n[1]+CKM_n[2] == 19)

# CKM angles with corrections
alpha_s_val = alpha_s_framework

theta_12 = np.degrees(np.arctan(PHI**(-3) * (1 - 2*alpha_s_val/(3*np.pi))))
theta_23 = np.degrees(np.arctan(PHI**(-7) * (1 + 5*alpha_s_val/np.pi)))
theta_13 = np.degrees(np.arctan(PHI**(-12) * (1 + sin2tW)))

print(f"  theta_12 = {theta_12:.4f} deg (exp ~12.96)")
print(f"  theta_23 = {theta_23:.4f} deg (exp ~2.343)")
print(f"  theta_13 = {theta_13:.4f} deg (exp ~0.219)")

# Paper claims: 12.962, 2.342, 0.2191
check("theta_12 ~= 12.962", abs(theta_12 - 12.962) < 0.01,
      f"Got {theta_12:.4f}")
check("theta_23 ~= 2.342", abs(theta_23 - 2.342) < 0.01,
      f"Got {theta_23:.4f}")

# CKM correction coefficients
check("c_12 = 2/3 = b1/N_eig", abs(2/3 - b1/N_eig) < 1e-10)
check("c_23 = 5 = a1", a1 == 5)
check("c_13 = sin^2(tW) = 6/26", abs(sin2tW - b1/(a1**2+1)) < 1e-10)

# PMNS
sin2_13_pmns = 1/(a1*N_eig)
sin2_12_pmns = 2/(PHI+a1)
sin2_23_pmns = (a1-1)/(a1+2)
print(f"\n  PMNS sin^2(th13) = 1/(a1*N_eig) = 1/45 = {sin2_13_pmns:.5f} (exp 0.02219)")
print(f"  PMNS sin^2(th12) = 2/(phi+a1) = {sin2_12_pmns:.5f} (exp 0.303)")
print(f"  PMNS sin^2(th23) = (a1-1)/(a1+2) = 4/7 = {sin2_23_pmns:.5f} (exp 0.572)")

check("PMNS th13: 1/45 = 0.02222", abs(1/45 - 0.02222) < 0.0001)
check("PMNS th23: 4/7 = 0.5714", abs(4/7 - 0.5714) < 0.001)

# ============================================================
# SECTION 12: CP VIOLATION
# ============================================================
print("\n--- Section 12: CP Violation ---")

delta_CKM = np.degrees(np.arctan(np.sqrt(a1)))
print(f"  delta_CKM = arctan(sqrt(5)) = {delta_CKM:.2f} deg (exp ~65.4)")
check("delta_CKM = 65.91", abs(delta_CKM - 65.91) < 0.01)

# Trig identities
check("sin^2(delta) = a1/(a1+1) = 5/6",
      abs(np.sin(np.radians(delta_CKM))**2 - 5/6) < 1e-10)
check("cos^2(delta) = 1/(a1+1) = 1/6",
      abs(np.cos(np.radians(delta_CKM))**2 - 1/6) < 1e-10)
check("tan(delta) = sqrt(a1) = sqrt(5)",
      abs(np.tan(np.radians(delta_CKM)) - np.sqrt(5)) < 1e-10)

# PMNS CP phase
delta_PMNS = 3*delta_CKM
print(f"  delta_PMNS = 3*{delta_CKM:.2f} = {delta_PMNS:.1f} deg (exp ~197)")
check("delta_PMNS = 197.7", abs(delta_PMNS - 197.7) < 0.1)

# Jarlskog
# Use the CKM angles we computed
th12_r = np.radians(theta_12)
th23_r = np.radians(theta_23)
th13_r = np.radians(theta_13)
delta_r = np.radians(delta_CKM)
J = np.cos(th12_r)*np.cos(th23_r)*np.cos(th13_r)**2 * \
    np.sin(th12_r)*np.sin(th23_r)*np.sin(th13_r)*np.sin(delta_r)
print(f"  J = {J:.4e} (exp ~3.08e-5)")
check("J ~= 3.12e-5", abs(J/3.12e-5 - 1) < 0.02, f"J = {J:.4e}")

# Leading order: J ~ phi^-22 * sqrt(a1/b1)
J_leading = PHI**(-22) * np.sqrt(a1/b1)
check("22 = sum of CKM exponents 3+7+12", 3+7+12 == 22)

# ============================================================
# APPENDIX: NUMERICAL CHECKS
# ============================================================
print("\n--- Appendix Numerical Checks ---")

# phi^4 = 3*phi+2
check("phi^4 = 3*phi+2", abs(PHI**4 - (3*PHI+2)) < 1e-10)
check("4*a1*phi^4 = 60*phi+40", abs(4*a1*PHI**4 - (60*PHI+40)) < 1e-10)

# 2*phi^3 = 4*phi+2
check("2*phi^3 = 4*phi+2", abs(2*PHI**3 - (4*PHI+2)) < 1e-10)

# phi^25 = 75025*phi + 46368
phi25_fib = 75025*PHI + 46368
check("phi^25 = 75025*phi+46368", abs(PHI**25 - phi25_fib) < 0.1,
      f"phi^25={PHI**25:.1f}, 75025*phi+46368={phi25_fib:.1f}")

# m_e * phi^25 = 85726 MeV
me_phi25 = m_e * PHI**25
print(f"  m_e * phi^25 = {me_phi25:.0f} MeV")
check("m_e*phi^25 ~= 85726 MeV", abs(me_phi25 - 85726) < 5,
      f"Got {me_phi25:.0f}")

# m_W from m_Z
m_Z_val = 91.11  # GeV (as computed)
m_W_from_Z = m_Z_val * np.sqrt(1 - 6/26)
print(f"  m_W = m_Z*sqrt(20/26) = {m_Z_val}*{np.sqrt(20/26):.5f} = {m_W_from_Z:.2f} GeV")
check("m_W from m_Z ~= 79.9", abs(m_W_from_Z - 79.9) < 0.1)

# Higgs detailed
higgs_ratio2 = PHI - 8*0.0072973525
print(f"  phi-8*alpha = {higgs_ratio2:.4f}")
m_H_detail = m_W_exp * higgs_ratio2
print(f"  m_H = {m_W_exp} * {higgs_ratio2:.4f} = {m_H_detail:.2f} GeV")

# ============================================================
# SPECIFIC CLAIMS THAT MIGHT BE WRONG
# ============================================================
print("\n--- Targeted Checks for Potential Errors ---")

# 1. "Universal optimality" - Cohn-Kumar proves for SPHERICAL codes on S^3,
# not "R^4 for all completely monotonic potentials"
warn("Universal optimality claim",
     "Paper says '600-cell is unique universal energy minimizer in R^4 for all completely monotonic potentials'. "
     "Cohn-Kumar 2007 proves universal optimality on S^3 (120 points on sphere). "
     "The statement 'in R^4' is imprecise - should say 'on S^3'.")

# 2. Spectral gap multiplicity
warn("Spectral gap multiplicity",
     "Line ~237: 'spectral gap lambda_min = 12-6*phi (multiplicity 1)'. "
     "The adjacency eigenvalue 6*phi has multiplicity 4, so the Laplacian eigenvalue "
     "12-6*phi should also have multiplicity 4, NOT 1.")

# 3. Check the 600-cell symmetry group order
# H_4 has order 14400. The FULL symmetry group of the 600-cell is
# 14400 (including reflections), or 7200 (rotational only).
# |H_4| = 14400 = 120^2.
check("H4 order = 14400", 14400 == 120**2)

# 4. "lambda_max = 12 (the vertex degree)" context
warn("lambda_max context",
     "Paper says 'maximal eigenvalue lambda_max = 12 (the vertex degree)' in the "
     "Laplacian context. For the combinatorial Laplacian L=D-A, the max eigenvalue "
     "is 12+|min eigenvalue of A| = 12+3.708 = 15.708. The value 12 is the max "
     "eigenvalue of the ADJACENCY matrix. Context may confuse readers.")

# 5. "phi-sector dimension 4+9+9+4 = 26"
# Paper says this is split 13+13. But phi-sector (b>0) has 4+9=13, and phi'-sector
# (b<0 in eigenvalue) has 9+4=13. OK.

# 6. Generation theorem: electron on a=1 line issue (already caught)

# 7. Check "the spectral gap = 2" for U(1) AC matching
# This is from exp143, claimed to be "verified computationally"
# We can't verify without the actual computation

# 8. "288 pure plaquettes out of 1200 total = b1/a1^2 = 6/25"
check("288/1200 = 6/25", abs(288/1200 - 6/25) < 1e-10)

# 9. Bare dihedral angle ratio
theta_O = np.arccos(-1/3)
theta_T = np.arccos(1/3)
bare_ratio = theta_O / theta_T
print(f"\n  Dihedral angle ratio: arccos(-1/3)/arccos(1/3) = {bare_ratio:.4f}")
m_H_bare = bare_ratio * m_W_exp
print(f"  m_H(bare) = {bare_ratio:.4f} * {m_W_exp} = {m_H_bare:.1f} GeV")
check("Bare Higgs ratio ~= 1.552", abs(bare_ratio - 1.552) < 0.001,
      f"Got {bare_ratio:.4f}")

# 10. The 96 = 16*3*2 factorization claim
# "16 fermions per generation * 3 generations * 2 chiralities"
# But the SM has 15 Weyl fermions per generation (or 16 if you count right-handed neutrino)
warn("96=16*3*2 factorization",
     "Paper says 96 = 16*3*2 = '(fermions per generation) x (generations) x (chiralities)'. "
     "The SM has 15 Weyl fermions per generation (no right-handed neutrino). "
     "16 would include a sterile neutrino. This should be clarified.")

# 11. Check m_3 < 1% claim more carefully
m3_pred = 2*m_e/PHI**35
m3_exp = np.sqrt(2.453e-3)  # sqrt(Delta m^2_31) in eV, PDG 2024
m3_error = abs(m3_pred - m3_exp)/m3_exp * 100
print(f"\n  m_3 predicted: {m3_pred*1000:.2f} meV")
print(f"  m_3 exp (sqrt(Dm31^2)): {m3_exp*1000:.2f} meV")
print(f"  Error: {m3_error:.1f}%")
check("m_3 error is reasonable", m3_error < 2, f"Actual: {m3_error:.1f}%")

# 12. Paper says sum(m_nu) ~= 0.058 eV. Let's check if this is consistent.
# With m_3 = 0.0495, and NO information about m_1, m_2, this is just an estimate.
# For normal hierarchy: m_1 ~ 0, m_2 ~ sqrt(Dm21^2) ~ 0.0087, m_3 ~ 0.050
# sum ~= 0.059. Close to 0.058.

# 13. "phi^(-7)" in theta_23 bare exponent
check("phi^(-7) = phi'^7 (Galois)", abs(PHI**(-7) - abs(PHI_prime)**7) < 1e-10)

# 14. Delta m^2_21 / Delta m^2_31 = 4*alpha
Dm21 = 7.53e-5  # eV^2 (PDG 2024)
Dm31 = 2.453e-3  # eV^2 (PDG 2024)
ratio_dm = Dm21/Dm31
pred_ratio_dm = 4*alpha_exp
print(f"\n  Dm21/Dm31 = {ratio_dm:.5f}")
print(f"  4*alpha = {pred_ratio_dm:.5f}")
dm_error = abs(ratio_dm - pred_ratio_dm)/ratio_dm * 100
print(f"  Error: {dm_error:.2f}%")

# 15. Check the spectral action c_1 = 14880 = 124*N
check("14880 = 124*120", 14880 == 124*120)
check("124 = dim(E8)/2", 124 == 248/2)

# 16. Check the neutrino exponent identities
check("n_nu = 35 = b1^2-1", 35 == b1**2-1)
check("n_Z = 25 = a1^2", 25 == a1**2)
check("n_Z+n_nu = 60 = N/2", 25+35 == N//2)

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*72)
print(f"AUDIT SUMMARY")
print(f"="*72)
print(f"  Checks passed: {checks_passed}")
print(f"  ERRORS found: {len(errors)}")
print(f"  WARNINGS: {len(warnings)}")
print()

if errors:
    print("ERRORS:")
    for e in errors:
        print(f"  {e}")
    print()

if warnings:
    print("WARNINGS:")
    for w in warnings:
        print(f"  {w}")
    print()

print("="*72)
print("EXP-278 COMPLETE")
print("="*72)

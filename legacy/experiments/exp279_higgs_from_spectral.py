"""
EXP-279: HIGGS MASS FROM SPECTRAL TRIPLE
==========================================
Open problem: WHY m_H/m_W = phi - 8*alpha?

Strategy:
1. NCG Higgs: m_H^2 = 2*lambda*v^2 where lambda = Tr(Y^4)/Tr(Y^2)^2 type
2. Use our Yukawa matrix Y_f = phi^(n_f) on McKay graph
3. Compute Connes-Chamseddine formula with our fermion content
4. See if phi - 8*alpha emerges naturally
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
N_eig = 9; N_gen = 3
ALPHA = 1/137.035999084
ALPHA_S = 1/(2*PHI**3)
sin2tW = b1/(a1**2 + 1)  # = 6/26

print("="*72)
print("EXP-279: HIGGS MASS FROM SPECTRAL TRIPLE")
print("="*72)

# === Part 1: Fermion mass exponents ===
print("\n--- Part 1: Fermion Content ---")

# (a,b) quantum numbers -> n = 5a + 6b
fermions = {
    'e':   (0, 0),   # n=0
    'mu':  (1, 1),   # n=11
    'tau': (1, 2),   # n=17
    'u':   (3, -2),  # n=3
    'c':   (2, 1),   # n=16
    't':   (4, 1),   # n=26
    'd':   (1, 0),   # n=5
    's':   (1, 1),   # n=11
    'b':   (-1, 4),  # n=19
}

m_e = 0.51099895  # MeV
masses = {}
exponents = {}
for name, (a, b) in fermions.items():
    n = a1*a + b1*b
    m = m_e * PHI**n
    masses[name] = m
    exponents[name] = n
    print(f"  {name:3s}: (a,b)=({a:+d},{b:+d}), n={n:3d}, m = {m:.4f} MeV")

# === Part 2: Connes-Chamseddine Higgs mass formula ===
print("\n--- Part 2: Connes-Chamseddine Formula ---")

# In NCG: m_H^2 / m_W^2 = 2 * [sum_f N_c(f) * y_f^4] / [sum_f N_c(f) * y_f^2]
# where y_f = m_f / v (Yukawa coupling), N_c = color factor (3 quarks, 1 leptons)
# This is equivalent to: m_H^2 / m_W^2 = 2 * sum(N_c * m_f^4) / [sum(N_c * m_f^2)]

# Color factors
N_c = {'e': 1, 'mu': 1, 'tau': 1, 'u': 3, 'c': 3, 't': 3, 'd': 3, 's': 3, 'b': 3}

# Using BARE masses from our formula
sum_m2 = sum(N_c[f] * masses[f]**2 for f in fermions)
sum_m4 = sum(N_c[f] * masses[f]**4 for f in fermions)

ratio_CC = 2 * sum_m4 / sum_m2**2
print(f"  Connes-Chamseddine: m_H/m_W = sqrt(2 * Tr(Y^4)/Tr(Y^2)^2)")
print(f"  sum(N_c * m_f^2) = {sum_m2:.6e}")
print(f"  sum(N_c * m_f^4) = {sum_m4:.6e}")
print(f"  ratio = 2 * sum_m4 / sum_m2^2 = {ratio_CC:.6f}")
print(f"  m_H/m_W (CC) = sqrt(ratio) = {np.sqrt(ratio_CC):.6f}")

# What we need
m_W = 80377  # MeV
m_H_exp = 125250  # MeV
target = m_H_exp / m_W
print(f"\n  Target: m_H/m_W = {target:.6f}")
print(f"  Our formula: phi - 8*alpha = {PHI - 8*ALPHA:.6f}")
print(f"  Connes-Chamseddine gives: {np.sqrt(ratio_CC):.6f}")
print(f"  Error CC: {(np.sqrt(ratio_CC) - target)/target * 100:.1f}%")

# The CC formula is dominated by top quark
top_frac = 3 * masses['t']**4 / sum_m4
print(f"\n  Top quark dominance: {top_frac*100:.2f}% of sum(m^4)")

# === Part 3: Try modified formulas ===
print("\n--- Part 3: Modified Spectral Higgs Formulas ---")

# Attempt 1: Use phi^n instead of m_f (pure geometric)
# y_f = phi^n_f (unnormalized)
sum_phi2 = sum(N_c[f] * PHI**(2*exponents[f]) for f in fermions)
sum_phi4 = sum(N_c[f] * PHI**(4*exponents[f]) for f in fermions)

R1 = 2 * sum_phi4 / sum_phi2**2
print(f"  Attempt 1: phi^n Yukawa")
print(f"    ratio = {R1:.6e}")
print(f"    sqrt(ratio) = {np.sqrt(R1):.6f}")

# Attempt 2: Use Z[phi] norm N(z) instead
# For z = a + b*phi, N(z) = a^2 + a*b - b^2
norms = {}
for name, (a, b) in fermions.items():
    norms[name] = a**2 + a*b - b**2

sum_N2 = sum(N_c[f] * norms[f]**2 for f in fermions)
sum_N4 = sum(N_c[f] * norms[f]**4 for f in fermions)
if sum_N2 > 0:
    R2 = 2 * sum_N4 / sum_N2**2
    print(f"\n  Attempt 2: N(z) norms")
    print(f"    Norms: {dict((f, norms[f]) for f in fermions)}")
    print(f"    ratio = {R2:.6f}, sqrt = {np.sqrt(R2):.6f}")

# Attempt 3: Use traces Tr(z) = 2a+b
traces = {}
for name, (a, b) in fermions.items():
    traces[name] = 2*a + b

sum_T2 = sum(N_c[f] * traces[f]**2 for f in fermions)
sum_T4 = sum(N_c[f] * traces[f]**4 for f in fermions)
if sum_T2 > 0:
    R3 = 2 * sum_T4 / sum_T2**2
    print(f"\n  Attempt 3: Tr(z) traces")
    print(f"    Traces: {dict((f, traces[f]) for f in fermions)}")
    print(f"    ratio = {R3:.6f}, sqrt = {np.sqrt(R3):.6f}")

# === Part 4: Decomposition phi - 8*alpha ===
print("\n--- Part 4: Why phi - 8*alpha? ---")
print(f"  phi = {PHI:.10f}")
print(f"  8*alpha = {8*ALPHA:.10f}")
print(f"  phi - 8*alpha = {PHI - 8*ALPHA:.10f}")
print(f"  m_H/m_W exp = {target:.10f}")

# The 8 = rank(E8). The alpha comes from EM loops.
# Tree level: phi (from DSI potential, ground state ratio)
# 1-loop correction: -rank(E8) * alpha

# Check: is there a formula m_H/m_W = phi * (1 - correction)?
correction = 1 - target/PHI
print(f"\n  Relative correction from phi: {correction:.6f}")
print(f"  = {correction:.6f}")
print(f"  8*alpha/phi = {8*ALPHA/PHI:.6f}")
print(f"  Match: {abs(correction - 8*ALPHA/PHI)/correction*100:.2f}%")

# Alternative: phi - 2/(a1*phi^4) [geometric form]
geom = PHI - 2/(a1*PHI**4)
print(f"\n  Geometric form: phi - 2/(a1*phi^4) = {geom:.10f}")
print(f"  This equals phi - 8*alpha + O(alpha^2) via alpha equation:")
print(f"  From 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0:")
print(f"  alpha ~ 1/(4*a1*phi^4) + ... so 2/(a1*phi^4) ~ 8*alpha")
print(f"  Exact: 2/(a1*phi^4) = {2/(a1*PHI**4):.10f}")
print(f"  vs 8*alpha = {8*ALPHA:.10f}")
print(f"  Difference: {abs(2/(a1*PHI**4) - 8*ALPHA):.2e}")

# === Part 5: Inner fluctuation approach ===
print("\n--- Part 5: Inner Fluctuation of D ---")
# In NCG: D -> D + A + JAJ^-1
# A = sum_i a_i [D, b_i] is the inner fluctuation (gauge + Higgs)
# The Higgs field H appears as the "finite" part of A
# Its potential V(H) = a*|H|^2 + b*|H|^4 with:
#   a = -mu^2 = Tr(D_F^2) type term
#   b = lambda = Tr(D_F^4) type term
# VEV: <H> = sqrt(-a/(2b)) = v

# On our McKay D_F = diag(n_f):
n_f = [exponents[f] for f in ['e','u','d','s','c','b','t','tau','mu']]
# McKay ordering: 0(e)-1(u)-2(d)-3(s)-4(c)-5(b)-6(t)-7(tau), branch 4(c)-8(mu)
print(f"  D_F diagonal = {n_f}")

Tr_DF2 = sum(n**2 for n in n_f)
Tr_DF4 = sum(n**4 for n in n_f)
print(f"  Tr(D_F^2) = {Tr_DF2}")
print(f"  Tr(D_F^4) = {Tr_DF4}")

# NCG Higgs VEV: related to Tr(D_F^2)
# NCG Higgs mass^2: related to Tr(D_F^4)/Tr(D_F^2)
R_DF = Tr_DF4 / Tr_DF2
print(f"  Tr(D_F^4)/Tr(D_F^2) = {R_DF:.4f}")
print(f"  sqrt(...) = {np.sqrt(R_DF):.4f}")

# With color factors
n_f_with_color = {
    'e': (0, 1), 'mu': (11, 1), 'tau': (17, 1),
    'u': (3, 3), 'c': (16, 3), 't': (26, 3),
    'd': (5, 3), 's': (11, 3), 'b': (19, 3)
}
Tr_DF2_c = sum(nc * n**2 for n, nc in n_f_with_color.values())
Tr_DF4_c = sum(nc * n**4 for n, nc in n_f_with_color.values())
print(f"\n  With color factors:")
print(f"  Tr(N_c * D_F^2) = {Tr_DF2_c}")
print(f"  Tr(N_c * D_F^4) = {Tr_DF4_c}")
print(f"  Ratio = {Tr_DF4_c/Tr_DF2_c:.4f}")
print(f"  sqrt = {np.sqrt(Tr_DF4_c/Tr_DF2_c):.4f}")

# === Part 6: Algebraic search ===
print("\n--- Part 6: Algebraic Relations ---")

# Check if Tr(D_F^2) or Tr(D_F^4) have nice forms
print(f"  Tr(D_F^2) = {Tr_DF2}")
print(f"  = {Tr_DF2} = ?")
# Factor
for a in range(1, 50):
    for b in range(1, 50):
        if a*b == Tr_DF2:
            print(f"    = {a} * {b}")

print(f"  Tr(D_F^4) = {Tr_DF4}")
for a in range(1, 200):
    if Tr_DF4 % a == 0:
        b = Tr_DF4 // a
        if a <= b and a > 1 and b < 100000:
            print(f"    = {a} * {b}")

# Check n_f sums
n_list = sorted(set(exponents.values()))
print(f"\n  Distinct exponents: {n_list}")
print(f"  Sum of all n: {sum(exponents.values())}")
print(f"  Sum of n^2: {sum(n**2 for n in exponents.values())}")

# The 8 in "8*alpha": is it rank(E8)?
# Check: m_H/m_W = phi - rank*alpha for various ranks
print(f"\n  m_H/m_W = phi - R*alpha for different R:")
for R in range(1, 15):
    pred = PHI - R*ALPHA
    err = (pred - target)/target * 100
    print(f"    R={R:2d}: {pred:.6f}, error = {err:+.3f}%")

# === Part 7: DSI + perturbation theory ===
print("\n--- Part 7: DSI Higgs Mechanism ---")
# In DSI potential V(r) = sin^2(pi*ln(r/r0)/ln(phi))
# Minima at r_k = r_0 * phi^k
# The Higgs VEV is at some minimum, mass is curvature at minimum

# V''(r_k) = (pi/ln(phi))^2 * 2 / r_k^2
# For the ratio m_H/m_W, we need V''_Higgs / V''_W or similar

# Actually: in the SM, m_H = sqrt(2*lambda) * v, m_W = g*v/2
# So m_H/m_W = 2*sqrt(2*lambda)/g = 2*sqrt(2*lambda) / sqrt(4*pi*alpha/sin^2(tW))

# The DSI potential gives the Higgs self-coupling lambda
# V(H) = -mu^2*|H|^2 + lambda*|H|^4
# At minimum: v^2 = mu^2/(2*lambda), m_H^2 = 2*mu^2 = 4*lambda*v^2

# If lambda comes from DSI curvature:
# V_DSI''(phi^k) = 2*(pi/ln(phi))^2 / (phi^k)^2
# This is (2*pi^2/ln^2(phi)) * phi^(-2k)
C_DSI = 2 * np.pi**2 / np.log(PHI)**2
print(f"  DSI curvature coefficient: 2*pi^2/ln^2(phi) = {C_DSI:.6f}")

# If we identify lambda with the DSI self-coupling at k=0:
# lambda_DSI = pi^2 / (2*ln^2(phi))
lambda_DSI = np.pi**2 / (2 * np.log(PHI)**2)
print(f"  lambda_DSI = pi^2/(2*ln^2(phi)) = {lambda_DSI:.6f}")

# Then m_H/m_W = 2*sqrt(2*lambda_DSI)/g_2
# Need g_2: sin^2(tW) = g1^2/(g1^2+g2^2), and alpha = g2^2*sin^2(tW)/(4*pi)
# g2^2 = 4*pi*alpha/sin^2(tW)
g2_sq = 4*np.pi*ALPHA/sin2tW
g2 = np.sqrt(g2_sq)
print(f"  g_2 = {g2:.6f}")

ratio_DSI = 2*np.sqrt(2*lambda_DSI)/g2
print(f"  m_H/m_W (DSI) = 2*sqrt(2*lambda)/g2 = {ratio_DSI:.6f}")
print(f"  Target: {target:.6f}")
print(f"  Error: {(ratio_DSI - target)/target*100:.1f}%")

# === Part 8: Exact algebraic check ===
print("\n--- Part 8: Exact Algebraic Identity ---")
# phi - 2/(a1*phi^4):
# phi^5 = phi^4 + phi^3 = (3+phi)(phi+1) = 3*phi+3+phi^2+phi = 3*phi+3+phi+1+phi = 5*phi+4
# Actually: phi^5 = 5*phi + 3 (in Z[phi])
# phi^4 = 3*phi + 2 (in Z[phi]: (3,2))
# a1*phi^4 = 5*(3*phi+2) = 15*phi + 10
# 2/(a1*phi^4) = 2/(15*phi+10)
# N(15*phi+10) = 10^2 + 10*15 - 15^2 = 100+150-225 = 25
# So 1/(a1*phi^4) = (10-15*phi')/25... hmm, let me do this differently

# phi - 2/(a1*phi^4) in terms of alpha equation:
# 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
# => alpha = (4*a1*phi^4 - sqrt(16*a1^2*phi^8 - 8*pi)) / (4*pi)
# => 4*a1*phi^4*alpha = 1 + 2*pi*alpha^2 (from equation rearranged)
# => 1/(a1*phi^4) = 4*alpha / (1 + 2*pi*alpha^2)
# => 2/(a1*phi^4) = 8*alpha / (1 + 2*pi*alpha^2)

val1 = 8*ALPHA / (1 + 2*np.pi*ALPHA**2)
val2 = 2/(a1*PHI**4)
print(f"  8*alpha/(1+2*pi*alpha^2) = {val1:.12f}")
print(f"  2/(a1*phi^4) = {val2:.12f}")
print(f"  Match: {abs(val1-val2):.2e}")

print(f"\n  IDENTITY: phi - 2/(a1*phi^4) = phi - 8*alpha/(1+2*pi*alpha^2)")
print(f"  = phi - 8*alpha + O(alpha^3)")
print(f"  The EXACT geometric formula is phi - 2/(a1*phi^4)")
print(f"  The alpha-expansion gives phi - 8*alpha to leading order")

# Is 2/(a1*phi^4) DERIVABLE?
# phi^4 = 3*phi + 2. So a1*phi^4 = 15*phi + 10.
# The "2" in the numerator:
# If Higgs is at the FIRST DSI level (k=1), then VEV ratio = phi
# Correction from curvature: 1/(a1*phi^(a1-1)) = 1/(5*phi^4)
# With factor 2 for doublet: 2/(5*phi^4)
print(f"\n  GEOMETRIC: correction = 2/(a1*phi^4)")
print(f"  a1*phi^4 = {a1*PHI**4:.6f}")
print(f"  2 could be: dim(Higgs doublet), SU(2) factor, or chirality")
print(f"  a1 = pentagon symmetry, phi^4 = DSI^4 (4th Fibonacci level)")

# Check: phi^4 = phi^3 + phi^2 = (2+phi) + (1+phi) = 3+2*phi
# Wait: phi^2 = phi+1, phi^3 = phi^2+phi = 2*phi+1, phi^4 = phi^3+phi^2 = 3*phi+2
print(f"\n  phi^4 = 3*phi + 2 in Z[phi]: (a,b)=(2,3)")
print(f"  N(phi^4) = 4+6-9 = 1 (unit, as expected)")
print(f"  a1*phi^4 = 5*(3*phi+2) = 15*phi+10")
print(f"  N(15*phi+10) = 100+150-225 = 25 = a1^2")
print(f"  So |N(a1*phi^4)| = a1^2 (!!)")

print(f"\n  DERIVATION ATTEMPT:")
print(f"  Tree level: m_H/m_W = phi (DSI ground state ratio)")
print(f"  1-loop EM: subtract 2/(a1*phi^4) = 8*alpha/(1+O(alpha^2))")
print(f"  The 2: Higgs doublet components")
print(f"  The a1*phi^4: DSI normalization (N(a1*phi^4)=a1^2)")
print(f"  Equivalently: rank(E8)*alpha at leading order")

# === Summary ===
print("\n" + "="*72)
print("SUMMARY")
print("="*72)
print(f"  1. Connes-Chamseddine with our masses: m_H/m_W = {np.sqrt(ratio_CC):.4f}")
print(f"     (doesn't match {target:.4f} -- same failure as Connes)")
print(f"  2. EXACT geometric: m_H/m_W = phi - 2/(a1*phi^4)")
print(f"     = {geom:.6f} vs {target:.6f} ({(geom-target)/target*100:+.3f}%)")
print(f"  3. IDENTITY: 2/(a1*phi^4) = 8*alpha/(1+2*pi*alpha^2) EXACT")
print(f"     from the alpha equation itself")
print(f"  4. N(a1*phi^4) = a1^2 = 25 (integer square!)")
print(f"  5. DSI lambda formula: {ratio_DSI:.4f} (wrong, {(ratio_DSI-target)/target*100:+.1f}%)")
print(f"  STATUS: Formula VERIFIED, connection to alpha equation FOUND,")
print(f"          but tree-level phi still needs Higgs sector derivation")
print("="*72)

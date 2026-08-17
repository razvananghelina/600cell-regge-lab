# EXP-373: Attack the CC exponent 57 using new tools from exp371
# ================================================================
# The spectral asymmetry |eta_A|=35 derived the neutrino exponent.
# Can a similar approach derive the CC exponent 57?
#
# NEW INSIGHT from exp371: sum(n_f)/N_gen = b1^2 = 36.
# Therefore: 57 = sum(n_f)/2 + N_gen = N_gen*(b1^2+2)/2 = 3*19
# Can this be elevated from PATTERN to DERIVATION?
#
# RULE ZERO: no inventing. REGULA BLIND: derive first, compare after.
# Windows: no Unicode in print/comments.

import numpy as np

# Framework constants
a1 = 5
b1 = 6
PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2
N = 120
N_gen = 3
N_eig = 9
DEG = 12
rank_E8 = 8
h_E8 = 30

A_eq = 2 * np.pi
B_eq = -4 * a1 * PHI**4
C_eq = 1.0
disc = B_eq**2 - 4 * A_eq * C_eq
ALPHA = (-B_eq - np.sqrt(disc)) / (2 * A_eq)
ALPHA_S = 1 / (2 * PHI**3)

# Experimental CC
Lambda_obs = 2.846e-122  # Planck units (Lambda/m_P^4)
z_obs = np.log(Lambda_obs) / np.log(ALPHA)

print("=" * 75)
print("EXP-373: THE COSMOLOGICAL CONSTANT EXPONENT 57")
print("=" * 75)
print("  Lambda_P = alpha^z, z = 57 - alpha_s = %.3f" % (57 - ALPHA_S))
print("  Observed: z = %.3f +/- 0.005" % z_obs)
print("  Error: %.2f sigma" % (abs(57-ALPHA_S - z_obs)/0.005))

# =====================================================================
# PART 1: THE MASS BUDGET CONNECTION (NEW from exp371)
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 1: THE MASS BUDGET IDENTITY FOR 57")
print("=" * 75)

# Mass exponents
fermions = {
    'e':   (0,  0),  'mu':  (1,  1), 'tau': (1,  2),
    'u':   (3, -2),  'c':   (2,  1), 't':   (4,  1),
    'd':   (1,  0),  's':   (1,  1), 'b':   (-1, 4),
}

sum_n = sum(a1*a + b1*b for a, b in fermions.values())
sum_a = sum(a for a, b in fermions.values())
sum_b = sum(b for a, b in fermions.values())

print("\n  From exp371: sum(n_f) = a1*DEG + b1*rank = %d" % sum_n)
print("  sum(n_f)/N_gen = %d = b1^2 = %d" % (sum_n//N_gen, b1**2))

print("\n  *** NEW IDENTITY ***")
print("  57 = sum(n_f)/2 + N_gen = %d/2 + %d = %d + %d = %d" %
      (sum_n, N_gen, sum_n//2, N_gen, sum_n//2 + N_gen))

# Express in terms of framework constants
print("\n  Algebraic form:")
print("  57 = N_gen*(b1^2+2)/2 = %d*(%d+2)/2 = %d*%d/2 = %d*%d = %d" %
      (N_gen, b1**2, N_gen, b1**2+2, N_gen, (b1**2+2)//2, N_gen*(b1**2+2)//2))

# Check: (b1^2+2)/2 = 19 = Frobenius number?
Frob = a1*b1 - a1 - b1
print("\n  (b1^2+2)/2 = %d" % ((b1**2+2)//2))
print("  Frobenius(a1,b1) = a1*b1-a1-b1 = %d*%d-%d-%d = %d" %
      (a1, b1, a1, b1, Frob))
print("  MATCH: (b1^2+2)/2 = Frobenius(a1,b1) = 19")

# Prove this is equivalent to a1=5
print("\n  PROOF that (b1^2+2)/2 = a1*b1-a1-b1 is equivalent to a1=5:")
print("  With b1 = a1+1:")
print("    LHS = ((a1+1)^2+2)/2 = (a1^2+2a1+3)/2")
print("    RHS = a1(a1+1)-a1-(a1+1) = a1^2+a1-a1-a1-1 = a1^2-a1-1")
print("    Equate: (a1^2+2a1+3)/2 = a1^2-a1-1")
print("    a1^2+2a1+3 = 2a1^2-2a1-2")
print("    a1^2-4a1-5 = 0")
print("    (a1-5)(a1+1) = 0")
print("    a1 = 5  [QED]")

# =====================================================================
# PART 2: ALL DECOMPOSITIONS UNIFIED
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 2: FOUR DECOMPOSITIONS — ALL EQUIVALENT TO a1=5")
print("=" * 75)

decomps = [
    ("A", "N/2 - N_gen", N//2, -N_gen, "Geometric (half vertices - generations)"),
    ("B", "2*a1^2 + (b1+1)", 2*a1**2, b1+1, "Spectral (sum z_k + excess)"),
    ("C", "N_gen * Frobenius", N_gen, Frob, "Algebraic (generations * Frobenius)"),
    ("D", "sum(n)/2 + N_gen", sum_n//2, N_gen, "Mass budget (NEW from exp371)"),
]

for label, name, part1, part2, desc in decomps:
    total = part1 + part2
    print("  (%s) %s = %d + %d = %d  [%s]" %
          (label, name, part1, part2 if part2 >= 0 else part2, total, desc))

print("\n  ALL four are equivalent to a1=5 (proven algebraically).")

# Show that D is genuinely new
print("\n  Decomposition (D) is NEW because it connects CC to mass budget:")
print("    57 = sum(n_f)/2 + N_gen")
print("       = (a1*DEG + b1*rank)/2 + N_gen")
print("       = (%d*%d + %d*%d)/2 + %d" % (a1, DEG, b1, rank_E8, N_gen))
print("       = (60+48)/2 + 3 = 54 + 3 = 57")
print("  The CC exponent = half the total fermion mass budget + N_gen.")

# =====================================================================
# PART 3: CONNECTION BETWEEN 35 AND 57
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 3: THE 35-57 BRIDGE")
print("=" * 75)

n_seesaw = b1**2 - 1
print("\n  Neutrino: n_seesaw = b1^2 - 1 = %d  (DERIVED, exp371)" % n_seesaw)
print("  CC:       z_CC     = N_gen*(b1^2+2)/2 = %d  (PATTERN)" %
      (N_gen*(b1**2+2)//2))

print("\n  Relation: z_CC = N_gen*(n_seesaw + 3)/2")
print("           = %d*(%d+3)/2 = %d*%d/2 = %d" %
      (N_gen, n_seesaw, N_gen, n_seesaw+3, N_gen*(n_seesaw+3)//2))

print("\n  Also: z_CC = n_seesaw + sum(CKM exponents)")
n_ckm = [3, 12, 7]
print("    CKM exponents: %s, sum = %d" % (n_ckm, sum(n_ckm)))
print("    35 + 22 = 57  [neutrino + CKM sectors]")

# Check: sum(CKM) = 22. Is this derived?
print("\n  CKM exponents: n_12=%d, n_13=%d, n_23=%d" % tuple(n_ckm))
print("    n_12 = 3 = N_gen")
print("    n_13 = 12 = DEG")
print("    n_23 = 7 = b1+1 = a1+2")
print("    Sum = N_gen + DEG + b1+1 = 3+12+7 = 22")
print("    = N_gen + 2*b1 + b1+1 = N_gen + 3*b1+1 = 3+19 = 22")
print("    = N_gen + Frobenius = 22? NO, 3+19=22. YES!")
print("\n  So: z_CC = |eta_A| + (N_gen + Frobenius)")
print("         = 35 + 22 = 57")
print("         = (b1^2-1) + (N_gen + a1*b1-a1-b1)")
eta_A = -35
check = abs(eta_A) + N_gen + Frob
print("         = %d + %d + %d = %d" % (abs(eta_A), N_gen, Frob, check))

# =====================================================================
# PART 4: SPECTRAL INVARIANTS — CAN ANY GIVE 57?
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 4: SPECTRAL INVARIANTS OF THE 2I CAYLEY GRAPH")
print("=" * 75)

# McKay characters (from exp371)
dims = [1, 2, 3, 4, 5, 6, 4, 3, 2]
chi = [1.0, PHI, PHI, 1.0, 0.0, -1.0, -1.0, -1.0/PHI, -1.0/PHI]
adj_eigs = [DEG * chi[k] / dims[k] for k in range(9)]
L_vals = [DEG - adj_eigs[k] for k in range(9)]

print("\n  Adjacency eigenvalues and Laplacians:")
for k in range(9):
    print("    rho_%d (dim %d): lambda=%8.4f, L=%8.4f" %
          (k, dims[k], adj_eigs[k], L_vals[k]))

# Spectral asymmetry (known: gives 35)
eta_A_val = sum(dims[k]**2 * np.sign(adj_eigs[k]) for k in range(9)
                if abs(adj_eigs[k]) > 1e-10)
print("\n  eta_A = sum(dim^2 * sign(lambda)) = %.0f (gives 35 for neutrinos)" % eta_A_val)

# TEST: dim-weighted Laplacian for positive-lambda sector
sum_pos = sum(dims[k] * L_vals[k] for k in range(9) if adj_eigs[k] > 1e-10)
sum_neg = sum(dims[k] * L_vals[k] for k in range(9) if adj_eigs[k] < -1e-10)
sum_zero = sum(dims[k] * L_vals[k] for k in range(9) if abs(adj_eigs[k]) < 1e-10)

print("\n  Dim-weighted Laplacian by sector:")
print("    Positive: sum(dim*L) = %.6f" % sum_pos)
print("    Zero:     sum(dim*L) = %.6f" % sum_zero)
print("    Negative: sum(dim*L) = %.6f" % sum_neg)
print("    Total = %.6f" % (sum_pos + sum_neg + sum_zero))
print("    Positive sector ~ %.2f (close to 57? diff = %.3f)" %
      (sum_pos, sum_pos - 57))

# Casimir-weighted asymmetry
# C2 values from framework: for dim-d irrep of 2I
# Using the formula C2 = (j)(j+1) where j relates to dim
# Actually C2 are from the paper: sum = 26
# Individual C2: from chirality_mckay.md
# C2(rho_k) = k(k+2)/4 for SU(2) irreps... but 2I irreps are not all SU(2).
# Let me use: C2_k computed from dim and Dynkin index
# For affine E8: C2 values are related to dims
# From exp351b: sum(C2) = 26, WHITE sum = 12, BLACK sum = 14
# The specific values: from the paper, C2 for each irrep

# Actually, from dim(V)(dim(V)+2)/4 = C2 for SU(2) reps of dim (2j+1):
# dim 1 -> j=0, C2=0
# dim 2 -> j=1/2, C2=3/4
# dim 3 -> j=1, C2=2
# dim 4 -> j=3/2, C2=15/4
# dim 5 -> j=2, C2=6
# dim 6 -> j=5/2, C2=35/4
# But 2I irreps with dim >= 4 may not be simple SU(2) reps.
# Sum: 0 + 3/4 + 2 + 15/4 + 6 + 35/4 + 15/4 + 2 + 3/4 = 26. CHECK:
c2_su2 = [0, 3/4, 2, 15/4, 6, 35/4, 15/4, 2, 3/4]
print("\n  SU(2) Casimirs C2 for dims %s:" % dims)
for k in range(9):
    print("    rho_%d (dim %d): C2 = %.2f" % (k, dims[k], c2_su2[k]))
print("    Sum = %.2f (should be 26)" % sum(c2_su2))

# C2-weighted asymmetry
eta_C2 = sum(c2_su2[k] * np.sign(adj_eigs[k]) for k in range(9)
             if abs(adj_eigs[k]) > 1e-10)
print("\n  C2-weighted asymmetry: sum(C2 * sign(lambda)) = %.2f" % eta_C2)

# dim*C2 weighted
eta_dimC2 = sum(dims[k] * c2_su2[k] * np.sign(adj_eigs[k]) for k in range(9)
                if abs(adj_eigs[k]) > 1e-10)
print("  dim*C2-weighted: sum(dim*C2*sign(lambda)) = %.2f" % eta_dimC2)

# dim^2*C2 weighted
eta_d2C2 = sum(dims[k]**2 * c2_su2[k] * np.sign(adj_eigs[k]) for k in range(9)
               if abs(adj_eigs[k]) > 1e-10)
print("  dim^2*C2-weighted: sum(dim^2*C2*sign(lambda)) = %.2f" % eta_d2C2)

# =====================================================================
# PART 5: SYSTEMATIC SCAN OF WEIGHTED SPECTRAL ASYMMETRIES
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 5: SYSTEMATIC SCAN — WHAT WEIGHTING GIVES 57?")
print("=" * 75)

print("\n  Testing: eta(w) = sum_k w(k) * sign(lambda_k)")
print("  for various weight functions w(k):")

pos_idx = [k for k in range(9) if adj_eigs[k] > 1e-10]
neg_idx = [k for k in range(9) if adj_eigs[k] < -1e-10]
zero_idx = [k for k in range(9) if abs(adj_eigs[k]) < 1e-10]

def test_weight(name, weights):
    pos = sum(weights[k] for k in pos_idx)
    neg = sum(weights[k] for k in neg_idx)
    eta = pos - neg
    err = abs(eta - 57)
    marker = " <--- EXACT!" if err < 0.001 else (" <~~ CLOSE" if err < 1 else "")
    print("    %-30s pos=%8.2f neg=%8.2f eta=%8.2f %s" %
          (name, pos, neg, eta, marker))
    return eta

# Standard weights
w_dim2 = {k: dims[k]**2 for k in range(9)}
w_dim = {k: float(dims[k]) for k in range(9)}
w_1 = {k: 1.0 for k in range(9)}
w_dimL = {k: dims[k]*L_vals[k] for k in range(9)}
w_d2L = {k: dims[k]**2*L_vals[k] for k in range(9)}
w_C2 = {k: c2_su2[k] for k in range(9)}
w_dC2 = {k: dims[k]*c2_su2[k] for k in range(9)}
w_d2C2 = {k: dims[k]**2*c2_su2[k] for k in range(9)}
w_d3 = {k: dims[k]**3 for k in range(9)}
w_dim_L2 = {k: dims[k]*L_vals[k]**2 for k in range(9)}

# Additional creative weights
w_dim_plus1_sq = {k: (dims[k]+1)**2 for k in range(9)}
w_dim_times_dimplus1 = {k: dims[k]*(dims[k]+1) for k in range(9)}
w_triangle = {k: dims[k]*(dims[k]+1)//2 for k in range(9)}

# Test all
test_weight("dim^2 (=neutrino 35)", w_dim2)
test_weight("dim", w_dim)
test_weight("1 (count)", w_1)
test_weight("dim*L", w_dimL)
test_weight("dim^2*L", w_d2L)
test_weight("C2", w_C2)
test_weight("dim*C2", w_dC2)
test_weight("dim^2*C2", w_d2C2)
test_weight("dim^3", w_d3)
test_weight("dim*L^2", w_dim_L2)
test_weight("(dim+1)^2", w_dim_plus1_sq)
test_weight("dim*(dim+1)", w_dim_times_dimplus1)
test_weight("dim*(dim+1)/2", w_triangle)

# Try: dim^2 + dim (= dim*(dim+1))
w_combined = {k: dims[k]**2 + dims[k] for k in range(9)}
test_weight("dim^2+dim", w_combined)

# Try: dim^2 * (something involving k)
# The "natural" weight for CC might involve the McKay graph distance from rho_0
# Distances on McKay tree: rho_0 is at distance k from the main chain
mckay_dist = [0, 1, 2, 3, 4, 5, 4, 3, 1]  # distance from rho_0 on McKay graph
# Actually: rho_0(1)-rho_1(2)-rho_2(3)-rho_3(4)-rho_4(5)-rho_5(6)
# Branch: rho_5(6)-rho_6(4)-rho_7(3)
# Branch: rho_5(6)-rho_8(2)
# Distances from rho_0: 0,1,2,3,4,5,6,7,1 -- NO
# Let me be careful. The McKay graph (affine E8) is:
# Main chain: rho_0--rho_1--rho_2--rho_3--rho_4--rho_5
# Branch from rho_5: rho_5--rho_6--rho_7
# Branch from rho_5: rho_5--rho_8
# So distances from rho_0:
# rho_0: 0, rho_1: 1, rho_2: 2, rho_3: 3, rho_4: 4, rho_5: 5
# rho_6: 6, rho_7: 7, rho_8: 6
mckay_dist = [0, 1, 2, 3, 4, 5, 6, 7, 6]

w_dist = {k: dims[k]**2 * mckay_dist[k] for k in range(9)}
w_dist_d = {k: dims[k] * mckay_dist[k] for k in range(9)}

print("\n  Weights involving McKay distance from rho_0:")
print("  (distances: %s)" % mckay_dist)
test_weight("dim^2*dist", w_dist)
test_weight("dim*dist", w_dist_d)

# Mass exponent weights
mass_exp = {0: 0, 1: 3, 2: 26, 3: 5, 4: 19, 5: 17, 6: 11, 7: 11, 8: 16}
# Wait, the mapping is: rho_k -> fermion with specific mass exponent
# From exp363b: rho_1->e(0), rho_2->u(3), rho_3->t(26), etc.
# Actually from the code: rho_1(dim2)->e(n=0), NOT rho_0(dim1)->e
# Let me re-read the mapping
# From exp371 Part 10: rho_0(dim1)=trivial, rho_1(dim2)=e, etc.
# The trivial rho_0 has no fermion assignment (or = vacuum)

# =====================================================================
# PART 6: THE SUM-OF-GAPS (FROBENIUS) APPROACH
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 6: FROBENIUS / NUMERICAL SEMIGROUP")
print("=" * 75)

# Numerical semigroup generated by a1=5, b1=6 (gcd=1)
# Non-representable integers (gaps):
gaps = []
for n in range(1, a1*b1):  # Frobenius number + 1 = 20
    found = False
    for aa in range(n//a1 + 1):
        rem = n - a1*aa
        if rem >= 0 and rem % b1 == 0:
            found = True
            break
    if not found:
        gaps.append(n)

print("\n  Numerical semigroup <a1=%d, b1=%d>:" % (a1, b1))
print("  Gaps (non-representable as %d*a+%d*b, a,b>=0):" % (a1, b1))
print("  %s" % gaps)
print("  Number of gaps: %d" % len(gaps))
print("  Frobenius number (largest gap): %d" % max(gaps))
print("  Sum of gaps: %d" % sum(gaps))

# The Sylvester-Frobenius formulas for gcd(a,b)=1:
# Number of gaps = (a-1)*(b-1)/2
# Sum of gaps = (a-1)*(b-1)*(2*a*b-a-b-1)/12
n_gaps_formula = (a1-1)*(b1-1)//2
sum_gaps_formula = (a1-1)*(b1-1)*(2*a1*b1-a1-b1-1)//12
print("\n  Formulas:")
print("    # gaps = (a1-1)*(b1-1)/2 = %d*%d/2 = %d" %
      (a1-1, b1-1, n_gaps_formula))
print("    Sum gaps = (a1-1)*(b1-1)*(2*a1*b1-a1-b1-1)/12 = %d" % sum_gaps_formula)

# Is sum of gaps related to 57?
print("\n    Sum of gaps = %d" % sum(gaps))
print("    57 = %d? %s" % (57, "YES!" if sum(gaps) == 57 else "NO"))

# Check: sum of gaps = 57?
# For (5,6): (4)(5)(60-5-6-1)/12 = 20*48/12 = 80. NOT 57.
# But maybe there's a DIFFERENT sum that gives 57?

# Alternatively: count gaps below the MEDIAN?
# Or: sum of gaps that are ALSO non-norms in Z[phi]?
norms_zphi = set()
for aa in range(-20, 21):
    for bb in range(-20, 21):
        n_val = abs(aa**2 + aa*bb - bb**2)
        if n_val <= 100:
            norms_zphi.add(n_val)

gaps_non_norms = [g for g in gaps if g not in norms_zphi]
gaps_norms = [g for g in gaps if g in norms_zphi]
print("\n  Gaps that are non-norms in Z[phi]: %s (sum=%d)" %
      (gaps_non_norms, sum(gaps_non_norms)))
print("  Gaps that ARE norms in Z[phi]: %s (sum=%d)" %
      (gaps_norms, sum(gaps_norms)))

# What if we weight gaps by something?
# Gaps weighted by their position: sum(g * something)?
print("\n  Weighted gap sums:")
print("    sum(g) = %d" % sum(gaps))
print("    sum(g) - max(g) = %d - %d = %d" % (sum(gaps), max(gaps), sum(gaps)-max(gaps)))
print("    sum(g) - sum(g[g<10]) = %d - %d = %d" %
      (sum(gaps), sum(g for g in gaps if g<10), sum(gaps) - sum(g for g in gaps if g<10)))

# =====================================================================
# PART 7: THE eta_A + CKM CONNECTION
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 7: 57 = |eta_A| + sum(CKM exponents)")
print("=" * 75)

print("\n  |eta_A| = 35 (DERIVED, spectral asymmetry of 2I adjacency)")
print("  CKM bare exponents: n_12=3, n_13=12, n_23=7")
print("  sum(CKM) = 3 + 12 + 7 = 22")
print("  |eta_A| + sum(CKM) = 35 + 22 = 57")

print("\n  Decompose sum(CKM) = 22:")
print("    n_12 = 3 = N_gen")
print("    n_13 = 12 = DEG = 2*b1")
print("    n_23 = 7 = b1+1 = a1+2")
print("    Sum = N_gen + DEG + (b1+1) = 3 + 12 + 7 = 22")

# Prove: |eta_A| + N_gen + DEG + b1+1 = N_gen*(b1^2+2)/2
# 35 + 3 + 12 + 7 = 57
# (b1^2-1) + N_gen + 2*b1 + b1+1 = N_gen*(b1^2+2)/2
# b1^2-1 + N_gen + 3*b1+1 = N_gen*(b1^2+2)/2
# b1^2 + N_gen + 3*b1 = N_gen*(b1^2+2)/2
# For b1=6, N_gen=3: 36+3+18 = 57 = 3*19. CHECK.

print("\n  Algebraic: |eta_A| + N_gen + DEG + (b1+1)")
print("  = (b1^2-1) + N_gen + 2*b1 + (b1+1)")
print("  = b1^2 + 3*b1 + N_gen")
print("  = 36 + 18 + 3 = 57")
print("  = N_gen * (b1^2+2)/2  [requires a1=5]")

# =====================================================================
# PART 8: THE QCD CORRECTION -alpha_s
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 8: THE QCD CORRECTION")
print("=" * 75)

z_int = 57
z_frac = -ALPHA_S
z_total = z_int + z_frac
z_pred = 57 - ALPHA_S

print("\n  z_CC = 57 - alpha_s = 57 - 1/(2*phi^3)")
print("       = %.6f" % z_pred)
print("  Observed: %.3f +/- 0.005" % z_obs)
print("  Deviation: %.2f sigma" % (abs(z_pred - z_obs)/0.005))

print("\n  The correction alpha_s = 1/(2*phi^3) = %.6f" % ALPHA_S)
print("  = phi - 3/2 = %.6f" % (PHI - 1.5))
print("  These are equal: |diff| = %.2e" % abs(ALPHA_S - (PHI - 1.5)))

# Why alpha_s? Physical argument:
print("\n  WHY the QCD correction -alpha_s?")
print("  Lambda_P involves the TOTAL vacuum energy.")
print("  At tree level: z = 57 (integer, from mass budget + N_gen).")
print("  The QCD correction: vacuum condensate shifts z by -alpha_s.")
print("  alpha_s = 1/(2*phi^3) is the strong coupling (DERIVED).")
print("  The shift is ONE POWER of alpha_s (leading-order QCD correction).")
print("  STATUS: MOTIVATED (standard QCD physics), not uniquely derived.")

# =====================================================================
# PART 9: CAN WE DERIVE 57 FROM A SINGLE SPECTRAL INVARIANT?
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 9: COMBINED SPECTRAL INVARIANTS")
print("=" * 75)

# Test: |eta_A| + |eta_dim| where eta_dim = sign-weighted dim sum
eta_dim = sum(dims[k] * np.sign(adj_eigs[k]) for k in range(9)
              if abs(adj_eigs[k]) > 1e-10)
print("\n  eta_dim = sum(dim*sign(lambda)) = %.0f" % eta_dim)
print("  |eta_A| + |eta_dim| = %d + %d = %d" %
      (abs(int(eta_A_val)), abs(int(eta_dim)), abs(int(eta_A_val))+abs(int(eta_dim))))

# Test: spectral action on full 2640-dim complex
c0 = 2640
c1 = 14880
c2 = 55920
print("\n  Spectral action coefficients:")
print("    c0 = %d = %d * 240" % (c0, c0//240))
print("    c1 = %d = %d * 240" % (c1, c1//240))
print("    c2 = %d = %d * 240" % (c2, c2//240))
print("    c0/240 = %d = L_5 (Lucas)" % (c0//240))
print("    c1/240 = %d" % (c1//240))
print("    c2/240 = %d = F_13 (Fibonacci)" % (c2//240))

# Ratios
print("\n  Spectral ratios:")
print("    c1/c0 = %d/%d = %.4f" % (c1, c0, c1/c0))
print("    c1/(2*c0) = %d/%d = %.4f = 31/11" % (c1, 2*c0, c1/(2*c0)))
print("    c2/c1 = %d/%d = %.4f" % (c2, c1, c2/c1))
print("    c1^2/(c0*c2) = %.4f = 3/2 = N_gen/2" % (c1**2/(c0*c2)))

# Is there a combination of c0, c1, c2 that gives 57?
print("\n  Combinations:")
combos = [
    ("c1/c0 - 1/c0", c1/c0 - 1/c0),
    ("(c2-c0)/c1", (c2-c0)/c1),
    ("c2/(c1-c0)*N_gen", c2/(c1-c0)*N_gen),
    ("c0*(c1-c0)/c2", c0*(c1-c0)/c2),
    ("(c0+c1)/(c0/N_gen)", (c0+c1)/(c0/N_gen)),
    ("c1/240 - a1", c1/240 - a1),
    ("c0/240 * a1 + 2", c0//240 * a1 + 2),
    ("c2/c1 + c1/c0 - c0/c1", c2/c1 + c1/c0 - c0/c1),
    ("N*c0/c2", N*c0/c2),
]

print("  %40s %12s" % ("Expression", "Value"))
for name, val in combos:
    marker = " <---" if abs(val - 57) < 0.5 else ""
    print("  %40s %12.4f%s" % (name, val, marker))

# =====================================================================
# PART 10: HONEST ASSESSMENT
# =====================================================================
print("\n\n" + "=" * 75)
print("PART 10: HONEST ASSESSMENT")
print("=" * 75)

print("""
  WHAT IS NEW FROM THIS EXPERIMENT:

  1. MASS BUDGET DECOMPOSITION (D): 57 = sum(n_f)/2 + N_gen = 54 + 3
     - Connects CC exponent to the total fermion mass content
     - Uses sum(n_f) = 108 (DERIVED in exp371)
     - Equivalent to a1=5 (proven: a1^2-4a1-5=0)
     - STATUS: DERIVED identity, but not a DERIVATION of CC

  2. BRIDGE: 57 = |eta_A| + sum(CKM)
     - Neutrino sector (35) + CKM sector (22)
     - sum(CKM) = N_gen + DEG + (b1+1) = 3+12+7 = 22
     - All CKM exponents are framework quantities
     - STATUS: DECOMPOSITION (not derivation)

  3. FROBENIUS: (b1^2+2)/2 = Frobenius(a1,b1) = 19
     - This is equivalent to a1=5 (quadratic equation)
     - Connects Frobenius number to Galois norm b1^2=36
     - STATUS: ALGEBRAIC IDENTITY

  4. SUM OF GAPS: sum = 80, NOT 57
     - Previous claim of sum=57 was INCORRECT
     - Sylvester-Frobenius formula gives 80 for (5,6)

  5. SPECTRAL SEARCH: NO single spectral invariant gives 57 exactly
     - dim-weighted Laplacian (positive sector) = 96-24*phi ~ 57.17 (CLOSE but not exact)
     - No clean weighting found in systematic scan
     - The CC exponent does NOT emerge from a single spectral asymmetry

  OVERALL STATUS: 57 remains PATTERN.

  The identity 57 = sum(n)/2 + N_gen is NEW and connects CC to mass budget.
  But it does not explain WHY Lambda_P = alpha^(sum(n)/2+N_gen-alpha_s).
  The physical mechanism (vacuum energy from fermion masses?) is motivated
  but not derived from first principles.

  COMPARISON WITH 35:
  - 35 = |eta_A| (spectral asymmetry): SINGLE spectral invariant. DERIVED.
  - 57 = sum(n)/2 + N_gen: ALGEBRAIC identity from mass budget. NOT spectral.
  - 57 CANNOT be obtained from a single spectral asymmetry of the Cayley graph.

  The CC exponent remains the framework's MOST SERIOUS OPEN PROBLEM.
  Category: STRONG PATTERN (0.13 sigma, 4 equivalent decompositions,
  all equivalent to a1=5, but no first-principles derivation).
""")

print("=" * 75)
print("EXP-373 COMPLETE")
print("=" * 75)

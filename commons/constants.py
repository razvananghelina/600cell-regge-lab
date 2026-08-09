"""
Fundamental constants derived from a1 = 5.
"""
import numpy as np

# Fundamental integer
a1 = 5
b1 = 6       # = a1 + 1
N = 120      # = 2 * a1!  (vertices of 600-cell, |2I|)
h = 30       # = a1 * b1  (Coxeter number of E8)
rank_E8 = 8
dim_E8 = 248 # = 2*N + rank_E8
N_gen = 3    # = number of generations (from McKay p4/p2)
degree = 12  # = 2 * b1 (vertex degree of 600-cell)
d_ST = 4     # legacy name: real dim(C^2/2I)=a1-1; spacetime use is OPEN

# Golden ratio
PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2  # = -1/PHI
SQRT5 = np.sqrt(5)

# Coupling constants
alpha_s = 1 / (2 * PHI**3)
sin2_tW = b1 / (a1**2 + 1)  # = 6/26

# Alpha from quadratic: 2*pi*a^2 - 4*a1*phi^4*a + 1 = 0
_disc = 16 * a1**2 * PHI**8 - 8 * np.pi
alpha_em = (4 * a1 * PHI**4 - np.sqrt(_disc)) / (4 * np.pi)
inv_alpha = 1 / alpha_em
ln_inv_alpha = np.log(inv_alpha)

# Seeley-DeWitt per root pair
c1_per_2N = dim_E8 // 4  # = 62

# ============================================================
# 600-CELL LAPLACIAN SPECTRUM
# Verified against numerical diagonalization: Tr(L) = 1440 = N*degree
# 9 distinct eigenvalues, one per irrep of 2I
# ============================================================

# Eigenvalues: exact algebraic forms and multiplicities
# Format: (eigenvalue, multiplicity, "exact_form", galois_type)
SPECTRUM_600CELL = [
    (0,              1,  "0",           "fixed"),
    (12 - 6*PHI,     4,  "12-6phi",     "broken"),  # Galois pair A (phys)
    (10 - 2*SQRT5,   9,  "10-2sqrt5",   "broken"),  # Galois pair B (phys)
    (9,             16,  "9",           "fixed"),
    (12,            25,  "12",          "fixed"),
    (14,            36,  "14",          "fixed"),
    (10 + 2*SQRT5,   9,  "10+2sqrt5",   "broken"),  # Galois pair B (dark)
    (15,            16,  "15",          "fixed"),
    (6 + 6*PHI,      4,  "6+6phi",      "broken"),  # Galois pair A (dark)
]

# Named eigenvalues for convenience
L0 = 0
L1 = 12 - 6*PHI       # = b1/phi^2      = 2.2918...  (Galois pair A, physical)
L2 = 10 - 2*SQRT5     # = 2*D'^2        = 5.5279...  (Galois pair B, physical)
L3 = 9
L4 = 12
L5 = 14
L6 = 10 + 2*SQRT5     # = 2*D^2         = 14.4721... (Galois pair B, dark)
L7 = 15
L8 = 6 + 6*PHI        # = b1*phi^2      = 15.7082... (Galois pair A, dark)

# ---- EIGENVALUE-TQFT IDENTITIES (discovered exp521) ----
# L1 = b1/phi^2,  L8 = b1*phi^2    (fusion norm scaled by quantum dim)
# L2 = 2*D'^2,    L6 = 2*D^2       (twice the TQFT total quantum dim squared!)
# These connect the 600-cell graph spectrum directly to the SU(2)_3 TQFT.

# Galois pair products (EXACT, rational = Galois norms)
GALOIS_NORM_A = L1 * L8   # = b1^2         = 36
GALOIS_NORM_B = L2 * L6   # = d_ST^2 * a1  = 80  (= 4*D^2*D'^2)

# Galois pair sums (EXACT, rational)
GALOIS_SUM_A = L1 + L8    # = 3*b1 = 18    (since 1/phi^2+phi^2 = 3)
GALOIS_SUM_B = L2 + L6    # = d_ST*a1 = 20 (= 2*(D^2+D'^2) = 4*a1)

# Galois pair ratios (EXACT)
GALOIS_RATIO_A = L8 / L1  # = phi^4
GALOIS_RATIO_B = L6 / L2  # = phi^2 = D^2/D'^2

# Eigenvalue differences (EXACT)
GALOIS_DIFF_A = L8 - L1   # = 6*sqrt5 = b1*sqrt(a1)
GALOIS_DIFF_B = L6 - L2   # = 4*sqrt5 = d_ST*sqrt(a1)

# Mode counts
N_GALOIS_BROKEN = 4 + 9 + 9 + 4   # = 26 = a1^2 + 1
N_GALOIS_FIXED = 1 + 16 + 25 + 36 + 16  # = 94

# TQFT quantum dimensions (SU(2)_3, k=3, k+2=a1=5)
D2_PHYS = a1 + SQRT5       # = 5+sqrt5 = 7.236... (D^2 physical)
D2_DARK = a1 - SQRT5       # = 5-sqrt5 = 2.764... (D'^2 dark)
DELTA_D2 = D2_PHYS - D2_DARK  # = 2*sqrt5

# Dark quantum dimension
d_dark = 1 / PHI           # = |sigma(phi)| = phi-1

# DM abundance (STRUCTURAL)
Omega_DM = b1 - d_dark     # = 7-phi = 5.382...

# ============================================================
# DIRAC SPECTRUM ON S^3/2I (POINCARÉ HOMOLOGY SPHERE)
# Computed via McKay recursion on extended E8 Dynkin diagram.
# Eigenvalue +(k+3/2) has mult = n_{rho_1, k+1}.
# Verified: exp522 (convergent at 50+ eigenvalues).
# ============================================================

# McKay graph adjacency (extended E8): 0-1-2-3-4-5-6-7, plus 5-8
MCKAY_EDGES = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(5,8)]
IRREP_DIMS_2I = [1, 2, 3, 4, 5, 6, 4, 2, 3]  # 9 irreps of 2I

# First Dirac eigenvalues on S^3/2I (exact):
# lambda_1 = 5/2 = a1/2  (mult 1, from V_2 = rho_1)
# lambda_2 = 25/2 = a1^2/2  (mult 1, from V_12 containing rho_1)
# Ratio lambda_2/lambda_1 = a1 = 5
DIRAC_LAMBDA_1 = a1 / 2        # = 5/2 = 2.5
DIRAC_LAMBDA_2 = a1**2 / 2     # = 25/2 = 12.5
DIRAC_GAP_RATIO = a1            # lambda_2/lambda_1 = 5

# ============================================================
# BOOTSTRAP RESULTS (exp523)
# ============================================================

# Alpha existence: dark sector has NO real EM coupling
# Alpha equation: 2*pi*x^2 - 4*a1*d^4*x + 1 = 0
# Real solutions require d >= d_crit = (pi/(2*a1^2))^{1/8}
# phi > d_crit > 1/phi: physical has EM, dark doesn't
ALPHA_COEFF = 4 * a1 * PHI**4         # = 70 + 30*sqrt5 = 137.082...
ALPHA_COEFF_DARK = 4 * a1 / PHI**4    # = 70 - 30*sqrt5 = 2.918... (disc < 0!)
ALPHA_COEFF_SUM = 4 * a1 * 7          # = 140 (Galois trace, L_4=7)
ALPHA_COEFF_PRODUCT = 16 * a1**2      # = 400 (Galois norm)

# Galois norms of couplings
GALOIS_NORM_OMEGA = 41   # = Omega * sigma(Omega) = b1^2+a1 = a1^2+3a1+1 (prime)
# Omega minimal polynomial: x^2 - 13x + 41 = 0, discriminant = a1

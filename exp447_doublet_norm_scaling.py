"""
exp447_doublet_norm_scaling.py
==============================
Doublet Wilson-line products and their Galois norms.

OBSERVATION: The SU(2)_L doublet products z_up * z_down have norms
|N| = 1, a_1, (4a_1-1)^2 across 3 generations.

GOALS:
1. Verify the diagonal (same-generation) doublet products
2. Compute the FULL 3x3 matrix of products z_{u_i} * z_{d_j}
3. Check for CKM connections
4. Look for sum rules on the product norms
5. Investigate the generation scaling

Author: Razvan-Constantin Anghelina
Date: 2026-03-23
"""

import numpy as np
from itertools import product as iter_product

# ============================================================
# PART 0: SETUP
# ============================================================
a1 = 5
phi = (1 + np.sqrt(a1)) / 2  # golden ratio
phi_conj = (1 - np.sqrt(a1)) / 2  # Galois conjugate

def galois_norm(a, b):
    """N(a + b*phi) = a^2 + ab - b^2"""
    return a**2 + a*b - b**2

def galois_trace(a, b):
    """Tr(a + b*phi) = 2a + b"""
    return 2*a + b

def z_value(a, b):
    """Numerical value of a + b*phi"""
    return a + b * phi

def z_conj(a, b):
    """Galois conjugate a + b*phi'"""
    return a + b * phi_conj

def multiply_zphi(a1, b1, a2, b2):
    """Multiply (a1+b1*phi)*(a2+b2*phi) = (a1*a2+b1*b2) + (a1*b2+a2*b1+b1*b2)*phi
    Using phi^2 = phi + 1"""
    a_new = a1*a2 + b1*b2
    b_new = a1*b2 + a2*b1 + b1*b2
    return a_new, b_new

# Wilson line quantum numbers (a, b) for each fermion
# m_f = m_e * phi^(5a + 6b + corrections)
fermions = {
    'e':   (0, 0),
    'mu':  (1, 1),
    'tau': (1, 2),
    'u':   (3, -2),
    'd':   (1, 0),
    's':   (1, 1),
    'c':   (2, 1),
    'b':   (-1, 4),
    't':   (4, 1),
}

up_quarks = ['u', 'c', 't']
down_quarks = ['d', 's', 'b']
leptons = ['e', 'mu', 'tau']

print("=" * 70)
print("EXP447: DOUBLET NORM SCALING")
print("=" * 70)

# ============================================================
# PART 1: INDIVIDUAL NORMS
# ============================================================
print("\n" + "=" * 70)
print("PART 1: INDIVIDUAL FERMION WILSON LINES")
print("=" * 70)

print(f"\n{'Fermion':<8} {'(a,b)':<12} {'z_f':<12} {'z_f_conj':<12} {'N(z_f)':<10} {'Tr(z_f)':<10}")
print("-" * 64)
for f in ['e', 'mu', 'tau', 'u', 'd', 's', 'c', 'b', 't']:
    a, b = fermions[f]
    z = z_value(a, b)
    zc = z_conj(a, b)
    N = galois_norm(a, b)
    T = galois_trace(a, b)
    print(f"{f:<8} ({a:+d},{b:+d})     {z:+10.4f}  {zc:+10.4f}  {N:+8d}  {T:+8d}")

# ============================================================
# PART 2: DIAGONAL DOUBLET PRODUCTS (same generation)
# ============================================================
print("\n" + "=" * 70)
print("PART 2: DIAGONAL DOUBLET PRODUCTS (same generation)")
print("=" * 70)

generations = [('u', 'd', 'Gen 1'), ('c', 's', 'Gen 2'), ('t', 'b', 'Gen 3')]

for up, down, label in generations:
    a_u, b_u = fermions[up]
    a_d, b_d = fermions[down]

    # Product in Z[phi]
    a_prod, b_prod = multiply_zphi(a_u, b_u, a_d, b_d)
    z_prod = z_value(a_prod, b_prod)
    N_prod = galois_norm(a_prod, b_prod)

    # Verify: N(product) = N(up) * N(down)
    N_up = galois_norm(a_u, b_u)
    N_down = galois_norm(a_d, b_d)

    print(f"\n{label}: z_{up} * z_{down} = ({a_u}+{b_u}phi)*({a_d}+{b_d}phi)")
    print(f"  = {a_prod} + {b_prod}*phi = {z_prod:.4f}")
    print(f"  N(product) = {N_prod}")
    print(f"  |N(up)|*|N(down)| = {abs(N_up)}*{abs(N_down)} = {abs(N_up)*abs(N_down)}")
    print(f"  Verify: N(product) = N(up)*N(down) = {N_up}*{N_down} = {N_up*N_down}  [{'OK' if N_prod == N_up*N_down else 'FAIL'}]")

print(f"\n--- SCALING PATTERN ---")
norms_diag = []
for up, down, label in generations:
    a_u, b_u = fermions[up]
    a_d, b_d = fermions[down]
    N = abs(galois_norm(a_u, b_u) * galois_norm(a_d, b_d))
    norms_diag.append(N)
    print(f"  {label}: |N| = {N}")

print(f"\n  Ratios: {norms_diag[0]} : {norms_diag[1]} : {norms_diag[2]}")
print(f"  Expected: 1 : a_1 : (4a_1-1)^2 = 1 : {a1} : {(4*a1-1)**2}")
print(f"  Match: {norms_diag == [1, a1, (4*a1-1)**2]}")

# ============================================================
# PART 3: FULL 3x3 PRODUCT MATRIX
# ============================================================
print("\n" + "=" * 70)
print("PART 3: FULL 3x3 PRODUCT MATRIX z_{up_i} * z_{down_j}")
print("=" * 70)

print("\n--- Wilson line products (a + b*phi representation) ---")
print(f"{'':>8}", end="")
for d in down_quarks:
    print(f"{'z_'+d:>20}", end="")
print()

product_matrix_zphi = {}
product_matrix_norm = np.zeros((3, 3), dtype=int)
product_matrix_trace = np.zeros((3, 3), dtype=int)

for i, up in enumerate(up_quarks):
    print(f"z_{up:>4}", end="  ")
    for j, down in enumerate(down_quarks):
        a_u, b_u = fermions[up]
        a_d, b_d = fermions[down]
        a_p, b_p = multiply_zphi(a_u, b_u, a_d, b_d)
        N_p = galois_norm(a_p, b_p)
        T_p = galois_trace(a_p, b_p)
        product_matrix_zphi[(up, down)] = (a_p, b_p)
        product_matrix_norm[i, j] = N_p
        product_matrix_trace[i, j] = T_p
        print(f"  ({a_p:+d}{b_p:+d}phi)={N_p:+4d}", end="")
    print()

print("\n--- Norm matrix N(z_up * z_down) ---")
print(f"{'':>8}", end="")
for d in down_quarks:
    print(f"{d:>8}", end="")
print()
for i, up in enumerate(up_quarks):
    print(f"  {up:>4}  ", end="")
    for j in range(3):
        print(f"{product_matrix_norm[i,j]:>8}", end="")
    print()

print("\n--- |Norm| matrix ---")
abs_norm = np.abs(product_matrix_norm)
print(f"{'':>8}", end="")
for d in down_quarks:
    print(f"{d:>8}", end="")
print()
for i, up in enumerate(up_quarks):
    print(f"  {up:>4}  ", end="")
    for j in range(3):
        print(f"{abs_norm[i,j]:>8}", end="")
    print()

# ============================================================
# PART 4: CKM CONNECTION TEST
# ============================================================
print("\n" + "=" * 70)
print("PART 4: CKM CONNECTION TEST")
print("=" * 70)

# CKM experimental values |V_ij|^2
V_CKM_sq = np.array([
    [0.97373**2, 0.2243**2, 0.00382**2],   # u row
    [0.221**2,   0.987**2,  0.0410**2],     # c row
    [0.0086**2,  0.0415**2, 1.014**2],      # t row (V_tb from unitarity)
])

# CKM exponents: V_ij ~ phi^(-n_ij)
n_CKM = np.array([
    [0,  3, 12],    # n_ud=0 (diagonal), n_us=3, n_ub=12
    [3,  0,  7],    # n_cd=3, n_cs=0, n_cb=7
    [12, 7,  0],    # n_td=12, n_ts=7, n_tb=0
])

print("\nDoes |N(z_up*z_down)| correlate with CKM?")
print("\n  CKM |V|^2 matrix:")
for i, up in enumerate(up_quarks):
    print(f"    {up}: ", [f"{V_CKM_sq[i,j]:.6f}" for j in range(3)])

print("\n  |Norm| matrix (flat):")
for i, up in enumerate(up_quarks):
    print(f"    {up}: ", [f"{abs_norm[i,j]}" for j in range(3)])

# Test: is there a function f such that |V_ij|^2 ~ f(|N_ij|)?
print("\n  Correlation test:")
x = abs_norm.flatten()
y = V_CKM_sq.flatten()
# Remove diagonal (self-mixing)
mask = ~np.eye(3, dtype=bool).flatten()
x_off = x[mask]
y_off = y[mask]
if len(x_off) > 2:
    corr = np.corrcoef(x_off, y_off)[0, 1]
    corr_log = np.corrcoef(np.log(x_off + 1), np.log(y_off + 1e-10))[0, 1]
    print(f"  Pearson (off-diag |N| vs |V|^2): {corr:.4f}")
    print(f"  Pearson (log|N| vs log|V|^2):    {corr_log:.4f}")

# Test: N_ij vs n_ij (CKM exponents)
print("\n  CKM exponents n_ij vs norm products:")
for i, up in enumerate(up_quarks):
    for j, down in enumerate(down_quarks):
        if i != j or True:  # show all
            print(f"    ({up},{down}): n_CKM={n_CKM[i,j]:>3}, |N(z_u*z_d)|={abs_norm[i,j]:>4}, ratio={abs_norm[i,j]/(n_CKM[i,j]+0.001):.3f}")

print("\n  VERDICT: Norm matrix is RANK-1 (outer product of individual |N|)")
print("  N(z_up*z_down) = N(z_up)*N(z_down) by multiplicativity")
print("  No CKM structure beyond individual norms")
print("  STATUS: NEGATIVE for direct CKM connection")

# ============================================================
# PART 5: SUM RULES ON PRODUCT MATRIX
# ============================================================
print("\n" + "=" * 70)
print("PART 5: SUM RULES ON PRODUCT MATRIX")
print("=" * 70)

# Row sums
print("\n--- Row sums of N(z_up*z_down) ---")
for i, up in enumerate(up_quarks):
    rs = sum(product_matrix_norm[i, :])
    print(f"  {up}: sum = {rs}")

# Column sums
print("\n--- Column sums ---")
for j, down in enumerate(down_quarks):
    cs = sum(product_matrix_norm[:, j])
    print(f"  {down}: sum = {cs}")

# Total
total = product_matrix_norm.sum()
print(f"\n  Grand total: {total}")

# Trace (diagonal sum)
diag_sum = sum(product_matrix_norm[i, i] for i in range(3))
print(f"  Trace (diagonal): {diag_sum}")

# Off-diagonal sum
off_diag = total - diag_sum
print(f"  Off-diagonal sum: {off_diag}")

# Determinant
det_N = int(round(np.linalg.det(product_matrix_norm.astype(float))))
print(f"  Determinant: {det_N}")
print(f"  (Expected 0 for rank-1 matrix)")

# Trace of |N| matrix
print(f"\n  Tr(|N|) = {sum(abs_norm[i,i] for i in range(3))}")
print(f"  = 1 + {a1} + {(4*a1-1)**2} = {1 + a1 + (4*a1-1)**2}")

# ============================================================
# PART 6: TRACE MATRIX
# ============================================================
print("\n" + "=" * 70)
print("PART 6: TRACE MATRIX Tr(z_up*z_down)")
print("=" * 70)

print(f"\n{'':>8}", end="")
for d in down_quarks:
    print(f"{d:>8}", end="")
print()
for i, up in enumerate(up_quarks):
    print(f"  {up:>4}  ", end="")
    for j in range(3):
        print(f"{product_matrix_trace[i,j]:>8}", end="")
    print()

trace_total = product_matrix_trace.sum()
trace_diag = sum(product_matrix_trace[i, i] for i in range(3))
print(f"\n  Grand total Tr: {trace_total}")
print(f"  Diagonal Tr: {trace_diag}")

# ============================================================
# PART 7: z_PLANCK CONNECTION
# ============================================================
print("\n" + "=" * 70)
print("PART 7: z_PLANCK = 4*phi^2 = 4 + 4*phi")
print("=" * 70)

z_Planck_a, z_Planck_b = 4, 4  # 4 + 4*phi
N_Planck = galois_norm(z_Planck_a, z_Planck_b)
Tr_Planck = galois_trace(z_Planck_a, z_Planck_b)

print(f"\n  z_Planck = 4 + 4*phi = {z_value(4, 4):.6f}")
print(f"  N(z_Planck) = {N_Planck} = (a_1-1)^2 = {(a1-1)**2}")
print(f"  Tr(z_Planck) = {Tr_Planck} = 2*b_1 = {2*(a1+1)}")

print(f"\n  TRIANGLE OF 16 = (a_1-1)^2:")
print(f"    (1) N(z_Planck) = {N_Planck}  [Hierarchy: m_e/M_P = alpha^(4phi^2)]")
print(f"    (2) Tr(A^2) = {(a1-1)**2}     [Higgs: m_H^2/m_W^2 = phi^2 - 16*alpha*phi]")
print(f"    (3) p_2 = {(a1-1)**2}          [McKay power sum]")
print(f"    (4) |N(diagonal)|_sum = {1+a1+(4*a1-1)**2} = 1+{a1}+{(4*a1-1)**2}")
print(f"        (NOT 16, but interesting: {1+a1+(4*a1-1)**2} = {1+a1+(4*a1-1)**2})")

# Check: Tr(|N|_diagonal) = 1 + 5 + 361 = 367 = ?
TrN = 1 + a1 + (4*a1-1)**2
print(f"\n  1 + a_1 + (4a_1-1)^2 = {TrN}")
print(f"  = 1 + a_1 + 16*a_1^2 - 8*a_1 + 1 = 16*a_1^2 - 7*a_1 + 2")
print(f"  = 16*{a1**2} - 7*{a1} + 2 = {16*a1**2 - 7*a1 + 2}")

# ============================================================
# PART 8: GENERATION-DEPENDENT INVARIANTS
# ============================================================
print("\n" + "=" * 70)
print("PART 8: GENERATION-DEPENDENT INVARIANTS")
print("=" * 70)

print("\n--- Individual |N| by sector ---")
print("  Up quarks:  ", [abs(galois_norm(*fermions[q])) for q in up_quarks])
print("  Down quarks:", [abs(galois_norm(*fermions[q])) for q in down_quarks])
print("  Leptons:    ", [abs(galois_norm(*fermions[q])) for q in leptons])

print("\n--- Up quark norm pattern ---")
N_up = [galois_norm(*fermions[q]) for q in up_quarks]
print(f"  N(u)={N_up[0]}, N(c)={N_up[1]}, N(t)={N_up[2]}")
print(f"  |N|: {abs(N_up[0])}, {abs(N_up[1])}, {abs(N_up[2])}")
print(f"  Pattern: 1, a_1, 4*a_1-1 = 1, {a1}, {4*a1-1}")

print("\n--- Down quark norm pattern ---")
N_down = [galois_norm(*fermions[q]) for q in down_quarks]
print(f"  N(d)={N_down[0]}, N(s)={N_down[1]}, N(b)={N_down[2]}")
print(f"  |N|: {abs(N_down[0])}, {abs(N_down[1])}, {abs(N_down[2])}")
print(f"  Pattern: 1, 1, 4*a_1-1 = 1, 1, {4*a1-1}")

print("\n--- Lepton norm pattern ---")
N_lep = [galois_norm(*fermions[q]) for q in leptons]
print(f"  N(e)={N_lep[0]}, N(mu)={N_lep[1]}, N(tau)={N_lep[2]}")
print(f"  Sum = {sum(N_lep)} (anomaly cancellation)")

# ============================================================
# PART 9: DEEPER INVARIANTS
# ============================================================
print("\n" + "=" * 70)
print("PART 9: DEEPER INVARIANTS")
print("=" * 70)

# Check: is there a Casimir-like invariant per generation?
print("\n--- Sum N(z_up) + N(z_down) per generation ---")
gen_pairs = [('u', 'd'), ('c', 's'), ('t', 'b')]
for up, down in gen_pairs:
    Nu = galois_norm(*fermions[up])
    Nd = galois_norm(*fermions[down])
    print(f"  ({up},{down}): N_up + N_down = {Nu} + {Nd} = {Nu+Nd}")

print("\n--- N(z_up) - N(z_down) per generation ---")
for up, down in gen_pairs:
    Nu = galois_norm(*fermions[up])
    Nd = galois_norm(*fermions[down])
    print(f"  ({up},{down}): N_up - N_down = {Nu} - {Nd} = {Nu-Nd}")

print("\n--- N(z_up) * N(z_down) per generation ---")
for up, down in gen_pairs:
    Nu = galois_norm(*fermions[up])
    Nd = galois_norm(*fermions[down])
    print(f"  ({up},{down}): N_up * N_down = {Nu} * {Nd} = {Nu*Nd}")

# z_Planck as sum
print(f"\n--- z_Planck decomposition ---")
sum_up = sum(z_value(*fermions[q]) for q in up_quarks)
sum_down = sum(z_value(*fermions[q]) for q in down_quarks)
sum_lep = sum(z_value(*fermions[q]) for q in leptons)
sum_all = sum(z_value(*fermions[q]) for q in fermions)
print(f"  Sum z (up quarks) = {sum_up:.6f}")
print(f"  Sum z (down quarks) = {sum_down:.6f}")
print(f"  Sum z (leptons) = {sum_lep:.6f}")
print(f"  Sum z (all) = {sum_all:.6f}")
print(f"  z_Planck = {z_value(4, 4):.6f}")
print(f"  Ratio sum_all/z_Planck = {sum_all/z_value(4, 4):.6f}")

# Sum of (a,b) components
sum_a = sum(fermions[f][0] for f in fermions)
sum_b = sum(fermions[f][1] for f in fermions)
print(f"\n  Sum(a) = {sum_a} = dim(G) = 12")
print(f"  Sum(b) = {sum_b} = rank(E8) = 8")
print(f"  Sum(a) + Sum(b) = {sum_a + sum_b} = dim(G)+rank(E8) = 20 = 4*a_1")

# Is 4*a_1 = 20 = sum(a)+sum(b) related to z_Planck?
# z_Planck = 4+4phi, components sum = 4+4 = 8
# But sum(a)+sum(b) = 12+8 = 20 = 4*a_1
print(f"  4*a_1 = {4*a1}")
print(f"  4*phi^2 = 4*(phi+1) = 4*phi + 4, trace = 4+4+4 = 12 = sum(a)? No.")
print(f"  z_Planck trace = {Tr_Planck} = 12 = sum(a) = dim(G)!")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
1. DIAGONAL SCALING: |N(z_up*z_down)| = 1 : {a1} : {(4*a1-1)**2}
   = 1 : a_1 : (4a_1-1)^2. VERIFIED.
   This is just |N_up|*|N_down| by multiplicativity.

2. FULL MATRIX: Rank 1 (outer product of |N_up| and |N_down|).
   NO independent CKM structure. STATUS: NEGATIVE for CKM.

3. TRIANGLE OF 16: N(z_Planck) = Tr(A^2) = p_2 = (a_1-1)^2.
   Three different objects sharing same value. STRUCTURAL.

4. TRACE CONNECTION: Tr(z_Planck) = 12 = sum(a) = dim(G).
   The sum of all a-quantum numbers = trace of Planck Wilson line.

5. SUM RULE: sum(a) + sum(b) = 12 + 8 = 20 = 4*a_1.
   This is dim(G) + rank(E8) = 4*a_1.

STATUS: PATTERN (scaling), NEGATIVE (CKM), STRUCTURAL (triangle 16)
""")

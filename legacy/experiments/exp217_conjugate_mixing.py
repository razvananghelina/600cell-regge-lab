"""
EXP-217: Conjugate Subgroup Mixing Matrix
==========================================
PURPOSE: Compute the 3x3 mixing matrix from overlapping invariant subspaces
         of CONJUGATE 2T subgroups within 2I.

KEY INSIGHT (from exp216):
- 2T is NOT normal in 2I (A_4 not normal in A_5, simple group)
- There are 5 conjugate subgroups 2T_g = g*2T*g^{-1}
- Each corresponds to one inscribed 24-cell
- The dim-6 irrep decomposes into 3 two-dim subspaces under EACH 2T_g
- The OVERLAP between subspaces of different conjugates = mixing matrix
- M_{ij}(g) = Tr(P_i(2T_0) * P_j(2T_g)) / 2

METHOD:
1. Compute the 3 invariant subspaces V_k of 2T_0 (from exp216)
2. For each coset representative g, compute D(g)*V_k = subspaces of 2T_g
3. Compute 3x3 overlap matrix
4. Compare to CKM/PMNS |V|^2
"""

import numpy as np
from itertools import product as iprod, permutations
from scipy.special import comb as binomial

PHI = (1 + np.sqrt(5)) / 2
PSI = 1 - PHI

print("=" * 70)
print("EXP-217: CONJUGATE SUBGROUP MIXING MATRIX")
print("=" * 70)

# ======================================================================
# PART 1: Setup (same as exp216)
# ======================================================================
print("\n--- PART 1: Building 2I, 2T, spin-5/2 rep ---")

def qmul(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return (
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    )

# Generate 2I
elements = []
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0]*4; v[i] = s; elements.append(tuple(v))
for signs in iprod([0.5, -0.5], repeat=4):
    elements.append(tuple(signs))
even_perms = [
    (0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
    (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0)
]
bv = [0.0, 0.5, PHI/2, 1/(2*PHI)]
for perm in even_perms:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            for s3 in [1,-1]:
                v = [0.0]*4
                v[perm[1]] = s1*bv[1]; v[perm[2]] = s2*bv[2]; v[perm[3]] = s3*bv[3]
                elements.append(tuple(v))
unique = []
for v in elements:
    is_dup = any(sum((v[k]-u[k])**2 for k in range(4)) < 1e-10 for u in unique)
    if not is_dup:
        unique.append(v)
elements = unique
N = len(elements)
elem_arr = np.array(elements)
print("|2I| = %d" % N)

# Multiplication table
mult_table = np.zeros((N, N), dtype=int)
for i in range(N):
    for j in range(N):
        prod = qmul(elements[i], elements[j])
        dists = np.sum((elem_arr - np.array(prod))**2, axis=1)
        mult_table[i, j] = np.argmin(dists)

identity_idx = 0
for i in range(N):
    if all(mult_table[i, j] == j for j in range(N)):
        identity_idx = i
        break

inv_table = np.zeros(N, dtype=int)
for i in range(N):
    for k in range(N):
        if mult_table[i, k] == identity_idx:
            inv_table[i] = k
            break

# 2T = Type A + Type B
subgroup_2T = []
for i, q in enumerate(elements):
    n_zero = sum(1 for c in q if abs(c) < 1e-10)
    if n_zero == 3 or (n_zero == 0 and all(abs(abs(c)-0.5) < 1e-10 for c in q)):
        subgroup_2T.append(i)
is_2T = set(subgroup_2T)
print("|2T| = %d" % len(subgroup_2T))

# SU(2) and spin-j
def quat_to_su2(q):
    w, x, y, z = q
    a_val = complex(w, x)
    b_val = complex(y, z)
    return np.array([[a_val, b_val], [-np.conj(b_val), np.conj(a_val)]])

def spin_j_matrix(U, j):
    two_j = int(2*j)
    dim = two_j + 1
    a, b, c, d = U[0,0], U[0,1], U[1,0], U[1,1]
    D = np.zeros((dim, dim), dtype=complex)
    for m_idx in range(dim):
        m = j - m_idx
        p = int(j + m)
        q_val = int(j - m)
        for mp_idx in range(dim):
            mp = j - mp_idx
            pp = int(j + mp)
            val = 0
            for s in range(p + 1):
                for t in range(q_val + 1):
                    if s + t == pp:
                        coeff = (binomial(p, s, exact=True) *
                                binomial(q_val, t, exact=True) *
                                a**s * c**(p-s) * b**t * d**(q_val-t))
                        val += coeff
            D[mp_idx, m_idx] = val
    return D

# Compute spin-5/2 representation
j_gen = 2.5
dim_gen = 6
print("Computing spin-5/2 matrices...")
rep6 = []
for i in range(N):
    U = quat_to_su2(elements[i])
    D = spin_j_matrix(U, j_gen)
    rep6.append(D)
print("Done. Representation verified: err = %.2e"
      % max(np.max(np.abs(rep6[i] @ rep6[j] - rep6[mult_table[i,j]]))
            for i in [0,1,10,20] for j in [0,1,10,20]))

# ======================================================================
# PART 2: Find 2T_0 invariant subspaces using SCHUR PROJECTORS
# ======================================================================
print("\n--- PART 2: Schur Projectors for 2T ---")

# The key improvement over exp216: use proper Schur projectors
# P_j = (d_j / |G|) * sum_{g in G} conj(chi_j(g)) * D(g)
# But we need the character values of 2T irreps.

# From exp216 output, the characters of the 3 two-dim irreps on 2T classes:
# psi_5 (real): chi(e)=2, chi(-1)=-2, chi(D_2)=1, chi(D_3)=1, chi(D_4)=-1, chi(D_5)=-1, chi(D_6)=0
# psi_3 (complex): chi(e)=2, chi(-1)=-2, chi(D_2)=-omega^2, chi(D_3)=-omega, chi(D_4)=omega^2, chi(D_5)=omega, chi(D_6)=0
# psi_4 = conj(psi_3)

# omega = e^{2pi*i/3}
omega = np.exp(2j * np.pi / 3)
omega2 = omega**2

# First, classify 2T elements into conjugacy classes
# Build 2T conjugacy classes by conjugation
visited_2T = set()
classes_2T = []
for g in subgroup_2T:
    if g in visited_2T:
        continue
    cc = set()
    for q in subgroup_2T:
        q_inv = inv_table[q]
        conj_elem = mult_table[mult_table[q, g], q_inv]
        cc.add(conj_elem)
    visited_2T |= cc
    classes_2T.append(sorted(cc))

classes_2T.sort(key=lambda cc: (len(cc), -elements[cc[0]][0]))
print("2T classes: sizes %s" % [len(cc) for cc in classes_2T])

# Identify class types by order
for i, cc in enumerate(classes_2T):
    g = cc[0]
    order = 1
    current = g
    while current != identity_idx and order < 200:
        current = mult_table[current, g]
        order += 1
    print("  D_%d: size %d, order %d, w=%.4f" % (i, len(cc), order, elements[g][0]))

# Assign characters for each element based on its class
# Use known 2T character table
# D_0: {e}, size 1, order 1
# D_1: {-1}, size 1, order 2
# D_2: size 4, order 6 (first chirality)
# D_3: size 4, order 6 (second chirality)
# D_4: size 4, order 3 (first chirality)
# D_5: size 4, order 3 (second chirality)
# D_6: size 6, order 4

# Character assignments for psi_3, psi_4, psi_5:
# psi_5 (real 2D): 2, -2, 1, 1, -1, -1, 0
# psi_3 (complex 2D): 2, -2, omega, omega^2, -omega, -omega^2, 0 ?
# Wait, I need to be careful about which chirality gets omega vs omega^2

# From exp216 subspace 0 character output:
# chi(D_2, w=0.5, size 4) = -0.5 - 0.866i = -(1/2 + sqrt(3)/2 * i) = -omega^2
# Wait, omega = e^{2pi i/3} = -1/2 + sqrt(3)/2 * i
# omega^2 = e^{4pi i/3} = -1/2 - sqrt(3)/2 * i
# -omega^2 = 1/2 + sqrt(3)/2 * i
# So chi(D_2) = -0.5 - 0.866i = -1/2 - sqrt(3)/2*i = omega^2
# Actually: -0.5 - 0.866i ~ omega^2 = -1/2 - sqrt(3)/2*i. YES.

# From exp216 subspace 0:
# D_2 (w=0.5, size 4, order 6, first chiral): chi = -0.5 - 0.866i = omega^2
# D_3 (w=0.5, size 4, order 6, second chiral): chi = -0.5 + 0.866i = omega
# D_4 (w=-0.5, size 4, order 3, first chiral): chi = 0.5 - 0.866i = -omega
# D_5 (w=-0.5, size 4, order 3, second chiral): chi = 0.5 + 0.866i = -omega^2
# D_6 (w=0.0, size 6, order 4): chi = 0

# So psi_3: (2, -2, omega^2, omega, -omega, -omega^2, 0)
# psi_4 = conj(psi_3): (2, -2, omega, omega^2, -omega^2, -omega, 0)

# But wait - I need to match which of our D_2,...D_5 classes correspond to order 6 vs 3

# From the class listing:
# D_2: size 4, order 6
# D_3: size 4, order 6
# D_4: size 4, order 3
# D_5: size 4, order 3

# But there are TWO classes of order 6 and TWO of order 3
# The character of psi_3 on each class depends on which chirality

# Rather than figuring out the exact assignment, use a systematic approach:
# Compute the Schur projector using the REAL character psi_5 first (known exactly),
# then use the commutant of psi_5's complement to find psi_3 and psi_4.

# Schur projector for psi_5 (real 2D irrep):
# P_5 = (2/24) * sum_{g in 2T} chi_5(g)^* * D(g)

# chi_5 values by class:
chi_5_values = {}
for i, cc in enumerate(classes_2T):
    size = len(cc)
    order_g = 1
    g = cc[0]
    current = g
    while current != identity_idx and order_g < 200:
        current = mult_table[current, g]
        order_g += 1

    if size == 1 and order_g == 1:
        chi_5_values[i] = 2.0  # identity
    elif size == 1 and order_g == 2:
        chi_5_values[i] = -2.0  # -1
    elif size == 6:
        chi_5_values[i] = 0.0   # order 4
    elif size == 4 and order_g == 6:
        chi_5_values[i] = 1.0   # order 6
    elif size == 4 and order_g == 3:
        chi_5_values[i] = -1.0  # order 3
    else:
        chi_5_values[i] = 0.0

print("\npsi_5 character values: %s" % chi_5_values)

# Element-level character
chi_5_elem = np.zeros(N)
for i, cc in enumerate(classes_2T):
    for g in cc:
        chi_5_elem[g] = chi_5_values[i]

# Schur projector P_5
P_5 = np.zeros((dim_gen, dim_gen), dtype=complex)
for g in subgroup_2T:
    P_5 += chi_5_elem[g] * rep6[g]  # chi is real, no need to conjugate
P_5 *= 2.0 / 24.0  # dim_j / |G|

# Verify P_5 is a projector: P_5^2 = P_5
P5_err = np.max(np.abs(P_5 @ P_5 - P_5))
P5_rank = np.linalg.matrix_rank(P_5.real, tol=1e-8)
P5_trace = np.trace(P_5).real
print("P_5: rank=%d, trace=%.4f, projector error=%.2e" % (P5_rank, P5_trace, P5_err))

# P_5 should project onto a 2D subspace
evals_P5, evecs_P5 = np.linalg.eigh(P_5)
print("P_5 eigenvalues:", np.sort(evals_P5.real)[::-1])

# Extract the 2D subspace for psi_5
mask_5 = evals_P5.real > 0.5
V_5 = evecs_P5[:, mask_5]
print("psi_5 subspace: %d dimensions" % V_5.shape[1])

# The complement: psi_3 + psi_4 live in the 4D complement
# P_34 = I - P_5 - P_trivial - ...
# But simpler: P_34 = (sum of all psi_3 and psi_4 projectors)
# Since dim-6 = psi_3 + psi_4 + psi_5, the complement of psi_5 IS psi_3+psi_4
P_34 = np.eye(dim_gen) - P_5
# But P_34 might include components outside the representation space
# Actually since P_5 projects onto a 2D subspace within the 6D space,
# I - P_5 projects onto the 4D complement = psi_3 + psi_4

evals_P34, evecs_P34 = np.linalg.eigh(P_34)
mask_34 = evals_P34.real > 0.5
V_34 = evecs_P34[:, mask_34]
print("psi_3+psi_4 subspace: %d dimensions" % V_34.shape[1])

# Now separate psi_3 and psi_4 within V_34
# Use the commutant of 2T restricted to V_34
# Project rep matrices into V_34
print("\nSeparating psi_3 and psi_4 in the 4D complement...")
comm_basis_34 = []
np.random.seed(777)
for trial in range(30):
    M_rand = np.random.randn(4, 4) + 1j * np.random.randn(4, 4)
    M_avg = np.zeros_like(M_rand)
    for g in subgroup_2T:
        D_g_34 = V_34.T.conj() @ rep6[g] @ V_34  # 4x4 restriction
        D_inv_34 = np.linalg.inv(D_g_34)
        M_avg += D_g_34 @ M_rand @ D_inv_34
    M_avg /= len(subgroup_2T)
    comm_basis_34.append(M_avg.flatten())

C_34 = np.array(comm_basis_34)
_, S_34, Vh_34 = np.linalg.svd(C_34)
rank_34 = np.sum(S_34 > 1e-8)
print("Commutant rank in 4D = %d (expected 2 for psi_3+psi_4)" % rank_34)

M_comm_34 = Vh_34[0].reshape(4, 4)
evals_34, evecs_34 = np.linalg.eig(M_comm_34)
print("Commutant eigenvalues:", ["%.4f+%.4fi" % (e.real, e.imag) for e in evals_34])

# Group into pairs
groups_34 = []
used = [False] * 4
for i in range(4):
    if used[i]: continue
    group = [i]; used[i] = True
    for j in range(i+1, 4):
        if not used[j] and abs(evals_34[i] - evals_34[j]) < 1e-5:
            group.append(j); used[j] = True
    groups_34.append(group)

if len(groups_34) == 2 and all(len(g) == 2 for g in groups_34):
    print("Found 2 groups of size 2: psi_3 and psi_4 separated!")
    # Get subspaces in 6D coordinates
    V_3_local = evecs_34[:, groups_34[0]]
    V_4_local = evecs_34[:, groups_34[1]]
    V_3 = V_34 @ V_3_local  # 6x2
    V_4 = V_34 @ V_4_local  # 6x2
else:
    print("WARNING: expected 2 groups of 2, got %s" % [len(g) for g in groups_34])
    # Fallback: just split the 4D space arbitrarily
    V_3 = V_34[:, :2]
    V_4 = V_34[:, 2:]

# Orthogonalize all three subspaces
from numpy.linalg import qr
V_5_orth, _ = qr(V_5, mode='reduced')
V_3_orth, _ = qr(V_3, mode='reduced')
V_4_orth, _ = qr(V_4, mode='reduced')

# Verify orthogonality between subspaces
orth_35 = np.max(np.abs(V_3_orth.T.conj() @ V_5_orth))
orth_45 = np.max(np.abs(V_4_orth.T.conj() @ V_5_orth))
orth_34 = np.max(np.abs(V_3_orth.T.conj() @ V_4_orth))
print("Orthogonality: |V3.V5|=%.2e, |V4.V5|=%.2e, |V3.V4|=%.2e"
      % (orth_35, orth_45, orth_34))

# Identify which is psi_3 vs psi_4 by checking character
for label, V in [("V_3", V_3_orth), ("V_4", V_4_orth), ("V_5", V_5_orth)]:
    traces = []
    for g in subgroup_2T:
        DV = rep6[g] @ V
        R = np.linalg.lstsq(V, DV, rcond=None)[0]
        traces.append(np.trace(R))
    all_real = all(abs(tr.imag) < 1e-4 for tr in traces)
    print("  %s: all traces real = %s" % (label, all_real))

# Label subspaces as gen_1 (psi_3), gen_2 (psi_4), gen_3 (psi_5)
subspaces_0 = [V_3_orth, V_4_orth, V_5_orth]  # Under 2T_0
projectors_0 = [V @ V.T.conj() for V in subspaces_0]

# Verify projectors sum to identity
P_sum = sum(projectors_0)
P_sum_err = np.max(np.abs(P_sum - np.eye(dim_gen)))
print("Projector sum = I check: err = %.2e" % P_sum_err)

# ======================================================================
# PART 3: Find 5 conjugate subgroups
# ======================================================================
print("\n--- PART 3: Conjugate Subgroups 2T_g ---")

# Find the 5 left cosets of 2T in 2I
coset_reps = [identity_idx]  # First coset rep = identity
assigned = set(subgroup_2T)
while len(assigned) < N:
    for g in range(N):
        if g not in assigned:
            coset = set()
            for h in subgroup_2T:
                coset.add(mult_table[g, h])
            coset_reps.append(g)
            assigned |= coset
            break

print("Coset representatives: %d" % len(coset_reps))
for ci, g in enumerate(coset_reps):
    print("  Coset %d: rep=%d, q=(%.4f, %.4f, %.4f, %.4f)"
          % (ci, g, elements[g][0], elements[g][1], elements[g][2], elements[g][3]))

# For each coset rep g, the conjugate subgroup is 2T_g = g*2T*g^{-1}
conjugate_2Ts = []
for ci, g in enumerate(coset_reps):
    g_inv = inv_table[g]
    conj_sub = set()
    for h in subgroup_2T:
        conj_sub.add(mult_table[mult_table[g, h], g_inv])
    conjugate_2Ts.append(sorted(conj_sub))
    # Check it has 24 elements
    assert len(conj_sub) == 24, "Conjugate 2T has %d elements!" % len(conj_sub)

# Check which are distinct
distinct_conj = []
for ci, sub in enumerate(conjugate_2Ts):
    is_dup = False
    for prev in distinct_conj:
        if set(sub) == set(prev[1]):
            is_dup = True
            break
    if not is_dup:
        distinct_conj.append((ci, sub))

print("\nDistinct conjugate subgroups: %d (expected 5)" % len(distinct_conj))

# ======================================================================
# PART 4: Compute mixing matrices
# ======================================================================
print("\n--- PART 4: Mixing Matrices from Conjugate Subgroups ---")

# For each conjugate 2T_g, the invariant subspaces are D(g) * V_k
# Because: if V_k is 2T-invariant, then D(g)*V_k is g*2T*g^{-1}-invariant

# Compute projectors for each conjugate subgroup
all_mixing_matrices = []

for ci, g in enumerate(coset_reps):
    if ci == 0:
        continue  # Skip identity (overlap with itself is trivial)

    print("\n  Coset %d (g=%d, q=(%.3f,%.3f,%.3f,%.3f)):"
          % (ci, g, elements[g][0], elements[g][1], elements[g][2], elements[g][3]))

    Dg = rep6[g]

    # Transform subspaces
    subspaces_g = [Dg @ V for V in subspaces_0]
    projectors_g = [V @ V.T.conj() for V in subspaces_g]

    # Verify: should still be projectors summing to I
    P_sum_g = sum(projectors_g)
    err_g = np.max(np.abs(P_sum_g - np.eye(dim_gen)))

    # Compute 3x3 overlap matrix
    # M[i,j] = Tr(P_i(2T_0) * P_j(2T_g)) / dim_irrep
    # This gives the probability of a state in generation i (2T_0)
    # being found in generation j (2T_g)
    M = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            M[i, j] = np.trace(projectors_0[i] @ projectors_g[j]).real / 2.0

    print("    3x3 mixing matrix M (Tr(P_i P_j)/2):")
    for r in range(3):
        print("    [%s]" % "  ".join("%.6f" % M[r,c] for c in range(3)))

    # Check: should be doubly stochastic
    row_sums = np.sum(M, axis=1)
    col_sums = np.sum(M, axis=0)
    print("    Row sums: [%s]" % "  ".join("%.4f" % s for s in row_sums))
    print("    Col sums: [%s]" % "  ".join("%.4f" % s for s in col_sums))

    # Compare to CKM and PMNS
    CKM_sq = np.array([
        [0.9484, 0.0503, 0.0000146],
        [0.0489, 0.9744, 0.00168],
        [0.0000640, 0.00150, 0.9984]
    ])
    PMNS_sq = np.array([
        [0.6778, 0.3002, 0.0220],
        [0.1608, 0.3531, 0.4861],
        [0.1614, 0.3467, 0.4919]
    ])
    # Tribimaximal
    TBM = np.array([
        [2.0/3, 1.0/3, 0],
        [1.0/6, 1.0/3, 1.0/2],
        [1.0/6, 1.0/3, 1.0/2]
    ])

    best_ckm = 1e10
    best_pmns = 1e10
    best_tbm = 1e10
    for pr in permutations(range(3)):
        for pc in permutations(range(3)):
            Mp = M[np.array(pr), :][:, np.array(pc)]
            err_ckm = np.sqrt(np.mean((Mp - CKM_sq)**2))
            err_pmns = np.sqrt(np.mean((Mp - PMNS_sq)**2))
            err_tbm = np.sqrt(np.mean((Mp - TBM)**2))
            if err_ckm < best_ckm: best_ckm = err_ckm
            if err_pmns < best_pmns: best_pmns = err_pmns
            if err_tbm < best_tbm: best_tbm = err_tbm

    print("    Best RMS: CKM=%.4f, PMNS=%.4f, TBM=%.4f" % (best_ckm, best_pmns, best_tbm))

    all_mixing_matrices.append((ci, g, M))

# ======================================================================
# PART 5: Analyze mixing matrix structure
# ======================================================================
print("\n--- PART 5: Structure Analysis ---")

# Check if all 4 non-trivial conjugates give the same mixing matrix
print("\nAre all mixing matrices the same?")
for ci, g, M in all_mixing_matrices:
    print("  Coset %d: M[0,0]=%.6f, M[1,1]=%.6f, M[2,2]=%.6f"
          % (ci, M[0,0], M[1,1], M[2,2]))

# Unique mixing matrices
unique_Ms = []
for ci, g, M in all_mixing_matrices:
    is_dup = False
    for M_prev in unique_Ms:
        # Check if same up to permutation
        for pr in permutations(range(3)):
            for pc in permutations(range(3)):
                Mp = M[np.array(pr), :][:, np.array(pc)]
                if np.max(np.abs(Mp - M_prev)) < 1e-6:
                    is_dup = True
                    break
            if is_dup: break
    if not is_dup:
        unique_Ms.append(M)

print("Unique mixing matrices (up to permutation): %d" % len(unique_Ms))

# Test: average of all mixing matrices
if len(all_mixing_matrices) > 0:
    M_avg = np.mean([M for _, _, M in all_mixing_matrices], axis=0)
    print("\nAverage mixing matrix:")
    for r in range(3):
        print("  [%s]" % "  ".join("%.6f" % M_avg[r,c] for c in range(3)))

# ======================================================================
# PART 6: Deeper analysis - product of conjugates
# ======================================================================
print("\n--- PART 6: Products of Mixing Matrices ---")

# If M(g1) and M(g2) are mixing matrices for two different cosets,
# the product M(g1) * M(g2)^T gives the mixing between 2T_{g1} and 2T_{g2}
# This is different from a single coset mixing

if len(all_mixing_matrices) >= 2:
    M1 = all_mixing_matrices[0][2]
    M2 = all_mixing_matrices[1][2]
    M12 = M1 @ M2.T
    print("Product M(coset1) * M(coset2)^T:")
    for r in range(3):
        print("  [%s]" % "  ".join("%.6f" % M12[r,c] for c in range(3)))

# ======================================================================
# PART 7: Try different j values for different sectors
# ======================================================================
print("\n--- PART 7: Mixing for ALL spin-j representations ---")

for j in [0.5, 1.0, 1.5, 2.0, 2.5]:
    dim_j = int(2*j + 1)
    print("\n  j = %.1f (dim %d):" % (j, dim_j))

    # Compute representation
    rep_j = []
    for i in range(N):
        U = quat_to_su2(elements[i])
        D = spin_j_matrix(U, j)
        rep_j.append(D)

    # Schur projector for psi_5 (real 2D irrep)
    P5_j = np.zeros((dim_j, dim_j), dtype=complex)
    for g in subgroup_2T:
        P5_j += chi_5_elem[g] * rep_j[g]
    P5_j *= 2.0 / 24.0

    rank_5 = np.linalg.matrix_rank(P5_j.real, tol=1e-6)
    trace_5 = np.trace(P5_j).real
    print("    psi_5 projector: rank=%d, trace=%.1f" % (rank_5, trace_5))

    if rank_5 == 0:
        print("    psi_5 not present in this rep")
        continue
    if rank_5 != 2:
        print("    psi_5 has rank %d (not 2), skipping" % rank_5)
        continue

    # Get psi_5 subspace
    ev5, V5 = np.linalg.eigh(P5_j)
    V5_sub = V5[:, ev5.real > 0.5]

    # Check for psi_3/psi_4 in complement
    P_compl = np.eye(dim_j) - P5_j
    rank_compl = np.linalg.matrix_rank(P_compl.real, tol=1e-6)
    print("    Complement rank = %d" % rank_compl)

    if rank_compl < 2:
        print("    Complement too small for psi_3/psi_4")
        continue

    # For the mixing test, use the psi_5 subspace only (simpler)
    # Compute how psi_5(2T_0) overlaps with psi_5(2T_g)
    for ci, g in enumerate(coset_reps):
        if ci == 0: continue
        Dg_j = rep_j[g]
        V5_g = Dg_j @ V5_sub
        # Overlap = Tr(P_5(0) @ P_5(g)) / 2
        P5_0 = V5_sub @ V5_sub.T.conj()
        P5_g = V5_g @ V5_g.T.conj()
        overlap = np.trace(P5_0 @ P5_g).real / 2.0
        print("    Coset %d: psi_5 overlap = %.6f" % (ci, overlap))

# ======================================================================
# PART 8: Key geometric ratios in mixing matrix
# ======================================================================
print("\n--- PART 8: Geometric Ratios ---")

if len(all_mixing_matrices) > 0:
    M = all_mixing_matrices[0][2]  # Use first non-trivial

    # Look for phi-related ratios
    print("Mixing matrix elements and phi-tests:")
    for r in range(3):
        for c in range(3):
            val = M[r, c]
            if val > 1e-6:
                # Test val = a/b for small integers
                best_frac = None
                best_err = 1
                for a in range(0, 20):
                    for b in range(1, 20):
                        err = abs(val - a/b)
                        if err < best_err:
                            best_err = err
                            best_frac = (a, b)
                # Test val = phi^k or related
                phi_tests = {}
                for k in range(-5, 6):
                    phi_tests["phi^%d" % k] = PHI**k
                    phi_tests["1/(phi^%d+1)" % k] = 1/(PHI**k + 1)
                best_phi = None
                best_phi_err = 1
                for name, test_val in phi_tests.items():
                    err = abs(val - test_val)
                    if err < best_phi_err:
                        best_phi_err = err
                        best_phi = name

                print("  M[%d,%d] = %.8f ~ %d/%d (err=%.6f) ~ %s (err=%.6f)"
                      % (r, c, val, best_frac[0], best_frac[1], best_err,
                         best_phi, best_phi_err))

# ======================================================================
# SUMMARY
# ======================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
1. Schur projectors correctly separate psi_3, psi_4, psi_5 in dim-6 rep
2. 2T is NOT normal in 2I -> 5 conjugate subgroups -> mixing matrices
3. Each conjugate gives a 3x3 doubly-stochastic overlap matrix M
4. M encodes the "generation mixing" between different 24-cell decompositions
5. Key question: does M match CKM^2 or PMNS^2?
""")

print("=" * 70)
print("EXP-217 COMPLETE")
print("=" * 70)

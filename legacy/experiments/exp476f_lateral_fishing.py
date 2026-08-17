"""
exp476f: Lateral Fishing on the Decomposition Graph
====================================================

Trying UNEXPECTED directions on the 40 E8/H4 decompositions:

A) Galois action: does phi->phi' permute the 40 decompositions?
B) Spectral determinant: det'(T), det'(E) -- framework numbers?
C) Power traces: Tr(T^k), Tr(E^k) -- combinatorial invariants
D) Heat kernel at special times (t=1/phi^2, t=1/a1, etc.)
E) The INTEGER structure: T has entries in {2, 2.4, 2.8, 3.2, 4}
   = (2/a1)*{5,6,7,8,10}. So (a1/2)*T is INTEGER. Determinants over Z?

Looking for SURPRISE framework numbers.

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: FISHING EXPEDITION
"""

import numpy as np
from itertools import permutations as iterperms
from math import cos, sin, pi, sqrt, factorial, log, exp
from collections import Counter

phi = (1 + sqrt(5)) / 2
phi_p = (1 - sqrt(5)) / 2
a1 = 5; b1 = 6; h = 30; N = 120; rank_E8 = 8
alpha = 1/137.035999084
alpha_s = 1/(2*phi**3)

# ============================================================
# BUILD 40 DECOMPOSITIONS (compact)
# ============================================================
print("Building 40 decompositions...")

E8_roots = []
for i in range(8):
    for j in range(i+1, 8):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = np.zeros(8); v[i] = si; v[j] = sj
                E8_roots.append(v)
for bits in range(256):
    signs = [(-1)**((bits >> k) & 1) for k in range(8)]
    if sum(1 for s in signs if s < 0) % 2 == 0:
        E8_roots.append(np.array(signs) / 2.0)
E8 = np.array(E8_roots)

simple = np.zeros((8, 8))
simple[0]=[1,-1,0,0,0,0,0,0]; simple[1]=[0,1,-1,0,0,0,0,0]
simple[2]=[0,0,1,-1,0,0,0,0]; simple[3]=[0,0,0,1,-1,0,0,0]
simple[4]=[0,0,0,0,1,-1,0,0]; simple[5]=[0,0,0,0,0,1,-1,0]
simple[6]=[0,0,0,0,0,1,1,0];  simple[7]=[-0.5]*8
refs = [np.eye(8)-2*np.outer(simple[i],simple[i])/np.dot(simple[i],simple[i]) for i in range(8)]

def make_cox(perm):
    C = np.eye(8)
    for i in perm: C = refs[i] @ C
    return C

def get_proj(C):
    if not np.allclose(np.linalg.matrix_power(C,30),np.eye(8),atol=1e-8): return None
    eigvals,eigvecs = np.linalg.eig(C)
    args = np.angle(eigvals)
    exps = [(round(args[i]*30/(2*pi))%30 or 30,i) for i in range(8)]
    h4_idx = [i for m,i in exps if m in {1,11,19,29}]
    if len(h4_idx)!=4: return None
    v_raw = eigvecs[:,h4_idx]
    basis = []; used = set()
    for k in range(4):
        if k in used: continue
        v = v_raw[:,k]
        if np.max(np.abs(np.imag(v)))>0.01:
            basis.extend([np.real(v),np.imag(v)])
            for k2 in range(k+1,4):
                if k2 not in used and np.allclose(v_raw[:,k2],np.conj(v)):
                    used.add(k2); break
        else: basis.append(np.real(v))
    if len(basis)!=4: return None
    B = np.array(basis).T; Q,_ = np.linalg.qr(B)
    return Q[:,:4]

# Also store the Coxeter elements themselves
projs = []; coxeters = []; perm_list = []
np.random.seed(42)
perms = list(iterperms(range(8))); np.random.shuffle(perms)
for perm in perms:
    C = make_cox(perm)
    P = get_proj(C)
    if P is None: continue
    if all(not np.allclose(np.linalg.svd(Po.T@P,compute_uv=False),np.ones(4),atol=1e-6) for Po in projs):
        projs.append(P)
        coxeters.append(C)
        perm_list.append(perm)
    if len(projs)>=40: break
assert len(projs)==40
print(f"Found 40/40")

# Build matrices
Proj8 = [P @ P.T for P in projs]
M_int = np.zeros((40,40), dtype=int)  # overlap/12 = m-level
for i in range(40):
    for j in range(i+1,40):
        svs = np.linalg.svd(projs[i].T@projs[j],compute_uv=False)
        m = round(h*np.sum(svs**2)) // 12
        M_int[i,j] = m; M_int[j,i] = m
    M_int[i,i] = 10

T = np.zeros((40,40))
for i in range(40):
    for j in range(40):
        T[i,j] = 2.0 * M_int[i,j] / a1

# Integer version: (a1/2)*T = M_int scaled
M_scaled = M_int  # entries in {5,6,7,8,10}

# ============================================================
# A) GALOIS ACTION ON THE 40 DECOMPOSITIONS
# ============================================================
print(f"\n{'='*70}")
print("A) GALOIS ACTION: phi -> phi'")
print("="*70)

# The Galois conjugation sends phi -> phi' = -1/phi
# On the Coxeter element level, this swaps exponents 1<->11 and 19<->29
# (since exp(2pi*i*1/30) is conjugated to exp(2pi*i*11/30) under sqrt(5)->-sqrt(5))
# The H4 eigenplane has exponents {1,11,19,29}.
# The Galois involution swaps 1<->11 and 19<->29 WITHIN the H4 set.
# So the H4 eigenplane is GALOIS-INVARIANT as a set of exponents.
# But the PROJECTOR might change because the eigenvectors change.

# Alternative approach: the Galois conjugation acts on R^8 by swapping
# the H4 and H4' eigenspaces (exponents {1,11,19,29} <-> {7,13,17,23}).
# This gives a COMPLEMENTARY projector: P' = I - P.

# Let's check: does P -> I-P map our 40 projectors to themselves?
# If P_i projects onto a 4-plane V_i, then I-P_i projects onto V_i^perp.
# The overlap: Tr((I-P_i)(I-P_j)) = Tr(I - P_i - P_j + P_i P_j)
# = 8 - 4 - 4 + Tr(P_i P_j) = Tr(P_i P_j)
# So T_complement = T! The Gram matrix is GALOIS-INVARIANT.

print(f"\n  Key insight: P -> P_complement = I - P gives")
print(f"  Tr(P_c_i @ P_c_j) = 8 - 4 - 4 + Tr(Pi@Pj) = Tr(Pi@Pj)")
print(f"  So the GRAM MATRIX T is Galois-invariant (algebraic identity).")

# But: does P -> I-P MAP the set of 40 projectors TO ITSELF?
# i.e., for each Pi, is I-Pi one of the 40 projectors?

print(f"\n  Does P -> I-P permute the 40 decompositions?")
galois_perm = []
for i in range(40):
    Pc = np.eye(8) - Proj8[i]  # complement projector
    found = -1
    for j in range(40):
        # Check if Pc == Proj8[j] (as subspaces)
        svs = np.linalg.svd(projs[j].T @ (np.eye(8) - Proj8[i]) @ projs[j], compute_uv=False)
        # If Pc projects onto the same subspace as Pj, then Pj^T @ Pc @ Pj should be identity-like
        # Actually: Tr(Pc @ Pj) = Tr((I-Pi)Pj) = 4 - Tr(PiPj) = 4 - 2*m_ij/5
        # For Pc = Pj: Tr(Pc@Pj) = 4. So 4 - 2*m/5 = 4 => m = 0. But m >= 5.
        # So Tr((I-Pi)Pj) = 4 - T[i,j]. For j=i: 4 - 4 = 0. For j!=i: 4 - {2,2.4,2.8,3.2} = {2,1.6,1.2,0.8}
        # None of these equal 4, so the complement is NEVER one of the 40!
        pass

    galois_perm.append(found)

print(f"  Tr((I-Pi)@Pj) = 4 - Tr(Pi@Pj) = 4 - T[i,j]")
print(f"  For Pj to equal I-Pi: need Tr((I-Pi)@Pj) = 4 (rank 4)")
print(f"  But 4 - T[i,j] = 4 - {{2, 2.4, 2.8, 3.2, 4}} = {{0, 0.8, 1.2, 1.6, 2}}")
print(f"  Maximum is 2 (at m=5), but need 4. IMPOSSIBLE.")
print(f"\n  *** The Galois complement I-P is NEVER one of the 40! ***")
print(f"  The 40 H4 projectors and 40 complement H4' projectors are DISJOINT sets.")
print(f"  This means the Galois sector is geometrically SEPARATE.")

# But: does the complement set form its OWN 40-element graph?
# And what's the CROSS overlap between the two sets?
print(f"\n  Cross-overlaps Tr(Pi @ (I-Pj)) = 4 - T[i,j]:")
cross_vals = set()
for i in range(40):
    for j in range(40):
        cross_vals.add(round(4 - T[i,j], 4))
print(f"  Values: {sorted(cross_vals)}")
print(f"  = {{0, 0.8, 1.2, 1.6, 2}} = (2/a1)*{{0, 2, 3, 4, 5}}")
print(f"  Galois map: m -> 10-m on overlap levels!")
print(f"  This is EXACTLY the CKM duality sigma(m) = 2*a1 - m!")

# ============================================================
# B) SPECTRAL DETERMINANTS
# ============================================================
print(f"\n{'='*70}")
print("B) SPECTRAL DETERMINANTS AND FRAMEWORK NUMBERS")
print("="*70)

T_eigvals = np.sort(np.linalg.eigvalsh(T))[::-1]
nonzero_T = [e for e in T_eigvals if abs(e) > 0.01]

log_det_T = sum(log(abs(e)) for e in nonzero_T)
print(f"\n  log det'(T) = {log_det_T:.6f}")
print(f"  log det'(T) / ln(1/alpha) = {log_det_T / log(1/alpha):.6f}")
print(f"  log det'(T) / ln(phi) = {log_det_T / log(phi):.6f}")

# Integer matrix: (a1/2)*T has integer entries. Its determinant is an integer!
T_int = (a1 * T / 2).astype(int)  # wait, T[i,j] = 2m/5, so 5T/2 = m
# Actually (a1/2)*T[i,j] = (5/2)*(2m/5) = m. So M_int = (a1/2)*T on diagonal is 10.
# M_int is the overlap/12 matrix with integer entries {5,6,7,8,10}

print(f"\n  Integer matrix M = (a1/2)*T has entries in {{5,6,7,8}} off-diag, 10 on diag")
M_eigvals = np.sort(np.linalg.eigvalsh(M_int.astype(float)))[::-1]
nonzero_M = [e for e in M_eigvals if abs(e) > 0.01]

log_det_M = sum(log(abs(e)) for e in nonzero_M)
det_M = np.prod([e for e in nonzero_M])
print(f"  log det'(M_int) = {log_det_M:.6f}")
print(f"  det'(M_int) = {det_M:.4e}")
print(f"  log det'(M_int) / ln(1/alpha) = {log_det_M / log(1/alpha):.6f}")

# E operator: integer matrix
E_int = np.zeros((40,40), dtype=int)
for i in range(40):
    for j in range(40):
        if i != j:
            E_int[i,j] = 10 - M_int[i,j]

E_eigvals = np.sort(np.linalg.eigvalsh(E_int.astype(float)))[::-1]
nonzero_E = [e for e in E_eigvals if abs(e) > 0.01]
log_det_E = sum(log(abs(e)) for e in nonzero_E)

print(f"\n  log det'(E) = {log_det_E:.6f}")
print(f"  log det'(E) / ln(1/alpha) = {log_det_E / log(1/alpha):.6f}")
print(f"  det'(E) = {np.prod([abs(e) for e in nonzero_E]):.4e}")

# ============================================================
# C) POWER TRACES
# ============================================================
print(f"\n{'='*70}")
print("C) POWER TRACES Tr(M^k)")
print("="*70)

M_float = M_int.astype(float)
for k in range(1, 9):
    Mk = np.linalg.matrix_power(M_float, k)
    tr = np.trace(Mk)
    print(f"  Tr(M^{k}) = {tr:.1f}", end="")
    # Check framework numbers
    for name, val in [('N', 120), ('h', 30), ('a1^2', 25), ('N^2', 14400),
                       ('40', 40), ('phi', phi), ('rank*a1', 40)]:
        if abs(val) > 0:
            r = tr / val
            if abs(r - round(r)) < 0.1 and 0 < abs(round(r)) < 10000:
                print(f"  = {int(round(r))}*{name}", end="")
    print()

# Tr(M) = sum of diagonal = 40*10 = 400
# Tr(M^2) = sum_ij M_ij^2

# ============================================================
# D) HEAT KERNEL AT SPECIAL TIMES
# ============================================================
print(f"\n{'='*70}")
print("D) HEAT KERNEL ON DECOMPOSITION GRAPH")
print("="*70)

# Use M_int as the "Hamiltonian" (it's positive, symmetric)
# K(t) = Tr(exp(-t * M_int))
# But M has large eigenvalues, so exp(-t*M) decays quickly.

# Better: use the graph Laplacian of one of the color classes
# Or: use L = 10*I - M_offdiag (a shifted Laplacian)
# L[i,j] = 10*delta[i,j] - m_ij*(1-delta[i,j]) -- this is the "exchanged" operator E with diagonal

# Actually use E (which has 0 diagonal and is the CKM-weighted adjacency)
# L_CKM = diag(row sums of E) - E (Laplacian)

row_sums_E = np.sum(E_int, axis=1)
L_CKM = np.diag(row_sums_E) - E_int.astype(float)
L_eigvals = np.sort(np.linalg.eigvalsh(L_CKM))

print(f"\n  CKM Laplacian L = D_E - E:")
print(f"  Smallest eigenvalues: {[f'{e:.4f}' for e in L_eigvals[:5]]}")
print(f"  Largest: {L_eigvals[-1]:.4f}")

# Spectral gap
lambda_1_CKM = L_eigvals[1] if abs(L_eigvals[0]) < 0.01 else L_eigvals[0]
print(f"  Spectral gap lambda_1 = {lambda_1_CKM:.6f}")

# Heat kernel at special times
special_times = {
    '1/phi^2': 1/phi**2,
    '1/a1': 1/a1,
    '1/b1': 1/b1,
    '1/h': 1/h,
    'alpha': alpha,
    'alpha_s': alpha_s,
    '1/N': 1/N,
    '1': 1.0,
    'ln(phi)': log(phi),
}

print(f"\n  Heat kernel K(t) = Tr(exp(-t*L_CKM)):")
for name, t in special_times.items():
    K = sum(exp(-t * e) for e in L_eigvals)
    print(f"    K({name:>8}) = {K:12.6f}", end="")
    # Check ratios
    for ref_name, ref_val in [('N', 120), ('40', 40), ('26', 26), ('14', 14),
                                ('a1^2+1', 26), ('2*n_23', 14), ('b1', 6)]:
        r = K / ref_val
        if abs(r - round(r, 1)) < 0.05:
            print(f"  ~ {r:.3f}*{ref_name}", end="")
    print()

# ============================================================
# E) THE Z[phi] STRUCTURE
# ============================================================
print(f"\n{'='*70}")
print("E) GALOIS NORM AND Z[phi] STRUCTURE")
print("="*70)

# T[i,j] = 2m/a1 where m in {5,6,7,8} for off-diagonal
# The values are: 2, 12/5, 14/5, 16/5, 4
# In Z[phi]: 2 = 2, 12/5 = ?, 14/5 = ?, 16/5 = ?, 4 = 4
# These are NOT in Z[phi] directly.

# But M_int has entries in Z.
# (a1/2)*T = M_int in Z.

# Galois norm of M_int entries:
# N(5) = 5, N(6) = 6, N(7) = 7, N(8) = 8, N(10) = 10
# These are just the numbers themselves (all rational).

# More interesting: the EIGENVALUES of T
print(f"\n  T eigenvalues and their Galois structure:")
for e in sorted(nonzero_T, reverse=True):
    # Check if eigenvalue is in Z[phi]
    # i.e., e = a + b*phi for some rationals a, b
    # phi = (1+sqrt(5))/2, so b*phi = b*(1+sqrt(5))/2
    # e = a + b/2 + b*sqrt(5)/2
    # So: rational part = a + b/2, irrational part = b*sqrt(5)/2
    # If e is rational: b = 0

    is_rational = abs(e - round(e, 6)) < 1e-4
    if is_rational:
        print(f"    {e:10.6f} -- rational ({round(e)})")
    else:
        # Try e = a + b*sqrt(5)
        # Check: e = c*phi + d for integers c, d
        # c*phi + d = c*(1+sqrt(5))/2 + d = (c/2 + d) + c*sqrt(5)/2
        for c in range(-20, 21):
            for d in range(-20, 21):
                if abs(e - (c*phi + d)) < 0.001:
                    norm = (c*phi + d) * (c*phi_p + d)
                    print(f"    {e:10.6f} = {c}*phi + {d}, N = {norm:.4f}")
                    break
            else:
                continue
            break
        else:
            print(f"    {e:10.6f} -- not in Z[phi] (small coefficients)")

# ============================================================
# F) THE CROSS-GRAM MATRIX (H4 x H4')
# ============================================================
print(f"\n{'='*70}")
print("F) H4 x H4' CROSS STRUCTURE")
print("="*70)

# From part A: Tr(Pi @ (I-Pj)) = 4 - T[i,j]
# This is the "cross-Gram matrix" between physical and Galois sectors

T_cross = 4 - T  # Cross overlaps
T_cross_eigvals = np.sort(np.linalg.eigvalsh(T_cross))[::-1]

# T_cross = 4*J - T. Eigenvalues: 4*40 - lambda(T) = 160 - lambda(T)... no.
# T_cross = 4*J40 - T where J40 is the 40x40 all-ones matrix
# J40 has eigenvalues 40 (x1) and 0 (x39)
# So this is not straightforward

print(f"  Cross Gram matrix: T_cross = 4*J - T")
print(f"  T_cross eigenvalues (top 5): {[f'{e:.4f}' for e in T_cross_eigvals[:5]]}")
print(f"  T_cross eigenvalues (bottom 5): {[f'{e:.4f}' for e in T_cross_eigvals[-5:]]}")

# T_cross diagonal: 4 - 4 = 0. Off-diagonal: 4 - {2, 2.4, 2.8, 3.2} = {2, 1.6, 1.2, 0.8}
# = (2/a1)*{5, 4, 3, 2} = (2/a1)*{a1, n_ut, n_12, n_db}
# So T_cross off-diagonal = (2/a1)*(2*a1 - m) = EXCHANGED multipliers scaled!

# T_cross = (2/a1) * E + 0*I = (2/a1)*E
# CHECK: E[i,j] = 10 - m for i!=j, 0 for i=i. T_cross[i,j] = 4 - 2m/5 = (10-m)*2/5 for i!=j, 0 for i=i.
# So T_cross = (2/a1)*E. YES!

print(f"\n  KEY: T_cross = (2/a1) * E")
print(f"  The cross-Gram matrix IS the CKM-weighted operator (up to scale)!")
print(f"  ker(T_cross) = ker(E) = 26-dim = a1^2+1")
print(f"\n  PHYSICAL MEANING:")
print(f"  The 26 'CKM-harmonic' modes are EXACTLY the modes where")
print(f"  the H4-H4' cross-overlap vanishes.")
print(f"  i.e., ker(E) = {{functions on decompositions with zero Galois cross-talk}}")

# ============================================================
# G) CHARACTERISTIC POLYNOMIAL OF M_int
# ============================================================
print(f"\n{'='*70}")
print("G) CHARACTERISTIC POLYNOMIAL")
print("="*70)

# M_int is 40x40 integer matrix. Its char poly has integer coefficients.
# det(M_int - lambda*I) = 0

# We can compute coefficients from eigenvalues
M_eigvals_all = np.sort(np.linalg.eigvalsh(M_int.astype(float)))[::-1]

# Elementary symmetric polynomials from eigenvalues
# e_1 = sum lambda_i = Tr(M) = 40*10 = 400
# e_2 = sum_{i<j} lambda_i*lambda_j = (Tr(M)^2 - Tr(M^2))/2
trM = np.trace(M_float)
trM2 = np.trace(M_float @ M_float)
trM3 = np.trace(np.linalg.matrix_power(M_float, 3))

print(f"  Tr(M) = {trM:.0f} = 40*10")
print(f"  Tr(M^2) = {trM2:.0f}")
print(f"  Tr(M^3) = {trM3:.0f}")

# Tr(M^2) = sum M[i,j]^2 = 40*100 (diagonal) + sum_{i!=j} m_ij^2
off_diag_sq = trM2 - 40*100
print(f"  sum_{i!=j} m_ij^2 = {off_diag_sq:.0f}")
print(f"    = 2*(176*25 + 296*36 + 184*49 + 124*64)")
expected = 2*(176*25 + 296*36 + 184*49 + 124*64)
print(f"    = {expected} (check: {abs(off_diag_sq - expected) < 1})")

# Average m^2
avg_m2 = off_diag_sq / (40*39)
print(f"  <m^2> = {avg_m2:.6f}")
print(f"  <m>^2 = {(sum(m*p for m,p in [(5,176),(6,296),(7,184),(8,124)])/780)**2:.6f}")

# ============================================================
# H) SURPRISE HUNT: RATIOS OF SPECTRAL DATA
# ============================================================
print(f"\n{'='*70}")
print("H) SURPRISE HUNT: RATIOS")
print("="*70)

fw = {'a1':5,'b1':6,'phi':phi,'phi^2':phi**2,'phi^3':phi**3,
      'N':120,'h':30,'rank':8,'n_12':3,'n_23':7,'n_13':12,
      'alpha':alpha,'alpha_s':alpha_s,'1/alpha':1/alpha,
      'sin2tW':6/26,'26':26,'14':14,'40':40,'pi':pi,
      'z_Planck':4*phi**2, 'z_CC':57-alpha_s}

# Compute various "natural" quantities
quantities = {}
quantities['det_T'] = np.prod(nonzero_T)
quantities['det_E'] = np.prod([abs(e) for e in nonzero_E])
quantities['log_det_T'] = log_det_T
quantities['log_det_E'] = log_det_E
quantities['Tr(E^2)'] = np.trace(E_int.astype(float) @ E_int.astype(float))
quantities['Tr(M^2)_offdiag'] = off_diag_sq
quantities['sum_nonzero_E_pos'] = sum(e for e in nonzero_E if e > 0)
quantities['sum_|nonzero_E|'] = sum(abs(e) for e in nonzero_E)
quantities['lambda_max_T'] = max(nonzero_T)
quantities['lambda_max_E'] = max(nonzero_E)
quantities['spectral_gap_L'] = lambda_1_CKM

print(f"\nQuantity ratios against framework numbers:")
for q_name, q_val in quantities.items():
    if abs(q_val) < 1e-10: continue
    hits = []
    for f_name, f_val in fw.items():
        if abs(f_val) < 1e-10: continue
        r = q_val / f_val
        if abs(r - round(r)) < 0.02 and 0 < abs(round(r)) <= 1000:
            hits.append(f"{int(round(r))}*{f_name}")
        elif abs(1/r - round(1/r)) < 0.02 and 0 < abs(round(1/r)) <= 100:
            hits.append(f"{f_name}/{int(round(1/r))}")
    if hits:
        print(f"  {q_name:>25} = {q_val:15.4f} : {', '.join(hits[:5])}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY: WHAT'S NEW AND SURPRISING?")
print("="*70)

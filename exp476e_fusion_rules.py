"""
exp476e: Selection Rules as Fusion Rules?
==========================================

From exp476b, the products A_m @ A_n have EXACT ZERO BLOCKS:
  A_8^2: zero on colors {5, 8}
  A_7^2: zero on color {8}
  A_7 @ A_8: zero on color {7}

These look like FUSION RULES. Our TQFT has SU(2)_3 (k+2=a1=5)
with 4 anyons {0, 1/2, 1, 3/2} and Fibonacci sub-theory.

QUESTION: Do the selection rules from E8 decomposition graph
MATCH the SU(2)_3 fusion rules under some mapping m <-> j?

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: INVESTIGATION
"""

import numpy as np
from itertools import permutations as iterperms
from math import cos, sin, pi, sqrt, factorial
from collections import Counter

phi = (1 + sqrt(5)) / 2
a1 = 5; b1 = 6; h = 30; N = 120; rank_E8 = 8

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

projs = []; np.random.seed(42)
perms = list(iterperms(range(8))); np.random.shuffle(perms)
for perm in perms:
    P = get_proj(make_cox(perm))
    if P is None: continue
    if all(not np.allclose(np.linalg.svd(Po.T@P,compute_uv=False),np.ones(4),atol=1e-6) for Po in projs):
        projs.append(P)
    if len(projs)>=40: break
assert len(projs)==40
print(f"Found 40/40")

# Build adjacency matrices
M_ov = np.zeros((40,40),dtype=int)
for i in range(40):
    for j in range(i+1,40):
        svs = np.linalg.svd(projs[i].T@projs[j],compute_uv=False)
        ov = round(h*np.sum(svs**2))
        M_ov[i,j]=ov; M_ov[j,i]=ov
    M_ov[i,i]=N

A = {}
for m in [5,6,7,8]:
    A[m] = ((M_ov==12*m) & ~np.eye(40,dtype=bool)).astype(int)

# ============================================================
# PART 1: COMPLETE FUSION TABLE
# ============================================================
print(f"\n{'='*70}")
print("PART 1: COMPLETE SELECTION RULE TABLE")
print("="*70)

colors = [0, 5, 6, 7, 8]  # 0 = diagonal (identity)
color_labels = {0: 'I', 5: 'A5', 6: 'A6', 7: 'A7', 8: 'A8'}

# For each product A_m @ A_n, determine which colors have ALL ZERO entries
fusion = {}
print(f"\n{'Product':>10} | {'I':>4} {'A5':>4} {'A6':>4} {'A7':>4} {'A8':>4} | Allowed channels")
print("-" * 65)

for m in [5,6,7,8]:
    for n in range(m, 9):
        if n > 8: break
        prod = A[m] @ A[n]

        channels = {}
        for c in colors:
            if c == 0:
                vals = [prod[i,i] for i in range(40)]
            else:
                vals = [prod[i,j] for i in range(40) for j in range(40) if A[c][i,j]==1]

            all_zero = all(v == 0 for v in vals)
            any_nonzero = any(v != 0 for v in vals)
            channels[c] = 0 if all_zero else 1

        allowed = [color_labels[c] for c in colors if channels[c] == 1]
        status = " ".join(f"{channels[c]:4d}" for c in colors)
        print(f"  A{m}*A{n:1d} | {status} | {', '.join(allowed)}")
        fusion[(m,n)] = channels
        if m != n:
            fusion[(n,m)] = channels

# ============================================================
# PART 2: SU(2)_3 FUSION RULES
# ============================================================
print(f"\n{'='*70}")
print("PART 2: SU(2)_3 FUSION RULES (k=3, k+2=a1=5)")
print("="*70)

# Representations: j = 0, 1, 2, 3 (using integer labeling 2j)
# Fusion: j1 x j2 = sum_{j=|j1-j2|, step 2}^{min(j1+j2, 2k-j1-j2)} j
k = 3

def su2k_fusion(j1, j2, k):
    """SU(2)_k fusion: returns list of j in product j1 x j2."""
    result = []
    for j in range(abs(j1-j2), min(j1+j2, 2*k-j1-j2)+1, 2):
        result.append(j)
    return result

print(f"\nSU(2)_{k} fusion table (integer labels 0,1,2,3 = spins 0,1/2,1,3/2):")
print(f"{'j1 x j2':>8} | Result")
print("-" * 35)
su2_fusion = {}
for j1 in range(k+1):
    for j2 in range(j1, k+1):
        result = su2k_fusion(j1, j2, k)
        print(f"  {j1} x {j2}   | {result}")
        su2_fusion[(j1,j2)] = set(result)
        su2_fusion[(j2,j1)] = set(result)

# ============================================================
# PART 3: TRY ALL 24 MAPPINGS m <-> j
# ============================================================
print(f"\n{'='*70}")
print("PART 3: COMPARING ALL 24 MAPPINGS")
print("="*70)

m_values = [5, 6, 7, 8]
j_values = [0, 1, 2, 3]

best_score = -1
best_mapping = None

print(f"\n{'Mapping':>30} | {'Match':>5} {'Mismatch':>8} {'Score':>6}")
print("-" * 60)

for perm in iterperms(j_values):
    mapping = dict(zip(m_values, perm))  # m -> j
    inv_mapping = {j: m for m, j in mapping.items()}  # j -> m

    matches = 0
    mismatches = 0

    for m1 in m_values:
        for m2 in range(m1, 9):
            if m2 > 8: break
            j1, j2 = mapping[m1], mapping[m2]
            su2_result = su2_fusion[(j1, j2)]

            # What the graph says: which colors are allowed?
            graph_allowed = set()
            for c in [5,6,7,8]:
                if fusion[(m1,m2)][c] == 1:
                    graph_allowed.add(mapping[c])

            # Diagonal (identity) corresponds to j=0 in fusion
            if fusion[(m1,m2)][0] == 1:
                graph_allowed.add(0)  # identity channel

            # Compare: does graph_allowed match su2_result?
            # Score by: #(correct zeros) + #(correct nonzeros)
            for j in j_values:
                if j in su2_result and j in graph_allowed:
                    matches += 1
                elif j not in su2_result and j not in graph_allowed:
                    matches += 1
                else:
                    mismatches += 1

    score = matches / (matches + mismatches) * 100
    map_str = f"5->{perm[0]}, 6->{perm[1]}, 7->{perm[2]}, 8->{perm[3]}"
    print(f"  {map_str:>28} | {matches:5d} {mismatches:8d} {score:5.1f}%")

    if score > best_score:
        best_score = score
        best_mapping = mapping

print(f"\n  Best mapping: {best_mapping} with score {best_score:.1f}%")

# ============================================================
# PART 4: DETAILED COMPARISON FOR BEST MAPPING
# ============================================================
print(f"\n{'='*70}")
print(f"PART 4: DETAILED COMPARISON (best mapping)")
print("="*70)

inv_best = {j: m for m, j in best_mapping.items()}
spin_names = {0: '0', 1: '1/2', 2: '1', 3: '3/2'}

print(f"\nMapping: m -> spin")
for m in m_values:
    j = best_mapping[m]
    print(f"  overlap-{m} (multiplier {2*a1-m}) <-> spin {spin_names[j]} (j={j})")

print(f"\nDetailed comparison:")
print(f"{'Product':>10} | {'SU(2)_3 prediction':>25} | {'Graph observation':>25} | {'Match?':>6}")
print("-" * 80)

total_match = 0
total_mismatch = 0

for m1 in m_values:
    for m2 in range(m1, 9):
        if m2 > 8: break
        j1, j2 = best_mapping[m1], best_mapping[m2]
        su2_pred = su2_fusion[(j1, j2)]
        su2_str = "{" + ",".join(spin_names[j] for j in sorted(su2_pred)) + "}"

        graph_obs = set()
        for c in [5,6,7,8]:
            if fusion[(m1,m2)][c] == 1:
                graph_obs.add(best_mapping[c])
        if fusion[(m1,m2)][0] == 1 and m1 == m2:
            graph_obs.add(0)  # diagonal = identity channel (only for m x m)
        graph_str = "{" + ",".join(spin_names[j] for j in sorted(graph_obs)) + "}"

        match = (su2_pred == graph_obs)
        if match:
            total_match += 1
            status = "YES"
        else:
            total_mismatch += 1
            # Show what's different
            extra = graph_obs - su2_pred
            missing = su2_pred - graph_obs
            status = f"NO (extra:{extra}, miss:{missing})"

        print(f"  A{m1}*A{m2} | {su2_str:>25} | {graph_str:>25} | {status}")

print(f"\n  Exact matches: {total_match}/{total_match+total_mismatch}")

# ============================================================
# PART 5: FIBONACCI SUB-THEORY
# ============================================================
print(f"\n{'='*70}")
print("PART 5: FIBONACCI FUSION (simpler test)")
print("="*70)

# Fibonacci: two objects, 1 and tau, with tau x tau = 1 + tau
# In SU(2)_3: tau = spin 1 (j=2 in integer labeling)
# Or: identify "even" (j=0,2) and "odd" (j=1,3) sectors

# Can we find a bipartition of {5,6,7,8} into two sets
# such that the "even x even" and "odd x odd" give "even",
# and "even x odd" gives "odd"?

print(f"\nTesting Z/2 grading on overlap colors:")
for partition in [({5,7},{6,8}), ({5,8},{6,7}), ({5,6},{7,8})]:
    even, odd = partition
    print(f"\n  Even={even}, Odd={odd}")
    consistent = True

    for m1 in m_values:
        for m2 in range(m1, 9):
            if m2 > 8: break
            # Product parity
            p1 = 'E' if m1 in even else 'O'
            p2 = 'E' if m2 in even else 'O'
            expected_parity = 'E' if p1 == p2 else 'O'

            # Actual: which parities appear in the product?
            actual_parities = set()
            for c in [5,6,7,8]:
                if fusion[(m1,m2)][c] == 1:
                    actual_parities.add('E' if c in even else 'O')

            if expected_parity not in actual_parities:
                consistent = False
            if len(actual_parities) == 1 and expected_parity in actual_parities:
                pass  # perfect
            elif len(actual_parities) > 1:
                pass  # not forbidden, just not pure

    print(f"  Z/2 consistent: {consistent}")

# ============================================================
# PART 6: THE "TRUNCATED" COMPARISON
# ============================================================
print(f"\n{'='*70}")
print("PART 6: ZERO-PATTERN ONLY (ignoring extra channels)")
print("="*70)

# A weaker test: do the FORBIDDEN channels match?
# i.e., does "graph says 0" imply "SU(2)_3 says 0"?

print(f"\nFor best mapping {best_mapping}:")
print(f"Checking: every graph zero is a SU(2)_3 zero?")

all_zeros_match = True
for m1 in m_values:
    for m2 in range(m1, 9):
        if m2 > 8: break
        j1, j2 = best_mapping[m1], best_mapping[m2]
        su2_pred = su2_fusion[(j1, j2)]

        for c in [5,6,7,8]:
            jc = best_mapping[c]
            graph_zero = (fusion[(m1,m2)][c] == 0)
            su2_zero = (jc not in su2_pred)

            if graph_zero and not su2_zero:
                print(f"  MISMATCH: A{m1}*A{m2} on color {c}: graph=0 but SU(2)_3 allows it")
                all_zeros_match = False
            elif not graph_zero and su2_zero:
                print(f"  EXTRA: A{m1}*A{m2} on color {c}: graph!=0 but SU(2)_3 forbids it")

if all_zeros_match:
    print(f"  ALL graph zeros are SU(2)_3 zeros: the graph is a REFINEMENT of fusion rules!")

# ============================================================
# PART 7: ALTERNATIVE -- VERLINDE RING STRUCTURE
# ============================================================
print(f"\n{'='*70}")
print("PART 7: VERLINDE RING MULTIPLICATION TABLE")
print("="*70)

# The Verlinde ring for SU(2)_3 has basis {e_0, e_1, e_2, e_3}
# with multiplication N_{ij}^k = fusion coefficients
# The S-matrix diagonalizes the fusion rules

print(f"\nSU(2)_3 fusion coefficients N_{{ij}}^k:")
for j1 in range(4):
    for j2 in range(j1, 4):
        result = su2k_fusion(j1, j2, k)
        n_vec = [0]*4
        for j in result:
            n_vec[j] = 1
        print(f"  N_{{{j1},{j2}}} = {n_vec}")

# The TQFT S-matrix for SU(2)_3
# S_{jj'} = sqrt(2/(k+2)) * sin(pi*(2j+1)*(2j'+1)/(k+2))
print(f"\nS-matrix (SU(2)_3, k+2={a1}):")
S_mat = np.zeros((4, 4))
for j1 in range(4):
    for j2 in range(4):
        S_mat[j1, j2] = sqrt(2/a1) * sin(pi*(j1+1)*(j2+1)/a1)
    s_str = " ".join(f"{S_mat[j1,j2]:+.4f}" for j2 in range(4))
    print(f"  S[{j1},:] = [{s_str}]")

# Verlinde formula: N_{ij}^k = sum_l S_{il} S_{jl} S_{kl}^* / S_{0l}
print(f"\nVerlinde formula verification:")
for j1 in range(4):
    for j2 in range(j1, 4):
        result_verlinde = []
        for j3 in range(4):
            N = sum(S_mat[j1,l] * S_mat[j2,l] * S_mat[j3,l] / S_mat[0,l] for l in range(4))
            N = round(N.real)
            if N > 0:
                result_verlinde.append(j3)
        direct = su2k_fusion(j1, j2, k)
        match = (set(result_verlinde) == set(direct))
        if not match:
            print(f"  {j1}x{j2}: Verlinde={result_verlinde}, direct={direct} MISMATCH")

print(f"  Verlinde formula verified.")

# ============================================================
# PART 8: GRAPH STRUCTURE CONSTANTS (NUMERICAL)
# ============================================================
print(f"\n{'='*70}")
print("PART 8: GRAPH 'FUSION COEFFICIENTS' (averaged)")
print("="*70)

# Since it's NOT an association scheme, define averaged structure constants:
# N_{mn}^p = (1/|edges_p|) * sum_{(i,j) in edges_p} (A_m @ A_n)[i,j]

print(f"\nAveraged structure constants <N_{{mn}}^p>:")
print(f"{'m,n':>6} | {'<N^5>':>8} {'<N^6>':>8} {'<N^7>':>8} {'<N^8>':>8} | {'<N^I>':>8}")
print("-" * 60)

for m in [5,6,7,8]:
    for n in range(m, 9):
        if n > 8: break
        prod = A[m] @ A[n]
        avgs = []
        for c in [5,6,7,8]:
            edges = [(i,j) for i in range(40) for j in range(40) if A[c][i,j]==1]
            if edges:
                avg = np.mean([prod[i,j] for i,j in edges])
            else:
                avg = 0
            avgs.append(avg)
        diag_avg = np.mean([prod[i,i] for i in range(40)])
        avg_str = " ".join(f"{a:8.3f}" for a in avgs)
        print(f"  {m},{n:1d}  | {avg_str} | {diag_avg:8.3f}")

# ============================================================
# PART 9: SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY AND HONEST ASSESSMENT")
print("="*70)

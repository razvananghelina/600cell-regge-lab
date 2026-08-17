"""
exp476i: Is "3 missing tangent directions" trivial or non-trivial?
===================================================================

CRITICAL QUESTION: Is the number 3 just codim(orbit) in Gr(4,8),
or is it genuinely related to N_gen?

TESTS:
  1. RANDOM 40 points on Gr(4,8) -- tangent rank should be 16
  2. A4/H2 pair (h=5, rank 4, Gr(2,4)=dim 4) -- how many missing?
  3. Random subsets of Coxeter orderings -- does 3 change?

If random gives 16 and A4/H2 gives a different number:
  => 3 is NON-TRIVIAL and specific to E8/H4
If random also gives 13:
  => 3 is a geometric artifact, NOT related to N_gen

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: CRITICAL VERIFICATION
"""

import numpy as np
from itertools import permutations as iterperms
from math import cos, sin, pi, sqrt

phi = (1+sqrt(5))/2; a1=5; b1=6; h=30; N=120

def sym_vec(M, n):
    """Vectorize nxn symmetric matrix."""
    v = []
    for i in range(n):
        for j in range(i, n):
            v.append(M[i,j]*sqrt(2) if i!=j else M[i,j])
    return np.array(v)

def get_tangent_rank(projectors, n):
    """Compute tangent rank of a set of rank-k projectors in R^n."""
    if len(projectors) < 2:
        return 0
    k = int(round(np.trace(projectors[0])))
    dim_gr = k * (n - k)
    dim_sym2 = n * (n + 1) // 2

    # Vectorize
    vecs = np.array([sym_vec(P, n) for P in projectors])

    # Full linear rank
    full_rank = np.linalg.matrix_rank(vecs, tol=0.01)

    # Tangent rank at first projector
    P0 = projectors[0]
    e_vals, e_vecs = np.linalg.eigh(P0)
    idx_1 = [i for i in range(n) if abs(e_vals[i]-1)<0.01]
    idx_0 = [i for i in range(n) if abs(e_vals[i])<0.01]
    V1 = e_vecs[:,idx_1]
    V0 = e_vecs[:,idx_0]

    # Tangent basis
    tang = []
    for a in range(k):
        for b in range(n-k):
            X = np.zeros((k, n-k)); X[a,b] = 1
            dP = V1 @ X @ V0.T + V0 @ X.T @ V1.T
            tang.append(sym_vec(dP, n))
    tang_mat = np.array(tang)

    # Center and project
    centered = vecs - sym_vec(P0, n)
    proj_tang = centered @ tang_mat.T
    tang_rank = np.linalg.matrix_rank(proj_tang, tol=0.01)

    return full_rank, tang_rank, dim_gr

# ============================================================
# TEST 1: RANDOM POINTS ON Gr(4,8)
# ============================================================
print("="*70)
print("TEST 1: RANDOM 40 POINTS ON Gr(4,8)")
print("="*70)

np.random.seed(123)
n_trials = 10
for trial in range(n_trials):
    random_projs = []
    for _ in range(40):
        # Random 4-plane in R^8: take random 8x4 matrix, orthogonalize
        M = np.random.randn(8, 4)
        Q, _ = np.linalg.qr(M)
        P = Q[:,:4] @ Q[:,:4].T
        random_projs.append(P)

    full_r, tang_r, dim_gr = get_tangent_rank(random_projs, 8)
    n_missing = dim_gr - tang_r
    print(f"  Trial {trial}: full_rank={full_r}, tang_rank={tang_r}/{dim_gr}, missing={n_missing}")

print(f"\n  Expected for random: tang_rank = {4*4} = dim(Gr(4,8)), missing = 0")

# ============================================================
# TEST 2: FEWER RANDOM POINTS
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: HOW MANY RANDOM POINTS NEEDED TO SPAN Gr(4,8)?")
print("="*70)

for n_pts in [5, 10, 16, 17, 20, 30, 40]:
    ranks = []
    for trial in range(5):
        rp = []
        for _ in range(n_pts):
            M = np.random.randn(8, 4)
            Q, _ = np.linalg.qr(M)
            rp.append(Q[:,:4] @ Q[:,:4].T)
        _, tr, dg = get_tangent_rank(rp, 8)
        ranks.append(tr)
    avg = np.mean(ranks)
    print(f"  {n_pts:3d} random points: avg tang_rank = {avg:.1f} / {dg}")

# ============================================================
# TEST 3: A4/H2 PAIR (rank 4, h=5)
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: A4/H2 -- rank 4, Coxeter number 5")
print("="*70)

# A4 root system: e_i - e_j for 1 <= i != j <= 5 in R^5
# But realized in R^4 as the hyperplane sum(x_i)=0 in R^5
# Simple roots: alpha_i = e_i - e_{i+1} for i=1,...,4

simple_A4 = np.zeros((4, 4))
# In R^4: use the standard A4 embedding
# alpha_1 = (1,-1,0,0), alpha_2 = (0,1,-1,0), alpha_3 = (0,0,1,-1)
# alpha_4 needs care. A4 in R^4:
# Actually, A4 has rank 4. Simple roots in R^4:
simple_A4[0] = [1,-1,0,0]
simple_A4[1] = [0,1,-1,0]
simple_A4[2] = [0,0,1,-1]
# For A4, alpha_4 = e_4 - e_5, but in R^4 we need a different embedding
# Standard: embed A4 in R^5 with constraint sum=0, then project to R^4
# Simpler: use the Cartan matrix
# A4 Cartan: tridiagonal with 2 on diagonal, -1 on off-diagonal
# Let's just use R^5 constrained to sum=0

# Actually, let me use R^5 and project to the 4D hyperplane sum=0
simple_A4_5 = np.zeros((4, 5))
for i in range(4):
    simple_A4_5[i, i] = 1
    simple_A4_5[i, i+1] = -1

print(f"A4 simple roots (in R^5):")
for i in range(4):
    print(f"  alpha_{i} = {simple_A4_5[i].tolist()}")

# A4 Cartan matrix
cartan_A4 = np.zeros((4,4), dtype=int)
for i in range(4):
    for j in range(4):
        cartan_A4[i,j] = int(round(2 * simple_A4_5[i] @ simple_A4_5[j] / (simple_A4_5[i] @ simple_A4_5[i])))
print(f"A4 Cartan matrix:\n{cartan_A4}")

# Reflections in R^5
refs_A4 = []
for i in range(4):
    a = simple_A4_5[i]
    refs_A4.append(np.eye(5) - 2*np.outer(a,a)/(a@a))

# Coxeter elements and H2 projectors
# A4 exponents: {1, 2, 3, 4}, h=5
# H2 exponents: {1, 4} (the Coxeter exponents of the dihedral group I2(5))
# H2' exponents: {2, 3}
# H2 projector: 2D eigenspace for exp 2*pi*i*1/5 and 2*pi*i*4/5

def make_cox_A4(perm):
    C = np.eye(5)
    for i in perm:
        C = refs_A4[i] @ C
    return C

def get_proj_A4(C, h_val=5):
    """Get H2 projector from A4 Coxeter element."""
    if not np.allclose(np.linalg.matrix_power(C, h_val), np.eye(5), atol=1e-8):
        return None
    ev, evec = np.linalg.eig(C)
    args = np.angle(ev)
    exps = [(round(args[i]*h_val/(2*pi)) % h_val or h_val, i) for i in range(5)]

    # H2 exponents: {1, 4} (= {1, h-1})
    h2_idx = [i for m, i in exps if m in {1, 4}]
    if len(h2_idx) != 2:
        return None

    v_raw = evec[:, h2_idx]
    basis = []
    used = set()
    for k in range(2):
        if k in used: continue
        v = v_raw[:, k]
        if np.max(np.abs(np.imag(v))) > 0.01:
            basis.extend([np.real(v), np.imag(v)])
            for k2 in range(k+1, 2):
                if k2 not in used and np.allclose(v_raw[:,k2], np.conj(v)):
                    used.add(k2); break
        else:
            basis.append(np.real(v))
    if len(basis) != 2:
        return None
    B = np.array(basis).T
    Q, _ = np.linalg.qr(B)
    return Q[:,:2]

# Find all distinct A4/H2 decompositions
projs_A4 = []
for perm in iterperms(range(4)):
    C = make_cox_A4(perm)
    P = get_proj_A4(C)
    if P is None: continue
    PP = P @ P.T  # 5x5 projector
    is_new = True
    for Po in projs_A4:
        svs = np.linalg.svd(Po[:,:2].T @ P, compute_uv=False) if Po.shape[1]==2 else [0]
        # Compare as projector matrices
        diff = np.linalg.norm(PP - Po)
        if diff < 0.01:
            is_new = False
            break
    if is_new:
        projs_A4.append(PP)

n_A4 = len(projs_A4)
print(f"\nA4/H2 decompositions found: {n_A4}")
print(f"  Compare: for E8/H4 we had 40 = rank(E8)*a1 = 8*5")
print(f"  For A4/H2: rank(A4)=4, a1_analog=?")
if n_A4 > 0 and n_A4 % 4 == 0:
    print(f"  {n_A4}/rank(A4) = {n_A4//4}")

if n_A4 >= 2:
    # Compute tangent rank for A4/H2
    # Gr(2,5) has dim 2*3 = 6
    full_r, tang_r, dim_gr = get_tangent_rank(projs_A4, 5)
    n_missing_A4 = dim_gr - tang_r
    print(f"\n  A4/H2 on Gr(2,5) (dim {dim_gr}):")
    print(f"  Full rank: {full_r}")
    print(f"  Tangent rank: {tang_r}")
    print(f"  Missing tangent directions: {n_missing_A4}")
    print(f"  Formula: dim(Gr) - tang_rank = {dim_gr} - {tang_r} = {n_missing_A4}")

    # Gram matrix rank
    T_A4 = np.array([[np.trace(projs_A4[i] @ projs_A4[j])
                       for j in range(n_A4)] for i in range(n_A4)])
    rank_T_A4 = np.linalg.matrix_rank(T_A4, tol=0.01)
    print(f"  rank(T_A4) = {rank_T_A4}")
    print(f"  ker(T_A4) = {n_A4} - {rank_T_A4} = {n_A4 - rank_T_A4}")

    # Overlap values
    overlaps_A4 = set()
    for i in range(n_A4):
        for j in range(i+1, n_A4):
            overlaps_A4.add(round(np.trace(projs_A4[i] @ projs_A4[j]), 4))
    print(f"  Overlap values: {sorted(overlaps_A4)}")

# ============================================================
# TEST 4: E8/H4 WITH ONLY 17 POINTS (just above dim Gr = 16)
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: E8/H4 SUBSETS -- HOW MANY POINTS NEEDED FOR RANK 13?")
print("="*70)

# Build the E8/H4 projectors (reuse infrastructure)
E8_roots = []
for i in range(8):
    for j in range(i+1,8):
        for si in [1,-1]:
            for sj in [1,-1]:
                v=np.zeros(8);v[i]=si;v[j]=sj;E8_roots.append(v)
for bits in range(256):
    signs=[(-1)**((bits>>k)&1) for k in range(8)]
    if sum(1 for s in signs if s<0)%2==0:
        E8_roots.append(np.array(signs)/2.0)
simple_E8=np.zeros((8,8))
simple_E8[0]=[1,-1,0,0,0,0,0,0];simple_E8[1]=[0,1,-1,0,0,0,0,0]
simple_E8[2]=[0,0,1,-1,0,0,0,0];simple_E8[3]=[0,0,0,1,-1,0,0,0]
simple_E8[4]=[0,0,0,0,1,-1,0,0];simple_E8[5]=[0,0,0,0,0,1,-1,0]
simple_E8[6]=[0,0,0,0,0,1,1,0]; simple_E8[7]=[-0.5]*8
refs_E8=[np.eye(8)-2*np.outer(simple_E8[i],simple_E8[i])/np.dot(simple_E8[i],simple_E8[i]) for i in range(8)]
def make_cox_E8(perm):
    C=np.eye(8)
    for i in perm: C=refs_E8[i]@C
    return C
def get_proj_E8(C):
    if not np.allclose(np.linalg.matrix_power(C,30),np.eye(8),atol=1e-8): return None
    ev,evec=np.linalg.eig(C);args=np.angle(ev)
    exps=[(round(args[i]*30/(2*pi))%30 or 30,i) for i in range(8)]
    h4=[i for m,i in exps if m in{1,11,19,29}]
    if len(h4)!=4: return None
    vr=evec[:,h4];b=[];used=set()
    for k in range(4):
        if k in used: continue
        v=vr[:,k]
        if np.max(np.abs(np.imag(v)))>0.01:
            b.extend([np.real(v),np.imag(v)])
            for k2 in range(k+1,4):
                if k2 not in used and np.allclose(vr[:,k2],np.conj(v)):
                    used.add(k2);break
        else: b.append(np.real(v))
    if len(b)!=4: return None
    B=np.array(b).T;Q,_=np.linalg.qr(B)
    return Q[:,:4]@(Q[:,:4]).T

projs_E8=[];np.random.seed(42)
perms_E8=list(iterperms(range(8)));np.random.shuffle(perms_E8)
for perm in perms_E8:
    P=get_proj_E8(make_cox_E8(perm))
    if P is None: continue
    is_new=True
    for Po in projs_E8:
        if np.linalg.norm(P-Po)<0.01: is_new=False; break
    if is_new: projs_E8.append(P)
    if len(projs_E8)>=40: break

print(f"Built {len(projs_E8)} E8/H4 projectors")

# Test subsets of different sizes
for n_pts in [5, 10, 14, 15, 16, 17, 20, 30, 40]:
    if n_pts > len(projs_E8): break
    # Try multiple random subsets
    ranks = []
    for trial in range(min(10, len(projs_E8))):
        np.random.seed(trial*1000 + n_pts)
        subset = [projs_E8[i] for i in np.random.choice(len(projs_E8), n_pts, replace=False)]
        _, tr, dg = get_tangent_rank(subset, 8)
        ranks.append(tr)
    avg = np.mean(ranks)
    min_r = min(ranks)
    max_r = max(ranks)
    print(f"  {n_pts:3d} E8/H4 points: tang_rank = {min_r:.0f}-{max_r:.0f} (avg {avg:.1f}) / {dg}  missing = {dg-avg:.1f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY: IS 3 TRIVIAL?")
print("="*70)

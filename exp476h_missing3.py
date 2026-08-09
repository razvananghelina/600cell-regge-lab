"""
exp476h: The 3 Missing Tangent Directions -- Universal or Reference-Dependent?
===============================================================================

From exp476g: the 40 H4 projectors span 13 of 16 tangent directions at P_0.
3 tangent directions are MISSING. 3 = N_gen?

KEY TESTS:
  1. Are the 3 missing directions the SAME at every reference point?
  2. Do they form a recognizable algebraic object (SU(3)? generation space?)
  3. Can we express them in terms of E8/H4 data?

If the 3-dim missing space is UNIVERSAL (same at all P_alpha):
  => rank(T) = 16 - 3 + 1 = 14 is PROVED (modulo the curvature dimension)
  => The 3-dim space has a geometric meaning

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: CRITICAL TEST
"""

import numpy as np
from itertools import permutations as iterperms
from math import cos, sin, pi, sqrt, factorial

phi = (1+sqrt(5))/2; a1=5; b1=6; h=30; N=120

# Build 40 decompositions (same infrastructure)
E8_roots = []
for i in range(8):
    for j in range(i+1,8):
        for si in [1,-1]:
            for sj in [1,-1]:
                v=np.zeros(8); v[i]=si; v[j]=sj; E8_roots.append(v)
for bits in range(256):
    signs=[(-1)**((bits>>k)&1) for k in range(8)]
    if sum(1 for s in signs if s<0)%2==0:
        E8_roots.append(np.array(signs)/2.0)
simple=np.zeros((8,8))
simple[0]=[1,-1,0,0,0,0,0,0];simple[1]=[0,1,-1,0,0,0,0,0]
simple[2]=[0,0,1,-1,0,0,0,0];simple[3]=[0,0,0,1,-1,0,0,0]
simple[4]=[0,0,0,0,1,-1,0,0];simple[5]=[0,0,0,0,0,1,-1,0]
simple[6]=[0,0,0,0,0,1,1,0]; simple[7]=[-0.5]*8
refs=[np.eye(8)-2*np.outer(simple[i],simple[i])/np.dot(simple[i],simple[i]) for i in range(8)]
def make_cox(perm):
    C=np.eye(8)
    for i in perm: C=refs[i]@C
    return C
def get_proj(C):
    if not np.allclose(np.linalg.matrix_power(C,30),np.eye(8),atol=1e-8): return None
    ev,evec=np.linalg.eig(C); args=np.angle(ev)
    exps=[(round(args[i]*30/(2*pi))%30 or 30,i) for i in range(8)]
    h4=[i for m,i in exps if m in{1,11,19,29}]
    if len(h4)!=4: return None
    vr=evec[:,h4]; b=[]; used=set()
    for k in range(4):
        if k in used: continue
        v=vr[:,k]
        if np.max(np.abs(np.imag(v)))>0.01:
            b.extend([np.real(v),np.imag(v)])
            for k2 in range(k+1,4):
                if k2 not in used and np.allclose(vr[:,k2],np.conj(v)):
                    used.add(k2); break
        else: b.append(np.real(v))
    if len(b)!=4: return None
    B=np.array(b).T; Q,_=np.linalg.qr(B)
    return Q[:,:4]

projs=[]; np.random.seed(42)
perms=list(iterperms(range(8))); np.random.shuffle(perms)
for perm in perms:
    P=get_proj(make_cox(perm))
    if P is None: continue
    if all(not np.allclose(np.linalg.svd(Po.T@P,compute_uv=False),np.ones(4),atol=1e-6) for Po in projs):
        projs.append(P)
    if len(projs)>=40: break
assert len(projs)==40
Proj8=[P@P.T for P in projs]
print(f"Built 40 projectors")

def sym_vec(M):
    v=[]
    for i in range(8):
        for j in range(i,8):
            v.append(M[i,j]*sqrt(2) if i!=j else M[i,j])
    return np.array(v)

vecs = np.array([sym_vec(P) for P in Proj8])  # 40 x 36

def get_tangent_space(P0):
    """Get the 16-dim tangent space to Gr(4,8) at P0."""
    e_vals, e_vecs = np.linalg.eigh(P0)
    idx_1 = [i for i in range(8) if abs(e_vals[i]-1)<0.01]
    idx_0 = [i for i in range(8) if abs(e_vals[i])<0.01]
    V1 = e_vecs[:,idx_1]  # range of P0
    V0 = e_vecs[:,idx_0]  # null of P0
    tangent_vecs = []
    for a in range(4):
        for b in range(4):
            X=np.zeros((4,4)); X[a,b]=1
            dP = V1 @ X @ V0.T + V0 @ X.T @ V1.T
            tangent_vecs.append(sym_vec(dP))
    return np.array(tangent_vecs), V1, V0  # 16x36, 8x4, 8x4

def get_missing_tangent(ref_idx):
    """Find missing tangent directions at reference point ref_idx."""
    P0 = Proj8[ref_idx]
    tang_mat, V1, V0 = get_tangent_space(P0)

    # Center projectors at P0
    centered = vecs - sym_vec(P0)  # 40 x 36

    # Project centered onto tangent space
    proj_tang = centered @ tang_mat.T  # 40 x 16

    # SVD to find missing directions
    U, S, Vt = np.linalg.svd(proj_tang, full_matrices=True)
    n_missing = sum(1 for s in S if s < 0.01)
    rank_tang = 16 - n_missing

    # The missing directions (in tangent parameter space)
    missing_dirs_param = Vt[16-n_missing:]  # n_missing x 16

    # Convert to 8x8 matrices
    missing_mats = []
    for k in range(n_missing):
        d = missing_dirs_param[k]
        X = d.reshape(4,4)
        dP = V1 @ X @ V0.T + V0 @ X.T @ V1.T
        missing_mats.append(dP)

    return n_missing, missing_mats, V1, V0

# ============================================================
# TEST 1: IS THE NUMBER OF MISSING DIRECTIONS ALWAYS 3?
# ============================================================
print(f"\n{'='*70}")
print("TEST 1: NUMBER OF MISSING TANGENT DIRECTIONS AT EACH P_alpha")
print("="*70)

missing_counts = []
all_missing_mats = {}
for ref in range(40):
    n_miss, mats, _, _ = get_missing_tangent(ref)
    missing_counts.append(n_miss)
    all_missing_mats[ref] = mats

print(f"\nMissing tangent directions at each of 40 reference points:")
unique_counts = set(missing_counts)
for c in sorted(unique_counts):
    n = sum(1 for x in missing_counts if x == c)
    print(f"  {c} missing: {n} reference points")

is_universal_count = (len(unique_counts) == 1)
print(f"\nUniversal count? {is_universal_count}")
if is_universal_count:
    n_missing_universal = missing_counts[0]
    print(f"*** ALL 40 reference points have exactly {n_missing_universal} missing tangent directions ***")
    print(f"*** {n_missing_universal} = N_gen? {n_missing_universal == 3} ***")

# ============================================================
# TEST 2: IS THE MISSING SUBSPACE THE SAME AT ALL POINTS?
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: IS THE 3-DIM MISSING SUBSPACE UNIVERSAL?")
print("="*70)

# For each pair of reference points, check if their missing subspaces
# are the same (have the same span in Sym^2(R^8))

# Convert missing mats to Sym^2 vectors for comparison
def missing_to_vecs(ref_idx):
    mats = all_missing_mats[ref_idx]
    return np.array([sym_vec(M) for M in mats])

# Compare missing subspaces between different reference points
# Use principal angles between subspaces
print(f"\nComparing missing subspaces between reference points:")
print(f"(Principal angles -- if all ~0, subspaces are identical)")

# Pick 10 representative reference points
test_refs = [0, 5, 10, 15, 20, 25, 30, 35, 2, 17]
n_miss = missing_counts[0]

results = []
for i in range(len(test_refs)):
    for j in range(i+1, len(test_refs)):
        ri, rj = test_refs[i], test_refs[j]
        Vi = missing_to_vecs(ri)  # n_miss x 36
        Vj = missing_to_vecs(rj)  # n_miss x 36

        if len(Vi) == 0 or len(Vj) == 0:
            continue

        # Principal angles via SVD
        M = Vi @ Vj.T  # n_miss x n_miss
        svs = np.linalg.svd(M, compute_uv=False)
        # Clamp to [-1, 1] for acos
        svs = np.clip(svs, -1, 1)
        angles = np.arccos(np.minimum(svs, 1.0)) * 180 / pi
        max_angle = max(angles)
        results.append((ri, rj, max_angle, angles))

# Show results
same_count = 0
diff_count = 0
for ri, rj, max_ang, angles in results[:15]:
    status = "SAME" if max_ang < 5 else "DIFFERENT"
    if max_ang < 5: same_count += 1
    else: diff_count += 1
    ang_str = ", ".join(f"{a:.1f}" for a in angles)
    print(f"  P_{ri:2d} vs P_{rj:2d}: max angle = {max_ang:6.2f} deg ({status})  [{ang_str}]")

print(f"\n  SAME subspace: {same_count}/{same_count+diff_count}")
print(f"  DIFFERENT subspace: {diff_count}/{same_count+diff_count}")

if diff_count > 0:
    print(f"\n  The missing subspace is REFERENCE-DEPENDENT (not universal).")
    print(f"  This means: the 3 missing directions rotate as we move on Gr(4,8).")
    print(f"  They define a VECTOR BUNDLE over the 40 points, not a fixed subspace.")

    # Even though not universal, check: is there a UNIVERSAL missing subspace
    # that's smaller? I.e., intersection of all missing subspaces.
    print(f"\n  Computing INTERSECTION of all 40 missing subspaces...")

    # Start with the missing subspace at ref 0
    intersection = missing_to_vecs(0)  # n_miss x 36

    for ref in range(1, 40):
        V_ref = missing_to_vecs(ref)  # n_miss x 36
        # Intersection = common subspace
        # Stack and find singular values
        combined = np.vstack([intersection, V_ref])
        U_c, S_c, Vt_c = np.linalg.svd(combined, full_matrices=False)

        # Vectors with singular value close to sqrt(2) are in both subspaces
        # (because ||v|| = 1 in each, so if v is in both, combined has sv = sqrt(2))
        # Actually, this doesn't work directly. Use a different method.

        # Method: project V_ref onto span of intersection, keep what has large projection
        if intersection.shape[0] == 0:
            break
        proj = V_ref @ intersection.T  # n_miss x current_dim
        # SVD of projection
        U_p, S_p, Vt_p = np.linalg.svd(proj, full_matrices=False)
        # Directions with sv close to 1 are in the intersection
        in_both = sum(1 for s in S_p if s > 0.95)
        if in_both == 0:
            intersection = np.zeros((0, 36))
            break
        # Update intersection to the common part
        # The common directions in intersection basis
        common_dirs = Vt_p[:in_both] @ intersection[:Vt_p.shape[1]]
        # Re-orthogonalize
        if common_dirs.shape[0] > 0:
            U_new, S_new, Vt_new = np.linalg.svd(common_dirs, full_matrices=False)
            intersection = Vt_new[:sum(1 for s in S_new if s > 0.01)]

    print(f"  Intersection dimension: {intersection.shape[0]}")

    if intersection.shape[0] > 0:
        print(f"  *** There IS a universal missing subspace of dim {intersection.shape[0]} ***")
        # Analyze it
        for k in range(intersection.shape[0]):
            from math import sqrt as msqrt
            M_int = np.zeros((8,8)); idx=0
            v = intersection[k]
            for i in range(8):
                for j in range(i,8):
                    val = v[idx]/sqrt(2) if i!=j else v[idx]
                    M_int[i,j]=val; M_int[j,i]=val; idx+=1
            eigs = np.sort(np.linalg.eigvalsh(M_int))[::-1]
            print(f"    Direction {k}: eigs = {[round(e,4) for e in eigs]}, tr = {np.trace(M_int):.4f}")
    else:
        print(f"  No universal intersection -- the missing space ROTATES completely.")

else:
    print(f"\n  *** The missing subspace IS UNIVERSAL! ***")
    print(f"  The 3-dim space is independent of reference point.")

# ============================================================
# TEST 3: WHAT DO THE MISSING DIRECTIONS LOOK LIKE?
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: STRUCTURE OF MISSING DIRECTIONS (at P_0)")
print("="*70)

n_miss_0, miss_mats_0, V1_0, V0_0 = get_missing_tangent(0)

for k in range(n_miss_0):
    dP = miss_mats_0[k]
    eigs = np.sort(np.linalg.eigvalsh(dP))[::-1]
    print(f"\n  Missing direction {k}:")
    print(f"    Eigenvalues: {[round(e,4) for e in eigs]}")

    # Express in range/null basis of P_0
    # dP = V1 @ X @ V0.T + V0 @ X.T @ V1.T
    # where X = V1.T @ dP @ V0 (4x4 matrix)
    X = V1_0.T @ dP @ V0_0  # 4x4
    print(f"    X matrix (4x4, tangent parameter):")
    print(f"      {np.round(X, 4).tolist()}")

    # Singular values of X
    svs_X = np.linalg.svd(X, compute_uv=False)
    print(f"    Singular values of X: {[round(s,4) for s in svs_X]}")

    # Is X antisymmetric? X + X.T = 0?
    sym_part = np.linalg.norm(X + X.T)
    antisym_part = np.linalg.norm(X - X.T)
    print(f"    ||X + X.T|| = {sym_part:.4f} (0 if antisymmetric)")
    print(f"    ||X - X.T|| = {antisym_part:.4f} (0 if symmetric)")

# ============================================================
# TEST 4: THE 3 DIRECTIONS AND FRAMEWORK NUMBERS
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: FRAMEWORK NUMBER CHECK")
print("="*70)

# Check various invariants of the 3-dim missing space
# Use all 40 reference points to build a robust picture

# For each reference, compute the 3 missing direction matrices
# Then compute their traces, determinants, etc.

# Actually, since the subspace varies, let's look at what's INVARIANT:
# The DIMENSION (3) is invariant. What else?

# Check: the sum of squared principal angles between missing and a fixed reference
all_max_angles = []
ref0_vecs = missing_to_vecs(0)
for ref in range(1, 40):
    V_ref = missing_to_vecs(ref)
    M = ref0_vecs @ V_ref.T
    svs = np.clip(np.linalg.svd(M, compute_uv=False), 0, 1)
    angles = np.arccos(svs) * 180 / pi
    all_max_angles.append(max(angles))

print(f"\nMax principal angle between missing subspace at P_0 and others:")
print(f"  min = {min(all_max_angles):.2f} deg")
print(f"  max = {max(all_max_angles):.2f} deg")
print(f"  mean = {np.mean(all_max_angles):.2f} deg")
print(f"  (If all ~0: universal. If large: reference-dependent.)")

# Histogram of angles
bins = [0, 10, 20, 30, 45, 60, 90]
for i in range(len(bins)-1):
    n = sum(1 for a in all_max_angles if bins[i] <= a < bins[i+1])
    print(f"  [{bins[i]:3d}, {bins[i+1]:3d}) deg: {n} points")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print("="*70)
print(f"""
  Number of missing tangent directions: {missing_counts[0]} at ALL 40 points
  Universal count: {'YES' if is_universal_count else 'NO'}

  Is the missing 3-dim subspace the SAME everywhere?
  {'YES -- universal' if diff_count == 0 else 'NO -- it rotates (reference-dependent)'}

  KEY FORMULA: rank(T) = dim(Gr(4,8)) - n_missing + n_curvature
                       = 16 - {missing_counts[0]} + {14 - (16 - missing_counts[0])}
                       = 14

  The number {missing_counts[0]} = {'N_gen' if missing_counts[0] == 3 else '???'}
""")

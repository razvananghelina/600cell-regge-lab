"""
exp476g: Proving rank(T) = 14 -- Analyzing the 22-dim complement
=================================================================

Strategy: instead of asking "why 14 in span?", ask "why 22 in complement?"

The 40 projectors span 14 dims of Sym^2(R^8) (dim 36).
The orthogonal complement has dim 22.
If we can characterize these 22 directions algebraically, we have a proof.

Key tools:
  - Each P_alpha is a RANK-4 PROJECTOR: P^2=P, P^T=P, Tr(P)=4
  - The 40 P's come from H4 eigenspaces of E8 Coxeter elements
  - The complement consists of all M in Sym^2(R^8) with Tr(M @ P_alpha)=0 for all alpha

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: PROOF ATTEMPT
"""

import numpy as np
from itertools import permutations as iterperms
from math import cos, sin, pi, sqrt, factorial
from collections import Counter

phi = (1+sqrt(5))/2; a1=5; b1=6; h=30; N=120; rank_E8=8

# Build 40 decompositions (compact, verified code)
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
E8=np.array(E8_roots)
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

# ============================================================
# PART 1: VECTORIZE AND FIND THE 14-DIM SPAN AND 22-DIM COMPLEMENT
# ============================================================
print(f"\n{'='*70}")
print("PART 1: THE 14-DIM SPAN AND 22-DIM COMPLEMENT")
print("="*70)

def sym_vec(M):
    v=[]
    for i in range(8):
        for j in range(i,8):
            v.append(M[i,j]*sqrt(2) if i!=j else M[i,j])
    return np.array(v)

def sym_unvec(v):
    M=np.zeros((8,8)); idx=0
    for i in range(8):
        for j in range(i,8):
            val=v[idx]/sqrt(2) if i!=j else v[idx]
            M[i,j]=val; M[j,i]=val; idx+=1
    return M

vecs = np.array([sym_vec(P) for P in Proj8])  # 40 x 36
U, S_vals, Vt = np.linalg.svd(vecs, full_matrices=True)

# The 14 span directions (rows of Vt[:14])
V_span = Vt[:14, :]   # 14 x 36
# The 22 complement directions
V_comp = Vt[14:, :]    # 22 x 36

print(f"Span dimension: {sum(1 for s in S_vals if s > 0.01)} = 14")
print(f"Complement dimension: {36 - sum(1 for s in S_vals if s > 0.01)} = 22")

# ============================================================
# PART 2: PROPERTIES OF THE 22 COMPLEMENT DIRECTIONS
# ============================================================
print(f"\n{'='*70}")
print("PART 2: THE 22 COMPLEMENT MATRICES")
print("="*70)

comp_mats = [sym_unvec(V_comp[k]) for k in range(22)]

# Basic properties
print(f"\nComplement basis matrices (8x8 symmetric, Frobenius-orthogonal to all P_alpha):")
for k in range(22):
    M = comp_mats[k]
    tr = np.trace(M)
    rk = np.linalg.matrix_rank(M, tol=0.01)
    eigs = np.sort(np.linalg.eigvalsh(M))[::-1]
    n_pos = sum(1 for e in eigs if e > 0.01)
    n_neg = sum(1 for e in eigs if e < -0.01)
    n_zero = 8 - n_pos - n_neg
    print(f"  C{k:2d}: tr={tr:+.4f}, rank={rk}, signature=({n_pos},{n_zero},{n_neg})")

# ============================================================
# PART 3: IS THE COMPLEMENT RELATED TO ANTISYMMETRIC MATRICES?
# ============================================================
print(f"\n{'='*70}")
print("PART 3: COMPLEMENT vs KNOWN SUBSPACES")
print("="*70)

# Sym^2(R^8) has dim 36. What natural subspaces of dim 22 exist?
# - Anti^2(R^8) has dim 28 (antisymmetric matrices), NOT inside Sym^2
# - But we can look at specific SYMMETRIC subspaces of dim 22

# Key decomposition under SO(8):
# Sym^2(R^8) = R (trace) + 35_v (traceless symmetric)
# Under smaller groups, 35_v decomposes further

# Check: is the complement entirely traceless?
comp_traces = [np.trace(M) for M in comp_mats]
print(f"\nTraces of complement matrices: {[round(t,4) for t in comp_traces]}")
n_nonzero_trace = sum(1 for t in comp_traces if abs(t) > 0.01)
print(f"Number with non-zero trace: {n_nonzero_trace}")

# Check: how many complement directions are diagonal?
diag_projections = []
for k in range(22):
    M = comp_mats[k]
    diag_part = np.diag(np.diag(M))
    offdiag_part = M - diag_part
    diag_norm = np.linalg.norm(sym_vec(diag_part))
    offdiag_norm = np.linalg.norm(sym_vec(offdiag_part))
    total_norm = np.linalg.norm(V_comp[k])
    diag_projections.append(diag_norm / total_norm if total_norm > 0 else 0)

print(f"\nDiagonal content of complement:")
print(f"  Fraction diagonal: {[round(d,3) for d in diag_projections]}")

# ============================================================
# PART 4: THE E8 ROOT STRUCTURE IN THE COMPLEMENT
# ============================================================
print(f"\n{'='*70}")
print("PART 4: E8 ROOT PRODUCTS IN SPAN vs COMPLEMENT")
print("="*70)

# Each E8 root r gives a rank-1 matrix r@r.T in Sym^2(R^8)
# These 240 rank-1 matrices span some subspace
# Question: how much of Sym^2 do E8 roots span?

root_vecs = np.array([sym_vec(np.outer(r,r)) for r in E8])  # 240 x 36
rank_roots = np.linalg.matrix_rank(root_vecs, tol=0.01)
print(f"\nE8 root outer products span {rank_roots} dims of Sym^2(R^8)")

# Project E8 root products onto span and complement
proj_on_span = root_vecs @ V_span.T  # 240 x 14
proj_on_comp = root_vecs @ V_comp.T  # 240 x 22

rank_on_span = np.linalg.matrix_rank(proj_on_span, tol=0.01)
rank_on_comp = np.linalg.matrix_rank(proj_on_comp, tol=0.01)
print(f"E8 roots project onto span: rank {rank_on_span} (of 14)")
print(f"E8 roots project onto complement: rank {rank_on_comp} (of 22)")

# ============================================================
# PART 5: SIMPLE ROOT PRODUCTS
# ============================================================
print(f"\n{'='*70}")
print("PART 5: SIMPLE ROOT PRODUCTS AND THE SPAN")
print("="*70)

# The 8 simple roots give 8 rank-1 matrices (diagonal of Sym^2)
# Plus C(8,2)=28 cross products (alpha_i alpha_j^T + alpha_j alpha_i^T)
# Total: 36 products (spanning all of Sym^2 if independent)

simple_products = []
for i in range(8):
    simple_products.append(sym_vec(np.outer(simple[i], simple[i])))
for i in range(8):
    for j in range(i+1, 8):
        cross = np.outer(simple[i], simple[j]) + np.outer(simple[j], simple[i])
        simple_products.append(sym_vec(cross))

simple_mat = np.array(simple_products)
rank_simple = np.linalg.matrix_rank(simple_mat, tol=0.01)
print(f"Simple root products (self + cross): {len(simple_products)} matrices, rank {rank_simple}")

# ============================================================
# PART 6: THE CRITICAL TEST -- QUADRATIC CONSTRAINT
# ============================================================
print(f"\n{'='*70}")
print("PART 6: QUADRATIC CONSTRAINTS P^2 = P")
print("="*70)

# Each P_alpha satisfies P^2 = P. When vectorized:
# vec(P^2) = vec(P) imposes constraints.
# But P^2 in the vectorized form involves the product of P with itself,
# which is a QUADRATIC relation. Not directly useful for linear rank.

# However: consider the LINEARIZATION at a reference projector P_0.
# If P = P_0 + epsilon*dP, then P^2 = P gives:
# P_0*dP + dP*P_0 = dP (to first order)
# This is a LINEAR constraint on dP.
# dP must satisfy: (P_0 - I/2)*dP + dP*(P_0 - I/2) = 0
# i.e., {P_0 - I/2, dP} = 0 (anticommutator)
# i.e., dP anticommutes with Q_0 = P_0 - I/2

# The space of symmetric dP anticommuting with Q_0 has dim = ?
# Q_0 has eigenvalues +1/2 (4 times) and -1/2 (4 times)
# In the +1/2 eigenspace (dim 4): dP maps within
# In the -1/2 eigenspace (dim 4): dP maps within
# Between eigenspaces: {Q_0, dP} = 0 means dP maps +1/2 to -1/2 and vice versa

# The tangent space to Gr(4,8) at P_0 consists of:
# dP such that P_0*dP + dP*P_0 = dP, Tr(dP) = 0, dP = dP^T
# This gives: dP = P_0*X*(I-P_0) + (I-P_0)*X^T*P_0 for arbitrary 4x4 matrix X
# dim = 4*4 = 16. This is the TANGENT SPACE to Gr(4,8).

# Now: the 40 projectors span 14 dims. So in the tangent space (dim 16),
# the 40 projectors span at most 14 tangent directions.
# The 2 "missing" tangent directions are where the 40 points are constant.

# Let's FIND these 2 missing tangent directions at P_0.

print(f"\nTangent space analysis at P_0:")
P0 = Proj8[0]
Q0 = P0 - np.eye(8)/2

# Compute tangent space: dP = P0*X*(I-P0) + (I-P0)*X^T*P0
# for general 4x4 matrix X (using the eigenvectors of P0)
e_vals, e_vecs = np.linalg.eigh(P0)
# Eigenvalues are 0 (4 times) and 1 (4 times)
idx_1 = [i for i in range(8) if abs(e_vals[i] - 1) < 0.01]
idx_0 = [i for i in range(8) if abs(e_vals[i]) < 0.01]
V1 = e_vecs[:, idx_1]  # 8x4, columns span Range(P0)
V0 = e_vecs[:, idx_0]  # 8x4, columns span Null(P0)

# Tangent vectors: dP = V1 @ X @ V0.T + V0 @ X.T @ V1.T
# for arbitrary 4x4 matrix X
# Vectorize this: 16 basis tangent vectors
tangent_vecs = []
for a in range(4):
    for b in range(4):
        X = np.zeros((4,4)); X[a,b] = 1
        dP = V1 @ X @ V0.T + V0 @ X.T @ V1.T
        tangent_vecs.append(sym_vec(dP))

tangent_mat = np.array(tangent_vecs)  # 16 x 36
rank_tangent = np.linalg.matrix_rank(tangent_mat, tol=0.01)
print(f"  Tangent space at P0: {len(tangent_vecs)} vectors, rank {rank_tangent}")

# Project the 40 projectors onto the tangent space
# First center at P0
centered_vecs = vecs - sym_vec(P0)  # 40 x 36 (centered at P0)
proj_tangent = centered_vecs @ tangent_mat.T  # 40 x 16
rank_proj_tangent = np.linalg.matrix_rank(proj_tangent, tol=0.01)
print(f"  Centered projectors projected on tangent: rank {rank_proj_tangent} out of 16")
print(f"  Missing tangent directions: {16 - rank_proj_tangent}")

# Find the missing directions
U_t, S_t, Vt_t = np.linalg.svd(proj_tangent, full_matrices=True)
print(f"  Singular values of tangent projection:")
for i, s in enumerate(S_t):
    print(f"    sv[{i:2d}] = {s:.6f}{'  <-- zero' if s < 0.01 else ''}")

# The missing tangent directions
n_missing = sum(1 for s in S_t if s < 0.01)
if n_missing > 0:
    # The missing directions in tangent parameter space
    missing_idx = [i for i, s in enumerate(S_t) if s < 0.01]
    print(f"\n  Missing tangent directions: {n_missing}")

    # Convert back to 8x8 matrices
    for idx in missing_idx:
        # The direction in the 16-dim tangent parameter space
        direction = Vt_t[idx]  # 16-dim vector
        # Convert to X matrix (4x4)
        X = direction.reshape(4, 4)
        dP = V1 @ X @ V0.T + V0 @ X.T @ V1.T
        print(f"\n  Missing direction {idx}:")
        print(f"    dP matrix (8x8):")
        # Eigenvalues of dP
        dP_eigs = np.sort(np.linalg.eigvalsh(dP))[::-1]
        print(f"    Eigenvalues: {[round(e,4) for e in dP_eigs]}")
        print(f"    Trace: {np.trace(dP):.6f}")
        print(f"    Frobenius norm: {np.linalg.norm(dP):.6f}")

        # Check if dP is related to E8 simple roots or Coxeter
        # Project onto simple root basis
        for k, alpha in enumerate(simple):
            overlap = np.trace(dP @ np.outer(alpha, alpha))
            if abs(overlap) > 0.01:
                print(f"    Overlap with alpha_{k} outer: {overlap:.4f}")

# ============================================================
# PART 7: NORMAL SPACE ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("PART 7: NORMAL (NON-TANGENT) COMPONENT OF COMPLEMENT")
print("="*70)

# The 22-dim complement splits into:
# - Tangent complement: directions in Gr(4,8) tangent space not hit by 40 points
# - Normal complement: directions orthogonal to Gr(4,8) at P0

# Normal space at P0: everything in Sym^2 NOT tangent to Gr(4,8)
# dim = 36 - 16 = 20

# Check: how much of the 22-dim complement is tangent vs normal?
proj_comp_on_tangent = V_comp @ tangent_mat.T  # 22 x 16
rank_comp_tangent = np.linalg.matrix_rank(proj_comp_on_tangent, tol=0.01)
print(f"  22-dim complement projected on tangent space: rank {rank_comp_tangent}")
print(f"  22-dim complement in normal space: {22 - rank_comp_tangent}")

# So: 22 = rank_comp_tangent (tangent) + (22-rank_comp_tangent) (normal)
# Expected: 20 normal + 2 tangent = 22
# This would confirm: the 40 points miss exactly 2 tangent directions

print(f"\n  DECOMPOSITION of 22-dim complement:")
print(f"    Tangent complement: {rank_comp_tangent} (missing Gr directions)")
print(f"    Normal complement: {22 - rank_comp_tangent} (Gr codimension)")
print(f"    Expected: 2 tangent + 20 normal = 22")
print(f"    Match: {rank_comp_tangent == 2 and (22-rank_comp_tangent) == 20}")

# ============================================================
# PART 8: WHAT ARE THE 2 MISSING TANGENT DIRECTIONS?
# ============================================================
if rank_comp_tangent >= 1:
    print(f"\n{'='*70}")
    print("PART 8: THE MISSING TANGENT DIRECTIONS")
    print("="*70)

    # Find the tangent directions that ARE in the complement
    # These are the directions along Gr(4,8) that the 40 points don't explore

    # SVD of the 22x16 projection
    U_ct, S_ct, Vt_ct = np.linalg.svd(proj_comp_on_tangent, full_matrices=True)

    print(f"\n  Singular values of complement-on-tangent projection:")
    for i in range(min(16, len(S_ct))):
        marker = " <-- tangent missing dir" if i < rank_comp_tangent else ""
        print(f"    sv[{i:2d}] = {S_ct[i]:.6f}{marker}")

    # The rank_comp_tangent strongest right-singular vectors of Vt_ct
    # are the tangent directions present in the complement
    # (i.e., the directions the 40 points DON'T span)

    for idx in range(min(rank_comp_tangent, 4)):
        tang_dir = Vt_ct[idx]  # 16-dim, in tangent parameter space
        X = tang_dir.reshape(4, 4)
        dP = V1 @ X @ V0.T + V0 @ X.T @ V1.T

        print(f"\n  Missing tangent direction {idx}:")
        dP_eigs = np.sort(np.linalg.eigvalsh(dP))[::-1]
        print(f"    dP eigenvalues: {[round(e,4) for e in dP_eigs]}")
        print(f"    dP trace: {np.trace(dP):.6f}")

        # Key test: is this direction related to the Coxeter element?
        C0 = make_cox(range(8))
        comm_C0 = np.linalg.norm(C0 @ dP - dP @ C0)
        anticomm_C0 = np.linalg.norm(C0 @ dP + dP @ C0)
        print(f"    ||[C0, dP]|| = {comm_C0:.4f}")
        print(f"    ||{{C0, dP}}|| = {anticomm_C0:.4f}")

        # Check against H4 projector
        comm_P0 = np.linalg.norm(P0 @ dP - dP @ P0)
        print(f"    ||[P0, dP]|| = {comm_P0:.4f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print("="*70)

"""
EXP-217b: Conjugate Subgroup Mixing with UNITARY Representation
================================================================
FIX: exp217 used non-unitary polynomial basis for spin-5/2.
     This version unitarizes via the invariant Hermitian form.

KEY COMPUTATION:
1. Unitarize spin-5/2 rep: D_u(g) = H^{1/2} D(g) H^{-1/2}
2. Schur projectors become proper Hermitian projectors (eigenvalues 0,1)
3. Mixing matrix = Tr(P_i(2T_0) @ P_j(2T_g)) / 2 is doubly stochastic
4. Compare to CKM/PMNS
"""

import numpy as np
from itertools import product as iprod, permutations
from scipy.special import comb as binomial

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-217b: UNITARY SPIN-5/2 AND CONJUGATE MIXING")
print("=" * 70)

# === Build 2I ===
def qmul(q1, q2):
    w1,x1,y1,z1 = q1; w2,x2,y2,z2 = q2
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
            w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)

elements = []
for i in range(4):
    for s in [1.0,-1.0]:
        v=[0.0]*4; v[i]=s; elements.append(tuple(v))
for signs in iprod([0.5,-0.5], repeat=4):
    elements.append(tuple(signs))
ep = [(0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
      (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0)]
bv = [0.0, 0.5, PHI/2, 1/(2*PHI)]
for perm in ep:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            for s3 in [1,-1]:
                v=[0.0]*4
                v[perm[1]]=s1*bv[1]; v[perm[2]]=s2*bv[2]; v[perm[3]]=s3*bv[3]
                elements.append(tuple(v))
unique = []
for v in elements:
    if not any(sum((v[k]-u[k])**2 for k in range(4))<1e-10 for u in unique):
        unique.append(v)
elements = unique
N = len(elements)
ea = np.array(elements)

# Mult table
mt = np.zeros((N,N), dtype=int)
for i in range(N):
    for j in range(N):
        p = qmul(elements[i], elements[j])
        mt[i,j] = np.argmin(np.sum((ea-np.array(p))**2, axis=1))
ide = next(i for i in range(N) if all(mt[i,j]==j for j in range(N)))
inv_t = np.zeros(N, dtype=int)
for i in range(N):
    for k in range(N):
        if mt[i,k]==ide: inv_t[i]=k; break

# 2T
sub2T = sorted([i for i,q in enumerate(elements) if
    sum(1 for c in q if abs(c)<1e-10)==3 or
    (sum(1 for c in q if abs(c)<1e-10)==0 and all(abs(abs(c)-0.5)<1e-10 for c in q))])
print("|2I|=%d, |2T|=%d, identity=%d" % (N, len(sub2T), ide))

# === Spin-j representation ===
def quat_to_su2(q):
    w,x,y,z = q
    return np.array([[complex(w,x), complex(y,z)],
                     [-complex(y,-z), complex(w,-x)]])

def spin_j_matrix(U, j):
    d = int(2*j+1)
    a,b,c,d_val = U[0,0],U[0,1],U[1,0],U[1,1]
    D = np.zeros((d,d), dtype=complex)
    for mi in range(d):
        m = j-mi; p=int(j+m); qv=int(j-m)
        for mpi in range(d):
            mp = j-mpi; pp=int(j+mp)
            val = 0
            for s in range(p+1):
                for t in range(qv+1):
                    if s+t==pp:
                        val += (binomial(p,s,exact=True)*binomial(qv,t,exact=True)*
                                a**s * c**(p-s) * b**t * d_val**(qv-t))
            D[mpi,mi] = val
    return D

# Compute spin-5/2 in polynomial basis
j = 2.5; dim = 6
print("\nComputing spin-5/2 representation...")
rep = [spin_j_matrix(quat_to_su2(elements[i]), j) for i in range(N)]

# === UNITARIZE ===
print("Unitarizing via invariant Hermitian form...")
H = np.zeros((dim,dim), dtype=complex)
for i in range(N):
    H += rep[i].T.conj() @ rep[i]
H /= N

# H should be positive definite
eH, vH = np.linalg.eigh(H)
print("H eigenvalues:", eH)
assert all(eH > 0), "H not positive definite!"

Hsqrt = vH @ np.diag(np.sqrt(eH)) @ vH.T.conj()
Hsqrt_inv = vH @ np.diag(1.0/np.sqrt(eH)) @ vH.T.conj()

rep_u = [Hsqrt @ rep[i] @ Hsqrt_inv for i in range(N)]

# Verify unitarity
max_err = max(np.max(np.abs(rep_u[i] @ rep_u[i].T.conj() - np.eye(dim)))
              for i in range(N))
print("Max unitarity error: %.2e" % max_err)

# Verify representation property
sample_err = max(np.max(np.abs(rep_u[i] @ rep_u[j] - rep_u[mt[i,j]]))
                 for i in [0,1,10,50] for j in [0,1,10,50])
print("Representation error: %.2e" % sample_err)

# === 2T conjugacy classes and characters ===
print("\n--- 2T Schur Projectors ---")

# Conjugacy classes of 2T
vis = set()
cl2T = []
for g in sub2T:
    if g in vis: continue
    cc = set()
    for q in sub2T:
        cc.add(mt[mt[q,g], inv_t[q]])
    vis |= cc
    cl2T.append(sorted(cc))
cl2T.sort(key=lambda c: (len(c), -elements[c[0]][0]))

# Identify class types and assign psi_5 characters
# psi_5 (real 2D): chi(e)=2, chi(-1)=-2, chi(order4)=0, chi(order6)=1, chi(order3)=-1
chi5 = np.zeros(N)
for ci, cc in enumerate(cl2T):
    g = cc[0]; ord_g = 1; cur = g
    while cur != ide and ord_g < 200: cur=mt[cur,g]; ord_g += 1
    sz = len(cc)
    if sz==1 and ord_g==1: val = 2.0
    elif sz==1 and ord_g==2: val = -2.0
    elif sz==6: val = 0.0
    elif sz==4 and ord_g==6: val = 1.0
    elif sz==4 and ord_g==3: val = -1.0
    else: val = 0.0
    for g2 in cc: chi5[g2] = val
    print("  Class %d: size %d, order %d, chi_5=%.1f" % (ci, sz, ord_g, val))

# Schur projector for psi_5
P5 = np.zeros((dim,dim), dtype=complex)
for g in sub2T:
    P5 += chi5[g] * rep_u[g]
P5 *= 2.0 / 24.0

# Verify projector properties with UNITARY representation
P5_err = np.max(np.abs(P5 @ P5 - P5))
P5_herm = np.max(np.abs(P5 - P5.T.conj()))
P5_trace = np.trace(P5).real
print("\nP_5: Hermitian err=%.2e, projector err=%.2e, trace=%.4f"
      % (P5_herm, P5_err, P5_trace))

evals_P5 = np.linalg.eigvalsh(P5)
print("P_5 eigenvalues:", np.sort(evals_P5)[::-1])

# Extract psi_5 subspace (eigenvalue 1)
ev5, V5_full = np.linalg.eigh(P5)
V5 = V5_full[:, ev5 > 0.5]
print("psi_5 subspace: %d dimensions" % V5.shape[1])

# Complement = psi_3 + psi_4 (4-dimensional)
V34 = V5_full[:, ev5 < 0.5]
print("psi_3+psi_4 subspace: %d dimensions" % V34.shape[1])

# Separate psi_3 and psi_4 using commutant in 4D complement
print("\nSeparating psi_3 and psi_4...")
np.random.seed(42)
comm_basis = []
for trial in range(30):
    M = np.random.randn(4,4) + 1j*np.random.randn(4,4)
    Mavg = np.zeros_like(M)
    for g in sub2T:
        Dg = V34.T.conj() @ rep_u[g] @ V34
        Mavg += Dg @ M @ np.linalg.inv(Dg)
    Mavg /= len(sub2T)
    comm_basis.append(Mavg.flatten())
_, S_c, Vh_c = np.linalg.svd(np.array(comm_basis))
rank_c = np.sum(S_c > 1e-8)
print("Commutant rank = %d" % rank_c)

# The commutant should be 2D (2 irreps of multiplicity 1)
# Use the first commutant basis vector to separate
M_c = Vh_c[0].reshape(4,4)
ev_c, vec_c = np.linalg.eig(M_c)
print("Commutant eigenvalues:", ["%.6f%+.6fi" % (e.real, e.imag) for e in ev_c])

# Group eigenvalues in pairs
groups = []
used = [False]*4
for i in range(4):
    if used[i]: continue
    grp = [i]; used[i] = True
    for jj in range(i+1, 4):
        if not used[jj] and abs(ev_c[i]-ev_c[jj]) < 1e-6:
            grp.append(jj); used[jj] = True
    groups.append(grp)
print("Groups:", [len(g) for g in groups])

if len(groups) == 2 and all(len(g)==2 for g in groups):
    V3_loc = vec_c[:, groups[0]]
    V4_loc = vec_c[:, groups[1]]
    V3 = V34 @ V3_loc
    V4 = V34 @ V4_loc
    # Orthogonalize within each subspace
    V3, _ = np.linalg.qr(V3, mode='reduced')
    V4, _ = np.linalg.qr(V4, mode='reduced')
    V5, _ = np.linalg.qr(V5, mode='reduced')
    separated = True
else:
    # Alternative: use TWO commutant matrices
    print("Trying combined commutant...")
    M_c2 = Vh_c[1].reshape(4,4) if rank_c >= 2 else M_c
    M_combo = M_c + 0.7123 * M_c2  # break degeneracy
    ev_combo, vec_combo = np.linalg.eig(M_combo)
    print("Combined eigenvalues:", ["%.6f%+.6fi" % (e.real,e.imag) for e in ev_combo])
    groups2 = []
    used2 = [False]*4
    for i in range(4):
        if used2[i]: continue
        grp = [i]; used2[i] = True
        for jj in range(i+1, 4):
            if not used2[jj] and abs(ev_combo[i]-ev_combo[jj]) < 1e-5:
                grp.append(jj); used2[jj] = True
        groups2.append(grp)
    print("Groups2:", [len(g) for g in groups2])
    if len(groups2) == 2 and all(len(g)==2 for g in groups2):
        V3_loc = vec_combo[:, groups2[0]]
        V4_loc = vec_combo[:, groups2[1]]
        V3 = V34 @ V3_loc
        V4 = V34 @ V4_loc
        V3, _ = np.linalg.qr(V3, mode='reduced')
        V4, _ = np.linalg.qr(V4, mode='reduced')
        V5, _ = np.linalg.qr(V5, mode='reduced')
        separated = True
    else:
        print("FAILED to separate psi_3 and psi_4")
        separated = False

if separated:
    # Verify orthogonality
    o35 = np.max(np.abs(V3.T.conj() @ V5))
    o45 = np.max(np.abs(V4.T.conj() @ V5))
    o34 = np.max(np.abs(V3.T.conj() @ V4))
    print("\nOrthogonality: |V3.V5|=%.2e, |V4.V5|=%.2e, |V3.V4|=%.2e" % (o35, o45, o34))

    # Build projectors
    P = [V3 @ V3.T.conj(), V4 @ V4.T.conj(), V5 @ V5.T.conj()]

    # Verify sum = I
    Psum_err = np.max(np.abs(P[0]+P[1]+P[2] - np.eye(dim)))
    print("Sum of projectors = I: err=%.2e" % Psum_err)

    # Verify each is rank-2 Hermitian projector
    for k in range(3):
        print("  P_%d: rank=%d, trace=%.4f, Hermitian err=%.2e, projector err=%.2e"
              % (k, np.linalg.matrix_rank(P[k], tol=1e-8),
                 np.trace(P[k]).real,
                 np.max(np.abs(P[k]-P[k].T.conj())),
                 np.max(np.abs(P[k]@P[k]-P[k]))))

    # Identify which subspace is real (psi_5) vs complex (psi_3/psi_4)
    labels = []
    for k, V in enumerate([V3, V4, V5]):
        traces = []
        for g in sub2T:
            DV = rep_u[g] @ V
            R = np.linalg.lstsq(V, DV, rcond=None)[0]
            traces.append(np.trace(R))
        is_real = all(abs(tr.imag) < 1e-4 for tr in traces)
        labels.append("REAL(psi_5)" if is_real else "COMPLEX(psi_3/4)")
        print("  Subspace %d: %s" % (k, labels[-1]))

    # === CONJUGATE SUBGROUP MIXING ===
    print("\n--- CONJUGATE SUBGROUP MIXING ---")

    # Find coset representatives
    coset_reps = [ide]
    assigned = set(sub2T)
    while len(assigned) < N:
        for g in range(N):
            if g not in assigned:
                coset = set(mt[g, h] for h in sub2T)
                coset_reps.append(g)
                assigned |= coset
                break

    print("Coset reps: %d" % len(coset_reps))

    # Verify conjugate subgroups are distinct
    conj_subs = []
    for g in coset_reps:
        gi = inv_t[g]
        cs = frozenset(mt[mt[g,h], gi] for h in sub2T)
        conj_subs.append(cs)
    n_distinct = len(set(conj_subs))
    print("Distinct conjugate subgroups: %d" % n_distinct)

    # Compute mixing matrix for each coset
    for ci, g in enumerate(coset_reps):
        if ci == 0: continue
        Dg = rep_u[g]

        # Conjugated subspaces
        V_conj = [Dg @ V for V in [V3, V4, V5]]
        P_conj = [V @ V.T.conj() for V in V_conj]

        # 3x3 mixing matrix
        M = np.zeros((3,3))
        for i in range(3):
            for jj in range(3):
                M[i,jj] = np.trace(P[i] @ P_conj[jj]).real / 2.0

        print("\n  Coset %d (g=%d, q=(%.3f,%.3f,%.3f,%.3f)):"
              % (ci, g, elements[g][0], elements[g][1], elements[g][2], elements[g][3]))
        for r in range(3):
            print("    [%s]" % "  ".join("%.6f" % M[r,c] for c in range(3)))
        rs = np.sum(M, axis=1)
        cs = np.sum(M, axis=0)
        print("    Row sums: [%.4f, %.4f, %.4f]" % tuple(rs))
        print("    Col sums: [%.4f, %.4f, %.4f]" % tuple(cs))

        # Compare to known matrices
        CKM2 = np.array([[0.9484,0.0503,0.0000146],[0.0489,0.9744,0.00168],[0.0001,0.0015,0.9984]])
        PMNS2 = np.array([[0.6778,0.3002,0.0220],[0.1608,0.3531,0.4861],[0.1614,0.3467,0.4919]])
        TBM = np.array([[2./3,1./3,0],[1./6,1./3,1./2],[1./6,1./3,1./2]])
        DEM = np.ones((3,3))/3

        best = {}
        for name, ref in [("CKM", CKM2), ("PMNS", PMNS2), ("TBM", TBM), ("DEM", DEM)]:
            best_e = 1e10
            for pr in permutations(range(3)):
                for pc in permutations(range(3)):
                    Mp = M[np.array(pr),:][:,np.array(pc)]
                    e = np.sqrt(np.mean((Mp-ref)**2))
                    if e < best_e: best_e = e
            best[name] = best_e
        print("    RMS errors: CKM=%.4f, PMNS=%.4f, TBM=%.4f, DEM=%.4f"
              % (best["CKM"], best["PMNS"], best["TBM"], best["DEM"]))

        # Look for exact ratios
        print("    Geometric analysis:")
        for r in range(3):
            for c in range(3):
                v = M[r,c]
                if v > 1e-6:
                    # Test simple fractions
                    bf = None; be = 1
                    for aa in range(0,21):
                        for bb in range(1,21):
                            e = abs(v - aa/bb)
                            if e < be: be = e; bf = (aa, bb)
                    # Test phi expressions
                    bp = None; bpe = 1
                    for a_t in range(-5,6):
                        for b_t in range(-5,6):
                            for c_t in range(1,10):
                                test = (a_t + b_t*PHI) / c_t
                                if test > -0.1 and abs(v-test) < bpe:
                                    bpe = abs(v-test)
                                    bp = "%d+%d*phi/%d" % (a_t,b_t,c_t) if b_t!=0 else "%d/%d" % (a_t,c_t)
                    print("      M[%d,%d]=%.8f ~ %d/%d (err=%.1e) ~ %s (err=%.1e)"
                          % (r,c,v, bf[0],bf[1],be, bp, bpe))

    # === ALSO TRY j=2 (dim 5, commutant rank 3) ===
    print("\n--- MIXING FOR j=2 (dim 5) ---")
    j2 = 2.0; d2 = 5
    rep2 = [spin_j_matrix(quat_to_su2(elements[i]), j2) for i in range(N)]
    H2 = sum(rep2[i].T.conj() @ rep2[i] for i in range(N)) / N
    eH2, vH2 = np.linalg.eigh(H2)
    Hs2 = vH2 @ np.diag(np.sqrt(eH2)) @ vH2.T.conj()
    Hs2i = vH2 @ np.diag(1./np.sqrt(eH2)) @ vH2.T.conj()
    rep2_u = [Hs2 @ rep2[i] @ Hs2i for i in range(N)]

    # Schur projector for psi_5 in j=2
    P5_j2 = np.zeros((d2,d2), dtype=complex)
    for g in sub2T:
        P5_j2 += chi5[g] * rep2_u[g]
    P5_j2 *= 2.0/24.0
    ev5_j2 = np.linalg.eigvalsh(P5_j2)
    print("P_5(j=2) eigenvalues:", np.sort(ev5_j2)[::-1])
    print("P_5(j=2) trace=%.4f, projector err=%.2e"
          % (np.trace(P5_j2).real, np.max(np.abs(P5_j2@P5_j2-P5_j2))))

    # From exp215: j=2 (dim 5) = rho_7, branching: psi_1 + psi_2 + psi_6 (1+1+3)
    # No 2D irreps! So no "generation" mixing in dim 5.
    # The commutant rank=3 from exp216 was for COMMUTANT of 2T, which has
    # 3 distinct irreps (psi_1, psi_2, psi_6) with dims 1, 1, 3
    print("  dim 5 = psi_1(1D) + psi_2(1D) + psi_6(3D) -> no 2D generation reps")
    print("  Generation mixing only in dim 6 (j=5/2)")

    # === ALSO: PAIR-WISE MIXING BETWEEN ALL 5 CONJUGATES ===
    print("\n--- PAIRWISE MIXING BETWEEN ALL 5 CONJUGATES ---")
    # Compute mixing for all pairs (ci, cj)
    mixing_matrices = {}
    for ci, gi in enumerate(coset_reps):
        Dgi = rep_u[gi]
        P_ci = [Dgi @ V @ (Dgi @ V).T.conj() for V in [V3, V4, V5]]
        # Actually we need orthogonalized versions
        V_ci = [Dgi @ V for V in [V3, V4, V5]]
        P_ci = [v @ v.T.conj() for v in V_ci]
        for cj, gj in enumerate(coset_reps):
            if cj <= ci: continue
            Dgj = rep_u[gj]
            V_cj = [Dgj @ V for V in [V3, V4, V5]]
            P_cj = [v @ v.T.conj() for v in V_cj]
            M_pair = np.zeros((3,3))
            for ii in range(3):
                for jj in range(3):
                    M_pair[ii,jj] = np.trace(P_ci[ii] @ P_cj[jj]).real / 2.0
            mixing_matrices[(ci,cj)] = M_pair
            rs = np.sum(M_pair, axis=1)
            print("  (%d,%d): diag=[%.4f,%.4f,%.4f] rows=[%.4f,%.4f,%.4f]"
                  % (ci,cj, M_pair[0,0],M_pair[1,1],M_pair[2,2], rs[0],rs[1],rs[2]))

    # Are all pairwise matrices the same?
    vals = list(mixing_matrices.values())
    all_same = all(np.max(np.abs(v - vals[0])) < 1e-10 for v in vals)
    print("All pairwise matrices identical: %s" % all_same)

    # Check if there's a pattern depending on "distance" between cosets
    if not all_same:
        print("\nDistinct pairwise mixing matrices:")
        seen = []
        for (ci,cj), M_pair in mixing_matrices.items():
            is_new = True
            for M_prev in seen:
                if np.max(np.abs(M_pair - M_prev)) < 1e-6:
                    is_new = False; break
                # Also check permutations
                for pr in permutations(range(3)):
                    for pc in permutations(range(3)):
                        if np.max(np.abs(M_pair[np.array(pr),:][:,np.array(pc)] - M_prev)) < 1e-6:
                            is_new = False; break
                    if not is_new: break
            if is_new:
                seen.append(M_pair)
                print("  (%d,%d):" % (ci,cj))
                for r in range(3):
                    print("    [%s]" % "  ".join("%.6f" % M_pair[r,c] for c in range(3)))
        print("Total distinct: %d" % len(seen))

print("\n" + "=" * 70)
print("EXP-217b COMPLETE")
print("=" * 70)

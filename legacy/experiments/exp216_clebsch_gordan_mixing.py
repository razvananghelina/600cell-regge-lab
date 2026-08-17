"""
EXP-216: Clebsch-Gordan Decomposition and Mixing Matrix
========================================================
PURPOSE: Compute the EXPLICIT Clebsch-Gordan decomposition of the dim-6
         irrep (rho_8) of 2I restricted to 2T. Extract the 3x3 mixing
         matrix from the three 2D subspaces within the 6D space.

KEY IDEA:
- From exp215: rho_8 (dim 6) = psi_3 + psi_4 + psi_5 under 2T
- 2I is a subgroup of SU(2). Irreps of 2I = restrictions of spin-j reps.
- dim 6 => j = 5/2 (spin-5/2 representation)
- Compute explicit 6x6 matrices for each 2I element via spin-5/2
- Restrict to 2T, find 3 two-dim invariant subspaces
- The change-of-basis matrix IS the Clebsch-Gordan mixing matrix
- Compare 3x3 structure to CKM/PMNS

METHOD:
1. Quaternion -> SU(2) matrix for each element
2. SU(2) -> spin-j matrix via polynomial action
3. Compute all irreps (j=0,1/2,1,...,5/2) and verify character match
4. For j=5/2: restrict to 2T, find invariant 2D subspaces
5. Extract mixing matrix from subspace orientations
"""

import numpy as np
from itertools import product as iprod
from scipy.special import comb as binomial

PHI = (1 + np.sqrt(5)) / 2
PSI = 1 - PHI
ALPHA = 7.2973525693e-3

print("=" * 70)
print("EXP-216: CLEBSCH-GORDAN DECOMPOSITION AND MIXING MATRIX")
print("=" * 70)

# ======================================================================
# PART 1: Construct 2I and 2T
# ======================================================================
print("\n--- PART 1: Groups 2I and 2T ---")

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
# Deduplicate
unique = []
ea = np.array(elements)
for v in elements:
    is_dup = False
    for u in unique:
        if sum((v[k]-u[k])**2 for k in range(4)) < 1e-10:
            is_dup = True; break
    if not is_dup:
        unique.append(v)
elements = unique
N = len(elements)
print("|2I| = %d" % N)

# Identify 2T
subgroup_2T = []
for i, q in enumerate(elements):
    n_zero = sum(1 for c in q if abs(c) < 1e-10)
    if n_zero == 3 or (n_zero == 0 and all(abs(abs(c)-0.5) < 1e-10 for c in q)):
        subgroup_2T.append(i)
print("|2T| = %d" % len(subgroup_2T))
is_2T = set(subgroup_2T)

# ======================================================================
# PART 2: Quaternion -> SU(2) -> Spin-j matrices
# ======================================================================
print("\n--- PART 2: Spin-j Representations ---")

def quat_to_su2(q):
    """Convert quaternion (w,x,y,z) to SU(2) matrix [[a,b],[-b*,a*]]"""
    w, x, y, z = q
    a = complex(w, x)
    b = complex(-y, -z)
    # Actually: q = w + xi + yj + zk -> U = [[w+xi, -y-zi],[y-zi, w-xi]]
    # So a = w+xi, b = -y-zi
    # Wait let me be precise. The standard map is:
    # q = a + bi + cj + dk maps to:
    #   [[a+bi, c+di],[-c+di, a-bi]]
    # Verify: det = (a+bi)(a-bi) - (c+di)(-c+di) = a^2+b^2+c^2+d^2 = 1
    a_val = complex(w, x)   # a + bi
    b_val = complex(y, z)   # c + di
    U = np.array([[a_val, b_val], [-np.conj(b_val), np.conj(a_val)]])
    return U

def spin_j_matrix(U, j):
    """Compute the spin-j representation matrix of SU(2) element U.
    Uses the polynomial action: z1 -> a*z1 + c*z2, z2 -> b*z1 + d*z2
    Basis: e_m = z1^(j+m) * z2^(j-m) / sqrt(C(2j, j-m)) for m = j, j-1, ..., -j
    """
    two_j = int(2*j)
    dim = two_j + 1
    a, b = U[0, 0], U[0, 1]
    c, d = U[1, 0], U[1, 1]

    D = np.zeros((dim, dim), dtype=complex)

    for m_idx in range(dim):  # row index, m = j - m_idx
        m = j - m_idx
        p = int(j + m)   # power of z1 in input basis
        q_val = int(j - m)   # power of z2 in input basis

        for mp_idx in range(dim):  # column index, m' = j - mp_idx
            mp = j - mp_idx
            pp = int(j + mp)   # power of z1 in output
            qp = int(j - mp)   # power of z2 in output

            # Coefficient of z1^pp * z2^qp in (a*z1 + c*z2)^p * (b*z1 + d*z2)^q
            val = 0
            for s in range(p + 1):
                for t in range(q_val + 1):
                    if s + t == pp:  # z1 power matches
                        # (p choose s) * a^s * c^(p-s) * (q choose t) * b^t * d^(q-t)
                        coeff = (binomial(p, s, exact=True) *
                                binomial(q_val, t, exact=True) *
                                a**s * c**(p-s) * b**t * d**(q_val-t))
                        val += coeff
            D[mp_idx, m_idx] = val

    return D

# Test: spin-0 should be trivial
U_test = quat_to_su2(elements[0])  # identity
D0 = spin_j_matrix(U_test, 0)
print("Spin-0 at identity: %s" % D0)

# Compute all spin-j representations for j = 0, 1/2, 1, 3/2, 2, 5/2
# and verify which ones give distinct irreps of 2I
j_values = [0, 0.5, 1, 1.5, 2, 2.5]
print("\nSpin-j representation analysis:")

# For each j, compute trace (= character) for representative of each class
# Use a few test elements to identify which j gives which irrep
test_elems = [0]  # identity
# Find element with w ~ 0.809 (nearest neighbor, order 10)
for i in range(N):
    if abs(elements[i][0] - 0.809) < 0.01:
        test_elems.append(i)
        break
# Find element with w ~ 0.5 (type B, order 6)
for i in range(N):
    if abs(elements[i][0] - 0.5) < 0.01 and all(abs(abs(c)-0.5) < 0.01 for c in elements[i]):
        test_elems.append(i)
        break

# Known eigenvalues: {12, 6phi, 6-6phi, 4phi, 4-4phi, 3, -3, 0, -2}
# With dims: {1, 2, 2, 3, 3, 4, 4, 5, 6}
# Spin-j gives dim 2j+1: j=0->1, j=1/2->2, j=1->3, j=3/2->4, j=2->5, j=5/2->6

for j in j_values:
    dim = int(2*j + 1)
    print("\n  j = %.1f (dim %d):" % (j, dim))

    # Compute character (trace) at all elements
    chars = np.zeros(N, dtype=complex)
    for i in range(N):
        U = quat_to_su2(elements[i])
        D = spin_j_matrix(U, j)
        chars[i] = np.trace(D)

    # Verify character at identity = dim
    print("    chi(e) = %.4f (should be %d)" % (chars[0].real, dim))

    # Character at nearest neighbor
    if len(test_elems) > 1:
        chi_nn = chars[test_elems[1]]
        # Expected: chi(g) = lambda * dim / 12
        print("    chi(nearest neighbor) = %.4f + %.4fi"
              % (chi_nn.real, chi_nn.imag))

    # Compute 600-cell adjacency eigenvalue from characters
    # lambda = (1/dim) * sum_{g in S} chi(g) where S = 12 nearest neighbors
    # S = elements at distance^2 = 2 - 2*cos(pi/5) from identity
    cos_pi5 = np.cos(np.pi/5)
    threshold = 2 - 2*cos_pi5 + 0.01
    nn_sum = 0
    nn_count = 0
    for i in range(N):
        d2 = sum((elements[i][k] - elements[0][k])**2 for k in range(4))
        if 0.01 < d2 < threshold:
            nn_sum += chars[i]
            nn_count += 1
    if nn_count > 0:
        lam = nn_sum / dim
        # Try a+b*phi decomposition
        best = None
        for a_try in range(-10, 13):
            for b_try in range(-10, 11):
                if abs(a_try + b_try*PHI - lam.real) < 0.01 and abs(lam.imag) < 0.01:
                    best = (a_try, b_try)
                    break
            if best: break
        label = "%d + %d*phi" % best if best else "%.4f" % lam.real
        print("    lambda = sum(chi_S)/dim = %.4f = %s (from %d neighbors)"
              % (lam.real, label, nn_count))

# ======================================================================
# PART 3: Focus on j=5/2 (dim 6) - THE generation representation
# ======================================================================
print("\n--- PART 3: j=5/2 (dim 6) Detailed Analysis ---")

j_gen = 2.5  # spin-5/2
dim_gen = 6

# Compute all 6x6 representation matrices
print("Computing 120 matrices of dim 6x6...")
rep_matrices = []
for i in range(N):
    U = quat_to_su2(elements[i])
    D = spin_j_matrix(U, j_gen)
    rep_matrices.append(D)

# Verify: these should form a representation
# Check D(g)*D(h) = D(g*h) for a sample
elem_arr = np.array(elements)
print("Representation check (sample)...")
errors = []
for trial in range(50):
    i = np.random.randint(N)
    j_idx = np.random.randint(N)
    # Find product g_i * g_j
    prod = qmul(elements[i], elements[j_idx])
    prod_arr = np.array(prod)
    dists = np.sum((elem_arr - prod_arr)**2, axis=1)
    k = np.argmin(dists)
    # Check D(i) @ D(j) = D(k)
    err = np.max(np.abs(rep_matrices[i] @ rep_matrices[j_idx] - rep_matrices[k]))
    errors.append(err)
print("  Max representation error: %.2e" % max(errors))

# ======================================================================
# PART 4: Restrict to 2T and find invariant subspaces
# ======================================================================
print("\n--- PART 4: 2T Invariant Subspaces ---")

# For the 2T subgroup, find the common invariant subspaces of all D(g), g in 2T
# Method: compute the projection operators for each irrep of 2T
# P_j = (dim_j / |2T|) * sum_{g in 2T} conj(chi_j(g)) * D(g)

# First, compute the sum of all 2T representation matrices
D_sum = np.zeros((dim_gen, dim_gen), dtype=complex)
for i in subgroup_2T:
    D_sum += rep_matrices[i]
D_avg = D_sum / len(subgroup_2T)  # = projection onto trivial component

# Check: D_avg should be rank 0 (no trivial component in spin-5/2 restricted to 2T)
# From branching: rho_8 = psi_3 + psi_4 + psi_5 (all 2D, no trivial)
print("D_avg eigenvalues:", np.sort(np.abs(np.linalg.eigvals(D_avg)))[::-1])

# Method 2: Use the commutant algebra
# The space of matrices commuting with all D(g), g in 2T, has dimension = sum(m_j^2)
# where m_j are the multiplicities of 2T irreps in the decomposition
# From branching: m_3 = m_4 = m_5 = 1, so commutant dim = 1+1+1 = 3

# Compute commutant: C = (1/|2T|) sum_g D(g) (x) conj(D(g))
# Or equivalently, find the space of 6x6 matrices M such that D(g) M D(g)^{-1} = M for all g in 2T

# Schur's lemma approach: average the projection
# For a random matrix M, the average (1/|2T|) sum_g D(g) M D(g)^{-1} projects onto commutant
print("\nFinding commutant algebra of 2T in dim-6 rep...")
np.random.seed(42)

# Use multiple random matrices to span the commutant
commutant_basis = []
for trial in range(20):
    M_rand = np.random.randn(dim_gen, dim_gen) + 1j * np.random.randn(dim_gen, dim_gen)
    M_avg = np.zeros_like(M_rand)
    for i in subgroup_2T:
        D = rep_matrices[i]
        D_inv = np.linalg.inv(D)
        M_avg += D @ M_rand @ D_inv
    M_avg /= len(subgroup_2T)
    commutant_basis.append(M_avg.flatten())

# Find rank of commutant
C_matrix = np.array(commutant_basis)
U_c, S_c, Vh_c = np.linalg.svd(C_matrix)
rank = np.sum(S_c > 1e-8)
print("Commutant rank = %d (expected 3 for three multiplicity-1 irreps)" % rank)

# Extract 3 independent commutant matrices
comm_matrices = []
for k in range(rank):
    M = Vh_c[k].reshape(dim_gen, dim_gen)
    comm_matrices.append(M)
    evals_M = np.linalg.eigvals(M)
    print("  Commutant basis %d eigenvalues: %s"
          % (k, ", ".join("%.3f+%.3fi" % (e.real, e.imag) for e in sorted(evals_M, key=lambda x: -abs(x)))))

# The commutant matrices should be simultaneously diagonalizable
# with eigenvalues grouped in pairs (since each irrep is 2D)
# Diagonalize the first commutant matrix
if rank >= 1:
    evals_c, evecs_c = np.linalg.eig(comm_matrices[0])
    print("\nFirst commutant matrix eigenvalues:")
    for e in sorted(evals_c, key=lambda x: -abs(x)):
        print("  %.6f + %.6fi" % (e.real, e.imag))

    # Group eigenvalues in pairs (should get 3 pairs for three 2D irreps)
    # Each pair of identical eigenvalues spans a 2D invariant subspace
    groups = []
    used = [False] * dim_gen
    for i in range(dim_gen):
        if used[i]:
            continue
        group = [i]
        used[i] = True
        for j in range(i+1, dim_gen):
            if not used[j] and abs(evals_c[i] - evals_c[j]) < 1e-6:
                group.append(j)
                used[j] = True
        groups.append(group)

    print("\nEigenvalue groups (= 2D invariant subspaces):")
    for g_idx, group in enumerate(groups):
        print("  Group %d: indices %s, eigenvalue = %.6f + %.6fi"
              % (g_idx, group, evals_c[group[0]].real, evals_c[group[0]].imag))

    # If we have exactly 3 groups of size 2, extract the 2D subspaces
    if len(groups) == 3 and all(len(g) == 2 for g in groups):
        print("\n  *** Found 3 two-dimensional invariant subspaces! ***")

        # The eigenvectors of the commutant matrix span the 2D subspaces
        subspaces = []
        for g_idx, group in enumerate(groups):
            V = evecs_c[:, group]  # 6x2 matrix spanning the 2D subspace
            subspaces.append(V)
            print("\n  Subspace %d (eigenvalue %.4f+%.4fi):" % (g_idx, evals_c[group[0]].real, evals_c[group[0]].imag))
            print("    Basis vectors (columns):")
            for row in range(dim_gen):
                print("    [%.4f+%.4fi  %.4f+%.4fi]"
                      % (V[row,0].real, V[row,0].imag, V[row,1].real, V[row,1].imag))

        # Verify: each subspace should be 2T-invariant
        print("\n  Invariance check:")
        for g_idx, V in enumerate(subspaces):
            max_err = 0
            for i in subgroup_2T:
                DV = rep_matrices[i] @ V  # Should be V @ R for some 2x2 matrix R
                # Check if DV is in the column span of V
                R = np.linalg.lstsq(V, DV, rcond=None)[0]
                err = np.max(np.abs(DV - V @ R))
                if err > max_err:
                    max_err = err
            print("    Subspace %d: max invariance error = %.2e %s"
                  % (g_idx, max_err, "OK" if max_err < 1e-8 else "FAIL"))

        # ======================================================================
        # PART 5: Identify which 2T irrep each subspace carries
        # ======================================================================
        print("\n--- PART 5: Irrep Identification ---")

        # For each subspace, compute the character of the 2T representation
        for g_idx, V in enumerate(subspaces):
            print("\n  Subspace %d:" % g_idx)
            traces = []
            for i in subgroup_2T:
                DV = rep_matrices[i] @ V
                R = np.linalg.lstsq(V, DV, rcond=None)[0]  # 2x2 rep matrix
                tr = np.trace(R)
                traces.append(tr)

            # Compute character on each 2T conjugacy class
            # Group 2T elements by conjugacy class
            # (using the fact that conjugate elements have same trace)
            class_traces = {}
            for idx_loc, i_glob in enumerate(subgroup_2T):
                tr = traces[idx_loc]
                found = False
                for key in class_traces:
                    if abs(tr - key) < 1e-4:
                        class_traces[key].append(i_glob)
                        found = True
                        break
                if not found:
                    class_traces[tr] = [i_glob]

            print("    Character values (grouped by class):")
            for tr_val, members in sorted(class_traces.items(), key=lambda x: -x[1].__len__()):
                # Representative element
                rep_elem = elements[members[0]]
                print("      chi = %7.4f + %7.4fi  (class size %d, rep w=%.3f)"
                      % (tr_val.real, tr_val.imag, len(members), rep_elem[0]))

            # Check: is this the real 2D rep (psi_5) or a complex one (psi_3/4)?
            all_real = all(abs(tr.imag) < 1e-4 for tr in traces)
            print("    All traces real: %s -> %s"
                  % (all_real, "psi_5 (REAL 2D)" if all_real else "psi_3 or psi_4 (COMPLEX 2D)"))

        # ======================================================================
        # PART 6: The 3x3 Mixing Matrix
        # ======================================================================
        print("\n--- PART 6: The 3x3 Mixing Matrix ---")

        # The full change-of-basis matrix from spin-5/2 basis to generation basis
        # is the 6x6 matrix [V_0 | V_1 | V_2]
        W = np.hstack(subspaces)  # 6x6 unitary matrix
        print("Change-of-basis matrix W (6x6):")
        print("  Columns: [gen_1a, gen_1b | gen_2a, gen_2b | gen_3a, gen_3b]")

        # Check unitarity
        WW = W.T.conj() @ W
        print("  W^dag W diag:", np.abs(np.diag(WW)))
        print("  Max off-diag:", np.max(np.abs(WW - np.diag(np.diag(WW)))))

        # The mixing information is in the OVERLAPS between subspaces
        # Project the spin-5/2 basis vectors onto each generation subspace
        print("\n  Projection of spin-5/2 basis |m> onto generation subspaces:")
        print("  (probability = |<m|gen_k>|^2, summed over 2 basis vectors per gen)")
        print("  %8s" % "m", end="")
        for g_idx in range(3):
            print("  Gen_%d  " % (g_idx+1), end="")
        print("  Total")

        m_values = [j_gen - k for k in range(dim_gen)]  # m = 5/2, 3/2, 1/2, -1/2, -3/2, -5/2
        projection_matrix = np.zeros((dim_gen, 3))

        for m_idx in range(dim_gen):
            m = m_values[m_idx]
            e_m = np.zeros(dim_gen)
            e_m[m_idx] = 1.0
            probs = []
            for g_idx, V in enumerate(subspaces):
                # Probability of |m> being in subspace g_idx
                proj = V.T.conj() @ e_m  # 2-vector of overlaps
                prob = np.sum(np.abs(proj)**2)
                probs.append(prob)
                projection_matrix[m_idx, g_idx] = prob.real
            total = sum(probs)
            print("  %6.1f:" % m, end="")
            for p in probs:
                print("  %.4f " % p.real, end="")
            print("  %.4f" % total.real)

        # The 3x3 matrix we want: how each PAIR of m-values maps to generations
        # Pair (m, -m) together (particle-antiparticle)
        print("\n  Paired (m, -m) projection sums:")
        n_pairs = dim_gen // 2
        for pair_idx in range(n_pairs):
            m_pos = m_values[pair_idx]
            m_neg = m_values[dim_gen - 1 - pair_idx]
            probs_pair = projection_matrix[pair_idx, :] + projection_matrix[dim_gen-1-pair_idx, :]
            print("  (m=%.1f, %.1f): [%s]"
                  % (m_pos, m_neg, "  ".join("%.4f" % p for p in probs_pair)))

        # Also compute the overlap matrix between j=5/2 and j=3/2 (dim 4 = branch)
        # This tests the relation between branch node and generation structure
        print("\n--- PART 7: Comparison with dim-4 (j=3/2, branch node) ---")

        j_branch = 1.5
        dim_branch = 4

        # Compute commutant for j=3/2 restricted to 2T
        rep_branch = []
        for i in range(N):
            U = quat_to_su2(elements[i])
            D = spin_j_matrix(U, j_branch)
            rep_branch.append(D)

        comm_basis_br = []
        for trial in range(20):
            M_rand = np.random.randn(dim_branch, dim_branch) + 1j * np.random.randn(dim_branch, dim_branch)
            M_avg = np.zeros_like(M_rand)
            for i in subgroup_2T:
                D = rep_branch[i]
                D_inv = np.linalg.inv(D)
                M_avg += D @ M_rand @ D_inv
            M_avg /= len(subgroup_2T)
            comm_basis_br.append(M_avg.flatten())

        C_br = np.array(comm_basis_br)
        _, S_br, Vh_br = np.linalg.svd(C_br)
        rank_br = np.sum(S_br > 1e-8)
        print("dim-4 commutant rank = %d (expected 2 for psi_3 + psi_4)" % rank_br)

        if rank_br >= 1:
            M_br = Vh_br[0].reshape(dim_branch, dim_branch)
            evals_br, evecs_br = np.linalg.eig(M_br)
            print("dim-4 commutant eigenvalues:", ["%.4f+%.4fi" % (e.real, e.imag) for e in evals_br])

            groups_br = []
            used_br = [False] * dim_branch
            for i in range(dim_branch):
                if used_br[i]: continue
                group = [i]; used_br[i] = True
                for j_v in range(i+1, dim_branch):
                    if not used_br[j_v] and abs(evals_br[i] - evals_br[j_v]) < 1e-6:
                        group.append(j_v); used_br[j_v] = True
                groups_br.append(group)

            print("Groups: %s" % groups_br)
            if len(groups_br) == 2:
                print("  Two 2D subspaces in dim-4 (as expected: psi_3 + psi_4)")

                # Check if each is complex 2D
                for g_idx, group in enumerate(groups_br):
                    V_br = evecs_br[:, group]
                    traces_br = []
                    for i in subgroup_2T:
                        DV = rep_branch[i] @ V_br
                        R = np.linalg.lstsq(V_br, DV, rcond=None)[0]
                        traces_br.append(np.trace(R))
                    all_real_br = all(abs(tr.imag) < 1e-4 for tr in traces_br)
                    print("    Subspace %d: all traces real = %s" % (g_idx, all_real_br))

        # ======================================================================
        # PART 8: CKM/PMNS Comparison
        # ======================================================================
        print("\n--- PART 8: CKM/PMNS Comparison ---")

        # The projection matrix gives the probability of each spin component
        # being in each generation subspace. Compare to CKM^2 and PMNS^2.

        print("\n6x3 Projection matrix (|<m|gen>|^2):")
        print("  m\\gen     Gen_1   Gen_2   Gen_3")
        for m_idx in range(dim_gen):
            print("  %5.1f:  " % m_values[m_idx], end="")
            for g in range(3):
                print("%.4f  " % projection_matrix[m_idx, g], end="")
            print()

        # Extract the 3x3 by pairing m values
        P33 = np.zeros((3, 3))
        for pair_idx in range(3):
            P33[pair_idx, :] = (projection_matrix[pair_idx, :] +
                               projection_matrix[5-pair_idx, :]) / 2.0

        print("\n3x3 matrix (averaged over m-pairs):")
        for r in range(3):
            print("  [%s]" % "  ".join("%.4f" % P33[r,c] for c in range(3)))

        # CKM |V|^2
        CKM_sq = np.array([
            [0.9484, 0.0503, 0.0000146],
            [0.0489, 0.9744, 0.00168],
            [0.0000640, 0.00150, 1.026]
        ])
        # PMNS |U|^2
        s12_sq, s23_sq, s13_sq = 0.307, 0.546, 0.0220
        c12_sq, c23_sq, c13_sq = 1-s12_sq, 1-s23_sq, 1-s13_sq
        PMNS_sq = np.array([
            [c12_sq*c13_sq, s12_sq*c13_sq, s13_sq],
            [0.327, 0.353, 0.320],
            [0.174, 0.340, 0.487]
        ])

        print("\nCKM |V|^2:")
        for r in range(3):
            print("  [%s]" % "  ".join("%.4f" % CKM_sq[r,c] for c in range(3)))

        print("\nPMNS |U|^2:")
        for r in range(3):
            print("  [%s]" % "  ".join("%.4f" % PMNS_sq[r,c] for c in range(3)))

        # Check if P33 is related by any permutation
        print("\nPermutation tests (P33 vs CKM/PMNS):")
        from itertools import permutations
        best_ckm_err = 1e10
        best_pmns_err = 1e10
        best_ckm_perm = None
        best_pmns_perm = None
        for perm_r in permutations(range(3)):
            for perm_c in permutations(range(3)):
                P_perm = P33[np.array(perm_r), :][:, np.array(perm_c)]
                err_ckm = np.max(np.abs(P_perm - CKM_sq))
                err_pmns = np.max(np.abs(P_perm - PMNS_sq))
                if err_ckm < best_ckm_err:
                    best_ckm_err = err_ckm
                    best_ckm_perm = (perm_r, perm_c)
                if err_pmns < best_pmns_err:
                    best_pmns_err = err_pmns
                    best_pmns_perm = (perm_r, perm_c)

        print("  Best CKM match: err=%.4f (perm_r=%s, perm_c=%s)"
              % (best_ckm_err, best_ckm_perm[0], best_ckm_perm[1]))
        print("  Best PMNS match: err=%.4f (perm_r=%s, perm_c=%s)"
              % (best_pmns_err, best_pmns_perm[0], best_pmns_perm[1]))

        # Democratic matrix for comparison
        dem_err = np.max(np.abs(P33 - 1.0/3.0))
        print("  Distance from democratic (1/3 everywhere): %.4f" % dem_err)

    else:
        print("  WARNING: Expected 3 groups of size 2, got %d groups: %s"
              % (len(groups), [len(g) for g in groups]))
        print("  Trying alternative decomposition...")

        # If commutant didn't separate properly, try using TWO commutant matrices
        if rank >= 2:
            # Simultaneously diagonalize two commutant matrices
            M1 = comm_matrices[0]
            M2 = comm_matrices[1]
            # Linear combination to break degeneracy
            M_combo = M1 + 0.3713 * M2  # irrational coefficient
            evals_combo2, evecs_combo2 = np.linalg.eig(M_combo)
            print("\n  Combined commutant eigenvalues:")
            for e in sorted(evals_combo2, key=lambda x: -abs(x)):
                print("    %.6f + %.6fi" % (e.real, e.imag))

            groups2 = []
            used2 = [False] * dim_gen
            for i in range(dim_gen):
                if used2[i]: continue
                group = [i]; used2[i] = True
                for j_v in range(i+1, dim_gen):
                    if not used2[j_v] and abs(evals_combo2[i] - evals_combo2[j_v]) < 1e-5:
                        group.append(j_v); used2[j_v] = True
                groups2.append(group)
            print("  Groups from combined: %s" % [len(g) for g in groups2])


else:
    print("  Commutant rank not 3. Branching structure may be different.")
    print("  Check: is rho_{5/2} irreducible for 2I?")

# ======================================================================
# PART 9: Additional spin-j analysis
# ======================================================================
print("\n--- PART 9: All Spin-j Branching Summary ---")

for j in [0, 0.5, 1, 1.5, 2, 2.5]:
    dim_j = int(2*j + 1)
    rep_j = []
    for i in range(N):
        U = quat_to_su2(elements[i])
        D = spin_j_matrix(U, j)
        rep_j.append(D)

    # Commutant rank for 2T
    comm_basis_j = []
    for trial in range(max(20, dim_j**2)):
        M_rand = np.random.randn(dim_j, dim_j) + 1j * np.random.randn(dim_j, dim_j)
        M_avg = np.zeros_like(M_rand)
        for i in subgroup_2T:
            D = rep_j[i]
            D_inv = np.linalg.inv(D)
            M_avg += D @ M_rand @ D_inv
        M_avg /= len(subgroup_2T)
        comm_basis_j.append(M_avg.flatten())

    C_j = np.array(comm_basis_j)
    _, S_j, _ = np.linalg.svd(C_j)
    rank_j = np.sum(S_j > 1e-8)

    # Character at identity element
    chi_e = dim_j
    # Character at nearest neighbor
    U_nn = quat_to_su2(elements[1] if abs(elements[1][0] - 0.809) < 0.01 else elements[2])
    D_nn = spin_j_matrix(U_nn, j)
    chi_nn = np.trace(D_nn).real

    print("  j=%.1f (dim=%d): commutant rank=%d, chi(nn)=%.4f"
          % (j, dim_j, rank_j, chi_nn))

print("\n" + "=" * 70)
print("EXP-216 COMPLETE")
print("=" * 70)

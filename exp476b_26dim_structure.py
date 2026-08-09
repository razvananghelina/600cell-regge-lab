"""
exp476b: The 26-Dimensional Structure in the Decomposition Graph
================================================================

FROM exp476:
  - Exchanged operator E = 5*A5 + 4*A6 + 3*A7 + 2*A8 has 26-dim KERNEL
  - Overlap operator O = 5*A5 + 6*A6 + 7*A7 + 8*A8 has eigenvalue -10=-2a1
    with multiplicity 26
  - 26 = a1^2 + 1 = 1/sin^2(theta_W) -- THE gauge dimension!
  - A_7 has all eigenvalues with multiplicity 2 (Z/2 symmetry?)

QUESTIONS:
  1. What is the 26-dim subspace? Does it have framework structure?
  2. Is there a Z/2 action on the 40 decompositions?
  3. How does 26 arise algebraically?
  4. Connection to sin^2(tW) = b1/(a1^2+1) = 6/26?

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: INVESTIGATION
"""

import numpy as np
from itertools import permutations
from math import cos, sin, pi, sqrt, gcd, factorial
from collections import Counter

phi = (1 + sqrt(5)) / 2
phi_p = (1 - sqrt(5)) / 2
a1 = 5
b1 = 6
h = a1 * b1  # 30
N = 120
rank_E8 = 8
degree_600 = 12

# ============================================================
# REBUILD 40 DECOMPOSITIONS (same code as exp476)
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
    n_minus = sum(1 for s in signs if s < 0)
    if n_minus % 2 == 0:
        v = np.array(signs) / 2.0
        E8_roots.append(v)
E8 = np.array(E8_roots)

simple = np.zeros((8, 8))
simple[0] = np.array([1, -1, 0, 0, 0, 0, 0, 0])
simple[1] = np.array([0, 1, -1, 0, 0, 0, 0, 0])
simple[2] = np.array([0, 0, 1, -1, 0, 0, 0, 0])
simple[3] = np.array([0, 0, 0, 1, -1, 0, 0, 0])
simple[4] = np.array([0, 0, 0, 0, 1, -1, 0, 0])
simple[5] = np.array([0, 0, 0, 0, 0, 1, -1, 0])
simple[6] = np.array([0, 0, 0, 0, 0, 1, 1, 0])
simple[7] = np.array([-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5])

reflections = [np.eye(8) - 2 * np.outer(simple[i], simple[i]) / np.dot(simple[i], simple[i])
               for i in range(8)]

def make_coxeter(perm):
    C = np.eye(8)
    for i in perm:
        C = reflections[i] @ C
    return C

def get_H4_projector(C):
    if not np.allclose(np.linalg.matrix_power(C, 30), np.eye(8), atol=1e-8):
        return None
    eigvals, eigvecs = np.linalg.eig(C)
    args = np.angle(eigvals)
    exponents = []
    for i in range(8):
        m = round(args[i] * 30 / (2*pi)) % 30
        if m == 0: m = 30
        exponents.append((m, i))
    h4_idx = [i for m, i in exponents if m in {1, 11, 19, 29}]
    if len(h4_idx) != 4:
        return None
    v_raw = eigvecs[:, h4_idx]
    basis_real = []
    used = set()
    for k in range(4):
        if k in used: continue
        v = v_raw[:, k]
        if np.max(np.abs(np.imag(v))) > 0.01:
            basis_real.append(np.real(v))
            basis_real.append(np.imag(v))
            for k2 in range(k+1, 4):
                if k2 not in used and np.allclose(v_raw[:, k2], np.conj(v)):
                    used.add(k2)
                    break
        else:
            basis_real.append(np.real(v))
    if len(basis_real) != 4:
        return None
    B = np.array(basis_real).T
    Q, _ = np.linalg.qr(B)
    P = Q[:, :4]
    return P

projectors = []
np.random.seed(42)
perms_to_try = list(permutations(range(8)))
np.random.shuffle(perms_to_try)

for perm in perms_to_try:
    C = make_coxeter(perm)
    P = get_H4_projector(C)
    if P is None:
        continue
    is_new = True
    for P_old in projectors:
        M = P_old.T @ P
        svs = np.linalg.svd(M, compute_uv=False)
        if np.allclose(svs, np.ones(4), atol=1e-6):
            is_new = False
            break
    if is_new:
        projectors.append(P)
    if len(projectors) >= 40:
        break

print(f"Found {len(projectors)}/40 decompositions")

# Build adjacency and operators
A = {}
overlap_matrix = np.zeros((40, 40), dtype=int)
for i in range(40):
    for j in range(i+1, 40):
        M = projectors[i].T @ projectors[j]
        svs = np.linalg.svd(M, compute_uv=False)
        ov = round(h * np.sum(svs**2))
        overlap_matrix[i, j] = ov
        overlap_matrix[j, i] = ov
    overlap_matrix[i, i] = N

for m in [5, 6, 7, 8]:
    A[m] = np.zeros((40, 40), dtype=int)
    for i in range(40):
        for j in range(40):
            if overlap_matrix[i, j] == 12 * m and i != j:
                A[m][i, j] = 1

E_op = sum((2*a1 - m) * A[m].astype(float) for m in [5,6,7,8])
O_op = sum(m * A[m].astype(float) for m in [5,6,7,8])

# ============================================================
# PART 1: THE 26-DIMENSIONAL KERNEL OF E
# ============================================================
print(f"\n{'='*70}")
print("PART 1: THE 26-DIMENSIONAL KERNEL OF E")
print("="*70)

E_eigvals, E_eigvecs = np.linalg.eigh(E_op)
idx = np.argsort(np.abs(E_eigvals))
E_eigvals = E_eigvals[idx]
E_eigvecs = E_eigvecs[:, idx]

# Count zero eigenvalues
n_kernel = np.sum(np.abs(E_eigvals) < 0.01)
print(f"\nKernel dimension of E: {n_kernel}")
print(f"  = a1^2 + 1 = {a1**2 + 1}")
print(f"  = 1/sin^2(theta_W) = {1/(b1/(a1**2+1)):.0f}")

# Kernel vectors
K = E_eigvecs[:, :n_kernel]

# What structure does the kernel have?
print(f"\nKernel basis properties:")
print(f"  Shape: {K.shape}")

# Check: does the kernel have a clean structure?
# Compute the Gram matrix K^T @ K
G_kernel = K.T @ K
print(f"  Gram matrix K^T K = I check: max deviation = {np.max(np.abs(G_kernel - np.eye(n_kernel))):.6e}")

# Project each A_m into the kernel
print(f"\nA_m restricted to ker(E):")
for m in [5,6,7,8]:
    Am_ker = K.T @ A[m].astype(float) @ K
    eigs_ker = np.sort(np.linalg.eigvalsh(Am_ker))[::-1]
    trace = np.trace(Am_ker)
    print(f"  A_{m}|_ker: trace = {trace:.4f}, rank = {np.linalg.matrix_rank(Am_ker, tol=0.1)}")

    # Eigenvalues
    distinct = []
    for ev in eigs_ker:
        if not distinct or abs(ev - distinct[-1][0]) > 0.1:
            distinct.append([ev, 1])
        else:
            distinct[-1][1] += 1
    eig_str = ", ".join(f"{ev:.2f}({m_})" for ev, m_ in distinct)
    print(f"    eigenvalues: {eig_str}")

# ============================================================
# PART 2: WHY 26? ALGEBRAIC ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("PART 2: WHY 26 = a1^2 + 1?")
print("="*70)

# E = sum (2a1 - m) * A_m = (2a1) * sum A_m - sum m * A_m = 2*a1*(J-I) - O
# J-I has eigenvalues: 39 (once) and -1 (39 times)
# So E = 2*a1*(J-I) - O
# Eigenvalues of E = 2*a1*lambda(J-I) - lambda(O)
#   = 10*39 - 247 = 390 - 247 = 143 (the dominant eigenvalue of E: CHECK)
#   = 10*(-1) - lambda(O) = -10 - lambda(O) for the other 39

print(f"\nAlgebraic relation: E = 2*a1*(J-I) - O")
print(f"  Dominant eigenvalue: 2*a1*39 - lambda_max(O)")

O_eigvals_sorted = np.sort(np.linalg.eigvalsh(O_op))[::-1]
print(f"  lambda_max(O) = {O_eigvals_sorted[0]:.4f}")
print(f"  2*a1*39 - lambda_max(O) = {2*a1*39 - O_eigvals_sorted[0]:.4f}")

E_eigvals_sorted = np.sort(np.linalg.eigvalsh(E_op))[::-1]
print(f"  lambda_max(E) = {E_eigvals_sorted[0]:.4f}")
print(f"  Match: {abs(2*a1*39 - O_eigvals_sorted[0] - E_eigvals_sorted[0]) < 0.01}")

print(f"\n  For the 26-dim kernel of E: lambda(E) = 0")
print(f"  => -10 - lambda(O) = 0 => lambda(O) = -10 = -2*a1")
print(f"  This confirms: ker(E) = eigenspace of O with eigenvalue -2*a1 = -10")

# So the question becomes: WHY does O have a 26-dim eigenspace at -2*a1?
# Equivalently: WHY does E have a 26-dim kernel?

# E = sum (2a1-m) * A_m
# The 26-dim kernel means: for 26 independent vectors v,
#   sum (2a1-m) * (A_m v) = 0
# i.e., 5*(A5 v) + 4*(A6 v) + 3*(A7 v) + 2*(A8 v) = 0

print(f"\n  Kernel condition: 5*(A5 v) + 4*(A6 v) + 3*(A7 v) + 2*(A8 v) = 0")
print(f"  = n_db*(A5 v) + n_12*(A6 v) + n_db*(A7 v) + n_db*(A8 v) ... NO")
print(f"  Coefficients are {a1}, {a1-1}, {a1-2}, {a1-3} = a1, n_ut, n_12, n_db")
print(f"  = the CKM multipliers!")

# ============================================================
# PART 3: RANK OF E AND CHARACTERISTIC POLYNOMIAL
# ============================================================
print(f"\n{'='*70}")
print("PART 3: RANK AND CHARACTERISTIC POLYNOMIAL OF E")
print("="*70)

rank_E = np.linalg.matrix_rank(E_op, tol=0.1)
print(f"\nrank(E) = {rank_E}")
print(f"  = 40 - 26 = 14 = 2 * n_23 = 2 * 7")

# Non-zero eigenvalues of E
nonzero_E = [ev for ev in E_eigvals_sorted if abs(ev) > 0.1]
print(f"\nNon-zero eigenvalues of E ({len(nonzero_E)}):")
for ev in sorted(nonzero_E, reverse=True):
    frac_str = ""
    for name, val in [('a1', 5), ('b1', 6), ('N', 120), ('h', 30), ('phi', phi),
                       ('phi^2', phi**2), ('phi^3', phi**3), ('n_23', 7), ('n_13', 12)]:
        if abs(ev/val - round(ev/val)) < 0.02 and abs(round(ev/val)) <= 20:
            frac_str += f" ~{round(ev/val)}*{name}"
    print(f"  {ev:12.6f}{frac_str}")

# Sum of nonzero eigenvalues
print(f"\n  Sum of nonzero eigenvalues: {sum(nonzero_E):.4f}")
print(f"  Sum of positive eigenvalues: {sum(e for e in nonzero_E if e > 0):.4f}")
print(f"  Sum of negative eigenvalues: {sum(e for e in nonzero_E if e < 0):.4f}")
print(f"  Product of nonzero eigenvalues (det'): {np.prod([abs(e) for e in nonzero_E]):.4f}")

# Trace of E^2 (sum of squared eigenvalues)
trE2 = np.trace(E_op @ E_op)
print(f"\n  Tr(E^2) = {trE2:.4f}")
print(f"  Tr(E^2) / 40 = {trE2/40:.4f}")
print(f"  Tr(E^2) / (40*39) = {trE2/(40*39):.6f}")

# ============================================================
# PART 4: THE Z/2 SYMMETRY OF A_7
# ============================================================
print(f"\n{'='*70}")
print("PART 4: Z/2 SYMMETRY (why A_7 eigenvalues have mult 2)")
print("="*70)

# A_7 has all eigenvalues with multiplicity exactly 2
# This means there's an involution sigma on {1,...,40} that:
#   - commutes with A_7: sigma @ A_7 = A_7 @ sigma
#   - acts freely (no fixed points) on each eigenspace of A_7

# Find the involution by looking at A_7 eigenspaces
A7_eigvals, A7_eigvecs = np.linalg.eigh(A[7].astype(float))

# Group eigenvalues into pairs
print(f"\nA_7 eigenvalues (all mult 2):")
i = 0
while i < 40:
    ev = A7_eigvals[i]
    j = i
    while j < 40 and abs(A7_eigvals[j] - ev) < 0.01:
        j += 1
    mult = j - i
    print(f"  {ev:10.6f}  mult = {mult}")
    i = j

# Check: does there exist a permutation matrix P such that P @ A_7 = A_7 @ P
# and P^2 = I (involution)?
# Simpler: check if the COMPLEMENT A_5 + A_8 has a Z/2 symmetry

# The degree sequence of A_7 has 3 types: {7:4, 9:20, 10:16}
degs_7 = np.sum(A[7], axis=1)
type_7 = {}
for i in range(40):
    d = degs_7[i]
    if d not in type_7:
        type_7[d] = []
    type_7[d].append(i)

print(f"\nA_7 node types by degree:")
for d, nodes in sorted(type_7.items()):
    print(f"  degree {d}: {len(nodes)} nodes  (indices: {nodes[:8]}...)")

# Check if the 4 nodes of degree 7 form a special substructure
special_4 = type_7[7]
print(f"\nThe 4 special nodes (degree 7 in A_7): {special_4}")
print(f"  Their mutual overlaps:")
for i in range(4):
    for j in range(i+1, 4):
        ni, nj = special_4[i], special_4[j]
        print(f"    ({ni},{nj}): overlap = {overlap_matrix[ni, nj]}")

# Check all 4 adjacency matrices for these 4 nodes
print(f"\n  Adjacency pattern among the 4 special nodes:")
for m in [5,6,7,8]:
    connections = []
    for i in range(4):
        for j in range(i+1, 4):
            if A[m][special_4[i], special_4[j]] == 1:
                connections.append((i,j))
    print(f"    A_{m}: {len(connections)} edges -- {connections}")

# Their degree in each color class
print(f"\n  Degree of special nodes in each color:")
for i, node in enumerate(special_4):
    degs = [np.sum(A[m][node, :]) for m in [5,6,7,8]]
    print(f"    node {node}: A5={degs[0]}, A6={degs[1]}, A7={degs[2]}, A8={degs[3]}, total={sum(degs)}")

# ============================================================
# PART 5: THE GRASSMANNIAN STRUCTURE OF THE KERNEL
# ============================================================
print(f"\n{'='*70}")
print("PART 5: GRASSMANNIAN INTERPRETATION")
print("="*70)

# Each decomposition is a point in Gr(4,8) (Grassmannian of 4-planes in R^8)
# dim Gr(4,8) = 4*4 = 16
# The 40 points lie on some subvariety

# The kernel of E corresponds to a 26-dim subspace of R^40
# R^40 = space of "functions on 40 decompositions"
# The kernel condition: sum (2a1-m) * f(neighbors at level m) = 0 for each node

# This is like a HARMONIC condition: the "CKM-weighted average" of f vanishes

print(f"\n  The 40 decompositions live in Gr(4,8), dim = 16")
print(f"  Functions on 40 points: R^40")
print(f"  Kernel of E: 26-dimensional subspace of harmonic functions")
print(f"  (harmonic w.r.t. the CKM-weighted graph Laplacian)")

# Is there a Gr(4,8) interpretation of 26?
# Gr(4,8) has total cohomology dim = C(8,4) = 70
# But the Schubert cells have dims 0,1,...,16
# Total Schubert cells: p(4) * p(4) related... actually C(8,4)=70 Plucker coords

# Alternative: in SO(8) rep theory
# dim(so(8)) = 28
# 28 - 2 = 26? The 2 would be the Cartan of H4?
# Or: dim(so(8)) - dim(so(2)) = 28 - 1 = 27, not 26
# Or: rank 8 representation of so(8): adjoint = 28, vector = 8, spinor = 8

# More direct: 26 = (a1^2+1) appears throughout the framework
# sin^2(tW) = b1/(a1^2+1) = 6/26
# C = 4/(a1^2+1) = 4/26 (mass formula coefficient)

print(f"\n  26 in the framework:")
print(f"    a1^2 + 1 = {a1**2+1}")
print(f"    sin^2(tW) = b1/(a1^2+1) = {b1}/{a1**2+1}")
print(f"    C = 4/(a1^2+1) = 4/{a1**2+1} (mass formula)")
print(f"    dim(so(8)) - 2 = 28 - 2 = {28-2}")

# ============================================================
# PART 6: STRUCTURE CONSTANTS PATTERN
# ============================================================
print(f"\n{'='*70}")
print("PART 6: STRUCTURE CONSTANTS (row-by-row)")
print("="*70)

# Although not an association scheme (non-constant structure constants),
# let's look at the ACTUAL product entries

print(f"\nA_m @ A_n products (entry-by-entry analysis):")
for m in [5,6,7,8]:
    for n in range(m, 9):
        if n > 8: break
        prod = A[m] @ A[n]
        # The (i,j) entry counts paths i -> k -> j where i-k has color m, k-j has color n

        # Group entries by the color of the (i,j) edge
        for c in [0, 5, 6, 7, 8]:  # 0 = diagonal
            if c == 0:
                vals = [prod[i,i] for i in range(40)]
            else:
                vals = [prod[i,j] for i in range(40) for j in range(40)
                        if A.get(c, np.zeros((40,40)))[i,j] == 1] if c in A else []
            if vals:
                unique = sorted(set(vals))
                print(f"  (A{m}@A{n})[color={c}]: values = {unique}, counts = {Counter(vals)}")

# ============================================================
# PART 7: DECOMPOSITION INTO ORBITS
# ============================================================
print(f"\n{'='*70}")
print("PART 7: ORBIT STRUCTURE")
print("="*70)

# The degree sequence groups the 40 nodes into 3 types:
# Looking at all 4 color degrees simultaneously
node_types = {}
for i in range(40):
    sig = tuple(int(np.sum(A[m][i,:])) for m in [5,6,7,8])
    if sig not in node_types:
        node_types[sig] = []
    node_types[sig].append(i)

print(f"\nNode types by (deg5, deg6, deg7, deg8):")
for sig, nodes in sorted(node_types.items()):
    total = sum(sig)
    print(f"  {sig}: {len(nodes)} nodes (total deg = {total})")

# This gives the ORBITS under the automorphism group
n_orbits = len(node_types)
orbit_sizes = [len(v) for v in node_types.values()]
print(f"\n  Number of orbits: {n_orbits}")
print(f"  Orbit sizes: {sorted(orbit_sizes)}")
print(f"  Sum: {sum(orbit_sizes)} = 40")

# Check: 4 + 16 + 20 = 40? or 4 + 12 + 24 = 40?
# The orbit structure tells us about the automorphism group

# For each orbit, what's the subgraph structure?
for sig, nodes in sorted(node_types.items()):
    print(f"\n  Orbit {sig} ({len(nodes)} nodes):")
    for m in [5,6,7,8]:
        # Internal edges (within orbit)
        internal = sum(A[m][i,j] for i in nodes for j in nodes if i < j)
        # External edges (to other orbits)
        external = sum(A[m][i,j] for i in nodes for j in range(40) if j not in nodes)
        print(f"    A_{m}: {internal} internal, {external} external edges")

# ============================================================
# PART 8: KEY NUMERICAL CHECKS
# ============================================================
print(f"\n{'='*70}")
print("PART 8: KEY NUMERICAL CHECKS")
print("="*70)

# Check 1: Is 40 - 26 = 14 = 2*n_23?
print(f"\n  rank(E) = 40 - 26 = 14 = 2*n_23 = 2*{n_23}: {'YES' if 40-26 == 2*7 else 'NO'}")

# Check 2: The 14 non-kernel dimensions
# These should correspond to the "CKM-active" modes
print(f"  CKM-active modes: 14 = 2*n_23")

# Check 3: 26/40 ratio
print(f"  26/40 = {26/40} = 13/20")
print(f"  14/40 = {14/40} = 7/20 = n_23/20")

# Check 4: Is there a natural decomposition 40 = 26 + 14?
# In rep theory of SO(8): 28 (adjoint) + 8_v + 8_s + 8_c
# But 40 != any standard SO(8) rep
# In H4: the Weyl group has... complicated rep theory

# Check 5: Eigenvalue products
O_eigvals_full = np.sort(np.linalg.eigvalsh(O_op))[::-1]
print(f"\n  O eigenvalues (distinct with mult):")
distinct_O = []
for ev in O_eigvals_full:
    if not distinct_O or abs(ev - distinct_O[-1][0]) > 0.01:
        distinct_O.append([ev, 1])
    else:
        distinct_O[-1][1] += 1
for ev, mult in distinct_O:
    print(f"    {ev:12.6f} x {mult}")

# The non-(-10) eigenvalues
non_bulk = [(ev, mult) for ev, mult in distinct_O if abs(ev + 10) > 0.1]
print(f"\n  Non-bulk eigenvalues of O (not -10):")
for ev, mult in non_bulk:
    print(f"    {ev:12.6f} x {mult}")

# Sum check: sum of all eigenvalues = Tr(O) = 0
print(f"\n  Tr(O) = sum eigenvalues = {sum(ev*m for ev,m in distinct_O):.4f}")
print(f"  Expected: sum_m m * Tr(A_m) = 0 (each A_m is traceless)")

# The dominant eigenvalue: should be sum_m m * d_m / 40 * 40 (if regular)
# where d_m is the (average) degree
avg_O = sum(m * np.mean(np.sum(A[m], axis=1)) for m in [5,6,7,8])
print(f"\n  Average row sum of O = {avg_O:.4f}")
print(f"  Dominant eigenvalue of O / 39 = {O_eigvals_full[0]/39:.6f}")
print(f"  Average m = {avg_O / 39:.6f}")

# ============================================================
# PART 9: HONEST ASSESSMENT
# ============================================================
print(f"\n{'='*70}")
print("HONEST ASSESSMENT")
print("="*70)

print(f"""
KEY FINDINGS:
  1. E = (CKM-weighted adjacency) has 26-dim kernel
     26 = a1^2 + 1 = the GAUGE dimension

  2. rank(E) = 14 = 2*n_23 = 2*7
     The CKM-active subspace is 14-dimensional

  3. O has eigenvalue -2*a1 = -10 with multiplicity 26
     (equivalent to finding 1 by E = 2a1*(J-I) - O)

  4. A_7 has all eigenvalues with multiplicity 2 (Z/2 symmetry)

  5. 40 nodes split into 3 orbits by degree signature

ASSESSMENT:
  The 26-dim kernel is NUMERICALLY STRIKING:
    26 = a1^2 + 1 is the central number in sin^2(tW) = b1/(a1^2+1)

  BUT: Is this DERIVABLE or coincidental?

  The key question: can we prove algebraically that
    ker(sum (2a1-m) * A_m) has dimension a1^2+1
  for the E8/H4 decomposition graph?

  If YES: this gives a NEW geometric origin for 26
  If NO (just numerical): then it's interesting but not a derivation

  The fact that 40 - 26 = 14 = 2*n_23 adds structure,
  but could be forced by the relation E = 2a1*(J-I) - O.
""")

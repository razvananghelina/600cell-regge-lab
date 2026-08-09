"""
exp343b: Derive Mass Formula Coefficients (5,6) from McKay Graph
================================================================
The mass formula n = 5a + 6b = a1*a + b1*b uses coefficients (a1, b1).
Can we derive WHY these are the coefficients from the McKay graph?

KEY INSIGHT FROM exp343a:
  z_mass = a1 + b1*phi has N(z_mass) = 19 = unique prime.
  The coefficients (5,6) ARE (a1, b1) by definition. But WHY a1 and b1?

Approaches:
  1. McKay Cayley matrix eigenvectors
  2. Perron-Frobenius analysis
  3. Graph Laplacian / distance functions
  4. The mass as an "evaluation" functional on the character ring
  5. Height function on the affine E8 Dynkin diagram

CATEGORY: EXPLORATORY
"""

import numpy as np
from math import factorial, gcd
import time

t0 = time.time()

# Framework
a1 = 5
b1 = a1 + 1  # 6
PHI = (1 + np.sqrt(5)) / 2
phi_prime = (1 - np.sqrt(5)) / 2
N = factorial(a1)  # 120
N_gen = 3
N_eig = 9
rank_E8 = 8

def norm_zphi(a, b):
    return a*a + a*b - b*b

def trace_zphi(a, b):
    return 2*a + b

print("=" * 72)
print("EXP343b: DERIVE MASS FORMULA FROM McKAY GRAPH")
print("=" * 72)

# ============================================================
# SECTION 1: THE McKAY GRAPH AND CAYLEY MATRIX
# ============================================================
print("\n" + "=" * 72)
print("SECTION 1: McKAY GRAPH STRUCTURE")
print("=" * 72)

# The McKay graph of 2I is the affine E8 Dynkin diagram.
# 9 nodes with dimensions (Dynkin labels): [1, 2, 3, 4, 5, 6, 4, 2, 3]
# The graph structure:
#
# Numbering (following standard affine E8):
# Node 0: dim 1 (extending node)
# Node 1: dim 2
# Node 2: dim 3
# Node 3: dim 4
# Node 4: dim 5
# Node 5: dim 6 (branch point)
# Node 6: dim 4
# Node 7: dim 2
# Node 8: dim 3
#
# Edges: 0-1, 1-2, 2-3, 3-4, 4-5, 5-6, 6-7, 5-8
# (8 edges total)
#
# Topology: main chain 0-1-2-3-4-5, with branches 5-6-7 and 5-8

dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])
node_labels = ['0(1)', '1(2)', '2(3)', '3(4)', '4(5)', '5(6)', '6(4)', '7(2)', '8(3)']

# Adjacency matrix of affine E8
A_mckay = np.zeros((9, 9), dtype=int)
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
for i, j in edges:
    A_mckay[i,j] = 1
    A_mckay[j,i] = 1

print(f"\n  Affine E8 adjacency matrix:")
print(f"  Node dims: {list(dims)}")
print(f"  Edges: {edges}")
print(f"  Branch point: node 5 (dim 6 = b1)")
print(f"  Legs from branch: (5,4,3,2,1), (4,2), (3)")
print(f"  Leg lengths: 5, 2, 1")

# McKay Cayley matrix: C_{ij} = number of times rho_j appears in rho_i x rho_fund
# For affine E8, the fundamental representation has dim 2.
# C_{ij} = A_{ij} (the adjacency matrix) since the McKay graph is the
# graph where edges represent tensor product decomposition.
# But the CAYLEY matrix used in spectral analysis includes dimension weighting.

# The McKay Cayley matrix (normalized): M_{ij} = A_{ij} * d_j / d_i
# This is NOT symmetric. The symmetric version is:
# M_sym = D^{-1/2} * A * D^{1/2} where D = diag(dims)

D = np.diag(dims.astype(float))
D_inv = np.diag(1.0 / dims.astype(float))
D_sqrt = np.diag(np.sqrt(dims.astype(float)))
D_inv_sqrt = np.diag(1.0 / np.sqrt(dims.astype(float)))

# Cayley matrix: C = D^{-1} * A * D (acts on characters)
C_cayley = D_inv @ A_mckay @ D

# Symmetric Cayley: S = D^{-1/2} * A * D^{1/2}
S_cayley = D_inv_sqrt @ A_mckay @ D_sqrt

# Also: the "normalized Laplacian" of McKay graph
# L = I - D^{-1} * A (random walk Laplacian)
L_rw = np.eye(9) - D_inv @ A_mckay.astype(float)

# Vertex Laplacian (degree - adjacency)
deg_mckay = np.diag(A_mckay.sum(axis=1).astype(float))
L_mckay = deg_mckay - A_mckay.astype(float)

print(f"\n  Cayley matrix C = D^{{-1}} A D:")
for i in range(9):
    row_str = " ".join(f"{C_cayley[i,j]:>6.3f}" for j in range(9))
    print(f"    [{row_str}]")

# ============================================================
# SECTION 2: EIGENANALYSIS
# ============================================================
print("\n" + "=" * 72)
print("SECTION 2: EIGENANALYSIS OF McKAY MATRICES")
print("=" * 72)

# Eigenvalues of A (adjacency)
eig_A, vec_A = np.linalg.eigh(A_mckay.astype(float))
idx_A = np.argsort(-eig_A)
eig_A = eig_A[idx_A]
vec_A = vec_A[:, idx_A]

print(f"\n  Adjacency matrix eigenvalues:")
for i, ev in enumerate(eig_A):
    print(f"    lambda_{i} = {ev:>10.6f}")

# Eigenvalues of Cayley matrix C
eig_C, vec_C_r = np.linalg.eig(C_cayley)
idx_C = np.argsort(-np.real(eig_C))
eig_C = eig_C[idx_C]
vec_C_r = vec_C_r[:, idx_C]

print(f"\n  Cayley matrix C eigenvalues:")
for i, ev in enumerate(eig_C):
    ev_r = np.real(ev)
    print(f"    mu_{i} = {ev_r:>10.6f}")

# Check: are these the McKay spectrum {0, +-1/phi, +-1, +-phi, +-2}?
mckay_spectrum = sorted([0, 1/PHI, -1/PHI, 1, -1, PHI, -PHI, 2, -2], reverse=True)
print(f"\n  Expected McKay spectrum: {[f'{x:.6f}' for x in mckay_spectrum]}")
print(f"  Computed:                {[f'{np.real(e):.6f}' for e in eig_C]}")
match = all(abs(np.real(eig_C[i]) - mckay_spectrum[i]) < 1e-6 for i in range(9))
print(f"  Match: {'YES' if match else 'NO'}")

# Perron vector (largest eigenvalue = 2)
print(f"\n  Perron eigenvector (lambda = 2):")
perron = vec_A[:, 0] if eig_A[0] > 1.5 else vec_A[:, -1]
# Normalize to match dims
perron = perron / perron[0] * dims[0]
# Check if proportional to dims
ratio = perron / dims
print(f"    Perron = [{', '.join(f'{p:.4f}' for p in perron)}]")
print(f"    Dims   = {list(dims)}")
is_perron_dims = np.allclose(ratio, ratio[0], atol=1e-4)
print(f"    Perron proportional to dims: {is_perron_dims}")

# Symmetric Cayley eigenanalysis
eig_S, vec_S = np.linalg.eigh(S_cayley)
idx_S = np.argsort(-eig_S)
eig_S = eig_S[idx_S]
vec_S = vec_S[:, idx_S]

print(f"\n  Symmetric Cayley S eigenvalues:")
for i, ev in enumerate(eig_S):
    print(f"    sigma_{i} = {ev:>10.6f}")

# ============================================================
# SECTION 3: CAN WE RECOVER (a1, b1) FROM EIGENVECTORS?
# ============================================================
print("\n" + "=" * 72)
print("SECTION 3: RECOVERING MASS COEFFICIENTS FROM EIGENVECTORS")
print("=" * 72)

# The mass formula n_k = a1*a_k + b1*b_k = 5*a_k + 6*b_k
# for the fermion at McKay node k.

# Fermion assignment (from exp323):
# Node 0 (dim 1): electron, n=0
# Node 1 (dim 2): ?, n=?
# Node 2 (dim 3): ?, n=?
# etc.
# The exact assignment depends on the fermion-to-node mapping.

# Known fermion data:
fermions = {
    'e':   {'ab': (0, 0),  'n': 0},
    'u':   {'ab': (3, -2), 'n': 3},
    'd':   {'ab': (1, 0),  'n': 5},
    'mu':  {'ab': (1, 1),  'n': 11},
    's':   {'ab': (1, 1),  'n': 11},
    'c':   {'ab': (2, 1),  'n': 16},
    'tau': {'ab': (1, 2),  'n': 17},
    'b':   {'ab': (-1, 4), 'n': 19},
    't':   {'ab': (4, 1),  'n': 26},
}

# McKay node assignment (from exp323, Solution 3):
# The main chain [0-1-2-3-4-5] maps to:
# w=0: e(n=0), w=1: ?(n=?), w=2: ?(n=?), w=3: ?(n=?), ...
# From memory: main chain FORCED: [N_gen, F(3), b1, 0, a1] = [3,2,6,0,5]
# This seems to be weights, not exponents.

# Let me use the weight assignment from exp323:
# McKay weights (graph distance or similar):
# Solution 3: weights = {0, 1, 2, 3, 3, 5, 6, 9}
# with mu-s degeneracy (both have same weight)

# Node-to-weight (Solution 3 from exp323):
# The fermion-to-node assignment:
# e -> node 0 (dim 1, weight 0)
# u -> node ? (weight 3)
# d -> node 1 (dim 2, weight 1)... this isn't quite right

# Let me approach differently. The mass exponents are:
mass_exponents = {'e': 0, 'u': 3, 'd': 5, 's': 11, 'mu': 11, 'c': 16, 'tau': 17, 'b': 19, 't': 26}

# The 9 exponents (with mu=s degeneracy): 0, 3, 5, 11, 11, 16, 17, 19, 26
exps_sorted = sorted(mass_exponents.values())
print(f"\n  Mass exponents: {exps_sorted}")
print(f"  Sum = {sum(exps_sorted)} = N_eig * vertex_degree = {N_eig * 12}")

# Can these 9 numbers be expressed as eigenvalues or components of
# some natural 9-dimensional object on the McKay graph?

# Test 1: Are exponents linear combinations of McKay eigenvalues?
print(f"\n  Test 1: Linear decomposition of exponents in McKay eigenvalues")
# Eigenvalues: 2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2
mckay_eigs = np.array([2, PHI, 1, 1/PHI, 0, -1/PHI, -1, -PHI, -2])
# Try: n_k = c1*mu_1 + c2*mu_2 + ... (linear in eigenvalues)
# This doesn't make sense directly since n_k is per-node, not per-eigenvalue.

# Test 2: Are exponents related to Cayley Laplacian eigenvalues?
# From exp342e: lambda_k = 12(1 - chi_k(10a)/d_k)
# The Cayley Laplacian on 2I regular representation.
# These eigenvalues are per-irrep, with multiplicity d_k^2.
# We already know from exp342c that "Cayley eigenvalues NOT proportional to mass"

# Test 3: The fundamental question - WHY 5 and 6?
print(f"\n  Test 3: WHY coefficients (5, 6) = (a1, a1+1)?")
print(f"  The mass formula n = a1*a + b1*b in the (a,b) coordinates of Z[phi].")
print(f"  z = a + b*phi is the algebraic integer, n = a1*a + b1*b is the 'height'.")
print(f"")

# In Z[phi], every element z = a + b*phi has:
# Trace: Tr(z) = z + z' = 2a + b
# Norm: N(z) = z*z' = a^2 + ab - b^2
# "a1-height": h(z) = a1*a + b1*b = 5a + 6b

# Can we express the height h in terms of Tr and N?
# h = 5a + 6b
# Tr = 2a + b => a = (Tr - b)/2, so h = 5*(Tr-b)/2 + 6b = 5*Tr/2 + 7*b/2
# This mixes Tr and b (not clean).

# Alternative: use the "evaluation map"
# If we evaluate z at phi -> b1/a1 = 6/5:
# z(6/5) = a + b*(6/5) = (5a + 6b)/5 = n/5
# So n = 5 * z(b1/a1)!

print(f"  KEY OBSERVATION:")
print(f"  n = a1 * z(b1/a1) where z(x) = a + b*x")
print(f"  i.e., n = 5 * z(6/5) = 5*(a + 6b/5) = 5a + 6b")
print(f"  The mass exponent is a1 times z evaluated at b1/a1 = 6/5!")

# But WHY evaluate at 6/5 and not at phi?
# If we evaluate at phi: z(phi) = a + b*phi (the actual algebraic integer)
# If we evaluate at b1/a1 = 6/5:
# z(6/5) = a + 6b/5 = n/5

# What is 6/5 in the framework?
print(f"\n  What is b1/a1 = 6/5 = 1.2?")
print(f"    phi = {PHI:.6f}")
print(f"    b1/a1 = {b1/a1:.6f}")
print(f"    phi - b1/a1 = {PHI - b1/a1:.6f}")
print(f"    b1/a1 = phi - (phi-1-1/a1) = phi - {PHI-1-1/a1:.6f}")
print(f"    b1/a1 = (a1+1)/a1 = 1 + 1/a1 = {1+1/a1:.6f}")
print(f"    phi^2 - b1/a1 = {PHI**2 - b1/a1:.6f}")

# Actually: n = a1*a + (a1+1)*b = a1*(a+b) + b = a1*Tr(z)/2 + b*(1+a1/2)
# Hmm, still not clean.

# Better: the evaluation at b1/a1 = 1+1/a1 is the RATIONAL APPROXIMATION of phi!
# phi ~ 1.618, 1+1/a1 = 1.200. Not a great approximation...
# But 1+1/a1 appears in the continued fraction: phi = [1; 1, 1, 1, ...]
# Convergents: 1/1, 2/1, 3/2, 5/3, 8/5, 13/8, 21/13, ...
# 6/5 is NOT a convergent of phi!

# Test 4: McKay graph and the Cartan matrix
print(f"\n  Test 4: Cartan Matrix Approach")
# The Cartan matrix of affine E8: C = 2I - A
C_cartan = 2 * np.eye(9) - A_mckay.astype(float)
print(f"  Cartan matrix C = 2I - A:")
for i in range(9):
    row_str = " ".join(f"{C_cartan[i,j]:>3.0f}" for j in range(9))
    print(f"    [{row_str}]")

# The null vector of C_cartan is the Dynkin marks = dims
null_C = dims.astype(float) / np.linalg.norm(dims)
check = C_cartan @ dims.astype(float)
print(f"  C * dims = [{', '.join(f'{x:.1f}' for x in check)}]")
print(f"  (should be zero: {'YES' if np.allclose(check, 0) else 'NO'})")

# The Cartan matrix inverse (pseudoinverse, since det=0):
# Could the mass exponents be related to the pseudoinverse?
C_pinv = np.linalg.pinv(C_cartan)
print(f"\n  Pseudoinverse of Cartan:")
for i in range(9):
    row_str = " ".join(f"{C_pinv[i,j]:>7.3f}" for j in range(9))
    print(f"    [{row_str}]")

# ============================================================
# SECTION 4: HEIGHT FUNCTION ON DYNKIN DIAGRAM
# ============================================================
print("\n" + "=" * 72)
print("SECTION 4: HEIGHT FUNCTION ON AFFINE E8")
print("=" * 72)

# On the affine E8 diagram, there's a natural "height" or "level" function.
# The Dynkin marks d_i satisfy: sum_j A_{ij} * d_j = 2 * d_i for affine case.
# This means the marks are the Perron eigenvector.

# The "imaginary root" delta = sum n_i * alpha_i where n_i = marks = dims
# For affine E8: delta = 1*a0 + 2*a1 + 3*a2 + 4*a3 + 5*a4 + 6*a5 + 4*a6 + 2*a7 + 3*a8
# This is the null root of the Cartan matrix.

# A natural height function h(node_i) could be:
# h_i = (graph distance from node 0) * (Dynkin mark d_i)
# or h_i = some weighted sum

graph_dist = np.zeros(9, dtype=int)
# BFS from node 0
from collections import deque
queue = deque([0])
visited = {0}
while queue:
    node = queue.popleft()
    for j in range(9):
        if A_mckay[node, j] and j not in visited:
            graph_dist[j] = graph_dist[node] + 1
            visited.add(j)
            queue.append(j)

print(f"\n  Graph distances from node 0:")
for i in range(9):
    print(f"    Node {i} (dim {dims[i]}): dist = {graph_dist[i]}")

# Height candidates
print(f"\n  Height function candidates:")
h1 = graph_dist * dims  # dist * dim
h2 = graph_dist + dims  # dist + dim
h3 = graph_dist * a1 + dims  # dist*a1 + dim
h4 = dims * a1 + graph_dist * b1  # dim*a1 + dist*b1
h5 = 2 * graph_dist * dims + dims - 1  # ad hoc

print(f"  {'Node':>6} {'dim':>4} {'dist':>5} {'d*dist':>7} {'d+dist':>7} {'d*5+dist':>8} {'5d+6dist':>9}")
for i in range(9):
    print(f"  {i:>6} {dims[i]:>4} {graph_dist[i]:>5} {h1[i]:>7} {h2[i]:>7} {dims[i]*a1+graph_dist[i]:>8} {dims[i]*a1+graph_dist[i]*b1:>9}")

# The sorted mass exponents are: 0, 3, 5, 11, 11, 16, 17, 19, 26
# Can any of these height functions reproduce them?
target_exps_sorted = [0, 3, 5, 11, 11, 16, 17, 19, 26]

# Try systematic: n_node = c1*dim + c2*dist for various (c1, c2)
print(f"\n  Systematic search: n = c1*dim + c2*dist")
print(f"  Target (sorted): {target_exps_sorted}")

best_match = None
best_c = None
best_rms = float('inf')

for c1 in range(-10, 15):
    for c2 in range(-10, 15):
        if c1 == 0 and c2 == 0:
            continue
        h = c1 * dims + c2 * graph_dist
        h_sorted = sorted(h)
        # Check if matches target
        if list(h_sorted) == target_exps_sorted:
            print(f"  *** EXACT MATCH: c1={c1}, c2={c2} ***")
            print(f"      h = {list(h)}")
            best_match = h
            best_c = (c1, c2)
            best_rms = 0
        else:
            diff = np.array(h_sorted) - np.array(target_exps_sorted)
            rms = np.sqrt(np.mean(diff**2))
            if rms < best_rms:
                best_rms = rms
                best_c = (c1, c2)
                best_match = list(h)

print(f"\n  Best match: c1={best_c[0]}, c2={best_c[1]}, RMS={best_rms:.4f}")
print(f"  h = {best_match}")
print(f"  sorted: {sorted(best_match)}")
print(f"  target: {target_exps_sorted}")

# ============================================================
# SECTION 5: EFFECTIVE RESISTANCE / GREEN'S FUNCTION
# ============================================================
print("\n" + "=" * 72)
print("SECTION 5: GRAPH-THEORETIC DISTANCES")
print("=" * 72)

# Effective resistance between node 0 and node i
# R_eff(0, i) = (e_0 - e_i)^T * L^+ * (e_0 - e_i)
# where L^+ is the pseudoinverse of the Laplacian

L_pinv = np.linalg.pinv(L_mckay)

eff_resistance = np.zeros(9)
for i in range(9):
    e = np.zeros(9)
    e[0] = 1
    e[i] = -1
    eff_resistance[i] = e @ L_pinv @ e

print(f"\n  Effective resistance from node 0:")
for i in range(9):
    print(f"    R(0,{i}) = {eff_resistance[i]:.6f}")

# Green's function
print("\n  Green's function G(0, i) = L_pinv[0,i]:")
for i in range(9):
    print(f"    G(0,{i}) = {L_pinv[0,i]:.6f}")

# Heat kernel at various times
print("\n  Heat kernel K(0, i, t) = exp(-t*L)[0,i]:")
for t in [0.1, 0.5, 1.0, 2.0]:
    K = np.real(np.linalg.matrix_power(np.eye(9) - 0.01*L_mckay, int(t/0.01)))
    # Better: matrix exponential
    eigL, vecL = np.linalg.eigh(L_mckay)
    K = vecL @ np.diag(np.exp(-t * eigL)) @ vecL.T
    heat_vals = K[0, :]
    print(f"    t={t:.1f}: [{', '.join(f'{x:.4f}' for x in heat_vals)}]")

# ============================================================
# SECTION 6: THE FUNDAMENTAL THEOREM ATTEMPT
# ============================================================
print("\n" + "=" * 72)
print("SECTION 6: WHY n = a1*a + b1*b (ALGEBRAIC ARGUMENT)")
print("=" * 72)

print("""
  The mass formula n = a1*a + b1*b arises from:

  1. Each fermion is characterized by z = a + b*phi in Z[phi].
     (This is derived: the mass is m = m_e * phi^n for integer n,
      and n decomposes as n = 5a + 6b for unique (a,b) given constraints.)

  2. The mass exponent n = a1*a + b1*b IS the evaluation:
     n = [z evaluated at phi_rational] * a1
     where phi_rational = b1/a1 = (a1+1)/a1.

  3. But there's a DEEPER way to see this.
     The ring Z[phi] has two embeddings into R:
       sigma_1: phi -> phi = (1+sqrt(5))/2   (identity)
       sigma_2: phi -> phi' = (1-sqrt(5))/2   (Galois conjugate)

     Under sigma_1: z -> z = a + b*phi (the actual algebraic integer)
     Under sigma_2: z -> z' = a + b*phi' = a - b/phi (Galois conjugate)

     The mass exponent n is NEITHER sigma_1(z) nor sigma_2(z).

  4. ALTERNATIVE: The mass formula comes from the INNER PRODUCT
     n = (w, z) in Z^2 where w = (a1, b1) is the "mass weight vector"
     and z = (a, b) is the fermion's Z[phi] coordinate.

     The weight vector w = (a1, b1) has N(w) = N(a1+b1*phi) = 19 = p.
     This is the UNIQUE prime from the uniqueness theorem.

  5. KEY QUESTION: Why is the mass weight vector w = (a1, a1+1)?

     ANSWER ATTEMPT: In the McKay graph (affine E8 Dynkin diagram),
     the Dynkin marks are d = (1, 2, 3, 4, 5, 6, 4, 2, 3).
     The sum of marks = 30 = h(E8) (Coxeter number).
     The highest mark = 6 = b1.
     The second-highest mark = 5 = a1.

     The mass weight vector (a1, b1) = (second-highest mark, highest mark)!
""")

# Verify: Dynkin marks sorted
marks = list(dims)
marks_sorted = sorted(marks, reverse=True)
print(f"  Dynkin marks sorted: {marks_sorted}")
print(f"  Top 2: ({marks_sorted[1]}, {marks_sorted[0]}) = ({a1}, {b1})? ", end="")
print("YES!" if marks_sorted[1] == a1 and marks_sorted[0] == b1 else "NO")

# ============================================================
# SECTION 7: THE COXETER ELEMENT APPROACH
# ============================================================
print("\n" + "=" * 72)
print("SECTION 7: COXETER ELEMENT AND EXPONENTS")
print("=" * 72)

# The Coxeter number of E8 is h = 30 = a1*b1.
# The exponents of E8 are: 1, 7, 11, 13, 17, 19, 23, 29
# These are the integers < 30 coprime to 30.

E8_exponents = [1, 7, 11, 13, 17, 19, 23, 29]
print(f"  E8 exponents: {E8_exponents}")
print(f"  Sum = {sum(E8_exponents)} = 120 = N = |2I|")
print(f"  Product = {np.prod(E8_exponents)}")
print(f"  Note: 19 appears as an E8 exponent!")
print(f"  And N(a1+b1*phi) = 19 = an E8 exponent.")
print(f"  Exponents coprime to h=30: {[k for k in range(1,30) if gcd(k,30)==1]}")

# The E8 exponents mod a1=5:
print(f"\n  E8 exponents mod a1 = {a1}:")
for e in E8_exponents:
    print(f"    {e} mod {a1} = {e % a1}")

# And mod b1=6:
print(f"\n  E8 exponents mod b1 = {b1}:")
for e in E8_exponents:
    print(f"    {e} mod {b1} = {e % b1}")

# Compare mass exponents with E8 exponents
print(f"\n  Mass exponents:  {sorted(mass_exponents.values())}")
print(f"  E8 exponents:    {E8_exponents}")
print(f"  Intersection: {sorted(set(mass_exponents.values()) & set(E8_exponents))}")
# 11, 17, 19 are in both lists!
intersect = sorted(set(mass_exponents.values()) & set(E8_exponents))
print(f"  Shared values: {intersect} ({len(intersect)} out of 8)")

# ============================================================
# SECTION 8: Z[phi] LATTICE AND QUADRATIC FORM
# ============================================================
print("\n" + "=" * 72)
print("SECTION 8: Z[phi] LATTICE INTERPRETATION")
print("=" * 72)

# The 9 fermion Z[phi] elements lie on a 2D lattice (a,b) in Z^2.
# The mass formula n = 5a + 6b defines a "height direction" (5,6).
# The level curves n = const are lines perpendicular to (5,6).

# The Z[phi] norm form: N(a+b*phi) = a^2 + ab - b^2
# This is an indefinite binary quadratic form of discriminant 5.
# It represents: {..., -19, -11, -9, -4, -1, 0, 1, 4, 5, 9, 11, 19, ...}

# The mass height function h(a,b) = 5a + 6b defines a foliation of Z^2.
# Each fermion sits on a specific h-level.

print(f"  Fermion positions in Z^2:")
print(f"  {'Name':>4} {'(a,b)':>8} {'n':>4} {'N(z)':>6} {'Tr(z)':>6}")
for name in ['e','u','d','s','mu','c','tau','b','t']:
    a, b = fermions[name]['ab']
    n = fermions[name]['n']
    nz = norm_zphi(a, b)
    tz = trace_zphi(a, b)
    print(f"  {name:>4} ({a:+2},{b:+2}) {n:>4} {nz:>6} {tz:>6}")

# The perpendicular direction to (5,6):
# (5,6) . (a,b) = 0 => 5a + 6b = 0 => a/b = -6/5
# Perpendicular direction: (-6, 5) or (6, -5)
print(f"\n  Mass direction: (5, 6)")
print(f"  Perpendicular: (-6, 5)")
print(f"  N(mass dir) = N(5+6*phi) = {norm_zphi(5,6)} = 19 (prime!)")
print(f"  N(perp dir) = N(-6+5*phi) = {norm_zphi(-6,5)}")
print(f"  Note: N(-6+5*phi) = 36 - 30 - 25 = {36-30-25}")

# The "perpendicular norm" distinguishes fermions at the same height n
print(f"\n  Perpendicular component p = -6*a + 5*b (mod gcd):")
for name in ['e','u','d','s','mu','c','tau','b','t']:
    a, b = fermions[name]['ab']
    n = fermions[name]['n']
    p = -6*a + 5*b
    print(f"  {name:>4}: n={n:>3}, p={p:>4}")

# The mu-s degeneracy: both have n=11 AND same (a,b)=(1,1)
# So they're at the SAME point in Z^2!

# ============================================================
# SECTION 9: THE a1, b1 AS DYNKIN MARKS
# ============================================================
print("\n" + "=" * 72)
print("SECTION 9: WHY a1 AND b1? - DYNKIN MARK ARGUMENT")
print("=" * 72)

print(f"""
  THEOREM ATTEMPT: n = a1*a + b1*b follows from the affine E8 structure.

  The affine E8 Dynkin diagram has marks d = (1,2,3,4,5,6,4,2,3).
  The two largest marks are d_max = 6 = b1 and d_sub = 5 = a1.
  These occur at the branch point (dim 6) and one step along the long leg (dim 5).

  In the representation ring R(2I) = Z[chi_1, ..., chi_9]:
  - Each irrep chi_k has dimension d_k (= Dynkin mark)
  - The fundamental representation chi_2 (dim 2) generates via tensor products
  - The mass exponent n_k should be the "level" of chi_k in this tower

  The tensor product tower:
    chi_2^0 = chi_1 (dim 1)
    chi_2^1 = chi_2 (dim 2)
    chi_2^2 = chi_1 + chi_3 (dims 1+3)
    chi_2^3 = 2*chi_2 + chi_4 (dims 2+2+4)
    chi_2^4 = chi_1 + 2*chi_3 + chi_5 (dims 1+3+3+5)
    chi_2^5 = 2*chi_2 + 2*chi_4 + chi_6 (dims 2+2+4+4+6)
    etc.

  The FIRST time chi_k appears in chi_2^n determines n_k(tower).
  For affine E8, this gives the "level" of each node.
""")

# Compute tensor product tower
# In the McKay graph, chi_2^n decomposition follows from the adjacency matrix.
# If we denote v_n = multiplicity vector at level n:
# v_0 = [1,0,0,0,0,0,0,0,0] (just chi_1)
# v_1 = A * v_0 (tensor with chi_2)
# v_n = A * v_{n-1}

v = np.zeros(9, dtype=float)
v[0] = 1  # start with trivial

first_appearance = -np.ones(9, dtype=int)
first_appearance[0] = 0

print(f"  Tensor product tower chi_2^n:")
print(f"  {'n':>3} {'decomposition':>45} {'new nodes':>15}")

for n in range(1, 20):
    v_new = A_mckay.astype(float) @ v
    new_nodes = []
    for i in range(9):
        if v_new[i] > 0.5 and first_appearance[i] < 0:
            first_appearance[i] = n
            new_nodes.append(i)
    v_str = ", ".join(f"{int(v_new[i])}*R{i}" for i in range(9) if v_new[i] > 0.5)
    new_str = ", ".join(f"R{i}(d={dims[i]})" for i in new_nodes)
    print(f"  {n:>3} {v_str:>45} {new_str:>15}")
    v = v_new
    if all(first_appearance >= 0):
        break

print(f"\n  First appearance level for each node:")
for i in range(9):
    print(f"    Node {i} (dim {dims[i]}): first at n={first_appearance[i]}")

print(f"\n  First appearance levels sorted: {sorted(first_appearance)}")
print(f"  Mass exponents sorted:          {sorted(mass_exponents.values())}")

# Are they related?
fa_sorted = sorted(first_appearance)
me_sorted = sorted(mass_exponents.values())
ratios = [me_sorted[i]/fa_sorted[i] if fa_sorted[i] > 0 else 0 for i in range(9)]
print(f"  Ratios (mass_exp / first_app): {[f'{r:.2f}' for r in ratios]}")

# ============================================================
# SECTION 10: WEIGHTED DISTANCE / RESISTANCE
# ============================================================
print("\n" + "=" * 72)
print("SECTION 10: WEIGHTED GRAPH DISTANCE")
print("=" * 72)

# Try: weighted distance where edge weight = product of endpoint dims
# or: weight = dim at target node
# or: weight = something involving phi

# Path weights: each edge (i,j) has weight w(i,j) = f(d_i, d_j)
print(f"  Weighted distance with various edge weights:")

# Weight 1: dim of target
def weighted_dist(start, weights_func):
    """Compute shortest weighted distance from start to all nodes."""
    import heapq
    dist = [float('inf')] * 9
    dist[start] = 0
    pq = [(0, start)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist[u]:
            continue
        for v in range(9):
            if A_mckay[u,v]:
                w = weights_func(u, v)
                if dist[u] + w < dist[v]:
                    dist[v] = dist[u] + w
                    heapq.heappush(pq, (dist[v], v))
    return dist

# Various weight functions
weight_fns = {
    "unweighted": lambda u, v: 1,
    "dim_target": lambda u, v: dims[v],
    "dim_sum": lambda u, v: dims[u] + dims[v],
    "dim_product": lambda u, v: dims[u] * dims[v],
    "a1*1+dim": lambda u, v: a1 + dims[v],
    "1/dim_target": lambda u, v: 1.0/dims[v],
}

for wname, wfn in weight_fns.items():
    dists = weighted_dist(0, wfn)
    dists_sorted = sorted(dists)
    print(f"\n  {wname}: {[f'{d:.1f}' for d in dists]}")
    print(f"   sorted: {[f'{d:.1f}' for d in dists_sorted]}")

# ============================================================
# SECTION 11: HONEST ASSESSMENT
# ============================================================
print("\n" + "=" * 72)
print("SECTION 11: HONEST ASSESSMENT")
print("=" * 72)

print(f"""
  WHAT WE FOUND:
  1. The mass coefficients (5,6) = (a1, b1) are the two largest Dynkin marks.
     This is a structural observation but not a derivation.

  2. The norm N(a1+b1*phi) = 19 = unique prime from exp342f.
     This connects mass formula to uniqueness theorem.

  3. n = a1 * z(b1/a1) where z(x) = a + bx.
     The mass exponent is a1 times the evaluation of z at b1/a1.

  4. The tensor product tower gives "first appearance" levels for each node.
     These are related to but NOT equal to the mass exponents.

  5. No simple weighted graph distance reproduces the mass exponents exactly.

  WHAT REMAINS:
  1. WHY the mass height direction is (a1, b1) and not some other vector.
     Best answer so far: (a1, b1) = (second-largest, largest) Dynkin marks.
     But this isn't a DERIVATION.

  2. WHY the mass formula is LINEAR (n = 5a+6b) rather than, say, quadratic.
     This is related to the phi^n scaling: masses are exponential in n,
     and n should be additive => linear in (a,b).

  3. The tensor product tower gives a DIFFERENT set of levels.
     The mass exponents are NOT the first-appearance levels.

  CATEGORY: PARTIALLY DERIVED
  - The fact that (a1,b1) = top Dynkin marks: STRUCTURAL OBSERVATION
  - N(a1+b1*phi) = 19 = unique prime: DERIVED CONNECTION
  - Full derivation of 5a+6b from McKay: OPEN
""")

elapsed = time.time() - t0
print(f"\nCompleted in {elapsed:.1f}s")

"""
exp342d: Exact Analysis of Character Formula Weights
=====================================================
The mass exponent n_k = sum_c w_c * chi_k(c) is EXACT with 9 weights.
Analyze whether these weights are framework constants.

Key result from exp342c: w(3a) = 4 EXACT, w(6a) = -2 EXACT, w(4a) = -13/4.
"""
import numpy as np
from numpy.linalg import eigh

PHI = (1 + np.sqrt(5)) / 2
PHIp = (1 - np.sqrt(5)) / 2
A1 = 5
B1 = 6
N_GEN = 3
N_EIG = 9
G_ORDER = 120

# Setup
dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3], dtype=float)
A = np.zeros((9, 9))
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
for i, j in edges:
    A[i,j] = 1
    A[j,i] = 1

evals_raw, evecs_raw = eigh(A)
idx = np.argsort(evals_raw)[::-1]
evals = evals_raw[idx]
evecs = evecs_raw[:, idx]
for j in range(9):
    if evecs[0, j] < 0:
        evecs[:, j] *= -1

# Character table
class_names = ['1a', '10a', '3a', '10b', '4a', '5b', '6a', '5a', '2a']
class_sizes = [1, 12, 20, 12, 30, 12, 20, 12, 1]

char_table = np.zeros((9, 9))
for c in range(9):
    norm = np.sqrt(G_ORDER / class_sizes[c])
    for k in range(9):
        char_table[k, c] = evecs[k, c] * norm

# Fermion mapping (Solution 3)
fermions = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
n_bare = np.array([0, 3, 5, 11, 11, 16, 17, 26, 19], dtype=float)
a_qn = np.array([0, 3, 1, 1, 1, 2, 1, 4, -1], dtype=float)
b_qn = np.array([0, -2, 0, 1, 1, 1, 2, 1, 4], dtype=float)

# ================================================================
# PART 46: Exact Weight Computation
# ================================================================
print("=" * 72)
print("PART 46: Exact Weight Computation w_c = (h_c/|G|)*sum_k n_k*chi_k(c)")
print("=" * 72)

# Compute numerators (before division by |G|)
print(f"\n  Numerator S_c = sum_k n_k * chi_k(c):")
for c in range(9):
    S = sum(n_bare[k] * char_table[k, c] for k in range(9))
    w = (class_sizes[c] / G_ORDER) * S
    # Check if S is in Z[phi]
    # Try S = p + q*phi
    for p in range(-500, 501):
        q = round((S - p) / PHI)
        if abs(p + q*PHI - S) < 0.001:
            Fval = A1*p + B1*q
            print(f"  {class_names[c]:>4}: S = {S:>+12.4f} = ({p:+4d}) + ({q:+4d})*phi, "
                  f"w = S*{class_sizes[c]}/{G_ORDER} = {w:>+10.6f}, "
                  f"5p+6q = {Fval}")
            break
    else:
        print(f"  {class_names[c]:>4}: S = {S:>+12.4f}, w = {w:>+10.6f}")

# ================================================================
# PART 47: Same for Quantum Numbers a and b separately
# ================================================================
print(f"\n{'='*72}")
print("PART 47: Character Expansion of a(rho) and b(rho) Separately")
print("=" * 72)

for qn_name, qn_vals in [('a', a_qn), ('b', b_qn)]:
    print(f"\n  --- Quantum number {qn_name} ---")
    print(f"  Values: {[int(v) for v in qn_vals]}")
    for c in range(9):
        S = sum(qn_vals[k] * char_table[k, c] for k in range(9))
        w = (class_sizes[c] / G_ORDER) * S
        for p in range(-200, 201):
            q = round((S - p) / PHI)
            if abs(p + q*PHI - S) < 0.001:
                print(f"  {class_names[c]:>4}: S_{qn_name} = {S:>+10.4f} = "
                      f"({p:+3d})+({q:+3d})*phi, "
                      f"w = {w:>+10.6f}")
                break

# ================================================================
# PART 48: Structure of the Numerators
# ================================================================
print(f"\n{'='*72}")
print("PART 48: Numerator Structure")
print("=" * 72)

print(f"\n  S_c = sum_k n_k * chi_k(c) for mass function n:")
sums_n = {}
for c in range(9):
    S = sum(n_bare[k] * char_table[k, c] for k in range(9))
    sums_n[class_names[c]] = S

# Rational class sums
print(f"\n  RATIONAL CLASSES:")
print(f"  S(1a)  = sum(n_k * d_k)     = {sums_n['1a']:.0f}")
print(f"  S(2a)  = sum(n_k * +/-d_k)  = {sums_n['2a']:.0f}")
print(f"  S(3a)  = sum(n_k * chi(3a)) = {sums_n['3a']:.0f}")
print(f"  S(4a)  = sum(n_k * chi(4a)) = {sums_n['4a']:.0f}")
print(f"  S(6a)  = sum(n_k * chi(6a)) = {sums_n['6a']:.0f}")

# Decompose S(1a)
S_A5 = sum(n_bare[k] * dims[k] for k in range(9) if char_table[k, 8] > 0)
S_spin = sum(n_bare[k] * dims[k] for k in range(9) if char_table[k, 8] < 0)
print(f"\n  S(1a) = {sums_n['1a']:.0f} = S_A5 + S_spin = {S_A5:.0f} + {S_spin:.0f}")
print(f"  S(2a) = {sums_n['2a']:.0f} = S_A5 - S_spin = {S_A5:.0f} - {S_spin:.0f}")
print(f"\n  S_A5 = {S_A5:.0f}, S_spin = {S_spin:.0f}")
print(f"  Ratio: {S_A5/S_spin:.6f}")

# Framework identification
print(f"\n  FRAMEWORK IDENTIFICATION:")
print(f"  S(3a)  = {sums_n['3a']:.0f} = |2T| = 24? YES" if abs(sums_n['3a']-24)<0.1 else
      f"  S(3a)  = {sums_n['3a']:.0f}")
print(f"  S(6a)  = {sums_n['6a']:.0f} = -12 = -vertex_degree? "
      f"{'YES' if abs(sums_n['6a']+12)<0.1 else 'NO'}")
print(f"  S(4a)  = {sums_n['4a']:.0f} = -13 = -(a1^2+1)/2? "
      f"{'YES' if abs(sums_n['4a']+13)<0.1 else 'NO'}")
print(f"  S(1a)  = {sums_n['1a']:.0f}")
print(f"  S(2a)  = {sums_n['2a']:.0f}")

# Irrational class sums
print(f"\n  IRRATIONAL CLASSES (Z[phi] decomposition p + q*phi):")
for cname in ['10a', '10b', '5a', '5b']:
    S = sums_n[cname]
    for p in range(-500, 501):
        q = round((S - p) / PHI)
        if abs(p + q*PHI - S) < 0.001:
            print(f"  S({cname:>3}) = {S:>+10.4f} = ({p:+4d}) + ({q:+4d})*phi")
            break

# ================================================================
# PART 49: Galois Pairs in the Weights
# ================================================================
print(f"\n{'='*72}")
print("PART 49: Galois Structure of Weights")
print("=" * 72)

# The weights come in Galois pairs:
# w(10a) <-> w(5b) under sigma: phi -> phi'
# w(10b) <-> w(5a) under sigma: phi -> phi'

# This is because sigma swaps classes 10a<->5b and 10b<->5a
# (the elements of order 10 and their negatives of order 5)

# For the sums:
print(f"\n  Galois pairs of numerator sums:")
for name1, name2 in [('10a','5b'), ('10b','5a')]:
    S1 = sums_n[name1]
    S2 = sums_n[name2]
    print(f"  S({name1}) + S({name2}) = {S1:.4f} + {S2:.4f} = {S1+S2:.4f}")
    print(f"  S({name1}) - S({name2}) = {S1:.4f} - {S2:.4f} = {S1-S2:.4f}")
    # Decompose sum and difference
    sm = S1 + S2
    df = S1 - S2
    for p in range(-200, 201):
        q = round((sm - p) / PHI)
        if abs(p + q*PHI - sm) < 0.001:
            print(f"    Sum = ({p:+d}) + ({q:+d})*phi", end="")
            if q == 0:
                print(f" = {p} (RATIONAL!)")
            else:
                print()
            break
    for p in range(-200, 201):
        q = round((df - p) / PHI)
        if abs(p + q*PHI - df) < 0.001:
            print(f"    Diff = ({p:+d}) + ({q:+d})*phi")
            break

# ================================================================
# PART 50: The Class Function as Spectral Action
# ================================================================
print(f"\n{'='*72}")
print("PART 50: Weights as Spectral Action Coefficients")
print("=" * 72)

# The mass exponent defines a class function F on 2I:
# F(g) = w_c(g)  for g in class c
# This function has specific values on each class:

print(f"\n  Mass class function F(g):")
for c in range(9):
    w = (class_sizes[c] / G_ORDER) * sums_n[class_names[c]]
    order_g = [1, 10, 3, 10, 4, 5, 6, 5, 2][c]
    print(f"  Class {class_names[c]:>4} (order {order_g:>2}, size {class_sizes[c]:>2}): "
          f"F = {w:>+10.6f}")

# The spectral action connection:
# In NCG, Tr(f(D/Lambda)) = sum_k f(lambda_k/Lambda)
# Here, if D is the Cayley Laplacian with gen set S:
# Tr(rho(f)) = sum_k f(lambda_k) * dim(rho_k)
# For our mass function: we need f such that f(lambda_k) = n_k
# This is possible since we have 9 distinct lambda_k values.

# On the Cayley graph with gen set S = all of 2I \ {e}:
# lambda_k = |G|-1 - (sum_{g!=e} chi_k(g))/dim_k
# = 119 - (sum_c |class_c|*chi_k(c))/dim_k + chi_k(e)/dim_k
# Actually, for the regular representation...

# Let's try a simpler approach: Cayley Laplacian with gen = 5-cycles (both classes)
print(f"\n  Cayley Laplacian (gen = all 5-cycles, |S|=24):")
print(f"  lambda_k = 24 - 12*(chi(5a)+chi(5b))/dim")
for k in range(9):
    c5a = char_table[k, 7]  # column 5a
    c5b = char_table[k, 5]  # column 5b
    lam = 24 - 12*(c5a + c5b)/dims[k]
    print(f"  {fermions[k]:>4} (rho_{k+1}, dim {int(dims[k])}): "
          f"chi(5a)+chi(5b)={c5a+c5b:>+8.4f}, lambda={lam:>10.4f}, n={int(n_bare[k])}")

# ================================================================
# PART 51: The Minimal Framework Formula
# ================================================================
print(f"\n{'='*72}")
print("PART 51: Can a 2-3 Parameter Formula Work?")
print("=" * 72)

# Test n_k = alpha * dim_k + beta * (dim_k mod 2) + gamma
X = np.column_stack([dims, dims % 2, np.ones(9)])
result = np.linalg.lstsq(X, n_bare, rcond=None)
pred = X @ result[0]
rms = np.sqrt(np.mean((n_bare - pred)**2))
print(f"\n  n = {result[0][0]:.3f}*dim + {result[0][1]:.3f}*is_odd + {result[0][2]:.3f}")
print(f"  RMS = {rms:.3f}")

# Test n_k = alpha * graph_distance + beta * dim_k
from collections import deque
def bfs(adj, s):
    d = [-1]*adj.shape[0]
    d[s] = 0
    q = deque([s])
    while q:
        u = q.popleft()
        for v in range(adj.shape[0]):
            if adj[u,v] > 0 and d[v] < 0:
                d[v] = d[u]+1
                q.append(v)
    return d
gdist = np.array(bfs(A, 0), dtype=float)

X2 = np.column_stack([gdist, dims])
result2 = np.linalg.lstsq(X2, n_bare, rcond=None)
pred2 = X2 @ result2[0]
rms2 = np.sqrt(np.mean((n_bare - pred2)**2))
print(f"\n  n = {result2[0][0]:.3f}*dist + {result2[0][1]:.3f}*dim")
print(f"  RMS = {rms2:.3f}")
for k in range(9):
    print(f"    {fermions[k]:>4}: pred={pred2[k]:.1f}, actual={int(n_bare[k])}, "
          f"err={n_bare[k]-pred2[k]:+.1f}")

# Test n_k = alpha * gdist + beta * cumsum_dim
def bfs_path(adj, start):
    n = adj.shape[0]
    paths = {start: [start]}
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in range(n):
            if adj[u,v] > 0 and v not in paths:
                paths[v] = paths[u] + [v]
                queue.append(v)
    return paths
paths = bfs_path(A, 0)
cumsum = np.array([sum(dims[node] for node in paths[k][1:]) for k in range(9)])

X3 = np.column_stack([gdist, cumsum])
result3 = np.linalg.lstsq(X3, n_bare, rcond=None)
pred3 = X3 @ result3[0]
rms3 = np.sqrt(np.mean((n_bare - pred3)**2))
print(f"\n  n = {result3[0][0]:.4f}*dist + {result3[0][1]:.4f}*cumsum")
print(f"  RMS = {rms3:.3f}")
for k in range(9):
    print(f"    {fermions[k]:>4}: pred={pred3[k]:.1f}, actual={int(n_bare[k])}, "
          f"err={n_bare[k]-pred3[k]:+.1f}")

# ================================================================
# PART 52: The CRUCIAL Test - Cayley eigenvalue as single predictor
# ================================================================
print(f"\n{'='*72}")
print("PART 52: n vs Cayley Eigenvalue lambda(10a) - Detailed")
print("=" * 72)

lam_10a = np.array([12*(1 - char_table[k,1]/dims[k]) for k in range(9)])

# Sort by lambda
order = np.argsort(lam_10a)
print(f"\n  Sorted by lambda(10a):")
print(f"  {'Fermion':>8} {'lambda':>10} {'n':>6} {'Delta_lam':>10} {'Delta_n':>8}")
prev_lam = 0
prev_n = 0
for idx_k in order:
    dl = lam_10a[idx_k] - prev_lam
    dn = n_bare[idx_k] - prev_n
    print(f"  {fermions[idx_k]:>8} {lam_10a[idx_k]:>10.4f} {int(n_bare[idx_k]):>6} "
          f"{dl:>10.4f} {dn:>+8.0f}")
    prev_lam = lam_10a[idx_k]
    prev_n = n_bare[idx_k]

# The ratios Delta_n / Delta_lambda
print(f"\n  Ratios Delta_n / Delta_lambda:")
prev_lam = 0
prev_n = 0
for idx_k in order:
    dl = lam_10a[idx_k] - prev_lam
    dn = n_bare[idx_k] - prev_n
    if dl > 0.01:
        ratio = dn / dl
        print(f"  {fermions[idx_k]:>8}: {dn:+.0f} / {dl:.4f} = {ratio:.4f}")
    prev_lam = lam_10a[idx_k]
    prev_n = n_bare[idx_k]

# Monotone interpolation using lambda
# Check if lambda ordering preserves n ordering
print(f"\n  Lambda ordering vs n ordering:")
lam_order = [fermions[i] for i in np.argsort(lam_10a)]
n_order = [fermions[i] for i in np.argsort(n_bare)]
print(f"  By lambda: {lam_order}")
print(f"  By n:      {n_order}")
print(f"  Same order: {lam_order == n_order}")

# ================================================================
# PART 53: Key Observation - Sums Over Symmetry Classes
# ================================================================
print(f"\n{'='*72}")
print("PART 53: Sum Rules for Mass Exponents")
print("=" * 72)

print(f"\n  Total: sum(n_k) = {sum(n_bare):.0f}")
print(f"  = 108 = 4*27 = 4*(a1^2+2) ? No. 108 = 9*12 = N_eig * degree")
print(f"  = {N_EIG} * 12 = N_eig * vertex_degree: "
      f"{'YES' if abs(sum(n_bare) - N_EIG*12) < 0.1 else 'NO'}")

print(f"\n  Weighted sums:")
print(f"  sum(n_k * d_k) = {sum(n_bare[k]*dims[k] for k in range(9)):.0f}")
print(f"    = 393 = 3 * 131")

print(f"\n  By A5/spin type:")
A5_nodes = [k for k in range(9) if char_table[k,8] > 0]
spin_nodes = [k for k in range(9) if char_table[k,8] < 0]
print(f"  A5 reps (nodes {A5_nodes}): "
      f"fermions = {[fermions[k] for k in A5_nodes]}")
print(f"  spin reps (nodes {spin_nodes}): "
      f"fermions = {[fermions[k] for k in spin_nodes]}")

sum_A5 = sum(n_bare[k] for k in A5_nodes)
sum_spin = sum(n_bare[k] for k in spin_nodes)
print(f"  sum_A5 = {sum_A5:.0f} (e+d+mu+tau+b = 0+5+11+17+19 = 52)")
print(f"  sum_spin = {sum_spin:.0f} (u+s+c+t = 3+11+16+26 = 56)")

print(f"\n  By Galois pairs (same dim at opposite ends):")
pairs = [(1,7), (2,8), (3,6)]  # (u,t), (d,b), (s,tau)
for i, j in pairs:
    print(f"  rho_{i+1},{j+1} ({fermions[i]},{fermions[j]}): "
          f"n = {int(n_bare[i])}+{int(n_bare[j])} = {int(n_bare[i]+n_bare[j])}, "
          f"dim = {int(dims[i])}={int(dims[j])}")
print(f"  Unpaired: rho_1({fermions[0]}, n={int(n_bare[0])}), "
      f"rho_5({fermions[4]}, n={int(n_bare[4])}), "
      f"rho_6({fermions[5]}, n={int(n_bare[5])})")

# Sum of Galois pair exponents
pair_sums = [int(n_bare[i]+n_bare[j]) for i,j in pairs]
print(f"\n  Pair sums: {pair_sums}")
print(f"  Mean: {np.mean(pair_sums):.2f}")
print(f"  The pair sums are {pair_sums[0]}, {pair_sums[1]}, {pair_sums[2]}")
print(f"  = 29, 24, 28")
print(f"  Sum of pair sums: {sum(pair_sums)} = 81 = 3^4")
print(f"  Product: {pair_sums[0]*pair_sums[1]*pair_sums[2]}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*72}")
print("SUMMARY: Character Formula Weight Analysis")
print("=" * 72)

print(f"""
The mass exponent n_k = sum_c w_c * chi_k(c) is EXACT with weights:

  RATIONAL CLASS WEIGHTS (exact integers/rationals):
    w(1a)  = 393/120 = 131/40
    w(2a)  = -3/120 = -1/40
    w(3a)  = 24*20/120 = 480/120 = 4         <-- |2T|/h_c * h_c/|G| = ...
    w(4a)  = -13*30/120 = -390/120 = -13/4   <-- -(a1^2+1)/4
    w(6a)  = -12*20/120 = -240/120 = -2      <-- -vertex_degree/b1 = -12/6

  The NUMERATOR SUMS are key:
    S(3a) = 24 = |2T| = 2 * 12 = 2 * vertex_degree  <-- FRAMEWORK
    S(6a) = -12 = -vertex_degree = -2*(a1+1)         <-- FRAMEWORK
    S(4a) = -13 = -(a1^2+1)/2                        <-- FRAMEWORK

  TOTAL: sum(n_k) = 108 = N_eig * vertex_degree = 9 * 12 <-- FRAMEWORK

  IRRATIONAL CLASS WEIGHTS (Z[phi]):
    S(10a) = 23 - 37*phi    (37 = b1^2 + 1)
    S(5b)  = -14 + 37*phi   (Galois conjugate!)
    S(10a) + S(5b) = 9 = N_eig            <-- FRAMEWORK
    S(10b) = -10 - 9*phi
    S(5a)  = -19 + 9*phi    (Galois conjugate!)
    S(10b) + S(5a) = -29 = -(h(E8) - 1)  <-- FRAMEWORK?

  GALOIS PAIRS on the graph:
    (u,t): n_u + n_t = 29, dims (2,2)
    (d,b): n_d + n_b = 24, dims (3,3)
    (s,tau): n_s + n_tau = 28, dims (4,4)
""")
print("Done.")

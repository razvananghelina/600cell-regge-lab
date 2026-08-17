"""
exp342c: Cayley Eigenvalues vs Mass Exponents
==============================================
Uses the fermion-node mapping from exp323 (Solution 3) and the
2I character table from exp342b to test whether Cayley Laplacian
eigenvalues determine mass exponents.

Key hypothesis: lambda_k(10a) = 12(1 - chi_k(10a)/dim_k) is related
to the mass exponent n_k by a simple function.
"""
import numpy as np
from numpy.linalg import eigh

PHI = (1 + np.sqrt(5)) / 2
PHIp = (1 - np.sqrt(5)) / 2
A1 = 5
B1 = 6

# ======================
# Setup: Character table
# ======================
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

# Eigenvalue -> class: 2->1a, phi->10a, 1->3a, 1/phi->10b, 0->4a,
#                       -1/phi->5b, -1->6a, -phi->5a, -2->2a
class_map = {0:'1a', 1:'10a', 2:'3a', 3:'10b', 4:'4a',
             5:'5b', 6:'6a', 7:'5a', 8:'2a'}
class_sizes = {0:1, 1:12, 2:20, 3:12, 4:30, 5:12, 6:20, 7:12, 8:1}

char_table = np.zeros((9, 9))
for c in range(9):
    norm = np.sqrt(120.0 / class_sizes[c])
    for k in range(9):
        char_table[k, c] = evecs[k, c] * norm

# ======================
# Fermion-node mapping (Solution 3 from exp323)
# ======================
# node: (fermion, mass_exponent, (a,b))
fermion_map = {
    0: ('e',   0,  (0, 0)),
    1: ('u',   3,  (3,-2)),
    2: ('d',   5,  (1, 0)),
    3: ('s',  11,  (1, 1)),    # mu/s degenerate, s here
    4: ('mu', 11,  (1, 1)),    # mu/s degenerate, mu here
    5: ('c',  16,  (2, 1)),
    6: ('tau',17,  (1, 2)),
    7: ('t',  26,  (4, 1)),
    8: ('b',  19, (-1, 4)),
}

fermions = [fermion_map[k][0] for k in range(9)]
n_bare = [fermion_map[k][1] for k in range(9)]
ab_qn = [fermion_map[k][2] for k in range(9)]

# ================================================================
# PART 34: Cayley Eigenvalues for Class 10a
# ================================================================
print("=" * 72)
print("PART 34: Cayley Laplacian Eigenvalues (gen set = class 10a)")
print("=" * 72)

# Column index for class 10a (eigenvalue = phi)
col_10a = 1  # second eigenvalue is phi

print(f"\n  {'Node':>6} {'Fermion':>8} {'dim':>4} {'chi(10a)':>10} "
      f"{'lambda':>10} {'n_bare':>8} {'n/lambda':>10}")
print("-" * 72)

lam_10a = []
for k in range(9):
    d = dims[k]
    chi = char_table[k, col_10a]
    lam = 12 * (1 - chi/d)
    lam_10a.append(lam)
    ratio = n_bare[k]/lam if abs(lam) > 0.01 else float('inf')
    print(f"  rho_{k+1:1d}  {fermions[k]:>8} {int(d):>4} {chi:>+10.4f} "
          f"{lam:>10.4f} {n_bare[k]:>8} {ratio:>10.3f}")

# ================================================================
# PART 35: Cayley Eigenvalues for Class 5a
# ================================================================
print(f"\n{'='*72}")
print("PART 35: Cayley Laplacian Eigenvalues (gen set = class 5a)")
print("=" * 72)

col_5a = 7  # eigenvalue = -phi

print(f"\n  {'Node':>6} {'Fermion':>8} {'dim':>4} {'chi(5a)':>10} "
      f"{'lambda':>10} {'n_bare':>8}")
print("-" * 72)

lam_5a = []
for k in range(9):
    d = dims[k]
    chi = char_table[k, col_5a]
    lam = 12 * (1 - chi/d)
    lam_5a.append(lam)
    print(f"  rho_{k+1:1d}  {fermions[k]:>8} {int(d):>4} {chi:>+10.4f} "
          f"{lam:>10.4f} {n_bare[k]:>8}")

# ================================================================
# PART 36: Combined: n = alpha*lam_10a + beta*lam_5a + gamma
# ================================================================
print(f"\n{'='*72}")
print("PART 36: Fit n = alpha*lam(10a) + beta*lam(5a) + gamma")
print("=" * 72)

X = np.column_stack([lam_10a, lam_5a, np.ones(9)])
n_arr = np.array(n_bare, dtype=float)

# Least squares
result = np.linalg.lstsq(X, n_arr, rcond=None)
params = result[0]
pred = X @ params
resid = n_arr - pred
rms = np.sqrt(np.mean(resid**2))

print(f"\n  alpha = {params[0]:.6f}")
print(f"  beta  = {params[1]:.6f}")
print(f"  gamma = {params[2]:.6f}")
print(f"  RMS   = {rms:.4f}")
print(f"\n  {'Fermion':>8} {'n_bare':>8} {'predicted':>10} {'residual':>10}")
for k in range(9):
    print(f"  {fermions[k]:>8} {n_bare[k]:>8} {pred[k]:>10.3f} {resid[k]:>+10.3f}")

# ================================================================
# PART 37: Try n = alpha*lam(10a) + beta*lam(5a)  (no constant)
# ================================================================
print(f"\n{'='*72}")
print("PART 37: Fit n = alpha*lam(10a) + beta*lam(5a) (zero-intercept)")
print("=" * 72)

X2 = np.column_stack([lam_10a, lam_5a])
result2 = np.linalg.lstsq(X2, n_arr, rcond=None)
params2 = result2[0]
pred2 = X2 @ params2
resid2 = n_arr - pred2
rms2 = np.sqrt(np.mean(resid2**2))

print(f"\n  alpha = {params2[0]:.6f}")
print(f"  beta  = {params2[1]:.6f}")
print(f"  RMS   = {rms2:.4f}")
print(f"\n  {'Fermion':>8} {'n_bare':>8} {'predicted':>10} {'residual':>10}")
for k in range(9):
    print(f"  {fermions[k]:>8} {n_bare[k]:>8} {pred2[k]:>10.3f} {resid2[k]:>+10.3f}")

# Check if alpha, beta are framework numbers
print(f"\n  alpha/beta = {params2[0]/params2[1]:.6f}")
print(f"  a1/b1 = {A1/B1:.6f}")
print(f"  phi = {PHI:.6f}")

# ================================================================
# PART 38: Use ALL 9 Cayley eigenvalues
# ================================================================
print(f"\n{'='*72}")
print("PART 38: All 9 Cayley Eigenvalues per Node")
print("=" * 72)

# For each class c, compute lambda_k(c) = |class_c| * (1 - chi_k(c)/dim_k)
class_cols = list(range(9))
class_names = ['1a', '10a', '3a', '10b', '4a', '5b', '6a', '5a', '2a']
class_sz = [1, 12, 20, 12, 30, 12, 20, 12, 1]

all_lam = np.zeros((9, 9))
for c in range(9):
    for k in range(9):
        all_lam[k, c] = class_sz[c] * (1 - char_table[k, c] / dims[k])

# Fit n = sum_c w_c * lambda_k(c)
result3 = np.linalg.lstsq(all_lam, n_arr, rcond=None)
params3 = result3[0]
pred3 = all_lam @ params3
resid3 = n_arr - pred3
rms3 = np.sqrt(np.mean(resid3**2))

print(f"\n  Weights for 9-class fit:")
for c in range(9):
    print(f"    w({class_names[c]:>4}) = {params3[c]:>+10.6f}")
print(f"  RMS = {rms3:.6f}")
print(f"\n  {'Fermion':>8} {'n_bare':>8} {'predicted':>10} {'residual':>10}")
for k in range(9):
    print(f"  {fermions[k]:>8} {n_bare[k]:>8} {pred3[k]:>10.4f} {resid3[k]:>+10.4f}")

# ================================================================
# PART 39: Character-based Mass Formula
# ================================================================
print(f"\n{'='*72}")
print("PART 39: Direct Character-Based Mass Formula")
print("=" * 72)

# Test: n_k = sum_c w_c * chi_k(c) for weight vector w
X4 = char_table.copy()
result4 = np.linalg.lstsq(X4, n_arr, rcond=None)
params4 = result4[0]
pred4 = X4 @ params4
resid4 = n_arr - pred4
rms4 = np.sqrt(np.mean(resid4**2))

print(f"\n  n_k = sum_c w_c * chi_k(c)")
print(f"  Weights:")
for c in range(9):
    print(f"    w({class_names[c]:>4}) = {params4[c]:>+10.6f}")
print(f"  RMS = {rms4:.6f}")
print(f"\n  {'Fermion':>8} {'n_bare':>8} {'predicted':>10} {'residual':>10}")
for k in range(9):
    print(f"  {fermions[k]:>8} {n_bare[k]:>8} {pred4[k]:>10.4f} {resid4[k]:>+10.4f}")

# Check if weights are framework numbers
print(f"\n  Checking if weights are framework constants:")
print(f"  w(1a) * 120 = {params4[0]*120:.4f}")
print(f"  w(10a) * 12 = {params4[1]*12:.4f}")
print(f"  w(5a) * 12 = {params4[7]*12:.4f}")

# ================================================================
# PART 40: The KEY - Z[phi] Quantum Numbers from Characters
# ================================================================
print(f"\n{'='*72}")
print("PART 40: Can (a,b) Be Read from Characters?")
print("=" * 72)

# Z[phi] decomposition of chi(5a): chi = p + q*phi
def to_zphi(x):
    best = None
    for p in range(-10, 11):
        q = round((x - p) / PHI)
        if abs(p + q*PHI - x) < 0.001:
            if best is None or abs(p)+abs(q) < abs(best[0])+abs(best[1]):
                best = (p, q)
    return best if best else (round(x), 0)

print(f"\n  Comparison of Z[phi] coordinates:")
print(f"  {'Node':>6} {'Fermion':>8} {'n':>4} {'(a,b)_mass':>12} "
      f"{'chi(5a)':>10} {'(p,q)_chi':>12} {'chi(10a)':>10} {'(p,q)_10a':>12}")
print("-" * 90)

col_10a_idx = 1
col_5a_idx = 7

for k in range(9):
    a, b = ab_qn[k]
    c5a = char_table[k, col_5a_idx]
    c10a = char_table[k, col_10a_idx]
    p5, q5 = to_zphi(c5a)
    p10, q10 = to_zphi(c10a)
    print(f"  rho_{k+1:1d}  {fermions[k]:>8} {n_bare[k]:>4}  ({a:+2d},{b:+2d})   "
          f"{c5a:>+10.4f}  ({p5:+2d},{q5:+2d})   "
          f"{c10a:>+10.4f}  ({p10:+2d},{q10:+2d})")

# ================================================================
# PART 41: The Cayley Laplacian in Z[phi]
# ================================================================
print(f"\n{'='*72}")
print("PART 41: Cayley Eigenvalues in Z[phi]")
print("=" * 72)

print(f"\n  lambda(10a) = 12 - 12*chi(10a)/dim  in Z[phi]:")
print(f"  {'Fermion':>8} {'lambda':>10} {'p+q*phi':>15} {'5p+6q':>8} {'n_bare':>8}")
print("-" * 60)

for k in range(9):
    lam = lam_10a[k]
    p, q = to_zphi(lam)
    F = A1*p + B1*q
    print(f"  {fermions[k]:>8} {lam:>10.4f}   ({p:+3d})+({q:+3d})*phi "
          f"{F:>8} {n_bare[k]:>8}")

# ================================================================
# PART 42: Dimension * (normalized chi difference)
# ================================================================
print(f"\n{'='*72}")
print("PART 42: Dimension-Weighted Character Differences")
print("=" * 72)

print(f"\n  Test: d_k * (chi_k(5a) - chi_k(5b))")
for k in range(9):
    d = dims[k]
    diff = char_table[k, col_5a_idx] - char_table[k, 5]  # 5b is col 5
    weighted = d * diff
    print(f"  {fermions[k]:>8}: dim={int(d)}, chi(5a)-chi(5b)={diff:+.4f}, "
          f"d*(diff)={weighted:+.4f}")

print(f"\n  Test: d_k^2 * (chi_k(5a) - chi_k(5b))^2 / 5")
for k in range(9):
    d = dims[k]
    diff = char_table[k, col_5a_idx] - char_table[k, 5]
    val = d**2 * diff**2 / A1
    print(f"  {fermions[k]:>8}: val = {val:.4f}  (n = {n_bare[k]})")

# ================================================================
# PART 43: Path-based Mass with Fermion Mapping
# ================================================================
print(f"\n{'='*72}")
print("PART 43: Path Metrics with Known Fermion Mapping")
print("=" * 72)

from collections import deque
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

print(f"\n  Graph distance and mass exponent:")
print(f"  {'Fermion':>8} {'dist':>6} {'n_bare':>8} {'cumsum_d':>10} {'log_phi(prod)':>14}")
for k in range(9):
    path = paths[k]
    dist = len(path) - 1
    cumsum = sum(dims[node] for node in path[1:])
    prod = 1
    for node in path:
        prod *= int(dims[node])
    logp = np.log(prod)/np.log(PHI) if prod > 1 else 0
    print(f"  {fermions[k]:>8} {dist:>6} {n_bare[k]:>8} {cumsum:>10.0f} {logp:>14.3f}")

# Best linear combination of metrics
print(f"\n  Fit: n = alpha*dist + beta*cumsum + gamma*log_phi(prod)")
dist_arr = np.array([len(paths[k])-1 for k in range(9)], dtype=float)
cumsum_arr = np.array([sum(dims[node] for node in paths[k][1:]) for k in range(9)])
prod_arr = np.array([np.log(max(1, np.prod([dims[node] for node in paths[k]])))/np.log(PHI)
                     for k in range(9)])

X5 = np.column_stack([dist_arr, cumsum_arr, prod_arr])
result5 = np.linalg.lstsq(X5, n_arr, rcond=None)
params5 = result5[0]
pred5 = X5 @ params5
rms5 = np.sqrt(np.mean((n_arr - pred5)**2))

print(f"  alpha={params5[0]:.4f}, beta={params5[1]:.4f}, gamma={params5[2]:.4f}")
print(f"  RMS = {rms5:.4f}")
for k in range(9):
    print(f"  {fermions[k]:>8}: pred={pred5[k]:>8.2f}, actual={n_bare[k]:>4}, "
          f"diff={n_bare[k]-pred5[k]:>+6.2f}")

# ================================================================
# PART 44: The SPECTRAL Mass Operator
# ================================================================
print(f"\n{'='*72}")
print("PART 44: Spectral Mass Operator - Exact Solution")
print("=" * 72)

# If M is a 9x9 matrix on the rep ring with M|k> = n_k|k>,
# can we write M as a polynomial in A (the adjacency matrix)?
# M = sum c_p * A^p

# Since A has 9 distinct eigenvalues, any function f(A) is determined
# by f(lambda_k) for each eigenvalue. We need:
# f(lambda_k) = n_k for each k.

# But the EIGENVECTORS of A are the character table columns (classes),
# not the fermion basis vectors (irreps). The fermion at node k is NOT
# an eigenvector of A.

# Instead, the mass operator in the NODE basis is simply:
M_node = np.diag(n_arr)
print(f"\n  M = diag(n_bare) in the node basis")
print(f"  M expressed as polynomial in A:")

# Express M as sum c_p * A^p for p=0,...,8
A_powers = [np.eye(9)]
for p in range(1, 9):
    A_powers.append(A_powers[-1] @ A)

# Solve: sum_p c_p * A^p = M_node
# Vectorize: reshape each A^p as 81-dim vector
nent = 81
A_mat = np.zeros((nent, 9))
for p in range(9):
    A_mat[:, p] = A_powers[p].flatten()
M_vec = M_node.flatten()

result6 = np.linalg.lstsq(A_mat, M_vec, rcond=None)
coeffs = result6[0]
M_recon = sum(coeffs[p] * A_powers[p] for p in range(9))
err = np.max(np.abs(M_recon - M_node))

print(f"  M = ", end="")
terms = []
for p in range(9):
    if abs(coeffs[p]) > 1e-6:
        terms.append(f"{coeffs[p]:+.4f}*A^{p}")
print(" + ".join(terms))
print(f"  Max reconstruction error: {err:.2e}")

# Check if coefficients are Z[phi] numbers
print(f"\n  Coefficients in Z[phi]?")
for p in range(9):
    if abs(coeffs[p]) > 1e-6:
        pz, qz = to_zphi(coeffs[p])
        actual = pz + qz * PHI
        print(f"    c_{p} = {coeffs[p]:>+10.6f} ~ ({pz:+d})+({qz:+d})*phi = "
              f"{actual:.6f} (err={abs(coeffs[p]-actual):.4e})")

# ================================================================
# PART 45: Character Inner Product with Mass
# ================================================================
print(f"\n{'='*72}")
print("PART 45: Mass Exponent as Class Function Inner Product")
print("=" * 72)

# If n_k is a class function on 2I (function of irreps), we can
# decompose it in the character basis:
# n_k = sum_c w_c * chi_k(c)
# This was done in Part 39. Let's interpret the weights.

# The inner product on class functions is:
# <f, g> = (1/|G|) sum_{g in G} f(g) conj(g(g))
# In terms of characters: <chi_i, chi_j> = delta_ij

# So w_c = <n, chi_c_bar> = (1/|G|) sum_k n_k * chi_k(c_bar) * |class_c|/dim_k ...
# Actually for functions on irreps: w_c = sum_k n_k * chi_k(c) * class_size_c / |G|

# Actually, since {chi_k(c)} is the character table with orthogonality:
# sum_k chi_k(c) * chi_k(c') = |G|/|class_c| * delta_{c,c'}
# The expansion n_k = sum_c w_c * chi_k(c) gives:
# w_c = sum_k n_k * chi_k(c) * |class_c| / |G|

# Let's compute this directly:
print(f"\n  n_k as linear combination of characters:")
print(f"  n_k = sum_c w_c * chi_k(c) where w_c = (|class_c|/|G|) * sum_k n_k * chi_k(c)")

for c in range(9):
    sz = class_sz[c]
    w = (sz / 120.0) * sum(n_bare[k] * char_table[k, c] for k in range(9))
    print(f"    w({class_names[c]:>4}) = {w:>+10.6f}")
    # Check against Part 39 result
    if abs(w - params4[c]) > 0.01:
        print(f"      WARNING: mismatch with lstsq result {params4[c]:.6f}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*72}")
print("SUMMARY OF KEY FINDINGS")
print("=" * 72)

print(f"""
FERMION-NODE MAPPING (Solution 3):
  e->rho_1(1), u->rho_2(2), d->rho_3(3), s->rho_4(4), mu->rho_5(5)
  c->rho_6(6), tau->rho_7(4), t->rho_8(2), b->rho_9(3)

CAYLEY EIGENVALUE (class 10a) VS MASS EXPONENT:
  Roughly monotone but NOT simple proportionality.

CHARACTER-BASED FIT QUALITY:
  RMS of 2-param Cayley fit: {rms2:.4f}
  RMS of 3-param Cayley fit: {rms:.4f}
  RMS of 9-param character fit: {rms4:.6f} (exact if 9 params for 9 points)
  RMS of path metric fit: {rms5:.4f}

The character table CAN encode mass exponents (trivially, with 9 params
for 9 data points). The question is whether the WEIGHTS are framework constants.
""")

print("Done.")

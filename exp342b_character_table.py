"""
exp342b: Extract 2I Character Table from Affine E8 Eigenvectors
================================================================
Key insight: eigenvectors of the affine E8 adjacency matrix ARE
the columns of the 2I character table (up to normalization).

Eigenvalue lambda_c = chi_2(class_c), so each eigenvector corresponds
to a conjugacy class, and its components give chi_k(class_c).

Tests whether mass exponents n = 5a + 6b arise from character values
on the 5-cycle classes of 2I.
"""
import numpy as np
from numpy.linalg import eigh

PHI = (1 + np.sqrt(5)) / 2
PHIp = (1 - np.sqrt(5)) / 2  # = -1/phi
A1 = 5
B1 = 6
G_ORDER = 120  # |2I|

# Affine E8 adjacency
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

# Fix sign: v[0] > 0 (since chi_trivial = 1 > 0 for all classes)
for j in range(9):
    if evecs[0, j] < 0:
        evecs[:, j] *= -1

# ================================================================
# PART 21: Eigenvalue -> Conjugacy Class Assignment
# ================================================================
print("=" * 72)
print("PART 21: Eigenvalue -> Conjugacy Class of 2I")
print("=" * 72)

# chi_2(class) = eigenvalue of A
# Classes of 2I: 1a(1), 2a(1), 3a(20), 4a(30), 5a(12), 5b(12),
#                6a(20), 10a(12), 10b(12)
# chi_2 values: 2, -2, 1, 0, -phi, -1/phi, -1, phi, 1/phi

class_data = [
    (2.0,    '1a',  1,  1),    # identity
    (-2.0,   '2a',  1,  2),    # central element -I
    (1.0,    '3a',  20, 3),    # order 3
    (0.0,    '4a',  30, 4),    # order 4
    (-PHI,   '5a',  12, 5),    # order 5
    (-1/PHI, '5b',  12, 5),    # order 5
    (-1.0,   '6a',  20, 6),    # order 6
    (PHI,    '10a', 12, 10),   # order 10
    (1/PHI,  '10b', 12, 10),   # order 10
]

# Match each computed eigenvalue to a class
col_to_class = {}
for j in range(9):
    ev = evals[j]
    best = min(class_data, key=lambda x: abs(x[0] - ev))
    col_to_class[j] = best
    print(f"  lambda_{j} = {ev:+.6f} -> class {best[1]} "
          f"(size {best[2]}, order {best[3]})")

# Verify total: sum of class sizes = |2I|
total = sum(cd[2] for cd in class_data)
print(f"\n  Sum of class sizes: {total} (should be {G_ORDER})")

# ================================================================
# PART 22: Extract Character Table
# ================================================================
print(f"\n{'='*72}")
print("PART 22: Full 2I Character Table")
print("=" * 72)

# chi_k(class_c) = v_c[k] * sqrt(|G| / |class_c|)
char_table = np.zeros((9, 9))
class_names = []
class_sizes = []

for c in range(9):
    _, name, size, _ = col_to_class[c]
    class_names.append(name)
    class_sizes.append(size)
    norm = np.sqrt(G_ORDER / size)
    for k in range(9):
        char_table[k, c] = evecs[k, c] * norm

# Round to nice values (characters of finite groups are algebraic integers)
# For 2I, characters are in Z[phi]
def to_zphi(x):
    """Decompose x = p + q*phi, return (p, q)"""
    # x = p + q*(1+sqrt5)/2 = (2p+q)/2 + q*sqrt5/2
    # So q = round(2*(x - round(x))/sqrt(5)*2) ... tricky
    # Better: q*phi = x - p, so q = (x-p)/phi for integer p
    # Try all small p:
    best = None
    for p in range(-10, 11):
        q_float = (x - p) / PHI
        q = round(q_float)
        if abs(p + q*PHI - x) < 0.001:
            if best is None or abs(p) + abs(q) < abs(best[0]) + abs(best[1]):
                best = (p, q)
    if best:
        return best
    return (round(x), 0)

print(f"\n  {'':>10}", end="")
for c in range(9):
    print(f" {class_names[c]:>7}", end="")
print(f"  (size)")

print(f"  {'':>10}", end="")
for c in range(9):
    print(f" {class_sizes[c]:>7}", end="")
print()
print("-" * 90)

for k in range(9):
    print(f"  rho_{k+1}({int(dims[k])})", end="")
    for c in range(9):
        val = char_table[k, c]
        p, q = to_zphi(val)
        if q == 0:
            print(f" {p:>7}", end="")
        elif q == 1:
            if p == 0:
                print(f" {'phi':>7}", end="")
            else:
                print(f" {p}+phi".rjust(7), end="")
        elif q == -1:
            if p == 0:
                print(f" {'-phi':>7}", end="")
            else:
                print(f" {p}-phi".rjust(7), end="")
        else:
            print(f" {val:>7.3f}", end="")
    print()

# ================================================================
# PART 23: Verify Key Properties
# ================================================================
print(f"\n{'='*72}")
print("PART 23: Verification of Character Table Properties")
print("=" * 72)

# 1. chi_k(1a) = dim(rho_k)
col_1a = class_names.index('1a')
print(f"\n  chi_k(1a) vs dims:")
for k in range(9):
    print(f"    rho_{k+1}: chi = {char_table[k, col_1a]:.4f}, "
          f"dim = {dims[k]:.0f}, "
          f"match = {abs(char_table[k, col_1a] - dims[k]) < 0.01}")

# 2. chi_k(2a) = +-dim (A5 reps: +, spin reps: -)
col_2a = class_names.index('2a')
print(f"\n  chi_k(2a) / dim_k (should be +1 for A5 reps, -1 for spin):")
is_A5_rep = []
for k in range(9):
    ratio = char_table[k, col_2a] / dims[k]
    is_a5 = ratio > 0
    is_A5_rep.append(is_a5)
    print(f"    rho_{k+1}(dim {int(dims[k])}): chi(2a)/d = {ratio:+.4f} "
          f"-> {'A5 rep' if is_a5 else 'spin rep'}")

# 3. For A5 reps: chi(5a) = chi(10a), for spin: chi(5a) = -chi(10a)
col_5a = class_names.index('5a')
col_10a = class_names.index('10a')
print(f"\n  chi(5a) vs chi(10a) (A5: equal, spin: opposite):")
for k in range(9):
    c5a = char_table[k, col_5a]
    c10a = char_table[k, col_10a]
    if is_A5_rep[k]:
        match = abs(c5a - c10a) < 0.01
        print(f"    rho_{k+1} (A5):  chi(5a)={c5a:+.4f}, chi(10a)={c10a:+.4f}, "
              f"equal={match}")
    else:
        match = abs(c5a + c10a) < 0.01
        print(f"    rho_{k+1} (spin): chi(5a)={c5a:+.4f}, chi(10a)={c10a:+.4f}, "
              f"sum=0: {match}")

# ================================================================
# PART 24: Characters on 5-Cycle Classes - Z[phi] Decomposition
# ================================================================
print(f"\n{'='*72}")
print("PART 24: Characters on 5a and 5b Classes in Z[phi]")
print("=" * 72)

col_5a = class_names.index('5a')
col_5b = class_names.index('5b')
col_10a = class_names.index('10a')
col_10b = class_names.index('10b')

print(f"\n  {'Irrep':>10} {'dim':>4} {'chi(5a)':>10} {'chi(5b)':>10} "
      f"{'chi(10a)':>10} {'chi(10b)':>10} {'Z[phi] of chi(5a)':>18}")
print("-" * 80)

zphi_5a = []
zphi_5b = []
for k in range(9):
    c5a = char_table[k, col_5a]
    c5b = char_table[k, col_5b]
    c10a = char_table[k, col_10a]
    c10b = char_table[k, col_10b]
    p, q = to_zphi(c5a)
    zphi_5a.append((p, q))
    p2, q2 = to_zphi(c5b)
    zphi_5b.append((p2, q2))
    print(f"  rho_{k+1}({int(dims[k]):1d})  {int(dims[k]):>4} {c5a:>+10.4f} {c5b:>+10.4f} "
          f"{c10a:>+10.4f} {c10b:>+10.4f}   ({p:+d}) + ({q:+d})*phi")

# ================================================================
# PART 25: Mass Functional on Z[phi] Characters
# ================================================================
print(f"\n{'='*72}")
print("PART 25: Mass Functional F = a1*p + b1*q on Z[phi] of chi(5a)")
print("=" * 72)

print(f"\n  If chi_k(5a) = p + q*phi, then F = {A1}*p + {B1}*q:")
print(f"\n  {'Irrep':>10} {'dim':>4} {'(p,q)':>10} {'F=5p+6q':>10}")
print("-" * 50)

F_values = []
for k in range(9):
    p, q = zphi_5a[k]
    F = A1*p + B1*q
    F_values.append(F)
    print(f"  rho_{k+1}({int(dims[k]):1d})  {int(dims[k]):>4}  ({p:+d},{q:+d})   {F:>+10d}")

# Mass exponents needed
n_bare = {'e':0, 'mu':11, 'tau':17, 'u':3, 'c':16, 't':26,
          'd':5, 's':11, 'b':19}
print(f"\n  F values from character table: {sorted(F_values)}")
print(f"  Mass exponents needed:         {sorted(n_bare.values())}")
print(f"\n  Do they match for ANY permutation?")

F_sorted = sorted(F_values)
n_sorted = sorted(n_bare.values())
print(f"  Direct match: {F_sorted == n_sorted}")

# ================================================================
# PART 26: Other Linear Functionals on chi(5a), chi(5b)
# ================================================================
print(f"\n{'='*72}")
print("PART 26: Search for Linear Functional Giving Mass Exponents")
print("=" * 72)

# Try F = alpha * chi(5a) + beta * chi(5b) for various alpha, beta
# The mass exponents to reproduce: {0, 3, 5, 11, 11, 16, 17, 19, 26}
target = sorted(n_bare.values())

print(f"\n  Testing F = alpha * chi(5a) + beta * chi(5b):")
best_rms = 1e10
best_params = None

for a_try in np.arange(-20, 21, 0.5):
    for b_try in np.arange(-20, 21, 0.5):
        F_test = []
        for k in range(9):
            c5a = char_table[k, col_5a]
            c5b = char_table[k, col_5b]
            F_test.append(a_try * c5a + b_try * c5b)
        F_test_sorted = sorted(F_test)
        rms = np.sqrt(np.mean([(f-t)**2 for f, t in zip(F_test_sorted, target)]))
        if rms < best_rms:
            best_rms = rms
            best_params = (a_try, b_try)
            best_F = list(F_test)

print(f"  Best linear: alpha={best_params[0]:.1f}, beta={best_params[1]:.1f}, "
      f"RMS={best_rms:.3f}")
vals_sorted = sorted(best_F)
print(f"  Values:  {[f'{v:.1f}' for v in vals_sorted]}")
print(f"  Target:  {target}")

# Try with offset: F = alpha * chi(5a) + beta * chi(5b) + gamma
print(f"\n  Testing F = alpha * chi(5a) + beta * chi(5b) + gamma:")
best_rms2 = 1e10
best_params2 = None

for a_try in np.arange(-15, 16, 1):
    for b_try in np.arange(-15, 16, 1):
        for g_try in np.arange(-5, 30, 1):
            F_test = []
            for k in range(9):
                c5a = char_table[k, col_5a]
                c5b = char_table[k, col_5b]
                F_test.append(a_try * c5a + b_try * c5b + g_try)
            F_test_sorted = sorted(F_test)
            rms = np.sqrt(np.mean([(f-t)**2 for f,t in zip(F_test_sorted, target)]))
            if rms < best_rms2:
                best_rms2 = rms
                best_params2 = (a_try, b_try, g_try)
                best_F2 = list(F_test)

a, b, g = best_params2
print(f"  Best: alpha={a:.0f}, beta={b:.0f}, gamma={g:.0f}, RMS={best_rms2:.3f}")
vals_sorted2 = sorted(best_F2)
print(f"  Values:  {[f'{v:.1f}' for v in vals_sorted2]}")
print(f"  Target:  {target}")

# ================================================================
# PART 27: Characters on ALL Classes - Full Functional Search
# ================================================================
print(f"\n{'='*72}")
print("PART 27: Single-Class Character -> Mass Exponent")
print("=" * 72)

# For each class c, check if alpha*chi_k(c) + beta can give mass exponents
for c in range(9):
    chi_c = char_table[:, c]
    chi_sorted = sorted(chi_c)
    target_sorted = sorted(n_bare.values())

    # Linear regression: n = alpha*chi + beta
    # Using sorted values
    chi_arr = np.array(sorted(chi_c))
    tgt_arr = np.array(target_sorted)

    if np.std(chi_arr) > 0.01:
        alpha_opt = np.sum((chi_arr - np.mean(chi_arr))*(tgt_arr - np.mean(tgt_arr))) / \
                    np.sum((chi_arr - np.mean(chi_arr))**2)
        beta_opt = np.mean(tgt_arr) - alpha_opt * np.mean(chi_arr)
        pred = alpha_opt * chi_arr + beta_opt
        rms = np.sqrt(np.mean((pred - tgt_arr)**2))
        print(f"  Class {class_names[c]:>4}: alpha={alpha_opt:>+8.3f}, "
              f"beta={beta_opt:>+8.3f}, RMS={rms:.3f}")
    else:
        print(f"  Class {class_names[c]:>4}: constant character, skip")

# ================================================================
# PART 28: Cumulative Dimension Sums Along Paths
# ================================================================
print(f"\n{'='*72}")
print("PART 28: Cumulative Dimension Sums from Node 0")
print("=" * 72)

# BFS from node 0
from collections import deque

def bfs_path(adj, start):
    """Return shortest path from start to each node"""
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

print(f"\n  Cumulative dim sums (including start):")
for j in range(9):
    path = paths[j]
    cum_sum = sum(dims[node] for node in path)
    path_dims = [int(dims[node]) for node in path]
    print(f"  rho_{j+1}(dim {int(dims[j])}): path={path}, "
          f"dims={path_dims}, cumsum={cum_sum:.0f}")

print(f"\n  Cumulative dim sums (excluding start):")
cum_excl = []
for j in range(9):
    path = paths[j]
    cum_sum = sum(dims[node] for node in path[1:])  # exclude start
    cum_excl.append(cum_sum)
    print(f"  rho_{j+1}: cumsum_excl = {cum_sum:.0f}")

print(f"\n  Sorted cumulative (excl): {sorted(cum_excl)}")
print(f"  Target mass exponents:    {sorted(n_bare.values())}")

# ================================================================
# PART 29: Product of Dims Along Path
# ================================================================
print(f"\n{'='*72}")
print("PART 29: Other Path Metrics from Node 0")
print("=" * 72)

# Product of dims along path
print(f"\n  Product of dims along path:")
for j in range(9):
    path = paths[j]
    prod = 1
    for node in path:
        prod *= int(dims[node])
    log_phi = np.log(prod) / np.log(PHI) if prod > 0 else 0
    print(f"  rho_{j+1}: prod={prod:>8d}, log_phi(prod)={log_phi:>8.3f}")

# Weighted path: edge weight = sqrt(d_i * d_j)
print(f"\n  Weighted path length (weight = sqrt(d_i*d_j)):")
for j in range(9):
    path = paths[j]
    wt_len = 0
    for step in range(len(path)-1):
        u, v = path[step], path[step+1]
        wt_len += np.sqrt(dims[u] * dims[v])
    print(f"  rho_{j+1}: weighted_path = {wt_len:.4f}")

# Weighted path: edge weight = (d_i + d_j)/2
print(f"\n  Weighted path length (weight = (d_i+d_j)/2):")
for j in range(9):
    path = paths[j]
    wt_len = 0
    for step in range(len(path)-1):
        u, v = path[step], path[step+1]
        wt_len += (dims[u] + dims[v]) / 2
    print(f"  rho_{j+1}: weighted_path = {wt_len:.4f}")

# ================================================================
# PART 30: The Key Test - chi(5a)*chi(5b) and Norm Map
# ================================================================
print(f"\n{'='*72}")
print("PART 30: Norm Map N(chi(5a)) = chi(5a)*chi(5b)")
print("=" * 72)

print(f"\n  {'Irrep':>10} {'chi(5a)':>10} {'chi(5b)':>10} "
      f"{'chi*chi_G':>10} {'chi+chi_G':>10} {'(chi-chi_G)^2':>12}")
for k in range(9):
    c5a = char_table[k, col_5a]
    c5b = char_table[k, col_5b]
    prod = c5a * c5b  # Norm
    trace = c5a + c5b  # Trace
    diff2 = (c5a - c5b)**2
    print(f"  rho_{k+1}({int(dims[k]):1d})  {c5a:>+10.4f} {c5b:>+10.4f} "
          f"{prod:>+10.4f} {trace:>+10.4f} {diff2:>12.4f}")

# ================================================================
# PART 31: Mass Operator M = a1*L_5a + b1*L_5b on Rep Ring
# ================================================================
print(f"\n{'='*72}")
print("PART 31: Two-Cayley-Laplacian Mass Operator")
print("=" * 72)

# Cayley Laplacian with generating set S_5a (12 elements):
# lambda_k^(5a) = 12 - 12*chi_k(5a)/dim_k
# Similarly for S_5b
print(f"\n  Cayley eigenvalues for generating sets 5a and 5b:")
print(f"  {'Irrep':>10} {'dim':>4} {'lam_5a':>10} {'lam_5b':>10} "
      f"{'5*l_5a+6*l_5b':>14} {'scaled':>10}")

for k in range(9):
    d = dims[k]
    c5a = char_table[k, col_5a]
    c5b = char_table[k, col_5b]
    lam_5a = 12 * (1 - c5a/d)
    lam_5b = 12 * (1 - c5b/d)
    combined = A1 * lam_5a + B1 * lam_5b
    # Scale to match n_t = 26: what factor gives 26 from max combined?
    print(f"  rho_{k+1}({int(d):1d})  {int(d):>4} {lam_5a:>+10.4f} {lam_5b:>+10.4f} "
          f"{combined:>14.4f}")

# Find the actual eigenvalues and scale
lam_5a_vals = [12*(1 - char_table[k, col_5a]/dims[k]) for k in range(9)]
lam_5b_vals = [12*(1 - char_table[k, col_5b]/dims[k]) for k in range(9)]
combined_vals = [A1*la + B1*lb for la, lb in zip(lam_5a_vals, lam_5b_vals)]

print(f"\n  Combined values sorted: "
      f"{[f'{v:.2f}' for v in sorted(combined_vals)]}")
print(f"  Target (sorted):        {sorted(n_bare.values())}")

# Scale to match endpoint
cmin, cmax = min(combined_vals), max(combined_vals)
if cmax > cmin:
    scaled = [(v - cmin)/(cmax - cmin) * 26 for v in combined_vals]
    print(f"\n  Linearly scaled to [0,26]:")
    print(f"  {[f'{v:.1f}' for v in sorted(scaled)]}")
    rms_scaled = np.sqrt(np.mean([(s-t)**2 for s,t in
                  zip(sorted(scaled), sorted(n_bare.values()))]))
    print(f"  RMS vs target: {rms_scaled:.2f}")

# ================================================================
# PART 32: Direct Test - chi_k(10a) as Mass Generator
# ================================================================
print(f"\n{'='*72}")
print("PART 32: chi_k(10a) and chi_k(10b) as Mass Sources")
print("=" * 72)

# On class 10a, chi_2 = phi. These are the "phi-carrying" elements.
# chi_k(10a) gives the character of each irrep on these elements.
print(f"\n  Characters on 10a and 10b classes:")
print(f"  {'Irrep':>10} {'chi(10a)':>10} {'chi(10b)':>10} "
      f"{'5*p+6*q':>10} {'(p,q) of chi(10a)':>20}")

for k in range(9):
    c10a = char_table[k, col_10a]
    c10b = char_table[k, col_10b]
    p, q = to_zphi(c10a)
    F = A1*p + B1*q
    print(f"  rho_{k+1}({int(dims[k]):1d})  {c10a:>+10.4f} {c10b:>+10.4f} "
          f"{F:>+10d}   ({p:+d},{q:+d})")

# ================================================================
# PART 33: The Galois Structure of Characters
# ================================================================
print(f"\n{'='*72}")
print("PART 33: Galois Action on Characters")
print("=" * 72)

# Galois sigma: phi -> phi' = -1/phi
# On characters: sigma(chi(5a)) = chi(5b) for A5 reps
#                sigma(chi(5a)) = -chi(5b) for spin reps? No...
# Actually: Galois acts on the VALUE, not on the class.
# sigma(chi_k(5a)) replaces phi->phi' in the algebraic expression.

print(f"\n  Galois: phi -> phi' = -1/phi")
print(f"  For each irrep, sigma(chi(5a)) should equal chi(5b) if A5 rep")
print(f"  and sigma(chi(10a)) should equal chi(10b)")
print()

for k in range(9):
    c5a = char_table[k, col_5a]
    c5b = char_table[k, col_5b]
    c10a = char_table[k, col_10a]
    c10b = char_table[k, col_10b]

    # Galois conjugate of chi(5a): replace phi->phi' in p+q*phi
    p, q = zphi_5a[k]
    sigma_5a = p + q * PHIp

    # Galois conjugate of chi(10a)
    p10, q10 = to_zphi(c10a)
    sigma_10a = p10 + q10 * PHIp

    print(f"  rho_{k+1}: sigma(chi(5a))={sigma_5a:+.4f} vs chi(5b)={c5b:+.4f} "
          f"{'MATCH' if abs(sigma_5a-c5b)<0.01 else 'NO'} | "
          f"sigma(chi(10a))={sigma_10a:+.4f} vs chi(10b)={c10b:+.4f} "
          f"{'MATCH' if abs(sigma_10a-c10b)<0.01 else 'NO'}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*72}")
print("SUMMARY")
print("=" * 72)

print(f"""
2I Character Table EXTRACTED from affine E8 eigenvectors.

Key findings:
1. Perron eigenvector (lambda=2) -> chi_k(1a) = dim(rho_k) EXACT
2. Lambda=-2 eigenvector -> chi_k(2a) distinguishes A5 vs spin reps
3. Characters decompose as p + q*phi in Z[phi]
4. Galois sigma acts on characters by phi -> phi'

Mass formula test results shown above.
""")

print("Done.")

"""
exp423_ryu_takayanagi.py
========================
Ryu-Takayanagi on the 600-cell: Holographic Newton's Constant.

Goal: Extract an effective Newton's constant G_eff from the Ryu-Takayanagi
formula S_EE = A/(4G) applied to the Cayley graph of 2I, using the TQFT
entanglement entropy from SU(2)_3 Chern-Simons theory (k+2 = a1 = 5).

Background:
  - exp422 showed kappa_Ollivier = 0 (Ricci-flat), so G cannot come from curvature.
  - The RT formula in tensor networks gives G_eff = 1/(4*ln(chi)) where chi is
    the bond dimension. The 600-cell has a natural TQFT (SU(2)_3) with known
    quantum dimensions.
  - If we connect the min-cut of the Cayley graph to the TQFT entanglement
    entropy, we might extract G from the framework.

Key formulas:
  D^2 = a1 + sqrt(a1) = 5 + sqrt(5) = 7.236
  TEE: gamma = ln(D) = 0.9895
  Tensor network RT: S(A) = |min-cut| * ln(chi)
  3D gravity (Witten): k = l/(4G), G_3D = l/(4k) = l/12

Framework constants:
  a_1 = 5, phi = (1+sqrt(5))/2, N = 120, degree = 12

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path
from itertools import permutations, product as cartesian_product
import time
import math

# =====================================================================
# CONSTANTS
# =====================================================================
PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2
SQRT5 = np.sqrt(5)
a_1 = 5
b_1 = 6
N_ORDER = 120        # |2I|
DEGREE = 12          # valency of Cayley graph
N_gen = 3
N_eig = 9           # distinct eigenvalues
TOL = 1e-8
EDGE_TOL = 1e-6

# Framework couplings
alpha_s = 1 / (2 * PHI**3)
sin2_tW = b_1 / (a_1**2 + 1)
# alpha from quadratic: 2*pi*a^2 - (N/b1)*phi^4*a + 1 = 0
_A_coef = 2 * np.pi
_B_coef = -(N_ORDER / b_1) * PHI**4
_C_coef = 1.0
_disc = _B_coef**2 - 4 * _A_coef * _C_coef
alpha_em = (-_B_coef - np.sqrt(_disc)) / (2 * _A_coef)

N_PASS = 0
N_FAIL = 0
N_PATTERN = 0


def check(condition, label, detail=""):
    global N_PASS, N_FAIL
    if condition:
        N_PASS += 1
        tag = "PASS"
    else:
        N_FAIL += 1
        tag = "FAIL"
    print(f"  [{tag}] {label}")
    if detail:
        print(f"         {detail}")
    return condition


def pattern(condition, label, detail=""):
    global N_PATTERN
    if condition:
        N_PATTERN += 1
    tag = "PATTERN" if condition else "NEGATIVE"
    print(f"  [{tag}] {label}")
    if detail:
        print(f"           {detail}")
    return condition


# =====================================================================
# PART 1: Build the 600-cell Cayley graph (reused from exp422)
# =====================================================================
print("=" * 72)
print("EXP423: RYU-TAKAYANAGI ON THE 600-CELL")
print("        Holographic Newton's Constant")
print("=" * 72)

print("\n--- PART 1: Construct 600-cell Cayley graph ---")


def quat_mult(q1, q2):
    """Hamilton product of two quaternions (w, x, y, z)."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])


def build_2I():
    """Construct all 120 unit quaternions of the binary icosahedral group 2I."""
    verts = set()

    def add_vert(v):
        arr = np.array(v, dtype=float)
        n = np.linalg.norm(arr)
        if n > 1e-12:
            arr = arr / n
        verts.add(tuple(np.round(arr, 10)))

    # Type A: 8 axis quaternions
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            add_vert(v)

    # Type B: 16 half-integer quaternions
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    add_vert([s0, s1, s2, s3])

    # Type C: 96 golden-ratio quaternions
    base = [0.0, 0.5, PHI / 2.0, 1.0 / (2.0 * PHI)]
    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4)
                  if p[i] > p[j])
        if inv % 2 == 0:
            even_perms.append(p)

    for perm in even_perms:
        coords = [base[perm[i]] for i in range(4)]
        nz_indices = [i for i in range(4) if abs(coords[i]) > 1e-12]
        for signs in cartesian_product([1, -1], repeat=len(nz_indices)):
            v = list(coords)
            for idx, s in zip(nz_indices, signs):
                v[idx] *= s
            add_vert(v)

    return np.array(sorted(verts))


t0 = time.time()
verts = build_2I()
N = len(verts)
check(N == 120, f"N = |2I| = {N}")

# Build adjacency matrix
dot_threshold = PHI / 2.0
dots = verts @ verts.T
np.clip(dots, -1.0, 1.0, out=dots)

A = np.zeros((N, N), dtype=float)
neighbors = [[] for _ in range(N)]
for i in range(N):
    for j in range(i + 1, N):
        if abs(dots[i, j] - dot_threshold) < EDGE_TOL:
            A[i, j] = 1.0
            A[j, i] = 1.0
            neighbors[i].append(j)
            neighbors[j].append(i)

degrees = np.sum(A, axis=1).astype(int)
num_edges = int(np.sum(A) / 2)
check(np.all(degrees == DEGREE), f"All vertices have degree {DEGREE}")
check(num_edges == 720, f"Number of edges = {num_edges}")

# Shortest-path distance matrix
print("\n  Computing shortest-path distances...")
dist_matrix = shortest_path(A, method='D', unweighted=True).astype(int)
diameter = np.max(dist_matrix)
check(diameter == 5, f"Diameter of Cayley graph = {diameter}")

# Edge list
edge_list = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            edge_list.append((i, j))
check(len(edge_list) == 720, f"Edge list has {len(edge_list)} edges")

t_build = time.time() - t0
print(f"  Build time: {t_build:.2f}s")

# Laplacian
D_mat = np.diag(degrees.astype(float))
L = D_mat - A


# =====================================================================
# PART 2: TQFT Entanglement Entropy (analytical, from exp404/exp405)
# =====================================================================
print("\n--- PART 2: TQFT Entanglement Entropy ---")

k = a_1 - 2   # CS level k = 3

# Quantum dimensions: d_j = sin(pi*(2j+1)/(k+2)) / sin(pi/(k+2))
print(f"  SU(2)_{k} Chern-Simons (k+2 = {k+2} = a1)")
print(f"\n  Quantum dimensions:")

q_dims = []
for j2 in range(k + 1):  # j = 0, 1/2, 1, 3/2
    j = j2 / 2
    d_j = math.sin(math.pi * (j2 + 1) / (k + 2)) / math.sin(math.pi / (k + 2))
    q_dims.append(d_j)
    label = ""
    if abs(d_j - 1) < 1e-10:
        label = " = 1"
    elif abs(d_j - PHI) < 1e-10:
        label = " = phi"
    print(f"    d(j={j}) = {d_j:.6f}{label}")

# Total quantum dimension D^2 = sum d_j^2
D_sq = sum(d**2 for d in q_dims)
D_total = math.sqrt(D_sq)
gamma_TEE = math.log(D_total)

print(f"\n  D^2 = sum(d_j^2) = {D_sq:.6f}")
print(f"  D^2 = a1 + sqrt(a1) = {a_1 + math.sqrt(a_1):.6f}")
check(abs(D_sq - (a_1 + math.sqrt(a_1))) < TOL,
      f"D^2 = a1 + sqrt(a1) = {D_sq:.6f}")

print(f"  D = {D_total:.6f}")
print(f"  TEE: gamma = ln(D) = {gamma_TEE:.6f}")

# TQFT entanglement entropy (Shannon entropy of probabilities p_j = d_j^2/D^2)
p_j = [d**2 / D_sq for d in q_dims]
S_TQFT = -sum(p * math.log(p) for p in p_j if p > 0)
print(f"\n  TQFT Shannon entropy:")
print(f"    p_j = d_j^2 / D^2:")
for j2 in range(k + 1):
    j = j2 / 2
    print(f"      p(j={j}) = {p_j[j2]:.6f}")
print(f"    S_TQFT = -sum(p_j * ln(p_j)) = {S_TQFT:.6f}")

# Galois conjugate sector (from exp405)
D_sq_prime = a_1 - math.sqrt(a_1)
D_prime = math.sqrt(D_sq_prime)
gamma_prime = math.log(D_prime)

print(f"\n  Galois conjugate sector:")
print(f"    D'^2 = a1 - sqrt(a1) = {D_sq_prime:.6f}")
print(f"    D' = {D_prime:.6f}")
print(f"    gamma' = ln(D') = {gamma_prime:.6f}")
print(f"    gamma - gamma' = {gamma_TEE - gamma_prime:.6f}")
print(f"    ln(phi) = {math.log(PHI):.6f}")
check(abs((gamma_TEE - gamma_prime) - math.log(PHI)) < TOL,
      "gamma - gamma' = ln(phi)")

# Heegaard entropy: for genus-1 Heegaard splitting of S^3
# S(Heegaard) = ln(D^2) = 2*gamma
S_Heegaard = math.log(D_sq)
print(f"\n  Heegaard entropy (genus-1 splitting of S^3):")
print(f"    S_Heegaard = ln(D^2) = 2*gamma = {S_Heegaard:.6f}")


# =====================================================================
# PART 3: Natural bipartitions of 120 vertices
# =====================================================================
print("\n--- PART 3: Natural Bipartitions ---")

# For each bipartition: compute |S|, complement |S^c|, cut size |dS|


def cut_size(partition_A, adj_matrix):
    """Count edges between partition A and its complement."""
    mask_A = np.zeros(adj_matrix.shape[0], dtype=bool)
    mask_A[list(partition_A)] = True
    mask_B = ~mask_A
    # Count edges from A to B
    return int(np.sum(adj_matrix[np.ix_(mask_A, mask_B)]))


bipartitions = {}

# 3a) Spectral (Fiedler vector): sign of eigenvector for lambda_1
print("\n  3a) Spectral bipartition (Fiedler vector):")
evals_L, evecs_L = eigh(L)
# lambda_0 = 0 (constant mode), lambda_1 = algebraic connectivity
idx_sorted = np.argsort(evals_L)
lambda_1 = evals_L[idx_sorted[1]]
fiedler = evecs_L[:, idx_sorted[1]]

# Fiedler bipartition: sign of Fiedler vector
A_fiedler = set(np.where(fiedler >= 0)[0])
B_fiedler = set(np.where(fiedler < 0)[0])
cut_fiedler = cut_size(A_fiedler, A)

print(f"    lambda_1 (algebraic connectivity) = {lambda_1:.6f}")
print(f"    Expected: 12 - 6*phi = {12 - 6*PHI:.6f}")
check(abs(lambda_1 - (12 - 6 * PHI)) < 0.01,
      f"lambda_1 = 12 - 6*phi = {lambda_1:.6f}")
print(f"    |A| = {len(A_fiedler)}, |B| = {len(B_fiedler)}")
print(f"    Cut size |dS| = {cut_fiedler}")
print(f"    Ratio |dS|/min(|A|,|B|) = {cut_fiedler / min(len(A_fiedler), len(B_fiedler)):.4f}")

bipartitions['spectral'] = {
    'A': A_fiedler, 'cut': cut_fiedler,
    'size_A': len(A_fiedler), 'size_B': len(B_fiedler)
}

# 3b) Geometric bipartition: cut by x_0=0 (first quaternion coordinate)
print("\n  3b) Geometric bipartition (x_0 = 0 hyperplane):")
A_geo = set(np.where(verts[:, 0] >= 0)[0])
B_geo = set(np.where(verts[:, 0] < 0)[0])
cut_geo = cut_size(A_geo, A)

print(f"    |A| = {len(A_geo)}, |B| = {len(B_geo)}")
print(f"    Cut size |dS| = {cut_geo}")
print(f"    Ratio |dS|/min(|A|,|B|) = {cut_geo / min(len(A_geo), len(B_geo)):.4f}")

bipartitions['geometric_x0'] = {
    'A': A_geo, 'cut': cut_geo,
    'size_A': len(A_geo), 'size_B': len(B_geo)
}

# 3c) Distance-based: vertices at d<=2 from vertex 0 vs rest
print("\n  3c) Distance-based bipartition (d <= 2 from vertex 0):")
A_dist = set(np.where(dist_matrix[0] <= 2)[0])
B_dist = set(range(N)) - A_dist
cut_dist = cut_size(A_dist, A)

print(f"    |A| (d<=2) = {len(A_dist)}, |B| (d>2) = {len(B_dist)}")
print(f"    Cut size |dS| = {cut_dist}")
# Distance shell sizes
for d in range(diameter + 1):
    shell = np.sum(dist_matrix[0] == d)
    print(f"      d={d}: {shell} vertices")

bipartitions['distance'] = {
    'A': A_dist, 'cut': cut_dist,
    'size_A': len(A_dist), 'size_B': len(B_dist)
}

# 3d) Antipodal bipartition: vertex i and its antipode (d=5) in same/different sets
print("\n  3d) Antipodal bipartition (top half vs bottom half in S^3):")
# Group vertices by distance from vertex 0, take d<=2 vs d>=3
# This gives a "hemisphere" cut
A_hemi = set(np.where(dist_matrix[0] <= 2)[0])
B_hemi = set(range(N)) - A_hemi
# Already same as 3c, so use a different cut: x_0 + x_1 >= 0
vals = verts[:, 0] + verts[:, 1]
A_diag = set(np.where(vals >= 0)[0])
B_diag = set(range(N)) - A_diag
cut_diag = cut_size(A_diag, A)

print(f"    (x_0 + x_1 >= 0 hyperplane)")
print(f"    |A| = {len(A_diag)}, |B| = {len(B_diag)}")
print(f"    Cut size |dS| = {cut_diag}")

bipartitions['diagonal'] = {
    'A': A_diag, 'cut': cut_diag,
    'size_A': len(A_diag), 'size_B': len(B_diag)
}

# 3e) Random: average over 100 random balanced (60-60) splits
print("\n  3e) Random balanced bipartitions (60-60):")
rng = np.random.RandomState(42)
random_cuts = []
for trial in range(100):
    perm = rng.permutation(N)
    A_rand = set(perm[:60])
    c = cut_size(A_rand, A)
    random_cuts.append(c)

random_cuts = np.array(random_cuts)
print(f"    Mean cut size = {random_cuts.mean():.1f}")
print(f"    Std  cut size = {random_cuts.std():.1f}")
print(f"    Min  cut size = {random_cuts.min()}")
print(f"    Max  cut size = {random_cuts.max()}")

bipartitions['random_mean'] = {
    'A': None, 'cut': random_cuts.mean(),
    'size_A': 60, 'size_B': 60
}

# Summary table
print(f"\n  Bipartition summary:")
print(f"  {'Name':<20} {'|A|':>5} {'|B|':>5} {'Cut |dS|':>10} {'|dS|/min':>10}")
print(f"  {'-'*20} {'-'*5} {'-'*5} {'-'*10} {'-'*10}")
for name, bp in bipartitions.items():
    mn = min(bp['size_A'], bp['size_B'])
    ratio = bp['cut'] / mn if mn > 0 else 0
    print(f"  {name:<20} {bp['size_A']:5d} {bp['size_B']:5d} {bp['cut']:10.1f} {ratio:10.4f}")


# =====================================================================
# PART 4: Edge connectivity, Cheeger constant
# =====================================================================
print("\n--- PART 4: Edge Connectivity and Cheeger Constant ---")

# Edge connectivity of a vertex-transitive graph = degree (by Whitney)
# For 12-regular vertex-transitive: edge connectivity = 12
print(f"  Edge connectivity (vertex-transitive) = degree = {DEGREE}")

# Cheeger constant: h = min over S (|dS| / min(|S|, |S^c|))
# with |S| <= N/2
# From spectral gap: lambda_1/2 <= h <= sqrt(2*k*lambda_1)
lambda_1_val = lambda_1
h_lower = lambda_1_val / 2
h_upper = math.sqrt(2 * DEGREE * lambda_1_val)

print(f"\n  Cheeger constant bounds (Cheeger inequality):")
print(f"    lambda_1 = {lambda_1_val:.6f}")
print(f"    lambda_1/2 = {h_lower:.6f} <= h <= sqrt(2*k*lambda_1) = {h_upper:.6f}")

# Compute Cheeger constant from our bipartitions
h_values = {}
for name, bp in bipartitions.items():
    if bp['A'] is not None:
        mn = min(bp['size_A'], bp['size_B'])
        if mn > 0:
            h_val = bp['cut'] / mn
            h_values[name] = h_val

best_h_name = min(h_values, key=h_values.get)
best_h = h_values[best_h_name]

print(f"\n  Cheeger estimates from bipartitions:")
for name, h_val in sorted(h_values.items(), key=lambda x: x[1]):
    print(f"    h({name}) = {h_val:.4f}")

print(f"\n  Best (lowest) Cheeger estimate: h >= {best_h:.4f} (from {best_h_name})")
check(h_lower <= best_h + TOL,
      f"Cheeger bound check: lambda_1/2 = {h_lower:.4f} <= h_est = {best_h:.4f}")

# Full Laplacian spectrum (needed later)
print(f"\n  Full Laplacian spectrum:")
evals_all = evals_L[idx_sorted]
unique_evals = np.unique(np.round(evals_all, 6))
print(f"    {len(unique_evals)} distinct eigenvalues:")
for i, ev in enumerate(unique_evals):
    mult = np.sum(np.abs(evals_all - ev) < 0.001)
    print(f"      lambda_{i} = {ev:.6f} (multiplicity {mult})")


# =====================================================================
# PART 5: Ryu-Takayanagi Tests (CORE)
# =====================================================================
print("\n--- PART 5: Ryu-Takayanagi Tests ---")
print("  For each bipartition with cut size A_cut,")
print("  test G_eff from RT formula S = A/(4G)")

# Use the balanced spectral cut as primary
# Also test geometric cut
test_partitions = {
    'spectral': bipartitions['spectral'],
    'geometric_x0': bipartitions['geometric_x0'],
    'distance': bipartitions['distance'],
    'diagonal': bipartitions['diagonal'],
}

print(f"\n  Framework reference values:")
z_Planck = 4 * PHI**2
print(f"    z_Planck = 4*phi^2 = {z_Planck:.6f}")
print(f"    alpha_em = {alpha_em:.6f}")
print(f"    alpha_s = {alpha_s:.6f}")
print(f"    1/alpha_em = {1/alpha_em:.4f}")
print(f"    gamma (TEE) = {gamma_TEE:.6f}")
print(f"    S_TQFT = {S_TQFT:.6f}")
print(f"    S_Heegaard = {S_Heegaard:.6f}")

for part_name, bp in test_partitions.items():
    A_cut = bp['cut']
    print(f"\n  === Bipartition: {part_name} (cut = {A_cut}) ===")

    # 5a) chi = D (total quantum dimension as bond dimension)
    G_5a = A_cut / (4 * S_TQFT)
    print(f"\n    5a) RT with S = S_TQFT:")
    print(f"        G_RT = A_cut / (4*S_TQFT) = {A_cut} / {4*S_TQFT:.4f} = {G_5a:.6f}")
    z_5a = math.log(G_5a) / math.log(alpha_em) if G_5a > 0 else float('nan')
    z_5a_inv = math.log(1/G_5a) / math.log(alpha_em) if G_5a > 0 else float('nan')
    print(f"        ln(G_RT)/ln(alpha) = {z_5a:.4f} (compare z_Planck = {z_Planck:.4f})")
    print(f"        ln(1/G_RT)/ln(alpha) = {z_5a_inv:.4f}")

    # 5b) chi = phi (Fibonacci bond dimension)
    G_5b = A_cut / (4 * gamma_TEE)
    print(f"\n    5b) RT with S = gamma (TEE):")
    print(f"        G_RT = A_cut / (4*gamma) = {A_cut} / {4*gamma_TEE:.4f} = {G_5b:.6f}")
    z_5b = math.log(G_5b) / math.log(alpha_em) if G_5b > 0 else float('nan')
    print(f"        ln(G_RT)/ln(alpha) = {z_5b:.4f}")

    # 5c) Weighted area: Area = A_cut * ln(phi) (each edge contributes ln(phi))
    Area_weighted = A_cut * math.log(PHI)
    G_5c = Area_weighted / (4 * S_TQFT)
    print(f"\n    5c) Weighted area (edge weight = ln(phi)):")
    print(f"        Area = A_cut*ln(phi) = {Area_weighted:.4f}")
    print(f"        G_RT = Area/(4*S_TQFT) = {G_5c:.6f}")
    z_5c = math.log(G_5c) / math.log(alpha_em) if G_5c > 0 else float('nan')
    print(f"        ln(G_RT)/ln(alpha) = {z_5c:.4f}")

    # 5d) Heegaard entropy as bulk entropy
    G_5d = A_cut / (4 * S_Heegaard)
    print(f"\n    5d) RT with S = S_Heegaard:")
    print(f"        G_RT = A_cut / (4*S_Heeg) = {G_5d:.6f}")
    z_5d = math.log(G_5d) / math.log(alpha_em) if G_5d > 0 else float('nan')
    print(f"        ln(G_RT)/ln(alpha) = {z_5d:.4f}")

    # 5e) Tensor network: 1/(4G) = ln(chi), chi = phi
    G_TN_phi = 1 / (4 * math.log(PHI))
    print(f"\n    5e) Tensor network (chi=phi):")
    print(f"        G_TN = 1/(4*ln(phi)) = {G_TN_phi:.6f}")
    print(f"        S_predicted = A_cut * ln(phi) = {A_cut * math.log(PHI):.4f}")
    print(f"        Compare S_TQFT = {S_TQFT:.4f}")

    # 5f) 3D gravity (Witten): G_3D = l/(4k), k=3
    G_3D_unit = 1.0 / (4 * k)        # l = 1
    G_3D_phi = 1.0 / (PHI * 4 * k)   # l = 1/phi (edge length)
    print(f"\n    5f) 3D gravity (Witten, k={k}):")
    print(f"        G_3D(l=1) = 1/(4k) = {G_3D_unit:.6f}")
    print(f"        G_3D(l=1/phi) = 1/(phi*4k) = {G_3D_phi:.6f}")

# Now the specific search: which A_cut * gamma = framework quantity?
print(f"\n\n  --- Cross-bipartition search for framework matches ---")
print(f"  Testing: A_cut * gamma = ? and A_cut / (4*S) -> z_Planck?")
print()
for part_name, bp in test_partitions.items():
    A_cut = bp['cut']
    prod_gamma = A_cut * gamma_TEE
    prod_ln_phi = A_cut * math.log(PHI)
    G_RT = A_cut / (4 * S_TQFT)

    print(f"  {part_name}: A_cut={A_cut}")
    print(f"    A_cut * gamma = {prod_gamma:.4f}")
    print(f"    A_cut * ln(phi) = {prod_ln_phi:.4f}")
    print(f"    A_cut / N = {A_cut / N_ORDER:.4f}")
    print(f"    A_cut / (2*|E|) = {A_cut / (2 * num_edges):.4f}")
    print(f"    A_cut / DEGREE = {A_cut / DEGREE:.4f}")

    # Check various framework quantities
    for target_name, target_val in [
        ("N", N_ORDER), ("2*N", 2*N_ORDER), ("|E|", num_edges),
        ("N*gamma", N_ORDER * gamma_TEE), ("N*ln(phi)", N_ORDER * math.log(PHI)),
        ("DEGREE*N_gen", DEGREE * N_gen), ("a1^2", a_1**2),
        ("N/b1", N_ORDER / b_1), ("2*pi", 2 * math.pi),
        ("a1*(a1+1)", a_1 * (a_1 + 1)), ("b1^2", b_1**2),
    ]:
        ratio = A_cut / target_val
        if 0.9 < ratio < 1.1:
            print(f"    *** A_cut / {target_name} = {ratio:.6f} (NEAR 1) ***")
    print()


# =====================================================================
# PART 6: Graph Entanglement Entropy (free bosonic field)
# =====================================================================
print("\n--- PART 6: Graph Entanglement Entropy (Bosonic) ---")
print("  Using Laplacian pseudoinverse as covariance matrix")

# Laplacian pseudoinverse: C = L^+ (remove zero mode)
# L has one zero eigenvalue (constant mode)
evals_L_all = evals_all.copy()
evecs_L_all = evecs_L[:, idx_sorted].copy()

# Pseudoinverse: replace 0 eigenvalue with 0 in inverse
evals_pinv = np.zeros_like(evals_L_all)
for i in range(len(evals_L_all)):
    if abs(evals_L_all[i]) > TOL:
        evals_pinv[i] = 1.0 / evals_L_all[i]

# C = V * diag(1/lambda) * V^T (pseudoinverse)
C_cov = evecs_L_all @ np.diag(evals_pinv) @ evecs_L_all.T

# Verify: L @ C_cov should be the projector P = I - (1/N)|1><1|
LC = L @ C_cov
P_expected = np.eye(N) - np.ones((N, N)) / N
print(f"  ||L*L^+ - P||_F = {np.linalg.norm(LC - P_expected):.2e}")
check(np.linalg.norm(LC - P_expected) < 1e-8,
      "L * L^+ = projector (pseudoinverse correct)")

# For a bipartition A, the entanglement entropy of a free bosonic field
# on the graph is:
# S_ent = sum over eigenvalues nu_i of C_A:
#   S = sum [ (nu+1/2)*ln(nu+1/2) - (nu-1/2)*ln(nu-1/2) ] for nu >= 1/2
# Actually, for the bosonic formula, use the restricted correlation matrix.
# The entropy comes from the symplectic eigenvalues.
# Simpler approach: use the relation S = Tr[ (C_A + 1/2) ln(C_A + 1/2)
#                                           - (C_A - 1/2) ln(C_A - 1/2) ]
# But C_A needs to be the restricted covariance of the ground state.
#
# For a free scalar on a graph, the vacuum correlation matrix is
# C_ij = <phi_i phi_j> = (L^+)_ij / 2 + (something from zero mode regularization)
#
# Simpler: use the entanglement entropy from the restricted Laplacian submatrix.
# For region A with |A| = n_A:
#   L_A = L restricted to A (n_A x n_A submatrix)
#   L_B = L restricted to B
#   The entanglement entropy is:
#     S = (1/2) * [ ln det(L_A) + ln det(L_B) - ln det(L/zero mode) ]
#   (for massless free field, with proper regularization)

# More tractable: use the numerical eigenvalue method
# The ground state of H = (1/2) sum L_ij phi_i phi_j has
# <phi_i phi_j> ~ (L^{-1/2})_ij and <pi_i pi_j> ~ (L^{1/2})_ij
# Use the matrix sqrt approach

print("\n  Computing graph entanglement entropy for each bipartition...")

# Use the formula based on restricted Laplacian eigenvalues
# Following Peschel (2003): entanglement from correlation matrices
# For free bosons on a lattice with Hamiltonian H = (1/2)(p^2 + phi L phi),
# the correlation matrices are X_ij = <phi_i phi_j> and P_ij = <pi_i pi_j>
# X = (1/2) L^{-1/2}, P = (1/2) L^{1/2} (in regularized form)

# Remove zero mode for regularization
# Use L_reg = L + epsilon * |0><0| with small epsilon
epsilon_reg = 0.01  # small regularization
L_reg = L.copy()
# Add mass term to remove zero mode
L_reg += epsilon_reg * np.eye(N)

evals_reg, evecs_reg = eigh(L_reg)
# X = (1/2) L_reg^{-1/2}
L_sqrt = evecs_reg @ np.diag(np.sqrt(evals_reg)) @ evecs_reg.T
L_sqrt_inv = evecs_reg @ np.diag(1.0 / np.sqrt(evals_reg)) @ evecs_reg.T

X_mat = 0.5 * L_sqrt_inv  # position-position correlator
P_mat = 0.5 * L_sqrt       # momentum-momentum correlator

# For region A, the entanglement entropy is computed from
# the symplectic eigenvalues of the restricted (X_A, P_A) pair
# nu_i = eigenvalues of sqrt(X_A @ P_A)
# S = sum_i [ (nu_i + 1/2)*ln(nu_i + 1/2) - (nu_i - 1/2)*ln(nu_i - 1/2) ]

S_graph_results = {}

for part_name in ['spectral', 'geometric_x0', 'distance', 'diagonal']:
    bp = test_partitions[part_name]
    A_set = sorted(list(bp['A']))
    n_A = len(A_set)

    X_A = X_mat[np.ix_(A_set, A_set)]
    P_A = P_mat[np.ix_(A_set, A_set)]

    # Symplectic eigenvalues
    M = X_A @ P_A
    evals_M = np.linalg.eigvalsh(M)
    # nu_i should be >= 1/4 (uncertainty principle)
    nu = np.sqrt(np.maximum(evals_M, 0))

    # Entanglement entropy
    S_ent = 0.0
    for nui in nu:
        if nui > 0.5 + 1e-12:
            S_ent += (nui + 0.5) * math.log(nui + 0.5) - (nui - 0.5) * math.log(nui - 0.5)
        elif nui > 0.5 - 1e-12:
            # nu ~ 1/2, minimal uncertainty
            S_ent += 0  # no contribution at minimum

    S_graph_results[part_name] = S_ent
    print(f"\n  {part_name}: |A| = {n_A}")
    print(f"    S_graph (bosonic) = {S_ent:.4f}")
    print(f"    S_graph / ln(phi) = {S_ent / math.log(PHI):.4f}")
    print(f"    S_graph / gamma = {S_ent / gamma_TEE:.4f}")
    print(f"    S_graph / S_TQFT = {S_ent / S_TQFT:.4f}")
    print(f"    S_graph / A_cut = {S_ent / bp['cut']:.6f} (entropy per cut edge)")

# Area law check: S ~ c * |cut|
print(f"\n  Area law check: S_graph vs cut size")
for part_name in ['spectral', 'geometric_x0', 'distance', 'diagonal']:
    bp = test_partitions[part_name]
    S = S_graph_results[part_name]
    if bp['cut'] > 0:
        ratio = S / bp['cut']
        print(f"    {part_name}: S/|cut| = {ratio:.6f}")

# Compare with gamma
print(f"\n  Compare S_graph with TQFT quantities:")
print(f"    TEE gamma = {gamma_TEE:.6f}")
print(f"    S_TQFT = {S_TQFT:.6f}")
print(f"    S_Heegaard = {S_Heegaard:.6f}")


# =====================================================================
# PART 7: Framework Connections (CORE ANALYSIS)
# =====================================================================
print("\n--- PART 7: Framework Connections ---")

# Use the spectral (Fiedler) bipartition as the canonical one
canonical = 'spectral'
bp_canon = test_partitions[canonical]
A_cut_canon = bp_canon['cut']
S_graph_canon = S_graph_results[canonical]

print(f"\n  Canonical bipartition: {canonical}")
print(f"  A_cut = {A_cut_canon}, S_graph = {S_graph_canon:.6f}")

# Extract G_RT using various entropy definitions
print(f"\n  7a) G_RT from graph entropy:")
G_from_graph = A_cut_canon / (4 * S_graph_canon)
print(f"    G_RT = A_cut / (4*S_graph) = {G_from_graph:.6f}")
z_graph = math.log(1 / G_from_graph) / math.log(1 / alpha_em)
print(f"    ln(1/G_RT) / ln(1/alpha) = {z_graph:.4f}")
print(f"    Compare z_Planck = 4*phi^2 = {z_Planck:.4f}")
pattern(abs(z_graph - z_Planck) / z_Planck < 0.05,
        f"z_graph = {z_graph:.4f} vs z_Planck = {z_Planck:.4f} "
        f"({abs(z_graph - z_Planck)/z_Planck*100:.1f}%)")

print(f"\n  7b) G_RT from TQFT entropy:")
G_from_TQFT = A_cut_canon / (4 * S_TQFT)
print(f"    G_RT = A_cut / (4*S_TQFT) = {G_from_TQFT:.6f}")
z_TQFT = math.log(1 / G_from_TQFT) / math.log(1 / alpha_em)
print(f"    ln(1/G_RT) / ln(1/alpha) = {z_TQFT:.4f}")
pattern(abs(z_TQFT - z_Planck) / z_Planck < 0.05,
        f"z_TQFT = {z_TQFT:.4f} vs z_Planck = {z_Planck:.4f} "
        f"({abs(z_TQFT - z_Planck)/z_Planck*100:.1f}%)")

print(f"\n  7c) G_RT from TEE gamma:")
G_from_TEE = A_cut_canon / (4 * gamma_TEE)
print(f"    G_RT = A_cut / (4*gamma) = {G_from_TEE:.6f}")
z_TEE = math.log(1 / G_from_TEE) / math.log(1 / alpha_em)
print(f"    ln(1/G_RT) / ln(1/alpha) = {z_TEE:.4f}")
pattern(abs(z_TEE - z_Planck) / z_Planck < 0.05,
        f"z_TEE = {z_TEE:.4f} vs z_Planck = {z_Planck:.4f} "
        f"({abs(z_TEE - z_Planck)/z_Planck*100:.1f}%)")

# 7d) 3D gravity: G_3D = 1/(4k) and S(cut) = A_cut * gamma
# In 3D, RT: S = L/(4G) where L is length of geodesic (= cut on boundary)
print(f"\n  7d) 3D gravity (k = {k}):")
G_3D = 1.0 / (4 * k)
print(f"    G_3D = 1/(4k) = {G_3D:.6f}")
print(f"    1/(4*G_3D) = 4*k = {4*k}")
S_predicted_3D = A_cut_canon * G_3D
print(f"    S = A_cut * G_3D = {A_cut_canon} * {G_3D:.4f} = {S_predicted_3D:.4f}")
print(f"    (This is NOT the right formula -- S should be A/(4G), not A*G)")
S_from_3D = A_cut_canon / (4 * G_3D)
print(f"    S = A_cut / (4*G_3D) = {S_from_3D:.4f}")
print(f"    Compare: actual graph S = {S_graph_canon:.4f}")

# 7e) Key test: A_balanced * gamma = framework quantity?
print(f"\n  7e) Products A_cut * (entropy per edge):")
prod_gamma = A_cut_canon * gamma_TEE
prod_ln_phi = A_cut_canon * math.log(PHI)
prod_ln2 = A_cut_canon * math.log(2)
print(f"    A_cut * gamma = {prod_gamma:.4f}")
print(f"    A_cut * ln(phi) = {prod_ln_phi:.4f}")
print(f"    A_cut * ln(2) = {prod_ln2:.4f}")

# Test against N, |E|, etc.
for name, val in [
    ("N = 120", N_ORDER),
    ("|E| = 720", num_edges),
    ("N*N_gen = 360", N_ORDER * N_gen),
    ("a1^2 = 25", a_1**2),
    ("b1^2 = 36", b_1**2),
    ("N/2 = 60", N_ORDER / 2),
    ("4*pi = 12.57", 4 * math.pi),
    ("8*pi^2 = 78.96", 8 * math.pi**2),
    ("2*pi^2 = 19.74 (Vol S^3)", 2 * math.pi**2),
]:
    if abs(prod_gamma) > TOL:
        ratio = prod_gamma / val
        if 0.8 < ratio < 1.2:
            print(f"    A_cut * gamma / ({name}) = {ratio:.6f}")

# 7f) Direct test: does G_eff = alpha^(z_Planck) for some A_cut?
print(f"\n  7f) What A_cut gives z_Planck = {z_Planck:.4f}?")
# G_eff = A_cut / (4*S) => 1/G_eff = 4*S/A_cut
# We want: 1/G_eff = alpha^(-z_Planck) = (1/alpha)^z_Planck
# => 4*S/A_cut = (1/alpha)^z_Planck
# => A_cut = 4*S * alpha^z_Planck
target_G_inv = (1 / alpha_em) ** z_Planck
print(f"    (1/alpha)^z_Planck = {target_G_inv:.4e}")

for S_name, S_val in [("S_TQFT", S_TQFT), ("gamma", gamma_TEE),
                       ("S_Heeg", S_Heegaard), ("ln(phi)", math.log(PHI))]:
    A_needed = 4 * S_val / target_G_inv
    print(f"    For S = {S_name}: A_cut needed = {A_needed:.4e} (vs actual {A_cut_canon})")
    print(f"      Ratio: actual/needed = {A_cut_canon / A_needed:.4e}")

# 7g) Tensor network interpretation
print(f"\n  7g) Tensor network: bond dimension chi")
print(f"    If each edge carries bond dim chi, entropy per edge = ln(chi)")
print(f"    S_total = |cut| * ln(chi)")
print(f"    From S_graph: ln(chi) = S_graph / A_cut = {S_graph_canon / A_cut_canon:.6f}")
chi_eff = math.exp(S_graph_canon / A_cut_canon)
print(f"    chi_eff = {chi_eff:.6f}")
print(f"    Compare: phi = {PHI:.6f}")
print(f"    Compare: D = {D_total:.6f}")
print(f"    Compare: phi^2 = {PHI**2:.6f}")
print(f"    Compare: 2 = {2.0}")

pattern(abs(chi_eff - PHI) / PHI < 0.05,
        f"chi_eff = {chi_eff:.4f} vs phi = {PHI:.4f} "
        f"({abs(chi_eff - PHI)/PHI*100:.1f}%)",
        "Bond dimension = golden ratio?")

# Also check if chi = phi gives a reasonable G
if abs(chi_eff - PHI) / PHI < 0.2:
    G_TN = 1 / (4 * math.log(PHI))
    z_TN = math.log(1 / G_TN) / math.log(1 / alpha_em)
    print(f"\n    If chi = phi: G_TN = 1/(4*ln(phi)) = {G_TN:.6f}")
    print(f"    z_TN = ln(1/G_TN)/ln(1/alpha) = {z_TN:.4f}")
    print(f"    This is NOT the Planck hierarchy z={z_Planck:.4f}")
    print(f"    G_TN is O(1), not exponentially small.")


# =====================================================================
# PART 8: Summary and Verdict
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: SUMMARY AND VERDICT")
print("=" * 72)

print(f"""
  INFRASTRUCTURE:
    - 600-cell Cayley graph: N=120, |E|=720, degree=12, diameter=5
    - TQFT (SU(2)_3): D^2 = a1+sqrt(a1) = {D_sq:.4f}
    - TEE: gamma = ln(D) = {gamma_TEE:.6f}
    - S_TQFT (Shannon) = {S_TQFT:.6f}
    - S_Heegaard = {S_Heegaard:.6f}

  BIPARTITIONS:""")
for name, bp in bipartitions.items():
    print(f"    {name}: |A|={bp['size_A']}, cut={bp['cut']:.0f}")

print(f"""
  CHEEGER CONSTANT:
    lambda_1 = {lambda_1_val:.6f} = 12 - 6*phi
    Cheeger bounds: {h_lower:.4f} <= h <= {h_upper:.4f}
    Best estimate: h ~ {best_h:.4f} ({best_h_name})

  GRAPH ENTANGLEMENT:""")
for name, S in S_graph_results.items():
    bp = test_partitions[name]
    ratio = S / bp['cut'] if bp['cut'] > 0 else 0
    print(f"    {name}: S = {S:.4f}, S/|cut| = {ratio:.6f}")

print(f"""
  RT EXTRACTIONS (spectral bipartition, A_cut = {A_cut_canon}):
    G_RT(S_TQFT)  = {A_cut_canon/(4*S_TQFT):.4f}  ->  z = {math.log(1/(A_cut_canon/(4*S_TQFT)))/math.log(1/alpha_em):.4f}
    G_RT(gamma)   = {A_cut_canon/(4*gamma_TEE):.4f}  ->  z = {math.log(1/(A_cut_canon/(4*gamma_TEE)))/math.log(1/alpha_em):.4f}
    G_RT(S_graph) = {G_from_graph:.4f}  ->  z = {z_graph:.4f}
    Target: z_Planck = 4*phi^2 = {z_Planck:.4f}

  TENSOR NETWORK:
    chi_eff = exp(S_graph/A_cut) = {chi_eff:.4f}
    Compare: phi = {PHI:.4f}, D = {D_total:.4f}

  3D GRAVITY (Witten):
    G_3D = 1/(4k) = {1/(4*k):.4f} (k = a1-2 = {k})
""")

# Final verdict
print(f"  VERDICT:")
print(f"  ========")
print()

# Check if ANY z is close to z_Planck
z_values = {
    "S_TQFT": math.log(1 / (A_cut_canon / (4*S_TQFT))) / math.log(1/alpha_em),
    "gamma": math.log(1 / (A_cut_canon / (4*gamma_TEE))) / math.log(1/alpha_em),
    "S_graph": z_graph,
}

any_match = any(abs(z - z_Planck) / z_Planck < 0.1 for z in z_values.values())

if any_match:
    best_z_name = min(z_values, key=lambda k: abs(z_values[k] - z_Planck))
    best_z = z_values[best_z_name]
    print(f"  PATTERN: z({best_z_name}) = {best_z:.4f} is within 10% of z_Planck = {z_Planck:.4f}")
    print(f"  The RT formula gives an effective G whose scale is CONSISTENT")
    print(f"  with the hierarchy, but the match is not exact.")
else:
    print(f"  NEGATIVE: No RT extraction gives z close to z_Planck = {z_Planck:.4f}")
    print(f"  The extracted z-values are:")
    for name, z in z_values.items():
        print(f"    z({name}) = {z:.4f} (off by {abs(z-z_Planck)/z_Planck*100:.1f}%)")

print(f"""
  REASON: The RT formula S = A/(4G) gives G of ORDER UNITY (G ~ 10-100).
  Newton's constant requires G ~ alpha^{{z_Planck}} ~ 10^{{-{z_Planck/math.log10(1/alpha_em):.0f}}}.
  The Cayley graph cut sizes (~ 100-400 edges) are TOO SMALL relative
  to the TQFT entropy (~ 1) to produce an exponentially small G.

  The RT formula on a FINITE graph with O(100) vertices CANNOT produce
  the hierarchy. The Planck hierarchy z = 4*phi^2 must come from the
  SPECTRAL ACTION (exp413) or another non-combinatorial mechanism.

  STRUCTURAL INSIGHTS (positive):
  1. The graph entanglement entropy S_graph obeys an AREA LAW:
     S ~ c * |cut|, with c = S/|cut| ~ {S_graph_canon/A_cut_canon:.4f} per edge.
  2. The effective bond dimension chi_eff = {chi_eff:.4f}.
  3. The TQFT entropy gamma = ln(D) = {gamma_TEE:.4f} is a topological invariant.
  4. The Cheeger constant satisfies spectral bounds.
  5. The 3D gravity G_3D = 1/12 from Witten's CS/gravity is a
     TOPOLOGICAL Newton's constant, not the dynamical one.

  NEGATIVE for deriving G from RT on the 600-cell.
  The hierarchy problem is NOT solved by holographic entanglement on
  a finite Cayley graph.
""")

# =====================================================================
# Final count
# =====================================================================
print("=" * 72)
print(f"  Infrastructure: {N_PASS} PASS, {N_FAIL} FAIL, {N_PATTERN} PATTERN")
print("=" * 72)

if N_FAIL == 0:
    print("  STATUS: ALL INFRASTRUCTURE TESTS PASS")
else:
    print(f"  STATUS: {N_FAIL} INFRASTRUCTURE FAILURES")

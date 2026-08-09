"""
exp110_spectral_anatomy.py
==========================
Spectral Anatomy of the 600-Cell: Eigenvalues, Association Scheme,
and the Search for Fermion Mass Patterns

MOTIVATION:
The 600-cell has EXACTLY 9 angular distance classes (shells of distinct
angular separation from any vertex). This means its adjacency matrix
has exactly 9 distinct eigenvalues. We have exactly 9 massive fermions
(e, mu, tau, u, d, s, c, b, t). Coincidence or structure?

APPROACH:
A) Build 600-cell and identify all 9 angular distance classes
B) Compute full adjacency spectrum (120 eigenvalues, 9 distinct)
C) Compute association scheme eigenmatrix (P-matrix)
D) Test: do eigenvalue ratios connect to mass ratios?
E) Compute Hodge spectrum (edge Laplacian)
F) Dirac operator spectrum (vertices + edges)

Category: DERIVAT (exact computation on known geometry)
"""

import numpy as np
from itertools import product, permutations
from collections import defaultdict, Counter

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP110: SPECTRAL ANATOMY OF THE 600-CELL")
print("=" * 70)

# ============================================================
# PART A: BUILD 600-CELL AND ANGULAR DISTANCE CLASSES
# ============================================================
print("\n--- PART A: 600-cell angular distance classes ---")

def build_600cell():
    phi = PHI
    iphi = 1.0 / phi
    verts = set()
    def add(v):
        n = np.sqrt(sum(x**2 for x in v))
        verts.add(tuple(round(x/n, 10) for x in v))
    for i in range(4):
        for s in [+1, -1]:
            v = [0., 0., 0., 0.]; v[i] = float(s); add(v)
    for ss in product([+1, -1], repeat=4):
        add([s * 0.5 for s in ss])
    base = [phi/2, 0.5, iphi/2, 0.0]
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            pv = [base[p[i]] for i in range(4)]
            nz = [i for i in range(4) if abs(pv[i]) > 1e-12]
            for ss in product([+1, -1], repeat=len(nz)):
                v = list(pv)
                for k, idx in enumerate(nz): v[idx] *= ss[k]
                add(v)
    return np.array([list(v) for v in verts])

V = build_600cell()
N = len(V)
print(f"  Vertices: {N}")

# Compute all pairwise dot products
dots = V @ V.T
np.clip(dots, -1, 1, out=dots)

# Find distinct angular distances (from vertex 0)
angles_from_0 = np.arccos(dots[0])
unique_angles = sorted(set(np.round(angles_from_0, 6)))

print(f"\n  Angular distance classes from vertex 0:")
print(f"  {'Class':>6} {'Angle (rad)':>12} {'Angle (deg)':>12} {'cos(angle)':>12} {'Count':>6} {'as pi/n':>10}")

class_info = []
for k, theta in enumerate(unique_angles):
    count = np.sum(np.abs(angles_from_0 - theta) < 0.001)
    cos_val = np.cos(theta)
    # Try to express as pi/n
    pi_frac = ""
    for n in range(1, 20):
        if abs(theta - np.pi/n) < 0.001:
            pi_frac = f"pi/{n}"
        elif abs(theta - 2*np.pi/n) < 0.001:
            pi_frac = f"2pi/{n}"
        elif abs(theta - 3*np.pi/n) < 0.001:
            pi_frac = f"3pi/{n}"
        elif abs(theta - 4*np.pi/n) < 0.001:
            pi_frac = f"4pi/{n}"
    class_info.append((theta, count, cos_val))
    print(f"  {k:6d} {theta:12.6f} {np.degrees(theta):12.4f} {cos_val:12.6f} {count:6d} {pi_frac:>10}")

print(f"\n  Total classes (including d=0): {len(unique_angles)}")
print(f"  Non-trivial classes: {len(unique_angles) - 1}")
print(f"  Sum of counts: {sum(c for _, c, _ in class_info)}")

# Build adjacency from distance class 1 (nearest neighbors)
adj = defaultdict(set)
edge_set = set()
theta_edge = unique_angles[1]  # first non-zero angle
for i in range(N):
    for j in range(i+1, N):
        if abs(np.arccos(dots[i, j]) - theta_edge) < 0.001:
            adj[i].add(j)
            adj[j].add(i)
            edge_set.add((i, j))
edges = list(edge_set)
print(f"\n  Edges (class 1, theta={np.degrees(theta_edge):.1f} deg): {len(edges)}")

# ============================================================
# PART B: FULL ADJACENCY SPECTRUM
# ============================================================
print("\n--- PART B: Full adjacency spectrum ---")

# Build adjacency matrix
A = np.zeros((N, N))
for i, j in edges:
    A[i, j] = 1.0
    A[j, i] = 1.0

# Eigenvalues
eig_A = np.linalg.eigvalsh(A)
eig_A_sorted = np.sort(eig_A)[::-1]

# Find distinct eigenvalues
eig_rounded = np.round(eig_A_sorted, 4)
unique_eigs, counts_eigs = np.unique(eig_rounded, return_counts=True)
unique_eigs = unique_eigs[::-1]
counts_eigs = counts_eigs[::-1]

print(f"\n  Distinct eigenvalues of adjacency matrix: {len(unique_eigs)}")
print(f"  {'#':>4} {'Eigenvalue':>14} {'Multiplicity':>14} {'Ratio to max':>14}")
for i, (ev, mult) in enumerate(zip(unique_eigs, counts_eigs)):
    ratio = ev / unique_eigs[0] if unique_eigs[0] != 0 else 0
    print(f"  {i:4d} {ev:14.6f} {mult:14d} {ratio:14.6f}")

print(f"\n  Sum of multiplicities: {sum(counts_eigs)} (expect {N})")
print(f"  Largest eigenvalue: {unique_eigs[0]:.6f} (expect 12 = degree)")

# Laplacian eigenvalues
L_eigs = 12.0 - eig_A_sorted  # L = D - A, degree = 12
L_unique = np.round(12.0 - unique_eigs, 4)
print(f"\n  Corresponding Laplacian eigenvalues (lambda = 12 - theta):")
for i, (lev, mult) in enumerate(zip(L_unique, counts_eigs)):
    print(f"    lambda_{i} = {lev:10.4f}  (mult = {mult})")

# ============================================================
# PART C: ASSOCIATION SCHEME STRUCTURE
# ============================================================
print("\n--- PART C: Association scheme ---")

# Build distance matrices D_k for each angular class
n_classes = len(unique_angles)
D = np.zeros((n_classes, N, N))
class_angles = [theta for theta, _, _ in class_info]
class_counts = [count for _, count, _ in class_info]

for i in range(N):
    for j in range(N):
        theta_ij = np.arccos(dots[i, j])
        for k, theta_k in enumerate(class_angles):
            if abs(theta_ij - theta_k) < 0.001:
                D[k, i, j] = 1.0
                break

# Verify: sum of all D_k = J (all-ones matrix)
J = np.sum(D, axis=0)
print(f"  Sum of all D_k = all-ones? max deviation: {np.max(np.abs(J - 1.0)):.10f}")

# Verify: each D_k is symmetric with correct row sums
for k in range(n_classes):
    row_sum = np.sum(D[k], axis=1)
    print(f"    D_{k}: row_sum = {row_sum[0]:.0f} (expect {class_counts[k]}), "
          f"symmetric: {np.allclose(D[k], D[k].T)}")

# Compute intersection numbers: p^k_{ij} = (D_i @ D_j)[0, :] projected onto D_k
# For the association scheme, D_i @ D_j = sum_k p^k_{ij} D_k
print(f"\n  Intersection numbers (D_1 @ D_k for adjacency class):")
D1 = D[1]  # adjacency matrix
for k in range(n_classes):
    product_row = D1[0] @ D[k]  # row 0 of D_1 @ D_k
    # Decompose into classes
    decomp = []
    for m in range(n_classes):
        coeff = np.dot(product_row, D[m][0]) / max(class_counts[m], 1)
        if abs(coeff) > 0.01:
            decomp.append(f"p^{m}_{{1,{k}}}={coeff:.1f}")
    print(f"    D_1 @ D_{k}: {', '.join(decomp)}")

# ============================================================
# PART D: EIGENVALUES AND PHYSICAL CONSTANTS
# ============================================================
print("\n--- PART D: Eigenvalue ratios and physical constants ---")

phi = PHI
alpha_em = 1.0 / 137.036
alpha_s = 1.0 / (2 * phi**3)
sin2_thetaW = 6.0 / 26.0

# Known: lambda_4 / lambda_1 = 2*phi^2 (from exp072)
# lambda_1 = smallest non-zero Laplacian eigenvalue
# Let's check all ratios

print(f"\n  Reference constants:")
print(f"    phi = {phi:.6f}")
print(f"    phi^2 = {phi**2:.6f}")
print(f"    2*phi^2 = {2*phi**2:.6f}")
print(f"    alpha_EM = {alpha_em:.6f}")
print(f"    alpha_s = {alpha_s:.6f}")
print(f"    sin^2(theta_W) = {sin2_thetaW:.6f}")

print(f"\n  Adjacency eigenvalue ratios (theta_i / theta_j):")
for i in range(len(unique_eigs)):
    for j in range(i+1, len(unique_eigs)):
        if abs(unique_eigs[j]) > 0.1:
            ratio = unique_eigs[i] / unique_eigs[j]
            # Check if ratio matches any known constant
            matches = []
            for name, val in [("phi", phi), ("phi^2", phi**2), ("2phi^2", 2*phi**2),
                              ("phi^3", phi**3), ("phi^4", phi**4),
                              ("1/phi", 1/phi), ("2", 2), ("3", 3), ("5", 5),
                              ("phi+1", phi+1), ("phi+2", phi+2),
                              ("2*phi", 2*phi), ("3*phi", 3*phi)]:
                if abs(ratio - val) < 0.01 or (val != 0 and abs(ratio/val - 1) < 0.01):
                    matches.append(f"{name}={val:.4f}")
                if abs(ratio + val) < 0.01 or (val != 0 and abs(ratio/val + 1) < 0.01):
                    matches.append(f"-{name}={-val:.4f}")
            if matches:
                print(f"    theta_{i}/theta_{j} = {ratio:10.6f} = {', '.join(matches)}")

# Check Laplacian ratios
print(f"\n  Laplacian eigenvalue ratios:")
nonzero_L = [(lev, mult) for lev, mult in zip(L_unique, counts_eigs) if abs(lev) > 0.01]
for i in range(len(nonzero_L)):
    for j in range(i+1, len(nonzero_L)):
        if abs(nonzero_L[j][0]) > 0.1:
            ratio = nonzero_L[i][0] / nonzero_L[j][0]
            matches = []
            for name, val in [("phi", phi), ("phi^2", phi**2), ("2phi^2", 2*phi**2),
                              ("1/phi", 1/phi), ("2", 2), ("3", 3), ("5", 5),
                              ("phi+1", phi+1), ("2*phi", 2*phi)]:
                if abs(ratio - val) < 0.02 or (val != 0 and abs(ratio/val - 1) < 0.01):
                    matches.append(f"{name}={val:.4f}")
            if matches:
                li, mi = nonzero_L[i]
                lj, mj = nonzero_L[j]
                print(f"    L({li:.2f},m={mi})/L({lj:.2f},m={mj}) = "
                      f"{ratio:10.6f} = {', '.join(matches)}")

# ============================================================
# PART E: 9 EIGENVALUES vs 9 FERMION MASSES
# ============================================================
print("\n--- PART E: 9 eigenvalues vs 9 fermion masses ---")

# Fermion masses (PDG 2024, running to pole masses where applicable)
fermion_masses = {
    'e':     0.000511,     # GeV
    'mu':    0.10566,      # GeV
    'tau':   1.777,        # GeV
    'u':     0.00216,      # GeV (MS-bar at 2 GeV)
    'd':     0.00467,      # GeV
    's':     0.0934,       # GeV
    'c':     1.27,         # GeV (MS-bar)
    'b':     4.18,         # GeV (MS-bar)
    't':     172.76,       # GeV (pole)
}

# Sort by mass
sorted_fermions = sorted(fermion_masses.items(), key=lambda x: x[1])
sorted_names = [name for name, _ in sorted_fermions]
sorted_masses = np.array([mass for _, mass in sorted_fermions])
log_masses = np.log(sorted_masses)

print(f"  9 fermion masses (sorted):")
for name, mass in sorted_fermions:
    print(f"    {name:>4s}: {mass:.6f} GeV  (log = {np.log(mass):.4f})")

# Sort eigenvalues (9 distinct)
eig_sorted = np.sort(unique_eigs)[::-1]  # descending
n_eig = len(eig_sorted)
print(f"\n  {n_eig} distinct adjacency eigenvalues:")
for i, ev in enumerate(eig_sorted):
    print(f"    theta_{i}: {ev:12.6f}  (mult = {counts_eigs[i]})")

# Test: can we map eigenvalues to log(masses) via linear transformation?
# log(m) = a * theta + b?
if n_eig == 9:
    print(f"\n  TESTING: log(m_i) = a * theta_i + b")
    from numpy.polynomial import polynomial as P
    coeffs = np.polyfit(eig_sorted, log_masses[::-1], 1)  # reverse mass order
    a_fit, b_fit = coeffs
    print(f"    Best fit: a = {a_fit:.6f}, b = {b_fit:.6f}")
    print(f"    {'Fermion':>8} {'Mass (GeV)':>12} {'log(m)':>10} {'theta':>10} "
          f"{'pred log(m)':>12} {'error %':>10}")
    for i, (name, mass) in enumerate(sorted_fermions[::-1]):
        pred = a_fit * eig_sorted[i] + b_fit
        actual = np.log(mass)
        err = (pred - actual) / abs(actual) * 100
        print(f"    {name:>8s} {mass:12.6f} {actual:10.4f} {eig_sorted[i]:10.4f} "
              f"{pred:12.4f} {err:10.1f}%")

    # Also test: log(m) = a * theta^2 + b
    print(f"\n  TESTING: log(m_i) = a * theta_i^2 + b")
    coeffs2 = np.polyfit(eig_sorted**2, log_masses[::-1], 1)
    a2, b2 = coeffs2
    print(f"    Best fit: a = {a2:.6f}, b = {b2:.6f}")
    rms2 = np.sqrt(np.mean((a2 * eig_sorted**2 + b2 - log_masses[::-1])**2))
    rms1 = np.sqrt(np.mean((a_fit * eig_sorted + b_fit - log_masses[::-1])**2))
    print(f"    RMS error (linear): {rms1:.4f}")
    print(f"    RMS error (quadratic): {rms2:.4f}")

elif n_eig > 9:
    print(f"\n  WARNING: {n_eig} eigenvalues > 9 fermions. Cannot do 1-to-1 mapping.")
    print(f"  Trying subset matching...")
elif n_eig < 9:
    print(f"\n  WARNING: {n_eig} eigenvalues < 9 fermions. Not enough for 1-to-1 mapping.")

# Test: eigenvalue ratios vs mass ratios
print(f"\n  EIGENVALUE RATIOS vs MASS RATIOS:")
# Ratio of consecutive eigenvalues
print(f"  Consecutive eigenvalue ratios:")
for i in range(len(eig_sorted)-1):
    if abs(eig_sorted[i+1]) > 0.01:
        ratio = eig_sorted[i] / eig_sorted[i+1]
        print(f"    theta_{i}/theta_{i+1} = {ratio:.6f}")
print(f"  Consecutive mass ratios (descending):")
for i in range(len(sorted_masses)-2, -1, -1):
    ratio = sorted_masses[i+1] / sorted_masses[i] if sorted_masses[i] > 0 else 0
    print(f"    m_{sorted_names[i+1]}/m_{sorted_names[i]} = {ratio:.6f}")

# ============================================================
# PART F: MULTIPLICITIES AND PARTICLE COUNTING
# ============================================================
print("\n--- PART F: Eigenvalue multiplicities ---")

print(f"\n  Multiplicity structure:")
mult_list = sorted(counts_eigs)
print(f"  Multiplicities (sorted): {list(mult_list)}")
print(f"  Sum: {sum(mult_list)} (= {N})")

# Check if multiplicities match any particle counting
print(f"\n  Known particle counts:")
print(f"    Quarks: 6 flavors x 3 colors x 2 chiralities = 36")
print(f"    Leptons: 3 charged x 2 chiralities = 6")
print(f"    Gauge bosons: 8 gluons + 3 weak + 1 photon = 12")
print(f"    Total SM particles (excl. Higgs, neutrinos): 54")
print(f"    Neutrinos (3): add 6 (L+R) = 60")
print(f"    With Higgs: 64")

# Check sums of multiplicities
print(f"\n  Partial sums of multiplicities:")
cumsum = np.cumsum(mult_list)
for i, cs in enumerate(cumsum):
    print(f"    Sum of {i+1} smallest multiplicities: {cs}")

# Specific multiplicity patterns
print(f"\n  Multiplicity decomposition:")
for ev, mult in zip(unique_eigs, counts_eigs):
    # Check if multiplicity is a sum of known numbers
    factors = []
    for a in range(1, mult+1):
        if mult % a == 0:
            factors.append(a)
    print(f"    theta = {ev:10.4f}, mult = {mult:4d}, divisors: {factors}")

# ============================================================
# PART G: EDGE LAPLACIAN (HODGE-1)
# ============================================================
print("\n--- PART G: Edge Laplacian (Hodge-1) spectrum ---")

# Build signed incidence matrix B (N_edges x N_vertices)
# For each edge (i,j) with i<j: B[e, i] = -1, B[e, j] = +1
n_edges = len(edges)
B = np.zeros((n_edges, N))
for e_idx, (i, j) in enumerate(edges):
    B[e_idx, i] = -1.0
    B[e_idx, j] = +1.0

# Vertex Laplacian: L_0 = B^T @ B
L0 = B.T @ B
L0_eigs = np.sort(np.linalg.eigvalsh(L0))
print(f"  Vertex Laplacian L_0 = B^T B: {N}x{N}")
print(f"  L_0 eigenvalues: min={L0_eigs[0]:.6f}, max={L0_eigs[-1]:.6f}")
print(f"  (Should match 12 - adj eigenvalues)")

# Edge Laplacian: L_1 = B @ B^T
L1 = B @ B.T
print(f"\n  Edge Laplacian L_1 = B B^T: {n_edges}x{n_edges}")
print(f"  Computing eigenvalues...")
L1_eigs = np.sort(np.linalg.eigvalsh(L1))
print(f"  L_1 eigenvalues: min={L1_eigs[0]:.6f}, max={L1_eigs[-1]:.6f}")

# Count distinct L1 eigenvalues
L1_unique = np.unique(np.round(L1_eigs, 3))
print(f"  Distinct L_1 eigenvalues: {len(L1_unique)}")

# Zero modes of L_1 (harmonic 1-forms = Betti number b_1)
n_zero_L1 = np.sum(np.abs(L1_eigs) < 0.01)
print(f"  Zero modes of L_1: {n_zero_L1} (expect b_1(S^3) = 0)")

# Non-zero L1 spectrum
L1_nonzero = L1_eigs[np.abs(L1_eigs) > 0.01]
L1_nz_unique, L1_nz_counts = np.unique(np.round(L1_nonzero, 3), return_counts=True)
print(f"\n  Non-zero L_1 eigenvalues (distinct): {len(L1_nz_unique)}")
print(f"  {'#':>4} {'Eigenvalue':>12} {'Multiplicity':>14}")
for i, (ev, mult) in enumerate(zip(L1_nz_unique, L1_nz_counts)):
    print(f"  {i:4d} {ev:12.4f} {mult:14d}")

# ============================================================
# PART H: SIMPLE DIRAC OPERATOR
# ============================================================
print("\n--- PART H: Dirac operator (vertices + edges) ---")

# Dirac operator D = [[0, B^T], [B, 0]]
# Size: (N + n_edges) x (N + n_edges) = 840 x 840
n_total = N + n_edges
D_dirac = np.zeros((n_total, n_total))
D_dirac[:N, N:] = B.T
D_dirac[N:, :N] = B

print(f"  Dirac operator size: {n_total}x{n_total}")
print(f"  Computing eigenvalues...")
D_eigs = np.sort(np.linalg.eigvalsh(D_dirac))
print(f"  Dirac eigenvalues: min={D_eigs[0]:.6f}, max={D_eigs[-1]:.6f}")

# Zero modes
n_zero_D = np.sum(np.abs(D_eigs) < 0.01)
print(f"  Zero modes: {n_zero_D}")
# Expected: dim(ker D) = b_0 + b_1 = 1 + 0 = 1 for vertices,
# but also need to account for edge contribution
# Actually: ker(D) = ker(B) intersect ker(B^T)
# Positive eigenvalues of D are sqrt(eigenvalues of L_0) and sqrt(eigenvalues of L_1)

# Positive Dirac eigenvalues
D_pos = D_eigs[D_eigs > 0.01]
D_pos_unique, D_pos_counts = np.unique(np.round(D_pos, 4), return_counts=True)
print(f"\n  Positive Dirac eigenvalues (distinct): {len(D_pos_unique)}")
print(f"  {'#':>4} {'Eigenvalue':>12} {'lambda^2':>12} {'Multiplicity':>14}")
for i, (ev, mult) in enumerate(zip(D_pos_unique, D_pos_counts)):
    print(f"  {i:4d} {ev:12.6f} {ev**2:12.4f} {mult:14d}")

# Check: D eigenvalues should be +/- sqrt(L_0 eigenvalues) union +/- sqrt(L_1 eigenvalues)
# minus the shared non-zero parts (which are the same for L_0 and L_1)
print(f"\n  Verification: D^2 eigenvalues vs L_0 and L_1")
D_sq_eigs = np.sort(D_eigs**2)
D_sq_unique = np.unique(np.round(D_sq_eigs, 2))
L0_unique_vals = np.unique(np.round(L0_eigs, 2))
L1_unique_vals = np.unique(np.round(L1_eigs, 2))
print(f"  D^2 distinct: {len(D_sq_unique)}")
print(f"  L_0 distinct: {len(L0_unique_vals)}")
print(f"  L_1 distinct: {len(L1_unique_vals)}")
union = np.unique(np.round(np.concatenate([L0_eigs, L1_eigs]), 2))
print(f"  L_0 union L_1 distinct: {len(union)}")

# ============================================================
# PART I: PHI-BASED PATTERNS IN SPECTRUM
# ============================================================
print("\n--- PART I: Golden ratio patterns in spectrum ---")

print(f"\n  Testing: are adjacency eigenvalues expressible as phi-polynomials?")
for ev, mult in zip(unique_eigs, counts_eigs):
    # Try ev = a + b*phi for small integers a, b
    best_match = None
    best_err = 999
    for a in range(-15, 16):
        for b in range(-10, 11):
            val = a + b * phi
            err = abs(ev - val)
            if err < best_err:
                best_err = err
                best_match = (a, b, val)
    a, b, val = best_match
    if best_err < 0.001:
        sign_b = "+" if b >= 0 else ""
        print(f"    theta = {ev:10.4f} = {a}{sign_b}{b}*phi = {val:.6f} "
              f"(mult={mult}, err={best_err:.6f}) EXACT")
    else:
        # Try a + b*phi + c*phi^2
        for a in range(-10, 11):
            for b in range(-10, 11):
                for c in range(-5, 6):
                    val = a + b*phi + c*phi**2
                    err = abs(ev - val)
                    if err < 0.001:
                        print(f"    theta = {ev:10.4f} = {a}+{b}*phi+{c}*phi^2 "
                              f"(mult={mult}) EXACT")
                        break
                else:
                    continue
                break
            else:
                continue
            break
        else:
            print(f"    theta = {ev:10.4f} (mult={mult}): no simple phi-expression found")

# ============================================================
# PART J: MASS FORMULA CONNECTION
# ============================================================
print("\n--- PART J: Connecting spectrum to mass formula ---")

# Our mass formula: m = m_e * phi^n, n = 5a + 6b
# The mass levels are n = 0, 3, 5, 11, 16, 17, 19, 26
# Can these n values be expressed using the eigenvalues?

mass_levels = {'e': 0, 'u': 3, 'd': 5, 'mu': 11, 's': 11, 'c': 16,
               'tau': 17, 'b': 19, 't': 26}

print(f"\n  Mass levels n = 5a + 6b:")
for name, n in sorted(mass_levels.items(), key=lambda x: x[1]):
    print(f"    {name:>4s}: n = {n}")

print(f"\n  Distinct n values: {sorted(set(mass_levels.values()))}")
print(f"  Number of distinct n: {len(set(mass_levels.values()))}")
print(f"  (mu and strange share n=11)")

# Check: 5 and 6 are fundamental parameters
# 5 = diameter of 600-cell graph
# 6 = decagons per vertex
# Can 5 and 6 be extracted from the eigenvalue spectrum?

# The Laplacian eigenvalues in order
print(f"\n  Laplacian spectrum (distinct, ascending):")
L_asc = sorted(L_unique)
for i, lev in enumerate(L_asc):
    mult_i = counts_eigs[list(L_unique).index(lev)] if lev in L_unique else "?"
    print(f"    lambda_{i} = {lev:.6f} (mult = {mult_i})")

# Ratio of specific eigenvalues
if len(L_asc) >= 5:
    print(f"\n  Key Laplacian ratios:")
    # lambda_4/lambda_1 should be 2*phi^2 (from exp072)
    for i in range(1, len(L_asc)):
        for j in range(1, i):
            if L_asc[j] > 0.01:
                ratio = L_asc[i] / L_asc[j]
                for name, val in [("2*phi^2", 2*phi**2), ("phi^2", phi**2),
                                  ("phi", phi), ("2", 2), ("3", 3), ("5", 5),
                                  ("6", 6), ("phi+1", phi+1), ("2*phi", 2*phi),
                                  ("phi^3", phi**3), ("5/3", 5/3), ("6/5", 6/5)]:
                    if abs(ratio/val - 1) < 0.005:
                        print(f"    L_{i}/L_{j} = {ratio:.6f} = {name} ({val:.6f}) "
                              f"err={abs(ratio-val)/val*100:.2f}%")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
600-CELL ASSOCIATION SCHEME:
  Angular distance classes: {len(unique_angles)} (including d=0)
  Non-trivial classes: {len(unique_angles) - 1}
  Class sizes: {[c for _, c, _ in class_info]}

ADJACENCY SPECTRUM:
  Distinct eigenvalues: {len(unique_eigs)}
  Multiplicities: {list(counts_eigs)}
  Largest: {unique_eigs[0]:.4f} (degree)
  Smallest: {unique_eigs[-1]:.4f}

HODGE SPECTRUM:
  L_0 (vertex): {len(np.unique(np.round(L0_eigs, 3)))} distinct eigenvalues
  L_1 (edge): {len(L1_nz_unique)} distinct non-zero eigenvalues
  Dirac zero modes: {n_zero_D}

COMPARISON WITH FERMIONS:
  Number of eigenvalue classes: {len(unique_eigs)}
  Number of fermion masses: 9 (or 8 distinct, mu/s share n=11)
""")

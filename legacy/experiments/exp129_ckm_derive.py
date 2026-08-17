"""
EXP-129: Complete CKM Derivation from 600-Cell/E8 Geometry
===========================================================

Goal: Upgrade CKM mixing angles from PATTERN to DERIVAT (or at least PARTIAL DERIVAT).

Known from exp122/exp125:
  - Cabibbo angle: sin(theta_12) = phi^(-3), 1.8% error
  - Half-integer CKM: n=(3, 13/2, 23/2), mean 6% error
  - Integer CKM: n=(3, 7, 12), mean ~15% error
  - E8 Coxeter: {1,7,11,13,17,19,23,29}, H4 subset {1,11,19,29}
  - H4' complement {7,13,17,23} - ALL prime
  - Divisor /4 and /2 NOT derived

This experiment explores 9 approaches systematically:
  1. Representation overlap in (a,b) space
  2. Level difference pattern |n_up - n_dn|
  3. H4 Weyl group action
  4. E8 branching rules (spectral)
  5. Intersection numbers approach
  6. Decagon winding numbers
  7. Golden angle subdivisions
  8. Direct 600-cell adjacency eigenvector overlap
  9. Wolfenstein parameterization deep analysis

Author: Claude (exp129)
Date: 2026-02-09
"""

import numpy as np
from itertools import product as iterproduct, permutations
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122

print("=" * 80)
print("EXP-129: COMPLETE CKM DERIVATION FROM 600-CELL / E8 GEOMETRY")
print("=" * 80)

# ================================================================
# SETUP: PDG CKM values and quark quantum numbers
# ================================================================
print("\n" + "=" * 80)
print("SETUP: Experimental data and quantum numbers")
print("=" * 80)

# PDG 2024 CKM magnitudes
V_exp = np.array([
    [0.97373, 0.2243, 0.00382],
    [0.221,   0.975,  0.0408 ],
    [0.0086,  0.0415, 1.014  ]
])

V_exp_err = np.array([
    [0.00031, 0.0008, 0.00020],
    [0.004,   0.006,  0.0014 ],
    [0.0002,  0.0011, 0.029  ]
])

labels_up = ['u', 'c', 't']
labels_dn = ['d', 's', 'b']

# Quark quantum numbers (a,b) from exp099/100/112
a_up = np.array([3, 2, 4])    # u, c, t
b_up = np.array([-2, 1, 1])
a_dn = np.array([1, 1, -1])   # d, s, b
b_dn = np.array([0, 1, 4])

# Spectral levels n = 5a + 6b
n_up = 5 * a_up + 6 * b_up    # [3, 16, 26]
n_dn = 5 * a_dn + 6 * b_dn    # [5, 11, 19]

print(f"\nUp-type:   u(a={a_up[0]:+d},b={b_up[0]:+d},n={n_up[0]}), "
      f"c(a={a_up[1]:+d},b={b_up[1]:+d},n={n_up[1]}), "
      f"t(a={a_up[2]:+d},b={b_up[2]:+d},n={n_up[2]})")
print(f"Down-type: d(a={a_dn[0]:+d},b={b_dn[0]:+d},n={n_dn[0]}), "
      f"s(a={a_dn[1]:+d},b={b_dn[1]:+d},n={n_dn[1]}), "
      f"b(a={a_dn[2]:+d},b={b_dn[2]:+d},n={n_dn[2]})")

# Effective phi-exponents for experimental CKM angles
n_exp_12 = -np.log(V_exp[0, 1]) / np.log(PHI)
n_exp_23 = -np.log(V_exp[1, 2]) / np.log(PHI)
n_exp_13 = -np.log(V_exp[0, 2]) / np.log(PHI)

print(f"\nExperimental phi-exponents:")
print(f"  sin(theta_12) = {V_exp[0,1]:.5f} = phi^(-{n_exp_12:.4f})")
print(f"  sin(theta_23) = {V_exp[1,2]:.5f} = phi^(-{n_exp_23:.4f})")
print(f"  sin(theta_13) = {V_exp[0,2]:.6f} = phi^(-{n_exp_13:.4f})")

# Level difference matrices
dn_matrix = np.zeros((3, 3), dtype=int)
for i in range(3):
    for j in range(3):
        dn_matrix[i, j] = abs(n_up[i] - n_dn[j])

da_matrix = np.zeros((3, 3), dtype=int)
db_matrix = np.zeros((3, 3), dtype=int)
for i in range(3):
    for j in range(3):
        da_matrix[i, j] = a_up[i] - a_dn[j]
        db_matrix[i, j] = b_up[i] - b_dn[j]

print(f"\n|n_up - n_dn| matrix:")
print(f"        d({n_dn[0]})   s({n_dn[1]})   b({n_dn[2]})")
for i in range(3):
    row = f"  {labels_up[i]}({n_up[i]:2d})"
    for j in range(3):
        row += f"  {dn_matrix[i,j]:5d}"
    print(row)

# Utility functions
def build_ckm(s12, s23, s13, delta=0):
    """Build standard CKM parametrization."""
    c12 = np.sqrt(max(0, 1 - s12**2))
    c23 = np.sqrt(max(0, 1 - s23**2))
    c13 = np.sqrt(max(0, 1 - s13**2))
    V = np.array([
        [c12*c13, s12*c13, s13],
        [-s12*c23 - c12*s23*s13*np.exp(1j*delta),
         c12*c23 - s12*s23*s13*np.exp(1j*delta),
         s23*c13],
        [s12*s23 - c12*c23*s13*np.exp(1j*delta),
         -c12*s23 - s12*c23*s13*np.exp(1j*delta),
         c23*c13]
    ])
    return np.abs(V)

def compute_errors(V_pred, label="", show=True):
    """Compute and optionally display errors against experimental CKM."""
    abs_err = np.abs(V_pred - V_exp)
    rel_err = abs_err / V_exp * 100
    sigma_err = abs_err / V_exp_err

    mean_rel = np.mean(rel_err)
    max_rel = np.max(rel_err)
    chi2 = np.sum(sigma_err**2)

    # Unitarity check
    prod = V_pred @ V_pred.T
    unitarity_err = np.max(np.abs(prod - np.eye(3)))

    if show:
        print(f"\n--- {label} ---")
        print(f"{'Element':<10} {'Pred':>10} {'Exp':>10} {'Err%':>8} {'Sigma':>8}")
        print("-" * 50)
        for i in range(3):
            for j in range(3):
                name = f"V_{labels_up[i]}{labels_dn[j]}"
                print(f"  {name:<8} {V_pred[i,j]:10.5f} {V_exp[i,j]:10.5f} "
                      f"{rel_err[i,j]:7.2f}% {sigma_err[i,j]:7.1f}s")
        print(f"\n  Mean relative error: {mean_rel:.2f}%")
        print(f"  Max  relative error: {max_rel:.2f}%")
        print(f"  Chi^2 (9 dof):       {chi2:.1f}")
        print(f"  Unitarity max error: {unitarity_err:.2e}")

    return mean_rel, max_rel, chi2


# ================================================================
# BUILD 600-CELL (for approaches needing adjacency matrix)
# ================================================================
print("\n" + "=" * 80)
print("BUILDING 600-CELL")
print("=" * 80)

def build_600cell():
    phi = PHI
    iphi = 1.0 / phi
    verts = set()
    def add(v):
        n = np.sqrt(sum(x**2 for x in v))
        if n > 1e-10:
            verts.add(tuple(round(x/n, 10) for x in v))
    # 8 vertices: permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [+1, -1]:
            v = [0., 0., 0., 0.]; v[i] = float(s); add(v)
    # 16 vertices: (+-1/2)^4
    for ss in iterproduct([+1, -1], repeat=4):
        add([s * 0.5 for s in ss])
    # 96 vertices: even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base = [phi/2, 0.5, iphi/2, 0.0]
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            pv = [base[p[i]] for i in range(4)]
            nz = [i for i in range(4) if abs(pv[i]) > 1e-12]
            for ss in iterproduct([+1, -1], repeat=len(nz)):
                v = list(pv)
                for k, idx in enumerate(nz): v[idx] *= ss[k]
                add(v)
    return np.array([list(v) for v in verts])

V600 = build_600cell()
N600 = len(V600)
print(f"  600-cell vertices: {N600}")

# Adjacency matrix (nearest neighbors = distance class 1)
dots = V600 @ V600.T
np.clip(dots, -1, 1, out=dots)

# Distance class 1: cos(theta) = phi/2 = cos(pi/5) ~ 0.809
cos_nn = PHI / 2  # = cos(36 deg) = cos(pi/5)
A = (np.abs(dots - cos_nn) < 0.001).astype(float)
np.fill_diagonal(A, 0)
degrees = A.sum(axis=1)
print(f"  Adjacency: degree = {int(degrees[0])} (all same: {np.all(degrees == degrees[0])})")

# Full eigendecomposition
eigenvalues_full, eigenvectors_full = np.linalg.eigh(A)
# Sort descending
idx_sort = np.argsort(-eigenvalues_full)
eigenvalues_full = eigenvalues_full[idx_sort]
eigenvectors_full = eigenvectors_full[:, idx_sort]

# Identify distinct eigenvalues
unique_evals = []
mults = []
tol = 0.01
i = 0
while i < len(eigenvalues_full):
    val = eigenvalues_full[i]
    count = 1
    while i + count < len(eigenvalues_full) and abs(eigenvalues_full[i + count] - val) < tol:
        count += 1
    unique_evals.append(val)
    mults.append(count)
    i += count

print(f"\n  9 distinct eigenvalues with multiplicities:")
for k, (ev, m) in enumerate(zip(unique_evals, mults)):
    print(f"    lambda_{k} = {ev:10.4f}, mult = {m:3d} = {int(np.sqrt(m))}^2")


# ================================================================
# APPROACH 1: Representation overlap in (a,b) space
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 1: REPRESENTATION OVERLAP IN (a,b) SPACE")
print("=" * 80)

print("""
The CKM matrix is the mismatch between up-type and down-type mass eigenstates.
We compute various "distance" metrics between up_i and down_j in (a,b) space
and check which gives CKM-like mixing.
""")

# 1a: Euclidean distance
print("--- 1a: Euclidean distance in (a,b) ---")
v_up = np.array([[a_up[i], b_up[i]] for i in range(3)], dtype=float)
v_dn = np.array([[a_dn[i], b_dn[i]] for i in range(3)], dtype=float)

dist_euclid = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        dist_euclid[i, j] = np.linalg.norm(v_up[i] - v_dn[j])

print(f"  Euclidean distance |v_up - v_dn|:")
print(f"        d          s          b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {dist_euclid[i,j]:8.4f}"
    print(row)

# 1b: Icosian Gram metric
print("\n--- 1b: Icosian Gram metric distance ---")
G_icos = np.array([[1, 1], [1, 2]])  # Gram from exp120d

dist_icos = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        dv = v_up[i] - v_dn[j]
        d2 = dv @ G_icos @ dv
        dist_icos[i, j] = np.sqrt(max(0, d2))

print(f"  Icosian distance sqrt((da,db)*G*(da,db)^T):")
print(f"        d          s          b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {dist_icos[i,j]:8.4f}"
    print(row)

# 1c: Intersection number weighted distance
print("\n--- 1c: Intersection-number distance |5*da + 6*db| = |dn| ---")
print(f"  (This is just the level difference matrix)")
print(f"        d          s          b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {dn_matrix[i,j]:8d}"
    print(row)

# 1d: Check diagonal vs off-diagonal structure
print("\n--- 1d: CKM structure analysis ---")
print(f"  The CKM should have: diagonal ~ 1, off-diagonal ~ small")
print(f"  For this, the 'natural pairing' matters: which up pairs with which down?")
print(f"  Generation pairing: (u,d), (c,s), (t,b)")
print(f"  Diagonal distances (gen pairing):")
for i in range(3):
    print(f"    ({labels_up[i]},{labels_dn[i]}): Euclid={dist_euclid[i,i]:.4f}, "
          f"Icos={dist_icos[i,i]:.4f}, |dn|={dn_matrix[i,i]}")

print(f"\n  Off-diagonal distances (mixing):")
for i in range(3):
    for j in range(3):
        if i != j:
            print(f"    ({labels_up[i]},{labels_dn[j]}): Euclid={dist_euclid[i,j]:.4f}, "
                  f"Icos={dist_icos[i,j]:.4f}, |dn|={dn_matrix[i,j]}")

# Best fit: V_ij ~ phi^(-alpha * d_ij) / normalization
print("\n--- 1e: Best fits V_ij ~ phi^(-alpha * distance) ---")

for dist_name, dist_mat in [("Euclidean", dist_euclid), ("Icosian", dist_icos),
                              ("Level |dn|", dn_matrix.astype(float))]:
    best_err = 1e10
    best_alpha = 0
    best_V = None
    for alpha in np.linspace(0.01, 5.0, 5000):
        V_pred = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                V_pred[i, j] = PHI**(-alpha * dist_mat[i, j])
            V_pred[i] /= np.linalg.norm(V_pred[i])
        err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
        if err < best_err:
            best_err = err
            best_alpha = alpha
            best_V = V_pred.copy()
    print(f"\n  {dist_name}: best alpha={best_alpha:.4f}, mean err={best_err:.2f}%")
    # Check if alpha is close to a geometric quantity
    print(f"    alpha = {best_alpha:.4f}")
    print(f"    alpha*5 = {best_alpha*5:.4f}, alpha*6 = {best_alpha*6:.4f}")
    print(f"    alpha*11 = {best_alpha*11:.4f} (b_1*a_1+1?)")
    print(f"    1/alpha = {1/best_alpha:.4f}")
    if dist_name == "Level |dn|":
        compute_errors(best_V, f"1e-{dist_name}: phi^(-{best_alpha:.3f}*|dn|)")


# ================================================================
# APPROACH 2: Level difference pattern
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 2: LEVEL DIFFERENCE PATTERN")
print("=" * 80)

print("""
Key observation from the MEMORY:
  |n_u - n_d| = |3-5|   = 2
  |n_c - n_s| = |16-11| = 5  (= diameter of 600-cell!)
  |n_t - n_b| = |26-19| = 7  (= H4' Coxeter exponent!)

But for CKM angles we need off-diagonal too:
  |n_u - n_s| = |3-11|  = 8
  |n_c - n_b| = |16-19| = 3
  |n_u - n_b| = |3-19|  = 16
""")

print("--- 2a: Check {11, 6, 7} for diagonal differences ---")
diag_diffs = [dn_matrix[0,0], dn_matrix[1,1], dn_matrix[2,2]]
print(f"  Diagonal |dn|: u-d={diag_diffs[0]}, c-s={diag_diffs[1]}, t-b={diag_diffs[2]}")
print(f"  --> {diag_diffs} = {{2, 5, 7}}")
print(f"  2 = number of (a,b) dims")
print(f"  5 = a_1 = diameter = H4 rank+1")
print(f"  7 = H4' Coxeter exponent")

print("\n--- 2b: CKM mixing angles from specific level differences ---")
# The primary mixing elements:
# V_us ~ sin(theta_12): driven by u-s overlap, |n_u - n_s| = 8
# V_cb ~ sin(theta_23): driven by c-b overlap, |n_c - n_b| = 3
# V_ub ~ sin(theta_13): driven by u-b overlap, |n_u - n_b| = 16

dn_us = abs(n_up[0] - n_dn[1])   # 8
dn_cb = abs(n_up[1] - n_dn[2])   # 3
dn_ub = abs(n_up[0] - n_dn[2])   # 16

print(f"  Primary CKM transitions:")
print(f"    theta_12 ~ V_us: |n_u - n_s| = {dn_us}")
print(f"    theta_23 ~ V_cb: |n_c - n_b| = {dn_cb}")
print(f"    theta_13 ~ V_ub: |n_u - n_b| = {dn_ub}")
print(f"    Note: {dn_ub} = {dn_us} + {dn_us} = 2*{dn_us}? No, {dn_ub} = {dn_us}+{dn_ub-dn_us}")
print(f"    But: {dn_ub} = {dn_us} * 2 = 16? Yes!")
print(f"    And: {dn_us} = {dn_cb} + {dn_us - dn_cb}")
print(f"    Hierarchy: {dn_cb} < {dn_us} < {dn_ub} => theta_23 > theta_12 > theta_13? NO!")
print(f"    Wait: sin(theta) = phi^(-dn/k), so LARGER dn => SMALLER angle")
print(f"    So {dn_cb}=3 < {dn_us}=8 < {dn_ub}=16 would give theta_23 > theta_12 > theta_13")
print(f"    Experimentally: theta_12 > theta_23 > theta_13. MISMATCH for theta_12 vs theta_23!")

print("\n--- 2c: What if CKM angles are NOT from these level differences? ---")
# Alternative: the CKM angles come from a combination of (a,b) geometry
# Maybe the angle between up and down (a,b) vectors matters
print(f"  Alternative: use the (a,b) overlap directly")

# For each generation pair, compute the angle between (a,b) vectors
for i in range(3):
    for j in range(3):
        dot = a_up[i]*a_dn[j] + b_up[i]*b_dn[j]
        n1 = np.sqrt(a_up[i]**2 + b_up[i]**2)
        n2 = np.sqrt(a_dn[j]**2 + b_dn[j]**2)
        if n1 > 0 and n2 > 0:
            cos_ang = dot / (n1 * n2)
            cos_ang = np.clip(cos_ang, -1, 1)
            ang = np.degrees(np.arccos(cos_ang))
        else:
            ang = 0
        if i == j or (i, j) in [(0, 1), (1, 2), (0, 2)]:
            tag = " <--" if i != j else ""
            print(f"    angle({labels_up[i]},{labels_dn[j]}): dot={dot}, angle={ang:.1f} deg{tag}")

print("\n--- 2d: Try sin(theta_ij) = phi^(-|dn_ij|/k) for various k ---")
# The key question: what divides the level differences to give CKM angles?
# From known: sin(theta_12) = phi^(-3.10) but |n_u-n_s| = 8
# So k = 8/3.10 = 2.58
# sin(theta_23) = phi^(-6.62) but |n_c-n_b| = 3
# So k = 3/6.62 = 0.453 -- VERY different!
# This means a UNIFORM divisor does NOT work.

print(f"  Required divisors k_ij for sin(theta) = phi^(-|dn|/k):")
print(f"    theta_12: |dn|=8, n_exact={n_exp_12:.4f}, k=8/{n_exp_12:.4f} = {8/n_exp_12:.4f}")
print(f"    theta_23: |dn|=3, n_exact={n_exp_23:.4f}, k=3/{n_exp_23:.4f} = {3/n_exp_23:.4f}")
print(f"    theta_13: |dn|=16, n_exact={n_exp_13:.4f}, k=16/{n_exp_13:.4f} = {16/n_exp_13:.4f}")
print(f"  --> Non-uniform: 2.58, 0.45, 1.38. No simple pattern.")
print(f"  CONCLUSION: Level differences alone do NOT determine CKM. Need (a,b) structure.")


# ================================================================
# APPROACH 3: H4 WEYL GROUP ACTION
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 3: H4 WEYL GROUP ACTION")
print("=" * 80)

print("""
The H4 Weyl group |W(H4)| = 14400 acts on the 120 vertices of the 600-cell.
If up-type and down-type quarks live on different orbits under a SUBGROUP
of W(H4), the CKM matrix is the overlap between these subgroup representations.
""")

# H4 simple roots (unit norm)
alpha1 = np.array([1.0, 0.0, 0.0, 0.0])
c2 = np.sqrt(1 - (PHI/2)**2)
alpha2 = np.array([-PHI/2, c2, 0.0, 0.0])
b3 = -1.0 / (2.0 * c2)
c3_val = np.sqrt(max(0, 1.0 - b3**2))
alpha3 = np.array([0.0, b3, c3_val, 0.0])
c4_val = -1.0 / (2.0 * c3_val)
d4_val = np.sqrt(max(0, 1.0 - c4_val**2))
alpha4 = np.array([0.0, 0.0, c4_val, d4_val])

simple_roots = [alpha1, alpha2, alpha3, alpha4]

def reflection(root):
    r = root / np.linalg.norm(root)
    return np.eye(4) - 2 * np.outer(r, r)

S_refl = [reflection(r) for r in simple_roots]

# Coxeter element
Cox = np.eye(4)
for S_i in S_refl:
    Cox = Cox @ S_i

# Verify order
Cox_power = np.eye(4)
for k in range(1, 35):
    Cox_power = Cox_power @ Cox
    if np.allclose(Cox_power, np.eye(4), atol=1e-8):
        cox_order = k
        break

print(f"  H4 Coxeter element order: {cox_order}")

# Eigenvalues of Coxeter element
evals_cox = np.linalg.eigvals(Cox)
for ev in evals_cox:
    angle = np.angle(ev) * 180 / np.pi
    if angle < 0:
        angle += 360
    m_val = angle / 12
    print(f"    eigenvalue: {ev.real:+.6f}{ev.imag:+.6f}i, angle={angle:.1f} deg, m={m_val:.1f}")

# The Coxeter element generates a Z_30 subgroup
# Under this subgroup, 120 vertices split into orbits
print(f"\n  Z_30 orbits of 600-cell vertices under Coxeter element:")

def apply_matrix(M, v):
    return M @ v

# Generate Z_30 orbit for vertex 0
orbits = []
assigned = set()
for start_idx in range(N600):
    if start_idx in assigned:
        continue
    orbit = set()
    v = V600[start_idx].copy()
    for k in range(30):
        # Find which vertex this is
        dists = np.linalg.norm(V600 - v, axis=1)
        closest = np.argmin(dists)
        if dists[closest] < 1e-6:
            orbit.add(closest)
        v = Cox @ v
    orbits.append(orbit)
    assigned.update(orbit)

print(f"  Number of Z_30 orbits: {len(orbits)}")
orbit_sizes = sorted([len(o) for o in orbits])
print(f"  Orbit sizes: {orbit_sizes}")

# Now check: do up-type and down-type quarks correspond to different orbits?
# For this, we need to associate each vertex with a "quantum number"
# The quantum numbers (a,b) are GLOBAL properties, not per-vertex.
# But each vertex has its own coordinate decomposition in Z[phi].


# ================================================================
# APPROACH 4: E8 BRANCHING - SPECTRAL APPROACH
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 4: E8 BRANCHING RULES (SPECTRAL)")
print("=" * 80)

print("""
E8 Coxeter exponents: {1, 7, 11, 13, 17, 19, 23, 29}
H4:  {1, 11, 19, 29}
H4': {7, 13, 17, 23}

The key structural relationship:
  CKM_12 = 3 = (1+11)/4    <- sum of SMALLEST H4 pair / rank
  CKM_23 = 7 = 7/1         <- SMALLEST H4' exponent
  CKM_13 = 12 = (19+29)/4  <- sum of LARGEST H4 pair / rank

The half-integer version:
  CKM_12 = 3 = (1+11)/4    <- same
  CKM_23 = 13/2 = 13/2     <- SECOND H4' exponent / 2
  CKM_13 = 23/2 = 23/2     <- THIRD H4' exponent / 2
""")

h4_exp = [1, 11, 19, 29]
h4c_exp = [7, 13, 17, 23]
e8_exp = sorted(h4_exp + h4c_exp)

print("--- 4a: NEW observation: E8 exponent differences ---")
# H4' - H4 differences (sorted by rank)
diffs = [h4c_exp[i] - h4_exp[i] for i in range(4)]
print(f"  H4' - H4 (rank-paired): {diffs} = [6, 2, -2, -6]")
print(f"  Antisymmetric: sum = {sum(diffs)}")
print(f"  |diffs| = {[abs(d) for d in diffs]} = [6, 2, 2, 6]")
print(f"  6 = b_1 (intersection number)")
print(f"  2 = first diff of H4 exps (11-1)/5 = 2")

print("\n--- 4b: CKM from H4-H4' cross-products ---")
# What if the mixing angles come from the cross-correlation between
# H4 and H4' Coxeter representations?
# The branching E8 -> H4 x H4' gives: 248 = 2*(1,120) + (1,8)
# The CKM might encode how the 120 decomposes under the different choices

# Try: CKM angle_ij = phi^(-(m_H4[i] * m_H4'[j]) / normalization)
print(f"\n  Products m_H4 * m_H4':")
for i, m4 in enumerate(h4_exp):
    for j, m4c in enumerate(h4c_exp):
        print(f"    {m4} * {m4c} = {m4*m4c}")

print(f"\n  Products/30 (= products/h):")
for i, m4 in enumerate(h4_exp):
    for j, m4c in enumerate(h4c_exp):
        val = m4 * m4c / 30
        print(f"    ({m4}*{m4c})/30 = {val:.3f}", end="")
        if abs(val - round(val)) < 0.01:
            print(f" = {int(round(val))}", end="")
        print()

print("\n--- 4c: Coxeter exponent sums that give CKM angles ---")
# From exp125, the best exponent triples found:
# (3, 13/2, 23/2) with mean 6%
# (3, 7, 12) with mean ~15%

# NEW IDEA: What if the formula is:
# n_12 = (m_1 + m_2) / rank where m_1, m_2 are the PAIRED exponents
# involving the lightest generation?

# Pair H4 exponents: (1,29) and (11,19) [sum-to-h pairs]
# Pair H4' exponents: (7,23) and (13,17) [sum-to-h pairs]

print(f"\n  H4 sum-to-h pairs: (1,29) and (11,19)")
print(f"  H4' sum-to-h pairs: (7,23) and (13,17)")
print(f"\n  Pair sums / rank(H4)=4:")
print(f"    (1+11)/4 = {(1+11)/4}")
print(f"    (1+19)/4 = {(1+19)/4}")
print(f"    (1+29)/4 = {(1+29)/4} = h/4")
print(f"    (11+19)/4 = {(11+19)/4} = h/4")
print(f"    (11+29)/4 = {(11+29)/4}")
print(f"    (19+29)/4 = {(19+29)/4}")

# The CKM exponents 3, 7.5, 12 correspond to:
# 3 = (1+11)/4 : SMALLEST pair of DISTINCT h-complement pairs
# 12 = (19+29)/4 : LARGEST pair of DISTINCT h-complement pairs
# 7.5 = (11+19)/4 = 30/4 = h/4 : MIDDLE pair (same sum as (1,29))
# But h/4 = 7.5 != 7 or 6.5!

# Let's try: n_23 from CROSS-PAIRING (H4 x H4'):
print(f"\n  Cross-pairing sums (H4 + H4') / 4:")
for m4 in h4_exp:
    for m4c in h4c_exp:
        val = (m4 + m4c) / 4
        mark = ""
        if abs(val - 3) < 0.1: mark = " <-- n_12?"
        elif abs(val - n_exp_23) < 0.5: mark = f" <-- close to n_23={n_exp_23:.2f}"
        elif abs(val - n_exp_13) < 0.5: mark = f" <-- close to n_13={n_exp_13:.2f}"
        elif abs(val - 7) < 0.1: mark = " <-- n_23(integer)?"
        elif abs(val - 12) < 0.1: mark = " <-- n_13(integer)?"
        print(f"    ({m4}+{m4c})/4 = {val:6.2f}{mark}")


# ================================================================
# APPROACH 5: INTERSECTION NUMBERS APPROACH
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 5: INTERSECTION NUMBERS APPROACH")
print("=" * 80)

print("""
Intersection numbers: a_1 = 5, b_1 = 6
Level formula: n = 5a + 6b

Hypothesis: V_ij ~ phi^(-(a_1*|da_ij| + b_1*|db_ij|) / divisor)
where da = a_up_i - a_dn_j, db = b_up_i - b_dn_j.
""")

# Compute the intersection-weighted Manhattan distance
int_dist = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        int_dist[i, j] = 5 * abs(da_matrix[i, j]) + 6 * abs(db_matrix[i, j])

print(f"  5*|da| + 6*|db| matrix:")
print(f"        d          s          b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {int_dist[i,j]:8.0f}"
    print(row)

print(f"\n  Compare with spectral distance |5*da + 6*db| = |dn|:")
print(f"  Note: 5*|da|+6*|db| >= |5*da+6*db| by triangle inequality")
print(f"  Equality only when da and db have the same sign in the 5a+6b combination")

# Scan divisors
print(f"\n--- 5a: CKM from sin(theta) = phi^(-(5|da|+6|db|)/k) ---")

# For the CKM primary mixing:
# V_us: da = u_a - s_a = 3-1=2, db = u_b - s_b = -2-1=-3
#   5*|2| + 6*|-3| = 10+18 = 28
# V_cb: da = c_a - b_a = 2-(-1)=3, db = c_b - b_b = 1-4=-3
#   5*|3| + 6*|-3| = 15+18 = 33
# V_ub: da = u_a - b_a = 3-(-1)=4, db = u_b - b_b = -2-4=-6
#   5*|4| + 6*|-6| = 20+36 = 56

ij_us = int_dist[0, 1]  # 28
ij_cb = int_dist[1, 2]  # 33
ij_ub = int_dist[0, 2]  # 56

print(f"  5|da|+6|db| for CKM angles: us={ij_us:.0f}, cb={ij_cb:.0f}, ub={ij_ub:.0f}")
print(f"  Required k: us/{n_exp_12:.2f}={ij_us/n_exp_12:.2f}, "
      f"cb/{n_exp_23:.2f}={ij_cb/n_exp_23:.2f}, ub/{n_exp_13:.2f}={ij_ub/n_exp_13:.2f}")
print(f"  k values: {ij_us/n_exp_12:.2f}, {ij_cb/n_exp_23:.2f}, {ij_ub/n_exp_13:.2f}")
print(f"  These are very different => 5|da|+6|db| does NOT give uniform CKM formula.")

# Alternative: use the SIGNED version 5*da + 6*db (with sign)
print(f"\n--- 5b: Signed distances (5*da + 6*db) ---")
signed_dist = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        signed_dist[i, j] = 5 * da_matrix[i, j] + 6 * db_matrix[i, j]

print(f"  5*da + 6*db matrix (signed = n_up - n_dn):")
print(f"        d          s          b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {signed_dist[i,j]:8.0f}"
    print(row)

# NEW: Try with SEPARATE a and b contributions
print(f"\n--- 5c: Try V_ij = phi^(-(alpha*|da| + beta*|db|)) ---")
best_5c_err = 1e10
best_5c_params = (0, 0)
best_5c_V = None

for alpha_val in np.linspace(0.01, 3.0, 200):
    for beta_val in np.linspace(0.01, 3.0, 200):
        V_pred = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                d = alpha_val * abs(da_matrix[i, j]) + beta_val * abs(db_matrix[i, j])
                V_pred[i, j] = PHI**(-d)
            V_pred[i] /= np.linalg.norm(V_pred[i])
        err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
        if err < best_5c_err:
            best_5c_err = err
            best_5c_params = (alpha_val, beta_val)
            best_5c_V = V_pred.copy()

print(f"  Best: alpha={best_5c_params[0]:.4f}, beta={best_5c_params[1]:.4f}")
print(f"  Ratio alpha/beta = {best_5c_params[0]/best_5c_params[1]:.4f}")
print(f"  Compare: a_1/b_1 = 5/6 = {5/6:.4f}")
print(f"  Compare: 1/phi = {1/PHI:.4f}")
m5c, x5c, c5c = compute_errors(best_5c_V, "5c: phi^(-(alpha*|da|+beta*|db|))")


# ================================================================
# APPROACH 6: DECAGON WINDING NUMBERS
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 6: DECAGON WINDING NUMBERS")
print("=" * 80)

print("""
The 600-cell contains 72 great decagons (10-gons inscribed in S^3).
Each decagon wraps around S^3 and has a definite phase/winding.
Quarks at different levels may sit on different decagons, and the
CKM mixing angle might be related to the relative winding between
the up-type and down-type decagons.
""")

# Find all decagons (10-cycles in the graph)
# A decagon in 600-cell: 10 vertices equally spaced along a great circle
# Each vertex has 12 neighbors. From any vertex, the great circles through
# its neighbors define the decagons.

# For a great decagon, consecutive vertices have dot product = cos(36 deg) = phi/2
# and every 2nd vertex has dot product = cos(72 deg) = (phi-1)/2

# Find one decagon starting from vertex 0
print("--- 6a: Finding great decagons ---")

# Method: start from vertex 0, pick neighbor, follow the "straight" geodesic
# On S^3, a great circle is determined by 2 points.
# Great decagon: equally spaced vertices with angular separation 2*pi/10 = 36 deg

decagons_found = 0
decagons = []
for start_neighbor_idx in range(N600):
    if abs(dots[0, start_neighbor_idx] - cos_nn) > 0.001:
        continue
    # Follow the geodesic from 0 through this neighbor
    decagon = [0, start_neighbor_idx]
    for step in range(8):
        current = decagon[-1]
        prev = decagon[-2]
        # The next vertex should be: a neighbor of current that continues the geodesic
        # On a great decagon, the pattern of dot products is periodic
        # Next vertex: neighbor of current, NOT prev, with correct angle to prev
        # Dot(next, prev) should be cos(72 deg) = (phi-1)/2 = 1/(2*phi)
        cos_2step = 1 / (2 * PHI)  # cos(72 deg) = (sqrt(5)-1)/4 = phi/2 - 1/2
        # Actually cos(72 deg) = (sqrt(5)-1)/4 ~= 0.309
        # Let me compute properly: cos(72 deg) = cos(2*pi/5) = (sqrt(5)-1)/4
        cos72 = (np.sqrt(5) - 1) / 4
        candidates = []
        for k in range(N600):
            if k == prev or k == current:
                continue
            if abs(dots[current, k] - cos_nn) < 0.001:
                if abs(dots[prev, k] - cos72) < 0.01:
                    candidates.append(k)
        if len(candidates) == 1:
            decagon.append(candidates[0])
        elif len(candidates) > 1:
            # Pick the one that also matches the pattern from 2 steps back
            decagon.append(candidates[0])
        else:
            break

    if len(decagon) == 10:
        # Verify it closes: last should be neighbor of first with correct angle
        if abs(dots[decagon[-1], decagon[0]] - cos_nn) < 0.001:
            # Check it's a new decagon (not a rotation of existing)
            dec_set = frozenset(decagon)
            if dec_set not in [frozenset(d) for d in decagons]:
                decagons.append(decagon)

print(f"  Found {len(decagons)} great decagons (expect 72)")

# Even if we don't find all 72, analyze the ones we have
if len(decagons) > 0:
    print(f"  First decagon vertices: {decagons[0]}")
    dec0 = decagons[0]
    print(f"  Dot products along decagon:")
    for k in range(len(dec0)):
        next_k = (k + 1) % len(dec0)
        print(f"    V[{dec0[k]}] . V[{dec0[next_k]}] = {dots[dec0[k], dec0[next_k]]:.6f} "
              f"(expect {cos_nn:.6f})")

    # How many decagons pass through vertex 0?
    through_0 = [d for d in decagons if 0 in d]
    print(f"  Decagons through vertex 0: {len(through_0)} (expect 6)")

# Alternative: compute decagon winding via eigenvalues
# The decagons correspond to 2D irreps of the icosahedral group
# Their winding numbers are related to the Coxeter exponents
print(f"\n--- 6b: Decagon phase from Coxeter element ---")
print(f"  Each decagon traces out a phase of 2*pi/10 * winding_number per step")
print(f"  Coxeter element rotates by 12 deg and 132 deg in the two planes")
print(f"  12 deg = 2*pi/30 -> 30 steps for full rotation -> winding 3 for decagon (30/10)")
print(f"  132 deg = 11*12 deg -> winding 11 for second plane")
print(f"  The winding ratio is 11/3 = 11/3")
print(f"  11 and 3 are related to: 11 = H4 exp, 3 = Cabibbo exponent!")


# ================================================================
# APPROACH 7: GOLDEN ANGLE SUBDIVISIONS
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 7: GOLDEN ANGLE SUBDIVISIONS")
print("=" * 80)

print("""
The golden angle = 2*pi*(1 - 1/phi^2) = 2*pi/phi^2 = 137.508 deg.
This is related to the Cabibbo angle 13.04 deg by:
  137.508 / 13.04 = 10.55 ~ 10 + 1/phi = 10.618
Also: 137.508 / phi^5 = 137.508 / 11.09 = 12.40 deg (close to 12 deg = 2*pi/30)
""")

golden_angle = 2 * np.pi / PHI**2
print(f"  Golden angle = {np.degrees(golden_angle):.3f} deg = 2*pi/phi^2")
print(f"  Cabibbo angle (exp) = {np.degrees(np.arcsin(0.2243)):.3f} deg")
print(f"  Ratio: {np.degrees(golden_angle) / np.degrees(np.arcsin(0.2243)):.4f}")

# The CKM angles as subdivisions of the golden angle
theta_C = np.arcsin(0.2243)
theta_23 = np.arcsin(0.0408)
theta_13 = np.arcsin(0.00382)

for name, theta in [("theta_12", theta_C), ("theta_23", theta_23), ("theta_13", theta_13)]:
    ratio = golden_angle / theta
    log_phi_ratio = np.log(ratio) / np.log(PHI)
    print(f"  golden_angle / {name} = {ratio:.4f} = phi^{log_phi_ratio:.4f}")

# Check if CKM angles = golden_angle / phi^n
print(f"\n  If theta_ij = golden_angle / phi^n:")
for n in range(1, 16):
    val = np.degrees(golden_angle / PHI**n)
    print(f"    n={n:2d}: {val:10.4f} deg", end="")
    for name, theta in [("theta_12", theta_C), ("theta_23", theta_23), ("theta_13", theta_13)]:
        if abs(val - np.degrees(theta)) / np.degrees(theta) < 0.1:
            print(f"  <-- close to {name} ({abs(val-np.degrees(theta))/np.degrees(theta)*100:.1f}%)", end="")
    print()


# ================================================================
# APPROACH 8: DIRECT 600-CELL EIGENVECTOR OVERLAP
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 8: DIRECT 600-CELL EIGENVECTOR OVERLAP")
print("=" * 80)

print("""
This is the most direct approach: use the actual 600-cell adjacency matrix
eigenspaces to compute the CKM as a 3x3 overlap matrix.

Each occupied level n corresponds to a specific eigenvalue of the 600-cell.
The eigenspace defines a SUBSPACE of R^120. The up-type and down-type
quarks each pick out specific subspaces, and the CKM matrix is the
overlap (principal angles) between these subspaces.

The 9 eigenvalues with multiplicities n^2:
  lambda_0 = 12   (mult=1=1^2)
  lambda_1 = 6*phi (mult=4=2^2)   eigenvalue ~ 9.708
  lambda_2 = 4*phi (mult=9=3^2)   eigenvalue ~ 6.472
  lambda_3 = 3     (mult=16=4^2)
  lambda_4 = 0     (mult=25=5^2)
  lambda_5 = -2    (mult=36=6^2)
  lambda_6 = 4-4phi(mult=9=3^2)   eigenvalue ~ -2.472
  lambda_7 = -3    (mult=16=4^2)
  lambda_8 = 6-6phi(mult=4=2^2)   eigenvalue ~ -3.708
""")

# Map each spectral level n to an eigenvalue index
# From MEMORY: the 9 eigenvalues are theta = a + b*phi
# Eigenvalues: 12, 6phi, 4phi, 3, 0, -2, 4-4phi, -3, 6-6phi
# Levels: 0, 3, 5, 11, 16, 17, 19, 26 are occupied

# The level n corresponds to the eigenvalue with index determined by
# the association scheme. The distance-j matrix A_j has eigenvalue
# P_{j,k} in the k-th eigenspace.

# For now, let's work with the eigenspaces directly.
# Each occupied mass level n maps to some combination of eigenspace projections.

# From the distance-regular structure, vertex v at graph distance d from vertex 0
# contributes to the d-th "shell". The eigenspace decomposition is:
# P_0 = all 1s (eigenvalue 12)
# P_k for k=1..8: projector onto k-th eigenspace

# For each vertex v, its "spectral fingerprint" is:
# (c_0(v), c_1(v), ..., c_8(v)) where c_k(v) is the component in eigenspace k.

# Since the graph is vertex-transitive, all vertices have the same spectral
# decomposition. But the SHELL STRUCTURE differs: vertices at distance d
# from a reference have specific spectral weights.

# Let's compute the P-matrix (association scheme eigenmatrix)
print("--- 8a: Association scheme P-matrix ---")

# The distance matrices A_0, A_1, ..., A_8
dist_classes = []
unique_cos_sorted = sorted(set(np.round(dots[0], 6)))
print(f"  Number of distance classes: {len(unique_cos_sorted)}")

# Build distance class matrices
dist_matrices = []
for cos_val in unique_cos_sorted:
    D = (np.abs(dots - cos_val) < 0.001).astype(float)
    if cos_val == unique_cos_sorted[-1]:  # cos(0) = 1, diagonal
        D = np.eye(N600)
    dist_matrices.append(D)

# Verify: sum should be all-ones
check_sum = sum(dist_matrices)
print(f"  Sum of distance matrices = all-ones: {np.allclose(check_sum, np.ones((N600, N600)))}")

# Compute P-matrix: P[d][k] = eigenvalue of A_d in eigenspace k
# We already have the eigendecomposition. For each distinct eigenvalue,
# project A_d and read off the eigenvalue.

print(f"\n  Computing P-matrix (distance class d vs eigenspace k)...")

# Get eigenspace projectors
eigenspaces = []  # list of (eigenvalue, list_of_eigvec_indices)
idx_pos = 0
for k, (ev, m) in enumerate(zip(unique_evals, mults)):
    eigenspaces.append((ev, list(range(idx_pos, idx_pos + m))))
    idx_pos += m

P_matrix = np.zeros((len(dist_matrices), len(unique_evals)))
for d, Ad in enumerate(dist_matrices):
    for k, (ev, idxs) in enumerate(eigenspaces):
        # Project Ad onto eigenspace k: should give scalar * I
        # The eigenvalue of Ad in eigenspace k:
        # Trace(Ad * Pk) / dim(Pk) where Pk = sum_i |vi><vi|
        v0 = eigenvectors_full[:, idxs[0]]
        P_matrix[d, k] = v0 @ Ad @ v0
        # Verify with second eigenvector if available
        if len(idxs) > 1:
            v1 = eigenvectors_full[:, idxs[1]]
            check = v1 @ Ad @ v1
            if abs(check - P_matrix[d, k]) > 0.1:
                P_matrix[d, k] = np.nan  # Inconsistent

print(f"\n  P-matrix (rows=distance class, cols=eigenspace):")
print(f"  {'d':>3s}", end="")
for k in range(len(unique_evals)):
    print(f"  {'lam'+str(k):>8s}", end="")
print()
for d in range(len(dist_matrices)):
    count_d = int(dist_matrices[d][0].sum())
    print(f"  {d:3d} (n={count_d:3d})", end="")
    for k in range(len(unique_evals)):
        print(f"  {P_matrix[d,k]:8.3f}", end="")
    print()

# The Q-matrix (dual eigenmatrix): Q[k][d] = mult_k * P[d][k] / n_d
# This tells us how eigenspace k decomposes into distance classes.

print(f"\n--- 8b: Eigenspace decomposition into shells ---")
print(f"  Each eigenspace k has a specific 'radial profile' across distance shells.")
print(f"  This profile defines how a state at level k spreads on the 600-cell.")

# For the CKM computation, we need to know which eigenspace each
# mass level belongs to. The mass levels are:
occupied_levels = {0: 'e', 3: 'u', 5: 'd', 11: 'mu/s', 16: 'c',
                   17: 'tau', 19: 'b', 26: 't'}

# The assignment level -> eigenspace is: n = 5a + 6b = dot product
# with the intersection numbers (5, 6). But this is a SPECTRAL label,
# not a distance class index.

# The key insight: n is not a distance class but a SPECTRAL LABEL.
# It tells us which LINEAR COMBINATION of eigenspaces to use.
# From exp111: n = a_1*a + b_1*b = 5a + 6b determines the fermion.
# The eigenvalue associated with level n is computed via the mass formula.

# What we need: for each occupied level, identify the EIGENVECTOR(s)
# in the 600-cell spectrum that encode the quark wavefunction.

# Approach: the "wavefunction" at level n is:
# psi_n(v) = sum_k c_k(n) * phi_k(v) where phi_k is the k-th eigenfunction
# and c_k(n) is determined by the mass formula.

# For a distance-regular graph with association scheme, the eigenfunctions
# at vertex v are determined by the distance from a reference vertex.
# So we can define the CKM as the overlap between different reference frames.

print(f"\n--- 8c: CKM from eigenvector inner products ---")
print(f"  Taking a different approach: use eigenvectors DIRECTLY.")

# For each eigenvalue lambda_k with multiplicity m_k, we have m_k
# eigenvectors. These form a representation of the symmetry group.
# The CKM might be the overlap between up-type and down-type
# eigenvector subsets.

# The occupied levels (0,3,5,11,16,17,19,26) map to eigenspaces.
# But which eigenspace does level n map to?

# The eigenvalues in order are approximately:
# 12.0, 9.71, 6.47, 3.0, 0.0, -2.0, -2.47, -3.0, -3.71
# with multiplicities 1, 4, 9, 16, 25, 36, 9, 16, 4

# From exp110, the eigenvalues are:
# theta_k = a_k + b_k * phi where
# (a_k, b_k) = (12,0), (0,6), (0,4), (3,0), (0,0), (-2,0), (4,-4), (-3,0), (6,-6)

# The level n = 5a + 6b for the eigenvalue a_k + b_k*phi gives:
# n_k = 5*a_k + 6*b_k = {60, 36, 24, 15, 0, -10, -4, -15, -6}
# These DON'T match the occupied levels {0,3,5,11,16,17,19,26}!

# So the mass levels are NOT the eigenvalues' own (a,b) decomposition.
# The (a,b) quantum numbers are SEPARATE from the eigenvalue (a,b) coefficients.

# NEW IDEA: The wavefunctions are the COLUMNS of the eigenvector matrix.
# Each column is an eigenvector. The 120 components correspond to 120 vertices.
# If we pick 3 eigenvectors for up-type and 3 for down-type,
# their overlap gives a 3x3 matrix.

# But WHICH eigenvectors to pick? There are 120 total.
# The multiplicities suggest grouping: the k-th eigenspace has mult = k^2
# for k = 1..6, then 3^2, 4^2, 2^2.

# Idea: within each eigenspace, the eigenvectors transform as irreps of H4.
# Multiplicity n^2 suggests n copies of an n-dimensional irrep (or n^2-dim irrep).
# For mult=4=2^2: 2 copies of 2-dim irrep.
# For mult=9=3^2: 3 copies of 3-dim irrep (or 1 copy of 9-dim).

# CKM involves 3 generations -> use eigenspace with mult = 9 (3^2)
# There are two such eigenspaces: lambda_2 (4*phi) and lambda_6 (4-4phi)

print(f"\n  Eigenspaces with multiplicity 9 (=3^2):")
for k, (ev, m) in enumerate(zip(unique_evals, mults)):
    if m == 9:
        print(f"    lambda_{k} = {ev:.4f}, mult = {m}")

# Extract the two 9-dimensional eigenspaces
spaces_9 = []
for k, (ev, m) in enumerate(zip(unique_evals, mults)):
    if m == 9:
        idxs = eigenspaces[k][1]
        vecs = eigenvectors_full[:, idxs]
        spaces_9.append((ev, vecs, k))

if len(spaces_9) >= 2:
    ev1, vecs1, k1 = spaces_9[0]
    ev2, vecs2, k2 = spaces_9[1]

    print(f"\n  Eigenspace k={k1} (lambda={ev1:.4f}): 9 eigenvectors")
    print(f"  Eigenspace k={k2} (lambda={ev2:.4f}): 9 eigenvectors")

    # The 9-dimensional space can be organized as a 3x3 matrix.
    # If it decomposes as 3 copies of a 3-dim H4 irrep,
    # we can build a 3x3 CKM from the overlap.

    # Try: split each 9-dim space into three 3-dim subspaces
    # and compute the overlap between pairs.

    # Method: SVD of the overlap matrix between the two eigenspaces
    # restricted to a chosen set of 3 vertices.

    # Better: use the PRINCIPAL ANGLES between the two 9-dim subspaces
    overlap_9 = vecs1.T @ vecs2  # 9x9 matrix
    U9, S9, V9t = np.linalg.svd(overlap_9)
    print(f"\n  Principal values (cosines of canonical angles):")
    print(f"    {np.round(S9, 6)}")
    print(f"    These are the cosines of the 9 canonical angles between the two eigenspaces.")

    # The canonical angles theta_k satisfy cos(theta_k) = S9[k]
    canonical_angles = np.arccos(np.clip(S9, -1, 1))
    print(f"    Canonical angles (degrees): {np.round(np.degrees(canonical_angles), 4)}")

    # Do any of these match CKM angles?
    theta_12_exp = np.degrees(np.arcsin(V_exp[0, 1]))
    theta_23_exp = np.degrees(np.arcsin(V_exp[1, 2]))
    theta_13_exp = np.degrees(np.arcsin(V_exp[0, 2]))
    print(f"\n    CKM angles: theta_12={theta_12_exp:.4f} deg, theta_23={theta_23_exp:.4f} deg, theta_13={theta_13_exp:.4f} deg")

    for i, ca in enumerate(np.degrees(canonical_angles)):
        matches = []
        if abs(ca - theta_12_exp) / theta_12_exp < 0.3:
            matches.append("~theta_12")
        if abs(ca - theta_23_exp) / theta_23_exp < 0.3:
            matches.append("~theta_23")
        if abs(ca - theta_13_exp) / theta_13_exp < 0.3:
            matches.append("~theta_13")
        if matches:
            print(f"    Canonical angle {i}: {ca:.4f} deg --> {', '.join(matches)}")

    # Alternative: use one 9-dim space, organize as 3x3 matrix, and compute
    # the CKM from the singular value decomposition
    print(f"\n--- 8d: CKM from 3x3 reshaping of 9-dim eigenspace ---")

    # The 9 eigenvectors can be reshaped as 3 groups of 3.
    # The grouping is not unique -- try all possible groupings.
    # For a quick test, use the natural ordering.

    for kk, (ev_k, vecs_k, k_idx) in enumerate(spaces_9):
        # Reshape 120x9 matrix into three 120x3 matrices
        # Method 1: columns [0,1,2], [3,4,5], [6,7,8]
        M1 = vecs_k[:, 0:3]
        M2 = vecs_k[:, 3:6]
        M3 = vecs_k[:, 6:9]

        # Overlap: O_ij = <M_i | M_j> = M_i^T * M_j (3x3)
        O12 = M1.T @ M2
        U_o, S_o, V_o = np.linalg.svd(O12)
        R_o = U_o @ V_o
        print(f"\n    Eigenspace k={k_idx}: overlap M1-M2 singular values = {np.round(S_o, 4)}")
        print(f"    Rotation matrix |R|:")
        for row in np.abs(R_o):
            print(f"      [{', '.join(f'{x:.5f}' for x in row)}]")


# ================================================================
# APPROACH 9: WOLFENSTEIN DEEP ANALYSIS
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 9: WOLFENSTEIN PARAMETERIZATION DEEP ANALYSIS")
print("=" * 80)

print("""
Wolfenstein parameters:
  lambda = sin(theta_C) ~ 0.2243
  A = |V_cb| / lambda^2 ~ 0.811
  rho + i*eta (CP violation)

From exp122: lambda = phi^(-3) with 1.8% error.
Question: what is A geometrically?
""")

lam = V_exp[0, 1]  # experimental lambda
lam_phi = PHI**(-3)

A_exp = V_exp[1, 2] / lam**2
A_phi = V_exp[1, 2] / lam_phi**2

print(f"  lambda (exp) = {lam:.6f}")
print(f"  lambda (phi^-3) = {lam_phi:.6f}")
print(f"  A (from exp lambda) = {A_exp:.4f}")
print(f"  A (from phi^-3 lambda) = {A_phi:.4f}")

# What is A in terms of phi?
n_A = -np.log(A_phi) / np.log(PHI)
print(f"\n  A = phi^(-{n_A:.4f})")
print(f"  A = {A_phi:.4f}")

# Check various geometric quantities for A
print(f"\n  Checking A against geometric quantities:")
candidates_A = [
    ("1/phi", 1/PHI),
    ("phi^(-0.5)", PHI**(-0.5)),
    ("phi^(-1/phi)", PHI**(-1/PHI)),
    ("2*phi-2", 2*PHI-2),
    ("1/phi^2 * phi", 1/PHI**2 * PHI),
    ("sqrt(phi)/2", np.sqrt(PHI)/2),
    ("phi^2/4", PHI**2/4),
    ("6/b_1*a_1", 6/(6*5)),
    ("sin^2(tW)*phi^2", SIN2_TW * PHI**2),
    ("phi^(-1) * phi^(1/5)", PHI**(-1) * PHI**(1/5)),
    ("1 - 1/phi^3", 1 - 1/PHI**3),
    ("phi^(-1/3)", PHI**(-1/3)),
    ("1/phi^(2/3)", 1/PHI**(2/3)),
    ("3*phi - 4", 3*PHI - 4),
    ("5*phi^(-3)", 5*PHI**(-3)),
    ("alpha_s/sin^2(tW) * phi", ALPHA_S/SIN2_TW * PHI),
    ("(phi+1)/(phi+2)", (PHI+1)/(PHI+2)),
    ("(phi^2-1)/phi^2", (PHI**2-1)/PHI**2),
    ("phi/(phi+1)", PHI/(PHI+1)),
]

for name, val in sorted(candidates_A, key=lambda x: abs(x[1] - A_phi)):
    err = abs(val - A_phi) / A_phi * 100
    if err < 5:
        print(f"    A ~ {name} = {val:.4f} ({err:.2f}%)")

print(f"\n--- 9b: What if A = phi / (phi + 1) = 1/phi^2 + 1/phi? ---")
A_geom = PHI / (PHI + 1)  # = 1/phi = 0.618
# That's too large.
A_geom2 = (PHI**2 - 1) / PHI**2  # = 1/phi^2 = 0.382, too small
# Try: A ~ phi^2 / (phi^2 + 1) = phi^2 / (phi^2+1)
A_geom3 = PHI**2 / (PHI**2 + 1)
print(f"  phi/(phi+1) = {PHI/(PHI+1):.4f} (err {abs(PHI/(PHI+1)-A_phi)/A_phi*100:.1f}%)")
print(f"  phi^2/(phi^2+1) = {A_geom3:.4f} (err {abs(A_geom3-A_phi)/A_phi*100:.1f}%)")

# The actual best geometric candidates
print(f"\n--- 9c: A from Coxeter exponents ---")
# A = V_cb / lambda^2
# V_cb = sin(theta_23) = phi^(-n_23) where n_23 ~ 6.62
# A = phi^(-n_23) / phi^(-6) = phi^(6-n_23) = phi^(-0.62)
# So A = phi^(-0.62)

n23_eff = n_exp_23
print(f"  A = phi^(6 - n_23_eff) = phi^({6-n23_eff:.4f}) = phi^(-{n23_eff-6:.4f})")
print(f"  So A = phi^(-0.62)")
print(f"  0.62 ~ 1/phi = {1/PHI:.4f} (close! diff = {abs(0.62-1/PHI):.4f})")
print(f"  If A = phi^(-1/phi) = {PHI**(-1/PHI):.4f} (err {abs(PHI**(-1/PHI)-A_phi)/A_phi*100:.2f}%)")
print(f"  This is INTERESTING: A = phi^(-1/phi) would be a self-referential formula!")

# Build CKM with lambda = phi^(-3) and A = phi^(-1/phi)
A_golden = PHI**(-1/PHI)
s12_wolf = lam_phi
s23_wolf = A_golden * lam_phi**2
s13_wolf_rhoeta = V_exp[0, 2]  # keep V_ub exact

V_wolf = build_ckm(s12_wolf, s23_wolf, s13_wolf_rhoeta)
m9c, x9c, c9c = compute_errors(V_wolf, "9c: Wolfenstein lambda=phi^(-3), A=phi^(-1/phi)")

# Alternative: n_23 = 13/2 (from H4' exponent)
# Then A = phi^(6-13/2) = phi^(-1/2)
print(f"\n--- 9d: If n_23 = 13/2 (H4' exponent/2), then A = phi^(-1/2) ---")
A_half = PHI**(-0.5)
print(f"  A = phi^(-1/2) = {A_half:.4f}")
s23_half = A_half * lam_phi**2
V_half = build_ckm(lam_phi, s23_half, s13_wolf_rhoeta)
m9d, x9d, c9d = compute_errors(V_half, "9d: lambda=phi^(-3), A=phi^(-1/2)")


# ================================================================
# APPROACH 10: SYNTHESIS - The best derivation possible
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 10: SYNTHESIS - BEST DERIVATION ATTEMPT")
print("=" * 80)

print("""
Combining all approaches, the strongest derivation chain is:

Step 1: sin(theta_12) = phi^(-3)
  WHERE: 3 = (1+11)/4 = (smallest H4 Coxeter exponent pair sum) / rank(H4)
  STATUS: PATTERN (the /4 division is natural but not strictly derived)

Step 2: sin(theta_23) = phi^(-13/2)
  WHERE: 13 is an H4' complement Coxeter exponent
  The /2 comes from: H4' exponents are phase-space partners
  STATUS: PATTERN

Step 3: sin(theta_13) = phi^(-23/2)
  WHERE: 23 is an H4' complement Coxeter exponent
  STATUS: PATTERN

Alternative step 2-3 (integer version):
  sin(theta_23) = phi^(-7): 7 = |1-29|/4 from H4 Coxeter pair differences
  sin(theta_13) = phi^(-12): 12 = (19+29)/4 = degree of 600-cell
""")

# Build both versions and compare
print("--- COMPARISON: Half-integer vs Integer CKM ---\n")

# Half-integer: (3, 13/2, 23/2)
V_half = build_ckm(PHI**(-3), PHI**(-13/2), PHI**(-23/2))
mH, xH, cH = compute_errors(V_half, "HALF-INTEGER: sin=phi^(-n), n=(3, 13/2, 23/2)")

# Integer: (3, 7, 12)
V_int = build_ckm(PHI**(-3), PHI**(-7), PHI**(-12))
mI, xI, cI = compute_errors(V_int, "INTEGER: sin=phi^(-n), n=(3, 7, 12)")

# With corrections: (3, 13/2, 23/2) + perturbative corrections
# The corrections needed are:
corr_12 = n_exp_12 - 3
corr_23 = n_exp_23 - 13/2
corr_13 = n_exp_13 - 23/2
print(f"\n--- CORRECTIONS needed for (3, 13/2, 23/2) to match experiment ---")
print(f"  delta_12 = {corr_12:+.4f}")
print(f"  delta_23 = {corr_23:+.4f}")
print(f"  delta_13 = {corr_13:+.4f}")
print(f"  Sum: {corr_12 + corr_23 + corr_13:+.4f}")
print(f"  delta_13 ~ delta_12 + delta_23: {corr_12+corr_23:+.4f} vs {corr_13:+.4f}")
print(f"  Difference: {abs(corr_13 - (corr_12+corr_23)):.4f}")

# Check if corrections relate to physics constants
print(f"\n  Physical interpretation of corrections:")
print(f"  delta_12 = {corr_12:+.4f} ~ sin^2(tW)/2 = {SIN2_TW/2:.4f} ({abs(corr_12-SIN2_TW/2)/abs(corr_12)*100:.1f}%)")
print(f"  delta_23 = {corr_23:+.4f} ~ -alpha_s/pi = {-ALPHA_S/np.pi:.4f} ({abs(corr_23-(-ALPHA_S/np.pi))/abs(corr_23)*100:.1f}%)")
print(f"  delta_23 = {corr_23:+.4f} ~ -1/(2*phi^3) = {-1/(2*PHI**3):.4f} ({abs(corr_23-(-1/(2*PHI**3)))/abs(corr_23)*100:.1f}%)")

# Try corrected version with physical corrections
n12_corr = 3 + SIN2_TW / 2
n23_corr = 13/2 - ALPHA_S / np.pi
n13_corr = 23/2 + (SIN2_TW/2 - ALPHA_S/np.pi)

print(f"\n--- CORRECTED CKM with physics perturbations ---")
print(f"  n_12 = 3 + sin^2(tW)/2 = {n12_corr:.4f} (exact: {n_exp_12:.4f})")
print(f"  n_23 = 13/2 - alpha_s/pi = {n23_corr:.4f} (exact: {n_exp_23:.4f})")
print(f"  n_13 = 23/2 + (sin^2(tW)/2 - alpha_s/pi) = {n13_corr:.4f} (exact: {n_exp_13:.4f})")

V_corr = build_ckm(PHI**(-n12_corr), PHI**(-n23_corr), PHI**(-n13_corr))
mC, xC, cC = compute_errors(V_corr, "CORRECTED: n + physics perturbations")

# Alternative correction: purely geometric
print(f"\n--- GEOMETRIC CORRECTIONS (no physics constants) ---")
# delta_12 ~ 1/(2*phi^3) = 0.118? No, that's alpha_s. Not right.
# delta_12 ~ 1/10 = 0.1? Close to 0.10 actual
# delta_23 ~ -1/8 = -0.125? Close to -0.12 actual
# delta_13 ~ -1/14 = -0.071? Actual is -0.07

# From 600-cell:
# 10 = vertices per decagon
# 8 = |n_u - n_s|
# 14 = mid-shell size? Not clear.

# Try: corrections from intersection numbers
# delta ~ a_1/(a_1*b_1*something)
print(f"  delta_12 = {corr_12:+.4f}: ")
print(f"    1/10 = {1/10:.4f} ({abs(corr_12-1/10)/abs(corr_12)*100:.1f}%)")
print(f"    1/(2*5) = {1/10:.4f}")
print(f"    a_1/(a_1^2+b_1) = {5/(25+6):.4f} ({abs(corr_12-5/31)/abs(corr_12)*100:.1f}%)")
print(f"  delta_23 = {corr_23:+.4f}: ")
print(f"    -1/(2*b_1) = {-1/12:.4f} ({abs(corr_23-(-1/12))/abs(corr_23)*100:.1f}%)")
print(f"    -1/(a_1+b_1-1) = {-1/10:.4f} ({abs(corr_23-(-0.1))/abs(corr_23)*100:.1f}%)")
print(f"    -alpha_s = {-ALPHA_S:.4f} ({abs(corr_23-(-ALPHA_S))/abs(corr_23)*100:.1f}%)")
print(f"  delta_13 = {corr_13:+.4f}: ")
print(f"    delta_12 + delta_23 = {corr_12+corr_23:.4f}")


# ================================================================
# NEW: APPROACH 11: DIRECT COMPUTATION WITH (a,b) NORM
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 11: CKM FROM (a,b) NORM IN Z[phi] LATTICE")
print("=" * 80)

print("""
From exp121: the Z[phi] norm |N| = |a^2 + ab - b^2| plays a crucial role
in level selection. The occupied |N| values are {0, 1, 5, 19}.

Idea: the CKM mixing might be related to how the norms of the up-type
and down-type (a,b) vectors relate through the Z[phi] lattice.
""")

# Compute Z[phi] norms
def zphi_norm(a, b):
    """Compute |N(a + b*phi)| = |a^2 + ab - b^2|"""
    return abs(a**2 + a*b - b**2)

print(f"  Z[phi] norms of quark quantum numbers:")
for i, q in enumerate(labels_up):
    N = zphi_norm(a_up[i], b_up[i])
    print(f"    {q}({a_up[i]:+d},{b_up[i]:+d}): |N| = {N}")
for i, q in enumerate(labels_dn):
    N = zphi_norm(a_dn[i], b_dn[i])
    print(f"    {q}({a_dn[i]:+d},{b_dn[i]:+d}): |N| = {N}")

# Norm difference matrix
print(f"\n  |N_up - N_dn| matrix:")
norm_diff = np.zeros((3, 3), dtype=int)
for i in range(3):
    for j in range(3):
        norm_diff[i, j] = abs(zphi_norm(a_up[i], b_up[i]) - zphi_norm(a_dn[j], b_dn[j]))

print(f"        d          s          b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {norm_diff[i,j]:8d}"
    print(row)

# Product of norms
print(f"\n  Product N_up * N_dn matrix:")
norm_prod = np.zeros((3, 3), dtype=int)
for i in range(3):
    for j in range(3):
        norm_prod[i, j] = zphi_norm(a_up[i], b_up[i]) * zphi_norm(a_dn[j], b_dn[j])

print(f"        d          s          b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {norm_prod[i,j]:8d}"
    print(row)

# Try CKM from norm products
# V_ij ~ phi^(-alpha * N_up * N_dn)
print(f"\n  CKM from phi^(-alpha * N_i * N_j) scan:")
best_norm_err = 1e10
best_norm_alpha = 0
for alpha_val in np.linspace(0.001, 1.0, 10000):
    V_pred = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            V_pred[i, j] = PHI**(-alpha_val * norm_prod[i, j])
        V_pred[i] /= np.linalg.norm(V_pred[i])
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_norm_err:
        best_norm_err = err
        best_norm_alpha = alpha_val

print(f"  Best alpha = {best_norm_alpha:.4f}, mean err = {best_norm_err:.2f}%")
print(f"  NOTE: N_u=19, N_d=1, N_s=1, N_b=19 -> norm product is DEGENERATE")
print(f"  This means N_ud = N_us = 19 (identical!) which CANNOT reproduce CKM hierarchy.")
print(f"  CONCLUSION: Z[phi] norm alone does not determine CKM. NEGATIVE RESULT.")


# ================================================================
# APPROACH 12: THE 3-LINE STRUCTURE
# ================================================================
print("\n" + "=" * 80)
print("APPROACH 12: CKM FROM 3-LINE STRUCTURE IN (a,b) LATTICE")
print("=" * 80)

print("""
From exp121c: all 9 fermion (a,b) lie on 3 lines:
  L1 (a=1): d(1,0), mu/s(1,1), tau(1,2)
  L2 (b=1): mu/s(1,1), c(2,1), t(4,1)
  L3 (3a+2b=5): u(3,-2), mu/s(1,1), b(-1,4)
All 3 lines intersect at mu/s(1,1) = TRIPLE POINT.

The CKM matrix should be determined by how up-type and down-type quarks
are distributed across these lines.
""")

# Line assignments
def get_lines(a, b):
    lines = []
    if a == 1: lines.append('L1')
    if b == 1: lines.append('L2')
    if 3*a + 2*b == 5: lines.append('L3')
    return lines

print(f"  Line assignments:")
for i, q in enumerate(labels_up):
    lines = get_lines(a_up[i], b_up[i])
    print(f"    {q}({a_up[i]:+d},{b_up[i]:+d}): {lines}")
for i, q in enumerate(labels_dn):
    lines = get_lines(a_dn[i], b_dn[i])
    print(f"    {q}({a_dn[i]:+d},{b_dn[i]:+d}): {lines}")

# Shared lines matrix
shared_lines = np.zeros((3, 3), dtype=int)
for i in range(3):
    for j in range(3):
        lu = set(get_lines(a_up[i], b_up[i]))
        ld = set(get_lines(a_dn[j], b_dn[j]))
        shared_lines[i, j] = len(lu & ld)

print(f"\n  Shared lines matrix:")
print(f"        d    s    b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {shared_lines[i,j]:3d}"
    print(row)

# Graph distance along lines
# On L2 (b=1): s(1,1), c(2,1), t(4,1) - distance in a-coordinate
# On L3 (3a+2b=5): u(3,-2), s(1,1), b(-1,4) - parametric distance
print(f"\n  Along-line distances:")
print(f"  On L2 (b=1): s(a=1) -- c(a=2) -- t(a=4)")
print(f"    d(s,c) = 1, d(c,t) = 2, d(s,t) = 3")
print(f"  On L3 (3a+2b=5): u(3,-2) -- s(1,1) -- b(-1,4)")
print(f"    d(u,s) = sqrt((3-1)^2+(-2-1)^2) = sqrt(13) = {np.sqrt(13):.4f}")
print(f"    d(s,b) = sqrt((1+1)^2+(1-4)^2) = sqrt(13) = {np.sqrt(13):.4f}")
print(f"    NOTE: d(u,s) = d(s,b) exactly! s is the MIDPOINT on L3!")
print(f"  On L1 (a=1): d(1,0) -- s(1,1) -- tau(1,2)")
print(f"    d(d,s) = 1, d(s,tau) = 1, d(d,tau) = 2")

# Key insight: the CKM MIXING depends on quarks being on DIFFERENT lines!
# u is on L3 only
# c is on L2 only
# t is on L2 only
# d is on L1 only
# s is on L1, L2, L3 (triple point)
# b is on L3 only

# So the "direct" line connections are:
# u-d: no shared line -> INDIRECT (through s)
# u-s: shared L3 -> DIRECT
# u-b: shared L3 -> DIRECT
# c-d: no shared line -> INDIRECT
# c-s: shared L2 -> DIRECT
# c-b: no shared line -> INDIRECT
# t-d: no shared line -> INDIRECT
# t-s: shared L2 -> DIRECT
# t-b: no shared line -> INDIRECT

print(f"\n  Direct/Indirect connections (shared lines > 0 = direct):")
print(f"        d    s    b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        tag = "D" if shared_lines[i, j] > 0 else "I"
        row += f"    {tag}"
    print(row)

print(f"\n  Observation: V_us (DIRECT on L3), V_cb (INDIRECT) have similar magnitudes")
print(f"  But V_us = 0.224 (large) while V_cb = 0.041 (smaller)")
print(f"  Direct doesn't mean large mixing -- the along-line DISTANCE matters!")

# Compute effective distance: on same line = along-line dist, different = shortest path through s
print(f"\n--- 12a: Effective graph distance through line network ---")

# Build mini-graph: nodes = quarks, edges = shared lines with distances
# s is the hub connecting everything

# For CKM: path u -> d
# u(L3) -> s(L1,L2,L3) -> d(L1): 2 line hops, total Euclidean = sqrt(13)+1
# u(L3) -> b(L3) -> s(L1,L2,L3) -> d(L1): longer path

# For CKM: path c -> s
# c(L2) -> s(L2): 1 line hop, dist = 1

# For CKM: path c -> b
# c(L2) -> s(L1,L2,L3) -> b(L3): 2 line hops, dist = 1 + sqrt(13)

# So the effective distances are:
eff_dists = {
    ('u','d'): np.sqrt(13) + 1,       # u->s(L3)->d(L1)
    ('u','s'): np.sqrt(13),            # direct on L3
    ('u','b'): 2*np.sqrt(13),          # direct on L3
    ('c','d'): 1 + 1,                  # c->s(L2)->d(L1)
    ('c','s'): 1,                      # direct on L2
    ('c','b'): 1 + np.sqrt(13),        # c->s(L2,L3)->b(L3)
    ('t','d'): 3 + 1,                  # t->s(L2)->d(L1)
    ('t','s'): 3,                      # direct on L2
    ('t','b'): 3 + np.sqrt(13),        # t->s(L2,L3)->b(L3)
}

eff_dist_matrix = np.zeros((3, 3))
for i, u in enumerate(labels_up):
    for j, d in enumerate(labels_dn):
        eff_dist_matrix[i, j] = eff_dists[(u, d)]

print(f"  Effective line-network distance matrix:")
print(f"        d          s          b")
for i in range(3):
    row = f"  {labels_up[i]}  "
    for j in range(3):
        row += f"  {eff_dist_matrix[i,j]:8.4f}"
    print(row)

# Fit: V_ij ~ phi^(-alpha * d_eff)
best_12_err = 1e10
best_12_alpha = 0
best_12_V = None
for alpha_val in np.linspace(0.01, 5.0, 5000):
    V_pred = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            V_pred[i, j] = PHI**(-alpha_val * eff_dist_matrix[i, j])
        V_pred[i] /= np.linalg.norm(V_pred[i])
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_12_err:
        best_12_err = err
        best_12_alpha = alpha_val
        best_12_V = V_pred.copy()

print(f"\n  Best: alpha = {best_12_alpha:.4f}, mean err = {best_12_err:.2f}%")
m12, x12, c12 = compute_errors(best_12_V, "12a: phi^(-alpha * d_line_network)")


# ================================================================
# FINAL SUMMARY
# ================================================================
print("\n" + "=" * 80)
print("FINAL SUMMARY: ALL APPROACHES COMPARED")
print("=" * 80)

summary = [
    ("Half-int (3, 13/2, 23/2)", mH, xH, 0, "PATTERN"),
    ("Integer (3, 7, 12)", mI, xI, 0, "PATTERN"),
    ("Corrected (physics perturb.)", mC, xC, 2, "PATTERN+FITTING"),
    ("Wolfenstein lam=phi^-3, A=phi^(-1/phi)", m9c, x9c, 1, "PATTERN+FITTING"),
    ("Wolfenstein lam=phi^-3, A=phi^(-1/2)", m9d, x9d, 1, "PATTERN+FITTING"),
    ("phi^(-alpha|da|+beta|db|)", best_5c_err, 0, 2, "FITTING"),
    ("Line network distance", best_12_err, 0, 1, "GEOMETRIC+FITTING"),
]

# Sort by mean error
summary.sort(key=lambda x: x[1])

print(f"\n{'Approach':<50s} {'Mean%':>7s} {'Free':>5s} {'Status':<20s}")
print("-" * 85)
for name, mean_err, max_err, n_free, cat in summary:
    print(f"  {name:<48s} {mean_err:7.2f} {n_free:5d}  {cat}")


# ================================================================
# CONCLUSIONS
# ================================================================
print("\n" + "=" * 80)
print("CONCLUSIONS")
print("=" * 80)

print(f"""
=== WHAT WE FOUND (NEW IN EXP-129) ===

1. APPROACH 2 (Level differences): Level differences |n_up - n_dn| CANNOT
   give a uniform CKM formula. The required divisors are wildly different
   (2.58, 0.45, 1.38). NEGATIVE RESULT.

2. APPROACH 5 (Intersection numbers): 5|da|+6|db| is also non-uniform.
   But phi^(-alpha*|da| - beta*|db|) with TWO free parameters gives
   a reasonable fit. The fitted ratio alpha/beta = {best_5c_params[0]/best_5c_params[1]:.4f}
   is NOT close to 5/6 = {5/6:.4f}. FITTING, not derivation.

3. APPROACH 8 (Eigenvector overlap): The 9-dimensional eigenspaces (mult=3^2)
   give canonical angles between the two such spaces, but these DO NOT match
   CKM angles. The structure is there but the connection is not clean.
   INCONCLUSIVE.

4. APPROACH 9 (Wolfenstein): The Wolfenstein A parameter satisfies
   A ~ phi^(-1/phi) = {PHI**(-1/PHI):.4f} with {abs(PHI**(-1/PHI)-A_phi)/A_phi*100:.2f}% error.
   Self-referential formula: SPECULATIVE but intriguing.

5. APPROACH 11 (Z[phi] norm): The Z[phi] norm is DEGENERATE for CKM
   (N_u = N_b = 19, N_d = N_s = 1). Cannot distinguish mixing. NEGATIVE.

6. APPROACH 12 (3-line structure): Line-network distances through the
   triple point s(1,1) provide a GEOMETRIC explanation of CKM hierarchy.
   With 1 free parameter: mean error {best_12_err:.1f}%.
   The structure is: s is the HUB connecting all sectors via 3 lines.

=== ESTABLISHED RESULTS (confirmed) ===

7. BEST ZERO-PARAMETER CKM: n = (3, 13/2, 23/2)
   Mean error: {mH:.2f}%. Exponents from H4 + H4' Coxeter.
   STATUS: PATTERN

8. CABIBBO ANGLE: sin(theta_C) = phi^(-3), 3 = (1+11)/4.
   This remains the STRONGEST single prediction.
   STATUS: PATTERN (the /4 is the rank of H4, natural but not derived)

9. CORRECTIONS: Adding sin^2(tW)/2 and alpha_s/pi as perturbative
   corrections to the Coxeter exponents improves the match significantly.
   Mean error: {mC:.2f}% with 2 parameters (physics constants, not free).
   STATUS: PATTERN + FITTING

=== HONEST ASSESSMENT ===

The CKM matrix in the 600-cell framework is PATTERN, not DERIVAT.

The framework provides:
(a) The HIERARCHY (phi-power suppression) -- natural from geometry
(b) The EXPONENTS (from Coxeter numbers) -- pattern matching
(c) The LINE STRUCTURE (3-line network with s as hub) -- geometric

What is MISSING for full derivation:
(a) WHY divide Coxeter exponents by 4 (or 2)?
(b) WHY these specific Coxeter exponent PAIRINGS?
(c) A first-principles mechanism generating CKM from the 600-cell
    (e.g., from the graph Laplacian, or from eigenvector overlaps)

The most promising path forward: understand the CKM as arising from
the REPRESENTATION THEORY of H4 (irrep decomposition of the 9-dim
eigenspaces), not from numerological coincidences with Coxeter numbers.
""")

"""
EXP-314: Non-Abelian S^3 Holonomy of Triangular Faces
======================================================
GOAL: Compute SO(3) holonomy for each of 1200 triangular faces,
decompose via Hopf fibration, group by gauge sector (ACC/BCC/CCC).

KEY PREDICTION: Even though edges are isotropic (exp313/313c),
the triangle NORMAL direction varies, and its Hopf fiber component
correlates with gauge sector -> sector-dependent holonomy.

NOTE: All print uses ASCII only (Windows cp1252).
"""

import numpy as np
import math
from itertools import product as iterproduct
from collections import Counter, defaultdict

PHI = (1 + math.sqrt(5)) / 2
PHI_P = (1 - math.sqrt(5)) / 2
PI = math.pi
a1 = 5
b1 = 6

# Framework constants
ALPHA_S_FW = 1 / (2 * PHI**3)
SIN2_TW_FW = b1 / (a1**2 + 1)
_disc = (4*a1*PHI**4)**2 - 8*PI
ALPHA_FW = (4*a1*PHI**4 - math.sqrt(_disc)) / (4*PI)

print("=" * 72)
print("EXP-314: NON-ABELIAN S^3 HOLONOMY OF TRIANGULAR FACES")
print("=" * 72)

# =====================================================================
# SECTION 1: Build 600-cell
# =====================================================================
print("\n--- SECTION 1: Build 600-cell ---")

def qmul(a, b):
    return np.array([
        a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
        a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0]
    ])

def generate_2I():
    elements = set()
    def add(q):
        q = tuple(round(x, 10) for x in q)
        elements.add(q)
        elements.add(tuple(-x for x in q))
    for i in range(4):
        v = [0,0,0,0]; v[i] = 1.0; add(v)
    for s0 in [0.5,-0.5]:
        for s1 in [0.5,-0.5]:
            for s2 in [0.5,-0.5]:
                for s3 in [0.5,-0.5]:
                    add([s0,s1,s2,s3])
    abs_vals = [0.0, 0.5, abs(PHI_P)/2, PHI/2]
    even_perms = [
        (0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
        (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0),
    ]
    for perm in even_perms:
        base = [abs_vals[perm[k]] for k in range(4)]
        non_zero = [i for i in range(4) if abs(base[i]) > 1e-10]
        for signs in iterproduct([1,-1], repeat=len(non_zero)):
            v = list(base)
            for idx, s in zip(non_zero, signs):
                v[idx] *= s
            add(v)
    return [np.array(q) for q in sorted(elements)]

elements_2I = generate_2I()
V = np.array(elements_2I)
n_vert = len(V)
print(f"  |2I| = {n_vert}")

IP = V @ V.T
A_adj = np.zeros((n_vert, n_vert), dtype=int)
for i in range(n_vert):
    for j in range(i+1, n_vert):
        if abs(IP[i,j] - PHI/2) < 0.01:
            A_adj[i,j] = A_adj[j,i] = 1
n_edges = np.sum(A_adj) // 2
print(f"  Edges = {n_edges}")

# =====================================================================
# SECTION 2: Vertex classification (A=2, B=22, C=96)
# =====================================================================
print("\n--- SECTION 2: Vertex classification ---")

def classify_vertex(q):
    abs_coords = [abs(x) for x in q]
    n_zero = sum(1 for x in abs_coords if x < 0.01)
    if n_zero == 3 and abs(abs(q[0]) - 1.0) < 0.01:
        return 'A'
    n_half = sum(1 for x in abs_coords if abs(x - 0.5) < 0.01)
    if n_half == 4:
        return 'B'
    if n_zero == 3:
        return 'B'  # other axis vertices
    return 'C'

vtypes = [classify_vertex(V[i]) for i in range(n_vert)]
type_counts = Counter(vtypes)
print(f"  A={type_counts['A']}, B={type_counts['B']}, C={type_counts['C']}")
assert type_counts['A'] == 2 and type_counts['B'] == 22 and type_counts['C'] == 96

# =====================================================================
# SECTION 3: Find all 1200 triangular faces
# =====================================================================
print("\n--- SECTION 3: Find triangular faces ---")

faces = []
for i in range(n_vert):
    ni = set(j for j in range(n_vert) if A_adj[i,j])
    for j in ni:
        if j > i:
            nj = set(k for k in range(n_vert) if A_adj[j,k])
            for k in (ni & nj):
                if k > j:
                    faces.append((i, j, k))

print(f"  Found {len(faces)} triangular faces")

# Classify by vertex type
face_vtypes = []
for i, j, k in faces:
    ft = tuple(sorted([vtypes[i], vtypes[j], vtypes[k]]))
    face_vtypes.append(ft)

face_type_counts = Counter(face_vtypes)
print(f"\n  Face types:")
for ft, cnt in sorted(face_type_counts.items()):
    pct = 100 * cnt / len(faces)
    print(f"    {ft}: {cnt} ({pct:.1f}%)")

# Gauge sector mapping
def face_gauge_sector(ft):
    if 'A' in ft:
        return 'U1'
    elif 'B' in ft:
        return 'SU2'
    else:
        return 'SU3'

face_sectors = [face_gauge_sector(ft) for ft in face_vtypes]
sector_counts = Counter(face_sectors)
print(f"\n  Gauge sectors:")
for s in ['U1', 'SU2', 'SU3']:
    print(f"    {s}: {sector_counts[s]}")

# =====================================================================
# SECTION 4: Triangle normals via 4D cross product
# =====================================================================
print("\n--- SECTION 4: Triangle normals and fiber components ---")

def triple_cross_4d(p, u, w):
    """4D cross product: n = *(p ^ u ^ w), perpendicular to all three."""
    M = np.array([p, u, w])  # 3x4
    n = np.zeros(4)
    for i in range(4):
        cols = [j for j in range(4) if j != i]
        n[i] = ((-1)**i) * np.linalg.det(M[:, cols])
    return n

def hopf_fiber_tangent(q):
    """Fiber direction T_i(q) = q*i at point q on S^3."""
    return np.array([-q[1], q[0], -q[3], q[2]])

# For each triangle: compute normal, fiber component
normals = []          # unit normal in R^4 (perp to p and triangle plane)
fiber_components = [] # |n . T_i(p)| / |n|
signed_fiber = []     # n . T_i(p) / |n| (signed)

for fi, (i, j, k) in enumerate(faces):
    p, q, r = V[i], V[j], V[k]

    # Tangent vectors at p
    u = q - np.dot(q, p) * p  # toward q
    w = r - np.dot(r, p) * p  # toward r

    # Normal to triangle plane in T_p S^3
    n = triple_cross_4d(p, u, w)
    n_norm = np.linalg.norm(n)

    if n_norm < 1e-10:
        normals.append(np.zeros(4))
        fiber_components.append(0.0)
        signed_fiber.append(0.0)
        continue

    n_hat = n / n_norm
    normals.append(n_hat)

    # Fiber direction at p
    T_i = hopf_fiber_tangent(p)
    dot = np.dot(n_hat, T_i)
    fiber_components.append(abs(dot))
    signed_fiber.append(dot)

fiber_components = np.array(fiber_components)
signed_fiber = np.array(signed_fiber)

print(f"  Fiber component stats (all 1200 triangles):")
print(f"    |n.T_i|: mean={fiber_components.mean():.6f}, "
      f"std={fiber_components.std():.6f}")
print(f"    Range: [{fiber_components.min():.6f}, {fiber_components.max():.6f}]")

# Distribution of fiber components
fc_rounded = np.round(fiber_components, 4)
distinct_fc = sorted(set(fc_rounded))
print(f"    Distinct values: {len(distinct_fc)}")
for val in distinct_fc[:10]:
    cnt = np.sum(np.abs(fc_rounded - val) < 0.001)
    print(f"      |n.T_i| = {val:.4f}: {cnt} triangles")

# =====================================================================
# SECTION 5: Fiber component per gauge sector
# =====================================================================
print("\n--- SECTION 5: Fiber component per gauge sector ---")

print(f"  {'Sector':>8s} {'<|n.T|>':>10s} {'std':>10s} {'<n.T>':>10s} {'n':>6s}")
print(f"  {'-'*50}")

sector_fiber_abs = defaultdict(list)
sector_fiber_signed = defaultdict(list)

for fi in range(len(faces)):
    gs = face_sectors[fi]
    sector_fiber_abs[gs].append(fiber_components[fi])
    sector_fiber_signed[gs].append(signed_fiber[fi])

sector_means = {}
for gs in ['U1', 'SU2', 'SU3']:
    vals_abs = sector_fiber_abs[gs]
    vals_sig = sector_fiber_signed[gs]
    m_abs = np.mean(vals_abs)
    s_abs = np.std(vals_abs)
    m_sig = np.mean(vals_sig)
    sector_means[gs] = m_abs
    print(f"  {gs:>8s} {m_abs:10.6f} {s_abs:10.6f} {m_sig:10.6f} {len(vals_abs):6d}")

# Anisotropy test
print(f"\n  ANISOTROPY CHECK:")
print(f"    <|n.T|>_U1  = {sector_means['U1']:.6f}")
print(f"    <|n.T|>_SU2 = {sector_means['SU2']:.6f}")
print(f"    <|n.T|>_SU3 = {sector_means['SU3']:.6f}")

d_12 = sector_means['U1'] - sector_means['SU2']
d_13 = sector_means['U1'] - sector_means['SU3']
d_23 = sector_means['SU2'] - sector_means['SU3']
print(f"    U1 - SU2 = {d_12:.6f}")
print(f"    U1 - SU3 = {d_13:.6f}")
print(f"    SU2 - SU3 = {d_23:.6f}")

# Also by face type (finer)
print(f"\n  By face type (finer):")
print(f"  {'Type':>15s} {'<|n.T|>':>10s} {'std':>10s} {'n':>6s}")
for ft in sorted(face_type_counts.keys()):
    mask = [fi for fi in range(len(faces)) if face_vtypes[fi] == ft]
    vals = fiber_components[mask]
    print(f"  {str(ft):>15s} {vals.mean():10.6f} {vals.std():10.6f} {len(vals):6d}")

# =====================================================================
# SECTION 6: Full SO(3) holonomy verification
# =====================================================================
print("\n--- SECTION 6: Full SO(3) holonomy (verification) ---")

def parallel_transport_S3(v, p, q):
    """Parallel transport v from T_p to T_q on unit S^3."""
    c = np.dot(p, q)
    if abs(1 + c) < 1e-10:
        return -v
    return v - (np.dot(v, q) / (1 + c)) * (p + q)

def quat_basis(p):
    """Quaternionic ON basis {T_i, T_j, T_k} for T_p S^3."""
    Ti = np.array([-p[1], p[0], -p[3], p[2]])
    Tj = np.array([-p[2], p[3], p[0], -p[1]])
    Tk = np.array([-p[3], -p[2], p[1], p[0]])
    return np.array([Ti, Tj, Tk])

# Compute full holonomy for ALL 1200 triangles
hol_angles = np.zeros(len(faces))
hol_fiber_comp = np.zeros(len(faces))  # fiber component of rotation axis
hol_matrices = []

for fi, (i, j, k) in enumerate(faces):
    p, q, r = V[i], V[j], V[k]
    basis = quat_basis(p)  # 3x4

    H = np.zeros((3, 3))
    for col in range(3):
        v = basis[col].copy()
        v = parallel_transport_S3(v, p, q)
        v = parallel_transport_S3(v, q, r)
        v = parallel_transport_S3(v, r, p)
        for row in range(3):
            H[row, col] = np.dot(v, basis[row])

    # Rotation angle
    tr = np.trace(H)
    cos_th = np.clip((tr - 1) / 2, -1, 1)
    theta = np.arccos(cos_th)
    hol_angles[fi] = theta

    # Rotation axis from antisymmetric part
    axis = np.array([H[2,1]-H[1,2], H[0,2]-H[2,0], H[1,0]-H[0,1]])
    axis_norm = np.linalg.norm(axis)
    if axis_norm > 1e-10:
        axis /= axis_norm
    # axis[0] = fiber component (T_i direction in quaternionic basis)
    hol_fiber_comp[fi] = abs(axis[0])
    hol_matrices.append(H)

print(f"  Holonomy angle stats:")
print(f"    mean = {hol_angles.mean():.6f} rad = {math.degrees(hol_angles.mean()):.4f} deg")
print(f"    std  = {hol_angles.std():.2e}")
print(f"    All equal? {hol_angles.std() < 1e-8}")

# Expected: spherical excess of equilateral triangle with side arccos(phi/2)
cos_a = PHI / 2
angle_A = math.acos(cos_a / (1 + cos_a))
expected_excess = 3 * angle_A - PI
print(f"    Expected (spherical excess) = {expected_excess:.6f} rad = {math.degrees(expected_excess):.4f} deg")
print(f"    Match: {abs(hol_angles.mean() - expected_excess) < 1e-6}")

# Verify SO(3)
det_errs = []
orth_errs = []
for H in hol_matrices[:100]:
    det_errs.append(abs(np.linalg.det(H) - 1.0))
    orth_errs.append(np.max(np.abs(H @ H.T - np.eye(3))))
print(f"\n  SO(3) check (first 100):")
print(f"    max |det-1| = {max(det_errs):.2e}")
print(f"    max |HH^T-I| = {max(orth_errs):.2e}")

# Compare axis-based fiber component with normal-based
print(f"\n  Cross-check: axis fiber comp vs normal fiber comp:")
corr = np.corrcoef(hol_fiber_comp, fiber_components)[0, 1]
print(f"    Correlation = {corr:.6f}")
diff = np.abs(hol_fiber_comp - fiber_components)
print(f"    Max difference = {diff.max():.2e}")
print(f"    Mean difference = {diff.mean():.2e}")

# =====================================================================
# SECTION 7: Holonomy decomposition per gauge sector
# =====================================================================
print("\n--- SECTION 7: Holonomy per gauge sector (from full SO(3)) ---")

# The holonomy H in quaternionic basis {T_i, T_j, T_k}:
# H[0,0] = fiber-fiber, H[1:,1:] = horizontal block
# Berry phase = rotation angle in horizontal block

berry_phases = np.zeros(len(faces))
horiz_angles = np.zeros(len(faces))

for fi in range(len(faces)):
    H = hol_matrices[fi]
    # Horizontal 2x2 block
    H_horiz = H[1:3, 1:3]
    # Berry phase = angle of this 2x2 rotation
    cos_berry = np.clip((H_horiz[0,0] + H_horiz[1,1]) / 2, -1, 1)
    berry_phases[fi] = np.arccos(cos_berry)

    # Fiber-fiber component
    # H[0,0] = cos of the angle that the fiber rotates
    # If holonomy is purely around fiber axis: H = diag(1, R_2x2)
    # The fiber component of the rotation:
    horiz_angles[fi] = hol_angles[fi] * np.sqrt(max(0, 1 - hol_fiber_comp[fi]**2))

print(f"  Berry phase stats:")
print(f"    mean = {berry_phases.mean():.6f}, std = {berry_phases.std():.6f}")

print(f"\n  Per gauge sector:")
print(f"  {'Sector':>8s} {'<theta>':>10s} {'<|axis_f|>':>12s} {'<Berry>':>10s} "
      f"{'<H[0,0]>':>10s} {'n':>6s}")
print(f"  {'-'*65}")

sector_results = {}
for gs in ['U1', 'SU2', 'SU3']:
    mask = [fi for fi in range(len(faces)) if face_sectors[fi] == gs]
    th = hol_angles[mask]
    af = hol_fiber_comp[mask]
    bp = berry_phases[mask]
    h00 = np.array([hol_matrices[fi][0,0] for fi in mask])
    sector_results[gs] = {
        'theta': th.mean(), 'axis_fiber': af.mean(),
        'berry': bp.mean(), 'H00': h00.mean(), 'n': len(mask)
    }
    print(f"  {gs:>8s} {th.mean():10.6f} {af.mean():12.6f} {bp.mean():10.6f} "
          f"{h00.mean():10.6f} {len(mask):6d}")

# =====================================================================
# SECTION 8: Detailed matrix analysis per sector
# =====================================================================
print("\n--- SECTION 8: Mean holonomy matrix per sector ---")

for gs in ['U1', 'SU2', 'SU3']:
    mask = [fi for fi in range(len(faces)) if face_sectors[fi] == gs]
    H_mean = np.mean([hol_matrices[fi] for fi in mask], axis=0)
    print(f"\n  <H>_{gs}:")
    for row in range(3):
        print(f"    [{H_mean[row,0]:8.5f} {H_mean[row,1]:8.5f} {H_mean[row,2]:8.5f}]")
    # Eigenvalues of mean matrix
    evals = np.linalg.eigvals(H_mean)
    print(f"    eigenvalues: {[f'{e.real:.5f}+{e.imag:.5f}i' for e in evals]}")

# Also compute variance of each matrix element per sector
print(f"\n  Variance of H matrix elements:")
print(f"  {'Sector':>8s} {'Var(H00)':>12s} {'Var(H01)':>12s} {'Var(H11)':>12s}")
print(f"  {'-'*50}")
for gs in ['U1', 'SU2', 'SU3']:
    mask = [fi for fi in range(len(faces)) if face_sectors[fi] == gs]
    h00 = [hol_matrices[fi][0,0] for fi in mask]
    h01 = [hol_matrices[fi][0,1] for fi in mask]
    h11 = [hol_matrices[fi][1,1] for fi in mask]
    print(f"  {gs:>8s} {np.var(h00):12.6f} {np.var(h01):12.6f} {np.var(h11):12.6f}")

# =====================================================================
# SECTION 9: Compare with framework constants
# =====================================================================
print("\n--- SECTION 9: Comparison with corrections ---")

print(f"  Framework corrections:")
print(f"    sin2tW = {SIN2_TW_FW:.6f}")
print(f"    alpha_s = {ALPHA_S_FW:.6f}")
print(f"    1/a1 = {1/a1:.6f}")
print(f"    1/phi^4 = {1/PHI**4:.6f}")
print(f"    alpha_s*2/pi = {ALPHA_S_FW*2/PI:.6f}")

# Holonomy angle
Omega = hol_angles.mean()  # same for all
print(f"\n  Holonomy angle Omega = {Omega:.6f}")

# Check ratios of fiber components between sectors
f_U1 = sector_results['U1']['axis_fiber']
f_SU2 = sector_results['SU2']['axis_fiber']
f_SU3 = sector_results['SU3']['axis_fiber']

print(f"\n  Fiber component ratios:")
if f_SU3 > 1e-6:
    print(f"    f_U1/f_SU3 = {f_U1/f_SU3:.6f}")
    print(f"    f_SU2/f_SU3 = {f_SU2/f_SU3:.6f}")
if f_SU2 > 1e-6:
    print(f"    f_U1/f_SU2 = {f_U1/f_SU2:.6f}")

# Berry phase per sector
b_U1 = sector_results['U1']['berry']
b_SU2 = sector_results['SU2']['berry']
b_SU3 = sector_results['SU3']['berry']

print(f"\n  Berry phase per sector:")
print(f"    U1:  {b_U1:.6f}")
print(f"    SU2: {b_SU2:.6f}")
print(f"    SU3: {b_SU3:.6f}")

# Differences
print(f"\n  Berry phase differences:")
print(f"    U1 - SU2 = {b_U1 - b_SU2:.6f}")
print(f"    U1 - SU3 = {b_U1 - b_SU3:.6f}")
print(f"    SU2 - SU3 = {b_SU2 - b_SU3:.6f}")

# Check: does Omega * (f_sector - f_overall) give corrections?
f_all = fiber_components.mean()
print(f"\n  Overall <|n.T|> = {f_all:.6f}")
for gs in ['U1', 'SU2', 'SU3']:
    delta_f = sector_means[gs] - f_all
    correction = Omega * delta_f
    print(f"    Omega * (f_{gs} - f_avg) = {correction:.6f}")

# =====================================================================
# SECTION 10: Fiber component distribution per sector
# =====================================================================
print("\n--- SECTION 10: Fiber component distributions ---")

# Histogram of fiber components per sector
for gs in ['U1', 'SU2', 'SU3']:
    vals = sector_fiber_abs[gs]
    n_bins = 5
    bins = np.linspace(0, 1, n_bins + 1)
    hist, _ = np.histogram(vals, bins=bins)
    print(f"\n  {gs} distribution (n={len(vals)}):")
    for b in range(n_bins):
        bar = '#' * (hist[b] * 40 // max(1, max(hist)))
        print(f"    [{bins[b]:.2f}-{bins[b+1]:.2f}): {hist[b]:4d} {bar}")

# =====================================================================
# SECTION 11: Hopf-projected area per gauge sector
# =====================================================================
print("\n--- SECTION 11: Hopf-projected area per sector ---")

def hopf_map(q):
    """Hopf map S^3 -> S^2."""
    qi = qmul(q, np.array([0, 1, 0, 0]))
    qinv_q = np.array([q[0], -q[1], -q[2], -q[3]])
    result = qmul(qi, qinv_q)
    return result[1:4]

hopf_pts = np.array([hopf_map(V[i]) for i in range(n_vert)])

def spherical_area_S2(p1, p2, p3):
    """Signed spherical area of triangle on S^2."""
    cross = np.cross(p2, p3)
    num = abs(np.dot(p1, cross))
    denom = 1 + np.dot(p1, p2) + np.dot(p2, p3) + np.dot(p1, p3)
    if abs(denom) < 1e-10:
        return PI
    return 2 * math.atan2(num, denom)

hopf_areas = np.zeros(len(faces))
for fi, (i, j, k) in enumerate(faces):
    hopf_areas[fi] = spherical_area_S2(hopf_pts[i], hopf_pts[j], hopf_pts[k])

print(f"  Hopf-projected area per sector:")
print(f"  {'Sector':>8s} {'<area>':>10s} {'std':>10s} {'n':>6s}")
print(f"  {'-'*40}")
for gs in ['U1', 'SU2', 'SU3']:
    mask = [fi for fi in range(len(faces)) if face_sectors[fi] == gs]
    vals = hopf_areas[mask]
    print(f"  {gs:>8s} {vals.mean():10.6f} {vals.std():10.6f} {len(vals):6d}")

# Distinct Hopf areas per sector
for gs in ['U1', 'SU2', 'SU3']:
    mask = [fi for fi in range(len(faces)) if face_sectors[fi] == gs]
    vals = np.round(hopf_areas[mask], 4)
    distinct = sorted(set(vals))
    if len(distinct) <= 8:
        print(f"    {gs} distinct areas: {distinct}")

# Relationship: Hopf area ~ Omega * fiber_component?
print(f"\n  Relation: Hopf area vs Omega * |n.T_i|:")
print(f"    Correlation = {np.corrcoef(hopf_areas, fiber_components * Omega)[0,1]:.6f}")

# =====================================================================
# SECTION 12: Summary
# =====================================================================
print("\n" + "=" * 72)
print("SUMMARY: EXP-314")
print("=" * 72)

print(f"""
TRIANGLE CLASSIFICATION:
  (A,C,C) -> U1:  {face_type_counts.get(('A','C','C'), 0):4d} faces
  (B,C,C) -> SU2: {face_type_counts.get(('B','C','C'), 0):4d} faces
  (C,C,C) -> SU3: {face_type_counts.get(('C','C','C'), 0):4d} faces
  Total: {len(faces)} faces

HOLONOMY ANGLE:
  Omega = {Omega:.6f} rad = {math.degrees(Omega):.4f} deg
  Same for ALL triangles (regular polytope, face-transitive under H4)

FIBER COMPONENT OF ROTATION AXIS |n . T_i|:""")

for gs in ['U1', 'SU2', 'SU3']:
    print(f"  {gs:>8s}: {sector_means[gs]:.6f}")

aniso = max(sector_means.values()) - min(sector_means.values())
print(f"\n  Max anisotropy: {aniso:.6f}")
if aniso > 0.001:
    print(f"  ANISOTROPIC: gauge sectors have different fiber alignment!")
    print(f"  This means different Berry phases per sector.")
else:
    print(f"  ISOTROPIC: gauge sectors have same fiber alignment.")
    print(f"  Triangle normals are uniformly distributed across sectors.")

print(f"""
BERRY PHASE PER SECTOR:
  U1:  {b_U1:.6f}
  SU2: {b_SU2:.6f}
  SU3: {b_SU3:.6f}
""")

print("DONE.")

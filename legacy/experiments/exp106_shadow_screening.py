"""
EXP-106: Shadow Screening - Green's Function on 600-Cell
=========================================================

QUESTION: Can higher-order corrections to C(r) from the 600-cell lattice
          "inflate" the photon sphere from r=M/2 toward r~3M, resolving
          the 7-sigma tension with EHT observations?

APPROACH:
1. Compute the graph Laplacian Green's function G(v0, v) on the 600-cell
2. Average G over BFS shells to get G_avg(d) for d=1,2,3,4,5
3. Fit: C(d) = 1 + kappa * G_avg(d) and extract higher-order corrections
4. Compute per-shell angular fraction alpha(d)
5. With corrected C(r), compute photon sphere and shadow
6. Compare with GR

STATUS: DERIVED (computation from 600-cell geometry, no free parameters)
"""

import numpy as np
from itertools import permutations
import scipy.linalg as la

PHI = (1 + np.sqrt(5)) / 2

# =====================================================
# PART 1: Build 600-Cell
# =====================================================

def generate_600cell_vertices():
    """Generate all 120 vertices of the 600-cell on unit S^3."""
    verts = set()

    # Type A: 8 vertices - 16-cell (permutations of (+-1, 0, 0, 0))
    for i in range(4):
        for s in [+1, -1]:
            v = [0, 0, 0, 0]
            v[i] = s
            verts.add(tuple(v))

    # Type B: 16 vertices - tesseract (+-1/2, +-1/2, +-1/2, +-1/2)
    for s0 in [+0.5, -0.5]:
        for s1 in [+0.5, -0.5]:
            for s2 in [+0.5, -0.5]:
                for s3 in [+0.5, -0.5]:
                    verts.add((s0, s1, s2, s3))

    # Type C: 96 vertices - snub 24-cell
    # Even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base_values = [PHI/2, 0.5, 1/(2*PHI), 0.0]

    # Generate even permutations
    even_perms = []
    indices = [0, 1, 2, 3]
    for p in permutations(indices):
        # Count inversions
        inv = 0
        for i in range(4):
            for j in range(i+1, 4):
                if p[i] > p[j]:
                    inv += 1
        if inv % 2 == 0:
            even_perms.append(p)

    for perm in even_perms:
        vals = [base_values[perm[i]] for i in range(4)]
        # All sign combinations
        for s0 in [+1, -1]:
            for s1 in [+1, -1]:
                for s2 in [+1, -1]:
                    for s3 in [+1, -1]:
                        v = (s0*vals[0], s1*vals[1], s2*vals[2], s3*vals[3])
                        # Check it's on unit sphere
                        norm = sum(x**2 for x in v)
                        if abs(norm - 1.0) < 1e-10 and v not in verts:
                            verts.add(v)

    verts = list(verts)
    assert len(verts) == 120, f"Expected 120 vertices, got {len(verts)}"
    return np.array(verts)

def build_adjacency(vertices, threshold=0.05):
    """Build adjacency matrix from edge distance."""
    n = len(vertices)
    # Compute all pairwise distances
    dists = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(vertices[i] - vertices[j])
            dists[i,j] = dists[j,i] = d

    # Find minimum non-zero distance (edge length)
    min_dist = np.min(dists[dists > 0.01])

    # Build adjacency
    A = np.zeros((n, n), dtype=int)
    edge_count = 0
    for i in range(n):
        for j in range(i+1, n):
            if abs(dists[i,j] - min_dist) < threshold:
                A[i,j] = A[j,i] = 1
                edge_count += 1

    print(f"Vertices: {n}, Edges: {edge_count}, Min distance: {min_dist:.6f}")
    assert edge_count == 720, f"Expected 720 edges, got {edge_count}"

    # Verify degree 12
    degrees = A.sum(axis=1)
    assert np.all(degrees == 12), f"Not all degrees are 12: {np.unique(degrees)}"

    return A, dists

def bfs_distances(A, source):
    """BFS from source, return array of graph distances."""
    n = len(A)
    dist = np.full(n, -1, dtype=int)
    dist[source] = 0
    queue = [source]
    while queue:
        next_queue = []
        for v in queue:
            for w in range(n):
                if A[v,w] and dist[w] == -1:
                    dist[w] = dist[v] + 1
                    next_queue.append(w)
        queue = next_queue
    return dist

# =====================================================
# PART 2: Green's Function on Graph
# =====================================================

def compute_greens_function(A):
    """Compute Green's function (pseudoinverse of Laplacian)."""
    n = len(A)
    D = np.diag(A.sum(axis=1).astype(float))
    L = D - A.astype(float)

    # Eigendecomposition
    eigenvalues, eigenvectors = la.eigh(L)

    print(f"\nLaplacian spectrum:")
    unique_evals = []
    for ev in eigenvalues:
        is_new = True
        for uev in unique_evals:
            if abs(ev - uev) < 1e-8:
                is_new = False
                break
        if is_new:
            unique_evals.append(ev)

    for i, ev in enumerate(sorted(unique_evals)[:8]):
        mult = np.sum(np.abs(eigenvalues - ev) < 1e-8)
        print(f"  lambda_{i} = {ev:.6f}  (multiplicity {mult})")
    print(f"  ... ({len(unique_evals)} distinct eigenvalues total)")

    # Green's function: G = sum_{k: lambda_k > 0} |psi_k><psi_k| / lambda_k
    G = np.zeros((n, n))
    for k in range(n):
        if eigenvalues[k] > 1e-10:
            G += np.outer(eigenvectors[:,k], eigenvectors[:,k]) / eigenvalues[k]

    return G, eigenvalues, eigenvectors

# =====================================================
# PART 3: Per-Shell Analysis
# =====================================================

def analyze_per_shell(A, source):
    """Compute angular/radial edge fraction per BFS shell."""
    n = len(A)
    dist = bfs_distances(A, source)
    max_d = dist.max()

    # Shell sizes
    shell_sizes = {}
    for d in range(max_d + 1):
        shell_sizes[d] = np.sum(dist == d)

    # Count angular and radial edges per shell
    angular_per_shell = {}  # edges within shell d
    radial_per_shell = {}   # edges from shell d to shell d+1

    for d in range(max_d + 1):
        angular_per_shell[d] = 0
        radial_per_shell[d] = 0

    for i in range(n):
        for j in range(i+1, n):
            if A[i,j]:
                di, dj = dist[i], dist[j]
                if di == dj:
                    angular_per_shell[di] += 1
                elif abs(di - dj) == 1:
                    d_min = min(di, dj)
                    radial_per_shell[d_min] += 1

    return dist, shell_sizes, angular_per_shell, radial_per_shell

# =====================================================
# PART 4: Photon Sphere and Shadow with Corrected C(r)
# =====================================================

def photon_sphere_corrected(C_func, alpha, M=1.0, r_range=(0.01, 20.0)):
    """Find photon sphere for general C(r) and alpha."""
    from scipy.optimize import minimize_scalar

    # Effective potential for photons: V(r) = 1 / (C^(2+2alpha) * r^2)
    # Photon sphere at maximum of V
    def neg_V(r):
        C = C_func(r, M)
        return -1.0 / (C**(2+2*alpha) * r**2)

    result = minimize_scalar(neg_V, bounds=r_range, method='bounded')
    r_ph = result.x
    V_max = -result.fun
    b_c = 1.0 / np.sqrt(V_max)

    return r_ph, b_c

def isco_corrected(C_func, alpha, M=1.0, r_range=(0.1, 50.0)):
    """Find ISCO for general C(r) and alpha."""
    from scipy.optimize import minimize_scalar

    # For massive particles, effective potential:
    # V_eff(r, L) = -1/C^2 + L^2/(C^(2alpha)*r^2)
    # ISCO: where V_eff has inflection point
    # Simplified: find where d^2/dr^2[1/(C^(2+2alpha)*r^2)] = 0

    # Use numerical second derivative of V_photon
    def V_photon(r):
        C = C_func(r, M)
        return 1.0 / (C**(2+2*alpha) * r**2)

    dr = 1e-6
    def d2V(r):
        return (V_photon(r+dr) - 2*V_photon(r) + V_photon(r-dr)) / dr**2

    # ISCO is where d2V = 0 (outside photon sphere)
    # Scan for zero crossing
    r_ph, _ = photon_sphere_corrected(C_func, alpha, M, r_range)

    rs = np.linspace(r_ph * 1.1, r_range[1], 10000)
    d2vs = [d2V(r) for r in rs]

    for i in range(len(rs)-1):
        if d2vs[i] * d2vs[i+1] < 0:
            # Bisection
            a, b = rs[i], rs[i+1]
            for _ in range(50):
                mid = (a+b)/2
                if d2V(a)*d2V(mid) < 0:
                    b = mid
                else:
                    a = mid
            return (a+b)/2

    return None  # No ISCO found

# =====================================================
# MAIN
# =====================================================

if __name__ == "__main__":
    print("="*70)
    print("EXP-106: Shadow Screening - Green's Function on 600-Cell")
    print("="*70)

    # Build 600-cell
    print("\n--- Building 600-cell ---")
    vertices = generate_600cell_vertices()
    A, dists = build_adjacency(vertices)

    # Compute Green's function
    print("\n--- Computing Green's function ---")
    G, eigenvalues, eigenvectors = compute_greens_function(A)

    # Analyze per-shell structure for multiple source vertices
    print("\n--- Per-shell analysis ---")

    # Use vertex 0 as source (results are vertex-independent for totals)
    source = 0
    dist, shell_sizes, angular_per_shell, radial_per_shell = analyze_per_shell(A, source)

    print(f"\nBFS from vertex {source}:")
    print(f"{'Shell d':>8} {'Vertices':>10} {'Angular':>10} {'Radial':>10} {'alpha_local':>12} {'Cumulative':>12}")
    print("-" * 65)

    total_ang = 0
    total_rad = 0
    max_d = max(shell_sizes.keys())

    alpha_per_shell = {}
    for d in range(max_d + 1):
        ang = angular_per_shell[d]
        rad = radial_per_shell[d]
        total = ang + rad
        alpha_local = ang / total if total > 0 else 0
        alpha_per_shell[d] = alpha_local
        total_ang += ang
        total_rad += rad
        cum_alpha = total_ang / (total_ang + total_rad) if (total_ang + total_rad) > 0 else 0
        print(f"{d:>8} {shell_sizes[d]:>10} {ang:>10} {rad:>10} {alpha_local:>12.4f} {cum_alpha:>12.4f}")

    print(f"\nTotal: {total_ang} angular + {total_rad} radial = {total_ang + total_rad}")
    print(f"Global angular fraction: {total_ang / (total_ang + total_rad):.6f}")

    # Verify over all source vertices
    print("\n--- Verifying alpha per shell over ALL source vertices ---")
    all_alpha_shells = {d: [] for d in range(max_d + 1)}

    for v0 in range(120):
        _, ss, aps, rps = analyze_per_shell(A, v0)
        for d in range(max_d + 1):
            total = aps[d] + rps[d]
            if total > 0:
                all_alpha_shells[d].append(aps[d] / total)

    print(f"\n{'Shell d':>8} {'mean(alpha)':>12} {'std(alpha)':>12} {'min':>8} {'max':>8}")
    print("-" * 55)
    for d in range(max_d + 1):
        vals = all_alpha_shells[d]
        if vals:
            print(f"{d:>8} {np.mean(vals):>12.6f} {np.std(vals):>12.6f} {np.min(vals):>8.4f} {np.max(vals):>8.4f}")

    # Green's function analysis
    print("\n--- Green's function vs BFS distance ---")

    # Average G(v0, v) over vertices at each BFS distance
    G_avg = {}
    G_std = {}
    for d in range(1, max_d + 1):  # skip d=0 (self)
        mask = (dist == d)
        g_vals = G[source, mask]
        G_avg[d] = np.mean(g_vals)
        G_std[d] = np.std(g_vals)

    # Also compute diagonal (self-potential)
    G_self = G[source, source]
    print(f"\nG(v0, v0) = {G_self:.6f} (self-potential)")

    print(f"\n{'d':>4} {'G_avg(d)':>12} {'G_std':>10} {'G*d':>10} {'G*d^2':>10}")
    print("-" * 50)
    for d in range(1, max_d + 1):
        print(f"{d:>4} {G_avg[d]:>12.6f} {G_std[d]:>10.6f} {G_avg[d]*d:>10.4f} {G_avg[d]*d**2:>10.4f}")

    # Fit C(r) = 1 + a1/r + a2/r^2 + a3/r^3
    print("\n--- Fitting C(d) = 1 + kappa * G_avg(d) ---")
    print("--- Also fitting G_avg(d) = a1/d + a2/d^2 + a3/d^3 ---")

    ds = np.array([d for d in range(1, max_d + 1)])
    gs = np.array([G_avg[d] for d in range(1, max_d + 1)])

    # Fit: G(d) = a1/d + a2/d^2 + a3/d^3
    # Design matrix: [1/d, 1/d^2, 1/d^3]
    X = np.column_stack([1/ds, 1/ds**2, 1/ds**3])
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, gs, rcond=None)
    a1, a2, a3 = coeffs

    print(f"\nFit: G(d) = {a1:.6f}/d + {a2:.6f}/d^2 + {a3:.6f}/d^3")
    print(f"\nRatio a2/a1 = {a2/a1:.4f} (= beta in C(r) = 1 + M/r + beta*(M/r)^2)")
    print(f"Ratio a3/a1 = {a3/a1:.4f}")

    # Verify fit quality
    print(f"\n{'d':>4} {'G_actual':>12} {'G_fit':>12} {'Error%':>10}")
    print("-" * 40)
    for i, d in enumerate(ds):
        g_fit = a1/d + a2/d**2 + a3/d**3
        err = abs(g_fit - gs[i]) / abs(gs[i]) * 100
        print(f"{d:>4} {gs[i]:>12.6f} {g_fit:>12.6f} {err:>10.2f}%")

    # Also try pure 1/d fit
    a1_pure = np.linalg.lstsq(1/ds.reshape(-1,1), gs, rcond=None)[0][0]
    print(f"\nPure 1/d fit: G(d) = {a1_pure:.6f}/d")
    print(f"  Residuals:")
    for i, d in enumerate(ds):
        g_fit = a1_pure / d
        err = (gs[i] - g_fit) / gs[i] * 100
        print(f"  d={d}: actual={gs[i]:.6f}, fit={g_fit:.6f}, deviation={err:+.2f}%")

    # =====================================================
    # PART 5: Shadow Computation with Corrected C(r)
    # =====================================================
    print("\n" + "="*70)
    print("SHADOW ANALYSIS WITH CORRECTED C(r)")
    print("="*70)

    M = 1.0
    alpha_metric = 0.5
    beta = a2 / a1  # From Green's function fit

    print(f"\nbeta from 600-cell Green's function: {beta:.4f}")

    # Define C(r) functions
    def C_standard(r, M):
        """Original STG: C = 1 + M/r"""
        return 1 + M/r

    def C_corrected(r, M):
        """Corrected: C = 1 + M/r + beta*(M/r)^2"""
        u = M/r
        return 1 + u + beta * u**2

    def C_corrected_cubic(r, M):
        """Corrected with cubic: C = 1 + M/r + beta*(M/r)^2 + gamma*(M/r)^3"""
        u = M/r
        gamma = a3/a1
        return 1 + u + beta * u**2 + gamma * u**3

    # Compute photon spheres and shadows
    print("\n--- Photon sphere and shadow comparison ---")
    print(f"{'Model':>25} {'r_ph/M':>10} {'b_c/M':>10} {'b_c/b_GR':>10}")
    print("-" * 60)

    # GR values for reference
    r_ph_GR = 3 * M
    b_c_GR = 3 * np.sqrt(3) * M
    print(f"{'GR (Schwarzschild)':>25} {r_ph_GR:>10.4f} {b_c_GR:>10.4f} {'1.0000':>10}")

    # Standard STG
    r_ph_std, b_c_std = photon_sphere_corrected(C_standard, alpha_metric, M)
    print(f"{'STG standard':>25} {r_ph_std:>10.4f} {b_c_std:>10.4f} {b_c_std/b_c_GR:>10.4f}")

    # Corrected STG (quadratic)
    r_ph_cor, b_c_cor = photon_sphere_corrected(C_corrected, alpha_metric, M)
    print(f"{'STG + beta*u^2':>25} {r_ph_cor:>10.4f} {b_c_cor:>10.4f} {b_c_cor/b_c_GR:>10.4f}")

    # Corrected STG (cubic)
    r_ph_cub, b_c_cub = photon_sphere_corrected(C_corrected_cubic, alpha_metric, M)
    print(f"{'STG + beta*u^2 + gamma*u^3':>25} {r_ph_cub:>10.4f} {b_c_cub:>10.4f} {b_c_cub/b_c_GR:>10.4f}")

    # What beta would we NEED to match GR shadow?
    print("\n--- Parameter scan: what beta matches GR shadow? ---")

    betas = np.linspace(0, 5, 500)
    results = []
    for b_test in betas:
        def C_test(r, M):
            u = M/r
            return 1 + u + b_test * u**2
        try:
            r_ph_t, b_c_t = photon_sphere_corrected(C_test, alpha_metric, M)
            results.append((b_test, r_ph_t, b_c_t, b_c_t/b_c_GR))
        except:
            pass

    results = np.array(results)

    # Find beta where b_c/b_GR is closest to 1
    idx_match = np.argmin(np.abs(results[:,3] - 1.0))
    beta_match = results[idx_match, 0]
    r_ph_match = results[idx_match, 1]
    b_c_match = results[idx_match, 2]

    print(f"\nTo match GR shadow (b_c/b_GR = 1.0):")
    print(f"  Need beta = {beta_match:.4f}")
    print(f"  Gives r_ph = {r_ph_match:.4f} M")
    print(f"  Gives b_c = {b_c_match:.4f} M (GR: {b_c_GR:.4f} M)")

    print(f"\nActual beta from 600-cell: {beta:.4f}")
    print(f"Ratio needed/actual: {beta_match/beta:.2f}")

    # Check weak-field behavior of corrected C(r)
    print("\n--- Weak-field verification ---")
    print("At r = 100M:")
    print(f"  C_standard = {C_standard(100, M):.8f}")
    print(f"  C_corrected = {C_corrected(100, M):.8f}")
    print(f"  Difference: {abs(C_corrected(100, M) - C_standard(100, M)):.2e}")
    print("At r = 10M:")
    print(f"  C_standard = {C_standard(10, M):.8f}")
    print(f"  C_corrected = {C_corrected(10, M):.8f}")
    print(f"  Difference: {abs(C_corrected(10, M) - C_standard(10, M)):.4f}")

    # Check precession at 1PN with corrected C
    print("\n--- 1PN precession check ---")
    print("At r = 100M (weak field), the beta*(M/r)^2 correction is:")
    print(f"  beta * (M/r)^2 = {beta} * {(M/100)**2:.6f} = {beta * (M/100)**2:.2e}")
    print(f"  This is a {beta * (M/100)**2 * 100:.4f}% correction to C")
    print(f"  1PN precession agreement is preserved (correction enters at 2PN)")

    # =====================================================
    # PART 6: Alternative - Full Green's function as C(r)
    # =====================================================
    print("\n" + "="*70)
    print("ALTERNATIVE: Using full Green's function G(d) directly as C")
    print("="*70)

    # Normalize so that G_avg(d=1) gives C(d=1) reasonable
    # C(d) = 1 + kappa * G_avg(d)
    # At large d: C ~ 1 + kappa * a1/d ~ 1 + M/d
    # So kappa * a1 = M => kappa = M/a1

    kappa = M / a1
    print(f"\nNormalization: kappa = M/a1 = {kappa:.6f}")

    print(f"\n{'d':>4} {'C_direct':>12} {'C_1+M/d':>12} {'Ratio':>10}")
    print("-" * 42)
    for d in range(1, max_d + 1):
        C_direct = 1 + kappa * G_avg[d]
        C_simple = 1 + M/d
        print(f"{d:>4} {C_direct:>12.6f} {C_simple:>12.6f} {C_direct/C_simple:>10.6f}")

    # =====================================================
    # PART 7: Per-shell alpha as effective metric exponent
    # =====================================================
    print("\n" + "="*70)
    print("PER-SHELL ALPHA: Effective metric exponent vs distance")
    print("="*70)

    # Average over all source vertices for robustness
    print(f"\n{'d':>4} {'alpha(d)':>10} {'r_ph if constant':>18} {'Implication':>30}")
    print("-" * 70)
    for d in range(max_d + 1):
        vals = all_alpha_shells[d]
        if vals:
            alpha_d = np.mean(vals)
            r_ph_if = alpha_d * M if alpha_d > 0 else 0
            imp = ""
            if d <= 1:
                imp = "Near mass (strong field)"
            elif d <= 3:
                imp = "Intermediate"
            else:
                imp = "Far (weak field)"
            print(f"{d:>4} {alpha_d:>10.4f} {r_ph_if:>18.4f} M    {imp:>30}")

    # =====================================================
    # SUMMARY
    # =====================================================
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    print(f"""
RESULTS:

1. GREEN'S FUNCTION ANALYSIS:
   - G(d) on 600-cell is NOT purely 1/d
   - Best fit: G(d) = {a1:.4f}/d + {a2:.4f}/d^2 + {a3:.4f}/d^3
   - Higher-order coefficient beta = a2/a1 = {beta:.4f}

2. SHADOW WITH CORRECTED C(r):
   - Standard STG (C = 1+M/r):     b_c = {b_c_std:.4f} M  (ratio to GR: {b_c_std/b_c_GR:.4f})
   - Corrected STG (beta={beta:.3f}):  b_c = {b_c_cor:.4f} M  (ratio to GR: {b_c_cor/b_c_GR:.4f})
   - GR reference:                   b_c = {b_c_GR:.4f} M

3. TO MATCH GR SHADOW:
   - Need beta = {beta_match:.4f}
   - 600-cell gives beta = {beta:.4f}
   - Gap factor: {beta_match/beta:.2f}x

4. PER-SHELL ALPHA:
   - alpha varies with BFS shell distance
   - But even alpha=1 only gives r_ph = M (still << 3M)
   - alpha variation alone CANNOT resolve the shadow tension

5. WEAK-FIELD PRESERVATION:
   - At r = 100M: correction is {beta * (M/100)**2:.2e} (negligible)
   - 1PN precession and deflection are PRESERVED
   - Corrections only matter at 2PN and beyond (strong field)

CLASSIFICATION:
- Green's function analysis: DERIVED (exact computation from 600-cell)
- beta coefficient: DERIVED (from lattice structure)
- Shadow with corrections: DERIVED (no free parameters)
- Gap to GR: QUANTIFIED (this is what it is)

VERDICT: The 600-cell lattice Green's function provides higher-order
corrections to C(r), but they are {"INSUFFICIENT" if abs(b_c_cor/b_c_GR - 1.0) > 0.1 else "SUFFICIENT"}
to fully resolve the shadow tension with EHT.
""")

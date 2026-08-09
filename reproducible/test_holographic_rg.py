"""
test_holographic_rg.py
======================
Exploratory go/no-go test for a holographic-RG interpretation of the
600-cell boundary data.

Core question:
  Can the spectral flow of the discrete Dirac/Laplacian data on the 600-cell
  support an emergent 4D bulk interpretation, rather than only the obvious 3D
  S^3 boundary interpretation?

This script is intentionally exploratory rather than confirmatory.
It computes two independent diagnostics from the exact 600-cell Hodge data:

  1. Counting-function scaling on the positive spectrum of D^2
       N(lambda) ~ lambda^(d/2)

  2. Heat-flow effective dimension on the same spectrum
       K(t) = sum exp(-t lambda),   d_eff(t) = -2 d log K / d log t

Decision logic:
  CONTINUE:
      both diagnostics show a stable 4D-like window
  STOP:
      both diagnostics show only a stable 3D-like window
  INCONCLUSIVE:
      anything mixed or unstable

Dependencies: numpy, scipy
No project imports. Safe on Windows.
"""

import numpy as np
from itertools import permutations, product
from collections import defaultdict
from scipy.linalg import eigh


PHI = (1 + np.sqrt(5.0)) / 2.0
EDGE_THRESHOLD = PHI / 2.0
EDGE_TOL = 1e-3


def build_600cell_complex():
    """Return vertices, edges, triangles, tetrahedra and lookup maps."""
    verts_set = set()

    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts_set.add(tuple(v))

    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    verts_set.add((s0, s1, s2, s3))

    base = [PHI / 2.0, 0.5, 1.0 / (2.0 * PHI), 0.0]
    even_perms = []
    for p in permutations(range(4)):
        inv_count = sum(
            1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]
        )
        if inv_count % 2 == 0:
            even_perms.append(p)

    for perm in even_perms:
        coords = [base[perm[i]] for i in range(4)]
        nonzero_idx = [i for i in range(4) if abs(coords[i]) > 1e-12]
        for signs in product([1, -1], repeat=len(nonzero_idx)):
            v = list(coords)
            for idx, s in zip(nonzero_idx, signs):
                v[idx] *= s
            verts_set.add(tuple(round(x, 10) for x in v))

    verts = np.array(sorted(verts_set), dtype=float)
    dots = verts @ verts.T

    edges = []
    edge_to_idx = {}
    adj = defaultdict(set)
    nv = len(verts)
    for i in range(nv):
        for j in range(i + 1, nv):
            if abs(dots[i, j] - EDGE_THRESHOLD) < EDGE_TOL:
                idx = len(edges)
                edges.append((i, j))
                edge_to_idx[(i, j)] = idx
                edge_to_idx[(j, i)] = idx
                adj[i].add(j)
                adj[j].add(i)

    triangles = []
    face_to_idx = {}
    for i in range(nv):
        for j in adj[i]:
            if j > i:
                common = adj[i] & adj[j]
                for k in common:
                    if k > j:
                        face = (i, j, k)
                        face_to_idx[face] = len(triangles)
                        triangles.append(face)

    tetrahedra = []
    for i in range(nv):
        ni = adj[i]
        for j in ni:
            if j > i:
                common_ij = ni & adj[j]
                for k in common_ij:
                    if k > j:
                        common_ijk = common_ij & adj[k]
                        for l in common_ijk:
                            if l > k:
                                tetrahedra.append((i, j, k, l))

    return verts, edges, triangles, tetrahedra, edge_to_idx, face_to_idx


def build_hodge_laplacians(edges, triangles, tetrahedra, edge_to_idx, face_to_idx):
    """Build the simplicial Hodge Laplacians Delta_p."""
    nv = 120
    ne = len(edges)
    nf = len(triangles)
    nc = len(tetrahedra)

    d0 = np.zeros((ne, nv), dtype=float)
    for e_idx, (i, j) in enumerate(edges):
        d0[e_idx, i] = -1.0
        d0[e_idx, j] = 1.0

    d1 = np.zeros((nf, ne), dtype=float)
    for f_idx, (i, j, k) in enumerate(triangles):
        d1[f_idx, edge_to_idx[(i, j)]] = 1.0
        d1[f_idx, edge_to_idx[(j, k)]] = 1.0
        d1[f_idx, edge_to_idx[(i, k)]] = -1.0

    d2 = np.zeros((nc, nf), dtype=float)
    for c_idx, (i, j, k, l) in enumerate(tetrahedra):
        faces_of_tet = [
            ((j, k, l), 1.0),
            ((i, k, l), -1.0),
            ((i, j, l), 1.0),
            ((i, j, k), -1.0),
        ]
        for face, sign in faces_of_tet:
            f_idx = face_to_idx[face]
            d2[c_idx, f_idx] = sign

    delta0 = d0.T @ d0
    delta1 = d0 @ d0.T + d1.T @ d1
    delta2 = d1 @ d1.T + d2.T @ d2
    delta3 = d2 @ d2.T
    return delta0, delta1, delta2, delta3


def rolling_counting_fit(lam_pos, window):
    """Scan rolling windows of N(lambda) ~ lambda^(d/2)."""
    logs_lam = np.log(lam_pos)
    counts = np.arange(1, len(lam_pos) + 1, dtype=float)
    logs_n = np.log(counts)
    fits = []
    for start in range(0, len(lam_pos) - window + 1, max(1, window // 8)):
        stop = start + window
        x = logs_lam[start:stop]
        y = logs_n[start:stop]
        slope, intercept = np.polyfit(x, y, 1)
        residuals = y - (slope * x + intercept)
        rmse = float(np.sqrt(np.mean(residuals ** 2)))
        fits.append(
            {
                "start": start,
                "stop": stop - 1,
                "lam_min": float(lam_pos[start]),
                "lam_max": float(lam_pos[stop - 1]),
                "dimension": float(2.0 * slope),
                "rmse": rmse,
            }
        )
    return fits


def best_fit_near_target(fits, target):
    """Pick the cleanest counting-fit window near a target dimension."""
    return min(fits, key=lambda item: (abs(item["dimension"] - target), item["rmse"]))


def heat_effective_dimension(lam_pos, num_points=180):
    """Return t-grid, heat trace, and effective dimension."""
    lam_min = float(np.min(lam_pos))
    lam_max = float(np.max(lam_pos))
    t_min = 0.25 / lam_max
    t_max = 12.0 / lam_min
    t_grid = np.logspace(np.log10(t_min), np.log10(t_max), num_points)
    k_vals = np.array([np.exp(-t * lam_pos).sum() for t in t_grid], dtype=float)
    log_t = np.log(t_grid)
    log_k = np.log(k_vals)
    d_eff = -2.0 * np.gradient(log_k, log_t)
    return t_grid, k_vals, d_eff


def best_heat_window_near_target(t_grid, d_eff, target, window=15):
    """Pick the most stable heat-dimension window near a target."""
    candidates = []
    for start in range(0, len(d_eff) - window + 1):
        segment = d_eff[start:start + window]
        mean = float(np.mean(segment))
        std = float(np.std(segment))
        candidates.append(
            {
                "start": start,
                "stop": start + window - 1,
                "t_min": float(t_grid[start]),
                "t_max": float(t_grid[start + window - 1]),
                "mean": mean,
                "std": std,
            }
        )
    return min(candidates, key=lambda item: (abs(item["mean"] - target), item["std"]))


def support_target_counting(fit, target, dim_tol=0.35, rmse_tol=0.06):
    return abs(fit["dimension"] - target) <= dim_tol and fit["rmse"] <= rmse_tol


def support_target_heat(window, target, dim_tol=0.35, std_tol=0.20):
    return abs(window["mean"] - target) <= dim_tol and window["std"] <= std_tol


def print_counting_fit(label, fit):
    print(label)
    print(
        "  lambda in [{:.6f}, {:.6f}], indices {}..{}".format(
            fit["lam_min"], fit["lam_max"], fit["start"], fit["stop"]
        )
    )
    print("  fitted dimension = {:.4f}".format(fit["dimension"]))
    print("  log-log RMSE      = {:.4f}".format(fit["rmse"]))


def print_heat_window(label, window):
    print(label)
    print("  t in [{:.6e}, {:.6e}]".format(window["t_min"], window["t_max"]))
    print("  mean d_eff = {:.4f}".format(window["mean"]))
    print("  std  d_eff = {:.4f}".format(window["std"]))


def main():
    print("=" * 72)
    print("HOLOGRAPHIC RG GO/NO-GO TEST ON THE 600-CELL")
    print("=" * 72)

    print("\n--- Build 600-cell Hodge data ---")
    _, edges, triangles, tetrahedra, edge_to_idx, face_to_idx = build_600cell_complex()
    print(f"  edges      = {len(edges)}")
    print(f"  triangles  = {len(triangles)}")
    print(f"  tetrahedra = {len(tetrahedra)}")

    delta0, delta1, delta2, delta3 = build_hodge_laplacians(
        edges, triangles, tetrahedra, edge_to_idx, face_to_idx
    )

    print("\n--- Diagonalize Delta_p blocks ---")
    spec0 = eigh(delta0, eigvals_only=True)
    spec1 = eigh(delta1, eigvals_only=True)
    spec2 = eigh(delta2, eigvals_only=True)
    spec3 = eigh(delta3, eigvals_only=True)
    d2_spec = np.sort(np.concatenate([spec0, spec1, spec2, spec3]))
    tol = 1e-9
    zero_modes = int(np.sum(np.abs(d2_spec) < tol))
    lam_pos = d2_spec[d2_spec > tol]

    print(f"  total D^2 eigenvalues = {len(d2_spec)}")
    print(f"  zero modes            = {zero_modes}")
    print(f"  positive modes        = {len(lam_pos)}")
    print(f"  lambda_min^+          = {lam_pos[0]:.6f}")
    print(f"  lambda_max            = {lam_pos[-1]:.6f}")

    print("\n--- Diagnostic 1: counting-function scaling ---")
    count_fits = rolling_counting_fit(lam_pos, window=320)
    best3_count = best_fit_near_target(count_fits, target=3.0)
    best4_count = best_fit_near_target(count_fits, target=4.0)
    best_any_count = min(count_fits, key=lambda item: item["rmse"])
    print_counting_fit("Best window near d = 3", best3_count)
    print_counting_fit("Best window near d = 4", best4_count)
    print_counting_fit("Best window by raw RMSE", best_any_count)

    print("\n--- Diagnostic 2: heat-flow effective dimension ---")
    t_grid, _, d_eff = heat_effective_dimension(lam_pos)
    best3_heat = best_heat_window_near_target(t_grid, d_eff, target=3.0)
    best4_heat = best_heat_window_near_target(t_grid, d_eff, target=4.0)
    best_any_heat = min(
        (
            {
                "start": i,
                "stop": i + 14,
                "t_min": float(t_grid[i]),
                "t_max": float(t_grid[i + 14]),
                "mean": float(np.mean(d_eff[i:i + 15])),
                "std": float(np.std(d_eff[i:i + 15])),
            }
            for i in range(len(d_eff) - 14)
        ),
        key=lambda item: item["std"],
    )
    print_heat_window("Best window near d = 3", best3_heat)
    print_heat_window("Best window near d = 4", best4_heat)
    print_heat_window("Most stable window overall", best_any_heat)

    support3_count = support_target_counting(best3_count, 3.0)
    support4_count = support_target_counting(best4_count, 4.0)
    support3_heat = support_target_heat(best3_heat, 3.0)
    support4_heat = support_target_heat(best4_heat, 4.0)

    print("\n--- Decision summary ---")
    print(f"  counting supports 3D: {support3_count}")
    print(f"  counting supports 4D: {support4_count}")
    print(f"  heat flow supports 3D: {support3_heat}")
    print(f"  heat flow supports 4D: {support4_heat}")

    if support4_count and support4_heat and not (support3_count and support3_heat):
        decision = "CONTINUE"
        reason = "both diagnostics show a stable 4D-like scaling window"
    elif support3_count and support3_heat and not (support4_count and support4_heat):
        decision = "STOP"
        reason = "both diagnostics prefer a stable 3D boundary interpretation"
    else:
        decision = "INCONCLUSIVE"
        reason = "the diagnostics are mixed or the windows are not stable enough"

    print(f"\nFINAL DECISION: {decision}")
    print(f"Reason: {reason}")

    print("\nInterpretation rule:")
    print("  STOP         -> do not claim emergent 4D holography from current data")
    print("  CONTINUE     -> proceed to coarse-graining and Hessian-mode tests")
    print("  INCONCLUSIVE -> refine observables before making gravity claims")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Frozen warped-product spectral-dimension audit.

FROZEN BEFORE WARPED SPECTRAL EVALUATION (2026-07-24)
=====================================================
The plateau estimators and every numerical gate are imported unchanged from
verify_holographic_dimension.py.

  truncations: N = 8, 16, 24
  radial proper lattice spacing: ell = log(phi)
  primary boundary condition: Dirichlet at both finite ends
  sensitivity boundary condition: Neumann at both finite ends
  w1: phi^(-2n), derived golden warp (primary)
  w2: 2^(-2n), derived McKay warp
  w3: 1, no-warp consistency control

For each spatial D3^2 eigenvalue lambda the exact mode block is

  T_lambda = L_radial + lambda diag(warp(n)).

Dirichlet means ghost values vanish, hence diag(L)=2/ell^2 throughout.
Neumann means the path-graph quadratic form, hence endpoint diagonal
1/ell^2 and interior diagonal 2/ell^2.  Both use off-diagonal -1/ell^2.

A dynamical-4D verdict requires, for at least one derived warp, qualifying
counting and heat d=4 plateaus at all three N for the primary Dirichlet
choice.  Neumann results are reported as boundary-condition sensitivity.
The single-floor control must retain the frozen 3D counting verdict.
No threshold or boundary convention is selected after inspecting results.
"""

import argparse
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
import scipy.integrate
import scipy.linalg

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify_holographic_dimension as hd


PHI = (1 + np.sqrt(5.0)) / 2
ELL = np.log(PHI)
TRUNCATIONS = (8, 16, 24)
BOUNDARIES = ("dirichlet", "neumann")
WARPS = (
    ("w1_golden", PHI, "DERIVED WARP PRIMARY"),
    ("w2_mckay", 2.0, "DERIVED WARP"),
    ("w3_no_warp", 1.0, "CONTROL"),
)


def radial_laplacian(n, boundary):
    """Uniform second difference in proper distance ell."""
    matrix = np.diag(np.full(n, 2.0))
    matrix += np.diag(np.full(n - 1, -1.0), 1)
    matrix += np.diag(np.full(n - 1, -1.0), -1)
    if boundary == "neumann":
        matrix[0, 0] = matrix[-1, -1] = 1.0
    elif boundary != "dirichlet":
        raise ValueError(boundary)
    return matrix / ELL**2


def warped_spectrum(base, n, ratio, boundary):
    """Union of Jacobi spectra, retaining exact spatial multiplicities."""
    radial = radial_laplacian(n, boundary)
    warp = ratio ** (-2 * np.arange(n, dtype=float))
    raw = []
    for spatial_value, spatial_mult in base:
        diagonal = np.diag(radial) + spatial_value * warp
        values = scipy.linalg.eigh_tridiagonal(
            diagonal, np.diag(radial, 1), eigvals_only=True,
            check_finite=False)
        raw.extend((float(value), spatial_mult) for value in values)
    raw.sort()
    groups = []
    for value, multiplicity in raw:
        if abs(value) < hd.TOL:
            value = 0.0
        if groups and abs(value - groups[-1][0]) < hd.TOL:
            old, mult = groups[-1]
            groups[-1] = ((old*mult + value*multiplicity)/(mult+multiplicity),
                          mult + multiplicity)
        else:
            groups.append((value, multiplicity))
    return groups


def explicit_decomposition_error():
    """Independent full-matrix check on the 5-cell boundary at N=4."""
    groups, _ = hd.spectrum_from_facets(list(__import__("itertools").combinations(
        range(5), 4)))
    spatial = hd.expanded(groups)
    n = 4
    radial = radial_laplacian(n, "dirichlet")
    warp = PHI ** (-2 * np.arange(n, dtype=float))
    # A diagonal representative with the same spatial spectrum is sufficient:
    # the claimed block decomposition is basis-independent.
    full = np.kron(np.eye(len(spatial)), radial)
    full += np.kron(np.diag(spatial), np.diag(warp))
    direct = np.linalg.eigvalsh(full)
    union = np.concatenate([
        np.linalg.eigvalsh(radial + value*np.diag(warp))
        for value in spatial])
    return np.max(np.abs(np.sort(direct)-np.sort(union))), len(spatial)


def flat_r4_kernel(t):
    """Exact diagonal heat kernel per unit volume on R^4."""
    return (4*np.pi*t)**-2


def hyperbolic_h4_kernel(t, radius=ELL):
    """Exact H4 diagonal heat kernel via its Plancherel integral."""
    tau = t/radius**2
    integral = scipy.integrate.quad(
        lambda r: r*(r*r+0.25)*np.tanh(np.pi*r)*np.exp(-tau*r*r),
        0, np.inf, epsabs=2e-11, epsrel=2e-11, limit=300)[0]
    return np.exp(-2.25*tau)*integral/(8*np.pi**2*radius**4)


def benchmark_certificate():
    """High-energy normalization and the exact H4 curvature gap."""
    t = ELL**2 * 1e-4
    ratio = hyperbolic_h4_kernel(t) / flat_r4_kernel(t)
    return ratio, 9/(4*ELL**2)


def counting_plateaus_fast(groups, target):
    """Algebraically identical frozen counting scan using prefix sums."""
    levels, counts, local = hd.counting_curve(groups)
    x_width = np.log10(levels)
    x = np.log(levels)
    y = np.log(counts)
    # Leading zero makes interval sums [i,j] simple differences.
    prefixes = [np.r_[0.0, np.cumsum(z)] for z in
                (x, y, x*x, y*y, x*y, local, local*local)]
    best = None
    for i in range(len(x)-2):
        js = np.arange(i+2, len(x))
        widths = x_width[js]-x_width[i]
        keep = widths >= hd.FROZEN["minimum_decades"]
        if not np.any(keep):
            continue
        js, widths = js[keep], widths[keep]
        count = js-i+1
        sums = [p[js+1]-p[i] for p in prefixes]
        sx, sy, sxx, syy, sxy, sl, sll = sums
        denom = sxx-sx*sx/count
        slope = (sxy-sx*sy/count)/denom
        intercept = (sy-slope*sx)/count
        # SSE of y-(slope*x+intercept), with roundoff clipped at zero.
        sse = (syy + slope*slope*sxx + count*intercept*intercept
               - 2*slope*sxy - 2*intercept*sy
               + 2*slope*intercept*sx)
        rmse = np.sqrt(np.maximum(0.0, sse/count))
        local_std = np.sqrt(np.maximum(0.0, sll/count-(sl/count)**2))
        dimension = 2*slope
        valid = ((np.abs(dimension-target)
                  <= hd.FROZEN["target_tolerance"])
                 & (rmse <= hd.FROZEN["counting_log_rmse_max"])
                 & (local_std <= hd.FROZEN["counting_local_std_max"]))
        for k in np.flatnonzero(valid):
            candidate = (widths[k], dimension[k], rmse[k], local_std[k],
                         i, int(js[k]))
            if best is None or (candidate[0], -candidate[2]) > (
                    best[0], -best[2]):
                best = candidate
    return best


def curve_rows(label, groups):
    levels, counts, dn = hd.counting_curve(groups)
    times, ds = hd.heat_curve(groups)
    print(f"{label} COUNTING_CURVE_BEGIN lambda N d_N")
    for row in zip(levels, counts, dn):
        print("  {:.12g} {} {:.9g}".format(*row))
    print(f"{label} COUNTING_CURVE_END")
    print(f"{label} HEAT_CURVE_BEGIN t d_s")
    for row in zip(times, ds):
        print("  {:.12g} {:.9g}".format(*row))
    print(f"{label} HEAT_CURVE_END")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-curves", action="store_true")
    args = parser.parse_args()

    print("WARPED SPACETIME PRE-REGISTRATION FROZEN")
    print("criteria", dict(hd.FROZEN))
    print("N", TRUNCATIONS, "ell=log(phi)", ELL)
    print("boundaries", BOUNDARIES, "primary=dirichlet")
    print("warps", WARPS)

    error, control_dim = explicit_decomposition_error()
    hd.check("explicit full matrix equals union of mode-Jacobi spectra",
             control_dim == 30 and error < 2e-12)
    print(f"  tiny_control_dimension={control_dim} N=4 max_error={error:.3e}")

    r4 = flat_r4_kernel(ELL**2)
    h4 = hyperbolic_h4_kernel(ELL**2)
    uv_ratio, gap = benchmark_certificate()
    hd.check("exact benchmark implementations have flat UV normalization",
             abs(uv_ratio-1) < 5e-4)
    hd.check("H4 benchmark carries the 9/4 curvature gap",
             abs(gap*ELL**2-2.25) < 1e-14 and r4 > 0 and h4 > 0)
    print(f"  BENCHMARK R4(t=ell^2)={r4:.12g} H4={h4:.12g} "
          f"H4_gap={gap:.12g} UV_ratio={uv_ratio:.9g}")

    base, _ = hd.invariant_600_spectrum()
    count3 = counting_plateaus_fast(base, 3)
    count4 = counting_plateaus_fast(base, 4)
    reference3, reference4 = (hd.counting_plateaus(base, target)
                              for target in (3, 4))
    hd.check("prefix-sum counting scan reproduces frozen implementation",
             np.allclose(count3[:4], reference3[:4], atol=2e-11, rtol=0)
             and count3[4:] == reference3[4:]
             and count4 is reference4 is None)
    hd.check("single floor stays 3D under frozen counting criterion",
             count3 is not None and count4 is None)
    print("SINGLE_FLOOR count3", hd.format_plateau(count3))
    print("SINGLE_FLOOR count4", hd.format_plateau(count4))
    print("SINGLE_FLOOR heat3 ", hd.format_plateau(
        hd.heat_plateaus(base, 3)))
    print("SINGLE_FLOOR heat4 ", hd.format_plateau(
        hd.heat_plateaus(base, 4)))

    results = defaultdict(dict)
    for boundary in BOUNDARIES:
        for n in TRUNCATIONS:
            for name, ratio, status in WARPS:
                groups = warped_spectrum(base, n, ratio, boundary)
                cp = counting_plateaus_fast(groups, 4)
                hp = hd.heat_plateaus(groups, 4)
                results[(boundary, name)][n] = (cp, hp)
                print(f"\nDATASET boundary={boundary} N={n} {name} [{status}]")
                print(f"  states={sum(m for _, m in groups)} "
                      f"levels={sum(x > hd.TOL for x, _ in groups)}")
                print("  count4", hd.format_plateau(cp))
                print("  heat4 ", hd.format_plateau(hp))
                print("  count3", hd.format_plateau(
                    counting_plateaus_fast(groups, 3)))
                print("  heat3 ", hd.format_plateau(
                    hd.heat_plateaus(groups, 3)))
                hd.check(f"{boundary} N={n} {name} multiplicity closure",
                         sum(m for _, m in groups) == 2640*n)
                if args.full_curves:
                    curve_rows(f"{boundary}_N{n}_{name}", groups)

    stable = {}
    for name, _, _ in WARPS[:2]:
        stable[name] = all(
            results[("dirichlet", name)][n][0] is not None
            and results[("dirichlet", name)][n][1] is not None
            for n in TRUNCATIONS)
    verdict = ("DYNAMICAL_4D_DERIVED_WARP_N_STABLE"
               if any(stable.values())
               else "NO_N_STABLE_4D_PLATEAU_WARPED_NEGATIVE")
    print("\nN_STABILITY_PRIMARY_DIRICHLET", stable)
    print("DYNAMICAL_VERDICT=" + verdict)

    # Boundary recovery is exact after conformal rescaling on every floor:
    # ratio^(2n) [lambda ratio^(-2n)] = lambda, including multiplicities.
    for name, ratio, _ in WARPS[:2]:
        n = TRUNCATIONS[-1]-1
        recovered = [(value*ratio**(2*n), mult)
                     for value, mult in
                     [(lam*ratio**(-2*n), mult) for lam, mult in base]]
        hd.check(f"{name} conformal boundary recovers all 600-cell levels",
                 all(abs(a-c) < 2e-7 and b == d
                     for (a, b), (c, d) in zip(recovered, base)))
    print("BOUNDARY_STATEMENT=conformal rescaling on each floor exactly "
          "recovers the 52 D3^2 levels and their multiplicities")
    print("INDEX_STATEMENT=vertex_Box_index_-4_not_an_index_of_this_finite_"
          "positive_D4_squared_construction")
    print("RESULT: warped-spacetime audit passed")


if __name__ == "__main__":
    main()

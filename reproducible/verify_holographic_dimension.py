#!/usr/bin/env python3
"""Calibrated spectral-dimension audit.

FROZEN BEFORE SPECTRAL EVALUATION (2026-07-24)
================================================
The constants and control roster in this block are the pre-registration.
They must not be tuned in response to the curves produced below.

Both diagnostics use eigenvalues of D^2.  Thus Weyl scaling is
N(Lambda) ~ Lambda^(d/2), and both reported dimensions include the factor 2:
    d_N = 2 d log N / d log Lambda,
    d_s = -2 d log Tr exp(-t D^2) / d log t.

A target plateau is a contiguous log-scale interval satisfying all of:
  * width >= 0.50 decade;
  * fitted/mean dimension within 0.35 of the integer target;
  * counting: RMSE of log N about its fitted line <= 0.08 and the standard
    deviation of local d_N values <= 0.35;
  * heat: standard deviation of d_s <= 0.35 and
    max |d d_s / d log10(t)| <= 1.00.

Full curves use a frozen 241-point logarithmic grid for heat and one point
per distinct positive eigenvalue for counting.  Heat times cover
[10^-2/lambda_max, 10^2/lambda_min].  Counting local slopes are the slopes
of frozen centered 7-level regressions (shrunk only at endpoints).

Controls:
  S^3: boundaries of the 5-cell, 16-cell, and 24-cell.
  T^3: periodic 4^3 Freudenthal triangulation.
  T^4: periodic 3^4 Freudenthal triangulation (genuinely four-dimensional).

Decision rule:
The 600-cell has a genuine holographic 4D-counting anomaly iff
(i) weighted counting has a registered d=4 plateau, (ii) every S^3 control
has a d=3 counting plateau or no d=4 plateau, and (iii) T^4 has a d=4
plateau in both estimators.  Every other outcome is artifact or explicitly
inconclusive.  Degeneracy-stripped curves are specificity diagnostics and
do not replace the registered weighted decision.
"""

from collections import defaultdict
import argparse
from itertools import combinations, permutations, product
from pathlib import Path
from types import MappingProxyType

import numpy as np
import scipy.linalg
from scipy.spatial import ConvexHull


FROZEN = MappingProxyType({
    "minimum_decades": 0.50,
    "target_tolerance": 0.35,
    "counting_log_rmse_max": 0.08,
    "counting_local_std_max": 0.35,
    "heat_std_max": 0.35,
    "heat_derivative_max": 1.00,
    "heat_grid_points": 241,
    "counting_local_levels": 7,
    "heat_low_over_lambda_max": 1e-2,
    "heat_high_over_lambda_min": 1e2,
})

CONTROL_ROSTER = (
    ("5-cell boundary", "S3"),
    ("16-cell boundary", "S3"),
    ("24-cell boundary", "S3"),
    ("Freudenthal T3 n=4", "S3"),
    ("Freudenthal T4 n=3", "D4"),
)

TOL = 2e-7
HERE = Path(__file__).resolve().parent


def cluster(values, tol=TOL):
    groups = []
    for x in sorted(float(v) for v in values):
        if abs(x) < tol:
            x = 0.0
        if groups and abs(x - groups[-1][0]) < tol:
            mean, n = groups[-1]
            groups[-1] = ((mean*n+x)/(n+1), n+1)
        else:
            groups.append((x, 1))
    return groups


def spectrum_from_facets(facets):
    """Full cochain D^2 spectrum of a finite pure simplicial complex."""
    facets = [tuple(sorted(f)) for f in facets]
    dim = len(facets[0])-1
    cells = []
    vertices = sorted({v for f in facets for v in f})
    cells.append([(v,) for v in vertices])
    for k in range(1, dim+1):
        cells.append(sorted({face for f in facets
                             for face in combinations(f, k+1)}))
    indices = [{c: i for i, c in enumerate(layer)} for layer in cells]
    boundaries = []
    for k in range(dim):
        b = np.zeros((len(cells[k]), len(cells[k+1])))
        for j, simplex in enumerate(cells[k+1]):
            for omit in range(k+2):
                face = simplex[:omit]+simplex[omit+1:]
                b[indices[k][face], j] = (-1)**omit
        boundaries.append(b)
    eigenvalues = []
    for k, layer in enumerate(cells):
        lap = np.zeros((len(layer), len(layer)))
        if k:
            lap += boundaries[k-1].T @ boundaries[k-1]
        if k < dim:
            lap += boundaries[k] @ boundaries[k].T
        eigenvalues.extend(scipy.linalg.eigvalsh(
            lap, driver="evd", check_finite=False))
    return cluster(eigenvalues), tuple(map(len, cells))


def regular_polytope_controls():
    # A regular simplex centered in R^5, represented by its five vertices.
    simplex = np.eye(5)-np.ones((5, 5))/5
    cross = np.vstack((np.eye(4), -np.eye(4)))
    cell24 = []
    for i, j in combinations(range(4), 2):
        for a, b in product((-1, 1), repeat=2):
            v = np.zeros(4)
            v[i], v[j] = a, b
            cell24.append(v)
    answer = {}
    datasets = (("5-cell boundary", simplex),
                ("16-cell boundary", cross),
                ("24-cell boundary", np.array(cell24)))
    for name, points in datasets:
        if name == "5-cell boundary":
            facets = list(combinations(range(5), 4))
        else:
            facets = sorted({tuple(sorted(f))
                             for f in ConvexHull(points).simplices})
        answer[name] = spectrum_from_facets(facets)
    return answer


def canonical_orbit(simplex):
    """Translation-orbit representative and translation of a sorted cell."""
    simplex = tuple(sorted(tuple(map(int, v)) for v in simplex))
    shift = np.array(simplex[0], dtype=int)
    representative = tuple(tuple(np.array(v)-shift) for v in simplex)
    return representative, shift


def freudenthal_orbits(dim):
    """Cell orbits of the standard translation-invariant triangulation."""
    origin = np.zeros(dim, dtype=int)
    tops = set()
    for p in permutations(range(dim)):
        vertices = [origin.copy()]
        current = origin.copy()
        for axis in p:
            current = current.copy()
            current[axis] += 1
            vertices.append(current)
        tops.add(canonical_orbit(vertices)[0])
    cells = [set() for _ in range(dim+1)]
    cells[dim] = tops
    for k in range(dim, 0, -1):
        for simplex in cells[k]:
            for omit in range(k+1):
                face = simplex[:omit]+simplex[omit+1:]
                cells[k-1].add(canonical_orbit(face)[0])
    return [sorted(layer) for layer in cells]


def torus_spectrum(dim, n):
    """Spectrum by Bloch blocks; avoids a 12,150-square T^4 matrix."""
    cells = freudenthal_orbits(dim)
    indices = [{c: i for i, c in enumerate(layer)} for layer in cells]
    all_values = []
    for momentum in product(range(n), repeat=dim):
        angle = 2*np.pi*np.array(momentum)/n
        boundaries = []
        for k in range(dim):
            b = np.zeros((len(cells[k]), len(cells[k+1])),
                         dtype=complex)
            for j, simplex in enumerate(cells[k+1]):
                for omit in range(k+2):
                    face = simplex[:omit]+simplex[omit+1:]
                    representative, shift = canonical_orbit(face)
                    phase = np.exp(1j*np.dot(angle, shift))
                    b[indices[k][representative], j] += (-1)**omit*phase
            boundaries.append(b)
        for k, layer in enumerate(cells):
            lap = np.zeros((len(layer), len(layer)), dtype=complex)
            if k:
                lap += boundaries[k-1].conj().T @ boundaries[k-1]
            if k < dim:
                lap += boundaries[k] @ boundaries[k].conj().T
            all_values.extend(scipy.linalg.eigvalsh(
                lap, driver="evd", check_finite=False))
    f_vector = tuple(n**dim*len(layer) for layer in cells)
    return cluster(all_values), f_vector


def invariant_600_spectrum():
    # Canonical 52-level table emitted and independently checked by
    # verify_invariant_spectrum.py.  Keeping it here avoids repeating its
    # expensive invariant block diagonalization in run_all.py.
    levels = (
        (0, 2), (.145898034, 8), (.381966011, 18), (.527864045, 12),
        (.697224362, 32), (1.074577082, 50), (1.145898034, 32),
        (1.481801301, 72), (1.763932023, 48), (1.933189667, 60),
        (2.208712153, 32), (2.291796068, 8), (2.381966011, 48),
        (2.618033989, 18), (2.807417596, 96), (2.82180592, 72),
        (3, 80), (3.44807057, 50), (3.556327176, 60),
        (3.585786438, 96), (3.697224362, 80), (4, 108),
        (4.302775638, 32), (4.381966011, 60), (4.618033989, 48),
        (4.763932023, 48), (5, 48), (5.527864045, 18), (6, 16),
        (6.236067977, 48), (6.381966011, 48), (6.414213562, 96),
        (6.477352348, 50), (6.618033989, 60), (6.696392779, 72),
        (6.791287847, 32), (6.854101966, 8), (7, 144),
        (7.302775638, 80), (7.443672824, 60), (7.854101966, 32),
        (8, 50), (8.192582404, 96), (8.618033989, 48), (9, 64),
        (9.066810333, 60), (9.236067977, 48), (9.472135955, 12),
        (12, 50), (14, 72), (14.472135955, 18), (15, 32),
        (15.708203932, 8),
    )
    return list(levels), (120, 720, 1200, 600)


def expanded(groups):
    return np.repeat([x for x, m in groups], [m for x, m in groups])


def counting_curve(groups):
    """Curve on distinct levels; N includes all weighted states at the level."""
    levels = np.array([x for x, _ in groups if x > TOL])
    mult = np.array([m for x, m in groups if x > TOL])
    counts = np.cumsum(mult)
    x = np.log10(levels)
    y = np.log(counts)
    local = np.empty(len(x))
    if len(x) < 3:
        return levels, counts, np.zeros(len(x))
    half = FROZEN["counting_local_levels"]//2
    for i in range(len(x)):
        lo, hi = max(0, i-half), min(len(x), i+half+1)
        if hi-lo < 3:
            lo, hi = max(0, hi-3), min(len(x), lo+3)
        local[i] = 2*np.polyfit(np.log(levels[lo:hi]),
                                y[lo:hi], 1)[0]
    return levels, counts, local


def heat_curve(groups):
    values = expanded([(x, m) for x, m in groups if x > TOL])
    zero_modes = sum(m for x, m in groups if x <= TOL)
    t = np.logspace(
        np.log10(FROZEN["heat_low_over_lambda_max"]/values.max()),
        np.log10(FROZEN["heat_high_over_lambda_min"]/values.min()),
        FROZEN["heat_grid_points"])
    weights = np.exp(-np.outer(t, values))
    # Analytic logarithmic derivative, avoiding finite-difference bias.
    ds = 2*t*(weights @ values)/(zero_modes+weights.sum(axis=1))
    return t, ds


def counting_plateaus(groups, target):
    levels, counts, local = counting_curve(groups)
    x, y = np.log10(levels), np.log(counts)
    found = []
    for i in range(len(x)):
        for j in range(i+2, len(x)):
            width = x[j]-x[i]
            if width < FROZEN["minimum_decades"]:
                continue
            slope, intercept = np.polyfit(np.log(levels[i:j+1]),
                                          y[i:j+1], 1)
            dimension = 2*slope
            rmse = np.sqrt(np.mean(
                (y[i:j+1]-(slope*np.log(levels[i:j+1])+intercept))**2))
            local_std = np.std(local[i:j+1])
            if (abs(dimension-target) <= FROZEN["target_tolerance"]
                    and rmse <= FROZEN["counting_log_rmse_max"]
                    and local_std <= FROZEN["counting_local_std_max"]):
                found.append((width, dimension, rmse, local_std, i, j))
    return max(found, default=None, key=lambda z: (z[0], -z[2]))


def heat_plateaus(groups, target):
    t, ds = heat_curve(groups)
    x = np.log10(t)
    derivative = np.gradient(ds, x)
    found = []
    for i in range(len(x)):
        for j in range(i+2, len(x)):
            width = x[j]-x[i]
            if width < FROZEN["minimum_decades"]:
                continue
            mean, std = np.mean(ds[i:j+1]), np.std(ds[i:j+1])
            maxder = np.max(abs(derivative[i:j+1]))
            if (abs(mean-target) <= FROZEN["target_tolerance"]
                    and std <= FROZEN["heat_std_max"]
                    and maxder <= FROZEN["heat_derivative_max"]):
                found.append((width, mean, std, maxder, i, j))
    return max(found, default=None, key=lambda z: (z[0], -z[2]))


def format_plateau(p):
    if p is None:
        return "NONE"
    return "width={:.3f} d={:.4f} residual/std={:.4f} stability={:.4f}".format(
        p[0], p[1], p[2], p[3])


def check(label, condition):
    print(f"[{'PASS' if condition else 'FAIL'}] {label}")
    if not condition:
        raise AssertionError(label)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", action="store_true",
                        help="omit the full machine-readable curve rows")
    args = parser.parse_args()
    print("HOLOGRAPHIC DIMENSION PRE-REGISTRATION FROZEN")
    print(dict(FROZEN))
    print(CONTROL_ROSTER)
    datasets = {}
    datasets["600-cell boundary"] = (invariant_600_spectrum(), "S3")
    datasets.update({name: (data, "S3")
                     for name, data in regular_polytope_controls().items()})
    datasets["Freudenthal T3 n=4"] = (torus_spectrum(3, 4), "S3")
    datasets["Freudenthal T4 n=3"] = (torus_spectrum(4, 3), "D4")

    results = {}
    for name, ((groups, f_vector), kind) in datasets.items():
        positive = [(x, m) for x, m in groups if x > TOL]
        stripped = ([(0.0, 1)] if any(x <= TOL for x, _ in groups) else []) + [
            (x, 1) for x, _ in positive]
        print(f"\nDATASET {name}: kind={kind} f={f_vector} "
              f"states={sum(m for _, m in groups)} "
              f"positive_levels={len(positive)}")
        entry = {}
        for weighting, spectrum in (("weighted", groups),
                                    ("stripped", stripped)):
            print(f"  {weighting}")
            entry[weighting] = {}
            for target in (3, 4):
                cp = counting_plateaus(spectrum, target)
                hp = heat_plateaus(spectrum, target)
                entry[weighting][target] = (cp, hp)
                print(f"    d={target} counting {format_plateau(cp)}")
                print(f"    d={target} heat     {format_plateau(hp)}")
        levels, counts, dcount = counting_curve(positive)
        t, ds = heat_curve(groups)
        if not args.summary:
            print("  COUNTING_CURVE_BEGIN lambda N d_N")
            for row in zip(levels, counts, dcount):
                print("  {:.12g} {} {:.9g}".format(*row))
            print("  COUNTING_CURVE_END")
            print("  HEAT_CURVE_BEGIN t d_s")
            for row in zip(t, ds):
                print("  {:.12g} {:.9g}".format(*row))
            print("  HEAT_CURVE_END")
        print(f"  EXTREMA d_N=[{dcount.min():.6g},{dcount.max():.6g}] "
              f"d_s=[{ds.min():.6g},{ds.max():.6g}]")
        heat_maxima = [(t[i], ds[i]) for i in range(1, len(ds)-1)
                       if ds[i] >= ds[i-1] and ds[i] > ds[i+1]]
        print("  HEAT_LOCAL_MAXIMA " + (
            ", ".join(f"(t={a:.6g},d={b:.6g})" for a, b in heat_maxima)
            if heat_maxima else "NONE"))
        results[name] = entry

    r600 = results["600-cell boundary"]["weighted"]
    s3_ok = all(
        results[name]["weighted"][3][0] is not None
        or results[name]["weighted"][4][0] is None
        for name, kind in CONTROL_ROSTER if kind == "S3")
    t4 = results["Freudenthal T4 n=3"]["weighted"][4]
    anomaly = r600[4][0] is not None and s3_ok and all(x is not None for x in t4)
    if anomaly:
        verdict = "HOLOGRAPHIC_4D_COUNTING_ANOMALY"
    elif r600[4][0] is not None and not all(x is not None for x in t4):
        verdict = "INCONCLUSIVE_CONTROL_TOO_SMALL_OR_CRITERION_UNCALIBRATED"
    else:
        verdict = "ARTIFACT_NO_CALIBRATED_4D_COUNTING_PLATEAU"
    print(f"\nCALIBRATED_VERDICT={verdict}")
    print(f"DECISION_COMPONENTS 600_count4={r600[4][0] is not None} "
          f"S3_specificity={s3_ok} T4_both4={all(x is not None for x in t4)}")

    check("all complexes have the expected total cochain dimension",
          sum(m for _, m in datasets["600-cell boundary"][0][0]) == 2640
          and sum(m for _, m in datasets["Freudenthal T3 n=4"][0][0])
          == sum(datasets["Freudenthal T3 n=4"][0][1])
          and sum(m for _, m in datasets["Freudenthal T4 n=3"][0][0])
          == sum(datasets["Freudenthal T4 n=3"][0][1]))
    check("600-cell gap is phi^-4",
          abs(min(x for x, _ in datasets["600-cell boundary"][0][0]
                  if x > TOL)-(7-3*np.sqrt(5))/2) < 2e-7)
    check("control Betti kernels are S3=2, T3=8, T4=16",
          sum(m for x, m in datasets["5-cell boundary"][0][0] if x == 0) == 2
          and sum(m for x, m in datasets["Freudenthal T3 n=4"][0][0]
                  if x == 0) == 8
          and sum(m for x, m in datasets["Freudenthal T4 n=3"][0][0]
                  if x == 0) == 16)
    print("RESULT: calibrated holographic-dimension audit passed")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Frozen finite tower-product spectral-dimension audit.

FROZEN BEFORE TOWER SPECTRAL EVALUATION (2026-07-24)
====================================================
This verifier imports, without alteration, the plateau criteria and curve
estimators from verify_holographic_dimension.py.  The registered tower data
are:

  truncations N = 8, 12, 16;
  w1: c_n = phi^(-n), registered golden-inflation scale candidate;
  w2: c_n = 2^(-n), registered McKay/PF scale candidate;
  w3: c_n = 1, flat path control only.

D_tower is the N by N real symmetric Jacobi matrix with c_n in entries
(n,n+1) and (n+1,n).  Product eigenvalues are formed only as pair sums:
spec(D4^2) = {lambda_3 + lambda_t}, with multiplicities multiplied.

A dynamical-4D claim requires a qualifying 4D plateau for w1 or w2 at every
registered N, while the single-floor spectrum retains its qualifying 3D
counting plateau.  A plateau only for w3 means that dimension comes from the
flat product, not the derived scale weights.  Otherwise the verdict is a
finite-size/negative result.  No threshold below is fitted to the output.
"""

import argparse
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify_holographic_dimension as hd


TRUNCATIONS = (8, 12, 16)
PHI = (1 + np.sqrt(5.0)) / 2
WEIGHTS = (
    ("w1_golden", PHI, "REGISTERED SCALE CANDIDATE"),
    ("w2_mckay_pf", 2.0, "REGISTERED SCALE CANDIDATE"),
    ("w3_uniform", 1.0, "FLAT CONTROL"),
)
ANCHOR_EXPONENTS = (5, 6, 25, 35)
RATIO_TOL = 1e-10


def tower_spectrum(n, ratio):
    d = np.zeros((n, n))
    for level in range(n - 1):
        d[level, level + 1] = d[level + 1, level] = ratio ** (-level)
    return hd.cluster(np.linalg.eigvalsh(d) ** 2, tol=1e-9), d


def product_spectrum(left, right):
    """Cluster pair sums without constructing a tensor-product matrix."""
    raw = [(a + b, ma * mb) for a, ma in left for b, mb in right]
    raw.sort()
    out = []
    for value, multiplicity in raw:
        if abs(value) < hd.TOL:
            value = 0.0
        if out and abs(value - out[-1][0]) < hd.TOL:
            old, mult = out[-1]
            out[-1] = ((old * mult + value * multiplicity)
                       / (mult + multiplicity), mult + multiplicity)
        else:
            out.append((value, multiplicity))
    return out


def exact_factorization_certificate():
    """Exact odd-product identity over Q(sqrt(5)) on a finite witness."""
    q = sp.Matrix([[1, -1, 2], [0, 1, 1]])
    z22, z33 = sp.zeros(2), sp.zeros(3)
    d3 = z22.row_join(q)
    d3 = d3.col_join(q.T.row_join(z33))
    gamma = sp.diag(1, 1, -1, -1, -1)
    c0, c1 = sp.symbols("c0 c1", real=True)
    dt = sp.Matrix([[0, c0, 0], [c0, 0, c1], [0, c1, 0]])
    d4 = sp.kronecker_product(d3, sp.eye(3)) + sp.kronecker_product(gamma, dt)
    rhs = (sp.kronecker_product(d3 * d3, sp.eye(3))
           + sp.kronecker_product(sp.eye(5), dt * dt))
    return gamma * d3 + d3 * gamma == sp.zeros(5), sp.simplify(d4*d4-rhs) == sp.zeros(15)


def anchor_hits(groups):
    positive = sorted(x for x, _ in groups if x > hd.TOL)
    hits = []
    comparisons = 0
    for i, a in enumerate(positive):
        for b in positive[i + 1:]:
            comparisons += 1
            exponent = np.log(b / a) / np.log(PHI)
            for target in ANCHOR_EXPONENTS:
                if abs(exponent - target) < RATIO_TOL:
                    hits.append((a, b, target))
    return comparisons, hits


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
    parser.add_argument("--full-curves", action="store_true",
                        help="print every counting and heat-curve row")
    args = parser.parse_args()

    print("TOWER SPACETIME PRE-REGISTRATION FROZEN")
    print("criteria", dict(hd.FROZEN))
    print("truncations", TRUNCATIONS)
    print("weights", [(name, ratio, status) for name, ratio, status in WEIGHTS])
    print("anchors", ANCHOR_EXPONENTS, "ratio tolerance", RATIO_TOL)

    anti, factor = exact_factorization_certificate()
    hd.check("exact form-parity anticommutation certificate", anti)
    hd.check("exact D4^2 tensor-sum factorization certificate", factor)

    base, _ = hd.invariant_600_spectrum()
    base_c3 = hd.counting_plateaus(base, 3)
    base_c4 = hd.counting_plateaus(base, 4)
    print("\nSINGLE_FLOOR")
    print("  count3", hd.format_plateau(base_c3))
    print("  count4", hd.format_plateau(base_c4))
    print("  heat3 ", hd.format_plateau(hd.heat_plateaus(base, 3)))
    print("  heat4 ", hd.format_plateau(hd.heat_plateaus(base, 4)))
    hd.check("single floor retains frozen 3D counting verdict",
             base_c3 is not None and base_c4 is None)

    results = defaultdict(dict)
    anchor_total = 0
    anchor_hits_total = []
    for n in TRUNCATIONS:
        for name, ratio, status in WEIGHTS:
            tower, matrix = tower_spectrum(n, ratio)
            product = product_spectrum(base, tower)
            tc1, th1 = hd.counting_plateaus(tower, 1), hd.heat_plateaus(tower, 1)
            pc4, ph4 = hd.counting_plateaus(product, 4), hd.heat_plateaus(product, 4)
            results[name][n] = (pc4, ph4)
            comparisons, hits = anchor_hits(tower)
            anchor_total += comparisons
            anchor_hits_total.extend((name, n, *hit) for hit in hits)
            print(f"\nDATASET N={n} {name} [{status}]")
            print(f"  tower_states={n} tower_positive_levels="
                  f"{sum(x > hd.TOL for x, _ in tower)} "
                  f"product_states={sum(m for _, m in product)} "
                  f"product_levels={sum(x > hd.TOL for x, _ in product)}")
            print("  tower count1", hd.format_plateau(tc1))
            print("  tower heat1 ", hd.format_plateau(th1))
            print("  product count4", hd.format_plateau(pc4))
            print("  product heat4 ", hd.format_plateau(ph4))
            print("  product count3",
                  hd.format_plateau(hd.counting_plateaus(product, 3)))
            print("  product heat3 ",
                  hd.format_plateau(hd.heat_plateaus(product, 3)))
            print(f"  anchor spectral-ratio comparisons={comparisons} hits={len(hits)}")
            hd.check(f"N={n} {name} tower trace D^2 equals 2 sum c_n^2",
                     abs(np.trace(matrix @ matrix)
                         - 2*sum(ratio**(-2*k) for k in range(n-1))) < 1e-10)
            hd.check(f"N={n} {name} product multiplicity closes",
                     sum(m for _, m in product) == 2640*n)
            if args.full_curves:
                curve_rows(f"N={n}_{name}_TOWER", tower)
                curve_rows(f"N={n}_{name}_PRODUCT", product)

    derived_stable = {
        name: all(results[name][n][0] is not None
                  and results[name][n][1] is not None for n in TRUNCATIONS)
        for name, _, _ in WEIGHTS[:2]
    }
    flat_stable = all(results["w3_uniform"][n][0] is not None
                      and results["w3_uniform"][n][1] is not None
                      for n in TRUNCATIONS)
    if any(derived_stable.values()):
        verdict = "DYNAMICAL_4D_DERIVED_WEIGHT_N_STABLE"
    elif flat_stable:
        verdict = "DIMENSION_FROM_FLAT_PRODUCT_NOT_DERIVED_WEIGHTS"
    else:
        verdict = "NO_N_STABLE_4D_PLATEAU_FINITE_SIZE_NEGATIVE"
    print("\nANCHOR_SEARCH")
    print(f"  registered exponents={ANCHOR_EXPONENTS}")
    print(f"  spectral-ratio comparisons={anchor_total}")
    print(f"  exact-tolerance hits={anchor_hits_total or 'NONE'}")
    print("  w1 hop ratios at separation k are tautologically phi^-k; "
          "k=5,6 occur when the truncation permits, while 25,35 do not. "
          "No spacing is distinguished from the other integer separations.")
    print(f"\nDYNAMICAL_VERDICT={verdict}")
    print("N_STABILITY derived", derived_stable, "flat", flat_stable)
    hd.check("registered anchor spectral search is negative",
             not anchor_hits_total)
    print("RESULT: tower-spacetime audit passed")


if __name__ == "__main__":
    main()

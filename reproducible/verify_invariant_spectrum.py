#!/usr/bin/env python3
"""Invariant Kaehler--Dirac spectrum audit for the full 600-cell complex.

ANTI-NUMEROLOGY FREEZE
----------------------
The registry below was written before computing the spectrum.  Only these
targets may be called registered matches.  Anything noticed after spectral
computation is explicitly an unregistered pattern and carries a
multiple-comparisons caveat.
"""

from collections import Counter, defaultdict
from itertools import permutations, product
from types import MappingProxyType

import numpy as np
import scipy.linalg


# Frozen pre-registration: do not extend from computed spectral output.
REGISTERED_TARGETS = MappingProxyType({
    "phi_definition": "(1 + sqrt(5))/2",
    "phi_powers": "phi**n for integer n",
    # Nine entries are retained; exponent 11 occurs twice in the assignment.
    "mass_exponents_n_5a_plus_6b": (0, 5, 3, 11, 16, 19, 26, 17, 11),
    "lambda_1": "1/(2*phi**2)",
    "alpha_inverse_observed": "137.036",
    "alpha_inverse_tree": "20*phi**4",
    "alpha_s": "1/(2*phi**3)",
    "sin2_theta_w": "6/26",
    "a_1": 5,
    "b_1": 6,
    "N": 120,
    "spectral_action_coefficients": (2640, 14880, 55920),
    "spectral_action_reduced": (11, 62, 233),
    "neutrino_conformal_s": "3/4",
    "cp_angle": "atan(sqrt(5))",
})


def permutation_sign(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p))
                     for j in range(i + 1, len(p))) % 2 else 1


def cluster(values, tol=2e-7):
    out = []
    for value in sorted(map(float, values)):
        if out and abs(value - out[-1][0]) < tol:
            old, count = out[-1]
            out[-1] = ((old * count + value) / (count + 1), count + 1)
        else:
            out.append((value, 1))
    return out


def build_complex():
    phi = (1 + np.sqrt(5.0)) / 2
    vertices = set()
    for i in range(4):
        for sign in (-1.0, 1.0):
            q = [0.0] * 4
            q[i] = sign
            vertices.add(tuple(q))
    vertices.update(product((-0.5, 0.5), repeat=4))
    base = [phi / 2, 0.5, 1 / (2 * phi), 0.0]
    even = [p for p in permutations(range(4)) if permutation_sign(p) == 1]
    for p in even:
        q = [base[p[i]] for i in range(4)]
        nonzero = [i for i, x in enumerate(q) if abs(x) > 1e-12]
        for signs in product((-1, 1), repeat=3):
            r = q[:]
            for i, sign in zip(nonzero, signs):
                r[i] *= sign
            vertices.add(tuple(round(x, 10) for x in r))
    vertices = np.array(sorted(vertices))
    dots = vertices @ vertices.T
    edges = [(i, j) for i in range(120) for j in range(i + 1, 120)
             if abs(dots[i, j] - phi / 2) < 1e-3]
    adj = defaultdict(set)
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    triangles = [(i, j, k) for i, j in edges for k in adj[i] & adj[j]
                 if j < k]
    tetrahedra = [(i, j, k, ell) for i, j, k in triangles
                  for ell in adj[i] & adj[j] & adj[k] if k < ell]
    cells = [[(i,) for i in range(120)], edges, triangles, tetrahedra]
    indices = [{cell: i for i, cell in enumerate(layer)} for layer in cells]
    d = []
    for degree in range(3):
        matrix = np.zeros((len(cells[degree + 1]), len(cells[degree])))
        for row, simplex in enumerate(cells[degree + 1]):
            for omit in range(degree + 2):
                face = simplex[:omit] + simplex[omit + 1:]
                matrix[row, indices[degree][face]] = (-1) ** omit
        d.append(matrix)
    delta = [
        d[0].T @ d[0],
        d[0] @ d[0].T + d[1].T @ d[1],
        d[1] @ d[1].T + d[2].T @ d[2],
        d[2] @ d[2].T,
    ]
    return phi, vertices, cells, indices, delta


def qmul(a, b):
    w, x, y, z = a
    W, X, Y, Z = b
    return np.array((w*W-x*X-y*Y-z*Z, w*X+x*W+y*Z-z*Y,
                     w*Y-x*Z+y*W+z*X, w*Z+x*Y-y*X+z*W))


def invariant_spectrum():
    phi, vertices, cells, indices, delta = build_complex()
    sqrt5 = np.sqrt(5.0)
    class_x = np.array([2, -2, 0, 1, -1, phi, -phi,
                        phi - 1, (1 - sqrt5) / 2])
    class_sizes = np.array([1, 1, 30, 20, 20, 12, 12, 12, 12])
    xp = np.array([2, -2, 0, 1, -1, (1-sqrt5)/2,
                   phi-1, -phi, phi])

    def sym_powers(t):
        answer = [np.ones(9), t]
        for _ in range(2, 6):
            answer.append(t * answer[-1] - answer[-2])
        return answer

    sx, sxp = sym_powers(class_x), sym_powers(xp)
    standard = [sx[0], class_x, xp, sx[2], sxp[2],
                class_x*xp, sx[3], sx[4], sx[5]]
    order = (0, 1, 3, 6, 7, 8, 5, 2, 4)
    chars = np.array([standard[i] for i in order])
    names = ("rho0", "rho1", "rho2", "rho3", "rho4",
             "rho5", "rho6", "rho7", "rho8")
    dims_ir = np.rint(chars[:, 0]).astype(int)

    weights = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256],
                       dtype=float)
    central_by_degree = [
        np.zeros((len(layer), len(layer))) for layer in cells
    ]
    for g in vertices:
        moved = np.array([qmul(g, v) for v in vertices])
        distances = ((moved[:, None] - vertices[None]) ** 2).sum(axis=2)
        p = distances.argmin(axis=1)
        conjugacy_class = int(np.argmin(abs(class_x - 2*g[0])))
        for degree, layer in enumerate(cells):
            for source, cell in enumerate(layer):
                image = [int(p[i]) for i in cell]
                ordering = sorted(range(len(image)), key=image.__getitem__)
                target = indices[degree][tuple(sorted(image))]
                central_by_degree[degree][target, source] += (
                    weights[conjugacy_class] * permutation_sign(ordering)
                )

    # A generic real central class sum. Its nine scalars distinguish all
    # irreps, and it never selects bases inside their multiplicity spaces.
    central_scalars = (
        chars @ (class_sizes * weights) / dims_ir
    )
    assert min(abs(a-b) for i, a in enumerate(central_scalars)
               for b in central_scalars[i+1:]) > 1e-5
    epsilon = 1e-4
    rows = []
    for degree in range(4):
        central = central_by_degree[degree]
        central = (central + central.T) / 2
        raw_delta = scipy.linalg.eigvalsh(delta[degree],
                                          driver="evd",
                                          check_finite=False)
        raw_joint = scipy.linalg.eigvalsh(delta[degree] + epsilon * central,
                                          driver="evd",
                                          check_finite=False)
        delta_clusters = cluster(raw_delta)
        joint_clusters = cluster(raw_joint)
        candidates = []
        for lam, _ in delta_clusters:
            for irrep, scalar in enumerate(central_scalars):
                candidates.append((lam + epsilon * scalar, lam, irrep))
        for value, multiplicity in joint_clusters:
            _, lam, irrep = min(candidates, key=lambda item:
                                abs(item[0] - value))
            residual = abs(value - (lam + epsilon*central_scalars[irrep]))
            assert residual < 2e-6, (degree, value, residual)
            rows.append((degree, lam, irrep, multiplicity))
    return phi, names, dims_ir, chars, rows


def main():
    print("INVARIANT SPECTRUM PRE-REGISTRATION FROZEN")
    print(f"registered target families: {len(REGISTERED_TARGETS)}")
    phi, names, dims_ir, chars, rows = invariant_spectrum()
    print("SPECTRUM_ROWS_BEGIN")
    for degree, lam, irrep, multiplicity in sorted(
            rows, key=lambda x: (x[1], x[2], x[0])):
        print(f"{lam:.12g} {multiplicity:4d} {names[irrep]:4s} C^{degree}")
    print("SPECTRUM_ROWS_END")
    total = sum(row[3] for row in rows)
    zeros = sum(row[3] for row in rows if abs(row[1]) < 1e-7)
    assert total == 2640
    assert zeros == 2
    positive = [(lam, mult) for _, lam, _, mult in rows if lam > 1e-7]
    distinct = sorted(set(round(lam, 10) for lam, _ in positive))
    # Independent exact finite-moment path.  These Laplacian blocks have
    # integral entries.  For symmetric Delta, Tr(Delta^2) is the sum of
    # squared entries, so this does not reuse numerical eigenvalues.
    _, _, _, _, exact_source = build_complex()
    delta_integer = [np.rint(block).astype(np.int64) for block in exact_source]
    exact_moment1 = sum(int(np.trace(block)) for block in delta_integer)
    exact_moment2_half = sum(int(np.sum(block * block))
                             for block in delta_integer) // 2
    assert exact_moment1 == 14880
    assert exact_moment2_half == 55920

    moment1 = sum(lam * mult for lam, mult in positive)
    moment2_half = sum(lam**2 * mult for lam, mult in positive) / 2
    print(f"DISTINCT_POSITIVE={len(distinct)}")
    print(f"TRACE_D2={moment1:.12f}")
    print(f"HALF_TRACE_D4={moment2_half:.12f}")
    print(f"EXACT_TRACE_D2={exact_moment1}")
    print(f"EXACT_HALF_TRACE_D4={exact_moment2_half}")
    for s in (0.75, 1.0, 2.0):
        value = sum(mult * lam**(-s) for lam, mult in positive)
        print(f"ZETA_D2({s})={value:.15g}")
    for t in (5.0, phi):
        value = zeros + sum(mult * np.exp(-t*lam)
                            for lam, mult in positive)
        print(f"THETA({t:.15g})={value:.15g}")

    # Pre-registered golden-ratio confrontation.  The comparison universe is
    # all unordered ratios of distinct positive eigenvalues.  A half-integer
    # exponent hit uses a deliberately severe 1e-10 absolute log threshold.
    exponent_hits = []
    ratio_comparisons = 0
    for i, numerator in enumerate(distinct):
        for denominator in distinct[:i]:
            ratio_comparisons += 1
            exponent = np.log(numerator / denominator) / np.log(phi)
            nearest_half = round(2*exponent) / 2
            if abs(exponent-nearest_half) < 1e-10:
                exponent_hits.append((numerator, denominator, nearest_half))
    registered_exponents = set(REGISTERED_TARGETS[
        "mass_exponents_n_5a_plus_6b"])
    found_registered = sorted({
        int(exponent) for _, _, exponent in exponent_hits
        if exponent == int(exponent) and int(exponent) in registered_exponents
    })
    gap = min(distinct)
    target_gap = 1 / (2*phi**2)
    print(f"RATIO_COMPARISONS={ratio_comparisons}")
    print(f"HALF_INTEGER_HITS={len(exponent_hits)}")
    for hit in exponent_hits:
        print("PHI_RATIO_HIT {:.12g}/{:.12g}=phi^{}".format(*hit))
    print(f"REGISTERED_MASS_EXPONENTS_FOUND={found_registered}")
    print(f"GAP={gap:.15g}")
    print(f"REGISTERED_LAMBDA1={target_gap:.15g}")
    print(f"GAP_OVER_LAMBDA1={gap/target_gap:.15g}")

    quaternionic = {1, 3, 5, 7}
    kramers_ok = all(mult % 2 == 0 for _, _, irrep, mult in rows
                     if irrep in quaternionic)
    assert kramers_ok
    assert abs(moment1 - 14880) < 2e-8
    assert abs(moment2_half - 55920) < 2e-7
    print("KRAMERS_EVEN_DEGENERACIES=True")
    print("RESULT: all invariant spectrum checks passed")


if __name__ == "__main__":
    main()

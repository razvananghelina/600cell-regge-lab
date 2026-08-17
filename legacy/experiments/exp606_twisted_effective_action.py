"""
exp606: Non-circular one-loop-style effective-action test on the twisted 600-cell.

Purpose
-------
exp605 was circular as a derivation of the alpha equation: it built a
functional whose Euler-Lagrange equation was already known in advance.
This script instead studies only quantities computed directly from the
gauge-twisted simplicial complex.

We introduce the canonical twisted operators

    d0(theta)         : twisted incidence on edges
    Delta0(theta)     = d0(theta)^* d0(theta)          on vertices
    B(theta)          = d0(theta) d0(theta)^*          exact 1-form sector
    Delta1(theta)     = B(theta) + C                   on edges
    Box1(theta)       = C - a1 B(theta)                signed 1-form contrast

where C = d1^* d1 is the coexact 1-form operator and the twist is the Hopf
fiber U(1) phase used in exp570/571.

Main question
-------------
Do natural spectral functionals built from these twisted operators produce a
non-circular alpha-like coefficient from their small-theta curvature?

We test the canonical logarithmic determinants

    W_O(theta) = (1/2) log det' O(theta)

for O in {Delta0, B, Delta1, |Box1|}, together with the gauge-fixed
combination

    W_gf(theta) = (1/2) log det' Delta1(theta) - log det' Delta0(theta).

For a one-loop-style effective action

    Gamma(theta) = S_YM(theta)/(4*pi*alpha) + W(theta),

the quadratic coefficient at theta = 0 would imply the candidate relation

    alpha^{-1}_cand = -2*pi * W''(0) / K_YM,

where S_YM(theta) = K_YM theta^2 + O(theta^4) with K_YM = 6.

This is exploratory but non-circular: W(theta) is computed from the twisted
operator family, not designed to reproduce the known alpha equation.
"""

from __future__ import annotations

from collections import defaultdict
import math
import sys

import numpy as np

sys.path.insert(0, ".")
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0
A1 = 5
B1 = 6
N = 120
TOL = 1e-9

ALPHA_INV_TREE = 4.0 * A1 * PHI**4
ALPHA_INV_FULL = 137.036188770018


def qmul(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.array(
        [
            p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
            p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
            p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
            p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0],
        ]
    )


def find_idx(v: np.ndarray, verts: np.ndarray, tol: float = 1e-6) -> int:
    dots = verts @ v
    idx = int(np.argmax(dots))
    return idx if dots[idx] > 1.0 - tol else -1


def find_fibration(verts: np.ndarray) -> list[list[int]]:
    target_w = PHI / 2.0
    for i in range(N):
        if abs(verts[i, 0] - target_w) >= 1e-6:
            continue

        g = verts[i]
        p = g.copy()
        ok = True
        for k in range(2, 11):
            p = qmul(p, g)
            if k == 5 and not np.allclose(p, [-1, 0, 0, 0], atol=1e-6):
                ok = False
                break
            if k == 10 and not np.allclose(p, [1, 0, 0, 0], atol=1e-6):
                ok = False
        if not ok:
            continue

        used: set[int] = set()
        fibers: list[list[int]] = []
        subgroup: list[int] = []
        pp = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_idx(pp, verts))
            pp = qmul(pp, g)

        for s in range(N):
            if s in used:
                continue
            fib: list[int] = []
            for si in subgroup:
                q = qmul(verts[s], verts[si])
                idx = find_idx(q, verts)
                if idx >= 0 and idx not in used:
                    fib.append(idx)
                    used.add(idx)
            if len(fib) == 10:
                fibers.append(fib)

        if len(fibers) == 12:
            return fibers

    raise RuntimeError("Could not find Hopf fibration")


def get_ordered_cycles(fibers: list[list[int]], adj: np.ndarray) -> list[list[int]]:
    cycles: list[list[int]] = []
    for fib in fibers:
        fib_set = set(fib)
        cycle = [fib[0]]
        visited = {fib[0]}
        while len(cycle) < 10:
            curr = cycle[-1]
            found = False
            for j in fib_set - visited:
                if adj[curr, j] > 0.5:
                    cycle.append(j)
                    visited.add(j)
                    found = True
                    break
            if not found:
                raise RuntimeError("Could not order Hopf fiber as a 10-cycle")
        cycles.append(cycle)
    return cycles


def build_complex(adj: np.ndarray):
    adj_list = defaultdict(set)
    edges: list[tuple[int, int]] = []
    edge_to_idx: dict[tuple[int, int], int] = {}

    for i in range(N):
        for j in range(i + 1, N):
            if adj[i, j] > 0.5:
                adj_list[i].add(j)
                adj_list[j].add(i)
                edge_to_idx[(i, j)] = len(edges)
                edges.append((i, j))

    triangles: list[tuple[int, int, int]] = []
    for i in range(N):
        for j in adj_list[i]:
            if j <= i:
                continue
            for k in adj_list[i] & adj_list[j]:
                if k > j:
                    triangles.append((i, j, k))

    return edges, edge_to_idx, triangles


def build_d1(
    triangles: list[tuple[int, int, int]],
    edge_to_idx: dict[tuple[int, int], int],
    n_edges: int,
) -> np.ndarray:
    d1 = np.zeros((len(triangles), n_edges), dtype=complex)
    for f_idx, (i, j, k) in enumerate(triangles):
        d1[f_idx, edge_to_idx[(i, j)]] = +1.0
        d1[f_idx, edge_to_idx[(j, k)]] = +1.0
        d1[f_idx, edge_to_idx[(i, k)]] = -1.0
    return d1


def build_fiber_forward(cycles: list[list[int]]) -> set[tuple[int, int]]:
    forward: set[tuple[int, int]] = set()
    for cyc in cycles:
        for k in range(10):
            forward.add((cyc[k], cyc[(k + 1) % 10]))
    return forward


def edge_phase(i: int, j: int, theta: float, fiber_forward: set[tuple[int, int]]) -> complex:
    if (i, j) in fiber_forward:
        return np.exp(1j * theta / 10.0)
    if (j, i) in fiber_forward:
        return np.exp(-1j * theta / 10.0)
    return 1.0 + 0.0j


def build_d0(
    theta: float,
    edges: list[tuple[int, int]],
    fiber_forward: set[tuple[int, int]],
) -> np.ndarray:
    d0 = np.zeros((len(edges), N), dtype=complex)
    for e_idx, (i, j) in enumerate(edges):
        d0[e_idx, i] = -edge_phase(i, j, theta, fiber_forward)
        d0[e_idx, j] = +1.0
    return d0


def sym_eigs(mat: np.ndarray) -> np.ndarray:
    herm = 0.5 * (mat + mat.conj().T)
    return np.linalg.eigvalsh(herm)


def logdet_prime_positive(mat: np.ndarray) -> float:
    eigs = sym_eigs(mat)
    nz = eigs[eigs > TOL]
    return float(np.sum(np.log(nz)))


def logdet_prime_abs(mat: np.ndarray) -> float:
    eigs = sym_eigs(mat)
    nz = eigs[np.abs(eigs) > TOL]
    return float(np.sum(np.log(np.abs(nz))))


def kernel_dim(mat: np.ndarray, use_abs: bool = False) -> int:
    eigs = sym_eigs(mat)
    if use_abs:
        return int(np.sum(np.abs(eigs) <= TOL))
    return int(np.sum(np.abs(eigs) <= TOL))


def small_theta_quadratic(
    theta_vals: np.ndarray, values: np.ndarray
) -> tuple[float, float, float]:
    x = theta_vals**2
    y = values - values[0]
    coeffs = np.polyfit(x, y, 2)
    c4, c2, _ = coeffs
    second = 2.0 * c2
    fit_vals = coeffs[0] * x**2 + coeffs[1] * x + coeffs[2]
    residual = float(np.max(np.abs(y - fit_vals)))
    return float(c2), float(second), residual


print("=" * 78)
print("EXP606: TWISTED EFFECTIVE-ACTION TEST")
print("=" * 78)

verts, adj, _ = build_600cell()
fibers = find_fibration(verts)
cycles = get_ordered_cycles(fibers, adj)
fiber_forward = build_fiber_forward(cycles)
edges, edge_to_idx, triangles = build_complex(adj)
d1 = build_d1(triangles, edge_to_idx, len(edges))
C = d1.conj().T @ d1

print("\nGeometry")
print(f"  Vertices:  {len(verts)}")
print(f"  Edges:     {len(edges)}")
print(f"  Triangles: {len(triangles)}")
print(f"  Fibers:    {len(fibers)} x 10")

theta_vals = np.array([0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.16])
sym_theta_vals = theta_vals

sym_data: dict[str, list[float]] = {
    "sym": [],
    "ym": [],
    "delta0_logdet": [],
    "b_logdet": [],
    "delta1_logdet": [],
    "box1_logdet_abs": [],
    "w_delta0": [],
    "w_b": [],
    "w_delta1": [],
    "w_box1": [],
    "w_gf": [],
}

for theta in sym_theta_vals:
    d0 = build_d0(theta, edges, fiber_forward)
    delta0 = d0.conj().T @ d0
    b_theta = d0 @ d0.conj().T
    delta1 = b_theta + C
    box1 = C - A1 * b_theta
    F = d1 @ d0
    s_ym = float(np.real(np.trace(F.conj().T @ F)))

    logdet_delta0 = logdet_prime_positive(delta0)
    logdet_b = logdet_prime_positive(b_theta)
    logdet_delta1 = logdet_prime_positive(delta1)
    logdet_box1 = logdet_prime_abs(box1)

    sym_data["sym"].append(theta)
    sym_data["ym"].append(s_ym)
    sym_data["delta0_logdet"].append(logdet_delta0)
    sym_data["b_logdet"].append(logdet_b)
    sym_data["delta1_logdet"].append(logdet_delta1)
    sym_data["box1_logdet_abs"].append(logdet_box1)
    sym_data["w_delta0"].append(0.5 * logdet_delta0)
    sym_data["w_b"].append(0.5 * logdet_b)
    sym_data["w_delta1"].append(0.5 * logdet_delta1)
    sym_data["w_box1"].append(0.5 * logdet_box1)
    sym_data["w_gf"].append(0.5 * logdet_delta1 - logdet_delta0)

ym_vals = np.array(sym_data["ym"])
ym_c2, ym_second, ym_resid = small_theta_quadratic(theta_vals, ym_vals)
print("\nYang-Mills control")
print(f"  K_YM from fit = {ym_c2:.12f}")
print(f"  S_YM''(0)     = {ym_second:.12f}")
print(f"  expected K_YM = {B1}")
print(f"  fit residual  = {ym_resid:.3e}")

theta_probe = 0.01
d0_0 = build_d0(0.0, edges, fiber_forward)
d0_probe = build_d0(theta_probe, edges, fiber_forward)
B_0 = d0_0 @ d0_0.conj().T
B_probe = d0_probe @ d0_probe.conj().T
Delta0_0 = d0_0.conj().T @ d0_0
Delta0_probe = d0_probe.conj().T @ d0_probe
Delta1_0 = B_0 + C
Delta1_probe = B_probe + C
Box1_0 = C - A1 * B_0
Box1_probe = C - A1 * B_probe

print("\nAnalyticity diagnostic at theta -> 0")
diag_rows = [
    ("Delta0", Delta0_0, Delta0_probe, False),
    ("B", B_0, B_probe, False),
    ("Delta1", Delta1_0, Delta1_probe, False),
    ("|Box1|", Box1_0, Box1_probe, True),
]
for label, mat0, matp, use_abs in diag_rows:
    eigs_p = sym_eigs(matp)
    smallest = (
        float(np.min(np.abs(eigs_p[np.abs(eigs_p) > TOL])))
        if use_abs
        else float(np.min(eigs_p[eigs_p > TOL]))
    )
    print(
        f"  {label:>6s}: ker(0)={kernel_dim(mat0, use_abs):4d}, "
        f"ker({theta_probe:g})={kernel_dim(matp, use_abs):4d}, "
        f"smallest nonzero @ {theta_probe:g} = {smallest:.6e}"
    )
print("  If the kernel changes under an infinitesimal twist, log det' picks up a")
print("  non-analytic log(theta) term and its raw quadratic curvature is not a")
print("  trustworthy effective-action coefficient.")

functionals = [
    ("W_Delta0", "w_delta0"),
    ("W_B", "w_b"),
    ("W_Delta1", "w_delta1"),
    ("W_|Box1|", "w_box1"),
    ("W_gf", "w_gf"),
]

print("\nSmall-theta curvature of computed functionals")
print(
    f"{'functional':>12s} {'status':>12s} {'W\"(0)':>16s} {'alpha^-1_cand':>18s}"
)

results: list[tuple[str, float, float]] = []
for label, key in functionals:
    vals = np.array(sym_data[key])
    _, second, resid = small_theta_quadratic(theta_vals, vals)
    if label in {"W_Delta0", "W_B", "W_gf"}:
        status = "nonanalytic"
        print(f"{label:>12s} {status:>12s} {second:16.9f} {'N/A':>18s}")
        print(f"  raw fit residual for {label}: {resid:.3e}")
        continue

    alpha_inv_cand = -2.0 * math.pi * second / B1
    status = "analytic"
    results.append((label, second, alpha_inv_cand))
    print(f"{label:>12s} {status:>12s} {second:16.9f} {alpha_inv_cand:18.9f}")
    print(
        f"  ratios for {label}: "
        f"cand/tree={alpha_inv_cand / ALPHA_INV_TREE:.6f}, "
        f"cand/full={alpha_inv_cand / ALPHA_INV_FULL:.6f}, "
        f"fit residual={resid:.3e}"
    )

print("\nSample values at distinguished angles")
angles = [
    ("0", 0.0),
    ("pi/5", math.pi / 5.0),
    ("2pi", 2.0 * math.pi),
]

angle_data: dict[float, dict[str, float]] = {}
for _, theta in angles:
    d0 = build_d0(theta, edges, fiber_forward)
    delta0 = d0.conj().T @ d0
    b_theta = d0 @ d0.conj().T
    delta1 = b_theta + C
    box1 = C - A1 * b_theta
    F = d1 @ d0
    angle_data[theta] = {
        "S_YM": float(np.real(np.trace(F.conj().T @ F))),
        "logdet_Delta0": logdet_prime_positive(delta0),
        "logdet_B": logdet_prime_positive(b_theta),
        "logdet_Delta1": logdet_prime_positive(delta1),
        "logdet_Box1_abs": logdet_prime_abs(box1),
    }

for label, theta in angles:
    vals = angle_data[theta]
    print(
        f"  theta={label:>4s}: "
        f"S_YM={vals['S_YM']:12.6f}, "
        f"logdet Delta0={vals['logdet_Delta0']:12.6f}, "
        f"logdet B={vals['logdet_B']:12.6f}, "
        f"logdet Delta1={vals['logdet_Delta1']:12.6f}, "
        f"logdet |Box1|={vals['logdet_Box1_abs']:12.6f}"
    )

print("\nInterpretation")
best = min(results, key=lambda item: abs(item[2] - ALPHA_INV_FULL))
print(
    f"  Closest analytic one-loop-style candidate among the tested canonical "
    f"functionals: "
    f"{best[0]} with alpha^-1_cand = {best[2]:.9f}"
)
print("  The mapping alpha^-1_cand = -2*pi*W''(0)/K_YM is the standard")
print("  curvature matching for Gamma(theta)=S_YM/(4*pi*alpha)+W(theta).")
print("  Here the only analytic determinant candidates are tiny compared with 137,")
print("  while the formally closer Delta0/B values are invalid because the twist")
print("  lifts zero modes and induces non-analytic log(theta) behavior.")
print("  So the naive determinant-based one-loop closure is insufficient and the")
print("  KK/topological bridge remains genuinely open.")

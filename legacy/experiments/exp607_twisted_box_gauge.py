"""
exp607: Twisted determinant directly on the gauge-projected Box sector.

Motivation
----------
exp606 showed that determinant functionals built from the full simplicial
Laplacians do not recover the fine-structure coefficient. That is not too
surprising: alpha is expected to live in the gauge-projected Box sector, not
in the raw geometry.

So this script targets exactly that sector. We build

    Box0(theta) = L_cross - a1 * L_fiber(theta)

where only the Hopf-fiber edges are twisted, and then project onto the
12-dimensional fiber-constant sector:

    Box_gauge(theta) = P0^* Box0(theta) P0.

Two spectral functionals are tested:

  1. raw   : W_raw(theta) = (1/2) log det' |Box_gauge(theta)|
  2. reduced: W_red(theta) = (1/2) log det |Box_gauge(theta)|_{mean-zero}

The reduced version removes the trivial base constant mode, so it remains
analytic if the raw prime-determinant develops a lifted near-zero mode.

As in exp606, we test the one-loop-style curvature matching

    Gamma(theta) = S_YM(theta)/(4*pi*alpha) + W(theta)

which implies the candidate coefficient

    alpha^{-1}_cand = -2*pi * W''(0) / K_YM

with K_YM = 6 from the exact simplicial Yang-Mills action.
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


def build_edge_lists(
    adj: np.ndarray, fibers: list[list[int]]
) -> tuple[list[tuple[int, int]], list[tuple[int, int]], set[tuple[int, int]]]:
    cycles = get_ordered_cycles(fibers, adj)
    fiber_forward: set[tuple[int, int]] = set()
    fiber_undirected: set[tuple[int, int]] = set()
    for cyc in cycles:
        for k in range(10):
            u = cyc[k]
            v = cyc[(k + 1) % 10]
            fiber_forward.add((u, v))
            fiber_undirected.add((u, v) if u < v else (v, u))

    fiber_edges: list[tuple[int, int]] = []
    cross_edges: list[tuple[int, int]] = []
    for i in range(N):
        for j in range(i + 1, N):
            if adj[i, j] <= 0.5:
                continue
            if (i, j) in fiber_undirected:
                fiber_edges.append((i, j))
            else:
                cross_edges.append((i, j))

    return fiber_edges, cross_edges, fiber_forward


def build_fiber_basis(fibers: list[list[int]]) -> np.ndarray:
    basis = np.zeros((N, len(fibers)))
    for fi, fiber in enumerate(fibers):
        basis[fiber, fi] = 1.0 / math.sqrt(len(fiber))
    return basis


def build_mean_zero_basis(nfib: int) -> np.ndarray:
    u = np.ones((nfib, 1)) / math.sqrt(nfib)
    proj = np.eye(nfib) - u @ u.T
    q, _ = np.linalg.qr(proj[:, :-1])
    return q[:, : nfib - 1]


def edge_phase(i: int, j: int, theta: float, fiber_forward: set[tuple[int, int]]) -> complex:
    if (i, j) in fiber_forward:
        return np.exp(1j * theta / 10.0)
    if (j, i) in fiber_forward:
        return np.exp(-1j * theta / 10.0)
    return 1.0 + 0.0j


def build_d0_fiber(
    theta: float, fiber_edges: list[tuple[int, int]], fiber_forward: set[tuple[int, int]]
) -> np.ndarray:
    d0 = np.zeros((len(fiber_edges), N), dtype=complex)
    for e_idx, (i, j) in enumerate(fiber_edges):
        d0[e_idx, i] = -edge_phase(i, j, theta, fiber_forward)
        d0[e_idx, j] = +1.0
    return d0


def build_d0_cross(cross_edges: list[tuple[int, int]]) -> np.ndarray:
    d0 = np.zeros((len(cross_edges), N), dtype=complex)
    for e_idx, (i, j) in enumerate(cross_edges):
        d0[e_idx, i] = -1.0
        d0[e_idx, j] = +1.0
    return d0


def build_d1(
    adj: np.ndarray, fiber_forward: set[tuple[int, int]]
) -> tuple[np.ndarray, list[tuple[int, int]], dict[tuple[int, int], int]]:
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

    d1 = np.zeros((len(triangles), len(edges)), dtype=complex)
    for f_idx, (i, j, k) in enumerate(triangles):
        d1[f_idx, edge_to_idx[(i, j)]] = +1.0
        d1[f_idx, edge_to_idx[(j, k)]] = +1.0
        d1[f_idx, edge_to_idx[(i, k)]] = -1.0
    return d1, edges, edge_to_idx


def build_full_d0(
    theta: float,
    all_edges: list[tuple[int, int]],
    fiber_forward: set[tuple[int, int]],
) -> np.ndarray:
    d0 = np.zeros((len(all_edges), N), dtype=complex)
    for e_idx, (i, j) in enumerate(all_edges):
        if (i, j) in fiber_forward or (j, i) in fiber_forward:
            d0[e_idx, i] = -edge_phase(i, j, theta, fiber_forward)
        else:
            d0[e_idx, i] = -1.0
        d0[e_idx, j] = +1.0
    return d0


def hermitian_eigs(mat: np.ndarray) -> np.ndarray:
    herm = 0.5 * (mat + mat.conj().T)
    return np.linalg.eigvalsh(herm)


def logdet_prime_abs(mat: np.ndarray) -> float:
    eigs = hermitian_eigs(mat)
    nz = eigs[np.abs(eigs) > TOL]
    return float(np.sum(np.log(np.abs(nz))))


def logdet_abs(mat: np.ndarray) -> float:
    eigs = hermitian_eigs(mat)
    return float(np.sum(np.log(np.abs(eigs))))


def smallest_nonzero_abs(mat: np.ndarray) -> float:
    eigs = hermitian_eigs(mat)
    nz = np.abs(eigs[np.abs(eigs) > TOL])
    return float(np.min(nz))


def kernel_dim(mat: np.ndarray) -> int:
    return int(np.sum(np.abs(hermitian_eigs(mat)) <= TOL))


def small_theta_quadratic(theta_vals: np.ndarray, values: np.ndarray) -> tuple[float, float, float]:
    x = theta_vals**2
    y = values - values[0]
    coeffs = np.polyfit(x, y, 2)
    fit_vals = coeffs[0] * x**2 + coeffs[1] * x + coeffs[2]
    second = 2.0 * coeffs[1]
    resid = float(np.max(np.abs(y - fit_vals)))
    return float(coeffs[1]), float(second), resid


print("=" * 78)
print("EXP607: TWISTED DETERMINANT ON BOX_GAUGE")
print("=" * 78)

verts, adj, _ = build_600cell()
fibers = find_fibration(verts)
fiber_edges, cross_edges, fiber_forward = build_edge_lists(adj, fibers)
basis = build_fiber_basis(fibers)
mean_zero = build_mean_zero_basis(len(fibers))
d0_cross = build_d0_cross(cross_edges)
l_cross = d0_cross.conj().T @ d0_cross

d1, all_edges, _ = build_d1(adj, fiber_forward)
theta_vals = np.array([0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.16])

raw_vals: list[float] = []
red_vals: list[float] = []
ym_vals: list[float] = []

print("\nGeometry")
print(f"  Fiber edges: {len(fiber_edges)}")
print(f"  Cross edges: {len(cross_edges)}")
print(f"  Gauge dimension: {basis.shape[1]}")
print(f"  Reduced gauge dimension: {mean_zero.shape[1]}")

for theta in theta_vals:
    d0_fiber = build_d0_fiber(theta, fiber_edges, fiber_forward)
    l_fiber = d0_fiber.conj().T @ d0_fiber
    box0 = l_cross - A1 * l_fiber
    box_gauge = basis.conj().T @ box0 @ basis
    box_red = mean_zero.conj().T @ box_gauge @ mean_zero

    raw_vals.append(0.5 * logdet_prime_abs(box_gauge))
    red_vals.append(0.5 * logdet_abs(box_red))

    d0_full = build_full_d0(theta, all_edges, fiber_forward)
    F = d1 @ d0_full
    ym_vals.append(float(np.real(np.trace(F.conj().T @ F))))

ym_vals_np = np.array(ym_vals)
ym_c2, ym_second, ym_resid = small_theta_quadratic(theta_vals, ym_vals_np)

print("\nYang-Mills control")
print(f"  K_YM from fit = {ym_c2:.12f}")
print(f"  S_YM''(0)     = {ym_second:.12f}")
print(f"  expected K_YM = {B1}")
print(f"  fit residual  = {ym_resid:.3e}")

theta_probe = 0.01
d0_fiber_0 = build_d0_fiber(0.0, fiber_edges, fiber_forward)
d0_fiber_p = build_d0_fiber(theta_probe, fiber_edges, fiber_forward)
box0_0 = l_cross - A1 * (d0_fiber_0.conj().T @ d0_fiber_0)
box0_p = l_cross - A1 * (d0_fiber_p.conj().T @ d0_fiber_p)
box_gauge_0 = basis.conj().T @ box0_0 @ basis
box_gauge_p = basis.conj().T @ box0_p @ basis
box_red_0 = mean_zero.conj().T @ box_gauge_0 @ mean_zero
box_red_p = mean_zero.conj().T @ box_gauge_p @ mean_zero

print("\nAnalyticity diagnostic")
print(
    f"  raw Box_gauge: ker(0)={kernel_dim(box_gauge_0):2d}, "
    f"ker({theta_probe:g})={kernel_dim(box_gauge_p):2d}, "
    f"smallest nonzero @ {theta_probe:g} = {smallest_nonzero_abs(box_gauge_p):.6e}"
)
print(
    f"  reduced Box_gauge: ker(0)={kernel_dim(box_red_0):2d}, "
    f"ker({theta_probe:g})={kernel_dim(box_red_p):2d}, "
    f"smallest nonzero @ {theta_probe:g} = {smallest_nonzero_abs(box_red_p):.6e}"
)

raw_vals_np = np.array(raw_vals)
red_vals_np = np.array(red_vals)
_, raw_second, raw_resid = small_theta_quadratic(theta_vals, raw_vals_np)
_, red_second, red_resid = small_theta_quadratic(theta_vals, red_vals_np)

alpha_inv_raw = -2.0 * math.pi * raw_second / B1
alpha_inv_red = -2.0 * math.pi * red_second / B1

print("\nCurvature matching for the twisted gauge determinant")
print(f"  raw W''(0)      = {raw_second:.12f}")
print(f"  raw fit residual= {raw_resid:.3e}")
print(f"  raw alpha^-1    = {alpha_inv_raw:.12f}")
print(f"  reduced W''(0)  = {red_second:.12f}")
print(f"  reduced residual= {red_resid:.3e}")
print(f"  reduced alpha^-1= {alpha_inv_red:.12f}")
print(f"  target tree     = {ALPHA_INV_TREE:.12f}")
print(f"  target full     = {ALPHA_INV_FULL:.12f}")

print("\nDistinguished angles")
for label, theta in [
    ("0", 0.0),
    ("pi/5", math.pi / 5.0),
    ("2pi", 2.0 * math.pi),
    ("20pi", 20.0 * math.pi),
]:
    d0_fiber = build_d0_fiber(theta, fiber_edges, fiber_forward)
    l_fiber = d0_fiber.conj().T @ d0_fiber
    box0 = l_cross - A1 * l_fiber
    box_gauge = basis.conj().T @ box0 @ basis
    box_red = mean_zero.conj().T @ box_gauge @ mean_zero
    eigs = hermitian_eigs(box_gauge)
    rounded = np.round(eigs, 10)
    print(
        f"  theta={label:>4s}: "
        f"logdet'|raw|={0.5 * logdet_prime_abs(box_gauge):12.6f}, "
        f"logdet|red|={0.5 * logdet_abs(box_red):12.6f}, "
        f"spec={rounded}"
    )

print("\nVerdict")
print("  This is the direct test on the intended sector: the fiber-constant")
print("  gauge projection of the twisted Box operator.")
print("  If the reduced determinant still fails to land near 1/alpha, then the")
print("  missing piece is probably not a simple spectral determinant of Box_gauge")
print("  alone, but a genuinely structural/topological ingredient beyond the")
print("  gauge-sector spectrum.")

r"""
exp569: Twisted edge operator on the 600-cell Hopf fibers.

Motivation
----------
The untwisted edge-space coexact operator C = d1^T d1 does not separate
cleanly into fiber/cross sectors. The exact lift

    Box1_lift = b1 d0 A_f d0^T - d0 A d0^T = d0 Box0 d0^T

contains no genuinely new spectral information: its nonzero spectrum equals
the spectrum of Box0 Delta0 on the vertex space.

To inject the Hopf holonomy directly into the 1-form sector, we instead build
a twisted incidence operator d_theta with total holonomy theta around each
decagonal Hopf fiber. Since each fiber has 10 edges, the edge holonomy is
theta/10, so the incidence row uses half-phases e^{\pm i theta/20}. This
gives two natural 1-form operators:

    B_theta   = d_theta d_theta^*              exact-sector energy
    Delta1    = B_theta + C                    positive twisted edge Laplacian
    Box1      = C - a1 B_theta                 signed exact/coexact contrast

The signed operator Box1(theta) is the main exploratory object: it allows
log det', eta, and spectral-flow diagnostics.
"""

from __future__ import annotations

from collections import Counter
import math
import sys

import numpy as np
from numpy.linalg import eigvalsh, norm

sys.path.insert(0, ".")
from commons import build_600cell


PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
TOL = 1e-8
N_SCAN = 81


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
    return idx if dots[idx] > 1 - tol else -1


def find_ordered_fibration(verts: np.ndarray, adj: np.ndarray) -> list[list[int]]:
    target_w = PHI / 2.0
    n = len(verts)
    for i in range(n):
        if abs(verts[i, 0] - target_w) >= 1e-6:
            continue

        g = verts[i]
        power = g.copy()
        ok = True
        for k in range(2, 11):
            power = qmul(power, g)
            if k == 5 and not np.allclose(power, [-1, 0, 0, 0], atol=1e-6):
                ok = False
                break
            if k == 10 and not np.allclose(power, [1, 0, 0, 0], atol=1e-6):
                ok = False
                break
        if not ok:
            continue

        subgroup = []
        p = np.array([1.0, 0, 0, 0])
        for _ in range(10):
            subgroup.append(find_idx(p, verts))
            p = qmul(p, g)

        used: set[int] = set()
        fibers: list[list[int]] = []
        for s in range(n):
            if s in used:
                continue
            fiber = []
            q = verts[s].copy()
            for _ in range(10):
                idx = find_idx(q, verts)
                if idx < 0 or idx in used:
                    break
                fiber.append(idx)
                used.add(idx)
                q = qmul(q, g)
            if len(fiber) == 10:
                fibers.append(fiber)

        if len(fibers) == 12 and len(used) == n:
            return fibers

    raise RuntimeError("Could not find ordered Hopf fibration")


def build_complex(verts: np.ndarray, adj: np.ndarray):
    n = len(verts)
    edges: list[tuple[int, int]] = []
    edge_to_idx: dict[tuple[int, int], int] = {}
    adj_list = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j] > 0.5:
                edge_to_idx[(i, j)] = len(edges)
                edges.append((i, j))
                adj_list[i].add(j)
                adj_list[j].add(i)

    triangles: list[tuple[int, int, int]] = []
    for i in range(n):
        for j in adj_list[i]:
            if j <= i:
                continue
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    triangles.append((i, j, k))

    return edges, edge_to_idx, triangles


def build_face_boundary(
    n_edges: int, triangles: list[tuple[int, int, int]], edge_to_idx: dict[tuple[int, int], int]
) -> np.ndarray:
    d1 = np.zeros((len(triangles), n_edges), dtype=float)
    for f_idx, (i, j, k) in enumerate(triangles):
        d1[f_idx, edge_to_idx[(i, j)]] = +1.0
        d1[f_idx, edge_to_idx[(j, k)]] = +1.0
        d1[f_idx, edge_to_idx[(i, k)]] = -1.0
    return d1


def build_fiber_orientation(
    fibers: list[list[int]], edge_to_idx: dict[tuple[int, int], int]
) -> tuple[np.ndarray, np.ndarray]:
    n_edges = len(edge_to_idx)
    edge_sign = np.zeros(n_edges, dtype=int)
    is_fiber = np.zeros(n_edges, dtype=bool)

    for fiber in fibers:
        for k in range(10):
            u = fiber[k]
            v = fiber[(k + 1) % 10]
            i, j = (u, v) if u < v else (v, u)
            e_idx = edge_to_idx[(i, j)]
            is_fiber[e_idx] = True
            edge_sign[e_idx] = +1 if (u, v) == (i, j) else -1

    return is_fiber, edge_sign


def build_twisted_d0(
    n_vertices: int,
    edges: list[tuple[int, int]],
    is_fiber: np.ndarray,
    edge_sign: np.ndarray,
    theta: float,
) -> np.ndarray:
    d0 = np.zeros((len(edges), n_vertices), dtype=complex)
    for e_idx, (i, j) in enumerate(edges):
        if is_fiber[e_idx]:
            sgn = edge_sign[e_idx]
            phase_minus = np.exp(-0.05j * sgn * theta)
            phase_plus = np.exp(+0.05j * sgn * theta)
            d0[e_idx, i] = -phase_minus
            d0[e_idx, j] = +phase_plus
        else:
            d0[e_idx, i] = -1.0
            d0[e_idx, j] = +1.0
    return d0


def spectral_summary(evals: np.ndarray) -> dict[str, float]:
    nonzero = evals[np.abs(evals) > 1e-7]
    pos = nonzero[nonzero > 0]
    neg = nonzero[nonzero < 0]
    out = {
        "kernel": float(len(evals) - len(nonzero)),
        "eta0": float(len(pos) - len(neg)),
        "logdet": float(np.sum(np.log(np.abs(nonzero)))) if len(nonzero) else float("-inf"),
        "gap_abs": float(np.min(np.abs(nonzero))) if len(nonzero) else 0.0,
        "tr2": float(np.sum(evals**2)),
    }
    return out


def second_derivative_at_zero(xs: np.ndarray, ys: np.ndarray) -> float:
    # Symmetric 3-point estimate using x = [0, h, 2h, ...]
    h = xs[1] - xs[0]
    return float((ys[2] - 2 * ys[1] + ys[0]) / (h * h))


print("=" * 72)
print("EXP569: TWISTED EDGE OPERATOR")
print("=" * 72)

verts, adj, lap = build_600cell()
edges, edge_to_idx, triangles = build_complex(verts, adj)
d1 = build_face_boundary(len(edges), triangles, edge_to_idx)
C = d1.T @ d1
fibers = find_ordered_fibration(verts, adj)
is_fiber, edge_sign = build_fiber_orientation(fibers, edge_to_idx)

print("\nGeometry")
print(f"  Vertices:  {len(verts)}")
print(f"  Edges:     {len(edges)}")
print(f"  Triangles: {len(triangles)}")
print(f"  Fiber edges: {np.sum(is_fiber)} (expect 120)")
print(f"  Cross edges: {len(edges) - np.sum(is_fiber)} (expect 600)")
print(f"  ||d1 d0(0)|| check will be nonzero once twisted away from 0; at theta=0 we recover the ordinary complex.")

d0_0 = build_twisted_d0(len(verts), edges, is_fiber, edge_sign, 0.0)
err_complex = norm(d1 @ d0_0)
B0 = d0_0 @ d0_0.conj().T
Box1_0 = C - a1 * B0
evals0 = np.sort(eigvalsh(Box1_0))

print("\nUntwisted Checks")
print(f"  ||d1 d0(0)|| = {err_complex:.2e}")
print(f"  ker(B0)      = {np.sum(np.abs(eigvalsh(B0)) < 1e-7)} (expect 601)")
print(f"  ker(C)       = {np.sum(np.abs(eigvalsh(C)) < 1e-7)} (expect 119)")
print(f"  ker(Box1(0)) = {np.sum(np.abs(evals0) < 1e-7)}")

abs_counter = Counter(round(float(x), 6) for x in np.sort(np.abs(evals0[np.abs(evals0) > 1e-7])))
print("  Distinct |eigenvalues| of Box1(0):")
for val, mult in sorted(abs_counter.items()):
    print(f"    {val:10.6f} : {mult:3d}")

thetas = np.linspace(0.0, 2 * np.pi, N_SCAN)
box_logdet = np.zeros(N_SCAN)
box_eta = np.zeros(N_SCAN)
box_gap = np.zeros(N_SCAN)
box_tr2 = np.zeros(N_SCAN)
lap_logdet = np.zeros(N_SCAN)
lap_gap = np.zeros(N_SCAN)
lap_tr2 = np.zeros(N_SCAN)
box_kernel = np.zeros(N_SCAN)

for idx, theta in enumerate(thetas):
    d0_theta = build_twisted_d0(len(verts), edges, is_fiber, edge_sign, theta)
    B_theta = d0_theta @ d0_theta.conj().T
    Delta1_theta = B_theta + C
    Box1_theta = C - a1 * B_theta

    evals_box = np.sort(eigvalsh(Box1_theta))
    evals_lap = np.sort(eigvalsh(Delta1_theta))
    s_box = spectral_summary(evals_box)
    s_lap = spectral_summary(evals_lap)

    box_logdet[idx] = s_box["logdet"]
    box_eta[idx] = s_box["eta0"]
    box_gap[idx] = s_box["gap_abs"]
    box_tr2[idx] = s_box["tr2"]
    box_kernel[idx] = s_box["kernel"]

    lap_logdet[idx] = s_lap["logdet"]
    lap_gap[idx] = s_lap["gap_abs"]
    lap_tr2[idx] = s_lap["tr2"]

print("\nSpecial Angles")
special_angles = [
    ("0", 0.0),
    ("pi/10", np.pi / 10),
    ("pi/5", np.pi / 5),
    ("2pi/5", 2 * np.pi / 5),
    ("pi", np.pi),
    ("2pi", 2 * np.pi),
]
for label, theta in special_angles:
    j = int(np.argmin(np.abs(thetas - theta)))
    print(
        f"  theta={label:>5s} ({thetas[j]:8.5f}) : "
        f"logdet Box1={box_logdet[j]:10.4f}, eta={box_eta[j]:6.0f}, "
        f"gap={box_gap[j]:9.6f}, ker={box_kernel[j]:5.0f}, "
        f"logdet Delta1={lap_logdet[j]:10.4f}"
    )

print("\nGlobal Scan")
imin_box = int(np.argmin(box_logdet))
imax_box = int(np.argmax(box_logdet))
imin_gap = int(np.argmin(box_gap))
imax_gap = int(np.argmax(box_gap))
imin_lap = int(np.argmin(lap_logdet))
imax_lap = int(np.argmax(lap_logdet))
print(f"  Box1 logdet min at theta={thetas[imin_box]:.6f}, value={box_logdet[imin_box]:.6f}")
print(f"  Box1 logdet max at theta={thetas[imax_box]:.6f}, value={box_logdet[imax_box]:.6f}")
print(f"  Box1 |gap| min   at theta={thetas[imin_gap]:.6f}, value={box_gap[imin_gap]:.6f}")
print(f"  Box1 |gap| max   at theta={thetas[imax_gap]:.6f}, value={box_gap[imax_gap]:.6f}")
print(f"  Delta1 logdet min at theta={thetas[imin_lap]:.6f}, value={lap_logdet[imin_lap]:.6f}")
print(f"  Delta1 logdet max at theta={thetas[imax_lap]:.6f}, value={lap_logdet[imax_lap]:.6f}")

eta_jumps = []
for i in range(N_SCAN - 1):
    if abs(box_eta[i + 1] - box_eta[i]) > 0.5:
        eta_jumps.append((thetas[i], thetas[i + 1], box_eta[i], box_eta[i + 1]))

print("\nEta / Spectral Flow")
if eta_jumps:
    for lo, hi, e0, e1 in eta_jumps:
        print(
            f"  jump near theta in [{lo:.5f}, {hi:.5f}] : "
            f"eta {e0:.0f} -> {e1:.0f}"
        )
else:
    print("  No eta jump on the scan grid.")

theta_small = np.linspace(0.0, 0.20, 21)
box_logdet_small = np.zeros_like(theta_small)
box_tr2_small = np.zeros_like(theta_small)
lap_logdet_small = np.zeros_like(theta_small)
for i, theta in enumerate(theta_small):
    d0_theta = build_twisted_d0(len(verts), edges, is_fiber, edge_sign, theta)
    B_theta = d0_theta @ d0_theta.conj().T
    Delta1_theta = B_theta + C
    Box1_theta = C - a1 * B_theta
    evals_box = np.sort(eigvalsh(Box1_theta))
    evals_lap = np.sort(eigvalsh(Delta1_theta))
    s_box = spectral_summary(evals_box)
    s_lap = spectral_summary(evals_lap)
    box_logdet_small[i] = s_box["logdet"]
    box_tr2_small[i] = s_box["tr2"]
    lap_logdet_small[i] = s_lap["logdet"]

curv_box_logdet = second_derivative_at_zero(theta_small, box_logdet_small)
curv_box_tr2 = second_derivative_at_zero(theta_small, box_tr2_small)
curv_lap_logdet = second_derivative_at_zero(theta_small, lap_logdet_small)

print("\nSmall-theta Curvature")
print(f"  d^2/dtheta^2 logdet Box1 |0  ~= {curv_box_logdet:.6f}")
print(f"  d^2/dtheta^2 Tr(Box1^2) |0   ~= {curv_box_tr2:.6f}")
print(f"  d^2/dtheta^2 logdet Delta1|0 ~= {curv_lap_logdet:.6f}")
print(f"  Ratios to 2pi: {curv_box_logdet / (2*np.pi):.6f}, {curv_box_tr2 / (2*np.pi):.6f}, {curv_lap_logdet / (2*np.pi):.6f}")

print("\nPeriodicity")
print(f"  Box1 logdet(2pi)-logdet(0)   = {box_logdet[-1] - box_logdet[0]:.6e}")
print(f"  Box1 gap(2pi)-gap(0)         = {box_gap[-1] - box_gap[0]:.6e}")
print(f"  Delta1 logdet(2pi)-logdet(0) = {lap_logdet[-1] - lap_logdet[0]:.6e}")

print("\nVerdict")
print("  The holonomy now enters directly through theta in the 1-form sector.")
print("  If a clean alpha-like coefficient exists here, it should show up through")
print("  periodicity, a special angle such as pi/a1, or a distinguished curvature")
print("  of log det' / eta / spectral gap. This script is a diagnostic for that.")
print("  The missing exact 2pi periodicity is itself informative: twisting only")
print("  d_theta while keeping the face operator untwisted is not a fully gauge-")
print("  consistent magnetic complex. The next refinement would be a compatible")
print("  twist of the 2-simplex sector so that the discrete connection closes on")
print("  faces, not only on fiber edges.")

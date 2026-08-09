"""
exp605: Variational functional for alpha from Box_gauge + Hopf holonomy.

Goal
----
exp604 showed that the full quadratic alpha equation closes naturally from:

    spectral coefficient A := 1/alpha_0  from Box_gauge
    topological coefficient H := 2*pi    from the closed Hopf fiber

This script asks a sharper question:

    Does there exist a simple mixed functional J(alpha) whose stationarity
    condition is exactly the alpha equation, and is the physical root a stable
    minimum on the weak-coupling domain 0 < alpha < 1?

We use the minimal positive-coupling functional

    J_H(alpha) = A*alpha - (H/2)*alpha^2 - ln(alpha),   alpha > 0

Then

    dJ/dalpha = A - H*alpha - 1/alpha
              = 0

is equivalent to

    H*alpha^2 - A*alpha + 1 = 0.

Interpretation
--------------
  + A*alpha         : spectral drive from Box_gauge
  - (H/2)*alpha^2   : topological holonomy cost
  - ln(alpha)       : positivity / normalization barrier

This does NOT prove that J is the unique physical action. It only shows that
the full alpha equation can be realized as the Euler-Lagrange condition of a
natural mixed functional with no fitted coefficients.
"""

from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, ".")
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0
A1 = 5
B1 = 6
CODATA_ALPHA_INV = 137.035999084
TOL = 1e-10


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


def find_hopf_fibration(verts: np.ndarray) -> list[list[int]]:
    target_w = PHI / 2.0
    nv = len(verts)
    for i in range(nv):
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

        used = set()
        fibers: list[list[int]] = []
        subgroup = []
        pp = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_idx(pp, verts))
            pp = qmul(pp, g)

        for s in range(nv):
            if s in used:
                continue
            fib = []
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


def build_box_gauge() -> np.ndarray:
    verts, adj, _ = build_600cell()
    fibers = find_hopf_fibration(verts)
    nv = len(verts)

    vertex_to_fiber: dict[int, int] = {}
    for fi, fiber in enumerate(fibers):
        for v in fiber:
            vertex_to_fiber[v] = fi

    a_fiber = adj * np.array(
        [
            [float(vertex_to_fiber[i] == vertex_to_fiber[j]) for j in range(nv)]
            for i in range(nv)
        ]
    )
    a_cross = adj - a_fiber
    l_fiber = np.diag(np.sum(a_fiber, axis=1)) - a_fiber
    l_cross = np.diag(np.sum(a_cross, axis=1)) - a_cross
    box0 = l_cross - A1 * l_fiber

    basis = np.zeros((nv, len(fibers)))
    for fi, fiber in enumerate(fibers):
        basis[fiber, fi] = 1.0 / math.sqrt(len(fiber))
    return basis.T @ box0 @ basis


def physical_root(A: float, H: float) -> float:
    disc = A * A - 4.0 * H
    return (A - math.sqrt(disc)) / (2.0 * H)


def large_root(A: float, H: float) -> float:
    disc = A * A - 4.0 * H
    return (A + math.sqrt(disc)) / (2.0 * H)


def J(alpha: float, A: float, H: float) -> float:
    return A * alpha - 0.5 * H * alpha * alpha - math.log(alpha)


def dJ(alpha: float, A: float, H: float) -> float:
    return A - H * alpha - 1.0 / alpha


def d2J(alpha: float, H: float) -> float:
    return 1.0 / (alpha * alpha) - H


print("=" * 78)
print("VARIATIONAL FUNCTIONAL FOR ALPHA")
print("=" * 78)

print("\nBuilding Box_gauge...")
box_gauge = build_box_gauge()
tr_box_g2 = float(np.sum(np.linalg.eigvalsh(box_gauge) ** 2))
lambda_fiber = 2.0 - PHI
A = tr_box_g2 / (12.0 * B1 * lambda_fiber**2)
theta_edge = math.acos(1.0 - lambda_fiber / 2.0)
H = (2 * A1) * theta_edge

alpha_phys = physical_root(A, H)
alpha_large = large_root(A, H)
alpha_inv = 1.0 / alpha_phys

print(f"  Spectral coefficient A = 1/alpha_0 = {A:.12f}")
print(f"  Topological coefficient H = 2*pi   = {H:.12f}")
print(f"  Physical root alpha                = {alpha_phys:.15f}")
print(f"  Physical root alpha^-1             = {alpha_inv:.12f}")
print(f"  CODATA alpha^-1                    = {CODATA_ALPHA_INV:.12f}")
print(f"  Large root                         = {alpha_large:.12f}")

print("\n" + "=" * 78)
print("1. STATIONARITY = ALPHA EQUATION")
print("=" * 78)

print(f"  dJ(alpha_phys)  = {dJ(alpha_phys, A, H):+.3e}")
print(f"  dJ(alpha_large) = {dJ(alpha_large, A, H):+.3e}")
print(f"  d2J(alpha_phys) = {d2J(alpha_phys, H):.12f}")
print(f"  d2J(alpha_large)= {d2J(alpha_large, H):.12f}")

print("\nInterpretation:")
print("  The weak-coupling root is a local minimum because d2J > 0 there.")
print("  The large root is a local maximum because d2J < 0 there.")

print("\n" + "=" * 78)
print("2. GLOBAL SHAPE ON THE PHYSICAL DOMAIN 0 < alpha < 1")
print("=" * 78)

grid = np.geomspace(1e-5, 1.0, 4000)
vals = np.array([J(x, A, H) for x in grid])
idx_min = int(np.argmin(vals))
alpha_grid_min = float(grid[idx_min])

print(f"  Grid minimum on (0,1] at alpha = {alpha_grid_min:.15f}")
print(f"  J(alpha_phys)              = {J(alpha_phys, A, H):.12f}")
print(f"  J(alpha_grid_min)          = {vals[idx_min]:.12f}")
print(f"  |alpha_grid_min-alpha_phys| = {abs(alpha_grid_min - alpha_phys):.3e}")

print("\nBoundary behavior:")
print(f"  J(1e-5) = {J(1e-5, A, H):.12f}")
print(f"  J(1.0)  = {J(1.0, A, H):.12f}")
print("  So on the physical weak-coupling interval, the smaller root is the unique")
print("  global minimum of this functional.")

print("\n" + "=" * 78)
print("3. SENSITIVITY TO THE HOLOMONY MULTIPLE n * theta_edge")
print("=" * 78)

print(f"{'n':>3s} {'H_n':>16s} {'alpha_inv(n)':>16s} {'J_min alpha':>16s} {'status':>16s}")
best_n = None
best_err = None
for n in range(1, 21):
    h_n = n * theta_edge
    disc = A * A - 4.0 * h_n
    if disc <= 0:
        continue
    alpha_n = physical_root(A, h_n)
    err = abs(1.0 / alpha_n - CODATA_ALPHA_INV) / CODATA_ALPHA_INV
    status = ""
    if n == 2 * A1:
        status = "<-- closed fiber"
    if best_err is None or err < best_err:
        best_err = err
        best_n = n
    print(f"{n:3d} {h_n:16.12f} {1.0/alpha_n:16.9f} {alpha_n:16.12f} {status:>16s}")

print(f"\n  Best n in [1,20] by alpha agreement: {best_n}")
print(f"  Closed-fiber n: {2*A1}")

print("\n" + "=" * 78)
print("4. WHAT THE FUNCTIONAL ADDS")
print("=" * 78)

print("  exp604 showed that the alpha equation closes from spectral + topological")
print("  coefficients. This script shows more:")
print("    - the same equation is the Euler-Lagrange condition of a simple mixed")
print("      functional J(alpha)")
print("    - the physical root is a stable minimum on the weak-coupling domain")
print("    - the larger algebraic root is automatically unstable")

print("\nCaveat:")
print("  J(alpha) is presently a mathematically natural candidate functional, not a")
print("  derived action from first principles. What is gained is variational form,")
print("  stability selection, and a cleaner bridge from 'equation' to 'principle'.")

"""
verify_gauge_continuum_map.py -- Low-harmonic S^2 map on the Hopf base.

Verifies that the 12 Hopf fibers define an icosahedron graph whose scalar
low-harmonic content matches the round S^2 exactly up to l = 2:

  (a) Hopf-base graph is the icosahedron: 12 vertices, degree 5, 30 edges
  (b) Base Laplacian spectrum is 0 + (5-sqrt(5))^3 + 6^5 + (5+sqrt(5))^3
  (c) Spectral embedding from the first 3-dimensional eigenspace gives a
      regular icosahedron on S^2
  (d) Sampled scalar harmonics have rank 9 = 1+3+5
  (e) l=1 and l=2 spaces are exact invariant eigenspaces
  (f) The sampled low-harmonic projector matches the spectral projector
  (g) A quadratic polynomial in the base Laplacian reproduces the continuum
      Casimir values l(l+1)=0,2,6 on the low-harmonic subspace

This is a controlled scalar continuum map for the gauge skeleton. It does not
claim a full vector-harmonic or Lie-bracket completion.
"""

import os
import sys
import math

import numpy as np
from numpy.linalg import eigh, matrix_rank, norm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0
TOL = 1e-8
PASS = 0
tests_run = 0
tests_pass = 0


def check(name, condition, detail=""):
    global PASS, tests_run, tests_pass
    tests_run += 1
    if condition:
        tests_pass += 1
        print(f"  [PASS] {name}")
    else:
        PASS = 1
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def qmul(p, q):
    return np.array(
        [
            p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
            p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
            p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
            p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0],
        ]
    )


def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = int(np.argmax(dots))
    return idx if dots[idx] > 1.0 - tol else -1


def find_hopf_fibration(verts):
    nv = len(verts)
    for i in range(nv):
        if abs(verts[i, 0] - PHI / 2.0) >= 1e-6:
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
        if not ok:
            continue

        subgroup = []
        power = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_idx(power, verts))
            power = qmul(power, g)

        used = set()
        fibers = []
        for s in range(nv):
            if s in used:
                continue
            fiber = []
            for si in subgroup:
                idx = find_idx(qmul(verts[s], verts[si]), verts)
                if idx >= 0 and idx not in used:
                    fiber.append(idx)
                    used.add(idx)
            if len(fiber) == 10:
                fibers.append(fiber)

        if len(fibers) == 12:
            return fibers

    return None


def build_base_graph(adj, fibers):
    vertex_to_fiber = {}
    for fi, fiber in enumerate(fibers):
        for v in fiber:
            vertex_to_fiber[v] = fi

    weighted = np.zeros((12, 12), dtype=float)
    nv = adj.shape[0]
    for i in range(nv):
        for j in range(i + 1, nv):
            if adj[i, j] > 0.5:
                fi = vertex_to_fiber[i]
                fj = vertex_to_fiber[j]
                if fi != fj:
                    weighted[fi, fj] += 1.0
                    weighted[fj, fi] += 1.0

    adjacency = (weighted > 0).astype(float)
    laplacian = np.diag(np.sum(adjacency, axis=1)) - adjacency
    return weighted, adjacency, laplacian


def harmonic_basis_from_embedding(xyz):
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    return np.column_stack(
        [
            np.ones(12),
            x,
            y,
            z,
            x * y,
            y * z,
            z * x,
            x * x - y * y,
            2.0 * z * z - x * x - y * y,
        ]
    )


def projector_from_columns(mat):
    gram = mat.T @ mat
    return mat @ np.linalg.inv(gram) @ mat.T


def subspace_action(op, basis):
    gram = basis.T @ basis
    action = np.linalg.solve(gram, basis.T @ op @ basis)
    residual = norm(op @ basis - basis @ action)
    return action, residual


print("Building Hopf-base gauge continuum map...")
verts, adj, _ = build_600cell()
fibers = find_hopf_fibration(verts)
check("Found 12 Hopf fibers", fibers is not None and len(fibers) == 12)

weighted, adjacency, laplacian = build_base_graph(adj, fibers)

check("Base graph has 12 vertices", adjacency.shape == (12, 12))
check(
    "Base graph is 5-regular",
    np.allclose(np.sum(adjacency, axis=1), 5.0, atol=1e-12),
    f"degrees = {np.sum(adjacency, axis=1).astype(int).tolist()}",
)
check("Base graph has 30 edges", int(np.sum(adjacency) // 2) == 30)
check(
    "Every adjacent fiber pair has weight 20",
    np.array_equal(np.unique(weighted[weighted > 0]), np.array([20.0])),
    f"nonzero weights = {np.unique(weighted[weighted > 0])}",
)

evals, evecs = eigh(laplacian)
expected = np.array(
    [0.0]
    + [5.0 - math.sqrt(5.0)] * 3
    + [6.0] * 5
    + [5.0 + math.sqrt(5.0)] * 3,
    dtype=float,
)
check(
    "Base Laplacian matches icosahedron spectrum",
    np.allclose(evals, expected, atol=1e-8),
    f"evals = {np.round(evals, 8)}",
)

xyz = evecs[:, 1:4]
row_norms = np.linalg.norm(xyz, axis=1)
check(
    "Spectral embedding row norms are constant",
    np.allclose(row_norms, 0.5, atol=1e-12),
    f"spread = [{row_norms.min():.12f}, {row_norms.max():.12f}]",
)
xyz = xyz / row_norms[:, None]

dots = []
for i in range(12):
    for j in range(i + 1, 12):
        dots.append(np.dot(xyz[i], xyz[j]))
unique_dots = np.unique(np.round(np.sort(dots), 8))
expected_dots = np.array([-1.0, -1.0 / math.sqrt(5.0), 1.0 / math.sqrt(5.0)])
check(
    "Embedding gives regular icosahedron dot products",
    np.allclose(unique_dots, np.round(expected_dots, 8), atol=1e-8),
    f"dots = {unique_dots}",
)

harmonics = harmonic_basis_from_embedding(xyz)
check("Sampled l<=2 harmonic space has rank 9", matrix_rank(harmonics, tol=1e-10) == 9)

l0 = harmonics[:, :1]
l1 = harmonics[:, 1:4]
l2 = harmonics[:, 4:]

action_1, residual_1 = subspace_action(laplacian, l1)
action_2, residual_2 = subspace_action(laplacian, l2)
check(
    "l=1 subspace is exact eigenspace",
    np.allclose(np.linalg.eigvalsh((action_1 + action_1.T) / 2.0), 5.0 - math.sqrt(5.0), atol=1e-8)
    and residual_1 < 1e-10,
    f"residual = {residual_1:.3e}",
)
check(
    "l=2 subspace is exact eigenspace",
    np.allclose(np.linalg.eigvalsh((action_2 + action_2.T) / 2.0), 6.0, atol=1e-8)
    and residual_2 < 1e-10,
    f"residual = {residual_2:.3e}",
)

low_mask = (
    (np.abs(evals) < TOL)
    | (np.abs(evals - (5.0 - math.sqrt(5.0))) < TOL)
    | (np.abs(evals - 6.0) < TOL)
)
proj_spec = evecs[:, low_mask] @ evecs[:, low_mask].T
proj_harm = projector_from_columns(harmonics)
check(
    "Sampled low harmonics equal spectral projector 1+3+5",
    norm(proj_spec - proj_harm) < 1e-10,
    f"||P_spec-P_harm|| = {norm(proj_spec - proj_harm):.3e}",
)
check(
    "Alias sector has dimension 3",
    int(np.sum(np.abs(evals - (5.0 + math.sqrt(5.0))) < TOL)) == 3,
)

alpha, beta = np.linalg.solve(
    np.array(
        [
            [(5.0 - math.sqrt(5.0)) ** 2, 5.0 - math.sqrt(5.0)],
            [36.0, 6.0],
        ],
        dtype=float,
    ),
    np.array([2.0, 6.0], dtype=float),
)
casimir = alpha * (laplacian @ laplacian) + beta * laplacian
res0 = norm(casimir @ l0 - 0.0 * l0)
res1 = norm(casimir @ l1 - 2.0 * l1)
res2 = norm(casimir @ l2 - 6.0 * l2)
check(
    "Quadratic polynomial reproduces continuum Casimir on l<=2",
    max(res0, res1, res2) < 1e-10,
    f"residuals = ({res0:.3e}, {res1:.3e}, {res2:.3e})",
)

print("\n" + "=" * 60)
print(f"TOTAL: {tests_pass}/{tests_run} tests PASSED")
print("=" * 60)
sys.exit(PASS)

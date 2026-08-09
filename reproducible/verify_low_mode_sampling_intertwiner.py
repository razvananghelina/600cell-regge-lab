"""Verify the exact l<=2 sampling/intertwining theorem on the Hopf base.

The theorem is finite and exact up to floating-point certification:

* the 12 Hopf fibers carry the icosahedron graph;
* evaluation E of real spherical polynomials of degrees 0, 1, and 2 is
  injective (rank 9);
* C_disc = p(L_base), with p(0)=0, p(5-sqrt(5))=2, p(6)=6,
  satisfies C_disc E = E C_cont on H_0 + H_1 + H_2;
* the Moore-Penrose reconstruction R gives R E = I and E R equal to the
  spectral projector onto 1+3+5.

This does not assert convergence, a vector-field map, or a gauge connection.
"""

import math
import os
import sys

import numpy as np
from numpy.linalg import eigh, matrix_rank, norm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0
TOL = 1.0e-9
failed = False
tests_run = 0
tests_passed = 0


def check(name, condition, detail=""):
    global failed, tests_run, tests_passed
    tests_run += 1
    if condition:
        tests_passed += 1
        print(f"  [PASS] {name}")
    else:
        failed = True
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0],
    ])


def find_idx(v, verts, tol=1.0e-6):
    dots = verts @ v
    idx = int(np.argmax(dots))
    return idx if dots[idx] > 1.0 - tol else -1


def find_fibers(verts):
    for candidate in verts:
        if abs(candidate[0] - PHI / 2.0) >= 1.0e-6:
            continue
        power = candidate.copy()
        valid = True
        for k in range(2, 11):
            power = qmul(power, candidate)
            if k == 5 and not np.allclose(power, [-1, 0, 0, 0], atol=1.0e-6):
                valid = False
                break
            if k == 10 and not np.allclose(power, [1, 0, 0, 0], atol=1.0e-6):
                valid = False
        if not valid:
            continue

        subgroup = []
        power = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_idx(power, verts))
            power = qmul(power, candidate)

        used = set()
        fibers = []
        for seed in range(len(verts)):
            if seed in used:
                continue
            fiber = []
            for element in subgroup:
                idx = find_idx(qmul(verts[seed], verts[element]), verts)
                if idx >= 0 and idx not in used:
                    fiber.append(idx)
                    used.add(idx)
            if len(fiber) == 10:
                fibers.append(fiber)
        if len(fibers) == 12:
            return fibers
    raise RuntimeError("No 12-fiber Hopf decomposition found")


def base_laplacian(adj, fibers):
    vertex_to_fiber = {}
    for fiber_index, fiber in enumerate(fibers):
        for vertex in fiber:
            vertex_to_fiber[vertex] = fiber_index
    base_adj = np.zeros((12, 12))
    for i in range(len(adj)):
        for j in range(i + 1, len(adj)):
            if adj[i, j] > 0.5:
                fi, fj = vertex_to_fiber[i], vertex_to_fiber[j]
                if fi != fj:
                    base_adj[fi, fj] = base_adj[fj, fi] = 1.0
    return np.diag(base_adj.sum(axis=1)) - base_adj


print("Verifying exact low-mode sampling/intertwining theorem...")
verts, adj, _ = build_600cell()
fibers = find_fibers(verts)
lap = base_laplacian(adj, fibers)
evals, evecs = eigh(lap)

# The first nonconstant eigenspace is an equivariant spectral embedding.
xyz = evecs[:, 1:4].copy()
xyz /= np.linalg.norm(xyz, axis=1)[:, None]
x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]

# Evaluation of a real polynomial basis for H_0 + H_1 + H_2.
evaluation = np.column_stack([
    np.ones(12), x, y, z, x*y, y*z, z*x,
    x*x-y*y, 2.0*z*z-x*x-y*y,
])
continuum_casimir = np.diag([0.0] + [2.0]*3 + [6.0]*5)

sqrt5 = math.sqrt(5.0)
quad = (3.0 - sqrt5) / (4.0 * sqrt5)
linear = (5.0 * sqrt5 - 9.0) / (2.0 * sqrt5)
discrete_casimir = quad * (lap @ lap) + linear * lap

check("evaluation has rank 9", matrix_rank(evaluation, tol=1.0e-10) == 9)
intertwining_error = norm(discrete_casimir @ evaluation - evaluation @ continuum_casimir)
check("C_disc E = E C_cont on l<=2", intertwining_error < TOL,
      f"residual = {intertwining_error:.3e}")

reconstruction = np.linalg.inv(evaluation.T @ evaluation) @ evaluation.T
left_error = norm(reconstruction @ evaluation - np.eye(9))
check("sampling has exact left reconstruction", left_error < TOL,
      f"||R E-I|| = {left_error:.3e}")

low_basis = evecs[:, :9]
spectral_projector = low_basis @ low_basis.T
projector_error = norm(evaluation @ reconstruction - spectral_projector)
check("E R is the 1+3+5 spectral projector", projector_error < TOL,
      f"||E R-P_low|| = {projector_error:.3e}")

alias = evecs[:, 9:]
alias_error = norm(reconstruction @ alias)
check("reconstruction annihilates the 3-dimensional alias sector", alias_error < TOL,
      f"||R P_alias|| = {alias_error:.3e}")

print("\nClassification:")
print("  DERIVED: exact finite sampling, reconstruction, and Casimir intertwining on l<=2")
print("  OPEN: vector modes, local gauge transport, Lie bracket, and refinement convergence")
print(f"\nTOTAL: {tests_passed}/{tests_run} tests PASSED")
sys.exit(1 if failed else 0)

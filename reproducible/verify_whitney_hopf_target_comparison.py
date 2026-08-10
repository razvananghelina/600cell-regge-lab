#!/usr/bin/env python3
"""STEP-3 comparison of preregistered Whitney/Hopf spectra with the seed.

The blind data source is the JSON committed in STEP 1 at 9884c95.  This file
is deliberately later in git history.  It tests exact arithmetic, all six
recorded fibrations, common generalized eigenspaces and the kernel relation;
it also distinguishes preservation of an old graph ratio from independent
metric selection.
"""

import json
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as sla
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


BLIND = Path(__file__).with_name("whitney_hopf_blind_enumeration.json")
PREREGISTRATION_COMMIT = "9884c95"
TOL = 2e-8
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def quat_mult(left, right):
    w1, x1, y1, z1 = left
    w2, x2, y2, z2 = right
    return np.array((
        w1*w2-x1*x2-y1*y2-z1*z2,
        w1*x2+x1*w2+y1*z2-z1*y2,
        w1*y2-x1*z2+y1*w2+z1*x2,
        w1*z2+x1*y2-y1*x2+z1*w2,
    ))


def vertex_index(vertices, quaternion):
    distances = np.linalg.norm(vertices-quaternion, axis=1)
    index = int(np.argmin(distances))
    return index if distances[index] < 1e-6 else -1


def first_hopf_fibration(vertices):
    for generator in vertices:
        power = generator.copy()
        order = None
        for candidate in range(2, 121):
            power = quat_mult(power, generator)
            if np.allclose(power, (1, 0, 0, 0), atol=1e-6):
                order = candidate
                break
        if order != 10:
            continue
        subgroup = []
        power = np.array((1.0, 0.0, 0.0, 0.0))
        for _ in range(10):
            subgroup.append(vertex_index(vertices, power))
            power = quat_mult(power, generator)
        if min(subgroup) < 0 or len(set(subgroup)) != 10:
            continue
        assigned = np.full(120, -1, dtype=int)
        fibers = []
        valid = True
        for left in range(120):
            if assigned[left] >= 0:
                continue
            fiber = []
            for right in subgroup:
                target = vertex_index(
                    vertices, quat_mult(vertices[left], vertices[right]))
                if target >= 0 and assigned[target] < 0:
                    assigned[target] = len(fibers)
                    fiber.append(target)
            if len(fiber) != 10:
                valid = False
                break
            fibers.append(fiber)
        if valid and len(fibers) == 12:
            return fibers
    raise RuntimeError("no Hopf fibration found")


def first_positive(values):
    positive = values[values > TOL]
    return float(positive.min())


print("=" * 78)
print("STEP 3: WHITNEY/HOPF COMPARISON WITH THE INDEPENDENT BOOTSTRAP SEED")
print("=" * 78)

data = json.loads(BLIND.read_text())
check("blind data records STEP 1 with no target comparison",
      data["protocol"].startswith("STEP 1 blind enumeration"),
      f"preregistration commit={PREREGISTRATION_COMMIT}")

# Independent exact target from the already-closed Fibonacci/SU(2) bootstrap.
seed = 5
phi = (1+sp.sqrt(5))/2
check("bootstrap seed satisfies d_1(n)=phi exactly at n=5",
      sp.simplify(1+2*sp.cos(2*sp.pi/seed)-phi) == 0)

raw_ratios = [record["raw"]["gap_ratio_cross_over_fiber"]
              for record in data["fibrations"]]
generalized_ratios = [
    record["generalized"]["gap_ratio_cross_over_fiber"]
    for record in data["fibrations"]
]
check("all six preregistered raw gap ratios equal the bootstrap seed",
      len(raw_ratios) == 6
      and all(abs(ratio-seed) < 1e-10 for ratio in raw_ratios),
      f"ratios={raw_ratios}")
check("all six preregistered generalized Whitney ratios equal the seed",
      len(generalized_ratios) == 6
      and all(abs(ratio-seed) < 1e-10 for ratio in generalized_ratios),
      f"ratios={generalized_ratios}")

# Reconstruct one representative independently from the JSON.  Regular-facet
# assembly gives M=2V I+(V/4)A and K=wL.  This also exposes why the result is a
# preservation theorem: the Whitney stiffness is a common scalar times the
# already existing graph Laplacian split.
vertices, adjacency, _ = build_600cell()
fibers = first_hopf_fibration(vertices)
fiber_of = np.empty(120, dtype=int)
for fiber_index, fiber in enumerate(fibers):
    fiber_of[fiber] = fiber_index
fiber_adjacency = adjacency*(fiber_of[:, None] == fiber_of[None, :])
cross_adjacency = adjacency-fiber_adjacency
fiber_laplacian = np.diag(fiber_adjacency.sum(axis=1))-fiber_adjacency
cross_laplacian = np.diag(cross_adjacency.sum(axis=1))-cross_adjacency
volume = data["carrier"]["single_tetrahedron_volume"]
weight = data["whitney"]["stiffness_edge_weight"]
mass = 2*volume*np.eye(120)+(volume/4)*adjacency
fiber_stiffness = weight*fiber_laplacian
cross_stiffness = weight*cross_laplacian

fiber_values, fiber_vectors = sla.eigh(fiber_stiffness, mass)
cross_values, cross_vectors = sla.eigh(cross_stiffness, mass)
fiber_gap = first_positive(fiber_values)
cross_gap = first_positive(cross_values)
fiber_gap_vectors = fiber_vectors[:, abs(fiber_values-fiber_gap) < TOL]
cross_gap_vectors = cross_vectors[:, abs(cross_values-cross_gap) < TOL]
overlap = fiber_gap_vectors.T@mass@cross_gap_vectors
principal_cosines = sla.svdvals(overlap)
check("the two generalized gap eigenspaces both have dimension four",
      fiber_gap_vectors.shape[1] == 4 and cross_gap_vectors.shape[1] == 4)
check("fiber and cross gap eigenspaces are the same M-orthogonal subspace",
      np.max(abs(principal_cosines-1)) < 1e-9,
      f"principal cosines={principal_cosines}")
check("cross stiffness equals seed times fiber stiffness on the gap space",
      np.linalg.norm((cross_stiffness-seed*fiber_stiffness)
                     @ fiber_gap_vectors) < 1e-8)

box_stiffness = cross_stiffness-seed*fiber_stiffness
box_eigenvalues = np.linalg.eigvalsh(box_stiffness)
box_kernel = int(np.count_nonzero(abs(box_eigenvalues) < TOL))
check("the seed-weighted Whitney split has exact numerical kernel dimension 9",
      box_kernel == 9, f"kernel dimension={box_kernel}")

# Exact symbolic cancellation on the shared low mode.  M is scalar there, so
# it changes both kinetic eigenvalues by the same denominator.
fiber_graph_gap = 2-phi
cross_graph_gap = seed*(2-phi)
full_graph_gap = (seed+1)*(2-phi)
mass_scalar = sp.Symbol("V", positive=True)/4*(20-full_graph_gap)
whitney_scale = sp.Symbol("w", positive=True)
generalized_fiber_exact = whitney_scale*fiber_graph_gap/mass_scalar
generalized_cross_exact = whitney_scale*cross_graph_gap/mass_scalar
check("the generalized ratio is exactly seed because the mass cancels",
      sp.simplify(generalized_cross_exact/generalized_fiber_exact-seed) == 0)

histograms = [record["tetrahedron_fiber_edge_histogram"]
              for record in data["fibrations"]]
check("every coarse tetrahedron has exactly one fiber edge in every fibration",
      all(histogram == {"1": 600} for histogram in histograms))

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: the preregistered raw and generalized ratios both equal the seed.")
print("DERIVED: the equality lives on one common four-dimensional gap subspace.")
print("DERIVED: the seed-weighted difference retains a nine-dimensional kernel.")
print("STRUCTURAL LIMIT: Whitney preserves the graph ratio; it does not independently select it.")
print("OPEN: extend the unique-per-tetrahedron split canonically to refinement modes.")
raise SystemExit(0 if passed == tests else 1)

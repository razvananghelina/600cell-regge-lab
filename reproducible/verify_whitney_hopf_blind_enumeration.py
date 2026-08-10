#!/usr/bin/env python3
"""Blind enumeration of the Hopf-split Whitney scalar kinetic spectra.

This STEP-1 certificate intentionally contains no comparison with any
bootstrap integer or proposed physical speed.  It constructs the geometry,
all six discrete Hopf fibrations, the exact-geometric lowest-order Whitney
mass/stiffness matrices on the unrefined 600-cell boundary, and writes the
complete observed spectral data used by the later comparison.
"""

from collections import Counter, defaultdict
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as sla

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_hopf_blind_enumeration.json")
TOL = 1e-8
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


def all_hopf_fibrations(vertices):
    fibrations = []
    seen = set()
    generators = []
    for index, generator in enumerate(vertices):
        power = generator.copy()
        for order in range(2, 121):
            power = quat_mult(power, generator)
            if np.allclose(power, (1, 0, 0, 0), atol=1e-6):
                if order == 10:
                    generators.append(index)
                break

    for generator_index in generators:
        generator = vertices[generator_index]
        subgroup = []
        power = np.array((1.0, 0.0, 0.0, 0.0))
        for _ in range(10):
            index = vertex_index(vertices, power)
            if index < 0:
                subgroup = []
                break
            subgroup.append(index)
            power = quat_mult(power, generator)
        if len(set(subgroup)) != 10:
            continue

        assigned = np.full(len(vertices), -1, dtype=int)
        fibers = []
        valid = True
        for left in range(len(vertices)):
            if assigned[left] >= 0:
                continue
            fiber = []
            fiber_index = len(fibers)
            for right in subgroup:
                target = vertex_index(
                    vertices, quat_mult(vertices[left], vertices[right]))
                if target >= 0 and assigned[target] < 0:
                    assigned[target] = fiber_index
                    fiber.append(target)
            if len(fiber) != 10:
                valid = False
                break
            fibers.append(tuple(sorted(fiber)))
        if not valid or len(fibers) != 12:
            continue
        signature = tuple(sorted(fibers))
        if signature not in seen:
            seen.add(signature)
            fibrations.append(signature)
    return sorted(fibrations)


def distinct_spectrum(eigenvalues, tolerance=2e-8):
    clusters = []
    for value in sorted(map(float, eigenvalues)):
        if not clusters or abs(value-clusters[-1][0]) > tolerance:
            clusters.append([value, 1])
        else:
            old, count = clusters[-1]
            clusters[-1][0] = (old*count+value)/(count+1)
            clusters[-1][1] += 1
    # Generalized eigensolvers differ at roughly 1e-12 inside degenerate
    # conjugate spectra.  Ten decimals is stricter than the declared 2e-8
    # clustering tolerance and makes the committed conjugacy invariant
    # deterministic without manufacturing extra digits.
    return [{"value": round(value, 10), "multiplicity": count}
            for value, count in clusters]


def first_positive(eigenvalues):
    positive = np.asarray(eigenvalues)[np.asarray(eigenvalues) > TOL]
    return float(np.min(positive))


print("=" * 78)
print("BLIND HOPF-SPLIT WHITNEY SPECTRAL ENUMERATION -- NO TARGET COMPARISON")
print("=" * 78)

vertices, adjacency, graph_laplacian = build_600cell()
neighbors = [set(np.flatnonzero(adjacency[index] > 0.5))
             for index in range(120)]
tetrahedra = []
for i in range(120):
    for j in sorted(vertex for vertex in neighbors[i] if vertex > i):
        common_ij = neighbors[i] & neighbors[j]
        for k in sorted(vertex for vertex in common_ij if vertex > j):
            common_ijk = common_ij & neighbors[k]
            for ell in sorted(vertex for vertex in common_ijk if vertex > k):
                tetrahedra.append((i, j, k, ell))

check("600-cell carrier has f0=120, f1=720, f3=600",
      len(vertices) == 120 and int(adjacency.sum()/2) == 720
      and len(tetrahedra) == 600)

# Assemble consistent Whitney zero-form mass and stiffness in the inherited
# piecewise-Euclidean facet metric.  The tetrahedra live affinely in R4.
mass = np.zeros((120, 120), dtype=float)
stiffness = np.zeros((120, 120), dtype=float)
tetra_volumes = []
for tetrahedron in tetrahedra:
    points = vertices[list(tetrahedron)]
    edge_jacobian = np.column_stack(
        [points[index]-points[0] for index in range(1, 4)])
    gram = edge_jacobian.T@edge_jacobian
    volume = float(np.sqrt(np.linalg.det(gram))/6)
    tetra_volumes.append(volume)
    inverse_gram = np.linalg.inv(gram)
    gradients = np.zeros((4, 4), dtype=float)
    gradients[1:] = (edge_jacobian@inverse_gram).T
    gradients[0] = -gradients[1:].sum(axis=0)
    local_mass = volume*(np.ones((4, 4))+np.eye(4))/20
    local_stiffness = volume*(gradients@gradients.T)
    for local_left, global_left in enumerate(tetrahedron):
        for local_right, global_right in enumerate(tetrahedron):
            mass[global_left, global_right] += local_mass[local_left, local_right]
            stiffness[global_left, global_right] += local_stiffness[local_left, local_right]

check("all 600 Euclidean facet volumes agree",
      max(tetra_volumes)-min(tetra_volumes) < 1e-10,
      f"volume={np.mean(tetra_volumes):.12f}")
check("Whitney mass is positive definite and stiffness has one zero mode",
      np.linalg.eigvalsh(mass)[0] > 0
      and np.count_nonzero(np.abs(np.linalg.eigvalsh(stiffness)) < TOL) == 1)
check("assembled stiffness has graph-Laplacian support",
      np.max(np.abs(stiffness[np.logical_and(adjacency < 0.5,
                                             ~np.eye(120, dtype=bool))])) < TOL)

# Record, but do not interpret, the symmetry-forced proportionality constants.
offdiag_edges = stiffness[np.triu(adjacency > 0.5, 1)]
stiffness_edge_weight = float(-np.mean(offdiag_edges))
stiffness_residual = float(np.max(np.abs(
    stiffness-stiffness_edge_weight*graph_laplacian)))
check("Whitney stiffness is one scalar times the graph Laplacian",
      stiffness_residual < 1e-8,
      f"edge weight={stiffness_edge_weight:.12f}")

fibrations = all_hopf_fibrations(vertices)
check("enumeration contains exactly six distinct Hopf fibrations",
      len(fibrations) == 6)

records = []
for fibration_index, fibers in enumerate(fibrations):
    fiber_of = np.empty(120, dtype=int)
    for fiber_index, fiber in enumerate(fibers):
        fiber_of[list(fiber)] = fiber_index
    fiber_adjacency = adjacency*(fiber_of[:, None] == fiber_of[None, :])
    cross_adjacency = adjacency-fiber_adjacency
    fiber_laplacian = (np.diag(fiber_adjacency.sum(axis=1))
                       - fiber_adjacency)
    cross_laplacian = (np.diag(cross_adjacency.sum(axis=1))
                       - cross_adjacency)
    fiber_stiffness = stiffness_edge_weight*fiber_laplacian
    cross_stiffness = stiffness_edge_weight*cross_laplacian

    raw_fiber = np.linalg.eigvalsh(fiber_stiffness)
    raw_cross = np.linalg.eigvalsh(cross_stiffness)
    generalized_fiber = sla.eigh(fiber_stiffness, mass,
                                 eigvals_only=True)
    generalized_cross = sla.eigh(cross_stiffness, mass,
                                 eigvals_only=True)
    generalized_full = sla.eigh(stiffness, mass, eigvals_only=True)
    raw_gap_fiber = first_positive(raw_fiber)
    raw_gap_cross = first_positive(raw_cross)
    generalized_gap_fiber = first_positive(generalized_fiber)
    generalized_gap_cross = first_positive(generalized_cross)

    # Count the number of already-classified old fiber edges in each coarse
    # tetrahedron.  This is relevant to whether the split has a unique local
    # continuation under subdivision, but no preferred answer is assumed.
    histogram = Counter()
    for tetrahedron in tetrahedra:
        number = sum(
            fiber_adjacency[left, right] > 0.5
            for left, right in combinations(tetrahedron, 2)
        )
        histogram[int(number)] += 1

    record = {
        "fibration_index": fibration_index,
        "fiber_edges": int(fiber_adjacency.sum()/2),
        "cross_edges": int(cross_adjacency.sum()/2),
        "tetrahedron_fiber_edge_histogram": {
            str(key): histogram[key] for key in sorted(histogram)
        },
        "raw": {
            "fiber_zero_multiplicity": int(np.count_nonzero(abs(raw_fiber) < TOL)),
            "cross_zero_multiplicity": int(np.count_nonzero(abs(raw_cross) < TOL)),
            "fiber_gap": round(raw_gap_fiber, 12),
            "cross_gap": round(raw_gap_cross, 12),
            "gap_ratio_cross_over_fiber": round(raw_gap_cross/raw_gap_fiber, 12),
            "fiber_spectrum": distinct_spectrum(raw_fiber),
            "cross_spectrum": distinct_spectrum(raw_cross),
        },
        "generalized": {
            "fiber_zero_multiplicity": int(np.count_nonzero(abs(generalized_fiber) < TOL)),
            "cross_zero_multiplicity": int(np.count_nonzero(abs(generalized_cross) < TOL)),
            "fiber_gap": round(generalized_gap_fiber, 12),
            "cross_gap": round(generalized_gap_cross, 12),
            "gap_ratio_cross_over_fiber": round(
                generalized_gap_cross/generalized_gap_fiber, 12),
            "fiber_spectrum": distinct_spectrum(generalized_fiber),
            "cross_spectrum": distinct_spectrum(generalized_cross),
            "full_spectrum": distinct_spectrum(generalized_full),
        },
    }
    records.append(record)

check("each split has 120 fiber and 600 cross edges",
      all(record["fiber_edges"] == 120 and record["cross_edges"] == 600
          for record in records))
check("all six raw spectral enumerations are identical",
      all(record["raw"] == records[0]["raw"] for record in records[1:]))
check("all six generalized spectral enumerations are identical",
      all(record["generalized"] == records[0]["generalized"]
          for record in records[1:]))
check("all six local fiber-edge histograms are identical",
      all(record["tetrahedron_fiber_edge_histogram"]
          == records[0]["tetrahedron_fiber_edge_histogram"]
          for record in records[1:]))

payload = {
    "protocol": "STEP 1 blind enumeration; no bootstrap/speed target comparison",
    "carrier": {
        "vertices": 120,
        "edges": 720,
        "tetrahedra": 600,
        "total_piecewise_euclidean_volume": round(sum(tetra_volumes), 12),
        "single_tetrahedron_volume": round(float(np.mean(tetra_volumes)), 12),
    },
    "whitney": {
        "stiffness_edge_weight": round(stiffness_edge_weight, 12),
        "stiffness_graph_laplacian_residual": round(stiffness_residual, 15),
        "mass_spectrum": distinct_spectrum(np.linalg.eigvalsh(mass)),
    },
    "number_of_fibrations": len(fibrations),
    "fibrations": records,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("blind enumeration JSON was written deterministically", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"BLIND_OUTPUT: {OUTPUT.name}")
print("DERIVED: complete base-level Whitney/Hopf split spectra are enumerated.")
print("DERIVED: all six fibrations have identical blind spectral data.")
print("OPEN BY PROTOCOL: no bootstrap integer or speed comparison is made here.")
raise SystemExit(0 if passed == tests else 1)

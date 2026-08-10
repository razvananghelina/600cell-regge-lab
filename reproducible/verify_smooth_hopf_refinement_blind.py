#!/usr/bin/env python3
"""Blind smooth-Hopf orthogonal split on projected 600-cell refinements.

Definition commit is recorded in smooth_hopf_refinement_preregistration.md.
This script contains no bootstrap integer and makes no target comparison.  It
writes geometric audits, canonical-mode Ritz data and low spectra.
"""

from collections import defaultdict
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as la
import scipy.sparse as sparse
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("smooth_hopf_refinement_blind.json")
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


def quat_conjugate(quaternion):
    result = quaternion.copy()
    result[1:] *= -1
    return result


def vertex_index(vertices, quaternion):
    distances = np.linalg.norm(vertices-quaternion, axis=1)
    index = int(np.argmin(distances))
    return index if distances[index] < 1e-6 else -1


def first_hopf_axis(vertices, adjacency):
    """Return the axis of the first valid right-coset decagon fibration."""
    identity = np.array((1.0, 0.0, 0.0, 0.0))
    for generator in vertices:
        power = generator.copy()
        order = None
        for candidate in range(2, 121):
            power = quat_mult(power, generator)
            if np.allclose(power, identity, atol=1e-6):
                order = candidate
                break
        if order != 10 or np.linalg.norm(generator[1:]) < 1e-10:
            continue
        subgroup = []
        power = identity.copy()
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
        if not valid or len(fibers) != 12:
            continue
        fiber_of = np.empty(120, dtype=int)
        for fiber_index, fiber in enumerate(fibers):
            fiber_of[fiber] = fiber_index
        fiber_adjacency = adjacency*(fiber_of[:, None] == fiber_of[None, :])
        axis = generator[1:]/np.linalg.norm(generator[1:])
        return generator, axis, fibers, fiber_adjacency
    raise RuntimeError("no valid order-ten Hopf fibration found")


def tetra_geometry(points):
    edge_jacobian = np.column_stack(
        [points[index]-points[0] for index in range(1, 4)])
    gram = edge_jacobian.T@edge_jacobian
    inverse_gram = np.linalg.inv(gram)
    volume = float(np.sqrt(np.linalg.det(gram))/6)
    gradients = np.zeros((4, 4), dtype=float)
    gradients[1:] = (edge_jacobian@inverse_gram).T
    gradients[0] = -gradients[1:].sum(axis=0)
    tangent = edge_jacobian@inverse_gram@edge_jacobian.T
    return volume, gradients, tangent


def append_local(rows, cols, data, indices, local):
    for local_row, global_row in enumerate(indices):
        for local_col, global_col in enumerate(indices):
            value = float(local[local_row, local_col])
            if abs(value) > 1e-15:
                rows.append(global_row)
                cols.append(global_col)
                data.append(value)


def sparse_matrix(rows, cols, data, dimension):
    result = sparse.coo_matrix(
        (data, (rows, cols)), shape=(dimension, dimension)).tocsr()
    result.sum_duplicates()
    return result


def assemble(positions, tetrahedra, axis):
    rows = defaultdict(list)
    cols = defaultdict(list)
    data = defaultdict(list)
    projection_norms = []
    vertical_spectra = []
    horizontal_spectra = []
    max_edge = 0.0

    pure_axis = np.r_[0.0, axis]
    for tetrahedron in tetrahedra:
        points = positions[list(tetrahedron)]
        volume, gradients, tangent = tetra_geometry(points)
        centroid = points.sum(axis=0)
        centroid /= np.linalg.norm(centroid)
        hopf_vector = quat_mult(centroid, pure_axis)
        projected = tangent@hopf_vector
        projected_norm = float(np.linalg.norm(projected))
        if projected_norm < 1e-12:
            raise RuntimeError("Hopf field degenerates in an element tangent plane")
        unit_vertical = projected/projected_norm
        vertical_projector = np.outer(unit_vertical, unit_vertical)
        horizontal_projector = tangent-vertical_projector

        projection_norms.append(projected_norm)
        vertical_spectra.append(np.linalg.eigvalsh(vertical_projector))
        horizontal_spectra.append(np.linalg.eigvalsh(horizontal_projector))
        for left, right in combinations(range(4), 2):
            max_edge = max(max_edge, float(np.linalg.norm(
                points[left]-points[right])))

        local_mass = volume*(np.ones((4, 4))+np.eye(4))/20
        local_vertical = volume*(gradients@vertical_projector@gradients.T)
        local_horizontal = volume*(gradients@horizontal_projector@gradients.T)
        local_full = volume*(gradients@gradients.T)
        for name, local in (("mass", local_mass),
                            ("vertical", local_vertical),
                            ("horizontal", local_horizontal),
                            ("full", local_full)):
            append_local(rows[name], cols[name], data[name], tetrahedron, local)

    matrices = {
        name: sparse_matrix(rows[name], cols[name], data[name], len(positions))
        for name in ("mass", "vertical", "horizontal", "full")
    }
    split_residual = float(spla.norm(
        matrices["vertical"]+matrices["horizontal"]-matrices["full"])
        / max(spla.norm(matrices["full"]), 1e-30))
    return matrices, {
        "minimum_projected_hopf_norm": min(projection_norms),
        "maximum_projected_hopf_norm": max(projection_norms),
        "maximum_chord_length": max_edge,
        "split_relative_residual": split_residual,
        "vertical_projector_eigenvalues": np.mean(vertical_spectra, axis=0),
        "horizontal_projector_eigenvalues": np.mean(horizontal_spectra, axis=0),
    }


def hopf_base_coordinates(positions, axis):
    pure_axis = np.r_[0.0, axis]
    values = []
    for q in positions:
        image = quat_mult(quat_mult(q, pure_axis), quat_conjugate(q))
        values.append(image[1:])
    return np.asarray(values)


def ritz_values(basis, operator, mass):
    gram = np.asarray(basis.T@(mass@basis))
    form = np.asarray(basis.T@(operator@basis))
    gram = (gram+gram.T)/2
    form = (form+form.T)/2
    return la.eigvalsh(form, gram)


def low_spectrum(operator, mass, count=24):
    initial = np.linspace(1.0, 2.0, operator.shape[0])
    initial /= np.linalg.norm(initial)
    values = spla.eigsh(
        operator, k=min(count, operator.shape[0]-2), M=mass,
        sigma=-1e-8, which="LM", return_eigenvectors=False,
        tol=2e-10, maxiter=30000, v0=initial)
    values.sort()
    return values


def rounded(values):
    return [round(float(value), 12) for value in values]


print("=" * 78)
print("BLIND SMOOTH-HOPF ORTHOGONAL REFINEMENT")
print("=" * 78)

vertices, adjacency, _ = build_600cell()
neighbors = [set(np.flatnonzero(adjacency[index] > 0.5))
             for index in range(120)]
coarse_tetrahedra = []
for i in range(120):
    for j in sorted(vertex for vertex in neighbors[i] if vertex > i):
        common_ij = neighbors[i] & neighbors[j]
        for k in sorted(vertex for vertex in common_ij if vertex > j):
            common_ijk = common_ij & neighbors[k]
            for ell in sorted(vertex for vertex in common_ijk if vertex > k):
                coarse_tetrahedra.append((i, j, k, ell))

triangles = sorted({face for tetrahedron in coarse_tetrahedra
                    for face in combinations(tetrahedron, 3)})
edges = sorted({edge for tetrahedron in coarse_tetrahedra
                for edge in combinations(tetrahedron, 2)})
coarse_cells = ([(index,) for index in range(120)], edges,
                triangles, coarse_tetrahedra)
fine_vertex_cells = [cell for layer in coarse_cells for cell in layer]
fine_vertex_index = {cell: index for index, cell in enumerate(fine_vertex_cells)}
raw_fine_positions = np.array([
    vertices[list(cell)].mean(axis=0) for cell in fine_vertex_cells
])
fine_positions = raw_fine_positions/np.linalg.norm(
    raw_fine_positions, axis=1)[:, None]
fine_tetrahedra = []
for tetrahedron in coarse_tetrahedra:
    for ordering in permutations(tetrahedron):
        flag = (
            (ordering[0],),
            tuple(sorted(ordering[:2])),
            tuple(sorted(ordering[:3])),
            tetrahedron,
        )
        fine_tetrahedra.append(tuple(fine_vertex_index[cell] for cell in flag))

generator, axis, fibers, fiber_adjacency = first_hopf_axis(vertices, adjacency)
local_fiber_counts = []
for tetrahedron in coarse_tetrahedra:
    local_fiber_counts.append(sum(
        fiber_adjacency[left, right] > 0.5
        for left, right in combinations(tetrahedron, 2)))

check("600-cell f-vector data and first barycentric carrier are exact",
      (len(vertices), len(edges), len(triangles), len(coarse_tetrahedra),
       len(fine_positions), len(fine_tetrahedra))
      == (120, 720, 1200, 600, 2640, 14400))
check("chosen smooth axis extends a discrete fibration with one edge per tetrahedron",
      len(fibers) == 12 and set(local_fiber_counts) == {1})
check("all projected fine nodes lie on the unit three-sphere",
      np.max(abs(np.linalg.norm(fine_positions, axis=1)-1)) < 2e-15)

coarse_matrices, coarse_audit = assemble(vertices, coarse_tetrahedra, axis)
fine_matrices, fine_audit = assemble(fine_positions, fine_tetrahedra, axis)

for label, audit in (("coarse", coarse_audit), ("fine", fine_audit)):
    check(f"{label} orthogonal projectors have ranks one and two",
          np.allclose(audit["vertical_projector_eigenvalues"], (0, 0, 0, 1),
                      atol=2e-10)
          and np.allclose(audit["horizontal_projector_eigenvalues"], (0, 0, 1, 1),
                          atol=2e-10),
          f"V={rounded(audit['vertical_projector_eigenvalues'])}, "
          f"H={rounded(audit['horizontal_projector_eigenvalues'])}")
    check(f"{label} vertical plus horizontal equals full stiffness",
          audit["split_relative_residual"] < 2e-13,
          f"relative residual={audit['split_relative_residual']:.3e}")
    check(f"{label} projected Hopf field does not degenerate",
          audit["minimum_projected_hopf_norm"] > 0.8,
          f"norm range={audit['minimum_projected_hopf_norm']:.6f}.."
          f"{audit['maximum_projected_hopf_norm']:.6f}")

payload = {
    "protocol": "smooth Hopf definition committed before blind execution; no bootstrap target",
    "hopf_generator": rounded(generator),
    "hopf_axis": rounded(axis),
    "dimensions": {
        "coarse_vertices": len(vertices),
        "coarse_tetrahedra": len(coarse_tetrahedra),
        "fine_vertices": len(fine_positions),
        "fine_tetrahedra": len(fine_tetrahedra),
    },
    "levels": {},
}

for label, positions, matrices, audit in (
        ("coarse", vertices, coarse_matrices, coarse_audit),
        ("fine", fine_positions, fine_matrices, fine_audit)):
    charged_basis = positions
    base_basis = hopf_base_coordinates(positions, axis)
    level = {
        "geometry": {
            "maximum_chord_length": round(audit["maximum_chord_length"], 12),
            "minimum_projected_hopf_norm": round(
                audit["minimum_projected_hopf_norm"], 12),
            "maximum_projected_hopf_norm": round(
                audit["maximum_projected_hopf_norm"], 12),
            "split_relative_residual": round(
                audit["split_relative_residual"], 15),
            "vertical_projector_eigenvalues": rounded(
                audit["vertical_projector_eigenvalues"]),
            "horizontal_projector_eigenvalues": rounded(
                audit["horizontal_projector_eigenvalues"]),
        },
        "charged_coordinate_ritz": {},
        "base_pullback_ritz": {},
        "low_spectra": {},
    }
    for name in ("vertical", "horizontal", "full"):
        level["charged_coordinate_ritz"][name] = rounded(ritz_values(
            charged_basis, matrices[name], matrices["mass"]))
        level["base_pullback_ritz"][name] = rounded(ritz_values(
            base_basis, matrices[name], matrices["mass"]))
        level["low_spectra"][name] = rounded(low_spectrum(
            matrices[name], matrices["mass"]))
    payload["levels"][label] = level

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("blind artifact written without a bootstrap target",
      "bootstrap target" in payload["protocol"]
      and "a1" not in json.dumps(payload).lower())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"BLIND_ARTIFACT={OUTPUT}")
print("NO continuum target comparison performed in this script.")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Blind two-level projected red refinement of the smooth Hopf split.

The 1-to-8 rule, central-octahedron diagonal and recorded observables are fixed
in smooth_hopf_red_refinement_preregistration.md.  No continuum target or
bootstrap integer is compared here.
"""

from collections import Counter, defaultdict
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as la
import scipy.sparse as sparse
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("smooth_hopf_red_refinement_blind.json")
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
        if valid and len(fibers) == 12:
            return generator, generator[1:]/np.linalg.norm(generator[1:])
    raise RuntimeError("no valid order-ten Hopf axis found")


def coarse_tetrahedra(adjacency):
    neighbors = [set(np.flatnonzero(adjacency[index] > 0.5))
                 for index in range(120)]
    result = []
    for i in range(120):
        for j in sorted(vertex for vertex in neighbors[i] if vertex > i):
            common_ij = neighbors[i] & neighbors[j]
            for k in sorted(vertex for vertex in common_ij if vertex > j):
                common_ijk = common_ij & neighbors[k]
                for ell in sorted(vertex for vertex in common_ijk if vertex > k):
                    result.append((i, j, k, ell))
    return result


def mesh_edges(tetrahedra):
    return sorted({tuple(sorted(edge)) for tet in tetrahedra
                   for edge in combinations(tet, 2)})


def mesh_faces(tetrahedra):
    return [tuple(sorted(face)) for tet in tetrahedra
            for face in combinations(tet, 3)]


def red_refine(positions, tetrahedra):
    new_positions = [point.copy() for point in positions]
    midpoint_index = {}
    for edge in mesh_edges(tetrahedra):
        midpoint = positions[edge[0]]+positions[edge[1]]
        midpoint /= np.linalg.norm(midpoint)
        midpoint_index[edge] = len(new_positions)
        new_positions.append(midpoint)

    refined = []
    selected_diagonal_lengths = []
    for tetrahedron in tetrahedra:
        a, b, c, d = sorted(tetrahedron)

        def midpoint(left, right):
            return midpoint_index[tuple(sorted((left, right)))]

        ab, ac, ad = midpoint(a, b), midpoint(a, c), midpoint(a, d)
        bc, bd, cd = midpoint(b, c), midpoint(b, d), midpoint(c, d)
        refined.extend((
            (a, ab, ac, ad),
            (b, ab, bc, bd),
            (c, ac, bc, cd),
            (d, ad, bd, cd),
        ))

        candidates = (
            (ab, cd, (ac, ad, bd, bc)),
            (ac, bd, (ab, ad, cd, bc)),
            (ad, bc, (ab, ac, cd, bd)),
        )
        ranked = []
        for left, right, cycle in candidates:
            length = float(np.linalg.norm(
                new_positions[left]-new_positions[right]))
            ranked.append((length, min(left, right), max(left, right),
                           left, right, cycle))
        length, _, _, left, right, cycle = min(ranked)
        selected_diagonal_lengths.append(length)
        for index in range(4):
            refined.append((left, right, cycle[index], cycle[(index+1) % 4]))

    return np.asarray(new_positions), refined, selected_diagonal_lengths


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


def mean_ratio_quality(points, volume):
    edge_square_sum = sum(
        float(np.dot(points[left]-points[right], points[left]-points[right]))
        for left, right in combinations(range(4), 2))
    return float(12*(3*volume)**(2/3)/edge_square_sum)


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
    qualities = []
    max_edge = 0.0
    pure_axis = np.r_[0.0, axis]

    for tetrahedron in tetrahedra:
        points = positions[list(tetrahedron)]
        volume, gradients, tangent = tetra_geometry(points)
        qualities.append(mean_ratio_quality(points, volume))
        centroid = points.sum(axis=0)
        centroid /= np.linalg.norm(centroid)
        hopf_vector = quat_mult(centroid, pure_axis)
        projected = tangent@hopf_vector
        projected_norm = float(np.linalg.norm(projected))
        if projected_norm < 1e-12:
            raise RuntimeError("projected Hopf field degenerates")
        unit_vertical = projected/projected_norm
        vertical_projector = np.outer(unit_vertical, unit_vertical)
        horizontal_projector = tangent-vertical_projector
        projection_norms.append(projected_norm)
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
    qualities = np.asarray(qualities)
    return matrices, {
        "maximum_chord_length": max_edge,
        "minimum_quality": float(qualities.min()),
        "median_quality": float(np.median(qualities)),
        "maximum_quality": float(qualities.max()),
        "minimum_projected_hopf_norm": min(projection_norms),
        "maximum_projected_hopf_norm": max(projection_norms),
        "split_relative_residual": split_residual,
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
    return la.eigvalsh((form+form.T)/2, (gram+gram.T)/2)


def low_spectrum(operator, mass, count=32):
    initial = np.linspace(1.0, 2.0, operator.shape[0])
    initial /= np.linalg.norm(initial)
    values = spla.eigsh(
        operator, k=min(count, operator.shape[0]-2), M=mass,
        sigma=-1e-8, which="LM", return_eigenvectors=False,
        tol=2e-10, maxiter=40000, v0=initial)
    values.sort()
    return values


def rounded(values):
    return [round(float(value), 12) for value in values]


print("=" * 78)
print("BLIND SMOOTH-HOPF PROJECTED RED REFINEMENT")
print("=" * 78)

vertices, adjacency, _ = build_600cell()
level0_tetrahedra = coarse_tetrahedra(adjacency)
level1_positions, level1_tetrahedra, diagonal1 = red_refine(
    vertices, level0_tetrahedra)
level2_positions, level2_tetrahedra, diagonal2 = red_refine(
    level1_positions, level1_tetrahedra)
generator, axis = first_hopf_axis(vertices, adjacency)

levels = (
    ("level0", vertices, level0_tetrahedra),
    ("level1", level1_positions, level1_tetrahedra),
    ("level2", level2_positions, level2_tetrahedra),
)
expected_sizes = ((120, 600), (840, 4800), (6480, 38400))
check("projected red-refinement sizes match the preregistration",
      tuple((len(positions), len(tetrahedra))
            for _, positions, tetrahedra in levels) == expected_sizes)
check("all vertices remain on the unit three-sphere",
      all(np.max(abs(np.linalg.norm(positions, axis=1)-1)) < 3e-15
          for _, positions, _ in levels))

topology = {}
for label, positions, tetrahedra in levels:
    edges = mesh_edges(tetrahedra)
    face_counts = Counter(mesh_faces(tetrahedra))
    topology[label] = {
        "vertices": len(positions),
        "edges": len(edges),
        "faces": len(face_counts),
        "tetrahedra": len(tetrahedra),
        "euler_characteristic": (
            len(positions)-len(edges)+len(face_counts)-len(tetrahedra)),
        "face_incidence_values": sorted(set(face_counts.values())),
    }
check("all three levels are closed Euler-zero tetrahedral three-complexes",
      all(item["euler_characteristic"] == 0
          and item["face_incidence_values"] == [2]
          for item in topology.values()))

payload = {
    "protocol": "projected red refinement fixed before blind execution; no bootstrap target",
    "hopf_generator": rounded(generator),
    "hopf_axis": rounded(axis),
    "topology": topology,
    "selected_central_diagonal_lengths": {
        "level1_min_mean_max": rounded((min(diagonal1), np.mean(diagonal1), max(diagonal1))),
        "level2_min_mean_max": rounded((min(diagonal2), np.mean(diagonal2), max(diagonal2))),
    },
    "levels": {},
}

for label, positions, tetrahedra in levels:
    matrices, audit = assemble(positions, tetrahedra, axis)
    check(f"{label} split reconstructs the full stiffness",
          audit["split_relative_residual"] < 2e-13,
          f"residual={audit['split_relative_residual']:.3e}")
    check(f"{label} elements and projected Hopf field are nondegenerate",
          audit["minimum_quality"] > 0
          and audit["minimum_projected_hopf_norm"] > 0,
          f"q_min={audit['minimum_quality']:.6f}, "
          f"|X|_min={audit['minimum_projected_hopf_norm']:.6f}")

    charged_basis = positions
    base_basis = hopf_base_coordinates(positions, axis)
    result = {
        "geometry": {
            key: round(float(value), 15)
            for key, value in audit.items()
        },
        "charged_coordinate_ritz": {},
        "base_pullback_ritz": {},
        "low_spectra": {},
    }
    for name in ("vertical", "horizontal", "full"):
        result["charged_coordinate_ritz"][name] = rounded(ritz_values(
            charged_basis, matrices[name], matrices["mass"]))
        result["base_pullback_ritz"][name] = rounded(ritz_values(
            base_basis, matrices[name], matrices["mass"]))
        result["low_spectra"][name] = rounded(low_spectrum(
            matrices[name], matrices["mass"]))
    payload["levels"][label] = result

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("blind artifact contains no bootstrap target",
      "no bootstrap target" in payload["protocol"]
      and "a1" not in json.dumps(payload).lower())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"BLIND_ARTIFACT={OUTPUT}")
print("NO shape-threshold or continuum-mode comparison performed.")
raise SystemExit(0 if passed == tests else 1)


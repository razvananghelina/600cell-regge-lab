#!/usr/bin/env python3
"""Target-blind first-refinement Hopf/Whitney tensor enumeration.

The local rule is fixed before inspecting the refined gaps.  Every coarse
tetrahedron has exactly one preregistered fiber edge.  Its single equal-edge
stiffness contribution defines a positive rank-one tangent tensor Q_f; the
complement Q_c=P_tangent-Q_f is positive.  These tensors, not hand-labelled
new edges, are integrated on all 24 barycentric children.

No bootstrap integer or proposed speed is referenced by this verifier.
"""

from collections import defaultdict
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_hopf_refinement_blind.json")
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


def tetra_geometry(points):
    edge_jacobian = np.column_stack(
        [points[index]-points[0] for index in range(1, 4)])
    gram = edge_jacobian.T@edge_jacobian
    inverse_gram = np.linalg.inv(gram)
    volume = float(np.sqrt(np.linalg.det(gram))/6)
    gradients = np.zeros((4, 4), dtype=float)
    gradients[1:] = (edge_jacobian@inverse_gram).T
    gradients[0] = -gradients[1:].sum(axis=0)
    tangent_projector = edge_jacobian@inverse_gram@edge_jacobian.T
    return volume, gradients, tangent_projector


def append_local(rows, cols, data, indices, local):
    for local_row, global_row in enumerate(indices):
        for local_col, global_col in enumerate(indices):
            value = float(local[local_row, local_col])
            if abs(value) > 1e-15:
                rows.append(global_row)
                cols.append(global_col)
                data.append(value)


def low_generalized_spectrum(operator, mass, maximum_k=384):
    """Return all zero modes and the first observed positive cluster.

    A small negative shift orders true zeros before positive eigenvalues while
    avoiding a singular factorization.  k is doubled until a positive value
    appears, so the returned zero multiplicity is not silently truncated.
    """
    dimension = operator.shape[0]
    for requested in (32, 64, 128, 256, maximum_k):
        requested = min(requested, dimension-2)
        values = spla.eigsh(operator, k=requested, M=mass, sigma=-1e-7,
                            which="LM", return_eigenvectors=False,
                            tol=2e-10, maxiter=20000)
        values.sort()
        positive = values[values > TOL]
        if len(positive):
            gap = float(positive[0])
            return values, int(np.count_nonzero(abs(values) <= TOL)), gap
        if requested == dimension-2:
            break
    raise RuntimeError("positive gap not reached within preregistered eigen window")


def cluster(values, tolerance=3e-7):
    clusters = []
    for value in values:
        value = float(value)
        if not clusters or abs(value-clusters[-1][0]) > tolerance:
            clusters.append([value, 1])
        else:
            old, count = clusters[-1]
            clusters[-1][0] = (old*count+value)/(count+1)
            clusters[-1][1] += 1
    return [{"value": round(value, 10), "multiplicity_in_window": count}
            for value, count in clusters]


print("=" * 78)
print("BLIND FIRST-REFINEMENT HOPF/WHITNEY TENSOR ENUMERATION")
print("=" * 78)

vertices, adjacency, _ = build_600cell()
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

triangles = sorted({face for tetrahedron in tetrahedra
                    for face in combinations(tetrahedron, 3)})
edges = sorted({edge for tetrahedron in tetrahedra
                for edge in combinations(tetrahedron, 2)})
coarse_cells = ([(index,) for index in range(120)], edges,
                triangles, tetrahedra)
fine_vertex_cells = [cell for layer in coarse_cells for cell in layer]
fine_vertex_index = {cell: index for index, cell in enumerate(fine_vertex_cells)}
fine_positions = np.array([
    vertices[list(cell)].mean(axis=0) for cell in fine_vertex_cells
])

fine_tetrahedra = []
fine_parent = []
for parent_index, tetrahedron in enumerate(tetrahedra):
    for ordering in permutations(tetrahedron):
        flag = (
            (ordering[0],),
            tuple(sorted(ordering[:2])),
            tuple(sorted(ordering[:3])),
            tetrahedron,
        )
        fine_tetrahedra.append(tuple(fine_vertex_index[cell] for cell in flag))
        fine_parent.append(parent_index)

check("coarse and fine scalar carriers have 120 and 2640 vertices",
      len(vertices) == 120 and len(fine_positions) == 2640)
check("first barycentric level has exactly 14400 top tetrahedra",
      len(fine_tetrahedra) == 14400)

fibers = first_hopf_fibration(vertices)
fiber_of = np.empty(120, dtype=int)
for fiber_index, fiber in enumerate(fibers):
    fiber_of[fiber] = fiber_index
fiber_adjacency = adjacency*(fiber_of[:, None] == fiber_of[None, :])

parent_tensors = []
nonzero_fiber_eigenvalues = []
minimum_cross_eigenvalues = []
for tetrahedron in tetrahedra:
    local_fiber_edges = [
        (local_left, local_right)
        for local_left, local_right in combinations(range(4), 2)
        if fiber_adjacency[tetrahedron[local_left],
                           tetrahedron[local_right]] > 0.5
    ]
    if len(local_fiber_edges) != 1:
        raise RuntimeError("local fiber edge is not unique")
    points = vertices[list(tetrahedron)]
    volume, gradients, tangent = tetra_geometry(points)
    full_local = volume*(gradients@gradients.T)
    local_left, local_right = local_fiber_edges[0]
    edge_weight = -full_local[local_left, local_right]
    edge_vector = points[local_left]-points[local_right]
    fiber_tensor = (edge_weight/volume)*np.outer(edge_vector, edge_vector)
    cross_tensor = tangent-fiber_tensor
    parent_tensors.append((fiber_tensor, cross_tensor))
    fiber_eigenvalues = np.linalg.eigvalsh(fiber_tensor)
    cross_eigenvalues = np.linalg.eigvalsh(cross_tensor)
    nonzero_fiber_eigenvalues.append(fiber_eigenvalues[-1])
    minimum_cross_eigenvalues.append(cross_eigenvalues[1])

check("local fiber tensor is rank one with one common positive eigenvalue",
      min(nonzero_fiber_eigenvalues) > 0
      and max(nonzero_fiber_eigenvalues)-min(nonzero_fiber_eigenvalues) < 1e-9,
      f"eigenvalue={np.mean(nonzero_fiber_eigenvalues):.12f}")
check("local cross tensor is positive definite on every tangent space",
      min(minimum_cross_eigenvalues) > 1e-10,
      f"minimum tangent eigenvalue={min(minimum_cross_eigenvalues):.12f}")

coarse_rows = defaultdict(list)
coarse_cols = defaultdict(list)
coarse_data = defaultdict(list)
for parent_index, tetrahedron in enumerate(tetrahedra):
    points = vertices[list(tetrahedron)]
    volume, gradients, _ = tetra_geometry(points)
    local_mass = volume*(np.ones((4, 4))+np.eye(4))/20
    fiber_tensor, cross_tensor = parent_tensors[parent_index]
    local_fiber = volume*(gradients@fiber_tensor@gradients.T)
    local_cross = volume*(gradients@cross_tensor@gradients.T)
    append_local(coarse_rows["mass"], coarse_cols["mass"],
                 coarse_data["mass"], tetrahedron, local_mass)
    append_local(coarse_rows["fiber"], coarse_cols["fiber"],
                 coarse_data["fiber"], tetrahedron, local_fiber)
    append_local(coarse_rows["cross"], coarse_cols["cross"],
                 coarse_data["cross"], tetrahedron, local_cross)

fine_rows = defaultdict(list)
fine_cols = defaultdict(list)
fine_data = defaultdict(list)
for fine_index, fine_tetrahedron in enumerate(fine_tetrahedra):
    points = fine_positions[list(fine_tetrahedron)]
    volume, gradients, _ = tetra_geometry(points)
    local_mass = volume*(np.ones((4, 4))+np.eye(4))/20
    fiber_tensor, cross_tensor = parent_tensors[fine_parent[fine_index]]
    local_fiber = volume*(gradients@fiber_tensor@gradients.T)
    local_cross = volume*(gradients@cross_tensor@gradients.T)
    append_local(fine_rows["mass"], fine_cols["mass"], fine_data["mass"],
                 fine_tetrahedron, local_mass)
    append_local(fine_rows["fiber"], fine_cols["fiber"], fine_data["fiber"],
                 fine_tetrahedron, local_fiber)
    append_local(fine_rows["cross"], fine_cols["cross"], fine_data["cross"],
                 fine_tetrahedron, local_cross)

def matrix(rows, cols, data, dimension):
    result = sparse.coo_matrix((data, (rows, cols)),
                               shape=(dimension, dimension)).tocsr()
    result.sum_duplicates()
    return result


coarse_mass = matrix(coarse_rows["mass"], coarse_cols["mass"],
                     coarse_data["mass"], 120)
coarse_fiber = matrix(coarse_rows["fiber"], coarse_cols["fiber"],
                      coarse_data["fiber"], 120)
coarse_cross = matrix(coarse_rows["cross"], coarse_cols["cross"],
                      coarse_data["cross"], 120)
fine_mass = matrix(fine_rows["mass"], fine_cols["mass"], fine_data["mass"],
                   2640)
fine_fiber = matrix(fine_rows["fiber"], fine_cols["fiber"],
                    fine_data["fiber"], 2640)
fine_cross = matrix(fine_rows["cross"], fine_cols["cross"],
                    fine_data["cross"], 2640)

# Canonical Whitney zero-form inclusion: evaluate each coarse barycentric
# hat function at every barycentre (fine vertex).
p_rows, p_cols, p_data = [], [], []
for fine_index, cell in enumerate(fine_vertex_cells):
    for coarse_vertex in cell:
        p_rows.append(fine_index)
        p_cols.append(coarse_vertex)
        p_data.append(1/len(cell))
inclusion = sparse.coo_matrix((p_data, (p_rows, p_cols)),
                              shape=(2640, 120)).tocsr()

def compression_residual(fine_matrix, coarse_matrix):
    difference = inclusion.T@fine_matrix@inclusion-coarse_matrix
    denominator = max(spla.norm(coarse_matrix), 1e-30)
    return float(spla.norm(difference)/denominator)


mass_residual = compression_residual(fine_mass, coarse_mass)
fiber_residual = compression_residual(fine_fiber, coarse_fiber)
cross_residual = compression_residual(fine_cross, coarse_cross)
check("mass and both tensor forms compress to coarse level",
      max(mass_residual, fiber_residual, cross_residual) < 2e-9,
      f"relative residuals M/F/C={mass_residual:.3e}/"
      f"{fiber_residual:.3e}/{cross_residual:.3e}")

full_row_residual = float(np.max(abs(np.asarray(
    (fine_fiber+fine_cross).sum(axis=1)).ravel())))
check("refined kinetic sum has zero row sums",
      full_row_residual < 2e-9,
      f"max row residual={full_row_residual:.3e}")

coarse_fiber_values, coarse_fiber_zeros, coarse_fiber_gap = (
    low_generalized_spectrum(coarse_fiber, coarse_mass, maximum_k=128))
coarse_cross_values, coarse_cross_zeros, coarse_cross_gap = (
    low_generalized_spectrum(coarse_cross, coarse_mass, maximum_k=128))
fine_fiber_values, fine_fiber_zeros, fine_fiber_gap = (
    low_generalized_spectrum(fine_fiber, fine_mass))
fine_cross_values, fine_cross_zeros, fine_cross_gap = (
    low_generalized_spectrum(fine_cross, fine_mass))

payload = {
    "protocol": "target-blind tensor extension fixed before refined comparison",
    "definition": {
        "local_fiber_edges_per_parent": 1,
        "fiber_tensor_rank": 1,
        "fiber_tensor_positive_eigenvalue": round(
            float(np.mean(nonzero_fiber_eigenvalues)), 12),
        "cross_tensor_minimum_tangent_eigenvalue": round(
            float(min(minimum_cross_eigenvalues)), 12),
    },
    "dimensions": {"coarse": 120, "fine": 2640,
                   "fine_tetrahedra": 14400},
    "compression_relative_residuals": {
        "mass": mass_residual, "fiber": fiber_residual,
        "cross": cross_residual,
    },
    "coarse": {
        "fiber_kernel_in_window": coarse_fiber_zeros,
        "cross_kernel_in_window": coarse_cross_zeros,
        "fiber_gap": coarse_fiber_gap,
        "cross_gap": coarse_cross_gap,
        "gap_ratio_cross_over_fiber": coarse_cross_gap/coarse_fiber_gap,
        "fiber_low_spectrum": cluster(coarse_fiber_values),
        "cross_low_spectrum": cluster(coarse_cross_values),
    },
    "fine": {
        "fiber_kernel_in_window": fine_fiber_zeros,
        "cross_kernel_in_window": fine_cross_zeros,
        "fiber_gap": fine_fiber_gap,
        "cross_gap": fine_cross_gap,
        "gap_ratio_cross_over_fiber": fine_cross_gap/fine_fiber_gap,
        "fiber_low_spectrum": cluster(fine_fiber_values),
        "cross_low_spectrum": cluster(fine_cross_values),
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("blind refined enumeration JSON was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"BLIND_OUTPUT: {OUTPUT.name}")
print("DERIVED: the unique local edge defines positive fiber/cross tensors.")
print("DERIVED: both refined tensor forms are exactly Galerkin-compatible.")
print("OPEN BY PROTOCOL: refined gaps are recorded without a target comparison.")
raise SystemExit(0 if passed == tests else 1)

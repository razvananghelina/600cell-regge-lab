#!/usr/bin/env python3
"""Complete 600-cell exact Whitney trace-stiffness quotient spectra.

Protocol commit 7ed6d49 froze the carriers, quotient pencil, LOBPCG settings,
dense calibration, residual gates, outputs, and paired comparison before the
complete extremal spectra were evaluated.
"""

from itertools import combinations, permutations
import gc
import json
from math import factorial
from pathlib import Path
import sys
import warnings

import numpy as np
from scipy import linalg
import scipy.sparse as sparse
from scipy.sparse.linalg import LinearOperator, lobpcg
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_trace_stiffness_full.json")
CONTROL_CERTIFICATE = Path(__file__).with_name(
    "whitney_trace_stiffness.json"
)
UNWEIGHTED_CERTIFICATE = Path(__file__).with_name(
    "whitney_stiffness_refinement.json"
)
PROTOCOL_COMMIT = "7ed6d49"
EXPECTED_CONTROL_PROTOCOL = "b9a4104"
EXPECTED_UNWEIGHTED_PROTOCOL = "03e0abc"
RANDOM_SEED = 60_020_260_811
BLOCK_SIZE = 5
LOBPCG_TOLERANCE = 1e-9
LOBPCG_MAXITER = 2000
RITZ_RESIDUAL_GATE = 1e-7
CALIBRATION_VALUE_RELATIVE_GATE = 5e-7
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


coordinate_bases = [list(combinations(range(3), degree))
                    for degree in range(4)]


def wedge_components(covectors, degree):
    if degree == 0:
        return sy.Matrix((1,))
    return sy.Matrix([
        sy.det(sy.Matrix([
            [covector[index] for index in basis]
            for covector in covectors
        ]))
        for basis in coordinate_bases[degree]
    ])


def local_whitney_mass(points, degree):
    affine = sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    )
    inverse = affine.inv()
    gradients = [-sum(
        (sy.Matrix(inverse.row(row)).T for row in range(3)),
        sy.zeros(3, 1),
    )]
    gradients.extend(sy.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det()) / 6
    moment = volume * (sy.ones(4, 4) + sy.eye(4)) / 20
    forms = list(combinations(range(4), degree + 1))
    coefficients = []
    for form in forms:
        coefficient = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            coefficient[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                coefficient[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * wedge_components(covectors, degree)
                )
        coefficients.append(coefficient)
    mass = sy.zeros(len(forms), len(forms))
    for row, left in enumerate(coefficients):
        for column, right in enumerate(coefficients):
            mass[row, column] = sy.simplify(sum(
                (left[basis, :] * moment * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
            ))
    return mass


def triangle_gram(points):
    affine = sy.Matrix.hstack(points[1] - points[0], points[2] - points[0])
    return sy.simplify(affine.T * affine)


def triangle_whitney_mass_from_gram(gram, degree):
    gradients = (
        sy.Matrix((-1, -1)),
        sy.Matrix((1, 0)),
        sy.Matrix((0, 1)),
    )
    face_bases = [list(combinations(range(2), value)) for value in range(3)]

    def face_wedge(covectors):
        if degree == 0:
            return sy.Matrix((1,))
        return sy.Matrix([
            sy.det(sy.Matrix([
                [covector[index] for index in basis]
                for covector in covectors
            ]))
            for basis in face_bases[degree]
        ])

    inverse_gram = sy.simplify(gram.inv())
    wedge_metric = (
        sy.ones(1, 1) if degree == 0 else
        inverse_gram if degree == 1 else
        sy.Matrix(((sy.simplify(1 / gram.det()),),))
    )
    area = sy.sqrt(sy.simplify(gram.det())) / 2
    moment = area * (sy.ones(3, 3) + sy.eye(3)) / 12
    forms = list(combinations(range(3), degree + 1))
    coefficients = []
    for form in forms:
        coefficient = sy.zeros(len(face_bases[degree]), 3)
        if degree == 0:
            coefficient[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                coefficient[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * face_wedge(covectors)
                )
        coefficients.append(coefficient)
    mass = sy.zeros(len(forms), len(forms))
    for row, left in enumerate(coefficients):
        for column, right in enumerate(coefficients):
            mass[row, column] = sy.simplify(sum(
                wedge_metric[left_basis, right_basis]
                * (left[left_basis, :] * moment
                   * right[right_basis, :].T)[0]
                for left_basis in range(wedge_metric.rows)
                for right_basis in range(wedge_metric.cols)
            ))
    return mass


def local_coboundaries():
    faces = [list(combinations(range(4), degree + 1))
             for degree in range(4)]
    indices = [{face: index for index, face in enumerate(layer)}
               for layer in faces]
    differentials = []
    for degree in range(3):
        matrix = sy.zeros(len(faces[degree + 1]), len(faces[degree]))
        for row, simplex in enumerate(faces[degree + 1]):
            for omitted in range(degree + 2):
                face = simplex[:omitted] + simplex[omitted + 1:]
                matrix[row, indices[degree][face]] = (-1) ** omitted
        differentials.append(matrix)
    return differentials


def local_geometry(points, differentials):
    masses = [local_whitney_mass(points, degree) for degree in range(4)]
    offsets = np.cumsum((0, 4, 6, 4, 1))
    metric = sy.diag(*masses)
    weak = sy.zeros(15, 15)
    for degree, differential in enumerate(differentials):
        low_start, low_stop = offsets[degree:degree + 2]
        high_start, high_stop = offsets[degree + 1:degree + 3]
        forward = masses[degree + 1] * differential
        weak[high_start:high_stop, low_start:low_stop] = forward
        weak[low_start:low_stop, high_start:high_stop] = forward.T
    eigenvalues = linalg.eigvalsh(
        np.asarray(weak, dtype=np.float64),
        np.asarray(metric, dtype=np.float64),
    )
    return masses, float(np.max(np.abs(eigenvalues)))


def make_child_types(reference_vertices):
    types = []
    for ordering in permutations(range(4)):
        types.append((
            reference_vertices[ordering[0]],
            sum((reference_vertices[index] for index in ordering[:2]),
                sy.zeros(3, 1)) / 2,
            sum((reference_vertices[index] for index in ordering[:3]),
                sy.zeros(3, 1)) / 3,
            sum(reference_vertices, sy.zeros(3, 1)) / 4,
        ))
    return tuple(types)


def build_control_levels(reference_vertices, child_types):
    base_top = tuple(combinations(range(5), 4))
    base_cells = all_simplices(base_top)
    base_types = tuple(0 for _ in base_top)
    base_type_points = (reference_vertices,)

    fine_vertex_cells = tuple(cell for layer in base_cells for cell in layer)
    fine_vertex_index = {cell: index for index, cell in enumerate(
        fine_vertex_cells
    )}
    fine_top = []
    fine_types = []
    ordering_to_type = {
        ordering: index for index, ordering in enumerate(permutations(range(4)))
    }
    for top in base_top:
        position = {vertex: index for index, vertex in enumerate(top)}
        for ordering in permutations(top):
            flag = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                top,
            )
            fine_top.append(tuple(fine_vertex_index[cell] for cell in flag))
            local_ordering = tuple(position[vertex] for vertex in ordering)
            fine_types.append(ordering_to_type[local_ordering])
    fine_top = tuple(fine_top)
    return (
        ("base", base_top, base_cells, base_types, base_type_points),
        ("first_barycentric", fine_top, all_simplices(fine_top),
         tuple(fine_types), child_types),
    )


def build_full_levels(reference_vertices, child_types):
    vertices, adjacency, _ = build_600cell()
    neighbours = tuple(
        frozenset(np.flatnonzero(adjacency[index]).tolist())
        for index in range(120)
    )
    edges = tuple(
        (left, right)
        for left in range(120)
        for right in sorted(neighbours[left])
        if left < right
    )
    triangles = tuple(
        (left, right, third)
        for left, right in edges
        for third in sorted(neighbours[left] & neighbours[right])
        if right < third
    )
    top = tuple(
        (first, second, third, fourth)
        for first, second, third in triangles
        for fourth in sorted(
            neighbours[first] & neighbours[second] & neighbours[third]
        )
        if third < fourth
    )
    base_cells = (
        tuple((index,) for index in range(120)), edges, triangles, top
    )
    base_types = tuple(0 for _ in top)

    fine_vertex_cells = tuple(cell for layer in base_cells for cell in layer)
    fine_vertex_index = {cell: index for index, cell in enumerate(
        fine_vertex_cells
    )}
    fine_top = []
    fine_types = []
    ordering_list = tuple(permutations(range(4)))
    ordering_to_type = {
        ordering: index for index, ordering in enumerate(ordering_list)
    }
    for tetrahedron in top:
        position = {vertex: index for index, vertex in enumerate(tetrahedron)}
        for ordering in permutations(tetrahedron):
            flag = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                tetrahedron,
            )
            fine_top.append(tuple(fine_vertex_index[cell] for cell in flag))
            local_ordering = tuple(position[vertex] for vertex in ordering)
            fine_types.append(ordering_to_type[local_ordering])
    fine_top = tuple(fine_top)
    return (
        ("base", top, base_cells, base_types, (reference_vertices,)),
        ("first_barycentric", fine_top, all_simplices(fine_top),
         tuple(fine_types), child_types),
    )


class DisjointSet:
    def __init__(self, size):
        self.parent = np.arange(size, dtype=np.int32)

    def find(self, item):
        item = int(item)
        root = item
        while self.parent[root] != root:
            root = int(self.parent[root])
        while self.parent[item] != item:
            parent = int(self.parent[item])
            self.parent[item] = root
            item = parent
        return root

    def union(self, left, right):
        left = self.find(left)
        right = self.find(right)
        if left != right:
            self.parent[right] = left


def precompute_face_grams(type_points):
    result = {}
    for type_index, points in enumerate(type_points):
        for positions in combinations(range(4), 3):
            result[(type_index, positions)] = triangle_gram(tuple(
                points[position] for position in positions
            ))
    return result


def build_jump_quotient(top_cells, cells, top_types, type_points,
                        degree, local_mass):
    top_count = len(top_cells)
    local_faces = list(combinations(range(4), degree + 1))
    local_count = len(local_faces)
    cell_indices = {cell: index
                    for index, cell in enumerate(cells[degree])}
    triangle_indices = {cell: index
                        for index, cell in enumerate(cells[2])}

    lookup = {}
    triangle_parents = [[] for _ in cells[2]]
    for top_index, top in enumerate(top_cells):
        for local_index, positions in enumerate(local_faces):
            cell = tuple(top[position] for position in positions)
            lookup[(top_index, cell_indices[cell])] = (
                top_index * local_count + local_index
            )
        for triangle in combinations(top, 3):
            triangle_parents[triangle_indices[tuple(triangle)]].append(
                top_index
            )

    face_count = len(cells[2])
    face_local_count = len(list(combinations(range(3), degree + 1)))
    row_count = face_count * face_local_count
    r_rows = np.empty(2 * row_count, dtype=np.int32)
    r_columns = np.empty(2 * row_count, dtype=np.int32)
    r_data = np.empty(2 * row_count, dtype=np.float64)
    h_rows = np.empty(face_count * face_local_count ** 2, dtype=np.int32)
    h_columns = np.empty_like(h_rows)
    h_data = np.empty(face_count * face_local_count ** 2, dtype=np.float64)
    pairs_by_simplex = [[] for _ in cells[degree]]
    face_grams = precompute_face_grams(type_points)
    face_mass_cache = {}
    face_metric_mismatches = 0
    face_metric_types = set()

    r_cursor = 0
    h_cursor = 0
    for triangle_index, triangle in enumerate(cells[2]):
        parents = sorted(triangle_parents[triangle_index])
        if len(parents) != 2:
            raise AssertionError("not a closed 3-complex")
        grams = []
        for parent in parents:
            position = {
                vertex: index for index, vertex in enumerate(top_cells[parent])
            }
            positions = tuple(position[vertex] for vertex in triangle)
            grams.append(face_grams[(top_types[parent], positions)])
        if grams[0] != grams[1]:
            face_metric_mismatches += 1
        gram_key = tuple(sy.srepr(value) for value in grams[0])
        face_metric_types.add(gram_key)
        if gram_key not in face_mass_cache:
            exact_mass = triangle_whitney_mass_from_gram(grams[0], degree)
            face_mass_cache[gram_key] = np.asarray(
                exact_mass, dtype=np.float64
            )
        face_mass = face_mass_cache[gram_key]

        for local_row, simplex in enumerate(
            combinations(triangle, degree + 1)
        ):
            global_index = cell_indices[tuple(simplex)]
            left = lookup[(parents[0], global_index)]
            right = lookup[(parents[1], global_index)]
            row = triangle_index * face_local_count + local_row
            r_rows[r_cursor:r_cursor + 2] = row
            r_columns[r_cursor:r_cursor + 2] = (left, right)
            r_data[r_cursor:r_cursor + 2] = (1.0, -1.0)
            r_cursor += 2
            pairs_by_simplex[global_index].append((row, left, right))

        rows = np.arange(
            triangle_index * face_local_count,
            (triangle_index + 1) * face_local_count,
            dtype=np.int32,
        )
        block_size = face_local_count ** 2
        h_rows[h_cursor:h_cursor + block_size] = np.repeat(
            rows, face_local_count
        )
        h_columns[h_cursor:h_cursor + block_size] = np.tile(
            rows, face_local_count
        )
        h_data[h_cursor:h_cursor + block_size] = face_mass.ravel()
        h_cursor += block_size

    jump = sparse.csr_matrix(
        (r_data, (r_rows, r_columns)),
        shape=(row_count, top_count * local_count),
    )
    face_metric = sparse.csr_matrix(
        (h_data, (h_rows, h_columns)), shape=(row_count, row_count)
    )

    rank = sum(
        len(set(node for _, left, right in pairs for node in (left, right))) - 1
        for pairs in pairs_by_simplex
    )
    v_nonzeros = sum(
        len(pairs)
        * (len(set(node for _, left, right in pairs
                   for node in (left, right))) - 1)
        for pairs in pairs_by_simplex
    )
    v_rows = np.empty(v_nonzeros, dtype=np.int32)
    v_columns = np.empty(v_nonzeros, dtype=np.int32)
    v_data = np.empty(v_nonzeros, dtype=np.float64)
    rank_offset = 0
    v_cursor = 0
    all_connected = True
    maximum_basis_residual = 0.0
    for pairs in pairs_by_simplex:
        row_ids = [row for row, _, _ in pairs]
        nodes = sorted(set(
            node for _, left, right in pairs for node in (left, right)
        ))
        node_index = {node: index for index, node in enumerate(nodes)}
        incidence = np.zeros((len(pairs), len(nodes)), dtype=np.float64)
        dsu = DisjointSet(len(nodes))
        for local_row, (_, left, right) in enumerate(pairs):
            incidence[local_row, node_index[left]] = 1.0
            incidence[local_row, node_index[right]] = -1.0
            dsu.union(node_index[left], node_index[right])
        all_connected &= len({
            dsu.find(index) for index in range(len(nodes))
        }) == 1
        local_rank = len(nodes) - 1
        basis, _ = np.linalg.qr(incidence[:, :-1], mode="reduced")
        maximum_basis_residual = max(
            maximum_basis_residual,
            float(np.max(np.abs(basis.T @ basis - np.eye(local_rank))))
        )
        block_size = len(pairs) * local_rank
        v_rows[v_cursor:v_cursor + block_size] = np.repeat(
            np.asarray(row_ids, dtype=np.int32), local_rank
        )
        v_columns[v_cursor:v_cursor + block_size] = np.tile(
            np.arange(rank_offset, rank_offset + local_rank, dtype=np.int32),
            len(pairs),
        )
        v_data[v_cursor:v_cursor + block_size] = basis.ravel()
        v_cursor += block_size
        rank_offset += local_rank
    row_image_basis = sparse.csr_matrix(
        (v_data, (v_rows, v_columns)), shape=(row_count, rank)
    )
    return jump, face_metric, row_image_basis, {
        "degree": degree,
        "local_dimension": top_count * local_count,
        "global_dimension": len(cells[degree]),
        "constraint_rows": row_count,
        "constraint_rank": rank,
        "constraint_redundancy": row_count - rank,
        "row_image_basis_nonzeros": int(row_image_basis.nnz),
        "all_occurrence_graphs_connected": bool(all_connected),
        "maximum_basis_orthonormality_residual": maximum_basis_residual,
        "face_metric_mismatches": face_metric_mismatches,
        "face_metric_type_count": len(face_metric_types),
    }


def quotient_operators(jump, face_metric, basis, local_mass, top_count):
    local_inverse = np.linalg.inv(np.asarray(local_mass, dtype=np.float64))
    local_count = local_inverse.shape[0]

    def mass_solve(array):
        array = np.asarray(array)
        one_vector = array.ndim == 1
        if one_vector:
            array = array[:, None]
        reshaped = array.reshape(top_count, local_count, array.shape[1])
        solved = np.einsum("ij,tjk->tik", local_inverse, reshaped)
        result = solved.reshape(top_count * local_count, array.shape[1])
        return result[:, 0] if one_vector else result

    def apply_metric(array):
        return np.asarray(basis.T @ (face_metric @ (basis @ array)))

    def apply_energy(array):
        rows = basis @ array
        weighted_rows = face_metric @ rows
        copies = jump.T @ weighted_rows
        solved = mass_solve(copies)
        returned_rows = jump @ solved
        return np.asarray(basis.T @ (face_metric @ returned_rows))

    shape = (basis.shape[1], basis.shape[1])
    energy = LinearOperator(
        shape, matvec=apply_energy, matmat=apply_energy,
        rmatvec=apply_energy, dtype=np.float64,
    )
    metric = LinearOperator(
        shape, matvec=apply_metric, matmat=apply_metric,
        rmatvec=apply_metric, dtype=np.float64,
    )
    return energy, metric, apply_energy, apply_metric


def solve_extrema(jump, face_metric, basis, local_mass, top_count, seed):
    energy, metric, apply_energy, apply_metric = quotient_operators(
        jump, face_metric, basis, local_mass, top_count
    )
    rng = np.random.default_rng(seed)
    initial = rng.standard_normal((basis.shape[1], BLOCK_SIZE))
    warning_messages = []
    solutions = {}
    for label, largest in (("smallest", False), ("largest", True)):
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            values, vectors = lobpcg(
                energy, initial.copy(), B=metric,
                largest=largest, tol=LOBPCG_TOLERANCE,
                maxiter=LOBPCG_MAXITER,
            )
        warning_messages.extend(str(item.message) for item in captured)
        order = np.argsort(values)
        values = values[order]
        vectors = vectors[:, order]
        residuals = []
        for index, value in enumerate(values):
            vector = vectors[:, index]
            left = np.asarray(apply_energy(vector)).ravel()
            right = value * np.asarray(apply_metric(vector)).ravel()
            residuals.append(float(
                np.linalg.norm(left - right)
                / max(1.0, np.linalg.norm(left), np.linalg.norm(right))
            ))
        solutions[label] = {
            "values": values.tolist(),
            "maximum_relative_ritz_residual": max(residuals),
        }
    return {
        "five_smallest_positive_eigenvalues": solutions["smallest"]["values"],
        "five_largest_positive_eigenvalues": solutions["largest"]["values"],
        "positive_gap": solutions["smallest"]["values"][0],
        "maximum_positive_eigenvalue": solutions["largest"]["values"][-1],
        "maximum_relative_ritz_residual": max(
            solutions["smallest"]["maximum_relative_ritz_residual"],
            solutions["largest"]["maximum_relative_ritz_residual"],
        ),
        "solver_warning_count": len(warning_messages),
        "solver_warnings": warning_messages,
    }


def audit_level(level_index, level, local_masses, local_norm,
                target_records=None):
    name, top_cells, cells, top_types, type_points = level
    print(f"\n-- {name} --")
    records = []
    for degree in range(3):
        jump, face_metric, basis, structure = build_jump_quotient(
            top_cells, cells, top_types, type_points,
            degree, local_masses[degree],
        )
        spectral = solve_extrema(
            jump, face_metric, basis, local_masses[degree], len(top_cells),
            RANDOM_SEED + 100 * level_index + degree,
        )
        record = {**structure, **spectral}
        record["a_over_gap"] = local_norm / record["positive_gap"]
        if target_records is not None:
            target = target_records[degree]
            record["target_gap"] = target["positive_gap"]
            record["target_maximum"] = target["maximum_eigenvalue"]
            record["gap_relative_error"] = abs(
                record["positive_gap"] / target["positive_gap"] - 1.0
            )
            record["maximum_relative_error"] = abs(
                record["maximum_positive_eigenvalue"]
                / target["maximum_eigenvalue"] - 1.0
            )
        records.append(record)
        print(
            f"degree {degree}: gap={record['positive_gap']:.12g}, "
            f"max={record['maximum_positive_eigenvalue']:.12g}, "
            f"residual={record['maximum_relative_ritz_residual']:.3e}"
        )
        del jump, face_metric, basis
        gc.collect()
    return {
        "level": name,
        "f_vector": list(map(len, cells)),
        "top_count": len(top_cells),
        "local_dirac_norm": local_norm,
        "degree_records": records,
    }


print("=" * 78)
print("COMPLETE 600-CELL EXACT WHITNEY TRACE STIFFNESS")
print("=" * 78)

control_certificate = json.loads(CONTROL_CERTIFICATE.read_text())
unweighted_certificate = json.loads(UNWEIGHTED_CERTIFICATE.read_text())
check("the dense control and paired unweighted certificates have frozen protocols",
      control_certificate["protocol_commit"] == EXPECTED_CONTROL_PROTOCOL
      and unweighted_certificate["protocol_commit"]
      == EXPECTED_UNWEIGHTED_PROTOCOL)

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
child_types = make_child_types(reference_vertices)
element_grams = [
    sy.simplify(sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    ).T * sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    ))
    for points in child_types
]
check("all 24 ordered barycentric child element metrics agree exactly",
      all(gram == element_grams[0] for gram in element_grams))

local_d = local_coboundaries()
coarse_masses, coarse_norm = local_geometry(reference_vertices, local_d)
fine_masses, fine_norm = local_geometry(child_types[0], local_d)

# Dense-control calibration through the new quotient representation.
control_levels = build_control_levels(reference_vertices, child_types)
control_audits = []
for level_index, level in enumerate(control_levels):
    target_level = control_certificate["levels"][level_index]
    target_records = target_level["trace_records"]
    masses = coarse_masses if level_index == 0 else fine_masses
    norm = coarse_norm if level_index == 0 else fine_norm
    control_audits.append(audit_level(
        level_index, level, masses, norm, target_records=target_records
    ))
maximum_calibration_value_error = max(
    record[key]
    for audit in control_audits for record in audit["degree_records"]
    for key in ("gap_relative_error", "maximum_relative_error")
)
maximum_calibration_residual = max(
    record["maximum_relative_ritz_residual"]
    for audit in control_audits for record in audit["degree_records"]
)
calibration_passed = (
    maximum_calibration_value_error < CALIBRATION_VALUE_RELATIVE_GATE
    and maximum_calibration_residual < RITZ_RESIDUAL_GATE
)
check("the quotient LOBPCG reproduces all dense control gaps and maxima",
      calibration_passed,
      f"max value error={maximum_calibration_value_error:.3e}, "
      f"max residual={maximum_calibration_residual:.3e}")
if not calibration_passed:
    raise SystemExit(1)

# Complete carriers.
full_levels = build_full_levels(reference_vertices, child_types)
check("the complete base and refined f-vectors are exact",
      full_levels[0][2] and list(map(len, full_levels[0][2]))
      == [120, 720, 1200, 600]
      and list(map(len, full_levels[1][2]))
      == [2640, 17040, 28800, 14400])
full_audits = []
for level_index, level in enumerate(full_levels):
    masses = coarse_masses if level_index == 0 else fine_masses
    norm = coarse_norm if level_index == 0 else fine_norm
    full_audits.append(audit_level(
        level_index, level, masses, norm
    ))

expected_rows = {
    "base": [3600, 3600, 1200],
    "first_barycentric": [86400, 86400, 28800],
}
expected_ranks = {
    "base": [2280, 2880, 1200],
    "first_barycentric": [54960, 69360, 28800],
}
check("all complete row and exact quotient-rank counts are reproduced",
      all(
          [record["constraint_rows"] for record in audit["degree_records"]]
          == expected_rows[audit["level"]]
          and [record["constraint_rank"] for record in audit["degree_records"]]
          == expected_ranks[audit["level"]]
          for audit in full_audits
      ))
check("all complete occurrence graphs are connected and face metrics agree",
      all(
          record["all_occurrence_graphs_connected"]
          and record["face_metric_mismatches"] == 0
          for audit in full_audits for record in audit["degree_records"]
      ))
maximum_basis_residual = max(
    record["maximum_basis_orthonormality_residual"]
    for audit in full_audits for record in audit["degree_records"]
)
check("all complete row-image bases remain orthonormal",
      maximum_basis_residual < 1e-11,
      f"maximum residual={maximum_basis_residual:.3e}")
maximum_full_ritz_residual = max(
    record["maximum_relative_ritz_residual"]
    for audit in full_audits for record in audit["degree_records"]
)
check("all complete extremal blocks meet the frozen Ritz residual gate",
      maximum_full_ritz_residual < RITZ_RESIDUAL_GATE,
      f"maximum residual={maximum_full_ritz_residual:.3e}")

trace_ratios = [
    full_audits[1]["degree_records"][degree]["a_over_gap"]
    / full_audits[0]["degree_records"][degree]["a_over_gap"]
    for degree in range(3)
]
trace_spread = max(trace_ratios) / min(trace_ratios)
unweighted_ratios = unweighted_certificate[
    "refinement_diagnostics"
]["degreewise_scale_ratios"]
unweighted_spread = max(unweighted_ratios) / min(unweighted_ratios)
exact_compatible = all(
    abs(value - 1.0) < CALIBRATION_VALUE_RELATIVE_GATE
    for value in trace_ratios
)
improved_balance = trace_spread < unweighted_spread
check("the complete compatibility label follows the frozen exact criterion",
      exact_compatible == all(
          abs(value - 1.0) < CALIBRATION_VALUE_RELATIVE_GATE
          for value in trace_ratios
      ))
check("the complete balance label follows the paired unweighted comparison",
      improved_balance == (trace_spread < unweighted_spread))

compatibility_verdict = (
    "DERIVED NUMERICAL: exact trace stiffness is first-step compatible on the complete 600-cell"
    if exact_compatible else
    "DERIVED NUMERICAL NEGATIVE: exact trace stiffness is not first-step compatible on the complete 600-cell"
)
balance_verdict = (
    "PATTERN: exact trace weighting improves complete degree balance"
    if improved_balance else
    "DERIVED NEGATIVE: exact trace weighting does not improve complete degree balance"
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "control_protocol_commit": EXPECTED_CONTROL_PROTOCOL,
    "unweighted_protocol_commit": EXPECTED_UNWEIGHTED_PROTOCOL,
    "phenomenological_target_used": False,
    "candidate_count": 1,
    "solver": {
        "method": "generalized block LOBPCG on exact row-image quotient",
        "seed": RANDOM_SEED,
        "block_size": BLOCK_SIZE,
        "tolerance": LOBPCG_TOLERANCE,
        "maximum_iterations": LOBPCG_MAXITER,
        "ritz_residual_gate": RITZ_RESIDUAL_GATE,
    },
    "calibration": {
        "audits": control_audits,
        "maximum_value_relative_error": maximum_calibration_value_error,
        "maximum_ritz_relative_residual": maximum_calibration_residual,
    },
    "full_audits": full_audits,
    "refinement_comparison": {
        "trace_degreewise_ratios": trace_ratios,
        "trace_ratio_spread": trace_spread,
        "unweighted_degreewise_ratios": unweighted_ratios,
        "unweighted_ratio_spread": unweighted_spread,
        "exact_first_step_compatible": exact_compatible,
        "degree_balance_improved": improved_balance,
    },
    "verdicts": [
        compatibility_verdict,
        balance_verdict,
        "OPEN: repeated-refinement law and overall dimensionless stiffness",
        "OPEN: chiral finite-stiffness realization and causal dynamics",
    ],
    "scope": (
        "Complete 600-cell at levels zero and one; numerical quotient "
        "extremal spectra calibrated on complete dense controls."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured complete trace-stiffness certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("TRACE_RATIOS=" + str([float(value) for value in trace_ratios]))
print(f"TRACE_SPREAD={trace_spread:.12g}")
print(f"UNWEIGHTED_SPREAD={unweighted_spread:.12g}")
print("COMPATIBILITY_VERDICT: " + compatibility_verdict)
print("BALANCE_VERDICT: " + balance_verdict)
raise SystemExit(0 if passed == tests else 1)

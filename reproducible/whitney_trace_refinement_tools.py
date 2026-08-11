"""Shared exact-geometry tools for iterated Whitney trace refinement tests."""

from itertools import combinations, permutations
from math import factorial
import warnings

import numpy as np
from scipy import linalg
import scipy.sparse as sparse
from scipy.sparse.linalg import LinearOperator, eigsh, lobpcg
import sympy as sy


BLOCK_SIZE = 5
LOBPCG_TOLERANCE = 1e-9
LOBPCG_MAXITER = 2000
LANCZOS_TOLERANCE = 1e-11
LANCZOS_MAXITER = 20_000


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


def tetrahedron_gram(points):
    affine = sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    )
    return sy.simplify(affine.T * affine)


def matrix_key(matrix):
    return tuple(sy.srepr(value) for value in matrix)


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


LOCAL_D = local_coboundaries()


def local_geometry(points):
    masses = [local_whitney_mass(points, degree) for degree in range(4)]
    offsets = np.cumsum((0, 4, 6, 4, 1))
    metric = sy.diag(*masses)
    weak = sy.zeros(15, 15)
    for degree, differential in enumerate(LOCAL_D):
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


def make_base_level(reference_vertices):
    top = tuple(combinations(range(5), 4))
    return {
        "level": 0,
        "name": "base",
        "top": top,
        "cells": all_simplices(top),
        "top_points": tuple(reference_vertices for _ in top),
    }


def barycentric_refine(level):
    cells = level["cells"]
    fine_vertex_cells = tuple(cell for layer in cells for cell in layer)
    fine_vertex_index = {cell: index for index, cell in enumerate(
        fine_vertex_cells
    )}
    fine_top = []
    fine_top_points = []
    for top, points in zip(level["top"], level["top_points"]):
        for ordering in permutations(range(4)):
            ordered_vertices = tuple(top[index] for index in ordering)
            flag = (
                (ordered_vertices[0],),
                tuple(sorted(ordered_vertices[:2])),
                tuple(sorted(ordered_vertices[:3])),
                top,
            )
            fine_top.append(tuple(fine_vertex_index[cell] for cell in flag))
            fine_top_points.append((
                points[ordering[0]],
                (points[ordering[0]] + points[ordering[1]]) / 2,
                (points[ordering[0]] + points[ordering[1]]
                 + points[ordering[2]]) / 3,
                sum(points, sy.zeros(3, 1)) / 4,
            ))
    fine_top = tuple(fine_top)
    next_level = level["level"] + 1
    return {
        "level": next_level,
        "name": f"barycentric_{next_level}",
        "top": fine_top,
        "cells": all_simplices(fine_top),
        "top_points": tuple(fine_top_points),
    }


def classify_element_types(level):
    key_to_type = {}
    representatives = []
    top_types = []
    for points in level["top_points"]:
        key = matrix_key(tetrahedron_gram(points))
        if key not in key_to_type:
            key_to_type[key] = len(representatives)
            representatives.append(points)
        top_types.append(key_to_type[key])
    geometries = [local_geometry(points) for points in representatives]
    masses = [geometry[0] for geometry in geometries]
    norms = [geometry[1] for geometry in geometries]
    return {
        "top_types": np.asarray(top_types, dtype=np.int32),
        "representatives": tuple(representatives),
        "masses": masses,
        "norms": norms,
        "maximum_dirac_norm": max(norms),
        "type_count": len(representatives),
    }


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


def build_jump_quotient(level, element_types, degree):
    top_cells = level["top"]
    cells = level["cells"]
    top_count = len(top_cells)
    top_types = element_types["top_types"]
    representatives = element_types["representatives"]
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

    face_grams = {}
    for type_index, points in enumerate(representatives):
        for positions in combinations(range(4), 3):
            face_grams[(type_index, positions)] = triangle_gram(tuple(
                points[position] for position in positions
            ))

    face_count = len(cells[2])
    face_local_count = len(list(combinations(range(3), degree + 1)))
    row_count = face_count * face_local_count
    r_rows = np.empty(2 * row_count, dtype=np.int32)
    r_columns = np.empty(2 * row_count, dtype=np.int32)
    r_data = np.empty(2 * row_count, dtype=np.float64)
    h_size = face_count * face_local_count ** 2
    h_rows = np.empty(h_size, dtype=np.int32)
    h_columns = np.empty(h_size, dtype=np.int32)
    h_data = np.empty(h_size, dtype=np.float64)
    pairs_by_simplex = [[] for _ in cells[degree]]
    face_mass_cache = {}
    face_metric_mismatches = 0
    face_metric_types = set()
    r_cursor = h_cursor = 0

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
            grams.append(face_grams[(int(top_types[parent]), positions)])
        if grams[0] != grams[1]:
            face_metric_mismatches += 1
        gram_key = matrix_key(grams[0])
        face_metric_types.add(gram_key)
        if gram_key not in face_mass_cache:
            face_mass_cache[gram_key] = np.asarray(
                triangle_whitney_mass_from_gram(grams[0], degree),
                dtype=np.float64,
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
    rank_offset = v_cursor = 0
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
    inverse_by_type = [
        np.linalg.inv(np.asarray(masses[degree], dtype=np.float64))
        for masses in element_types["masses"]
    ]
    inverse_by_top = np.asarray([
        inverse_by_type[int(type_index)] for type_index in top_types
    ])
    return jump, face_metric, row_image_basis, inverse_by_top, {
        "degree": degree,
        "local_dimension": top_count * local_count,
        "global_dimension": len(cells[degree]),
        "constraint_rows": row_count,
        "constraint_rank": rank,
        "constraint_redundancy": row_count - rank,
        "all_occurrence_graphs_connected": bool(all_connected),
        "maximum_basis_orthonormality_residual": maximum_basis_residual,
        "face_metric_mismatches": face_metric_mismatches,
        "face_metric_type_count": len(face_metric_types),
        "row_image_basis_nonzeros": int(row_image_basis.nnz),
    }


def mass_solve(array, inverse_by_top):
    array = np.asarray(array)
    one_vector = array.ndim == 1
    if one_vector:
        array = array[:, None]
    top_count, local_count, _ = inverse_by_top.shape
    reshaped = array.reshape(top_count, local_count, array.shape[1])
    solved = np.einsum("tij,tjk->tik", inverse_by_top, reshaped)
    result = solved.reshape(top_count * local_count, array.shape[1])
    return result[:, 0] if one_vector else result


def quotient_operators(jump, face_metric, basis, inverse_by_top):
    def apply_metric(array):
        return np.asarray(basis.T @ (face_metric @ (basis @ array)))

    def apply_energy(array):
        rows = basis @ array
        weighted = face_metric @ rows
        copies = jump.T @ weighted
        solved = mass_solve(copies, inverse_by_top)
        returned = jump @ solved
        return np.asarray(basis.T @ (face_metric @ returned))

    shape = (basis.shape[1], basis.shape[1])
    return (
        LinearOperator(
            shape, matvec=apply_energy, matmat=apply_energy,
            rmatvec=apply_energy, dtype=np.float64,
        ),
        LinearOperator(
            shape, matvec=apply_metric, matmat=apply_metric,
            rmatvec=apply_metric, dtype=np.float64,
        ),
        apply_energy,
        apply_metric,
    )


def solve_generalized(jump, face_metric, basis, inverse_by_top, seed):
    energy, metric, apply_energy, apply_metric = quotient_operators(
        jump, face_metric, basis, inverse_by_top
    )
    rng = np.random.default_rng(seed)
    initial = rng.standard_normal((basis.shape[1], BLOCK_SIZE))
    warnings_seen = []
    solutions = {}
    for label, largest in (("smallest", False), ("largest", True)):
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            values, vectors = lobpcg(
                energy, initial.copy(), B=metric, largest=largest,
                tol=LOBPCG_TOLERANCE, maxiter=LOBPCG_MAXITER,
            )
        warnings_seen.extend(str(item.message) for item in captured)
        order = np.argsort(values)
        values = values[order]
        vectors = vectors[:, order]
        residuals = []
        for index, value in enumerate(values):
            left = np.asarray(apply_energy(vectors[:, index])).ravel()
            right = value * np.asarray(
                apply_metric(vectors[:, index])
            ).ravel()
            residuals.append(float(
                np.linalg.norm(left - right)
                / max(1.0, np.linalg.norm(left), np.linalg.norm(right))
            ))
        solutions[label] = (values, max(residuals))
    return {
        "five_smallest_positive_eigenvalues": solutions["smallest"][0].tolist(),
        "five_largest_positive_eigenvalues": solutions["largest"][0].tolist(),
        "positive_gap": float(solutions["smallest"][0][0]),
        "maximum_positive_eigenvalue": float(solutions["largest"][0][-1]),
        "maximum_relative_ritz_residual": max(
            solutions["smallest"][1], solutions["largest"][1]
        ),
        "solver_method": "generalized_block_lobpcg",
        "solver_warning_count": len(warnings_seen),
        "solver_warnings": warnings_seen,
    }


def solve_degree_two(jump, face_metric, inverse_by_top, seed):
    weights = np.asarray(face_metric.diagonal(), dtype=np.float64)
    square_root = np.sqrt(weights)

    def apply(array):
        array = np.asarray(array)
        one_vector = array.ndim == 1
        if one_vector:
            array = array[:, None]
        weighted = square_root[:, None] * array
        copies = jump.T @ weighted
        solved = mass_solve(copies, inverse_by_top)
        returned = np.asarray(jump @ solved)
        result = square_root[:, None] * returned
        return result[:, 0] if one_vector else result

    operator = LinearOperator(
        (jump.shape[0], jump.shape[0]),
        matvec=apply, matmat=apply, rmatvec=apply, dtype=np.float64,
    )
    rng = np.random.default_rng(seed)
    initial = rng.standard_normal(jump.shape[0])
    initial /= np.linalg.norm(initial)
    solutions = {}
    for label, which in (("smallest", "SA"), ("largest", "LA")):
        values, vectors = eigsh(
            operator, k=BLOCK_SIZE, which=which, v0=initial,
            tol=LANCZOS_TOLERANCE, maxiter=LANCZOS_MAXITER,
        )
        order = np.argsort(values)
        values = values[order]
        vectors = vectors[:, order]
        residuals = []
        for index, value in enumerate(values):
            left = np.asarray(apply(vectors[:, index])).ravel()
            right = value * vectors[:, index]
            residuals.append(float(
                np.linalg.norm(left - right)
                / max(1.0, np.linalg.norm(left), np.linalg.norm(right))
            ))
        solutions[label] = (values, max(residuals))
    return {
        "five_smallest_positive_eigenvalues": solutions["smallest"][0].tolist(),
        "five_largest_positive_eigenvalues": solutions["largest"][0].tolist(),
        "positive_gap": float(solutions["smallest"][0][0]),
        "maximum_positive_eigenvalue": float(solutions["largest"][0][-1]),
        "maximum_relative_ritz_residual": max(
            solutions["smallest"][1], solutions["largest"][1]
        ),
        "solver_method": "symmetric_lanczos_full_row_degree_two",
        "solver_warning_count": 0,
        "solver_warnings": [],
    }


def solve_explicit_generalized(jump, face_metric, basis, inverse_by_top, seed):
    """Explicit sparse quotient frozen after the level-two LOBPCG failure."""
    inverse_blocks = [sparse.csr_matrix(block) for block in inverse_by_top]
    mass_inverse = sparse.block_diag(inverse_blocks, format="csr")
    weighted_basis = (face_metric @ basis).tocsr()
    pullback = (jump.T @ weighted_basis).tocsr()
    metric = (basis.T @ weighted_basis).tocsr()
    energy = (pullback.T @ (mass_inverse @ pullback)).tocsr()
    metric.eliminate_zeros()
    energy.eliminate_zeros()

    def symmetry_residual(matrix):
        difference = (matrix - matrix.T).tocsr()
        if difference.nnz == 0:
            return 0.0
        return float(np.max(np.abs(difference.data)))

    energy_symmetry_residual = symmetry_residual(energy)
    metric_symmetry_residual = symmetry_residual(metric)
    # Sparse multiplication may sum identical floating terms in different
    # orders.  After recording the residual, remove that roundoff asymmetry.
    energy = ((energy + energy.T) * 0.5).tocsr()
    metric = ((metric + metric.T) * 0.5).tocsr()

    rng = np.random.default_rng(seed)
    initial = rng.standard_normal(energy.shape[0])
    initial /= np.linalg.norm(initial)
    solutions = {}
    for label, options in (
        ("smallest", {"sigma": 0.0, "which": "LM"}),
        ("largest", {"which": "LA"}),
    ):
        values, vectors = eigsh(
            energy, k=BLOCK_SIZE, M=metric, v0=initial,
            tol=LANCZOS_TOLERANCE, maxiter=LANCZOS_MAXITER,
            **options,
        )
        order = np.argsort(values)
        values = values[order]
        vectors = vectors[:, order]
        residuals = []
        for index, value in enumerate(values):
            left = np.asarray(energy @ vectors[:, index]).ravel()
            right = value * np.asarray(
                metric @ vectors[:, index]
            ).ravel()
            residuals.append(float(
                np.linalg.norm(left - right)
                / max(1.0, np.linalg.norm(left), np.linalg.norm(right))
            ))
        solutions[label] = (values, max(residuals))
    return {
        "five_smallest_positive_eigenvalues": solutions["smallest"][0].tolist(),
        "five_largest_positive_eigenvalues": solutions["largest"][0].tolist(),
        "positive_gap": float(solutions["smallest"][0][0]),
        "maximum_positive_eigenvalue": float(solutions["largest"][0][-1]),
        "maximum_relative_ritz_residual": max(
            solutions["smallest"][1], solutions["largest"][1]
        ),
        "solver_method": "explicit_sparse_generalized_lanczos",
        "solver_warning_count": 0,
        "solver_warnings": [],
        "explicit_energy_nonzeros": int(energy.nnz),
        "explicit_metric_nonzeros": int(metric.nnz),
        "energy_symmetry_residual_before_symmetrization": (
            energy_symmetry_residual
        ),
        "metric_symmetry_residual_before_symmetrization": (
            metric_symmetry_residual
        ),
    }


def audit_level(level, seed_base):
    element_types = classify_element_types(level)
    records = []
    for degree in range(3):
        jump, face_metric, basis, inverse_by_top, structure = (
            build_jump_quotient(level, element_types, degree)
        )
        if degree == 2:
            spectral = solve_degree_two(
                jump, face_metric, inverse_by_top, seed_base + degree
            )
        else:
            spectral = solve_explicit_generalized(
                jump, face_metric, basis, inverse_by_top, seed_base + degree
            )
        record = {**structure, **spectral}
        record["a_over_gap"] = (
            element_types["maximum_dirac_norm"] / record["positive_gap"]
        )
        records.append(record)
        print(
            f"degree {degree}: gap={record['positive_gap']:.12g}, "
            f"max={record['maximum_positive_eigenvalue']:.12g}, "
            f"residual={record['maximum_relative_ritz_residual']:.3e}"
        )
    return {
        "level": level["level"],
        "name": level["name"],
        "f_vector": list(map(len, level["cells"])),
        "top_count": len(level["top"]),
        "duplicated_dimension": 15 * len(level["top"]),
        "element_metric_type_count": element_types["type_count"],
        "element_type_dirac_norms": element_types["norms"],
        "maximum_local_dirac_norm": element_types["maximum_dirac_norm"],
        "degree_records": records,
    }

#!/usr/bin/env python3
"""Exact local constrained-pencil representation of Whitney spectral data.

Protocol commit ba6035e froze the element carrier, canonical equality
constraint, exact pencil equivalence, small calibration and dynamical scope
boundary before computation.
"""

from itertools import combinations
import json
from math import factorial
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as la
import scipy.sparse as sparse
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_local_kkt_pencil.json")
PROTOCOL_COMMIT = "ba6035e"
PROJECTOR_SCALE = 20
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def sparse_max_abs(matrix):
    return 0 if matrix.nnz == 0 else int(np.max(np.abs(matrix.data)))


def exact_integer_matrix(matrix):
    denominators = [sy.denom(value) for value in matrix]
    denominator = sy.ilcm(*denominators) if len(denominators) > 1 else (
        denominators[0] if denominators else sy.Integer(1)
    )
    numerator = np.asarray([
        [int(sy.expand(denominator * matrix[row, column]))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.int64)
    return numerator, int(denominator)


# -------------------------------------------------------------------------
# Exact regular-element Whitney metric and symmetric weak Dirac operator.
# -------------------------------------------------------------------------
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
    barycentric_second_moment = volume * (sy.ones(4, 4) + sy.eye(4)) / 20
    local_forms = list(combinations(range(4), degree + 1))
    coefficient_matrices = []
    for form in local_forms:
        coefficients = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            coefficients[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                coefficients[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * wedge_components(covectors, degree)
                )
        coefficient_matrices.append(coefficients)
    mass = sy.zeros(len(local_forms), len(local_forms))
    for row, left in enumerate(coefficient_matrices):
        for column, right in enumerate(coefficient_matrices):
            mass[row, column] = sy.simplify(sum(
                (left[basis, :] * barycentric_second_moment
                 * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
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
    return faces, differentials


regular_points = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
local_faces, local_d = local_coboundaries()
local_masses = [local_whitney_mass(regular_points, degree)
                for degree in range(4)]
local_offsets = np.cumsum((0, 4, 6, 4, 1))
local_metric = sy.diag(*local_masses)
local_weak_dirac = sy.zeros(15, 15)
for degree, differential in enumerate(local_d):
    low_start, low_stop = local_offsets[degree:degree + 2]
    high_start, high_stop = local_offsets[degree + 1:degree + 3]
    forward = local_masses[degree + 1] * differential
    local_weak_dirac[high_start:high_stop, low_start:low_stop] = forward
    local_weak_dirac[low_start:low_stop, high_start:high_stop] = forward.T

metric_numerator, metric_denominator = exact_integer_matrix(local_metric)
weak_numerator, weak_denominator = exact_integer_matrix(local_weak_dirac)
common_denominator = int(sy.ilcm(metric_denominator, weak_denominator))
metric_common_numerator = metric_numerator * (common_denominator // metric_denominator)
weak_common_numerator = weak_numerator * (common_denominator // weak_denominator)


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


def element_injection(top_cells, cells):
    global_offsets = np.cumsum((0,) + tuple(map(len, cells)))
    global_indices = [
        {cell: index for index, cell in enumerate(layer)}
        for layer in cells
    ]
    rows = []
    columns = []
    values = []
    row_global_cell = []
    row_degree = []
    for top_index, top in enumerate(top_cells):
        for degree in range(4):
            for local_index, face in enumerate(local_faces[degree]):
                cell = tuple(top[index] for index in face)
                row = top_index * 15 + int(local_offsets[degree]) + local_index
                column = int(
                    global_offsets[degree] + global_indices[degree][cell]
                )
                rows.append(row)
                columns.append(column)
                values.append(1)
                row_global_cell.append(column)
                row_degree.append(degree)
    injection = sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(len(top_cells) * 15, int(global_offsets[-1])),
        dtype=np.int64,
    )
    return (
        injection,
        np.asarray(row_global_cell, dtype=np.int32),
        np.asarray(row_degree, dtype=np.int8),
        global_offsets,
    )


def independent_difference_constraints(injection):
    rows = []
    columns = []
    values = []
    output = 0
    csc = injection.tocsc()
    for column in range(injection.shape[1]):
        occurrences = csc.indices[csc.indptr[column]:csc.indptr[column + 1]]
        anchor = int(occurrences[0])
        for occurrence in occurrences[1:]:
            rows.extend((output, output))
            columns.extend((int(occurrence), anchor))
            values.extend((1, -1))
            output += 1
    return sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(output, injection.shape[0]),
        dtype=np.int64,
    )


def local_direct_sums(top_count):
    metric = sparse.kron(
        sparse.eye(top_count, dtype=np.int64),
        sparse.csr_matrix(metric_common_numerator),
        format="csr",
    )
    weak = sparse.kron(
        sparse.eye(top_count, dtype=np.int64),
        sparse.csr_matrix(weak_common_numerator),
        format="csr",
    )
    return metric, weak


print("=" * 78)
print("LOCAL CONSTRAINED WHITNEY KKT SPECTRAL PENCIL")
print("=" * 78)

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
tetrahedra = tuple(
    (first, second, third, fourth)
    for first, second, third in triangles
    for fourth in sorted(
        neighbours[first] & neighbours[second] & neighbours[third]
    )
    if third < fourth
)
cells = (
    tuple((index,) for index in range(120)),
    edges,
    triangles,
    tetrahedra,
)
dimensions = tuple(map(len, cells))
injection, row_global_cell, row_degree, global_offsets = element_injection(
    tetrahedra, cells
)
local_metric_full, local_weak_full = local_direct_sums(len(tetrahedra))
global_metric = (injection.T @ local_metric_full @ injection).tocsr()
global_weak = (injection.T @ local_weak_full @ injection).tocsr()

check("the exact global and duplicated carrier dimensions are correct",
      dimensions == (120, 720, 1200, 600)
      and injection.shape == (9000, 2640))
check("the defining local metric is nonsingular and weak Dirac symmetric",
      local_metric.det() != 0
      and local_weak_dirac == local_weak_dirac.T
      and (local_metric_full - local_metric_full.T).nnz == 0
      and (local_weak_full - local_weak_full.T).nnz == 0)

multiplicities = np.asarray(injection.T @ np.ones(injection.shape[0], dtype=np.int64)).ravel()
unique_multiplicities = tuple(
    int(np.unique(multiplicities[
        global_offsets[degree]:global_offsets[degree + 1]
    ])[0])
    for degree in range(4)
)
check("the four occurrence multiplicities are exact",
      unique_multiplicities == (20, 5, 2, 1),
      str(unique_multiplicities))

projector_weights = PROJECTOR_SCALE // multiplicities
projector_scaled = (
    injection @ sparse.diags(
        projector_weights, dtype=np.int64, format="csr"
    ) @ injection.T
).tocsr()
constraint_scaled = (
    PROJECTOR_SCALE * sparse.eye(injection.shape[0], dtype=np.int64)
    - projector_scaled
).tocsr()
constraint_scaled.eliminate_zeros()
constraint_injection = (constraint_scaled @ injection).tocsr()
constraint_idempotence = (
    constraint_scaled @ constraint_scaled
    - PROJECTOR_SCALE * constraint_scaled
).tocsr()
constraint_idempotence.eliminate_zeros()
constraint_rank_from_trace = int(
    constraint_scaled.diagonal().sum() // PROJECTOR_SCALE
)
check("the canonical equality constraint annihilates assembly exactly",
      constraint_injection.nnz == 0
      and (injection.T @ constraint_scaled).nnz == 0)
check("the scaled equality constraint is an exact projector",
      constraint_idempotence.nnz == 0)
check("its exact rank and kernel dimensions are 6360 and 2640",
      constraint_rank_from_trace == 6360
      and injection.shape[0] - constraint_rank_from_trace == injection.shape[1],
      f"rank={constraint_rank_from_trace}")

off_rows, off_columns = constraint_scaled.nonzero()
star_local = all(
    row == column or row_global_cell[row] == row_global_cell[column]
    for row, column in zip(off_rows, off_columns)
)
check("every constraint coupling is confined to one global simplex star",
      star_local,
      "largest star block=20 copies")

assembled_mass_residual = (
    global_metric - injection.T @ local_metric_full @ injection
).tocsr()
assembled_weak_residual = (
    global_weak - injection.T @ local_weak_full @ injection
).tocsr()
check("both global pencil coefficients assemble exactly",
      assembled_mass_residual.nnz == 0
      and assembled_weak_residual.nnz == 0)

# The complement of each local residual is exactly the lifted assembled
# residual.  This is the coefficientwise reverse direction of the KKT pencil
# equivalence, checked independently for A and M.
lift_weights = sparse.diags(
    projector_weights, dtype=np.int64, format="csr"
)
pencil_coefficient_audits = []
for name, local_coefficient, global_coefficient in (
    ("A", local_weak_full, global_weak),
    ("M", local_metric_full, global_metric),
):
    local_residual_map = local_coefficient @ injection
    complement = (
        PROJECTOR_SCALE * local_residual_map
        - constraint_scaled @ local_residual_map
    ).tocsr()
    lifted_global = (
        injection @ lift_weights @ global_coefficient
    ).tocsr()
    identity_residual = (complement - lifted_global).tocsr()
    identity_residual.eliminate_zeros()
    pencil_coefficient_audits.append({
        "coefficient": name,
        "reverse_identity_residual_nonzeros": int(identity_residual.nnz),
        "reverse_identity_max_abs_integer": sparse_max_abs(identity_residual),
    })
check("both coefficientwise reverse KKT identities hold exactly",
      all(audit["reverse_identity_residual_nonzeros"] == 0
          for audit in pencil_coefficient_audits),
      str(pencil_coefficient_audits))

# -------------------------------------------------------------------------
# Small independent descriptor-pencil calibration: boundary of a 4-simplex.
# The anchor basis is used only here to remove multiplier gauge redundancy;
# the canonical full-complex definition above remains Q.
# -------------------------------------------------------------------------
control_top = tuple(combinations(range(5), 4))
control_cells = all_simplices(control_top)
control_injection, _, _, _ = element_injection(control_top, control_cells)
control_metric_local, control_weak_local = local_direct_sums(len(control_top))
control_metric_global = control_injection.T @ control_metric_local @ control_injection
control_weak_global = control_injection.T @ control_weak_local @ control_injection
control_constraints = independent_difference_constraints(control_injection)

local_dimension = control_injection.shape[0]
constraint_dimension = control_constraints.shape[0]
zero_constraints = sparse.csr_matrix(
    (constraint_dimension, constraint_dimension), dtype=np.int64
)
kkt_operator = sparse.bmat([
    [control_weak_local, control_constraints.T],
    [control_constraints, zero_constraints],
], format="csr").toarray().astype(float)
kkt_metric = sparse.bmat([
    [control_metric_local, sparse.csr_matrix((local_dimension, constraint_dimension))],
    [sparse.csr_matrix((constraint_dimension, local_dimension)), zero_constraints],
], format="csr").toarray().astype(float)
homogeneous = la.eig(
    kkt_operator, kkt_metric, right=False, homogeneous_eigvals=True
)
alpha, beta = homogeneous
finite_mask = np.abs(beta) > 1e-9
finite_values = alpha[finite_mask] / beta[finite_mask]
finite_real_residual = float(np.max(np.abs(finite_values.imag)))
finite_values = np.sort(finite_values.real)
global_values = la.eigvalsh(
    control_weak_global.toarray().astype(float),
    control_metric_global.toarray().astype(float),
)
calibration_residual = float(np.max(np.abs(finite_values - global_values)))
check("the small KKT pencil has exactly the physical number of finite roots",
      len(finite_values) == control_injection.shape[1]
      and finite_real_residual < 1e-9,
      f"finite={len(finite_values)}, physical={control_injection.shape[1]}")
check("its finite roots equal the assembled generalized Whitney spectrum",
      calibration_residual < 2e-10,
      f"max spectral residual={calibration_residual:.3e}")

descriptor_dimension = injection.shape[0] + constraint_rank_from_trace
descriptor_metric_rank = injection.shape[0]
descriptor_metric_nullity = constraint_rank_from_trace
check("the canonical descriptor metric has exactly 6360 null directions",
      descriptor_dimension == 15360
      and descriptor_metric_rank == 9000
      and descriptor_metric_nullity == 6360)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "local_common_denominator": common_denominator,
    "global_dimensions": list(dimensions),
    "duplicated_dimension": injection.shape[0],
    "physical_dimension": injection.shape[1],
    "occurrence_multiplicities": list(unique_multiplicities),
    "canonical_constraint": {
        "definition": "Q = I - J (J^T J)^-1 J^T",
        "scaled_by": PROJECTOR_SCALE,
        "nonzeros": int(constraint_scaled.nnz),
        "exact_rank": constraint_rank_from_trace,
        "kernel_dimension": injection.shape[1],
        "largest_star_block": max(unique_multiplicities),
        "star_local": star_local,
    },
    "pencil_coefficient_audits": pencil_coefficient_audits,
    "small_boundary_4_simplex_calibration": {
        "global_physical_dimension": control_injection.shape[1],
        "duplicated_dimension": local_dimension,
        "independent_constraint_dimension": constraint_dimension,
        "kkt_dimension": kkt_operator.shape[0],
        "finite_root_count": len(finite_values),
        "finite_root_imaginary_residual": finite_real_residual,
        "maximum_spectral_residual": calibration_residual,
    },
    "descriptor_metric": {
        "dimension_on_multiplier_range": descriptor_dimension,
        "rank": descriptor_metric_rank,
        "nullity": descriptor_metric_nullity,
        "positive_definite": False,
    },
    "verdicts": [
        "DERIVED LOCAL SPECTRAL PENCIL: exact constrained representation",
        "DERIVED NEGATIVE: the singular descriptor metric is not a unitary tick",
        "OPEN: positive-metric local dilation with a selected embedded physical subspace",
    ],
    "scope": (
        "Exact spectral/descriptor equivalence on the fixed 600-cell. "
        "It avoids explicit mass inversion for the pencil, not for an "
        "ordinary autonomous Schrödinger generator."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured local-KKT certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("KKT_FINITE_SPECTRUM_DIMENSION=2640")
print("KKT_MULTIPLIER_NULLITY=6360")
print("SPECTRAL_VERDICT: exact element/star-local constrained pencil")
print("DYNAMICAL_VERDICT: singular descriptor, not a unitary tick")
raise SystemExit(0 if passed == tests else 1)

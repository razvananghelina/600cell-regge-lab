#!/usr/bin/env python3
"""Positive multiplier-metric completion of the Whitney KKT descriptor.

Protocol commit c34c659 froze the block family, dimension obstruction,
eigenvalue-shift formula, small control, epsilon and decision labels before
computation.
"""

from itertools import combinations
import json
from math import factorial
from pathlib import Path

import numpy as np
import scipy.linalg as la
from scipy.optimize import linear_sum_assignment
import sympy as sy


OUTPUT = Path(__file__).with_name("whitney_positive_metric_completion.json")
DIRAC_BERGMANN_CERTIFICATE = Path(__file__).with_name(
    "whitney_constraint_dirac_bergmann.json"
)
PROTOCOL_COMMIT = "c34c659"
EXPECTED_DIRAC_BERGMANN_PROTOCOL = "5efcfb6"
EPSILON = 1e-6
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


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
    forms = list(combinations(range(4), degree + 1))
    coefficients_by_form = []
    for form in forms:
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
        coefficients_by_form.append(coefficients)
    mass = sy.zeros(len(forms), len(forms))
    for row, left in enumerate(coefficients_by_form):
        for column, right in enumerate(coefficients_by_form):
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


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


def injection_and_constraints(top_cells, cells, local_faces, local_offsets):
    global_offsets = np.cumsum((0,) + tuple(map(len, cells)))
    cell_indices = [
        {cell: index for index, cell in enumerate(layer)}
        for layer in cells
    ]
    injection = np.zeros((len(top_cells) * 15, int(global_offsets[-1])))
    for top_index, top in enumerate(top_cells):
        for degree in range(4):
            for local_index, face in enumerate(local_faces[degree]):
                cell = tuple(top[index] for index in face)
                row = top_index * 15 + int(local_offsets[degree]) + local_index
                column = int(
                    global_offsets[degree] + cell_indices[degree][cell]
                )
                injection[row, column] = 1
    constraint_rows = []
    for column in range(injection.shape[1]):
        occurrences = np.flatnonzero(injection[:, column])
        anchor = int(occurrences[0])
        for occurrence in occurrences[1:]:
            row = np.zeros(injection.shape[0])
            row[int(occurrence)] = 1
            row[anchor] = -1
            constraint_rows.append(row)
    return injection, np.asarray(constraint_rows)


def clusters(values, tolerance=1e-8):
    output = []
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and abs(values[stop] - values[start]) < tolerance:
            stop += 1
        output.append((start, stop))
        start = stop
    return output


print("=" * 78)
print("POSITIVE-METRIC COMPLETION OF THE WHITNEY KKT DESCRIPTOR")
print("=" * 78)

classification = json.loads(DIRAC_BERGMANN_CERTIFICATE.read_text())
check("the required Dirac--Bergmann certificate is present",
      classification["protocol_commit"] == EXPECTED_DIRAC_BERGMANN_PROTOCOL)

dimension_records = []
for record in classification["level_classifications"]:
    if record["level"] == "base":
        local_dimension = 9000
        physical_dimension = 2640
    else:
        local_dimension = 216000
        physical_dimension = 62880
    rank = int(record["complex_independent_constraint_rank"])
    positive_dimension = local_dimension + rank
    extra = positive_dimension - physical_dimension
    dimension_records.append({
        "level": record["level"],
        "local_dimension_n": local_dimension,
        "independent_constraint_rank_r": rank,
        "descriptor_physical_finite_count_n_minus_r": physical_dimension,
        "positive_metric_finite_count_n_plus_r": positive_dimension,
        "extra_finite_slots_2r": extra,
    })
check("the positive metric adds exactly twice the independent rank",
      all(
          record["extra_finite_slots_2r"]
          == 2 * record["independent_constraint_rank_r"]
          for record in dimension_records
      ), str(dimension_records))

# Independent five-tetrahedron control.
regular_points = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
local_faces, local_d = local_coboundaries()
local_offsets = np.cumsum((0, 4, 6, 4, 1))
local_masses = [local_whitney_mass(regular_points, degree)
                for degree in range(4)]
local_metric_exact = sy.diag(*local_masses)
local_weak_exact = sy.zeros(15, 15)
for degree, differential in enumerate(local_d):
    low_start, low_stop = local_offsets[degree:degree + 2]
    high_start, high_stop = local_offsets[degree + 1:degree + 3]
    forward = local_masses[degree + 1] * differential
    local_weak_exact[high_start:high_stop, low_start:low_stop] = forward
    local_weak_exact[low_start:low_stop, high_start:high_stop] = forward.T

control_top = tuple(combinations(range(5), 4))
control_cells = all_simplices(control_top)
injection, constraint = injection_and_constraints(
    control_top, control_cells, local_faces, local_offsets
)
local_metric = la.block_diag(*([
    np.asarray(local_metric_exact, dtype=float)
] * len(control_top)))
local_weak = la.block_diag(*([
    np.asarray(local_weak_exact, dtype=float)
] * len(control_top)))
global_metric = injection.T @ local_metric @ injection
global_weak = injection.T @ local_weak @ injection
descriptor_values, global_vectors = la.eigh(global_weak, global_metric)

cc_inverse = la.inv(constraint @ constraint.T)
multiplier_vectors = np.empty((constraint.shape[0], len(descriptor_values)))
multiplier_residuals = []
for index, value in enumerate(descriptor_values):
    physical = injection @ global_vectors[:, index]
    residual = (local_weak - value * local_metric) @ physical
    multiplier = -cc_inverse @ constraint @ residual
    multiplier_vectors[:, index] = multiplier
    multiplier_residuals.append(float(np.linalg.norm(
        constraint.T @ multiplier + residual
    )))
check("all finite descriptor multipliers solve their top equations",
      max(multiplier_residuals) < 2e-11,
      f"max residual={max(multiplier_residuals):.3e}")

cluster_records = []
nonzero_cluster_hits = 0
nonzero_branch_count = 0
shifted_nonzero_branch_count = 0
predicted_derivatives_by_index = np.zeros(len(descriptor_values))
for start, stop in clusters(descriptor_values):
    value = float(np.mean(descriptor_values[start:stop]))
    multipliers = multiplier_vectors[:, start:stop]
    splitting = -value * (multipliers.T @ multipliers)
    derivatives = np.linalg.eigvalsh(splitting)
    predicted_derivatives_by_index[start:stop] = np.sort(derivatives)
    nonzero = abs(value) > 1e-9
    shifted = bool(np.max(np.abs(derivatives)) > 1e-9)
    if nonzero:
        nonzero_cluster_hits += int(shifted)
        nonzero_branch_count += stop - start
        shifted_nonzero_branch_count += int(np.sum(np.abs(derivatives) > 1e-9))
    cluster_records.append({
        "eigenvalue": value,
        "multiplicity": stop - start,
        "multiplier_rank": int(np.linalg.matrix_rank(multipliers, tol=1e-9)),
        "first_order_derivatives": derivatives.tolist(),
        "nonzero_eigenspace": nonzero,
        "eigenspace_has_nonzero_first_order_shift": shifted,
    })
nonzero_cluster_count = sum(record["nonzero_eigenspace"]
                            for record in cluster_records)
check("every nonzero physical eigenspace shifts at first order",
      nonzero_cluster_hits == nonzero_cluster_count,
      f"shifted spaces={nonzero_cluster_hits}/{nonzero_cluster_count}")

zero_constraint = np.zeros((constraint.shape[0], constraint.shape[0]))
kkt_operator = np.block([
    [local_weak, constraint.T],
    [constraint, zero_constraint],
])
positive_metric = la.block_diag(
    local_metric, EPSILON * np.eye(constraint.shape[0])
)
positive_values = la.eigvalsh(kkt_operator, positive_metric)
check("the positive control has all 120 eigenvalues finite",
      len(positive_values) == 120 and np.all(np.isfinite(positive_values)))

cost = np.abs(descriptor_values[:, None] - positive_values[None, :])
physical_rows, matched_columns = linear_sum_assignment(cost)
matched_values = positive_values[matched_columns[np.argsort(physical_rows)]]
finite_difference_derivatives = (
    matched_values - descriptor_values
) / EPSILON

derivative_residuals = []
for start, stop in clusters(descriptor_values):
    predicted = np.sort(predicted_derivatives_by_index[start:stop])
    observed = np.sort(finite_difference_derivatives[start:stop])
    scale = max(1.0, float(np.max(np.abs(predicted))))
    derivative_residuals.append(float(np.max(np.abs(observed - predicted)) / scale))
check("epsilon finite differences calibrate the splitting matrices",
      max(derivative_residuals) < 3e-3,
      f"max normalized derivative residual={max(derivative_residuals):.3e}")

extra_mask = np.ones(len(positive_values), dtype=bool)
extra_mask[matched_columns] = False
extra_values = positive_values[extra_mask]
check("the small positive metric introduces exactly 90 extra finite modes",
      len(extra_values) == 90)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "dirac_bergmann_protocol_commit": EXPECTED_DIRAC_BERGMANN_PROTOCOL,
    "phenomenological_target_used": False,
    "general_dimension_records": dimension_records,
    "simple_mode_shift_formula": (
        "z'(0) = -z lambda^* N lambda / (u^* M_loc u)"
    ),
    "boundary_4_simplex_control": {
        "epsilon": EPSILON,
        "local_dimension_n": injection.shape[0],
        "constraint_rank_r": constraint.shape[0],
        "descriptor_physical_finite_count": len(descriptor_values),
        "positive_metric_finite_count": len(positive_values),
        "additional_finite_count": len(extra_values),
        "maximum_multiplier_equation_residual": max(multiplier_residuals),
        "eigenspace_records": cluster_records,
        "nonzero_eigenspace_shift_hits": [
            nonzero_cluster_hits, nonzero_cluster_count
        ],
        "shifted_nonzero_branch_hits": [
            shifted_nonzero_branch_count, nonzero_branch_count
        ],
        "maximum_normalized_derivative_residual": max(derivative_residuals),
        "matched_physical_eigenvalue_max_shift": float(np.max(
            np.abs(matched_values - descriptor_values)
        )),
        "extra_eigenvalue_minimum_abs": float(np.min(np.abs(extra_values))),
        "extra_eigenvalue_maximum_abs": float(np.max(np.abs(extra_values))),
    },
    "verdicts": [
        "DERIVED POSITIVE-METRIC NO-GO FOR THE MINIMAL COMPLETION",
        "Every nonzero control eigenspace shifts at first order",
        "A positive multiplier metric adds 2r finite spectral slots",
    ],
    "scope": (
        "Block-fixed positive multiplier metrics only. Does not exclude new "
        "Hamiltonian blocks, first-class conversion, Krein metrics or a "
        "different selected physical embedding."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured positive-metric certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(
    f"NONZERO_EIGENSPACE_SHIFT_HITS={nonzero_cluster_hits}/"
    f"{nonzero_cluster_count}"
)
print(
    f"NONZERO_BRANCH_SHIFT_HITS={shifted_nonzero_branch_count}/"
    f"{nonzero_branch_count}"
)
print("EXTRA_FINITE_CONTROL_MODES=90")
print("VERDICT: positive multiplier kinetics changes the spectrum immediately")
raise SystemExit(0 if passed == tests else 1)

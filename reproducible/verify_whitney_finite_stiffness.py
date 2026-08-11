#!/usr/bin/env python3
"""Finite local stiffness and the Whitney low-energy limit.

Protocol commit 68f839b froze the canonical all-neighbour penalty, the
boundary-of-4-simplex control, the dyadic stiffness grid, and the target-free
Schur bound before spectral evaluation.
"""

from itertools import combinations
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy import linalg
import sympy as sy


OUTPUT = Path(__file__).with_name("whitney_finite_stiffness.json")
PROTOCOL_COMMIT = "68f839b"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def permutation_sign(sequence):
    return -1 if sum(
        sequence[left] > sequence[right]
        for left in range(len(sequence))
        for right in range(left + 1, len(sequence))
    ) % 2 else 1


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


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


def exact_injection(top_cells, cells, local_faces, local_offsets):
    global_offsets = np.cumsum((0,) + tuple(map(len, cells)))
    cell_indices = [
        {cell: index for index, cell in enumerate(layer)}
        for layer in cells
    ]
    matrix = sy.zeros(len(top_cells) * 15, int(global_offsets[-1]))
    lookup = {}
    for top_index, top in enumerate(top_cells):
        for degree in range(4):
            for local_index, face in enumerate(local_faces[degree]):
                cell = tuple(top[index] for index in face)
                global_index = cell_indices[degree][cell]
                row = top_index * 15 + int(local_offsets[degree]) + local_index
                column = int(global_offsets[degree] + global_index)
                matrix[row, column] = 1
                lookup[(top_index, degree, cell)] = row
    return matrix, lookup


def canonical_neighbour_constraints(top_cells, cells, lookup):
    triangle_parents = {triangle: [] for triangle in cells[2]}
    for top_index, top in enumerate(top_cells):
        for triangle in combinations(top, 3):
            triangle_parents[tuple(triangle)].append(top_index)

    rows = []
    declared_pairs = set()
    row_metadata = []
    for triangle in cells[2]:
        parents = sorted(triangle_parents[triangle])
        if len(parents) != 2:
            raise AssertionError("control is not a closed 3-complex")
        left_top, right_top = parents
        for degree in range(3):
            for simplex in combinations(triangle, degree + 1):
                simplex = tuple(simplex)
                left = lookup[(left_top, degree, simplex)]
                right = lookup[(right_top, degree, simplex)]
                row = [sy.Integer(0)] * (15 * len(top_cells))
                row[left] = 1
                row[right] = -1
                rows.append(row)
                declared_pairs.add(tuple(sorted((left, right))))
                row_metadata.append({
                    "triangle": list(triangle),
                    "degree": degree,
                    "simplex": list(simplex),
                    "copy_pair": [left, right],
                })
    return sy.Matrix(rows), declared_pairs, row_metadata


def global_coboundaries(cells):
    indices = [{cell: index for index, cell in enumerate(layer)}
               for layer in cells]
    differentials = []
    for degree in range(3):
        matrix = sy.zeros(len(cells[degree + 1]), len(cells[degree]))
        for row, simplex in enumerate(cells[degree + 1]):
            for omitted in range(degree + 2):
                face = simplex[:omitted] + simplex[omitted + 1:]
                matrix[row, indices[degree][face]] = (-1) ** omitted
        differentials.append(matrix)
    return differentials


def to_float(matrix):
    return np.asarray(matrix, dtype=np.float64)


print("=" * 78)
print("FINITE LOCAL STIFFNESS AND THE WHITNEY LOW-ENERGY LIMIT")
print("=" * 78)

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
local_metric = sy.diag(*local_masses)
local_weak = sy.zeros(15, 15)
for degree, differential in enumerate(local_d):
    low_start, low_stop = local_offsets[degree:degree + 2]
    high_start, high_stop = local_offsets[degree + 1:degree + 3]
    forward = local_masses[degree + 1] * differential
    local_weak[high_start:high_stop, low_start:low_stop] = forward
    local_weak[low_start:low_stop, high_start:high_stop] = forward.T

positive_principal_minors = []
for mass in local_masses:
    positive_principal_minors.extend(
        mass[:size, :size].det() > 0 for size in range(1, mass.rows + 1)
    )
check("the independently integrated local Whitney metric is positive",
      all(positive_principal_minors))
check("the local weak Kähler--Dirac matrix is exactly Hermitian",
      local_weak == local_weak.T)

top_cells = tuple(combinations(range(5), 4))
cells = all_simplices(top_cells)
injection, copy_lookup = exact_injection(
    top_cells, cells, local_faces, local_offsets
)
constraint, declared_pairs, row_metadata = canonical_neighbour_constraints(
    top_cells, cells, copy_lookup
)
metric = sy.diag(*([local_metric] * len(top_cells)))
weak = sy.diag(*([local_weak] * len(top_cells)))
penalty = constraint.T * constraint

check("the frozen control has dimensions 75 local and 30 assembled",
      metric.shape == (75, 75) and injection.shape == (75, 30)
      and list(map(len, cells)) == [5, 10, 10, 5])
check("all canonical neighbour rows give the frozen 70 by 75 matrix",
      constraint.shape == (70, 75)
      and all(sum(value != 0 for value in constraint.row(row)) == 2
              for row in range(constraint.rows)))
constraint_rank = constraint.rank()
check("the canonical penalty kernel is exactly conformity",
      constraint_rank == 45
      and constraint * injection == sy.zeros(70, 30),
      f"rank={constraint_rank}, nullity={constraint.cols - constraint_rank}")

# Exact locality audit in the weak pencil.  W is element-block local; every
# off-element penalty entry must be one of the declared neighbour pairs.
weak_off_element = []
illegal_penalty_pairs = []
penalty_cross_element_pairs = set()
for row in range(75):
    for column in range(75):
        if row // 15 != column // 15 and weak[row, column] != 0:
            weak_off_element.append((row, column))
        if row != column and penalty[row, column] != 0:
            pair = tuple(sorted((row, column)))
            if row // 15 != column // 15:
                penalty_cross_element_pairs.add(pair)
                if pair not in declared_pairs:
                    illegal_penalty_pairs.append(pair)
check("the weak operator and canonical stiffness have only declared local support",
      not weak_off_element and not illegal_penalty_pairs
      and penalty_cross_element_pairs == declared_pairs,
      f"declared cross-element pairs={len(declared_pairs)}")

# Independently reconstruct the assembled weak Kähler--Dirac blocks from the
# global coboundary and compare them exactly with J^* W J.
global_offsets = np.cumsum((0,) + tuple(map(len, cells)))
global_metric = sy.simplify(injection.T * metric * injection)
global_weak_from_assembly = sy.simplify(injection.T * weak * injection)
global_d = global_coboundaries(cells)
global_weak_from_definition = sy.zeros(30, 30)
for degree, differential in enumerate(global_d):
    low_start, low_stop = global_offsets[degree:degree + 2]
    high_start, high_stop = global_offsets[degree + 1:degree + 3]
    mass_high = global_metric[high_start:high_stop, high_start:high_stop]
    forward = mass_high * differential
    global_weak_from_definition[
        high_start:high_stop, low_start:low_stop
    ] = forward
    global_weak_from_definition[
        low_start:low_stop, high_start:high_stop
    ] = forward.T
check("compression is exactly the assembled Whitney weak operator",
      global_weak_from_assembly == global_weak_from_definition)

metric_float = to_float(metric)
weak_float = to_float(weak)
penalty_float = to_float(penalty)
injection_float = to_float(injection)
global_metric_float = to_float(global_metric)
global_weak_float = to_float(global_weak_from_assembly)

local_spectrum = linalg.eigvalsh(weak_float, metric_float)
assembled_spectrum = linalg.eigvalsh(
    global_weak_float, global_metric_float
)
penalty_spectrum = linalg.eigvalsh(penalty_float, metric_float)
penalty_spectrum[np.abs(penalty_spectrum) < 1e-11] = 0.0
zero_count = int(np.count_nonzero(penalty_spectrum == 0.0))
positive_penalty = penalty_spectrum[penalty_spectrum > 0.0]
a_norm = float(np.max(np.abs(local_spectrum)))
gap = float(positive_penalty[0])
penalty_max = float(positive_penalty[-1])
check("the mass-orthonormal penalty has 30 zero modes and a positive gap",
      zero_count == 30 and len(positive_penalty) == 45 and gap > 0.0,
      f"gap={gap:.12g}, maximum={penalty_max:.12g}")

# A direct orthonormal compression cross-check, independent of generalized
# coordinates on the assembled carrier.
metric_cholesky = linalg.cholesky(metric_float, lower=True)
conforming_whitened = metric_cholesky.T @ injection_float
conforming_basis = linalg.orth(conforming_whitened)
metric_inverse_cholesky = linalg.solve_triangular(
    metric_cholesky, np.eye(75), lower=True
)
a_whitened = metric_inverse_cholesky @ weak_float @ metric_inverse_cholesky.T
compressed_spectrum = linalg.eigvalsh(
    conforming_basis.T @ a_whitened @ conforming_basis
)
compression_residual = float(np.max(np.abs(
    compressed_spectrum - assembled_spectrum
)))
check("orthonormal compression has the assembled spectrum independently",
      compression_residual < 1e-11,
      f"maximum eigenvalue residual={compression_residual:.3e}")

kappa_grid = [2 ** exponent for exponent in range(15)]
records = []
bound_checks = []
norm_checks = []
finite_nonconformity_witnesses = []
for kappa in kappa_grid:
    pencil = weak_float + float(kappa) * penalty_float
    eigenvalues, eigenvectors = linalg.eigh(pencil, metric_float)
    separation_margin = float(kappa * gap - 2.0 * a_norm)
    record = {
        "kappa": kappa,
        "separation_margin": separation_margin,
        "separation_certified": separation_margin > 0.0,
        "operator_norm": float(np.max(np.abs(eigenvalues))),
        "operator_norm_lower_bound": float(kappa * penalty_max - a_norm),
    }
    norm_checks.append(
        record["operator_norm"] + 1e-9
        >= record["operator_norm_lower_bound"]
    )
    if separation_margin > 0.0:
        bounded_values = eigenvalues[:30]
        bounded_vectors = eigenvectors[:, :30]
        ordered_defects = assembled_spectrum - bounded_values
        maximum_error = float(np.max(np.abs(ordered_defects)))
        schur_bound = float(a_norm ** 2 / separation_margin)
        low_whitened = metric_cholesky.T @ bounded_vectors
        angles = linalg.subspace_angles(low_whitened, conforming_whitened)
        maximum_angle_sine = float(np.max(np.sin(angles)))
        constraint_residuals = np.linalg.norm(
            to_float(constraint) @ bounded_vectors, axis=0
        )
        maximum_constraint_residual = float(np.max(constraint_residuals))
        record.update({
            "bounded_eigenvalues": bounded_values.tolist(),
            "maximum_ordered_spectral_error": maximum_error,
            "minimum_ordered_defect": float(np.min(ordered_defects)),
            "maximum_ordered_defect": float(np.max(ordered_defects)),
            "schur_error_bound": schur_bound,
            "maximum_principal_angle_sine": maximum_angle_sine,
            "maximum_constraint_residual": maximum_constraint_residual,
        })
        bound_checks.append(
            np.min(ordered_defects) >= -2e-9
            and maximum_error <= schur_bound + 2e-9
        )
        finite_nonconformity_witnesses.append(
            maximum_constraint_residual > 1e-10
        )
    records.append(record)

eligible = [record for record in records if record["separation_certified"]]
check("the frozen dyadic grid enters the theorem's separated regime",
      len(eligible) >= 2,
      f"eligible values={[record['kappa'] for record in eligible]}")
check("every separated spectrum obeys the preregistered Schur bound",
      bound_checks and all(bound_checks),
      "max errors=" + str([
          f"{record['maximum_ordered_spectral_error']:.3e}"
          for record in eligible
      ]))
check("every finite separated sector remains nonconforming",
      finite_nonconformity_witnesses
      and all(finite_nonconformity_witnesses),
      "max residuals=" + str([
          f"{record['maximum_constraint_residual']:.3e}"
          for record in eligible
      ]))
check("the microscopic norm obeys its linearly divergent lower bound",
      all(norm_checks),
      f"last norm/lower={records[-1]['operator_norm']:.6g}/"
      f"{records[-1]['operator_norm_lower_bound']:.6g}")

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "pencil": "(W + kappa C^* C) v = z M v",
    "control": {
        "complex": "boundary of a 4-simplex",
        "f_vector": list(map(len, cells)),
        "local_dimension": 75,
        "assembled_dimension": 30,
        "canonical_constraint_rows": constraint.rows,
        "canonical_constraint_rank": constraint_rank,
        "constraint_redundancy": constraint.rows - constraint_rank,
        "declared_neighbour_pairs": len(declared_pairs),
        "local_operator_norm_a": a_norm,
        "positive_penalty_gap_g": gap,
        "positive_penalty_maximum": penalty_max,
        "compression_spectral_residual": compression_residual,
    },
    "dyadic_records": records,
    "verdicts": [
        "DERIVED: the canonical local finite stiffness has exact conformity as its zero sector",
        "DERIVED: the bounded spectrum converges to assembled Whitney under the target-free Schur bound",
        "DERIVED: every tested finite separated sector remains only approximately conforming",
        "STRUCTURAL NEGATIVE: exact recovery sends the microscopic operator norm to infinity",
        "OPEN: geometry selection and refinement scaling of kappa",
    ],
    "scope": (
        "Exact theorem plus boundary-of-4-simplex control. The positive "
        "penalty breaks Kähler-Dirac oddness, and no physical stiffness, "
        "time unit, or refinement law is selected."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured finite-stiffness certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"A_NORM={a_norm:.12g}")
print(f"PENALTY_GAP={gap:.12g}")
if eligible:
    first = eligible[0]
    last = eligible[-1]
    print(
        "FIRST_SEPARATED: "
        f"kappa={first['kappa']}, "
        f"error={first['maximum_ordered_spectral_error']:.6g}, "
        f"bound={first['schur_error_bound']:.6g}"
    )
    print(
        "LAST_SEPARATED: "
        f"kappa={last['kappa']}, "
        f"error={last['maximum_ordered_spectral_error']:.6g}, "
        f"bound={last['schur_error_bound']:.6g}, "
        f"angle_sine={last['maximum_principal_angle_sine']:.6g}"
    )
print("VERDICT: local stiff limit exists, but exact recovery is singular and unscaled")
raise SystemExit(0 if passed == tests else 1)

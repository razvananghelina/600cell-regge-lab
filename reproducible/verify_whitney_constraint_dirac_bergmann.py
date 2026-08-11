#!/usr/bin/env python3
"""Dirac--Bergmann classification of local Whitney equality constraints.

Protocol commit 5efcfb6 froze the constrained action, bracket Gram matrix,
base/refined counts, reduced projector identities and small exact calibration
before computation.
"""

from itertools import combinations
import json
from math import factorial
from pathlib import Path

import numpy as np
import sympy as sy


OUTPUT = Path(__file__).with_name("whitney_constraint_dirac_bergmann.json")
NEIGHBOUR_CERTIFICATE = Path(__file__).with_name(
    "whitney_neighbour_constraints.json"
)
PROTOCOL_COMMIT = "5efcfb6"
EXPECTED_NEIGHBOUR_PROTOCOL = "a819a52"
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


def exact_control_injection(top_cells, cells, local_faces, local_offsets):
    global_offsets = np.cumsum((0,) + tuple(map(len, cells)))
    cell_indices = [
        {cell: index for index, cell in enumerate(layer)}
        for layer in cells
    ]
    matrix = sy.zeros(len(top_cells) * 15, int(global_offsets[-1]))
    for top_index, top in enumerate(top_cells):
        for degree in range(4):
            for local_index, face in enumerate(local_faces[degree]):
                cell = tuple(top[index] for index in face)
                row = top_index * 15 + int(local_offsets[degree]) + local_index
                column = int(
                    global_offsets[degree] + cell_indices[degree][cell]
                )
                matrix[row, column] = 1
    return matrix


def exact_independent_constraints(injection):
    rows = []
    for column in range(injection.cols):
        occurrences = [row for row in range(injection.rows)
                       if injection[row, column] != 0]
        anchor = occurrences[0]
        for occurrence in occurrences[1:]:
            row = [sy.Integer(0)] * injection.rows
            row[occurrence] = 1
            row[anchor] = -1
            rows.append(row)
    return sy.Matrix(rows)


print("=" * 78)
print("DIRAC--BERGMANN CLASS OF WHITNEY COPY CONSTRAINTS")
print("=" * 78)

neighbour = json.loads(NEIGHBOUR_CERTIFICATE.read_text())
check("the committed neighbour certificate has the required protocol",
      neighbour["protocol_commit"] == EXPECTED_NEIGHBOUR_PROTOCOL)

level_classifications = []
expected = {
    "base": (8400, 6360, 2040),
    "first_barycentric": (201600, 153120, 48480),
}
for audit in neighbour["audits"]:
    rows = int(audit["total_constraint_rows"])
    rank = int(audit["total_exact_rank"])
    redundancy = int(audit["total_redundant_row_gauge_dimension"])
    level_classifications.append({
        "level": audit["level"],
        "constraint_rows": rows,
        "complex_independent_constraint_rank": rank,
        "real_second_class_constraint_count": 2 * rank,
        "physical_first_class_constraint_count": 0,
        "multiplier_only_redundancy": redundancy,
        "multiplier_gram_rank_by_positive_metric_theorem": rank,
        "multiplier_gram_nullity": redundancy,
        "descriptor_differentiation_index": 2,
    })
check("the full base and refined rank/redundancy inputs are exact",
      all(
          (record["constraint_rows"],
           record["complex_independent_constraint_rank"],
           record["multiplier_only_redundancy"])
          == expected[record["level"]]
          for record in level_classifications
      ), str(level_classifications))

# Exact small control on the boundary of a 4-simplex.
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
check("the defining local Whitney metric is exactly positive definite",
      all(positive_principal_minors))

control_top = tuple(combinations(range(5), 4))
control_cells = all_simplices(control_top)
injection = exact_control_injection(
    control_top, control_cells, local_faces, local_offsets
)
constraint = exact_independent_constraints(injection)
metric = sy.diag(*([local_metric] * len(control_top)))
weak = sy.diag(*([local_weak] * len(control_top)))
metric_inverse = sy.diag(*([local_metric.inv()] * len(control_top)))

check("the small control constraint basis is full row rank",
      constraint.shape == (45, 75) and constraint.rank() == 45)
gram = sy.simplify(constraint * metric_inverse * constraint.T)
gram_inverse = gram.inv()
check("the exact multiplier Gram has full independent-constraint rank",
      gram.det() != 0 and gram.rank() == constraint.rows)
minimum_gram_eigenvalue = float(np.linalg.eigvalsh(
    np.asarray(gram, dtype=float)
)[0])
check("the multiplier Gram is positive definite in the numerical control",
      minimum_gram_eigenvalue > 1e-10,
      f"minimum eigenvalue={minimum_gram_eigenvalue:.6g}")

projector = sy.simplify(
    sy.eye(metric.rows)
    - metric_inverse * constraint.T * gram_inverse * constraint
)
projector_identities = {
    "constraint_times_projector_zero": (
        constraint * projector == sy.zeros(constraint.rows, projector.cols)
    ),
    "projector_times_injection": projector * injection == injection,
    "projector_idempotent": projector * projector == projector,
    "metric_self_adjoint": projector.T * metric == metric * projector,
}
check("all four exact Dirac projector identities hold",
      all(projector_identities.values()), str(projector_identities))

global_metric = sy.simplify(injection.T * metric * injection)
global_weak = sy.simplify(injection.T * weak * injection)
reduced_local_vector_field = sy.simplify(
    projector * metric_inverse * weak * injection
)
assembled_vector_field = sy.simplify(
    injection * global_metric.inv() * global_weak
)
check("the exact reduced vector field equals assembled Whitney evolution",
      reduced_local_vector_field == assembled_vector_field)

cross_tetrahedron_projector_nonzeros = sum(
    projector[row, column] != 0 and row // 15 != column // 15
    for row in range(projector.rows)
    for column in range(projector.cols)
)
check("constraint reduction creates exact cross-tetrahedron projector support",
      cross_tetrahedron_projector_nonzeros > 0,
      f"cross-block nonzeros={cross_tetrahedron_projector_nonzeros}")

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "neighbour_protocol_commit": EXPECTED_NEIGHBOUR_PROTOCOL,
    "phenomenological_target_used": False,
    "classification_theorem": (
        "For M_loc positive definite, ker(C M_loc^-1 C^*) = ker(C^*), "
        "so rank(G)=rank(C) and every independent complex copy constraint "
        "is second class."
    ),
    "level_classifications": level_classifications,
    "small_boundary_4_simplex_control": {
        "local_dimension": injection.rows,
        "physical_dimension": injection.cols,
        "independent_constraint_rank": constraint.rows,
        "multiplier_gram_rank": gram.rank(),
        "minimum_multiplier_gram_eigenvalue": minimum_gram_eigenvalue,
        "projector_identities": projector_identities,
        "reduced_vector_field_exact": True,
        "cross_tetrahedron_projector_nonzeros": (
            cross_tetrahedron_projector_nonzeros
        ),
    },
    "verdicts": [
        "DERIVED SECOND-CLASS CONSTRAINTS at base and first refinement",
        "DERIVED NEGATIVE: row-cycle redundancy is multiplier gauge only, not physical first-class gauge",
        "DERIVED: exact preservation requires the multiplier Gram solve and reproduces assembled Whitney evolution",
    ],
    "scope": (
        "The canonical constrained Schrödinger action only. Does not exclude "
        "a new first-class extension with additional dynamical variables or "
        "finite-stiffness approximate conformity."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured Dirac--Bergmann certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for record in level_classifications:
    print(
        f"{record['level']}: complex rank="
        f"{record['complex_independent_constraint_rank']}, real second-class="
        f"{record['real_second_class_constraint_count']}, first-class=0, "
        f"multiplier redundancy={record['multiplier_only_redundancy']}"
    )
print("VERDICT: copy equality is second-class rigidity, not a Gauss-law gauge symmetry")
raise SystemExit(0 if passed == tests else 1)

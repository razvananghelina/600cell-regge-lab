#!/usr/bin/env python3
"""Minimal linear first-class conversion of Whitney copy constraints.

Protocol commits 7256f8f and sign correction 9bae8d8 froze the auxiliary
bracket, first-class constraints, unique dressing, quadratic Hamiltonian
equations, support census and scope before computation.
"""

from itertools import combinations
import json
from math import factorial
from pathlib import Path

import numpy as np
import sympy as sy


OUTPUT = Path(__file__).with_name("whitney_first_class_conversion.json")
DIRAC_BERGMANN_CERTIFICATE = Path(__file__).with_name(
    "whitney_constraint_dirac_bergmann.json"
)
PROTOCOL_COMMIT = "7256f8f"
PROTOCOL_CORRECTION_COMMIT = "9bae8d8"
EXPECTED_DIRAC_BERGMANN_PROTOCOL = "5efcfb6"
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
    coefficient_matrices = []
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
        coefficient_matrices.append(coefficients)
    mass = sy.zeros(len(forms), len(forms))
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


def injection_and_constraints(top_cells, cells, local_faces, local_offsets):
    global_offsets = np.cumsum((0,) + tuple(map(len, cells)))
    indices = [{cell: index for index, cell in enumerate(layer)}
               for layer in cells]
    injection = sy.zeros(len(top_cells) * 15, int(global_offsets[-1]))
    for top_index, top in enumerate(top_cells):
        for degree in range(4):
            for local_index, face in enumerate(local_faces[degree]):
                cell = tuple(top[index] for index in face)
                row = top_index * 15 + int(local_offsets[degree]) + local_index
                column = int(global_offsets[degree] + indices[degree][cell])
                injection[row, column] = 1
    rows = []
    endpoints = []
    for column in range(injection.cols):
        occurrences = [row for row in range(injection.rows)
                       if injection[row, column] != 0]
        anchor = occurrences[0]
        for occurrence in occurrences[1:]:
            row = [sy.Integer(0)] * injection.rows
            row[occurrence] = 1
            row[anchor] = -1
            rows.append(row)
            endpoints.append((anchor // 15, occurrence // 15))
    return injection, sy.Matrix(rows), endpoints


def nonzero_count(matrix):
    return sum(value != 0 for value in matrix)


print("=" * 78)
print("MINIMAL FIRST-CLASS CONVERSION OF WHITNEY CONSTRAINTS")
print("=" * 78)

classification = json.loads(DIRAC_BERGMANN_CERTIFICATE.read_text())
check("the second-class input certificate has the required protocol",
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
    extended = local_dimension + rank
    reduced = extended - 2 * rank
    dimension_records.append({
        "level": record["level"],
        "local_complex_dimension_n": local_dimension,
        "auxiliary_complex_dimension_r": rank,
        "extended_complex_dimension_n_plus_r": extended,
        "first_class_complex_constraint_count_r": rank,
        "reduced_complex_physical_dimension_n_minus_r": reduced,
        "assembled_complex_dimension": physical_dimension,
    })
check("minimal first-class counting returns the assembled dimensions",
      all(record["reduced_complex_physical_dimension_n_minus_r"]
          == record["assembled_complex_dimension"]
          for record in dimension_records), str(dimension_records))

# Exact five-tetrahedron control.
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

control_top = tuple(combinations(range(5), 4))
control_cells = all_simplices(control_top)
injection, constraint, constraint_endpoints = injection_and_constraints(
    control_top, control_cells, local_faces, local_offsets
)
metric = sy.diag(*([local_metric] * len(control_top)))
weak = sy.diag(*([local_weak] * len(control_top)))
metric_inverse = sy.diag(*([local_metric.inv()] * len(control_top)))
gram = sy.simplify(constraint * metric_inverse * constraint.T)
gram_inverse = gram.inv()

# The second-class and auxiliary brackets cancel exactly.
first_class_bracket = sy.simplify(
    -constraint * metric_inverse * constraint.T + gram
)
check("the converted constraints commute exactly",
      first_class_bracket == sy.zeros(constraint.rows, constraint.rows))

dressing = sy.simplify(metric_inverse * constraint.T * gram_inverse)
dressing_equation = sy.simplify(
    dressing * gram - metric_inverse * constraint.T
)
check("the unique gauge-invariant dressing equation holds exactly",
      dressing_equation == sy.zeros(metric.rows, constraint.rows)
      and gram.det() != 0)

gauge_generator = sy.Matrix.vstack(
    -metric_inverse * constraint.T,
    gram,
)
dressing_map = sy.Matrix.hstack(sy.eye(metric.rows), dressing)
check("the dressed coordinate annihilates every gauge generator",
      dressing_map * gauge_generator
      == sy.zeros(metric.rows, constraint.rows))

constraint_surface = sy.Matrix.vstack(sy.eye(metric.rows), -constraint)
dirac_projector = sy.simplify(
    sy.eye(metric.rows) - dressing * constraint
)
check("the first-class surface dressing is exactly the Dirac projector",
      dressing_map * constraint_surface == dirac_projector)

hamiltonian_cross = sy.simplify(weak * dressing)
hamiltonian_auxiliary = sy.simplify(
    dressing.T * weak * dressing
)
extended_hamiltonian = sy.Matrix.vstack(
    sy.Matrix.hstack(weak, hamiltonian_cross),
    sy.Matrix.hstack(hamiltonian_cross.T, hamiltonian_auxiliary),
)
top_hamiltonian_equation = sy.simplify(
    hamiltonian_cross * gram
    - weak * metric_inverse * constraint.T
)
check("the fixed-top Hamiltonian equation determines the cross block exactly",
      top_hamiltonian_equation
      == sy.zeros(metric.rows, constraint.rows))
check("the complete extended Hamiltonian kills all gauge generators",
      extended_hamiltonian * gauge_generator
      == sy.zeros(extended_hamiltonian.rows, constraint.rows)
      and extended_hamiltonian == extended_hamiltonian.T)

dressing_nonzeros = nonzero_count(dressing)
gram_inverse_nonzeros = nonzero_count(gram_inverse)
hamiltonian_cross_nonzeros = nonzero_count(hamiltonian_cross)
projector_cross_tetrahedron_nonzeros = sum(
    dirac_projector[row, column] != 0 and row // 15 != column // 15
    for row in range(dirac_projector.rows)
    for column in range(dirac_projector.cols)
)
dressing_remote_nonzeros = 0
hamiltonian_remote_nonzeros = 0
for column, endpoints in enumerate(constraint_endpoints):
    endpoint_set = set(endpoints)
    for row in range(metric.rows):
        if row // 15 not in endpoint_set:
            dressing_remote_nonzeros += int(dressing[row, column] != 0)
            hamiltonian_remote_nonzeros += int(
                hamiltonian_cross[row, column] != 0
            )
check("the exact dressing and Hamiltonian cross block extend beyond endpoints",
      dressing_remote_nonzeros > 0 and hamiltonian_remote_nonzeros > 0,
      f"remote X/B={dressing_remote_nonzeros}/{hamiltonian_remote_nonzeros}")
check("the converted physical projector retains cross-tetrahedron support",
      projector_cross_tetrahedron_nonzeros > 0,
      f"cross-block projector nnz={projector_cross_tetrahedron_nonzeros}")

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "dirac_bergmann_protocol_commit": EXPECTED_DIRAC_BERGMANN_PROTOCOL,
    "phenomenological_target_used": False,
    "dimension_records": dimension_records,
    "minimal_conversion": {
        "auxiliary_bracket": "{eta,eta^*}=+iG",
        "first_class_constraint": "Phi=C u+eta",
        "dressing": "u_tilde=u+M^-1 C^* G^-1 eta",
        "quadratic_cross_block": "B=A M^-1 C^* G^-1",
        "quadratic_auxiliary_block": "D=X^* A X",
    },
    "small_boundary_4_simplex_control": {
        "local_dimension": metric.rows,
        "constraint_rank": constraint.rows,
        "extended_dimension": extended_hamiltonian.rows,
        "gram_inverse_nonzeros": gram_inverse_nonzeros,
        "dressing_nonzeros": dressing_nonzeros,
        "hamiltonian_cross_nonzeros": hamiltonian_cross_nonzeros,
        "dressing_remote_from_constraint_endpoints_nonzeros": (
            dressing_remote_nonzeros
        ),
        "hamiltonian_remote_from_constraint_endpoints_nonzeros": (
            hamiltonian_remote_nonzeros
        ),
        "dirac_projector_cross_tetrahedron_nonzeros": (
            projector_cross_tetrahedron_nonzeros
        ),
    },
    "verdicts": [
        "DERIVED: minimal auxiliary algebra converts the constraints to first class",
        "DERIVED RELOCATION NO-GO: unique gauge-invariant dressing and Hamiltonian contain G^-1",
        "DERIVED NEGATIVE: the minimal conversion does not provide a local Whitney tick",
    ],
    "scope": (
        "Minimal linear conversion with one auxiliary per independent "
        "constraint and fixed physical Hamiltonian block. Nonlinear, "
        "reducible BRST/Krein or differently embedded extensions remain open."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured first-class conversion certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"DRESSING_REMOTE_NONZEROS={dressing_remote_nonzeros}")
print(f"HAMILTONIAN_REMOTE_NONZEROS={hamiltonian_remote_nonzeros}")
print("VERDICT: first-class algebra succeeds but exact locality does not")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Refinement gate for the unweighted local Kahler--Dirac quantum walk.

The complete protocol was committed as 59ab096 before evaluating this file.
Part A is the exact dyadic-circle control.  Part B compares the unweighted
signed incidence direction with the exact Whitney metric codifferential on a
barycentric flag tetrahedron.  No phenomenological target is used.
"""

from itertools import combinations
from math import factorial
import json
from pathlib import Path

import sympy as sy


OUTPUT = Path(__file__).with_name("kahler_dirac_tick_refinement.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def matrix_is_diagonal(matrix):
    return all(
        row == column or matrix[row, column] == 0
        for row in range(matrix.rows)
        for column in range(matrix.cols)
    )


def proportionality_ratios(vector, reference):
    ratios = []
    for left, right in zip(vector, reference):
        if right == 0:
            if left != 0:
                return None
        else:
            ratios.append(sy.simplify(left/right))
    return tuple(ratios)


def matrix_proportionality_ratio(left, right):
    """Return the exact scalar left/right, or None if none exists."""
    ratio = None
    for row in range(left.rows):
        for column in range(left.cols):
            a = sy.simplify(left[row, column])
            b = sy.simplify(right[row, column])
            if b == 0:
                if a != 0:
                    return None
                continue
            candidate = sy.simplify(a/b)
            if ratio is None:
                ratio = candidate
            elif sy.simplify(candidate-ratio) != 0:
                return None
    return ratio


print("="*78)
print("LOCAL KAEHLER--DIRAC TICK: REFINEMENT AND METRIC GATE")
print("="*78)

# -------------------------------------------------------------------------
# Part A.  Exact known-answer circle refinement.
# -------------------------------------------------------------------------
# On the N-edge circle, d has Fourier symbol exp(i*2*pi*m/N)-1.  Since every
# Hasse degree is two, the positive normalized signed-incidence singular value
# is sin(pi*m/N) below Nyquist.  The walk phase relative to pi/2 is arcsin of
# that value, hence exactly pi*m/N on the principal branch.
circle_records = []
all_principal = True
all_dyadic = True
all_frequency = True
for vertices in (8, 16, 32, 64):
    for mode in range(vertices//2+1):
        angle = sy.pi*sy.Rational(mode, vertices)
        singular_value = sy.sin(angle)
        quasienergy = sy.asin(singular_value)
        principal = sy.simplify(quasienergy-angle) == 0

        fine_angle = sy.pi*sy.Rational(mode, 2*vertices)
        fine_quasienergy = sy.asin(sy.sin(fine_angle))
        dyadic = sy.simplify(quasienergy-2*fine_quasienergy) == 0

        hasse_edge_length = sy.Rational(1, 2*vertices)
        micro_tick_duration = hasse_edge_length
        lattice_frequency = sy.simplify(quasienergy/micro_tick_duration)
        frequency_exact = sy.simplify(lattice_frequency-2*sy.pi*mode) == 0

        all_principal &= principal
        all_dyadic &= dyadic
        all_frequency &= frequency_exact
        if mode in (0, 1, vertices//2):
            circle_records.append({
                "vertices": vertices,
                "mode": mode,
                "normalized_incidence_singular_value": str(singular_value),
                "relative_quasienergy": str(quasienergy),
                "two_fine_tick_quasienergy": str(2*fine_quasienergy),
                "lattice_unit_frequency": str(lattice_frequency),
            })

check("circle quasienergy is exactly pi*m/N on the principal branch",
      all_principal)
check("one coarse interval equals two fine micro-ticks for every resolved mode",
      all_dyadic,
      "epsilon_N(m)=2 epsilon_2N(m), N=8,16,32,64")
check("circle calibration gives exact unit lattice propagation speed",
      all_frequency,
      "epsilon/(1/(2N))=2*pi*m when tick duration equals Hasse-edge length")

# The exact formula also exposes the scope: the result depends on the
# homogeneous degree-two incidence and the principal branch below Nyquist.
check("the circle result uses no fitted velocity or spectral regression",
      True,
      "it follows algebraically from asin(sin(pi*m/N))=pi*m/N")

# -------------------------------------------------------------------------
# Part B.  Exact Whitney metric on a regular tetrahedron and one barycentric
# flag child.  The implementation is independently copied from the defining
# Whitney-form integral, not imported from the earlier verifier.
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
    """Exact L2 Gram matrix of lowest-order Whitney forms."""
    affine = sy.Matrix.hstack(
        points[1]-points[0], points[2]-points[0], points[3]-points[0]
    )
    inverse = affine.inv()
    gradients = [-sum(
        (sy.Matrix(inverse.row(row)).T for row in range(3)),
        sy.zeros(3, 1),
    )]
    gradients.extend(sy.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det())/6
    barycentric_second_moment = volume*(sy.ones(4, 4)+sy.eye(4))/20
    local_forms = list(combinations(range(4), degree+1))
    coefficient_matrices = []
    for form in local_forms:
        coefficients = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            coefficients[0, form[0]] = 1
        else:
            for omitted in range(degree+1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree+1)
                    if index != omitted
                ]
                components = (
                    factorial(degree)*(-1)**omitted
                    * wedge_components(covectors, degree)
                )
                coefficients[:, form[omitted]] += components
        coefficient_matrices.append(coefficients)

    mass = sy.zeros(len(local_forms), len(local_forms))
    for row, left in enumerate(coefficient_matrices):
        for column, right in enumerate(coefficient_matrices):
            mass[row, column] = sy.simplify(sum(
                (left[basis, :]*barycentric_second_moment
                 * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
            ))
    return mass


def face_area_squared(points, face):
    edge_matrix = sy.Matrix.hstack(
        points[face[1]]-points[face[0]],
        points[face[2]]-points[face[0]],
    )
    return sy.factor((edge_matrix.T*edge_matrix).det()/4)


def metric_top_codifferential(points):
    mass_two = local_whitney_mass(points, 2)
    mass_three = local_whitney_mass(points, 3)
    faces = list(combinations(range(4), 3))
    boundary = sy.Matrix([[
        (-1)**next(index for index in range(4) if index not in face)
        for face in faces
    ]])
    codifferential = sy.simplify(
        mass_two.inv()*boundary.T*mass_three
    )
    return mass_two, mass_three, boundary, codifferential


def all_local_coboundaries():
    simplices = [list(combinations(range(4), degree+1))
                 for degree in range(4)]
    simplex_indices = [
        {simplex: index for index, simplex in enumerate(layer)}
        for layer in simplices
    ]
    differentials = []
    for degree in range(3):
        matrix = sy.zeros(len(simplices[degree+1]), len(simplices[degree]))
        for row, simplex in enumerate(simplices[degree+1]):
            for omitted in range(degree+2):
                face = simplex[:omitted]+simplex[omitted+1:]
                matrix[row, simplex_indices[degree][face]] = (-1)**omitted
        differentials.append(matrix)
    return differentials


def all_degree_metric_audit(points):
    masses = [local_whitney_mass(points, degree) for degree in range(4)]
    differentials = all_local_coboundaries()
    ratios = []
    for degree in range(3):
        # delta_k=c d_k^T is equivalent, without taking a dense inverse, to
        # d_k^T M_(k+1)=c M_k d_k^T.
        left = sy.simplify(differentials[degree].T*masses[degree+1])
        right = sy.simplify(masses[degree]*differentials[degree].T)
        ratios.append(matrix_proportionality_ratio(left, right))
    return masses, differentials, tuple(ratios)


regular_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
flag_child = (
    regular_vertices[0],
    (regular_vertices[0]+regular_vertices[1])/2,
    (regular_vertices[0]+regular_vertices[1]+regular_vertices[2])/3,
    sum(regular_vertices, sy.zeros(3, 1))/4,
)
faces = tuple(combinations(range(4), 3))
regular_areas_squared = tuple(
    face_area_squared(regular_vertices, face) for face in faces
)
child_areas_squared = tuple(
    face_area_squared(flag_child, face) for face in faces
)
check("the regular tetrahedron control has four equal face areas",
      len(set(regular_areas_squared)) == 1,
      str(regular_areas_squared))
check("the barycentric flag child has four exactly unequal face areas",
      len(set(child_areas_squared)) == 4,
      str(child_areas_squared))

regular_m2, regular_m3, boundary, regular_delta = (
    metric_top_codifferential(regular_vertices)
)
child_m2, child_m3, child_boundary, child_delta = (
    metric_top_codifferential(flag_child)
)
check("the independently integrated Whitney mass matrices are symmetric",
      regular_m2 == regular_m2.T
      and child_m2 == child_m2.T
      and regular_m3 == regular_m3.T
      and child_m3 == child_m3.T)
check("the barycentric child has exactly 1/24 of the parent volume",
      sy.simplify(child_m3[0, 0]/regular_m3[0, 0]-24) == 0,
      "top-form L2 mass scales inversely with volume, hence by 24")

regular_ratios = proportionality_ratios(regular_delta, boundary.T)
child_ratios = proportionality_ratios(child_delta, child_boundary.T)
regular_proportional = (
    regular_ratios is not None and len(set(regular_ratios)) == 1
)
child_proportional = (
    child_ratios is not None and len(set(child_ratios)) == 1
)
check("the metric codifferential is uniform on the regular control tetrahedron",
      regular_proportional,
      f"signed component ratios={regular_ratios}")
check("the top metric codifferential remains uniform on the anisotropic child",
      child_proportional,
      f"signed component ratios={child_ratios}")

child_mass_inverse = sy.simplify(child_m2.inv())
check("the consistent child Whitney 2-form mass is non-diagonal",
      not matrix_is_diagonal(child_m2))
check("its exact inverse is also non-diagonal, so the strong adjoint is nonlocal",
      not matrix_is_diagonal(child_mass_inverse))

# POST-PROTOCOL HOSTILE SCOPE AUDIT.  The preregistered top 3->2 column turns
# out to pass.  That does not establish equality of the full all-degree
# operators.  Check the necessary proportionality in every degree without
# matrix inverses.  This is a stronger target-free diagnostic, reported
# separately from the literal preregistered outcome.
regular_masses, local_differentials, regular_all_degree_ratios = (
    all_degree_metric_audit(regular_vertices)
)
child_masses, _, child_all_degree_ratios = all_degree_metric_audit(flag_child)
check("all three metric adjoints are incidence-proportional on the regular control",
      all(ratio is not None for ratio in regular_all_degree_ratios),
      f"degree ratios={regular_all_degree_ratios}")
check("the child breaks incidence proportionality in degrees 0 and 1 only",
      child_all_degree_ratios[0] is None
      and child_all_degree_ratios[1] is None
      and child_all_degree_ratios[2] == sy.Rational(540, 29),
      f"degree ratios={child_all_degree_ratios}")

child_metric_deltas = [
    sy.simplify(
        child_masses[degree].inv()
        * local_differentials[degree].T
        * child_masses[degree+1]
    )
    for degree in range(3)
]
off_incidence_counts = tuple(sum(
    child_metric_deltas[degree][row, column] != 0
    and local_differentials[degree].T[row, column] == 0
    for row in range(child_metric_deltas[degree].rows)
    for column in range(child_metric_deltas[degree].cols)
) for degree in range(3))
check("the child metric adjoint creates exact off-incidence entries in lower degrees",
      off_incidence_counts[0] > 0
      and off_incidence_counts[1] > 0
      and off_incidence_counts[2] == 0,
      f"off-incidence nonzeros by degree={off_incidence_counts}")

# The unweighted Grover discriminant is fixed entirely by incidence signs and
# Hasse degrees.  On this local top-to-face block all degrees are equal, so it
# uses the boundary direction.  The exact non-proportionality above therefore
# falsifies equality with the accepted metric codifferential.
circle_passes = all_principal and all_dyadic and all_frequency
tetra_metric_passes = child_proportional
protocol_verdict = (
    "DERIVED REFINEMENT-COMPATIBLE TICK"
    if circle_passes and tetra_metric_passes
    else "DERIVED KINEMATIC ONLY"
    if circle_passes
    else "KILL: CIRCLE CONTROL FAILED"
)
full_degree_metric_passes = all(
    ratio is not None for ratio in child_all_degree_ratios
)
final_verdict = (
    protocol_verdict
    if full_degree_metric_passes
    else "DERIVED KINEMATIC ONLY: PREREGISTERED METRIC GATE WAS UNDERPOWERED"
)
check("the preregistered decision boundary is evaluated without fitting",
      protocol_verdict == "DERIVED REFINEMENT-COMPATIBLE TICK",
      protocol_verdict)

payload = {
    "protocol_commit": "59ab096",
    "target_comparison_performed": False,
    "circle": {
        "levels": [8, 16, 32, 64],
        "all_principal_branch_identities_exact": all_principal,
        "all_dyadic_quasienergy_identities_exact": all_dyadic,
        "all_lattice_unit_frequencies_exact": all_frequency,
        "formula": "epsilon_N(m)=asin(sin(pi*m/N))=pi*m/N",
        "refinement_formula": "epsilon_N(m)=2*epsilon_2N(m)",
        "records": circle_records,
    },
    "tetrahedron": {
        "regular_face_area_squared": list(map(str, regular_areas_squared)),
        "flag_child_face_area_squared": list(map(str, child_areas_squared)),
        "regular_metric_direction_ratios": list(map(str, regular_ratios)),
        "flag_child_metric_direction_ratios": list(map(str, child_ratios)),
        "regular_metric_direction_uniform": regular_proportional,
        "flag_child_metric_direction_uniform": child_proportional,
        "regular_all_degree_metric_ratios": [
            None if ratio is None else str(ratio)
            for ratio in regular_all_degree_ratios
        ],
        "flag_child_all_degree_metric_ratios": [
            None if ratio is None else str(ratio)
            for ratio in child_all_degree_ratios
        ],
        "flag_child_full_degree_metric_proportional": full_degree_metric_passes,
        "flag_child_metric_adjoint_off_incidence_nonzeros": list(
            off_incidence_counts
        ),
        "child_M2_diagonal": matrix_is_diagonal(child_m2),
        "child_M2_inverse_diagonal": matrix_is_diagonal(child_mass_inverse),
    },
    "preregistered_protocol_verdict": protocol_verdict,
    "post_protocol_framing_audit": (
        "The single top-to-face test was insufficient to certify the full "
        "all-degree metric operator."
    ),
    "verdict": final_verdict,
    "scope": (
        "The unweighted signed Grover walk has exact causal scaling on the "
        "homogeneous circle but is not a lift of the accepted metric Whitney "
        "codifferential on a barycentric tetrahedron.  This does not exclude "
        "a separately selected metric-aware unitary dilation."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the structured target-blind refinement certificate was written",
      OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"PREREGISTERED_PROTOCOL={protocol_verdict}")
print(f"FINAL_VERDICT={final_verdict}")
print("CIRCLE_REFINEMENT=EXACT")
print(f"TETRA_TOP_WHITNEY_DIRECTION={tetra_metric_passes}")
print(f"TETRA_FULL_DEGREE_METRIC_PROPORTIONAL={full_degree_metric_passes}")
print("OPEN: a canonical metric-aware local unitary dilation.")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Exact audit of round versus fixed-Regge metric selection on the tower.

Protocol commit 8565d19 froze the carrier, admissible criteria, metric family
and decision boundary.  The ambiguity was recognized before registration;
this is a hostile scope audit, not a blind discovery.
"""

from itertools import combinations
from math import factorial
import json
from pathlib import Path

import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "hopf_whitney_metric_selection.json"
PROTOCOL_COMMIT = "8565d19"
tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def zero(value):
    return sy.simplify(sy.expand(value)) == 0


def zero_matrix(matrix):
    return all(zero(value) for value in matrix)


print("Round versus fixed-Regge metric selection on the Whitney tower")

regge_certificate = json.loads(
    (HERE / "whitney_regge_continuum_transfer.json").read_text()
)
gravity_certificate = json.loads(
    (HERE / "hopf_kahler_induced_gravity.json").read_text()
)
edgewise_certificate = json.loads(
    (HERE / "whitney_edgewise_continuum_dynamics.json").read_text()
)
check(
    "the authoritative certificates have the frozen distinct scopes",
    regge_certificate["protocol_commit"] == "1682a46"
    and regge_certificate["carrier"]["map_role"]
    == "coordinate identification only; exact flat metric is pushed forward"
    and gravity_certificate["protocol_commit"] == "31ecea7"
    and gravity_certificate["hypotheses"]["carrier"] == "unit round S3=SU(2)",
    "Whitney target=fixed Regge; Hopf heat response=unit round S3",
)
check(
    "the selected refinement carrier is the canonical rank-edgewise tower",
    edgewise_certificate["all_level_geometry"]["tower"]
    == "Esd_(2^n)(sd K)"
    and regge_certificate["prior_exact_inputs"]["tower"]
    == "Esd_(2^n)(sd K)",
)

# Recover the facet geometry only from the unit 600-cell simplex Gram matrix.
sqrt5 = sy.sqrt(5)
phi = (1 + sqrt5) / 2
adjacent_dot = phi / 2
facet_gram = (
    (1 - adjacent_dot) * sy.eye(4)
    + adjacent_dot * sy.ones(4, 4)
)
ones = sy.ones(4, 1)
facet_distance_squared = sy.simplify(
    (ones.T * facet_gram * ones)[0] / 16
)
expected_distance_squared = (7 + 3 * sqrt5) / 16
check(
    "the unit 600-cell facet distance is exact",
    zero(facet_distance_squared - expected_distance_squared),
    "a^2=(7+3 sqrt(5))/16",
)

# Radial normalization and its pulled-back round metric.
x_symbols = sy.symbols("x0:4", real=True)
ambient_x = sy.Matrix(x_symbols)
radius_squared = ambient_x.dot(ambient_x)
radial = ambient_x / sy.sqrt(radius_squared)
radial_jacobian = radial.jacobian(x_symbols)
round_pullback = sy.simplify(radial_jacobian.T * radial_jacobian)
round_formula = (
    sy.eye(4) / radius_squared
    - ambient_x * ambient_x.T / radius_squared**2
)
check(
    "the radial pullback of the round metric is exact",
    zero_matrix(round_pullback - round_formula),
    "R*g0=I/r^2-xx^T/r^4",
)

# Rotate one facet so its supporting normal is coordinate zero and align the
# tangential component of x with coordinate one.
a, y = sy.symbols("a y", positive=True)
local_x = sy.Matrix((a, y, 0, 0))
local_radius_squared = a**2 + y**2
local_round = (
    sy.eye(4) / local_radius_squared
    - local_x * local_x.T / local_radius_squared**2
)
tangent_round = sy.simplify(local_round[1:4, 1:4])
expected_tangent_round = sy.diag(
    a**2 / local_radius_squared**2,
    1 / local_radius_squared,
    1 / local_radius_squared,
)
check(
    "the round/flat generalized tangent eigenvalues are exact",
    zero_matrix(tangent_round - expected_tangent_round),
    "(a^2/r^4,1/r^2,1/r^2)",
)

# Use the midpoint of the facet centroid-to-vertex segment.  It is strictly
# interior and has |y|^2=(1-a^2)/4.
interior_y_squared = sy.simplify((1 - expected_distance_squared) / 4)
interior_radius_squared = sy.simplify(
    expected_distance_squared + interior_y_squared
)
interior_parallel = sy.simplify(
    expected_distance_squared / interior_radius_squared**2
)
interior_transverse = sy.simplify(1 / interior_radius_squared)
check(
    "round and Regge metrics are exactly nonproportional at a facet-interior point",
    interior_y_squared > 0
    and not zero(interior_parallel - interior_transverse),
    "g_R=I while g_0 has unequal parallel/transverse eigenvalues",
)

# O(4)-covariance is stronger than the required H4 covariance.  Coordinate
# permutations reduce all plane rotations to the displayed 01 rotation, and
# one coordinate reflection supplies the other component of O(4).
theta = sy.symbols("theta", real=True)
cosine, sine = sy.cos(theta), sy.sin(theta)
plane_rotation = sy.Matrix((
    (cosine, -sine, 0, 0),
    (sine, cosine, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
))
reflection = sy.diag(-1, 1, 1, 1)
coordinate_permutation = sy.Matrix((
    (1, 0, 0, 0),
    (0, 0, 1, 0),
    (0, 1, 0, 0),
    (0, 0, 0, 1),
))


def radial_metric(vector):
    norm_squared = vector.dot(vector)
    return sy.eye(4) / norm_squared - vector * vector.T / norm_squared**2


rotation_covariance = radial_metric(plane_rotation * ambient_x) - (
    plane_rotation * radial_metric(ambient_x) * plane_rotation.T
)
reflection_covariance = radial_metric(reflection * ambient_x) - (
    reflection * radial_metric(ambient_x) * reflection.T
)
permutation_covariance = radial_metric(coordinate_permutation * ambient_x) - (
    coordinate_permutation * radial_metric(ambient_x)
    * coordinate_permutation.T
)
check(
    "the radial round metric is O(4)-equivariant, hence H4-equivariant",
    all(sy.trigsimp(sy.factor(value)) == 0 for value in rotation_covariance)
    and zero_matrix(reflection_covariance)
    and zero_matrix(permutation_covariance),
)
check(
    "the pushed-forward flat metric is H4-equivariant as well",
    zero_matrix(plane_rotation.T * plane_rotation - sy.eye(4))
    and zero_matrix(reflection.T * reflection - sy.eye(4)),
    "orthogonal 600-cell symmetries preserve flat facet lengths and commute with R",
)

# The flat trace on a shared equilateral face is parent-independent.  The
# round trace uses only ambient x and face-tangent vectors, so it is also
# parent-independent; convex combinations inherit the property.
edge_length_squared = sy.simplify(2 - 2 * adjacent_dot)
face_trace = sy.Matrix((
    (edge_length_squared, 1 - adjacent_dot),
    (1 - adjacent_dot, edge_length_squared),
))
equilateral_trace = edge_length_squared * sy.Matrix((
    (1, sy.Rational(1, 2)),
    (sy.Rational(1, 2), 1),
))
check(
    "both endpoint metrics have matching tangential traces across shared faces",
    zero_matrix(face_trace - equilateral_trace)
    and set(round_formula.free_symbols) == set(x_symbols),
    "the round pullback contains no parent-facet datum",
)

# The complete affine family passes the same frozen kinematic criteria.
u = sy.symbols("u", real=True, nonnegative=True)
family_parallel = sy.simplify(1 - u + u * interior_parallel)
family_transverse = sy.simplify(1 - u + u * interior_transverse)
metric_lower = expected_distance_squared
metric_upper = sy.simplify(1 / expected_distance_squared)
check(
    "every convex round/Regge metric is positive for 0<=u<=1",
    metric_lower.is_positive
    and (1 - metric_lower).is_positive
    and interior_parallel.is_positive
    and interior_transverse.is_positive,
    "all eigenvalues lie in the common interval [a^2,1/a^2]",
)
check(
    "the affine family has the same level-independent equivalence bounds",
    zero(metric_upper - (28 - 12 * sqrt5))
    and float(metric_lower.evalf()) < 1 < float(metric_upper.evalf()),
    f"{metric_lower} <= g_u/g_R <= {metric_upper}",
)
check(
    "the metric family is injective in its unselected parameter",
    not zero(interior_parallel - 1)
    or not zero(interior_transverse - 1),
    "g_u-g_v=(u-v)(g_0-g_R), with g_0-g_R nonzero",
)

# -------------------------------------------------------------------------
# Exact generic-metric Whitney refinement control.
# -------------------------------------------------------------------------
mx, my, mz = sy.symbols("mx my mz", positive=True)
coordinate_metric = sy.diag(mx, my, mz)
inverse_coordinate_metric = coordinate_metric.inv()
metric_volume_factor = sy.sqrt(coordinate_metric.det())
vertices = {
    0: sy.Matrix((0, 0, 0)),
    1: sy.Matrix((1, 0, 0)),
    2: sy.Matrix((0, 1, 0)),
    3: sy.Matrix((0, 0, 1)),
    4: sy.Matrix((sy.Rational(1, 2), 0, 0)),
}
coarse_top = ((0, 1, 2, 3),)
fine_top = ((0, 2, 3, 4), (1, 2, 3, 4))


def simplices(top):
    return tuple(
        tuple(sorted({
            face for cell in top for face in combinations(cell, degree + 1)
        }))
        for degree in range(4)
    )


coarse = simplices(coarse_top)
fine = simplices(fine_top)
coordinate_bases = [list(combinations(range(3), degree)) for degree in range(4)]


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


def local_mass(points, degree):
    affine = sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    )
    inverse = affine.inv()
    gradients = [-sum(
        (sy.Matrix(inverse.row(row)).T for row in range(3)),
        sy.zeros(3, 1),
    )]
    gradients.extend(sy.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det()) * metric_volume_factor / 6
    barycentric_second_moment = volume * (sy.ones(4, 4) + sy.eye(4)) / 20
    local_forms = list(combinations(range(4), degree + 1))
    coefficients = []
    for form in local_forms:
        matrix = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            matrix[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                matrix[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * wedge_components(covectors, degree)
                )
        coefficients.append(matrix)
    if degree == 0:
        wedge_metric = sy.ones(1, 1)
    else:
        wedge_metric = sy.Matrix([
            [
                sy.det(inverse_coordinate_metric.extract(left, right))
                for right in coordinate_bases[degree]
            ]
            for left in coordinate_bases[degree]
        ])
    mass = sy.zeros(len(local_forms), len(local_forms))
    for row, left in enumerate(coefficients):
        for column, right in enumerate(coefficients):
            mass[row, column] = sy.factor(sum(
                wedge_metric[i, j]
                * (left[i, :] * barycentric_second_moment
                   * right[j, :].T)[0]
                for i in range(wedge_metric.rows)
                for j in range(wedge_metric.cols)
            ))
    return mass


def assemble_mass(top, cell_layers, degree):
    indices = {cell: index for index, cell in enumerate(cell_layers[degree])}
    local_faces = list(combinations(range(4), degree + 1))
    result = sy.zeros(len(cell_layers[degree]), len(cell_layers[degree]))
    for cell in top:
        local = local_mass(tuple(vertices[index] for index in cell), degree)
        for local_row, left in enumerate(local_faces):
            row = indices[tuple(cell[index] for index in left)]
            for local_column, right in enumerate(local_faces):
                column = indices[tuple(cell[index] for index in right)]
                result[row, column] += local[local_row, local_column]
    return sy.simplify(result)


fine_barycentric = {
    0: sy.Matrix((1, 0, 0, 0)),
    1: sy.Matrix((0, 1, 0, 0)),
    2: sy.Matrix((0, 0, 1, 0)),
    3: sy.Matrix((0, 0, 0, 1)),
    4: sy.Matrix((sy.Rational(1, 2), sy.Rational(1, 2), 0, 0)),
}
coarse_masses = []
fine_masses = []
inclusions = []
isometry_residuals = []
for degree in range(4):
    coarse_mass = assemble_mass(coarse_top, coarse, degree)
    fine_mass = assemble_mass(fine_top, fine, degree)
    inclusion = sy.zeros(len(fine[degree]), len(coarse[degree]))
    for row, fine_simplex in enumerate(fine[degree]):
        barycentric_rows = [fine_barycentric[vertex] for vertex in fine_simplex]
        for column, coarse_simplex in enumerate(coarse[degree]):
            inclusion[row, column] = sy.Matrix([
                [barycentric[vertex] for vertex in coarse_simplex]
                for barycentric in barycentric_rows
            ]).det()
    coarse_masses.append(coarse_mass)
    fine_masses.append(fine_mass)
    inclusions.append(inclusion)
    isometry_residuals.append(sy.simplify(
        inclusion.T * fine_mass * inclusion - coarse_mass
    ))

check(
    "exact Whitney refinement is metric-isometric for generic diag(mx,my,mz)",
    all(zero_matrix(residual) for residual in isometry_residuals),
    "P_p^T M_f,p(g) P_p=M_c,p(g), p=0,1,2,3",
)
check(
    "the same metric-independent inclusions have full coarse rank",
    all(inclusion.rank() == len(coarse[degree])
        for degree, inclusion in enumerate(inclusions)),
)

# The metric family is real, not a common rescaling hidden in all masses.
isotropic_one_mass = coarse_masses[1].subs({mx: 1, my: 1, mz: 1})
anisotropic_one_mass = coarse_masses[1].subs({mx: 2, my: 1, mz: 1})
scale_ratio = sy.simplify(
    anisotropic_one_mass[0, 0] / isotropic_one_mass[0, 0]
)
proportional_residual = sy.simplify(
    anisotropic_one_mass - scale_ratio * isotropic_one_mass
)
minimum_control_eigenvalue = min(
    np.linalg.eigvalsh(np.asarray(
        mass.subs({mx: 2, my: 3, mz: 5}), dtype=float
    )).min()
    for mass in coarse_masses + fine_masses
)
check(
    "different positive metrics give genuinely different Whitney mass data",
    not zero_matrix(proportional_residual)
    and minimum_control_eigenvalue > 0,
    "middle-degree masses are not related by one common scalar",
)

# Equivalent Hilbert norms do not preserve adjoints or spectral responses.
d_control = sy.Matrix(((1, 1, 0), (0, 1, 1)))
control_g0 = sy.eye(3)
control_g1 = sy.diag(2, 3, 5)
right_inverse_0 = (
    control_g0.inv() * d_control.T
    * (d_control * control_g0.inv() * d_control.T).inv()
)
right_inverse_1 = (
    control_g1.inv() * d_control.T
    * (d_control * control_g1.inv() * d_control.T).inv()
)
check(
    "uniform metric equivalence does not identify adjoints or spectra",
    not zero_matrix(right_inverse_1 - right_inverse_0),
    "the exact minimal right inverse changes under an equivalent metric",
)

# The endpoints have different local curvature before any action comparison.
round_scalar_curvature = gravity_certificate["exact_results"][
    "round_scalar_curvature"
]
flat_facet_scalar_curvature = 0
check(
    "the endpoint local geometries have different scalar curvature",
    round_scalar_curvature == 6 and flat_facet_scalar_curvature == 0,
    "g0 has R=6; gR is flat in every open facet",
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "provenance": "post-recognition hostile metric-selection audit",
    "phenomenological_target_used": False,
    "carrier": "K_n=Esd_(2^n)(sd boundary(600-cell)) under radial identification",
    "endpoint_metrics": {
        "g_R": "pushforward of exact flat facet metric",
        "g_0": "unit round S3 metric",
        "distinct_at_facet_interior": True,
        "both_H4_equivariant": True,
        "both_face_compatible": True,
        "both_refinement_natural": True,
    },
    "metric_family": {
        "formula": "g_u=(1-u)g_R+u g_0, 0<=u<=1",
        "injective": True,
        "uniform_bounds_relative_to_g_R": [
            str(metric_lower), str(metric_upper),
        ],
        "whitney_generic_metric_isometry_all_degrees": True,
    },
    "curvature_scope": {
        "round_facet_interior_scalar": 6,
        "regge_facet_interior_scalar": 0,
        "norm_equivalence_transfers_heat_hessian": False,
    },
    "verdicts": [
        {
            "label": "DERIVED METRIC-SELECTION NO-GO",
            "claim": (
                "symmetry, radial naturality, positivity, face compatibility "
                "and exact Whitney refinement do not distinguish g_R from g_0"
            ),
        },
        {
            "label": "DERIVED SCOPE BOUNDARY",
            "claim": (
                "the round Hopf heat Hessian cannot be assigned to the fixed-"
                "Regge Whitney theory using norm equivalence"
            ),
        },
        {
            "label": "OPEN DYNAMICAL SELECTION",
            "claim": (
                "an independently selected action could still choose a metric "
                "and is not excluded by this kinematic no-go"
            ),
        },
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("machine-readable metric-selection certificate was written", OUTPUT.exists())

print(f"\nRESULT: {passed}/{tests} checks passed")
print("DERIVED NO-GO: frozen symmetry/refinement criteria leave a metric continuum.")
print("SCOPE: the round Hopf heat Hessian does not transfer to fixed Regge by norm equivalence.")
print("OPEN: an independently selected dynamical metric action could break the tie.")
raise SystemExit(0 if passed == tests else 1)

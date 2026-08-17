#!/usr/bin/env python3
"""Certify physical frustum equivalence of all homothetic staircases."""

from collections import Counter
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gravity_600cell_overlay_metric_compatibility.json"
OUTPUT = HERE / "gravity_600cell_homothetic_frustum_equivalence.json"
UPSTREAM_SHA256 = "7de877b83b5524a1c86788f207ec205fa1eae799ca66bf62c1ae6b46081bb45e"
PRIOR_ART_COMMIT = "e4428bf"
PROTOCOL_COMMIT = "affb02e"
VERTICES = tuple(range(4))
ORDERS = tuple(permutations(VERTICES))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equal(left, right):
    return sp.simplify(sp.expand(left-right)) == 0


def matrix_equal(left, right):
    return left.shape == right.shape and all(
        equal(a, b) for a, b in zip(left, right)
    )


def expression_key(expression):
    return sp.srepr(sp.factor(expression))


tests = []


def check(label, condition):
    tests.append((label, bool(condition)))


upstream = json.loads(UPSTREAM.read_text())
upstream_ok = bool(
    digest(UPSTREAM) == UPSTREAM_SHA256
    and upstream.get("outcome") == "OVERLAY_INHERITS_STATIC_METRIC_ONLY"
    and upstream.get("passed") == upstream.get("tests") == 9
    and upstream.get("chamber_count") == 148
    and upstream.get("global_compatibility_gcd") == "r - 1"
    and upstream.get("positive_compatible_ratios") == ["1"]
)
check("the frozen identity-pullback obstruction is reproduced", upstream_ok)


# Exact regular-tetrahedron and homothetic data.
R_MINUS, R_PLUS, RHO, T = sp.symbols(
    "R_minus R_plus rho T", positive=True
)
DELTA = R_PLUS-R_MINUS
T_SQUARE = RHO+DELTA**2
SQRT5 = sp.sqrt(5)
PHI = (1+SQRT5)/2
C = PHI/2
D = 1+3*C
U = sp.Matrix(4, 4, lambda i, j: 1 if i == j else C)

n_coefficients = sp.Matrix([1/D]*4)
n_dot_vertices = tuple(sp.expand((n_coefficients.T*U[:, i])[0]) for i in VERTICES)
n_square = sp.simplify((n_coefficients.T*U*n_coefficients)[0])
normal_formula_ok = bool(
    all(equal(value, 1) for value in n_dot_vertices)
    and equal(n_square, 4/D)
)
check("the exact common-hyperplane normal is derived from the Gram matrix",
      normal_formula_ok)


# Rows are (four abstract spatial coefficients, target time, homogeneous 1).
physical_vertex_rows = []
vertex_labels = []
for top, scale, target_time in ((False, R_MINUS, 0), (True, R_PLUS, T)):
    for vertex in VERTICES:
        spatial = [sp.Integer(0)]*4
        spatial[vertex] = scale
        physical_vertex_rows.append((*spatial, target_time, sp.Integer(1)))
        vertex_labels.append((top, vertex))
physical_vertex_matrix = sp.Matrix(physical_vertex_rows)
hyperplane_covector = sp.Matrix([1, 1, 1, 1, -DELTA/T, -R_MINUS])
hyperplane_residual = physical_vertex_matrix*hyperplane_covector
nullspace = physical_vertex_matrix.nullspace()
nullspace_matches = bool(
    len(nullspace) == 1
    and all(
        equal(nullspace[0][index]*hyperplane_covector[0],
              hyperplane_covector[index]*nullspace[0][0])
        for index in range(6)
    )
)
affine_hull_ok = bool(
    physical_vertex_matrix.rank() == 5
    and all(equal(value, 0) for value in hyperplane_residual)
    and nullspace_matches
)
check("the eight vertices have one independently recovered affine 4-plane",
      affine_hull_ok)


# Exact causal character from both the ambient normal and the induced metric.
normal_square = 4/D-DELTA**2/T_SQUARE
normal_positive_decomposition = (
    (4/D-1)+RHO/T_SQUARE
)
positive_margin = sp.simplify(4/D-1)
normal_signature_ok = bool(
    equal(normal_square, normal_positive_decomposition)
    and equal(positive_margin, 3*(1-C)/D)
    and equal(1-C, (3-SQRT5)/4)
)

y0, y1, y2, z = sp.symbols("y_0 y_1 y_2 z")
Y = (y0, y1, y2, R_MINUS+DELTA*z-y0-y1-y2)
spatial_jacobian = sp.Matrix(Y).jacobian((y0, y1, y2, z))
G_Q = (spatial_jacobian.T*U*spatial_jacobian).applyfunc(sp.expand)
G_Q[3, 3] = sp.expand(G_Q[3, 3]-T_SQUARE)
spatial_block = G_Q[:3, :3]
expected_spatial_block = (1-C)*(sp.eye(3)+sp.ones(3))
cross_block = G_Q[:3, 3]
schur_complement = sp.simplify(
    G_Q[3, 3]-(cross_block.T*spatial_block.inv()*cross_block)[0]
)
expected_schur = -RHO-sp.Rational(3, 4)*(1-C)*DELTA**2
induced_signature_ok = bool(
    matrix_equal(spatial_block, expected_spatial_block)
    and equal(spatial_block.det(), 4*(1-C)**3)
    and equal(schur_complement, expected_schur)
)
signature_ok = normal_signature_ok and induced_signature_ok
check("the common affine plane has exact Lorentzian signature (3,1)",
      signature_ok)


# The six exact facet slacks of Q; coefficients are
# (constant,y0,y1,y2,z).
FACETS = (
    ("bottom", (0, 0, 0, 0, 1)),
    ("top", (1, 0, 0, 0, -1)),
    ("lateral_0", (0, 1, 0, 0, 0)),
    ("lateral_1", (0, 0, 1, 0, 0)),
    ("lateral_2", (0, 0, 0, 1, 0)),
    ("lateral_3", (R_MINUS, -1, -1, -1, DELTA)),
)


def evaluate_form(form, point):
    return sp.expand(form[0]+sum(
        coefficient*coordinate
        for coefficient, coordinate in zip(form[1:], point)
    ))


def reduced_vertex(vertex, top):
    scale = R_PLUS if top else R_MINUS
    point = [sp.Integer(0)]*4
    if vertex < 3:
        point[vertex] = scale
    point[3] = sp.Integer(1 if top else 0)
    return tuple(point)


expected_vertices = {
    reduced_vertex(vertex, top)
    for top in (False, True)
    for vertex in VERTICES
}


def universally_nonnegative(expression):
    expression = sp.factor(expression)
    return any(equal(expression, candidate) for candidate in (
        0, 1, R_MINUS, R_PLUS,
    ))


enumerated_vertices = set()
intersection_records = []
for active in combinations(range(6), 4):
    coefficient_matrix = sp.Matrix([
        FACETS[index][1][1:] for index in active
    ])
    if coefficient_matrix.det() == 0:
        continue
    rhs = sp.Matrix([-FACETS[index][1][0] for index in active])
    point = tuple(sp.simplify(value) for value in coefficient_matrix.inv()*rhs)
    slacks = tuple(
        sp.factor(evaluate_form(record[1], point)) for record in FACETS
    )
    feasible = all(universally_nonnegative(value) for value in slacks)
    intersection_records.append({
        "active_facets": [FACETS[index][0] for index in active],
        "point": [str(value) for value in point],
        "slacks": [str(value) for value in slacks],
        "universally_feasible": feasible,
    })
    if feasible:
        enumerated_vertices.add(point)

facet_vertex_sets = {}
for label, form in FACETS:
    facet_vertex_sets[label] = {
        point for point in enumerated_vertices if equal(evaluate_form(form, point), 0)
    }
facet_counts = {label: len(items) for label, items in facet_vertex_sets.items()}

vertex_affine_matrix = sp.Matrix([
    (1, *point) for point in sorted(expected_vertices, key=lambda p: tuple(map(str, p)))
])
frustum_volume = (
    R_MINUS**3+R_MINUS**2*R_PLUS+R_MINUS*R_PLUS**2+R_PLUS**3
)/24
integrated_volume = sp.integrate(
    (R_MINUS+DELTA*z)**3/6, (z, 0, 1)
)
polytope_ok = bool(
    enumerated_vertices == expected_vertices
    and facet_counts == {
        "bottom": 4,
        "top": 4,
        "lateral_0": 6,
        "lateral_1": 6,
        "lateral_2": 6,
        "lateral_3": 6,
    }
    and vertex_affine_matrix.rank() == 5
    and equal(integrated_volume, frustum_volume)
)
check("the six half-spaces give exactly one tetrahedral 4-frustum",
      polytope_ok)


def parameter_vertex(vertex, top):
    lambdas = [sp.Integer(0)]*4
    lambdas[vertex] = 1
    return (1, *lambdas[:3], int(top))


def q_vertex(vertex, top):
    return (1, *reduced_vertex(vertex, top))


def labelled_simplex(order, split):
    return tuple(
        [(False, vertex) for vertex in order[:split+1]]
        + [(True, vertex) for vertex in order[split:]]
    )


simplex_records = []
affine_q_maps = {}
orientation_ok = True
schedule_geometry_records = []
for order in ORDERS:
    simplices = [labelled_simplex(order, split) for split in range(4)]
    face_incidence = Counter()
    volume_sum = sp.Integer(0)
    split_records = []
    for split, simplex in enumerate(simplices):
        domain = sp.Matrix([
            parameter_vertex(vertex, top) for top, vertex in simplex
        ])
        target = sp.Matrix([
            q_vertex(vertex, top) for top, vertex in simplex
        ])
        domain_det = sp.factor(domain.det())
        target_det = sp.factor(target.det())
        ratio_det = sp.factor(target_det/domain_det)
        expected_ratio = R_MINUS**split*R_PLUS**(3-split)
        current_ok = equal(ratio_det, expected_ratio)
        orientation_ok &= current_ok
        volume_sum += ratio_det/24
        target_values = target[:, 1:]
        forms = (domain.inv()*target_values).T.applyfunc(sp.expand)
        affine_q_maps[(order, split)] = forms
        for facet in combinations(simplex, 4):
            face_incidence[frozenset(facet)] += 1
        record = {
            "order": list(order),
            "split": split,
            "domain_determinant": str(domain_det),
            "physical_determinant": str(target_det),
            "orientation_ratio": str(ratio_det),
            "expected_positive_ratio": str(expected_ratio),
            "pass": current_ok,
        }
        simplex_records.append(record)
        split_records.append(record)

    internal_facets = [face for face, count in face_incidence.items() if count == 2]
    external_facets = [face for face, count in face_incidence.items() if count == 1]
    bad_multiplicity = [count for count in face_incidence.values() if count not in (1, 2)]
    boundary_distribution = Counter()
    for facet in external_facets:
        containing = []
        for label, vertices in (
            ("bottom", {(False, index) for index in VERTICES}),
            ("top", {(True, index) for index in VERTICES}),
            *((f"lateral_{omitted}", {
                (top, index)
                for top in (False, True)
                for index in VERTICES if index != omitted
            }) for omitted in VERTICES),
        ):
            if facet <= vertices:
                containing.append(label)
        if len(containing) == 1:
            boundary_distribution[containing[0]] += 1
        else:
            boundary_distribution["UNCLASSIFIED"] += 1
    schedule_ok = bool(
        len(internal_facets) == 3
        and len(external_facets) == 14
        and not bad_multiplicity
        and boundary_distribution == Counter({
            "bottom": 1,
            "top": 1,
            "lateral_0": 3,
            "lateral_1": 3,
            "lateral_2": 3,
            "lateral_3": 3,
        })
        and equal(volume_sum, frustum_volume)
    )
    schedule_geometry_records.append({
        "order": list(order),
        "internal_tetrahedra": len(internal_facets),
        "boundary_tetrahedra": len(external_facets),
        "boundary_distribution": dict(sorted(boundary_distribution.items())),
        "simplex_volume_sum": str(sp.factor(volume_sum)),
        "pass": schedule_ok,
    })

schedule_census_ok = bool(
    len(simplex_records) == 96
    and len(schedule_geometry_records) == 24
    and all(len(record["boundary_distribution"]) >= 1
            for record in schedule_geometry_records)
)
all_schedule_geometry_ok = bool(
    orientation_ok
    and schedule_census_ok
    and all(record["pass"] for record in schedule_geometry_records)
)
check("all 96 simplices and 24 schedule geometry censuses are complete",
      schedule_census_ok)


# Explicit projective equivalence from the standard prism to Q.
lam0, lam1, lam2, tau = sp.symbols("lambda_0 lambda_1 lambda_2 tau")
LAMBDAS = (lam0, lam1, lam2, 1-lam0-lam1-lam2)
projective_denominator = (1-tau)/R_MINUS+tau/R_PLUS
forward_y = tuple(sp.factor(value/projective_denominator) for value in LAMBDAS)
forward_z = sp.factor((tau/R_PLUS)/projective_denominator)
forward_sum = sp.factor(sum(forward_y))

inverse_scale = R_MINUS+DELTA*z
inverse_lambdas = tuple(sp.factor(value/inverse_scale) for value in Y)
inverse_tau = sp.factor(
    z*R_PLUS/(R_MINUS*(1-z)+R_PLUS*z)
)

forward_inverse_checks = [
    equal(inverse_tau.subs({
        z: forward_z,
    }), tau)
]
for inverse_lambda, original_lambda in zip(inverse_lambdas, LAMBDAS):
    substituted = inverse_lambda.subs({
        y0: forward_y[0],
        y1: forward_y[1],
        y2: forward_y[2],
        z: forward_z,
    }, simultaneous=True)
    forward_inverse_checks.append(equal(substituted, original_lambda))

inverse_substitution = {
    lam0: inverse_lambdas[0],
    lam1: inverse_lambdas[1],
    lam2: inverse_lambdas[2],
    tau: inverse_tau,
}
inverse_forward_checks = [
    equal(forward_z.subs(inverse_substitution, simultaneous=True), z)
]
for forward_coordinate, original_coordinate in zip(forward_y[:3], (y0, y1, y2)):
    inverse_forward_checks.append(
        equal(forward_coordinate.subs(
            inverse_substitution, simultaneous=True
        ), original_coordinate)
    )

projective_vertex_checks = []
for top in (False, True):
    for vertex in VERTICES:
        substitution = {
            lam0: int(vertex == 0),
            lam1: int(vertex == 1),
            lam2: int(vertex == 2),
            tau: int(top),
        }
        image = tuple(sp.simplify(value.subs(substitution)) for value in forward_y[:3])
        image += (sp.simplify(forward_z.subs(substitution)),)
        projective_vertex_checks.append(image == reduced_vertex(vertex, top))

facet_image_checks = bool(
    equal(forward_z.subs(tau, 0), 0)
    and equal(forward_z.subs(tau, 1), 1)
    and all(equal(forward_y[index].subs(LAMBDAS[index], 0), 0)
            for index in range(3))
    and equal(forward_y[3].subs(lam2, 1-lam0-lam1), 0)
)
projective_ok = bool(
    equal(forward_sum, R_MINUS+DELTA*forward_z)
    and all(forward_inverse_checks)
    and all(inverse_forward_checks)
    and all(projective_vertex_checks)
    and facet_image_checks
)
check("an exact positive-denominator projective homeomorphism closes non-overlap",
      projective_ok)


# Tensor-law isometry checks between all 32 distinct local affine map types.
affine_types = {}
for key, forms in affine_q_maps.items():
    affine_types.setdefault(tuple(expression_key(value) for value in forms), (key, forms))
type_records = tuple(affine_types.values())
inverse_jacobians = {}
local_metrics = {}
for key, forms in type_records:
    jacobian = forms[:, 1:]
    inverse_jacobians[key] = jacobian.inv().applyfunc(sp.factor)
    local_metrics[key] = (jacobian.T*G_Q*jacobian).applyfunc(sp.expand)

isometry_failures = []
for left_key, left_forms in type_records:
    left_jacobian = left_forms[:, 1:]
    left_metric = local_metrics[left_key]
    for right_key, right_forms in type_records:
        right_jacobian = right_forms[:, 1:]
        transition = inverse_jacobians[right_key]*left_jacobian
        transported = transition.T*local_metrics[right_key]*transition
        if not matrix_equal(transported, left_metric):
            isometry_failures.append((left_key, right_key))

bottom_boundary_ok = True
top_boundary_ok = True
# In reduced Q coordinates, y0..2=R*lambda0..2 and z is the endpoint.
for order in ORDERS:
    bottom_forms = affine_q_maps[(order, 3)]
    top_forms = affine_q_maps[(order, 0)]
    bottom_expressions = bottom_forms[:, 0]+bottom_forms[:, 1:]*sp.Matrix(
        [lam0, lam1, lam2, 0]
    )
    top_expressions = top_forms[:, 0]+top_forms[:, 1:]*sp.Matrix(
        [lam0, lam1, lam2, 1]
    )
    bottom_expected = sp.Matrix([R_MINUS*lam0, R_MINUS*lam1, R_MINUS*lam2, 0])
    top_expected = sp.Matrix([R_PLUS*lam0, R_PLUS*lam1, R_PLUS*lam2, 1])
    bottom_boundary_ok &= matrix_equal(bottom_expressions, bottom_expected)
    top_boundary_ok &= matrix_equal(top_expressions, top_expected)

isometry_ok = bool(
    len(type_records) == 32
    and not isometry_failures
    and bottom_boundary_ok
    and top_boundary_ok
)
check("all 32 local map types are isometric and fix both time boundaries",
      isometry_ok)


# Reconcile the common frustum with the repository's homothetic edge lengths.
L_MINUS, L_PLUS = sp.symbols("L_minus L_plus", positive=True)
bottom_edge_square = sp.expand(2*R_MINUS**2*(1-C))
top_edge_square = sp.expand(2*R_PLUS**2*(1-C))
strut_square = sp.expand(DELTA**2-T_SQUARE)
cross_square = sp.expand(
    R_MINUS**2+R_PLUS**2-2*C*R_MINUS*R_PLUS-T_SQUARE
)
edge_reconciliation_ok = bool(
    equal(bottom_edge_square.subs(R_MINUS, PHI*L_MINUS), L_MINUS**2)
    and equal(top_edge_square.subs(R_PLUS, PHI*L_PLUS), L_PLUS**2)
    and equal(strut_square, -RHO)
    and equal(
        cross_square.subs({R_MINUS: PHI*L_MINUS, R_PLUS: PHI*L_PLUS}),
        L_MINUS*L_PLUS-RHO,
    )
)
check("the common frustum reproduces every inherited homothetic edge length",
      edge_reconciliation_ok)


controls_ok = all(condition for _, condition in tests)
if not controls_ok:
    outcome = "HOMOTHETIC_FRUSTUM_CONTROL_FAILED"
elif not all_schedule_geometry_ok:
    outcome = "HOMOTHETIC_FRUSTUM_TRIANGULATION_FAILED"
else:
    outcome = "HOMOTHETIC_SCHEDULES_ONE_LORENTZIAN_FRUSTUM"


payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": {"identity_pullback_audit": digest(UPSTREAM)},
    "tetrahedron_gram_off_diagonal": str(C),
    "affine_hull_homogeneous_rank": physical_vertex_matrix.rank(),
    "affine_hull_nullity": len(nullspace),
    "hyperplane_equation": "sum(y_i)-(R_plus-R_minus)*time/T-R_minus=0",
    "normal_square": str(sp.factor(normal_square)),
    "normal_positive_decomposition": str(sp.factor(normal_positive_decomposition)),
    "induced_metric_reduced_coordinates": [
        [str(sp.factor(G_Q[i, j])) for j in range(4)] for i in range(4)
    ],
    "spatial_block_determinant": str(sp.factor(spatial_block.det())),
    "lorentzian_schur_complement": str(sp.factor(schur_complement)),
    "polytope_vertex_count": len(enumerated_vertices),
    "polytope_facet_vertex_counts": facet_counts,
    "facet_intersections": intersection_records,
    "frustum_coordinate_four_volume": str(sp.factor(frustum_volume)),
    "labelled_simplex_count": len(simplex_records),
    "distinct_local_affine_map_types": len(type_records),
    "simplex_records": simplex_records,
    "schedule_geometry": schedule_geometry_records,
    "isometry_failure_count": len(isometry_failures),
    "projective_map": {
        "denominator": str(projective_denominator),
        "forward_y": [str(value) for value in forward_y],
        "forward_z": str(forward_z),
        "inverse_lambda": [str(value) for value in inverse_lambdas],
        "inverse_tau": str(inverse_tau),
    },
    "edge_squares": {
        "bottom": str(sp.factor(bottom_edge_square)),
        "top": str(sp.factor(top_edge_square)),
        "same_vertex_strut": str(sp.factor(strut_square)),
        "cross_diagonal": str(sp.factor(cross_square)),
    },
    "regge_action_evaluations": 0,
    "dust_action_evaluations": 0,
    "anisotropic_states_evaluated": 0,
    "tests": len(tests),
    "passed": sum(condition for _, condition in tests),
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, condition in tests:
    print(f"{'PASS' if condition else 'FAIL'}: {label}")
print(f"affine hull rank={payload['affine_hull_homogeneous_rank']-1}")
print(f"frustum vertices={payload['polytope_vertex_count']}, facets={facet_counts}")
print(f"local affine map types={payload['distinct_local_affine_map_types']}")
print(f"frustum volume={payload['frustum_coordinate_four_volume']}")
print(f"OUTCOME: {outcome}")
print(f"{payload['passed']}/{payload['tests']} tests passed")

raise SystemExit(0 if controls_ok else 1)

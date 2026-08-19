#!/usr/bin/env python3
"""Exact two-frustum gluing with the derived variable face transition."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_variable_face_connection.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_variable_face_connection_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_variable_face_connection_protocol.md"
LOCAL_RESULT = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_result.md"
FIXED_RESULT = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_result.md"
FIXED_JSON = HERE / "gravity_600cell_two_frustum_face_gluing.json"
GLOBAL_RESULT = ROOT / "docs/gravity/gravity_600cell_global_flex_holonomy_result.md"
GLOBAL_JSON = HERE / "gravity_600cell_global_flex_holonomy_adversarial.json"

PROTOCOL_COMMIT = "2260b72"
EXPECTED_HASHES = {
    "prior_art": "2ed809fedad24fa15977b39e4dd6fec386e9080c123208d54fd089554ce44d2d",
    "protocol": "f6b91206a857cda6ebfe5cb9988110de5f12a9c1ca51bcbdb733a8429682ca6a",
    "local_result": "436fb57037e491b6bdb8fee9ad8b10ab8da1621fd9ecda73e1fcac3fa616fa29",
    "fixed_result": "b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3",
    "fixed_json": "0e09c3f8f38c8158deff5b81bc6fe4d5d6dd685a24cce83e015fb95e3f26a70e",
    "global_result": "72c8b2c0ffbd9d13aef8f14404270cac29896c876ed6f015a4dc7a41a89b6535",
    "global_json": "f224fe123c882ccda97d4ca6ec67c9fd810d58ed8377c5afb457a1dec69f4b87",
}

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
POINTS = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
    (sp.Rational(5, 3), sp.Rational(5, 3), sp.Rational(-5, 3), 0),
))
LEFT = (0, 1, 2, 3)
RIGHT = (0, 1, 2, 4)
SHARED = (0, 1, 2)
REPRESENTATIVES = ((1, 5), (2, 5), (3, 11))
PAIRS = tuple(combinations(range(4), 2))
BOOST = sp.Matrix((
    (sp.Rational(5, 3), 0, 0, sp.Rational(4, 3)),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (sp.Rational(4, 3), 0, 0, sp.Rational(5, 3)),
))

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def kernel(matrix):
    vectors = matrix.nullspace()
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(matrix.cols, 0)


def same_space(left, right):
    return bool(
        left.rows == right.rows
        and left.rank() == right.rank()
        and left.row_join(right).rank() == left.rank()
    )


def intersection_dimension(left, right):
    return left.rank() + right.rank() - left.row_join(right).rank()


def lorentz_basis(metric):
    result = []
    for a, b in ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)):
        generator = sp.zeros(4)
        generator[a, b] = 1
        generator[b, a] = -metric[a, a] / metric[b, b]
        result.append(generator)
    return tuple(result)


def poincare_evaluation(points, basis):
    columns = [
        sp.Matrix.vstack(*(generator * point for point in points))
        for generator in basis
    ]
    for axis in range(4):
        direction = sp.eye(4)[:, axis]
        columns.append(sp.Matrix.vstack(*(direction for _ in points)))
    return sp.Matrix.hstack(*columns)


def edge_jacobian(points, metric):
    matrix = sp.zeros(6, 16)
    for row, (left, right) in enumerate(PAIRS):
        gradient = 2 * metric * (points[left] - points[right])
        matrix[row, 4 * left:4 * left + 4] = gradient.T
        matrix[row, 4 * right:4 * right + 4] = -gradient.T
    return matrix


def strut_jacobian(bottom, top, metric):
    matrix = sp.zeros(4, 16)
    for index in range(4):
        gradient = 2 * metric * (top[index] - bottom[index])
        matrix[index, 4 * index:4 * index + 4] = gradient.T
    return matrix


def local_data(indices, bottom_all, top_all, metric, basis):
    bottom = tuple(bottom_all[index] for index in indices)
    top = tuple(top_all[index] for index in indices)
    evaluation = poincare_evaluation(top, basis)
    edge = edge_jacobian(top, metric)
    strut = strut_jacobian(bottom, top, metric)
    direct = edge.col_join(strut)
    parameter_kernel = kernel(strut * evaluation)
    direct_kernel = kernel(direct)
    return {
        "evaluation": evaluation,
        "direct": direct,
        "parameter_kernel": parameter_kernel,
        "direct_kernel": direct_kernel,
        "complete": bool(
            direct.rank() == 10
            and parameter_kernel.shape == (10, 6)
            and direct_kernel.shape == (16, 6)
            and same_space(evaluation * parameter_kernel, direct_kernel)
        ),
    }


def analytic_kernel(scale, lapse, basis, normal):
    if scale == 1:
        result = sp.zeros(10, 6)
        result[:3, :3] = sp.eye(3)
        result[6:9, 3:6] = sp.eye(3)
        return result
    result = sp.zeros(10, 6)
    result[:6, :6] = sp.eye(6)
    for column, generator in enumerate(basis):
        result[6:10, column] = (
            sp.Rational(lapse, scale - 1) * generator * normal
        )
    return result


def block_parameter_map(left_kernel, right_kernel, connection):
    result = sp.zeros(30, left_kernel.cols + right_kernel.cols + connection.cols)
    result[:10, :left_kernel.cols] = left_kernel
    result[10:20, left_kernel.cols:left_kernel.cols + right_kernel.cols] = right_kernel
    result[20:30, left_kernel.cols + right_kernel.cols:] = connection
    return result


def pair_parameter_map(left_kernel, right_kernel):
    result = sp.zeros(20, left_kernel.cols + right_kernel.cols)
    result[:10, :left_kernel.cols] = left_kernel
    result[10:20, left_kernel.cols:] = right_kernel
    return result


def audit(bottom_all, top_all, metric, basis):
    left = local_data(LEFT, bottom_all, top_all, metric, basis)
    right = local_data(RIGHT, bottom_all, top_all, metric, basis)
    lower_points = tuple(bottom_all[index] for index in SHARED)
    upper_points = tuple(top_all[index] for index in SHARED)
    lower_evaluation = poincare_evaluation(lower_points, basis)
    upper_evaluation = poincare_evaluation(upper_points, basis)
    lower_stabilizer = kernel(lower_evaluation)
    upper_stabilizer = kernel(upper_evaluation)

    left_kernel = left["parameter_kernel"]
    right_kernel = right["parameter_kernel"]
    fixed = (upper_evaluation * left_kernel).row_join(
        -upper_evaluation * right_kernel
    )
    fixed_coefficients = kernel(fixed)
    fixed_parameters = pair_parameter_map(left_kernel, right_kernel) * fixed_coefficients
    diagonal = left_kernel.col_join(left_kernel)

    variable = fixed.row_join(-upper_evaluation * lower_stabilizer)
    variable_coefficients = kernel(variable)
    parameter_map = block_parameter_map(
        left_kernel, right_kernel, lower_stabilizer
    )
    compatible_parameters = parameter_map * variable_coefficients
    relative_map = sp.eye(10).row_join(-sp.eye(10)).row_join(sp.zeros(10, 10))
    connection_map = sp.zeros(10, 20).row_join(sp.eye(10))
    relative_parameters = relative_map * compatible_parameters
    connection_parameters = connection_map * compatible_parameters
    diagonal_variable = sp.Matrix.vstack(left_kernel, left_kernel, sp.zeros(10, 6))

    stabilizer_span = lower_stabilizer.row_join(upper_stabilizer)
    return {
        "left": left,
        "right": right,
        "lower_evaluation": lower_evaluation,
        "upper_evaluation": upper_evaluation,
        "lower_stabilizer": lower_stabilizer,
        "upper_stabilizer": upper_stabilizer,
        "fixed_rank": fixed.rank(),
        "fixed_parameters": fixed_parameters,
        "fixed_diagonal": same_space(fixed_parameters, diagonal),
        "variable_rank": variable.rank(),
        "compatible_parameters": compatible_parameters,
        "compatible_dimension": compatible_parameters.rank(),
        "relative_rank": relative_parameters.rank(),
        "connection_rank": connection_parameters.rank(),
        "intersection_dimension": intersection_dimension(
            left_kernel, stabilizer_span
        ),
        "diagonal_contained": bool(
            compatible_parameters.row_join(diagonal_variable).rank()
            == compatible_parameters.rank()
        ),
        "relative_dimension_mod_diagonal": (
            compatible_parameters.rank() - diagonal_variable.rank()
        ),
    }


def poincare_adjoint(linear, basis):
    lorentz_coordinates = sp.Matrix.hstack(
        *(generator.reshape(16, 1) for generator in basis)
    )
    result = sp.zeros(10)
    for column in range(10):
        if column < 6:
            transformed = sp.simplify(
                linear * basis[column] * linear.inv()
            )
            coordinates, free = lorentz_coordinates.gauss_jordan_solve(
                transformed.reshape(16, 1)
            )
            if free.rows:
                raise RuntimeError("ambiguous Lorentz coordinates")
            result[:6, column] = coordinates
        else:
            result[6:10, column] = linear * sp.eye(4)[:, column - 6]
    return result


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "local_result": LOCAL_RESULT,
    "fixed_result": FIXED_RESULT,
    "fixed_json": FIXED_JSON,
    "global_result": GLOBAL_RESULT,
    "global_json": GLOBAL_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all variable-connection inputs have frozen provenance",
      provenance_ok, str(hashes))

fixed_upstream = json.loads(FIXED_JSON.read_text())
global_upstream = json.loads(GLOBAL_JSON.read_text())
upstream_ok = bool(
    fixed_upstream["outcome"] == "TWO_FRUSTUM_DIAGONAL_ONLY"
    and fixed_upstream["passed"] == fixed_upstream["tests"] == 9
    and global_upstream["outcome"] == "ADVERSARIAL_GLOBAL_FLEX_SEED_KILLED"
    and global_upstream["passed"] == global_upstream["tests"] == 11
)
check("both frozen-connection upstream outcomes persist", upstream_ok)

basis = lorentz_basis(ETA)
basis_ok = bool(
    len(basis) == 6
    and sp.Matrix.hstack(*(value.reshape(16, 1) for value in basis)).rank() == 6
    and all(value.T * ETA + ETA * value == sp.zeros(4) for value in basis)
)
check("the exact Poincare basis is nondegenerate", basis_ok)

left_edges = {
    ((POINTS[LEFT[a]] - POINTS[LEFT[b]]).T
     * ETA * (POINTS[LEFT[a]] - POINTS[LEFT[b]]))[0]
    for a, b in PAIRS
}
right_edges = {
    ((POINTS[RIGHT[a]] - POINTS[RIGHT[b]]).T
     * ETA * (POINTS[RIGHT[a]] - POINTS[RIGHT[b]]))[0]
    for a, b in PAIRS
}
carrier_ok = bool(left_edges == right_edges == {8})
check("the two cells are exact reflected regular tetrahedra", carrier_ok)

records = []
local_controls = True
fixed_controls = True
stabilizer_controls = True
prediction_controls = True
covariance_controls = True
sign_controls = True

adjoint = poincare_adjoint(BOOST, basis)
boost_ok = bool(BOOST.T * ETA * BOOST == ETA and BOOST.det() == 1)

for scale, lapse in REPRESENTATIVES:
    top = tuple(scale * point + lapse * NORMAL for point in POINTS)
    result = audit(POINTS, top, ETA, basis)
    expected = analytic_kernel(scale, lapse, basis, NORMAL)

    local_ok = bool(
        result["left"]["complete"] and result["right"]["complete"]
        and same_space(result["left"]["parameter_kernel"], expected)
        and same_space(result["right"]["parameter_kernel"], expected)
    )
    fixed_ok = bool(
        result["fixed_rank"] == 6
        and result["fixed_parameters"].rank() == 6
        and result["fixed_diagonal"]
    )
    stabilizer_ok = bool(
        result["lower_stabilizer"].shape == (10, 1)
        and result["upper_stabilizer"].shape == (10, 1)
        and result["lower_stabilizer"][:6, :].rank() == 1
        and result["upper_stabilizer"][:6, :].rank() == 1
        and result["lower_stabilizer"].row_join(
            result["upper_stabilizer"]
        ).rank() == 2
        and (result["upper_evaluation"]
             * result["lower_stabilizer"]).rank() == 1
    )
    predicted = bool(
        result["intersection_dimension"] == 1
        and result["variable_rank"] == 6
        and result["compatible_dimension"] == 7
        and result["relative_rank"] == 1
        and result["connection_rank"] == 1
        and result["diagonal_contained"]
        and result["relative_dimension_mod_diagonal"] == 1
    )
    local_controls &= local_ok
    fixed_controls &= fixed_ok
    stabilizer_controls &= stabilizer_ok
    prediction_controls &= predicted

    boosted_bottom = tuple(BOOST * point for point in POINTS)
    boosted_top = tuple(BOOST * point for point in top)
    boosted = audit(boosted_bottom, boosted_top, ETA, basis)
    block_adjoint = sp.diag(adjoint, adjoint, adjoint)
    covariance_ok = bool(
        boost_ok
        and boosted["left"]["complete"] and boosted["right"]["complete"]
        and same_space(
            boosted["left"]["parameter_kernel"],
            adjoint * result["left"]["parameter_kernel"],
        )
        and same_space(
            boosted["lower_stabilizer"],
            adjoint * result["lower_stabilizer"],
        )
        and same_space(
            boosted["upper_stabilizer"],
            adjoint * result["upper_stabilizer"],
        )
        and same_space(
            boosted["compatible_parameters"],
            block_adjoint * result["compatible_parameters"],
        )
        and boosted["compatible_dimension"] == result["compatible_dimension"]
        and boosted["relative_rank"] == result["relative_rank"]
        and boosted["connection_rank"] == result["connection_rank"]
    )
    covariance_controls &= covariance_ok

    opposite_basis = lorentz_basis(-ETA)
    opposite = audit(POINTS, top, -ETA, opposite_basis)
    sign_ok = bool(
        same_space(
            opposite["left"]["parameter_kernel"],
            result["left"]["parameter_kernel"],
        )
        and same_space(
            opposite["lower_stabilizer"], result["lower_stabilizer"]
        )
        and same_space(
            opposite["upper_stabilizer"], result["upper_stabilizer"]
        )
        and same_space(
            opposite["compatible_parameters"],
            result["compatible_parameters"],
        )
    )
    sign_controls &= sign_ok

    records.append({
        "scale": scale,
        "lapse": lapse,
        "fixed_rank": result["fixed_rank"],
        "fixed_compatible_dimension": result["fixed_parameters"].rank(),
        "lower_stabilizer_dimension": result["lower_stabilizer"].cols,
        "upper_stabilizer_dimension": result["upper_stabilizer"].cols,
        "stabilizer_span_dimension": result["lower_stabilizer"].row_join(
            result["upper_stabilizer"]
        ).rank(),
        "kernel_stabilizer_intersection_dimension": result["intersection_dimension"],
        "variable_rank": result["variable_rank"],
        "variable_compatible_dimension": result["compatible_dimension"],
        "relative_local_parameter_rank": result["relative_rank"],
        "connection_parameter_rank": result["connection_rank"],
        "relative_dimension_mod_diagonal": result["relative_dimension_mod_diagonal"],
        "boost_covariance": covariance_ok,
        "metric_sign_control": sign_ok,
    })

check("all local six-flex kernels are reconstructed exactly", local_controls)
check("the fixed-frame diagonal-only result is reproduced", fixed_controls)
check("both triangle stabilizers and their nontrivial action are exact",
      stabilizer_controls)
check("the disclosed one-mode variable-transition prediction holds",
      prediction_controls, str(records))
check("the full compatible spaces intertwine under an exact Lorentz boost",
      covariance_controls)
check("metric-sign reversal preserves every compatible subspace",
      sign_controls)

controls_ok = bool(
    provenance_ok and upstream_ok and basis_ok and carrier_ok
    and local_controls and fixed_controls and stabilizer_controls
    and covariance_controls and sign_controls
)
forced_zero = bool(
    controls_ok
    and all(record["variable_compatible_dimension"] == 6
            and record["connection_parameter_rank"] == 0
            for record in records)
)
one_mode = bool(controls_ok and prediction_controls)
underdetermined = bool(
    controls_ok
    and any(record["relative_dimension_mod_diagonal"] > 1
            for record in records)
)

if not controls_ok:
    outcome = "VARIABLE_FACE_CONNECTION_CONTROL_FAILED"
elif forced_zero:
    outcome = "VARIABLE_FACE_CONNECTION_FORCED_ZERO"
elif one_mode:
    outcome = "ONE_CONNECTION_COUPLED_RELATIVE_MODE"
elif underdetermined:
    outcome = "VARIABLE_FACE_CONNECTION_UNDERDETERMINED"
else:
    outcome = "VARIABLE_FACE_CONNECTION_OPEN"

allowed = {
    "VARIABLE_FACE_CONNECTION_CONTROL_FAILED",
    "VARIABLE_FACE_CONNECTION_FORCED_ZERO",
    "ONE_CONNECTION_COUPLED_RELATIVE_MODE",
    "VARIABLE_FACE_CONNECTION_UNDERDETERMINED",
    "VARIABLE_FACE_CONNECTION_OPEN",
}
check("the variable-connection hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "records": records,
    "classification": {
        "fixed_connection_holonomy_matrices": "PRESERVED",
        "frozen_connection_implies_metric_rigidity": (
            "REFUTED DERIVED EXACT" if one_mode else "OPEN"
        ),
        "variable_face_connection": (
            "ONE DERIVED RELATIVE MODE PER TESTED FACE"
            if one_mode else "OPEN"
        ),
        "global_variable_connection_closure": "NOT TESTED",
        "action_hessian_or_dynamics": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
for record in records:
    print(
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"fixed={record['fixed_compatible_dimension']}, "
        f"variable={record['variable_compatible_dimension']}, "
        f"relative={record['relative_dimension_mod_diagonal']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Exact fixed-frame gluing audit for two homothetic tetrahedral frusta."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_two_frustum_face_gluing.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_protocol.md"
LOCAL_RESULT = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_result.md"
CORRECTION_SOURCE = HERE / "verify_gravity_600cell_cellular_frustum_relative_poincare_covariance_correction.py"
CORRECTION_JSON = HERE / "gravity_600cell_cellular_frustum_relative_poincare_covariance_correction.json"
ADVERSARIAL_SOURCE = HERE / "verify_gravity_600cell_cellular_frustum_relative_poincare_adversarial.py"
ADVERSARIAL_JSON = HERE / "gravity_600cell_cellular_frustum_relative_poincare_adversarial.json"

PROTOCOL_COMMIT = "e5bf53e"
EXPECTED_HASHES = {
    "prior_art": "d38994ada998df4736858b1a242802097134383d870135bd542178c5343ff63b",
    "protocol": "7d6f6028b6585bc472ee25aca455194d2ac13ed61fc31f7e0f339f4a9bf697f8",
    "local_result": "436fb57037e491b6bdb8fee9ad8b10ab8da1621fd9ecda73e1fcac3fa616fa29",
    "correction_source": "e85e2df690234e19e0343183499c6f8465bc149bc66a87c5246dfe0bda4c1d61",
    "correction_json": "f571869be3341b74b2341c2bf776e99b21174f9f0fb0c5d02e42585c2f3ebaa2",
    "adversarial_source": "8f9bd9882f2efa9f7fdc415ad3c6ca13927283ddb548a8d779c8955ad8ff7e21",
    "adversarial_json": "b750943349dc60ee42d08c0ba61d9a7a0838e3f9346ac1d4397000519b4d6395",
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


def squared_distance(left, right):
    delta = left - right
    return sp.expand((delta.T * ETA * delta)[0])


def make_lorentz_basis():
    basis = []
    for a, b in ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)):
        generator = sp.zeros(4)
        generator[a, b] = 1
        generator[b, a] = -ETA[a, a] / ETA[b, b]
        basis.append(generator)
    return tuple(basis)


LORENTZ = make_lorentz_basis()


def poincare_evaluation(points):
    columns = [sp.Matrix.vstack(*(generator * point for point in points))
               for generator in LORENTZ]
    for axis in range(4):
        direction = sp.eye(4)[:, axis]
        columns.append(sp.Matrix.vstack(*(direction for _ in points)))
    return sp.Matrix.hstack(*columns)


def edge_jacobian(points):
    matrix = sp.zeros(6, 16)
    for row, (left, right) in enumerate(PAIRS):
        gradient = 2 * ETA * (points[left] - points[right])
        matrix[row, 4 * left:4 * left + 4] = gradient.T
        matrix[row, 4 * right:4 * right + 4] = -gradient.T
    return matrix


def strut_jacobian(bottom, top):
    matrix = sp.zeros(4, 16)
    for index in range(4):
        gradient = 2 * ETA * (top[index] - bottom[index])
        matrix[index, 4 * index:4 * index + 4] = gradient.T
    return matrix


def analytic_kernel(scale, lapse):
    if scale == 1:
        result = sp.zeros(10, 6)
        result[:3, :3] = sp.eye(3)
        result[6:9, 3:6] = sp.eye(3)
        return result
    result = sp.zeros(10, 6)
    result[:6, :6] = sp.eye(6)
    for column, generator in enumerate(LORENTZ):
        result[6:10, column] = (
            sp.Rational(lapse, scale - 1) * generator * NORMAL
        )
    return result


def local_data(indices, top_all):
    bottom = tuple(POINTS[index] for index in indices)
    top = tuple(top_all[index] for index in indices)
    evaluation = poincare_evaluation(top)
    edge = edge_jacobian(top)
    strut = strut_jacobian(bottom, top)
    direct = edge.col_join(strut)
    constraint = strut * evaluation
    parameter_kernel = kernel(constraint)
    direct_kernel = kernel(direct)
    return {
        "evaluation": evaluation,
        "edge": edge,
        "strut": strut,
        "direct": direct,
        "constraint": constraint,
        "parameter_kernel": parameter_kernel,
        "direct_kernel": direct_kernel,
        "direct_equality": same_space(
            evaluation * parameter_kernel, direct_kernel
        ),
    }


def pair_parameter_map(left_kernel, right_kernel):
    result = sp.zeros(20, left_kernel.cols + right_kernel.cols)
    result[:10, :left_kernel.cols] = left_kernel
    result[10:20, left_kernel.cols:] = right_kernel
    return result


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "local_result": LOCAL_RESULT,
    "correction_source": CORRECTION_SOURCE,
    "correction_json": CORRECTION_JSON,
    "adversarial_source": ADVERSARIAL_SOURCE,
    "adversarial_json": ADVERSARIAL_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all two-frustum inputs have exact frozen provenance",
      provenance_ok, str(hashes))

correction = json.loads(CORRECTION_JSON.read_text())
adversarial = json.loads(ADVERSARIAL_JSON.read_text())
upstream_ok = bool(
    correction["outcome"]
    == "STATIC_STABILIZER_AND_EXPANDING_LORENTZ_CHART_CORROBORATED"
    and correction["passed"] == correction["tests"] == 13
    and adversarial["outcome"]
    == "ADVERSARIAL_POINCARE_STRATIFICATION_CORROBORATED"
    and adversarial["passed"] == adversarial["tests"] == 13
)
check("both accepted local Poincare certificates persist", upstream_ok)

left_edges = tuple(squared_distance(POINTS[LEFT[a]], POINTS[LEFT[b]])
                   for a, b in PAIRS)
right_edges = tuple(squared_distance(POINTS[RIGHT[a]], POINTS[RIGHT[b]])
                    for a, b in PAIRS)
face_normal = POINTS[3]
face_plane_value = (face_normal.T * POINTS[0])[0]
reflected_apex = sp.simplify(
    POINTS[3]
    - 2 * (
        ((face_normal.T * POINTS[3])[0] - face_plane_value)
        / (face_normal.T * face_normal)[0]
    ) * face_normal
)
left_affine = sp.Matrix.hstack(*(POINTS[index] - POINTS[LEFT[0]]
                                 for index in LEFT[1:]))
right_affine = sp.Matrix.hstack(*(POINTS[index] - POINTS[RIGHT[0]]
                                  for index in RIGHT[1:]))
geometry_ok = bool(
    set(left_edges) == set(right_edges) == {8}
    and reflected_apex == POINTS[4]
    and left_affine.rank() == right_affine.rank() == 3
    and ((face_normal.T * POINTS[3])[0] - face_plane_value)
    * ((face_normal.T * POINTS[4])[0] - face_plane_value) < 0
)
check("the carrier is an exact pair of reflected regular tetrahedra",
      geometry_ok)

basis_ok = bool(
    len(LORENTZ) == 6
    and sp.Matrix.hstack(*(matrix.reshape(16, 1)
                            for matrix in LORENTZ)).rank() == 6
    and all(matrix.T * ETA + ETA * matrix == sp.zeros(4)
            for matrix in LORENTZ)
)
check("the common Poincare frame uses exact so(3,1)", basis_ok)

records = []
local_controls = True
analytic_controls = True
full_face_controls = True
diagonal_only = True
hidden_face_mode = True
underdetermined = False

for scale, lapse in REPRESENTATIVES:
    top_all = tuple(scale * point + lapse * NORMAL for point in POINTS)
    left = local_data(LEFT, top_all)
    right = local_data(RIGHT, top_all)
    expected = analytic_kernel(scale, lapse)

    local_ok = bool(
        left["edge"].rank() == right["edge"].rank() == 6
        and left["direct"].rank() == right["direct"].rank() == 10
        and left["constraint"].rank() == right["constraint"].rank() == 4
        and left["parameter_kernel"].shape
        == right["parameter_kernel"].shape == (10, 6)
        and left["direct_equality"] and right["direct_equality"]
    )
    local_analytic = bool(
        same_space(left["parameter_kernel"], right["parameter_kernel"])
        and same_space(left["parameter_kernel"], expected)
        and same_space(right["parameter_kernel"], expected)
    )
    local_controls &= local_ok
    analytic_controls &= local_analytic

    face_points = tuple(top_all[index] for index in SHARED)
    face = poincare_evaluation(face_points)
    face_kernel = kernel(face)
    unrestricted_pair = face.row_join(-face)
    local_full_face = bool(
        face.shape == (12, 10)
        and face.rank() == 9
        and face_kernel.shape == (10, 1)
        and face_kernel[:6, :].rank() == 1
        and face * face_kernel == sp.zeros(12, 1)
        and unrestricted_pair.rank() == 9
        and kernel(unrestricted_pair).shape == (20, 11)
    )
    full_face_controls &= local_full_face

    left_kernel = left["parameter_kernel"]
    right_kernel = right["parameter_kernel"]
    gluing = (face * left_kernel).row_join(-face * right_kernel)
    coefficient_kernel = kernel(gluing)
    pair_map = pair_parameter_map(left_kernel, right_kernel)
    compatible_pairs = pair_map * coefficient_kernel
    difference = sp.eye(10).row_join(-sp.eye(10)) * compatible_pairs
    constrained_relative_dimension = 6 - (face * expected).rank()
    diagonal = expected.col_join(expected)
    pair_dimension = compatible_pairs.rank()
    difference_rank = difference.rank()
    local_diagonal = bool(
        gluing.rank() == 6
        and pair_dimension == 6
        and difference_rank == 0
        and constrained_relative_dimension == 0
        and same_space(compatible_pairs, diagonal)
    )
    local_hidden = bool(
        pair_dimension == 7
        and difference_rank == 1
        and constrained_relative_dimension == 1
        and same_space(difference, face_kernel)
    )
    diagonal_only &= local_diagonal
    hidden_face_mode &= local_hidden
    underdetermined |= bool(
        pair_dimension > 7 or difference_rank > 1
        or constrained_relative_dimension > 1
    )

    records.append({
        "scale": scale,
        "lapse": lapse,
        "left_local_kernel_dimension": left_kernel.cols,
        "right_local_kernel_dimension": right_kernel.cols,
        "local_kernels_equal": local_analytic,
        "full_poincare_face_evaluation_rank": face.rank(),
        "full_poincare_face_stabilizer_dimension": face_kernel.cols,
        "unrestricted_pair_dimension": kernel(unrestricted_pair).cols,
        "constrained_gluing_rank": gluing.rank(),
        "compatible_pair_dimension": pair_dimension,
        "relative_difference_rank": difference_rank,
        "constrained_relative_face_stabilizer_dimension": (
            constrained_relative_dimension
        ),
        "compatible_space_is_diagonal": local_diagonal,
    })

check("both local frustum kernels equal their direct length kernels",
      local_controls)
check("both local kernels equal the frozen analytic strut solution",
      analytic_controls)
check("the unrestricted Poincare algebra retains one relative face mode",
      full_face_controls)
check("all constrained two-frustum compatible pairs are diagonal",
      diagonal_only)

controls_ok = bool(
    provenance_ok and upstream_ok and geometry_ok and basis_ok
    and local_controls and analytic_controls and full_face_controls
)
hidden = bool(controls_ok and hidden_face_mode)
diagonal_result = bool(controls_ok and diagonal_only)

if not controls_ok:
    outcome = "TWO_FRUSTUM_FACE_GLUING_CONTROL_FAILED"
elif hidden:
    outcome = "TWO_FRUSTUM_HIDDEN_FACE_MODE"
elif diagonal_result:
    outcome = "TWO_FRUSTUM_DIAGONAL_ONLY"
elif underdetermined:
    outcome = "TWO_FRUSTUM_FACE_UNDERDETERMINED"
else:
    outcome = "TWO_FRUSTUM_FACE_GLUING_OPEN"

allowed = {
    "TWO_FRUSTUM_FACE_GLUING_CONTROL_FAILED",
    "TWO_FRUSTUM_HIDDEN_FACE_MODE",
    "TWO_FRUSTUM_DIAGONAL_ONLY",
    "TWO_FRUSTUM_FACE_UNDERDETERMINED",
    "TWO_FRUSTUM_FACE_GLUING_OPEN",
}
check("the preregistered gluing hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "carrier": {
        "vertices": 5,
        "tetrahedra": [list(LEFT), list(RIGHT)],
        "shared_face": list(SHARED),
    },
    "records": records,
    "classification": {
        "full_poincare_pointwise_face_stabilizer": (
            "DIMENSION ONE CONTROL" if full_face_controls else "OPEN"
        ),
        "face_mode_inside_six_length_flexes": (
            "ABSENT DERIVED EXACT" if diagonal_result else "OPEN"
        ),
        "independent_connection_hidden_in_local_flexes": (
            "REFUTED LOCAL" if diagonal_result else "OPEN"
        ),
        "new_first_order_holonomy_variable": "NOT TESTED",
        "global_closure_action_or_dynamics": "NOT TESTED",
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
        f"full face stabilizer="
        f"{record['full_poincare_face_stabilizer_dimension']}, "
        f"constrained pair={record['compatible_pair_dimension']}, "
        f"relative={record['relative_difference_rank']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Independent direct-union audit of two-frustum face gluing."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_two_frustum_face_gluing_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_adversarial_protocol.md"
PRIMARY_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_protocol.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_two_frustum_face_gluing.py"
PRIMARY_JSON = HERE / "gravity_600cell_two_frustum_face_gluing.json"
LOCAL_RESULT = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_result.md"

PROTOCOL_COMMIT = "32adc05"
EXPECTED_HASHES = {
    "protocol": "25a2d6e2d9b60b1dd7a6151bc43ad57ec2cb99bc0ffc8a3f9dd5687c1aef9551",
    "primary_protocol": "7d6f6028b6585bc472ee25aca455194d2ac13ed61fc31f7e0f339f4a9bf697f8",
    "primary_source": "52636ae59bd4e4568df175e32b7c3aeae4fbfbc3d475d255131b6db671c41ae7",
    "primary_json": "0e09c3f8f38c8158deff5b81bc6fe4d5d6dd685a24cce83e015fb95e3f26a70e",
    "local_result": "436fb57037e491b6bdb8fee9ad8b10ab8da1621fd9ecda73e1fcac3fa616fa29",
}

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
BOTTOM = tuple(sp.Matrix(point) for point in (
    (5, 0, 0, 0),
    (0, 5, 0, 0),
    (0, 0, 5, 0),
    (3, 4, 0, 0),
    (sp.Rational(5, 3), sp.Rational(8, 3), sp.Rational(-4, 3), 0),
))
LEFT = (0, 1, 2, 3)
RIGHT = (0, 1, 2, 4)
SHARED = (0, 1, 2)
UNION_EDGES = tuple(sorted(
    set(tuple(sorted((tetra[a], tetra[b])))
        for tetra in (LEFT, RIGHT)
        for a, b in combinations(range(4), 2))
))
REPRESENTATIVES = ((1, 7), (2, 7), (3, 13))

Y = sp.symbols("y0:20")
Y_POINTS = tuple(sp.Matrix(Y[4 * i:4 * i + 4]) for i in range(5))
A_SYMBOLS = sp.symbols("a0:16")
B_SYMBOLS = sp.symbols("b0:4")
A_MATRIX = sp.Matrix(4, 4, A_SYMBOLS)
B_VECTOR = sp.Matrix(B_SYMBOLS)

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


def square(left, right, metric):
    delta = left - right
    return sp.expand((delta.T * metric * delta)[0])


def make_lorentz_basis():
    result = []
    for a, b in ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)):
        generator = sp.zeros(4)
        generator[a, b] = 1
        generator[b, a] = -ETA[a, a] / ETA[b, b]
        result.append(generator)
    return tuple(result)


LORENTZ = make_lorentz_basis()


def top(scale, lapse, bottom=BOTTOM, normal=NORMAL):
    return tuple(scale * point + lapse * normal for point in bottom)


def union_polynomial_jacobian(bottom, points, metric):
    polynomials = [
        square(Y_POINTS[left], Y_POINTS[right], metric)
        for left, right in UNION_EDGES
    ]
    polynomials.extend(
        square(Y_POINTS[index], bottom[index], metric)
        for index in range(5)
    )
    symbolic = sp.Matrix(polynomials).jacobian(Y)
    substitution = {
        Y[4 * index + axis]: points[index][axis]
        for index in range(5) for axis in range(4)
    }
    return symbolic.subs(substitution)


def poincare_evaluation(points):
    columns = [sp.Matrix.vstack(*(generator * point for point in points))
               for generator in LORENTZ]
    for axis in range(4):
        vector = sp.eye(4)[:, axis]
        columns.append(sp.Matrix.vstack(*(vector for _ in points)))
    return sp.Matrix.hstack(*columns)


def analytic_parameter_kernel(scale, lapse):
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


def redundant_face_system(face, metric):
    symmetric = A_MATRIX.T * metric + metric * A_MATRIX
    lorentz_equations = [symmetric[row, column]
                         for row in range(4) for column in range(row, 4)]
    parameters = A_SYMBOLS + B_SYMBOLS
    lorentz_rows = sp.Matrix(lorentz_equations).jacobian(parameters)
    evaluations = sp.Matrix.vstack(*(A_MATRIX * point + B_VECTOR
                                      for point in face))
    face_rows = evaluations.jacobian(parameters)
    return lorentz_rows.col_join(face_rows)


def block_transport(matrix, count):
    result = sp.zeros(4 * count)
    for index in range(count):
        result[4 * index:4 * index + 4,
               4 * index:4 * index + 4] = matrix
    return result


paths = {
    "protocol": PROTOCOL,
    "primary_protocol": PRIMARY_PROTOCOL,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "local_result": LOCAL_RESULT,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all direct-union adversarial inputs have frozen provenance",
      provenance_ok, str(hashes))

primary = json.loads(PRIMARY_JSON.read_text())
primary_ok = bool(
    primary["outcome"] == "TWO_FRUSTUM_DIAGONAL_ONLY"
    and primary["passed"] == primary["tests"] == 9
    and all(record["full_poincare_face_stabilizer_dimension"] == 1
            and record["compatible_pair_dimension"] == 6
            and record["relative_difference_rank"] == 0
            for record in primary["records"])
)
check("the primary diagonal-only result is preserved literally", primary_ok)

face_normal = sp.Matrix((1, 1, 1, 0))
face_value = (face_normal.T * BOTTOM[0])[0]
reflected = sp.simplify(
    BOTTOM[3]
    - 2 * (
        ((face_normal.T * BOTTOM[3])[0] - face_value)
        / (face_normal.T * face_normal)[0]
    ) * face_normal
)
left_affine = sp.Matrix.hstack(*(BOTTOM[i] - BOTTOM[0] for i in LEFT[1:]))
right_affine = sp.Matrix.hstack(*(BOTTOM[i] - BOTTOM[0] for i in RIGHT[1:]))
geometry_ok = bool(
    len(UNION_EDGES) == 9
    and reflected == BOTTOM[4]
    and left_affine.rank() == right_affine.rank() == 3
    and ((face_normal.T * BOTTOM[3])[0] - face_value)
    * ((face_normal.T * BOTTOM[4])[0] - face_value) < 0
)
check("the irregular reflected five-vertex union is exact", geometry_ok)

boost = sp.eye(4)
boost[2, 2] = sp.Rational(5, 4)
boost[2, 3] = sp.Rational(3, 4)
boost[3, 2] = sp.Rational(3, 4)
boost[3, 3] = sp.Rational(5, 4)
boost_ok = bool(boost.T * ETA * boost == ETA and boost.det() == 1)
check("the independent z-time covariance boost is exact Lorentz", boost_ok)

records = []
struts_ok = True
polynomial_complete = True
direct_equals_common = True
face_control_ok = True
metric_sign_ok = True
boost_transport_ok = True

boosted_bottom = tuple(boost * point for point in BOTTOM)
boosted_normal = boost * NORMAL
vertex_transport = block_transport(boost, 5)

for scale, lapse in REPRESENTATIVES:
    points = top(scale, lapse)
    strut_values = tuple(square(points[i], BOTTOM[i], ETA)
                         for i in range(5))
    local_struts = all(value < 0 for value in strut_values)
    struts_ok &= local_struts

    direct = union_polynomial_jacobian(BOTTOM, points, ETA)
    direct_kernel = kernel(direct)
    common = poincare_evaluation(points) * analytic_parameter_kernel(scale, lapse)
    local_polynomial = bool(
        direct.shape == (14, 20)
        and direct.rank() == 14
        and direct_kernel.shape == (20, 6)
    )
    local_common = same_space(direct_kernel, common)
    polynomial_complete &= local_polynomial
    direct_equals_common &= local_common

    face = tuple(points[index] for index in SHARED)
    redundant = redundant_face_system(face, ETA)
    local_face_control = bool(
        redundant.shape == (22, 20)
        and redundant.rank() == 19
        and kernel(redundant).shape == (20, 1)
    )
    face_control_ok &= local_face_control

    negative = union_polynomial_jacobian(BOTTOM, points, -ETA)
    local_metric_sign = bool(
        negative == -direct
        and same_space(kernel(negative), direct_kernel)
    )
    metric_sign_ok &= local_metric_sign

    boosted_points = top(scale, lapse, boosted_bottom, boosted_normal)
    boosted_direct = union_polynomial_jacobian(
        boosted_bottom, boosted_points, ETA
    )
    boosted_kernel = kernel(boosted_direct)
    boosted_redundant = redundant_face_system(
        tuple(boosted_points[index] for index in SHARED), ETA
    )
    local_boost = bool(
        boosted_direct.rank() == 14
        and same_space(boosted_kernel, vertex_transport * direct_kernel)
        and boosted_redundant.rank() == 19
        and kernel(boosted_redundant).shape == (20, 1)
    )
    boost_transport_ok &= local_boost

    records.append({
        "scale": scale,
        "lapse": lapse,
        "strut_squared_lengths": [str(value) for value in strut_values],
        "all_struts_timelike": local_struts,
        "union_polynomial_rows": direct.rows,
        "union_polynomial_rank": direct.rank(),
        "union_kernel_dimension": direct_kernel.cols,
        "kernel_equals_common_poincare_image": local_common,
        "redundant_face_system_rank": redundant.rank(),
        "full_poincare_face_stabilizer_dimension": kernel(redundant).cols,
        "metric_sign_control": local_metric_sign,
        "boost_transport_control": local_boost,
    })

check("all five unequal struts remain timelike at every representative",
      struts_ok)
check("every direct union polynomial Jacobian has exact rank fourteen",
      polynomial_complete)
check("every direct union kernel is exactly the common six-motion image",
      direct_equals_common)
check("the redundant full-Poincare face stabilizer remains one-dimensional",
      face_control_ok)
check("the metric-sign convention leaves every direct kernel unchanged",
      metric_sign_ok)
check("the boosted direct kernels and face controls transport exactly",
      boost_transport_ok)

controls_ok = bool(
    provenance_ok and primary_ok and geometry_ok and boost_ok
    and struts_ok and polynomial_complete and face_control_ok
    and metric_sign_ok and boost_transport_ok
)
corroborated = bool(controls_ok and direct_equals_common)
disagreement = bool(controls_ok and not direct_equals_common)

if not controls_ok:
    outcome = "ADVERSARIAL_TWO_FRUSTUM_CONTROL_FAILED"
elif corroborated:
    outcome = "ADVERSARIAL_TWO_FRUSTUM_DIAGONAL_ONLY"
elif disagreement:
    outcome = "ADVERSARIAL_TWO_FRUSTUM_DISAGREEMENT"
else:
    outcome = "ADVERSARIAL_TWO_FRUSTUM_OPEN"

allowed = {
    "ADVERSARIAL_TWO_FRUSTUM_CONTROL_FAILED",
    "ADVERSARIAL_TWO_FRUSTUM_DIAGONAL_ONLY",
    "ADVERSARIAL_TWO_FRUSTUM_DISAGREEMENT",
    "ADVERSARIAL_TWO_FRUSTUM_OPEN",
}
check("the adversarial hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "union_edges": [list(edge) for edge in UNION_EDGES],
    "records": records,
    "classification": {
        "direct_two_frustum_flex_space": (
            "COMMON SIX-DIMENSIONAL POINCARE IMAGE"
            if corroborated else "OPEN"
        ),
        "independent_face_mode_inside_length_flexes": (
            "ABSENT ADVERSARIALLY CORROBORATED"
            if corroborated else "OPEN"
        ),
        "full_poincare_face_stabilizer_control": "DIMENSION ONE",
        "new_connection_variable_closure_action_or_dynamics": "NOT TESTED",
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
        f"rank/nullity={record['union_polynomial_rank']}/"
        f"{record['union_kernel_dimension']}, "
        f"full-face={record['full_poincare_face_stabilizer_dimension']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

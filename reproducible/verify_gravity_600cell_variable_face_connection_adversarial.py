#!/usr/bin/env python3
"""Independent polynomial audit of variable face-transition gluing."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_variable_face_connection_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_variable_face_connection_adversarial_protocol.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_variable_face_connection_prior_art.md"
PRIMARY_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_variable_face_connection_protocol.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_variable_face_connection.py"
PRIMARY_JSON = HERE / "gravity_600cell_variable_face_connection.json"
FIXED_RESULT = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_result.md"
FIXED_ADV_JSON = HERE / "gravity_600cell_two_frustum_face_gluing_adversarial.json"

PROTOCOL_COMMIT = "3867f97"
EXPECTED_HASHES = {
    "protocol": "94c86b5ed80b4e3bac23b139cb3d9c9f7f708d59ea2a240c21908db8cbf5e34f",
    "prior_art": "2ed809fedad24fa15977b39e4dd6fec386e9080c123208d54fd089554ce44d2d",
    "primary_protocol": "f6b91206a857cda6ebfe5cb9988110de5f12a9c1ca51bcbdb733a8429682ca6a",
    "primary_source": "69a5d7479a5df427cead76f82db31fe62a9190c28c967f699c846881634fb0f6",
    "primary_json": "001212016553d006862e68edc4f780f37ca1476110b6e0aed3e987f52a43b5e3",
    "fixed_result": "b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3",
    "fixed_adv_json": "0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542",
}

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
LEFT = (0, 1, 2, 3)
RIGHT = (0, 1, 2, 4)
SHARED = (0, 1, 2)
PAIRS = tuple(combinations(range(4), 2))
REPRESENTATIVES = ((1, 7), (2, 7), (4, 13))
CARRIERS = {
    "A": tuple(sp.Matrix(point) for point in (
        (0, 0, 0, 0),
        (2, 0, 0, 0),
        (1, 3, 0, 0),
        (1, -1, 4, 0),
        (1, -1, -4, 0),
    )),
    "B": tuple(sp.Matrix(point) for point in (
        (-1, 1, 0, 0),
        (3, 1, 0, 0),
        (0, 4, 0, 0),
        (2, -2, 7, 0),
        (2, -2, -7, 0),
    )),
}

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


def squared_distance(left, right, metric):
    difference = left - right
    return sp.expand((difference.T * metric * difference)[0])


def affine_volume(points, indices):
    base = points[indices[0]][:3, :]
    columns = [points[index][:3, :] - base for index in indices[1:]]
    return sp.Matrix.hstack(*columns).det()


def carrier_control(points):
    left_lengths = {
        squared_distance(points[LEFT[a]], points[LEFT[b]], ETA)
        for a, b in PAIRS
    }
    right_lengths = {
        squared_distance(points[RIGHT[a]], points[RIGHT[b]], ETA)
        for a, b in PAIRS
    }
    return bool(
        affine_volume(points, LEFT) != 0
        and affine_volume(points, RIGHT) != 0
        and len(left_lengths) > 1
        and left_lengths == right_lengths
        and all(points[4][axis] == points[3][axis] for axis in (0, 1, 3))
        and points[4][2] == -points[3][2]
    )


def polynomial_audit(points, scale, lapse, metric):
    left_variables = sp.symbols("left_0:16")
    right_variables = sp.symbols("right_0:16")
    matrix_variables = sp.symbols("matrix_0:16")
    translation_variables = sp.symbols("translation_0:4")
    variables = (
        list(left_variables) + list(right_variables)
        + list(matrix_variables) + list(translation_variables)
    )
    connection_variables = list(matrix_variables) + list(translation_variables)

    left = tuple(sp.Matrix(left_variables[4 * index:4 * index + 4])
                 for index in range(4))
    right = tuple(sp.Matrix(right_variables[4 * index:4 * index + 4])
                  for index in range(4))
    matrix = sp.Matrix(4, 4, matrix_variables)
    translation = sp.Matrix(translation_variables)
    top_all = tuple(scale * point + lapse * NORMAL for point in points)
    left_background = tuple(top_all[index] for index in LEFT)
    right_background = tuple(top_all[index] for index in RIGHT)

    equations = []
    for coordinates, background_indices in ((left, LEFT), (right, RIGHT)):
        for first, second in PAIRS:
            equations.append(squared_distance(
                coordinates[first], coordinates[second], metric
            ))
        for local_index, global_index in enumerate(background_indices):
            equations.append(squared_distance(
                coordinates[local_index], points[global_index], metric
            ))

    lorentz = matrix.T * metric + metric * matrix
    lorentz_equations = [lorentz[row, column]
                         for row in range(4)
                         for column in range(row, 4)]
    equations.extend(lorentz_equations)

    lower_fixing = []
    for global_index in SHARED:
        lower_fixing.extend(matrix * points[global_index] + translation)
    equations.extend(lower_fixing)

    matching = []
    for local_index in range(3):
        matching.extend(
            left[local_index]
            - ((sp.eye(4) + matrix) * right[local_index] + translation)
        )
    equations.extend(matching)

    substitutions = {}
    for coordinates, background in ((left, left_background),
                                     (right, right_background)):
        for point, value in zip(coordinates, background):
            substitutions.update(zip(point, value))
    substitutions.update({value: 0 for value in connection_variables})

    jacobian = sp.Matrix(equations).jacobian(variables).subs(substitutions)
    local_left = jacobian[:10, :16]
    local_right = jacobian[10:20, 16:32]
    connection_system = jacobian[20:42, 32:52]

    upper_fixing = []
    for global_index in SHARED:
        upper_fixing.extend(matrix * top_all[global_index] + translation)
    upper_fixing_jacobian = sp.Matrix(upper_fixing).jacobian(
        connection_variables
    ).subs(substitutions)
    fixed_face_connection = connection_system.col_join(upper_fixing_jacobian)

    frozen_rows = sp.zeros(20, 52)
    frozen_rows[:, 32:52] = sp.eye(20)
    frozen_system = jacobian.col_join(frozen_rows)

    full_kernel_vectors = jacobian.nullspace()
    full_kernel = (
        sp.Matrix.hstack(*full_kernel_vectors)
        if full_kernel_vectors else sp.zeros(52, 0)
    )
    shared_difference = sp.zeros(12, 52)
    for index in range(12):
        shared_difference[index, index] = 1
        shared_difference[index, 16 + index] = -1

    return {
        "local_left_rank": local_left.rank(),
        "local_right_rank": local_right.rank(),
        "connection_rank": connection_system.rank(),
        "connection_nullity": 20 - connection_system.rank(),
        "fixed_face_connection_rank": fixed_face_connection.rank(),
        "fixed_face_connection_nullity": 20 - fixed_face_connection.rank(),
        "full_rank": jacobian.rank(),
        "full_nullity": 52 - jacobian.rank(),
        "frozen_rank": frozen_system.rank(),
        "frozen_nullity": 52 - frozen_system.rank(),
        "kernel_connection_projection_rank": full_kernel[32:52, :].rank(),
        "kernel_shared_difference_rank": (shared_difference * full_kernel).rank(),
    }


paths = {
    "protocol": PROTOCOL,
    "prior_art": PRIOR_ART,
    "primary_protocol": PRIMARY_PROTOCOL,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "fixed_result": FIXED_RESULT,
    "fixed_adv_json": FIXED_ADV_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all adversarial variable-face inputs have frozen provenance",
      provenance_ok, str(hashes))

primary = json.loads(PRIMARY_JSON.read_text())
fixed = json.loads(FIXED_ADV_JSON.read_text())
upstream_ok = bool(
    primary["outcome"] == "ONE_CONNECTION_COUPLED_RELATIVE_MODE"
    and primary["passed"] == primary["tests"] == 11
    and fixed["outcome"] == "ADVERSARIAL_TWO_FRUSTUM_DIAGONAL_ONLY"
    and fixed["passed"] == fixed["tests"] == 11
)
check("the primary correction and old fixed-frame control persist", upstream_ok)

geometry_ok = all(carrier_control(points) for points in CARRIERS.values())
check("both rational carriers are nondegenerate irregular reflected pairs",
      geometry_ok)

records = []
local_controls = True
stabilizer_controls = True
frozen_controls = True
decisive_controls = True
sign_controls = True

for carrier_name, points in CARRIERS.items():
    for scale, lapse in REPRESENTATIVES:
        result = polynomial_audit(points, scale, lapse, ETA)
        opposite = polynomial_audit(points, scale, lapse, -ETA)
        local_ok = bool(
            result["local_left_rank"] == result["local_right_rank"] == 10
        )
        stabilizer_ok = bool(
            result["connection_rank"] == 19
            and result["connection_nullity"] == 1
            and result["fixed_face_connection_rank"] == 20
            and result["fixed_face_connection_nullity"] == 0
        )
        frozen_ok = bool(
            result["frozen_rank"] == 46
            and result["frozen_nullity"] == 6
        )
        decisive = bool(
            result["full_rank"] == 45
            and result["full_nullity"] == 7
            and result["kernel_connection_projection_rank"] == 1
            and result["kernel_shared_difference_rank"] == 1
        )
        sign_ok = result == opposite
        local_controls &= local_ok
        stabilizer_controls &= stabilizer_ok
        frozen_controls &= frozen_ok
        decisive_controls &= decisive
        sign_controls &= sign_ok
        records.append({
            "carrier": carrier_name,
            "scale": scale,
            "lapse": lapse,
            **result,
            "metric_sign_control": sign_ok,
        })

check("all twelve direct local polynomial Jacobians have rank/nullity 10/6",
      local_controls)
check("the lower-face mode is 1D and fixing the upper face kills it",
      stabilizer_controls)
check("freezing all connection entries restores exact nullity six",
      frozen_controls)
check("all six variable systems have nullity seven and rank-one projections",
      decisive_controls, str(records))
check("metric-sign reversal preserves every adversarial rank", sign_controls)

controls_ok = bool(
    provenance_ok and upstream_ok and geometry_ok and local_controls
    and stabilizer_controls and frozen_controls and sign_controls
)
one_mode = bool(controls_ok and decisive_controls)
forced_zero = bool(
    controls_ok
    and all(record["full_nullity"] == 6
            and record["kernel_connection_projection_rank"] == 0
            for record in records)
)
underdetermined = bool(
    controls_ok
    and any(record["full_nullity"] > 7
            or record["kernel_connection_projection_rank"] > 1
            or record["kernel_shared_difference_rank"] > 1
            for record in records)
)

if not controls_ok:
    outcome = "ADVERSARIAL_VARIABLE_FACE_CONTROL_FAILED"
elif one_mode:
    outcome = "ADVERSARIAL_ONE_CONNECTION_MODE"
elif forced_zero:
    outcome = "ADVERSARIAL_CONNECTION_FORCED_ZERO"
elif underdetermined:
    outcome = "ADVERSARIAL_VARIABLE_FACE_UNDERDETERMINED"
else:
    outcome = "ADVERSARIAL_VARIABLE_FACE_OPEN"

allowed = {
    "ADVERSARIAL_VARIABLE_FACE_CONTROL_FAILED",
    "ADVERSARIAL_ONE_CONNECTION_MODE",
    "ADVERSARIAL_CONNECTION_FORCED_ZERO",
    "ADVERSARIAL_VARIABLE_FACE_UNDERDETERMINED",
    "ADVERSARIAL_VARIABLE_FACE_OPEN",
}
check("the adversarial variable-face hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "records": records,
    "classification": {
        "one_connection_coupled_relative_mode": (
            "ADVERSARIALLY CORROBORATED DERIVED EXACT" if one_mode else "OPEN"
        ),
        "fixed_connection_holonomy_implies_metric_rigidity": (
            "REFUTED" if one_mode else "OPEN"
        ),
        "complete_variable_connection_closure": "NOT TESTED",
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
        f"carrier={record['carrier']} "
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"frozen={record['frozen_nullity']}, variable={record['full_nullity']}, "
        f"connection={record['kernel_connection_projection_rank']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

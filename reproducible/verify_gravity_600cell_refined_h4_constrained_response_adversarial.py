#!/usr/bin/env python3
"""Direct-scalar-action adversarial audit of the constrained H4 response."""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
CURVATURE = HERE / "gravity_600cell_refined_local_curvature_mass.json"
FEASIBILITY = HERE / "gravity_600cell_refined_canonical_map_feasibility.json"
NULL_COUPLING = HERE / "gravity_600cell_refined_h4_null_coupling_adversarial.json"
PRIMARY = HERE / "gravity_600cell_refined_h4_constrained_response_corrected.json"
PRIMARY_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_primary_result.md"
)
PRIOR_ART = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_hessian_prior_art.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_adversarial_protocol.md"
)
CELL600 = ROOT / "commons/cell600.py"
OUTPUT = HERE / "gravity_600cell_refined_h4_constrained_response_adversarial.json"

PROTOCOL_COMMIT = "a9ee74b"
EXPECTED_HASHES = {
    "action_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "curvature": "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "feasibility": "ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e",
    "null_coupling": "5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9",
    "primary": "85adea23f6a19153f61f3ed066137a5e40ab77b8901d4cc81cfc4f864e0bc093",
    "primary_result": "fd07977bdb2e45bf3170d1ba98919690e57bc8ee476ce48506859bbffc0253ad",
    "prior_art": "222f31862e911e03a1a7740696618948e370e43164812120120d85e834f0f639",
    "protocol": "b0e4df01e5af9d23418d6dde89dc94cb9722f34594bd1e3a298b1b871596a856",
    "cell600": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}

PAIR4 = tuple(combinations(range(4), 2))
BOUNDARY_VARIABLES = (
    tuple(("old",) + pair for pair in PAIR4)
    + tuple(("new",) + pair for pair in PAIR4)
)
INTERNAL_VARIABLES = (
    tuple(("cross",) + pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
VARIABLES = BOUNDARY_VARIABLES + INTERNAL_VARIABLES
LOCAL_TRIANGLES = np.asarray(tuple(combinations(range(5), 3)), dtype=np.int8)
TAU_TEXT = "0.0102"
ACTION_DPS = 180
REPEAT_DPS = 220
STEP_TEXTS = ("1e-10", "5e-11", "2.5e-11", "1.25e-11", "6.25e-12")
STATIONARY_STEP_TEXTS = ("1e-15", "5e-16")
REPEAT_INDICES = (0, 23)
BOUNDARY_PIVOT = 3
BOUNDARY_ALT_PIVOT = 9
INTERNAL_PIVOT = 9
INTERNAL_ALT_PIVOT = 6

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=70):
    return mp.nstr(value, digits, strip_zeros=False)


def serialize_matrix(matrix, digits=65):
    return [
        [mp_text(matrix[row, column], digits) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def matrix_from_text(rows):
    return mp.matrix([[mp.mpf(value) for value in row] for row in rows])


def matrix_max(matrix):
    if not matrix.rows or not matrix.cols:
        return mp.mpf(0)
    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows)
        for column in range(matrix.cols)
    )


def matrix_difference(left, right):
    return matrix_max(left - right)


def real_symmetric(matrix):
    return mp.matrix([
        [
            mp.re(matrix[row, column] + matrix[column, row]) / 2
            for column in range(matrix.cols)
        ]
        for row in range(matrix.rows)
    ])


def submatrix(matrix, rows, columns):
    return mp.matrix([
        [matrix[row, column] for column in columns]
        for row in rows
    ])


def load_action_definitions():
    tree = ast.parse(ACTION_SOURCE.read_text(), filename=str(ACTION_SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "mp": mp,
        "np": np,
        "json": json,
        "Path": Path,
        "sha256": sha256,
        "combinations": combinations,
        "permutations": permutations,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "HERE": HERE,
        "ROOT": ROOT,
        "PAIR4": PAIR4,
        "LOCAL_TRIANGLES": LOCAL_TRIANGLES,
        "TAU_TEXT": TAU_TEXT,
        "VARIABLES": VARIABLES,
        "INTERNAL_VARIABLES": INTERNAL_VARIABLES,
        "FD_STEP_TEXTS": STATIONARY_STEP_TEXTS,
        "FD_GATE_TEXT": "1e-60",
        "EXPECTED_F": (2640, 17040, 28800, 14400),
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(ACTION_SOURCE), "exec"), namespace)
    return namespace


def kernel_basis(vector, pivot):
    size = len(vector)
    columns = tuple(index for index in range(size) if index != pivot)
    if vector[pivot] == 0:
        raise ZeroDivisionError("kernel-basis pivot is zero")
    result = mp.matrix(size, size - 1)
    for column, index in enumerate(columns):
        result[index, column] = 1
        result[pivot, column] = -vector[index] / vector[pivot]
    return result, columns


def reversal_matrix():
    result = mp.matrix(12, 12)
    for index in range(6):
        result[index, 6 + index] = 1
        result[6 + index, index] = 1
    return result


def analytic_tangent(base):
    tau = mp.mpf(TAU_TEXT)
    return mp.matrix(
        [-tau * tau / base[("cross",) + pair] for pair in PAIR4]
        + [mp.mpf(1)] * 4
    )


def direct_slice(boundary_basis, internal_basis):
    result = mp.matrix(22, 20)
    for row in range(12):
        for column in range(11):
            result[row, column] = boundary_basis[row, column]
    for row in range(10):
        for column in range(9):
            result[12 + row, 11 + column] = internal_basis[row, column]
    return result


def shifted_coordinates(base, direction, step):
    result = dict(base)
    for index, key in enumerate(VARIABLES):
        result[key] *= mp.exp(step * direction[index])
    return result


def new_branch_tracker():
    return {
        "maximum_angle_identity_residual": mp.mpf(0),
        "maximum_imaginary_curvature": mp.mpf(0),
        "maximum_relative_imaginary_action": mp.mpf(0),
        "minimum_angle_argument": mp.inf,
        "evaluations": 0,
    }


def update_branch_tracker(tracker, record):
    action = record["action"]
    tracker["maximum_angle_identity_residual"] = max(
        tracker["maximum_angle_identity_residual"],
        record["maximum_angle_identity_residual"],
    )
    tracker["maximum_imaginary_curvature"] = max(
        tracker["maximum_imaginary_curvature"],
        record["maximum_imaginary_curvature"],
    )
    tracker["maximum_relative_imaginary_action"] = max(
        tracker["maximum_relative_imaginary_action"],
        abs(mp.im(action)) / max(mp.mpf(1), abs(mp.re(action))),
    )
    tracker["minimum_angle_argument"] = min(
        tracker["minimum_angle_argument"], record["minimum_angle_argument"]
    )
    tracker["evaluations"] += 1


def scalar_schedule_action(combinatorics, geometry, coordinates, masses, tracker):
    """Evaluate only the scalar action; never construct or read its gradient."""
    angle_lookup = {}
    maximum_identity = mp.mpf(0)
    minimum_argument = mp.inf
    for simplex_record in combinatorics["simplex_types"]:
        states = simplex_record["states"]
        squared = actions["simplex_squared"](states, coordinates)
        angles, identity, argument = actions["angle_record"](squared)
        maximum_identity = max(maximum_identity, identity)
        minimum_argument = min(minimum_argument, argument)
        for local_triangle, angle in angles.items():
            triangle_states = actions["canonical_states"](
                states[index] for index in local_triangle
            )
            angle_lookup[states, triangle_states] = angle

    gravitational_sum = mp.mpc(0)
    maximum_imaginary_curvature = mp.mpf(0)
    for triangle_record in combinatorics["triangle_types"]:
        states = triangle_record["states"]
        boundary = len({layer for _, layer in states}) == 1
        curvature = mp.pi if boundary else 2 * mp.pi
        for contribution in triangle_record["contributions"]:
            curvature += contribution["multiplicity"] * angle_lookup[
                contribution["simplex"], contribution["triangle"]
            ]
        values = []
        for left, right in combinations(range(3), 2):
            key, jacobian = actions["edge_coordinate"](states[left], states[right])
            values.append(jacobian * coordinates[key])
        x, y, z = values
        area_square = (2 * (x * y + x * z + y * z) - x * x - y * y - z * z) / 16
        area = mp.sqrt(mp.mpc(area_square))
        gravitational_sum += triangle_record["count"] * area * curvature
        maximum_imaginary_curvature = max(
            maximum_imaginary_curvature, abs(mp.im(curvature))
        )

    gravitational = -mp.j * gravitational_sum
    dust = -8 * mp.pi * mp.fsum(
        masses[rank] * mp.sqrt(coordinates["rho", rank])
        for rank in range(4)
    )
    record = {
        "action": gravitational + dust,
        "maximum_angle_identity_residual": maximum_identity,
        "minimum_angle_argument": minimum_argument,
        "maximum_imaginary_curvature": maximum_imaginary_curvature,
    }
    update_branch_tracker(tracker, record)
    return record["action"]


def extrapolation(values):
    richardson = tuple(
        (4 * values[index + 1] - values[index]) / 3
        for index in range(4)
    )
    sixth = tuple(
        (16 * richardson[index + 1] - richardson[index]) / 15
        for index in range(3)
    )
    eighth = tuple(
        (64 * sixth[index + 1] - sixth[index]) / 63
        for index in range(2)
    )
    return richardson, sixth, eighth


def action_second_ladder(schedule, geometry, base, masses, direction, centre, tracker):
    values = []
    for step_text in STEP_TEXTS:
        step = mp.mpf(step_text)
        upper = scalar_schedule_action(
            schedule, geometry, shifted_coordinates(base, direction, step),
            masses, tracker,
        )
        lower = scalar_schedule_action(
            schedule, geometry, shifted_coordinates(base, direction, -step),
            masses, tracker,
        )
        values.append((upper - 2 * centre + lower) / (step * step))
    richardson, sixth, eighth = extrapolation(tuple(values))
    return {
        "centred": tuple(values),
        "richardson": richardson,
        "sixth": sixth,
        "eighth": eighth,
    }


def reconstruct_matrix(diagonal, pair_values, level):
    result = mp.matrix(20, 20)
    for index in range(20):
        result[index, index] = diagonal[index][level]
    for (left, right), values in pair_values.items():
        value = (
            values[level] - diagonal[left][level] - diagonal[right][level]
        ) / 2
        result[left, right] = result[right, left] = value
    return result


def reduce_matrix(matrix):
    boundary = tuple(range(11))
    internal = tuple(range(11, 20))
    a = submatrix(matrix, boundary, boundary)
    b = submatrix(matrix, boundary, internal)
    c = submatrix(matrix, internal, internal)
    y = mp.matrix(9, 11)
    right = -b.T
    for column in range(11):
        solution = mp.lu_solve(c, mp.matrix([right[row, column] for row in range(9)]))
        for row in range(9):
            y[row, column] = solution[row]
    residual = c * y + b.T
    response = real_symmetric(a + b * y)
    return response, y, c, matrix_max(residual) / max(mp.mpf(1), matrix_max(b))


def direct_reconstruction(schedule, geometry, base, masses, boundary_basis,
                          internal_basis, dps, label):
    with mp.workdps(dps):
        tracker = new_branch_tracker()
        centre = scalar_schedule_action(schedule, geometry, base, masses, tracker)
        slice_matrix = direct_slice(boundary_basis, internal_basis)
        diagonal = []
        pair_values = {}
        max_displacement = mp.mpf(0)
        max_ladder_imaginary = mp.mpf(0)
        max_final_difference = mp.mpf(0)

        for index in range(20):
            direction = slice_matrix[:, index]
            max_displacement = max(
                max_displacement,
                mp.mpf(STEP_TEXTS[0]) * matrix_max(direction),
            )
            ladder = action_second_ladder(
                schedule, geometry, base, masses, direction, centre, tracker
            )
            diagonal.append(ladder["eighth"])
            max_ladder_imaginary = max(
                max_ladder_imaginary,
                *(abs(mp.im(value)) for value in ladder["eighth"]),
            )
            max_final_difference = max(
                max_final_difference,
                abs(ladder["eighth"][1] - ladder["eighth"][0]),
            )

        pair_list = tuple(combinations(range(20), 2))
        for pair_number, (left, right) in enumerate(pair_list, start=1):
            direction = slice_matrix[:, left] + slice_matrix[:, right]
            max_displacement = max(
                max_displacement,
                mp.mpf(STEP_TEXTS[0]) * matrix_max(direction),
            )
            ladder = action_second_ladder(
                schedule, geometry, base, masses, direction, centre, tracker
            )
            pair_values[left, right] = ladder["eighth"]
            max_ladder_imaginary = max(
                max_ladder_imaginary,
                *(abs(mp.im(value)) for value in ladder["eighth"]),
            )
            max_final_difference = max(
                max_final_difference,
                abs(ladder["eighth"][1] - ladder["eighth"][0]),
            )
            if pair_number in (50, 100, 150, 190):
                print(
                    f"[INFO] {label}: {20 + pair_number}/210 directions",
                    flush=True,
                )

        raw0 = reconstruct_matrix(diagonal, pair_values, 0)
        raw1 = reconstruct_matrix(diagonal, pair_values, 1)
        g0 = real_symmetric(raw0)
        g1 = real_symmetric(raw1)
        e_g = (
            100 * matrix_difference(raw1, raw0)
            + mp.mpf("1e-80") * max(mp.mpf(1), matrix_max(raw1))
        )
        k0, y0, c0, solve0 = reduce_matrix(g0)
        k1, y1, c1, solve1 = reduce_matrix(g1)
        e_k = (
            100 * matrix_difference(k1, k0)
            + mp.mpf("1e-60") * max(mp.mpf(1), matrix_max(k1))
        )
        eigenvalues, _ = mp.eigsy(c1)
        return {
            "restricted_second_variation": g1,
            "response": k1,
            "lift_coefficients": y1,
            "internal_eigenvalues": tuple(eigenvalues[row] for row in range(9)),
            "second_variation_envelope": e_g,
            "response_envelope": e_k,
            "solve_residual": max(solve0, solve1),
            "maximum_coordinate_displacement": max_displacement,
            "maximum_ladder_imaginary": max_ladder_imaginary,
            "maximum_eighth_level_difference": max_final_difference,
            "branch": tracker,
        }


def polynomial_value(vector, matrix):
    quadratic = (vector.T * matrix * vector)[0] / 2
    return (
        quadratic
        + 11 * mp.fsum(value ** 4 for value in vector)
        + 13 * mp.fsum(value ** 6 for value in vector)
        + 17 * mp.fsum(value ** 8 for value in vector)
        + 19 * mp.fsum(value ** 10 for value in vector)
    )


def polynomial_control():
    with mp.workdps(ACTION_DPS):
        expected = mp.matrix([[7, 2, -1], [2, 5, 3], [-1, 3, 11]])
        diagonal = []
        pairs = {}
        directions = [mp.eye(3)[:, index] for index in range(3)]

        def ladder(direction):
            values = []
            for step_text in STEP_TEXTS:
                step = mp.mpf(step_text)
                values.append(
                    (
                        polynomial_value(step * direction, expected)
                        - 2 * polynomial_value(mp.matrix(3, 1), expected)
                        + polynomial_value(-step * direction, expected)
                    ) / (step * step)
                )
            return extrapolation(tuple(values))[2]

        for direction in directions:
            diagonal.append(ladder(direction))
        for left, right in combinations(range(3), 2):
            pairs[left, right] = ladder(directions[left] + directions[right])

        matrices = []
        wrong = None
        for level in range(2):
            result = mp.matrix(3, 3)
            bad = mp.matrix(3, 3)
            for index in range(3):
                result[index, index] = bad[index, index] = diagonal[index][level]
            for (left, right), values in pairs.items():
                numerator = values[level] - diagonal[left][level] - diagonal[right][level]
                result[left, right] = result[right, left] = numerator / 2
                bad[left, right] = bad[right, left] = numerator
            matrices.append(result)
            if level == 1:
                wrong = bad
        envelope = (
            100 * matrix_difference(matrices[1], matrices[0])
            + mp.mpf("1e-80") * max(mp.mpf(1), matrix_max(matrices[1]))
        )
        return {
            "expected": expected,
            "computed": matrices[1],
            "envelope": envelope,
            "error": matrix_difference(matrices[1], expected),
            "wrong_factor_error": matrix_difference(wrong, expected),
        }


def transformed_envelope(matrix, envelope):
    maximum_column_sum = max(
        mp.fsum(abs(matrix[row, column]) for row in range(matrix.rows))
        for column in range(matrix.cols)
    )
    return maximum_column_sum ** 2 * envelope


def classify(matrices, envelopes):
    representatives = []
    memberships = []
    assignments = []
    for index, matrix in enumerate(matrices):
        assigned = None
        for class_index, representative in enumerate(representatives):
            if matrix_difference(matrix, matrices[representative]) <= (
                envelopes[index] + envelopes[representative]
            ):
                assigned = class_index
                break
        if assigned is None:
            assigned = len(representatives)
            representatives.append(index)
            memberships.append([])
        memberships[assigned].append(index)
        assignments.append(assigned)
    return representatives, memberships, assignments


print("=" * 78)
print("DIRECT-ACTION ADVERSARIAL CONSTRAINED H4 RESPONSE")
print("=" * 78)

paths = {
    "action_source": ACTION_SOURCE,
    "curvature": CURVATURE,
    "feasibility": FEASIBILITY,
    "null_coupling": NULL_COUPLING,
    "primary": PRIMARY,
    "primary_result": PRIMARY_RESULT,
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "cell600": CELL600,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all direct-action, geometry, protocol and frozen-comparison inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES and PROTOCOL_COMMIT == "a9ee74b",
)

curvature = json.loads(CURVATURE.read_text())
null_coupling = json.loads(NULL_COUPLING.read_text())
upstream_ok = check(
    "the mass and compatibility inputs carry their accepted scoped outcomes",
    curvature["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and null_coupling["outcome"]
        == "ADVERSARIAL_REFINED_H4_NULL_COUPLING_CORROBORATED"
    and null_coupling["compatibility"]["rank"] == 1,
)

actions = load_action_definitions()
definitions_ok = check(
    "the scalar route loads definitions only and does not call the gradient/Hessian routine",
    {
        "tetrahedra_from_adjacency",
        "barycentric_chambers",
        "all_simplices",
        "schedule_combinatorics",
        "exact_geometry",
        "base_coordinates",
        "angle_record",
        "simplex_squared",
        "edge_coordinate",
    } <= set(actions)
    and "OUTPUT" not in actions,
)

_, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
_, top, colours = actions["barycentric_chambers"](coarse_top)
spatial_cells = actions["all_simplices"](tuple(map(tuple, top)))
orders = tuple(permutations(range(4)))
combinatorics = tuple(
    actions["schedule_combinatorics"](top, colours, order) for order in orders
)
topology_ok = check(
    "the exact carrier and complete 24-schedule census are reconstructed",
    tuple(len(layer) for layer in spatial_cells) == (2640, 17040, 28800, 14400)
    and len(combinatorics) == 24
    and all(
        record["pentachora"] == record["distinct_pentachora"] == 57600
        and record["triangles"] == 149280
        and record["boundary_triangles"] == 57600
        and record["mixed_triangle_types"] == 0
        and len(record["simplex_types"]) == 4
        and len(record["triangle_types"]) == 28
        for record in combinatorics
    ),
)

mp.mp.dps = REPEAT_DPS + 10
geometry180 = actions["exact_geometry"](ACTION_DPS)
geometry220 = actions["exact_geometry"](REPEAT_DPS)
geometry180["mass"] = mp.mpf(0)
geometry220["mass"] = mp.mpf(0)
base180 = actions["base_coordinates"](geometry180)
base220 = actions["base_coordinates"](geometry220)
masses = tuple(
    mp.mpf(value) for value in curvature["selected_rank_matter"]["total_masses"]
)
c = mp.matrix([
    mp.mpf(value)
    for value in null_coupling["compatibility"]["adversarial_row"]
])
n = analytic_tangent(base220)
p, p_columns = kernel_basis(c, BOUNDARY_PIVOT)
p_alt, p_alt_columns = kernel_basis(c, BOUNDARY_ALT_PIVOT)
q, q_columns = kernel_basis(n, INTERNAL_PIVOT)
q_alt, q_alt_columns = kernel_basis(n, INTERNAL_ALT_PIVOT)
p_change = submatrix(p_alt, p_columns, tuple(range(11)))
q_change = submatrix(q_alt, q_columns, tuple(range(9)))
reversal = reversal_matrix()
rp = reversal * p
t_reversal = submatrix(rp, p_columns, tuple(range(11)))
basis_error = max(
    matrix_max(c.T * p),
    matrix_max(c.T * p_alt),
    matrix_max(n.T * q),
    matrix_max(n.T * q_alt),
    matrix_difference(p * p_change, p_alt),
    matrix_difference(q * q_change, q_alt),
    matrix_difference(p * t_reversal, rp),
)
bases_ok = check(
    "the algebraic compatibility, quotient and reversal bases satisfy exact identities",
    p.rows == 12 and p.cols == 11
    and q.rows == 10 and q.cols == 9
    and basis_error < mp.mpf("1e-180"),
    f"max identity error={mp_text(basis_error, 8)}",
)

stationarity_tracker = new_branch_tracker()
stationarity_residuals = []
with mp.workdps(ACTION_DPS):
    for schedule_index, schedule in enumerate(combinatorics):
        centre = scalar_schedule_action(
            schedule, geometry180, base180, masses, stationarity_tracker
        )
        for internal_index in range(10):
            direction = mp.matrix(22, 1)
            direction[12 + internal_index] = 1
            derivatives = []
            for step_text in STATIONARY_STEP_TEXTS:
                step = mp.mpf(step_text)
                upper = scalar_schedule_action(
                    schedule, geometry180,
                    shifted_coordinates(base180, direction, step),
                    masses, stationarity_tracker,
                )
                lower = scalar_schedule_action(
                    schedule, geometry180,
                    shifted_coordinates(base180, direction, -step),
                    masses, stationarity_tracker,
                )
                derivatives.append((upper - lower) / (2 * step))
            stationarity_residuals.append((4 * derivatives[1] - derivatives[0]) / 3)
        if schedule_index in (5, 11, 17, 23):
            print(f"[INFO] direct stationarity: {schedule_index + 1}/24 schedules", flush=True)
maximum_stationarity_residual = max(abs(value) for value in stationarity_residuals)
stationarity_ok = check(
    "all ten internal equations vanish directly from the scalar action on all schedules",
    len(stationarity_residuals) == 240
    and maximum_stationarity_residual < mp.mpf("1e-60"),
    f"max residual={mp_text(maximum_stationarity_residual, 8)}",
)

polynomial = polynomial_control()
polynomial_ok = check(
    "the known polynomial Hessian is recovered inside the frozen direct envelope",
    polynomial["error"] <= polynomial["envelope"],
    f"error/envelope={mp_text(polynomial['error'] / polynomial['envelope'], 8)}",
)
wrong_factor_ok = check(
    "the known polynomial control rejects omission of the polarization factor 1/2",
    polynomial["wrong_factor_error"] > mp.mpf("1e6") * polynomial["envelope"],
    f"wrong/envelope={mp_text(polynomial['wrong_factor_error'] / polynomial['envelope'], 8)}",
)

direct_records = []
print("[INFO] reconstructing 24 complete 20x20 forms from scalar action", flush=True)
for index, schedule in enumerate(combinatorics):
    record = direct_reconstruction(
        schedule, geometry180, base180, masses, p, q, ACTION_DPS,
        f"schedule {index + 1}/24 {orders[index]}",
    )
    direct_records.append(record)
    print(
        f"[INFO] direct schedules completed: {index + 1}/24; "
        f"e_K={mp_text(record['response_envelope'], 5)}",
        flush=True,
    )

branch_trackers = [stationarity_tracker] + [record["branch"] for record in direct_records]
branch_ok = check(
    "every direct scalar-action evaluation remains on the same real Lorentzian branch",
    max(item["maximum_angle_identity_residual"] for item in branch_trackers)
        < mp.mpf("1e-80")
    and max(item["maximum_imaginary_curvature"] for item in branch_trackers)
        < mp.mpf("1e-80")
    and max(item["maximum_relative_imaginary_action"] for item in branch_trackers)
        < mp.mpf("1e-80")
    and min(item["minimum_angle_argument"] for item in branch_trackers)
        > mp.mpf("1e-8"),
)

direct_numerics_ok = check(
    "all 24 direct internal blocks, solves, displacements and imaginary errors pass",
    all(
        min(record["internal_eigenvalues"])
            > 100 * record["second_variation_envelope"]
        and record["solve_residual"] < mp.mpf("1e-80")
        and record["maximum_coordinate_displacement"] < mp.mpf("2e-5")
        and record["maximum_ladder_imaginary"]
            <= record["second_variation_envelope"]
        for record in direct_records
    ),
    f"min eigen={mp_text(min(min(r['internal_eigenvalues']) for r in direct_records), 8)}, "
    f"max solve={mp_text(max(r['solve_residual'] for r in direct_records), 8)}",
)

order_index = {order: index for index, order in enumerate(orders)}
reversal_differences = []
reversal_envelopes = []
canonical_matrices = []
canonical_envelopes = []
for index, order in enumerate(orders):
    reverse_index = order_index[tuple(reversed(order))]
    reversed_response = (
        t_reversal.T * direct_records[reverse_index]["response"] * t_reversal
    )
    reversal_differences.append(
        matrix_difference(direct_records[index]["response"], reversed_response)
    )
    reversal_envelopes.append(
        direct_records[index]["response_envelope"]
        + transformed_envelope(
            t_reversal, direct_records[reverse_index]["response_envelope"]
        )
    )
    if order <= tuple(reversed(order)):
        canonical_matrices.append(direct_records[index]["response"])
        canonical_envelopes.append(direct_records[index]["response_envelope"])
    else:
        canonical_matrices.append(
            t_reversal.T * direct_records[index]["response"] * t_reversal
        )
        canonical_envelopes.append(
            transformed_envelope(t_reversal, direct_records[index]["response_envelope"])
        )

time_reversal_covariant = all(
    difference <= envelope
    for difference, envelope in zip(reversal_differences, reversal_envelopes)
)
representatives, memberships, assignments = classify(
    canonical_matrices, canonical_envelopes
)
print("[INFO] complete direct schedule classes:", flush=True)
for class_index, members in enumerate(memberships):
    print(
        f"[INFO] class {class_index}: "
        + ", ".join(str(orders[index]) for index in members),
        flush=True,
    )
class_census_ok = check(
    "the direct time-reversal test and complete target-free class census are resolved",
    len(assignments) == 24 and len(memberships) >= 1,
    f"reversal={time_reversal_covariant}, classes={len(memberships)}",
)

repeat_records = {}
for index in REPEAT_INDICES:
    repeat_records[index] = direct_reconstruction(
        combinatorics[index], geometry220, base220, masses, p, q, REPEAT_DPS,
        f"precision repeat {index} {orders[index]}",
    )
    print(f"[INFO] precision repeat completed: schedule {index}", flush=True)
repeat_differences = {
    index: matrix_difference(
        direct_records[index]["response"], repeat_records[index]["response"]
    )
    for index in REPEAT_INDICES
}
repeat_gates = {
    index: 10 * (
        direct_records[index]["response_envelope"]
        + repeat_records[index]["response_envelope"]
    )
    for index in REPEAT_INDICES
}
precision_ok = check(
    "the two complete 220-digit reconstructions agree with the 180-digit route",
    all(repeat_differences[index] <= repeat_gates[index] for index in REPEAT_INDICES),
    f"max fraction={mp_text(max(
        repeat_differences[i] / repeat_gates[i] for i in REPEAT_INDICES
    ), 8)}",
)

# The primary JSON is deliberately loaded only after all direct matrices and
# their target-free class census have been completed.
primary = json.loads(PRIMARY.read_text())
primary_upstream_ok = check(
    "the delayed comparison target is exactly the frozen clean primary result",
    primary["outcome"] == "REFINED_H4_CONSTRAINED_RESPONSE_SINGLE_SCHEDULE_CLASS"
    and primary["tests"] == {"passed": 19, "total": 19}
    and primary["census"]["schedule_count"] == 24
    and primary["census"]["class_count"] == 1,
)
primary_differences = []
primary_gates = []
for index in range(24):
    primary_record = primary["census"]["schedules"][index]
    primary_matrix = matrix_from_text(primary_record["primary_response"])
    primary_envelope = mp.mpf(primary_record["response_envelopes"]["primary"])
    primary_differences.append(
        matrix_difference(direct_records[index]["response"], primary_matrix)
    )
    primary_gates.append(
        10 * (direct_records[index]["response_envelope"] + primary_envelope)
    )
primary_matches = [
    difference <= gate
    for difference, gate in zip(primary_differences, primary_gates)
]
primary_comparison_ok = check(
    "all 24 direct scalar-action matrices match the analytic-gradient primary matrices",
    all(primary_matches),
    f"matches={sum(primary_matches)}/24, max fraction={mp_text(max(
        difference / gate for difference, gate in zip(primary_differences, primary_gates)
    ), 8)}",
)

corruption_size = mp.mpf("1e-6") * max(mp.mpf(1), matrix_max(canonical_matrices[0]))
corrupted = mp.matrix(canonical_matrices[0])
corrupted[0, 0] += corruption_size
corrupt_representatives, corrupt_memberships, _ = classify(
    [canonical_matrices[0], corrupted],
    [canonical_envelopes[0], canonical_envelopes[0]],
)
corrupt_primary_difference = matrix_difference(
    corrupted,
    matrix_from_text(primary["census"]["schedules"][0]["primary_response"]),
)
corruption_ok = check(
    "a resolved one-entry corruption splits the class and fails primary comparison",
    len(corrupt_representatives) == len(corrupt_memberships) == 2
    and corrupt_primary_difference > primary_gates[0],
    f"corruption/gate={mp_text(corrupt_primary_difference / primary_gates[0], 8)}",
)

naive_reversal_pass_count = 0
for index, order in enumerate(orders):
    reverse_index = order_index[tuple(reversed(order))]
    difference = matrix_difference(
        direct_records[index]["response"], direct_records[reverse_index]["response"]
    )
    envelope = (
        direct_records[index]["response_envelope"]
        + direct_records[reverse_index]["response_envelope"]
    )
    naive_reversal_pass_count += int(difference <= envelope)
naive_convention_ok = check(
    "the alternative untransformed reversal convention is reported without selection",
    0 <= naive_reversal_pass_count <= 24,
    f"naive passes={naive_reversal_pass_count}/24",
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    definitions_ok,
    topology_ok,
    bases_ok,
    stationarity_ok,
    polynomial_ok,
    wrong_factor_ok,
    branch_ok,
    direct_numerics_ok,
    class_census_ok,
    precision_ok,
    primary_upstream_ok,
    corruption_ok,
    naive_convention_ok,
))
if not controls_ok:
    outcome = "ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED"
elif (
    not time_reversal_covariant
    or len(memberships) != 1
    or not primary_comparison_ok
):
    outcome = "ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_DISAGREEMENT"
else:
    outcome = "ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_CORROBORATED"

expected_corroboration = (
    controls_ok
    and time_reversal_covariant
    and len(memberships) == 1
    and primary_comparison_ok
)
outcome_ok = check(
    "the frozen outcome follows mechanically from controls and direct comparisons",
    (outcome.endswith("_CORROBORATED")) == expected_corroboration
    and (outcome.endswith("_DISAGREEMENT")) == (
        controls_ok and not expected_corroboration
    )
    and (outcome.endswith("_CONTROL_FAILED")) == (not controls_ok),
    outcome,
)

all_trackers = branch_trackers + [repeat_records[index]["branch"] for index in REPEAT_INDICES]
artifact = {
    "title": "Direct-action adversarial constrained H4 response",
    "date": "2026-08-21",
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "carrier": "P(sd K_600)",
        "coordinates": "total-orbit log squared-edge coordinates",
        "schedule_count": 24,
        "direct_slice_dimension": 20,
        "boundary_compatibility_dimension": 11,
        "internal_quotient_dimension": 9,
        "directions_per_complete_reconstruction": 210,
        "action_decimal_precision": ACTION_DPS,
        "repeat_decimal_precision": REPEAT_DPS,
        "steps": STEP_TEXTS,
        "repeat_indices": REPEAT_INDICES,
        "primary_loaded_after_direct_census": True,
        "analytic_gradient_or_hessian_used": False,
    },
    "controls": {
        "maximum_stationarity_residual": mp_text(maximum_stationarity_residual),
        "basis_identity_error": mp_text(basis_error),
        "polynomial_error": mp_text(polynomial["error"]),
        "polynomial_envelope": mp_text(polynomial["envelope"]),
        "wrong_polarization_error": mp_text(polynomial["wrong_factor_error"]),
        "precision_repeat": {
            str(index): {
                "difference": mp_text(repeat_differences[index]),
                "gate": mp_text(repeat_gates[index]),
                "response_envelope_220": mp_text(
                    repeat_records[index]["response_envelope"]
                ),
            }
            for index in REPEAT_INDICES
        },
        "corruption_size": mp_text(corruption_size),
        "corruption_class_count": len(corrupt_memberships),
        "corruption_primary_difference": mp_text(corrupt_primary_difference),
        "naive_untransformed_reversal_pass_count": naive_reversal_pass_count,
    },
    "branch": {
        "evaluation_count": sum(item["evaluations"] for item in all_trackers),
        "maximum_angle_identity_residual": mp_text(max(
            item["maximum_angle_identity_residual"] for item in all_trackers
        )),
        "maximum_imaginary_curvature": mp_text(max(
            item["maximum_imaginary_curvature"] for item in all_trackers
        )),
        "maximum_relative_imaginary_action": mp_text(max(
            item["maximum_relative_imaginary_action"] for item in all_trackers
        )),
        "minimum_angle_argument": mp_text(min(
            item["minimum_angle_argument"] for item in all_trackers
        )),
    },
    "census": {
        "class_count": len(memberships),
        "classes": [
            {
                "indices": members,
                "orders": [orders[index] for index in members],
            }
            for members in memberships
        ],
        "time_reversal_covariant": time_reversal_covariant,
        "maximum_reversal_difference": mp_text(max(reversal_differences)),
        "maximum_reversal_fraction": mp_text(max(
            difference / envelope
            for difference, envelope in zip(reversal_differences, reversal_envelopes)
        )),
        "primary_match_count": sum(primary_matches),
        "maximum_primary_difference": mp_text(max(primary_differences)),
        "maximum_primary_fraction": mp_text(max(
            difference / gate
            for difference, gate in zip(primary_differences, primary_gates)
        )),
        "schedules": [
            {
                "order": orders[index],
                "canonical_class": assignments[index],
                "restricted_second_variation": serialize_matrix(
                    record["restricted_second_variation"]
                ),
                "response": serialize_matrix(record["response"]),
                "lift_coefficients": serialize_matrix(record["lift_coefficients"]),
                "internal_eigenvalues": [
                    mp_text(value) for value in record["internal_eigenvalues"]
                ],
                "second_variation_envelope": mp_text(
                    record["second_variation_envelope"]
                ),
                "response_envelope": mp_text(record["response_envelope"]),
                "solve_residual": mp_text(record["solve_residual"]),
                "maximum_coordinate_displacement": mp_text(
                    record["maximum_coordinate_displacement"]
                ),
                "maximum_ladder_imaginary": mp_text(
                    record["maximum_ladder_imaginary"]
                ),
                "primary_difference": mp_text(primary_differences[index]),
                "primary_gate": mp_text(primary_gates[index]),
                "primary_match": primary_matches[index],
            }
            for index, record in enumerate(direct_records)
        ],
    },
    "scope": {
        "nonlinear_constraint_surface_computed": False,
        "nonhomogeneous_operator_or_spectrum_computed": False,
        "root_search_or_deferred_census_executed": False,
        "continuum_or_particle_target_loaded": False,
        "physical_constant_extracted": False,
    },
    "status_labels": {
        "direct_H4_response": (
            "DERIVED COMPUTATIONAL" if outcome.endswith("_CORROBORATED") else "OPEN"
        ),
        "external_novelty": "OPEN",
        "nonhomogeneous_propagation": "OPEN",
        "tick_c_G_planck_particles": "NOT ESTABLISHED",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Direct classes: {len(memberships)}")
print(f"Primary matches: {sum(primary_matches)}/24")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")

if passed != tests:
    raise SystemExit(1)


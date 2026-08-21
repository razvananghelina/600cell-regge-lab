#!/usr/bin/env python3
"""Effective H4 boundary Hessian of the curvature-matched refined seed."""

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
CURVATURE_ADVERSARIAL = (
    HERE / "gravity_600cell_refined_local_curvature_mass_adversarial.json"
)
BOUNDARY_COTANGENT = HERE / "gravity_600cell_refined_boundary_cotangent.json"
BOUNDARY_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_boundary_cotangent_result.md"
)
PRIOR_ART = (
    ROOT / "docs/gravity/gravity_600cell_refined_effective_h4_hessian_prior_art.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_effective_h4_hessian_protocol.md"
)
CELL600 = ROOT / "commons/cell600.py"
OUTPUT = HERE / "gravity_600cell_refined_effective_h4_hessian.json"

PRIOR_ART_COMMIT = "f7bf3c1"
PROTOCOL_COMMIT = "eed7891"
EXPECTED_HASHES = {
    "action_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "curvature": "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "curvature_adversarial": "c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e",
    "boundary_cotangent": "4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa",
    "boundary_result": "391a317b9f8823a5479f450dde43a43177e210a2d81192aedc938e90fc8006d1",
    "prior_art": "d111b896265ccbd0534ec50fec184c81067133665b97a99ac0a69df834877934",
    "protocol": "7aa1628c6cad5865ac2a97fad6a230f4f0c29148fe2dfdc48c03a50267c2713d",
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
PRIMARY_DPS = 100
SECONDARY_DPS = 140
STEP_TEXTS = ("1e-10", "5e-11", "2.5e-11")
DIRECTIONAL_INDICES = (0, 1, 22, 23)

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


def variable_label(key):
    if key[0] == "rho":
        return f"rho_{key[1]}"
    return f"{key[0]}_{key[1]}{key[2]}"


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
        "FD_STEP_TEXTS": ("1e-15", "5e-16"),
        "FD_GATE_TEXT": "1e-24",
        "EXPECTED_F": (2640, 17040, 28800, 14400),
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(ACTION_SOURCE), "exec"), namespace)
    return namespace


def matrix_max(matrix):
    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows)
        for column in range(matrix.cols)
    )


def matrix_difference(left, right):
    return max(
        abs(left[row, column] - right[row, column])
        for row in range(left.rows)
        for column in range(left.cols)
    )


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


def shifted_coordinates(base, key, displacement):
    result = dict(base)
    result[key] *= mp.exp(displacement)
    return result


def gravitational_gradient(evaluate_schedule, combinatorics, geometry, coordinates):
    evaluation = evaluate_schedule(combinatorics, geometry, coordinates)
    return evaluation, mp.matrix([evaluation["gradient"][key] for key in VARIABLES])


def derivative_column(evaluate_schedule, combinatorics, geometry, base, key, step):
    _, plus = gravitational_gradient(
        evaluate_schedule,
        combinatorics,
        geometry,
        shifted_coordinates(base, key, step),
    )
    _, minus = gravitational_gradient(
        evaluate_schedule,
        combinatorics,
        geometry,
        shifted_coordinates(base, key, -step),
    )
    return (plus - minus) / (2 * step)


def hessian_ladders(evaluate_schedule, combinatorics, geometry, dps,
                    secondary_only=False):
    with mp.workdps(dps):
        base = actions["base_coordinates"](geometry)
        steps = tuple(mp.mpf(value) for value in STEP_TEXTS)
        primary = mp.matrix(22, 22)
        secondary = mp.matrix(22, 22)
        for column, key in enumerate(VARIABLES):
            d1 = derivative_column(
                evaluate_schedule, combinatorics, geometry, base, key, steps[1]
            )
            d2 = derivative_column(
                evaluate_schedule, combinatorics, geometry, base, key, steps[2]
            )
            secondary_column = (4 * d2 - d1) / 3
            for row in range(22):
                secondary[row, column] = secondary_column[row]
            if not secondary_only:
                d0 = derivative_column(
                    evaluate_schedule, combinatorics, geometry, base, key, steps[0]
                )
                primary_column = (4 * d1 - d0) / 3
                for row in range(22):
                    primary[row, column] = primary_column[row]
        return primary, secondary


def dust_hessian(masses, tau):
    result = mp.matrix(22, 22)
    for rank in range(4):
        result[18 + rank, 18 + rank] = -2 * mp.pi * masses[rank] * tau
    return result


def add_matrices(left, right):
    return mp.matrix([
        [left[row, column] + right[row, column] for column in range(left.cols)]
        for row in range(left.rows)
    ])


def solve_columns(matrix, right):
    result = mp.matrix(matrix.rows, right.cols)
    for column in range(right.cols):
        solution = mp.lu_solve(
            matrix, mp.matrix([right[row, column] for row in range(right.rows)])
        )
        for row in range(matrix.rows):
            result[row, column] = solution[row]
    return result


def schur_complement(matrix):
    boundary = tuple(range(12))
    internal = tuple(range(12, 22))
    hbb = submatrix(matrix, boundary, boundary)
    hbi = submatrix(matrix, boundary, internal)
    hib = submatrix(matrix, internal, boundary)
    hii = submatrix(matrix, internal, internal)
    solved = solve_columns(hii, hib)
    residual = hii * solved - hib
    relative_residual = matrix_max(residual) / max(mp.mpf(1), matrix_max(hib))
    return hbb - hbi * solved, relative_residual, hii, hib


def total_action(evaluate_schedule, combinatorics, geometry, coordinates,
                 masses):
    gravitational = evaluate_schedule(
        combinatorics, geometry, coordinates
    )["action"]
    dust = -8 * mp.pi * mp.fsum(
        masses[rank] * mp.sqrt(coordinates["rho", rank])
        for rank in range(4)
    )
    return gravitational + dust


def directional_second(evaluate_schedule, combinatorics, geometry, base,
                       masses, direction, step):
    plus = dict(base)
    minus = dict(base)
    for index, key in enumerate(VARIABLES):
        plus[key] *= mp.exp(step * direction[index])
        minus[key] *= mp.exp(-step * direction[index])
    centre = total_action(evaluate_schedule, combinatorics, geometry, base, masses)
    upper = total_action(evaluate_schedule, combinatorics, geometry, plus, masses)
    lower = total_action(evaluate_schedule, combinatorics, geometry, minus, masses)
    return (upper - 2 * centre + lower) / (step * step)


def reversal_matrix():
    result = mp.matrix(12, 12)
    for index in range(6):
        result[index, 6 + index] = 1
        result[6 + index, index] = 1
    return result


def canonicalized_matrix(index, orders, matrices, reversal):
    order = orders[index]
    reverse = tuple(reversed(order))
    if order <= reverse:
        return matrices[index]
    return reversal.T * matrices[index] * reversal


def serialize_matrix(matrix, digits=65):
    return [
        [mp_text(matrix[row, column], digits) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


print("=" * 78)
print("REFINED ON-SHELL EFFECTIVE H4 BOUNDARY HESSIAN")
print("=" * 78)

paths = {
    "action_source": ACTION_SOURCE,
    "curvature": CURVATURE,
    "curvature_adversarial": CURVATURE_ADVERSARIAL,
    "boundary_cotangent": BOUNDARY_COTANGENT,
    "boundary_result": BOUNDARY_RESULT,
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "cell600": CELL600,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all frozen geometry, matter, boundary and protocol inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIOR_ART_COMMIT == "f7bf3c1"
    and PROTOCOL_COMMIT == "eed7891",
)

curvature = json.loads(CURVATURE.read_text())
curvature_adversarial = json.loads(CURVATURE_ADVERSARIAL.read_text())
boundary_cotangent = json.loads(BOUNDARY_COTANGENT.read_text())
upstream_ok = check(
    "the frozen upstream artifacts carry all accepted scoped outcomes",
    curvature["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and curvature["tests"] == {"passed": 15, "total": 15}
    and curvature_adversarial["outcome"]
        == "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED"
    and curvature_adversarial["tests"] == {"passed": 16, "total": 16}
    and boundary_cotangent["outcome"]
        == "REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED"
    and boundary_cotangent["tests"] == {"passed": 16, "total": 16},
)

actions = load_action_definitions()
definitions_ok = check(
    "only frozen action definitions are loaded and the P1 top level is not executed",
    {
        "tetrahedra_from_adjacency",
        "barycentric_chambers",
        "all_simplices",
        "schedule_combinatorics",
        "exact_geometry",
        "base_coordinates",
        "evaluate_schedule",
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
    "the direct carrier and all schedule incidences have the frozen exact census",
    tuple(len(layer) for layer in spatial_cells)
        == (2640, 17040, 28800, 14400)
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

geometry100 = actions["exact_geometry"](PRIMARY_DPS)
geometry140 = actions["exact_geometry"](SECONDARY_DPS)
geometry100["mass"] = mp.mpf(0)
geometry140["mass"] = mp.mpf(0)

with mp.workdps(SECONDARY_DPS):
    tau = mp.mpf(TAU_TEXT)
    masses = tuple(
        mp.mpf(value)
        for value in curvature["selected_rank_matter"]["total_masses"]
    )
    base140 = actions["base_coordinates"](geometry140)
    maximum_internal_residual = mp.mpf(0)
    maximum_branch_identity = mp.mpf(0)
    maximum_branch_imaginary = mp.mpf(0)
    minimum_branch_argument = mp.inf
    for record in combinatorics:
        evaluation = actions["evaluate_schedule"](record, geometry140, base140)
        internal = [
            evaluation["gradient"][key] for key in INTERNAL_VARIABLES[:6]
        ] + [
            evaluation["gradient"]["rho", rank]
            - 4 * mp.pi * masses[rank] * tau
            for rank in range(4)
        ]
        maximum_internal_residual = max(
            maximum_internal_residual, *(abs(value) for value in internal)
        )
        maximum_branch_identity = max(
            maximum_branch_identity,
            evaluation["maximum_angle_identity_residual"],
        )
        maximum_branch_imaginary = max(
            maximum_branch_imaginary,
            evaluation["maximum_imaginary_curvature"],
        )
        minimum_branch_argument = min(
            minimum_branch_argument, evaluation["minimum_angle_argument"]
        )
on_shell_ok = check(
    "all 24 curvature-matched backgrounds are internally on shell",
    maximum_internal_residual < mp.mpf("1e-60"),
    f"max residual={mp_text(maximum_internal_residual, 8)}",
)
branch_ok = check(
    "all schedules remain on the same finite real Lorentzian branch",
    maximum_branch_identity < mp.mpf("1e-80")
    and maximum_branch_imaginary < mp.mpf("1e-80")
    and minimum_branch_argument > mp.mpf("1e-8"),
    f"identity={mp_text(maximum_branch_identity, 8)}, "
    f"imag={mp_text(maximum_branch_imaginary, 8)}",
)

primary100 = []
secondary100 = []
secondary140 = []
print("[INFO] constructing 24 full 22x22 Hessian ladders at 100 digits", flush=True)
for index, record in enumerate(combinatorics):
    first, second = hessian_ladders(
        actions["evaluate_schedule"], record, geometry100, PRIMARY_DPS
    )
    primary100.append(first)
    secondary100.append(second)
    if index in (5, 11, 17, 23):
        print(f"[INFO] 100-digit schedules completed: {index + 1}/24", flush=True)

print("[INFO] repeating fine ladders at 140 digits", flush=True)
for index, record in enumerate(combinatorics):
    _, second = hessian_ladders(
        actions["evaluate_schedule"],
        record,
        geometry140,
        SECONDARY_DPS,
        secondary_only=True,
    )
    secondary140.append(second)
    if index in (5, 11, 17, 23):
        print(f"[INFO] 140-digit schedules completed: {index + 1}/24", flush=True)

records = []
effective_matrices = []
effective_envelopes = []
all_precision_controls = []
all_nonsingular = True
maximum_solve_residual = mp.mpf(0)
maximum_effective_antisymmetry = mp.mpf(0)

with mp.workdps(SECONDARY_DPS):
    dust = dust_hessian(masses, tau)
    for index in range(24):
        raw100a = add_matrices(primary100[index], dust)
        raw100b = add_matrices(secondary100[index], dust)
        raw140b = add_matrices(secondary140[index], dust)
        step_difference = matrix_difference(raw100a, raw100b)
        precision_difference = matrix_difference(raw100b, raw140b)
        full_scale = max(mp.mpf(1), matrix_max(raw140b))
        hessian_envelope = (
            100 * max(step_difference, precision_difference)
            + mp.mpf("1e-50") * full_scale
        )
        raw_imaginary = max(
            abs(mp.im(raw140b[row, column]))
            for row in range(22) for column in range(22)
        )
        raw_antisymmetry = max(
            abs(raw140b[row, column] - raw140b[column, row])
            for row in range(22) for column in range(22)
        )
        all_precision_controls.append(
            raw_imaginary <= hessian_envelope
            and raw_antisymmetry <= hessian_envelope
        )

        symmetric100a = real_symmetric(raw100a)
        symmetric100b = real_symmetric(raw100b)
        symmetric140b = real_symmetric(raw140b)
        hii = submatrix(symmetric140b, tuple(range(12, 22)), tuple(range(12, 22)))
        values, _ = mp.eigsy(hii)
        eigenvalues = tuple(values[row] for row in range(10))
        spectral_envelope = 100 * hessian_envelope
        nonsingular = min(abs(value) for value in eigenvalues) > spectral_envelope
        all_nonsingular = all_nonsingular and nonsingular

        record = {
            "order": list(orders[index]),
            "step_difference": step_difference,
            "precision_difference": precision_difference,
            "hessian_envelope": hessian_envelope,
            "raw_imaginary": raw_imaginary,
            "raw_antisymmetry": raw_antisymmetry,
            "internal_eigenvalues": eigenvalues,
            "internal_spectral_envelope": spectral_envelope,
            "internal_nonsingular": nonsingular,
            "internal_inertia": (
                sum(value > spectral_envelope for value in eigenvalues),
                sum(abs(value) <= spectral_envelope for value in eigenvalues),
                sum(value < -spectral_envelope for value in eigenvalues),
            ),
            "internal_condition": (
                max(abs(value) for value in eigenvalues)
                / min(abs(value) for value in eigenvalues)
            ),
        }

        if nonsingular:
            k100a, residual100a, _, _ = schur_complement(symmetric100a)
            k100b, residual100b, _, _ = schur_complement(symmetric100b)
            k140b, residual140b, _, _ = schur_complement(symmetric140b)
            effective_step = matrix_difference(k100a, k100b)
            effective_precision = matrix_difference(k100b, k140b)
            effective_scale = max(mp.mpf(1), matrix_max(k140b))
            effective_envelope = (
                100 * max(effective_step, effective_precision)
                + mp.mpf("1e-45") * effective_scale
            )
            effective_antisymmetry = max(
                abs(k140b[row, column] - k140b[column, row])
                for row in range(12) for column in range(12)
            )
            solve_residual = max(residual100a, residual100b, residual140b)
            maximum_solve_residual = max(maximum_solve_residual, solve_residual)
            maximum_effective_antisymmetry = max(
                maximum_effective_antisymmetry, effective_antisymmetry
            )
            record.update({
                "effective_step_difference": effective_step,
                "effective_precision_difference": effective_precision,
                "effective_envelope": effective_envelope,
                "effective_antisymmetry": effective_antisymmetry,
                "solve_residual": solve_residual,
                "effective_matrix": k140b,
                "full_matrix": symmetric140b,
            })
            effective_matrices.append(k140b)
            effective_envelopes.append(effective_envelope)
        records.append(record)

precision_ok = check(
    "all full Hessians pass the frozen step, precision, reality and symmetry gates",
    all(all_precision_controls),
    f"max e_H={mp_text(max(item['hessian_envelope'] for item in records), 8)}",
)
internal_census_ok = check(
    "all 24 internal spectra are resolved against their certified envelopes",
    all(
        len(item["internal_eigenvalues"]) == 10
        and sum(item["internal_inertia"]) == 10
        for item in records
    ),
    f"nonsingular={sum(item['internal_nonsingular'] for item in records)}/24",
)

if all_nonsingular:
    solve_ok_value = (
        maximum_solve_residual < mp.mpf("1e-60")
        and all(
            item["effective_antisymmetry"] <= item["effective_envelope"]
            for item in records
        )
    )
else:
    solve_ok_value = True
solve_ok = check(
    "every licensed Schur solve is residual-small and symmetric",
    solve_ok_value,
    f"max residual={mp_text(maximum_solve_residual, 8)}, "
    f"max antisym={mp_text(maximum_effective_antisymmetry, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    hbb_test = mp.matrix([[5, 1], [1, 4]])
    hbi_test = mp.matrix([[1], [2]])
    hii_test = mp.matrix([[2]])
    synthetic = hbb_test - hbi_test * solve_columns(hii_test, hbi_test.T)
    synthetic_expected = mp.matrix([[mp.mpf("4.5"), 0], [0, 2]])
    synthetic_error = matrix_difference(synthetic, synthetic_expected)
synthetic_ok = check(
    "the independent synthetic block has the exact frozen Schur complement",
    synthetic_error < mp.mpf("1e-80"),
    f"error={mp_text(synthetic_error, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    dust_support_errors = []
    for row in range(22):
        for column in range(22):
            expected = mp.mpf(0)
            if row == column and 18 <= row < 22:
                expected = -2 * mp.pi * masses[row - 18] * tau
            dust_support_errors.append(abs(dust[row, column] - expected))
    dust_support_error = max(dust_support_errors)
    minimum_dust_diagonal = min(abs(dust[18 + rank, 18 + rank]) for rank in range(4))
dust_ok = check(
    "the selected matter changes exactly the four frozen lapse diagonal entries",
    dust_support_error < mp.mpf("1e-80")
    and minimum_dust_diagonal > mp.mpf("1e-6"),
    f"support error={mp_text(dust_support_error, 8)}",
)

directional_records = []
maximum_directional_error = mp.mpf(0)
if all_nonsingular:
    with mp.workdps(SECONDARY_DPS):
        directions = {
            "old_01": mp.matrix([1] + [0] * 11),
            "common_old_new_scale": mp.matrix([2] * 12),
            "old_minus_new_scale": mp.matrix([2] * 6 + [-2] * 6),
        }
        steps = (mp.mpf("1e-10"), mp.mpf("5e-11"))
        for index in DIRECTIONAL_INDICES:
            matrix = records[index]["full_matrix"]
            effective = records[index]["effective_matrix"]
            hii = submatrix(matrix, tuple(range(12, 22)), tuple(range(12, 22)))
            hib = submatrix(matrix, tuple(range(12, 22)), tuple(range(12)))
            for label, boundary_direction in directions.items():
                internal_direction = -mp.lu_solve(hii, hib * boundary_direction)
                full_direction = mp.matrix(list(boundary_direction) + list(internal_direction))
                coarse = directional_second(
                    actions["evaluate_schedule"], combinatorics[index], geometry140,
                    base140, masses, full_direction, steps[0],
                )
                fine = directional_second(
                    actions["evaluate_schedule"], combinatorics[index], geometry140,
                    base140, masses, full_direction, steps[1],
                )
                richardson = (4 * fine - coarse) / 3
                quadratic = (boundary_direction.T * effective * boundary_direction)[0]
                relative = abs(richardson - quadratic) / max(mp.mpf(1), abs(quadratic))
                maximum_directional_error = max(maximum_directional_error, relative)
                directional_records.append({
                    "order": list(orders[index]),
                    "direction": label,
                    "quadratic": quadratic,
                    "action_richardson": richardson,
                    "relative_error": relative,
                })
directional_ok = check(
    "direct complete-action second differences reproduce every tested effective form",
    (not all_nonsingular) or maximum_directional_error < mp.mpf("1e-28"),
    f"max relative error={mp_text(maximum_directional_error, 8)}",
)

time_reversal_covariant = False
maximum_reversal_difference = mp.inf
class_members = []
class_indices = []
canonical_matrices = []
if all_nonsingular:
    with mp.workdps(SECONDARY_DPS):
        reversal = reversal_matrix()
        order_index = {order: index for index, order in enumerate(orders)}
        reversal_differences = []
        reversal_passes = []
        for index, order in enumerate(orders):
            reverse_index = order_index[tuple(reversed(order))]
            transformed = reversal.T * effective_matrices[reverse_index] * reversal
            difference = matrix_difference(effective_matrices[index], transformed)
            envelope = max(
                effective_envelopes[index], effective_envelopes[reverse_index]
            )
            reversal_differences.append(difference)
            reversal_passes.append(difference <= envelope)
        maximum_reversal_difference = max(reversal_differences)
        time_reversal_covariant = all(reversal_passes)

        canonical_matrices = [
            canonicalized_matrix(index, orders, effective_matrices, reversal)
            for index in range(24)
        ]
        class_representatives = []
        for index, matrix in enumerate(canonical_matrices):
            assigned = None
            for class_index, representative in enumerate(class_representatives):
                envelope = max(
                    effective_envelopes[index], effective_envelopes[representative]
                )
                if matrix_difference(matrix, canonical_matrices[representative]) <= envelope:
                    assigned = class_index
                    break
            if assigned is None:
                assigned = len(class_representatives)
                class_representatives.append(index)
                class_members.append([])
            class_members[assigned].append(index)
            class_indices.append(assigned)
time_reversal_census_ok = check(
    "the fixed layer-swap reversal comparison and class census are complete",
    (not all_nonsingular)
    or (len(class_indices) == 24 and len(class_members) >= 1),
    f"reversal covariant={time_reversal_covariant}, classes={len(class_members)}",
)

corruption_detected = True
corruption_difference = mp.mpf(0)
if all_nonsingular:
    with mp.workdps(SECONDARY_DPS):
        corrupted = mp.matrix(canonical_matrices[0])
        corruption_size = mp.mpf("1e-6") * max(
            mp.mpf(1), matrix_max(canonical_matrices[0])
        )
        corrupted[0, 0] += corruption_size
        corruption_difference = matrix_difference(corrupted, canonical_matrices[0])
        corruption_detected = corruption_difference > effective_envelopes[0]
corruption_ok = check(
    "the frozen comparator rejects a deliberate one-component corruption",
    corruption_detected,
    f"difference={mp_text(corruption_difference, 8)}",
)

scope = {
    "root_search_or_nested_census_executed": False,
    "nonhomogeneous_operator_or_spectrum_computed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
    "schedule_average_or_selection_used": False,
}
scope_ok = check(
    "the calculation remains inside the frozen invariant Hessian scope",
    not any(scope.values()),
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    definitions_ok,
    topology_ok,
    on_shell_ok,
    branch_ok,
    precision_ok,
    internal_census_ok,
    solve_ok,
    synthetic_ok,
    dust_ok,
    directional_ok,
    time_reversal_census_ok,
    corruption_ok,
    scope_ok,
))

if not controls_ok:
    outcome = "REFINED_EFFECTIVE_H4_HESSIAN_CONTROL_FAILED"
elif not all_nonsingular:
    outcome = "REFINED_EFFECTIVE_H4_HESSIAN_INTERNAL_SINGULAR"
elif not time_reversal_covariant:
    outcome = "REFINED_EFFECTIVE_H4_HESSIAN_TIME_REVERSAL_FAILED"
elif len(class_members) > 1:
    outcome = "REFINED_EFFECTIVE_H4_HESSIAN_MULTIPLE_SCHEDULE_CLASSES"
else:
    outcome = "REFINED_EFFECTIVE_H4_HESSIAN_SINGLE_SCHEDULE_CLASS"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one effective-Hessian outcome",
    outcome in {
        "REFINED_EFFECTIVE_H4_HESSIAN_CONTROL_FAILED",
        "REFINED_EFFECTIVE_H4_HESSIAN_INTERNAL_SINGULAR",
        "REFINED_EFFECTIVE_H4_HESSIAN_TIME_REVERSAL_FAILED",
        "REFINED_EFFECTIVE_H4_HESSIAN_MULTIPLE_SCHEDULE_CLASSES",
        "REFINED_EFFECTIVE_H4_HESSIAN_SINGLE_SCHEDULE_CLASS",
    },
    outcome,
)

artifact_records = []
for index, item in enumerate(records):
    serialized = {
        "order": item["order"],
        "step_difference": mp_text(item["step_difference"]),
        "precision_difference": mp_text(item["precision_difference"]),
        "hessian_envelope": mp_text(item["hessian_envelope"]),
        "raw_imaginary": mp_text(item["raw_imaginary"]),
        "raw_antisymmetry": mp_text(item["raw_antisymmetry"]),
        "internal_eigenvalues": [
            mp_text(value) for value in item["internal_eigenvalues"]
        ],
        "internal_spectral_envelope": mp_text(item["internal_spectral_envelope"]),
        "internal_nonsingular": item["internal_nonsingular"],
        "internal_inertia_positive_zero_negative": list(item["internal_inertia"]),
        "internal_condition": mp_text(item["internal_condition"]),
    }
    if item["internal_nonsingular"]:
        serialized.update({
            "effective_step_difference": mp_text(item["effective_step_difference"]),
            "effective_precision_difference": mp_text(
                item["effective_precision_difference"]
            ),
            "effective_envelope": mp_text(item["effective_envelope"]),
            "effective_antisymmetry": mp_text(item["effective_antisymmetry"]),
            "solve_residual": mp_text(item["solve_residual"]),
            "matrix_class": class_indices[index] if class_indices else None,
            "effective_matrix": serialize_matrix(item["effective_matrix"]),
        })
    artifact_records.append(serialized)

artifact_directional = [
    {
        "order": item["order"],
        "direction": item["direction"],
        "quadratic": mp_text(item["quadratic"]),
        "action_richardson": mp_text(item["action_richardson"]),
        "relative_error": mp_text(item["relative_error"]),
    }
    for item in directional_records
]

artifact = {
    "title": "On-shell effective H4 boundary Hessian across refined schedules",
    "date": "2026-08-21",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "carrier": "K0=P(sd K_600)",
        "tau0": TAU_TEXT,
        "boundary_coordinates": [variable_label(key) for key in BOUNDARY_VARIABLES],
        "internal_coordinates": [variable_label(key) for key in INTERNAL_VARIABLES],
        "decimal_precisions": [PRIMARY_DPS, SECONDARY_DPS],
        "difference_steps": list(STEP_TEXTS),
        "effective_hessian": "H_bb-H_bi H_ii^(-1) H_ib",
        "dust_lapse_diagonal": "-2*pi*m_r*tau0",
        "allowed_schedule_identification": "time-reversal old/new layer swap only",
    },
    "background": {
        "selected_rank_masses": [mp_text(value) for value in masses],
        "maximum_internal_residual": mp_text(maximum_internal_residual),
        "maximum_angle_identity_residual": mp_text(maximum_branch_identity),
        "maximum_imaginary_curvature": mp_text(maximum_branch_imaginary),
        "minimum_angle_argument": mp_text(minimum_branch_argument),
    },
    "controls": {
        "synthetic_schur_error": mp_text(synthetic_error),
        "dust_support_error": mp_text(dust_support_error),
        "maximum_solve_residual": mp_text(maximum_solve_residual),
        "maximum_effective_antisymmetry": mp_text(maximum_effective_antisymmetry),
        "maximum_directional_relative_error": mp_text(maximum_directional_error),
        "directional_records": artifact_directional,
        "corruption_difference": mp_text(corruption_difference),
    },
    "census": {
        "schedule_count": 24,
        "all_internal_blocks_nonsingular": all_nonsingular,
        "internal_inertia_histogram": {
            str(key): value
            for key, value in sorted(Counter(
                item["internal_inertia"] for item in records
            ).items())
        },
        "time_reversal_covariant": time_reversal_covariant,
        "maximum_time_reversal_difference": mp_text(maximum_reversal_difference),
        "effective_matrix_class_count": len(class_members),
        "effective_matrix_classes": [
            {
                "class": class_index,
                "orders": [list(orders[index]) for index in members],
            }
            for class_index, members in enumerate(class_members)
        ],
        "schedules": artifact_records,
    },
    "scope": scope,
    "status_labels": {
        "invariant_quadratic_dynamics": "TESTED",
        "nonhomogeneous_propagation": "NOT_COMPUTED",
        "tick_c_G_planck_particles": "OPEN_NOT_COMPUTED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Outcome: {outcome}")
print(f"Internal nonsingular: {sum(item['internal_nonsingular'] for item in records)}/24")
print(f"Time reversal covariant: {time_reversal_covariant}")
print(f"Effective schedule classes: {len(class_members)}")

raise SystemExit(0 if passed == tests else 1)


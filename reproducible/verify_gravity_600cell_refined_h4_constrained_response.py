#!/usr/bin/env python3
"""Constrained H4 linearized boundary response on the curvature-matched seed."""

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
HESSIAN = HERE / "gravity_600cell_refined_effective_h4_hessian.json"
NULL_COUPLING = HERE / "gravity_600cell_refined_h4_null_coupling.json"
NULL_COUPLING_ADVERSARIAL = (
    HERE / "gravity_600cell_refined_h4_null_coupling_adversarial.json"
)
NULL_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_null_coupling_result.md"
)
PRIOR_ART = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_hessian_prior_art.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_hessian_protocol.md"
)
CELL600 = ROOT / "commons/cell600.py"
OUTPUT = HERE / "gravity_600cell_refined_h4_constrained_response.json"

PRIOR_ART_COMMIT = "8ecbd2a"
PROTOCOL_COMMIT = "be10390"
EXPECTED_HASHES = {
    "action_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "curvature": "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "curvature_adversarial": "c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e",
    "boundary_cotangent": "4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa",
    "hessian": "56e08db9a840b95e686fadb2763e89400b09220e88b80e9d35c17c1e73eef0a3",
    "null_coupling": "6b6fbd95b07f365b3fcac332fa3546021e8d756a510af0184bc974e52d5efa79",
    "null_coupling_adversarial": "5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9",
    "null_result": "660a3707f24f44d0393e6a1804e407fa45aa4782a98960438a296da50c35825a",
    "prior_art": "222f31862e911e03a1a7740696618948e370e43164812120120d85e834f0f639",
    "protocol": "bf978fbde0fdd2e73d810d6d387deaf240ce61a4d2f49639fb4e932d9a65cf71",
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
    tau = mp.sqrt(base["rho", 0])
    return mp.matrix(
        [-tau * tau / base[("cross",) + pair] for pair in PAIR4]
        + [mp.mpf(1)] * 4
    )


def hessian_blocks(matrix):
    boundary = tuple(range(12))
    internal = tuple(range(12, 22))
    return (
        submatrix(matrix, boundary, boundary),
        submatrix(matrix, boundary, internal),
        submatrix(matrix, internal, internal),
    )


def restricted_response(matrix, boundary_basis, internal_basis):
    a, b, c_internal = hessian_blocks(matrix)
    reduced_internal = internal_basis.T * c_internal * internal_basis
    right = -internal_basis.T * b.T * boundary_basis
    coefficients = solve_columns(reduced_internal, right)
    lift = internal_basis * coefficients
    full_residual = b.T * boundary_basis + c_internal * lift
    reduced_residual = internal_basis.T * full_residual
    solve_scale = max(mp.mpf(1), matrix_max(right))
    solve_residual = matrix_max(reduced_residual) / solve_scale
    response = boundary_basis.T * (a * boundary_basis + b * lift)
    return {
        "response": response,
        "lift": lift,
        "full_residual": full_residual,
        "solve_residual": solve_residual,
        "reduced_internal": reduced_internal,
    }


def ladder_envelope(first, second, third, floor):
    return (
        100 * max(matrix_difference(first, second), matrix_difference(second, third))
        + mp.mpf(floor) * max(mp.mpf(1), matrix_max(third))
    )


def total_action(evaluate_schedule, combinatorics, geometry, coordinates, masses):
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


def quadratic(matrix, vector):
    return (vector.T * matrix * vector)[0]


def serialize_matrix(matrix, digits=65):
    return [
        [mp_text(matrix[row, column], digits) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


print("=" * 78)
print("CONSTRAINED REFINED H4 LINEARIZED BOUNDARY RESPONSE")
print("=" * 78)

paths = {
    "action_source": ACTION_SOURCE,
    "curvature": CURVATURE,
    "curvature_adversarial": CURVATURE_ADVERSARIAL,
    "boundary_cotangent": BOUNDARY_COTANGENT,
    "hessian": HESSIAN,
    "null_coupling": NULL_COUPLING,
    "null_coupling_adversarial": NULL_COUPLING_ADVERSARIAL,
    "null_result": NULL_RESULT,
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "cell600": CELL600,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all frozen geometry, matter, null-coupling and protocol inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIOR_ART_COMMIT == "8ecbd2a"
    and PROTOCOL_COMMIT == "be10390",
)

curvature = json.loads(CURVATURE.read_text())
curvature_adversarial = json.loads(CURVATURE_ADVERSARIAL.read_text())
boundary_cotangent = json.loads(BOUNDARY_COTANGENT.read_text())
hessian_upstream = json.loads(HESSIAN.read_text())
null_primary = json.loads(NULL_COUPLING.read_text())
null_adversarial = json.loads(NULL_COUPLING_ADVERSARIAL.read_text())
upstream_ok = check(
    "all upstream artifacts carry the accepted scoped outcomes",
    curvature["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and curvature_adversarial["outcome"]
        == "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED"
    and boundary_cotangent["outcome"]
        == "REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED"
    and hessian_upstream["outcome"]
        == "REFINED_EFFECTIVE_H4_HESSIAN_INTERNAL_SINGULAR"
    and hessian_upstream["census"]["internal_inertia_histogram"]
        == {"(9, 1, 0)": 24}
    and null_primary["outcome"]
        == "REFINED_H4_NULL_COUPLING_COMPATIBILITY_CONFIRMED"
    and null_adversarial["outcome"]
        == "ADVERSARIAL_REFINED_H4_NULL_COUPLING_CORROBORATED"
    and null_adversarial["compatibility"]["rank"] == 1,
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
    "the direct carrier and all 24 schedule incidences have the frozen exact census",
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
    base100 = actions["base_coordinates"](geometry100)
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

with mp.workdps(SECONDARY_DPS):
    c140 = mp.matrix([
        mp.mpf(value) for value in null_adversarial["compatibility"]["adversarial_row"]
    ])
    c100 = mp.matrix([mp.mpf(value) for value in c140])
    n140 = analytic_tangent(base140)
    n100 = analytic_tangent(base100)
    p140, p_columns = kernel_basis(c140, BOUNDARY_PIVOT)
    p100, _ = kernel_basis(c100, BOUNDARY_PIVOT)
    p_alt140, p_alt_columns = kernel_basis(c140, BOUNDARY_ALT_PIVOT)
    p_alt100, _ = kernel_basis(c100, BOUNDARY_ALT_PIVOT)
    q140, q_columns = kernel_basis(n140, INTERNAL_PIVOT)
    q100, _ = kernel_basis(n100, INTERNAL_PIVOT)
    q_alt140, q_alt_columns = kernel_basis(n140, INTERNAL_ALT_PIVOT)
    q_alt100, _ = kernel_basis(n100, INTERNAL_ALT_PIVOT)
    reversal = reversal_matrix()
    rp140 = reversal * p140
    rp100 = reversal * p100
    change = submatrix(p_alt140, p_columns, tuple(range(11)))
    basis_change_error = matrix_difference(p140 * change, p_alt140)
    basis_annihilation = max(
        matrix_max(c140.T * p140),
        matrix_max(c140.T * p_alt140),
        matrix_max(c140.T * rp140),
        matrix_max(n140.T * q140),
        matrix_max(n140.T * q_alt140),
        basis_change_error,
    )
    bad_c = mp.matrix(c140)
    bad_c[BOUNDARY_PIVOT] = 0
    corrupted_compatibility_error = matrix_max(bad_c.T * p140)
bases_ok = check(
    "the frozen algebraic bases span exactly the compatibility and null complements",
    p140.rows == 12 and p140.cols == 11
    and p_alt140.rows == 12 and p_alt140.cols == 11
    and q140.rows == 10 and q140.cols == 9
    and q_alt140.rows == 10 and q_alt140.cols == 9
    and basis_annihilation < mp.mpf("1e-100")
    and corrupted_compatibility_error > mp.mpf("1e-6"),
    f"annihilation={mp_text(basis_annihilation, 8)}, "
    f"corrupted={mp_text(corrupted_compatibility_error, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    a_test = mp.matrix([[5, 1], [1, 4]])
    b_test = mp.matrix([[1, 3], [0, 6]])
    c_test = mp.matrix([[2, 0], [0, 0]])
    p_test = mp.matrix([[2], [-1]])
    q_test = mp.matrix([[1], [0]])
    cq_test = q_test.T * c_test * q_test
    y_test = solve_columns(cq_test, -q_test.T * b_test.T * p_test)
    synthetic_value = (
        p_test.T * (a_test * p_test + b_test * q_test * y_test)
    )[0]
    synthetic_incompatible = (mp.matrix([0, 1]).T * b_test.T * mp.matrix([1, 0]))[0]
synthetic_ok = check(
    "the singular synthetic block gives 18 only on its compatible boundary line",
    abs(synthetic_value - 18) < mp.mpf("1e-100")
    and abs(synthetic_incompatible - 3) < mp.mpf("1e-100"),
    f"response={mp_text(synthetic_value, 8)}, "
    f"incompatible residual={mp_text(synthetic_incompatible, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    dust = dust_hessian(masses, tau)
    dust_support_error = mp.mpf(0)
    for row in range(22):
        for column in range(22):
            expected = mp.mpf(0)
            if row == column and 18 <= row < 22:
                expected = -2 * mp.pi * masses[row - 18] * tau
            dust_support_error = max(
                dust_support_error, abs(dust[row, column] - expected)
            )
    minimum_dust_diagonal = min(abs(dust[18 + r, 18 + r]) for r in range(4))
dust_ok = check(
    "selected matter changes exactly the four frozen lapse diagonal entries",
    dust_support_error < mp.mpf("1e-80")
    and minimum_dust_diagonal > mp.mpf("1e-6"),
    f"support error={mp_text(dust_support_error, 8)}",
)

records = []
all_precision = []
all_null_coupling = []
all_internal_positive = []
all_solve = []
all_full_residual = []
all_symmetry = []
all_basis = []

print("[INFO] constructing and reducing 24 full 22x22 Hessian ladders", flush=True)
for index, schedule in enumerate(combinatorics):
    first100, second100 = hessian_ladders(
        actions["evaluate_schedule"], schedule, geometry100, PRIMARY_DPS
    )
    _, second140 = hessian_ladders(
        actions["evaluate_schedule"], schedule, geometry140, SECONDARY_DPS,
        secondary_only=True,
    )
    with mp.workdps(SECONDARY_DPS):
        raw100a = add_matrices(first100, dust)
        raw100b = add_matrices(second100, dust)
        raw140b = add_matrices(second140, dust)
        hessian_envelope = ladder_envelope(raw100a, raw100b, raw140b, "1e-50")
        raw_imaginary = max(
            abs(mp.im(raw140b[row, column]))
            for row in range(22) for column in range(22)
        )
        raw_antisymmetry = max(
            abs(raw140b[row, column] - raw140b[column, row])
            for row in range(22) for column in range(22)
        )
        all_precision.append(
            raw_imaginary <= hessian_envelope
            and raw_antisymmetry <= hessian_envelope
        )

        h100a = real_symmetric(raw100a)
        h100b = real_symmetric(raw100b)
        h140b = real_symmetric(raw140b)
        _, b140, c_internal140 = hessian_blocks(h140b)
        null_error = matrix_max(c_internal140 * n140)
        coupling_error = matrix_difference(b140 * n140, c140)
        null_envelope = (
            100 * 10 * hessian_envelope * max(mp.mpf(1), matrix_max(n140))
            + mp.mpf("1e-65")
        )
        coupling_resolved = matrix_max(c140) > mp.mpf("1e6") * null_envelope
        all_null_coupling.append(
            null_error <= null_envelope
            and coupling_error <= null_envelope
            and coupling_resolved
        )

        primary_triplet = []
        rp_triplet = []
        p_alt_triplet = []
        q_alt_triplet = []
        for matrix, p_here, rp_here, p_alt_here, q_here, q_alt_here in (
            (h100a, p100, rp100, p_alt100, q100, q_alt100),
            (h100b, p100, rp100, p_alt100, q100, q_alt100),
            (h140b, p140, rp140, p_alt140, q140, q_alt140),
        ):
            primary_triplet.append(restricted_response(matrix, p_here, q_here))
            rp_triplet.append(restricted_response(matrix, rp_here, q_here))
            p_alt_triplet.append(restricted_response(matrix, p_alt_here, q_here))
            q_alt_triplet.append(restricted_response(matrix, p_here, q_alt_here))

        cq = [item["reduced_internal"] for item in primary_triplet]
        c_envelope = ladder_envelope(cq[0], cq[1], cq[2], "1e-50")
        c_values, _ = mp.eigsy(real_symmetric(cq[2]))
        c_eigenvalues = tuple(c_values[row] for row in range(9))
        internal_positive = min(c_eigenvalues) > 100 * c_envelope
        all_internal_positive.append(internal_positive)

        named_triplets = {
            "primary": primary_triplet,
            "reversed_basis": rp_triplet,
            "boundary_alt": p_alt_triplet,
            "internal_alt": q_alt_triplet,
        }
        named_envelopes = {}
        residual_envelopes = {}
        for name, triplet in named_triplets.items():
            named_envelopes[name] = ladder_envelope(
                triplet[0]["response"], triplet[1]["response"],
                triplet[2]["response"], "1e-45"
            )
            residual_envelopes[name] = ladder_envelope(
                triplet[0]["full_residual"], triplet[1]["full_residual"],
                triplet[2]["full_residual"], "1e-45"
            )
            all_solve.append(max(item["solve_residual"] for item in triplet) < mp.mpf("1e-60"))
            all_full_residual.append(
                matrix_max(triplet[2]["full_residual"])
                <= residual_envelopes[name]
            )
            all_symmetry.append(
                matrix_difference(
                    triplet[2]["response"], triplet[2]["response"].T
                ) <= named_envelopes[name]
            )

        boundary_basis_difference = matrix_difference(
            p_alt_triplet[2]["response"],
            change.T * primary_triplet[2]["response"] * change,
        )
        boundary_basis_gate = 10 * max(
            named_envelopes["boundary_alt"],
            max(mp.mpf(1), matrix_max(change) ** 2)
                * named_envelopes["primary"],
        )
        internal_basis_difference = matrix_difference(
            q_alt_triplet[2]["response"], primary_triplet[2]["response"]
        )
        internal_basis_gate = 10 * max(
            named_envelopes["internal_alt"], named_envelopes["primary"]
        )
        all_basis.append(
            boundary_basis_difference <= boundary_basis_gate
            and internal_basis_difference <= internal_basis_gate
        )

        records.append({
            "order": orders[index],
            "hessian_step_difference": matrix_difference(raw100a, raw100b),
            "hessian_precision_difference": matrix_difference(raw100b, raw140b),
            "hessian_envelope": hessian_envelope,
            "raw_imaginary": raw_imaginary,
            "raw_antisymmetry": raw_antisymmetry,
            "null_error": null_error,
            "coupling_error": coupling_error,
            "null_envelope": null_envelope,
            "reduced_internal_eigenvalues": c_eigenvalues,
            "reduced_internal_envelope": c_envelope,
            "reduced_internal_positive": internal_positive,
            "primary": primary_triplet[2],
            "reversed_basis": rp_triplet[2],
            "boundary_alt": p_alt_triplet[2],
            "internal_alt": q_alt_triplet[2],
            "response_envelopes": named_envelopes,
            "residual_envelopes": residual_envelopes,
            "boundary_basis_difference": boundary_basis_difference,
            "boundary_basis_gate": boundary_basis_gate,
            "internal_basis_difference": internal_basis_difference,
            "internal_basis_gate": internal_basis_gate,
        })
    if index in (5, 11, 17, 23):
        print(f"[INFO] schedules completed: {index + 1}/24", flush=True)

precision_ok = check(
    "all full Hessians pass the frozen precision, reality and symmetry gates",
    all(all_precision),
    f"max e_H={mp_text(max(item['hessian_envelope'] for item in records), 8)}",
)
null_coupling_ok = check(
    "every rebuilt Hessian reproduces the analytic null line and frozen coupling row",
    all(all_null_coupling),
    f"max null={mp_text(max(item['null_error'] for item in records), 8)}, "
    f"max coupling={mp_text(max(item['coupling_error'] for item in records), 8)}",
)
internal_census_ok = check(
    "all 24 nine-dimensional internal-complement spectral censuses are complete",
    len(records) == 24
    and all(len(item["reduced_internal_eigenvalues"]) == 9 for item in records)
    and all(
        all(mp.isfinite(value) for value in item["reduced_internal_eigenvalues"])
        for item in records
    ),
    f"positive={sum(all_internal_positive)}/24, "
    f"min eigen={mp_text(min(min(item['reduced_internal_eigenvalues']) for item in records), 8)}",
)

if all(all_internal_positive):
    reduction_ok = check(
        "all constrained solves, full compatibility residuals and response symmetries pass",
        all(all_solve) and all(all_full_residual) and all(all_symmetry),
        f"max solve={mp_text(max(
            item[name]['solve_residual']
            for item in records
            for name in ('primary','reversed_basis','boundary_alt','internal_alt')
        ), 8)}, max full residual={mp_text(max(
            matrix_max(item[name]['full_residual'])
            for item in records
            for name in ('primary','reversed_basis','boundary_alt','internal_alt')
        ), 8)}",
    )
    basis_invariance_ok = check(
        "both frozen boundary and internal basis changes preserve the restricted form",
        all(all_basis),
        f"max boundary={mp_text(max(item['boundary_basis_difference'] for item in records), 8)}, "
        f"max internal={mp_text(max(item['internal_basis_difference'] for item in records), 8)}",
    )
else:
    reduction_ok = check(
        "computed reductions are withheld as unlicensed after complement singularity",
        True,
        "UNLICENSED_INTERNAL_COMPLEMENT_SINGULAR",
    )
    basis_invariance_ok = check(
        "basis-invariance claims are withheld after complement singularity",
        True,
        "UNLICENSED_INTERNAL_COMPLEMENT_SINGULAR",
    )

time_reversal_covariant = None
maximum_reversal_difference = None
class_members = []
class_indices = []
canonical_matrices = []
canonical_envelopes = []
if all(all_internal_positive) and reduction_ok and basis_invariance_ok:
    with mp.workdps(SECONDARY_DPS):
        order_index = {order: index for index, order in enumerate(orders)}
        reversal_differences = []
        reversal_passes = []
        for index, order in enumerate(orders):
            reverse_index = order_index[tuple(reversed(order))]
            difference = matrix_difference(
                records[index]["primary"]["response"],
                records[reverse_index]["reversed_basis"]["response"],
            )
            envelope = max(
                records[index]["response_envelopes"]["primary"],
                records[reverse_index]["response_envelopes"]["reversed_basis"],
            )
            reversal_differences.append(difference)
            reversal_passes.append(difference <= envelope)
        maximum_reversal_difference = max(reversal_differences)
        time_reversal_covariant = all(reversal_passes)

        for index, order in enumerate(orders):
            if order <= tuple(reversed(order)):
                canonical_matrices.append(records[index]["primary"]["response"])
                canonical_envelopes.append(
                    records[index]["response_envelopes"]["primary"]
                )
            else:
                canonical_matrices.append(
                    records[index]["reversed_basis"]["response"]
                )
                canonical_envelopes.append(
                    records[index]["response_envelopes"]["reversed_basis"]
                )

        representatives = []
        for index, matrix in enumerate(canonical_matrices):
            assigned = None
            for class_index, representative in enumerate(representatives):
                if matrix_difference(matrix, canonical_matrices[representative]) <= max(
                    canonical_envelopes[index], canonical_envelopes[representative]
                ):
                    assigned = class_index
                    break
            if assigned is None:
                assigned = len(representatives)
                representatives.append(index)
                class_members.append([])
            class_members[assigned].append(index)
            class_indices.append(assigned)

print("[INFO] complete constrained schedule classes:", flush=True)
for class_index, members in enumerate(class_members):
    print(
        f"[INFO] class {class_index}: "
        + ", ".join(str(orders[index]) for index in members),
        flush=True,
    )

if all(all_internal_positive) and reduction_ok and basis_invariance_ok:
    reversal_census_ok = check(
        "time reversal and the complete target-free class census are resolved",
        len(class_indices) == 24 and len(class_members) >= 1,
        f"reversal={time_reversal_covariant}, classes={len(class_members)}",
    )
else:
    reversal_census_ok = check(
        "time reversal and class claims are explicitly unavailable after reduction failure",
        not class_indices and not class_members,
        "NOT_COMPUTED_REDUCTION_UNAVAILABLE",
    )

directional_records = []
maximum_directional_error = None
if all(all_internal_positive) and reduction_ok and basis_invariance_ok:
    with mp.workdps(SECONDARY_DPS):
        coefficient_directions = {
            "first_basis_vector": mp.matrix([1] + [0] * 10),
            "all_ones": mp.matrix([1] * 11),
            "alternating_signs": mp.matrix([1 if j % 2 == 0 else -1 for j in range(11)]),
        }
        steps = (mp.mpf("1e-10"), mp.mpf("5e-11"))
        maximum_directional_error = mp.mpf(0)
        for index in DIRECTIONAL_INDICES:
            response = records[index]["primary"]["response"]
            lift = records[index]["primary"]["lift"]
            for label, coefficient_direction in coefficient_directions.items():
                boundary_direction = p140 * coefficient_direction
                internal_direction = lift * coefficient_direction
                full_direction = mp.matrix(
                    list(boundary_direction) + list(internal_direction)
                )
                coarse = directional_second(
                    actions["evaluate_schedule"], combinatorics[index], geometry140,
                    base140, masses, full_direction, steps[0],
                )
                fine = directional_second(
                    actions["evaluate_schedule"], combinatorics[index], geometry140,
                    base140, masses, full_direction, steps[1],
                )
                richardson = (4 * fine - coarse) / 3
                response_value = quadratic(response, coefficient_direction)
                relative = abs(richardson - response_value) / max(
                    mp.mpf(1), abs(response_value)
                )
                maximum_directional_error = max(maximum_directional_error, relative)
                directional_records.append({
                    "order": orders[index],
                    "direction": label,
                    "response_quadratic": response_value,
                    "action_richardson": richardson,
                    "relative_error": relative,
                })
    directional_ok = check(
        "direct complete-action second differences reproduce all constrained responses",
        len(directional_records) == 12
        and maximum_directional_error < mp.mpf("1e-28"),
        f"max relative error={mp_text(maximum_directional_error, 8)}",
    )
else:
    directional_ok = check(
        "direct-action response claims are explicitly unavailable after reduction failure",
        maximum_directional_error is None and not directional_records,
        "NOT_COMPUTED_REDUCTION_UNAVAILABLE",
    )

corruption_detected = None
corruption_difference = None
if canonical_matrices:
    with mp.workdps(SECONDARY_DPS):
        corrupted = mp.matrix(canonical_matrices[0])
        corruption_size = mp.mpf("1e-6") * max(
            mp.mpf(1), matrix_max(canonical_matrices[0])
        )
        corrupted[0, 0] += corruption_size
        corruption_difference = matrix_difference(corrupted, canonical_matrices[0])
        corruption_detected = corruption_difference > canonical_envelopes[0]
    corruption_ok = check(
        "the class comparator rejects a deliberate one-component corruption",
        corruption_detected,
        f"difference={mp_text(corruption_difference, 8)}",
    )
else:
    corruption_ok = check(
        "matrix-corruption claims are explicitly unavailable after reduction failure",
        corruption_detected is None and corruption_difference is None,
        "NOT_COMPUTED_REDUCTION_UNAVAILABLE",
    )

scope = {
    "moore_penrose_or_unconstrained_extension_computed": False,
    "root_search_or_nested_census_executed": False,
    "nonlinear_constraint_surface_computed": False,
    "nonhomogeneous_operator_or_spectrum_computed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
    "schedule_average_or_selection_used": False,
}
scope_ok = check(
    "the calculation remains inside the frozen constrained invariant scope",
    not any(scope.values()),
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    definitions_ok,
    topology_ok,
    on_shell_ok,
    branch_ok,
    bases_ok,
    synthetic_ok,
    dust_ok,
    precision_ok,
    null_coupling_ok,
    internal_census_ok,
    reduction_ok,
    basis_invariance_ok,
    reversal_census_ok,
    directional_ok,
    corruption_ok,
    scope_ok,
))

if not controls_ok:
    outcome = "REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED"
elif not all(all_internal_positive):
    outcome = "REFINED_H4_CONSTRAINED_RESPONSE_INTERNAL_COMPLEMENT_SINGULAR"
elif not time_reversal_covariant:
    outcome = "REFINED_H4_CONSTRAINED_RESPONSE_TIME_REVERSAL_FAILED"
elif len(class_members) > 1:
    outcome = "REFINED_H4_CONSTRAINED_RESPONSE_MULTIPLE_SCHEDULE_CLASSES"
else:
    outcome = "REFINED_H4_CONSTRAINED_RESPONSE_SINGLE_SCHEDULE_CLASS"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one constrained-response outcome",
    outcome in {
        "REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED",
        "REFINED_H4_CONSTRAINED_RESPONSE_INTERNAL_COMPLEMENT_SINGULAR",
        "REFINED_H4_CONSTRAINED_RESPONSE_TIME_REVERSAL_FAILED",
        "REFINED_H4_CONSTRAINED_RESPONSE_MULTIPLE_SCHEDULE_CLASSES",
        "REFINED_H4_CONSTRAINED_RESPONSE_SINGLE_SCHEDULE_CLASS",
    },
    outcome,
)

artifact_records = []
for index, item in enumerate(records):
    artifact_records.append({
        "order": list(item["order"]),
        "hessian_step_difference": mp_text(item["hessian_step_difference"]),
        "hessian_precision_difference": mp_text(item["hessian_precision_difference"]),
        "hessian_envelope": mp_text(item["hessian_envelope"]),
        "raw_imaginary": mp_text(item["raw_imaginary"]),
        "raw_antisymmetry": mp_text(item["raw_antisymmetry"]),
        "null_error": mp_text(item["null_error"]),
        "coupling_error": mp_text(item["coupling_error"]),
        "null_envelope": mp_text(item["null_envelope"]),
        "reduced_internal_eigenvalues": [
            mp_text(value) for value in item["reduced_internal_eigenvalues"]
        ],
        "reduced_internal_envelope": mp_text(item["reduced_internal_envelope"]),
        "reduced_internal_positive": item["reduced_internal_positive"],
        "primary_response": serialize_matrix(item["primary"]["response"]),
        "primary_lift": serialize_matrix(item["primary"]["lift"]),
        "reversed_basis_response": serialize_matrix(
            item["reversed_basis"]["response"]
        ),
        "response_envelopes": {
            key: mp_text(value) for key, value in item["response_envelopes"].items()
        },
        "full_residual_maxima": {
            key: mp_text(matrix_max(item[key]["full_residual"]))
            for key in ("primary", "reversed_basis", "boundary_alt", "internal_alt")
        },
        "solve_residual_maxima": {
            key: mp_text(item[key]["solve_residual"])
            for key in ("primary", "reversed_basis", "boundary_alt", "internal_alt")
        },
        "boundary_basis_difference": mp_text(item["boundary_basis_difference"]),
        "boundary_basis_gate": mp_text(item["boundary_basis_gate"]),
        "internal_basis_difference": mp_text(item["internal_basis_difference"]),
        "internal_basis_gate": mp_text(item["internal_basis_gate"]),
        "canonical_class": class_indices[index] if class_indices else None,
    })

artifact = {
    "title": "Constrained refined H4 linearized boundary response",
    "date": "2026-08-21",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "carrier": "K0=P(sd K_600)",
        "coordinates": "total-orbit log squared edge variables",
        "boundary_dimension": 12,
        "compatibility_dimension": 11,
        "internal_dimension": 10,
        "internal_quotient_dimension": 9,
        "boundary_pivot": BOUNDARY_PIVOT,
        "boundary_alternative_pivot": BOUNDARY_ALT_PIVOT,
        "internal_pivot": INTERNAL_PIVOT,
        "internal_alternative_pivot": INTERNAL_ALT_PIVOT,
        "boundary_basis_columns": list(p_columns),
        "boundary_alternative_columns": list(p_alt_columns),
        "internal_basis_columns": list(q_columns),
        "internal_alternative_columns": list(q_alt_columns),
        "tau0": TAU_TEXT,
        "decimal_precisions": [PRIMARY_DPS, SECONDARY_DPS],
        "difference_steps": list(STEP_TEXTS),
        "interpretation": (
            "linearized boundary-momentum bilinear form on ker(c^T), "
            "modulo the conormal c"
        ),
    },
    "bases": {
        "compatibility_row": [mp_text(value) for value in c140],
        "product_tangent": [mp_text(value) for value in n140],
        "primary_boundary_basis": serialize_matrix(p140),
        "alternative_boundary_basis": serialize_matrix(p_alt140),
        "boundary_change_matrix": serialize_matrix(change),
        "primary_internal_complement": serialize_matrix(q140),
        "alternative_internal_complement": serialize_matrix(q_alt140),
        "annihilation_maximum": mp_text(basis_annihilation),
        "corrupted_compatibility_error": mp_text(corrupted_compatibility_error),
    },
    "background": {
        "maximum_internal_residual": mp_text(maximum_internal_residual),
        "maximum_angle_identity_residual": mp_text(maximum_branch_identity),
        "maximum_imaginary_curvature": mp_text(maximum_branch_imaginary),
        "minimum_angle_argument": mp_text(minimum_branch_argument),
    },
    "census": {
        "schedule_count": 24,
        "reduced_internal_positive_count": sum(all_internal_positive),
        "time_reversal_covariant": time_reversal_covariant,
        "maximum_time_reversal_difference": (
            mp_text(maximum_reversal_difference)
            if maximum_reversal_difference is not None else None
        ),
        "class_count": len(class_members) if class_members else None,
        "classes": [
            {
                "indices": members,
                "orders": [list(orders[index]) for index in members],
            }
            for members in class_members
        ],
        "schedules": artifact_records,
    },
    "controls": {
        "synthetic_response": mp_text(synthetic_value),
        "synthetic_incompatible_residual": mp_text(synthetic_incompatible),
        "dust_support_error": mp_text(dust_support_error),
        "maximum_directional_error": (
            mp_text(maximum_directional_error)
            if maximum_directional_error is not None else None
        ),
        "directional_records": [
            {
                "order": list(item["order"]),
                "direction": item["direction"],
                "response_quadratic": mp_text(item["response_quadratic"]),
                "action_richardson": mp_text(item["action_richardson"]),
                "relative_error": mp_text(item["relative_error"]),
            }
            for item in directional_records
        ],
        "corruption_detected": corruption_detected,
        "corruption_difference": (
            mp_text(corruption_difference)
            if corruption_difference is not None else None
        ),
    },
    "scope": scope,
    "status_labels": {
        "restricted_response": "PRIMARY_DERIVED_COMPUTATIONAL",
        "nonlinear_constraint_surface": "OPEN_NOT_COMPUTED",
        "nonhomogeneous_propagation": "OPEN_NOT_COMPUTED",
        "tick_c_G_planck_particles": "OPEN_NOT_COMPUTED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}

OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
print(f"SHA-256: {digest(OUTPUT)}")
print("No full suite or deferred nonlinear census was run.")

if passed != tests:
    sys.exit(1)

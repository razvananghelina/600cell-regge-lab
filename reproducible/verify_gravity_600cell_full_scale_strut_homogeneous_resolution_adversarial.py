#!/usr/bin/env python3
"""Direct-matrix adversarial replication of the homogeneous weak-pole line."""

import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

from flint import acb, acb_mat, ctx
import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial_protocol.md"
PRIMARY_RESULT = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_result.md"
PRIMARY_RESOLUTION = HERE / "gravity_600cell_full_scale_strut_homogeneous_resolution.json"
ACTION_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
ACTION_INPUT = HERE / "gravity_600cell_dust_hyperbolic_lapse_alignment.json"
HOMOGENEOUS_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
CARRIER_INPUT = HERE / "gravity_600cell_full_scale_strut_carrier.json"
INTERSECTION_INPUT = HERE / "gravity_600cell_full_scale_strut_canonical_intersection.json"
OUTPUT = HERE / "gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial.json"

PROTOCOL_COMMIT = "1cca153"
PRIMARY_RESULT_COMMIT = "3ee5c55"
EXPECTED_HASHES = {
    "protocol": "4b5fd36881c3d67e795599d7e00c44c9fd69a680a8ec84ec98a9e97e8c8dbff1",
    "primary_result": "0d56d6b0122f8a059bbbdaf4f32623cf1b5ae05461947467c1f18bd25e26f9fb",
    "primary_resolution": "70d7583756acdbee77893f98d57054ab074d9353a86247840cc1eb2c7b6be931",
    "action_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "action": "a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff",
    "homogeneous": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "carrier": "6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d",
    "intersection": "b29cc33a9effeb2087fb6133359ee747d100d203778586372a7ceeebc2e4f070",
}
INPUTS = {
    "protocol": PROTOCOL,
    "primary_result": PRIMARY_RESULT,
    "primary_resolution": PRIMARY_RESOLUTION,
    "action_source": ACTION_SOURCE,
    "action": ACTION_INPUT,
    "homogeneous": HOMOGENEOUS_INPUT,
    "carrier": CARRIER_INPUT,
    "intersection": INTERSECTION_INPUT,
}

VERTICES = 120
P160_DPS = 160
P160_BALL_DPS = 140
P160_DIGITS = 145
P160_STEPS = {
    "operational_primary": mp.mpf("1e-40"),
    "operational_shadow": mp.mpf("1e-30"),
    "validation_primary": mp.mpf("3e-40"),
    "validation_shadow": mp.mpf("3e-30"),
}

tests = passed = 0


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
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=70):
    if isinstance(value, mp.mpc):
        return {
            "real": mp.nstr(mp.re(value), digits),
            "imaginary": mp.nstr(mp.im(value), digits),
        }
    return mp.nstr(value, digits)


def arb_mid_to_mp(value):
    mantissa, exponent = value.mid().man_exp()
    return mp.mpf(int(mantissa)) * mp.power(2, int(exponent))


def acb_mid_to_mp(value):
    return mp.mpc(arb_mid_to_mp(value.real), arb_mid_to_mp(value.imag))


def mp_to_acb(value):
    return acb(mp.nstr(mp.re(value), P160_DIGITS), mp.nstr(mp.im(value), P160_DIGITS))


def mp_matrix_to_acb(matrix):
    return acb_mat(
        matrix.rows,
        matrix.cols,
        [mp_to_acb(matrix[row, column])
         for row in range(matrix.rows) for column in range(matrix.cols)],
    )


def acb_matrix_midpoint(matrix):
    return mp.matrix([
        [acb_mid_to_mp(matrix[row, column]) for column in range(matrix.ncols())]
        for row in range(matrix.nrows())
    ])


def acb_subcolumns(matrix, columns):
    return acb_mat(
        matrix.nrows(), len(columns),
        [matrix[row, column] for row in range(matrix.nrows()) for column in columns],
    )


def gram_determinant_record(matrix):
    determinant = (matrix.conjugate().transpose() * matrix).det()
    return {
        "ball": str(determinant),
        "contains_zero": bool(determinant.contains(0)),
        "abs_lower": str(determinant.abs_lower()),
        "abs_upper": str(determinant.abs_upper()),
    }


def vector_norm(vector):
    return mp.sqrt(mp.fsum(abs(vector[index]) ** 2 for index in range(len(vector))))


def matrix_norm(matrix):
    return mp.sqrt(mp.fsum(
        abs(matrix[row, column]) ** 2
        for row in range(matrix.rows) for column in range(matrix.cols)
    ))


def normalized_residual(matrix, vector):
    return vector_norm(matrix * vector) / (matrix_norm(matrix) * vector_norm(vector))


def normalize_phase(vector):
    norm = vector_norm(vector)
    vector = vector / norm
    pivot = max(range(len(vector)), key=lambda index: abs(vector[index]))
    if abs(vector[pivot]):
        vector *= mp.conj(vector[pivot]) / abs(vector[pivot])
    return vector


def projector_distance(left, right):
    left = normalize_phase(left.copy())
    right = normalize_phase(right.copy())
    return mp.sqrt(mp.fsum(
        abs(left[row] * mp.conj(left[column])
            - right[row] * mp.conj(right[column])) ** 2
        for row in range(len(left)) for column in range(len(left))
    ))


def deterministic_normal_null(matrix):
    """Fix the last entry to one and solve normal equations for all earlier entries."""
    dependent = matrix[:, :matrix.cols - 1]
    rhs = -matrix[:, matrix.cols - 1]
    normal = dependent.H * dependent
    solution = mp.lu_solve(normal, dependent.H * rhs)
    vector = mp.matrix(matrix.cols, 1)
    for index in range(matrix.cols - 1):
        vector[index] = solution[index]
    vector[matrix.cols - 1] = 1
    return normalize_phase(vector)


def configure_action_precision(old):
    mp.mp.dps = P160_DPS
    ctx.dps = P160_BALL_DPS
    m_star = mp.mpf(10)
    zeta = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
    r0 = 4 * m_star / (3 * mp.pi)
    l0 = zeta * r0
    l0_square = l0**2
    epsilon3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
    mass = (90 / mp.pi) * epsilon3 * l0
    rho0 = mp.mpf("0.0102") ** 2
    functions = (
        old["group_and_index_data"], old["high_precision_pattern_cache"],
        old["assemble_full_representative_kernels"], old["area_data"],
    )
    for function in functions:
        namespace = function.__globals__
        namespace["L0_SQUARE"] = l0_square
        namespace["MASS"] = mass
        namespace["RHO0"] = rho0
        namespace["DERIVATIVE_STEPS"] = P160_STEPS
        namespace["VARIANTS"] = tuple(P160_STEPS)
        namespace["ARITHMETIC_FLOOR"] = mp.mpf("1e-130")
    return {"L0_square": l0_square, "mass": mass, "rho0": rho0}


def direct_carrier_entries(parity, index_data, weak_positions, carrier, tick, constants):
    l0_square = constants["L0_square"]
    lam = mp.exp(mp.mpf(tick["solutions"][parity]["state"][0]))
    rho = index_data["rho"]
    q_diagonal = lam * l0_square - rho

    # Raw exact carrier responses, independently restated from the frozen
    # scale+strut formula rather than imported from the precision verifier.
    source_scale = (l0_square / (8 * q_diagonal)) * (
        -16 * rho / (l0_square * (lam - 1) ** 2)
    )
    target_scale = (l0_square / (8 * q_diagonal)) * (
        8 + 16 * rho / (l0_square * (lam - 1) ** 2)
    )
    source_strut = rho / ((lam - 1) * q_diagonal)
    target_strut = -lam * source_strut

    orbit_order = tuple(range(30, 65)) + tuple(range(65, 95))
    orbit_position = {orbit: position for position, orbit in enumerate(orbit_order)}
    pole_edges = [
        tuple(map(int, edge))
        for position in weak_positions
        for edge in index_data["orbit_edges"][30 + position]
    ]
    logical_to_column = {edge[0]: column for column, edge in enumerate(pole_edges)}
    frozen = carrier["parities"][parity]
    internal_edges = [tuple(map(int, edge)) for edge in frozen["internal_edge_order"]]
    final_edges = [tuple(map(int, edge)) for edge in frozen["final_edge_order"]]
    entries = []

    def append(row, column, value):
        if value:
            entries.append((row, column, value))

    for edge in internal_edges:
        global_index = index_data["edge_to_index"][edge]
        orbit, group = divmod(global_index, 24)
        row = 24 * orbit_position[orbit] + group
        lower, upper = edge
        if upper == lower + VERTICES:
            append(row, VERTICES + logical_to_column[lower], mp.mpf(1))
        else:
            target = upper - VERTICES
            append(row, logical_to_column[lower], source_scale)
            append(row, logical_to_column[target], target_scale)
            append(row, VERTICES + logical_to_column[lower], source_strut)
            append(row, VERTICES + logical_to_column[target], target_strut)
    for edge in final_edges:
        global_index = index_data["edge_to_index"][edge]
        orbit, group = divmod(global_index, 24)
        row = 24 * orbit_position[orbit] + group
        append(row, logical_to_column[edge[0] - VERTICES], 1 / lam)
        append(row, logical_to_column[edge[1] - VERTICES], 1 / lam)
    return entries, logical_to_column


def project_entries(entries, sector):
    dimension = int(sector["dimension"])
    basis = sector["basis"]
    result = mp.matrix(65 * dimension, 10 * dimension)
    for row, column, value in entries:
        row_type, row_group = divmod(row, 24)
        column_type, column_group = divmod(column, 24)
        for left in range(dimension):
            factor = mp.conj(basis[row_group, left]) * value
            for right in range(dimension):
                result[row_type * dimension + left,
                       column_type * dimension + right] += (
                    factor * basis[column_group, right]
                )
    return result


def canonical_lift(block, weak_positions, old):
    dimension = 1
    internal = old["expanded_types"](30, 65, dimension)
    new = old["expanded_types"](65, 95, dimension)
    old_types = old["expanded_types"](0, 30, dimension)
    j_matrix = mp.matrix(65, 65)
    k_xx = old["mp_submatrix"](block, internal, internal)
    k_xn = old["mp_submatrix"](block, internal, new)
    k_ox = old["mp_submatrix"](block, old_types, internal)
    k_on = old["mp_submatrix"](block, old_types, new)
    for row in range(35):
        for column in range(35):
            j_matrix[row, column] = k_xx[row, column]
        for column in range(30):
            j_matrix[row, 35 + column] = k_xn[row, column]
    for row in range(30):
        for column in range(35):
            j_matrix[35 + row, column] = -k_ox[row, column]
        for column in range(30):
            j_matrix[35 + row, 35 + column] = -k_on[row, column]

    weak = list(weak_positions)
    strong = [index for index in range(65) if index not in set(weak)]
    a = old["mp_submatrix"](j_matrix, strong, strong)
    b = old["mp_submatrix"](j_matrix, strong, weak)
    a_ball = mp_matrix_to_acb(a)
    b_ball = mp_matrix_to_acb(b)
    response = -(a_ball.solve(b_ball))
    lift = acb_mat(65, 5)
    for ordered_row, original_row in enumerate(strong):
        for column in range(5):
            lift[original_row, column] = response[ordered_row, column]
    for column, original_row in enumerate(weak):
        lift[original_row, column] = 1
    return lift, a_ball.det()


def intersection_matrices(carrier_matrix, lift, scale_norm, strut_norm):
    carrier_ball = mp_matrix_to_acb(carrier_matrix)
    d_matrix = acb_mat(65, 10)
    k_matrix = acb_mat(65, 15)
    scale_ball = mp_to_acb(scale_norm)
    strut_ball = mp_to_acb(strut_norm)
    for row in range(65):
        for column in range(5):
            scale = carrier_ball[row, column] * scale_ball
            strut = carrier_ball[row, 5 + column] * strut_ball
            canonical = lift[row, column] * strut_ball
            d_matrix[row, column] = scale
            d_matrix[row, 5 + column] = strut - canonical
            k_matrix[row, column] = scale
            k_matrix[row, 5 + column] = strut
            k_matrix[row, 10 + column] = -canonical
    return d_matrix, k_matrix


def all_deleted_gram_records(matrix):
    records = []
    for removed in range(matrix.ncols()):
        columns = [index for index in range(matrix.ncols()) if index != removed]
        records.append({"removed": removed, **gram_determinant_record(
            acb_subcolumns(matrix, columns)
        )})
    return records


def synthetic_controls():
    planted = acb_mat(65, 10)
    for index in range(9):
        planted[index, index] = 1
        planted[index, 9] = 1
    identity = acb_mat(65, 10)
    for index in range(10):
        identity[index, index] = 1

    def one_line(matrix):
        full = gram_determinant_record(matrix)
        deleted = all_deleted_gram_records(matrix)
        return bool(full["contains_zero"] and all(
            not record["contains_zero"] for record in deleted
        ))

    return one_line(planted) and not one_line(identity)


print("=" * 78)
print("ADVERSARIAL DIRECT-MATRIX HOMOGENEOUS LINE REPLICATION")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
primary = json.loads(PRIMARY_RESOLUTION.read_text())
homogeneous = json.loads(HOMOGENEOUS_INPUT.read_text())
carrier = json.loads(CARRIER_INPUT.read_text())
intersection = json.loads(INTERSECTION_INPUT.read_text())
action_before = ACTION_INPUT.read_bytes()
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and primary["outcome"] == "HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE"
    and primary["passed"] == primary["tests"] == 10
    and homogeneous["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and homogeneous["passed"] == homogeneous["tests"] == 7
    and intersection["outcome"] == "FULL_SCALE_STRUT_CANONICAL_NUMERICALLY_OPEN"
)
check("all adversarial inputs retain frozen provenance", provenance_ok)
check("rank-nine and rank-ten synthetic controls separate", synthetic_controls())

# Re-execute the independent full-action builder while suppressing its own report.
captured = io.StringIO()
original_exit = sys.exit


def audited_exit(code=0):
    if code not in (None, 0):
        raise SystemExit(code)


try:
    sys.exit = audited_exit
    with contextlib.redirect_stdout(captured):
        old = runpy.run_path(str(ACTION_SOURCE))
finally:
    sys.exit = original_exit

action_rebuild_ok = bool(
    old["tests"] == old["passed"] == 14
    and old["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
    and ACTION_INPUT.read_bytes() == action_before
)
check("the frozen full-action builder reproduces without mutating its artifact", action_rebuild_ok)

constants = configure_action_precision(old)
tick = old["tick"]
records = {}
rank_ok = True
line_ok = True
control_ok = True
d_vectors = {}

for parity in ("even", "odd"):
    print(f"[{parity}] rebuilding the homogeneous full matrices", flush=True)
    model = old["models"][parity]
    state = tick["solutions"][parity]["state"]
    index_data = old["group_and_index_data"](model, state)
    geometry = old["prepare_geometry"](model, index_data)
    weak_positions = [
        position for position in range(35)
        if index_data["edge_kind"][24 * (30 + position)] == "pole"
    ]
    sectors, sector_control = old["high_precision_sector_bases"](index_data)
    homogeneous_sectors = [
        sector for sector in sectors
        if int(sector["dimension"]) == 1
        and abs(sector["constant_overlap"] - 1) < mp.mpf("1e-80")
    ]
    sector_selection_ok = len(homogeneous_sectors) == 1 and len(weak_positions) == 5
    control_ok &= sector_selection_ok
    if not sector_selection_ok:
        records[parity] = {"sector_selection_failed": True}
        continue
    sector = homogeneous_sectors[0]

    kind_values = {
        "old": constants["L0_square"],
        "internal": mp.exp(mp.mpf(state[0])) * constants["L0_square"] - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * mp.mpf(state[0])) * constants["L0_square"],
    }
    pattern_cache, branch_control = old["high_precision_pattern_cache"](
        geometry["patterns"], kind_values
    )
    kernels, kernel_control = old["assemble_full_representative_kernels"](
        index_data, geometry, pattern_cache
    )
    entries, logical_columns = direct_carrier_entries(
        parity, index_data, weak_positions, carrier, tick, constants
    )
    carrier_matrix = project_entries(entries, sector)
    block = old["project_full_kernel"](kernels["operational_primary"], sector)
    lift, strong_determinant = canonical_lift(block, weak_positions, old)
    scales = intersection["parities"][parity]["sectors"][6]["scaling"]
    scale_norm = mp.mpf(scales["scale"])
    strut_norm = mp.mpf(scales["strut"])
    d_ball, k_ball = intersection_matrices(
        carrier_matrix, lift, scale_norm, strut_norm
    )
    d_matrix = acb_matrix_midpoint(d_ball)
    k_matrix = acb_matrix_midpoint(k_ball)

    d_deleted = all_deleted_gram_records(d_ball)
    k_deleted = all_deleted_gram_records(k_ball)
    parity_rank_ok = bool(
        not strong_determinant.contains(0)
        and all(not item["contains_zero"] for item in d_deleted)
        and all(not item["contains_zero"] for item in k_deleted)
    )
    rank_ok &= parity_rank_ok

    d_vector = deterministic_normal_null(d_matrix)
    k_vector = deterministic_normal_null(k_matrix)
    d_vectors[parity] = d_vector
    d_residual = normalized_residual(d_matrix, d_vector)
    k_residual = normalized_residual(k_matrix, k_vector)
    scale_values = [d_vector[index] for index in range(5)]
    strut_values = [d_vector[5 + index] for index in range(5)]
    scale_mean = mp.fsum(scale_values) / 5
    strut_mean = mp.fsum(strut_values) / 5
    scale_relative_spread = max(
        abs(value - scale_mean) for value in scale_values
    ) / abs(scale_mean)
    strut_relative_spread = max(
        abs(value - strut_mean) for value in strut_values
    ) / abs(strut_mean)
    joined = mp.matrix(list(d_vector[:5]) + list(d_vector[5:]) + list(d_vector[5:]))
    joined_distance = projector_distance(k_vector, joined)
    direct_ratio = scale_norm * scale_mean / (strut_norm * strut_mean)
    primary_ratio = mp.mpf(primary["bridges"][parity]["analytic_sigma_over_c"])
    ratio_error = abs(direct_ratio - primary_ratio)

    lam = mp.mpf(primary["bridges"][parity]["lambda"])
    missing_lambda_ratio = primary_ratio / lam
    wrong_sign_ratio = -primary_ratio

    def physical_vector(ratio):
        result = mp.matrix(10, 1)
        for index in range(5):
            result[index] = ratio / scale_norm
            result[5 + index] = 1 / strut_norm
        return normalize_phase(result)

    missing_lambda_residual = normalized_residual(
        d_matrix, physical_vector(missing_lambda_ratio)
    )
    wrong_sign_residual = normalized_residual(
        d_matrix, physical_vector(wrong_sign_ratio)
    )
    parity_line_ok = bool(
        d_residual < mp.mpf("1e-30")
        and k_residual < mp.mpf("1e-30")
        and scale_relative_spread < mp.mpf("1e-30")
        and strut_relative_spread < mp.mpf("1e-30")
        and joined_distance < mp.mpf("1e-30")
        and ratio_error < mp.mpf("1e-30")
    )
    parity_control_ok = bool(
        len(entries) == 4440
        and set(logical_columns) == set(range(VERTICES))
        and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and branch_control["entry_pass"]
        and missing_lambda_residual > mp.mpf("1e-10")
        and wrong_sign_residual > mp.mpf("1e-3")
    )
    line_ok &= parity_line_ok
    control_ok &= parity_control_ok
    records[parity] = {
        "D_deleted_gram_determinants": d_deleted,
        "K_deleted_gram_determinants": k_deleted,
        "strong_solve_determinant": str(strong_determinant),
        "D_normalized_residual": mp_text(d_residual),
        "K_normalized_residual": mp_text(k_residual),
        "D_scale_relative_spread": mp_text(scale_relative_spread),
        "D_strut_relative_spread": mp_text(strut_relative_spread),
        "K_joined_D_projector_distance": mp_text(joined_distance),
        "direct_sigma_over_c": mp_text(direct_ratio),
        "primary_sigma_over_c": mp_text(primary_ratio),
        "ratio_absolute_error": mp_text(ratio_error),
        "missing_lambda_normalized_residual": mp_text(missing_lambda_residual),
        "wrong_sign_normalized_residual": mp_text(wrong_sign_residual),
        "entry_count": len(entries),
        "maximum_kernel_imaginary": mp_text(kernel_control["maximum_imaginary"]),
        "rank_ok": parity_rank_ok,
        "line_ok": parity_line_ok,
        "controls_ok": parity_control_ok,
    }

parity_distance = (
    projector_distance(d_vectors["even"], d_vectors["odd"])
    if set(d_vectors) == {"even", "odd"} else mp.inf
)
line_ok &= parity_distance < mp.mpf("1e-60")

check("both direct P160 matrix reconstructions pass hostile controls", control_ok)
check("all 50 frozen-choice deleted Gram determinants exclude zero", rank_ok)
check("normal-equation D/K lines agree with each other and the primary line", line_ok,
      f"parity projector distance={mp_text(parity_distance, 12)}")

if not provenance_ok or not action_rebuild_ok or not control_ok:
    outcome = "HOMOGENEOUS_ADVERSARIAL_CONTROL_FAILED"
elif not rank_ok:
    outcome = "HOMOGENEOUS_ADVERSARIAL_RANK_DISAGREEMENT"
elif not line_ok:
    outcome = "HOMOGENEOUS_ADVERSARIAL_LINE_DISAGREEMENT"
else:
    outcome = "HOMOGENEOUS_WEAK_POLE_LINE_REPLICATED"

allowed = {
    "HOMOGENEOUS_ADVERSARIAL_CONTROL_FAILED",
    "HOMOGENEOUS_ADVERSARIAL_RANK_DISAGREEMENT",
    "HOMOGENEOUS_ADVERSARIAL_LINE_DISAGREEMENT",
    "HOMOGENEOUS_WEAK_POLE_LINE_REPLICATED",
}
check("the preregistered adversarial hierarchy assigns one verdict", outcome in allowed, outcome)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "primary_result_commit": PRIMARY_RESULT_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "method": {
        "precision": "P160/Arb P140",
        "candidate_extraction": "fixed-last-component normal equations",
        "rank_lower_bound": "all single-column-deleted Gram determinants",
        "symbolic_primary_generator_used_to_construct_candidate": False,
        "stored_P100_SVD_candidate_used": False,
    },
    "parities": records,
    "even_odd_D_projector_distance": mp_text(parity_distance),
    "classification": {
        "homogeneous_weak_pole_line": (
            "ONE DIMENSION; ADVERSARIALLY REPLICATED"
            if outcome == "HOMOGENEOUS_WEAK_POLE_LINE_REPLICATED" else "OPEN"
        ),
        "omitted_pole_equation": "NOT EVALUATED",
        "gauge_or_physical": "OPEN",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"TOTAL: {passed}/{tests} tests PASSED")
print(f"Artifact: {OUTPUT.name}")
if passed != tests:
    raise SystemExit(1)


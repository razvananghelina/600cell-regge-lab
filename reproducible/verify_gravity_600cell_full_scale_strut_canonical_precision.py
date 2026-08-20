#!/usr/bin/env python3
"""Multiprecision resolver for the complete carrier/action intersection."""

from collections import Counter
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
PRIOR = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_canonical_precision_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_canonical_precision_protocol.md"
FIRST_RESULT = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_canonical_intersection_first_result.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_canonical_intersection.py"
PRIMARY_INPUT = HERE / "gravity_600cell_full_scale_strut_canonical_intersection.json"
FIRST_FAILURE = HERE / "gravity_600cell_full_scale_strut_canonical_intersection_first_failure.json"
CARRIER_INPUT = HERE / "gravity_600cell_full_scale_strut_carrier.json"
SYMBOLIC_INPUT = HERE / "gravity_600cell_full_scale_strut_symbolic_gap_resolution.json"
ACTION_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
ACTION_INPUT = HERE / "gravity_600cell_dust_hyperbolic_lapse_alignment.json"
HOMOGENEOUS_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
OUTPUT = HERE / "gravity_600cell_full_scale_strut_canonical_precision.json"

PRIOR_COMMIT = "03939f8"
PROTOCOL_COMMIT = "f011db5"
EXPECTED_HASHES = {
    "prior": "fd2e230fdc0c0f7aaa771a4781973d0b476b5758d8d82f5b84ba471b391a722c",
    "protocol": "f044c0738fc7f507b89b1bc3658836ba5fa7a1d34f00f533bf821146663686b0",
    "first_result": "971b1eddcb09e4a35a72ba5d1c359b710ba74673815824c1029ddb68e897bfc5",
    "primary_source": "a2d5390d39c725a5fb586fefce9da34cede3a1fb84bbe36791f8b0599b3eae42",
    "primary": "b29cc33a9effeb2087fb6133359ee747d100d203778586372a7ceeebc2e4f070",
    "first_failure": "6423c3efc03ba6107a82c1b0d813e0226ccf757d242cc3ecc0522003095e97d5",
    "symbolic": "ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179",
    "action_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "action": "a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff",
    "homogeneous": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
}
INPUTS = {
    "prior": PRIOR,
    "protocol": PROTOCOL,
    "first_result": FIRST_RESULT,
    "primary_source": PRIMARY_SOURCE,
    "primary": PRIMARY_INPUT,
    "first_failure": FIRST_FAILURE,
    "symbolic": SYMBOLIC_INPUT,
    "action_source": ACTION_SOURCE,
    "action": ACTION_INPUT,
    "homogeneous": HOMOGENEOUS_INPUT,
}

LEVELS = {
    "P100": {
        "dps": 100,
        "ball_dps": 80,
        "digits": 85,
        "steps": {
            "operational_primary": "1e-20",
            "operational_shadow": "1e-15",
            "validation_primary": "3e-20",
            "validation_shadow": "3e-15",
        },
        "floor": "1e-70",
    },
    "P160": {
        "dps": 160,
        "ball_dps": 140,
        "digits": 145,
        "steps": {
            "operational_primary": "1e-40",
            "operational_shadow": "1e-30",
            "validation_primary": "3e-40",
            "validation_shadow": "3e-30",
        },
        "floor": "1e-130",
    },
}
VERTICES = 120
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


def mp_string(value, digits=50):
    if isinstance(value, mp.mpc):
        return {
            "real": mp.nstr(mp.re(value), digits),
            "imaginary": mp.nstr(mp.im(value), digits),
        }
    return mp.nstr(value, digits)


def mp_norm(vector):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in vector))


def arb_mid_to_mp(value):
    mantissa, exponent = value.mid().man_exp()
    return mp.mpf(int(mantissa)) * mp.power(2, int(exponent))


def arb_upper_to_mp(value):
    mantissa, exponent = value.upper().man_exp()
    return mp.mpf(int(mantissa)) * mp.power(2, int(exponent))


def acb_mid_to_mp(value):
    return mp.mpc(
        arb_mid_to_mp(value.real),
        arb_mid_to_mp(value.imag),
    )


def acb_radius_to_mp(value):
    return mp.sqrt(
        arb_upper_to_mp(value.real.rad()) ** 2
        + arb_upper_to_mp(value.imag.rad()) ** 2
    )


def mp_to_acb(value, digits):
    return acb(
        mp.nstr(mp.re(value), digits),
        mp.nstr(mp.im(value), digits),
    )


def mp_matrix_to_acb(matrix, digits):
    return acb_mat(
        matrix.rows,
        matrix.cols,
        [mp_to_acb(matrix[row, column], digits)
         for row in range(matrix.rows) for column in range(matrix.cols)],
    )


def acb_matrix_midpoint(matrix):
    return mp.matrix([
        [acb_mid_to_mp(matrix[row, column]) for column in range(matrix.ncols())]
        for row in range(matrix.nrows())
    ])


def acb_matrix_radius(matrix):
    return mp.sqrt(mp.fsum(
        acb_radius_to_mp(matrix[row, column]) ** 2
        for row in range(matrix.nrows()) for column in range(matrix.ncols())
    ))


def acb_subcolumns(matrix, columns):
    return acb_mat(
        matrix.nrows(), len(columns),
        [matrix[row, column] for row in range(matrix.nrows()) for column in columns],
    )


def acb_det_record(matrix):
    gram = matrix.conjugate().transpose() * matrix
    determinant = gram.det()
    return {
        "ball": str(determinant),
        "contains_zero": bool(determinant.contains(0)),
        "abs_lower": str(determinant.abs_lower()),
        "abs_upper": str(determinant.abs_upper()),
        "radius": str(determinant.rad()),
    }


def gram_spectrum(matrix, vectors=False):
    gram = matrix.H * matrix
    values, eigenvectors = mp.eighe(gram)
    tolerance = mp.power(10, -(mp.mp.dps - 20)) * max(
        mp.mpf(1), max(abs(values[index]) for index in range(len(values)))
    )
    negative_bad = any(values[index] < -tolerance for index in range(len(values)))
    singulars = [
        mp.sqrt(max(mp.mpf(0), values[index])) for index in range(len(values))
    ]
    return singulars, eigenvectors if vectors else None, negative_bad


def configure_precision(old, level):
    specification = LEVELS[level]
    mp.mp.dps = specification["dps"]
    ctx.dps = specification["ball_dps"]
    steps = {name: mp.mpf(value) for name, value in specification["steps"].items()}

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
        namespace["DERIVATIVE_STEPS"] = steps
        namespace["VARIANTS"] = tuple(steps)
        namespace["ARITHMETIC_FLOOR"] = mp.mpf(specification["floor"])
    old["L0_SQUARE"] = l0_square
    old["MASS"] = mass
    old["RHO0"] = rho0
    old["DERIVATIVE_STEPS"] = steps
    old["VARIANTS"] = tuple(steps)
    return {
        "L0_square": l0_square,
        "mass": mass,
        "rho0": rho0,
        "steps": steps,
        "digits": specification["digits"],
    }


def carrier_sparse(parity, index_data, weak_positions, committed, constants, corrupt=False):
    l0_square = constants["L0_square"]
    state = primary["parities"][parity]
    del state  # provenance only; geometry state comes from the accepted action input
    lam = mp.exp(mp.mpf(action_tick["solutions"][parity]["state"][0]))
    rho = index_data["rho"]
    q_diag = lam * l0_square - rho
    a_raw = -16 * rho / (l0_square * (lam - 1) ** 2)
    b_raw = 8 + 16 * rho / (l0_square * (lam - 1) ** 2)
    scale_factor = l0_square / (8 * q_diag)
    kappa = rho / ((lam - 1) * q_diag)

    orbit_order = tuple(range(30, 65)) + tuple(range(65, 95))
    orbit_position = {orbit: position for position, orbit in enumerate(orbit_order)}
    weak_types = [30 + position for position in weak_positions]
    pole_edges = [
        tuple(map(int, edge))
        for orbit in weak_types for edge in index_data["orbit_edges"][orbit]
    ]
    logical_to_column = {edge[0]: column for column, edge in enumerate(pole_edges)}
    record = committed["parities"][parity]
    internal_edges = [tuple(map(int, edge)) for edge in record["internal_edge_order"]]
    final_edges = [tuple(map(int, edge)) for edge in record["final_edge_order"]]
    first_diagonal = min(
        edge for edge in internal_edges
        if edge[0] < VERTICES <= edge[1] and edge[1] != edge[0] + VERTICES
    )
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
            source_scale = scale_factor * a_raw
            target_scale = scale_factor * b_raw
            source_strut = kappa
            target_strut = -lam * kappa
            if corrupt and edge == first_diagonal:
                source_scale, target_scale = target_scale, source_scale
                source_strut, target_strut = target_strut, source_strut
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
    return entries, logical_to_column, first_diagonal


def project_sparse(entries, sector):
    dimension = sector["dimension"]
    basis = sector["basis"]
    result = mp.matrix(65 * dimension, 10 * dimension)
    for row, column, value in entries:
        row_type, row_group = divmod(row, 24)
        column_type, column_group = divmod(column, 24)
        for left in range(dimension):
            factor = mp.conj(basis[row_group, left]) * value
            for right in range(dimension):
                result[
                    row_type * dimension + left,
                    column_type * dimension + right,
                ] += factor * basis[column_group, right]
    return result


def canonical_lift_ball(block, dimension, weak_positions, old, digits):
    internal = old["expanded_types"](30, 65, dimension)
    new = old["expanded_types"](65, 95, dimension)
    old_types = old["expanded_types"](0, 30, dimension)
    xd = 35 * dimension
    nd = 30 * dimension
    j_matrix = mp.matrix(xd + nd, xd + nd)
    k_xx = old["mp_submatrix"](block, internal, internal)
    k_xn = old["mp_submatrix"](block, internal, new)
    k_ox = old["mp_submatrix"](block, old_types, internal)
    k_on = old["mp_submatrix"](block, old_types, new)
    for row in range(xd):
        for column in range(xd):
            j_matrix[row, column] = k_xx[row, column]
        for column in range(nd):
            j_matrix[row, xd + column] = k_xn[row, column]
    for row in range(nd):
        for column in range(xd):
            j_matrix[xd + row, column] = -k_ox[row, column]
        for column in range(nd):
            j_matrix[xd + row, xd + column] = -k_on[row, column]

    weak = [
        position * dimension + component
        for position in weak_positions for component in range(dimension)
    ]
    weak_set = set(weak)
    strong = [index for index in range(xd + nd) if index not in weak_set]
    a_mp = old["mp_submatrix"](j_matrix, strong, strong)
    b_mp = old["mp_submatrix"](j_matrix, strong, weak)
    a_ball = mp_matrix_to_acb(a_mp, digits)
    b_ball = mp_matrix_to_acb(b_mp, digits)
    solve_ball = a_ball.solve(b_ball)
    ordered = acb_mat(xd + nd, len(weak))
    for row in range(len(strong)):
        for column in range(len(weak)):
            ordered[row, column] = -solve_ball[row, column]
    for index in range(len(weak)):
        ordered[len(strong) + index, index] = 1
    order = strong + weak
    lift = acb_mat(xd + nd, len(weak))
    for ordered_row, original_row in enumerate(order):
        for column in range(len(weak)):
            lift[original_row, column] = ordered[ordered_row, column]
    return lift, a_ball.det(), weak


def scaled_intersection_matrices(g_mp, c_ball, scale_value, strut_value, digits):
    rows = g_mp.rows
    half = g_mp.cols // 2
    g_ball = mp_matrix_to_acb(g_mp, digits)
    d_ball = acb_mat(rows, 2 * half)
    k_ball = acb_mat(rows, 3 * half)
    scale_ball = mp_to_acb(scale_value, digits)
    strut_ball = mp_to_acb(strut_value, digits)
    for row in range(rows):
        for column in range(half):
            gs = g_ball[row, column] * scale_ball
            gt = g_ball[row, half + column] * strut_ball
            cc = c_ball[row, column] * strut_ball
            d_ball[row, column] = gs
            d_ball[row, half + column] = gt - cc
            k_ball[row, column] = gs
            k_ball[row, half + column] = gt
            k_ball[row, 2 * half + column] = -cc
    return d_ball, k_ball


def freeze_candidate(vector):
    index = max(range(len(vector)), key=lambda item: abs(vector[item]))
    if abs(vector[index]):
        phase = mp.conj(vector[index]) / abs(vector[index])
        vector = vector * phase
    norm = mp_norm(vector)
    vector = vector / norm
    return mp.matrix([
        mp.mpc(
            mp.mpf(mp.nstr(mp.re(value), 70)),
            mp.mpf(mp.nstr(mp.im(value), 70)),
        )
        for value in vector
    ])


def synthetic_controls(digits):
    full = acb_mat(6, 4)
    for index in range(4):
        full[index, index] = 1
    planted = acb_mat(6, 4)
    for index in range(3):
        planted[index, index] = 1
    full_record = acb_det_record(full)
    planted_record = acb_det_record(planted)
    minor = acb_subcolumns(planted, (0, 1, 2))
    minor_record = acb_det_record(minor)
    return bool(
        not full_record["contains_zero"]
        and planted_record["contains_zero"]
        and not minor_record["contains_zero"]
    ), {
        "full": full_record,
        "planted": planted_record,
        "planted_rank_minor": minor_record,
        "digits": digits,
    }


print("=" * 78)
print("MULTIPRECISION COMPLETE-CARRIER INTERSECTION RESOLVER")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
primary = json.loads(PRIMARY_INPUT.read_text())
first_failure = json.loads(FIRST_FAILURE.read_text())
symbolic = json.loads(SYMBOLIC_INPUT.read_text())
action_before = json.loads(ACTION_INPUT.read_text())
homogeneous_input = json.loads(HOMOGENEOUS_INPUT.read_text())
carrier_input = json.loads(CARRIER_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and primary["outcome"] == "FULL_SCALE_STRUT_CANONICAL_NUMERICALLY_OPEN"
    and primary["passed"] == primary["tests"] == 13
    and not primary["all_resolved"]
    and all(
        not sector["resolved"]
        for parity in primary["parities"].values()
        for sector in parity["sectors"]
    )
    and first_failure["outcome"] == "FULL_SCALE_STRUT_CANONICAL_CONTROL_FAILED"
    and symbolic["outcome"] == "FULL_SCALE_STRUT_GAP_REAL_RESOLVED"
    and homogeneous_input["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and digest(CARRIER_INPUT) == primary["input_sha256"]["carrier"]
)
check("all precision-resolver inputs retain frozen provenance", provenance_ok)

print("loading the frozen action construction", flush=True)
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
action_ok = bool(
    old["tests"] == old["passed"] == 14
    and old["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
    and digest(ACTION_INPUT) == EXPECTED_HASHES["action"]
    and json.loads(ACTION_INPUT.read_text()) == action_before
)
check("the frozen action construction reproduces before precision changes", action_ok)

models = old["models"]
action_tick = old["tick"]
level_records = {}
level_live = {}
frozen_candidates = {}
all_level_controls = True

for level in ("P100", "P160"):
    specification = LEVELS[level]
    constants = configure_precision(old, level)
    synthetic_ok, synthetic_record = synthetic_controls(specification["digits"])
    check(f"{level}: interval determinant synthetic controls pass", synthetic_ok)
    all_level_controls &= synthetic_ok
    level_records[level] = {"synthetic": synthetic_record, "parities": {}}
    level_live[level] = {}

    for parity in ("even", "odd"):
        print(f"[{level}/{parity}] rebuilding geometry and action kernels", flush=True)
        model = models[parity]
        state = action_tick["solutions"][parity]["state"]
        index_data = old["group_and_index_data"](model, state)
        geometry = old["prepare_geometry"](model, index_data)
        weak_positions = [
            position for position in range(35)
            if index_data["edge_kind"][24 * (30 + position)] == "pole"
        ]
        sectors, sector_control = old["high_precision_sector_bases"](index_data)
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
        entries, logical_to_column, first_diagonal = carrier_sparse(
            parity, index_data, weak_positions, carrier_input, constants
        )
        corrupt_entries, _, _ = carrier_sparse(
            parity, index_data, weak_positions, carrier_input, constants, corrupt=True
        )
        geometry_ok = bool(
            set(logical_to_column) == set(range(VERTICES))
            and len(entries) == 4440
            and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
            and branch_control["entry_pass"]
            and branch_control["base_negative_counts"] == Counter({1: 2400})
            and branch_control["displaced_negative_counts"] == Counter({1: 1600})
            and kernel_control["maximum_imaginary"] < mp.mpf(specification["floor"])
        )
        check(f"{level}/{parity}: multiprecision geometry controls pass", geometry_ok)
        all_level_controls &= geometry_ok

        parity_records = []
        level_live[level][parity] = {}
        for sector_index, sector in enumerate(sectors):
            dimension = int(sector["dimension"])
            print(
                f"[{level}/{parity}] sector {sector_index + 1}/7 d={dimension}",
                flush=True,
            )
            g_mp = project_sparse(entries, sector)
            g_corrupt_mp = project_sparse(corrupt_entries, sector)
            block = old["project_full_kernel"](
                kernels["operational_primary"], sector
            )
            c_ball, strong_det, weak_rows = canonical_lift_ball(
                block, dimension, weak_positions, old, specification["digits"]
            )
            scaling = primary["parities"][parity]["sectors"][sector_index]["scaling"]
            scale_value = mp.mpf(scaling["scale"])
            strut_value = mp.mpf(scaling["strut"])
            d_ball, k_ball = scaled_intersection_matrices(
                g_mp, c_ball, scale_value, strut_value, specification["digits"]
            )
            d_mp = acb_matrix_midpoint(d_ball)
            k_mp = acb_matrix_midpoint(k_ball)
            d_singulars, d_vectors, d_negative = gram_spectrum(d_mp, vectors=True)
            k_singulars, _, k_negative = gram_spectrum(k_mp)
            d_det = acb_det_record(d_ball)
            k_det = acb_det_record(k_ball)
            d_radius = acb_matrix_radius(d_ball)
            k_radius = acb_matrix_radius(k_ball)
            homogeneous = abs(sector["constant_overlap"] - 1) < mp.mpf("1e-50")

            d_minor_records = []
            k_minor_records = []
            if homogeneous:
                for removed in range(d_ball.ncols()):
                    columns = [index for index in range(d_ball.ncols()) if index != removed]
                    d_minor_records.append({
                        "removed": removed,
                        **acb_det_record(acb_subcolumns(d_ball, columns)),
                    })
                for removed in range(k_ball.ncols()):
                    columns = [index for index in range(k_ball.ncols()) if index != removed]
                    k_minor_records.append({
                        "removed": removed,
                        **acb_det_record(acb_subcolumns(k_ball, columns)),
                    })
                if level == "P100":
                    frozen_candidates[parity] = freeze_candidate(d_vectors[:, 0])

            c_mid = acb_matrix_midpoint(c_ball)
            half = 5 * dimension
            d_corrupt = mp.matrix(d_mp.rows, d_mp.cols)
            for row in range(d_mp.rows):
                for column in range(half):
                    d_corrupt[row, column] = g_corrupt_mp[row, column] * scale_value
                    d_corrupt[row, half + column] = (
                        g_corrupt_mp[row, half + column] - c_mid[row, column]
                    ) * strut_value
            corrupt_singulars, _, corrupt_negative = gram_spectrum(d_corrupt)
            corruption_change = max(
                abs(corrupt_singulars[index] - d_singulars[index])
                for index in range(len(d_singulars))
            )

            conjugate_singulars, _, conjugate_negative = gram_spectrum(
                d_mp.apply(mp.conj)
            )
            conjugation_error = max(
                abs(conjugate_singulars[index] - d_singulars[index])
                for index in range(len(d_singulars))
            )
            sector_control_ok = bool(
                not d_negative and not k_negative and not corrupt_negative
                and not conjugate_negative
                and not strong_det.contains(0)
                and conjugation_error < mp.power(10, -(specification["dps"] - 30))
            )
            all_level_controls &= sector_control_ok
            parity_records.append({
                "sector_index": sector_index,
                "dimension": dimension,
                "constant_overlap": mp_string(sector["constant_overlap"]),
                "homogeneous": homogeneous,
                "strong_solve_determinant": str(strong_det),
                "D_determinant": d_det,
                "K_determinant": k_det,
                "D_radius": mp_string(d_radius),
                "K_radius": mp_string(k_radius),
                "D_singular_values_ascending": [mp_string(value) for value in d_singulars],
                "K_singular_values_ascending": [mp_string(value) for value in k_singulars],
                "D_single_column_deleted_minors": d_minor_records,
                "K_single_column_deleted_minors": k_minor_records,
                "source_target_corruption_max_singular_change": mp_string(corruption_change),
                "conjugation_singular_error": mp_string(conjugation_error),
                "controls_pass": sector_control_ok,
            })
            level_live[level][parity][sector_index] = {
                "dimension": dimension,
                "homogeneous": homogeneous,
                "D": d_mp,
                "K": k_mp,
                "D_singulars": d_singulars,
                "K_singulars": k_singulars,
                "D_radius": d_radius,
                "K_radius": k_radius,
                "D_det": d_det,
                "K_det": k_det,
                "D_minors": d_minor_records,
                "K_minors": k_minor_records,
                "corruption_change": corruption_change,
                "weak_rows": weak_rows,
            }
        level_records[level]["parities"][parity] = parity_records

check("both precision levels retain every geometry and arithmetic control", all_level_controls)

# Cross-level classification.
nonhomogeneous_full_rank = True
nonhomogeneous_disagreement = False
homogeneous_status = {}
no_refit_records = {}
minimum_nonhomogeneous_margin = mp.inf
maximum_nonhomogeneous_relative_change = mp.mpf(0)

for parity in ("even", "odd"):
    homogeneous_indices = [
        index for index, record in level_live["P160"][parity].items()
        if record["homogeneous"]
    ]
    homogeneous_index_ok = len(homogeneous_indices) == 1
    check(f"{parity}: exactly one constant-overlap sector is identified", homogeneous_index_ok)
    all_level_controls &= homogeneous_index_ok
    homogeneous_index = homogeneous_indices[0]

    for sector_index, high in level_live["P160"][parity].items():
        low = level_live["P100"][parity][sector_index]
        if high["homogeneous"]:
            continue
        d_stable = abs(high["D_singulars"][0] - low["D_singulars"][0]) / high["D_singulars"][0]
        k_stable = abs(high["K_singulars"][0] - low["K_singulars"][0]) / high["K_singulars"][0]
        stability = max(d_stable, k_stable)
        maximum_nonhomogeneous_relative_change = max(
            maximum_nonhomogeneous_relative_change, stability
        )
        margin = min(
            high["D_singulars"][0] / max(high["D_radius"], mp.mpf("1e-150")),
            high["K_singulars"][0] / max(high["K_radius"], mp.mpf("1e-150")),
        )
        minimum_nonhomogeneous_margin = min(minimum_nonhomogeneous_margin, margin)
        full_rank = bool(
            not low["D_det"]["contains_zero"]
            and not low["K_det"]["contains_zero"]
            and not high["D_det"]["contains_zero"]
            and not high["K_det"]["contains_zero"]
            and stability < mp.mpf("1e-20")
            and margin > mp.mpf("1e20")
        )
        nonhomogeneous_full_rank &= full_rank
        binary_minimum = mp.mpf(
            primary["parities"][parity]["sectors"][sector_index]
            ["D_singular_values"]["operational_primary"]["gesvd"][-1]
        )
        nonhomogeneous_disagreement |= bool(
            abs(high["D_singulars"][0] - binary_minimum)
            > mp.mpf("1e-5") * max(high["D_singulars"][0], binary_minimum)
        )

    high = level_live["P160"][parity][homogeneous_index]
    low = level_live["P100"][parity][homogeneous_index]
    candidate = frozen_candidates[parity]
    half = 5
    joined_candidate = mp.matrix(15, 1)
    for index in range(half):
        joined_candidate[index] = candidate[index]
        joined_candidate[half + index] = candidate[half + index]
        joined_candidate[2 * half + index] = candidate[half + index]
    d_residual = mp_norm(high["D"] * candidate)
    k_residual = mp_norm(high["K"] * joined_candidate)
    next_d = high["D_singulars"][1]
    next_k = high["K_singulars"][1]
    d_minor_rank = any(not record["contains_zero"] for record in high["D_minors"])
    k_minor_rank = any(not record["contains_zero"] for record in high["K_minors"])
    d_decay = low["D_singulars"][0] / max(high["D_singulars"][0], mp.mpf("1e-160"))
    k_decay = low["K_singulars"][0] / max(high["K_singulars"][0], mp.mpf("1e-160"))
    d_next_stability = abs(next_d - low["D_singulars"][1]) / next_d
    k_next_stability = abs(next_k - low["K_singulars"][1]) / next_k
    midpoint_small_counts = {
        "D": sum(value < mp.mpf("1e-50") for value in high["D_singulars"]),
        "K": sum(value < mp.mpf("1e-50") for value in high["K_singulars"]),
    }
    one_kernel = bool(
        low["D_det"]["contains_zero"] and low["K_det"]["contains_zero"]
        and high["D_det"]["contains_zero"] and high["K_det"]["contains_zero"]
        and d_minor_rank and k_minor_rank
        and midpoint_small_counts == {"D": 1, "K": 1}
        and next_d > mp.mpf("1e-8") and next_k > mp.mpf("1e-8")
        and d_decay > mp.mpf("1e20") and k_decay > mp.mpf("1e20")
        and d_next_stability < mp.mpf("1e-20")
        and k_next_stability < mp.mpf("1e-20")
        and d_residual < mp.mpf("1e-50")
        and k_residual < mp.mpf("1e-50")
        and d_residual < mp.mpf("1e-40") * next_d
        and k_residual < mp.mpf("1e-40") * next_k
    )
    zero_kernel = bool(
        not low["D_det"]["contains_zero"]
        and not low["K_det"]["contains_zero"]
        and not high["D_det"]["contains_zero"]
        and not high["K_det"]["contains_zero"]
        and high["D_singulars"][0] / max(high["D_radius"], mp.mpf("1e-150")) > mp.mpf("1e20")
        and high["K_singulars"][0] / max(high["K_radius"], mp.mpf("1e-150")) > mp.mpf("1e20")
    )
    homogeneous_status[parity] = (
        "ONE" if one_kernel else "ZERO" if zero_kernel else "OPEN"
    )
    no_refit_records[parity] = {
        "sector_index": homogeneous_index,
        "candidate_P100_rounded_70_digits": [mp_string(value, 75) for value in candidate],
        "D_P160_residual": mp_string(d_residual),
        "K_P160_residual": mp_string(k_residual),
        "D_next_singular": mp_string(next_d),
        "K_next_singular": mp_string(next_k),
        "D_decay_P100_over_P160": mp_string(d_decay),
        "K_decay_P100_over_P160": mp_string(k_decay),
        "D_next_relative_change": mp_string(d_next_stability),
        "K_next_relative_change": mp_string(k_next_stability),
        "P160_small_counts_below_1e-50": midpoint_small_counts,
        "D_rank_at_least_9_minor": d_minor_rank,
        "K_rank_at_least_14_minor": k_minor_rank,
        "classification": homogeneous_status[parity],
    }

candidate_even = frozen_candidates["even"]
candidate_odd = frozen_candidates["odd"]
projector_difference = mp.sqrt(mp.fsum(
    abs(
        candidate_even[row] * mp.conj(candidate_even[column])
        - candidate_odd[row] * mp.conj(candidate_odd[column])
    ) ** 2
    for row in range(10) for column in range(10)
))
parity_candidate_ok = projector_difference < mp.mpf("1e-30")
check(
    "the two frozen homogeneous candidate projectors agree",
    parity_candidate_ok,
    f"difference={mp.nstr(projector_difference, 8)}",
)
all_level_controls &= parity_candidate_ok

source_target_margins = [
    record["corruption_change"]
    / max(record["D_radius"], mp.mpf("1e-150"))
    for parity in ("even", "odd")
    for record in level_live["P160"][parity].values()
    if not record["homogeneous"]
]
source_target_ok = bool(source_target_margins and max(source_target_margins) > mp.mpf("1e20"))
check(
    "the source/target reversal is resolved above P160 uncertainty",
    source_target_ok,
    f"maximum margin={mp.nstr(max(source_target_margins), 8)}",
)
all_level_controls &= source_target_ok

# Hostile pole deletion on the no-refit candidate.
pole_deletion_records = {}
pole_deletion_ok = True
for parity in ("even", "odd"):
    index = no_refit_records[parity]["sector_index"]
    high = level_live["P160"][parity][index]
    candidate = frozen_candidates[parity]
    corrupted = high["D"].copy()
    # The first weak row and first strut column contain the literal identity.
    pole_row = high["weak_rows"][0]
    corrupted[pole_row, 5] -= mp.mpf(
        primary["parities"][parity]["sectors"][index]["scaling"]["strut"]
    )
    residual = mp_norm(corrupted * candidate)
    baseline = mp.mpf(no_refit_records[parity]["D_P160_residual"])
    ratio = residual / max(baseline, mp.mpf("1e-150"))
    corrupted_det = acb_det_record(mp_matrix_to_acb(corrupted, LEVELS["P160"]["digits"]))
    determinant_changed = corrupted_det != high["D_det"]
    passed_control = bool(ratio > mp.mpf("1e20") and determinant_changed)
    pole_deletion_ok &= passed_control
    pole_deletion_records[parity] = {
        "residual": mp_string(residual),
        "baseline_ratio": mp_string(ratio),
        "determinant": corrupted_det,
        "determinant_changed": determinant_changed,
        "passed": passed_control,
    }
check("pole-identity deletion destroys the homogeneous candidate", pole_deletion_ok)
all_level_controls &= pole_deletion_ok

if not all_level_controls:
    outcome = "FULL_SCALE_STRUT_CANONICAL_PRECISION_CONTROL_FAILED"
elif nonhomogeneous_disagreement:
    outcome = "FULL_SCALE_STRUT_CANONICAL_PRECISION_DISAGREEMENT"
elif not nonhomogeneous_full_rank:
    outcome = "FULL_SCALE_STRUT_CANONICAL_NONHOMOGENEOUS_OPEN"
elif any(value == "OPEN" for value in homogeneous_status.values()):
    outcome = "FULL_SCALE_STRUT_CANONICAL_HOMOGENEOUS_OPEN"
elif set(homogeneous_status.values()) == {"ZERO"}:
    outcome = "FULL_SCALE_STRUT_CANONICAL_ZERO_INTERSECTION_RESOLVED"
elif set(homogeneous_status.values()) == {"ONE"}:
    outcome = "FULL_SCALE_STRUT_CANONICAL_ONE_HOMOGENEOUS_RESOLVED"
else:
    outcome = "FULL_SCALE_STRUT_CANONICAL_PRECISION_DISAGREEMENT"

allowed = {
    "FULL_SCALE_STRUT_CANONICAL_PRECISION_CONTROL_FAILED",
    "FULL_SCALE_STRUT_CANONICAL_PRECISION_DISAGREEMENT",
    "FULL_SCALE_STRUT_CANONICAL_NONHOMOGENEOUS_OPEN",
    "FULL_SCALE_STRUT_CANONICAL_HOMOGENEOUS_OPEN",
    "FULL_SCALE_STRUT_CANONICAL_ZERO_INTERSECTION_RESOLVED",
    "FULL_SCALE_STRUT_CANONICAL_ONE_HOMOGENEOUS_RESOLVED",
}
check("the preregistered precision hierarchy assigns one verdict", outcome in allowed, outcome)

payload = {
    "prior_commit": PRIOR_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "levels": level_records,
    "nonhomogeneous": {
        "all_full_rank": nonhomogeneous_full_rank,
        "binary_disagreement": nonhomogeneous_disagreement,
        "minimum_singular_over_ball_radius": mp_string(minimum_nonhomogeneous_margin),
        "maximum_P100_P160_relative_change": mp_string(maximum_nonhomogeneous_relative_change),
        "maximum_source_target_corruption_margin": mp_string(max(source_target_margins)),
    },
    "homogeneous_status": homogeneous_status,
    "no_refit_validation": no_refit_records,
    "parity_candidate_projector_difference": mp_string(projector_difference),
    "pole_deletion": pole_deletion_records,
    "classification": {
        "nonhomogeneous_canonical_intersection": (
            "ZERO; PRIMARY MULTIPRECISION" if nonhomogeneous_full_rank else "OPEN"
        ),
        "homogeneous_canonical_intersection": homogeneous_status,
        "gauge_or_physical": "NOT CLASSIFIED",
        "propagation_tick_c_G_planck_mass": "NOT EVALUATED",
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
if passed != tests or outcome == "FULL_SCALE_STRUT_CANONICAL_PRECISION_CONTROL_FAILED":
    raise SystemExit(1)

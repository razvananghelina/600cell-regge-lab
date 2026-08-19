#!/usr/bin/env python3
"""Freeze the corrected geometric 120-column strut carrier before targets."""

from collections import Counter
import contextlib
import hashlib
import io
from itertools import combinations
import json
import math
from pathlib import Path
import runpy

import mpmath as mp
import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_corrected_strut_carrier.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_corrected_strut_carrier_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_corrected_strut_carrier_protocol.md"
PROTOCOL_CORRECTION = (
    ROOT / "docs/gravity/gravity_600cell_corrected_strut_carrier_protocol_correction.md"
)
LOCAL_SOURCE = HERE / "verify_gravity_600cell_local_data_lift.py"
LOCAL_ARTIFACT = HERE / "gravity_600cell_local_data_lift.json"
LOCAL_RESULT = ROOT / "docs/gravity/gravity_600cell_local_data_lift_result.md"
FRUSTUM_SOURCE = HERE / "verify_gravity_600cell_homothetic_frustum_equivalence.py"
FRUSTUM_ARTIFACT = HERE / "gravity_600cell_homothetic_frustum_equivalence.json"
FRUSTUM_RESULT = ROOT / "docs/gravity/gravity_600cell_homothetic_frustum_equivalence_result.md"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
TICK_ARTIFACT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"

PROTOCOL_COMMIT = "5e46f63"
CORRECTION_COMMIT = "67f4e0d"
EXPECTED_HASHES = {
    "prior_art": "e0064e73d161f7ba64b5a5c0c14ace4276cd0ce21d22fcb04b329506132064ca",
    "protocol": "73bf0666a43e61ebb0c0f425362396b032ab1b8026aafdc8d18c00337e4351e8",
    "protocol_correction": "c72884e9ceeb9020bcd4750bb00da9d3769bf105aec8d6dcda1eb7587e8dcd76",
    "local_source": "4389861a4b64d043325e0661ae9c2340f61e5c8eb50399c9fd2083a360dadbc1",
    "local_artifact": "0a569e48189c56bc081efcee33f7826fedd52afb93b6135ddb2fec385b56fbdf",
    "local_result": "646972a19450f1734ef522cb0b9693cc809b19d7895eb21823b20332a958d56d",
    "frustum_source": "99f47f0cfc70d2c0784d002cc08898e29f28a53e51930e6683c95629af128587",
    "frustum_artifact": "7e7c23efaf24a2c99a68f3b302b9ef575e0f777ef46f73ccaea9f99e1ecd58dc",
    "frustum_result": "b63808c260f12711ab25bdc72414f36c6c0f89f9420619c448c875d2dac7b093",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "tick_artifact": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
}
INPUT_PATHS = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "protocol_correction": PROTOCOL_CORRECTION,
    "local_source": LOCAL_SOURCE,
    "local_artifact": LOCAL_ARTIFACT,
    "local_result": LOCAL_RESULT,
    "frustum_source": FRUSTUM_SOURCE,
    "frustum_artifact": FRUSTUM_ARTIFACT,
    "frustum_result": FRUSTUM_RESULT,
    "geometry_source": GEOMETRY_SOURCE,
    "tick_artifact": TICK_ARTIFACT,
}

mp.mp.dps = 100
VERTEX_COUNT = 120
INTERNAL_COUNT = 840
FINAL_COUNT = 720
ROW_COUNT = INTERNAL_COUNT + FINAL_COUNT

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


def mp_string(value, digits=70):
    return mp.nstr(value, digits)


def symbolic_trapezoid_control():
    lam, a_square, a_parallel, edge_square = sy.symbols(
        "lambda a_square a_parallel edge_square"
    )
    delta = lam - 1
    strut_u = a_square
    strut_v = a_square + 2 * delta * a_parallel + delta**2 * edge_square
    diagonal = a_square + 2 * lam * a_parallel + lam**2 * edge_square
    variables = (a_square, a_parallel)
    differential = sy.simplify(
        delta * sy.Matrix([sy.diff(diagonal, item) for item in variables])
        + sy.Matrix([sy.diff(strut_u, item) for item in variables])
        - lam * sy.Matrix([sy.diff(strut_v, item) for item in variables])
    )
    natural_jacobian = sy.Matrix([
        [sy.diff(strut_u, item) for item in variables],
        [sy.diff(strut_v, item) for item in variables],
    ])
    generic_det = sy.factor(natural_jacobian.det())
    static_rank = natural_jacobian.subs(lam, 1).rank()
    return bool(
        differential == sy.zeros(2, 1)
        and sy.simplify(generic_det - 2 * (lam - 1)) == 0
        and static_rank == 1
    ), {
        "differential_residual": [str(sy.factor(value)) for value in differential],
        "natural_jacobian_determinant": str(generic_det),
        "static_natural_rank": static_rank,
        "derived_cross_variation": "(-delta_s_u+lambda*delta_s_v)/(lambda-1)",
    }


ETA = sy.diag(1, 1, 1, -1)
NORMAL = sy.Matrix((0, 0, 0, 1))
CANONICAL = tuple(sy.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
))
LOCAL_PAIRS = tuple(combinations(range(4), 2))


def exact_local_geometry(scale, lapse):
    top = tuple(scale * point + lapse * NORMAL for point in CANONICAL)
    jacobian = sy.zeros(10, 16)
    for row, (left, right) in enumerate(LOCAL_PAIRS):
        covector = 2 * (top[left] - top[right]).T * ETA
        jacobian[row, 4 * left:4 * left + 4] = covector
        jacobian[row, 4 * right:4 * right + 4] = -covector
    for vertex in range(4):
        covector = 2 * (top[vertex] - CANONICAL[vertex]).T * ETA
        jacobian[6 + vertex, 4 * vertex:4 * vertex + 4] = covector
    pivots = tuple(jacobian.rref()[1])
    right_inverse = sy.zeros(16, 10)
    inverse = jacobian[:, pivots].inv()
    for row, column in enumerate(pivots):
        right_inverse[column, :] = inverse[row, :]
    kernel_vectors = jacobian.nullspace()
    kernel = sy.Matrix.hstack(*kernel_vectors)
    data = sy.zeros(10, 8)
    for row, (left, right) in enumerate(LOCAL_PAIRS):
        data[row, left] = 8 * scale
        data[row, right] = 8 * scale
    for vertex in range(4):
        data[6 + vertex, 4 + vertex] = 1
    return top, jacobian, right_inverse, kernel, data, pivots


def rational_block_controls(local_artifact):
    records = {
        (int(record["scale"]), int(record["lapse"])): record
        for record in local_artifact["records"]
        if record["name"] == "baseline"
    }
    details = []
    all_ok = set(records) == {(2, 5), (3, 11)}
    for (scale, lapse), record in sorted(records.items()):
        top, jacobian, right_inverse, kernel, data, pivots = exact_local_geometry(
            sy.Integer(scale), sy.Integer(lapse)
        )
        block = sy.Matrix([
            [sy.Rational(value) for value in row] for row in record["block"]
        ])
        physical = sy.simplify(kernel * block + right_inverse * data)
        natural_ok = bool(jacobian * physical == data)
        cross_ok = True
        cross_checks = 0
        for lower in range(4):
            for upper in range(4):
                if lower == upper:
                    continue
                covector = 2 * (top[upper] - CANONICAL[lower]).T * ETA
                response = sy.simplify(
                    covector * physical[4 * upper:4 * upper + 4, 4:8]
                )
                expected = sy.zeros(1, 4)
                expected[0, lower] = -sy.Rational(1, scale - 1)
                expected[0, upper] += sy.Rational(scale, scale - 1)
                cross_ok &= bool(response == expected)
                cross_checks += 4
        pole_ok = True
        for vertex in range(4):
            covector = 2 * (top[vertex] - CANONICAL[vertex]).T * ETA
            response = sy.simplify(
                covector * physical[4 * vertex:4 * vertex + 4, 4:8]
            )
            expected = sy.zeros(1, 4)
            expected[0, vertex] = 1
            pole_ok &= bool(response == expected)
        record_ok = bool(
            natural_ok and cross_ok and pole_ok and len(pivots) == 10
            and kernel.shape == (16, 6)
        )
        all_ok &= record_ok
        details.append({
            "scale": scale,
            "lapse": lapse,
            "pivot_columns": list(pivots),
            "natural_data_exact": natural_ok,
            "ordered_cross_coefficient_checks": cross_checks,
            "cross_formula_exact": cross_ok,
            "pole_identity_exact": pole_ok,
        })
    return all_ok, details


def logical_edge(edge):
    return tuple(vertex % VERTEX_COUNT for vertex in edge)


def map_edge(edge, permutation):
    mapped = []
    for vertex in edge:
        sheet, logical = divmod(int(vertex), VERTEX_COUNT)
        mapped.append(int(permutation[logical]) + sheet * VERTEX_COUNT)
    return tuple(sorted(mapped))


def role_rows(internal_edges, final_edges):
    rows = {}
    for edge in internal_edges:
        lower, upper = edge
        if upper == lower + VERTEX_COUNT:
            rows[edge] = {lower: "one"}
        else:
            rows[edge] = {
                lower: "source_kappa",
                upper - VERTEX_COUNT: "target_minus_lambda_kappa",
            }
    for edge in final_edges:
        rows[edge] = {}
    return rows


def equivariance_control(rows, stabilizer):
    failures = []
    for group_index, permutation in enumerate(stabilizer):
        for edge, roles in rows.items():
            mapped_edge = map_edge(edge, permutation)
            mapped_roles = {
                int(permutation[column]): role for column, role in roles.items()
            }
            if mapped_edge not in rows or rows[mapped_edge] != mapped_roles:
                failures.append((group_index, edge, mapped_edge))
                break
    return not failures, failures[:3]


def corrupted_control(rows, stabilizer, first_diagonal):
    corrupted = {edge: dict(roles) for edge, roles in rows.items()}
    source = first_diagonal[0]
    del corrupted[first_diagonal][source]
    equivariant, failures = equivariance_control(corrupted, stabilizer)
    collective_failed = len(corrupted[first_diagonal]) == 1
    return bool(not equivariant and collective_failed), failures


def build_numeric_carrier(internal_edges, final_edges, lam, rho, q_diag):
    kappa = rho / ((lam - 1) * q_diag)
    matrix = mp.matrix(ROW_COUNT, VERTEX_COUNT)
    row_records = []
    for row, edge in enumerate(internal_edges):
        lower, upper = edge
        if upper == lower + VERTEX_COUNT:
            matrix[row, lower] = 1
            coefficients = [{"column": lower, "role": "one", "value": "1"}]
            kind = "pole"
        else:
            target = upper - VERTEX_COUNT
            matrix[row, lower] = kappa
            matrix[row, target] = -lam * kappa
            coefficients = [
                {"column": lower, "role": "source_kappa", "value": mp_string(kappa)},
                {"column": target, "role": "target_minus_lambda_kappa", "value": mp_string(-lam * kappa)},
            ]
            kind = "diagonal"
        row_records.append({"edge": list(edge), "kind": kind, "coefficients": coefficients})
    for offset, edge in enumerate(final_edges):
        row_records.append({"edge": list(edge), "kind": "new_boundary", "coefficients": []})
    return matrix, row_records, kappa


def householder_complement_mp(size):
    uniform = mp.matrix([1 / mp.sqrt(size) for _ in range(size)])
    e0 = mp.matrix(size, 1)
    e0[0] = 1
    vector = e0 - uniform
    norm_square = (vector.T * vector)[0]
    householder = mp.eye(size) - 2 * (vector * vector.T) / norm_square
    complement = householder[:, 1:size]
    return uniform, complement


def spectrum_census(matrix):
    gram = matrix.T * matrix
    eigenvalues_mp = list(mp.eigsy(gram, eigvals_only=True))
    eigenvalues_mp = sorted((max(mp.mpf(0), value) for value in eigenvalues_mp), reverse=True)
    singular_mp = [mp.sqrt(value) for value in eigenvalues_mp]

    uniform, complement = householder_complement_mp(VERTEX_COUNT)
    uniform_gain = mp.sqrt((uniform.T * gram * uniform)[0])
    coupling = complement.T * gram * uniform
    coupling_norm = mp.sqrt((coupling.T * coupling)[0])
    restricted_gram = complement.T * gram * complement
    restricted_eigen_mp = list(mp.eigsy(restricted_gram, eigvals_only=True))
    restricted_eigen_mp = sorted(
        (max(mp.mpf(0), value) for value in restricted_eigen_mp), reverse=True
    )
    restricted_gains_mp = [mp.sqrt(value) for value in restricted_eigen_mp]

    numpy_matrix = np.array(matrix.tolist(), dtype=float)
    gram64 = numpy_matrix.T @ numpy_matrix
    eigen64 = np.linalg.eigvalsh(gram64)[::-1]
    singular64 = np.sqrt(np.maximum(0.0, eigen64))
    uniform64 = np.full(VERTEX_COUNT, 1 / math.sqrt(VERTEX_COUNT))
    e064 = np.zeros(VERTEX_COUNT)
    e064[0] = 1
    vector64 = e064 - uniform64
    householder64 = np.eye(VERTEX_COUNT) - 2 * np.outer(vector64, vector64) / np.dot(vector64, vector64)
    complement64 = householder64[:, 1:]
    restricted64 = np.linalg.eigvalsh(complement64.T @ gram64 @ complement64)[::-1]
    restricted_gains64 = np.sqrt(np.maximum(0.0, restricted64))
    uniform_gain64 = float(np.linalg.norm(numpy_matrix @ uniform64))
    coupling64 = float(np.linalg.norm(complement64.T @ gram64 @ uniform64))

    singular_mp_float = np.asarray([float(value) for value in singular_mp])
    restricted_mp_float = np.asarray([float(value) for value in restricted_gains_mp])
    full_error = float(np.max(
        np.abs(singular64 - singular_mp_float)
        / np.maximum(1.0, np.abs(singular_mp_float))
    ))
    restricted_error = float(np.max(
        np.abs(restricted_gains64 - restricted_mp_float)
        / np.maximum(1.0, np.abs(restricted_mp_float))
    ))
    threshold64 = float(
        np.finfo(float).eps * max(matrix.rows, matrix.cols) * singular64[0]
    )
    threshold_mp = mp.mpf("1e-80") * singular_mp[0]
    rank64 = int(np.sum(singular64 > threshold64))
    rank_mp = sum(value > threshold_mp for value in singular_mp)
    return {
        "singular_values_mp": [mp_string(value) for value in singular_mp],
        "singular_values_binary64": [f"{value:.17e}" for value in singular64],
        "restricted_complement_gains_mp": [mp_string(value) for value in restricted_gains_mp],
        "restricted_complement_gains_binary64": [f"{value:.17e}" for value in restricted_gains64],
        "uniform_gain_mp": mp_string(uniform_gain),
        "uniform_gain_binary64": f"{uniform_gain64:.17e}",
        "uniform_complement_coupling_mp": mp_string(coupling_norm),
        "uniform_complement_coupling_binary64": f"{coupling64:.17e}",
        "condition_number_mp": mp_string(singular_mp[0] / singular_mp[-1]),
        "condition_number_binary64": f"{singular64[0] / singular64[-1]:.17e}",
        "rank_mp": rank_mp,
        "rank_binary64": rank64,
        "threshold_mp": mp_string(threshold_mp),
        "threshold_binary64": f"{threshold64:.17e}",
        "maximum_full_relative_discrepancy": f"{full_error:.17e}",
        "maximum_restricted_relative_discrepancy": f"{restricted_error:.17e}",
        "agreement": bool(
            rank_mp == rank64 == VERTEX_COUNT
            and full_error < 1e-10 and restricted_error < 1e-10
            and abs(uniform_gain64 - float(uniform_gain)) / max(1.0, float(uniform_gain)) < 1e-10
            and abs(coupling64 - float(coupling_norm)) / max(1.0, float(coupling_norm)) < 1e-10
        ),
    }


print("=" * 78)
print("TARGET-BLIND CORRECTED 600-CELL STRUT CARRIER")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUT_PATHS.items()}
local_artifact = json.loads(LOCAL_ARTIFACT.read_text())
frustum_artifact = json.loads(FRUSTUM_ARTIFACT.read_text())
tick = json.loads(TICK_ARTIFACT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and local_artifact.get("outcome") == "LOCAL_STAR_DATA_LIFT_DERIVED"
    and local_artifact.get("passed") == local_artifact.get("tests") == 13
    and frustum_artifact.get("outcome") == "HOMOTHETIC_SCHEDULES_ONE_LORENTZIAN_FRUSTUM"
    and frustum_artifact.get("passed") == frustum_artifact.get("tests")
    and tick.get("outcome") == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
)
check("all target-blind inputs retain exact frozen provenance", provenance_ok, str(hashes))

symbolic_ok, symbolic_record = symbolic_trapezoid_control()
check("the planar-trapezoid differential identity is exact and rejects the static inversion", symbolic_ok, str(symbolic_record))

rational_ok, rational_records = rational_block_controls(local_artifact)
check("both independent rational face-gluing blocks reproduce every cross and pole response", rational_ok, str(rational_records))

with contextlib.redirect_stdout(io.StringIO()):
    geometry_namespace = runpy.run_path(str(GEOMETRY_SOURCE))
geometry_import_ok = geometry_namespace["tests"] == geometry_namespace["passed"] == 43
check("the directly reconstructed one-slab geometry retains all 43 controls", geometry_import_ok)

models = geometry_namespace["models"]
geometry_counts_ok = True
geometry_records = {}
for parity in ("even", "odd"):
    model = models[parity]
    internal_edges = tuple(sorted(tuple(map(int, edge)) for edge in model["internal_edges"]))
    final_edges = tuple(sorted(tuple(map(int, edge)) for edge in model["new_edges"]))
    poles = tuple(edge for edge in internal_edges if edge[1] == edge[0] + VERTEX_COUNT)
    diagonals = tuple(edge for edge in internal_edges if edge not in poles)
    count_ok = bool(
        len(model["old_edges"]) == 720
        and len(final_edges) == 720
        and len(diagonals) == 720
        and len(poles) == 120
        and len(model["stabilizer"]) == 24
        and all(0 <= edge[0] < 120 <= edge[1] < 240 for edge in internal_edges)
    )
    geometry_counts_ok &= count_ok
    geometry_records[parity] = {
        "internal_edges": internal_edges,
        "final_edges": final_edges,
        "poles": poles,
        "diagonals": diagonals,
        "stabilizer": tuple(model["stabilizer"]),
    }
check("both schedules have the exact 720+120 internal and 720 final edge census", geometry_counts_ok)

states_equal = tick["solutions"]["even"]["state"] == tick["solutions"]["odd"]["state"]
s_value, r_value = map(mp.mpf, tick["solutions"]["even"]["state"])
M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
RHO0 = mp.mpf("0.0102")**2
lam = mp.exp(s_value)
rho = RHO0 * mp.exp(r_value)
q_diag = lam * L0_SQUARE - rho
background_ok = bool(states_equal and lam > 0 and lam != 1 and rho > 0 and q_diag > 0)
check(
    "the accepted non-static background lies in the disclosed carrier domain",
    background_ok,
    f"lambda={mp_string(lam,40)}, rho={mp_string(rho,40)}, q_diag={mp_string(q_diag,40)}",
)

parity_records = {}
all_incidence_ok = True
all_equivariance_ok = True
all_collective_ok = True
all_corruption_ok = True
all_numeric_ok = True
for parity in ("even", "odd"):
    item = geometry_records[parity]
    internal_edges = item["internal_edges"]
    final_edges = item["final_edges"]
    rows = role_rows(internal_edges, final_edges)
    support_histogram = Counter(len(roles) for roles in rows.values())
    column_support = Counter()
    for roles in rows.values():
        for column in roles:
            column_support[column] += 1
    pole_rows = [edge for edge in internal_edges if edge[1] == edge[0] + VERTEX_COUNT]
    pole_identity = all(rows[edge] == {edge[0]: "one"} for edge in pole_rows)
    incidence_ok = bool(
        support_histogram == Counter({2: 720, 1: 120, 0: 720})
        and set(column_support) == set(range(VERTEX_COUNT))
        and set(column_support.values()) == {13}
        and pole_identity
    )
    all_incidence_ok &= incidence_ok

    equivariant, equivariance_failures = equivariance_control(rows, item["stabilizer"])
    all_equivariance_ok &= equivariant

    first_diagonal = min(item["diagonals"])
    corruption_ok, corruption_failures = corrupted_control(
        rows, item["stabilizer"], first_diagonal
    )
    all_corruption_ok &= corruption_ok

    matrix, serialized_rows, kappa = build_numeric_carrier(
        internal_edges, final_edges, lam, rho, q_diag
    )
    collective = matrix * mp.matrix([1 for _ in range(VERTEX_COUNT)])
    expected_collective = mp.matrix(ROW_COUNT, 1)
    for row, edge in enumerate(internal_edges):
        expected_collective[row] = 1 if edge[1] == edge[0] + VERTEX_COUNT else -rho / q_diag
    collective_error = max(abs(collective[row] - expected_collective[row]) for row in range(ROW_COUNT))
    collective_ok = bool(collective_error < mp.mpf("1e-90"))
    all_collective_ok &= collective_ok

    print(f"[{parity}] target-blind 120-column intrinsic spectrum", flush=True)
    census = spectrum_census(matrix)
    all_numeric_ok &= census["agreement"]
    parity_records[parity] = {
        "internal_edge_order": [list(edge) for edge in internal_edges],
        "final_edge_order": [list(edge) for edge in final_edges],
        "rows": serialized_rows,
        "support_histogram": {str(key): value for key, value in sorted(support_histogram.items())},
        "column_support_histogram": {str(key): value for key, value in sorted(Counter(column_support.values()).items())},
        "literal_pole_identity_rank": VERTEX_COUNT if pole_identity else 0,
        "equivariance": {
            "group_elements_checked": len(item["stabilizer"]),
            "exact": equivariant,
            "first_failures": [list(map(str, failure)) for failure in equivariance_failures],
        },
        "collective": {
            "maximum_error": mp_string(collective_error),
            "diagonal_value": mp_string(-rho / q_diag),
            "exact_by_coefficient_identity": True,
        },
        "corruption": {
            "edge": list(first_diagonal),
            "rejected": corruption_ok,
            "first_equivariance_failures": [list(map(str, failure)) for failure in corruption_failures],
        },
        "kappa": mp_string(kappa),
        "spectrum": census,
    }

check("every carrier row and column has the exact preregistered incidence support", all_incidence_ok)
check("all 24 stabilizer elements intertwine both corrected carriers exactly", all_equivariance_ok)
check("both corrected carriers sum to the analytic collective lapse column", all_collective_ok)
check("the one-row endpoint corruption breaks collective identity and equivariance", all_corruption_ok)
check("100-digit and binary64 target-blind carrier censuses agree", all_numeric_ok)

even_singular = np.asarray([
    float(value) for value in parity_records["even"]["spectrum"]["singular_values_mp"]
])
odd_singular = np.asarray([
    float(value) for value in parity_records["odd"]["spectrum"]["singular_values_mp"]
])
parity_distance = float(np.max(
    np.abs(even_singular - odd_singular) / np.maximum(1.0, np.abs(even_singular))
))
finite_census_ok = bool(np.isfinite(parity_distance))
check("the complete target-blind parity census is finite and serialized", finite_census_ok, f"ordered-spectrum distance={parity_distance:.3e}")

firewall_ok = set(hashes) == set(INPUT_PATHS) and all(path.exists() for path in INPUT_PATHS.values())
check("the target firewall loaded only the eleven preregistered geometric inputs", firewall_ok, str(sorted(hashes)))

controls_ok = bool(
    provenance_ok and symbolic_ok and rational_ok and geometry_import_ok
    and geometry_counts_ok and background_ok and all_incidence_ok
    and all_equivariance_ok and all_collective_ok and all_corruption_ok
    and finite_census_ok and firewall_ok
)
if not controls_ok:
    outcome = "CORRECTED_STRUT_CARRIER_CONTROL_FAILED"
elif not all_numeric_ok:
    outcome = "CORRECTED_STRUT_CARRIER_NUMERICALLY_OPEN"
else:
    outcome = "CORRECTED_STRUT_CARRIER_FROZEN"

payload = {
    "prior_art_commit": "f419238",
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": CORRECTION_COMMIT,
    "input_sha256": hashes,
    "loaded_inputs_only": sorted(hashes),
    "target_artifacts_loaded": False,
    "candidate_count_per_parity": 1,
    "symbolic_trapezoid": symbolic_record,
    "rational_face_gluing_controls": rational_records,
    "background": {
        "state": [str(s_value), str(r_value)],
        "lambda": mp_string(lam),
        "rho": mp_string(rho),
        "L0_square": mp_string(L0_SQUARE),
        "q_diagonal": mp_string(q_diag),
        "kappa_formula": "rho/((lambda-1)*(lambda*L0_square-rho))",
    },
    "parities": parity_records,
    "ordered_parity_spectrum_distance": f"{parity_distance:.17e}",
    "classification": {
        "carrier_formula": "DERIVED EXACT FOR lambda != 1",
        "accepted_background_matrix": "DERIVED COMPUTATIONAL TARGET-BLIND",
        "dynamic_alignment": "NOT EVALUATED",
        "gauge_or_physical_interpretation": "OPEN",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"TOTAL: {passed}/{tests} tests PASSED")
if outcome not in {
    "CORRECTED_STRUT_CARRIER_FROZEN",
    "CORRECTED_STRUT_CARRIER_NUMERICALLY_OPEN",
} or passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Target-disclosed exact audit of the full scale--strut carrier."""

from collections import Counter, deque
import contextlib
from hashlib import sha256
import io
from itertools import combinations
import json
import math
from pathlib import Path
import runpy

import mpmath as mp
import numpy as np
from scipy import linalg as sla
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_full_scale_strut_carrier.json"
DISCLOSURE = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_carrier_exploratory_disclosure.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_carrier_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_carrier_protocol.md"
LOCAL_SOURCE = HERE / "verify_gravity_600cell_local_data_lift.py"
LOCAL_JSON = HERE / "gravity_600cell_local_data_lift.json"
LOCAL_RESULT = ROOT / "docs/gravity/gravity_600cell_local_data_lift_result.md"
STRUT_SOURCE = HERE / "verify_gravity_600cell_corrected_strut_carrier.py"
STRUT_JSON = HERE / "gravity_600cell_corrected_strut_carrier.json"
TWO_SOURCE = HERE / "verify_gravity_600cell_two_frustum_face_gluing_adversarial.py"
TWO_JSON = HERE / "gravity_600cell_two_frustum_face_gluing_adversarial.json"
TWO_RESULT = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_result.md"
FACE_SOURCE = HERE / "verify_gravity_600cell_canonical_data_admissibility.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
TICK_JSON = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"

PROTOCOL_COMMIT = "19dd489"
EXPECTED_HASHES = {
    "disclosure": "e3dba59118e35cc2370beec4b081cc18fdfd1753fda74a6ca2b0a013d86bd473",
    "prior_art": "3fc6c3e75ad92c3c20bb420d97e26e73fdd62b69a9aac44162fa95c74c29219a",
    "protocol": "3017d856f838dcc885a916c9b3eca265b2ae949a134fe176b94872466cad32fb",
    "local_source": "4389861a4b64d043325e0661ae9c2340f61e5c8eb50399c9fd2083a360dadbc1",
    "local_json": "0a569e48189c56bc081efcee33f7826fedd52afb93b6135ddb2fec385b56fbdf",
    "local_result": "646972a19450f1734ef522cb0b9693cc809b19d7895eb21823b20332a958d56d",
    "strut_source": "80f0a17960adee496fe7d51678ea99849280ecd3fca6254efc8acd3753aad348",
    "strut_json": "e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7",
    "two_source": "b7a1f63e193aad50783929c8448ce99c18f1b50dc8e5ea27e3ed1102ec9dfa26",
    "two_json": "0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542",
    "two_result": "b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3",
    "face_source": "4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "tick_json": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
}
INPUTS = {
    "disclosure": DISCLOSURE,
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "local_source": LOCAL_SOURCE,
    "local_json": LOCAL_JSON,
    "local_result": LOCAL_RESULT,
    "strut_source": STRUT_SOURCE,
    "strut_json": STRUT_JSON,
    "two_source": TWO_SOURCE,
    "two_json": TWO_JSON,
    "two_result": TWO_RESULT,
    "face_source": FACE_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "tick_json": TICK_JSON,
}

VERTEX_COUNT = 120
EDGE_COUNT = 720
CELL_FLEX_COUNT = 3600
LOCAL_UNKNOWN_COUNT = 48
CONSTANT_COLUMN = 48
ROW_COUNT = 1560
DATA_COUNT = 240
PRIME = 1000003
NEW_REPRESENTATIVES = ((4, 7), (5, 13), (7, 17))
OLD_REPRESENTATIVES = ((2, 5), (3, 11))

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
CANONICAL = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
))

mp.mp.dps = 100
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


def mp_string(value, digits=70):
    return mp.nstr(value, digits)


def rational(value):
    value = sp.cancel(value)
    if value and value.is_Rational is not True:
        raise TypeError(f"non-rational coefficient: {value}")
    return value


def add_value(row, column, value):
    if not value:
        return
    updated = row.get(column, 0) + value
    if updated:
        row[column] = updated
    elif column in row:
        del row[column]


def boundary_data_row(original, scale, edges):
    """Convert global upper-edge/strut columns to 240 vertex data."""
    result = {}
    strut_start = CELL_FLEX_COUNT + EDGE_COUNT
    edge_factor = sp.Integer(8) * sp.Integer(scale)
    for column, raw in original.items():
        value = rational(raw)
        if column < CELL_FLEX_COUNT:
            continue
        if column < strut_start:
            left, right = edges[column - CELL_FLEX_COUNT]
            add_value(result, left, value * edge_factor)
            add_value(result, right, value * edge_factor)
        else:
            add_value(result, VERTEX_COUNT + column - strut_start, value)
    return result


def universal_affine_rows(built, scale, tetrahedra, edges):
    """Substitute one universal 6x8 block into all 600-cell face rows."""
    constraints = []
    metadata = []
    for global_row, original in enumerate(built["rows"]):
        by_data = {}
        for data_column, value in boundary_data_row(original, scale, edges).items():
            by_data.setdefault(data_column, {})[CONSTANT_COLUMN] = value
        for flex_column, raw in original.items():
            if flex_column >= CELL_FLEX_COUNT:
                continue
            coefficient = rational(raw)
            cell_index, local_flex = divmod(flex_column, 6)
            tetrahedron = tetrahedra[cell_index]
            for local_vertex, global_vertex in enumerate(tetrahedron):
                row = by_data.setdefault(global_vertex, {})
                add_value(row, 8 * local_flex + local_vertex, coefficient)
                row = by_data.setdefault(VERTEX_COUNT + global_vertex, {})
                add_value(row, 8 * local_flex + 4 + local_vertex, coefficient)
        for data_column in sorted(by_data):
            row = by_data[data_column]
            if row:
                constraints.append(row)
                metadata.append((global_row, data_column))
    return constraints, metadata


def rational_mod(value, prime=PRIME):
    value = sp.Rational(value)
    numerator = int(value.p) % prime
    denominator = int(value.q) % prime
    if denominator == 0:
        raise ZeroDivisionError("denominator vanishes modulo selection prime")
    return numerator * pow(denominator, -1, prime) % prime


def select_independent_rows(rows, unknown_count=LOCAL_UNKNOWN_COUNT):
    """Use modular elimination only to select independent coefficient rows."""
    pivots = {}
    selected = []
    for row_index, original in enumerate(rows):
        row = {
            column: rational_mod(value)
            for column, value in original.items()
            if column < unknown_count and rational_mod(value)
        }
        while row:
            pivot = min(row)
            if pivot not in pivots:
                inverse = pow(row[pivot], -1, PRIME)
                row = {
                    column: value * inverse % PRIME
                    for column, value in row.items()
                    if value * inverse % PRIME
                }
                pivots[pivot] = row
                selected.append(row_index)
                break
            factor = row[pivot]
            for column, value in pivots[pivot].items():
                updated = (row.get(column, 0) - factor * value) % PRIME
                if updated:
                    row[column] = updated
                elif column in row:
                    del row[column]
        if len(selected) == unknown_count:
            break
    return selected, tuple(sorted(pivots))


def exact_selected_solve(rows):
    selected, modular_pivots = select_independent_rows(rows)
    if len(selected) != LOCAL_UNKNOWN_COUNT:
        return {
            "selected": selected,
            "modular_pivots": modular_pivots,
            "matrix": None,
            "determinant": 0,
            "solution": None,
            "residual_count": None,
            "first_residual": None,
        }
    matrix = sp.Matrix([
        [rational(rows[index].get(column, 0)) for column in range(LOCAL_UNKNOWN_COUNT)]
        for index in selected
    ])
    rhs = sp.Matrix([
        -rational(rows[index].get(CONSTANT_COLUMN, 0)) for index in selected
    ])
    determinant = matrix.det(method="domain-ge")
    solution = matrix.inv(method="DM") * rhs if determinant else None
    residual_count = 0
    first_residual = None
    if solution is not None:
        for row_index, row in enumerate(rows):
            residual = rational(
                sum(
                    value * solution[column]
                    for column, value in row.items()
                    if column < LOCAL_UNKNOWN_COUNT
                ) + row.get(CONSTANT_COLUMN, 0)
            )
            if residual:
                residual_count += 1
                if first_residual is None:
                    first_residual = (row_index, residual)
    return {
        "selected": selected,
        "modular_pivots": modular_pivots,
        "matrix": matrix,
        "determinant": determinant,
        "solution": solution,
        "residual_count": residual_count,
        "first_residual": first_residual,
    }


def data_matrix(scale):
    result = sp.zeros(10, 8)
    for row, (left, right) in enumerate(combinations(range(4), 2)):
        result[row, left] = 8 * scale
        result[row, right] = 8 * scale
    for vertex in range(4):
        result[6 + vertex, 4 + vertex] = 1
    return result


def cross_jacobian(scale, lapse):
    top = tuple(scale * point + lapse * NORMAL for point in CANONICAL)
    ordered = tuple((left, right) for left in range(4) for right in range(4) if left != right)
    result = sp.zeros(12, 16)
    for row, (lower, upper) in enumerate(ordered):
        covector = 2 * (top[upper] - CANONICAL[lower]).T * ETA
        result[row, 4 * upper:4 * upper + 4] = covector
    return ordered, result


def expected_cross(scale, lapse):
    scale = sp.Rational(scale)
    lapse = sp.Rational(lapse)
    denominator = (scale - 1) ** 2
    A = 6 - 2 * lapse**2 / denominator
    B = 2 + 2 * lapse**2 / denominator
    C = -1 / (scale - 1)
    D = scale / (scale - 1)
    ordered, _ = cross_jacobian(scale, lapse)
    result = sp.zeros(12, 8)
    for row, (source, target) in enumerate(ordered):
        result[row, source] = A
        result[row, target] = B
        result[row, 4 + source] = C
        result[row, 4 + target] = D
    return result, (A, B, C, D)


def block_digest(block):
    hasher = sha256()
    for row in range(block.rows):
        for column in range(block.cols):
            value = sp.Rational(block[row, column])
            hasher.update(f"{row},{column},{int(value.p)},{int(value.q)}\n".encode())
    return hasher.hexdigest()


def physical_response(built, solution, scale):
    block = sp.Matrix(6, 8, lambda row, column: solution[8 * row + column])
    geometry = built["local_geometry"]
    physical = sp.simplify(
        geometry["kernel"] * block
        + geometry["right_inverse"] * data_matrix(sp.Rational(scale))
    )
    ordered, cross = cross_jacobian(sp.Rational(scale), sp.Rational(built["_lapse"]))
    return block, physical, ordered, sp.simplify(cross * physical)


def local_underdetermination_control(face_namespace):
    scale, lapse = map(sp.Integer, NEW_REPRESENTATIVES[0])
    geometry = face_namespace["local_length_geometry"](
        scale, lapse, ETA, CANONICAL
    )
    natural = geometry["jacobian"]
    _, cross = cross_jacobian(scale, lapse)
    complete = natural.col_join(cross)
    left_vectors = complete.T.nullspace()
    left = sp.Matrix.hstack(*left_vectors)
    A, B, C, D = sp.symbols("A B C D")
    response = sp.zeros(12, 8)
    ordered = tuple((i, j) for i in range(4) for j in range(4) if i != j)
    for row, (source, target) in enumerate(ordered):
        response[row, source] = A
        response[row, target] = B
        response[row, 4 + source] = C
        response[row, 4 + target] = D
    rhs = data_matrix(scale).col_join(response)
    residual = sp.simplify(left.T * rhs)
    coefficient_rows = []
    variables = (A, B, C, D)
    for value in residual:
        value = sp.expand(value)
        if not value:
            continue
        coefficients = [sp.diff(value, variable) for variable in variables]
        constant = value.subs({variable: 0 for variable in variables})
        coefficient_rows.append(coefficients + [constant])
    observed = sp.Matrix(coefficient_rows)
    expected = sp.Matrix([
        [1, 1, 0, 0, -8],
        [0, 0, 1, 1, -1],
    ])
    same = bool(
        complete.rank() == 16
        and left.shape == (22, 6)
        and observed.rank() == 2
        and observed.col_join(expected).rank() == 2
    )
    return same, {
        "representative": [int(scale), int(lapse)],
        "complete_length_jacobian_rank": complete.rank(),
        "left_nullity": left.cols,
        "compatibility_rank": observed.rank(),
        "compatibility_rowspace": [
            [str(value) for value in row] for row in expected.tolist()
        ],
    }


def graph_controls(edges):
    neighbours = [set() for _ in range(VERTEX_COUNT)]
    for left, right in edges:
        neighbours[left].add(right)
        neighbours[right].add(left)
    seen = {0}
    queue = deque([0])
    while queue:
        vertex = queue.popleft()
        for target in neighbours[vertex]:
            if target not in seen:
                seen.add(target)
                queue.append(target)
    triangle = any(
        right in neighbours[left]
        for vertex in range(VERTEX_COUNT)
        for left in neighbours[vertex]
        for right in neighbours[vertex]
        if left < right
    )
    return bool(
        len(edges) == EDGE_COUNT
        and len(seen) == VERTEX_COUNT
        and set(map(len, neighbours)) == {12}
        and triangle
    )


def modular_matrix_rank(rows, column_count, prime=PRIME):
    pivots = {}
    for original in rows:
        row = {column: value % prime for column, value in original.items() if value % prime}
        while row:
            pivot = min(row)
            if pivot not in pivots:
                inverse = pow(row[pivot], -1, prime)
                pivots[pivot] = {
                    column: value * inverse % prime
                    for column, value in row.items()
                    if value * inverse % prime
                }
                break
            factor = row[pivot]
            for column, value in pivots[pivot].items():
                updated = (row.get(column, 0) - factor * value) % prime
                if updated:
                    row[column] = updated
                elif column in row:
                    del row[column]
    return len(pivots)


def map_edge(edge, action):
    mapped = []
    for vertex in edge:
        mapped.append(int(action[vertex] if vertex < VERTEX_COUNT else action[vertex - VERTEX_COUNT] + VERTEX_COUNT))
    return tuple(sorted(mapped))


def map_column(column, action):
    if column < VERTEX_COUNT:
        return int(action[column])
    return VERTEX_COUNT + int(action[column - VERTEX_COUNT])


def role_equivariance(rows, stabilizer):
    failures = []
    for group_index, action in enumerate(stabilizer):
        for edge, roles in rows.items():
            mapped_edge = map_edge(edge, action)
            mapped_roles = {map_column(column, action): role for column, role in roles.items()}
            if mapped_edge not in rows or rows[mapped_edge] != mapped_roles:
                failures.append((group_index, edge, mapped_edge))
                break
    return not failures, failures[:3]


def accepted_carrier(internal_edges, final_edges, lam, rho, L0_square):
    q_diag = lam * L0_square - rho
    A = -16 * rho / (L0_square * (lam - 1) ** 2)
    B = 8 + 16 * rho / (L0_square * (lam - 1) ** 2)
    scale_factor = L0_square / (8 * q_diag)
    kappa = rho / ((lam - 1) * q_diag)
    matrix = np.zeros((ROW_COUNT, DATA_COUNT), dtype=float)
    roles = {}
    for row, edge in enumerate(internal_edges):
        lower, upper = edge
        if upper == lower + VERTEX_COUNT:
            matrix[row, VERTEX_COUNT + lower] = 1.0
            roles[edge] = {VERTEX_COUNT + lower: "pole_one"}
        else:
            target = upper - VERTEX_COUNT
            matrix[row, lower] = float(scale_factor * A)
            matrix[row, target] = float(scale_factor * B)
            matrix[row, VERTEX_COUNT + lower] = float(kappa)
            matrix[row, VERTEX_COUNT + target] = float(-lam * kappa)
            roles[edge] = {
                lower: "scale_source_A",
                target: "scale_target_B",
                VERTEX_COUNT + lower: "strut_source_C",
                VERTEX_COUNT + target: "strut_target_D",
            }
    for offset, edge in enumerate(final_edges):
        left = edge[0] - VERTEX_COUNT
        right = edge[1] - VERTEX_COUNT
        row = len(internal_edges) + offset
        matrix[row, left] = float(1 / lam)
        matrix[row, right] = float(1 / lam)
        roles[edge] = {left: "upper_scale", right: "upper_scale"}
    return matrix, roles, {
        "A": A,
        "B": B,
        "scale_factor": scale_factor,
        "kappa": kappa,
        "q_diag": q_diag,
    }


def spectral_diagnostic(matrix):
    direct = sla.svdvals(matrix, overwrite_a=False, check_finite=True)
    gram_values = np.linalg.eigvalsh(matrix.T @ matrix)[::-1]
    minimum_gram = float(np.min(gram_values))
    gram = np.sqrt(np.maximum(0.0, gram_values))
    discrepancy = float(np.max(np.abs(direct - gram) / np.maximum(1.0, np.abs(direct))))
    finite = bool(np.all(np.isfinite(direct)) and np.all(np.isfinite(gram)))
    return {
        "direct_singular_values": [f"{value:.17e}" for value in direct],
        "gram_singular_values": [f"{value:.17e}" for value in gram],
        "direct_condition_number": f"{direct[0] / direct[-1]:.17e}",
        "gram_condition_number": f"{gram[0] / gram[-1]:.17e}" if gram[-1] else "inf",
        "minimum_raw_gram_eigenvalue": f"{minimum_gram:.17e}",
        "maximum_relative_discrepancy": f"{discrepancy:.17e}",
        "finite": finite,
        "agreement": bool(finite and discrepancy < 1e-8 and gram[-1] > 0),
    }


print("=" * 78)
print("TARGET-DISCLOSED FULL 600-CELL SCALE--STRUT CARRIER AUDIT")
print("=" * 78)

hashes_before = {name: digest(path) for name, path in INPUTS.items()}
provenance_ok = hashes_before == EXPECTED_HASHES
check("all target-disclosed inputs retain frozen provenance", provenance_ok, str(hashes_before))

local_artifact = json.loads(LOCAL_JSON.read_text())
strut_artifact = json.loads(STRUT_JSON.read_text())
two_artifact = json.loads(TWO_JSON.read_text())
tick = json.loads(TICK_JSON.read_text())
upstream_ok = bool(
    local_artifact.get("outcome") == "LOCAL_STAR_DATA_LIFT_DERIVED"
    and local_artifact.get("passed") == local_artifact.get("tests") == 13
    and strut_artifact.get("outcome") == "CORRECTED_STRUT_CARRIER_FROZEN"
    and strut_artifact.get("passed") == strut_artifact.get("tests")
    and two_artifact.get("outcome") == "ADVERSARIAL_TWO_FRUSTUM_DIAGONAL_ONLY"
    and two_artifact.get("passed") == two_artifact.get("tests") == 11
    and all(record["full_poincare_face_stabilizer_dimension"] == 1 for record in two_artifact["records"])
    and tick.get("outcome") == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
)
check("all frozen local, strut, face and background controls persist", upstream_ok)

print("reconstructing the frozen exact face-equation namespace", flush=True)
with contextlib.redirect_stdout(io.StringIO()):
    face_namespace = runpy.run_path(str(FACE_SOURCE))
face_source_ok = bool(
    face_namespace["tests"] == face_namespace["passed"]
    and digest(FACE_SOURCE) == EXPECTED_HASHES["face_source"]
)
check("the exact face-equation source reruns without changing its frozen bytes", face_source_ok)

complex_data = face_namespace["complex_data"]
tetrahedra = tuple(complex_data["tetrahedra"])
edges = tuple(sorted(complex_data["edges"]))
complex_ok = bool(
    len(complex_data["vertices"]) == 120
    and len(edges) == 720
    and len(complex_data["triangles"]) == 1200
    and len(tetrahedra) == 600
    and graph_controls(edges)
)
check("the exact 600-cell carrier has f=(120,720,1200,600) and an odd cycle", complex_ok)

local_ok, local_record = local_underdetermination_control(face_namespace)
check("one frustum imposes exactly A+B=8 and C+D=1, not local uniqueness", local_ok, str(local_record))

old_records = {
    (int(record["scale"]), int(record["lapse"])): record
    for record in local_artifact["records"] if record["name"] == "baseline"
}
old_ok = set(old_records) == set(OLD_REPRESENTATIVES)
old_details = []
for representative in OLD_REPRESENTATIVES:
    scale, lapse = map(sp.Integer, representative)
    built = face_namespace["baseline_builds"][representative]
    built["_lapse"] = lapse
    block = sp.Matrix([[sp.Rational(value) for value in row] for row in old_records[representative]["block"]])
    physical = sp.simplify(
        built["local_geometry"]["kernel"] * block
        + built["local_geometry"]["right_inverse"] * data_matrix(scale)
    )
    _, cross = cross_jacobian(scale, lapse)
    observed = sp.simplify(cross * physical)
    expected, coefficients = expected_cross(scale, lapse)
    match = bool(observed == expected)
    old_ok &= match
    old_details.append({
        "representative": list(representative),
        "match": match,
        "coefficients": [str(value) for value in coefficients],
        "block_sha256": block_digest(block),
    })
check("the two old rational blocks reproduce the disclosed physical response", old_ok, str(old_details))

new_records = []
new_controls_ok = True
finite_formula_agrees = True
corrupted_rejected = True
for scale_integer, lapse_integer in NEW_REPRESENTATIVES:
    print(f"building new exact global representative ({scale_integer},{lapse_integer})", flush=True)
    scale = sp.Integer(scale_integer)
    lapse = sp.Integer(lapse_integer)
    built = face_namespace["build_global"](
        complex_data, scale, lapse,
        face_namespace["ETA"], face_namespace["CANONICAL_BASE"],
    )
    built["_lapse"] = lapse
    local_controls = bool(
        built["local_geometry"]["controls"]
        and built["transition_control"]
        and built["inverse_control"]
        and built["face_control"]
        and built["local_fixed_ranks"] == {5: 1200}
    )
    rows, metadata = universal_affine_rows(built, scale, tetrahedra, edges)
    solved = exact_selected_solve(rows)
    exact_controls = bool(
        len(rows) == 51320
        and len(solved["selected"]) == LOCAL_UNKNOWN_COUNT
        and solved["modular_pivots"] == tuple(range(LOCAL_UNKNOWN_COUNT))
        and solved["determinant"] != 0
        and solved["solution"] is not None
        and solved["residual_count"] == 0
        and solved["first_residual"] is None
    )
    match = False
    mismatch_count = None
    corrupted_count = None
    block = None
    observed = None
    coefficients = None
    if exact_controls:
        block, physical, ordered, observed = physical_response(built, solved["solution"], scale)
        expected, coefficients = expected_cross(scale, lapse)
        mismatch_count = sum(observed[row, column] != expected[row, column] for row in range(12) for column in range(8))
        corrupted = expected.copy()
        for row, (_, target) in enumerate(ordered):
            corrupted[row, 4 + target] += 1
        corrupted_count = sum(observed[row, column] != corrupted[row, column] for row in range(12) for column in range(8))
        match = mismatch_count == 0
        corrupted_rejected &= corrupted_count > 0
    new_controls_ok &= local_controls and exact_controls
    finite_formula_agrees &= match
    determinant = solved["determinant"]
    new_records.append({
        "representative": [scale_integer, lapse_integer],
        "nondegenerate": bool(
            scale != 1 and lapse != 0
            and (scale - 1) ** 2 - 3 * lapse**2 != 0
        ),
        "global_face_controls": local_controls,
        "affine_constraint_count": len(rows),
        "modular_selection_prime": PRIME,
        "selected_row_indices": solved["selected"],
        "selected_metadata": [list(metadata[index]) for index in solved["selected"]],
        "modular_pivots": list(solved["modular_pivots"]),
        "exact_determinant": str(determinant),
        "exact_determinant_sha256": sha256(str(determinant).encode()).hexdigest(),
        "full_residual_nonzero_count": solved["residual_count"],
        "first_residual": None if solved["first_residual"] is None else [solved["first_residual"][0], str(solved["first_residual"][1])],
        "block_sha256": block_digest(block) if block is not None else None,
        "candidate_coefficients": [str(value) for value in coefficients] if coefficients is not None else None,
        "cross_response": [[str(observed[row, column]) for column in range(8)] for row in range(12)] if observed is not None else None,
        "candidate_mismatch_count": mismatch_count,
        "corrupted_D_plus_one_mismatch_count": corrupted_count,
    })

check("all three new complete global systems have exact rank-48 solutions and zero residual", new_controls_ok)
check("the finite candidate comparison and D+1 corruption are both classified", bool(corrupted_rejected and all(record["candidate_mismatch_count"] is not None for record in new_records)), f"candidate_agrees={finite_formula_agrees}")

print("reconstructing both frozen staircase schedules", flush=True)
with contextlib.redirect_stdout(io.StringIO()):
    geometry_namespace = runpy.run_path(str(GEOMETRY_SOURCE))
models = geometry_namespace["models"]
schedule_ok = bool(
    geometry_namespace["tests"] == geometry_namespace["passed"] == 43
    and all(len(models[parity]["stabilizer"]) == 24 for parity in ("even", "odd"))
)
check("both staircase schedules retain their exact order-24 stabilizers", schedule_ok)

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
check("the accepted curved slab lies in the disclosed nonstatic carrier domain", background_ok, f"lambda={mp_string(lam,40)}, rho={mp_string(rho,40)}")

parity_records = {}
all_incidence_ok = True
all_rank_ok = True
all_equivariance_ok = True
all_collective_ok = True
all_corruption_ok = True
all_numeric_finite = True
all_numeric_agreement = True
for parity in ("even", "odd"):
    corrected = strut_artifact["parities"][parity]
    internal_edges = tuple(tuple(map(int, edge)) for edge in corrected["internal_edge_order"])
    final_edges = tuple(tuple(map(int, edge)) for edge in corrected["final_edge_order"])
    model = models[parity]
    model_internal = tuple(sorted(tuple(map(int, edge)) for edge in model["internal_edges"]))
    model_final = tuple(sorted(tuple(map(int, edge)) for edge in model["new_edges"]))
    schedule_edges_match = internal_edges == model_internal and final_edges == model_final
    matrix, roles, coefficients = accepted_carrier(internal_edges, final_edges, lam, rho, L0_SQUARE)

    support_histogram = Counter(map(len, roles.values()))
    column_supports = Counter()
    for row_roles in roles.values():
        for column in row_roles:
            column_supports[column] += 1
    scale_supports = Counter(column_supports[column] for column in range(VERTEX_COUNT))
    strut_supports = Counter(column_supports[VERTEX_COUNT + column] for column in range(VERTEX_COUNT))
    incidence_ok = bool(
        schedule_edges_match
        and len(internal_edges) == 840
        and len(final_edges) == 720
        and support_histogram == Counter({4: 720, 1: 120, 2: 720})
        and scale_supports == Counter({24: 120})
        and strut_supports == Counter({13: 120})
    )
    all_incidence_ok &= incidence_ok

    upper_graph = tuple(sorted((edge[0] - VERTEX_COUNT, edge[1] - VERTEX_COUNT) for edge in final_edges))
    incidence_rows = [{left: 1, right: 1} for left, right in upper_graph]
    incidence_rank = modular_matrix_rank(incidence_rows, VERTEX_COUNT)
    poles = tuple(edge for edge in internal_edges if edge[1] == edge[0] + VERTEX_COUNT)
    exact_rank = bool(
        graph_controls(upper_graph)
        and incidence_rank == VERTEX_COUNT
        and len(poles) == VERTEX_COUNT
        and all(roles[edge] == {VERTEX_COUNT + edge[0]: "pole_one"} for edge in poles)
    )
    all_rank_ok &= exact_rank

    equivariant, equivariance_failures = role_equivariance(roles, model["stabilizer"])
    all_equivariance_ok &= equivariant

    scale_diagonal = coefficients["scale_factor"] * (coefficients["A"] + coefficients["B"])
    strut_diagonal = coefficients["kappa"] * (1 - lam)
    collective_errors = {
        "scale_diagonal": abs(scale_diagonal - L0_SQUARE / q_diag),
        "scale_upper": abs(2 / lam - 2 / lam),
        "strut_diagonal": abs(strut_diagonal + rho / q_diag),
        "strut_pole": mp.mpf(0),
    }
    collective_ok = max(collective_errors.values()) < mp.mpf("1e-90")
    all_collective_ok &= collective_ok

    first_diagonal = min(edge for edge in internal_edges if edge not in poles)
    corrupted_roles = {edge: dict(value) for edge, value in roles.items()}
    del corrupted_roles[first_diagonal][first_diagonal[0]]
    corrupted_equivariant, corrupted_failures = role_equivariance(corrupted_roles, model["stabilizer"])
    corruption_ok = not corrupted_equivariant and len(corrupted_roles[first_diagonal]) == 3
    all_corruption_ok &= corruption_ok

    diagnostic = spectral_diagnostic(matrix)
    all_numeric_finite &= diagnostic["finite"]
    all_numeric_agreement &= diagnostic["agreement"]
    parity_records[parity] = {
        "internal_edge_order": [list(edge) for edge in internal_edges],
        "final_edge_order": [list(edge) for edge in final_edges],
        "support_histogram": {str(key): value for key, value in sorted(support_histogram.items())},
        "scale_column_support_histogram": {str(key): value for key, value in sorted(scale_supports.items())},
        "strut_column_support_histogram": {str(key): value for key, value in sorted(strut_supports.items())},
        "upper_graph_connected_nonbipartite": graph_controls(upper_graph),
        "upper_unsigned_incidence_modular_rank": incidence_rank,
        "literal_pole_identity_rank": VERTEX_COUNT if exact_rank else 0,
        "complete_exact_rank": DATA_COUNT if exact_rank else "OPEN",
        "equivariance": {
            "group_elements_checked": len(model["stabilizer"]),
            "exact": equivariant,
            "first_failures": [list(map(str, value)) for value in equivariance_failures],
        },
        "collective_errors": {key: mp_string(value) for key, value in collective_errors.items()},
        "corruption": {
            "edge": list(first_diagonal),
            "rejected": corruption_ok,
            "first_failures": [list(map(str, value)) for value in corrupted_failures],
        },
        "coefficient_values": {key: mp_string(value) for key, value in coefficients.items()},
        "spectral_diagnostic": diagnostic,
    }

check("both accepted-background carriers have the exact preregistered supports", all_incidence_ok)
check("pole identity plus the nonbipartite upper graph proves exact rank 240", all_rank_ok)
check("all 24 schedule symmetries preserve both full carrier coefficient roles", all_equivariance_ok)
check("both homogeneous scale and homogeneous strut derivatives are exact", all_collective_ok)
check("a one-endpoint corruption is rejected by schedule equivariance", all_corruption_ok)
check("both target-blind spectral diagnostics are finite and classified", all_numeric_finite, f"agreement={all_numeric_agreement}")

hashes_after = {name: digest(path) for name, path in INPUTS.items()}
firewall_ok = bool(
    hashes_after == hashes_before
    and set(hashes_before) == set(INPUTS)
    and not any("hessian" in str(path).lower() or "strong" in str(path).lower() for path in INPUTS.values())
)
check("the target firewall loaded no Hessian, strong-equation or sector-target artifact", firewall_ok)

controls_ok = bool(
    provenance_ok and upstream_ok and face_source_ok and complex_ok and local_ok
    and old_ok and new_controls_ok and corrupted_rejected and schedule_ok
    and background_ok
)
curved_ok = bool(
    all_incidence_ok and all_rank_ok and all_equivariance_ok
    and all_collective_ok and all_corruption_ok
)
if not controls_ok or not firewall_ok:
    outcome = "FULL_SCALE_STRUT_CONTROL_FAILED"
elif not finite_formula_agrees:
    outcome = "FULL_SCALE_STRUT_FINITE_DISAGREEMENT"
elif not curved_ok:
    outcome = "FULL_SCALE_STRUT_CURVED_CARRIER_FAILED"
elif not all_numeric_finite or not all_numeric_agreement:
    outcome = "FULL_SCALE_STRUT_NUMERICALLY_OPEN"
else:
    outcome = "FULL_SCALE_STRUT_FINITE_CONTROLS_CORROBORATE"

allowed = {
    "FULL_SCALE_STRUT_CONTROL_FAILED",
    "FULL_SCALE_STRUT_FINITE_DISAGREEMENT",
    "FULL_SCALE_STRUT_CURVED_CARRIER_FAILED",
    "FULL_SCALE_STRUT_FINITE_CONTROLS_CORROBORATE",
    "FULL_SCALE_STRUT_NUMERICALLY_OPEN",
}
check("the preregistered target-disclosed hierarchy assigns one outcome", outcome in allowed, outcome)

payload = {
    "disclosure_and_prior_art_commit": "bfce559",
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes_before,
    "source_sha256": digest(Path(__file__)),
    "loaded_inputs_only": sorted(hashes_before),
    "action_hessian_or_sector_target_loaded": False,
    "candidate_count_per_parity": 1,
    "local_underdetermination": local_record,
    "old_controls": old_details,
    "new_exact_global_controls": new_records,
    "finite_formula_agrees": finite_formula_agrees,
    "background": {
        "state": [str(s_value), str(r_value)],
        "lambda": mp_string(lam),
        "rho": mp_string(rho),
        "L0_square": mp_string(L0_SQUARE),
        "q_diagonal": mp_string(q_diag),
    },
    "parities": parity_records,
    "classification": {
        "generic_formula": "OPEN PENDING SYMBOLIC ADVERSARIAL REPLICATION",
        "three_new_rational_controls": "DERIVED EXACT",
        "accepted_background_matrix": "TARGET-BLIND KINEMATIC CANDIDATE",
        "complete_rank": "DERIVED EXACT BY BLOCK-ROW ARGUMENT" if curved_ok else "OPEN",
        "gauge_constraint_or_physical_interpretation": "NOT EVALUATED",
        "tick_c_G_Planck_particle_mass": "NOT EVALUATED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"TOTAL: {passed}/{tests} tests PASSED")
if passed != tests or outcome == "FULL_SCALE_STRUT_CONTROL_FAILED":
    raise SystemExit(1)

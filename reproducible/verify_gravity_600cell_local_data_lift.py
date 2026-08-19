#!/usr/bin/env python3
"""Derive the exact universal 6x8 local canonical-data lift."""

from collections import Counter, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import runpy

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_local_data_lift.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_local_data_lift_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_local_data_lift_protocol.md"
GLOBAL_LIFT_SOURCE = HERE / "verify_gravity_600cell_rational_data_lift.py"
GLOBAL_LIFT_JSON = HERE / "gravity_600cell_rational_data_lift.json"
GLOBAL_LIFT_RESULT = ROOT / "docs/gravity/gravity_600cell_rational_data_lift_result.md"
ADMISSIBILITY_SOURCE = HERE / "verify_gravity_600cell_canonical_data_admissibility.py"
REFUTED_SOURCE = HERE / "verify_gravity_600cell_canonical_data_carrier.py"
REFUTED_JSON = HERE / "gravity_600cell_canonical_data_carrier.json"

PROTOCOL_COMMIT = "dd302d8"
EXPECTED_HASHES = {
    "prior_art": "91a98568926afa6c556b143050a081354bcce28bc02da673d8863aa1aadc8aa7",
    "protocol": "a8e6f4f60c3d3688c4d41a789c16f8290d2bbe6b7f6b1395dcdd2751df31407b",
    "global_lift_source": "65a097cd11dea830fd16bad988cd6d1b88ce4b84e0700b7dfaf7477a2c198ecb",
    "global_lift_json": "1b6ac46a0ea4889f476cc71d51ac464c27caa6d4b6a9b2f2d74ff93da77b123f",
    "global_lift_result": "c69e367fed93498705a30058134000bd77be4845589e7942a90b059e53aa3ecc",
    "admissibility_source": "4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e",
    "refuted_source": "3adb80448e19fd99f0b8ec205497f11325b0a6f7a72c9a2785f9b65778707750",
    "refuted_json": "4065950aaac4180ec1cdd0b82f7a8bc403b2969c50d26cf14cc28592085cb2c5",
}

FLEX_COUNT = 3600
EDGE_COUNT = 720
VERTEX_COUNT = 120
DATA_COUNT = 240
LOCAL_UNKNOWNS = 48
CONSTANT_COLUMN = 48

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


def add_scaled(target, source, factor):
    if not factor:
        return
    for column, value in source.items():
        add_value(target, column, factor * value)


def scaled(source, factor):
    return {
        column: factor * value
        for column, value in source.items()
        if factor * value
    }


def graph_controls(vertex_count, edges, tetrahedra):
    neighbours = [set() for _ in range(vertex_count)]
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
    triangle = tuple(tetrahedra[0][:3])
    triangle_edges = {
        tuple(sorted((triangle[index], triangle[(index + 1) % 3])))
        for index in range(3)
    }
    return bool(
        len(seen) == vertex_count
        and triangle_edges <= set(edges)
        and all(len(value) == 12 for value in neighbours)
    )


def boundary_data_row(original, scale, edges, corrupt=False):
    data = {}
    strut_start = FLEX_COUNT + EDGE_COUNT
    edge_factor = sp.Integer(8) * sp.Integer(scale)
    for column, raw_value in original.items():
        value = rational(raw_value)
        if column < FLEX_COUNT:
            continue
        if column < strut_start:
            edge_index = column - FLEX_COUNT
            left, right = edges[edge_index]
            add_value(data, left, value * edge_factor)
            if not (corrupt and edge_index == 0):
                add_value(data, right, value * edge_factor)
        else:
            vertex = column - strut_start
            add_value(data, VERTEX_COUNT + vertex, value)
    return data


def universal_affine_rows(built, scale, tetrahedra, edges, corrupt=False):
    """Substitute one universal 6x8 block into all global residuals."""
    constraints = []
    metadata = []
    for global_row, original in enumerate(built["rows"]):
        by_data = {}
        boundary = boundary_data_row(original, scale, edges, corrupt=corrupt)
        for data_column, value in boundary.items():
            by_data.setdefault(data_column, {})[CONSTANT_COLUMN] = value
        for flex_column, raw_value in original.items():
            if flex_column >= FLEX_COUNT:
                continue
            coefficient = rational(raw_value)
            cell_index, local_flex = divmod(flex_column, 6)
            tetrahedron = tetrahedra[cell_index]
            for local_vertex, global_vertex in enumerate(tetrahedron):
                scale_row = by_data.setdefault(global_vertex, {})
                add_value(
                    scale_row,
                    8 * local_flex + local_vertex,
                    coefficient,
                )
                strut_row = by_data.setdefault(
                    VERTEX_COUNT + global_vertex, {}
                )
                add_value(
                    strut_row,
                    8 * local_flex + 4 + local_vertex,
                    coefficient,
                )
        for data_column in sorted(by_data):
            row = by_data[data_column]
            if row:
                constraints.append(row)
                metadata.append((global_row, data_column))
    return constraints, metadata


def solve_affine(rows, metadata, unknown_count):
    """Exact sparse affine elimination, independent of the global solver."""
    pivots = {}
    pivot_metadata = {}
    obstructions = []
    maximum_width = 0
    for row_index, original in enumerate(rows):
        row = dict(original)
        while True:
            unknowns = [column for column in row if column < unknown_count]
            if not unknowns:
                if row.get(CONSTANT_COLUMN, 0):
                    obstructions.append((metadata[row_index], dict(row)))
                break
            pivot = min(unknowns)
            if pivot not in pivots:
                inverse = sp.Integer(1) / row[pivot]
                row = scaled(row, inverse)
                pivots[pivot] = row
                pivot_metadata[pivot] = metadata[row_index]
                maximum_width = max(maximum_width, len(row))
                break
            add_scaled(row, pivots[pivot], -row[pivot])
    full_rank = len(pivots) == unknown_count and set(pivots) == set(range(unknown_count))
    solution = None
    if full_rank and not obstructions:
        values = [None] * unknown_count
        for pivot in range(unknown_count - 1, -1, -1):
            row = pivots[pivot]
            value = -row.get(CONSTANT_COLUMN, 0)
            for column, coefficient in row.items():
                if column == pivot or column == CONSTANT_COLUMN:
                    continue
                if column <= pivot or values[column] is None:
                    raise RuntimeError("affine pivots are not echelon ordered")
                value -= coefficient * values[column]
            values[pivot] = rational(value)
        solution = values
    return {
        "pivots": pivots,
        "pivot_metadata": pivot_metadata,
        "rank": len(pivots),
        "full_rank": full_rank,
        "obstructions": obstructions,
        "maximum_width": maximum_width,
        "solution": solution,
    }


def block_matrix(solution):
    return sp.Matrix(6, 8, lambda row, column: solution[8 * row + column])


def block_record(block):
    return [[str(block[row, column]) for column in range(8)] for row in range(6)]


def block_digest(block):
    hasher = sha256()
    for row in range(6):
        for column in range(8):
            value = sp.Rational(block[row, column])
            hasher.update(f"{row},{column},{int(value.p)},{int(value.q)}\n".encode())
    return hasher.hexdigest()


def local_star_lift(block, tetrahedra):
    lift = []
    for tetrahedron in tetrahedra:
        for local_flex in range(6):
            row = {}
            for local_vertex, global_vertex in enumerate(tetrahedron):
                add_value(row, global_vertex, block[local_flex, local_vertex])
                add_value(
                    row,
                    VERTEX_COUNT + global_vertex,
                    block[local_flex, 4 + local_vertex],
                )
            lift.append(row)
    return lift


def direct_residual(built, scale, edges, lift):
    count = 0
    first = None
    maximum_width = 0
    for row_index, original in enumerate(built["rows"]):
        residual = boundary_data_row(original, scale, edges, corrupt=False)
        for flex_column, raw_value in original.items():
            if flex_column < FLEX_COUNT:
                add_scaled(residual, lift[flex_column], rational(raw_value))
        if residual:
            count += 1
            maximum_width = max(maximum_width, len(residual))
            if first is None:
                first = (row_index, residual)
    return count, first, maximum_width


def exact_item_record(item):
    if item is None:
        return None
    location, row = item
    return {
        "location": list(location) if isinstance(location, tuple) else location,
        "coefficients": {
            str(column): str(value) for column, value in sorted(row.items())
        },
    }


def support_controls(lift, tetrahedra):
    stars = [set() for _ in range(VERTEX_COUNT)]
    for cell_index, tetrahedron in enumerate(tetrahedra):
        for vertex in tetrahedron:
            stars[vertex].add(cell_index)
    flex_exact = True
    data_cells = [set() for _ in range(DATA_COUNT)]
    data_flex = [0] * DATA_COUNT
    for flex_index, row in enumerate(lift):
        tetrahedron = tetrahedra[flex_index // 6]
        expected = set(tetrahedron) | {
            VERTEX_COUNT + vertex for vertex in tetrahedron
        }
        flex_exact &= set(row) == expected and len(row) == 8
        for data_column in row:
            data_cells[data_column].add(flex_index // 6)
            data_flex[data_column] += 1
    data_exact = all(
        data_cells[vertex] == stars[vertex]
        and data_cells[VERTEX_COUNT + vertex] == stars[vertex]
        and data_flex[vertex] == 120
        and data_flex[VERTEX_COUNT + vertex] == 120
        for vertex in range(VERTEX_COUNT)
    )
    return {
        "exact_flex_support": flex_exact,
        "exact_data_star_support": data_exact,
        "vertex_star_histogram": {
            str(key): value for key, value in sorted(Counter(map(len, stars)).items())
        },
        "data_cell_support_histogram": {
            str(key): value
            for key, value in sorted(Counter(map(len, data_cells)).items())
        },
        "data_flex_support_histogram": {
            str(key): value
            for key, value in sorted(Counter(data_flex).items())
        },
    }


def local_data_block(scale):
    data = sp.zeros(10, 8)
    factor = sp.Integer(8) * sp.Integer(scale)
    for row, (left, right) in enumerate(combinations(range(4), 2)):
        data[row, left] = factor
        data[row, right] = factor
    for vertex in range(4):
        data[6 + vertex, 4 + vertex] = 1
    return data


def radial_physical_block(canonical, normal):
    result = sp.zeros(16, 8)
    for vertex in range(4):
        result[4 * vertex : 4 * vertex + 4, vertex] = canonical[vertex]
        result[4 * vertex : 4 * vertex + 4, 4 + vertex] = normal
    return result


def construction_ok(built):
    return bool(
        built["local_geometry"]["controls"]
        and built["transition_control"]
        and built["inverse_control"]
        and built["face_control"]
        and built["local_fixed_ranks"] == {5: 1200}
    )


def solve_construction(name, scale, lapse, built, tetrahedra, edges, namespace, use_alternate=False):
    print(f"universal local lift {name} ({scale},{lapse})", flush=True)
    rows, metadata = universal_affine_rows(
        built, scale, tetrahedra, edges, corrupt=False
    )
    exact = solve_affine(rows, metadata, LOCAL_UNKNOWNS)
    bad_rows, bad_metadata = universal_affine_rows(
        built, scale, tetrahedra, edges, corrupt=True
    )
    corrupted = solve_affine(bad_rows, bad_metadata, LOCAL_UNKNOWNS)

    block = None
    lift = None
    residual = (None, None, None)
    support = None
    local_reconciliation = None
    physical = None
    if exact["solution"] is not None:
        block = block_matrix(exact["solution"])
        lift = local_star_lift(block, tetrahedra)
        residual = direct_residual(built, scale, edges, lift)
        support = support_controls(lift, tetrahedra)

        geometry = built["local_geometry"]
        kernel = geometry["kernel"]
        right_inverse = geometry[
            "alternate_right_inverse" if use_alternate else "right_inverse"
        ]
        data = local_data_block(scale)
        radial_seed = radial_physical_block(
            namespace["CANONICAL_BASE"] if name != "odd_relabelling" else (
                namespace["CANONICAL_BASE"][1],
                namespace["CANONICAL_BASE"][0],
                namespace["CANONICAL_BASE"][2],
                namespace["CANONICAL_BASE"][3],
            ),
            namespace["NORMAL"],
        )
        if geometry["jacobian"] * radial_seed == data:
            radial_convention_factor = 1
        elif geometry["jacobian"] * radial_seed == -data:
            radial_convention_factor = -1
        else:
            raise RuntimeError("the old radial ansatz does not reproduce the declared data up to metric sign")
        radial = radial_convention_factor * radial_seed
        old_block = namespace["solve_full_column"](
            kernel, radial - right_inverse * data
        )
        physical = kernel * block + right_inverse * data
        delta = block - old_block
        local_reconciliation = {
            "local_data_from_old_radial_is_exact": bool(
                geometry["jacobian"] * radial == data
            ),
            "radial_metric_convention_factor": radial_convention_factor,
            "new_local_data_is_exact": bool(
                geometry["jacobian"] * physical == data
            ),
            "new_differs_from_old": bool(block != old_block),
            "difference_is_kernel_correction": bool(
                physical - radial == kernel * delta
                and geometry["jacobian"] * (physical - radial)
                == sp.zeros(10, 8)
            ),
            "old_block": block_record(old_block),
            "difference_rank": delta.rank(),
            "difference_nonzero_entries": sum(
                delta[row, column] != 0
                for row in range(6) for column in range(8)
            ),
        }

    record = {
        "name": name,
        "scale": scale,
        "lapse": lapse,
        "construction_controls": construction_ok(built),
        "affine_constraint_count": len(rows),
        "affine_rank": exact["rank"],
        "affine_maximum_width": exact["maximum_width"],
        "candidate_consistent": not exact["obstructions"],
        "candidate_obstruction_count": len(exact["obstructions"]),
        "candidate_first_obstruction": exact_item_record(
            exact["obstructions"][0] if exact["obstructions"] else None
        ),
        "block_all_nonzero": bool(
            block is not None
            and all(block[row, column] != 0 for row in range(6) for column in range(8))
        ),
        "block_sha256": block_digest(block) if block is not None else None,
        "block": block_record(block) if block is not None else None,
        "direct_residual_nonzero_rows": residual[0],
        "direct_first_residual": exact_item_record(residual[1]),
        "direct_residual_maximum_width": residual[2],
        "support": support,
        "corrupted_affine_rank": corrupted["rank"],
        "corrupted_obstruction_count": len(corrupted["obstructions"]),
        "corrupted_first_obstruction": exact_item_record(
            corrupted["obstructions"][0]
            if corrupted["obstructions"] else None
        ),
        "old_block_reconciliation": local_reconciliation,
    }
    return record, block, physical


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "global_lift_source": GLOBAL_LIFT_SOURCE,
    "global_lift_json": GLOBAL_LIFT_JSON,
    "global_lift_result": GLOBAL_LIFT_RESULT,
    "admissibility_source": ADMISSIBILITY_SOURCE,
    "refuted_source": REFUTED_SOURCE,
    "refuted_json": REFUTED_JSON,
}
hashes_before = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes_before == EXPECTED_HASHES
check("all universal-local-lift inputs have frozen provenance", provenance_ok, str(hashes_before))

global_lift = json.loads(GLOBAL_LIFT_JSON.read_text())
refuted = json.loads(REFUTED_JSON.read_text())
upstream_ok = bool(
    global_lift["outcome"] == "RATIONAL_VERTEX_STRUT_DATA_LIFT_DERIVED"
    and global_lift["passed"] == global_lift["tests"] == 12
    and refuted["outcome"] == "CANONICAL_DATA_VERTEX_CARRIER_REFUTED"
    and refuted["passed"] == refuted["tests"] == 10
)
check("the exact global positive and old local negative remain frozen", upstream_ok)

# Solver control with the exact solution (1,2) and one inconsistent affine row.
synthetic_rows = [
    {0: 1, CONSTANT_COLUMN: -1},
    {1: 1, CONSTANT_COLUMN: -2},
    {0: 1, 1: 1, CONSTANT_COLUMN: -3},
]
synthetic_metadata = [(0, 0), (1, 0), (2, 0)]
synthetic = solve_affine(synthetic_rows, synthetic_metadata, 2)
synthetic_bad = solve_affine(
    synthetic_rows + [{0: 1, CONSTANT_COLUMN: -2}],
    synthetic_metadata + [(3, 0)],
    2,
)
synthetic_ok = bool(
    synthetic["rank"] == 2
    and synthetic["solution"] == [1, 2]
    and not synthetic["obstructions"]
    and synthetic_bad["obstructions"]
)
check("the independent 48-style affine solver passes its known control", synthetic_ok)

print("reconstructing complete exact face equations", flush=True)
namespace = runpy.run_path(str(ADMISSIBILITY_SOURCE))
hashes_after = {name: digest(path) for name, path in paths.items()}
source_reproduction_ok = hashes_after == EXPECTED_HASHES
check("the complete exact face-equation source preserves every frozen byte", source_reproduction_ok)

complex_data = namespace["complex_data"]
edges = tuple(sorted(complex_data["edges"]))
tetrahedra = tuple(complex_data["tetrahedra"])
vertex_stars = Counter(
    vertex for tetrahedron in tetrahedra for vertex in tetrahedron
)
combinatorics_ok = bool(
    len(complex_data["vertices"]) == VERTEX_COUNT
    and len(edges) == EDGE_COUNT
    and len(complex_data["triangles"]) == 1200
    and len(tetrahedra) == 600
    and set(vertex_stars.values()) == {20}
    and graph_controls(VERTEX_COUNT, edges, tetrahedra)
)
check("the 600-cell graph and its 20-cell vertex stars are exact", combinatorics_ok)

constructions = []
for scale, lapse in namespace["REPRESENTATIVES"]:
    constructions.append((
        "baseline", scale, lapse,
        namespace["baseline_builds"][(scale, lapse)], False,
    ))
for scale, lapse in namespace["REPRESENTATIVES"]:
    constructions.append((
        "alternate_right_inverse", scale, lapse,
        namespace["build_global"](
            complex_data, scale, lapse, namespace["ETA"],
            namespace["CANONICAL_BASE"], use_alternate=True,
        ), True,
    ))
scale, lapse = namespace["REPRESENTATIVES"][0]
constructions.append((
    "reverse_faces", scale, lapse,
    namespace["build_global"](
        complex_data, scale, lapse, namespace["ETA"],
        namespace["CANONICAL_BASE"], reverse_orientation=True,
    ), False,
))
odd_canonical = (
    namespace["CANONICAL_BASE"][1], namespace["CANONICAL_BASE"][0],
    namespace["CANONICAL_BASE"][2], namespace["CANONICAL_BASE"][3],
)
constructions.append((
    "odd_relabelling", scale, lapse,
    namespace["build_global"](
        complex_data, scale, lapse, namespace["ETA"], odd_canonical,
    ), False,
))
constructions.append((
    "metric_sign", scale, lapse,
    namespace["build_global"](
        complex_data, scale, lapse, -namespace["ETA"],
        namespace["CANONICAL_BASE"],
    ), False,
))

records = []
blocks = {}
physical_responses = {}
for name, scale, lapse, built, use_alternate in constructions:
    record, block, physical = solve_construction(
        name, scale, lapse, built, tetrahedra, edges, namespace,
        use_alternate=use_alternate,
    )
    records.append(record)
    blocks[(name, scale, lapse)] = block
    physical_responses[(name, scale, lapse)] = physical

construction_controls_ok = all(
    record["construction_controls"] for record in records
)
check("all seven local-block constructions retain every exact face control", construction_controls_ok)

rank_consistency_decisions = {
    (record["affine_rank"] == LOCAL_UNKNOWNS, record["candidate_consistent"])
    for record in records
}
decision_agreement_ok = len(rank_consistency_decisions) == 1
check("all constructions agree on rank-48 universal-block consistency", decision_agreement_ok, str(rank_consistency_decisions))
candidate_positive = bool(
    decision_agreement_ok and rank_consistency_decisions == {(True, True)}
)

blocks_ok = all(
    record["affine_rank"] == LOCAL_UNKNOWNS
    and record["candidate_consistent"]
    and record["block_all_nonzero"]
    and record["direct_residual_nonzero_rows"] == 0
    and record["direct_first_residual"] is None
    for record in records
)
check(
    "every accepted unique all-nonzero 6x8 block vanishes on all 6000 rows",
    not candidate_positive or blocks_ok,
)

supports_ok = all(
    record["support"]["exact_flex_support"]
    and record["support"]["exact_data_star_support"]
    and record["support"]["vertex_star_histogram"] == {"20": 120}
    and record["support"]["data_cell_support_histogram"] == {"20": 240}
    and record["support"]["data_flex_support_histogram"] == {"120": 240}
    for record in records if record["support"] is not None
) and all(record["support"] is not None for record in records)
check(
    "every accepted support set equals its exact four-vertex cell or 20-cell star",
    not candidate_positive or supports_ok,
)

corruption_ok = all(
    record["corrupted_obstruction_count"] > 0
    and record["corrupted_first_obstruction"] is not None
    for record in records
)
check("the one-row incidence corruption is affinely inconsistent throughout", corruption_ok)

reconciliation_ok = all(
    value is not None
    and value["local_data_from_old_radial_is_exact"]
    and value["new_local_data_is_exact"]
    and value["new_differs_from_old"]
    and value["difference_is_kernel_correction"]
    for value in (record["old_block_reconciliation"] for record in records)
)
check(
    "every accepted new block is a length-preserving Poincare repair of the old block",
    not candidate_positive or reconciliation_ok,
)

graph_physical_ok = all(
    physical_responses[("baseline", scale, lapse)]
    == physical_responses[("alternate_right_inverse", scale, lapse)]
    and blocks[("baseline", scale, lapse)]
    != blocks[("alternate_right_inverse", scale, lapse)]
    for scale, lapse in namespace["REPRESENTATIVES"]
)
check(
    "accepted physical responses agree across right-inverse graphs while flex blocks differ",
    not candidate_positive or graph_physical_ok,
)

base_controls_ok = bool(
    provenance_ok and upstream_ok and synthetic_ok and source_reproduction_ok
    and combinatorics_ok and construction_controls_ok and corruption_ok
)
positive_decision = rank_consistency_decisions == {(True, True)}
positive_controls_ok = bool(
    not positive_decision
    or (blocks_ok and supports_ok and reconciliation_ok and graph_physical_ok)
)

if not base_controls_ok or not positive_controls_ok:
    outcome = "LOCAL_DATA_LIFT_CONTROL_FAILED"
elif not decision_agreement_ok:
    outcome = "LOCAL_DATA_LIFT_DISAGREEMENT_OPEN"
elif not positive_decision:
    outcome = "UNIVERSAL_LOCAL_DATA_LIFT_REFUTED"
else:
    outcome = "LOCAL_STAR_DATA_LIFT_DERIVED"

allowed = {
    "LOCAL_DATA_LIFT_CONTROL_FAILED",
    "LOCAL_DATA_LIFT_DISAGREEMENT_OPEN",
    "UNIVERSAL_LOCAL_DATA_LIFT_REFUTED",
    "LOCAL_STAR_DATA_LIFT_DERIVED",
}
check("the preregistered local-lift hierarchy assigns one outcome", outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes_before,
    "source_sha256": digest(Path(__file__)),
    "records": records,
    "classification": {
        "universal_local_6x8_block": (
            "DERIVED" if outcome == "LOCAL_STAR_DATA_LIFT_DERIVED" else "OPEN"
        ),
        "exact_vertex_star_support": (
            "DERIVED" if outcome == "LOCAL_STAR_DATA_LIFT_DERIVED" else "OPEN"
        ),
        "old_local_formula": "DERIVED NEGATIVE; RECONCILED BY POINCARE CORRECTION",
        "action_or_physical_dynamics": "NOT EVALUATED",
        "external_novelty": "OPEN",
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
        f"{record['name']} ({record['scale']},{record['lapse']}): "
        f"rank={record['affine_rank']}, "
        f"consistent={record['candidate_consistent']}, "
        f"all_nonzero={record['block_all_nonzero']}, "
        f"residual_rows={record['direct_residual_nonzero_rows']}, "
        f"bad_obstructions={record['corrupted_obstruction_count']}, "
        f"block={record['block_sha256']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

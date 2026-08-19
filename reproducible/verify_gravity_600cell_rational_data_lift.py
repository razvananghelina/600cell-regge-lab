#!/usr/bin/env python3
"""Construct and verify the exact rational global canonical-data lift."""

from collections import Counter, deque
from hashlib import sha256
import json
from pathlib import Path
import runpy

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_rational_data_lift.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_rational_data_lift_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_rational_data_lift_protocol.md"
IDENTIFICATION_SOURCE = (
    HERE / "verify_gravity_600cell_canonical_data_space_identification.py"
)
IDENTIFICATION_JSON = (
    HERE / "gravity_600cell_canonical_data_space_identification.json"
)
ADMISSIBILITY_SOURCE = HERE / "verify_gravity_600cell_canonical_data_admissibility.py"
ADMISSIBILITY_JSON = HERE / "gravity_600cell_canonical_data_admissibility.json"
REFUTED_SOURCE = HERE / "verify_gravity_600cell_canonical_data_carrier.py"
REFUTED_JSON = HERE / "gravity_600cell_canonical_data_carrier.json"

PROTOCOL_COMMIT = "1cbe2e1"
EXPECTED_HASHES = {
    "prior_art": "5a1a695a4526a70320f9316771f015d3c14bffdf3bca3d115a52fa8e9fe73d27",
    "protocol": "264c0ea6dff510481f1c4b09221596cca7fe3239924a312675b681d0a5497398",
    "identification_source": "21022ab2014f5c95bd3e2f06bf1137f713533df465d152a40c069f107dd153ae",
    "identification_json": "3db0b9ce8c90cba9de3fbbff818129388d79a98e0483a0ca3ae53b2e4d271434",
    "admissibility_source": "4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e",
    "admissibility_json": "fa45c80739ca0dda4f82c9da98a4b22f4d8a18c182a40696a2a22d1d26ec89a1",
    "refuted_source": "3adb80448e19fd99f0b8ec205497f11325b0a6f7a72c9a2785f9b65778707750",
    "refuted_json": "4065950aaac4180ec1cdd0b82f7a8bc403b2969c50d26cf14cc28592085cb2c5",
}

FLEX_COUNT = 3600
EDGE_COUNT = 720
VERTEX_COUNT = 120
DATA_COUNT = 240

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
        raise TypeError(f"non-rational exact coefficient: {value}")
    return value


def add_scaled(target, source, factor):
    if not factor:
        return
    for column, value in source.items():
        updated = target.get(column, 0) + factor * value
        if updated:
            target[column] = updated
        elif column in target:
            del target[column]


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


def incidence_rows(edges, corrupt=False):
    rows = []
    for edge_index, (left, right) in enumerate(edges):
        row = {left: sp.Integer(1), right: sp.Integer(1)}
        if corrupt and edge_index == 0:
            row.pop(right)
        rows.append(row)
    return rows


def sparse_matrix(rows, columns):
    entries = {
        (row_index, column): value
        for row_index, row in enumerate(rows)
        for column, value in row.items()
    }
    return sp.SparseMatrix(len(rows), columns, entries)


def split_exact_rows(built, scale, edges, corrupt=False):
    """Build F and [E U,S] directly, without modular composition."""
    flex_rows = []
    data_rows = []
    strut_start = FLEX_COUNT + EDGE_COUNT
    edge_factor = sp.Integer(8) * sp.Integer(scale)
    for original in built["rows"]:
        flex = {}
        data = {}
        for column, raw_value in original.items():
            value = rational(raw_value)
            if column < FLEX_COUNT:
                flex[column] = value
            elif column < strut_start:
                edge_index = column - FLEX_COUNT
                left, right = edges[edge_index]
                add_scaled(data, {left: sp.Integer(1)}, value * edge_factor)
                if not (corrupt and edge_index == 0):
                    add_scaled(data, {right: sp.Integer(1)}, value * edge_factor)
            else:
                vertex = column - strut_start
                add_scaled(
                    data,
                    {VERTEX_COUNT + vertex: sp.Integer(1)},
                    value,
                )
        flex_rows.append(flex)
        data_rows.append(data)
    return flex_rows, data_rows


def exact_eliminate(flex_rows, data_rows, bad_rows, flex_count):
    """Echelonize over Q, with pivots selected only from the flex block."""
    pivots = {}
    candidate_obstructions = []
    bad_obstructions = []
    maximum_widths = {"flex": 0, "candidate": 0, "corrupted": 0}
    pivot_sources = []

    for source_index, (source_flex, source_data, source_bad) in enumerate(
        zip(flex_rows, data_rows, bad_rows)
    ):
        flex = dict(source_flex)
        data = dict(source_data)
        bad = dict(source_bad)
        while flex:
            pivot = min(flex)
            if pivot not in pivots:
                inverse = sp.Integer(1) / flex[pivot]
                flex = scaled(flex, inverse)
                data = scaled(data, inverse)
                bad = scaled(bad, inverse)
                pivots[pivot] = (flex, data, bad, source_index)
                pivot_sources.append(source_index)
                maximum_widths["flex"] = max(
                    maximum_widths["flex"], len(flex)
                )
                maximum_widths["candidate"] = max(
                    maximum_widths["candidate"], len(data)
                )
                maximum_widths["corrupted"] = max(
                    maximum_widths["corrupted"], len(bad)
                )
                break
            factor = -flex[pivot]
            pivot_flex, pivot_data, pivot_bad, _ = pivots[pivot]
            add_scaled(flex, pivot_flex, factor)
            add_scaled(data, pivot_data, factor)
            add_scaled(bad, pivot_bad, factor)
        else:
            if data:
                candidate_obstructions.append((source_index, data))
            if bad:
                bad_obstructions.append((source_index, bad))

    full_rank = len(pivots) == flex_count and set(pivots) == set(range(flex_count))
    return {
        "pivots": pivots,
        "pivot_sources": pivot_sources,
        "full_rank": full_rank,
        "candidate_obstructions": candidate_obstructions,
        "bad_obstructions": bad_obstructions,
        "maximum_widths": maximum_widths,
    }


def back_substitute(pivots, flex_count):
    lift = [None] * flex_count
    for pivot in range(flex_count - 1, -1, -1):
        flex, data, _, _ = pivots[pivot]
        expression = scaled(data, -sp.Integer(1))
        for column, coefficient in flex.items():
            if column == pivot:
                continue
            if column <= pivot or lift[column] is None:
                raise RuntimeError("non-echelon pivot encountered")
            add_scaled(expression, lift[column], -coefficient)
        lift[pivot] = expression
    return lift


def direct_residual(flex_rows, data_rows, lift):
    nonzero_rows = 0
    first = None
    maximum_width = 0
    for row_index, (flex, data) in enumerate(zip(flex_rows, data_rows)):
        residual = dict(data)
        for column, coefficient in flex.items():
            add_scaled(residual, lift[column], coefficient)
        if residual:
            nonzero_rows += 1
            maximum_width = max(maximum_width, len(residual))
            if first is None:
                first = (row_index, residual)
    return nonzero_rows, first, maximum_width


def exact_dict_record(item):
    if item is None:
        return None
    row_index, row = item
    return {
        "row": row_index,
        "coefficients": {
            str(column): str(value) for column, value in sorted(row.items())
        },
    }


def histogram(values):
    return {
        str(key): value
        for key, value in sorted(Counter(values).items())
    }


def lift_statistics(lift, data_count):
    hasher = sha256()
    nonzero = 0
    flex_supports = []
    data_cells = [set() for _ in range(data_count)]
    data_flex_counts = [0] * data_count
    max_num_bits = 0
    max_den_bits = 0
    for flex_index, row in enumerate(lift):
        flex_supports.append(len(row))
        for data_index, raw_value in sorted(row.items()):
            value = sp.Rational(raw_value)
            numerator = int(value.p)
            denominator = int(value.q)
            hasher.update(
                f"{flex_index},{data_index},{numerator},{denominator}\n".encode()
            )
            nonzero += 1
            data_flex_counts[data_index] += 1
            data_cells[data_index].add(flex_index // 6)
            max_num_bits = max(max_num_bits, abs(numerator).bit_length())
            max_den_bits = max(max_den_bits, denominator.bit_length())
    cell_supports = [len(value) for value in data_cells]
    return {
        "sha256": hasher.hexdigest(),
        "nonzero_coefficients": nonzero,
        "flex_row_support_histogram": histogram(flex_supports),
        "data_flex_support_histogram": histogram(data_flex_counts),
        "data_cell_support_histogram": histogram(cell_supports),
        "vertex_scale_cell_support_range": [
            min(cell_supports[:VERTEX_COUNT]),
            max(cell_supports[:VERTEX_COUNT]),
        ],
        "strut_cell_support_range": [
            min(cell_supports[VERTEX_COUNT:]),
            max(cell_supports[VERTEX_COUNT:]),
        ],
        "maximum_numerator_bits": max_num_bits,
        "maximum_denominator_bits": max_den_bits,
    }


def construction_ok(built):
    return bool(
        built["local_geometry"]["controls"]
        and built["transition_control"]
        and built["inverse_control"]
        and built["face_control"]
        and built["local_fixed_ranks"] == {5: 1200}
    )


def solve_construction(name, scale, lapse, built, edges):
    print(f"exact rational lift {name} ({scale},{lapse})", flush=True)
    flex_rows, data_rows = split_exact_rows(
        built, scale, edges, corrupt=False
    )
    bad_flex_rows, bad_rows = split_exact_rows(
        built, scale, edges, corrupt=True
    )
    if bad_flex_rows != flex_rows:
        raise RuntimeError("corruption changed the flex block")
    elimination = exact_eliminate(
        flex_rows, data_rows, bad_rows, FLEX_COUNT
    )
    candidate_consistent = not elimination["candidate_obstructions"]
    corrupted_rejected = bool(elimination["bad_obstructions"])
    lift = None
    residual_count = None
    residual_first = None
    residual_width = None
    statistics = None
    if elimination["full_rank"] and candidate_consistent:
        lift = back_substitute(elimination["pivots"], FLEX_COUNT)
        residual_count, residual_first, residual_width = direct_residual(
            flex_rows, data_rows, lift
        )
        statistics = lift_statistics(lift, DATA_COUNT)
    return {
        "name": name,
        "scale": scale,
        "lapse": lapse,
        "construction_controls": construction_ok(built),
        "equation_rows": len(flex_rows),
        "exact_flex_pivots": len(elimination["pivots"]),
        "pivot_source_rows": len(elimination["pivot_sources"]),
        "maximum_pivot_widths": elimination["maximum_widths"],
        "candidate_consistent": candidate_consistent,
        "candidate_obstruction_count": len(
            elimination["candidate_obstructions"]
        ),
        "candidate_first_obstruction": exact_dict_record(
            elimination["candidate_obstructions"][0]
            if elimination["candidate_obstructions"] else None
        ),
        "direct_residual_nonzero_rows": residual_count,
        "direct_residual_maximum_width": residual_width,
        "direct_first_residual": exact_dict_record(residual_first),
        "corrupted_rejected": corrupted_rejected,
        "corrupted_obstruction_count": len(
            elimination["bad_obstructions"]
        ),
        "corrupted_first_obstruction": exact_dict_record(
            elimination["bad_obstructions"][0]
            if elimination["bad_obstructions"] else None
        ),
        "lift_statistics": statistics,
    }


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "identification_source": IDENTIFICATION_SOURCE,
    "identification_json": IDENTIFICATION_JSON,
    "admissibility_source": ADMISSIBILITY_SOURCE,
    "admissibility_json": ADMISSIBILITY_JSON,
    "refuted_source": REFUTED_SOURCE,
    "refuted_json": REFUTED_JSON,
}
hashes_before = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes_before == EXPECTED_HASHES
check("all rational-lift inputs have frozen provenance", provenance_ok, str(hashes_before))

identification = json.loads(IDENTIFICATION_JSON.read_text())
admissibility_before = json.loads(ADMISSIBILITY_JSON.read_text())
refuted = json.loads(REFUTED_JSON.read_text())
upstream_artifacts_ok = bool(
    identification["outcome"]
    == "CANONICAL_DATA_MODULAR_VERTEX_SCALE_PLUS_ARBITRARY_STRUTS"
    and identification["passed"] == identification["tests"] == 11
    and admissibility_before["outcome"]
    == "CANONICAL_DATA_INTERMEDIATE_MODULAR_OPEN"
    and admissibility_before["passed"] == admissibility_before["tests"] == 10
    and refuted["outcome"] == "CANONICAL_DATA_VERTEX_CARRIER_REFUTED"
    and refuted["passed"] == refuted["tests"] == 10
)
check("the modular positive and old local-lift negative are both retained", upstream_artifacts_ok)

# A solver-only construction known exactly in advance.  The third candidate
# equation is dependent; changing its datum from 3 to 4 must expose one
# zero-flex obstruction.
synthetic_flex = [
    {0: sp.Integer(1)},
    {1: sp.Integer(1)},
    {0: sp.Integer(1), 1: sp.Integer(1)},
]
synthetic_data = [
    {0: sp.Integer(1)},
    {0: sp.Integer(2)},
    {0: sp.Integer(3)},
]
synthetic_bad = [
    {0: sp.Integer(1)},
    {0: sp.Integer(2)},
    {0: sp.Integer(4)},
]
synthetic = exact_eliminate(
    synthetic_flex, synthetic_data, synthetic_bad, 2
)
synthetic_lift = back_substitute(synthetic["pivots"], 2)
synthetic_residual = direct_residual(
    synthetic_flex, synthetic_data, synthetic_lift
)
synthetic_ok = bool(
    synthetic["full_rank"]
    and not synthetic["candidate_obstructions"]
    and synthetic["bad_obstructions"]
    and synthetic_lift == [{0: -1}, {0: -2}]
    and synthetic_residual[0] == 0
)
check("the independent exact solver accepts a known lift and rejects its corruption", synthetic_ok)

print("reconstructing exact complete compatibility equations", flush=True)
namespace = runpy.run_path(str(ADMISSIBILITY_SOURCE))
hashes_after = {name: digest(path) for name, path in paths.items()}
admissibility_after = json.loads(ADMISSIBILITY_JSON.read_text())
upstream_reproduction_ok = bool(
    hashes_after == EXPECTED_HASHES
    and admissibility_after == admissibility_before
)
check("the exact upstream equation source reproduces byte-for-byte", upstream_reproduction_ok)

complex_data = namespace["complex_data"]
edges = tuple(sorted(complex_data["edges"]))
tetrahedra = tuple(complex_data["tetrahedra"])
incidence = incidence_rows(edges, corrupt=False)
incidence_bad = incidence_rows(edges, corrupt=True)
joint = []
for unsigned, corrupted in zip(incidence, incidence_bad):
    row = dict(unsigned)
    for column, value in corrupted.items():
        row[VERTEX_COUNT + column] = value
    joint.append(row)
incidence_ranks = {
    "unsigned": sparse_matrix(incidence, VERTEX_COUNT).rank(),
    "corrupted": sparse_matrix(incidence_bad, VERTEX_COUNT).rank(),
    "joint": sparse_matrix(joint, 2 * VERTEX_COUNT).rank(),
}
combinatorics_ok = bool(
    len(complex_data["vertices"]) == VERTEX_COUNT
    and len(edges) == EDGE_COUNT
    and len(complex_data["triangles"]) == 1200
    and len(tetrahedra) == 600
    and graph_controls(VERTEX_COUNT, edges, tetrahedra)
    and incidence_ranks
    == {"unsigned": 120, "corrupted": 120, "joint": 121}
)
check("the exact graph and 120,120,121 incidence controls pass", combinatorics_ok, str(incidence_ranks))

constructions = []
for scale, lapse in namespace["REPRESENTATIVES"]:
    constructions.append((
        "baseline",
        scale,
        lapse,
        namespace["baseline_builds"][(scale, lapse)],
    ))
for scale, lapse in namespace["REPRESENTATIVES"]:
    constructions.append((
        "alternate_right_inverse",
        scale,
        lapse,
        namespace["build_global"](
            complex_data,
            scale,
            lapse,
            namespace["ETA"],
            namespace["CANONICAL_BASE"],
            use_alternate=True,
        ),
    ))
scale, lapse = namespace["REPRESENTATIVES"][0]
constructions.append((
    "reverse_faces",
    scale,
    lapse,
    namespace["build_global"](
        complex_data,
        scale,
        lapse,
        namespace["ETA"],
        namespace["CANONICAL_BASE"],
        reverse_orientation=True,
    ),
))
odd_canonical = (
    namespace["CANONICAL_BASE"][1],
    namespace["CANONICAL_BASE"][0],
    namespace["CANONICAL_BASE"][2],
    namespace["CANONICAL_BASE"][3],
)
constructions.append((
    "odd_relabelling",
    scale,
    lapse,
    namespace["build_global"](
        complex_data,
        scale,
        lapse,
        namespace["ETA"],
        odd_canonical,
    ),
))
constructions.append((
    "metric_sign",
    scale,
    lapse,
    namespace["build_global"](
        complex_data,
        scale,
        lapse,
        -namespace["ETA"],
        namespace["CANONICAL_BASE"],
    ),
))

records = [
    solve_construction(name, scale, lapse, built, edges)
    for name, scale, lapse, built in constructions
]

construction_controls_ok = all(
    record["construction_controls"] and record["equation_rows"] == 6000
    for record in records
)
check("all seven exact constructions retain every geometry and face control", construction_controls_ok)

full_pivots_ok = all(
    record["exact_flex_pivots"] == FLEX_COUNT
    and record["pivot_source_rows"] == FLEX_COUNT
    for record in records
)
check("every rational flex block has exactly 3600 deterministic pivots", full_pivots_ok)

candidate_decisions = {record["candidate_consistent"] for record in records}
decision_agreement_ok = len(candidate_decisions) == 1
check("all exact constructions agree on candidate consistency", decision_agreement_ok, str(candidate_decisions))

direct_residuals_ok = all(
    (not record["candidate_consistent"])
    or (
        record["direct_residual_nonzero_rows"] == 0
        and record["direct_first_residual"] is None
        and record["lift_statistics"] is not None
    )
    for record in records
)
check("every accepted lift vanishes on every original rational equation", direct_residuals_ok)

corruption_ok = all(
    record["corrupted_rejected"]
    and record["corrupted_obstruction_count"] > 0
    and record["corrupted_first_obstruction"] is not None
    for record in records
)
check("the distinct one-row-corrupted image is rejected exactly throughout", corruption_ok)

serialization_ok = all(
    (not record["candidate_consistent"])
    or (
        len(record["lift_statistics"]["sha256"]) == 64
        and record["lift_statistics"]["nonzero_coefficients"] > 0
    )
    for record in records
)
check("every exact lift has a complete canonical coefficient digest", serialization_ok)

controls_ok = bool(
    provenance_ok
    and upstream_artifacts_ok
    and synthetic_ok
    and upstream_reproduction_ok
    and combinatorics_ok
    and construction_controls_ok
    and full_pivots_ok
    and direct_residuals_ok
    and corruption_ok
    and serialization_ok
)

candidate_consistent = bool(
    decision_agreement_ok and next(iter(candidate_decisions))
)
if not controls_ok:
    outcome = "RATIONAL_DATA_LIFT_CONTROL_FAILED"
elif not decision_agreement_ok:
    outcome = "RATIONAL_DATA_LIFT_DISAGREEMENT_OPEN"
elif not candidate_consistent:
    outcome = "RATIONAL_VERTEX_STRUT_DATA_LIFT_REFUTED"
else:
    outcome = "RATIONAL_VERTEX_STRUT_DATA_LIFT_DERIVED"

allowed = {
    "RATIONAL_DATA_LIFT_CONTROL_FAILED",
    "RATIONAL_DATA_LIFT_DISAGREEMENT_OPEN",
    "RATIONAL_VERTEX_STRUT_DATA_LIFT_REFUTED",
    "RATIONAL_VERTEX_STRUT_DATA_LIFT_DERIVED",
}
check("the preregistered rational-lift hierarchy assigns one outcome", outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes_before,
    "source_sha256": digest(Path(__file__)),
    "first_sorted_edge": list(edges[0]),
    "incidence_ranks_over_Q": incidence_ranks,
    "records": records,
    "classification": {
        "rational_boundary_data_lift": (
            "DERIVED"
            if outcome == "RATIONAL_VERTEX_STRUT_DATA_LIFT_DERIVED"
            else "OPEN"
        ),
        "old_local_lift": "DERIVED NEGATIVE",
        "action_or_symplectic_dynamics": "NOT EVALUATED",
        "physical_tick_speed_G_or_Planck_scale": "NOT EVALUATED",
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
    stats = record["lift_statistics"] or {}
    print(
        f"{record['name']} ({record['scale']},{record['lapse']}): "
        f"pivots={record['exact_flex_pivots']}, "
        f"candidate={record['candidate_consistent']}, "
        f"residual_rows={record['direct_residual_nonzero_rows']}, "
        f"corrupted={record['corrupted_rejected']}, "
        f"lift_nnz={stats.get('nonzero_coefficients')}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

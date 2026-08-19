#!/usr/bin/env python3
"""Test the disclosed unsigned-incidence plus arbitrary-strut data carrier."""

from collections import deque
from hashlib import sha256
import json
from pathlib import Path
import runpy

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_canonical_data_space_identification.json"
PRIOR_ART = (
    ROOT
    / "docs/gravity/gravity_600cell_canonical_data_space_identification_prior_art.md"
)
PROTOCOL = (
    ROOT
    / "docs/gravity/gravity_600cell_canonical_data_space_identification_protocol.md"
)
DEVIATION = (
    ROOT
    / "docs/gravity/gravity_600cell_canonical_data_projection_protocol_deviation.md"
)
ADMISSIBILITY_PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_canonical_data_admissibility_protocol.md"
)
ADMISSIBILITY_SOURCE = HERE / "verify_gravity_600cell_canonical_data_admissibility.py"
ADMISSIBILITY_JSON = HERE / "gravity_600cell_canonical_data_admissibility.json"
PROJECTION_SOURCE = HERE / "verify_gravity_600cell_canonical_data_projection.py"
PROJECTION_JSON = HERE / "gravity_600cell_canonical_data_projection.json"

PROTOCOL_COMMIT = "8576c84"
EXPECTED_HASHES = {
    "prior_art": "f648178c79e5b23fa586d34614bc7e030d199880c5225ee3a735809a9f563f49",
    "protocol": "b5b011abcfd2f98a228505360c7fe856906ac1a203263e6bbaedbf8f0abc3209",
    "projection_deviation": "1ffd2c7d29537cada652167c5d73ff0093eecbefd59eff84eed095c383a9f809",
    "admissibility_protocol": "8db29cb9af699da660b969988eeb76c5e605e67c5ec65716795ada2e34674185",
    "admissibility_source": "4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e",
    "admissibility_json": "fa45c80739ca0dda4f82c9da98a4b22f4d8a18c182a40696a2a22d1d26ec89a1",
    "projection_source": "6c59003278d91058fc89efbbc9ffa6d105785001909093555b3ee13bf4be5a1c",
    "projection_json": "f011ef9848a6139408a9f8495a12e0d8e0050e04f39aa5c00ca88c02dde26beb",
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


def add_value(row, column, value):
    value = sp.cancel(value)
    if not value:
        return
    updated = sp.cancel(row.get(column, 0) + value)
    if updated:
        row[column] = updated
    elif column in row:
        del row[column]


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


def incidence_rows(vertex_count, edges, corrupt=False):
    rows = []
    for edge_index, (left, right) in enumerate(edges):
        row = {left: 1, right: 1}
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


def compose_rows(built, scale, edges, fixed_count, edge_count, corrupt=False):
    result = []
    for original in built["rows"]:
        row = {
            column: value
            for column, value in original.items()
            if column < fixed_count
        }
        for column, value in original.items():
            if fixed_count <= column < fixed_count + edge_count:
                edge_index = column - fixed_count
                left, right = edges[edge_index]
                add_value(row, fixed_count + left, value * 8 * scale)
                if not (corrupt and edge_index == 0):
                    add_value(row, fixed_count + right, value * 8 * scale)
        result.append(row)
    return result


def restricted_rows(built, fixed_count, edge_count, include_edges, include_struts):
    result = []
    strut_start = fixed_count + edge_count
    for original in built["rows"]:
        row = {}
        for column, value in original.items():
            if column < fixed_count:
                row[column] = value
            elif column < strut_start and include_edges:
                row[column] = value
            elif column >= strut_start and include_struts:
                target = (
                    column
                    if include_edges
                    else fixed_count + column - strut_start
                )
                row[target] = value
        result.append(row)
    return result


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "projection_deviation": DEVIATION,
    "admissibility_protocol": ADMISSIBILITY_PROTOCOL,
    "admissibility_source": ADMISSIBILITY_SOURCE,
    "admissibility_json": ADMISSIBILITY_JSON,
    "projection_source": PROJECTION_SOURCE,
    "projection_json": PROJECTION_JSON,
}
hashes_before = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes_before == EXPECTED_HASHES
check("all data-space identification inputs have frozen provenance", provenance_ok, str(hashes_before))

projection = json.loads(PROJECTION_JSON.read_text())
expected_dimensions = {
    "kernel": 240,
    "edge_only": 120,
    "strut_only": 120,
    "edge_projection": 120,
    "strut_projection": 120,
}
projection_ok = bool(
    projection["outcome"] == "CANONICAL_DATA_PROJECTION_STABLE_MODULAR_OPEN"
    and projection["passed"] == projection["tests"] == 11
    and projection["target_comparison"] == "NOT PERFORMED"
    and len(projection["records"]) == 7
    and all(
        prime_record["dimensions"] == expected_dimensions
        for record in projection["records"]
        for prime_record in record["prime_records"].values()
    )
)
check("the frozen 240=120+120 projection artifact is reproduced", projection_ok)

print("reconstructing complete compatibility matrices", flush=True)
namespace = runpy.run_path(str(ADMISSIBILITY_SOURCE))
hashes_after = {name: digest(path) for name, path in paths.items()}
admissibility = json.loads(ADMISSIBILITY_JSON.read_text())
upstream_ok = bool(
    hashes_after == EXPECTED_HASHES
    and admissibility["outcome"] == "CANONICAL_DATA_INTERMEDIATE_MODULAR_OPEN"
    and admissibility["passed"] == admissibility["tests"] == 10
)
check("the upstream compatibility source reproduces byte-for-byte", upstream_ok)

complex_data = namespace["complex_data"]
edges = tuple(sorted(complex_data["edges"]))
tetrahedra = tuple(complex_data["tetrahedra"])
vertex_count = namespace["VERTEX_COUNT"]
edge_count = namespace["EDGE_COUNT"]
fixed_count = namespace["CELL_COLUMNS"]
primes = tuple(namespace["PRIMES"])
modular_rank = namespace["modular_rank"]

combinatorics_ok = bool(
    vertex_count == 120
    and edge_count == 720
    and graph_controls(vertex_count, edges, tetrahedra)
)
check("the 600-cell graph is connected, 12-regular and contains a triangle", combinatorics_ok)

unsigned_rows = incidence_rows(vertex_count, edges)
corrupt_rows = incidence_rows(vertex_count, edges, corrupt=True)
joint_rows = []
for unsigned, corrupt in zip(unsigned_rows, corrupt_rows):
    joint = dict(unsigned)
    for column, value in corrupt.items():
        joint[vertex_count + column] = value
    joint_rows.append(joint)

unsigned_matrix = sparse_matrix(unsigned_rows, vertex_count)
corrupt_matrix = sparse_matrix(corrupt_rows, vertex_count)
joint_matrix = sparse_matrix(joint_rows, 2 * vertex_count)
rational_incidence_ranks = {
    "unsigned": unsigned_matrix.rank(),
    "corrupted": corrupt_matrix.rank(),
    "joint": joint_matrix.rank(),
}
rational_incidence_ok = rational_incidence_ranks == {
    "unsigned": 120,
    "corrupted": 120,
    "joint": 121,
}
check("the exact rational incidence and corruption ranks are 120,120,121", rational_incidence_ok, str(rational_incidence_ranks))

modular_incidence_ranks = {}
for prime in primes:
    modular_incidence_ranks[str(prime)] = {
        "unsigned": modular_rank(unsigned_rows, vertex_count, prime)[0],
        "corrupted": modular_rank(corrupt_rows, vertex_count, prime)[0],
        "joint": modular_rank(joint_rows, 2 * vertex_count, prime)[0],
    }
modular_incidence_ok = all(
    ranks == {"unsigned": 120, "corrupted": 120, "joint": 121}
    for ranks in modular_incidence_ranks.values()
)
check("both primes reproduce the 120,120,121 incidence ranks", modular_incidence_ok, str(modular_incidence_ranks))


def construction_ok(built):
    return bool(
        built["local_geometry"]["controls"]
        and built["transition_control"]
        and built["inverse_control"]
        and built["face_control"]
    )


def identify_record(name, scale, lapse, built):
    print(f"testing data-space candidate {name} ({scale},{lapse})", flush=True)
    fixed_rows = restricted_rows(
        built, fixed_count, edge_count, include_edges=False, include_struts=False
    )
    fixed_edge_rows = restricted_rows(
        built, fixed_count, edge_count, include_edges=True, include_struts=False
    )
    fixed_strut_rows = restricted_rows(
        built, fixed_count, edge_count, include_edges=False, include_struts=True
    )
    candidate_rows = compose_rows(
        built, scale, edges, fixed_count, edge_count, corrupt=False
    )
    corrupted_candidate_rows = compose_rows(
        built, scale, edges, fixed_count, edge_count, corrupt=True
    )
    prime_records = {}
    for prime in primes:
        fixed_rank, fixed_width = modular_rank(fixed_rows, fixed_count, prime)
        edge_rank, edge_width = modular_rank(
            fixed_edge_rows, fixed_count + edge_count, prime
        )
        strut_rank, strut_width = modular_rank(
            fixed_strut_rows, fixed_count + vertex_count, prime
        )
        candidate_rank, candidate_width = modular_rank(
            candidate_rows, fixed_count + vertex_count, prime
        )
        corrupted_rank, corrupted_width = modular_rank(
            corrupted_candidate_rows, fixed_count + vertex_count, prime
        )
        prime_records[str(prime)] = {
            "ranks": {
                "fixed": fixed_rank,
                "fixed_edge": edge_rank,
                "fixed_strut": strut_rank,
                "candidate": candidate_rank,
                "corrupted_candidate": corrupted_rank,
            },
            "candidate_included": candidate_rank == fixed_rank,
            "strut_ambient_included": strut_rank == fixed_rank,
            "corrupted_candidate_included": corrupted_rank == fixed_rank,
            "maximum_elimination_widths": {
                "fixed": fixed_width,
                "fixed_edge": edge_width,
                "fixed_strut": strut_width,
                "candidate": candidate_width,
                "corrupted_candidate": corrupted_width,
            },
        }
    return {
        "name": name,
        "scale": scale,
        "lapse": lapse,
        "construction_controls": construction_ok(built),
        "prime_records": prime_records,
    }


records = []
for scale, lapse in namespace["REPRESENTATIVES"]:
    records.append(
        identify_record(
            "baseline",
            scale,
            lapse,
            namespace["baseline_builds"][(scale, lapse)],
        )
    )

for scale, lapse in namespace["REPRESENTATIVES"]:
    alternate = namespace["build_global"](
        complex_data,
        scale,
        lapse,
        namespace["ETA"],
        namespace["CANONICAL_BASE"],
        use_alternate=True,
    )
    records.append(
        identify_record("alternate_right_inverse", scale, lapse, alternate)
    )

scale, lapse = namespace["REPRESENTATIVES"][0]
reverse = namespace["build_global"](
    complex_data,
    scale,
    lapse,
    namespace["ETA"],
    namespace["CANONICAL_BASE"],
    reverse_orientation=True,
)
records.append(identify_record("reverse_faces", scale, lapse, reverse))

odd_canonical = (
    namespace["CANONICAL_BASE"][1],
    namespace["CANONICAL_BASE"][0],
    namespace["CANONICAL_BASE"][2],
    namespace["CANONICAL_BASE"][3],
)
odd = namespace["build_global"](
    complex_data,
    scale,
    lapse,
    namespace["ETA"],
    odd_canonical,
)
records.append(identify_record("odd_relabelling", scale, lapse, odd))

metric_sign = namespace["build_global"](
    complex_data,
    scale,
    lapse,
    -namespace["ETA"],
    namespace["CANONICAL_BASE"],
)
records.append(identify_record("metric_sign", scale, lapse, metric_sign))

construction_controls_ok = all(record["construction_controls"] for record in records)
check("all seven constructions retain every local and face control", construction_controls_ok)

frozen_rank_controls_ok = all(
    prime_record["ranks"]["fixed"] == 3600
    and prime_record["ranks"]["fixed_edge"] == 4200
    and prime_record["ranks"]["fixed_strut"] == 3600
    for record in records
    for prime_record in record["prime_records"].values()
)
check("all fixed, edge and strut ranks reproduce the frozen census", frozen_rank_controls_ok)

decision_tuples = {
    (
        prime_record["candidate_included"],
        prime_record["strut_ambient_included"],
        prime_record["corrupted_candidate_included"],
    )
    for record in records
    for prime_record in record["prime_records"].values()
}
decision_agreement_ok = len(decision_tuples) == 1
check("both primes and all constructions give one inclusion decision", decision_agreement_ok, str(decision_tuples))

candidate_included = bool(
    decision_agreement_ok and next(iter(decision_tuples))[0]
)
strut_included = bool(
    decision_agreement_ok and next(iter(decision_tuples))[1]
)
corrupted_included = bool(
    decision_agreement_ok and next(iter(decision_tuples))[2]
)
conditional_corruption_ok = bool(
    not candidate_included or not corrupted_included
)
check("an accepted candidate excludes the distinct corrupted image", conditional_corruption_ok)

controls_ok = bool(
    provenance_ok
    and projection_ok
    and upstream_ok
    and combinatorics_ok
    and rational_incidence_ok
    and modular_incidence_ok
    and construction_controls_ok
    and frozen_rank_controls_ok
    and conditional_corruption_ok
)

if not controls_ok:
    outcome = "CANONICAL_DATA_SPACE_CONTROL_FAILED"
elif not decision_agreement_ok:
    outcome = "CANONICAL_DATA_SPACE_DISAGREEMENT_OPEN"
elif not candidate_included:
    outcome = "CANONICAL_DATA_VERTEX_SCALE_IMAGE_REFUTED"
elif candidate_included and strut_included and not corrupted_included:
    outcome = "CANONICAL_DATA_MODULAR_VERTEX_SCALE_PLUS_ARBITRARY_STRUTS"
else:
    outcome = "CANONICAL_DATA_SPACE_CONTROL_FAILED"

allowed = {
    "CANONICAL_DATA_SPACE_CONTROL_FAILED",
    "CANONICAL_DATA_SPACE_DISAGREEMENT_OPEN",
    "CANONICAL_DATA_VERTEX_SCALE_IMAGE_REFUTED",
    "CANONICAL_DATA_MODULAR_VERTEX_SCALE_PLUS_ARBITRARY_STRUTS",
}
check("the preregistered data-space hierarchy assigns one outcome", outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes_before,
    "source_sha256": digest(Path(__file__)),
    "first_sorted_edge": list(edges[0]),
    "rational_incidence_ranks": rational_incidence_ranks,
    "modular_incidence_ranks": modular_incidence_ranks,
    "records": records,
    "classification": {
        "modular_data_space": (
            "unsigned vertex-edge scale image plus arbitrary struts"
            if outcome
            == "CANONICAL_DATA_MODULAR_VERTEX_SCALE_PLUS_ARBITRARY_STRUTS"
            else "OPEN"
        ),
        "rational_equality": "OPEN",
        "global_cell_flex_lift": "OPEN",
        "action_or_physics": "NOT EVALUATED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
for record in records:
    decisions = {
        prime: (
            value["candidate_included"],
            value["strut_ambient_included"],
            value["corrupted_candidate_included"],
        )
        for prime, value in record["prime_records"].items()
    }
    print(
        f"{record['name']} ({record['scale']},{record['lapse']}): {decisions}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

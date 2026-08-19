#!/usr/bin/env python3
"""Classify the frozen 240-dimensional canonical-data kernel exactly."""

from hashlib import sha256
import json
from pathlib import Path
import runpy

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_canonical_data_carrier.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_canonical_data_carrier_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_canonical_data_carrier_protocol.md"
CLASSIFIER_CORRECTION = (
    ROOT
    / "docs/gravity/gravity_600cell_canonical_data_carrier_classifier_correction.md"
)
ADMISSIBILITY_PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_canonical_data_admissibility_protocol.md"
)
ADMISSIBILITY_SOURCE = HERE / "verify_gravity_600cell_canonical_data_admissibility.py"
ADMISSIBILITY_JSON = HERE / "gravity_600cell_canonical_data_admissibility.json"

PROTOCOL_COMMIT = "fbabe33"
EXPECTED_HASHES = {
    "prior_art": "fd7158f80af48fadc88c121c6001e258d2bccab480d8624b709f0ee142d145af",
    "protocol": "e5b1217cb49da317e6c334ed8b5458b6139b31087f740b31e6498dd7e74f4bb3",
    "classifier_correction": "f448aa5858cb4e2d15d17f50a61da933bdb3e54796157cc08aae43558f5feded",
    "admissibility_protocol": "8db29cb9af699da660b969988eeb76c5e605e67c5ec65716795ada2e34674185",
    "admissibility_source": "4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e",
    "admissibility_json": "fa45c80739ca0dda4f82c9da98a4b22f4d8a18c182a40696a2a22d1d26ec89a1",
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


def map_residuals(matrix_rows, coordinate_map):
    nonzero_rows = 0
    maximum_support = 0
    first_nonzero = None
    for row_number, row in enumerate(matrix_rows):
        residual = {}
        for coordinate, coefficient in row.items():
            for carrier_column, value in coordinate_map[coordinate].items():
                add_value(residual, carrier_column, coefficient * value)
        maximum_support = max(maximum_support, len(residual))
        if residual:
            nonzero_rows += 1
            if first_nonzero is None:
                first_nonzero = {
                    "row": row_number,
                    "residual": {str(key): str(value) for key, value in residual.items()},
                }
    return {
        "nonzero_rows": nonzero_rows,
        "maximum_residual_support": maximum_support,
        "first_nonzero": first_nonzero,
    }


def column_sum(coordinate_map, carrier_columns):
    result = {}
    selected = set(carrier_columns)
    for coordinate, row in enumerate(coordinate_map):
        value = sp.cancel(sum(entry for key, entry in row.items() if key in selected))
        if value:
            result[coordinate] = value
    return result


def build_carrier(namespace, built, scale, lapse, right_inverse):
    canonical = namespace["CANONICAL_BASE"]
    normal = namespace["NORMAL"]
    cell_columns = namespace["CELL_COLUMNS"]
    edge_count = namespace["EDGE_COUNT"]
    vertex_count = namespace["VERTEX_COUNT"]
    total_columns = namespace["TOTAL_COLUMNS"]
    complex_data = namespace["complex_data"]
    tetrahedra = complex_data["tetrahedra"]
    edges = tuple(sorted(complex_data["edges"]))
    local = built["local_geometry"]
    jacobian = local["jacobian"]
    local_kernel = local["kernel"]

    displacement = sp.zeros(16, 8)
    for local_vertex in range(4):
        displacement[
            4 * local_vertex : 4 * local_vertex + 4, local_vertex
        ] = canonical[local_vertex]
        displacement[
            4 * local_vertex : 4 * local_vertex + 4, 4 + local_vertex
        ] = normal
    local_data = jacobian * displacement
    residual = displacement - right_inverse * local_data
    local_flex = namespace["solve_full_column"](local_kernel, residual)
    decomposition_ok = bool(
        displacement == right_inverse * local_data + local_kernel * local_flex
    )

    formula_ok = True
    for row, (left, right) in enumerate(local["edge_pairs"]):
        expected = sp.zeros(1, 8)
        expected[0, left] = 8 * scale
        expected[0, right] = 8 * scale
        formula_ok &= bool(local_data[row, :] == expected)
    for local_vertex in range(4):
        expected = sp.zeros(1, 8)
        expected[0, local_vertex] = 6 * (scale - 1)
        expected[0, 4 + local_vertex] = -2 * lapse
        formula_ok &= bool(local_data[6 + local_vertex, :] == expected)

    coordinate_map = [dict() for _ in range(total_columns)]
    for tetrahedron_index, tetrahedron in enumerate(tetrahedra):
        for flex_row in range(6):
            target = coordinate_map[6 * tetrahedron_index + flex_row]
            for local_vertex, global_vertex in enumerate(tetrahedron):
                add_value(target, global_vertex, local_flex[flex_row, local_vertex])
                add_value(
                    target,
                    vertex_count + global_vertex,
                    local_flex[flex_row, 4 + local_vertex],
                )
    for edge_index, (left, right) in enumerate(edges):
        target = coordinate_map[cell_columns + edge_index]
        add_value(target, left, 8 * scale)
        add_value(target, right, 8 * scale)
    for vertex in range(vertex_count):
        target = coordinate_map[cell_columns + edge_count + vertex]
        add_value(target, vertex, 6 * (scale - 1))
        add_value(target, vertex_count + vertex, -2 * lapse)

    wrong_difference = [dict(row) for row in coordinate_map]
    for edge_index, (left, right) in enumerate(edges):
        target = wrong_difference[cell_columns + edge_index]
        target.clear()
        add_value(target, left, 8 * scale)
        add_value(target, right, -8 * scale)

    wrong_lapse = [dict(row) for row in coordinate_map]
    first_strut = wrong_lapse[cell_columns + edge_count]
    first_strut.pop(vertex_count, None)

    incidence_rows = []
    for left, right in edges:
        incidence_rows.append({left: 1, right: 1})
    incidence_ranks = {
        str(prime): namespace["modular_rank"](
            incidence_rows, vertex_count, prime
        )[0]
        for prime in namespace["PRIMES"]
    }
    data_rank = (
        2 * vertex_count
        if all(value == vertex_count for value in incidence_ranks.values())
        and scale != 0
        and lapse != 0
        else None
    )
    return {
        "coordinate_map": coordinate_map,
        "wrong_difference": wrong_difference,
        "wrong_lapse": wrong_lapse,
        "local_flex": local_flex,
        "local_data": local_data,
        "formula_ok": formula_ok,
        "decomposition_ok": decomposition_ok,
        "incidence_ranks": incidence_ranks,
        "data_rank": data_rank,
    }


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "classifier_correction": CLASSIFIER_CORRECTION,
    "admissibility_protocol": ADMISSIBILITY_PROTOCOL,
    "admissibility_source": ADMISSIBILITY_SOURCE,
    "admissibility_json": ADMISSIBILITY_JSON,
}
hashes_before = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes_before == EXPECTED_HASHES
check("all carrier-classification inputs have frozen provenance", provenance_ok, str(hashes_before))

print("reconstructing the frozen augmented matrices", flush=True)
namespace = runpy.run_path(str(ADMISSIBILITY_SOURCE))
hashes_after = {name: digest(path) for name, path in paths.items()}
artifact_reproduced = hashes_after == EXPECTED_HASHES
frozen = json.loads(ADMISSIBILITY_JSON.read_text())
frozen_ok = bool(
    artifact_reproduced
    and frozen["outcome"] == "CANONICAL_DATA_INTERMEDIATE_MODULAR_OPEN"
    and frozen["passed"] == frozen["tests"] == 10
    and all(
        set(record["fixed_ranks"].values()) == {3600}
        and set(record["augmented_ranks"].values()) == {4200}
        and set(record["modular_nullities"].values()) == {240}
        for record in frozen["records"]
    )
)
check("the target-blind 4200/4440 census reproduces byte-for-byte", frozen_ok)

records = []
carrier_controls = True
inclusion_ok = True
data_rank_ok = True
alternate_ok = True
difference_control_ok = True
lapse_control_ok = True
constant_control_ok = True

for scale, lapse in namespace["REPRESENTATIVES"]:
    built = namespace["baseline_builds"][(scale, lapse)]
    carrier = build_carrier(
        namespace,
        built,
        scale,
        lapse,
        built["local_geometry"]["right_inverse"],
    )
    inclusion = map_residuals(built["rows"], carrier["coordinate_map"])
    wrong_difference = map_residuals(
        built["rows"], carrier["wrong_difference"]
    )
    wrong_lapse = map_residuals(built["rows"], carrier["wrong_lapse"])

    print(
        f"rebuilding alternate graph at (lambda,tau)=({scale},{lapse})",
        flush=True,
    )
    alternate_built = namespace["build_global"](
        namespace["complex_data"],
        scale,
        lapse,
        namespace["ETA"],
        namespace["CANONICAL_BASE"],
        use_alternate=True,
    )
    alternate = build_carrier(
        namespace,
        alternate_built,
        scale,
        lapse,
        alternate_built["local_geometry"]["alternate_right_inverse"],
    )
    alternate_inclusion = map_residuals(
        alternate_built["rows"], alternate["coordinate_map"]
    )

    homothetic, homothetic_data_rank = namespace["homothetic_vectors"](
        built["local_geometry"],
        built["local_geometry"]["right_inverse"],
        len(namespace["tetrahedra"]),
    )
    homothetic_by_name = dict(homothetic)
    sigma_constant = column_sum(
        carrier["coordinate_map"], range(namespace["VERTEX_COUNT"])
    )
    nu_constant = column_sum(
        carrier["coordinate_map"],
        range(namespace["VERTEX_COUNT"], 2 * namespace["VERTEX_COUNT"]),
    )
    constants_match = bool(
        homothetic_data_rank == 2
        and sigma_constant == homothetic_by_name["scale"]
        and nu_constant == homothetic_by_name["lapse"]
    )

    carrier_controls &= bool(
        carrier["formula_ok"]
        and carrier["decomposition_ok"]
        and alternate["formula_ok"]
        and alternate["decomposition_ok"]
    )
    inclusion_ok &= inclusion["nonzero_rows"] == 0
    data_rank_ok &= carrier["data_rank"] == 240
    alternate_ok &= bool(
        alternate["data_rank"] == 240
        and (
            (inclusion["nonzero_rows"] == 0)
            == (alternate_inclusion["nonzero_rows"] == 0)
        )
    )
    difference_control_ok &= wrong_difference["nonzero_rows"] > 0
    lapse_control_ok &= wrong_lapse["nonzero_rows"] > 0
    constant_control_ok &= constants_match
    records.append({
        "scale": scale,
        "lapse": lapse,
        "local_formula_derived": carrier["formula_ok"],
        "local_decomposition_exact": carrier["decomposition_ok"],
        "unsigned_incidence_modular_ranks": carrier["incidence_ranks"],
        "rational_data_rank": carrier["data_rank"],
        "complete_inclusion": inclusion,
        "alternate_graph_inclusion": alternate_inclusion,
        "wrong_endpoint_difference": wrong_difference,
        "wrong_first_lapse": wrong_lapse,
        "constant_columns_match_frozen_controls": constants_match,
    })

check("the local radial/lapse formulas are derived and decomposed exactly", carrier_controls)
check("all complete rational face rows annihilate the 240-column carrier", inclusion_ok)
check("the data projection has exact rational rank 240", data_rank_ok)
check("the inclusion decision is independent of the right-inverse graph", alternate_ok)
check("the endpoint-difference carrier is rejected exactly", difference_control_ok)
check("deleting one lapse datum is rejected exactly", lapse_control_ok)
check("the two constant columns reproduce the frozen homothetic controls", constant_control_ok)

controls_ok = bool(
    provenance_ok
    and frozen_ok
    and carrier_controls
    and data_rank_ok
    and alternate_ok
    and difference_control_ok
    and lapse_control_ok
    and constant_control_ok
)
exact_exhaustion = bool(controls_ok and inclusion_ok)

if not controls_ok:
    outcome = "CANONICAL_DATA_CARRIER_CONTROL_FAILED"
elif exact_exhaustion:
    outcome = "CANONICAL_DATA_EXACTLY_VERTEX_SCALE_PLUS_LAPSE"
elif not inclusion_ok:
    outcome = "CANONICAL_DATA_VERTEX_CARRIER_REFUTED"
else:
    outcome = "CANONICAL_DATA_CARRIER_OPEN"

allowed = {
    "CANONICAL_DATA_CARRIER_CONTROL_FAILED",
    "CANONICAL_DATA_EXACTLY_VERTEX_SCALE_PLUS_LAPSE",
    "CANONICAL_DATA_VERTEX_CARRIER_PROPER_SUBSPACE",
    "CANONICAL_DATA_VERTEX_CARRIER_REFUTED",
    "CANONICAL_DATA_CARRIER_OPEN",
}
check("the preregistered carrier hierarchy assigns one outcome", outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes_before,
    "frozen_artifact_reproduced": artifact_reproduced,
    "records": records,
    "classification": {
        "rational_kernel_dimension": 240 if exact_exhaustion else "OPEN",
        "kernel_carrier": (
            "Q^120 vertex scale plus Q^120 vertex lapse"
            if exact_exhaustion
            else "OPEN"
        ),
        "independent_upper_edge_shape_directions": (
            0 if exact_exhaustion else "OPEN"
        ),
        "action_or_dynamics": "NOT EVALUATED",
        "tensor_modes_or_speed": "NOT EVALUATED",
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
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"data_rank={record['rational_data_rank']}, "
        f"nonzero_rows={record['complete_inclusion']['nonzero_rows']}, "
        f"wrong_difference={record['wrong_endpoint_difference']['nonzero_rows']}, "
        f"wrong_lapse={record['wrong_first_lapse']['nonzero_rows']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

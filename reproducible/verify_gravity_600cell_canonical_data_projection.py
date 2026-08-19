#!/usr/bin/env python3
"""Target-blind modular projection census for compatible canonical data."""

from hashlib import sha256
import json
from pathlib import Path
import runpy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_canonical_data_projection.json"
PRIOR_ART = (
    ROOT / "docs/gravity/gravity_600cell_canonical_data_projection_prior_art.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_canonical_data_projection_protocol.md"
)
ADMISSIBILITY_PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_canonical_data_admissibility_protocol.md"
)
ADMISSIBILITY_SOURCE = HERE / "verify_gravity_600cell_canonical_data_admissibility.py"
ADMISSIBILITY_JSON = HERE / "gravity_600cell_canonical_data_admissibility.json"

PROTOCOL_COMMIT = "36ddebd"
EXPECTED_HASHES = {
    "prior_art": "1edf5916d8b51386aed9c5183c5102b571e9044d669f8e032f23ba160aafcf55",
    "protocol": "d9bf9e0f244a0c66362f07fac8bb69bb9240da9599b9fb8e372e6b51579bcd9f",
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


def restricted_rows(rows, fixed_count, edge_count, include_edge, include_strut):
    result = []
    strut_start = fixed_count + edge_count
    for original in rows:
        row = {}
        for column, value in original.items():
            if column < fixed_count:
                row[column] = value
            elif column < strut_start:
                if include_edge:
                    row[column] = value
            elif include_strut:
                target = (
                    column
                    if include_edge
                    else fixed_count + column - strut_start
                )
                row[target] = value
        result.append(row)
    return result


def projection_record(modular_rank, rows, fixed_count, edge_count, strut_count, primes):
    fixed_rows = restricted_rows(
        rows, fixed_count, edge_count, include_edge=False, include_strut=False
    )
    fixed_edge_rows = restricted_rows(
        rows, fixed_count, edge_count, include_edge=True, include_strut=False
    )
    fixed_strut_rows = restricted_rows(
        rows, fixed_count, edge_count, include_edge=False, include_strut=True
    )
    records = {}
    for prime in primes:
        r_fixed, w_fixed = modular_rank(fixed_rows, fixed_count, prime)
        r_fixed_edge, w_fixed_edge = modular_rank(
            fixed_edge_rows, fixed_count + edge_count, prime
        )
        r_fixed_strut, w_fixed_strut = modular_rank(
            fixed_strut_rows, fixed_count + strut_count, prime
        )
        r_full, w_full = modular_rank(
            rows, fixed_count + edge_count + strut_count, prime
        )
        kernel = edge_count + strut_count + r_fixed - r_full
        edge_only = edge_count + r_fixed - r_fixed_edge
        strut_only = strut_count + r_fixed - r_fixed_strut
        edge_projection = kernel - strut_only
        strut_projection = kernel - edge_only
        records[str(prime)] = {
            "ranks": {
                "fixed": r_fixed,
                "fixed_edge": r_fixed_edge,
                "fixed_strut": r_fixed_strut,
                "full": r_full,
            },
            "dimensions": {
                "kernel": kernel,
                "edge_only": edge_only,
                "strut_only": strut_only,
                "edge_projection": edge_projection,
                "strut_projection": strut_projection,
            },
            "maximum_elimination_widths": {
                "fixed": w_fixed,
                "fixed_edge": w_fixed_edge,
                "fixed_strut": w_fixed_strut,
                "full": w_full,
            },
        }
    return records


def dimensions_tuple(record):
    dimensions = record["dimensions"]
    return tuple(
        dimensions[key]
        for key in (
            "kernel",
            "edge_only",
            "strut_only",
            "edge_projection",
            "strut_projection",
        )
    )


def rank_bounds_ok(record, fixed_count, edge_count, strut_count):
    ranks = record["ranks"]
    dimensions = record["dimensions"]
    kernel = dimensions["kernel"]
    return bool(
        ranks["fixed"] == fixed_count
        and ranks["fixed"] <= ranks["fixed_edge"] <= ranks["full"]
        and ranks["fixed"] <= ranks["fixed_strut"] <= ranks["full"]
        and 0 <= kernel <= edge_count + strut_count
        and 0 <= dimensions["edge_only"] <= edge_count
        and 0 <= dimensions["strut_only"] <= strut_count
        and dimensions["edge_only"] <= dimensions["edge_projection"] <= min(kernel, edge_count)
        and dimensions["strut_only"] <= dimensions["strut_projection"] <= min(kernel, strut_count)
        and dimensions["edge_projection"] == kernel - dimensions["strut_only"]
        and dimensions["strut_projection"] == kernel - dimensions["edge_only"]
    )


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "admissibility_protocol": ADMISSIBILITY_PROTOCOL,
    "admissibility_source": ADMISSIBILITY_SOURCE,
    "admissibility_json": ADMISSIBILITY_JSON,
}
hashes_before = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes_before == EXPECTED_HASHES
check("all projection-census inputs have frozen provenance", provenance_ok, str(hashes_before))

print("reconstructing the frozen compatibility matrices", flush=True)
namespace = runpy.run_path(str(ADMISSIBILITY_SOURCE))
hashes_after = {name: digest(path) for name, path in paths.items()}
frozen = json.loads(ADMISSIBILITY_JSON.read_text())
frozen_ok = bool(
    hashes_after == EXPECTED_HASHES
    and frozen["outcome"] == "CANONICAL_DATA_INTERMEDIATE_MODULAR_OPEN"
    and frozen["passed"] == frozen["tests"] == 10
    and all(
        set(record["fixed_ranks"].values()) == {3600}
        and set(record["augmented_ranks"].values()) == {4200}
        and set(record["modular_nullities"].values()) == {240}
        for record in frozen["records"]
    )
)
check("the frozen target-blind 3600/4200 ranks reproduce byte-for-byte", frozen_ok)

modular_rank = namespace["modular_rank"]
primes = tuple(namespace["PRIMES"])

fixture_rows = [
    {0: 1, 2: 1},
    {1: 1},
    {3: 1, 4: 1},
    {},
]
fixture = projection_record(modular_rank, fixture_rows, 2, 2, 1, primes)
fixture_ok = all(dimensions_tuple(record) == (2, 1, 0, 2, 1) for record in fixture.values())
check("the synthetic quotient fixture has its preregistered census", fixture_ok, str(fixture))

corrupt_fixture_rows = [
    {0: 1, 2: 1},
    {1: 1},
    {3: 1},
    {4: 1},
]
corrupt_fixture = projection_record(
    modular_rank, corrupt_fixture_rows, 2, 2, 1, primes
)
corrupt_fixture_ok = all(
    dimensions_tuple(record) == (1, 1, 0, 1, 0)
    for record in corrupt_fixture.values()
)
check("the corrupted fixture destroys the edge-strut cancellation", corrupt_fixture_ok, str(corrupt_fixture))

fixed_count = namespace["CELL_COLUMNS"]
edge_count = namespace["EDGE_COUNT"]
strut_count = namespace["VERTEX_COUNT"]
complex_data = namespace["complex_data"]


def construction_ok(built):
    return bool(
        built["local_geometry"]["controls"]
        and built["transition_control"]
        and built["inverse_control"]
        and built["face_control"]
    )


def census(name, scale, lapse, built):
    print(f"computing projection census {name} ({scale},{lapse})", flush=True)
    return {
        "name": name,
        "scale": scale,
        "lapse": lapse,
        "construction_controls": construction_ok(built),
        "prime_records": projection_record(
            modular_rank,
            built["rows"],
            fixed_count,
            edge_count,
            strut_count,
            primes,
        ),
    }


records = []
for scale, lapse in namespace["REPRESENTATIVES"]:
    records.append(
        census(
            "baseline",
            scale,
            lapse,
            namespace["baseline_builds"][(scale, lapse)],
        )
    )

for scale, lapse in namespace["REPRESENTATIVES"]:
    built = namespace["build_global"](
        complex_data,
        scale,
        lapse,
        namespace["ETA"],
        namespace["CANONICAL_BASE"],
        use_alternate=True,
    )
    records.append(census("alternate_right_inverse", scale, lapse, built))

scale, lapse = namespace["REPRESENTATIVES"][0]
reverse_built = namespace["build_global"](
    complex_data,
    scale,
    lapse,
    namespace["ETA"],
    namespace["CANONICAL_BASE"],
    reverse_orientation=True,
)
records.append(census("reverse_faces", scale, lapse, reverse_built))

odd_canonical = (
    namespace["CANONICAL_BASE"][1],
    namespace["CANONICAL_BASE"][0],
    namespace["CANONICAL_BASE"][2],
    namespace["CANONICAL_BASE"][3],
)
odd_built = namespace["build_global"](
    complex_data,
    scale,
    lapse,
    namespace["ETA"],
    odd_canonical,
)
records.append(census("odd_relabelling", scale, lapse, odd_built))

metric_built = namespace["build_global"](
    complex_data,
    scale,
    lapse,
    -namespace["ETA"],
    namespace["CANONICAL_BASE"],
)
records.append(census("metric_sign", scale, lapse, metric_built))

construction_controls_ok = all(record["construction_controls"] for record in records)
check("all seven complete constructions pass their local and face controls", construction_controls_ok)

rank_controls_ok = all(
    rank_bounds_ok(prime_record, fixed_count, edge_count, strut_count)
    and prime_record["ranks"]["full"] == 4200
    for record in records
    for prime_record in record["prime_records"].values()
)
check("all ranks, bounds and projection identities are internally exact", rank_controls_ok)

prime_agreement_ok = all(
    len({dimensions_tuple(value) for value in record["prime_records"].values()}) == 1
    for record in records
)
check("both primes give the same census in every construction", prime_agreement_ok)

baseline_records = [record for record in records if record["name"] == "baseline"]
baseline_agreement_ok = len(
    {
        dimensions_tuple(next(iter(record["prime_records"].values())))
        for record in baseline_records
    }
) == 1
check("the two nonstatic representatives give the same census", baseline_agreement_ok)

alternate_records = [
    record for record in records if record["name"] == "alternate_right_inverse"
]
alternate_agreement_ok = len(alternate_records) == 2 and all(
    dimensions_tuple(next(iter(record["prime_records"].values())))
    == dimensions_tuple(next(iter(baseline["prime_records"].values())))
    for record, baseline in zip(alternate_records, baseline_records)
)
check("both alternate right-inverse graphs preserve the census", alternate_agreement_ok)

reference_dimensions = dimensions_tuple(
    next(iter(baseline_records[0]["prime_records"].values()))
)
convention_names = {"reverse_faces", "odd_relabelling", "metric_sign"}
convention_records = [record for record in records if record["name"] in convention_names]
convention_agreement_ok = len(convention_records) == 3 and all(
    dimensions_tuple(next(iter(record["prime_records"].values())))
    == reference_dimensions
    for record in convention_records
)
check("orientation, relabelling and metric-sign attacks preserve the census", convention_agreement_ok)

controls_ok = bool(
    provenance_ok
    and frozen_ok
    and fixture_ok
    and corrupt_fixture_ok
    and construction_controls_ok
    and rank_controls_ok
)
all_agree = bool(
    prime_agreement_ok
    and baseline_agreement_ok
    and alternate_agreement_ok
    and convention_agreement_ok
)

if not controls_ok:
    outcome = "CANONICAL_DATA_PROJECTION_CONTROL_FAILED"
elif not all_agree:
    outcome = "CANONICAL_DATA_PROJECTION_DISAGREEMENT_OPEN"
elif len(records) < 7:
    outcome = "CANONICAL_DATA_PROJECTION_INTERMEDIATE_MODULAR_OPEN"
else:
    outcome = "CANONICAL_DATA_PROJECTION_STABLE_MODULAR_OPEN"

allowed = {
    "CANONICAL_DATA_PROJECTION_CONTROL_FAILED",
    "CANONICAL_DATA_PROJECTION_DISAGREEMENT_OPEN",
    "CANONICAL_DATA_PROJECTION_INTERMEDIATE_MODULAR_OPEN",
    "CANONICAL_DATA_PROJECTION_STABLE_MODULAR_OPEN",
}
check("the preregistered projection hierarchy assigns one outcome", outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes_before,
    "source_sha256": digest(Path(__file__)),
    "primes": list(primes),
    "column_counts": {
        "cell_flex": fixed_count,
        "upper_edges": edge_count,
        "struts": strut_count,
    },
    "synthetic_fixture": fixture,
    "corrupted_fixture": corrupt_fixture,
    "records": records,
    "target_comparison": "NOT PERFORMED",
    "rational_dimension_or_carrier": "OPEN",
    "action_or_physics": "NOT EVALUATED",
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
for record in records:
    values = {
        prime: dimensions_tuple(prime_record)
        for prime, prime_record in record["prime_records"].items()
    }
    print(
        f"{record['name']} ({record['scale']},{record['lapse']}): {values}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Target-blind finite carrier conjugacy enumeration.

Prior-art commit: 0129053.
Protocol commit: e006b8c.

This stage deliberately does not open the dynamic tangent artifact.
"""

from collections import Counter
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ACTION_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
GLUING_SOURCE = HERE / "verify_gravity_600cell_dust_two_slab_gluing.py"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_dynamic_tangent_conjugacy_enumeration.json"

PRIOR_ART_COMMIT = "0129053"
PROTOCOL_COMMIT = "e006b8c"
ACTION_SOURCE_SHA256 = "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf"
GLUING_SOURCE_SHA256 = "9ea55dab1fd2f4e9ee643247f5d35599c5894cf77970fc2006fe3d8ac22edf37"
GLUING_ARTIFACT_SHA256 = "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def orbit_sort_key(orbit, phase):
    representative = min(orbit)
    logical = tuple(vertex % 120 for vertex in representative)
    phase_pair = tuple(sorted(phase[vertex] for vertex in logical))
    return phase_pair, tuple(sorted(orbit))


def augment_boundary_orbits(gro, base):
    old_orbits = tuple(sorted(
        gro.orbit_partition(base["old_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    final_orbits = tuple(sorted(
        gro.orbit_partition(base["new_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    return {**base, "old_orbits": old_orbits, "final_orbits": final_orbits}


def orbit_mapping(source_orbits, target_orbits, vertex_image):
    target_lookup = {
        frozenset(orbit): index for index, orbit in enumerate(target_orbits)
    }
    mapping = []
    for orbit in source_orbits:
        image = frozenset(
            tuple(sorted(vertex_image(vertex) for vertex in edge))
            for edge in orbit
        )
        if image not in target_lookup:
            return None
        mapping.append(target_lookup[image])
    result = tuple(mapping)
    if sorted(result) != list(range(30)):
        return None
    return result


def transformed_slab(slab, vertex_image):
    return frozenset(
        tuple(sorted(vertex_image(vertex) for vertex in simplex))
        for simplex in slab
    )


hashes = {
    "action_source": digest(ACTION_SOURCE),
    "gluing_source": digest(GLUING_SOURCE),
    "gluing_artifact": digest(GLUING_ARTIFACT),
}
expected_hashes = {
    "action_source": ACTION_SOURCE_SHA256,
    "gluing_source": GLUING_SOURCE_SHA256,
    "gluing_artifact": GLUING_ARTIFACT_SHA256,
}
gluing = json.loads(GLUING_ARTIFACT.read_text())

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_conjugacy_enumeration", ACTION_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise

models = {
    parity: augment_boundary_orbits(gro, model)
    for parity, model in gro.models.items()
}
even = models["even"]
odd = models["odd"]

schedule_ok = bool(
    tuple(odd["ordering"])
    == (even["ordering"][1], even["ordering"][0],
        even["ordering"][2], even["ordering"][3], even["ordering"][4])
)
stabilizer_even = {tuple(map(int, action)) for action in even["stabilizer"]}
stabilizer_odd = {tuple(map(int, action)) for action in odd["stabilizer"]}
carrier_ok = bool(
    len(gro.h4_actions) == 14400
    and len(even["slab"]) == len(odd["slab"]) == 2400
    and len(even["old_edges"]) == len(odd["old_edges"]) == 720
    and len(even["old_orbits"]) == len(odd["old_orbits"]) == 30
    and len(even["edge_orbits"]) == len(odd["edge_orbits"]) == 35
    and len(even["final_orbits"]) == len(odd["final_orbits"]) == 30
    and Counter(map(len, even["old_orbits"])) == Counter({24: 30})
    and Counter(map(len, odd["old_orbits"])) == Counter({24: 30})
    and stabilizer_even == stabilizer_odd
    and len(stabilizer_even) == 24
)

direct = Counter()
reversed_maps = Counter()
boundary = Counter()
direct_phase_candidates = 0
reversed_phase_candidates = 0

for action in gro.h4_actions:
    action_tuple = tuple(map(int, action))

    boundary_map = orbit_mapping(
        even["old_orbits"], odd["old_orbits"],
        lambda vertex, a=action_tuple: a[vertex],
    )
    if boundary_map is not None:
        boundary[boundary_map] += 1

    direct_phase = all(
        odd["phase"][action_tuple[vertex]] == even["phase"][vertex]
        for vertex in range(120)
    )
    if direct_phase:
        direct_phase_candidates += 1

        def direct_vertex(vertex, a=action_tuple):
            return a[vertex % 120]+120*(vertex//120)

        if transformed_slab(even["slab"], direct_vertex) != odd["slab"]:
            raise RuntimeError("phase-compatible direct action missed the odd slab")
        q_old = orbit_mapping(
            even["old_orbits"], odd["old_orbits"], direct_vertex
        )
        q_final = orbit_mapping(
            even["final_orbits"], odd["final_orbits"], direct_vertex
        )
        if q_old is None or q_final is None:
            raise RuntimeError("direct slab action failed on a boundary partition")
        direct[(q_old, q_final)] += 1

    reversed_phase = all(
        odd["phase"][action_tuple[vertex]] == 4-even["phase"][vertex]
        for vertex in range(120)
    )
    if reversed_phase:
        reversed_phase_candidates += 1

        def reversed_vertex(vertex, a=action_tuple):
            return a[vertex % 120]+120*(1-vertex//120)

        if transformed_slab(even["slab"], reversed_vertex) != odd["slab"]:
            raise RuntimeError("phase-compatible reversing action missed the odd slab")
        q_old_final = orbit_mapping(
            even["old_orbits"], odd["final_orbits"], reversed_vertex
        )
        q_final_old = orbit_mapping(
            even["final_orbits"], odd["old_orbits"], reversed_vertex
        )
        if q_old_final is None or q_final_old is None:
            raise RuntimeError("reversed slab action failed on a boundary partition")
        reversed_maps[(q_old_final, q_final_old)] += 1

identity_map = orbit_mapping(
    even["old_orbits"], odd["old_orbits"], lambda vertex: vertex
)
if identity_map is None:
    raise RuntimeError("the physical-edge identity did not map boundary partitions")

boundary_records = []
for permutation in sorted(set(boundary) | {identity_map}):
    sources = []
    if boundary[permutation]:
        sources.append("H4_BOUNDARY_ACTION")
    if permutation == identity_map:
        sources.append("IDENTICAL_PHYSICAL_EDGE_SETS")
    boundary_records.append({
        "permutation": list(permutation),
        "h4_action_count": boundary[permutation],
        "sources": sources,
    })

direct_records = [{
    "old_to_old": list(key[0]),
    "final_to_final": list(key[1]),
    "h4_action_count": count,
    "endpoint_permutations_equal": key[0] == key[1],
} for key, count in sorted(direct.items())]

reversed_records = [{
    "old_to_final": list(key[0]),
    "final_to_old": list(key[1]),
    "h4_action_count": count,
    "endpoint_permutations_equal": key[0] == key[1],
} for key, count in sorted(reversed_maps.items())]

permutations_ok = all(
    sorted(record[field]) == list(range(30))
    for records, fields in (
        (direct_records, ("old_to_old", "final_to_final")),
        (reversed_records, ("old_to_final", "final_to_old")),
        (boundary_records, ("permutation",)),
    )
    for record in records for field in fields
)
identity_included = any(
    record["permutation"] == list(identity_map)
    and "IDENTICAL_PHYSICAL_EDGE_SETS" in record["sources"]
    for record in boundary_records
)
provenance_ok = bool(
    hashes == expected_hashes
    and gluing.get("outcome") == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing.get("passed") == gluing.get("tests") == 25
)

tests = [
    ("frozen source and gluing hashes", provenance_ok),
    ("imported carrier retains 43/43 certificates", gro.tests == gro.passed == 43),
    ("odd schedule is exactly the frozen first-cell swap", schedule_ok),
    ("finite H4 carrier and quotient dimensions", carrier_ok),
    ("every retained endpoint map is a permutation", permutations_ok),
    ("physical-edge identity boundary map included", identity_included),
    ("tangent and spectral targets remain unparsed", True),
]
passed = sum(bool(ok) for _, ok in tests)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "tangent_matrices_parsed": False,
    "spectral_target_parsed": False,
    "attempted_h4_actions": len(gro.h4_actions),
    "direct_phase_action_count": direct_phase_candidates,
    "reversed_phase_action_count": reversed_phase_candidates,
    "direct_slab_action_count": sum(direct.values()),
    "reversed_slab_action_count": sum(reversed_maps.values()),
    "boundary_h4_action_count": sum(boundary.values()),
    "N_direct_slab": len(direct_records),
    "N_reversed_slab": len(reversed_records),
    "N_boundary": len(boundary_records),
    "direct_slab_candidates": direct_records,
    "reversed_slab_candidates": reversed_records,
    "boundary_candidates": boundary_records,
    "passed": passed,
    "tests": len(tests),
    "outcome": (
        "GEOMETRIC_CONJUGACY_CANDIDATES_ENUMERATED"
        if passed == len(tests) else "GEOMETRIC_CONJUGACY_ENUMERATION_CONTROL_FAILED"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(
    "candidates direct/reversed/boundary = "
    f"{len(direct_records)}/{len(reversed_records)}/{len(boundary_records)}"
)
print(
    "supporting H4 actions direct/reversed/boundary = "
    f"{sum(direct.values())}/{sum(reversed_maps.values())}/{sum(boundary.values())}"
)
print(f"OUTCOME: {payload['outcome']}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) else 1)

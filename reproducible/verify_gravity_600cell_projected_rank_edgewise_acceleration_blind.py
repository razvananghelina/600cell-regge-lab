#!/usr/bin/env python3
"""Blind homogeneous Regge calibration on the canonical projected carrier.

Prior-art commit 1921519 and protocol commit ebe3889 precede execution.  This
stage reconstructs two coefficients but deliberately loads no continuum or
projected-red coefficient target.
"""

import ast
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_acceleration_blind.json"
PRIOR_ART_COMMIT = "1921519"
PROTOCOL_COMMIT = "ebe3889"
INPUT_HASHES = {
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "reproducible/verify_gravity_600cell_projected_rank_edgewise_carrier.py":
        "50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23",
    "reproducible/gravity_600cell_projected_rank_edgewise_carrier.json":
        "b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84",
    "reproducible/verify_gravity_600cell_projected_refinement_acceleration_blind.py":
        "e88111adaeb333abf80b68e06e23d7840ef14399238ada9d0f3cd722d7934e50",
}
ACTION_SOURCE = ROOT / "reproducible/verify_gravity_600cell_projected_refinement_acceleration_blind.py"
ETAS_G = (0.04, 0.02, 0.01, 0.005)
ETAS_F = ETAS_G[:3]
A_SENTINELS = (0.0, -1.0, -2.0, -3.0)
PRIMARY_DERIVATIVE_STEP = 2e-5
SECONDARY_DERIVATIVE_STEP = 1e-5
PRIMARY_LAPSE_DERIVATIVE_STEP = 4e-3
SECONDARY_LAPSE_DERIVATIVE_STEP = 8e-3
EXPECTED = {
    "projected_barycentric": {
        "f_vector": (2640, 17040, 28800, 14400),
        "volume": 19.147932918312847,
    },
    "projected_rank_edgewise_2": {
        "f_vector": (19680, 134880, 230400, 115200),
        "volume": 19.583480465413963,
    },
}
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return condition


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def tetrahedra_from_adjacency(adjacency):
    neighbours = [set(np.flatnonzero(adjacency[index] > 0.5))
                  for index in range(len(adjacency))]
    result = []
    for first in range(len(adjacency)):
        for second in sorted(value for value in neighbours[first]
                             if value > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(value for value in common_two
                                if value > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(value for value in common_three
                                     if value > third):
                    result.append((first, second, third, fourth))
    return tuple(result)


def all_simplices(top):
    return tuple(
        tuple(sorted({
            tuple(sorted(face))
            for tetrahedron in top
            for face in combinations(tetrahedron, degree+1)
        }))
        for degree in range(4)
    )


def projected_barycentric(vertices, coarse_cells, coarse_top):
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    indices = {cell: index for index, cell in enumerate(vertex_cells)}
    raw = np.asarray([vertices[list(cell)].mean(axis=0)
                      for cell in vertex_cells])
    positions = raw/np.linalg.norm(raw, axis=1)[:, None]
    top = []
    for tetrahedron in coarse_top:
        for ordering in permutations(tetrahedron):
            flags = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                tetrahedron,
            )
            top.append(tuple(indices[flag] for flag in flags))
    return positions, tuple(top)


DIRECT_CHILDREN = (
    ("v0", "m01", "m02", "m03"),
    ("v1", "m01", "m12", "m13"),
    ("v2", "m02", "m12", "m23"),
    ("v3", "m03", "m13", "m23"),
    ("m01", "m02", "m03", "m13"),
    ("m01", "m02", "m12", "m13"),
    ("m02", "m03", "m13", "m23"),
    ("m02", "m12", "m13", "m23"),
)


def direct_rank_edgewise_two(base_positions, base_top):
    top_keys = []
    vertex_keys = set()
    for chamber in base_top:
        names = {f"v{index}": (chamber[index], chamber[index])
                 for index in range(4)}
        for left, right in combinations(range(4), 2):
            names[f"m{left}{right}"] = tuple(sorted(
                (chamber[left], chamber[right])
            ))
        for child in DIRECT_CHILDREN:
            keys = tuple(sorted(names[name] for name in child))
            top_keys.append(keys)
            vertex_keys.update(keys)
    ordered_keys = tuple(sorted(vertex_keys))
    index = {key: position for position, key in enumerate(ordered_keys)}
    raw = np.asarray([base_positions[list(key)].mean(axis=0)
                      for key in ordered_keys])
    positions = raw/np.linalg.norm(raw, axis=1)[:, None]
    top = tuple(tuple(sorted(index[key] for key in tetrahedron))
                for tetrahedron in top_keys)
    return positions, top


def topology_record(positions, top):
    cells = all_simplices(top)
    face_counts = Counter(tuple(sorted(face)) for tetrahedron in top
                          for face in combinations(tetrahedron, 3))
    f_vector = tuple(len(layer) for layer in cells)
    return {
        "f_vector": list(f_vector),
        "euler_characteristic": sum(
            (-1)**degree*value for degree, value in enumerate(f_vector)
        ),
        "face_incidence_values": sorted(set(face_counts.values())),
        "unique_top": len(set(top)) == len(top),
        "unit_norm_residual": float(np.max(np.abs(
            np.linalg.norm(positions, axis=1)-1
        ))),
    }


def load_frozen_action_definitions():
    """Load definitions only, never the old verifier's executable body."""
    wanted = {
        "mesh_edges", "mesh_faces_with_counts", "triangle_area",
        "tetra_volume", "richardson_four", "richardson_three",
        "affine_root", "affine_residual", "quadratic_dynamic_root",
        "quadratic_residual", "static_lapse_residual", "CellularMesh",
        "mesh_record", "evaluate_mesh",
    }
    tree = ast.parse(ACTION_SOURCE.read_text(), filename=str(ACTION_SOURCE))
    selected = [
        node for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.ClassDef))
        and node.name in wanted
    ]
    found = {node.name for node in selected}
    if found != wanted:
        raise RuntimeError(f"frozen action definitions missing: {wanted-found}")
    namespace = {
        "np": np,
        "Counter": Counter,
        "combinations": combinations,
        "hashlib": hashlib,
        "ETAS_G": ETAS_G,
        "ETAS_F": ETAS_F,
        "A_SENTINELS": A_SENTINELS,
        "PRIMARY_DERIVATIVE_STEP": PRIMARY_DERIVATIVE_STEP,
        "SECONDARY_DERIVATIVE_STEP": SECONDARY_DERIVATIVE_STEP,
        "PRIMARY_LAPSE_DERIVATIVE_STEP": PRIMARY_LAPSE_DERIVATIVE_STEP,
        "SECONDARY_LAPSE_DERIVATIVE_STEP": SECONDARY_LAPSE_DERIVATIVE_STEP,
    }
    module = ast.Module(body=selected, type_ignores=[])
    exec(compile(module, str(ACTION_SOURCE), "exec"), namespace)
    return namespace


print("="*78)
print("BLIND REGGE ACCELERATION ON THE CANONICAL PROJECTED CARRIER")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "the preregistered carrier and action inputs have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "1921519"
    and PROTOCOL_COMMIT == "ebe3889",
    str(actual_hashes),
)

vertices, adjacency, _ = build_600cell()
vertices = vertices/np.linalg.norm(vertices, axis=1)[:, None]
coarse_top = tetrahedra_from_adjacency(adjacency)
coarse_cells = all_simplices(coarse_top)
base_positions, base_top = projected_barycentric(
    vertices, coarse_cells, coarse_top
)
fine_positions, fine_top = direct_rank_edgewise_two(
    base_positions, base_top
)

carriers = {
    "projected_barycentric": (base_positions, base_top),
    "projected_rank_edgewise_2": (fine_positions, fine_top),
}
topologies = {name: topology_record(*carrier)
              for name, carrier in carriers.items()}
for name in carriers:
    expected = EXPECTED[name]
    topology = topologies[name]
    check(
        f"{name} reproduces the frozen closed carrier",
        tuple(topology["f_vector"]) == expected["f_vector"]
        and topology["euler_characteristic"] == 0
        and topology["face_incidence_values"] == [2]
        and topology["unique_top"]
        and topology["unit_norm_residual"] < 2e-12,
        str(topology),
    )

action = load_frozen_action_definitions()
check(
    "only the frozen action definitions were loaded",
    "CellularMesh" in action and "evaluate_mesh" in action
    and "OUTPUT" not in action and "VARIANTS" not in action,
    "the old verifier executable body was not run",
)

records = {}
internal = {}
for name, (positions, top) in carriers.items():
    record, ok = action["evaluate_mesh"](name, positions, top)
    records[name] = record
    internal[name] = bool(ok)
    expected_volume = EXPECTED[name]["volume"]
    check(
        f"{name} reproduces its frozen chordal volume",
        abs(record["volume_bar"]-expected_volume) < 5e-11,
        f"V={record['volume_bar']:.15f}",
    )
    check(
        f"{name} passes every frozen direct-action control",
        internal[name],
        f"a={record['coefficient_audit']['coefficient']:.12g}, "
        f"a_lapse={record['coefficient_audit']['lapse_coefficient']:.12g}",
    )

coefficients_finite = all(np.isfinite(
    record["coefficient_audit"]["coefficient"]
) for record in records.values())
check(
    "both blind coefficients are finite and target-free",
    coefficients_finite,
    str({name: record["coefficient_audit"]["coefficient"]
         for name, record in records.items()}),
)

all_internal = bool(
    provenance_ok and all(internal.values()) and coefficients_finite
    and all(abs(records[name]["volume_bar"]-EXPECTED[name]["volume"])
            < 5e-11 for name in records)
)
outcome = (
    "CANONICAL_CARRIER_ACCELERATION_COEFFICIENTS_DERIVED"
    if all_internal else
    "CANONICAL_CARRIER_ACCELERATION_INTERNAL_FAILURE"
)
check(
    "the frozen hierarchy assigns exactly one blind outcome",
    outcome in {
        "CANONICAL_CARRIER_ACCELERATION_COEFFICIENTS_DERIVED",
        "CANONICAL_CARRIER_ACCELERATION_INTERNAL_FAILURE",
    },
    outcome,
)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "continuum_target_loaded": False,
    "projected_red_coefficients_loaded": False,
    "carrier_builder": "direct four-corner plus four-around-m02--m13 rank split",
    "action_source_mode": "AST definitions only from exact frozen source",
    "topology": topologies,
    "levels": records,
    "blind_coefficients": {
        name: record["coefficient_audit"]["coefficient"]
        for name, record in records.items()
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("="*78)
for name, value in artifact["blind_coefficients"].items():
    print(f"{name}: {value:.12g}")
print(f"Outcome: {outcome}")
print(f"{passed}/{tests} checks passed")
sys.exit(0 if passed == tests and all_internal else 1)

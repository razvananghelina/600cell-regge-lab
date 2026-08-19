#!/usr/bin/env python3
"""Independent direct-split audit of the projected rank-edgewise carrier."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


PRIMARY = HERE / "gravity_600cell_projected_rank_edgewise_carrier.json"
OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_carrier_adversarial.json"
PRIMARY_SHA256 = "b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84"
SOURCE_SHA256 = "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f"
INDEPENDENCE_COMMIT = "88b12a0"
PROTOCOL_COMMIT = "3840156"

tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return condition


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def parity(permutation):
    return sum(
        permutation[left] > permutation[right]
        for left, right in combinations(range(len(permutation)), 2)
    ) % 2


def all_simplices(top):
    return tuple(
        tuple(sorted({
            tuple(sorted(face))
            for tetrahedron in top
            for face in combinations(tetrahedron, degree+1)
        }))
        for degree in range(4)
    )


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


def face_barycenter(face):
    return tuple(Fraction(int(vertex in face), len(face))
                 for vertex in range(4))


def midpoint(left, right):
    return tuple((a+b)/2 for a, b in zip(left, right))


def direct_children(points):
    names = {f"v{index}": point for index, point in enumerate(points)}
    for left, right in combinations(range(4), 2):
        names[f"m{left}{right}"] = midpoint(points[left], points[right])
    return tuple(tuple(sorted(names[name] for name in child))
                 for child in DIRECT_CHILDREN)


def permute_point(point, permutation):
    result = [Fraction(0) for _ in point]
    for old, new in enumerate(permutation):
        result[new] = point[old]
    return tuple(result)


def permute_complex(top, permutation):
    return frozenset(tuple(sorted(
        permute_point(point, permutation) for point in tetrahedron
    )) for tetrahedron in top)


def local_direct_union():
    result = set()
    for ordering in permutations(range(4)):
        flags = tuple(tuple(sorted(ordering[:rank+1])) for rank in range(4))
        points = tuple(face_barycenter(flag) for flag in flags)
        result.update(direct_children(points))
    return frozenset(result)


def direct_red_parent():
    vertices = tuple(tuple(Fraction(int(axis == vertex)) for axis in range(4))
                     for vertex in range(4))
    return frozenset(direct_children(vertices))


def networkx_tetrahedra(adjacency):
    graph = nx.from_numpy_array(adjacency)
    cliques = [tuple(sorted(clique)) for clique in nx.find_cliques(graph)]
    return tuple(sorted(cliques))


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


def direct_global_split(base_positions, base_top):
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
    return ordered_keys, positions, top


def face_incidence(top):
    return Counter(
        tuple(sorted(face))
        for tetrahedron in top
        for face in combinations(tetrahedron, 3)
    )


def cayley_geometry(positions, top):
    tetrahedra = np.asarray(top, dtype=np.int32)
    qualities = []
    volumes = []
    maximum_chord = 0.0
    maximum_sag = 0.0
    chunk = 12000
    for start in range(0, len(tetrahedra), chunk):
        points = positions[tetrahedra[start:start+chunk]]
        difference = points[:, :, None, :]-points[:, None, :, :]
        distance_squared = np.einsum("nija,nija->nij", difference, difference)
        cayley = np.ones((len(points), 5, 5), dtype=float)
        cayley[:, 0, 0] = 0.0
        cayley[:, 1:, 1:] = distance_squared
        determinants = np.linalg.det(cayley)
        local_volumes = np.sqrt(np.maximum(determinants/288.0, 0.0))
        squared_sum = np.asarray([
            distance_squared[:, left, right]
            for left, right in combinations(range(4), 2)
        ]).sum(axis=0)
        local_quality = 12.0*np.power(3.0*local_volumes, 2.0/3.0)/squared_sum
        maximum_chord = max(maximum_chord, float(np.sqrt(
            max(distance_squared[:, left, right].max()
                for left, right in combinations(range(4), 2))
        )))
        maximum_sag = max(maximum_sag, float(np.max(
            1.0-np.linalg.norm(points.mean(axis=1), axis=1)
        )))
        qualities.append(local_quality)
        volumes.append(local_volumes)
    qualities = np.concatenate(qualities)
    volumes = np.concatenate(volumes)
    total = float(volumes.sum())
    target = float(2*np.pi**2)
    return {
        "tetrahedra": int(len(top)),
        "quality_minimum": float(qualities.min()),
        "quality_median": float(np.median(qualities)),
        "quality_maximum": float(qualities.max()),
        "volume_minimum": float(volumes.min()),
        "volume_maximum": float(volumes.max()),
        "volume_total": total,
        "volume_target_2pi2": target,
        "volume_absolute_error": abs(total-target),
        "maximum_chord_length": maximum_chord,
        "maximum_centroid_radial_sag": maximum_sag,
    }


def regular_volume_control():
    points = np.asarray((
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ))
    affine = points[1:]-points[0]
    gram_volume = np.sqrt(np.linalg.det(affine@affine.T))/6.0
    difference = points[:, None, :]-points[None, :, :]
    distance_squared = np.einsum("ija,ija->ij", difference, difference)
    cayley = np.ones((5, 5))
    cayley[0, 0] = 0.0
    cayley[1:, 1:] = distance_squared
    cayley_volume = np.sqrt(np.linalg.det(cayley)/288.0)
    relative = abs(cayley_volume-gram_volume)/gram_volume
    return gram_volume, cayley_volume, relative


def canonical_top_digest(top):
    array = np.asarray(sorted(top), dtype="<i4")
    return sha256(array.tobytes()).hexdigest()


def coordinate_digest(positions):
    return sha256(np.asarray(positions, dtype="<f8").tobytes()).hexdigest()


print("="*78)
print("ADVERSARIAL DIRECT-SPLIT PROJECTED CARRIER AUDIT")
print("="*78)

primary_hash = digest(PRIMARY)
source_hash = digest(ROOT/"commons/cell600.py")
primary = json.loads(PRIMARY.read_text())
provenance_ok = check(
    "the audit has exact frozen primary and source provenance",
    primary_hash == PRIMARY_SHA256 and source_hash == SOURCE_SHA256
    and primary["outcome"] == "PROJECTED_RANK_EDGEWISE_CARRIER_DERIVED"
    and primary["tests"] == primary["passed"] == 16,
    f"primary={primary_hash}, source={source_hash}",
)

vertices, adjacency, _ = build_600cell()
coarse_top = networkx_tetrahedra(adjacency)
coarse_size_set = {len(tetrahedron) for tetrahedron in coarse_top}
source_ok = check(
    "NetworkX independently finds exactly 600 maximal tetrahedra",
    len(coarse_top) == 600 and coarse_size_set == {4},
    f"count={len(coarse_top)}, sizes={coarse_size_set}",
)
coarse_cells = all_simplices(coarse_top)
base_positions, base_top = projected_barycentric(
    vertices, coarse_cells, coarse_top
)
base_cells = all_simplices(base_top)
base_f = tuple(len(layer) for layer in base_cells)
base_ok = check(
    "the independent containment flags reproduce P(sd K)",
    tuple(len(layer) for layer in coarse_cells) == (120, 720, 1200, 600)
    and base_f == (2640, 17040, 28800, 14400)
    and set(face_incidence(base_top).values()) == {2},
    f"base f={base_f}",
)

local_union = local_direct_union()
local_failures = [
    permutation for permutation in permutations(range(4))
    if permute_complex(local_union, permutation) != local_union
]
local_positive_ok = check(
    "the direct rank-selected split has 192 children and full S4 invariance",
    len(local_union) == 192 and not local_failures,
    f"top={len(local_union)}, failures={len(local_failures)}",
)

red = direct_red_parent()
even = [permutation for permutation in permutations(range(4))
        if parity(permutation) == 0]
red_preserving = [permutation for permutation in even
                  if permute_complex(red, permutation) == red]
red_control_ok = check(
    "the unranked fixed-diagonal red split fails full A4 invariance",
    0 < len(red_preserving) < 12,
    f"A4 preserving={len(red_preserving)}, failing={12-len(red_preserving)}",
)

gram_volume, cayley_volume, volume_relative = regular_volume_control()
volume_control_ok = check(
    "Cayley-Menger reproduces a regular tetrahedron volume",
    volume_relative < 1e-13,
    f"Gram={gram_volume:.16g}, Cayley={cayley_volume:.16g}, rel={volume_relative:.3e}",
)

fine_keys, fine_positions, fine_top = direct_global_split(
    base_positions, base_top
)
fine_cells = all_simplices(fine_top)
fine_f = tuple(len(layer) for layer in fine_cells)
fine_incidence = face_incidence(fine_top)
fine_topology_ok = check(
    "the direct global split has the frozen closed f-vector",
    fine_f == (19680, 134880, 230400, 115200)
    and len(set(fine_top)) == len(fine_top)
    and set(fine_incidence.values()) == {2}
    and sum((-1)**degree*value for degree, value in enumerate(fine_f)) == 0,
    f"f={fine_f}, keys={len(fine_keys)}",
)

base_geometry = cayley_geometry(base_positions, base_top)
fine_geometry = cayley_geometry(fine_positions, fine_top)
fields = (
    "quality_minimum", "quality_median", "quality_maximum",
    "volume_minimum", "volume_maximum", "volume_total",
    "volume_target_2pi2", "volume_absolute_error",
    "maximum_chord_length", "maximum_centroid_radial_sag",
)
comparison = {}
comparison_ok = True
for label, actual in (("projected_barycentric", base_geometry),
                      ("projected_rank_edgewise_2", fine_geometry)):
    frozen = primary["levels"][label]["geometry"]
    comparison[label] = {}
    for field in fields:
        difference = abs(actual[field]-frozen[field])
        tolerance = 5e-11*max(1.0, abs(frozen[field]))
        field_ok = difference <= tolerance
        comparison_ok &= field_ok
        comparison[label][field] = {
            "actual": actual[field],
            "primary": frozen[field],
            "absolute_difference": difference,
            "tolerance": tolerance,
            "passes": bool(field_ok),
        }
geometry_comparison_ok = check(
    "all 20 independent Cayley-Menger diagnostics match the primary artifact",
    comparison_ok,
    f"max normalized difference={max(v['absolute_difference']/v['tolerance'] for level in comparison.values() for v in level.values()):.3e}",
)

inequalities_ok = (
    base_geometry["volume_minimum"] > 0
    and base_geometry["quality_minimum"] > 0.25
    and fine_geometry["volume_minimum"] > 0
    and fine_geometry["quality_minimum"] > 0.25
    and fine_geometry["quality_minimum"]
        >= 0.5*base_geometry["quality_minimum"]
    and fine_geometry["maximum_chord_length"]
        < base_geometry["maximum_chord_length"]
    and fine_geometry["maximum_centroid_radial_sag"]
        < base_geometry["maximum_centroid_radial_sag"]
)
geometry_gate_ok = check(
    "the direct carrier independently passes every finite-geometry inequality",
    inequalities_ok,
    f"q={base_geometry['quality_minimum']:.6f}->{fine_geometry['quality_minimum']:.6f}, "
    f"h={base_geometry['maximum_chord_length']:.6f}->{fine_geometry['maximum_chord_length']:.6f}",
)

controls_ok = provenance_ok and source_ok and base_ok and local_positive_ok \
    and red_control_ok and volume_control_ok
actual_ok = fine_topology_ok and geometry_comparison_ok and geometry_gate_ok
if not controls_ok:
    outcome = "ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_CONTROL_FAILED"
elif not actual_ok:
    outcome = "ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_DISAGREEMENT_OPEN"
else:
    outcome = "ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_CORROBORATED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one adversarial outcome",
    outcome in {
        "ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_CONTROL_FAILED",
        "ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_DISAGREEMENT_OPEN",
        "ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_CORROBORATED",
    },
    outcome,
)

payload = {
    "classification": "DERIVED_COMPUTATIONAL_INDEPENDENT_CARRIER_AUDIT",
    "independence_commit": INDEPENDENCE_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "primary_sha256": primary_hash,
    "source_sha256": source_hash,
    "source_enumerator": "networkx.find_cliques",
    "split_constructor": "direct four corners plus four around m02--m13",
    "volume_constructor": "Cayley-Menger",
    "local_direct_union_tetrahedra": len(local_union),
    "local_S4_failures": len(local_failures),
    "unranked_red_A4_preserving": len(red_preserving),
    "unranked_red_A4_failing": 12-len(red_preserving),
    "base_f_vector": list(base_f),
    "fine_f_vector": list(fine_f),
    "fine_topology_sha256": canonical_top_digest(fine_top),
    "fine_coordinates_sha256": coordinate_digest(fine_positions),
    "base_geometry": base_geometry,
    "fine_geometry": fine_geometry,
    "geometry_comparison": comparison,
    "independence_scope": (
        "direct split, NetworkX cliques and Cayley-Menger geometry; "
        "does not independently reconstruct the full quaternionic H4 action"
    ),
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print("="*78)
print(f"Direct fine topology SHA-256: {payload['fine_topology_sha256']}")
print(f"Direct fine coordinate SHA-256: {payload['fine_coordinates_sha256']}")
print(f"Outcome: {outcome}")
print(f"{passed}/{tests} checks passed")
raise SystemExit(0 if passed == tests and outcome_ok else 1)

#!/usr/bin/env python3
"""Canonical projected rank-edgewise carrier for 600-cell Regge refinement.

The prior-art gate is commit b361eca.  Protocol f00a0bd, corrected before
execution in 251997a, fixes the object, counts, thresholds and outcome
hierarchy.  This verifier contains no Lorentzian action or physics target.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_carrier.json"
PRIOR_ART_COMMIT = "b361eca"
PROTOCOL_COMMIT = "f00a0bd"
PROTOCOL_CORRECTION_COMMIT = "251997a"
INPUT_HASHES = {
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "reproducible/verify_whitney_rank_edgewise_refinement.py":
        "371e28210fcf54f87acba114d26e6ffa8e72246842ef75000b832b4a6847e5dc",
    "reproducible/whitney_rank_edgewise_refinement.json":
        "af0c615a104ee7e22e0f003bde249d90a1b661b5e3cabec242a5763cc44aa77f",
    "reproducible/verify_smooth_hopf_refinement_blind.py":
        "2c0d0cb2ee1a9b1b6e0d3df7c35031f07e49a16850ecd02c8536b5b371f92b8a",
    "reproducible/smooth_hopf_refinement_blind.json":
        "7258a4755ac32af9d32d2415c09bf65f7b1c4a064475a0e2da304f43eb362ba8",
}

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


def all_simplices(top):
    return tuple(
        tuple(sorted({
            tuple(sorted(face))
            for tetrahedron in top
            for face in combinations(tetrahedron, degree + 1)
        }))
        for degree in range(4)
    )


def parity(permutation):
    return sum(
        permutation[left] > permutation[right]
        for left, right in combinations(range(len(permutation)), 2)
    ) % 2


def quaternion_multiply(left, right):
    w1, x1, y1, z1 = np.asarray(left)
    w2, x2, y2, z2 = np.asarray(right)
    return np.array((
        w1*w2-x1*x2-y1*y2-z1*z2,
        w1*x2+x1*w2+y1*z2-z1*y2,
        w1*y2-x1*z2+y1*w2+z1*x2,
        w1*z2+x1*y2-y1*x2+z1*w2,
    ))


def multiplication_table(vertices):
    table = np.empty((len(vertices), len(vertices)), dtype=np.int16)
    maximum_residual = 0.0
    for left, quaternion in enumerate(vertices):
        products = np.asarray([
            quaternion_multiply(quaternion, right) for right in vertices
        ])
        matches = np.argmax(products @ vertices.T, axis=1)
        residuals = np.linalg.norm(products-vertices[matches], axis=1)
        table[left] = matches
        maximum_residual = max(maximum_residual, float(residuals.max()))
    return table, maximum_residual


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


PAIRINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def canonical_pairing(pairs):
    return tuple(sorted(tuple(sorted(pair)) for pair in pairs))


def pairing_action(pairing, permutation):
    return canonical_pairing(
        (permutation[left], permutation[right]) for left, right in pairing
    )


def cell_stabilizers(vertices, table, tetrahedron):
    conjugates = []
    for value in vertices:
        target = value.copy()
        target[1:] *= -1
        index = int(np.argmax(vertices @ target))
        conjugates.append(index)
    conjugates = np.asarray(conjugates, dtype=np.int16)
    inverse = conjugates
    positions = {vertex: index for index, vertex in enumerate(tetrahedron)}
    target_set = set(tetrahedron)
    rotational = set()
    reflected = set()
    for left in range(len(vertices)):
        for right in range(len(vertices)):
            images = tuple(int(table[left, table[vertex, inverse[right]]])
                           for vertex in tetrahedron)
            if set(images) == target_set:
                rotational.add(tuple(positions[value] for value in images))
            conjugated_images = tuple(
                int(table[left, table[conjugates[vertex], inverse[right]]])
                for vertex in tetrahedron
            )
            if set(conjugated_images) == target_set:
                reflected.add(tuple(positions[value]
                                    for value in conjugated_images))
    return rotational, rotational | reflected, conjugates


def weak_compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in weak_compositions(total-first, parts-1):
            yield (first,) + rest


def edgewise_facets(k, dimension=3):
    """Independent exact Edelsbrunner--Grayson color-scheme enumeration."""
    width = dimension + 1
    facets = set()
    for counts in weak_compositions(k*width, width):
        sequence = tuple(
            color for color, count in enumerate(counts)
            for _ in range(count)
        )
        rows = tuple(sequence[row*width:(row+1)*width]
                     for row in range(k))
        columns = tuple(tuple(rows[row][column] for row in range(k))
                        for column in range(width))
        if len(set(columns)) != width:
            continue
        facets.add(tuple(
            tuple(column.count(color) for color in range(width))
            for column in columns
        ))
    return tuple(sorted(facets))


def determinant_three(columns):
    matrix = np.asarray(columns, dtype=np.int64).T
    return int(round(np.linalg.det(matrix)))


def face_barycenter(face):
    return tuple(Fraction(int(vertex in face), len(face))
                 for vertex in range(4))


def combine_fractional(numerator, rank_vertices, k):
    return tuple(sum(
        Fraction(numerator[rank], k)*rank_vertices[rank][axis]
        for rank in range(4)
    ) for axis in range(4))


def local_rank_edgewise_union(local_facets):
    top = set()
    for ordering in permutations(range(4)):
        flags = tuple(tuple(sorted(ordering[:rank+1])) for rank in range(4))
        rank_vertices = tuple(face_barycenter(face) for face in flags)
        for facet in local_facets:
            points = tuple(combine_fractional(point, rank_vertices, 2)
                           for point in facet)
            top.add(tuple(sorted(points)))
    return frozenset(top)


def permute_fractional_point(point, permutation):
    result = [Fraction(0) for _ in point]
    for old, new in enumerate(permutation):
        result[new] = point[old]
    return tuple(result)


def permute_fractional_complex(top, permutation):
    return frozenset(tuple(sorted(
        permute_fractional_point(point, permutation)
        for point in tetrahedron
    )) for tetrahedron in top)


def projected_barycentric_carrier(vertices, coarse_cells, coarse_top):
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    index = {cell: position for position, cell in enumerate(vertex_cells)}
    raw = np.asarray([vertices[list(cell)].mean(axis=0)
                      for cell in vertex_cells])
    positions = raw/np.linalg.norm(raw, axis=1)[:, None]
    top = []
    for tetrahedron in coarse_top:
        for ordering in permutations(tetrahedron):
            flag = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                tetrahedron,
            )
            top.append(tuple(index[cell] for cell in flag))
    return vertex_cells, index, positions, tuple(top)


def projected_edgewise_two(base_positions, base_top, local_facets):
    keys = {}
    top = []

    def point_key(chamber, numerator):
        entries = []
        for rank, weight in enumerate(numerator):
            entries.extend([chamber[rank]]*weight)
        return tuple(sorted(entries))

    for chamber in base_top:
        for facet in local_facets:
            child = []
            for numerator in facet:
                key = point_key(chamber, numerator)
                if key not in keys:
                    keys[key] = len(keys)
                child.append(keys[key])
            top.append(tuple(sorted(child)))

    ordered_keys = [None]*len(keys)
    for key, index in keys.items():
        ordered_keys[index] = key
    raw = np.asarray([
        base_positions[list(key)].mean(axis=0) for key in ordered_keys
    ])
    positions = raw/np.linalg.norm(raw, axis=1)[:, None]
    return tuple(ordered_keys), positions, tuple(top)


def geometry_diagnostics(positions, top):
    tetrahedra = np.asarray(top, dtype=np.int32)
    qualities = []
    volumes = []
    maximum_chord = 0.0
    maximum_sag = 0.0
    edge_pairs = tuple(combinations(range(4), 2))
    chunk = 20000
    for start in range(0, len(tetrahedra), chunk):
        points = positions[tetrahedra[start:start+chunk]]
        affine = points[:, 1:, :] - points[:, :1, :]
        gram = np.einsum("nia,nja->nij", affine, affine)
        determinants = np.linalg.det(gram)
        local_volumes = np.sqrt(np.maximum(determinants, 0.0))/6.0
        squared_sum = np.zeros(len(points))
        for left, right in edge_pairs:
            squared = np.einsum(
                "ni,ni->n", points[:, left]-points[:, right],
                points[:, left]-points[:, right]
            )
            squared_sum += squared
            maximum_chord = max(maximum_chord,
                                float(np.sqrt(squared.max())))
        local_quality = 12.0*np.power(3.0*local_volumes, 2.0/3.0)/squared_sum
        centroids = points.mean(axis=1)
        maximum_sag = max(maximum_sag, float(
            np.max(1.0-np.linalg.norm(centroids, axis=1))))
        qualities.append(local_quality)
        volumes.append(local_volumes)
    qualities = np.concatenate(qualities)
    volumes = np.concatenate(volumes)
    total_volume = float(volumes.sum())
    target_volume = 2.0*np.pi**2
    return {
        "tetrahedra": int(len(top)),
        "quality_minimum": float(qualities.min()),
        "quality_median": float(np.median(qualities)),
        "quality_maximum": float(qualities.max()),
        "volume_minimum": float(volumes.min()),
        "volume_maximum": float(volumes.max()),
        "volume_total": total_volume,
        "volume_target_2pi2": float(target_volume),
        "volume_absolute_error": abs(total_volume-target_volume),
        "maximum_chord_length": maximum_chord,
        "maximum_centroid_radial_sag": maximum_sag,
    }


def face_incidence(top):
    return Counter(
        tuple(sorted(face))
        for tetrahedron in top
        for face in combinations(tetrahedron, 3)
    )


def equivariance_audit(vertices, table, conjugates, vertex_cells,
                       cell_index, positions, base_top):
    top_set = set(base_top)
    maximum_coordinate_residual = 0.0
    topology_failures = 0
    actions = []
    for group_index in range(len(vertices)):
        actions.append(("left", group_index, table[group_index].astype(int)))
        actions.append(("right", group_index, table[:, group_index].astype(int)))
    actions.append(("conjugation", -1, conjugates.astype(int)))

    for kind, group_index, permutation in actions:
        mapped_cells = []
        valid = True
        for cell in vertex_cells:
            image = tuple(sorted(int(permutation[vertex]) for vertex in cell))
            if image not in cell_index:
                valid = False
                break
            mapped_cells.append(cell_index[image])
        if not valid:
            topology_failures += 1
            continue
        mapped_cells = np.asarray(mapped_cells, dtype=np.int32)
        if any(tuple(mapped_cells[list(tetrahedron)]) not in top_set
               for tetrahedron in base_top):
            topology_failures += 1

        if kind == "left":
            transformed = np.asarray([
                quaternion_multiply(vertices[group_index], point)
                for point in positions
            ])
        elif kind == "right":
            transformed = np.asarray([
                quaternion_multiply(point, vertices[group_index])
                for point in positions
            ])
        else:
            transformed = positions.copy()
            transformed[:, 1:] *= -1
        maximum_coordinate_residual = max(
            maximum_coordinate_residual,
            float(np.linalg.norm(
                transformed-positions[mapped_cells], axis=1
            ).max()),
        )
    return {
        "actions_tested": len(actions),
        "topology_failures": topology_failures,
        "maximum_coordinate_residual": maximum_coordinate_residual,
    }


print("="*78)
print("CANONICAL PROJECTED RANK-EDGEWISE 600-CELL CARRIER")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "all preregistered inputs have exact frozen provenance",
    actual_hashes == INPUT_HASHES,
    str(actual_hashes),
)

vertices, adjacency, _ = build_600cell()
coarse_top = tetrahedra_from_adjacency(adjacency)
coarse_cells = all_simplices(coarse_top)
coarse_faces = face_incidence(coarse_top)
source_f = tuple(len(layer) for layer in coarse_cells)
source_ok = check(
    "the source is the closed 600-cell boundary",
    source_f == (120, 720, 1200, 600)
    and np.all(adjacency.sum(axis=1) == 12)
    and set(coarse_faces.values()) == {2}
    and sum((-1)**degree*value for degree, value in enumerate(source_f)) == 0,
    f"f={source_f}, norm residual={np.max(abs(np.linalg.norm(vertices, axis=1)-1)):.3e}",
)

table, multiplication_residual = multiplication_table(vertices)
identity = int(np.argmax(vertices @ np.array((1.0, 0.0, 0.0, 0.0))))
group_ok = check(
    "the 120 vertices close as a quaternion group",
    multiplication_residual < 2e-8
    and all(len(set(row)) == 120 for row in table)
    and np.all(table[identity] == np.arange(120))
    and np.all(table[:, identity] == np.arange(120)),
    f"identity={identity}, max product residual={multiplication_residual:.3e}",
)

rotational, full_stabilizer, conjugates = cell_stabilizers(
    vertices, table, coarse_top[0]
)
even_permutations = set(permutation for permutation in permutations(range(4))
                        if parity(permutation) == 0)
stabilizer_ok = check(
    "the explicit 600-cell tetrahedron stabilizer induces A4 and S4",
    rotational == even_permutations
    and full_stabilizer == set(permutations(range(4))),
    f"rotational={len(rotational)}, full={len(full_stabilizer)}",
)

pairing_orbits = [
    {pairing_action(pairing, permutation) for permutation in rotational}
    for pairing in PAIRINGS
]
fixed_pairings = [
    pairing for pairing in PAIRINGS
    if all(pairing_action(pairing, permutation) == pairing
           for permutation in rotational)
]
red_nogo_ok = check(
    "A4 has no fixed central-octahedron diagonal",
    {len(orbit) for orbit in pairing_orbits} == {3} and not fixed_pairings,
    f"orbit sizes={[len(orbit) for orbit in pairing_orbits]}, fixed={fixed_pairings}",
)

old_rank = json.loads((HERE/"whitney_rank_edgewise_refinement.json").read_text())
old_hopf = json.loads((HERE/"smooth_hopf_refinement_blind.json").read_text())
frozen_control_ok = check(
    "the frozen local rank-edgewise and projected-barycentric controls agree",
    old_rank["direct_edgewise_variants"] == 3
    and old_rank["direct_variant_A4_orbit_sizes"] == [3, 3, 3]
    and old_hopf["dimensions"]["fine_vertices"] == 2640
    and old_hopf["dimensions"]["fine_tetrahedra"] == 14400,
)

vertex_cells, cell_index, base_positions, base_top = (
    projected_barycentric_carrier(vertices, coarse_cells, coarse_top)
)
base_cells = all_simplices(base_top)
base_incidence = face_incidence(base_top)
base_f = tuple(len(layer) for layer in base_cells)
base_top_unique = len(set(base_top)) == len(base_top)
base_topology_ok = check(
    "P(sd K) has the complete closed canonical f-vector",
    base_f == (2640, 17040, 28800, 14400)
    and base_top_unique and set(base_incidence.values()) == {2}
    and sum((-1)**degree*value for degree, value in enumerate(base_f)) == 0
    and np.max(abs(np.linalg.norm(base_positions, axis=1)-1)) < 2e-12,
    f"f={base_f}, unique top={base_top_unique}",
)

local_facets = edgewise_facets(2)
local_vertices = {point for facet in local_facets for point in facet}
local_volume_ratios = set()
for facet in local_facets:
    points = [np.asarray(point[1:], dtype=np.int64) for point in facet]
    columns = [points[index]-points[0] for index in range(1, 4)]
    local_volume_ratios.add(Fraction(abs(determinant_three(columns)), 8))
local_ok = check(
    "the independent k=2 color schemes give eight equal-volume children",
    len(local_facets) == 8 and len(local_vertices) == 10
    and local_volume_ratios == {Fraction(1, 8)},
    f"top={len(local_facets)}, vertices={len(local_vertices)}, volumes={local_volume_ratios}",
)

local_union = local_rank_edgewise_union(local_facets)
local_failures = [
    permutation for permutation in permutations(range(4))
    if permute_fractional_complex(local_union, permutation) != local_union
]
local_equivariance_ok = check(
    "rank-barycentric Esd2 has 192 children and full local S4 invariance",
    len(local_union) == 192 and not local_failures,
    f"top={len(local_union)}, failures={len(local_failures)}",
)

fine_keys, fine_positions, fine_top = projected_edgewise_two(
    base_positions, base_top, local_facets
)
fine_cells = all_simplices(fine_top)
fine_incidence = face_incidence(fine_top)
fine_f = tuple(len(layer) for layer in fine_cells)
fine_top_unique = len(set(fine_top)) == len(fine_top)
fine_merge_ok = check(
    "global Esd2 merges exactly old vertices and shared edge midpoints",
    len(fine_keys) == base_f[0]+base_f[1]
    and sum(left == right for left, right in fine_keys) == base_f[0],
    f"keys={len(fine_keys)}, old={sum(left == right for left, right in fine_keys)}",
)
fine_topology_ok = check(
    "projected Esd2(sd K) has the preregistered closed f-vector",
    fine_f == (19680, 134880, 230400, 115200)
    and fine_top_unique and set(fine_incidence.values()) == {2}
    and sum((-1)**degree*value for degree, value in enumerate(fine_f)) == 0
    and np.max(abs(np.linalg.norm(fine_positions, axis=1)-1)) < 2e-12,
    f"f={fine_f}, unique top={fine_top_unique}",
)

equivariance = equivariance_audit(
    vertices, table, conjugates, vertex_cells, cell_index,
    base_positions, base_top,
)
equivariance_ok = check(
    "left/right 2I actions and conjugation preserve the projected carrier",
    equivariance["actions_tested"] == 241
    and equivariance["topology_failures"] == 0
    and equivariance["maximum_coordinate_residual"] < 2e-8,
    str(equivariance),
)

base_geometry = geometry_diagnostics(base_positions, base_top)
fine_geometry = geometry_diagnostics(fine_positions, fine_top)
base_control_ok = check(
    "the projected barycentric geometry reproduces its frozen chord control",
    abs(base_geometry["maximum_chord_length"]-0.385707678423) < 2e-12,
    f"h={base_geometry['maximum_chord_length']:.12f}",
)
base_geometry_ok = check(
    "the projected barycentric carrier passes its finite shape gate",
    base_geometry["volume_minimum"] > 0
    and base_geometry["quality_minimum"] > 0.25,
    f"q_min={base_geometry['quality_minimum']:.6f}, "
    f"sag={base_geometry['maximum_centroid_radial_sag']:.6e}",
)
fine_geometry_ok = check(
    "the projected edgewise carrier passes the frozen geometry gates",
    fine_geometry["volume_minimum"] > 0
    and fine_geometry["quality_minimum"] > 0.25
    and fine_geometry["quality_minimum"]
        >= 0.5*base_geometry["quality_minimum"]
    and fine_geometry["maximum_chord_length"]
        < base_geometry["maximum_chord_length"]
    and fine_geometry["maximum_centroid_radial_sag"]
        < base_geometry["maximum_centroid_radial_sag"],
    f"q_min={fine_geometry['quality_minimum']:.6f}, "
    f"h={fine_geometry['maximum_chord_length']:.6f}, "
    f"sag={fine_geometry['maximum_centroid_radial_sag']:.6e}",
)

volume_label = (
    "IMPROVES" if fine_geometry["volume_absolute_error"]
    < base_geometry["volume_absolute_error"] else "DOES_NOT_IMPROVE"
)
print(f"[INFO] chordal-volume comparison: {volume_label}")
print(f"       errors={base_geometry['volume_absolute_error']:.6e} -> "
      f"{fine_geometry['volume_absolute_error']:.6e}")

control_ok = provenance_ok and source_ok and group_ok and frozen_control_ok \
    and base_control_ok
canonicity_ok = stabilizer_ok and red_nogo_ok and base_topology_ok \
    and local_ok and local_equivariance_ok and fine_merge_ok \
    and fine_topology_ok and equivariance_ok
geometry_ok = base_geometry_ok and fine_geometry_ok
if not control_ok:
    outcome = "PROJECTED_RANK_EDGEWISE_CONTROL_FAILED"
elif not canonicity_ok:
    outcome = "PROJECTED_RANK_EDGEWISE_CANONICITY_FAILED"
elif not geometry_ok:
    outcome = "PROJECTED_RANK_EDGEWISE_FINITE_GEOMETRY_FAILED"
else:
    outcome = "PROJECTED_RANK_EDGEWISE_CARRIER_DERIVED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one carrier outcome",
    outcome in {
        "PROJECTED_RANK_EDGEWISE_CONTROL_FAILED",
        "PROJECTED_RANK_EDGEWISE_CANONICITY_FAILED",
        "PROJECTED_RANK_EDGEWISE_FINITE_GEOMETRY_FAILED",
        "PROJECTED_RANK_EDGEWISE_CARRIER_DERIVED",
    },
    outcome,
)

payload = {
    "classification": "DERIVED_COMPUTATIONAL_CARRIER_WITH_STRUCTURAL_ROUND_BACKGROUND",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "input_sha256": actual_hashes,
    "source_f_vector": list(source_f),
    "quaternion_group": {
        "order": 120,
        "maximum_product_residual": multiplication_residual,
        "rotational_cell_stabilizer_permutations": len(rotational),
        "full_cell_stabilizer_permutations": len(full_stabilizer),
    },
    "red_diagonal_no_go": {
        "pairing_orbit_sizes": [len(orbit) for orbit in pairing_orbits],
        "fixed_pairings": len(fixed_pairings),
        "conclusion": "NO_H4_EQUIVARIANT_RED_DIAGONAL_ASSIGNMENT",
    },
    "local_edgewise": {
        "k": 2,
        "one_rank_chamber_tetrahedra": len(local_facets),
        "one_rank_chamber_vertices": len(local_vertices),
        "parent_union_tetrahedra": len(local_union),
        "S4_failures": len(local_failures),
    },
    "levels": {
        "projected_barycentric": {
            "f_vector": list(base_f),
            "geometry": base_geometry,
        },
        "projected_rank_edgewise_2": {
            "f_vector": list(fine_f),
            "geometry": fine_geometry,
        },
    },
    "equivariance": equivariance,
    "volume_comparison": volume_label,
    "round_background_status": "STRUCTURAL",
    "infinite_projected_tower_status": "OPEN",
    "lorentzian_action_loaded": False,
    "continuum_dispersion_target_loaded": False,
    "particle_target_loaded": False,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print("="*78)
print(f"Base f={base_f}; fine f={fine_f}")
print(f"Base q_min={base_geometry['quality_minimum']:.9f}; "
      f"fine q_min={fine_geometry['quality_minimum']:.9f}")
print(f"Base h={base_geometry['maximum_chord_length']:.9f}; "
      f"fine h={fine_geometry['maximum_chord_length']:.9f}")
print(f"Volume trend: {volume_label}")
print(f"Outcome: {outcome}")
print(f"{passed}/{tests} checks passed")
raise SystemExit(0 if passed == tests and outcome_ok else 1)

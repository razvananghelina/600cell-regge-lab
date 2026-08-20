#!/usr/bin/env python3
"""Exact combinatorial feasibility census for refined canonical Regge maps.

Prior-art commit 883b4e7 and protocol commit b0d42b8 precede this
implementation.  This verifier computes no action derivative or spectrum.
"""

from collections import Counter
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


OUTPUT = HERE / "gravity_600cell_refined_canonical_map_feasibility.json"
PRIOR_ART_COMMIT = "883b4e7"
PROTOCOL_COMMIT = "b0d42b8"
INPUT_HASHES = {
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "reproducible/verify_gravity_600cell_projected_rank_edgewise_carrier.py":
        "50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23",
    "reproducible/gravity_600cell_projected_rank_edgewise_carrier.json":
        "b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84",
    "reproducible/verify_gravity_600cell_projected_rank_edgewise_balanced_slab.py":
        "f59b8fc89106b42077eca281ff3d956a5a5d6fb4be70b73465133035b1ce0f57",
    "reproducible/gravity_600cell_projected_rank_edgewise_balanced_slab.json":
        "0a9e9e796cd671c82f2e428bfa21ba63ccb07fe76867e4553979c3c54b22a0d5",
}
EXPECTED = {
    "projected_barycentric": (2640, 17040, 28800, 14400),
    "projected_rank_edgewise_2": (19680, 134880, 230400, 115200),
}
PAIR5 = np.asarray(tuple(combinations(range(5), 2)), dtype=np.int64)
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
tests = 0
passed = 0


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
    colours = np.asarray([len(cell)-1 for cell in vertex_cells], dtype=np.int8)
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
    return positions, np.asarray(top, dtype=np.int64), colours


def direct_rank_edgewise_two(base_positions, base_top, base_colours):
    top_keys = []
    vertex_keys = set()
    for chamber in base_top:
        names = {f"v{index}": (int(chamber[index]), int(chamber[index]))
                 for index in range(4)}
        for left, right in combinations(range(4), 2):
            names[f"m{left}{right}"] = tuple(sorted(
                (int(chamber[left]), int(chamber[right]))
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
    colours = np.asarray([
        (int(base_colours[key[0]])+int(base_colours[key[1]])) % 4
        for key in ordered_keys
    ], dtype=np.int8)
    top = np.asarray([
        tuple(sorted(index[key] for key in tetrahedron))
        for tetrahedron in top_keys
    ], dtype=np.int64)
    return positions, top, colours


def spatial_record(top, colours):
    cells = all_simplices(tuple(map(tuple, top)))
    face_counts = Counter(
        tuple(sorted(face))
        for tetrahedron in top
        for face in combinations(map(int, tetrahedron), 3)
    )
    colour_words = np.sort(colours[top], axis=1)
    proper = bool(np.all(colour_words == np.arange(4, dtype=np.int8)))
    return {
        "f_vector": tuple(len(layer) for layer in cells),
        "euler_characteristic": sum(
            (-1)**degree*len(layer) for degree, layer in enumerate(cells)
        ),
        "face_incidence_values": tuple(sorted(set(face_counts.values()))),
        "proper_four_colouring": proper,
        "colour_class_sizes": tuple(
            int(np.count_nonzero(colours == colour)) for colour in range(4)
        ),
        "edges": np.asarray(cells[1], dtype=np.int64),
    }


def encode_edges(edges, slab_vertices):
    edges = np.sort(np.asarray(edges, dtype=np.int64), axis=1)
    return edges[:, 0]*np.int64(slab_vertices)+edges[:, 1]


def staircase_pentachora(top, colours, order, vertex_count):
    position = np.empty(4, dtype=np.int8)
    for rank, colour in enumerate(order):
        position[colour] = rank
    local_order = np.argsort(position[colours[top]], axis=1)
    ordered = np.take_along_axis(top, local_order, axis=1)
    bottom = ordered
    upper = ordered+vertex_count
    pieces = (
        np.column_stack((bottom[:, 0], upper[:, 0], upper[:, 1],
                         upper[:, 2], upper[:, 3])),
        np.column_stack((bottom[:, 0], bottom[:, 1], upper[:, 1],
                         upper[:, 2], upper[:, 3])),
        np.column_stack((bottom[:, 0], bottom[:, 1], bottom[:, 2],
                         upper[:, 2], upper[:, 3])),
        np.column_stack((bottom[:, 0], bottom[:, 1], bottom[:, 2],
                         bottom[:, 3], upper[:, 3])),
    )
    return np.concatenate(pieces, axis=0), position


def schedule_record(top, colours, spatial_edges, order):
    vertex_count = len(colours)
    slab_vertices = 2*vertex_count
    pentachora, position = staircase_pentachora(
        top, colours, order, vertex_count)
    canonical_pentachora = np.sort(pentachora, axis=1)
    distinct_pentachora = len(np.unique(canonical_pentachora, axis=0))

    local_edges = pentachora[:, PAIR5].reshape(-1, 2)
    total_codes = np.unique(encode_edges(local_edges, slab_vertices))
    boundary_edges = np.concatenate((
        spatial_edges,
        spatial_edges+vertex_count,
    ), axis=0)
    boundary_codes = np.unique(encode_edges(boundary_edges, slab_vertices))
    internal_codes = np.setdiff1d(total_codes, boundary_codes,
                                  assume_unique=True)

    edge_positions = position[colours[spatial_edges]]
    lower_first = edge_positions[:, 0] < edge_positions[:, 1]
    lower = np.where(lower_first, spatial_edges[:, 0], spatial_edges[:, 1])
    upper = np.where(lower_first, spatial_edges[:, 1], spatial_edges[:, 0])
    cross_edges = np.column_stack((lower, upper+vertex_count))
    cross_codes = encode_edges(cross_edges, slab_vertices)
    vertical_edges = np.column_stack((
        np.arange(vertex_count, dtype=np.int64),
        np.arange(vertex_count, dtype=np.int64)+vertex_count,
    ))
    vertical_codes = encode_edges(vertical_edges, slab_vertices)
    formula_internal = np.sort(np.concatenate((vertical_codes, cross_codes)))
    formula_total = np.sort(np.concatenate((boundary_codes, formula_internal)))

    return {
        "order": tuple(map(int, order)),
        "position": tuple(map(int, position)),
        "pentachora": int(len(pentachora)),
        "distinct_pentachora": int(distinct_pentachora),
        "total_edges": int(len(total_codes)),
        "boundary_edges": int(len(boundary_codes)),
        "internal_edges": int(len(internal_codes)),
        "direct_formula_edges_equal": bool(
            np.array_equal(total_codes, formula_total)
            and np.array_equal(internal_codes, formula_internal)
        ),
        "internal_digest": sha256(formula_internal.tobytes()).hexdigest(),
        "cross_codes": cross_codes,
    }


def storage_record(dimension):
    byte_count = int(dimension)*int(dimension)*8
    return {
        "dimension": int(dimension),
        "float64_bytes": byte_count,
        "float64_decimal_gb": byte_count/1e9,
        "float64_gib": byte_count/(1024**3),
    }


def level_census(name, top, colours, spatial):
    print(f"[INFO] enumerating 24 staircase schedules on {name}", flush=True)
    vertex_count, edge_count, _, tetrahedron_count = spatial["f_vector"]
    spatial_edges = spatial["edges"]
    colour_pair_counts = Counter(
        tuple(sorted(map(int, colours[edge]))) for edge in spatial_edges
    )
    orders = tuple(permutations(range(4)))
    schedules = [
        schedule_record(top, colours, spatial_edges, order)
        for order in orders
    ]

    direct_formula_ok = all(
        record["direct_formula_edges_equal"]
        and record["pentachora"] == 4*tetrahedron_count
        and record["distinct_pentachora"] == 4*tetrahedron_count
        and record["boundary_edges"] == 2*edge_count
        and record["internal_edges"] == vertex_count+edge_count
        and record["total_edges"] == vertex_count+3*edge_count
        for record in schedules
    )

    distances = []
    distance_formula_ok = True
    for left, right in combinations(range(len(schedules)), 2):
        direct = int(np.count_nonzero(
            schedules[left]["cross_codes"] != schedules[right]["cross_codes"]
        ))
        pos_left = schedules[left]["position"]
        pos_right = schedules[right]["position"]
        formula = sum(
            population
            for (a, b), population in colour_pair_counts.items()
            if ((pos_left[a] < pos_left[b]) != (pos_right[a] < pos_right[b]))
        )
        distances.append(direct)
        distance_formula_ok &= direct == formula

    all_cross = np.concatenate([
        record["cross_codes"] for record in schedules
    ])
    _, multiplicities = np.unique(all_cross, return_counts=True)
    cross_union = int(len(multiplicities))
    cross_intersection = int(np.count_nonzero(multiplicities == len(schedules)))

    phase_dimension = 2*edge_count
    legendre_dimension = vertex_count+2*edge_count
    total_dimension = vertex_count+3*edge_count
    distance_histogram = Counter(distances)
    serial_schedules = [{
        key: value for key, value in record.items()
        if key not in {"position", "cross_codes"}
    } for record in schedules]
    for record, source in zip(serial_schedules, schedules):
        record["position"] = list(source["position"])

    return {
        "spatial_f_vector": list(spatial["f_vector"]),
        "colour_class_sizes": list(spatial["colour_class_sizes"]),
        "colour_pair_edge_populations": {
            f"{a}-{b}": int(value)
            for (a, b), value in sorted(colour_pair_counts.items())
        },
        "schedule_count": len(schedules),
        "distinct_internal_edge_sets": len({
            record["internal_digest"] for record in schedules
        }),
        "cross_edge_union": cross_union,
        "cross_edge_intersection": cross_intersection,
        "expected_cross_edge_union": 2*edge_count,
        "expected_cross_edge_intersection": 0,
        "pairwise_schedule_comparisons": len(distances),
        "weighted_inversion_distance_histogram": {
            str(key): int(value)
            for key, value in sorted(distance_histogram.items())
        },
        "minimum_positive_weighted_distance": min(
            value for value in distances if value > 0),
        "maximum_weighted_distance": max(distances),
        "distance_formula_all_pairs": bool(distance_formula_ok),
        "direct_formula_all_schedules": bool(direct_formula_ok),
        "dimensions": {
            "slab_hessian": storage_record(total_dimension),
            "pre_legendre_jacobian": storage_record(legendre_dimension),
            "boundary_phase_map": storage_record(phase_dimension),
        },
        "pentachora": 4*tetrahedron_count,
        "hessian_upper_incidence_bound": 55*4*tetrahedron_count,
        "schedules": serial_schedules,
    }


print("="*78)
print("REFINED CANONICAL MAP FEASIBILITY CENSUS")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "the frozen carrier and slab inputs have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "883b4e7"
    and PROTOCOL_COMMIT == "b0d42b8",
    str(actual_hashes),
)

carrier_artifact = json.loads(
    (HERE/"gravity_600cell_projected_rank_edgewise_carrier.json").read_text())
slab_artifact = json.loads(
    (HERE/"gravity_600cell_projected_rank_edgewise_balanced_slab.json").read_text())
upstream_ok = check(
    "the frozen upstream artifacts carry accepted carrier/existence controls",
    carrier_artifact["passed"] == carrier_artifact["tests"] == 16
    and carrier_artifact["outcome"]
    == "PROJECTED_RANK_EDGEWISE_CARRIER_DERIVED"
    and slab_artifact["tests"]["passed"]
    == slab_artifact["tests"]["total"] == 15
    and slab_artifact["selection"]["existence_passes"]
    and slab_artifact["selection"]["ordered_slab_alternatives"] == 24,
)

vertices, adjacency, _ = build_600cell()
vertices = vertices/np.linalg.norm(vertices, axis=1)[:, None]
coarse_top = tetrahedra_from_adjacency(adjacency)
coarse_cells = all_simplices(coarse_top)
base_positions, base_top, base_colours = projected_barycentric(
    vertices, coarse_cells, coarse_top)
fine_positions, fine_top, fine_colours = direct_rank_edgewise_two(
    base_positions, base_top, base_colours)
del fine_positions

carriers = {
    "projected_barycentric": (base_top, base_colours),
    "projected_rank_edgewise_2": (fine_top, fine_colours),
}
spatial = {
    name: spatial_record(top, colours)
    for name, (top, colours) in carriers.items()
}
spatial_ok = check(
    "both spatial carriers reproduce topology and proper four-colour controls",
    all(
        record["f_vector"] == EXPECTED[name]
        and record["euler_characteristic"] == 0
        and record["face_incidence_values"] == (2,)
        and record["proper_four_colouring"]
        for name, record in spatial.items()
    ),
    str({name: {key: value for key, value in record.items() if key != "edges"}
         for name, record in spatial.items()}),
)

levels = {
    name: level_census(name, *carriers[name], spatial[name])
    for name in carriers
}
schedule_controls_ok = check(
    "all 48 explicit schedule slabs satisfy the frozen simplex and edge formulas",
    all(
        level["schedule_count"] == 24
        and level["direct_formula_all_schedules"]
        and level["pairwise_schedule_comparisons"] == 276
        and level["distance_formula_all_pairs"]
        for level in levels.values()
    ),
)

boundary_controls_ok = check(
    "the upstream slabs have exactly the two certified spatial boundaries",
    all(
        slab_artifact["levels"][name]["slab"]["boundary_exact"]
        and slab_artifact["levels"][name]["slab"]["bad_facet_incidences"] == 0
        and slab_artifact["levels"][name]["slab"]["bottom_boundary"]
            == EXPECTED[name][3]
        and slab_artifact["levels"][name]["slab"]["top_boundary"]
            == EXPECTED[name][3]
        for name in levels
    ),
)

distinct_schedule_ok = check(
    "all 24 schedules have distinct internal diagonal sets at both levels",
    all(
        level["distinct_internal_edge_sets"] == 24
        and level["cross_edge_union"] == level["expected_cross_edge_union"]
        and level["cross_edge_intersection"]
            == level["expected_cross_edge_intersection"]
        for level in levels.values()
    ),
    str({name: {
        "distinct": level["distinct_internal_edge_sets"],
        "union": level["cross_edge_union"],
        "intersection": level["cross_edge_intersection"],
    } for name, level in levels.items()}),
)

size_ledger_ok = check(
    "the dense sizes and local sparse-support bounds are exact integer ledgers",
    all(
        level["dimensions"]["slab_hessian"]["dimension"]
            == level["spatial_f_vector"][0]+3*level["spatial_f_vector"][1]
        and level["dimensions"]["pre_legendre_jacobian"]["dimension"]
            == level["spatial_f_vector"][0]+2*level["spatial_f_vector"][1]
        and level["dimensions"]["boundary_phase_map"]["dimension"]
            == 2*level["spatial_f_vector"][1]
        and level["hessian_upper_incidence_bound"] == 55*level["pentachora"]
        for level in levels.values()
    ),
)

controls_ok = bool(
    provenance_ok and upstream_ok and spatial_ok and schedule_controls_ok
    and boundary_controls_ok and distinct_schedule_ok and size_ledger_ok
)
if not controls_ok:
    outcome = "REFINED_MAP_FEASIBILITY_CONTROL_FAILED"
elif all(level["distinct_internal_edge_sets"] == 1 for level in levels.values()):
    outcome = "REFINED_MAP_SINGLE_TEMPORAL_CARRIER"
else:
    outcome = "REFINED_MAP_SCHEDULE_ELIMINATION_REQUIRED"

check(
    "the frozen outcome hierarchy assigns the feasibility verdict",
    outcome in {
        "REFINED_MAP_FEASIBILITY_CONTROL_FAILED",
        "REFINED_MAP_SINGLE_TEMPORAL_CARRIER",
        "REFINED_MAP_SCHEDULE_ELIMINATION_REQUIRED",
    },
    outcome,
)

artifact = {
    "title": "Refined canonical map feasibility census",
    "date": "2026-08-20",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_hashes": actual_hashes,
    "definitions": {
        "schedule_orders": "all permutations of four colour classes",
        "time_orientation": "old layer 0 to new layer 1",
        "continuum_target_loaded": False,
        "spectrum_computed": False,
        "action_hessian_computed": False,
        "coarse_fine_transport_selected": False,
    },
    "levels": levels,
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
for name, level in levels.items():
    print(
        f"{name}: schedules={level['distinct_internal_edge_sets']}, "
        f"phase_dim={level['dimensions']['boundary_phase_map']['dimension']}, "
        f"phase_dense={level['dimensions']['boundary_phase_map']['float64_gib']:.3f} GiB, "
        f"upper_support_bound={level['hessian_upper_incidence_bound']}",
        flush=True,
    )
print(f"OUTCOME: {outcome}")
print(f"RESULT: {passed}/{tests} checks passed")

if not controls_ok:
    raise SystemExit(1)

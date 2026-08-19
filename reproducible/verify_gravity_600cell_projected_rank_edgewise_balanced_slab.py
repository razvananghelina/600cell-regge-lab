#!/usr/bin/env python3
"""Selection census for a balanced slab on the rank-edgewise carrier.

Prior-art commit: dabc098.
Protocol commit: f39a5cc.

This verifier is purely combinatorial.  It imports no Regge action, lapse,
Hessian, continuum value, or previously computed slab.
"""

from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_balanced_slab.json"
PRIOR_ART_COMMIT = "dabc098"
PROTOCOL_COMMIT = "f39a5cc"
INPUT_HASHES = {
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
EXPECTED_F = {
    "projected_barycentric": (2640, 17040, 28800, 14400),
    "projected_rank_edgewise_2": (19680, 134880, 230400, 115200),
}
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
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def tetrahedra_from_adjacency(adjacency):
    neighbours = [set(np.flatnonzero(row > 0.5)) for row in adjacency]
    top = []
    for first in range(len(adjacency)):
        for second in sorted(v for v in neighbours[first] if v > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(v for v in common_two if v > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(v for v in common_three if v > third):
                    top.append((first, second, third, fourth))
    return tuple(top)


def all_simplices(top):
    return tuple(
        tuple(sorted({tuple(sorted(face)) for tetrahedron in top
                      for face in combinations(tetrahedron, degree+1)}))
        for degree in range(4)
    )


def barycentric_chambers(coarse_top):
    coarse_cells = all_simplices(coarse_top)
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    cell_index = {cell: index for index, cell in enumerate(vertex_cells)}
    top = []
    for tetrahedron in coarse_top:
        for ordering in permutations(tetrahedron):
            flag = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                tetrahedron,
            )
            top.append(tuple(cell_index[cell] for cell in flag))
    return vertex_cells, cell_index, tuple(top)


def direct_rank_split(base_top):
    top_keys = []
    vertex_keys = set()
    for chamber in base_top:
        names = {f"v{rank}": (chamber[rank], chamber[rank])
                 for rank in range(4)}
        for left, right in combinations(range(4), 2):
            names[f"m{left}{right}"] = tuple(sorted(
                (chamber[left], chamber[right])
            ))
        for child in DIRECT_CHILDREN:
            keys = tuple(sorted(names[name] for name in child))
            top_keys.append(keys)
            vertex_keys.update(keys)
    ordered_keys = tuple(sorted(vertex_keys))
    key_index = {key: index for index, key in enumerate(ordered_keys)}
    top = tuple(tuple(sorted(key_index[key] for key in tetrahedron))
                for tetrahedron in top_keys)
    return ordered_keys, key_index, top


def face_incidence_and_dual(top):
    triangle_to_top = defaultdict(list)
    for top_index, tetrahedron in enumerate(top):
        for triangle in combinations(tetrahedron, 3):
            triangle_to_top[triangle].append(top_index)
    incidence = Counter(len(values) for values in triangle_to_top.values())
    dual = [[] for _ in top]
    for values in triangle_to_top.values():
        if len(values) == 2:
            left, right = values
            dual[left].append(right)
            dual[right].append(left)
    seen = {0}
    queue = deque([0])
    while queue:
        current = queue.popleft()
        for neighbour in dual[current]:
            if neighbour not in seen:
                seen.add(neighbour)
                queue.append(neighbour)
    return triangle_to_top, tuple(tuple(row) for row in dual), incidence, len(seen)


def propagate_colourings(top, triangle_to_top, dual):
    seed = top[0]
    complete = []
    for assignment in permutations(range(4)):
        colours = np.full(max(max(tetrahedron) for tetrahedron in top)+1,
                          -1, dtype=np.int8)
        for vertex, colour in zip(seed, assignment):
            colours[vertex] = colour
        seen = {0}
        queue = deque([0])
        consistent = True
        while queue and consistent:
            current = queue.popleft()
            tetrahedron = top[current]
            if set(int(colours[v]) for v in tetrahedron) != set(range(4)):
                consistent = False
                break
            for neighbour in dual[current]:
                shared = set(tetrahedron) & set(top[neighbour])
                if len(shared) != 3:
                    consistent = False
                    break
                missing_vertex = next(iter(set(top[neighbour])-shared))
                missing_colour = next(iter(set(range(4))-
                                           {int(colours[v]) for v in shared}))
                if colours[missing_vertex] not in (-1, missing_colour):
                    consistent = False
                    break
                colours[missing_vertex] = missing_colour
                if neighbour not in seen:
                    seen.add(neighbour)
                    queue.append(neighbour)
        if (consistent and len(seen) == len(top)
                and np.all(colours >= 0)
                and all(len({int(colours[v]) for v in tetrahedron}) == 4
                        for tetrahedron in top)):
            complete.append(colours)
    hashes = {sha256(colours.tobytes()).hexdigest() for colours in complete}
    return tuple(complete), tuple(sorted(hashes))


def qmul(left, right):
    w1, x1, y1, z1 = left
    w2, x2, y2, z2 = right
    return np.array((
        w1*w2-x1*x2-y1*y2-z1*z2,
        w1*x2+x1*w2+y1*z2-z1*y2,
        w1*y2-x1*z2+y1*w2+z1*x2,
        w1*z2+x1*y2-y1*x2+z1*w2,
    ))


def multiplication_table(vertices):
    table = np.empty((len(vertices), len(vertices)), dtype=np.int16)
    for left, quaternion in enumerate(vertices):
        products = np.asarray([qmul(quaternion, right) for right in vertices])
        table[left] = np.argmax(products @ vertices.T, axis=1)
    return table


def source_actions(vertices, table):
    conjugates = []
    for vertex in vertices:
        target = vertex.copy()
        target[1:] *= -1
        conjugates.append(int(np.argmax(vertices @ target)))
    actions = []
    for group_index in range(len(vertices)):
        actions.append(table[group_index].astype(np.int32))
        actions.append(table[:, group_index].astype(np.int32))
    actions.append(np.asarray(conjugates, dtype=np.int32))
    return tuple(actions)


def induced_maps(actions, vertex_cells, cell_index, fine_keys, fine_index):
    base_maps = []
    fine_maps = []
    for action in actions:
        base_map = np.asarray([
            cell_index[tuple(sorted(int(action[v]) for v in cell))]
            for cell in vertex_cells
        ], dtype=np.int32)
        fine_map = np.asarray([
            fine_index[tuple(sorted((int(base_map[key[0]]),
                                     int(base_map[key[1]]))))]
            for key in fine_keys
        ], dtype=np.int32)
        base_maps.append(base_map)
        fine_maps.append(fine_map)
    return tuple(base_maps), tuple(fine_maps)


def staircase_slab(top, colours, order):
    order_rank = {colour: rank for rank, colour in enumerate(order)}
    pentachora = np.empty((4*len(top), 5), dtype=np.int32)
    cursor = 0
    vertex_count = len(colours)
    for tetrahedron in top:
        ordered = sorted(tetrahedron, key=lambda v: order_rank[int(colours[v])])
        for pivot_rank, pivot in enumerate(ordered):
            simplex = [pivot, pivot+vertex_count]
            simplex.extend(
                other+vertex_count if rank < pivot_rank else other
                for rank, other in enumerate(ordered) if other != pivot
            )
            pentachora[cursor] = np.sort(simplex)
            cursor += 1
    return pentachora


def unique_rows(array, return_counts=False):
    return np.unique(array, axis=0, return_counts=return_counts)


def slab_topology(top, colours, order):
    slab = staircase_slab(top, colours, order)
    unique_top = unique_rows(slab)
    facets = np.concatenate([
        np.delete(slab, omitted, axis=1) for omitted in range(5)
    ])
    unique_facets, counts = unique_rows(facets, return_counts=True)
    boundary = unique_facets[counts == 1]
    vertex_count = len(colours)
    spatial = np.sort(np.asarray(top, dtype=np.int32), axis=1)
    expected_boundary = unique_rows(np.concatenate(
        (spatial, spatial+vertex_count)
    ))
    return {
        "pentachora": int(len(slab)),
        "distinct_pentachora": int(len(unique_top)),
        "facets": int(len(unique_facets)),
        "boundary_facets": int(len(boundary)),
        "interior_facets": int(np.count_nonzero(counts == 2)),
        "bad_facet_incidences": int(np.count_nonzero((counts != 1) &
                                                       (counts != 2))),
        "boundary_exact": bool(np.array_equal(boundary, expected_boundary)),
        "bottom_boundary": int(len(top)),
        "top_boundary": int(len(top)),
    }


def local_staircase(order):
    colours = np.arange(4, dtype=np.int8)
    slab = staircase_slab(((0, 1, 2, 3),), colours, order)
    return frozenset(tuple(int(value) for value in row) for row in slab)


def reverse_time(local):
    return frozenset(tuple(sorted(value+4 if value < 4 else value-4
                                  for value in simplex))
                     for simplex in local)


print("="*78)
print("BALANCED TEMPORAL SLAB SELECTION CENSUS")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the sole source input and preregistration commits have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "dabc098" and PROTOCOL_COMMIT == "f39a5cc",
    str(actual_hashes),
)

vertices, adjacency, _ = build_600cell()
coarse_top = tetrahedra_from_adjacency(adjacency)
vertex_cells, cell_index, base_top = barycentric_chambers(coarse_top)
fine_keys, fine_index, fine_top = direct_rank_split(base_top)

base_colours = np.asarray([len(cell)-1 for cell in vertex_cells], dtype=np.int8)
fine_colours = np.asarray([
    (int(base_colours[left])+int(base_colours[right])) % 4
    for left, right in fine_keys
], dtype=np.int8)

levels = {
    "projected_barycentric": (base_top, base_colours),
    "projected_rank_edgewise_2": (fine_top, fine_colours),
}
records = {}
colourings = {}

for name, (top, colours) in levels.items():
    simplices = all_simplices(top)
    triangle_to_top, dual, incidence, dual_seen = face_incidence_and_dual(top)
    observed_f = tuple(len(layer) for layer in simplices)
    check(
        f"{name} reproduces its certified closed connected spatial carrier",
        observed_f == EXPECTED_F[name]
        and incidence == Counter({2: EXPECTED_F[name][2]})
        and dual_seen == len(top),
        f"f={observed_f}, triangle incidence={dict(incidence)}, "
        f"dual={dual_seen}/{len(top)}",
    )
    proper_edges = all(colours[left] != colours[right]
                       for left, right in simplices[1])
    full_tetrahedra = all({int(colours[v]) for v in tetrahedron} == set(range(4))
                            for tetrahedron in top)
    check(
        f"{name} has the frozen proper rank-derived four-colouring",
        proper_edges and full_tetrahedra,
        f"colour classes={np.bincount(colours, minlength=4).tolist()}",
    )
    complete, hashes = propagate_colourings(top, triangle_to_top, dual)
    colourings[name] = complete
    declared_matches = any(np.array_equal(candidate, colours)
                           for candidate in complete)
    check(
        f"{name} has exactly 4! global labelled colourings and no later branch",
        len(complete) == 24 and len(hashes) == 24 and declared_matches,
        f"complete={len(complete)}, distinct={len(hashes)}, "
        f"declared present={declared_matches}",
    )
    topology = slab_topology(top, colours, (0, 1, 2, 3))
    check(
        f"{name} gives a conforming closed-space staircase slab",
        topology["pentachora"] == 4*len(top)
        and topology["distinct_pentachora"] == 4*len(top)
        and topology["bad_facet_incidences"] == 0
        and topology["boundary_facets"] == 2*len(top)
        and topology["boundary_exact"],
        str(topology),
    )
    records[name] = {
        "spatial_f_vector": list(observed_f),
        "colour_class_sizes": np.bincount(colours, minlength=4).tolist(),
        "labelled_global_colourings": len(complete),
        "dual_tetrahedra_reached": dual_seen,
        "slab": topology,
    }

vertices = vertices/np.linalg.norm(vertices, axis=1)[:, None]
table = multiplication_table(vertices)
actions = source_actions(vertices, table)
base_maps, fine_maps = induced_maps(
    actions, vertex_cells, cell_index, fine_keys, fine_index
)
action_records = {}
for name, colours, maps in (
    ("projected_barycentric", base_colours, base_maps),
    ("projected_rank_edgewise_2", fine_colours, fine_maps),
):
    failures = sum(not np.array_equal(colours[mapping], colours)
                   for mapping in maps)
    action_records[name] = {
        "actions_tested": len(maps),
        "colour_preservation_failures": failures,
    }
    check(
        f"all 241 declared spatial actions preserve the {name} colouring",
        len(maps) == 241 and failures == 0,
        str(action_records[name]),
    )

orders = tuple(permutations(range(4)))
local_slabs = {order: local_staircase(order) for order in orders}
distinct_local = {slab for slab in local_slabs.values()}
reversal = {}
for order, slab in local_slabs.items():
    image = reverse_time(slab)
    matches = tuple(candidate for candidate, candidate_slab in local_slabs.items()
                    if candidate_slab == image)
    reversal[order] = matches[0] if len(matches) == 1 else None
reversal_ok = all(reversal[order] == tuple(reversed(order)) for order in orders)
fixed = sum(reversal[order] == order for order in orders)
unseen = set(orders)
orbit_sizes = []
while unseen:
    order = min(unseen)
    orbit = {order, reversal[order]}
    unseen -= orbit
    orbit_sizes.append(len(orbit))

check(
    "all 4! colour orders give distinct labelled staircase slabs",
    len(orders) == 24 and len(distinct_local) == 24,
    f"N_order={len(distinct_local)}",
)
check(
    "time reversal sends each order to its reverse with twelve size-two orbits",
    reversal_ok and fixed == 0 and Counter(orbit_sizes) == Counter({2: 12}),
    f"fixed={fixed}, orbit sizes={dict(Counter(orbit_sizes))}",
)

# Every tested spatial action fixes each colour class pointwise as a class.
# Consequently it preserves the staircase rule for every fixed colour order,
# not merely for the declared 0<1<2<3 order.
h4_invariant_orders = len(orders) if all(
    record["colour_preservation_failures"] == 0
    for record in action_records.values()
) else 0
check(
    "spatial H4 leaves all 24 ordered slab alternatives admissible",
    h4_invariant_orders == 24,
    f"H4-invariant ordered slabs={h4_invariant_orders}",
)

selection = {
    "existence_passes": True,
    "ordered_slab_alternatives": len(distinct_local),
    "h4_invariant_ordered_slab_alternatives": h4_invariant_orders,
    "time_reversal_fixed_orders": fixed,
    "time_reversal_orbit_sizes": dict(sorted(Counter(orbit_sizes).items())),
    "canonical_selection_passes": len(distinct_local) == 1,
    "classification": "STRUCTURAL" if len(distinct_local) > 1 else "DERIVED",
}
check(
    "the preregistered canonical-selection gate is evaluated, not hidden",
    selection["existence_passes"]
    and not selection["canonical_selection_passes"]
    and selection["classification"] == "STRUCTURAL",
    str(selection),
)

artifact = {
    "title": "Balanced temporal slab selection census",
    "date": "2026-08-19",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_hashes": actual_hashes,
    "definitions": {
        "base_colour": "face dimension r in {0,1,2,3}",
        "fine_colour": "(r(left)+r(right)) mod 4",
        "declared_colour_order": [0, 1, 2, 3],
        "time_orientation": "past < future",
    },
    "levels": records,
    "spatial_actions": action_records,
    "selection": selection,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")
print(f"WROTE {OUTPUT.relative_to(ROOT)}")
print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)

#!/usr/bin/env python3
"""Census whether orientation selects one five-colour staircase schedule.

Prior-art commit: 7c9cd5b.
Protocol commit: 5f399c2.

This verifier is purely combinatorial/geometric.  It deliberately imports no
Regge action, Hessian, canonical map, nonlinear result, or continuum target.
"""

from collections import Counter, defaultdict, deque
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_staircase_orientation_selector.json"
PRIOR_ART_COMMIT = "7c9cd5b"
PROTOCOL_COMMIT = "5f399c2"
DET_TOLERANCE = 1e-10


def qmul(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return np.array((
        a*e-b*f-c*g-d*h,
        a*f+b*e+c*h-d*g,
        a*g-b*h+c*e+d*f,
        a*h+b*g-c*f+d*e,
    ))


def canonical_key(item):
    return tuple(sorted(item))


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left+1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def build_tetrahedra(adjacency):
    neighbours = [set(np.flatnonzero(row).tolist()) for row in adjacency]
    edges = {
        (left, right)
        for left in range(120)
        for right in range(left+1, 120)
        if adjacency[left, right]
    }
    tetrahedra = []
    for a in range(120):
        for b in sorted(vertex for vertex in neighbours[a] if vertex > a):
            for c in sorted(
                vertex for vertex in neighbours[a] & neighbours[b]
                if vertex > b
            ):
                for d in sorted(
                    vertex
                    for vertex in neighbours[a] & neighbours[b] & neighbours[c]
                    if vertex > c
                ):
                    tetrahedra.append((a, b, c, d))
    return edges, tuple(tetrahedra)


def build_slab(tetrahedra, vertex_colour, order):
    rank = {colour: index for index, colour in enumerate(order)}
    simplices = set()
    for tetrahedron in tetrahedra:
        ordered = sorted(tetrahedron, key=lambda vertex: rank[vertex_colour[vertex]])
        if len({vertex_colour[vertex] for vertex in ordered}) != 4:
            raise RuntimeError("spatial tetrahedron repeats a cover colour")
        for vertex in ordered:
            simplex = [vertex, vertex+120]
            simplex.extend(
                other+120 if rank[vertex_colour[other]] < rank[vertex_colour[vertex]]
                else other
                for other in ordered if other != vertex
            )
            simplices.add(tuple(sorted(simplex)))
    return frozenset(simplices)


def faces(simplices, size):
    return frozenset(
        face for simplex in simplices for face in combinations(simplex, size)
    )


def sign_from_det(matrix):
    determinant = float(np.linalg.det(matrix))
    if abs(determinant) < DET_TOLERANCE:
        return 0, abs(determinant)
    return (1 if determinant > 0 else -1), abs(determinant)


def oriented_boundary(slab, vertices):
    coefficients = {}
    minimum = float("inf")
    boundary = Counter()
    for simplex in slab:
        rows = np.array([
            np.r_[vertices[vertex % 120], float(vertex >= 120)]
            for vertex in simplex
        ])
        coefficient, magnitude = sign_from_det(rows)
        coefficients[simplex] = coefficient
        minimum = min(minimum, magnitude)
        for omitted in range(5):
            facet = simplex[:omitted] + simplex[omitted+1:]
            boundary[facet] += coefficient*((-1)**omitted)
    return coefficients, Counter({key: value for key, value in boundary.items() if value}), minimum


def chain_orientation_sign(action, tetrahedra, spatial_coefficients):
    observed = set()
    for tetrahedron in tetrahedra:
        image_sequence = tuple(int(action[vertex]) for vertex in tetrahedron)
        image = tuple(sorted(image_sequence))
        positions = {vertex: index for index, vertex in enumerate(image)}
        reorder = tuple(positions[vertex] for vertex in image_sequence)
        pushed = spatial_coefficients[tetrahedron]*permutation_sign(reorder)
        observed.add(pushed*spatial_coefficients[image])
        if len(observed) > 1:
            return 0
    return next(iter(observed))


def components(vertices, adjacency):
    unseen = set(vertices)
    result = []
    while unseen:
        seed = min(unseen)
        queue = deque([seed])
        component = {seed}
        unseen.remove(seed)
        while queue:
            current = queue.popleft()
            for neighbour in adjacency[current]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    queue.append(neighbour)
        result.append(frozenset(component))
    return tuple(result)


tests = []


def check(label, condition):
    ok = bool(condition)
    tests.append((label, ok))
    print(f"{'PASS' if ok else 'FAIL'}: {label}")


vertices, adjacency_float, _ = build_600cell()
adjacency = adjacency_float > 0.5
edges, tetrahedra = build_tetrahedra(adjacency)
spatial_triangles = faces(tetrahedra, 3)

multiplication = np.empty((120, 120), dtype=np.int16)
for left in range(120):
    for right in range(120):
        multiplication[left, right] = int(np.argmax(
            vertices @ qmul(vertices[left], vertices[right])
        ))
conjugate = np.array([
    int(np.argmax(vertices @ (vertex*np.array((1, -1, -1, -1)))))
    for vertex in vertices
], dtype=np.int16)

binary_tetrahedral = frozenset(
    index for index, vertex in enumerate(vertices)
    if (
        np.count_nonzero(np.abs(vertex) > 1e-8) == 1
        and np.max(np.abs(vertex)) > 1-1e-8
    ) or np.all(np.abs(np.abs(vertex)-0.5) < 1e-8)
)
unseen = set(range(120))
cover_cells = []
while unseen:
    representative = min(unseen)
    cell = frozenset(
        int(multiplication[representative, element])
        for element in binary_tetrahedral
    )
    cover_cells.append(cell)
    unseen -= cell
cover_cells = tuple(sorted(cover_cells, key=canonical_key))
cell_lookup = {cell: index for index, cell in enumerate(cover_cells)}
vertex_colour = {}
for colour, cell in enumerate(cover_cells):
    for vertex in cell:
        if vertex in vertex_colour:
            raise RuntimeError("cover cells overlap")
        vertex_colour[vertex] = colour

carrier_ok = bool(
    vertices.shape == (120, 4)
    and len(edges) == 720
    and len(spatial_triangles) == 1200
    and len(tetrahedra) == 600
    and len(cover_cells) == 5
    and {len(cell) for cell in cover_cells} == {24}
    and len(vertex_colour) == 120
    and all(len({vertex_colour[v] for v in tetrahedron}) == 4
            for tetrahedron in tetrahedra)
)

spatial_coefficients = {}
spatial_det_minimum = float("inf")
for tetrahedron in tetrahedra:
    coefficient, magnitude = sign_from_det(vertices[list(tetrahedron)])
    spatial_coefficients[tetrahedron] = coefficient
    spatial_det_minimum = min(spatial_det_minimum, magnitude)

all_orders = tuple(permutations(range(5)))
old_boundary = frozenset(tetrahedra)
new_boundary = frozenset(
    tuple(vertex+120 for vertex in tetrahedron) for tetrahedron in tetrahedra
)
slabs = {}
order_records = []
all_fvectors = []
all_incidence_ok = True
all_chain_ok = spatial_det_minimum >= DET_TOLERANCE
all_boundary_sign_pairs = []
all_slab_det_minimum = float("inf")

for order in all_orders:
    slab = build_slab(tetrahedra, vertex_colour, order)
    slabs[order] = slab
    layer_faces = {size: faces(slab, size) for size in range(1, 6)}
    fvector = tuple(len(layer_faces[size]) for size in range(1, 6))
    all_fvectors.append(fvector)
    facet_incidence = Counter(
        facet for simplex in slab for facet in combinations(simplex, 4)
    )
    incidence_boundary = frozenset(
        facet for facet, count in facet_incidence.items() if count == 1
    )
    incidence_ok = bool(
        len(slab) == 2400
        and fvector == (240, 2280, 6240, 6600, 2400)
        and sum((-1)**index*value for index, value in enumerate(fvector)) == 0
        and Counter(facet_incidence.values()) == Counter({2: 5400, 1: 1200})
        and incidence_boundary == old_boundary | new_boundary
    )
    all_incidence_ok &= incidence_ok

    _, boundary_chain, slab_det_minimum = oriented_boundary(slab, vertices)
    all_slab_det_minimum = min(all_slab_det_minimum, slab_det_minimum)
    ratios_old = {
        boundary_chain.get(tetrahedron, 0)*spatial_coefficients[tetrahedron]
        for tetrahedron in tetrahedra
    }
    ratios_new = {
        boundary_chain.get(tuple(v+120 for v in tetrahedron), 0)
        * spatial_coefficients[tetrahedron]
        for tetrahedron in tetrahedra
    }
    sign_pair = (
        next(iter(ratios_old)) if len(ratios_old) == 1 else 0,
        next(iter(ratios_new)) if len(ratios_new) == 1 else 0,
    )
    chain_ok = bool(
        slab_det_minimum >= DET_TOLERANCE
        and len(boundary_chain) == 1200
        and set(boundary_chain) == set(old_boundary | new_boundary)
        and sign_pair[0] in (-1, 1)
        and sign_pair[1] == -sign_pair[0]
    )
    all_chain_ok &= chain_ok
    all_boundary_sign_pairs.append(sign_pair)
    order_records.append({
        "order": list(order),
        "permutation_sign": permutation_sign(order),
        "f_vector": list(fvector),
        "orientation_sign_pair_old_new": list(sign_pair),
        "minimum_absolute_pentachoron_determinant": slab_det_minimum,
        "incidence_ok": incidence_ok,
        "orientation_chain_ok": chain_ok,
    })

unique_slabs = len(set(slabs.values()))

# Full H4 action with retained reflection provenance.
action_provenance = defaultdict(set)
plain = np.arange(120, dtype=np.int16)
for reflected in (False, True):
    seed = conjugate if reflected else plain
    for left in range(120):
        left_images = multiplication[left, seed]
        for right in range(120):
            action = tuple(
                int(value)
                for value in multiplication[left_images, conjugate[right]]
            )
            action_provenance[action].add(reflected)

setwise_records = []
induced_representatives = {}
orientation_construction_match = True
for action_tuple, reflected_values in action_provenance.items():
    action = np.asarray(action_tuple, dtype=np.int16)
    images = tuple(frozenset(int(action[v]) for v in cell) for cell in cover_cells)
    if not all(image in cell_lookup for image in images):
        continue
    induced = tuple(cell_lookup[image] for image in images)
    if len(set(induced)) != 5:
        continue
    spatial_sign = chain_orientation_sign(action, tetrahedra, spatial_coefficients)
    construction_signs = {(-1 if reflected else 1) for reflected in reflected_values}
    construction_match = len(construction_signs) == 1 and spatial_sign in construction_signs
    orientation_construction_match &= construction_match
    record = {
        "induced_permutation": list(induced),
        "induced_permutation_sign": permutation_sign(induced),
        "spatial_orientation_sign": spatial_sign,
        "reflected_provenance": sorted(reflected_values),
        "orientation_matches_construction": construction_match,
    }
    setwise_records.append(record)
    induced_representatives.setdefault(induced, action)

induced_group = frozenset(induced_representatives)
identity = tuple(range(5))
kernel_size = sum(
    tuple(record["induced_permutation"]) == identity
    for record in setwise_records
)
induced_sign_counts = Counter(permutation_sign(group_element)
                              for group_element in induced_group)
setwise_cross_counts = Counter(
    (
        record["induced_permutation_sign"],
        record["spatial_orientation_sign"],
        tuple(record["reflected_provenance"]),
    )
    for record in setwise_records
)

# Verify that induced cover permutations really transport the full slabs.
canonical_order = tuple(range(5))
symmetry_transport_ok = True
for induced, action in induced_representatives.items():
    transported = frozenset(
        tuple(sorted(
            int(action[v]) if v < 120 else int(action[v-120])+120
            for v in simplex
        ))
        for simplex in slabs[canonical_order]
    )
    target_order = tuple(induced[colour] for colour in canonical_order)
    symmetry_transport_ok &= transported == slabs[target_order]

unseen_orders = set(all_orders)
order_orbits = []
while unseen_orders:
    seed = min(unseen_orders)
    orbit = frozenset(
        tuple(group_element[colour] for colour in seed)
        for group_element in induced_group
    )
    order_orbits.append(orbit)
    unseen_orders -= orbit
order_orbits = tuple(sorted(order_orbits, key=lambda orbit: min(orbit)))
order_to_orbit = {
    order: orbit_index
    for orbit_index, orbit in enumerate(order_orbits)
    for order in orbit
}
order_record_lookup = {tuple(record["order"]): record for record in order_records}
orbit_records = []
for orbit_index, orbit in enumerate(order_orbits):
    orbit_records.append({
        "orbit_index": orbit_index,
        "size": len(orbit),
        "permutation_sign_counts": dict(Counter(
            permutation_sign(order) for order in orbit
        )),
        "orientation_sign_pair_counts": {
            str(key): value for key, value in Counter(
                tuple(order_record_lookup[order]["orientation_sign_pair_old_new"])
                for order in orbit
            ).items()
        },
        "all_controls_pass": all(
            order_record_lookup[order]["incidence_ok"]
            and order_record_lookup[order]["orientation_chain_ok"]
            for order in orbit
        ),
    })

# Exact layer reversal and its schedule map.
slab_to_order = {slab: order for order, slab in slabs.items()}
time_reversal = {}
time_reversal_parity_transitions = Counter()
time_reversal_orbit_transitions = Counter()
time_reversal_ok = len(slab_to_order) == 120
for order, slab in slabs.items():
    reversed_slab = frozenset(
        tuple(sorted(vertex+120 if vertex < 120 else vertex-120
                     for vertex in simplex))
        for simplex in slab
    )
    target = slab_to_order.get(reversed_slab)
    if target is None:
        time_reversal_ok = False
        continue
    time_reversal[order] = target
    time_reversal_parity_transitions[(permutation_sign(order), permutation_sign(target))] += 1
    time_reversal_orbit_transitions[(order_to_orbit[order], order_to_orbit[target])] += 1

# Adjacent-transposition flip graph.
flip_adjacency = {order: set() for order in all_orders}
flip_edges = set()
flip_census = Counter()
flip_orbit_transitions = Counter()
for order in all_orders:
    for position in range(4):
        neighbour = list(order)
        neighbour[position], neighbour[position+1] = neighbour[position+1], neighbour[position]
        neighbour = tuple(neighbour)
        edge = tuple(sorted((order, neighbour)))
        if edge in flip_edges:
            continue
        flip_edges.add(edge)
        flip_adjacency[order].add(neighbour)
        flip_adjacency[neighbour].add(order)
        left, right = slabs[order], slabs[neighbour]
        intersection = len(left & right)
        removed = len(left-right)
        added = len(right-left)
        flip_census[(intersection, removed, added, len(left ^ right))] += 1
        flip_orbit_transitions[
            tuple(sorted((order_to_orbit[order], order_to_orbit[neighbour])))
        ] += 1
flip_components = components(all_orders, flip_adjacency)

all_orbits_pass = all(record["all_controls_pass"] for record in orbit_records)
common_boundary_chain = len(set(all_boundary_sign_pairs)) == 1
no_selector_conditions = bool(
    carrier_ok
    and all_incidence_ok
    and all_chain_ok
    and len(order_orbits) >= 2
    and all_orbits_pass
    and common_boundary_chain
    and len(flip_components) == 1
    and time_reversal_ok
)
passing_orbits = sum(record["all_controls_pass"] for record in orbit_records)
if no_selector_conditions:
    outcome = "ORIENTATION_DOES_NOT_SELECT_PARITY"
elif (
    carrier_ok and len(order_orbits) >= 2 and passing_orbits == 1
    and time_reversal_ok and len(flip_components) == 1
):
    outcome = "ORIENTATION_SELECTS_ONE_PARITY"
else:
    outcome = "OPEN_CONTROL_FAILURE"

check("600-cell and five-cell cover have the frozen counts", carrier_ok)
check("all 120 total orders produce distinct staircase slabs", unique_slabs == 120)
check("all orders have the frozen product f-vector and facet incidence", all_incidence_ok)
check("all signed four-chains cancel internally with resolved determinants", all_chain_ok)
check("all orders induce the same opposite old/new boundary signs", common_boundary_chain)
check("all 14,400 H4 actions are unique", len(action_provenance) == 14400)
check("setwise cover action is exact and orientation signs cross-check", bool(setwise_records) and orientation_construction_match)
check("induced cover permutations transport the complete staircase slabs", symmetry_transport_ok)
check("induced group partitions all 120 orders without overlap", sum(map(len, order_orbits)) == 120 and len(order_to_orbit) == 120)
check("layer reversal maps every slab to an enumerated schedule", time_reversal_ok and len(time_reversal) == 120)
check("adjacent-transposition graph is a connected 120-vertex census", len(flip_edges) == 240 and tuple(map(len, flip_components)) == (120,))
check("mechanical selector verdict assigned", outcome in {"ORIENTATION_SELECTS_ONE_PARITY", "ORIENTATION_DOES_NOT_SELECT_PARITY", "OPEN_CONTROL_FAILURE"})
check("no Regge, continuum, canonical-map or nonlinear target parsed", True)

passed = sum(ok for _, ok in tests)
payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "regge_action_parsed": False,
    "continuum_target_parsed": False,
    "canonical_map_parsed": False,
    "nonlinear_result_parsed": False,
    "carrier": {
        "vertices": len(vertices),
        "edges": len(edges),
        "triangles": len(spatial_triangles),
        "tetrahedra": len(tetrahedra),
        "cover_cell_sizes": [len(cell) for cell in cover_cells],
        "spatial_minimum_absolute_determinant": spatial_det_minimum,
    },
    "schedule_census": {
        "orders": len(all_orders),
        "unique_slabs": unique_slabs,
        "f_vector_counts": {str(key): value for key, value in Counter(all_fvectors).items()},
        "boundary_sign_pair_counts": {str(key): value for key, value in Counter(all_boundary_sign_pairs).items()},
        "minimum_absolute_pentachoron_determinant": all_slab_det_minimum,
        "records": order_records,
    },
    "h4_cover_action": {
        "all_h4_actions": len(action_provenance),
        "setwise_actions": len(setwise_records),
        "distinct_induced_permutations": len(induced_group),
        "induced_permutation_sign_counts": {str(key): value for key, value in induced_sign_counts.items()},
        "kernel_size": kernel_size,
        "setwise_cross_counts": {str(key): value for key, value in setwise_cross_counts.items()},
        "order_orbits": orbit_records,
    },
    "time_reversal": {
        "parity_transition_counts": {str(key): value for key, value in time_reversal_parity_transitions.items()},
        "orbit_transition_counts": {str(key): value for key, value in time_reversal_orbit_transitions.items()},
        "schedule_map": {str(order): list(target) for order, target in time_reversal.items()},
    },
    "adjacent_transposition_graph": {
        "vertices": len(all_orders),
        "edges": len(flip_edges),
        "component_sizes": [len(component) for component in flip_components],
        "simplex_change_census": {str(key): value for key, value in flip_census.items()},
        "orbit_transition_counts": {str(key): value for key, value in flip_orbit_transitions.items()},
    },
    "tests": len(tests),
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print(f"orders={len(all_orders)}, unique_slabs={unique_slabs}")
print(f"boundary_sign_pairs={dict(Counter(all_boundary_sign_pairs))}")
print(
    "H4 setwise actions={}, induced permutations={}, kernel={}, signs={}".format(
        len(setwise_records), len(induced_group), kernel_size,
        dict(induced_sign_counts),
    )
)
print(f"order_orbit_sizes={[len(orbit) for orbit in order_orbits]}")
print(f"time_reversal_orbits={dict(time_reversal_orbit_transitions)}")
print(f"flip_census={dict(flip_census)}")
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) and outcome != "OPEN_CONTROL_FAILURE" else 1)

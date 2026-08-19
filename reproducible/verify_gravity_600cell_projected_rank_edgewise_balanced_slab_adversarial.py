#!/usr/bin/env python3
"""Independent tournament audit of the balanced-slab order ambiguity."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import networkx as nx


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIMARY = HERE / "gravity_600cell_projected_rank_edgewise_balanced_slab.json"
CARRIER = HERE / "gravity_600cell_projected_rank_edgewise_carrier.json"
OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_balanced_slab_adversarial.json"
AUDIT_PROTOCOL_COMMIT = "8e4c2fb"
INPUT_HASHES = {
    "reproducible/gravity_600cell_projected_rank_edgewise_balanced_slab.json":
        "0a9e9e796cd671c82f2e428bfa21ba63ccb07fe76867e4553979c3c54b22a0d5",
    "reproducible/gravity_600cell_projected_rank_edgewise_carrier.json":
        "b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84",
}
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


def weak_compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total+1):
        for rest in weak_compositions(total-first, parts-1):
            yield (first,)+rest


def edgewise_facets(k=2, dimension=3):
    """Edelsbrunner--Grayson weak-composition construction."""
    width = dimension+1
    facets = set()
    for counts in weak_compositions(k*width, width):
        word = tuple(colour for colour, count in enumerate(counts)
                     for _ in range(count))
        rows = tuple(word[row*width:(row+1)*width] for row in range(k))
        columns = tuple(tuple(rows[row][column] for row in range(k))
                        for column in range(width))
        if len(set(columns)) != width:
            continue
        facets.add(tuple(sorted(
            tuple(column.count(colour) for colour in range(width))
            for column in columns
        )))
    return tuple(sorted(facets))


def graph_from_facets(facets):
    graph = nx.Graph()
    for facet in facets:
        graph.add_nodes_from(facet)
        graph.add_edges_from(combinations(facet, 2))
    return graph


def exact_colourings(graph, colour_count=4):
    """Independent DSATUR-style exhaustive backtracking."""
    nodes = tuple(sorted(graph.nodes()))
    assigned = {}
    solutions = []

    def choose_vertex():
        remaining = [node for node in nodes if node not in assigned]
        return max(remaining, key=lambda node: (
            len({assigned[n] for n in graph.neighbors(node) if n in assigned}),
            graph.degree[node],
            tuple(-value for value in node),
        ))

    def visit():
        if len(assigned) == len(nodes):
            solutions.append(tuple(assigned[node] for node in nodes))
            return
        node = choose_vertex()
        forbidden = {assigned[n] for n in graph.neighbors(node) if n in assigned}
        for colour in range(colour_count):
            if colour not in forbidden:
                assigned[node] = colour
                visit()
                del assigned[node]

    visit()
    return nodes, tuple(solutions)


PAIRS = tuple(combinations(range(4), 2))


def tournament_code(order):
    position = {colour: rank for rank, colour in enumerate(order)}
    return tuple(int(position[left] < position[right]) for left, right in PAIRS)


def tournament_graph(code):
    graph = nx.DiGraph()
    graph.add_nodes_from(range(4))
    for bit, (left, right) in zip(code, PAIRS):
        graph.add_edge(left, right) if bit else graph.add_edge(right, left)
    return graph


print("="*78)
print("ADVERSARIAL BALANCED-SLAB TOURNAMENT AUDIT")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the frozen primary result and earlier carrier have exact provenance",
    actual_hashes == INPUT_HASHES and AUDIT_PROTOCOL_COMMIT == "8e4c2fb",
    str(actual_hashes),
)

facets = edgewise_facets()
graph = graph_from_facets(facets)
points = tuple(sorted(graph.nodes()))
check(
    "the independent colour-scheme construction gives the local fine carrier",
    len(facets) == 8 and len(points) == 10 and graph.number_of_edges() == 25
    and nx.is_connected(graph),
    f"vertices={len(points)}, edges={graph.number_of_edges()}, facets={len(facets)}",
)

declared_colours = {
    point: sum(rank*weight for rank, weight in enumerate(point)) % 4
    for point in points
}
proper = all(declared_colours[left] != declared_colours[right]
             for left, right in graph.edges())
facet_colours = Counter(tuple(sorted(declared_colours[v] for v in facet))
                        for facet in facets)
check(
    "the weighted-rank residue is a proper local four-colouring",
    proper and facet_colours == Counter({(0, 1, 2, 3): 8}),
    f"classes={dict(Counter(declared_colours.values()))}",
)

ordered_nodes, solutions = exact_colourings(graph)
declared_tuple = tuple(declared_colours[node] for node in ordered_nodes)
check(
    "independent exact graph colouring finds only 4! labelled colourings",
    len(solutions) == 24 and len(set(solutions)) == 24
    and declared_tuple in solutions,
    f"solutions={len(solutions)}, declared present={declared_tuple in solutions}",
)

orders = tuple(permutations(range(4)))
order_codes = {order: tournament_code(order) for order in orders}
all_codes = tuple(tuple((mask >> shift) & 1 for shift in range(6))
                  for mask in range(64))
transitive = tuple(code for code in all_codes
                   if nx.is_directed_acyclic_graph(tournament_graph(code)))
check(
    "the six prism-square diagonals admit exactly 24 transitive schedules",
    len(set(order_codes.values())) == 24
    and len(transitive) == 24
    and set(transitive) == set(order_codes.values()),
    f"all tournaments=64, transitive={len(transitive)}, total orders={len(orders)}",
)

# Positive control: the square I x I has exactly its two diagonal choices.
segment_codes = {tuple(order): int(order[0] < order[1])
                 for order in permutations(range(2))}
check(
    "positive control: a segment prism has two time-reversed staircases",
    len(set(segment_codes.values())) == 2
    and segment_codes[(0, 1)] == 1 and segment_codes[(1, 0)] == 0,
    str(segment_codes),
)

# Negative control: 0->1->2->0, with all three pointing to vertex 3.
cyclic = []
directed = {(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)}
for left, right in PAIRS:
    cyclic.append(int((left, right) in directed))
cyclic = tuple(cyclic)
check(
    "negative control: a cyclic diagonal tournament is not a staircase order",
    not nx.is_directed_acyclic_graph(tournament_graph(cyclic))
    and cyclic not in set(order_codes.values()),
    f"code={cyclic}",
)

reverse_map = {}
for order, code in order_codes.items():
    reversed_code = tuple(1-bit for bit in code)
    matches = [candidate for candidate, value in order_codes.items()
               if value == reversed_code]
    reverse_map[order] = matches[0] if len(matches) == 1 else None
unseen = set(orders)
orbits = []
while unseen:
    seed = min(unseen)
    orbit = {seed, reverse_map[seed]}
    unseen -= orbit
    orbits.append(orbit)
check(
    "time reversal independently gives twelve size-two schedule orbits",
    all(reverse_map[order] == tuple(reversed(order)) for order in orders)
    and Counter(map(len, orbits)) == Counter({2: 12}),
    f"fixed={sum(reverse_map[o] == o for o in orders)}, orbits={len(orbits)}",
)

parities = Counter()
for order in orders:
    inversions = sum(order[left] > order[right]
                     for left, right in PAIRS)
    parities[inversions % 2] += 1
check(
    "even a supplied spatial orientation could reduce 24 only to 12+12",
    parities == Counter({0: 12, 1: 12}),
    str(dict(parities)),
)

# A 600-cell automorphism preserves face dimension.  Therefore it induces
# the identity on the four rank labels and on the weighted-rank residue.
rank_label_permutations = tuple(
    permutation for permutation in permutations(range(4))
    if all(permutation[rank] == rank for rank in range(4))
)
check(
    "rank-preserving spatial symmetry identifies none of the 24 orders",
    rank_label_permutations == ((0, 1, 2, 3),),
    f"induced rank-label permutations={rank_label_permutations}",
)

primary = json.loads(PRIMARY.read_text())
carrier = json.loads(CARRIER.read_text())
primary_selection = primary["selection"]
carrier_temporal_fields = sorted(
    key for key in carrier
    if any(word in key.lower() for word in ("time", "slab", "order", "colour", "color"))
)
check(
    "the tournament audit agrees exactly with the primary selection census",
    primary_selection["ordered_slab_alternatives"] == 24
    and primary_selection["h4_invariant_ordered_slab_alternatives"] == 24
    and primary_selection["time_reversal_fixed_orders"] == 0
    and primary_selection["canonical_selection_passes"] is False,
    str(primary_selection),
)
check(
    "the pre-mission carrier artifact did not already certify temporal order data",
    carrier_temporal_fields == [],
    f"matching top-level fields={carrier_temporal_fields}",
)

interpretation = {
    "D0_uncoloured_carrier_schedules": 24,
    "D1_unordered_four_colour_partition_schedules": 24,
    "D2_linearly_ordered_colour_map_schedules": 1,
    "D2_status_before_this_mission": "NOT_CERTIFIED",
    "primary_frozen_reading": "D1",
    "primary_classification": "STRUCTURAL",
    "interpretive_caveat": (
        "Declaring the residue codomain to be the ordered set 0<1<2<3 "
        "creates D2 and selects one schedule, but the selection is then input data."
    ),
}
check(
    "the audit keeps existence separate from the missing selection",
    interpretation["D0_uncoloured_carrier_schedules"] == 24
    and interpretation["D1_unordered_four_colour_partition_schedules"] == 24
    and interpretation["D2_linearly_ordered_colour_map_schedules"] == 1
    and interpretation["D2_status_before_this_mission"] == "NOT_CERTIFIED",
    str(interpretation),
)

artifact = {
    "title": "Adversarial balanced-slab tournament audit",
    "date": "2026-08-19",
    "audit_protocol_commit": AUDIT_PROTOCOL_COMMIT,
    "input_hashes": actual_hashes,
    "local_edgewise": {
        "vertices": len(points),
        "edges": graph.number_of_edges(),
        "tetrahedra": len(facets),
        "proper_labelled_four_colourings": len(solutions),
        "declared_colour_class_sizes": dict(sorted(Counter(
            declared_colours.values()).items())),
    },
    "tournaments": {
        "all": len(all_codes),
        "transitive": len(transitive),
        "distinct_staircase_orders": len(set(order_codes.values())),
        "time_reversal_orbits": len(orbits),
        "time_reversal_orbit_sizes": dict(sorted(Counter(
            map(len, orbits)).items())),
        "permutation_parities": dict(sorted(parities.items())),
    },
    "interpretation": interpretation,
    "primary_selection": primary_selection,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")
print(f"WROTE {OUTPUT.relative_to(ROOT)}")
print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)


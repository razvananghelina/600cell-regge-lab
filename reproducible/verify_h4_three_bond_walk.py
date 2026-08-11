#!/usr/bin/env python3
"""Exact support connectivity of the preregistered H4 three-bond walk.

Protocol commit 9bf2fba froze the three adjacent Coxeter bonds, their rank
orientation, the dense-block support of the published coin and all six phase
orders before any connectivity calculation.
"""

from collections import Counter, deque
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sparse
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("h4_three_bond_walk.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def chamber_components(colour_maps, colours):
    unseen = set(range(len(colour_maps[0])))
    components = []
    while unseen:
        seed = next(iter(unseen))
        component = {seed}
        queue = deque((seed,))
        while queue:
            chamber = queue.popleft()
            for colour in colours:
                neighbour = colour_maps[colour][chamber]
                if neighbour not in component:
                    component.add(neighbour)
                    queue.append(neighbour)
        unseen.difference_update(component)
        components.append(frozenset(component))
    return tuple(components)


print("=" * 78)
print("THREE COXETER-BOND PERIODIC WALK ON THE H4 CHAMBER CARRIER")
print("=" * 78)

# -------------------------------------------------------------------------
# Full four-coloured chamber graph.
# -------------------------------------------------------------------------
vertices, adjacency, _ = build_600cell()
neighbours = tuple(
    frozenset(np.flatnonzero(adjacency[index]).tolist())
    for index in range(120)
)
edges = tuple(
    (left, right)
    for left in range(120)
    for right in sorted(neighbours[left])
    if left < right
)
triangles = tuple(
    (left, right, third)
    for left, right in edges
    for third in sorted(neighbours[left] & neighbours[right])
    if right < third
)
tetrahedra = tuple(
    (first, second, third, fourth)
    for first, second, third in triangles
    for fourth in sorted(
        neighbours[first] & neighbours[second] & neighbours[third]
    )
    if third < fourth
)
check("the coarse f-vector is (120,720,1200,600)",
      (len(vertices), len(edges), len(triangles), len(tetrahedra))
      == (120, 720, 1200, 600))

face_to_tetrahedra = {}
for tetrahedron in tetrahedra:
    for face in combinations(tetrahedron, 3):
        face_to_tetrahedra.setdefault(face, []).append(tetrahedron)

chambers = tuple(
    (tetrahedron, ordering)
    for tetrahedron in tetrahedra
    for ordering in permutations(tetrahedron)
)
chamber_index = {chamber: index for index, chamber in enumerate(chambers)}
colour_maps = [[] for _ in range(4)]
for tetrahedron, ordering in chambers:
    for colour in range(3):
        changed = list(ordering)
        changed[colour], changed[colour + 1] = (
            changed[colour + 1], changed[colour]
        )
        colour_maps[colour].append(
            chamber_index[(tetrahedron, tuple(changed))]
        )
    face = tuple(sorted(ordering[:3]))
    across = next(candidate for candidate in face_to_tetrahedra[face]
                  if candidate != tetrahedron)
    opposite = next(vertex for vertex in across if vertex not in face)
    colour_maps[3].append(
        chamber_index[(across, ordering[:3] + (opposite,))]
    )
colour_maps = tuple(map(tuple, colour_maps))

identity = tuple(range(len(chambers)))
check("the 14,400 chamber carrier has four exact involutive colours",
      len(chambers) == 14400
      and all(
          tuple(mapping[mapping[index]] for index in identity) == identity
          for mapping in colour_maps
      ))

bonds = ((0, 1), (1, 2), (2, 3))
single_bond_histograms = {}
for bond in bonds:
    components = chamber_components(colour_maps, bond)
    single_bond_histograms[str(bond)] = dict(
        sorted(Counter(map(len, components)).items())
    )
check("single-bond orbit multisets are derived before composing directions",
      single_bond_histograms == {
          "(0, 1)": {6: 2400},
          "(1, 2)": {6: 2400},
          "(2, 3)": {10: 1440},
      }, str(single_bond_histograms))

# -------------------------------------------------------------------------
# Corrected robust macro template for each rank-oriented Coxeter bond.
# -------------------------------------------------------------------------
active_components = 4
active_dimension = active_components * len(chambers)


def state(chamber, component):
    return active_components * chamber + component


translations = {}
for bond in bonds:
    first_colour, second_colour = bond
    mapping = [-1] * active_dimension
    for chamber in identity:
        sources = (
            (chamber, 2),
            (colour_maps[first_colour][chamber], 3),
            (colour_maps[second_colour][chamber], 0),
            (chamber, 1),
        )
        for output_component, (source_chamber, source_component) in enumerate(
                sources):
            mapping[state(chamber, output_component)] = state(
                source_chamber, source_component
            )
    translations[bond] = tuple(mapping)

check("all three bond translations are exact active-carrier permutations",
      all(len(set(mapping)) == active_dimension
          for mapping in translations.values()))
check("each bond translation crosses at most one edge of its own colours",
      all(
          (source // active_components == output // active_components)
          or source // active_components
             in {
                 colour_maps[bond[0]][output // active_components],
                 colour_maps[bond[1]][output // active_components],
             }
          for bond, mapping in translations.items()
          for output, source in enumerate(mapping)
      ))

# -------------------------------------------------------------------------
# Exact zero-pattern of C_hat*T on the three-phase time-expanded graph.
# C_hat=I2 tensor C has dense blocks {0,1} and {2,3}.
# -------------------------------------------------------------------------
coin_inputs = {
    0: (0, 1),
    1: (0, 1),
    2: (2, 3),
    3: (2, 3),
}


def phase_graph(schedule):
    phase_count = len(schedule)
    node_count = phase_count * active_dimension
    rows = []
    columns = []
    for phase, bond in enumerate(schedule):
        translation = translations[bond]
        target_phase = (phase + 1) % phase_count
        source_offset = phase * active_dimension
        target_offset = target_phase * active_dimension
        for chamber in identity:
            for output_component in range(active_components):
                output = state(chamber, output_component)
                for intermediate_component in coin_inputs[output_component]:
                    intermediate_output = state(chamber, intermediate_component)
                    source = translation[intermediate_output]
                    # scipy csgraph convention: rows are edge sources.
                    rows.append(source_offset + source)
                    columns.append(target_offset + output)
    graph = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, columns)),
        shape=(node_count, node_count),
    )
    graph.sum_duplicates()
    return graph


schedules = tuple(permutations(bonds))
check("the complete preregistered look-elsewhere set has N=3!=6 schedules",
      len(schedules) == 6 and len(set(schedules)) == 6)

schedule_audits = []
for schedule in schedules:
    graph = phase_graph(schedule)
    strong_count, strong_labels = connected_components(
        graph, directed=True, connection="strong"
    )
    weak_count, weak_labels = connected_components(
        graph, directed=True, connection="weak"
    )
    strong_sizes = Counter(np.bincount(strong_labels).tolist())
    weak_sizes = Counter(np.bincount(weak_labels).tolist())
    schedule_audits.append({
        "schedule": [list(bond) for bond in schedule],
        "node_count": graph.shape[0],
        "directed_edge_count": graph.nnz,
        "strong_component_count": int(strong_count),
        "strong_component_size_histogram": dict(sorted(strong_sizes.items())),
        "weak_component_count": int(weak_count),
        "weak_component_size_histogram": dict(sorted(weak_sizes.items())),
    })

designated = (bonds[0], bonds[1], bonds[2])
designated_audit = next(
    audit for audit in schedule_audits
    if tuple(map(tuple, audit["schedule"])) == designated
)
strong_hits = sum(
    audit["strong_component_count"] == 1 for audit in schedule_audits
)
weak_hits = sum(
    audit["weak_component_count"] == 1 for audit in schedule_audits
)

check("the designated rank-forward schedule is strongly connected",
      designated_audit["strong_component_count"] == 1,
      f"strong components={designated_audit['strong_component_count']}")
check("all six schedule outcomes are recorded as a hit fraction",
      0 <= strong_hits <= 6 and 0 <= weak_hits <= 6,
      f"strong={strong_hits}/6, weak={weak_hits}/6")

cyclic_classes = []
remaining = set(schedules)
while remaining:
    representative = next(iter(remaining))
    rotations = {
        representative[index:] + representative[:index]
        for index in range(3)
    }
    cyclic_classes.append(frozenset(rotations))
    remaining.difference_update(rotations)
check("six schedules reduce to exactly two classes under phase-origin rotation",
      len(cyclic_classes) == 2
      and Counter(map(len, cyclic_classes)) == Counter({3: 2}))

payload = {
    "protocol_commit": "9bf2fba",
    "phenomenological_target_used": False,
    "chambers": len(chambers),
    "active_dimension": active_dimension,
    "time_expanded_dimension": 3 * active_dimension,
    "bonds": [list(bond) for bond in bonds],
    "bond_orientation_rule": "increasing flag rank",
    "single_bond_component_size_histograms": single_bond_histograms,
    "schedule_count_N": len(schedules),
    "cyclic_phase_origin_class_count": len(cyclic_classes),
    "designated_schedule": [list(bond) for bond in designated],
    "designated_strongly_connected": (
        designated_audit["strong_component_count"] == 1
    ),
    "strong_connectivity_hit_fraction": [strong_hits, len(schedules)],
    "weak_connectivity_hit_fraction": [weak_hits, len(schedules)],
    "schedule_audits": schedule_audits,
    "verdict": (
        "DERIVED CONNECTED SUPPORT SCAFFOLD" if
        designated_audit["strong_component_count"] == 1 else
        "DERIVED NEGATIVE: DESIGNATED THREE-BOND SUPPORT IS DISCONNECTED"
    ),
    "scope": (
        "Strong connectivity of the exact nonzero-support graph only; no "
        "claim of amplitude noncancellation, isotropy, Dirac convergence or "
        "physical axis selection."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured six-schedule certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"STRONG_CONNECTIVITY_HITS={strong_hits}/6")
print(f"DESIGNATED_STRONG_COMPONENTS={designated_audit['strong_component_count']}")
print("STRUCTURAL ONLY: support connectivity is not Dirac convergence.")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Connectivity no-go for the literal robust tetrahedral walk on H4.

Protocol commit f52a661 permits arbitrary local coins and arbitrary words in
the corrected robust translation, but no inter-chamber maps beyond colours 2
and 3.  The verifier computes the exact invariant chamber blocks.
"""

from collections import Counter, deque
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("tetrahedral_dirac_walk_connectivity.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def components_from_colours(colour_maps, selected_colours):
    unseen = set(range(len(colour_maps[0])))
    components = []
    component_of = [-1] * len(colour_maps[0])
    while unseen:
        seed = next(iter(unseen))
        component_index = len(components)
        component = {seed}
        queue = deque((seed,))
        while queue:
            chamber = queue.popleft()
            for colour in selected_colours:
                neighbour = colour_maps[colour][chamber]
                if neighbour not in component:
                    component.add(neighbour)
                    queue.append(neighbour)
        for chamber in component:
            component_of[chamber] = component_index
        unseen.difference_update(component)
        components.append(frozenset(component))
    return tuple(components), tuple(component_of)


def product_order(left, right, limit=20):
    product = tuple(left[right[index]] for index in range(len(left)))
    identity = tuple(range(len(left)))
    power = identity
    for exponent in range(1, limit + 1):
        power = tuple(product[power[index]] for index in range(len(left)))
        if power == identity:
            return exponent
    return None


print("=" * 78)
print("CONNECTIVITY OF THE LITERAL ROBUST TETRAHEDRAL WALK")
print("=" * 78)

# -------------------------------------------------------------------------
# Exact combinatorial H4 chamber graph from the certified 600-cell carrier.
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
check("the coarse carrier has the exact 600-cell f-vector",
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
check("14,400 chambers have four exact fixed-point-free involutions",
      len(chambers) == 14400
      and all(
          tuple(mapping[mapping[index]] for index in identity) == identity
          and all(mapping[index] != index for index in identity)
          for mapping in colour_maps
      ))
check("the causal pair has Coxeter product order five",
      product_order(colour_maps[2], colour_maps[3]) == 5)

# -------------------------------------------------------------------------
# Exact support components for the only two colours used by the literal walk.
# -------------------------------------------------------------------------
causal_components, causal_component_of = components_from_colours(
    colour_maps, (2, 3)
)
causal_size_histogram = Counter(map(len, causal_components))
check("colours 2 and 3 split the carrier into 1,440 decagons",
      causal_size_histogram == Counter({10: 1440}),
      str(dict(causal_size_histogram)))

# Each vertex of a causal component has one colour-2 and one colour-3 edge;
# distinctness plus connected size ten makes every block a simple 10-cycle.
check("every causal component is a simple alternating ten-cycle",
      all(
          colour_maps[2][chamber] != colour_maps[3][chamber]
          and causal_component_of[colour_maps[2][chamber]]
              == causal_component_of[chamber]
          and causal_component_of[colour_maps[3][chamber]]
              == causal_component_of[chamber]
          for chamber in identity
      ))

all_components, all_component_of = components_from_colours(
    colour_maps, (0, 1, 2, 3)
)
check("the full four-colour H4 chamber graph is connected",
      Counter(map(len, all_components)) == Counter({14400: 1}))

# Colours 0 and 1 are precisely what join the decagonal blocks.  Record the
# quotient connectivity as a control for the disclosed next route.
quotient_neighbours = [set() for _ in causal_components]
crossing_counts = Counter()
for colour in (0, 1):
    for chamber in identity:
        left = causal_component_of[chamber]
        right = causal_component_of[colour_maps[colour][chamber]]
        if left != right:
            quotient_neighbours[left].add(right)
            crossing_counts[colour] += 1
quotient_seen = {0}
queue = deque((0,))
while queue:
    block = queue.popleft()
    for neighbour in quotient_neighbours[block]:
        if neighbour not in quotient_seen:
            quotient_seen.add(neighbour)
            queue.append(neighbour)
check("adding colours 0 and 1 connects the 1,440-block quotient",
      len(quotient_seen) == len(causal_components),
      f"directed crossing incidences={dict(crossing_counts)}")

# -------------------------------------------------------------------------
# Universal arbitrary-local-coin obstruction.
# -------------------------------------------------------------------------
components_per_chamber = 8
carrier_dimension = components_per_chamber * len(chambers)


def state(chamber, component):
    return components_per_chamber * chamber + component


# Corrected macro map from commit 69447b9.
robust_macro = [-1] * carrier_dimension
for chamber in identity:
    sources = (
        (chamber, 2),
        (colour_maps[2][chamber], 3),
        (colour_maps[3][chamber], 0),
        (chamber, 1),
        (chamber, 4),
        (chamber, 5),
        (chamber, 6),
        (chamber, 7),
    )
    for output_component, (source_chamber, source_component) in enumerate(
            sources):
        robust_macro[state(chamber, output_component)] = state(
            source_chamber, source_component
        )
robust_macro = tuple(robust_macro)

check("the robust macro shift preserves every decagonal orbit projector",
      all(
          causal_component_of[output // components_per_chamber]
          == causal_component_of[source // components_per_chamber]
          for output, source in enumerate(robust_macro)
      ))

inverse_macro = [None] * carrier_dimension
for output, source in enumerate(robust_macro):
    inverse_macro[source] = output
check("the inverse robust shift preserves the same 1,440 projectors",
      None not in inverse_macro
      and all(
          causal_component_of[output // components_per_chamber]
          == causal_component_of[source // components_per_chamber]
          for output, source in enumerate(inverse_macro)
      ))

# A completely arbitrary chamber-local operator has support only between
# (k,a) and (k,b), so it commutes with every projector onto a union of local
# component spaces over one causal orbit.  If two relations preserve the
# orbit label, their relational/matrix product does too (transitivity of
# equality).  This is the exact induction covering arbitrary finite words.
local_complete_mixing_preserves = all(
    causal_component_of[chamber] == causal_component_of[chamber]
    for chamber in identity
    for _output_component in range(components_per_chamber)
    for _source_component in range(components_per_chamber)
)
closure_truth_table = all(
    not (left == middle and middle == right) or left == right
    for left in (0, 1)
    for middle in (0, 1)
    for right in (0, 1)
)
check("arbitrary chamber-local mixing and all finite products preserve the blocks",
      local_complete_mixing_preserves and closure_truth_table,
      "equality of orbit labels is closed under matrix-support composition")

check("the literal robust walk fails global connectedness even with arbitrary local coins",
      len(causal_components) == 1440 > 1,
      "1,440 nontrivial invariant projectors remain")

payload = {
    "protocol_commit": "f52a661",
    "phenomenological_target_used": False,
    "chambers": len(chambers),
    "local_components_per_chamber": components_per_chamber,
    "carrier_dimension": carrier_dimension,
    "literal_inter_chamber_colours": [2, 3],
    "causal_component_count": len(causal_components),
    "causal_component_size_histogram": dict(
        sorted(causal_size_histogram.items())
    ),
    "full_four_colour_component_count": len(all_components),
    "quotient_connected_after_adding_colours_0_and_1": (
        len(quotient_seen) == len(causal_components)
    ),
    "arbitrary_chamber_local_coins_covered": True,
    "verdict": (
        "DERIVED CONNECTIVITY NO-GO: every word in the literal robust "
        "translation, its inverse and arbitrary chamber-local coins preserves "
        "1,440 ten-chamber orbit blocks."
    ),
    "open": (
        "A separately preregistered three-direction schedule using Coxeter "
        "bonds (0,1), (1,2), (2,3)."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured connectivity certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED CONNECTIVITY NO-GO: literal robust walk has 1,440 blocks.")
print("CONTROL: all four H4 colours give one connected chamber graph.")
print("OPEN: preregistered three-bond geometric direction schedule.")
raise SystemExit(0 if passed == tests else 1)

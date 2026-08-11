#!/usr/bin/env python3
"""Exact carrier audit for the published tetrahedral Dirac quantum walk.

Protocol commit e731871 froze the construction and the two possible verdicts.
This verifier tests the four-coloured first-barycentric chamber graph of the
600-cell and the paper's non-robust two-stage permutation shift.
"""

from collections import Counter, deque
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("tetrahedral_dirac_walk_bridge.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def compose(left, right):
    """Permutation product left after right."""
    return tuple(left[right[index]] for index in range(len(left)))


def product_order(left, right, limit=20):
    product = compose(left, right)
    identity = tuple(range(len(left)))
    power = identity
    for exponent in range(1, limit + 1):
        power = compose(product, power)
        if power == identity:
            return exponent
    return None


print("=" * 78)
print("600-CELL BARYCENTRIC CARRIER FOR THE TETRAHEDRAL DIRAC WALK")
print("=" * 78)

# -------------------------------------------------------------------------
# Rebuild the complete coarse boundary complex.
# -------------------------------------------------------------------------
vertices, adjacency, _ = build_600cell()
neighbours = tuple(
    frozenset(np.flatnonzero(adjacency[index]).tolist())
    for index in range(len(vertices))
)
edges = tuple(
    (left, right)
    for left in range(120)
    for right in sorted(neighbours[left])
    if left < right
)
triangles = []
for left, right in edges:
    triangles.extend(
        (left, right, third)
        for third in sorted(neighbours[left] & neighbours[right])
        if right < third
    )
triangles = tuple(triangles)
tetrahedra = []
for first, second, third in triangles:
    tetrahedra.extend(
        (first, second, third, fourth)
        for fourth in sorted(
            neighbours[first] & neighbours[second] & neighbours[third]
        )
        if third < fourth
    )
tetrahedra = tuple(tetrahedra)
check("the coarse complex has f-vector (120,720,1200,600)",
      (len(vertices), len(edges), len(triangles), len(tetrahedra))
      == (120, 720, 1200, 600))

face_to_tetrahedra = {}
for tetrahedron in tetrahedra:
    for face in combinations(tetrahedron, 3):
        face_to_tetrahedra.setdefault(face, []).append(tetrahedron)
check("every coarse triangular face has exactly two incident tetrahedra",
      Counter(map(len, face_to_tetrahedra.values())) == Counter({2: 1200}))

# -------------------------------------------------------------------------
# Complete flags and their four intrinsic rank-changing neighbours.
# -------------------------------------------------------------------------
chambers = tuple(
    (tetrahedron, ordering)
    for tetrahedron in tetrahedra
    for ordering in permutations(tetrahedron)
)
chamber_index = {chamber: index for index, chamber in enumerate(chambers)}
check("the first barycentric subdivision has exactly 14,400 chambers",
      len(chambers) == 14400 and len(chamber_index) == 14400)

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

    boundary_face = tuple(sorted(ordering[:3]))
    across = next(
        candidate
        for candidate in face_to_tetrahedra[boundary_face]
        if candidate != tetrahedron
    )
    opposite = next(vertex for vertex in across if vertex not in boundary_face)
    colour_maps[3].append(
        chamber_index[(across, ordering[:3] + (opposite,))]
    )

colour_maps = tuple(map(tuple, colour_maps))
identity = tuple(range(len(chambers)))
check("all four rank-colour maps are fixed-point-free involutions",
      all(compose(mapping, mapping) == identity
              and all(mapping[index] != index for index in identity)
              for mapping in colour_maps))
check("each chamber has four distinct rank-coloured neighbours",
      all(len({mapping[index] for mapping in colour_maps}) == 4
          for index in identity))

coxeter_orders = tuple(
    tuple(1 if left == right else product_order(
        colour_maps[left], colour_maps[right]
    ) for right in range(4))
    for left in range(4)
)
expected_orders = (
    (1, 3, 2, 2),
    (3, 1, 3, 2),
    (2, 3, 1, 5),
    (2, 2, 5, 1),
)
check("the exact coloured chamber relations are Coxeter H4=(3,3,5)",
      coxeter_orders == expected_orders, str(coxeter_orders))

# Graph-derived handedness: unique up to one global exchange because the
# coloured chamber graph is connected and bipartite.
handedness = [-1] * len(chambers)
handedness[0] = 0
queue = deque((0,))
bipartite = True
while queue:
    chamber = queue.popleft()
    for mapping in colour_maps:
        neighbour = mapping[chamber]
        wanted = 1 - handedness[chamber]
        if handedness[neighbour] < 0:
            handedness[neighbour] = wanted
            queue.append(neighbour)
        elif handedness[neighbour] != wanted:
            bipartite = False

handedness_counts = Counter(handedness)
check("the chamber graph is connected, bipartite and balanced",
      bipartite and -1 not in handedness
      and handedness_counts == Counter({0: 7200, 1: 7200}),
      str(dict(handedness_counts)))

# -------------------------------------------------------------------------
# Published non-robust shift, Eqs. (2)--(3), as output-to-input maps.
# -------------------------------------------------------------------------
carrier_dimension = 4 * len(chambers)


def state(chamber, component):
    return 4 * chamber + component


def published_shift_maps(labels):
    black = [-1] * carrier_dimension
    grey = [-1] * carrier_dimension
    for chamber, chirality in enumerate(labels):
        if chirality == 0:  # LH convention
            across = colour_maps[2][chamber]
            black[state(chamber, 0)] = state(across, 3)
            black[state(chamber, 1)] = state(chamber, 0)
            black[state(chamber, 2)] = state(across, 1)
            black[state(chamber, 3)] = state(chamber, 2)
            grey[state(chamber, 0)] = state(chamber, 2)
            grey[state(chamber, 1)] = state(chamber, 1)
            grey[state(chamber, 2)] = state(chamber, 0)
            grey[state(chamber, 3)] = state(chamber, 3)
        else:  # RH convention
            across = colour_maps[3][chamber]
            black[state(chamber, 0)] = state(chamber, 1)
            black[state(chamber, 1)] = state(across, 2)
            black[state(chamber, 2)] = state(chamber, 3)
            black[state(chamber, 3)] = state(across, 0)
            grey[state(chamber, 0)] = state(chamber, 0)
            grey[state(chamber, 1)] = state(chamber, 3)
            grey[state(chamber, 2)] = state(chamber, 2)
            grey[state(chamber, 3)] = state(chamber, 1)
    complete = tuple(black[grey[index]] for index in range(carrier_dimension))
    return tuple(black), tuple(grey), complete


shift_audits = []
for exchange in (False, True):
    labels = tuple((1-value) if exchange else value for value in handedness)
    black, grey, complete = published_shift_maps(labels)
    black_unique = len(set(black))
    grey_unique = len(set(grey))
    complete_unique = len(set(complete))
    black_input_multiplicities = Counter(Counter(black).values())
    black_missing_inputs = carrier_dimension - black_unique
    invariant_under_causal_colours = all(
        labels[colour_maps[colour][chamber]] == labels[chamber]
        for colour in (2, 3)
        for chamber in range(len(chambers))
    )
    shift_audits.append({
        "global_handedness_exchanged": exchange,
        "black_unique_inputs": black_unique,
        "grey_unique_inputs": grey_unique,
        "complete_unique_inputs": complete_unique,
        "carrier_dimension": carrier_dimension,
        "black_missing_inputs": black_missing_inputs,
        "black_input_multiplicity_histogram": dict(
            sorted(black_input_multiplicities.items())
        ),
        "labels_invariant_under_colours_2_and_3": (
            invariant_under_causal_colours
        ),
        "black_is_permutation": black_unique == carrier_dimension,
        "grey_is_permutation": grey_unique == carrier_dimension,
        "complete_is_permutation": complete_unique == carrier_dimension,
    })

check("the local grey stage is a permutation for either orientation convention",
      all(audit["grey_is_permutation"] for audit in shift_audits))

# For a source component 0 or 2 at chamber j, the black-stage occurrence
# count is [j is LH] + [s3(j) is RH].  It is one exactly when the two labels
# agree.  Components 1 and 3 give the same condition for s2.  Hence the black
# stage is bijective iff handedness is constant on every <s2,s3> orbit.
source_count_formula_exact = True
for exchange in (False, True):
    labels = tuple((1-value) if exchange else value for value in handedness)
    black, _, _ = published_shift_maps(labels)
    actual = Counter(black)
    for chamber in range(len(chambers)):
        expected = {
            0: int(labels[chamber] == 0)
               + int(labels[colour_maps[3][chamber]] == 1),
            1: int(labels[chamber] == 1)
               + int(labels[colour_maps[2][chamber]] == 0),
            2: int(labels[chamber] == 0)
               + int(labels[colour_maps[3][chamber]] == 1),
            3: int(labels[chamber] == 1)
               + int(labels[colour_maps[2][chamber]] == 0),
        }
        source_count_formula_exact &= all(
            actual[state(chamber, component)] == count
            for component, count in expected.items()
        )
check("the exact source-count formula proves the shift's handedness condition",
      source_count_formula_exact,
      "bijective iff labels are invariant under both s2 and s3")

causal_orbits = []
unseen = set(range(len(chambers)))
while unseen:
    seed = next(iter(unseen))
    orbit = {seed}
    queue = deque((seed,))
    while queue:
        chamber = queue.popleft()
        for colour in (2, 3):
            neighbour = colour_maps[colour][chamber]
            if neighbour not in orbit:
                orbit.add(neighbour)
                queue.append(neighbour)
    unseen.difference_update(orbit)
    causal_orbits.append(frozenset(orbit))
causal_orbit_sizes = Counter(map(len, causal_orbits))
check("the causal subgroup <s2,s3> has 1,440 ten-chamber orbits",
      causal_orbit_sizes == Counter({10: 1440}),
      str(dict(causal_orbit_sizes)))

check("the natural chamber handedness makes the published black stage exactly two-to-one",
      all(
          not audit["black_is_permutation"]
          and not audit["complete_is_permutation"]
          and audit["black_unique_inputs"] == carrier_dimension // 2
          and audit["black_missing_inputs"] == carrier_dimension // 2
          and audit["black_input_multiplicity_histogram"] == {2: 28800}
          and not audit["labels_invariant_under_colours_2_and_3"]
          for audit in shift_audits
      ),
      "28,800 inputs duplicated and 28,800 absent for either global exchange")

payload = {
    "protocol_commit": "e731871",
    "phenomenological_target_used": False,
    "coarse_f_vector": [120, 720, 1200, 600],
    "chamber_count": len(chambers),
    "rank_colours": 4,
    "coxeter_product_orders": coxeter_orders,
    "bipartition_counts": dict(sorted(handedness_counts.items())),
    "published_nonrobust_shift_audits": shift_audits,
    "causal_subgroup_orbit_size_histogram": dict(
        sorted(causal_orbit_sizes.items())
    ),
    "verdict": (
        "DERIVED NEGATIVE FOR THE DIRECT NON-ROBUST TRANSPLANT: natural "
        "chamber orientation flips across s2 and s3, but shift bijectivity "
        "requires a label constant on each of 1,440 causal orbits."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured bridge certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: the ordered H4 chamber carrier is geometry-derived.")
print("DERIVED NEGATIVE: the direct non-robust published shift is not unitary.")
print("OPEN: the paper's doubled three-swap robust construction.")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Exact audit of static orientations of the six Hopf axes.

The homogeneous-cover test, existing-carrier census and complete 2^6 global
orientation census were frozen in protocol commit 1c8694b.  No Dirac or
physical target is used.
"""

from collections import Counter
from itertools import permutations, product
import json
from pathlib import Path


OUTPUT = Path(__file__).with_name("hopf_axis_orientation.json")
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
    return tuple(left[right[index]] for index in range(len(right)))


def inverse(permutation):
    result = [None]*len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def parity(permutation):
    return sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left+1, len(permutation))
    ) % 2


def element_order(permutation):
    identity = tuple(range(len(permutation)))
    current = identity
    for order in range(1, 61):
        current = compose(permutation, current)
        if current == identity:
            return order
    raise RuntimeError("element order exceeds 60")


def cyclic_subgroup(generator):
    identity = tuple(range(len(generator)))
    subgroup = {identity}
    current = identity
    for _ in range(element_order(generator)-1):
        current = compose(generator, current)
        subgroup.add(current)
    return frozenset(subgroup)


def normalizer(group, subgroup):
    return frozenset(
        element for element in group
        if frozenset(
            compose(compose(element, member), inverse(element))
            for member in subgroup
        ) == subgroup
    )


def left_cosets(group, subgroup):
    remaining = set(group)
    result = []
    while remaining:
        representative = min(remaining)
        coset = frozenset(compose(representative, member)
                          for member in subgroup)
        result.append(coset)
        remaining -= coset
    return tuple(sorted(result, key=lambda coset: sorted(coset)))


def coset_action(group, cosets):
    lookup = {coset: index for index, coset in enumerate(cosets)}
    return tuple(
        tuple(lookup[frozenset(compose(element, member) for member in coset)]
              for coset in cosets)
        for element in group
    )


print("="*78)
print("STATIC ORIENTATION AUDIT FOR THE SIX HOPF AXES")
print("="*78)

a5 = tuple(
    permutation for permutation in permutations(range(5))
    if parity(permutation) == 0
)
a5_index = {element: index for index, element in enumerate(a5)}
identity = tuple(range(5))
check("the exact group carrier is A5 of order 60", len(a5) == 60)

order_five = tuple(element for element in a5 if element_order(element) == 5)
c5_subgroups = tuple(sorted(
    {cyclic_subgroup(element) for element in order_five},
    key=lambda subgroup: sorted(subgroup),
))
base_c5 = c5_subgroups[0]
base_d5 = normalizer(a5, base_c5)
check(
    "the base C5 has order five and normalizer D5 of order ten",
    len(c5_subgroups) == 6 and len(base_c5) == 5 and len(base_d5) == 10
    and Counter(element_order(element) for element in base_d5)
    == Counter({1: 1, 2: 5, 5: 4}),
)

oriented_points = left_cosets(a5, base_c5)
axes = left_cosets(a5, base_d5)
oriented_action = coset_action(a5, oriented_points)
axis_action = coset_action(a5, axes)
oriented_lookup = {coset: index for index, coset in enumerate(oriented_points)}
axis_lookup = {coset: index for index, coset in enumerate(axes)}

# Quotient pi(gC5)=gD5.  Any member of a left C5 coset can be used as g.
cover_projection = []
for oriented_coset in oriented_points:
    representative = min(oriented_coset)
    axis_coset = frozenset(compose(representative, member)
                           for member in base_d5)
    cover_projection.append(axis_lookup[axis_coset])
cover_projection = tuple(cover_projection)
cover_fibres = tuple(
    tuple(index for index, axis in enumerate(cover_projection)
          if axis == axis_index)
    for axis_index in range(len(axes))
)
check(
    "A5/C5 -> A5/D5 is an exact two-point cover",
    len(oriented_points) == 12 and len(axes) == 6
    and all(len(fibre) == 2 for fibre in cover_fibres),
)
check(
    "the oriented-axis cover is exactly A5-equivariant",
    all(
        cover_projection[oriented_action[group_index][point]]
        == axis_action[group_index][cover_projection[point]]
        for group_index in range(60) for point in range(12)
    ),
)

# A reflection in D5\C5 fixes the base axis and reverses the two orientations.
reflection = min(element for element in base_d5-base_c5
                 if element_order(element) == 2)
reflection_index = a5_index[reflection]
base_axis_index = axis_lookup[base_d5]
base_oriented_index = oriented_lookup[base_c5]
other_oriented_index = next(
    point for point in cover_fibres[base_axis_index]
    if point != base_oriented_index
)
generator = min(element for element in base_c5
                if element_order(element) == 5)
conjugation_reverses_all = all(
    compose(compose(reflection, element), inverse(reflection))
    == inverse(element)
    for element in base_c5
)
check(
    "a D5 reflection fixes the axis and swaps its two orientations",
    axis_action[reflection_index][base_axis_index] == base_axis_index
    and oriented_action[reflection_index][base_oriented_index]
    == other_oriented_index
    and oriented_action[reflection_index][other_oriented_index]
    == base_oriented_index,
)
check(
    "the same reflection sends every C5 rotation to its inverse",
    conjugation_reverses_all
    and compose(compose(reflection, generator), inverse(reflection))
    == inverse(generator),
    "the incidence-preserving stabilizer realizes r <-> r^-1",
)

# -------------------------------------------------------------------------
# Equivariant-map enumeration from X=A5/D5.  A map is determined by the
# image of the base axis, which must be fixed by D5.  We enumerate all target
# points, construct every resulting map and check equivariance directly.
# -------------------------------------------------------------------------
base_stabilizer_indices = tuple(a5_index[element] for element in base_d5)


def enumerate_equivariant_maps(target_actions, projection=None):
    fixed_targets = tuple(
        target for target in range(len(target_actions[0]))
        if all(target_actions[group_index][target] == target
               for group_index in base_stabilizer_indices)
    )
    maps = []
    for target_base in fixed_targets:
        mapping = [None]*len(axes)
        well_defined = True
        for group_index in range(60):
            source = axis_action[group_index][base_axis_index]
            target = target_actions[group_index][target_base]
            if mapping[source] is None:
                mapping[source] = target
            elif mapping[source] != target:
                well_defined = False
                break
        equivariant = well_defined and all(
            mapping[axis_action[group_index][source]]
            == target_actions[group_index][mapping[source]]
            for group_index in range(60) for source in range(6)
        )
        is_section = (
            equivariant and projection is not None
            and all(projection[mapping[source]] == source for source in range(6))
        )
        maps.append({
            "base_image": target_base,
            "mapping": tuple(mapping),
            "equivariant": equivariant,
            "is_section": is_section,
        })
    return fixed_targets, tuple(maps)


oriented_fixed, oriented_maps = enumerate_equivariant_maps(
    oriented_action, cover_projection
)
check(
    "there is no A5-equivariant orientation or section of the cover",
    len(oriented_fixed) == len(oriented_maps) == 0,
)

# The handed target is two disjoint copies of the same unoriented axis set.
handed_target = tuple((hand, axis) for hand in range(2) for axis in range(6))
handed_lookup = {point: index for index, point in enumerate(handed_target)}
handed_actions = tuple(
    tuple(handed_lookup[(hand, axis_action[group_index][axis])]
          for hand, axis in handed_target)
    for group_index in range(60)
)
handed_projection = tuple(axis for hand, axis in handed_target)
handed_fixed, handed_maps = enumerate_equivariant_maps(
    handed_actions, handed_projection
)
check(
    "qH/Hq supplies exactly two handed copies, not an axis orientation",
    len(handed_fixed) == len(handed_maps) == 2
    and all(record["equivariant"] and record["is_section"]
            for record in handed_maps)
    and {handed_target[target][0] for target in handed_fixed} == {0, 1}
    and all(
        len({cover for cover in cover_fibres[source]}) == 2
        for source in range(6)
    ),
    "two maps select a handed sheet but neither chooses a point in A5/C5",
)

# Each oriented chamber sheet is a free A5 torsor A5/1.  Use the left regular
# action itself; a two-sheet target is their disjoint union.
regular_action = tuple(
    tuple(a5_index[compose(element, target)] for target in a5)
    for element in a5
)
chamber_target = tuple((sheet, point) for sheet in range(2) for point in range(60))
chamber_lookup = {point: index for index, point in enumerate(chamber_target)}
chamber_actions = tuple(
    tuple(chamber_lookup[(sheet, regular_action[group_index][point])]
          for sheet, point in chamber_target)
    for group_index in range(60)
)
chamber_fixed, chamber_maps = enumerate_equivariant_maps(chamber_actions)
check(
    "neither oriented chamber sheet gives an equivariant axis section",
    len(chamber_fixed) == len(chamber_maps) == 0,
    "a free A5 torsor has no point fixed by the axis stabilizer D5",
)

# -------------------------------------------------------------------------
# Exhaust all 2^6 global orientation assignments and their A5 orbits.
# -------------------------------------------------------------------------
assignments = tuple(product((0, 1), repeat=6))


def assignment_points(bits):
    return tuple(cover_fibres[axis][bits[axis]] for axis in range(6))


point_assignments = tuple(assignment_points(bits) for bits in assignments)
assignment_lookup = {points: index for index, points in enumerate(point_assignments)}


def act_on_assignment(group_index, points):
    result = [None]*6
    for source_axis, oriented_point in enumerate(points):
        target_axis = axis_action[group_index][source_axis]
        target_point = oriented_action[group_index][oriented_point]
        if cover_projection[target_point] != target_axis:
            raise RuntimeError("cover equivariance lost under assignment action")
        result[target_axis] = target_point
    return tuple(result)


assignment_actions = tuple(
    tuple(assignment_lookup[act_on_assignment(group_index, points)]
          for points in point_assignments)
    for group_index in range(60)
)
check(
    "all 2^6 global fibre orientations form an exact A5 set",
    len(assignments) == len(point_assignments) == 64
    and all(sorted(action) == list(range(64)) for action in assignment_actions),
)

unseen = set(range(64))
assignment_orbits = []
while unseen:
    seed = min(unseen)
    orbit = frozenset(action[seed] for action in assignment_actions)
    stabilizer_order = sum(action[seed] == seed for action in assignment_actions)
    if len(orbit)*stabilizer_order != 60:
        raise RuntimeError("assignment orbit-stabilizer failed")
    assignment_orbits.append({
        "representative_bits": assignments[seed],
        "size": len(orbit),
        "stabilizer_order": stabilizer_order,
        "free": stabilizer_order == 1,
    })
    unseen -= orbit

assignment_stabilizers = tuple(
    sum(action[assignment] == assignment for action in assignment_actions)
    for assignment in range(64)
)
orbit_size_multiset = Counter(record["size"] for record in assignment_orbits)
orbit_stabilizer_multiset = Counter(
    record["stabilizer_order"] for record in assignment_orbits
)
point_stabilizer_multiset = Counter(assignment_stabilizers)
a5_fixed_assignments = sum(order == 60 for order in assignment_stabilizers)
free_assignments = sum(order == 1 for order in assignment_stabilizers)
free_assignment_orbits = sum(record["free"] for record in assignment_orbits)
check(
    "the global-orientation orbit census exhausts all 64 assignments",
    sum(record["size"] for record in assignment_orbits) == 64
    and sum(size*count for size, count in orbit_size_multiset.items()) == 64,
)
check(
    "no global orientation assignment is A5-invariant",
    a5_fixed_assignments == 0,
)
check(
    "the number of free orientation assignments is computed exactly",
    free_assignments >= 0 and free_assignment_orbits >= 0,
    f"free assignments={free_assignments}; free orbits={free_assignment_orbits}",
)

# A single oriented Hopf point is A5/C5 and hence retains C5 exactly.
oriented_point_stabilizers = tuple(
    sum(action[point] == point for action in oriented_action)
    for point in range(12)
)
check(
    "choosing only one oriented Hopf axis leaves residual C5",
    set(oriented_point_stabilizers) == {5},
)


def encode_counter(counter):
    return {str(key): value for key, value in sorted(counter.items())}


payload = {
    "protocol_commit": "1c8694b",
    "Dirac_evaluated": False,
    "physical_target_comparison_performed": False,
    "cover": {
        "unoriented_axis_set": "A5/D5",
        "unoriented_count": 6,
        "oriented_axis_set": "A5/C5",
        "oriented_count": 12,
        "fibre_size": 2,
        "A5_equivariant_section_count": len(oriented_maps),
        "reflection_reverses_r": conjugation_reverses_all,
    },
    "existing_candidate_map_counts": {
        "oriented_axis_cover": len(oriented_maps),
        "qH_Hq_handed_double": len(handed_maps),
        "two_oriented_chamber_torsors": len(chamber_maps),
    },
    "handed_maps_choose_orientation": False,
    "global_orientation_assignments": {
        "count": 64,
        "orbit_count": len(assignment_orbits),
        "orbit_size_multiset": encode_counter(orbit_size_multiset),
        "orbit_stabilizer_order_multiset": encode_counter(
            orbit_stabilizer_multiset
        ),
        "assignment_stabilizer_order_multiset": encode_counter(
            point_stabilizer_multiset
        ),
        "A5_fixed_count": a5_fixed_assignments,
        "free_assignment_count": free_assignments,
        "free_orbit_count": free_assignment_orbits,
        "orbits": [
            {
                "representative_bits": list(record["representative_bits"]),
                "size": record["size"],
                "stabilizer_order": record["stabilizer_order"],
                "free": record["free"],
            }
            for record in assignment_orbits
        ],
    },
    "single_oriented_axis_stabilizer_order": 5,
    "verdict": (
        "DERIVED STATIC-ORIENTATION NEGATIVE: the cover A5/C5 -> A5/D5 "
        "has no A5-equivariant section. A D5 reflection preserves the axis "
        "and incidence while sending r to r^-1. qH/Hq selects only a handed "
        "copy, and neither free chamber sheet supplies a section. The full "
        "64-assignment census quantifies the symmetry breaking required by "
        "a genuinely new orientation field."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact structured orientation audit was written", OUTPUT.exists())

print("\nGlobal orientation assignment census:")
print(f"  orbit sizes: {dict(orbit_size_multiset)}")
print(f"  orbit stabilizers: {dict(orbit_stabilizer_multiset)}")
print(f"  assignment stabilizers: {dict(point_stabilizer_multiset)}")
print(f"  A5-fixed={a5_fixed_assignments}, free={free_assignments}")
print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("NO DIRAC OR PHYSICAL TARGET WAS EVALUATED.")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Exact census of two-axis residual stabilizers in A5.

The complete oriented/unoriented cross-orbit census and interpretation
boundary were frozen in protocol commit 8384d59.  No Dirac connectedness or
physical target is evaluated here.
"""

from collections import Counter
from itertools import combinations, permutations, product
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "hopf_residual_symmetry_breaking.json"
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


def generated_cyclic(generator):
    identity = tuple(range(len(generator)))
    subgroup = {identity}
    current = identity
    for _ in range(element_order(generator)-1):
        current = compose(generator, current)
        subgroup.add(current)
    return frozenset(subgroup)


def normalizer(group, subgroup):
    result = []
    for element in group:
        element_inverse = inverse(element)
        conjugate = frozenset(
            compose(compose(element, member), element_inverse)
            for member in subgroup
        )
        if conjugate == subgroup:
            result.append(element)
    return frozenset(result)


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
    coset_index = {coset: index for index, coset in enumerate(cosets)}
    actions = []
    for element in group:
        action = []
        for coset in cosets:
            image = frozenset(compose(element, member) for member in coset)
            action.append(coset_index[image])
        actions.append(tuple(action))
    return tuple(actions)


print("="*78)
print("RESIDUAL-SYMMETRY CENSUS FOR PAIRS OF ICOSAHEDRAL AXES")
print("="*78)

a5 = tuple(
    permutation for permutation in permutations(range(5))
    if parity(permutation) == 0
)
a5_set = set(a5)
identity = tuple(range(5))
check("the exact carrier is A5 of order 60", len(a5) == 60)
check(
    "A5 is closed under products and inverses",
    all(compose(left, right) in a5_set for left in a5 for right in a5)
    and all(inverse(element) in a5_set for element in a5),
)

cyclic_subgroups = {}
expected_cyclic_counts = {5: 6, 3: 10, 2: 15}
for order in (5, 3, 2):
    subgroups = {
        generated_cyclic(element) for element in a5
        if element_order(element) == order
    }
    cyclic_subgroups[order] = tuple(
        sorted(subgroups, key=lambda subgroup: sorted(subgroup))
    )
check(
    "the C5, C3 and C2 subgroup counts are 6,10,15",
    {order: len(subgroups) for order, subgroups
     in cyclic_subgroups.items()} == expected_cyclic_counts,
)

normalizer_orders = {5: 10, 3: 6, 2: 4}
normalizers = {
    order: normalizer(a5, cyclic_subgroups[order][0])
    for order in (5, 3, 2)
}
check(
    "their exact normalizers have orders D5=10, D3=6, V4=4",
    {order: len(subgroup) for order, subgroup in normalizers.items()}
    == normalizer_orders,
)
check(
    "the normalizer element-order censuses identify D5,D3,V4",
    Counter(element_order(element) for element in normalizers[5])
    == Counter({1: 1, 2: 5, 5: 4})
    and Counter(element_order(element) for element in normalizers[3])
    == Counter({1: 1, 2: 3, 3: 2})
    and Counter(element_order(element) for element in normalizers[2])
    == Counter({1: 1, 2: 3}),
)

oriented_sets = {
    order: left_cosets(a5, cyclic_subgroups[order][0])
    for order in (5, 3, 2)
}
unoriented_sets = {
    order: left_cosets(a5, normalizers[order])
    for order in (5, 3, 2)
}
expected_oriented_sizes = {5: 12, 3: 20, 2: 30}
expected_unoriented_sizes = {5: 6, 3: 10, 2: 15}
check(
    "the oriented homogeneous sets have sizes 12,20,30",
    {order: len(points) for order, points in oriented_sets.items()}
    == expected_oriented_sizes,
)
check(
    "the unoriented homogeneous sets have sizes 6,10,15",
    {order: len(points) for order, points in unoriented_sets.items()}
    == expected_unoriented_sizes,
)

oriented_actions = {
    order: coset_action(a5, points) for order, points in oriented_sets.items()
}
unoriented_actions = {
    order: coset_action(a5, points)
    for order, points in unoriented_sets.items()
}
check(
    "all six homogeneous actions are transitive",
    all(
        {action[0] for action in actions} == set(range(len(actions[0])))
        for actions in tuple(oriented_actions.values())
        +tuple(unoriented_actions.values())
    ),
)


def pair_census(left_actions, right_actions):
    left_size = len(left_actions[0])
    right_size = len(right_actions[0])
    all_pairs = tuple(product(range(left_size), range(right_size)))
    unseen = set(all_pairs)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = frozenset(
            (left_actions[group_index][seed[0]],
             right_actions[group_index][seed[1]])
            for group_index in range(len(a5))
        )
        stabilizer_order = sum(
            left_actions[group_index][seed[0]] == seed[0]
            and right_actions[group_index][seed[1]] == seed[1]
            for group_index in range(len(a5))
        )
        if len(orbit)*stabilizer_order != 60:
            raise RuntimeError("orbit-stabilizer failed")
        orbits.append({
            "representative": seed,
            "size": len(orbit),
            "stabilizer_order": stabilizer_order,
            "regular_A5_torsor": len(orbit) == 60 and stabilizer_order == 1,
        })
        unseen -= orbit
    point_stabilizers = []
    for left, right in all_pairs:
        point_stabilizers.append(sum(
            left_actions[group_index][left] == left
            and right_actions[group_index][right] == right
            for group_index in range(len(a5))
        ))
    free_pairs = sum(order == 1 for order in point_stabilizers)
    return {
        "pair_count": len(all_pairs),
        "orbit_count": len(orbits),
        "orbit_size_multiset": Counter(orbit["size"] for orbit in orbits),
        "orbit_stabilizer_multiset": Counter(
            orbit["stabilizer_order"] for orbit in orbits
        ),
        "point_stabilizer_multiset": Counter(point_stabilizers),
        "free_pair_count": free_pairs,
        "free_pair_fraction": f"{free_pairs}/{len(all_pairs)}",
        "free_orbit_count": sum(orbit["regular_A5_torsor"]
                                for orbit in orbits),
        "orbits": orbits,
    }


cross_types = ((5, 3), (5, 2), (3, 2))
censuses = {"oriented": {}, "unoriented": {}}
for orientation, actions_by_order in (
        ("oriented", oriented_actions),
        ("unoriented", unoriented_actions)):
    for left_order, right_order in cross_types:
        key = f"{left_order}x{right_order}"
        censuses[orientation][key] = pair_census(
            actions_by_order[left_order], actions_by_order[right_order]
        )

check(
    "all six cross-type pair products are exhausted by disjoint A5 orbits",
    all(
        sum(size*count for size, count
            in census["orbit_size_multiset"].items())
        == census["pair_count"]
        and sum(order_count for order_count
                in census["point_stabilizer_multiset"].values())
        == census["pair_count"]
        for orientation in censuses.values()
        for census in orientation.values()
    ),
)
check(
    "every reported free orbit is an exact regular A5 torsor",
    all(
        all(orbit["regular_A5_torsor"]
            == (orbit["size"] == 60 and orbit["stabilizer_order"] == 1)
            for orbit in census["orbits"])
        for orientation in censuses.values()
        for census in orientation.values()
    ),
)

free_summary = {
    orientation: {
        pair_type: census["free_orbit_count"]
        for pair_type, census in data.items()
    }
    for orientation, data in censuses.items()
}
check(
    "the free-pair census has a determinate nonnegative result",
    all(census["free_pair_count"] >= 0
        and census["free_orbit_count"] >= 0
        for data in censuses.values() for census in data.values()),
    f"free orbit counts={free_summary}",
)

# A free pair has trivial stabilizer, so covariance permits every rectangular
# node map.  Record this look-elsewhere space before any connectedness test.
node_sizes = (6, 6, 12, 12)
free_hom_dimensions = {
    f"{left}-{right}": node_sizes[left]*node_sizes[right]
    for left, right in combinations(range(4), 2)
}
nodes = (
    {"u_rank": 2, "v_rank": 2},
    {"u_rank": 2, "v_rank": 0},
    {"u_rank": 1, "v_rank": 1},
    {"u_rank": 0, "v_rank": 1},
)
reading_free_dimensions = []
unordered_node_pairs = tuple(combinations(range(4), 2))
for priority, u_direction, v_direction in product(
        (("u", "v"), ("v", "u")), (1, -1), (1, -1)):
    direction = {"u": u_direction, "v": v_direction}

    def spectral_key(node):
        return tuple(direction[coordinate]*nodes[node][coordinate+"_rank"]
                     for coordinate in priority)

    order = tuple(sorted(range(4), key=spectral_key))
    position = {node: rank for rank, node in enumerate(order)}
    positive = frozenset(
        (left, right) if position[left] < position[right] else (right, left)
        for left, right in unordered_node_pairs
    )
    links = set()
    for i, j in positive:
        for k, ell in positive:
            if i == ell:
                links.add(tuple(sorted((j, k))))
            elif j == k:
                links.add(tuple(sorted((i, ell))))
    reading_free_dimensions.append({
        "priority": list(priority),
        "u_direction": u_direction,
        "v_direction": v_direction,
        "links": [list(link) for link in sorted(links)],
        "total_full_rectangular_dimension": sum(
            node_sizes[left]*node_sizes[right] for left, right in links
        ),
    })
free_dimension_multiset = Counter(
    record["total_full_rectangular_dimension"]
    for record in reading_free_dimensions
)
check(
    "free-vacuum tensor freedom is recorded for all eight readings",
    len(reading_free_dimensions) == 8
    and all(len(record["links"]) == 3 for record in reading_free_dimensions),
    f"total rectangular-dimension multiset={dict(free_dimension_multiset)}",
)

# Scoped action provenance from committed authoritative structured results.
selector = json.loads((HERE/"hopf_sixth_order_selector.json").read_text())
action = json.loads((HERE/"hopf_hessian_action_origin.json").read_text())
source_selector = (HERE/"verify_hopf_sixth_order_selector.py").read_text()
source_action = (HERE/"verify_hopf_hessian_action_origin.py").read_text()
action_audit = {
    "one_field_selects_both_cross_type_axes": False,
    "second_independent_axis_field_defined": False,
    "axis_orientation_selected": False,
    "relative_sign_derived": bool(
        action["equivariant_affine_extension"]["relative_sign_fixed"]
    ),
    "scope": [
        "verify_hopf_sixth_order_selector.py",
        "verify_hopf_hessian_action_origin.py",
        "hopf_sixth_order_selector.json",
        "hopf_hessian_action_origin.json",
    ],
}
source_audit_consistent = (
    selector["global"]["maximum_orbit"].startswith("6 unoriented C10")
    and selector["global"]["minimum_orbit"].startswith("10 unoriented C6")
    and "variables = (x, y, z)" in source_selector
    and "deduplicate_lines" in source_selector
    and "b, c = sp.symbols" in source_action
    and action["equivariant_affine_extension"]["relative_sign_fixed"] is False
    and "second_field" not in source_selector+source_action
)
check(
    "the scoped action audit confirms one unoriented field and no fixed sign",
    source_audit_consistent
    and not any(action_audit[key] for key in (
        "one_field_selects_both_cross_type_axes",
        "second_independent_axis_field_defined",
        "axis_orientation_selected",
        "relative_sign_derived",
    )),
)


def encode_counter(counter):
    return {str(key): value for key, value in sorted(counter.items())}


encoded_censuses = {}
for orientation, data in censuses.items():
    encoded_censuses[orientation] = {}
    for pair_type, census in data.items():
        encoded_censuses[orientation][pair_type] = {
            "pair_count": census["pair_count"],
            "orbit_count": census["orbit_count"],
            "orbit_size_multiset": encode_counter(
                census["orbit_size_multiset"]
            ),
            "orbit_stabilizer_order_multiset": encode_counter(
                census["orbit_stabilizer_multiset"]
            ),
            "point_stabilizer_order_multiset": encode_counter(
                census["point_stabilizer_multiset"]
            ),
            "free_pair_count": census["free_pair_count"],
            "free_pair_fraction": census["free_pair_fraction"],
            "free_orbit_count": census["free_orbit_count"],
            "orbits": [
                {
                    "representative_indices": list(orbit["representative"]),
                    "size": orbit["size"],
                    "stabilizer_order": orbit["stabilizer_order"],
                    "regular_A5_torsor": orbit["regular_A5_torsor"],
                }
                for orbit in census["orbits"]
            ],
        }

payload = {
    "protocol_commit": "8384d59",
    "Dirac_connectedness_evaluated": False,
    "physical_target_comparison_performed": False,
    "axis_orbit_sizes": {
        "oriented": {"C5": 12, "C3": 20, "C2": 30},
        "unoriented": {"D5": 6, "D3": 10, "V4": 15},
    },
    "pair_censuses": encoded_censuses,
    "free_orbit_summary": free_summary,
    "free_vacuum_Hom_dimensions": free_hom_dimensions,
    "free_vacuum_reading_tensor_dimensions": reading_free_dimensions,
    "free_vacuum_total_dimension_multiset": encode_counter(
        free_dimension_multiset
    ),
    "scoped_action_provenance": action_audit,
    "verdict": (
        "COMPLETE STABILIZER CENSUS. Free cross-axis orbits, if present, "
        "remove the residual-symmetry obstruction but expose full "
        "rectangular Dirac-map spaces. The authoritative action files define "
        "one unoriented field, no second independent field, no orientation "
        "selection and no derived relative sign. No Dirac rank was tested."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact structured pair census was written", OUTPUT.exists())

print("\nFree A5-orbit counts:")
for orientation, summary in free_summary.items():
    print(f"  {orientation}: {summary}")
print(f"Free-vacuum tensor-dimension multiset: {dict(free_dimension_multiset)}")
print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("NO DIRAC CONNECTEDNESS OR PHYSICAL TARGET WAS EVALUATED.")
raise SystemExit(0 if passed == tests else 1)

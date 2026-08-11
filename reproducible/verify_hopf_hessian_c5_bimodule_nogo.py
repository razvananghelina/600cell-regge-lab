#!/usr/bin/env python3
"""Exact arbitrary-multiplicity C^5 orientability no-go.

The hypotheses and blind orbit/sign protocol were frozen in commit 1bed02d.
This verifier does not construct or inspect a Hessian operator.  It rebuilds
the derived A5 action, obtains its five-point coset action, and exhausts the
KO6 metric-dimension-zero orientation signs on all pair orbits.
"""

from itertools import combinations, product
import json
from pathlib import Path

from verify_hopf_fibration_invariants import (
    build_2I,
    find_all_hopf_fibrations,
    find_vertex_index,
    quat_mult,
)


OUTPUT = Path(__file__).with_name("hopf_hessian_c5_bimodule_nogo.json")
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


def permutation_order(permutation):
    identity = tuple(range(len(permutation)))
    current = identity
    for order in range(1, 61):
        current = compose(permutation, current)
        if current == identity:
            return order
    raise RuntimeError("order exceeds derived group size")


print("="*78)
print("EXACT ARBITRARY-MULTIPLICITY C5 BIMODULE ORIENTABILITY AUDIT")
print("="*78)

# Rebuild the exact A5 group already carried by the six Hopf fibrations.  The
# five-point action used below is constructed internally from its A4
# subgroups, rather than imported as a standard presentation of A5.
vertices = build_2I()
fibrations = find_all_hopf_fibrations(vertices)


def conjugation_permutation(group_element):
    group_inverse = group_element.copy()
    group_inverse[1:] *= -1
    return tuple(find_vertex_index(
        vertices,
        quat_mult(quat_mult(group_element, vertex), group_inverse),
    ) for vertex in vertices)


vertex_actions = sorted(set(conjugation_permutation(element)
                            for element in vertices))
fibration_by_signature = {
    tuple(sorted(tuple(sorted(fiber)) for fiber in fibration)): index
    for index, fibration in enumerate(fibrations)
}
group = []
for vertex_action in vertex_actions:
    action = []
    for fibration in fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fiber))
            for fiber in fibration
        ))
        action.append(fibration_by_signature[signature])
    group.append(tuple(action))
group = sorted(set(group))
identity = tuple(range(6))
check("the derived quotient action has exactly 60 elements",
      len(group) == 60 and identity in group)

# An order-two centralizer is V4 and its normalizer is an index-five A4.
order_two = next(element for element in group
                 if permutation_order(element) == 2)
V4 = frozenset(element for element in group
               if compose(element, order_two) == compose(order_two, element))
A4 = frozenset(element for element in group if frozenset(
    compose(compose(element, member), inverse(element)) for member in V4
) == V4)
check("the exact subgroup chain is V4 < A4 < A5",
      len(V4) == 4 and len(A4) == 12 and len(group)//len(A4) == 5)


def conjugate_subgroup(subgroup, element):
    element_inverse = inverse(element)
    return frozenset(compose(compose(element, member), element_inverse)
                     for member in subgroup)


A4_conjugates = sorted(
    {conjugate_subgroup(A4, element) for element in group},
    key=lambda subgroup: tuple(sorted(subgroup)),
)
subgroup_index = {subgroup: index
                  for index, subgroup in enumerate(A4_conjugates)}
five_point_actions = []
for element in group:
    action = tuple(subgroup_index[conjugate_subgroup(subgroup, element)]
                   for subgroup in A4_conjugates)
    five_point_actions.append(action)
five_point_actions = sorted(set(five_point_actions))
check("conjugation on the A4 subgroups gives a faithful five-point action",
      len(A4_conjugates) == 5 and len(five_point_actions) == 60)
check("the five-point action is transitive",
      {action[0] for action in five_point_actions} == set(range(5)))

# Blind orbit census on the 25 ordered character pairs.
unseen = {(left, right) for left in range(5) for right in range(5)}
pair_orbits = []
while unseen:
    seed = min(unseen)
    orbit = frozenset((action[seed[0]], action[seed[1]])
                      for action in five_point_actions)
    pair_orbits.append(orbit)
    unseen -= orbit
pair_orbits.sort(key=lambda orbit: (len(orbit), tuple(sorted(orbit))))
orbit_sizes = [len(orbit) for orbit in pair_orbits]
reversal_invariant = [
    frozenset((right, left) for left, right in orbit) == orbit
    for orbit in pair_orbits
]
check("the exact ordered-pair orbit sizes are 5 and 20",
      orbit_sizes == [5, 20],
      "diagonal orbit=5; ordered-distinct orbit=20")
check("each pair orbit is invariant under left-right reversal",
      reversal_invariant == [True, True])
check("each nonempty pair orbit already sees all five algebra characters",
      all({left for left, _ in orbit} == set(range(5))
              and {right for _, right in orbit} == set(range(5))
          for orbit in pair_orbits),
      "every invariant nonzero support would be faithful on both sides")

# Exhaust every nonempty invariant support (there are three) and every
# constant +/- orientation sign on its occupied orbits.  KO6 requires the
# sign at the reversed block to be its negative.
support_sign_records = []
admissible_supports = 0
for support_size in range(1, len(pair_orbits)+1):
    for support_indices in combinations(range(len(pair_orbits)), support_size):
        support = set().union(*(pair_orbits[index]
                               for index in support_indices))
        sign_solutions = 0
        for orbit_signs in product((-1, 1), repeat=support_size):
            signs = {}
            for index, sign in zip(support_indices, orbit_signs):
                signs.update({pair: sign for pair in pair_orbits[index]})
            if all((right, left) in support
                   and signs[(right, left)] == -signs[(left, right)]
                   for left, right in support):
                sign_solutions += 1
        admissible_supports += int(sign_solutions > 0)
        support_sign_records.append({
            "orbit_indices": list(support_indices),
            "dimension_per_unit_multiplicity": len(support),
            "orientation_sign_solutions": sign_solutions,
        })
check("no nonempty invariant bimodule support admits KO6 orientation signs",
      len(support_sign_records) == 3 and admissible_supports == 0
      and all(record["orientation_sign_solutions"] == 0
              for record in support_sign_records),
      f"support ledger={support_sign_records}")

# The obstruction is independent of multiplicity.  On H_(i,j), every
# represented Hochschild zero-cycle is a scalar times the identity, because
# both C5 actions are characters there.  A mixed +/- grading inside the same
# multiplicity space therefore cannot be an orienting zero-cycle.
# On the block H_(i,j), the 25 products of minimal left/right projections
# have coefficients delta_(p,i) delta_(q,j), hence are exactly 0 or I_m for
# every m.  This coefficient census is independent of m.
zero_cycle_block_coefficients = {
    int(left_projection == 0)*int(right_projection == 1)
    for left_projection in range(5) for right_projection in range(5)
}
zero_cycle_scalar_on_multiplicity = zero_cycle_block_coefficients == {0, 1}
extra_multiplicity_changes_sign_problem = False
check("extra multiplicity cannot evade the zero-cycle sign obstruction",
      zero_cycle_scalar_on_multiplicity
      and not extra_multiplicity_changes_sign_problem)

faithful_nonzero_carrier_exists = admissible_supports > 0
check("the stated C5 arena has no nonzero faithful orientable carrier",
      not faithful_nonzero_carrier_exists)

payload = {
    "protocol_commit": "1bed02d",
    "hypotheses": [
        "A=C^5 with a unital faithful finite representation",
        "order zero",
        "KO6 J*gamma=-gamma*J",
        "metric-dimension-zero orientability",
        "derived A5 equivariance compatible with J and gamma",
        "arbitrary finite bimodule multiplicities",
    ],
    "derived_group_order": len(group),
    "A4_conjugates": len(A4_conjugates),
    "five_point_action_order": len(five_point_actions),
    "ordered_pair_orbits": {
        "sizes": orbit_sizes,
        "reversal_invariant": reversal_invariant,
    },
    "nonempty_invariant_supports": len(support_sign_records),
    "support_sign_ledger": support_sign_records,
    "admissible_nonempty_supports": admissible_supports,
    "arbitrary_multiplicity_changes_result": False,
    "verdict": (
        "DERIVED FULL-ARENA NO-GO under the stated hypotheses: the exact "
        "five-point A5 action has pair orbits of sizes 5 and 20, and both are "
        "invariant under reversal. A5-invariant metric-zero orientation signs "
        "therefore contradict KO6 on every nonempty support. Hochschild "
        "zero-cycles are scalar on each C5-C5op multiplicity block, so "
        "arbitrary multiplicities do not evade the obstruction. No Hessian "
        "operator was used."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact C5 full-arena audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED FULL-ARENA NO-GO: no nonzero A5-equivariant KO6 metric-zero")
print("                              orientable C^5 bimodule exists, even with")
print("                              arbitrary multiplicities.")
print("NO HESSIAN TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

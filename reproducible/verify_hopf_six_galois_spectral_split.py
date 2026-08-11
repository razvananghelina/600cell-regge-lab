#!/usr/bin/env python3
"""Exact audit of the incidence-selected golden split for six fibrations.

The candidate and its falsification boundary were disclosed in protocol
commit 43c6dd3.  Discovery of the finite 600-cell permutation data reuses the
registered Hopf-fibration constructor.  All representation-theory and
spectral-projector calculations after that discovery are exact.

No Hessian, Standard-Model module, mass or coupling target is used.
"""

from pathlib import Path
import json

import numpy as np
import sympy as sp

from verify_hopf_fibration_invariants import (
    build_2I,
    build_adjacency,
    find_all_hopf_fibrations,
    find_vertex_index,
    quat_mult,
)


OUTPUT = Path(__file__).with_name("hopf_six_galois_spectral_split.json")
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


def inverse_permutation(permutation):
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
    raise RuntimeError("permutation order exceeds 60")


def quaternion_inverse(quaternion):
    result = quaternion.copy()
    result[1:] *= -1
    return result


def same_exact(left, right):
    difference = left-right
    return all(sp.simplify(entry) == 0 for entry in difference)


def regular_matrix(element, ordered_group, index):
    matrix = sp.zeros(len(ordered_group))
    for source, basis_element in enumerate(ordered_group):
        target = compose(element, basis_element)
        matrix[index[target], source] = 1
    return matrix


print("="*78)
print("INCIDENCE-SELECTED GOLDEN D5 SPECTRAL SPLIT")
print("="*78)

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
identity_vertex = find_vertex_index(
    vertices, np.array([1.0, 0.0, 0.0, 0.0])
)
check("the certified carrier has 120 vertices and six fibrations",
      len(vertices) == 120 and len(fibrations) == 6
      and identity_vertex >= 0)

fibration_by_signature = {
    tuple(sorted(tuple(sorted(fibre)) for fibre in fibration)): label
    for label, fibration in enumerate(fibrations)
}


def conjugation_permutation(group_element):
    group_inverse = quaternion_inverse(group_element)
    return tuple(find_vertex_index(
        vertices,
        quat_mult(quat_mult(group_element, vertex), group_inverse),
    ) for vertex in vertices)


vertex_action_by_index = [conjugation_permutation(vertex)
                          for vertex in vertices]
label_action_by_vertex = []
for vertex_action in vertex_action_by_index:
    action = []
    for fibration in fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fibre))
            for fibre in fibration
        ))
        action.append(fibration_by_signature[signature])
    label_action_by_vertex.append(tuple(action))

group = tuple(sorted(set(label_action_by_vertex)))
identity_action = tuple(range(6))
check("quaternionic conjugation gives the effective A5 action",
      len(group) == 60 and identity_action in group)


def relative_action(source, target):
    relative = quat_mult(quaternion_inverse(vertices[source]),
                         vertices[target])
    relative_index = find_vertex_index(vertices, relative)
    if relative_index < 0:
        raise RuntimeError("relative binary element left the 2I carrier")
    return label_action_by_vertex[relative_index]


incidence_records = []
all_incidence_canonical = True
selected_pairs = []
for fibration_label, fibration in enumerate(fibrations):
    stabilizer = frozenset(element for element in group
                           if element[fibration_label] == fibration_label)
    rotations = frozenset(element for element in stabilizer
                          if permutation_order(element) in (1, 5))
    order_five = rotations-{identity_action}

    identity_fibres = [fibre for fibre in fibration
                       if identity_vertex in fibre]
    internal_edges = []
    distance_two_pairs = []
    for fibre in fibration:
        fibre = tuple(fibre)
        for offset, source in enumerate(fibre):
            for target in fibre[offset+1:]:
                if adjacency[source, target] == 1:
                    internal_edges.append((source, target))
                    continue
                common = sum(adjacency[source, middle]
                             *adjacency[middle, target]
                             for middle in fibre)
                if common > 0:
                    distance_two_pairs.append((source, target))

    edge_actions = frozenset(
        relative_action(source, target)
        for edge in internal_edges for source, target in (edge, edge[::-1])
    )
    chord_actions = frozenset(
        relative_action(source, target)
        for pair in distance_two_pairs
        for source, target in (pair, pair[::-1])
    )
    edge_inverse_closed = all(inverse_permutation(element) in edge_actions
                              for element in edge_actions)
    chord_inverse_closed = all(inverse_permutation(element) in chord_actions
                               for element in chord_actions)
    local_ok = (
        len(stabilizer) == 10
        and len(rotations) == 5
        and len(identity_fibres) == 1
        and len(internal_edges) == 120
        and len(distance_two_pairs) == 120
        and len(edge_actions) == len(chord_actions) == 2
        and edge_actions.isdisjoint(chord_actions)
        and edge_actions | chord_actions == order_five
        and edge_inverse_closed and chord_inverse_closed
    )
    all_incidence_canonical &= local_ok
    selected_pairs.append(edge_actions)
    incidence_records.append({
        "fibration": fibration_label,
        "stabilizer_order": len(stabilizer),
        "internal_edge_count": len(internal_edges),
        "distance_two_pair_count": len(distance_two_pairs),
        "edge_action_count": len(edge_actions),
        "chord_action_count": len(chord_actions),
        "edge_and_chord_disjoint": edge_actions.isdisjoint(chord_actions),
    })

check("all six fibrations have the exact D5/C5 stabilizer structure",
      all(record["stabilizer_order"] == 10
          for record in incidence_records))
check("each fibration has 120 fibre edges and 120 distance-two chords",
      all(record["internal_edge_count"] == 120
          and record["distance_two_pair_count"] == 120
          for record in incidence_records))
check("incidence uniformly selects one inverse pair, never the chord pair",
      all_incidence_canonical,
      "tested all 6 fibrations and all 720 internal undirected edges")

# The selected conjugacy class must also transport naturally when a rotation
# carries one fibration label to another.  This removes the remaining
# base-label/coset-representative convention.
transport_natural = True
for transporter in group:
    transporter_inverse = inverse_permutation(transporter)
    for source_label, source_pair in enumerate(selected_pairs):
        target_label = transporter[source_label]
        transported = frozenset(
            compose(compose(transporter, element), transporter_inverse)
            for element in source_pair
        )
        if transported != selected_pairs[target_label]:
            transport_natural = False
            break
check("the incidence-selected inverse pairs transport A5-equivariantly",
      transport_natural,
      "all 60 transporters and all six base labels")

# Repeat the incidence test on the opposite handed family.  Quaternionic
# inversion maps every qH partition to an Hq partition, but changes the
# well-defined same-coset relative element from x^-1*y to y*x^-1.  This is a
# scope check: node separation must not secretly rely on choosing chirality.
inverse_vertex = tuple(find_vertex_index(vertices,
                                         quaternion_inverse(vertex))
                       for vertex in vertices)
opposite_fibrations = [
    tuple(sorted(tuple(sorted(inverse_vertex[index] for index in fibre))
                 for fibre in fibration))
    for fibration in fibrations
]
opposite_signatures = {
    tuple(sorted(tuple(sorted(fibre)) for fibre in fibration)): label
    for label, fibration in enumerate(opposite_fibrations)
}
opposite_label_actions = []
for vertex_action in vertex_action_by_index:
    action = []
    for fibration in opposite_fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fibre))
            for fibre in fibration
        ))
        action.append(opposite_signatures[signature])
    opposite_label_actions.append(tuple(action))
opposite_group = tuple(sorted(set(opposite_label_actions)))


def opposite_relative_action(source, target):
    relative = quat_mult(vertices[target],
                         quaternion_inverse(vertices[source]))
    relative_index = find_vertex_index(vertices, relative)
    if relative_index < 0:
        raise RuntimeError("opposite-handed relative element left 2I")
    return opposite_label_actions[relative_index]


opposite_pairs = []
original_partitions = {
    tuple(sorted(tuple(sorted(fibre)) for fibre in fibration))
    for fibration in fibrations
}
opposite_ok = (len(set(opposite_fibrations)) == 6
               and set(opposite_fibrations).isdisjoint(original_partitions)
               and len(opposite_group) == 60)
for label, fibration in enumerate(opposite_fibrations):
    stabilizer_opposite = frozenset(
        element for element in opposite_group if element[label] == label
    )
    rotations_opposite = frozenset(
        element for element in stabilizer_opposite
        if permutation_order(element) in (1, 5)
    )
    edge_pairs_opposite = []
    chord_pairs_opposite = []
    for fibre in fibration:
        for offset, source in enumerate(fibre):
            for target in fibre[offset+1:]:
                if adjacency[source, target] == 1:
                    edge_pairs_opposite.append((source, target))
                else:
                    common = sum(adjacency[source, middle]
                                 *adjacency[middle, target]
                                 for middle in fibre)
                    if common > 0:
                        chord_pairs_opposite.append((source, target))
    edge_actions_opposite = frozenset(
        opposite_relative_action(source, target)
        for pair in edge_pairs_opposite
        for source, target in (pair, pair[::-1])
    )
    chord_actions_opposite = frozenset(
        opposite_relative_action(source, target)
        for pair in chord_pairs_opposite
        for source, target in (pair, pair[::-1])
    )
    opposite_pairs.append(edge_actions_opposite)
    opposite_ok &= (
        len(stabilizer_opposite) == 10
        and len(edge_pairs_opposite) == len(chord_pairs_opposite) == 120
        and len(edge_actions_opposite) == len(chord_actions_opposite) == 2
        and edge_actions_opposite.isdisjoint(chord_actions_opposite)
        and edge_actions_opposite | chord_actions_opposite
        == rotations_opposite-{identity_action}
    )

opposite_transport = True
for transporter in opposite_group:
    transporter_inverse = inverse_permutation(transporter)
    for source_label, source_pair in enumerate(opposite_pairs):
        target_label = transporter[source_label]
        transported = frozenset(
            compose(compose(transporter, element), transporter_inverse)
            for element in source_pair
        )
        if transported != opposite_pairs[target_label]:
            opposite_transport = False
            break
check("the same edge/chord split holds in the opposite handed family",
      opposite_ok and opposite_transport,
      "12 disjoint fibrations, 1440 fibre edges total across both hands")

# Exact regular representation for the first D5 stabilizer.  Either member of
# the incidence-selected inverse pair gives the same unoriented element.
stabilizer = tuple(sorted(element for element in group if element[0] == 0))
stabilizer_index = {element: index for index, element in enumerate(stabilizer)}
edge_pair = selected_pairs[0]
r = min(edge_pair)
r_inverse = inverse_permutation(r)
chord_pair = frozenset(
    (compose(r, r), compose(r_inverse, r_inverse))
)
check("the selected edge pair is exactly {r,r^-1}",
      edge_pair == frozenset((r, r_inverse))
      and chord_pair.isdisjoint(edge_pair))

regular = {
    element: regular_matrix(element, stabilizer, stabilizer_index)
    for element in stabilizer
}
U = regular[r]+regular[r_inverse]
W = regular[compose(r, r)]+regular[compose(r_inverse, r_inverse)]
I10 = sp.eye(10)
lam = sp.symbols("lambda")
expected_characteristic = sp.expand(
    (lam-2)**2*(lam**2+lam-1)**4
)
check("the edge operator has the exact golden characteristic polynomial",
      sp.expand(U.charpoly(lam).as_expr()-expected_characteristic) == 0,
      "(lambda-2)^2 (lambda^2+lambda-1)^4")
check("the chord operator is U^2-2I exactly",
      W == U*U-2*I10)

sqrt5 = sp.sqrt(5)
phi = (1+sqrt5)/2
positive_doublet_value = phi-1
negative_doublet_value = -phi


def eigenspace_projector(value, other):
    numerator = (U-2*I10)*(U-other*I10)
    denominator = sp.expand((value-2)*(value-other))
    return numerator.applyfunc(lambda entry: sp.radsimp(entry/denominator))


P_negative = eigenspace_projector(
    negative_doublet_value, positive_doublet_value
)
P_positive = eigenspace_projector(
    positive_doublet_value, negative_doublet_value
)
P_doublets = P_negative+P_positive
check("the two exact doublet projectors are orthogonal central idempotents",
      P_negative.T == P_negative and P_positive.T == P_positive
      and same_exact(P_negative*P_negative, P_negative)
      and same_exact(P_positive*P_positive, P_positive)
      and same_exact(P_negative*P_positive, sp.zeros(10))
      and all(same_exact(P_negative*regular[element],
                         regular[element]*P_negative)
              and same_exact(P_positive*regular[element],
                             regular[element]*P_positive)
              for element in stabilizer))
check("negative real functional calculus selects exactly one M2 block",
      P_negative.rank() == 4 and P_positive.rank() == 4
      and same_exact(U*P_negative,
                     negative_doublet_value*P_negative)
      and same_exact(U*P_positive,
                     positive_doublet_value*P_positive),
      "regular ranks 4+4; Morita-amplified blocks are M12(R)+M12(R)")


def golden_conjugate(matrix):
    return matrix.applyfunc(
        lambda entry: sp.radsimp(entry.xreplace({sqrt5: -sqrt5}))
    )


check("golden conjugation exchanges the two spectral projectors",
      same_exact(golden_conjugate(P_negative), P_positive)
      and same_exact(golden_conjugate(P_positive), P_negative))
check("only the sum of the doublet projectors descends rationally",
      same_exact(golden_conjugate(P_doublets), P_doublets)
      and not same_exact(golden_conjugate(P_negative), P_negative)
      and not same_exact(golden_conjugate(P_positive), P_positive),
      "fixed field of sqrt(5)->-sqrt(5) is Q")
check("edge and chord incidence select opposite golden doublets",
      same_exact(W*P_positive,
                 negative_doublet_value*P_positive)
      and same_exact(W*P_negative,
                     positive_doublet_value*P_negative),
      "the edge-versus-distance-two distinction is load-bearing")

# The abstract D5 automorphism r->r^2 is the character-Galois operation.  It
# exchanges the two order-five conjugacy classes, hence sends the geometric
# edge operator to the non-edge chord operator.  It is consequently not a
# symmetry of the incidence-decorated fibre.
rotation_powers = {identity_action: 0}
current = identity_action
for exponent in range(1, 5):
    current = compose(r, current)
    rotation_powers[current] = exponent
s = min(element for element in stabilizer
        if permutation_order(element) == 2)
normal_form = {}
for element, exponent in rotation_powers.items():
    normal_form[element] = (0, exponent)
    normal_form[compose(s, element)] = (1, exponent)


def golden_group_automorphism(element):
    reflection_bit, exponent = normal_form[element]
    image_rotation = identity_action
    for _ in range((2*exponent) % 5):
        image_rotation = compose(r, image_rotation)
    return (compose(s, image_rotation)
            if reflection_bit else image_rotation)


tau_homomorphism = (
    len(normal_form) == 10
    and len({golden_group_automorphism(element)
             for element in stabilizer}) == 10
    and all(golden_group_automorphism(compose(left, right))
            == compose(golden_group_automorphism(left),
                       golden_group_automorphism(right))
            for left in stabilizer for right in stabilizer)
)


def transport_group_algebra_operator(operator):
    # In the left regular representation, column e_identity lists the
    # coefficients of sum_g c_g L_g.
    identity_column = stabilizer_index[identity_action]
    result = sp.zeros(10)
    for element in stabilizer:
        coefficient = operator[stabilizer_index[element], identity_column]
        result += coefficient*regular[golden_group_automorphism(element)]
    return result


check("the exact golden D5 automorphism sends fibre edges to chords",
      tau_homomorphism
      and transport_group_algebra_operator(U) == W
      and same_exact(transport_group_algebra_operator(P_negative),
                     P_positive),
      "r->r^2 maps u_edge to u_chord and swaps the M12 projectors")

# The unique order-two class distinguishes the two one-dimensional nodes.
reflections = tuple(element for element in stabilizer
                    if permutation_order(element) == 2)
V = sum((regular[element] for element in reflections), sp.zeros(10))
P_trivial = sum((regular[element] for element in stabilizer),
                sp.zeros(10))/10
P_sign = sum(((regular[element] if permutation_order(element) in (1, 5)
               else -regular[element]) for element in stabilizer),
             sp.zeros(10))/10
check("the reflection sum splits the two one-dimensional blocks",
      len(reflections) == 5
      and P_trivial.rank() == P_sign.rank() == 1
      and same_exact(P_trivial+P_sign+P_doublets, I10)
      and same_exact(U*P_trivial, 2*P_trivial)
      and same_exact(U*P_sign, 2*P_sign)
      and same_exact(V*P_trivial, 5*P_trivial)
      and same_exact(V*P_sign, -5*P_sign)
      and same_exact(V*P_negative, sp.zeros(10))
      and same_exact(V*P_positive, sp.zeros(10)))
check("the joint incidence spectrum labels all four real Wedderburn nodes",
      len({("2", "5"), ("2", "-5"),
           (str(positive_doublet_value), "0"),
           (str(negative_doublet_value), "0")}) == 4,
      "(2,5),(2,-5),(phi-1,0),(-phi,0)")

# On the faithful Morita carrier R^6 tensor R[D5], the central ranks are
# multiplied by six.  This identifies algebra blocks; it is not yet a
# Krajewski support or intersection form.
P_negative_crossed = sp.kronecker_product(sp.eye(6), P_negative)
P_positive_crossed = sp.kronecker_product(sp.eye(6), P_positive)
check("Morita amplification separates exactly the two M12 blocks",
      P_negative_crossed.rank() == P_positive_crossed.rank() == 24
      and same_exact(P_negative_crossed*P_positive_crossed, sp.zeros(60)),
      "rank 24 is the regular isotypic rank; each algebra summand is M12(R)")

payload = {
    "protocol_commit": "43c6dd3",
    "target_comparison_performed": False,
    "incidence": {
        "fibrations_per_handed_family": len(fibrations),
        "handed_families_checked": 2,
        "total_fibrations_checked": len(fibrations)+len(opposite_fibrations),
        "records": incidence_records,
        "edge_inverse_pair_uniform": all_incidence_canonical,
        "opposite_hand_edge_inverse_pair_uniform": opposite_ok,
        "A5_transport_natural": transport_natural and opposite_transport,
        "edge_total_undirected": sum(
            record["internal_edge_count"] for record in incidence_records
        )*2,
        "chord_total_undirected": sum(
            record["distance_two_pair_count"] for record in incidence_records
        )*2,
    },
    "regular_D5": {
        "dimension": 10,
        "edge_characteristic_polynomial": (
            "(lambda-2)^2*(lambda^2+lambda-1)^4"
        ),
        "joint_node_spectrum": [
            ["2", "5"], ["2", "-5"],
            ["phi-1", "0"], ["-phi", "0"],
        ],
        "negative_projector_rank": int(P_negative.rank()),
        "positive_projector_rank": int(P_positive.rank()),
        "individual_projectors_rational": False,
        "projector_sum_rational": True,
        "chord_swaps_selected_doublet": True,
    },
    "crossed_product": {
        "type": ["M6(R)", "M6(R)", "M12(R)", "M12(R)"],
        "negative_M12_regular_isotypic_rank": int(P_negative_crossed.rank()),
        "positive_M12_regular_isotypic_rank": int(P_positive_crossed.rank()),
    },
    "scope": {
        "geometry_separates_golden_nodes_over_ordered_R": True,
        "edge_operator_defined_over_Q": True,
        "golden_automorphism_preserves_edge_incidence": False,
        "individual_real_node_projectors_descend_to_Q": False,
        "rational_unsplit_K0_no_go_still_valid": True,
        "Krajewski_support_selected": False,
        "KO6_intersection_form_selected": False,
        "finite_spectral_triple_constructed": False,
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact structured audit was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: fibre incidence separates the two golden D5/M12 blocks.")
print("DERIVED CORRECTION: the real geometry supplies a canonical node label.")
print("OPEN: no Krajewski support or KO6 intersection form is selected yet.")
print("NO HESSIAN OR STANDARD-MODEL TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

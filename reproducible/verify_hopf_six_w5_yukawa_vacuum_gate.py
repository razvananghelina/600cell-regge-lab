#!/usr/bin/env python3
"""Exact maximal-span gate for a W5 simplex vacuum on the 936 carrier.

The field, complete Phase-1 census and Phase-2 decision boundary were frozen
in commits b9623c4 and c148821 before any vacuum commutant was evaluated.
All calculations below are exact.  No physical target is used.
"""

from collections import Counter
from itertools import combinations, permutations, product
import json
from pathlib import Path

import numpy as np
import sympy as sp


OUTPUT = Path(__file__).with_name("hopf_six_w5_yukawa_vacuum_gate.json")
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
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left+1, len(permutation))
    )
    return inversions % 2


def permutation_order(permutation):
    identity = tuple(range(len(permutation)))
    current = identity
    for order in range(1, 61):
        current = compose(permutation, current)
        if current == identity:
            return order
    raise RuntimeError("permutation order exceeds group order")


print("="*78)
print("W5 SIMPLEX-VACUUM MAXIMAL DIRAC-SPAN GATE")
print("="*78)

# -------------------------------------------------------------------------
# Independent exact combinatorial model of A5 and its six Sylow-5 subgroups.
# The six points are the same abstract A5/D5 orbit as the derived fibrations,
# without using rotation eigenvectors or imported six-label permutations.
# -------------------------------------------------------------------------
a5 = tuple(
    permutation for permutation in permutations(range(5))
    if parity(permutation) == 0
)
identity = tuple(range(5))
a5_set = set(a5)
check("the even permutations form the exact group A5", len(a5) == 60)
check(
    "the permutation model is closed under products and inverses",
    all(compose(left, right) in a5_set for left in a5 for right in a5)
    and all(inverse(element) in a5_set for element in a5),
)

order_five = tuple(element for element in a5
                   if permutation_order(element) == 5)
sylow5 = set()
for generator in order_five:
    subgroup = {identity}
    current = identity
    for _ in range(4):
        current = compose(generator, current)
        subgroup.add(current)
    sylow5.add(frozenset(subgroup))
sylow5 = tuple(sorted(sylow5, key=lambda subgroup: sorted(subgroup)))
check(
    "A5 has exactly six Sylow-5 subgroups",
    len(order_five) == 24 and len(sylow5) == 6
    and all(len(subgroup) == 5 for subgroup in sylow5),
)

sylow_index = {subgroup: index for index, subgroup in enumerate(sylow5)}
six_point_action = []
for group_element in a5:
    group_inverse = inverse(group_element)
    action = []
    for subgroup in sylow5:
        conjugate = frozenset(
            compose(compose(group_element, element), group_inverse)
            for element in subgroup
        )
        action.append(sylow_index[conjugate])
    six_point_action.append(tuple(action))
check(
    "conjugation gives a faithful transitive six-point action",
    len(set(six_point_action)) == 60
    and {action[0] for action in six_point_action} == set(range(6)),
)

stabilizers = tuple(
    tuple(a5[index] for index, action in enumerate(six_point_action)
          if action[point] == point)
    for point in range(6)
)
check(
    "all six simplex vacua have stabilizer D5 of order ten",
    all(len(stabilizer) == 10
        and Counter(permutation_order(element) for element in stabilizer)
        == Counter({1: 1, 2: 5, 5: 4})
        for stabilizer in stabilizers),
    "orbit-stabilizer: 60/10=6; order census 1+5+4",
)

# -------------------------------------------------------------------------
# Exact characters and restriction of the four A5 node modules to the
# stabilizer D5.  Class conventions agree with edge/chord spectral labels.
# -------------------------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1+sqrt5)/2
phi_bar = (1-sqrt5)/2
beta = phi-1
alpha = -phi

a5_names = ("1", "3", "3p", "4", "5")
a5_class_sizes = (1, 15, 20, 12, 12)
a5_characters = {
    "1": (1, 1, 1, 1, 1),
    "3": (3, -1, 0, phi, phi_bar),
    "3p": (3, -1, 0, phi_bar, phi),
    "4": (4, 0, 1, -1, -1),
    "5": (5, 1, -1, 0, 0),
}
d5_names = (
    "trivial", "reflection_sign", "positive_doublet", "negative_doublet"
)
d5_dimensions = (1, 1, 2, 2)
d5_class_sizes = (1, 5, 2, 2)
d5_characters = {
    "trivial": (1, 1, 1, 1),
    "reflection_sign": (1, -1, 1, 1),
    "positive_doublet": (2, 0, beta, alpha),
    "negative_doublet": (2, 0, alpha, beta),
}
restricted_a5_irreps = {
    name: (character[0], character[1], character[3], character[4])
    for name, character in a5_characters.items()
}


def d5_inner(left, right):
    return sp.simplify(sp.Rational(1, 10)*sum(
        size*a*b for size, a, b in zip(d5_class_sizes, left, right)
    ))


node_a5_irreps = (
    ("1", "5"),
    ("3", "3p"),
    ("3", "4", "5"),
    ("3p", "4", "5"),
)
restriction_rows = []
for irreps in node_a5_irreps:
    restricted_character = tuple(sp.simplify(sum(
        restricted_a5_irreps[name][class_index] for name in irreps
    )) for class_index in range(4))
    restriction_rows.append(tuple(
        int(d5_inner(restricted_character, d5_characters[d5_name]))
        for d5_name in d5_names
    ))
restriction_rows = tuple(restriction_rows)
expected_restrictions = (
    (2, 0, 1, 1),
    (0, 2, 1, 1),
    (1, 1, 3, 2),
    (1, 1, 2, 3),
)
check(
    "all four node restrictions to the vacuum D5 are exact",
    restriction_rows == expected_restrictions,
    f"multiplicity rows={restriction_rows}",
)

node_sizes = tuple(sum(
    multiplicity*dimension
    for multiplicity, dimension in zip(row, d5_dimensions)
) for row in restriction_rows)
check("the restricted models retain dimensions 6,6,12,12",
      node_sizes == (6, 6, 12, 12))

hom_d5 = sp.Matrix([
    [sum(a*b for a, b in zip(left, right))
     for right in restriction_rows]
    for left in restriction_rows
])
hom_a5 = sp.Matrix([
    [2, 0, 1, 1],
    [0, 2, 1, 1],
    [1, 1, 3, 2],
    [1, 1, 2, 3],
])
w5_couplings = sp.Matrix([
    [4, 2, 6, 6],
    [2, 4, 6, 6],
    [6, 6, 12, 12],
    [6, 6, 12, 12],
])
check(
    "Frobenius evaluation gives Hom_D5 = Hom_A5 plus W5 couplings",
    hom_d5 == hom_a5+w5_couplings,
    f"Hom_D5={hom_d5.tolist()}",
)

# Exact block models for full Hom_D5 spaces.  Each node is ordered by D5
# irrep, then multiplicity copy, then irrep coordinate.
node_irrep_slices = []
for node, multiplicities in enumerate(restriction_rows):
    offset = 0
    slices = {}
    for irrep, (multiplicity, dimension) in enumerate(
            zip(multiplicities, d5_dimensions)):
        slices[irrep] = (offset, multiplicity, dimension)
        offset += multiplicity*dimension
    if offset != node_sizes[node]:
        raise RuntimeError("restricted-model dimension mismatch")
    node_irrep_slices.append(slices)


def hom_d5_basis(target, source):
    matrices = []
    labels = []
    for irrep, dimension in enumerate(d5_dimensions):
        source_start, source_mult, _ = node_irrep_slices[source][irrep]
        target_start, target_mult, _ = node_irrep_slices[target][irrep]
        for target_copy in range(target_mult):
            for source_copy in range(source_mult):
                matrix = np.zeros(
                    (node_sizes[target], node_sizes[source]), dtype=int
                )
                for coordinate in range(dimension):
                    row = target_start+target_copy*dimension+coordinate
                    column = source_start+source_copy*dimension+coordinate
                    matrix[row, column] = 1
                matrices.append(matrix)
                labels.append((irrep, target_copy, source_copy))
    return tuple(matrices), tuple(labels)


check(
    "explicit matrix units span every full D5 intertwiner space",
    all(
        len(hom_d5_basis(target, source)[0]) == int(hom_d5[target, source])
        for source in range(4) for target in range(4)
    ),
)

# -------------------------------------------------------------------------
# Eight spectral readings and maximal legal Hom_D5 spans at one vacuum.
# The other five vacua are A5-conjugate, and the negative signed orbit has the
# same stabilizer, so ranks and commutants are identical on all twelve.
# -------------------------------------------------------------------------
nodes = (
    {"name": "trivial", "u_rank": 2, "v_rank": 2},
    {"name": "reflection_sign", "u_rank": 2, "v_rank": 0},
    {"name": "positive_doublet", "u_rank": 1, "v_rank": 1},
    {"name": "negative_doublet", "u_rank": 0, "v_rank": 1},
)
unordered_pairs = tuple(combinations(range(4), 2))
readings = []
for priority, u_direction, v_direction in product(
        (("u", "v"), ("v", "u")), (1, -1), (1, -1)):
    direction = {"u": u_direction, "v": v_direction}

    def spectral_key(node):
        return tuple(direction[coordinate]*nodes[node][coordinate+"_rank"]
                     for coordinate in priority)

    order = tuple(sorted(range(4), key=spectral_key))
    position = {node: rank for rank, node in enumerate(order)}
    positive_cells = frozenset(
        (left, right) if position[left] < position[right] else (right, left)
        for left, right in unordered_pairs
    )
    legal_blocks = []
    for source in positive_cells:
        i, j = source
        for second_positive in positive_cells:
            k, ell = second_positive
            if i != ell and j != k:
                continue
            if i == ell:
                link = tuple(sorted((j, k)))
                side = "right"
            else:
                link = tuple(sorted((i, ell)))
                side = "left"
            legal_blocks.append((link, side))
    links = tuple(sorted({link for link, _ in legal_blocks}))
    readings.append({
        "priority": priority,
        "u_direction": u_direction,
        "v_direction": v_direction,
        "order": order,
        "links": links,
        "legal_blocks": tuple(legal_blocks),
    })

check(
    "all eight legal supports are independently reconstructed",
    len(readings) == 8
    and len({record["order"] for record in readings}) == 8
    and all(len(record["links"]) == 3
            and len(record["legal_blocks"]) == 8 for record in readings),
)

variable_offsets = []
variable_count = 0
for size in node_sizes:
    variable_offsets.append(variable_count)
    variable_count += size*size
check("the full real algebra again has 360 coefficients",
      variable_count == 360)


def variable(node, row, column):
    return variable_offsets[node]+row*node_sizes[node]+column


def maximal_constraint(links):
    rows = []
    for source, target in links:
        matrices, _ = hom_d5_basis(target, source)
        for T in matrices:
            for row in range(node_sizes[target]):
                for column in range(node_sizes[source]):
                    equation = {}
                    for inner in range(node_sizes[target]):
                        coefficient = int(T[inner, column])
                        if coefficient:
                            index = variable(target, row, inner)
                            equation[index] = equation.get(index, 0)+coefficient
                    for inner in range(node_sizes[source]):
                        coefficient = int(T[row, inner])
                        if coefficient:
                            index = variable(source, inner, column)
                            equation[index] = equation.get(index, 0)-coefficient
                    if equation:
                        rows.append(equation)
            TT = T.T
            for row in range(node_sizes[source]):
                for column in range(node_sizes[target]):
                    equation = {}
                    for inner in range(node_sizes[source]):
                        coefficient = int(TT[inner, column])
                        if coefficient:
                            index = variable(source, row, inner)
                            equation[index] = equation.get(index, 0)+coefficient
                    for inner in range(node_sizes[target]):
                        coefficient = int(TT[row, inner])
                        if coefficient:
                            index = variable(target, inner, column)
                            equation[index] = equation.get(index, 0)-coefficient
                    if equation:
                        rows.append(equation)
    constraint = sp.MutableSparseMatrix(len(rows), variable_count, {})
    for row_index, equation in enumerate(rows):
        for column, coefficient in equation.items():
            constraint[row_index, column] = coefficient
    return constraint


def residual_irrep_projector(irrep):
    vector = sp.zeros(variable_count, 1)
    for node in range(4):
        start, multiplicity, dimension = node_irrep_slices[node][irrep]
        for coordinate in range(multiplicity*dimension):
            absolute = start+coordinate
            vector[variable(node, absolute, absolute)] = 1
    return vector


residual_projectors = tuple(
    residual_irrep_projector(irrep) for irrep in range(4)
)
global_scalar = sp.zeros(variable_count, 1)
for node, size in enumerate(node_sizes):
    for coordinate in range(size):
        global_scalar[variable(node, coordinate, coordinate)] = 1
projector_supports = tuple(
    {index for index in range(variable_count) if projector[index] != 0}
    for projector in residual_projectors
)
check(
    "four nonzero orthogonal residual-D5 projectors resolve the identity",
    all(projector != sp.zeros(variable_count, 1)
        and projector != global_scalar for projector in residual_projectors)
    and all(projector_supports[left].isdisjoint(projector_supports[right])
            for left in range(4) for right in range(left+1, 4))
    and sum(residual_projectors, sp.zeros(variable_count, 1)) == global_scalar,
)

cache = {}
records = []
all_projectors_commute = True
for reading_index, record in enumerate(readings):
    links = record["links"]
    if links not in cache:
        constraint = maximal_constraint(links)
        rank = int(constraint.rank())
        projector_checks = tuple(
            constraint*projector == sp.zeros(constraint.rows, 1)
            for projector in residual_projectors
        )
        cache[links] = {
            "rows": int(constraint.rows),
            "rank": rank,
            "commutant_dimension": variable_count-rank,
            "projector_checks": projector_checks,
        }
    result = cache[links]
    all_projectors_commute &= all(result["projector_checks"])
    records.append({
        "reading_index": reading_index,
        "priority": list(record["priority"]),
        "u_direction": record["u_direction"],
        "v_direction": record["v_direction"],
        "node_order": [nodes[node]["name"] for node in record["order"]],
        "links": [list(link) for link in links],
        "link_Hom_D5_dimensions": {
            f"{source}-{target}": int(hom_d5[source, target])
            for source, target in links
        },
        "maximal_constraint_rows": result["rows"],
        "maximal_constraint_rank": result["rank"],
        "maximal_commutant_dimension": result["commutant_dimension"],
        "all_four_residual_projectors_commute": all(
            result["projector_checks"]
        ),
    })

check(
    "all four residual-D5 projectors commute on every maximal legal span",
    all_projectors_commute,
)

dimension_multiset = Counter(
    record["maximal_commutant_dimension"] for record in records
)
rank_multiset = Counter(record["maximal_constraint_rank"] for record in records)
check(
    "the complete affine W5 vacuum span is nonconnected for all readings",
    all(record["maximal_commutant_dimension"] >= 4 for record in records),
    f"commutant dimensions={dict(dimension_multiset)}; "
    f"ranks={dict(rank_multiset)}",
)
check(
    "the connected hit fraction is exactly 0/8 readings and 0/12 vacua",
    sum(record["maximal_commutant_dimension"] == 1
        for record in records) == 0,
)

payload = {
    "protocol_commit": "b9623c4",
    "blind_enumeration_commit": "c148821",
    "physical_target_comparison_performed": False,
    "vacuum_orbits": {
        "positive_points": 6,
        "negative_points": 6,
        "stabilizer": "D5",
        "stabilizer_order": 10,
    },
    "node_restriction_multiplicities_D5": [list(row)
                                            for row in restriction_rows],
    "Hom_A5": [[int(hom_a5[row, column]) for column in range(4)]
               for row in range(4)],
    "W5_linear_couplings": [
        [int(w5_couplings[row, column]) for column in range(4)]
        for row in range(4)
    ],
    "maximal_affine_vacuum_Hom_D5": [
        [int(hom_d5[row, column]) for column in range(4)]
        for row in range(4)
    ],
    "reading_count": len(records),
    "distinct_legal_link_sets": len(cache),
    "connected_readings": sum(
        record["maximal_commutant_dimension"] == 1 for record in records
    ),
    "connected_signed_vacua": 0,
    "maximal_commutant_dimension_multiset": {
        str(key): value for key, value in sorted(dimension_multiset.items())
    },
    "maximal_constraint_rank_multiset": {
        str(key): value for key, value in sorted(rank_multiset.items())
    },
    "readings": records,
    "verdict": (
        "DERIVED LINEAR-FIELD NO-GO: a simplex vacuum leaves D5 unbroken. "
        "The complete affine evaluation span is contained in Hom_D5 and "
        "therefore commutes with four nontrivial residual-isotypic algebra "
        "projectors. Even granting every Hom_D5 tensor on every legal link, "
        "all eight readings and all twelve signed vacua are nonconnected."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact structured vacuum gate was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED LINEAR-FIELD NO-GO: connected readings 0/8; signed vacua 0/12.")
print("NO MATTER OR STANDARD-MODEL TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

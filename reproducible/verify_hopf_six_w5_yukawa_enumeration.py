#!/usr/bin/env python3
"""Blind exact census of W5-covariant inter-node Dirac couplings.

The complete census and the ban on vacuum/connectedness evaluation were
frozen in protocol commit b9623c4.  This STEP 1 uses only exact characters
and central first-order supports.  No vacuum matrix or physical target is
constructed.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import sympy as sp


OUTPUT = Path(__file__).with_name("hopf_six_w5_yukawa_enumeration.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("="*78)
print("BLIND W5-COVARIANT DIRAC-COUPLING CENSUS")
print("="*78)

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


def a5_inner(left, right):
    return sp.simplify(sp.Rational(1, 60)*sum(
        size*a*b for size, a, b in zip(a5_class_sizes, left, right)
    ))


check(
    "the exact A5 characters are orthonormal",
    all(
        a5_inner(a5_characters[left], a5_characters[right])
        == (1 if left == right else 0)
        for left in a5_names for right in a5_names
    ),
)

d5_names = (
    "trivial", "reflection_sign", "positive_doublet", "negative_doublet"
)
d5_class_sizes = (1, 5, 2, 2)
d5_characters = {
    "trivial": (1, 1, 1, 1),
    "reflection_sign": (1, -1, 1, 1),
    "positive_doublet": (2, 0, beta, alpha),
    "negative_doublet": (2, 0, alpha, beta),
}
restricted_a5 = {
    name: (character[0], character[1], character[3], character[4])
    for name, character in a5_characters.items()
}


def d5_inner(left, right):
    return sp.simplify(sp.Rational(1, 10)*sum(
        size*a*b for size, a, b in zip(d5_class_sizes, left, right)
    ))


induction_rows = tuple(
    tuple(int(d5_inner(d5_characters[d5_name], restricted_a5[a5_name]))
          for a5_name in a5_names)
    for d5_name in d5_names
)
expected_rows = (
    (1, 0, 0, 0, 1),
    (0, 1, 1, 0, 0),
    (0, 1, 0, 1, 1),
    (0, 0, 1, 1, 1),
)
check(
    "the four induced node modules are reconstructed exactly",
    induction_rows == expected_rows,
    "V0=1+5; V1=3+3'; V2=3+4+5; V3=3'+4+5",
)


def sum_characters(multiplicities):
    return tuple(sp.simplify(sum(
        multiplicity*a5_characters[name][class_index]
        for name, multiplicity in zip(a5_names, multiplicities)
    )) for class_index in range(5))


node_characters = tuple(sum_characters(row) for row in induction_rows)
node_dimensions = tuple(int(character[0]) for character in node_characters)
check("the exact node dimensions are 6,6,12,12",
      node_dimensions == (6, 6, 12, 12))

w5_character = a5_characters["5"]
coupling_matrix = sp.zeros(4)
for left in range(4):
    for right in range(4):
        product_character = tuple(
            node_characters[left][class_index]
            *node_characters[right][class_index]
            for class_index in range(5)
        )
        multiplicity = sp.simplify(a5_inner(w5_character, product_character))
        if not multiplicity.is_Integer or multiplicity < 0:
            raise RuntimeError("covariant coupling multiplicity is invalid")
        coupling_matrix[left, right] = multiplicity

check(
    "all W5-covariant coupling multiplicities are nonnegative integers",
    all(coupling_matrix[row, column].is_Integer
        and coupling_matrix[row, column] >= 0
        for row in range(4) for column in range(4)),
)
check(
    "the complete coupling matrix is symmetric",
    coupling_matrix == coupling_matrix.T,
)

offdiagonal_pairs = tuple(combinations(range(4), 2))
offdiagonal_multiplicities = tuple(
    int(coupling_matrix[left, right]) for left, right in offdiagonal_pairs
)
offdiagonal_multiset = Counter(offdiagonal_multiplicities)

# Reconstruct the eight derived gradings and their complete central
# first-order link supports.  This is support enumeration only: no tensor is
# evaluated at a field value.
nodes = (
    {"name": "trivial", "u_rank": 2, "v_rank": 2},
    {"name": "reflection_sign", "u_rank": 2, "v_rank": 0},
    {"name": "positive_doublet", "u_rank": 1, "v_rank": 1},
    {"name": "negative_doublet", "u_rank": 0, "v_rank": 1},
)
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
        for left, right in offdiagonal_pairs
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
                changed_side = "right"
            else:
                link = tuple(sorted((i, ell)))
                changed_side = "left"
            legal_blocks.append((link, changed_side))
    links = tuple(sorted({link for link, _ in legal_blocks}))
    multiplicities = tuple(
        int(coupling_matrix[left, right]) for left, right in links
    )
    if any(value == 0 for value in multiplicities):
        status = "zero_link"
    elif all(value == 1 for value in multiplicities):
        status = "all_unique"
    else:
        status = "ambiguous"
    readings.append({
        "priority": priority,
        "u_direction": u_direction,
        "v_direction": v_direction,
        "order": order,
        "positive_cell_count": len(positive_cells),
        "legal_block_count": len(legal_blocks),
        "links": links,
        "multiplicities": multiplicities,
        "status": status,
        "both_sides": all(
            {link for link, side in legal_blocks if side == tensor_side}
            == set(links)
            for tensor_side in ("left", "right")
        ),
    })

check(
    "all eight derived readings have the complete legal support",
    len(readings) == 8
    and len({record["order"] for record in readings}) == 8
    and all(record["positive_cell_count"] == 6
            and record["legal_block_count"] == 8
            and len(record["links"]) == 3
            and record["both_sides"] for record in readings),
)

status_counts = Counter(record["status"] for record in readings)
check(
    "the reading classification is exhaustive",
    sum(status_counts.values()) == 8
    and set(status_counts) <= {"all_unique", "zero_link", "ambiguous"},
)

payload = {
    "protocol_commit": "b9623c4",
    "phase": "STEP 1 target-blind covariant-coupling census",
    "vacuum_evaluated": False,
    "connectedness_evaluated": False,
    "physical_target_comparison_performed": False,
    "node_order": list(d5_names),
    "node_dimensions": list(node_dimensions),
    "W5_coupling_multiplicity_matrix": [
        [int(coupling_matrix[row, column]) for column in range(4)]
        for row in range(4)
    ],
    "offdiagonal_pairs": [list(pair) for pair in offdiagonal_pairs],
    "offdiagonal_multiplicities": list(offdiagonal_multiplicities),
    "offdiagonal_multiplicity_multiset": {
        str(key): value for key, value in sorted(offdiagonal_multiset.items())
    },
    "reading_status_counts": dict(sorted(status_counts.items())),
    "readings": [
        {
            "reading_index": index,
            "priority": list(record["priority"]),
            "u_direction": record["u_direction"],
            "v_direction": record["v_direction"],
            "node_order": [nodes[node]["name"] for node in record["order"]],
            "links": [list(link) for link in record["links"]],
            "W5_coupling_multiplicities": list(record["multiplicities"]),
            "classification": record["status"],
        }
        for index, record in enumerate(readings)
    ],
    "verdict": (
        "TARGET-BLIND ENUMERATION ONLY. The complete W5-covariant coupling "
        "multiplicity matrix and all eight legal-link readings were frozen "
        "before evaluating a Hopf vacuum or a Dirac commutant."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the target-blind structured census was written", OUTPUT.exists())

print("\nComplete W5 coupling matrix:")
print(coupling_matrix)
print(f"Off-diagonal multiset: {dict(offdiagonal_multiset)}")
print(f"Reading classification: {dict(status_counts)}")
print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("STEP 1 ONLY: NO VACUUM OR CONNECTEDNESS WAS EVALUATED.")
raise SystemExit(0 if passed == tests else 1)

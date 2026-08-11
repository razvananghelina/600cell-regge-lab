#!/usr/bin/env python3
"""Universal connectedness bound for A5-equivariant first-order Dirac data.

The hypothesis, complete candidate space and kill boundary were frozen in
protocol commit c2c32df.  The calculation uses exact Q(sqrt(5)) character
arithmetic and exact rational linear algebra.  No physical target is used.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
import sympy as sp


OUTPUT = Path(__file__).with_name(
    "hopf_six_equivariant_dirac_universal.json"
)
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
print("UNIVERSAL A5-EQUIVARIANT DIRAC CONNECTEDNESS BOUND")
print("="*78)

# -------------------------------------------------------------------------
# Exact real representation theory.  A5 class order is
# (1A,2A,3A,5A,5B).  Squaring exchanges the two order-five classes.
# -------------------------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1+sqrt5)/2
phi_bar = (1-sqrt5)/2
beta = phi-1
alpha = -phi

a5_names = ("1", "3", "3p", "4", "5")
a5_dims = (1, 3, 3, 4, 5)
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
    "the exact A5 character table is orthonormal",
    all(
        a5_inner(a5_characters[left], a5_characters[right])
        == (1 if left == right else 0)
        for left in a5_names for right in a5_names
    ),
)

square_class = (0, 0, 2, 4, 3)
frobenius_schur = {}
for name, character in a5_characters.items():
    frobenius_schur[name] = sp.simplify(sp.Rational(1, 60)*sum(
        size*character[square_class[class_index]]
        for class_index, size in enumerate(a5_class_sizes)
    ))
check(
    "all five A5 irreps are of real type",
    all(value == 1 for value in frobenius_schur.values()),
    f"Frobenius-Schur indicators={frobenius_schur}",
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


check(
    "the exact D5 character table is orthonormal",
    all(
        d5_inner(d5_characters[left], d5_characters[right])
        == (1 if left == right else 0)
        for left in d5_names for right in d5_names
    ),
)

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
    "Frobenius reciprocity reconstructs all four induced modules",
    induction_rows == expected_rows,
    "V0=1+5; V1=3+3'; V2=3+4+5; V3=3'+4+5",
)

hom_gram = sp.Matrix([
    [sum(a*b for a, b in zip(left, right)) for right in induction_rows]
    for left in induction_rows
])
expected_gram = sp.Matrix([
    [2, 0, 1, 1],
    [0, 2, 1, 1],
    [1, 1, 3, 2],
    [1, 1, 2, 3],
])
check(
    "the complete real equivariant Hom Gram matrix is exact",
    hom_gram == expected_gram,
    f"Gram={hom_gram.tolist()}",
)

# Fixed multiplicity-one orthogonal sum models.  Because every A5 irrep is
# real type, identity maps on common summands form a real Hom basis.
node_irreps = (
    ("1", "5"),
    ("3", "3p"),
    ("3", "4", "5"),
    ("3p", "4", "5"),
)
irrep_dimensions = dict(zip(a5_names, a5_dims))
node_sizes = tuple(sum(irrep_dimensions[name] for name in irreps)
                   for irreps in node_irreps)
node_slices = []
for irreps in node_irreps:
    offset = 0
    slices = {}
    for irrep in irreps:
        dimension = irrep_dimensions[irrep]
        slices[irrep] = (offset, offset+dimension)
        offset += dimension
    node_slices.append(slices)
check("the real node sizes are 6,6,12,12", node_sizes == (6, 6, 12, 12))


def hom_basis(target, source):
    """Exact real basis of Hom_A5(V_source,V_target)."""
    matrices = []
    labels = []
    for irrep in a5_names:
        if irrep not in node_slices[source] or irrep not in node_slices[target]:
            continue
        source_start, source_stop = node_slices[source][irrep]
        target_start, target_stop = node_slices[target][irrep]
        if source_stop-source_start != target_stop-target_start:
            raise RuntimeError("common irreducible dimensions disagree")
        matrix = np.zeros((node_sizes[target], node_sizes[source]), dtype=int)
        for offset in range(source_stop-source_start):
            matrix[target_start+offset, source_start+offset] = 1
        matrices.append(matrix)
        labels.append(irrep)
    return tuple(matrices), tuple(labels)


all_hom_bases_correct = True
for source in range(4):
    for target in range(4):
        matrices, _ = hom_basis(target, source)
        all_hom_bases_correct &= len(matrices) == int(hom_gram[target, source])
check(
    "identity maps on common real irreps span every equivariant Hom space",
    all_hom_bases_correct,
)

# -------------------------------------------------------------------------
# Reconstruct all eight gradings and every centrally legal odd cell block.
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
            target = (ell, k)  # J-transpose of second_positive
            if i == ell:
                link = tuple(sorted((j, k)))
                changed_side = "right"
            else:
                link = tuple(sorted((i, ell)))
                changed_side = "left"
            legal_blocks.append({
                "source": source,
                "target": target,
                "link": link,
                "changed_side": changed_side,
            })
    links = tuple(sorted({block["link"] for block in legal_blocks}))
    readings.append({
        "priority": priority,
        "u_direction": u_direction,
        "v_direction": v_direction,
        "order": order,
        "positive_cells": positive_cells,
        "legal_blocks": tuple(legal_blocks),
        "links": links,
    })

check(
    "all eight spectral readings and all legal odd positions are reconstructed",
    len(readings) == 8
    and len({record["order"] for record in readings}) == 8
    and all(len(record["positive_cells"]) == 6
            and len(record["legal_blocks"]) == 8
            and len(record["links"]) == 3
            for record in readings),
    "each reading: 6 positive cells, 8 legal odd blocks, 3 node links",
)
check(
    "every link occurs on both tensor sides",
    all(
        {block["link"] for block in record["legal_blocks"]
         if block["changed_side"] == side} == set(record["links"])
        for record in readings for side in ("left", "right")
    ),
)
check(
    "the ambiguous V2-V3 Hom space is absent from every legal support",
    all((2, 3) not in record["links"] for record in readings),
)

# -------------------------------------------------------------------------
# Maximal commutator constraints.  For every legal link, impose equations for
# an entire basis of Hom_A5, rather than for one selected linear combination.
# Any individual equivariant D imposes a subset/linear combination of these
# equations, so its algebra commutant contains the maximal common commutant.
# -------------------------------------------------------------------------
variable_offsets = []
variable_count = 0
for size in node_sizes:
    variable_offsets.append(variable_count)
    variable_count += size*size
check("the represented real algebra has 360 matrix coefficients",
      variable_count == 360)


def variable(node, row, column):
    return variable_offsets[node]+row*node_sizes[node]+column


def maximal_constraint(links):
    rows = []
    hom_labels = {}
    for source, target in links:
        matrices, labels = hom_basis(target, source)
        hom_labels[(source, target)] = labels
        for T in matrices:
            # A_target*T - T*A_source = 0.
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

            # The self-adjoint Dirac block also imposes
            # A_source*T^t - T^t*A_target = 0.
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
    return constraint, hom_labels


def occurrence_components(links):
    occurrences = tuple(
        (node, irrep) for node, irreps in enumerate(node_irreps)
        for irrep in irreps
    )
    graph = {occurrence: set() for occurrence in occurrences}
    for source, target in links:
        _, labels = hom_basis(target, source)
        for irrep in labels:
            left = (source, irrep)
            right = (target, irrep)
            graph[left].add(right)
            graph[right].add(left)
    components = []
    unseen = set(occurrences)
    while unseen:
        seed = min(unseen)
        stack = [seed]
        component = set()
        while stack:
            current = stack.pop()
            if current in component:
                continue
            component.add(current)
            unseen.discard(current)
            stack.extend(graph[current]-component)
        components.append(frozenset(component))
    return tuple(components)


def component_projector(component):
    vector = sp.zeros(variable_count, 1)
    for node, irrep in component:
        start, stop = node_slices[node][irrep]
        for coordinate in range(start, stop):
            vector[variable(node, coordinate, coordinate)] = 1
    return vector


cache = {}
records = []
all_witnesses_exact = True
for index, record in enumerate(readings):
    links = record["links"]
    if links not in cache:
        constraint, hom_labels = maximal_constraint(links)
        rank = int(constraint.rank())
        kernel_dimension = variable_count-rank
        components = occurrence_components(links)
        witnesses = tuple(component_projector(component)
                          for component in components)
        scalar = sum(witnesses, sp.zeros(variable_count, 1))
        witness_exact = all(constraint*witness == sp.zeros(constraint.rows, 1)
                            for witness in witnesses)
        non_scalar_witness = next(
            (witness for witness in witnesses if witness != scalar), None
        )
        witness_exact &= (
            non_scalar_witness is not None
            and constraint*non_scalar_witness
            == sp.zeros(constraint.rows, 1)
        )
        cache[links] = {
            "constraint_rows": int(constraint.rows),
            "rank": rank,
            "kernel_dimension": kernel_dimension,
            "hom_labels": hom_labels,
            "components": components,
            "component_lower_bound": len(components),
            "witness_exact": witness_exact,
            "non_scalar_witness_support": [
                [node, irrep]
                for node, irrep in sorted(next(
                    component for component, witness
                    in zip(components, witnesses)
                    if witness == non_scalar_witness
                ))
            ],
        }
    result = cache[links]
    all_witnesses_exact &= result["witness_exact"]
    records.append({
        "reading_index": index,
        "priority": list(record["priority"]),
        "u_direction": record["u_direction"],
        "v_direction": record["v_direction"],
        "node_order": [nodes[node]["name"] for node in record["order"]],
        "links": [list(link) for link in links],
        "link_Hom_dimensions": {
            f"{source}-{target}": int(hom_gram[source, target])
            for source, target in links
        },
        "maximal_constraint_rows": result["constraint_rows"],
        "maximal_constraint_rank": result["rank"],
        "maximal_commutant_dimension": result["kernel_dimension"],
        "irrep_component_lower_bound": result["component_lower_bound"],
        "explicit_non_scalar_witness_support": (
            result["non_scalar_witness_support"]
        ),
        "witness_checked_exactly": result["witness_exact"],
    })

check(
    "every reading has an exact non-scalar maximal-commutant witness",
    all_witnesses_exact,
)

dimension_multiset = Counter(
    record["maximal_commutant_dimension"] for record in records
)
rank_multiset = Counter(record["maximal_constraint_rank"] for record in records)
check(
    "even the maximal equivariant constraint system is never connected",
    all(record["maximal_commutant_dimension"] > 1 for record in records),
    f"commutant dimensions={dict(dimension_multiset)}; "
    f"ranks={dict(rank_multiset)}",
)
check(
    "all edge-first readings reproduce the previous maximal bounds",
    Counter(
        record["maximal_commutant_dimension"] for record in records
        if record["priority"] == ["u", "v"]
    ) == Counter({109: 2, 141: 2}),
)
check(
    "every reflection-first reading contains the forbidden V0-V1 link",
    all(
        [0, 1] in record["links"]
        and record["link_Hom_dimensions"]["0-1"] == 0
        for record in records if record["priority"] == ["v", "u"]
    ),
)
check(
    "the universal equivariant connectedness hit fraction is 0/8",
    sum(record["maximal_commutant_dimension"] == 1
        for record in records) == 0,
)

payload = {
    "protocol_commit": "c2c32df",
    "physical_target_comparison_performed": False,
    "hypothesis": (
        "fixed 936-state H_off; one of eight derived gradings; full-matrix "
        "first order; odd scalar D natural under the derived A5 action"
    ),
    "node_sizes": list(node_sizes),
    "Hom_Gram": [
        [int(hom_gram[row, column]) for column in range(4)]
        for row in range(4)
    ],
    "reading_count": len(records),
    "distinct_legal_link_sets": len(cache),
    "connected_maximal_spans": sum(
        record["maximal_commutant_dimension"] == 1 for record in records
    ),
    "maximal_commutant_dimension_multiset": {
        str(key): value for key, value in sorted(dimension_multiset.items())
    },
    "maximal_constraint_rank_multiset": {
        str(key): value for key, value in sorted(rank_multiset.items())
    },
    "readings": records,
    "verdict": (
        "DERIVED UNIVERSAL EQUIVARIANT NO-GO: for all eight derived "
        "gradings, the common algebra commutant of the complete legal "
        "A5-equivariant first-order Hom span is non-scalar. Every individual "
        "equivariant D has a commutant containing this one, so unequal "
        "magnitudes, zeros and arbitrary Hom-space combinations cannot "
        "restore connectedness. A continuation requires independently "
        "selected A5-breaking data or a different carrier/algebra."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact structured universal audit was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED UNIVERSAL EQUIVARIANT NO-GO: connected maximal spans 0/8.")
print("NO HESSIAN OR STANDARD-MODEL TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

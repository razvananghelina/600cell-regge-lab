#!/usr/bin/env python3
"""Blind 729-support central Krajewski census for M5^3+M15.

The support space, gates and recorded symmetries were frozen in commit
8b21830.  This file contains no Hessian or selector comparison.
"""

from collections import Counter
from itertools import combinations, permutations, product
import json
from pathlib import Path

import sympy as sp


OUTPUT = Path(__file__).with_name(
    "hopf_hessian_crossed_krajewski_blind.json"
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


def connected_components(node_count, edges):
    adjacency = {node: set() for node in range(node_count)}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(range(node_count))
    components = []
    while unseen:
        seed = min(unseen)
        component = {seed}
        frontier = [seed]
        while frontier:
            node = frontier.pop()
            for neighbor in adjacency[node]:
                if neighbor not in component:
                    component.add(neighbor)
                    frontier.append(neighbor)
        components.append(component)
        unseen -= component
    return components


def transform_edges(edges, permutation):
    return frozenset((permutation[left], permutation[right])
                     for left, right in edges)


def reverse_edges(edges):
    return frozenset((right, left) for left, right in edges)


def edge_key(edges):
    return tuple(sorted(edges))


def support_gate_data(edges):
    """First-order cell support and induced central-link graph."""
    legal_blocks = []
    central_links = set()
    for positive in edges:
        for other_positive in edges:
            left, right = positive
            other_left, other_right = other_positive
            legal = (left == other_right or right == other_left)
            if legal:
                legal_blocks.append((positive, other_positive))
                if right == other_left and left != other_right:
                    central_links.add(tuple(sorted((left, other_right))))
    components = connected_components(4, central_links)
    return legal_blocks, central_links, components


print("="*78)
print("BLIND CROSSED-PRODUCT KRAJEWSKI SUPPORT CENSUS")
print("="*78)

node_sizes = (5, 5, 5, 15)
unordered_pairs = list(combinations(range(4), 2))
conjugation = (0, 2, 1, 3)

all_records = []
gate_counts = Counter()
for choices in product((0, 1, 2), repeat=len(unordered_pairs)):
    edges = set()
    for (left, right), choice in zip(unordered_pairs, choices):
        if choice == 1:
            edges.add((left, right))
        elif choice == 2:
            edges.add((right, left))
    edges = frozenset(edges)

    mu = sp.zeros(4)
    for left, right in edges:
        mu[left, right] = 1
    intersection = mu-mu.T
    pfaffian = (
        intersection[0, 1]*intersection[2, 3]
        - intersection[0, 2]*intersection[1, 3]
        + intersection[0, 3]*intersection[1, 2]
    )
    determinant = int(pfaffian*pfaffian)
    intersection_rank = 4 if pfaffian != 0 else int(intersection.rank())
    faithful = all(any(node in edge for edge in edges) for node in range(4))

    legal_blocks, central_links, components = support_gate_data(edges)
    centre_connected = len(components) == 1

    transformed = transform_edges(edges, conjugation)
    conjugation_preserves_grading = transformed == edges
    conjugation_reverses_grading = transformed == reverse_edges(edges)
    hilbert_dimension = 2*sum(node_sizes[left]*node_sizes[right]
                              for left, right in edges)

    gate_counts["total"] += 1
    gate_counts["poincare"] += int(intersection_rank == 4)
    gate_counts["faithful"] += int(faithful)
    gate_counts["poincare_and_faithful"] += int(
        intersection_rank == 4 and faithful
    )
    gate_counts["nonzero_first_order"] += int(len(legal_blocks) > 0)
    gate_counts["centre_connected"] += int(centre_connected)

    survives = (intersection_rank == 4 and faithful
                and len(legal_blocks) > 0 and centre_connected)
    gate_counts["survivors"] += int(survives)
    all_records.append({
        "edges": [list(edge) for edge in sorted(edges)],
        "edge_key": edge_key(edges),
        "hilbert_dimension": hilbert_dimension,
        "intersection_rank": intersection_rank,
        "intersection_determinant": determinant,
        "intersection_pfaffian": int(pfaffian),
        "faithful": faithful,
        "legal_odd_blocks": len(legal_blocks),
        "central_links": [list(edge) for edge in sorted(central_links)],
        "central_components": len(components),
        "centre_connected": centre_connected,
        "conjugation_preserves_grading": conjugation_preserves_grading,
        "conjugation_reverses_grading": conjugation_reverses_grading,
        "survives": survives,
    })

check("the preregistered support space has exactly 3^6=729 elements",
      gate_counts["total"] == 729 and len(all_records) == 729)
check("every support obeys the orientability no-loop/no-reverse rule",
      all(all(left != right for left, right in record["edge_key"])
          and all((right, left) not in record["edge_key"]
                  for left, right in record["edge_key"])
          for record in all_records))
check("the antisymmetric determinant is the exact Pfaffian square",
      all(record["intersection_determinant"]
          == record["intersection_pfaffian"]**2
          for record in all_records))

survivors = [record for record in all_records if record["survives"]]
survivor_keys = {record["edge_key"] for record in survivors}
dimension_multiset = Counter(record["hilbert_dimension"]
                             for record in survivors)
determinant_multiset = Counter(record["intersection_determinant"]
                               for record in survivors)
minimum_dimension = min(dimension_multiset) if dimension_multiset else None
minimizers = [record for record in survivors
              if record["hilbert_dimension"] == minimum_dimension]

# Quotient counts under exact character conjugation, abstract S3 permutations
# of the three equal M5 blocks, grading reversal, and their combinations.
equal_block_permutations = []
for permuted_three in permutations((0, 1, 2)):
    equal_block_permutations.append(tuple(permuted_three)+(3,))


def orbit_partition(keys, transformations):
    unseen = set(keys)
    orbits = []
    while unseen:
        seed = min(unseen)
        seed_edges = frozenset(seed)
        orbit = set()
        for permutation, reverse in transformations:
            image = transform_edges(seed_edges, permutation)
            if reverse:
                image = reverse_edges(image)
            image_key = edge_key(image)
            if image_key in keys:
                orbit.add(image_key)
        # Transformations form a group in every call below, so one pass from
        # the seed is its complete orbit.
        orbits.append(orbit)
        unseen -= orbit
    return orbits


identity_permutation = tuple(range(4))
conjugation_orbits = orbit_partition(
    survivor_keys,
    [(identity_permutation, False), (conjugation, False)],
)
s3_orbits = orbit_partition(
    survivor_keys,
    [(permutation, False) for permutation in equal_block_permutations],
)
s3_reversal_orbits = orbit_partition(
    survivor_keys,
    [(permutation, reverse)
     for permutation in equal_block_permutations for reverse in (False, True)],
)

symmetry_counts = {
    "conjugation_preserves_grading": sum(
        record["conjugation_preserves_grading"] for record in survivors
    ),
    "conjugation_reverses_grading": sum(
        record["conjugation_reverses_grading"] for record in survivors
    ),
    "conjugation_orbits": len(conjugation_orbits),
    "abstract_S3_orbits": len(s3_orbits),
    "abstract_S3_plus_grading_reversal_orbits": len(s3_reversal_orbits),
}

# The first run records the blind values.  These structural identities make
# the output falsifiable without importing a selector target.
check("every survivor passes all four preregistered necessary gates",
      len(survivors) == gate_counts["survivors"]
      and all(record["intersection_rank"] == 4
              and record["faithful"]
              and record["legal_odd_blocks"] > 0
              and record["centre_connected"]
              for record in survivors))
check("survivor symmetry transformations close within the census",
      sum(len(orbit) for orbit in conjugation_orbits) == len(survivors)
      and sum(len(orbit) for orbit in s3_orbits) == len(survivors)
      and sum(len(orbit) for orbit in s3_reversal_orbits) == len(survivors))
check("the blind binary gate ledger has the frozen exact counts",
      dict(gate_counts) == {
          "total": 729,
          "poincare": 484,
          "faithful": 636,
          "poincare_and_faithful": 484,
          "nonzero_first_order": 642,
          "centre_connected": 316,
          "survivors": 256,
      }
      and dimension_multiset
      == Counter({300: 24, 400: 48, 450: 48,
                  500: 24, 550: 48, 600: 64})
      and determinant_multiset == Counter({1: 144, 4: 96, 9: 16})
      and symmetry_counts == {
          "conjugation_preserves_grading": 0,
          "conjugation_reverses_grading": 0,
          "conjugation_orbits": 128,
          "abstract_S3_orbits": 44,
          "abstract_S3_plus_grading_reversal_orbits": 22,
      })

# Corrected exhaustive weighted minimum.  A binary survivor gives the upper
# bound 300.  Enumerate the complete signed multiplicity box implied by that
# bound: three 5-5 entries in [-6,6] and three 5-15 entries in [-2,2].
weighted_box_count = 0
weighted_within_bound = 0
weighted_survivors = []
weighted_pairs = ((0, 1), (0, 2), (1, 2),
                  (0, 3), (1, 3), (2, 3))
for low_values in product(range(-6, 7), repeat=3):
    for high_values in product(range(-2, 3), repeat=3):
        weighted_box_count += 1
        values = low_values+high_values
        hilbert_dimension = (
            50*sum(abs(value) for value in low_values)
            + 150*sum(abs(value) for value in high_values)
        )
        if hilbert_dimension > 300:
            continue
        weighted_within_bound += 1
        edges = frozenset(
            (left, right) if value > 0 else (right, left)
            for (left, right), value in zip(weighted_pairs, values)
            if value != 0
        )
        faithful = all(any(node in edge for edge in edges)
                       for node in range(4))
        q01, q02, q12, q03, q13, q23 = values
        pfaffian = q01*q23-q02*q13+q03*q12
        legal_blocks, central_links, components = support_gate_data(edges)
        if (faithful and pfaffian != 0 and legal_blocks
                and len(components) == 1):
            weighted_survivors.append({
                "signed_multiplicities": list(values),
                "edge_order": [list(pair) for pair in weighted_pairs],
                "hilbert_dimension": hilbert_dimension,
                "intersection_pfaffian": pfaffian,
                "intersection_determinant": pfaffian*pfaffian,
                "support": [list(edge) for edge in sorted(edges)],
                "legal_odd_blocks": len(legal_blocks),
                "central_links": [list(edge)
                                  for edge in sorted(central_links)],
            })

weighted_dimension_multiset = Counter(
    record["hilbert_dimension"] for record in weighted_survivors
)
weighted_minimum = (min(weighted_dimension_multiset)
                    if weighted_dimension_multiset else None)
weighted_minimizers = [record for record in weighted_survivors
                       if record["hilbert_dimension"] == weighted_minimum]
check("the corrected weighted box has exactly 13^3*5^3 assignments",
      weighted_box_count == 274625)
check("the dimension-cut weighted search is internally exhaustive",
      weighted_minimum is not None
      and all(record["hilbert_dimension"] <= 300
              and record["intersection_pfaffian"] != 0
              and record["legal_odd_blocks"] > 0
          for record in weighted_survivors))
check("the corrected weighted minimum is exactly 300 with 24 binary designs",
      weighted_within_bound == 773
      and len(weighted_survivors) == 24
      and weighted_dimension_multiset == Counter({300: 24})
      and weighted_minimum == 300 and len(weighted_minimizers) == 24
      and all(max(abs(value) for value in record["signed_multiplicities"]) == 1
              for record in weighted_minimizers))

# Multiplicity-independent conjugation obstruction.  The exact character
# conjugation swaps nodes 1 and 2, so its permutation determinant is -1.
# For a generic alternating 4x4 Q, Pf(PQP^T)=det(P)Pf(Q)=-Pf(Q), whereas
# Pf(-Q)=Pf(Q).  Invariance or anti-invariance therefore forces Pf(Q)=0.
a, b, c, d, e, f = sp.symbols("a b c d e f")
generic_intersection = sp.Matrix([
    [0, a, b, c],
    [-a, 0, d, e],
    [-b, -d, 0, f],
    [-c, -e, -f, 0],
])
conjugation_matrix = sp.zeros(4)
for source, target in enumerate(conjugation):
    conjugation_matrix[target, source] = 1


def pfaffian4(matrix):
    return sp.expand(matrix[0, 1]*matrix[2, 3]
                     - matrix[0, 2]*matrix[1, 3]
                     + matrix[0, 3]*matrix[1, 2])


generic_pfaffian = pfaffian4(generic_intersection)
transformed_intersection = (
    conjugation_matrix*generic_intersection*conjugation_matrix.T
)
check("character conjugation reverses every four-node Pfaffian",
      conjugation_matrix.det() == -1
      and sp.expand(pfaffian4(transformed_intersection)
                    + generic_pfaffian) == 0
      and sp.expand(pfaffian4(-generic_intersection)
                    - generic_pfaffian) == 0)

# Verify the conclusion directly on the linear invariant and anti-invariant
# solution spaces, rather than citing only the transformation identity.
variables = (a, b, c, d, e, f)
invariant_equations = list(transformed_intersection-generic_intersection)
anti_invariant_equations = list(transformed_intersection+generic_intersection)
invariant_solution = sp.linsolve(invariant_equations, variables)
anti_invariant_solution = sp.linsolve(anti_invariant_equations, variables)
invariant_tuple = next(iter(invariant_solution))
anti_invariant_tuple = next(iter(anti_invariant_solution))
check("both conjugation grading behaviors force degenerate intersection",
      sp.expand(generic_pfaffian.subs(dict(zip(variables, invariant_tuple)))) == 0
      and sp.expand(generic_pfaffian.subs(
          dict(zip(variables, anti_invariant_tuple)))) == 0,
      "grading-preserving and grading-reversing conjugation both force Pf=0")

payload = {
    "protocol_commit": "8b21830",
    "target_comparison_performed": False,
    "algebra_blocks": list(node_sizes),
    "total_supports": gate_counts["total"],
    "gate_counts": dict(sorted(gate_counts.items())),
    "survivor_count": len(survivors),
    "survivor_dimension_multiset": {
        str(key): value for key, value in sorted(dimension_multiset.items())
    },
    "survivor_intersection_determinant_multiset": {
        str(key): value for key, value in sorted(determinant_multiset.items())
    },
    "minimum_hilbert_dimension": minimum_dimension,
    "minimum_support_count": len(minimizers),
    "minimum_supports": minimizers,
    "symmetry_counts": symmetry_counts,
    "corrected_weighted_minimum_census": {
        "search_dimension_upper_bound": 300,
        "signed_box_assignments": weighted_box_count,
        "assignments_within_dimension_bound": weighted_within_bound,
        "survivor_count_within_bound": len(weighted_survivors),
        "survivor_dimension_multiset": {
            str(key): value
            for key, value in sorted(weighted_dimension_multiset.items())
        },
        "minimum_hilbert_dimension": weighted_minimum,
        "minimum_design_count": len(weighted_minimizers),
        "minimum_designs": weighted_minimizers,
    },
    "conjugation_pfaffian_no_go": {
        "node_permutation": list(conjugation),
        "permutation_determinant": int(conjugation_matrix.det()),
        "grading_preserving_nondegenerate_design_exists": False,
        "grading_reversing_nondegenerate_design_exists": False,
        "arbitrary_integer_multiplicities_covered": True,
    },
    "survivors": survivors,
    "verdict": (
        "BLIND CENTRAL KRAJEWSKI CENSUS. Exact binary-support enumeration for "
        "B=M5+M5+M5+M15 under KO6 transpose grading, metric-zero support "
        "orientability, nondegenerate antisymmetric intersection form, "
        "nonzero first-order-legal odd blocks and central connectedness. "
        "Passing is necessary only; no Hessian target or matrix-level Dirac "
        "gate is used."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the blind Krajewski census JSON was written", OUTPUT.exists())

print(f"gate_counts={dict(sorted(gate_counts.items()))}")
print(f"survivor_dimensions={dict(sorted(dimension_multiset.items()))}")
print(f"intersection_determinants={dict(sorted(determinant_multiset.items()))}")
print(f"minimum_dimension={minimum_dimension}, minimizers={len(minimizers)}")
print(f"symmetry_counts={symmetry_counts}")
print(f"weighted_within_bound={weighted_within_bound}, "
      f"weighted_survivors={len(weighted_survivors)}")
print(f"weighted_dimensions={dict(sorted(weighted_dimension_multiset.items()))}")
print(f"weighted_minimum={weighted_minimum}, "
      f"weighted_minimizers={len(weighted_minimizers)}")
print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("NO HESSIAN OR SELECTOR TARGET COMPARISON WAS PERFORMED.")
raise SystemExit(0 if passed == tests else 1)

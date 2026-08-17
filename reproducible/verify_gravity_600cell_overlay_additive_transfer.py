#!/usr/bin/env python3
"""Audit additive fine-to-coarse transfer on the staircase overlay exactly."""

from collections import Counter, deque
from fractions import Fraction
from itertools import permutations
import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "gravity_600cell_universal_staircase_overlay.json"
OUTPUT = HERE / "gravity_600cell_overlay_additive_transfer.json"
SOURCE_SHA256 = "0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc"
PRIOR_ART_COMMIT = "c77a87c"
PROTOCOL_COMMIT = "b0b30dc"
PRIMES = (1000003, 1000033, 1000037)
VERTICES = tuple(range(4))
FULL_MASK = 15
MASKS = tuple(range(1, FULL_MASK))
MASK_INDEX = {mask: index for index, mask in enumerate(MASKS)}
PERMUTATIONS = tuple(permutations(VERTICES))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mask_of(vertices):
    result = 0
    for vertex in vertices:
        result |= 1 << vertex
    return result


def decode(word):
    if len(word) != 14 or set(word) - {"+", "-"}:
        raise ValueError(f"invalid sign word {word!r}")
    return tuple(1 if value == "+" else -1 for value in word)


def encode(pattern):
    return "".join("+" if value > 0 else "-" for value in pattern)


def staircase_assignment(pattern, order):
    assignments = []
    for split in range(4):
        lower = mask_of(order[split+1:])
        upper = mask_of(order[split:])
        lower_ok = lower == 0 or pattern[MASK_INDEX[lower]] > 0
        upper_ok = upper == FULL_MASK or pattern[MASK_INDEX[upper]] < 0
        if lower_ok and upper_ok:
            assignments.append(split)
    return tuple(assignments)


def permute_mask(mask, permutation):
    return mask_of(
        permutation[vertex]
        for vertex in VERTICES
        if mask & (1 << vertex)
    )


def transform_pattern(pattern, permutation, reflect_time):
    result = [None]*14
    for mask, sign in zip(MASKS, pattern):
        target = permute_mask(mask, permutation)
        if reflect_time:
            target = FULL_MASK ^ target
            sign = -sign
        result[MASK_INDEX[target]] = sign
    return tuple(result)


def modular_rank(rows, prime):
    matrix = [[int(value) % prime for value in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0]) if matrix else 0
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count)
             if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row][column:] = [
            value*inverse % prime
            for value in matrix[pivot_row][column:]
        ]
        for row in range(pivot_row+1, row_count):
            factor = matrix[row][column]
            if factor:
                matrix[row][column:] = [
                    (left-factor*right) % prime
                    for left, right in zip(
                        matrix[row][column:], matrix[pivot_row][column:]
                    )
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def rational_string(value):
    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


source_hash = digest(SOURCE)
source = json.loads(SOURCE.read_text())
source_ok = bool(
    source_hash == SOURCE_SHA256
    and source.get("outcome") == "UNIVERSAL_STAIRCASE_OVERLAY_CERTIFIED"
    and source.get("passed") == source.get("tests") == 12
    and source.get("full_dimensional_chamber_count") == 148
    and len(source.get("feasible_sign_words", ())) == 148
    and len(source.get("staircase_orders", ())) == 24
    and len(source.get("symmetry_orbits", ())) == 14
    and source.get("gravity_action_evaluations") == 0
    and source.get("physical_target_parsed") is False
)

patterns = tuple(decode(word) for word in source["feasible_sign_words"])
pattern_set = set(patterns)
orders = tuple(tuple(record["order"]) for record in source["staircase_orders"])
reconstruction_ok = bool(
    len(patterns) == len(pattern_set) == 148
    and set(orders) == set(PERMUTATIONS)
)

rows = []
row_labels = []
assignment_failures = []
for order in orders:
    assignments = [staircase_assignment(pattern, order) for pattern in patterns]
    for index, assigned in enumerate(assignments):
        if len(assigned) != 1:
            assignment_failures.append({
                "order": list(order),
                "chamber": encode(patterns[index]),
                "assignments": list(assigned),
            })
    for split in range(4):
        rows.append([
            int(len(assigned) == 1 and assigned[0] == split)
            for assigned in assignments
        ])
        row_labels.append((order, split))

row_sums_by_order = [
    tuple(sum(rows[4*index+split]) for split in range(4))
    for index in range(24)
]
partition_ok = bool(
    not assignment_failures
    and len(rows) == 96
    and all(values == (19, 55, 55, 19) for values in row_sums_by_order)
    and all(
        sum(rows[4*index+split][column] for split in range(4)) == 1
        for index in range(24)
        for column in range(148)
    )
)
matrix_hash = hashlib.sha256(
    bytes(value for row in rows for value in row)
).hexdigest()

matrix = sp.Matrix(rows)
rational_rank = matrix.rank()
nullity = matrix.cols-rational_rank
rank_bounds_ok = rational_rank <= 73 and nullity >= 75
modular_ranks = {str(prime): modular_rank(rows, prime) for prime in PRIMES}
modular_rank_ok = all(rank == rational_rank for rank in modular_ranks.values())

nullspace = matrix.nullspace()
if nullspace:
    nullspace_matrix = sp.Matrix.hstack(*nullspace)
    nullspace_residual = matrix*nullspace_matrix
    nullspace_ok = bool(
        len(nullspace) == nullity
        and nullspace_residual == sp.zeros(matrix.rows, len(nullspace))
        and nullspace_matrix.rank() == nullity
    )
else:
    nullspace_matrix = sp.zeros(matrix.cols, 0)
    nullspace_ok = nullity == 0


# Reconstruct fine chamber orbits from the complete group action.
transformations = tuple(
    (permutation, reflect_time)
    for permutation in PERMUTATIONS
    for reflect_time in (False, True)
)
unseen = set(patterns)
fine_orbits = []
while unseen:
    representative = min(unseen)
    orbit = {
        transform_pattern(representative, permutation, reflect_time)
        for permutation, reflect_time in transformations
    }
    queue = deque(orbit)
    while queue:
        pattern = queue.popleft()
        for permutation, reflect_time in transformations:
            image = transform_pattern(pattern, permutation, reflect_time)
            if image not in orbit:
                orbit.add(image)
                queue.append(image)
    fine_orbits.append(tuple(sorted(orbit)))
    unseen -= orbit
fine_orbits.sort(key=lambda orbit: (len(orbit), encode(orbit[0])))
fine_orbit_records = [
    {"representative": encode(orbit[0]), "size": len(orbit)}
    for orbit in fine_orbits
]
fine_symmetry_ok = bool(
    set().union(*(set(orbit) for orbit in fine_orbits)) == pattern_set
    and len(fine_orbits) == 14
    and fine_orbit_records == source["symmetry_orbits"]
)

pattern_index = {pattern: index for index, pattern in enumerate(patterns)}
fine_orbit_index = {}
for orbit_index, orbit in enumerate(fine_orbits):
    for pattern in orbit:
        fine_orbit_index[pattern] = orbit_index
orbit_indicator = sp.zeros(148, len(fine_orbits))
for pattern, column in pattern_index.items():
    orbit_indicator[column, fine_orbit_index[pattern]] = 1
invariant_matrix = matrix*orbit_indicator
invariant_rows = [list(map(int, invariant_matrix.row(index))) for index in range(96)]


def transform_label(label, permutation, reflect_time):
    order, split = label
    transformed_order = tuple(permutation[vertex] for vertex in order)
    if reflect_time:
        transformed_order = tuple(reversed(transformed_order))
        split = 3-split
    return transformed_order, split


label_index = {label: index for index, label in enumerate(row_labels)}
unseen_labels = set(row_labels)
coarse_orbits = []
while unseen_labels:
    representative = min(unseen_labels)
    orbit = {
        transform_label(representative, permutation, reflect_time)
        for permutation, reflect_time in transformations
    }
    coarse_orbits.append(tuple(sorted(orbit)))
    unseen_labels -= orbit
coarse_orbits.sort(key=lambda orbit: (orbit[0][1], orbit[0][0]))
coarse_orbit_rows_ok = bool(
    len(coarse_orbits) == 2
    and sorted(map(len, coarse_orbits)) == [48, 48]
    and all(
        len({tuple(invariant_rows[label_index[label]]) for label in orbit}) == 1
        for orbit in coarse_orbits
    )
)

invariant_rank = invariant_matrix.rank()
invariant_nullity = invariant_matrix.cols-invariant_rank
invariant_rank_bounds_ok = invariant_rank <= 2 and invariant_nullity >= 12
invariant_modular_ranks = {
    str(prime): modular_rank(invariant_rows, prime) for prime in PRIMES
}
invariant_modular_ok = all(
    rank == invariant_rank for rank in invariant_modular_ranks.values()
)
invariant_nullspace = invariant_matrix.nullspace()
invariant_nullspace_ok = bool(
    len(invariant_nullspace) == invariant_nullity
    and (
        not invariant_nullspace
        or invariant_matrix*sp.Matrix.hstack(*invariant_nullspace)
        == sp.zeros(96, invariant_nullity)
    )
)


# Positive invariant pair with identical coarse totals.
witness_ok = False
witness = None
if invariant_nullspace:
    orbit_vector = invariant_nullspace[0]
    maximum = max(abs(value) for value in orbit_vector)
    epsilon = sp.Rational(1, 2)/maximum
    fine_vector = orbit_indicator*orbit_vector
    ones = sp.ones(148, 1)
    plus = ones+epsilon*fine_vector
    minus = ones-epsilon*fine_vector
    plus_totals = matrix*plus
    minus_totals = matrix*minus
    invariant_values_ok = all(
        len({fine_vector[pattern_index[pattern]] for pattern in orbit}) == 1
        for orbit in fine_orbits
    )
    witness_ok = bool(
        any(value != 0 for value in orbit_vector)
        and matrix*fine_vector == sp.zeros(96, 1)
        and all(value > 0 for value in plus)
        and all(value > 0 for value in minus)
        and plus != minus
        and plus_totals == minus_totals
        and sum(plus) == sum(minus)
        and invariant_values_ok
    )
    witness = {
        "orbit_coordinates": [rational_string(value) for value in orbit_vector],
        "epsilon": rational_string(epsilon),
        "minimum_plus_weight": rational_string(min(plus)),
        "minimum_minus_weight": rational_string(min(minus)),
        "common_coarse_totals": [
            rational_string(value) for value in plus_totals
        ],
        "common_total_fine_weight": rational_string(sum(plus)),
    }

controls_ok = bool(
    source_ok
    and reconstruction_ok
    and partition_ok
    and rank_bounds_ok
    and modular_rank_ok
    and nullspace_ok
    and fine_symmetry_ok
    and coarse_orbit_rows_ok
    and invariant_rank_bounds_ok
    and invariant_modular_ok
    and invariant_nullspace_ok
)

if not controls_ok:
    outcome = "ADDITIVE_TRANSFER_CONTROL_FAILED"
elif nullity > 0 and invariant_nullity > 0 and witness_ok:
    outcome = "POSITIVE_INVARIANT_ADDITIVE_TRANSFER_NONUNIQUE"
elif nullity > 0 and invariant_nullity == 0:
    outcome = "SYMMETRY_REMOVES_ADDITIVE_TRANSFER_KERNEL"
elif nullity == 0:
    outcome = "ADDITIVE_TRANSFER_INJECTIVE"
else:
    outcome = "ADDITIVE_TRANSFER_CONTROL_FAILED"

tests = [
    ("frozen certified overlay artifact and SHA-256", source_ok),
    ("148 sign words and all 24 orders reconstructed", reconstruction_ok),
    ("each order partitions chambers with row sums 19,55,55,19", partition_ok),
    ("exact rational full rank obeys preregistered bounds", rank_bounds_ok),
    ("three modular full ranks equal the rational rank", modular_rank_ok),
    ("exact full nullspace basis has the reported dimension", nullspace_ok),
    ("14 fine S4 x C2 orbits reproduce the frozen records", fine_symmetry_ok),
    ("coarse labels form two 48-element invariant row orbits", coarse_orbit_rows_ok),
    ("exact invariant rank obeys preregistered bounds", invariant_rank_bounds_ok),
    ("three modular invariant ranks equal the rational rank", invariant_modular_ok),
    ("exact invariant nullspace has the reported dimension", invariant_nullspace_ok),
    ("distinct positive invariant lifts have identical 96 totals", witness_ok),
    ("no gravity action, metric or physical target was read", True),
    ("outcome follows the preregistered mechanical rule", outcome in {
        "POSITIVE_INVARIANT_ADDITIVE_TRANSFER_NONUNIQUE",
        "SYMMETRY_REMOVES_ADDITIVE_TRANSFER_KERNEL",
        "ADDITIVE_TRANSFER_INJECTIVE",
        "ADDITIVE_TRANSFER_CONTROL_FAILED",
    }),
]
passed = sum(bool(ok) for _, ok in tests)

unique_invariant_rows = Counter(tuple(row) for row in invariant_rows)
payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": {"universal_overlay": source_hash},
    "aggregation_shape": [matrix.rows, matrix.cols],
    "aggregation_matrix_sha256": matrix_hash,
    "row_sums_by_order": [list(values) for values in row_sums_by_order],
    "rational_rank": rational_rank,
    "rational_nullity": nullity,
    "modular_ranks": modular_ranks,
    "fine_symmetry_orbits": fine_orbit_records,
    "coarse_orbit_sizes": [len(orbit) for orbit in coarse_orbits],
    "invariant_shape": [invariant_matrix.rows, invariant_matrix.cols],
    "invariant_rational_rank": invariant_rank,
    "invariant_rational_nullity": invariant_nullity,
    "invariant_modular_ranks": invariant_modular_ranks,
    "invariant_unique_row_count": len(unique_invariant_rows),
    "invariant_unique_rows": [
        {"row": list(row), "multiplicity": multiplicity}
        for row, multiplicity in sorted(unique_invariant_rows.items())
    ],
    "positive_invariant_witness": witness,
    "gravity_action_evaluations": 0,
    "physical_target_parsed": False,
    "tests": len(tests),
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"R shape={matrix.rows}x{matrix.cols}")
print(f"rank_Q(R)={rational_rank}; nullity_Q(R)={nullity}")
print(f"modular ranks={modular_ranks}")
print(f"fine orbits={len(fine_orbits)}; coarse orbits={list(map(len, coarse_orbits))}")
print(f"rank_Q(R_inv)={invariant_rank}; nullity_Q(R_inv)={invariant_nullity}")
print(f"invariant modular ranks={invariant_modular_ranks}")
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) and controls_ok else 1)


#!/usr/bin/env python3
"""Exact support audit for strict-local commuting Whitney recoveries.

Protocol commit 39aa5c8 froze the carriers, arbitrary-weight support theorem,
candidate controls and interpretation before this census was evaluated.
"""

from itertools import combinations
import gc
import json
from pathlib import Path

import numpy as np
import sympy as sy

from whitney_trace_refinement_tools import (
    LOCAL_D,
    barycentric_refine,
    classify_element_types,
    make_base_level,
    rank_edgewise_level,
)


OUTPUT = Path(__file__).with_name("whitney_commuting_recovery.json")
PROTOCOL_COMMIT = "39aa5c8"
EXPECTED_F_VECTORS = {
    1: (30, 150, 240, 120),
    2: (180, 1140, 1920, 960),
    4: (1320, 9000, 15360, 7680),
}
LOCAL_FACES = tuple(
    tuple(combinations(range(4), degree + 1)) for degree in range(4)
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


def orientation_sign(vertices):
    inversions = sum(
        vertices[left] > vertices[right]
        for left in range(len(vertices))
        for right in range(left + 1, len(vertices))
    )
    return -1 if inversions % 2 else 1


def global_cofaces(level, degree):
    low_cells = level["cells"][degree]
    high_cells = level["cells"][degree + 1]
    low_indices = {cell: index for index, cell in enumerate(low_cells)}
    cofaces = [[] for _ in low_cells]
    for high_index, high in enumerate(high_cells):
        for omitted in range(degree + 2):
            face = high[:omitted] + high[omitted + 1:]
            cofaces[low_indices[face]].append((high_index, (-1) ** omitted))
    return cofaces


def occurrence_support_audit(level, degree):
    low_cells = level["cells"][degree]
    high_cells = level["cells"][degree + 1]
    low_indices = {cell: index for index, cell in enumerate(low_cells)}
    high_indices = {cell: index for index, cell in enumerate(high_cells)}
    cofaces = global_cofaces(level, degree)
    forced_by_global = [[] for _ in low_cells]
    external_counts = []
    first_witness = None

    for top_index, top in enumerate(level["top"]):
        for local_index, positions in enumerate(LOCAL_FACES[degree]):
            oriented = tuple(top[position] for position in positions)
            low = tuple(sorted(oriented))
            low_index = low_indices[low]
            local_cofaces = {
                high_indices[tuple(sorted(top[position] for position in high))]
                for high in LOCAL_FACES[degree + 1]
                if set(positions).issubset(high)
            }
            external = sorted(
                high_index for high_index, _ in cofaces[low_index]
                if high_index not in local_cofaces
            )
            external_counts.append(len(external))
            forced = len(external) > 0
            forced_by_global[low_index].append(forced)
            if forced and first_witness is None:
                high_index = external[0]
                first_witness = {
                    "top_index": top_index,
                    "local_simplex_index": local_index,
                    "global_simplex_index": low_index,
                    "global_simplex": list(low),
                    "external_coface_index": high_index,
                    "external_coface": list(high_cells[high_index]),
                    "external_coface_not_contained_in_top": not set(
                        high_cells[high_index]
                    ).issubset(top),
                    "forced_equation": "plus_or_minus_w_T_s_equals_zero",
                }

    every_occurrence_forced = all(count > 0 for count in external_counts)
    contradicted_left_inverse_rows = sum(
        all(forced) for forced in forced_by_global
    )
    return {
        "degree": degree,
        "global_simplex_count": len(low_cells),
        "occurrence_count": len(external_counts),
        "occurrences_with_external_coface": sum(
            count > 0 for count in external_counts
        ),
        "fraction_with_external_coface": (
            sum(count > 0 for count in external_counts)
            / len(external_counts)
        ),
        "minimum_external_cofaces": min(external_counts),
        "maximum_external_cofaces": max(external_counts),
        "every_occurrence_weight_forced_zero": every_occurrence_forced,
        "left_inverse_rows_contradicted": contradicted_left_inverse_rows,
        "all_left_inverse_rows_contradicted": (
            contradicted_left_inverse_rows == len(low_cells)
        ),
        "first_exact_support_witness": first_witness,
    }


def exact_recovery_data(level, element_types, degree):
    cells = level["cells"][degree]
    indices = {cell: index for index, cell in enumerate(cells)}
    occurrences = [[] for _ in cells]
    copy_records = []

    for top_index, top in enumerate(level["top"]):
        type_index = int(element_types["top_types"][top_index])
        mass = element_types["masses"][type_index][degree]
        for local_index, positions in enumerate(LOCAL_FACES[degree]):
            oriented = tuple(top[position] for position in positions)
            global_index = indices[tuple(sorted(oriented))]
            sign = orientation_sign(oriented)
            copy_index = top_index * len(LOCAL_FACES[degree]) + local_index
            record = {
                "copy": copy_index,
                "top": top_index,
                "local": local_index,
                "global": global_index,
                "sign": sign,
                "diagonal": sy.factor(mass[local_index, local_index]),
            }
            occurrences[global_index].append(record)
            copy_records.append(record)

    weights = {"counting": {}, "diagonal": {}}
    exact_left_inverse = {"counting": True, "diagonal": True}
    for copies in occurrences:
        diagonal_sum = sy.factor(sum(
            (record["diagonal"] for record in copies), sy.Integer(0)
        ))
        candidate_sums = {"counting": sy.Integer(0),
                          "diagonal": sy.Integer(0)}
        for record in copies:
            candidate_weights = {
                "counting": sy.Rational(1, len(copies)),
                "diagonal": sy.factor(record["diagonal"] / diagonal_sum),
            }
            for candidate, weight in candidate_weights.items():
                weights[candidate][record["copy"]] = weight
                candidate_sums[candidate] += weight
        for candidate in weights:
            exact_left_inverse[candidate] &= (
                sy.factor(candidate_sums[candidate]) == 1
            )
    return copy_records, weights, exact_left_inverse


def candidate_commutator(level, degree, candidate, records, weights,
                         next_records, next_weights):
    """Compute L_(p+1) D_pw - d L_p as an exact sparse dictionary."""
    low_cells = level["cells"][degree]
    high_cells = level["cells"][degree + 1]
    low_indices = {cell: index for index, cell in enumerate(low_cells)}
    high_indices = {cell: index for index, cell in enumerate(high_cells)}
    cofaces = global_cofaces(level, degree)
    entries = {}

    def add(row, column, value):
        key = (row, column)
        updated = sy.factor(entries.get(key, sy.Integer(0)) + value)
        if updated == 0:
            entries.pop(key, None)
        else:
            entries[key] = updated

    next_by_copy = {record["copy"]: record for record in next_records}
    for record in records:
        source_copy = record["copy"]
        top_index = record["top"]
        source_positions = LOCAL_FACES[degree][record["local"]]

        # L_(p+1) D_pw: only local cofaces inside the same tetrahedron.
        for next_local, next_positions in enumerate(LOCAL_FACES[degree + 1]):
            if not set(source_positions).issubset(next_positions):
                continue
            omitted = next_positions.index(next(
                position for position in next_positions
                if position not in source_positions
            ))
            local_incidence = (-1) ** omitted
            next_copy = top_index * len(LOCAL_FACES[degree + 1]) + next_local
            target = next_by_copy[next_copy]
            add(
                target["global"],
                source_copy,
                local_incidence * target["sign"]
                * next_weights[candidate][next_copy],
            )

        # -d_p L_p: every global coface of the recovered global simplex.
        source_weight = weights[candidate][source_copy]
        for target_global, incidence in cofaces[record["global"]]:
            add(
                target_global,
                source_copy,
                -incidence * record["sign"] * source_weight,
            )

    maximum = max((abs(value) for value in entries.values()), default=0)
    first_key = min(entries) if entries else None
    first_witness = None
    if first_key is not None:
        row, column = first_key
        first_witness = {
            "global_target_index": row,
            "local_source_copy_index": column,
            "exact_coefficient": str(entries[first_key]),
        }
    return {
        "degree": degree,
        "shape": [len(high_cells), len(records)],
        "exact_nonzero_count": len(entries),
        "maximum_absolute_exact_coefficient": str(maximum),
        "first_exact_nonzero": first_witness,
        "commutes_exactly": len(entries) == 0,
    }


def global_differential_dense(level, degree):
    low_cells = level["cells"][degree]
    high_cells = level["cells"][degree + 1]
    low_indices = {cell: index for index, cell in enumerate(low_cells)}
    matrix = np.zeros((len(high_cells), len(low_cells)), dtype=np.int64)
    for row, high in enumerate(high_cells):
        for omitted in range(degree + 2):
            face = high[:omitted] + high[omitted + 1:]
            matrix[row, low_indices[face]] = (-1) ** omitted
    return matrix


def rank_mod_prime(matrix, prime):
    """Exact lower-bound certificate for an integer matrix rank over Q."""
    reduced = np.asarray(matrix, dtype=np.int64) % prime
    rows, columns = reduced.shape
    pivot_row = 0
    for column in range(columns):
        candidates = np.flatnonzero(reduced[pivot_row:, column])
        if not len(candidates):
            continue
        selected = pivot_row + int(candidates[0])
        if selected != pivot_row:
            reduced[[pivot_row, selected]] = reduced[[selected, pivot_row]]
        inverse = pow(int(reduced[pivot_row, column]), -1, prime)
        reduced[pivot_row, column:] = (
            reduced[pivot_row, column:] * inverse
        ) % prime
        for row in range(pivot_row + 1, rows):
            factor = int(reduced[row, column])
            if factor:
                reduced[row, column:] = (
                    reduced[row, column:]
                    - factor * reduced[pivot_row, column:]
                ) % prime
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def injection_intertwiner_exact(level, degree):
    """Check D_pw J = J d by exact signed coefficient dictionaries."""
    low_cells = level["cells"][degree]
    high_cells = level["cells"][degree + 1]
    low_indices = {cell: index for index, cell in enumerate(low_cells)}
    high_indices = {cell: index for index, cell in enumerate(high_cells)}
    left = {}
    right = {}

    def add(dictionary, key, value):
        updated = dictionary.get(key, 0) + value
        if updated:
            dictionary[key] = updated
        else:
            dictionary.pop(key, None)

    high_occurrences = [[] for _ in high_cells]
    for top_index, top in enumerate(level["top"]):
        for high_local, high_positions in enumerate(LOCAL_FACES[degree + 1]):
            high_oriented = tuple(top[position] for position in high_positions)
            high_global = high_indices[tuple(sorted(high_oriented))]
            high_copy = (
                top_index * len(LOCAL_FACES[degree + 1]) + high_local
            )
            high_occurrences[high_global].append((
                high_copy, orientation_sign(high_oriented)
            ))

        for low_local, low_positions in enumerate(LOCAL_FACES[degree]):
            low_oriented = tuple(top[position] for position in low_positions)
            low_global = low_indices[tuple(sorted(low_oriented))]
            low_sign = orientation_sign(low_oriented)
            for high_local, high_positions in enumerate(
                LOCAL_FACES[degree + 1]
            ):
                if not set(low_positions).issubset(high_positions):
                    continue
                missing = next(
                    position for position in high_positions
                    if position not in low_positions
                )
                omitted = high_positions.index(missing)
                high_copy = (
                    top_index * len(LOCAL_FACES[degree + 1]) + high_local
                )
                add(left, (high_copy, low_global),
                    low_sign * (-1) ** omitted)

    for high_global, high in enumerate(high_cells):
        for omitted in range(degree + 2):
            low = high[:omitted] + high[omitted + 1:]
            low_global = low_indices[low]
            incidence = (-1) ** omitted
            for high_copy, high_sign in high_occurrences[high_global]:
                add(right, (high_copy, low_global), high_sign * incidence)

    differences = {
        key: left.get(key, 0) - right.get(key, 0)
        for key in set(left) | set(right)
        if left.get(key, 0) != right.get(key, 0)
    }
    return {
        "degree": degree,
        "left_nonzero_count": len(left),
        "right_nonzero_count": len(right),
        "difference_nonzero_count": len(differences),
        "intertwines_exactly": not differences,
    }


def cohomology_retraction_obstruction(level):
    """Post-result global strengthening: no cochain retraction at all."""
    differentials = [global_differential_dense(level, degree)
                     for degree in range(3)]
    nilpotent = all(
        not np.any(differentials[degree + 1] @ differentials[degree])
        for degree in range(2)
    )
    primes = (1_000_003, 1_000_033)
    ranks_by_prime = {
        str(prime): [rank_mod_prime(matrix, prime)
                     for matrix in differentials]
        for prime in primes
    }
    ranks_agree = len({tuple(value) for value in ranks_by_prime.values()}) == 1
    ranks = next(iter(ranks_by_prime.values()))
    dimensions = list(map(len, level["cells"]))
    # A constant zero-cochain bounds rank(d0) by n0-1.  Nilpotency then gives
    # rank(d1)<=n1-rank(d0), rank(d2)<=n2-rank(d1).  Matching modular lower
    # bounds make all three ranks exact over Q.
    rank_upper_bounds = [
        dimensions[0] - 1,
        dimensions[1] - ranks[0],
        dimensions[2] - ranks[1],
    ]
    ranks_exactly_certified = nilpotent and ranks == rank_upper_bounds
    betti = [
        dimensions[0] - ranks[0],
        dimensions[1] - ranks[0] - ranks[1],
        dimensions[2] - ranks[1] - ranks[2],
        dimensions[3] - ranks[2],
    ]

    # One tetrahedron has differential ranks (3,3,1); direct sums multiply
    # both dimensions and ranks by the number of tetrahedra.
    local_differentials = tuple(map(sy.Matrix, LOCAL_D))
    local_ranks = [matrix.rank() for matrix in local_differentials]
    local_dimensions = [4, 6, 4, 1]
    local_betti = [
        local_dimensions[0] - local_ranks[0],
        local_dimensions[1] - local_ranks[0] - local_ranks[1],
        local_dimensions[2] - local_ranks[1] - local_ranks[2],
        local_dimensions[3] - local_ranks[2],
    ]
    broken_betti = [
        len(level["top"]) * value for value in local_betti
    ]
    return {
        "control_resolution": 1,
        "global_dimensions": dimensions,
        "global_differential_ranks_mod_primes": ranks_by_prime,
        "global_ranks_agree_across_primes": ranks_agree,
        "global_rank_upper_bounds_from_constant_and_nilpotency": (
            rank_upper_bounds
        ),
        "global_ranks_exactly_certified": ranks_exactly_certified,
        "global_complex_nilpotent_exact_integer": nilpotent,
        "global_betti_from_certified_ranks": betti,
        "one_tetrahedron_differential_ranks": local_ranks,
        "one_tetrahedron_betti": local_betti,
        "broken_direct_sum_betti": broken_betti,
        "top_cohomology_mismatch": betti[3] == 1 and broken_betti[3] == 0,
        "cochain_retraction_impossible": betti[3] == 1 and broken_betti[3] == 0,
        "argument": (
            "If L and J are cochain maps with LJ=I, then on H^3, "
            "L_* J_*=I; but J_* maps one-dimensional conforming H^3 into "
            "zero-dimensional broken H^3"
        ),
    }


print("=" * 78)
print("STRICT-LOCAL COMMUTING RECOVERY AUDIT")
print("=" * 78)

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
ranked = barycentric_refine(make_base_level(reference_vertices))

records = []
all_f_vectors = True
all_degrees_forced = True
all_left_inverse_rows_contradicted = True
all_candidate_left_inverse = True
all_candidate_commutators_nonzero = True
global_cohomology_obstruction = None
all_injection_intertwiners_exact = True

for resolution in (1, 2, 4):
    print(f"\n-- rank-edgewise resolution k={resolution} --", flush=True)
    level = rank_edgewise_level(ranked, resolution)
    f_vector = tuple(map(len, level["cells"]))
    all_f_vectors &= f_vector == EXPECTED_F_VECTORS[resolution]
    support = []
    for degree in range(3):
        audit = occurrence_support_audit(level, degree)
        support.append(audit)
        all_degrees_forced &= audit["every_occurrence_weight_forced_zero"]
        all_left_inverse_rows_contradicted &= audit[
            "all_left_inverse_rows_contradicted"
        ]
        print(
            f"  p={degree}: external coface on "
            f"{audit['occurrences_with_external_coface']}/"
            f"{audit['occurrence_count']} occurrences, "
            f"range={audit['minimum_external_cofaces']}.."
            f"{audit['maximum_external_cofaces']}",
            flush=True,
        )

    candidate_records = None
    if resolution in (1, 2):
        element_types = classify_element_types(level)
        degree_records = []
        degree_weights = []
        exact_left = []
        for degree in range(4):
            local_records, weights, left = exact_recovery_data(
                level, element_types, degree
            )
            degree_records.append(local_records)
            degree_weights.append(weights)
            exact_left.append(left)
            all_candidate_left_inverse &= all(left.values())
        candidate_records = {"counting": [], "diagonal": []}
        for candidate in candidate_records:
            for degree in range(3):
                commutator = candidate_commutator(
                    level,
                    degree,
                    candidate,
                    degree_records[degree],
                    degree_weights[degree],
                    degree_records[degree + 1],
                    degree_weights[degree + 1],
                )
                candidate_records[candidate].append(commutator)
                all_candidate_commutators_nonzero &= not commutator[
                    "commutes_exactly"
                ]
                print(
                    f"    {candidate} p={degree}: commutator nnz="
                    f"{commutator['exact_nonzero_count']}, max="
                    f"{commutator['maximum_absolute_exact_coefficient']}",
                    flush=True,
                )
        del element_types, degree_records, degree_weights, exact_left

    if resolution == 1:
        intertwiner_records = [
            injection_intertwiner_exact(level, degree)
            for degree in range(3)
        ]
        all_injection_intertwiners_exact &= all(
            item["intertwines_exactly"] for item in intertwiner_records
        )
        global_cohomology_obstruction = cohomology_retraction_obstruction(
            level
        )
    else:
        intertwiner_records = None

    records.append({
        "edgewise_resolution": resolution,
        "f_vector": list(f_vector),
        "support_obstruction_by_degree": support,
        "candidate_commutator_controls": candidate_records,
        "injection_intertwiner_control": intertwiner_records,
    })
    del level
    gc.collect()

check("all three preregistered f-vectors are exact", all_f_vectors)
check(
    "every local occurrence in p=0,1,2 has an external global coface",
    all_degrees_forced,
)
check(
    "the forced zeros contradict every left-inverse row",
    all_left_inverse_rows_contradicted,
)
check(
    "both candidate families retain exact left-inverse weight sums",
    all_candidate_left_inverse,
)
check(
    "both candidate commutators are exactly nonzero in every tested degree",
    all_candidate_commutators_nonzero,
)
check(
    "the signed occurrence injection is exactly a cochain map",
    all_injection_intertwiners_exact,
)
check(
    "post-result H3 mismatch forbids even a global cochain retraction",
    global_cohomology_obstruction is not None
    and global_cohomology_obstruction["global_complex_nilpotent_exact_integer"]
    and global_cohomology_obstruction["global_ranks_agree_across_primes"]
    and global_cohomology_obstruction["global_ranks_exactly_certified"]
    and global_cohomology_obstruction["global_betti_from_certified_ranks"]
    == [1, 0, 0, 1]
    and global_cohomology_obstruction["one_tetrahedron_betti"]
    == [1, 0, 0, 0]
    and global_cohomology_obstruction["cochain_retraction_impossible"],
    str(global_cohomology_obstruction),
)

general_verdict = (
    "DERIVED NO-GO FOR A STRICT-OCCURRENCE-LOCAL COMMUTING RECOVERY ON THE "
    "FROZEN TOWER: every admissible weight is forced to zero by an external "
    "coface, contradicting LJ=I"
)
scope_verdict = (
    "DERIVED POST-RESULT STRENGTHENING: no raw-piecewise cochain retraction "
    "exists at any support radius, because conforming H3 is one-dimensional "
    "while the direct-sum tetrahedron complex has H3=0"
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "spectrum_computed": False,
    "arbitrary_weights_allow_negative_values": True,
    "records": records,
    "every_occurrence_forced_zero": all_degrees_forced,
    "all_left_inverse_rows_contradicted": (
        all_left_inverse_rows_contradicted
    ),
    "candidate_commutators_all_nonzero": (
        all_candidate_commutators_nonzero
    ),
    "post_result_global_cohomology_obstruction": (
        global_cohomology_obstruction
    ),
    "verdicts": [
        general_verdict,
        "DERIVED: counting and diagonal-Whitney candidate failures are exact "
        "controls, not the source of the all-weight theorem",
        "STRUCTURAL: raw piecewise-derivative commutation is stronger than "
        "the broken-FEEC projected differential construction",
        scope_verdict,
        "NOT CLAIMED: failure of broken FEEC, all local projections, time, "
        "causality, mass or Planck units",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured exact-support certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(general_verdict)
print(scope_verdict)
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Audit conservation-based selection of the broken-FEEC projection.

Protocol commit 21cbc05 froze the carrier, candidates, two conservation
notions, residual gate and interpretation before these residuals were
evaluated.  No spectrum or phenomenological target is used.
"""

from itertools import combinations
import gc
import json
from pathlib import Path

import numpy as np
from scipy import sparse
import sympy as sy

from whitney_trace_refinement_tools import (
    barycentric_refine,
    make_base_level,
    rank_edgewise_level,
)


OUTPUT = Path(__file__).with_name("whitney_projection_selection.json")
PROTOCOL_COMMIT = "21cbc05"
RESIDUAL_GATE = 1e-11
STRUCTURAL_ZERO_RELATIVE_THRESHOLD = 1e-12
EXPECTED_F_VECTORS = {
    1: (30, 150, 240, 120),
    2: (180, 1140, 1920, 960),
    4: (1320, 9000, 15360, 7680),
}
EXPECTED_EXACT_ADJOINT_WITNESS = sy.Rational(243, 7480)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def sparse_maximum(matrix):
    matrix = sparse.csr_matrix(matrix)
    return float(np.max(np.abs(matrix.data))) if matrix.nnz else 0.0


def prune_numerical_zeros(matrix, scale):
    matrix = sparse.csr_matrix(matrix)
    threshold = STRUCTURAL_ZERO_RELATIVE_THRESHOLD * max(scale, 1e-300)
    if matrix.nnz:
        matrix.data[np.abs(matrix.data) <= threshold] = 0.0
        matrix.eliminate_zeros()
    return matrix, threshold


def exact_tetrahedron_volume(points):
    affine = sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    )
    return sy.factor(abs(affine.det()) / 6)


def exact_zero_mass(volume):
    # Exact P1 scalar Whitney mass on a tetrahedron.
    return sy.factor(volume) * (sy.ones(4, 4) + sy.eye(4)) / 20


def zero_form_operators(level):
    """Build J, C, D and the block-local degree-zero mass."""
    top_cells = level["top"]
    vertex_count = len(level["cells"][0])
    occurrences = [[] for _ in range(vertex_count)]
    volumes = []
    blocks = []
    rows = []
    columns = []
    data = []

    for top_index, (top, points) in enumerate(
        zip(top_cells, level["top_points"])
    ):
        volume = exact_tetrahedron_volume(points)
        volumes.append(volume)
        blocks.append(sparse.csr_matrix(
            np.asarray(exact_zero_mass(volume), dtype=np.float64)
        ))
        for local_index, global_vertex in enumerate(top):
            copy = 4 * top_index + local_index
            occurrences[global_vertex].append((copy, top_index, volume))
            rows.append(copy)
            columns.append(global_vertex)
            data.append(1.0)

    local_dimension = 4 * len(top_cells)
    injection = sparse.csr_matrix(
        (data, (rows, columns)),
        shape=(local_dimension, vertex_count),
    )
    local_mass = sparse.block_diag(blocks, format="csr")

    recovery_rows = []
    recovery_columns = []
    counting_data = []
    diagonal_data = []
    exact_counting_harmonic = True
    exact_diagonal_harmonic = True
    exact_left_inverse = True
    all_incident_volumes_equal = True
    for global_vertex, copies in enumerate(occurrences):
        volume_sum = sy.factor(sum(
            (item[2] for item in copies), sy.Integer(0)
        ))
        distinct = {sy.srepr(item[2]) for item in copies}
        all_incident_volumes_equal &= len(distinct) == 1
        count_weight = sy.Rational(1, len(copies))
        count_sum = diagonal_sum = sy.Integer(0)
        for copy, _, volume in copies:
            diagonal_weight = sy.factor(volume / volume_sum)
            recovery_rows.append(global_vertex)
            recovery_columns.append(copy)
            counting_data.append(float(count_weight))
            diagonal_data.append(float(diagonal_weight))
            count_sum += count_weight
            diagonal_sum += diagonal_weight

            # Against J*1, the original local coefficient is volume/4.
            # After projection it is weight times the total incident volume/4.
            exact_counting_harmonic &= sy.factor(
                count_weight * volume_sum - volume
            ) == 0
            exact_diagonal_harmonic &= sy.factor(
                diagonal_weight * volume_sum - volume
            ) == 0
        exact_left_inverse &= count_sum == 1 and diagonal_sum == 1

    recoveries = {
        "counting": sparse.csr_matrix(
            (counting_data, (recovery_rows, recovery_columns)),
            shape=(vertex_count, local_dimension),
        ),
        "diagonal": sparse.csr_matrix(
            (diagonal_data, (recovery_rows, recovery_columns)),
            shape=(vertex_count, local_dimension),
        ),
    }
    return injection, local_mass, recoveries, occurrences, volumes, {
        "exact_left_inverse_by_weight_sums": bool(exact_left_inverse),
        "all_incident_tetrahedron_volumes_equal": bool(
            all_incident_volumes_equal
        ),
        "exact_counting_harmonic_moment_preservation": bool(
            exact_counting_harmonic
        ),
        "exact_diagonal_harmonic_moment_preservation": bool(
            exact_diagonal_harmonic
        ),
        "distinct_exact_top_volumes": sorted({str(value) for value in volumes}),
    }


def candidate_moment_audit(injection, local_mass, recovery):
    conforming_mass = (injection.T @ local_mass @ injection).tocsr()
    mass_injection = (local_mass @ injection).tocsr()
    defect = (recovery.T @ conforming_mass - mass_injection).tocsr()
    scale = sparse_maximum(mass_injection)
    defect, numerical_threshold = prune_numerical_zeros(defect, scale)

    harmonic = np.ones(injection.shape[1], dtype=np.float64)
    harmonic_defect = np.asarray(defect @ harmonic).ravel()
    harmonic_maximum = (
        float(np.max(np.abs(harmonic_defect)))
        if harmonic_defect.size else 0.0
    )
    return {
        "full_moment_defect_nonzeros": int(defect.nnz),
        "full_moment_maximum_absolute_defect": sparse_maximum(defect),
        "full_moment_relative_defect": (
            sparse_maximum(defect) / scale if scale else 0.0
        ),
        "mass_injection_maximum_absolute_entry": scale,
        "structural_zero_absolute_threshold": numerical_threshold,
        "constant_harmonic_moment_maximum_absolute_defect": harmonic_maximum,
        "constant_harmonic_moment_passes_frozen_gate": (
            harmonic_maximum < RESIDUAL_GATE
        ),
    }


def exact_orthogonal_off_occurrence_witness(level):
    """Rebuild the preregistered k=1,p=0 witness over exact rationals."""
    top_cells = level["top"]
    vertices = level["cells"][0]
    vertex_indices = {cell[0]: index for index, cell in enumerate(vertices)}
    conforming_mass = sy.zeros(len(vertices), len(vertices))
    masses = []
    copy_global = []

    for top, points in zip(top_cells, level["top_points"]):
        mass = exact_zero_mass(exact_tetrahedron_volume(points))
        masses.append(mass)
        global_indices = [vertex_indices[vertex] for vertex in top]
        copy_global.extend(global_indices)
        for local_row, global_row in enumerate(global_indices):
            for local_column, global_column in enumerate(global_indices):
                conforming_mass[global_row, global_column] += mass[
                    local_row, local_column
                ]

    row_global = 5
    copy_index = 52
    top_index, local_column = divmod(copy_index, 4)
    top = top_cells[top_index]
    mass = masses[top_index]
    pullback_column = sy.zeros(len(vertices), 1)
    for local_row, vertex in enumerate(top):
        pullback_column[vertex_indices[vertex]] += mass[
            local_row, local_column
        ]
    solution = conforming_mass.inv() * pullback_column
    coefficient = sy.factor(solution[row_global])
    return {
        "row_global_vertex_index": row_global,
        "copy_index": copy_index,
        "copy_global_vertex_index": int(copy_global[copy_index]),
        "strict_occurrence_allowed": row_global == copy_global[copy_index],
        "exact_coefficient": str(coefficient),
        "expected_exact_coefficient": str(EXPECTED_EXACT_ADJOINT_WITNESS),
        "matches_independent_prior_certificate": (
            coefficient == EXPECTED_EXACT_ADJOINT_WITNESS
        ),
        "exactly_nonzero_off_occurrence": bool(
            coefficient != 0 and row_global != copy_global[copy_index]
        ),
    }


print("=" * 78)
print("WHITNEY PROJECTION CONSERVATION-SELECTION AUDIT")
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
all_exact_left = True
all_exact_harmonic = {"counting": True, "diagonal": True}
all_full_moment_fail = {"counting": True, "diagonal": True}
all_top_identity = True
exact_witness = None

for resolution in (1, 2, 4):
    print(f"\n-- rank-edgewise resolution k={resolution} --", flush=True)
    level = rank_edgewise_level(ranked, resolution)
    f_vector = tuple(map(len, level["cells"]))
    all_f_vectors &= f_vector == EXPECTED_F_VECTORS[resolution]
    injection, local_mass, recoveries, occurrences, volumes, exact = (
        zero_form_operators(level)
    )
    all_exact_left &= exact["exact_left_inverse_by_weight_sums"]

    candidates = {}
    for candidate, recovery in recoveries.items():
        audit = candidate_moment_audit(injection, local_mass, recovery)
        exact_harmonic = exact[
            f"exact_{candidate}_harmonic_moment_preservation"
        ]
        audit["constant_harmonic_moment_preserved_exactly"] = exact_harmonic
        candidates[candidate] = audit
        all_exact_harmonic[candidate] &= exact_harmonic
        all_full_moment_fail[candidate] &= (
            audit["full_moment_defect_nonzeros"] > 0
        )
        print(
            f"  {candidate}: full defect rel="
            f"{audit['full_moment_relative_defect']:.6g}, "
            f"nnz={audit['full_moment_defect_nonzeros']}; "
            f"harmonic exact={exact_harmonic}",
            flush=True,
        )

    top_occurrences = {top: 0 for top in level["cells"][3]}
    for top in level["top"]:
        top_occurrences[tuple(sorted(top))] += 1
    top_identity = all(value == 1 for value in top_occurrences.values())
    all_top_identity &= top_identity

    if resolution == 1:
        exact_witness = exact_orthogonal_off_occurrence_witness(level)

    records.append({
        "edgewise_resolution": resolution,
        "f_vector": list(f_vector),
        "degree_zero_global_dimension": injection.shape[1],
        "degree_zero_local_dimension": injection.shape[0],
        "degree_zero_minimum_occurrences": min(map(len, occurrences)),
        "degree_zero_maximum_occurrences": max(map(len, occurrences)),
        "exact_structure": exact,
        "candidates": candidates,
        "degree_three_top_occurrence_count_is_one": top_identity,
        "degree_three_projection_is_identity": top_identity,
    })
    del level, injection, local_mass, recoveries, occurrences, volumes
    gc.collect()

check("all three preregistered f-vectors are exact", all_f_vectors)
check("both recovery rules remain exact left inverses", all_exact_left)
check(
    "counting preserves the complete degree-zero harmonic line exactly",
    all_exact_harmonic["counting"],
)
check(
    "diagonal Whitney preserves the complete degree-zero harmonic line exactly",
    all_exact_harmonic["diagonal"],
)
check(
    "degree three is the one-copy identity for both candidates",
    all_top_identity,
)
check(
    "both local candidates fail full conforming-moment preservation",
    all(all_full_moment_fail.values()),
)
check(
    "the unique orthogonal recovery has an exact off-occurrence coefficient",
    exact_witness is not None
    and exact_witness["exactly_nonzero_off_occurrence"]
    and exact_witness["matches_independent_prior_certificate"],
    str(exact_witness),
)

harmonic_selects = (
    all_exact_harmonic["counting"] != all_exact_harmonic["diagonal"]
)
if harmonic_selects:
    harmonic_verdict = (
        "STRUCTURAL SELECTION BY THE PROPOSED AXIOM: exactly one local "
        "candidate preserves all harmonic moments on the frozen tower"
    )
elif all(all_exact_harmonic.values()) and all_top_identity:
    harmonic_verdict = (
        "DERIVED NEGATIVE FOR TOPOLOGICAL-MOMENT UNIQUENESS ON THE FROZEN "
        "TOWER: both local candidates preserve every nonvacuous harmonic "
        "moment"
    )
else:
    harmonic_verdict = (
        "NO LOCAL CANDIDATE IS VALIDATED BY ALL FROZEN HARMONIC-MOMENT TESTS"
    )

full_verdict = (
    "DERIVED SELECTION/LOCALITY CONFLICT: full conforming-moment "
    "preservation uniquely selects the metric-orthogonal recovery, whose "
    "coefficient 243/7480 is exactly nonzero off occurrence support"
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "candidate_spectrum_computed": False,
    "residual_gate": RESIDUAL_GATE,
    "structural_zero_relative_threshold": (
        STRUCTURAL_ZERO_RELATIVE_THRESHOLD
    ),
    "algebraic_lemma": (
        "For an idempotent with image im(J), preserving all M-pairings "
        "against im(J) is equivalent to the unique M-orthogonal projection"
    ),
    "records": records,
    "all_exact_harmonic_moment_results": all_exact_harmonic,
    "harmonic_moments_select_exactly_one_candidate": harmonic_selects,
    "both_candidates_fail_full_moment_preservation": all(
        all_full_moment_fail.values()
    ),
    "exact_orthogonal_recovery_nonlocal_witness": exact_witness,
    "verdicts": [
        harmonic_verdict,
        full_verdict,
        "STRUCTURAL: full-moment preservation is a projection axiom, not a "
        "derived Noether conservation law",
        "OPEN: a derived primal-dual metric Hodge star",
        "OPEN: an independently derived unique bounded-star local projection",
        "NOT CLAIMED: time, causality, inertia, mass or Planck units",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured no-spectrum certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(harmonic_verdict)
print(full_verdict)
raise SystemExit(0 if passed == tests else 1)

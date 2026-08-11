#!/usr/bin/env python3
"""Canonicity audit for a three-dimensional Whitney flux recovery.

Protocol commit 64af02c froze the carriers, three recovery formulas,
locality notion, tolerances and decision labels before these comparisons were
evaluated.  No low spectrum or phenomenological target is used here.
"""

from itertools import combinations
import gc
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
from scipy.sparse.linalg import splu
import sympy as sy

from whitney_trace_refinement_tools import (
    barycentric_refine,
    classify_element_types,
    make_base_level,
    rank_edgewise_level,
)


OUTPUT = Path(__file__).with_name("whitney_3d_flux_canonicity.json")
PROTOCOL_COMMIT = "64af02c"
SUPPORT_RELATIVE_THRESHOLD = 1e-11
SOLVE_RESIDUAL_GATE = 1e-10
EXPECTED_F_VECTORS = {
    1: (30, 150, 240, 120),
    2: (180, 1140, 1920, 960),
}
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


def occurrence_data(level, element_types, degree):
    top_cells = level["top"]
    cells = level["cells"][degree]
    local_faces = tuple(combinations(range(4), degree + 1))
    local_count = len(local_faces)
    cell_indices = {cell: index for index, cell in enumerate(cells)}
    occurrences = [[] for _ in cells]
    copy_global = np.empty(len(top_cells) * local_count, dtype=np.int32)
    copy_sign = np.empty(len(top_cells) * local_count, dtype=np.int8)
    j_rows = np.empty(len(top_cells) * local_count, dtype=np.int32)
    j_columns = np.empty(len(top_cells) * local_count, dtype=np.int32)
    j_data = np.empty(len(top_cells) * local_count, dtype=np.float64)

    cursor = 0
    for top_index, top in enumerate(top_cells):
        type_index = int(element_types["top_types"][top_index])
        local_mass = element_types["masses"][type_index][degree]
        for local_index, positions in enumerate(local_faces):
            oriented = tuple(top[position] for position in positions)
            cell = tuple(sorted(oriented))
            global_index = cell_indices[cell]
            sign = orientation_sign(oriented)
            copy_index = top_index * local_count + local_index
            diagonal = sy.factor(local_mass[local_index, local_index])
            occurrences[global_index].append({
                "copy": copy_index,
                "top": top_index,
                "local": local_index,
                "sign": sign,
                "diagonal": diagonal,
            })
            copy_global[copy_index] = global_index
            copy_sign[copy_index] = sign
            j_rows[cursor] = copy_index
            j_columns[cursor] = global_index
            j_data[cursor] = sign
            cursor += 1

    injection = sparse.csr_matrix(
        (j_data, (j_rows, j_columns)),
        shape=(len(top_cells) * local_count, len(cells)),
    )
    return local_faces, occurrences, copy_global, copy_sign, injection


def local_metric(element_types, degree):
    blocks = [
        sparse.csr_matrix(
            np.asarray(
                element_types["masses"][int(type_index)][degree],
                dtype=np.float64,
            )
        )
        for type_index in element_types["top_types"]
    ]
    return sparse.block_diag(blocks, format="csr")


def exact_local_recoveries(occurrences, global_dimension, local_dimension):
    c_rows = []
    c_columns = []
    c_data = []
    d_data = []
    differing_rows = 0
    differing_coefficients = 0
    maximum_difference = sy.Integer(0)
    first_witness = None
    all_positive = True
    exact_left_inverse = True

    for global_index, copies in enumerate(occurrences):
        count = len(copies)
        diagonal_sum = sy.factor(sum(
            (copy["diagonal"] for copy in copies), sy.Integer(0)
        ))
        counting_weight = sy.Rational(1, count)
        row_differs = False
        counting_sum = sy.Integer(0)
        diagonal_weight_sum = sy.Integer(0)
        for copy in copies:
            diagonal_weight = sy.factor(copy["diagonal"] / diagonal_sum)
            all_positive &= bool(copy["diagonal"] > 0)
            counting_sum += counting_weight
            diagonal_weight_sum += diagonal_weight
            difference = sy.factor(diagonal_weight - counting_weight)
            if difference != 0:
                row_differs = True
                differing_coefficients += 1
                if abs(float(difference)) > abs(float(maximum_difference)):
                    maximum_difference = difference
                if first_witness is None:
                    first_witness = {
                        "global_simplex_index": global_index,
                        "occurrence_count": count,
                        "top_index": copy["top"],
                        "local_simplex_index": copy["local"],
                        "local_diagonal_mass": str(copy["diagonal"]),
                        "counting_weight": str(counting_weight),
                        "diagonal_weight": str(diagonal_weight),
                        "signed_weight_difference": str(difference),
                    }
            sign = copy["sign"]
            c_rows.append(global_index)
            c_columns.append(copy["copy"])
            c_data.append(float(sign * counting_weight))
            d_data.append(float(sign * diagonal_weight))
        differing_rows += int(row_differs)
        exact_left_inverse &= (
            sy.factor(counting_sum) == 1
            and sy.factor(diagonal_weight_sum) == 1
        )

    counting = sparse.csr_matrix(
        (c_data, (c_rows, c_columns)),
        shape=(global_dimension, local_dimension),
    )
    diagonal = sparse.csr_matrix(
        (d_data, (c_rows, c_columns)),
        shape=(global_dimension, local_dimension),
    )
    return counting, diagonal, {
        "exact_left_inverse_by_weight_sums": bool(exact_left_inverse),
        "all_diagonal_masses_positive": bool(all_positive),
        "differing_rows": differing_rows,
        "differing_coefficients": differing_coefficients,
        "maximum_signed_exact_weight_difference": str(maximum_difference),
        "maximum_absolute_weight_difference": abs(float(maximum_difference)),
        "first_exact_witness": first_witness,
    }


def simplex_adjacency(level, degree):
    cells = level["cells"][degree]
    indices = {cell: index for index, cell in enumerate(cells)}
    rows = []
    columns = []
    for top in level["top"]:
        incident = [
            indices[tuple(sorted(top[position] for position in face))]
            for face in combinations(range(4), degree + 1)
        ]
        for left in incident:
            for right in incident:
                if left != right:
                    rows.append(left)
                    columns.append(right)
    adjacency = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, columns)),
        shape=(len(cells), len(cells)),
    )
    adjacency.data[:] = 1
    adjacency.eliminate_zeros()
    return adjacency


def metric_adjoint_audit(level, degree, injection, local_mass, copy_global):
    conforming_mass = (injection.T @ (local_mass @ injection)).tocsc()
    pullback = (injection.T @ local_mass).tocsc()
    factor = splu(conforming_mass)
    # SuperLU accepts the moderate dense multi-right-hand-side controls here.
    recovery = factor.solve(pullback.toarray())
    maximum = float(np.max(np.abs(recovery)))
    threshold = SUPPORT_RELATIVE_THRESHOLD * maximum

    left_inverse = recovery @ injection
    identity = np.eye(left_inverse.shape[0])
    left_residual = float(np.max(np.abs(left_inverse - identity)))
    solve_residual = conforming_mass @ recovery - pullback
    solve_scale = max(1.0, float(np.max(np.abs(pullback.data))))
    solve_relative_residual = float(
        np.max(np.abs(solve_residual)) / solve_scale
    )

    significant = np.abs(recovery) > threshold
    row_indices = np.arange(recovery.shape[0], dtype=np.int32)[:, None]
    strict_allowed = row_indices == copy_global[None, :]
    outside = significant & ~strict_allowed
    outside_count = int(np.count_nonzero(outside))
    outside_maximum = float(np.max(np.abs(recovery[outside]))) \
        if outside_count else 0.0
    maximum_outside_witness = None
    if outside_count:
        masked = np.where(outside, np.abs(recovery), -1.0)
        witness_row, witness_copy = np.unravel_index(
            int(np.argmax(masked)), masked.shape
        )
        maximum_outside_witness = {
            "row_global_simplex_index": int(witness_row),
            "copy_index": int(witness_copy),
            "copy_global_simplex_index": int(copy_global[witness_copy]),
            "numerical_coefficient": float(
                recovery[witness_row, witness_copy]
            ),
        }

    collapsed = np.zeros(
        (recovery.shape[0], recovery.shape[0]), dtype=bool
    )
    for copy_index, global_index in enumerate(copy_global):
        collapsed[:, int(global_index)] |= significant[:, copy_index]
    adjacency = simplex_adjacency(level, degree)
    distances = csgraph.shortest_path(
        adjacency, directed=False, unweighted=True
    )
    reached = collapsed & np.isfinite(distances)
    maximum_distance = int(np.max(distances[reached])) if np.any(reached) else 0
    graph_diameter = int(np.max(distances[np.isfinite(distances)]))

    del recovery, significant, outside, collapsed, distances
    gc.collect()
    return {
        "solve_relative_residual": solve_relative_residual,
        "left_inverse_maximum_residual": left_residual,
        "maximum_absolute_coefficient": maximum,
        "relative_support_threshold": SUPPORT_RELATIVE_THRESHOLD,
        "absolute_support_threshold": threshold,
        "significant_coefficients_outside_strict_occurrence": outside_count,
        "maximum_absolute_outside_coefficient": outside_maximum,
        "maximum_outside_witness": maximum_outside_witness,
        "maximum_simplex_graph_distance_reached": maximum_distance,
        "simplex_graph_diameter": graph_diameter,
        "strict_occurrence_local_numerically": outside_count == 0,
    }


def exact_metric_adjoint_witness(level, element_types, degree, witness):
    """Certify one post-result off-support coefficient over the rationals."""
    top_cells = level["top"]
    cells = level["cells"][degree]
    local_faces = tuple(combinations(range(4), degree + 1))
    local_count = len(local_faces)
    cell_indices = {cell: index for index, cell in enumerate(cells)}
    conforming_mass = sy.zeros(len(cells), len(cells))

    for top_index, top in enumerate(top_cells):
        type_index = int(element_types["top_types"][top_index])
        mass = element_types["masses"][type_index][degree]
        global_indices = []
        signs = []
        for positions in local_faces:
            oriented = tuple(top[position] for position in positions)
            global_indices.append(cell_indices[tuple(sorted(oriented))])
            signs.append(orientation_sign(oriented))
        for local_row, global_row in enumerate(global_indices):
            for local_column, global_column in enumerate(global_indices):
                conforming_mass[global_row, global_column] += (
                    signs[local_row] * mass[local_row, local_column]
                    * signs[local_column]
                )

    copy_index = witness["copy_index"]
    top_index, local_column = divmod(copy_index, local_count)
    top = top_cells[top_index]
    type_index = int(element_types["top_types"][top_index])
    mass = element_types["masses"][type_index][degree]
    pullback_column = sy.zeros(len(cells), 1)
    for local_row, positions in enumerate(local_faces):
        oriented = tuple(top[position] for position in positions)
        global_row = cell_indices[tuple(sorted(oriented))]
        pullback_column[global_row] += (
            orientation_sign(oriented) * mass[local_row, local_column]
        )
    solution = conforming_mass.inv() * pullback_column
    coefficient = sy.factor(
        solution[witness["row_global_simplex_index"]]
    )
    return {
        **witness,
        "exact_coefficient": str(coefficient),
        "exactly_nonzero": coefficient != 0,
        "absolute_numerical_disagreement": abs(
            float(coefficient) - witness["numerical_coefficient"]
        ),
    }


print("=" * 78)
print("THREE-DIMENSIONAL WHITNEY FLUX CANONICITY AUDIT")
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
all_positive = True
all_numerical_residuals = True
degree_three_control = True
exact_adjoint_nonlocal_witness = None

for resolution in (1, 2):
    print(f"\n-- rank-edgewise resolution k={resolution} --", flush=True)
    level = rank_edgewise_level(ranked, resolution)
    f_vector = tuple(map(len, level["cells"]))
    all_f_vectors &= f_vector == EXPECTED_F_VECTORS[resolution]
    element_types = classify_element_types(level)
    print(
        f"f-vector={f_vector}, element types={element_types['type_count']}",
        flush=True,
    )
    level_record = {
        "edgewise_resolution": resolution,
        "f_vector": list(f_vector),
        "top_count": len(level["top"]),
        "element_type_count": int(element_types["type_count"]),
        "degrees": [],
    }
    for degree in range(4):
        print(f"  degree p={degree}", flush=True)
        _, occurrences, copy_global, _, injection = occurrence_data(
            level, element_types, degree
        )
        local_mass = local_metric(element_types, degree)
        counting, diagonal, comparison = exact_local_recoveries(
            occurrences, injection.shape[1], injection.shape[0]
        )
        counting_residual = float(np.max(np.abs(
            (counting @ injection).toarray() - np.eye(injection.shape[1])
        )))
        diagonal_residual = float(np.max(np.abs(
            (diagonal @ injection).toarray() - np.eye(injection.shape[1])
        )))
        all_exact_left &= comparison["exact_left_inverse_by_weight_sums"]
        all_positive &= comparison["all_diagonal_masses_positive"]

        adjoint = metric_adjoint_audit(
            level, degree, injection, local_mass, copy_global
        )
        if resolution == 1 and degree == 0:
            exact_adjoint_nonlocal_witness = exact_metric_adjoint_witness(
                level,
                element_types,
                degree,
                adjoint["maximum_outside_witness"],
            )
            adjoint["post_result_exact_nonlocal_witness"] = (
                exact_adjoint_nonlocal_witness
            )
        all_numerical_residuals &= (
            adjoint["solve_relative_residual"] < SOLVE_RESIDUAL_GATE
            and adjoint["left_inverse_maximum_residual"]
            < SOLVE_RESIDUAL_GATE
        )
        if degree == 3:
            degree_three_control &= (
                comparison["differing_rows"] == 0
                and adjoint[
                    "significant_coefficients_outside_strict_occurrence"
                ] == 0
            )
        degree_record = {
            "degree": degree,
            "global_dimension": injection.shape[1],
            "local_dimension": injection.shape[0],
            "minimum_occurrence_count": min(map(len, occurrences)),
            "maximum_occurrence_count": max(map(len, occurrences)),
            "counting_left_inverse_numerical_residual": counting_residual,
            "diagonal_left_inverse_numerical_residual": diagonal_residual,
            "counting_vs_diagonal": comparison,
            "metric_adjoint": adjoint,
        }
        level_record["degrees"].append(degree_record)
        print(
            "    C!=D rows={rows}; adjoint off-support={outside}; "
            "radius={radius}/{diameter}".format(
                rows=comparison["differing_rows"],
                outside=adjoint[
                    "significant_coefficients_outside_strict_occurrence"
                ],
                radius=adjoint["maximum_simplex_graph_distance_reached"],
                diameter=adjoint["simplex_graph_diameter"],
            ),
            flush=True,
        )
        del local_mass, counting, diagonal, injection
        gc.collect()
    records.append(level_record)
    del level, element_types
    gc.collect()

check("both frozen exact f-vectors are recovered", all_f_vectors)
check("all exact local Whitney diagonal masses are positive", all_positive)
check("both local recovery formulas are exact left inverses", all_exact_left)
check("all floating reconstructions of the two exact identities are stable",
      all(
          degree["counting_left_inverse_numerical_residual"] < 1e-14
          and degree["diagonal_left_inverse_numerical_residual"] < 1e-14
          for level in records for degree in level["degrees"]
      ))
check("all metric-adjoint solves and left inverses pass the frozen residual",
      all_numerical_residuals)
check("degree three is the frozen one-copy identity control",
      degree_three_control)
check("one metric-adjoint off-support coefficient is exactly nonzero",
      exact_adjoint_nonlocal_witness is not None
      and exact_adjoint_nonlocal_witness["exactly_nonzero"]
      and exact_adjoint_nonlocal_witness[
          "absolute_numerical_disagreement"
      ] < 1e-12,
      str(exact_adjoint_nonlocal_witness))

propagating = [
    degree for level in records for degree in level["degrees"]
    if degree["degree"] < 3
]
local_ambiguity = any(
    degree["counting_vs_diagonal"]["differing_rows"] > 0
    for degree in propagating
)
adjoint_nonlocal = any(
    degree["metric_adjoint"][
        "significant_coefficients_outside_strict_occurrence"
    ] > 0
    for degree in propagating
)
if local_ambiguity:
    uniqueness_verdict = (
        "DERIVED NEGATIVE FOR UNIQUENESS: counting and diagonal-Whitney "
        "recoveries are distinct exact local natural left inverses"
    )
else:
    uniqueness_verdict = (
        "OPEN, NOT UNIQUE: the two frozen local recoveries agree on these "
        "controls, which is not a uniqueness theorem"
    )
if adjoint_nonlocal:
    adjoint_verdict = (
        "DERIVED LOCALITY CONFLICT: the unique exact-metric adjoint recovery "
        "has an exactly certified coefficient outside the occurrence star"
    )
else:
    adjoint_verdict = (
        "NO NUMERICAL LOCALITY CONFLICT DETECTED on the frozen controls"
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "low_spectrum_computed": False,
    "support_relative_threshold": SUPPORT_RELATIVE_THRESHOLD,
    "solve_residual_gate": SOLVE_RESIDUAL_GATE,
    "records": records,
    "local_ambiguity_detected": local_ambiguity,
    "metric_adjoint_nonlocal_detected": adjoint_nonlocal,
    "post_result_exact_adjoint_nonlocal_witness": (
        exact_adjoint_nonlocal_witness
    ),
    "verdicts": [
        uniqueness_verdict,
        adjoint_verdict,
        "OPEN: a uniquely selected bounded-star three-dimensional flux",
        "NOT CLAIMED: continuum physics, time, mass, c or Planck units",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured no-target certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(uniqueness_verdict)
print(adjoint_verdict)
raise SystemExit(0 if passed == tests else 1)

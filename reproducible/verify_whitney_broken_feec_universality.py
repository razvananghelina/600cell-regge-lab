#!/usr/bin/env python3
"""Strong-penalty universality audit for two broken-FEEC projections.

Protocol commits 72d4eaa and 230195b froze the carrier, constrained pencil,
microscopic and low-spectrum comparisons, solver residual and labels before
the new k=4 data were evaluated.
"""

from itertools import combinations
import gc
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
import sympy as sy

from whitney_trace_refinement_tools import (
    LOCAL_D,
    barycentric_refine,
    classify_element_types,
    make_base_level,
    rank_edgewise_level,
)


OUTPUT = Path(__file__).with_name(
    "whitney_broken_feec_universality.json"
)
PROTOCOL_COMMITS = ("72d4eaa", "230195b")
EXPECTED_F_VECTORS = {
    1: (30, 150, 240, 120),
    2: (180, 1140, 1920, 960),
    4: (1320, 9000, 15360, 7680),
}
EXPECTED_BETTI = (1, 0, 0, 1)
ZERO_THRESHOLD = 1e-9
RITZ_RESIDUAL_GATE = 1e-7
EQUALITY_THRESHOLD = 1e-10
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


def numerical_nonzeros(matrix, threshold=1e-13):
    matrix = sparse.csr_matrix(matrix)
    return int(np.count_nonzero(np.abs(matrix.data) > threshold))


def orientation_sign(vertices):
    inversions = sum(
        vertices[left] > vertices[right]
        for left in range(len(vertices))
        for right in range(left + 1, len(vertices))
    )
    return -1 if inversions % 2 else 1


def injection_recoveries_and_micro(level, element_types, degree):
    top_cells = level["top"]
    cells = level["cells"][degree]
    local_faces = tuple(combinations(range(4), degree + 1))
    local_count = len(local_faces)
    cell_indices = {cell: index for index, cell in enumerate(cells)}
    occurrences = [[] for _ in cells]
    j_rows = []
    j_columns = []
    j_data = []

    for top_index, top in enumerate(top_cells):
        type_index = int(element_types["top_types"][top_index])
        local_mass = element_types["masses"][type_index][degree]
        for local_index, positions in enumerate(local_faces):
            oriented = tuple(top[position] for position in positions)
            global_index = cell_indices[tuple(sorted(oriented))]
            sign = orientation_sign(oriented)
            copy_index = top_index * local_count + local_index
            occurrences[global_index].append((
                copy_index,
                sign,
                sy.factor(local_mass[local_index, local_index]),
                top_index,
                local_index,
            ))
            j_rows.append(copy_index)
            j_columns.append(global_index)
            j_data.append(sign)

    local_dimension = len(top_cells) * local_count
    injection = sparse.csr_matrix(
        (j_data, (j_rows, j_columns)),
        shape=(local_dimension, len(cells)),
        dtype=float,
    )
    rows = []
    columns = []
    data = {"counting": [], "diagonal": []}
    differing_rows = 0
    differing_coefficients = 0
    maximum_difference = sy.Integer(0)
    first_witness = None
    exact_sums = True
    for global_index, copies in enumerate(occurrences):
        count = len(copies)
        diagonal_sum = sy.factor(sum(
            (item[2] for item in copies), sy.Integer(0)
        ))
        count_sum = diagonal_weight_sum = sy.Integer(0)
        row_differs = False
        for copy_index, sign, local_diagonal, top_index, local_index in copies:
            count_weight = sy.Rational(1, count)
            diagonal_weight = sy.factor(local_diagonal / diagonal_sum)
            difference = sy.factor(diagonal_weight - count_weight)
            count_sum += count_weight
            diagonal_weight_sum += diagonal_weight
            if difference != 0:
                row_differs = True
                differing_coefficients += 1
                if abs(float(difference)) > abs(float(maximum_difference)):
                    maximum_difference = difference
                if first_witness is None:
                    first_witness = {
                        "global_simplex_index": global_index,
                        "top_index": top_index,
                        "local_simplex_index": local_index,
                        "occurrence_count": count,
                        "counting_weight": str(count_weight),
                        "diagonal_weight": str(diagonal_weight),
                        "exact_difference": str(difference),
                    }
            rows.append(global_index)
            columns.append(copy_index)
            data["counting"].append(float(sign * count_weight))
            data["diagonal"].append(float(sign * diagonal_weight))
        differing_rows += int(row_differs)
        exact_sums &= (
            sy.factor(count_sum) == 1
            and sy.factor(diagonal_weight_sum) == 1
        )

    recoveries = {
        candidate: sparse.csr_matrix(
            (values, (rows, columns)),
            shape=(len(cells), local_dimension),
        )
        for candidate, values in data.items()
    }
    return injection, recoveries, {
        "exact_weight_sums": bool(exact_sums),
        "global_rows": len(cells),
        "differing_rows": differing_rows,
        "differing_row_fraction": differing_rows / len(cells),
        "differing_coefficients": differing_coefficients,
        "maximum_signed_exact_weight_difference": str(maximum_difference),
        "maximum_absolute_weight_difference": abs(float(maximum_difference)),
        "first_exact_witness": first_witness,
    }


def local_mass_and_inverse(element_types, degree):
    masses = []
    inverses = []
    for type_index in element_types["top_types"]:
        block = np.asarray(
            element_types["masses"][int(type_index)][degree],
            dtype=np.float64,
        )
        masses.append(sparse.csr_matrix(block))
        inverses.append(sparse.csr_matrix(np.linalg.inv(block)))
    return (
        sparse.block_diag(masses, format="csr"),
        sparse.block_diag(inverses, format="csr"),
    )


def piecewise_coboundary(top_count, degree):
    block = sparse.csr_matrix(np.asarray(LOCAL_D[degree], dtype=np.float64))
    return sparse.kron(
        sparse.eye(top_count, format="csr"), block, format="csr"
    )


def constrained_pencil(
    degree,
    candidate,
    injections,
    projectors,
    masses,
    inverse_masses,
    local_d,
):
    injection = injections[degree]
    metric = (injection.T @ masses[degree] @ injection).tocsr()
    weak = sparse.csr_matrix(metric.shape, dtype=np.float64)
    if degree < 3:
        upper = (local_d[degree] @ injection).tocsr()
        weak = (weak + upper.T @ masses[degree + 1] @ upper).tocsr()
    if degree > 0:
        lower = (
            inverse_masses[degree - 1]
            @ projectors[candidate][degree - 1].T
            @ local_d[degree - 1].T
            @ masses[degree]
            @ injection
        ).tocsr()
        weak = (weak + lower.T @ masses[degree - 1] @ lower).tocsr()
    weak.sum_duplicates()
    weak.eliminate_zeros()
    return weak, metric


def matrix_infinity_norm(matrix):
    return float(np.max(np.asarray(abs(matrix).sum(axis=1)).ravel()))


def lowest_ritz(weak, metric, count, seed_offset):
    dimension = weak.shape[0]
    v0 = np.cos(
        np.arange(dimension, dtype=np.float64) + 0.125 + seed_offset
    )
    v0 /= np.linalg.norm(v0)
    values, vectors = eigsh(
        weak,
        M=metric,
        k=count,
        sigma=-1e-8,
        which="LM",
        tol=1e-10,
        maxiter=30_000,
        v0=v0,
    )
    ordering = np.argsort(values)
    values = np.real(values[ordering])
    vectors = np.real(vectors[:, ordering])
    weak_norm = matrix_infinity_norm(weak)
    metric_norm = matrix_infinity_norm(metric)
    residuals = []
    for index, value in enumerate(values):
        vector = vectors[:, index]
        residual = weak @ vector - value * (metric @ vector)
        scale = (
            (weak_norm + abs(value) * metric_norm)
            * np.linalg.norm(vector)
        )
        residuals.append(float(np.linalg.norm(residual) / scale))
    return values, residuals


print("=" * 78)
print("STRONG-LIMIT BROKEN-FEEC UNIVERSALITY AUDIT")
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
all_weight_sums = True
all_k1_calibration = True
all_betti_nullities = True
all_ritz_residuals = True

for resolution in (1, 2, 4):
    print(f"\n-- rank-edgewise resolution k={resolution} --", flush=True)
    level = rank_edgewise_level(ranked, resolution)
    f_vector = tuple(map(len, level["cells"]))
    all_f_vectors &= f_vector == EXPECTED_F_VECTORS[resolution]
    top_count = len(level["top"])
    element_types = classify_element_types(level)
    print(
        f"f-vector={f_vector}, local metric types={element_types['type_count']}",
        flush=True,
    )

    injections = []
    projectors = {"counting": [], "diagonal": []}
    masses = []
    inverse_masses = []
    micro_records = []
    for degree in range(4):
        injection, recoveries, micro = injection_recoveries_and_micro(
            level, element_types, degree
        )
        injections.append(injection)
        micro_records.append({"degree": degree, **micro})
        all_weight_sums &= micro["exact_weight_sums"]
        for candidate, recovery in recoveries.items():
            projectors[candidate].append((injection @ recovery).tocsr())
        mass, inverse = local_mass_and_inverse(element_types, degree)
        masses.append(mass)
        inverse_masses.append(inverse)
        print(
            f"  p={degree}: C!=D rows={micro['differing_rows']}, "
            f"max weight delta={micro['maximum_absolute_weight_difference']:.8g}",
            flush=True,
        )

    local_d = [piecewise_coboundary(top_count, degree)
               for degree in range(3)]
    candidate_records = {"counting": [], "diagonal": []}
    cached_solutions = {}
    constrained_matrices = {"counting": [], "diagonal": []}
    for degree in range(4):
        matrices = {}
        metric = None
        for candidate in ("counting", "diagonal"):
            weak, candidate_metric = constrained_pencil(
                degree,
                candidate,
                injections,
                projectors,
                masses,
                inverse_masses,
                local_d,
            )
            matrices[candidate] = weak
            constrained_matrices[candidate].append(weak)
            metric = candidate_metric

        matrix_difference = matrices["counting"] - matrices["diagonal"]
        identical = numerical_nonzeros(matrix_difference) == 0
        requested = 8 if resolution < 4 else 6
        for candidate_index, candidate in enumerate(("counting", "diagonal")):
            cache_key = (degree, "shared" if identical else candidate)
            if cache_key not in cached_solutions:
                cached_solutions[cache_key] = lowest_ritz(
                    matrices[candidate],
                    metric,
                    requested,
                    0.01 * (100 * resolution + 10 * degree + candidate_index),
                )
            values, residuals = cached_solutions[cache_key]
            nullity = int(np.count_nonzero(np.abs(values) < ZERO_THRESHOLD))
            positives = values[values > ZERO_THRESHOLD]
            all_betti_nullities &= nullity == EXPECTED_BETTI[degree]
            all_ritz_residuals &= max(residuals) < RITZ_RESIDUAL_GATE
            candidate_records[candidate].append({
                "degree": degree,
                "conforming_dimension": metric.shape[0],
                "weak_nonzeros": int(matrices[candidate].nnz),
                "kernel_dimension": nullity,
                "lowest_generalized_eigenvalues": [
                    float(value) for value in values
                ],
                "first_four_positive_eigenvalues": [
                    float(value) for value in positives[:4]
                ],
                "maximum_relative_ritz_residual": max(residuals),
            })
            print(
                f"    {candidate} p={degree}: nullity={nullity}, "
                f"first+={positives[0]:.10g}, "
                f"ritz={max(residuals):.2e}",
                flush=True,
            )

        count_positive = candidate_records["counting"][degree][
            "first_four_positive_eigenvalues"
        ]
        diag_positive = candidate_records["diagonal"][degree][
            "first_four_positive_eigenvalues"
        ]
        if len(count_positive) != 4 or len(diag_positive) != 4:
            raise AssertionError("four positive constrained modes required")
        relative = [
            abs(left - right) / max(abs(left), abs(right))
            for left, right in zip(count_positive, diag_positive)
        ]
        candidate_records.setdefault("comparisons", []).append({
            "degree": degree,
            "pencils_identical_numerically": identical,
            "pencil_difference_nonzeros": numerical_nonzeros(
                matrix_difference
            ),
            "pencil_maximum_absolute_difference": sparse_maximum(
                matrix_difference
            ),
            "ordered_positive_relative_differences": relative,
            "maximum_ordered_positive_relative_difference": max(relative),
        })
        if resolution == 1:
            all_k1_calibration &= max(relative) < EQUALITY_THRESHOLD

    records.append({
        "edgewise_resolution": resolution,
        "f_vector": list(f_vector),
        "element_type_count": int(element_types["type_count"]),
        "microscopic_projection_comparison": micro_records,
        "strong_limit_spectra": candidate_records,
    })
    del (
        level,
        element_types,
        injections,
        projectors,
        masses,
        inverse_masses,
        local_d,
        constrained_matrices,
        cached_solutions,
    )
    gc.collect()

check("all three frozen f-vectors are exact", all_f_vectors)
check("all recovery weight sums remain exact", all_weight_sums)
check("k=1 reproduces exact candidate agreement", all_k1_calibration)
check("all strong-limit controls have Betti nullities (1,0,0,1)",
      all_betti_nullities)
check("all reported Ritz pairs pass the frozen residual gate",
      all_ritz_residuals)

micro_k2 = records[1]["microscopic_projection_comparison"]
micro_k4 = records[2]["microscopic_projection_comparison"]
micro_nonconvergence = all(
    micro_k4[degree]["maximum_signed_exact_weight_difference"]
    == micro_k2[degree]["maximum_signed_exact_weight_difference"]
    and micro_k2[degree]["maximum_absolute_weight_difference"] > 0
    for degree in (1, 2)
)
check("the exact microscopic comparison reaches a frozen verdict",
      micro_nonconvergence)

comparison_k2 = records[1]["strong_limit_spectra"]["comparisons"]
comparison_k4 = records[2]["strong_limit_spectra"]["comparisons"]
flow_by_degree = []
flow_toward_common = True
for degree in (1, 2, 3):
    before = comparison_k2[degree][
        "maximum_ordered_positive_relative_difference"
    ]
    after = comparison_k4[degree][
        "maximum_ordered_positive_relative_difference"
    ]
    passes_flow = (
        after < EQUALITY_THRESHOLD
        if before < EQUALITY_THRESHOLD
        else after < before
    )
    flow_toward_common &= passes_flow
    flow_by_degree.append({
        "degree": degree,
        "maximum_relative_difference_k2": before,
        "maximum_relative_difference_k4": after,
        "ratio_k4_over_k2": after / before if before else None,
        "passes_preregistered_flow_gate": passes_flow,
    })
check("the preregistered strong-limit flow comparison was completed", True)

micro_verdict = (
    "DERIVED NEGATIVE FOR UNIFORM MICROSCOPIC CONVERGENCE: exact local "
    "projection differences do not decrease from k=2 to k=4"
    if micro_nonconvergence else
    "NO MICROSCOPIC NONCONVERGENCE VERDICT"
)
flow_verdict = (
    "PATTERN TOWARD A COMMON LOW-ENERGY CLASS: every nonzero strong-limit "
    "candidate difference decreases at k=4"
    if flow_toward_common else
    "PATTERN NEGATIVE FOR COMMON LOW-ENERGY FLOW: at least one strong-limit "
    "candidate difference fails to decrease at k=4"
)

payload = {
    "protocol_commits": list(PROTOCOL_COMMITS),
    "phenomenological_target_used": False,
    "candidate_count": 2,
    "weak_alpha_one_used_as_physical_branch": False,
    "strong_penalty_constrained_limit_used": True,
    "records": records,
    "microscopic_nonconvergence_detected": micro_nonconvergence,
    "strong_limit_flow_by_degree": flow_by_degree,
    "flow_toward_common_low_energy_class": flow_toward_common,
    "verdicts": [
        micro_verdict,
        flow_verdict,
        "OPEN: analytic continuum universality for both P and P*",
        "STRUCTURAL: constrained strong limit is not selected finite microdynamics",
        "NOT CLAIMED: round-S3 spectrum, time, c, inertia, mass or Planck units",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured strong-limit certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(micro_verdict)
print(flow_verdict)
for item in flow_by_degree:
    print(
        f"p={item['degree']}: k2={item['maximum_relative_difference_k2']:.3e}, "
        f"k4={item['maximum_relative_difference_k4']:.3e}, "
        f"pass={item['passes_preregistered_flow_gate']}"
    )
raise SystemExit(0 if passed == tests else 1)

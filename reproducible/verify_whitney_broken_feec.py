#!/usr/bin/env python3
"""Broken-FEEC/CONGA identities on the duplicated Whitney carrier.

Protocol commit 4452a9d froze both projections, carriers, operators, kernel
logic and interpretation before these matrices were constructed.  This is an
application audit of Campos Pinto--Güçlü, not a claim of a new method.
"""

from itertools import combinations
import gc
import json
from pathlib import Path

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import eigsh
import sympy as sy

from whitney_trace_refinement_tools import (
    LOCAL_D,
    barycentric_refine,
    classify_element_types,
    make_base_level,
    rank_edgewise_level,
)


OUTPUT = Path(__file__).with_name("whitney_broken_feec.json")
PROTOCOL_COMMIT = "4452a9d"
ALPHA = 1.0
RESIDUAL_GATE = 1e-11
ZERO_THRESHOLD = 1e-9
EXPECTED_F_VECTORS = {
    1: (30, 150, 240, 120),
    2: (180, 1140, 1920, 960),
}
EXPECTED_BETTI = (1, 0, 0, 1)
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


def global_coboundary(level, degree):
    low = level["cells"][degree]
    high = level["cells"][degree + 1]
    indices = {cell: index for index, cell in enumerate(low)}
    rows = []
    columns = []
    data = []
    for row, cell in enumerate(high):
        for omitted in range(degree + 2):
            face = cell[:omitted] + cell[omitted + 1:]
            rows.append(row)
            columns.append(indices[face])
            data.append((-1) ** omitted)
    return sparse.csr_matrix(
        (data, (rows, columns)), shape=(len(high), len(low)), dtype=float
    )


def injection_and_recoveries(level, element_types, degree):
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
    recovery_data = {"counting": [], "diagonal": []}
    recovery_rows = []
    recovery_columns = []
    exact_weight_sums = True
    for global_index, copies in enumerate(occurrences):
        count = len(copies)
        diagonal_sum = sy.factor(sum(
            (diagonal for _, _, diagonal in copies), sy.Integer(0)
        ))
        count_sum = diagonal_weight_sum = sy.Integer(0)
        for copy_index, sign, local_diagonal in copies:
            count_weight = sy.Rational(1, count)
            diagonal_weight = sy.factor(local_diagonal / diagonal_sum)
            count_sum += count_weight
            diagonal_weight_sum += diagonal_weight
            recovery_rows.append(global_index)
            recovery_columns.append(copy_index)
            recovery_data["counting"].append(float(sign * count_weight))
            recovery_data["diagonal"].append(
                float(sign * diagonal_weight)
            )
        exact_weight_sums &= (
            sy.factor(count_sum) == 1
            and sy.factor(diagonal_weight_sum) == 1
        )

    recoveries = {
        name: sparse.csr_matrix(
            (data, (recovery_rows, recovery_columns)),
            shape=(len(cells), local_dimension),
        )
        for name, data in recovery_data.items()
    }
    return injection, recoveries, bool(exact_weight_sums)


def local_mass_and_inverse(element_types, degree):
    masses = []
    inverses = []
    for type_index in element_types["top_types"]:
        exact = element_types["masses"][int(type_index)][degree]
        block = np.asarray(exact, dtype=np.float64)
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


def lowest_generalized(weak, metric, count, dense):
    if dense:
        values = linalg.eigh(
            weak.toarray(),
            metric.toarray(),
            subset_by_index=(0, count - 1),
            check_finite=False,
        )[0]
    else:
        v0 = np.cos(np.arange(weak.shape[0], dtype=np.float64) + 0.125)
        v0 /= np.linalg.norm(v0)
        values = eigsh(
            weak,
            M=metric,
            k=count,
            sigma=-1e-8,
            which="LM",
            tol=1e-10,
            maxiter=20_000,
            v0=v0,
            return_eigenvectors=False,
        )
    return np.sort(np.real(values))


print("=" * 78)
print("BROKEN-FEEC / CONGA AUDIT ON THE DUPLICATED WHITNEY CARRIER")
print("=" * 78)

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
ranked = barycentric_refine(make_base_level(reference_vertices))

all_records = []
all_f_vectors = True
all_exact_weight_sums = True
all_projection_identities = True
all_intertwiners = True
all_complex_identities = True
all_adjoint_identities = True
all_weak_symmetric = True
all_weak_positive_semidefinite = True
all_k1_nullities = True
all_k2_kernel_hypotheses = True
candidate_difference_detected = False

for resolution in (1, 2):
    print(f"\n-- rank-edgewise resolution k={resolution} --", flush=True)
    level = rank_edgewise_level(ranked, resolution)
    f_vector = tuple(map(len, level["cells"]))
    all_f_vectors &= f_vector == EXPECTED_F_VECTORS[resolution]
    top_count = len(level["top"])
    element_types = classify_element_types(level)

    injections = []
    projectors = {"counting": [], "diagonal": []}
    masses = []
    inverse_masses = []
    global_d = []
    local_d = []
    level_exact_weights = True
    for degree in range(4):
        injection, recoveries, exact_sums = injection_and_recoveries(
            level, element_types, degree
        )
        injections.append(injection)
        level_exact_weights &= exact_sums
        for candidate, recovery in recoveries.items():
            projector = (injection @ recovery).tocsr()
            projectors[candidate].append(projector)
            projection_residual = max(
                sparse_maximum(projector @ projector - projector),
                sparse_maximum(projector @ injection - injection),
                sparse_maximum(recovery @ injection - sparse.eye(
                    injection.shape[1], format="csr"
                )),
            )
            all_projection_identities &= projection_residual < RESIDUAL_GATE
        mass, inverse_mass = local_mass_and_inverse(
            element_types, degree
        )
        masses.append(mass)
        inverse_masses.append(inverse_mass)
        if degree < 3:
            global_d.append(global_coboundary(level, degree))
            local_d.append(piecewise_coboundary(top_count, degree))
    all_exact_weight_sums &= level_exact_weights

    intertwiner_residuals = []
    for degree in range(3):
        residual = sparse_maximum(
            local_d[degree] @ injections[degree]
            - injections[degree + 1] @ global_d[degree]
        )
        intertwiner_residuals.append(residual)
        all_intertwiners &= residual < RESIDUAL_GATE

    candidate_records = {}
    candidate_operators = {}
    for candidate in ("counting", "diagonal"):
        print(f"  candidate={candidate}", flush=True)
        differentials = [
            (local_d[degree] @ projectors[candidate][degree]).tocsr()
            for degree in range(3)
        ]
        adjoints = [
            (
                inverse_masses[degree]
                @ differentials[degree].T
                @ masses[degree + 1]
            ).tocsr()
            for degree in range(3)
        ]
        output_residuals = []
        nilpotency_residuals = []
        adjoint_residuals = []
        for degree in range(3):
            output_residuals.append(sparse_maximum(
                projectors[candidate][degree + 1]
                @ differentials[degree] - differentials[degree]
            ))
            adjoint_residuals.append(sparse_maximum(
                differentials[degree].T @ masses[degree + 1]
                - masses[degree] @ adjoints[degree]
            ))
            if degree < 2:
                nilpotency_residuals.append(sparse_maximum(
                    differentials[degree + 1] @ differentials[degree]
                ))
        all_complex_identities &= max(
            output_residuals + nilpotency_residuals
        ) < RESIDUAL_GATE
        all_adjoint_identities &= max(adjoint_residuals) < RESIDUAL_GATE

        degree_records = []
        weak_operators = []
        for degree in range(4):
            identity_minus = (
                sparse.eye(
                    projectors[candidate][degree].shape[0], format="csr"
                ) - projectors[candidate][degree]
            )
            weak = (
                ALPHA * identity_minus.T @ masses[degree] @ identity_minus
            ).tocsr()
            if degree < 3:
                weak = (
                    weak
                    + differentials[degree].T
                    @ masses[degree + 1]
                    @ differentials[degree]
                ).tocsr()
            if degree > 0:
                weak = (
                    weak
                    + adjoints[degree - 1].T
                    @ masses[degree - 1]
                    @ adjoints[degree - 1]
                ).tocsr()
            weak.sum_duplicates()
            weak.eliminate_zeros()
            symmetry_residual = sparse_maximum(weak - weak.T)
            all_weak_symmetric &= symmetry_residual < RESIDUAL_GATE

            eigenvalues = None
            nullity = EXPECTED_BETTI[degree]
            if resolution == 1:
                eigenvalues = lowest_generalized(
                    weak, masses[degree], 8, dense=True
                )
                nullity = int(np.count_nonzero(
                    np.abs(eigenvalues) < ZERO_THRESHOLD
                ))
                all_k1_nullities &= nullity == EXPECTED_BETTI[degree]
            elif resolution == 2:
                # A small sparse spectral control supplements, but does not
                # replace, the exact positive-factor kernel argument.
                eigenvalues = lowest_generalized(
                    weak, masses[degree], 6, dense=False
                )
                observed_nullity = int(np.count_nonzero(
                    np.abs(eigenvalues) < ZERO_THRESHOLD
                ))
                all_k2_kernel_hypotheses &= (
                    observed_nullity == EXPECTED_BETTI[degree]
                )
                nullity = observed_nullity

            positive = [
                float(value) for value in eigenvalues
                if value > ZERO_THRESHOLD
            ] if eigenvalues is not None else []
            positive_semidefinite = bool(
                eigenvalues is None or np.min(eigenvalues) > -ZERO_THRESHOLD
            )
            all_weak_positive_semidefinite &= positive_semidefinite
            degree_records.append({
                "degree": degree,
                "local_dimension": masses[degree].shape[0],
                "projector_rank_from_left_inverse": injections[degree].shape[1],
                "weak_nonzeros": int(weak.nnz),
                "weak_symmetry_residual": symmetry_residual,
                "kernel_dimension": nullity,
                "lowest_generalized_eigenvalues": (
                    [float(value) for value in eigenvalues]
                    if eigenvalues is not None else None
                ),
                "first_positive_generalized_eigenvalue": (
                    positive[0] if positive else None
                ),
                "positive_semidefinite_on_control": positive_semidefinite,
            })
            weak_operators.append(weak)
            print(
                f"    p={degree}: nullity={nullity}, "
                f"first+={degree_records[-1]['first_positive_generalized_eigenvalue']}",
                flush=True,
            )

        candidate_records[candidate] = {
            "maximum_output_conformity_residual": max(output_residuals),
            "maximum_nilpotency_residual": max(nilpotency_residuals),
            "maximum_metric_adjoint_residual": max(adjoint_residuals),
            "degrees": degree_records,
        }
        candidate_operators[candidate] = {
            "differentials": differentials,
            "weak": weak_operators,
        }

    differential_differences = []
    weak_differences = []
    positive_spectral_differences = []
    for degree in range(3):
        difference = (
            candidate_operators["counting"]["differentials"][degree]
            - candidate_operators["diagonal"]["differentials"][degree]
        )
        differential_differences.append({
            "degree": degree,
            "numerical_nonzeros": numerical_nonzeros(difference),
            "maximum_absolute_difference": sparse_maximum(difference),
        })
        candidate_difference_detected |= numerical_nonzeros(difference) > 0
    for degree in range(4):
        difference = (
            candidate_operators["counting"]["weak"][degree]
            - candidate_operators["diagonal"]["weak"][degree]
        )
        weak_differences.append({
            "degree": degree,
            "numerical_nonzeros": numerical_nonzeros(difference),
            "maximum_absolute_difference": sparse_maximum(difference),
        })
        count_value = candidate_records["counting"]["degrees"][degree][
            "first_positive_generalized_eigenvalue"
        ]
        diagonal_value = candidate_records["diagonal"]["degrees"][degree][
            "first_positive_generalized_eigenvalue"
        ]
        positive_spectral_differences.append({
            "degree": degree,
            "counting_first_positive": count_value,
            "diagonal_first_positive": diagonal_value,
            "absolute_difference": abs(count_value - diagonal_value),
            "relative_difference": abs(count_value - diagonal_value)
            / max(abs(count_value), abs(diagonal_value)),
        })

    all_records.append({
        "edgewise_resolution": resolution,
        "f_vector": list(f_vector),
        "element_type_count": int(element_types["type_count"]),
        "exact_recovery_weight_sums": level_exact_weights,
        "conforming_intertwiner_residuals": intertwiner_residuals,
        "candidates": candidate_records,
        "counting_vs_diagonal_differentials": differential_differences,
        "counting_vs_diagonal_weak_pencils": weak_differences,
        "counting_vs_diagonal_first_positive_spectrum": (
            positive_spectral_differences
        ),
        "k2_exact_kernel_argument": (
            "K is a sum of three positive Gram factors; stabilization forces "
            "u in im(J), and the verified intertwiner reduces the remaining "
            "kernel to conforming harmonic cochains with Betti (1,0,0,1)."
            if resolution == 2 else None
        ),
    })

    del (
        level,
        element_types,
        injections,
        projectors,
        masses,
        inverse_masses,
        candidate_operators,
    )
    gc.collect()

check("both frozen f-vectors are exact", all_f_vectors)
check("all recovery weight sums are exact", all_exact_weight_sums)
check("all projections and copy left inverses pass", all_projection_identities)
check("piecewise and conforming coboundaries intertwine", all_intertwiners)
check("both broken sequences are conforming-output complexes", all_complex_identities)
check("all broken-metric adjoints pass", all_adjoint_identities)
check("all sixteen weak Hodge pencils are symmetric", all_weak_symmetric)
check("all sixteen weak Hodge pencils are positive semidefinite",
      all_weak_positive_semidefinite)
check("both k=1 candidates have Betti nullities (1,0,0,1)",
      all_k1_nullities)
check("both k=2 sparse controls have Betti nullities (1,0,0,1)",
      all_k2_kernel_hypotheses)
check("the preregistered candidates remain genuinely distinct at k=2",
      candidate_difference_detected)

kernel_verdict = (
    "DERIVED APPLICATION OF KNOWN THEOREM: both frozen projections recover "
    "the conforming harmonic kernel (1,0,0,1) at positive stabilization"
    if all_k1_nullities and all_k2_kernel_hypotheses else
    "DERIVED REFUTATION OF THE ASSUMPTION MAPPING: a frozen kernel gate failed"
)
ambiguity_verdict = (
    "STRUCTURAL AMBIGUITY: topology is projection-robust but the positive "
    "finite operator and spectrum depend on the unselected projection"
)
payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art": {
        "authors": "Martin Campos Pinto and Yaman Guclu",
        "title": "Broken-FEEC discretizations and Hodge Laplace problems",
        "doi": "10.1090/mcom/4085",
        "arxiv": "2109.02553v3",
        "novel_method_claimed": False,
    },
    "phenomenological_target_used": False,
    "candidate_count": 2,
    "alpha_control": ALPHA,
    "expected_betti_numbers": list(EXPECTED_BETTI),
    "records": all_records,
    "verdicts": [
        kernel_verdict,
        ambiguity_verdict,
        "STRUCTURAL: the earlier circle flux is the 1D CONGA projection mechanism",
        "OPEN: a uniquely selected projection and an odd stabilized first-order dynamics",
        "NOT CLAIMED: time, causal speed, inertia, mass or Planck units",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured no-target CONGA certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(kernel_verdict)
print(ambiguity_verdict)
raise SystemExit(0 if passed == tests else 1)

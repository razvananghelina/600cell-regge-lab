#!/usr/bin/env python3
"""Exact assembly-invariance gate for element-local Whitney dynamics.

Protocol commit 162ce61 froze the 9000-dimensional duplicated carrier,
all-degree invariance tests and projected weak identity before computation.
"""

from itertools import combinations
import json
from math import factorial
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sparse
from scipy.sparse.csgraph import structural_rank
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_element_local_assembly.json")
PROTOCOL_COMMIT = "162ce61"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def exact_integer_matrix(matrix):
    denominators = [sy.denom(value) for value in matrix]
    if not denominators:
        denominator = sy.Integer(1)
    elif len(denominators) == 1:
        denominator = denominators[0]
    else:
        denominator = sy.ilcm(*denominators)
    numerator = np.asarray(
        [[int(sy.expand(denominator * matrix[row, column]))
          for column in range(matrix.cols)]
         for row in range(matrix.rows)],
        dtype=np.int64,
    )
    return numerator, int(denominator)


def permutation_sign(sequence):
    return -1 if sum(
        sequence[left] > sequence[right]
        for left in range(len(sequence))
        for right in range(left + 1, len(sequence))
    ) % 2 else 1


# -------------------------------------------------------------------------
# Exact defining Whitney-form integral on a regular affine tetrahedron.
# -------------------------------------------------------------------------
coordinate_bases = [list(combinations(range(3), degree))
                    for degree in range(4)]


def wedge_components(covectors, degree):
    if degree == 0:
        return sy.Matrix((1,))
    return sy.Matrix([
        sy.det(sy.Matrix([
            [covector[index] for index in basis]
            for covector in covectors
        ]))
        for basis in coordinate_bases[degree]
    ])


def local_whitney_mass(points, degree):
    affine = sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    )
    inverse = affine.inv()
    gradients = [-sum(
        (sy.Matrix(inverse.row(row)).T for row in range(3)),
        sy.zeros(3, 1),
    )]
    gradients.extend(sy.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det()) / 6
    barycentric_second_moment = volume * (sy.ones(4, 4) + sy.eye(4)) / 20
    local_forms = list(combinations(range(4), degree + 1))
    coefficient_matrices = []
    for form in local_forms:
        coefficients = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            coefficients[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                components = (
                    factorial(degree) * (-1) ** omitted
                    * wedge_components(covectors, degree)
                )
                coefficients[:, form[omitted]] += components
        coefficient_matrices.append(coefficients)

    mass = sy.zeros(len(local_forms), len(local_forms))
    for row, left in enumerate(coefficient_matrices):
        for column, right in enumerate(coefficient_matrices):
            mass[row, column] = sy.simplify(sum(
                (left[basis, :] * barycentric_second_moment
                 * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
            ))
    return mass


def local_coboundaries():
    local_faces = [list(combinations(range(4), degree + 1))
                   for degree in range(4)]
    indices = [
        {simplex: index for index, simplex in enumerate(layer)}
        for layer in local_faces
    ]
    differentials = []
    for degree in range(3):
        matrix = sy.zeros(len(local_faces[degree + 1]),
                          len(local_faces[degree]))
        for row, simplex in enumerate(local_faces[degree + 1]):
            for omitted in range(degree + 2):
                face = simplex[:omitted] + simplex[omitted + 1:]
                matrix[row, indices[degree][face]] = (-1) ** omitted
        differentials.append(matrix)
    return local_faces, differentials


regular_points = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
local_faces, local_d = local_coboundaries()
local_masses = [local_whitney_mass(regular_points, degree)
                for degree in range(4)]
local_deltas = [
    sy.simplify(
        local_masses[degree].inv()
        * local_d[degree].T
        * local_masses[degree + 1]
    )
    for degree in range(3)
]
mass_integer = [exact_integer_matrix(matrix) for matrix in local_masses]
delta_integer = [exact_integer_matrix(matrix) for matrix in local_deltas]


def proportionality_ratio(left, right):
    ratio = None
    for row in range(left.rows):
        for column in range(left.cols):
            a = sy.simplify(left[row, column])
            b = sy.simplify(right[row, column])
            if b == 0:
                if a != 0:
                    return None
            else:
                candidate = sy.simplify(a / b)
                if ratio is None:
                    ratio = candidate
                elif sy.simplify(candidate - ratio) != 0:
                    return None
    return ratio


local_delta_ratios = [
    proportionality_ratio(local_deltas[degree], local_d[degree].T)
    for degree in range(3)
]


# -------------------------------------------------------------------------
# Complete 600-cell simplicial complex and global integer coboundaries.
# -------------------------------------------------------------------------
vertices, adjacency, _ = build_600cell()
neighbours = tuple(
    frozenset(np.flatnonzero(adjacency[index]).tolist())
    for index in range(120)
)
edges = tuple(
    (left, right)
    for left in range(120)
    for right in sorted(neighbours[left])
    if left < right
)
triangles = tuple(
    (left, right, third)
    for left, right in edges
    for third in sorted(neighbours[left] & neighbours[right])
    if right < third
)
tetrahedra = tuple(
    (first, second, third, fourth)
    for first, second, third in triangles
    for fourth in sorted(
        neighbours[first] & neighbours[second] & neighbours[third]
    )
    if third < fourth
)
cells = (
    tuple((index,) for index in range(120)),
    edges,
    triangles,
    tetrahedra,
)
dimensions = tuple(map(len, cells))
indices = [
    {simplex: index for index, simplex in enumerate(layer)}
    for layer in cells
]
global_d = []
for degree in range(3):
    rows = []
    columns = []
    values = []
    for row, simplex in enumerate(cells[degree + 1]):
        for omitted in range(degree + 2):
            face = simplex[:omitted] + simplex[omitted + 1:]
            rows.append(row)
            columns.append(indices[degree][face])
            values.append((-1) ** omitted)
    global_d.append(sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(dimensions[degree + 1], dimensions[degree]),
        dtype=np.int64,
    ))


# Restriction/duplication matrices.  The global cells and every tetrahedron
# are stored in increasing vertex order, so all restriction orientation signs
# are +1; compute them anyway rather than assuming that convention.
restrictions = []
occurrence_groups = []
occurrence_metadata = []
for degree in range(4):
    rows = []
    columns = []
    values = []
    metadata = []
    group = [[] for _ in range(dimensions[degree])]
    for tetrahedron_index, tetrahedron in enumerate(tetrahedra):
        for local_index, local_face in enumerate(local_faces[degree]):
            image = tuple(tetrahedron[index] for index in local_face)
            sorted_image = tuple(sorted(image))
            order = sorted(range(len(image)), key=image.__getitem__)
            sign = permutation_sign(order)
            global_index = indices[degree][sorted_image]
            row = tetrahedron_index * len(local_faces[degree]) + local_index
            rows.append(row)
            columns.append(global_index)
            values.append(sign)
            metadata.append((tetrahedron_index, local_index, global_index, sign))
            group[global_index].append(row)
    restriction = sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(len(tetrahedra) * len(local_faces[degree]), dimensions[degree]),
        dtype=np.int64,
    )
    restrictions.append(restriction)
    occurrence_groups.append(group)
    occurrence_metadata.append(metadata)


def equality_difference_matrix(groups, row_count):
    rows = []
    columns = []
    values = []
    metadata = []
    output_row = 0
    for global_index, occurrences in enumerate(groups):
        anchor = occurrences[0]
        for occurrence in occurrences[1:]:
            rows.extend((output_row, output_row))
            columns.extend((occurrence, anchor))
            values.extend((1, -1))
            metadata.append((global_index, occurrence, anchor))
            output_row += 1
    return sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(output_row, row_count),
        dtype=np.int64,
    ), metadata


differences = []
difference_metadata = []
for degree in range(4):
    matrix, metadata = equality_difference_matrix(
        occurrence_groups[degree], restrictions[degree].shape[0]
    )
    differences.append(matrix)
    difference_metadata.append(metadata)


def modular_rank(matrix, prime=1000003):
    """Sparse exact rank certificate over F_prime."""
    matrix = matrix.tocsr()
    pivots = {}
    for row_index in range(matrix.shape[0]):
        row = {
            int(column): int(value) % prime
            for column, value in zip(
                matrix.indices[matrix.indptr[row_index]:matrix.indptr[row_index + 1]],
                matrix.data[matrix.indptr[row_index]:matrix.indptr[row_index + 1]],
            )
            if int(value) % prime
        }
        while row:
            pivot = min(row)
            if pivot not in pivots:
                inverse = pow(row[pivot], -1, prime)
                row = {
                    column: (value * inverse) % prime
                    for column, value in row.items()
                    if (value * inverse) % prime
                }
                pivots[pivot] = row
                break
            factor = row[pivot]
            pivot_row = pivots[pivot]
            for column, value in pivot_row.items():
                updated = (row.get(column, 0) - factor * value) % prime
                if updated:
                    row[column] = updated
                elif column in row:
                    del row[column]
    return len(pivots)


print("=" * 78)
print("ELEMENT-LOCAL WHITNEY ASSEMBLY INVARIANCE")
print("=" * 78)

check("the global 600-cell f-vector is exact",
      dimensions == (120, 720, 1200, 600))
duplicated_dimensions = tuple(matrix.shape[0] for matrix in restrictions)
check("the element-local carrier has dimensions (2400,3600,2400,600)",
      duplicated_dimensions == (2400, 3600, 2400, 600)
      and sum(duplicated_dimensions) == 9000)
occurrence_histograms = []
for groups in occurrence_groups:
    values, counts = np.unique(list(map(len, groups)), return_counts=True)
    occurrence_histograms.append({int(value): int(count)
                                  for value, count in zip(values, counts)})
check("all simplex occurrence multiplicities are the exact regular values",
      occurrence_histograms == [
          {20: 120}, {5: 720}, {2: 1200}, {1: 600}
      ], str(occurrence_histograms))
check("the defining integral gives exact symmetric nonsingular local masses",
      all(matrix == matrix.T and matrix.det() != 0
          for matrix in local_masses))
check("on a regular element every local metric adjoint is incidence-proportional",
      tuple(local_delta_ratios) == (
          sy.Rational(5, 4), sy.Rational(5, 2), sy.Rational(15, 4)
      ), f"ratios={local_delta_ratios}")

upward_audits = []
leakage_audits = []
weak_audits = []
equality_glue_audits = []

assembled_mass_integer = []
for degree in range(4):
    local_mass_num, local_mass_denominator = mass_integer[degree]
    local_mass_global_num = sparse.kron(
        sparse.eye(len(tetrahedra), dtype=np.int64),
        sparse.csr_matrix(local_mass_num),
        format="csr",
    )
    assembled_mass_num = (
        restrictions[degree].T
        @ local_mass_global_num
        @ restrictions[degree]
    ).tocsr()
    assembled_mass_integer.append((assembled_mass_num, local_mass_denominator))

for degree in range(3):
    local_d_num = np.asarray(local_d[degree], dtype=np.int64)
    local_d_global = sparse.kron(
        sparse.eye(len(tetrahedra), dtype=np.int64),
        sparse.csr_matrix(local_d_num),
        format="csr",
    )
    upward_residual = (
        local_d_global @ restrictions[degree]
        - restrictions[degree + 1] @ global_d[degree]
    )
    upward_audits.append({
        "degree": degree,
        "residual_nonzeros": int(upward_residual.nnz),
    })

    delta_num, delta_denominator = delta_integer[degree]
    local_delta_global_num = sparse.kron(
        sparse.eye(len(tetrahedra), dtype=np.int64),
        sparse.csr_matrix(delta_num),
        format="csr",
    )
    leaked_num = (
        differences[degree]
        @ local_delta_global_num
        @ restrictions[degree + 1]
    ).tocsr()
    leaked_num.eliminate_zeros()
    s_rank = int(structural_rank(leaked_num))
    mod_rank = int(modular_rank(leaked_num))
    row, column = leaked_num.nonzero()
    witness = None
    if leaked_num.nnz:
        witness_row = int(row[0])
        witness_column = int(column[0])
        numerator = int(leaked_num[witness_row, witness_column])
        global_lower, occurrence, anchor = difference_metadata[degree][witness_row]
        witness = {
            "lower_global_simplex_index": global_lower,
            "lower_global_simplex": list(cells[degree][global_lower]),
            "higher_global_simplex_index": witness_column,
            "higher_global_simplex": list(cells[degree + 1][witness_column]),
            "compared_occurrence_rows": [occurrence, anchor],
            "exact_difference": f"{numerator}/{delta_denominator}",
        }
    off_incidence = int(sum(
        delta_num[row_index, column_index] != 0
        and local_d_num.T[row_index, column_index] == 0
        for row_index in range(delta_num.shape[0])
        for column_index in range(delta_num.shape[1])
    ))
    leakage_audits.append({
        "degree": degree,
        "local_delta_denominator": delta_denominator,
        "local_delta_off_incidence_nonzeros": off_incidence,
        "leakage_shape": list(leaked_num.shape),
        "leakage_exact_nonzeros": int(leaked_num.nnz),
        "leakage_max_abs_numerator": (
            int(np.max(np.abs(leaked_num.data))) if leaked_num.nnz else 0
        ),
        "structural_rank_upper_bound": s_rank,
        "modular_rank_lower_bound": mod_rank,
        "rank_certified_exact": bool(s_rank == mod_rank),
        "exact_rank": mod_rank if s_rank == mod_rank else None,
        "exact_witness": witness,
        "conforming_subspace_invariant": bool(leaked_num.nnz == 0),
    })

    # Clear the sole mass denominator.  Since M_p delta_p=d_p^T M_(p+1)
    # locally, both assembled sides use the numerator of M_(p+1).
    next_mass_num, next_mass_denominator = mass_integer[degree + 1]
    local_weak_num = sparse.kron(
        sparse.eye(len(tetrahedra), dtype=np.int64),
        sparse.csr_matrix(local_d_num.T @ next_mass_num),
        format="csr",
    )
    assembled_next_mass_num = assembled_mass_integer[degree + 1][0]
    weak_left = (
        restrictions[degree].T
        @ local_weak_num
        @ restrictions[degree + 1]
    )
    weak_right = global_d[degree].T @ assembled_next_mass_num
    weak_residual = (weak_left - weak_right).tocsr()
    weak_residual.eliminate_zeros()
    weak_audits.append({
        "degree": degree,
        "cleared_mass_denominator": next_mass_denominator,
        "residual_nonzeros": int(weak_residual.nnz),
        "assembled_mass_nonzeros": int(assembled_next_mass_num.nnz),
    })

    # POST-PROTOCOL HOSTILE AUDIT.  The Euclidean equality projector is a
    # product of small per-simplex Grover averages and is therefore the
    # cheapest reversible glue candidate.  Its compressed downward map is
    # X_eq=(J^T J)^-1 J^T delta_loc J.  Test exact equality with the global
    # Whitney adjoint without forming an inverse: M_p X_eq=d^T M_(p+1).
    current_mass_num, current_mass_denominator = assembled_mass_integer[degree]
    multiplicities = np.asarray(
        restrictions[degree].power(2).sum(axis=0)
    ).ravel().astype(np.int64)
    if len(set(multiplicities.tolist())) != 1:
        raise RuntimeError("base 600-cell degree multiplicities are not uniform")
    multiplicity = int(multiplicities[0])
    equality_sum_num = (
        restrictions[degree].T
        @ local_delta_global_num
        @ restrictions[degree + 1]
    )
    equality_left = (
        current_mass_num @ equality_sum_num * next_mass_denominator
    )
    equality_right = (
        global_d[degree].T @ assembled_next_mass_num
        * (current_mass_denominator * delta_denominator * multiplicity)
    )
    equality_residual = (equality_left - equality_right).tocsr()
    equality_residual.eliminate_zeros()
    witness = None
    if equality_residual.nnz:
        witness_row, witness_column = equality_residual.nonzero()
        row = int(witness_row[0])
        column = int(witness_column[0])
        witness = {
            "lower_global_simplex": list(cells[degree][row]),
            "higher_global_simplex": list(cells[degree + 1][column]),
            "cleared_integer_residual": int(equality_residual[row, column]),
        }
    equality_glue_audits.append({
        "degree": degree,
        "uniform_copy_multiplicity": multiplicity,
        "metric_identity_residual_nonzeros": int(equality_residual.nnz),
        "metric_identity_max_abs_integer": (
            int(np.max(np.abs(equality_residual.data)))
            if equality_residual.nnz else 0
        ),
        "exact_witness": witness,
        "equals_global_whitney_adjoint": bool(equality_residual.nnz == 0),
    })

check("all three upward restrictions intertwine exactly",
      all(audit["residual_nonzeros"] == 0 for audit in upward_audits),
      str(upward_audits))
leakage_count = sum(not audit["conforming_subspace_invariant"]
                    for audit in leakage_audits)
check("all three element-local codifferentials leak out of conformity exactly",
      leakage_count == 3
      and all(audit["exact_witness"] is not None
              for audit in leakage_audits),
      f"leakage degrees={leakage_count}/3")
check("all three leakage ranks are certified by matching exact bounds",
      all(audit["rank_certified_exact"] for audit in leakage_audits),
      str([(audit["exact_rank"], audit["leakage_exact_nonzeros"])
           for audit in leakage_audits]))
check("all three metric-projected weak assembly identities hold exactly",
      all(audit["residual_nonzeros"] == 0 for audit in weak_audits),
      str(weak_audits))
check("simple equality/Grover averaging fails the Whitney metric in all degrees",
      all(not audit["equals_global_whitney_adjoint"]
          for audit in equality_glue_audits),
      str([(audit["degree"], audit["metric_identity_residual_nonzeros"])
           for audit in equality_glue_audits]))

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "global_dimensions": list(dimensions),
    "duplicated_dimensions": list(duplicated_dimensions),
    "duplicated_total_dimension": sum(duplicated_dimensions),
    "occurrence_histograms": occurrence_histograms,
    "local_mass_matrices": [
        [[str(matrix[row, column]) for column in range(matrix.cols)]
         for row in range(matrix.rows)]
        for matrix in local_masses
    ],
    "local_codifferential_matrices": [
        [[str(matrix[row, column]) for column in range(matrix.cols)]
         for row in range(matrix.rows)]
        for matrix in local_deltas
    ],
    "local_codifferential_incidence_ratios": list(map(str, local_delta_ratios)),
    "upward_intertwining": upward_audits,
    "downward_assembly_leakage": leakage_audits,
    "invariant_degree_hit_fraction": [3 - leakage_count, 3],
    "weak_metric_projection_identity": weak_audits,
    "post_protocol_euclidean_equality_glue_audit": equality_glue_audits,
    "verdicts": [
        "DERIVED ASSEMBLY-LEAKAGE NO-GO: 3/3 downward degrees leak",
        "DERIVED PROJECTED FACTORIZATION: 3/3 weak identities hold exactly",
        "DERIVED NEGATIVE: simple equality/Grover glue misses Whitney in 3/3 degrees",
    ],
    "scope": (
        "Direct-sum regular-element Whitney generator on the base 600-cell. "
        "It does not exclude reversible multistep glue dynamics or enlarged "
        "discontinuous-Galerkin flux carriers."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured element-assembly certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for audit in leakage_audits:
    print(
        "degree {degree}: leakage nnz={leakage_exact_nonzeros}, "
        "rank={exact_rank}, witness={exact_witness}".format(**audit)
    )
print("INVARIANT_DOWNWARD_DEGREES=0/3")
print("PROJECTED_WEAK_IDENTITIES=3/3")
print("VERDICT: exact metric glue is necessary and contains the unresolved step.")
raise SystemExit(0 if passed == tests else 1)

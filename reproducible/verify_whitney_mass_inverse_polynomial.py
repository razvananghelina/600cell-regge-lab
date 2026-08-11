#!/usr/bin/env python3
"""Exact minimal-polynomial and locality audit for assembled Whitney masses.

Protocol commit 3323174 froze the four matrices, exact modular certificates,
inverse-polynomial construction and locality diagnostics before computation.
"""

from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sparse
from scipy.sparse.csgraph import connected_components, shortest_path
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_mass_inverse_polynomial.json")
PROTOCOL_COMMIT = "3323174"
PROBE_SEED = 600_20260811
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def permutation_sign(sequence):
    return -1 if sum(
        sequence[left] > sequence[right]
        for left in range(len(sequence))
        for right in range(left + 1, len(sequence))
    ) % 2 else 1


def exact_integer_matrix(matrix):
    denominators = [sy.denom(value) for value in matrix]
    denominator = sy.ilcm(*denominators) if len(denominators) > 1 else (
        denominators[0] if denominators else sy.Integer(1)
    )
    numerator = np.asarray([
        [int(sy.expand(denominator * matrix[row, column]))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.int64)
    return numerator, int(denominator)


# -------------------------------------------------------------------------
# Derive the local Whitney masses from their defining affine-form integral.
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
                coefficients[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * wedge_components(covectors, degree)
                )
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


regular_points = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
local_faces = [list(combinations(range(4), degree + 1))
               for degree in range(4)]
local_masses = [local_whitney_mass(regular_points, degree)
                for degree in range(4)]
local_mass_integer = [exact_integer_matrix(matrix) for matrix in local_masses]


# -------------------------------------------------------------------------
# Assemble the exact 600-cell mass blocks.
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

restrictions = []
for degree in range(4):
    rows = []
    columns = []
    values = []
    for tetrahedron_index, tetrahedron in enumerate(tetrahedra):
        for local_index, local_face in enumerate(local_faces[degree]):
            image = tuple(tetrahedron[index] for index in local_face)
            sorted_image = tuple(sorted(image))
            order = sorted(range(len(image)), key=image.__getitem__)
            rows.append(tetrahedron_index * len(local_faces[degree]) + local_index)
            columns.append(indices[degree][sorted_image])
            values.append(permutation_sign(order))
    restrictions.append(sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(len(tetrahedra) * len(local_faces[degree]), dimensions[degree]),
        dtype=np.int64,
    ))


primitive_masses = []
assembly_audits = []
for degree in range(4):
    local_numerator, local_denominator = local_mass_integer[degree]
    duplicated = sparse.kron(
        sparse.eye(len(tetrahedra), dtype=np.int64),
        sparse.csr_matrix(local_numerator),
        format="csr",
    )
    assembled = (
        restrictions[degree].T @ duplicated @ restrictions[degree]
    ).tocsr()
    common_divisor = 0
    for value in assembled.data:
        common_divisor = gcd(common_divisor, abs(int(value)))
    primitive = assembled.copy().astype(np.int64).tocsr()
    primitive.data //= common_divisor
    primitive.eliminate_zeros()
    primitive_masses.append(primitive)
    assembly_audits.append({
        "degree": degree,
        "dimension": dimensions[degree],
        "local_denominator": local_denominator,
        "assembled_numerator_gcd": common_divisor,
        "primitive_nonzeros": int(primitive.nnz),
        "symmetric": bool((primitive - primitive.T).nnz == 0),
    })


# -------------------------------------------------------------------------
# Exact modular polynomial machinery.
# -------------------------------------------------------------------------
def berlekamp_massey(sequence, prime):
    """Connection polynomial C with sum_i C_i s_(n-i)=0 modulo prime."""
    connection = [1]
    previous = [1]
    length = 0
    shift = 1
    previous_discrepancy = 1
    for index in range(len(sequence)):
        discrepancy = sequence[index] % prime
        for offset in range(1, length + 1):
            discrepancy = (
                discrepancy + connection[offset] * sequence[index - offset]
            ) % prime
        if discrepancy == 0:
            shift += 1
            continue
        old_connection = connection[:]
        scale = discrepancy * pow(previous_discrepancy, -1, prime) % prime
        required = len(previous) + shift
        if len(connection) < required:
            connection.extend([0] * (required - len(connection)))
        for offset, coefficient in enumerate(previous):
            connection[offset + shift] = (
                connection[offset + shift] - scale * coefficient
            ) % prime
        if 2 * length <= index:
            length = index + 1 - length
            previous = old_connection
            previous_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    connection = connection[:length + 1]
    # x^L + C1 x^(L-1) + ... + CL, returned low coefficient first.
    return list(reversed(connection))


def scalar_krylov_polynomial(matrix, probe, prime):
    vector = np.asarray(probe % prime, dtype=np.int64)
    initial = vector.copy()
    sequence = []
    for _ in range(2 * matrix.shape[0] + 16):
        sequence.append(int(initial @ vector % prime))
        vector = np.asarray(matrix @ vector, dtype=np.int64) % prime
    polynomial = berlekamp_massey(sequence, prime)
    return polynomial, len(polynomial) - 1


def centered_crt(residues, primes):
    value = 0
    modulus = 1
    for residue, prime in zip(residues, primes):
        step = ((int(residue) - value) % prime) * pow(modulus, -1, prime) % prime
        value += modulus * step
        modulus *= prime
    if value > modulus // 2:
        value -= modulus
    return value, modulus


def evaluate_polynomial_mod(matrix, coefficients, prime):
    dimension = matrix.shape[0]
    result = np.zeros((dimension, dimension), dtype=np.int64)
    diagonal = np.diag_indices(dimension)
    result[diagonal] = int(coefficients[-1] % prime)
    for coefficient in reversed(coefficients[:-1]):
        result = np.asarray(matrix @ result, dtype=np.int64) % prime
        result[diagonal] = (result[diagonal] + int(coefficient % prime)) % prime
    return result


def polynomial_bound(coefficients, row_sum_bound):
    return sum(
        abs(int(coefficient)) * int(row_sum_bound) ** power
        for power, coefficient in enumerate(coefficients)
    )


def next_recorded_prime(start):
    return int(sy.nextprime(start))


def exact_polynomial_audit(matrix, degree):
    dimension = matrix.shape[0]
    row_sum_bound = int(np.max(np.asarray(abs(matrix).sum(axis=1)).ravel()))
    rng = np.random.default_rng(PROBE_SEED + degree)
    probe = rng.integers(1, 1_000_000, size=dimension, dtype=np.int64)

    records = []
    prime = 1_000_000
    for _ in range(12):
        prime = next_recorded_prime(prime)
        polynomial, complexity = scalar_krylov_polynomial(matrix, probe, prime)
        records.append((prime, polynomial, complexity))

    maximum_complexity = max(record[2] for record in records)
    selected = [record for record in records
                if record[2] == maximum_complexity]
    coefficient_bound = (1 + row_sum_bound) ** maximum_complexity
    reconstruction_modulus = 1
    for prime, _polynomial, _complexity in selected:
        reconstruction_modulus *= prime
    while reconstruction_modulus <= 2 * coefficient_bound:
        prime = next_recorded_prime(prime)
        polynomial, complexity = scalar_krylov_polynomial(matrix, probe, prime)
        records.append((prime, polynomial, complexity))
        if complexity > maximum_complexity:
            maximum_complexity = complexity
            selected = [(prime, polynomial, complexity)]
            coefficient_bound = (1 + row_sum_bound) ** maximum_complexity
            reconstruction_modulus = prime
        elif complexity == maximum_complexity:
            selected.append((prime, polynomial, complexity))
            reconstruction_modulus *= prime

    selected_primes = [record[0] for record in selected]
    coefficients = []
    for index in range(maximum_complexity + 1):
        value, modulus = centered_crt(
            [record[1][index] for record in selected], selected_primes
        )
        coefficients.append(int(value))
    reconstruction_modulus = modulus
    coefficients_within_bound = all(
        abs(value) <= coefficient_bound for value in coefficients
    )

    residual_bound = polynomial_bound(coefficients, row_sum_bound)
    inverse_coefficients = [-coefficients[index]
                            for index in range(1, len(coefficients))]
    inverse_bound = polynomial_bound(inverse_coefficients, row_sum_bound)
    required_bound = max(residual_bound, inverse_bound)
    validation_modulus = 1
    validation_primes = []
    annihilates_whole_matrix = True
    inverse_nonzero_support = np.zeros((dimension, dimension), dtype=bool)
    validation_prime = 2_000_000
    while validation_modulus <= 2 * required_bound or len(validation_primes) < 2:
        validation_prime = next_recorded_prime(validation_prime)
        validation_primes.append(validation_prime)
        polynomial_value = evaluate_polynomial_mod(
            matrix, coefficients, validation_prime
        )
        annihilates_whole_matrix &= not np.any(polynomial_value)
        inverse_value = evaluate_polynomial_mod(
            matrix, inverse_coefficients, validation_prime
        )
        inverse_nonzero_support |= inverse_value != 0
        validation_modulus *= validation_prime

    support_graph = matrix.copy().tolil()
    support_graph.setdiag(0)
    support_graph = support_graph.tocsr()
    support_graph.eliminate_zeros()
    undirected = ((support_graph != 0) + (support_graph.T != 0)).astype(np.int8)
    component_count, component_labels = connected_components(
        undirected, directed=False
    )
    distances = shortest_path(undirected, directed=False, unweighted=True)
    finite_distances = distances[np.isfinite(distances)]
    diameter = int(np.max(finite_distances)) if finite_distances.size else 0
    inverse_distances = distances[
        inverse_nonzero_support & np.isfinite(distances)
    ]
    maximum_inverse_distance = (
        int(np.max(inverse_distances)) if inverse_distances.size else 0
    )
    cross_component_nonzeros = int(sum(
        inverse_nonzero_support[row, column]
        and component_labels[row] != component_labels[column]
        for row in range(dimension)
        for column in range(dimension)
    ))
    constant_nonzero = coefficients[0] != 0
    minimality_certified = bool(
        annihilates_whole_matrix
        and coefficients[-1] == 1
        and coefficients_within_bound
    )
    variable = sy.symbols("x")
    exact_polynomial = sum(
        sy.Integer(coefficient) * variable ** power
        for power, coefficient in enumerate(coefficients)
    )
    return {
        "degree": degree,
        "dimension": dimension,
        "primitive_nonzeros": int(matrix.nnz),
        "maximum_absolute_row_sum": row_sum_bound,
        "probe_seed": PROBE_SEED + degree,
        "probe_prime_records": [
            {"prime": record[0], "linear_complexity": record[2]}
            for record in records
        ],
        "selected_reconstruction_primes": selected_primes,
        "reconstruction_modulus": str(reconstruction_modulus),
        "coefficient_bound": str(coefficient_bound),
        "minimal_polynomial_coefficients_low_to_high": coefficients,
        "minimal_polynomial_factorization_over_Z": str(
            sy.factor(exact_polynomial)
        ),
        "minimal_polynomial_degree": maximum_complexity,
        "whole_matrix_residual_bound": str(residual_bound),
        "validation_primes": validation_primes,
        "validation_modulus": str(validation_modulus),
        "whole_matrix_annihilation_certified": annihilates_whole_matrix,
        "minimality_certified": minimality_certified,
        "constant_coefficient_nonzero": constant_nonzero,
        "inverse_polynomial_coefficients_low_to_high_numerator": (
            inverse_coefficients
        ),
        "inverse_polynomial_denominator": coefficients[0],
        "inverse_polynomial_degree": maximum_complexity - 1,
        "inverse_entry_bound_numerator": str(inverse_bound),
        "inverse_support_certified_exact": bool(
            validation_modulus > 2 * inverse_bound
        ),
        "inverse_exact_nonzeros": int(inverse_nonzero_support.sum()),
        "inverse_nonzero_fraction": float(
            inverse_nonzero_support.mean()
        ),
        "support_graph_components": int(component_count),
        "support_graph_diameter": diameter,
        "maximum_nonzero_inverse_distance": maximum_inverse_distance,
        "inverse_reaches_component_diameter": bool(
            maximum_inverse_distance == diameter
        ),
        "cross_component_inverse_nonzeros": cross_component_nonzeros,
        "inverse_is_one_step_local": bool(maximum_inverse_distance <= 1),
        "cayley_hamilton_degree_ratio": (
            float((maximum_complexity - 1) / (dimension - 1))
            if dimension > 1 else 0.0
        ),
    }


print("=" * 78)
print("EXACT WHITNEY MASS-INVERSE POLYNOMIAL LOCALITY")
print("=" * 78)

check("the independently rebuilt f-vector is exact",
      dimensions == (120, 720, 1200, 600))
check("all four defining local masses are symmetric and nonsingular",
      all(matrix == matrix.T and matrix.det() != 0
          for matrix in local_masses))
check("all four primitive assembled masses are symmetric",
      all(audit["symmetric"] for audit in assembly_audits),
      str([(audit["dimension"], audit["primitive_nonzeros"])
           for audit in assembly_audits]))

mass_laplacian_models = (
    20 * sparse.eye(dimensions[0], dtype=np.int64)
    - global_d[0].T @ global_d[0],
    50 * sparse.eye(dimensions[1], dtype=np.int64)
    - 3 * global_d[1].T @ global_d[1],
    20 * sparse.eye(dimensions[2], dtype=np.int64)
    - global_d[2].T @ global_d[2],
    sparse.eye(dimensions[3], dtype=np.int64),
)
mass_laplacian_residual_nonzeros = [
    int((primitive_masses[degree] - mass_laplacian_models[degree]).nnz)
    for degree in range(4)
]
check("all primitive masses are exact affine polynomials of upper Laplacians",
      mass_laplacian_residual_nonzeros == [0, 0, 0, 0],
      f"residual nnz={mass_laplacian_residual_nonzeros}")

audits = []
for degree, matrix in enumerate(primitive_masses):
    print(f"  auditing degree {degree}, dimension {matrix.shape[0]} ...")
    audits.append(exact_polynomial_audit(matrix, degree))

check("all four whole-matrix annihilation certificates are exact",
      all(audit["whole_matrix_annihilation_certified"]
          for audit in audits))
check("all four minimal-polynomial degrees have matching lower and upper bounds",
      all(audit["minimality_certified"] for audit in audits),
      str([audit["minimal_polynomial_degree"] for audit in audits]))
check("all four primitive masses are exactly invertible",
      all(audit["constant_coefficient_nonzero"] for audit in audits))
check("all four inverse support censuses are exact",
      all(audit["inverse_support_certified_exact"] for audit in audits))
check("no inverse couples distinct support-graph components",
      all(audit["cross_component_inverse_nonzeros"] == 0
          for audit in audits))

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "f_vector": list(dimensions),
    "assembly_audits": assembly_audits,
    "primitive_mass_upper_laplacian_identities": [
        "B0 = 20 I - d0^T d0",
        "B1 = 50 I - 3 d1^T d1",
        "B2 = 20 I - d2^T d2",
        "B3 = I",
    ],
    "primitive_mass_upper_laplacian_residual_nonzeros": (
        mass_laplacian_residual_nonzeros
    ),
    "mass_inverse_audits": audits,
    "verdicts": [
        "DERIVED exact minimal and inverse polynomials on the fixed complex",
        "Locality verdict is determined degree by degree from exact inverse support",
        "No bounded-depth refinement claim is made",
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured mass-polynomial certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for audit in audits:
    print(
        "p={degree}: N={dimension}, minpoly={minimal_polynomial_degree}, "
        "invdeg={inverse_polynomial_degree}, components={support_graph_components}, "
        "diam={support_graph_diameter}, invdist={maximum_nonzero_inverse_distance}, "
        "invnnz={inverse_exact_nonzeros}, fraction={inverse_nonzero_fraction:.6f}, "
        "one_step={inverse_is_one_step_local}".format(**audit)
    )
raise SystemExit(0 if passed == tests else 1)

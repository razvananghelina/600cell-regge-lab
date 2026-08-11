#!/usr/bin/env python3
"""Exact first-refinement lower bounds for Whitney inverse complexity.

Protocol commit 366fe4a froze the barycentric carrier, exact Whitney masses,
modular probes, primes, sequence lengths and decision labels before values
were computed.
"""

from itertools import combinations, permutations
import json
from math import factorial, gcd
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sparse
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_mass_inverse_refinement.json")
PROTOCOL_COMMIT = "366fe4a"
PROBE_SEED = 600_20260811
PRIMES = (1000003, 1000033, 1000037)
COARSE_MINIMAL_DEGREES = (9, 22, 27, 1)
SEQUENCE_LENGTHS = tuple(4 * degree + 32 for degree in COARSE_MINIMAL_DEGREES)
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
    denominator = sy.ilcm(*denominators) if len(denominators) > 1 else (
        denominators[0] if denominators else sy.Integer(1)
    )
    numerator = np.asarray([
        [int(sy.expand(denominator * matrix[row, column]))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.int64)
    return numerator, int(denominator)


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


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


def assemble_primitive_mass(top_cells, cells, form_degree, local_mass):
    local_faces = list(combinations(range(4), form_degree + 1))
    cell_indices = {cell: index
                    for index, cell in enumerate(cells[form_degree])}
    numerator, denominator = exact_integer_matrix(local_mass)
    rows = []
    columns = []
    data = []
    for top in top_cells:
        global_faces = [tuple(top[index] for index in face)
                        for face in local_faces]
        global_indices = [cell_indices[face] for face in global_faces]
        for local_row, global_row in enumerate(global_indices):
            for local_column, global_column in enumerate(global_indices):
                value = int(numerator[local_row, local_column])
                if value:
                    rows.append(global_row)
                    columns.append(global_column)
                    data.append(value)
    assembled = sparse.coo_matrix(
        (data, (rows, columns)),
        shape=(len(cells[form_degree]), len(cells[form_degree])),
        dtype=np.int64,
    ).tocsr()
    assembled.sum_duplicates()
    common_divisor = 0
    for value in assembled.data:
        common_divisor = gcd(common_divisor, abs(int(value)))
    primitive = assembled.copy()
    primitive.data //= common_divisor
    primitive.eliminate_zeros()
    return primitive, {
        "local_denominator": denominator,
        "assembled_numerator_gcd": common_divisor,
        "dimension": primitive.shape[0],
        "nonzeros": int(primitive.nnz),
        "maximum_absolute_entry": int(np.max(np.abs(primitive.data))),
        "maximum_absolute_row_sum": int(np.max(
            np.asarray(abs(primitive).sum(axis=1)).ravel()
        )),
    }


def berlekamp_massey(sequence, prime):
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
    return length


def krylov_complexity(matrix, form_degree, prime, sequence_length):
    rng = np.random.default_rng(PROBE_SEED + form_degree)
    probe = rng.integers(1, 1_000_000, size=matrix.shape[0], dtype=np.int64)
    initial = np.asarray(probe % prime, dtype=np.int64)
    vector = initial.copy()
    sequence = []
    for _ in range(sequence_length):
        sequence.append(int(initial @ vector % prime))
        vector = np.asarray(matrix @ vector, dtype=np.int64) % prime
    return berlekamp_massey(sequence, prime)


print("=" * 78)
print("WHITNEY INVERSE COMPLEXITY UNDER FIRST BARYCENTRIC REFINEMENT")
print("=" * 78)

vertices, adjacency, _ = build_600cell()
neighbours = tuple(
    frozenset(np.flatnonzero(adjacency[index]).tolist())
    for index in range(120)
)
coarse_edges = tuple(
    (left, right)
    for left in range(120)
    for right in sorted(neighbours[left])
    if left < right
)
coarse_triangles = tuple(
    (left, right, third)
    for left, right in coarse_edges
    for third in sorted(neighbours[left] & neighbours[right])
    if right < third
)
coarse_tetrahedra = tuple(
    (first, second, third, fourth)
    for first, second, third in coarse_triangles
    for fourth in sorted(
        neighbours[first] & neighbours[second] & neighbours[third]
    )
    if third < fourth
)
coarse_top = coarse_tetrahedra
coarse_cells = (
    tuple((index,) for index in range(120)),
    coarse_edges,
    coarse_triangles,
    coarse_tetrahedra,
)
coarse_dimensions = tuple(map(len, coarse_cells))
check("the coarse f-vector is exact",
      coarse_dimensions == (120, 720, 1200, 600))

fine_vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
fine_vertex_index = {cell: index for index, cell in enumerate(fine_vertex_cells)}
fine_top = []
for tetrahedron in coarse_tetrahedra:
    for ordering in permutations(tetrahedron):
        flag = (
            (ordering[0],),
            tuple(sorted(ordering[:2])),
            tuple(sorted(ordering[:3])),
            tetrahedron,
        )
        fine_top.append(tuple(fine_vertex_index[cell] for cell in flag))
fine_top = tuple(fine_top)
fine_cells = all_simplices(fine_top)
fine_dimensions = tuple(map(len, fine_cells))
check("the complete fine f-vector and top count are exact",
      fine_dimensions == (2640, 17040, 28800, 14400)
      and len(fine_top) == 14400,
      f"fine f-vector={fine_dimensions}")

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
coarse_local_masses = [
    local_whitney_mass(reference_vertices, degree) for degree in range(4)
]

child_mass_families = [[] for _ in range(4)]
for ordering in permutations(range(4)):
    child = (
        reference_vertices[ordering[0]],
        sum((reference_vertices[index] for index in ordering[:2]),
            sy.zeros(3, 1)) / 2,
        sum((reference_vertices[index] for index in ordering[:3]),
            sy.zeros(3, 1)) / 3,
        sum(reference_vertices, sy.zeros(3, 1)) / 4,
    )
    for degree in range(4):
        child_mass_families[degree].append(
            local_whitney_mass(child, degree)
        )
child_local_masses = [family[0] for family in child_mass_families]
check("all 24 flag children have identical rank-ordered local masses",
      all(matrix == family[0]
          for family in child_mass_families for matrix in family))
check("all coarse and child local masses are exact nonsingular Gram blocks",
      all(matrix == matrix.T and matrix.det() != 0
          for matrix in coarse_local_masses + child_local_masses))

coarse_masses = []
fine_masses = []
coarse_assembly = []
fine_assembly = []
for degree in range(4):
    coarse_mass, coarse_audit = assemble_primitive_mass(
        coarse_top, coarse_cells, degree, coarse_local_masses[degree]
    )
    fine_mass, fine_audit = assemble_primitive_mass(
        fine_top, fine_cells, degree, child_local_masses[degree]
    )
    coarse_masses.append(coarse_mass)
    fine_masses.append(fine_mass)
    coarse_assembly.append(coarse_audit)
    fine_assembly.append(fine_audit)

check("all eight primitive assembled masses are exactly symmetric",
      all((matrix - matrix.T).nnz == 0
          for matrix in coarse_masses + fine_masses))
check("the refined primitive top-form mass remains the identity",
      (fine_masses[3] - sparse.eye(
          fine_dimensions[3], dtype=np.int64
      )).nnz == 0)
check("all three recorded moduli are prime",
      all(sy.isprime(prime) for prime in PRIMES), str(PRIMES))

coarse_complexities = []
fine_complexities = []
for degree in range(4):
    print(
        f"  degree {degree}: coarse N={coarse_dimensions[degree]}, "
        f"fine N={fine_dimensions[degree]}, "
        f"sequence={SEQUENCE_LENGTHS[degree]}"
    )
    coarse_complexities.append([
        krylov_complexity(
            coarse_masses[degree], degree, prime, SEQUENCE_LENGTHS[degree]
        )
        for prime in PRIMES
    ])
    fine_complexities.append([
        krylov_complexity(
            fine_masses[degree], degree, prime, SEQUENCE_LENGTHS[degree]
        )
        for prime in PRIMES
    ])

check("the modular estimator recovers every exact coarse degree",
      all(complexity == COARSE_MINIMAL_DEGREES[degree]
          for degree in range(4)
          for complexity in coarse_complexities[degree]),
      f"coarse complexities={coarse_complexities}")
check("the top-form complexity remains one at refinement",
      fine_complexities[3] == [1, 1, 1],
      f"top complexities={fine_complexities[3]}")

fine_lower_bounds = tuple(map(max, fine_complexities))
growth = tuple(
    fine_lower_bounds[degree] > COARSE_MINIMAL_DEGREES[degree]
    for degree in range(4)
)
records = []
for degree in range(4):
    status = (
        "DERIVED DEGREE GROWTH AT LEVEL 1"
        if growth[degree]
        else "DERIVED NEGATIVE FOR GROWTH DETECTION"
    )
    records.append({
        "degree": degree,
        "coarse_dimension": coarse_dimensions[degree],
        "fine_dimension": fine_dimensions[degree],
        "sequence_length": SEQUENCE_LENGTHS[degree],
        "coarse_exact_minimal_degree": COARSE_MINIMAL_DEGREES[degree],
        "coarse_complexities_by_prime": coarse_complexities[degree],
        "fine_complexities_by_prime": fine_complexities[degree],
        "fine_certified_minimal_degree_lower_bound": fine_lower_bounds[degree],
        "fine_inverse_degree_lower_bound": fine_lower_bounds[degree] - 1,
        "lower_bound_ratio_to_coarse_exact_degree": float(
            fine_lower_bounds[degree] / COARSE_MINIMAL_DEGREES[degree]
        ),
        "growth_detected": growth[degree],
        "status": status,
    })

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "probe_seed_rule": "60020260811 + form_degree",
    "primes": list(PRIMES),
    "coarse_f_vector": list(coarse_dimensions),
    "fine_f_vector": list(fine_dimensions),
    "coarse_assembly": coarse_assembly,
    "fine_assembly": fine_assembly,
    "degree_records": records,
    "verdict": (
        "DERIVED LEVEL-1 GROWTH in degrees "
        + str([degree for degree in range(4) if growth[degree]])
    ),
    "scope": (
        "Exact deterministic lower bounds at one barycentric refinement; "
        "not full fine minimal polynomials and not a continuum scaling law."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured refinement lower-bound certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for record in records:
    print(
        "p={degree}: coarse exact={coarse_exact_minimal_degree}, "
        "fine complexities={fine_complexities_by_prime}, "
        "fine lower bound={fine_certified_minimal_degree_lower_bound}, "
        "ratio={lower_bound_ratio_to_coarse_exact_degree:.6f}, "
        "{status}".format(**record)
    )
raise SystemExit(0 if passed == tests else 1)

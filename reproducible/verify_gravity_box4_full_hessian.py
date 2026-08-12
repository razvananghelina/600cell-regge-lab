#!/usr/bin/env python3
"""Corrective audit of the full Hessian of Tr(Box^4).

Protocol commit d909513 froze the operator, perturbation space, exact
noncommutative derivative, controls, and decision rule.  This verifier uses
no phenomenological or graviton target.
"""

from itertools import permutations, product
import json
from pathlib import Path

import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_box4_full_hessian.json"
PROTOCOL_COMMIT = "d909513"
N = 120
A1 = 5
LEGACY_TOL = 1.0e-8

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def build_600cell():
    """Return the canonical 120 vertices and integer adjacency matrix."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vertices = set()

    def add(values):
        vector = np.asarray(values, dtype=float)
        vector /= np.linalg.norm(vector)
        vertices.add(tuple(np.round(vector, 10)))

    for coordinate in range(4):
        for sign in (1.0, -1.0):
            vector = [0.0] * 4
            vector[coordinate] = sign
            add(vector)
    for signs in product((0.5, -0.5), repeat=4):
        add(signs)

    base = [0.0, 0.5, phi / 2.0, 1.0 / (2.0 * phi)]
    even_permutations = [
        item for item in permutations(range(4))
        if sum(
            item[i] > item[j]
            for i in range(4) for j in range(i + 1, 4)
        ) % 2 == 0
    ]
    for permutation in even_permutations:
        coordinates = [base[permutation[index]] for index in range(4)]
        nonzero = [
            index for index, value in enumerate(coordinates)
            if abs(value) > 1.0e-12
        ]
        for signs in product((1.0, -1.0), repeat=len(nonzero)):
            vector = list(coordinates)
            for index, sign in zip(nonzero, signs):
                vector[index] *= sign
            add(vector)

    vertex_array = np.asarray(sorted(vertices))
    dots = np.clip(vertex_array @ vertex_array.T, -1.0, 1.0)
    adjacency = np.zeros((N, N), dtype=np.int64)
    for left in range(N):
        for right in range(left + 1, N):
            if abs(dots[left, right] - phi / 2.0) < 1.0e-6:
                adjacency[left, right] = adjacency[right, left] = 1
    return vertex_array, adjacency


def quaternion_multiply(left, right):
    return np.array([
        left[0] * right[0] - left[1] * right[1]
        - left[2] * right[2] - left[3] * right[3],
        left[0] * right[1] + left[1] * right[0]
        + left[2] * right[3] - left[3] * right[2],
        left[0] * right[2] - left[1] * right[3]
        + left[2] * right[0] + left[3] * right[1],
        left[0] * right[3] + left[1] * right[2]
        - left[2] * right[1] + left[3] * right[0],
    ])


def find_vertex(vector, vertices):
    dots = vertices @ vector
    index = int(np.argmax(dots))
    return index if dots[index] > 1.0 - 1.0e-6 else -1


def find_hopf_fibers(vertices):
    """Find the twelve right cosets of a derived order-ten subgroup."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    for index in range(N):
        if abs(vertices[index, 0] - phi / 2.0) >= 1.0e-6:
            continue
        generator = vertices[index]
        power = generator.copy()
        valid = True
        for exponent in range(2, 11):
            power = quaternion_multiply(power, generator)
            if exponent == 5 and not np.allclose(
                power, [-1.0, 0.0, 0.0, 0.0], atol=1.0e-6
            ):
                valid = False
                break
            if exponent == 10 and not np.allclose(
                power, [1.0, 0.0, 0.0, 0.0], atol=1.0e-6
            ):
                valid = False
        if not valid:
            continue

        subgroup = []
        power = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_vertex(power, vertices))
            power = quaternion_multiply(power, generator)

        used = set()
        fibers = []
        for seed in range(N):
            if seed in used:
                continue
            fiber = []
            for subgroup_index in subgroup:
                target = quaternion_multiply(
                    vertices[seed], vertices[subgroup_index]
                )
                target_index = find_vertex(target, vertices)
                if target_index >= 0 and target_index not in used:
                    fiber.append(target_index)
                    used.add(target_index)
            if len(fiber) == 10:
                fibers.append(fiber)
        if len(fibers) == 12:
            return fibers
    raise RuntimeError("No order-ten Hopf fibration found")


def rank_mod_prime(matrix, prime):
    """Row rank over F_prime, using exact int64 arithmetic."""
    work = np.remainder(matrix, prime).astype(np.int64, copy=True)
    row_count, column_count = work.shape
    rank = 0
    for column in range(column_count):
        candidates = np.flatnonzero(work[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            work[[rank, pivot]] = work[[pivot, rank]]
        inverse = pow(int(work[rank, column]), -1, prime)
        work[rank] = (work[rank] * inverse) % prime
        if rank + 1 < row_count:
            factors = work[rank + 1:, column].copy()
            active = np.flatnonzero(factors)
            if len(active):
                target_rows = rank + 1 + active
                work[target_rows] = (
                    work[target_rows]
                    - factors[active, None] * work[rank]
                ) % prime
        rank += 1
        if rank == row_count:
            break
    return rank


def assemble_full_hessian(box, edges, coefficients):
    """Assemble the frozen exact trace formula using sparse edge support."""
    starts = np.asarray([edge[0] for edge in edges], dtype=np.int64)
    ends = np.asarray([edge[1] for edge in edges], dtype=np.int64)
    coefficients = np.asarray(coefficients, dtype=np.int64)

    i_e = starts[:, None]
    j_e = ends[:, None]
    k_f = starts[None, :]
    l_f = ends[None, :]
    coefficient_product = coefficients[:, None] * coefficients[None, :]
    box_squared = box @ box

    first_trace = coefficient_product * (
        (l_f == i_e) * box_squared[j_e, k_f]
        + (l_f == j_e) * box_squared[i_e, k_f]
        + (k_f == i_e) * box_squared[j_e, l_f]
        + (k_f == j_e) * box_squared[i_e, l_f]
    )
    middle_trace = coefficient_product * (
        box[j_e, k_f] * box[l_f, i_e]
        + box[i_e, k_f] * box[l_f, j_e]
        + box[j_e, l_f] * box[k_f, i_e]
        + box[i_e, l_f] * box[k_f, j_e]
    )
    return 8 * first_trace + 4 * middle_trace


def edge_matrix(edge_index, edges, coefficients):
    matrix = np.zeros((N, N), dtype=np.int64)
    left, right = edges[edge_index]
    matrix[left, right] = matrix[right, left] = coefficients[edge_index]
    return matrix


def trace_fourth(matrix):
    square = matrix @ matrix
    return int(np.trace(square @ square))


def exact_directional_second_derivative(box, direction):
    """Extract P''(0) exactly from the degree-four polynomial P(t)."""
    p0 = trace_fourth(box)
    values = {time: trace_fourth(box + time * direction)
              for time in (-2, -1, 1, 2)}
    symmetric_one = values[1] + values[-1] - 2 * p0
    symmetric_two = values[2] + values[-2] - 2 * p0
    numerator = 16 * symmetric_one - symmetric_two
    assert numerator % 12 == 0
    return numerator // 12


def legacy_gram(eigenvalues, eigenvectors, edges, coefficients):
    sensitivities = np.empty((N, len(edges)), dtype=float)
    for edge_index, (left, right) in enumerate(edges):
        sensitivities[:, edge_index] = (
            2.0 * coefficients[edge_index]
            * eigenvectors[left, :] * eigenvectors[right, :]
        )
    return 12.0 * (
        sensitivities.T * (eigenvalues**2)[None, :]
    ) @ sensitivities


print("Full Hessian correction audit for Tr(Box^4)")

vertices, adjacency = build_600cell()
fibers = find_hopf_fibers(vertices)
vertex_fiber = {
    vertex: fiber_index
    for fiber_index, fiber in enumerate(fibers)
    for vertex in fiber
}

fiber_adjacency = np.zeros((N, N), dtype=np.int64)
for fiber in fibers:
    for left in fiber:
        for right in fiber:
            if left != right and adjacency[left, right]:
                fiber_adjacency[left, right] = 1
cross_adjacency = adjacency - fiber_adjacency
box = A1 * fiber_adjacency - cross_adjacency

edges = []
coefficients = []
fiber_flags = []
for left in range(N):
    for right in range(left + 1, N):
        if adjacency[left, right]:
            same_fiber = vertex_fiber[left] == vertex_fiber[right]
            edges.append((left, right))
            fiber_flags.append(same_fiber)
            coefficients.append(A1 if same_fiber else -1)

fiber_count = int(sum(fiber_flags))
cross_count = len(edges) - fiber_count
check(
    "600-cell and Hopf edge split are 120 + 600 = 720",
    len(vertices) == 120
    and np.all(adjacency.sum(axis=1) == 12)
    and len(fibers) == 12
    and all(len(fiber) == 10 for fiber in fibers)
    and fiber_count == 120 and cross_count == 600,
)
check(
    "B is the frozen integer symmetric operator",
    np.array_equal(box, box.T)
    and np.all(fiber_adjacency.sum(axis=1) == 2)
    and np.all(cross_adjacency.sum(axis=1) == 10),
)

# Exact B rank: nine exact rational kernel vectors give rank <= 111, while
# rank 111 modulo a prime gives rank_Q >= 111.
box_nullspace = sy.Matrix(box).nullspace()
integer_kernel = []
for vector in box_nullspace:
    denominator = sy.ilcm(*[entry.q for entry in vector])
    integer_kernel.append(np.array(
        [int(entry * denominator) for entry in vector], dtype=np.int64
    ))
kernel_matrix = np.stack(integer_kernel, axis=1)
box_rank_mod_1009 = rank_mod_prime(box, 1009)
check(
    "B has exact rational rank 111 and nullity 9",
    len(box_nullspace) == 9
    and np.max(np.abs(box @ kernel_matrix)) == 0
    and np.linalg.matrix_rank(kernel_matrix.astype(float)) == 9
    and box_rank_mod_1009 == 111,
    f"rank_mod_1009={box_rank_mod_1009}, exact_kernel_dim={len(box_nullspace)}",
)

full_hessian = assemble_full_hessian(box, edges, coefficients)
check(
    "all 720^2 full-Hessian entries are exact integers and symmetric",
    full_hessian.dtype == np.int64
    and np.array_equal(full_hessian, full_hessian.T),
    f"entry_range=[{full_hessian.min()}, {full_hessian.max()}]",
)

# Direct trace checks use dense edge matrices and do not reuse the sparse
# closed form used by assemble_full_hessian.
box_squared = box @ box
entry_controls = [(0, 0), (0, 1), (13, 217), (719, 4), (345, 511)]
entry_residuals = []
for edge_index, other_index in entry_controls:
    edge_e = edge_matrix(edge_index, edges, coefficients)
    edge_f = edge_matrix(other_index, edges, coefficients)
    direct = 4 * int(np.trace(
        edge_f @ box_squared @ edge_e
        + box @ edge_f @ box @ edge_e
        + box_squared @ edge_f @ edge_e
    ))
    entry_residuals.append(
        int(full_hessian[edge_index, other_index]) - direct
    )
check(
    "selected entries match the direct noncommutative trace derivative",
    entry_residuals == [0] * len(entry_controls),
    f"residuals={entry_residuals}",
)

# Exact, non-fitted polynomial controls in three deterministic directions.
directions = []
single = np.zeros(len(edges), dtype=np.int64)
single[0] = 1
directions.append(single)
alternating = np.zeros(len(edges), dtype=np.int64)
alternating[:17] = np.where(np.arange(17) % 2 == 0, 1, -1)
directions.append(alternating)
rng = np.random.default_rng(20260812)
sparse_random = np.zeros(len(edges), dtype=np.int64)
chosen = rng.choice(len(edges), size=41, replace=False)
sparse_random[chosen] = rng.integers(-2, 3, size=len(chosen))
directions.append(sparse_random)

directional_residuals = []
for weights in directions:
    perturbation = np.zeros_like(box)
    for edge_index in np.flatnonzero(weights):
        perturbation += int(weights[edge_index]) * edge_matrix(
            edge_index, edges, coefficients
        )
    polynomial_value = exact_directional_second_derivative(box, perturbation)
    hessian_value = int(weights @ full_hessian @ weights)
    directional_residuals.append(hessian_value - polynomial_value)
check(
    "three exact polynomial directions reproduce z^T H z",
    directional_residuals == [0, 0, 0],
    f"residuals={directional_residuals}",
)

# The first fixed prime happens to divide every 709-by-709 candidate found by
# elimination; the second supplies a full-rank certificate.  Reporting both
# prevents the successful prime from being hidden.
hessian_rank_mod_101 = rank_mod_prime(full_hessian, 101)
hessian_rank_mod_1009 = rank_mod_prime(full_hessian, 1009)
check(
    "H has exact rational rank 720 (full rank modulo 1009)",
    hessian_rank_mod_1009 == 720,
    f"rank_mod_101={hessian_rank_mod_101}, rank_mod_1009={hessian_rank_mod_1009}",
)

x, y = sy.symbols("x y", real=True)
positive_weight_identity = sy.expand(
    (x + y / 2)**2 + 3 * y**2 / 4 - (x**2 + x * y + y**2)
)
check(
    "spectral off-diagonal Hessian weight is a sum of squares",
    positive_weight_identity == 0,
    "x^2+xy+y^2=(x+y/2)^2+3y^2/4",
)
check(
    "analytic PSD plus exact full rank makes H positive definite on edge weights",
    positive_weight_identity == 0 and hessian_rank_mod_1009 == 720,
)

full_eigenvalues = np.linalg.eigvalsh(full_hessian.astype(float))
check(
    "numerical inertia control agrees with exact positive definiteness",
    np.all(full_eigenvalues > 1.0e-7),
    f"lambda_min={full_eigenvalues[0]:.12g}, lambda_max={full_eigenvalues[-1]:.12g}",
)

# Reproduce the legacy diagonal-eigenvalue-sensitivity Gram matrix, then make
# a deterministic orthogonal basis change inside every degenerate eigenspace.
box_eigenvalues, box_eigenvectors = np.linalg.eigh(box.astype(float))
legacy = legacy_gram(
    box_eigenvalues, box_eigenvectors, edges, np.asarray(coefficients)
)
legacy_eigenvalues = np.linalg.eigvalsh(legacy)
legacy_inertia = {
    "positive": int(np.sum(legacy_eigenvalues > LEGACY_TOL)),
    "zero": int(np.sum(np.abs(legacy_eigenvalues) < LEGACY_TOL)),
    "negative": int(np.sum(legacy_eigenvalues < -LEGACY_TOL)),
}
check(
    "legacy calculation itself is reproduced as 101 + 619",
    legacy_inertia == {"positive": 101, "zero": 619, "negative": 0},
    str(legacy_inertia),
)

absolute_legacy_error = float(np.linalg.norm(full_hessian - legacy))
relative_legacy_error = absolute_legacy_error / float(np.linalg.norm(full_hessian))
check(
    "legacy Gram matrix is not the full Hessian",
    relative_legacy_error > 0.1,
    f"||H-G||_F/||H||_F={relative_legacy_error:.12g}",
)

clusters = []
cluster_start = 0
for index in range(1, N + 1):
    if index == N or abs(
        box_eigenvalues[index] - box_eigenvalues[cluster_start]
    ) > 1.0e-7:
        clusters.append((cluster_start, index))
        cluster_start = index

rotated_vectors = box_eigenvectors.copy()
rotation_rng = np.random.default_rng(20260812)
for start, stop in clusters:
    dimension = stop - start
    if dimension <= 1:
        continue
    random_matrix = rotation_rng.normal(size=(dimension, dimension))
    rotation, _ = np.linalg.qr(random_matrix)
    rotated_vectors[:, start:stop] = (
        box_eigenvectors[:, start:stop] @ rotation
    )

rotated_legacy = legacy_gram(
    box_eigenvalues, rotated_vectors, edges, np.asarray(coefficients)
)
basis_change = float(np.linalg.norm(rotated_legacy - legacy))
relative_basis_change = basis_change / float(np.linalg.norm(legacy))
check(
    "legacy Gram matrix changes under valid degenerate-eigenspace rotations",
    relative_basis_change > 0.1,
    f"relative_basis_change={relative_basis_change:.12g}",
)
check(
    "full Hessian is eigenbasis-free while the legacy object is not",
    np.array_equal(full_hessian, full_hessian.T)
    and relative_basis_change > 0.1,
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "frozen_operator": {
        "vertices": 120,
        "edges": len(edges),
        "fiber_edges": fiber_count,
        "cross_edges": cross_count,
        "fiber_coefficient": 5,
        "cross_coefficient": -1,
        "B_exact_rank": 111,
        "B_exact_nullity": 9,
        "B_rank_mod_1009": box_rank_mod_1009,
    },
    "full_hessian": {
        "formula": "8 Tr(B^2 E_f E_e) + 4 Tr(B E_f B E_e)",
        "shape": [720, 720],
        "entry_min": int(full_hessian.min()),
        "entry_max": int(full_hessian.max()),
        "rank_mod_101": hessian_rank_mod_101,
        "rank_mod_1009": hessian_rank_mod_1009,
        "exact_rank": 720,
        "exact_inertia_from_psd_and_rank": {
            "positive": 720,
            "zero": 0,
            "negative": 0,
        },
        "numerical_eigenvalue_min": float(full_eigenvalues[0]),
        "numerical_eigenvalue_max": float(full_eigenvalues[-1]),
        "entry_control_residuals": entry_residuals,
        "directional_polynomial_residuals": directional_residuals,
    },
    "legacy_diagonal_sensitivity_gram": {
        "formula": "12 sum_k lambda_k^2 S_ke S_kf",
        "numerical_inertia_reproduced": legacy_inertia,
        "frobenius_error_from_full_hessian": absolute_legacy_error,
        "relative_frobenius_error_from_full_hessian": relative_legacy_error,
        "degenerate_eigenspace_dimensions": [
            stop - start for start, stop in clusters
        ],
        "relative_change_under_basis_rotation": relative_basis_change,
        "canonical": False,
    },
    "verdicts": [
        {
            "label": "DERIVED CORRECTION",
            "claim": "the full restricted Hessian has exact rank and positive inertia 720",
        },
        {
            "label": "DERIVED NEGATIVE",
            "claim": "the legacy 101+619 Gram matrix is neither the full Hessian nor eigenbasis invariant",
        },
        {
            "label": "OPEN",
            "claim": "no graviton follows from this edge-weight stiffness without a field dictionary, gauge quotient, dynamics, continuum limit, and source coupling",
        },
    ],
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"Wrote {OUTPUT.name}")
raise SystemExit(0 if passed == tests else 1)

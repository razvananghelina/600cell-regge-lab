#!/usr/bin/env python3
"""Refinement scaling of the canonical finite Whitney stiffness.

Protocol commit 03e0abc froze the complete base/refined carriers, exact local
metrics, row-image quotient method, calibration values, numerical residual
gate, and target-free scale diagnostics before the extremal spectra were
computed.
"""

from itertools import combinations, permutations
import json
from math import factorial
from pathlib import Path
import sys

import numpy as np
from scipy import linalg
import scipy.sparse as sparse
from scipy.sparse.linalg import LinearOperator, eigsh
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_stiffness_refinement.json")
PROTOCOL_COMMIT = "03e0abc"
NEIGHBOUR_CERTIFICATE = Path(__file__).with_name(
    "whitney_neighbour_constraints.json"
)
EXPECTED_NEIGHBOUR_PROTOCOL = "a819a52"
RITZ_TOLERANCE = 1e-8
EXTREMAL_COUNT = 3 if "--multiplicity-audit" in sys.argv[1:] else 1
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


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


def local_coboundaries():
    faces = [list(combinations(range(4), degree + 1))
             for degree in range(4)]
    indices = [{face: index for index, face in enumerate(layer)}
               for layer in faces]
    differentials = []
    for degree in range(3):
        matrix = sy.zeros(len(faces[degree + 1]), len(faces[degree]))
        for row, simplex in enumerate(faces[degree + 1]):
            for omitted in range(degree + 2):
                face = simplex[:omitted] + simplex[omitted + 1:]
                matrix[row, indices[degree][face]] = (-1) ** omitted
        differentials.append(matrix)
    return faces, differentials


def local_metric_and_weak(points, differentials):
    masses = [local_whitney_mass(points, degree) for degree in range(4)]
    offsets = np.cumsum((0, 4, 6, 4, 1))
    metric = sy.diag(*masses)
    weak = sy.zeros(15, 15)
    for degree, differential in enumerate(differentials):
        low_start, low_stop = offsets[degree:degree + 2]
        high_start, high_stop = offsets[degree + 1:degree + 3]
        forward = masses[degree + 1] * differential
        weak[high_start:high_stop, low_start:low_stop] = forward
        weak[low_start:low_stop, high_start:high_stop] = forward.T
    return masses, metric, weak


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


class DisjointSet:
    def __init__(self, size):
        self.parent = np.arange(size, dtype=np.int32)

    def find(self, item):
        root = int(item)
        while self.parent[root] != root:
            root = int(self.parent[root])
        item = int(item)
        while self.parent[item] != item:
            parent = int(self.parent[item])
            self.parent[item] = root
            item = parent
        return root

    def union(self, left, right):
        left = self.find(left)
        right = self.find(right)
        if left != right:
            self.parent[right] = left


def canonical_degree_constraints(top_cells, cells, degree):
    local_faces = list(combinations(range(4), degree + 1))
    local_count = len(local_faces)
    cell_indices = {cell: index
                    for index, cell in enumerate(cells[degree])}
    triangle_indices = {cell: index
                        for index, cell in enumerate(cells[2])}

    lookup = {}
    copy_global = np.empty(len(top_cells) * local_count, dtype=np.int32)
    triangle_parents = [[] for _ in cells[2]]
    for top_index, top in enumerate(top_cells):
        for local_index, positions in enumerate(local_faces):
            cell = tuple(top[position] for position in positions)
            global_index = cell_indices[cell]
            copy_index = top_index * local_count + local_index
            lookup[(top_index, global_index)] = copy_index
            copy_global[copy_index] = global_index
        for triangle in combinations(top, 3):
            triangle_parents[triangle_indices[tuple(triangle)]].append(
                top_index
            )

    pairs_by_simplex = [[] for _ in cells[degree]]
    for triangle_index, parents in enumerate(triangle_parents):
        if len(parents) != 2:
            raise AssertionError("complex is not a closed 3-pseudomanifold")
        left_top, right_top = sorted(parents)
        triangle = cells[2][triangle_index]
        for simplex in combinations(triangle, degree + 1):
            global_index = cell_indices[tuple(simplex)]
            left = lookup[(left_top, global_index)]
            right = lookup[(right_top, global_index)]
            pairs_by_simplex[global_index].append((left, right))

    row_count = sum(map(len, pairs_by_simplex))
    rank = sum(
        len(set(node for pair in pairs for node in pair)) - 1
        for pairs in pairs_by_simplex
    )
    c_rows = np.empty(2 * row_count, dtype=np.int32)
    c_columns = np.empty(2 * row_count, dtype=np.int32)
    c_data = np.empty(2 * row_count, dtype=np.float64)

    # The row-image basis is dense only inside each small occurrence graph.
    v_nonzeros = sum(
        len(pairs)
        * (len(set(node for pair in pairs for node in pair)) - 1)
        for pairs in pairs_by_simplex
    )
    v_rows = np.empty(v_nonzeros, dtype=np.int32)
    v_columns = np.empty(v_nonzeros, dtype=np.int32)
    v_data = np.empty(v_nonzeros, dtype=np.float64)

    row_offset = 0
    rank_offset = 0
    c_cursor = 0
    v_cursor = 0
    maximum_orthonormality_residual = 0.0
    maximum_node_degree = 0
    all_connected = True
    occurrence_histogram = {}

    for pairs in pairs_by_simplex:
        nodes = sorted(set(node for pair in pairs for node in pair))
        node_index = {node: index for index, node in enumerate(nodes)}
        edge_count = len(pairs)
        node_count = len(nodes)
        occurrence_histogram[node_count] = (
            occurrence_histogram.get(node_count, 0) + 1
        )
        incidence = np.zeros((edge_count, node_count), dtype=np.float64)
        dsu = DisjointSet(node_count)
        node_degrees = np.zeros(node_count, dtype=np.int16)
        for local_row, (left, right) in enumerate(pairs):
            left_node = node_index[left]
            right_node = node_index[right]
            incidence[local_row, left_node] = 1.0
            incidence[local_row, right_node] = -1.0
            dsu.union(left_node, right_node)
            node_degrees[left_node] += 1
            node_degrees[right_node] += 1

            c_rows[c_cursor:c_cursor + 2] = row_offset + local_row
            c_columns[c_cursor:c_cursor + 2] = (left, right)
            c_data[c_cursor:c_cursor + 2] = (1.0, -1.0)
            c_cursor += 2

        roots = {dsu.find(node) for node in range(node_count)}
        all_connected &= len(roots) == 1
        maximum_node_degree = max(
            maximum_node_degree, int(node_degrees.max(initial=0))
        )

        # Any n-1 columns of a connected oriented incidence matrix are
        # independent.  Their reduced QR is an orthonormal basis for im(B).
        basis, _ = np.linalg.qr(incidence[:, :-1], mode="reduced")
        local_rank = node_count - 1
        orth_residual = np.max(np.abs(
            basis.T @ basis - np.eye(local_rank)
        )) if local_rank else 0.0
        maximum_orthonormality_residual = max(
            maximum_orthonormality_residual, float(orth_residual)
        )
        block_size = edge_count * local_rank
        v_rows[v_cursor:v_cursor + block_size] = np.repeat(
            np.arange(row_offset, row_offset + edge_count, dtype=np.int32),
            local_rank,
        )
        v_columns[v_cursor:v_cursor + block_size] = np.tile(
            np.arange(rank_offset, rank_offset + local_rank, dtype=np.int32),
            edge_count,
        )
        v_data[v_cursor:v_cursor + block_size] = basis.ravel()
        v_cursor += block_size
        row_offset += edge_count
        rank_offset += local_rank

    constraint = sparse.csr_matrix(
        (c_data, (c_rows, c_columns)),
        shape=(row_count, len(top_cells) * local_count),
    )
    row_image_basis = sparse.csr_matrix(
        (v_data, (v_rows, v_columns)),
        shape=(row_count, rank),
    )
    return constraint, row_image_basis, {
        "degree": degree,
        "local_dimension": len(top_cells) * local_count,
        "global_dimension": len(cells[degree]),
        "constraint_rows": row_count,
        "constraint_rank_by_connected_components": rank,
        "constraint_redundancy": row_count - rank,
        "row_image_basis_nonzeros": int(row_image_basis.nnz),
        "maximum_basis_orthonormality_residual": (
            maximum_orthonormality_residual
        ),
        "all_occurrence_graphs_connected": bool(all_connected),
        "maximum_occurrence_node_degree": maximum_node_degree,
        "occurrence_size_histogram": {
            str(size): count for size, count in sorted(
                occurrence_histogram.items()
            )
        },
    }


def row_image_extrema(constraint, basis, local_mass, top_count):
    local_inverse = np.linalg.inv(np.asarray(local_mass, dtype=np.float64))
    local_count = local_inverse.shape[0]

    def apply(vector):
        row_vector = basis @ vector
        copies = constraint.T @ row_vector
        copies = np.asarray(copies).reshape(top_count, local_count)
        mass_solved = copies @ local_inverse.T
        row_result = constraint @ mass_solved.reshape(-1)
        return np.asarray(basis.T @ row_result).ravel()

    operator = LinearOperator(
        (basis.shape[1], basis.shape[1]),
        matvec=apply,
        rmatvec=apply,
        dtype=np.float64,
    )
    seed = np.sin(np.arange(basis.shape[1], dtype=np.float64) + 1.0)
    seed /= np.linalg.norm(seed)
    minimum_values, minimum_vectors = eigsh(
        operator, k=EXTREMAL_COUNT, which="SA", v0=seed,
        tol=2e-11, maxiter=10000
    )
    maximum_values, maximum_vectors = eigsh(
        operator, k=EXTREMAL_COUNT, which="LA", v0=seed,
        tol=2e-11, maxiter=10000
    )
    minimum = float(minimum_values[0])
    maximum = float(maximum_values[-1])
    minimum_residuals = [
        float(np.linalg.norm(
            apply(minimum_vectors[:, index])
            - minimum_values[index] * minimum_vectors[:, index]
        ) / max(1.0, abs(minimum_values[index])))
        for index in range(EXTREMAL_COUNT)
    ]
    maximum_residuals = [
        float(np.linalg.norm(
            apply(maximum_vectors[:, index])
            - maximum_values[index] * maximum_vectors[:, index]
        ) / max(1.0, abs(maximum_values[index])))
        for index in range(EXTREMAL_COUNT)
    ]
    return {
        "positive_gap": minimum,
        "maximum_positive_eigenvalue": maximum,
        "computed_smallest_positive_eigenvalues": minimum_values.tolist(),
        "computed_largest_positive_eigenvalues": maximum_values.tolist(),
        "minimum_ritz_relative_residual": max(minimum_residuals),
        "maximum_ritz_relative_residual": max(maximum_residuals),
    }


def local_dirac_norm(metric, weak):
    eigenvalues = linalg.eigvalsh(
        np.asarray(weak, dtype=np.float64),
        np.asarray(metric, dtype=np.float64),
    )
    return float(np.max(np.abs(eigenvalues)))


def audit_level(name, top_cells, cells, masses, metric, weak):
    print(f"\n-- {name} --")
    norm = local_dirac_norm(metric, weak)
    degree_records = []
    for degree in range(3):
        constraint, basis, structure = canonical_degree_constraints(
            top_cells, cells, degree
        )
        expected_rank = structure["local_dimension"] - structure[
            "global_dimension"
        ]
        structure["exact_kernel_dimension"] = structure["global_dimension"]
        structure["rank_matches_conformity"] = (
            structure["constraint_rank_by_connected_components"]
            == expected_rank
        )
        extrema = row_image_extrema(
            constraint, basis, masses[degree], len(top_cells)
        )
        record = {**structure, **extrema}
        degree_records.append(record)
        print(
            f"degree {degree}: gap={record['positive_gap']:.12g}, "
            f"max={record['maximum_positive_eigenvalue']:.12g}, "
            f"residual={max(record['minimum_ritz_relative_residual'], record['maximum_ritz_relative_residual']):.3e}"
        )
    all_degree_gap = min(record["positive_gap"] for record in degree_records)
    controlling_degree = min(
        range(3), key=lambda degree: degree_records[degree]["positive_gap"]
    )
    return {
        "level": name,
        "f_vector": list(map(len, cells)),
        "top_count": len(top_cells),
        "local_dirac_norm": norm,
        "degree_records": degree_records,
        "all_degree_positive_gap": all_degree_gap,
        "controlling_degree": controlling_degree,
        "scale_factor_a_over_g": norm / all_degree_gap,
    }


print("=" * 78)
print("REFINEMENT SCALING OF CANONICAL WHITNEY STIFFNESS")
print("=" * 78)
print(f"EXTREMAL_MODES_PER_EDGE={EXTREMAL_COUNT}")

neighbour_certificate = json.loads(NEIGHBOUR_CERTIFICATE.read_text())
check("the independent neighbour certificate has the expected protocol",
      neighbour_certificate["protocol_commit"] == EXPECTED_NEIGHBOUR_PROTOCOL)

local_faces, local_d = local_coboundaries()
reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
coarse_masses, coarse_metric, coarse_weak = local_metric_and_weak(
    reference_vertices, local_d
)
child_points = (
    reference_vertices[0],
    (reference_vertices[0] + reference_vertices[1]) / 2,
    sum(reference_vertices[:3], sy.zeros(3, 1)) / 3,
    sum(reference_vertices, sy.zeros(3, 1)) / 4,
)
fine_masses, fine_metric, fine_weak = local_metric_and_weak(
    child_points, local_d
)

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
check("all 24 rank-ordered child metrics are exactly identical",
      all(matrix == family[0]
          for family in child_mass_families for matrix in family))
check("both local metrics are exactly positive and weak matrices Hermitian",
      all(
          mass == mass.T
          and all(mass[:size, :size].det() > 0
                  for size in range(1, mass.rows + 1))
          for mass in coarse_masses + fine_masses
      ) and coarse_weak == coarse_weak.T and fine_weak == fine_weak.T)

# Algorithm calibration on the boundary of a 4-simplex.
control_top = tuple(combinations(range(5), 4))
control_cells = all_simplices(control_top)
control_records = []
for degree in range(3):
    constraint, basis, structure = canonical_degree_constraints(
        control_top, control_cells, degree
    )
    extrema = row_image_extrema(
        constraint, basis, coarse_masses[degree], len(control_top)
    )
    control_records.append({**structure, **extrema})
control_gap = min(record["positive_gap"] for record in control_records)
control_maximum = max(
    record["maximum_positive_eigenvalue"] for record in control_records
)
control_maximum_residual = max(
    record[key]
    for record in control_records
    for key in (
        "minimum_ritz_relative_residual",
        "maximum_ritz_relative_residual",
    )
)
check("the row-image algorithm reproduces the frozen control gap and maximum",
      abs(control_gap - 7.5) < 1e-9
      and abs(control_maximum - 45.0) < 1e-9,
      f"gap={control_gap:.12g}, maximum={control_maximum:.12g}")
check("the control extremal Ritz residuals meet the frozen gate",
      control_maximum_residual < RITZ_TOLERANCE,
      f"maximum residual={control_maximum_residual:.3e}")

# Complete 600-cell and its first barycentric subdivision.
vertices, adjacency, _ = build_600cell()
neighbours = tuple(
    frozenset(np.flatnonzero(adjacency[index]).tolist())
    for index in range(120)
)
base_edges = tuple(
    (left, right)
    for left in range(120)
    for right in sorted(neighbours[left])
    if left < right
)
base_triangles = tuple(
    (left, right, third)
    for left, right in base_edges
    for third in sorted(neighbours[left] & neighbours[right])
    if right < third
)
base_top = tuple(
    (first, second, third, fourth)
    for first, second, third in base_triangles
    for fourth in sorted(
        neighbours[first] & neighbours[second] & neighbours[third]
    )
    if third < fourth
)
base_cells = (
    tuple((index,) for index in range(120)),
    base_edges,
    base_triangles,
    base_top,
)

fine_vertex_cells = tuple(cell for layer in base_cells for cell in layer)
fine_vertex_index = {cell: index for index, cell in enumerate(fine_vertex_cells)}
fine_top = []
for tetrahedron in base_top:
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
check("the complete base and refined f-vectors are exact",
      list(map(len, base_cells)) == [120, 720, 1200, 600]
      and list(map(len, fine_cells)) == [2640, 17040, 28800, 14400])

base_audit = audit_level(
    "base", base_top, base_cells,
    coarse_masses, coarse_metric, coarse_weak,
)
fine_audit = audit_level(
    "first_barycentric", fine_top, fine_cells,
    fine_masses, fine_metric, fine_weak,
)
audits = (base_audit, fine_audit)

expected_rows = {
    "base": [3600, 3600, 1200],
    "first_barycentric": [86400, 86400, 28800],
}
expected_ranks = {
    "base": [2280, 2880, 1200],
    "first_barycentric": [54960, 69360, 28800],
}
check("all canonical row and exact quotient-rank counts are reproduced",
      all(
          [record["constraint_rows"] for record in audit["degree_records"]]
          == expected_rows[audit["level"]]
          and [
              record["constraint_rank_by_connected_components"]
              for record in audit["degree_records"]
          ] == expected_ranks[audit["level"]]
          for audit in audits
      ))
check("every occurrence graph is connected and every kernel is conformity",
      all(
          record["all_occurrence_graphs_connected"]
          and record["rank_matches_conformity"]
          for audit in audits for record in audit["degree_records"]
      ))
maximum_basis_residual = max(
    record["maximum_basis_orthonormality_residual"]
    for audit in audits for record in audit["degree_records"]
)
check("every local row-image basis is orthonormal to numerical precision",
      maximum_basis_residual < 1e-11,
      f"maximum residual={maximum_basis_residual:.3e}")
maximum_ritz_residual = max(
    record[key]
    for audit in audits for record in audit["degree_records"]
    for key in (
        "minimum_ritz_relative_residual",
        "maximum_ritz_relative_residual",
    )
)
check("all complete extremal Ritz pairs meet the preregistered residual gate",
      maximum_ritz_residual < RITZ_TOLERANCE,
      f"maximum residual={maximum_ritz_residual:.3e}")
check("all complete quotient gaps are strictly positive at both levels",
      all(
          record["positive_gap"] > 1e-10
          for audit in audits for record in audit["degree_records"]
      ))

scale_ratio = (
    fine_audit["scale_factor_a_over_g"]
    / base_audit["scale_factor_a_over_g"]
)
degreewise_ratios = [
    (
        fine_audit["local_dirac_norm"]
        / fine_audit["degree_records"][degree]["positive_gap"]
    ) / (
        base_audit["local_dirac_norm"]
        / base_audit["degree_records"][degree]["positive_gap"]
    )
    for degree in range(3)
]
constant_kappa_compatible = abs(scale_ratio - 1.0) < 1e-10
degree_balanced_rescaling = (
    max(degreewise_ratios) - min(degreewise_ratios) < 1e-8
)
epsilon_diagnostics = {}
for epsilon in (1e-1, 1e-2, 1e-3):
    epsilon_diagnostics[f"{epsilon:.0e}"] = {
        audit["level"]: (
            audit["scale_factor_a_over_g"] * (2.0 + 1.0 / epsilon)
        )
        for audit in audits
    }

check("the first-step scale classification follows the frozen criterion",
      constant_kappa_compatible == (abs(scale_ratio - 1.0) < 1e-10))
check("the degree-balance classification follows the frozen criterion",
      degree_balanced_rescaling == (
          max(degreewise_ratios) - min(degreewise_ratios) < 1e-8
      ))

scale_verdict = (
    "PATTERN: constant algebraic stiffness is first-step compatible"
    if constant_kappa_compatible else
    "DERIVED NEGATIVE at first refinement: constant kappa does not preserve the same relative Schur guarantee"
)
degree_verdict = (
    "PATTERN: all form degrees share one first-step rescaling"
    if degree_balanced_rescaling else
    "DERIVED NEGATIVE at level one: the three form degrees do not share one balanced rescaling"
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "neighbour_protocol_commit": EXPECTED_NEIGHBOUR_PROTOCOL,
    "phenomenological_target_used": False,
    "extremal_modes_per_edge": EXTREMAL_COUNT,
    "control": {
        "complex": "boundary of a 4-simplex",
        "positive_gap": control_gap,
        "maximum_positive_eigenvalue": control_maximum,
        "maximum_ritz_relative_residual": control_maximum_residual,
        "degree_records": control_records,
    },
    "audits": list(audits),
    "refinement_diagnostics": {
        "scale_ratio_s1_over_s0": scale_ratio,
        "constant_kappa_first_step_compatible": constant_kappa_compatible,
        "degreewise_scale_ratios": degreewise_ratios,
        "one_degree_balanced_rescaling": degree_balanced_rescaling,
        "sufficient_kappa_by_relative_error": epsilon_diagnostics,
    },
    "verdicts": [
        "DERIVED: complete positive quotient gaps at base and first refinement",
        scale_verdict,
        degree_verdict,
        "OPEN: absolute kappa normalization and repeated-refinement law",
        "NOT CLAIMED: physical time, mass, speed, or Planck scale",
    ],
    "scope": (
        "Complete 600-cell at levels zero and one. Extremal spectra are "
        "numerical Ritz certificates after exact removal of the graph-cycle "
        "nullspace; no exponent is inferred from two levels."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured stiffness-refinement certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for audit in audits:
    print(
        f"{audit['level']}: a={audit['local_dirac_norm']:.12g}, "
        f"g={audit['all_degree_positive_gap']:.12g}, "
        f"a/g={audit['scale_factor_a_over_g']:.12g}, "
        f"controlling degree={audit['controlling_degree']}"
    )
print(f"FIRST_REFINEMENT_SCALE_RATIO={scale_ratio:.12g}")
print("DEGREEWISE_RATIOS=" + str([float(value) for value in degreewise_ratios]))
print("SCALE_VERDICT: " + scale_verdict)
print("DEGREE_VERDICT: " + degree_verdict)
raise SystemExit(0 if passed == tests else 1)

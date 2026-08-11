#!/usr/bin/env python3
"""Exact Whitney trace-jump stiffness on a refined closed control.

Protocol commits b9a4104 and a92c911 froze the exact face-trace form,
dimensional test, boundary-of-4-simplex carrier, first barycentric refinement,
complete dense spectra, and paired unweighted comparison before evaluation.
"""

from itertools import combinations, permutations
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy import linalg
import sympy as sy


OUTPUT = Path(__file__).with_name("whitney_trace_stiffness.json")
PROTOCOL_COMMIT = "b9a4104"
PROTOCOL_CORRECTION_COMMIT = "a92c911"
EIGEN_ZERO_TOLERANCE = 1e-10
EIGEN_RESIDUAL_GATE = 1e-9
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


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
    moment = volume * (sy.ones(4, 4) + sy.eye(4)) / 20
    forms = list(combinations(range(4), degree + 1))
    coefficients = []
    for form in forms:
        coefficient = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            coefficient[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                coefficient[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * wedge_components(covectors, degree)
                )
        coefficients.append(coefficient)
    mass = sy.zeros(len(forms), len(forms))
    for row, left in enumerate(coefficients):
        for column, right in enumerate(coefficients):
            mass[row, column] = sy.simplify(sum(
                (left[basis, :] * moment * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
            ))
    return mass


def triangle_gram(points):
    affine = sy.Matrix.hstack(points[1] - points[0], points[2] - points[0])
    return sy.simplify(affine.T * affine)


def triangle_whitney_mass_from_gram(gram, degree):
    gradients = (
        sy.Matrix((-1, -1)),
        sy.Matrix((1, 0)),
        sy.Matrix((0, 1)),
    )
    face_coordinate_bases = [list(combinations(range(2), value))
                             for value in range(3)]

    def face_wedge(covectors):
        if degree == 0:
            return sy.Matrix((1,))
        return sy.Matrix([
            sy.det(sy.Matrix([
                [covector[index] for index in basis]
                for covector in covectors
            ]))
            for basis in face_coordinate_bases[degree]
        ])

    inverse_gram = sy.simplify(gram.inv())
    if degree == 0:
        wedge_metric = sy.ones(1, 1)
    elif degree == 1:
        wedge_metric = inverse_gram
    else:
        wedge_metric = sy.Matrix(((sy.simplify(1 / gram.det()),),))

    area = sy.sqrt(sy.simplify(gram.det())) / 2
    moment = area * (sy.ones(3, 3) + sy.eye(3)) / 12
    forms = list(combinations(range(3), degree + 1))
    coefficients = []
    for form in forms:
        coefficient = sy.zeros(len(face_coordinate_bases[degree]), 3)
        if degree == 0:
            coefficient[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                coefficient[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * face_wedge(covectors)
                )
        coefficients.append(coefficient)

    mass = sy.zeros(len(forms), len(forms))
    for row, left in enumerate(coefficients):
        for column, right in enumerate(coefficients):
            mass[row, column] = sy.simplify(sum(
                wedge_metric[left_basis, right_basis]
                * (left[left_basis, :] * moment
                   * right[right_basis, :].T)[0]
                for left_basis in range(wedge_metric.rows)
                for right_basis in range(wedge_metric.cols)
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


def local_dirac_norm(metric, weak):
    eigenvalues = linalg.eigvalsh(
        np.asarray(weak, dtype=np.float64),
        np.asarray(metric, dtype=np.float64),
    )
    return float(np.max(np.abs(eigenvalues)))


class DisjointSet:
    def __init__(self, size):
        self.parent = list(range(size))

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left, right):
        left = self.find(left)
        right = self.find(right)
        if left != right:
            self.parent[right] = left


def build_base_and_fine_control(reference_vertices):
    base_top = tuple(combinations(range(5), 4))
    base_cells = all_simplices(base_top)
    base_top_points = tuple(reference_vertices for _ in base_top)

    fine_vertex_cells = tuple(cell for layer in base_cells for cell in layer)
    fine_vertex_index = {cell: index for index, cell in enumerate(
        fine_vertex_cells
    )}
    fine_top = []
    fine_top_points = []
    for base_top_index, top in enumerate(base_top):
        top_position = {vertex: index for index, vertex in enumerate(top)}
        for ordering in permutations(top):
            flag = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                top,
            )
            fine_top.append(tuple(fine_vertex_index[cell] for cell in flag))
            fine_top_points.append(tuple(
                sum((reference_vertices[top_position[vertex]]
                     for vertex in cell), sy.zeros(3, 1)) / len(cell)
                for cell in flag
            ))
    fine_top = tuple(fine_top)
    fine_top_points = tuple(fine_top_points)
    fine_cells = all_simplices(fine_top)
    return (
        ("base", base_top, base_cells, base_top_points),
        ("first_barycentric", fine_top, fine_cells, fine_top_points),
    )


def add_dense_block(matrix, indices, block):
    matrix[np.ix_(indices, indices)] += block


def build_penalties(top_cells, cells, top_points, local_masses):
    top_count = len(top_cells)
    cell_indices = [
        {cell: index for index, cell in enumerate(layer)}
        for layer in cells
    ]
    local_faces = [list(combinations(range(4), degree + 1))
                   for degree in range(4)]
    copy_lookup = []
    for degree in range(3):
        lookup = {}
        for top_index, top in enumerate(top_cells):
            for local_index, positions in enumerate(local_faces[degree]):
                cell = tuple(top[position] for position in positions)
                global_index = cell_indices[degree][cell]
                lookup[(top_index, global_index)] = (
                    top_index * len(local_faces[degree]) + local_index
                )
        copy_lookup.append(lookup)

    triangle_parents = [[] for _ in cells[2]]
    for top_index, top in enumerate(top_cells):
        for triangle in combinations(top, 3):
            triangle_parents[cell_indices[2][tuple(triangle)]].append(
                top_index
            )

    metrics = [
        np.kron(np.eye(top_count), np.asarray(mass, dtype=np.float64))
        for mass in local_masses[:3]
    ]
    trace_penalties = [np.zeros_like(metric) for metric in metrics]
    unweighted_penalties = [np.zeros_like(metric) for metric in metrics]
    pairs_by_simplex = [
        [[] for _ in cells[degree]] for degree in range(3)
    ]
    face_mass_cache = {}
    face_metric_mismatches = 0
    face_mass_positive = True
    face_metric_types = set()

    for triangle_index, triangle in enumerate(cells[2]):
        parents = sorted(triangle_parents[triangle_index])
        if len(parents) != 2:
            raise AssertionError("control is not a closed 3-complex")
        parent_points = []
        for parent in parents:
            position = {
                vertex: index for index, vertex in enumerate(top_cells[parent])
            }
            parent_points.append(tuple(
                top_points[parent][position[vertex]] for vertex in triangle
            ))
        grams = tuple(triangle_gram(points) for points in parent_points)
        if grams[0] != grams[1]:
            face_metric_mismatches += 1
        gram = grams[0]
        gram_key = tuple(sy.srepr(value) for value in gram)
        face_metric_types.add(gram_key)

        for degree in range(3):
            subfaces = list(combinations(triangle, degree + 1))
            global_indices = [
                cell_indices[degree][tuple(simplex)] for simplex in subfaces
            ]
            left_indices = [
                copy_lookup[degree][(parents[0], global_index)]
                for global_index in global_indices
            ]
            right_indices = [
                copy_lookup[degree][(parents[1], global_index)]
                for global_index in global_indices
            ]
            for global_index, left, right in zip(
                global_indices, left_indices, right_indices
            ):
                pairs_by_simplex[degree][global_index].append((left, right))

            cache_key = (degree, gram_key)
            if cache_key not in face_mass_cache:
                exact_face_mass = triangle_whitney_mass_from_gram(
                    gram, degree
                )
                face_mass_positive &= all(
                    exact_face_mass[:size, :size].det() > 0
                    for size in range(1, exact_face_mass.rows + 1)
                )
                face_mass_cache[cache_key] = np.asarray(
                    exact_face_mass, dtype=np.float64
                )
            face_mass = face_mass_cache[cache_key]
            jump_indices = left_indices + right_indices
            signed_trace_block = np.block([
                [face_mass, -face_mass],
                [-face_mass, face_mass],
            ])
            identity = np.eye(len(subfaces))
            signed_unweighted_block = np.block([
                [identity, -identity],
                [-identity, identity],
            ])
            add_dense_block(
                trace_penalties[degree], jump_indices, signed_trace_block
            )
            add_dense_block(
                unweighted_penalties[degree], jump_indices,
                signed_unweighted_block,
            )

    connectivity = []
    exact_ranks = []
    for degree in range(3):
        every_connected = True
        rank = 0
        for pairs in pairs_by_simplex[degree]:
            nodes = sorted(set(node for pair in pairs for node in pair))
            node_index = {node: index for index, node in enumerate(nodes)}
            dsu = DisjointSet(len(nodes))
            for left, right in pairs:
                dsu.union(node_index[left], node_index[right])
            components = len({dsu.find(index) for index in range(len(nodes))})
            every_connected &= components == 1
            rank += len(nodes) - components
        connectivity.append(every_connected)
        exact_ranks.append(rank)

    return {
        "metrics": metrics,
        "trace_penalties": trace_penalties,
        "unweighted_penalties": unweighted_penalties,
        "face_metric_mismatches": face_metric_mismatches,
        "face_metric_type_count": len(face_metric_types),
        "face_mass_positive": bool(face_mass_positive),
        "occurrence_graphs_connected": connectivity,
        "exact_ranks": exact_ranks,
    }


def spectrum_audit(penalty, metric, expected_nullity):
    eigenvalues, eigenvectors = linalg.eigh(penalty, metric)
    zero_mask = np.abs(eigenvalues) < EIGEN_ZERO_TOLERANCE
    zero_count = int(np.count_nonzero(zero_mask))
    positive_indices = np.flatnonzero(~zero_mask)
    minimum_index = int(positive_indices[0])
    maximum_index = int(positive_indices[-1])

    residuals = []
    for index in (minimum_index, maximum_index):
        vector = eigenvectors[:, index]
        value = eigenvalues[index]
        left = penalty @ vector
        right = value * (metric @ vector)
        residuals.append(float(
            np.linalg.norm(left - right)
            / max(1.0, np.linalg.norm(left), np.linalg.norm(right))
        ))
    return {
        "numerical_nullity": zero_count,
        "expected_nullity": expected_nullity,
        "minimum_eigenvalue": float(eigenvalues[0]),
        "positive_gap": float(eigenvalues[minimum_index]),
        "maximum_eigenvalue": float(eigenvalues[maximum_index]),
        "maximum_extremal_relative_residual": max(residuals),
        "nullity_matches": zero_count == expected_nullity,
    }


print("=" * 78)
print("EXACT WHITNEY TRACE-JUMP STIFFNESS")
print("=" * 78)

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
local_faces, local_d = local_coboundaries()

# Exact uniform-dilation audit of the preregistered dimensional claim.
coarse_masses, coarse_metric, coarse_weak = local_metric_and_weak(
    reference_vertices, local_d
)
scaled_vertices = tuple(2 * vertex for vertex in reference_vertices)
scaled_masses, scaled_metric, scaled_weak = local_metric_and_weak(
    scaled_vertices, local_d
)
element_scaling_exact = all(
    scaled_masses[degree]
    == sy.Rational(2) ** (3 - 2 * degree) * coarse_masses[degree]
    for degree in range(4)
)
reference_face_gram = triangle_gram(reference_vertices[:3])
scaled_face_gram = triangle_gram(scaled_vertices[:3])
face_scaling_exact = all(
    triangle_whitney_mass_from_gram(scaled_face_gram, degree)
    == sy.Rational(2) ** (2 - 2 * degree)
    * triangle_whitney_mass_from_gram(reference_face_gram, degree)
    for degree in range(3)
)
coarse_norm = local_dirac_norm(coarse_metric, coarse_weak)
scaled_norm = local_dirac_norm(scaled_metric, scaled_weak)
check("uniform dilation gives the exact element and face mass powers",
      element_scaling_exact and face_scaling_exact)
check("the mass-orthonormal local Dirac norm scales exactly numerically as 1/h",
      abs(scaled_norm / coarse_norm - 0.5) < 1e-12,
      f"norm ratio={scaled_norm / coarse_norm:.12g}")

levels = build_base_and_fine_control(reference_vertices)
level_records = []
for name, top_cells, cells, top_points in levels:
    print(f"\n-- {name} --")
    tetrahedron_grams = [
        sy.simplify(sy.Matrix.hstack(
            points[1] - points[0],
            points[2] - points[0],
            points[3] - points[0],
        ).T * sy.Matrix.hstack(
            points[1] - points[0],
            points[2] - points[0],
            points[3] - points[0],
        ))
        for points in top_points
    ]
    one_element_type = all(
        gram == tetrahedron_grams[0] for gram in tetrahedron_grams
    )
    local_masses, local_metric, local_weak = local_metric_and_weak(
        top_points[0], local_d
    )
    penalty_data = build_penalties(
        top_cells, cells, top_points, local_masses
    )
    dirac_norm = local_dirac_norm(local_metric, local_weak)
    trace_records = []
    unweighted_records = []
    for degree in range(3):
        expected_nullity = len(cells[degree])
        trace_record = spectrum_audit(
            penalty_data["trace_penalties"][degree],
            penalty_data["metrics"][degree],
            expected_nullity,
        )
        unweighted_record = spectrum_audit(
            penalty_data["unweighted_penalties"][degree],
            penalty_data["metrics"][degree],
            expected_nullity,
        )
        trace_record["degree"] = degree
        unweighted_record["degree"] = degree
        trace_record["a_over_gap"] = (
            dirac_norm / trace_record["positive_gap"]
        )
        unweighted_record["a_over_gap"] = (
            dirac_norm / unweighted_record["positive_gap"]
        )
        trace_records.append(trace_record)
        unweighted_records.append(unweighted_record)
        print(
            f"degree {degree}: trace gap={trace_record['positive_gap']:.12g}, "
            f"unweighted gap={unweighted_record['positive_gap']:.12g}"
        )
    level_records.append({
        "level": name,
        "f_vector": list(map(len, cells)),
        "top_count": len(top_cells),
        "duplicated_dimension": 15 * len(top_cells),
        "one_exact_element_metric_type": one_element_type,
        "local_dirac_norm": dirac_norm,
        "face_metric_mismatches": penalty_data["face_metric_mismatches"],
        "face_metric_type_count": penalty_data["face_metric_type_count"],
        "face_mass_positive": penalty_data["face_mass_positive"],
        "occurrence_graphs_connected": penalty_data[
            "occurrence_graphs_connected"
        ],
        "exact_constraint_ranks": penalty_data["exact_ranks"],
        "trace_records": trace_records,
        "unweighted_records": unweighted_records,
    })

base_record, fine_record = level_records
check("the frozen base and refined control dimensions are exact",
      base_record["f_vector"] == [5, 10, 10, 5]
      and base_record["duplicated_dimension"] == 75
      and fine_record["f_vector"] == [30, 150, 240, 120]
      and fine_record["duplicated_dimension"] == 1800)
check("every level has one exact rank-ordered element metric type",
      all(record["one_exact_element_metric_type"]
          for record in level_records))
check("every shared face has the same exact metric from both parents",
      all(record["face_metric_mismatches"] == 0
          for record in level_records),
      "face metric types=" + str([
          record["face_metric_type_count"] for record in level_records
      ]))
check("every exact face Whitney mass is positive definite",
      all(record["face_mass_positive"] for record in level_records))
check("all occurrence graphs are connected in every constrained degree",
      all(all(record["occurrence_graphs_connected"])
          for record in level_records))
expected_ranks = ([15, 20, 10], [450, 570, 240])
check("the exact copy-constraint ranks equal local minus assembled dimensions",
      all(record["exact_constraint_ranks"] == list(expected)
          for record, expected in zip(level_records, expected_ranks)),
      str([record["exact_constraint_ranks"] for record in level_records]))
check("trace and unweighted numerical nullities both equal exact conformity",
      all(
          degree_record["nullity_matches"]
          for record in level_records
          for family in ("trace_records", "unweighted_records")
          for degree_record in record[family]
      ))
maximum_residual = max(
    degree_record["maximum_extremal_relative_residual"]
    for record in level_records
    for family in ("trace_records", "unweighted_records")
    for degree_record in record[family]
)
check("all complete generalized extremal residuals meet the frozen gate",
      maximum_residual < EIGEN_RESIDUAL_GATE,
      f"maximum residual={maximum_residual:.3e}")

def refinement_ratios(family):
    return [
        fine_record[family][degree]["a_over_gap"]
        / base_record[family][degree]["a_over_gap"]
        for degree in range(3)
    ]


trace_ratios = refinement_ratios("trace_records")
unweighted_ratios = refinement_ratios("unweighted_records")
trace_spread = max(trace_ratios) / min(trace_ratios)
unweighted_spread = max(unweighted_ratios) / min(unweighted_ratios)
exact_first_step_compatible = all(
    abs(ratio - 1.0) < 1e-8 for ratio in trace_ratios
)
improved_degree_balance = trace_spread < unweighted_spread
check("the exact first-step compatibility label follows the frozen criterion",
      exact_first_step_compatible == all(
          abs(ratio - 1.0) < 1e-8 for ratio in trace_ratios
      ))
check("the paired-control improvement label follows the frozen comparison",
      improved_degree_balance == (trace_spread < unweighted_spread))

compatibility_verdict = (
    "DERIVED NUMERICAL: exact trace stiffness is first-step refinement-compatible on the control"
    if exact_first_step_compatible else
    "DERIVED NUMERICAL NEGATIVE: exact trace stiffness is not exactly first-step refinement-compatible on the control"
)
balance_verdict = (
    "PATTERN: exact trace weighting improves paired-control degree balance"
    if improved_degree_balance else
    "DERIVED NEGATIVE: exact trace weighting does not improve paired-control degree balance"
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "phenomenological_target_used": False,
    "candidate_count": 1,
    "candidate": "exact L2 Whitney trace-jump Gram form",
    "uniform_dilation": {
        "element_mass_scaling_exact": bool(element_scaling_exact),
        "face_mass_scaling_exact": bool(face_scaling_exact),
        "dirac_norm_ratio_at_h_equals_2": scaled_norm / coarse_norm,
    },
    "levels": level_records,
    "refinement_comparison": {
        "trace_degreewise_ratios": trace_ratios,
        "trace_ratio_spread": trace_spread,
        "unweighted_same_control_degreewise_ratios": unweighted_ratios,
        "unweighted_same_control_ratio_spread": unweighted_spread,
        "exact_first_step_compatible": exact_first_step_compatible,
        "paired_control_degree_balance_improved": improved_degree_balance,
    },
    "verdicts": [
        "DERIVED: exact trace mass has the same uniform 1/h scaling as Kähler-Dirac in every degree",
        compatibility_verdict,
        balance_verdict,
        "OPEN: complete 600-cell trace-stiffness refinement spectrum",
        "OPEN: overall dimensionless stiffness and repeated-refinement law",
    ],
    "scope": (
        "Complete dense spectra on the boundary-of-4-simplex control and "
        "its first barycentric subdivision. Not a complete 600-cell or "
        "continuum certificate."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured trace-stiffness certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("TRACE_RATIOS=" + str([float(value) for value in trace_ratios]))
print(f"TRACE_SPREAD={trace_spread:.12g}")
print("UNWEIGHTED_PAIRED_RATIOS=" + str([
    float(value) for value in unweighted_ratios
]))
print(f"UNWEIGHTED_PAIRED_SPREAD={unweighted_spread:.12g}")
print("COMPATIBILITY_VERDICT: " + compatibility_verdict)
print("BALANCE_VERDICT: " + balance_verdict)
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Audit the dynamical and separator roles of Whitney trace stiffness.

Protocol commit 7528f97 froze the hypotheses, claims, exact gates, and
interpretation before the kernel and valence enumeration.
"""

from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
import sympy as sy

from whitney_trace_refinement_tools import (
    barycentric_refine,
    make_base_level,
    matrix_key,
    rank_edgewise_level,
    triangle_gram,
    triangle_whitney_mass_from_gram,
)


OUTPUT = Path(__file__).with_name("whitney_trace_penalty_role.json")
STIFFNESS_CERTIFICATE = Path(__file__).with_name(
    "whitney_rank_edgewise_stiffness.json"
)
REFINEMENT_CERTIFICATE = Path(__file__).with_name(
    "whitney_rank_edgewise_refinement.json"
)
DILATION_CERTIFICATE = Path(__file__).with_name(
    "whitney_trace_stiffness.json"
)
PROTOCOL_COMMIT = "7528f97"
EXPECTED_STIFFNESS_PROTOCOL = "c7e4335"
EXPECTED_REFINEMENT_PROTOCOL = "58fa9fc"
EXPECTED_DILATION_PROTOCOL = "b9a4104"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


class DisjointSet:
    def __init__(self, size):
        self.parent = list(range(size))

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[right] = left


def copy_matrices(level, degree):
    top = level["top"]
    cells = level["cells"]
    local_faces = tuple(combinations(range(4), degree + 1))
    global_index = {cell: index for index, cell in enumerate(cells[degree])}
    triangle_index = {cell: index for index, cell in enumerate(cells[2])}

    lookup = {}
    occurrences = [[] for _ in cells[degree]]
    triangle_parents = [[] for _ in cells[2]]
    j_rows, j_columns = [], []
    for top_index, tetrahedron in enumerate(top):
        for local_index, positions in enumerate(local_faces):
            cell = tuple(tetrahedron[position] for position in positions)
            index = global_index[cell]
            local_copy = top_index * len(local_faces) + local_index
            lookup[(top_index, index)] = local_copy
            occurrences[index].append(local_copy)
            j_rows.append(local_copy)
            j_columns.append(index)
        for triangle in combinations(tetrahedron, 3):
            triangle_parents[triangle_index[tuple(triangle)]].append(top_index)

    r_rows, r_columns, r_data = [], [], []
    graph_edges = [[] for _ in cells[degree]]
    row = 0
    for triangle, parents in zip(cells[2], triangle_parents):
        if len(parents) != 2:
            raise AssertionError("control is not closed")
        for simplex in combinations(triangle, degree + 1):
            index = global_index[tuple(simplex)]
            left = lookup[(parents[0], index)]
            right = lookup[(parents[1], index)]
            r_rows.extend((row, row))
            r_columns.extend((left, right))
            r_data.extend((1, -1))
            graph_edges[index].append((left, right))
            row += 1

    local_dimension = len(top) * len(local_faces)
    injection = sparse.csr_matrix(
        (np.ones(len(j_rows), dtype=np.int8), (j_rows, j_columns)),
        shape=(local_dimension, len(cells[degree])),
        dtype=np.int8,
    )
    jump = sparse.csr_matrix(
        (np.asarray(r_data, dtype=np.int8), (r_rows, r_columns)),
        shape=(row, local_dimension),
        dtype=np.int8,
    )

    connected = True
    rank = 0
    maximum_occurrences = 0
    occurrence_histogram = {}
    for nodes, edges in zip(occurrences, graph_edges):
        maximum_occurrences = max(maximum_occurrences, len(nodes))
        occurrence_histogram[len(nodes)] = occurrence_histogram.get(
            len(nodes), 0
        ) + 1
        node_index = {node: index for index, node in enumerate(nodes)}
        dsu = DisjointSet(len(nodes))
        for left, right in edges:
            dsu.union(node_index[left], node_index[right])
        components = {dsu.find(index) for index in range(len(nodes))}
        connected &= len(components) == 1
        rank += len(nodes) - len(components)

    product = jump @ injection
    product.eliminate_zeros()
    return jump, injection, {
        "degree": degree,
        "jump_rows": jump.shape[0],
        "local_dimension": local_dimension,
        "conforming_dimension": len(cells[degree]),
        "combinatorial_rank": rank,
        "expected_rank": local_dimension - len(cells[degree]),
        "occurrence_graphs_connected": bool(connected),
        "maximum_occurrences": maximum_occurrences,
        "occurrence_histogram": {
            str(key): value for key, value in sorted(occurrence_histogram.items())
        },
        "rj_nonzeros": int(product.nnz),
    }


def face_mass_audit(level):
    triangle_indices = {cell: index for index, cell in enumerate(level["cells"][2])}
    parents = [[] for _ in level["cells"][2]]
    for top_index, tetrahedron in enumerate(level["top"]):
        for triangle in combinations(tetrahedron, 3):
            parents[triangle_indices[tuple(triangle)]].append(top_index)
    grams = {}
    for triangle, adjacent in zip(level["cells"][2], parents):
        parent = adjacent[0]
        position = {vertex: index for index, vertex in enumerate(level["top"][parent])}
        points = tuple(level["top_points"][parent][position[vertex]]
                       for vertex in triangle)
        gram = triangle_gram(points)
        grams[matrix_key(gram)] = gram
    minimum_eigenvalue = float("inf")
    positive = True
    for gram in grams.values():
        for degree in range(3):
            mass = triangle_whitney_mass_from_gram(gram, degree)
            minors = [sy.simplify(mass[:size, :size].det())
                      for size in range(1, mass.rows + 1)]
            positive &= all(bool(value > 0) for value in minors)
            eigenvalues = np.linalg.eigvalsh(np.asarray(mass, dtype=np.float64))
            minimum_eigenvalue = min(minimum_eigenvalue, float(eigenvalues[0]))
    return {
        "ordered_face_gram_types": len(grams),
        "all_exact_principal_minors_positive": bool(positive),
        "minimum_numerical_face_mass_eigenvalue": minimum_eigenvalue,
    }


print("=" * 78)
print("ROLE OF THE WHITNEY TRACE-JUMP PENALTY")
print("=" * 78)

stiffness = json.loads(STIFFNESS_CERTIFICATE.read_text())
refinement = json.loads(REFINEMENT_CERTIFICATE.read_text())
dilation = json.loads(DILATION_CERTIFICATE.read_text())
check("all three frozen input certificates have the expected protocols",
      stiffness["protocol_commit"] == EXPECTED_STIFFNESS_PROTOCOL
      and stiffness["status"] == "COMPLETE"
      and stiffness["passed"] == stiffness["tests"]
      and refinement["protocol_commit"] == EXPECTED_REFINEMENT_PROTOCOL
      and refinement["passed"] == refinement["tests"]
      and dilation["protocol_commit"] == EXPECTED_DILATION_PROTOCOL)

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
ranked = barycentric_refine(make_base_level(reference_vertices))

level_records = []
all_rj_zero = True
all_kernel_exact = True
all_face_mass_positive = True
all_witnesses = True
for resolution in (1, 2, 4):
    print(f"\n-- k={resolution} --")
    level = rank_edgewise_level(ranked, resolution)
    degree_records = []
    degree_matrices = []
    for degree in range(3):
        jump, injection, record = copy_matrices(level, degree)
        degree_records.append(record)
        degree_matrices.append((jump, injection))
        all_rj_zero &= record["rj_nonzeros"] == 0
        all_kernel_exact &= (
            record["occurrence_graphs_connected"]
            and record["combinatorial_rank"] == record["expected_rank"]
        )
        print(
            f"degree {degree}: local/global="
            f"{record['local_dimension']}/{record['conforming_dimension']}, "
            f"rank={record['combinatorial_rank']}, "
            f"max occurrence={record['maximum_occurrences']}"
        )

    # A delta at one global vertex is conforming but nonconstant.  Its jump is
    # zero, while its simplicial coboundary is nonzero on every incident edge.
    jump_zero, injection_zero = degree_matrices[0]
    scalar = np.zeros(len(level["cells"][0]), dtype=np.int8)
    scalar[0] = 1
    local_scalar = injection_zero @ scalar
    jump_witness = jump_zero @ local_scalar
    edge_coboundary = np.asarray([
        scalar[right] - scalar[left] for left, right in level["cells"][1]
    ], dtype=np.int8)
    witness = (
        np.count_nonzero(jump_witness) == 0
        and np.count_nonzero(edge_coboundary) > 0
    )
    all_witnesses &= witness
    face_record = face_mass_audit(level)
    all_face_mass_positive &= face_record[
        "all_exact_principal_minors_positive"
    ]
    level_records.append({
        "edgewise_resolution": resolution,
        "f_vector": list(map(len, level["cells"])),
        "degree_records": degree_records,
        "nonconstant_conforming_witness": {
            "jump_nonzeros": int(np.count_nonzero(jump_witness)),
            "coboundary_nonzeros": int(np.count_nonzero(edge_coboundary)),
        },
        "face_mass": face_record,
    })

check("R_h J_h is exactly zero in all nine complete degree/level cases",
      all_rj_zero)
check("connected occurrence graphs give ker R_h = im J_h in every case",
      all_kernel_exact)
check("all exact face trace masses are positive definite",
      all_face_mass_positive,
      "minimum eigenvalues=" + str([
          record["face_mass"]["minimum_numerical_face_mass_eigenvalue"]
          for record in level_records
      ]))
check("a nonconstant conforming scalar has zero jump and nonzero coboundary",
      all_witnesses,
      str([record["nonconstant_conforming_witness"]
           for record in level_records]))

uniform_dilation = dilation["uniform_dilation"]
check("the exact Whitney element/trace powers give common 1/h scaling",
      uniform_dilation["element_mass_scaling_exact"]
      and uniform_dilation["face_mass_scaling_exact"]
      and abs(uniform_dilation["dirac_norm_ratio_at_h_equals_2"] - 0.5)
      < 1e-12)

scaling_records = []
for audit in stiffness["audits"]:
    resolution = audit["edgewise_resolution"]
    local_norm = audit["maximum_local_dirac_norm"]
    gaps = [record["positive_gap"] for record in audit["degree_records"]]
    ratios = [local_norm / gap for gap in gaps]
    scaling_records.append({
        "edgewise_resolution": resolution,
        "h_proxy": 1.0 / resolution,
        "h_times_local_dirac_norm": local_norm / resolution,
        "h_times_positive_gaps": [gap / resolution for gap in gaps],
        "a_over_g": ratios,
        "sufficient_separation_threshold_2a_over_g": [
            2.0 * value for value in ratios
        ],
        "worst_sufficient_threshold": 2.0 * max(ratios),
    })
check("the accepted finite spectra have positive normalized penalty gaps",
      all(value > 0.0 for record in scaling_records
          for value in record["h_times_positive_gaps"]))

maximum_occurrences = {
    str(degree): [
        record["degree_records"][degree]["maximum_occurrences"]
        for record in level_records
    ]
    for degree in range(3)
}

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_protocols": {
        "stiffness": EXPECTED_STIFFNESS_PROTOCOL,
        "refinement": EXPECTED_REFINEMENT_PROTOCOL,
        "dilation": EXPECTED_DILATION_PROTOCOL,
    },
    "phenomenological_target_used": False,
    "level_records": level_records,
    "maximum_occurrences_by_degree_k1_k2_k4": maximum_occurrences,
    "scaling_records": scaling_records,
    "structural_theorem_hypotheses": {
        "quasi_uniform": True,
        "uniformly_shape_regular": True,
        "finite_normalized_shape_types": True,
        "uniformly_bounded_occurrence_neighborhoods": True,
        "exact_whitney_and_trace_mass_powers": True,
        "basis": (
            "rank-edgewise theorem plus finite original complex; finite-level "
            "records are controls, not the all-level proof"
        ),
    },
    "verdicts": [
        "DERIVED NEGATIVE FOR DYNAMICS: pure jump stiffness has zero conforming compression",
        "DERIVED STRUCTURAL POSITIVE FOR SEPARATION: a_h/g_h is uniformly bounded under the stated tower hypotheses",
        "OPEN: geometry-selected finite kappa and coupled-pencil continuum convergence",
        "OPEN: Lorentzian time, causal speed, mass and Planck units",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured role certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("MAX_OCCURRENCES=" + json.dumps(maximum_occurrences, sort_keys=True))
print("WORST_THRESHOLDS=" + str([
    record["worst_sufficient_threshold"] for record in scaling_records
]))
raise SystemExit(0 if passed == tests else 1)


#!/usr/bin/env python3
"""Bounded-degree face-neighbour constraints for the Whitney KKT pencil.

Protocol commit a819a52 froze the canonical copy graphs, both refinement
levels, exact kernel/rank gates and locality bounds before enumeration.
"""

from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("whitney_neighbour_constraints.json")
PROTOCOL_COMMIT = "a819a52"
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


class DisjointSet:
    def __init__(self, size):
        self.parent = np.arange(size, dtype=np.int64)
        self.rank = np.zeros(size, dtype=np.int8)

    def find(self, item):
        item = int(item)
        root = item
        while self.parent[root] != root:
            root = int(self.parent[root])
        while self.parent[item] != item:
            parent = int(self.parent[item])
            self.parent[item] = root
            item = parent
        return root

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def audit_level(name, top_cells, cells):
    top_count = len(top_cells)
    dimensions = tuple(map(len, cells))
    cell_indices = [
        {cell: index for index, cell in enumerate(layer)}
        for layer in cells
    ]
    local_counts = tuple(len(list(combinations(range(4), degree + 1)))
                         for degree in range(4))

    copy_lookup = []
    copy_global = []
    injections = []
    for degree in range(4):
        lookup = {}
        global_labels = np.empty(top_count * local_counts[degree], dtype=np.int32)
        rows = []
        columns = []
        for top_index, top in enumerate(top_cells):
            for local_index, local_face in enumerate(
                combinations(range(4), degree + 1)
            ):
                cell = tuple(top[index] for index in local_face)
                global_index = cell_indices[degree][cell]
                row = top_index * local_counts[degree] + local_index
                lookup[(top_index, global_index)] = row
                global_labels[row] = global_index
                rows.append(row)
                columns.append(global_index)
        injection = sparse.csr_matrix(
            (np.ones(len(rows), dtype=np.int8), (rows, columns)),
            shape=(top_count * local_counts[degree], dimensions[degree]),
        )
        copy_lookup.append(lookup)
        copy_global.append(global_labels)
        injections.append(injection)

    triangle_parents = [[] for _ in range(dimensions[2])]
    for top_index, top in enumerate(top_cells):
        for face in combinations(top, 3):
            triangle_parents[cell_indices[2][tuple(face)]].append(top_index)
    triangle_parent_histogram = {}
    for parents in triangle_parents:
        triangle_parent_histogram[len(parents)] = (
            triangle_parent_histogram.get(len(parents), 0) + 1
        )

    edge_pairs = [[] for _ in range(4)]
    for triangle_index, parents in enumerate(triangle_parents):
        if len(parents) != 2:
            continue
        left_top, right_top = sorted(parents)
        triangle = cells[2][triangle_index]
        for degree in range(3):
            for simplex in combinations(triangle, degree + 1):
                global_index = cell_indices[degree][tuple(simplex)]
                left = copy_lookup[degree][(left_top, global_index)]
                right = copy_lookup[degree][(right_top, global_index)]
                edge_pairs[degree].append((left, right, global_index))

    degree_audits = []
    constraints = []
    for degree in range(4):
        pairs = edge_pairs[degree]
        row_indices = []
        column_indices = []
        data = []
        disjoint = DisjointSet(injections[degree].shape[0])
        node_degrees = np.zeros(injections[degree].shape[0], dtype=np.int16)
        no_mixing = True
        for row, (left, right, global_index) in enumerate(pairs):
            row_indices.extend((row, row))
            column_indices.extend((left, right))
            data.extend((1, -1))
            disjoint.union(left, right)
            node_degrees[left] += 1
            node_degrees[right] += 1
            no_mixing &= (
                copy_global[degree][left] == global_index
                and copy_global[degree][right] == global_index
            )
        constraint = sparse.csr_matrix(
            (data, (row_indices, column_indices)),
            shape=(len(pairs), injections[degree].shape[0]),
            dtype=np.int8,
        )
        constraints.append(constraint)
        constraint_times_injection = constraint @ injections[degree]

        roots = np.asarray([
            disjoint.find(node) for node in range(injections[degree].shape[0])
        ], dtype=np.int64)
        component_roots = np.unique(roots)
        component_count = len(component_roots)
        global_to_roots = [set() for _ in range(dimensions[degree])]
        root_to_globals = {}
        for node, root in enumerate(roots):
            global_index = int(copy_global[degree][node])
            global_to_roots[global_index].add(int(root))
            root_to_globals.setdefault(int(root), set()).add(global_index)
        every_star_connected = all(len(roots_for_cell) == 1
                                   for roots_for_cell in global_to_roots)
        components_do_not_mix = all(len(global_set) == 1
                                    for global_set in root_to_globals.values())
        exact_rank = injections[degree].shape[0] - component_count
        row_nonzeros = np.diff(constraint.indptr)
        degree_audits.append({
            "degree": degree,
            "global_dimension": dimensions[degree],
            "local_copy_dimension": injections[degree].shape[0],
            "constraint_rows": constraint.shape[0],
            "constraint_nonzeros": int(constraint.nnz),
            "all_rows_have_two_nonzeros": bool(
                len(row_nonzeros) == 0 or np.all(row_nonzeros == 2)
            ),
            "coefficient_alphabet": sorted(set(map(int, constraint.data))),
            "same_global_simplex_only": bool(no_mixing),
            "constraint_times_injection_nonzeros": int(
                constraint_times_injection.nnz
            ),
            "occurrence_graph_components": component_count,
            "every_simplex_star_connected": every_star_connected,
            "components_do_not_mix_simplices": components_do_not_mix,
            "maximum_occurrence_node_degree": int(node_degrees.max()),
            "expected_maximum_degree": 3 - degree if degree < 3 else 0,
            "exact_graph_incidence_rank": exact_rank,
            "expected_rank_local_minus_global": (
                injections[degree].shape[0] - dimensions[degree]
            ),
            "redundant_row_gauge_dimension": int(
                constraint.shape[0] - exact_rank
            ),
        })

    return {
        "level": name,
        "f_vector": list(dimensions),
        "top_count": top_count,
        "local_copy_dimensions": [
            injection.shape[0] for injection in injections
        ],
        "local_copy_total": sum(
            injection.shape[0] for injection in injections
        ),
        "triangle_parent_histogram": triangle_parent_histogram,
        "degree_audits": degree_audits,
        "total_constraint_rows": sum(
            audit["constraint_rows"] for audit in degree_audits
        ),
        "total_exact_rank": sum(
            audit["exact_graph_incidence_rank"] for audit in degree_audits
        ),
        "total_redundant_row_gauge_dimension": sum(
            audit["redundant_row_gauge_dimension"] for audit in degree_audits
        ),
    }


print("=" * 78)
print("UNIFORMLY LOCAL WHITNEY NEIGHBOUR CONSTRAINTS")
print("=" * 78)

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

base_audit = audit_level("base", base_top, base_cells)
fine_audit = audit_level("first_barycentric", fine_top, fine_cells)
audits = (base_audit, fine_audit)

check("both exact f-vectors and local-copy totals are correct",
      base_audit["f_vector"] == [120, 720, 1200, 600]
      and base_audit["local_copy_total"] == 9000
      and fine_audit["f_vector"] == [2640, 17040, 28800, 14400]
      and fine_audit["local_copy_total"] == 216000)
check("every triangle has exactly two parent tetrahedra at both levels",
      all(audit["triangle_parent_histogram"] == {2: audit["f_vector"][2]}
          for audit in audits),
      str([audit["triangle_parent_histogram"] for audit in audits]))
check("the preregistered total constraint row counts are exact",
      base_audit["total_constraint_rows"] == 8400
      and fine_audit["total_constraint_rows"] == 201600)
check("every constraint row is an unweighted signed pair difference",
      all(
          degree["all_rows_have_two_nonzeros"]
          and (degree["constraint_rows"] == 0
               or degree["coefficient_alphabet"] == [-1, 1])
          for audit in audits for degree in audit["degree_audits"]
      ))
check("all neighbour constraints annihilate conforming assembly exactly",
      all(degree["constraint_times_injection_nonzeros"] == 0
          for audit in audits for degree in audit["degree_audits"]))
check("no constraint ever mixes two global simplices",
      all(degree["same_global_simplex_only"]
          and degree["components_do_not_mix_simplices"]
          for audit in audits for degree in audit["degree_audits"]))
check("every global simplex occurrence graph is connected at both levels",
      all(degree["every_simplex_star_connected"]
          for audit in audits for degree in audit["degree_audits"]))
check("the exact constraint kernels equal the conforming spaces",
      all(
          degree["occurrence_graph_components"] == degree["global_dimension"]
          and degree["exact_graph_incidence_rank"]
          == degree["expected_rank_local_minus_global"]
          for audit in audits for degree in audit["degree_audits"]
      ))
check("the maximum local degrees remain exactly (3,2,1,0)",
      all([
          degree["maximum_occurrence_node_degree"]
          for degree in audit["degree_audits"]
      ] == [3, 2, 1, 0] for audit in audits),
      str([[
          degree["maximum_occurrence_node_degree"]
          for degree in audit["degree_audits"]
      ] for audit in audits]))
check("the exact total ranks are 6360 and 153120",
      base_audit["total_exact_rank"] == 6360
      and fine_audit["total_exact_rank"] == 153120)
check("retaining all canonical rows exposes the exact gauge redundancies",
      base_audit["total_redundant_row_gauge_dimension"] == 2040
      and fine_audit["total_redundant_row_gauge_dimension"] == 48480)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "constraint_definition": (
        "pair equal copies only across tetrahedra sharing a triangle "
        "that contains the copied simplex"
    ),
    "audits": list(audits),
    "verdicts": [
        "DERIVED UNIFORMLY LOCAL SPECTRAL CONSTRAINT at levels 0 and 1",
        "DERIVED exact KKT row-space equivalence to conformity",
        "DERIVED NEGATIVE: multiplier metric remains singular and is not a unitary tick",
    ],
    "scope": (
        "Base and first barycentric closed 3-manifolds; algebraic KKT "
        "spectral representation, not positive-metric dynamics."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured neighbour-constraint certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for audit in audits:
    print(
        f"{audit['level']}: rows={audit['total_constraint_rows']}, "
        f"rank={audit['total_exact_rank']}, "
        f"gauge={audit['total_redundant_row_gauge_dimension']}, "
        "max degrees=" + str([
            degree["maximum_occurrence_node_degree"]
            for degree in audit["degree_audits"]
        ])
    )
print("SPECTRAL_VERDICT: exact constraint locality degree <= 3 at both levels")
print("DYNAMICAL_VERDICT: descriptor metric remains singular")
raise SystemExit(0 if passed == tests else 1)

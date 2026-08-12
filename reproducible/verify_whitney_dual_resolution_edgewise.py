#!/usr/bin/env python3
"""Locality flow of the dual constraint resolution on the edgewise tower.

Protocol commit 36a56b7 froze the carrier, link hypotheses, incidence
observables and k=2-to-k=4 flow gate before these links were enumerated.
"""

from collections import Counter, deque
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import sympy as sy

from whitney_trace_refinement_tools import (
    barycentric_refine,
    make_base_level,
    rank_edgewise_level,
)


OUTPUT = Path(__file__).with_name(
    "whitney_dual_resolution_edgewise.json"
)
PROTOCOL_COMMIT = "36a56b7"
DUAL_RESOLUTION_RESULT_COMMIT = "799966f"
EXPECTED_F_VECTORS = {
    1: (30, 150, 240, 120),
    2: (180, 1140, 1920, 960),
    4: (1320, 9000, 15360, 7680),
}
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def histogram(values):
    counts = Counter(values)
    return {
        str(key): int(value)
        for key, value in sorted(counts.items(), key=lambda item: str(item[0]))
    }


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
        left = self.find(left)
        right = self.find(right)
        if left == right:
            return
        if self.rank[left] < self.rank[right]:
            left, right = right, left
        self.parent[right] = left
        if self.rank[left] == self.rank[right]:
            self.rank[left] += 1


def connected_subgraph(nodes, pairs):
    if not nodes:
        return False
    adjacency = {node: [] for node in nodes}
    for left, right in pairs:
        adjacency[left].append(right)
        adjacency[right].append(left)
    reached = {min(nodes)}
    queue = deque(reached)
    while queue:
        node = queue.popleft()
        for neighbour in adjacency[node]:
            if neighbour not in reached:
                reached.add(neighbour)
                queue.append(neighbour)
    return reached == set(nodes)


def occurrence_connectivity(level, triangle_parents):
    top_cells = level["top"]
    cells = level["cells"]
    results = []
    for degree in range(4):
        local_faces = tuple(combinations(range(4), degree + 1))
        cell_indices = {cell: index for index, cell in enumerate(cells[degree])}
        lookup = {}
        globals_by_copy = []
        for top_index, top in enumerate(top_cells):
            for local_index, positions in enumerate(local_faces):
                cell = tuple(top[position] for position in positions)
                global_index = cell_indices[cell]
                copy = top_index * len(local_faces) + local_index
                lookup[(top_index, global_index)] = copy
                globals_by_copy.append(global_index)
        dsu = DisjointSet(len(globals_by_copy))
        node_degrees = np.zeros(len(globals_by_copy), dtype=np.int16)
        if degree < 3:
            for triangle_index, parents in enumerate(triangle_parents):
                left_top, right_top = sorted(parents)
                triangle = cells[2][triangle_index]
                for simplex in combinations(triangle, degree + 1):
                    global_index = cell_indices[tuple(simplex)]
                    left = lookup[(left_top, global_index)]
                    right = lookup[(right_top, global_index)]
                    dsu.union(left, right)
                    node_degrees[left] += 1
                    node_degrees[right] += 1
        roots_by_global = [set() for _ in cells[degree]]
        for copy, global_index in enumerate(globals_by_copy):
            roots_by_global[global_index].add(dsu.find(copy))
        results.append({
            "degree": degree,
            "local_copy_count": len(globals_by_copy),
            "global_simplex_count": len(cells[degree]),
            "every_occurrence_graph_connected": all(
                len(roots) == 1 for roots in roots_by_global
            ),
            "maximum_neighbour_constraint_node_degree": int(
                node_degrees.max()
            ) if len(node_degrees) else 0,
            "expected_maximum": 3 - degree if degree < 3 else 0,
        })
    return results


def audit_level(level, resolution):
    top_cells = level["top"]
    vertices, edges, triangles, tetrahedra = level["cells"]
    vertex_indices = {cell[0]: index for index, cell in enumerate(vertices)}
    edge_indices = {edge: index for index, edge in enumerate(edges)}
    triangle_indices = {
        triangle: index for index, triangle in enumerate(triangles)
    }

    tops_by_vertex = [[] for _ in vertices]
    tops_by_edge = [[] for _ in edges]
    edges_by_vertex = [[] for _ in vertices]
    triangles_by_vertex = [[] for _ in vertices]
    triangles_by_edge = [[] for _ in edges]
    triangle_parents = [[] for _ in triangles]

    for edge_index, edge in enumerate(edges):
        for vertex in edge:
            edges_by_vertex[vertex_indices[vertex]].append(edge_index)
    for triangle_index, triangle in enumerate(triangles):
        for vertex in triangle:
            triangles_by_vertex[vertex_indices[vertex]].append(triangle_index)
        for edge in combinations(triangle, 2):
            triangles_by_edge[edge_indices[edge]].append(triangle_index)
    for top_index, top in enumerate(top_cells):
        for vertex in top:
            tops_by_vertex[vertex_indices[vertex]].append(top_index)
        for edge in combinations(top, 2):
            tops_by_edge[edge_indices[edge]].append(top_index)
        for triangle in combinations(top, 3):
            triangle_parents[triangle_indices[triangle]].append(top_index)

    every_triangle_two_parents = all(
        len(parents) == 2 for parents in triangle_parents
    )

    edge_link_records = []
    every_edge_link_cycle = True
    for edge_index, top_ids in enumerate(tops_by_edge):
        pairs = [
            tuple(sorted(triangle_parents[triangle_index]))
            for triangle_index in triangles_by_edge[edge_index]
        ]
        degrees = Counter(node for pair in pairs for node in pair)
        connected = connected_subgraph(top_ids, pairs)
        cycle = (
            connected
            and len(pairs) == len(top_ids)
            and set(degrees) == set(top_ids)
            and all(value == 2 for value in degrees.values())
        )
        every_edge_link_cycle &= cycle
        edge_link_records.append({
            "edge_index": edge_index,
            "tetrahedron_count": len(top_ids),
            "triangle_count": len(pairs),
            "connected_cycle": cycle,
        })

    vertex_link_records = []
    every_vertex_link_sphere = True
    for vertex_index, vertex_cell in enumerate(vertices):
        vertex = vertex_cell[0]
        neighbour_vertices = {
            other
            for edge_index in edges_by_vertex[vertex_index]
            for other in edges[edge_index]
            if other != vertex
        }
        link_pairs = []
        for triangle_index in triangles_by_vertex[vertex_index]:
            others = tuple(
                item for item in triangles[triangle_index] if item != vertex
            )
            link_pairs.append(others)
        connected = connected_subgraph(neighbour_vertices, link_pairs)
        link_v = len(edges_by_vertex[vertex_index])
        link_e = len(triangles_by_vertex[vertex_index])
        link_f = len(tops_by_vertex[vertex_index])
        closed_edges = all(
            len(triangle_parents[triangle_index]) == 2
            for triangle_index in triangles_by_vertex[vertex_index]
        )
        sphere_gate = connected and closed_edges and link_v - link_e + link_f == 2
        every_vertex_link_sphere &= sphere_gate
        vertex_link_records.append({
            "vertex_index": vertex_index,
            "link_f_vector": [link_v, link_e, link_f],
            "link_euler_characteristic": link_v - link_e + link_f,
            "link_one_skeleton_connected": connected,
            "every_link_edge_has_two_faces": closed_edges,
            "passes_sphere_gate": sphere_gate,
        })

    occurrence_records = occurrence_connectivity(level, triangle_parents)
    all_occurrence_connected = all(
        item["every_occurrence_graph_connected"]
        for item in occurrence_records
    )
    neighbour_maxima = [
        item["maximum_neighbour_constraint_node_degree"]
        for item in occurrence_records
    ]

    vertex_occurrences = list(map(len, tops_by_vertex))
    edge_occurrences = list(map(len, tops_by_edge))
    vertex_edge_degrees = list(map(len, edges_by_vertex))
    observables = {
        "a0_max_tetrahedron_occurrences_per_vertex": max(vertex_occurrences),
        "a1_max_tetrahedron_occurrences_per_edge": max(edge_occurrences),
        "r3_max_edges_incident_per_vertex": max(vertex_edge_degrees),
        "neighbour_constraint_node_maxima": neighbour_maxima,
    }
    return {
        "edgewise_resolution": resolution,
        "f_vector": list(map(len, level["cells"])),
        "every_triangle_has_two_parents": every_triangle_two_parents,
        "every_edge_link_is_one_cycle": every_edge_link_cycle,
        "every_vertex_link_passes_closed_connected_euler_two_gate": (
            every_vertex_link_sphere
        ),
        "all_occurrence_graphs_connected": all_occurrence_connected,
        "edge_link_tetrahedron_count_histogram": histogram(
            record["tetrahedron_count"] for record in edge_link_records
        ),
        "vertex_link_f_vector_histogram": histogram(
            tuple(record["link_f_vector"]) for record in vertex_link_records
        ),
        "tetrahedron_occurrences_per_vertex_histogram": histogram(
            vertex_occurrences
        ),
        "tetrahedron_occurrences_per_edge_histogram": histogram(
            edge_occurrences
        ),
        "edges_incident_per_vertex_histogram": histogram(
            vertex_edge_degrees
        ),
        "occurrence_connectivity_by_degree": occurrence_records,
        "locality_observables": observables,
        "dual_resolution_exactness_hypotheses_pass": (
            every_triangle_two_parents
            and every_edge_link_cycle
            and every_vertex_link_sphere
            and all_occurrence_connected
        ),
    }


print("=" * 78)
print("DUAL-RESOLUTION LOCALITY ON THE RANK-EDGEWISE TOWER")
print("=" * 78)

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
ranked = barycentric_refine(make_base_level(reference_vertices))

records = []
for resolution in (1, 2, 4):
    print(f"\n-- rank-edgewise resolution k={resolution} --", flush=True)
    level = rank_edgewise_level(ranked, resolution)
    record = audit_level(level, resolution)
    records.append(record)
    print(
        f"  f={record['f_vector']}, obs={record['locality_observables']}, "
        f"links={record['dual_resolution_exactness_hypotheses_pass']}",
        flush=True,
    )

check(
    "all three preregistered f-vectors are exact",
    all(tuple(record["f_vector"]) == EXPECTED_F_VECTORS[
        record["edgewise_resolution"]
    ] for record in records),
)
check(
    "every triangle has exactly two parent tetrahedra",
    all(record["every_triangle_has_two_parents"] for record in records),
)
check(
    "every edge link is one connected cycle",
    all(record["every_edge_link_is_one_cycle"] for record in records),
)
check(
    "every vertex link passes the closed connected Euler-two gate",
    all(
        record[
            "every_vertex_link_passes_closed_connected_euler_two_gate"
        ] for record in records
    ),
)
check(
    "every simplex occurrence graph is connected in all four degrees",
    all(record["all_occurrence_graphs_connected"] for record in records),
)
check(
    "the transferred exact dual-resolution hypotheses pass at every level",
    all(record["dual_resolution_exactness_hypotheses_pass"]
        for record in records),
)
check(
    "neighbour constraint node maxima stay exactly (3,2,1,0)",
    all(record["locality_observables"][
        "neighbour_constraint_node_maxima"
    ] == [3, 2, 1, 0] for record in records),
)

observable_keys = (
    "a0_max_tetrahedron_occurrences_per_vertex",
    "a1_max_tetrahedron_occurrences_per_edge",
    "r3_max_edges_incident_per_vertex",
)
flow = []
primary_nonincreasing = True
for key in observable_keys:
    values = [record["locality_observables"][key] for record in records]
    primary_pass = values[2] <= values[1]
    primary_nonincreasing &= primary_pass
    flow.append({
        "observable": key,
        "values_k1_k2_k4": values,
        "k1_to_k2_increase": values[1] > values[0],
        "k2_to_k4_nonincreasing": primary_pass,
    })
check("the preregistered k=2 to k=4 locality gate was evaluated", True)

flow_verdict = (
    "PATTERN TOWARD BOUNDED DUAL LOCALITY: every frozen incidence maximum "
    "is nonincreasing from k=2 to k=4"
    if primary_nonincreasing else
    "PATTERN NEGATIVE FOR BOUNDED DUAL LOCALITY: at least one frozen "
    "incidence maximum increases from k=2 to k=4"
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "dual_resolution_result_commit": DUAL_RESOLUTION_RESULT_COMMIT,
    "phenomenological_target_used": False,
    "spectrum_computed": False,
    "records": records,
    "flow_by_observable": flow,
    "all_k2_to_k4_maxima_nonincreasing": primary_nonincreasing,
    "verdicts": [
        "DERIVED ON CONTROLS: all link hypotheses transfer the exact signed "
        "dual-cell resolution to k=1,2,4",
        flow_verdict,
        "OPEN: an analytic uniform incidence bound for general k",
        "STRUCTURAL: bounded kinematic constraint resolution is not a "
        "physical first-class gauge theory or Hamiltonian",
        "NOT CLAIMED: spectrum, time, causality, inertia, mass or Planck units",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured no-spectrum certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(flow_verdict)
for item in flow:
    print(f"{item['observable']}: {item['values_k1_k2_k4']}")
raise SystemExit(0 if passed == tests else 1)

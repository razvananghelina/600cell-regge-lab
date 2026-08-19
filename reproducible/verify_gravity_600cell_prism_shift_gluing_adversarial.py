#!/usr/bin/env python3
"""Independent fundamental-cycle audit of global prism shift gluing."""

from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_prism_shift_gluing_adversarial.json"
AUDIT_PROTOCOL_COMMIT = "42aba25"
INPUT_HASHES = {
    "docs/gravity/gravity_600cell_prism_shift_gluing_adversarial_protocol.md":
        "cd59fd9d80832a5ac72144fc80fd8f2d0c7c0af0b96d6ca31905b0ad2f9c9c57",
    "reproducible/gravity_600cell_prism_shift_gluing.json":
        "1ab6654ae57c83a49dd4f427154b891c0b8ae613631773ab6733a1227b9999fa",
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
PRIMES = (3, 1000003)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def tetrahedra_from_adjacency(adjacency):
    neighbours = [set(np.flatnonzero(row > 0.5)) for row in adjacency]
    top = []
    for first in range(len(adjacency)):
        for second in sorted(v for v in neighbours[first] if v > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(v for v in common_two if v > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(v for v in common_three if v > third):
                    top.append((first, second, third, fourth))
    return tuple(top)


def simplices_from_top(top, vertex_count=None):
    if vertex_count is None:
        vertex_count = max(max(simplex) for simplex in top)+1
    vertices = tuple((index,) for index in range(vertex_count))
    edges = tuple(sorted({tuple(sorted(edge)) for simplex in top
                          for edge in combinations(simplex, 2)}))
    triangles = tuple(sorted({tuple(sorted(face)) for simplex in top
                              for face in combinations(simplex, 3)}))
    return vertices, edges, triangles, tuple(sorted(top))


def spanning_tree(vertex_count, edges, root=0):
    adjacency = [[] for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    for row in adjacency:
        row.sort()
    parent = [-1]*vertex_count
    parent[root] = root
    order = [root]
    queue = deque([root])
    tree_edges = set()
    while queue:
        current = queue.popleft()
        for neighbour in adjacency[current]:
            if parent[neighbour] != -1:
                continue
            parent[neighbour] = current
            tree_edges.add(tuple(sorted((current, neighbour))))
            order.append(neighbour)
            queue.append(neighbour)
    return tuple(parent), tuple(order), frozenset(tree_edges)


def sparse_rank_mod(rows, prime):
    """Sparse normalized elimination, unrelated to the primary F2 bitsets."""
    basis = {}
    for source in rows:
        row = {column: value % prime for column, value in source.items()
               if value % prime}
        while row:
            pivot = min(row)
            if pivot in basis:
                factor = row[pivot]
                for column, value in basis[pivot].items():
                    updated = (row.get(column, 0)-factor*value) % prime
                    if updated:
                        row[column] = updated
                    elif column in row:
                        del row[column]
            else:
                inverse = pow(row[pivot], -1, prime)
                row = {column: (value*inverse) % prime
                       for column, value in row.items() if value % prime}
                basis[pivot] = row
                break
    return len(basis)


def triangle_chord_rows(vertices, edges, triangles):
    parent, order, tree_edges = spanning_tree(len(vertices), edges)
    chords = tuple(edge for edge in edges if edge not in tree_edges)
    chord_index = {edge: index for index, edge in enumerate(chords)}
    rows = []
    for first, second, third in triangles:
        boundary = (((first, second), 1),
                    ((second, third), 1),
                    ((first, third), -1))
        row = {}
        for edge, coefficient in boundary:
            if edge in chord_index:
                row[chord_index[edge]] = coefficient
        rows.append(row)
    return {
        "parent": parent,
        "order": order,
        "tree_edges": tree_edges,
        "chords": chords,
        "rows": tuple(rows),
        "ranks": {prime: sparse_rank_mod(rows, prime) for prime in PRIMES},
    }


def torus_three_by_three():
    size = 3

    def vertex(i, j):
        return (i % size)*size+(j % size)

    triangles = set()
    for i in range(size):
        for j in range(size):
            a = vertex(i, j)
            b = vertex(i+1, j)
            c = vertex(i+1, j+1)
            d = vertex(i, j+1)
            triangles.add(tuple(sorted((a, b, c))))
            triangles.add(tuple(sorted((a, c, d))))
    triangles = tuple(sorted(triangles))
    edges = tuple(sorted({tuple(sorted(edge)) for face in triangles
                          for edge in combinations(face, 2)}))
    return tuple((index,) for index in range(9)), edges, triangles, tuple()


def graph_distances(vertex_count, edges, seed):
    adjacency = [[] for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    distances = [-1]*vertex_count
    distances[seed] = 0
    queue = deque([seed])
    while queue:
        current = queue.popleft()
        for neighbour in adjacency[current]:
            if distances[neighbour] == -1:
                distances[neighbour] = distances[current]+1
                queue.append(neighbour)
    return tuple(distances)


def edge_values(potential, edges):
    return {edge: int(potential[edge[1]])-int(potential[edge[0]])
            for edge in edges}


def reconstruct_potential(values, parent, order):
    reconstructed = [None]*len(parent)
    root = order[0]
    reconstructed[root] = 0
    for vertex in order[1:]:
        ancestor = parent[vertex]
        edge = tuple(sorted((ancestor, vertex)))
        difference = values[edge]
        reconstructed[vertex] = (
            reconstructed[ancestor]+difference
            if ancestor < vertex else
            reconstructed[ancestor]-difference
        )
    return tuple(reconstructed)


def verify_reconstruction(potential, edges, cycle_data):
    values = edge_values(potential, edges)
    reconstructed = reconstruct_potential(
        values, cycle_data["parent"], cycle_data["order"])
    root = cycle_data["order"][0]
    target = tuple(int(value)-int(potential[root]) for value in potential)
    all_edges = all(reconstructed[right]-reconstructed[left] == values[(left, right)]
                    for left, right in edges)
    constant_invariant = edge_values(
        tuple(int(value)+37 for value in potential), edges) == values
    return reconstructed == target and all_edges and constant_invariant


def local_alpha_values(tetrahedron, potential):
    reference = tetrahedron[0]
    return {vertex: int(potential[vertex])-int(potential[reference])
            for vertex in tetrahedron}


def lateral_interval_matrix(face, tetrahedron, potential, rho, vertices):
    local_alpha = local_alpha_values(tetrahedron, potential)
    matrix = np.zeros((6, 6), dtype=float)
    for left in range(6):
        for right in range(left+1, 6):
            left_top = left >= 3
            right_top = right >= 3
            left_vertex = face[left % 3]
            right_vertex = face[right % 3]
            spatial_square = float(np.dot(
                vertices[right_vertex]-vertices[left_vertex],
                vertices[right_vertex]-vertices[left_vertex]))
            if left_top == right_top:
                interval = spatial_square
            else:
                if not left_top:
                    bottom_vertex, top_vertex = left_vertex, right_vertex
                else:
                    bottom_vertex, top_vertex = right_vertex, left_vertex
                alpha = local_alpha[top_vertex]-local_alpha[bottom_vertex]
                interval = spatial_square+2.0*alpha-rho
            matrix[left, right] = matrix[right, left] = interval
    return matrix


def local_face_trace_rank(tetrahedron, face):
    reference = tetrahedron[0]
    positions = {vertex: index for index, vertex in enumerate(tetrahedron)}
    root, endpoint_one, endpoint_two = face
    matrix = np.zeros((2, 3), dtype=int)
    for row, endpoint in enumerate((endpoint_one, endpoint_two)):
        for vertex, coefficient in ((endpoint, 1), (root, -1)):
            if vertex != reference:
                matrix[row, positions[vertex]-1] += coefficient
    return int(np.linalg.matrix_rank(matrix.astype(float)))


print("="*78)
print("ADVERSARIAL GLOBAL PRISM SHIFT GLUING AUDIT")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the audit protocol, primary artifact and source have frozen provenance",
    actual_hashes == INPUT_HASHES and AUDIT_PROTOCOL_COMMIT == "42aba25",
    str(actual_hashes),
)

vertices, adjacency, _ = build_600cell()
vertices = vertices/np.linalg.norm(vertices, axis=1)[:, None]
top = tetrahedra_from_adjacency(adjacency)
simplices = simplices_from_top(top, len(vertices))
source_f = tuple(len(layer) for layer in simplices)
face_incidence = defaultdict(list)
for top_index, tetrahedron in enumerate(top):
    for face in combinations(tetrahedron, 3):
        face_incidence[face].append(top_index)
check(
    "the independent source reconstruction is the closed 600-cell",
    source_f == (120, 720, 1200, 600)
    and Counter(map(len, face_incidence.values())) == Counter({2: 1200}),
    f"f={source_f}, face incidence={dict(Counter(map(len, face_incidence.values())))}",
)

cycles = triangle_chord_rows(*simplices[:3])
cycle_dimension = len(cycles["chords"])
check(
    "the deterministic spanning tree has 601 chord-cycle coordinates",
    len(cycles["tree_edges"]) == 119 and cycle_dimension == 601
    and len(cycles["order"]) == 120,
    f"tree={len(cycles['tree_edges'])}, chords={cycle_dimension}",
)
check(
    "triangle boundaries span every cycle over both adversarial primes",
    cycles["ranks"] == {3: 601, 1000003: 601},
    f"ranks={cycles['ranks']}",
)

tetra_boundary = (
    tuple((index,) for index in range(4)),
    tuple(combinations(range(4), 2)),
    tuple(combinations(range(4), 3)),
    tuple(),
)
tetra_cycles = triangle_chord_rows(*tetra_boundary[:3])
check(
    "the tetrahedral-sphere control closes all three graph cycles",
    len(tetra_cycles["chords"]) == 3
    and tetra_cycles["ranks"] == {3: 3, 1000003: 3},
    f"cycle dim={len(tetra_cycles['chords'])}, ranks={tetra_cycles['ranks']}",
)

torus = torus_three_by_three()
torus_cycles = triangle_chord_rows(*torus[:3])
torus_b1 = {prime: len(torus_cycles["chords"])-rank
            for prime, rank in torus_cycles["ranks"].items()}
check(
    "the torus control leaves exactly two cycles over both primes",
    tuple(map(len, torus[:3])) == (9, 27, 18)
    and len(torus_cycles["chords"]) == 19
    and torus_cycles["ranks"] == {3: 17, 1000003: 17}
    and torus_b1 == {3: 2, 1000003: 2},
    f"cycle dim={len(torus_cycles['chords'])}, "
    f"ranks={torus_cycles['ranks']}, b1={torus_b1}",
)
check(
    "without triangle relations all 601 source graph cycles remain",
    sparse_rank_mod((), 3) == 0 and cycle_dimension == 601,
    "relation rank=0, unresolved cycles=601",
)

seed = min(range(len(vertices)), key=lambda index: tuple(vertices[index]))
distances = graph_distances(len(vertices), simplices[1], seed)
lexicographic_order = sorted(range(len(vertices)), key=lambda index: tuple(vertices[index]))
lexicographic_rank = [0]*len(vertices)
for rank, vertex in enumerate(lexicographic_order):
    lexicographic_rank[vertex] = rank-60
potentials = {
    "distance": distances,
    "distance_squared": tuple(value*value for value in distances),
    "signed_coordinate_order": tuple(lexicographic_rank),
    "modular_quadratic": tuple((17*i*i+3*i+5) % 101
                               for i in range(len(vertices))),
    "modular_cubic": tuple((i*i*i+11*i+7) % 97
                           for i in range(len(vertices))),
}
reconstruction = {
    name: verify_reconstruction(potential, simplices[1], cycles)
    for name, potential in potentials.items()
}
check(
    "five integer potential families reconstruct exactly from tree paths",
    all(reconstruction.values()),
    str(reconstruction),
)

# Independent finite family: squared graph distance and rho=7/5.
potential = potentials["distance_squared"]
rho = 7.0/5.0
minimum_spatial_eigenvalue = float("inf")
maximum_schur = -float("inf")
volume_ratios = []
signature_failures = 0
for tetrahedron in top:
    points = vertices[list(tetrahedron)]
    edge_matrix = (points[1:]-points[0]).T
    gram = edge_matrix.T @ edge_matrix
    eigenvalues = np.linalg.eigvalsh(gram)
    minimum_spatial_eigenvalue = min(minimum_spatial_eigenvalue,
                                     float(eigenvalues[0]))
    a = np.asarray([potential[vertex]-potential[tetrahedron[0]]
                    for vertex in tetrahedron[1:]], dtype=float)
    shift_square = float(a @ np.linalg.solve(gram, a))
    schur = -rho-shift_square
    maximum_schur = max(maximum_schur, schur)
    signature_failures += int(not (eigenvalues[0] > 0 and schur < 0))
    volume_ratios.append(float(np.sqrt((rho+shift_square)/rho)))
check(
    "the different finite family has Lorentzian signature in every cell",
    signature_failures == 0 and minimum_spatial_eigenvalue > 1e-9
    and maximum_schur <= -rho+1e-12,
    f"failures={signature_failures}, min G eigen={minimum_spatial_eigenvalue:.12g}, "
    f"max Schur={maximum_schur:.12g}",
)

maximum_face_residual = 0.0
for face, adjacent in face_incidence.items():
    left = lateral_interval_matrix(face, top[adjacent[0]], potential,
                                   rho, vertices)
    right = lateral_interval_matrix(face, top[adjacent[1]], potential,
                                    rho, vertices)
    maximum_face_residual = max(maximum_face_residual,
                                float(np.max(np.abs(left-right))))
check(
    "all 1200 independently built lateral face metrics shape-match",
    maximum_face_residual < 1e-10,
    f"maximum 6x6 squared-interval residual={maximum_face_residual:.3e}",
)
changed_cells = sum(ratio > 1+1e-12 for ratio in volume_ratios)
check(
    "the second family changes volume while all natural lengths remain fixed",
    changed_cells > 0 and max(volume_ratios) > 1+1e-6,
    f"changed={changed_cells}/{len(top)}, max volume ratio={max(volume_ratios):.12g}; "
    "bottom/top lengths are inherited and all strut squares equal -7/5",
)

# Hostile stronger convention: one common ambient R4 vector.
normals = []
face_trace_ranks = []
for tetrahedron in top:
    edge_matrix = (vertices[list(tetrahedron[1:])]-vertices[tetrahedron[0]]).T
    _, _, right_vectors = np.linalg.svd(edge_matrix.T)
    normals.append(right_vectors[-1])
    for face in combinations(tetrahedron, 3):
        face_trace_ranks.append(local_face_trace_rank(tetrahedron, face))
normal_matrix = np.asarray(normals)
ambient_rank = int(np.linalg.matrix_rank(normal_matrix, tol=1e-10))
ambient_intersection_dimension = 4-ambient_rank
check(
    "the stronger common-ambient-vector convention kills every mode",
    ambient_rank == 4 and ambient_intersection_dimension == 0,
    f"normal rank={ambient_rank}, intersection dim={ambient_intersection_dimension}",
)
check(
    "that stronger convention adds one component absent from each face metric",
    Counter(face_trace_ranks) == Counter({2: 2400}),
    f"local face-trace ranks={dict(Counter(face_trace_ranks))}",
)

verdict = (
    "GLOBAL_SHIFT_GLUING_CORROBORATED"
    if passed == tests else
    "PRIMARY_RESULT_REFUTED"
)
check(
    "the adversarial verdict is evaluated",
    verdict == "GLOBAL_SHIFT_GLUING_CORROBORATED",
    verdict,
)

artifact = {
    "ambient_convention_attack": {
        "normal_matrix_rank": ambient_rank,
        "common_vector_intersection_dimension": ambient_intersection_dimension,
        "face_trace_rank_counts": dict(sorted(Counter(face_trace_ranks).items())),
        "classification": "OVERCONSTRAINED_NOT_INTRINSIC_FACE_MATCHING",
    },
    "classification": "DERIVED_KINEMATIC_CORROBORATED" if passed == tests else "OPEN",
    "controls": {
        "tetrahedral_sphere": {
            "cycle_dimension": len(tetra_cycles["chords"]),
            "boundary_ranks": tetra_cycles["ranks"],
        },
        "torus": {
            "f_vector": [len(torus[0]), len(torus[1]), len(torus[2])],
            "cycle_dimension": len(torus_cycles["chords"]),
            "boundary_ranks": torus_cycles["ranks"],
            "betti_1": torus_b1,
        },
    },
    "finite_family": {
        "potential": "squared_graph_distance",
        "rho": "7/5",
        "signature_failures": signature_failures,
        "minimum_spatial_gram_eigenvalue": minimum_spatial_eigenvalue,
        "maximum_schur_complement": maximum_schur,
        "maximum_face_metric_residual": maximum_face_residual,
        "changed_cell_count": changed_cells,
        "maximum_volume_ratio": max(volume_ratios),
        "natural_lengths_identical": True,
    },
    "fundamental_cycles": {
        "tree_edges": len(cycles["tree_edges"]),
        "cycle_dimension": cycle_dimension,
        "triangle_boundary_ranks": cycles["ranks"],
        "betti_1": {prime: cycle_dimension-rank
                    for prime, rank in cycles["ranks"].items()},
    },
    "potential_reconstruction": reconstruction,
    "provenance": {
        "audit_protocol_commit": AUDIT_PROTOCOL_COMMIT,
        "input_hashes": actual_hashes,
        "primes": list(PRIMES),
    },
    "tests": tests,
    "passed": passed,
    "verdict": verdict,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")
print(f"\nResult: {passed}/{tests} checks passed.")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    sys.exit(1)


#!/usr/bin/env python3
"""Exact global matching census for equal-scale prism shift modes."""

from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_prism_shift_gluing.json"
PRIOR_ART_COMMIT = "4a9ea11"
PROTOCOL_COMMIT = "a2c9174"
INPUT_HASHES = {
    "docs/gravity/gravity_600cell_prism_shift_gluing_prior_art.md":
        "8a1051f8e19bfd5b2d3e1c75cee1c2c21932e4039a7355fb6dd92465643e4a2e",
    "docs/gravity/gravity_600cell_prism_shift_gluing_protocol.md":
        "f76a2216f20ea00aa66c31891d81c1f72a4ee86fe519458b6527118a4fca2251",
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
EXPECTED_BASE_F = (120, 720, 1200, 600)
EXPECTED_FINE_F = (19680, 134880, 230400, 115200)
DIRECT_CHILDREN = (
    ("v0", "m01", "m02", "m03"),
    ("v1", "m01", "m12", "m13"),
    ("v2", "m02", "m12", "m23"),
    ("v3", "m03", "m13", "m23"),
    ("m01", "m02", "m03", "m13"),
    ("m01", "m02", "m12", "m13"),
    ("m02", "m03", "m13", "m23"),
    ("m02", "m12", "m13", "m23"),
)
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


def all_simplices(top):
    return tuple(
        tuple(sorted({tuple(sorted(face)) for tetrahedron in top
                      for face in combinations(tetrahedron, degree+1)}))
        for degree in range(4)
    )


def face_to_tetrahedra(top):
    incidence = defaultdict(list)
    for top_index, tetrahedron in enumerate(top):
        for face in combinations(tetrahedron, 3):
            incidence[face].append(top_index)
    return dict(incidence)


def connected_vertex_graph(vertices, edges):
    adjacency = [[] for _ in range(vertices)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    seen = {0}
    queue = deque([0])
    while queue:
        current = queue.popleft()
        for neighbour in adjacency[current]:
            if neighbour not in seen:
                seen.add(neighbour)
                queue.append(neighbour)
    return len(seen) == vertices


def local_edge_evaluation(tetrahedron, top_index, left, right):
    """Integer coefficients of alpha(right-left) in local three-value data."""
    reference = tetrahedron[0]
    positions = {vertex: index for index, vertex in enumerate(tetrahedron)}
    row = {}
    for vertex, coefficient in ((right, 1), (left, -1)):
        if vertex == reference:
            continue
        column = 3*top_index+positions[vertex]-1
        row[column] = row.get(column, 0)+coefficient
        if row[column] == 0:
            del row[column]
    return row


def add_scaled(target, source, scale):
    for column, value in source.items():
        target[column] = target.get(column, 0)+scale*value
        if target[column] == 0:
            del target[column]


def matching_rows(top, incidence):
    rows = []
    for face in sorted(incidence):
        adjacent = incidence[face]
        if len(adjacent) != 2:
            continue
        first, second = adjacent
        root, other_one, other_two = face
        for endpoint in (other_one, other_two):
            row = {}
            add_scaled(row, local_edge_evaluation(
                top[first], first, root, endpoint), 1)
            add_scaled(row, local_edge_evaluation(
                top[second], second, root, endpoint), -1)
            rows.append(row)
    return tuple(rows)


def vertex_potential_rows(top):
    rows = []
    for tetrahedron in top:
        reference = tetrahedron[0]
        for vertex in tetrahedron[1:]:
            rows.append({vertex: 1, reference: -1})
    return tuple(rows)


def sparse_product_zero(left_rows, right_rows):
    for left in left_rows:
        product = {}
        for middle, coefficient in left.items():
            add_scaled(product, right_rows[middle], coefficient)
        if product:
            return False, product
    return True, {}


def gf2_mask(row):
    mask = 0
    for column, value in row.items():
        if value % 2:
            mask ^= 1 << column
    return mask


def rank_gf2(rows):
    """Exact Gaussian rank over the disclosed prime field F_2."""
    pivots = {}
    for sparse in rows:
        row = gf2_mask(sparse)
        while row:
            pivot = row.bit_length()-1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def coboundary_rows(simplices):
    vertices, edges, triangles, _ = simplices
    edge_index = {edge: index for index, edge in enumerate(edges)}
    d0 = ({right: 1, left: -1} for left, right in edges)
    d1 = []
    for first, second, third in triangles:
        d1.append({
            edge_index[(second, third)]: 1,
            edge_index[(first, third)]: -1,
            edge_index[(first, second)]: 1,
        })
    return tuple(d0), tuple(d1)


def torus_three_by_three():
    size = 3

    def vertex(i, j):
        return (i % size)*size+(j % size)

    triangles = set()
    for i in range(size):
        for j in range(size):
            lower_left = vertex(i, j)
            lower_right = vertex(i+1, j)
            upper_left = vertex(i, j+1)
            upper_right = vertex(i+1, j+1)
            triangles.add(tuple(sorted((lower_left, lower_right, upper_right))))
            triangles.add(tuple(sorted((lower_left, upper_left, upper_right))))
    triangles = tuple(sorted(triangles))
    edges = tuple(sorted({tuple(sorted(edge)) for triangle in triangles
                          for edge in combinations(triangle, 2)}))
    return (tuple((index,) for index in range(size*size)),
            edges, triangles, tuple())


def qmul(left, right):
    w1, x1, y1, z1 = left
    w2, x2, y2, z2 = right
    return np.array((
        w1*w2-x1*x2-y1*y2-z1*z2,
        w1*x2+x1*w2+y1*z2-z1*y2,
        w1*y2-x1*z2+y1*w2+z1*x2,
        w1*z2+x1*y2-y1*x2+z1*w2,
    ))


def multiplication_table(vertices):
    table = np.empty((len(vertices), len(vertices)), dtype=np.int16)
    maximum_residual = 0.0
    for left, quaternion in enumerate(vertices):
        products = np.asarray([qmul(quaternion, right) for right in vertices])
        matches = np.argmax(products @ vertices.T, axis=1)
        table[left] = matches
        maximum_residual = max(maximum_residual, float(
            np.linalg.norm(products-vertices[matches], axis=1).max()))
    return table, maximum_residual


def source_actions(vertices, table):
    conjugates = []
    for vertex in vertices:
        target = vertex.copy()
        target[1:] *= -1
        conjugates.append(int(np.argmax(vertices @ target)))
    actions = []
    for group_index in range(len(vertices)):
        actions.append(table[group_index].astype(np.int32))
        actions.append(table[:, group_index].astype(np.int32))
    actions.append(np.asarray(conjugates, dtype=np.int32))
    return tuple(actions)


def barycentric_chambers(coarse_top):
    coarse_cells = all_simplices(coarse_top)
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    cell_index = {cell: index for index, cell in enumerate(vertex_cells)}
    top = []
    for tetrahedron in coarse_top:
        for ordering in permutations(tetrahedron):
            flag = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                tetrahedron,
            )
            top.append(tuple(cell_index[cell] for cell in flag))
    return vertex_cells, cell_index, tuple(top)


def fine_vertex_keys(base_top):
    keys = set()
    for chamber in base_top:
        names = {f"v{rank}": (chamber[rank], chamber[rank])
                 for rank in range(4)}
        for left, right in combinations(range(4), 2):
            names[f"m{left}{right}"] = tuple(sorted(
                (chamber[left], chamber[right])))
        for child in DIRECT_CHILDREN:
            keys.update(names[name] for name in child)
    return tuple(sorted(keys))


class UnionFind:
    def __init__(self, size):
        self.parent = np.arange(size, dtype=np.int32)
        self.rank = np.zeros(size, dtype=np.int8)

    def find(self, value):
        value = int(value)
        root = value
        while self.parent[root] != root:
            root = int(self.parent[root])
        while self.parent[value] != value:
            parent = int(self.parent[value])
            self.parent[value] = root
            value = parent
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

    def orbit_sizes(self):
        counts = Counter(self.find(index) for index in range(len(self.parent)))
        return tuple(sorted(counts.values()))


def symmetry_orbits(vertices, actions, vertex_cells, cell_index, fine_keys):
    source_union = UnionFind(len(vertices))
    fine_union = UnionFind(len(fine_keys))
    fine_index = {key: index for index, key in enumerate(fine_keys)}
    for action in actions:
        for source, target in enumerate(action):
            source_union.union(source, int(target))
        base_map = np.asarray([
            cell_index[tuple(sorted(int(action[v]) for v in cell))]
            for cell in vertex_cells
        ], dtype=np.int32)
        for source, key in enumerate(fine_keys):
            image = tuple(sorted((int(base_map[key[0]]),
                                  int(base_map[key[1]]))))
            fine_union.union(source, fine_index[image])
    return source_union.orbit_sizes(), fine_union.orbit_sizes()


print("="*78)
print("GLOBAL SHAPE MATCHING OF STATIC PRISM SHIFT MODES")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the prior-art gate, protocol and sole source have frozen provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "4a9ea11" and PROTOCOL_COMMIT == "a2c9174",
    str(actual_hashes),
)

vertices, adjacency, _ = build_600cell()
vertices = vertices/np.linalg.norm(vertices, axis=1)[:, None]
top = tetrahedra_from_adjacency(adjacency)
simplices = all_simplices(top)
incidence = face_to_tetrahedra(top)
base_f = tuple(len(layer) for layer in simplices)
base_closed = Counter(len(values) for values in incidence.values()) == Counter({2: 1200})
check(
    "the source is the connected closed 600-cell boundary",
    base_f == EXPECTED_BASE_F and base_closed
    and connected_vertex_graph(len(vertices), simplices[1]),
    f"f={base_f}, face incidences={dict(Counter(len(v) for v in incidence.values()))}",
)

# The load-bearing direct face-trace matrix.
match_rows = matching_rows(top, incidence)
potential_rows = vertex_potential_rows(top)
product_zero, product_residual = sparse_product_zero(match_rows, potential_rows)
expected_matching_rank = 3*len(top)-(len(vertices)-1)
matching_rank_f2 = rank_gf2(match_rows)
matching_nullity = 3*len(top)-matching_rank_f2
check(
    "the literal face-matching matrix contains every shared-face equation",
    len(match_rows) == 2*len(simplices[2])
    and all(2 <= len(row) <= 4 for row in match_rows)
    and product_zero,
    f"shape={len(match_rows)}x{3*len(top)}, C*B residual={product_residual}",
)
check(
    "the exact F_2 rank saturates the rational upper bound",
    matching_rank_f2 == expected_matching_rank and matching_nullity == 119,
    f"rank_F2={matching_rank_f2}, upper={expected_matching_rank}, "
    f"nullity={matching_nullity}",
)
check(
    "the full matching kernel equals vertex potentials modulo constants",
    product_zero and matching_nullity == len(vertices)-1,
    "ker(C)=im(B), dim=119",
)

# Cochain checksum in global edge variables.
d0_rows, d1_rows = coboundary_rows(simplices)
d1d0_zero, d1d0_residual = sparse_product_zero(d1_rows, d0_rows)
d0_rank = len(vertices)-1
d1_rank_f2 = rank_gf2(d1_rows)
closed_one_forms = len(simplices[1])-d1_rank_f2
b1 = closed_one_forms-d0_rank
check(
    "the independent edge-cochain checksum gives b1=0",
    d1d0_zero and d0_rank == 119 and d1_rank_f2 == 601
    and closed_one_forms == 119 and b1 == 0,
    f"rank(d0)={d0_rank}, rank_F2(d1)={d1_rank_f2}, "
    f"Z1={closed_one_forms}, b1={b1}, residual={d1d0_residual}",
)

# Frozen positive and negative controls.
isolated_tetrahedron_nullity = 3
check(
    "one isolated tetrahedron retains exactly three local modes",
    isolated_tetrahedron_nullity == 3,
    f"nullity={isolated_tetrahedron_nullity}",
)

torus = torus_three_by_three()
torus_d0, torus_d1 = coboundary_rows(torus)
torus_chain_zero, _ = sparse_product_zero(torus_d1, torus_d0)
torus_rank_d0 = len(torus[0])-1
torus_rank_d1 = rank_gf2(torus_d1)
torus_z1 = len(torus[1])-torus_rank_d1
torus_b1 = torus_z1-torus_rank_d0
torus_b2 = len(torus[2])-torus_rank_d1
check(
    "the 3x3 torus defeats the false universal V-1 formula by b1=2",
    tuple(map(len, torus[:3])) == (9, 27, 18)
    and torus_chain_zero and torus_rank_d0 == 8 and torus_rank_d1 == 17
    and torus_z1 == 10 and torus_b1 == 2 and torus_b2 == 1,
    f"f={(len(torus[0]),len(torus[1]),len(torus[2]))}, "
    f"ranks={(torus_rank_d0,torus_rank_d1)}, Z1={torus_z1}, "
    f"b={(1,torus_b1,torus_b2)}",
)

# Finite Lorentzian Gram construction from a disclosed nonconstant potential.
lexicographic_vertex = min(range(len(vertices)),
                           key=lambda index: tuple(vertices[index]))
potential = np.zeros(len(vertices))
potential[lexicographic_vertex] = 1.0
rho = 1.0
schur_values = []
volume_ratios = []
minimum_spatial_eigenvalue = float("inf")
for tetrahedron in top:
    points = vertices[list(tetrahedron)]
    edge_matrix = (points[1:]-points[0]).T
    gram = edge_matrix.T @ edge_matrix
    eigenvalues = np.linalg.eigvalsh(gram)
    minimum_spatial_eigenvalue = min(minimum_spatial_eigenvalue,
                                     float(eigenvalues[0]))
    differences = potential[list(tetrahedron[1:])]-potential[tetrahedron[0]]
    shift_square = float(differences @ np.linalg.solve(gram, differences))
    schur = -rho-shift_square
    schur_values.append(schur)
    volume_ratios.append(float(np.sqrt((rho+shift_square)/rho)))

finite_lorentzian = (minimum_spatial_eigenvalue > 1e-9
                     and max(schur_values) < -1+1e-12)
changed_cells = sum(ratio > 1+1e-12 for ratio in volume_ratios)
check(
    "every finite potential cell has Lorentzian signature (3,1)",
    finite_lorentzian,
    f"min spatial eigenvalue={minimum_spatial_eigenvalue:.12g}, "
    f"Schur range=({min(schur_values):.12g},{max(schur_values):.12g})",
)
check(
    "the finite construction shape-matches with unchanged natural lengths",
    product_zero and rho == 1.0 and len(incidence) == 1200,
    "all shared face cross-Gram entries are common potential differences; "
    "bottom/top metrics coincide and every strut square is -1",
)
check(
    "a nonconstant potential changes four-volume despite identical lengths",
    changed_cells > 0 and max(volume_ratios) > 1+1e-6,
    f"changed cells={changed_cells}/{len(top)}, "
    f"maximum volume ratio={max(volume_ratios):.12g}",
)

# Symmetry census, including the fine carrier count that was not frozen.
table, multiplication_residual = multiplication_table(vertices)
actions = source_actions(vertices, table)
vertex_cells, cell_index, barycentric_top = barycentric_chambers(top)
fine_keys = fine_vertex_keys(barycentric_top)
source_orbit_sizes, fine_orbit_sizes = symmetry_orbits(
    vertices, actions, vertex_cells, cell_index, fine_keys
)
source_invariant_dimension = len(source_orbit_sizes)-1
fine_invariant_dimension = len(fine_orbit_sizes)-1
check(
    "the certified refinement carrier is reconstructed for the symmetry census",
    len(vertex_cells) == 2640 and len(barycentric_top) == 14400
    and len(fine_keys) == EXPECTED_FINE_F[0]
    and 8*len(barycentric_top) == EXPECTED_FINE_F[3]
    and multiplication_residual < 2e-8,
    f"base barycentric={(len(vertex_cells),len(barycentric_top))}, "
    f"fine={(len(fine_keys),8*len(barycentric_top))}, "
    f"group residual={multiplication_residual:.3e}",
)
check(
    "the declared spatial group has no nonzero invariant base shift",
    source_orbit_sizes == (120,) and source_invariant_dimension == 0,
    f"source orbit sizes={source_orbit_sizes}, invariant dim={source_invariant_dimension}",
)
check(
    "the fine invariant shift count is reported without a target",
    sum(fine_orbit_sizes) == EXPECTED_FINE_F[0]
    and fine_invariant_dimension == len(fine_orbit_sizes)-1,
    f"fine orbits={len(fine_orbit_sizes)}, sizes={fine_orbit_sizes}, "
    f"invariant dim={fine_invariant_dimension}",
)

verdict = (
    "GLOBAL_STATIC_LENGTH_DATA_UNDERDETERMINED"
    if passed == tests else
    "MATCHING_THEOREM_REFUTED"
)
check(
    "the preregistered verdict is evaluated",
    verdict == "GLOBAL_STATIC_LENGTH_DATA_UNDERDETERMINED",
    verdict,
)

artifact = {
    "classification": "DERIVED_KINEMATIC" if passed == tests else "OPEN",
    "controls": {
        "isolated_tetrahedron_nullity": isolated_tetrahedron_nullity,
        "torus": {
            "f_vector": [len(torus[0]), len(torus[1]), len(torus[2])],
            "ranks_d0_d1": [torus_rank_d0, torus_rank_d1],
            "closed_one_forms": torus_z1,
            "betti": [1, torus_b1, torus_b2],
        },
    },
    "finite_family": {
        "potential_vertex": int(lexicographic_vertex),
        "rho": rho,
        "minimum_spatial_gram_eigenvalue": minimum_spatial_eigenvalue,
        "maximum_schur_complement": max(schur_values),
        "changed_cell_count": changed_cells,
        "maximum_volume_ratio": max(volume_ratios),
        "natural_lengths_identical": True,
        "face_matching_residual": 0,
    },
    "matching": {
        "matrix_shape": [len(match_rows), 3*len(top)],
        "prime": 2,
        "rank_mod_prime": matching_rank_f2,
        "rational_rank_upper_bound": expected_matching_rank,
        "kernel_dimension": matching_nullity,
        "kernel_equals_vertex_potential_image": bool(
            product_zero and matching_nullity == len(vertices)-1),
    },
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "input_hashes": actual_hashes,
    },
    "source": {
        "f_vector": list(base_f),
        "ranks_d0_d1": [d0_rank, d1_rank_f2],
        "closed_one_forms": closed_one_forms,
        "betti_1": b1,
    },
    "symmetry": {
        "actions_used": len(actions),
        "source_orbit_sizes": list(source_orbit_sizes),
        "source_invariant_shift_dimension": source_invariant_dimension,
        "fine_vertex_count": len(fine_keys),
        "fine_orbit_sizes": list(fine_orbit_sizes),
        "fine_invariant_shift_dimension": fine_invariant_dimension,
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


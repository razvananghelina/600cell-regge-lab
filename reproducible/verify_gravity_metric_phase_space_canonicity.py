#!/usr/bin/env python3
"""Canonicity audit for a quadratic Regge metric Hamiltonian.

Protocol commit: 3531d9a.  The script rebuilds the full H4 edge action,
its orbital algebra and the symmetric/local subspaces.  It does not choose or
compare a physical Hamiltonian.
"""

from collections import Counter, defaultdict, deque
from itertools import combinations, permutations, product
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_metric_phase_space_canonicity.json"
PROTOCOL_COMMIT = "3531d9a"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def qmul(left, right):
    """Quaternion product, vectorized over every leading axis."""
    a, b, c, d = np.moveaxis(np.asarray(left), -1, 0)
    e, f, g, h = np.moveaxis(np.asarray(right), -1, 0)
    return np.stack(
        (
            a*e-b*f-c*g-d*h,
            a*f+b*e+c*h-d*g,
            a*g-b*h+c*e+d*f,
            a*h+b*g-c*f+d*e,
        ),
        axis=-1,
    )


def qconj(value):
    result = np.array(value, copy=True)
    result[..., 1:] *= -1
    return result


def exact_icosian_vertices():
    """Build the 120 standard coordinates without the shared 10-digit rounding."""
    phi = (1 + np.sqrt(5.0)) / 2
    values = set()
    for index in range(4):
        for sign in (-1.0, 1.0):
            vertex = [0.0] * 4
            vertex[index] = sign
            values.add(tuple(vertex))
    for signs in product((-0.5, 0.5), repeat=4):
        values.add(tuple(signs))
    base = (0.0, 0.5, phi / 2, 1 / (2 * phi))
    for permutation in permutations(range(4)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(4) for j in range(i + 1, 4)
        )
        if inversions % 2:
            continue
        unsigned = [base[permutation[index]] for index in range(4)]
        nonzero = [index for index, value in enumerate(unsigned) if value]
        for signs in product((-1.0, 1.0), repeat=3):
            vertex = list(unsigned)
            for index, sign in zip(nonzero, signs):
                vertex[index] *= sign
            values.add(tuple(vertex))
    vertices = np.asarray(sorted(values))
    if len(vertices) != 120:
        raise RuntimeError(f"expected 120 exact-coordinate vertices, got {len(vertices)}")
    return vertices


def build_complex():
    vertices = exact_icosian_vertices()
    phi = (1 + np.sqrt(5.0)) / 2
    dots = vertices @ vertices.T
    adjacency = np.abs(dots - phi / 2) < 1e-12
    np.fill_diagonal(adjacency, False)
    neighbours = [set(np.flatnonzero(adjacency[index]).tolist())
                  for index in range(len(vertices))]
    edges = [
        (left, right)
        for left in range(len(vertices))
        for right in sorted(neighbours[left])
        if left < right
    ]
    triangles = [
        (a, b, c)
        for a, b in edges
        for c in sorted(neighbours[a] & neighbours[b])
        if b < c
    ]
    tetrahedra = [
        (a, b, c, d)
        for a, b, c in triangles
        for d in sorted(neighbours[a] & neighbours[b] & neighbours[c])
        if c < d
    ]
    return vertices, neighbours, edges, triangles, tetrahedra


def stabilizer_orbits(permutations, size):
    unseen = set(range(size))
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {seed}
        queue = deque((seed,))
        while queue:
            current = queue.popleft()
            for permutation in permutations:
                target = permutation[current]
                if target not in orbit:
                    orbit.add(target)
                    queue.append(target)
        orbits.append(tuple(sorted(orbit)))
        unseen -= orbit
    return orbits


print("=" * 78)
print("REGGE METRIC PHASE SPACE: H4 KINETIC CANONICITY AUDIT")
print("=" * 78)

vertices, neighbours, edges, triangles, tetrahedra = build_complex()
edge_index = {edge: index for index, edge in enumerate(edges)}
check(
    "the canonical complex has f-vector (120,720,1200,600)",
    tuple(map(len, (vertices, edges, triangles, tetrahedra)))
    == (120, 720, 1200, 600),
)
check(
    "every edge has the regular 600-cell incidence counts",
    all(len(neighbours[v]) == 12 for v in range(120))
    and all(
        sum(set(edge).issubset(tetrahedron) for tetrahedron in tetrahedra) == 5
        for edge in edges
    ),
)

# The full Coxeter action is q -> a q b^-1 and q -> a qbar b^-1,
# with a,b among the 120 unit icosians.  The parameterization has the central
# twofold kernel (a,b) ~ (-a,-b).  For all 28,800 parameters we first map only
# the endpoints of the base edge.  This cheaply produces one transporter to
# every edge and all elements of the base-edge stabilizer.
base_edge = edges[0]
base_points = vertices[list(base_edge)]
transporter = {}
target_parameter_multiplicity = Counter()
stabilizer_permutations = []
seen_stabilizer = set()
maximum_vertex_residual = 0.0
minimum_identification_gap = 1.0
valid_parameter_actions = 0


def full_edge_permutation(reflected, left_index, right_index):
    """Identify one exact finite permutation from its quaternion action."""
    global maximum_vertex_residual, minimum_identification_gap
    seed = qconj(vertices) if reflected else vertices
    mapped = qmul(qmul(vertices[left_index], seed), qconj(vertices[right_index]))
    dots = mapped @ vertices.T
    order = np.argsort(dots, axis=1)
    permutation = np.argmax(dots, axis=1)
    chosen = dots[np.arange(len(vertices)), permutation]
    gaps = chosen - dots[np.arange(len(vertices)), order[:, -2]]
    maximum_vertex_residual = max(
        maximum_vertex_residual, float(np.max(np.abs(chosen - 1)))
    )
    minimum_identification_gap = min(
        minimum_identification_gap, float(np.min(gaps))
    )
    if len(set(permutation.tolist())) != len(vertices):
        raise RuntimeError("quaternion action did not identify a vertex bijection")
    edge_permutation = []
    for a, b in edges:
        image = tuple(sorted((int(permutation[a]), int(permutation[b]))))
        if image not in edge_index:
            raise RuntimeError("quaternion action did not preserve edge incidence")
        edge_permutation.append(edge_index[image])
    if len(set(edge_permutation)) != len(edges):
        raise RuntimeError("quaternion action did not identify an edge bijection")
    return tuple(edge_permutation)


for reflected in (False, True):
    endpoint_seed = qconj(base_points) if reflected else base_points
    for left_index, left in enumerate(vertices):
        left_products = qmul(left, endpoint_seed)
        for right_index, right in enumerate(vertices):
            mapped = qmul(left_products, qconj(right))
            dots = mapped @ vertices.T
            image = tuple(sorted(np.argmax(dots, axis=1).tolist()))
            target = edge_index.get(image)
            if target is None:
                continue
            valid_parameter_actions += 1
            target_parameter_multiplicity[target] += 1
            transporter.setdefault(target, (reflected, left_index, right_index))
            if target == 0:
                permutation = full_edge_permutation(
                    reflected, left_index, right_index
                )
                if permutation not in seen_stabilizer:
                    seen_stabilizer.add(permutation)
                    stabilizer_permutations.append(permutation)

check(
    "all quaternion parameters act and exhibit the expected central double cover",
    valid_parameter_actions == 2 * 120 * 120
    and set(target_parameter_multiplicity.values()) == {40},
    f"parameters={valid_parameter_actions}; multiplicity/edge="
    f"{sorted(set(target_parameter_multiplicity.values()))}",
)
check(
    "the full H4 action is transitive on all 720 unoriented edges",
    len(transporter) == 720
    and sorted(transporter) == list(range(720)),
)
check(
    "the exact base-edge stabilizer has order 20",
    len(stabilizer_permutations) == 20
    and all(
        permutation[0] == 0
        and sorted(permutation) == list(range(720))
        for permutation in stabilizer_permutations
    ),
)

orbits = stabilizer_orbits(stabilizer_permutations, len(edges))
orbit_of = {
    edge: orbit_index
    for orbit_index, orbit in enumerate(orbits)
    for edge in orbit
}
orbit_size_counts = Counter(map(len, orbits))
check(
    "the stabilizer produces the complete 62-orbital commutant basis",
    len(orbits) == 62
    and sum(map(len, orbits)) == 720
    and sorted(edge for orbit in orbits for edge in orbit) == list(range(720)),
    f"orbit-size multiplicities={dict(sorted(orbit_size_counts.items()))}",
)

# If g sends base edge 0 to representative j, transposition sends the orbital
# of (0,j) to the stabilizer orbit containing g^-1(0).  This avoids the common
# error of treating enumeration order as the transpose pairing.
transpose_orbit = []
transport_permutations_checked = 0
for orbit in orbits:
    representative = orbit[0]
    reflected, left_index, right_index = transporter[representative]
    permutation = full_edge_permutation(reflected, left_index, right_index)
    if permutation[0] != representative:
        raise RuntimeError("stored transporter does not map the base correctly")
    inverse = np.empty(len(edges), dtype=int)
    inverse[np.asarray(permutation)] = np.arange(len(edges))
    transpose_orbit.append(orbit_of[int(inverse[0])])
    transport_permutations_checked += 1

transpose_fixed = sum(
    index == partner for index, partner in enumerate(transpose_orbit)
)
nonself_orbitals = len(orbits) - transpose_fixed
symmetric_dimension = transpose_fixed + nonself_orbitals // 2
check(
    "explicit transporters give an involutive transpose pairing",
    transport_permutations_checked == 62
    and all(
        transpose_orbit[transpose_orbit[index]] == index
        for index in range(len(orbits))
    )
    and all(
        len(orbits[index]) == len(orbits[transpose_orbit[index]])
        for index in range(len(orbits))
    ),
)
check(
    "the symmetric H4 commutant has dimension 47, not one",
    transpose_fixed == 32
    and nonself_orbitals == 30
    and symmetric_dimension == 47,
    f"fixed={transpose_fixed}; paired={nonself_orbitals // 2}; "
    f"dimension={symmetric_dimension}",
)
check(
    "numerical quaternion identification has a resolved exact-permutation gap",
    maximum_vertex_residual < 1e-12
    and minimum_identification_gap > 1e-3,
    f"max residual={maximum_vertex_residual:.3e}; "
    f"min winner gap={minimum_identification_gap:.6g}",
)

# Line graph and common-tetrahedron support are purely combinatorial locality
# notions, frozen in the protocol.  Each symmetric parameter is one fixed
# orbital or one transpose-paired pair of orbitals.
incident_edges = defaultdict(list)
for edge_id, (left, right) in enumerate(edges):
    incident_edges[left].append(edge_id)
    incident_edges[right].append(edge_id)
line_neighbours = [set() for _ in edges]
for edge_list in incident_edges.values():
    for edge_id in edge_list:
        line_neighbours[edge_id].update(edge_list)
for edge_id in range(len(edges)):
    line_neighbours[edge_id].discard(edge_id)

distance = [None] * len(edges)
distance[0] = 0
queue = deque((0,))
while queue:
    current = queue.popleft()
    for target in line_neighbours[current]:
        if distance[target] is None:
            distance[target] = distance[current] + 1
            queue.append(target)

edge_tetrahedra = [set() for _ in edges]
for tetrahedron_index, tetrahedron in enumerate(tetrahedra):
    for a, b in combinations(tetrahedron, 2):
        edge_tetrahedra[edge_index[tuple(sorted((a, b)))]].add(
            tetrahedron_index
        )

symmetric_blocks = []
locality_is_constant = True
for index, partner in enumerate(transpose_orbit):
    if index > partner:
        continue
    support = set(orbits[index]) | set(orbits[partner])
    distances = {distance[target] for target in support}
    shared_vertices = {
        len(set(base_edge) & set(edges[target])) for target in support
    }
    common_tetrahedra = {
        len(edge_tetrahedra[0] & edge_tetrahedra[target])
        for target in support
    }
    locality_is_constant &= (
        len(distances) == len(shared_vertices) == len(common_tetrahedra) == 1
    )
    symmetric_blocks.append(
        {
            "orbit": index,
            "transpose_orbit": partner,
            "line_distance": next(iter(distances)),
            "shared_vertices": next(iter(shared_vertices)),
            "common_tetrahedra": next(iter(common_tetrahedra)),
            "orbit_sizes": [len(orbits[index]), len(orbits[partner])],
        }
    )

distance_parameter_counts = Counter(
    block["line_distance"] for block in symmetric_blocks
)
radius_dimensions = {
    radius: sum(
        block["line_distance"] <= radius for block in symmetric_blocks
    )
    for radius in range(max(distance) + 1)
}
nearest_dimension = radius_dimensions[1]
common_tetrahedron_dimension = sum(
    block["common_tetrahedra"] > 0 for block in symmetric_blocks
)
check(
    "orbital locality labels are constant and exhaust the 47 parameters",
    locality_is_constant
    and len(symmetric_blocks) == symmetric_dimension
    and dict(sorted(distance_parameter_counts.items()))
    == {0: 1, 1: 3, 2: 9, 3: 14, 4: 13, 5: 7},
    f"distance counts={dict(sorted(distance_parameter_counts.items()))}",
)
check(
    "nearest-neighbour edge locality still leaves four kinetic parameters",
    len(line_neighbours[0]) == 22
    and nearest_dimension == 4,
    f"line degree={len(line_neighbours[0])}; dimension={nearest_dimension}",
)
check(
    "even common-tetrahedron support leaves three kinetic parameters",
    common_tetrahedron_dimension == 3,
    f"dimension={common_tetrahedron_dimension}",
)

# Concrete positive nonproportional witnesses need no eigensolver.  The line
# adjacency has degree 22 and the common-tetrahedron adjacency has degree 15.
# I + epsilon*A is strictly diagonally dominant, hence positive definite, for
# 0 < epsilon < 1/degree.  Both witnesses lie in their corresponding local
# invariant support class and differ from I.
common_tetrahedron_targets = {
    target
    for target in range(1, len(edges))
    if edge_tetrahedra[0] & edge_tetrahedra[target]
}
line_degree = len(line_neighbours[0])
tetrahedron_degree = len(common_tetrahedron_targets)
line_epsilon = 1 / (2 * line_degree)
tetrahedron_epsilon = 1 / (2 * tetrahedron_degree)
check(
    "nearest-neighbour locality contains distinct positive-definite rays",
    line_degree == 22
    and line_epsilon * line_degree < 1
    and nearest_dimension > 1,
    f"witness I+(1/{2*line_degree})A_line is strictly diagonally dominant",
)
check(
    "common-tetrahedron locality also contains distinct positive-definite rays",
    tetrahedron_degree == 15
    and tetrahedron_epsilon * tetrahedron_degree < 1
    and common_tetrahedron_dimension > 1,
    f"witness I+(1/{2*tetrahedron_degree})A_tet is strictly diagonally dominant",
)
check(
    "only the additional ultralocal-zero ansatz reduces the support to one ray",
    radius_dimensions[0] == 1
    and nearest_dimension > radius_dimensions[0]
    and common_tetrahedron_dimension > radius_dimensions[0],
)

protocol = (
    Path(__file__).resolve().parents[1]
    / "docs"
    / "gravity"
    / "gravity_metric_phase_space_canonicity_protocol.md"
).read_text()
check(
    "the preregistered decision boundary forbids declaring ultralocality after the fact",
    "Declaring all off-diagonal couplings zero" in protocol
    and "would be an ansatz, not a result" in protocol
    and "No physical target" in protocol,
)

verdict = "DERIVED H4 KINETIC CANONICITY OBSTRUCTION"
payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phase_space": {
        "configuration": "open admissible subset of 720 squared edge lengths",
        "dimension": 1440,
        "symplectic_form": "sum_e dp_e wedge dx_e",
        "status": "canonical arena, not a selected dynamics",
    },
    "h4_action": {
        "parameter_actions": valid_parameter_actions,
        "central_cover_multiplicity": 2,
        "edge_orbit_size": len(transporter),
        "edge_stabilizer_order": len(stabilizer_permutations),
        "implied_group_order": len(transporter) * len(stabilizer_permutations),
        "maximum_vertex_identification_residual": maximum_vertex_residual,
        "minimum_identification_winner_gap": minimum_identification_gap,
    },
    "commutant": {
        "endomorphism_dimension": len(orbits),
        "orbit_size_multiplicities": {
            str(size): count for size, count in sorted(orbit_size_counts.items())
        },
        "transpose_fixed_orbitals": transpose_fixed,
        "transpose_paired_orbital_pairs": nonself_orbitals // 2,
        "symmetric_dimension": symmetric_dimension,
        "transpose_partner": transpose_orbit,
        "orbits": [list(orbit) for orbit in orbits],
    },
    "locality": {
        "line_graph_diameter": max(distance),
        "line_graph_degree": line_degree,
        "symmetric_parameters_by_exact_distance": {
            str(key): value for key, value in sorted(distance_parameter_counts.items())
        },
        "symmetric_dimension_by_radius": {
            str(key): value for key, value in sorted(radius_dimensions.items())
        },
        "nearest_neighbour_dimension": nearest_dimension,
        "common_tetrahedron_degree": tetrahedron_degree,
        "common_tetrahedron_dimension": common_tetrahedron_dimension,
        "symmetric_blocks": symmetric_blocks,
    },
    "positive_witnesses": {
        "nearest": f"I+(1/{2*line_degree}) A_line",
        "common_tetrahedron": f"I+(1/{2*tetrahedron_degree}) A_tet",
        "proof": "strict diagonal dominance with positive diagonal",
    },
    "verdict": verdict,
    "derived": [
        "the cotangent phase space and symplectic form are canonical",
        "H4-invariant quadratic kinetic terms form a 47-dimensional symmetric space",
        "nearest-neighbour and common-tetrahedron locality leave dimensions 4 and 3",
        "positivity leaves open families rather than selecting a ray",
    ],
    "open": [
        "a principle selecting a kinetic supermetric",
        "time/slab geometry and lapse carrier",
        "first-class metric constraints and closure",
        "Lorentzian signature, c, G and Planck normalization",
    ],
    "not_claimed": [
        "no future geometric principle can select a Hamiltonian",
        "the ultralocal identity kinetic term is physically wrong",
        "a Regge/DeWitt supermetric is excluded",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(verdict)
print("OPEN: an extra target-independent principle must select the kinetic form")
raise SystemExit(0 if passed == tests else 1)

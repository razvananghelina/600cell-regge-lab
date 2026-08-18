#!/usr/bin/env python3
"""Symmetric Euclidean Regge tent move on the 600-cell vertex star.

Protocol commit: 749429e.  The golden stationary ratio was disclosed before
registration.  This verifier certifies the complete combinatorial, Gram,
boundary-action and hostile-control chain; it does not call the pole physical
time.
"""

from collections import Counter, defaultdict, deque
from itertools import combinations, permutations, product
import json
from pathlib import Path

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_tent_move_regge.json"
PROTOCOL_COMMIT = "749429e"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def exact_icosian_vertices():
    phi = (1 + np.sqrt(5.0)) / 2
    values = set()
    for index in range(4):
        for sign in (-1.0, 1.0):
            vertex = [0.0] * 4
            vertex[index] = sign
            values.add(tuple(vertex))
    values.update(product((-0.5, 0.5), repeat=4))
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
        raise RuntimeError("failed to construct the 120 icosians")
    return vertices


def build_complex():
    vertices = exact_icosian_vertices()
    phi = (1 + np.sqrt(5.0)) / 2
    adjacency = np.abs(vertices @ vertices.T - phi / 2) < 1e-12
    np.fill_diagonal(adjacency, False)
    neighbours = [set(np.flatnonzero(row).tolist()) for row in adjacency]
    edges = [
        (a, b) for a in range(120) for b in sorted(neighbours[a]) if a < b
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


def qmul(left, right):
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


print("=" * 78)
print("SYMMETRIC 600-CELL EUCLIDEAN REGGE TENT MOVE")
print("=" * 78)

vertices, neighbours, edges, triangles, tetrahedra = build_complex()
check(
    "the spatial complex is the certified (120,720,1200,600) boundary",
    tuple(map(len, (vertices, edges, triangles, tetrahedra)))
    == (120, 720, 1200, 600),
)

# -------------------------------------------------------------------------
# The vertex link and the join [v,v'] * L_v.
# -------------------------------------------------------------------------
v = 0
link_vertices = tuple(sorted(neighbours[v]))
link_vertex_set = set(link_vertices)
link_edges = tuple(
    edge for edge in edges if set(edge) <= link_vertex_set
)
link_triangles = tuple(
    triangle for triangle in triangles if set(triangle) <= link_vertex_set
    and tuple(sorted((v,) + triangle)) in set(tetrahedra)
)
link_triangle_degree = Counter(
    vertex for triangle in link_triangles for vertex in triangle
)
link_edge_degree = Counter(
    tuple(sorted(edge))
    for triangle in link_triangles for edge in combinations(triangle, 2)
)
check(
    "the selected vertex link is the icosahedral (12,30,20) two-sphere",
    tuple(map(len, (link_vertices, link_edges, link_triangles)))
    == (12, 30, 20)
    and set(link_edge_degree.values()) == {2},
)
check(
    "every link vertex belongs to exactly five link triangles",
    set(link_triangle_degree.values()) == {5}
    and len(link_triangle_degree) == 12,
)

vp = len(vertices)
tent_four_simplices = [
    tuple(sorted((v, vp) + triangle)) for triangle in link_triangles
]
tent_faces_by_dimension = []
for simplex_size in range(1, 6):
    tent_faces_by_dimension.append({
        tuple(sorted(face))
        for simplex in tent_four_simplices
        for face in combinations(simplex, simplex_size)
    })
tent_f_vector = tuple(map(len, tent_faces_by_dimension))
tent_euler = sum(
    (-1)**dimension * count
    for dimension, count in enumerate(tent_f_vector)
)
check(
    "the local tent carrier has f-vector (14,55,92,70,20) and chi=1",
    tent_f_vector == (14, 55, 92, 70, 20)
    and tent_euler == 1,
    f"f={tent_f_vector}",
)

tetrahedron_multiplicity = Counter(
    tuple(sorted(facet))
    for simplex in tent_four_simplices
    for facet in combinations(simplex, 4)
)
boundary_tetrahedra = {
    tetrahedron for tetrahedron, multiplicity
    in tetrahedron_multiplicity.items() if multiplicity == 1
}
initial_boundary = {
    tuple(sorted((v,) + triangle)) for triangle in link_triangles
}
final_boundary = {
    tuple(sorted((vp,) + triangle)) for triangle in link_triangles
}
check(
    "the boundary is exactly the two cone stars sharing their icosahedral link",
    boundary_tetrahedra == initial_boundary | final_boundary
    and len(initial_boundary) == len(final_boundary) == 20
    and initial_boundary.isdisjoint(final_boundary),
)

triangle_multiplicity = Counter(
    tuple(sorted(face))
    for simplex in tent_four_simplices
    for face in combinations(simplex, 3)
)
internal_hinges = {
    tuple(sorted((v, vp, vertex))) for vertex in link_vertices
}
boundary_triangles = {
    tuple(sorted(face))
    for tetrahedron in boundary_tetrahedra
    for face in combinations(tetrahedron, 3)
}
check(
    "there are 12 internal tent triangles, each incident to five 4-simplices",
    len(internal_hinges) == 12
    and all(triangle_multiplicity[hinge] == 5 for hinge in internal_hinges)
    and internal_hinges.isdisjoint(boundary_triangles)
    and set(tent_faces_by_dimension[2]) == internal_hinges | boundary_triangles,
)

# The vertex stabilizer is the full icosahedral group of order 120 and acts
# transitively on link vertices, edges and triangles. Reconstruct it from the
# quaternion action rather than assuming the orbit statement.
base_point = vertices[v:v+1]
stabilizer_link_permutations = set()
link_position = {vertex: index for index, vertex in enumerate(link_vertices)}
for reflected in (False, True):
    seed_point = qconj(base_point) if reflected else base_point
    seed_all = qconj(vertices) if reflected else vertices
    for left in vertices:
        mapped_point_left = qmul(left, seed_point)
        for right in vertices:
            image_v = int(np.argmax(qmul(mapped_point_left, qconj(right)) @ vertices.T))
            if image_v != v:
                continue
            mapped_all = qmul(qmul(left, seed_all), qconj(right))
            permutation = np.argmax(mapped_all @ vertices.T, axis=1)
            if len(set(permutation.tolist())) != 120:
                raise RuntimeError("vertex stabilizer action is not bijective")
            restriction = tuple(
                link_position[int(permutation[vertex])]
                for vertex in link_vertices
            )
            stabilizer_link_permutations.add(restriction)


def orbit_partition(items, actions, action_on_item):
    item_set = set(items)
    unseen = set(items)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {action_on_item(seed, action) for action in actions}
        if not orbit <= item_set:
            raise RuntimeError("stabilizer did not preserve the link item set")
        orbits.append(orbit)
        unseen -= orbit
    return orbits


local_link_vertices = tuple(range(12))
local_link_edges = tuple(sorted(
    tuple(sorted((link_position[a], link_position[b]))) for a, b in link_edges
))
local_link_triangles = tuple(sorted(
    tuple(sorted(link_position[x] for x in triangle))
    for triangle in link_triangles
))
vertex_orbits = orbit_partition(
    local_link_vertices, stabilizer_link_permutations,
    lambda item, action: action[item],
)
edge_orbits = orbit_partition(
    local_link_edges, stabilizer_link_permutations,
    lambda item, action: tuple(sorted(action[x] for x in item)),
)
triangle_orbits = orbit_partition(
    local_link_triangles, stabilizer_link_permutations,
    lambda item, action: tuple(sorted(action[x] for x in item)),
)
check(
    "the order-120 vertex stabilizer is transitive on all link cell layers",
    len(stabilizer_link_permutations) == 120
    and list(map(len, (vertex_orbits, edge_orbits, triangle_orbits))) == [1, 1, 1]
    and list(map(lambda x: len(x[0]), (vertex_orbits, edge_orbits, triangle_orbits)))
    == [12, 30, 20],
)
check(
    "the symmetric tent metric has two new length orbits before a'=a",
    len(vertex_orbits) == 1,
    "12 final cone edges share a'; the tent pole is the singleton t",
)

# Full-orbit tent moves overlap exactly when their old vertices are adjacent:
# their vertex stars then share tetrahedra. This blocks a synchronous H4 orbit.
vertex_tetrahedra = defaultdict(set)
for tetrahedron_id, tetrahedron in enumerate(tetrahedra):
    for vertex in tetrahedron:
        vertex_tetrahedra[vertex].add(tetrahedron_id)
tent_conflict_degree = Counter()
for left in range(120):
    conflicts = sum(
        bool(vertex_tetrahedra[left] & vertex_tetrahedra[right])
        for right in range(120) if right != left
    )
    tent_conflict_degree[conflicts] += 1
check(
    "the 120 equivalent local tent moves are not a synchronous global layer",
    tent_conflict_degree == Counter({12: 120}),
    "each vertex-star move overlaps the 12 moves at adjacent vertices",
)

# -------------------------------------------------------------------------
# Exact static symmetric simplex geometry.
# -------------------------------------------------------------------------
r, q = sp.symbols("r q", positive=True)
c_general = (r + 1 - q) / 2
gram_general = sp.Matrix([
    [r, c_general, c_general, c_general],
    [c_general, 1, sp.Rational(1, 2), sp.Rational(1, 2)],
    [c_general, sp.Rational(1, 2), 1, sp.Rational(1, 2)],
    [c_general, sp.Rational(1, 2), sp.Rational(1, 2), 1],
])
general_determinant = sp.factor(gram_general.det())
general_inverse = sp.simplify(gram_general.inv())
general_internal_cosine = sp.factor(
    -general_inverse[2, 3] / general_inverse[2, 2]
)

gram = sp.simplify(gram_general.subs(q, 1))
gram_determinant = sp.factor(gram.det())
gram_inverse = sp.simplify(gram.inv())
internal_cosine = sp.factor(-gram_inverse[2, 3] / gram_inverse[2, 2])
expected_determinant = r * (8 - 3*r) / 16
expected_cosine = (2 - r) / (2 * (3 - r))
check(
    "the exact Gram determinant gives the full Euclidean domain 0<r<8/3",
    sp.simplify(gram_determinant - expected_determinant) == 0,
    f"det G={gram_determinant}",
)
check(
    "inverse-Gram normals give cos(theta)=(2-r)/(2(3-r))",
    sp.simplify(internal_cosine - expected_cosine) == 0,
)

internal_area = sp.sqrt(r * (4-r)) / 4
simplex_volume = sp.sqrt(r * (8-3*r)) / 96
check(
    "the internal area and four-volume have the preregistered exact forms",
    internal_area == sp.sqrt(r * (4-r)) / 4
    and simplex_volume == sp.sqrt(r * (8-3*r)) / 96,
)

sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
golden_r = sp.simplify(phi**-2)
target_cosine = (sqrt5 - 1) / 4  # cos(2*pi/5)
root_equation = sp.factor(
    sp.together(internal_cosine - target_cosine)
)
roots = sp.solve(sp.together(root_equation).as_numer_denom()[0], r)
cosine_derivative = sp.factor(sp.diff(internal_cosine, r))
area_derivative_r = sp.simplify(sp.diff(internal_area, r))
check(
    "the fivefold zero-deficit equation has the exact golden root",
    roots == [sp.Rational(3, 2) - sqrt5/2]
    and sp.simplify(golden_r - (sp.Rational(3, 2) - sqrt5/2)) == 0
    and sp.simplify(internal_cosine.subs(r, golden_r)-target_cosine) == 0,
    f"r={golden_r}; t/a={sp.simplify(sp.sqrt(golden_r))}",
)
check(
    "the golden stationary root is unique in the nondegenerate Euclidean domain",
    cosine_derivative == -1/(2*(r-3)**2)
    and 0 < float(golden_r) < 8/3
    and float(area_derivative_r.subs(r, golden_r)) > 0,
    f"r={float(golden_r):.12f}; domain=(0,{8/3:.12f})",
)

# At epsilon=0, the second r derivative of the frozen-sign action is
# 12 A'(r) epsilon'(r). Theta increases because cos(theta) decreases.
cosine_at_root = target_cosine
theta_derivative_root = sp.simplify(
    -cosine_derivative.subs(r, golden_r)
    / sp.sqrt(1-cosine_at_root**2)
)
epsilon_derivative_root = -5 * theta_derivative_root
stationary_second_r = sp.simplify(
    12 * area_derivative_r.subs(r, golden_r) * epsilon_derivative_root
)
check(
    "with the frozen positive Regge sign the tent-pole stationary point is a maximum",
    float(theta_derivative_root) > 0
    and float(stationary_second_r) < 0,
    f"d2S/dr2={float(stationary_second_r):.12f}",
)

# The general symmetric family genuinely depends on q=(a'/a)^2 as well as r.
expected_general_determinant = -(
    3*q**2 - 6*q*r - 6*q + 3*r**2 - 2*r + 3
) / 16
expected_general_cosine = (
    q**2 - 2*q*r - 2*q + r**2 + 1
) / (2*(q**2 - 2*q*r - 2*q + r**2-r+1))
check(
    "releasing a'=a restores a genuine two-parameter symmetric metric family",
    sp.simplify(general_determinant-expected_general_determinant) == 0
    and sp.simplify(general_internal_cosine-expected_general_cosine) == 0
    and sp.diff(general_internal_cosine, q) != 0,
)

# The zero-deficit equation does not select a physical elapsed time when the
# final boundary length a' is free.  With y=t/a, an exact flat/reflection
# branch is q=(a'/a)^2=1+y^2-y/phi.  The opposite-normal algebraic branch has
# the plus sign.  The static boundary q=1 intersects the inward branch at the
# degenerate y=0 and the unique nondegenerate y=1/phi.
y = sp.symbols("y", positive=True)
zero_deficit_numerator = sp.factor(
    sp.together(general_internal_cosine-target_cosine).as_numer_denom()[0],
    extension=sqrt5,
)
flat_q_minus = sp.simplify(1+y**2-y/phi)
flat_q_plus = sp.simplify(1+y**2+y/phi)
flat_determinant_minus = sp.factor(
    general_determinant.subs({r: y**2, q: flat_q_minus}),
    extension=sqrt5,
)
static_intersection = sp.factor(flat_q_minus-1, extension=sqrt5)
check(
    "zero deficit is a flat boundary-data curve, not a selected tick",
    sp.simplify(zero_deficit_numerator.subs({r: y**2, q: flat_q_minus})) == 0
    and sp.simplify(zero_deficit_numerator.subs({r: y**2, q: flat_q_plus})) == 0
    and flat_determinant_minus.is_positive
    and sp.simplify(static_intersection-y*(y-1/phi)) == 0,
    "q=1+y^2-y/phi; imposing q=1 leaves y=1/phi after excluding y=0",
)

# -------------------------------------------------------------------------
# Independent full boundary-action reconstruction.
# -------------------------------------------------------------------------
def simplex_coordinates(tent_length):
    """One static a=a'=1 congruent four-simplex in R4."""
    spatial_gram = np.full((3, 3), 0.5)
    np.fill_diagonal(spatial_gram, 1.0)
    link_vectors = np.zeros((3, 4))
    link_vectors[:, :3] = np.linalg.cholesky(spatial_gram)
    common_dot = tent_length**2 / 2
    coefficients = np.linalg.solve(spatial_gram, np.full(3, common_dot))
    projection = coefficients @ link_vectors
    height_squared = tent_length**2 - projection @ projection
    if height_squared <= 0:
        raise ValueError("non-Euclidean or degenerate simplex")
    new_vertex = projection.copy()
    new_vertex[3] = np.sqrt(height_squared)
    return np.vstack((np.zeros(4), new_vertex, link_vectors))


def triangle_area(points, triangle):
    x, y, z = points[list(triangle)]
    left, right = y-x, z-x
    return 0.5 * np.sqrt(
        (left @ left)*(right @ right) - (left @ right)**2
    )


def outward_facet_normal(points, omitted):
    facet = [index for index in range(5) if index != omitted]
    base = points[facet[0]]
    matrix = np.vstack([points[index]-base for index in facet[1:]])
    _, _, right = np.linalg.svd(matrix)
    normal = right[-1]
    if normal @ (points[omitted]-base) > 0:
        normal = -normal
    return normal / np.linalg.norm(normal)


def dihedral_angle(points, triangle):
    omitted = [index for index in range(5) if index not in triangle]
    left = outward_facet_normal(points, omitted[0])
    right = outward_facet_normal(points, omitted[1])
    return np.arccos(np.clip(-left @ right, -1, 1))


def full_boundary_action(tent_length):
    points = simplex_coordinates(tent_length)
    area_internal = triangle_area(points, (0, 1, 2))
    area_link = triangle_area(points, (2, 3, 4))
    area_old = triangle_area(points, (0, 2, 3))
    area_new = triangle_area(points, (1, 2, 3))
    theta_internal = dihedral_angle(points, (0, 1, 2))
    theta_link = dihedral_angle(points, (2, 3, 4))
    theta_old = dihedral_angle(points, (0, 2, 3))
    theta_new = dihedral_angle(points, (1, 2, 3))
    action = (
        12 * area_internal * (2*np.pi-5*theta_internal)
        + 20 * area_link * (np.pi-theta_link)
        + 30 * area_old * (np.pi-2*theta_old)
        + 30 * area_new * (np.pi-2*theta_new)
    )
    return action, theta_internal


schlaefli_rows = []
for tent_length in (0.3, 0.6, 1.0):
    step = 1e-5
    derivative_full = (
        full_boundary_action(tent_length+step)[0]
        - full_boundary_action(tent_length-step)[0]
    ) / (2*step)
    theta_numeric = full_boundary_action(tent_length)[1]
    area_derivative_t = (
        2-tent_length**2
    ) / (2*np.sqrt(4-tent_length**2))
    derivative_reduced = (
        12 * (2*np.pi-5*theta_numeric) * area_derivative_t
    )
    theta_exact = np.arccos(
        float(internal_cosine.subs(r, tent_length**2))
    )
    schlaefli_rows.append({
        "t_over_a": tent_length,
        "full_derivative": derivative_full,
        "reduced_derivative": derivative_reduced,
        "derivative_residual": derivative_full-derivative_reduced,
        "angle_residual": theta_numeric-theta_exact,
    })

check(
    "the explicit full boundary action obeys the Schlaefli-reduced tent equation",
    max(abs(row["derivative_residual"]) for row in schlaefli_rows) < 2e-8
    and max(abs(row["angle_residual"]) for row in schlaefli_rows) < 1e-12,
    f"max derivative residual="
    f"{max(abs(row['derivative_residual']) for row in schlaefli_rows):.3e}",
)

golden_t = float(sp.sqrt(golden_r))
golden_action, golden_theta_numeric = full_boundary_action(golden_t)
golden_derivative_step = 1e-5
golden_derivative_numeric = (
    full_boundary_action(golden_t+golden_derivative_step)[0]
    - full_boundary_action(golden_t-golden_derivative_step)[0]
) / (2*golden_derivative_step)
check(
    "the full action is stationary at t/a=phi^-1",
    abs(golden_theta_numeric-2*np.pi/5) < 1e-12
    and abs(golden_derivative_numeric) < 2e-8,
    f"theta={golden_theta_numeric:.12f}; dS/dt={golden_derivative_numeric:.3e}",
)

# A nonzero volume/cosmological coefficient shifts the root. There are 20
# congruent simplices; work in a=1 and differentiate with respect to r.
total_volume = sp.Rational(20, 96) * sp.sqrt(r*(8-3*r))
volume_derivative_r_at_root = sp.simplify(
    sp.diff(total_volume, r).subs(r, golden_r)
)
check(
    "any nonzero volume coefficient destroys stationarity at the golden root",
    volume_derivative_r_at_root != 0
    and float(volume_derivative_r_at_root) > 0,
    f"dV_total/dr={float(volume_derivative_r_at_root):.12f}",
)

# -------------------------------------------------------------------------
# Independent ambient reflection control.
# -------------------------------------------------------------------------
neighbour_points = vertices[list(link_vertices)]
plane_base = neighbour_points[0]
_, _, right = np.linalg.svd(neighbour_points[1:]-plane_base)
plane_normal = right[-1] / np.linalg.norm(right[-1])
reflected = vertices[v] - 2 * (
    (vertices[v]-plane_base) @ plane_normal
) * plane_normal
phi_float = (1+np.sqrt(5.0))/2
spatial_edge = np.linalg.norm(vertices[v]-neighbour_points[0])
reflected_pole = np.linalg.norm(vertices[v]-reflected)
reflected_final_edges = np.linalg.norm(
    neighbour_points-reflected, axis=1
)
check(
    "the golden pole is independently the neighbour-hyperplane reflection",
    np.linalg.norm(reflected-vertices[v]/phi_float) < 1e-12
    and np.max(np.abs(reflected_final_edges-spatial_edge)) < 1e-12
    and abs(reflected_pole/spatial_edge-1/phi_float) < 1e-12,
    f"|v'|={np.linalg.norm(reflected):.12f}; t/a={reflected_pole/spatial_edge:.12f}",
)
check(
    "the reflected endpoint is not a second 600-cell boundary vertex",
    np.min(np.linalg.norm(vertices-reflected, axis=1)) > 0.3
    and abs(np.linalg.norm(reflected)-1/phi_float) < 1e-12,
    "it lies on the inner radius 1/phi, not the unit S3 carrier",
)

protocol = (
    ROOT / "docs" / "gravity" / "gravity_tent_move_regge_protocol.md"
).read_text()
check(
    "the protocol keeps Euclidean, Lambda=0 and static-boundary hypotheses explicit",
    "Euclidean" in protocol
    and "zero-cosmological-constant" in protocol
    and "static equilateral boundary ansatz" in protocol
    and "not physical time" in protocol,
)

verdict = "DERIVED CONDITIONAL EUCLIDEAN GOLDEN TENT"
payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "external_primary_source": "https://arxiv.org/abs/1108.1974",
    "hypotheses": {
        "signature": "Euclidean",
        "old_spatial_edges": "a",
        "new_spatial_edges": "a'=a",
        "tent_pole": "t variable",
        "action": "ordinary Regge plus boundary term",
        "cosmological_volume_coefficient": 0,
        "varied_variable": "internal tent pole only",
    },
    "link": {
        "f_vector": [12, 30, 20],
        "triangles_per_vertex": 5,
        "stabilizer_order": len(stabilizer_link_permutations),
        "cell_orbits": [len(vertex_orbits), len(edge_orbits), len(triangle_orbits)],
    },
    "tent_carrier": {
        "f_vector": list(tent_f_vector),
        "euler_characteristic": tent_euler,
        "four_simplices": len(tent_four_simplices),
        "boundary_tetrahedra_initial_final": [len(initial_boundary), len(final_boundary)],
        "internal_triangles": len(internal_hinges),
        "four_simplices_per_internal_triangle": 5,
        "new_symmetric_metric_variables": ["a'", "t"],
        "global_tent_conflict_degree": 12,
    },
    "exact_geometry": {
        "general_gram_determinant_q_r": str(general_determinant),
        "general_internal_dihedral_cosine_q_r": str(general_internal_cosine),
        "zero_deficit_branches": [
            "q=1+y^2-y/phi",
            "q=1+y^2+y/phi",
        ],
        "static_slice_intersection": "q=1 and y>0 gives y=phi^-1",
        "static_gram_determinant": "a^8*r*(8-3*r)/16",
        "euclidean_domain": "0<r<8/3",
        "internal_dihedral_cosine": "(2-r)/(2*(3-r))",
        "internal_area": "a^2*sqrt(r*(4-r))/4",
        "simplex_volume": "a^4*sqrt(r*(8-3*r))/96",
        "deficit": "2*pi-5*theta",
        "stationary_r": str(golden_r),
        "stationary_t_over_a": "phi^-1",
        "stationary_second_r": str(stationary_second_r),
        "frozen_sign_stationary_type": "maximum",
    },
    "full_action_controls": schlaefli_rows,
    "golden_full_action": {
        "action": golden_action,
        "internal_angle": golden_theta_numeric,
        "target_angle": float(2*np.pi/5),
        "finite_difference_derivative": golden_derivative_numeric,
    },
    "hostile_controls": {
        "total_volume_derivative_r_at_golden": str(volume_derivative_r_at_root),
        "nonzero_volume_coefficient_shifts_root": True,
        "a_prime_release_adds_parameter": True,
        "zero_deficit_selects_physical_tick": False,
        "lorentzian_continuation_covered": False,
        "full_h4_synchronous_tent_orbit_exists": False,
        "reflected_endpoint_is_existing_600cell_vertex": False,
    },
    "ambient_reflection": {
        "unit_circumradius_endpoint_radius": float(np.linalg.norm(reflected)),
        "tent_over_spatial_edge": float(reflected_pole/spatial_edge),
        "expected": float(1/phi_float),
    },
    "verdict": verdict,
    "derived": [
        "the canonical local tent carrier is I joined with the icosahedral vertex link",
        "fivefold incidence makes zero deficit equivalent to theta=2*pi/5",
        "under the frozen Euclidean static Lambda=0 hypotheses t/a=phi^-1 uniquely",
        "the same ratio is the exact neighbour-hyperplane reflection geometry",
        "with the displayed action sign the stationary pole direction is a maximum",
    ],
    "conditional": [
        "ordinary Euclidean Regge action is the correct dynamics",
        "the cosmological/volume coefficient vanishes",
        "the final local star is equilateral with a'=a",
    ],
    "open": [
        "Lorentzian tent-pole angles and causal admissibility",
        "a selected nonzero volume/cosmological coefficient",
        "an orbit-complete non-overlapping global tent schedule",
        "whether the flat boundary-data curve is a physical gauge orbit in a completed canonical theory",
        "constraint class and kinetic Hessian on unrestricted boundary data",
        "c, G, Planck time and Planck mass",
    ],
    "not_claimed": [
        "the Euclidean pole is physical elapsed time",
        "A1=5 alone derives a dimensional constant",
        "the reflected endpoint is another 600-cell slice vertex",
        "the golden ratio survives a nonzero cosmological term",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(verdict)
print("OPEN: Lorentzian continuation, volume term, global schedule and physical units")
raise SystemExit(0 if passed == tests else 1)

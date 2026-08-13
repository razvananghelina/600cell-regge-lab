#!/usr/bin/env python3
"""Symmetry-reduced Lorentzian Regge action on the global tent slab.

Protocol commit: 06a1c6a.  This verifier constructs and validates the
35-variable invariant action, including full-versus-orbit controls and its
linearization at regular data.  It deliberately performs no root search.
"""

from collections import Counter, defaultdict
import cmath
from itertools import combinations
import json
import math
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_global_regge_orbits.json"
PROTOCOL_COMMIT = "06a1c6a"
CORRECTION_COMMIT = "08b638c"
TRACE_CORRECTION_COMMIT = "4a28c41"
TRACE_STEP_CORRECTION_COMMIT = "074a82b"
UPSTREAM_COMMIT = "d439b07"
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
    a, b, c, d = left
    e, f, g, h = right
    return np.array((
        a*e-b*f-c*g-d*h,
        a*f+b*e+c*h-d*g,
        a*g-b*h+c*e+d*f,
        a*h+b*g-c*f+d*e,
    ))


def canonical_key(item):
    return tuple(sorted(item))


def build_tetrahedra(adjacency):
    neighbours = [set(np.flatnonzero(row).tolist()) for row in adjacency]
    tetrahedra = []
    for a in range(120):
        for b in sorted(vertex for vertex in neighbours[a] if vertex > a):
            for c in sorted(
                vertex for vertex in neighbours[a] & neighbours[b]
                if vertex > b
            ):
                for d in sorted(
                    vertex
                    for vertex in neighbours[a] & neighbours[b] & neighbours[c]
                    if vertex > c
                ):
                    tetrahedra.append((a, b, c, d))
    return neighbours, tuple(tetrahedra)


def build_slab(tetrahedra, phase):
    simplices = set()
    for tetrahedron in tetrahedra:
        ordered = sorted(tetrahedron, key=phase.__getitem__)
        if len({phase[vertex] for vertex in ordered}) != 4:
            raise RuntimeError("repeated phase color in spatial tetrahedron")
        for vertex in ordered:
            simplex = [vertex, vertex + 120]
            simplex.extend(
                other + 120 if phase[other] < phase[vertex] else other
                for other in ordered if other != vertex
            )
            simplices.add(tuple(sorted(simplex)))
    return frozenset(simplices)


def extended_image(action, vertex):
    return int(action[vertex] if vertex < 120 else action[vertex - 120] + 120)


def orbit_partition(items, stabilizer):
    item_set = set(items)
    unseen = set(items)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = frozenset(
            tuple(sorted(extended_image(action, vertex) for vertex in seed))
            for action in stabilizer
        )
        if not orbit <= item_set:
            raise RuntimeError("stabilizer left the simplex layer")
        orbits.append(orbit)
        unseen -= orbit
    return tuple(orbits)


def log_minus(value):
    """Principal logarithm with the negative-real cut approached below."""
    scale = max(1.0, abs(value))
    if abs(value.imag) < 2e-13*scale:
        value = complex(value.real, 0.0)
        if value.real < 0:
            return complex(math.log(-value.real), -math.pi)
    return cmath.log(value)


def signed_volume_square(squared, local_vertices):
    vertices = list(local_vertices)
    dimension = len(vertices) - 1
    if dimension == 0:
        return 1.0
    base = vertices[0]
    others = vertices[1:]
    gram = np.array([
        [
            (squared[base, left] + squared[base, right]
             - squared[left, right]) / 2
            for right in others
        ]
        for left in others
    ])
    return float(np.linalg.det(gram)/(math.factorial(dimension)**2))


def angle_data(squared):
    """Return the ten corrected theta+ angles for one Lorentzian 4-simplex."""
    gram = np.array([
        [
            (squared[0, left] + squared[0, right] - squared[left, right])/2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    inverse = np.linalg.inv(gram)
    simplex_volume_square = signed_volume_square(squared, range(5))
    facet_volume_squares = {
        omitted: signed_volume_square(
            squared, [vertex for vertex in range(5) if vertex != omitted]
        )
        for omitted in range(5)
    }
    angles = {}
    arguments = {}
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge_vertices = [
            vertex for vertex in range(5)
            if vertex not in (omitted_a, omitted_b)
        ]
        hinge_volume_square = signed_volume_square(squared, hinge_vertices)
        gram_derivative = np.zeros((4, 4))
        opposite_edge = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                hit_0_left = {0, left} == opposite_edge
                hit_0_right = {0, right} == opposite_edge
                hit_pair = left != right and {left, right} == opposite_edge
                gram_derivative[left-1, right-1] = (
                    float(hit_0_left) + float(hit_0_right) - float(hit_pair)
                )/2
        volume_derivative = simplex_volume_square*np.trace(
            inverse @ gram_derivative
        )
        denominator = (
            cmath.sqrt(complex(facet_volume_squares[omitted_a]))
            * cmath.sqrt(complex(facet_volume_squares[omitted_b]))
        )
        cosine = 16*volume_derivative/denominator
        sine = -(4/3)*(
            cmath.sqrt(complex(hinge_volume_square))
            * cmath.sqrt(complex(simplex_volume_square))
        )/denominator
        if abs(cosine*cosine + sine*sine - 1) > 3e-9:
            raise RuntimeError("complex sine/cosine identity failed")
        argument = cosine + 1j*sine
        angles[tuple(hinge_vertices)] = -1j*log_minus(argument)
        arguments[tuple(hinge_vertices)] = argument
    eigenvalues = np.linalg.eigvalsh(gram)
    return angles, arguments, eigenvalues


def triangle_area_square(edge_values):
    x, y, z = edge_values
    return (2*(x*y + x*z + y*z) - x*x - y*y - z*z)/16


def triangle_area_partials(edge_values):
    x, y, z = edge_values
    return ((y+z-x)/8, (x+z-y)/8, (x+y-z)/8)


def relative_error(left, right):
    return abs(left-right)/max(1.0, abs(left), abs(right))


def construct_model(vertices, adjacency, tetrahedra, multiplication, conjugate,
                    h4_actions, cover, ordering):
    phase = {
        vertex: phase_index
        for phase_index, cell in enumerate(ordering)
        for vertex in cell
    }
    slab = build_slab(tetrahedra, phase)
    stabilizer = tuple(
        action for action in h4_actions
        if all(
            frozenset(int(action[vertex]) for vertex in cell) == cell
            for cell in ordering
        )
    )
    if len(stabilizer) != 24:
        raise RuntimeError("ordered schedule stabilizer is not order 24")
    if not all(
        frozenset(
            tuple(sorted(extended_image(action, vertex) for vertex in simplex))
            for simplex in slab
        ) == slab
        for action in stabilizer
    ):
        raise RuntimeError("ordered stabilizer does not preserve the slab")

    faces = {
        size: frozenset(
            tuple(sorted(face))
            for simplex in slab for face in combinations(simplex, size)
        )
        for size in (2, 3, 4)
    }
    old_edges = frozenset(
        tuple(sorted((left, right)))
        for left in range(120) for right in range(left + 1, 120)
        if adjacency[left, right]
    )
    new_edges = frozenset(
        (left+120, right+120) for left, right in old_edges
    )
    internal_edges = faces[2] - old_edges - new_edges
    old_triangles = frozenset(
        tuple(sorted(face))
        for tetrahedron in tetrahedra for face in combinations(tetrahedron, 3)
    )
    new_triangles = frozenset(
        tuple(vertex+120 for vertex in triangle) for triangle in old_triangles
    )
    boundary_triangles = old_triangles | new_triangles
    internal_triangles = faces[3] - boundary_triangles
    facet_multiplicity = Counter(
        tuple(sorted(facet))
        for simplex in slab for facet in combinations(simplex, 4)
    )
    internal_facets = frozenset(
        facet for facet, multiplicity in facet_multiplicity.items()
        if multiplicity == 2
    )

    edge_orbits_unsorted = orbit_partition(internal_edges, stabilizer)

    def edge_orbit_key(orbit):
        representative = min(orbit)
        if representative[1] - representative[0] == 120:
            return (1, phase[representative[0]], canonical_key(orbit))
        logical = tuple(vertex % 120 for vertex in representative)
        phase_pair = tuple(sorted(phase[vertex] for vertex in logical))
        return (0, phase_pair, canonical_key(orbit))

    edge_orbits = tuple(sorted(edge_orbits_unsorted, key=edge_orbit_key))
    diagonal_orbits = tuple(
        orbit for orbit in edge_orbits
        if all(right-left != 120 for left, right in orbit)
    )
    pole_orbits = tuple(
        orbit for orbit in edge_orbits
        if all(right-left == 120 for left, right in orbit)
    )
    edge_orbits = diagonal_orbits + pole_orbits
    edge_to_variable = {
        edge: variable
        for variable, orbit in enumerate(edge_orbits) for edge in orbit
    }
    edge_jacobian = {
        edge: (-1.0 if edge[1]-edge[0] == 120 else 1.0)
        for edge in internal_edges
    }
    individual_edges = tuple(sorted(internal_edges))
    individual_edge_position = {
        edge: position for position, edge in enumerate(individual_edges)
    }

    triangle_orbits = orbit_partition(faces[3], stabilizer)
    triangle_orbits = tuple(sorted(triangle_orbits, key=canonical_key))
    triangle_to_orbit = {
        triangle: orbit_index
        for orbit_index, orbit in enumerate(triangle_orbits)
        for triangle in orbit
    }
    internal_triangle_orbits = tuple(
        orbit for orbit in triangle_orbits if min(orbit) in internal_triangles
    )
    simplex_orbits = orbit_partition(slab, stabilizer)
    simplex_orbits = tuple(sorted(simplex_orbits, key=canonical_key))
    internal_facet_orbits = orbit_partition(internal_facets, stabilizer)

    return {
        "phase": phase,
        "cover": cover,
        "ordering": ordering,
        "slab": slab,
        "stabilizer": stabilizer,
        "faces": faces,
        "old_edges": old_edges,
        "new_edges": new_edges,
        "internal_edges": internal_edges,
        "boundary_triangles": boundary_triangles,
        "internal_triangles": internal_triangles,
        "internal_facets": internal_facets,
        "edge_orbits": edge_orbits,
        "diagonal_orbits": diagonal_orbits,
        "pole_orbits": pole_orbits,
        "edge_to_variable": edge_to_variable,
        "edge_jacobian": edge_jacobian,
        "individual_edges": individual_edges,
        "individual_edge_position": individual_edge_position,
        "triangle_orbits": triangle_orbits,
        "triangle_to_orbit": triangle_to_orbit,
        "internal_triangle_orbits": internal_triangle_orbits,
        "simplex_orbits": simplex_orbits,
        "internal_facet_orbits": internal_facet_orbits,
    }


def edge_square(model, edge, variables):
    edge = tuple(sorted(edge))
    if edge in model["edge_to_variable"]:
        variable = model["edge_to_variable"][edge]
        return model["edge_jacobian"][edge]*variables[variable]
    if edge in model["old_edges"] or edge in model["new_edges"]:
        return 1.0
    raise ValueError(f"edge absent from slab: {edge}")


def simplex_squared(model, simplex, variables):
    squared = np.zeros((5, 5))
    for left, right in combinations(range(5), 2):
        value = edge_square(model, (simplex[left], simplex[right]), variables)
        squared[left, right] = squared[right, left] = value
    return squared


def triangle_area_data(model, triangle, variables):
    edges = tuple(tuple(sorted(edge)) for edge in combinations(triangle, 2))
    values = tuple(edge_square(model, edge, variables) for edge in edges)
    area_square = triangle_area_square(values)
    derivative_by_edge = {}
    for edge, partial in zip(edges, triangle_area_partials(values)):
        if edge in model["edge_to_variable"]:
            derivative_by_edge[edge] = partial*model["edge_jacobian"][edge]
    return area_square, derivative_by_edge


def hinge_constant(model, triangle):
    return math.pi if triangle in model["boundary_triangles"] else 2*math.pi


def reduced_evaluation(model, variables, diagnostics=False):
    triangle_orbits = model["triangle_orbits"]
    curvature = np.array([
        hinge_constant(model, min(orbit)) for orbit in triangle_orbits
    ], dtype=complex)
    minimum_argument = math.inf
    minimum_absolute_gram_eigenvalue = math.inf
    negative_counts = Counter()
    angle_records = []
    for simplex_orbit in model["simplex_orbits"]:
        simplex = min(simplex_orbit)
        squared = simplex_squared(model, simplex, variables)
        angles, arguments, eigenvalues = angle_data(squared)
        negative_counts[int(np.sum(eigenvalues < -1e-10))] += 1
        minimum_absolute_gram_eigenvalue = min(
            minimum_absolute_gram_eigenvalue,
            float(np.min(np.abs(eigenvalues))),
        )
        for local_hinge, angle in angles.items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            orbit_index = model["triangle_to_orbit"][triangle]
            curvature[orbit_index] += angle
            minimum_argument = min(minimum_argument, abs(arguments[local_hinge]))
            if diagnostics:
                angle_records.append((simplex, triangle, angle))

    action_sum = 0j
    gradient = np.zeros(35, dtype=complex)
    for orbit_index, orbit in enumerate(triangle_orbits):
        triangle = min(orbit)
        area_square, derivative_by_edge = triangle_area_data(
            model, triangle, variables
        )
        root_area = cmath.sqrt(complex(area_square))
        action_sum += 24*root_area*curvature[orbit_index]
        for edge, area_square_derivative in derivative_by_edge.items():
            variable = model["edge_to_variable"][edge]
            gradient[variable] += (
                -1j*24*curvature[orbit_index]
                * area_square_derivative/(2*root_area)
            )
    return -1j*action_sum, gradient, {
        "curvature": curvature,
        "minimum_argument": minimum_argument,
        "minimum_absolute_gram_eigenvalue": minimum_absolute_gram_eigenvalue,
        "negative_counts": negative_counts,
        "angle_records": angle_records,
    }


def full_evaluation(model, variables):
    angle_incidence = defaultdict(list)
    minimum_argument = math.inf
    negative_counts = Counter()
    for simplex in model["slab"]:
        squared = simplex_squared(model, simplex, variables)
        angles, arguments, eigenvalues = angle_data(squared)
        negative_counts[int(np.sum(eigenvalues < -1e-10))] += 1
        for local_hinge, angle in angles.items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            angle_incidence[triangle].append(angle)
            minimum_argument = min(minimum_argument, abs(arguments[local_hinge]))
    curvature = {
        triangle: hinge_constant(model, triangle) + sum(angles)
        for triangle, angles in angle_incidence.items()
    }
    action_sum = 0j
    individual_gradient = np.zeros(840, dtype=complex)
    for triangle, triangle_curvature in curvature.items():
        area_square, derivative_by_edge = triangle_area_data(
            model, triangle, variables
        )
        root_area = cmath.sqrt(complex(area_square))
        action_sum += root_area*triangle_curvature
        for edge, area_square_derivative in derivative_by_edge.items():
            position = model["individual_edge_position"][edge]
            individual_gradient[position] += (
                -1j*triangle_curvature*area_square_derivative/(2*root_area)
            )
    restricted_gradient = np.array([
        sum(
            individual_gradient[model["individual_edge_position"][edge]]
            for edge in orbit
        )
        for orbit in model["edge_orbits"]
    ])
    orbit_curvature = np.empty(len(model["triangle_orbits"]), dtype=complex)
    maximum_orbit_spread = 0.0
    for orbit_index, orbit in enumerate(model["triangle_orbits"]):
        values = [curvature[triangle] for triangle in orbit]
        orbit_curvature[orbit_index] = sum(values)/len(values)
        maximum_orbit_spread = max(
            maximum_orbit_spread,
            max(abs(value-orbit_curvature[orbit_index]) for value in values),
        )
    maximum_edge_gradient_spread = 0.0
    for orbit in model["edge_orbits"]:
        values = [
            individual_gradient[model["individual_edge_position"][edge]]
            for edge in orbit
        ]
        mean = sum(values)/len(values)
        maximum_edge_gradient_spread = max(
            maximum_edge_gradient_spread,
            max(abs(value-mean) for value in values),
        )
    return -1j*action_sum, restricted_gradient, {
        "curvature": orbit_curvature,
        "minimum_argument": minimum_argument,
        "negative_counts": negative_counts,
        "maximum_orbit_curvature_spread": maximum_orbit_spread,
        "maximum_edge_gradient_spread": maximum_edge_gradient_spread,
        "individual_gradient": individual_gradient,
    }


def control_vectors():
    regular = np.r_[np.ones(30), np.full(5, 0.25)]
    upward = np.r_[
        np.array([1 + (index+1)/1000 for index in range(30)]),
        np.array([0.25 + (index+1)/1000 for index in range(5)]),
    ]
    downward = np.r_[
        np.array([1 - (index+1)/2000 for index in range(30)]),
        np.array([0.25 + (5-index)/1500 for index in range(5)]),
    ]
    return {"R0": regular, "R1": upward, "R2": downward}


print("=" * 78)
print("GLOBAL LORENTZIAN REGGE ORBIT REDUCTION")
print("=" * 78)

vertices, adjacency_float, _ = build_600cell()
adjacency = adjacency_float > 0.5
neighbours, tetrahedra = build_tetrahedra(adjacency)

multiplication = np.empty((120, 120), dtype=np.int16)
for left in range(120):
    for right in range(120):
        multiplication[left, right] = int(np.argmax(
            vertices @ qmul(vertices[left], vertices[right])
        ))
conjugate = np.array([
    int(np.argmax(vertices @ (vertex*np.array((1, -1, -1, -1)))))
    for vertex in vertices
], dtype=np.int16)

binary_tetrahedral = frozenset(
    index for index, vertex in enumerate(vertices)
    if (
        np.count_nonzero(np.abs(vertex) > 1e-8) == 1
        and np.max(np.abs(vertex)) > 1-1e-8
    ) or np.all(np.abs(np.abs(vertex)-0.5) < 1e-8)
)
unseen = set(range(120))
cover_cells = []
while unseen:
    representative = min(unseen)
    cell = frozenset(
        int(multiplication[representative, element])
        for element in binary_tetrahedral
    )
    cover_cells.append(cell)
    unseen -= cell
cover_cells = tuple(sorted(cover_cells, key=canonical_key))
cover = frozenset(cover_cells)

h4_actions_set = set()
plain = np.arange(120, dtype=np.int16)
for reflected in (False, True):
    seed = conjugate if reflected else plain
    for left in range(120):
        left_images = multiplication[left, seed]
        for right in range(120):
            h4_actions_set.add(tuple(
                int(value)
                for value in multiplication[left_images, conjugate[right]]
            ))
h4_actions = tuple(np.asarray(action, dtype=np.int16)
                   for action in h4_actions_set)

even_ordering = cover_cells
odd_ordering = (
    cover_cells[1], cover_cells[0], cover_cells[2], cover_cells[3], cover_cells[4]
)
models = {
    "even": construct_model(
        vertices, adjacency, tetrahedra, multiplication, conjugate,
        h4_actions, cover, even_ordering,
    ),
    "odd": construct_model(
        vertices, adjacency, tetrahedra, multiplication, conjugate,
        h4_actions, cover, odd_ordering,
    ),
}

for parity, model in models.items():
    check(
        f"{parity}: ordered slab stabilizer has exact order 24",
        len(model["stabilizer"]) == 24,
    )
    edge_orbit_sizes = Counter(map(len, model["edge_orbits"]))
    triangle_orbit_sizes = Counter(map(len, model["internal_triangle_orbits"]))
    facet_orbit_sizes = Counter(map(len, model["internal_facet_orbits"]))
    simplex_orbit_sizes = Counter(map(len, model["simplex_orbits"]))
    check(
        f"{parity}: 840 internal edges reduce to 30 diagonal plus 5 pole orbits",
        len(model["internal_edges"]) == 840
        and len(model["diagonal_orbits"]) == 30
        and len(model["pole_orbits"]) == 5
        and edge_orbit_sizes == Counter({24: 35}),
    )
    diagonal_phase_pairs = Counter()
    for orbit in model["diagonal_orbits"]:
        representative = min(orbit)
        logical = tuple(vertex % 120 for vertex in representative)
        phase_pair = tuple(sorted(model["phase"][vertex] for vertex in logical))
        diagonal_phase_pairs[phase_pair] += 1
    check(
        f"{parity}: every phase pair has three distinct diagonal orbits",
        diagonal_phase_pairs == Counter({pair: 3 for pair in combinations(range(5), 2)}),
    )
    check(
        f"{parity}: internal hinge/facet/simplex orbit counts reproduce 160/225/100",
        len(model["internal_triangles"]) == 3840
        and len(model["internal_triangle_orbits"]) == 160
        and triangle_orbit_sizes == Counter({24: 160})
        and len(model["internal_facets"]) == 5400
        and len(model["internal_facet_orbits"]) == 225
        and facet_orbit_sizes == Counter({24: 225})
        and len(model["slab"]) == 2400
        and len(model["simplex_orbits"]) == 100
        and simplex_orbit_sizes == Counter({24: 100}),
    )

controls = control_vectors()
control_results = {parity: {} for parity in models}
for parity, model in models.items():
    for control_name, variables in controls.items():
        reduced_action, reduced_gradient, reduced_data = reduced_evaluation(
            model, variables
        )
        full_action, full_gradient, full_data = full_evaluation(model, variables)
        action_error = relative_error(reduced_action, full_action)
        gradient_error = max(
            relative_error(left, right)
            for left, right in zip(reduced_gradient, full_gradient)
        )
        curvature_error = max(
            relative_error(left, right)
            for left, right in zip(
                reduced_data["curvature"], full_data["curvature"]
            )
        )
        real_residual = max(
            abs(reduced_action.imag),
            float(np.max(np.abs(reduced_gradient.imag))),
        )
        causal = (
            reduced_data["negative_counts"] == Counter({1: 100})
            and full_data["negative_counts"] == Counter({1: 2400})
        )
        check(
            f"{parity}/{control_name}: full and 100-orbit actions/gradients agree",
            action_error < 2e-8
            and gradient_error < 2e-8
            and curvature_error < 2e-9
            and full_data["maximum_orbit_curvature_spread"] < 2e-9
            and full_data["maximum_edge_gradient_spread"] < 2e-9,
            f"action={action_error:.2e}, gradient={gradient_error:.2e}, "
            f"curvature={curvature_error:.2e}",
        )
        check(
            f"{parity}/{control_name}: all simplices remain on a real Lorentzian branch",
            causal and real_residual < 2e-8
            and reduced_data["minimum_argument"] > 1e-5,
            f"real residual={real_residual:.2e}, "
            f"min |angle argument|={reduced_data['minimum_argument']:.3e}",
        )
        control_results[parity][control_name] = {
            "variables": variables,
            "action": reduced_action,
            "gradient": reduced_gradient,
            "minimum_argument": reduced_data["minimum_argument"],
            "full_action_error": action_error,
            "full_gradient_error": gradient_error,
            "full_curvature_error": curvature_error,
            "real_residual": real_residual,
        }

# Direct complete-action differences and per-simplex Schlaefli controls.
finite_difference_step = 1e-6
schlaefli_step = 2e-6
hessian_step = 5e-4
direct_errors = {}
schlaefli_errors = {}
regular_perturbed_gradients = {}
regular_perturbed_actions = {}
for parity, model in models.items():
    direct_errors[parity] = {}
    schlaefli_errors[parity] = {}
    regular_perturbed_gradients[parity] = {}
    regular_perturbed_actions[parity] = {}
    for control_name, variables in controls.items():
        analytic_gradient = control_results[parity][control_name]["gradient"]
        direct_gradient = np.empty(35, dtype=complex)
        difference_step = (
            hessian_step if control_name == "R0" else finite_difference_step
        )
        for variable in range(35):
            plus = variables.copy()
            minus = variables.copy()
            plus[variable] += difference_step
            minus[variable] -= difference_step
            plus_action, plus_gradient, _ = reduced_evaluation(model, plus)
            minus_action, minus_gradient, _ = reduced_evaluation(model, minus)
            direct_gradient[variable] = (
                plus_action-minus_action
            )/(2*difference_step)
            if control_name == "R0":
                regular_perturbed_gradients[parity][variable] = (
                    plus_gradient, minus_gradient
                )
                regular_perturbed_actions[parity][variable] = (
                    plus_action, minus_action
                )
        maximum_direct_error = max(
            relative_error(left, right)
            for left, right in zip(direct_gradient, analytic_gradient)
        )
        direct_errors[parity][control_name] = maximum_direct_error
        check(
            f"{parity}/{control_name}: direct complete-action differences reproduce 35 derivatives",
            maximum_direct_error < 2e-5,
            f"max relative error={maximum_direct_error:.3e}",
        )

        schlaefli_residuals = []
        for simplex_orbit in model["simplex_orbits"]:
            simplex = min(simplex_orbit)
            variables_in_simplex = sorted({
                model["edge_to_variable"][tuple(sorted(edge))]
                for edge in combinations(simplex, 2)
                if tuple(sorted(edge)) in model["edge_to_variable"]
            })
            variable = variables_in_simplex[0]
            plus = variables.copy()
            minus = variables.copy()
            plus[variable] += schlaefli_step
            minus[variable] -= schlaefli_step
            base_squared = simplex_squared(model, simplex, variables)
            plus_angles = angle_data(simplex_squared(model, simplex, plus))[0]
            minus_angles = angle_data(simplex_squared(model, simplex, minus))[0]
            residual = 0j
            for local_hinge in plus_angles:
                area_square = signed_volume_square(base_squared, local_hinge)
                angle_derivative = (
                    plus_angles[local_hinge]-minus_angles[local_hinge]
                )/(2*schlaefli_step)
                residual += cmath.sqrt(complex(area_square))*angle_derivative
            schlaefli_residuals.append(residual)
        maximum_schlaefli = max(map(abs, schlaefli_residuals))
        schlaefli_errors[parity][control_name] = maximum_schlaefli
        check(
            f"{parity}/{control_name}: corrected angles satisfy per-simplex Schlaefli",
            maximum_schlaefli < 2e-7,
            f"100 orbit representatives; max residual={maximum_schlaefli:.3e}",
        )

# The regular local no-go must appear in every phase pole orbit.
rho_regular = 0.25
regular_cosine = (2+rho_regular)/(2*(3+rho_regular))
regular_deficit = 2*math.pi - 5*math.acos(regular_cosine)
regular_pole_gradients = {}
for parity in models:
    pole_gradient = control_results[parity]["R0"]["gradient"][30:].real
    regular_pole_gradients[parity] = pole_gradient
    check(
        f"{parity}: all five regular pole-orbit derivatives reproduce the nonzero local obstruction",
        regular_deficit > 0.1
        and np.min(np.abs(pole_gradient)) > 1e-3,
        f"deficit={regular_deficit:.12f}, gradients={pole_gradient.tolist()}",
    )

# Centered Hessian of the restricted action at the nonstationary regular point.
# Reuse the already evaluated R0 perturbations from the direct-action control.
hessian_results = {}
for parity, model in models.items():
    hessian = np.empty((35, 35), dtype=complex)
    for variable in range(35):
        plus_gradient, minus_gradient = regular_perturbed_gradients[parity][variable]
        hessian[:, variable] = (
            plus_gradient-minus_gradient
        )/(2*hessian_step)
    imaginary_residual = float(np.max(np.abs(hessian.imag)))
    real_hessian = hessian.real
    symmetry_residual = float(
        np.linalg.norm(real_hessian-real_hessian.T)
        / max(1.0, np.linalg.norm(real_hessian))
    )
    singular_values = np.linalg.svd(real_hessian, compute_uv=False)
    eigenvalues = np.linalg.eigvalsh(real_hessian)
    determinant_sign, log_abs_determinant = np.linalg.slogdet(real_hessian)
    ranks = {
        str(threshold): int(np.sum(
            singular_values > threshold*singular_values[0]
        ))
        for threshold in (1e-7, 1e-9, 1e-11)
    }
    hessian_results[parity] = {
        "matrix": real_hessian,
        "imaginary_residual": imaginary_residual,
        "symmetry_residual": symmetry_residual,
        "singular_values": singular_values,
        "eigenvalues": eigenvalues,
        "trace": float(np.trace(real_hessian)),
        "frobenius_norm": float(np.linalg.norm(real_hessian)),
        "determinant_sign": float(determinant_sign),
        "log_abs_determinant": float(log_abs_determinant),
        "ranks": ranks,
    }
    check(
        f"{parity}: the 35x35 regular Hessian is real and symmetric",
        imaginary_residual < 2e-8 and symmetry_residual < 2e-5,
        f"imag={imaginary_residual:.2e}, symmetry={symmetry_residual:.2e}",
    )
    check(
        f"{parity}: regular Hessian rank is stable at frozen thresholds",
        len(set(ranks.values())) == 1,
        f"ranks={ranks}, smin/smax={singular_values[-1]/singular_values[0]:.3e}",
    )

# Independent trace from second differences of the complete reduced action.
# The plus/minus actions are the R0 evaluations already made above.
trace_step = hessian_step
direct_hessian_traces = {}
for parity, model in models.items():
    base_action = control_results[parity]["R0"]["action"]
    diagonal_second_derivatives = []
    for variable in range(35):
        plus_action, minus_action = regular_perturbed_actions[parity][variable]
        diagonal_second_derivatives.append(
            (plus_action-2*base_action+minus_action)/(trace_step**2)
        )
    direct_trace = sum(diagonal_second_derivatives)
    direct_hessian_traces[parity] = direct_trace
    trace_error = relative_error(direct_trace, hessian_results[parity]["trace"])
    check(
        f"{parity}: direct complete-action second differences reproduce Hessian trace",
        trace_error < 2e-5 and abs(direct_trace.imag) < 2e-4,
        f"direct={direct_trace.real:.12g}, Hessian={hessian_results[parity]['trace']:.12g}, "
        f"relative error={trace_error:.3e}",
    )

even_singular = hessian_results["even"]["singular_values"]
odd_singular = hessian_results["odd"]["singular_values"]
disclosed_small_step_singular = {
    "even": np.array([
        69.33508179112337, 65.68762441608266, 64.09592396755176,
        52.143230985556855, 49.01114435191091, 44.43474100230955,
        44.20036866790382, 44.17055408197872, 41.316118106439,
        39.726931482502486, 38.88002467311924, 38.19121901429888,
        36.92341806609214, 35.848997086036526, 34.78408725469383,
        34.346616289477026, 32.253226858138795, 29.03920216915895,
        28.058106986827003, 25.835928282565224, 24.903771334251537,
        24.759274310328994, 22.400645531715107, 21.356913983306313,
        20.597980933787657, 20.54895721116826, 19.534393624597975,
        19.360621123796278, 17.431480235506807, 16.953329069321374,
        3.414459213553449, 3.076800977522834, 2.7125882147528757,
        1.9551218633407204, 1.2536263934325031,
    ]),
    "odd": np.array([
        69.39912303761076, 65.22606922964981, 64.42124520835846,
        52.09549846854919, 48.15012946494844, 45.56132919980328,
        44.442265352887446, 44.02729634473625, 42.136268727525206,
        39.21988082005066, 38.60111221901721, 38.15382836816019,
        36.398831210643245, 35.64262027697267, 35.326980073787375,
        34.34423700219703, 32.300828372803075, 29.159710079808214,
        28.018475712033897, 25.98658421095778, 25.707311574263215,
        23.707495357323268, 22.449100602377914, 21.37306631661759,
        20.584320038681668, 20.382090447005112, 19.732057264600655,
        19.26458231150661, 17.10492473053301, 17.087336056421254,
        3.6441434524815843, 2.954333414226676, 2.856940896648947,
        1.9735923859184643, 1.0863340874169345,
    ]),
}
step_stability_errors = {}
for parity, current in (("even", even_singular), ("odd", odd_singular)):
    previous = disclosed_small_step_singular[parity]
    step_error = float(np.max(np.abs(current-previous))/previous[0])
    step_stability_errors[parity] = step_error
    check(
        f"{parity}: Hessian spectrum is stable from h=2e-5 to h=5e-4",
        step_error < 2e-4,
        f"relative max change={step_error:.3e}",
    )
parity_spectrum_error = float(
    np.max(np.abs(even_singular-odd_singular))
    / max(1.0, even_singular[0], odd_singular[0])
)
check(
    "the two phase-parity classes have decisively separated regular Hessian spectra",
    parity_spectrum_error > 1e-3,
    f"relative max spectrum difference={parity_spectrum_error:.3e}",
)

verdict = (
    "DERIVED COMPUTATIONAL: the honest fixed-symmetry global Regge problem "
    "has 35 variables (30 diagonal and 5 pole orbits).  Its 100-simplex-orbit "
    "action reproduces the unreduced 2400-simplex action and all restricted "
    "gradients at three causal controls.  Regular data remain nonstationary, "
    "and the two phase parities have distinct quadratic actions.  No root was "
    "searched in this protocol."
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "correction_commit": CORRECTION_COMMIT,
    "trace_resource_correction_commit": TRACE_CORRECTION_COMMIT,
    "trace_step_correction_commit": TRACE_STEP_CORRECTION_COMMIT,
    "disclosed_original_run": {
        "passed": 38,
        "tests": 39,
        "failed_assertion": "phase-parity Hessian spectra are equal",
        "observed_relative_separation": 0.016233464461549153,
    },
    "disclosed_second_failed_run": {
        "passed": 39,
        "tests": 41,
        "failed_assertions": "two direct trace imaginary residuals at h=2e-5",
        "even_imaginary_residual": 0.012098596436229125,
        "odd_imaginary_residual": 0.010138901928121979,
    },
    "upstream_schedule_commit": UPSTREAM_COMMIT,
    "variables": {
        "internal_edges": 840,
        "edge_orbits": 35,
        "diagonal_orbits": 30,
        "pole_orbits": 5,
    },
    "orbit_counts": {
        "internal_triangles": [3840, 160],
        "internal_tetrahedra": [5400, 225],
        "four_simplices": [2400, 100],
        "common_orbit_size": 24,
    },
    "controls": {
        parity: {
            name: {
                "action_real": float(record["action"].real),
                "action_imaginary": float(record["action"].imag),
                "gradient_norm": float(np.linalg.norm(record["gradient"].real)),
                "gradient_imaginary_residual": float(np.max(np.abs(record["gradient"].imag))),
                "minimum_angle_argument": record["minimum_argument"],
                "full_action_relative_error": record["full_action_error"],
                "full_gradient_relative_error": record["full_gradient_error"],
                "full_curvature_relative_error": record["full_curvature_error"],
                "direct_gradient_relative_error": direct_errors[parity][name],
                "schlaefli_residual": schlaefli_errors[parity][name],
            }
            for name, record in parity_records.items()
        }
        for parity, parity_records in control_results.items()
    },
    "regular": {
        "rho": rho_regular,
        "pole_deficit": regular_deficit,
        "pole_orbit_gradients": {
            parity: values.tolist()
            for parity, values in regular_pole_gradients.items()
        },
        "stationary": False,
    },
    "hessian": {
        parity: {
            "imaginary_residual": result["imaginary_residual"],
            "symmetry_residual": result["symmetry_residual"],
            "singular_values": result["singular_values"].tolist(),
            "eigenvalues": result["eigenvalues"].tolist(),
            "trace": result["trace"],
            "direct_action_trace": {
                "real": float(direct_hessian_traces[parity].real),
                "imaginary": float(direct_hessian_traces[parity].imag),
            },
            "frobenius_norm": result["frobenius_norm"],
            "determinant_sign": result["determinant_sign"],
            "log_abs_determinant": result["log_abs_determinant"],
            "ranks": result["ranks"],
            "smallest_over_largest": float(
                result["singular_values"][-1]/result["singular_values"][0]
            ),
        }
        for parity, result in hessian_results.items()
    },
    "phase_parity_hessian_spectrum_relative_error": parity_spectrum_error,
    "hessian_step_spectrum_stability_errors": step_stability_errors,
    "labels": {
        "35_variable_orbit_reduction": "DERIVED COMPUTATIONAL",
        "regular_global_stationarity": "DERIVED NEGATIVE",
        "phase_parity_quadratic_equivalence": "DERIVED NEGATIVE",
        "root_in_35_variable_fixed_subspace": "NOT SEARCHED",
        "root_in_full_840_variable_space": "OPEN",
        "physical_phase_parity": "OPEN",
    },
    "verdict": verdict,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(verdict)
print("OPEN: preregistered root search in 35 variables; full 840-variable roots")
if __name__ == "__main__":
    raise SystemExit(0 if passed == tests else 1)

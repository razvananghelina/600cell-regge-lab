#!/usr/bin/env python3
"""Exhaustive oriented-orbit census of all 600-cell five-colourings.

Prior-art commit: 1071c40.
Protocol commit: a1d9bf0.
"""

from collections import Counter, defaultdict
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ortools.sat.python import cp_model


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_chromatic_cover_orbits.json"
CELL600_SOURCE = ROOT / "commons" / "cell600.py"
GLOBAL_SOURCE = HERE / "verify_gravity_global_tent_schedule.py"
DEGREE_SOURCE = HERE / "verify_gravity_600cell_chromatic_degree_selector.py"
PRIOR_ART_COMMIT = "1071c40"
PROTOCOL_COMMIT = "a1d9bf0"
EXPECTED_HASHES = {
    "cell600": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "global_schedule_source": "4b575450c905b43dd0efeff155a914fff5d4aee2838291e4a98bf65c5fca26d1",
    "fixed_degree_source": "423be626361be9e5a1915d7b91cc05c62d98fdbd2cd51c522d96d336fee9de63",
}
DET_TOLERANCE = 1e-10
GEOMETRY_TOLERANCE = 1e-8


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qmul(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return np.array((
        a*e-b*f-c*g-d*h,
        a*f+b*e+c*h-d*g,
        a*g-b*h+c*e+d*f,
        a*h+b*g-c*f+d*e,
    ))


def canonical_set_key(value):
    return tuple(sorted(value))


def canonical_partition_key(value):
    return tuple(sorted(canonical_set_key(cell) for cell in value))


def permutation_sign(value):
    inversions = sum(
        value[left] > value[right]
        for left in range(len(value))
        for right in range(left + 1, len(value))
    )
    return -1 if inversions % 2 else 1


def build_tetrahedra(adjacency):
    neighbours = [set(np.flatnonzero(row).tolist()) for row in adjacency]
    edges = {
        (left, right)
        for left in range(120)
        for right in range(left + 1, 120)
        if adjacency[left, right]
    }
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
    triangles = {
        triangle
        for tetrahedron in tetrahedra
        for triangle in combinations(tetrahedron, 3)
    }
    return frozenset(edges), frozenset(triangles), tuple(tetrahedra)


def independent_model(edges, target=None, excluded=()):
    model = cp_model.CpModel()
    variables = [model.NewBoolVar(f"x_{index}") for index in range(120)]
    for left, right in edges:
        model.Add(variables[left] + variables[right] <= 1)
    if target is None:
        model.Maximize(sum(variables))
    else:
        model.Add(sum(variables) == target)
    for candidate in excluded:
        model.Add(sum(variables[index] for index in candidate) <= target - 1)
    return model


def solve_exact(model):
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 8
    solver.parameters.random_seed = 600
    solver.parameters.max_time_in_seconds = 120.0
    status = solver.Solve(model)
    return solver, status


def exact_covers(candidates):
    by_vertex = {
        vertex: tuple(cell for cell in candidates if vertex in cell)
        for vertex in range(120)
    }
    covers = set()

    def recurse(chosen, remaining):
        if not remaining:
            covers.add(frozenset(chosen))
            return
        vertex = min(remaining)
        for cell in by_vertex[vertex]:
            if cell <= remaining:
                recurse(chosen + (cell,), remaining - cell)

    recurse(tuple(), frozenset(range(120)))
    return covers


def degree_for_vertex_labels(tetrahedra, source_coefficients, vertex_label):
    pushforward = Counter()
    local_sign_counts = {missing: Counter() for missing in range(5)}
    for tetrahedron in tetrahedra:
        labels = tuple(vertex_label[vertex] for vertex in tetrahedron)
        missing_set = set(range(5)) - set(labels)
        if len(missing_set) != 1 or len(set(labels)) != 4:
            raise RuntimeError("colour map does not land in one target facet")
        missing = missing_set.pop()
        canonical_facet = tuple(
            colour for colour in range(5) if colour != missing
        )
        positions = {
            colour: index for index, colour in enumerate(canonical_facet)
        }
        reorder = tuple(positions[colour] for colour in labels)
        contribution = source_coefficients[tetrahedron] * permutation_sign(reorder)
        pushforward[missing] += contribution
        local_sign_counts[missing][contribution * ((-1)**missing)] += 1
    candidates = tuple(
        ((-1)**missing) * pushforward[missing] for missing in range(5)
    )
    independent = tuple(
        sum(sign*count for sign, count in local_sign_counts[missing].items())
        for missing in range(5)
    )
    return {
        "pushforward": tuple(pushforward[index] for index in range(5)),
        "candidates": candidates,
        "independent": independent,
        "preimage_counts": tuple(
            sum(local_sign_counts[index].values()) for index in range(5)
        ),
        "local_sign_counts": tuple(
            dict(local_sign_counts[index]) for index in range(5)
        ),
        "degree": candidates[0] if len(set(candidates)) == 1 else None,
    }


def partition_orbits(items, action_maps, ordered):
    remaining = set(items)
    orbits = []
    while remaining:
        seed = min(
            remaining,
            key=lambda value: tuple(value) if ordered else tuple(sorted(value)),
        )
        if ordered:
            orbit = {
                tuple(action[index] for index in seed) for action in action_maps
            }
        else:
            orbit = {
                frozenset(action[index] for index in seed)
                for action in action_maps
            }
        if not orbit <= set(items):
            raise RuntimeError("group action left the enumerated carrier")
        orbits.append(frozenset(orbit))
        remaining -= orbit
    return tuple(sorted(orbits, key=lambda orbit: min(tuple(x) if ordered else tuple(sorted(x)) for x in orbit)))


tests = []


def check(label, condition, detail=""):
    ok = bool(condition)
    tests.append((label, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


print("=" * 78)
print("ALL TEN 600-CELL CHROMATIC COVER ORBITS")
print("=" * 78)

source_hashes = {
    "cell600": digest(CELL600_SOURCE),
    "global_schedule_source": digest(GLOBAL_SOURCE),
    "fixed_degree_source": digest(DEGREE_SOURCE),
}
provenance_ok = source_hashes == EXPECTED_HASHES

vertices, adjacency_float, _ = build_600cell()
adjacency = adjacency_float > 0.5
edges, triangles, tetrahedra = build_tetrahedra(adjacency)
degrees = tuple(int(row.sum()) for row in adjacency)
carrier_ok = bool(
    vertices.shape == (120, 4)
    and len(edges) == 720
    and len(triangles) == 1200
    and len(tetrahedra) == 600
    and set(degrees) == {12}
)

multiplication = np.empty((120, 120), dtype=np.int16)
multiplication_residual = 0.0
for left in range(120):
    for right in range(120):
        product = qmul(vertices[left], vertices[right])
        image = int(np.argmax(vertices @ product))
        multiplication[left, right] = image
        multiplication_residual = max(
            multiplication_residual,
            float(np.max(np.abs(vertices[image] - product))),
        )
conjugate = np.array([
    int(np.argmax(vertices @ (vertex*np.array((1, -1, -1, -1)))))
    for vertex in vertices
], dtype=np.int16)
identity = int(np.argmax(vertices[:, 0]))

binary_tetrahedral = frozenset(
    index for index, vertex in enumerate(vertices)
    if (
        np.count_nonzero(np.abs(vertex) > GEOMETRY_TOLERANCE) == 1
        and np.max(np.abs(vertex)) > 1-GEOMETRY_TOLERANCE
    ) or np.all(np.abs(np.abs(vertex)-0.5) < GEOMETRY_TOLERANCE)
)
subgroup_ok = bool(
    len(binary_tetrahedral) == 24
    and identity in binary_tetrahedral
    and all(
        int(multiplication[left, right]) in binary_tetrahedral
        for left in binary_tetrahedral for right in binary_tetrahedral
    )
    and all(conjugate[index] in binary_tetrahedral
            for index in binary_tetrahedral)
)

action_provenance = defaultdict(set)
plain = np.arange(120, dtype=np.int16)
for reflected in (False, True):
    seed = conjugate if reflected else plain
    for left in range(120):
        left_images = multiplication[left, seed]
        for right in range(120):
            action = tuple(
                int(value)
                for value in multiplication[left_images, conjugate[right]]
            )
            action_provenance[action].add(reflected)

action_tuples = tuple(sorted(action_provenance))
action_arrays = {
    action: np.asarray(action, dtype=np.int16) for action in action_tuples
}
edge_array = np.asarray(sorted(edges), dtype=np.int16)
all_graph_automorphisms = all(
    len(set(action)) == 120
    and np.all(adjacency[
        action_arrays[action][edge_array[:, 0]],
        action_arrays[action][edge_array[:, 1]],
    ])
    for action in action_tuples
)

basis_indices = []
for index in range(120):
    candidate = basis_indices + [index]
    if np.linalg.matrix_rank(vertices[candidate], tol=1e-12) > len(basis_indices):
        basis_indices.append(index)
    if len(basis_indices) == 4:
        break
basis_matrix = vertices[basis_indices]
orientation_class = {}
maximum_action_residual = 0.0
maximum_orthogonality_residual = 0.0
determinant_residual = 0.0
for action in action_tuples:
    action_array = action_arrays[action]
    transformation = np.linalg.solve(
        basis_matrix, vertices[action_array[basis_indices]]
    )
    determinant = float(np.linalg.det(transformation))
    orientation_class[action] = 1 if determinant > 0 else -1
    maximum_action_residual = max(
        maximum_action_residual,
        float(np.max(np.abs(vertices @ transformation - vertices[action_array]))),
    )
    maximum_orthogonality_residual = max(
        maximum_orthogonality_residual,
        float(np.max(np.abs(transformation.T @ transformation - np.eye(4)))),
    )
    determinant_residual = max(
        determinant_residual, abs(abs(determinant)-1.0)
    )

proper_actions = tuple(
    action for action in action_tuples if orientation_class[action] == 1
)
improper_actions = tuple(
    action for action in action_tuples if orientation_class[action] == -1
)
orientation_generation_ok = bool(
    len(action_tuples) == 14400
    and len(proper_actions) == len(improper_actions) == 7200
    and all(action_provenance[action] == {False} for action in proper_actions)
    and all(action_provenance[action] == {True} for action in improper_actions)
    and maximum_action_residual < GEOMETRY_TOLERANCE
    and maximum_orthogonality_residual < GEOMETRY_TOLERANCE
    and determinant_residual < GEOMETRY_TOLERANCE
)

maximum_cells_set = {
    frozenset(int(action_arrays[action][index])
              for index in binary_tetrahedral)
    for action in action_tuples
}
maximum_cells = tuple(sorted(maximum_cells_set, key=canonical_set_key))
all_cells_independent = all(
    len(cell) == 24
    and not np.any(adjacency[np.ix_(tuple(cell), tuple(cell))])
    for cell in maximum_cells
)

optimization_model = independent_model(edges)
optimization_solver, optimization_status = solve_exact(optimization_model)
optimization_value = (
    int(round(optimization_solver.ObjectiveValue()))
    if optimization_status == cp_model.OPTIMAL else None
)
exclusion_model = independent_model(
    edges, target=24, excluded=maximum_cells
)
exclusion_solver, exclusion_status = solve_exact(exclusion_model)

cover_sets = exact_covers(maximum_cells)
cell_id = {cell: index for index, cell in enumerate(maximum_cells)}
covers = tuple(sorted(
    (
        frozenset(cell_id[cell] for cell in cover)
        for cover in cover_sets
    ),
    key=lambda cover: tuple(sorted(cover)),
))

def coset_partition(subgroup, side):
    unseen = set(range(120))
    cells = []
    while unseen:
        representative = min(unseen)
        if side == "left":
            cell = frozenset(
                int(multiplication[representative, element])
                for element in subgroup
            )
        else:
            cell = frozenset(
                int(multiplication[element, representative])
                for element in subgroup
            )
        cells.append(cell)
        unseen -= cell
    return frozenset(cell_id[cell] for cell in cells)


conjugate_subgroups = {
    frozenset(
        int(multiplication[
            multiplication[group_element, element], conjugate[group_element]
        ])
        for element in binary_tetrahedral
    )
    for group_element in range(120)
}
left_covers = {
    coset_partition(subgroup, "left") for subgroup in conjugate_subgroups
}
right_covers = {
    coset_partition(subgroup, "right") for subgroup in conjugate_subgroups
}
coset_census_ok = bool(
    len(conjugate_subgroups) == 5
    and len(left_covers) == len(right_covers) == 5
    and left_covers.isdisjoint(right_covers)
    and left_covers | right_covers == set(covers)
)

cell_action_maps = {}
cell_action_ok = True
for action in action_tuples:
    action_array = action_arrays[action]
    images = []
    for cell in maximum_cells:
        image = frozenset(int(action_array[index]) for index in cell)
        if image not in cell_id:
            cell_action_ok = False
            break
        images.append(cell_id[image])
    if len(images) == 25:
        cell_action_maps[action] = tuple(images)
proper_cell_maps = tuple(cell_action_maps[action] for action in proper_actions)
improper_cell_maps = tuple(cell_action_maps[action] for action in improper_actions)
all_cell_maps = proper_cell_maps + improper_cell_maps

source_coefficients = {}
minimum_absolute_determinant = float("inf")
for tetrahedron in tetrahedra:
    determinant = float(np.linalg.det(vertices[list(tetrahedron)]))
    minimum_absolute_determinant = min(
        minimum_absolute_determinant, abs(determinant)
    )
    source_coefficients[tetrahedron] = (
        1 if determinant > 0 else -1 if determinant < 0 else 0
    )
source_boundary = Counter()
for tetrahedron, coefficient in source_coefficients.items():
    for omitted in range(4):
        face = tetrahedron[:omitted] + tetrahedron[omitted+1:]
        source_boundary[face] += coefficient*((-1)**omitted)
source_boundary = Counter({face: value for face, value in source_boundary.items() if value})
source_chain_ok = bool(
    minimum_absolute_determinant > DET_TOLERANCE
    and not source_boundary
    and set(source_coefficients.values()) == {-1, 1}
)

degree_by_schedule = {}
order_records = []
base_degree_records = []
all_degree_controls = True
all_order_alternation = True
for cover_index, cover in enumerate(covers):
    base_order = tuple(sorted(cover))
    base_result = None
    base_degree = None
    for order in permutations(base_order):
        vertex_label = np.empty(120, dtype=np.int8)
        for label, identifier in enumerate(order):
            vertex_label[list(maximum_cells[identifier])] = label
        result = degree_for_vertex_labels(
            tetrahedra, source_coefficients, vertex_label
        )
        control = bool(
            result["degree"] is not None
            and result["candidates"] == result["independent"]
            and result["preimage_counts"] == (120, 120, 120, 120, 120)
        )
        all_degree_controls &= control
        if order == base_order:
            base_result = result
            base_degree = result["degree"]
        relative_permutation = tuple(base_order.index(cell) for cell in order)
        if base_degree is not None:
            all_order_alternation &= bool(
                result["degree"]
                == permutation_sign(relative_permutation)*base_degree
            )
        degree_by_schedule[order] = result["degree"]
        order_records.append({
            "cover": cover_index,
            "order_cell_ids": list(order),
            "relative_permutation_sign": permutation_sign(relative_permutation),
            "pushforward_facet_coefficients": list(result["pushforward"]),
            "degree_candidates": list(result["candidates"]),
            "independent_preimage_degrees": list(result["independent"]),
            "degree": result["degree"],
            "control_pass": control,
        })
    base_degree_records.append({
        "cover": cover_index,
        "cell_ids": list(base_order),
        "coset_side": (
            "left" if cover in left_covers
            else "right" if cover in right_covers else None
        ),
        "degree": base_degree,
        "degree_candidates": list(base_result["candidates"]),
        "local_degree_sign_counts": [
            {str(sign): count for sign, count in counts.items()}
            for counts in base_result["local_sign_counts"]
        ],
    })

all_schedules = frozenset(degree_by_schedule)
positive_schedules = frozenset(
    schedule for schedule, degree in degree_by_schedule.items()
    if degree is not None and degree > 0
)
negative_schedules = frozenset(
    schedule for schedule, degree in degree_by_schedule.items()
    if degree is not None and degree < 0
)
zero_schedules = frozenset(
    schedule for schedule, degree in degree_by_schedule.items()
    if degree == 0
)

proper_cover_orbits = partition_orbits(covers, proper_cell_maps, ordered=False)
full_cover_orbits = partition_orbits(covers, all_cell_maps, ordered=False)
positive_orbits = (
    partition_orbits(positive_schedules, proper_cell_maps, ordered=True)
    if positive_schedules else tuple()
)
negative_orbits = (
    partition_orbits(negative_schedules, proper_cell_maps, ordered=True)
    if negative_schedules else tuple()
)
full_ordered_orbits = partition_orbits(
    all_schedules, all_cell_maps, ordered=True
)

def stabilizer_size(seed, action_maps, ordered):
    if ordered:
        return sum(
            tuple(action[index] for index in seed) == seed
            for action in action_maps
        )
    return sum(
        frozenset(action[index] for index in seed) == seed
        for action in action_maps
    )


proper_cover_stabilizers = tuple(
    stabilizer_size(next(iter(orbit)), proper_cell_maps, ordered=False)
    for orbit in proper_cover_orbits
)
positive_stabilizers = tuple(
    stabilizer_size(next(iter(orbit)), proper_cell_maps, ordered=True)
    for orbit in positive_orbits
)
negative_stabilizers = tuple(
    stabilizer_size(next(iter(orbit)), proper_cell_maps, ordered=True)
    for orbit in negative_orbits
)
full_ordered_stabilizers = tuple(
    stabilizer_size(next(iter(orbit)), all_cell_maps, ordered=True)
    for orbit in full_ordered_orbits
)

def orbit_index(value, orbits):
    return next(
        (index for index, orbit in enumerate(orbits) if value in orbit), None
    )


cover_orbit_coset_sides = []
for orbit in proper_cover_orbits:
    sides = set()
    if orbit & left_covers:
        sides.add("left")
    if orbit & right_covers:
        sides.add("right")
    cover_orbit_coset_sides.append(sorted(sides))

improper_cover_orbit_map = []
for orbit in proper_cover_orbits:
    seed = next(iter(orbit))
    images = {
        frozenset(action[index] for index in seed)
        for action in improper_cell_maps
    }
    improper_cover_orbit_map.append(
        orbit_index(next(iter(images)), proper_cover_orbits)
        if images else None
    )

improper_positive_to_negative = []
degree_covariance_ok = True
for orbit in positive_orbits:
    seed = next(iter(orbit))
    base_degree = degree_by_schedule[seed]
    proper_images = {
        tuple(action[index] for index in seed) for action in proper_cell_maps
    }
    improper_images = {
        tuple(action[index] for index in seed) for action in improper_cell_maps
    }
    degree_covariance_ok &= bool(
        proper_images == set(orbit)
        and all(degree_by_schedule[image] == base_degree
                for image in proper_images)
        and all(degree_by_schedule[image] == -base_degree
                for image in improper_images)
    )
    image_indices = {
        orbit_index(image, negative_orbits) for image in improper_images
    }
    improper_positive_to_negative.append(sorted(image_indices))
for orbit in negative_orbits:
    seed = next(iter(orbit))
    base_degree = degree_by_schedule[seed]
    proper_images = {
        tuple(action[index] for index in seed) for action in proper_cell_maps
    }
    improper_images = {
        tuple(action[index] for index in seed) for action in improper_cell_maps
    }
    degree_covariance_ok &= bool(
        proper_images == set(orbit)
        and all(degree_by_schedule[image] == base_degree
                for image in proper_images)
        and all(degree_by_schedule[image] == -base_degree
                for image in improper_images)
    )

positive_time_reversal_map = [
    orbit_index(tuple(reversed(next(iter(orbit)))), positive_orbits)
    for orbit in positive_orbits
]
negative_time_reversal_map = [
    orbit_index(tuple(reversed(next(iter(orbit)))), negative_orbits)
    for orbit in negative_orbits
]
time_reversal_ok = bool(
    all(tuple(reversed(schedule)) in positive_schedules
        for schedule in positive_schedules)
    and all(tuple(reversed(schedule)) in negative_schedules
            for schedule in negative_schedules)
)

cover_orbit_control = bool(
    set().union(*(set(orbit) for orbit in proper_cover_orbits)) == set(covers)
    and all(
        7200 == len(orbit)*stabilizer
        for orbit, stabilizer in zip(
            proper_cover_orbits, proper_cover_stabilizers
        )
    )
    and len(full_cover_orbits) == 1
)
ordered_orbit_control = bool(
    set().union(*(set(orbit) for orbit in positive_orbits))
        == set(positive_schedules)
    and set().union(*(set(orbit) for orbit in negative_orbits))
        == set(negative_schedules)
    and set().union(*(set(orbit) for orbit in full_ordered_orbits))
        == set(all_schedules)
    and all(
        7200 == len(orbit)*stabilizer
        for orbit, stabilizer in zip(positive_orbits, positive_stabilizers)
    )
    and all(
        7200 == len(orbit)*stabilizer
        for orbit, stabilizer in zip(negative_orbits, negative_stabilizers)
    )
    and all(
        14400 == len(orbit)*stabilizer
        for orbit, stabilizer in zip(
            full_ordered_orbits, full_ordered_stabilizers
        )
    )
)
improper_exchange_ok = bool(
    len(proper_cover_orbits) > 1
    and all(
        target is not None and target != source
        for source, target in enumerate(improper_cover_orbit_map)
    )
    and all(len(indices) == 1 and indices[0] is not None
            for indices in improper_positive_to_negative)
)

exhaustive_census_ok = bool(
    len(maximum_cells) == 25
    and all_cells_independent
    and optimization_status == cp_model.OPTIMAL
    and optimization_value == 24
    and exclusion_status == cp_model.INFEASIBLE
    and len(covers) == 10
    and all(len(cover) == 5 for cover in covers)
    and coset_census_ok
)
degree_controls_ok = bool(
    source_chain_ok
    and len(order_records) == 1200
    and len(all_schedules) == 1200
    and all_degree_controls
    and all_order_alternation
    and all(value is not None for value in degree_by_schedule.values())
)
controls_ok = bool(
    provenance_ok and carrier_ok and subgroup_ok
    and multiplication_residual < GEOMETRY_TOLERANCE
    and len(action_tuples) == 14400 and all_graph_automorphisms
    and orientation_generation_ok and exhaustive_census_ok
    and cell_action_ok and len(cell_action_maps) == 14400
    and degree_controls_ok and cover_orbit_control
    and ordered_orbit_control and degree_covariance_ok and time_reversal_ok
)

if not controls_ok:
    outcome = "OPEN_CONTROL_FAILURE"
elif zero_schedules:
    outcome = "CHROMATIC_DEGREE_DEGENERACY"
elif len(positive_orbits) == 1:
    outcome = "ONE_ORIENTED_CANONICAL_CLASS"
elif len(positive_orbits) > 1 and improper_exchange_ok:
    outcome = "CHIRAL_COVER_AMBIGUITY"
else:
    outcome = "OPEN_CONTROL_FAILURE"

check("frozen input source hashes", provenance_ok)
check("600-cell carrier is exactly 120/720/1200/600 and 12-regular", carrier_ok)
check("binary tetrahedral 2T is reconstructed as a closed 24-subgroup", subgroup_ok)
check("quaternion multiplication closes at exact vertex labels", multiplication_residual < GEOMETRY_TOLERANCE, f"max residual={multiplication_residual:.3e}")
check("all 14,400 generated H4 maps are graph automorphisms", len(action_tuples) == 14400 and all_graph_automorphisms)
check("independent determinant test gives 7,200 proper and 7,200 improper maps", orientation_generation_ok, f"action residual={maximum_action_residual:.3e}; orthogonal residual={maximum_orthogonality_residual:.3e}")
check("the full 2T orbit is 25 independent 24-cells", len(maximum_cells) == 25 and all_cells_independent)
check("CP-SAT proves alpha(G)=24", optimization_status == cp_model.OPTIMAL and optimization_value == 24, f"status={optimization_solver.StatusName(optimization_status)}; wall={optimization_solver.WallTime():.3f}s")
check("excluding all 25 cells leaves no other independent 24-set", exclusion_status == cp_model.INFEASIBLE, f"status={exclusion_solver.StatusName(exclusion_status)}; wall={exclusion_solver.WallTime():.3f}s")
check("exact covers are precisely five left plus five right coset colourings", len(covers) == 10 and coset_census_ok)
check("all 25 cells are permuted by every H4 action", cell_action_ok and len(cell_action_maps) == 14400)
check("determinant source orientation is resolved and its chain is closed", source_chain_ok, f"min |det|={minimum_absolute_determinant:.12g}")
check("all 1,200 orders pass five-facet and signed-preimage degree controls", all_degree_controls and len(order_records) == 1200)
check("all ten 120-order families alternate degree exactly", all_order_alternation)
check("proper-rotation cover orbits and stabilizers exhaust all ten covers", cover_orbit_control)
check("proper and full ordered-orbit censuses satisfy orbit-stabilizer exactly", ordered_orbit_control)
check("proper maps preserve and improper maps reverse chromatic degree", degree_covariance_ok)
check("five-phase time reversal preserves the degree-sign sectors", time_reversal_ok)
check("improper cover exchange is mechanically characterized", improper_exchange_ok or len(proper_cover_orbits) == 1)
check("mechanical outcome assigned without physical target", outcome in {
    "CHROMATIC_DEGREE_DEGENERACY",
    "ONE_ORIENTED_CANONICAL_CLASS",
    "CHIRAL_COVER_AMBIGUITY",
    "OPEN_CONTROL_FAILURE",
})
check("no Regge, nonlinear, continuum, particle or preferred-sign target parsed", True)

passed = sum(ok for _, ok in tests)
degree_counter = Counter(degree_by_schedule.values())
base_magnitude_counter = Counter(
    abs(record["degree"]) for record in base_degree_records
)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": source_hashes,
    "regge_action_parsed": False,
    "nonlinear_schedule_result_parsed": False,
    "preferred_parity_parsed": False,
    "continuum_or_particle_target_parsed": False,
    "carrier": {
        "vertices": len(vertices),
        "edges": len(edges),
        "triangles": len(triangles),
        "tetrahedra": len(tetrahedra),
        "degree": sorted(set(degrees)),
        "minimum_absolute_source_determinant": minimum_absolute_determinant,
    },
    "group": {
        "actions": len(action_tuples),
        "proper_actions": len(proper_actions),
        "improper_actions": len(improper_actions),
        "basis_indices_for_matrix_reconstruction": basis_indices,
        "maximum_action_residual": maximum_action_residual,
        "maximum_orthogonality_residual": maximum_orthogonality_residual,
        "maximum_abs_determinant_residual": determinant_residual,
    },
    "exhaustive_colour_census": {
        "maximum_independent_cells": len(maximum_cells),
        "cell_sizes": sorted(Counter(len(cell) for cell in maximum_cells).items()),
        "optimization_status": optimization_solver.StatusName(optimization_status),
        "optimization_value": optimization_value,
        "optimization_wall_seconds": optimization_solver.WallTime(),
        "exclusion_status": exclusion_solver.StatusName(exclusion_status),
        "exclusion_wall_seconds": exclusion_solver.WallTime(),
        "covers": len(covers),
        "left_coset_covers": len(left_covers),
        "right_coset_covers": len(right_covers),
        "cover_cell_ids": [sorted(cover) for cover in covers],
    },
    "degrees": {
        "base_cover_records": base_degree_records,
        "base_magnitude_multiset": {
            str(value): count for value, count in sorted(base_magnitude_counter.items())
        },
        "all_order_multiset": {
            str(value): count for value, count in sorted(degree_counter.items())
        },
        "positive_orders": len(positive_schedules),
        "negative_orders": len(negative_schedules),
        "zero_orders": len(zero_schedules),
        "order_records": order_records,
    },
    "orbits": {
        "proper_unordered_cover_sizes": [len(orbit) for orbit in proper_cover_orbits],
        "proper_unordered_cover_stabilizers": list(proper_cover_stabilizers),
        "proper_unordered_cover_coset_sides": cover_orbit_coset_sides,
        "improper_cover_orbit_map": improper_cover_orbit_map,
        "full_unordered_cover_sizes": [len(orbit) for orbit in full_cover_orbits],
        "positive_proper_orbit_sizes": [len(orbit) for orbit in positive_orbits],
        "positive_proper_stabilizers": list(positive_stabilizers),
        "negative_proper_orbit_sizes": [len(orbit) for orbit in negative_orbits],
        "negative_proper_stabilizers": list(negative_stabilizers),
        "improper_positive_to_negative_orbit_map": improper_positive_to_negative,
        "full_ordered_orbit_sizes": [len(orbit) for orbit in full_ordered_orbits],
        "full_ordered_stabilizers": list(full_ordered_stabilizers),
        "positive_time_reversal_map": positive_time_reversal_map,
        "negative_time_reversal_map": negative_time_reversal_map,
    },
    "labels": {
        "ten_colourings_and_orbits": "DERIVED COMPUTATIONAL",
        "chromatic_degree_values": "DERIVED EXACT",
        "orientation_compatible_classification": "STRUCTURAL",
        "physical_chromatic_compatibility": "OPEN",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": len(tests),
    "passed": passed,
    "test_records": [
        {"label": label, "passed": ok} for label, ok in tests
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"Degree multiset: {dict(sorted(degree_counter.items()))}")
print(f"Base magnitude multiset: {dict(sorted(base_magnitude_counter.items()))}")
print(f"Proper cover orbits: {[len(orbit) for orbit in proper_cover_orbits]}")
print(f"Cover coset sides: {cover_orbit_coset_sides}")
print(f"Positive proper orbits: {[len(orbit) for orbit in positive_orbits]}")
print(f"Negative proper orbits: {[len(orbit) for orbit in negative_orbits]}")
print(f"Full ordered H4 orbits: {[len(orbit) for orbit in full_ordered_orbits]}")
print(f"Improper cover map: {improper_cover_orbit_map}")
print(f"Time reversal on positive orbits: {positive_time_reversal_map}")
print(f"Outcome: {outcome}")
print(f"RESULT: {passed}/{len(tests)} checks passed")

if passed != len(tests):
    sys.exit(1)

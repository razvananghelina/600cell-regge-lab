#!/usr/bin/env python3
"""Exact finite audit of phased global tent schedules on the 600-cell.

Protocol commit: af50841.  The values 24, 25, 10 and the preliminary slab
f-vector were disclosed as exploratory observations before that commit.  This
script reproduces them, exhausts all 10*5! ordered schedules, and keeps the
combinatorial/kinematic conclusion separate from Regge stationarity.
"""

from collections import Counter
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ortools.sat.python import cp_model
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_global_tent_schedule.json"
PROTOCOL_COMMIT = "af50841"
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


def canonical_set_key(value):
    return tuple(sorted(value))


def canonical_partition_key(value):
    return tuple(sorted(canonical_set_key(cell) for cell in value))


def permutation_parity(value):
    inversions = sum(
        value[i] > value[j]
        for i in range(len(value)) for j in range(i + 1, len(value))
    )
    return inversions % 2


def build_tetrahedra(adjacency):
    neighbours = [set(np.flatnonzero(row).tolist()) for row in adjacency]
    edges = []
    triangles = []
    tetrahedra = []
    for a in range(120):
        for b in sorted(vertex for vertex in neighbours[a] if vertex > a):
            edges.append((a, b))
            for c in sorted(
                vertex for vertex in neighbours[a] & neighbours[b]
                if vertex > b
            ):
                triangles.append((a, b, c))
                for d in sorted(
                    vertex
                    for vertex in neighbours[a] & neighbours[b] & neighbours[c]
                    if vertex > c
                ):
                    tetrahedra.append((a, b, c, d))
    return neighbours, tuple(edges), tuple(triangles), tuple(tetrahedra)


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
    return model, variables


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


def build_slab(tetrahedra, phase):
    four_simplices = set()
    missing_colors = Counter()
    for tetrahedron in tetrahedra:
        ordered = sorted(tetrahedron, key=phase.__getitem__)
        colors = {phase[vertex] for vertex in tetrahedron}
        if len(colors) != 4:
            raise RuntimeError("a spatial tetrahedron repeats a phase color")
        missing_colors[(set(range(5)) - colors).pop()] += 1
        for vertex in ordered:
            simplex = [vertex, vertex + 120]
            simplex.extend(
                other + 120 if phase[other] < phase[vertex] else other
                for other in ordered if other != vertex
            )
            if len(set(simplex)) != 5:
                raise RuntimeError("degenerate staircase 4-simplex")
            four_simplices.add(tuple(sorted(simplex)))
    return four_simplices, missing_colors


print("=" * 78)
print("GLOBAL FIVE-PHASE 600-CELL TENT SCHEDULE")
print("=" * 78)

vertices, adjacency_float, _ = build_600cell()
adjacency = adjacency_float > 0.5
neighbours, edges, triangles, tetrahedra = build_tetrahedra(adjacency)

check(
    "the spatial carrier is the (120,720,1200,600) 600-cell boundary",
    tuple(map(len, (vertices, edges, triangles, tetrahedra)))
    == (120, 720, 1200, 600),
)
check(
    "the chamber graph is 12-regular",
    {len(row) for row in neighbours} == {12},
)

# Local tent carriers share a 600-cell tetrahedron exactly for adjacent poles.
vertex_stars = [set() for _ in range(120)]
for tetrahedron_id, tetrahedron in enumerate(tetrahedra):
    for vertex in tetrahedron:
        vertex_stars[vertex].add(tetrahedron_id)
conflict = np.zeros((120, 120), dtype=bool)
for left in range(120):
    for right in range(left + 1, 120):
        conflict[left, right] = conflict[right, left] = bool(
            vertex_stars[left] & vertex_stars[right]
        )
check(
    "tent-move conflict is exactly 600-cell edge adjacency",
    np.array_equal(conflict, adjacency),
)

# Quaternion multiplication table and the complete order-14400 H4 action.
multiplication = np.empty((120, 120), dtype=np.int16)
closure_residual = 0.0
for left in range(120):
    for right in range(120):
        product = qmul(vertices[left], vertices[right])
        image = int(np.argmax(vertices @ product))
        multiplication[left, right] = image
        closure_residual = max(
            closure_residual, float(np.linalg.norm(product - vertices[image]))
        )

conjugate = np.empty(120, dtype=np.int16)
negative = np.empty(120, dtype=np.int16)
for index, vertex in enumerate(vertices):
    conjugated = vertex * np.array((1.0, -1.0, -1.0, -1.0))
    conjugate[index] = int(np.argmax(vertices @ conjugated))
    negative[index] = int(np.argmax(vertices @ (-vertex)))

identity = int(np.argmax(vertices[:, 0]))
check(
    "the 120 vertices close as unit quaternions with controlled residual",
    closure_residual < 5e-9
    and np.all(multiplication[identity] == np.arange(120))
    and np.all(multiplication[:, identity] == np.arange(120)),
    f"max closure residual={closure_residual:.3e}",
)
check(
    "quaternion conjugation supplies every group inverse",
    all(
        multiplication[index, conjugate[index]] == identity
        and multiplication[conjugate[index], index] == identity
        for index in range(120)
    ),
)

h4_permutations = set()
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
            h4_permutations.add(action)

h4_arrays = tuple(np.asarray(action, dtype=np.int16) for action in h4_permutations)
edge_array = np.asarray(edges, dtype=np.int16)
all_automorphisms = all(
    len(set(action.tolist())) == 120
    and np.all(adjacency[action[edge_array[:, 0]], action[edge_array[:, 1]]])
    for action in h4_arrays
)
check(
    "left-right quaternion maps construct 14400 distinct H4 graph automorphisms",
    len(h4_arrays) == 14400 and all_automorphisms,
)

# The standard Hurwitz units: Q8 axes plus all 16 half-coordinate units.
binary_tetrahedral = frozenset(
    index for index, vertex in enumerate(vertices)
    if (
        np.count_nonzero(np.abs(vertex) > 1e-8) == 1
        and np.max(np.abs(vertex)) > 1 - 1e-8
    ) or np.all(np.abs(np.abs(vertex) - 0.5) < 1e-8)
)
subgroup_closed = all(
    int(multiplication[left, right]) in binary_tetrahedral
    for left in binary_tetrahedral for right in binary_tetrahedral
)
subgroup_inverse = all(conjugate[index] in binary_tetrahedral
                       for index in binary_tetrahedral)
check(
    "the 24 Hurwitz units form the binary tetrahedral subgroup 2T",
    len(binary_tetrahedral) == 24
    and identity in binary_tetrahedral
    and subgroup_closed and subgroup_inverse,
)
check(
    "2T is an independent 24-set in the tent-conflict graph",
    not np.any(adjacency[np.ix_(
        tuple(binary_tetrahedral), tuple(binary_tetrahedral)
    )]),
)

# Orbit census of maximum independent sets.
maximal_candidates = {
    frozenset(int(action[index]) for index in binary_tetrahedral)
    for action in h4_arrays
}
all_orbit_sets_independent = all(
    len(cell) == 24
    and not np.any(adjacency[np.ix_(tuple(cell), tuple(cell))])
    for cell in maximal_candidates
)
check(
    "the H4 orbit of 2T contains 25 independent 24-sets",
    len(maximal_candidates) == 25 and all_orbit_sets_independent,
)

independent_optimization, optimization_variables = independent_model(edges)
optimization_solver, optimization_status = solve_exact(independent_optimization)
optimization_value = (
    int(round(optimization_solver.ObjectiveValue()))
    if optimization_status == cp_model.OPTIMAL else None
)
check(
    "CP-SAT proves the independence number alpha(G)=24",
    optimization_status == cp_model.OPTIMAL and optimization_value == 24,
    f"status={optimization_solver.StatusName(optimization_status)}, "
    f"objective={optimization_value}, wall={optimization_solver.WallTime():.3f}s",
)

exclusion_model, _ = independent_model(
    edges, target=24, excluded=maximal_candidates
)
exclusion_solver, exclusion_status = solve_exact(exclusion_model)
check(
    "excluding the H4 orbit leaves no further independent 24-set",
    exclusion_status == cp_model.INFEASIBLE,
    f"status={exclusion_solver.StatusName(exclusion_status)}, "
    f"wall={exclusion_solver.WallTime():.3f}s",
)

dot_values = tuple(sorted(set(np.round((vertices @ vertices.T).ravel(), 8))))
largest_nonedge_dot = max(
    float(vertices[left] @ vertices[right])
    for left in range(120) for right in range(left + 1, 120)
    if not adjacency[left, right]
)
check(
    "the alpha bound is independently recognizable as the 4D kissing bound",
    np.isclose(largest_nonedge_dot, 0.5, atol=1e-8)
    and len(dot_values) == 9,
    "an independent set is a spherical code with pairwise dot <=1/2; "
    "Musin's theorem k(4)=24 supplies an external analytic upper bound",
)

intersection_histogram = Counter(
    len(left & right)
    for position, left in enumerate(maximal_candidates)
    for right in tuple(maximal_candidates)[position + 1:]
)
membership_counts = Counter(
    sum(vertex in cell for cell in maximal_candidates)
    for vertex in range(120)
)
check(
    "the 25 maximum sets form the (0,6)-intersection design",
    intersection_histogram == Counter({0: 100, 6: 200})
    and membership_counts == Counter({5: 120}),
    f"intersections={dict(sorted(intersection_histogram.items()))}",
)

covers = exact_covers(tuple(maximal_candidates))
cover_membership = Counter(
    sum(cell in cover for cover in covers) for cell in maximal_candidates
)
check(
    "the 25 maximum sets have exactly ten five-cell exact covers",
    len(covers) == 10
    and all(len(cover) == 5 for cover in covers)
    and cover_membership == Counter({2: 25}),
)

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
    return frozenset(cells)


conjugate_subgroups = {
    frozenset(
        int(multiplication[multiplication[group_element, element],
                           conjugate[group_element]])
        for element in binary_tetrahedral
    )
    for group_element in range(120)
}
coset_covers = {
    coset_partition(subgroup, side)
    for subgroup in conjugate_subgroups for side in ("left", "right")
}
check(
    "all ten covers are precisely left/right coset covers of five conjugate 2T subgroups",
    len(conjugate_subgroups) == 5
    and len(coset_covers) == 10
    and coset_covers == covers,
)

# Symmetry orbit and phase-order census.
ordered_covers = tuple(sorted(covers, key=canonical_partition_key))
base_cover = ordered_covers[0]
cover_orbit = set()
cover_stabilizer = []
for action in h4_arrays:
    image = frozenset(
        frozenset(int(action[index]) for index in cell)
        for cell in base_cover
    )
    cover_orbit.add(image)
    if image == base_cover:
        cover_stabilizer.append(action)

base_cells = tuple(sorted(base_cover, key=canonical_set_key))
base_cell_position = {cell: position for position, cell in enumerate(base_cells)}
induced_color_actions = set()
for action in cover_stabilizer:
    induced_color_actions.add(tuple(
        base_cell_position[frozenset(int(action[index]) for index in cell)]
        for cell in base_cells
    ))
all_even_color_actions = {
    action for action in permutations(range(5))
    if permutation_parity(action) == 0
}
pointwise_color_stabilizer = sum(
    all(
        frozenset(int(action[index]) for index in cell) == cell
        for cell in base_cells
    )
    for action in cover_stabilizer
)
check(
    "the ten unordered covers form one H4 orbit with stabilizer 1440",
    cover_orbit == covers and len(cover_stabilizer) == 1440,
)
check(
    "a cover stabilizer induces exactly A5 on its five colors",
    induced_color_actions == all_even_color_actions
    and len(induced_color_actions) == 60
    and pointwise_color_stabilizer == 24,
    f"induced group={len(induced_color_actions)}, kernel={pointwise_color_stabilizer}",
)

base_order_even = base_cells
base_order_odd = (
    base_cells[1], base_cells[0], base_cells[2], base_cells[3], base_cells[4]
)

def ordered_orbit(seed):
    return {
        tuple(
            frozenset(int(action[index]) for index in cell)
            for cell in seed
        )
        for action in h4_arrays
    }


even_schedule_orbit = ordered_orbit(base_order_even)
odd_schedule_orbit = ordered_orbit(base_order_odd)
all_ordered_schedules = {
    ordering
    for cover in covers
    for ordering in permutations(tuple(cover))
}
check(
    "the 1200 ordered schedules split into two H4 parity orbits of 600",
    len(all_ordered_schedules) == 1200
    and len(even_schedule_orbit) == len(odd_schedule_orbit) == 600
    and even_schedule_orbit.isdisjoint(odd_schedule_orbit)
    and even_schedule_orbit | odd_schedule_orbit == all_ordered_schedules,
)
check(
    "time reversal preserves, rather than exchanges, both parity orbits",
    tuple(reversed(base_order_even)) in even_schedule_orbit
    and tuple(reversed(base_order_odd)) in odd_schedule_orbit,
)

# Exhaust every cover and every phase ordering.
old_boundary = set(tetrahedra)
new_boundary = {
    tuple(vertex + 120 for vertex in tetrahedron) for tetrahedron in tetrahedra
}
expected_boundary = old_boundary | new_boundary
expected_f_vector = (240, 2280, 6240, 6600, 2400)
expected_edge_types = Counter({
    "old": 720, "new": 720, "diagonal": 720, "pole": 120,
})
schedule_count = 0
all_schedule_checks = True
schedule_summaries = set()
slab_digests = set()
representative_slab = None
representative_by_parity = {}

for cover in ordered_covers:
    cells = tuple(sorted(cover, key=canonical_set_key))
    for ordering in permutations(cells):
        schedule_count += 1
        phase = {
            vertex: phase_index
            for phase_index, cell in enumerate(ordering)
            for vertex in cell
        }
        phases_independent = all(
            len(cell) == 24
            and not np.any(adjacency[np.ix_(tuple(cell), tuple(cell))])
            for cell in ordering
        )
        slab, missing_colors = build_slab(tetrahedra, phase)
        if representative_slab is None:
            representative_slab = slab
        parity = 0 if ordering in even_schedule_orbit else 1
        if parity not in representative_by_parity:
            representative_by_parity[parity] = slab

        # These checks are performed on every one of the 1200 schedules.
        # The more allocation-heavy face census below is performed once in
        # each of the two already proved H4 orbits.  Since H4 acts separately
        # on old and new vertex copies and the staircase construction is
        # equivariant, this is an exact symmetry reduction, not sampling.
        local_ok = (
            phases_independent
            and len(slab) == 2400
            and missing_colors == Counter({color: 120 for color in range(5)})
        )
        all_schedule_checks &= local_ok
        simplex_array = np.asarray(sorted(slab), dtype=np.uint16)
        slab_digests.add(hashlib.sha256(simplex_array.tobytes()).hexdigest())

for parity, slab in representative_by_parity.items():
    faces_by_size = tuple(
        {
            tuple(sorted(face))
            for simplex in slab for face in combinations(simplex, size)
        }
        for size in range(1, 6)
    )
    f_vector = tuple(len(layer) for layer in faces_by_size)
    euler = sum((-1)**dimension * count
                for dimension, count in enumerate(f_vector))
    facet_multiplicity = Counter(
        tuple(sorted(facet))
        for simplex in slab for facet in combinations(simplex, 4)
    )
    multiplicity_histogram = Counter(facet_multiplicity.values())
    boundary = {
        facet for facet, multiplicity in facet_multiplicity.items()
        if multiplicity == 1
    }
    edge_types = Counter()
    for left, right in faces_by_size[1]:
        if right < 120:
            edge_types["old"] += 1
        elif left >= 120:
            edge_types["new"] += 1
        elif right - left == 120:
            edge_types["pole"] += 1
        else:
            edge_types["diagonal"] += 1
    orbit_representative_ok = (
        f_vector == expected_f_vector
        and euler == 0
        and multiplicity_histogram == Counter({2: 5400, 1: 1200})
        and boundary == expected_boundary
        and edge_types == expected_edge_types
    )
    all_schedule_checks &= orbit_representative_ok
    schedule_summaries.add((
        len(slab), f_vector, euler,
        tuple(sorted(multiplicity_histogram.items())),
        tuple(sorted(edge_types.items())),
    ))

check(
    "all 10*5!=1200 phase schedules are nonoverlapping and balanced",
    schedule_count == 1200 and all_schedule_checks,
)
check(
    "both exhaustive H4 order orbits give the same product slab and exact boundaries",
    len(schedule_summaries) == 1
    and len(representative_by_parity) == 2
    and next(iter(schedule_summaries))[1] == expected_f_vector,
    f"f={expected_f_vector}, chi=0, boundary=600+600",
)
check(
    "the 1200 labelled phase orders give 1200 distinct staircase triangulations",
    len(slab_digests) == 1200,
)

# Exact Lorentzian metric family, plus a direct check on all 2400 simplices.
rho = sp.symbols("rho", positive=True)
link_gram = sp.Matrix([
    [1, sp.Rational(1, 2), sp.Rational(1, 2)],
    [sp.Rational(1, 2), 1, sp.Rational(1, 2)],
    [sp.Rational(1, 2), sp.Rational(1, 2), 1],
])
pole_link = sp.Matrix([[-rho/2, -rho/2, -rho/2]])
regular_gram = pole_link.col_join(link_gram)
regular_gram = sp.Matrix.hstack(sp.Matrix([[-rho], [-rho/2], [-rho/2], [-rho/2]]),
                                regular_gram)
regular_determinant = sp.factor(regular_gram.det())
schur_complement = sp.factor(
    -rho - (pole_link * link_gram.inv() * pole_link.T)[0]
)
check(
    "the regular slab simplex has Lorentzian inertia for every rho>0",
    regular_determinant == -rho*(3*rho + 8)/16
    and schur_complement == -rho*(3*rho + 8)/8
    and all(value > 0 for value in link_gram.eigenvals()),
    f"det G={regular_determinant}; Schur={schur_complement}",
)

rho_numeric = 0.25
expected_determinant = -rho_numeric*(8 + 3*rho_numeric)/16
max_determinant_residual = 0.0
inertia_histogram = Counter()
pole_pair_histogram = Counter()
for simplex in representative_slab:
    simplex = tuple(simplex)
    pole_pairs = [
        (left, right)
        for left, right in combinations(simplex, 2)
        if abs(left - right) == 120
    ]
    pole_pair_histogram[len(pole_pairs)] += 1
    squared = np.ones((5, 5))
    np.fill_diagonal(squared, 0.0)
    for left_position, right_position in combinations(range(5), 2):
        if abs(simplex[left_position] - simplex[right_position]) == 120:
            squared[left_position, right_position] = -rho_numeric
            squared[right_position, left_position] = -rho_numeric
    gram = np.empty((4, 4))
    for left in range(1, 5):
        for right in range(1, 5):
            gram[left-1, right-1] = (
                squared[0, left] + squared[0, right] - squared[left, right]
            ) / 2
    eigenvalues = np.linalg.eigvalsh(gram)
    inertia_histogram[(
        int(np.sum(eigenvalues < -1e-10)),
        int(np.sum(eigenvalues > 1e-10)),
    )] += 1
    max_determinant_residual = max(
        max_determinant_residual,
        abs(float(np.linalg.det(gram)) - expected_determinant),
    )

check(
    "all 2400 staircase 4-simplices carry the same coherent Lorentzian metric",
    pole_pair_histogram == Counter({1: 2400})
    and inertia_histogram == Counter({(1, 3): 2400})
    and max_determinant_residual < 1e-12,
    f"inertia={dict(inertia_histogram)}, max det residual={max_determinant_residual:.3e}",
)

regular_cosine = (2 + rho_numeric)/(2*(3 + rho_numeric))
regular_deficit = 2*np.pi - 5*np.arccos(regular_cosine)
check(
    "the coherent regular Lorentzian slab is not a stationary vacuum tick",
    regular_deficit > 0.1,
    f"each regular pole has positive deficit {regular_deficit:.12f} rad",
)

chi_lower_bound = int(np.ceil(120/optimization_value))
chromatic_number = 5 if chi_lower_bound == 5 and len(covers) > 0 else None
check(
    "alpha=24 and an exact cover prove the chromatic number chi(G)=5",
    chromatic_number == 5,
    "five phases are forced for simultaneous nonoverlapping local moves",
)

verdict = (
    "DERIVED COMPUTATIONAL: the 600-cell tent-conflict graph has chromatic "
    "number five; its 25 maximum independent 24-sets form ten H4-equivalent "
    "unordered schedules.  All 1200 ordered schedules triangulate a coherent "
    "Lorentzian S3 product slab, but split into two H4 parity orbits and the "
    "regular metric is not Regge-stationary."
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "provenance": "frozen after disclosed exploratory observations",
    "spatial_f_vector": [120, 720, 1200, 600],
    "tent_conflict_degree": 12,
    "h4_action_order": len(h4_arrays),
    "binary_tetrahedral_order": len(binary_tetrahedral),
    "independence_number": optimization_value,
    "maximum_independent_sets": len(maximal_candidates),
    "maximum_set_intersections": dict(sorted(intersection_histogram.items())),
    "maximum_sets_per_vertex": 5,
    "unordered_exact_covers": len(covers),
    "covers_per_maximum_set": 2,
    "cover_h4_orbits": 1,
    "cover_stabilizer_order": len(cover_stabilizer),
    "induced_color_group": "A5",
    "induced_color_group_order": len(induced_color_actions),
    "ordered_schedules": schedule_count,
    "ordered_schedule_h4_orbits": [600, 600],
    "time_reversal_exchanges_ordered_orbits": False,
    "chromatic_number": chromatic_number,
    "global_slab_f_vector": list(expected_f_vector),
    "global_slab_euler_characteristic": 0,
    "global_slab_boundary_tetrahedra": [600, 600],
    "global_slab_edge_types": dict(expected_edge_types),
    "global_slab_four_simplices": 2400,
    "lorentzian_metric_family": {
        "spacelike_squared_length": "a^2",
        "pole_squared_length": "-rho*a^2, rho>0",
        "gram_determinant": str(regular_determinant),
        "signature": "(-,+,+,+)",
        "regular_rho_control": rho_numeric,
        "regular_pole_deficit": regular_deficit,
        "regge_stationary": False,
    },
    "external_analytic_control": {
        "statement": "an independent set is a 4D kissing configuration; k(4)=24",
        "reference": "O. R. Musin, The kissing number in four dimensions, Annals of Mathematics 168 (2008), 1-32",
        "url": "https://annals.math.princeton.edu/2008/168-1/p01",
    },
    "labels": {
        "minimal_five_phase_schedule": "DERIVED COMPUTATIONAL",
        "one_unordered_schedule_up_to_H4": "DERIVED COMPUTATIONAL",
        "two_ordered_phase_parity_orbits": "DERIVED COMPUTATIONAL",
        "phased_foliation_as_physics": "STRUCTURAL",
        "relation_to_a1_5": "NOT TESTED IN THIS PROTOCOL",
        "global_Regge_stationarity": "OPEN",
        "phase_order_action_independence": "OPEN",
        "physical_tick_scale_c_G_Planck_units": "OPEN",
    },
    "verdict": verdict,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(verdict)
print("OPEN: full Regge stationarity, order independence and physical units")
raise SystemExit(0 if passed == tests else 1)

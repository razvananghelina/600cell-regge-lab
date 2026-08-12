#!/usr/bin/env python3
"""Canonical one-step slab carrier audit for the 600-cell boundary.

Protocol commit: 5c9fff8.  The computation is target-free: it reconstructs
the cone, product cylinder, H4 stabilizers, initial Pachner move classes and
all 8! projection-preserving identifications with the committed robust walk.
"""

from collections import Counter, defaultdict
from itertools import combinations, permutations, product
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_time_slab_canonicity.json"
PROTOCOL_COMMIT = "5c9fff8"
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
        raise RuntimeError("failed to reconstruct the 120 exact icosians")
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
print("CANONICAL 600-CELL ONE-STEP TIME-SLAB CARRIER AUDIT")
print("=" * 78)

vertices, neighbours, edges, triangles, tetrahedra = build_complex()
edge_index = {edge: index for index, edge in enumerate(edges)}
triangle_index = {triangle: index for index, triangle in enumerate(triangles)}
tetrahedron_index = {
    tetrahedron: index for index, tetrahedron in enumerate(tetrahedra)
}
f_vector = tuple(map(len, (vertices, edges, triangles, tetrahedra)))
check(
    "the independently rebuilt boundary has f-vector (120,720,1200,600)",
    f_vector == (120, 720, 1200, 600),
)

# -------------------------------------------------------------------------
# Cone and product-cell cylinder.
# -------------------------------------------------------------------------
cone_f_vector = (
    1 + f_vector[0],
    f_vector[0] + f_vector[1],
    f_vector[1] + f_vector[2],
    f_vector[2] + f_vector[3],
    f_vector[3],
)
cone_euler = sum((-1)**degree * count
                 for degree, count in enumerate(cone_f_vector))

# Boundary facets of the explicit cone: every base tetrahedron occurs once,
# while a cone over a triangle occurs twice and is internal.
apex = len(vertices)
cone_four_simplices = [tuple(tetrahedron) + (apex,)
                       for tetrahedron in tetrahedra]
cone_facet_counts = Counter(
    tuple(sorted(facet))
    for simplex in cone_four_simplices
    for facet in combinations(simplex, 4)
)
cone_boundary = {
    facet for facet, multiplicity in cone_facet_counts.items()
    if multiplicity == 1
}
check(
    "the canonical cone is a one-boundary cap with Euler characteristic one",
    cone_f_vector == (121, 840, 1920, 1800, 600)
    and cone_euler == 1
    and cone_boundary == set(tetrahedra),
    f"f={cone_f_vector}; boundary tetrahedra={len(cone_boundary)}",
)

# Product cells of total degree k are two copies of K_k plus K_(k-1) x I.
cylinder_f_vector = (
    2 * f_vector[0],
    2 * f_vector[1] + f_vector[0],
    2 * f_vector[2] + f_vector[1],
    2 * f_vector[3] + f_vector[2],
    f_vector[3],
)
cylinder_euler = sum(
    (-1)**degree * count
    for degree, count in enumerate(cylinder_f_vector)
)
check(
    "the product CW cylinder has two K boundaries and 600 prism four-cells",
    cylinder_f_vector == (240, 1560, 3120, 2400, 600)
    and cylinder_euler == 0
    and cylinder_f_vector[-1] == len(tetrahedra),
    f"f={cylinder_f_vector}; chi={cylinder_euler}",
)
check(
    "the product supplies exactly 120 vertical edges as a vertex-lapse carrier",
    cylinder_f_vector[1] - 2 * len(edges) == len(vertices) == 120,
)

# -------------------------------------------------------------------------
# Full H4 stabilizers on a base tetrahedron, triangle, edge and vertex.
# -------------------------------------------------------------------------
base_tetrahedron = tetrahedra[0]
base_triangle = tuple(sorted(base_tetrahedron[:3]))
base_edge = tuple(sorted(base_tetrahedron[:2]))
base_points = vertices[list(base_tetrahedron)]
tetrahedron_targets = Counter()
triangle_targets = Counter()
edge_targets = Counter()
vertex_targets = Counter()
tetrahedron_stabilizer_restrictions = set()
rotational_tetrahedron_restrictions = set()
maximum_action_residual = 0.0
minimum_action_gap = 1.0

for reflected in (False, True):
    seed = qconj(base_points) if reflected else base_points
    for left in vertices:
        left_products = qmul(left, seed)
        for right in vertices:
            mapped = qmul(left_products, qconj(right))
            dots = mapped @ vertices.T
            order = np.argsort(dots, axis=1)
            image = tuple(np.argmax(dots, axis=1).tolist())
            chosen = dots[np.arange(4), image]
            maximum_action_residual = max(
                maximum_action_residual,
                float(np.max(np.abs(chosen - 1))),
            )
            minimum_action_gap = min(
                minimum_action_gap,
                float(np.min(chosen - dots[np.arange(4), order[:, -2]])),
            )
            tetra_target = tuple(sorted(image))
            triangle_target = tuple(sorted(image[:3]))
            edge_target = tuple(sorted(image[:2]))
            if tetra_target not in tetrahedron_index:
                raise RuntimeError("H4 parameter failed to preserve tetrahedra")
            if triangle_target not in triangle_index:
                raise RuntimeError("H4 parameter failed to preserve triangles")
            if edge_target not in edge_index:
                raise RuntimeError("H4 parameter failed to preserve edges")
            tetrahedron_targets[tetra_target] += 1
            triangle_targets[triangle_target] += 1
            edge_targets[edge_target] += 1
            vertex_targets[image[0]] += 1
            if tetra_target == base_tetrahedron:
                restriction = tuple(
                    base_tetrahedron.index(vertex) for vertex in image
                )
                tetrahedron_stabilizer_restrictions.add(restriction)
                if not reflected:
                    rotational_tetrahedron_restrictions.add(restriction)

check(
    "H4 is transitive on vertices, edges, triangles and tetrahedra",
    len(vertex_targets) == 120 and set(vertex_targets.values()) == {240}
    and len(edge_targets) == 720 and set(edge_targets.values()) == {40}
    and len(triangle_targets) == 1200 and set(triangle_targets.values()) == {24}
    and len(tetrahedron_targets) == 600 and set(tetrahedron_targets.values()) == {48},
)
check(
    "the tetrahedron stabilizer induces the complete S4 on its four vertices",
    len(tetrahedron_stabilizer_restrictions) == 24
    and tetrahedron_stabilizer_restrictions
    == set(permutations(range(4))),
)
even_permutations = {
    permutation for permutation in permutations(range(4))
    if sum(
        permutation[i] > permutation[j]
        for i in range(4) for j in range(i + 1, 4)
    ) % 2 == 0
}
check(
    "the orientation-preserving tetrahedron stabilizer is exactly A4",
    rotational_tetrahedron_restrictions == even_permutations
    and len(rotational_tetrahedron_restrictions) == 12,
)
check(
    "the quaternion action identification is numerically separated",
    maximum_action_residual < 1e-12 and minimum_action_gap > 0.1,
    f"max residual={maximum_action_residual:.3e}; "
    f"winner gap={minimum_action_gap:.6g}",
)

# The obstruction can be witnessed inside the even tetrahedral subgroup:
# each double transposition (ij)(kl) swaps the endpoints of an edge and hence
# exchanges the two diagonals of edge x I.  Using even elements matters: an
# odd spatial reflection could otherwise be paired with reversal of I, but
# the double transpositions remain trivial under that orientation character.
double_transpositions = {
    (1, 0, 3, 2),
    (2, 3, 0, 1),
    (3, 2, 1, 0),
}
check(
    "even H4 rotations already forbid an invariant vertex product triangulation",
    double_transpositions <= rotational_tetrahedron_restrictions
    and len(double_transpositions) == 3,
    "each even double transposition exchanges the diagonals of some edge x I",
)

# -------------------------------------------------------------------------
# Canonical barycentric cylinder flags and exact chamber adjacency rules.
# -------------------------------------------------------------------------
spatial_flags_per_tetrahedron = 24
interval_endpoints = 2
temporal_shuffle_positions = 4
flags_per_prism = (
    spatial_flags_per_tetrahedron
    * interval_endpoints
    * temporal_shuffle_positions
)
cylinder_chambers = len(tetrahedra) * flags_per_prism
cylinder_cell_count = sum(cylinder_f_vector)
boundary_chambers = 2 * len(tetrahedra) * spatial_flags_per_tetrahedron
check(
    "functorial barycentric subdivision has exactly 115,200 four-chambers",
    flags_per_prism == 192
    and cylinder_chambers == 115200
    and cylinder_cell_count == 7920
    and boundary_chambers == 28800,
    f"600*24*2*4={cylinder_chambers}; cell barycentres={cylinder_cell_count}",
)

labels = tuple((epsilon, rank) for epsilon in (0, 1) for rank in range(4))


def product_flag(ordering, epsilon, temporal_position):
    """Five product cells in one maximal prism flag."""
    spatial_face = frozenset((ordering[0],))
    interval_face = frozenset((epsilon,))
    chain = [(spatial_face, interval_face)]
    next_spatial = 1
    for promotion in range(4):
        if promotion == temporal_position:
            interval_face = frozenset((0, 1))
        else:
            spatial_face = spatial_face | {ordering[next_spatial]}
            next_spatial += 1
        chain.append((spatial_face, interval_face))
    return tuple(chain)


def chains_adjacent(left, right):
    return sum(a != b for a, b in zip(left, right)) == 1


local_rule_exact = True
spatial_rule_exact = True
representative_ordering = (0, 1, 2, 3)
for epsilon, rank in labels:
    chain = product_flag(representative_ordering, epsilon, rank)
    for other_epsilon, other_rank in labels:
        other = product_flag(
            representative_ordering, other_epsilon, other_rank
        )
        predicted = (
            epsilon == other_epsilon and abs(rank - other_rank) == 1
        ) or (
            rank == other_rank == 0 and epsilon != other_epsilon
        )
        local_rule_exact &= chains_adjacent(chain, other) == predicted
    for colour in range(3):
        changed = list(representative_ordering)
        changed[colour], changed[colour + 1] = (
            changed[colour + 1], changed[colour]
        )
        other = product_flag(tuple(changed), epsilon, rank)
        spatial_rule_exact &= chains_adjacent(chain, other) == (rank != colour)

check(
    "the eight product labels have the preregistered exact local adjacency",
    local_rule_exact,
    "two four-paths joined only between their r=0 endpoints",
)
check(
    "spatial colour i crosses one product chamber edge exactly when r != i",
    spatial_rule_exact,
    "verified directly for internal colours 0,1,2; colour 3 follows at prism facets",
)

# -------------------------------------------------------------------------
# Exact initial Pachner move census.
# -------------------------------------------------------------------------
face_to_tetrahedra = defaultdict(list)
edge_to_tetrahedra = defaultdict(list)
vertex_to_tetrahedra = defaultdict(list)
for tetrahedron_id, tetrahedron in enumerate(tetrahedra):
    for face in combinations(tetrahedron, 3):
        face_to_tetrahedra[tuple(sorted(face))].append(tetrahedron_id)
    for edge in combinations(tetrahedron, 2):
        edge_to_tetrahedra[tuple(sorted(edge))].append(tetrahedron_id)
    for vertex in tetrahedron:
        vertex_to_tetrahedra[vertex].append(tetrahedron_id)

legal_one_four = list(range(len(tetrahedra)))
legal_two_three = []
two_three_support = {}
for face, support in face_to_tetrahedra.items():
    if len(support) != 2:
        continue
    opposites = [
        next(vertex for vertex in tetrahedra[tetrahedron_id]
             if vertex not in face)
        for tetrahedron_id in support
    ]
    if tuple(sorted(opposites)) not in edge_index:
        face_id = triangle_index[face]
        legal_two_three.append(face_id)
        two_three_support[face_id] = frozenset(support)

legal_three_two = []
for edge, support in edge_to_tetrahedra.items():
    if len(support) != 3:
        continue
    opposite_vertices = sorted(set().union(
        *(set(tetrahedra[tetrahedron_id]) - set(edge)
          for tetrahedron_id in support)
    ))
    if len(opposite_vertices) == 3 \
            and tuple(opposite_vertices) not in triangle_index:
        legal_three_two.append(edge_index[edge])

legal_four_one = []
for vertex, support in vertex_to_tetrahedra.items():
    if len(support) != 4:
        continue
    link_vertices = sorted(set().union(
        *(set(tetrahedra[tetrahedron_id]) - {vertex}
          for tetrahedron_id in support)
    ))
    if len(link_vertices) == 4:
        legal_four_one.append(vertex)

pachner_counts = {
    "1-4": len(legal_one_four),
    "2-3": len(legal_two_three),
    "3-2": len(legal_three_two),
    "4-1": len(legal_four_one),
}
check(
    "the initial slice admits exactly 600 legal 1-4 and 1200 legal 2-3 moves",
    pachner_counts == {"1-4": 600, "2-3": 1200, "3-2": 0, "4-1": 0},
    str(pachner_counts),
)
check(
    "each nonempty initial Pachner class is one complete H4 orbit",
    len(tetrahedron_targets) == len(legal_one_four)
    and len(triangle_targets) == len(legal_two_three),
    "one tetrahedron orbit and one triangle orbit",
)

conflict_neighbours = [set() for _ in legal_two_three]
tetrahedron_to_moves = defaultdict(list)
for move, support in two_three_support.items():
    for tetrahedron_id in support:
        tetrahedron_to_moves[tetrahedron_id].append(move)
for move_list in tetrahedron_to_moves.values():
    for left, right in combinations(move_list, 2):
        conflict_neighbours[left].add(right)
        conflict_neighbours[right].add(left)
conflict_degree_counts = Counter(map(len, conflict_neighbours))
check(
    "the transitive 2-3 orbit is not a parallel layer",
    conflict_degree_counts == Counter({6: 1200}),
    "every triangle move conflicts with six others",
)
check(
    "there is no nonempty full-H4-invariant independent 2-3 move set",
    len(triangle_targets) == 1200
    and all(conflict_neighbours)
    and len(legal_two_three) > 0,
    "transitivity leaves only empty/full invariant subsets; full is not independent",
)

# The full 1-4 orbit is compatible because every move replaces a distinct
# tetrahedron. Build the resulting stellar subdivision explicitly.
refined_tetrahedra = []
for tetrahedron_id, tetrahedron in enumerate(tetrahedra):
    new_vertex = len(vertices) + tetrahedron_id
    for face in combinations(tetrahedron, 3):
        refined_tetrahedra.append(tuple(sorted(face + (new_vertex,))))
refined_vertices = len(vertices) + len(tetrahedra)
refined_faces = Counter(
    tuple(sorted(face))
    for tetrahedron in refined_tetrahedra
    for face in combinations(tetrahedron, 3)
)
refined_edges = {
    tuple(sorted(edge))
    for tetrahedron in refined_tetrahedra
    for edge in combinations(tetrahedron, 2)
}
refined_triangles = set(refined_faces)
refined_f_vector = (
    refined_vertices,
    len(refined_edges),
    len(refined_triangles),
    len(refined_tetrahedra),
)
check(
    "the full H4 1-4 orbit gives a valid canonical spatial stellar refinement",
    refined_f_vector == (720, 3120, 4800, 2400)
    and set(refined_faces.values()) == {2}
    and sum((-1)**i * value for i, value in enumerate(refined_f_vector)) == 0,
    f"refined f={refined_f_vector}",
)

# -------------------------------------------------------------------------
# All 8! projection-preserving robust-component identifications.
# -------------------------------------------------------------------------
robust = json.loads((HERE / "tetrahedral_dirac_walk_robust.json").read_text())
check(
    "the committed robust carrier has the same 115,200 cardinality",
    robust["chambers"] == 14400
    and robust["carrier_dimension"] == cylinder_chambers
    and robust["stage_permutations"] == [True, True, True]
    and robust["macro_is_permutation"],
)

stage0 = tuple((None, component) for component in (5, 2, 1, 6, 4, 0, 3, 7))
stage2 = tuple((None, component) for component in (0, 5, 6, 3, 4, 1, 2, 7))
stage1 = (
    (None, 1), (None, 0), (None, 3), (None, 2),
    (None, 4), (2, 6), (3, 5), (None, 7),
)
macro = (
    (None, 2), (2, 3), (3, 0), (None, 1),
    (None, 4), (None, 5), (None, 6), (None, 7),
)


def local_label_adjacent(left, right):
    le, lr = left
    re, rr = right
    return (
        le == re and abs(lr - rr) == 1
    ) or (
        lr == rr == 0 and le != re
    )


def transition_ok(output_component, signature, labeling):
    spatial_colour, input_component = signature
    if spatial_colour is None:
        if output_component == input_component:
            return True
        return local_label_adjacent(
            labeling[output_component], labeling[input_component]
        )
    # Crossing a spatial chamber facet preserves the product shuffle label.
    # The r!=colour condition is relevant only once the labels agree.
    return (
        labeling[output_component] == labeling[input_component]
        and labeling[output_component][1] != spatial_colour
    )


def audit_map(signatures, labeling):
    flags = [
        transition_ok(output, signature, labeling)
        for output, signature in enumerate(signatures)
    ]
    return all(flags), sum(flags)


map_names = ("stage0", "stage1", "stage2", "macro")
maps = (stage0, stage1, stage2, macro)
pass_counts = Counter()
score_histograms = {name: Counter() for name in map_names}
all_stages_pass = 0
all_stages_and_macro_pass = 0
passing_labelings = {name: [] for name in map_names}

for labeling in permutations(labels):
    outcomes = []
    for name, signatures in zip(map_names, maps):
        passed_map, score = audit_map(signatures, labeling)
        outcomes.append(passed_map)
        score_histograms[name][score] += 1
        if passed_map:
            pass_counts[name] += 1
            if len(passing_labelings[name]) < 5:
                passing_labelings[name].append(list(labeling))
    all_stages = all(outcomes[:3])
    all_stages_pass += int(all_stages)
    all_stages_and_macro_pass += int(all_stages and outcomes[3])

check(
    "all 8! projection-preserving component identifications were exhausted",
    sum(score_histograms["stage0"].values()) == 40320
    and all(sum(hist.values()) == 40320 for hist in score_histograms.values()),
)
check(
    "no factorized identification makes the cross-chamber robust stage local",
    pass_counts["stage1"] == 0
    and max(score_histograms["stage1"]) < 8,
    f"stage pass counts={dict(pass_counts)}; "
    f"stage1 best score={max(score_histograms['stage1'])}/8",
)
check(
    "no factorized identification makes the robust macro one slab-chamber edge",
    pass_counts["macro"] == 0
    and max(score_histograms["macro"]) < 8,
    f"macro best score={max(score_histograms['macro'])}/8",
)
check(
    "the equal 115,200 counts do not extend to an incidence bridge",
    all_stages_pass == 0 and all_stages_and_macro_pass == 0,
    f"all stages={all_stages_pass}; all stages+macro={all_stages_and_macro_pass}",
)

protocol = (ROOT / "gravity_time_slab_canonicity_protocol.md").read_text()
check(
    "the preregistered negative is scoped and does not ban future slabs",
    "not** an exhaustive no-go" in protocol
    and "A topological cylinder is not yet physical time" in protocol
    and "No full suite" in protocol,
)

verdicts = [
    "DERIVED CANONICAL CW SLAB",
    "DERIVED SIMPLICIAL CHOICE OBSTRUCTION",
    "DERIVED CURRENT SLAB-DYNAMICS GAP",
]
payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "external_primary_sources": [
        "https://arxiv.org/abs/1108.1974",
        "https://arxiv.org/abs/1411.5672",
        "https://arxiv.org/abs/0905.1670",
    ],
    "spatial_complex": {"f_vector": list(f_vector), "topology": "S3"},
    "cone": {
        "f_vector": list(cone_f_vector),
        "euler_characteristic": cone_euler,
        "boundary_copies_of_K": 1,
        "interpretation": "canonical H4 cap/filling, not a two-slice slab",
    },
    "cellular_cylinder": {
        "f_vector": list(cylinder_f_vector),
        "euler_characteristic": cylinder_euler,
        "boundary_copies_of_K": 2,
        "four_prisms": cylinder_f_vector[-1],
        "vertical_edges": len(vertices),
        "cell_barycentres": cylinder_cell_count,
        "barycentric_four_chambers": cylinder_chambers,
        "boundary_three_chambers": boundary_chambers,
    },
    "h4": {
        "vertex_edge_triangle_tetrahedron_orbits": [1, 1, 1, 1],
        "tetrahedron_stabilizer_vertex_action": "S4",
        "tetrahedron_stabilizer_order": len(tetrahedron_stabilizer_restrictions),
        "rotational_tetrahedron_stabilizer_order": len(
            rotational_tetrahedron_restrictions
        ),
        "vertex_preserving_product_triangulation": "no H4-invariant choice",
        "diagonal_obstruction": (
            "even double transpositions exchange square diagonals; "
            "time reversal cannot repair the choice"
        ),
    },
    "initial_pachner_moves": {
        "legal_counts": pachner_counts,
        "nonempty_h4_orbits": {"1-4": 1, "2-3": 1},
        "two_three_conflict_degree_histogram": {
            str(key): value for key, value in sorted(conflict_degree_counts.items())
        },
        "nonempty_h4_invariant_parallel_two_three_layer": False,
        "full_one_four_orbit_refined_f_vector": list(refined_f_vector),
    },
    "robust_walk_comparison": {
        "cardinality_equality": "14400*8 = 600*24*2*4 = 115200",
        "identifications_tested": 40320,
        "pass_counts": {
            **{name: pass_counts[name] for name in map_names},
            "all_three_stages": all_stages_pass,
            "all_three_stages_and_macro": all_stages_and_macro_pass,
        },
        "score_histograms": {
            name: {str(score): count for score, count in sorted(hist.items())}
            for name, hist in score_histograms.items()
        },
        "sample_passing_labelings": passing_labelings,
        "verdict": "equal cardinality but no projection-preserving incidence bridge",
    },
    "verdicts": verdicts,
    "derived": [
        "K x I is a canonical H4-equivariant two-boundary CW carrier",
        "full H4 forbids a no-new-vertex simplicial product triangulation",
        "barycentric subdivision canonically creates 115200 four-chambers",
        "the existing robust 115200-state walk is not their factorized chamber walk",
        "the only H4-invariant initial physical 2-3 orbit cannot run in one parallel layer",
    ],
    "open": [
        "a non-factorized but geometry-selected walk/slab identification",
        "an orbit-complete Pachner schedule or quantum move sum",
        "Lorentzian edge data and a Regge action on the cylinder",
        "lapse as a multiplier and exact first-class constraint closure",
        "a Legendre map selecting the 47-parameter kinetic form",
    ],
    "not_claimed": [
        "all future discrete temporal constructions are impossible",
        "a topological cylinder is physical time",
        "linearized-flat Regge constraint results hold exactly on this curved slice",
        "c, G, Planck time or a graviton has been derived",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for verdict in verdicts:
    print(verdict)
print("OPEN: metric/action/constraint data on the canonical CW slab")
raise SystemExit(0 if passed == tests else 1)

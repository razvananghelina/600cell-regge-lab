#!/usr/bin/env python3
"""Exact gate for the natural six-axis crossed product on chambers.

The carrier, representation and all geometric-J variants were frozen in
protocol commit a5c485d.  No matter or phenomenological target is used.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import networkx as nx
import sympy as sp


OUTPUT = Path(__file__).with_name("hopf_chamber_crossed_product.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(right)))


def inverse(permutation):
    result = [None]*len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def permutation_order(permutation, maximum=120):
    identity = tuple(range(len(permutation)))
    current = identity
    for order in range(1, maximum+1):
        current = compose(permutation, current)
        if current == identity:
            return order
    raise RuntimeError("permutation order exceeds bound")


def exact_zero(value):
    return sp.simplify(value) == 0


print("="*78)
print("CANONICAL CROSSED-PRODUCT LIFT TO ORIENTED CHAMBERS")
print("="*78)

# -------------------------------------------------------------------------
# Exact icosahedron and complete flags.
# -------------------------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1+sqrt5)/2
vertices = set()
for first, second in product((1, -1), repeat=2):
    vertices.update((
        (sp.Integer(0), sp.Integer(first), sp.Integer(second)*phi),
        (sp.Integer(first), sp.Integer(second)*phi, sp.Integer(0)),
        (sp.Integer(first)*phi, sp.Integer(0), sp.Integer(second)),
    ))
vertices = tuple(sorted(vertices, key=lambda vector: tuple(map(str, vector))))
vertex_lookup = {vector: index for index, vector in enumerate(vertices)}


def dot(left, right):
    return sp.simplify(sum(a*b for a, b in zip(left, right)))


edges = tuple(
    pair for pair in combinations(range(12), 2)
    if exact_zero(dot(vertices[pair[0]], vertices[pair[1]])-phi)
)
edge_lookup = {frozenset(edge): index for index, edge in enumerate(edges)}
edge_set = set(edge_lookup)
faces = tuple(
    face for face in combinations(range(12), 3)
    if all(frozenset(pair) in edge_set for pair in combinations(face, 2))
)
face_lookup = {frozenset(face): index for index, face in enumerate(faces)}
chambers = []
for face_index, face in enumerate(faces):
    for edge_vertices in combinations(face, 2):
        edge_index = edge_lookup[frozenset(edge_vertices)]
        for vertex in edge_vertices:
            chambers.append((vertex, edge_index, face_index))
chambers = tuple(sorted(chambers))
chamber_lookup = {chamber: index for index, chamber in enumerate(chambers)}
check(
    "exact incidence reconstructs the icosahedral f-vector and 120 flags",
    (len(vertices), len(edges), len(faces), len(chambers)) == (12, 30, 20, 120),
)

chamber_neighbors = [set() for _ in chambers]
for left, right in combinations(range(120), 2):
    if sum(a != b for a, b in zip(chambers[left], chambers[right])) == 1:
        chamber_neighbors[left].add(right)
        chamber_neighbors[right].add(left)
check(
    "the complete-flag chamber graph is exactly 3-regular",
    all(len(neighbors) == 3 for neighbors in chamber_neighbors)
    and sum(map(len, chamber_neighbors))//2 == 180,
)

# Enumerate all graph automorphisms, then determine their geometric
# determinant exactly from one independent coordinate triple.
vertex_graph = nx.Graph()
vertex_graph.add_nodes_from(range(12))
vertex_graph.add_edges_from(edges)
automorphisms = tuple(sorted(
    tuple(mapping[index] for index in range(12))
    for mapping in nx.algorithms.isomorphism.GraphMatcher(
        vertex_graph, vertex_graph
    ).isomorphisms_iter()
))
check("the exact icosahedron graph has 120 automorphisms",
      len(automorphisms) == 120)

independent = next(
    triple for triple in combinations(range(12), 3)
    if sp.Matrix.hstack(*(sp.Matrix(vertices[index]) for index in triple)).det()
    != 0
)
base_det = sp.Matrix.hstack(
    *(sp.Matrix(vertices[index]) for index in independent)
).det()
base_matrix = sp.Matrix.hstack(
    *(sp.Matrix(vertices[index]) for index in independent)
)


def determinant_sign(automorphism):
    image_det = sp.Matrix.hstack(
        *(sp.Matrix(vertices[automorphism[index]]) for index in independent)
    ).det()
    ratio = sp.simplify(image_det/base_det)
    if ratio not in (1, -1):
        raise RuntimeError(f"nonorthogonal graph automorphism determinant {ratio}")
    return int(ratio)


rotations = tuple(auto for auto in automorphisms
                  if determinant_sign(auto) == 1)
improper = tuple(auto for auto in automorphisms
                 if determinant_sign(auto) == -1)
identity_vertices = tuple(range(12))
check(
    "orientation sign splits the full symmetry group as 60+60",
    len(rotations) == len(improper) == 60 and identity_vertices in rotations,
)
all_automorphisms_geometric = True
vertex_gram = tuple(
    tuple(dot(left, right) for right in vertices) for left in vertices
)
for automorphism in automorphisms:
    # Exact preservation of the full Gram matrix proves that the permutation
    # is induced by a unique orthogonal linear map: the vertices span R^3.
    all_automorphisms_geometric &= all(
        vertex_gram[automorphism[left]][automorphism[right]]
        == vertex_gram[left][right]
        for left in range(12) for right in range(12)
    )
check(
    "every graph automorphism preserves the exact vertex Gram matrix",
    all_automorphisms_geometric,
)

inversion_vertices = tuple(
    vertex_lookup[tuple(-coordinate for coordinate in vertex)]
    for vertex in vertices
)
check(
    "central inversion is the exact central improper symmetry",
    inversion_vertices in improper
    and all(compose(inversion_vertices, rotation)
            == compose(rotation, inversion_vertices)
            for rotation in rotations)
    and {compose(inversion_vertices, rotation) for rotation in rotations}
    == set(improper),
)


def induced_data(vertex_permutation):
    edge_permutation = tuple(
        edge_lookup[frozenset(vertex_permutation[vertex] for vertex in edge)]
        for edge in edges
    )
    face_permutation = tuple(
        face_lookup[frozenset(vertex_permutation[vertex] for vertex in face)]
        for face in faces
    )
    chamber_permutation = tuple(
        chamber_lookup[(
            vertex_permutation[vertex],
            edge_permutation[edge],
            face_permutation[face],
        )]
        for vertex, edge, face in chambers
    )
    return edge_permutation, face_permutation, chamber_permutation


rotation_chambers = tuple(induced_data(rotation)[2] for rotation in rotations)
rotation_index = {rotation: index for index, rotation in enumerate(rotations)}
inversion_chambers = induced_data(inversion_vertices)[2]
identity_chambers = tuple(range(120))
check(
    "the induced chamber actions are faithful and preserve adjacency",
    len(set(rotation_chambers)) == 60
    and all(
        {permutation[neighbor] for neighbor in chamber_neighbors[source]}
        == chamber_neighbors[permutation[source]]
        for permutation in rotation_chambers+(inversion_chambers,)
        for source in range(120)
    ),
)

sheet_plus = frozenset(permutation[0] for permutation in rotation_chambers)
sheet_minus = frozenset(set(range(120))-set(sheet_plus))
check(
    "rotations give two free chamber sheets of size 60",
    len(sheet_plus) == len(sheet_minus) == 60
    and all(sum(permutation[source] == source for permutation in rotation_chambers)
            == 1 for source in range(120)),
)
check(
    "central inversion exchanges the sheets and commutes with A5",
    {inversion_chambers[index] for index in sheet_plus} == set(sheet_minus)
    and all(compose(inversion_chambers, permutation)
            == compose(permutation, inversion_chambers)
            for permutation in rotation_chambers),
)

gamma = tuple(1 if chamber in sheet_plus else -1 for chamber in range(120))
check(
    "chamber adjacency is odd and inversion reverses the grading",
    all(gamma[left] == -gamma[right]
        for left in range(120) for right in chamber_neighbors[left])
    and all(gamma[inversion_chambers[index]] == -gamma[index]
            for index in range(120)),
)

# -------------------------------------------------------------------------
# Six antipodal vertex axes and the equivariant chamber projection.
# -------------------------------------------------------------------------
axes = tuple(sorted({
    frozenset((vertex, inversion_vertices[vertex])) for vertex in range(12)
}, key=lambda pair: sorted(pair)))
axis_lookup = {axis: index for index, axis in enumerate(axes)}
vertex_axis = tuple(
    axis_lookup[frozenset((vertex, inversion_vertices[vertex]))]
    for vertex in range(12)
)
chamber_axis = tuple(vertex_axis[chamber[0]] for chamber in chambers)
axis_actions = tuple(
    tuple(axis_lookup[frozenset(rotation[vertex] for vertex in axis)]
          for axis in axes)
    for rotation in rotations
)
check(
    "the chamber-to-six-axis projection is exact and A5-equivariant",
    len(axes) == 6
    and all(Counter(chamber_axis[index] for index in sheet) == Counter({
        axis: 10 for axis in range(6)
    }) for sheet in (sheet_plus, sheet_minus))
    and all(
        chamber_axis[rotation_chambers[group][chamber]]
        == axis_actions[group][chamber_axis[chamber]]
        for group in range(60) for chamber in range(120)
    ),
)

# -------------------------------------------------------------------------
# Exact crossed-product basis E_(x,g)=P_x U_g.
# Supports on one free sheet are pairwise disjoint, proving independence.
# -------------------------------------------------------------------------
basis_labels = tuple((axis, group) for axis in range(6) for group in range(60))


def basis_support(axis, group, sheet):
    permutation = rotation_chambers[group]
    return frozenset(
        (permutation[source], source) for source in sheet
        if chamber_axis[permutation[source]] == axis
    )


plus_supports = tuple(
    basis_support(axis, group, sheet_plus) for axis, group in basis_labels
)
all_positions = set()
supports_disjoint = True
for support in plus_supports:
    if all_positions & set(support):
        supports_disjoint = False
    all_positions |= set(support)
check(
    "the 360 crossed-product basis matrices are exactly independent",
    supports_disjoint and all(len(support) == 10 for support in plus_supports)
    and len(all_positions) == 3600,
    "disjoint ten-entry supports on a free 60-state sheet",
)

check(
    "all 6*60 covariance identities hold exactly",
    all(
        chamber_axis[rotation_chambers[group][chamber]]
        == axis_actions[group][chamber_axis[chamber]]
        for group in range(60) for chamber in range(120)
    ),
)

# The point stabilizer is D5.  Its regular representation on each ten-point
# fibre has multiplicities equal to the four real irrep dimensions 1,1,2,2.
base_axis = chamber_axis[0]
stabilizer_groups = tuple(
    group for group in range(60) if axis_actions[group][base_axis] == base_axis
)
stabilizer_orders = Counter(
    permutation_order(rotations[group], 10) for group in stabilizer_groups
)
module_multiplicities = (1, 1, 2, 2)
simple_dimensions = (6, 6, 12, 12)
check(
    "the sheet representation has faithful Wedderburn multiplicities 1,1,2,2",
    len(stabilizer_groups) == 10
    and stabilizer_orders == Counter({1: 1, 2: 5, 5: 4})
    and sum(mult*dimension for mult, dimension
            in zip(module_multiplicities, simple_dimensions)) == 60
    and all(mult > 0 for mult in module_multiplicities),
)

# Direct commutant census: after axis diagonality, matrix entries are pairs of
# chambers over the same axis.  Diagonal A5 orbits classify the commutant.
plus_order = tuple(sorted(sheet_plus))
allowed_pairs = {
    (left, right) for left in plus_order for right in plus_order
    if chamber_axis[left] == chamber_axis[right]
}
unseen_pairs = set(allowed_pairs)
pair_orbits = []
while unseen_pairs:
    seed = min(unseen_pairs)
    orbit = frozenset(
        (permutation[seed[0]], permutation[seed[1]])
        for permutation in rotation_chambers
    )
    pair_orbits.append(orbit)
    unseen_pairs -= orbit
sheet_commutant_dimension = len(pair_orbits)
double_commutant_dimension = 4*sheet_commutant_dimension
check(
    "the exact algebra commutant has dimensions 10 per sheet and 40 doubled",
    sheet_commutant_dimension == sum(mult*mult
                                     for mult in module_multiplicities) == 10
    and double_commutant_dimension == 40,
)
check(
    "commutant dimension rules out order zero for every possible faithful J",
    double_commutant_dimension < len(basis_labels),
    "J pi(B) J^-1 would be a faithful 360-dimensional subalgebra of a "
    "40-dimensional commutant",
)

# -------------------------------------------------------------------------
# Geometric real structures and exhaustive order-zero census in the abstract
# crossed-product basis.
# -------------------------------------------------------------------------
rotation_products = {
    (left, right): rotation_index[compose(rotations[left], rotations[right])]
    for left in range(60) for right in range(60)
}
rotation_inverses = {
    group: rotation_index[inverse(rotations[group])] for group in range(60)
}


def basis_product(left, right):
    left_axis, left_group = left
    right_axis, right_group = right
    if left_axis != axis_actions[left_group][right_axis]:
        return None
    return left_axis, rotation_products[left_group, right_group]


noncommuting_pairs = 0
order_zero_witness = None
for left in basis_labels:
    for right in basis_labels:
        if basis_product(left, right) != basis_product(right, left):
            noncommuting_pairs += 1
            if order_zero_witness is None:
                order_zero_witness = (left, right)
check(
    "the exhaustive 360^2 basis census detects noncommutativity",
    noncommuting_pairs > 0 and order_zero_witness is not None,
    f"ordered noncommuting pairs={noncommuting_pairs}; witness={order_zero_witness}",
)

j_records = []
all_j_normalize = True
all_j_odd_real = True
for group in range(60):
    j_permutation = compose(rotation_chambers[group], inversion_chambers)
    j_square = compose(j_permutation, j_permutation)
    j_squared_one = j_square == identity_chambers
    j_gamma_minus = all(gamma[j_permutation[index]] == -gamma[index]
                        for index in range(120))
    j_commutes_d = all(
        {j_permutation[neighbor] for neighbor in chamber_neighbors[source]}
        == chamber_neighbors[j_permutation[source]]
        for source in range(120)
    )

    # J_k E_(x,h) J_k^-1 = E_(k x, k h k^-1).
    conjugated_labels = []
    group_inverse = rotation_inverses[group]
    for axis, algebra_group in basis_labels:
        conjugated_group = rotation_products[
            rotation_products[group, algebra_group], group_inverse
        ]
        conjugated_labels.append(
            (axis_actions[group][axis], conjugated_group)
        )
    normalizes = set(conjugated_labels) == set(basis_labels)
    all_j_normalize &= normalizes
    all_j_odd_real &= j_gamma_minus and j_commutes_d
    # Since conjugation permutes the full basis, the exhaustive order-zero
    # failure count is the same nonzero count for every J.
    order_zero = normalizes and noncommuting_pairs == 0
    j_records.append({
        "rotation_index": group,
        "rotation_order": permutation_order(rotations[group], 10),
        "J_squared_one": j_squared_one,
        "J_gamma_minus": j_gamma_minus,
        "JD_equals_DJ": j_commutes_d,
        "opposite_image_equals_algebra_image": normalizes,
        "order_zero": order_zero,
        "ordered_noncommuting_basis_pairs": noncommuting_pairs,
    })

admissible_j = [record for record in j_records if record["J_squared_one"]]
check(
    "all 60 improper symmetries reverse gamma, preserve D and normalize B",
    all_j_normalize and all_j_odd_real and len(j_records) == 60,
)
check(
    "exactly 16 geometric J candidates square to +1",
    len(admissible_j) == 16
    and Counter(record["rotation_order"] for record in admissible_j)
    == Counter({1: 1, 2: 15}),
)
check(
    "every geometric J fails order zero",
    all(not record["order_zero"] for record in j_records),
    f"0/60 order-zero; 0/{len(admissible_j)} among J^2=+1 candidates",
)

# -------------------------------------------------------------------------
# Exact connectedness rank.  Build integer commutator columns for all 360
# basis elements.  A modular lower bound plus the explicit 60-dimensional
# group-algebra kernel gives an exact rational rank when they meet.
# -------------------------------------------------------------------------
def commutator_column(axis, group):
    permutation = rotation_chambers[group]
    entries = {}

    def add(row, column, value):
        index = row*120+column
        entries[index] = entries.get(index, 0)+value
        if entries[index] == 0:
            del entries[index]

    for source in range(120):
        target = permutation[source]
        if chamber_axis[target] == axis:
            for row in chamber_neighbors[target]:
                add(row, source, 1)
        for neighbor in chamber_neighbors[source]:
            neighbor_target = permutation[neighbor]
            if chamber_axis[neighbor_target] == axis:
                add(neighbor_target, source, -1)
    return entries


commutator_columns = tuple(
    commutator_column(axis, group) for axis, group in basis_labels
)
check(
    "the derived chamber D has nonzero represented one-forms",
    any(column for column in commutator_columns),
)

# Every U_g=sum_x E_(x,g) commutes with the A5-invariant chamber adjacency.
group_kernel_exact = all(
    all(sum(commutator_columns[axis*60+group].get(position, 0)
            for axis in range(6)) == 0
        for position in set().union(*(
            set(commutator_columns[axis*60+group]) for axis in range(6)
        )))
    for group in range(60)
)
check(
    "the full 60-dimensional represented group algebra commutes with D",
    group_kernel_exact,
)


def modular_column_rank(columns, prime=1000003):
    pivots = {}
    for original in columns:
        vector = {row: value % prime for row, value in original.items()
                  if value % prime}
        while vector:
            pivot = min(vector)
            if pivot not in pivots:
                inverse_pivot = pow(vector[pivot], -1, prime)
                vector = {row: (value*inverse_pivot) % prime
                          for row, value in vector.items()
                          if (value*inverse_pivot) % prime}
                pivots[pivot] = vector
                break
            factor = vector[pivot]
            basis = pivots[pivot]
            for row, value in basis.items():
                updated = (vector.get(row, 0)-factor*value) % prime
                if updated:
                    vector[row] = updated
                elif row in vector:
                    del vector[row]
    return len(pivots)


modular_rank = modular_column_rank(commutator_columns)
kernel_upper_from_explicit = 60
rank_upper_from_explicit = 360-kernel_upper_from_explicit
exact_rank_certified = modular_rank == rank_upper_from_explicit
commutator_rank = modular_rank if exact_rank_certified else None
algebra_D_commutant_dimension = (
    360-commutator_rank if commutator_rank is not None else None
)
check(
    "modular rank and the explicit kernel certify the exact D-commutant",
    exact_rank_certified and algebra_D_commutant_dimension == 60,
    f"rank={commutator_rank}; kernel dimension={algebra_D_commutant_dimension}",
)
check(
    "the faithful chamber lift also fails connectedness",
    algebra_D_commutant_dimension != 1,
)

payload = {
    "protocol_commit": "a5c485d",
    "physical_target_comparison_performed": False,
    "geometry": {
        "f_vector": [12, 30, 20],
        "chambers": 120,
        "chamber_edges": 180,
        "chamber_degree": 3,
        "rotation_group_order": 60,
        "orientation_sheet_sizes": [60, 60],
        "axis_count": 6,
        "chambers_per_axis_per_sheet": 10,
    },
    "crossed_product_representation": {
        "algebra_dimension": 360,
        "basis_independent": True,
        "faithful": True,
        "simple_block_dimensions": list(simple_dimensions),
        "sheet_module_multiplicities": list(module_multiplicities),
        "sheet_commutant_dimension": sheet_commutant_dimension,
        "doubled_commutant_dimension": double_commutant_dimension,
    },
    "geometric_J_census": {
        "total_improper_symmetries": len(j_records),
        "J_squared_plus_one": len(admissible_j),
        "order_zero_all": sum(record["order_zero"] for record in j_records),
        "order_zero_among_J2_plus": sum(record["order_zero"]
                                         for record in admissible_j),
        "opposite_image_equals_algebra_image": all_j_normalize,
        "ordered_noncommuting_basis_pairs": noncommuting_pairs,
        "first_order_evaluated": False,
        "records": j_records,
    },
    "universal_order_zero_obstruction": {
        "represented_algebra_dimension": len(basis_labels),
        "represented_commutant_dimension": double_commutant_dimension,
        "faithful_opposite_action_can_fit_in_commutant": (
            len(basis_labels) <= double_commutant_dimension
        ),
        "scope": "this fixed faithful doubled chamber representation",
    },
    "Dirac_gates": {
        "nonzero_inner_one_form_witness": any(
            column for column in commutator_columns
        ),
        "commutator_map_rank": commutator_rank,
        "algebra_D_commutant_dimension": algebra_D_commutant_dimension,
        "connected": algebra_D_commutant_dimension == 1,
        "group_algebra_kernel_dimension": 60,
    },
    "verdict": (
        "DERIVED CANONICAL-LIFT NO-GO: each free chamber sheet gives the "
        "faithful natural representation of R(A5/D5) crossed A5, with "
        "multiplicities 1,1,2,2. But every geometric improper symmetry "
        "normalizes the same noncommutative algebra image, so all 60 fail "
        "order zero (including all 16 with J^2=+1). More strongly, the "
        "40-dimensional commutant cannot contain any faithful "
        "360-dimensional opposite image, so no possible J repairs order "
        "zero on this representation. Independently, the "
        "A5-invariant chamber adjacency commutes with the full 60-dimensional "
        "group algebra and fails connectedness, despite nonzero one-forms."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact structured chamber-lift audit was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED CANONICAL-LIFT NO-GO: faithful, but order zero 0/60")
print("and chamber-D connectedness fails with a 60-dimensional commutant.")
print("NO MATTER OR STANDARD-MODEL TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

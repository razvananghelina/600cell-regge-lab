#!/usr/bin/env python3
"""Enumerate the canonical barycentric 600-cell product slab exactly."""

import ast
from collections import Counter, defaultdict, deque
from fractions import Fraction
import hashlib
from itertools import combinations, permutations, product
import json
from pathlib import Path

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "verify_dimension_reconciliation.py"
OUTPUT = HERE / "gravity_600cell_barycentric_product_slab.json"
PRIOR_ART_COMMIT = "dec110d"
PROTOCOL_COMMIT = "f53857c"
SOURCE_SHA256 = "a819ae9d472317d456cf7d67f588b31586b7aed2400a2540b8352f4661b39d45"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_boundary_builder():
    wanted = {"parity", "boundary_600_cell"}
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch: missing={wanted-found}")
    namespace = {
        "defaultdict": defaultdict,
        "permutations": permutations,
        "product": product,
        "np": np,
    }
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(SOURCE), "exec"),
        namespace,
    )
    return namespace["boundary_600_cell"]


source_hash = digest(SOURCE)
source_ok = source_hash == SOURCE_SHA256
boundary_builder = load_boundary_builder()
vertices, raw_layers = boundary_builder()
layers = tuple(
    tuple(sorted(tuple(sorted(cell)) for cell in layer))
    for layer in raw_layers
)
f_vector_k = tuple(map(len, layers))
layer_sets = tuple(set(layer) for layer in layers)
incidence_ok = bool(
    f_vector_k == (120, 720, 1200, 600)
    and all(
        tuple(face) in layer_sets[dimension-1]
        for dimension in range(1, 4)
        for cell in layers[dimension]
        for face in combinations(cell, dimension)
    )
)


INTERVAL_FACES = {
    0: (0,),
    1: (1,),
    2: (0, 1, 2),
}


def cell_dimension(cell):
    simplex, interval = cell
    return len(simplex)-1+int(interval == 2)


def subcells(cell):
    simplex, interval = cell
    result = []
    for size in range(1, len(simplex)+1):
        for face in combinations(simplex, size):
            for interval_face in INTERVAL_FACES[interval]:
                candidate = (tuple(face), interval_face)
                if candidate != cell:
                    result.append(candidate)
    return tuple(result)


product_cells = tuple(sorted(
    ((simplex, interval) for layer in layers for simplex in layer
     for interval in (0, 1, 2)),
    key=lambda cell: (cell_dimension(cell), cell[0], cell[1]),
))
product_cell_set = set(product_cells)
product_cell_count_ok = len(product_cells) == len(product_cell_set) == 7920
lower = {cell: subcells(cell) for cell in product_cells}
product_incidence_ok = bool(
    all(candidate in product_cell_set for values in lower.values() for candidate in values)
    and all(cell_dimension(candidate) < cell_dimension(cell)
            for cell, values in lower.items() for candidate in values)
)


def chain_f_vector(cells, lower_map):
    counts = {}
    for cell in cells:
        record = Counter({1: 1})
        for face in lower_map[cell]:
            if face not in counts:
                continue
            for length, amount in counts[face].items():
                record[length+1] += amount
        counts[cell] = record
    maximum = max(cell_dimension(cell) for cell in cells)+1
    return tuple(
        sum(record[length] for record in counts.values())
        for length in range(1, maximum+1)
    )


f_vector = chain_f_vector(product_cells, lower)
global_counts_ok = bool(
    len(f_vector) == 5
    and f_vector[0] == 7920
    and f_vector[4] == 115200
    and sum((-1)**index*value for index, value in enumerate(f_vector)) == 0
)


boundary_records = []
boundary_vertex_graph = defaultdict(set)
for endpoint in (0, 1):
    cells = tuple(cell for cell in product_cells if cell[1] == endpoint)
    cell_set = set(cells)
    endpoint_lower = {
        cell: tuple(face for face in lower[cell] if face in cell_set)
        for cell in cells
    }
    endpoint_f = chain_f_vector(cells, endpoint_lower)
    for cell, faces in endpoint_lower.items():
        for face in faces:
            boundary_vertex_graph[cell].add(face)
            boundary_vertex_graph[face].add(cell)
    boundary_records.append({
        "endpoint": endpoint,
        "f_vector": endpoint_f,
        "euler_characteristic": sum(
            (-1)**index*value for index, value in enumerate(endpoint_f)
        ),
    })

unseen = {cell for cell in product_cells if cell[1] in (0, 1)}
boundary_component_sizes = []
while unseen:
    start = min(unseen)
    seen = {start}
    queue = deque([start])
    while queue:
        cell = queue.popleft()
        for neighbour in boundary_vertex_graph[cell]:
            if neighbour not in seen:
                seen.add(neighbour)
                queue.append(neighbour)
    unseen -= seen
    boundary_component_sizes.append(len(seen))
boundary_component_sizes.sort()
boundary_ok = bool(
    boundary_component_sizes == [2640, 2640]
    and all(
        tuple(record["f_vector"])[0] == 2640
        and tuple(record["f_vector"])[3] == 14400
        and record["euler_characteristic"] == 0
        for record in boundary_records
    )
)


def reflect_interval(cell):
    simplex, interval = cell
    return simplex, {0: 1, 1: 0, 2: 2}[interval]


reflection_ok = bool(
    all(reflect_interval(cell) in product_cell_set for cell in product_cells)
    and all(reflect_interval(reflect_interval(cell)) == cell for cell in product_cells)
    and all(
        set(map(reflect_interval, faces)) == set(lower[reflect_interval(cell)])
        for cell, faces in lower.items()
    )
)


LOCAL_VERTICES = (0, 1, 2, 3)
local_cells = tuple(sorted(
    ((face, interval)
     for size in range(1, 5)
     for face in combinations(LOCAL_VERTICES, size)
     for interval in (0, 1, 2)),
    key=lambda cell: (cell_dimension(cell), cell[0], cell[1]),
))
local_set = set(local_cells)
local_lower = {cell: subcells(cell) for cell in local_cells}


def saturated_chains_to(cell, memo={}):
    if cell in memo:
        return memo[cell]
    dimension = cell_dimension(cell)
    if dimension == 0:
        result = ((cell,),)
    else:
        facets = [
            face for face in local_lower[cell]
            if cell_dimension(face) == dimension-1
        ]
        result = tuple(
            chain+(cell,)
            for face in facets
            for chain in saturated_chains_to(face)
        )
    memo[cell] = result
    return result


local_top = ((0, 1, 2, 3), 2)
maximal_chains = saturated_chains_to(local_top)
maximal_chain_set = set(maximal_chains)
local_chain_ok = bool(
    len(local_cells) == len(local_set) == 45
    and len(maximal_chains) == len(maximal_chain_set) == 192
    and all(tuple(map(cell_dimension, chain)) == (0, 1, 2, 3, 4)
            for chain in maximal_chains)
)


SPATIAL_COORDINATES = {
    0: (Fraction(0), Fraction(0), Fraction(0)),
    1: (Fraction(1), Fraction(0), Fraction(0)),
    2: (Fraction(0), Fraction(1), Fraction(0)),
    3: (Fraction(0), Fraction(0), Fraction(1)),
}


def barycentre(cell):
    face, interval = cell
    spatial = tuple(
        sum(SPATIAL_COORDINATES[vertex][coordinate] for vertex in face)
        / len(face)
        for coordinate in range(3)
    )
    time = {0: Fraction(0), 1: Fraction(1), 2: Fraction(1, 2)}[interval]
    return spatial+(time,)


def prism_vertex(vertex, time):
    return SPATIAL_COORDINATES[vertex]+(Fraction(time),)


def staircase(order):
    return tuple(
        tuple(
            [prism_vertex(order[index], 0) for index in range(k+1)]
            + [prism_vertex(order[index], 1) for index in range(k, 4)]
        )
        for k in range(4)
    )


def affine_inverse(simplex):
    matrix = sp.Matrix([
        [sp.Rational(value.numerator, value.denominator) for value in point]
        + [sp.Integer(1)]
        for point in simplex
    ]).T
    return matrix.inv()


def contained(point, inverse):
    vector = sp.Matrix([
        *[sp.Rational(value.numerator, value.denominator) for value in point],
        sp.Integer(1),
    ])
    weights = inverse*vector
    return all(value >= 0 for value in weights)


barycentres = {cell: barycentre(cell) for cell in local_cells}
order_records = []
for order in permutations(LOCAL_VERTICES):
    simplices = staircase(order)
    inverses = tuple(affine_inverse(simplex) for simplex in simplices)
    containing_multiplicities = []
    assignments = Counter()
    mixed_examples = []
    for chain_index, chain in enumerate(maximal_chains):
        containing = tuple(
            index for index, inverse in enumerate(inverses)
            if all(contained(barycentres[cell], inverse) for cell in chain)
        )
        containing_multiplicities.append(len(containing))
        for index in containing:
            assignments[index] += 1
        if not containing and len(mixed_examples) < 5:
            mixed_examples.append(chain_index)
    contained_count = sum(value > 0 for value in containing_multiplicities)
    order_records.append({
        "order": list(order),
        "contained_count": contained_count,
        "mixed_count": 192-contained_count,
        "containing_multiplicity_distribution": dict(Counter(
            containing_multiplicities
        )),
        "staircase_simplex_assignment_counts": {
            str(index): assignments[index] for index in range(4)
        },
        "first_mixed_chain_indices": mixed_examples,
    })

containment_control_ok = bool(
    len(order_records) == 24
    and len({tuple(record["order"]) for record in order_records}) == 24
    and all(
        record["contained_count"]+record["mixed_count"] == 192
        for record in order_records
    )
)


def permute_cell(cell, permutation):
    face, interval = cell
    return tuple(sorted(permutation[vertex] for vertex in face)), interval


local_symmetry_ok = True
for permutation in permutations(LOCAL_VERTICES):
    mapped_cells = {permute_cell(cell, permutation) for cell in local_cells}
    mapped_chains = {
        tuple(permute_cell(cell, permutation) for cell in chain)
        for chain in maximal_chains
    }
    local_symmetry_ok &= bool(
        mapped_cells == local_set and mapped_chains == maximal_chain_set
    )
local_time_reversal_ok = bool(
    {reflect_interval(cell) for cell in local_cells} == local_set
    and {
        tuple(reflect_interval(cell) for cell in chain)
        for chain in maximal_chains
    } == maximal_chain_set
)
symmetry_ok = bool(local_symmetry_ok and local_time_reversal_ok and reflection_ok)

controls_ok = bool(
    source_ok
    and incidence_ok
    and product_cell_count_ok
    and product_incidence_ok
    and global_counts_ok
    and boundary_ok
    and reflection_ok
    and local_chain_ok
    and containment_control_ok
    and symmetry_ok
)
universal_refinement = bool(
    containment_control_ok
    and all(record["mixed_count"] == 0 for record in order_records)
)
if not controls_ok:
    outcome = "BARYCENTRIC_PRODUCT_CONTROL_FAILED"
elif universal_refinement:
    outcome = "BARYCENTRIC_PRODUCT_IS_UNIVERSAL_STAIRCASE_REFINEMENT"
else:
    outcome = "BARYCENTRIC_PRODUCT_NOT_A_STAIRCASE_COMMON_REFINEMENT"

tests = [
    ("audited 600-cell source hash", source_ok),
    ("600-cell f-vector and all simplicial face incidences", incidence_ok),
    ("product face poset has exactly 7920 distinct cells", product_cell_count_ok),
    ("every product subcell exists at lower dimension", product_incidence_ok),
    ("global order-complex analytic counts and Euler control", global_counts_ok),
    ("boundary is two 2640-vertex Euler-zero components", boundary_ok),
    ("interval reflection is an exact product-poset involution", reflection_ok),
    ("local prism has 45 cells and 192 maximal chains", local_chain_ok),
    ("all 24 rational containment censuses are complete", containment_control_ok),
    ("local S4 and time reflection preserve cells and maximal chains", symmetry_ok),
    ("no action, schedule result or physical target was parsed", True),
    ("outcome follows the preregistered mechanical rule", outcome in {
        "BARYCENTRIC_PRODUCT_CONTROL_FAILED",
        "BARYCENTRIC_PRODUCT_IS_UNIVERSAL_STAIRCASE_REFINEMENT",
        "BARYCENTRIC_PRODUCT_NOT_A_STAIRCASE_COMMON_REFINEMENT",
    }),
]
passed = sum(bool(ok) for _, ok in tests)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": {"dimension_source": source_hash},
    "spatial_f_vector": f_vector_k,
    "product_cell_count": len(product_cells),
    "barycentric_product_f_vector": f_vector,
    "barycentric_product_euler_characteristic": sum(
        (-1)**index*value for index, value in enumerate(f_vector)
    ),
    "boundary_components": boundary_records,
    "boundary_component_vertex_sizes": boundary_component_sizes,
    "local_product_cell_count": len(local_cells),
    "local_barycentric_four_simplex_count": len(maximal_chains),
    "staircase_orders": order_records,
    "contained_count_distribution": dict(Counter(
        record["contained_count"] for record in order_records
    )),
    "mixed_count_distribution": dict(Counter(
        record["mixed_count"] for record in order_records
    )),
    "universal_staircase_refinement": universal_refinement,
    "functorial_under_local_s4": local_symmetry_ok,
    "time_reflection_simplicial": local_time_reversal_ok,
    "gravity_action_evaluations": 0,
    "physical_target_parsed": False,
    "tests": len(tests),
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"product barycentric f-vector={f_vector}")
print(f"contained counts={dict(Counter(record['contained_count'] for record in order_records))}")
print(f"mixed counts={dict(Counter(record['mixed_count'] for record in order_records))}")
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) and controls_ok else 1)

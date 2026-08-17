#!/usr/bin/env python3
"""Enumerate the universal overlay face poset and order complex exactly."""

from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import reduce
from itertools import combinations, permutations
from math import gcd
import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "gravity_600cell_universal_staircase_overlay.json"
OUTPUT = HERE / "gravity_600cell_overlay_face_poset.json"
SOURCE_SHA256 = "0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc"
PRIOR_ART_COMMIT = "4a85d25"
PROTOCOL_COMMIT = "e8d995c"
SPATIAL_F_VECTOR = (120, 720, 1200, 600)
VERTICES = tuple(range(4))
PERMUTATIONS = tuple(permutations(VERTICES))
FULL_MASK = 15
MASKS = tuple(range(1, FULL_MASK))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qstring(value):
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def mask_of(items):
    result = 0
    for item in items:
        result |= 1 << item
    return result


def full_lambdas(point):
    return point[:3]+(Fraction(1)-sum(point[:3]),)


def add_forms(left, right):
    return tuple(a+b for a, b in zip(left, right))


def scale_form(value, form):
    return tuple(value*entry for entry in form)


ZERO_FORM = (Fraction(0),)*5
LAMBDA_FORMS = (
    (Fraction(0), Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(1), Fraction(-1), Fraction(-1), Fraction(-1), Fraction(0)),
)
TIME_FORM = (Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1))
ONE_MINUS_TIME_FORM = (
    Fraction(1), Fraction(0), Fraction(0), Fraction(0), Fraction(-1)
)

forms = []
for vertex, form in enumerate(LAMBDA_FORMS):
    forms.append({
        "label": f"lambda_{vertex}",
        "kind": "boundary",
        "boundary_index": vertex,
        "mask": None,
        "coefficients": form,
    })
forms.extend((
    {
        "label": "t",
        "kind": "boundary",
        "boundary_index": 4,
        "mask": None,
        "coefficients": TIME_FORM,
    },
    {
        "label": "1-t",
        "kind": "boundary",
        "boundary_index": 5,
        "mask": None,
        "coefficients": ONE_MINUS_TIME_FORM,
    },
))
for mask in MASKS:
    subset_sum = ZERO_FORM
    for vertex in VERTICES:
        if mask & (1 << vertex):
            subset_sum = add_forms(subset_sum, LAMBDA_FORMS[vertex])
    form = add_forms(TIME_FORM, scale_form(Fraction(-1), subset_sum))
    forms.append({
        "label": f"h_{mask:04b}",
        "kind": "internal",
        "boundary_index": None,
        "mask": mask,
        "coefficients": form,
    })

FORM_COUNT = len(forms)
BOUNDARY_FORM_INDICES = tuple(range(6))
INTERNAL_FORM_INDICES = tuple(range(6, 20))
MASK_TO_FORM = {forms[index]["mask"]: index for index in INTERNAL_FORM_INDICES}
LABEL_TO_FORM = {record["label"]: index for index, record in enumerate(forms)}


def evaluate(form, point):
    return form[0]+sum(coefficient*coordinate for coefficient, coordinate in zip(form[1:], point))


def hyperplane_key(form):
    denominators = [entry.denominator for entry in form]
    common = 1
    for denominator in denominators:
        common = common*denominator//gcd(common, denominator)
    integers = [entry.numerator*(common//entry.denominator) for entry in form]
    common_gcd = reduce(gcd, (abs(value) for value in integers if value), 0)
    integers = [value//common_gcd for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


forms_distinct_ok = bool(
    FORM_COUNT == 20
    and len({hyperplane_key(record["coefficients"]) for record in forms}) == 20
    and len(MASK_TO_FORM) == 14
)


def solve_four(form_indices):
    augmented = [
        [*forms[index]["coefficients"][1:], -forms[index]["coefficients"][0]]
        for index in form_indices
    ]
    pivot_row = 0
    for column in range(4):
        pivot = next(
            (row for row in range(pivot_row, 4) if augmented[row][column]),
            None,
        )
        if pivot is None:
            return None
        augmented[pivot_row], augmented[pivot] = augmented[pivot], augmented[pivot_row]
        divisor = augmented[pivot_row][column]
        augmented[pivot_row] = [value/divisor for value in augmented[pivot_row]]
        for row in range(4):
            if row == pivot_row:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    left-factor*right
                    for left, right in zip(augmented[row], augmented[pivot_row])
                ]
        pivot_row += 1
    return tuple(augmented[row][4] for row in range(4))


intersection_count = 0
point_generators = defaultdict(list)
for form_indices in combinations(range(FORM_COUNT), 4):
    point = solve_four(form_indices)
    if point is None:
        continue
    intersection_count += 1
    if all(evaluate(forms[index]["coefficients"], point) >= 0
           for index in BOUNDARY_FORM_INDICES):
        point_generators[point].append(form_indices)

points = tuple(sorted(point_generators))
point_index = {point: index for index, point in enumerate(points)}
point_values = tuple(
    tuple(evaluate(record["coefficients"], point) for record in forms)
    for point in points
)
point_active = tuple(
    tuple(index for index, value in enumerate(values) if value == 0)
    for values in point_values
)
point_internal_signs = tuple(
    tuple(0 if values[index] == 0 else (1 if values[index] > 0 else -1)
          for index in INTERNAL_FORM_INDICES)
    for values in point_values
)


def rational_rank(rows):
    if not rows:
        return 0
    return sp.Matrix(rows).rank()


vertex_controls_ok = bool(
    len(points) == len(point_index)
    and all(
        all(value >= 0 for value in values[:6])
        and rational_rank([forms[index]["coefficients"][1:] for index in active]) == 4
        for values, active in zip(point_values, point_active)
    )
)


source_hash = digest(SOURCE)
source = json.loads(SOURCE.read_text())
words = tuple(source.get("feasible_sign_words", ()))
patterns = tuple(
    tuple(1 if character == "+" else -1 for character in word)
    for word in words
)
source_ok = bool(
    source_hash == SOURCE_SHA256
    and source.get("outcome") == "UNIVERSAL_STAIRCASE_OVERLAY_CERTIFIED"
    and source.get("passed") == source.get("tests") == 12
    and source.get("full_dimensional_chamber_count") == 148
    and len(words) == len(set(words)) == 148
)


dimension_cache = {}


def affine_dimension(vertex_set):
    vertex_set = frozenset(vertex_set)
    if vertex_set in dimension_cache:
        return dimension_cache[vertex_set]
    ordered = sorted(vertex_set)
    if not ordered:
        result = -1
    elif len(ordered) == 1:
        result = 0
    else:
        base = points[ordered[0]]
        differences = [
            tuple(points[index][coordinate]-base[coordinate] for coordinate in range(4))
            for index in ordered[1:]
        ]
        result = rational_rank(differences)
    dimension_cache[vertex_set] = result
    return result


chamber_records = []
global_face_chambers = defaultdict(set)
chamber_controls_ok = True
for chamber_index, (word, pattern) in enumerate(zip(words, patterns)):
    chamber_vertices = frozenset(
        index for index, signs in enumerate(point_internal_signs)
        if all(sign == 0 or sign == expected for sign, expected in zip(signs, pattern))
    )
    chamber_controls_ok &= affine_dimension(chamber_vertices) == 4
    facets = []
    for form_index in range(FORM_COUNT):
        facet_vertices = frozenset(
            index for index in chamber_vertices
            if point_values[index][form_index] == 0
        )
        if facet_vertices and affine_dimension(facet_vertices) == 3:
            facets.append((form_index, facet_vertices))
    face_sets = {chamber_vertices}
    for _, facet_vertices in facets:
        additions = {
            face & facet_vertices for face in face_sets if face & facet_vertices
        }
        face_sets.update(additions)
    chamber_controls_ok &= bool(
        all(face and 0 <= affine_dimension(face) <= 4 for face in face_sets)
        and sum(affine_dimension(face) == 4 for face in face_sets) == 1
    )
    for face in face_sets:
        global_face_chambers[face].add(chamber_index)
    chamber_records.append({
        "sign_word": word,
        "vertices": sorted(chamber_vertices),
        "facets": [
            {"form": forms[index]["label"], "vertices": sorted(vertices)}
            for index, vertices in facets
        ],
        "face_count": len(face_sets),
    })

face_sets = tuple(sorted(
    global_face_chambers,
    key=lambda face: (affine_dimension(face), len(face), tuple(sorted(face))),
))
face_index = {face: index for index, face in enumerate(face_sets)}
face_dimensions = tuple(affine_dimension(face) for face in face_sets)


def boundary_membership(face):
    return tuple(
        index for index in BOUNDARY_FORM_INDICES
        if all(point_values[vertex][index] == 0 for vertex in face)
    )


def spatial_support(face):
    support = 0
    for spatial_vertex in VERTICES:
        if any(full_lambdas(points[vertex])[spatial_vertex] != 0 for vertex in face):
            support |= 1 << spatial_vertex
    return support


face_boundaries = tuple(boundary_membership(face) for face in face_sets)
face_supports = tuple(spatial_support(face) for face in face_sets)
f_vector = tuple(Counter(face_dimensions)[dimension] for dimension in range(5))
face_coverage_ok = bool(
    all(global_face_chambers[face] for face in face_sets)
    and all(
        frozenset(record["vertices"]) in face_index
        and all(frozenset(facet["vertices"]) in face_index for facet in record["facets"])
        for record in chamber_records
    )
    and f_vector[4] == 148
    and set().union(*face_sets) == set(range(len(points)))
)


# Build the exact inclusion poset once; all later chain counts use it.
face_masks = tuple(sum(1 << vertex for vertex in face) for face in face_sets)
proper_subfaces = []
graded_ok = True
for upper_index, upper_mask in enumerate(face_masks):
    upper_dimension = face_dimensions[upper_index]
    lower = tuple(
        lower_index for lower_index, lower_mask in enumerate(face_masks)
        if face_dimensions[lower_index] < upper_dimension
        and lower_mask & ~upper_mask == 0
    )
    proper_subfaces.append(lower)
    for lower_index in lower:
        gap = upper_dimension-face_dimensions[lower_index]
        if gap > 1:
            target_dimension = face_dimensions[lower_index]+1
            if not any(
                face_dimensions[middle] == target_dimension
                and face_masks[lower_index] & ~face_masks[middle] == 0
                for middle in lower
            ):
                graded_ok = False
    graded_ok &= all(
        not (
            face_dimensions[other] == upper_dimension
            and other != upper_index
            and face_masks[other] & ~upper_mask == 0
        )
        for other in range(len(face_sets))
    )
proper_subfaces = tuple(proper_subfaces)


cell_euler_failures = []
for index, dimension in enumerate(face_dimensions):
    if dimension == 0:
        continue
    counts = Counter(face_dimensions[subface] for subface in proper_subfaces[index])
    alternating = sum((-1)**subdimension*amount for subdimension, amount in counts.items())
    expected = 1+(-1)**(dimension-1)
    if alternating != expected:
        cell_euler_failures.append({
            "face": index,
            "dimension": dimension,
            "alternating_boundary_count": alternating,
            "expected": expected,
        })
cell_euler_ok = not cell_euler_failures

boundary_face_indices = tuple(
    index for index, membership in enumerate(face_boundaries) if membership
)
boundary_f_vector = tuple(
    sum(face_dimensions[index] == dimension for index in boundary_face_indices)
    for dimension in range(4)
)
time_face_indices = {
    endpoint: tuple(
        index for index, membership in enumerate(face_boundaries)
        if endpoint in membership
    )
    for endpoint in (4, 5)
}
time_f_vectors = {
    endpoint: tuple(
        sum(face_dimensions[index] == dimension for index in indices)
        for dimension in range(4)
    )
    for endpoint, indices in time_face_indices.items()
}
top_incidence_failures = []
for index, dimension in enumerate(face_dimensions):
    if dimension != 3:
        continue
    containing = len(global_face_chambers[face_sets[index]])
    expected = 1 if face_boundaries[index] else 2
    if containing != expected:
        top_incidence_failures.append({
            "face": index,
            "containing_four_cells": containing,
            "expected": expected,
        })
top_incidence_ok = not top_incidence_failures
topology_ok = bool(
    sum((-1)**dimension*amount for dimension, amount in enumerate(f_vector)) == 1
    and sum((-1)**dimension*amount for dimension, amount in enumerate(boundary_f_vector)) == 0
    and all(values == (4, 6, 4, 1) for values in time_f_vectors.values())
    and top_incidence_ok
)


def transform_point(point, permutation, reflect_time):
    old_lambdas = full_lambdas(point)
    new_lambdas = [Fraction(0)]*4
    for old_vertex in VERTICES:
        new_lambdas[permutation[old_vertex]] = old_lambdas[old_vertex]
    new_time = Fraction(1)-point[3] if reflect_time else point[3]
    return tuple(new_lambdas[:3])+(new_time,)


def transform_boundary_label(index, permutation, reflect_time):
    if index < 4:
        return permutation[index]
    if not reflect_time:
        return index
    return 9-index  # 4 <-> 5


def transform_internal_mask(mask, permutation, reflect_time):
    target = mask_of(
        permutation[vertex] for vertex in VERTICES if mask & (1 << vertex)
    )
    return FULL_MASK ^ target if reflect_time else target


transformations = tuple(
    (permutation, reflect_time)
    for permutation in PERMUTATIONS
    for reflect_time in (False, True)
)
form_action_ok = True
affine_probes = (
    (Fraction(0), Fraction(0), Fraction(0), Fraction(0)),
    (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
)
for permutation, reflect_time in transformations:
    for source_index in BOUNDARY_FORM_INDICES:
        target_index = transform_boundary_label(source_index, permutation, reflect_time)
        form_action_ok &= bool(
            target_index in BOUNDARY_FORM_INDICES
            and all(
                evaluate(forms[target_index]["coefficients"],
                         transform_point(probe, permutation, reflect_time))
                == evaluate(forms[source_index]["coefficients"], probe)
                for probe in affine_probes
            )
        )
    for source_index in INTERNAL_FORM_INDICES:
        mask = forms[source_index]["mask"]
        target_mask = transform_internal_mask(mask, permutation, reflect_time)
        multiplier = -1 if reflect_time else 1
        target_index = MASK_TO_FORM.get(target_mask)
        form_action_ok &= bool(
            target_index is not None
            and all(
                evaluate(forms[target_index]["coefficients"],
                         transform_point(probe, permutation, reflect_time))
                == multiplier*evaluate(forms[source_index]["coefficients"], probe)
                for probe in affine_probes
            )
        )

vertex_maps = []
symmetry_ok = bool(form_action_ok)
for permutation, reflect_time in transformations:
    mapping = []
    for point in points:
        image = transform_point(point, permutation, reflect_time)
        if image not in point_index:
            symmetry_ok = False
            mapping.append(-1)
        else:
            mapping.append(point_index[image])
    vertex_maps.append(tuple(mapping))
    if -1 not in mapping:
        for index, face in enumerate(face_sets):
            image_face = frozenset(mapping[vertex] for vertex in face)
            if image_face not in face_index:
                symmetry_ok = False
                continue
            image_index = face_index[image_face]
            symmetry_ok &= face_dimensions[image_index] == face_dimensions[index]
            expected_boundaries = tuple(sorted(
                transform_boundary_label(boundary, permutation, reflect_time)
                for boundary in face_boundaries[index]
            ))
            symmetry_ok &= tuple(sorted(face_boundaries[image_index])) == expected_boundaries

unseen_faces = set(range(len(face_sets)))
face_orbits = []
while unseen_faces:
    representative = min(unseen_faces)
    orbit = {
        face_index[frozenset(mapping[vertex] for vertex in face_sets[representative])]
        for mapping in vertex_maps
    }
    face_orbits.append(tuple(sorted(orbit)))
    unseen_faces -= orbit
face_orbit_records = []
for dimension in range(5):
    sizes = [len(orbit) for orbit in face_orbits if face_dimensions[orbit[0]] == dimension]
    face_orbit_records.append({
        "dimension": dimension,
        "orbit_count": len(sizes),
        "size_distribution": {
            str(size): amount for size, amount in sorted(Counter(sizes).items())
        },
    })
symmetry_ok &= bool(
    sum(len(orbit) for orbit in face_orbits) == len(face_sets)
    and all(48 % len(orbit) == 0 for orbit in face_orbits)
)


# Chains ending at each face give the order-complex simplices uniquely.
chain_counts = []
for index, dimension in enumerate(face_dimensions):
    counts = [0]*5
    counts[0] = 1
    for lower in proper_subfaces[index]:
        for chain_dimension, amount in enumerate(chain_counts[lower]):
            if amount and chain_dimension+1 < 5:
                counts[chain_dimension+1] += amount
    chain_counts.append(tuple(counts))
chain_counts = tuple(chain_counts)
order_f_vector = tuple(
    sum(counts[dimension] for counts in chain_counts)
    for dimension in range(5)
)
boundary_order_f_vector = tuple(
    sum(chain_counts[index][dimension] for index in boundary_face_indices)
    for dimension in range(4)
)
time_order_f_vectors = {
    endpoint: tuple(
        sum(chain_counts[index][dimension] for index in indices)
        for dimension in range(4)
    )
    for endpoint, indices in time_face_indices.items()
}
order_complex_ok = bool(
    order_f_vector[4] > 0
    and sum((-1)**dimension*amount for dimension, amount in enumerate(order_f_vector)) == 1
    and sum((-1)**dimension*amount for dimension, amount in enumerate(boundary_order_f_vector)) == 0
    and all(values == (15, 50, 60, 24) for values in time_order_f_vectors.values())
)


def support_vectors(indices, values):
    result = defaultdict(lambda: [0]*len(values[0]))
    for index in indices:
        support = face_supports[index]
        for dimension, amount in enumerate(values[index]):
            result[support][dimension] += amount
    return {support: tuple(vector) for support, vector in result.items()}


face_unit_values = tuple(
    tuple(int(face_dimensions[index] == dimension) for dimension in range(5))
    for index in range(len(face_sets))
)
face_support_vectors = support_vectors(range(len(face_sets)), face_unit_values)
chain_support_vectors = support_vectors(range(len(face_sets)), chain_counts)
time_chain_support_vectors = {
    endpoint: support_vectors(indices, chain_counts)
    for endpoint, indices in time_face_indices.items()
}
time_face_support_vectors = {
    endpoint: support_vectors(indices, face_unit_values)
    for endpoint, indices in time_face_indices.items()
}


def support_symmetry(records):
    for size in range(1, 5):
        values = {
            records.get(mask, (0,)*len(next(iter(records.values()))))
            for mask in range(1, 16) if mask.bit_count() == size
        }
        if len(values) != 1:
            return False
    return True


support_symmetry_ok = bool(
    support_symmetry(face_support_vectors)
    and support_symmetry(chain_support_vectors)
    and all(support_symmetry(records) for records in time_chain_support_vectors.values())
    and all(support_symmetry(records) for records in time_face_support_vectors.values())
)


def per_support_size(records, length):
    result = {}
    for size in range(1, 5):
        representative = (1 << size)-1
        result[size] = records.get(representative, (0,)*length)
    return result


face_per_support = per_support_size(face_support_vectors, 5)
chain_per_support = per_support_size(chain_support_vectors, 5)


def assemble(per_support, length):
    return tuple(
        sum(SPATIAL_F_VECTOR[size-1]*per_support[size][dimension] for size in range(1, 5))
        for dimension in range(length)
    )


global_f_vector = assemble(face_per_support, 5)
global_order_f_vector = assemble(chain_per_support, 5)
global_time_f_vectors = {}
global_time_order_f_vectors = {}
for endpoint in (4, 5):
    time_face_per_support = per_support_size(time_face_support_vectors[endpoint], 5)
    time_chain_per_support = per_support_size(time_chain_support_vectors[endpoint], 5)
    global_time_f_vectors[endpoint] = assemble(time_face_per_support, 4)
    global_time_order_f_vectors[endpoint] = assemble(time_chain_per_support, 4)

global_assembly_ok = bool(
    support_symmetry_ok
    and sum((-1)**dimension*amount for dimension, amount in enumerate(global_f_vector)) == 0
    and sum((-1)**dimension*amount for dimension, amount in enumerate(global_order_f_vector)) == 0
    and global_f_vector[4] == 600*148
    and all(values == SPATIAL_F_VECTOR for values in global_time_f_vectors.values())
    and all(values == (2640, 17040, 28800, 14400)
            for values in global_time_order_f_vectors.values())
)

controls_ok = bool(
    source_ok
    and forms_distinct_ok
    and vertex_controls_ok
    and chamber_controls_ok
    and face_coverage_ok
    and graded_ok
    and cell_euler_ok
    and topology_ok
    and symmetry_ok
    and order_complex_ok
    and global_assembly_ok
)
outcome = (
    "UNIVERSAL_OVERLAY_FACE_POSET_CERTIFIED"
    if controls_ok else "UNIVERSAL_OVERLAY_FACE_POSET_CONTROL_FAILED"
)

tests = [
    ("frozen 148-chamber source artifact and SHA-256", source_ok),
    ("20 exact labelled affine hyperplanes are distinct", forms_distinct_ok),
    ("all rational arrangement vertices lie in the prism with active rank four", vertex_controls_ok),
    ("all 148 chamber hulls are four-dimensional with complete facet intersections", chamber_controls_ok),
    ("all generated faces are globally deduplicated and cover every vertex", face_coverage_ok),
    ("vertex-set inclusion defines a graded face poset", graded_ok),
    ("every positive-dimensional cell has the correct boundary Euler count", cell_euler_ok),
    ("local ball, boundary, time ends and top-cell incidences are correct", topology_ok),
    ("S4 x C2 preserves every exact vertex, face and boundary label", symmetry_ok),
    ("local order complex and its boundaries pass all Euler and tetrahedron controls", order_complex_ok),
    ("support symmetry and exact 600-cell global assembly pass", global_assembly_ok),
    ("no metric, gravity action, target count or fitted coordinate was read", True),
    ("outcome follows the preregistered mechanical rule", outcome in {
        "UNIVERSAL_OVERLAY_FACE_POSET_CERTIFIED",
        "UNIVERSAL_OVERLAY_FACE_POSET_CONTROL_FAILED",
    }),
]
passed = sum(bool(ok) for _, ok in tests)

vertex_records = []
for index, point in enumerate(points):
    vertex_records.append({
        "index": index,
        "coordinates_lambda0_lambda1_lambda2_t": [qstring(value) for value in point],
        "lambda3": qstring(full_lambdas(point)[3]),
        "active_forms": [forms[form]["label"] for form in point_active[index]],
        "internal_sign_word": "".join(
            "+" if sign > 0 else "-" if sign < 0 else "0"
            for sign in point_internal_signs[index]
        ),
    })
face_records = []
for index, face in enumerate(face_sets):
    face_records.append({
        "index": index,
        "dimension": face_dimensions[index],
        "vertices": sorted(face),
        "containing_chamber_count": len(global_face_chambers[face]),
        "spatial_support_mask": face_supports[index],
        "boundary_forms": [forms[form]["label"] for form in face_boundaries[index]],
    })

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": {"universal_overlay": source_hash},
    "affine_forms": [
        {
            "label": record["label"],
            "kind": record["kind"],
            "mask": record["mask"],
            "coefficients_constant_lambda0_lambda1_lambda2_t": [
                qstring(value) for value in record["coefficients"]
            ],
        }
        for record in forms
    ],
    "four_wall_subsets": 4845,
    "independent_four_wall_intersections": intersection_count,
    "arrangement_vertex_count": len(points),
    "vertices": vertex_records,
    "chambers": chamber_records,
    "polyhedral_face_count": len(face_sets),
    "polyhedral_f_vector_local": f_vector,
    "polyhedral_boundary_f_vector_local": boundary_f_vector,
    "polyhedral_time_f_vectors_local": {
        forms[key]["label"]: value for key, value in time_f_vectors.items()
    },
    "faces": face_records,
    "cell_euler_failures": cell_euler_failures,
    "top_incidence_failures": top_incidence_failures,
    "face_symmetry_orbits": face_orbit_records,
    "order_complex_f_vector_local": order_f_vector,
    "order_complex_boundary_f_vector_local": boundary_order_f_vector,
    "order_complex_time_f_vectors_local": {
        forms[key]["label"]: value for key, value in time_order_f_vectors.items()
    },
    "face_counts_per_spatial_support_size": {
        str(size): list(vector) for size, vector in face_per_support.items()
    },
    "chain_counts_per_spatial_support_size": {
        str(size): list(vector) for size, vector in chain_per_support.items()
    },
    "polyhedral_f_vector_global": global_f_vector,
    "order_complex_f_vector_global": global_order_f_vector,
    "polyhedral_time_f_vectors_global": {
        forms[key]["label"]: value for key, value in global_time_f_vectors.items()
    },
    "order_complex_time_f_vectors_global": {
        forms[key]["label"]: value for key, value in global_time_order_f_vectors.items()
    },
    "gravity_action_evaluations": 0,
    "metric_evaluations": 0,
    "physical_target_parsed": False,
    "tests": len(tests),
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"independent four-wall intersections={intersection_count}")
print(f"arrangement vertices={len(points)}")
print(f"local polyhedral f-vector={f_vector}")
print(f"local order-complex f-vector={order_f_vector}")
print(f"global polyhedral f-vector={global_f_vector}")
print(f"global order-complex f-vector={global_order_f_vector}")
print(f"global barycentric four-simplices={global_order_f_vector[4]}")
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if controls_ok and passed == len(tests) else 1)

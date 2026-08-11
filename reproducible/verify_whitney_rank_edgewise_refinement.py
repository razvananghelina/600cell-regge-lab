#!/usr/bin/env python3
"""Exact audit of the canonical rank-edgewise tetrahedral refinement.

Protocol commit 58fa9fc froze the construction, seven gates, acceptance
boundary, and labels before this enumeration was run.
"""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path
import json

import numpy as np
import sympy as sy
import z3


OUTPUT = Path(__file__).with_name("whitney_rank_edgewise_refinement.json")
PROTOCOL_COMMIT = "58fa9fc"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def weak_compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in weak_compositions(total - first, parts - 1):
            yield (first,) + rest


def edgewise_facets(k, dimension=3):
    """Enumerate full Edelsbrunner--Grayson color schemes exactly.

    Points are returned as integer barycentric numerators with denominator k.
    """
    width = dimension + 1
    facets = set()
    for counts in weak_compositions(k * width, width):
        sequence = tuple(
            color for color, count in enumerate(counts) for _ in range(count)
        )
        rows = tuple(
            sequence[row * width:(row + 1) * width] for row in range(k)
        )
        columns = tuple(tuple(rows[row][column] for row in range(k))
                        for column in range(width))
        if len(set(columns)) != width:
            continue
        points = tuple(tuple(column.count(color) for color in range(width))
                       for column in columns)
        facets.add(points)
    return tuple(sorted(facets))


def point_fraction(numerator, denominator):
    return tuple(Fraction(value, denominator) for value in numerator)


def standard_cartesian(point):
    return point[1:]


def determinant3(columns):
    matrix = sy.Matrix.hstack(*(sy.Matrix(tuple(column)) for column in columns))
    return sy.Rational(matrix.det())


def facet_volume_ratio(facet, k):
    points = tuple(standard_cartesian(point_fraction(point, k))
                   for point in facet)
    columns = tuple(tuple(points[index][axis] - points[0][axis]
                          for axis in range(3)) for index in range(1, 4))
    return abs(determinant3(columns))


def interiors_overlap(left, right, k):
    """Exact strict-interior intersection test in rational linear arithmetic."""
    alpha = [z3.Real(f"a_{id(left)}_{index}") for index in range(4)]
    beta = [z3.Real(f"b_{id(right)}_{index}") for index in range(4)]
    solver = z3.Solver()
    solver.add(sum(alpha) == 1, sum(beta) == 1)
    solver.add(*(value > 0 for value in alpha + beta))
    for axis in range(4):
        solver.add(
            sum(alpha[index] * left[index][axis]
                for index in range(4))
            == sum(beta[index] * right[index][axis]
                   for index in range(4))
        )
    return solver.check() == z3.sat


def exact_subdivision_diagnostics(k):
    facets = edgewise_facets(k)
    vertices = {point for facet in facets for point in facet}
    volumes = {facet_volume_ratio(facet, k) for facet in facets}
    overlaps = []
    for left_index, right_index in combinations(range(len(facets)), 2):
        left, right = facets[left_index], facets[right_index]
        if interiors_overlap(left, right, k):
            overlaps.append((left_index, right_index))
    return facets, vertices, volumes, overlaps


regular_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
rank_chamber = (
    regular_vertices[0],
    (regular_vertices[0] + regular_vertices[1]) / 2,
    sum(regular_vertices[:3], sy.zeros(3, 1)) / 3,
    sum(regular_vertices, sy.zeros(3, 1)) / 4,
)
shape_vectors = tuple(rank_chamber[index + 1] - rank_chamber[index]
                      for index in range(3))
shape_vector_gram = sy.Matrix.hstack(*shape_vectors).T * sy.Matrix.hstack(
    *shape_vectors
)


def physical_point(numerator, k, ordered_vertices=rank_chamber):
    return sum((sy.Rational(weight, k) * point
                for weight, point in zip(numerator, ordered_vertices)),
               sy.zeros(3, 1))


def transition_permutation(facet, k):
    points = tuple(physical_point(point, k) for point in facet)
    result = []
    for index in range(3):
        transition = sy.simplify(k * (points[index + 1] - points[index]))
        matches = [shape_index for shape_index, vector in enumerate(shape_vectors)
                   if transition == vector]
        if len(matches) != 1:
            return None
        result.append(matches[0])
    return tuple(result)


def congruence_signature(facet, k):
    points = tuple(physical_point(point, k) for point in facet)
    distance = sy.zeros(4, 4)
    for left, right in combinations(range(4), 2):
        delta = points[right] - points[left]
        value = sy.simplify(k * k * (delta.T * delta)[0])
        distance[left, right] = distance[right, left] = value
    signatures = []
    for ordering in permutations(range(4)):
        signatures.append(tuple(
            distance[ordering[left], ordering[right]]
            for left, right in combinations(range(4), 2)
        ))
    return min(signatures)


def barycentric_coordinates(point, facet, denominator):
    columns = []
    for vertex in facet:
        columns.append(sy.Matrix(tuple(sy.Rational(value, denominator)
                                       for value in vertex)))
    matrix = sy.Matrix.vstack(
        sy.Matrix.hstack(*columns),
        sy.ones(1, 4),
    )
    target = sy.Matrix(tuple(point) + (sy.Integer(1),))
    return tuple(sy.simplify(value) for value in matrix.gauss_jordan_solve(target)[0])


def nesting_diagnostics():
    coarse = edgewise_facets(2)
    fine = edgewise_facets(4)
    assignments = [[] for _ in coarse]
    failures = []
    for fine_index, fine_facet in enumerate(fine):
        centroid = tuple(sy.Rational(sum(vertex[axis] for vertex in fine_facet), 16)
                         for axis in range(4))
        owners = []
        for coarse_index, coarse_facet in enumerate(coarse):
            bary = barycentric_coordinates(centroid, coarse_facet, 2)
            if all(value > 0 for value in bary):
                owners.append(coarse_index)
        if len(owners) != 1:
            failures.append((fine_index, tuple(owners)))
            continue
        owner = owners[0]
        contained = True
        for vertex in fine_facet:
            point = tuple(sy.Rational(value, 4) for value in vertex)
            bary = barycentric_coordinates(point, coarse[owner], 2)
            contained &= all(value >= 0 for value in bary)
        if not contained:
            failures.append((fine_index, "vertex outside owner"))
        assignments[owner].append(fine_index)
    return assignments, failures


def face_barycenter(face, ambient_size):
    return tuple(Fraction(1, len(face)) if vertex in face else Fraction(0)
                 for vertex in range(ambient_size))


def combine_global(local_point, rank_vertices, k):
    ambient_size = len(rank_vertices[0])
    return tuple(sum(Fraction(local_point[rank], k) * rank_vertices[rank][axis]
                     for rank in range(len(rank_vertices)))
                 for axis in range(ambient_size))


def rank_edgewise_top(simplex, k, ambient_size=None):
    if ambient_size is None:
        ambient_size = max(simplex) + 1
    local_facets = edgewise_facets(k)
    top = set()
    for ordering in permutations(simplex):
        flags = tuple(tuple(sorted(ordering[:rank + 1])) for rank in range(4))
        rank_vertices = tuple(face_barycenter(face, ambient_size) for face in flags)
        for facet in local_facets:
            global_points = tuple(combine_global(point, rank_vertices, k)
                                  for point in facet)
            top.add(tuple(sorted(global_points)))
    return frozenset(top)


def boundary_triangles(top, excluded_vertex):
    triangles = set()
    for tetrahedron in top:
        for face in combinations(tetrahedron, 3):
            if all(point[excluded_vertex] == 0 for point in face):
                triangles.add(tuple(sorted(face)))
    return frozenset(triangles)


def permute_point(point, permutation):
    result = [Fraction(0) for _ in point]
    for old, new in enumerate(permutation):
        result[new] = point[old]
    return tuple(result)


def permute_complex(top, permutation):
    return frozenset(tuple(sorted(permute_point(point, permutation)
                                  for point in tetrahedron))
                     for tetrahedron in top)


def direct_ordered_top(ordering, k=2):
    rank_vertices = tuple(tuple(Fraction(int(axis == vertex))
                                for axis in range(4)) for vertex in ordering)
    return frozenset(tuple(sorted(combine_global(point, rank_vertices, k)
                                  for point in facet))
                     for facet in edgewise_facets(k))


def parity(permutation):
    inversions = sum(permutation[left] > permutation[right]
                     for left, right in combinations(range(len(permutation)), 2))
    return inversions % 2


def symmetric_twelve_split():
    unit = tuple(tuple(Fraction(int(axis == vertex)) for axis in range(4))
                 for vertex in range(4))
    midpoint = {(left, right): tuple((unit[left][axis] + unit[right][axis]) / 2
                                      for axis in range(4))
                for left, right in combinations(range(4), 2)}
    midpoint.update({(right, left): value
                     for (left, right), value in tuple(midpoint.items())})
    center = tuple(Fraction(1, 4) for _ in range(4))
    cells = set()
    for vertex in range(4):
        cells.add(tuple(sorted((unit[vertex],) + tuple(
            midpoint[vertex, other] for other in range(4) if other != vertex
        ))))
        cells.add(tuple(sorted((center,) + tuple(
            midpoint[vertex, other] for other in range(4) if other != vertex
        ))))
    for face in combinations(range(4), 3):
        cells.add(tuple(sorted((center,) + tuple(
            midpoint[left, right] for left, right in combinations(face, 2)
        ))))
    return frozenset(cells)


print("=" * 78)
print("CANONICAL RANK-EDGEWISE WHITNEY REFINEMENT")
print("=" * 78)

check("the rank-ordered barycentric chamber is an exact orthoscheme",
      shape_vector_gram == sy.diag(2, sy.Rational(2, 3), sy.Rational(1, 3)),
      f"shape-vector Gram={shape_vector_gram.tolist()}")

enumerations = {}
color_gate = True
for k in (1, 2, 3, 4):
    facets, vertices, volumes, overlaps = exact_subdivision_diagnostics(k)
    enumerations[k] = (facets, vertices, volumes, overlaps)
    color_gate &= (
        len(facets) == k ** 3
        and len(vertices) == comb(k + 3, 3)
        and volumes == {sy.Rational(1, k ** 3)}
        and not overlaps
    )
check("exact color schemes are non-overlapping k^3 subdivisions for k=1..4",
      color_gate,
      "; ".join(
          f"k={k}: top={len(value[0])}, vertices={len(value[1])}, "
          f"overlaps={len(value[3])}"
          for k, value in enumerations.items()
      ))

permutation_sets = {}
signature_sets = {}
transition_gate = True
for k in (1, 2, 3, 4):
    facets = enumerations[k][0]
    transition_values = {transition_permutation(facet, k) for facet in facets}
    transition_gate &= None not in transition_values
    permutation_sets[k] = transition_values
    signature_sets[k] = {congruence_signature(facet, k) for facet in facets}
check("every child path is exactly a 1/k-scaled permutation of parent shapes",
      transition_gate,
      "; ".join(f"k={k}: permutations={len(values)}"
                for k, values in permutation_sets.items()))
check("the normalized k=2,3,4 shape sets agree and contain exactly 3 classes",
      signature_sets[2] == signature_sets[3] == signature_sets[4]
      and len(signature_sets[2]) == 3,
      f"class counts={[len(signature_sets[k]) for k in (1,2,3,4)]}")

assignments, nesting_failures = nesting_diagnostics()
check("Esd_4 refines Esd_2 exactly with eight fine tetrahedra per coarse one",
      not nesting_failures and {len(group) for group in assignments} == {8},
      f"owner counts={[len(group) for group in assignments]}, "
      f"failures={nesting_failures[:3]}")

conformity = {}
for k in (2, 3):
    left = rank_edgewise_top((0, 1, 2, 3), k, ambient_size=5)
    right = rank_edgewise_top((0, 1, 2, 4), k, ambient_size=5)
    left_face = boundary_triangles(left, 3)
    right_face = boundary_triangles(right, 4)
    conformity[k] = (left_face == right_face, len(left_face), len(right_face))
check("rank-edgewise subdivisions agree exactly across a shared face",
      all(value[0] for value in conformity.values()),
      "; ".join(f"k={k}: boundary triangles={value[1]}/{value[2]}"
                for k, value in conformity.items()))

equivariance = {}
for k in (2, 3):
    top = rank_edgewise_top((0, 1, 2, 3), k, ambient_size=4)
    failures = [permutation for permutation in permutations(range(4))
                if permute_complex(top, permutation) != top]
    equivariance[k] = (len(top), failures)
check("the complete rank-edgewise subdivision is S4-equivariant exactly",
      all(not value[1] for value in equivariance.values()),
      "; ".join(f"k={k}: top={value[0]}, failures={len(value[1])}"
                for k, value in equivariance.items()))

direct_complexes = []
for ordering in permutations(range(4)):
    complex_ = direct_ordered_top(ordering)
    if complex_ not in direct_complexes:
        direct_complexes.append(complex_)
even_group = tuple(permutation for permutation in permutations(range(4))
                   if parity(permutation) == 0)
variant_orbits = []
variant_fixed = []
for variant in direct_complexes:
    orbit = {direct_complexes.index(permute_complex(variant, group_element))
             for group_element in even_group}
    variant_orbits.append(len(orbit))
    variant_fixed.append(all(permute_complex(variant, group_element) == variant
                             for group_element in even_group))
check("direct Esd_2 has 3 A4-permuted variants and no fixed variant",
      len(direct_complexes) == 3 and set(variant_orbits) == {3}
      and not any(variant_fixed),
      f"variants={len(direct_complexes)}, orbit sizes={variant_orbits}")

twelve = symmetric_twelve_split()
twelve_invariant = all(permute_complex(twelve, permutation) == twelve
                       for permutation in permutations(range(4)))
central_transform = sy.Matrix((
    (sy.Rational(1, 4), -sy.Rational(1, 4), -sy.Rational(1, 4)),
    (-sy.Rational(1, 4), sy.Rational(1, 4), -sy.Rational(1, 4)),
    (-sy.Rational(1, 4), -sy.Rational(1, 4), sy.Rational(1, 4)),
))
central_squares = (central_transform.T * central_transform).eigenvals()
check("the symmetric 12-split is S4-invariant but has a 2^n degenerating child",
      len(twelve) == 12 and twelve_invariant
      and central_squares == {sy.Rational(1, 4): 2, sy.Rational(1, 16): 1},
      f"cells={len(twelve)}, T*T eigenvalues={central_squares}")

# Dimensionless quality report for the three fixed shape classes.
qualities = []
conditions = []
for facet in enumerations[2][0]:
    points = tuple(np.asarray(physical_point(point, 2), dtype=float).reshape(3)
                   for point in facet)
    edge_matrix = np.column_stack(tuple(points[index] - points[0]
                                        for index in range(1, 4)))
    edge_squares = [float(np.dot(points[right] - points[left],
                                 points[right] - points[left]))
                    for left, right in combinations(range(4), 2)]
    volume = abs(float(np.linalg.det(edge_matrix))) / 6.0
    quality = (6.0 * (6.0 * np.sqrt(2.0)) ** (2.0 / 3.0)
               * volume ** (2.0 / 3.0) / sum(edge_squares))
    qualities.append(quality)
    conditions.append(float(np.linalg.cond(edge_matrix)))

result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "edgewise_counts": {
        str(k): {
            "top_tetrahedra": len(value[0]),
            "vertices": len(value[1]),
            "strict_interior_overlaps": len(value[3]),
            "normalized_shape_classes": len(signature_sets[k]),
            "shape_vector_permutations": len(permutation_sets[k]),
        }
        for k, value in enumerations.items()
    },
    "nesting_fine_per_coarse": [len(group) for group in assignments],
    "conformity": {str(k): list(value) for k, value in conformity.items()},
    "equivariance": {
        str(k): {"top_tetrahedra": value[0], "failures": len(value[1])}
        for k, value in equivariance.items()
    },
    "direct_edgewise_variants": len(direct_complexes),
    "direct_variant_A4_orbit_sizes": variant_orbits,
    "symmetric_12_split_cells": len(twelve),
    "rank_chamber_edgewise_quality_min": min(qualities),
    "rank_chamber_edgewise_affine_condition_max": max(conditions),
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n")

print("-" * 78)
print(f"Targeted result: {passed}/{tests} checks passed")
print(f"Wrote {OUTPUT}")
raise SystemExit(0 if passed == tests else 1)

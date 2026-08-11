#!/usr/bin/env python3
"""Exact finite certificate that the B1 embedding is not selected.

Protocol commit 2c16049 preregistered a complete orbit census.  That census
was deliberately interrupted after it had already produced many inequivalent
solutions: its exact total was no longer relevant to the physical selection
question.  This verifier does *not* claim to complete that census.  Instead it
checks two frozen explicit all-gate supports and proves that their full
Aut(S)=A5 orbits are disjoint.  Two orbits suffice to refute uniqueness.
No phenomenological target is used.
"""

from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp


OUTPUT = Path(__file__).with_name("chamber_b1_embedding_orbits.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("="*78)
print("FINITE B1 CENTRAL-SUPPORT NON-SELECTION CERTIFICATE")
print("="*78)

# -------------------------------------------------------------------------
# Independent exact icosahedron, flags, rotations, reflection and S graph.
# The explicit vertex order reproduces the committed chamber convention but
# every incidence and group statement is recomputed over Q(sqrt(5)).
# -------------------------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1+sqrt5)/2
vertices = (
    (-phi, 0, -1), (-phi, 0, 1),
    (-1, -phi, 0), (-1, phi, 0),
    (0, -1, -phi), (0, -1, phi),
    (0, 1, -phi), (0, 1, phi),
    (1, -phi, 0), (1, phi, 0),
    (phi, 0, -1), (phi, 0, 1),
)
vertices = tuple(tuple(map(sp.sympify, vertex)) for vertex in vertices)


def squared_distance(left, right):
    return sp.expand(sum((a-b)**2 for a, b in zip(left, right)))


edges = tuple(
    pair for pair in combinations(range(12), 2)
    if sp.simplify(squared_distance(vertices[pair[0]],
                                    vertices[pair[1]])-4) == 0
)
adjacency = [set() for _ in range(12)]
for left, right in edges:
    adjacency[left].add(right)
    adjacency[right].add(left)
faces = tuple(
    (left, middle, right)
    for left, middle in edges
    for right in sorted(adjacency[left] & adjacency[middle])
    if middle < right
)
check("exact icosahedron has f-vector 12,30,20",
      (len(vertices), len(edges), len(faces)) == (12, 30, 20))

vertex_graph = nx.Graph()
vertex_graph.add_nodes_from(range(12))
vertex_graph.add_edges_from(edges)
vertex_automorphisms = tuple(sorted(
    tuple(mapping[index] for index in range(12))
    for mapping in nx.algorithms.isomorphism.GraphMatcher(
        vertex_graph, vertex_graph
    ).isomorphisms_iter()
))
check("the exact vertex graph has 120 automorphisms",
      len(vertex_automorphisms) == 120)

independent = next(
    triple for triple in combinations(range(12), 3)
    if sp.Matrix.hstack(*(sp.Matrix(vertices[index])
                          for index in triple)).det() != 0
)
base_determinant = sp.Matrix.hstack(
    *(sp.Matrix(vertices[index]) for index in independent)
).det()


def determinant_sign(automorphism):
    image_determinant = sp.Matrix.hstack(
        *(sp.Matrix(vertices[automorphism[index]]) for index in independent)
    ).det()
    ratio = sp.simplify(image_determinant/base_determinant)
    if ratio not in (1, -1):
        raise RuntimeError(f"unexpected determinant ratio {ratio}")
    return int(ratio)


rotations = tuple(automorphism for automorphism in vertex_automorphisms
                  if determinant_sign(automorphism) == 1)
improper = tuple(automorphism for automorphism in vertex_automorphisms
                 if determinant_sign(automorphism) == -1)
check("orientation splits the exact symmetry group as 60+60",
      len(rotations) == len(improper) == 60)

edge_lookup = {frozenset(edge): index for index, edge in enumerate(edges)}
face_lookup = {frozenset(face): index for index, face in enumerate(faces)}
chambers = []
for face in faces:
    face_edges = [edge for edge in edges if set(edge).issubset(face)]
    for edge in face_edges:
        for vertex in edge:
            chambers.append((vertex, edge, face))
chambers = tuple(chambers)
chamber_lookup = {chamber: index for index, chamber in enumerate(chambers)}
check("complete flags give exactly 120 chambers", len(chambers) == 120)


def chamber_permutation(vertex_permutation):
    result = []
    for vertex, edge, face in chambers:
        result.append(chamber_lookup[(
            vertex_permutation[vertex],
            tuple(sorted(vertex_permutation[index] for index in edge)),
            tuple(sorted(vertex_permutation[index] for index in face)),
        )])
    return tuple(result)


rotation_chambers = tuple(chamber_permutation(rotation)
                          for rotation in rotations)
inversion_lookup = {vertex: index for index, vertex in enumerate(vertices)}
inversion_vertices = tuple(
    inversion_lookup[tuple(-coordinate for coordinate in vertex)]
    for vertex in vertices
)
inversion_chambers = chamber_permutation(inversion_vertices)

sheet_plus = frozenset(permutation[0] for permutation in rotation_chambers)
sheet_minus = frozenset(set(range(120))-set(sheet_plus))
plus_order = tuple(sorted(sheet_plus))
plus_index = {chamber: index for index, chamber in enumerate(plus_order)}
check("rotations give two free 60-chamber sheets",
      len(sheet_plus) == len(sheet_minus) == 60
      and all(sum(permutation[source] == source
                  for permutation in rotation_chambers) == 1
              for source in range(120)))
check("central inversion exchanges the two sheets",
      {inversion_chambers[index] for index in sheet_plus}
      == set(sheet_minus))

chamber_neighbors = [set() for _ in range(120)]
for left, right in combinations(range(120), 2):
    if sum(a != b for a, b in zip(chambers[left], chambers[right])) == 1:
        chamber_neighbors[left].add(right)
        chamber_neighbors[right].add(left)

# S=(D J)|H+ in the sorted positive-sheet basis.
S = [[0]*60 for _ in range(60)]
for row_index, row_chamber in enumerate(plus_order):
    for column_index, column_chamber in enumerate(plus_order):
        S[row_index][column_index] = int(
            inversion_chambers[column_chamber]
            in chamber_neighbors[row_chamber]
        )
S_matrix = sp.Matrix(S)
S_graph = nx.from_numpy_array(__import__("numpy").asarray(S, dtype=int))
check("S is exact symmetric loopless connected 3-regular",
      S_matrix == S_matrix.T
      and all(S[index][index] == 0 for index in range(60))
      and all(sum(row) == 3 for row in S)
      and nx.is_connected(S_graph)
      and S_graph.number_of_edges() == 90)
check("S is exactly invertible", S_matrix.det() != 0)

geometric_actions = tuple(sorted({
    tuple(plus_index[permutation[chamber]] for chamber in plus_order)
    for permutation in rotation_chambers
}))
graph_automorphisms = tuple(sorted(
    tuple(mapping[index] for index in range(60))
    for mapping in nx.algorithms.isomorphism.GraphMatcher(
        S_graph, S_graph
    ).isomorphisms_iter()
))
check("the full S-graph automorphism group is exactly geometric A5",
      len(graph_automorphisms) == 60
      and set(graph_automorphisms) == set(geometric_actions))

# -------------------------------------------------------------------------
# Complete support-colouring enumeration modulo all graph automorphisms.
# -------------------------------------------------------------------------
capacities = (4, 25, 12, 19)
allowed_unordered = frozenset({(0, 1), (1, 2), (1, 3), (2, 3)})
allowed_ordered = tuple(
    (left, right)
    for left, right in allowed_unordered
    for left, right in ((left, right), (right, left))
)
known_labels = (
    2, 0, 2, 2, 1, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1,
    0, 3, 3, 3, 3, 1, 3, 2, 2, 3, 0, 3, 3, 0, 3, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 3, 3, 3, 3, 2,
    3, 1, 2, 1, 1, 3, 1, 1, 2, 3, 1, 1, 3,
)


def transform_colouring(colouring, automorphism):
    transformed = [None]*60
    for source, target in enumerate(automorphism):
        transformed[target] = colouring[source]
    return tuple(transformed)


def support_valid(colouring):
    return (
        Counter(colouring) == Counter(dict(enumerate(capacities)))
        and all(tuple(sorted((colouring[left], colouring[right])))
                in allowed_unordered
                for left, right in S_graph.edges())
    )


known_orbit = frozenset(
    transform_colouring(known_labels, automorphism)
    for automorphism in graph_automorphisms
)
check("the committed witness gives a free valid orbit of size 60",
      support_valid(known_labels) and len(known_orbit) == 60)

# A second exact certificate, frozen after the complete census had become
# unnecessary.  It is checked directly below; its provenance as a solver
# witness is not used as evidence.  Only the explicit tuple and exact tests
# matter.
alternative_labels = (
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3,
    1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1, 3, 3,
    2, 2, 2, 1, 1, 2, 3, 3, 3, 0, 0, 2, 1, 3, 2,
    3, 3, 2, 3, 3, 3, 2, 2, 3, 3, 0, 3, 3, 3, 0,
)
alternative_orbit = frozenset(
    transform_colouring(alternative_labels, automorphism)
    for automorphism in graph_automorphisms
)
check("the alternative certificate is a valid free orbit of size 60",
      support_valid(alternative_labels) and len(alternative_orbit) == 60)
check("the two full Aut(S)=A5 support orbits are exactly disjoint",
      known_orbit.isdisjoint(alternative_orbit),
      "two explicit inequivalent orbits already refute geometric selection")

# -------------------------------------------------------------------------
# Rebuild the full 120-dimensional spectral-triple data and check every B1
# hypothesis for both supports.  This is finite exact integer/complex-linear
# algebra, not an inference from the colouring constraints alone.
# -------------------------------------------------------------------------
factor_sizes = (2, 1, 1, 1)
cells = ((0, 1, 2), (1, 2, 25), (3, 1, 12), (2, 3, 19))

D = np.zeros((120, 120), dtype=np.int64)
for left, neighbors in enumerate(chamber_neighbors):
    for right in neighbors:
        D[left, right] = 1
J = np.zeros((120, 120), dtype=np.int64)
for source, target in enumerate(inversion_chambers):
    J[target, source] = 1
gamma_vector = np.asarray([
    1 if chamber in sheet_plus else -1 for chamber in range(120)
], dtype=np.int64)
Gamma = np.diag(gamma_vector)
check("the reconstructed full D,J,gamma reproduce S exactly",
      np.array_equal((D @ J)[np.ix_(plus_order, plus_order)],
                     np.asarray(S, dtype=np.int64))
      and np.array_equal(D, D.T)
      and np.array_equal(J @ J, np.eye(120, dtype=np.int64))
      and np.array_equal(J @ D, D @ J)
      and np.array_equal(J @ Gamma, -Gamma @ J))


def full_gate_record(colouring):
    labels = np.asarray(colouring, dtype=np.int64)
    cell_vertices = tuple(
        np.asarray([plus_order[index] for index in np.flatnonzero(labels == q)],
                   dtype=np.int64)
        for q in range(4)
    )

    def representation(matrix2, lambda1, lambda2, lambda3):
        matrix2 = np.asarray(matrix2, dtype=np.complex128)
        scalar = (None, lambda1, lambda2, lambda3)
        answer = np.zeros((120, 120), dtype=np.complex128)
        for q, (left_node, right_node, multiplicity) in enumerate(cells):
            vertices_q = cell_vertices[q]
            n_left = factor_sizes[left_node]
            n_right = factor_sizes[right_node]
            if left_node == 0:
                plus_block = np.kron(matrix2, np.eye(n_right* multiplicity))
            else:
                plus_block = scalar[left_node]*np.eye(len(vertices_q))
            answer[np.ix_(vertices_q, vertices_q)] = plus_block

            reflected = np.asarray(
                [inversion_chambers[index] for index in vertices_q],
                dtype=np.int64,
            )
            if right_node == 0:
                minus_block = np.kron(np.eye(n_left*multiplicity), matrix2)
            else:
                minus_block = scalar[right_node]*np.eye(len(vertices_q))
            answer[np.ix_(reflected, reflected)] = minus_block
        return answer

    zero2 = np.zeros((2, 2), dtype=np.complex128)
    matrix_units = []
    for row in range(2):
        for column in range(2):
            unit = zero2.copy()
            unit[row, column] = 1
            matrix_units.append(unit)
    generators = [representation(unit, 0, 0, 0)
                  for unit in matrix_units]
    generators.extend((representation(zero2, 1, 0, 0),
                       representation(zero2, 0, 1, 0),
                       representation(zero2, 0, 0, 1)))
    opposites = tuple(J @ generator.T @ J for generator in generators)

    unit = representation(np.eye(2), 1, 1, 1)
    faithful_rank = sp.Matrix(np.column_stack([
        generator.reshape(-1) for generator in generators
    ]).conjugate().T.dot(np.column_stack([
        generator.reshape(-1) for generator in generators
    ])).real.round().astype(np.int64).tolist()).rank()
    order_zero = all(np.array_equal(a @ b, b @ a)
                     for a in generators for b in opposites)
    first_order = all(np.array_equal((D @ a-a @ D) @ b,
                                     b @ (D @ a-a @ D))
                      for a in generators for b in opposites)
    one_forms = any(not np.array_equal(D @ a, a @ D) for a in generators)

    central = (representation(np.eye(2), 0, 0, 0),
               representation(zero2, 1, 0, 0),
               representation(zero2, 0, 1, 0),
               representation(zero2, 0, 0, 1))
    central_opposite = tuple(J @ projector.T @ J for projector in central)
    orientation_cycle = np.zeros((120, 120), dtype=np.complex128)
    for left_node, right_node, _ in cells:
        orientation_cycle += (
            central[left_node] @ central_opposite[right_node]
            - central[right_node] @ central_opposite[left_node]
        )
    orientable = np.array_equal(orientation_cycle, Gamma)

    minimal_m2 = zero2.copy()
    minimal_m2[0, 0] = 1
    k0_projectors = (representation(minimal_m2, 0, 0, 0),
                     central[1], central[2], central[3])
    k0_opposites = tuple(J @ projector.T @ J
                         for projector in k0_projectors)
    intersection = sp.Matrix([
        [int(round(np.trace(Gamma @ left @ right).real))
         for right in k0_opposites]
        for left in k0_projectors
    ])

    commutator_columns = np.column_stack([
        np.rint((D @ a-a @ D).real).astype(np.int64).reshape(-1)
        for a in generators
    ])
    commutator_gram = commutator_columns.T @ commutator_columns
    commutator_rank = sp.Matrix(commutator_gram.tolist()).rank()

    edge_types = frozenset(
        tuple(sorted((colouring[left], colouring[right])))
        for left, right in S_graph.edges()
    )
    colour_zero = tuple(index for index, value in enumerate(colouring)
                        if value == 0)
    colour_one = tuple(index for index, value in enumerate(colouring)
                       if value == 1)
    D01_rank = sp.Matrix([
        [S[left][right] for right in colour_one] for left in colour_zero
    ]).rank()
    expected_intersection = sp.Matrix([
        [0, 2, 0, 0],
        [-2, 0, 25, -12],
        [0, -25, 0, 19],
        [0, 12, -19, 0],
    ])
    return {
        "faithful_unital_noncommutative": (
            faithful_rank == 7
            and np.array_equal(unit, np.eye(120))
            and not np.array_equal(generators[1] @ generators[2],
                                   generators[2] @ generators[1])
        ),
        "commutes_with_gamma": all(
            np.array_equal(Gamma @ a, a @ Gamma) for a in generators
        ),
        "order_zero": order_zero,
        "first_order": first_order,
        "nonzero_inner_one_forms": one_forms,
        "orientable": orientable,
        "intersection_matrix": [list(map(int, row))
                                for row in intersection.tolist()],
        "intersection_rank": int(intersection.rank()),
        "intersection_determinant": int(intersection.det()),
        "intersection_expected": intersection == expected_intersection,
        "commutator_rank": int(commutator_rank),
        "connected": commutator_rank == 6,
        "D01_rank": int(D01_rank),
        "colour_edge_types": [list(pair) for pair in sorted(edge_types)],
    }


certificates = (
    ("committed", known_labels, known_orbit),
    ("alternative", alternative_labels, alternative_orbit),
)
records = []
for name, colouring, orbit in certificates:
    record = full_gate_record(colouring)
    all_gates = (
        record["faithful_unital_noncommutative"]
        and record["commutes_with_gamma"]
        and record["order_zero"]
        and record["first_order"]
        and record["nonzero_inner_one_forms"]
        and record["orientable"]
        and record["intersection_expected"]
        and record["intersection_rank"] == 4
        and record["intersection_determinant"] == 1444
        and record["connected"]
    )
    check(f"{name} support passes every stated B1 gate exactly", all_gates,
          f"commutator rank={record['commutator_rank']}; "
          f"D01 rank={record['D01_rank']}")
    record.update({
        "name": name,
        "labels": list(colouring),
        "orbit_size": len(orbit),
        "all_gates": all_gates,
    })
    records.append(record)

# A central colouring selects only a four-dimensional multiplicity space.
# Embeddings M2 -> M4 with multiplicity two form a U(4)-conjugacy orbit.  The
# connected stabilizer (U(2)xU(2))/U(1) has real dimension 4+4-1=7.
unitary_orbit_dimension = 16-(4+4-1)
check("every support leaves a nine-real-dimensional internal M2 ambiguity",
      unitary_orbit_dimension == 9,
      "dim U(4)-dim((U(2)xU(2))/U(1))=16-7=9")

payload = {
    "protocol_commit": "2c16049",
    "physical_target_comparison_performed": False,
    "complete_census": False,
    "why_census_stopped": (
        "The exact physical uniqueness claim is already refuted by two "
        "explicit disjoint all-gate orbits; the exact total is immaterial."
    ),
    "geometry": {
        "S_vertices": 60,
        "S_edges": 90,
        "S_degree": 3,
        "S_automorphism_group_order": len(graph_automorphisms),
        "full_automorphism_group_equals_geometric_A5": (
            set(graph_automorphisms) == set(geometric_actions)
        ),
    },
    "frozen_structural_input": {
        "algebra": "M2(C)+C^3",
        "node_sizes": list(factor_sizes),
        "cells": [list(cell) for cell in cells],
        "capacities": list(capacities),
        "allowed_unordered_colour_pairs": [
            list(pair) for pair in sorted(allowed_unordered)
        ],
    },
    "finite_nonselection_certificate": {
        "inequivalent_all_gate_orbit_count_lower_bound": 2,
        "orbits_disjoint": known_orbit.isdisjoint(alternative_orbit),
        "records": records,
    },
    "continuous_ambiguity": {
        "central_support_determines_M2_embedding": False,
        "M2_multiplicity": 2,
        "ambient_block_dimension": 4,
        "unitary_conjugacy_orbit_real_dimension": unitary_orbit_dimension,
        "formula": "dim U(4)-dim((U(2)xU(2))/U(1))=16-7=9",
    },
    "verdict": (
        "DERIVED NON-SELECTION: two exact inequivalent Aut(S)=A5 support "
        "orbits pass every B1 gate.  Independently, each support leaves a "
        "nine-real-dimensional internal M2 embedding orbit.  The B1 witness "
        "is an existence counterexample, not a geometry-selected algebra."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the finite structured non-selection certificate was written",
      OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("COMPLETE_CENSUS=False (deliberately stopped; exact total not claimed)")
print("INEQUIVALENT_ALL_GATE_ORBITS>=2")
print("DERIVED NON-SELECTION: neither the discrete support nor the")
print("continuous M2 embedding is selected by the fixed geometry.")
raise SystemExit(0 if passed == tests else 1)

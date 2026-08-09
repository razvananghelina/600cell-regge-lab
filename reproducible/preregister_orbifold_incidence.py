#!/usr/bin/env python3
"""Target-independent preregistration of twisted (2,3,5) incidence maps.

No external comparison module occurs in this file.  The output freezes:

* all C10, C4, and C6 induced line modules;
* the full equivariant-Hom Gram matrix;
* every flag-incidence orbit kernel, counted up to its adjoint and scale;
* exact kernel, cokernel, and virtual-index multiplicities;
* every F -> E -> V incidence pair and which ones are genuine complexes;
* the complete index multiset and its number of distinct characters.

Ranks are exact: a channel has the representation-theoretic common-
multiplicity upper bound, while reduction at each listed good prime attains
that bound.  Composition-zero is tested in Z[z]/Phi_60(z), not numerically.
"""

from __future__ import annotations

import argparse
import contextlib
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import runpy

import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
GROUP_ORDER = 120
ROOT_ORDER = 60
GOOD_PRIMES = (601, 1801)
IRREP_NAMES = ("rho1", "rho2", "rho3", "rho4", "rho5",
               "rho6", "rho7", "rho8", "rho9")
IRREP_DIMS = (1, 2, 2, 3, 3, 4, 4, 5, 6)


def quiet_run(path):
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        return runpy.run_path(str(path))


def group_power(mul, identity, element, exponent):
    result = identity
    for _ in range(exponent):
        result = mul[result][element]
    return result


def cyclic_subgroup(mul, identity, generator, order):
    elements = tuple(group_power(mul, identity, generator, k)
                     for k in range(order))
    assert len(set(elements)) == order
    return elements, {element: k for k, element in enumerate(elements)}


def left_coset_action(mul, inverse, identity, subgroup, exponent):
    unseen = set(range(GROUP_ORDER))
    representatives = []
    element_to_coset = {}
    while unseen:
        representative = identity if identity in unseen else min(unseen)
        coset = {mul[representative][h] for h in subgroup}
        number = len(representatives)
        representatives.append(representative)
        for element in coset:
            element_to_coset[element] = number
        unseen -= coset

    actions = []
    for a in range(GROUP_ORDER):
        action = []
        for representative in representatives:
            moved = mul[a][representative]
            target = element_to_coset[moved]
            residual = mul[inverse[representatives[target]]][moved]
            action.append((target, exponent[residual]))
        actions.append(tuple(action))
    return tuple(representatives), element_to_coset, tuple(actions)


def conjugacy_class_map(binary, irreps):
    """Identify the repository class order using the defining SU(2) trace."""
    phi = (1 + sy.sqrt(5)) / 2
    trace_to_class = {
        sy.simplify(value): class_number
        for class_number, value in enumerate(irreps[1])
    }
    answer = {}
    for element, quaternion in enumerate(binary["group"]):
        a, b = quaternion[0]
        trace = sy.Rational(2*a.numerator, a.denominator)
        trace += sy.Rational(2*b.numerator, b.denominator) * phi
        answer[element] = trace_to_class[sy.simplify(trace)]

    mul = binary["mul"]
    inverse = binary["inverse"]
    for element in range(GROUP_ORDER):
        conjugates = {
            mul[mul[g][element]][inverse[g]] for g in range(GROUP_ORDER)
        }
        assert {answer[g] for g in conjugates} == {answer[element]}
    return answer


def induction_multiplicities(subgroup, order, harmonic, class_of, irreps):
    values = []
    for character in irreps:
        inner = sum(
            character[class_of[element]]
            * sy.exp(-2 * sy.pi * sy.I * harmonic * k / order)
            for k, element in enumerate(subgroup)
        ) / order
        inner = sy.simplify(sy.expand_complex(inner))
        assert inner.is_Integer
        values.append(int(inner))
    return tuple(values)


def finite_character_table(prime, root):
    inv = lambda value: pow(int(value), -1, prime)
    phi = (pow(root, 6, prime) + pow(root, -6, prime)) % prime
    x = (2, -2, 0, 1, -1, phi, -phi, phi-1, 1-phi)
    xp = (2, -2, 0, 1, -1, 1-phi, phi-1, -phi, phi)
    x = tuple(value % prime for value in x)
    xp = tuple(value % prime for value in xp)

    def symmetric_power(power, defining):
        if power == 0:
            return (1,) * 9
        if power == 1:
            return defining
        older = (1,) * 9
        old = defining
        for _ in range(2, power+1):
            new = tuple((v*a-b) % prime
                        for v, a, b in zip(defining, old, older))
            older, old = old, new
        return old

    table = (
        symmetric_power(0, x),
        symmetric_power(1, x),
        symmetric_power(1, xp),
        symmetric_power(2, x),
        symmetric_power(2, xp),
        tuple((a*b) % prime for a, b in zip(x, xp)),
        symmetric_power(3, x),
        symmetric_power(4, x),
        symmetric_power(5, x),
    )
    assert inv(120)
    return table


def matrix_rank_mod(matrix, prime):
    work = np.asarray(matrix, dtype=np.int64).copy() % prime
    rows, columns = work.shape
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows)
                      if work[row, column] % prime), None)
        if pivot is None:
            continue
        if pivot != rank:
            work[[rank, pivot]] = work[[pivot, rank]]
        work[rank] = (work[rank] * pow(int(work[rank, column]), -1, prime)) % prime
        for row in range(rows):
            if row != rank and work[row, column]:
                work[row] = (work[row] - work[row, column]*work[rank]) % prime
        rank += 1
        if rank == rows:
            break
    return rank


def representation_matrix(module, element, prime, root):
    dimension = GROUP_ORDER // module["order"]
    matrix = np.zeros((dimension, dimension), dtype=np.int64)
    step = ROOT_ORDER // module["order"]
    for column, (row, exponent) in enumerate(module["actions"][element]):
        matrix[row, column] = pow(
            root, module["harmonic"] * exponent * step, prime
        )
    return matrix


def central_projectors(module, prime, root, class_of, finite_irreps):
    dimension = GROUP_ORDER // module["order"]
    actions = [representation_matrix(module, g, prime, root)
               for g in range(GROUP_ORDER)]
    projectors = []
    for irrep_dimension, character in zip(IRREP_DIMS, finite_irreps):
        projector = np.zeros((dimension, dimension), dtype=np.int64)
        for element, action in enumerate(actions):
            projector += character[class_of[element]] * action
            projector %= prime
        projector *= irrep_dimension * pow(GROUP_ORDER, -1, prime)
        projector %= prime
        projectors.append(projector)
    return tuple(projectors)


def incidence_exponents(source, target, source_coset):
    """Unique orbit kernel supported on the fixed incident flag orbit."""
    rows = GROUP_ORDER // target["order"]
    columns = GROUP_ORDER // source["order"]
    exponents = np.full((rows, columns), -1, dtype=np.int16)
    multiplicities = np.zeros((rows, columns), dtype=np.int16)
    source_step = ROOT_ORDER // source["order"]
    target_step = ROOT_ORDER // target["order"]
    for group_element in range(GROUP_ORDER):
        row, target_exponent = target["actions"][group_element][0]
        column, source_exponent = source["actions"][group_element][source_coset]
        value = (
            target["harmonic"] * target_exponent * target_step
            - source["harmonic"] * source_exponent * source_step
        ) % ROOT_ORDER
        if exponents[row, column] >= 0:
            assert exponents[row, column] == value
        else:
            exponents[row, column] = value
        multiplicities[row, column] += 1
    nonzero = exponents >= 0
    # An incident flag is stabilized exactly by the binary center.
    assert np.all(multiplicities[nonzero] == 2)
    assert np.count_nonzero(nonzero) == 60
    return exponents


def double_cosets(mul, identity, left_subgroup, right_subgroup):
    unseen = set(range(GROUP_ORDER))
    records = []
    while unseen:
        representative = identity if identity in unseen else min(unseen)
        members = {
            mul[mul[left][representative]][right]
            for left in left_subgroup for right in right_subgroup
        }
        records.append({"representative": representative,
                        "members": frozenset(members)})
        unseen -= members
    return records


def icosahedral_cell_certificate(binary, bases, identity, center, x, y, z):
    """Find incidence cosets from the 12-vertex combinatorial icosahedron."""
    group = binary["group"]
    mul = binary["mul"]
    q_mul = binary["q_mul"]
    q_conj = binary["q_conj"]
    zp_add = binary["zp_add"]
    zp_sub = binary["zp_sub"]
    zp_mul = binary["zp_mul"]
    zero = binary["zero"]

    def vector_add(left, right):
        return tuple(zp_add(a, b) for a, b in zip(left, right))

    def vector_dot(left, right):
        answer = zero
        for a, b in zip(left, right):
            answer = zp_add(answer, zp_mul(a, b))
        return answer

    def vector_cross(left, right):
        return (
            zp_sub(zp_mul(left[1], right[2]), zp_mul(left[2], right[1])),
            zp_sub(zp_mul(left[2], right[0]), zp_mul(left[0], right[2])),
            zp_sub(zp_mul(left[0], right[1]), zp_mul(left[1], right[0])),
        )

    def rotate_exact(quaternion, vector):
        pure = (zero,) + tuple(vector)
        return q_mul(q_mul(quaternion, pure), q_conj(quaternion))[1:]

    def zp_expression(value):
        a, b = value
        return sy.Rational(a.numerator, a.denominator) + sy.Rational(
            b.numerator, b.denominator
        ) * (1 + sy.sqrt(5))/2

    def positive(value):
        expression = sy.simplify(zp_expression(value))
        assert expression.is_positive is not None
        return bool(expression.is_positive)

    z_axis = group[z][1:]
    vertex_keys = sorted({rotate_exact(quaternion, z_axis)
                          for quaternion in group})
    assert len(vertex_keys) == 12
    vertices = tuple(vertex_keys)
    vertex_index = {vertex: number for number, vertex in enumerate(vertices)}
    vertex_actions = []
    for quaternion in group:
        permutation = [vertex_index[rotate_exact(quaternion, vertex)]
                       for vertex in vertices]
        assert len(set(permutation)) == 12
        vertex_actions.append(tuple(permutation))

    distinct_dot_values = {
        vector_dot(vertices[i], vertices[j])
        for i in range(12) for j in range(i+1, 12)
    }
    distinct_dots = sorted(distinct_dot_values,
                           key=lambda value: float(zp_expression(value)))
    edge_dot = distinct_dots[-1]
    edges = tuple((i, j) for i in range(12) for j in range(i+1, 12)
                  if vector_dot(vertices[i], vertices[j]) == edge_dot)
    assert len(edges) == 30
    edge_set = set(edges)
    faces = tuple((i, j, k) for i in range(12) for j in range(i+1, 12)
                  for k in range(j+1, 12)
                  if (i, j) in edge_set and (i, k) in edge_set
                  and (j, k) in edge_set)
    assert len(faces) == 20

    edge_index = {edge: number for number, edge in enumerate(edges)}
    face_index = {face: number for number, face in enumerate(faces)}

    def move_cell(cell_type, cell, element):
        permutation = vertex_actions[element]
        if cell_type == "V":
            return permutation[cell]
        moved = tuple(sorted(permutation[vertex] for vertex in cell))
        if cell_type == "E":
            assert moved in edge_index
        else:
            assert moved in face_index
        return moved

    # Fix the base cells without reading the incidence relation: choose the
    # exact cell-centre ray parallel to the oriented imaginary part of each
    # generator.  Parallelism and orientation are both checked in Q(phi).
    axes = {"V": group[z][1:], "E": group[x][1:], "F": group[y][1:]}
    cell_lists = {"V": tuple(range(12)), "E": edges, "F": faces}

    def centre(cell_type, cell):
        if cell_type == "V":
            return vertices[cell]
        answer = (zero, zero, zero)
        for vertex in cell:
            answer = vector_add(answer, vertices[vertex])
        return answer

    base_cells = {}
    for cell_type, cells in cell_lists.items():
        aligned = [cell for cell in cells
                   if vector_cross(centre(cell_type, cell), axes[cell_type])
                   == (zero, zero, zero)
                   and positive(vector_dot(centre(cell_type, cell), axes[cell_type]))]
        assert len(aligned) == 1
        base_cells[cell_type] = aligned[0]

    # Stabilizers are checked from the actual permutation action, not inferred
    # from group orders.
    for cell_type, base in base_cells.items():
        stabilizer = {
            element for element in range(GROUP_ORDER)
            if move_cell(cell_type, base, element) == base
        }
        assert stabilizer == set(bases[cell_type]["subgroup"])

    def incident(target_type, target_cell, source_type, source_cell):
        pair = {target_type, source_type}
        if pair == {"V", "E"}:
            vertex = target_cell if target_type == "V" else source_cell
            edge = target_cell if target_type == "E" else source_cell
            return vertex in edge
        if pair == {"V", "F"}:
            vertex = target_cell if target_type == "V" else source_cell
            face = target_cell if target_type == "F" else source_cell
            return vertex in face
        if pair == {"E", "F"}:
            edge = target_cell if target_type == "E" else source_cell
            face = target_cell if target_type == "F" else source_cell
            return set(edge) <= set(face)
        raise AssertionError(pair)

    relation_records = {}
    pair_types = (("E", "V"), ("F", "V"), ("F", "E"))
    expected = {("E", "V"): (6, 20),
                ("F", "V"): (4, 30),
                ("F", "E"): (10, 12)}
    for source_type, target_type in pair_types:
        double_coset_records = double_cosets(
            mul, identity, bases[target_type]["subgroup"],
            bases[source_type]["subgroup"]
        )
        expected_count, expected_size = expected[(source_type, target_type)]
        assert len(double_coset_records) == expected_count
        assert all(len(record["members"]) == expected_size
                   for record in double_coset_records)
        intersection_orders = []
        for record in double_coset_records:
            representative = record["representative"]
            conjugate_source = {
                mul[mul[representative][element]][binary["inverse"][representative]]
                for element in bases[source_type]["subgroup"]
            }
            intersection = set(bases[target_type]["subgroup"]) & conjugate_source
            assert intersection == {identity, center}
            intersection_orders.append(len(intersection))
        incidence_cosets = []
        orbit_incidence_counts = []
        for number, record in enumerate(double_coset_records):
            moved_source = move_cell(
                source_type, base_cells[source_type], record["representative"]
            )
            if incident(target_type, base_cells[target_type],
                        source_type, moved_source):
                incidence_cosets.append(number)

            pair_orbit = set()
            for group_element in range(GROUP_ORDER):
                target_cell = move_cell(
                    target_type, base_cells[target_type], group_element
                )
                source_cell = move_cell(source_type, moved_source, group_element)
                pair_orbit.add((target_cell, source_cell))
            assert len(pair_orbit) == 60
            orbit_incidence_counts.append(sum(
                incident(target_type, target_cell, source_type, source_cell)
                for target_cell, source_cell in pair_orbit
            ))
        assert len(incidence_cosets) == 1
        incidence_number = incidence_cosets[0]
        assert orbit_incidence_counts[incidence_number] == 60
        assert all(count == 0 for number, count in enumerate(orbit_incidence_counts)
                   if number != incidence_number)
        representative = double_coset_records[incidence_number]["representative"]
        source_coset = bases[source_type]["element_to_coset"][representative]
        relation_records[(source_type, target_type)] = {
            "source_type": source_type,
            "target_type": target_type,
            "double_coset_count": len(double_coset_records),
            "double_coset_sizes": [len(record["members"])
                                    for record in double_coset_records],
            "intersection_orders": intersection_orders,
            "incidence_double_coset_index": incidence_number,
            "incidence_double_coset_representative": representative,
            "incidence_source_coset": source_coset,
            "pair_orbit_incidence_counts": orbit_incidence_counts,
        }

    return {
        "construction": (
            "quaternion orbit of the oriented C10 axis; graph edges are the "
            "maximal distinct vertex-dot-product pairs; faces are 3-cliques"
        ),
        "normalized_vertex_dot_products_exact": [
            str(sy.simplify(zp_expression(value)
                            / zp_expression(vector_dot(vertices[0], vertices[0]))))
            for value in distinct_dots
        ],
        "cell_counts": {"vertices": len(vertices), "edges": len(edges),
                        "faces": len(faces)},
        "base_cells": {"V": base_cells["V"],
                       "E": list(base_cells["E"]),
                       "F": list(base_cells["F"])},
        "base_cells_form_flag": (
            base_cells["V"] in base_cells["E"]
            and set(base_cells["E"]) <= set(base_cells["F"])
        ),
        "relations": relation_records,
    }


def finite_matrix(exponents, prime, root):
    matrix = np.zeros(exponents.shape, dtype=np.int64)
    for exponent in range(ROOT_ORDER):
        matrix[exponents == exponent] = pow(root, exponent, prime)
    return matrix


def exact_channel_ranks(exponents, source, target, modular_data,
                        group_generators):
    upper_bounds = tuple(min(a, b) for a, b in zip(
        source["irrep_multiplicities"], target["irrep_multiplicities"]
    ))
    prime_results = []
    for prime, root, projectors in modular_data:
        matrix = finite_matrix(exponents, prime, root)
        for generator in group_generators:
            source_action = representation_matrix(source, generator, prime, root)
            target_action = representation_matrix(target, generator, prime, root)
            assert np.array_equal(
                (target_action @ matrix) % prime,
                (matrix @ source_action) % prime,
            )
        ranks = []
        for irrep_dimension, projector in zip(IRREP_DIMS, projectors):
            channel_rank = matrix_rank_mod((matrix @ projector) % prime, prime)
            assert channel_rank % irrep_dimension == 0
            ranks.append(channel_rank // irrep_dimension)
        prime_results.append(tuple(ranks))
    # Each modular rank is a lower bound on characteristic-zero rank.  The
    # common multiplicity is an upper bound.  Equality is an exact sandwich.
    assert all(result == upper_bounds for result in prime_results)
    return upper_bounds, tuple(prime_results)


def reduced_root_vectors():
    variable = sy.symbols("z")
    cyclotomic = sy.Poly(sy.cyclotomic_poly(ROOT_ORDER, variable), variable,
                         domain=sy.ZZ)
    degree = cyclotomic.degree()
    vectors = []
    for exponent in range(ROOT_ORDER):
        remainder = sy.Poly(variable**exponent, variable, domain=sy.ZZ).rem(cyclotomic)
        vectors.append(tuple(int(remainder.nth(k)) for k in range(degree)))
    return np.asarray(vectors, dtype=np.int64)


def exact_composition(first, second, root_vectors):
    # first: E x F; second: V x E.  Every entry is a sum of 60th roots.
    rows, middle = second.shape
    assert first.shape[0] == middle
    columns = first.shape[1]
    reduced = np.zeros((rows, columns, root_vectors.shape[1]), dtype=np.int64)
    for row in range(rows):
        for column in range(columns):
            for k in range(middle):
                if second[row, k] >= 0 and first[k, column] >= 0:
                    exponent = (int(second[row, k]) + int(first[k, column])) % ROOT_ORDER
                    reduced[row, column] += root_vectors[exponent]
    nonzero_entries = int(np.count_nonzero(np.any(reduced != 0, axis=2)))
    digest = hashlib.sha256(reduced.tobytes()).hexdigest()
    return nonzero_entries == 0, nonzero_entries, digest


def character_counter(vectors):
    counter = Counter(tuple(vector) for vector in vectors)
    return [
        {"multiplicity": multiplicity, "irrep_multiplicities": list(vector)}
        for vector, multiplicity in sorted(counter.items())
    ]


def enumerate_preregistration():
    binary = quiet_run(HERE / "verify_nonnormal_c10_selection.py")
    harmonic = quiet_run(HERE / "verify_fibonacci_nonbinary_dynamics.py")
    mul = binary["mul"]
    inverse = binary["inverse"]
    identity = binary["identity"]
    center = binary["central_minus"]
    element_order = binary["element_order"]
    irreps = harmonic["irreps_2i"]
    class_of = conjugacy_class_map(binary, irreps)

    triangles = []
    for x in range(GROUP_ORDER):
        if element_order(x) != 4:
            continue
        for y in range(GROUP_ORDER):
            if element_order(y) != 6:
                continue
            z = mul[inverse[mul[x][y]]][center]
            if element_order(z) == 10:
                triangles.append((x, y, z))
    x, y, z = min(triangles)
    assert len(triangles) == 120
    assert group_power(mul, identity, x, 2) == center
    assert group_power(mul, identity, y, 3) == center
    assert group_power(mul, identity, z, 5) == center
    assert mul[mul[x][y]][z] == center

    specifications = (("V", 10, z), ("E", 4, x), ("F", 6, y))
    bases = {}
    modules = {}
    module_records = []
    for cell_type, order, generator in specifications:
        subgroup, subgroup_exponent = cyclic_subgroup(
            mul, identity, generator, order
        )
        representatives, element_to_coset, actions = left_coset_action(
            mul, inverse, identity, subgroup, subgroup_exponent
        )
        bases[cell_type] = {
            "cell_type": cell_type,
            "order": order,
            "generator": generator,
            "subgroup": subgroup,
            "subgroup_exponent": subgroup_exponent,
            "representatives": representatives,
            "element_to_coset": element_to_coset,
            "actions": actions,
        }
        for harmonic_number in range(order):
            module_id = f"{cell_type}{harmonic_number}"
            multiplicities = induction_multiplicities(
                subgroup, order, harmonic_number, class_of, irreps
            )
            module = dict(
                bases[cell_type], id=module_id, harmonic=harmonic_number,
                irrep_multiplicities=multiplicities
            )
            modules[module_id] = module
            module_records.append({
                "id": module_id,
                "cell_type": cell_type,
                "stabilizer": f"C{order}",
                "harmonic": harmonic_number,
                "dimension": GROUP_ORDER // order,
                "central_parity": "even" if harmonic_number % 2 == 0 else "odd",
                "irrep_multiplicities": list(multiplicities),
            })

    geometry = icosahedral_cell_certificate(
        binary, bases, identity, center, x, y, z
    )

    induction_matrix = sy.Matrix([
        record["irrep_multiplicities"] for record in module_records
    ])
    gram = induction_matrix * induction_matrix.T
    diagonal_histogram = Counter(int(gram[i, i]) for i in range(20))
    off_diagonal_histogram = Counter(
        int(gram[i, j]) for i in range(20) for j in range(20) if i != j
    )
    assert induction_matrix.rank() == 9
    assert diagonal_histogram == Counter({3: 8, 4: 2, 7: 4,
                                          8: 2, 15: 2, 16: 2})
    assert off_diagonal_histogram == Counter({
        0: 200, 2: 32, 3: 8, 4: 60, 6: 48,
        7: 4, 10: 24, 14: 2, 15: 2,
    })
    assert off_diagonal_histogram[1] == 0

    modular_module_data = {}
    modular_metadata = []
    for prime in GOOD_PRIMES:
        assert sy.isprime(prime) and (prime-1) % ROOT_ORDER == 0
        primitive = int(sy.primitive_root(prime))
        root = pow(primitive, (prime-1)//ROOT_ORDER, prime)
        assert pow(root, ROOT_ORDER, prime) == 1
        assert all(pow(root, ROOT_ORDER//q, prime) != 1 for q in (2, 3, 5))
        finite_irreps = finite_character_table(prime, root)
        modular_metadata.append({"prime": prime, "primitive_60th_root": root})
        for module_id, module in modules.items():
            projectors = central_projectors(
                module, prime, root, class_of, finite_irreps
            )
            projector_multiplicities = tuple(
                matrix_rank_mod(projector, prime) // dimension
                for projector, dimension in zip(projectors, IRREP_DIMS)
            )
            assert projector_multiplicities == module["irrep_multiplicities"]
            modular_module_data.setdefault(module_id, []).append(
                (prime, root, projectors)
            )

    # Boundary orientation: F -> E -> V.  V-F is the direct flag relation.
    pair_types = (("E", "V", "edge has endpoint vertex"),
                  ("F", "E", "face has boundary edge"),
                  ("F", "V", "face has vertex"))
    operator_records = []
    exponent_matrices = {}
    operator_lookup = {}
    for source_type, target_type, relation in pair_types:
        geometric_relation = geometry["relations"][(source_type, target_type)]
        incidence_source_coset = geometric_relation["incidence_source_coset"]
        source_order = bases[source_type]["order"]
        target_order = bases[target_type]["order"]
        for source_harmonic in range(source_order):
            for target_harmonic in range(target_order):
                # The incident-flag stabilizer is the center C2.  Its two
                # restrictions agree exactly when the harmonic parities do.
                if source_harmonic % 2 != target_harmonic % 2:
                    continue
                source = modules[f"{source_type}{source_harmonic}"]
                target = modules[f"{target_type}{target_harmonic}"]
                exponents = incidence_exponents(
                    source, target, incidence_source_coset
                )
                channel_ranks, ranks_by_prime = exact_channel_ranks(
                    exponents, source, target,
                    modular_module_data[source["id"]], (x, y)
                )
                kernel = tuple(a-r for a, r in zip(
                    source["irrep_multiplicities"], channel_ranks
                ))
                cokernel = tuple(a-r for a, r in zip(
                    target["irrep_multiplicities"], channel_ranks
                ))
                index = tuple(a-b for a, b in zip(
                    source["irrep_multiplicities"],
                    target["irrep_multiplicities"]
                ))
                operator_id = f"I{len(operator_records):03d}"
                full_hom_dimension = int(
                    gram[[record["id"] for record in module_records].index(source["id"]),
                         [record["id"] for record in module_records].index(target["id"])]
                )
                record = {
                    "id": operator_id,
                    "source": source["id"],
                    "target": target["id"],
                    "adjoint_direction": f"{target['id']}->{source['id']}",
                    "full_equivariant_hom_dimension": full_hom_dimension,
                    "geometric_selector": relation,
                    "incidence_double_coset_index": geometric_relation[
                        "incidence_double_coset_index"
                    ],
                    "incidence_double_coset_representative": geometric_relation[
                        "incidence_double_coset_representative"
                    ],
                    "double_coset_basis_support": [
                        1 if number == geometric_relation["incidence_double_coset_index"]
                        else 0
                        for number in range(geometric_relation["double_coset_count"])
                    ],
                    "incident_flag_pairs": 60,
                    "incident_flag_stabilizer": "C2={1,-1}",
                    "incidence_supported_hom_dimension": 1,
                    "alternative_reading_changes_hom_point": False,
                    "canonicity_reason": (
                        "support on exactly this flag-incidence orbit plus "
                        "2I-equivariance leaves one line; changing the base "
                        "fiber trivialization rescales the whole operator only"
                    ),
                    "normalization": "coefficient +1 at the fixed incident flag",
                    "matrix_shape": [GROUP_ORDER//target["order"],
                                     GROUP_ORDER//source["order"]],
                    "matrix_nonzero_entries": 60,
                    "channel_ranks": list(channel_ranks),
                    "channel_upper_bounds": list(channel_ranks),
                    "modular_channel_rank_certificates": [list(x) for x in ranks_by_prime],
                    "matrix_rank": sum(r*d for r, d in zip(channel_ranks, IRREP_DIMS)),
                    "kernel_dimension": sum(r*d for r, d in zip(kernel, IRREP_DIMS)),
                    "cokernel_dimension": sum(r*d for r, d in zip(cokernel, IRREP_DIMS)),
                    "kernel_irrep_multiplicities": list(kernel),
                    "cokernel_irrep_multiplicities": list(cokernel),
                    "index_irrep_multiplicities": list(index),
                    "adjoint_index_irrep_multiplicities": [-x for x in index],
                }
                operator_records.append(record)
                exponent_matrices[operator_id] = exponents
                operator_lookup[(source["id"], target["id"])] = record

    assert len(operator_records) == 20 + 12 + 30 == 62

    root_vectors = reduced_root_vectors()
    complex_candidates = []
    short_complexes = []
    for face_harmonic in range(6):
        for edge_harmonic in range(4):
            for vertex_harmonic in range(10):
                first = operator_lookup.get((f"F{face_harmonic}",
                                             f"E{edge_harmonic}"))
                second = operator_lookup.get((f"E{edge_harmonic}",
                                              f"V{vertex_harmonic}"))
                if first is None or second is None:
                    continue
                is_zero, nonzero_entries, digest = exact_composition(
                    exponent_matrices[first["id"]],
                    exponent_matrices[second["id"]],
                    root_vectors,
                )
                index = tuple(
                    modules[f"F{face_harmonic}"]["irrep_multiplicities"][i]
                    - modules[f"E{edge_harmonic}"]["irrep_multiplicities"][i]
                    + modules[f"V{vertex_harmonic}"]["irrep_multiplicities"][i]
                    for i in range(9)
                )
                candidate = {
                    "modules": [f"F{face_harmonic}", f"E{edge_harmonic}",
                                f"V{vertex_harmonic}"],
                    "maps": [first["id"], second["id"]],
                    "composition_zero_exactly": is_zero,
                    "composition_nonzero_matrix_entries": nonzero_entries,
                    "composition_reduced_matrix_sha256": digest,
                    "euler_index_irrep_multiplicities": list(index),
                }
                complex_candidates.append(candidate)
                if is_zero:
                    rank_first = first["channel_ranks"]
                    rank_second = second["channel_ranks"]
                    h2 = first["kernel_irrep_multiplicities"]
                    h0 = second["cokernel_irrep_multiplicities"]
                    edge_mult = modules[f"E{edge_harmonic}"]["irrep_multiplicities"]
                    h1 = tuple(edge_mult[i] - rank_first[i] - rank_second[i]
                               for i in range(9))
                    assert all(value >= 0 for value in h1)
                    assert tuple(h2[i]-h1[i]+h0[i] for i in range(9)) == index
                    short_complexes.append({
                        "id": f"K{len(short_complexes):03d}",
                        **candidate,
                        "H2_irrep_multiplicities": list(h2),
                        "H1_irrep_multiplicities": list(h1),
                        "H0_irrep_multiplicities": list(h0),
                        "H2_dimension": sum(a*b for a, b in zip(h2, IRREP_DIMS)),
                        "H1_dimension": sum(a*b for a, b in zip(h1, IRREP_DIMS)),
                        "H0_dimension": sum(a*b for a, b in zip(h0, IRREP_DIMS)),
                    })
    assert len(complex_candidates) == 60

    object_indices = [record["index_irrep_multiplicities"]
                      for record in operator_records]
    object_indices += [record["euler_index_irrep_multiplicities"]
                       for record in short_complexes]
    directed_indices = object_indices + [tuple(-x for x in index)
                                         for index in object_indices]
    kernel_characters = [record["kernel_irrep_multiplicities"]
                         for record in operator_records]
    cokernel_characters = [record["cokernel_irrep_multiplicities"]
                           for record in operator_records]
    middle_characters = [record["H1_irrep_multiplicities"]
                         for record in short_complexes]

    return {
        "schema": "binary-235-incidence-preregistration-v1",
        "status": "TARGET_INDEPENDENT_NO_EXTERNAL_MODULE_COMPARISON",
        "scope_and_hypotheses": {
            "geometry": (
                "the 2I action on the vertex-edge-face flags of the "
                "icosahedral boundary, with stabilizers C10,C4,C6"
            ),
            "line_modules": "every complex character of each cyclic stabilizer",
            "operator": (
                "the unique-up-to-scale equivariant kernel supported exactly "
                "on one of the three actual flag-incidence relations"
            ),
            "adjointness": "the reverse incidence operator is the Hermitian adjoint",
            "normalization": "one fixed incident-flag coefficient is +1",
            "short_complex": (
                "a normalized F-to-E incidence followed by normalized E-to-V "
                "incidence whose composition vanishes exactly"
            ),
            "not_counted_as_canonical": [
                "an arbitrary linear combination of distinct double-coset kernels",
                "an arbitrary choice of Schur-channel coefficients",
                "a nonincident double-coset relation without additional cell data",
                "a direct sum selected after an external comparison",
            ],
            "counting_unit": (
                "one incidence map together with its forced adjoint; one short "
                "complex together with its adjoint cochain complex"
            ),
        },
        "group_certificate": {
            "group": "2I",
            "order": GROUP_ORDER,
            "triangle_generators_by_exact_group_index": {"x_C4": x,
                                                          "y_C6": y,
                                                          "z_C10": z},
            "relations": "x^2=y^3=z^5=xyz=-1",
            "compatible_triangle_flags": len(triangles),
        },
        "icosahedral_cell_certificate": {
            key: value for key, value in geometry.items() if key != "relations"
        } | {
            "relations": [geometry["relations"][key]
                          for key in (("E", "V"), ("F", "V"), ("F", "E"))]
        },
        "irrep_order": list(IRREP_NAMES),
        "irrep_dimensions": list(IRREP_DIMS),
        "modular_exact_rank_certificate": modular_metadata,
        "induction_matrix_rank": induction_matrix.rank(),
        "induction_relation_lattice_rank": 20 - induction_matrix.rank(),
        "modules": module_records,
        "hom_gram_matrix": [list(map(int, gram.row(i))) for i in range(20)],
        "hom_diagonal_histogram": {str(k): v for k, v in sorted(diagonal_histogram.items())},
        "hom_off_diagonal_ordered_histogram": {
            str(k): v for k, v in sorted(off_diagonal_histogram.items())
        },
        "off_diagonal_pairs_with_hom_dimension_one": off_diagonal_histogram[1],
        "operators": operator_records,
        "complex_candidates": complex_candidates,
        "short_complexes": short_complexes,
        "preregistered_index_multiset_up_to_adjoint": character_counter(object_indices),
        "preregistered_index_multiset_with_adjoints": character_counter(directed_indices),
        "operator_kernel_character_multiset": character_counter(kernel_characters),
        "operator_cokernel_character_multiset": character_counter(cokernel_characters),
        "short_complex_middle_cohomology_character_multiset": character_counter(
            middle_characters
        ),
        "counts": {
            "modules": len(module_records),
            "canonical_incidence_operator_adjoint_pairs": len(operator_records),
            "canonical_incidence_operators_with_adjoints": 2*len(operator_records),
            "composable_incidence_pairs_tested": len(complex_candidates),
            "short_complex_adjoint_pairs": len(short_complexes),
            "short_complexes_with_adjoints": 2*len(short_complexes),
            "N_canonical_objects_up_to_adjoint": len(object_indices),
            "N_canonical_objects_with_adjoints": len(directed_indices),
            "distinct_index_characters_up_to_adjoint": len({tuple(x) for x in object_indices}),
            "distinct_index_characters_with_adjoints": len({tuple(x) for x in directed_indices}),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    data = enumerate_preregistration()
    payload = json.dumps(data, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(payload, end="")
    else:
        arguments.output.write_text(payload, encoding="utf-8")
    counts = data["counts"]
    print(
        f"operators={counts['canonical_incidence_operator_adjoint_pairs']} "
        f"complexes={counts['short_complex_adjoint_pairs']} "
        f"N={counts['N_canonical_objects_up_to_adjoint']} "
        f"distinct={counts['distinct_index_characters_up_to_adjoint']}",
        file=__import__("sys").stderr,
    )


if __name__ == "__main__":
    main()

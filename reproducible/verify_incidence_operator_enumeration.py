#!/usr/bin/env python3
"""Target-independent verifier for the corrected (2,3,5) incidence census.

This supersedes the three-operator calculation in commit 36bd682.  That
calculation averaged a twisted kernel over both full stabilizers, so the sum
factorized and killed every nontrivial character.  A Mackey orbit kernel is
instead a covariant section over (H x K)/C2.  It is nonzero precisely when the
two characters have the same central parity.

No external comparison module is defined or loaded here.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

from preregister_orbifold_incidence import enumerate_preregistration, quiet_run


HERE = Path(__file__).resolve().parent
DATA_PATH = HERE / "orbifold_incidence_preregistered.json"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def literal_right_coset_crosscheck(c6_axis_sign=1):
    """Reproduce the old coset labels with exact combinatorial cells.

    This uses the same first order-10/4/6 generators and right-coset ordering
    as commit 36bd682, but it never diagonalizes a rotation matrix.  Vertices
    are an exact quaternion orbit, edges are maximal-dot-product pairs, and
    faces are graph triangles.
    """
    binary = quiet_run(HERE / "verify_nonnormal_c10_selection.py")
    group = binary["group"]
    mul = binary["mul"]
    inverse = binary["inverse"]
    identity = binary["identity"]
    order = binary["element_order"]
    q_mul = binary["q_mul"]
    q_conj = binary["q_conj"]
    zp_add = binary["zp_add"]
    zp_sub = binary["zp_sub"]
    zp_mul = binary["zp_mul"]
    zero = binary["zero"]

    generators = {
        "C10": next(g for g in range(120) if order(g) == 10),
        "C4": next(g for g in range(120) if order(g) == 4),
        "C6": next(g for g in range(120) if order(g) == 6),
    }

    def subgroup(generator):
        answer = []
        current = identity
        while current not in answer:
            answer.append(current)
            current = mul[current][generator]
        return tuple(answer)

    subgroups = {name: subgroup(generator)
                 for name, generator in generators.items()}

    def right_representatives(subgroup_elements):
        representatives = []
        seen = set()
        for element in range(120):
            if element in seen:
                continue
            representatives.append(element)
            seen |= {mul[h][element] for h in subgroup_elements}
        return tuple(representatives)

    representatives = {name: right_representatives(elements)
                       for name, elements in subgroups.items()}

    def rotate(quaternion, vector):
        pure = (zero,) + tuple(vector)
        return q_mul(q_mul(quaternion, pure), q_conj(quaternion))[1:]

    def add(left, right):
        return tuple(zp_add(a, b) for a, b in zip(left, right))

    def dot(left, right):
        answer = zero
        for a, b in zip(left, right):
            answer = zp_add(answer, zp_mul(a, b))
        return answer

    def cross(left, right):
        return (
            zp_sub(zp_mul(left[1], right[2]), zp_mul(left[2], right[1])),
            zp_sub(zp_mul(left[2], right[0]), zp_mul(left[0], right[2])),
            zp_sub(zp_mul(left[0], right[1]), zp_mul(left[1], right[0])),
        )

    def positive(value):
        a, b = value
        return float(a) + float(b)*(1 + 5**0.5)/2 > 0

    points = {}
    for name, generator in generators.items():
        axis = group[generator][1:]
        if name == "C6" and c6_axis_sign == -1:
            axis = tuple(binary["zp_neg"](coordinate) for coordinate in axis)
        points[name] = tuple(
            rotate(group[inverse[representative]], axis)
            for representative in representatives[name]
        )
        assert len(set(points[name])) == 120 // len(subgroups[name])

    vertices = points["C10"]
    dot_values = {dot(vertices[i], vertices[j])
                  for i in range(12) for j in range(i+1, 12)}
    edge_dot = max(dot_values,
                   key=lambda value: float(value[0])
                   + float(value[1])*(1 + 5**0.5)/2)
    edges = tuple((i, j) for i in range(12) for j in range(i+1, 12)
                  if dot(vertices[i], vertices[j]) == edge_dot)
    edge_set = set(edges)
    faces = tuple((i, j, k) for i in range(12) for j in range(i+1, 12)
                  for k in range(j+1, 12)
                  if (i, j) in edge_set and (i, k) in edge_set
                  and (j, k) in edge_set)
    assert (len(vertices), len(edges), len(faces)) == (12, 30, 20)

    def cell_centres(cells):
        answer = []
        for cell in cells:
            total = (zero, zero, zero)
            for vertex in cell:
                total = add(total, vertices[vertex])
            answer.append(total)
        return tuple(answer)

    edge_centres = cell_centres(edges)
    face_centres = cell_centres(faces)

    def identify_cells(point_vectors, cells, centres):
        identified = []
        for point in point_vectors:
            hits = [cell for cell, centre in zip(cells, centres)
                    if cross(point, centre) == (zero, zero, zero)
                    and positive(dot(point, centre))]
            assert len(hits) == 1
            identified.append(hits[0])
        return tuple(identified)

    cells = {
        "C10": tuple(range(12)),
        "C4": identify_cells(points["C4"], edges, edge_centres),
        "C6": identify_cells(points["C6"], faces, face_centres),
    }

    def double_coset_table(left, right):
        table = {}
        seen = set()
        number = 0
        for representative in range(120):
            if representative in seen:
                continue
            members = {mul[mul[h][representative]][k]
                       for h in left for k in right}
            for member in members:
                table[member] = number
            seen |= members
            number += 1
        return table, number

    def is_incident(left_type, left_cell, right_type, right_cell):
        if {left_type, right_type} == {"C10", "C4"}:
            vertex = left_cell if left_type == "C10" else right_cell
            edge = left_cell if left_type == "C4" else right_cell
            return vertex in edge
        if {left_type, right_type} == {"C10", "C6"}:
            vertex = left_cell if left_type == "C10" else right_cell
            face = left_cell if left_type == "C6" else right_cell
            return vertex in face
        edge = left_cell if left_type == "C4" else right_cell
        face = left_cell if left_type == "C6" else right_cell
        return set(edge) <= set(face)

    output = {}
    for left_type, right_type in (("C10", "C4"),
                                  ("C10", "C6"),
                                  ("C4", "C6")):
        table, number = double_coset_table(
            subgroups[left_type], subgroups[right_type]
        )
        tally = {index: [0, 0] for index in range(number)}
        for i, left_rep in enumerate(representatives[left_type]):
            for j, right_rep in enumerate(representatives[right_type]):
                # Right-coset convention: invariant of Hx,Ky is H*x*y^-1*K.
                index = table[mul[left_rep][inverse[right_rep]]]
                tally[index][0] += 1
                tally[index][1] += int(is_incident(
                    left_type, cells[left_type][i],
                    right_type, cells[right_type][j]
                ))
        pure = [index for index, (total, incident) in tally.items()
                if total == incident and total]
        assert len(pure) == 1
        assert sorted(incident for _, incident in tally.values())
        output[(left_type, right_type)] = (number, pure[0], tally)
    return output


print("=" * 78)
print("CORRECTED TARGET-INDEPENDENT (2,3,5) INCIDENCE ENUMERATION")
print("=" * 78)

committed = json.loads(DATA_PATH.read_text(encoding="utf-8"))
regenerated = enumerate_preregistration()
check("[DERIVED] committed census regenerates byte-for-object exactly",
      regenerated == committed)
check("[STRUCTURAL] census declares no external-module comparison",
      committed["status"] == "TARGET_INDEPENDENT_NO_EXTERNAL_MODULE_COMPARISON")

geometry = committed["icosahedral_cell_certificate"]
check("[DERIVED] exact Q(phi) cell construction has f-vector (12,30,20)",
      geometry["cell_counts"] == {"vertices": 12, "edges": 30, "faces": 20})
check("[DERIVED] exact vertex dot products are the icosahedral values",
      geometry["normalized_vertex_dot_products_exact"]
      == ["-1", "-sqrt(5)/5", "sqrt(5)/5"])

relations = geometry["relations"]
relation_signature = [
    (record["source_type"], record["target_type"],
     record["double_coset_count"], record["double_coset_sizes"])
    for record in relations
]
check("[DERIVED] cross-cell double-coset counts and sizes are exact",
      relation_signature == [
          ("E", "V", 6, [20]*6),
          ("F", "V", 4, [30]*4),
          ("F", "E", 10, [12]*10),
      ])
check("[DERIVED] each relation has one pure 60-pair incidence orbit",
      all(sorted(record["pair_orbit_incidence_counts"])
              == [0]*(record["double_coset_count"]-1) + [60]
              for record in relations))
check("[DERIVED] every cross-cell double-coset intersection is the center C2",
      all(record["intersection_orders"] == [2]*record["double_coset_count"]
          for record in relations))

literal_natural = literal_right_coset_crosscheck()
natural_indices = {
    pair: (count, incidence_index)
    for pair, (count, incidence_index, _) in literal_natural.items()
}
literal_eigen_sign = literal_right_coset_crosscheck(c6_axis_sign=-1)
eigen_sign_indices = {
    pair: (count, incidence_index)
    for pair, (count, incidence_index, _) in literal_eigen_sign.items()
}
check("[DERIVED] exact combinatorial cells expose the axis-sign label convention",
      natural_indices == {
          ("C10", "C4"): (6, 0),
          ("C10", "C6"): (4, 1),
          ("C4", "C6"): (10, 8),
      } and eigen_sign_indices == {
          ("C10", "C4"): (6, 0),
          ("C10", "C6"): (4, 3),
          ("C4", "C6"): (10, 6),
      }, f"oriented-axis labels={natural_indices}; C6-antipode labels={eigen_sign_indices}")

check("[DERIVED] all 20 induced rows span the rank-nine representation ring",
      committed["induction_matrix_rank"] == 9
      and committed["induction_relation_lattice_rank"] == 11)
check("[DERIVED] exact Hom diagonal histogram matches the Gram census",
      committed["hom_diagonal_histogram"]
      == {"3": 8, "4": 2, "7": 4, "8": 2, "15": 2, "16": 2})
check("[DERIVED] exact ordered off-diagonal Hom histogram has no dimension one",
      committed["hom_off_diagonal_ordered_histogram"] == {
          "0": 200, "2": 32, "3": 8, "4": 60, "6": 48,
          "7": 4, "10": 24, "14": 2, "15": 2,
      } and committed["off_diagonal_pairs_with_hom_dimension_one"] == 0)

operators = committed["operators"]
by_relation = Counter((record["source"][0], record["target"][0])
                      for record in operators)
check("[DERIVED] central parity leaves 20+30+12=62 incidence-map lines",
      by_relation == Counter({("E", "V"): 20,
                              ("F", "V"): 30,
                              ("F", "E"): 12})
      and len(operators) == 62)
check("[DERIVED] every map has support only on its incidence double coset",
      all(sum(record["double_coset_basis_support"]) == 1
          and record["double_coset_basis_support"][
              record["incidence_double_coset_index"]
          ] == 1
          and record["incidence_supported_hom_dimension"] == 1
          for record in operators))
check("[DERIVED] two good-prime rank certificates attain every exact upper bound",
      all(len(record["modular_channel_rank_certificates"]) == 2
          and all(certificate == record["channel_ranks"]
                  for certificate in record["modular_channel_rank_certificates"])
          and record["channel_ranks"] == record["channel_upper_bounds"]
          for record in operators))

dims = committed["irrep_dimensions"]
check("[DERIVED] all kernel and cokernel dimensions match their exact characters",
      all(sum(a*b for a, b in zip(record["kernel_irrep_multiplicities"], dims))
              == record["kernel_dimension"]
              and sum(a*b for a, b in zip(record["cokernel_irrep_multiplicities"], dims))
              == record["cokernel_dimension"]
              and record["matrix_rank"] + record["kernel_dimension"]
              == record["matrix_shape"][1]
              and record["matrix_rank"] + record["cokernel_dimension"]
              == record["matrix_shape"][0]
              for record in operators))

complexes = committed["short_complexes"]
check("[DERIVED] exactly one of 60 normalized boundary pairs is a complex",
      committed["counts"]["composable_incidence_pairs_tested"] == 60
      and len(complexes) == 1
      and complexes[0]["modules"] == ["F0", "E2", "V0"]
      and complexes[0]["composition_zero_exactly"])
check("[DERIVED] the cellular complex has H2=1, H1=0, H0=1",
      (complexes[0]["H2_irrep_multiplicities"],
       complexes[0]["H1_irrep_multiplicities"],
       complexes[0]["H0_irrep_multiplicities"])
      == ([1, 0, 0, 0, 0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0, 0, 0, 0, 0],
          [1, 0, 0, 0, 0, 0, 0, 0, 0]))

counts = committed["counts"]
check("[DERIVED] preregistered look-elsewhere census is N=63, 28 indices",
      counts["N_canonical_objects_up_to_adjoint"] == 63
      and counts["distinct_index_characters_up_to_adjoint"] == 28)

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT_OLD_THREE_OPERATOR_ENUMERATION=REFUTED")
print("VERDICT_CANONICAL_INCIDENCE_MAP_PAIRS=62")
print("VERDICT_CANONICAL_SHORT_COMPLEX_PAIRS=1")
print("NO_EXTERNAL_MODULE_COMPARISON_PERFORMED")
if passed != tests:
    raise SystemExit(1)

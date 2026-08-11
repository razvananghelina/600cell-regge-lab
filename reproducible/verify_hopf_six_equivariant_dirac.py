#!/usr/bin/env python3
"""Exact canonical equivariant-Dirac gate on the 936-state Hopf carrier.

The character calculation, eight-reading eligibility rule, normalized rook
operator and all acceptance boundaries were frozen in protocol commit
705edeb.  No Hessian or particle target is used.
"""

from collections import Counter, deque
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sparse
import sympy as sp


OUTPUT = Path(__file__).with_name("hopf_six_equivariant_dirac.json")
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
print("CANONICAL A5-EQUIVARIANT DIRAC GATE ON THE SIX-FIBRATION CARRIER")
print("="*78)

sqrt5 = sp.sqrt(5)
phi = (1+sqrt5)/2
phi_bar = (1-sqrt5)/2
beta = phi-1
alpha = -phi

# Exact class restrictions.  A5 order is (1A,2A,3A,5A,5B); the selected D5
# classes are (1, reflections, edge rotations, distance-two rotations).
a5_names = ("1", "3", "3p", "4", "5")
a5_dims = (1, 3, 3, 4, 5)
a5_characters = {
    "1": (1, 1, 1, 1, 1),
    "3": (3, -1, 0, phi, phi_bar),
    "3p": (3, -1, 0, phi_bar, phi),
    "4": (4, 0, 1, -1, -1),
    "5": (5, 1, -1, 0, 0),
}
d5_names = ("trivial", "reflection_sign",
            "positive_doublet", "negative_doublet")
d5_class_sizes = (1, 5, 2, 2)
d5_characters = {
    "trivial": (1, 1, 1, 1),
    "reflection_sign": (1, -1, 1, 1),
    "positive_doublet": (2, 0, beta, alpha),
    "negative_doublet": (2, 0, alpha, beta),
}
restricted_a5 = {
    name: (character[0], character[1], character[3], character[4])
    for name, character in a5_characters.items()
}


def exact_inner(left, right):
    return sp.simplify(sp.Rational(1, 10)*sum(
        size*a*b for size, a, b in zip(d5_class_sizes, left, right)
    ))


d5_orthonormal = all(
    exact_inner(d5_characters[left], d5_characters[right])
    == (1 if left == right else 0)
    for left in d5_names for right in d5_names
)
check("the four exact D5 characters are orthonormal", d5_orthonormal)

induction_multiplicities = {}
for d5_name in d5_names:
    induction_multiplicities[d5_name] = tuple(
        int(exact_inner(d5_characters[d5_name], restricted_a5[a5_name]))
        for a5_name in a5_names
    )
expected_inductions = {
    "trivial": (1, 0, 0, 0, 1),
    "reflection_sign": (0, 1, 1, 0, 0),
    "positive_doublet": (0, 1, 0, 1, 1),
    "negative_doublet": (0, 0, 1, 1, 1),
}
check("Frobenius reciprocity gives the four exact induced modules",
      induction_multiplicities == expected_inductions,
      "V0=1+5; V1=3+3'; V2=3+4+5; V3=3'+4+5")
check("the induced dimensions are exactly 6,6,12,12",
      tuple(sum(mult*dim for mult, dim in zip(
          induction_multiplicities[name], a5_dims
      )) for name in d5_names) == (6, 6, 12, 12))

induction_rows = tuple(induction_multiplicities[name] for name in d5_names)
hom_gram = sp.Matrix([
    [sum(a*b for a, b in zip(left, right))
     for right in induction_rows]
    for left in induction_rows
])
expected_hom_gram = sp.Matrix([
    [2, 0, 1, 1],
    [0, 2, 1, 1],
    [1, 1, 3, 2],
    [1, 1, 2, 3],
])
check("the complete induced-module Hom Gram matrix is exact",
      hom_gram == expected_hom_gram,
      f"Gram={hom_gram.tolist()}")

# Rebuild the eight carrier readings and their legal odd blocks.
nodes = (
    {"name": "trivial", "size": 6, "u_rank": 2, "v_rank": 2},
    {"name": "reflection_sign", "size": 6, "u_rank": 2, "v_rank": 0},
    {"name": "positive_doublet", "size": 12, "u_rank": 1, "v_rank": 1},
    {"name": "negative_doublet", "size": 12, "u_rank": 0, "v_rank": 1},
)
node_sizes = tuple(node["size"] for node in nodes)
unordered_pairs = tuple(combinations(range(4), 2))
reading_records = []
for priority, u_direction, v_direction in product(
        (("u", "v"), ("v", "u")), (1, -1), (1, -1)):
    direction = {"u": u_direction, "v": v_direction}

    def key(node_index):
        return tuple(direction[coordinate]*nodes[node_index][coordinate+"_rank"]
                     for coordinate in priority)

    order = tuple(sorted(range(4), key=key))
    position = {node: rank for rank, node in enumerate(order)}
    positive_cells = frozenset(
        (left, right) if position[left] < position[right]
        else (right, left)
        for left, right in unordered_pairs
    )
    legal_blocks = []
    links = set()
    for first in positive_cells:
        for second in positive_cells:
            i, j = first
            k, ell = second
            if i == ell or j == k:
                target = (ell, k)  # J transpose of the second positive cell
                if i == ell:
                    link = tuple(sorted((j, k)))
                    changed_side = "right"
                else:
                    link = tuple(sorted((i, ell)))
                    changed_side = "left"
                links.add(link)
                legal_blocks.append({
                    "source": first,
                    "target": target,
                    "link": link,
                    "changed_side": changed_side,
                })
    link_hom_dimensions = {
        link: int(hom_gram[link[0], link[1]]) for link in links
    }
    eligible = (len(links) == 3 and len(legal_blocks) == 8
                and all(value == 1 for value in link_hom_dimensions.values()))
    reading_records.append({
        "priority": priority,
        "u_direction": u_direction,
        "v_direction": v_direction,
        "order": order,
        "positive_cells": positive_cells,
        "legal_blocks": legal_blocks,
        "links": frozenset(links),
        "link_hom_dimensions": link_hom_dimensions,
        "eligible": eligible,
    })

eligible_readings = [record for record in reading_records if record["eligible"]]
ineligible_readings = [record for record in reading_records
                       if not record["eligible"]]
check("the canonical-equivariant eligibility hit fraction is exactly 4/8",
      len(eligible_readings) == len(ineligible_readings) == 4
      and all(record["priority"] == ("u", "v")
              for record in eligible_readings)
      and all(record["priority"] == ("v", "u")
              for record in ineligible_readings),
      "edge-first 4/4 eligible; reflection-first 0/4 eligible")
check("every ineligible reading is killed by the zero Hom V0<->V1",
      all((0, 1) in record["links"]
          and record["link_hom_dimensions"][(0, 1)] == 0
          for record in ineligible_readings))
check("the two-dimensional V2<->V3 Hom is never silently fitted",
      all((2, 3) not in record["links"] for record in reading_records),
      "no reading requires choosing a line in the dimension-two Hom")

# Fixed orthogonal direct-sum models for the four induced modules.
irrep_dimensions = dict(zip(a5_names, a5_dims))
node_irreps = (
    ("1", "5"),
    ("3", "3p"),
    ("3", "4", "5"),
    ("3p", "4", "5"),
)
node_slices = []
for irreps in node_irreps:
    offset = 0
    slices = {}
    for irrep in irreps:
        dimension = irrep_dimensions[irrep]
        slices[irrep] = (offset, offset+dimension)
        offset += dimension
    node_slices.append(slices)


def normalized_intertwiner(target, source):
    common = set(node_slices[target]) & set(node_slices[source])
    if len(common) != 1:
        raise ValueError("intertwiner is not unique")
    irrep = next(iter(common))
    target_start, target_end = node_slices[target][irrep]
    source_start, source_end = node_slices[source][irrep]
    matrix = np.zeros((node_sizes[target], node_sizes[source]), dtype=np.int8)
    dimension = target_end-target_start
    if dimension != source_end-source_start:
        raise RuntimeError("common irrep dimensions disagree")
    for offset in range(dimension):
        matrix[target_start+offset, source_start+offset] = 1
    return matrix, irrep


unique_intertwiners = {}
for left, right in unordered_pairs:
    if hom_gram[left, right] == 1:
        forward, irrep = normalized_intertwiner(right, left)
        backward, reverse_irrep = normalized_intertwiner(left, right)
        if reverse_irrep != irrep or not np.array_equal(backward, forward.T):
            raise RuntimeError("partial-isometry adjoint mismatch")
        unique_intertwiners[(right, left)] = forward
        unique_intertwiners[(left, right)] = backward
check("all four one-dimensional Hom lines have canonical partial isometries",
      len(unique_intertwiners) == 8
      and all(np.array_equal(matrix@matrix.T@matrix, matrix)
              for matrix in unique_intertwiners.values()))

# Full 936-state ordered-cell carrier.
cells = tuple((left, right) for left in range(4) for right in range(4)
              if left != right)
cell_offsets = {}
dimension = 0
for cell in cells:
    cell_offsets[cell] = dimension
    dimension += node_sizes[cell[0]]*node_sizes[cell[1]]
check("the explicit ordered-cell carrier has dimension 936", dimension == 936)


def basis_index(cell, left_coordinate, right_coordinate):
    return (cell_offsets[cell]
            + left_coordinate*node_sizes[cell[1]]+right_coordinate)


j_rows = []
j_cols = []
gamma_values_by_reading = []
for cell in cells:
    left, right = cell
    for a in range(node_sizes[left]):
        for b in range(node_sizes[right]):
            j_cols.append(basis_index(cell, a, b))
            j_rows.append(basis_index((right, left), b, a))
J = sparse.coo_matrix((np.ones(len(j_rows), dtype=np.int8),
                       (j_rows, j_cols)), shape=(dimension, dimension)).tocsr()
I_H = sparse.eye(dimension, dtype=np.int8, format="csr")
check("the explicit cell-transpose real structure squares to one",
      (J@J-I_H).nnz == 0)


def build_dirac(record, link_signs):
    rows = []
    cols = []
    values = []
    for block in record["legal_blocks"]:
        source = block["source"]
        target = block["target"]
        sign = link_signs[block["link"]]
        if block["changed_side"] == "left":
            # source=(i,j), target=(ell,j)
            source_left, shared_right = source
            target_left, target_right = target
            if target_right != shared_right:
                raise RuntimeError("left-changing block lost its shared side")
            T = unique_intertwiners[(target_left, source_left)]
            for target_a, source_a in zip(*np.nonzero(T)):
                for b in range(node_sizes[shared_right]):
                    rows.append(basis_index(target, target_a, b))
                    cols.append(basis_index(source, source_a, b))
                    values.append(sign)
        else:
            # source=(i,j), target=(i,k)
            shared_left, source_right = source
            target_left, target_right = target
            if target_left != shared_left:
                raise RuntimeError("right-changing block lost its shared side")
            T = unique_intertwiners[(target_right, source_right)]
            for target_b, source_b in zip(*np.nonzero(T)):
                for a in range(node_sizes[shared_left]):
                    rows.append(basis_index(target, a, target_b))
                    cols.append(basis_index(source, a, source_b))
                    values.append(sign)
    half = sparse.coo_matrix((np.asarray(values, dtype=np.int8), (rows, cols)),
                             shape=(dimension, dimension)).tocsr()
    return (half+half.T).tocsr(), half


def exact_first_order_tensor_support(record, link_signs, half):
    """Compare every nonzero odd block with its exact Kronecker formula."""
    expected_nnz = 0
    for block in record["legal_blocks"]:
        source = block["source"]
        target = block["target"]
        source_dimension = node_sizes[source[0]]*node_sizes[source[1]]
        target_dimension = node_sizes[target[0]]*node_sizes[target[1]]
        source_slice = slice(cell_offsets[source],
                             cell_offsets[source]+source_dimension)
        target_slice = slice(cell_offsets[target],
                             cell_offsets[target]+target_dimension)
        actual = half[target_slice, source_slice]
        sign = link_signs[block["link"]]
        if block["changed_side"] == "left":
            T = unique_intertwiners[(target[0], source[0])]
            expected = sparse.kron(
                sparse.csr_matrix(T),
                sparse.eye(node_sizes[source[1]], dtype=np.int8),
                format="csr",
            )
        else:
            T = unique_intertwiners[(target[1], source[1])]
            expected = sparse.kron(
                sparse.eye(node_sizes[source[0]], dtype=np.int8),
                sparse.csr_matrix(T),
                format="csr",
            )
        expected = sign*expected
        expected_nnz += expected.nnz
        if (actual-expected).nnz != 0:
            return False
    # Legal blocks have disjoint source/target cell rectangles, so this also
    # excludes any unlicensed nonzero position.
    return half.nnz == expected_nnz


def grading(record):
    positive = record["positive_cells"]
    diagonal = np.empty(dimension, dtype=np.int8)
    for cell in cells:
        start = cell_offsets[cell]
        stop = start+node_sizes[cell[0]]*node_sizes[cell[1]]
        diagonal[start:stop] = 1 if cell in positive else -1
    return sparse.diags(diagonal, dtype=np.int8, format="csr"), diagonal


def central_link_constraint_rank(links):
    variable_offsets = []
    offset = 0
    for size in node_sizes:
        variable_offsets.append(offset)
        offset += size*size
    rows = []

    def variable(node, row, col):
        return variable_offsets[node]+row*node_sizes[node]+col

    for source, target in sorted(links):
        T = unique_intertwiners[(target, source)].astype(int)
        # A_target*T-T*A_source=0.
        for row in range(node_sizes[target]):
            for col in range(node_sizes[source]):
                equation = {}
                for inner in range(node_sizes[target]):
                    coefficient = int(T[inner, col])
                    if coefficient:
                        index = variable(target, row, inner)
                        equation[index] = equation.get(index, 0)+coefficient
                for inner in range(node_sizes[source]):
                    coefficient = int(T[row, inner])
                    if coefficient:
                        index = variable(source, inner, col)
                        equation[index] = equation.get(index, 0)-coefficient
                if equation:
                    rows.append(equation)
        # The adjoint D block adds A_source*T^T-T^T*A_target=0.
        TT = T.T
        for row in range(node_sizes[source]):
            for col in range(node_sizes[target]):
                equation = {}
                for inner in range(node_sizes[source]):
                    coefficient = int(TT[inner, col])
                    if coefficient:
                        index = variable(source, row, inner)
                        equation[index] = equation.get(index, 0)+coefficient
                for inner in range(node_sizes[target]):
                    coefficient = int(TT[row, inner])
                    if coefficient:
                        index = variable(target, inner, col)
                        equation[index] = equation.get(index, 0)-coefficient
                if equation:
                    rows.append(equation)

    constraint = sp.MutableSparseMatrix(len(rows), 360, {})
    for row_index, equation in enumerate(rows):
        for column, coefficient in equation.items():
            constraint[row_index, column] = coefficient
    rank = int(constraint.rank())
    return rank, 360-rank, len(rows)


candidate_records = []
all_matrix_gates = True
all_gauge_equivalent = True
for reading_index, record in enumerate(eligible_readings):
    Gamma, gamma_diagonal = grading(record)
    links = tuple(sorted(record["links"]))
    constraint_rank, commutant_dimension, constraint_rows = (
        central_link_constraint_rank(links)
    )
    sign_diracs = {}
    for signs in product((1, -1), repeat=len(links)):
        link_signs = dict(zip(links, signs))
        D, half = build_dirac(record, link_signs)
        self_adjoint = (D-D.T).nnz == 0
        odd = (Gamma@D+D@Gamma).nnz == 0
        j_real = (J@D-D@J).nnz == 0
        first_order_support = (
            all(block["changed_side"] in ("left", "right")
                and int(hom_gram[block["link"][0],
                                 block["link"][1]]) == 1
                for block in record["legal_blocks"])
            and {block["link"] for block in record["legal_blocks"]
                 if block["changed_side"] == "left"} == set(links)
            and {block["link"] for block in record["legal_blocks"]
                 if block["changed_side"] == "right"} == set(links)
            and exact_first_order_tensor_support(record, link_signs, half)
        )
        nonzero_forms = constraint_rank > 0 and D.nnz > 0
        matrix_gates = (self_adjoint and odd and j_real
                        and first_order_support and nonzero_forms)
        all_matrix_gates &= matrix_gates
        sign_diracs[signs] = D
        candidate_records.append({
            "reading_index": reading_index,
            "priority": list(record["priority"]),
            "u_direction": record["u_direction"],
            "v_direction": record["v_direction"],
            "node_order": [nodes[node]["name"] for node in record["order"]],
            "links": [list(link) for link in links],
            "signs": list(signs),
            "D_nnz": int(D.nnz),
            "self_adjoint": self_adjoint,
            "odd": odd,
            "J_real": j_real,
            "first_order_by_exact_tensor_support": first_order_support,
            "commutator_map_rank": constraint_rank,
            "algebra_commutant_dimension": commutant_dimension,
            "connected": commutant_dimension == 1,
            "nonzero_one_form_witness": nonzero_forms,
            "constraint_rows": constraint_rows,
        })

    # Tree signs are removable by a representation-preserving cellwise gauge
    # U_(i,j)=g_i*g_j.  Construct and check all eight matrices explicitly.
    all_positive = tuple(1 for _ in links)
    reference = sign_diracs[all_positive]
    for signs, candidate in sign_diracs.items():
        edge_sign = dict(zip(links, signs))
        graph = {node: [] for node in range(4)}
        for left, right in links:
            graph[left].append((right, edge_sign[(left, right)]))
            graph[right].append((left, edge_sign[(left, right)]))
        gauges = {0: 1}
        queue = deque([0])
        while queue:
            node = queue.popleft()
            for neighbor, sign in graph[node]:
                if neighbor not in gauges:
                    gauges[neighbor] = gauges[node]*sign
                    queue.append(neighbor)
        gauge_diagonal = np.empty(dimension, dtype=np.int8)
        for cell in cells:
            start = cell_offsets[cell]
            stop = start+node_sizes[cell[0]]*node_sizes[cell[1]]
            gauge_diagonal[start:stop] = gauges[cell[0]]*gauges[cell[1]]
        Gauge = sparse.diags(gauge_diagonal, dtype=np.int8, format="csr")
        if (Gauge@reference@Gauge-candidate).nnz != 0:
            all_gauge_equivalent = False

check("all 32 eligible normalized rook operators pass matrix-level gates",
      len(candidate_records) == 32 and all_matrix_gates,
      "4 readings * 8 signs; every odd block equals T tensor I or I tensor T")
check("all link-sign choices are exact cellwise orthogonal gauges",
      all_gauge_equivalent)

commutant_dimensions = Counter(
    record["algebra_commutant_dimension"] for record in candidate_records
)
commutator_ranks = Counter(
    record["commutator_map_rank"] for record in candidate_records
)
connected_count = sum(record["connected"] for record in candidate_records)
check("every canonical equivariant rook operator fails connectedness",
      connected_count == 0
      and all(record["algebra_commutant_dimension"] > 1
              for record in candidate_records),
      f"commutant dimensions={dict(commutant_dimensions)}; "
      f"commutator ranks={dict(commutator_ranks)}")
check("every failed candidate nevertheless has a nonzero one-form witness",
      all(record["nonzero_one_form_witness"]
          for record in candidate_records))

payload = {
    "protocol_commit": "705edeb",
    "target_comparison_performed": False,
    "induced_A5_modules": {
        name: {a5_name: multiplicity for a5_name, multiplicity in zip(
            a5_names, induction_multiplicities[name]
        ) if multiplicity}
        for name in d5_names
    },
    "Hom_Gram": [[int(hom_gram[row, col]) for col in range(4)]
                  for row in range(4)],
    "reading_census": {
        "total": len(reading_records),
        "eligible_unique_Hom": len(eligible_readings),
        "ineligible": len(ineligible_readings),
        "hit_fraction": "4/8",
        "eligible_priority": "u_edge then v_ref",
        "ineligible_priority": "v_ref then u_edge",
    },
    "normalized_rook_census": {
        "candidate_count": len(candidate_records),
        "sign_choices_per_reading": 8,
        "all_signs_gauge_equivalent": all_gauge_equivalent,
        "matrix_gate_pass_count": sum(
            record["self_adjoint"] and record["odd"] and record["J_real"]
            and record["first_order_by_exact_tensor_support"]
            for record in candidate_records
        ),
        "connected_count": connected_count,
        "commutant_dimension_multiset": {
            str(key): value for key, value in sorted(commutant_dimensions.items())
        },
        "commutator_rank_multiset": {
            str(key): value for key, value in sorted(commutator_ranks.items())
        },
    },
    "candidates": candidate_records,
    "verdict": (
        "PATTERN: exactly 4/8 spectral readings, namely all edge-first "
        "readings, have one-dimensional equivariant Hom spaces on every "
        "required link. Their 32 normalized sign variants are gauge-"
        "equivalent and pass self-adjointness, KO6 reality, oddness, first "
        "order and nonzero forms, but every one has a non-scalar algebra "
        "commutant. DERIVED SCOPED KILL: no canonical A5-equivariant rook "
        "operator on this carrier is connected."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact structured equivariant-Dirac audit was written",
      OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("PATTERN: 4/8 readings have unique equivariant links (edge-first only).")
print("DERIVED SCOPED KILL: all normalized equivariant rook D fail connectedness.")
print("NO HESSIAN OR STANDARD-MODEL TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)

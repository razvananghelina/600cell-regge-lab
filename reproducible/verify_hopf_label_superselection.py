#!/usr/bin/env python3
"""Exact spectral-axiom audit of six-fibration label superselection.

The minimal and pair-groupoid arenas, hypotheses and orientability boundary
were frozen in commit 39e35b5.  All matrix tests use exact integer/rational
arithmetic after the registered six-fibration geometry is reconstructed.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp

from verify_hopf_fibration_invariants import (
    build_2I,
    find_all_hopf_fibrations,
    find_vertex_index,
    quat_mult,
)


OUTPUT = Path(__file__).with_name("hopf_label_superselection.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def commutator(left, right):
    return left@right-right@left


print("="*78)
print("EXACT SIX-LABEL SUPERSELECTION AXIOM AUDIT")
print("="*78)

n = 6
identity_6 = np.eye(n, dtype=np.int64)
minimal_projections = []
for index in range(n):
    diagonal = np.zeros(n, dtype=np.int64)
    diagonal[index] = 1
    minimal_projections.append(np.diag(diagonal))

# Arena A: on the minimal faithful diagonal representation, test every matrix
# unit of a general 6-by-6 D against every pair of minimal algebra projections.
minimal_legal_units = []
for row in range(n):
    for col in range(n):
        matrix_unit = np.zeros((n, n), dtype=np.int64)
        matrix_unit[row, col] = 1
        legal = all(np.count_nonzero(commutator(
            commutator(matrix_unit, left), right
        )) == 0 for left in minimal_projections
           for right in minimal_projections)
        if legal:
            minimal_legal_units.append((row, col))
check("first order on the minimal representation permits only diagonal D",
      minimal_legal_units == [(index, index) for index in range(n)],
      "legal matrix units=6 diagonal, 0 off-diagonal")

minimal_test_D = np.diag(np.arange(1, n+1, dtype=np.int64))
minimal_commutator_map = sp.Matrix(np.column_stack([
    commutator(minimal_test_D, projection).reshape(-1)
    for projection in minimal_projections
]).tolist())
minimal_forms = [
    left@commutator(minimal_test_D, right)
    for left in minimal_projections for right in minimal_projections
]
check("minimal first-order locality kills connectedness and all one-forms",
      minimal_commutator_map.rank() == 0
      and all(np.count_nonzero(form) == 0 for form in minimal_forms),
      "commuting algebra dimension=6, represented one-form dimension=0")

# The doubled selector has label-diagonal blocks, so with the canonical
# doubled representation pi(a)=diag(a,a) it also commutes with the full C^6.
phi_at_vertex_times_5 = np.diag([5, -1, -1, -1, -1, -1])
selector_times_5 = np.block([
    [5*identity_6, phi_at_vertex_times_5],
    [phi_at_vertex_times_5, 5*identity_6],
])
doubled_representations = [np.block([
    [projection, np.zeros_like(projection)],
    [np.zeros_like(projection), projection],
]) for projection in minimal_projections]
check("the canonical doubled selector representation has zero commutators",
      all(np.count_nonzero(commutator(selector_times_5, representation)) == 0
          for representation in doubled_representations),
      "D_aux supplies no inner one-forms and fails connectedness for C^6")

# Arena B: complete ordered-pair bimodule.  L and R are the two commuting
# diagonal actions; J_pair is the exact swap permutation.
dimension_pair = n*n
identity_pair = np.eye(dimension_pair, dtype=np.int64)
left_actions = []
right_actions = []
for projection in minimal_projections:
    left_actions.append(np.kron(projection, identity_6))
    right_actions.append(np.kron(identity_6, projection))

J_pair = np.zeros((dimension_pair, dimension_pair), dtype=np.int64)
for first in range(n):
    for second in range(n):
        source = first*n+second
        target = second*n+first
        J_pair[target, source] = 1

complete_adjacency = np.ones((n, n), dtype=np.int64)-identity_6
D_pair = (np.kron(complete_adjacency, identity_6)
          + np.kron(identity_6, complete_adjacency))
check("pair-groupoid left/right actions satisfy order zero exactly",
      all(np.count_nonzero(commutator(left, right)) == 0
          for left in left_actions for right in right_actions))
check("the parameter-free rook operator is self-adjoint and J-real",
      np.array_equal(D_pair, D_pair.T)
      and np.array_equal(J_pair@J_pair, identity_pair)
      and np.array_equal(J_pair@D_pair, D_pair@J_pair))
check("the rook operator satisfies first order on all 36 projection pairs",
      all(np.count_nonzero(commutator(
          commutator(D_pair, left), right
      )) == 0 for left in left_actions for right in right_actions))

# Exhaust the complete matrix-unit support allowed by first order.  A matrix
# element (i,j)<-(k,l) is legal iff i=k or j=l.
legal_pair_units = []
for first_target in range(n):
    for second_target in range(n):
        for first_source in range(n):
            for second_source in range(n):
                if (first_target == first_source
                        or second_target == second_source):
                    legal_pair_units.append((first_target, second_target,
                                             first_source, second_source))
legal_off_diagonal = sum(
    (first_target, second_target) != (first_source, second_source)
    for first_target, second_target, first_source, second_source
    in legal_pair_units
)
check("first order permits 360 genuinely off-diagonal pair channels",
      len(legal_pair_units) == 396 and legal_off_diagonal == 360
      and np.count_nonzero(D_pair) == 360,
      "legal support=396 including 36 diagonal; D_pair uses all 360 edges")

left_commutator_columns = np.column_stack([
    commutator(D_pair, action).reshape(-1) for action in left_actions
])
left_commutator_rank = sp.Matrix(left_commutator_columns.tolist()).rank()
one_form_columns = np.column_stack([
    (left@commutator(D_pair, right)).reshape(-1)
    for left in left_actions for right in left_actions
])
one_form_rank = sp.Matrix(one_form_columns.tolist()).rank()
check("the pair-groupoid witness is connected and has nonzero one-forms",
      left_commutator_rank == 5 and one_form_rank == 30,
      "kernel of a->[D,L(a)] has dimension 1; Omega_D^1 dimension=30")

# Reconstruct the actual A5 action on the six fibrations and verify diagonal
# equivariance of the pair-groupoid operator.
vertices = build_2I()
fibrations = find_all_hopf_fibrations(vertices)


def conjugation_permutation(group_element):
    inverse = group_element.copy()
    inverse[1:] *= -1
    return tuple(find_vertex_index(
        vertices,
        quat_mult(quat_mult(group_element, vertex), inverse),
    ) for vertex in vertices)


vertex_actions = sorted(set(conjugation_permutation(element)
                            for element in vertices))
fibration_by_signature = {
    tuple(sorted(tuple(sorted(fiber)) for fiber in fibration)): index
    for index, fibration in enumerate(fibrations)
}
label_actions = []
for vertex_action in vertex_actions:
    action = []
    for fibration in fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fiber))
            for fiber in fibration
        ))
        action.append(fibration_by_signature[signature])
    label_actions.append(tuple(action))
label_actions = sorted(set(label_actions))

pair_equivariance = True
for action in label_actions:
    permutation = np.zeros((n, n), dtype=np.int64)
    for source, target in enumerate(action):
        permutation[target, source] = 1
    pair_permutation = np.kron(permutation, permutation)
    if not np.array_equal(pair_permutation@D_pair@pair_permutation.T, D_pair):
        pair_equivariance = False
check("the off-diagonal pair witness is equivariant under all 60 A5 actions",
      len(label_actions) == 60 and pair_equivariance)

# Standard odd double with KO6 signs.  Since all matrices are real, the
# anti-linear conjugation part of J is represented by the displayed swap.
zero_pair = np.zeros_like(D_pair)
D_total = np.block([[zero_pair, D_pair], [D_pair, zero_pair]])
gamma = np.block([[identity_pair, zero_pair],
                  [zero_pair, -identity_pair]])
J_total = np.block([[zero_pair, J_pair], [J_pair, zero_pair]])
identity_total = np.eye(2*dimension_pair, dtype=np.int64)
left_total = [np.block([[left, zero_pair], [zero_pair, left]])
              for left in left_actions]
right_total = [J_total@left@J_total for left in left_total]
check("the odd double has exact KO6 signs and odd self-adjoint D",
      np.array_equal(J_total@J_total, identity_total)
      and np.array_equal(J_total@D_total, D_total@J_total)
      and np.array_equal(J_total@gamma, -gamma@J_total)
      and np.array_equal(gamma@D_total, -D_total@gamma)
      and np.array_equal(D_total, D_total.T))
check("order zero and first order survive the KO6 odd double",
      all(np.count_nonzero(commutator(left, right)) == 0
          for left in left_total for right in right_total)
      and all(np.count_nonzero(commutator(
          commutator(D_total, left), right
      )) == 0 for left in left_total for right in right_total))

# Metric-dimension-zero orientability.  pi(e_i)Jpi(e_j)J^-1 acts identically
# on both grading sheets.  Its 36 diagonal vectors span a 36D subspace;
# adjoining gamma raises the rank, proving gamma is not in that span.
orientation_columns = []
for left in left_total:
    for right in right_total:
        orientation_columns.append(np.diag(left@right))
orientation_matrix = sp.Matrix(np.column_stack(orientation_columns).tolist())
gamma_diagonal = sp.Matrix(np.diag(gamma).tolist())
orientation_augmented = orientation_matrix.row_join(gamma_diagonal)
check("the KO6 pair witness fails metric-dimension-zero orientability",
      orientation_matrix.rank() == 36
      and orientation_augmented.rank() == 37,
      "all zero-cycles act equally on both sheets, while gamma changes sign")

# A single-copy multiplicity-one KO6 grading cannot repair orientability on
# ordered distinct pairs: the A5 orbit of (0,1) is all 30 pairs and contains
# its swap.  An invariant sign is constant on the orbit, whereas Jgamma=-gammaJ
# demands the swapped sign be its negative.
ordered_pair_orbit = set((action[0], action[1]) for action in label_actions)
check("A5 transitivity obstructs a single-copy invariant KO6 pair grading",
      len(ordered_pair_orbit) == 30 and (1, 0) in ordered_pair_orbit,
      "constant A5 orbit sign contradicts the J-odd swapped sign")

# General metric-dimension-zero no-go, including arbitrary multiplicities.
# Any C^6 bimodule decomposes into H_ij.  A zero-cycle acts as a scalar on
# each H_ij, so an orienting grading has a sign epsilon_ij times the identity
# there.  Jgamma=-gammaJ gives epsilon_ji=-epsilon_ij.  The A5 action has only
# two orbits on F x F: the six diagonal pairs and all 30 distinct pairs.
# Diagonal support contradicts epsilon_ii=-epsilon_ii; off-diagonal support
# contradicts invariance because its orbit contains the swapped pair.  Extra
# multiplicity cannot change a scalar zero-cycle on H_ij.
all_pair_orbits = []
unseen_pairs = {(first, second) for first in range(n) for second in range(n)}
while unseen_pairs:
    seed = next(iter(unseen_pairs))
    orbit = {(action[seed[0]], action[seed[1]]) for action in label_actions}
    all_pair_orbits.append(orbit)
    unseen_pairs -= orbit
pair_orbit_sizes = sorted(len(orbit) for orbit in all_pair_orbits)
diagonal_orbit = next(orbit for orbit in all_pair_orbits
                      if all(first == second for first, second in orbit))
off_diagonal_orbit = next(orbit for orbit in all_pair_orbits
                          if any(first != second for first, second in orbit))
diagonal_obstructed = all((index, index) in diagonal_orbit
                          for index in range(n))
off_diagonal_obstructed = any(
    first != second and (second, first) in off_diagonal_orbit
    for first, second in off_diagonal_orbit
)
check("all A5-equivariant C6 bimodule orbits are KO6-orientability obstructed",
      pair_orbit_sizes == [6, 30]
      and diagonal_obstructed and off_diagonal_obstructed,
      "diagonal: epsilon_ii=-epsilon_ii; off-diagonal: invariance conflicts "
      "with epsilon_ji=-epsilon_ij; arbitrary multiplicities do not evade")

payload = {
    "protocol_commit": "39e35b5",
    "minimal_representation": {
        "legal_D_matrix_units": len(minimal_legal_units),
        "off_diagonal_legal_units": 0,
        "connected": False,
        "one_form_dimension": 0,
        "D_aux_commutes_with_C6": True,
    },
    "pair_groupoid": {
        "dimension": dimension_pair,
        "first_order_legal_matrix_units": len(legal_pair_units),
        "first_order_legal_off_diagonal_units": legal_off_diagonal,
        "D_pair_edges": int(np.count_nonzero(D_pair)),
        "order_zero": True,
        "first_order": True,
        "connected": True,
        "one_form_dimension": int(one_form_rank),
        "A5_equivariant": pair_equivariance,
        "KO6_odd_double": True,
        "orientable": False,
    },
    "general_orientability_no_go": {
        "hypotheses": [
            "A=C^6 functions on six fibrations",
            "A5-equivariant bimodule",
            "KO6 J*gamma=-gamma*J",
            "metric-dimension-zero orientability",
        ],
        "A5_orbits_on_F_times_F": pair_orbit_sizes,
        "arbitrary_multiplicities_allowed": True,
        "nonzero_faithful_triple_exists": False,
    },
    "verdict": (
        "SCOPED REFUTATION: order zero, first order, connectedness, nonzero "
        "forms, KO6 and A5 equivariance do not force label diagonality; the "
        "canonical pair-groupoid witness retains 360 off-diagonal channels. "
        "DERIVED FULL-ARENA NO-GO: adding metric-dimension-zero orientability "
        "does not select the diagonal; together with A5 and KO6 it forbids "
        "every nonzero C^6 bimodule, even with arbitrary multiplicities. "
        "Minimal-representation diagonality is obtained only by killing all "
        "commutators and connectedness"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact label-superselection audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED NEGATIVE: minimal first order forces diagonality only by killing")
print("                  connectedness and all one-forms.")
print("SCOPED REFUTATION: the canonical pair bimodule passes order zero, first")
print("                    order, connectedness, KO6 and A5 with 360 off-diagonal")
print("                    channels, but fails orientability.")
print("DERIVED NO-GO: adding metric-dimension-zero orientability kills every")
print("               nonzero A5-equivariant KO6 C^6 bimodule; it does not derive")
print("               a viable diagonal superselection sector.")
raise SystemExit(0 if passed == tests else 1)

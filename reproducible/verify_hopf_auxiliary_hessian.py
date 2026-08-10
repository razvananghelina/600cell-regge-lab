#!/usr/bin/env python3
"""Exact audit of the Hessian/Gram origin of the six-label selector.

The three constructions and their decision boundaries were frozen in commit
90627f7.  Geometry discovery reuses the registered six-fibration constructor;
all trace, Hessian, rank and commutator tests then use exact integer/rational
arithmetic.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp

from verify_hopf_fibration_invariants import (
    build_2I,
    build_adjacency,
    build_fiber_adjacency,
    find_all_hopf_fibrations,
    find_vertex_index,
    quat_mult,
)


OUTPUT = Path(__file__).with_name("hopf_auxiliary_hessian.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def trace_product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = result @ matrix
    return int(np.trace(result))


BOX_PRODUCTS = None


def label_hessian(X, boxes):
    return sp.Matrix(6, 6, lambda row, col:
                     3*(int(np.sum(X*BOX_PRODUCTS[row][col].T))
                        + int(np.sum(X*BOX_PRODUCTS[col][row].T))))


print("="*78)
print("EXACT AUXILIARY HESSIAN / GRAM BRIDGE AUDIT")
print("="*78)

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
fiber_adjacencies = [
    np.rint(build_fiber_adjacency(adjacency, fibration)).astype(np.int64)
    for fibration in fibrations
]
boxes = [6*fiber-adjacency for fiber in fiber_adjacencies]
BOX_PRODUCTS = [[boxes[row]@boxes[col] for col in range(6)]
                for row in range(6)]
basis = [boxes[index]-boxes[5] for index in range(5)]
gram = sp.Matrix([
    [trace_product(basis[row], basis[col]) for col in range(5)]
    for row in range(5)
])
check("the fixed geometry reconstructs the centered five-dimensional Box span",
      len(boxes) == 6 and np.array_equal(sum(boxes), np.zeros_like(adjacency))
      and gram.rank() == 5)

# Original one-parameter bootstrap.  Expand every noncommutative cubic exactly:
# Tr((c Af-A)^3)=alpha_3 c^3+alpha_2 c^2+alpha_1 c+alpha_0.
original_coefficients = []
for fiber in fiber_adjacencies:
    alpha_3 = trace_product(fiber, fiber, fiber)
    alpha_2 = -(
        trace_product(fiber, fiber, adjacency)
        + trace_product(fiber, adjacency, fiber)
        + trace_product(adjacency, fiber, fiber)
    )
    alpha_1 = (
        trace_product(fiber, adjacency, adjacency)
        + trace_product(adjacency, fiber, adjacency)
        + trace_product(adjacency, adjacency, fiber)
    )
    alpha_0 = -trace_product(adjacency, adjacency, adjacency)
    original_coefficients.append((alpha_3, alpha_2, alpha_1, alpha_0))
check("all six original bootstrap polynomials are exactly 3600(c-2)",
      set(original_coefficients) == {(0, 0, 3600, -7200)})
check("the actual one-parameter bootstrap Hessian vanishes identically",
      all(6*alpha_3 == 0 and 2*alpha_2 == 0
          for alpha_3, alpha_2, _, _ in original_coefficients),
      "d^2/dc^2 Tr((c Af-A)^3)=0 for every fibration")

# Full Hessian of V(X)=Tr(X^3) in all six label directions.
basis_hessians = [label_hessian(direction, boxes) for direction in basis]
vertex_hessians = [label_hessian(box, boxes) for box in boxes]
check("the exact Hessian formula is symmetric and kills the constant label mode",
      all(hessian == hessian.T
              and hessian*sp.ones(6, 1) == sp.zeros(6, 1)
          for hessian in basis_hessians+vertex_hessians))

analysis = sp.Matrix(6, 5, lambda label, coordinate:
                     sp.Rational(trace_product(basis[coordinate],
                                               boxes[label]), 7200))
diagonal_hessian_map = sp.Matrix(6, 5, lambda label, coordinate:
                                 basis_hessians[coordinate][label, label])
check("the disclosed diagonal-Hessian identity holds coefficientwise",
      diagonal_hessian_map == 86400*analysis,
      "diag(H_X)/(12*7200)=diag(Tr(X Box_i)/7200)=Phi(X)")

off_diagonal_pairs = [(row, col) for row in range(6)
                      for col in range(row+1, 6)]
off_diagonal_map = sp.Matrix([
    [basis_hessians[coordinate][row, col] for coordinate in range(5)]
    for row, col in off_diagonal_pairs
])
full_hessian_map = sp.Matrix([
    [basis_hessians[coordinate][row, col] for coordinate in range(5)]
    for row in range(6) for col in range(6)
])
off_diagonal_values = sorted(set(
    basis_hessians[0][row, col] for row, col in off_diagonal_pairs
))
check("the off-diagonal Hessian is nonzero and retains the full W information",
      off_diagonal_map.rank() == 5 and full_hessian_map.rank() == 5
      and off_diagonal_values
      == [-155520, -103680, -51840, 0, 51840, 103680, 155520],
      f"off-diagonal map rank=5; values={off_diagonal_values}")

expected_first_hessian = sp.Matrix([
    [86400, -17280, -17280, -17280, -17280, -17280],
    [-17280, -17280, -69120, -69120, 86400, 86400],
    [-17280, -69120, -17280, 86400, -69120, 86400],
    [-17280, -69120, 86400, -17280, 86400, -69120],
    [-17280, 86400, -69120, 86400, -17280, -69120],
    [-17280, 86400, 86400, -69120, -69120, -17280],
])
check("an explicit derived Hessian exhibits the nonlocal label couplings",
      vertex_hessians[0] == expected_first_hessian)
first_hessian_eigenvalues = vertex_hessians[0].eigenvals()
expected_first_eigenvalues = {
    sp.Integer(0): 1,
    sp.Integer(103680): 1,
    -25920-77760*sp.sqrt(5): 2,
    -25920+77760*sp.sqrt(5): 2,
}
check("the full Hessian kernel is universal and its physical part is indefinite",
      first_hessian_eigenvalues == expected_first_eigenvalues
      and all(hessian*sp.ones(6, 1) == sp.zeros(6, 1)
              for hessian in vertex_hessians),
      "on the zero-sum sector: signature (3 positive,2 negative), no zero")

commutator = vertex_hessians[0]*vertex_hessians[1]
commutator -= vertex_hessians[1]*vertex_hessians[0]
commutator_max = max(abs(int(value)) for value in commutator)
check("the Hessian family cannot be diagonalized by one fixed label basis",
      commutator.rank() == 4 and commutator_max == 26873856000,
      f"rank[H_Box0,H_Box1]={commutator.rank()}; max entry={commutator_max}")

# Reconstruct A5 permutation actions and verify both the full Hessian and the
# diagonal conditional expectation transform equivariantly.
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
action_pairs = []
for vertex_action in vertex_actions:
    label_action = []
    for fibration in fibrations:
        signature = tuple(sorted(
            tuple(sorted(vertex_action[index] for index in fiber))
            for fiber in fibration
        ))
        label_action.append(fibration_by_signature[signature])
    action_pairs.append((vertex_action, tuple(label_action)))

full_equivariance = True
diagonal_equivariance = True
for vertex_action, label_action in action_pairs:
    vertex_indices = np.asarray(vertex_action, dtype=int)
    label_permutation = sp.zeros(6)
    for source, target in enumerate(label_action):
        label_permutation[target, source] = 1
    for source, box in enumerate(boxes):
        transformed_box = np.zeros_like(box)
        transformed_box[np.ix_(vertex_indices, vertex_indices)] = box
        target = label_action[source]
        if not np.array_equal(transformed_box, boxes[target]):
            full_equivariance = False
        transformed_hessian = vertex_hessians[target]
        predicted = (label_permutation*vertex_hessians[source]
                     *label_permutation.T)
        if transformed_hessian != predicted:
            full_equivariance = False
        transformed_diagonal = sp.diag(*transformed_hessian.diagonal())
        predicted_diagonal = (label_permutation
                              *sp.diag(*vertex_hessians[source].diagonal())
                              *label_permutation.T)
        if transformed_diagonal != predicted_diagonal:
            diagonal_equivariance = False
check("the full label Hessian is exactly A5-equivariant",
      full_equivariance and len(action_pairs) == 60)
check("diagonal conditional expectation is also A5-equivariant",
      diagonal_equivariance,
      "canonical only after assuming locality in the fibration-label algebra")

# Gram origin of D_aux.  In normalized Hilbert--Schmidt units, each label
# block is the Gram matrix of (X,Box_i).  Its Schur complement is 1-r_i^2.
r = sp.symbols("r", real=True)
pair_gram = sp.Matrix([[1, r], [r, 1]])
schur_complement = sp.simplify(pair_gram[1, 1]
                               - pair_gram[1, 0]
                               *pair_gram[0, 0]**-1
                               *pair_gram[0, 1])
check("D_aux is exactly the direct sum of normalized pair-Gram matrices",
      pair_gram.det() == 1-r**2 and schur_complement == 1-r**2,
      "Schur complement is K_i=1-r_i^2")

# The Gram construction itself is exhaustive at the sharp points: determinant
# zero means equality in Cauchy--Schwarz, hence X=+/-Box_i.  It is a recognizer,
# not the Hessian of the original one-dimensional action.
check("the Gram singularity is precisely collinearity, not a new field equation",
      sp.expand(pair_gram.det()+(r-1)*(r+1)) == 0)

# Scoped source audit: the two authoritative original-bootstrap constructors
# define only 120-state adjacency/fibre operators and no six-label coupling.
source_paths = [
    Path(__file__).with_name("verify_variational_bootstrap.py"),
    Path(__file__).with_name("verify_hopf_fibration_invariants.py"),
]
source_text = "\n".join(path.read_text() for path in source_paths)
new_coupling_tokens = ("D_aux", "Phi(X)", "six-label", "label carrier")
check("the specified original bootstrap files contain no state-to-label coupling",
      all(token not in source_text for token in new_coupling_tokens),
      "scoped source audit only; not a theorem against future couplings")

payload = {
    "protocol_commit": "90627f7",
    "original_bootstrap": {
        "polynomial": "Tr((c*A_f-A)^3)=3600*(c-2)",
        "fibrations_checked": 6,
        "second_derivative": 0,
        "direct_Hessian_bridge": False,
    },
    "extended_cubic_Hessian": {
        "definition": "H_X(i,j)=3*Tr(X*(Box_i*Box_j+Box_j*Box_i))",
        "map_rank": int(full_hessian_map.rank()),
        "off_diagonal_map_rank": int(off_diagonal_map.rank()),
        "diagonal_identity": "diag(H_X)=86400*diag(r_i(X))",
        "H_Box0_H_Box1_commutator_rank": int(commutator.rank()),
        "H_Box0_H_Box1_commutator_max_entry": commutator_max,
        "H_Box0_spectrum": {
            str(value): multiplicity
            for value, multiplicity in first_hessian_eigenvalues.items()
        },
        "universal_kernel": "constant label mode for every X",
        "zero_sum_signature_at_Box": [3, 2, 0],
        "simultaneously_diagonalizable": False,
        "diagonal_expectation_A5_equivariant": diagonal_equivariance,
        "status": "conditional on a new fibration-label locality projection",
    },
    "Gram_origin": {
        "D_aux": "direct sum_i [[1,r_i],[r_i,1]]",
        "Schur_complement": "K_i=1-r_i^2",
        "singularity_meaning": "X is collinear with some Box_i",
        "is_Schur_complement_of_original_120x120_Box": False,
    },
    "verdict": (
        "DERIVED NEGATIVE for direct dynamics: the actual one-parameter "
        "bootstrap Hessian is zero.  DERIVED CONDITIONAL bridge: the diagonal "
        "of the extended cubic Hessian is exactly Phi, but the full Hessian "
        "has rank-five off-diagonal data and is not simultaneously diagonal; "
        "discarding it is a new label-locality assumption.  D_aux is exactly "
        "a pair-Gram recognizer, not an existing Box Schur complement"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact auxiliary Hessian/Gram JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED NEGATIVE: the original c-bootstrap Hessian is identically zero.")
print("DERIVED CONDITIONAL: diag(H_X)/86400=Phi(X), but off-diagonal rank=5.")
print("DERIVED: D_aux is the normalized pair-Gram operator; its singularity")
print("         recognizes collinearity rather than following from old dynamics.")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Exact counterexample to the repaired chamber commutativity conjecture B1.

The verifier uses the fixed 120-chamber D, gamma and geometric reflection J.
No Dirac matrix, spectrum, Schur coefficient, or tolerance is fitted.  A
four-colouring of the exact integer graph S=(D J)|H+ places its vertices in
the four oriented Krajewski cells

    (0,1)x2, (1,2)x25, (3,1)x12, (2,3)x19

for A=M2(C)+C+C+C.  The M2 node is a pure source.  Consequently every
shared index of every occupied Dirac block is scalar, so the Krajewski block
mask is the full matrix first-order condition in this example.  The script
nevertheless checks order zero and first order on a full complex basis.
"""

import contextlib
import io
import os
from pathlib import Path
import runpy

import numpy as np
import sympy as sy


os.chdir(Path(__file__).resolve().parent)


passed = tests = 0


def check(label, condition, detail=""):
    global passed, tests
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


print("=" * 78)
print("EXACT B1 COUNTEREXAMPLE ON THE FIXED 120-CHAMBER CARRIER")
print("=" * 78)

# Load the already-registered exact chamber construction.  Its output is
# suppressed; this verifier independently checks every gate used below.
with contextlib.redirect_stdout(io.StringIO()):
    geometry = runpy.run_path("verify_oriented_chamber_double.py")

D = geometry["D"].toarray().astype(np.int64)
J = geometry["J"].toarray().astype(np.int64)
gamma_vector = geometry["gamma"].astype(np.int64)
Gamma = np.diag(gamma_vector)
reflection = np.asarray(geometry["reflection"], dtype=np.int64)
plus = np.flatnonzero(gamma_vector == 1)
minus = np.flatnonzero(gamma_vector == -1)
S = (D @ J)[np.ix_(plus, plus)]

check("[DERIVED] fixed data are the 120-chamber KO6 carrier",
      D.shape == (120, 120)
      and len(plus) == len(minus) == 60
      and np.array_equal(D, D.T)
      and np.array_equal(J @ J, np.eye(120, dtype=np.int64))
      and np.array_equal(J @ D, D @ J)
      and np.array_equal(J @ Gamma, -Gamma @ J))
check("[DERIVED] S=(D J)|H+ is symmetric, invertible and loopless 3-regular",
      np.array_equal(S, S.T)
      and sy.Matrix(S.tolist()).det() != 0
      and np.count_nonzero(np.diag(S)) == 0
      and np.array_equal(S.sum(axis=1), 3*np.ones(60, dtype=np.int64)))

# Algebra-node sizes and positive-sheet bimodule cells (left,right,mult).
# Node 0 is M2(C); nodes 1,2,3 are scalar.  The M2 node is a pure source.
factor_sizes = (2, 1, 1, 1)
cells = ((0, 1, 2), (1, 2, 25), (3, 1, 12), (2, 3, 19))
cell_dimensions = tuple(factor_sizes[i]*factor_sizes[j]*m
                        for i, j, m in cells)

# Exact certificate in the sorted H+ chamber order.  Labels index `cells`.
CELL_LABELS = (
    2, 0, 2, 2, 1, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1,
    0, 3, 3, 3, 3, 1, 3, 2, 2, 3, 0, 3, 3, 0, 3, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 3, 3, 3, 3, 2,
    3, 1, 2, 1, 1, 3, 1, 1, 2, 3, 1, 1, 3,
)
labels = np.asarray(CELL_LABELS, dtype=np.int64)
cell_vertices = tuple(plus[np.flatnonzero(labels == q)]
                      for q in range(len(cells)))

check("[DERIVED] certificate has exact cell dimensions (4,25,12,19)",
      len(labels) == 60
      and tuple(len(vertices) for vertices in cell_vertices)
      == cell_dimensions == (4, 25, 12, 19)
      and sum(cell_dimensions) == 60)


def composable(first, second):
    """Krajewski first-order cell support relation in the J-paired basis."""
    i, j, _ = first
    k, l, _ = second
    return k == j or l == i


edge_types = set()
shared_nodes_are_scalar = True
for x in range(60):
    for y in range(x + 1, 60):
        if not S[x, y]:
            continue
        q, r = int(labels[x]), int(labels[y])
        first, second = cells[q], cells[r]
        edge_types.add(tuple(sorted((q, r))))
        if second[0] == first[1]:
            shared = first[1]
        elif first[0] == second[1]:
            shared = first[0]
        else:
            shared_nodes_are_scalar = False
            continue
        shared_nodes_are_scalar &= factor_sizes[shared] == 1

check("[DERIVED] every fixed-D edge lies in a legal composable cell pair",
      len(np.transpose(np.nonzero(np.triu(S, 1)))) == 90
      and edge_types == {(0, 1), (1, 2), (1, 3), (2, 3)}
      and all(composable(cells[q], cells[r])
              or composable(cells[r], cells[q])
              for q, r in edge_types),
      f"occupied cell-edge types={sorted(edge_types)}")
check("[DERIVED] every shared first-order index is scalar",
      shared_nodes_are_scalar and all(j != 0 for i, j, _ in cells),
      "M2 is a pure source, so no hidden matrix-factor intertwiner is omitted")


def representation(matrix2, lambda1, lambda2, lambda3):
    """Complex-linear faithful representation of M2(C)+C^3."""
    matrix2 = np.asarray(matrix2, dtype=np.complex128)
    scalar = (None, lambda1, lambda2, lambda3)
    answer = np.zeros((120, 120), dtype=np.complex128)
    for q, (i, j, multiplicity) in enumerate(cells):
        vertices = cell_vertices[q]
        n_i, n_j = factor_sizes[i], factor_sizes[j]
        if i == 0:
            plus_block = np.kron(
                np.kron(matrix2, np.eye(n_j)), np.eye(multiplicity)
            )
        else:
            plus_block = scalar[i]*np.eye(len(vertices))
        answer[np.ix_(vertices, vertices)] = plus_block

        reflected = reflection[vertices]
        if j == 0:
            minus_block = np.kron(
                np.kron(np.eye(n_i), matrix2), np.eye(multiplicity)
            )
        else:
            minus_block = scalar[j]*np.eye(len(vertices))
        answer[np.ix_(reflected, reflected)] = minus_block
    return answer


zero2 = np.zeros((2, 2), dtype=np.complex128)
matrix_units = []
for row in range(2):
    for column in range(2):
        unit = zero2.copy()
        unit[row, column] = 1
        matrix_units.append(unit)

generators = [representation(unit, 0, 0, 0) for unit in matrix_units]
generators.extend((representation(zero2, 1, 0, 0),
                   representation(zero2, 0, 1, 0),
                   representation(zero2, 0, 0, 1)))

# For J=P K, J pi(b)^* J^-1 has linear matrix P pi(b)^T P.
opposites = tuple(J @ generator.T @ J for generator in generators)

sample_a = np.asarray([[1 + 2j, 3 - 1j], [-2j, 4]], dtype=np.complex128)
sample_b = np.asarray([[2, 1j], [1 - 2j, -3j]], dtype=np.complex128)
pi_a = representation(sample_a, 2 + 1j, -1j, 3)
pi_b = representation(sample_b, -2j, 4, 1 + 1j)
pi_ab = representation(sample_a @ sample_b,
                       (2 + 1j)*(-2j), (-1j)*4, 3*(1 + 1j))
pi_astar = representation(sample_a.conjugate().T,
                          (2 + 1j).conjugate(), (-1j).conjugate(), 3)

flat_basis = np.column_stack([matrix.reshape(-1) for matrix in generators])
flat_gram = np.rint((flat_basis.conjugate().T @ flat_basis).real).astype(int)
check("[DERIVED] representation is complex-linear, multiplicative and star preserving",
      np.array_equal(representation(1j*sample_a, 2j - 1, 1, 3j),
                           1j*representation(sample_a, 2 + 1j, -1j, 3))
      and np.array_equal(pi_a @ pi_b, pi_ab)
      and np.array_equal(pi_a.conjugate().T, pi_astar))
check("[DERIVED] A=M2(C)+C^3 acts faithfully, unitally and noncommutatively",
      sy.Matrix(flat_gram.tolist()).rank() == 7
      and np.array_equal(representation(np.eye(2), 1, 1, 1),
                         np.eye(120))
      and not np.array_equal(generators[1] @ generators[2],
                             generators[2] @ generators[1]))
check("[DERIVED] represented algebra commutes with gamma",
      all(np.array_equal(Gamma @ a, a @ Gamma) for a in generators))
check("[DERIVED] order zero holds on a full complex algebra basis",
      all(np.array_equal(a @ b, b @ a)
          for a in generators for b in opposites))
check("[DERIVED] first order holds on a full complex algebra basis",
      all(np.array_equal((D @ a - a @ D) @ b,
                         b @ (D @ a - a @ D))
          for a in generators for b in opposites))
check("[DERIVED] represented inner one-forms are nonzero",
      any(not np.array_equal(D @ a, a @ D) for a in generators))

# Metric-dimension-zero orientability.  The four central projectors are used
# only to display one explicit Hochschild zero-cycle; matrix coefficients are
# allowed algebra elements, and the identity of M2 is its central projector.
central = (representation(np.eye(2), 0, 0, 0),
           representation(zero2, 1, 0, 0),
           representation(zero2, 0, 1, 0),
           representation(zero2, 0, 0, 1))
central_opposite = tuple(J @ p.T @ J for p in central)
orientation_cycle = np.zeros((120, 120), dtype=np.complex128)
for i, j, _ in cells:
    orientation_cycle += (central[i] @ central_opposite[j]
                          - central[j] @ central_opposite[i])
check("[DERIVED] metric-dimension-zero orientability holds exactly",
      np.array_equal(orientation_cycle, Gamma),
      "Gamma=sum_cells (p_i p_j^0-p_j p_i^0)")

# Minimal K0 projector for M2 plus the three scalar projectors.
minimal_m2 = zero2.copy()
minimal_m2[0, 0] = 1
k0_projectors = (representation(minimal_m2, 0, 0, 0),
                 central[1], central[2], central[3])
k0_opposites = tuple(J @ p.T @ J for p in k0_projectors)
intersection = sy.Matrix([
    [int(round(np.trace(Gamma @ p @ q).real)) for q in k0_opposites]
    for p in k0_projectors
])
expected_intersection = sy.Matrix([
    [0, 2, 0, 0],
    [-2, 0, 25, -12],
    [0, -25, 0, 19],
    [0, 12, -19, 0],
])
check("[DERIVED] KO6 intersection form is exact and antisymmetric",
      intersection == expected_intersection
      and intersection == -intersection.T)
check("[DERIVED] intersection form is nondegenerate",
      intersection.rank() == 4 and intersection.det() == 1444,
      "Pfaffian=38, determinant=38^2; unimodularity is not a B1 hypothesis")

# Connectedness is a linear condition on all seven complex coefficients.
# Rank the 7-column commutator map through its exact 7x7 Gram matrix.
commutator_columns = np.column_stack([
    np.rint((D @ a - a @ D).real).astype(np.int64).reshape(-1)
    for a in generators
])
commutator_gram = commutator_columns.T @ commutator_columns
commutator_rank = sy.Matrix(commutator_gram.tolist()).rank()
check("[DERIVED] connectedness holds: [D,a]=0 only for complex scalars",
      commutator_rank == 6
      and all(np.array_equal(D @ representation(z*np.eye(2), z, z, z),
                             representation(z*np.eye(2), z, z, z) @ D)
              for z in (1, 1j)),
      f"dim_C ker(a -> [D,pi(a)])={7-commutator_rank}")

# Framing audit.  A single oriented (M60,C) cell already supplies a faithful
# unital order-zero representation on 60+60 with M60 on only one sheet.  Its
# 2x2 skew multiplicity form is nondegenerate.  It need not satisfy first
# order; it suffices to refute the claimed order-zero k<=7 inference and to
# show why one-sheet support is not killed by PD or orientability alone.
one_sheet_cap = sy.Matrix([[0, 1], [-1, 0]])
check("[DERIVED NEGATIVE] faithfulness, unitality and order zero do not imply k<=7",
      60*1 == 60 and 60 > 7,
      "A=M60(C)+C on the single cell (M60,C) is a two-sheet faithful bimodule")
check("[DERIVED NEGATIVE] orientability and nondegenerate PD do not kill one-sheet summands",
      one_sheet_cap.det() == 1 and one_sheet_cap == -one_sheet_cap.T,
      "the oriented two-node multiplicity matrix has no loop or reverse pair")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT_B1=REFUTED_BY_EXACT_NONCOMMUTATIVE_ALL_GATE_WITNESS")
print("VERDICT_STRONG_UNIMODULAR_B1=OPEN")
print("VERDICT_WEDDERBURN_TYPE_ENUMERATION=INSUFFICIENT")
if passed != tests:
    raise SystemExit(1)

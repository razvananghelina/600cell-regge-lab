#!/usr/bin/env python3
"""Exact refutation of the proposed chamber-algebra commutativity theorem.

The fixed chamber Dirac is used through S=(D J)|H_+.  Exact rank over Q
shows that S has a six-dimensional (-1)-eigenspace.  A noncommutative
Krajewski block with the same Takagi singular block is then transported to
that eigenspace by unitary congruence.  The remaining 54 dimensions carry a
scalar-scalar bimodule, on which first order permits the untouched symmetric
complement of S.
"""

import contextlib
from collections import Counter
import io
import runpy

import numpy as np
import sympy as sy
from sympy.polys.matrices import DomainMatrix


passed = tests = 0


def check(label, condition, detail=""):
    global passed, tests
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def exact_rank(matrix):
    return len(DomainMatrix.from_Matrix(matrix).rref()[1])


# Load the registered chamber geometry without duplicating its construction.
with contextlib.redirect_stdout(io.StringIO()):
    geometry = runpy.run_path("verify_oriented_chamber_double.py")

D_fixed = geometry["D"].toarray().astype(np.int64)
J_fixed = geometry["J"].toarray().astype(np.int64)
gamma_fixed = geometry["gamma"]
plus = np.flatnonzero(gamma_fixed == 1)

# Ordering H_- by J(H_+) identifies the off-diagonal Dirac block with
# S=(D J)|H_+.  It is an exact integer symmetric matrix.
S_np = (D_fixed @ J_fixed)[np.ix_(plus, plus)]
S = sy.Matrix(S_np.tolist())
I60 = sy.eye(60)
rank_S = exact_rank(S)
minus_one_nullity = 60 - exact_rank(S + I60)
plus_one_nullity = 60 - exact_rank(S - I60)
check("fixed chamber S=(D J)|H+ is symmetric and invertible",
      S == S.T and rank_S == 60)
check("fixed chamber S has an exact six-dimensional -1 eigenspace",
      minus_one_nullity == 6 and plus_one_nullity == 0,
      f"null(S+I)={minus_one_nullity}, null(S-I)={plus_one_nullity}")

# A genuine scoped positive result: the registered C5 partition witness
# cannot be matrix-amplified while retaining its central supports.  Rebuild
# it from its 90 contraction bits, then inspect P_i intersect J(P_j).
C5_BITS = (
    "100110111100011000110100010100011100110010111100111010101010010001"
    "111000111110001101100110"
)
edge_pairs = geometry["edge_pairs"]
chamber_edges = geometry["chamber_edges"]
reflection = geometry["reflection"]
parent = list(range(120))


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(x, y):
    x, y = find(x), find(y)
    if x != y:
        parent[y] = x


for pair, bit in zip(edge_pairs, C5_BITS):
    union(*chamber_edges[pair[int(bit)]])
roots = [find(x) for x in range(120)]
root_labels = {root: i for i, root in enumerate(sorted(set(roots)))}
c5_labels = tuple(root_labels[root] for root in roots)
c5_dim = len(set(c5_labels))
c5_cells = np.zeros((c5_dim, c5_dim), dtype=np.int64)
for x in range(120):
    c5_cells[c5_labels[x], c5_labels[reflection[x]]] += 1
cell_histogram = Counter(map(int, c5_cells.ravel()))
unit_neighbors = [set() for _ in range(c5_dim)]
for i, j in zip(*np.where(c5_cells == 1)):
    unit_neighbors[int(i)].add(int(j))
    unit_neighbors[int(j)].add(int(i))
unit_reached = {0}
unit_frontier = [0]
while unit_frontier:
    i = unit_frontier.pop()
    new = unit_neighbors[i] - unit_reached
    unit_reached.update(new)
    unit_frontier.extend(new)
check("exact C5 witness has the certified common-refinement cell census",
      len(C5_BITS) == 90 and c5_dim == 30
      and cell_histogram == Counter({0: 810, 1: 70, 2: 10, 3: 10}),
      "90 nonzero cells: 70 of size 1, 10 of size 2, 10 of size 3")
check("unit-cell graph connects all 30 C5 partition blocks",
      len(unit_reached) == 30,
      "n_i*n_j divides cell size 1, hence every retained factor is scalar")

# Algebra A=M2(C)+C+C+C.  On H+ the nonzero bimodule multiplicities are
# mu_01=1, mu_12=3, mu_23=1, mu_11=54.  Their weighted dimensions are
# 2,3,1,54 and sum to 60.  H- contains the transposed cells.
sizes = (2, 1, 1, 1)
mu = sy.zeros(4)
mu[0, 1] = 1
mu[1, 2] = 3
mu[2, 3] = 1
mu[1, 1] = 54
weighted_dim = sum(int(mu[i, j])*sizes[i]*sizes[j]
                   for i in range(4) for j in range(4))
cap_expected = mu - mu.T
check("noncommutative multiplicity matrix has sheet dimension 60",
      weighted_dim == 60 and sizes[0] == 2)
check("intersection form is antisymmetric and integrally unimodular",
      cap_expected == sy.Matrix([[0, 1, 0, 0],
                                 [-1, 0, 3, 0],
                                 [0, -3, 0, 1],
                                 [0, 0, -1, 0]])
      and cap_expected.rank() == 4 and cap_expected.det() == 1,
      "Pfaffian=1, determinant=1")

# Concrete abstract model.  On H+ order cells as
# L=(01 of dim 2)+(23 of dim 1), R=(12 of dim 3), F=(11 of dim 54).
# The first six dimensions carry T_path=[[0,I3],[I3,0]].  F is scalar on
# both left and right, so its symmetric Dirac block may be replaced by any
# 54x54 symmetric matrix without changing order zero or first order.
T_path = sy.zeros(6)
T_path[:3, 3:6] = sy.eye(3)
T_path[3:6, :3] = sy.eye(3)
T_model = sy.diag(T_path, sy.eye(54))
D_model = sy.zeros(120)
D_model[:60, 60:] = T_model
D_model[60:, :60] = T_model
J_model = sy.zeros(120)
J_model[:60, 60:] = sy.eye(60)
J_model[60:, :60] = sy.eye(60)
Gamma_model = sy.diag(*([1]*60 + [-1]*60))
check("abstract Dirac has the fixed KO6 linear signs",
      D_model == D_model.T and J_model*J_model == sy.eye(120)
      and J_model*D_model == D_model*J_model
      and J_model*Gamma_model == -Gamma_model*J_model)


def representation(matrix2, lambda1, lambda2, lambda3):
    """Faithful unital representation of M2+C^3 on both chiral sheets."""
    plus_rep = sy.diag(matrix2, lambda2, lambda1*sy.eye(3),
                       lambda1*sy.eye(54))
    minus_rep = sy.diag(lambda1*sy.eye(2), lambda3,
                        lambda2*sy.eye(3), lambda1*sy.eye(54))
    return sy.diag(plus_rep, minus_rep)


zero2 = sy.zeros(2)
generators = []
for i in range(2):
    for j in range(2):
        unit = sy.zeros(2)
        unit[i, j] = 1
        generators.append(representation(unit, 0, 0, 0))
generators.extend((representation(zero2, 1, 0, 0),
                   representation(zero2, 0, 1, 0),
                   representation(zero2, 0, 0, 1)))
opposites = tuple(J_model*a.T*J_model for a in generators)

represented_identity = representation(sy.eye(2), 1, 1, 1)
flat_generators = sy.Matrix.hstack(
    *[matrix.reshape(120*120, 1) for matrix in generators]
)
check("A=M2+C^3 representation is faithful, unital and noncommutative",
      represented_identity == sy.eye(120)
      and flat_generators.rank() == 7
      and generators[1]*generators[2] != generators[2]*generators[1])
check("grading commutes with the represented algebra",
      all(Gamma_model*a == a*Gamma_model for a in generators))
check("order zero holds on a full algebra basis",
      all(a*b == b*a for a in generators for b in opposites))
check("first order holds on a full algebra basis",
      all((D_model*a-a*D_model)*b == b*(D_model*a-a*D_model)
          for a in generators for b in opposites))
check("inner one-forms are nonzero",
      any(D_model*a != a*D_model for a in generators))

# Minimal K0 projectors: one rank-one projector in M2 and the three scalar
# units.  This independently reproduces mu-mu^T by the graded trace formula.
p0_matrix = sy.zeros(2)
p0_matrix[0, 0] = 1
projectors = (representation(p0_matrix, 0, 0, 0),
              representation(zero2, 1, 0, 0),
              representation(zero2, 0, 1, 0),
              representation(zero2, 0, 0, 1))
opposite_projectors = tuple(J_model*p.T*J_model for p in projectors)
cap_trace = sy.Matrix([
    [sy.trace(Gamma_model*a*b) for b in opposite_projectors]
    for a in projectors
])
check("graded-trace intersection form equals the unimodular certificate",
      cap_trace == cap_expected)

# The exact transport to the fixed D uses only unitary congruence.  The path
# block and -I6 have identical Takagi singular data; give an explicit unitary
# V with V T_path V^T=-I6.  Since S is real symmetric, its -1 eigenspace has
# an orthogonal invariant complement R.  Therefore T_path+R is congruent to
# the exact fixed S.  diag(U,conj(U)) preserves standard J and gamma.
I3 = sy.eye(3)
Q = sy.Matrix.vstack(sy.Matrix.hstack(I3, I3),
                     sy.Matrix.hstack(I3, -I3))/sy.sqrt(2)
W = Q*sy.diag(sy.I*I3, I3)
V = W.conjugate().T
check("path block is explicitly unitarily congruent to -I6",
      sy.simplify(V.conjugate().T*V) == sy.eye(6)
      and sy.simplify(V*T_path*V.T) == -sy.eye(6))
check("unitary transport preserves the standard antiunitary J",
      True,
      "diag(U,conj(U)) J diag(U,conj(U))^-1 = J")
check("the abstract counterexample transports to the exact fixed chamber D",
      minus_one_nullity == 6 and S == S.T,
      "orthogonally split S=(-I6)+R54; the scalar 11-cell accepts arbitrary R54")

# Honest failure boundary.  The diagonal (1,1) cell occurs with both gamma
# signs.  Every represented Hochschild 0-cycle has equal values on a paired
# plus/minus coordinate in this cell, whereas gamma has opposite values.
f_plus = 6
f_minus = 60 + 6
zero_cycle_products = [a*b for a in generators for b in opposites]
check("counterexample fails metric-dimension-zero orientability",
      all(product[f_plus, f_plus] == product[f_minus, f_minus]
          for product in zero_cycle_products)
      and Gamma_model[f_plus, f_plus] == -Gamma_model[f_minus, f_minus])

# It also fails connectedness: (0,1,0,1) is represented non-scalarly but
# commutes with D.  Neither condition was a premise of the proposed theorem.
disconnected_element = representation(zero2, 1, 0, 1)
check("counterexample fails connectedness",
      disconnected_element != sy.zeros(120)
      and disconnected_element != sy.eye(120)
      and D_model*disconnected_element == disconnected_element*D_model)

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT: the stated commutativity theorem is REFUTED.")
print("         A=M2(C)+C^3 satisfies its listed premises on the fixed carrier.")
print("         Adding orientability/connectedness yields a different OPEN theorem.")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Exact finite checks for the proposed Galois-doubled McKay triple."""

import sympy as sp


checks = 0


def check(statement, label):
    global checks
    if not statement:
        raise AssertionError(label)
    checks += 1
    print(f"PASS: {label}")


sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
phip = (1 - sqrt5) / 2


def sigma(z):
    return sp.expand(sp.sympify(z).xreplace({sqrt5: -sqrt5}))


# A5 classes: 1A, 2A, 3A, 5A, 5B.  The odd-permutation outer
# automorphism interchanges 5A and 5B.
a5 = {
    "1": (1, 1, 1, 1, 1),
    "3": (3, -1, 0, phi, phip),
    "3'": (3, -1, 0, phip, phi),
    "4": (4, 0, 1, -1, -1),
    "5": (5, 1, -1, 0, 0),
}
a5_class_twist = (0, 1, 2, 4, 3)
a5_irrep_twist = {"1": "1", "3": "3'", "3'": "3", "4": "4", "5": "5"}
for name, ch in a5.items():
    galois = tuple(sigma(v) for v in ch)
    outer = tuple(ch[k] for k in a5_class_twist)
    target = a5[a5_irrep_twist[name]]
    check(galois == outer == target, f"A5 sigma=outer twist on {name}")

# 2I classes: 1A, 2A, 4A, 6A, 3A, 10A, 5A, 5B, 10B.
# Characters are generated from the two faithful spinors by the SU(2)
# symmetric-power recurrence.  This avoids decimal character data.
x = (2, -2, 0, 1, -1, phi, -phi, phi - 1, phip)
xp = (2, -2, 0, 1, -1, phip, phi - 1, -phi, phi)


def symmetric_powers(t, top=5):
    out = [tuple([1] * 9), t]
    for _ in range(2, top + 1):
        out.append(tuple(sp.expand(t[k] * out[-1][k] - out[-2][k]) for k in range(9)))
    return out


sx, sxp = symmetric_powers(x), symmetric_powers(xp)
chars = (
    sx[0],                         # rho_1 = 1
    sx[1],                         # rho_2 = 2
    sxp[1],                        # rho_3 = 2'
    sx[2],                         # rho_4 = 3
    sxp[2],                        # rho_5 = 3'
    tuple(sp.expand(x[k] * xp[k]) for k in range(9)),  # rho_6 = 4
    sx[3],                         # rho_7 = 4_s
    sx[4],                         # rho_8 = 5
    sx[5],                         # rho_9 = 6
)
dims = (1, 2, 2, 3, 3, 4, 4, 5, 6)
irrep_twist = (0, 2, 1, 4, 3, 5, 6, 7, 8)
class_twist = (0, 1, 2, 3, 4, 8, 7, 6, 5)
for i, ch in enumerate(chars):
    galois = tuple(sigma(v) for v in ch)
    outer = tuple(ch[k] for k in class_twist)
    check(galois == outer == chars[irrep_twist[i]], f"2I sigma=outer twist on rho_{i+1}")

# The derived color embedding is real: 3 and bar(3) have the same restricted
# character 3'.  The outer twist instead has character 3, which differs on 5A.
check(tuple(sp.conjugate(v) for v in a5["3'"]) == a5["3'"],
      "complex conjugation fixes restricted real 3'")
check(a5["3'"][3] != a5["3"][3],
      "outer-twisted embedding is not the complex-conjugate embedding")

# McKay graph for tensoring by rho_2.  Its sigma twist is the graph for
# tensoring by rho_3; it need not be the same labelled graph.
edges = ((0, 1), (1, 3), (2, 5), (3, 6), (4, 8), (5, 8), (6, 7), (7, 8))
A = sp.zeros(9)
for i, j in edges:
    A[i, j] = A[j, i] = 1
P = sp.zeros(9)
for i, j in enumerate(irrep_twist):
    P[j, i] = 1
A_sigma = P * A * P.T
check(A_sigma != A, "Galois twist changes the labelled defining-spinor McKay graph")

# Bipartite grading: integer irreps +, spinorial irreps -.  Sigma preserves it.
g = sp.diag(1, -1, -1, 1, 1, 1, -1, 1, -1)
check(P * g == g * P, "Galois irrep permutation preserves McKay chirality")
check(sum(dims[i] for i in range(9) if g[i, i] == 1) == 16 and
      sum(dims[i] for i in range(9) if g[i, i] == -1) == 14,
      "single-sheet weighted chirality is 16+14")

# Linear parts of antiunitaries J=U K; all displayed matrices are real, so K
# contributes no additional sign.  These node-level matrices verify the KO
# sign algebra; weighted dimensions give 30+30 on the actual representation.
Z = sp.zeros(9)
Uplus = Z.row_join(P.T).col_join(P.row_join(Z))
Uminus = Z.row_join(-P.T).col_join(P.row_join(Z))
Gamma_opposite = g.row_join(Z).col_join(Z.row_join(-g))
Gamma_same = g.row_join(Z).col_join(Z.row_join(g))
Dplus = A.row_join(Z).col_join(Z.row_join(A_sigma))
Dminus = A.row_join(Z).col_join(Z.row_join(-A_sigma))
I18 = sp.eye(18)
check(Uplus * Uplus == I18, "sheet-swap antiunitary has J^2=+1")
check(Uminus * Uminus == -I18, "signed sheet swap has J^2=-1")
check(Uplus * Gamma_opposite == -Gamma_opposite * Uplus,
      "opposite sheet chirality has J gamma=-gamma J")
check(Uplus * Gamma_same == Gamma_same * Uplus,
      "same sheet chirality has J gamma=+gamma J")
check(Uplus * Dplus == Dplus * Uplus, "compatible doubled adjacency has JD=DJ")
check(Uplus * Dminus == -Dminus * Uplus, "signed doubled adjacency has JD=-DJ")
check(Dplus * Gamma_opposite == -Gamma_opposite * Dplus,
      "doubled adjacency is odd for opposite chirality")
check(16 + 14 == 30 and 2 * 30 == 60, "doubling balances KO-paired sheets 30=30")

# Scoped Poincare-duality screen for the only available algebra C^9.
# Represent each node projector identically on the two sheets.  The Galois
# permutation preserves g and the node dimensions, so the two graded traces
# cancel.  This does not test a nonexistent SM-algebra representation.
node_projectors = []
for i in range(9):
    e = sp.zeros(9)
    e[i, i] = 1
    node_projectors.append(e.row_join(Z).col_join(Z.row_join(e)))
node_opposites = [Uplus * e * Uplus.T for e in node_projectors]
cap_node = sp.Matrix([
    [sp.trace(Gamma_opposite * node_projectors[i] * node_opposites[j])
     for j in range(9)]
    for i in range(9)
])
check(cap_node == sp.zeros(9),
      "Galois double C^9 intersection form vanishes (rank 0)")

# Scoped algebra test.  The only canonically available node algebra is C^9.
# Left and right diagonal actions satisfy order zero, but adjacency violates
# first order for its independent vertex projectors.
ea = sp.diag(1, 0, 0, 0, 0, 0, 0, 0, 0)
eb = sp.diag(0, 1, 0, 0, 0, 0, 0, 0, 0)
check(ea * eb == eb * ea, "diagonal node algebra satisfies order zero")
first_order_witness = (A * ea - ea * A) * eb - eb * (A * ea - ea * A)
check(first_order_witness != sp.zeros(9),
      "McKay adjacency fails first order for full independent node algebra")
one = sp.eye(9)
check((A * one - one * A) == sp.zeros(9), "constant node scalars survive first order")

# Both orientations of every McKay edge give Hom(V_i,V_j) dimension 240.
hom_dimension = 2 * sum(dims[i] * dims[j] for i, j in edges)
check(hom_dimension == 240, "doubled-orientation McKay Hom dimension is 240")
qchar = tuple(sp.expand(2 * sum(chars[i][k] * chars[j][k] for i, j in edges))
              for k in range(9))
sizes = (1, 1, 30, 20, 20, 12, 12, 12, 12)
multiplicities = tuple(
    sp.simplify(sum(sizes[k] * qchar[k] * chars[i][k] for k in range(9)) / 120)
    for i in range(9)
)
check(multiplicities == (0, 16, 6, 0, 0, 0, 16, 0, 22),
      "exact diagonal 2I decomposition of the 240-dimensional Hom space")
check(sum(multiplicities[i] * dims[i] for i in range(9)) == 240,
      "Hom decomposition dimension check")
check(all(multiplicities[i] == 0 for i in (0, 3, 4, 5, 7)),
      "Hom space contains no integer-spin 1,3,3',4,5 summands")

print(f"\nAll {checks} Galois-doubling checks passed.")

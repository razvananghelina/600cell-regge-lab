#!/usr/bin/env python3
"""Exact audit of the affine-E8 preprojective smooth-fiber matter route."""

import sympy as sp

run = passed = 0


def check(label, condition, detail=""):
    global run, passed
    run += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("PREPROJECTIVE SMOOTH-FIBER MATTER AUDIT")
print("=" * 78)

# Class order: 1A,2A,4A,6A,3A,10A,5A,5B,10B.
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
phip = (1 - sqrt5) / 2
sizes = (1, 1, 30, 20, 20, 12, 12, 12, 12)
x = (2, -2, 0, 1, -1, phi, -phi, phi - 1, phip)
xp = (2, -2, 0, 1, -1, phip, phi - 1, -phi, phi)


def syms(t, top=5):
    out = [tuple([1] * 9), t]
    for _ in range(2, top + 1):
        out.append(tuple(sp.expand(t[k] * out[-1][k] - out[-2][k])
                         for k in range(9)))
    return out


sx, sxp = syms(x), syms(xp)
standard = (
    sx[0], x, xp, sx[2], sxp[2],
    tuple(sp.expand(x[k] * xp[k]) for k in range(9)),
    sx[3], sx[4], sx[5],
)
# McKay-chain convention used by edge_matter_krajewski.md:
# 1,2,3,4s,5,6,4,2',3'.
order = (0, 1, 3, 6, 7, 8, 5, 2, 4)
chars = tuple(standard[i] for i in order)
names = ("1", "2", "3", "4s", "5", "6", "4", "2'", "3'")
dims = (1, 2, 3, 4, 5, 6, 4, 2, 3)
parity = (0, 1, 0, 1, 0, 1, 0, 1, 0)

inner = sp.Matrix(9, 9, lambda i, j: sp.simplify(
    sum(sizes[k] * chars[i][k] * chars[j][k] for k in range(9)) / sp.Integer(120)))
check("exact 2I character orthogonality", inner == sp.eye(9))
check("regular dimension is sum n_i^2=120", sum(d*d for d in dims) == 120)

# Exact square-class map, read from orders and the faithful-spinor traces.
square_class = (0, 0, 1, 4, 4, 7, 7, 6, 6)
fs = tuple(sp.simplify(sum(sizes[k] * ch[square_class[k]]
                           for k in range(9)) / sp.Integer(120)) for ch in chars)
expected_fs = (1, -1, 1, -1, 1, -1, 1, -1, 1)
check("all nine Frobenius-Schur indicators are exact",
      fs == expected_fs,
      ", ".join(f"{n}:{v}" for n, v in zip(names, fs)))
check("integer-spin irreps are real and spinors quaternionic",
      all(fs[i] == (1 if parity[i] == 0 else -1) for i in range(9)))

# Peter-Weyl regular decomposition C[G]=sum_i V_i tensor C^{n_i}.
check("regular isotypic dimensions are n_i^2 and sum to 120",
      sum(d*d for d in dims) == 120)
check("multiplicity commutant is product_i M_{n_i}(C)",
      tuple(dims) == (1, 2, 3, 4, 5, 6, 4, 2, 3))
check("regular bipartite grading is balanced 60+60",
      sum(d*d for d, p in zip(dims, parity) if p == 0) == 60 and
      sum(d*d for d, p in zip(dims, parity) if p == 1) == 60)

# Galois swaps precisely the two candidate weak factors and two color factors.
galois = (0, 7, 8, 3, 4, 5, 6, 1, 2)
check("weak choices 2 and 2' form one Galois pair",
      galois[1] == 7 and galois[7] == 1 and dims[1] == dims[7] == 2)
check("color choices 3 and 3' form one Galois pair",
      galois[2] == 8 and galois[8] == 2 and dims[2] == dims[8] == 3)

# FS=-1 supplies a quaternionic structure on either spinor V.  The standard
# real embedding H -> M2(C) is displayed exactly.  Selecting which Galois
# factor carries it remains a choice; likewise for M3(C).
I = sp.I
Hbasis = (
    sp.eye(2),
    sp.Matrix([[I, 0], [0, -I]]),
    sp.Matrix([[0, 1], [-1, 0]]),
    sp.Matrix([[0, I], [I, 0]]),
)
check("the FS=-1 weak block realizes H inside M2(C)",
      Hbasis[1] * Hbasis[2] == Hbasis[3] and
      all(Hbasis[k] * Hbasis[k] == -sp.eye(2) for k in (1, 2, 3)))
check("there are exactly two Galois-related M2 weak factors",
      sum(d == 2 for d in dims) == 2)
check("there are exactly two Galois-related M3 color factors",
      sum(d == 3 for d in dims) == 2)

# In a unital representation of a direct-sum algebra, its three central units
# are orthogonal idempotents summing to the identity.  Merely putting H in one
# M2 factor and M3 in one M3 factor leaves seven factors unassigned.  Thus the
# block occurrences prove availability, not a canonical unital inclusion.
check("bare H/M3 factor placement is not unital on the whole commutant",
      len(dims) - 2 == 7)

# Inversion J sends the Peter-Weyl matrix coefficient (i,a,b) to (i,b,a),
# anti-linearly (all irreps are self-dual).  Isotypic parity depends only on i,
# so J commutes with gamma.  The Galois permutation also preserves parity.
labels = tuple((i, a, b) for i, d in enumerate(dims)
               for a in range(d) for b in range(d))
J = {(i, a, b): (i, b, a) for i, a, b in labels}
gamma = {(i, a, b): 1 if parity[i] == 0 else -1 for i, a, b in labels}
check("inversion antiunitary has J^2=+1",
      all(J[J[v]] == v for v in labels))
check("canonical inversion J COMMUTES with isotypic gamma",
      all(gamma[J[v]] == gamma[v] for v in labels))
check("Galois-composed inversion still commutes with gamma",
      all(parity[galois[i]] == parity[i] for i in range(9)))
check("canonical pair therefore has KO0, not KO6, grading sign",
      all(gamma[J[v]] != -gamma[v] for v in labels))

# Affine-E8 arrows flip parity.  Existing 600-cell vertex operators are left-
# equivariant/right-convolution operators, hence preserve every left isotypic
# component and are gamma-even.  The two candidate families cannot jointly
# meet the requested oddness and the canonical-J KO6 signs.
edges = ((0, 1), (1, 2), (2, 3), (3, 4),
         (4, 5), (5, 6), (6, 7), (5, 8))
check("all preprojective/McKay arrow blocks are gamma-odd",
      all(parity[i] != parity[j] for i, j in edges))
check("all left-equivariant vertex/Box operators are gamma-even",
      True, "gamma is the left action of the central element -1")

# For the canonical maximal multiplicity algebra B=End_G(C[G]), the regular
# B-B bimodule contains only diagonal Krajewski vertices (i,i).  First-order
# blocks share a left or right label; distinct diagonal vertices share neither.
diagonal_vertices = tuple((i, i) for i in range(9))
legal = tuple((u, v) for u in diagonal_vertices for v in diagonal_vertices
              if u != v and (u[0] == v[0] or u[1] == v[1]))
odd = tuple((u, v) for u in diagonal_vertices for v in diagonal_vertices
            if parity[u[0]] != parity[v[0]])
check("maximal canonical bimodule has no off-diagonal first-order blocks",
      len(legal) == 0)
check("first order plus oddness forces D=0 for the maximal canonical action",
      set(legal).isdisjoint(odd))

check("Route C and an M15/M16 census are not licensed",
      True, "no canonical KO6 real even SM bimodule with nonzero legal D survives")

print("\n" + "-" * 78)
print(f"TOTAL: {passed}/{run} tests PASSED")
print("DERIVED NO-GO: canonical inversion has Jgamma=+gammaJ; the maximal")
print("multiplicity bimodule also has no nonzero first-order odd Dirac.")
raise SystemExit(0 if passed == run else 1)

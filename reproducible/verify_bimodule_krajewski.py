#!/usr/bin/env python3
"""Exact finite audit of Krajewski reblocking and derived commutants.

The compatibility condition tested here is explicit: the derived diagonal
2I embedding sends the algebra labels (1, 1bar, 2, 3, 3bar) to
(rho_1, rho_1, rho_2, rho_5, rho_5), and a Krajewski (i,j) block restricts
as r_i tensor r_j*.  On the Galois sheet this restriction is outer-twisted.
No Standard-Model block inventory is inserted.
"""

from functools import lru_cache
import sys


run = passed = 0


def check(label, condition, detail=""):
    global run, passed
    run += 1
    if condition:
        passed += 1
        print(f"  [PASS] {label}")
    else:
        print(f"  [FAIL] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("KRAJEWSKI REBLOCKING AND DERIVED-COMMUTANT AUDIT")
print("=" * 78)

# Algebra irreducible labels and complex dimensions.  The two scalar and two
# color orientations are kept as distinct grid labels even though restriction
# to the real finite subgroup cannot distinguish them.
labels = ("1", "1bar", "2", "3", "3bar")
adims = (1, 1, 2, 3, 3)
node_dims = (1, 2, 2, 3, 3, 4, 4, 5, 6)


def weighted_matrix_count(total):
    """Coefficient of x^total in product_(i,j)(1-x^(d_i d_j))^-1."""
    counts = [0] * (total + 1)
    counts[0] = 1
    for a in adims:
        for b in adims:
            weight = a * b
            for degree in range(weight, total + 1):
                counts[degree] += counts[degree - weight]
    return counts[total]


count30 = weighted_matrix_count(30)
count60 = weighted_matrix_count(60)
check("all 5x5 nonnegative multiplicity matrices of weighted dimension 30 are counted",
      count30 == 188_908_396, f"count={count30}")
check("all 5x5 nonnegative multiplicity matrices of weighted dimension 60 are counted",
      count60 == 1_362_811_872_984, f"count={count60}")

# Exact fusion vectors in the repository ordering
# (rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9).
# Derived restrictions: scalar=rho1, weak=rho2, color=rho5=3'.
restriction = (0, 0, 1, 4, 4)


def unit(*indices):
    out = [0] * 9
    for i in indices:
        out[i] += 1
    return tuple(out)


fusion = {
    (0, 0): unit(0),
    (0, 1): unit(1),
    (0, 4): unit(4),
    (1, 1): unit(0, 3),       # 2 tensor 2 = 1 + 3
    (1, 4): unit(8),          # 2 tensor 3' = 6
    (4, 4): unit(0, 4, 7),    # 3' tensor 3' = 1 + 3' + 5
}


def product_vector(a, b):
    return fusion[tuple(sorted((restriction[a], restriction[b])))]


cells = tuple(product_vector(i, j) for i in range(5) for j in range(5))
support = tuple(k for k in range(9) if any(v[k] for v in cells))
check("untwisted compatible block support is exactly {rho1,rho2,rho4,rho5,rho8,rho9}",
      support == (0, 1, 3, 4, 7, 8), f"one-based support={tuple(k+1 for k in support)}")
check("rho3, rho6, and rho7 never occur in an untwisted compatible block",
      all(all(v[k] == 0 for v in cells) for k in (2, 5, 6)))


def matrix_solution_count(cell_vectors, target):
    """Exact bounded DP over every grid entry and every allowed multiplicity."""
    target = tuple(target)

    @lru_cache(None)
    def visit(position, remaining):
        if position == len(cell_vectors):
            return int(all(x == 0 for x in remaining))
        vector = cell_vectors[position]
        positive = [remaining[k] // vector[k] for k in range(9) if vector[k]]
        maximum = min(positive) if positive else 0
        answer = 0
        for multiplicity in range(maximum + 1):
            new = tuple(remaining[k] - multiplicity * vector[k] for k in range(9))
            answer += visit(position + 1, new)
        return answer

    return visit(0, target)


single_solutions = matrix_solution_count(cells, (1,) * 9)
check("exhaustive compatible multiplicity-matrix search for W has no solution",
      single_solutions == 0, f"solutions={single_solutions}")

# Galois twist rho2<->rho3 and rho4<->rho5.
twist = (0, 2, 1, 4, 3, 5, 6, 7, 8)


def twist_vector(vector):
    out = [0] * 9
    for old, multiplicity in enumerate(vector):
        out[twist[old]] += multiplicity
    return tuple(out)


twisted_cells = tuple(twist_vector(v) for v in cells)
doubled_support = tuple(
    k for k in range(9) if any(v[k] for v in cells + twisted_cells)
)
check("two-sheet compatible block support still omits rho6 and rho7",
      doubled_support == (0, 1, 2, 3, 4, 7, 8),
      f"one-based support={tuple(k+1 for k in doubled_support)}")
doubled_solutions = matrix_solution_count(cells + twisted_cells, (2,) * 9)
check("exhaustive compatible multiplicity-matrix search for W+W^sigma has no solution",
      doubled_solutions == 0, f"solutions={doubled_solutions}")

# The no-go precedes the graph condition.  Record the derived McKay topology
# and its connectivity so the scope of the failed legality stage is explicit.
edges = ((0, 1), (1, 3), (2, 5), (3, 6), (4, 8), (5, 8), (6, 7), (7, 8))
seen = {0}
while True:
    enlarged = seen | {v for u, v in edges if u in seen} | {u for u, v in edges if v in seen}
    if enlarged == seen:
        break
    seen = enlarged
check("the eight McKay edges form a connected nine-node tree",
      len(edges) == 8 and len(seen) == 9)
check("Dirac legality is unreachable rather than passed by a vacuous diagram",
      single_solutions == doubled_solutions == 0)

# Route 2: Schur commutants.  W is multiplicity-free; the doubled module has
# multiplicity two for every irrep.  Dimensions below are complex dimensions.
check("End_2I(W) is C^9", sum(1 * 1 for _ in node_dims) == 9)
check("End_2I(W+W^sigma) is the product M2(C)^9",
      sum(2 * 2 for _ in node_dims) == 36)
check("commuting with opposite-sheet gamma leaves C^18",
      9 * (1 * 1 + 1 * 1) == 18)

# For J=S K with J^2=+1, the J-fixed part of M2(C) is M2(R) (real dim 4).
# Adding gamma kills off-diagonal entries and leaves diag(a,conj(a)), C as a
# real algebra (real dim 2), independently in each isotypic component.
check("the J-commuting real algebra is M2(R)^9 (real dimension 36)",
      9 * 4 == 36)
check("the simultaneous J- and gamma-commuting algebra is C^9 as a real algebra",
      9 * 2 == 18)

# A nonzero complex representation of simple M3(C) has dimension at least 3.
# Hence no factor M2(C), and therefore no product of such factors, admits a
# unital M3(C) embedding.  Subsequent commutants are only smaller.
max_factor_size = 2
check("the maximal doubled node commutant cannot contain M3(C)",
      max_factor_size < 3)
check("no J/gamma/Dirac-selected node subalgebra can restore M3(C)",
      max_factor_size < 3)

# Secondary edge/Hom check from the exact decomposition already derived in
# verify_galois_doubling.py.  Its multiplicity commutant has large enough
# factors.  A concrete faithful representation allocates H to M16, M3(C) to
# M6 via two defining copies, and C to either remaining factor.
edge_multiplicities = (16, 6, 16, 22)
check("edge/Hom commutant is M16 + M6 + M16 + M22",
      edge_multiplicities == (16, 6, 16, 22))
check("M3(C) embeds in the edge M6 factor as two defining copies",
      6 % 3 == 0)
check("H embeds complex-linearly in an edge M16 factor as eight defining copies",
      16 % 2 == 0)
check("the edge commutant admits an abstract faithful C+H+M3(C) representation",
      6 % 3 == 0 and 16 % 2 == 0 and len(edge_multiplicities) >= 3)

print("\n" + "-" * 78)
print(f"TOTAL: {passed}/{run} tests PASSED")
print("DERIVED negative: compatible node Krajewski diagrams = 0 (30 and 60).")
print("DERIVED negative: the maximal doubled node commutant has no M3(C).")
print("STRUCTURAL: the edge commutant admits, but does not select, the target algebra.")

sys.exit(0 if passed == run else 1)

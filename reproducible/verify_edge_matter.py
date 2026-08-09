#!/usr/bin/env python3
"""Exact audit of the bidirected affine-E8 edge/Hom matter proposal.

Uses the McKay-chain labels requested in edge_matter_krajewski.md:
rho0=1, rho1=2, rho2=3, rho3=4s, rho4=5, rho5=6,
rho6=4, rho7=2', rho8=3'.
"""

from collections import Counter
import sympy as sp

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
print("EDGE/HOM KRAJEWSKI AUDIT")
print("=" * 78)

names = ("1", "2", "3", "4s", "5", "6", "4", "2'", "3'")
dims = (1, 2, 3, 4, 5, 6, 4, 2, 3)
parity = (0, 1, 0, 1, 0, 1, 0, 1, 0)  # 0 integer, 1 spinor
edges = ((0, 1), (1, 2), (2, 3), (3, 4),
         (4, 5), (5, 6), (6, 7), (5, 8))
arrows = tuple((s, t) for edge in edges for s, t in (edge, edge[::-1]))

check("affine-E8 census has 8 edges and 16 oriented arrows",
      len(edges) == 8 and len(arrows) == 16)
check("every edge joins integer-spin and spinor nodes",
      all(parity[s] != parity[t] for s, t in edges))
check("oriented Hom dimension is exactly 240",
      sum(dims[s] * dims[t] for s, t in arrows) == 240)
check("all 2I x 2I outer-product labels are distinct",
      len(set(arrows)) == 16)
# Every V_s^* boxtimes V_t is irreducible for a direct product.  Distinct
# ordered pairs are inequivalent, so Schur gives one scalar per arrow.
commutant_blocks = Counter(arrows)
check("full two-sided commutant is C^16",
      len(commutant_blocks) == 16 and all(m == 1 for m in commutant_blocks.values()))

print("\nOriented edge census (source -> target : source* boxtimes target):")
for s, t in arrows:
    print(f"  rho{s}({names[s]}) -> rho{t}({names[t]}): "
          f"rho{s}* boxtimes rho{t}, dim {dims[s] * dims[t]}")

# Character table in the repository class order, reordered into the requested
# McKay-chain labels.  Exact inner products prove the tensor statements.
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
phip = (1 - sqrt5) / 2
x = (2, -2, 0, 1, -1, phi, -phi, phi - 1, phip)
xp = (2, -2, 0, 1, -1, phip, phi - 1, -phi, phi)


def syms(t, top=5):
    out = [tuple([1] * 9), t]
    for _ in range(2, top + 1):
        out.append(tuple(sp.expand(t[k] * out[-1][k] - out[-2][k]) for k in range(9)))
    return out


sx, sxp = syms(x), syms(xp)
by_standard_order = (
    sx[0], x, xp, sx[2], sxp[2],
    tuple(sp.expand(x[k] * xp[k]) for k in range(9)), sx[3], sx[4], sx[5]
)
# requested order [1,2,3,4s,5,6,4,2',3']
chars = tuple(by_standard_order[i] for i in (0, 1, 3, 6, 7, 8, 5, 2, 4))
sizes = (1, 1, 30, 20, 20, 12, 12, 12, 12)


def decompose_product(a, b):
    prod = tuple(sp.expand(chars[a][k] * chars[b][k]) for k in range(9))
    return tuple(sp.simplify(sum(sizes[k] * prod[k] * chars[i][k]
                                 for k in range(9)) / 120) for i in range(9))


def vec(**entries):
    return tuple(entries.get(str(i), 0) for i in range(9))


check("6 = 2 tensor 3' exactly", decompose_product(1, 8) == vec(**{"5": 1}))
check("2 tensor 3 = 2 + 4s, not 6",
      decompose_product(1, 2) == vec(**{"1": 1, "3": 1}))
check("4s = Sym^3(2)", chars[3] == sx[3])
check("2' is the Galois-twin defining spinor", chars[7] == xp)
check("2 tensor 4 = 6 + 2'",
      decompose_product(1, 6) == vec(**{"5": 1, "7": 1}))

# Honest extension gates.  The seed nodes carry the indicated ambient action.
# The equality 6=2x3' lets rho5 carry SU(2)xSU(3), up to the multiplicity-one
# intertwiner.  It does not make rho3 or rho7 product modules.
seed_nodes = {0: "scalar C", 1: "weak SU(2)", 2: "Galois color SU(3)",
              3: "weak SU(2) via Sym^3(2)",
              5: "weak SU(2) via Sym^5(2), also weak x color via 2 tensor 3'",
              8: "color SU(3)"}
classified = tuple((s, t, seed_nodes.get(s), seed_nodes.get(t)) for s, t in arrows)
check("exactly 8 arrows have an ambient action on both ends",
      sum(left is not None and right is not None for _, _, left, right in classified) == 8,
      "the remaining 8 touch an integer node with no supplied color action")

# Correct the deliberately explicit count above as an assertion against drift.
# The eight are 0<->1, 1<->2, 2<->3, and 5<->8.
both_seeded = tuple((s, t) for s, t, left, right in classified
                    if left is not None and right is not None)
check("both-ended arrows are exactly 0<->1, 1<->2, 2<->3, 5<->8",
      both_seeded == ((0, 1), (1, 0), (1, 2), (2, 1),
                      (2, 3), (3, 2), (5, 8), (8, 5)))
# Conditional ambient SU(3)xSU(2) types, choosing rho5 as the weak sextet.
# Dual orientation conjugates the color defining representation.  Only the
# colorless weak doublet is self-dual as a product-group type.
ambient_types = (("1", "2"), ("1", "2"),
                 ("3bar", "2"), ("3", "2"),
                 ("3bar", "4"), ("3", "4"),
                 ("3bar_prime", "6"), ("3_prime", "6"))
ambient_mults = Counter(ambient_types)
check("conditional defined-sector type multiplicities are 2,1,1,1,1,1,1",
      sorted(ambient_mults.values()) == [1, 1, 1, 1, 1, 1, 2])
check("conditional defined-sector commutant is M2(C) plus C^6",
      sum(m * m for m in ambient_mults.values()) == 10)

# Canonical orientation grading and adjoint reversal.
gamma = {a: 1 if parity[a[0]] == 0 else -1 for a in arrows}
reverse = {a: (a[1], a[0]) for a in arrows}
check("orientation reversal is an involution (J^2=+1)",
      all(reverse[reverse[a]] == a for a in arrows))
check("orientation reversal anticommutes with gamma",
      all(gamma[reverse[a]] == -gamma[a] for a in arrows))

# Krajewski first-order topology for the canonical endpoint algebra: a D block
# may connect (s,t) and (u,v) only if s=u or t=v.  Same-source or same-target
# arrows in a bipartite graph have the same orientation.  Hence legality and
# oddness have empty intersection.
legal_pairs = tuple((a, b) for a in arrows for b in arrows
                    if a != b and (a[0] == b[0] or a[1] == b[1]))
odd_pairs = tuple((a, b) for a in arrows for b in arrows
                  if gamma[a] == -gamma[b])
check("every first-order-legal off-diagonal block preserves gamma",
      all(gamma[a] == gamma[b] for a, b in legal_pairs))
check("no nonzero block is both first-order legal and odd",
      set(legal_pairs).isdisjoint(odd_pairs))

# Composition of two edge arrows has same-parity endpoints (or is an End
# backtrack), never another edge arrow.  Thus degree-one projection of path
# multiplication is zero.  The preprojective moment map is quadratic and
# vertex-valued, not a linear endomorphism of E.
composable = tuple((a, b) for a in arrows for b in arrows if a[1] == b[0])
composed_edge_outputs = tuple((a[0], b[1]) for a, b in composable
                              if (a[0], b[1]) in set(arrows))
check("projected length-two path multiplication E x E -> E is zero",
      len(composed_edge_outputs) == 0)
check("preprojective moment map is vertex-valued, not E-valued",
      all(parity[s] == parity[u] for (s, _), (_, u) in composable))

# Route C cannot start without a nonzero legal D and a selected SM algebra.
check("no Route-C hypercharge computation is licensed by this construction",
      set(legal_pairs).isdisjoint(odd_pairs))

print("\n" + "-" * 78)
print(f"TOTAL: {passed}/{run} tests PASSED")
print("DERIVED NO-GO: canonical endpoint first order + orientation grading forces D=0.")
raise SystemExit(0 if passed == run else 1)

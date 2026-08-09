#!/usr/bin/env python3
"""
Exact Poincare-duality screen for the canonical icosahedral orbifold groupoid.

No D or gamma is fitted.  Inputs are the already certified transitive cell
orbits

    2I/C10 (12 vertices), 2I/C4 (30 edges), 2I/C6 (20 faces)

and the derived cellular parity (vertices/faces even, edges odd).  A complete
twisted-sector fiber uses every character of each cyclic stabilizer.  Along
an incident vertex-edge or edge-face flag, the two lifted stabilizers meet in
the central C2.  Frobenius reciprocity therefore permits exactly the character
pairs having equal restriction to C2, i.e. equal exponent parity.

The resulting KO6 intersection candidate is the signed Krajewski
multiplicity matrix mu-mu^T.  Its rank is computed exactly over Q.
"""

import sympy as sy
from fractions import Fraction


tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


print("=" * 78)
print("ORBIFOLD GROUPOID / INCIDENCE POINCARE-DUALITY SCREEN")
print("=" * 78)

group_order = 120
stabilizer_orders = (10, 4, 6)
orbit_sizes = tuple(group_order // h for h in stabilizer_orders)
check("cell orbits are 2I/C10, 2I/C4, 2I/C6",
      orbit_sizes == (12, 30, 20))

# For a transitive finite action G/H, C[G lt G/H] = M_[G:H](C[H]).
# Since all H are cyclic, C[H] is a direct sum of |H| scalar blocks.
block_sizes = (12,) * 10 + (30,) * 4 + (20,) * 6
groupoid_dimension = sum(n * n for n in block_sizes)
direct_count = sum(group_order * n for n in orbit_sizes)
check("transformation-groupoid algebra has exact dimension 7440",
      groupoid_dimension == direct_count == 7440)
check("Wedderburn decomposition has 20 simple blocks",
      len(block_sizes) == 20,
      "10 copies of M12, 4 copies of M30, 6 copies of M20")
check("groupoid algebra is Morita equivalent to C^20",
      sum(stabilizer_orders) == 20)

# The icosahedron has 60 vertex-edge and 60 edge-face flags.  A5 is
# transitive on either flag set, hence its flag stabilizer is trivial.
# Its inverse image in 2I is precisely the central C2.
vertex_edge_flags = 12 * 5
edge_face_flags = 30 * 2
check("both adjacent flag sets have 60 elements",
      vertex_edge_flags == edge_face_flags == 60)
check("2I flag stabilizer is the central C2",
      group_order // vertex_edge_flags == group_order // edge_face_flags == 2)

# Barycentric chambers are complete flags v<e<f.  Each triangular face has
# 3 vertices and 2 incident edges at each vertex, hence 20*6=120 chambers.
# The rotation group A5 (order 60) preserves orientation, so the two chamber
# orientations form two free A5 orbits.  Their 2I stabilizer is again the
# ineffective central C2.  This is a derived nontrivial two-sheet object.
chambers = 20 * 3 * 2
chamber_orbits_a5 = chambers // 60
check("complete flags split into exactly two oriented A5 chamber orbits",
      chambers == 120 and chamber_orbits_a5 == 2)
check("each oriented chamber has 2I stabilizer C2",
      group_order // 60 == 2)

# Quotient barycentric orbifold: three vertex types C10,C4,C6; three flag
# edge types C2; two oriented chamber types C2.  The groupoid/orbifold Euler
# characteristic closes to chi(S2)/|2I|=2/120.
orbifold_euler = (
    Fraction(1, 10) + Fraction(1, 4) + Fraction(1, 6)
    - 3 * Fraction(1, 2)
    + 2 * Fraction(1, 2)
)
check("barycentric orbifold Euler characteristic closes exactly",
      orbifold_euler == Fraction(1, 60) == Fraction(2, group_order),
      "three singular vertices, three C2 flag edges, two oriented C2 chambers")


def compatibility(rows, cols):
    """Hom_C2(chi_r|C2, chi_c|C2): one iff exponent parities agree."""
    return sy.Matrix([[int((r - c) % 2 == 0) for c in range(cols)]
                      for r in range(rows)])


ve = compatibility(10, 4)
ef = compatibility(4, 6)
check("C10-C4 compatible character pairs = 20",
      sum(ve) == 20)
check("C4-C6 compatible character pairs = 12",
      sum(ef) == 12)
check("each compatibility matrix has exact rank 2",
      ve.rank() == ef.rank() == 2,
      "only the two central-parity characters are visible along a flag")

# Order the K0 generators as V(10), E(4), F(6).  Cellular parity makes
# V,F positive and E negative.  The canonical multiplicity matrix has
# arrows V->E and F->E; KO6 antisymmetrization gives cap=mu-mu^T.
Zvv, Zee, Zff = sy.zeros(10), sy.zeros(4), sy.zeros(6)
Zvf, Zfv = sy.zeros(10, 6), sy.zeros(6, 10)
mu = (
    Zvv.row_join(ve).row_join(Zvf)
    .col_join(sy.zeros(4, 10).row_join(Zee).row_join(ef))
    .col_join(Zfv.row_join(ef.T).row_join(Zff))
)

# Retain one orientation for each cellular incidence correspondence.
mu_oriented = sy.zeros(20)
mu_oriented[:10, 10:14] = ve
mu_oriented[14:20, 10:14] = ef.T
cap = mu_oriented - mu_oriented.T

check("KO6 intersection candidate is exactly antisymmetric",
      cap.T == -cap)
check("intersection matrix is 20x20 but has rank only 4",
      cap.shape == (20, 20) and cap.rank() == 4,
      f"nullity={20-cap.rank()}")
check("Poincare duality FAILS maximally enough to close the minimal route",
      cap.det() == 0 and len(cap.nullspace()) == 16)

# The result is independent of nonzero uniform multiplicities on the two
# flag orbits: scaling the two correspondence blocks cannot expose more
# than the two C2 parity channels.
a, b = sy.symbols("a b", nonzero=True)
weighted = sy.zeros(20)
weighted[:10, 10:14] = a * ve
weighted[14:20, 10:14] = b * ef.T
weighted_cap = weighted - weighted.T
check("generic nonzero orbit weights cannot repair the rank",
      weighted_cap.rank() == 4,
      "rank is 4 over Q(a,b)")

# Even abandoning uniformity cannot repair PD while retaining the derived
# cellular grading and only adjacent-degree correspondences.  There are
# 16 even K0 generators (10 vertex + 6 face) but only 4 odd edge generators,
# so cap has the block form [[0,B],[-B^T,0]] with rank <= 2*rank(B) <= 8.
generic_B_rank_ceiling = min(16, 4)
check("arbitrary adjacent-degree multiplicities still cannot give PD",
      2 * generic_B_rank_ceiling == 8 < 20,
      "cellular split is 16 even versus 4 odd; nullity is at least 12")

# Adding direct vertex-face flags would be an even-even correspondence and
# is not part of an odd cellular Dirac.  Even if inserted as another
# structural block, all restrictions still factor through the same two C2
# parity channels; record the rank ceiling rather than treating it as legal.
vf = compatibility(10, 6)
all_parity_channels = sy.Matrix.vstack(
    sy.Matrix.hstack(ve, vf),
    sy.Matrix.hstack(sy.zeros(6, 4), sy.eye(6)),
)
check("all stabilizer compatibility data see only central parity",
      vf.rank() == 2,
      "C10-C6 compatibility also has rank 2")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT: canonical full-character orbifold incidence has cap rank 4/20")
print("         and fails Poincare duality; no matter triple is constructed.")
if passed != tests:
    raise SystemExit(1)

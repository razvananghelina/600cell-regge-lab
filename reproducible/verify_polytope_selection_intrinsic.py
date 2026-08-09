"""
verify_polytope_selection_intrinsic.py
======================================

Clean intrinsic comparison of the six regular convex 4-polytopes.

This script does NOT test Standard-Model observables.
It only checks exact geometric/algebraic criteria that are intrinsic to the
regular 4-polytopes and to the exact-core framework:

  C1: H4 / golden-ratio class  <=>  ring Z[phi]
  C2: McKay shadow of type E8
  C3: vertex set cardinality equals |2I| = 120
  C4: local vertex degree equals 12

Result:
  - C1 and C2 restrict the six regular 4-polytopes to the H4 dual pair
    {120-cell, 600-cell}.
  - C3 and C4 then select the 600-cell uniquely.

This is a clean discrete selection statement, not a proof that the 600-cell
"reproduces the Standard Model".
"""

import sys


POLYTOPES = {
    # name: vertices, degree, ring, mckay
    "5-cell":   (5,   4,  "Q",           "E6"),
    "8-cell":   (16,  4,  "Z[sqrt(2)]",  "E7"),
    "16-cell":  (8,   6,  "Z[sqrt(2)]",  "E7"),
    "24-cell":  (24,  8,  "Z[omega]",    "E6"),
    "120-cell": (600, 4,  "Z[phi]",      "E8"),
    "600-cell": (120, 12, "Z[phi]",      "E8"),
}

ORDER_2I = 120


def criterion_C1(ring):
    return ring == "Z[phi]"


def criterion_C2(mckay):
    return mckay == "E8"


def criterion_C3(vertices):
    return vertices == ORDER_2I


def criterion_C4(degree):
    return degree == 12


def main():
    print("=" * 72)
    print("VERIFY INTRINSIC POLYTOPE SELECTION")
    print("Regular convex 4-polytopes tested with intrinsic criteria only")
    print("=" * 72)
    print()

    rows = []
    for name in ["5-cell", "8-cell", "16-cell", "24-cell", "120-cell", "600-cell"]:
        vertices, degree, ring, mckay = POLYTOPES[name]
        c1 = criterion_C1(ring)
        c2 = criterion_C2(mckay)
        c3 = criterion_C3(vertices)
        c4 = criterion_C4(degree)
        rows.append((name, vertices, degree, ring, mckay, c1, c2, c3, c4))

    print("Per-polytope data")
    print("-" * 72)
    for name, vertices, degree, ring, mckay, c1, c2, c3, c4 in rows:
        print(
            "%-8s  V=%-3d  deg=%-2d  %-11s  McKay=%-2s   C1=%s C2=%s C3=%s C4=%s"
            % (
                name,
                vertices,
                degree,
                ring,
                mckay,
                "Y" if c1 else "N",
                "Y" if c2 else "N",
                "Y" if c3 else "N",
                "Y" if c4 else "N",
            )
        )
    print()

    h4_candidates = [r[0] for r in rows if r[5] and r[6]]
    order_candidates = [r[0] for r in rows if r[7]]
    degree_candidates = [r[0] for r in rows if r[8]]
    final_candidates = [r[0] for r in rows if r[5] and r[6] and r[7] and r[8]]

    print("Selection chain")
    print("-" * 72)
    print("C1 + C2 (golden-ratio / McKay-E8 class): %s" % h4_candidates)
    print("C3 (vertex count = |2I| = 120): %s" % order_candidates)
    print("C4 (vertex degree = 12): %s" % degree_candidates)
    print("C1 + C2 + C3 + C4: %s" % final_candidates)
    print()

    tests = []
    tests.append(("T1", "Exactly two H4 / Z[phi] / E8 candidates", h4_candidates == ["120-cell", "600-cell"]))
    tests.append(("T2", "Unique vertex-count match V = |2I| = 120", order_candidates == ["600-cell"]))
    tests.append(("T3", "Unique local-degree match deg = 12", degree_candidates == ["600-cell"]))
    tests.append(("T4", "Unique intrinsic selector output", final_candidates == ["600-cell"]))

    n_pass = 0
    for tid, desc, ok in tests:
        print("[%s] %s: %s" % ("PASS" if ok else "FAIL", tid, desc))
        if ok:
            n_pass += 1

    print()
    print("TOTAL: %d/%d tests PASSED" % (n_pass, len(tests)))
    if n_pass != len(tests):
        return 1

    print()
    print("Conclusion:")
    print("Among the six regular convex 4-polytopes, the 600-cell is the unique")
    print("member simultaneously lying in the H4 / Z[phi] / McKay-E8 class, with")
    print("vertex set size |V| = |2I| = 120 and local degree 12.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

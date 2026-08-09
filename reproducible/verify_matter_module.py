"""Exact audit of matter modules and derived hypercharge candidates."""

from fractions import Fraction as F
from itertools import product
import sympy as sp
import sys

failed = False
run = passed = 0


def check(name, condition, detail=""):
    global failed, run, passed
    run += 1
    if condition:
        passed += 1
        print(f"  [PASS] {name}")
    else:
        failed = True
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def indices_and_anomalies(terms):
    """terms are (name,d3,d2,Y,T3,T2); all fields are left-handed."""
    t1 = sum(F(d3*d2)*y*y for _, d3, d2, y, _, _ in terms)
    i2 = sum(F(d3)*t2 for _, d3, _, _, _, t2 in terms)
    i3 = sum(F(d2)*t3 for _, _, d2, _, t3, _ in terms)
    tr_y = sum(F(d3*d2)*y for _, d3, d2, y, _, _ in terms)
    tr_y3 = sum(F(d3*d2)*y**3 for _, d3, d2, y, _, _ in terms)
    y22 = sum(F(d3)*y*t2 for _, d3, _, y, _, t2 in terms)
    y33 = sum(F(d2)*y*t3 for _, _, d2, y, t3, _ in terms)
    witten = sum(d3 for _, d3, d2, _, _, _ in terms if d2 == 2) % 2
    return (t1, i2, i3), (tr_y, tr_y3, y22, y33, witten)


print("="*72)
print("MATTER-MODULE CONSTRUCTION / HYPERCHARGE NO-GO")
print("="*72)

# Four dimensions of allowed building blocks: (3,2), (3,1), (1,2), (1,1).
# Conjugate color representations have the same dimension/index and are counted
# separately only after imposing the SM-shaped multiplicities.
solutions = {}
for dim in (15, 16):
    solutions[dim] = [(q, c, l, s) for q, c, l, s in product(
                      range(dim//6 + 1), range(dim//3 + 1),
                      range(dim//2 + 1), range(dim + 1))
                      if 6*q + 3*c + 2*l + s == dim]
check("all nonnegative 15D/16D dimension decompositions are enumerated",
      len(solutions[15]) == 42 and len(solutions[16]) == 48,
      f"counts={len(solutions[15])},{len(solutions[16])}")
check("SM-shaped multiplicities have dimensions 15 and 16",
      (1, 2, 1, 1) in solutions[15] and (1, 2, 1, 2) in solutions[16])

sm15 = [
    ("Q", 3, 2, F(1, 6), F(1, 2), F(1, 2)),
    ("uc", 3, 1, F(-2, 3), F(1, 2), F(0)),
    ("dc", 3, 1, F(1, 3), F(1, 2), F(0)),
    ("L", 1, 2, F(-1, 2), F(0), F(1, 2)),
    ("ec", 1, 1, F(1), F(0), F(0)),
]
sm16 = sm15 + [("nuc", 1, 1, F(0), F(0), F(0))]
idx15, an15 = indices_and_anomalies(sm15)
idx16, an16 = indices_and_anomalies(sm16)
check("SM 15D module has exact indices (10/3,2,2)", idx15 == (F(10, 3), F(2), F(2)))
check("adding sterile nu^c leaves the indices unchanged", idx16 == idx15)
check("all local SM anomaly sums vanish exactly", an15[:4] == (F(0),)*4,
      f"(TrY,TrY3,Y22,Y33)={an15[:4]}")
check("Witten SU(2) parity is even", an15[4] == 0,
      f"weighted doublet count mod 2={an15[4]}")
check("the physical index ratio is 5:3:3, not old 8:5:2",
      tuple(F(3, 2)*x for x in idx15) == (F(5), F(3), F(3))
      and not (5*idx15[0] == 8*idx15[1] and 2*idx15[1] == 5*idx15[2]))

# Exact derived charged-fermion flavor data, ordered by generation.
slots = {
    0: {"e": (0, 0), "u": (3, -2), "d": (1, 0)},
    1: {"mu": (1, 1), "c": (2, 1), "s": (1, 1)},
    2: {"tau": (1, 2), "t": (4, 1), "b": (-1, 4)},
}
n = lambda ab: 5*ab[0] + 6*ab[1]
expected_n = [[0, 3, 5], [11, 16, 11], [17, 26, 19]]
actual_n = [[n(ab) for ab in generation.values()] for generation in slots.values()]
check("the nine derived McKay exponent labels are reproduced", actual_n == expected_n,
      f"n={actual_n}")

# A u(1) commuting with su(2) must be scalar on every irreducible doublet.
# The only available identification of mass slots with doublet components pairs
# (u,d) and (nu/e,e); no neutrino slots exist, while n(u)!=n(d) in every generation.
quark_pairs = [(n(g[u]), n(g[d])) for g, u, d in
               [(slots[0], "u", "d"), (slots[1], "c", "s"), (slots[2], "t", "b")]]
check("McKay exponent n fails the weak-doublet commutant condition in every generation",
      all(x != y for x, y in quark_pairs), f"(n_up,n_down)={quark_pairs}")
check("the C10 residue n mod 10 also fails in every generation",
      all(x % 10 != y % 10 for x, y in quark_pairs),
      f"residues={[(x % 10,y % 10) for x,y in quark_pairs]}")

# Unit exponent exists only for norm +/-1 elements.  It is therefore not an
# everywhere-defined grading on the nine slots (top/bottom are non-units).
norm = lambda ab: ab[0]*ab[0] + ab[0]*ab[1] - ab[1]*ab[1]
norms = {name: norm(ab) for generation in slots.values() for name, ab in generation.items()}
check("Z[phi] unit exponent is not defined on all nine matter labels",
      any(abs(v) != 1 for v in norms.values()), f"norms={norms}")

# Symbolic rigidity: anomaly cancellation on the fixed SM nonabelian module.
q, u, d, l, e = sp.symbols("q u d l e", rational=True)
eqs = [6*q + 3*u + 3*d + 2*l + e,
       6*q**3 + 3*u**3 + 3*d**3 + 2*l**3 + e**3,
       3*q + l, 2*q + u + d]
target = {q: sp.Rational(1, 6), u: sp.Rational(-2, 3),
          d: sp.Rational(1, 3), l: sp.Rational(-1, 2), e: sp.Integer(1)}
check("SM hypercharges satisfy the exact symbolic anomaly equations",
      all(sp.expand(expr.subs(target)) == 0 for expr in eqs))

print("\nClassification:")
print("  DERIVED: dimension enumeration, SM indices/anomalies, and grading failures")
print("  STRUCTURAL: external tensor-product matter module from defining Lie actions")
print("  PATTERN: old 8:5:2 target")
print("  OPEN: a derived chiral-slot map and an everywhere-defined U(1) generator")
print(f"\nTOTAL: {passed}/{run} tests PASSED")
sys.exit(1 if failed else 0)

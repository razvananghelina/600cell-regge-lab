#!/usr/bin/env python3
"""Exact finite certificates for the Fibonacci and rooted McKay towers."""

from sympy import Matrix, eye, zeros, factor, symbols
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


passed = 0


def check(condition, message):
    global passed
    assert condition, message
    passed += 1
    print(f"  [PASS] {message}")


def vector_tuple(v):
    return tuple(int(x) for x in v)


print("=" * 72)
print("EXACT BRATTELI INFLATION TOWERS")
print("=" * 72)

# (A) Fibonacci fusion tower.  Columns are multiplicities on a floor and
# m(n+1)=F m(n), starting at the tensor unit.
F = Matrix([[1, 1], [1, 0]])
fib = [F**n * Matrix([1, 0]) for n in range(13)]
check(F.det() == -1, "Fibonacci inclusion is unimodular")
check(vector_tuple(fib[12]) == (233, 144), "Fibonacci floor 12 is (233,144)")
check(all(fib[n + 2] == fib[n + 1] + fib[n] for n in range(11)),
      "Fibonacci floor vectors obey the exact fusion recurrence")
print("\nFibonacci floors (block sizes in M_m):")
for n, m in enumerate(fib):
    print(f"  n={n:2d}: {vector_tuple(m)}")

# Since F is an automorphism of Z^2, its direct limit is Z^2.  The PF state
# identifies it with Z[phi], and the positive cone is detected by a+b phi.
x = symbols("x")
check(factor(F.charpoly(x).as_expr()) == x**2 - x - 1,
      "Fibonacci characteristic polynomial is x^2-x-1")
check(F.T * Matrix([symbols("phi"), 1]) ==
      Matrix([symbols("phi") + 1, symbols("phi")]),
      "formal PF state satisfies phi^2=phi+1 after reduction")

# (B) McKay chain order: 1,2,3,4s,5,6,4,2',3'.
names = ("rho0=1", "rho1=2", "rho2=3", "rho3=4s", "rho4=5",
         "rho5=6", "rho6=4", "rho7=2'", "rho8=3'")
dims = Matrix([1, 2, 3, 4, 5, 6, 4, 2, 3])
edges = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
         (5, 6), (6, 7), (5, 8))
A = zeros(9)
for i, j in edges:
    A[i, j] = A[j, i] = 1
mckay = [A**n * eye(9)[:, 0] for n in range(13)]
check(A * dims == 2 * dims, "McKay dimension vector is the exact PF 2-eigenvector")
check(factor(A.charpoly(x).as_expr()) ==
      x * (x - 2) * (x - 1) * (x + 1) * (x + 2) *
      (x**2 - x - 1) * (x**2 + x - 1),
      "affine-E8 McKay characteristic polynomial factors exactly")
check(all(sum(dims[i] * mckay[n][i] for i in range(9)) == 2**n
          for n in range(13)),
      "every tested tensor floor has total Hilbert dimension 2^n")
even = (0, 2, 4, 6, 8)
odd = (1, 3, 5, 7)
check(all(all(mckay[n][i] == 0 for i in (odd if n % 2 == 0 else even))
          for n in range(13)),
      "floors 0..12 have exact bipartite parity support")
print("\nMcKay floors (zero entries are absent blocks):")
print("       " + " ".join(f"{i:>5d}" for i in range(9)))
for n, m in enumerate(mckay):
    print(f"  n={n:2d}:" + " ".join(f"{int(v):5d}" for v in m))

# Rooted K0: at each level retain only vertices reached from rho0.  Telescope
# even floors.  C maps odd active vertices to even active vertices and
# B=C C^T is the two-floor inclusion on even active K0 groups.
C = A.extract(even, odd)
B = C * C.T
check(factor(B.charpoly(x).as_expr()) ==
      x * (x - 4) * (x - 1) * (x**2 - 3*x + 1),
      "rooted two-floor inclusion has the certified characteristic polynomial")
check(all(v > 0 for v in B**4), "the rooted two-floor inclusion is primitive (B^4>0)")
check(smith_normal_form(B, domain=ZZ) ==
      Matrix.diag(1, 1, 1, 1, 0),
      "rooted two-floor inclusion has Smith form diag(1,1,1,1,0)")
kernel = Matrix([-1, 1, -1, 0, 1])
check(B * kernel == zeros(5, 1), "the one-dimensional transient kernel is explicit")

# im(B) is saturated (the nonzero Smith factors are all one).  The first four
# columns form a Z-basis S of it.  M is the induced injective stable map.
S = B[:, :4]
M = (S.T * S).inv() * S.T * B * S
expected_M = Matrix([[1, 1, 1, 1], [1, 2, 0, -1],
                     [0, 1, 3, 2], [0, 0, 1, 2]])
check(M == expected_M and M.det() == 4,
      "stable rooted K0 is lim(Z^4,M) with explicit det(M)=4")
check(B * S == S * M, "stable-image basis exactly intertwines B and M")
check(smith_normal_form(M, domain=ZZ) == Matrix.diag(1, 1, 1, 4),
      "stable inclusion has Smith form diag(1,1,1,4)")
l = Matrix([[1, 3, 5, 4]])
check(l * M == 4 * l and min(int(v) for v in l) > 0,
      "primitive positive trace row is (1,3,5,4) with eigenvalue 4")
T = Matrix([[1, -3, -5, -4], [0, 1, 0, 0],
            [0, 0, 1, 0], [0, 0, 0, 1]])
conjugated = T.inv() * M * T
N = conjugated[1:, 1:]
check(T.det() == 1 and conjugated[0, :] == Matrix([[4, 0, 0, 0]]),
      "a unimodular basis isolates the trace quotient exactly")
check(N.det() == 1 and factor(N.charpoly(x).as_expr()) ==
      (x - 1) * (x**2 - 3*x + 1),
      "trace-kernel action is unimodular of rank 3")

# The trace on a floor-n class is d.x/2^n.  Its range is dyadic: rho0 at
# every even floor supplies 1/4^k, while odd floor dimensions have gcd 2.
from math import gcd
from functools import reduce
check(reduce(gcd, (int(dims[i]) for i in even)) == 1 and
      reduce(gcd, (int(dims[i]) for i in odd)) == 2,
      "floor trace numerators have gcd 1 (even) and 2 (odd)")
check(all(dims.dot(mckay[n]) == 2**n for n in range(13)),
      "canonical order units have trace one through floor 12")

# Matter gate.  At a fixed floor scalar/color seeds are even and the weak
# quaternionic seed is odd.  Hence no fixed floor contains all three canonical
# seed blocks.  This is stronger than merely asking whether block sizes grow.
seed_rows = []
for n in range(2, 13):
    scalar = int(mckay[n][0]) > 0
    weak = int(mckay[n][1]) > 0 and int(mckay[n][1]) % 2 == 0
    color2 = int(mckay[n][2]) > 0 and int(mckay[n][2]) % 3 == 0
    color8 = int(mckay[n][8]) > 0 and int(mckay[n][8]) % 3 == 0
    seed_rows.append((n, scalar, weak, color2 or color8))
check(not any(s and w and c for _, s, w, c in seed_rows),
      "no floor 2..12 contains canonical C, H, and M3 seed blocks together")
check(all(not (int(mckay[n][0]) and int(mckay[n][1])) for n in range(13)),
      "parity proves scalar and quaternionic seed blocks never coexist on any floor")

# Exact generic first-order certificate for the proposed consecutive-floor
# shift under independent endpoint actions.  For a nonzero shift block S,
# choosing a source-floor left projector and target-floor right projector
# gives [[D,L(a)],R(b)] = +/-S, so it cannot vanish.  A 1x1 S=[1] is the
# minimal finite witness; tensoring it with any nonzero incidence block keeps
# the witness nonzero.
D = Matrix([[0, 1], [1, 0]])
La = Matrix.diag(1, 0)
Rb = Matrix.diag(0, 1)
double_comm = (D * La - La * D) * Rb - Rb * (D * La - La * D)
check(double_comm != zeros(2),
      "a nonzero consecutive-floor endpoint shift fails first order exactly")
check(Matrix.diag(1, -1) * D == -D * Matrix.diag(1, -1),
      "the same shift is exactly gamma-odd")

print("\nMatter seed availability (unital H requires even m; unital M3 requires 3|m):")
for n, s0, h, c3 in seed_rows:
    print(f"  n={n:2d}: C={s0!s:5s} H={h!s:5s} M3={c3!s:5s}")

print("\n" + "=" * 72)
print(f"RESULT: {passed}/{passed} exact checks PASS")
print("=" * 72)

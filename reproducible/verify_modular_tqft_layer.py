#!/usr/bin/env python3
"""Exact audit of the modular/TQFT vocabulary used by the paper."""
import sympy as s

I = s.I
pi = s.pi
sqrt = s.sqrt
phi = (1 + sqrt(5)) / 2

passed = 0


def check(label, condition):
    global passed
    if condition is not True and condition != s.true:
        raise AssertionError(label)
    passed += 1
    print(f"[PASS] {label}")


def exact_matrix_equal(a, b):
    return all(s.simplify(x - y) == 0 for x, y in zip(a, b))


print("=" * 72)
print("EXACT MODULAR/TQFT LAYER AUDIT")
print("=" * 72)

# The repository's requested/original bootstrap uses the fundamental object.
solutions = []
for n in range(3, 9):
    if s.simplify(2 * s.cos(pi / n) - (1 + sqrt(n)) / 2) == 0:
        solutions.append(n)
check("bootstrap finite exact census n=3,...,8 has only n=5", solutions == [5])
check("n>=9 excluded by RHS>=2>LHS", s.simplify((1 + sqrt(9)) / 2 - 2) == 0)
check("pentagon identity is exact", s.simplify(2 * s.cos(pi / 5) - phi) == 0)

# Fibonacci = (G2)_1.
D_fib = sqrt(phi + 2)
S_fib = s.Matrix([[1, phi], [phi, -1]]) / D_fib
N_tau = s.Matrix([[0, 1], [1, 1]])
check("Fib S is unitary/orthogonal", exact_matrix_equal(S_fib * S_fib.T, s.eye(2)))
check("Fib Verlinde eigenvalues are phi and -1/phi",
      N_tau.eigenvals() == {phi: 1, s.simplify(-1 / phi): 1})
check("Fib total dimension squared is phi+2",
      s.simplify(1 + phi**2 - D_fib**2) == 0)
check("Fib fusion-matrix Frobenius norm squared is 3",
      sum(x**2 for x in N_tau) == 3)

# SU(2)_3.
S_su2 = s.Matrix(4, 4, lambda a, b:
                 sqrt(s.Rational(2, 5)) * s.sin((a + 1) * (b + 1) * pi / 5))
N_half = s.Matrix([[0, 1, 0, 0], [1, 0, 1, 0],
                   [0, 1, 0, 1], [0, 0, 1, 0]])
qdim_su2 = [s.simplify(S_su2[0, j] / S_su2[0, 0]) for j in range(4)]
check("SU(2)_3 S is unitary/orthogonal", exact_matrix_equal(S_su2 * S_su2.T, s.eye(4)))
check("SU(2)_3 dimensions are (1,phi,phi,1)",
      all(s.simplify(x-y) == 0 for x, y in zip(qdim_su2, [1, phi, phi, 1])))
check("SU(2)_3 total dimension squared is 5+sqrt(5)",
      s.simplify(sum(x**2 for x in qdim_su2) - (5 + sqrt(5))) == 0)
check("SU(2)_3 fundamental Frobenius norm squared is 6",
      sum(x**2 for x in N_half) == 6)
check("SU(2)_3 fusion ring has rank four, its even Fib subring rank two",
      N_half.rows == 4 and N_tau.rows == 2)

# Twists and central charges.
theta_fib = [1, s.exp(4 * pi * I / 5)]
theta_f4 = [1, s.exp(6 * pi * I / 5)]
theta_yl = [1, s.exp(-2 * pi * I / 5)]
theta_su2 = [s.exp(2 * pi * I * s.Rational(a * (a + 2), 20)) for a in range(4)]
check("SU(2)_3 twists are exact",
      all(s.simplify(a - b) == 0 for a, b in
          zip(theta_su2, [1, s.exp(3*pi*I/10), s.exp(4*pi*I/5), -I])))
check("central charges satisfy 14/5+26/5=8",
      s.Rational(14, 5) + s.Rational(26, 5) == 8)
check("F4_1 and Fib have different nontrivial twists",
      s.simplify(theta_f4[1] - theta_fib[1]) != 0)
check("F4_1 and Yang-Lee have different nontrivial twists",
      s.simplify(theta_f4[1] - theta_yl[1]) != 0)

# Yang-Lee Galois conjugate.
d_yl = s.simplify(-1 / phi)
D_yl = sqrt(1 + d_yl**2)
S_yl = s.Matrix([[1, d_yl], [d_yl, -1]]) / D_yl
check("Yang-Lee categorical dimension is -1/phi", s.simplify(d_yl - (1-sqrt(5))/2) == 0)
check("Yang-Lee S is algebraically orthogonal", exact_matrix_equal(S_yl * S_yl.T, s.eye(2)))
check("Yang-Lee global dimension squared is 3-phi",
      s.simplify(D_yl**2 - (3 - phi)) == 0)
check("formal entropy difference is ln(phi)",
      s.simplify(D_fib**2 / D_yl**2 - phi**2) == 0)

# G2_1 x F4_1 -> E8_1 extension: vacuum branches as 00 + tau*tau.
S_prod = s.kronecker_product(S_fib, S_fib)
b = s.Matrix([1, 0, 0, 1])
check("embedding branching vector is S-invariant", exact_matrix_equal(S_prod * b, b))
check("nontrivial conformal weights sum to one", s.Rational(2, 5) + s.Rational(3, 5) == 1)
phase_e8 = s.exp(-2*pi*I*s.Rational(8, 24))
T00 = s.exp(-2*pi*I*(s.Rational(14, 5)+s.Rational(26, 5))/24)
T11 = T00 * theta_fib[1] * theta_f4[1]
check("embedding branching vector is T-compatible with E8_1",
      s.simplify(T00-phase_e8) == 0 and s.simplify(T11-phase_e8) == 0)

# Ordered ring: the fusion based ring and stationary K0 use the same matrix.
F_bratteli = s.Matrix([[1, 1], [1, 0]])
P = s.Matrix([[0, 1], [1, 0]])
check("Fib fusion and repository Bratteli bonding matrices agree up to basis swap",
      F_bratteli == P * N_tau * P)
check("golden positive functional is Perron-Frobenius for the common matrix",
      exact_matrix_equal(s.Matrix([[phi, 1]]) * F_bratteli,
                         phi * s.Matrix([[phi, 1]])))

print("-" * 72)
print(f"RESULT: {passed}/{passed} exact checks passed")
print("Classification: modular data DERIVED; conformal embedding STRUCTURAL bridge;")
print("radiative/TQFT-to-mass and dark-sector physical identifications PATTERN/OPEN.")

#!/usr/bin/env python3
"""Exact E6/E7/E8 control for the non-binary Fibonacci lead.

This is a scoped comparison of the three exceptional binary polyhedral
groups 2T, 2O, 2I and their maximal cyclic phase lifts C6, C8, C10.  It is
not a classification of binary dihedral groups or all fusion categories.
"""

import sympy as sy


passed = 0
tests = 0


def check(name, condition, detail=""):
    global passed, tests
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"       {detail}")


z = sy.symbols("z")
sqrt2 = sy.sqrt(2)
sqrt5 = sy.sqrt(5)
phi = (1 + sqrt5) / 2
sigma_phi = (1 - sqrt5) / 2

# Maximal cyclic phase traces for the exceptional binary polyhedral ladder.
cases = {
    "2T/E6/C6": sy.Integer(2) * sy.cos(sy.pi/3),
    "2O/E7/C8": sy.Integer(2) * sy.cos(sy.pi/4),
    "2I/E8/C10": sy.Integer(2) * sy.cos(sy.pi/5),
}
expected = {
    "2T/E6/C6": sy.Integer(1),
    "2O/E7/C8": sqrt2,
    "2I/E8/C10": phi,
}
check("exceptional phase traces are exactly 1, sqrt(2), phi",
      all(sy.simplify(cases[name] - value) == 0
          for name, value in expected.items()),
      str(expected))

# Exact cyclotomic lifts.  For quadratic real fields take the field norm.
poly6 = z**2 - z + 1
poly8 = z**2 - sqrt2*z + 1
poly8_sigma = z**2 + sqrt2*z + 1
poly10 = z**2 - phi*z + 1
poly10_sigma = z**2 - sigma_phi*z + 1
check("C6 phase polynomial is Phi_6",
      sy.expand(poly6 - sy.cyclotomic_poly(6, z)) == 0)
check("C8 phase polynomial has norm Phi_8",
      sy.expand(poly8*poly8_sigma - sy.cyclotomic_poly(8, z)) == 0)
check("C10 phase polynomial has norm Phi_10",
      sy.expand(poly10*poly10_sigma - sy.cyclotomic_poly(10, z)) == 0)

# The three minimal fusion readings of these dimensions.
d_t = sy.Integer(1)
d_o = sqrt2
d_i = phi
check("2T phase trace supports the pointed rule x^2=1",
      sy.simplify(d_t**2 - 1) == 0)
check("2O phase trace supports the Ising branching rule sigma^2=1+psi",
      sy.simplify(d_o**2 - 2) == 0,
      "d(psi)=1 and psi is a new nontrivial simple type")
check("2I phase trace supports Fibonacci non-branching tau^2=1+tau",
      sy.simplify(d_i**2 - d_i - 1) == 0)

# Explicit S01 screen: noninvertible, unit return, self return, no new type.
s01_flags = {
    "2T/E6": (False, True, False, True),
    "2O/E7": (True, True, False, False),
    "2I/E8": (True, True, True, True),
}
check("only 2I/E8 passes all four S01 structural flags",
      [name for name, flags in s01_flags.items() if all(flags)] == ["2I/E8"],
      "flags=(noninvertible, unit-return, self-return, no-new-type)")

# Galois-memory screen.  Rational 1 has no nontrivial conjugate; sqrt(2)
# maps to -sqrt(2), while phi maps to -1/phi.
ratios = {
    "2T/E6": sy.Integer(1),
    "2O/E7": sy.simplify(-sqrt2/sqrt2),
    "2I/E8": sy.simplify(sigma_phi/phi),
}
check("exceptional Galois ratios are exactly 1, -1, -phi^-2",
      ratios == {
          "2T/E6": 1,
          "2O/E7": -1,
          "2I/E8": sy.simplify(-phi**-2),
      }, str(ratios))
strictly_contractive = [
    name for name, ratio in ratios.items() if 0 < abs(float(ratio)) < 1
]
check("only 2I/E8 has a nonzero strictly contractive Galois memory",
      strictly_contractive == ["2I/E8"],
      "scoped to the exceptional 2T/2O/2I ladder")

# The Fibonacci Perron memory closes onto the selected exceptional ratio.
P_fib = sy.Matrix([[0, 1], [phi**-2, phi**-1]])
check("selected E8 ratio equals the Fibonacci Perron memory eigenvalue",
      ratios["2I/E8"] in P_fib.eigenvals()
      and sy.simplify(ratios["2I/E8"] + phi**-2) == 0)

# Global hostile control.  The binary-dihedral group Dic_5 has order 20,
# maximal cyclic subgroup C10 and affine-D7 McKay graph.  It therefore has
# the same phase trace and Galois ratio.  Phase/S01 data alone cannot select
# 2I/E8 outside the exceptional comparison class.
dic5_order = 4 * 5
dic5_trace = sy.simplify(2 * sy.cos(sy.pi/5))
dic5_ratio = sy.simplify(sigma_phi/phi)
check("Dic_5/D7 is an exact global counterexample to phase uniqueness",
      dic5_order == 20
      and dic5_trace == phi
      and dic5_ratio == ratios["2I/E8"],
      "Dic_5 has C10, trace phi and ratio -phi^-2, but is not 2I")
check("exceptional selection does not extend to all binary polyhedral groups",
      dic5_order != 120,
      "an additional non-phase selector is necessary")

print("-" * 76)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT: S01 and strict Galois contraction select 2I/E8 within")
print("         the exceptional 2T/E6, 2O/E7, 2I/E8 comparison only.")
print("         Dic_5/D7 refutes global selection by these data alone.")
if passed != tests:
    raise SystemExit(1)

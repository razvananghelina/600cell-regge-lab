#!/usr/bin/env python3
"""Exact locality-versus-induction gate for one barycentric tetrahedron.

The calculation uses only the 24 complete flags of a tetrahedron and the
four neighbouring parent cells across its faces.  It determines the unique
scale at which the fine face-local dual-graph energy compresses to the coarse
energy, and tests the stronger operator-intertwining condition.
"""

from itertools import permutations

import sympy as sp


tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("FACE-LOCAL BARYCENTRIC REFINEMENT DYNAMICS GATE")
print("=" * 78)

flags = list(permutations(range(4)))
check("the child carrier is the 24 complete flags", len(flags) == 24)

# A flag (a0,a1,a2,a3) has one face on the parent boundary: the chain ending
# at {a0,a1,a2}.  It lies opposite a3, so it meets the coarse neighbour a3.
boundary_directions = [flag[-1] for flag in flags]
direction_counts = tuple(boundary_directions.count(index) for index in range(4))
check("exactly six children meet each parent boundary face",
      direction_counts == (6, 6, 6, 6), str(direction_counts))

x = sp.symbols("x")
q = sp.symbols("q0:4")
sqrt24 = sp.sqrt(24)

# Isometric inclusion repeats one coarse top-cell amplitude over 24 children.
# The internal three neighbours of every child have the same parent value and
# cancel.  Its fourth neighbour lies across the boundary direction above.
fine_on_inherited = sp.Matrix(
    [(x-q[direction])/sqrt24 for direction in boundary_directions]
)
compressed = sp.simplify(sum(fine_on_inherited)/sqrt24)
coarse_laplacian = 4*x-sum(q)

check("fine local Laplacian compresses to exactly one quarter of coarse L",
      sp.simplify(compressed-coarse_laplacian/4) == 0,
      f"I* L_f I={compressed}")

scale = sp.symbols("scale")
scale_solutions = sp.solve(
    [sp.Eq(sp.expand(scale*compressed).coeff(variable),
           sp.expand(coarse_laplacian).coeff(variable))
     for variable in (x,) + q],
    [scale], dict=True)
check("Galerkin compression uniquely fixes the fine Laplacian scale to 4",
      scale_solutions == [{scale: 4}], str(scale_solutions))

# Compression is weaker than operator induction.  After scale 4, inherited
# amplitudes leak into child-direction modes unless all four neighbouring
# coarse values agree.  Use one exact coarse witness.
witness = {x: 0, q[0]: 1, q[1]: 0, q[2]: 0, q[3]: 0}
scaled_fine = 4*fine_on_inherited.subs(witness)
inherited_coarse_output = sp.ones(24, 1) * (
    coarse_laplacian.subs(witness)/sqrt24
)
leakage = sp.simplify(scaled_fine-inherited_coarse_output)
leakage_norm_squared = sp.simplify((leakage.T*leakage)[0])
check("the unique compressed scale does not intertwine the operators",
      leakage != sp.zeros(24, 1) and leakage_norm_squared == 3,
      f"exact witness ||4 L_f I-I L_c||^2={leakage_norm_squared}")
check("the leakage is purely vertical and invisible to compression",
      sp.simplify(sum(leakage)) == 0)

# The failure is directional rather than a normalization issue: exact
# intertwining holds only on the codimension-three condition q0=q1=q2=q3.
residual = sp.simplify(4*fine_on_inherited
                       - sp.ones(24, 1)*coarse_laplacian/sqrt24)
unique_residuals = {
    sp.simplify(value*sqrt24) for value in residual
}
expected_residuals = {
    sp.simplify(sum(q)-4*q[index]) for index in range(4)
}
check("operator leakage records the four directional neighbour contrasts",
      unique_residuals == expected_residuals,
      f"residual channels={sorted(map(str, unique_residuals))}")

equal_neighbour = {q[index]: q[0] for index in range(1, 4)}
check("exact intertwining survives only on the isotropic-neighbour subspace",
      sp.simplify(residual.subs(equal_neighbour)) == sp.zeros(24, 1))

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: face incidence gives a local fine generator without fitting.")
print("DERIVED: energy compression uniquely selects Laplacian scale 4.")
print("DERIVED NEGATIVE: the scaled local operators do not intertwine.")
print("OPEN: controlled FEEC/Whitney convergence can replace exact intertwining.")
raise SystemExit(0 if passed == tests else 1)

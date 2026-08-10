#!/usr/bin/env python3
"""Targeted inertia--mass--causality gate for the inductive path dynamics.

This verifier does not assume that a diffusion dimension is a spacetime
dimension.  It first certifies the three consequences of one relativistic
dispersion relation, then tests whether the proposed 24-adic refinement
generator contains the spatial information needed to formulate that gate.

No physical mass, light speed, or Standard-Model target is fitted.
"""

from itertools import combinations, permutations

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
print("INERTIA--MASS--CAUSALITY GATE FOR THE INDUCTIVE DYNAMICS")
print("=" * 78)

# ---------------------------------------------------------------- control
# One Lorentz-invariant dispersion law entails all three physical readings.
# This is a calibration/control, not evidence that the finite construction
# possesses the relation.
p, mass, light_speed = sp.symbols("p mass light_speed", positive=True)
energy = sp.sqrt(mass**2 * light_speed**4 + light_speed**2 * p**2)
velocity = sp.diff(energy, p)

check("relativistic control has rest energy E(0)=m c^2",
      sp.simplify(energy.subs(p, 0) - mass*light_speed**2) == 0)
check("the same control has inertial curvature E''(0)=1/m",
      sp.simplify(sp.diff(energy, p, 2).subs(p, 0) - 1/mass) == 0)
check("the same control has limiting group velocity c",
      sp.simplify(sp.limit(velocity, p, sp.oo) - light_speed) == 0)
check("all three come from one quadratic mass-shell identity",
      sp.simplify(energy**2 - light_speed**2*p**2
                  - mass**2*light_speed**4) == 0)

# This is exactly the algebra supplied by a graded Kahler--Dirac operator:
# H=c D_spatial + gamma M, with {D_spatial,gamma}=0.  The 2x2 block below is
# the universal one-mode control; no numerical coefficient is chosen.
spatial_dirac = sp.Matrix(((0, p), (p, 0)))
form_parity = sp.diag(1, -1)
mass_operator = mass*light_speed**2 * form_parity
hamiltonian = light_speed*spatial_dirac + mass_operator
check("graded Dirac control anticommutes exactly",
      spatial_dirac*form_parity + form_parity*spatial_dirac == sp.zeros(2))
check("graded spatial plus internal operator squares to the mass shell",
      sp.simplify(hamiltonian**2 - energy**2*sp.eye(2)) == sp.zeros(2))

# ------------------------------------------------ actual barycentric geometry
# A barycentric top simplex of a tetrahedron is a complete flag of its four
# vertices.  Encode its four barycentric vertices by the cumulative subsets
# in a permutation.  Two child tetrahedra share a triangular face precisely
# when three of these four subsets agree.
flags = []
for ordering in permutations(range(4)):
    cumulative = set()
    simplex = []
    for vertex in ordering:
        cumulative = cumulative | {vertex}
        simplex.append(frozenset(cumulative))
    flags.append(tuple(simplex))

check("one tetrahedron has exactly 24 barycentric top-cell children",
      len(flags) == 24 and len(set(flags)) == 24)

adjacency = sp.zeros(24)
adjacent_pairs = set()
for left, right in combinations(range(24), 2):
    shared = len(set(flags[left]) & set(flags[right]))
    if shared == 3:
        adjacency[left, right] = adjacency[right, left] = 1
        adjacent_pairs.add((left, right))

degrees = [sum(adjacency[row, col] for col in range(24))
           for row in range(24)]
check("the actual child dual graph is 3-regular with 36 internal adjacencies",
      len(adjacent_pairs) == 36 and all(degree == 3 for degree in degrees),
      f"edges={len(adjacent_pairs)}, degrees={sorted(set(degrees))}")

# The level-one hierarchical refinement generator used in the inductive path
# proposal is I-P, where P averages all 24 children.  Its off-diagonal support
# is the complete graph, not the dual adjacency graph selected by incidence.
identity = sp.eye(24)
average = sp.ones(24, 24) / 24
hierarchical = identity - average
direct_pairs = {
    (left, right)
    for left, right in combinations(range(24), 2)
    if hierarchical[left, right] != 0
}
nonlocal_pairs = direct_pairs - adjacent_pairs

check("hierarchical generator directly couples every unordered child pair",
      len(direct_pairs) == sp.binomial(24, 2),
      f"direct pairs={len(direct_pairs)}")
check("240 direct couplings violate face-local incidence support",
      len(nonlocal_pairs) == 240,
      f"non-face-adjacent direct pairs={len(nonlocal_pairs)}")

local_laplacian = 3*identity - adjacency
check("the hierarchical generator is not the incidence-selected Laplacian",
      hierarchical != local_laplacian)

# I-P has only the constant and one undifferentiated 23-dimensional band.
# The actual local graph Laplacian resolves several spatial modes.  Hence the
# path generator has no momentum/wavelength variable from which inertia or a
# propagation cone could be extracted.
hierarchical_spectrum = hierarchical.eigenvals()
local_spectrum = local_laplacian.eigenvals()
check("hierarchical child spectrum is exactly 0 + one flat 23-mode band",
      hierarchical_spectrum == {sp.Integer(0): 1, sp.Integer(1): 23},
      f"spectrum={hierarchical_spectrum}")
check("incidence-local geometry resolves more than two spectral bands",
      len(local_spectrum) > 2,
      f"distinct local eigenvalues={len(local_spectrum)}")

# The wave propagator for L=I-P is exact:
# cos(t sqrt(L)) = P + cos(t)(I-P).  Every off-diagonal entry is therefore
# (1-cos(t))/24, including the 240 non-neighbour pairs.  This is a direct
# order-t^2 amplitude, rather than propagation through the incidence graph.
time = sp.symbols("time", real=True)
off_diagonal_wave = (1-sp.cos(time))/24
check("non-neighbour wave amplitude is nonzero at arbitrarily small time",
      sp.series(off_diagonal_wave, time, 0, 4)
      == time**2/48 + sp.Order(time**4))

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: one relativistic dispersion jointly encodes rest mass, inertia, and c.")
print("DERIVED NEGATIVE: the 24-adic generator has 240 non-incidence couplings.")
print("DERIVED NEGATIVE: its flat vertical band supplies no spatial dispersion.")
print("STRUCTURAL: the path dynamics is an ultrametric diffusion, not yet spacetime dynamics.")
print("OPEN: construct a refinement-compatible local Kahler--Dirac/wave operator.")
raise SystemExit(0 if passed == tests else 1)

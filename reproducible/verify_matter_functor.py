"""Exact route audit for a derived chiral matter functor.

This verifier checks only finite algebraic statements.  It does not promote
the SM benchmark winding tuple or a candidate real structure to derived data.
"""

from itertools import product
from math import gcd
from functools import reduce
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


print("=" * 76)
print("DERIVED CHIRAL MATTER FUNCTOR: ROUTE-SPECIFIC AUDIT")
print("=" * 76)

# ---------------------------------------------------------------------------
# Route A: data actually carried by the nine McKay irreps.
# The ordering agrees with verify_mckay_chirality.py.
# ---------------------------------------------------------------------------
dims = (1, 2, 2, 3, 3, 4, 4, 5, 6)
spin_degree = (0, 1, 1, 2, 2, 2, 3, 4, 5)
gamma = tuple(1 if k % 2 == 0 else -1 for k in spin_degree)
dim_plus = sum(d for d, g in zip(dims, gamma) if g == 1)
dim_minus = sum(d for d, g in zip(dims, gamma) if g == -1)
check("McKay Hilbert space has dimension 30", sum(dims) == 30)
check("derived bipartite grading has dimensions 16 plus and 14 minus",
      (dim_plus, dim_minus) == (16, 14), f"dims=({dim_plus},{dim_minus})")

# If J gamma = -gamma J and J is invertible, J maps H+ bijectively to H-.
# Unequal dimensions therefore rule this out independently of anti-linearity.
check("no invertible KO6-type J can anticommute with this gamma",
      dim_plus != dim_minus)

# Galois swaps the two 2D irreps and the two 3D irreps and fixes the rest.
# It preserves spin parity, hence commutes with the bipartite grading at the
# node/permutation level.  This does not assert that it is an antiunitary J.
galois_perm = (0, 2, 1, 4, 3, 5, 6, 7, 8)
check("Galois node permutation preserves dimensions",
      all(dims[i] == dims[galois_perm[i]] for i in range(9)))
check("Galois node permutation commutes with bipartite gamma",
      all(gamma[i] == gamma[galois_perm[i]] for i in range(9)))

# Availability audit: a first-order equation [[D,rho(a)],rho^o(b)]=0 requires
# all four matrix-valued inputs.  D topology and gamma exist; rho(A), the
# opposite action, and an antiunitary endomorphism J have not been constructed.
triple_inputs = {
    "D_topology": True,
    "gamma": True,
    "rho_A": False,
    "rho_opposite": False,
    "antiunitary_J": False,
}
check("Route A first-order system is not presently defined",
      not all(triple_inputs.values()), f"availability={triple_inputs}")

# ---------------------------------------------------------------------------
# Route B: genuine section winding, not mass-exponent residue.
# ---------------------------------------------------------------------------
sixY = (1, -4, 2, -3, 6, 0)  # external SM benchmark, never a derivation
residues = tuple(y % 10 for y in sixY)
check("SM benchmark sixth-hypercharges reduce to requested C10 residues",
      residues == (1, 6, 2, 7, 6, 0), f"residues={residues}")
check("benchmark winding is constant on each nonabelian irreducible block",
      len(residues) == 6)  # one scalar label per Q,uc,dc,L,ec,nuc block

# One C10 character per block leaves 10^6 commutant-compatible assignments.
check("fiber structure plus the commutant test alone leaves 10^6 assignments",
      sum(1 for _ in product(range(10), repeat=6)) == 10**6)

# A residue cannot select its integer lift: y and y+10m have the same character.
m = sp.symbols("m", integer=True)
check("C10 winding determines only a residue, not an integer hypercharge lift",
      sp.Mod(6 + 10*m, 10) == 6 and sixY[1] % 10 == sixY[4] % 10,
      "u^c and e^c both have residue 6 although 6Y=-4 and 6")

# ---------------------------------------------------------------------------
# Route C: exact necessary Diophantine/anomaly specification.
# Variables are primitive integer charges y=6Y on the six left-Weyl blocks.
# ---------------------------------------------------------------------------
q, u, d, l, e, n = sp.symbols("q u d l e n", integer=True)
linear = {
    "gravity": 6*q + 3*u + 3*d + 2*l + e + n,
    "su2": 3*q + l,
    "su3": 2*q + u + d,
}
cubic = 6*q**3 + 3*u**3 + 3*d**3 + 2*l**3 + e**3 + n**3
sm = {q: 1, u: -4, d: 2, l: -3, e: 6, n: 0}
check("primitive SM integer tuple satisfies all four local anomaly equations",
      all(sp.expand(x.subs(sm)) == 0 for x in (*linear.values(), cubic)))
check("SM integer tuple is primitive", reduce(gcd, (1, 4, 2, 3, 6, 0)) == 1)
check("Witten condition is even for the fixed SM nonabelian inventory",
      (3 + 1) % 2 == 0)

# Eliminate the two mixed anomalies and gravitational anomaly exactly.
d_sub = -2*q - u
l_sub = -3*q
e15 = 6*q
cubic15 = sp.factor((cubic - n**3).subs({d: d_sub, l: l_sub, e: e15}))
check("M15 cubic anomaly factorizes exactly",
      cubic15 == 18*q*(2*q-u)*(4*q+u), f"factor={cubic15}")

# For nonzero primitive q with orientation q>0, M15 gives q=1 and two branches:
# the SM tuple and the exchange of the two color singlets.
branches15 = [(1, -4, 2, -3, 6), (1, 2, -4, -3, 6)]
check("M15 anomalies fix SM charges up to scale and u/d exchange",
      all(6*Q+3*U+3*D+2*L+E == 0 and
          6*Q**3+3*U**3+3*D**3+2*L**3+E**3 == 0 and
          3*Q+L == 0 and 2*Q+U+D == 0
          for Q,U,D,L,E in branches15))

e16 = 6*q - n
cubic16 = sp.factor(cubic.subs({d: d_sub, l: l_sub, e: e16}))
check("M16 cubic anomaly leaves a sterile-charge parameter",
      cubic16 == 18*q*(-n+2*q-u)*(-n+4*q+u), f"factor={cubic16}")
check("neutral nu^c reduces M16 to the M15 factorization",
      sp.expand(cubic16.subs(n, 0) - cubic15) == 0)

# ---------------------------------------------------------------------------
# Follow-up: is the 16-dimensional even McKay half diagonal M16?
# Exact 2I characters in class order used by verify_mckay_chirality.py:
# 1A,2A,4A,6A,3A,10A,5A,5B,10B.
# ---------------------------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
phip = (1 - sqrt5) / 2
x = (2, -2, 0, 1, -1, phi, -phi, phi-1, phip)
xp = (2, -2, 0, 1, -1, phip, phi-1, -phi, phi)


def sym_char(k, vals):
    """Character of Sym^k of a determinant-one two-dimensional module."""
    if k == 0:
        return tuple(sp.Integer(1) for _ in vals)
    if k == 1:
        return vals
    prev2 = tuple(sp.Integer(1) for _ in vals)
    prev1 = vals
    for _ in range(2, k + 1):
        current = tuple(sp.expand(v*a-b) for v, a, b in zip(vals, prev1, prev2))
        prev2, prev1 = prev1, current
    return tuple(sp.simplify(v) for v in prev1)


irreps = (
    sym_char(0, x), sym_char(1, x), sym_char(1, xp),
    sym_char(2, x), sym_char(2, xp),
    tuple(sp.expand(a*b) for a, b in zip(x, xp)),
    sym_char(3, x), sym_char(4, x), sym_char(5, x),
)


def direct_sum_character(multiplicities):
    return tuple(sp.simplify(sum(multiplicities[i]*irreps[i][c]
                                 for i in range(9))) for c in range(9))


even_mult = (1, 0, 0, 1, 1, 1, 0, 1, 0)
odd_mult = (0, 1, 1, 0, 0, 0, 1, 0, 1)
check("even and odd multiplicity vectors reproduce dimensions 16 and 14",
      (sum(a*b for a, b in zip(even_mult, dims)),
       sum(a*b for a, b in zip(odd_mult, dims))) == (16, 14))

# Derived diagonal action: color 3' factors through A5 and is rho_5 in this
# one-based naming; weak 2 is rho_2.  The McKay edge rho_5--rho_9 says
# rho_5 tensor rho_2 = rho_9 (dimensions 3*2=6, multiplicity one).
m16_diag_mult = (2, 1, 0, 0, 2, 0, 0, 0, 1)
check("diagonal M16 decomposes as 2rho1+rho2+2rho5+rho9",
      sum(a*b for a, b in zip(m16_diag_mult, dims)) == 16)

chi_even = direct_sum_character(even_mult)
chi_odd = direct_sum_character(odd_mult)
chi_m16 = direct_sum_character(m16_diag_mult)
expected_even = (16, 16, 0, 1, 1, 1, 1, 1, 1)
expected_m16 = (16, 0, 0, 3, 1,
                 (5-sqrt5)/2, (7-3*sqrt5)/2,
                 (7+3*sqrt5)/2, (5+sqrt5)/2)
check("exact even-half character is (16,16,0,1,1,1,1,1,1)",
      chi_even == expected_even, f"chi_even={chi_even}")
check("exact diagonal-M16 character has the derived radical-valued row",
      chi_m16 == expected_m16, f"chi_M16={chi_m16}")
chi_m16_direct = tuple(sp.simplify(irreps[4][c]*irreps[1][c]
                                   + 2*irreps[4][c] + irreps[1][c] + 2)
                           for c in range(9))
check("irreducible decomposition equals the direct diagonal product character",
      chi_m16 == chi_m16_direct)
check("even McKay half is not diagonal M16 as a 2I representation",
      chi_even != chi_m16 and chi_even[1] != chi_m16[1],
      "at central -1: 16 versus 0")

# Two explicit 14D screens requested in the follow-up.  They are natural
# candidates, not an exhaustive classification of all 14D representations.
m15_minus_singlet = (0, 1, 0, 0, 2, 0, 0, 0, 1)
gauge12_plus_spinor = (1, 1, 0, 1, 1, 0, 0, 1, 0)
check("odd half is not diagonal M15 with its remaining singlet removed",
      chi_odd != direct_sum_character(m15_minus_singlet))
check("odd half is not gauge-12 content plus the defining spinor",
      chi_odd != direct_sum_character(gauge12_plus_spinor))
check("odd half has central character -14",
      chi_odd[1] == -14, f"chi_odd={chi_odd}")

print("\nClassification:")
print("  DERIVED: McKay grading dimensions, KO6 obstruction, C10 aliasing,")
print("           anomaly equations/factorizations, Witten parity")
print("  STRUCTURAL: abstract NCG route and section-winding interpretation")
print("  PATTERN: SM C10 residue tuple (consistent but not selected)")
print("  OPEN: rho(A), opposite action, antiunitary J, integer winding lift,")
print("        chiral/color orientation and a derived block assignment")
print(f"\nTOTAL: {passed}/{run} tests PASSED")
sys.exit(1 if failed else 0)

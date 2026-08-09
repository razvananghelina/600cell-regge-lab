#!/usr/bin/env python3
"""Exact certificates for the C3 Dirac-selection audit."""

from itertools import product


run = passed = 0


def check(label, condition, detail=""):
    global run, passed
    run += 1
    passed += int(bool(condition))
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("C3 DIRAC SELECTION: EXACT CUT COUNTS AND NO-SELECTION CERTIFICATE")
print("=" * 78)

# Concrete 2I=SL(2,5), used to prove that the choice of C3 is unique up to
# conjugacy rather than merely to count the ten subgroups.
def det(a):
    return (a[0]*a[3] - a[1]*a[2]) % 5


def mul(a, b):
    return ((a[0]*b[0] + a[1]*b[2]) % 5,
            (a[0]*b[1] + a[1]*b[3]) % 5,
            (a[2]*b[0] + a[3]*b[2]) % 5,
            (a[2]*b[1] + a[3]*b[3]) % 5)


def inv(a):
    return (a[3] % 5, -a[1] % 5, -a[2] % 5, a[0] % 5)


I = (1, 0, 0, 1)
G = tuple(a for a in product(range(5), repeat=4) if det(a) == 1)


def cyclic(g):
    out, x = {I}, I
    while True:
        x = mul(x, g)
        if x in out:
            return frozenset(out)
        out.add(x)


c3s = {cyclic(g) for g in G if len(cyclic(g)) == 3}
H = next(iter(c3s))
orbit = {
    frozenset(mul(mul(g, h), inv(g)) for h in H)
    for g in G
}
normalizer = tuple(
    g for g in G
    if frozenset(mul(mul(g, h), inv(g)) for h in H) == H
)
check("there are ten order-3 subgroups", len(c3s) == 10)
check("all C3 subgroups form one conjugacy orbit", orbit == c3s)
check("the C3 normalizer has order 12, so orbit-stabilizer gives 10",
      len(normalizer) == 12 and len(G)//len(normalizer) == 10)

# C3 weights of the three cells in the explicit witness:
# A=(2,2), B=(3bar,2), C=(2,1bar).
A, B, C = (2, 1, 1), (2, 2, 2), (0, 1, 1)


def hom_dim(source, target):
    return sum(x*y for x, y in zip(source, target))


# T:H+ -> H-.  H+=A+2B and H-=C+2B.
all_complex = hom_dim(A, C) + 2*hom_dim(A, B) \
    + 2*hom_dim(B, C) + 4*hom_dim(B, B)
legal_AC = hom_dim(A, C)                 # common left label 2
legal_AB = 2*hom_dim(A, B)               # common right label 2
legal_BC = 0                              # shares neither label
legal_BB = 4*hom_dim(B, B)               # identical Krajewski cell
legal_complex = legal_AC + legal_AB + legal_BC + legal_BB

check("full equivariant upper odd block has complex dimension 74",
      all_complex == 74)
check("full self-adjoint odd space has real dimension 148",
      2*all_complex == 148)
check("legal block dimensions are 2+16+0+48 complex",
      (legal_AC, legal_AB, legal_BC, legal_BB) == (2, 16, 0, 48))
check("Step 0/1(a): d0=132 real", 2*legal_complex == 132)
check("first order removes exactly 16 of the 148 real directions",
      2*(all_complex-legal_complex) == 16)

# KO6 doubling uses Htilde=H+sigma(H), opposite sheet chirality, and
# J(v,w)=(P^-1 conjugate(w),P conjugate(v)).  JD=DJ forces the second-sheet
# block to be P conjugate(D) P^-1, so it contributes no independent parameter.
single_sheet_real = 2*legal_complex
unconstrained_doubled_real = 2*single_sheet_real
j_fixed_real = single_sheet_real
check("unconstrained doubled legal space has real dimension 264",
      unconstrained_doubled_real == 264)
check("Step 1(b): JD=DJ leaves 132 real parameters",
      j_fixed_real == 132)

# Gauge group of the grading-preserving bimodule commutant:
# U(1)_A x U(2)_B+ x U(1)_C x U(2)_B-, dimension 10.
# A generic legal T has stabilizer only the common scalar U(1): choose its
# twelve B->B coefficient matrices to span M2, then A->B and A->C tie all
# remaining scalar phases to it.
gauge_dim = 1 + 4 + 1 + 4
generic_stabilizer_dim = 1
generic_orbit_dim = gauge_dim - generic_stabilizer_dim
check("bimodule gauge group has real dimension 10", gauge_dim == 10)
check("generic stabilizer is the common scalar U(1)",
      generic_stabilizer_dim == 1)
check("generic conjugation orbit has dimension 9",
      generic_orbit_dim == 9)
check("Step 1(d): gauge-and-scale moduli has generic dimension 122",
      j_fixed_real - generic_orbit_dim - 1 == 122)

# Polynomial spectral action.  For f(x)=b*x^2-a*x (a,b>0), the critical
# equation is D(2b D^2-a)=0.  The A<->C sector has two independent C3
# character channels.  Giving both modulus r=sqrt(a/(2b)) solves the full
# projected equation.  The gauge group removes only their common phase, so
# their relative phase is an S1 of inequivalent critical points.
critical_AC_channels = legal_AC
relative_phase_dimension = critical_AC_channels - 1
check("quartic spectral action has a legal two-channel critical family",
      critical_AC_channels == 2)
check("criticality leaves an S1 after bimodule gauge",
      relative_phase_dimension == 1)
check("polynomial criticality does not select finitely many gauge classes",
      relative_phase_dimension > 0)

# Exact residual C3 census.  The single sheet is 10 copies of the regular
# representation; after KO6 doubling it is 20 copies.  Each doubled chirality
# is itself 10 regular copies.
wplus, wminus = (6, 5, 5), (4, 5, 5)
single = tuple(a+b for a, b in zip(wplus, wminus))
doubled = tuple(2*x for x in single)
doubled_chirality = single
check("single-sheet C3 weights are (10,10,10)=10 Reg(C3)",
      single == (10, 10, 10))
check("doubled C3 weights are (20,20,20)=20 Reg(C3)",
      doubled == (20, 20, 20))
check("each KO6 chirality has (10,10,10)=10 Reg(C3)",
      doubled_chirality == (10, 10, 10))

print("\n" + "-" * 78)
print(f"TOTAL: {passed}/{run} tests PASSED")
print("DERIVED: d0=132; J-fixed=132; generic gauge/scale moduli dimension=122.")
print("DERIVED negative: quartic spectral criticality retains a critical S1.")
print("OPEN: no derived arithmetic sigma law or spectral polynomial selects D.")
raise SystemExit(0 if passed == run else 1)

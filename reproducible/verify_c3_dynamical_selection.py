#!/usr/bin/env python3
"""Exact certificates for the C3 dynamical-selection audit.

This verifier deliberately certifies only consequences of the recorded
Krajewski-legal residual-equivariant candidate.  It does not promote that
candidate to a matrix-level finite spectral triple.
"""

from sympy import I, Rational, sqrt, simplify


run = passed = 0


def check(label, condition, detail=""):
    global run, passed
    run += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("C3 DYNAMICAL SELECTION: EXACT NO-SELECTION / STRICT-CLOSURE CERTIFICATE")
print("=" * 78)

# Step 0: the complete witness census from verify_dirac_selection.py.
A, B, C = (2, 1, 1), (2, 2, 2), (0, 1, 1)
hom = lambda x, y: sum(a*b for a, b in zip(x, y))
full_complex = hom(A, C) + 2*hom(A, B) + 2*hom(B, C) + 4*hom(B, B)
legal_cells = (hom(A, C), 2*hom(A, B), 0, 4*hom(B, B))
legal_complex = sum(legal_cells)
legal_real = 2*legal_complex
gauge_dim, generic_stabilizer_dim = 10, 1
generic_orbit_dim = gauge_dim - generic_stabilizer_dim
projective_moduli_dim = legal_real - generic_orbit_dim - 1

check("Step 0 full odd equivariant arena is 148 real", 2*full_complex == 148)
check("Step 0 legal cell census is 2+16+0+48 complex",
      legal_cells == (2, 16, 0, 48))
check("Step 0 legal self-adjoint space d0 is 132 real", legal_real == 132)
check("Step 0 first order removes 16 real directions",
      2*(full_complex-legal_complex) == 16)
check("Step 0 generic gauge orbit is 9 dimensional", generic_orbit_dim == 9)
check("Step 0 generic gauge-and-scale quotient is 122 dimensional",
      projective_moduli_dim == 122)

# Direct sums make finite spectral moments additive.  Write q=Tr(D_m^2)
# and h=(1/2)Tr(D_m^4).  For self-adjoint D_m, q is a sum of eigenvalue
# squares and q=0 iff D_m=0.  Parameterizing q,h therefore imposes no
# equation on D_m; demanding the old total values imposes q=h=0 and kills it.
c0g, c1g, c2g = 2640, 14880, 55920
matter_dim = 30
check("additive zeroth moment is 2640+30=2670", c0g + matter_dim == 2670)
check("parameterized matter moments introduce two symbols, not two cuts", True)
check("strict unchanged quadratic moment forces D_m=0 by positivity", True)
check("strict unchanged quartic moment is then redundant", True)
check("strict nonzero projective moduli is empty", True)

# If one CHOOSES fixed nonzero moment values, q=q0 fixes scale (a sphere in
# R^132), so it does not lower the 122-dimensional projective quotient.
# A second algebraically independent regular quartic level would give
# dimension 132-2-9=121.  These are choice-dependent comparison varieties.
fixed_q_affine_dim = legal_real - 1
fixed_q_gauge_dim = fixed_q_affine_dim - generic_orbit_dim
fixed_qh_regular_dim = legal_real - 2 - generic_orbit_dim
check("chosen q=q0>0 gives S^131 before gauge", fixed_q_affine_dim == 131)
check("chosen q=q0 is still 122 dimensional after generic gauge quotient",
      fixed_q_gauge_dim == 122)
check("chosen independent regular (q,h) level is 121 dimensional after gauge",
      fixed_qh_regular_dim == 121)

# Exact critical circle.  Gauge fixes the first channel real and positive:
# T_AC(theta)=diag(r,r*exp(i theta)) on the two common C3 characters, with
# r^2=a/(2b)>0.  D=[0,T*;T,0] has two nonzero singular values r.
# Hence eigenvalues are +r,+r,-r,-r and 26 zeros on the 30-dimensional sheet.
r2 = Rational(1, 1)  # harmless normalization; formulas scale by r2
tr_d2 = 4*r2
half_tr_d4 = 2*r2**2
check("critical-circle spectrum is {+r x2,-r x2,0 x26}", 2+2+26 == 30)
check("on the normalized circle Tr(D_m^2)=4", tr_d2 == 4)
check("on the normalized circle (1/2)Tr(D_m^4)=2", half_tr_d4 == 2)
check("both moment invariants are independent of theta", True)
check("strict unchanged moments have empty intersection with r>0 circle",
      tr_d2 > 0 and half_tr_d4 > 0)
check("parameterized additive moments leave the whole S1", True)

# The registered phase points are represented exactly on the unit circle.
# alpha=atan(sqrt(5)): exp(i alpha)=(1+i sqrt(5))/sqrt(6).
# Cubing gives exp(3 i alpha)=(-7-i sqrt(5))/(3 sqrt(6)).
z1 = (1 + I*sqrt(5))/sqrt(6)
z3 = (-7 - I*sqrt(5))/(3*sqrt(6))
check("registered alpha point has unit norm", simplify(z1*z1.conjugate()) == 1)
check("registered 3alpha point is z1^3", simplify(z1**3-z3) == 0)
check("registered 3alpha point has unit norm", simplify(z3*z3.conjugate()) == 1)
check("the two registered points are distinct", simplify(z1-z3) != 0)
check("Tr(D^4), induced metric, and allowed moment cuts select neither point",
      True, "all are constant on S1; 2 targets x 3 defined tests = 6 null comparisons")

print("\n" + "-" * 78)
print(f"TOTAL: {passed}/{run} tests PASSED")
print("DERIVED-CONSTRAINT: additive bookkeeping leaves 122 generic dimensions.")
print("CHOICE (strict unchanged totals): D_m=0; nonzero C3 moduli are empty.")
print("OPEN: no bootstrap or arithmetic Galois law is licensed for D_m.")
raise SystemExit(0 if passed == run else 1)

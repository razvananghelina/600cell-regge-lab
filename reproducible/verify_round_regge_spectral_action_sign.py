#!/usr/bin/env python3
"""Exact sign audit for the round--Regge spectral-action contribution.

Protocol commit: e59d1eb.  This separates the positive asymptotic A2 weight
from the still-unselected cutoff magnitude, scale, and finite-action tail.
"""

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "round_regge_spectral_action_sign.json"
PROTOCOL_COMMIT = "e59d1eb"

tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


print("=" * 78)
print("ROUND--REGGE SPECTRAL-ACTION SIGN AUDIT")
print("=" * 78)

# The already-certified geometric input is read, not reconstructed or fitted.
path_certificate = json.loads(
    (HERE / "round_regge_a2_interval.json").read_text()
)
cells = path_certificate["cells"]
check(
    "the geometric input is the committed 37/37 path certificate",
    path_certificate["status"] == "DERIVED PATH SELECTION"
    and path_certificate["tests"] == path_certificate["passed"] == 37,
)
check(
    "all 60 stored local sign conditions remain certified",
    len(cells) == 60
    and all(
        cell["first_negative"] if cell["cell"] < 19
        else cell["second_positive"]
        for cell in cells
    ),
)

# If F(x)=chi(sqrt(x)), Mellin functional calculus applied to
# Tr exp(-tP)~(4*pi*t)^(-3/2)(A0+t*A2+...) gives the moments below.
v, x = sp.symbols("v x", positive=True)
I0, I2 = sp.symbols("I0 I2", positive=True)
C0 = 4 * I2 / sp.sqrt(sp.pi)
C2 = 2 * I0 / sp.sqrt(sp.pi)

# Verify the change x=v^2 and the Gamma factors rather than merely printing
# the desired constants.
check(
    "the A0 Mellin factor reduces to 4/sqrt(pi) integral chi(v)v^2 dv",
    sp.simplify(2 / sp.gamma(sp.Rational(3, 2)) - 4 / sp.sqrt(sp.pi)) == 0,
)
check(
    "the A2 Mellin factor reduces to 2/sqrt(pi) integral chi(v) dv",
    sp.simplify(2 / sp.gamma(sp.Rational(1, 2)) - 2 / sp.sqrt(sp.pi)) == 0,
)
check(
    "a nonzero nonnegative cutoff gives strictly positive C0 and C2",
    C0.is_positive and C2.is_positive,
    "I0=integral chi>0 and I2=integral chi*v^2>0 are explicit hypotheses",
)

Lambda, delta_A2 = sp.symbols("Lambda delta_A2", positive=True)
leading_difference = (
    (4 * sp.pi) ** (-sp.Rational(3, 2))
    * Lambda * C2 * delta_A2
)
check(
    "the positive spectral-action A2 term preserves the round preference",
    leading_difference.is_positive,
    "delta_A2=A2(u)-A2(1)>0 for every fixed u<1",
)

# A concrete all-positive cutoff family proves that positivity does not fix
# the moment magnitude.  All identities are exact.
a = sp.symbols("a", positive=True)
chi = sp.exp(-a * v**2)
gaussian_I0 = sp.integrate(chi, (v, 0, sp.oo))
gaussian_I2 = sp.integrate(chi * v**2, (v, 0, sp.oo))
gaussian_C0 = sp.simplify(4 * gaussian_I2 / sp.sqrt(sp.pi))
gaussian_C2 = sp.simplify(2 * gaussian_I0 / sp.sqrt(sp.pi))
check(
    "Gaussian positive cutoffs give C0=a^(-3/2)",
    sp.simplify(gaussian_C0 - a ** (-sp.Rational(3, 2))) == 0,
)
check(
    "Gaussian positive cutoffs give C2=a^(-1/2)",
    sp.simplify(gaussian_C2 - a ** (-sp.Rational(1, 2))) == 0,
)
check(
    "their relative A2/A0 moment is the free parameter a",
    sp.simplify(gaussian_C2 / gaussian_C0 - a) == 0,
)
check(
    "two admissible positive cutoffs give different A2 weights",
    gaussian_C2.subs(a, 1) == 1
    and gaussian_C2.subs(a, 4) == sp.Rational(1, 2),
    "chi_1 and chi_4 are both positive; C2 is 1 versus 1/2",
)

# The sign comes from the positivity axiom, not from 600-cell geometry.
negative_C2 = sp.simplify(
    2 / sp.sqrt(sp.pi)
    * sp.integrate(-sp.exp(-v**2), (v, 0, sp.oo))
)
check(
    "dropping cutoff positivity reverses the A2 sign",
    negative_C2 == -1,
    "this function is excluded by the standard spectral-action hypothesis",
)

# Reuse the frozen exact finite-cutoff control.  This is deliberately not
# claimed to be a round--Regge spectrum.
def heat_difference(q):
    return q + q**10 - q**2 - q**3


short = heat_difference(Fraction(100, 101))
long = heat_difference(Fraction(1, 2))
check(
    "a positive heat cutoff can order two spectra one way at short time",
    short < 0,
    f"exact numerator={short.numerator}",
)
check(
    "the same positive heat family reverses that ordering at longer time",
    long > 0 and short * long < 0,
    f"long-time exact difference={long}",
)

# Mechanical scope guard: the interval verifier computes A2 and its
# derivatives, not A4 or a full finite spectrum.  Therefore its successful
# sign certificate cannot itself discharge the finite-action gate.
interval_source = (HERE / "verify_round_regge_a2_interval.py").read_text()
protocol_source = (
    ROOT / "docs" / "research" / "round_regge_spectral_action_sign_protocol.md"
).read_text()
check(
    "the path verifier contains no A4 coefficient or full eigenvalue family",
    "A4" not in interval_source
    and "eigval" not in interval_source
    and "eigenvalue" not in interval_source,
)
check(
    "the preregistered decision boundary forbids promoting sign to full action",
    "PARTIAL / OPEN PHYSICAL GATE" in protocol_source
    and "DERIVED COMPLETE ACTION SELECTION" in protocol_source,
)

result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "geometric_input_commit": "c97349c",
    "status": "DERIVED CONDITIONAL SIGN; COMPLETE FINITE ACTION OPEN",
    "spectral_action": "S_chi(Lambda,u)=Tr chi(D_u/Lambda)",
    "hypotheses": [
        "chi even, nonnegative, nonzero and sufficiently regular/decaying",
        "complete ordinary de Rham operator with frozen transmittal/conic domain",
        "equal-volume round--Regge path",
        "standard Mellin asymptotic expansion",
    ],
    "asymptotic_formula": (
        "(4*pi)^(-3/2)[Lambda^3*C0*A0+Lambda*C2*A2+...]"
    ),
    "C0": "4/sqrt(pi) integral_0^infinity chi(v)*v^2 dv > 0",
    "C2": "2/sqrt(pi) integral_0^infinity chi(v) dv > 0",
    "gaussian_family": {
        "chi_a": "exp(-a*v^2), a>0",
        "C0": "a^(-3/2)",
        "C2": "a^(-1/2)",
        "C2_over_C0": "a",
    },
    "derived": [
        "standard positive cutoff cannot reverse the asymptotic A2 sign",
        "at each fixed u<1 the leading shape term prefers round",
    ],
    "open": [
        "selection of chi and Lambda",
        "uniform finite-cutoff remainder and higher singular coefficients",
        "complete-action minimum on the whole path",
        "Newton/Planck normalization and Lorentzian dynamics",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(result["status"])
raise SystemExit(0 if passed == tests else 1)

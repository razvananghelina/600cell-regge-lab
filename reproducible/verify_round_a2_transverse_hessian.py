#!/usr/bin/env python3
"""Exact transverse conformal-Hessian audit for ordinary de Rham A2.

Protocol commit: b366a72.  The calculation is target-free and distinguishes
the smooth full-metric statement from the fixed 600-cell Regge parameter
space.
"""

import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "round_a2_transverse_hessian.json"
PROTOCOL_COMMIT = "b366a72"

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
print("ROUND S3 ORDINARY-DE RHAM A2: TRANSVERSE CONFORMAL HESSIAN")
print("=" * 78)

# Freeze the previously certified positive homogeneous direction without
# re-running that calculation.  Only its exact machine-readable output is an
# input to the new sign/inertia decision.
hopf = json.loads((HERE / "hopf_kahler_induced_gravity.json").read_text())
hopf_exact = hopf["exact_results"]
check(
    "the committed homogeneous Hopf certificate supplies a positive direction",
    hopf_exact["normalized_heat_a2_hessian_on_stf"] == "8/3 I_5"
    and hopf_exact["heat_hessian_inertia"] == [5, 0, 0]
    and hopf_exact["ordinary_derham_heat_a2_curvature_multiplier"] == "-2/3",
)

# Exact mean-zero conformal expansion.  M is kept until after differentiation
# so the absence of a first variation is checked rather than silently built in.
epsilon = sp.symbols("epsilon", real=True)
V, M, F, lam = sp.symbols("V M F lambda", positive=True)
M = sp.symbols("M", real=True)
gradient_energy = lam * F

numerator = 8 * epsilon**2 * gradient_energy + 6 * (
    V + 2 * epsilon * M + epsilon**2 * F
)
volume_polynomial = V + 6 * epsilon * M + 15 * epsilon**2 * F
yamabe = numerator / volume_polynomial ** sp.Rational(1, 3)

yamabe_at_zero = sp.simplify(yamabe.subs(epsilon, 0))
yamabe_first = sp.simplify(sp.diff(yamabe, epsilon).subs(epsilon, 0))
yamabe_second_mean_zero = sp.simplify(
    sp.diff(yamabe, epsilon, 2).subs({epsilon: 0, M: 0})
)
expected_yamabe_second = 16 * (lam - 3) * F / V ** sp.Rational(1, 3)

check(
    "the normalized Einstein functional has the exact round value",
    yamabe_at_zero == 6 * V ** sp.Rational(2, 3),
)
check(
    "the mean-zero conformal path is stationary at round",
    sp.simplify(yamabe_first.subs(M, 0)) == 0,
)
check(
    "the independently expanded normalized Einstein Hessian is exact",
    sp.simplify(yamabe_second_mean_zero - expected_yamabe_second) == 0,
    "delta^2 Y=16(lambda-3)V^(-1/3)F",
)

# The ordinary full-exterior de Rham multiplier is -2/3.  Use it only after
# the geometric expansion has been derived.
ordinary_multiplier = -sp.Rational(2, 3)
a2_hessian = sp.simplify(ordinary_multiplier * yamabe_second_mean_zero)
expected_a2_hessian = (
    -sp.Rational(32, 3) * (lam - 3) * F / V ** sp.Rational(1, 3)
)
check(
    "the ordinary de Rham A2 conformal Hessian follows with no fitted sign",
    sp.simplify(a2_hessian - expected_a2_hessian) == 0,
    "delta^2 A2_hat=-(32/3)(lambda-3)V^(-1/3)F",
)

# Independent explicit l=2 control on the unit S3 subset of R4.
x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4", real=True)
coordinates = (x1, x2, x3, x4)
f = x1**2 - x2**2
euclidean_laplacian = sp.simplify(sum(sp.diff(f, x, 2) for x in coordinates))
euler_degree = sp.simplify(sum(x * sp.diff(f, x) for x in coordinates))
check(
    "f=x1^2-x2^2 is a homogeneous harmonic polynomial of degree two",
    euclidean_laplacian == 0 and euler_degree == 2 * f,
)

# Exact rotation-invariant S3 moments in four ambient dimensions.
# E[x_i^4]=3/(n(n+2)), E[x_i^2 x_j^2]=1/(n(n+2)).
ambient_dimension = sp.Integer(4)
mean_x4 = sp.Rational(3, 1) / (
    ambient_dimension * (ambient_dimension + 2)
)
mean_x2y2 = sp.Rational(1, 1) / (
    ambient_dimension * (ambient_dimension + 2)
)
mean_f2 = sp.simplify(2 * mean_x4 - 2 * mean_x2y2)
check(
    "the explicit l=2 harmonic has zero mean by coordinate symmetry",
    f.xreplace({x1: x2, x2: x1}) == -f,
)
check(
    "its exact S3 squared norm is integral(f^2)=V/6",
    mean_f2 == sp.Rational(1, 6),
)

# Direct tangential-gradient average, independent of inserting the desired
# eigenvalue: |grad_S f|^2=|grad_R4 f|^2-(x.grad_R4 f)^2 on the unit sphere.
mean_x1x2 = sp.Rational(2, ambient_dimension)
mean_tangent_gradient_sq = sp.simplify(4 * mean_x1x2 - 4 * mean_f2)
lambda_explicit = sp.simplify(mean_tangent_gradient_sq / mean_f2)
check(
    "direct tangential-gradient moments give the S3 eigenvalue lambda=8",
    mean_tangent_gradient_sq == sp.Rational(4, 3) and lambda_explicit == 8,
)
check(
    "lambda=8 agrees with the exact scalar-harmonic law l(l+2) at l=2",
    lambda_explicit == 2 * (2 + 2),
)

# Evaluate all frozen sectors.  F>0 and V>0 remain symbolic, so SymPy can
# decide signs without floating-point tolerance.
hessian_l1 = sp.simplify(a2_hessian.subs(lam, 3))
hessian_l2 = sp.simplify(a2_hessian.subs(lam, lambda_explicit))
check(
    "the l=1 conformal sector is an exact zero mode",
    hessian_l1 == 0,
    "this is the Möbius/diffeomorphism conformal sector",
)
check(
    "the first non-gauge l=2 conformal A2 Hessian is strictly negative",
    hessian_l2.is_negative,
    f"delta^2 A2_hat={hessian_l2}",
)

# A direct gauge discriminator avoids relying only on the harmonic label.
# For g'=u^4 g, R'=u^-5(8 Delta u+6u) in the positive-Laplacian convention,
# so delta R=8(lambda-3)f.  Diffeomorphisms of constant R give delta R=0;
# scale gives a constant delta R.  The l=2 variation is therefore neither.
delta_scalar_l2_factor = sp.simplify(8 * (lambda_explicit - 3))
check(
    "the l=2 conformal direction is neither round diffeomorphism nor scale",
    delta_scalar_l2_factor == 40,
    "delta R=40 f is nonzero and nonconstant, whereas L_X(6)=0",
)

# An exactly volume-preserving representative exists by global rescaling.
# If I(epsilon)=integral u^6, multiplying the metric by
# c=(V/I)^(2/3) makes c^(3/2) I=V.  The normalized functional is scale
# invariant, hence has the same Hessian on that equal-volume representative.
I = sp.symbols("I", positive=True)
scale_factor = (V / I) ** sp.Rational(2, 3)
check(
    "global rescaling gives an exact equal-volume representative",
    sp.simplify(scale_factor ** sp.Rational(3, 2) * I - V) == 0,
)

# Formal hostile sign control: the opposite curvature multiplier would make
# the l=2 direction positive, showing where the result enters.  It is not the
# ordinary de Rham coefficient.
control_hessian_l2 = sp.simplify(
    yamabe_second_mean_zero.subs(lam, lambda_explicit)
)
check(
    "the formal +1 curvature-multiplier control reverses the l=2 sign",
    control_hessian_l2.is_positive,
)

# Scope/provenance guards.
protocol = (
    HERE.parent / "docs" / "research" / "round_a2_transverse_hessian_protocol.md"
).read_text()
check(
    "the preregistered protocol forbids projection onto fixed Regge edge space",
    "smooth unit round three-sphere" in protocol
    and "discretization/refinement argument" in protocol
    and "finite fixed 600-cell edge-length space" in protocol,
)
check(
    "the decision boundary was frozen before this result",
    "DERIVED SMOOTH SADDLE" in protocol
    and "REFUTED HOSTILE PREDICTION" in protocol
    and "NO TRANSVERSE-HESSIAN RESULT" in protocol,
)

opposite_signs = (
    hopf_exact["heat_hessian_inertia"] == [5, 0, 0]
    and hessian_l2.is_negative
)
verdict = (
    "DERIVED SMOOTH SADDLE"
    if opposite_signs
    else "REFUTED HOSTILE PREDICTION"
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "provenance": "preregistered hostile transverse-Hessian audit",
    "physical_target_used": False,
    "scope": {
        "proved": "smooth metric space at unit round S3",
        "not_proved": "fixed finite 600-cell Regge edge-length space",
        "operator": "ordinary full-exterior de Rham D=d+d*",
        "functional": "A2_hat=Vol^(-1/3)*[-(2/3) integral R dVol]",
    },
    "exact_results": {
        "yamabe_hessian": "16(lambda-3)V^(-1/3)F",
        "ordinary_derham_A2_hessian": "-(32/3)(lambda-3)V^(-1/3)F",
        "l1_hessian": "0",
        "l2_lambda": 8,
        "l2_norm": "F=V/6",
        "l2_A2_hessian": "-(160/3)V^(-1/3)F = -(80/9)V^(2/3)",
        "homogeneous_hopf_hessian": "positive, inertia (5,0,0)",
        "smooth_full_hessian": "indefinite",
    },
    "verdict": verdict,
    "derived": [
        "round is a saddle of ordinary de Rham A2_hat on smooth metric space",
        "positive homogeneous Hopf and negative non-gauge conformal directions coexist",
        "A2 alone cannot be promoted to a complete smooth Euclidean vacuum selector",
    ],
    "open": [
        "transverse Hessian in the finite 600-cell Regge edge-length space",
        "higher spectral terms and complete finite-cutoff action",
        "Lorentzian contour/dynamics and Newton/Planck normalization",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(verdict)
print("SCOPE: smooth full-metric statement; fixed finite Regge transverse gate OPEN")
raise SystemExit(0 if passed == tests else 1)

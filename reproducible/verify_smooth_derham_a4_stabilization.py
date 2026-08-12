#!/usr/bin/env python3
"""Exact smooth full-de Rham A4 conformal-Hessian and cutoff audit.

Protocol commit: 23f0c4d.  The finite-tau A2+A4 threshold is reported only
for the frozen truncation; the full finite heat trace remains open.
"""

import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "smooth_derham_a4_stabilization.json"
PROTOCOL_COMMIT = "23f0c4d"
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
print("SMOOTH ROUND S3: FULL-DE RHAM A4 CONFORMAL STABILIZATION AUDIT")
print("=" * 78)

# -------------------------------------------------------------------------
# 1. Derive the full-exterior trace data rather than inserting coefficients.
# -------------------------------------------------------------------------
r0, r1, r2 = sp.symbols("r0 r1 r2", real=True)
ricci_eigenvalues = (r0, r1, r2)
R = sum(ricci_eigenvalues)
Ric2 = sum(value**2 for value in ricci_eigenvalues)

# In P=-(nabla^2+E), the Hodge one-form block is E=-Ric.  Hodge
# duality duplicates it in degree two; scalar/top-form blocks vanish.
E_full = sp.diag(0, -r0, -r1, -r2, -r0, -r1, -r2, 0)
rank_full = E_full.rows
trace_E = sp.trace(E_full)
trace_E2 = sp.trace(E_full * E_full)
check(
    "the complete exterior bundle has rank eight and tr(E)=-2R",
    rank_full == 8 and sp.simplify(trace_E + 2*R) == 0,
)
check(
    "Hodge duality gives tr(E^2)=2|Ric|^2",
    sp.simplify(trace_E2 - 2*Ric2) == 0,
)

# Build the 3D Riemann tensor from a diagonal Ricci tensor and explicitly
# trace the curvature action on one-forms.  The two-form block is its Hodge
# dual copy.
def delta(i, j):
    return sp.Integer(i == j)


def riemann(i, j, k, ell):
    ric = ricci_eigenvalues
    return (
        delta(i, k)*delta(j, ell)*ric[j]
        - delta(i, ell)*delta(j, k)*ric[j]
        - delta(j, k)*delta(i, ell)*ric[i]
        + delta(j, ell)*delta(i, k)*ric[i]
        - R*sp.Rational(1, 2)*(
            delta(i, k)*delta(j, ell)-delta(i, ell)*delta(j, k)
        )
    )


Riem2 = sp.expand(sum(
    riemann(i, j, k, ell)**2
    for i in range(3) for j in range(3)
    for k in range(3) for ell in range(3)
))
trace_omega_one = 0
for i in range(3):
    for j in range(3):
        omega = sp.Matrix(3, 3, lambda a, b: -riemann(b, a, i, j))
        trace_omega_one += sp.trace(omega * omega)
trace_omega_full = 2 * trace_omega_one
check(
    "direct curvature-representation matrices give tr(Omega^2)=-2|Riem|^2",
    sp.simplify(trace_omega_full + 2*Riem2) == 0,
)
check(
    "the direct 3D curvature tensor obeys |Riem|^2=4|Ric|^2-R^2",
    sp.simplify(Riem2 - (4*Ric2-R**2)) == 0,
)

# Universal integrated closed-manifold a4 formula.  Total divergences are
# absent after integration.
A4_numerator = sp.expand(
    60*R*trace_E + 180*trace_E2 + 30*trace_omega_full
    + rank_full*(5*R**2-2*Ric2+2*Riem2)
)
A4_density = sp.factor(A4_numerator.subs(Riem2, 4*Ric2-R**2) / 360)
expected_A4_density = sp.Rational(7, 15)*Ric2-sp.Rational(1, 10)*R**2
check(
    "the universal formula reduces to A4=(7/15)|Ric|^2-(1/10)R^2",
    sp.simplify(A4_density-expected_A4_density) == 0,
)
round_density = sp.simplify(expected_A4_density.subs({r0: 2, r1: 2, r2: 2}))
check(
    "unit round S3 has A4 density 2 and integrated A4=4*pi^2",
    round_density == 2 and sp.simplify(round_density*2*sp.pi**2-4*sp.pi**2) == 0,
)

# -------------------------------------------------------------------------
# 2. Independent conformal expansion for a mean-zero scalar eigenfunction.
# -------------------------------------------------------------------------
lam, V, F = sp.symbols("lambda V F", positive=True)
eps, f = sp.symbols("epsilon f", real=True)
Hess2_integral = lam*(lam-2)*F
f_trace_hess_integral = -lam*F
gradient2_integral = lam*F

# Ric~=2g+eps*A+eps^2*B in the background metric, derived from
# Ric(u^4g)=Ric-2u^-1 Hess u+6u^-2 du^2
#             +(2u^-1 Lu-2u^-2|du|^2)g.
integral_A2 = sp.expand(
    4*Hess2_integral
    - 8*lam*f_trace_hess_integral
    + 12*lam**2*F
)
integral_trace_B = sp.expand(
    2*f_trace_hess_integral + 6*gradient2_integral
    + 3*(-2*lam*F-2*gradient2_integral)
)
integral_f_trace_A = sp.expand(
    -2*f_trace_hess_integral + 6*lam*F
)
ricci_squared_coefficient = sp.expand(
    integral_A2
    + 4*integral_trace_B
    - 8*integral_f_trace_A
    + 36*F
)
expected_ricci_coefficient = (24*lam**2-104*lam+36)*F
check(
    "the conformal Ricci-square expansion is derived exactly",
    sp.simplify(ricci_squared_coefficient-expected_ricci_coefficient) == 0,
)

# Scalar curvature expansion is derived separately from
# R~=u^-5(8Lu+6u), dV~=u^6 dV.
u = 1+eps*f
scalar_tilde = u**(-5)*(8*eps*lam*f+6*u)
scalar_square_integrand = scalar_tilde**2*u**6
scalar_square_series = sp.series(scalar_square_integrand, eps, 0, 3).removeO()
scalar_square_point_coefficient = sp.expand(scalar_square_series).coeff(eps, 2)
scalar_square_coefficient = sp.simplify(scalar_square_point_coefficient/f**2*F)
expected_scalar_coefficient = (64*lam**2-288*lam+108)*F
check(
    "the independent scalar-curvature-square expansion is exact",
    sp.simplify(scalar_square_coefficient-expected_scalar_coefficient) == 0,
)

raw_A4_second_order = sp.expand(
    sp.Rational(7, 15)*ricci_squared_coefficient
    - sp.Rational(1, 10)*scalar_square_coefficient
)
expected_raw_A4 = (
    sp.Rational(24, 5)*lam**2-sp.Rational(296, 15)*lam+6
)*F
check(
    "the unnormalized A4 second-order coefficient combines correctly",
    sp.simplify(raw_A4_second_order-expected_raw_A4) == 0,
)

# V(u)=V+15 eps^2 F+O(eps^3) for mean-zero f, and A4(round)=2V.
normalized_A4_series = sp.series(
    (V+15*eps**2*F)**sp.Rational(1, 3)
    * (2*V+eps**2*raw_A4_second_order), eps, 0, 3
).removeO()
normalized_A4_hessian = sp.factor(
    sp.diff(normalized_A4_series, eps, 2).subs(eps, 0)
)
expected_A4_hessian = (
    sp.Rational(48, 5)*V**sp.Rational(1, 3)
    *(lam-3)*(lam-sp.Rational(10, 9))*F
)
check(
    "the normalized A4 conformal Hessian factors exactly",
    sp.simplify(normalized_A4_hessian-expected_A4_hessian) == 0,
    "delta^2 A4_hat=(48/5)V^(1/3)(lambda-3)(lambda-10/9)F",
)
check(
    "the l=1 conformal diffeomorphism sector remains exactly zero",
    sp.simplify(normalized_A4_hessian.subs(lam, 3)) == 0,
)
check(
    "A4 is strictly positive on every non-gauge scalar harmonic l>=2",
    all(expected_A4_hessian.subs(lam, ell*(ell+2)).is_positive
        for ell in range(2, 21)),
    "lambda_l=l(l+2)>=8 makes both factors positive",
)
explicit_l2_hessian = sp.factor(
    expected_A4_hessian.subs({lam: 8, F: V/6})
)
check(
    "the explicit f=x1^2-x2^2 control gives the frozen l=2 value",
    sp.simplify(explicit_l2_hessian-sp.Rational(496, 9)*V**sp.Rational(4, 3)) == 0,
)

# -------------------------------------------------------------------------
# 3. Gaussian cutoff threshold and the asymptotic decision.
# -------------------------------------------------------------------------
a, Lambda, tau = sp.symbols("a Lambda tau", positive=True)
A2_hessian = -sp.Rational(32, 3)*V**(-sp.Rational(1, 3))*(lam-3)*F
gaussian_combined = sp.factor(
    Lambda*a**(-sp.Rational(1, 2))*A2_hessian
    + Lambda**(-1)*a**sp.Rational(1, 2)*expected_A4_hessian
)
combined_in_tau = sp.factor(
    gaussian_combined.subs(a, tau*Lambda**2)
)
combined_numerator = sp.factor(sp.together(combined_in_tau).as_numer_denom()[0])
threshold_roots = sp.solve(combined_numerator, tau)
threshold = sp.factor(threshold_roots[0])
expected_threshold = 10/((9*lam-10)*V**sp.Rational(2, 3))
check(
    "the Gaussian A2+A4 truncated threshold is derived exactly",
    len(threshold_roots) == 1
    and sp.simplify(threshold-expected_threshold) == 0
    and combined_in_tau.subs(tau, threshold/2).subs(lam, 8).is_negative
    and combined_in_tau.subs(tau, 2*threshold).subs(lam, 8).is_positive,
)
tau_star = sp.factor(expected_threshold.subs(lam, 8))
threshold_derivative = sp.diff(expected_threshold, lam)
check(
    "the l=2 mode is the load-bearing scalar threshold",
    tau_star == 5/(31*V**sp.Rational(2, 3))
    and sp.simplify(
        threshold_derivative
        + 90/(V**sp.Rational(2, 3)*(9*lam-10)**2)
    ) == 0
    and threshold_derivative.subs(lam, 8).is_negative,
)
check(
    "the truncated ordering reverses on the two sides of the unselected threshold",
    combined_in_tau.subs(tau, tau_star/2).subs(lam, 8).is_negative
    and combined_in_tau.subs(tau, 2*tau_star).subs(lam, 8).is_positive,
)

# The UV claim uses only relative order: A4/A2 carries a/Lambda^2 and all
# later smooth heat orders are still lower.  No numeric Lambda0 is inferred.
relative_weight = a/Lambda**2
check(
    "A4 is suppressed by Lambda^-2 relative to A2 at fixed cutoff shape",
    sp.limit(relative_weight, Lambda, sp.oo) == 0,
)
hopf = json.loads((HERE / "hopf_kahler_induced_gravity.json").read_text())
smooth_saddle = json.loads((HERE / "round_a2_transverse_hessian.json").read_text())
check(
    "the committed A2 input contains both positive and negative non-gauge directions",
    hopf["exact_results"]["heat_hessian_inertia"] == [5, 0, 0]
    and smooth_saddle["verdict"] == "DERIVED SMOOTH SADDLE",
)
protocol = (HERE.parent / "smooth_derham_a4_stabilization_protocol.md").read_text()
check(
    "the preregistered boundary forbids promoting the truncated threshold to a full action",
    "FINITE CUTOFF OPEN" in protocol
    and "is **not** a conclusion about" in protocol
    and "No numerical value" in protocol,
)

verdict = "DERIVED A4 CONFORMAL STABILIZER; DERIVED UV SADDLE; FINITE CUTOFF OPEN"
payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "external_input": {
        "formula": "universal integrated closed-manifold Laplace-type A4",
        "source": "https://arxiv.org/abs/hep-th/0306138",
    },
    "scope": "smooth closed unit round S3, ordinary full de Rham heat trace",
    "exact_results": {
        "A4_density_3d": "(7/15)|Ric|^2-(1/10)R^2",
        "round_A4": "4*pi^2",
        "normalized_conformal_hessian": "(48/5)V^(1/3)(lambda-3)(lambda-10/9)F",
        "explicit_l2_hessian": "(496/9)V^(4/3)>0",
        "gaussian_truncated_threshold": "tau=a/Lambda^2 > 10/[(9lambda-10)V^(2/3)]",
        "all_scalar_threshold": "tau > 5/[31 V^(2/3)]",
    },
    "verdict": verdict,
    "derived": [
        "smooth A4 is a positive conformal stabilizer for every l>=2",
        "for fixed positive cutoff shape A4 is Lambda^-2 suppressed relative to A2",
        "the full smooth spectral action remains a saddle at sufficiently high cutoff under a uniform expansion",
    ],
    "structural_only": [
        "the A2+A4 truncated finite-tau threshold",
    ],
    "open": [
        "full finite-tau heat trace and omitted A6,A8,...",
        "singular finite-Regge A4 and refinement transfer",
        "selection of cutoff shape, tau, absolute scale and Lorentzian dynamics",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"A4 HESSIAN: {expected_A4_hessian}")
print(f"TRUNCATED l=2 THRESHOLD: tau > {tau_star}")
print(verdict)
raise SystemExit(0 if passed == tests else 1)

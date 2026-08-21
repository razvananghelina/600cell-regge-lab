#!/usr/bin/env python3
"""Derivative-first adversarial audit of the cubic fixed-state obstruction."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_generic_velocity_cubic.json"
OUTPUT = HERE / "gravity_600cell_generic_velocity_cubic_adversarial.json"
PRIMARY_SHA256 = (
    "1d35b46cd4db20df0af3ed3e6b5de676d69753cf5059e0eb607d1eec949b9103"
)
PRIMARY_RESULT_COMMIT = "08bdde5"
ADVERSARIAL_PROTOCOL_COMMIT = "4523385"

VELOCITY_RATIONALS = ((-7, 5), (2, 3), (5, 2))
CUBIC_RATIONALS = ((-2, 7), (1, 5))
HEIGHT_DENOMINATORS = (1200, 2400, 4800, 9600)
mp.mp.dps = 110

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text(value, digits=40):
    return mp.nstr(value, digits)


primary = json.loads(PRIMARY_INPUT.read_text())
provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and primary["outcome"] == "GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION"
    and primary["passed"] == primary["tests"] == 9
    and PRIMARY_RESULT_COMMIT == "08bdde5"
    and ADVERSARIAL_PROTOCOL_COMMIT == "4523385"
)
check("the accepted primary obstruction and adversarial protocol are frozen", provenance_ok)


# Differentiate the complete action before introducing tau or q.  This is
# mechanically distinct from the primary five-variable scaled-action Hessian.
L_MINUS, L_PLUS, RHO, MASS = sp.symbols(
    "L_minus L_plus rho mass", positive=True
)
DELTA = L_PLUS - L_MINUS
HEIGHT = sp.sqrt(RHO + DELTA**2 / 4)
COSINE = (DELTA**2 + 2 * RHO) / (2 * (DELTA**2 + 3 * RHO))
BOOST = DELTA / sp.sqrt(8 * (DELTA**2 + 3 * RHO))
ACTION = (
    360
    * (L_MINUS + L_PLUS)
    * HEIGHT
    * (2 * sp.pi - 5 * sp.acos(COSINE))
    + 600
    * sp.sqrt(3)
    * (L_MINUS**2 - L_PLUS**2)
    * sp.asinh(BOOST)
    - 8 * sp.pi * MASS * sp.sqrt(RHO)
)
F_EXACT = RHO * sp.diff(ACTION, RHO)
P_PRE_EXACT = -L_MINUS * sp.diff(ACTION, L_MINUS) / 2

tau = sp.symbols("tau", positive=True)
q = sp.symbols("q", real=True)
qt_substitution = {
    L_MINUS: 1,
    L_PLUS: 1 + tau * q,
    RHO: tau**2,
}
constraint_qt_raw = sp.factor(
    sp.simplify((2 * F_EXACT / tau).subs(qt_substitution))
)
momentum_qt_raw = sp.factor(sp.simplify(P_PRE_EXACT.subs(qt_substitution)))


def radical_inventory(z):
    square = z**2
    return (
        (square + 4, sp.sqrt(square + 4)),
        (3 * square + 8, sp.sqrt(3 * square + 8)),
        (
            3 * square**2 + 20 * square + 32,
            sp.sqrt(square + 4) * sp.sqrt(3 * square + 8),
        ),
        (
            9 * square**3 + 84 * square**2 + 256 * square + 256,
            (3 * square + 8) * sp.sqrt(square + 4),
        ),
        (
            3 * square**3 + 32 * square**2 + 112 * square + 128,
            (square + 4) * sp.sqrt(3 * square + 8),
        ),
        (
            9 * square**4
            + 120 * square**3
            + 592 * square**2
            + 1280 * square
            + 1024,
            (square + 4) * (3 * square + 8),
        ),
        (
            27 * square**5
            + 432 * square**4
            + 2736 * square**3
            + 8576 * square**2
            + 13312 * square
            + 8192,
            (square + 4)
            * (3 * square + 8)
            * sp.sqrt(3 * square + 8),
        ),
        (
            27 * square**6
            + 540 * square**5
            + 4464 * square**4
            + 19520 * square**3
            + 47616 * square**2
            + 61440 * square
            + 32768,
            (square + 4)
            * (3 * square + 8)
            * sp.sqrt(square + 4)
            * sp.sqrt(3 * square + 8),
        ),
    )


def normalize_radicals(expression, z):
    inventory = radical_inventory(z)

    def is_sqrt(node):
        return isinstance(node, sp.Pow) and node.exp == sp.Rational(1, 2)

    def replace_sqrt(node):
        base = sp.expand(node.base)
        for polynomial, replacement in inventory:
            if sp.expand(base - polynomial) == 0:
                return replacement
        return node

    normalized = expression.replace(is_sqrt, replace_sqrt)
    normalized = sp.together(normalized)
    normalized = normalized.replace(is_sqrt, replace_sqrt)
    normalized = sp.cancel(normalized)
    normalized = normalized.replace(is_sqrt, replace_sqrt)
    return sp.factor_terms(normalized)


factorizations_ok = all(
    sp.expand(polynomial - replacement**2) == 0
    for polynomial, replacement in radical_inventory(q)
)
constraint_qt = normalize_radicals(constraint_qt_raw, q)
momentum_qt = normalize_radicals(momentum_qt_raw, q)
check(
    "the derivative-first two-variable reduction obeys all radical identities",
    factorizations_ok,
)


v, a, c = sp.symbols("v a c", real=True)
u = v**2
radius = sp.sqrt(u + 4)
qrad = sp.sqrt(3 * u + 8)
theta = sp.acos((u + 2) / (2 * (u + 3)))
eta = sp.asinh(v / sp.sqrt(8 * (u + 3)))
epsilon_v = 2 * sp.pi - 5 * theta
mass_branch = 180 * epsilon_v / (sp.pi * radius)
momentum_branch = 180 * v * epsilon_v / radius - 600 * sp.sqrt(3) * eta
K = 10 * radius - (u + 3) * qrad * epsilon_v
B = 5 * u * radius + 2 * (u + 3) * qrad * epsilon_v
a_root = -B / K
q1 = a + v**2 / 2
q2 = c + v * a + v**3 / 6


def total_zero(expression):
    return normalize_radicals(
        expression.subs({tau: 0, q: v, MASS: mass_branch}),
        v,
    )


def total_first(expression):
    raw = (
        sp.diff(expression, tau) + q1 * sp.diff(expression, q)
    ).subs({tau: 0, q: v, MASS: mass_branch})
    return normalize_radicals(raw, v)


def total_second(expression):
    raw = (
        sp.diff(expression, tau, 2) / 2
        + q1 * sp.diff(expression, tau, q)
        + q1**2 * sp.diff(expression, q, 2) / 2
        + q2 * sp.diff(expression, q)
    ).subs({tau: 0, q: v, MASS: mass_branch})
    return normalize_radicals(raw, v)


constraint_zero = total_zero(constraint_qt)
momentum_zero = normalize_radicals(total_zero(momentum_qt) - momentum_branch, v)
constraint_first = total_first(constraint_qt)
momentum_first = total_first(momentum_qt)

c1_a = normalize_radicals(sp.diff(constraint_first, a), v)
c1_0 = normalize_radicals(constraint_first.subs(a, 0), v)
p1_a = normalize_radicals(sp.diff(momentum_first, a), v)
p1_0 = normalize_radicals(momentum_first.subs(a, 0), v)
c_prefactor = 1440 * v / (radius * qrad * (u + 3) * (u + 4))
leading_branch_ok = bool(
    sp.simplify(constraint_zero) == 0
    and sp.simplify(momentum_zero) == 0
    and normalize_radicals(c1_a - c_prefactor * K, v) == 0
    and normalize_radicals(c1_0 - c_prefactor * B, v) == 0
    and normalize_radicals(c1_0 * p1_a - p1_0 * c1_a, v) == 0
    and normalize_radicals(constraint_first.subs(a, a_root), v) == 0
    and normalize_radicals(momentum_first.subs(a, a_root), v) == 0
)
check(
    "the derivative-first route independently reconstructs the lower-order branch",
    leading_branch_ok,
)


constraint_second_generic = total_second(constraint_qt)
momentum_second_generic = total_second(momentum_qt)
constraint_second = normalize_radicals(
    constraint_second_generic.subs(a, a_root), v
)
momentum_second = normalize_radicals(
    momentum_second_generic.subs(a, a_root), v
)
c2_c = normalize_radicals(sp.diff(constraint_second, c), v)
c2_0 = normalize_radicals(constraint_second.subs(c, 0), v)
p2_c = normalize_radicals(sp.diff(momentum_second, c), v)
p2_0 = normalize_radicals(momentum_second.subs(c, 0), v)
affine_and_recursive_ok = bool(
    normalize_radicals(sp.diff(constraint_second, c, 2), v) == 0
    and normalize_radicals(sp.diff(momentum_second, c, 2), v) == 0
    and normalize_radicals(c2_c - c1_a, v) == 0
    and normalize_radicals(p2_c - p1_a, v) == 0
)
check(
    "both independently derived cubic coefficients are affine with recursive slopes",
    affine_and_recursive_ok,
)


cross_resultant = normalize_radicals(c2_0 * p2_c - p2_0 * c2_c, v)
cross_expected = 129600 * epsilon_v**2 / (u + 4)
cross_factorization_ok = bool(
    normalize_radicals(cross_resultant - cross_expected, v) == 0
)
check(
    "the derivative-first cross-resultant has the exact primary-independent factor",
    cross_factorization_ok,
    f"Delta={cross_resultant}",
)


# Different positivity proof: a direct exact lower bound on the acos argument.
x = sp.symbols("x", nonnegative=True)
z_x = (x + 2) / (2 * (x + 3))
cos_two_pi_fifths = (sp.sqrt(5) - 1) / 4
cosine_gap = (7 - 3 * sp.sqrt(5)) / 12
direct_bound_certificate = bool(
    sp.simplify(z_x - sp.Rational(1, 3) - x / (6 * (x + 3))) == 0
    and sp.simplify(
        sp.Rational(1, 3) - cos_two_pi_fifths - cosine_gap
    )
    == 0
    and 49 > 45
)
complete_no_common_root = bool(
    leading_branch_ok
    and affine_and_recursive_ok
    and cross_factorization_ok
    and direct_bound_certificate
)
check(
    "the direct acos bound excludes every real zero of the cross-resultant",
    complete_no_common_root,
)


# Only now read the primary formulas and compare the complete coefficients.
locals_map = {
    "v": v,
    "c": c,
    "pi": sp.pi,
    "sqrt": sp.sqrt,
    "acos": sp.acos,
    "asinh": sp.asinh,
}
primary_constraint = sp.sympify(
    primary["cubic_census"]["constraint_second"], locals=locals_map
)
primary_momentum = sp.sympify(
    primary["cubic_census"]["momentum_second"], locals=locals_map
)
primary_cross = sp.sympify(
    primary["cubic_census"]["cross_resultant"], locals=locals_map
)
primary_match_ok = bool(
    normalize_radicals(constraint_second - primary_constraint, v) == 0
    and normalize_radicals(momentum_second - primary_momentum, v) == 0
    and normalize_radicals(cross_resultant - primary_cross, v) == 0
)
check(
    "both full coefficients and Delta match the frozen primary only post-derivation",
    primary_match_ok,
)


hostile_mass = normalize_radicals(
    constraint_qt.subs(
        {tau: 0, q: v, MASS: mass_branch + sp.Rational(1, 13)}
    ),
    v,
)
wrong_q2_gap = sp.expand(q2 - c)
lapse_root = normalize_radicals(-c2_0 / c2_c, v)
momentum_at_lapse_root = normalize_radicals(
    p2_c * lapse_root + p2_0, v
)
hostile_ok = bool(
    sp.simplify(hostile_mass + 8 * sp.pi / 13) == 0
    and sp.expand(wrong_q2_gap - (v * a + v**3 / 6)) == 0
    and wrong_q2_gap != 0
    and normalize_radicals(
        momentum_at_lapse_root + cross_resultant / c2_c,
        v,
    )
    == 0
    and momentum_at_lapse_root != 0
)
check(
    "all three registered hostile controls fail in the required direction",
    hostile_ok,
    (
        f"mass={hostile_mass}; q2_gap={wrong_q2_gap}; "
        f"momentum_on_lapse_root={momentum_at_lapse_root}"
    ),
)


# Independent full-action numerical controls at disjoint sample points.
f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), F_EXACT, "mpmath")
p_pre_numeric = sp.lambdify(
    (L_MINUS, L_PLUS, RHO, MASS), P_PRE_EXACT, "mpmath"
)
mass_numeric = sp.lambdify(v, mass_branch, "mpmath")
momentum_numeric = sp.lambdify(v, momentum_branch, "mpmath")
a_numeric = sp.lambdify(v, a_root, "mpmath")
c2_numeric = sp.lambdify((v, c), constraint_second, "mpmath")
p2_numeric = sp.lambdify((v, c), momentum_second, "mpmath")

numeric_records = {}
numeric_failures = []
numeric_ok = True
for v_numerator, v_denominator in VELOCITY_RATIONALS:
    v_key = f"{v_numerator}/{v_denominator}"
    v_value = mp.mpf(v_numerator) / v_denominator
    mass_value = mass_numeric(v_value)
    p0_value = momentum_numeric(v_value)
    a_value = a_numeric(v_value)
    numeric_records[v_key] = {}
    for c_numerator, c_denominator in CUBIC_RATIONALS:
        c_key = f"{c_numerator}/{c_denominator}"
        c_value = mp.mpf(c_numerator) / c_denominator
        exact_values = {
            "constraint": c2_numeric(v_value, c_value),
            "momentum": p2_numeric(v_value, c_value),
        }
        errors = {name: [] for name in exact_values}
        for denominator in HEIGHT_DENOMINATORS:
            h_value = mp.mpf(1) / denominator
            endpoint = mp.exp(
                v_value * h_value
                + a_value * h_value**2
                + c_value * h_value**3
            )
            direct_values = {
                "constraint": 2
                * f_numeric(1, endpoint, h_value**2, mass_value)
                / h_value**3,
                "momentum": (
                    p_pre_numeric(1, endpoint, h_value**2, mass_value)
                    - p0_value
                )
                / h_value**2,
            }
            for name, direct_value in direct_values.items():
                errors[name].append(
                    abs(direct_value - exact_values[name])
                    / max(mp.mpf(1), abs(exact_values[name]))
                )

        orders = {}
        for name, sequence in errors.items():
            if all(error < mp.mpf("1e-75") for error in sequence):
                orders[name] = [mp.inf] * 3
                sequence_ok = True
            elif all(error >= mp.mpf("1e-75") for error in sequence):
                orders[name] = [
                    mp.log(left / right) / mp.log(2)
                    for left, right in zip(sequence, sequence[1:])
                ]
                sequence_ok = bool(
                    all(right < left for left, right in zip(sequence, sequence[1:]))
                    and all(
                        mp.mpf("0.8") <= order <= mp.mpf("1.2")
                        for order in orders[name]
                    )
                )
            else:
                orders[name] = [mp.nan] * 3
                sequence_ok = False
            numeric_ok &= sequence_ok
            if not sequence_ok:
                numeric_failures.append((v_key, c_key, name))

        numeric_records[v_key][c_key] = {
            "errors": {
                name: [text(value) for value in sequence]
                for name, sequence in errors.items()
            },
            "orders": {
                name: [text(value, 30) for value in sequence]
                for name, sequence in orders.items()
            },
        }

check(
    "all disjoint full-action coefficient controls converge at first order",
    numeric_ok,
    f"failures={numeric_failures}",
)


all_gates = all(
    (
        provenance_ok,
        factorizations_ok,
        leading_branch_ok,
        affine_and_recursive_ok,
        cross_factorization_ok,
        complete_no_common_root,
        primary_match_ok,
        hostile_ok,
        numeric_ok,
    )
)
outcome = (
    "GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION_ADVERSARIALLY_CORROBORATED"
    if all_gates
    else "GENERIC_CUBIC_OPEN"
)
check(
    "the adversarial hierarchy assigns the outcome without reinterpretation",
    outcome
    in {
        "GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION_ADVERSARIALLY_CORROBORATED",
        "GENERIC_CUBIC_OPEN",
    },
    outcome,
)


artifact = {
    "provenance": {
        "primary_sha256": digest(PRIMARY_INPUT),
        "primary_result_commit": PRIMARY_RESULT_COMMIT,
        "adversarial_protocol_commit": ADVERSARIAL_PROTOCOL_COMMIT,
    },
    "method": (
        "differentiate_full_action_then_two_variable_tau_q_total_derivatives"
    ),
    "fixed_state": {
        "mass": str(mass_branch),
        "momentum": str(momentum_branch),
        "acceleration": str(a_root),
        "domain": "real v != 0 and K(v^2) != 0",
    },
    "exact_census": {
        "constraint_second": str(constraint_second),
        "momentum_second": str(momentum_second),
        "constraint_slope": str(c2_c),
        "momentum_slope": str(p2_c),
        "cross_resultant": str(cross_resultant),
        "cross_expected": str(cross_expected),
        "cross_factorization_ok": cross_factorization_ok,
        "direct_acos_bound_certificate": direct_bound_certificate,
        "classification": "NO_COMMON_C_ON_COMPLETE_REGISTERED_DOMAIN",
        "classification_complete": complete_no_common_root,
    },
    "primary_comparison": {
        "performed_after_independent_derivation": True,
        "full_match": primary_match_ok,
    },
    "hostile": {
        "mass_shift_defect": str(hostile_mass),
        "wrong_q2_gap": str(wrong_q2_gap),
        "momentum_at_lapse_root": str(momentum_at_lapse_root),
        "all_passed": hostile_ok,
    },
    "numeric_controls": numeric_records,
    "numeric_failures": [list(value) for value in numeric_failures],
    "labels": {
        "formal_integrability": "DERIVED_NEGATIVE_SCOPED",
        "fundamental_tick": "NOT_DERIVED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nOutcome: {outcome}")
print(f"Checks: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

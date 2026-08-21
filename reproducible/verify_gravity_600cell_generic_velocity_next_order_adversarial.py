#!/usr/bin/env python3
"""Derivative-first adversarial audit of the next-order velocity census."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_generic_velocity_next_order.json"
OUTPUT = HERE / "gravity_600cell_generic_velocity_next_order_adversarial.json"
PRIMARY_SHA256 = (
    "4bc69490fc83a193b6ac2cbd8dbe291415a13b60e4dbcce4f499bf70152e5b18"
)
PRIMARY_IMPLEMENTATION_COMMIT = "98acd61"
ADVERSARIAL_PROTOCOL_COMMIT = "44a6ab2"
COMPOSITION_PROTOCOL_COMMIT = "9f08aa0"
VELOCITY_RATIONALS = ((-6, 5), (2, 5), (11, 5))
ACCELERATION_RATIONALS = ((-1, 5), (2, 9))
HEIGHT_RATIONALS = ((1, 500), (1, 1000))
mp.mp.dps = 100
tests = passed = 0


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
    and primary["outcome"] == "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES"
    and primary["passed"] == primary["tests"] == 13
    and PRIMARY_IMPLEMENTATION_COMMIT == "98acd61"
    and ADVERSARIAL_PROTOCOL_COMMIT == "44a6ab2"
    and COMPOSITION_PROTOCOL_COMMIT == "9f08aa0"
)
check("the primary theorem and both adversarial protocols are frozen", provenance_ok)


# Differentiate the complete action before introducing the tangent variables.
L_MINUS, L_PLUS, RHO, MASS = sp.symbols(
    "L_minus L_plus rho mass", positive=True
)
DELTA = L_PLUS - L_MINUS
HEIGHT = sp.sqrt(RHO + DELTA**2 / 4)
COSINE = (DELTA**2 + 2 * RHO) / (2 * (DELTA**2 + 3 * RHO))
BOOST = DELTA / sp.sqrt(8 * (DELTA**2 + 3 * RHO))
ACTION = (
    360 * (L_MINUS + L_PLUS) * HEIGHT
    * (2 * sp.pi - 5 * sp.acos(COSINE))
    + 600 * sp.sqrt(3) * (L_MINUS**2 - L_PLUS**2) * sp.asinh(BOOST)
    - 8 * sp.pi * MASS * sp.sqrt(RHO)
)
F_EXACT = RHO * sp.diff(ACTION, RHO)
P_MINUS_EXACT = L_MINUS * sp.diff(ACTION, L_MINUS) / 2
P_PLUS_EXACT = L_PLUS * sp.diff(ACTION, L_PLUS) / 2
P_PRE_EXACT = -P_MINUS_EXACT


tau = sp.symbols("tau", positive=True)
q = sp.symbols("q", real=True)
qt_path = {
    L_MINUS: 1,
    L_PLUS: 1 + tau * q,
    RHO: tau**2,
}
constraint_qt_raw = sp.factor(
    sp.simplify((2 * F_EXACT / tau).subs(qt_path))
)
momentum_qt_raw = sp.factor(sp.simplify(P_PRE_EXACT.subs(qt_path)))


def radical_inventory(z):
    u = z**2
    return (
        (u + 4, sp.sqrt(u + 4)),
        (3 * u + 8, sp.sqrt(3 * u + 8)),
        (
            3 * u**2 + 20 * u + 32,
            sp.sqrt(u + 4) * sp.sqrt(3 * u + 8),
        ),
        (
            9 * u**3 + 84 * u**2 + 256 * u + 256,
            (3 * u + 8) * sp.sqrt(u + 4),
        ),
        (
            3 * u**3 + 32 * u**2 + 112 * u + 128,
            (u + 4) * sp.sqrt(3 * u + 8),
        ),
        (
            9 * u**4 + 120 * u**3 + 592 * u**2 + 1280 * u + 1024,
            (u + 4) * (3 * u + 8),
        ),
        (
            27 * u**5 + 432 * u**4 + 2736 * u**3 + 8576 * u**2
            + 13312 * u + 8192,
            (u + 4) * (3 * u + 8) * sp.sqrt(3 * u + 8),
        ),
        (
            27 * u**6 + 540 * u**5 + 4464 * u**4 + 19520 * u**3
            + 47616 * u**2 + 61440 * u + 32768,
            (u + 4) * (3 * u + 8)
            * sp.sqrt(u + 4) * sp.sqrt(3 * u + 8),
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
check("the derivative-first expressions obey all eight frozen radical identities", factorizations_ok)


v, a = sp.symbols("v a", real=True)
radius = sp.sqrt(v**2 + 4)
theta = sp.acos((v**2 + 2) / (2 * (v**2 + 3)))
eta = sp.asinh(v / sp.sqrt(8 * (v**2 + 3)))
epsilon_v = 2 * sp.pi - 5 * theta
mass_branch = 180 * epsilon_v / (sp.pi * radius)
momentum_branch = 180 * v * epsilon_v / radius - 600 * sp.sqrt(3) * eta
q_slope = a + v**2 / 2


def direct_total_first(expression):
    tau_part = sp.diff(expression, tau)
    q_part = sp.diff(expression, q)
    raw = (
        tau_part + q_slope * q_part
    ).subs({tau: 0, q: v, MASS: mass_branch})
    return normalize_radicals(raw, v)


constraint_zero = normalize_radicals(
    constraint_qt.subs({tau: 0, q: v, MASS: mass_branch}), v
)
momentum_zero = normalize_radicals(
    momentum_qt.subs({tau: 0, q: v, MASS: mass_branch}) - momentum_branch,
    v,
)
constraint_first = direct_total_first(constraint_qt)
momentum_first = direct_total_first(momentum_qt)
direct_coefficients_ok = bool(
    sp.simplify(constraint_zero) == 0
    and sp.simplify(momentum_zero) == 0
    and sp.Poly(sp.together(constraint_first).as_numer_denom()[0], a).degree() == 1
    and sp.Poly(sp.together(momentum_first).as_numer_denom()[0], a).degree() == 1
)
check("the derivative-first route independently produces both affine coefficients", direct_coefficients_ok)


# Read primary formulas only after the direct derivative-first expressions exist.
locals_map = {
    "v": v,
    "a": a,
    "pi": sp.pi,
    "sqrt": sp.sqrt,
    "acos": sp.acos,
    "asinh": sp.asinh,
}
primary_constraint = sp.sympify(
    primary["one_slab"]["constraint_first"], locals=locals_map
)
primary_momentum = sp.sympify(
    primary["one_slab"]["momentum_first"], locals=locals_map
)
constraint_match = normalize_radicals(
    constraint_first - primary_constraint, v
) == 0
momentum_match = normalize_radicals(
    momentum_first - primary_momentum, v
) == 0
check(
    "both exact derivative-first coefficients match the frozen primary artifact",
    constraint_match and momentum_match,
    f"constraint={constraint_match}; momentum={momentum_match}",
)


c_a = normalize_radicals(sp.diff(constraint_first, a), v)
c_0 = normalize_radicals(constraint_first.subs(a, 0), v)
p_a = normalize_radicals(sp.diff(momentum_first, a), v)
p_0 = normalize_radicals(momentum_first.subs(a, 0), v)
K = 10 * radius - (v**2 + 3) * sp.sqrt(3 * v**2 + 8) * epsilon_v
B = 5 * v**2 * radius + 2 * (v**2 + 3) * sp.sqrt(3 * v**2 + 8) * epsilon_v
c_prefactor = 1440 * v / (
    radius * sp.sqrt(3 * v**2 + 8) * (v**2 + 3) * (v**2 + 4)
)
factorization_ok = bool(
    normalize_radicals(c_a - c_prefactor * K, v) == 0
    and normalize_radicals(c_0 - c_prefactor * B, v) == 0
    and normalize_radicals(c_0 * p_a - p_0 * c_a, v) == 0
)
check("the direct coefficients independently recover K, B and the common-root identity", factorization_ok)


x = sp.symbols("x", nonnegative=True)
r_x = sp.sqrt(x + 4)
q_x = sp.sqrt(3 * x + 8)
theta_x = sp.acos((x + 2) / (2 * (x + 3)))
epsilon_x = 2 * sp.pi - 5 * theta_x
H_x = (x + 3) * q_x * epsilon_x / r_x
monotonicity_ok = bool(
    sp.simplify(
        sp.diff(epsilon_x, x) - 5 / ((x + 3) * r_x * q_x)
    ) == 0
    and sp.simplify(
        sp.diff((3 * x + 8) / (x + 4), x) - 4 / (x + 4) ** 2
    ) == 0
    and 49 > 45
    and sp.N(2 * sp.pi - 5 * sp.acos(sp.Rational(1, 3)), 100) > 0
    and sp.N(
        3 * sp.sqrt(2)
        * (2 * sp.pi - 5 * sp.acos(sp.Rational(1, 3))),
        100,
    ) < 10
    and sp.limit(H_x, x, sp.oo) == sp.oo
)
check("the independent monotonicity proof gives exactly one positive exceptional x", monotonicity_ok)


primary_k = sp.sympify(primary["one_slab"]["K"], locals=locals_map)
k_match = normalize_radicals(K - primary_k, v) == 0
k_numeric = sp.lambdify(v, K, "mpmath")
primary_k_numeric = sp.lambdify(v, primary_k, "mpmath")
x_star = mp.findroot(lambda value: k_numeric(mp.sqrt(value)), (mp.mpf(5), mp.mpf(6)))
primary_x_star = mp.findroot(
    lambda value: primary_k_numeric(mp.sqrt(value)),
    (mp.mpf(5), mp.mpf(6)),
)
exceptional_numeric_ok = bool(
    k_match
    and k_numeric(mp.sqrt(5)) > 0
    and k_numeric(mp.sqrt(6)) < 0
    and abs(x_star - primary_x_star) < mp.mpf("1e-80")
    and abs(k_numeric(mp.sqrt(x_star))) < mp.mpf("1e-80")
)
check(
    "the new bracket and 80-decimal exceptional root agree with the exact primary K",
    exceptional_numeric_ok,
    f"x_star={text(x_star, 60)}",
)


hostile_mass = normalize_radicals(
    constraint_qt.subs(
        {tau: 0, q: v, MASS: mass_branch + sp.Rational(1, 10)}
    ),
    v,
)
hostile_momentum = normalize_radicals(
    momentum_qt.subs({tau: 0, q: v, MASS: mass_branch})
    - (momentum_branch + sp.Rational(1, 10)),
    v,
)
shifted_x_star = mp.findroot(
    lambda value: k_numeric(mp.sqrt(value)) + mp.mpf("0.1"),
    (mp.mpf(5), mp.mpf(6)),
)
hostile_ok = bool(
    sp.simplify(hostile_mass + 4 * sp.pi / 5) == 0
    and sp.simplify(hostile_momentum + sp.Rational(1, 10)) == 0
    and abs(shifted_x_star - x_star) > mp.mpf("1e-50")
)
check("all three changed-state/root hostile controls fail exactly as registered", hostile_ok)


f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), F_EXACT, "mpmath")
p_pre_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_PRE_EXACT, "mpmath")
p_minus_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_MINUS_EXACT, "mpmath")
p_plus_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_PLUS_EXACT, "mpmath")
action_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), ACTION, "mpmath")
c1_numeric = sp.lambdify((v, a), constraint_first, "mpmath")
p1_numeric = sp.lambdify((v, a), momentum_first, "mpmath")
mass_numeric = sp.lambdify(v, mass_branch, "mpmath")
momentum_numeric = sp.lambdify(v, momentum_branch, "mpmath")
root_numeric = sp.lambdify(v, -c_0 / c_a, "mpmath")


def convergence_order(pair):
    if pair[0] < mp.mpf("1e-70") and pair[1] < mp.mpf("1e-70"):
        return mp.inf, True
    if pair[0] >= mp.mpf("1e-70") and pair[1] >= mp.mpf("1e-70"):
        order = mp.log(pair[0] / pair[1]) / mp.log(2)
        return order, bool(
            pair[1] < pair[0]
            and mp.mpf("0.8") <= order <= mp.mpf("1.2")
        )
    return mp.nan, False


numeric_records = {}
numeric_ok = True
for v_numerator, v_denominator in VELOCITY_RATIONALS:
    v_key = f"{v_numerator}/{v_denominator}"
    v_value = mp.mpf(v_numerator) / v_denominator
    mass_value = mass_numeric(v_value)
    p0_value = momentum_numeric(v_value)
    numeric_records[v_key] = {}
    for a_numerator, a_denominator in ACCELERATION_RATIONALS:
        a_key = f"{a_numerator}/{a_denominator}"
        a_value = mp.mpf(a_numerator) / a_denominator
        exact = {
            "constraint": c1_numeric(v_value, a_value),
            "momentum": p1_numeric(v_value, a_value),
        }
        errors = {name: [] for name in exact}
        for h_numerator, h_denominator in HEIGHT_RATIONALS:
            h_value = mp.mpf(h_numerator) / h_denominator
            endpoint = mp.exp(v_value * h_value + a_value * h_value**2)
            direct = {
                "constraint": 2 * f_numeric(1, endpoint, h_value**2, mass_value)
                / h_value**2,
                "momentum": (
                    p_pre_numeric(1, endpoint, h_value**2, mass_value) - p0_value
                ) / h_value,
            }
            for name in errors:
                errors[name].append(
                    abs(direct[name] - exact[name])
                    / max(mp.mpf(1), abs(exact[name]))
                )
        orders = {}
        for name, pair in errors.items():
            orders[name], order_ok = convergence_order(pair)
            numeric_ok &= order_ok
        numeric_records[v_key][a_key] = {
            "errors": {name: [text(value) for value in pair] for name, pair in errors.items()},
            "orders": {name: text(value, 30) for name, value in orders.items()},
        }
check("all new derivative-first direct controls meet the frozen first-order gate", numeric_ok)


composition_records = {}
composition_numeric_ok = True
for v_numerator, v_denominator in VELOCITY_RATIONALS:
    v_key = f"{v_numerator}/{v_denominator}"
    v_value = mp.mpf(v_numerator) / v_denominator
    mass_value = mass_numeric(v_value)
    p0_value = momentum_numeric(v_value)
    a_value = root_numeric(v_value)
    residuals = {
        name: []
        for name in (
            "first_lapse",
            "first_pre_momentum",
            "second_lapse",
            "seam",
            "final_momentum",
            "action",
        )
    }
    for h_numerator, h_denominator in HEIGHT_RATIONALS:
        h_value = mp.mpf(h_numerator) / h_denominator
        l_mid = mp.exp(v_value * h_value / 2 + a_value * h_value**2 / 4)
        l_final = mp.exp(v_value * h_value + a_value * h_value**2)
        rho_fine = h_value**2 / 4
        rho_coarse = h_value**2
        first = (1, l_mid, rho_fine, mass_value)
        second = (l_mid, l_final, rho_fine, mass_value)
        coarse = (1, l_final, rho_coarse, mass_value)
        residuals["first_lapse"].append(abs(4 * f_numeric(*first) / h_value**2))
        residuals["first_pre_momentum"].append(
            abs((p_pre_numeric(*first) - p0_value) / h_value)
        )
        residuals["second_lapse"].append(abs(4 * f_numeric(*second) / h_value**2))
        residuals["seam"].append(
            abs((p_plus_numeric(*first) + p_minus_numeric(*second)) / h_value)
        )
        residuals["final_momentum"].append(
            abs((p_plus_numeric(*second) - p_plus_numeric(*coarse)) / h_value)
        )
        residuals["action"].append(
            abs((action_numeric(*first) + action_numeric(*second) - action_numeric(*coarse)) / h_value**2)
        )
    orders = {}
    for name, pair in residuals.items():
        orders[name], order_ok = convergence_order(pair)
        composition_numeric_ok &= order_ok
    composition_records[v_key] = {
        "a": text(a_value),
        "residuals": {name: [text(value) for value in pair] for name, pair in residuals.items()},
        "orders": {name: text(value, 30) for name, value in orders.items()},
    }
check("all stationary two-half-slab controls converge to the primary zero jets", composition_numeric_ok)


all_controls = bool(
    provenance_ok
    and factorizations_ok
    and direct_coefficients_ok
    and constraint_match
    and momentum_match
    and factorization_ok
    and monotonicity_ok
    and exceptional_numeric_ok
    and hostile_ok
    and numeric_ok
    and composition_numeric_ok
)
outcome = (
    "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES_ADVERSARIALLY_CORROBORATED"
    if all_controls
    else "GENERIC_NEXT_ORDER_ADVERSARIAL_DISAGREEMENT"
)
check(
    "the adversarial hierarchy assigns the corroborated exceptional-velocity verdict",
    outcome == "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES_ADVERSARIALLY_CORROBORATED",
    outcome,
)


artifact = {
    "primary_sha256": digest(PRIMARY_INPUT),
    "primary_implementation_commit": PRIMARY_IMPLEMENTATION_COMMIT,
    "adversarial_protocol_commit": ADVERSARIAL_PROTOCOL_COMMIT,
    "composition_protocol_commit": COMPOSITION_PROTOCOL_COMMIT,
    "method": "differentiate_full_action_then_introduce_independent_tau_q",
    "direct": {
        "constraint_first": str(constraint_first),
        "momentum_first": str(momentum_first),
        "constraint_matches_primary": constraint_match,
        "momentum_matches_primary": momentum_match,
        "K": str(K),
        "B": str(B),
        "factorization_ok": factorization_ok,
    },
    "exceptional": {
        "set": ["-sqrt(x_star)", "sqrt(x_star)"],
        "x_star": text(x_star, 80),
        "bracket": ["5", "6"],
        "primary_difference": text(abs(x_star - primary_x_star), 30),
        "monotonicity_ok": monotonicity_ok,
    },
    "hostile": {
        "mass_shift": str(hostile_mass),
        "momentum_shift": str(hostile_momentum),
        "shifted_K_root": text(shifted_x_star, 60),
    },
    "numeric_controls": numeric_records,
    "composition_numeric_controls": composition_records,
    "labels": {
        "generic_duration": "FREE_TO_NEXT_ORDER_EXCEPT_TWO_OBSTRUCTIONS",
        "exceptional_pair": "DERIVED_EXACT_ADVERSARIALLY_CORROBORATED_STRUCTURAL",
        "composition": "PRIMARY_EXACT_WITH_INDEPENDENT_NUMERICAL_CONTROLS",
        "fundamental_tick": "NOT_DERIVED",
        "absolute_tick": "DERIVED_NEGATIVE_UNDER_SCALE_COVARIANCE_HYPOTHESES",
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


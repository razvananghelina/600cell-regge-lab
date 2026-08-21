#!/usr/bin/env python3
"""Exact-derivative-first adversarial generic-velocity audit."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_generic_velocity_composition.json"
OUTPUT = HERE / "gravity_600cell_generic_velocity_composition_adversarial.json"
PRIMARY_SHA256 = (
    "8ded36f1fa00307fcb23369c25290c9f5bd701709762d6a865437c2507eabfc9"
)
PRIMARY_IMPLEMENTATION_COMMIT = "a5ebcea"
ADVERSARIAL_PROTOCOL_COMMIT = "7ed16ee"
RADICAL_PROTOCOL_COMMIT = "f424f31"
VELOCITY_TEXTS = ("-0.7", "0.3", "1.3")
S_RATIONALS = ((3, 4), (1, 3))
E_RATIONALS = ((1, 300), (1, 600))
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
    and primary["outcome"] == "GENERIC_VELOCITY_LEADING_REPARAMETRIZATION"
    and primary["passed"] == primary["tests"] == 11
    and PRIMARY_IMPLEMENTATION_COMMIT == "a5ebcea"
    and ADVERSARIAL_PROTOCOL_COMMIT == "7ed16ee"
    and RADICAL_PROTOCOL_COMMIT == "f424f31"
)
check("the corrected primary result and adversarial protocols are frozen", provenance_ok)


# Reconstruct and differentiate the unexpanded action before any limit.
L_MINUS, L_PLUS, RHO, MASS = sp.symbols(
    "L_minus L_plus rho mass", positive=True
)
DELTA = L_PLUS - L_MINUS
HEIGHT = sp.sqrt(RHO + DELTA**2 / 4)
COSINE = (DELTA**2 + 2 * RHO) / (2 * (DELTA**2 + 3 * RHO))
BOOST = DELTA / sp.sqrt(8 * (DELTA**2 + 3 * RHO))
ACTION = (
    360 * (L_MINUS + L_PLUS) * HEIGHT * (2 * sp.pi - 5 * sp.acos(COSINE))
    + 600 * sp.sqrt(3) * (L_MINUS**2 - L_PLUS**2) * sp.asinh(BOOST)
    - 8 * sp.pi * MASS * sp.sqrt(RHO)
)
F_EXACT = RHO * sp.diff(ACTION, RHO)
P_PRE_EXACT = -L_MINUS * sp.diff(ACTION, L_MINUS) / 2

tau = sp.symbols("tau", positive=True)
v, mu, s, e = sp.symbols("v mu s e", real=True)
linear_path = {
    L_MINUS: 1,
    L_PLUS: 1 + v * tau,
    RHO: tau**2,
    MASS: mu,
}
direct_action = sp.factor(
    sp.limit(ACTION.subs(linear_path) / tau, tau, 0, dir="+")
)
direct_constraint_raw = sp.factor(
    sp.limit(2 * F_EXACT.subs(linear_path) / tau, tau, 0, dir="+")
)
direct_momentum_raw = sp.factor(
    sp.limit(P_PRE_EXACT.subs(linear_path), tau, 0, dir="+")
)

tangent_ok = bool(
    sp.limit((sp.exp(s * v * e) - 1) / (s * e), e, 0, dir="+") == v
    and sp.simplify((s * e) / (s * e) - 1) == 0
    and not any(expression.has(s) for expression in (
        direct_action, direct_constraint_raw, direct_momentum_raw
    ))
)
check("the direct linear path is the same interval-independent exponential tangent", tangent_ok)


u = v**2
radical_data = (
    (3 * u**2 + 20 * u + 32, sp.sqrt(u + 4) * sp.sqrt(3 * u + 8)),
    (
        3 * u**4 + 38 * u**3 + 179 * u**2 + 372 * u + 288,
        (u + 3) * sp.sqrt(u + 4) * sp.sqrt(3 * u + 8),
    ),
    (
        9 * u**4 + 120 * u**3 + 592 * u**2 + 1280 * u + 1024,
        (u + 4) * (3 * u + 8),
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
        27 * u**4 + 360 * u**3 + 1776 * u**2 + 3840 * u + 3072,
        sp.sqrt(3) * (u + 4) * (3 * u + 8),
    ),
    (
        3 * u**5 + 50 * u**4 + 331 * u**3 + 1088 * u**2
        + 1776 * u + 1152,
        (u + 3) * (u + 4) * sp.sqrt(3 * u + 8),
    ),
)
factorization_ok = all(
    sp.factor(polynomial - replacement**2) == 0
    for polynomial, replacement in radical_data
)


def normalize_positive_radicals(expression):
    def is_sqrt(node):
        return isinstance(node, sp.Pow) and node.exp == sp.Rational(1, 2)

    def replace_sqrt(node):
        base = sp.factor(node.base)
        for polynomial, replacement in radical_data:
            if sp.expand(base - polynomial) == 0:
                return replacement
        return node

    return expression.replace(is_sqrt, replace_sqrt)


direct_constraint = sp.factor(normalize_positive_radicals(direct_constraint_raw))
direct_momentum = sp.factor(normalize_positive_radicals(direct_momentum_raw))
check("all preregistered positive radical normalizations are exact", factorization_ok)


# Parse primary formulas only after the independent direct limits exist.
locals_map = {
    "v": v,
    "mu": mu,
    "pi": sp.pi,
    "sqrt": sp.sqrt,
    "acos": sp.acos,
    "asinh": sp.asinh,
}
primary_action = sp.sympify(primary["leading"]["action"], locals=locals_map)
primary_constraint = sp.sympify(primary["leading"]["constraint"], locals=locals_map)
primary_momentum = sp.sympify(primary["leading"]["momentum"], locals=locals_map)
primary_mass = sp.sympify(primary["leading"]["mass_branch"], locals=locals_map)

action_match = sp.simplify(sp.powsimp(direct_action - primary_action, force=True)) == 0
constraint_match = sp.simplify(
    sp.powsimp(direct_constraint - primary_constraint, force=True)
) == 0
momentum_match = sp.simplify(
    sp.powsimp(direct_momentum - primary_momentum, force=True)
) == 0
direct_matches = bool(action_match and constraint_match and momentum_match)
check(
    "the exact derivative-first limits match all three frozen primary formulas",
    direct_matches,
    f"action={action_match}; constraint={constraint_match}; momentum={momentum_match}",
)


mu_coefficient = sp.diff(direct_constraint, mu)
direct_mass = sp.factor(
    -direct_constraint.subs(mu, 0) / mu_coefficient
)
mass_match = sp.simplify(sp.powsimp(direct_mass - primary_mass, force=True)) == 0
cosine = (v**2 + 2) / (2 * (v**2 + 3))
positive_domain = bool(
    sp.simplify(cosine - sp.Rational(1, 3) - v**2 / (6 * (v**2 + 3))) == 0
    and sp.simplify(sp.Rational(1, 2) - cosine - 1 / (2 * (v**2 + 3))) == 0
)
branch_ok = bool(
    mu_coefficient == -8 * sp.pi
    and sp.simplify(direct_constraint.subs(mu, direct_mass)) == 0
    and mass_match
    and positive_domain
)
check("the direct constraint independently has the same unique positive mass branch", branch_ok)


epsilon_static = 2 * sp.pi - 5 * sp.acos(sp.Rational(1, 3))
symmetry_ok = bool(
    sp.simplify(direct_mass.subs(v, -v) - direct_mass) == 0
    and sp.simplify(direct_momentum.subs(v, -v) + direct_momentum) == 0
    and sp.simplify(direct_mass.subs(v, 0) - 90 * epsilon_static / sp.pi) == 0
    and sp.simplify(direct_momentum.subs(v, 0)) == 0
)
check("the derivative-first branch passes static and time-reversal controls", symmetry_ok)


qroot = sp.sqrt(3 * v**2 + 8)
omitted_boundary_defect = -3600 * v**2 / ((v**2 + 3) * qroot)
mass_shift_defect = sp.simplify(
    direct_constraint.subs(mu, direct_mass + sp.Rational(1, 10))
)
hostile_ok = bool(
    all(
        abs(sp.N(omitted_boundary_defect.subs(v, sp.Rational(value)), 80)) > 0
        for value in (-sp.Rational(7, 10), sp.Rational(3, 10), sp.Rational(13, 10))
    )
    and sp.simplify(mass_shift_defect + 4 * sp.pi / 5) == 0
)
check("both deleted-term and shifted-mass hostile controls fail as required", hostile_ok)


# New 100-decimal controls on the literal exponential path.
action_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), ACTION, "mpmath")
f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), F_EXACT, "mpmath")
p_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_PRE_EXACT, "mpmath")


def exact_numeric(v_value):
    radius = mp.sqrt(v_value**2 + 4)
    theta = mp.acos((v_value**2 + 2) / (2 * (v_value**2 + 3)))
    eta = mp.asinh(v_value / mp.sqrt(8 * (v_value**2 + 3)))
    epsilon = 2 * mp.pi - 5 * theta
    mass = 180 * epsilon / (mp.pi * radius)
    action = 360 * radius * epsilon - 1200 * mp.sqrt(3) * v_value * eta - 8 * mp.pi * mass
    momentum = 180 * v_value * epsilon / radius - 600 * mp.sqrt(3) * eta
    return mass, action, mp.mpf(0), momentum


numeric_records = {}
numeric_ok = True
for v_text in VELOCITY_TEXTS:
    v_value = mp.mpf(v_text)
    mass_value, action_limit, constraint_limit, momentum_limit = exact_numeric(v_value)
    numeric_records[v_text] = {}
    for s_numerator, s_denominator in S_RATIONALS:
        s_text = f"{s_numerator}/{s_denominator}"
        s_value = mp.mpf(s_numerator) / s_denominator
        errors = {name: [] for name in ("action", "constraint", "momentum")}
        for e_numerator, e_denominator in E_RATIONALS:
            e_value = mp.mpf(e_numerator) / e_denominator
            tau_value = s_value * e_value
            point = (1, mp.exp(v_value * tau_value), tau_value**2, mass_value)
            direct = {
                "action": action_numeric(*point) / tau_value,
                "constraint": 2 * f_numeric(*point) / tau_value,
                "momentum": p_numeric(*point),
            }
            limits = {
                "action": action_limit,
                "constraint": constraint_limit,
                "momentum": momentum_limit,
            }
            for name in errors:
                errors[name].append(
                    abs(direct[name] - limits[name]) / max(mp.mpf(1), abs(limits[name]))
                )
        orders = {}
        for name, pair in errors.items():
            if pair[0] < mp.mpf("1e-70") and pair[1] < mp.mpf("1e-70"):
                orders[name] = mp.inf
                order_ok = True
            elif pair[0] >= mp.mpf("1e-70") and pair[1] >= mp.mpf("1e-70"):
                orders[name] = mp.log(pair[0] / pair[1]) / mp.log(2)
                order_ok = bool(
                    pair[1] < pair[0]
                    and mp.mpf("0.8") <= orders[name] <= mp.mpf("1.2")
                )
            else:
                orders[name] = mp.nan
                order_ok = False
            numeric_ok &= order_ok
        numeric_records[v_text][s_text] = {
            "errors": {name: [text(value, 40) for value in pair] for name, pair in errors.items()},
            "orders": {name: text(value, 30) for name, value in orders.items()},
        }
check("all new direct controls converge at the frozen first order", numeric_ok)


all_controls = bool(
    provenance_ok
    and tangent_ok
    and factorization_ok
    and direct_matches
    and branch_ok
    and symmetry_ok
    and hostile_ok
    and numeric_ok
)
outcome = (
    "GENERIC_VELOCITY_LEADING_REPARAMETRIZATION_ADVERSARIALLY_CORROBORATED"
    if all_controls
    else "GENERIC_VELOCITY_ADVERSARIAL_DISAGREEMENT"
)
check(
    "the adversarial hierarchy assigns the corroborated generic-velocity verdict",
    outcome == "GENERIC_VELOCITY_LEADING_REPARAMETRIZATION_ADVERSARIALLY_CORROBORATED",
    outcome,
)

artifact = {
    "primary_input_sha256": digest(PRIMARY_INPUT),
    "primary_implementation_commit": PRIMARY_IMPLEMENTATION_COMMIT,
    "adversarial_protocol_commit": ADVERSARIAL_PROTOCOL_COMMIT,
    "radical_protocol_commit": RADICAL_PROTOCOL_COMMIT,
    "method": "differentiate_full_action_before_direct_tangent_limit",
    "direct": {
        "action": str(direct_action),
        "constraint": str(direct_constraint),
        "momentum": str(direct_momentum),
        "mass_branch": str(direct_mass),
        "matches_primary": {
            "action": action_match,
            "constraint": constraint_match,
            "momentum": momentum_match,
            "mass": mass_match,
        },
    },
    "hostile": {
        "omitted_boundary_defect": str(omitted_boundary_defect),
        "mass_shift_defect": str(mass_shift_defect),
    },
    "numeric_controls": numeric_records,
    "labels": {
        "generic_velocity_leading_reparametrization": "DERIVED_EXACT_ADVERSARIALLY_CORROBORATED_STRUCTURAL",
        "next_order_composition": "OPEN",
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

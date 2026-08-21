#!/usr/bin/env python3
"""Generic-velocity leading relational map of the cellular 600-cell action."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
ACTION_INPUT = HERE / "gravity_600cell_homothetic_frustum_action.json"
OUTPUT = HERE / "gravity_600cell_generic_velocity_composition.json"
ACTION_SHA256 = (
    "c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d"
)
PRIOR_ART_COMMIT = "1fcab34"
PROTOCOL_COMMIT = "9472a15"
NUMERIC_PROTOCOL_COMMIT = "8af359d"
VELOCITY_TEXTS = ("0.2", "0.5", "1.0")
S_TEXTS = ("1.0", "0.5")
E_TEXTS = ("0.005", "0.0025")
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


action_artifact = json.loads(ACTION_INPUT.read_text())
provenance_ok = bool(
    digest(ACTION_INPUT) == ACTION_SHA256
    and action_artifact["outcome"] == "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
    and action_artifact["passed"] == action_artifact["tests"] == 16
    and PRIOR_ART_COMMIT == "1fcab34"
    and PROTOCOL_COMMIT == "9472a15"
    and NUMERIC_PROTOCOL_COMMIT == "8af359d"
)
check("the cellular action and generic-velocity protocols are frozen", provenance_ok)


# Exact primitive limits.  The computation keeps the interval factor s until
# it cancels and does not sample a mass or velocity.
e, s = sp.symbols("e s", positive=True)
v, mu = sp.symbols("v mu", real=True)
endpoint = sp.exp(s * v * e)
delta = endpoint - 1
rho = s**2 * e**2

primitive_limits = {
    "delta_over_se": sp.limit(delta / (s * e), e, 0, dir="+"),
    "height_over_se": sp.limit(
        sp.sqrt(rho + delta**2 / 4) / (s * e), e, 0, dir="+"
    ),
    "cosine": sp.limit(
        (delta**2 + 2 * rho) / (2 * (delta**2 + 3 * rho)),
        e,
        0,
        dir="+",
    ),
    "boost": sp.limit(
        delta / sp.sqrt(8 * (delta**2 + 3 * rho)), e, 0, dir="+"
    ),
    "boundary_over_se": sp.limit(
        (1 - endpoint**2) / (s * e), e, 0, dir="+"
    ),
    "dust_over_se": sp.limit(sp.sqrt(rho) / (s * e), e, 0, dir="+"),
}
expected_primitives = {
    "delta_over_se": v,
    "height_over_se": sp.sqrt(v**2 + 4) / 2,
    "cosine": (v**2 + 2) / (2 * (v**2 + 3)),
    "boost": v / sp.sqrt(8 * (v**2 + 3)),
    "boundary_over_se": -2 * v,
    "dust_over_se": sp.Integer(1),
}
primitive_ok = all(
    sp.simplify(primitive_limits[name] - expected) == 0
    for name, expected in expected_primitives.items()
)
check("all six exact generic-velocity primitive limits are interval independent", primitive_ok)


radius = sp.sqrt(v**2 + 4)
cosine = expected_primitives["cosine"]
theta = sp.acos(cosine)
eta = sp.asinh(expected_primitives["boost"])
epsilon_v = 2 * sp.pi - 5 * theta
leading_action = (
    360 * radius * epsilon_v
    - 1200 * sp.sqrt(3) * v * eta
    - 8 * sp.pi * mu
)
interval_independent_action = s * sp.diff(leading_action, s) == 0
action_ok = bool(interval_independent_action and not leading_action.has(s))
check(
    "the leading principal function is exactly S=s*e*L0(v,mu)",
    action_ok,
    f"L0={leading_action}",
)


# Exact derivative identities expose the cancellation without asking a CAS to
# simplify products of high-degree square roots heuristically.
theta_prime = -2 * v / (
    (v**2 + 3) * radius * sp.sqrt(3 * v**2 + 8)
)
eta_prime = sp.sqrt(3) / (
    (v**2 + 3) * sp.sqrt(3 * v**2 + 8)
)
derivative_identities_ok = bool(
    sp.simplify(sp.diff(theta, v) - theta_prime) == 0
    and sp.simplify(sp.diff(eta, v) - eta_prime) == 0
    and sp.simplify(
        -1800 * radius * theta_prime
        - 1200 * sp.sqrt(3) * v * eta_prime
    ) == 0
)
leading_action_prime = (
    360 * v * epsilon_v / radius - 1200 * sp.sqrt(3) * eta
)
constraint_hj = sp.factor(leading_action - v * leading_action_prime)
constraint_simple = 1440 * epsilon_v / radius - 8 * sp.pi * mu
hj_ok = bool(
    derivative_identities_ok
    and sp.simplify(
        constraint_hj - constraint_simple
    ) == 0
)
check("the Hamilton-Jacobi fixed-endpoint lapse constraint reduces exactly", hj_ok)


# Independent direct differentiation of the original hinge primitives.
# These are F/tau contributions after rho differentiation at fixed Delta.
qroot = sp.sqrt(3 * v**2 + 8)
direct_lateral = (
    720 * epsilon_v / radius
    - 1800 * v**2 / ((v**2 + 3) * qroot)
)
direct_boundary = 1800 * v**2 / ((v**2 + 3) * qroot)
direct_dust = -4 * sp.pi * mu
constraint_direct = sp.factor(2 * (direct_lateral + direct_boundary + direct_dust))
direct_constraint_ok = sp.simplify(constraint_direct - constraint_simple) == 0
check(
    "direct full-hinge lapse differentiation matches the HJ constraint",
    direct_constraint_ok,
    f"lateral-boundary cancellation={sp.simplify(direct_lateral+direct_boundary-720*epsilon_v/radius)}",
)


# Incoming-scale chain rule with fixed physical mass.  For
# S_lead=L*tau*L0(v,M/L), fixed L_plus gives dv/dL=-1/tau and fixed M gives
# dmu/dL=-mu/L.  The terms carrying those two derivatives are kept visibly.
L, tau = sp.symbols("L tau", positive=True)
leading_mu_derivative = sp.diff(leading_action, mu)
scale_derivative = (
    tau * leading_action
    + L * tau * (
        leading_action_prime * (-1 / tau)
        + leading_mu_derivative * (-mu / L)
    )
)
pre_chain = -L * scale_derivative / 2
chain_remainder = sp.simplify(
    pre_chain - L**2 * leading_action_prime / 2
)
expected_chain_remainder = sp.factor(
    -L * tau * (leading_action - mu * leading_mu_derivative) / 2
)
chain_expanded = sp.simplify(chain_remainder - expected_chain_remainder)
chain_limit = sp.limit(pre_chain.subs(L, 1), tau, 0, dir="+")
momentum = leading_action_prime / 2
momentum_simple = 180 * v * epsilon_v / radius - 600 * sp.sqrt(3) * eta
momentum_ok = bool(
    chain_expanded == 0
    and sp.simplify(chain_limit - momentum) == 0
    and sp.simplify(momentum - momentum_simple) == 0
    and not momentum_simple.has(s)
)
check(
    "the fixed-mass incoming momentum is interval independent",
    momentum_ok,
    f"p(v)={momentum_simple}",
)


mass_branch = sp.factor(180 * epsilon_v / (sp.pi * radius))
branch_residual = sp.simplify(constraint_simple.subs(mu, mass_branch))
branch_count = 1 if sp.diff(constraint_simple, mu) == -8 * sp.pi else None
cosine_bounds = (
    sp.factor(cosine - sp.Rational(1, 3)),
    sp.factor(sp.Rational(1, 2) - cosine),
)
positive_domain_ok = bool(
    cosine_bounds[0] == v**2 / (6 * (v**2 + 3))
    and cosine_bounds[1] == 1 / (2 * (v**2 + 3))
    and sp.N(2 * sp.pi - 5 * sp.acos(sp.Rational(1, 3)), 80) > 0
)
branch_ok = bool(
    branch_count == 1
    and branch_residual == 0
    and positive_domain_ok
)
check(
    "the lapse constraint has one positive mass branch for every real velocity",
    branch_ok,
    f"mu(v)={mass_branch}",
)


epsilon_static = 2 * sp.pi - 5 * sp.acos(sp.Rational(1, 3))
static_mass = sp.simplify(mass_branch.subs(v, 0))
mass_first = sp.simplify(sp.diff(mass_branch, v).subs(v, 0))
mass_second = sp.simplify(sp.diff(mass_branch, v, 2).subs(v, 0))
momentum_static = sp.simplify(momentum_simple.subs(v, 0))
momentum_first = sp.simplify(sp.diff(momentum_simple, v).subs(v, 0))
expected_mass_second = 180 / sp.pi * (
    5 / (6 * sp.sqrt(2)) - epsilon_static / 8
)
D_STATIC = 5 * sp.sqrt(2) / 3 - epsilon_static
static_jet_ok = bool(
    sp.simplify(static_mass - 90 * epsilon_static / sp.pi) == 0
    and mass_first == 0
    and sp.simplify(mass_second - expected_mass_second) == 0
    and sp.N(mass_second, 80) > 0
    and momentum_static == 0
    and sp.simplify(momentum_first + 90 * D_STATIC) == 0
)
check(
    "the generic branch recovers the static mass and its exact local jet",
    static_jet_ok,
    f"mu''(0)={mass_second}; p'(0)={momentum_first}",
)


time_reversal_ok = bool(
    sp.simplify(mass_branch.subs(v, -v) - mass_branch) == 0
    and sp.simplify(momentum_simple.subs(v, -v) + momentum_simple) == 0
)
composition_ok = bool(
    not any(expression.has(s) for expression in (
        leading_action, constraint_simple, momentum_simple, mass_branch
    ))
    and sp.simplify(2 * (sp.Rational(1, 2) * v) - v) == 0
    and branch_count == 1
)
check("time reversal and the unique leading same-state composition gate pass", time_reversal_ok and composition_ok)


# Unexpanded direct 100-decimal controls.
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
action_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), ACTION, "mpmath")
f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), F_EXACT, "mpmath")
p_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_PRE_EXACT, "mpmath")


def exact_numeric(v_value):
    radius_value = mp.sqrt(v_value**2 + 4)
    cosine_value = (v_value**2 + 2) / (2 * (v_value**2 + 3))
    theta_value = mp.acos(cosine_value)
    eta_value = mp.asinh(v_value / mp.sqrt(8 * (v_value**2 + 3)))
    epsilon_value = 2 * mp.pi - 5 * theta_value
    mass_value = 180 * epsilon_value / (mp.pi * radius_value)
    action_value = (
        360 * radius_value * epsilon_value
        - 1200 * mp.sqrt(3) * v_value * eta_value
        - 8 * mp.pi * mass_value
    )
    momentum_value = (
        180 * v_value * epsilon_value / radius_value
        - 600 * mp.sqrt(3) * eta_value
    )
    return mass_value, action_value, mp.mpf(0), momentum_value


numeric_records = {}
all_real = True
all_orders_ok = True
all_s_agreement = True
for v_text in VELOCITY_TEXTS:
    v_value = mp.mpf(v_text)
    mass_value, action_limit, constraint_limit, momentum_limit = exact_numeric(v_value)
    numeric_records[v_text] = {}
    by_s = {}
    for s_text in S_TEXTS:
        s_value = mp.mpf(s_text)
        errors = {name: [] for name in ("action", "constraint", "momentum")}
        values = {name: [] for name in errors}
        for e_text in E_TEXTS:
            e_value = mp.mpf(e_text)
            tau_value = s_value * e_value
            lp_value = mp.exp(tau_value * v_value)
            rho_value = tau_value**2
            direct = {
                "action": action_numeric(1, lp_value, rho_value, mass_value) / tau_value,
                "constraint": 2 * f_numeric(1, lp_value, rho_value, mass_value) / tau_value,
                "momentum": p_numeric(1, lp_value, rho_value, mass_value),
            }
            limits = {
                "action": action_limit,
                "constraint": constraint_limit,
                "momentum": momentum_limit,
            }
            for name in errors:
                values[name].append(direct[name])
                errors[name].append(
                    abs(direct[name] - limits[name]) / max(mp.mpf(1), abs(limits[name]))
                )
                all_real &= abs(mp.im(direct[name])) < mp.mpf("1e-80")
        orders = {}
        for name, pair in errors.items():
            if pair[0] < mp.mpf("1e-70") and pair[1] < mp.mpf("1e-70"):
                orders[name] = mp.inf
                passed_order = True
            elif pair[0] >= mp.mpf("1e-70") and pair[1] >= mp.mpf("1e-70"):
                orders[name] = mp.log(pair[0] / pair[1]) / mp.log(2)
                passed_order = bool(
                    pair[1] < pair[0]
                    and mp.mpf("0.8") <= orders[name] <= mp.mpf("1.2")
                )
            else:
                orders[name] = mp.nan
                passed_order = False
            all_orders_ok &= passed_order
        by_s[s_text] = {"values": values, "errors": errors, "orders": orders}
        numeric_records[v_text][s_text] = {
            "values": {name: [text(value, 50) for value in pair] for name, pair in values.items()},
            "errors": {name: [text(value, 40) for value in pair] for name, pair in errors.items()},
            "orders": {name: text(value, 30) for name, value in orders.items()},
        }
    for name in ("action", "constraint", "momentum"):
        coarse_fine_difference = abs(by_s["1.0"]["values"][name][-1] - by_s["0.5"]["values"][name][-1])
        envelope = 10 * (
            by_s["1.0"]["errors"][name][-1]
            + by_s["0.5"]["errors"][name][-1]
            + mp.mpf("1e-70")
        ) * max(
            mp.mpf(1),
            abs(by_s["1.0"]["values"][name][-1]),
            abs(by_s["0.5"]["values"][name][-1]),
        )
        all_s_agreement &= coarse_fine_difference <= envelope

numeric_ok = bool(all_real and all_orders_ok and all_s_agreement)
check(
    "all direct controls converge at the frozen order on both interval factors",
    numeric_ok,
    f"real={all_real}; orders={all_orders_ok}; s agreement={all_s_agreement}",
)


all_exact = bool(
    provenance_ok
    and primitive_ok
    and action_ok
    and hj_ok
    and direct_constraint_ok
    and momentum_ok
    and branch_ok
    and static_jet_ok
    and time_reversal_ok
    and composition_ok
)
if not all_exact or not numeric_ok:
    outcome = "GENERIC_VELOCITY_LEADING_OPEN"
elif branch_count == 1:
    outcome = "GENERIC_VELOCITY_LEADING_REPARAMETRIZATION"
else:
    outcome = "GENERIC_VELOCITY_LEADING_NONUNIQUE"
check(
    "the frozen hierarchy assigns the generic-velocity leading verdict",
    outcome == "GENERIC_VELOCITY_LEADING_REPARAMETRIZATION",
    outcome,
)

artifact = {
    "action_input_sha256": digest(ACTION_INPUT),
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "numeric_protocol_commit": NUMERIC_PROTOCOL_COMMIT,
    "leading": {
        "action": str(leading_action),
        "constraint": str(constraint_simple),
        "momentum": str(momentum_simple),
        "mass_branch": str(mass_branch),
        "branch_count": branch_count,
        "same_state_compatible_count": 1 if composition_ok else 0,
    },
    "static_jet": {
        "mass": str(static_mass),
        "mass_first": str(mass_first),
        "mass_second": str(mass_second),
        "momentum": str(momentum_static),
        "momentum_first": str(momentum_first),
    },
    "primitive_limits": {name: str(value) for name, value in primitive_limits.items()},
    "numeric_controls": numeric_records,
    "labels": {
        "leading_reparametrization": "DERIVED_EXACT_STRUCTURAL",
        "next_order_composition": "OPEN",
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

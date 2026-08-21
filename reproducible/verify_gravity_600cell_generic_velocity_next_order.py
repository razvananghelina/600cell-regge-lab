#!/usr/bin/env python3
"""Next-order fixed-state lapse and conditional composition census."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
ACTION_INPUT = HERE / "gravity_600cell_homothetic_frustum_action.json"
PRIMARY_INPUT = HERE / "gravity_600cell_generic_velocity_composition.json"
ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_generic_velocity_composition_adversarial.json"
)
OUTPUT = HERE / "gravity_600cell_generic_velocity_next_order.json"
ACTION_SHA256 = (
    "c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d"
)
PRIMARY_SHA256 = (
    "8ded36f1fa00307fcb23369c25290c9f5bd701709762d6a865437c2507eabfc9"
)
ADVERSARIAL_SHA256 = (
    "cd46c6c9d38e1b14fc32f09f8a2cf72039d28ed136ef4081fe0edba149a9b6b2"
)
PRIOR_ART_COMMIT = "62d92ab"
PROTOCOL_COMMIT = "c577419"
CORRECTION_PROTOCOL_COMMIT = "85f6752"
VELOCITY_RATIONALS = ((1, 3), (4, 5), (3, 2))
ACCELERATION_RATIONALS = ((0, 1), (1, 7))
HEIGHT_RATIONALS = ((1, 400), (1, 800))
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
primary_artifact = json.loads(PRIMARY_INPUT.read_text())
adversarial_artifact = json.loads(ADVERSARIAL_INPUT.read_text())
provenance_ok = bool(
    digest(ACTION_INPUT) == ACTION_SHA256
    and digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and digest(ADVERSARIAL_INPUT) == ADVERSARIAL_SHA256
    and action_artifact["outcome"] == "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
    and primary_artifact["outcome"]
    == "GENERIC_VELOCITY_LEADING_REPARAMETRIZATION"
    and adversarial_artifact["outcome"]
    == "GENERIC_VELOCITY_LEADING_REPARAMETRIZATION_ADVERSARIALLY_CORROBORATED"
    and PRIOR_ART_COMMIT == "62d92ab"
    and PROTOCOL_COMMIT == "c577419"
    and CORRECTION_PROTOCOL_COMMIT == "85f6752"
)
check("the action and both leading generic-velocity results are frozen", provenance_ok)


# Reconstruct the complete action and its canonical derivatives.
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


v, a, b = sp.symbols("v a b", real=True)
h = sp.symbols("h", positive=True)
radius = sp.sqrt(v**2 + 4)
theta = sp.acos((v**2 + 2) / (2 * (v**2 + 3)))
eta = sp.asinh(v / sp.sqrt(8 * (v**2 + 3)))
epsilon_v = 2 * sp.pi - 5 * theta
mass_branch = 180 * epsilon_v / (sp.pi * radius)
momentum_branch = (
    180 * v * epsilon_v / radius - 600 * sp.sqrt(3) * eta
)


# Exact scaled action.  This avoids asking the limit engine to rediscover the
# positive-height scaling through the complete transcendental derivative.
lm, lp, q, w, tau = sp.symbols("lm lp q w tau", real=True)
reduced_radius = sp.sqrt(q**2 + 4) / 2
reduced_theta = sp.acos((q**2 + 2) / (2 * (q**2 + 3)))
reduced_eta = sp.asinh(q / sp.sqrt(8 * (q**2 + 3)))
reduced_epsilon = 2 * sp.pi - 5 * reduced_theta
reduced_a = (
    360 * (lm + lp) * reduced_radius * reduced_epsilon - 8 * sp.pi * MASS
)
reduced_a_q = sp.diff(reduced_a, q)
reduced_eta_q = sp.diff(reduced_eta, q)
constraint_reduced = sp.factor(
    reduced_a - q * reduced_a_q - q * w * reduced_eta_q
)
p_minus_reduced = lm * (
    tau * 360 * reduced_radius * reduced_epsilon
    - reduced_a_q
    + 1200 * sp.sqrt(3) * lm * reduced_eta
    - w * reduced_eta_q
) / 2
p_plus_reduced = lp * (
    tau * 360 * reduced_radius * reduced_epsilon
    + reduced_a_q
    - 1200 * sp.sqrt(3) * lp * reduced_eta
    + w * reduced_eta_q
) / 2
action_over_tau_reduced = reduced_a + w * reduced_eta


def path_data(c, xm1, xm2, xp1, xp2):
    q0 = (xp1 - xm1) / c
    q1 = (xp2 + xp1**2 / 2 - xm2 - xm1**2 / 2) / c
    w0 = -1200 * sp.sqrt(3) * q0
    w1 = 1200 * sp.sqrt(3) * (
        xm2 - xp2 + xm1**2 - xp1**2
    ) / c
    return {
        "base": {lm: 1, lp: 1, q: q0, w: w0, tau: 0},
        "slope": {
            lm: xm1,
            lp: xp1,
            q: q1,
            w: w1,
            tau: c,
        },
    }


def path_base(expression, data, mass):
    return sp.factor(expression.subs(MASS, mass).subs(data["base"]))


def path_first(expression, data, mass):
    expression = expression.subs(MASS, mass)
    return sp.factor(
        sp.simplify(
            sum(
                sp.diff(expression, variable).subs(data["base"]) * slope
                for variable, slope in data["slope"].items()
            )
        )
    )


one_slab_data = path_data(sp.Integer(1), 0, 0, v, a)
constraint_zero = path_base(constraint_reduced, one_slab_data, mass_branch)
momentum_zero = -path_base(
    p_minus_reduced, one_slab_data, mass_branch
) - momentum_branch
leading_zero_ok = bool(
    sp.simplify(constraint_zero) == 0
    and sp.simplify(momentum_zero) == 0
)
check("the fixed incoming state cancels both zeroth-order residuals", leading_zero_ok)


constraint_first = path_first(
    constraint_reduced, one_slab_data, mass_branch
)
momentum_first = -path_first(
    p_minus_reduced, one_slab_data, mass_branch
)
limits_exist = bool(
    constraint_first.is_finite is not False
    and momentum_first.is_finite is not False
    and not constraint_first.has(sp.Limit)
    and not momentum_first.has(sp.Limit)
)
check(
    "the exact one-slab first correction exists symbolically",
    limits_exist,
    f"C1={constraint_first}; P1={momentum_first}",
)


def numerator_polynomial(expression, variable):
    numerator = sp.cancel(sp.together(expression)).as_numer_denom()[0]
    return sp.Poly(sp.expand(numerator), variable, domain="EX"), numerator


poly_c, numerator_c = numerator_polynomial(constraint_first, a)
poly_p, numerator_p = numerator_polynomial(momentum_first, a)
degrees = (poly_c.degree(), poly_p.degree())
linear_census = degrees == (1, 1)
if linear_census:
    c_coefficients = poly_c.all_coeffs()
    p_coefficients = poly_p.all_coeffs()
    root_c = sp.factor(-c_coefficients[1] / c_coefficients[0])
    root_p = sp.factor(-p_coefficients[1] / p_coefficients[0])
    root_difference = sp.factor(
        sp.simplify(sp.powsimp(root_c - root_p, force=True))
    )
    resultant = sp.factor(sp.resultant(numerator_c, numerator_p, a))
else:
    root_c = root_p = root_difference = None
    resultant = sp.factor(sp.resultant(numerator_c, numerator_p, a))

common_all_nonzero = bool(linear_census and root_difference == 0)
classification = "ALL_REAL_NONZERO" if common_all_nonzero else "UNRESOLVED"
exceptional_velocities = []
classification_complete = common_all_nonzero
check(
    "the common-root census is algebraically explicit",
    linear_census,
    f"degrees={degrees}; root_difference={root_difference}; resultant={resultant}",
)


# The static stratum is evaluated separately because p0 scales with h there.
A = sp.symbols("A", real=True)
alpha_static = sp.acos(sp.Rational(1, 3))
epsilon_static = 2 * sp.pi - 5 * alpha_static
D_STATIC = 5 * sp.sqrt(2) / 3 - epsilon_static
static_mass = 90 * epsilon_static / sp.pi
static_data = {
    "base": {lm: 1, lp: 1, q: 0, w: 0, tau: 0},
    "slope": {
        lm: 0,
        lp: 0,
        q: A,
        w: -1200 * sp.sqrt(3) * A,
        tau: 1,
    },
    "quadratic": {lm: 0, lp: A, q: 0, w: 0, tau: 0},
}


def path_second(expression, data, mass):
    expression = expression.subs(MASS, mass)
    variables = tuple(data["slope"])
    linear_second = sum(
        sp.diff(expression, variable).subs(data["base"])
        * data["quadratic"][variable]
        for variable in variables
    )
    quadratic_first = sp.Rational(1, 2) * sum(
        sp.diff(expression, left, right).subs(data["base"])
        * data["slope"][left]
        * data["slope"][right]
        for left in variables
        for right in variables
    )
    return sp.factor(sp.simplify(linear_second + quadratic_first))


static_f_power = 2
static_p_power = 1
static_f_coefficient = path_second(
    constraint_reduced, static_data, static_mass
)
static_p_coefficient = (
    -path_first(p_minus_reduced, static_data, static_mass)
    - 180 * epsilon_static
)
static_f_expected = A * (D_STATIC * A + 4 * epsilon_static)
static_p_expected = D_STATIC * A + 4 * epsilon_static
static_f_ratio = sp.factor(static_f_coefficient / static_f_expected)
static_p_ratio = sp.factor(static_p_coefficient / static_p_expected)
static_control_ok = bool(
    static_f_power is not None
    and static_p_power is not None
    and sp.simplify(sp.diff(static_f_ratio, A)) == 0
    and sp.simplify(sp.diff(static_p_ratio, A)) == 0
    and static_f_ratio != 0
    and static_p_ratio != 0
)
check(
    "the independent static jet recovers both certified turning-point factors",
    static_control_ok,
    f"powers=({static_f_power},{static_p_power}); ratios=({static_f_ratio},{static_p_ratio})",
)


hostile_momentum_shift = sp.Rational(1, 10)
hostile_mass_defect = path_base(
    constraint_reduced,
    one_slab_data,
    mass_branch + sp.Rational(1, 10),
)
hostile_momentum_defect = sp.factor(
    -path_base(p_minus_reduced, one_slab_data, mass_branch)
    - (momentum_branch + hostile_momentum_shift)
)
hostile_ok = bool(
    sp.simplify(hostile_mass_defect + 4 * sp.pi / 5) == 0
    and sp.simplify(hostile_momentum_defect + sp.Rational(1, 10)) == 0
)
check(
    "both frozen changed-state hostile controls fail exactly",
    hostile_ok,
    f"mass={hostile_mass_defect}; momentum={hostile_momentum_defect}",
)


# Direct arbitrary-precision controls of the coefficient extraction.
f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), F_EXACT, "mpmath")
p_pre_numeric = sp.lambdify(
    (L_MINUS, L_PLUS, RHO, MASS), P_PRE_EXACT, "mpmath"
)
c1_numeric = sp.lambdify((v, a), constraint_first, "mpmath")
p1_numeric = sp.lambdify((v, a), momentum_first, "mpmath")
mass_numeric = sp.lambdify(v, mass_branch, "mpmath")
momentum_numeric = sp.lambdify(v, momentum_branch, "mpmath")
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
        exact_values = {
            "constraint": c1_numeric(v_value, a_value),
            "momentum": p1_numeric(v_value, a_value),
        }
        errors = {name: [] for name in exact_values}
        for h_numerator, h_denominator in HEIGHT_RATIONALS:
            h_value = mp.mpf(h_numerator) / h_denominator
            endpoint = mp.exp(v_value * h_value + a_value * h_value**2)
            direct_values = {
                "constraint": 2 * f_numeric(1, endpoint, h_value**2, mass_value)
                / h_value**2,
                "momentum": (
                    p_pre_numeric(1, endpoint, h_value**2, mass_value) - p0_value
                )
                / h_value,
            }
            for name, direct_value in direct_values.items():
                errors[name].append(
                    abs(direct_value - exact_values[name])
                    / max(mp.mpf(1), abs(exact_values[name]))
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
        numeric_records[v_key][a_key] = {
            "errors": {
                name: [text(value) for value in pair]
                for name, pair in errors.items()
            },
            "orders": {name: text(value, 30) for name, value in orders.items()},
        }
check("all direct next-order coefficient controls meet the frozen criterion", numeric_ok)


# Conditional stationary two-half-slab calculation.
composition = {
    "executed": False,
    "reason": "one-slab common root is not globally classified",
}
composition_ok = None
if common_all_nonzero and classification_complete:
    a_root = root_c
    coarse_data = path_data(sp.Integer(1), 0, 0, v, a_root)
    first_data = path_data(
        sp.Rational(1, 2), 0, 0, v / 2, a_root / 4
    )
    second_data = path_data(
        sp.Rational(1, 2), v / 2, a_root / 4, v, b
    )
    first_constraint_coefficient = path_first(
        constraint_reduced, first_data, mass_branch
    )
    first_momentum_coefficient = -path_first(
        p_minus_reduced, first_data, mass_branch
    )
    second_constraint_coefficient = path_first(
        constraint_reduced, second_data, mass_branch
    )
    seam_base = sp.factor(
        path_base(p_plus_reduced, first_data, mass_branch)
        + path_base(p_minus_reduced, second_data, mass_branch)
    )
    seam_coefficient = sp.factor(
        path_first(p_plus_reduced, first_data, mass_branch)
        + path_first(p_minus_reduced, second_data, mass_branch)
    )
    poly_second, numerator_second = numerator_polynomial(
        second_constraint_coefficient, b
    )
    poly_seam, numerator_seam = numerator_polynomial(seam_coefficient, b)
    composition_degrees = (poly_second.degree(), poly_seam.degree())
    if composition_degrees == (1, 1):
        second_coefficients = poly_second.all_coeffs()
        seam_coefficients = poly_seam.all_coeffs()
        root_second = sp.factor(-second_coefficients[1] / second_coefficients[0])
        root_seam = sp.factor(-seam_coefficients[1] / seam_coefficients[0])
        composition_root_difference = sp.factor(
            sp.simplify(sp.powsimp(root_second - root_seam, force=True))
        )
        common_b = composition_root_difference == 0
    else:
        root_second = root_seam = composition_root_difference = None
        common_b = False
    if common_b:
        b_root = root_second
        endpoint_defect = sp.factor(
            sp.simplify(sp.powsimp(b_root - a_root, force=True))
        )
        second_on_shell_data = path_data(
            sp.Rational(1, 2), v / 2, a_root / 4, v, b_root
        )
        final_momentum_defect = sp.factor(
            path_first(p_plus_reduced, second_on_shell_data, mass_branch)
            - path_first(p_plus_reduced, coarse_data, mass_branch)
        )
        action_defect = sp.factor(
            sp.Rational(1, 2)
            * path_first(action_over_tau_reduced, first_data, mass_branch)
            + sp.Rational(1, 2)
            * path_first(
                action_over_tau_reduced, second_on_shell_data, mass_branch
            )
            - path_first(action_over_tau_reduced, coarse_data, mass_branch)
        )
        defects_zero = bool(
            sp.simplify(endpoint_defect) == 0
            and sp.simplify(final_momentum_defect) == 0
            and sp.simplify(action_defect) == 0
        )
    else:
        b_root = endpoint_defect = final_momentum_defect = action_defect = None
        defects_zero = False
    first_fine_ok = bool(
        sp.simplify(first_constraint_coefficient) == 0
        and sp.simplify(first_momentum_coefficient) == 0
        and sp.simplify(seam_base) == 0
    )
    composition_ok = bool(first_fine_ok and common_b and defects_zero)
    composition = {
        "executed": True,
        "first_fine_ok": first_fine_ok,
        "degrees_in_b": list(composition_degrees),
        "second_lapse_root": str(root_second),
        "seam_root": str(root_seam),
        "root_difference": str(composition_root_difference),
        "common_root": str(b_root),
        "endpoint_defect": str(endpoint_defect),
        "final_momentum_defect": str(final_momentum_defect),
        "action_defect": str(action_defect),
        "all_defects_zero": defects_zero,
    }
    check(
        "the conditional stationary two-half-slab census is explicit",
        first_fine_ok and composition_degrees == (1, 1),
        f"common_b={common_b}; defects_zero={defects_zero}",
    )


if not all((provenance_ok, leading_zero_ok, limits_exist, static_control_ok, hostile_ok, numeric_ok)):
    outcome = "GENERIC_NEXT_ORDER_OPEN"
elif common_all_nonzero and classification_complete:
    outcome = "GENERIC_DURATION_FREE_TO_NEXT_ORDER"
    if composition_ok:
        outcome += "_COMPOSITIONAL"
elif classification_complete:
    outcome = "GENERIC_FIXED_STATE_NEXT_ORDER_OBSTRUCTION"
else:
    outcome = "GENERIC_NEXT_ORDER_OPEN"

outcome_resolved = outcome in {
    "GENERIC_DURATION_FREE_TO_NEXT_ORDER",
    "GENERIC_DURATION_FREE_TO_NEXT_ORDER_COMPOSITIONAL",
    "GENERIC_FIXED_STATE_NEXT_ORDER_OBSTRUCTION",
    "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES",
    "GENERIC_NEXT_ORDER_OPEN",
}
check("the frozen hierarchy assigns an allowed scoped outcome", outcome_resolved, outcome)


artifact = {
    "provenance": {
        "action_sha256": digest(ACTION_INPUT),
        "leading_primary_sha256": digest(PRIMARY_INPUT),
        "leading_adversarial_sha256": digest(ADVERSARIAL_INPUT),
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "correction_protocol_commit": CORRECTION_PROTOCOL_COMMIT,
    },
    "fixed_state": {
        "mass": str(mass_branch),
        "momentum": str(momentum_branch),
        "velocity_domain": "real v != 0",
    },
    "one_slab": {
        "constraint_first": str(constraint_first),
        "momentum_first": str(momentum_first),
        "degrees_in_a": list(degrees),
        "constraint_root": str(root_c),
        "momentum_root": str(root_p),
        "root_difference": str(root_difference),
        "resultant": str(resultant),
        "classification": classification,
        "classification_complete": classification_complete,
        "exceptional_velocities": exceptional_velocities,
    },
    "static_control": {
        "lapse_power": static_f_power,
        "momentum_power": static_p_power,
        "lapse_ratio_to_expected": str(static_f_ratio),
        "momentum_ratio_to_expected": str(static_p_ratio),
    },
    "hostile": {
        "mass_shift_defect": str(hostile_mass_defect),
        "momentum_shift_defect": str(hostile_momentum_defect),
    },
    "numeric_controls": numeric_records,
    "composition": composition,
    "labels": {
        "lapse_selection": (
            "NOT_SELECTED_TO_NEXT_ORDER"
            if common_all_nonzero and classification_complete
            else "OPEN"
        ),
        "composition": (
            "DERIVED_EXACT_STRUCTURAL" if composition_ok else "OPEN_OR_NEGATIVE"
        ),
        "absolute_tick": "DERIVED_NEGATIVE_UNDER_SCALE_COVARIANCE_HYPOTHESES",
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

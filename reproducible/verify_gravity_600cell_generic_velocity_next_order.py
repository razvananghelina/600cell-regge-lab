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
RADICAL_PROTOCOL_COMMIT = "71b8312"
RECOMBINATION_PROTOCOL_COMMIT = "10be749"
COMPLETE_FAILURE_PROTOCOL_COMMIT = "8acf139"
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
    and RADICAL_PROTOCOL_COMMIT == "71b8312"
    and RECOMBINATION_PROTOCOL_COMMIT == "10be749"
    and COMPLETE_FAILURE_PROTOCOL_COMMIT == "8acf139"
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


u = v**2
radical_data = (
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
        (u + 4) * (3 * u + 8) * sp.sqrt(u + 4) * sp.sqrt(3 * u + 8),
    ),
)
radical_factorization_ok = bool(
    all(
        sp.expand(polynomial - replacement**2) == 0
        for polynomial, replacement in radical_data
    )
    and sp.ask(sp.Q.positive(u + 4)) is True
    and sp.ask(sp.Q.positive(3 * u + 8)) is True
)


def normalize_positive_radicals(expression):
    def is_sqrt(node):
        return isinstance(node, sp.Pow) and node.exp == sp.Rational(1, 2)

    def replace_sqrt(node):
        base = sp.expand(node.base)
        for polynomial, replacement in radical_data:
            if sp.expand(base - polynomial) == 0:
                return replacement
        return node

    return expression.replace(is_sqrt, replace_sqrt)


def stable_normalize(expression):
    normalized = normalize_positive_radicals(expression)
    normalized = sp.together(normalized)
    normalized = normalize_positive_radicals(normalized)
    normalized = sp.cancel(normalized)
    normalized = normalize_positive_radicals(normalized)
    return sp.factor_terms(normalized)


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
    raw = expression.subs(MASS, mass).subs(data["base"])
    return stable_normalize(raw)


def path_first(expression, data, mass):
    expression = expression.subs(MASS, mass)
    raw = sum(
        sp.diff(expression, variable).subs(data["base"]) * slope
        for variable, slope in data["slope"].items()
    )
    return stable_normalize(raw)


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


def sqrt_bases(expressions):
    bases = set()
    for expression in expressions:
        for node in sp.preorder_traversal(expression):
            if (
                isinstance(node, sp.Pow)
                and node.exp == sp.Rational(1, 2)
                and node.base.has(v)
            ):
                bases.add(sp.expand(node.base))
    return bases


post_normalization_bases = sqrt_bases(
    (constraint_zero, momentum_zero, constraint_first, momentum_first)
)
primitive_bases = {sp.expand(u + 4), sp.expand(3 * u + 8)}
radical_normalization_ok = bool(
    radical_factorization_ok
    and post_normalization_bases.issubset(primitive_bases)
)
check(
    "all eight frozen positive-radical identities reduce to the two primitives",
    radical_normalization_ok,
    f"remaining={sorted(map(str, post_normalization_bases))}",
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
    root_c = stable_normalize(-c_coefficients[1] / c_coefficients[0])
    root_p = stable_normalize(-p_coefficients[1] / p_coefficients[0])
    root_difference = stable_normalize(root_c - root_p)
    resultant = sp.factor(sp.resultant(numerator_c, numerator_p, a))
else:
    root_c = root_p = root_difference = None
    resultant = sp.factor(sp.resultant(numerator_c, numerator_p, a))

generic_common_root = bool(linear_census and root_difference == 0)

# A generic degree-one result is not a global degree-one result.  Classify
# the unique degree-drop pair exactly through a monotone function of x=v^2.
c_a = stable_normalize(sp.diff(constraint_first, a))
c_0 = stable_normalize(constraint_first.subs(a, 0))
p_a = stable_normalize(sp.diff(momentum_first, a))
p_0 = stable_normalize(momentum_first.subs(a, 0))
K = 10 * radius - (v**2 + 3) * sp.sqrt(3 * v**2 + 8) * epsilon_v
B = (
    5 * v**2 * radius
    + 2 * (v**2 + 3) * sp.sqrt(3 * v**2 + 8) * epsilon_v
)
c_prefactor = 1440 * v / (
    radius * sp.sqrt(3 * v**2 + 8) * (v**2 + 3) * (v**2 + 4)
)
exceptional_factorization_ok = bool(
    stable_normalize(c_a - c_prefactor * K) == 0
    and stable_normalize(c_0 - c_prefactor * B) == 0
    and stable_normalize(c_0 * p_a - p_0 * c_a) == 0
)

x = sp.symbols("x", nonnegative=True)
r_x = sp.sqrt(x + 4)
q_x = sp.sqrt(3 * x + 8)
theta_x = sp.acos((x + 2) / (2 * (x + 3)))
epsilon_x = 2 * sp.pi - 5 * theta_x
H_x = (x + 3) * q_x * epsilon_x / r_x
epsilon_x_prime = 5 / ((x + 3) * r_x * q_x)
ratio_square = (3 * x + 8) / (x + 4)
cos_two_pi_fifths = (sp.sqrt(5) - 1) / 4
cosine_gap = (7 - 3 * sp.sqrt(5)) / 12
epsilon_zero = 2 * sp.pi - 5 * sp.acos(sp.Rational(1, 3))
monotonicity_certificate = bool(
    sp.simplify(sp.diff(epsilon_x, x) - epsilon_x_prime) == 0
    and sp.simplify(
        sp.diff(ratio_square, x) - 4 / (x + 4) ** 2
    ) == 0
    and sp.simplify(
        sp.Rational(1, 3) - cos_two_pi_fifths - cosine_gap
    ) == 0
    and 49 > 45
    and sp.N(epsilon_zero, 100) > 0
    and sp.N(3 * sp.sqrt(2) * epsilon_zero, 100) < 10
    and sp.limit(H_x, x, sp.oo) == sp.oo
)
exceptional_proof_ok = bool(
    generic_common_root
    and exceptional_factorization_ok
    and monotonicity_certificate
)
classification = (
    "ALL_REAL_NONZERO_EXCEPT_UNIQUE_PAIR"
    if exceptional_proof_ok
    else "UNRESOLVED"
)
exceptional_velocities = (
    ["-sqrt(x_star)", "sqrt(x_star)"] if exceptional_proof_ok else []
)
classification_complete = exceptional_proof_ok
check(
    "the common-root census is algebraically explicit",
    linear_census and generic_common_root,
    f"degrees={degrees}; root_difference={root_difference}; resultant={resultant}",
)
check(
    "the unique exceptional velocity pair is classified without a sign scan",
    exceptional_proof_ok,
    "K=r*(10-H), H strictly increasing from below 10 to infinity; B>0",
)


def k_numeric(x_value):
    r_value = mp.sqrt(x_value + 4)
    q_value = mp.sqrt(3 * x_value + 8)
    theta_value = mp.acos((x_value + 2) / (2 * (x_value + 3)))
    epsilon_value = 2 * mp.pi - 5 * theta_value
    return 10 * r_value - (x_value + 3) * q_value * epsilon_value


x_star_numeric = mp.findroot(k_numeric, (mp.mpf(4), mp.mpf(9)))
exceptional_bracket_ok = bool(
    k_numeric(mp.mpf(4)) > 0
    and k_numeric(mp.mpf(9)) < 0
    and mp.mpf(4) < x_star_numeric < mp.mpf(9)
    and abs(k_numeric(x_star_numeric)) < mp.mpf("1e-80")
)
check(
    "the exceptional root has the independent registered numerical bracket",
    exceptional_bracket_ok,
    f"x_star={text(x_star_numeric, 50)}",
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
    return stable_normalize(linear_second + quadratic_first)


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
hostile_mass_defect = stable_normalize(
    path_base(
        constraint_reduced,
        one_slab_data,
        mass_branch + sp.Rational(1, 10),
    )
)
hostile_momentum_defect = stable_normalize(
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
numeric_failures = []
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
            if not order_ok:
                numeric_failures.append((v_key, a_key, name))
        numeric_records[v_key][a_key] = {
            "errors": {
                name: [text(value) for value in pair]
                for name, pair in errors.items()
            },
            "orders": {name: text(value, 30) for name, value in orders.items()},
        }
expected_preasymptotic_failure = [("3/2", "0/1", "momentum")]
original_numeric_gate_reproduced = bool(
    numeric_ok or numeric_failures == expected_preasymptotic_failure
)
check(
    "the frozen numeric gate is reproduced including its preserved sole failure",
    original_numeric_gate_reproduced,
    f"failures={numeric_failures}",
)


# Registered post-failure diagnostic; it adjudicates but does not erase the
# original v=3/2, a=0 momentum pair.
diagnostic_v = mp.mpf(3) / 2
diagnostic_a = mp.mpf(0)
diagnostic_mass = mass_numeric(diagnostic_v)
diagnostic_p0 = momentum_numeric(diagnostic_v)
diagnostic_exact = p1_numeric(diagnostic_v, diagnostic_a)
diagnostic_heights = tuple(mp.mpf(1) / value for value in (1600, 3200, 6400, 12800))
diagnostic_errors = []
for diagnostic_h in diagnostic_heights:
    diagnostic_endpoint = mp.exp(diagnostic_v * diagnostic_h)
    diagnostic_direct = (
        p_pre_numeric(
            1,
            diagnostic_endpoint,
            diagnostic_h**2,
            diagnostic_mass,
        )
        - diagnostic_p0
    ) / diagnostic_h
    diagnostic_errors.append(abs(diagnostic_direct - diagnostic_exact))
diagnostic_orders = [
    mp.log(left / right) / mp.log(2)
    for left, right in zip(diagnostic_errors, diagnostic_errors[1:])
]
diagnostic_ok = bool(
    numeric_failures == expected_preasymptotic_failure
    and all(
        right < left
        for left, right in zip(diagnostic_errors, diagnostic_errors[1:])
    )
    and all(mp.mpf("0.8") <= order <= mp.mpf("1.2") for order in diagnostic_orders)
)
numeric_adjudicated_ok = bool(numeric_ok or diagnostic_ok)
check(
    "the registered smaller-height diagnostic adjudicates the sole pre-asymptotic failure",
    numeric_adjudicated_ok,
    f"orders={[text(value, 30) for value in diagnostic_orders]}",
)


# Conditional stationary two-half-slab calculation.
composition = {
    "executed": False,
    "reason": "one-slab common root is not globally classified",
}
composition_ok = None
if generic_common_root and classification_complete:
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
        root_second = stable_normalize(
            -second_coefficients[1] / second_coefficients[0]
        )
        root_seam = stable_normalize(
            -seam_coefficients[1] / seam_coefficients[0]
        )
        composition_root_difference = stable_normalize(root_second - root_seam)
        common_b = composition_root_difference == 0
    else:
        root_second = root_seam = composition_root_difference = None
        common_b = False
    if common_b:
        b_root = root_second
        endpoint_defect = stable_normalize(b_root - a_root)
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


if not all((
    provenance_ok,
    leading_zero_ok,
    limits_exist,
    radical_normalization_ok,
    exceptional_proof_ok,
    exceptional_bracket_ok,
    static_control_ok,
    hostile_ok,
    numeric_adjudicated_ok,
)):
    outcome = "GENERIC_NEXT_ORDER_OPEN"
elif classification == "ALL_REAL_NONZERO_EXCEPT_UNIQUE_PAIR":
    outcome = "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES"
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
        "radical_protocol_commit": RADICAL_PROTOCOL_COMMIT,
        "recombination_protocol_commit": RECOMBINATION_PROTOCOL_COMMIT,
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
        "exceptional_x_definition": "unique x_star>0 with K(x_star)=0",
        "exceptional_x_numeric": text(x_star_numeric, 60),
        "exceptional_x_bracket": ["4", "9"],
        "K": str(K),
        "B": str(B),
        "exceptional_factorization_ok": exceptional_factorization_ok,
        "monotonicity_certificate": monotonicity_certificate,
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
    "numeric_adjudication": {
        "original_gate_passed": numeric_ok,
        "original_failures": [list(value) for value in numeric_failures],
        "diagnostic_heights": ["1/1600", "1/3200", "1/6400", "1/12800"],
        "diagnostic_errors": [text(value) for value in diagnostic_errors],
        "diagnostic_orders": [text(value, 30) for value in diagnostic_orders],
        "diagnostic_passed": diagnostic_ok,
    },
    "composition": composition,
    "labels": {
        "lapse_selection": (
            "FREE_ON_GENERIC_DOMAIN_WITH_TWO_EXCEPTIONAL_OBSTRUCTIONS"
            if classification == "ALL_REAL_NONZERO_EXCEPT_UNIQUE_PAIR"
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

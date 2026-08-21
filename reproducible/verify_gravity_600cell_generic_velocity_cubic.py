#!/usr/bin/env python3
"""Exact cubic endpoint-jet formal-integrability gate for the 600-cell slab."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
ACTION_INPUT = HERE / "gravity_600cell_homothetic_frustum_action.json"
NEXT_ORDER_INPUT = HERE / "gravity_600cell_generic_velocity_next_order.json"
NEXT_ORDER_ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_generic_velocity_next_order_adversarial.json"
)
OUTPUT = HERE / "gravity_600cell_generic_velocity_cubic.json"

ACTION_SHA256 = (
    "c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d"
)
NEXT_ORDER_SHA256 = (
    "4bc69490fc83a193b6ac2cbd8dbe291415a13b60e4dbcce4f499bf70152e5b18"
)
NEXT_ORDER_ADVERSARIAL_SHA256 = (
    "3ab16e6d19b527590b3dce6e8b3caa093efb6cc504a2a7824362ffc529a83a05"
)
PRIOR_ART_COMMIT = "70e7ca2"
PROTOCOL_COMMIT = "d2efdf4"
CLASSIFICATION_PROTOCOL_COMMIT = "3de74e7"

VELOCITY_RATIONALS = ((1, 2), (3, 2), (3, 1))
CUBIC_RATIONALS = ((0, 1), (1, 11))
HEIGHT_DENOMINATORS = (1000, 2000, 4000, 8000)
mp.mp.dps = 100

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


action_artifact = json.loads(ACTION_INPUT.read_text())
next_order_artifact = json.loads(NEXT_ORDER_INPUT.read_text())
next_order_adversarial_artifact = json.loads(
    NEXT_ORDER_ADVERSARIAL_INPUT.read_text()
)
provenance_ok = bool(
    digest(ACTION_INPUT) == ACTION_SHA256
    and digest(NEXT_ORDER_INPUT) == NEXT_ORDER_SHA256
    and digest(NEXT_ORDER_ADVERSARIAL_INPUT)
    == NEXT_ORDER_ADVERSARIAL_SHA256
    and action_artifact["outcome"] == "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
    and next_order_artifact["outcome"]
    == "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES"
    and next_order_adversarial_artifact["outcome"]
    == "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES_ADVERSARIALLY_CORROBORATED"
    and PRIOR_ART_COMMIT == "70e7ca2"
    and PROTOCOL_COMMIT == "d2efdf4"
    and CLASSIFICATION_PROTOCOL_COMMIT == "3de74e7"
)
check("the action and both accepted next-order artifacts are frozen", provenance_ok)


# Reconstruct the full action for hostile and arbitrary-precision controls.
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


v, a, c = sp.symbols("v a c", real=True)
h = sp.symbols("h", positive=True)
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


# Frozen positive radical inventory from the accepted next-order route.
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
        27 * u**5
        + 432 * u**4
        + 2736 * u**3
        + 8576 * u**2
        + 13312 * u
        + 8192,
        (u + 4) * (3 * u + 8) * sp.sqrt(3 * u + 8),
    ),
    (
        27 * u**6
        + 540 * u**5
        + 4464 * u**4
        + 19520 * u**3
        + 47616 * u**2
        + 61440 * u
        + 32768,
        (u + 4)
        * (3 * u + 8)
        * sp.sqrt(u + 4)
        * sp.sqrt(3 * u + 8),
    ),
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


# Exact scaled expressions.  In the accepted convention constraint_reduced is
# 2F/h and -p_minus_reduced is p_pre.
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

variables = (lm, lp, q, w, tau)
base = {lm: 1, lp: 1, q: v, w: -1200 * sp.sqrt(3) * v, tau: 0}
slope = {
    lm: 0,
    lp: v,
    q: a + v**2 / 2,
    w: -1200 * sp.sqrt(3) * (a + v**2),
    tau: 1,
}
quadratic = {
    lm: 0,
    lp: a + v**2 / 2,
    q: c + v * a + v**3 / 6,
    w: -1200 * sp.sqrt(3) * (c + 2 * v * a + 2 * v**3 / 3),
    tau: 0,
}


def path_base_raw(expression, mass):
    return expression.subs(MASS, mass).subs(base)


def path_first_raw(expression, mass):
    expression = expression.subs(MASS, mass)
    return sum(
        sp.diff(expression, variable).subs(base) * slope[variable]
        for variable in variables
    )


def path_second_raw(expression, mass):
    expression = expression.subs(MASS, mass)
    linear_second = sum(
        sp.diff(expression, variable).subs(base) * quadratic[variable]
        for variable in variables
    )
    quadratic_first = sp.Rational(1, 2) * sum(
        sp.diff(expression, left, right).subs(base)
        * slope[left]
        * slope[right]
        for left in variables
        for right in variables
    )
    return linear_second + quadratic_first


constraint_zero = stable_normalize(path_base_raw(constraint_reduced, mass_branch))
momentum_zero = stable_normalize(
    -path_base_raw(p_minus_reduced, mass_branch) - momentum_branch
)
constraint_first = stable_normalize(path_first_raw(constraint_reduced, mass_branch))
momentum_first = stable_normalize(-path_first_raw(p_minus_reduced, mass_branch))

c1_a = stable_normalize(sp.diff(constraint_first, a))
c1_0 = stable_normalize(constraint_first.subs(a, 0))
p1_a = stable_normalize(sp.diff(momentum_first, a))
p1_0 = stable_normalize(momentum_first.subs(a, 0))
c_prefactor = 1440 * v / (
    radius * qrad * (u + 3) * (u + 4)
)
leading_controls_ok = bool(
    sp.simplify(constraint_zero) == 0
    and sp.simplify(momentum_zero) == 0
    and stable_normalize(c1_a - c_prefactor * K) == 0
    and stable_normalize(c1_0 - c_prefactor * B) == 0
    and stable_normalize(c1_0 * p1_a - p1_0 * c1_a) == 0
    and stable_normalize(constraint_first.subs(a, a_root)) == 0
    and stable_normalize(momentum_first.subs(a, a_root)) == 0
)
check(
    "the accepted zeroth and first-order branch is reconstructed exactly",
    leading_controls_ok,
)


# Derive the full generic second path coefficient first, then insert the
# already-frozen lower-order solution as required by the protocol.
constraint_second_generic = path_second_raw(constraint_reduced, mass_branch)
momentum_second_generic = -path_second_raw(p_minus_reduced, mass_branch)
constraint_second = stable_normalize(
    constraint_second_generic.subs(a, a_root)
)
momentum_second = stable_normalize(momentum_second_generic.subs(a, a_root))

coefficients_exist = bool(
    not constraint_second.has(sp.Limit)
    and not momentum_second.has(sp.Limit)
    and constraint_second.is_finite is not False
    and momentum_second.is_finite is not False
)
check(
    "the exact on-shell cubic lapse and momentum coefficients exist",
    coefficients_exist,
)


primitive_bases = {sp.expand(u + 4), sp.expand(3 * u + 8)}
remaining_sqrt_bases = sqrt_bases(
    (constraint_first, momentum_first, constraint_second, momentum_second)
)
radical_inventory_ok = remaining_sqrt_bases.issubset(primitive_bases)
check(
    "the frozen positive-radical inventory suffices at cubic order",
    radical_inventory_ok,
    f"remaining={sorted(map(str, remaining_sqrt_bases))}",
)


c2_c = stable_normalize(sp.diff(constraint_second, c))
c2_0 = stable_normalize(constraint_second.subs(c, 0))
p2_c = stable_normalize(sp.diff(momentum_second, c))
p2_0 = stable_normalize(momentum_second.subs(c, 0))
linear_in_c = bool(
    stable_normalize(sp.diff(constraint_second, c, 2)) == 0
    and stable_normalize(sp.diff(momentum_second, c, 2)) == 0
    and c2_c != 0
)
coefficient_recursion_ok = bool(
    linear_in_c
    and stable_normalize(c2_c - c1_a) == 0
    and stable_normalize(p2_c - p1_a) == 0
)
check(
    "both cubic equations are affine in c with the registered recursive slopes",
    coefficient_recursion_ok,
)


constraint_root = stable_normalize(-c2_0 / c2_c) if linear_in_c else None
momentum_root = stable_normalize(-p2_0 / p2_c) if p2_c != 0 else None
cross_resultant = stable_normalize(c2_0 * p2_c - p2_0 * c2_c)
generic_common_root = bool(
    linear_in_c
    and p2_c != 0
    and cross_resultant == 0
    and stable_normalize(constraint_root - momentum_root) == 0
)

# The first OPEN run exposed a simple exact positive cross-resultant.  Its
# classification proof was frozen in commit 3de74e7 before this logic changed.
cross_expected = 129600 * epsilon_v**2 / (u + 4)
cross_factorization_ok = bool(
    stable_normalize(cross_resultant - cross_expected) == 0
)
x = sp.symbols("x", nonnegative=True)
theta_x = sp.acos((x + 2) / (2 * (x + 3)))
epsilon_x = 2 * sp.pi - 5 * theta_x
epsilon_x_prime = 5 / (
    (x + 3) * sp.sqrt(x + 4) * sp.sqrt(3 * x + 8)
)
cos_two_pi_fifths = (sp.sqrt(5) - 1) / 4
cosine_gap = (7 - 3 * sp.sqrt(5)) / 12
epsilon_positivity_certificate = bool(
    sp.simplify(sp.diff(epsilon_x, x) - epsilon_x_prime) == 0
    and sp.simplify(
        sp.Rational(1, 3) - cos_two_pi_fifths - cosine_gap
    )
    == 0
    and 49 > 45
)

# C2's slope is nonzero on the registered domain because it is the
# positive-denominator prefactor times v*K.  All inverse-function arguments
# are real and nonsingular there.
branch_argument = (u + 2) / (2 * (u + 3))
inverse_branch_ok = bool(
    sp.simplify(branch_argument - sp.Rational(1, 3)) == u / (6 * (u + 3))
    and sp.simplify(sp.Rational(1, 2) - branch_argument)
    == sp.Rational(1, 2) / (u + 3)
)
complete_domain = bool(
    coefficient_recursion_ok
    and cross_factorization_ok
    and epsilon_positivity_certificate
    and inverse_branch_ok
)
classification = (
    "NO_COMMON_C_ON_COMPLETE_REGISTERED_DOMAIN"
    if complete_domain
    else "UNRESOLVED"
)

check(
    "the positive cross-resultant excludes a common c on the complete domain",
    bool(
        complete_domain
        and not generic_common_root
        and classification == "NO_COMMON_C_ON_COMPLETE_REGISTERED_DOMAIN"
        and cross_resultant != 0
        and c2_c != 0
    ),
    (
        f"factorization={cross_factorization_ok}; "
        f"positivity={epsilon_positivity_certificate}; "
        f"classification={classification}"
    ),
)

root_parity_data = {
    "constraint_root_under_v_to_minus_v": str(
        stable_normalize(constraint_root.subs(v, -v))
    ),
    "momentum_root_under_v_to_minus_v": str(
        stable_normalize(momentum_root.subs(v, -v))
        if momentum_root is not None
        else None
    ),
}


hostile_mass_defect = stable_normalize(
    path_base_raw(constraint_reduced, mass_branch + sp.Rational(1, 10))
)
hostile_acceleration_defect = stable_normalize(
    constraint_first.subs(a, a_root + sp.Rational(1, 10))
)
hostile_ok = bool(
    sp.simplify(hostile_mass_defect + 4 * sp.pi / 5) == 0
    and stable_normalize(hostile_acceleration_defect - c1_a / 10) == 0
    and c1_a != 0
)
check(
    "the frozen changed-mass and changed-acceleration controls fail exactly",
    hostile_ok,
    f"mass={hostile_mass_defect}; acceleration={hostile_acceleration_defect}",
)


# Direct arbitrary-precision controls are evaluated only after exact
# extraction and the symbolic census above.
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
numeric_ok = True
numeric_failures = []
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
            if all(error < mp.mpf("1e-70") for error in sequence):
                orders[name] = [mp.inf] * 3
                sequence_ok = True
            elif all(error >= mp.mpf("1e-70") for error in sequence):
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
    "all preregistered direct cubic-coefficient controls converge at first order",
    numeric_ok,
    f"failures={numeric_failures}",
)


if not all(
    (
        provenance_ok,
        leading_controls_ok,
        coefficients_exist,
        radical_inventory_ok,
        coefficient_recursion_ok,
        hostile_ok,
        numeric_ok,
    )
):
    outcome = "GENERIC_CUBIC_OPEN"
elif complete_domain and not generic_common_root:
    outcome = "GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION"
else:
    outcome = "GENERIC_CUBIC_OPEN"

allowed_outcomes = {
    "GENERIC_DURATION_FREE_TO_CUBIC_ORDER",
    "GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION",
    "GENERIC_CUBIC_EXCEPTIONAL_STRATA",
    "GENERIC_CUBIC_OPEN",
}
check(
    "the preregistered hierarchy assigns an allowed scoped outcome",
    outcome in allowed_outcomes,
    outcome,
)


artifact = {
    "provenance": {
        "action_sha256": digest(ACTION_INPUT),
        "next_order_primary_sha256": digest(NEXT_ORDER_INPUT),
        "next_order_adversarial_sha256": digest(NEXT_ORDER_ADVERSARIAL_INPUT),
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "classification_protocol_commit": CLASSIFICATION_PROTOCOL_COMMIT,
    },
    "fixed_state": {
        "mass": str(mass_branch),
        "momentum": str(momentum_branch),
        "acceleration": str(a_root),
        "domain": "real v != 0 and K(v^2) != 0",
    },
    "controls": {
        "leading_reconstructed": leading_controls_ok,
        "hostile_mass_defect": str(hostile_mass_defect),
        "hostile_acceleration_defect": str(hostile_acceleration_defect),
        "radical_inventory_sufficient": radical_inventory_ok,
        "remaining_sqrt_bases": sorted(map(str, remaining_sqrt_bases)),
    },
    "cubic_census": {
        "constraint_second": str(constraint_second),
        "momentum_second": str(momentum_second),
        "degrees_in_c": [1, 1] if coefficient_recursion_ok else None,
        "constraint_leading_coefficient": str(c2_c),
        "momentum_leading_coefficient": str(p2_c),
        "constraint_constant_coefficient": str(c2_0),
        "momentum_constant_coefficient": str(p2_0),
        "constraint_root": str(constraint_root),
        "momentum_root": str(momentum_root),
        "cross_resultant": str(cross_resultant),
        "cross_resultant_expected": str(cross_expected),
        "cross_factorization_ok": cross_factorization_ok,
        "epsilon_positivity_certificate": epsilon_positivity_certificate,
        "generic_common_root": generic_common_root,
        "classification": classification,
        "classification_complete": complete_domain,
        "inverse_branch_certificate": inverse_branch_ok,
        "root_parity_data": root_parity_data,
        "excluded_degree_drop_loci": ["v=0", "K(v^2)=0"],
    },
    "numeric_controls": numeric_records,
    "numeric_failures": [list(value) for value in numeric_failures],
    "labels": {
        "formal_integrability": (
            "DERIVED_NEGATIVE_SCOPED"
            if complete_domain and not generic_common_root
            else "OPEN"
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

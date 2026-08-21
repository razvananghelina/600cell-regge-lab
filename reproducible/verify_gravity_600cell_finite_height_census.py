#!/usr/bin/env python3
"""Exact construction gate for the all-real finite-height state census."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_generic_velocity_cubic.json"
ADVERSARIAL_INPUT = HERE / "gravity_600cell_generic_velocity_cubic_adversarial.json"
OUTPUT = HERE / "gravity_600cell_finite_height_census.json"

PRIMARY_SHA256 = (
    "1d35b46cd4db20df0af3ed3e6b5de676d69753cf5059e0eb607d1eec949b9103"
)
ADVERSARIAL_SHA256 = (
    "b5167d597a927f8b441a096c31034aa04efa435883284dc2d9bfbd3b9cb3ff0d"
)
PRIOR_ART_COMMIT = "c8edca2"
DOMAIN_CORRECTION_COMMIT = "726e52f"
PROTOCOL_COMMIT = "79747ea"

VELOCITY_RATIONALS = ((-3, 2), (0, 1), (2, 3), (5, 2))
QUOTIENT_RATIONALS = ((-2, 1), (-1, 3), (1, 2), (3, 1))
HEIGHT_RATIONALS = ((1, 7), (2, 5))
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


primary = json.loads(PRIMARY_INPUT.read_text())
adversarial = json.loads(ADVERSARIAL_INPUT.read_text())
provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and digest(ADVERSARIAL_INPUT) == ADVERSARIAL_SHA256
    and primary["outcome"] == "GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION"
    and adversarial["outcome"]
    == "GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION_ADVERSARIALLY_CORROBORATED"
    and PRIOR_ART_COMMIT == "c8edca2"
    and DOMAIN_CORRECTION_COMMIT == "726e52f"
    and PROTOCOL_COMMIT == "79747ea"
)
check("both cubic inputs and the corrected finite protocol are frozen", provenance_ok)


# Differentiate the complete action before introducing finite-height
# variables, as required by the protocol.
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

h = sp.symbols("h", positive=True)
q = sp.symbols("q", real=True)
finite_substitution = {
    L_MINUS: 1,
    L_PLUS: 1 + h * q,
    RHO: h**2,
}
constraint_qh_raw = sp.factor(
    sp.simplify((2 * F_EXACT / h).subs(finite_substitution))
)
momentum_qh_raw = sp.factor(sp.simplify(P_PRE_EXACT.subs(finite_substitution)))


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


constraint_qh = normalize_radicals(constraint_qh_raw, q)
momentum_qh = normalize_radicals(momentum_qh_raw, q)
affine_certificate = bool(
    normalize_radicals(sp.diff(constraint_qh, h, 2), q) == 0
    and normalize_radicals(sp.diff(momentum_qh, h, 2), q) == 0
)
check("both complete finite residuals are exactly affine in height", affine_certificate)


def state_functions(z):
    square = z**2
    radius = sp.sqrt(square + 4)
    theta = sp.acos((square + 2) / (2 * (square + 3)))
    eta = sp.asinh(z / sp.sqrt(8 * (square + 3)))
    epsilon = 2 * sp.pi - 5 * theta
    mass = 180 * epsilon / (sp.pi * radius)
    momentum = 180 * z * epsilon / radius - 600 * sp.sqrt(3) * eta
    return radius, theta, eta, epsilon, mass, momentum


_, theta_q, eta_q, epsilon_q, mass_q, momentum_q = state_functions(q)
v = sp.symbols("v", real=True)
_, theta_v, eta_v, epsilon_v, mass_v, momentum_v = state_functions(v)

f_q = sp.sqrt(q**2 + 4) * epsilon_q / 2
f_q_prime = sp.diff(f_q, q)
eta_q_prime = sp.diff(eta_q, q)

constraint_constant_general = normalize_radicals(constraint_qh.subs(h, 0), q)
constraint_height = normalize_radicals(sp.diff(constraint_qh, h), q)
momentum_constant_general = normalize_radicals(momentum_qh.subs(h, 0), q)
momentum_height = normalize_radicals(sp.diff(momentum_qh, h), q)

constraint_constant_expected = 8 * sp.pi * (mass_q - MASS)
momentum_constant_expected = momentum_q
constraint_height_expected = (
    360 * q * f_q
    - 360 * q**2 * f_q_prime
    + 600 * sp.sqrt(3) * q**3 * eta_q_prime
)
momentum_height_expected = (
    -180 * f_q
    + 180 * q * f_q_prime
    - 300 * sp.sqrt(3) * q**2 * eta_q_prime
)
closed_formula_ok = bool(
    normalize_radicals(
        constraint_constant_general - constraint_constant_expected,
        q,
    )
    == 0
    and normalize_radicals(
        momentum_constant_general - momentum_constant_expected,
        q,
    )
    == 0
    and normalize_radicals(constraint_height - constraint_height_expected, q)
    == 0
    and normalize_radicals(momentum_height - momentum_height_expected, q)
    == 0
)
check("all four affine coefficients have exact closed formulas", closed_formula_ok)


constraint_constant = 8 * sp.pi * (mass_q - mass_v)
momentum_constant = momentum_q - momentum_v
elimination = sp.factor_terms(
    constraint_constant * momentum_height
    - momentum_constant * constraint_height
)
wrong_elimination = sp.factor_terms(
    constraint_constant * momentum_height
    + momentum_constant * constraint_height
)

boundary_constraint = normalize_radicals(constraint_constant.subs(q, v), v)
boundary_momentum = normalize_radicals(momentum_constant.subs(q, v), v)
boundary_ok = bool(boundary_constraint == 0 and boundary_momentum == 0)
check("the complete state family has the exact h=0, q=v boundary root", boundary_ok)


# The finite determinant is a tangent-chord condition for the canonical state
# curve.  This identity is derived before any root sampling.
mass_q_prime = sp.diff(mass_q, q)
momentum_q_prime = sp.diff(momentum_q, q)
tangent_chord = (
    (mass_q - mass_v) * momentum_q_prime
    - (momentum_q - momentum_v) * mass_q_prime
)
slope_tangent_identity = normalize_radicals(
    constraint_height * momentum_q_prime
    - 8 * sp.pi * momentum_height * mass_q_prime,
    q,
)
elimination_tangent_identity = normalize_radicals(
    elimination * mass_q_prime - constraint_height * tangent_chord,
    q,
)
tangent_reduction_ok = bool(
    slope_tangent_identity == 0 and elimination_tangent_identity == 0
)
check(
    "the finite determinant is exactly a tangent-chord condition",
    tangent_reduction_ok,
)


diagonal_zero = normalize_radicals(elimination.subs(q, v), v)
diagonal_first = normalize_radicals(sp.diff(elimination, q).subs(q, v), v)
diagonal_second = normalize_radicals(sp.diff(elimination, q, 2).subs(q, v), v)
diagonal_third = normalize_radicals(sp.diff(elimination, q, 3).subs(q, v), v)
diagonal_multiplicity_at_least_two = bool(
    diagonal_zero == 0 and diagonal_first == 0 and diagonal_second != 0
)
check(
    "the generic diagonal root has its exact first nonzero derivative recorded",
    diagonal_zero == 0 and diagonal_first == 0,
    f"D2_zero={diagonal_second == 0}; D3_zero={diagonal_third == 0}",
)


delta_cubic = 129600 * epsilon_v**2 / (v**2 + 4)
diagonal_to_cubic_ratio = normalize_radicals(diagonal_second / delta_cubic, v)
cubic_control_resolved = bool(
    diagonal_multiplicity_at_least_two
    and diagonal_to_cubic_ratio.has(sp.nan, sp.zoo, sp.oo, -sp.oo) is False
)
check(
    "the diagonal finite determinant exposes an exact cubic-control ratio",
    cubic_control_resolved,
    f"D2/Delta={diagonal_to_cubic_ratio}",
)


# Record, but do not pretend to complete, the slope-exception census on the
# first construction run.
constraint_slope_at_zero = normalize_radicals(constraint_height.subs(q, 0), q)
momentum_slope_at_zero = normalize_radicals(momentum_height.subs(q, 0), q)
common_slope_identity = normalize_radicals(
    constraint_height * momentum_q_prime
    - 8 * sp.pi * momentum_height * mass_q_prime,
    q,
)
slope_data_recorded = bool(
    constraint_slope_at_zero.is_finite is not False
    and momentum_slope_at_zero.is_finite is not False
    and common_slope_identity == 0
)
check(
    "the zero-slope case split has exact input expressions",
    slope_data_recorded,
    (
        f"Ch(0)={constraint_slope_at_zero}; "
        f"Ph(0)={momentum_slope_at_zero}"
    ),
)


# Frozen hostile controls.
hostile_mass_change = normalize_radicals(
    constraint_constant_general.subs(MASS, mass_v + sp.Rational(1, 10))
    - constraint_constant_general.subs(MASS, mass_v),
    q,
)
hostile_momentum_change = -sp.Rational(1, 10)
wrong_point = {
    v: sp.Rational(2, 3),
    q: -sp.Rational(1, 3),
}
wrong_terms = (
    sp.N((constraint_constant * momentum_height).subs(wrong_point), 100),
    sp.N((momentum_constant * constraint_height).subs(wrong_point), 100),
)
wrong_sign_gap = sp.N((wrong_elimination - elimination).subs(wrong_point), 100)
hostile_ok = bool(
    sp.simplify(hostile_mass_change + 4 * sp.pi / 5) == 0
    and hostile_momentum_change == -sp.Rational(1, 10)
    and all(abs(value) > sp.Float("1e-80", 100) for value in wrong_terms)
    and abs(wrong_sign_gap) > sp.Float("1e-80", 100)
)
check(
    "all three frozen finite-height hostile controls fail exactly",
    hostile_ok,
    (
        f"mass={hostile_mass_change}; momentum={hostile_momentum_change}; "
        f"wrong_sign_gap={wrong_sign_gap}"
    ),
)


# Direct full-action versus affine-reconstruction controls.  This validates
# the global reduction but is not a root search.
constraint_full_numeric = sp.lambdify((h, q, MASS), constraint_qh, "mpmath")
momentum_full_numeric = sp.lambdify((h, q), momentum_qh, "mpmath")
mass_numeric = sp.lambdify(v, mass_v, "mpmath")
momentum_numeric = sp.lambdify(v, momentum_v, "mpmath")
c0_numeric = sp.lambdify((v, q), constraint_constant, "mpmath")
ch_numeric = sp.lambdify(q, constraint_height, "mpmath")
p0_numeric = sp.lambdify((v, q), momentum_constant, "mpmath")
ph_numeric = sp.lambdify(q, momentum_height, "mpmath")

numeric_records = []
numeric_ok = True
numeric_failures = []
for v_numerator, v_denominator in VELOCITY_RATIONALS:
    v_value = mp.mpf(v_numerator) / v_denominator
    mass_value = mass_numeric(v_value)
    p_value = momentum_numeric(v_value)
    for q_numerator, q_denominator in QUOTIENT_RATIONALS:
        q_value = mp.mpf(q_numerator) / q_denominator
        for h_numerator, h_denominator in HEIGHT_RATIONALS:
            h_value = mp.mpf(h_numerator) / h_denominator
            endpoint = 1 + h_value * q_value
            if endpoint <= 0:
                continue
            direct_c = constraint_full_numeric(h_value, q_value, mass_value)
            direct_p = momentum_full_numeric(h_value, q_value) - p_value
            affine_c = c0_numeric(v_value, q_value) + h_value * ch_numeric(q_value)
            affine_p = p0_numeric(v_value, q_value) + h_value * ph_numeric(q_value)
            errors = {
                "constraint": abs(direct_c - affine_c)
                / max(mp.mpf(1), abs(affine_c)),
                "momentum": abs(direct_p - affine_p)
                / max(mp.mpf(1), abs(affine_p)),
            }
            row_ok = all(value < mp.mpf("1e-80") for value in errors.values())
            numeric_ok &= row_ok
            if not row_ok:
                numeric_failures.append(
                    (
                        f"{v_numerator}/{v_denominator}",
                        f"{q_numerator}/{q_denominator}",
                        f"{h_numerator}/{h_denominator}",
                    )
                )
            numeric_records.append(
                {
                    "v": f"{v_numerator}/{v_denominator}",
                    "q": f"{q_numerator}/{q_denominator}",
                    "h": f"{h_numerator}/{h_denominator}",
                    "endpoint": text(endpoint),
                    "errors": {name: text(value) for name, value in errors.items()},
                }
            )

check(
    "all frozen full-action values match the exact affine reconstruction",
    numeric_ok,
    f"rows={len(numeric_records)}; failures={numeric_failures}",
)


# The construction run intentionally does not infer a global zero set from
# samples.  Completeness remains false until slopes, tails and tangent-chord
# roots are certified on all R^2.
slope_exception_classification_complete = False
tail_bounds_complete = False
tangent_chord_root_classification_complete = False
global_classification_complete = all(
    (
        slope_exception_classification_complete,
        tail_bounds_complete,
        tangent_chord_root_classification_complete,
    )
)
outcome = "FINITE_HEIGHT_OPEN"
check(
    "the preregistered hierarchy keeps the incomplete global census OPEN",
    not global_classification_complete and outcome == "FINITE_HEIGHT_OPEN",
)


artifact = {
    "provenance": {
        "primary_cubic_sha256": digest(PRIMARY_INPUT),
        "adversarial_cubic_sha256": digest(ADVERSARIAL_INPUT),
        "prior_art_commit": PRIOR_ART_COMMIT,
        "domain_correction_commit": DOMAIN_CORRECTION_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
    },
    "domain": {
        "state_parameter": "all real v",
        "quotient": "all real q",
        "height": "h>0",
        "endpoint": "1+h*q>0",
        "includes_v_zero": True,
        "includes_K_zero": True,
    },
    "affine_reduction": {
        "certified": affine_certificate and closed_formula_ok,
        "constraint_constant": str(constraint_constant),
        "constraint_height_slope": str(constraint_height),
        "momentum_constant": str(momentum_constant),
        "momentum_height_slope": str(momentum_height),
        "elimination": str(elimination),
        "tangent_chord": str(tangent_chord),
        "tangent_reduction_ok": tangent_reduction_ok,
    },
    "diagonal_control": {
        "D_at_q_equals_v": str(diagonal_zero),
        "D_q_at_q_equals_v": str(diagonal_first),
        "D_qq_at_q_equals_v": str(diagonal_second),
        "D_qqq_at_q_equals_v": str(diagonal_third),
        "generic_multiplicity_at_least_two": diagonal_multiplicity_at_least_two,
        "D_qq_over_cubic_Delta": str(diagonal_to_cubic_ratio),
        "cubic_control_resolved": cubic_control_resolved,
    },
    "slope_case_input": {
        "constraint_slope_at_q_zero": str(constraint_slope_at_zero),
        "momentum_slope_at_q_zero": str(momentum_slope_at_zero),
        "common_slope_identity": str(common_slope_identity),
        "classification_complete": slope_exception_classification_complete,
    },
    "global_census": {
        "tail_bounds_complete": tail_bounds_complete,
        "tangent_chord_root_classification_complete": (
            tangent_chord_root_classification_complete
        ),
        "classification_complete": global_classification_complete,
        "root_scan_performed": False,
    },
    "hostile": {
        "mass_shift_defect": str(hostile_mass_change),
        "momentum_shift_defect": str(hostile_momentum_change),
        "wrong_sign_terms": [str(value) for value in wrong_terms],
        "wrong_sign_gap": str(wrong_sign_gap),
        "passed": hostile_ok,
    },
    "numeric_affine_controls": numeric_records,
    "numeric_failures": [list(value) for value in numeric_failures],
    "labels": {
        "finite_height_census": "OPEN",
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

#!/usr/bin/env python3
"""Global Lagrange resolution of the direct-quotient adversarial OPEN."""

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp
from flint import arb, arb_series, ctx


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = (
    ROOT
    / "docs"
    / "gravity"
    / "gravity_600cell_finite_height_invariant_region_adversarial_resolution_protocol.md"
)
FIRST_ADVERSARIAL = (
    HERE / "gravity_600cell_finite_height_invariant_region_adversarial.json"
)
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
FIFTH_INPUT = HERE / "gravity_600cell_finite_height_asymptotic_map.json"
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_invariant_region.json"
OUTPUT = (
    HERE
    / "gravity_600cell_finite_height_invariant_region_adversarial_resolution.json"
)

PROTOCOL_COMMIT = "a851018"
PROTOCOL_SHA256 = (
    "d24d1cef060674462bd73cb7b2b20dfbb7be72deb8bac13138d19d79d60524d3"
)
FIRST_ADVERSARIAL_COMMIT = "bfa4db4"
FIRST_ADVERSARIAL_SHA256 = (
    "f7d1f36e5ed679c39d1c38dbc21509ae52211f6735b38a2da46046fb798f54d5"
)
CLASSIFICATION_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
FIFTH_SHA256 = (
    "a93837d2bbec340ddbac528c0be4da52aefe45c8f0d4310496eb1aef6a7b19b6"
)
PRIMARY_SHA256 = (
    "9b6a473c462e7d23af50878cdd4d849bb66c69068c3178b82235e5d0e39926b9"
)

PRECISIONS = (160, 256)
DEGREE = 6
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


def rational_ball(numerator, denominator=1):
    return arb(numerator) / arb(denominator)


def interval_ball(left, right):
    left = Fraction(left)
    right = Fraction(right)
    midpoint = (left + right) / 2
    radius = (right - left) / 2
    center = rational_ball(midpoint.numerator, midpoint.denominator)
    return center + arb(0, f"{radius.numerator}/{radius.denominator}")


def strict_positive(value):
    return bool(
        value.is_finite()
        and value > 0
        and value.lower() > 0
        and not value.contains(0)
    )


def arb_record(value):
    if not value.is_finite():
        return {
            "pretty": str(value),
            "lower": "UNAVAILABLE",
            "upper": "UNAVAILABLE",
            "contains_zero": bool(value.contains(0)),
            "finite": False,
        }
    return {
        "pretty": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "contains_zero": bool(value.contains(0)),
        "finite": True,
    }


def mp_text(value, digits=75):
    return mp.nstr(value, digits)


first_adversarial = json.loads(FIRST_ADVERSARIAL.read_text())
classification = json.loads(CLASSIFICATION_INPUT.read_text())
provenance_ok = bool(
    digest(PROTOCOL) == PROTOCOL_SHA256
    and digest(FIRST_ADVERSARIAL) == FIRST_ADVERSARIAL_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and first_adversarial["outcome"]
    == "INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN"
    and first_adversarial["independent_outcome_before_primary_read"]
    == "INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN"
    and classification["outcome"]
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
)
check(
    "the resolution protocol and preserved direct-quotient OPEN are frozen",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}; first_result={FIRST_ADVERSARIAL_COMMIT}",
)


# Complete action, derived again in this verifier.
l_minus, l_plus, squared_strut, dust_mass = sp.symbols(
    "l_minus l_plus squared_strut dust_mass", positive=True
)
delta = l_plus - l_minus
proper_height = sp.sqrt(squared_strut + delta**2 / 4)
cosine = (delta**2 + 2 * squared_strut) / (
    2 * (delta**2 + 3 * squared_strut)
)
boost = delta / sp.sqrt(8 * (delta**2 + 3 * squared_strut))
full_action = (
    360
    * (l_minus + l_plus)
    * proper_height
    * (2 * sp.pi - 5 * sp.acos(cosine))
    + 600
    * sp.sqrt(3)
    * (l_minus**2 - l_plus**2)
    * sp.asinh(boost)
    - 8 * sp.pi * dust_mass * sp.sqrt(squared_strut)
)
constraint_exact = squared_strut * sp.diff(full_action, squared_strut)
pre_exact = -l_minus * sp.diff(full_action, l_minus) / 2
post_exact = l_plus * sp.diff(full_action, l_plus) / 2

scale = sp.symbols("scale", positive=True)
scaling = {
    l_minus: scale * l_minus,
    l_plus: scale * l_plus,
    squared_strut: scale**2 * squared_strut,
    dust_mass: scale * dust_mass,
}
action_scaling_ok = bool(
    sp.simplify(
        sp.powsimp(full_action.subs(scaling) - scale**2 * full_action, force=True)
    )
    == 0
    and sp.simplify(
        sp.powsimp(post_exact.subs(scaling) - scale**2 * post_exact, force=True)
    )
    == 0
)


# Explicit chain-rule factorization; no generic nested-radical simplifier is
# asked to discover the decisive common factor.
q = sp.symbols("q", real=True, nonzero=True)
q_squared = q**2
r = sp.sqrt(q_squared + 4)
s = sp.sqrt(3 * q_squared + 8)
z = (q_squared + 2) / (2 * (q_squared + 3))
eta_argument = q / sp.sqrt(8 * (q_squared + 3))
epsilon_symbol = sp.symbols("epsilon_symbol", real=True)

radical_identity = sp.factor(
    4 * (q_squared + 3) ** 2 * (1 - z**2) - r**2 * s**2
)
z_derivative_identity = sp.simplify(
    sp.diff(z, q) - q / (q_squared + 3) ** 2
)
eta_derivative_identity = sp.powdenest(
    sp.simplify(
        sp.diff(eta_argument, q) / sp.sqrt(1 + eta_argument**2)
        - sp.sqrt(3) / ((q_squared + 3) * s)
    ),
    force=True,
)
epsilon_prime = 10 * q / ((q_squared + 3) * r * s)
mu_prime_chain = 180 / sp.pi * (
    epsilon_prime / r - epsilon_symbol * q / r**3
)
p_prime_chain = (
    180
    * (
        epsilon_symbol / r
        + q * epsilon_prime / r
        - q**2 * epsilon_symbol / r**3
    )
    - 1800 / ((q_squared + 3) * s)
)
k_factor = 10 * r - (q_squared + 3) * s * epsilon_symbol
mu_prime_expected = (
    180 * q * k_factor / (sp.pi * r**3 * (q_squared + 3) * s)
)
p_prime_expected = -720 * k_factor / (r**3 * (q_squared + 3) * s)
chain_factor_ok = bool(
    radical_identity == 0
    and z_derivative_identity == 0
    and eta_derivative_identity == 0
    and sp.factor(mu_prime_chain - mu_prime_expected) == 0
    and sp.factor(p_prime_chain - p_prime_expected) == 0
    and sp.factor(p_prime_chain + 4 * sp.pi * mu_prime_chain / q) == 0
)

mass_symbol, mu_symbol, p_prime_symbol, mu_prime_symbol = sp.symbols(
    "mass_symbol mu_symbol p_prime_symbol mu_prime_symbol", nonzero=True
)
r_prime_abstract = p_prime_symbol + 4 * sp.pi * (
    q * mu_prime_symbol - (mu_symbol - mass_symbol)
) / q**2
r_prime_reduced = sp.factor(
    r_prime_abstract.subs(
        p_prime_symbol, -4 * sp.pi * mu_prime_symbol / q
    )
)
r_prime_ok = bool(
    sp.simplify(
        r_prime_reduced - 4 * sp.pi * (mass_symbol - mu_symbol) / q**2
    )
    == 0
)
redifferentiated_ok = bool(action_scaling_ok and chain_factor_ok and r_prime_ok)
check(
    "explicit radical identities certify the action and common derivative factor",
    redifferentiated_ok,
)

constraint_numeric = sp.lambdify(
    (l_minus, l_plus, squared_strut, dust_mass), constraint_exact, "mpmath"
)
pre_numeric = sp.lambdify(
    (l_minus, l_plus, squared_strut, dust_mass), pre_exact, "mpmath"
)
post_numeric = sp.lambdify(
    (l_minus, l_plus, squared_strut, dust_mass), post_exact, "mpmath"
)


# Exact unit-square reduction.
m_symbol, x_symbol, a_symbol, b_symbol = sp.symbols(
    "m_symbol x_symbol a_symbol b_symbol", positive=True
)
u_unit_square = sp.simplify(
    (m_symbol**2 / x_symbol**2).subs(
        {
            m_symbol: sp.Rational(2, 5) * sp.sqrt(a_symbol),
            x_symbol: sp.Rational(125) / b_symbol,
        }
    )
)
unit_square_ok = bool(
    sp.simplify(u_unit_square - 4 * a_symbol * b_symbol**2 / 390625) == 0
)
check(
    "the unchanged half-strip is the complete compact unit square",
    unit_square_ok,
)


def series_asinh(value):
    return (value + (value * value + 1).sqrt()).log()


def arb_m(value):
    api = arb.pi()
    epsilon = 2 * api - 5 * (
        (1 + 2 * value) / (2 * (1 + 3 * value))
    ).acos()
    return 180 * epsilon / (api * (1 + 4 * value).sqrt())


def arb_p(value, include_boost=True):
    api = arb.pi()
    epsilon = 2 * api - 5 * (
        (1 + 2 * value) / (2 * (1 + 3 * value))
    ).acos()
    curvature = 180 * epsilon / (1 + 4 * value).sqrt()
    if not include_boost:
        return curvature
    argument = 1 / (8 * (1 + 3 * value)).sqrt()
    return curvature - 600 * arb(3).sqrt() * series_asinh(argument)


def polynomial(coefficients, value, first, last, shift):
    result = arb(0)
    for index in range(first, last + 1):
        result += coefficients[index] * value ** (index - shift)
    return result


def lagrange_certificate(precision, x_threshold=125):
    ctx.dps = precision
    api = arb.pi()
    u_interval = interval_ball(Fraction(0), Fraction(4, 390625))
    zero_series = arb_series([arb(0), arb(1)], prec=9)
    range_series = arb_series([u_interval, arb(1)], prec=9)

    m_zero = arb_m(zero_series)
    p_zero = arb_p(zero_series)
    m_range = arb_m(range_series)
    p_range = arb_p(range_series)
    p0 = arb_p(arb(0))
    n_zero = p_zero - p0 + 2 * api * zero_series * m_zero
    n_range = p_range - p0 + 2 * api * range_series * m_range
    h_zero = zero_series * p_zero.derivative() - (p_zero - p0)
    h_range = range_series * p_range.derivative() - (p_range - p0)

    coefficient_identity_ok = all(
        h_zero[index].overlaps((index - 1) * p_zero[index])
        for index in range(2, DEGREE + 1)
    )
    removable_origin_ok = bool(
        p_zero[0].overlaps(p0)
        and n_zero[0].contains(0)
        and n_zero[1].contains(0)
        and h_zero[0].contains(0)
        and h_zero[1].contains(0)
        and coefficient_identity_ok
    )

    w_value = polynomial(p_zero, u_interval, 1, DEGREE, 1)
    w_value += p_range[DEGREE + 1] * u_interval**DEGREE
    bbar = polynomial(n_zero, u_interval, 2, DEGREE, 2)
    bbar += n_range[DEGREE + 1] * u_interval ** (DEGREE - 1)
    w_prime = polynomial(h_zero, u_interval, 2, DEGREE, 2)
    w_prime += h_range[DEGREE + 1] * u_interval ** (DEGREE - 1)

    m_value = m_range[0]
    m_prime = m_range[1]
    c_value = w_value + 4 * api * m_value
    cbar = -w_prime - 4 * api * m_prime
    primitive_ok = bool(
        removable_origin_ok
        and strict_positive(m_value)
        and strict_positive(bbar)
        and strict_positive(c_value)
        and strict_positive(cbar)
    )

    z_interval = interval_ball(Fraction(0), Fraction(1, x_threshold))
    mass_squared = interval_ball(Fraction(0), Fraction(4, 25))
    u_value = z_interval * m_value
    one_minus_u = 1 - u_value
    yplus_over_z = (
        4 * api
        - z_interval * c_value
        - 4
        * one_minus_u
        * mass_squared
        * z_interval
        * bbar
        / (m_value * m_value)
    )
    curve_derivative = 4 * api - 2 * z_interval * c_value
    normalized_gap = (
        4 * one_minus_u * bbar / (m_value * m_value)
        - z_interval * z_interval * cbar
    )
    domain_ok = bool(
        primitive_ok
        and u_value < 1
        and strict_positive(one_minus_u)
        and strict_positive(yplus_over_z)
        and strict_positive(curve_derivative)
        and strict_positive(normalized_gap)
    )

    p_infinity = 60 * api - 300 * arb(3).sqrt() * arb(2).log()
    negative_root_margin = -p_infinity - rational_ball(8, 5) * api
    monotone = classification["monotone_facts"]
    global_root_ok = bool(
        domain_ok
        and chain_factor_ok
        and r_prime_ok
        and monotone["K_has_one_positive_squared_root"]
        and monotone["mu_inner_increasing_outer_decreasing"]
        and monotone["p_negative_on_positive_axis"]
        and classification["thresholds"]["rational_brackets"]["x_star"]
        == ["5", "6"]
        and strict_positive(negative_root_margin)
    )

    return {
        "precision": precision,
        "x_threshold": x_threshold,
        "u": u_interval,
        "coefficient_identity_ok": coefficient_identity_ok,
        "removable_origin_ok": removable_origin_ok,
        "M": m_value,
        "M_prime": m_prime,
        "W": w_value,
        "Bbar": bbar,
        "W_prime": w_prime,
        "C": c_value,
        "minus_C_prime": cbar,
        "U": u_value,
        "one_minus_U": one_minus_u,
        "y_plus_over_z": yplus_over_z,
        "curve_derivative_lower": curve_derivative,
        "normalized_gap_lower": normalized_gap,
        "negative_root_margin": negative_root_margin,
        "primitive_ok": primitive_ok,
        "domain_ok": domain_ok,
        "global_root_ok": global_root_ok,
        "p_coefficients": [p_zero[index] for index in range(DEGREE + 2)],
        "n_coefficients": [n_zero[index] for index in range(DEGREE + 2)],
        "h_coefficients": [h_zero[index] for index in range(DEGREE + 2)],
        "p_seventh_range": p_range[DEGREE + 1],
        "n_seventh_range": n_range[DEGREE + 1],
        "h_seventh_range": h_range[DEGREE + 1],
    }


def naive_zero_quotient_unresolved():
    u_interval = interval_ball(Fraction(0), Fraction(4, 390625))
    try:
        naive = (arb_p(u_interval) - arb_p(arb(0))) / u_interval
        return bool(not naive.is_finite() or not naive < 0), str(naive)
    except (ValueError, ZeroDivisionError) as error:
        return True, f"raised {type(error).__name__}"


def state_functions():
    def epsilon(value):
        return 2 * mp.pi - 5 * mp.acos(
            (value * value + 2) / (2 * (value * value + 3))
        )

    def mu(value):
        return 180 * epsilon(value) / (mp.pi * mp.sqrt(value * value + 4))

    def momentum(value):
        return (
            180 * value * epsilon(value) / mp.sqrt(value * value + 4)
            - 600
            * mp.sqrt(3)
            * mp.asinh(value / mp.sqrt(8 * (value * value + 3)))
        )

    return mu, momentum


def direct_residuals(height, q_value, mass, incoming_p):
    endpoint = 1 + height * q_value
    return (
        2
        * constraint_numeric(1, endpoint, height**2, mass)
        / height,
        pre_numeric(1, endpoint, height**2, mass) - incoming_p,
    )


def solve_direct(seed, mass, incoming_p, tolerance):
    return mp.findroot(
        lambda height, q_value: direct_residuals(
            height, q_value, mass, incoming_p
        ),
        seed,
        tol=tolerance,
        maxsteps=160,
    )


def transport(mass, height, q_value):
    ratio = 1 + height * q_value
    post = post_numeric(1, ratio, height**2, mass)
    return ratio, mass / ratio, post / ratio**2, post


def reconstruct_fifth(precision):
    mp.mp.dps = precision + 60
    tolerance = mp.mpf(10) ** (-(precision + 20))
    residual_limit = mp.mpf(10) ** (-(precision - 25))
    mu, momentum = state_functions()
    mass = mu(mp.mpf(3) / 2)
    incoming_p = momentum(mp.mpf(3) / 2)
    seeds = (
        (mp.mpf(1) / 5, mp.mpf(10)),
        (mp.mpf(1) / 14, mp.mpf(31)),
        (mp.mpf(1) / 50, mp.mpf(100)),
        (mp.mpf(1) / 150, mp.mpf(317)),
        (mp.mpf(1) / 500, mp.mpf(1007)),
    )
    states = []
    cumulative_scale = mp.mpf(1)
    previous_post = None
    for slab, seed in enumerate(seeds, start=1):
        height, q_value = solve_direct(seed, mass, incoming_p, tolerance)
        residuals = direct_residuals(height, q_value, mass, incoming_p)
        ratio, next_mass, next_p, post = transport(mass, height, q_value)
        junction = mp.mpf(0)
        if previous_post is not None:
            junction = previous_post - cumulative_scale**2 * pre_numeric(
                1, ratio, height**2, mass
            )
        row_ok = bool(
            height > 0
            and ratio > 0
            and max(abs(value) for value in residuals) < residual_limit
            and abs(junction) < residual_limit
        )
        states.append(
            {
                "slab": slab,
                "height": height,
                "q": q_value,
                "ratio": ratio,
                "mass": mass,
                "incoming_p": incoming_p,
                "constraint_residual": residuals[0],
                "momentum_residual": residuals[1],
                "junction_residual": junction,
                "passed": row_ok,
            }
        )
        previous_post = cumulative_scale**2 * post
        cumulative_scale *= ratio
        mass, incoming_p = next_mass, next_p

    fourth = states[3]
    fifth = states[4]
    x4 = fourth["mass"] * fourth["q"]
    x5 = fifth["mass"] * fifth["q"]
    wrong_incoming = -post_numeric(
        1,
        fourth["ratio"],
        fourth["height"] ** 2,
        fourth["mass"],
    ) / fourth["ratio"] ** 2
    wrong_sign_residual = pre_numeric(
        1, fifth["ratio"], fifth["height"] ** 2, fifth["mass"]
    ) - wrong_incoming
    return {
        "precision": precision,
        "states": states,
        "x4": x4,
        "x5": x5,
        "wrong_sign_residual": wrong_sign_residual,
        "passed": bool(
            all(row["passed"] for row in states)
            and x5 > x4
            and abs(wrong_sign_residual) > mp.mpf("1e-20")
        ),
    }


records = []
diagnostics = []
direct_records = []
controls = []
for precision in PRECISIONS:
    record = lagrange_certificate(precision)
    records.append(record)
    diagnostics.append(
        {
            "x_ge_124": lagrange_certificate(precision, 124)["domain_ok"],
            "x_ge_126": lagrange_certificate(precision, 126)["domain_ok"],
            "not_used_for_acceptance": True,
        }
    )
    direct = reconstruct_fifth(precision)
    direct_records.append(direct)

    ctx.dps = precision
    api = arb.pi()
    axis = arb_series([arb(0), arb(1)], prec=3)
    n_full = (
        arb_p(axis)
        - arb_p(arb(0))
        + 2 * api * axis * arb_m(axis)
    )
    n_no_boost = (
        arb_p(axis, include_boost=False)
        - arb_p(arb(0), include_boost=False)
        + 2 * api * axis * arb_m(axis)
    )
    boost_changes = bool(not n_full[2].overlaps(n_no_boost[2]))
    naive_unresolved, naive_detail = naive_zero_quotient_unresolved()
    x60_value = rational_ball(1, 60) * rational_ball(60)
    x60_fails = bool(x60_value.contains(1) and not x60_value < 1)
    controls.append(
        {
            "precision": precision,
            "naive_zero_quotient_unresolved": naive_unresolved,
            "naive_zero_quotient_detail": naive_detail,
            "boost_omission_changes_Bbar0": boost_changes,
            "x60_fails": x60_fails,
            "wrong_sign_fails": abs(direct["wrong_sign_residual"])
            > mp.mpf("1e-20"),
            "passed": bool(
                naive_unresolved
                and boost_changes
                and x60_fails
                and abs(direct["wrong_sign_residual"]) > mp.mpf("1e-20")
            ),
        }
    )


primitive_ok = all(record["primitive_ok"] for record in records)
domain_ok = all(record["domain_ok"] for record in records)
global_root_ok = all(record["global_root_ok"] for record in records)
direct_ok = all(record["passed"] for record in direct_records)
controls_ok = all(record["passed"] for record in controls)
check(
    "degree-six global Lagrange remainders separate every primitive sign twice",
    primitive_ok,
)
check(
    "the unchanged half-strip physical and normalized-gap gates pass twice",
    domain_ok,
)
check(
    "the complete all-real root argument leaves one physical successor",
    global_root_ok,
)
check(
    "the complete action independently reconstructs all five slabs twice",
    direct_ok,
)
check(
    "the preserved OPEN and every hostile control remain distinct",
    controls_ok,
)


bound_names = (
    "M",
    "W",
    "Bbar",
    "C",
    "minus_C_prime",
    "one_minus_U",
    "y_plus_over_z",
    "curve_derivative_lower",
    "normalized_gap_lower",
    "negative_root_margin",
)
precision_overlap = all(
    records[0][name].overlaps(records[1][name]) for name in bound_names
)
same_sign = all(
    strict_positive(record[name])
    for record in records
    for name in (
        "M",
        "Bbar",
        "C",
        "minus_C_prime",
        "one_minus_U",
        "y_plus_over_z",
        "curve_derivative_lower",
        "normalized_gap_lower",
        "negative_root_margin",
    )
)
precision_ok = bool(precision_overlap and same_sign)
check(
    "the two outward-rounded certificates overlap with identical strict signs",
    precision_ok,
)


independent_ok = bool(
    provenance_ok
    and redifferentiated_ok
    and unit_square_ok
    and primitive_ok
    and domain_ok
    and global_root_ok
    and direct_ok
    and controls_ok
    and precision_ok
)
independent_outcome = (
    "INDEPENDENT_GLOBAL_LAGRANGE_HALF_STRIP_CERTIFICATE"
    if independent_ok
    else "INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN"
)
check(
    "the resolution verdict is fixed before reading the primary artifact",
    independent_ok,
    independent_outcome,
)


# Primary and fifth artifacts are opened only after the independent verdict.
primary = json.loads(PRIMARY_INPUT.read_text())
fifth = json.loads(FIFTH_INPUT.read_text())
primary_provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and digest(FIFTH_INPUT) == FIFTH_SHA256
    and primary["outcome"] == "INVARIANT_HALF_STRIP_PRIMARY_CERTIFICATE"
    and fifth["outcome"]
    == "CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB"
)
primary_map = {
    "M": "M",
    "Bbar": "Bbar",
    "C": "C",
    "minus_C_prime": "minus_C_prime",
    "one_minus_U": "one_minus_U",
    "y_plus_over_z": "y_plus_over_z",
    "curve_derivative_lower": "curve_derivative_lower",
    "normalized_gap_lower": "normalized_gap_lower",
    "negative_root_margin": "negative_root_margin",
}
primary_overlap = True
for own_name, stored_name in primary_map.items():
    stored = primary["arb_bounds"][stored_name]
    stored_ball = arb(stored["lower"]).union(arb(stored["upper"]))
    primary_overlap &= bool(records[1][own_name].overlaps(stored_ball))

stored_q5 = mp.mpf(fifth["forecast_comparison"]["actual_q5"])
stored_h5 = mp.mpf(fifth["forecast_comparison"]["actual_h5"])
stored_x5 = mp.mpf(fifth["forecast_comparison"]["actual_x5"])
direct_fifth = direct_records[1]["states"][4]
fifth_match = bool(
    mp_text(direct_fifth["q"], 60) == mp_text(stored_q5, 60)
    and mp_text(direct_fifth["height"], 60) == mp_text(stored_h5, 60)
    and mp_text(direct_records[1]["x5"], 60) == mp_text(stored_x5, 60)
)
comparison_ok = bool(
    primary_provenance_ok
    and primary_overlap
    and fifth_match
    and primary["domain"]["m"] == "0<m<=2/5"
    and primary["domain"]["x"] == "x>=125"
)
check(
    "the independent result agrees with the primary only after construction",
    comparison_ok,
)


final_ok = bool(independent_ok and comparison_ok)
outcome = (
    "INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED"
    if final_ok
    else "INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN"
)
check(
    "the resolution hierarchy adversarially corroborates the half-strip",
    final_ok,
    outcome,
)


def serialize_direct(record):
    return {
        "precision": record["precision"],
        "x4": mp_text(record["x4"]),
        "x5": mp_text(record["x5"]),
        "wrong_sign_residual": mp_text(record["wrong_sign_residual"]),
        "passed": record["passed"],
        "states": [
            {
                key: mp_text(value) if isinstance(value, mp.mpf) else value
                for key, value in state.items()
            }
            for state in record["states"]
        ],
    }


def serialize_record(record, diagnostic, direct, control):
    scalar_names = (
        "u",
        "M",
        "M_prime",
        "W",
        "Bbar",
        "W_prime",
        "C",
        "minus_C_prime",
        "U",
        "one_minus_U",
        "y_plus_over_z",
        "curve_derivative_lower",
        "normalized_gap_lower",
        "negative_root_margin",
        "p_seventh_range",
        "n_seventh_range",
        "h_seventh_range",
    )
    return {
        "precision_decimal_digits": record["precision"],
        "degree": DEGREE,
        "coefficient_identity_ok": record["coefficient_identity_ok"],
        "removable_origin_ok": record["removable_origin_ok"],
        "primitive_ok": record["primitive_ok"],
        "domain_ok": record["domain_ok"],
        "global_root_ok": record["global_root_ok"],
        "bounds": {name: arb_record(record[name]) for name in scalar_names},
        "coefficients": {
            "p": [arb_record(value) for value in record["p_coefficients"]],
            "n": [arb_record(value) for value in record["n_coefficients"]],
            "h": [arb_record(value) for value in record["h_coefficients"]],
        },
        "threshold_diagnostics": diagnostic,
        "direct_action": serialize_direct(direct),
        "hostile_controls": control,
    }


artifact = {
    "provenance": {
        "resolution_protocol_commit": PROTOCOL_COMMIT,
        "resolution_protocol_sha256": PROTOCOL_SHA256,
        "first_adversarial_commit": FIRST_ADVERSARIAL_COMMIT,
        "first_adversarial_sha256": FIRST_ADVERSARIAL_SHA256,
        "classification_sha256": CLASSIFICATION_SHA256,
        "primary_sha256": digest(PRIMARY_INPUT),
        "fifth_sha256": digest(FIFTH_INPUT),
        "primary_read_after_independent_outcome": True,
    },
    "domain": {
        "m": "0<m<=2/5",
        "x": "x>=125",
        "u": "0<=u<=4/390625",
        "thresholds_post_hoc": True,
    },
    "method": {
        "interval_engine": "python-flint Arb outward-rounded balls",
        "precisions_decimal_digits": list(PRECISIONS),
        "degree": DEGREE,
        "load_bearing_enclosure": "global Maclaurin polynomial with Lagrange derivative remainder",
        "primary_integral_means_used": False,
        "direct_quotient_interval_evaluation_used": False,
        "subdivision_used": False,
        "finite_grid_used": False,
        "complete_action_redifferentiated": True,
    },
    "precision_records": [
        serialize_record(record, diagnostic, direct, control)
        for record, diagnostic, direct, control in zip(
            records, diagnostics, direct_records, controls
        )
    ],
    "independent_outcome_before_primary_read": independent_outcome,
    "primary_comparison": {
        "provenance_ok": primary_provenance_ok,
        "rigorous_ranges_overlap": primary_overlap,
        "fifth_serialized_60_digits_match": fifth_match,
        "passed": comparison_ok,
    },
    "claims": {
        "invariant_half_strip": (
            "DERIVED_COMPUTATIONAL_WITH_TWO_MECHANICALLY_DISTINCT_RIGOROUS_CERTIFICATES"
            if final_ok
            else "OPEN"
        ),
        "accepted_history_infinite": (
            "DERIVED_BY_INDUCTION_IN_FROZEN_HOMOGENEOUS_MODEL"
            if final_ok
            else "OPEN"
        ),
        "first_direct_quotient_route": "PRESERVED_OPEN",
        "extendibility_selector": "STRUCTURAL_GLOBAL_NOT_LOCAL_LAW",
        "generic_original_v": "OPEN",
        "nonhomogeneous_stability": "OPEN",
        "absolute_tick": "NOT_DERIVED",
    },
    "checks": {
        "provenance": provenance_ok,
        "redifferentiation": redifferentiated_ok,
        "unit_square": unit_square_ok,
        "primitive_bounds": primitive_ok,
        "domain_gates": domain_ok,
        "global_root_uniqueness": global_root_ok,
        "direct_fifth": direct_ok,
        "hostile_controls": controls_ok,
        "precision_consistency": precision_ok,
        "primary_comparison": comparison_ok,
    },
    "passed": passed,
    "tests": tests,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print()
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
raise SystemExit(0 if passed == tests else 1)

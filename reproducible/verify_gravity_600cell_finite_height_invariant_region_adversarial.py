#!/usr/bin/env python3
"""Direct-quotient adversarial certificate for the invariant half-strip."""

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
    / "gravity_600cell_finite_height_invariant_region_adversarial_protocol.md"
)
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
FIFTH_INPUT = HERE / "gravity_600cell_finite_height_asymptotic_map.json"
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_invariant_region.json"
OUTPUT = HERE / "gravity_600cell_finite_height_invariant_region_adversarial.json"

PROTOCOL_COMMIT = "26ef9c3"
PROTOCOL_SHA256 = (
    "8075d5ea7ff2be80eaab5e6c7243d1c5c282114b097bebbba3eec500bf186b5b"
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
LEAF_COUNT = 64
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
    try:
        return bool(
            value.is_finite()
            and value > 0
            and value.lower() > 0
            and not value.contains(0)
        )
    except (ValueError, ZeroDivisionError):
        return False


def strict_negative(value):
    try:
        return bool(
            value.is_finite()
            and value < 0
            and value.upper() < 0
            and not value.contains(0)
        )
    except (ValueError, ZeroDivisionError):
        return False


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
        "finite": bool(value.is_finite()),
    }


def hull(values):
    result = values[0]
    for value in values[1:]:
        result = result.union(value)
    return result


def mp_text(value, digits=75):
    return mp.nstr(value, digits)


# Only accepted pre-invariant inputs are read here.  The primary invariant
# artifact is deliberately not read until the independent verdict exists.
classification = json.loads(CLASSIFICATION_INPUT.read_text())
provenance_ok = bool(
    digest(PROTOCOL) == PROTOCOL_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and classification["outcome"]
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
)
check(
    "the frozen adversarial protocol and root-classification input are exact",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}; primary invariant artifact remains unread",
)


# Re-differentiate the complete one-slab action with symbols and names that
# are independent of the primary invariant verifier.
l0, l1, rho, dust = sp.symbols("l0 l1 rho dust", positive=True)
delta = l1 - l0
proper_height = sp.sqrt(rho + delta**2 / 4)
dihedral_cosine = (delta**2 + 2 * rho) / (2 * (delta**2 + 3 * rho))
boost_argument = delta / sp.sqrt(8 * (delta**2 + 3 * rho))
action = (
    360
    * (l0 + l1)
    * proper_height
    * (2 * sp.pi - 5 * sp.acos(dihedral_cosine))
    + 600
    * sp.sqrt(3)
    * (l0**2 - l1**2)
    * sp.asinh(boost_argument)
    - 8 * sp.pi * dust * sp.sqrt(rho)
)
constraint_exact = rho * sp.diff(action, rho)
pre_exact = -l0 * sp.diff(action, l0) / 2
post_exact = l1 * sp.diff(action, l1) / 2

scale = sp.symbols("scale", positive=True)
scaling = {l0: scale * l0, l1: scale * l1, rho: scale**2 * rho, dust: scale * dust}
action_scaling_ok = bool(
    sp.simplify(sp.powsimp(action.subs(scaling) - scale**2 * action, force=True))
    == 0
    and sp.simplify(
        sp.powsimp(post_exact.subs(scaling) - scale**2 * post_exact, force=True)
    )
    == 0
)

q = sp.symbols("q", real=True, nonzero=True)
q2 = q**2
eps_q = 2 * sp.pi - 5 * sp.acos((q2 + 2) / (2 * (q2 + 3)))
mu_q = 180 * eps_q / (sp.pi * sp.sqrt(q2 + 4))
p_q = (
    180 * q * eps_q / sp.sqrt(q2 + 4)
    - 600 * sp.sqrt(3) * sp.asinh(q / sp.sqrt(8 * (q2 + 3)))
)
mu_prime_direct = sp.diff(mu_q, q)
p_prime_direct = sp.diff(p_q, q)
common_derivative_ok = bool(
    sp.simplify(
        sp.trigsimp(p_prime_direct + 4 * sp.pi * mu_prime_direct / q)
    )
    == 0
)

mass_symbol, mu_symbol, pi_symbol = sp.symbols("mass_symbol mu_symbol pi_symbol")
r_function = p_q - pi_symbol + 4 * sp.pi * (mu_q - mass_symbol) / q
r_prime_reduced = sp.simplify(
    sp.diff(r_function, q).subs(
        sp.diff(mu_q, q), -q * sp.diff(p_q, q) / (4 * sp.pi)
    )
)
r_prime_ok = bool(
    sp.simplify(r_prime_reduced - 4 * sp.pi * (mass_symbol - mu_q) / q**2)
    == 0
)
redifferentiated_ok = bool(action_scaling_ok and common_derivative_ok and r_prime_ok)
check(
    "the complete action and dual root derivative are independently redifferentiated",
    redifferentiated_ok,
)

constraint_numeric = sp.lambdify((l0, l1, rho, dust), constraint_exact, "mpmath")
pre_numeric = sp.lambdify((l0, l1, rho, dust), pre_exact, "mpmath")
post_numeric = sp.lambdify((l0, l1, rho, dust), post_exact, "mpmath")


# The unit-square coordinates are used only to derive the complete rational
# u range.  The 64-leaf certificate itself is one-dimensional in u.
m_symbol, x_symbol, a_symbol, b_symbol = sp.symbols(
    "m_symbol x_symbol a_symbol b_symbol", positive=True
)
u_from_mx = m_symbol**2 / x_symbol**2
unit_square_substitution = {
    m_symbol: sp.Rational(2, 5) * sp.sqrt(a_symbol),
    x_symbol: sp.Rational(125, 1) / b_symbol,
}
u_unit_square = sp.simplify(u_from_mx.subs(unit_square_substitution))
unit_square_ok = bool(
    sp.simplify(u_unit_square - 4 * a_symbol * b_symbol**2 / 390625) == 0
)
check(
    "the half-strip compactifies to the frozen unit square and u interval",
    unit_square_ok,
    "a=(5m/2)^2, b=125/x, 0<=u<=4/390625",
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


def make_leaf(index, precision):
    ctx.dps = precision
    u_max = Fraction(4, 390625)
    left = index * u_max / LEAF_COUNT
    right = (index + 1) * u_max / LEAF_COUNT
    u_ball = interval_ball(left, right)
    api = arb.pi()
    p0 = arb_p(arb(0))
    jet = arb_series([u_ball, arb(1)], prec=5)
    m_jet = arb_m(jet)
    p_jet = arb_p(jet)
    n_jet = p_jet - p0 + 2 * api * jet * m_jet

    m_value = m_jet[0]
    m_prime = m_jet[1]
    construction = "direct_quotient"
    if index:
        p_value = p_jet[0]
        p_prime = p_jet[1]
        denominator = u_ball * u_ball
        w_value = (p_value - p0) / u_ball
        bbar = (p_value - p0 + 2 * api * u_ball * m_value) / denominator
        w_prime = (u_ball * p_prime - (p_value - p0)) / denominator
    else:
        construction = "axis_lagrange_taylor"
        zero = arb_series([arb(0), arb(1)], prec=5)
        m_zero = arb_m(zero)
        p_zero = arb_p(zero)
        n_zero = p_zero - p0 + 2 * api * zero * m_zero

        p1_zero = p_zero[1]
        p2_zero = 2 * p_zero[2]
        p3_zero = 6 * p_zero[3]
        n2_zero = 2 * n_zero[2]
        n3_zero = 6 * n_zero[3]
        p3_range = 6 * p_jet[3]
        p4_range = 24 * p_jet[4]
        n4_range = 24 * n_jet[4]
        u_squared = u_ball * u_ball

        w_value = (
            p1_zero + p2_zero * u_ball / 2 + p3_range * u_squared / 6
        )
        bbar = (
            n2_zero / 2
            + n3_zero * u_ball / 6
            + n4_range * u_squared / 24
        )
        w_prime_center = p2_zero / 2 + p3_zero * u_ball / 3
        remainder = (
            arb(0, 1)
            * rational_ball(5, 24)
            * p4_range.abs_upper()
            * u_squared
        )
        w_prime = w_prime_center + remainder

    c_value = w_value + 4 * api * m_value
    cbar = -w_prime - 4 * api * m_prime
    return {
        "index": index,
        "left": left,
        "right": right,
        "u": u_ball,
        "construction": construction,
        "M": m_value,
        "W": w_value,
        "Bbar": bbar,
        "C": c_value,
        "minus_C_prime": cbar,
    }


def evaluate_domain(leaves, x_threshold):
    api = arb.pi()
    z_interval = interval_ball(Fraction(0), Fraction(1, x_threshold))
    mass_squared = interval_ball(Fraction(0), Fraction(4, 25))
    gate_rows = []
    for leaf in leaves:
        m_value = leaf["M"]
        bbar = leaf["Bbar"]
        c_value = leaf["C"]
        cbar = leaf["minus_C_prime"]
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
        gap = (
            4 * one_minus_u * bbar / (m_value * m_value)
            - z_interval * z_interval * cbar
        )
        row_ok = bool(
            strict_positive(m_value)
            and strict_positive(bbar)
            and strict_positive(c_value)
            and strict_positive(cbar)
            and u_value < 1
            and strict_positive(one_minus_u)
            and strict_positive(yplus_over_z)
            and strict_positive(curve_derivative)
            and strict_positive(gap)
        )
        gate_rows.append(
            {
                "index": leaf["index"],
                "U": u_value,
                "one_minus_U": one_minus_u,
                "y_plus_over_z": yplus_over_z,
                "curve_derivative_lower": curve_derivative,
                "normalized_gap_lower": gap,
                "passed": row_ok,
            }
        )
    return gate_rows


def direct_quotient_zero_control():
    interval = interval_ball(Fraction(0), Fraction(4, 390625 * LEAF_COUNT))
    try:
        naive = (arb_p(interval) - arb_p(arb(0))) / interval
        unresolved = bool(
            not naive.is_finite() or naive.contains(0) or not strict_negative(naive)
        )
        return unresolved, str(naive)
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
        2 * constraint_numeric(1, endpoint, height**2, mass) / height,
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
    mass0 = mu(mp.mpf(3) / 2)
    pi0 = momentum(mp.mpf(3) / 2)

    seeds = (
        (mp.mpf(1) / 5, mp.mpf(10)),
        (mp.mpf(1) / 14, mp.mpf(31)),
        (mp.mpf(1) / 50, mp.mpf(100)),
        (mp.mpf(1) / 150, mp.mpf(317)),
        (mp.mpf(1) / 500, mp.mpf(1007)),
    )
    states = []
    mass = mass0
    incoming_p = pi0
    cumulative_scale = mp.mpf(1)
    previous_post = None
    for index, seed in enumerate(seeds, start=1):
        height, q_value = solve_direct(seed, mass, incoming_p, tolerance)
        residuals = direct_residuals(height, q_value, mass, incoming_p)
        ratio, next_mass, next_p, post = transport(mass, height, q_value)
        junction = mp.mpf(0)
        if previous_post is not None:
            junction = previous_post - cumulative_scale**2 * pre_numeric(
                1, ratio, height**2, mass
            )
        physical = bool(height > 0 and ratio > 0)
        row_ok = bool(
            physical
            and max(abs(value) for value in residuals) < residual_limit
            and abs(junction) < residual_limit
        )
        states.append(
            {
                "slab": index,
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

    # Reverse the outgoing fourth momentum before the fifth solve, but keep
    # the accepted fifth geometry fixed.  It must fail the pre-equation.
    wrong_incoming = -states[3]["ratio"] ** (-2) * post_numeric(
        1,
        states[3]["ratio"],
        states[3]["height"] ** 2,
        states[3]["mass"],
    )
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


def serialize_leaf(leaf):
    return {
        "index": leaf["index"],
        "left": str(leaf["left"]),
        "right": str(leaf["right"]),
        "construction": leaf["construction"],
        **{
            name: arb_record(leaf[name])
            for name in ("M", "W", "Bbar", "C", "minus_C_prime")
        },
    }


def serialize_gate(row):
    return {
        "index": row["index"],
        "passed": row["passed"],
        **{
            name: arb_record(row[name])
            for name in (
                "U",
                "one_minus_U",
                "y_plus_over_z",
                "curve_derivative_lower",
                "normalized_gap_lower",
            )
        },
    }


def serialize_fifth(record):
    return {
        "precision": record["precision"],
        "x4": mp_text(record["x4"]),
        "x5": mp_text(record["x5"]),
        "wrong_sign_residual": mp_text(record["wrong_sign_residual"]),
        "passed": record["passed"],
        "states": [
            {
                key: (
                    mp_text(value)
                    if isinstance(value, mp.mpf)
                    else value
                )
                for key, value in state.items()
            }
            for state in record["states"]
        ],
    }


precision_records = []
all_primitive_ok = True
all_domain_ok = True
all_global_root_ok = True
all_direct_ok = True
all_controls_ok = True

monotone = classification["monotone_facts"]
classification_logic_ok = bool(
    monotone["K_has_one_positive_squared_root"]
    and monotone["mu_inner_increasing_outer_decreasing"]
    and monotone["p_negative_on_positive_axis"]
    and classification["thresholds"]["rational_brackets"]["x_star"] == ["5", "6"]
)

for precision in PRECISIONS:
    ctx.dps = precision
    leaves = [make_leaf(index, precision) for index in range(LEAF_COUNT)]
    domain = evaluate_domain(leaves, 125)
    diagnostic_124 = evaluate_domain(leaves, 124)
    diagnostic_126 = evaluate_domain(leaves, 126)
    primitive_ok = all(
        strict_positive(leaf[name])
        for leaf in leaves
        for name in ("M", "Bbar", "C", "minus_C_prime")
    )
    domain_ok = bool(primitive_ok and all(row["passed"] for row in domain))

    api = arb.pi()
    p_infinity = 60 * api - 300 * arb(3).sqrt() * arb(2).log()
    negative_root_margin = -p_infinity - rational_ball(8, 5) * api
    global_root_ok = bool(
        classification_logic_ok
        and common_derivative_ok
        and r_prime_ok
        and strict_positive(negative_root_margin)
        and domain_ok
    )

    naive_unresolved, naive_detail = direct_quotient_zero_control()
    fifth_record = reconstruct_fifth(precision)

    # Symbolic boost-omission control at u=0, evaluated independently of the
    # interval quotient rows.
    axis = arb_series([arb(0), arb(1)], prec=3)
    p0_full = arb_p(arb(0))
    p0_no_boost = arb_p(arb(0), include_boost=False)
    n_full = arb_p(axis) - p0_full + 2 * api * axis * arb_m(axis)
    n_no_boost = (
        arb_p(axis, include_boost=False)
        - p0_no_boost
        + 2 * api * axis * arb_m(axis)
    )
    bbar0_full = n_full[2]
    bbar0_no_boost = n_no_boost[2]
    boost_omission_changes = bool(
        not bbar0_full.overlaps(bbar0_no_boost)
        and abs(bbar0_full - bbar0_no_boost).lower() > 0
    )
    boundary_x60 = rational_ball(1, 60) * rational_ball(60)
    x60_fails = bool(boundary_x60.contains(1) and not boundary_x60 < 1)
    controls_ok = bool(
        naive_unresolved
        and boost_omission_changes
        and x60_fails
        and abs(fifth_record["wrong_sign_residual"]) > mp.mpf("1e-20")
    )

    coefficient_hulls = {
        name: hull([leaf[name] for leaf in leaves])
        for name in ("M", "W", "Bbar", "C", "minus_C_prime")
    }
    gate_hulls = {
        name: hull([row[name] for row in domain])
        for name in (
            "U",
            "one_minus_U",
            "y_plus_over_z",
            "curve_derivative_lower",
            "normalized_gap_lower",
        )
    }
    gate_hulls["negative_root_margin"] = negative_root_margin

    precision_records.append(
        {
            "precision": precision,
            "leaves": leaves,
            "domain": domain,
            "diagnostic_124": diagnostic_124,
            "diagnostic_126": diagnostic_126,
            "coefficient_hulls": coefficient_hulls,
            "gate_hulls": gate_hulls,
            "primitive_ok": primitive_ok,
            "domain_ok": domain_ok,
            "global_root_ok": global_root_ok,
            "naive_zero_quotient_unresolved": naive_unresolved,
            "naive_zero_quotient_detail": naive_detail,
            "boost_omission_changes": boost_omission_changes,
            "x60_fails": x60_fails,
            "fifth": fifth_record,
            "controls_ok": controls_ok,
        }
    )
    all_primitive_ok &= primitive_ok
    all_domain_ok &= domain_ok
    all_global_root_ok &= global_root_ok
    all_direct_ok &= fifth_record["passed"]
    all_controls_ok &= controls_ok


check(
    "all 64 direct-quotient leaves separate every primitive sign at both precisions",
    all_primitive_ok,
    "the first leaf alone uses the frozen Taylor-Lagrange axis formula",
)
check(
    "the half-strip physical, momentum, monotonicity and gap gates pass twice",
    all_domain_ok,
    "precisions=160,256; no adaptive subdivision or threshold fallback",
)
check(
    "the independent all-real R-prime argument leaves one physical successor",
    all_global_root_ok,
)
check(
    "the complete redifferentiated action reconstructs the fifth successor twice",
    all_direct_ok,
)
check(
    "all frozen hostile controls distinguish the accepted construction",
    all_controls_ok,
)


low, high = precision_records
overlap_names = tuple(low["coefficient_hulls"]) + tuple(low["gate_hulls"])
precision_overlap = all(
    (
        low["coefficient_hulls"].get(name, low["gate_hulls"].get(name))
    ).overlaps(
        high["coefficient_hulls"].get(name, high["gate_hulls"].get(name))
    )
    for name in overlap_names
)
same_sign = all(
    strict_positive(record["coefficient_hulls"][name])
    for record in precision_records
    for name in ("M", "Bbar", "C", "minus_C_prime")
) and all(
    strict_positive(record["gate_hulls"][name])
    for record in precision_records
    for name in (
        "one_minus_U",
        "y_plus_over_z",
        "curve_derivative_lower",
        "normalized_gap_lower",
        "negative_root_margin",
    )
)
precision_consistency_ok = bool(precision_overlap and same_sign)
check(
    "the 160- and 256-digit rigorous ranges overlap with identical strict signs",
    precision_consistency_ok,
)


independent_certificate_ok = bool(
    provenance_ok
    and redifferentiated_ok
    and unit_square_ok
    and all_primitive_ok
    and all_domain_ok
    and all_global_root_ok
    and all_direct_ok
    and all_controls_ok
    and precision_consistency_ok
)
strict_refutation = any(
    strict_negative(leaf[name])
    for record in precision_records
    for leaf in record["leaves"]
    for name in ("M", "Bbar", "C", "minus_C_prime")
) or any(
    strict_negative(row[name])
    for record in precision_records
    for row in record["domain"]
    for name in (
        "one_minus_U",
        "y_plus_over_z",
        "curve_derivative_lower",
        "normalized_gap_lower",
    )
)
if independent_certificate_ok:
    independent_outcome = "INDEPENDENT_HALF_STRIP_CERTIFICATE"
elif strict_refutation or not all_direct_ok:
    independent_outcome = "INVARIANT_HALF_STRIP_ADVERSARIAL_REFUTED"
else:
    independent_outcome = "INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN"

check(
    "an independent verdict is fixed before the primary artifact is opened",
    independent_outcome == "INDEPENDENT_HALF_STRIP_CERTIFICATE",
    independent_outcome,
)


# Only now compare with the primary invariant artifact and accepted fifth
# serialization.  Neither file contributed a coefficient or sign above.
primary = json.loads(PRIMARY_INPUT.read_text())
fifth = json.loads(FIFTH_INPUT.read_text())
primary_provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and digest(FIFTH_INPUT) == FIFTH_SHA256
    and primary["outcome"] == "INVARIANT_HALF_STRIP_PRIMARY_CERTIFICATE"
    and fifth["outcome"]
    == "CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB"
)

primary_name_map = {
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
for own_name, primary_name in primary_name_map.items():
    source = (
        high["coefficient_hulls"]
        if own_name in high["coefficient_hulls"]
        else high["gate_hulls"]
    )
    stored = primary["arb_bounds"][primary_name]
    stored_ball = arb(stored["lower"]).union(arb(stored["upper"]))
    primary_overlap &= bool(source[own_name].overlaps(stored_ball))

stored_q5 = mp.mpf(fifth["forecast_comparison"]["actual_q5"])
stored_h5 = mp.mpf(fifth["forecast_comparison"]["actual_h5"])
stored_x5 = mp.mpf(fifth["forecast_comparison"]["actual_x5"])
direct_fifth = high["fifth"]["states"][4]
fifth_serialized_match = bool(
    mp_text(direct_fifth["q"], 60) == mp_text(stored_q5, 60)
    and mp_text(direct_fifth["height"], 60) == mp_text(stored_h5, 60)
    and mp_text(high["fifth"]["x5"], 60) == mp_text(stored_x5, 60)
)
comparison_ok = bool(
    primary_provenance_ok
    and primary_overlap
    and fifth_serialized_match
    and primary["domain"]["m"] == "0<m<=2/5"
    and primary["domain"]["x"] == "x>=125"
)
check(
    "the independent certificate agrees with primary bounds only after construction",
    comparison_ok,
)


final_ok = bool(independent_certificate_ok and comparison_ok)
if final_ok:
    outcome = "INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED"
elif independent_outcome == "INVARIANT_HALF_STRIP_ADVERSARIAL_REFUTED":
    outcome = independent_outcome
else:
    outcome = "INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN"
check(
    "the adversarial hierarchy corroborates the invariant half-strip",
    outcome == "INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED",
    outcome,
)


def serialize_precision(record):
    return {
        "precision_decimal_digits": record["precision"],
        "leaf_count": LEAF_COUNT,
        "primitive_ok": record["primitive_ok"],
        "domain_ok": record["domain_ok"],
        "global_root_ok": record["global_root_ok"],
        "controls_ok": record["controls_ok"],
        "naive_zero_quotient_unresolved": record[
            "naive_zero_quotient_unresolved"
        ],
        "naive_zero_quotient_detail": record["naive_zero_quotient_detail"],
        "boost_omission_changes": record["boost_omission_changes"],
        "x60_fails": record["x60_fails"],
        "coefficient_hulls": {
            name: arb_record(value)
            for name, value in record["coefficient_hulls"].items()
        },
        "gate_hulls": {
            name: arb_record(value)
            for name, value in record["gate_hulls"].items()
        },
        "leaves": [serialize_leaf(leaf) for leaf in record["leaves"]],
        "domain_gates": [serialize_gate(row) for row in record["domain"]],
        "threshold_diagnostics": {
            "x_ge_124_all_gates": all(
                row["passed"] for row in record["diagnostic_124"]
            ),
            "x_ge_126_all_gates": all(
                row["passed"] for row in record["diagnostic_126"]
            ),
            "not_used_for_acceptance": True,
        },
        "fifth_direct_action": serialize_fifth(record["fifth"]),
    }


artifact = {
    "provenance": {
        "adversarial_protocol_commit": PROTOCOL_COMMIT,
        "adversarial_protocol_sha256": PROTOCOL_SHA256,
        "classification_sha256": CLASSIFICATION_SHA256,
        "primary_sha256": digest(PRIMARY_INPUT),
        "fifth_sha256": digest(FIFTH_INPUT),
        "primary_read_after_independent_outcome": True,
    },
    "domain": {
        "m": "0<m<=2/5",
        "x": "x>=125",
        "unit_square": "a=(5m/2)^2, b=125/x",
        "u": "4*a*b^2/390625",
        "leaf_partition": "exactly 64 equal rational u intervals",
        "thresholds_post_hoc": True,
    },
    "method": {
        "interval_engine": "python-flint Arb outward-rounded balls",
        "precisions_decimal_digits": list(PRECISIONS),
        "away_from_axis": "direct quotient formulas on 63 leaves",
        "axis": "frozen Taylor-Lagrange formulas on first leaf only",
        "adaptive_subdivision": False,
        "primary_integral_means_used": False,
        "complete_action_redifferentiated": True,
    },
    "precision_records": [serialize_precision(record) for record in precision_records],
    "independent_outcome_before_primary_read": independent_outcome,
    "primary_comparison": {
        "provenance_ok": primary_provenance_ok,
        "rigorous_ranges_overlap": primary_overlap,
        "fifth_serialized_60_digits_match": fifth_serialized_match,
        "passed": comparison_ok,
    },
    "claims": {
        "invariant_half_strip": (
            "DERIVED_COMPUTATIONAL_WITH_TWO_RIGOROUS_INTERVAL_CERTIFICATES"
            if final_ok
            else "OPEN"
        ),
        "accepted_history_infinite": (
            "DERIVED_BY_INDUCTION_IN_FROZEN_HOMOGENEOUS_MODEL"
            if final_ok
            else "OPEN"
        ),
        "extendibility_selector": "STRUCTURAL_GLOBAL_NOT_LOCAL_LAW",
        "generic_original_v": "OPEN",
        "nonhomogeneous_stability": "OPEN",
        "absolute_tick": "NOT_DERIVED",
    },
    "checks": {
        "provenance": provenance_ok,
        "redifferentiation": redifferentiated_ok,
        "unit_square": unit_square_ok,
        "primitive_bounds": all_primitive_ok,
        "domain_gates": all_domain_ok,
        "global_root_uniqueness": all_global_root_ok,
        "direct_fifth": all_direct_ok,
        "hostile_controls": all_controls_ok,
        "precision_consistency": precision_consistency_ok,
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

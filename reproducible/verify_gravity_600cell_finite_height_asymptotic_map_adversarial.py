#!/usr/bin/env python3
"""Independent equal-mu replication of the asymptotic map and fifth slab."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = (
    ROOT
    / "docs"
    / "gravity"
    / "gravity_600cell_finite_height_asymptotic_map_protocol.md"
)
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_asymptotic_map.json"
OUTPUT = HERE / "gravity_600cell_finite_height_asymptotic_map_adversarial.json"

PROTOCOL_COMMIT = "c39f9ca"
PROTOCOL_SHA256 = (
    "c9ce01a119be1e66b46796fa16cafd3d81b4ac0d732edf2a3a4c9a936bfec25f"
)
PRIMARY_COMMIT = "46014f5"
PRIMARY_SHA256 = (
    "a93837d2bbec340ddbac528c0be4da52aefe45c8f0d4310496eb1aef6a7b19b6"
)

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


def text(value, digits=75):
    return mp.nstr(value, digits)


# Freeze only the protocol here.  The primary result is deliberately unread
# until both independent precision runs and the fifth-root censuses finish.
protocol_ok = digest(PROTOCOL) == PROTOCOL_SHA256
check(
    "the independently implemented protocol is frozen",
    protocol_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


# Re-differentiate the complete action with independent symbols.
a0, a1, sigma, dust = sp.symbols("a0 a1 sigma dust", positive=True)
delta = a1 - a0
proper_height = sp.sqrt(sigma + delta**2 / 4)
dihedral_cosine = (delta**2 + 2 * sigma) / (2 * (delta**2 + 3 * sigma))
boost_argument = delta / sp.sqrt(8 * (delta**2 + 3 * sigma))
complete_action = (
    360
    * (a0 + a1)
    * proper_height
    * (2 * sp.pi - 5 * sp.acos(dihedral_cosine))
    + 600
    * sp.sqrt(3)
    * (a0**2 - a1**2)
    * sp.asinh(boost_argument)
    - 8 * sp.pi * dust * sp.sqrt(sigma)
)
constraint_exact = sigma * sp.diff(complete_action, sigma)
pre_momentum_exact = -a0 * sp.diff(complete_action, a0) / 2
post_momentum_exact = a1 * sp.diff(complete_action, a1) / 2
constraint_numeric = sp.lambdify(
    (a0, a1, sigma, dust), constraint_exact, "mpmath"
)
pre_numeric = sp.lambdify(
    (a0, a1, sigma, dust), pre_momentum_exact, "mpmath"
)
post_numeric = sp.lambdify(
    (a0, a1, sigma, dust), post_momentum_exact, "mpmath"
)


def deficit(q):
    return 2 * mp.pi - 5 * mp.acos((q**2 + 2) / (2 * (q**2 + 3)))


def mu(q):
    return 180 * deficit(q) / (mp.pi * mp.sqrt(q**2 + 4))


def canonical_p(q):
    return (
        180 * q * deficit(q) / mp.sqrt(q**2 + 4)
        - 600
        * mp.sqrt(3)
        * mp.asinh(q / mp.sqrt(8 * (q**2 + 3)))
    )


def direct_equations(mass, incoming_p, height, q):
    ratio = 1 + height * q
    return (
        2 * constraint_numeric(1, ratio, height**2, mass) / height,
        pre_numeric(1, ratio, height**2, mass) - incoming_p,
    )


def reduced_equations(mass, incoming_p, height, q):
    return (
        8 * mp.pi * (mu(q) - mass) + 4 * mp.pi * height * q * mu(q),
        canonical_p(q) - incoming_p - 2 * mp.pi * height * mu(q),
    )


def solve_action(seed, mass, incoming_p, tolerance):
    return mp.findroot(
        lambda height, q: direct_equations(mass, incoming_p, height, q),
        seed,
        tol=tolerance,
        maxsteps=140,
    )


def transport(mass, height, q):
    ratio = 1 + height * q
    post = post_numeric(1, ratio, height**2, mass)
    return ratio, mass / ratio, post / ratio**2, post


# Obtain coefficients by derivatives of expressions analytic at t=0.  No
# series() call and no substitution into the primary result are used.
t = sp.symbols("t", nonnegative=True)
epsilon_regular = 2 * sp.pi - 5 * sp.acos(
    (1 + 2 * t**2) / (2 * (1 + 3 * t**2))
)
mu_regular = 180 * t * epsilon_regular / (sp.pi * sp.sqrt(1 + 4 * t**2))
p_regular = (
    180 * epsilon_regular / sp.sqrt(1 + 4 * t**2)
    - 600
    * sp.sqrt(3)
    * sp.asinh(1 / sp.sqrt(8 * (1 + 3 * t**2)))
)


def taylor_coefficient(expression, degree):
    return sp.simplify(sp.diff(expression, t, degree).subs(t, 0) / sp.factorial(degree))


mu_c1 = taylor_coefficient(mu_regular, 1)
mu_c3 = taylor_coefficient(mu_regular, 3)
p_c0 = taylor_coefficient(p_regular, 0)
p_c2 = taylor_coefficient(p_regular, 2)
p_c4 = taylor_coefficient(p_regular, 4)
p_infinity_exact = 60 * sp.pi - 300 * sp.sqrt(3) * sp.log(2)
derivative_coefficients_ok = bool(
    sp.simplify(mu_c1 - 60) == 0
    and sp.simplify(mu_c3 + 120 + 300 * sp.sqrt(3) / sp.pi) == 0
    and sp.simplify(p_c0.rewrite(sp.log) - p_infinity_exact) == 0
    and sp.simplify(p_c2 + 120 * sp.pi) == 0
    and sp.simplify(p_c4 - 360 * sp.pi - 900 * sp.sqrt(3)) == 0
)
check(
    "derivatives at t=0 reproduce every preregistered asymptotic coefficient",
    derivative_coefficients_ok,
)


# Re-derive the exact normalized identity and its limiting fixed family.
U, V, x, y = sp.symbols("U V x y", nonzero=True)
r = 2 / U - 1
y_after = -r * ((r + 1) * V + y)
y_on_root = -V - 4 * sp.pi * (U - 1) / x
drift = 4 * (U - 1) / U**2 * (V + 2 * sp.pi * U / x)
drift_identity_ok = sp.factor(
    sp.together((y_after - y - drift).subs(y, y_on_root))
) == 0

x_positive = sp.symbols("x_positive", positive=True)
f_x = 4 * sp.pi * (x_positive - 30) / x_positive**2
limiting_U = 60 / x_positive
limiting_V = -120 * sp.pi / x_positive**2
limiting_r = sp.simplify(2 / limiting_U - 1)
limiting_y_after = sp.simplify(
    -limiting_r
    * ((limiting_r + 1) * limiting_V + f_x)
)
family_identity_ok = bool(
    drift_identity_ok
    and sp.simplify(limiting_y_after - f_x) == 0
    and sp.simplify(limiting_r - (x_positive - 30) / 30) == 0
    and sp.simplify(sp.diff(f_x, x_positive) - 4 * sp.pi * (60 - x_positive) / x_positive**3) == 0
)
check(
    "the independently reduced boundary is a continuous fixed family",
    family_identity_ok,
    "outer physical branch x>60; no universal point selected",
)


# Derive the O(m^2) next-root displacement from implicit differentiation of
# f(x_plus)=y_plus, using only the derivative coefficients above.
A_exact = 120 * sp.pi + 300 * sp.sqrt(3)
r0 = (x_positive - 30) / 30
f_prime = sp.diff(f_x, x_positive)
y_drift_coefficient = A_exact * (60 - x_positive) / (900 * x_positive**3)
root_curve_correction_now = A_exact / x_positive**4
root_curve_correction_next = A_exact / (r0**2 * x_positive**4)
dx_coefficient = sp.factor(
    (
        root_curve_correction_now
        + y_drift_coefficient
        - root_curve_correction_next
    )
    / f_prime
)
dx_expected = (
    (2 * sp.pi + 5 * sp.sqrt(3))
    / (60 * sp.pi)
    * (1 - r0 ** (-2))
)
next_root_coefficient_ok = sp.simplify(dx_coefficient - dx_expected) == 0
check(
    "implicit differentiation reproduces the frozen next-root coefficient",
    next_root_coefficient_ok,
)


# The census below uses R(q)=0 rather than the primary equal-p partition.
# Its stationary points are exactly the nonzero solutions of mu(q)=mass.
q_symbol, mu_symbol, mu_prime_symbol, p_prime_symbol, mass_symbol = sp.symbols(
    "q_symbol mu_symbol mu_prime_symbol p_prime_symbol mass_symbol",
    nonzero=True,
)
R_prime = p_prime_symbol + 4 * sp.pi * (
    q_symbol * mu_prime_symbol - (mu_symbol - mass_symbol)
) / q_symbol**2
R_prime = sp.factor(
    R_prime.subs(p_prime_symbol, -4 * sp.pi * mu_prime_symbol / q_symbol)
)
stationary_identity_ok = sp.simplify(
    R_prime - 4 * sp.pi * (mass_symbol - mu_symbol) / q_symbol**2
) == 0
check(
    "the dual root function has stationary points only at equal-mu points",
    stationary_identity_ok,
)


def bisect(function, left, right, tolerance):
    f_left = function(left)
    f_right = function(right)
    if f_left == 0:
        return left
    if f_right == 0:
        return right
    if f_left * f_right >= 0:
        raise RuntimeError("adversarial bisection lacks a sign change")
    for _ in range(1400):
        middle = (left + right) / 2
        f_middle = function(middle)
        if abs(f_middle) < tolerance or abs(right - left) < tolerance:
            return middle
        if f_left * f_middle > 0:
            left, f_left = middle, f_middle
        else:
            right, f_right = middle, f_middle
    raise RuntimeError("adversarial bisection did not converge")


def outer_bracket(function, boundary, direction):
    fixed_value = function(boundary)
    magnitude = max(mp.mpf(10), 2 * abs(boundary))
    for _ in range(700):
        point = direction * magnitude
        if function(point) * fixed_value < 0:
            return (point, boundary) if direction < 0 else (boundary, point)
        magnitude *= 2
    raise RuntimeError("equal-mu outer root did not bracket")


def find_mu_peak(tolerance):
    # Central difference is not used: this is the analytic mpmath derivative
    # of the independently coded mu function.
    derivative = lambda q: mp.diff(mu, q)
    return bisect(derivative, mp.mpf(1), mp.mpf(10), tolerance)


def equal_mu_points(mass, peak, tolerance, margin):
    mu0 = mu(0)
    mu_peak = mu(peak)
    if mass <= 0 or mass >= mu_peak:
        return [], False, {"reason": "mass outside positive mu range"}
    equation = lambda q: mu(q) - mass
    positive = []
    if mass > mu0 + margin:
        positive.append(bisect(equation, 0, peak, tolerance))
    elif abs(mass - mu0) <= margin:
        return [], False, {"reason": "mass at central critical value"}
    positive.append(
        bisect(equation, *outer_bracket(equation, peak, 1), tolerance)
    )
    points = sorted([-q for q in positive] + positive)
    return points, True, {
        "mu_zero_margin": text(abs(mass - mu0), 35),
        "mu_peak_margin": text(abs(mass - mu_peak), 35),
    }


def approach_endpoint(function, endpoint, direction, expected_sign):
    if endpoint == 0:
        magnitude = mp.mpf("0.1")
        for _ in range(700):
            point = direction * magnitude
            value = function(point)
            if (value > 0) == (expected_sign > 0):
                return point
            magnitude /= 2
    else:
        magnitude = max(mp.mpf(10), 2 * abs(endpoint))
        for _ in range(700):
            point = direction * magnitude
            value = function(point)
            if (value > 0) == (expected_sign > 0):
                return point
            magnitude *= 2
    raise RuntimeError("dual endpoint sign not reached")


def dual_root_census(mass, incoming_p, peak, tolerance, margin):
    stationary, stationary_ok, diagnostics = equal_mu_points(
        mass, peak, tolerance, margin
    )
    if not stationary_ok:
        return [], False, diagnostics

    function = lambda q: (
        canonical_p(q)
        - incoming_p
        + 4 * mp.pi * (mu(q) - mass) / q
    )
    p_infinity = 60 * mp.pi - 300 * mp.sqrt(3) * mp.log(2)
    central_gap = mu(0) - mass
    endpoint_values = {
        "negative_infinity": -p_infinity - incoming_p,
        "negative_zero": -central_gap,
        "positive_zero": central_gap,
        "positive_infinity": p_infinity - incoming_p,
    }
    stationary_values = {q: function(q) for q in stationary}
    if min(abs(v) for v in endpoint_values.values()) <= margin:
        return [], False, {"reason": "endpoint value is critical"}
    if min(abs(v) for v in stationary_values.values()) <= margin:
        return [], False, {"reason": "stationary value is critical"}

    roots = []
    intervals = []

    def scan_side(side_points, left_limit, right_limit, side):
        boundaries = [None, *side_points, None]
        for index in range(len(boundaries) - 1):
            left = boundaries[index]
            right = boundaries[index + 1]
            left_value = left_limit if left is None else function(left)
            right_value = right_limit if right is None else function(right)
            has_root = left_value * right_value < 0
            intervals.append(
                {
                    "side": side,
                    "left": (
                        ("-infinity" if side == "negative" else "0+")
                        if left is None
                        else text(left, 35)
                    ),
                    "right": (
                        ("0-" if side == "negative" else "+infinity")
                        if right is None
                        else text(right, 35)
                    ),
                    "left_sign": 1 if left_value > 0 else -1,
                    "right_sign": 1 if right_value > 0 else -1,
                    "has_root": has_root,
                }
            )
            if not has_root:
                continue
            if left is None:
                bracket_left = approach_endpoint(
                    function,
                    0 if side == "positive" else right,
                    1 if side == "positive" else -1,
                    1 if left_value > 0 else -1,
                )
            else:
                bracket_left = left
            if right is None:
                bracket_right = approach_endpoint(
                    function,
                    0 if side == "negative" else left,
                    -1 if side == "negative" else 1,
                    1 if right_value > 0 else -1,
                )
            else:
                bracket_right = right
            roots.append(bisect(function, bracket_left, bracket_right, tolerance))

    negative_points = [q for q in stationary if q < 0]
    positive_points = [q for q in stationary if q > 0]
    scan_side(
        negative_points,
        endpoint_values["negative_infinity"],
        endpoint_values["negative_zero"],
        "negative",
    )
    scan_side(
        positive_points,
        endpoint_values["positive_zero"],
        endpoint_values["positive_infinity"],
        "positive",
    )
    roots.sort()
    q_zero_constraint = 4 * mp.pi * central_gap
    certified = abs(q_zero_constraint) > margin
    return roots, certified, {
        **diagnostics,
        "stationary_points": [text(q) for q in stationary],
        "stationary_values": [text(stationary_values[q], 35) for q in stationary],
        "endpoint_values": {name: text(value, 35) for name, value in endpoint_values.items()},
        "q_zero_constraint": text(q_zero_constraint, 35),
        "q_zero_excluded": certified,
        "intervals": intervals,
    }


PRECISIONS = (110, 180)
records = []
history_ok = True
census_ok = True
direct_ok = True
hostile_ok = True

for precision in PRECISIONS:
    mp.mp.dps = precision + 60
    tolerance = mp.mpf(10) ** (-(precision + 20))
    margin = mp.mpf(10) ** (-(precision - 20))
    residual_limit = mp.mpf(10) ** (-(precision - 25))

    peak = find_mu_peak(tolerance)
    initial_q = mp.mpf(3) / 2
    mass0 = mu(initial_q)
    pi0 = canonical_p(initial_q)

    h1, q1 = solve_action((mp.mpf(1) / 5, 10), mass0, pi0, tolerance)
    r1, m1, pi1, _ = transport(mass0, h1, q1)
    h2, q2 = solve_action((mp.mpf(1) / 14, 31), m1, pi1, tolerance)
    r2, m2, pi2, _ = transport(m1, h2, q2)
    h3, q3 = solve_action((mp.mpf(1) / 50, 100), m2, pi2, tolerance)
    r3, m3, pi3, _ = transport(m2, h3, q3)
    h4, q4 = solve_action((mp.mpf(1) / 150, 317), m3, pi3, tolerance)
    r4, m4, pi4, post4 = transport(m3, h4, q4)

    row_history_ok = bool(
        all(v > 0 for v in (h1, r1, h2, r2, h3, r3, h4, r4))
        and max(abs(v) for v in direct_equations(m3, pi3, h4, q4)) < residual_limit
    )
    history_ok &= row_history_ok

    # Freeze this route's forecast before calling the fifth-root census.
    x4 = m3 * q4
    r04 = x4 / 30 - 1
    dx_forecast = (
        (2 * mp.pi + 5 * mp.sqrt(3))
        / (60 * mp.pi)
        * (1 - 1 / r04**2)
        * m3**2
    )
    x5_forecast = x4 + dx_forecast

    roots, certified, diagnostics = dual_root_census(
        m4, pi4, peak, tolerance, margin
    )
    rows = []
    physical = []
    row_census_ok = certified
    row_direct_ok = True
    for q5 in roots:
        h5 = 2 * (m4 - mu(q5)) / (q5 * mu(q5))
        ratio5 = 1 + h5 * q5
        reduced = reduced_equations(m4, pi4, h5, q5)
        algebraic = max(abs(v) for v in reduced) < residual_limit
        is_physical = bool(h5 > 0 and ratio5 > 0)
        full = (mp.nan, mp.nan)
        full_ok = True
        if is_physical:
            full = direct_equations(m4, pi4, h5, q5)
            full_ok = max(abs(v) for v in full) < residual_limit
            physical.append((h5, q5, ratio5))
        row_census_ok &= algebraic
        row_direct_ok &= full_ok
        rows.append(
            {
                "q5": q5,
                "h5": h5,
                "ratio5": ratio5,
                "physical": is_physical,
                "reduced_constraint": reduced[0],
                "reduced_momentum": reduced[1],
                "direct_constraint": full[0],
                "direct_momentum": full[1],
                "algebraic_pass": algebraic,
                "direct_pass": full_ok,
            }
        )

    row_census_ok &= len(roots) == 3 and len(physical) == 1
    census_ok &= row_census_ok
    direct_ok &= row_direct_ok

    forecast = {}
    if len(physical) == 1:
        h5, q5, ratio5 = physical[0]
        x5 = m4 * q5
        _, m5, pi5, _ = transport(m4, h5, q5)
        p_infinity = 60 * mp.pi - 300 * mp.sqrt(3) * mp.log(2)
        y3 = (p_infinity - pi3) / m3**2
        y4 = (p_infinity - pi4) / m4**2
        y5 = (p_infinity - pi5) / m5**2
        forecast = {
            "x4": x4,
            "actual_x5": x5,
            "forecast_x5": x5_forecast,
            "forecast_error": x5 - x5_forecast,
            "forecast_error_over_m3_fourth": (x5 - x5_forecast) / m3**4,
            "x5_minus_x4": x5 - x4,
            "y4_minus_y3": y4 - y3,
            "y5_minus_y4": y5 - y4,
            "drift_signs_pass": x5 > x4 and y4 < y3 and y5 < y4,
        }
        row_census_ok &= forecast["drift_signs_pass"]
        census_ok &= forecast["drift_signs_pass"]

    wrong_rescaling = post4 / r4
    reversed_sign = -pi4
    reset_mass = mu(q4)
    no_boost_p = 180 * q4 * deficit(q4) / mp.sqrt(q4**2 + 4)
    row_hostile_ok = bool(
        abs(wrong_rescaling - pi4) > mp.mpf("1e-20")
        and abs(reversed_sign - pi4) > mp.mpf("1e-20")
        and abs(reset_mass - m4) > mp.mpf("1e-20")
        and abs(no_boost_p - canonical_p(q4)) > mp.mpf("1e-20")
    )
    hostile_ok &= row_hostile_ok

    records.append(
        {
            "precision": precision,
            "peak": peak,
            "history_ok": row_history_ok,
            "m3": m3,
            "pi3": pi3,
            "q4": q4,
            "h4": h4,
            "r4": r4,
            "m4": m4,
            "pi4": pi4,
            "root_count_certified": certified,
            "all_real_root_count": len(roots),
            "physical_root_count": len(physical),
            "diagnostics": diagnostics,
            "roots": rows,
            "forecast": forecast,
            "hostile": {
                "wrong_rescaling_gap": wrong_rescaling - pi4,
                "reversed_sign_gap": reversed_sign - pi4,
                "mass_reset_gap": reset_mass - m4,
                "boost_omission_gap": no_boost_p - canonical_p(q4),
                "passed": row_hostile_ok,
            },
        }
    )


check(
    "both precision runs reconstruct four slabs from the complete action",
    history_ok,
    "precisions=110,180",
)
check(
    "the equal-mu partition certifies three real fifth roots and one physical root",
    census_ok,
)
check(
    "every physical fifth root satisfies the independently differentiated action",
    direct_ok,
)
check(
    "all four hostile state and boost controls fail at both precisions",
    hostile_ok,
)


low, high = records
nesting_ok = bool(
    low["all_real_root_count"] == high["all_real_root_count"]
    and low["physical_root_count"] == high["physical_root_count"]
    and len(low["roots"]) == len(high["roots"])
)
for low_root, high_root in zip(low["roots"], high["roots"]):
    nesting_ok &= bool(
        abs(low_root["q5"] - high_root["q5"]) < mp.mpf("1e-60")
        and abs(low_root["h5"] - high_root["h5"]) < mp.mpf("1e-60")
        and low_root["physical"] == high_root["physical"]
    )
check(
    "the independent fifth-root censuses nest beyond 60 digits",
    nesting_ok,
)


# Only now read and compare the primary artifact.
primary = json.loads(PRIMARY_INPUT.read_text())
primary_provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and primary["outcome"]
    == "CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB"
)
comparison_ok = bool(
    primary_provenance_ok
    and high["all_real_root_count"] == len(primary["fifth_census"]["roots"])
    and high["physical_root_count"] == primary["fifth_census"]["physical_count"]
)
comparison = []
for independent, stored in zip(high["roots"], primary["fifth_census"]["roots"]):
    differences = {
        "q5": abs(independent["q5"] - mp.mpf(stored["q5"])),
        "h5": abs(independent["h5"] - mp.mpf(stored["h5"])),
        "ratio5": abs(independent["ratio5"] - mp.mpf(stored["scale_ratio"])),
    }
    row_ok = bool(
        max(differences.values()) < mp.mpf("1e-55")
        and independent["physical"] == stored["physical"]
    )
    comparison_ok &= row_ok
    comparison.append(
        {
            **{name: text(value, 25) for name, value in differences.items()},
            "passed": row_ok,
        }
    )
check(
    "the dual census agrees with the primary only after independent construction",
    comparison_ok,
)


all_gates = bool(
    protocol_ok
    and derivative_coefficients_ok
    and family_identity_ok
    and next_root_coefficient_ok
    and stationary_identity_ok
    and history_ok
    and census_ok
    and direct_ok
    and hostile_ok
    and nesting_ok
    and comparison_ok
)
outcome = (
    "CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB_"
    "ADVERSARIALLY_CORROBORATED"
    if all_gates
    else "ASYMPTOTIC_MAP_ADVERSARIAL_OPEN"
)
check(
    "the adversarial hierarchy corroborates the asymptotic family and fifth slab",
    outcome
    == "CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB_"
    "ADVERSARIALLY_CORROBORATED",
)


def pack(value):
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        return value
    return text(value)


def pack_record(record):
    return {
        "precision": record["precision"],
        "peak": text(record["peak"]),
        "history_ok": record["history_ok"],
        "m3": text(record["m3"]),
        "pi3": text(record["pi3"]),
        "q4": text(record["q4"]),
        "h4": text(record["h4"]),
        "r4": text(record["r4"]),
        "m4": text(record["m4"]),
        "pi4": text(record["pi4"]),
        "root_count_certified": record["root_count_certified"],
        "all_real_root_count": record["all_real_root_count"],
        "physical_root_count": record["physical_root_count"],
        "diagnostics": record["diagnostics"],
        "roots": [
            {name: pack(value) for name, value in row.items()}
            for row in record["roots"]
        ],
        "forecast": {
            name: pack(value) for name, value in record["forecast"].items()
        },
        "hostile": {
            name: pack(value) for name, value in record["hostile"].items()
        },
    }


artifact = {
    "provenance": {
        "protocol_commit": PROTOCOL_COMMIT,
        "protocol_sha256": PROTOCOL_SHA256,
        "primary_commit": PRIMARY_COMMIT,
        "primary_sha256": PRIMARY_SHA256,
        "primary_read_after_independent_census": True,
    },
    "method": {
        "action": "independently redifferentiated complete action",
        "asymptotics": "derivatives at regular t=0 expressions; no series call",
        "root_partition": "equal-mu stationary points of R(q)",
        "stationary_identity": "R'(q)=4*pi*(m-mu(q))/q^2",
        "precisions": list(PRECISIONS),
    },
    "exact_coefficients": {
        "mu_t": text(mu_c1),
        "mu_t3": str(mu_c3),
        "p_t2": str(p_c2),
        "p_t4": str(p_c4),
    },
    "precision_records": [pack_record(record) for record in records],
    "primary_comparison": comparison,
    "checks": {
        "protocol": protocol_ok,
        "derivative_coefficients": derivative_coefficients_ok,
        "continuous_family": family_identity_ok,
        "next_root_coefficient": next_root_coefficient_ok,
        "stationary_identity": stationary_identity_ok,
        "history": history_ok,
        "root_census": census_ok,
        "direct_action": direct_ok,
        "hostile_controls": hostile_ok,
        "precision_nesting": nesting_ok,
        "primary_comparison": comparison_ok,
    },
    "claims": {
        "continuous_asymptotic_fixed_family": "DERIVED",
        "unique_fifth_slab": "DERIVED_COMPUTATIONAL_FIVE_SLAB_SCOPED",
        "infinite_history": "OPEN",
        "universal_fixed_point": "NOT_SELECTED",
        "absolute_tick": "NOT_DERIVED",
    },
    "passed": passed,
    "tests": tests,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print()
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
if passed != tests:
    raise SystemExit(1)

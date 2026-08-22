#!/usr/bin/env python3
"""Scale-free finite-height map and out-of-sample fifth-slab census."""

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
FOURTH_INPUT = HERE / "gravity_600cell_finite_height_fourth_slab.json"
FOURTH_ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_finite_height_fourth_slab_adversarial.json"
)
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
OUTPUT = HERE / "gravity_600cell_finite_height_asymptotic_map.json"

PROTOCOL_COMMIT = "c39f9ca"
PROTOCOL_SHA256 = (
    "c9ce01a119be1e66b46796fa16cafd3d81b4ac0d732edf2a3a4c9a936bfec25f"
)
FOURTH_SHA256 = (
    "cf322cf0d60668d8f3f58e251425c9ad6bf43b112f22f9f3aebbc28f86212468"
)
FOURTH_ADVERSARIAL_SHA256 = (
    "ac1ed19fd72549cf7cd054107d921e2819580704391fbc294a55e106a8f7a1bd"
)
CLASSIFICATION_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)

mp.mp.dps = 180
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


fourth = json.loads(FOURTH_INPUT.read_text())
fourth_adversarial = json.loads(FOURTH_ADVERSARIAL_INPUT.read_text())
classification = json.loads(CLASSIFICATION_INPUT.read_text())
provenance_ok = bool(
    digest(PROTOCOL) == PROTOCOL_SHA256
    and digest(FOURTH_INPUT) == FOURTH_SHA256
    and digest(FOURTH_ADVERSARIAL_INPUT) == FOURTH_ADVERSARIAL_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and fourth["outcome"] == "SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB"
    and fourth_adversarial["outcome"]
    == "SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB_ADVERSARIALLY_CORROBORATED"
)
check(
    "the asymptotic protocol and accepted four-slab history are frozen",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


# Abstract exact scale-free map.
U_ABS, V_ABS, X_ABS = sp.symbols("U V x", nonzero=True)
Y_ROOT = -V_ABS - 4 * sp.pi * (U_ABS - 1) / X_ABS
R_ABS = 2 / U_ABS - 1
Y_PLUS = -R_ABS * ((R_ABS + 1) * V_ABS + Y_ROOT)
DRIFT_TARGET = (
    4
    * (U_ABS - 1)
    / U_ABS**2
    * (V_ABS + 2 * sp.pi * U_ABS / X_ABS)
)
exact_drift_residual = sp.factor(
    sp.together(Y_PLUS - Y_ROOT - DRIFT_TARGET)
)
exact_map_ok = exact_drift_residual == 0
check(
    "the exact canonical relation gives the preregistered scale-free drift identity",
    exact_map_ok,
)


# Derive asymptotic coefficients through an independent variable t=1/q.
t = sp.symbols("t", positive=True)
q_t = 1 / t
epsilon_t = 2 * sp.pi - 5 * sp.acos(
    (q_t**2 + 2) / (2 * (q_t**2 + 3))
)
mu_t = 180 * epsilon_t / (sp.pi * sp.sqrt(q_t**2 + 4))
p_t = (
    180 * q_t * epsilon_t / sp.sqrt(q_t**2 + 4)
    - 600
    * sp.sqrt(3)
    * sp.asinh(q_t / sp.sqrt(8 * (q_t**2 + 3)))
)
p_limit = sp.limit(p_t, t, 0, dir="+")
p_infinity_symbolic = 60 * sp.pi - 300 * sp.sqrt(3) * sp.log(2)
p_limit_ok = sp.simplify(
    p_limit.rewrite(sp.log) - p_infinity_symbolic
) == 0

mu_series = sp.series(mu_t, t, 0, 6).removeO().expand()
p_delta_series = sp.series(p_t - p_limit, t, 0, 6).removeO().expand()
b_mu = -120 - 300 * sp.sqrt(3) / sp.pi
d_p = 360 * sp.pi + 900 * sp.sqrt(3)
series_ok = bool(
    p_limit_ok
    and sp.simplify(mu_series.coeff(t, 1) - 60) == 0
    and sp.simplify(mu_series.coeff(t, 3) - b_mu) == 0
    and sp.simplify(p_delta_series.coeff(t, 2) + 120 * sp.pi) == 0
    and sp.simplify(p_delta_series.coeff(t, 4) - d_p) == 0
)
check(
    "the t=1/q expansion gives every frozen U and V coefficient",
    series_ok,
)


# Root curve, leading drift, and continuous limiting physical family.
x = sp.symbols("x", positive=True)
A = 120 * sp.pi + 300 * sp.sqrt(3)
f = 4 * sp.pi * (x - 30) / x**2
f_prime = sp.factor(sp.diff(f, x))
root_m2 = sp.factor(-(d_p + 4 * sp.pi * b_mu))
bracket_m2 = sp.factor(d_p + 2 * sp.pi * b_mu)
leading_drift = sp.factor(
    (x * (60 - x) / 900) * bracket_m2 / x**4
)
expected_drift = A * (60 - x) / (900 * x**3)
family_ok = bool(
    sp.simplify(root_m2 - A) == 0
    and sp.simplify(bracket_m2 - A) == 0
    and sp.simplify(leading_drift - expected_drift) == 0
    and sp.simplify(f_prime - 4 * sp.pi * (60 - x) / x**3) == 0
)
check(
    "the compactified physical branch is a continuous fixed family, not one point",
    family_ok,
    "physical boundary: x>60, 0<y<pi/30",
)


# Derive the next physical root drift by coefficient matching.
r0 = (x - 30) / 30
k_from_matching = sp.factor(
    (
        A / x**4
        + A * (60 - x) / (900 * x**3)
        - A / (r0**2 * x**4)
    )
    / f_prime
)
k_target = (
    (2 * sp.pi + 5 * sp.sqrt(3))
    / (60 * sp.pi)
    * (1 - r0 ** (-2))
)
root_drift_ok = sp.simplify(k_from_matching - k_target) == 0
check(
    "coefficient matching gives the frozen next-root drift",
    root_drift_ok,
)


# Complete action and numerical canonical reconstruction.
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
P_POST_EXACT = L_PLUS * sp.diff(ACTION, L_PLUS) / 2
f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), F_EXACT, "mpmath")
pre_numeric = sp.lambdify(
    (L_MINUS, L_PLUS, RHO, MASS), P_PRE_EXACT, "mpmath"
)
post_numeric = sp.lambdify(
    (L_MINUS, L_PLUS, RHO, MASS), P_POST_EXACT, "mpmath"
)


def epsilon(value):
    return 2 * mp.pi - 5 * mp.acos(
        (value**2 + 2) / (2 * (value**2 + 3))
    )


def mu(value):
    return 180 * epsilon(value) / (mp.pi * mp.sqrt(value**2 + 4))


def momentum(value):
    return (
        180 * value * epsilon(value) / mp.sqrt(value**2 + 4)
        - 600
        * mp.sqrt(3)
        * mp.asinh(value / mp.sqrt(8 * (value**2 + 3)))
    )


def momentum_no_boost(value):
    return 180 * value * epsilon(value) / mp.sqrt(value**2 + 4)


def direct_residuals(mass, incoming_p, height, slope):
    endpoint = 1 + height * slope
    return (
        2 * f_numeric(1, endpoint, height**2, mass) / height,
        pre_numeric(1, endpoint, height**2, mass) - incoming_p,
    )


def reduced_residuals(mass, incoming_p, height, slope):
    return (
        8 * mp.pi * (mu(slope) - mass)
        + 4 * mp.pi * height * slope * mu(slope),
        momentum(slope)
        - incoming_p
        - 2 * mp.pi * height * mu(slope),
    )


def solve_direct(seed, mass, incoming_p):
    return mp.findroot(
        lambda height, slope: direct_residuals(
            mass, incoming_p, height, slope
        ),
        seed,
        tol=mp.mpf("1e-150"),
        maxsteps=120,
    )


def advance_state(mass, height, slope):
    ratio = 1 + height * slope
    p_post = post_numeric(1, ratio, height**2, mass)
    next_mass = mass / ratio
    next_p = p_post / ratio**2
    closed_p = (
        momentum(slope)
        + 2 * mp.pi * height * mu(slope) / ratio
    )
    return ratio, next_mass, next_p, p_post, abs(next_p - closed_p)


v = mp.mpf(3) / 2
m0 = mu(v)
pi0 = momentum(v)
h1, q1 = solve_direct((mp.mpf(1) / 5, mp.mpf(10)), m0, pi0)
r1, m1, pi1, _, error1 = advance_state(m0, h1, q1)
h2, q2 = solve_direct((mp.mpf(1) / 14, mp.mpf(31)), m1, pi1)
r2, m2, pi2, _, error2 = advance_state(m1, h2, q2)
h3, q3 = solve_direct((mp.mpf(1) / 50, mp.mpf(100)), m2, pi2)
r3, m3, pi3, _, error3 = advance_state(m2, h3, q3)
h4, q4 = solve_direct((mp.mpf(1) / 150, mp.mpf(317)), m3, pi3)
r4, m4, pi4, p_post4, error4 = advance_state(m3, h4, q4)

stored_fourth = next(
    root for root in fourth["census"]["roots"] if root["physical"]
)
history_ok = bool(
    all(value > 0 for value in (h1, r1, h2, r2, h3, r3, h4, r4))
    and max(error1, error2, error3, error4) < mp.mpf("1e-145")
    and abs(h4 - mp.mpf(stored_fourth["h4"])) < mp.mpf("1e-60")
    and abs(q4 - mp.mpf(stored_fourth["q4"])) < mp.mpf("1e-60")
    and max(abs(value) for value in direct_residuals(m3, pi3, h4, q4))
    < mp.mpf("1e-125")
)
check(
    "the complete action reconstructs the accepted four-slab history and fifth input",
    history_ok,
    f"m4={text(m4, 25)}, pi4={text(pi4, 25)}",
)


# Frozen fifth-root forecast, evaluated before the actual fifth root.
p_infinity = 60 * mp.pi - 300 * mp.sqrt(3) * mp.log(2)
x4 = m3 * q4
y3 = (p_infinity - pi3) / m3**2
y4 = (p_infinity - pi4) / m4**2
r04 = x4 / 30 - 1
coefficient4 = (
    (2 * mp.pi + 5 * mp.sqrt(3))
    / (60 * mp.pi)
    * (1 - 1 / r04**2)
)
forecast_dx = coefficient4 * m3**2
forecast_x5 = x4 + forecast_dx
forecast_q5 = forecast_x5 / m4
forecast_frozen_ok = bool(
    abs(
        forecast_dx
        - mp.mpf(
            "0.011186371823854069885088912079005399125421359177972"
        )
    )
    < mp.mpf("1e-50")
)
check(
    "the committed fourth state reproduces the frozen fifth-root forecast",
    forecast_frozen_ok,
    f"x5_forecast={text(forecast_x5, 30)}",
)


# Global equal-p root census.
v_star = mp.mpf(classification["thresholds"]["v_star"])
p_star = momentum(v_star)
SIGN_TOL = mp.mpf("1e-110")
ROOT_TOL = mp.mpf("1e-140")


def topology_sign(value):
    if value > SIGN_TOL:
        return 1
    if value < -SIGN_TOL:
        return -1
    return 0


def bisect(function, left, right):
    f_left = function(left)
    f_right = function(right)
    if f_left == 0:
        return left
    if f_right == 0:
        return right
    if f_left * f_right >= 0:
        raise RuntimeError("fifth-slab bisection lacks a sign change")
    for _ in range(1200):
        middle = (left + right) / 2
        f_middle = function(middle)
        if (
            f_middle == 0
            or abs(f_middle) < ROOT_TOL
            or abs(right - left) < ROOT_TOL
        ):
            return middle
        if f_left * f_middle > 0:
            left, f_left = middle, f_middle
        else:
            right, f_right = middle, f_middle
    raise RuntimeError("fifth-slab bisection did not converge")


def outer_bracket(function, boundary, direction):
    magnitude = max(mp.mpf(10), 2 * abs(boundary))
    fixed = function(boundary)
    for _ in range(600):
        point = direction * magnitude
        value = function(point)
        if value * fixed < 0:
            return (point, boundary) if direction < 0 else (boundary, point)
        magnitude *= 2
    raise RuntimeError("fifth-slab outer root did not bracket")


def stationary_points(incoming_p):
    levels = (p_star, p_infinity, mp.mpf(0), -p_infinity, -p_star)
    margins = [abs(incoming_p - level) for level in levels]
    if min(margins) < SIGN_TOL:
        return [], False, margins
    function = lambda value: momentum(value) - incoming_p
    roots = []
    intervals = [
        (None, -v_star, -p_infinity, -p_star),
        (-v_star, mp.mpf(0), -p_star, mp.mpf(0)),
        (mp.mpf(0), v_star, mp.mpf(0), p_star),
        (v_star, None, p_star, p_infinity),
    ]
    for left, right, p_left, p_right in intervals:
        if not (min(p_left, p_right) < incoming_p < max(p_left, p_right)):
            continue
        bracket = (
            outer_bracket(function, right, -1)
            if left is None
            else outer_bracket(function, left, 1)
            if right is None
            else (left, right)
        )
        roots.append(bisect(function, *bracket))
    return sorted(roots), True, margins


def tail_point(function, direction, expected_sign, boundary):
    magnitude = max(mp.mpf(10), 2 * abs(boundary))
    for _ in range(600):
        point = direction * magnitude
        if topology_sign(function(point)) == expected_sign:
            return point
        magnitude *= 2
    raise RuntimeError("fifth-slab tail sign not reached")


def enumerate_roots(mass, incoming_p):
    stationary, stationary_ok, margins = stationary_points(incoming_p)
    if not stationary_ok:
        return [], False, {"reason": "critical momentum level"}
    function = lambda value: (
        4 * mp.pi * (mu(value) - mass)
        + value * (momentum(value) - incoming_p)
    )
    stationary_values = [function(value) for value in stationary]
    if any(abs(value) < SIGN_TOL for value in stationary_values):
        return [], False, {"reason": "stationary elimination root"}
    left_coefficient = -p_infinity - incoming_p
    right_coefficient = p_infinity - incoming_p
    left_sign = -topology_sign(left_coefficient)
    right_sign = topology_sign(right_coefficient)
    if left_sign == 0 or right_sign == 0:
        return [], False, {"reason": "cancelled linear tail"}

    boundaries = [None, *stationary, None]
    roots = []
    intervals = []
    for index in range(len(boundaries) - 1):
        left = boundaries[index]
        right = boundaries[index + 1]
        sign_left = left_sign if left is None else topology_sign(function(left))
        sign_right = right_sign if right is None else topology_sign(function(right))
        has_root = sign_left != sign_right
        intervals.append(
            {
                "left": "-infinity" if left is None else text(left, 35),
                "right": "+infinity" if right is None else text(right, 35),
                "left_sign": sign_left,
                "right_sign": sign_right,
                "has_root": has_root,
            }
        )
        if not has_root:
            continue
        bracket_left = (
            tail_point(function, -1, left_sign, right)
            if left is None
            else left
        )
        bracket_right = (
            tail_point(function, 1, right_sign, left)
            if right is None
            else right
        )
        roots.append(bisect(function, bracket_left, bracket_right))
    return roots, True, {
        "critical_level_margins": [text(value, 35) for value in margins],
        "stationary_points": [text(value) for value in stationary],
        "stationary_values": [text(value, 35) for value in stationary_values],
        "left_tail_coefficient": text(left_coefficient, 35),
        "right_tail_coefficient": text(right_coefficient, 35),
        "q_zero_value": text(function(0), 35),
        "q_zero_is_root": abs(function(0)) < SIGN_TOL,
        "intervals": intervals,
    }


roots, census_certified, diagnostics = enumerate_roots(m4, pi4)
root_rows = []
physical_roots = []
all_algebraic_ok = census_certified
all_direct_ok = True
for q5 in roots:
    h5 = (momentum(q5) - pi4) / (2 * mp.pi * mu(q5))
    r5 = 1 + h5 * q5
    reduced = reduced_residuals(m4, pi4, h5, q5)
    algebraic_ok = max(abs(value) for value in reduced) < mp.mpf("1e-110")
    physical = bool(h5 > 0 and r5 > 0)
    direct = (mp.nan, mp.nan)
    direct_ok = True
    junction = mp.nan
    if physical:
        direct = direct_residuals(m4, pi4, h5, q5)
        junction = pi4 - pre_numeric(1, r5, h5**2, m4)
        direct_ok = bool(
            max(abs(value) for value in direct) < mp.mpf("1e-100")
            and abs(junction) < mp.mpf("1e-100")
        )
        physical_roots.append((h5, q5, r5))
    all_algebraic_ok &= algebraic_ok
    all_direct_ok &= direct_ok
    root_rows.append(
        {
            "q5": text(q5),
            "h5": text(h5),
            "scale_ratio": text(r5),
            "constraint_residual": text(reduced[0], 35),
            "momentum_residual": text(reduced[1], 35),
            "direct_constraint_residual": text(direct[0], 35),
            "direct_momentum_residual": text(direct[1], 35),
            "junction_residual": text(junction, 35),
            "physical": physical,
            "algebraic_pass": algebraic_ok,
            "direct_pass": direct_ok,
        }
    )

check(
    "the equal-p partition classifies every real fifth-slab root",
    census_certified and all_algebraic_ok,
    f"all_real_roots={len(roots)}",
)
check(
    "every physical fifth root passes the complete action and junction",
    all_direct_ok,
)


forecast = {}
drift_signs_ok = False
if len(physical_roots) == 1:
    h5, q5, r5 = physical_roots[0]
    x5 = m4 * q5
    r5_out, m5, pi5, _, error5 = advance_state(m4, h5, q5)
    y5 = (p_infinity - pi5) / m5**2
    forecast_error = x5 - forecast_x5
    drift_signs_ok = bool(
        x5 > x4
        and y4 < y3
        and y5 < y4
        and error5 < mp.mpf("1e-145")
    )
    forecast = {
        "actual_x5": text(x5),
        "actual_q5": text(q5),
        "actual_h5": text(h5),
        "actual_r5": text(r5),
        "actual_y5": text(y5),
        "x5_minus_x4": text(x5 - x4),
        "y4_minus_y3": text(y4 - y3),
        "y5_minus_y4": text(y5 - y4),
        "forecast_x5": text(forecast_x5),
        "forecast_q5": text(forecast_q5),
        "x_forecast_error": text(forecast_error),
        "x_forecast_error_over_m3_fourth": text(forecast_error / m3**4),
    }

check(
    "the out-of-sample fifth root has the preregistered asymptotic drift signs",
    drift_signs_ok,
    f"physical_roots={len(physical_roots)}",
)


wrong_scale = p_post4 / r4
wrong_sign = -pi4
reset_mass = mu(q4)
boost_omission = momentum_no_boost(q4) - momentum(q4)
hostile_ok = bool(
    abs(wrong_scale - pi4) > mp.mpf("1e-20")
    and abs(wrong_sign - pi4) > mp.mpf("1e-20")
    and abs(reset_mass - m4) > mp.mpf("1e-20")
    and abs(boost_omission) > mp.mpf("1e-20")
)
check(
    "state-scaling, sign, mass-reset and boost-omission controls all fail",
    hostile_ok,
)


all_gates = bool(
    provenance_ok
    and exact_map_ok
    and series_ok
    and family_ok
    and root_drift_ok
    and history_ok
    and forecast_frozen_ok
    and census_certified
    and all_algebraic_ok
    and all_direct_ok
    and hostile_ok
)
if not all_gates:
    outcome = "ASYMPTOTIC_MAP_OPEN"
elif len(physical_roots) == 0:
    outcome = "ASYMPTOTIC_FAMILY_DERIVED_BUT_HISTORY_STOPS"
elif len(physical_roots) > 1:
    outcome = "FIFTH_SLAB_BRANCHES"
elif drift_signs_ok:
    outcome = "CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB"
else:
    outcome = "ASYMPTOTIC_MAP_REFUTED"

artifact = {
    "protocol": {
        "commit": PROTOCOL_COMMIT,
        "sha256": PROTOCOL_SHA256,
        "post_hoc_variable_disclosure": True,
        "fifth_slab_out_of_sample": True,
    },
    "symbolic_map": {
        "p_infinity": "60*pi-300*sqrt(3)*log(2)",
        "exact_root": "4*pi*(U-1)+x*(V+y)=0",
        "exact_r": "2/U-1",
        "exact_y_plus": "-r*((r+1)*V+y)",
        "exact_drift": "4*(U-1)/U^2*(V+2*pi*U/x)",
        "U_m0": "60/x",
        "U_m2": "(-120-300*sqrt(3)/pi)/x^3",
        "V_m0": "-120*pi/x^2",
        "V_m2": "(360*pi+900*sqrt(3))/x^4",
        "limiting_curve": "y=4*pi*(x-30)/x^2",
        "physical_boundary": "x>60",
        "limiting_map": "y_plus=y for a continuous family",
        "y_drift_m2": "A*(60-x)/(900*x^3)",
        "x_drift_m2": (
            "(2*pi+5*sqrt(3))/(60*pi)*(1-r0^-2)"
        ),
        "A": "120*pi+300*sqrt(3)",
    },
    "history": {
        "m3": text(m3),
        "pi3": text(pi3),
        "q4": text(q4),
        "h4": text(h4),
        "r4": text(r4),
        "x4": text(x4),
        "y3": text(y3),
        "m4": text(m4),
        "pi4": text(pi4),
        "y4": text(y4),
    },
    "frozen_forecast": {
        "coefficient": text(coefficient4),
        "x5_minus_x4": text(forecast_dx),
        "x5": text(forecast_x5),
        "q5": text(forecast_q5),
    },
    "fifth_census": {
        "certified": census_certified,
        "diagnostics": diagnostics,
        "roots": root_rows,
        "physical_count": len(physical_roots),
    },
    "forecast_comparison": forecast,
    "checks": {
        "provenance": provenance_ok,
        "exact_map": exact_map_ok,
        "series": series_ok,
        "continuous_family": family_ok,
        "next_root_drift": root_drift_ok,
        "four_slab_reconstruction": history_ok,
        "forecast_reproduction": forecast_frozen_ok,
        "root_census": census_certified and all_algebraic_ok,
        "direct_action": all_direct_ok,
        "drift_signs": drift_signs_ok,
        "hostile_controls": hostile_ok,
    },
    "claims": {
        "continuous_asymptotic_fixed_family": "DERIVED"
        if all_gates
        else "OPEN",
        "unique_fifth_slab": "DERIVED_COMPUTATIONAL_FIVE_SLAB_SCOPED"
        if len(physical_roots) == 1 and all_gates
        else "OPEN",
        "infinite_history": "OPEN",
        "universal_scale_ratio": "NOT_SELECTED_BY_LIMITING_MAP"
        if all_gates
        else "OPEN",
        "fundamental_tick": "NO",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print()
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
if passed != tests or outcome == "ASYMPTOTIC_MAP_OPEN":
    raise SystemExit(1)


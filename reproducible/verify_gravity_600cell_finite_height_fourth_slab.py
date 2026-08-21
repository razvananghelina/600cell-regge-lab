#!/usr/bin/env python3
"""Complete fourth-slab census along the unique three-slab history."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
THIRD_INPUT = HERE / "gravity_600cell_finite_height_third_slab.json"
THIRD_ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_finite_height_third_slab_adversarial.json"
)
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
OUTPUT = HERE / "gravity_600cell_finite_height_fourth_slab.json"

THIRD_SHA256 = (
    "6b0e92d031aa891fdc3e1b2045c35bd135a955bb1374c92f015dcd5727d3d8fc"
)
THIRD_ADVERSARIAL_SHA256 = (
    "df689f5360ace94d2212e1d71c799ed4e8019457d2702e989bf045ea566abda8"
)
CLASSIFICATION_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
PROTOCOL_COMMIT = "54c4554"

mp.mp.dps = 170
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


third = json.loads(THIRD_INPUT.read_text())
third_adversarial = json.loads(THIRD_ADVERSARIAL_INPUT.read_text())
classification = json.loads(CLASSIFICATION_INPUT.read_text())
provenance_ok = bool(
    digest(THIRD_INPUT) == THIRD_SHA256
    and digest(THIRD_ADVERSARIAL_INPUT) == THIRD_ADVERSARIAL_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and third["outcome"] == "ONE_SECOND_BRANCH_EXTENDS_UNIQUELY"
    and third_adversarial["outcome"]
    == "ONE_SECOND_BRANCH_EXTENDS_UNIQUELY_ADVERSARIALLY_CORROBORATED"
)
check(
    "the unique three-slab history and fourth-slab protocol are frozen",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


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
        tol=mp.mpf("1e-145"),
        maxsteps=100,
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


# Reconstruct all three accepted slabs from the complete action.
v = mp.mpf(3) / 2
mass0 = mu(v)
pi0 = momentum(v)
h1, q1 = solve_direct((mp.mpf(1) / 5, mp.mpf(10)), mass0, pi0)
r1, m1, pi1, _, error1 = advance_state(mass0, h1, q1)
h2, q2 = solve_direct((mp.mpf(1) / 14, mp.mpf(31)), m1, pi1)
r2, m2, pi2, _, error2 = advance_state(m1, h2, q2)
h3, q3 = solve_direct((mp.mpf(1) / 50, mp.mpf(100)), m2, pi2)
r3, m3, pi3, p_post3, error3 = advance_state(m2, h3, q3)

stored_third = next(
    root
    for root in third["third_slab_census"]["B"]["roots"]
    if root["physical"]
)
history_ok = bool(
    all(value > 0 for value in (h1, r1, h2, r2, h3, r3))
    and max(
        max(abs(value) for value in direct_residuals(mass0, pi0, h1, q1)),
        max(abs(value) for value in direct_residuals(m1, pi1, h2, q2)),
        max(abs(value) for value in direct_residuals(m2, pi2, h3, q3)),
    )
    < mp.mpf("1e-120")
    and abs(h3 - mp.mpf(stored_third["h3"])) < mp.mpf("1e-55")
    and abs(q3 - mp.mpf(stored_third["q3"])) < mp.mpf("1e-55")
)
check(
    "the complete action reconstructs the unique accepted three-slab history",
    history_ok,
)
recurrence_ok = max(error1, error2, error3) < mp.mpf("1e-140")
check(
    "the full-action and closed outgoing-state recurrences agree through slab three",
    recurrence_ok,
    f"third_error={text(error3, 20)}",
)


v_star = mp.mpf(classification["thresholds"]["v_star"])
p_star = momentum(v_star)
p_infinity = 60 * mp.pi - 300 * mp.sqrt(3) * mp.log(2)
SIGN_TOL = mp.mpf("1e-105")
ROOT_TOL = mp.mpf("1e-130")


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
        raise RuntimeError("fourth-slab bisection lacks a sign change")
    for _ in range(1000):
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
    raise RuntimeError("fourth-slab bisection did not converge")


def outer_bracket(function, boundary, direction):
    magnitude = max(mp.mpf(10), 2 * abs(boundary))
    fixed = function(boundary)
    for _ in range(500):
        point = direction * magnitude
        value = function(point)
        if value * fixed < 0:
            return (point, boundary) if direction < 0 else (boundary, point)
        magnitude *= 2
    raise RuntimeError("fourth-slab outer root did not bracket")


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
    for _ in range(500):
        point = direction * magnitude
        if topology_sign(function(point)) == expected_sign:
            return point
        magnitude *= 2
    raise RuntimeError("fourth-slab tail sign not reached")


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


roots, census_certified, diagnostics = enumerate_roots(m3, pi3)
root_rows = []
physical_count = 0
all_algebraic_ok = census_certified
all_direct_ok = True
for q4 in roots:
    h4 = (momentum(q4) - pi3) / (2 * mp.pi * mu(q4))
    r4 = 1 + h4 * q4
    reduced = reduced_residuals(m3, pi3, h4, q4)
    algebraic_ok = max(abs(value) for value in reduced) < mp.mpf("1e-100")
    physical = bool(h4 > 0 and r4 > 0)
    direct = (mp.nan, mp.nan)
    direct_ok = True
    junction = mp.nan
    if physical:
        direct = direct_residuals(m3, pi3, h4, q4)
        junction = pi3 - pre_numeric(1, r4, h4**2, m3)
        direct_ok = bool(
            max(abs(value) for value in direct) < mp.mpf("1e-90")
            and abs(junction) < mp.mpf("1e-90")
        )
        physical_count += 1
    all_algebraic_ok &= algebraic_ok
    all_direct_ok &= direct_ok
    root_rows.append(
        {
            "q4": text(q4),
            "h4": text(h4),
            "scale_ratio": text(r4),
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
    "the equal-p monotone proof classifies every real fourth-slab root",
    census_certified and all_algebraic_ok,
    f"all_real_roots={len(roots)}",
)
check(
    "every physical fourth slab passes the complete action and junction",
    all_direct_ok,
)


wrong_scale = p_post3 / r3
wrong_sign = -pi3
reset_mass = mu(q3)
hostile_ok = bool(
    abs(wrong_scale - pi3) > mp.mpf("1e-20")
    and abs(wrong_sign - pi3) > mp.mpf("1e-20")
    and abs(reset_mass - m3) > mp.mpf("1e-20")
)
check(
    "wrong scaling, sign and mass reset all change the fourth incoming state",
    hostile_ok,
)


if provenance_ok and history_ok and recurrence_ok and census_certified and all_algebraic_ok and all_direct_ok and hostile_ok:
    if physical_count == 0:
        outcome = "SURVIVING_HISTORY_STOPS_BEFORE_FOURTH_SLAB"
    elif physical_count == 1:
        outcome = "SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB"
    else:
        outcome = "SURVIVING_HISTORY_BRANCHES_AT_FOURTH_SLAB"
else:
    outcome = "FOURTH_SLAB_EXTENDIBILITY_OPEN"
check(
    "the preregistered hierarchy classifies fourth-slab extendibility",
    outcome != "FOURTH_SLAB_EXTENDIBILITY_OPEN",
    f"physical_roots={physical_count}",
)


artifact = {
    "provenance": {
        "third_slab_sha256": THIRD_SHA256,
        "third_slab_adversarial_sha256": THIRD_ADVERSARIAL_SHA256,
        "classification_sha256": CLASSIFICATION_SHA256,
        "protocol_commit": PROTOCOL_COMMIT,
    },
    "history": {
        "v": "3/2",
        "slab_1": {"h": text(h1), "q": text(q1), "ratio": text(r1)},
        "slab_2_B": {"h": text(h2), "q": text(q2), "ratio": text(r2)},
        "slab_3": {"h": text(h3), "q": text(q3), "ratio": text(r3)},
    },
    "fourth_incoming_state": {
        "m3": text(m3),
        "pi3": text(pi3),
        "recurrence_error": text(error3, 35),
    },
    "census": {
        "all_real_root_count": len(roots),
        "physical_root_count": physical_count,
        "root_count_certified": census_certified,
        "diagnostics": diagnostics,
        "roots": root_rows,
    },
    "hostile_controls": {
        "wrong_scale_gap": text(wrong_scale - pi3, 35),
        "wrong_sign_gap": text(wrong_sign - pi3, 35),
        "mass_reset_gap": text(reset_mass - m3, 35),
        "passed": hostile_ok,
    },
    "interpretation": {
        "label": "DERIVED COMPUTATIONAL, FOUR-SLAB SCOPED / STRUCTURAL",
        "indefinite_history": "OPEN",
        "fundamental_tick": False,
    },
    "tests": tests,
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
raise SystemExit(0 if passed == tests else 1)

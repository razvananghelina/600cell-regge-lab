#!/usr/bin/env python3
"""Equal-mu replication of the unique fourth-slab continuation."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
THIRD_ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_finite_height_third_slab_adversarial.json"
)
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_fourth_slab.json"
OUTPUT = HERE / "gravity_600cell_finite_height_fourth_slab_adversarial.json"

THIRD_ADVERSARIAL_SHA256 = (
    "df689f5360ace94d2212e1d71c799ed4e8019457d2702e989bf045ea566abda8"
)
PRIMARY_SHA256 = (
    "cf322cf0d60668d8f3f58e251425c9ad6bf43b112f22f9f3aebbc28f86212468"
)
PROTOCOL_COMMIT = "530b2f0"
PRIMARY_ARTIFACT_COMMIT = "7601c8f"

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


# Keep the primary fourth-slab artifact unread until the dual construction is
# complete.
third_adversarial = json.loads(THIRD_ADVERSARIAL_INPUT.read_text())
provenance_ok = bool(
    digest(THIRD_ADVERSARIAL_INPUT) == THIRD_ADVERSARIAL_SHA256
    and third_adversarial["outcome"]
    == "ONE_SECOND_BRANCH_EXTENDS_UNIQUELY_ADVERSARIALLY_CORROBORATED"
)
check(
    "the independently accepted three-slab history and dual protocol are frozen",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


q_s, u_s, up_s, pp_s, m_s = sp.symbols("q u up pp m", real=True)
R_prime = pp_s + 4 * sp.pi * (
    up_s * q_s - (u_s - m_s)
) / q_s**2
R_prime = sp.factor(R_prime.subs(pp_s, -4 * sp.pi * up_s / q_s))
dual_identity_ok = sp.simplify(
    R_prime - 4 * sp.pi * (m_s - u_s) / q_s**2
) == 0
check(
    "the dual fourth-root derivative is controlled only by equal-mu points",
    dual_identity_ok,
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


def k_function(x_value):
    return (
        10 * mp.sqrt(x_value + 4)
        - (x_value + 3)
        * mp.sqrt(3 * x_value + 8)
        * (
            2 * mp.pi
            - 5 * mp.acos((x_value + 2) / (2 * (x_value + 3)))
        )
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


def solve_direct(seed, mass, incoming_p, tolerance):
    return mp.findroot(
        lambda height, slope: direct_residuals(
            mass, incoming_p, height, slope
        ),
        seed,
        tol=tolerance,
        maxsteps=120,
    )


def advance(mass, height, slope):
    ratio = 1 + height * slope
    p_post = post_numeric(1, ratio, height**2, mass)
    return ratio, mass / ratio, p_post / ratio**2, p_post


def bisect(function, left, right, tolerance):
    f_left = function(left)
    f_right = function(right)
    if f_left == 0:
        return left
    if f_right == 0:
        return right
    if f_left * f_right >= 0:
        raise RuntimeError("dual fourth bisection lacks a sign change")
    for _ in range(1200):
        middle = (left + right) / 2
        f_middle = function(middle)
        if (
            f_middle == 0
            or abs(f_middle) < tolerance
            or abs(right - left) < tolerance
        ):
            return middle
        if f_left * f_middle > 0:
            left, f_left = middle, f_middle
        else:
            right, f_right = middle, f_middle
    raise RuntimeError("dual fourth bisection did not converge")


def outer_bracket(function, boundary, direction):
    fixed = function(boundary)
    magnitude = max(mp.mpf(10), 2 * abs(boundary))
    for _ in range(600):
        point = direction * magnitude
        value = function(point)
        if value * fixed < 0:
            return (point, boundary) if direction < 0 else (boundary, point)
        magnitude *= 2
    raise RuntimeError("dual fourth outer root did not bracket")


def equal_mu_points(mass, v_star, tolerance, margin):
    mu_zero = mu(0)
    mu_star = mu(v_star)
    if (
        mass <= 0
        or mass > mu_star
        or abs(mass - mu_zero) < margin
        or abs(mass - mu_star) < margin
    ):
        return [], False, {}
    function = lambda value: mu(value) - mass
    positive = []
    if mass > mu_zero:
        positive.append(bisect(function, 0, v_star, tolerance))
    positive.append(
        bisect(function, *outer_bracket(function, v_star, 1), tolerance)
    )
    points = sorted([-value for value in positive] + positive)
    return points, True, {
        "mu_zero_margin": text(abs(mass - mu_zero), 35),
        "mu_star_margin": text(abs(mass - mu_star), 35),
    }


def limit_point(function, direction, expected_sign, boundary, near_zero):
    magnitude = (
        min(mp.mpf("0.1"), abs(boundary) / 2)
        if near_zero
        else max(mp.mpf(10), 2 * abs(boundary))
    )
    for _ in range(600):
        point = direction * magnitude
        value = function(point)
        if (value > 0) == (expected_sign > 0):
            return point
        magnitude = magnitude / 2 if near_zero else magnitude * 2
    raise RuntimeError("dual fourth limit sign not reached")


def dual_census(mass, incoming_p, v_star, tolerance, margin):
    stationary, stationary_ok, diagnostics = equal_mu_points(
        mass, v_star, tolerance, margin
    )
    if not stationary_ok:
        return [], False, diagnostics
    function = lambda value: (
        momentum(value)
        - incoming_p
        + 4 * mp.pi * (mu(value) - mass) / value
    )
    p_infinity = 60 * mp.pi - 300 * mp.sqrt(3) * mp.log(2)
    mu_zero_gap = mu(0) - mass
    limit_sources = {
        "negative_infinity": -p_infinity - incoming_p,
        "negative_zero": -mu_zero_gap,
        "positive_zero": mu_zero_gap,
        "positive_infinity": p_infinity - incoming_p,
    }
    stationary_values = {value: function(value) for value in stationary}
    if (
        min(abs(value) for value in limit_sources.values()) < margin
        or min(abs(value) for value in stationary_values.values()) < margin
    ):
        return [], False, {"reason": "exceptional dual sign"}

    roots = []
    intervals = []

    def add_side(points, limit_signs, side):
        boundaries = [None, *points, None]
        for index in range(len(boundaries) - 1):
            left = boundaries[index]
            right = boundaries[index + 1]
            left_sign = (
                limit_signs[0]
                if left is None
                else (1 if function(left) > 0 else -1)
            )
            right_sign = (
                limit_signs[1]
                if right is None
                else (1 if function(right) > 0 else -1)
            )
            has_root = left_sign != right_sign
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
                    "left_sign": left_sign,
                    "right_sign": right_sign,
                    "has_root": has_root,
                }
            )
            if not has_root:
                continue
            if left is None:
                bracket_left = limit_point(
                    function,
                    -1 if side == "negative" else 1,
                    left_sign,
                    right,
                    near_zero=(side == "positive"),
                )
            else:
                bracket_left = left
            if right is None:
                bracket_right = limit_point(
                    function,
                    -1 if side == "negative" else 1,
                    right_sign,
                    left,
                    near_zero=(side == "negative"),
                )
            else:
                bracket_right = right
            roots.append(bisect(function, bracket_left, bracket_right, tolerance))

    negative = [value for value in stationary if value < 0]
    positive = [value for value in stationary if value > 0]
    add_side(
        negative,
        (
            1 if limit_sources["negative_infinity"] > 0 else -1,
            1 if limit_sources["negative_zero"] > 0 else -1,
        ),
        "negative",
    )
    add_side(
        positive,
        (
            1 if limit_sources["positive_zero"] > 0 else -1,
            1 if limit_sources["positive_infinity"] > 0 else -1,
        ),
        "positive",
    )
    roots.sort()
    q_zero_constraint = 8 * mp.pi * mu_zero_gap
    q_zero_excluded = abs(q_zero_constraint) > margin
    return roots, q_zero_excluded, {
        **diagnostics,
        "stationary_points": [text(value) for value in stationary],
        "stationary_values": [
            text(stationary_values[value], 35) for value in stationary
        ],
        "one_sided_limit_sources": {
            name: text(value, 35) for name, value in limit_sources.items()
        },
        "q_zero_constraint": text(q_zero_constraint, 35),
        "q_zero_excluded": q_zero_excluded,
        "intervals": intervals,
    }


TARGET_PRECISIONS = (110, 180)
precision_records = []
history_ok = True
census_ok = True
hostile_ok = True

for target_precision in TARGET_PRECISIONS:
    mp.mp.dps = target_precision + 60
    tolerance = mp.mpf(10) ** (-(target_precision + 20))
    margin = mp.mpf(10) ** (-(target_precision - 20))
    v_star = mp.sqrt(mp.findroot(k_function, (mp.mpf(5), mp.mpf(6))))

    v = mp.mpf(3) / 2
    mass0 = mu(v)
    pi0 = momentum(v)
    h1, q1 = solve_direct((mp.mpf(1) / 5, mp.mpf(10)), mass0, pi0, tolerance)
    r1, m1, pi1, _ = advance(mass0, h1, q1)
    h2, q2 = solve_direct((mp.mpf(1) / 14, mp.mpf(31)), m1, pi1, tolerance)
    r2, m2, pi2, _ = advance(m1, h2, q2)
    h3, q3 = solve_direct((mp.mpf(1) / 50, mp.mpf(100)), m2, pi2, tolerance)
    r3, m3, pi3, p_post3 = advance(m2, h3, q3)
    history_ok &= bool(
        all(value > 0 for value in (h1, r1, h2, r2, h3, r3))
        and max(abs(value) for value in direct_residuals(m2, pi2, h3, q3))
        < mp.mpf("1e-90")
    )

    roots, certified, diagnostics = dual_census(
        m3, pi3, v_star, tolerance, margin
    )
    root_rows = []
    physical_count = 0
    row_census_ok = certified
    for q4 in roots:
        h4 = 2 * (m3 - mu(q4)) / (q4 * mu(q4))
        r4 = 1 + h4 * q4
        reduced = reduced_residuals(m3, pi3, h4, q4)
        algebraic_ok = max(abs(value) for value in reduced) < mp.mpf("1e-90")
        physical = bool(h4 > 0 and r4 > 0)
        direct = (mp.nan, mp.nan)
        direct_ok = True
        if physical:
            direct = direct_residuals(m3, pi3, h4, q4)
            direct_ok = max(abs(value) for value in direct) < mp.mpf("1e-90")
            physical_count += 1
        row_census_ok &= algebraic_ok and direct_ok
        root_rows.append(
            {
                "q4": q4,
                "h4": h4,
                "ratio": r4,
                "constraint_residual": reduced[0],
                "momentum_residual": reduced[1],
                "direct_constraint_residual": direct[0],
                "direct_momentum_residual": direct[1],
                "physical": physical,
                "algebraic_pass": algebraic_ok,
                "direct_pass": direct_ok,
            }
        )
    census_ok &= row_census_ok

    wrong_scale = p_post3 / r3
    wrong_sign = -pi3
    reset_mass = mu(q3)
    row_hostile_ok = bool(
        abs(wrong_scale - pi3) > mp.mpf("1e-20")
        and abs(wrong_sign - pi3) > mp.mpf("1e-20")
        and abs(reset_mass - m3) > mp.mpf("1e-20")
    )
    hostile_ok &= row_hostile_ok
    precision_records.append(
        {
            "target_precision": target_precision,
            "history": {
                "h1": h1,
                "q1": q1,
                "h2": h2,
                "q2": q2,
                "h3": h3,
                "q3": q3,
                "r3": r3,
            },
            "m3": m3,
            "pi3": pi3,
            "root_count_certified": certified,
            "all_real_root_count": len(roots),
            "physical_root_count": physical_count,
            "diagnostics": diagnostics,
            "roots": root_rows,
            "hostile": {
                "wrong_scale_gap": wrong_scale - pi3,
                "wrong_sign_gap": wrong_sign - pi3,
                "mass_reset_gap": reset_mass - m3,
                "passed": row_hostile_ok,
            },
        }
    )

check(
    "both precision runs reconstruct the unique three-slab history from the action",
    history_ok,
    "precisions=110,180",
)
check(
    "the equal-mu proof certifies every algebraic and physical fourth root",
    census_ok,
)
check(
    "all outgoing-state hostile controls fail at both precisions",
    hostile_ok,
)


first, second = precision_records
nesting_ok = bool(
    first["all_real_root_count"] == second["all_real_root_count"]
    and first["physical_root_count"] == second["physical_root_count"]
    and len(first["roots"]) == len(second["roots"])
)
for first_root, second_root in zip(first["roots"], second["roots"]):
    nesting_ok &= bool(
        abs(first_root["q4"] - second_root["q4"]) < mp.mpf("1e-60")
        and abs(first_root["h4"] - second_root["h4"]) < mp.mpf("1e-60")
        and first_root["physical"] == second_root["physical"]
    )
check(
    "the dual fourth-root census nests beyond 60 digits",
    nesting_ok,
)


independent_unique = bool(
    second["all_real_root_count"] == 3
    and second["physical_root_count"] == 1
)
check(
    "the independent dual proof finds one physical fourth slab",
    independent_unique,
    f"all/physical={second['all_real_root_count']}/{second['physical_root_count']}",
)


# Only now compare with the primary artifact.
primary = json.loads(PRIMARY_INPUT.read_text())
primary_provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and primary["outcome"] == "SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB"
    and primary["passed"] == primary["tests"] == 7
)
comparison_ok = bool(
    primary_provenance_ok
    and second["all_real_root_count"] == primary["census"]["all_real_root_count"]
    and second["physical_root_count"] == primary["census"]["physical_root_count"]
    and len(second["roots"]) == len(primary["census"]["roots"])
)
comparison = []
for independent_root, primary_root in zip(
    second["roots"], primary["census"]["roots"]
):
    differences = {
        "q4": abs(independent_root["q4"] - mp.mpf(primary_root["q4"])),
        "h4": abs(independent_root["h4"] - mp.mpf(primary_root["h4"])),
        "ratio": abs(
            independent_root["ratio"] - mp.mpf(primary_root["scale_ratio"])
        ),
    }
    row_ok = bool(
        max(differences.values()) < mp.mpf("1e-55")
        and independent_root["physical"] == primary_root["physical"]
    )
    comparison_ok &= row_ok
    comparison.append(
        {
            **{name: text(value, 20) for name, value in differences.items()},
            "passed": row_ok,
        }
    )
state_differences = {
    "m3": abs(second["m3"] - mp.mpf(primary["fourth_incoming_state"]["m3"])),
    "pi3": abs(second["pi3"] - mp.mpf(primary["fourth_incoming_state"]["pi3"])),
}
comparison_ok &= max(state_differences.values()) < mp.mpf("1e-55")
check(
    "the dual census agrees with the equal-p primary only after construction",
    comparison_ok,
)


outcome = (
    "SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB_"
    "ADVERSARIALLY_CORROBORATED"
    if provenance_ok
    and dual_identity_ok
    and history_ok
    and census_ok
    and hostile_ok
    and nesting_ok
    and independent_unique
    and comparison_ok
    else "FOURTH_SLAB_EXTENDIBILITY_ADVERSARIAL_OPEN"
)
check(
    "the adversarial hierarchy corroborates the unique fourth slab",
    outcome
    == "SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB_"
    "ADVERSARIALLY_CORROBORATED",
)


def pack(value):
    return value if isinstance(value, bool) else text(value)


def pack_record(row):
    return {
        "target_precision": row["target_precision"],
        "history": {name: text(value) for name, value in row["history"].items()},
        "m3": text(row["m3"]),
        "pi3": text(row["pi3"]),
        "root_count_certified": row["root_count_certified"],
        "all_real_root_count": row["all_real_root_count"],
        "physical_root_count": row["physical_root_count"],
        "diagnostics": row["diagnostics"],
        "roots": [
            {name: pack(value) for name, value in root.items()}
            for root in row["roots"]
        ],
        "hostile": {name: pack(value) for name, value in row["hostile"].items()},
    }


artifact = {
    "provenance": {
        "third_slab_adversarial_sha256": THIRD_ADVERSARIAL_SHA256,
        "primary_fourth_slab_sha256": PRIMARY_SHA256,
        "primary_artifact_commit": PRIMARY_ARTIFACT_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
    },
    "method": {
        "elimination": "constraint-first R4(q)",
        "decisive_identity": "R4'(q)=4*pi*(m3-mu(q))/q^2",
        "primary_E4_used": False,
        "primary_artifact_read_after_construction": True,
    },
    "precision_records": [pack_record(row) for row in precision_records],
    "primary_comparison": {
        "root_differences": comparison,
        "state_differences": {
            name: text(value, 20) for name, value in state_differences.items()
        },
        "passed": comparison_ok,
    },
    "interpretation": {
        "label": (
            "DERIVED COMPUTATIONAL, FOUR-SLAB SCOPED / STRUCTURAL / "
            "ADVERSARIALLY CORROBORATED"
        ),
        "indefinite_history": "OPEN",
        "asymptotic_self_similarity": "OPEN",
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

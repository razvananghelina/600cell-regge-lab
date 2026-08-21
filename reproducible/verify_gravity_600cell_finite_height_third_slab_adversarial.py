#!/usr/bin/env python3
"""Equal-mu replication of the third-slab future-extendibility census."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
DIRECT_COMPOSITION_INPUT = (
    HERE / "gravity_600cell_finite_height_composition_adversarial.json"
)
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_third_slab.json"
OUTPUT = HERE / "gravity_600cell_finite_height_third_slab_adversarial.json"

DIRECT_COMPOSITION_SHA256 = (
    "d50e87f736e51585596aa1d7778238febaf7422840d668499878d8bd917f99e9"
)
PRIMARY_SHA256 = (
    "6b0e92d031aa891fdc3e1b2045c35bd135a955bb1374c92f015dcd5727d3d8fc"
)
ADVERSARIAL_PROTOCOL_COMMIT = "2f1c6f0"
PRIMARY_ARTIFACT_COMMIT = "a0abf0a"
ADVERSARIAL_CORRECTION_COMMIT = "e1898f6"

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


def text(value, digits=70):
    return mp.nstr(value, digits)


# The primary third-slab artifact remains unread until the independent dual
# census is complete.
direct_composition = json.loads(DIRECT_COMPOSITION_INPUT.read_text())
direct_provenance_ok = bool(
    digest(DIRECT_COMPOSITION_INPUT) == DIRECT_COMPOSITION_SHA256
    and direct_composition["outcome"]
    == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED"
)
check(
    "the direct-action two-branch input and dual protocol are frozen",
    direct_provenance_ok,
    f"protocol={ADVERSARIAL_PROTOCOL_COMMIT}",
)


# Derive the dual decisive identity without constructing the primary E(q).
q_symbol, u_symbol, up_symbol, pp_symbol, m_symbol = sp.symbols(
    "q u up pp m", real=True
)
pi_symbol = sp.pi
r_prime = (
    pp_symbol
    + 4
    * pi_symbol
    * (
        up_symbol * q_symbol - (u_symbol - m_symbol)
    )
    / q_symbol**2
)
r_prime_reduced = sp.factor(
    r_prime.subs(pp_symbol, -4 * pi_symbol * up_symbol / q_symbol)
)
dual_identity_ok = sp.simplify(
    r_prime_reduced
    - 4 * pi_symbol * (m_symbol - u_symbol) / q_symbol**2
) == 0
check(
    "constraint-first elimination gives R'(q)=4*pi*(m-mu(q))/q^2",
    dual_identity_ok,
)


# Complete action reconstructed independently of both third-slab artifacts.
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
        maxsteps=100,
    )


def raw_bisect(function, left, right, tolerance):
    f_left = function(left)
    f_right = function(right)
    if f_left == 0:
        return left
    if f_right == 0:
        return right
    if f_left * f_right >= 0:
        raise RuntimeError("dual bisection lacks a sign change")
    for _ in range(1000):
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
    raise RuntimeError("dual bisection did not converge")


def outer_bracket(function, boundary, direction):
    fixed = function(boundary)
    magnitude = max(mp.mpf(10), 2 * abs(boundary))
    for _ in range(400):
        point = direction * magnitude
        value = function(point)
        if value * fixed < 0:
            return (point, boundary) if direction < 0 else (boundary, point)
        magnitude *= 2
    raise RuntimeError("dual outer root did not bracket")


def solve_equal_mu(mass, v_star, tolerance, margin):
    mu_zero = mu(0)
    mu_star = mu(v_star)
    exceptional = bool(
        abs(mass) < margin
        or abs(mass - mu_zero) < margin
        or abs(mass - mu_star) < margin
    )
    if exceptional or mass <= 0 or mass > mu_star:
        return [], False, {
            "reason": "exceptional or inadmissible equal-mu level"
        }

    f = lambda value: mu(value) - mass
    positive = []
    if mass > mu_zero:
        positive.append(raw_bisect(f, mp.mpf(0), v_star, tolerance))
    bracket = outer_bracket(f, v_star, 1)
    positive.append(raw_bisect(f, *bracket, tolerance))
    stationary = sorted([-value for value in positive] + positive)
    return stationary, True, {
        "mu_zero_margin": text(abs(mass - mu_zero), 30),
        "mu_star_margin": text(abs(mass - mu_star), 30),
    }


def finite_point_for_limit(
    function, direction, target_sign, boundary, near_zero=False
):
    if near_zero:
        magnitude = min(mp.mpf("0.1"), abs(boundary) / 2)
        for _ in range(500):
            point = direction * magnitude
            value = function(point)
            if (value > 0) == (target_sign > 0):
                return point
            magnitude /= 2
    else:
        magnitude = max(mp.mpf(10), 2 * abs(boundary))
        for _ in range(500):
            point = direction * magnitude
            value = function(point)
            if (value > 0) == (target_sign > 0):
                return point
            magnitude *= 2
    raise RuntimeError("dual one-sided limit was not reached")


def dual_census(mass, incoming_p, v_star, tolerance, margin):
    stationary, stationary_ok, diagnostics = solve_equal_mu(
        mass, v_star, tolerance, margin
    )
    if not stationary_ok:
        return [], False, diagnostics

    function = lambda value: (
        momentum(value)
        - incoming_p
        + 4 * mp.pi * (mu(value) - mass) / value
    )
    mu_zero_gap = mu(0) - mass
    p_infinity = 60 * mp.pi - 300 * mp.sqrt(3) * mp.log(2)
    limit_values = {
        "negative_infinity": -p_infinity - incoming_p,
        "negative_zero_sign_source": -mu_zero_gap,
        "positive_zero_sign_source": mu_zero_gap,
        "positive_infinity": p_infinity - incoming_p,
    }
    if any(abs(value) < margin for value in limit_values.values()):
        return [], False, {
            **diagnostics,
            "reason": "dual one-sided limit is exceptional",
        }
    stationary_values = {value: function(value) for value in stationary}
    if any(abs(value) < margin for value in stationary_values.values()):
        return [], False, {
            **diagnostics,
            "reason": "dual stationary point is itself a root",
        }

    negative_stationary = [value for value in stationary if value < 0]
    positive_stationary = [value for value in stationary if value > 0]
    roots = []
    intervals = []

    def add_side(boundaries, limit_signs, side):
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
                        else text(left, 30)
                    ),
                    "right": (
                        ("0-" if side == "negative" else "+infinity")
                        if right is None
                        else text(right, 30)
                    ),
                    "left_sign": left_sign,
                    "right_sign": right_sign,
                    "has_root": has_root,
                }
            )
            if not has_root:
                continue
            if left is None:
                if side == "negative":
                    finite_left = finite_point_for_limit(
                        function,
                        -1,
                        left_sign,
                        right,
                        near_zero=False,
                    )
                else:
                    finite_left = finite_point_for_limit(
                        function,
                        1,
                        left_sign,
                        right,
                        near_zero=True,
                    )
                bracket_left = finite_left
            else:
                bracket_left = left
            if right is None:
                if side == "negative":
                    finite_right = finite_point_for_limit(
                        function,
                        -1,
                        right_sign,
                        left,
                        near_zero=True,
                    )
                else:
                    finite_right = finite_point_for_limit(
                        function,
                        1,
                        right_sign,
                        left,
                        near_zero=False,
                    )
                bracket_right = finite_right
            else:
                bracket_right = right
            roots.append(
                raw_bisect(function, bracket_left, bracket_right, tolerance)
            )

    negative_limits = (
        1 if limit_values["negative_infinity"] > 0 else -1,
        1 if limit_values["negative_zero_sign_source"] > 0 else -1,
    )
    positive_limits = (
        1 if limit_values["positive_zero_sign_source"] > 0 else -1,
        1 if limit_values["positive_infinity"] > 0 else -1,
    )
    add_side([None, *negative_stationary, None], negative_limits, "negative")
    add_side([None, *positive_stationary, None], positive_limits, "positive")
    roots.sort()

    q_zero_constraint = 8 * mp.pi * (mu(0) - mass)
    q_zero_excluded = abs(q_zero_constraint) > margin
    return roots, q_zero_excluded, {
        **diagnostics,
        "stationary_points": [text(value) for value in stationary],
        "stationary_values": [
            text(stationary_values[value], 30) for value in stationary
        ],
        "one_sided_limit_sources": {
            name: text(value, 30) for name, value in limit_values.items()
        },
        "q_zero_constraint": text(q_zero_constraint, 30),
        "q_zero_excluded": q_zero_excluded,
        "intervals": intervals,
    }


TARGET_PRECISIONS = (100, 160)
precision_records = []
construction_ok = True
census_ok = True
hostile_ok = True

for target_precision in TARGET_PRECISIONS:
    mp.mp.dps = target_precision + 60
    tolerance = mp.mpf(10) ** (-(target_precision + 20))
    margin = mp.mpf(10) ** (-(target_precision - 20))
    x_star = mp.findroot(k_function, (mp.mpf(5), mp.mpf(6)))
    v_star = mp.sqrt(x_star)

    v = mp.mpf(3) / 2
    mass0 = mu(v)
    pi0 = momentum(v)
    h1, q1 = solve_direct(
        (mp.mpf(1) / 5, mp.mpf(10)), mass0, pi0, tolerance
    )
    r1 = 1 + h1 * q1
    p_post1 = post_numeric(1, r1, h1**2, mass0)
    m1 = mass0 / r1
    pi1 = p_post1 / r1**2

    branches = {}
    for branch, seed in {
        "A": (mp.mpf(7), mp.mpf(1) / 50),
        "B": (mp.mpf(1) / 14, mp.mpf(31)),
    }.items():
        h2, q2 = solve_direct(seed, m1, pi1, tolerance)
        r2 = 1 + h2 * q2
        p_post2 = post_numeric(1, r2, h2**2, m1)
        m2 = m1 / r2
        pi2 = p_post2 / r2**2
        second_residuals = direct_residuals(m1, pi1, h2, q2)
        branch_construction_ok = bool(
            h2 > 0
            and r2 > 0
            and max(abs(value) for value in second_residuals)
            < mp.mpf("1e-90")
        )
        construction_ok &= branch_construction_ok

        roots, certified, diagnostics = dual_census(
            m2, pi2, v_star, tolerance, margin
        )
        root_rows = []
        physical_count = 0
        branch_direct_ok = certified
        for q3 in roots:
            h3 = 2 * (m2 - mu(q3)) / (q3 * mu(q3))
            r3 = 1 + h3 * q3
            reduced = reduced_residuals(m2, pi2, h3, q3)
            physical = bool(h3 > 0 and r3 > 0)
            algebraic_ok = max(abs(value) for value in reduced) < mp.mpf(
                "1e-90"
            )
            direct = (mp.nan, mp.nan)
            direct_ok = True
            if physical:
                direct = direct_residuals(m2, pi2, h3, q3)
                direct_ok = max(abs(value) for value in direct) < mp.mpf(
                    "1e-90"
                )
            branch_direct_ok &= algebraic_ok and direct_ok
            physical_count += int(physical)
            root_rows.append(
                {
                    "q3": q3,
                    "h3": h3,
                    "ratio": r3,
                    "constraint_residual": reduced[0],
                    "momentum_residual": reduced[1],
                    "direct_constraint_residual": direct[0],
                    "direct_momentum_residual": direct[1],
                    "physical": physical,
                    "algebraic_pass": algebraic_ok,
                    "direct_pass": direct_ok,
                }
            )
        census_ok &= branch_direct_ok

        wrong_scale = p_post2 / r2
        wrong_sign = -pi2
        reset_mass = mu(q2)
        branch_hostile_ok = bool(
            abs(wrong_scale - pi2) > mp.mpf("1e-20")
            and abs(wrong_sign - pi2) > mp.mpf("1e-20")
            and abs(reset_mass - m2) > mp.mpf("1e-20")
        )
        hostile_ok &= branch_hostile_ok
        branches[branch] = {
            "h2": h2,
            "q2": q2,
            "r2": r2,
            "m2": m2,
            "pi2": pi2,
            "root_count_certified": certified,
            "all_real_root_count": len(roots),
            "physical_root_count": physical_count,
            "diagnostics": diagnostics,
            "roots": root_rows,
            "hostile": {
                "wrong_scale_gap": wrong_scale - pi2,
                "wrong_sign_gap": wrong_sign - pi2,
                "mass_reset_gap": reset_mass - m2,
                "passed": branch_hostile_ok,
            },
        }
    precision_records.append(
        {
            "target_precision": target_precision,
            "v_star": v_star,
            "h1": h1,
            "q1": q1,
            "m1": m1,
            "pi1": pi1,
            "branches": branches,
        }
    )

check(
    "both precision runs reconstruct the first and second slabs from the full action",
    construction_ok,
    "precisions=100,160",
)
check(
    "the equal-mu monotone proof certifies every direct third root",
    census_ok,
    "q=0 and all four one-sided limits included",
)
check(
    "wrong scaling, sign and mass reset fail for both branches at both precisions",
    hostile_ok,
)


precision_nesting_ok = True
first_run, second_run = precision_records
for branch in ("A", "B"):
    first_branch = first_run["branches"][branch]
    second_branch = second_run["branches"][branch]
    precision_nesting_ok &= bool(
        first_branch["all_real_root_count"]
        == second_branch["all_real_root_count"]
        and first_branch["physical_root_count"]
        == second_branch["physical_root_count"]
        and len(first_branch["roots"]) == len(second_branch["roots"])
    )
    for first_root, second_root in zip(
        first_branch["roots"], second_branch["roots"]
    ):
        precision_nesting_ok &= bool(
            abs(first_root["q3"] - second_root["q3"]) < mp.mpf("1e-60")
            and abs(first_root["h3"] - second_root["h3"])
            < mp.mpf("1e-60")
            and first_root["physical"] == second_root["physical"]
        )
check(
    "the dual root censuses and physical labels nest beyond 60 digits",
    precision_nesting_ok,
)


independent_counts = {
    branch: second_run["branches"][branch]["physical_root_count"]
    for branch in ("A", "B")
}
independent_outcome_ok = independent_counts == {"A": 0, "B": 1}
check(
    "the independent dual census leaves only one B-branch continuation",
    independent_outcome_ok,
    f"physical_counts={independent_counts}",
)


# Only now read the primary third-slab artifact and compare every root.
primary = json.loads(PRIMARY_INPUT.read_text())
primary_provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and primary["outcome"] == "ONE_SECOND_BRANCH_EXTENDS_UNIQUELY"
    and primary["passed"] == primary["tests"] == 8
)
comparison_ok = primary_provenance_ok
comparison = {}
for branch in ("A", "B"):
    independent_branch = second_run["branches"][branch]
    primary_branch = primary["third_slab_census"][branch]
    row_ok = bool(
        independent_branch["all_real_root_count"]
        == primary_branch["all_real_root_count"]
        and independent_branch["physical_root_count"]
        == primary_branch["physical_root_count"]
        and len(independent_branch["roots"]) == len(primary_branch["roots"])
    )
    root_differences = []
    for independent_root, primary_root in zip(
        independent_branch["roots"], primary_branch["roots"]
    ):
        differences = {
            "q3": abs(
                independent_root["q3"] - mp.mpf(primary_root["q3"])
            ),
            "h3": abs(
                independent_root["h3"] - mp.mpf(primary_root["h3"])
            ),
            "ratio": abs(
                independent_root["ratio"]
                - mp.mpf(primary_root["scale_ratio"])
            ),
        }
        this_ok = bool(
            max(differences.values()) < mp.mpf("1e-55")
            and independent_root["physical"] == primary_root["physical"]
        )
        row_ok &= this_ok
        root_differences.append(
            {
                **{name: text(value, 20) for name, value in differences.items()},
                "passed": this_ok,
            }
        )
    comparison_ok &= row_ok
    comparison[branch] = {
        "root_differences": root_differences,
        "passed": row_ok,
    }
check(
    "the equal-mu census agrees with the equal-p primary only after construction",
    comparison_ok,
)


outcome = (
    "ONE_SECOND_BRANCH_EXTENDS_UNIQUELY_ADVERSARIALLY_CORROBORATED"
    if direct_provenance_ok
    and dual_identity_ok
    and construction_ok
    and census_ok
    and hostile_ok
    and precision_nesting_ok
    and independent_outcome_ok
    and comparison_ok
    else "THIRD_SLAB_EXTENDIBILITY_ADVERSARIAL_OPEN"
)
check(
    "the adversarial hierarchy corroborates the three-slab branch distinction",
    outcome
    == "ONE_SECOND_BRANCH_EXTENDS_UNIQUELY_ADVERSARIALLY_CORROBORATED",
)


def pack_number(value):
    return value if isinstance(value, bool) else text(value)


def pack_precision(row):
    return {
        "target_precision": row["target_precision"],
        "v_star": text(row["v_star"]),
        "h1": text(row["h1"]),
        "q1": text(row["q1"]),
        "m1": text(row["m1"]),
        "pi1": text(row["pi1"]),
        "branches": {
            branch: {
                "h2": text(branch_row["h2"]),
                "q2": text(branch_row["q2"]),
                "r2": text(branch_row["r2"]),
                "m2": text(branch_row["m2"]),
                "pi2": text(branch_row["pi2"]),
                "root_count_certified": branch_row["root_count_certified"],
                "all_real_root_count": branch_row["all_real_root_count"],
                "physical_root_count": branch_row["physical_root_count"],
                "diagnostics": branch_row["diagnostics"],
                "roots": [
                    {
                        name: pack_number(value)
                        for name, value in root.items()
                    }
                    for root in branch_row["roots"]
                ],
                "hostile": {
                    name: pack_number(value)
                    for name, value in branch_row["hostile"].items()
                },
            }
            for branch, branch_row in row["branches"].items()
        },
    }


artifact = {
    "provenance": {
        "direct_composition_sha256": DIRECT_COMPOSITION_SHA256,
        "primary_third_slab_sha256": PRIMARY_SHA256,
        "primary_artifact_commit": PRIMARY_ARTIFACT_COMMIT,
        "adversarial_protocol_commit": ADVERSARIAL_PROTOCOL_COMMIT,
        "adversarial_correction_commit": ADVERSARIAL_CORRECTION_COMMIT,
    },
    "method": {
        "elimination": "constraint-first R(q)",
        "decisive_identity": "R'(q)=4*pi*(m-mu(q))/q^2",
        "primary_E_used": False,
        "primary_artifact_read_after_construction": True,
    },
    "precision_records": [pack_precision(row) for row in precision_records],
    "physical_counts": independent_counts,
    "primary_comparison": comparison,
    "interpretation": {
        "label": (
            "DERIVED COMPUTATIONAL, THREE-SLAB SCOPED / STRUCTURAL / "
            "ADVERSARIALLY CORROBORATED"
        ),
        "future_extendibility_is_physical_selector": False,
        "deterministic_tick": False,
        "indefinite_extension": "OPEN",
    },
    "tests": tests,
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
raise SystemExit(0 if passed == tests else 1)

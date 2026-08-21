#!/usr/bin/env python3
"""Complete third-slab extendibility census for both second branches."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
COMPOSITION_INPUT = HERE / "gravity_600cell_finite_height_composition.json"
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
SELECTOR_INPUT = HERE / "gravity_600cell_finite_height_selector_audit.json"
SELECTOR_ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_finite_height_selector_audit_adversarial.json"
)
OUTPUT = HERE / "gravity_600cell_finite_height_third_slab.json"

COMPOSITION_SHA256 = (
    "d4e36141863bd2ae515b96eeeff4f50eb087016cca8cfb6f4b1e3355d6fba447"
)
CLASSIFICATION_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
SELECTOR_SHA256 = (
    "956cd655b8b3a5106029fb852df74b85bb59f922a4984542bc2e089f54799676"
)
SELECTOR_ADVERSARIAL_SHA256 = (
    "1fe11f006cd928dc5418c3171154a4b0e26db79225e69845f2edd9566f820f0e"
)
PROTOCOL_COMMIT = "0f31fe8"
BISECTION_CORRECTION_COMMIT = "9c40419"

mp.mp.dps = 150
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


composition = json.loads(COMPOSITION_INPUT.read_text())
classification = json.loads(CLASSIFICATION_INPUT.read_text())
selector = json.loads(SELECTOR_INPUT.read_text())
selector_adversarial = json.loads(SELECTOR_ADVERSARIAL_INPUT.read_text())
provenance_ok = bool(
    digest(COMPOSITION_INPUT) == COMPOSITION_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and digest(SELECTOR_INPUT) == SELECTOR_SHA256
    and digest(SELECTOR_ADVERSARIAL_INPUT) == SELECTOR_ADVERSARIAL_SHA256
    and composition["outcome"] == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE"
    and selector["outcome"]
    == "STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH"
    and selector_adversarial["outcome"]
    == "STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH_"
    "ADVERSARIALLY_CORROBORATED"
)
check(
    "the two branches, selector negatives and third-slab protocol are frozen",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


# Reconstruct the complete action and both canonical momenta before using a
# closed outgoing-state recurrence.
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

swapped_action = ACTION.xreplace({L_MINUS: L_PLUS, L_PLUS: L_MINUS})
swapped_pre = P_PRE_EXACT.xreplace({L_MINUS: L_PLUS, L_PLUS: L_MINUS})
action_reversal_ok = bool(
    sp.simplify(swapped_action - ACTION) == 0
    and sp.simplify(swapped_pre + P_POST_EXACT) == 0
)
check(
    "the complete action is endpoint symmetric and reverses pre/post momentum sign",
    action_reversal_ok,
)


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


def reduced_residuals(mass, incoming_p, height, slope):
    return (
        8 * mp.pi * (mu(slope) - mass)
        + 4 * mp.pi * height * slope * mu(slope),
        momentum(slope)
        - incoming_p
        - 2 * mp.pi * height * mu(slope),
    )


def direct_residuals(mass, incoming_p, height, slope):
    endpoint = 1 + height * slope
    return (
        2 * f_numeric(1, endpoint, height**2, mass) / height,
        pre_numeric(1, endpoint, height**2, mass) - incoming_p,
    )


def solve_direct(seed, mass, incoming_p):
    return mp.findroot(
        lambda height, slope: direct_residuals(
            mass, incoming_p, height, slope
        ),
        seed,
        tol=mp.mpf("1e-130"),
        maxsteps=100,
    )


# Reconstruct the first and both second slabs at working precision; the input
# artifact supplies provenance and branch labels, not truncated solve values.
v = mp.mpf(3) / 2
mass0 = mu(v)
pi0 = momentum(v)
h1, q1 = solve_direct((mp.mpf(1) / 5, mp.mpf(10)), mass0, pi0)
r1 = 1 + h1 * q1
p_post1 = post_numeric(1, r1, h1**2, mass0)
m1 = mass0 / r1
pi1 = p_post1 / r1**2

second_seeds = {
    "A": (mp.mpf(7), mp.mpf(1) / 50),
    "B": (mp.mpf(1) / 14, mp.mpf(31)),
}
second_roots = {}
reconstruction_ok = True
for branch, seed in second_seeds.items():
    h2, q2 = solve_direct(seed, m1, pi1)
    r2 = 1 + h2 * q2
    residuals = direct_residuals(m1, pi1, h2, q2)
    reconstruction_ok &= bool(
        h2 > 0
        and r2 > 0
        and max(abs(value) for value in residuals) < mp.mpf("1e-120")
    )
    second_roots[branch] = {"h": h2, "q": q2, "ratio": r2}

frozen_roots = [
    root
    for root in next(
        row for row in composition["composition"] if row["v"] == "1.5"
    )["roots"]
    if root["physical"]
]
for branch, frozen in zip(("A", "B"), frozen_roots):
    reconstruction_ok &= bool(
        abs(second_roots[branch]["h"] - mp.mpf(frozen["h2"]))
        < mp.mpf("1e-55")
        and abs(second_roots[branch]["q"] - mp.mpf(frozen["q2"]))
        < mp.mpf("1e-55")
    )
check(
    "the first slab and both frozen second branches reconstruct from the full action",
    reconstruction_ok,
)


# Derive and test the exact normalized outgoing recurrence.  Endpoint
# reversal gives p_pre(reverse)=-p_post/r^2, while even mu and odd p evaluate
# the reversed normalized slab at (-q,h/r,m/r).
recurrence_ok = True
outgoing_states = {}
for branch, root in second_roots.items():
    h2 = root["h"]
    q2 = root["q"]
    r2 = root["ratio"]
    p_post2 = post_numeric(1, r2, h2**2, m1)
    m2 = m1 / r2
    pi2_direct = p_post2 / r2**2
    pi2_closed = momentum(q2) + 2 * mp.pi * h2 * mu(q2) / r2
    recurrence_error = abs(pi2_direct - pi2_closed)
    recurrence_ok &= recurrence_error < mp.mpf("1e-120")
    outgoing_states[branch] = {
        "m": m2,
        "pi": pi2_direct,
        "p_post": p_post2,
        "recurrence_error": recurrence_error,
    }
check(
    "the full action confirms the closed outgoing-state recurrence on both branches",
    recurrence_ok,
)


v_star = mp.mpf(classification["thresholds"]["v_star"])
p_star = momentum(v_star)
p_infinity = 60 * mp.pi - 300 * mp.sqrt(3) * mp.log(2)
SIGN_TOL = mp.mpf("1e-90")
ROOT_TOL = mp.mpf("1e-115")


def sign(value):
    if value > SIGN_TOL:
        return 1
    if value < -SIGN_TOL:
        return -1
    return 0


def bisect_root(function, left, right):
    f_left = function(left)
    f_right = function(right)
    if f_left == 0:
        return left
    if f_right == 0:
        return right
    if f_left * f_right >= 0:
        raise RuntimeError("bisection interval has no certified sign change")
    for _ in range(700):
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
    raise RuntimeError("bisection did not converge")


def outer_bracket(function, boundary, direction):
    magnitude = max(mp.mpf(10), 2 * abs(boundary))
    fixed_value = function(boundary)
    for _ in range(300):
        tail_point = direction * magnitude
        tail_value = function(tail_point)
        if sign(tail_value) != 0 and sign(tail_value) != sign(fixed_value):
            return (
                (tail_point, boundary)
                if direction < 0
                else (boundary, tail_point)
            )
        magnitude *= 2
    raise RuntimeError("failed to bracket an outer monotone root")


def stationary_points(incoming_p):
    critical_levels = (
        p_star,
        p_infinity,
        mp.mpf(0),
        -p_infinity,
        -p_star,
    )
    exceptional = any(
        abs(incoming_p - level) < SIGN_TOL for level in critical_levels
    )
    if exceptional:
        return [], False, {"reason": "incoming momentum hits a critical level"}

    roots = []
    f = lambda value: momentum(value) - incoming_p
    intervals = [
        (None, -v_star, -p_infinity, -p_star),
        (-v_star, mp.mpf(0), -p_star, mp.mpf(0)),
        (mp.mpf(0), v_star, mp.mpf(0), p_star),
        (v_star, None, p_star, p_infinity),
    ]
    for left, right, p_left, p_right in intervals:
        low = min(p_left, p_right)
        high = max(p_left, p_right)
        if not (incoming_p > low and incoming_p < high):
            continue
        if left is None:
            bracket = outer_bracket(f, right, -1)
        elif right is None:
            bracket = outer_bracket(f, left, 1)
        else:
            bracket = (left, right)
        roots.append(bisect_root(f, *bracket))
    roots.sort()
    return roots, True, {
        "critical_level_margins": [
            text(abs(incoming_p - level), 30) for level in critical_levels
        ]
    }


def tail_point_for_sign(function, direction, expected_sign, boundary):
    magnitude = max(mp.mpf(10), 2 * abs(boundary))
    for _ in range(300):
        point = direction * magnitude
        if sign(function(point)) == expected_sign:
            return point
        magnitude *= 2
    raise RuntimeError("failed to certify elimination tail")


def enumerate_elimination_roots(mass, incoming_p):
    stationary, stationary_ok, diagnostics = stationary_points(incoming_p)
    if not stationary_ok:
        return [], False, diagnostics

    function = lambda value: (
        4 * mp.pi * (mu(value) - mass)
        + value * (momentum(value) - incoming_p)
    )
    stationary_values = [function(value) for value in stationary]
    if any(abs(value) < SIGN_TOL for value in stationary_values):
        return [], False, {
            **diagnostics,
            "reason": "stationary elimination root requires exceptional treatment",
        }

    right_coefficient = p_infinity - incoming_p
    left_coefficient = -p_infinity - incoming_p
    right_tail_sign = sign(right_coefficient)
    left_tail_sign = -sign(left_coefficient)
    if right_tail_sign == 0 or left_tail_sign == 0:
        return [], False, {
            **diagnostics,
            "reason": "linear elimination tail cancels",
        }

    boundaries = [None, *stationary, None]
    roots = []
    interval_certificates = []
    for index in range(len(boundaries) - 1):
        left = boundaries[index]
        right = boundaries[index + 1]
        left_sign = (
            left_tail_sign if left is None else sign(function(left))
        )
        right_sign = (
            right_tail_sign if right is None else sign(function(right))
        )
        has_root = left_sign != right_sign
        interval_certificates.append(
            {
                "left": "-infinity" if left is None else text(left, 30),
                "right": "+infinity" if right is None else text(right, 30),
                "left_sign": left_sign,
                "right_sign": right_sign,
                "has_root": has_root,
            }
        )
        if not has_root:
            continue
        if left is None:
            finite_left = tail_point_for_sign(
                function, -1, left_tail_sign, right
            )
            bracket = (finite_left, right)
        elif right is None:
            finite_right = tail_point_for_sign(
                function, 1, right_tail_sign, left
            )
            bracket = (left, finite_right)
        else:
            bracket = (left, right)
        roots.append(bisect_root(function, *bracket))

    q_zero_value = function(mp.mpf(0))
    roots.sort()
    deduplicated = []
    for root in roots:
        if not deduplicated or abs(root - deduplicated[-1]) > mp.mpf("1e-80"):
            deduplicated.append(root)
    return deduplicated, True, {
        **diagnostics,
        "stationary_points": [text(value) for value in stationary],
        "stationary_values": [text(value, 30) for value in stationary_values],
        "left_tail_coefficient": text(left_coefficient, 30),
        "right_tail_coefficient": text(right_coefficient, 30),
        "q_zero_value": text(q_zero_value, 30),
        "q_zero_is_root": abs(q_zero_value) < SIGN_TOL,
        "intervals": interval_certificates,
    }


branch_census = {}
census_ok = True
direct_checks_ok = True
for branch, state in outgoing_states.items():
    roots, certified, diagnostics = enumerate_elimination_roots(
        state["m"], state["pi"]
    )
    census_ok &= certified
    root_rows = []
    physical_count = 0
    for q3 in roots:
        h3 = (
            momentum(q3) - state["pi"]
        ) / (2 * mp.pi * mu(q3))
        r3 = 1 + h3 * q3
        determinant = 8 * mp.pi**2 * h3 * mu(q3) ** 2
        reduced = reduced_residuals(state["m"], state["pi"], h3, q3)
        physical = bool(h3 > 0 and r3 > 0 and determinant > 0)
        direct = (mp.nan, mp.nan)
        junction = mp.nan
        direct_pass = True
        if physical:
            direct = direct_residuals(state["m"], state["pi"], h3, q3)
            pi3_pre = pre_numeric(1, r3, h3**2, state["m"])
            junction = state["pi"] - pi3_pre
            direct_pass = bool(
                max(abs(value) for value in direct) < mp.mpf("1e-90")
                and abs(junction) < mp.mpf("1e-90")
            )
            direct_checks_ok &= direct_pass
            physical_count += 1
        root_rows.append(
            {
                "q3": text(q3),
                "h3": text(h3),
                "scale_ratio": text(r3),
                "determinant": text(determinant),
                "constraint_residual": text(reduced[0], 30),
                "momentum_residual": text(reduced[1], 30),
                "direct_constraint_residual": (
                    "not_evaluated" if not physical else text(direct[0], 30)
                ),
                "direct_momentum_residual": (
                    "not_evaluated" if not physical else text(direct[1], 30)
                ),
                "junction_residual": (
                    "not_evaluated" if not physical else text(junction, 30)
                ),
                "physical": physical,
                "direct_pass": direct_pass,
            }
        )
    branch_census[branch] = {
        "outgoing_state": {
            "m2": text(state["m"]),
            "pi2": text(state["pi"]),
            "recurrence_error": text(state["recurrence_error"], 30),
        },
        "all_real_root_count": len(roots),
        "physical_root_count": physical_count,
        "root_count_certified": certified,
        "diagnostics": diagnostics,
        "roots": root_rows,
    }

check(
    "both third-slab elimination equations have complete monotone all-real censuses",
    census_ok,
    "finite plotting boxes are not used",
)
check(
    "every physical third slab passes the complete action and shared-slice junction",
    direct_checks_ok,
)


hostile_controls = {}
hostile_ok = True
for branch, root in second_roots.items():
    r2 = root["ratio"]
    p_post2 = outgoing_states[branch]["p_post"]
    correct_pi = outgoing_states[branch]["pi"]
    wrong_scale = p_post2 / r2
    wrong_sign = -correct_pi
    reset_mass = mu(root["q"])
    row_ok = bool(
        abs(wrong_scale - correct_pi) > mp.mpf("1e-20")
        and abs(wrong_sign - correct_pi) > mp.mpf("1e-20")
        and abs(reset_mass - outgoing_states[branch]["m"])
        > mp.mpf("1e-20")
    )
    hostile_ok &= row_ok
    hostile_controls[branch] = {
        "wrong_scale_gap": text(wrong_scale - correct_pi, 30),
        "wrong_sign_gap": text(wrong_sign - correct_pi, 30),
        "mass_reset_gap": text(reset_mass - outgoing_states[branch]["m"], 30),
        "passed": row_ok,
    }
check(
    "wrong scaling, sign and mass-reset conventions all change both outgoing states",
    hostile_ok,
)


counts = {
    branch: branch_census[branch]["physical_root_count"]
    for branch in ("A", "B")
}
if census_ok and direct_checks_ok and recurrence_ok and hostile_ok:
    extending = [branch for branch, count in counts.items() if count > 0]
    if len(extending) == 2:
        outcome = "BOTH_SECOND_BRANCHES_EXTEND"
    elif len(extending) == 0:
        outcome = "NEITHER_SECOND_BRANCH_EXTENDS"
    elif counts[extending[0]] == 1:
        outcome = "ONE_SECOND_BRANCH_EXTENDS_UNIQUELY"
    else:
        outcome = "ONE_SECOND_BRANCH_EXTENDS_NONUNIQUELY"
else:
    outcome = "THIRD_SLAB_EXTENDIBILITY_OPEN"

check(
    "the preregistered hierarchy classifies third-slab future extendibility",
    outcome != "THIRD_SLAB_EXTENDIBILITY_OPEN",
    f"physical_counts={counts}",
)


artifact = {
    "provenance": {
        "composition_sha256": COMPOSITION_SHA256,
        "classification_sha256": CLASSIFICATION_SHA256,
        "selector_sha256": SELECTOR_SHA256,
        "selector_adversarial_sha256": SELECTOR_ADVERSARIAL_SHA256,
        "protocol_commit": PROTOCOL_COMMIT,
        "bisection_correction_commit": BISECTION_CORRECTION_COMMIT,
    },
    "exact_certificates": {
        "action_endpoint_symmetry": action_reversal_ok,
        "outgoing_mass": "m_next=m/r",
        "outgoing_momentum": "pi_next=p(q)+2*pi*h*mu(q)/r",
        "elimination_derivative": "E'(q)=p(q)-pi",
        "local_determinant": "8*pi^2*h*mu(q)^2",
    },
    "second_branches": {
        branch: {
            "h2": text(root["h"]),
            "q2": text(root["q"]),
            "scale_ratio": text(root["ratio"]),
        }
        for branch, root in second_roots.items()
    },
    "third_slab_census": branch_census,
    "physical_counts": counts,
    "hostile_controls": hostile_controls,
    "interpretation": {
        "label": "DERIVED COMPUTATIONAL, THREE-SLAB SCOPED / STRUCTURAL",
        "future_extendibility_is_physical_selector": False,
        "deterministic_tick": False,
        "global_history": "OPEN",
    },
    "tests": tests,
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
raise SystemExit(0 if passed == tests else 1)

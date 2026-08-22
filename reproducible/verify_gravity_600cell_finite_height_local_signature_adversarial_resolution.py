#!/usr/bin/env python3
"""Resolve the local-signature adversarial OPEN using exact p monotonicity."""

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp
from flint import arb, ctx


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = (
    ROOT
    / "docs"
    / "gravity"
    / "gravity_600cell_finite_height_local_signature_adversarial_resolution_protocol.md"
)
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
INVARIANT_INPUT = (
    HERE
    / "gravity_600cell_finite_height_invariant_region_adversarial_resolution.json"
)
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_local_signature.json"
FIRST_OPEN_INPUT = (
    HERE / "gravity_600cell_finite_height_local_signature_adversarial.json"
)
OUTPUT = (
    HERE
    / "gravity_600cell_finite_height_local_signature_adversarial_resolution.json"
)

PROTOCOL_COMMIT = "55c2ef0"
PROTOCOL_SHA256 = (
    "1bc25245ffae3a009242ed76816767128e928551030a016f66abdb8c7a6a73bc"
)
CLASSIFICATION_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
INVARIANT_SHA256 = (
    "813e05bd66b47cc3ae1cd35d0a2eddb9c645a850d84abeaad37d15b14a6a380f"
)
PRIMARY_SHA256 = (
    "9f524cc22df8cfb5083f372481b3efd19868252b85551d56378327eea7a6d613"
)
FIRST_OPEN_SHA256 = (
    "139dcee2e9ee021c131aae1090433fe16bd70c9f2b10ec52d32b0c5ebd7748a7"
)

ARB_DPS = 240
MAX_BISECTION_STEPS = 420
MIN_AMBIGUITY_STEP = 240

ctx.dps = ARB_DPS
IDENTITY_WIDTH_LIMIT = arb("1e-110")
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


def rational_ball(value):
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def interval_ball(left, right):
    left = Fraction(left)
    right = Fraction(right)
    midpoint = (left + right) / 2
    radius = (right - left) / 2
    return rational_ball(midpoint) + arb(
        0, f"{radius.numerator}/{radius.denominator}"
    )


def strict_sign(value):
    if value.lower() > 0:
        return 1
    if value.upper() < 0:
        return -1
    return 0


def arb_record(value):
    return {
        "pretty": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "contains_zero": bool(value.contains(0)),
        "sign": strict_sign(value),
    }


def series_asinh(value):
    return (value + (value * value + 1).sqrt()).log()


API = arb.pi()


def arb_epsilon(value):
    square = value * value
    return 2 * API - 5 * (
        (square + 2) / (2 * (square + 3))
    ).acos()


def arb_mu(value):
    return 180 * arb_epsilon(value) / (
        API * (value * value + 4).sqrt()
    )


def arb_p(value):
    square = value * value
    return (
        180 * value * arb_epsilon(value) / (square + 4).sqrt()
        - 600
        * arb(3).sqrt()
        * series_asinh(value / (8 * (square + 3)).sqrt())
    )


def arb_E(mass, pi_value, q_value):
    return 4 * API * (arb_mu(q_value) - mass) + q_value * (
        arb_p(q_value) - pi_value
    )


def bisect_from_signs(function, left, right):
    left = Fraction(left)
    right = Fraction(right)
    original_left = left
    original_right = right
    left_value = function(rational_ball(left))
    right_value = function(rational_ball(right))
    left_sign = strict_sign(left_value)
    right_sign = strict_sign(right_value)
    if left_sign * right_sign >= 0:
        return interval_ball(left, right), False, {
            "initial_left": str(original_left),
            "initial_right": str(original_right),
            "left_value": arb_record(left_value),
            "right_value": arb_record(right_value),
            "reason": "endpoint signs are not strictly opposite",
        }

    ambiguity_step = None
    steps = 0
    for step in range(MAX_BISECTION_STEPS):
        middle = (left + right) / 2
        middle_value = function(rational_ball(middle))
        middle_sign = strict_sign(middle_value)
        if middle_sign == 0:
            ambiguity_step = step
            break
        if left_sign * middle_sign < 0:
            right = middle
            right_sign = middle_sign
        else:
            left = middle
            left_sign = middle_sign
        steps = step + 1

    final_left_value = function(rational_ball(left))
    final_right_value = function(rational_ball(right))
    final_left_sign = strict_sign(final_left_value)
    final_right_sign = strict_sign(final_right_value)
    ambiguity_ok = bool(
        ambiguity_step is None or ambiguity_step >= MIN_AMBIGUITY_STEP
    )
    certified = bool(
        ambiguity_ok and final_left_sign * final_right_sign < 0
    )
    enclosure = interval_ball(left, right)
    return enclosure, certified, {
        "initial_left": str(original_left),
        "initial_right": str(original_right),
        "initial_left_value": arb_record(left_value),
        "initial_right_value": arb_record(right_value),
        "steps": steps,
        "ambiguity_step": ambiguity_step,
        "final_left": str(left),
        "final_right": str(right),
        "final_width_fraction": str(right - left),
        "final_left_value": arb_record(final_left_value),
        "final_right_value": arb_record(final_right_value),
        "certified": certified,
    }


classification = json.loads(CLASSIFICATION_INPUT.read_text())
invariant = json.loads(INVARIANT_INPUT.read_text())
primary = json.loads(PRIMARY_INPUT.read_text())
first_open = json.loads(FIRST_OPEN_INPUT.read_text())
provenance_ok = bool(
    digest(PROTOCOL) == PROTOCOL_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and digest(INVARIANT_INPUT) == INVARIANT_SHA256
    and digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and digest(FIRST_OPEN_INPUT) == FIRST_OPEN_SHA256
    and classification["outcome"]
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
    and invariant["outcome"]
    == "INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED"
    and primary["outcome"] == "LOCAL_SIGNATURE_PRIMARY_CERTIFIED"
    and first_open["outcome"] == "LOCAL_SIGNATURE_ADVERSARIAL_OPEN"
)
check(
    "the resolution protocol and all frozen inputs are exact",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


# Reconstruct the accepted p' factorization algebraically.
q_symbol = sp.symbols("q", real=True)
x_symbol = q_symbol**2
r_symbol = sp.sqrt(x_symbol + 4)
s_symbol = sp.sqrt(3 * x_symbol + 8)
epsilon_symbol = sp.symbols("epsilon", real=True)
epsilon_prime = 10 * q_symbol / (
    (x_symbol + 3) * r_symbol * s_symbol
)
p_prime_chain = (
    180
    * (
        epsilon_symbol / r_symbol
        + q_symbol * epsilon_prime / r_symbol
        - q_symbol**2 * epsilon_symbol / r_symbol**3
    )
    - 1800 / ((x_symbol + 3) * s_symbol)
)
K_symbol = 10 * r_symbol - (x_symbol + 3) * s_symbol * epsilon_symbol
p_prime_expected = -720 * K_symbol / (
    r_symbol**3 * (x_symbol + 3) * s_symbol
)
factorization_ok = sp.factor(p_prime_chain - p_prime_expected) == 0
check(
    "the exact K-factor formula for p prime is independently reconstructed",
    factorization_ok,
)


# Frozen polynomial controls before gravity.
poly_interval, poly_ok, poly_record = bisect_from_signs(
    lambda value: value * value - 2,
    Fraction(1),
    Fraction(2),
)
_, negative_ok, negative_record = bisect_from_signs(
    lambda value: value * value - 2,
    Fraction(2),
    Fraction(3),
)
polynomial_controls_ok = bool(poly_ok and not negative_ok)
check(
    "sign bisection accepts sqrt(2) on [1,2] and rejects [2,3]",
    polynomial_controls_ok,
)


STATIONARY_BRACKETS = {
    "root": ("diag", (5, 6)),
    "root/c0": ((1, 2), (17, 18)),
    "root/c0/c0": ((-1, 0),),
    "root/c0/c1": ((1, 2), (55, 56)),
    "root/c0/c1/c0": ((1, 2), (177, 178)),
}
ROOT_BRACKETS = {
    "root": ("diag", (9, 10)),
    "root/c0": ((0, 1), (9, 10), (31, 32)),
    "root/c0/c0": ((-3, -2), (0, 1)),
    "root/c0/c1": ((-1, 0), (31, 32), (99, 100)),
    "root/c0/c1/c0": ((-1, 0), (99, 100), (316, 317)),
}
EXPECTED_COUNTS = {
    "root": (2, 1),
    "root/c0": (3, 2),
    "root/c0/c0": (2, 0),
    "root/c0/c1": (3, 1),
    "root/c0/c1/c0": (3, 1),
}
EXPECTED_TERMINALS = {
    "root/c0": None,
    "root/c0/c0": "DEAD",
    "root/c0/c1": None,
    "root/c0/c1/c0": None,
    "root/c0/c1/c0/c0": "ENTERED_D",
}

V0 = arb(3) / 2
M0 = arb_mu(V0)
PI0 = arb_p(V0)
P_INFINITY = 60 * API - 300 * arb(3).sqrt() * arb(2).log()
V_STAR = arb(classification["thresholds"]["v_star"]) + arb(0, "1e-80")
NEG_V_STAR = -V_STAR
P_STAR = arb_p(V_STAR)


def positive_level_count(level):
    signs = {
        "minus_p_star": strict_sign(level - P_STAR),
        "minus_p_infinity": strict_sign(level - P_INFINITY),
        "level": strict_sign(level),
    }
    if 0 in signs.values():
        return None, signs
    if signs["minus_p_star"] < 0 or signs["level"] > 0:
        return 0, signs
    if signs["minus_p_infinity"] < 0:
        return 2, signs
    return 1, signs


def stationary_monotonic_region(left, right):
    left_ball = rational_ball(left)
    right_ball = rational_ball(right)
    if right_ball.upper() < NEG_V_STAR.lower():
        return "left_increasing"
    if (
        left_ball.lower() > NEG_V_STAR.upper()
        and right_ball.upper() < V_STAR.lower()
    ):
        return "central_decreasing"
    if left_ball.lower() > V_STAR.upper():
        return "right_increasing"
    return None


def p_gap_monotone_certificate(pi_value, left, right):
    evaluations = [
        ("left", arb_p(rational_ball(left)) - pi_value),
        ("right", arb_p(rational_ball(right)) - pi_value),
    ]
    left_ball = rational_ball(left)
    right_ball = rational_ball(right)
    if (
        left_ball.lower() < NEG_V_STAR.lower()
        and right_ball.upper() > NEG_V_STAR.upper()
    ):
        evaluations.append(("minus_v_star", arb_p(NEG_V_STAR) - pi_value))
    if (
        left_ball.lower() < V_STAR.lower()
        and right_ball.upper() > V_STAR.upper()
    ):
        evaluations.append(("v_star", arb_p(V_STAR) - pi_value))
    signs = [strict_sign(value) for _, value in evaluations]
    certified = bool(signs and signs[0] != 0 and all(s == signs[0] for s in signs))
    return signs[0] if certified else 0, certified, [
        {"point": name, "value": arb_record(value)}
        for name, value in evaluations
    ]


all_states_ok = True
all_gates_ok = True
all_identity_ok = True
state_records = {}
transition_records = []
d_entries = []


def certify_state(mass, pi_value, path):
    global all_states_ok, all_gates_ok, all_identity_ok

    pos_count, pos_signs = positive_level_count(pi_value)
    neg_count, neg_signs = positive_level_count(-pi_value)
    expected_stationary_count = (
        None
        if pos_count is None or neg_count is None
        else pos_count + neg_count
    )

    stationary_intervals = []
    stationary_records = []
    diagonal_index = None
    for index, bracket in enumerate(STATIONARY_BRACKETS[path]):
        if bracket == "diag":
            stationary_intervals.append(V0)
            stationary_records.append({"exact_diagonal": True})
            diagonal_index = index
            continue
        region = stationary_monotonic_region(
            Fraction(bracket[0]), Fraction(bracket[1])
        )
        enclosure, sign_ok, record = bisect_from_signs(
            lambda value: arb_p(value) - pi_value,
            Fraction(bracket[0]),
            Fraction(bracket[1]),
        )
        record["monotonic_region"] = region
        record["monotonicity_certified"] = region is not None
        stationary_intervals.append(enclosure)
        stationary_records.append(record)
        all_states_ok &= bool(sign_ok and region is not None)
    all_states_ok &= expected_stationary_count == len(stationary_intervals)
    all_states_ok &= all(
        stationary_intervals[index].upper()
        < stationary_intervals[index + 1].lower()
        for index in range(len(stationary_intervals) - 1)
    )

    stationary_values = []
    stationary_signs = []
    for index, enclosure in enumerate(stationary_intervals):
        value = arb_E(mass, pi_value, enclosure)
        stationary_values.append(value)
        if index == diagonal_index:
            stationary_signs.append(0)
        else:
            stationary_signs.append(strict_sign(value))
            all_states_ok &= strict_sign(value) != 0

    left_coefficient = -P_INFINITY - pi_value
    right_coefficient = P_INFINITY - pi_value
    left_tail_sign = -strict_sign(left_coefficient)
    right_tail_sign = strict_sign(right_coefficient)
    zero_value = 4 * API * (arb_mu(arb(0)) - mass)
    all_states_ok &= bool(
        left_tail_sign != 0
        and right_tail_sign != 0
        and strict_sign(zero_value) != 0
    )
    boundary_signs = [left_tail_sign, *stationary_signs, right_tail_sign]
    reversal_slots = [
        index
        for index in range(len(boundary_signs) - 1)
        if boundary_signs[index] * boundary_signs[index + 1] < 0
    ]
    predicted_count = len(reversal_slots) + int(diagonal_index is not None)

    root_intervals = []
    root_records = []
    diagonal_root_index = None
    for index, bracket in enumerate(ROOT_BRACKETS[path]):
        if bracket == "diag":
            root_intervals.append(V0)
            root_records.append({"exact_diagonal": True})
            diagonal_root_index = index
            continue
        derivative_sign, monotone_ok, monotone_record = (
            p_gap_monotone_certificate(
                pi_value, Fraction(bracket[0]), Fraction(bracket[1])
            )
        )
        enclosure, sign_ok, record = bisect_from_signs(
            lambda value: arb_E(mass, pi_value, value),
            Fraction(bracket[0]),
            Fraction(bracket[1]),
        )
        record["E_q_sign"] = derivative_sign
        record["monotone_segments"] = monotone_record
        record["monotonicity_certified"] = monotone_ok
        root_intervals.append(enclosure)
        root_records.append(record)
        all_states_ok &= bool(sign_ok and monotone_ok)
    all_states_ok &= predicted_count == len(root_intervals)
    all_states_ok &= all(
        root_intervals[index].upper() < root_intervals[index + 1].lower()
        for index in range(len(root_intervals) - 1)
    )

    used_slots = set()
    for index, enclosure in enumerate(root_intervals):
        if index == diagonal_root_index:
            all_states_ok &= diagonal_index is not None
            continue
        matching = []
        for slot in reversal_slots:
            left = None if slot == 0 else stationary_intervals[slot - 1]
            right = (
                None
                if slot == len(stationary_intervals)
                else stationary_intervals[slot]
            )
            if (
                (left is None or enclosure.lower() > left.upper())
                and (right is None or enclosure.upper() < right.lower())
            ):
                matching.append(slot)
        all_states_ok &= len(matching) == 1
        if len(matching) == 1:
            all_states_ok &= matching[0] not in used_slots
            used_slots.add(matching[0])
    all_states_ok &= used_slots == set(reversal_slots)

    root_rows = []
    physical = []
    for index, enclosure in enumerate(root_intervals):
        if index == diagonal_root_index:
            root_rows.append(
                {
                    "q": arb_record(enclosure),
                    "diagonal": True,
                    "physical": False,
                }
            )
            continue
        mu_q = arb_mu(enclosure)
        p_gap = arb_p(enclosure) - pi_value
        height = p_gap / (2 * API * mu_q)
        ratio = 2 * mass / mu_q - 1
        endpoint_gap = mu_q - 2 * mass
        signs = (
            strict_sign(height),
            strict_sign(ratio),
            strict_sign(endpoint_gap),
            strict_sign(p_gap),
        )
        all_gates_ok &= all(sign != 0 for sign in signs)
        is_physical = signs[0] > 0 and signs[1] > 0
        root_rows.append(
            {
                "q": arb_record(enclosure),
                "diagonal": False,
                "h": arb_record(height),
                "r": arb_record(ratio),
                "mu_minus_2m": arb_record(endpoint_gap),
                "E_q": arb_record(p_gap),
                "physical": is_physical,
            }
        )
        if is_physical:
            physical.append((enclosure, height, ratio))

    expected = EXPECTED_COUNTS[path]
    all_states_ok &= bool(
        len(root_intervals) == expected[0] and len(physical) == expected[1]
    )

    children = []
    for index, (q_value, height, ratio) in enumerate(physical):
        transition_path = f"{path}/c{index}"
        all_states_ok &= transition_path in EXPECTED_TERMINALS
        next_mass = mass / ratio
        p_gap = arb_p(q_value) - pi_value
        next_pi = arb_p(q_value) + p_gap / ratio
        identity = ratio - (1 + height * q_value)
        identity_width = identity.upper() - identity.lower()
        identity_ok = bool(
            identity.contains(0)
            and identity_width.upper() < IDENTITY_WIDTH_LIMIT.lower()
        )
        all_identity_ok &= identity_ok
        terminal = EXPECTED_TERMINALS.get(transition_path)
        transition_records.append(
            {
                "path": transition_path,
                "identity": arb_record(identity),
                "identity_width": arb_record(identity_width),
                "passed": identity_ok,
            }
        )
        child = {
            "path": transition_path,
            "terminal": terminal,
            "m_plus": arb_record(next_mass),
            "pi_plus": arb_record(next_pi),
        }
        if terminal == "ENTERED_D":
            x_value = mass * q_value
            entry_ok = bool(
                strict_sign(mass) > 0
                and strict_sign(arb(2) / 5 - mass) > 0
                and strict_sign(q_value) > 0
                and strict_sign(x_value - 125) > 0
            )
            hostile_rejected = strict_sign(126 - x_value) > 0
            all_states_ok &= entry_ok and hostile_rejected
            d_entries.append(
                {
                    "path": transition_path,
                    "mass": arb_record(mass),
                    "q": arb_record(q_value),
                    "x": arb_record(x_value),
                    "two_fifths_minus_m": arb_record(arb(2) / 5 - mass),
                    "x_minus_125": arb_record(x_value - 125),
                    "126_minus_x": arb_record(126 - x_value),
                    "passed": entry_ok,
                    "hostile_126_rejected": hostile_rejected,
                }
            )
        elif transition_path in STATIONARY_BRACKETS:
            nested = certify_state(next_mass, next_pi, transition_path)
            child["next"] = nested
            if terminal == "DEAD":
                all_states_ok &= nested["physical_count"] == 0
        children.append(child)

    record = {
        "path": path,
        "mass": arb_record(mass),
        "pi": arb_record(pi_value),
        "p_level_signs": {
            "positive_axis": pos_signs,
            "negative_axis": neg_signs,
            "expected_stationary_count": expected_stationary_count,
        },
        "stationary_bisections": stationary_records,
        "stationary_intervals": [arb_record(v) for v in stationary_intervals],
        "stationary_E": [arb_record(v) for v in stationary_values],
        "left_tail_coefficient": arb_record(left_coefficient),
        "right_tail_coefficient": arb_record(right_coefficient),
        "E_at_zero": arb_record(zero_value),
        "predicted_root_count": predicted_count,
        "root_bisections": root_records,
        "roots": root_rows,
        "all_real_count": len(root_intervals),
        "physical_count": len(physical),
        "children": children,
    }
    state_records[path] = record
    return record


tree = certify_state(M0, PI0, "root")
check(
    "the monotone-factor proof resolves all five fixed-bracket gravity states",
    all_states_ok and len(state_records) == 5,
    f"states={len(state_records)}",
)
check(
    "all non-diagonal physical and endpoint gates are strict",
    all_gates_ok,
)
check(
    "the resolved tree has the frozen DEAD and ENTERED_D terminals",
    all_states_ok and len(d_entries) == 1,
    f"D_entries={len(d_entries)}",
)
check(
    "the endpoint identity is narrow and contains zero on all physical edges",
    all_identity_ok and len(transition_records) == 5,
    f"edges={len(transition_records)}",
)


primary_states = {}


def collect_primary(state):
    primary_states[state["path"]] = state
    for child in state["children"]:
        if child.get("next") is not None:
            collect_primary(child["next"])


collect_primary(primary["tree"])
primary_comparison_ok = set(primary_states) == set(state_records)
comparison_rows = []
for path in sorted(state_records):
    left_state = state_records[path]
    right_state = primary_states.get(path, {})
    row_ok = bool(
        left_state["all_real_count"] == right_state.get("all_real_count")
        and left_state["physical_count"] == right_state.get("physical_count")
        and len(left_state["roots"]) == len(right_state.get("roots", []))
    )
    if row_ok:
        for left, right in zip(left_state["roots"], right_state["roots"]):
            row_ok &= arb(left["q"]["pretty"]).overlaps(
                arb(right["q"]["pretty"])
            )
            row_ok &= left["physical"] == right["physical"]
            if not left["diagonal"]:
                row_ok &= all(
                    left[key]["sign"] == right[key]["sign"]
                    for key in ("h", "r", "mu_minus_2m", "E_q")
                )
    primary_comparison_ok &= row_ok
    comparison_rows.append({"path": path, "passed": row_ok})

first_initial = first_open["tree"]
first_endpoint_comparison_ok = bool(
    first_initial["stationary_bisections"][1]["left_value"]["sign"]
    == tree["stationary_bisections"][1]["initial_left_value"]["sign"]
    and first_initial["stationary_bisections"][1]["right_value"]["sign"]
    == tree["stationary_bisections"][1]["initial_right_value"]["sign"]
    and first_initial["root_bisections"][1]["left_value"]["sign"]
    == tree["root_bisections"][1]["initial_left_value"]["sign"]
    and first_initial["root_bisections"][1]["right_value"]["sign"]
    == tree["root_bisections"][1]["initial_right_value"]["sign"]
)
check(
    "the resolved object agrees with the preserved OPEN endpoints and primary tree",
    primary_comparison_ok and first_endpoint_comparison_ok,
)


hostile_ok = bool(
    len(d_entries) == 1 and d_entries[0]["hostile_126_rejected"]
)
check(
    "the hostile m*q>126 entry claim remains rejected",
    hostile_ok,
)


complete = bool(
    provenance_ok
    and factorization_ok
    and polynomial_controls_ok
    and all_states_ok
    and all_gates_ok
    and all_identity_ok
    and primary_comparison_ok
    and first_endpoint_comparison_ok
    and hostile_ok
)
outcome = (
    "LOCAL_SIGNATURE_ADVERSARIAL_DISAGREEMENT_RESOLVED"
    if complete
    else "LOCAL_SIGNATURE_ADVERSARIAL_DISAGREEMENT_OPEN"
)
check(
    "the exact monotone-factor method resolves the adversarial disagreement",
    outcome == "LOCAL_SIGNATURE_ADVERSARIAL_DISAGREEMENT_RESOLVED",
    outcome,
)


artifact = {
    "provenance": {
        "protocol_commit": PROTOCOL_COMMIT,
        "protocol_sha256": PROTOCOL_SHA256,
        "classification_sha256": CLASSIFICATION_SHA256,
        "invariant_sha256": INVARIANT_SHA256,
        "primary_sha256": PRIMARY_SHA256,
        "first_open_sha256": FIRST_OPEN_SHA256,
    },
    "method": {
        "interval_engine": "python-flint Arb exact endpoint signs",
        "uniqueness_method": "accepted exact K-factor monotonicity of p",
        "arb_decimal_digits": ARB_DPS,
        "maximum_bisection_steps": MAX_BISECTION_STEPS,
        "interval_newton_used": False,
        "raw_wide_derivative_interval_used": False,
        "discovery_root_seeds_used": False,
        "explicit_epsilon_computed": False,
    },
    "polynomial_controls": {
        "positive": poly_record,
        "positive_interval": arb_record(poly_interval),
        "negative": negative_record,
    },
    "tree": tree,
    "transitions": transition_records,
    "D_entries": d_entries,
    "post_construction_primary_comparison": comparison_rows,
    "preserved_open_endpoint_comparison": first_endpoint_comparison_ok,
    "claims": {
        "local_signature": (
            "ADVERSARIALLY_CORROBORATED" if complete else "OPEN"
        ),
        "explicit_radius": "NOT_COMPUTED",
        "global_incoming_basin": "OPEN",
        "physical_selection_rule": "NOT_DERIVED",
        "nonhomogeneous_physics": "NOT_TESTED",
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


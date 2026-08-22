#!/usr/bin/env python3
"""Primary Arb certificate for local stability of the v=3/2 branch tree."""

import gzip
import hashlib
import json
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
    / "gravity_600cell_finite_height_local_signature_protocol.md"
)
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
INVARIANT_INPUT = (
    HERE
    / "gravity_600cell_finite_height_invariant_region_adversarial_resolution.json"
)
DISCOVERY_INPUT = (
    HERE / "gravity_600cell_finite_height_incoming_basin_discovery.json.gz"
)
OUTPUT = HERE / "gravity_600cell_finite_height_local_signature.json"

PROTOCOL_COMMIT = "6512791"
PROTOCOL_SHA256 = (
    "c0a8c7767db9689f6d37236696c151601ce3664aa824dab098a5daab9181bfe7"
)
CLASSIFICATION_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
INVARIANT_SHA256 = (
    "813e05bd66b47cc3ae1cd35d0a2eddb9c645a850d84abeaad37d15b14a6a380f"
)
DISCOVERY_SHA256 = (
    "f492f50cfcaa8e171fb6faa21524d824b4d11b3701b7d635ce483500aaffeb8d"
)

ARB_DPS = 220
MP_DPS = 140
SEED_RADIUS = "1e-38"

ctx.dps = ARB_DPS
mp.mp.dps = MP_DPS
ACTION_LIMIT = mp.mpf("1e-75")
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


def strict_sign(value):
    if value.lower() > 0:
        return 1
    if value.upper() < 0:
        return -1
    return 0


def strict_subset(inner, outer):
    return bool(
        inner.lower() > outer.lower() and inner.upper() < outer.upper()
    )


def disjoint(left, right):
    return bool(left.upper() < right.lower() or right.upper() < left.lower())


def arb_record(value):
    return {
        "pretty": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "contains_zero": bool(value.contains(0)),
        "sign": strict_sign(value),
    }


def mp_midpoint(value):
    return mp.mpf(str(value.mid().lower()))


def seed_ball(value):
    return arb(str(value)) + arb(0, SEED_RADIUS)


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


def arb_p_prime(value):
    jet = arb_series([value, arb(1)], prec=2)
    return arb_p(jet)[1]


def arb_E(mass, pi_value, q_value):
    return 4 * API * (arb_mu(q_value) - mass) + q_value * (
        arb_p(q_value) - pi_value
    )


def interval_newton(function, derivative, seed):
    enclosure = seed_ball(seed)
    history = []
    ok = True
    # One strict interval-Newton inclusion already proves existence and
    # uniqueness.  Further iterations can stall at the outward-rounding floor
    # and are not additional mathematical evidence.
    for _ in range(1):
        derivative_ball = derivative(enclosure)
        if strict_sign(derivative_ball) == 0:
            ok = False
            history.append(
                {
                    "input": arb_record(enclosure),
                    "derivative": arb_record(derivative_ball),
                    "reason": "derivative contains zero",
                }
            )
            break
        center = enclosure.mid()
        image = center - function(center) / derivative_ball
        included = strict_subset(image, enclosure)
        history.append(
            {
                "input": arb_record(enclosure),
                "derivative": arb_record(derivative_ball),
                "image": arb_record(image),
                "strictly_interior": included,
            }
        )
        if not included:
            ok = False
            break
        enclosure = image
    return enclosure, ok, history


classification = json.loads(CLASSIFICATION_INPUT.read_text())
invariant = json.loads(INVARIANT_INPUT.read_text())
with gzip.open(DISCOVERY_INPUT, "rt") as stream:
    discovery = json.load(stream)

provenance_ok = bool(
    digest(PROTOCOL) == PROTOCOL_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and digest(INVARIANT_INPUT) == INVARIANT_SHA256
    and digest(DISCOVERY_INPUT) == DISCOVERY_SHA256
    and classification["outcome"]
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
    and invariant["outcome"]
    == "INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED"
    and discovery["outcome"]
    == "INCOMING_BASIN_CANDIDATE_SKELETON_FROZEN"
)
check(
    "the local protocol and all accepted inputs are frozen",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


# Symbolic identities are checked independently of the numerical tree.
q_symbol, m_symbol, pi_symbol, v_symbol = sp.symbols(
    "q m pi v", real=True
)


def sym_epsilon(value):
    square = value**2
    return 2 * sp.pi - 5 * sp.acos(
        (square + 2) / (2 * (square + 3))
    )


def sym_mu(value):
    return 180 * sym_epsilon(value) / (
        sp.pi * sp.sqrt(value**2 + 4)
    )


def sym_p(value):
    square = value**2
    return (
        180 * value * sym_epsilon(value) / sp.sqrt(square + 4)
        - 600
        * sp.sqrt(3)
        * sp.asinh(value / sp.sqrt(8 * (square + 3)))
    )


E_symbol = 4 * sp.pi * (sym_mu(q_symbol) - m_symbol) + q_symbol * (
    sym_p(q_symbol) - pi_symbol
)
derivative_identity = sp.simplify(
    sp.diff(E_symbol, q_symbol) - (sym_p(q_symbol) - pi_symbol)
)
diagonal_identity = sp.simplify(
    E_symbol.subs(
        {
            m_symbol: sym_mu(v_symbol),
            pi_symbol: sym_p(v_symbol),
            q_symbol: v_symbol,
        }
    )
)
diagonal_derivative_identity = sp.simplify(
    (sym_p(q_symbol) - pi_symbol).subs(
        {pi_symbol: sym_p(v_symbol), q_symbol: v_symbol}
    )
)
symbolic_ok = bool(
    derivative_identity == 0
    and diagonal_identity == 0
    and diagonal_derivative_identity == 0
)
check(
    "the derivative and persistent diagonal identities hold symbolically",
    symbolic_ok,
)


V0 = arb(3) / 2
M0 = arb_mu(V0)
PI0 = arb_p(V0)
P_PRIME_V0 = arb_p_prime(V0)
diagonal_nondegenerate = strict_sign(P_PRIME_V0) != 0
check(
    "the exact initial diagonal tangency is nondegenerate",
    diagonal_nondegenerate,
    str(P_PRIME_V0),
)


thresholds = classification["thresholds"]
V_STAR = seed_ball(thresholds["v_star"])
P_STAR = arb_p(V_STAR)
P_INFINITY = 60 * API - 300 * arb(3).sqrt() * arb(2).log()


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

all_records = []
all_certificates_ok = True
all_gates_ok = True
all_tree_ok = True
d_entry_records = []
transition_records = []


def certify_state(stored_state, mass, pi_value, path, initial=False):
    global all_certificates_ok, all_gates_ok, all_tree_ok

    positive_count, positive_levels = positive_level_count(pi_value)
    negative_count, negative_levels = positive_level_count(-pi_value)
    expected_stationary_count = (
        None
        if positive_count is None or negative_count is None
        else positive_count + negative_count
    )
    stationary_seeds = stored_state["diagnostics"]["stationary_points"]
    stationary_balls = []
    stationary_newton = []
    stationary_ok = expected_stationary_count == len(stationary_seeds)
    for seed in stationary_seeds:
        enclosure, root_ok, history = interval_newton(
            lambda value: arb_p(value) - pi_value,
            arb_p_prime,
            seed,
        )
        stationary_balls.append(enclosure)
        stationary_newton.append(history)
        stationary_ok &= root_ok and strict_sign(arb_p_prime(enclosure)) != 0
    stationary_ok &= all(
        stationary_balls[index].upper()
        < stationary_balls[index + 1].lower()
        for index in range(len(stationary_balls) - 1)
    )

    diagonal_index = None
    if initial:
        matches = [
            index
            for index, seed in enumerate(stationary_seeds)
            if abs(mp.mpf(seed) - mp.mpf("1.5")) < mp.mpf("1e-80")
        ]
        if len(matches) == 1:
            diagonal_index = matches[0]
            diagonal_stationary = stationary_balls[diagonal_index]
            stationary_ok &= bool(
                diagonal_stationary.lower() <= V0.lower()
                and diagonal_stationary.upper() >= V0.upper()
            )
        else:
            stationary_ok = False

    stationary_values = []
    stationary_signs = []
    for index, enclosure in enumerate(stationary_balls):
        value = arb_E(mass, pi_value, enclosure)
        stationary_values.append(value)
        if index == diagonal_index:
            stationary_signs.append(0)
        else:
            stationary_signs.append(strict_sign(value))
            stationary_ok &= strict_sign(value) != 0

    right_tail = P_INFINITY - pi_value
    left_coefficient = -P_INFINITY - pi_value
    left_tail_sign = -strict_sign(left_coefficient)
    right_tail_sign = strict_sign(right_tail)
    central = 4 * API * (arb_mu(arb(0)) - mass)
    tail_origin_ok = bool(
        left_tail_sign != 0
        and right_tail_sign != 0
        and strict_sign(central) != 0
    )

    boundary_signs = [left_tail_sign, *stationary_signs, right_tail_sign]
    reversal_slots = []
    for index in range(len(boundary_signs) - 1):
        if boundary_signs[index] * boundary_signs[index + 1] < 0:
            reversal_slots.append(index)
    predicted_root_count = len(reversal_slots) + int(diagonal_index is not None)

    root_seeds = stored_state["all_real_roots"]
    root_balls = []
    root_newton = []
    diagonal_root_index = None
    roots_ok = predicted_root_count == len(root_seeds)
    for index, seed in enumerate(root_seeds):
        if initial and abs(mp.mpf(seed) - mp.mpf("1.5")) < mp.mpf("1e-80"):
            enclosure = V0
            history = [{"exact_diagonal": True}]
            root_ok = diagonal_root_index is None
            diagonal_root_index = index
        else:
            enclosure, root_ok, history = interval_newton(
                lambda value: arb_E(mass, pi_value, value),
                lambda value: arb_p(value) - pi_value,
                seed,
            )
            root_ok &= strict_sign(arb_p(enclosure) - pi_value) != 0
        root_balls.append(enclosure)
        root_newton.append(history)
        roots_ok &= root_ok
    roots_ok &= all(
        root_balls[index].upper() < root_balls[index + 1].lower()
        for index in range(len(root_balls) - 1)
    )

    used_slots = set()
    for index, enclosure in enumerate(root_balls):
        if index == diagonal_root_index:
            roots_ok &= diagonal_index is not None
            continue
        matching_slots = []
        for slot in reversal_slots:
            left = None if slot == 0 else stationary_balls[slot - 1]
            right = (
                None
                if slot == len(stationary_balls)
                else stationary_balls[slot]
            )
            inside = bool(
                (left is None or enclosure.lower() > left.upper())
                and (right is None or enclosure.upper() < right.lower())
            )
            if inside:
                matching_slots.append(slot)
        roots_ok &= len(matching_slots) == 1
        if len(matching_slots) == 1:
            roots_ok &= matching_slots[0] not in used_slots
            used_slots.add(matching_slots[0])
    roots_ok &= used_slots == set(reversal_slots)

    root_rows = []
    physical = []
    root_census = stored_state["root_census"]
    roots_ok &= len(root_census) == len(root_balls)
    for index, enclosure in enumerate(root_balls):
        if index == diagonal_root_index:
            row = {
                "q": arb_record(enclosure),
                "diagonal_zero_height": True,
                "h": arb_record(arb(0)),
                "r": arb_record(arb(1)),
                "mu_minus_2m": arb_record(arb_mu(enclosure) - 2 * mass),
                "E_q": arb_record(arb(0)),
                "physical": False,
            }
            expected_physical = False
        else:
            mu_q = arb_mu(enclosure)
            p_gap = arb_p(enclosure) - pi_value
            height = p_gap / (2 * API * mu_q)
            ratio = 2 * mass / mu_q - 1
            endpoint_gap = mu_q - 2 * mass
            derivative = p_gap
            signs = [
                strict_sign(height),
                strict_sign(ratio),
                strict_sign(endpoint_gap),
                strict_sign(derivative),
            ]
            gates_resolved = all(value != 0 for value in signs)
            all_gates_ok &= gates_resolved
            is_physical = bool(signs[0] > 0 and signs[1] > 0)
            expected_physical = bool(root_census[index]["physical"])
            all_gates_ok &= is_physical == expected_physical
            row = {
                "q": arb_record(enclosure),
                "diagonal_zero_height": False,
                "h": arb_record(height),
                "r": arb_record(ratio),
                "mu_minus_2m": arb_record(endpoint_gap),
                "E_q": arb_record(derivative),
                "physical": is_physical,
                "expected_physical": expected_physical,
            }
            if is_physical:
                physical.append((enclosure, height, ratio, row))
        root_rows.append(row)

    expected_counts = EXPECTED_COUNTS[path]
    state_count_ok = bool(
        len(root_balls) == expected_counts[0]
        and len(physical) == expected_counts[1]
        and stored_state["all_real_count"] == expected_counts[0]
        and stored_state["physical_count"] == expected_counts[1]
    )
    all_tree_ok &= state_count_ok
    all_certificates_ok &= stationary_ok and tail_origin_ok and roots_ok

    children = stored_state["children"]
    all_tree_ok &= len(children) == len(physical)
    child_records = []
    for index, (q_ball, height, ratio, _) in enumerate(physical):
        if index >= len(children):
            break
        child = children[index]
        transition_path = f"{path}/c{index}"
        all_tree_ok &= (
            transition_path in EXPECTED_TERMINALS
            and child["terminal"] == EXPECTED_TERMINALS[transition_path]
        )
        next_mass = mass / ratio
        p_gap = arb_p(q_ball) - pi_value
        next_pi = arb_p(q_ball) + p_gap / ratio
        transition = {
            "path": transition_path,
            "m_plus": arb_record(next_mass),
            "pi_plus": arb_record(next_pi),
            "terminal": child["terminal"],
        }
        transition_records.append(
            {
                "path": f"{path}/c{index}",
                "mass": mass,
                "pi": pi_value,
                "q": q_ball,
                "h": height,
                "r": ratio,
                "m_plus": next_mass,
                "pi_plus": next_pi,
            }
        )
        if child["terminal"] == "ENTERED_D":
            x_value = mass * q_ball
            d_ok = bool(
                strict_sign(mass) > 0
                and strict_sign(arb(2) / 5 - mass) > 0
                and strict_sign(q_ball) > 0
                and strict_sign(x_value - 125) > 0
            )
            hostile_126_rejected = strict_sign(126 - x_value) > 0
            all_tree_ok &= d_ok and hostile_126_rejected
            d_entry_records.append(
                {
                    "path": f"{path}/c{index}",
                    "mass": arb_record(mass),
                    "q": arb_record(q_ball),
                    "x": arb_record(x_value),
                    "two_fifths_minus_m": arb_record(arb(2) / 5 - mass),
                    "x_minus_125": arb_record(x_value - 125),
                    "126_minus_x": arb_record(126 - x_value),
                    "passed": d_ok,
                    "hostile_126_rejected": hostile_126_rejected,
                }
            )
        if child.get("next") is not None:
            nested = certify_state(
                child["next"],
                next_mass,
                next_pi,
                transition_path,
                initial=False,
            )
            transition["next"] = nested
            if child["terminal"] == "DEAD":
                all_tree_ok &= nested["physical_count"] == 0
        child_records.append(transition)

    record = {
        "path": path,
        "mass": arb_record(mass),
        "pi": arb_record(pi_value),
        "p_level_signs": {
            "positive_axis": positive_levels,
            "negative_axis": negative_levels,
            "expected_stationary_count": expected_stationary_count,
        },
        "stationary_points": [arb_record(value) for value in stationary_balls],
        "stationary_newton": stationary_newton,
        "stationary_E": [arb_record(value) for value in stationary_values],
        "left_tail_coefficient": arb_record(left_coefficient),
        "right_tail_coefficient": arb_record(right_tail),
        "E_at_zero": arb_record(central),
        "predicted_root_count": predicted_root_count,
        "root_newton": root_newton,
        "roots": root_rows,
        "all_real_count": len(root_balls),
        "physical_count": len(physical),
        "children": child_records,
        "certificate_passed": bool(
            stationary_ok and tail_origin_ok and roots_ok and state_count_ok
        ),
    }
    all_records.append(record)
    return record


stored_reference = discovery["post_construction_controls"][
    "known_representatives"
]["v_3_over_2"]
root_record = certify_state(stored_reference, M0, PI0, "root", initial=True)

check(
    "every visited state has a complete stationary-and-tail root certificate",
    all_certificates_ok and len(all_records) == 5,
    f"states={len(all_records)}",
)
check(
    "every non-diagonal real root has strict physical-gate signs",
    all_gates_ok,
)
check(
    "the ordered recursive tree is one, then two, with DEAD and ENTERED_D terminals",
    all_tree_ok and len(d_entry_records) == 1,
    f"D_entries={len(d_entry_records)}",
)


# The stored decimal values are controls only: every Newton result must remain
# inside its preregistered discovery-seed ball.
seed_control_ok = True
for state in all_records:
    for history in state["stationary_newton"] + state["root_newton"]:
        if history and not history[0].get("exact_diagonal"):
            seed_control_ok &= bool(history[0].get("strictly_interior", False))
check(
    "the certified roots reproduce the delayed discovery seeds",
    seed_control_ok,
)


# Redifferentiate the complete action only after the primary tree exists.
L0, L1, RHO, MASS = sp.symbols("L0 L1 RHO MASS", positive=True)
DELTA = L1 - L0
HEIGHT = sp.sqrt(RHO + DELTA**2 / 4)
COSINE = (DELTA**2 + 2 * RHO) / (2 * (DELTA**2 + 3 * RHO))
BOOST = DELTA / sp.sqrt(8 * (DELTA**2 + 3 * RHO))
ACTION = (
    360
    * (L0 + L1)
    * HEIGHT
    * (2 * sp.pi - 5 * sp.acos(COSINE))
    + 600
    * sp.sqrt(3)
    * (L0**2 - L1**2)
    * sp.asinh(BOOST)
    - 8 * sp.pi * MASS * sp.sqrt(RHO)
)
CONSTRAINT = RHO * sp.diff(ACTION, RHO)
PRE = -L0 * sp.diff(ACTION, L0) / 2
POST = L1 * sp.diff(ACTION, L1) / 2
constraint_numeric = sp.lambdify((L0, L1, RHO, MASS), CONSTRAINT, "mpmath")
pre_numeric = sp.lambdify((L0, L1, RHO, MASS), PRE, "mpmath")
post_numeric = sp.lambdify((L0, L1, RHO, MASS), POST, "mpmath")

action_rows = []
action_ok = True
for transition in transition_records:
    mass = mp_midpoint(transition["mass"])
    pi_value = mp_midpoint(transition["pi"])
    q_value = mp_midpoint(transition["q"])
    height = mp_midpoint(transition["h"])
    ratio = mp_midpoint(transition["r"])
    next_pi = mp_midpoint(transition["pi_plus"])
    constraint_residual = 2 * constraint_numeric(
        1, ratio, height**2, mass
    ) / height
    pre_residual = pre_numeric(1, ratio, height**2, mass) - pi_value
    post_residual = post_numeric(1, ratio, height**2, mass) / ratio**2 - next_pi
    row_ok = max(
        abs(constraint_residual), abs(pre_residual), abs(post_residual)
    ) < ACTION_LIMIT
    action_ok &= row_ok
    action_rows.append(
        {
            "path": transition["path"],
            "constraint": mp.nstr(constraint_residual, 90),
            "pre": mp.nstr(pre_residual, 90),
            "post": mp.nstr(post_residual, 90),
            "passed": row_ok,
        }
    )
check(
    "every physical edge satisfies the redifferentiated complete action",
    action_ok and len(action_rows) == 5,
    f"edges={len(action_rows)}",
)


first_transition = transition_records[0]
correct_pi = first_transition["pi_plus"]
wrong_scaled_pi = correct_pi * first_transition["r"]
wrong_reversed_pi = -correct_pi
hostile_momentum_ok = bool(
    disjoint(correct_pi, wrong_scaled_pi)
    and disjoint(correct_pi, wrong_reversed_pi)
)
hostile_d_ok = bool(
    len(d_entry_records) == 1
    and d_entry_records[0]["hostile_126_rejected"]
)
check(
    "wrong post-momentum conventions and the hostile x>126 gate are rejected",
    hostile_momentum_ok and hostile_d_ok,
)


local_lemma_hypotheses = bool(
    provenance_ok
    and symbolic_ok
    and diagonal_nondegenerate
    and all_certificates_ok
    and all_gates_ok
    and all_tree_ok
    and seed_control_ok
    and action_ok
    and hostile_momentum_ok
    and hostile_d_ok
)
outcome = (
    "LOCAL_SIGNATURE_PRIMARY_CERTIFIED"
    if local_lemma_hypotheses
    else "LOCAL_SIGNATURE_OPEN"
)
check(
    "the complete local-constancy hypotheses imply an unspecified open neighbourhood",
    outcome == "LOCAL_SIGNATURE_PRIMARY_CERTIFIED",
    outcome,
)


artifact = {
    "provenance": {
        "protocol_commit": PROTOCOL_COMMIT,
        "protocol_sha256": PROTOCOL_SHA256,
        "classification_sha256": CLASSIFICATION_SHA256,
        "invariant_sha256": INVARIANT_SHA256,
        "discovery_sha256": DISCOVERY_SHA256,
    },
    "method": {
        "interval_engine": "python-flint Arb outward-rounded balls",
        "arb_decimal_digits": ARB_DPS,
        "mpmath_control_digits": MP_DPS,
        "seed_radius": SEED_RADIUS,
        "explicit_epsilon_computed": False,
        "grid_or_radius_search_used": False,
        "root_completeness_method": "stationary partition plus analytic tails",
    },
    "symbolic": {
        "E_q_equals_p_minus_pi": derivative_identity == 0,
        "diagonal_E_identity": diagonal_identity == 0,
        "diagonal_derivative_identity": diagonal_derivative_identity == 0,
        "p_prime_at_three_halves": arb_record(P_PRIME_V0),
    },
    "tree": root_record,
    "D_entries": d_entry_records,
    "controls": {
        "discovery_seed_newton_inclusions": seed_control_ok,
        "action_rows": action_rows,
        "wrong_momentum_rejected": hostile_momentum_ok,
        "hostile_x_126_rejected": hostile_d_ok,
    },
    "logical_implication": {
        "persistent_nondegenerate_diagonal_tangency": bool(
            symbolic_ok and diagonal_nondegenerate
        ),
        "all_other_roots_simple": all_certificates_ok,
        "no_local_root_birth_or_escape": all_certificates_ok,
        "strict_physical_and_terminal_gates": bool(
            all_gates_ok and all_tree_ok
        ),
        "conclusion": (
            "there exists an unspecified epsilon>0 with constant ordered "
            "physical tree on the incoming curve"
        ),
    },
    "claims": {
        "local_signature": "PRIMARY_CERTIFICATE_PENDING_ADVERSARIAL_REPLICATION",
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

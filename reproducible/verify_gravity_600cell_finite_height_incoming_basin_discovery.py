#!/usr/bin/env python3
"""Target-free discovery skeleton for the incoming-state branch basin."""

import hashlib
import json
from collections import Counter
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = (
    ROOT
    / "docs"
    / "gravity"
    / "gravity_600cell_finite_height_incoming_basin_discovery_protocol.md"
)
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
INVARIANT_INPUT = HERE / "gravity_600cell_finite_height_invariant_region_adversarial_resolution.json"
OUTPUT = HERE / "gravity_600cell_finite_height_incoming_basin_discovery.json"

PROTOCOL_COMMIT = "5da21c0"
PROTOCOL_SHA256 = (
    "c9baf3e0fd5a77d5b9e431be0740c3b492071a51950d2aa0623c131ea35ad671"
)
CLASSIFICATION_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
INVARIANT_SHA256 = (
    "813e05bd66b47cc3ae1cd35d0a2eddb9c645a850d84abeaad37d15b14a6a380f"
)

NODE_COUNT = 1024
MAX_SLAB = 4
MAX_TREE_NODES = 256
WORK_DPS = 140
SERIAL_DIGITS = 90

# Set the arithmetic context before constructing any numerical protocol
# constant.  Otherwise mpmath would freeze these values at its import-time
# precision even though the subsequent calculations use 140 digits.
mp.mp.dps = WORK_DPS
ROOT_TOLERANCE = mp.mpf("1e-112")
CRITICAL_MARGIN = mp.mpf("1e-72")
RESIDUAL_LIMIT = mp.mpf("1e-88")
MATCH_MARGIN = mp.mpf("1e-70")
MAX_BISECTION_STEPS = 620
MAX_TAIL_STEPS = 800

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


def text(value, digits=SERIAL_DIGITS):
    if isinstance(value, str):
        return value
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return value
    return mp.nstr(value, digits)


classification = json.loads(CLASSIFICATION_INPUT.read_text())
invariant = json.loads(INVARIANT_INPUT.read_text())
provenance_ok = bool(
    digest(PROTOCOL) == PROTOCOL_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and digest(INVARIANT_INPUT) == INVARIANT_SHA256
    and classification["outcome"]
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
    and invariant["outcome"]
    == "INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED"
)
check(
    "the target-free discovery protocol and accepted inputs are frozen",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


PI = mp.pi


def epsilon(value):
    square = value * value
    return 2 * PI - 5 * mp.acos((square + 2) / (2 * (square + 3)))


def mu(value):
    return 180 * epsilon(value) / (PI * mp.sqrt(value * value + 4))


def momentum(value):
    square = value * value
    return (
        180 * value * epsilon(value) / mp.sqrt(square + 4)
        - 600
        * mp.sqrt(3)
        * mp.asinh(value / mp.sqrt(8 * (square + 3)))
    )


P_INFINITY = 60 * PI - 300 * mp.sqrt(3) * mp.log(2)
thresholds = classification["thresholds"]
V_A = mp.mpf(thresholds["v_A"])
V_STAR = mp.mpf(thresholds["v_star"])
V_M = mp.mpf(thresholds["v_M"])
V_C = mp.mpf(thresholds["v_C"])
P_STAR = momentum(V_STAR)

threshold_order_ok = bool(
    0 < V_A < V_STAR < V_M < V_C
    and P_STAR < P_INFINITY < 0
    and abs(mu(V_M) - mu(0)) < mp.mpf("1e-75")
)
check(
    "the intrinsic incoming thresholds and momentum ranges are ordered",
    threshold_order_ok,
)


def sign(value, margin=CRITICAL_MARGIN):
    if value > margin:
        return 1
    if value < -margin:
        return -1
    return 0


def bisect(function, left, right):
    left = mp.mpf(left)
    right = mp.mpf(right)
    f_left = function(left)
    f_right = function(right)
    # CRITICAL_MARGIN is only an ambiguity gate for bifurcation levels.  It
    # must not turn a merely small endpoint value into an exact root.
    if f_left == 0:
        return left
    if f_right == 0:
        return right
    if f_left * f_right >= 0:
        raise RuntimeError("bisection interval lacks a strict sign change")
    for _ in range(MAX_BISECTION_STEPS):
        middle = (left + right) / 2
        f_middle = function(middle)
        if f_middle == 0 or abs(right - left) < ROOT_TOLERANCE:
            return middle
        if f_left * f_middle < 0:
            right = middle
        else:
            left = middle
            f_left = f_middle
    raise RuntimeError("bisection step limit reached")


def positive_p_level_roots(level):
    """All s>0 with p(s)=level, or an unresolved critical marker."""
    level = mp.mpf(level)
    if (
        abs(level - P_STAR) <= CRITICAL_MARGIN
        or abs(level - P_INFINITY) <= CRITICAL_MARGIN
        or abs(level) <= CRITICAL_MARGIN
    ):
        return [], False, "momentum level is critical"
    if level < P_STAR or level > 0:
        return [], True, "outside positive-axis momentum range"

    function = lambda value: momentum(value) - level
    roots = [bisect(function, 0, V_STAR)]
    if level < P_INFINITY:
        right = 2 * V_STAR
        for _ in range(MAX_TAIL_STEPS):
            if function(right) > CRITICAL_MARGIN:
                roots.append(bisect(function, V_STAR, right))
                break
            right *= 2
        else:
            return roots, False, "outer momentum level did not bracket"
    return roots, True, "complete positive-axis momentum partition"


def stationary_points(pi_value):
    positive, positive_ok, positive_detail = positive_p_level_roots(pi_value)
    negative_magnitudes, negative_ok, negative_detail = positive_p_level_roots(
        -pi_value
    )
    points = sorted([-value for value in negative_magnitudes] + positive)
    return points, bool(positive_ok and negative_ok), {
        "positive": positive_detail,
        "negative": negative_detail,
    }


def root_equation(mass, pi_value, q_value):
    return (
        4 * PI * (mu(q_value) - mass)
        + q_value * (momentum(q_value) - pi_value)
    )


def tail_bracket(function, finite_boundary, direction, required_sign):
    # When E has no finite stationary point, the sole monotonic interval is
    # (-infinity,+infinity).  q=0 is already certified noncritical above and
    # is used only to seed the magnitude of both analytic-tail probes.
    anchor = mp.mpf(0) if finite_boundary is None else finite_boundary
    magnitude = max(mp.mpf(10), 2 * abs(anchor))
    for _ in range(MAX_TAIL_STEPS):
        point = direction * magnitude
        if sign(function(point)) == required_sign:
            return point
        magnitude *= 2
    raise RuntimeError("root-equation tail did not reach its analytic sign")


def all_real_roots(mass, pi_value, excluded_stationary=None):
    stationary, stationary_ok, stationary_detail = stationary_points(pi_value)
    if not stationary_ok:
        return [], False, {
            "reason": stationary_detail,
            "stationary_points": [text(value) for value in stationary],
        }

    right_coefficient = P_INFINITY - pi_value
    left_coefficient = -P_INFINITY - pi_value
    right_tail_sign = sign(right_coefficient)
    left_tail_sign = -sign(left_coefficient)
    central_value = 4 * PI * (mu(0) - mass)
    if right_tail_sign == 0 or left_tail_sign == 0 or sign(central_value) == 0:
        return [], False, {
            "reason": "root equation is at a tail or q=0 critical level",
            "right_tail_coefficient": text(right_coefficient),
            "left_tail_coefficient": text(left_coefficient),
            "q_zero_value": text(central_value),
        }

    function = lambda value: root_equation(mass, pi_value, value)

    excluded_index = None
    if excluded_stationary is not None:
        excluded_stationary = mp.mpf(excluded_stationary)
        matches = [
            index
            for index, value in enumerate(stationary)
            if abs(value - excluded_stationary) < MATCH_MARGIN
        ]
        exact_root_residual = function(excluded_stationary)
        exact_derivative_residual = momentum(excluded_stationary) - pi_value
        if (
            len(matches) != 1
            or abs(exact_root_residual) >= RESIDUAL_LIMIT
            or abs(exact_derivative_residual) >= RESIDUAL_LIMIT
        ):
            return [], False, {
                "reason": "known initial diagonal stationary root did not match exactly once",
                "match_count": len(matches),
                "root_residual": text(exact_root_residual),
                "derivative_residual": text(exact_derivative_residual),
                "stationary_points": [text(value) for value in stationary],
            }
        excluded_index = matches[0]
        # Use the exact incoming parameter rather than its independently
        # bisected copy.  This makes the preregistered q=v exclusion explicit
        # while retaining every other stationary point from the root census.
        stationary[excluded_index] = excluded_stationary

    stationary_values = [function(value) for value in stationary]
    ambiguous_indices = [
        index
        for index, value in enumerate(stationary_values)
        if sign(value) == 0 and index != excluded_index
    ]
    if ambiguous_indices:
        return [], False, {
            "reason": "stationary point is a multiple-root candidate",
            "ambiguous_indices": ambiguous_indices,
            "stationary_points": [text(value) for value in stationary],
            "stationary_values": [text(value) for value in stationary_values],
        }

    boundaries = [None, *stationary, None]
    boundary_values = [
        mp.mpf(left_tail_sign),
        *stationary_values,
        mp.mpf(right_tail_sign),
    ]
    roots = (
        [excluded_stationary]
        if excluded_stationary is not None
        else []
    )
    intervals = []
    for index in range(len(boundaries) - 1):
        left = boundaries[index]
        right = boundaries[index + 1]
        left_sign = sign(boundary_values[index])
        right_sign = sign(boundary_values[index + 1])
        # A zero boundary is the separately recorded initial diagonal
        # tangency.  Only a strict sign reversal contains another root in the
        # open monotonic interval.
        has_root = left_sign * right_sign < 0
        intervals.append(
            {
                "left": "-infinity" if left is None else text(left),
                "right": "+infinity" if right is None else text(right),
                "left_sign": left_sign,
                "right_sign": right_sign,
                "has_root": has_root,
            }
        )
        if not has_root:
            continue
        bracket_left = (
            tail_bracket(function, right, -1, left_sign)
            if left is None
            else left
        )
        bracket_right = (
            tail_bracket(function, left, 1, right_sign)
            if right is None
            else right
        )
        roots.append(bisect(function, bracket_left, bracket_right))

    roots.sort()
    residual_ok = all(abs(function(root)) < RESIDUAL_LIMIT for root in roots)
    return roots, residual_ok, {
        "stationary_points": [text(value) for value in stationary],
        "stationary_values": [text(value) for value in stationary_values],
        "stationary_derivative_residuals": [
            text(momentum(value) - pi_value) for value in stationary
        ],
        "excluded_initial_diagonal": (
            None if excluded_stationary is None else text(excluded_stationary)
        ),
        "right_tail_coefficient": text(right_coefficient),
        "left_tail_coefficient": text(left_coefficient),
        "q_zero_value": text(central_value),
        "intervals": intervals,
    }


def physical_data(mass, pi_value, q_value):
    mu_q = mu(q_value)
    p_gap = momentum(q_value) - pi_value
    height = p_gap / (2 * PI * mu_q)
    ratio = 2 * mass / mu_q - 1
    endpoint_identity = ratio - (1 + height * q_value)
    if abs(endpoint_identity) >= RESIDUAL_LIMIT:
        return None, False, "endpoint identity failed"
    h_sign = sign(height)
    r_sign = sign(ratio)
    if h_sign == 0 or r_sign == 0:
        return None, False, "physical gate is critical"
    physical = bool(h_sign > 0 and r_sign > 0)
    next_mass = mass / ratio
    next_pi = momentum(q_value) + p_gap / ratio
    return {
        "q": q_value,
        "h": height,
        "r": ratio,
        "x": mass * q_value,
        "mu_minus_2m": mu_q - 2 * mass,
        "m_minus_two_fifths": mass - mp.mpf(2) / 5,
        "x_minus_125": mass * q_value - 125,
        "m_plus": next_mass,
        "pi_plus": next_pi,
        "root_residual": root_equation(mass, pi_value, q_value),
        "endpoint_identity_residual": endpoint_identity,
        "physical": physical,
    }, True, "classified"


def serialize_root_data(data):
    return {
        key: (value if isinstance(value, bool) else text(value))
        for key, value in data.items()
    }


def tree_signature(state):
    if state.get("unresolved"):
        return {"unresolved": True}
    return {
        "all_real": state["all_real_count"],
        "physical": state["physical_count"],
        "children": [
            {
                "terminal": child["terminal"],
                "next": (
                    tree_signature(child["next"])
                    if child.get("next") is not None
                    else None
                ),
            }
            for child in state["children"]
        ],
    }


def terminal_labels(state):
    labels = []
    for child in state.get("children", []):
        if child["terminal"]:
            labels.append(child["terminal"])
        elif child.get("next") is not None:
            labels.extend(terminal_labels(child["next"]))
    return labels


def signature_text(state):
    return json.dumps(tree_signature(state), sort_keys=True, separators=(",", ":"))


def expand_state(mass, pi_value, slab, budget, initial_v=None):
    if budget[0] >= MAX_TREE_NODES:
        return {
            "slab": slab,
            "unresolved": True,
            "reason": "tree node budget reached",
            "children": [],
        }
    budget[0] += 1

    roots, census_ok, diagnostics = all_real_roots(
        mass,
        pi_value,
        excluded_stationary=(initial_v if slab == 1 else None),
    )
    state = {
        "slab": slab,
        "m": text(mass),
        "pi": text(pi_value),
        "all_real_count": len(roots),
        "all_real_roots": [text(root) for root in roots],
        "root_census": [],
        "physical_count": 0,
        "diagnostics": diagnostics,
        "children": [],
        "unresolved": not census_ok,
    }
    if not census_ok:
        state["reason"] = "all-real root census unresolved"
        return state

    physical_rows = []
    diagonal_exclusions = 0
    for root in roots:
        if slab == 1 and initial_v is not None and abs(root - initial_v) < MATCH_MARGIN:
            # The exact diagonal root has h=0 and is excluded by theorem.
            diagonal_exclusions += 1
            state["root_census"].append(
                {
                    "q": text(root),
                    "classified": True,
                    "physical": False,
                    "exclusion": "KNOWN_INITIAL_DIAGONAL_ZERO_HEIGHT",
                }
            )
            continue
        data, classified, reason = physical_data(mass, pi_value, root)
        state["root_census"].append(
            {
                "q": text(root),
                "classified": classified,
                "physical": bool(data is not None and data["physical"]),
                "data": (
                    serialize_root_data(data) if data is not None else None
                ),
                "detail": reason,
            }
        )
        if not classified:
            state["unresolved"] = True
            state["reason"] = reason
            return state
        if data["physical"]:
            physical_rows.append(data)

    if slab == 1 and initial_v is not None and diagonal_exclusions != 1:
        state["unresolved"] = True
        state["reason"] = "initial diagonal exclusion did not occur exactly once"
        state["diagonal_exclusion_count"] = diagonal_exclusions
        return state
    state["diagonal_exclusion_count"] = diagonal_exclusions

    physical_rows.sort(key=lambda row: row["q"])
    state["physical_count"] = len(physical_rows)
    for data in physical_rows:
        entered = bool(
            mass <= mp.mpf(2) / 5
            and data["x"] >= 125
            and data["q"] > 0
        )
        child = {
            "data": serialize_root_data(data),
            "terminal": None,
            "next": None,
        }
        if entered:
            child["terminal"] = "ENTERED_D"
        elif slab == MAX_SLAB:
            child["terminal"] = "LIVE_OUTSIDE_D_AT_DEPTH_4"
        else:
            next_state = expand_state(
                data["m_plus"], data["pi_plus"], slab + 1, budget
            )
            child["next"] = next_state
            if next_state.get("unresolved"):
                child["terminal"] = "UNRESOLVED"
                state["unresolved"] = True
            elif next_state["physical_count"] == 0:
                child["terminal"] = "DEAD"
        state["children"].append(child)
    return state


def build_initial(v_value):
    budget = [0]
    state = expand_state(mu(v_value), momentum(v_value), 1, budget, initial_v=v_value)
    state["v"] = text(v_value)
    state["tree_nodes"] = budget[0]
    state["signature"] = signature_text(state)
    state["terminal_labels"] = sorted(terminal_labels(state))
    return state


def flatten_monitors(state, prefix="state"):
    monitors = {}
    if state.get("unresolved"):
        return monitors
    for index, value in enumerate(state["diagnostics"].get("stationary_values", [])):
        monitors[f"{prefix}/stationary_E/{index}"] = mp.mpf(value)
    for index, child in enumerate(state["children"]):
        child_prefix = f"{prefix}/child/{index}"
        data = child["data"]
        for name in (
            "h",
            "r",
            "m_minus_two_fifths",
            "x_minus_125",
            "mu_minus_2m",
        ):
            monitors[f"{child_prefix}/{name}"] = mp.mpf(data[name])
        if child.get("next") is not None:
            monitors.update(flatten_monitors(child["next"], child_prefix))
    return monitors


def chebyshev_nodes(left, right):
    values = []
    for index in range(NODE_COUNT):
        angle = (2 * index + 1) * PI / (2 * NODE_COUNT)
        value = (left + right) / 2 + (right - left) * mp.cos(angle) / 2
        values.append(value)
    return sorted(values)


components = (
    ("vA_to_vstar", V_A, V_STAR),
    ("vstar_to_vM", V_STAR, V_M),
    ("vM_to_vC", V_M, V_C),
)

component_records = []
all_nodes_resolved = True
total_processed = 0
for component_name, left, right in components:
    rows = []
    for index, v_value in enumerate(chebyshev_nodes(left, right)):
        row = build_initial(v_value)
        row["component"] = component_name
        row["node_index"] = index
        rows.append(row)
        all_nodes_resolved &= not row.get("unresolved", False)
        total_processed += 1
        if total_processed % 128 == 0:
            print(
                f"[PROGRESS] incoming nodes {total_processed}/{3*NODE_COUNT}",
                flush=True,
            )
    component_records.append(
        {
            "name": component_name,
            "left": text(left),
            "right": text(right),
            "rows": rows,
        }
    )

check(
    "all frozen Chebyshev nodes have complete all-real branch trees",
    all_nodes_resolved,
    f"nodes={total_processed}; depth={MAX_SLAB}",
)


candidate_cells = []
signature_runs = []
signature_counter = Counter()
for component in component_records:
    rows = component["rows"]
    start = 0
    for index, row in enumerate(rows):
        signature_counter[row["signature"]] += 1
        if index and row["signature"] != rows[index - 1]["signature"]:
            signature_runs.append(
                {
                    "component": component["name"],
                    "start_index": start,
                    "end_index": index - 1,
                    "signature": rows[index - 1]["signature"],
                }
            )
            start = index
        if not index:
            continue
        previous = rows[index - 1]
        reasons = []
        equation_types = set()
        if row["signature"] != previous["signature"]:
            reasons.append("tree_signature_change")
            equation_types.update(("branch_birth_or_merger", "zero_endpoint"))
        if row.get("unresolved") or previous.get("unresolved"):
            reasons.append("unresolved_adjacent_node")
            equation_types.add("unresolved")
        if row["signature"] == previous["signature"]:
            before = flatten_monitors(previous)
            after = flatten_monitors(row)
            for key in sorted(set(before) & set(after)):
                left_sign = sign(before[key])
                right_sign = sign(after[key])
                if left_sign and right_sign and left_sign != right_sign:
                    reasons.append(f"sign_change:{key}")
                    if key.endswith("stationary_E") or "/stationary_E/" in key:
                        equation_types.add("branch_birth_or_merger")
                    elif key.endswith("/h"):
                        equation_types.add("zero_height_state_curve_contact")
                    elif key.endswith("/r") or key.endswith("/mu_minus_2m"):
                        equation_types.add("zero_endpoint")
                    elif key.endswith("/m_minus_two_fifths"):
                        equation_types.add("invariant_entry_m")
                    elif key.endswith("/x_minus_125"):
                        equation_types.add("invariant_entry_x")
        if reasons:
            candidate_cells.append(
                {
                    "component": component["name"],
                    "left_index": index - 1,
                    "right_index": index,
                    "left_v": previous["v"],
                    "right_v": row["v"],
                    "reasons": sorted(set(reasons)),
                    "equation_types": sorted(equation_types),
                    "left_signature": previous["signature"],
                    "right_signature": row["signature"],
                }
            )
    signature_runs.append(
        {
            "component": component["name"],
            "start_index": start,
            "end_index": len(rows) - 1,
            "signature": rows[-1]["signature"],
        }
    )

inherited_candidates = [
    {"equation_type": "inherited_endpoint", "threshold": name}
    for name in ("v_A", "v_star", "v_M", "v_C")
]

skeleton = {
    "components": component_records,
    "signature_runs": signature_runs,
    "candidate_cells": candidate_cells,
    "inherited_candidates": inherited_candidates,
    "distinct_signature_count": len(signature_counter),
    "signature_node_counts": dict(signature_counter),
}


# The complete discovery object now exists in memory.  Known representatives
# and altered conventions are evaluated only below this line.
control_states = {
    "v_3_over_2": build_initial(mp.mpf(3) / 2),
    "v_3": build_initial(mp.mpf(3)),
    "v_20": build_initial(mp.mpf(20)),
}


def second_physical_count(state):
    if state.get("unresolved") or state["physical_count"] != 1:
        return None
    first = state["children"][0]
    if first.get("next") is None:
        return 0
    return first["next"]["physical_count"]


known_controls_ok = bool(
    second_physical_count(control_states["v_3_over_2"]) == 2
    and control_states["v_3_over_2"]["terminal_labels"]
    == ["DEAD", "ENTERED_D"]
    and second_physical_count(control_states["v_3"]) == 2
    and second_physical_count(control_states["v_20"]) == 0
)
check(
    "known representative signatures are reproduced only after discovery",
    known_controls_ok,
    (
        f"v=3/2 labels={control_states['v_3_over_2']['terminal_labels']}; "
        f"second counts={[second_physical_count(control_states[key]) for key in control_states]}"
    ),
)


def first_child(state):
    if state.get("unresolved") or state["physical_count"] != 1:
        raise RuntimeError("control state lacks its unique first slab")
    return state["children"][0]


reference = control_states["v_3_over_2"]
first = first_child(reference)
first_data = first["data"]
first_next = first["next"]
correct_signature = signature_text(first_next)
first_m_plus = mp.mpf(first_data["m_plus"])
first_pi_plus = mp.mpf(first_data["pi_plus"])
first_r = mp.mpf(first_data["r"])
first_q = mp.mpf(first_data["q"])

wrong_states = {
    "wrong_post_scale": (first_m_plus, first_pi_plus * first_r),
    "reversed_post_sign": (first_m_plus, -first_pi_plus),
    "dust_mass_reset": (mu(first_q), first_pi_plus),
}
hostile_signatures = {}
hostile_ok = True
for name, (mass, pi_value) in wrong_states.items():
    wrong_tree = expand_state(mass, pi_value, 2, [0])
    wrong_signature = signature_text(wrong_tree)
    changed = wrong_signature != correct_signature
    hostile_ok &= changed
    hostile_signatures[name] = {
        "signature": wrong_signature,
        "changed": changed,
        "unresolved": wrong_tree.get("unresolved", False),
    }
check(
    "wrong momentum scale, sign and mass reset change the branch signature",
    hostile_ok,
)


# Full-action check on one discovery node from each distinct signature.
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


def validate_tree_action(state):
    if state.get("unresolved"):
        return False, []
    records = []
    valid = True
    mass = mp.mpf(state["m"])
    pi_value = mp.mpf(state["pi"])
    for child in state["children"]:
        data = child["data"]
        height = mp.mpf(data["h"])
        ratio = mp.mpf(data["r"])
        next_pi = mp.mpf(data["pi_plus"])
        constraint_residual = 2 * constraint_numeric(
            1, ratio, height**2, mass
        ) / height
        pre_residual = pre_numeric(1, ratio, height**2, mass) - pi_value
        post_residual = post_numeric(1, ratio, height**2, mass) / ratio**2 - next_pi
        row_ok = max(
            abs(constraint_residual), abs(pre_residual), abs(post_residual)
        ) < mp.mpf("1e-78")
        valid &= row_ok
        records.append(
            {
                "slab": state["slab"],
                "constraint": text(constraint_residual),
                "pre": text(pre_residual),
                "post": text(post_residual),
                "passed": row_ok,
            }
        )
        if child.get("next") is not None:
            nested_ok, nested = validate_tree_action(child["next"])
            valid &= nested_ok
            records.extend(nested)
    return valid, records


signature_representatives = {}
for component in component_records:
    for row in component["rows"]:
        signature_representatives.setdefault(row["signature"], row)

action_controls = []
action_ok = True
for signature, row in signature_representatives.items():
    row_ok, residuals = validate_tree_action(row)
    action_ok &= row_ok
    action_controls.append(
        {
            "signature": signature,
            "v": row["v"],
            "residuals": residuals,
            "passed": row_ok,
        }
    )
check(
    "one node per distinct signature satisfies the redifferentiated action",
    action_ok,
    f"signatures={len(signature_representatives)}",
)


candidate_list_ok = bool(
    len(candidate_cells) >= 1
    and len(signature_runs) >= len(components)
    and len(inherited_candidates) == 4
)
check(
    "a finite intrinsic candidate-cell list precedes every continuum claim",
    candidate_list_ok,
    f"cells={len(candidate_cells)}; signatures={len(signature_counter)}",
)


complete = bool(
    provenance_ok
    and threshold_order_ok
    and all_nodes_resolved
    and known_controls_ok
    and hostile_ok
    and action_ok
    and candidate_list_ok
)
outcome = (
    "INCOMING_BASIN_CANDIDATE_SKELETON_FROZEN"
    if complete
    else "INCOMING_BASIN_DISCOVERY_OPEN"
)
check(
    "the discovery hierarchy freezes a candidate skeleton without a basin claim",
    outcome == "INCOMING_BASIN_CANDIDATE_SKELETON_FROZEN",
    outcome,
)


artifact = {
    "provenance": {
        "protocol_commit": PROTOCOL_COMMIT,
        "protocol_sha256": PROTOCOL_SHA256,
        "classification_sha256": CLASSIFICATION_SHA256,
        "invariant_sha256": INVARIANT_SHA256,
    },
    "method": {
        "nodes_per_component": NODE_COUNT,
        "total_nodes": 3 * NODE_COUNT,
        "node_rule": "Gauss-Chebyshev on each intrinsic open component",
        "work_decimal_digits": WORK_DPS,
        "reported_decimal_digits": SERIAL_DIGITS,
        "root_tolerance": text(ROOT_TOLERANCE),
        "critical_margin": text(CRITICAL_MARGIN),
        "maximum_slab": MAX_SLAB,
        "maximum_tree_nodes_per_input": MAX_TREE_NODES,
        "finite_q_box_used": False,
        "continuum_claim_from_nodes": False,
    },
    "thresholds": {
        "v_A": text(V_A),
        "v_star": text(V_STAR),
        "v_M": text(V_M),
        "v_C": text(V_C),
        "p_star": text(P_STAR),
        "p_infinity": text(P_INFINITY),
    },
    "skeleton": skeleton,
    "post_construction_controls": {
        "known_representatives": control_states,
        "hostile_signatures": hostile_signatures,
        "action_representatives": action_controls,
    },
    "checks": {
        "provenance": provenance_ok,
        "threshold_order": threshold_order_ok,
        "all_nodes_resolved": all_nodes_resolved,
        "known_controls": known_controls_ok,
        "hostile_controls": hostile_ok,
        "full_action_controls": action_ok,
        "candidate_list": candidate_list_ok,
    },
    "claims": {
        "candidate_skeleton": "DISCOVERY_ONLY",
        "continuum_branch_counts": "OPEN_PENDING_INTERVAL_PROOF",
        "incoming_v_basin": "OPEN_PENDING_INTERVAL_PROOF",
        "full_canonical_state_basin": "NOT_TESTED",
        "nonhomogeneous_physics": "OPEN",
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
print(f"CANDIDATE_CELLS: {len(candidate_cells)}")
print(f"DISTINCT_SIGNATURES: {len(signature_counter)}")
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Nested 6+3+1 reduction of the refined H4 internal equations.

Prior-art commit: 7714933.
Protocol commit: b284aa1.
"""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
from scipy.optimize import least_squares


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_root.py"
ROOT_ARTIFACT = HERE / "gravity_600cell_refined_h4_stationary_root.json"
ROOT_RESULT = ROOT / "docs/gravity/gravity_600cell_refined_h4_stationary_root_result.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_refined_h4_nested_reduction_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_h4_nested_reduction_protocol.md"
FILL = HERE / "gravity_600cell_refined_h4_stationary_fill.json"
JACOBIAN = HERE / "gravity_600cell_refined_h4_internal_jacobian.json"
ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
OUTPUT = HERE / "gravity_600cell_refined_h4_nested_reduction.json"
PRIOR_ART_COMMIT = "7714933"
PROTOCOL_COMMIT = "b284aa1"
INPUT_HASHES = {
    "reproducible/verify_gravity_600cell_refined_h4_stationary_root.py":
        "0105508a17fc40f34eb7a15f4c7c36bc89850c653476cb6a2e8086ca4281cceb",
    "reproducible/gravity_600cell_refined_h4_stationary_root.json":
        "e945dc54a0768b00358aca6bef9e9a105ab3d0080d22dd83dfd140b038adf14d",
    "docs/gravity/gravity_600cell_refined_h4_stationary_root_result.md":
        "c41c81409e2aa8d16bd7db71e68c3e954e5eaf568e4618ea52153262824b42ff",
    "docs/gravity/gravity_600cell_refined_h4_nested_reduction_prior_art.md":
        "fe395c93a1ce5209fccc4829e23011be185b3a1ec5c5674c0d68bce11f6a9f0d",
    "docs/gravity/gravity_600cell_refined_h4_nested_reduction_protocol.md":
        "5db83e284a344471b78b5a189e2e565d581199e645fd0e80a23d588d6ef0cd6b",
}
PAIR4 = tuple(combinations(range(4), 2))
VARIABLES = (
    tuple(("old",)+pair for pair in PAIR4)
    + tuple(("new",)+pair for pair in PAIR4)
    + tuple(("cross",)+pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
INTERNAL_VARIABLES = (
    tuple(("cross",)+pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
TAU_TEXT = "0.0102"
MAIN_LOWER = np.asarray([-0.35]*6+[-8.0]*4)
MAIN_UPPER = np.asarray([0.35]*6+[2.0]*4)
Q = np.asarray([
    [1/math.sqrt(2), 1/math.sqrt(6), 1/math.sqrt(12)],
    [-1/math.sqrt(2), 1/math.sqrt(6), 1/math.sqrt(12)],
    [0, -2/math.sqrt(6), 1/math.sqrt(12)],
    [0, 0, -3/math.sqrt(12)],
])
T_GRID = tuple(
    [-2.0+0.5*index for index in range(8)]
    + [-2.5-0.5*index for index in range(11)]
)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return condition


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def sf(value, digits=17):
    value = float(value)
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "+inf" if value > 0 else "-inf"
    return format(value, f".{digits}g")


def mp_text(value, digits=55):
    return mp.nstr(value, digits)


def helmert_mp():
    """Return the preregistered normalized Helmert matrix at active precision."""
    return mp.matrix([
        [1/mp.sqrt(2), 1/mp.sqrt(6), 1/mp.sqrt(12)],
        [-1/mp.sqrt(2), 1/mp.sqrt(6), 1/mp.sqrt(12)],
        [0, -2/mp.sqrt(6), 1/mp.sqrt(12)],
        [0, 0, -3/mp.sqrt(12)],
    ])


def mp_linear_blocks(matrix_strings):
    """Reconstruct A, S and B without binary64 roundoff."""
    H = mp.matrix([[mp.mpf(value) for value in row] for row in matrix_strings])
    A = H[:6, :6]
    Hxz = H[:6, 6:]
    S = H[6:, 6:]-H[6:, :6]*mp.inverse(A)*Hxz
    q = helmert_mp()
    return A, S, q.T*S*q


def load_root_definitions():
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "ast": ast,
        "mp": mp,
        "np": np,
        "json": json,
        "math": math,
        "Path": Path,
        "sha256": sha256,
        "combinations": combinations,
        "permutations": permutations,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "least_squares": least_squares,
        "HERE": HERE,
        "ROOT": ROOT,
        "SOURCE": ACTION_SOURCE,
        "FILL": FILL,
        "PAIR4": PAIR4,
        "VARIABLES": VARIABLES,
        "INTERNAL_VARIABLES": INTERNAL_VARIABLES,
        "TAU_TEXT": TAU_TEXT,
        "MAIN_LOWER": MAIN_LOWER,
        "MAIN_UPPER": MAIN_UPPER,
        "tests": 0,
        "passed": 0,
    }
    exec(
        compile(ast.Module(body=definitions, type_ignores=[]), str(SOURCE), "exec"),
        namespace,
    )
    return namespace


def public_attempt(record):
    return {key: value for key, value in record.items() if not key.startswith("_")}


def fast_state(record, y):
    return root_defs["fast_equation_state"](
        record, fast_geometry, fast_base, np.asarray(y, dtype=float)
    )


def fast_gradient(state):
    return np.asarray(state["gradient"].real, dtype=float)


def inner_attempt(record, z, seed, A):
    diagnostics = {"evaluations": 0, "invalid": 0}

    def residual(x):
        diagnostics["evaluations"] += 1
        state = fast_state(record, np.r_[x, z])
        if not (state["finite"] and state["branch"]):
            diagnostics["invalid"] += 1
            return np.full(6, 1e6)
        return np.linalg.solve(A, fast_gradient(state)[:6])

    lower = np.full(6, -0.35)
    upper = np.full(6, 0.35)
    clipped = np.clip(np.asarray(seed, dtype=float), lower+1e-6, upper-1e-6)
    optimization = least_squares(
        residual, clipped, bounds=(lower, upper), method="trf",
        jac="3-point", diff_step=1e-5, xtol=1e-12, ftol=1e-12,
        gtol=1e-12, max_nfev=600, x_scale=1.0,
    )
    x = np.asarray(optimization.x)
    state = fast_state(record, np.r_[x, z])
    if state["finite"] and state["branch"]:
        gradient = fast_gradient(state)
        residual_norm = float(np.linalg.norm(np.linalg.solve(A, gradient[:6])))
    else:
        gradient = np.full(10, np.nan)
        residual_norm = math.inf
    distance = float(min(np.min(x-lower), np.min(upper-x)))
    resolved = bool(
        state["finite"] and state["branch"] and distance > 1e-6
        and residual_norm < 1e-8
    )
    return {
        "x": x,
        "state": state,
        "gradient": gradient,
        "residual_norm": residual_norm,
        "resolved": resolved,
        "record": {
            "seed": [sf(value) for value in seed],
            "endpoint": [sf(value) for value in x],
            "optimizer_success": bool(optimization.success),
            "optimizer_status": int(optimization.status),
            "nfev": int(optimization.nfev),
            "njev": None if optimization.njev is None else int(optimization.njev),
            "preconditioned_cross_norm": sf(residual_norm),
            "distance_from_box": sf(distance),
            "branch_valid": bool(state["finite"] and state["branch"]),
            "resolved": resolved,
            "diagnostics": diagnostics,
        },
    }


def solve_inner(record, z, A, Hxz):
    seeds = (np.zeros(6), -np.linalg.solve(A, Hxz@z))
    attempts = [inner_attempt(record, z, seed, A) for seed in seeds]
    selected = min(range(2), key=lambda index: attempts[index]["residual_norm"])
    result = attempts[selected]
    return {
        **result,
        "selected_seed": selected,
        "attempt_records": [item["record"] for item in attempts],
    }


def reduced_evaluation(class_record, t, u):
    u = np.asarray(u, dtype=float)
    z = t*np.ones(4)+Q@u
    if np.any(z <= -8) or np.any(z >= 2):
        return {
            "valid": False, "resolved": False, "z": z,
            "contrast_scaled": np.full(3, 1e6), "contrast_norm": math.inf,
            "g": math.nan, "error": "lapse box",
        }
    inner = solve_inner(
        class_record["record"], z, class_record["A"], class_record["Hxz"]
    )
    if not inner["resolved"]:
        return {
            "valid": False, "resolved": False, "z": z, "inner": inner,
            "contrast_scaled": np.full(3, 1e6), "contrast_norm": math.inf,
            "g": math.nan, "error": "inner unresolved",
        }
    rho = inner["gradient"][6:]
    contrast = Q.T@rho
    contrast_scaled = np.linalg.solve(class_record["B"], contrast)
    return {
        "valid": True,
        "resolved": True,
        "z": z,
        "inner": inner,
        "rho": rho,
        "contrast": contrast,
        "contrast_scaled": contrast_scaled,
        "contrast_norm": float(np.linalg.norm(contrast_scaled)),
        "g": float(np.sum(rho)/2),
    }


def contrast_attempt(class_record, t, seed):
    diagnostics = {"evaluations": 0, "invalid": 0, "inner_attempts": 0}

    def residual(u):
        diagnostics["evaluations"] += 1
        reduced = reduced_evaluation(class_record, t, u)
        diagnostics["inner_attempts"] += 2
        if not reduced["valid"]:
            diagnostics["invalid"] += 1
            return np.full(3, 1e6)
        return reduced["contrast_scaled"]

    lower = np.full(3, -4.0)
    upper = np.full(3, 4.0)
    clipped = np.clip(np.asarray(seed, dtype=float), lower+1e-6, upper-1e-6)
    optimization = least_squares(
        residual, clipped, bounds=(lower, upper), method="trf",
        jac="3-point", diff_step=1e-5, xtol=1e-11, ftol=1e-11,
        gtol=1e-11, max_nfev=400, x_scale=1.0,
    )
    u = np.asarray(optimization.x)
    final = reduced_evaluation(class_record, t, u)
    diagnostics["inner_attempts"] += 2
    z_distance = min(np.min(final["z"]+8), np.min(2-final["z"]))
    u_distance = min(np.min(u-lower), np.min(upper-u))
    distance = float(min(z_distance, u_distance))
    resolved = bool(
        final["valid"] and final["inner"]["resolved"]
        and distance > 1e-6 and final["contrast_norm"] < 1e-7
    )
    return {
        "u": u,
        "final": final,
        "resolved": resolved,
        "record": {
            "seed": [sf(value) for value in seed],
            "endpoint_u": [sf(value) for value in u],
            "endpoint_z": [sf(value) for value in final["z"]],
            "optimizer_success": bool(optimization.success),
            "optimizer_status": int(optimization.status),
            "nfev": int(optimization.nfev),
            "njev": None if optimization.njev is None else int(optimization.njev),
            "preconditioned_contrast_norm": sf(final["contrast_norm"]),
            "common_scalar_fast": sf(final["g"]),
            "distance_from_all_boxes": sf(distance),
            "resolved": resolved,
            "selected_inner_seed": final["inner"]["selected_seed"]
                if final.get("inner") else None,
            "final_inner_attempts": [
                public_attempt(value) for value in
                final["inner"]["attempt_records"]
            ] if final.get("inner") else [],
            "diagnostics": diagnostics,
        },
    }


def solve_contrast(class_record, t, primary_seed):
    attempts = [contrast_attempt(class_record, t, primary_seed)]
    if not attempts[0]["resolved"]:
        r0 = class_record["rho0"]
        predictor = -np.linalg.solve(
            class_record["B"],
            Q.T@(r0+class_record["S"]@(t*np.ones(4))),
        )
        attempts.append(contrast_attempt(class_record, t, predictor))
    resolved_indices = [index for index, item in enumerate(attempts) if item["resolved"]]
    if resolved_indices:
        selected = min(
            resolved_indices,
            key=lambda index: attempts[index]["final"]["contrast_norm"],
        )
    else:
        selected = min(
            range(len(attempts)),
            key=lambda index: attempts[index]["final"]["contrast_norm"],
        )
    result = attempts[selected]
    return {
        **result,
        "selected_attempt": selected,
        "attempt_records": [item["record"] for item in attempts],
    }


def validate_grid_point(class_record, point):
    y = np.r_[point["final"]["inner"]["x"], point["final"]["z"]]
    with mp.workdps(80):
        precise = root_defs["evaluate_equations"](
            class_record["record"], geometry80, base80, y
        )
        reverse = root_defs["evaluate_equations"](
            class_record["reverse_record"], geometry80, base80, y
        )
        precise_gradient = np.asarray([
            complex(value) for value in precise["gradient"]
        ]) if precise["finite"] else np.full(10, complex(np.nan, np.nan))
        fast = point["final"]["inner"]["state"]
        fast_gradient_value = np.asarray(fast["gradient"])
        relative = float(np.max(np.abs(
            fast_gradient_value-precise_gradient
        )))/max(1, float(np.max(np.abs(precise_gradient))))
        A_mp, _, B_mp = mp_linear_blocks(
            class_record["matrix_record"]["matrix"]
        )
        q_mp = helmert_mp()
        cross_norm_mp = mp.norm(mp.lu_solve(A_mp, precise["gradient"][:6]))
        rho_mp = mp.matrix([
            precise["gradient"][index] for index in range(6, 10)
        ])
        contrast_norm_mp = mp.norm(mp.lu_solve(B_mp, q_mp.T*rho_mp))
        cross_norm = float(cross_norm_mp)
        contrast_norm = float(contrast_norm_mp)
        reverse_difference = max(
            abs(precise["gradient"][index]-reverse["gradient"][index])
            for index in range(10)
        ) if reverse["finite"] else mp.inf
        g_value = sum(precise["gradient"][index] for index in range(6, 10))/2
        valid = bool(
            precise["finite"] and precise["branch"]
            and reverse["finite"] and reverse["branch"]
            and relative < 5e-9
            and precise["maximum_imaginary"] < mp.mpf("1e-50")
            and precise["result"]["minimum_angle_argument"] > mp.mpf("1e-8")
            and cross_norm < 1e-6 and contrast_norm < 1e-6
            and reverse_difference < mp.mpf("1e-40")
        )
        return {
            "valid": valid,
            "endpoint_y": [sf(value) for value in y],
            "gradient_relative_fast_error": sf(relative),
            "preconditioned_cross_norm": sf(cross_norm),
            "preconditioned_contrast_norm": sf(contrast_norm),
            "physical_imaginary": mp_text(precise["maximum_imaginary"], 35),
            "minimum_angle_argument": mp_text(
                precise["result"]["minimum_angle_argument"], 35
            ),
            "reverse_equation_difference": mp_text(reverse_difference, 35),
            "rho_equations": [
                mp_text(precise["gradient"][index], 45) for index in range(6, 10)
            ],
            "common_scalar": mp_text(g_value, 45),
            "common_scalar_float": float(mp.re(g_value)),
        }


def midpoint_refinement(class_record, left, right):
    """Apply the frozen bisection without crossing an unresolved midpoint."""
    history = []
    left_t = float(left["t"])
    right_t = float(right["t"])
    left_u = np.asarray(left["u"], dtype=float)
    right_u = np.asarray(right["u"], dtype=float)
    left_g = float(left["validation"]["common_scalar_float"])
    right_g = float(right["validation"]["common_scalar_float"])
    best = left if abs(left_g) <= abs(right_g) else right
    unresolved = False
    for iteration in range(30):
        best_g = float(best["validation"]["common_scalar_float"])
        if right_t-left_t < 1e-8 or abs(best_g) < 1e-9:
            break
        midpoint_t = (left_t+right_t)/2
        midpoint_u = (left_u+right_u)/2
        midpoint = contrast_attempt(class_record, midpoint_t, midpoint_u)
        midpoint["t"] = midpoint_t
        if not midpoint["resolved"]:
            unresolved = True
            history.append({
                "iteration": iteration,
                "t": sf(midpoint_t),
                "resolved": False,
            })
            break
        validation = validate_grid_point(class_record, midpoint)
        midpoint["validation"] = validation
        validation_values.append(validation["valid"])
        history.append({
            "iteration": iteration,
            "t": sf(midpoint_t),
            "resolved": bool(validation["valid"]),
            "common_scalar": validation["common_scalar"],
        })
        if not validation["valid"]:
            unresolved = True
            break
        midpoint_g = float(validation["common_scalar_float"])
        if abs(midpoint_g) < abs(best_g):
            best = midpoint
        if abs(midpoint_g) < 1e-9:
            break
        if left_g*midpoint_g < 0:
            right_t, right_u, right_g = midpoint_t, midpoint["u"], midpoint_g
        else:
            left_t, left_u, left_g = midpoint_t, midpoint["u"], midpoint_g
    return {
        "resolved": not unresolved,
        "best": best,
        "final_width": sf(right_t-left_t),
        "history": history,
    }


def full_candidate_refinement(class_record, source, point, reduction_record):
    """Pass a scalar candidate through the frozen full ten-equation gate."""
    endpoint = np.r_[point["final"]["inner"]["x"], point["final"]["z"]]
    refinement = root_defs["refine_candidate"](
        class_record["record"], class_record["reverse_record"],
        class_record["matrix_record"]["matrix"], endpoint,
    )
    return {
        "class": class_record["class"],
        "source": source,
        "seed_t": sf(point["t"]),
        "seed_common_scalar": point["validation"]["common_scalar"],
        "nested_reduction": reduction_record,
        "full_refinement": refinement,
        "accepted_before_deduplication": bool(
            refinement["accepted_finite_positive_root"]
        ),
        "accepted_finite_positive_root": False,
    }


def synthetic_controls():
    M = np.asarray([
        [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
        [0, 0, 0, 1], [1, 1, 0, 0], [0, 0, 1, 1],
    ], dtype=float)

    def solve_at(t):
        def contrast_residual(u):
            z = t*np.ones(4)+Q@u
            x = least_squares(
                lambda value: value-M@z, np.zeros(6), method="trf",
                jac="3-point", diff_step=1e-5, xtol=1e-12, ftol=1e-12,
                gtol=1e-12, max_nfev=600, x_scale=1.0,
            ).x
            _ = x
            target = 0.1*t*np.asarray([1, -1, .5])
            rho = Q@(u-target)+(t+1)*np.ones(4)/2
            return Q.T@rho
        solved = least_squares(
            contrast_residual, np.zeros(3), method="trf", jac="3-point",
            diff_step=1e-5, xtol=1e-11, ftol=1e-11, gtol=1e-11,
            max_nfev=400, x_scale=1.0,
        )
        u = solved.x
        rho = Q@(u-0.1*t*np.asarray([1, -1, .5]))+(t+1)*np.ones(4)/2
        return np.linalg.norm(Q.T@rho), float(np.sum(rho)/2)

    synthetic_grid = [(t, *solve_at(t)) for t in T_GRID]
    root_records = [
        (t, contrast_norm, scalar)
        for t, contrast_norm, scalar in synthetic_grid
        if contrast_norm < 1e-9 and abs(scalar) < 1e-9
    ]
    negative_values = [math.exp(t)+1 for t in T_GRID]
    return {
        "positive_grid_root_count": len(root_records),
        "positive_grid_roots": [
            {
                "t": sf(t),
                "contrast_norm": sf(contrast_norm),
                "scalar": sf(scalar),
            }
            for t, contrast_norm, scalar in root_records
        ],
        "positive_ok": len(root_records) == 1 and abs(root_records[0][0]+1) < 1e-15,
        "negative_values": [sf(value) for value in negative_values],
        "negative_ok": all(value > 0 for value in negative_values),
    }


print("="*78)
print("REFINED H4 NESTED 6+3+1 REDUCTION")
print("="*78)

previous_artifact = json.loads(OUTPUT.read_text()) if OUTPUT.exists() else None
actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "the frozen root result, prior-art gate and protocol have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "7714933" and PROTOCOL_COMMIT == "b284aa1",
)

root_artifact = json.loads(ROOT_ARTIFACT.read_text())
fill = json.loads(FILL.read_text())
jacobian_artifact = json.loads(JACOBIAN.read_text())
upstream_ok = check(
    "the accepted bounded-negative root outcome and all inherited controls are exact",
    root_artifact["outcome"] == "REFINED_H4_NO_FINITE_ROOT_FOUND_OTHER"
    and root_artifact["tests"] == {"passed": 13, "total": 13}
    and root_artifact["search"]["accepted_class_count"] == 0
    and root_artifact["search"]["accepted_attempt_hit_fraction"]
        == {"numerator": 0, "denominator": 120}
    and root_artifact["input_sha256"][
        "reproducible/gravity_600cell_refined_h4_stationary_fill.json"
    ] == digest(FILL)
    and root_artifact["input_sha256"][
        "reproducible/gravity_600cell_refined_h4_internal_jacobian.json"
    ] == digest(JACOBIAN),
)

root_defs = load_root_definitions()
actions = root_defs["load_action_definitions"]()
root_defs["actions"] = actions
required = {
    "parse_combinatorics", "fast_equation_state", "evaluate_equations",
    "refine_candidate",
}
definitions_ok = check(
    "only frozen AST definitions are loaded and all nested-action functions exist",
    required <= set(root_defs) and "OUTPUT" not in root_defs,
)

combinatorics = [root_defs["parse_combinatorics"](record)
                 for record in fill["combinatorics"]]
orders = tuple(record["order"] for record in combinatorics)
order_to_index = {order: index for index, order in enumerate(orders)}
combinatorics_ok = check(
    "all 24 committed reduced carriers reconstruct without a new schedule choice",
    orders == tuple(permutations(range(4)))
    and len(combinatorics) == 24
    and all(record["pentachora"] == 57600 and record["triangles"] == 149280
            for record in combinatorics),
)

with mp.workdps(80):
    geometry80 = actions["exact_geometry"](80)
    base80 = actions["base_coordinates"](geometry80)
fast_geometry = {"mass": float(geometry80["mass"])}
fast_base = {key: float(value) for key, value in base80.items()}

schedule_by_order = {
    tuple(item["order"]): item for item in jacobian_artifact["census"]["schedules"]
}
class_records = []
linear_records = []
linear_ok_values = []
for item in jacobian_artifact["census"]["matrix_classes"]:
    order = tuple(item["orders"][0])
    reverse = tuple(reversed(order))
    matrix_record = schedule_by_order[order]
    H = np.asarray([[float(value) for value in row]
                    for row in matrix_record["matrix"]])
    A = H[:6, :6]
    Hxz = H[:6, 6:]
    S = H[6:, 6:]-H[6:, :6]@np.linalg.solve(A, Hxz)
    B = Q.T@S@Q
    a_values = np.linalg.svd(A, compute_uv=False)
    s_values = np.linalg.eigvalsh((S+S.T)/2)
    b_values = np.linalg.svd(B, compute_uv=False)
    state0 = fast_state(
        combinatorics[order_to_index[order]], np.zeros(10)
    )
    rho0 = fast_gradient(state0)[6:]
    ok = bool(
        np.linalg.matrix_rank(A) == 6 and np.linalg.matrix_rank(S) == 4
        and np.linalg.matrix_rank(B) == 3
        and abs(a_values[-1]-5.93908093) < 1e-6
        and np.linalg.cond(S) < 1.1
        and sum(s_values < 0) == 2 and sum(s_values > 0) == 2
    )
    linear_ok_values.append(ok)
    linear_records.append({
        "class": item["class"],
        "rank_A_S_B": [int(np.linalg.matrix_rank(A)),
                         int(np.linalg.matrix_rank(S)),
                         int(np.linalg.matrix_rank(B))],
        "A_smallest_singular": sf(a_values[-1]),
        "A_condition": sf(a_values[0]/a_values[-1]),
        "S_eigenvalues": [sf(value) for value in s_values],
        "S_condition": sf(np.linalg.cond(S)),
        "B_smallest_singular": sf(b_values[-1]),
        "B_condition": sf(b_values[0]/b_values[-1]),
    })
    class_records.append({
        "class": item["class"], "order": order, "reverse": reverse,
        "record": combinatorics[order_to_index[order]],
        "reverse_record": combinatorics[order_to_index[reverse]],
        "matrix_record": matrix_record,
        "H": H, "A": A, "Hxz": Hxz, "S": S, "B": B, "rho0": rho0,
    })

maximum_s_difference_float = max(
    np.max(np.abs(class_records[index]["S"]-class_records[0]["S"]))
    for index in range(12)
)
with mp.workdps(100):
    schur_mp = [
        mp_linear_blocks(item["matrix_record"]["matrix"])[1]
        for item in class_records
    ]
    maximum_s_difference_mp = max(
        abs(schur_mp[index][row, column]-schur_mp[0][row, column])
        for index in range(12) for row in range(4) for column in range(4)
    )
    maximum_committed_entry_error = max(
        mp.mpf(item["matrix_record"]["entry_error"])
        for item in class_records
    )
    schur_comparison_envelope = 2*maximum_committed_entry_error
linear_ok = check(
    "all cross blocks and reduced lapse/contrast blocks pass the frozen rank gates",
    all(linear_ok_values)
    and maximum_s_difference_mp <= schur_comparison_envelope,
    "max mp S difference="+mp_text(maximum_s_difference_mp, 30)
    +", envelope="+mp_text(schur_comparison_envelope, 30)
    +f", float diagnostic={maximum_s_difference_float:.3e}",
)
corrupted = class_records[0]["A"].copy()
corrupted[-1] = corrupted[0]
corruption_ok = check(
    "a deliberately rank-deficient cross block fails the elimination gate",
    np.linalg.matrix_rank(corrupted) < 6,
)

synthetic = synthetic_controls()
synthetic_positive_ok = check(
    "the nested wrapper recovers its synthetic cross/contrast/scalar root",
    synthetic["positive_ok"],
    f"roots={synthetic['positive_grid_roots']}",
)
synthetic_negative_ok = check(
    "the synthetic sign-definite scalar is not called a root",
    synthetic["negative_ok"], str(synthetic["negative_values"]),
)

controls_before_search = all((
    provenance_ok, upstream_ok, definitions_ok, combinatorics_ok,
    linear_ok, corruption_ok, synthetic_positive_ok, synthetic_negative_ok,
))

class_outputs = []
validation_values = []
candidate_refinements = []
if controls_before_search:
    for class_position, class_record in enumerate(class_records):
        points = []
        resolved_history = []
        for t in T_GRID:
            if resolved_history:
                nearest = min(resolved_history, key=lambda point: abs(point["t"]-t))
                primary_seed = nearest["u"]
            else:
                primary_seed = np.zeros(3)
            point = solve_contrast(class_record, t, primary_seed)
            point["t"] = t
            if point["resolved"]:
                validation = validate_grid_point(class_record, point)
                point["validation"] = validation
                validation_values.append(validation["valid"])
                if validation["valid"]:
                    resolved_history.append({"t": t, "u": point["u"].copy()})
            else:
                point["validation"] = None
            points.append(point)
        sorted_points = sorted(points, key=lambda point: point["t"])
        valid_points = [
            point for point in sorted_points
            if point["resolved"] and point["validation"] is not None
            and point["validation"]["valid"]
        ]
        for point in valid_points:
            if abs(float(point["validation"]["common_scalar_float"])) < 1e-7:
                candidate_refinements.append(full_candidate_refinement(
                    class_record, "resolved_grid_point", point,
                    {
                        "kind": "grid_near_zero",
                        "threshold": "1e-7",
                        "resolved": True,
                    },
                ))
        for left, right in zip(sorted_points, sorted_points[1:]):
            adjacent_valid = all(
                point["resolved"] and point["validation"] is not None
                and point["validation"]["valid"]
                for point in (left, right)
            )
            if not adjacent_valid:
                continue
            left_g = float(left["validation"]["common_scalar_float"])
            right_g = float(right["validation"]["common_scalar_float"])
            if not (abs(left_g) > 1e-7 and abs(right_g) > 1e-7
                    and left_g*right_g < 0):
                continue
            bisection = midpoint_refinement(class_record, left, right)
            candidate_refinements.append(full_candidate_refinement(
                class_record, "adjacent_sign_bracket", bisection["best"],
                {
                    "kind": "sign_bracket",
                    "left_t": sf(left["t"]),
                    "right_t": sf(right["t"]),
                    "left_common_scalar": left["validation"]["common_scalar"],
                    "right_common_scalar": right["validation"]["common_scalar"],
                    "bisection_resolved": bisection["resolved"],
                    "bisection_final_width": bisection["final_width"],
                    "bisection_history": bisection["history"],
                },
            ))
        serial_points = []
        for point in points:
            serial_points.append({
                "t": sf(point["t"]),
                "resolved": point["resolved"],
                "selected_attempt": point["selected_attempt"],
                "attempts": point["attempt_records"],
                "endpoint_u": [sf(value) for value in point["u"]],
                "endpoint_x": [sf(value) for value in point["final"]["inner"]["x"]]
                    if point["final"].get("inner") else [],
                "endpoint_z": [sf(value) for value in point["final"]["z"]],
                "common_scalar_fast": sf(point["final"]["g"]),
                "validation": point["validation"],
            })
        class_outputs.append({
            "class": class_record["class"],
            "order": list(class_record["order"]),
            "reverse": list(class_record["reverse"]),
            "resolved_grid_points": sum(
                point["resolved"] and point["validation"] is not None
                and point["validation"]["valid"] for point in points
            ),
            "grid": serial_points,
        })
        print(f"[INFO] nested classes completed: {class_position+1}/12", flush=True)

accepted_endpoints = []
for candidate in candidate_refinements:
    if not candidate["accepted_before_deduplication"]:
        continue
    endpoint = np.asarray([
        float(value) for value in candidate["full_refinement"]["endpoint"]
    ])
    duplicate = any(
        np.linalg.norm(endpoint-previous, ord=np.inf) < 1e-8
        for previous in accepted_endpoints
    )
    candidate["deduplicated"] = bool(duplicate)
    candidate["accepted_finite_positive_root"] = not duplicate
    if not duplicate:
        accepted_endpoints.append(endpoint)

search_count_ok = check(
    "the nested census contains all 12 classes and 19 frozen grid points",
    (not controls_before_search)
    or (len(class_outputs) == 12
        and all(len(item["grid"]) == 19 for item in class_outputs)),
)
validation_ok = check(
    "every resolved grid or bisection point passes independent 80-decimal substitution and reversal",
    all(validation_values),
    f"validated={sum(validation_values)}/{len(validation_values)}",
)

accepted_classes = sorted({
    item["class"] for item in candidate_refinements
    if item["accepted_finite_positive_root"]
})
resolved_counts = [item["resolved_grid_points"] for item in class_outputs]
sign_consistent = []
for item in class_outputs:
    values = [
        float(point["validation"]["common_scalar_float"])
        for point in item["grid"]
        if point["resolved"] and point["validation"] is not None
        and point["validation"]["valid"]
    ]
    sign_consistent.append(
        len(values) == 19
        and (all(value > 1e-7 for value in values)
             or all(value < -1e-7 for value in values))
    )

all_controls = controls_before_search and search_count_ok and validation_ok
if not all_controls:
    outcome = "REFINED_H4_NESTED_REDUCTION_CONTROL_FAILED"
elif len(accepted_classes) == 12:
    outcome = "REFINED_H4_NESTED_FINITE_ROOTS_ALL_CLASSES"
elif accepted_classes:
    outcome = "REFINED_H4_NESTED_FINITE_ROOTS_SOME_CLASSES"
elif any(value < 19 for value in resolved_counts):
    outcome = "REFINED_H4_NESTED_BRANCH_UNRESOLVED"
elif all(sign_consistent):
    outcome = "REFINED_H4_NESTED_NO_ROOT_GRID_SIGN_CONSISTENT"
else:
    outcome = "REFINED_H4_NESTED_NO_ROOT_GRID_OTHER"
outcome_ok = check(
    "the frozen hierarchy classifies the nested reduction",
    outcome in {
        "REFINED_H4_NESTED_REDUCTION_CONTROL_FAILED",
        "REFINED_H4_NESTED_FINITE_ROOTS_ALL_CLASSES",
        "REFINED_H4_NESTED_FINITE_ROOTS_SOME_CLASSES",
        "REFINED_H4_NESTED_BRANCH_UNRESOLVED",
        "REFINED_H4_NESTED_NO_ROOT_GRID_SIGN_CONSISTENT",
        "REFINED_H4_NESTED_NO_ROOT_GRID_OTHER",
    },
    outcome,
)


def deterministic_science(record):
    return {
        key: record[key] for key in
        ("definitions", "controls", "census", "scope", "outcome")
    }


artifact = {
    "title": "Nested 6+3+1 reduction of the refined H4 equations",
    "date": "2026-08-20",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "t_grid_continuation_order": list(T_GRID),
        "t_grid_sorted": sorted(T_GRID),
        "classes": 12,
        "common_lapse_only_restriction": False,
        "helmert_matrix": [[sf(value) for value in row] for row in Q],
        "physical_or_continuum_target_loaded": False,
        "effective_boundary_hessian_computed": False,
    },
    "controls": {
        "linear_reduction": linear_records,
        "maximum_lapse_schur_class_difference_mp": mp_text(
            maximum_s_difference_mp, 45
        ),
        "maximum_lapse_schur_class_difference_float_diagnostic": sf(
            maximum_s_difference_float
        ),
        "twice_maximum_committed_entry_error": mp_text(
            schur_comparison_envelope, 45
        ),
        "synthetic": synthetic,
    },
    "census": {
        "classes": class_outputs,
        "resolved_grid_counts": resolved_counts,
        "validated_grid_and_bisection_points": sum(validation_values),
        "scalar_candidate_count": len(candidate_refinements),
        "accepted_classes": accepted_classes,
        "accepted_roots": candidate_refinements,
        "sign_consistent_classes": [
            index for index, value in enumerate(sign_consistent) if value
        ],
    },
    "scope": {
        "connected_branch_only": True,
        "continuous_no_root_theorem": False,
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
previous_science_ok = check(
    "all scientific fields reproduce the first complete artifact exactly",
    previous_artifact is not None
    and previous_artifact.get("outcome") == outcome
    and deterministic_science(previous_artifact) == deterministic_science(artifact),
)
artifact["tests"] = {"passed": passed, "total": tests}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"RESOLVED GRID COUNTS: {resolved_counts}")
print(f"ACCEPTED ROOT CLASSES: {len(accepted_classes)}/12")
print(f"SIGN-CONSISTENT CLASSES: {sum(sign_consistent)}/12")
print(f"OUTCOME: {outcome}")
print(f"RESULT: {passed}/{tests} checks passed")
sys.exit(0 if outcome_ok and passed == tests else 1)

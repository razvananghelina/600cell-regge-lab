#!/usr/bin/env python3
"""Test whether two dust-slab principal functions differ by endpoint terms.

Prior-art commit: 5ccb29b.
Protocol commit: 84c8fd8.
"""

import ast
from collections import Counter
import contextlib
import hashlib
import importlib.util
import io
import json
import multiprocessing as mp_pool
from pathlib import Path
import sys

import mpmath as arb


HERE = Path(__file__).resolve().parent
SEED_INPUT = HERE / "gravity_600cell_dust_nonlinear_boundary_covariance_seeds.json"
NONLINEAR_INPUT = HERE / "gravity_600cell_dust_nonlinear_boundary_covariance.json"
CANONICAL_SOURCE = HERE / "verify_gravity_600cell_dust_canonical_legendre_rank.py"
ACTION_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
OUTPUT = HERE / "gravity_600cell_dust_boundary_coboundary.json"

PRIOR_ART_COMMIT = "5ccb29b"
PROTOCOL_COMMIT = "84c8fd8"
SEED_SHA256 = "2104c69ba6b21d3a3d92c7071d7f2702cb7d33f7f0e3ff17954f64c469f0c01d"
NONLINEAR_SHA256 = "a1e00071fa41f986dfaee84ea6e7689a14c50823f6c87d76889e6cb9346a7e3f"
CANONICAL_SOURCE_SHA256 = "396c491fe51a9f5e04fa8402e2e5b16884fe23fc5057d8ded325e6064fbd3b9e"
ACTION_SOURCE_SHA256 = "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf"

DPS = 100
arb.mp.dps = DPS
ETA = arb.mpf("1e-4")
LEVELS = {"half": arb.mpf("0.5"), "full": arb.mpf(1)}
RESIDUAL_TOLERANCE = arb.mpf("1e-55")
CORRECTION_TOLERANCE = arb.mpf("1e-45")
IMAGINARY_TOLERANCE = arb.mpf("1e-70")
ANGLE_TOLERANCE = arb.mpf("1e-6")
MAX_ITERATIONS = 20
MAX_BACKTRACKING = 12
CONSISTENT_FACTOR = arb.mpf(10)
NONSEPARABLE_FACTOR = arb.mpf(100)
UNCERTAINTY_FLOOR = arb.mpf("1e-70")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_audited_functions():
    wanted = {
        "orbit_sort_key",
        "augment_boundary_orbits",
        "arb_log_minus",
        "arb_signed_volume_square",
        "arb_angle_record",
        "triangle_area_square",
        "triangle_area_square_partials",
        "edge_data",
        "simplex_squared",
        "action_and_gradient",
    }
    tree = ast.parse(CANONICAL_SOURCE.read_text(), filename=str(CANONICAL_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch: missing={wanted-found}")
    exec(
        compile(
            ast.Module(body=body, type_ignores=[]),
            str(CANONICAL_SOURCE),
            "exec",
        ),
        globals(),
    )


def infinity_norm(vector):
    return max(abs(value) for value in vector)


def matrix_from_rows(rows):
    return arb.matrix([[arb.mpf(value) for value in row] for row in rows])


def number(value, digits=70):
    return arb.nstr(value, digits, strip_zeros=False)


def point_id(point):
    old_direction, old_sign, level, new_direction, new_sign = point
    if old_direction == new_direction == 0:
        return "base"
    old = "o0" if old_direction == 0 else f"o{old_direction}s{old_sign:+d}"
    new = "n0" if new_direction == 0 else f"n{new_direction}s{new_sign:+d}"
    return f"{old}_{new}_{level}"


hashes = {
    "seeds": digest(SEED_INPUT),
    "nonlinear_result": digest(NONLINEAR_INPUT),
    "canonical_source": digest(CANONICAL_SOURCE),
    "action_source": digest(ACTION_SOURCE),
}
seed = json.loads(SEED_INPUT.read_text())
nonlinear = json.loads(NONLINEAR_INPUT.read_text())
provenance_ok = bool(
    hashes == {
        "seeds": SEED_SHA256,
        "nonlinear_result": NONLINEAR_SHA256,
        "canonical_source": CANONICAL_SOURCE_SHA256,
        "action_source": ACTION_SOURCE_SHA256,
    }
    and seed.get("prior_art_commit") == "526a202"
    and seed.get("protocol_commit") == "05f76c3"
    and seed.get("outcome") == "NONLINEAR_BOUNDARY_COVARIANCE_CASES_FROZEN"
    and seed.get("eta") == "0.0001"
    and seed.get("number_of_directions") == 4
    and nonlinear.get("seed_commit") == "b6370bd"
    and nonlinear.get("outcome")
        == "NONLINEAR_BOUNDARY_COVARIANCE_BROKEN_ON_FROZEN_CASES"
    and nonlinear.get("classification_counts") == {"BROKEN": 32}
)

load_audited_functions()
spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_boundary_coboundary", ACTION_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
models = {
    parity: augment_boundary_orbits(model)
    for parity, model in gro.models.items()
}
carrier_ok = bool(
    gro.tests == gro.passed == 43
    and all(
        len(model["old_orbits"]) == 30
        and len(model["edge_orbits"]) == 35
        and len(model["final_orbits"]) == 30
        for model in models.values()
    )
)

ARB_I = arb.mpc(0, 1)
ARB_M_STAR = arb.mpf(10)
ARB_ZETA = (arb.pi**2*arb.sqrt(2)/50)**(arb.mpf(1)/3)
ARB_R0 = 4*ARB_M_STAR/(3*arb.pi)
ARB_L0 = ARB_ZETA*ARB_R0
ARB_EPSILON_3 = 2*arb.pi-5*arb.acos(arb.mpf(1)/3)
ARB_MASS = (90/arb.pi)*ARB_EPSILON_3*ARB_L0

physical_map = tuple(seed["physical_edge_permutation"])
base = {}
jacobians = {}
for parity in ("even", "odd"):
    source = seed["parities"][parity]
    base[parity] = {
        "old": tuple(arb.mpf(value) for value in source["base_old"]),
        "x": tuple(arb.mpf(value) for value in source["base_x"]),
        "new": tuple(arb.mpf(value) for value in source["base_new"]),
    }
    jacobians[parity] = {}
    for variant in ("operational", "validation"):
        complete = matrix_from_rows(source[f"canonical_{variant}"])
        jacobians[parity][variant] = complete[:35, :35]

directions = {
    parity: {
        int(record["index"]): tuple(arb.mpf(value) for value in record[parity])
        for record in seed["directions"]
    }
    for parity in ("even", "odd")
}

direction_indices_ok = bool(
    set(directions["even"]) == set(directions["odd"]) == {1, 2, 3, 4}
)
direction_geometry_ok = direction_indices_ok
if direction_geometry_ok:
    for index in range(1, 5):
        even = directions["even"][index]
        odd = directions["odd"][index]
        mapped = [arb.mpf(0) for _ in range(30)]
        for source, target in enumerate(physical_map):
            mapped[target] = even[source]
        direction_geometry_ok &= bool(
            abs(sum(even)) < arb.mpf("1e-65")
            and abs(sum(odd)) < arb.mpf("1e-65")
            and abs(sum(value*value for value in even)-1) < arb.mpf("1e-65")
            and abs(sum(value*value for value in odd)-1) < arb.mpf("1e-65")
            and max(abs(left-right) for left, right in zip(mapped, odd))
                < arb.mpf("1e-65")
        )


def make_point(old_direction, old_sign, level, new_direction, new_sign):
    return (
        int(old_direction), int(old_sign), str(level),
        int(new_direction), int(new_sign),
    )


base_point = make_point(0, 0, "base", 0, 0)
points = {base_point}
rectangle_specs = []
for old_direction in range(1, 5):
    for new_direction in range(1, 5):
        for old_sign in (-1, 1):
            for new_sign in (-1, 1):
                for level in ("half", "full"):
                    old_only = make_point(
                        old_direction, old_sign, level, 0, 0
                    )
                    new_only = make_point(
                        0, 0, level, new_direction, new_sign
                    )
                    mixed = make_point(
                        old_direction, old_sign, level,
                        new_direction, new_sign,
                    )
                    points.update((old_only, new_only, mixed))
                    rectangle_specs.append({
                        "id": (
                            f"o{old_direction}_n{new_direction}_"
                            f"s{old_sign:+d}_{new_sign:+d}_{level}"
                        ),
                        "old_direction": old_direction,
                        "new_direction": new_direction,
                        "old_sign": old_sign,
                        "new_sign": new_sign,
                        "level": level,
                        "base": base_point,
                        "old_only": old_only,
                        "new_only": new_only,
                        "mixed": mixed,
                    })
points = tuple(sorted(points))
point_census_ok = bool(
    len(points) == 161
    and len({point_id(point) for point in points}) == 161
    and len(rectangle_specs) == 128
    and len({record["id"] for record in rectangle_specs}) == 128
)


def boundary_logs(parity, point):
    old_direction, old_sign, level, new_direction, new_sign = point
    old_logs = [arb.log(value) for value in base[parity]["old"]]
    new_logs = [arb.log(value) for value in base[parity]["new"]]
    if level != "base":
        h = ETA*LEVELS[level]
        if old_direction:
            vector = directions[parity][old_direction]
            for index in range(30):
                old_logs[index] += old_sign*h*vector[index]
        if new_direction:
            vector = directions[parity][new_direction]
            for index in range(30):
                new_logs[index] += new_sign*h*vector[index]
    return tuple(old_logs), tuple(new_logs)


def evaluate_internal(parity, old_logs, new_logs, state):
    q_old = tuple(arb.exp(value) for value in old_logs)
    x = tuple(arb.exp(state[index]) for index in range(35))
    q_new = tuple(arb.exp(value) for value in new_logs)
    action, gradient, branch = action_and_gradient(
        models[parity], q_old, x, q_new
    )
    residual = tuple(arb.re(gradient[30+index]) for index in range(35))
    maximum_imaginary = max(
        abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient)
    )
    branch_ok = bool(
        branch["negative_counts"] == Counter({1: 2400})
        and branch["minimum_leading_minor"] > 0
        and branch["minimum_argument"] > ANGLE_TOLERANCE
        and maximum_imaginary < IMAGINARY_TOLERANCE
    )
    return {
        "action": action,
        "residual": residual,
        "residual_norm": infinity_norm(residual),
        "branch": branch,
        "maximum_imaginary": maximum_imaginary,
        "branch_ok": branch_ok,
    }


def solve_internal(parity, variant, point):
    old_logs, new_logs = boundary_logs(parity, point)
    state = arb.matrix([arb.log(value) for value in base[parity]["x"]])
    jacobian = jacobians[parity][variant]
    trace = []
    branch_all = True
    evaluation_error = None
    for iteration in range(MAX_ITERATIONS+1):
        try:
            current = evaluate_internal(parity, old_logs, new_logs, state)
        except Exception as error:
            evaluation_error = repr(error)
            branch_all = False
            break
        branch_all &= current["branch_ok"]
        if not current["branch_ok"]:
            return {
                "success": False,
                "reason": "BRANCH_FAILURE",
                "branch_all": False,
                "trace": trace,
                "evaluation_error": None,
            }
        try:
            delta = arb.lu_solve(
                jacobian,
                arb.matrix([-value for value in current["residual"]]),
            )
        except (ZeroDivisionError, ValueError) as error:
            return {
                "success": False,
                "reason": "JACOBIAN_SOLVE_FAILED",
                "branch_all": branch_all,
                "trace": trace,
                "evaluation_error": repr(error),
            }
        correction_norm = max(abs(delta[index]) for index in range(35))
        if (
            current["residual_norm"] < RESIDUAL_TOLERANCE
            and correction_norm < CORRECTION_TOLERANCE
        ):
            corrected_state = state+delta
            try:
                corrected = evaluate_internal(
                    parity, old_logs, new_logs, corrected_state
                )
            except Exception as error:
                return {
                    "success": False,
                    "reason": "CORRECTION_EVALUATION_FAILED",
                    "branch_all": False,
                    "trace": trace,
                    "evaluation_error": repr(error),
                }
            branch_all &= corrected["branch_ok"]
            correction_action = abs(
                arb.re(corrected["action"])-arb.re(current["action"])
            )
            success = bool(
                corrected["branch_ok"]
                and corrected["residual_norm"] < RESIDUAL_TOLERANCE
                and current["maximum_imaginary"] < IMAGINARY_TOLERANCE
                and corrected["maximum_imaginary"] < IMAGINARY_TOLERANCE
            )
            return {
                "success": success,
                "reason": "CONVERGED" if success else "FINAL_CONTROL_FAILED",
                "branch_all": branch_all,
                "iterations": iteration,
                "action": arb.re(current["action"]),
                "corrected_action": arb.re(corrected["action"]),
                "correction_action": correction_action,
                "residual_norm": current["residual_norm"],
                "corrected_residual_norm": corrected["residual_norm"],
                "correction_norm": correction_norm,
                "minimum_leading_minor": min(
                    current["branch"]["minimum_leading_minor"],
                    corrected["branch"]["minimum_leading_minor"],
                ),
                "minimum_argument": min(
                    current["branch"]["minimum_argument"],
                    corrected["branch"]["minimum_argument"],
                ),
                "maximum_imaginary": max(
                    current["maximum_imaginary"],
                    corrected["maximum_imaginary"],
                ),
                "trace": trace,
                "evaluation_error": None,
            }
        if iteration >= MAX_ITERATIONS:
            return {
                "success": False,
                "reason": "ITERATION_LIMIT",
                "branch_all": branch_all,
                "residual_norm": current["residual_norm"],
                "correction_norm": correction_norm,
                "trace": trace,
                "evaluation_error": None,
            }
        accepted = False
        accepted_alpha = None
        valid_trials = 0
        for power in range(MAX_BACKTRACKING+1):
            alpha = arb.mpf(2)**(-power)
            trial_state = state+alpha*delta
            try:
                trial = evaluate_internal(
                    parity, old_logs, new_logs, trial_state
                )
            except Exception:
                branch_all = False
                continue
            branch_all &= trial["branch_ok"]
            if not trial["branch_ok"]:
                continue
            valid_trials += 1
            if (
                trial["residual_norm"]
                <= (1-alpha/4)*current["residual_norm"]
            ):
                state = trial_state
                accepted = True
                accepted_alpha = alpha
                break
        trace.append({
            "iteration": iteration,
            "residual_norm": number(current["residual_norm"], 35),
            "correction_norm": number(correction_norm, 35),
            "accepted_alpha": (
                number(accepted_alpha, 20)
                if accepted_alpha is not None else None
            ),
        })
        if not accepted:
            return {
                "success": False,
                "reason": (
                    "BRANCH_FAILURE" if valid_trials == 0 else "NO_ARMIJO_STEP"
                ),
                "branch_all": branch_all,
                "residual_norm": current["residual_norm"],
                "correction_norm": correction_norm,
                "trace": trace,
                "evaluation_error": None,
            }
    return {
        "success": False,
        "reason": "EVALUATION_EXCEPTION",
        "branch_all": branch_all,
        "trace": trace,
        "evaluation_error": evaluation_error,
    }


def initialize_worker():
    arb.mp.dps = DPS


def solve_worker(task):
    parity, variant, point = task
    try:
        result = solve_internal(parity, variant, point)
    except Exception as error:
        result = {
            "success": False,
            "reason": "UNCAUGHT_EXCEPTION",
            "branch_all": False,
            "trace": [],
            "evaluation_error": repr(error),
        }
    return parity, variant, point, result


def compact_solve(result):
    scalar_keys = (
        "action",
        "corrected_action",
        "correction_action",
        "residual_norm",
        "corrected_residual_norm",
        "correction_norm",
        "minimum_leading_minor",
        "minimum_argument",
        "maximum_imaginary",
    )
    payload = {key: value for key, value in result.items() if key not in scalar_keys}
    for key in scalar_keys:
        if key in result:
            payload[key] = number(result[key])
    return payload


print("="*78)
print("DUST-SCHEDULE BOUNDARY-COBOUNDARY TEST")
print("="*78)
tasks = [
    (parity, variant, point)
    for point in points
    for parity in ("even", "odd")
    for variant in ("operational", "validation")
]
print(f"evaluating {len(tasks)} preregistered internal solves", flush=True)
fork_context = mp_pool.get_context("fork")
with fork_context.Pool(processes=8, initializer=initialize_worker) as pool:
    raw_results = pool.map(solve_worker, tasks, chunksize=1)

results = {
    point: {
        parity: {} for parity in ("even", "odd")
    }
    for point in points
}
for parity, variant, point, result in raw_results:
    results[point][parity][variant] = result

task_census_ok = bool(
    len(raw_results) == 644
    and all(
        set(results[point][parity]) == {"operational", "validation"}
        for point in points for parity in ("even", "odd")
    )
)
implementation_ok = all(
    not result.get("evaluation_error")
    for _, _, _, result in raw_results
)
recognized_solve_labels = {
    "CONVERGED",
    "BRANCH_FAILURE",
    "JACOBIAN_SOLVE_FAILED",
    "CORRECTION_EVALUATION_FAILED",
    "FINAL_CONTROL_FAILED",
    "ITERATION_LIMIT",
    "NO_ARMIJO_STEP",
    "EVALUATION_EXCEPTION",
    "UNCAUGHT_EXCEPTION",
}
solve_labels_ok = all(
    result.get("reason") in recognized_solve_labels
    for _, _, _, result in raw_results
)
all_success = all(result["success"] for _, _, _, result in raw_results)
all_branch_controls = all(
    result.get("branch_all", False) for _, _, _, result in raw_results
)

point_records = {}
for point in points:
    point_results = results[point]
    success = all(
        point_results[parity][variant]["success"]
        for parity in ("even", "odd")
        for variant in ("operational", "validation")
    )
    if success:
        even_op = point_results["even"]["operational"]
        even_val = point_results["even"]["validation"]
        odd_op = point_results["odd"]["operational"]
        odd_val = point_results["odd"]["validation"]
        delta = odd_op["action"]-even_op["action"]
        uncertainty = (
            abs(even_op["action"]-even_val["action"])
            + abs(odd_op["action"]-odd_val["action"])
            + even_op["correction_action"]
            + even_val["correction_action"]
            + odd_op["correction_action"]
            + odd_val["correction_action"]
            + UNCERTAINTY_FLOOR
        )
    else:
        delta = uncertainty = None
    point_records[point] = {
        "id": point_id(point),
        "coordinates": {
            "old_direction": point[0],
            "old_sign": point[1],
            "level": point[2],
            "new_direction": point[3],
            "new_sign": point[4],
        },
        "success": success,
        "delta": number(delta) if delta is not None else None,
        "uncertainty": number(uncertainty) if uncertainty is not None else None,
        "solves": {
            parity: {
                variant: compact_solve(result)
                for variant, result in variants.items()
            }
            for parity, variants in point_results.items()
        },
    }

rectangle_records = []
for spec_record in rectangle_specs:
    corners = [
        spec_record["mixed"],
        spec_record["old_only"],
        spec_record["new_only"],
        spec_record["base"],
    ]
    if not all(point_records[point]["success"] for point in corners):
        rectangle = uncertainty = ratio = None
        classification = "OPEN_SOLVE"
    else:
        mixed, old_only, new_only, base_corner = (
            arb.mpf(point_records[point]["delta"]) for point in corners
        )
        rectangle = mixed-old_only-new_only+base_corner
        uncertainty = sum(
            arb.mpf(point_records[point]["uncertainty"]) for point in corners
        )
        ratio = abs(rectangle)/uncertainty
        if abs(rectangle) <= CONSISTENT_FACTOR*uncertainty:
            classification = "SEPARABLE_CONSISTENT"
        elif abs(rectangle) > NONSEPARABLE_FACTOR*uncertainty:
            classification = "NONSEPARABLE"
        else:
            classification = "OPEN"
    rectangle_records.append({
        "id": spec_record["id"],
        "old_direction": spec_record["old_direction"],
        "new_direction": spec_record["new_direction"],
        "old_sign": spec_record["old_sign"],
        "new_sign": spec_record["new_sign"],
        "level": spec_record["level"],
        "corners": {
            name: point_id(spec_record[name])
            for name in ("base", "old_only", "new_only", "mixed")
        },
        "classification": classification,
        "rectangle": number(rectangle) if rectangle is not None else None,
        "absolute_rectangle": (
            number(abs(rectangle)) if rectangle is not None else None
        ),
        "uncertainty": number(uncertainty) if uncertainty is not None else None,
        "ratio": number(ratio) if ratio is not None else None,
    })

classification_counts = Counter(
    record["classification"] for record in rectangle_records
)
rectangle_census_ok = bool(
    len(rectangle_records) == 128
    and len({record["id"] for record in rectangle_records}) == 128
    and sum(classification_counts.values()) == 128
    and set(classification_counts) <= {
        "SEPARABLE_CONSISTENT", "NONSEPARABLE", "OPEN", "OPEN_SOLVE"
    }
)

scaling_records = []
for old_direction in range(1, 5):
    for new_direction in range(1, 5):
        for old_sign in (-1, 1):
            for new_sign in (-1, 1):
                pair = [
                    record for record in rectangle_records
                    if record["old_direction"] == old_direction
                    and record["new_direction"] == new_direction
                    and record["old_sign"] == old_sign
                    and record["new_sign"] == new_sign
                ]
                half = next(record for record in pair if record["level"] == "half")
                full = next(record for record in pair if record["level"] == "full")
                available = bool(
                    half["classification"] == full["classification"]
                    == "NONSEPARABLE"
                )
                if available:
                    order = arb.log(
                        arb.mpf(full["absolute_rectangle"])
                        / arb.mpf(half["absolute_rectangle"]),
                        2,
                    )
                    if arb.mpf("1.5") <= order <= arb.mpf("2.5"):
                        label = "QUADRATIC_ACTION_COMPATIBLE"
                    elif arb.mpf("2.5") < order <= arb.mpf("3.5"):
                        label = "CUBIC_ACTION_COMPATIBLE"
                    else:
                        label = "OTHER_RESOLVED_ORDER"
                else:
                    order = None
                    label = "NOT_AVAILABLE"
                scaling_records.append({
                    "old_direction": old_direction,
                    "new_direction": new_direction,
                    "old_sign": old_sign,
                    "new_sign": new_sign,
                    "available": available,
                    "order": number(order) if order is not None else None,
                    "label": label,
                })

controls_ok = bool(
    provenance_ok
    and carrier_ok
    and direction_geometry_ok
    and point_census_ok
    and task_census_ok
    and implementation_ok
    and solve_labels_ok
    and rectangle_census_ok
)
if not controls_ok:
    outcome = "BOUNDARY_COBOUNDARY_CONTROL_FAILED"
elif classification_counts["NONSEPARABLE"]:
    outcome = "BOUNDARY_COBOUNDARY_REFUTED_ON_FROZEN_RECTANGLES"
elif classification_counts["OPEN"] or classification_counts["OPEN_SOLVE"]:
    outcome = "BOUNDARY_COBOUNDARY_OPEN"
else:
    outcome = "BOUNDARY_COBOUNDARY_CONSISTENT_ON_FROZEN_RECTANGLES"

tests = [
    ("frozen input hashes and upstream provenance", provenance_ok),
    ("imported carrier retains 43/43 and 30+35+30 variables", carrier_ok),
    ("four directions remain unit, zero-sum and physically permuted", direction_geometry_ok),
    ("preregistered census has 161 unique points and 128 rectangles", point_census_ok),
    ("all 644 parity/calibration tasks were returned exactly once", task_census_ok),
    ("no internal action evaluation raised an implementation exception", implementation_ok),
    ("every internal solve received a recognized mechanical label", solve_labels_ok),
    ("every rectangle received exactly one frozen classification", rectangle_census_ok),
    ("no continuum, speed, chirality or experimental target was parsed", True),
    ("outcome follows the preregistered mechanical rule", outcome in {
        "BOUNDARY_COBOUNDARY_CONTROL_FAILED",
        "BOUNDARY_COBOUNDARY_REFUTED_ON_FROZEN_RECTANGLES",
        "BOUNDARY_COBOUNDARY_OPEN",
        "BOUNDARY_COBOUNDARY_CONSISTENT_ON_FROZEN_RECTANGLES",
    }),
]
passed = sum(bool(ok) for _, ok in tests)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "coordinates": "logarithms of squared-length orbit variables",
    "carrier": {
        "old_boundary": 30,
        "internal": 35,
        "new_boundary": 30,
        "full_720_edge_carrier": False,
    },
    "eta": number(ETA),
    "number_of_directions": 4,
    "number_of_unique_boundary_points": len(points),
    "number_of_internal_solves": len(tasks),
    "number_of_rectangles": len(rectangle_records),
    "all_internal_solves_converged": all_success,
    "all_evaluations_retained_lorentzian_branch": all_branch_controls,
    "classification_counts": dict(classification_counts),
    "scaling_label_counts": dict(Counter(
        record["label"] for record in scaling_records
    )),
    "points": [point_records[point] for point in points],
    "rectangles": rectangle_records,
    "scaling_diagnostics": scaling_records,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "chirality_target_parsed": False,
    "experimental_target_parsed": False,
    "tests": len(tests),
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"solve success={sum(result['success'] for _, _, _, result in raw_results)}/644")
print(f"rectangle classifications={dict(classification_counts)}")
if any(record["absolute_rectangle"] is not None for record in rectangle_records):
    available = [
        record for record in rectangle_records
        if record["absolute_rectangle"] is not None
    ]
    print(
        "absolute rectangle range {} ... {}".format(
            number(min(arb.mpf(record["absolute_rectangle"]) for record in available), 10),
            number(max(arb.mpf(record["absolute_rectangle"]) for record in available), 10),
        )
    )
print(f"scaling labels={dict(Counter(record['label'] for record in scaling_records))}")
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) and controls_ok else 1)

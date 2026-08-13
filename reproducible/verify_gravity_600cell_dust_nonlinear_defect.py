#!/usr/bin/env python3
"""Raw-residual defect correction for nonlinear relative-phase cases.

Protocol commit: 2d695c6.  Every solve starts from the stored final state of
the frozen first nonlinear run.  Fixed-Hessian steps are accepted only when
the physical transverse equation norm decreases; candidates retain the same
100-decimal complete-action validation.
"""

import contextlib
import importlib.util
import io
import json
import math
import multiprocessing as mp
from pathlib import Path
import sys

import mpmath as arb
import numpy as np


HERE = Path(__file__).resolve().parent
PRECISION_INPUT = HERE / "gravity_600cell_dust_gauge_quotient_precision.json"
FIRST_INPUT = HERE / "gravity_600cell_dust_nonlinear_relative.json"
OUTPUT = HERE / "gravity_600cell_dust_nonlinear_defect.json"
PROTOCOL_COMMIT = "2d695c6"
DIAGNOSTIC_COMMIT = "daf5a64"
FIRST_RESULT_COMMIT = "b4a7828"
PRECISION_RESULT_COMMIT = "29a779f"
T_GRID = (-0.10, -0.05, 0.0, 0.05, 0.10)
DAMPING = tuple(2.0**(-power) for power in range(11))
ACTION_STEPS = (2.0e-5, 1.0e-5, 5.0e-6)
SCALAR_FLOOR = 1.0e-10
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


precision_input = json.loads(PRECISION_INPUT.read_text())
first_input = json.loads(FIRST_INPUT.read_text())
check(
    "the precision quotient and first nonlinear run retain frozen provenance",
    precision_input["protocol_commit"] == "da34272"
    and first_input["protocol_commit"] == "80f8de7"
    and precision_input["passed"] == precision_input["tests"] == 14
    and first_input["passed"] == first_input["tests"] == 10
    and first_input["eta"] == 1.0e-4
    and tuple(first_input["collective_scan"]) == T_GRID
    and first_input["hit_fractions"]["signed_cases"] == {
        "hits": 0, "total": 16
    },
)


print("Loading the independently certified complete-action implementation...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_nonlinear_defect",
    HERE / "verify_gravity_600cell_published_dust_control.py",
)
dust = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = dust
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(dust)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check(
    "the imported complete-action control retains all 14 certificates",
    dust.tests == dust.passed == 14,
)


quotient_basis = np.array(precision_input["quotient_basis"], dtype=float)
direction_input = first_input["direction_enumeration"]
check(
    "the stored enumeration retains four rank-four directions and sixteen cases",
    direction_input["count"] == 4
    and direction_input["absolute_rank_above_1e-8"] == 4
    and len(direction_input["directions"]) == 4
    and all(len(record["cases"]) == 8 for record in first_input["parities"].values())
    and sum(
        len(record["cases"]) for record in first_input["parities"].values()
    ) == 16,
)


def path_internal(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    return np.concatenate((np.full(30, diagonal), np.full(5, rho)))


def path_tangent(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    tangent = np.concatenate((
        np.full(30, -rho/diagonal), np.ones(5)
    ))
    return tangent/np.linalg.norm(tangent)


def complete_variables(t, z, boundary_log):
    internal = path_internal(t)*np.exp(quotient_basis@z)
    final = dust.L0_SQUARE*np.exp(boundary_log)
    return np.concatenate((internal, final))


def equation_evaluation(model, t, z, boundary_log):
    variables = complete_variables(t, z, boundary_log)
    try:
        _, gradient, _, data = dust.total_reduced_evaluation(
            model, variables, dust.old_values
        )
        equation_complex = variables[:35]*gradient[:35]/24.0
        equation = equation_complex.real
        finite = bool(
            np.all(np.isfinite(equation_complex.real))
            and np.all(np.isfinite(equation_complex.imag))
        )
        branch = bool(
            dict(data["negative_counts"]) == {1: 100}
            and float(data["minimum_gram"]) > 1e-8
            and float(data["minimum_argument"]) > 1e-6
        )
        return {
            "variables": variables,
            "equation": equation,
            "maximum_imaginary": float(np.max(np.abs(equation_complex.imag))),
            "minimum_gram": float(data["minimum_gram"]),
            "minimum_argument": float(data["minimum_argument"]),
            "lorentzian": dict(data["negative_counts"]) == {1: 100},
            "branch": branch,
            "finite": finite,
        }
    except Exception as error:
        return {
            "variables": variables,
            "equation": np.full(35, np.nan),
            "maximum_imaginary": math.inf,
            "minimum_gram": -math.inf,
            "minimum_argument": -math.inf,
            "lorentzian": False,
            "branch": False,
            "finite": False,
            "error": repr(error),
        }


def solve_defect(model, quotient_hessian, t, initial_z, boundary_log):
    z = np.array(initial_z, dtype=float, copy=True)
    history = []
    trials = []
    minimum_gram = math.inf
    minimum_argument = math.inf
    maximum_imaginary = 0.0
    accepted_iterations = 0
    stop_reason = "ITERATION_LIMIT"

    for iteration in range(51):
        state = equation_evaluation(model, t, z, boundary_log)
        minimum_gram = min(minimum_gram, state["minimum_gram"])
        minimum_argument = min(minimum_argument, state["minimum_argument"])
        maximum_imaginary = max(
            maximum_imaginary, state["maximum_imaginary"]
        )
        if not state["finite"] or not state["branch"]:
            stop_reason = "CURRENT_STATE_INVALID"
            break
        transverse = quotient_basis.T@state["equation"]
        correction = np.linalg.solve(quotient_hessian, transverse)
        transverse_norm = float(np.linalg.norm(transverse))
        correction_norm = float(np.linalg.norm(correction))
        history.append({
            "iteration": iteration,
            "transverse_norm": transverse_norm,
            "correction_norm": correction_norm,
            "accepted_damping": None,
        })
        if transverse_norm < 1e-12 and correction_norm < 1e-5:
            stop_reason = "TRANSVERSE_GATES_PASSED"
            break
        if iteration == 50:
            break

        accepted = False
        for damping in DAMPING:
            trial_z = z-damping*correction
            trial = equation_evaluation(model, t, trial_z, boundary_log)
            minimum_gram = min(minimum_gram, trial["minimum_gram"])
            minimum_argument = min(
                minimum_argument, trial["minimum_argument"]
            )
            maximum_imaginary = max(
                maximum_imaginary, trial["maximum_imaginary"]
            )
            trial_record = {
                "iteration": iteration,
                "damping": damping,
                "finite": trial["finite"],
                "branch": trial["branch"],
                "transverse_norm": math.inf,
            }
            if trial["finite"] and trial["branch"]:
                trial_transverse = quotient_basis.T@trial["equation"]
                trial_record["transverse_norm"] = float(
                    np.linalg.norm(trial_transverse)
                )
                if trial_record["transverse_norm"] < transverse_norm:
                    z = trial_z
                    history[-1]["accepted_damping"] = damping
                    accepted_iterations += 1
                    accepted = True
            trials.append(trial_record)
            if accepted:
                break
        if not accepted:
            stop_reason = "NO_RAW_RESIDUAL_DECREASING_STEP"
            break

    final = equation_evaluation(model, t, z, boundary_log)
    minimum_gram = min(minimum_gram, final["minimum_gram"])
    minimum_argument = min(minimum_argument, final["minimum_argument"])
    maximum_imaginary = max(maximum_imaginary, final["maximum_imaginary"])
    if final["finite"]:
        transverse = quotient_basis.T@final["equation"]
        correction = np.linalg.solve(quotient_hessian, transverse)
        transverse_norm = float(np.linalg.norm(transverse))
        correction_norm = float(np.linalg.norm(correction))
        scalar = float(path_tangent(t)@final["equation"])
        full_norm = float(np.linalg.norm(final["equation"]))
    else:
        transverse = np.full(34, np.nan)
        correction = np.full(34, np.nan)
        transverse_norm = correction_norm = full_norm = math.inf
        scalar = math.nan
    resolved = bool(
        final["finite"] and final["branch"]
        and transverse_norm < 1e-12
        and correction_norm < 1e-5
    )
    return {
        "t": float(t),
        "z": z,
        "outcome": (
            "DEFECT_TRANSVERSE_SOLVED" if resolved
            else "DEFECT_TRANSVERSE_UNRESOLVED"
        ),
        "stop_reason": stop_reason,
        "accepted_iterations": accepted_iterations,
        "history": history,
        "trials": trials,
        "variables": final["variables"],
        "equation": final["equation"],
        "transverse": transverse,
        "correction": correction,
        "transverse_norm": transverse_norm,
        "correction_norm": correction_norm,
        "scalar": scalar,
        "full_equation_norm": full_norm,
        "final_branch": final["branch"],
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
    }


def bisect_scalar(model, quotient_hessian, left, right, boundary_log):
    history = []
    current_left, current_right = left, right
    if current_left["scalar"]*current_right["scalar"] >= 0:
        return None, history, False
    for _iteration in range(30):
        midpoint_t = (current_left["t"]+current_right["t"])/2.0
        initial_z = (current_left["z"]+current_right["z"])/2.0
        midpoint = solve_defect(
            model, quotient_hessian, midpoint_t, initial_z, boundary_log
        )
        history.append(midpoint)
        if midpoint["outcome"] != "DEFECT_TRANSVERSE_SOLVED":
            return None, history, False
        if (
            abs(midpoint["scalar"]) <= SCALAR_FLOOR
            or current_right["t"]-current_left["t"] <= 1e-8
        ):
            return midpoint, history, True
        if current_left["scalar"]*midpoint["scalar"] < 0:
            current_right = midpoint
        else:
            current_left = midpoint
    return history[-1], history, True


def compact_solve(record):
    return {
        "t": record["t"],
        "z": record["z"].tolist(),
        "outcome": record["outcome"],
        "stop_reason": record["stop_reason"],
        "accepted_iterations": record["accepted_iterations"],
        "history": record["history"],
        "trials": record["trials"],
        "variables": record["variables"].tolist(),
        "equation": record["equation"].tolist(),
        "transverse": record["transverse"].tolist(),
        "correction": record["correction"].tolist(),
        "transverse_norm": record["transverse_norm"],
        "correction_norm": record["correction_norm"],
        "scalar": record["scalar"],
        "full_equation_norm": record["full_equation_norm"],
        "final_branch": record["final_branch"],
        "minimum_gram": record["minimum_gram"],
        "minimum_argument": record["minimum_argument"],
        "maximum_imaginary": record["maximum_imaginary"],
    }


def configure_arb_precision():
    arb.mp.dps = 100
    dust.ARB_I = arb.mpc(0, 1)
    dust.ARB_TAU = arb.mpf("0.0102")
    dust.ARB_M_STAR = arb.mpf(10)
    dust.ARB_ZETA = (
        arb.pi**2*arb.sqrt(2)/50
    )**(arb.mpf(1)/3)
    dust.ARB_R0 = 4*dust.ARB_M_STAR/(3*arb.pi)
    dust.ARB_L0 = dust.ARB_ZETA*dust.ARB_R0
    dust.ARB_L0_SQUARE = dust.ARB_L0**2
    dust.ARB_EPSILON_3 = 2*arb.pi-5*arb.acos(arb.mpf(1)/3)
    dust.ARB_MASS = (90/arb.pi)*dust.ARB_EPSILON_3*dust.ARB_L0
    dust.ARB_TAU_SQUARE = dust.ARB_TAU**2
    dust.ARB_SLANT_SQUARE = dust.ARB_L0_SQUARE-dust.ARB_TAU_SQUARE
    dust.ARB_OLD_VALUES[:] = [dust.ARB_L0_SQUARE]*30


_ARB_MODEL = None
_ARB_CANDIDATES = None


def initialize_candidate_worker(model, candidate_strings):
    global _ARB_MODEL, _ARB_CANDIDATES
    _ARB_MODEL = model
    configure_arb_precision()
    _ARB_CANDIDATES = [
        [arb.mpf(value) for value in candidate]
        for candidate in candidate_strings
    ]


def candidate_action_worker(point):
    candidate_index, coordinate, logarithmic_step = point
    variables = list(_ARB_CANDIDATES[candidate_index])
    variables[coordinate] *= arb.exp(arb.mpf(str(logarithmic_step)))
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


def audit_candidates(model, quotient_hessian, candidates):
    if not candidates:
        return []
    points = []
    branch_records = []
    for candidate_index, candidate in enumerate(candidates):
        variables = candidate["variables"]
        branch_pass = True
        minimum_gram = math.inf
        minimum_argument = math.inf
        for step in ACTION_STEPS:
            for coordinate in range(35):
                for sign in (1.0, -1.0):
                    logarithmic_step = sign*step
                    points.append((candidate_index, coordinate, logarithmic_step))
                    trial = variables.copy()
                    trial[coordinate] *= math.exp(logarithmic_step)
                    _, _, _, data = dust.total_reduced_evaluation(
                        model, trial, dust.old_values
                    )
                    minimum_gram = min(
                        minimum_gram, float(data["minimum_gram"])
                    )
                    minimum_argument = min(
                        minimum_argument, float(data["minimum_argument"])
                    )
                    branch_pass &= bool(
                        dict(data["negative_counts"]) == {1: 100}
                        and float(data["minimum_gram"]) > 1e-8
                        and float(data["minimum_argument"]) > 1e-6
                    )
        branch_records.append({
            "pass": branch_pass,
            "minimum_gram": minimum_gram,
            "minimum_argument": minimum_argument,
        })

    candidate_strings = [
        [f"{value:.17g}" for value in candidate["variables"]]
        for candidate in candidates
    ]
    fork_context = mp.get_context("fork")
    print(
        f"Validating {len(candidates)} defect candidates with "
        f"{len(points)} action points...",
        flush=True,
    )
    with fork_context.Pool(
        processes=8,
        initializer=initialize_candidate_worker,
        initargs=(model, candidate_strings),
    ) as pool:
        action_values = pool.map(candidate_action_worker, points, chunksize=1)

    audits = []
    cursor = 0
    for candidate_index, candidate in enumerate(candidates):
        derivatives = {}
        maximum_imaginary = arb.mpf(0)
        for step in ACTION_STEPS:
            step_arb = arb.mpf(str(step))
            row = []
            for _coordinate in range(35):
                plus, minus = action_values[cursor:cursor+2]
                cursor += 2
                row.append((plus-minus)/(48*step_arb))
                maximum_imaginary = max(
                    maximum_imaginary,
                    abs(arb.im(plus)), abs(arb.im(minus)),
                )
            derivatives[step] = row
        row1, row2, row3 = (derivatives[step] for step in ACTION_STEPS)
        richardson12 = [(4*b-a)/3 for a, b in zip(row1, row2)]
        richardson23 = [(4*c-b)/3 for b, c in zip(row2, row3)]
        sixth_order = [
            (16*r23-r12)/15
            for r12, r23 in zip(richardson12, richardson23)
        ]
        equation = np.array([float(arb.re(value)) for value in sixth_order])
        richardson23_float = np.array([
            float(arb.re(value)) for value in richardson23
        ])
        error_row = equation-richardson23_float
        equation_norm = float(np.linalg.norm(equation))
        equation_error = float(np.linalg.norm(error_row))
        transverse_correction = float(np.linalg.norm(
            np.linalg.solve(quotient_hessian, quotient_basis.T@equation)
        ))
        transverse_error = float(np.linalg.norm(
            np.linalg.solve(quotient_hessian, quotient_basis.T@error_row)
        ))
        tangent = path_tangent(candidate["t"])
        scalar = float(tangent@equation)
        scalar_error = float(abs(tangent@error_row))
        binary_agreement = float(np.linalg.norm(
            equation-candidate["equation"]
        ))
        maximum_imaginary = float(max(
            maximum_imaginary,
            max(abs(arb.im(value)) for value in sixth_order),
        ))
        gates = {
            "equation_zero_consistent": equation_norm <= 10*max(equation_error, 1e-30),
            "transverse_zero_consistent": transverse_correction <= 10*max(transverse_error, 1e-30),
            "scalar_zero_consistent": abs(scalar) <= 10*max(scalar_error, 1e-30),
            "transverse_error_below_1e-5": transverse_error < 1e-5,
            "transverse_correction_below_1e-5": transverse_correction < 1e-5,
            "all_210_branches_pass": branch_records[candidate_index]["pass"],
            "imaginary_below_1e-80": maximum_imaginary < 1e-80,
            "binary_action_agreement_below_1e-8": binary_agreement < 1e-8,
        }
        audits.append({
            "label": (
                "NONLINEAR_STATIONARY_CANDIDATE" if all(gates.values())
                else "ACTION_ONLY_CANDIDATE_FAILED"
            ),
            "gates": gates,
            "derivatives": derivatives,
            "richardson12": richardson12,
            "richardson23": richardson23,
            "sixth_order": sixth_order,
            "equation": equation,
            "error_row": error_row,
            "equation_norm": equation_norm,
            "equation_error": equation_error,
            "transverse_correction": transverse_correction,
            "transverse_error": transverse_error,
            "scalar": scalar,
            "scalar_error": scalar_error,
            "binary_agreement": binary_agreement,
            "maximum_imaginary": maximum_imaginary,
            "branch": branch_records[candidate_index],
        })
    if cursor != len(action_values):
        raise RuntimeError("defect candidate action cursor mismatch")
    return audits


results = {}
for parity, model in dust.bl.models.items():
    print(f"Applying raw-residual defect correction: {parity} parity...", flush=True)
    quotient_hessian = np.array(
        precision_input["parities"][parity]["quotient"]["matrix"],
        dtype=float,
    )
    quotient_hessian = (quotient_hessian+quotient_hessian.T)/2.0
    first_cases = first_input["parities"][parity]["cases"]
    cases = []
    all_candidates = []

    for first_case in first_cases:
        boundary_log = np.array(
            first_case["boundary_log_vector"], dtype=float
        )
        grid = {}
        for t in T_GRID:
            key = f"{t:+.2f}"
            stored = first_case["grid"][key]
            grid[t] = solve_defect(
                model,
                quotient_hessian,
                t,
                np.array(stored["z"], dtype=float),
                boundary_log,
            )
        ordered = [grid[t] for t in T_GRID]
        candidates = [
            record for record in ordered
            if record["outcome"] == "DEFECT_TRANSVERSE_SOLVED"
            and abs(record["scalar"]) <= SCALAR_FLOOR
        ]
        bisections = []
        intervals_resolved = True
        for left, right in zip(ordered[:-1], ordered[1:]):
            if (
                left["outcome"] != "DEFECT_TRANSVERSE_SOLVED"
                or right["outcome"] != "DEFECT_TRANSVERSE_SOLVED"
            ):
                intervals_resolved = False
                continue
            if (
                abs(left["scalar"]) > SCALAR_FLOOR
                and abs(right["scalar"]) > SCALAR_FLOOR
                and left["scalar"]*right["scalar"] < 0
            ):
                candidate, history, resolved = bisect_scalar(
                    model, quotient_hessian, left, right, boundary_log
                )
                bisections.append({
                    "interval": [left["t"], right["t"]],
                    "resolved": resolved,
                    "history": history,
                    "candidate": candidate,
                })
                intervals_resolved &= resolved
                if candidate is not None:
                    candidates.append(candidate)
        candidates.sort(key=lambda record: record["t"])
        deduplicated = []
        for candidate in candidates:
            if not deduplicated or abs(candidate["t"]-deduplicated[-1]["t"]) >= 1e-6:
                deduplicated.append(candidate)
            elif abs(candidate["scalar"]) < abs(deduplicated[-1]["scalar"]):
                deduplicated[-1] = candidate
        candidates = deduplicated
        candidate_indices = []
        for candidate in candidates:
            candidate_indices.append(len(all_candidates))
            all_candidates.append(candidate)
        cases.append({
            "direction_index": int(first_case["direction_index"]),
            "sign": int(first_case["sign"]),
            "boundary_log": boundary_log,
            "grid": grid,
            "bisections": bisections,
            "candidate_indices": candidate_indices,
            "all_grid_resolved": all(
                record["outcome"] == "DEFECT_TRANSVERSE_SOLVED"
                for record in ordered
            ),
            "all_intervals_resolved": intervals_resolved,
        })

    audits = audit_candidates(model, quotient_hessian, all_candidates)
    for case in cases:
        case["candidate_records"] = [
            all_candidates[index] for index in case["candidate_indices"]
        ]
        case["audits"] = [
            audits[index] for index in case["candidate_indices"]
        ]
        case["passing_candidates"] = sum(
            audit["label"] == "NONLINEAR_STATIONARY_CANDIDATE"
            for audit in case["audits"]
        )
        case["localization_resolved"] = bool(
            case["all_grid_resolved"] and case["all_intervals_resolved"]
        )
        if case["passing_candidates"] > 0:
            case["outcome"] = "NONLINEAR_CONTINUATION_FOUND_IN_SCAN"
        elif case["localization_resolved"]:
            case["outcome"] = "NO_NONLINEAR_CONTINUATION_IN_FROZEN_SCAN"
        else:
            case["outcome"] = "NONLINEAR_CONTINUATION_NUMERICALLY_UNRESOLVED"

    results[parity] = {
        "quotient_hessian": quotient_hessian,
        "cases": cases,
        "all_candidates": all_candidates,
        "audits": audits,
    }
    check(
        f"{parity}: all eight stored cases and forty stored starts were reused",
        len(cases) == 8
        and sum(len(case["grid"]) for case in cases) == 40
        and {(case["direction_index"], case["sign"]) for case in cases}
        == {(index, sign) for index in range(1, 5) for sign in (-1, 1)},
    )
    check(
        f"{parity}: every localized defect candidate has an action audit",
        len(all_candidates) == len(audits),
        f"candidates={len(audits)}, passing={sum(a['label']=='NONLINEAR_STATIONARY_CANDIDATE' for a in audits)}",
    )
    check(
        f"{parity}: all eight corrected outcomes follow the frozen rules",
        all(case["outcome"] in {
            "NONLINEAR_CONTINUATION_FOUND_IN_SCAN",
            "NO_NONLINEAR_CONTINUATION_IN_FROZEN_SCAN",
            "NONLINEAR_CONTINUATION_NUMERICALLY_UNRESOLVED",
        } for case in cases),
    )


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 80),
        "imaginary": arb.nstr(arb.im(value), 80),
    }


def serialize_audit(record):
    return {
        "label": record["label"],
        "gates": record["gates"],
        "derivatives": {
            f"{step:.1e}": [
                serialize_complex(value)
                for value in record["derivatives"][step]
            ] for step in ACTION_STEPS
        },
        "richardson12": [
            serialize_complex(value) for value in record["richardson12"]
        ],
        "richardson23": [
            serialize_complex(value) for value in record["richardson23"]
        ],
        "sixth_order": [
            serialize_complex(value) for value in record["sixth_order"]
        ],
        "sixth_order_float": record["equation"].tolist(),
        "error_row_float": record["error_row"].tolist(),
        "equation_norm": record["equation_norm"],
        "equation_error_norm": record["equation_error"],
        "transverse_correction_norm": record["transverse_correction"],
        "transverse_error_norm": record["transverse_error"],
        "collective_scalar": record["scalar"],
        "collective_scalar_error": record["scalar_error"],
        "binary_vs_action_equation_norm": record["binary_agreement"],
        "maximum_action_or_derivative_imaginary": record["maximum_imaginary"],
        "minimum_absolute_gram_eigenvalue": record["branch"]["minimum_gram"],
        "minimum_angle_argument_modulus": record["branch"]["minimum_argument"],
        "all_210_branches_pass": record["branch"]["pass"],
    }


def serialize_case(case):
    return {
        "direction_index": case["direction_index"],
        "sign": case["sign"],
        "outcome": case["outcome"],
        "passing_candidate_count": case["passing_candidates"],
        "all_grid_resolved": case["all_grid_resolved"],
        "all_intervals_resolved": case["all_intervals_resolved"],
        "localization_resolved": case["localization_resolved"],
        "boundary_log_vector": case["boundary_log"].tolist(),
        "grid": {
            f"{t:+.2f}": compact_solve(record)
            for t, record in sorted(case["grid"].items())
        },
        "bisections": [
            {
                "interval": record["interval"],
                "resolved": record["resolved"],
                "history": [compact_solve(item) for item in record["history"]],
                "candidate": None if record["candidate"] is None else compact_solve(record["candidate"]),
            }
            for record in case["bisections"]
        ],
        "candidates": [
            {
                "solve": compact_solve(candidate),
                "action_only": serialize_audit(audit),
            }
            for candidate, audit in zip(
                case["candidate_records"], case["audits"]
            )
        ],
    }


signed_hits = sum(
    case["outcome"] == "NONLINEAR_CONTINUATION_FOUND_IN_SCAN"
    for record in results.values() for case in record["cases"]
)
pair_hits = 0
contrast_hits = 0
for direction_index in range(1, 5):
    contrast_hit = False
    for parity, record in results.items():
        pair = [
            case for case in record["cases"]
            if case["direction_index"] == direction_index
        ]
        hit = any(
            case["outcome"] == "NONLINEAR_CONTINUATION_FOUND_IN_SCAN"
            for case in pair
        )
        pair_hits += hit
        contrast_hit |= hit
    contrast_hits += contrast_hit

check(
    "corrected hit fractions cover 16 cases, 8 pairs and 4 contrasts",
    0 <= signed_hits <= 16
    and 0 <= pair_hits <= 8
    and 0 <= contrast_hits <= 4,
    f"signed={signed_hits}/16, pairs={pair_hits}/8, contrasts={contrast_hits}/4",
)

if signed_hits > 0:
    verdict = (
        "DERIVED COMPUTATIONAL LOCAL WITH PATTERN-INFORMED METHOD: at least "
        "one preregistered relative-phase case has a complete action-validated "
        "nonlinear continuation; robustness is given by the hit fractions."
    )
elif all(
    case["outcome"] == "NO_NONLINEAR_CONTINUATION_IN_FROZEN_SCAN"
    for record in results.values() for case in record["cases"]
):
    verdict = (
        "DERIVED COMPUTATIONAL NEGATIVE WITHIN SCOPE: all sixteen scans "
        "resolved and none continued at eta=1e-4 within |t|<=0.1."
    )
else:
    verdict = (
        "OPEN NUMERICALLY: no action-validated continuation was found and at "
        "least one raw-residual corrected scan remained unresolved."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "post_result_diagnostic_commit": DIAGNOSTIC_COMMIT,
    "first_nonlinear_result_commit": FIRST_RESULT_COMMIT,
    "precision_response_result_commit": PRECISION_RESULT_COMMIT,
    "eta": first_input["eta"],
    "collective_scan": T_GRID,
    "damping_sequence": DAMPING,
    "action_precision_digits": 100,
    "action_derivative_steps": ACTION_STEPS,
    "parities": {
        parity: {
            "cases": [serialize_case(case) for case in record["cases"]]
        }
        for parity, record in results.items()
    },
    "hit_fractions": {
        "signed_cases": {"hits": int(signed_hits), "total": 16},
        "direction_parity_pairs": {"hits": int(pair_hits), "total": 8},
        "phase_contrasts": {"hits": int(contrast_hits), "total": 4},
    },
    "verdict": verdict,
    "claim_boundary": {
        "method_provenance": "PATTERN-INFORMED BY ONE-STEP DIAGNOSTIC",
        "nonlinear_relative_continuation": "DERIVED COMPUTATIONAL IF VALIDATED",
        "amplitude_scaling": "NOT TESTED",
        "all_29_boundary_directions": "NOT TESTED",
        "second_slab": "NOT TESTED",
        "full_840_edge_carrier": "NOT TESTED",
        "clock_speed_limit_planck_scale": "NOT DERIVED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)

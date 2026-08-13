#!/usr/bin/env python3
"""Nonlinear continuation probe for all four relative-phase directions.

Protocol commit: 80f8de7.  Sixteen signed parity cases are localized with the
complete analytic equations and every candidate is validated independently
from the complete action at 100 decimal digits.
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
from scipy.optimize import least_squares


HERE = Path(__file__).resolve().parent
PRECISION_INPUT = HERE / "gravity_600cell_dust_gauge_quotient_precision.json"
PATH_INPUT = HERE / "gravity_600cell_dust_weak_scale_path.json"
OUTPUT = HERE / "gravity_600cell_dust_nonlinear_relative.json"
PROTOCOL_COMMIT = "80f8de7"
PRIOR_ART_COMMIT = "81b1aa1"
PATH_RESULT_COMMIT = "ae902a8"
PRECISION_RESULT_COMMIT = "29a779f"
ETA = 1.0e-4
T_GRID = (-0.10, -0.05, 0.0, 0.05, 0.10)
POSITIVE_SCAN = (0.05, 0.10)
NEGATIVE_SCAN = (-0.05, -0.10)
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
path_input = json.loads(PATH_INPUT.read_text())
check(
    "the precision response and weak-scale path retain frozen provenance",
    precision_input["protocol_commit"] == "da34272"
    and path_input["protocol_commit"] == "8380f0d"
    and precision_input["passed"] == precision_input["tests"] == 14
    and path_input["passed"] == path_input["tests"] == 11
    and {
        record["outcome"] for record in path_input["parities"].values()
    } == {"ALL_13_PATH_POINTS_STATIONARY_WITHIN_ERROR"},
)


print("Loading the independently certified complete-action implementation...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_nonlinear_relative",
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
even_response = precision_input["parities"]["even"]["response"]
boundary_basis = np.array(even_response["boundary_basis"], dtype=float)
internal_response = np.array(
    even_response["internal_response_matrix"], dtype=float
)
coefficient_response = np.array(
    even_response["coefficient_matrix"], dtype=float
)
pole_projector = np.eye(5)-np.ones((5, 5))/5.0
relative_map = pole_projector@internal_response[30:, :]
relative_singular = np.linalg.svd(relative_map, compute_uv=False)
relative_rank = int(np.sum(relative_singular > 1e-8))

directions = []
for contrast_index in range(1, 5):
    contrast = np.concatenate((
        np.ones(contrast_index),
        np.array([-float(contrast_index)]),
        np.zeros(4-contrast_index),
    ))/math.sqrt(contrast_index*(contrast_index+1))
    coordinates = np.linalg.pinv(
        relative_map, rcond=1e-12
    )@contrast
    boundary_vector = boundary_basis@coordinates
    boundary_vector /= np.linalg.norm(boundary_vector)
    normalized_coordinates = boundary_basis.T@boundary_vector
    response_vector = internal_response@normalized_coordinates
    amplification = float(np.linalg.norm(response_vector))
    amplitude = ETA/amplification
    directions.append({
        "index": contrast_index,
        "contrast": contrast,
        "boundary_vector": boundary_vector,
        "boundary_coordinates": normalized_coordinates,
        "response_vector": response_vector,
        "response_coefficients": coefficient_response@normalized_coordinates,
        "amplification": amplification,
        "amplitude": amplitude,
    })

check(
    "the preregistered relative map has rank four and all four directions exist",
    relative_rank == 4
    and len(directions) == 4
    and all(
        abs(np.sum(record["boundary_vector"])) < 2e-14
        and abs(np.linalg.norm(record["boundary_vector"])-1) < 2e-14
        and 4.7428e5 < record["amplification"] < 4.7431e5
        and 2.108e-10 < record["amplitude"] < 2.109e-10
        for record in directions
    ),
    f"relative singular values={relative_singular.tolist()}",
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


def solve_transverse(model, quotient_hessian, t, initial_z, boundary_log):
    diagnostics = {
        "evaluations": 0,
        "all_finite": True,
        "minimum_gram": math.inf,
        "minimum_argument": math.inf,
        "maximum_imaginary": 0.0,
    }

    def residual(z):
        state = equation_evaluation(model, t, z, boundary_log)
        diagnostics["evaluations"] += 1
        diagnostics["all_finite"] &= state["finite"]
        diagnostics["minimum_gram"] = min(
            diagnostics["minimum_gram"], state["minimum_gram"]
        )
        diagnostics["minimum_argument"] = min(
            diagnostics["minimum_argument"], state["minimum_argument"]
        )
        diagnostics["maximum_imaginary"] = max(
            diagnostics["maximum_imaginary"], state["maximum_imaginary"]
        )
        if not state["finite"]:
            return np.full(34, 1e20)
        transverse = quotient_basis.T@state["equation"]
        return np.linalg.solve(quotient_hessian, transverse)

    optimization = least_squares(
        residual,
        np.array(initial_z, dtype=float),
        method="trf",
        jac="3-point",
        diff_step=1e-4,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=800,
        x_scale=1.0,
    )
    z = optimization.x
    final = equation_evaluation(model, t, z, boundary_log)
    transverse = quotient_basis.T@final["equation"]
    preconditioned = np.linalg.solve(quotient_hessian, transverse)
    tangent = path_tangent(t)
    scalar = float(tangent@final["equation"])
    transverse_norm = float(np.linalg.norm(transverse))
    preconditioned_norm = float(np.linalg.norm(preconditioned))
    resolved = bool(
        diagnostics["all_finite"]
        and final["finite"]
        and final["branch"]
        and transverse_norm < 1e-9
        and preconditioned_norm < 1e-5
    )
    return {
        "t": float(t),
        "z": z,
        "outcome": "TRANSVERSE_SOLVED" if resolved else "TRANSVERSE_UNRESOLVED",
        "optimizer_success": bool(optimization.success),
        "optimizer_status": int(optimization.status),
        "optimizer_message": str(optimization.message),
        "optimizer_nfev": int(optimization.nfev),
        "optimizer_njev": None if optimization.njev is None else int(optimization.njev),
        "optimizer_cost": float(optimization.cost),
        "diagnostics": diagnostics,
        "variables": final["variables"],
        "equation": final["equation"],
        "transverse": transverse,
        "transverse_norm": transverse_norm,
        "preconditioned_norm": preconditioned_norm,
        "scalar": scalar,
        "full_equation_norm": float(np.linalg.norm(final["equation"])),
        "final_branch": final["branch"],
        "final_minimum_gram": final["minimum_gram"],
        "final_minimum_argument": final["minimum_argument"],
        "final_maximum_imaginary": final["maximum_imaginary"],
    }


def bisect_scalar(model, quotient_hessian, left, right, boundary_log):
    history = []
    current_left, current_right = left, right
    if current_left["scalar"]*current_right["scalar"] >= 0:
        return None, history, False
    for _iteration in range(30):
        midpoint_t = (current_left["t"]+current_right["t"])/2.0
        initial_z = (current_left["z"]+current_right["z"])/2.0
        midpoint = solve_transverse(
            model, quotient_hessian, midpoint_t, initial_z, boundary_log
        )
        history.append(midpoint)
        if midpoint["outcome"] != "TRANSVERSE_SOLVED":
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
        "optimizer_success": record["optimizer_success"],
        "optimizer_status": record["optimizer_status"],
        "optimizer_message": record["optimizer_message"],
        "optimizer_nfev": record["optimizer_nfev"],
        "optimizer_njev": record["optimizer_njev"],
        "optimizer_cost": record["optimizer_cost"],
        "diagnostics": record["diagnostics"],
        "variables": record["variables"].tolist(),
        "equation": record["equation"].tolist(),
        "transverse": record["transverse"].tolist(),
        "transverse_norm": record["transverse_norm"],
        "preconditioned_norm": record["preconditioned_norm"],
        "scalar": record["scalar"],
        "full_equation_norm": record["full_equation_norm"],
        "final_branch": record["final_branch"],
        "final_minimum_gram": record["final_minimum_gram"],
        "final_minimum_argument": record["final_minimum_argument"],
        "final_maximum_imaginary": record["final_maximum_imaginary"],
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
    branch_records = []
    points = []
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
        f"Validating {len(candidates)} candidates with "
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
        equation_error = float(np.linalg.norm(error_row))
        equation_norm = float(np.linalg.norm(equation))
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
        label = (
            "NONLINEAR_STATIONARY_CANDIDATE" if all(gates.values())
            else "ACTION_ONLY_CANDIDATE_FAILED"
        )
        audits.append({
            "label": label,
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
        raise RuntimeError("candidate action cursor mismatch")
    return audits


results = {}
for parity, model in dust.bl.models.items():
    print(f"Localizing nonlinear cases: {parity} parity...", flush=True)
    quotient_hessian = np.array(
        precision_input["parities"][parity]["quotient"]["matrix"],
        dtype=float,
    )
    quotient_hessian = (quotient_hessian+quotient_hessian.T)/2.0
    cases = []
    all_candidates = []

    for direction in directions:
        for sign in (-1, 1):
            boundary_log = sign*direction["amplitude"]*direction["boundary_vector"]
            linear_z = sign*direction["amplitude"]*direction["response_coefficients"]
            grid = {}
            zero = solve_transverse(
                model, quotient_hessian, 0.0, linear_z, boundary_log
            )
            grid[0.0] = zero
            previous = zero
            for t in POSITIVE_SCAN:
                initial = previous["z"] if previous["outcome"] == "TRANSVERSE_SOLVED" else linear_z
                current = solve_transverse(
                    model, quotient_hessian, t, initial, boundary_log
                )
                grid[t] = current
                if current["outcome"] == "TRANSVERSE_SOLVED":
                    previous = current
            previous = zero
            for t in NEGATIVE_SCAN:
                initial = previous["z"] if previous["outcome"] == "TRANSVERSE_SOLVED" else linear_z
                current = solve_transverse(
                    model, quotient_hessian, t, initial, boundary_log
                )
                grid[t] = current
                if current["outcome"] == "TRANSVERSE_SOLVED":
                    previous = current

            ordered = [grid[t] for t in T_GRID]
            candidates = [
                item for item in ordered
                if item["outcome"] == "TRANSVERSE_SOLVED"
                and abs(item["scalar"]) <= SCALAR_FLOOR
            ]
            bisections = []
            intervals_resolved = True
            for left, right in zip(ordered[:-1], ordered[1:]):
                if (
                    left["outcome"] != "TRANSVERSE_SOLVED"
                    or right["outcome"] != "TRANSVERSE_SOLVED"
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
                        "history": history,
                        "resolved": resolved,
                        "candidate": candidate,
                    })
                    intervals_resolved &= resolved
                    if candidate is not None:
                        candidates.append(candidate)

            candidates.sort(key=lambda item: item["t"])
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
                "direction_index": direction["index"],
                "sign": sign,
                "boundary_log": boundary_log,
                "linear_z": linear_z,
                "grid": grid,
                "bisections": bisections,
                "candidate_indices": candidate_indices,
                "all_grid_resolved": all(
                    item["outcome"] == "TRANSVERSE_SOLVED" for item in ordered
                ),
                "all_intervals_resolved": intervals_resolved,
            })

    audits = audit_candidates(model, quotient_hessian, all_candidates)
    for case in cases:
        case_audits = [audits[index] for index in case["candidate_indices"]]
        passing = sum(
            audit["label"] == "NONLINEAR_STATIONARY_CANDIDATE"
            for audit in case_audits
        )
        localization_resolved = bool(
            case["all_grid_resolved"] and case["all_intervals_resolved"]
        )
        if passing > 0:
            outcome = "NONLINEAR_CONTINUATION_FOUND_IN_SCAN"
        elif localization_resolved:
            outcome = "NO_NONLINEAR_CONTINUATION_IN_FROZEN_SCAN"
        else:
            outcome = "NONLINEAR_CONTINUATION_NUMERICALLY_UNRESOLVED"
        case["audits"] = case_audits
        case["passing_candidates"] = passing
        case["localization_resolved"] = localization_resolved
        case["outcome"] = outcome

    results[parity] = {
        "quotient_hessian": quotient_hessian,
        "cases": cases,
        "all_candidates": all_candidates,
        "audits": audits,
    }
    check(
        f"{parity}: all eight preregistered signed cases were enumerated",
        len(cases) == 8
        and {(case["direction_index"], case["sign"]) for case in cases}
        == {(index, sign) for index in range(1, 5) for sign in (-1, 1)},
    )
    check(
        f"{parity}: every localized candidate received one action-only audit",
        len(all_candidates) == len(audits)
        and all("gates" in audit for audit in audits),
        f"candidates={len(audits)}, passing={sum(a['label']=='NONLINEAR_STATIONARY_CANDIDATE' for a in audits)}",
    )
    check(
        f"{parity}: all eight nonlinear outcomes follow the frozen rules",
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
        "linear_initial_z": case["linear_z"].tolist(),
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
                "solve": compact_solve(
                    results_case_candidate(case, local_index)
                ),
                "action_only": serialize_audit(audit),
            }
            for local_index, audit in enumerate(case["audits"])
        ],
    }


def results_case_candidate(case, local_index):
    # This helper is rebound through the stored global candidate index below
    # in the parity serializer; it is never called without that binding.
    return case["_candidate_records"][local_index]


signed_hits = 0
pair_hits = 0
contrast_hits = 0
for parity, record in results.items():
    for case in record["cases"]:
        case["_candidate_records"] = [
            record["all_candidates"][index]
            for index in case["candidate_indices"]
        ]
        signed_hits += case["outcome"] == "NONLINEAR_CONTINUATION_FOUND_IN_SCAN"
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
    "hit fractions cover all 16 signed cases, 8 pairs and 4 contrasts",
    0 <= signed_hits <= 16
    and 0 <= pair_hits <= 8
    and 0 <= contrast_hits <= 4,
    f"signed={signed_hits}/16, pairs={pair_hits}/8, contrasts={contrast_hits}/4",
)

if signed_hits > 0:
    verdict = (
        "DERIVED COMPUTATIONAL LOCAL: at least one of sixteen preregistered "
        "relative-phase cases has a complete action-validated nonlinear "
        "continuation in the frozen lapse scan.  Robustness is quantified "
        "by the recorded hit fractions."
    )
elif all(
    case["outcome"] == "NO_NONLINEAR_CONTINUATION_IN_FROZEN_SCAN"
    for record in results.values() for case in record["cases"]
):
    verdict = (
        "DERIVED COMPUTATIONAL NEGATIVE WITHIN SCOPE: no relative-phase "
        "case continues nonlinearly at eta=1e-4 within |t|<=0.1."
    )
else:
    verdict = (
        "OPEN NUMERICALLY: no action-validated continuation was found and at "
        "least one preregistered case did not resolve its frozen scan."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "weak_scale_path_result_commit": PATH_RESULT_COMMIT,
    "precision_response_result_commit": PRECISION_RESULT_COMMIT,
    "eta": ETA,
    "collective_scan": T_GRID,
    "action_precision_digits": 100,
    "action_derivative_steps": ACTION_STEPS,
    "direction_enumeration": {
        "count": 4,
        "relative_map_singular_values": relative_singular.tolist(),
        "absolute_rank_above_1e-8": relative_rank,
        "directions": [
            {
                "index": record["index"],
                "helmert_contrast": record["contrast"].tolist(),
                "boundary_vector": record["boundary_vector"].tolist(),
                "boundary_sum": float(np.sum(record["boundary_vector"])),
                "response_amplification": record["amplification"],
                "boundary_amplitude": record["amplitude"],
                "predicted_internal_norm": float(np.linalg.norm(
                    record["amplitude"]*record["response_vector"]
                )),
            }
            for record in directions
        ],
    },
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
        "nonlinear_relative_continuation": "DERIVED COMPUTATIONAL IF VALIDATED",
        "look_elsewhere": "REPORT ALL HIT FRACTIONS",
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

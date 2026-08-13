#!/usr/bin/env python3
"""Audit the mechanism behind frozen-Hessian nonlinear stagnation.

Protocol commit: e53dcaf.  This verifier does not iterate to a root.  It
computes displaced local Jacobians at all 80 frozen states, tries exactly one
local-Jacobian step where its inverse action is resolved, and independently
audits 48 deterministic states from the complete action at 100 decimals.
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
DEFECT_INPUT = HERE / "gravity_600cell_dust_nonlinear_defect.json"
OUTPUT = HERE / "gravity_600cell_dust_stagnation.json"
PROTOCOL_COMMIT = "e53dcaf"
PRIOR_ART_COMMIT = "53968b0"
DEFECT_RESULT_COMMIT = "8f651cd"
JACOBIAN_STEPS = (5.0e-4, 2.5e-4, 1.25e-4)
ACTION_STEPS = (2.0e-5, 1.0e-5, 5.0e-6)
DAMPING = tuple(2.0**(-power) for power in range(11))
T_GRID = (-0.10, -0.05, 0.0, 0.05, 0.10)
ACTION_T = (-0.10, 0.0, 0.10)
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
defect_input = json.loads(DEFECT_INPUT.read_text())
check(
    "all three frozen inputs retain their registered provenance",
    precision_input["protocol_commit"] == "da34272"
    and first_input["protocol_commit"] == "80f8de7"
    and defect_input["protocol_commit"] == "2d695c6"
    and precision_input["passed"] == precision_input["tests"] == 14
    and first_input["passed"] == first_input["tests"] == 10
    and defect_input["passed"] == defect_input["tests"] == 10
    and defect_input["hit_fractions"]["signed_cases"]
    == {"hits": 0, "total": 16},
)


print("Loading the independently certified complete-action control...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_stagnation",
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
check(
    "the frozen complement has the expected orthonormal 35 by 34 shape",
    quotient_basis.shape == (35, 34)
    and np.linalg.norm(quotient_basis.T@quotient_basis-np.eye(34), 2) < 1e-12,
)


def path_internal(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    return np.concatenate((np.full(30, diagonal), np.full(5, rho)))


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
            "transverse": quotient_basis.T@equation,
            "maximum_imaginary": float(
                np.max(np.abs(equation_complex.imag))
            ),
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
            "transverse": np.full(34, np.nan),
            "maximum_imaginary": math.inf,
            "minimum_gram": -math.inf,
            "minimum_argument": -math.inf,
            "lorentzian": False,
            "branch": False,
            "finite": False,
            "error": repr(error),
        }


def local_jacobian_audit(
    model, quotient_hessian, background_soft, t, z, boundary_log
):
    base = equation_evaluation(model, t, z, boundary_log)
    matrices = []
    displaced_branch = bool(base["branch"] and base["finite"])
    minimum_gram = base["minimum_gram"]
    minimum_argument = base["minimum_argument"]
    maximum_imaginary = base["maximum_imaginary"]
    displaced_count = 0

    for step in JACOBIAN_STEPS:
        matrix = np.zeros((34, 34))
        for coordinate in range(34):
            displacement = np.zeros(34)
            displacement[coordinate] = step
            plus = equation_evaluation(
                model, t, z+displacement, boundary_log
            )
            minus = equation_evaluation(
                model, t, z-displacement, boundary_log
            )
            displaced_count += 2
            for state in (plus, minus):
                displaced_branch &= bool(state["branch"] and state["finite"])
                minimum_gram = min(minimum_gram, state["minimum_gram"])
                minimum_argument = min(
                    minimum_argument, state["minimum_argument"]
                )
                maximum_imaginary = max(
                    maximum_imaginary, state["maximum_imaginary"]
                )
            matrix[:, coordinate] = (
                plus["transverse"]-minus["transverse"]
            )/(2.0*step)
        matrices.append(matrix)

    j1, j2, j3 = matrices
    r12 = (4.0*j2-j1)/3.0
    r23 = (4.0*j3-j2)/3.0
    j6_raw = (16.0*r23-r12)/15.0
    delta_j = j6_raw-r23
    symmetry_defect = float(
        np.linalg.norm((j6_raw-j6_raw.T)/2.0, 2)
    )
    error_norm = float(np.linalg.norm(delta_j, 2))
    jacobian = (j6_raw+j6_raw.T)/2.0
    jacobian23 = (r23+r23.T)/2.0
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    singular_values23 = np.linalg.svd(jacobian23, compute_uv=False)
    eigenvalues, eigenvectors = np.linalg.eigh(jacobian)
    soft_indices = np.argsort(np.abs(eigenvalues))[:4]
    local_soft = eigenvectors[:, soft_indices]
    principal_cosines = np.linalg.svd(
        background_soft.T@local_soft, compute_uv=False
    )

    residual = base["transverse"]
    residual_norm = float(np.linalg.norm(residual))
    fixed_step = np.linalg.solve(quotient_hessian, residual)
    eta_fixed = float(np.linalg.norm(
        residual-jacobian@fixed_step
    )/max(residual_norm, 1e-30))
    mu_fixed = float(
        residual@(jacobian@fixed_step)/max(residual_norm**2, 1e-60)
    )
    delta_mu_fixed = float(
        np.linalg.norm(delta_j@fixed_step)/max(residual_norm, 1e-30)
    )

    inverse_available = bool(
        singular_values[-1] > 1e-12
        and singular_values23[-1] > 1e-12
    )
    local_step = None
    local_step23 = None
    step_change = None
    step_model_error = None
    if inverse_available:
        local_step = np.linalg.solve(jacobian, residual)
        local_step23 = np.linalg.solve(jacobian23, residual)
        step_change = float(
            np.linalg.norm(local_step-local_step23)
            / max(np.linalg.norm(local_step), 1e-30)
        )
        step_model_error = float(
            np.linalg.norm(delta_j@local_step)/max(residual_norm, 1e-30)
        )

    symmetry_resolved = bool(
        symmetry_defect <= 10.0*max(error_norm, 1e-30)
    )
    jacobian_resolved = bool(
        displaced_branch
        and symmetry_resolved
        and inverse_available
        and step_change is not None
        and step_model_error is not None
        and step_change <= 0.1
        and step_model_error <= 0.1
    )

    if not jacobian_resolved:
        fixed_direction = "FIXED_DIRECTION_NOT_CLASSIFIED"
    elif mu_fixed+10.0*delta_mu_fixed < 0.0:
        fixed_direction = "ROBUST_FIXED_NONDESCENT"
    elif mu_fixed-10.0*delta_mu_fixed > 0.0:
        fixed_direction = "ROBUST_FIXED_DESCENT"
    else:
        fixed_direction = "FIXED_DIRECTION_UNRESOLVED"

    local_trials = []
    accepted = None
    if jacobian_resolved:
        for damping in DAMPING:
            trial = equation_evaluation(
                model, t, z-damping*local_step, boundary_log
            )
            trial_norm = float(np.linalg.norm(trial["transverse"]))
            local_trials.append({
                "damping": damping,
                "branch": bool(trial["branch"] and trial["finite"]),
                "transverse_norm": trial_norm,
                "minimum_gram": trial["minimum_gram"],
                "minimum_argument": trial["minimum_argument"],
                "maximum_imaginary": trial["maximum_imaginary"],
            })
            if (
                trial["branch"] and trial["finite"]
                and trial_norm < residual_norm
            ):
                accepted = {
                    "damping": damping,
                    "transverse_norm": trial_norm,
                    "reduction_factor": trial_norm/max(residual_norm, 1e-30),
                    "z": z-damping*local_step,
                    "variables": trial["variables"],
                    "equation": trial["equation"],
                }
                break

    strong_local_descent = bool(
        accepted is not None and accepted["reduction_factor"] < 0.5
    )
    if not jacobian_resolved:
        local_descent = "LOCAL_JACOBIAN_NOT_TRIED"
        mechanism = "PRECISION_LIMITED_LOCAL_MODEL"
    elif strong_local_descent:
        local_descent = "LOCAL_JACOBIAN_DESCENT"
        if fixed_direction == "ROBUST_FIXED_NONDESCENT":
            mechanism = "FIXED_MODEL_MISMATCH"
        else:
            mechanism = "RESOLVED_BUT_NOT_FIXED_MISMATCH"
    else:
        local_descent = "NO_STRONG_LOCAL_JACOBIAN_DESCENT"
        mechanism = "RESOLVED_BUT_NOT_FIXED_MISMATCH"

    return {
        "t": float(t),
        "z": z,
        "base_variables": base["variables"],
        "base_equation": base["equation"],
        "base_transverse": residual,
        "base_transverse_norm": residual_norm,
        "jacobian_steps": matrices,
        "richardson12": r12,
        "richardson23": r23,
        "jacobian6_raw": j6_raw,
        "jacobian": jacobian,
        "jacobian_error_matrix": delta_j,
        "jacobian_error_norm": error_norm,
        "symmetry_defect": symmetry_defect,
        "symmetry_resolved": symmetry_resolved,
        "singular_values": singular_values,
        "singular_values23": singular_values23,
        "eigenvalues": eigenvalues,
        "soft_indices": soft_indices,
        "principal_cosines": principal_cosines,
        "fixed_step": fixed_step,
        "fixed_forcing_factor": eta_fixed,
        "fixed_descent_mu": mu_fixed,
        "fixed_descent_uncertainty": delta_mu_fixed,
        "fixed_direction": fixed_direction,
        "local_step": local_step,
        "local_step23": local_step23,
        "step_change": step_change,
        "step_model_error": step_model_error,
        "jacobian_label": (
            "JACOBIAN_STEP_RESOLVED" if jacobian_resolved
            else "LOCAL_JACOBIAN_UNRESOLVED"
        ),
        "local_trials": local_trials,
        "accepted_local_trial": accepted,
        "local_descent": local_descent,
        "mechanism": mechanism,
        "displaced_count": displaced_count,
        "all_displaced_branches_pass": displaced_branch,
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
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
_ARB_ANCHORS = None


def initialize_action_worker(model, anchor_strings):
    global _ARB_MODEL, _ARB_ANCHORS
    _ARB_MODEL = model
    configure_arb_precision()
    _ARB_ANCHORS = [
        [arb.mpf(value) for value in anchor]
        for anchor in anchor_strings
    ]


def action_worker(point):
    anchor_index, coordinate, logarithmic_step = point
    variables = list(_ARB_ANCHORS[anchor_index])
    variables[coordinate] *= arb.exp(arb.mpf(str(logarithmic_step)))
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


def audit_action_anchors(model, anchors):
    points = []
    branch_records = []
    for anchor_index, anchor in enumerate(anchors):
        variables = anchor["variables"]
        branch_pass = True
        minimum_gram = math.inf
        minimum_argument = math.inf
        for step in ACTION_STEPS:
            for coordinate in range(35):
                for sign in (1.0, -1.0):
                    logarithmic_step = sign*step
                    points.append((
                        anchor_index, coordinate, logarithmic_step
                    ))
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

    anchor_strings = [
        [f"{value:.17g}" for value in anchor["variables"]]
        for anchor in anchors
    ]
    print(
        f"Validating {len(anchors)} deterministic action anchors with "
        f"{len(points)} action points...",
        flush=True,
    )
    fork_context = mp.get_context("fork")
    with fork_context.Pool(
        processes=8,
        initializer=initialize_action_worker,
        initargs=(model, anchor_strings),
    ) as pool:
        action_values = pool.map(action_worker, points, chunksize=1)

    audits = []
    cursor = 0
    for anchor_index, anchor in enumerate(anchors):
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
        action_transverse = quotient_basis.T@equation
        error_transverse = quotient_basis.T@error_row
        action_norm = float(np.linalg.norm(action_transverse))
        error_norm = float(np.linalg.norm(error_transverse))
        floor = max(error_norm, 1e-30)
        if action_norm <= 10.0*floor:
            label = "ACTION_ZERO_CONSISTENT"
        elif action_norm > 100.0*floor:
            label = "ACTION_RESOLVED_NONZERO"
        else:
            label = "ACTION_RESIDUAL_UNRESOLVED"
        binary_discrepancy = float(np.linalg.norm(
            anchor["equation"]-equation
        ))
        transverse_discrepancy = float(np.linalg.norm(
            anchor["transverse"]-action_transverse
        ))
        maximum_imaginary = float(max(
            maximum_imaginary,
            max(abs(arb.im(value)) for value in sixth_order),
        ))
        audits.append({
            "parity": anchor["parity"],
            "direction_index": anchor["direction_index"],
            "sign": anchor["sign"],
            "t": anchor["t"],
            "label": label,
            "derivatives": derivatives,
            "richardson12": richardson12,
            "richardson23": richardson23,
            "sixth_order": sixth_order,
            "equation": equation,
            "error_row": error_row,
            "action_transverse": action_transverse,
            "action_transverse_error": error_transverse,
            "action_transverse_norm": action_norm,
            "action_transverse_error_norm": error_norm,
            "binary_equation_discrepancy": binary_discrepancy,
            "binary_transverse_discrepancy": transverse_discrepancy,
            "discrepancy_over_action_error": (
                transverse_discrepancy/floor
            ),
            "analytic_action_agreement_below_1e-8": (
                binary_discrepancy < 1e-8
            ),
            "maximum_imaginary": maximum_imaginary,
            "imaginary_below_1e-80": maximum_imaginary < 1e-80,
            "branch": branch_records[anchor_index],
        })
    if cursor != len(action_values):
        raise RuntimeError("stagnation action cursor mismatch")
    return audits


def array_list(value):
    return np.asarray(value).tolist()


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 80),
        "imaginary": arb.nstr(arb.im(value), 80),
    }


def serialize_local(record):
    accepted = record["accepted_local_trial"]
    return {
        "t": record["t"],
        "z": array_list(record["z"]),
        "base_variables": array_list(record["base_variables"]),
        "base_equation": array_list(record["base_equation"]),
        "base_transverse": array_list(record["base_transverse"]),
        "base_transverse_norm": record["base_transverse_norm"],
        "jacobian_steps": [
            array_list(matrix) for matrix in record["jacobian_steps"]
        ],
        "richardson12": array_list(record["richardson12"]),
        "richardson23": array_list(record["richardson23"]),
        "jacobian6_raw": array_list(record["jacobian6_raw"]),
        "jacobian": array_list(record["jacobian"]),
        "jacobian_error_matrix": array_list(
            record["jacobian_error_matrix"]
        ),
        "jacobian_error_norm": record["jacobian_error_norm"],
        "symmetry_defect": record["symmetry_defect"],
        "symmetry_resolved": record["symmetry_resolved"],
        "singular_values": array_list(record["singular_values"]),
        "singular_values_richardson23": array_list(
            record["singular_values23"]
        ),
        "eigenvalues": array_list(record["eigenvalues"]),
        "soft_indices": array_list(record["soft_indices"]),
        "soft_principal_cosines_with_background": array_list(
            record["principal_cosines"]
        ),
        "fixed_step": array_list(record["fixed_step"]),
        "fixed_forcing_factor": record["fixed_forcing_factor"],
        "fixed_descent_mu": record["fixed_descent_mu"],
        "fixed_descent_uncertainty": record["fixed_descent_uncertainty"],
        "fixed_direction": record["fixed_direction"],
        "local_step": array_list(record["local_step"]),
        "local_step_richardson23": array_list(record["local_step23"]),
        "step_change": record["step_change"],
        "step_model_error": record["step_model_error"],
        "jacobian_label": record["jacobian_label"],
        "local_trials": record["local_trials"],
        "accepted_local_trial": None if accepted is None else {
            "damping": accepted["damping"],
            "transverse_norm": accepted["transverse_norm"],
            "reduction_factor": accepted["reduction_factor"],
            "z": array_list(accepted["z"]),
            "variables": array_list(accepted["variables"]),
            "equation": array_list(accepted["equation"]),
        },
        "local_descent": record["local_descent"],
        "mechanism": record["mechanism"],
        "displaced_count": record["displaced_count"],
        "all_displaced_branches_pass": record[
            "all_displaced_branches_pass"
        ],
        "minimum_gram": record["minimum_gram"],
        "minimum_argument": record["minimum_argument"],
        "maximum_imaginary": record["maximum_imaginary"],
    }


def serialize_action(record):
    return {
        "parity": record["parity"],
        "direction_index": record["direction_index"],
        "sign": record["sign"],
        "t": record["t"],
        "label": record["label"],
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
        "sixth_order_float": array_list(record["equation"]),
        "error_row_float": array_list(record["error_row"]),
        "action_transverse": array_list(record["action_transverse"]),
        "action_transverse_error": array_list(
            record["action_transverse_error"]
        ),
        "action_transverse_norm": record["action_transverse_norm"],
        "action_transverse_error_norm": record[
            "action_transverse_error_norm"
        ],
        "binary_equation_discrepancy": record[
            "binary_equation_discrepancy"
        ],
        "binary_transverse_discrepancy": record[
            "binary_transverse_discrepancy"
        ],
        "discrepancy_over_action_error": record[
            "discrepancy_over_action_error"
        ],
        "analytic_action_agreement_below_1e-8": record[
            "analytic_action_agreement_below_1e-8"
        ],
        "maximum_imaginary": record["maximum_imaginary"],
        "imaginary_below_1e-80": record["imaginary_below_1e-80"],
        "branch": record["branch"],
    }


results = {}
all_local = []
all_anchors = []
maximum_stored_equation_difference = 0.0
for parity, model in dust.bl.models.items():
    print(f"Computing all 40 local Jacobians: {parity} parity...", flush=True)
    quotient_hessian = np.array(
        precision_input["parities"][parity]["quotient"]["matrix"],
        dtype=float,
    )
    quotient_hessian = (quotient_hessian+quotient_hessian.T)/2.0
    h_eigenvalues, h_eigenvectors = np.linalg.eigh(quotient_hessian)
    h_soft_indices = np.argsort(np.abs(h_eigenvalues))[:4]
    background_soft = h_eigenvectors[:, h_soft_indices]
    cases = []
    defect_cases = defect_input["parities"][parity]["cases"]
    for case_number, defect_case in enumerate(defect_cases, 1):
        boundary_log = np.array(
            defect_case["boundary_log_vector"], dtype=float
        )
        states = {}
        for t in T_GRID:
            key = f"{t:+.2f}"
            stored = defect_case["grid"][key]
            z = np.array(stored["z"], dtype=float)
            record = local_jacobian_audit(
                model, quotient_hessian, background_soft,
                t, z, boundary_log,
            )
            maximum_stored_equation_difference = max(
                maximum_stored_equation_difference,
                float(np.linalg.norm(
                    record["base_equation"]
                    - np.array(stored["equation"], dtype=float)
                )),
            )
            record["previous_stop_reason"] = stored["stop_reason"]
            states[t] = record
            all_local.append(record)
            if t in ACTION_T:
                all_anchors.append({
                    "parity": parity,
                    "direction_index": int(defect_case["direction_index"]),
                    "sign": int(defect_case["sign"]),
                    "t": float(t),
                    "variables": record["base_variables"],
                    "equation": record["base_equation"],
                    "transverse": record["base_transverse"],
                })
        cases.append({
            "direction_index": int(defect_case["direction_index"]),
            "sign": int(defect_case["sign"]),
            "boundary_log_vector": boundary_log,
            "states": states,
        })
        print(
            f"  {parity}: completed case {case_number}/8",
            flush=True,
        )
    results[parity] = {
        "quotient_hessian": quotient_hessian,
        "background_eigenvalues": h_eigenvalues,
        "background_soft_indices": h_soft_indices,
        "cases": cases,
    }
    parity_records = [
        state for case in cases for state in case["states"].values()
    ]
    check(
        f"{parity}: all forty frozen states have complete local audits",
        len(parity_records) == 40
        and all(record["displaced_count"] == 204 for record in parity_records),
    )
    check(
        f"{parity}: every local state received exactly one mechanism label",
        all(record["mechanism"] in {
            "FIXED_MODEL_MISMATCH",
            "RESOLVED_BUT_NOT_FIXED_MISMATCH",
            "PRECISION_LIMITED_LOCAL_MODEL",
        } for record in parity_records),
    )


check(
    "all eighty recomputed base equations match their stored states",
    maximum_stored_equation_difference < 1e-12,
    f"maximum difference={maximum_stored_equation_difference:.3e}",
)
check(
    "all 16 cases contribute centre and both endpoints to 48 anchors",
    len(all_anchors) == 48
    and {
        (anchor["parity"], anchor["direction_index"], anchor["sign"], anchor["t"])
        for anchor in all_anchors
    } == {
        (parity, direction, sign, t)
        for parity in ("even", "odd")
        for direction in range(1, 5)
        for sign in (-1, 1)
        for t in ACTION_T
    },
)


all_action_audits = []
for parity, model in dust.bl.models.items():
    anchors = [
        anchor for anchor in all_anchors if anchor["parity"] == parity
    ]
    all_action_audits.extend(audit_action_anchors(model, anchors))

check(
    "all 48 complete-action anchors retain all 210 branch points",
    len(all_action_audits) == 48
    and all(audit["branch"]["pass"] for audit in all_action_audits),
)
check(
    "all complete-action anchors retain the frozen reality gate",
    all(audit["imaginary_below_1e-80"] for audit in all_action_audits),
    f"maximum={max(audit['maximum_imaginary'] for audit in all_action_audits):.3e}",
)
check(
    "all binary analytic rows agree with action rows below 1e-8",
    all(
        audit["analytic_action_agreement_below_1e-8"]
        for audit in all_action_audits
    ),
    f"maximum={max(audit['binary_equation_discrepancy'] for audit in all_action_audits):.3e}",
)


resolved_count = sum(
    record["jacobian_label"] == "JACOBIAN_STEP_RESOLVED"
    for record in all_local
)
mismatch_count = sum(
    record["mechanism"] == "FIXED_MODEL_MISMATCH"
    for record in all_local
)
if resolved_count >= 60 and mismatch_count >= 0.75*resolved_count:
    aggregate_outcome = "FIXED_MODEL_MISMATCH_DOMINANT"
elif resolved_count < 60:
    aggregate_outcome = "LOCAL_MODEL_PRECISION_LIMITED"
else:
    aggregate_outcome = "MIXED_OR_OTHER_STAGNATION"

mechanism_counts = {
    label: sum(record["mechanism"] == label for record in all_local)
    for label in (
        "FIXED_MODEL_MISMATCH",
        "RESOLVED_BUT_NOT_FIXED_MISMATCH",
        "PRECISION_LIMITED_LOCAL_MODEL",
    )
}
action_counts = {
    label: sum(audit["label"] == label for audit in all_action_audits)
    for label in (
        "ACTION_ZERO_CONSISTENT",
        "ACTION_RESOLVED_NONZERO",
        "ACTION_RESIDUAL_UNRESOLVED",
    )
}
check(
    "the aggregate outcome follows the frozen 60/80 and 3/4 rules",
    aggregate_outcome in {
        "FIXED_MODEL_MISMATCH_DOMINANT",
        "LOCAL_MODEL_PRECISION_LIMITED",
        "MIXED_OR_OTHER_STAGNATION",
    }
    and sum(mechanism_counts.values()) == 80
    and sum(action_counts.values()) == 48,
    f"resolved={resolved_count}/80, mechanisms={mechanism_counts}, action={action_counts}",
)


payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "frozen_defect_result_commit": DEFECT_RESULT_COMMIT,
    "jacobian_steps": JACOBIAN_STEPS,
    "action_steps": ACTION_STEPS,
    "action_precision_digits": 100,
    "action_anchor_t": ACTION_T,
    "maximum_stored_equation_difference": maximum_stored_equation_difference,
    "parities": {
        parity: {
            "quotient_hessian": array_list(record["quotient_hessian"]),
            "background_eigenvalues": array_list(
                record["background_eigenvalues"]
            ),
            "background_soft_indices": array_list(
                record["background_soft_indices"]
            ),
            "cases": [
                {
                    "direction_index": case["direction_index"],
                    "sign": case["sign"],
                    "boundary_log_vector": array_list(
                        case["boundary_log_vector"]
                    ),
                    "states": {
                        f"{t:+.2f}": {
                            **serialize_local(state),
                            "previous_stop_reason": state[
                                "previous_stop_reason"
                            ],
                        }
                        for t, state in sorted(case["states"].items())
                    },
                }
                for case in record["cases"]
            ],
        }
        for parity, record in results.items()
    },
    "action_anchors": [
        serialize_action(record) for record in all_action_audits
    ],
    "resolved_local_jacobians": {
        "count": int(resolved_count), "total": 80
    },
    "mechanism_counts": mechanism_counts,
    "action_residual_counts": action_counts,
    "aggregate_outcome": aggregate_outcome,
    "claim_boundary": {
        "finite_diagnostics": "DERIVED COMPUTATIONAL IF CHECKS PASS",
        "method_provenance": "PATTERN-INFORMED POST-RESULT DIAGNOSTIC",
        "stationary_root": "OPEN",
        "root_nonexistence": "OPEN AND NOT TESTED",
        "second_slab": "NOT TESTED",
        "full_840_edge_carrier": "NOT TESTED",
        "clock_speed_limit_planck_scale": "NOT DERIVED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(
    f"{aggregate_outcome}: resolved local Jacobians {resolved_count}/80; "
    f"mechanisms {mechanism_counts}; action anchors {action_counts}."
)
raise SystemExit(0 if passed == tests else 1)

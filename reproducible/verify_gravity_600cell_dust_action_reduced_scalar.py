#!/usr/bin/env python3
"""Complete-action transverse solve and Lyapunov--Schmidt scalar audit.

Protocol commit: 17f9560.  The physical residual is reconstructed only from
the complete action.  Binary local Jacobians are recomputed at accepted
iterates and used solely as preconditioners.  Solver and validation action
windows are disjoint.
"""

import contextlib
import importlib.util
import io
import json
import math
import multiprocessing as mp
from pathlib import Path
import pickle
import sys

import mpmath as arb
import numpy as np


HERE = Path(__file__).resolve().parent
PRECISION_INPUT = HERE / "gravity_600cell_dust_gauge_quotient_precision.json"
DEFECT_INPUT = HERE / "gravity_600cell_dust_nonlinear_defect.json"
STAGNATION_INPUT = HERE / "gravity_600cell_dust_stagnation.json"
OUTPUT = HERE / "gravity_600cell_dust_action_reduced_scalar.json"
CHECKPOINT = HERE / ".gravity_600cell_dust_action_reduced_scalar.checkpoint.pkl"
PROTOCOL_COMMIT = "17f9560"
PRIOR_ART_COMMIT = "8be130e"
STAGNATION_RESULT_COMMIT = "1d66278"
SOLVER_STEPS = (6.0e-5, 3.0e-5, 1.5e-5)
VALIDATION_STEPS = (2.0e-5, 1.0e-5, 5.0e-6)
JACOBIAN_STEPS = (5.0e-4, 2.5e-4, 1.25e-4)
DAMPING = tuple(2.0**(-power) for power in range(11))
T_GRID = (-0.10, -0.05, 0.0, 0.05, 0.10)
MAX_ACCEPTED_ITERATIONS = 6
tests = passed = 0


def checkpoint_metadata():
    """Bind an operational checkpoint to the frozen numerical protocol."""
    return {
        "schema": 1,
        "protocol_commit": PROTOCOL_COMMIT,
        "prior_art_commit": PRIOR_ART_COMMIT,
        "stagnation_result_commit": STAGNATION_RESULT_COMMIT,
        "solver_steps": SOLVER_STEPS,
        "validation_steps": VALIDATION_STEPS,
        "jacobian_steps": JACOBIAN_STEPS,
        "damping": DAMPING,
        "t_grid": T_GRID,
        "max_accepted_iterations": MAX_ACCEPTED_ITERATIONS,
    }


def load_checkpoint():
    if not CHECKPOINT.exists():
        return {}
    with CHECKPOINT.open("rb") as stream:
        payload = pickle.load(stream)
    if payload.get("metadata") != checkpoint_metadata():
        raise RuntimeError(
            f"incompatible operational checkpoint: {CHECKPOINT}"
        )
    completed = payload.get("completed_parities")
    if not isinstance(completed, dict):
        raise RuntimeError(f"malformed operational checkpoint: {CHECKPOINT}")
    print(
        "Loaded operational checkpoint for "
        f"{sorted(completed)}; frozen mathematics unchanged.",
        flush=True,
    )
    return completed


def save_checkpoint(completed_parities):
    payload = {
        "metadata": checkpoint_metadata(),
        "completed_parities": completed_parities,
    }
    temporary = CHECKPOINT.with_suffix(CHECKPOINT.suffix+".tmp")
    with temporary.open("wb") as stream:
        pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)
        stream.flush()
    temporary.replace(CHECKPOINT)
    print(
        "Saved operational checkpoint for "
        f"{sorted(completed_parities)}.",
        flush=True,
    )


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


precision_input = json.loads(PRECISION_INPUT.read_text())
defect_input = json.loads(DEFECT_INPUT.read_text())
stagnation_input = json.loads(STAGNATION_INPUT.read_text())
check(
    "all frozen inputs retain their registered provenance",
    precision_input["protocol_commit"] == "da34272"
    and defect_input["protocol_commit"] == "2d695c6"
    and stagnation_input["protocol_commit"] == "e53dcaf"
    and precision_input["passed"] == precision_input["tests"] == 14
    and defect_input["passed"] == defect_input["tests"] == 10
    and stagnation_input["passed"] == stagnation_input["tests"] == 13
    and stagnation_input["aggregate_outcome"]
    == "MIXED_OR_OTHER_STAGNATION",
)


print("Loading the independently certified complete-action control...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_action_reduced_scalar",
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
    "the exact complement retains its orthonormal 35 by 34 carrier",
    quotient_basis.shape == (35, 34)
    and np.linalg.norm(quotient_basis.T@quotient_basis-np.eye(34), 2) < 1e-12,
)


def path_internal(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    return np.concatenate((np.full(30, diagonal), np.full(5, rho)))


def path_tangent(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    raw = np.concatenate((
        np.full(30, -rho/diagonal), np.ones(5)
    ))
    return raw/np.linalg.norm(raw)


def complete_variables(t, z, boundary_log):
    internal = path_internal(t)*np.exp(quotient_basis@z)
    final = dust.L0_SQUARE*np.exp(boundary_log)
    return np.concatenate((internal, final))


def binary_equation_evaluation(model, t, z, boundary_log):
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
            "branch": branch,
            "finite": finite,
        }
    except Exception as error:
        return {
            "variables": variables,
            "equation": None,
            "transverse": None,
            "maximum_imaginary": None,
            "minimum_gram": None,
            "minimum_argument": None,
            "branch": False,
            "finite": False,
            "error": repr(error),
        }


def local_jacobian(model, t, z, boundary_log):
    base = binary_equation_evaluation(model, t, z, boundary_log)
    matrices = []
    all_branches = bool(base["finite"] and base["branch"])
    minimum_gram = base["minimum_gram"]
    minimum_argument = base["minimum_argument"]
    maximum_imaginary = base["maximum_imaginary"]
    for step in JACOBIAN_STEPS:
        matrix = np.zeros((34, 34))
        for coordinate in range(34):
            displacement = np.zeros(34)
            displacement[coordinate] = step
            plus = binary_equation_evaluation(
                model, t, z+displacement, boundary_log
            )
            minus = binary_equation_evaluation(
                model, t, z-displacement, boundary_log
            )
            for state in (plus, minus):
                all_branches &= bool(state["finite"] and state["branch"])
                if state["minimum_gram"] is not None:
                    minimum_gram = min(minimum_gram, state["minimum_gram"])
                    minimum_argument = min(
                        minimum_argument, state["minimum_argument"]
                    )
                    maximum_imaginary = max(
                        maximum_imaginary, state["maximum_imaginary"]
                    )
            if plus["transverse"] is None or minus["transverse"] is None:
                matrix[:, coordinate] = np.nan
            else:
                matrix[:, coordinate] = (
                    plus["transverse"]-minus["transverse"]
                )/(2.0*step)
        matrices.append(matrix)
    j1, j2, j3 = matrices
    r12 = (4.0*j2-j1)/3.0
    r23_raw = (4.0*j3-j2)/3.0
    j6_raw = (16.0*r23_raw-r12)/15.0
    j6 = (j6_raw+j6_raw.T)/2.0
    j23 = (r23_raw+r23_raw.T)/2.0
    delta_j = j6-j23
    singular6 = np.linalg.svd(j6, compute_uv=False)
    singular23 = np.linalg.svd(j23, compute_uv=False)
    eigenvalues = np.linalg.eigvalsh(j6)
    return {
        "base": base,
        "j6": j6,
        "j23": j23,
        "delta_j": delta_j,
        "singular_values": singular6,
        "singular_values23": singular23,
        "eigenvalues": eigenvalues,
        "all_204_branches_pass": all_branches,
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
    }


def precondition(jacobian, action_row):
    transverse = action_row["transverse"]
    transverse_error = action_row["transverse_error"]
    j6 = jacobian["j6"]
    j23 = jacobian["j23"]
    p6 = np.linalg.solve(j6, transverse)
    p23 = np.linalg.solve(j23, transverse)
    p_error = np.linalg.solve(j6, transverse_error)
    p_norm = float(np.linalg.norm(p6))
    p_error_norm = float(np.linalg.norm(p_error))
    step_change = float(
        np.linalg.norm(p6-p23)/max(p_norm, 1e-30)
    )
    model_error = float(
        np.linalg.norm(jacobian["delta_j"]@p6)
        / max(action_row["transverse_norm"], 1e-30)
    )
    usable = bool(
        jacobian["all_204_branches_pass"]
        and jacobian["singular_values"][-1] > 1e-12
        and jacobian["singular_values23"][-1] > 1e-12
        and step_change <= 0.1
        and model_error <= 0.1
    )
    zero_consistent = bool(
        usable
        and action_row["transverse_norm"]
        <= 10.0*max(action_row["transverse_error_norm"], 1e-30)
        and p_norm <= 10.0*max(p_error_norm, 1e-30)
        and p_error_norm < 1e-5
    )
    return {
        "p6": p6,
        "p23": p23,
        "p_error": p_error,
        "p_norm": p_norm,
        "p_error_norm": p_error_norm,
        "step_change": step_change,
        "model_error": model_error,
        "usable": usable,
        "zero_consistent": zero_consistent,
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
_ARB_VARIABLES = None


def initialize_action_worker(model, variable_strings):
    global _ARB_MODEL, _ARB_VARIABLES
    _ARB_MODEL = model
    configure_arb_precision()
    _ARB_VARIABLES = [
        [arb.mpf(value) for value in variables]
        for variables in variable_strings
    ]


def action_worker(point):
    state_index, coordinate, logarithmic_step = point
    variables = list(_ARB_VARIABLES[state_index])
    variables[coordinate] *= arb.exp(arb.mpf(str(logarithmic_step)))
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


def action_rows(model, states, steps, phase):
    if not states:
        return []
    points = []
    branch_records = []
    for state_index, state in enumerate(states):
        variables = complete_variables(
            state["t"], state["z"], state["boundary_log"]
        )
        state["variables_for_action"] = variables
        branch_pass = True
        minimum_gram = math.inf
        minimum_argument = math.inf
        for step in steps:
            for coordinate in range(35):
                for sign in (1.0, -1.0):
                    logarithmic_step = sign*step
                    points.append((
                        state_index, coordinate, logarithmic_step
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

    variable_strings = [
        [f"{value:.17g}" for value in state["variables_for_action"]]
        for state in states
    ]
    print(
        f"{phase}: {len(states)} rows, {len(points)} action points...",
        flush=True,
    )
    fork_context = mp.get_context("fork")
    with fork_context.Pool(
        processes=8,
        initializer=initialize_action_worker,
        initargs=(model, variable_strings),
    ) as pool:
        action_values = pool.map(action_worker, points, chunksize=1)

    rows = []
    cursor = 0
    for state_index, state in enumerate(states):
        derivatives = {}
        maximum_imaginary = arb.mpf(0)
        for step in steps:
            step_arb = arb.mpf(str(step))
            row = []
            for _coordinate in range(35):
                plus, minus = action_values[cursor:cursor+2]
                cursor += 2
                derivative = (plus-minus)/(48*step_arb)
                row.append(derivative)
                maximum_imaginary = max(
                    maximum_imaginary,
                    abs(arb.im(plus)), abs(arb.im(minus)),
                    abs(arb.im(derivative)),
                )
            derivatives[step] = row
        row1, row2, row3 = (derivatives[step] for step in steps)
        richardson12 = [(4*b-a)/3 for a, b in zip(row1, row2)]
        richardson23 = [(4*c-b)/3 for b, c in zip(row2, row3)]
        sixth_order = [
            (16*r23-r12)/15
            for r12, r23 in zip(richardson12, richardson23)
        ]
        error_arb = [
            value-r23 for value, r23 in zip(sixth_order, richardson23)
        ]
        equation = np.array([float(arb.re(value)) for value in sixth_order])
        error_row = np.array([float(arb.re(value)) for value in error_arb])
        transverse = quotient_basis.T@equation
        transverse_error = quotient_basis.T@error_row
        tangent = path_tangent(state["t"])
        scalar = float(tangent@equation)
        scalar_error = float(tangent@error_row)
        maximum_imaginary = float(max(
            maximum_imaginary,
            max(abs(arb.im(value)) for value in sixth_order),
            max(abs(arb.im(value)) for value in error_arb),
        ))
        rows.append({
            "phase": phase,
            "steps": steps,
            "sixth_order_arb": sixth_order,
            "error_arb": error_arb,
            "equation": equation,
            "error_row": error_row,
            "transverse": transverse,
            "transverse_error": transverse_error,
            "transverse_norm": float(np.linalg.norm(transverse)),
            "transverse_error_norm": float(
                np.linalg.norm(transverse_error)
            ),
            "scalar": scalar,
            "scalar_error": scalar_error,
            "maximum_imaginary": maximum_imaginary,
            "imaginary_below_1e-80": maximum_imaginary < 1e-80,
            "branch": branch_records[state_index],
        })
    if cursor != len(action_values):
        raise RuntimeError("action row cursor mismatch")
    return rows


def robust_descent(current_row, trial_row):
    return bool(
        trial_row["branch"]["pass"]
        and trial_row["imaginary_below_1e-80"]
        and trial_row["transverse_norm"]
        + 10.0*trial_row["transverse_error_norm"]
        < current_row["transverse_norm"]
        - 10.0*current_row["transverse_error_norm"]
    )


def jacobian_summary(record, response):
    return {
        "j6": record["j6"].tolist(),
        "j23": record["j23"].tolist(),
        "delta_j": record["delta_j"].tolist(),
        "singular_values": record["singular_values"].tolist(),
        "singular_values23": record["singular_values23"].tolist(),
        "eigenvalues": record["eigenvalues"].tolist(),
        "all_204_branches_pass": record["all_204_branches_pass"],
        "minimum_gram": record["minimum_gram"],
        "minimum_argument": record["minimum_argument"],
        "maximum_imaginary": record["maximum_imaginary"],
        "response": response_summary(response),
    }


def response_summary(response):
    return {
        "p6": response["p6"].tolist(),
        "p23": response["p23"].tolist(),
        "p_error": response["p_error"].tolist(),
        "p_norm": response["p_norm"],
        "p_error_norm": response["p_error_norm"],
        "step_change": response["step_change"],
        "model_error": response["model_error"],
        "usable": response["usable"],
        "zero_consistent": response["zero_consistent"],
    }


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 80),
        "imaginary": arb.nstr(arb.im(value), 80),
    }


def action_row_summary(row):
    return {
        "phase": row["phase"],
        "steps": row["steps"],
        "sixth_order": [
            serialize_complex(value) for value in row["sixth_order_arb"]
        ],
        "error": [
            serialize_complex(value) for value in row["error_arb"]
        ],
        "sixth_order_float": row["equation"].tolist(),
        "error_float": row["error_row"].tolist(),
        "transverse": row["transverse"].tolist(),
        "transverse_error": row["transverse_error"].tolist(),
        "transverse_norm": row["transverse_norm"],
        "transverse_error_norm": row["transverse_error_norm"],
        "collective_scalar": row["scalar"],
        "collective_scalar_error": row["scalar_error"],
        "maximum_imaginary": row["maximum_imaginary"],
        "imaginary_below_1e-80": row["imaginary_below_1e-80"],
        "branch": row["branch"],
    }


def solve_transverse_batch(model, states, batch_label):
    for state in states:
        state["accepted_iterations"] = 0
        state["solver_history"] = []
        state["solver_outcome"] = None
        state["current_row"] = None
        state["final_jacobian"] = None
        state["final_response"] = None

    initial_rows = action_rows(
        model, states, SOLVER_STEPS, f"{batch_label} solver initial"
    )
    for state, row in zip(states, initial_rows):
        state["current_row"] = row

    while any(state["solver_outcome"] is None for state in states):
        pending = [
            state for state in states if state["solver_outcome"] is None
        ]
        step_ready = []
        for index, state in enumerate(pending, 1):
            jacobian = local_jacobian(
                model, state["t"], state["z"], state["boundary_log"]
            )
            response = precondition(jacobian, state["current_row"])
            state["final_jacobian"] = jacobian
            state["final_response"] = response
            state["solver_history"].append({
                "kind": "ITERATE",
                "accepted_iterations": state["accepted_iterations"],
                "z": state["z"].copy(),
                "action_row": state["current_row"],
                "jacobian": jacobian,
                "response": response,
            })
            if response["zero_consistent"]:
                state["solver_outcome"] = (
                    "SOLVER_TRANSVERSE_ZERO_CONSISTENT"
                )
            elif not response["usable"]:
                state["solver_outcome"] = "LOCAL_PRECONDITIONER_UNRESOLVED"
            elif state["accepted_iterations"] >= MAX_ACCEPTED_ITERATIONS:
                state["solver_outcome"] = "SOLVER_ITERATION_LIMIT"
            else:
                step_ready.append(state)
            print(
                f"{batch_label}: local preconditioner {index}/{len(pending)}",
                flush=True,
            )

        remaining = list(step_ready)
        for damping in DAMPING:
            if not remaining:
                break
            trial_states = []
            for state in remaining:
                trial_states.append({
                    **state,
                    "z": state["z"]-damping*state["final_response"]["p6"],
                })
            trial_rows = action_rows(
                model, trial_states, SOLVER_STEPS,
                f"{batch_label} damping {damping:g}",
            )
            next_remaining = []
            for state, trial_state, trial_row in zip(
                remaining, trial_states, trial_rows
            ):
                accepted = robust_descent(state["current_row"], trial_row)
                provisional_zero = bool(
                    trial_row["branch"]["pass"]
                    and trial_row["imaginary_below_1e-80"]
                    and trial_row["transverse_norm"]
                    <= 10.0*max(
                        trial_row["transverse_error_norm"], 1e-30
                    )
                )
                zero_jacobian = None
                zero_response = None
                if not accepted and provisional_zero:
                    zero_jacobian = local_jacobian(
                        model, trial_state["t"], trial_state["z"],
                        trial_state["boundary_log"],
                    )
                    zero_response = precondition(zero_jacobian, trial_row)
                    accepted = zero_response["zero_consistent"]
                state["solver_history"].append({
                    "kind": "TRIAL",
                    "accepted_iterations": state["accepted_iterations"],
                    "damping": damping,
                    "z": trial_state["z"].copy(),
                    "action_row": trial_row,
                    "robust_descent": robust_descent(
                        state["current_row"], trial_row
                    ),
                    "provisional_zero": provisional_zero,
                    "zero_jacobian": zero_jacobian,
                    "zero_response": zero_response,
                    "accepted": accepted,
                })
                if accepted:
                    state["z"] = trial_state["z"]
                    state["current_row"] = trial_row
                    state["accepted_iterations"] += 1
                    if zero_response is not None:
                        state["final_jacobian"] = zero_jacobian
                        state["final_response"] = zero_response
                        state["solver_outcome"] = (
                            "SOLVER_TRANSVERSE_ZERO_CONSISTENT"
                        )
                else:
                    next_remaining.append(state)
            remaining = next_remaining
        for state in remaining:
            state["solver_outcome"] = "NO_ROBUST_ACTION_DESCENT"

    return states


def validate_transverse_batch(model, states, batch_label):
    candidates = [
        state for state in states
        if state["solver_outcome"] == "SOLVER_TRANSVERSE_ZERO_CONSISTENT"
    ]
    rows = action_rows(
        model, candidates, VALIDATION_STEPS,
        f"{batch_label} independent validation",
    )
    for state in states:
        state["validation"] = None
        state["scalar_label"] = None
    for state, row in zip(candidates, rows):
        response = precondition(state["final_jacobian"], row)
        binary = binary_equation_evaluation(
            model, state["t"], state["z"], state["boundary_log"]
        )
        binary_agreement = float(np.linalg.norm(
            binary["equation"]-row["equation"]
        ))
        gates = {
            "transverse_zero_consistent": (
                row["transverse_norm"]
                <= 10.0*max(row["transverse_error_norm"], 1e-30)
            ),
            "preconditioned_zero_consistent": (
                response["p_norm"]
                <= 10.0*max(response["p_error_norm"], 1e-30)
            ),
            "preconditioned_error_below_1e-5": (
                response["p_error_norm"] < 1e-5
            ),
            "all_210_branches_pass": row["branch"]["pass"],
            "imaginary_below_1e-80": row["imaginary_below_1e-80"],
            "binary_action_agreement_below_1e-8": binary_agreement < 1e-8,
        }
        transverse_validated = all(gates.values())
        scalar_floor = max(abs(row["scalar_error"]), 1e-30)
        if not transverse_validated:
            scalar_label = "SCALAR_NOT_CLASSIFIED"
        elif abs(row["scalar"]) <= 10.0*scalar_floor:
            scalar_label = "REDUCED_SCALAR_ZERO_CONSISTENT"
        elif abs(row["scalar"]) > 100.0*scalar_floor:
            scalar_label = "REDUCED_SCALAR_RESOLVED_NONZERO"
        else:
            scalar_label = "REDUCED_SCALAR_UNRESOLVED"
        state["validation"] = {
            "action_row": row,
            "response": response,
            "binary_equation": binary["equation"],
            "binary_action_agreement": binary_agreement,
            "gates": gates,
            "transverse_validated": transverse_validated,
            "scalar_label": scalar_label,
            "scalar_sign": (
                1 if row["scalar"] > 0 else -1 if row["scalar"] < 0 else 0
            ),
            "scalar_over_error": abs(row["scalar"])/scalar_floor,
        }
        state["scalar_label"] = scalar_label
    return states


def make_grid_states(parity):
    states = []
    for case in defect_input["parities"][parity]["cases"]:
        boundary_log = np.array(case["boundary_log_vector"], dtype=float)
        for t in T_GRID:
            stored = case["grid"][f"{t:+.2f}"]
            states.append({
                "parity": parity,
                "direction_index": int(case["direction_index"]),
                "sign": int(case["sign"]),
                "t": float(t),
                "boundary_log": boundary_log.copy(),
                "z": np.array(stored["z"], dtype=float),
                "origin": "GRID",
            })
    return states


def state_key(state):
    return (
        state["parity"], state["direction_index"],
        state["sign"], state["t"],
    )


def serialize_jacobian(record, response):
    if record is None or response is None:
        return None
    return jacobian_summary(record, response)


def serialize_history_item(item):
    result = {
        "kind": item["kind"],
        "accepted_iterations": item["accepted_iterations"],
        "z": item["z"].tolist(),
        "action_row": action_row_summary(item["action_row"]),
    }
    if item["kind"] == "ITERATE":
        result["jacobian"] = jacobian_summary(
            item["jacobian"], item["response"]
        )
    else:
        result.update({
            "damping": item["damping"],
            "robust_descent": item["robust_descent"],
            "provisional_zero": item["provisional_zero"],
            "accepted": item["accepted"],
            "zero_jacobian": serialize_jacobian(
                item["zero_jacobian"], item["zero_response"]
            ),
        })
    return result


def serialize_validation(record):
    if record is None:
        return None
    return {
        "action_row": action_row_summary(record["action_row"]),
        "response": response_summary(record["response"]),
        "binary_equation": record["binary_equation"].tolist(),
        "binary_action_agreement": record["binary_action_agreement"],
        "gates": record["gates"],
        "transverse_validated": record["transverse_validated"],
        "scalar_label": record["scalar_label"],
        "scalar_sign": record["scalar_sign"],
        "scalar_over_error": record["scalar_over_error"],
    }


def serialize_state(state):
    return {
        "parity": state["parity"],
        "direction_index": state["direction_index"],
        "sign": state["sign"],
        "t": state["t"],
        "origin": state["origin"],
        "boundary_log_vector": state["boundary_log"].tolist(),
        "final_z": state["z"].tolist(),
        "accepted_iterations": state["accepted_iterations"],
        "solver_outcome": state["solver_outcome"],
        "solver_history": [
            serialize_history_item(item) for item in state["solver_history"]
        ],
        "final_solver_action_row": action_row_summary(
            state["current_row"]
        ),
        "final_jacobian": serialize_jacobian(
            state["final_jacobian"], state["final_response"]
        ),
        "validation": serialize_validation(state["validation"]),
        "scalar_label": state["scalar_label"],
    }


parity_results = load_checkpoint()
all_grid_states = []
for parity, model in dust.bl.models.items():
    resumed = parity in parity_results
    if resumed:
        states = parity_results[parity]["grid_states"]
    else:
        states = make_grid_states(parity)
    check(
        f"{parity}: exactly forty preregistered grid states were loaded",
        len(states) == 40
        and len({state_key(state) for state in states}) == 40,
    )
    if not resumed:
        solve_transverse_batch(model, states, f"{parity} grid")
        validate_transverse_batch(model, states, f"{parity} grid")
    check(
        f"{parity}: every grid state has a frozen solver outcome",
        all(state["solver_outcome"] in {
            "SOLVER_TRANSVERSE_ZERO_CONSISTENT",
            "LOCAL_PRECONDITIONER_UNRESOLVED",
            "SOLVER_ITERATION_LIMIT",
            "NO_ROBUST_ACTION_DESCENT",
        } for state in states),
    )
    check(
        f"{parity}: every solver-zero state has independent validation",
        all(
            (state["validation"] is not None)
            == (state["solver_outcome"]
                == "SOLVER_TRANSVERSE_ZERO_CONSISTENT")
            for state in states
        ),
    )
    if not resumed:
        parity_results[parity] = {"grid_states": states, "bisections": []}
        save_checkpoint(parity_results)
    all_grid_states.extend(states)


# Frozen adjacent sign brackets.  The expected positive pattern may produce
# none; the general implementation remains deterministic if a sign differs.
all_bisection_states = []
for parity, model in dust.bl.models.items():
    grid_states = parity_results[parity]["grid_states"]
    for direction_index in range(1, 5):
        for sign in (-1, 1):
            ordered = sorted([
                state for state in grid_states
                if state["direction_index"] == direction_index
                and state["sign"] == sign
            ], key=lambda state: state["t"])
            for left, right in zip(ordered[:-1], ordered[1:]):
                left_v = left["validation"]
                right_v = right["validation"]
                if not (
                    left_v is not None and right_v is not None
                    and left_v["transverse_validated"]
                    and right_v["transverse_validated"]
                    and left_v["scalar_label"]
                    == "REDUCED_SCALAR_RESOLVED_NONZERO"
                    and right_v["scalar_label"]
                    == "REDUCED_SCALAR_RESOLVED_NONZERO"
                    and left_v["scalar_sign"]*right_v["scalar_sign"] < 0
                ):
                    continue
                record = {
                    "direction_index": direction_index,
                    "sign": sign,
                    "initial_interval": [left["t"], right["t"]],
                    "history": [],
                    "outcome": None,
                }
                current_left, current_right = left, right
                for _iteration in range(30):
                    midpoint = {
                        "parity": parity,
                        "direction_index": direction_index,
                        "sign": sign,
                        "t": (current_left["t"]+current_right["t"])/2.0,
                        "boundary_log": left["boundary_log"].copy(),
                        "z": (current_left["z"]+current_right["z"])/2.0,
                        "origin": "BISECTION",
                    }
                    solve_transverse_batch(
                        model, [midpoint],
                        f"{parity} bisection d{direction_index} s{sign}",
                    )
                    validate_transverse_batch(
                        model, [midpoint],
                        f"{parity} bisection validation d{direction_index} s{sign}",
                    )
                    record["history"].append(midpoint)
                    all_bisection_states.append(midpoint)
                    validation = midpoint["validation"]
                    if (
                        validation is None
                        or not validation["transverse_validated"]
                        or validation["scalar_label"]
                        == "REDUCED_SCALAR_UNRESOLVED"
                    ):
                        record["outcome"] = "BISECTION_NUMERICALLY_UNRESOLVED"
                        break
                    if validation["scalar_label"] == (
                        "REDUCED_SCALAR_ZERO_CONSISTENT"
                    ):
                        record["outcome"] = "BISECTION_STATIONARY_HIT"
                        break
                    if current_right["t"]-current_left["t"] < 1e-10:
                        record["outcome"] = (
                            "BRACKET_LOCALIZED_NOT_ZERO_VALIDATED"
                        )
                        break
                    if (
                        current_left["validation"]["scalar_sign"]
                        * validation["scalar_sign"] < 0
                    ):
                        current_right = midpoint
                    else:
                        current_left = midpoint
                if record["outcome"] is None:
                    record["outcome"] = (
                        "BRACKET_LOCALIZED_NOT_ZERO_VALIDATED"
                    )
                parity_results[parity]["bisections"].append(record)


case_records = []
for parity, record in parity_results.items():
    for direction_index in range(1, 5):
        for sign in (-1, 1):
            states = sorted([
                state for state in record["grid_states"]
                if state["direction_index"] == direction_index
                and state["sign"] == sign
            ], key=lambda state: state["t"])
            bisections = [
                item for item in record["bisections"]
                if item["direction_index"] == direction_index
                and item["sign"] == sign
            ]
            grid_hit = any(
                state["scalar_label"] == "REDUCED_SCALAR_ZERO_CONSISTENT"
                and state["validation"] is not None
                and state["validation"]["transverse_validated"]
                for state in states
            )
            bisection_hit = any(
                item["outcome"] == "BISECTION_STATIONARY_HIT"
                for item in bisections
            )
            all_grid_resolved = all(
                state["validation"] is not None
                and state["validation"]["transverse_validated"]
                and state["scalar_label"]
                == "REDUCED_SCALAR_RESOLVED_NONZERO"
                for state in states
            )
            all_brackets_resolved = all(
                item["outcome"] in {
                    "BISECTION_STATIONARY_HIT",
                    "BRACKET_LOCALIZED_NOT_ZERO_VALIDATED",
                } for item in bisections
            )
            if grid_hit or bisection_hit:
                outcome = "NONLINEAR_STATIONARY_CONTINUATION_FOUND"
            elif all_grid_resolved and all_brackets_resolved:
                outcome = "NO_STATIONARY_POINT_IN_FROZEN_GRID_AND_BRACKETS"
            else:
                outcome = "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED"
            case_records.append({
                "parity": parity,
                "direction_index": direction_index,
                "sign": sign,
                "grid_states": states,
                "bisections": bisections,
                "outcome": outcome,
            })


signed_hits = sum(
    case["outcome"] == "NONLINEAR_STATIONARY_CONTINUATION_FOUND"
    for case in case_records
)
pair_hits = 0
contrast_hits = 0
for direction_index in range(1, 5):
    contrast_hit = False
    for parity in ("even", "odd"):
        hit = any(
            case["outcome"] == "NONLINEAR_STATIONARY_CONTINUATION_FOUND"
            for case in case_records
            if case["direction_index"] == direction_index
            and case["parity"] == parity
        )
        pair_hits += hit
        contrast_hit |= hit
    contrast_hits += contrast_hit


all_grid_validated = all(
    state["validation"] is not None
    and state["validation"]["transverse_validated"]
    for state in all_grid_states
)
all_grid_nonzero = all(
    state["scalar_label"] == "REDUCED_SCALAR_RESOLVED_NONZERO"
    for state in all_grid_states
)
grid_signs = {
    state["validation"]["scalar_sign"]
    for state in all_grid_states
    if state["validation"] is not None
    and state["scalar_label"] == "REDUCED_SCALAR_RESOLVED_NONZERO"
}
all_cases_no_hit = all(
    case["outcome"] == "NO_STATIONARY_POINT_IN_FROZEN_GRID_AND_BRACKETS"
    for case in case_records
)
if signed_hits > 0:
    global_outcome = "NONLINEAR_STATIONARY_CONTINUATION_FOUND"
elif (
    all_grid_validated and all_grid_nonzero
    and len(grid_signs) == 1 and all_cases_no_hit
):
    global_outcome = "SIGN_DEFINITE_REDUCED_SCALAR_ON_FROZEN_GRID"
elif all_cases_no_hit:
    global_outcome = "NO_HIT_MIXED_REDUCED_SIGNS_ON_FROZEN_SCAN"
else:
    global_outcome = "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED"


check(
    "all sixteen signed cases receive exactly one frozen outcome",
    len(case_records) == 16
    and all(case["outcome"] in {
        "NONLINEAR_STATIONARY_CONTINUATION_FOUND",
        "NO_STATIONARY_POINT_IN_FROZEN_GRID_AND_BRACKETS",
        "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED",
    } for case in case_records),
)
check(
    "hit fractions cover 16 cases, 8 pairs and 4 contrasts",
    0 <= signed_hits <= 16
    and 0 <= pair_hits <= 8
    and 0 <= contrast_hits <= 4,
    f"signed={signed_hits}/16, pairs={pair_hits}/8, contrasts={contrast_hits}/4",
)
check(
    "the global label follows the preregistered hierarchy",
    global_outcome in {
        "NONLINEAR_STATIONARY_CONTINUATION_FOUND",
        "SIGN_DEFINITE_REDUCED_SCALAR_ON_FROZEN_GRID",
        "NO_HIT_MIXED_REDUCED_SIGNS_ON_FROZEN_SCAN",
        "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED",
    },
    f"{global_outcome}; grid signs={sorted(grid_signs)}",
)


payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "stagnation_result_commit": STAGNATION_RESULT_COMMIT,
    "solver_steps": SOLVER_STEPS,
    "validation_steps": VALIDATION_STEPS,
    "jacobian_steps": JACOBIAN_STEPS,
    "action_precision_digits": 100,
    "eta": defect_input["eta"],
    "grid": T_GRID,
    "cases": [
        {
            "parity": case["parity"],
            "direction_index": case["direction_index"],
            "sign": case["sign"],
            "outcome": case["outcome"],
            "grid_states": [
                serialize_state(state) for state in case["grid_states"]
            ],
            "bisections": [
                {
                    "direction_index": item["direction_index"],
                    "sign": item["sign"],
                    "initial_interval": item["initial_interval"],
                    "outcome": item["outcome"],
                    "history": [
                        serialize_state(state) for state in item["history"]
                    ],
                }
                for item in case["bisections"]
            ],
        }
        for case in case_records
    ],
    "attempt_count": {
        "grid_states": 80,
        "mechanically_forced_bisection_midpoints": len(all_bisection_states),
        "total": 80+len(all_bisection_states),
    },
    "hit_fractions": {
        "signed_cases": {"hits": int(signed_hits), "total": 16},
        "direction_parity_pairs": {"hits": int(pair_hits), "total": 8},
        "phase_contrasts": {"hits": int(contrast_hits), "total": 4},
    },
    "grid_summary": {
        "transverse_validated": sum(
            state["validation"] is not None
            and state["validation"]["transverse_validated"]
            for state in all_grid_states
        ),
        "scalar_zero_consistent": sum(
            state["scalar_label"] == "REDUCED_SCALAR_ZERO_CONSISTENT"
            for state in all_grid_states
        ),
        "scalar_resolved_nonzero": sum(
            state["scalar_label"] == "REDUCED_SCALAR_RESOLVED_NONZERO"
            for state in all_grid_states
        ),
        "scalar_unresolved": sum(
            state["scalar_label"] == "REDUCED_SCALAR_UNRESOLVED"
            for state in all_grid_states
        ),
        "resolved_nonzero_signs": sorted(grid_signs),
    },
    "global_outcome": global_outcome,
    "claim_boundary": {
        "finite_grid_result": "DERIVED COMPUTATIONAL IF CHECKS PASS",
        "method_and_positive_sign_provenance": "PATTERN-INFORMED",
        "continuous_interval_no_root": "NOT ESTABLISHED",
        "amplitude_scaling_and_quadratic_law": "NOT TESTED",
        "second_slab_and_full_carrier": "NOT TESTED",
        "clock_speed_limit_planck_scale": "NOT DERIVED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
CHECKPOINT.unlink(missing_ok=True)

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(
    f"{global_outcome}; grid {payload['grid_summary']}; "
    f"attempts {payload['attempt_count']}."
)
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Preregistered arbitrary-precision complete-action solver repair.

Prior-art commit: 6b7f9e4
Protocol commit: 4b6b10c

Only the complete Regge-plus-dust action supplies the physical residual.
The binary analytic Jacobian is used solely as a proposal generator and
temporary norm preconditioner.
"""

import contextlib
import hashlib
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
OLD_RESULT = HERE / "gravity_600cell_dust_action_reduced_scalar.json"
OUTPUT = HERE / "gravity_600cell_dust_action_solver_repair.json"
CHECKPOINT = HERE / ".gravity_600cell_dust_action_solver_repair.checkpoint.pkl"

PRIOR_ART_COMMIT = "6b7f9e4"
PROTOCOL_COMMIT = "4b6b10c"
FROZEN_RESULT_COMMIT = "64a13f6"
OLD_PROTOCOL_COMMIT = "17f9560"

OPERATIONAL_STEPS = ("1e-20", "1e-15")
VALIDATION_STEPS = ("3e-20", "3e-15")
BRANCH_AUDIT_STEP = 1.0e-6
JACOBIAN_STEPS = (5.0e-4, 2.5e-4, 1.25e-4)
DAMPING = tuple(2.0**(-power) for power in range(11))
T_GRID = (-0.10, -0.05, 0.0, 0.05, 0.10)
MAX_ACCEPTED_ITERATIONS = 12
PRECISION_DIGITS = 100
ARITHMETIC_FLOOR = 1.0e-60

tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


def checkpoint_metadata():
    return {
        "schema": 1,
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "frozen_result_commit": FROZEN_RESULT_COMMIT,
        "old_protocol_commit": OLD_PROTOCOL_COMMIT,
        "operational_steps": OPERATIONAL_STEPS,
        "validation_steps": VALIDATION_STEPS,
        "branch_audit_step": BRANCH_AUDIT_STEP,
        "jacobian_steps": JACOBIAN_STEPS,
        "damping": DAMPING,
        "t_grid": T_GRID,
        "max_accepted_iterations": MAX_ACCEPTED_ITERATIONS,
        "precision_digits": PRECISION_DIGITS,
        "arithmetic_floor": ARITHMETIC_FLOOR,
    }


def load_checkpoint():
    if not CHECKPOINT.exists():
        return {"completed_cases": {}, "bisections": []}
    with CHECKPOINT.open("rb") as stream:
        payload = pickle.load(stream)
    if payload.get("metadata") != checkpoint_metadata():
        raise RuntimeError(f"incompatible operational checkpoint: {CHECKPOINT}")
    completed = payload.get("completed_cases")
    bisections = payload.get("bisections")
    if not isinstance(completed, dict) or not isinstance(bisections, list):
        raise RuntimeError(f"malformed operational checkpoint: {CHECKPOINT}")
    print(
        f"Loaded checkpoint: {len(completed)}/16 cases, "
        f"{len(bisections)} bracket records.",
        flush=True,
    )
    return {"completed_cases": completed, "bisections": bisections}


def save_checkpoint(progress):
    payload = {"metadata": checkpoint_metadata(), **progress}
    temporary = CHECKPOINT.with_suffix(CHECKPOINT.suffix + ".tmp")
    with temporary.open("wb") as stream:
        pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)
        stream.flush()
    temporary.replace(CHECKPOINT)
    print(
        f"Saved checkpoint: {len(progress['completed_cases'])}/16 cases, "
        f"{len(progress['bisections'])} bracket records.",
        flush=True,
    )


precision_input = json.loads(PRECISION_INPUT.read_text())
old_result = json.loads(OLD_RESULT.read_text())
check(
    "the frozen carrier and numerical-boundary result retain provenance",
    precision_input["protocol_commit"] == "da34272"
    and precision_input["tests"] == precision_input["passed"] == 14
    and old_result["protocol_commit"] == OLD_PROTOCOL_COMMIT
    and old_result["tests"] == old_result["passed"] == 12
    and old_result["global_outcome"]
    == "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED",
)


print("Loading the independently certified complete-action control...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_solver_repair",
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


def configure_arb_precision():
    arb.mp.dps = PRECISION_DIGITS
    dust.ARB_I = arb.mpc(0, 1)
    dust.ARB_TAU = arb.mpf("0.0102")
    dust.ARB_M_STAR = arb.mpf(10)
    dust.ARB_ZETA = (
        arb.pi**2 * arb.sqrt(2) / 50
    ) ** (arb.mpf(1) / 3)
    dust.ARB_R0 = 4 * dust.ARB_M_STAR / (3 * arb.pi)
    dust.ARB_L0 = dust.ARB_ZETA * dust.ARB_R0
    dust.ARB_L0_SQUARE = dust.ARB_L0**2
    dust.ARB_EPSILON_3 = 2 * arb.pi - 5 * arb.acos(arb.mpf(1) / 3)
    dust.ARB_MASS = (
        (90 / arb.pi) * dust.ARB_EPSILON_3 * dust.ARB_L0
    )
    dust.ARB_TAU_SQUARE = dust.ARB_TAU**2
    dust.ARB_SLANT_SQUARE = dust.ARB_L0_SQUARE - dust.ARB_TAU_SQUARE
    # Mutate the list bound as arb_edge_square's default argument.
    dust.ARB_OLD_VALUES[:] = [dust.ARB_L0_SQUARE] * 30


configure_arb_precision()

quotient_basis_float = np.array(precision_input["quotient_basis"], dtype=float)
quotient_basis_arb = [
    [arb.mpf(str(value)) for value in row]
    for row in precision_input["quotient_basis"]
]
check(
    "the exact frozen complement retains its 35 by 34 carrier",
    quotient_basis_float.shape == (35, 34)
    and np.linalg.norm(
        quotient_basis_float.T @ quotient_basis_float - np.eye(34), 2
    ) < 1e-12,
)
check(
    "the derivative pairs are disjoint and have frozen comparable scales",
    not set(OPERATIONAL_STEPS) & set(VALIDATION_STEPS)
    and OPERATIONAL_STEPS == ("1e-20", "1e-15")
    and VALIDATION_STEPS == ("3e-20", "3e-15"),
)


def mp_string(value, digits=80):
    return arb.nstr(value, digits)


def mp_norm(values):
    return arb.sqrt(sum(arb.mpf(value) ** 2 for value in values))


def path_tangent_arb(t):
    t_arb = arb.mpf(str(t))
    rho = dust.ARB_TAU_SQUARE * arb.exp(t_arb)
    diagonal = dust.ARB_L0_SQUARE - rho
    raw = [-rho / diagonal] * 30 + [arb.mpf(1)] * 5
    norm = mp_norm(raw)
    return [value / norm for value in raw]


def complete_variables_arb(t, z, boundary_log):
    t_arb = arb.mpf(str(t))
    rho = dust.ARB_TAU_SQUARE * arb.exp(t_arb)
    diagonal = dust.ARB_L0_SQUARE - rho
    base = [diagonal] * 30 + [rho] * 5
    qz = [
        sum(quotient_basis_arb[row][column] * z[column]
            for column in range(34))
        for row in range(35)
    ]
    internal = [base[row] * arb.exp(qz[row]) for row in range(35)]
    final = [
        dust.ARB_L0_SQUARE * arb.exp(value) for value in boundary_log
    ]
    return internal + final


def base_variables_arb():
    return (
        [dust.ARB_SLANT_SQUARE] * 30
        + [dust.ARB_TAU_SQUARE] * 5
        + [dust.ARB_L0_SQUARE] * 30
    )


def state_variables_arb(state):
    if state.get("control_variables") is not None:
        return list(state["control_variables"])
    return complete_variables_arb(
        state["t"], state["z"], state["boundary_log"]
    )


def branch_audit(model, variables):
    variables_float = np.array([float(value) for value in variables])
    branch_pass = True
    minimum_gram = math.inf
    minimum_argument = math.inf
    maximum_imaginary = 0.0
    audit_vectors = [variables_float]
    for coordinate in range(35):
        for sign in (1.0, -1.0):
            trial = variables_float.copy()
            trial[coordinate] *= math.exp(sign * BRANCH_AUDIT_STEP)
            audit_vectors.append(trial)
    for trial in audit_vectors:
        try:
            _, _, _, data = dust.total_reduced_evaluation(
                model, trial, dust.old_values
            )
            minimum_gram = min(minimum_gram, float(data["minimum_gram"]))
            minimum_argument = min(
                minimum_argument, float(data["minimum_argument"])
            )
            branch_pass &= bool(
                dict(data["negative_counts"]) == {1: 100}
                and float(data["minimum_gram"]) > 1e-8
                and float(data["minimum_argument"]) > 1e-6
            )
        except Exception:
            branch_pass = False
    return {
        "pass": bool(branch_pass),
        "audit_count": len(audit_vectors),
        "step": BRANCH_AUDIT_STEP,
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
    }


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
    variables[coordinate] *= arb.exp(arb.mpf(logarithmic_step))
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


def project_row(equation, error, tangent):
    equation_real = [arb.re(value) for value in equation]
    error_real = [arb.re(value) for value in error]
    transverse = [
        sum(quotient_basis_arb[row][column] * equation_real[row]
            for row in range(35))
        for column in range(34)
    ]
    transverse_error = [
        sum(quotient_basis_arb[row][column] * error_real[row]
            for row in range(35))
        for column in range(34)
    ]
    scalar = sum(tangent[row] * equation_real[row] for row in range(35))
    scalar_error = sum(tangent[row] * error_real[row] for row in range(35))
    return transverse, transverse_error, scalar, scalar_error


def action_rows(model, states, steps, phase):
    if not states:
        return []
    variables_by_state = [state_variables_arb(state) for state in states]
    branches = [branch_audit(model, variables) for variables in variables_by_state]
    points = []
    for state_index in range(len(states)):
        for step in steps:
            for coordinate in range(35):
                points.append((state_index, coordinate, step))
                points.append((state_index, coordinate, "-" + step))
    variable_strings = [
        [mp_string(value, 100) for value in variables]
        for variables in variables_by_state
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
    primary_step = arb.mpf(steps[0])
    shadow_step = arb.mpf(steps[1])
    for state_index, state in enumerate(states):
        derivatives = {}
        maximum_imaginary = arb.mpf(0)
        for step_string, step_arb in (
            (steps[0], primary_step), (steps[1], shadow_step)
        ):
            row = []
            for _coordinate in range(35):
                plus, minus = action_values[cursor:cursor + 2]
                cursor += 2
                derivative = (plus - minus) / (48 * step_arb)
                row.append(derivative)
                maximum_imaginary = max(
                    maximum_imaginary,
                    abs(arb.im(plus)), abs(arb.im(minus)),
                    abs(arb.im(derivative)),
                )
            derivatives[step_string] = row
        primary = derivatives[steps[0]]
        shadow = derivatives[steps[1]]
        error = [left - right for left, right in zip(primary, shadow)]
        tangent = (
            state.get("control_tangent")
            if state.get("control_tangent") is not None
            else path_tangent_arb(state["t"])
        )
        transverse, transverse_error, scalar, scalar_error = project_row(
            primary, error, tangent
        )
        maximum_imaginary = max(
            maximum_imaginary,
            max(abs(arb.im(value)) for value in error),
        )
        row = {
            "phase": phase,
            "steps": steps,
            "primary_arb": primary,
            "shadow_arb": shadow,
            "error_arb": error,
            "transverse_arb": transverse,
            "transverse_error_arb": transverse_error,
            "transverse": np.array([float(value) for value in transverse]),
            "transverse_error": np.array(
                [float(value) for value in transverse_error]
            ),
            "transverse_norm": float(mp_norm(transverse)),
            "transverse_error_norm": float(mp_norm(transverse_error)),
            "scalar_arb": scalar,
            "scalar_error_arb": scalar_error,
            "scalar": float(scalar),
            "scalar_error": float(scalar_error),
            "maximum_imaginary": float(maximum_imaginary),
            "imaginary_below_1e-70": maximum_imaginary < arb.mpf("1e-70"),
            "branch": branches[state_index],
        }
        rows.append(row)
    if cursor != len(action_values):
        raise RuntimeError("action row cursor mismatch")
    return rows


def complete_variables_float(t, z, boundary_log):
    z_float = np.array([float(value) for value in z], dtype=float)
    boundary_float = np.array(
        [float(value) for value in boundary_log], dtype=float
    )
    rho = dust.TAU_SQUARE * math.exp(float(t))
    diagonal = dust.L0_SQUARE - rho
    base = np.concatenate((np.full(30, diagonal), np.full(5, rho)))
    internal = base * np.exp(quotient_basis_float @ z_float)
    final = dust.L0_SQUARE * np.exp(boundary_float)
    return np.concatenate((internal, final))


def binary_equation_evaluation(model, t, z, boundary_log):
    variables = complete_variables_float(t, z, boundary_log)
    try:
        _, gradient, _, data = dust.total_reduced_evaluation(
            model, variables, dust.old_values
        )
        equation_complex = variables[:35] * gradient[:35] / 24.0
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
            "equation": equation,
            "transverse": quotient_basis_float.T @ equation,
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
    z_float = np.array([float(value) for value in z])
    for step in JACOBIAN_STEPS:
        matrix = np.zeros((34, 34))
        for coordinate in range(34):
            displacement = np.zeros(34)
            displacement[coordinate] = step
            plus = binary_equation_evaluation(
                model, t, z_float + displacement, boundary_log
            )
            minus = binary_equation_evaluation(
                model, t, z_float - displacement, boundary_log
            )
            for record in (plus, minus):
                all_branches &= bool(record["finite"] and record["branch"])
                if record["minimum_gram"] is not None:
                    minimum_gram = min(minimum_gram, record["minimum_gram"])
                    minimum_argument = min(
                        minimum_argument, record["minimum_argument"]
                    )
                    maximum_imaginary = max(
                        maximum_imaginary, record["maximum_imaginary"]
                    )
            if plus["transverse"] is None or minus["transverse"] is None:
                matrix[:, coordinate] = np.nan
            else:
                matrix[:, coordinate] = (
                    plus["transverse"] - minus["transverse"]
                ) / (2.0 * step)
        matrices.append(matrix)
    j1, j2, j3 = matrices
    r12 = (4.0 * j2 - j1) / 3.0
    r23_raw = (4.0 * j3 - j2) / 3.0
    j6_raw = (16.0 * r23_raw - r12) / 15.0
    j6 = (j6_raw + j6_raw.T) / 2.0
    j23 = (r23_raw + r23_raw.T) / 2.0
    delta_j = j6 - j23
    singular6 = np.linalg.svd(j6, compute_uv=False)
    singular23 = np.linalg.svd(j23, compute_uv=False)
    eigenvalues = np.linalg.eigvalsh(j6)
    return {
        "j6": j6,
        "j23": j23,
        "delta_j": delta_j,
        "singular_values": singular6,
        "singular_values23": singular23,
        "eigenvalues": eigenvalues,
        "all_204_branches_pass": bool(all_branches),
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
    }


def precondition(jacobian, row):
    transverse = row["transverse"]
    transverse_error = row["transverse_error"]
    j6 = jacobian["j6"]
    j23 = jacobian["j23"]
    p = np.linalg.solve(j6, transverse)
    p23 = np.linalg.solve(j23, transverse)
    p_error = np.linalg.solve(j6, transverse_error)
    p_norm = float(np.linalg.norm(p))
    p_error_norm = float(np.linalg.norm(p_error))
    step_change = float(
        np.linalg.norm(p - p23) / max(p_norm, ARITHMETIC_FLOOR)
    )
    model_error = float(
        np.linalg.norm(jacobian["delta_j"] @ p)
        / max(row["transverse_norm"], ARITHMETIC_FLOOR)
    )
    usable = bool(
        jacobian["all_204_branches_pass"]
        and jacobian["singular_values"][-1] > 1e-12
        and jacobian["singular_values23"][-1] > 1e-12
        and step_change <= 0.1
        and model_error <= 0.1
    )
    row_valid = bool(
        row["branch"]["pass"] and row["imaginary_below_1e-70"]
    )
    propagated_accurate = bool(p_error_norm < 1e-5)
    zero_consistent = bool(
        row_valid and usable and propagated_accurate
        and row["transverse_norm"]
        <= 10.0 * max(row["transverse_error_norm"], ARITHMETIC_FLOOR)
        and p_norm <= 10.0 * max(p_error_norm, ARITHMETIC_FLOOR)
    )
    if p_norm > 10.0 * max(p_error_norm, ARITHMETIC_FLOOR):
        merit_kind = "NATURAL_PRECONDITIONED"
        merit = p_norm
        merit_error = p_error_norm
    else:
        merit_kind = "RAW_TRANSVERSE"
        merit = row["transverse_norm"]
        merit_error = row["transverse_error_norm"]
    return {
        "p": p,
        "p23": p23,
        "p_error": p_error,
        "p_norm": p_norm,
        "p_error_norm": p_error_norm,
        "step_change": step_change,
        "model_error": model_error,
        "usable": usable,
        "row_valid": row_valid,
        "propagated_accurate": propagated_accurate,
        "zero_consistent": zero_consistent,
        "merit_kind": merit_kind,
        "merit": merit,
        "merit_error": merit_error,
    }


def trial_merit(jacobian, row, merit_kind):
    if merit_kind == "NATURAL_PRECONDITIONED":
        value = np.linalg.solve(jacobian["j6"], row["transverse"])
        error = np.linalg.solve(
            jacobian["j6"], row["transverse_error"]
        )
        return float(np.linalg.norm(value)), float(np.linalg.norm(error))
    if merit_kind == "RAW_TRANSVERSE":
        return row["transverse_norm"], row["transverse_error_norm"]
    raise ValueError(f"unknown merit kind {merit_kind}")


def calibration_reference(parity):
    record = dust.records[parity]
    variables = base_variables_arb()
    return [
        variables[index] * record["arb_total_gradient"][index] / 24
        for index in range(35)
    ]


def row_agreement(left, right):
    component_difference = max(
        abs(arb.re(a - b))
        for a, b in zip(left["primary_arb"], right["primary_arb"])
    )
    component_proxy = max(
        abs(arb.re(a)) + abs(arb.re(b))
        for a, b in zip(left["error_arb"], right["error_arb"])
    )
    transverse_difference = float(np.linalg.norm(
        left["transverse"] - right["transverse"]
    ))
    transverse_proxy = (
        left["transverse_error_norm"] + right["transverse_error_norm"]
    )
    return {
        "component_difference": float(component_difference),
        "component_proxy": float(component_proxy),
        "component_pass": bool(
            component_difference
            <= 10 * max(component_proxy, arb.mpf("1e-60"))
        ),
        "transverse_difference": transverse_difference,
        "transverse_proxy": transverse_proxy,
        "transverse_pass": bool(
            transverse_difference
            <= 10 * max(transverse_proxy, ARITHMETIC_FLOOR)
        ),
    }


def calibrate_derivatives():
    records = {}
    for parity, model in dust.bl.models.items():
        state = {
            "control_variables": base_variables_arb(),
            "control_tangent": [arb.mpf(1)] + [arb.mpf(0)] * 34,
        }
        operational = action_rows(
            model, [state], OPERATIONAL_STEPS,
            f"{parity} derivative calibration operational",
        )[0]
        validation = action_rows(
            model, [state], VALIDATION_STEPS,
            f"{parity} derivative calibration validation",
        )[0]
        reference = calibration_reference(parity)
        operational_reference_error = max(
            abs(arb.re(left - right))
            for left, right in zip(operational["primary_arb"], reference)
        )
        validation_reference_error = max(
            abs(arb.re(left - right))
            for left, right in zip(validation["primary_arb"], reference)
        )
        agreement = row_agreement(operational, validation)
        gates = {
            "operational_branch": operational["branch"]["pass"],
            "validation_branch": validation["branch"]["pass"],
            "operational_reference_below_1e-10": (
                operational_reference_error < arb.mpf("1e-10")
            ),
            "validation_reference_below_1e-10": (
                validation_reference_error < arb.mpf("1e-10")
            ),
            "primary_rows_agree": agreement["component_pass"],
            "operational_imaginary_below_1e-70": (
                operational["imaginary_below_1e-70"]
            ),
            "validation_imaginary_below_1e-70": (
                validation["imaginary_below_1e-70"]
            ),
        }
        records[parity] = {
            "operational": operational,
            "validation": validation,
            "operational_reference_error": float(
                operational_reference_error
            ),
            "validation_reference_error": float(validation_reference_error),
            "agreement": agreement,
            "gates": gates,
            "pass": all(gates.values()),
        }
    return records


def decimal_list(values):
    return [arb.mpf(str(value)) for value in values]


def state_key(state):
    return (
        state["parity"], state["direction_index"],
        state["sign"], float(state["t"]),
    )


def case_key(parity, direction_index, sign):
    return f"{parity}:d{int(direction_index)}:s{int(sign):+d}"


def make_frozen_cases():
    cases = []
    for old_case in old_result["cases"]:
        states = []
        for old_state in old_case["grid_states"]:
            states.append({
                "parity": old_state["parity"],
                "direction_index": int(old_state["direction_index"]),
                "sign": int(old_state["sign"]),
                "t": float(old_state["t"]),
                "origin": "GRID_REPAIR_FROM_64a13f6",
                "boundary_log": decimal_list(
                    old_state["boundary_log_vector"]
                ),
                "z": decimal_list(old_state["final_z"]),
            })
        cases.append({
            "parity": old_case["parity"],
            "direction_index": int(old_case["direction_index"]),
            "sign": int(old_case["sign"]),
            "states": sorted(states, key=lambda state: state["t"]),
        })
    cases.sort(key=lambda item: (
        0 if item["parity"] == "even" else 1,
        item["direction_index"], item["sign"],
    ))
    return cases


def solve_transverse_batch(model, states, batch_label):
    for state in states:
        state["accepted_iterations"] = 0
        state["solver_history"] = []
        state["solver_outcome"] = None
        state["current_row"] = None
        state["final_jacobian"] = None
        state["final_response"] = None

    initial_rows = action_rows(
        model, states, OPERATIONAL_STEPS, f"{batch_label} initial"
    )
    for state, row in zip(states, initial_rows):
        state["current_row"] = row

    while any(state["solver_outcome"] is None for state in states):
        pending = [
            state for state in states if state["solver_outcome"] is None
        ]
        step_ready = []
        for index, state in enumerate(pending, 1):
            try:
                jacobian = local_jacobian(
                    model, state["t"], state["z"], state["boundary_log"]
                )
                response = precondition(jacobian, state["current_row"])
            except Exception as error:
                state["solver_history"].append({
                    "kind": "PRECONDITIONER_EXCEPTION",
                    "accepted_iterations": state["accepted_iterations"],
                    "z": list(state["z"]),
                    "error": repr(error),
                })
                state["solver_outcome"] = "LOCAL_PRECONDITIONER_UNRESOLVED"
                continue
            state["final_jacobian"] = jacobian
            state["final_response"] = response
            state["solver_history"].append({
                "kind": "ITERATE",
                "accepted_iterations": state["accepted_iterations"],
                "z": list(state["z"]),
                "action_row": state["current_row"],
                "jacobian": jacobian,
                "response": response,
            })
            if not response["row_valid"]:
                state["solver_outcome"] = "ACTION_ROW_ACCURACY_UNRESOLVED"
            elif not response["usable"]:
                state["solver_outcome"] = "LOCAL_PRECONDITIONER_UNRESOLVED"
            elif not response["propagated_accurate"]:
                state["solver_outcome"] = "PROPAGATED_ACCURACY_UNRESOLVED"
            elif response["zero_consistent"]:
                state["solver_outcome"] = (
                    "SOLVER_TRANSVERSE_ZERO_CONSISTENT"
                )
            elif state["accepted_iterations"] >= MAX_ACCEPTED_ITERATIONS:
                state["solver_outcome"] = "SOLVER_ITERATION_LIMIT"
            else:
                step_ready.append(state)
            print(
                f"{batch_label}: preconditioner {index}/{len(pending)}",
                flush=True,
            )

        remaining = list(step_ready)
        for damping in DAMPING:
            if not remaining:
                break
            damping_arb = arb.mpf(str(damping))
            trial_states = []
            for state in remaining:
                trial_states.append({
                    **state,
                    "z": [
                        value - damping_arb * arb.mpf(str(step))
                        for value, step in zip(
                            state["z"], state["final_response"]["p"]
                        )
                    ],
                })
            trial_rows = action_rows(
                model, trial_states, OPERATIONAL_STEPS,
                f"{batch_label} damping {damping:g}",
            )
            next_remaining = []
            for state, trial_state, trial_row in zip(
                remaining, trial_states, trial_rows
            ):
                trial_value, trial_error = trial_merit(
                    state["final_jacobian"], trial_row,
                    state["final_response"]["merit_kind"],
                )
                current_value = state["final_response"]["merit"]
                current_error = state["final_response"]["merit_error"]
                robust_descent = bool(
                    trial_row["branch"]["pass"]
                    and trial_row["imaginary_below_1e-70"]
                    and trial_value + 10.0 * trial_error
                    < current_value - 10.0 * current_error
                )
                state["solver_history"].append({
                    "kind": "TRIAL",
                    "accepted_iterations": state["accepted_iterations"],
                    "damping": damping,
                    "z": list(trial_state["z"]),
                    "action_row": trial_row,
                    "merit_kind": state["final_response"]["merit_kind"],
                    "current_merit": current_value,
                    "current_merit_error": current_error,
                    "trial_merit": trial_value,
                    "trial_merit_error": trial_error,
                    "accepted": robust_descent,
                })
                if robust_descent:
                    state["z"] = trial_state["z"]
                    state["current_row"] = trial_row
                    state["accepted_iterations"] += 1
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
    validation_rows = action_rows(
        model, candidates, VALIDATION_STEPS,
        f"{batch_label} independent validation",
    )
    for state in states:
        state["validation"] = None
        state["scalar_label"] = None
    for state, row in zip(candidates, validation_rows):
        response = precondition(state["final_jacobian"], row)
        agreement = row_agreement(state["current_row"], row)
        binary = binary_equation_evaluation(
            model, state["t"], state["z"], state["boundary_log"]
        )
        binary_difference = (
            float(np.linalg.norm(
                binary["equation"] - np.array([
                    float(arb.re(value)) for value in row["primary_arb"]
                ])
            ))
            if binary["equation"] is not None else math.inf
        )
        gates = {
            "transverse_zero_consistent": (
                row["transverse_norm"]
                <= 10.0 * max(
                    row["transverse_error_norm"], ARITHMETIC_FLOOR
                )
            ),
            "preconditioned_zero_consistent": (
                response["p_norm"]
                <= 10.0 * max(
                    response["p_error_norm"], ARITHMETIC_FLOOR
                )
            ),
            "preconditioned_error_below_1e-5": (
                response["p_error_norm"] < 1e-5
            ),
            "validation_row_valid": response["row_valid"],
            "preconditioner_usable": response["usable"],
            "full_rows_agree": agreement["component_pass"],
            "transverse_rows_agree": agreement["transverse_pass"],
        }
        transverse_validated = all(gates.values())
        scalar_floor = max(abs(row["scalar_error"]), ARITHMETIC_FLOOR)
        if not transverse_validated:
            scalar_label = "SCALAR_NOT_CLASSIFIED"
        elif abs(row["scalar"]) <= 10.0 * scalar_floor:
            scalar_label = "REDUCED_SCALAR_ZERO_CONSISTENT"
        elif abs(row["scalar"]) > 100.0 * scalar_floor:
            scalar_label = "REDUCED_SCALAR_RESOLVED_NONZERO"
        else:
            scalar_label = "REDUCED_SCALAR_UNRESOLVED"
        state["validation"] = {
            "action_row": row,
            "response": response,
            "agreement": agreement,
            "binary_equation": binary["equation"],
            "binary_action_difference": binary_difference,
            "gates": gates,
            "transverse_validated": transverse_validated,
            "scalar_label": scalar_label,
            "scalar_sign": (
                1 if row["scalar"] > 0 else -1 if row["scalar"] < 0 else 0
            ),
            "scalar_over_error": abs(row["scalar"]) / scalar_floor,
        }
        state["scalar_label"] = scalar_label
    return states


def serialize_complex(value):
    return {
        "real": mp_string(arb.re(value)),
        "imaginary": mp_string(arb.im(value)),
    }


def action_row_summary(row, compact=False):
    result = {
        "phase": row["phase"],
        "steps": list(row["steps"]),
        "transverse": [mp_string(value) for value in row["transverse_arb"]],
        "transverse_error": [
            mp_string(value) for value in row["transverse_error_arb"]
        ],
        "transverse_norm": row["transverse_norm"],
        "transverse_error_norm": row["transverse_error_norm"],
        "collective_scalar": mp_string(row["scalar_arb"]),
        "collective_scalar_error": mp_string(row["scalar_error_arb"]),
        "maximum_imaginary": row["maximum_imaginary"],
        "imaginary_below_1e-70": row["imaginary_below_1e-70"],
        "branch": row["branch"],
    }
    if not compact:
        result.update({
            "primary": [
                serialize_complex(value) for value in row["primary_arb"]
            ],
            "shadow": [
                serialize_complex(value) for value in row["shadow_arb"]
            ],
            "stability_proxy": [
                serialize_complex(value) for value in row["error_arb"]
            ],
        })
    return result


def response_summary(response):
    return {
        "p": response["p"].tolist(),
        "p23": response["p23"].tolist(),
        "p_error": response["p_error"].tolist(),
        "p_norm": response["p_norm"],
        "p_error_norm": response["p_error_norm"],
        "step_change": response["step_change"],
        "model_error": response["model_error"],
        "usable": response["usable"],
        "row_valid": response["row_valid"],
        "propagated_accurate": response["propagated_accurate"],
        "zero_consistent": response["zero_consistent"],
        "merit_kind": response["merit_kind"],
        "merit": response["merit"],
        "merit_error": response["merit_error"],
    }


def matrix_hash(matrix):
    return hashlib.sha256(np.asarray(matrix, dtype="<f8").tobytes()).hexdigest()


def jacobian_summary(jacobian, response):
    if jacobian is None or response is None:
        return None
    return {
        "j6_sha256": matrix_hash(jacobian["j6"]),
        "j23_sha256": matrix_hash(jacobian["j23"]),
        "delta_j_sha256": matrix_hash(jacobian["delta_j"]),
        "singular_values": jacobian["singular_values"].tolist(),
        "singular_values23": jacobian["singular_values23"].tolist(),
        "eigenvalues": jacobian["eigenvalues"].tolist(),
        "all_204_branches_pass": jacobian["all_204_branches_pass"],
        "minimum_gram": jacobian["minimum_gram"],
        "minimum_argument": jacobian["minimum_argument"],
        "maximum_imaginary": jacobian["maximum_imaginary"],
        "response": response_summary(response),
    }


def history_summary(item):
    result = {
        "kind": item["kind"],
        "accepted_iterations": item["accepted_iterations"],
        "z": [mp_string(value) for value in item["z"]],
    }
    if item["kind"] == "PRECONDITIONER_EXCEPTION":
        result["error"] = item["error"]
    elif item["kind"] == "ITERATE":
        result.update({
            "action_row": action_row_summary(item["action_row"]),
            "jacobian": jacobian_summary(
                item["jacobian"], item["response"]
            ),
        })
    else:
        result.update({
            "damping": item["damping"],
            "action_row": action_row_summary(
                item["action_row"], compact=True
            ),
            "merit_kind": item["merit_kind"],
            "current_merit": item["current_merit"],
            "current_merit_error": item["current_merit_error"],
            "trial_merit": item["trial_merit"],
            "trial_merit_error": item["trial_merit_error"],
            "accepted": item["accepted"],
        })
    return result


def validation_summary(validation):
    if validation is None:
        return None
    return {
        "action_row": action_row_summary(validation["action_row"]),
        "response": response_summary(validation["response"]),
        "operational_validation_agreement": validation["agreement"],
        "binary_equation": (
            validation["binary_equation"].tolist()
            if validation["binary_equation"] is not None else None
        ),
        "binary_action_difference": validation["binary_action_difference"],
        "gates": validation["gates"],
        "transverse_validated": validation["transverse_validated"],
        "scalar_label": validation["scalar_label"],
        "scalar_sign": validation["scalar_sign"],
        "scalar_over_error": validation["scalar_over_error"],
    }


def state_summary(state):
    return {
        "parity": state["parity"],
        "direction_index": state["direction_index"],
        "sign": state["sign"],
        "t": state["t"],
        "origin": state["origin"],
        "boundary_log_vector": [
            mp_string(value) for value in state["boundary_log"]
        ],
        "final_z": [mp_string(value) for value in state["z"]],
        "accepted_iterations": state["accepted_iterations"],
        "solver_outcome": state["solver_outcome"],
        "solver_history": [
            history_summary(item) for item in state["solver_history"]
        ],
        "final_solver_action_row": action_row_summary(state["current_row"]),
        "final_jacobian": jacobian_summary(
            state["final_jacobian"], state["final_response"]
        ),
        "validation": validation_summary(state["validation"]),
        "scalar_label": state["scalar_label"],
    }


def calibration_summary(records):
    return {
        parity: {
            "operational": action_row_summary(record["operational"]),
            "validation": action_row_summary(record["validation"]),
            "operational_reference_error": record[
                "operational_reference_error"
            ],
            "validation_reference_error": record[
                "validation_reference_error"
            ],
            "agreement": record["agreement"],
            "gates": record["gates"],
            "pass": record["pass"],
        }
        for parity, record in records.items()
    }


print("Running the preregistered derivative calibration...", flush=True)
calibration = calibrate_derivatives()
calibration_pass = bool(
    set(calibration) == {"even", "odd"}
    and all(record["pass"] for record in calibration.values())
)
check(
    "the calibration gate was evaluated for both parities before target rows",
    set(calibration) == {"even", "odd"}
    and all(set(record["gates"]) == {
        "operational_branch", "validation_branch",
        "operational_reference_below_1e-10",
        "validation_reference_below_1e-10",
        "primary_rows_agree",
        "operational_imaginary_below_1e-70",
        "validation_imaginary_below_1e-70",
    } for record in calibration.values()),
    f"calibration_pass={calibration_pass}",
)

if not calibration_pass:
    payload = {
        "protocol_commit": PROTOCOL_COMMIT,
        "prior_art_commit": PRIOR_ART_COMMIT,
        "frozen_result_commit": FROZEN_RESULT_COMMIT,
        "operational_steps": OPERATIONAL_STEPS,
        "validation_steps": VALIDATION_STEPS,
        "action_precision_digits": PRECISION_DIGITS,
        "calibration": calibration_summary(calibration),
        "target_rows_evaluated": 0,
        "global_outcome": "DERIVATIVE_CALIBRATION_FAILED",
        "claim_boundary": "OPEN NUMERICALLY; no target state was evaluated",
        "tests": tests,
        "passed": passed,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    CHECKPOINT.unlink(missing_ok=True)
    print("-" * 78)
    print(f"RESULT: {passed}/{tests} implementation checks passed")
    print("DERIVATIVE_CALIBRATION_FAILED; no target row evaluated.")
    raise SystemExit(0 if passed == tests else 1)


frozen_cases = make_frozen_cases()
all_frozen_states = [state for case in frozen_cases for state in case["states"]]
check(
    "exactly sixteen frozen cases and eighty unique grid states were loaded",
    len(frozen_cases) == 16
    and len(all_frozen_states) == 80
    and len({state_key(state) for state in all_frozen_states}) == 80
    and {state["t"] for state in all_frozen_states} == set(T_GRID),
)

progress = load_checkpoint()
for case in frozen_cases:
    key = case_key(case["parity"], case["direction_index"], case["sign"])
    if key in progress["completed_cases"]:
        case["states"] = progress["completed_cases"][key]
        continue
    model = dust.bl.models[case["parity"]]
    label = (
        f"{case['parity']} d{case['direction_index']} "
        f"s{case['sign']:+d}"
    )
    solve_transverse_batch(model, case["states"], label)
    validate_transverse_batch(model, case["states"], label)
    progress["completed_cases"][key] = case["states"]
    save_checkpoint(progress)

check(
    "all sixteen signed cases completed the frozen grid repair",
    len(progress["completed_cases"]) == 16
    and all(len(states) == 5 for states in progress["completed_cases"].values()),
)


def eligible_bracket(left, right):
    return bool(
        left["validation"] is not None
        and right["validation"] is not None
        and left["validation"]["transverse_validated"]
        and right["validation"]["transverse_validated"]
        and left["scalar_label"] == "REDUCED_SCALAR_RESOLVED_NONZERO"
        and right["scalar_label"] == "REDUCED_SCALAR_RESOLVED_NONZERO"
        and left["validation"]["scalar_sign"]
        * right["validation"]["scalar_sign"] < 0
    )


def bracket_key(parity, direction_index, sign, left_t, right_t):
    return (
        parity, int(direction_index), int(sign),
        float(left_t), float(right_t),
    )


def find_bracket_record(key):
    for record in progress["bisections"]:
        if tuple(record["key"]) == key:
            return record
    return None


all_bisection_states = []
for case in frozen_cases:
    model = dust.bl.models[case["parity"]]
    states = sorted(case["states"], key=lambda state: state["t"])
    for initial_left, initial_right in zip(states[:-1], states[1:]):
        if not eligible_bracket(initial_left, initial_right):
            continue
        key = bracket_key(
            case["parity"], case["direction_index"], case["sign"],
            initial_left["t"], initial_right["t"],
        )
        record = find_bracket_record(key)
        if record is None:
            record = {
                "key": key,
                "parity": case["parity"],
                "direction_index": case["direction_index"],
                "sign": case["sign"],
                "initial_interval": [initial_left["t"], initial_right["t"]],
                "left_state": initial_left,
                "right_state": initial_right,
                "history": [],
                "outcome": None,
            }
            progress["bisections"].append(record)
            save_checkpoint(progress)
        while record["outcome"] is None and len(record["history"]) < 30:
            left = record["left_state"]
            right = record["right_state"]
            midpoint = {
                "parity": case["parity"],
                "direction_index": case["direction_index"],
                "sign": case["sign"],
                "t": (left["t"] + right["t"]) / 2.0,
                "origin": "BISECTION_REPAIR",
                "boundary_log": list(initial_left["boundary_log"]),
                "z": [
                    (a + b) / 2 for a, b in zip(left["z"], right["z"])
                ],
            }
            label = (
                f"{case['parity']} bisection d{case['direction_index']} "
                f"s{case['sign']:+d} n{len(record['history']) + 1}"
            )
            solve_transverse_batch(model, [midpoint], label)
            validate_transverse_batch(model, [midpoint], label)
            record["history"].append(midpoint)
            all_bisection_states.append(midpoint)
            validation = midpoint["validation"]
            if (
                validation is None
                or not validation["transverse_validated"]
                or validation["scalar_label"] in {
                    "SCALAR_NOT_CLASSIFIED", "REDUCED_SCALAR_UNRESOLVED",
                }
            ):
                record["outcome"] = "BISECTION_NUMERICALLY_UNRESOLVED"
            elif validation["scalar_label"] == (
                "REDUCED_SCALAR_ZERO_CONSISTENT"
            ):
                record["outcome"] = "BISECTION_STATIONARY_HIT"
            elif right["t"] - left["t"] < 1e-10:
                record["outcome"] = "BRACKET_LOCALIZED_NOT_ZERO_VALIDATED"
            elif (
                left["validation"]["scalar_sign"]
                * validation["scalar_sign"] < 0
            ):
                record["right_state"] = midpoint
            else:
                record["left_state"] = midpoint
            save_checkpoint(progress)
        if record["outcome"] is None:
            record["outcome"] = "BRACKET_LOCALIZED_NOT_ZERO_VALIDATED"
            save_checkpoint(progress)


case_records = []
for case in frozen_cases:
    relevant_bisections = [
        record for record in progress["bisections"]
        if record["parity"] == case["parity"]
        and record["direction_index"] == case["direction_index"]
        and record["sign"] == case["sign"]
    ]
    states = sorted(case["states"], key=lambda state: state["t"])
    grid_hit = any(
        state["validation"] is not None
        and state["validation"]["transverse_validated"]
        and state["scalar_label"] == "REDUCED_SCALAR_ZERO_CONSISTENT"
        for state in states
    )
    bisection_hit = any(
        record["outcome"] == "BISECTION_STATIONARY_HIT"
        for record in relevant_bisections
    )
    all_grid_resolved = all(
        state["validation"] is not None
        and state["validation"]["transverse_validated"]
        and state["scalar_label"] == "REDUCED_SCALAR_RESOLVED_NONZERO"
        for state in states
    )
    all_brackets_resolved = all(
        record["outcome"] in {
            "BISECTION_STATIONARY_HIT",
            "BRACKET_LOCALIZED_NOT_ZERO_VALIDATED",
        }
        for record in relevant_bisections
    )
    if grid_hit or bisection_hit:
        outcome = "NONLINEAR_STATIONARY_CONTINUATION_FOUND"
    elif all_grid_resolved and all_brackets_resolved:
        outcome = "NO_STATIONARY_POINT_IN_FROZEN_GRID_AND_BRACKETS"
    else:
        outcome = "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED"
    case_records.append({
        **case,
        "bisections": relevant_bisections,
        "outcome": outcome,
    })


all_grid_states = [state for case in case_records for state in case["states"]]
allowed_solver_outcomes = {
    "SOLVER_TRANSVERSE_ZERO_CONSISTENT",
    "ACTION_ROW_ACCURACY_UNRESOLVED",
    "LOCAL_PRECONDITIONER_UNRESOLVED",
    "PROPAGATED_ACCURACY_UNRESOLVED",
    "SOLVER_ITERATION_LIMIT",
    "NO_ROBUST_ACTION_DESCENT",
}
check(
    "every grid state has one frozen solver outcome and validation iff eligible",
    all(state["solver_outcome"] in allowed_solver_outcomes
        for state in all_grid_states)
    and all(
        (state["validation"] is not None)
        == (state["solver_outcome"]
            == "SOLVER_TRANSVERSE_ZERO_CONSISTENT")
        for state in all_grid_states
    ),
)

allowed_case_outcomes = {
    "NONLINEAR_STATIONARY_CONTINUATION_FOUND",
    "NO_STATIONARY_POINT_IN_FROZEN_GRID_AND_BRACKETS",
    "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED",
}
check(
    "all sixteen signed cases receive exactly one preregistered outcome",
    len(case_records) == 16
    and all(case["outcome"] in allowed_case_outcomes for case in case_records),
)

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
        pair_hits += int(hit)
        contrast_hit |= hit
    contrast_hits += int(contrast_hit)
check(
    "hit fractions retain the frozen 16/8/4 look-elsewhere accounting",
    0 <= signed_hits <= 16
    and 0 <= pair_hits <= 8
    and 0 <= contrast_hits <= 4,
    f"signed={signed_hits}/16, pairs={pair_hits}/8, contrasts={contrast_hits}/4",
)

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
    "the global label follows the preregistered outcome hierarchy",
    global_outcome in {
        "NONLINEAR_STATIONARY_CONTINUATION_FOUND",
        "SIGN_DEFINITE_REDUCED_SCALAR_ON_FROZEN_GRID",
        "NO_HIT_MIXED_REDUCED_SIGNS_ON_FROZEN_SCAN",
        "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED",
    },
    f"{global_outcome}; resolved signs={sorted(grid_signs)}",
)


def bisection_summary(record):
    return {
        "parity": record["parity"],
        "direction_index": record["direction_index"],
        "sign": record["sign"],
        "initial_interval": record["initial_interval"],
        "outcome": record["outcome"],
        "history": [state_summary(state) for state in record["history"]],
    }


payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "frozen_result_commit": FROZEN_RESULT_COMMIT,
    "old_protocol_commit": OLD_PROTOCOL_COMMIT,
    "operational_steps": OPERATIONAL_STEPS,
    "validation_steps": VALIDATION_STEPS,
    "branch_audit_step": BRANCH_AUDIT_STEP,
    "jacobian_steps": JACOBIAN_STEPS,
    "action_precision_digits": PRECISION_DIGITS,
    "calibration": calibration_summary(calibration),
    "target_rows_evaluated": "ALL_ROWS_RECORDED_IN_CASE_HISTORIES",
    "cases": [
        {
            "parity": case["parity"],
            "direction_index": case["direction_index"],
            "sign": case["sign"],
            "outcome": case["outcome"],
            "grid_states": [state_summary(state) for state in case["states"]],
            "bisections": [
                bisection_summary(record) for record in case["bisections"]
            ],
        }
        for case in case_records
    ],
    "attempt_count": {
        "grid_states": 80,
        "mechanically_forced_bisection_midpoints": sum(
            len(record["history"]) for record in progress["bisections"]
        ),
        "total": 80 + sum(
            len(record["history"]) for record in progress["bisections"]
        ),
    },
    "hit_fractions": {
        "signed_cases": {"hits": int(signed_hits), "total": 16},
        "direction_parity_pairs": {"hits": int(pair_hits), "total": 8},
        "phase_contrasts": {"hits": int(contrast_hits), "total": 4},
    },
    "grid_summary": {
        "solver_outcomes": {
            outcome: sum(
                state["solver_outcome"] == outcome for state in all_grid_states
            )
            for outcome in sorted(allowed_solver_outcomes)
        },
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
        "finite_grid_result": "DERIVED COMPUTATIONAL IF VALIDATED",
        "positive_sign_provenance": "PATTERN-INFORMED",
        "continuous_interval_no_root": "NOT ESTABLISHED",
        "second_slab_and_full_carrier": "NOT TESTED",
        "physical_time_mass_speed_planck_units": "NOT DERIVED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
CHECKPOINT.unlink(missing_ok=True)

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(
    f"{global_outcome}; grid={payload['grid_summary']}; "
    f"attempts={payload['attempt_count']}",
)
raise SystemExit(0 if passed == tests else 1)


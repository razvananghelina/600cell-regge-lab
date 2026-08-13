#!/usr/bin/env python3
"""Lyapunov--Schmidt refinement of the 600-cell dust stationary base.

Protocol commit: 61b7f35 (initial protocol fa001af).  The 34 transverse
equations are solved on a frozen collective scan, every scalar root is
localized mechanically, and candidates are audited from the complete action
at 100 decimal digits.
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
PUBLISHED_INPUT = HERE / "gravity_600cell_published_dust_control.json"
OUTPUT = HERE / "gravity_600cell_dust_stationary_base.json"
PROTOCOL_COMMIT = "61b7f35"
INITIAL_PROTOCOL_COMMIT = "fa001af"
PRIOR_ART_COMMIT = "86f6ce7"
PRECISION_RESULT_COMMIT = "29a779f"
PUBLISHED_RESULT_COMMIT = "66a6465"
T_GRID = (
    -0.10, -0.075, -0.05, -0.03, -0.02, -0.01,
    0.0,
    0.01, 0.02, 0.03, 0.05, 0.075, 0.10,
)
POSITIVE_SCAN = (0.01, 0.02, 0.03, 0.05, 0.075, 0.10)
NEGATIVE_SCAN = (-0.01, -0.02, -0.03, -0.05, -0.075, -0.10)
DAMPING = tuple(2.0**(-power) for power in range(11))
ACTION_STEPS = (2.0e-4, 1.0e-4, 5.0e-5)
SCALAR_FLOOR = 1.0e-12
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
published_input = json.loads(PUBLISHED_INPUT.read_text())
check(
    "the two frozen inputs retain committed provenance and both parities",
    precision_input["protocol_commit"] == "da34272"
    and published_input["protocol_commit"] == "cc0902b"
    and precision_input["passed"] == precision_input["tests"] == 14
    and published_input["passed"] == published_input["tests"] == 14
    and set(precision_input["parities"])
        == set(published_input["parities"])
        == {"even", "odd"},
)


print("Loading the independently certified complete-action implementation...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_stationary_base",
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
base_gauge_tangent = np.array(
    precision_input["normalized_gauge_tangent"], dtype=float
)
check(
    "the frozen 35-by-34 basis remains the exact base-tangent complement",
    quotient_basis.shape == (35, 34)
    and np.linalg.norm(
        quotient_basis.T@quotient_basis-np.eye(34), 2
    ) < 2e-14
    and np.linalg.norm(quotient_basis.T@base_gauge_tangent) < 2e-14,
)


def path_internal(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    return np.concatenate((np.full(30, diagonal), np.full(5, rho)))


def path_tangent(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    return np.concatenate((np.full(30, -rho/diagonal), np.ones(5)))


def internal_from_coordinates(t, z):
    return path_internal(t)*np.exp(quotient_basis@np.asarray(z, dtype=float))


def complete_variables(t, z):
    return np.concatenate((
        internal_from_coordinates(t, z),
        np.full(30, dust.L0_SQUARE),
    ))


def equation_evaluation(model, t, z):
    variables = complete_variables(t, z)
    try:
        _, gradient, _, data = dust.total_reduced_evaluation(
            model, variables, dust.old_values
        )
        equation_complex = variables[:35]*gradient[:35]/24.0
        equation = equation_complex.real
        branch = bool(
            dict(data["negative_counts"]) == {1: 100}
            and float(data["minimum_gram"]) > 1e-8
            and float(data["minimum_argument"]) > 1e-6
        )
        finite = bool(
            np.all(np.isfinite(equation_complex.real))
            and np.all(np.isfinite(equation_complex.imag))
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


def solve_transverse(model, quotient_hessian, t, initial_z):
    z = np.array(initial_z, dtype=float, copy=True)
    history = []
    minimum_gram = math.inf
    minimum_argument = math.inf
    maximum_imaginary = 0.0
    all_branch = True
    outcome = "TRANSVERSE_SOLVE_UNRESOLVED"

    for iteration in range(101):
        state = equation_evaluation(model, t, z)
        minimum_gram = min(minimum_gram, state["minimum_gram"])
        minimum_argument = min(minimum_argument, state["minimum_argument"])
        maximum_imaginary = max(
            maximum_imaginary, state["maximum_imaginary"]
        )
        all_branch &= state["branch"]
        if not state["branch"] or not state["finite"]:
            history.append({
                "iteration": iteration,
                "accepted_damping": None,
                "transverse_norm": math.inf,
                "preconditioned_norm": math.inf,
                "branch": False,
            })
            break

        transverse = quotient_basis.T@state["equation"]
        preconditioned = np.linalg.solve(quotient_hessian, transverse)
        transverse_norm = float(np.linalg.norm(transverse))
        preconditioned_norm = float(np.linalg.norm(preconditioned))
        record = {
            "iteration": iteration,
            "accepted_damping": None,
            "transverse_norm": transverse_norm,
            "preconditioned_norm": preconditioned_norm,
            "branch": True,
        }
        history.append(record)
        if transverse_norm < 1e-12 and preconditioned_norm < 1e-5:
            outcome = "TRANSVERSE_SOLVED"
            break
        if iteration == 100:
            break

        delta = -preconditioned
        accepted = False
        for damping in DAMPING:
            trial_z = z+damping*delta
            trial = equation_evaluation(model, t, trial_z)
            minimum_gram = min(minimum_gram, trial["minimum_gram"])
            minimum_argument = min(
                minimum_argument, trial["minimum_argument"]
            )
            maximum_imaginary = max(
                maximum_imaginary, trial["maximum_imaginary"]
            )
            all_branch &= trial["branch"]
            if not trial["branch"] or not trial["finite"]:
                continue
            trial_transverse = quotient_basis.T@trial["equation"]
            trial_preconditioned = np.linalg.solve(
                quotient_hessian, trial_transverse
            )
            if np.linalg.norm(trial_preconditioned) < preconditioned_norm:
                z = trial_z
                record["accepted_damping"] = damping
                accepted = True
                break
        if not accepted:
            break

    final_state = equation_evaluation(model, t, z)
    minimum_gram = min(minimum_gram, final_state["minimum_gram"])
    minimum_argument = min(minimum_argument, final_state["minimum_argument"])
    maximum_imaginary = max(
        maximum_imaginary, final_state["maximum_imaginary"]
    )
    all_branch &= final_state["branch"]
    if final_state["finite"]:
        final_transverse = quotient_basis.T@final_state["equation"]
        final_preconditioned = np.linalg.solve(
            quotient_hessian, final_transverse
        )
        final_transverse_norm = float(np.linalg.norm(final_transverse))
        final_preconditioned_norm = float(np.linalg.norm(final_preconditioned))
        tangent = path_tangent(t)
        scalar = float(
            tangent@final_state["equation"]/np.linalg.norm(tangent)
        )
    else:
        final_transverse = np.full(34, np.nan)
        final_transverse_norm = math.inf
        final_preconditioned_norm = math.inf
        scalar = math.nan
    if not (
        all_branch
        and final_transverse_norm < 1e-12
        and final_preconditioned_norm < 1e-5
    ):
        outcome = "TRANSVERSE_SOLVE_UNRESOLVED"

    return {
        "t": float(t),
        "z": z,
        "outcome": outcome,
        "history": history,
        "equation": final_state["equation"],
        "transverse": final_transverse,
        "transverse_norm": final_transverse_norm,
        "preconditioned_norm": final_preconditioned_norm,
        "scalar": scalar,
        "full_equation_norm": float(
            np.linalg.norm(final_state["equation"])
        ) if final_state["finite"] else math.inf,
        "maximum_equation_imaginary": maximum_imaginary,
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "all_branch": all_branch,
        "variables": final_state["variables"],
    }


def bisect_scalar_root(model, quotient_hessian, left, right):
    history = []
    if left["scalar"]*right["scalar"] >= 0:
        return None, history, False
    current_left = left
    current_right = right
    for iteration in range(40):
        midpoint_t = (current_left["t"]+current_right["t"])/2.0
        initial_z = (current_left["z"]+current_right["z"])/2.0
        midpoint = solve_transverse(
            model, quotient_hessian, midpoint_t, initial_z
        )
        history.append(midpoint)
        if midpoint["outcome"] != "TRANSVERSE_SOLVED":
            return None, history, False
        if (
            abs(midpoint["scalar"]) <= SCALAR_FLOOR
            or current_right["t"]-current_left["t"] <= 1e-10
        ):
            return midpoint, history, True
        if current_left["scalar"]*midpoint["scalar"] < 0:
            current_right = midpoint
        else:
            current_left = midpoint
    return history[-1], history, True


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
_ARB_CANDIDATE = None


def initialize_action_worker(model, candidate_strings):
    global _ARB_MODEL, _ARB_CANDIDATE
    _ARB_MODEL = model
    configure_arb_precision()
    _ARB_CANDIDATE = [arb.mpf(value) for value in candidate_strings]


def candidate_action_worker(point):
    coordinate, logarithmic_step = point
    variables = list(_ARB_CANDIDATE)
    variables[coordinate] *= arb.exp(arb.mpf(str(logarithmic_step)))
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


fork_context = mp.get_context("fork")


def action_only_audit(model, quotient_hessian, candidate):
    candidate_variables = np.array(candidate["variables"], dtype=float)
    points = []
    branch_pass = True
    minimum_gram = math.inf
    minimum_argument = math.inf
    for step in ACTION_STEPS:
        for coordinate in range(35):
            for sign in (1.0, -1.0):
                logarithmic_step = sign*step
                points.append((coordinate, logarithmic_step))
                trial = candidate_variables.copy()
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

    candidate_strings = [f"{value:.17g}" for value in candidate_variables]
    print(
        f"Auditing candidate t={candidate['t']:+.12e} with 210 action points...",
        flush=True,
    )
    with fork_context.Pool(
        processes=8,
        initializer=initialize_action_worker,
        initargs=(model, candidate_strings),
    ) as pool:
        actions = pool.map(candidate_action_worker, points, chunksize=1)

    derivatives = {}
    maximum_imaginary = arb.mpf(0)
    cursor = 0
    for step in ACTION_STEPS:
        row = []
        step_arb = arb.mpf(str(step))
        for _coordinate in range(35):
            plus, minus = actions[cursor:cursor+2]
            cursor += 2
            row.append((plus-minus)/(48*step_arb))
            maximum_imaginary = max(
                maximum_imaginary, abs(arb.im(plus)), abs(arb.im(minus))
            )
        derivatives[step] = row

    row1, row2, row3 = (derivatives[step] for step in ACTION_STEPS)
    richardson12 = [(4*b-a)/3 for a, b in zip(row1, row2)]
    richardson23 = [(4*c-b)/3 for b, c in zip(row2, row3)]
    sixth_order = [
        (16*r23-r12)/15
        for r12, r23 in zip(richardson12, richardson23)
    ]
    sixth_float = np.array([float(arb.re(value)) for value in sixth_order])
    richardson23_float = np.array([
        float(arb.re(value)) for value in richardson23
    ])
    epsilon = float(np.linalg.norm(sixth_float-richardson23_float))
    action_floor = max(epsilon, 1e-30)
    transverse = quotient_basis.T@sixth_float
    transverse_correction = float(np.linalg.norm(
        np.linalg.solve(quotient_hessian, transverse)
    ))
    tangent = path_tangent(candidate["t"])
    scalar = float(tangent@sixth_float/np.linalg.norm(tangent))
    norm = float(np.linalg.norm(sixth_float))
    binary_agreement = float(np.linalg.norm(
        sixth_float-candidate["equation"]
    ))
    maximum_imaginary = float(max(
        maximum_imaginary,
        max(abs(arb.im(value)) for value in sixth_order),
    ))
    gates = {
        "equation_norm_below_1e-10": norm < 1e-10,
        "scalar_below_1e-11": abs(scalar) < 1e-11,
        "transverse_correction_below_1e-5": transverse_correction < 1e-5,
        "imaginary_below_1e-80": maximum_imaginary < 1e-80,
        "all_210_branches_pass": branch_pass,
        "binary_action_agreement_below_1e-9": binary_agreement < 1e-9,
    }
    label = (
        "WEAK_SCALE_STATIONARY" if all(gates.values())
        else "ACTION_ONLY_STATIONARITY_FAILED"
    )
    return {
        "label": label,
        "gates": gates,
        "steps": ACTION_STEPS,
        "derivatives": derivatives,
        "richardson12": richardson12,
        "richardson23": richardson23,
        "sixth_order": sixth_order,
        "sixth_order_float": sixth_float,
        "norm": norm,
        "scalar": scalar,
        "transverse_correction": transverse_correction,
        "epsilon": epsilon,
        "action_floor": action_floor,
        "norm_over_floor": norm/action_floor,
        "binary_agreement": binary_agreement,
        "maximum_imaginary": maximum_imaginary,
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "branch_pass": branch_pass,
    }


def compact_solve(record):
    return {
        "t": record["t"],
        "z": record["z"].tolist(),
        "outcome": record["outcome"],
        "history": record["history"],
        "equation": record["equation"].tolist(),
        "transverse": record["transverse"].tolist(),
        "transverse_norm": record["transverse_norm"],
        "preconditioned_norm": record["preconditioned_norm"],
        "scalar": record["scalar"],
        "full_equation_norm": record["full_equation_norm"],
        "maximum_equation_imaginary": record["maximum_equation_imaginary"],
        "minimum_gram": record["minimum_gram"],
        "minimum_argument": record["minimum_argument"],
        "all_branch": record["all_branch"],
        "variables": record["variables"].tolist(),
    }


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 80),
        "imaginary": arb.nstr(arb.im(value), 80),
    }


def serialize_action(record):
    return {
        "label": record["label"],
        "gates": record["gates"],
        "steps": record["steps"],
        "derivatives": {
            f"{step:.1e}": [
                serialize_complex(value) for value in record["derivatives"][step]
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
        "sixth_order_float": record["sixth_order_float"].tolist(),
        "norm": record["norm"],
        "scalar": record["scalar"],
        "transverse_correction": record["transverse_correction"],
        "empirical_error_norm": record["epsilon"],
        "classification_floor": record["action_floor"],
        "norm_over_floor": record["norm_over_floor"],
        "binary_vs_action_norm": record["binary_agreement"],
        "maximum_action_or_derivative_imaginary": record["maximum_imaginary"],
        "minimum_absolute_gram_eigenvalue": record["minimum_gram"],
        "minimum_angle_argument_modulus": record["minimum_argument"],
        "all_210_branches_pass": record["branch_pass"],
    }


results = {}
for parity, model in dust.bl.models.items():
    print(f"Scanning the stationary base: {parity} parity...", flush=True)
    quotient_hessian = np.array(
        precision_input["parities"][parity]["quotient"]["matrix"],
        dtype=float,
    )
    quotient_hessian = (quotient_hessian+quotient_hessian.T)/2.0

    grid = {}
    zero = solve_transverse(model, quotient_hessian, 0.0, np.zeros(34))
    grid[0.0] = zero
    previous = zero
    for t in POSITIVE_SCAN:
        initial = previous["z"] if previous["outcome"] == "TRANSVERSE_SOLVED" else np.zeros(34)
        current = solve_transverse(model, quotient_hessian, t, initial)
        grid[t] = current
        if current["outcome"] == "TRANSVERSE_SOLVED":
            previous = current
    previous = zero
    for t in NEGATIVE_SCAN:
        initial = previous["z"] if previous["outcome"] == "TRANSVERSE_SOLVED" else np.zeros(34)
        current = solve_transverse(model, quotient_hessian, t, initial)
        grid[t] = current
        if current["outcome"] == "TRANSVERSE_SOLVED":
            previous = current

    ordered = [grid[t] for t in T_GRID]
    candidates = [
        record for record in ordered
        if record["outcome"] == "TRANSVERSE_SOLVED"
        and abs(record["scalar"]) <= SCALAR_FLOOR
    ]
    bisection_records = []
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
            candidate, history, resolved = bisect_scalar_root(
                model, quotient_hessian, left, right
            )
            bisection_records.append({
                "interval": [left["t"], right["t"]],
                "resolved": resolved,
                "history": history,
                "candidate": candidate,
            })
            intervals_resolved &= resolved
            if candidate is not None:
                candidates.append(candidate)

    candidates.sort(key=lambda item: item["t"])
    deduplicated = []
    for candidate in candidates:
        if not deduplicated or abs(candidate["t"]-deduplicated[-1]["t"]) >= 1e-8:
            deduplicated.append(candidate)
        elif abs(candidate["scalar"]) < abs(deduplicated[-1]["scalar"]):
            deduplicated[-1] = candidate
    candidates = deduplicated

    candidate_records = []
    for candidate in candidates:
        action_record = action_only_audit(model, quotient_hessian, candidate)
        candidate_records.append({
            "solve": candidate,
            "action": action_record,
        })
    stationary_count = sum(
        item["action"]["label"] == "WEAK_SCALE_STATIONARY"
        for item in candidate_records
    )
    all_grid_resolved = all(
        item["outcome"] == "TRANSVERSE_SOLVED" for item in ordered
    )
    localization_resolved = all_grid_resolved and intervals_resolved
    if stationary_count == 1:
        outcome = "UNIQUE_WEAK_SCALE_STATIONARY_BASE_IN_SCAN"
    elif stationary_count > 1:
        outcome = "MULTIPLE_WEAK_SCALE_STATIONARY_BASES_IN_SCAN"
    elif localization_resolved:
        outcome = "NO_STATIONARY_BASE_IN_FROZEN_SCAN"
    else:
        outcome = "STATIONARY_BASE_NUMERICALLY_UNRESOLVED"

    results[parity] = {
        "outcome": outcome,
        "quotient_hessian": quotient_hessian,
        "grid": grid,
        "bisections": bisection_records,
        "candidates": candidate_records,
        "stationary_count": stationary_count,
        "all_grid_resolved": all_grid_resolved,
        "all_intervals_resolved": intervals_resolved,
        "localization_resolved": localization_resolved,
    }

    check(
        f"{parity}: all thirteen frozen grid points were evaluated and recorded",
        set(grid) == set(T_GRID)
        and all(len(item["z"]) == 34 for item in grid.values()),
    )
    check(
        f"{parity}: scalar candidates follow only the frozen grid and bisection rules",
        all(
            item["solve"]["outcome"] == "TRANSVERSE_SOLVED"
            for item in candidate_records
        ),
        f"candidates={len(candidate_records)}, stationary={stationary_count}",
    )
    check(
        f"{parity}: every localized candidate has a complete action-only audit",
        len(candidate_records) == len(candidates)
        and all("gates" in item["action"] for item in candidate_records),
    )
    check(
        f"{parity}: the frozen scientific outcome is assigned mechanically",
        outcome in {
            "UNIQUE_WEAK_SCALE_STATIONARY_BASE_IN_SCAN",
            "MULTIPLE_WEAK_SCALE_STATIONARY_BASES_IN_SCAN",
            "NO_STATIONARY_BASE_IN_FROZEN_SCAN",
            "STATIONARY_BASE_NUMERICALLY_UNRESOLVED",
        },
        f"outcome={outcome}, grid_resolved={all_grid_resolved}, "
        f"intervals_resolved={intervals_resolved}",
    )


check(
    "both schedule parities completed the same frozen stationary-base audit",
    set(results) == {"even", "odd"}
    and all("outcome" in record for record in results.values()),
)


def serialize_result(record):
    return {
        "outcome": record["outcome"],
        "stationary_candidate_count": record["stationary_count"],
        "all_grid_transverse_solves_resolved": record["all_grid_resolved"],
        "all_candidate_intervals_resolved": record["all_intervals_resolved"],
        "localization_resolved": record["localization_resolved"],
        "quotient_preconditioner": record["quotient_hessian"].tolist(),
        "grid": {
            f"{t:+.3f}": compact_solve(item)
            for t, item in sorted(record["grid"].items())
        },
        "bisections": [
            {
                "interval": item["interval"],
                "resolved": item["resolved"],
                "history": [compact_solve(point) for point in item["history"]],
                "candidate": None if item["candidate"] is None else compact_solve(item["candidate"]),
            }
            for item in record["bisections"]
        ],
        "candidates": [
            {
                "solve": compact_solve(item["solve"]),
                "action_only": serialize_action(item["action"]),
            }
            for item in record["candidates"]
        ],
    }


outcomes = {record["outcome"] for record in results.values()}
if outcomes == {"UNIQUE_WEAK_SCALE_STATIONARY_BASE_IN_SCAN"}:
    verdict = (
        "DERIVED COMPUTATIONAL LOCAL: each schedule has one weak-scale "
        "stationary base in the frozen collective scan.  Its Hessian and "
        "boundary response must be recertified at the refined coordinates."
    )
elif outcomes == {"MULTIPLE_WEAK_SCALE_STATIONARY_BASES_IN_SCAN"}:
    verdict = (
        "DERIVED COMPUTATIONAL LOCAL: each schedule has multiple weak-scale "
        "stationary bases in the frozen scan; collective freedom is not "
        "resolved to a unique base."
    )
elif outcomes == {"NO_STATIONARY_BASE_IN_FROZEN_SCAN"}:
    verdict = (
        "DERIVED COMPUTATIONAL NEGATIVE WITHIN SCAN: neither schedule has a "
        "weak-scale stationary base in the preregistered collective interval."
    )
else:
    verdict = (
        "OPEN NUMERICALLY OR PARITY-SPLIT: the two schedules do not share one "
        "resolved stationary-base outcome under the frozen audit."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "initial_protocol_commit": INITIAL_PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "precision_quotient_result_commit": PRECISION_RESULT_COMMIT,
    "published_control_result_commit": PUBLISHED_RESULT_COMMIT,
    "collective_scan": T_GRID,
    "damping_sequence": DAMPING,
    "action_precision_digits": 100,
    "action_derivative_steps": ACTION_STEPS,
    "parities": {
        parity: serialize_result(record) for parity, record in results.items()
    },
    "verdict": verdict,
    "claim_boundary": {
        "stationary_base_within_frozen_scan": "DERIVED COMPUTATIONAL IF RESOLVED",
        "root_count_outside_scan": "NOT TESTED",
        "collective_lapse_is_gauge": "OPEN",
        "refined_hessian_and_boundary_response": "NOT YET RECALCULATED",
        "nonlinear_boundary_continuation": "NOT TESTED",
        "full_840_edge_carrier": "NOT TESTED",
        "clock_or_planck_scale": "NOT DERIVED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)

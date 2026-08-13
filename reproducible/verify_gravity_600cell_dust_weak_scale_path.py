#!/usr/bin/env python3
"""100-decimal weak-scale gradient audit on the collective dust path.

Protocol commit: 8380f0d.  All thirteen points and both schedule parities are
tested directly from the complete action with smaller, preregistered steps;
no binary64 transverse root solver is used.
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
FIRST_INPUT = HERE / "gravity_600cell_dust_stationary_base.json"
OUTPUT = HERE / "gravity_600cell_dust_weak_scale_path.json"
PROTOCOL_COMMIT = "8380f0d"
COORDINATE_CORRECTION_COMMIT = "ad5f0ad"
PRECISION_RESULT_COMMIT = "29a779f"
T_GRID = (
    -0.10, -0.075, -0.05, -0.03, -0.02, -0.01,
    0.0,
    0.01, 0.02, 0.03, 0.05, 0.075, 0.10,
)
ACTION_STEPS = (2.0e-5, 1.0e-5, 5.0e-6)
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
    "the precision quotient and first stationary audit retain provenance",
    precision_input["protocol_commit"] == "da34272"
    and first_input["protocol_commit"] == "61b7f35"
    and precision_input["passed"] == precision_input["tests"] == 14
    and first_input["passed"] == first_input["tests"] == 12
    and set(precision_input["parities"])
        == set(first_input["parities"])
        == {"even", "odd"},
)


print("Loading the independently certified complete-action implementation...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_weak_path",
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
    "the frozen quotient basis remains orthonormal",
    quotient_basis.shape == (35, 34)
    and np.linalg.norm(
        quotient_basis.T@quotient_basis-np.eye(34), 2
    ) < 2e-14,
)


def float_path_variables(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    return np.concatenate((
        np.full(30, diagonal),
        np.full(5, rho),
        np.full(30, dust.L0_SQUARE),
    ))


def path_tangent(t):
    rho = dust.TAU_SQUARE*math.exp(float(t))
    diagonal = dust.L0_SQUARE-rho
    tangent = np.concatenate((
        np.full(30, -rho/diagonal), np.ones(5)
    ))
    return tangent/np.linalg.norm(tangent)


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


def arb_path_variables(t):
    t_arb = arb.mpf(str(t))
    rho = dust.ARB_TAU_SQUARE*arb.exp(t_arb)
    diagonal = dust.ARB_L0_SQUARE-rho
    return (
        [diagonal]*30
        + [rho]*5
        + [dust.ARB_L0_SQUARE]*30
    )


_ARB_MODEL = None


def initialize_action_worker(model):
    global _ARB_MODEL
    _ARB_MODEL = model
    configure_arb_precision()


def action_worker(point):
    t, coordinate, logarithmic_step = point
    variables = arb_path_variables(t)
    variables[coordinate] *= arb.exp(arb.mpf(str(logarithmic_step)))
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


fork_context = mp.get_context("fork")
results = {}
for parity, model in dust.bl.models.items():
    print(
        f"Auditing 13 weak-scale path points at 100 decimals: {parity} parity...",
        flush=True,
    )
    quotient_hessian = np.array(
        precision_input["parities"][parity]["quotient"]["matrix"],
        dtype=float,
    )
    quotient_hessian = (quotient_hessian+quotient_hessian.T)/2.0

    points = []
    branches = {}
    for t in T_GRID:
        branch_pass = True
        minimum_gram = math.inf
        minimum_argument = math.inf
        base = float_path_variables(t)
        for step in ACTION_STEPS:
            for coordinate in range(35):
                for sign in (1.0, -1.0):
                    logarithmic_step = sign*step
                    points.append((t, coordinate, logarithmic_step))
                    trial = base.copy()
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
        branches[t] = {
            "pass": branch_pass,
            "minimum_gram": minimum_gram,
            "minimum_argument": minimum_argument,
        }

    with fork_context.Pool(
        processes=8,
        initializer=initialize_action_worker,
        initargs=(model,),
    ) as pool:
        action_values = pool.map(action_worker, points, chunksize=1)

    cursor = 0
    path_records = {}
    for t in T_GRID:
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
        epsilon = float(np.linalg.norm(error_row))
        equation_norm = float(np.linalg.norm(equation))
        transverse_correction = float(np.linalg.norm(
            np.linalg.solve(quotient_hessian, quotient_basis.T@equation)
        ))
        transverse_error = float(np.linalg.norm(
            np.linalg.solve(quotient_hessian, quotient_basis.T@error_row)
        ))
        tangent = path_tangent(t)
        scalar = float(tangent@equation)
        scalar_error = float(abs(tangent@error_row))
        equation_floor = max(epsilon, 1e-30)
        transverse_floor = max(transverse_error, 1e-30)
        scalar_floor = max(scalar_error, 1e-30)
        maximum_imaginary = float(max(
            maximum_imaginary,
            max(abs(arb.im(value)) for value in sixth_order),
        ))
        numerical_pass = bool(
            branches[t]["pass"] and maximum_imaginary < 1e-80
        )
        adequacy_pass = bool(
            numerical_pass
            and transverse_error < 1e-5
            and transverse_correction < 1e-5
            and equation_norm < 1e-10
            and abs(scalar) < 1e-11
        )
        zero_consistent = bool(
            equation_norm <= 10*equation_floor
            and transverse_correction <= 10*transverse_floor
            and abs(scalar) <= 10*scalar_floor
        )
        resolved_nonzero = bool(
            equation_norm > 100*equation_floor
            or transverse_correction > 100*transverse_floor
            or abs(scalar) > 100*scalar_floor
        )
        resolved_nonstationary = bool(
            numerical_pass
            and transverse_correction > 1e-5
            and transverse_correction > 100*transverse_floor
        )
        if adequacy_pass and zero_consistent:
            label = "STATIONARY_WITHIN_ACTION_ERROR_AND_WEAK_SCALE"
        elif adequacy_pass and resolved_nonzero:
            label = "RESOLVED_SMALL_NONZERO_GRADIENT"
        elif resolved_nonstationary:
            label = "RESOLVED_NONSTATIONARY_ON_WEAK_SCALE"
        else:
            label = "WEAK_SCALE_PATH_NUMERICALLY_UNRESOLVED"

        path_records[t] = {
            "label": label,
            "derivatives": derivatives,
            "richardson12": richardson12,
            "richardson23": richardson23,
            "sixth_order": sixth_order,
            "equation": equation,
            "error_row": error_row,
            "equation_norm": equation_norm,
            "epsilon": epsilon,
            "equation_floor": equation_floor,
            "transverse_correction": transverse_correction,
            "transverse_error": transverse_error,
            "transverse_floor": transverse_floor,
            "scalar": scalar,
            "scalar_error": scalar_error,
            "scalar_floor": scalar_floor,
            "maximum_imaginary": maximum_imaginary,
            "branch": branches[t],
            "adequacy_pass": adequacy_pass,
            "zero_consistent": zero_consistent,
            "resolved_nonzero": resolved_nonzero,
        }

    labels = {record["label"] for record in path_records.values()}
    if labels == {"STATIONARY_WITHIN_ACTION_ERROR_AND_WEAK_SCALE"}:
        outcome = "ALL_13_PATH_POINTS_STATIONARY_WITHIN_ERROR"
    elif labels <= {
        "STATIONARY_WITHIN_ACTION_ERROR_AND_WEAK_SCALE",
        "RESOLVED_SMALL_NONZERO_GRADIENT",
    }:
        outcome = "ALL_13_PATH_POINTS_WEAK_SCALE_ADEQUATE_SOME_NONZERO"
    elif "RESOLVED_NONSTATIONARY_ON_WEAK_SCALE" in labels:
        outcome = "PATH_HAS_RESOLVED_NONSTATIONARY_POINT"
    else:
        outcome = "PATH_AUDIT_NUMERICALLY_UNRESOLVED"

    results[parity] = {
        "outcome": outcome,
        "quotient_hessian": quotient_hessian,
        "path": path_records,
    }
    check(
        f"{parity}: all 2730 action points were evaluated on certified branches",
        cursor == len(action_values)
        and len(action_values) == 2730
        and all(record["branch"]["pass"] for record in path_records.values()),
    )
    check(
        f"{parity}: all thirteen point labels follow the frozen gates",
        len(path_records) == 13
        and all(record["label"] in {
            "STATIONARY_WITHIN_ACTION_ERROR_AND_WEAK_SCALE",
            "RESOLVED_SMALL_NONZERO_GRADIENT",
            "RESOLVED_NONSTATIONARY_ON_WEAK_SCALE",
            "WEAK_SCALE_PATH_NUMERICALLY_UNRESOLVED",
        } for record in path_records.values()),
        f"labels={sorted(labels)}",
    )
    check(
        f"{parity}: the path outcome is assigned mechanically",
        outcome in {
            "ALL_13_PATH_POINTS_STATIONARY_WITHIN_ERROR",
            "ALL_13_PATH_POINTS_WEAK_SCALE_ADEQUATE_SOME_NONZERO",
            "PATH_HAS_RESOLVED_NONSTATIONARY_POINT",
            "PATH_AUDIT_NUMERICALLY_UNRESOLVED",
        },
        f"outcome={outcome}",
    )


parity_records = {}
for t in T_GRID:
    even = results["even"]["path"][t]
    odd = results["odd"]["path"][t]
    difference = float(np.linalg.norm(even["equation"]-odd["equation"]))
    envelope = even["epsilon"]+odd["epsilon"]
    floor = max(envelope, 1e-30)
    if difference <= 10*floor:
        label = "PARITY_AGREES_WITHIN_ACTION_ERROR"
    elif difference > 100*floor:
        label = "PARITY_RESOLVED_DIFFERENT"
    else:
        label = "PARITY_COMPARISON_UNRESOLVED"
    parity_records[t] = {
        "label": label,
        "difference_norm": difference,
        "combined_empirical_envelope": envelope,
        "difference_over_envelope": difference/floor,
    }

check(
    "all thirteen parity comparisons follow the frozen empirical envelopes",
    len(parity_records) == 13
    and all(record["label"] in {
        "PARITY_AGREES_WITHIN_ACTION_ERROR",
        "PARITY_RESOLVED_DIFFERENT",
        "PARITY_COMPARISON_UNRESOLVED",
    } for record in parity_records.values()),
)


old_even_errors = {
    float(candidate["solve"]["t"]): candidate["action_only"][
        "empirical_error_norm"
    ]
    for candidate in first_input["parities"]["even"]["candidates"]
}
error_ratios = np.array([
    old_even_errors[t]/results["even"]["path"][t]["epsilon"]
    for t in T_GRID
])
median_error_improvement = float(np.median(error_ratios))
check(
    "the smaller-step audit improves the shared empirical error by over 1e3",
    median_error_improvement > 1e3,
    f"median improvement={median_error_improvement:.6e}",
)


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 80),
        "imaginary": arb.nstr(arb.im(value), 80),
    }


def serialize_point(record):
    return {
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
        "sixth_order_float": record["equation"].tolist(),
        "error_row_float": record["error_row"].tolist(),
        "equation_norm": record["equation_norm"],
        "empirical_error_norm": record["epsilon"],
        "equation_floor": record["equation_floor"],
        "transverse_correction_norm": record["transverse_correction"],
        "transverse_error_norm": record["transverse_error"],
        "transverse_floor": record["transverse_floor"],
        "collective_scalar": record["scalar"],
        "collective_scalar_error": record["scalar_error"],
        "collective_scalar_floor": record["scalar_floor"],
        "maximum_action_or_derivative_imaginary": record[
            "maximum_imaginary"
        ],
        "minimum_absolute_gram_eigenvalue": record["branch"]["minimum_gram"],
        "minimum_angle_argument_modulus": record["branch"][
            "minimum_argument"
        ],
        "all_210_branches_pass": record["branch"]["pass"],
        "adequacy_pass": record["adequacy_pass"],
        "zero_consistent": record["zero_consistent"],
        "resolved_nonzero": record["resolved_nonzero"],
    }


outcomes = {record["outcome"] for record in results.values()}
parity_labels = {record["label"] for record in parity_records.values()}
if (
    outcomes == {"ALL_13_PATH_POINTS_STATIONARY_WITHIN_ERROR"}
    and parity_labels == {"PARITY_AGREES_WITHIN_ACTION_ERROR"}
):
    verdict = (
        "DERIVED COMPUTATIONAL ON THE FROZEN GRID: all thirteen collective "
        "points in both schedules are stationary within the smaller-step "
        "action error and weak-mode scale.  A continuous stationary family "
        "remains PATTERN, not an analytic interval theorem."
    )
elif outcomes <= {
    "ALL_13_PATH_POINTS_STATIONARY_WITHIN_ERROR",
    "ALL_13_PATH_POINTS_WEAK_SCALE_ADEQUATE_SOME_NONZERO",
}:
    verdict = (
        "DERIVED COMPUTATIONAL WEAK-SCALE ADEQUACY ON THE FROZEN GRID, but "
        "at least one point has a resolved small nonzero gradient or parity "
        "comparison not consistent with zero."
    )
elif "PATH_HAS_RESOLVED_NONSTATIONARY_POINT" in outcomes:
    verdict = (
        "DERIVED COMPUTATIONAL NEGATIVE: the collective path contains a "
        "point nonstationary on the four-soft-mode scale."
    )
else:
    verdict = (
        "OPEN NUMERICALLY: at least one collective point or parity remains "
        "unresolved under the frozen smaller-step action audit."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "coordinate_correction_commit": COORDINATE_CORRECTION_COMMIT,
    "precision_quotient_result_commit": PRECISION_RESULT_COMMIT,
    "precision_digits": 100,
    "collective_scan": T_GRID,
    "action_derivative_steps": ACTION_STEPS,
    "parities": {
        parity: {
            "outcome": record["outcome"],
            "quotient_preconditioner": record["quotient_hessian"].tolist(),
            "path": {
                f"{t:+.3f}": serialize_point(point)
                for t, point in record["path"].items()
            },
        }
        for parity, record in results.items()
    },
    "parity_comparisons": {
        f"{t:+.3f}": record for t, record in parity_records.items()
    },
    "old_to_new_even_error_ratios": error_ratios.tolist(),
    "median_old_to_new_error_improvement": median_error_improvement,
    "verdict": verdict,
    "claim_boundary": {
        "frozen_13_point_grid": "DERIVED COMPUTATIONAL IF RESOLVED",
        "continuous_stationary_family": "PATTERN ONLY",
        "analytic_gauge_identity": "NOT PROVED",
        "boundary_deformed_lapse_freedom": "NOT TESTED",
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

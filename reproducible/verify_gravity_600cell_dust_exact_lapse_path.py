#!/usr/bin/env python3
"""Direct 100-decimal audit of the exact collective lapse path.

Protocol commit: 515a509.  This does not alter the upstream FIVE_STIFF label;
it avoids the measured double-precision Schur-lift error by evaluating the
published path q=l0^2-rho itself.
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
SCHUR_INPUT = HERE / "gravity_600cell_dust_lapse_schur.json"
OUTPUT = HERE / "gravity_600cell_dust_exact_lapse_path.json"
PROTOCOL_COMMIT = "515a509"
SCHUR_RESULT_COMMIT = "dc927a5"
PRIOR_ART_COMMIT = "0882934"
PATH_POINTS = (
    -0.1, -0.03, -0.01, -0.003, -0.001,
    0.0,
    0.001, 0.003, 0.01, 0.03, 0.1,
)
DERIVATIVE_STEPS = (1.0e-2, 5.0e-3, 2.5e-3)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


schur_input = json.loads(SCHUR_INPUT.read_text())
check(
    "the frozen four-mode input is the committed 80-decimal Schur result",
    schur_input["protocol_commit"] == "5c0372a"
    and schur_input["precision_digits"] == 80
    and schur_input["passed"] == schur_input["tests"] == 13
    and {
        item["outcome"] for item in schur_input["parities"].values()
    } == {"FIVE_STIFF"},
)


print("Loading the independently certified complete-action implementation...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_exact_lapse",
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
    dust.ARB_MASS = (
        (90/arb.pi)*dust.ARB_EPSILON_3*dust.ARB_L0
    )
    dust.ARB_TAU_SQUARE = dust.ARB_TAU**2
    dust.ARB_SLANT_SQUARE = dust.ARB_L0_SQUARE-dust.ARB_TAU_SQUARE
    # arb_edge_square captured this list as a default argument at definition;
    # mutate it in place rather than merely rebinding the module name.
    dust.ARB_OLD_VALUES[:] = [dust.ARB_L0_SQUARE]*30


configure_arb_precision()


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


def action_worker(t):
    return dust.arb_action_components(_ARB_MODEL, arb_path_variables(t))[2]


fork_context = mp.get_context("fork")
results = {}
all_action_points = sorted(set(
    PATH_POINTS
    + tuple(step for step in DERIVATIVE_STEPS)
    + tuple(-step for step in DERIVATIVE_STEPS)
))

for parity, model in dust.bl.models.items():
    print(f"Evaluating exact lapse path at 100 decimals: {parity} parity...", flush=True)
    with fork_context.Pool(
        processes=8,
        initializer=initialize_action_worker,
        initargs=(model,),
    ) as pool:
        action_values = pool.map(action_worker, all_action_points, chunksize=1)
    actions = dict(zip(all_action_points, action_values))
    base_action = actions[0.0]

    path_action_records = {}
    maximum_action_error = 0.0
    maximum_action_imaginary = arb.mpf(0)
    for t in PATH_POINTS:
        value = actions[t]
        difference = value-base_action
        error = float(
            abs(difference)/max(arb.mpf(1), abs(value), abs(base_action))
        )
        maximum_action_error = max(maximum_action_error, error)
        maximum_action_imaginary = max(
            maximum_action_imaginary, abs(arb.im(value))
        )
        path_action_records[f"{t:+.3f}"] = {
            "action": value,
            "difference": difference,
            "normalized_difference": error,
        }

    first = []
    second = []
    for step in DERIVATIVE_STEPS:
        step_arb = arb.mpf(str(step))
        plus, minus = actions[step], actions[-step]
        first.append((plus-minus)/(2*step_arb))
        second.append((plus-2*base_action+minus)/(step_arb**2))
    first_coarse = (4*first[1]-first[0])/3
    first_fine = (4*first[2]-first[1])/3
    first_sixth = (16*first_fine-first_coarse)/15
    second_coarse = (4*second[1]-second[0])/3
    second_fine = (4*second[2]-second[1])/3
    second_sixth = (16*second_fine-second_coarse)/15

    # Full transverse stationarity and branch audit at the eleven frozen path
    # points uses the already certified analytic gradient implementation.
    residual_records = {}
    maximum_residual = 0.0
    branch_pass = True
    minimum_gram = math.inf
    minimum_argument = math.inf
    for t in PATH_POINTS:
        rho = dust.TAU_SQUARE*math.exp(t)
        diagonal = dust.L0_SQUARE-rho
        variables = np.concatenate((
            np.full(30, diagonal),
            np.full(5, rho),
            np.full(30, dust.L0_SQUARE),
        ))
        _, gradient, _, data = dust.total_reduced_evaluation(
            model, variables, dust.old_values
        )
        residual = gradient[:35]/24.0
        maximum = float(np.max(np.abs(residual)))
        maximum_residual = max(maximum_residual, maximum)
        minimum_gram = min(minimum_gram, float(data["minimum_gram"]))
        minimum_argument = min(
            minimum_argument, float(data["minimum_argument"])
        )
        lorentzian = dict(data["negative_counts"]) == {1: 100}
        branch_pass &= bool(
            lorentzian
            and data["minimum_gram"] > 1.0e-8
            and data["minimum_argument"] > 1.0e-6
        )
        residual_records[f"{t:+.3f}"] = {
            "rho": rho,
            "diagonal": diagonal,
            "residuals_real": residual.real,
            "maximum_absolute": maximum,
            "maximum_imaginary": float(np.max(np.abs(residual.imag))),
            "minimum_gram": float(data["minimum_gram"]),
            "minimum_argument": float(data["minimum_argument"]),
            "lorentzian": lorentzian,
        }

    schur_record = schur_input["parities"][parity]
    relative = np.array(schur_record["relative_eigenvalues"], dtype=float)
    epsilon = float(schur_record["empirical_error_norm"])
    relative_pass = bool(np.min(np.abs(relative)) > 100.0*epsilon)
    action_constant = maximum_action_error <= 1.0e-50
    tangent_null = abs(second_sixth) <= arb.mpf("1e-40")
    family_stationary = maximum_residual <= 1.0e-7
    precision_pass = maximum_action_imaginary < arb.mpf("1e-80")

    positive = bool(
        action_constant
        and tangent_null
        and family_stationary
        and branch_pass
        and precision_pass
        and relative_pass
    )
    tangential_only = bool(
        action_constant and tangent_null and not family_stationary
    )
    not_null = bool(
        branch_pass
        and precision_pass
        and (not action_constant or not tangent_null)
    )
    if positive:
        outcome = "ONE_COLLECTIVE_LAPSE_NULL_FOUR_PSEUDOCONSTRAINT_STIFF"
    elif tangential_only:
        outcome = "TANGENTIAL_NULL_NOT_SOLUTION_FAMILY"
    elif not_null:
        outcome = "COLLECTIVE_DIRECTION_NOT_NULL"
    else:
        outcome = "NUMERICALLY_UNRESOLVED"

    results[parity] = {
        "path_actions": path_action_records,
        "maximum_action_error": maximum_action_error,
        "maximum_action_imaginary": maximum_action_imaginary,
        "first_centered": first,
        "first_coarse": first_coarse,
        "first_fine": first_fine,
        "first_sixth": first_sixth,
        "second_centered": second,
        "second_coarse": second_coarse,
        "second_fine": second_fine,
        "second_sixth": second_sixth,
        "residuals": residual_records,
        "maximum_residual": maximum_residual,
        "branch_pass": branch_pass,
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "relative_eigenvalues": relative,
        "relative_epsilon": epsilon,
        "relative_pass": relative_pass,
        "action_constant": action_constant,
        "tangent_null": bool(tangent_null),
        "family_stationary": family_stationary,
        "precision_pass": precision_pass,
        "outcome": outcome,
    }

    check(
        f"{parity}: all eleven exact-path geometries remain Lorentzian and off branch boundaries",
        branch_pass,
        f"min Gram={minimum_gram:.3e}, min argument={minimum_argument:.3e}",
    )
    check(
        f"{parity}: the 100-decimal path actions are finite and real",
        precision_pass
        and all(math.isfinite(item["normalized_difference"])
                for item in path_action_records.values()),
        f"max action imaginary={float(maximum_action_imaginary):.3e}",
    )
    check(
        f"{parity}: all eleven full residual vectors were evaluated",
        len(residual_records) == 11
        and all(math.isfinite(item["maximum_absolute"])
                for item in residual_records.values()),
        f"maximum per-edge residual={maximum_residual:.3e}",
    )
    check(
        f"{parity}: the exact-path outcome is assigned mechanically",
        outcome in {
            "ONE_COLLECTIVE_LAPSE_NULL_FOUR_PSEUDOCONSTRAINT_STIFF",
            "TANGENTIAL_NULL_NOT_SOLUTION_FAMILY",
            "COLLECTIVE_DIRECTION_NOT_NULL",
            "NUMERICALLY_UNRESOLVED",
        },
        f"outcome={outcome}, action error={maximum_action_error:.3e}, "
        f"curvature={arb.nstr(second_sixth, 8)}, relative min={np.min(np.abs(relative)):.3e}",
    )


check(
    "both schedule parities completed the same frozen exact-path audit",
    set(results) == {"even", "odd"}
    and all(len(record["path_actions"]) == 11 for record in results.values()),
)


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 80),
        "imaginary": arb.nstr(arb.im(value), 80),
    }


def serialize_result(record):
    return {
        "outcome": record["outcome"],
        "path_actions": {
            t: {
                "action": serialize_complex(item["action"]),
                "difference_from_base": serialize_complex(item["difference"]),
                "normalized_difference": item["normalized_difference"],
            }
            for t, item in record["path_actions"].items()
        },
        "maximum_normalized_action_difference": record["maximum_action_error"],
        "maximum_action_imaginary": arb.nstr(
            record["maximum_action_imaginary"], 80
        ),
        "first_derivative": {
            "centered": [serialize_complex(value) for value in record["first_centered"]],
            "richardson_coarse": serialize_complex(record["first_coarse"]),
            "richardson_fine": serialize_complex(record["first_fine"]),
            "sixth_order": serialize_complex(record["first_sixth"]),
        },
        "second_derivative": {
            "centered": [serialize_complex(value) for value in record["second_centered"]],
            "richardson_coarse": serialize_complex(record["second_coarse"]),
            "richardson_fine": serialize_complex(record["second_fine"]),
            "sixth_order": serialize_complex(record["second_sixth"]),
        },
        "full_residuals": {
            t: {
                **{key: value for key, value in item.items()
                   if key != "residuals_real"},
                "residuals_real": item["residuals_real"].tolist(),
            }
            for t, item in record["residuals"].items()
        },
        "maximum_per_edge_residual": record["maximum_residual"],
        "branch_pass": record["branch_pass"],
        "minimum_absolute_gram_eigenvalue": record["minimum_gram"],
        "minimum_angle_argument_modulus": record["minimum_argument"],
        "recorded_relative_eigenvalues": record["relative_eigenvalues"].tolist(),
        "recorded_relative_empirical_error": record["relative_epsilon"],
        "relative_modes_pass_separation": record["relative_pass"],
        "gates": {
            "action_constant_1e-50": record["action_constant"],
            "exact_path_curvature_1e-40": record["tangent_null"],
            "full_family_stationary_1e-7": record["family_stationary"],
            "precision_imaginary_1e-80": record["precision_pass"],
        },
    }


outcomes = {record["outcome"] for record in results.values()}
if outcomes == {
    "ONE_COLLECTIVE_LAPSE_NULL_FOUR_PSEUDOCONSTRAINT_STIFF"
}:
    verdict = (
        "DERIVED COMPUTATIONAL: the exact time-symmetric path is a stationary "
        "collective lapse family in both parities, while four relative "
        "phase-lapse pseudo-constraint curvatures remain resolved."
    )
elif "COLLECTIVE_DIRECTION_NOT_NULL" in outcomes:
    verdict = (
        "DERIVED COMPUTATIONAL NEGATIVE: the exact published lapse path is "
        "not action-null under the frozen 100-decimal gate."
    )
else:
    verdict = (
        "OPEN: the exact-path audit does not support the complete scoped "
        "1+4 classification in both parities."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "schur_result_commit": SCHUR_RESULT_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "precision_digits": 100,
    "path_points": PATH_POINTS,
    "derivative_steps": DERIVATIVE_STEPS,
    "parities": {
        parity: serialize_result(record) for parity, record in results.items()
    },
    "verdict": verdict,
    "claim_boundary": {
        "collective_lapse_family": "DERIVED COMPUTATIONAL IF POSITIVE",
        "four_relative_pseudo_constraints": "KNOWN MECHANISM, EXPLICITLY COMPUTED",
        "exact_analytic_identity": "NOT PROVED",
        "full_carrier_or_continuum": "NOT TESTED",
        "physical_clock_or_tick_selection": "REFUTED ON THIS FAMILY",
        "multi_tick_evolution": "NOT TESTED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)

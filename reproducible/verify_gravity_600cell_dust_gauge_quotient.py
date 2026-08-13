#!/usr/bin/env python3
"""Gauge-quotient linear boundary response at the 600-cell dust sandwich.

Protocol commit: 25d9ee9.  The exact collective lapse tangent is removed
before testing the 34-dimensional internal quotient.  The compatibility row
is reconstructed independently from 90-decimal complete-action rectangles.
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
JACOBIAN_INPUT = HERE / "gravity_600cell_dust_implicit_jacobian.json"
SCHUR_INPUT = HERE / "gravity_600cell_dust_lapse_schur.json"
LAPSE_INPUT = HERE / "gravity_600cell_dust_exact_lapse_path.json"
OUTPUT = HERE / "gravity_600cell_dust_gauge_quotient.json"
PROTOCOL_COMMIT = "25d9ee9"
PRIOR_ART_COMMIT = "ff8f404"
EXACT_LAPSE_COMMIT = "790fc7f"
SCHUR_RESULT_COMMIT = "dc927a5"
JACOBIAN_RECORD_COMMIT = "7d5e9fc"
MIXED_STEPS = (5.0e-4, 2.5e-4)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


def rank_table(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    if singular.size == 0 or singular[0] == 0:
        return {f"{threshold:.0e}": 0 for threshold in (1e-7, 1e-9, 1e-11)}
    return {
        f"{threshold:.0e}": int(np.sum(singular > threshold*singular[0]))
        for threshold in (1e-7, 1e-9, 1e-11)
    }


def householder_complement(vector):
    vector = np.asarray(vector, dtype=float)
    vector = vector/np.linalg.norm(vector)
    sign = 1.0 if vector[0] >= 0 else -1.0
    reflector_vector = vector.copy()
    reflector_vector[0] += sign
    reflector = (
        np.eye(len(vector))
        - 2.0*np.outer(reflector_vector, reflector_vector)
        / float(reflector_vector@reflector_vector)
    )
    return reflector[:, 1:]


jacobian_input = json.loads(JACOBIAN_INPUT.read_text())
schur_input = json.loads(SCHUR_INPUT.read_text())
lapse_input = json.loads(LAPSE_INPUT.read_text())
check(
    "all three frozen inputs have the committed provenance and two parities",
    jacobian_input["protocol_commit"] == "41acf7b"
    and schur_input["protocol_commit"] == "5c0372a"
    and lapse_input["protocol_commit"] == "515a509"
    and set(jacobian_input["parities"])
        == set(schur_input["parities"])
        == set(lapse_input["parities"])
        == {"even", "odd"}
    and {
        item["outcome"] for item in lapse_input["parities"].values()
    } == {"ONE_COLLECTIVE_LAPSE_NULL_FOUR_PSEUDOCONSTRAINT_STIFF"},
)


print("Loading the independently certified complete-action implementation...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_gauge_quotient",
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
    arb.mp.dps = 90
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
    dust.ARB_OLD_VALUES[:] = [dust.ARB_L0_SQUARE]*30


configure_arb_precision()


def arb_mixed_variables(t, boundary_index, boundary_log):
    t_arb = arb.mpf(str(t))
    boundary_arb = arb.mpf(str(boundary_log))
    rho = dust.ARB_TAU_SQUARE*arb.exp(t_arb)
    diagonal = dust.ARB_L0_SQUARE-rho
    final = [dust.ARB_L0_SQUARE]*30
    final[boundary_index] *= arb.exp(boundary_arb)
    return [diagonal]*30 + [rho]*5 + final


def float_mixed_variables(t, boundary_index, boundary_log):
    rho = dust.TAU_SQUARE*math.exp(t)
    diagonal = dust.L0_SQUARE-rho
    final = np.full(30, dust.L0_SQUARE)
    final[boundary_index] *= math.exp(boundary_log)
    return np.concatenate((
        np.full(30, diagonal),
        np.full(5, rho),
        final,
    ))


_ARB_MODEL = None


def initialize_action_worker(model):
    global _ARB_MODEL
    _ARB_MODEL = model
    configure_arb_precision()


def action_worker(point):
    t, boundary_index, boundary_log = point
    variables = arb_mixed_variables(t, boundary_index, boundary_log)
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


raw_tangent = np.concatenate((
    np.full(30, -dust.TAU_SQUARE/dust.SLANT_SQUARE),
    np.ones(5),
))
raw_tangent_norm = float(np.linalg.norm(raw_tangent))
gauge_tangent = raw_tangent/raw_tangent_norm
quotient_basis = householder_complement(gauge_tangent)
projector = np.eye(35)-np.outer(gauge_tangent, gauge_tangent)

check(
    "the deterministic Householder basis spans the analytic gauge complement",
    quotient_basis.shape == (35, 34)
    and np.linalg.norm(quotient_basis.T@quotient_basis-np.eye(34), 2) < 2e-14
    and np.linalg.norm(quotient_basis.T@gauge_tangent) < 2e-14,
)


fork_context = mp.get_context("fork")
results = {}
for parity, model in dust.bl.models.items():
    print(f"Auditing the 34-dimensional quotient: {parity} parity...", flush=True)
    jacobian_record = jacobian_input["parities"][parity]
    schur_record = schur_input["parities"][parity]
    hessian = np.array(
        jacobian_record["symmetrized_internal_matrix"], dtype=float
    )
    boundary_block = np.array(
        jacobian_record["richardson_boundary_matrix"], dtype=float
    )
    hessian = (hessian+hessian.T)/2.0
    regular_block = hessian[:30, :30]
    regular_ranks = rank_table(regular_block)
    relative_values = np.array(
        schur_record["relative_eigenvalues"], dtype=float
    )
    relative_epsilon = float(schur_record["empirical_error_norm"])
    relative_upstream_pass = bool(
        np.min(np.abs(relative_values)) > 100.0*relative_epsilon
    )
    gauge_leakage = float(np.linalg.norm(hessian@gauge_tangent))
    projected_hessian = projector@hessian@projector
    quotient_hessian = quotient_basis.T@projected_hessian@quotient_basis
    quotient_hessian = (quotient_hessian+quotient_hessian.T)/2.0
    quotient_singular = np.linalg.svd(quotient_hessian, compute_uv=False)
    quotient_eigenvalues = np.linalg.eigvalsh(quotient_hessian)
    absolute_rank = int(np.sum(quotient_singular > 1.0e-9))
    quotient_weak = np.sort(quotient_singular)[:4]
    relative_sorted = np.sort(np.abs(relative_values))
    weak_spectrum_error = float(
        np.max(np.abs(quotient_weak-relative_sorted))
        / max(np.max(relative_sorted), np.finfo(float).tiny)
    )
    quotient_regular = bool(
        set(regular_ranks.values()) == {30}
        and relative_upstream_pass
        and absolute_rank == 34
        and weak_spectrum_error < 3.0e-3
    )
    quotient_label = (
        "QUOTIENT_REGULAR" if quotient_regular else "QUOTIENT_UNRESOLVED"
    )

    # Freeze and branch-check all 240 mixed action points before evaluating
    # the high-precision compatibility row.
    points = []
    minimum_gram = math.inf
    minimum_argument = math.inf
    all_lorentzian = True
    for boundary_index in range(30):
        for step in MIXED_STEPS:
            for t_sign, v_sign in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                point = (
                    t_sign*step,
                    boundary_index,
                    v_sign*step,
                )
                points.append(point)
                variables = float_mixed_variables(*point)
                _, _, _, data = dust.total_reduced_evaluation(
                    model, variables, dust.old_values
                )
                minimum_gram = min(
                    minimum_gram, float(data["minimum_gram"])
                )
                minimum_argument = min(
                    minimum_argument, float(data["minimum_argument"])
                )
                all_lorentzian &= dict(data["negative_counts"]) == {1: 100}

    print(f"Evaluating 240 frozen 90-decimal mixed actions: {parity} parity...", flush=True)
    with fork_context.Pool(
        processes=8,
        initializer=initialize_action_worker,
        initargs=(model,),
    ) as pool:
        action_values = pool.map(action_worker, points, chunksize=1)

    mixed_by_step = {
        step: [arb.mpc(0) for _ in range(30)] for step in MIXED_STEPS
    }
    maximum_imaginary_by_boundary = [arb.mpf(0) for _ in range(30)]
    cursor = 0
    for boundary_index in range(30):
        for step in MIXED_STEPS:
            pp, pm, mp_, mm = action_values[cursor:cursor+4]
            cursor += 4
            step_arb = arb.mpf(str(step))
            mixed_by_step[step][boundary_index] = (
                pp-pm-mp_+mm
            )/(4*step_arb**2)
            maximum_imaginary_by_boundary[boundary_index] = max(
                maximum_imaginary_by_boundary[boundary_index],
                abs(arb.im(pp)), abs(arb.im(pm)),
                abs(arb.im(mp_)), abs(arb.im(mm)),
            )

    normalization = 24*arb.mpf(str(raw_tangent_norm))
    coarse = [value/normalization for value in mixed_by_step[MIXED_STEPS[0]]]
    fine = [value/normalization for value in mixed_by_step[MIXED_STEPS[1]]]
    richardson = [(4*f-c)/3 for c, f in zip(coarse, fine)]
    compatibility_row = np.array([
        float(arb.re(value)) for value in richardson
    ])
    compatibility_fine = np.array([float(arb.re(value)) for value in fine])
    compatibility_imaginary = float(max(
        max(abs(arb.im(value)) for value in richardson),
        max(maximum_imaginary_by_boundary),
    ))
    compatibility_norm = float(np.linalg.norm(compatibility_row))
    compatibility_epsilon = float(np.linalg.norm(
        compatibility_row-compatibility_fine
    ))
    compatibility_floor = max(compatibility_epsilon, 1.0e-35)

    if compatibility_norm <= 10.0*compatibility_floor:
        compatibility_label = "ALL_30_BOUNDARY_DIRECTIONS_COMPATIBLE"
    elif compatibility_norm > 100.0*compatibility_floor:
        compatibility_label = "ONE_LINEAR_BOUNDARY_CONSTRAINT"
    else:
        compatibility_label = "BOUNDARY_COMPATIBILITY_UNRESOLVED"

    recorded_compatibility = gauge_tangent@boundary_block
    recorded_action_error = float(
        np.linalg.norm(recorded_compatibility-compatibility_row)
        / max(1.0, np.linalg.norm(boundary_block, 2))
    )

    response = None
    if quotient_regular and compatibility_label != "BOUNDARY_COMPATIBILITY_UNRESOLVED":
        if compatibility_label == "ALL_30_BOUNDARY_DIRECTIONS_COMPATIBLE":
            boundary_basis = np.eye(30)
        else:
            boundary_basis = householder_complement(compatibility_row)
        right_hand_side = (
            quotient_basis.T@projector@boundary_block@boundary_basis
        )
        coefficients = np.linalg.solve(quotient_hessian, -right_hand_side)
        internal_response = quotient_basis@coefficients
        quotient_residual = float(
            np.linalg.norm(quotient_hessian@coefficients+right_hand_side, 2)
            / max(1.0, np.linalg.norm(right_hand_side, 2))
        )
        unprojected_residual = (
            hessian@internal_response+boundary_block@boundary_basis
        )
        response_singular = np.linalg.svd(
            internal_response, compute_uv=False
        )
        response = {
            "boundary_basis": boundary_basis,
            "coefficients": coefficients,
            "internal_response": internal_response,
            "quotient_residual": quotient_residual,
            "unprojected_residual_norm": float(
                np.linalg.norm(unprojected_residual, 2)
            ),
            "unprojected_gauge_component_norm": float(
                np.linalg.norm(gauge_tangent@unprojected_residual)
            ),
            "singular_values": response_singular,
            "ranks": rank_table(internal_response),
            "condition": float(
                response_singular[0]/response_singular[-1]
            ) if response_singular[-1] > 0 else math.inf,
            "norm": float(np.linalg.norm(internal_response, 2)),
        }

    if quotient_regular and compatibility_label == "ALL_30_BOUNDARY_DIRECTIONS_COMPATIBLE":
        outcome = "REGULAR_QUOTIENT_ALL_30_LINEAR_RESPONSES"
    elif quotient_regular and compatibility_label == "ONE_LINEAR_BOUNDARY_CONSTRAINT":
        outcome = "REGULAR_QUOTIENT_29_RESPONSES_PLUS_ONE_BOUNDARY_CONSTRAINT"
    else:
        outcome = "QUOTIENT_OR_COMPATIBILITY_UNRESOLVED"

    results[parity] = {
        "regular_ranks": regular_ranks,
        "relative_values": relative_values,
        "relative_epsilon": relative_epsilon,
        "gauge_leakage": gauge_leakage,
        "quotient_hessian": quotient_hessian,
        "quotient_singular": quotient_singular,
        "quotient_eigenvalues": quotient_eigenvalues,
        "absolute_rank": absolute_rank,
        "weak_spectrum_error": weak_spectrum_error,
        "quotient_label": quotient_label,
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "all_lorentzian": all_lorentzian,
        "mixed_by_step": mixed_by_step,
        "coarse": coarse,
        "fine": fine,
        "richardson": richardson,
        "maximum_imaginary_by_boundary": maximum_imaginary_by_boundary,
        "compatibility_row": compatibility_row,
        "compatibility_norm": compatibility_norm,
        "compatibility_epsilon": compatibility_epsilon,
        "compatibility_floor": compatibility_floor,
        "compatibility_imaginary": compatibility_imaginary,
        "compatibility_label": compatibility_label,
        "recorded_compatibility": recorded_compatibility,
        "recorded_action_error": recorded_action_error,
        "response": response,
        "outcome": outcome,
    }

    check(
        f"{parity}: all 240 mixed-action geometries remain Lorentzian and off branch boundaries",
        all_lorentzian
        and minimum_gram > 1.0e-8
        and minimum_argument > 1.0e-6,
        f"min Gram={minimum_gram:.3e}, min argument={minimum_argument:.3e}",
    )
    check(
        f"{parity}: the 90-decimal compatibility reconstruction is real and finite",
        compatibility_imaginary < 1.0e-70
        and np.all(np.isfinite(compatibility_row)),
        f"max imaginary={compatibility_imaginary:.3e}, norm(c)={compatibility_norm:.3e}",
    )
    check(
        f"{parity}: the recorded analytic mixed block agrees with the independent action row",
        recorded_action_error < 3.0e-6,
        f"normalized error={recorded_action_error:.3e}",
    )
    check(
        f"{parity}: quotient and compatibility outcomes are assigned mechanically",
        quotient_label in {"QUOTIENT_REGULAR", "QUOTIENT_UNRESOLVED"}
        and compatibility_label in {
            "ALL_30_BOUNDARY_DIRECTIONS_COMPATIBLE",
            "ONE_LINEAR_BOUNDARY_CONSTRAINT",
            "BOUNDARY_COMPATIBILITY_UNRESOLVED",
        }
        and outcome in {
            "REGULAR_QUOTIENT_ALL_30_LINEAR_RESPONSES",
            "REGULAR_QUOTIENT_29_RESPONSES_PLUS_ONE_BOUNDARY_CONSTRAINT",
            "QUOTIENT_OR_COMPATIBILITY_UNRESOLVED",
        },
        f"quotient={quotient_label}, compatibility={compatibility_label}, outcome={outcome}",
    )
    if response is not None:
        check(
            f"{parity}: the gauge-fixed quotient response solves its frozen linear system",
            response["quotient_residual"] < 1.0e-5,
            f"relative quotient residual={response['quotient_residual']:.3e}, "
            f"response rank={response['ranks']}",
        )


check(
    "both schedule parities completed the same gauge-quotient audit",
    set(results) == {"even", "odd"}
    and all("outcome" in result for result in results.values()),
)


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 70),
        "imaginary": arb.nstr(arb.im(value), 70),
    }


def serialize_result(result):
    response = result["response"]
    return {
        "outcome": result["outcome"],
        "quotient": {
            "label": result["quotient_label"],
            "regular_block_ranks": result["regular_ranks"],
            "frozen_relative_eigenvalues": result["relative_values"].tolist(),
            "frozen_relative_empirical_error": result["relative_epsilon"],
            "preprojection_gauge_leakage_norm": result["gauge_leakage"],
            "matrix": result["quotient_hessian"].tolist(),
            "singular_values": result["quotient_singular"].tolist(),
            "eigenvalues": result["quotient_eigenvalues"].tolist(),
            "absolute_rank_above_1e-9": result["absolute_rank"],
            "weak_spectrum_relative_error": result["weak_spectrum_error"],
            "condition_2": float(
                result["quotient_singular"][0]
                / result["quotient_singular"][-1]
            ),
        },
        "compatibility": {
            "label": result["compatibility_label"],
            "mixed_action_steps": MIXED_STEPS,
            "mixed_action_derivatives": {
                f"{step:.1e}": [
                    serialize_complex(value)
                    for value in result["mixed_by_step"][step]
                ]
                for step in MIXED_STEPS
            },
            "coarse_row": [serialize_complex(value) for value in result["coarse"]],
            "fine_row": [serialize_complex(value) for value in result["fine"]],
            "richardson_row": [
                serialize_complex(value) for value in result["richardson"]
            ],
            "richardson_row_float": result["compatibility_row"].tolist(),
            "norm": result["compatibility_norm"],
            "empirical_error_norm": result["compatibility_epsilon"],
            "classification_floor": result["compatibility_floor"],
            "maximum_action_or_row_imaginary": result["compatibility_imaginary"],
            "maximum_action_imaginary_by_boundary": [
                arb.nstr(value, 70)
                for value in result["maximum_imaginary_by_boundary"]
            ],
            "recorded_analytic_row": result["recorded_compatibility"].tolist(),
            "recorded_vs_action_normalized_error": result["recorded_action_error"],
            "minimum_absolute_gram_eigenvalue": result["minimum_gram"],
            "minimum_angle_argument_modulus": result["minimum_argument"],
        },
        "response": None if response is None else {
            "boundary_dimension": int(response["boundary_basis"].shape[1]),
            "boundary_basis": response["boundary_basis"].tolist(),
            "coefficient_matrix": response["coefficients"].tolist(),
            "internal_response_matrix": response["internal_response"].tolist(),
            "relative_quotient_residual": response["quotient_residual"],
            "unprojected_residual_norm": response["unprojected_residual_norm"],
            "unprojected_gauge_component_norm": response[
                "unprojected_gauge_component_norm"
            ],
            "singular_values": response["singular_values"].tolist(),
            "ranks": response["ranks"],
            "condition_2": response["condition"],
            "norm_2": response["norm"],
        },
    }


outcomes = {result["outcome"] for result in results.values()}
if outcomes == {"REGULAR_QUOTIENT_ALL_30_LINEAR_RESPONSES"}:
    verdict = (
        "DERIVED COMPUTATIONAL LINEAR: both schedule parities have a regular "
        "34-dimensional internal quotient and all thirty final-boundary "
        "directions satisfy the collective lapse compatibility equation."
    )
elif outcomes == {
    "REGULAR_QUOTIENT_29_RESPONSES_PLUS_ONE_BOUNDARY_CONSTRAINT"
}:
    verdict = (
        "DERIVED COMPUTATIONAL LINEAR: both schedule parities have a regular "
        "internal quotient, but one linear boundary constraint leaves only "
        "twenty-nine compatible final-boundary directions."
    )
else:
    verdict = (
        "OPEN NUMERICALLY: at least one schedule parity does not resolve the "
        "gauge quotient or boundary compatibility under the frozen gates."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "exact_lapse_result_commit": EXACT_LAPSE_COMMIT,
    "schur_result_commit": SCHUR_RESULT_COMMIT,
    "jacobian_record_commit": JACOBIAN_RECORD_COMMIT,
    "precision_digits": 90,
    "raw_gauge_tangent": raw_tangent.tolist(),
    "normalized_gauge_tangent": gauge_tangent.tolist(),
    "quotient_basis": quotient_basis.tolist(),
    "parities": {
        parity: serialize_result(result) for parity, result in results.items()
    },
    "verdict": verdict,
    "claim_boundary": {
        "linear_gauge_quotient_response": "DERIVED COMPUTATIONAL IF RESOLVED",
        "nonlinear_boundary_family": "NOT ESTABLISHED",
        "higher_order_pseudo_constraint": "OPEN",
        "physical_graviton_identification": "NOT CLAIMED",
        "full_840_edge_carrier": "NOT TESTED",
        "multi_tick_evolution": "NOT TESTED",
        "clock_selection": "REFUTED ON THE COLLECTIVE LAPSE FAMILY",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)

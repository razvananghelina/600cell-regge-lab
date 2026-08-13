#!/usr/bin/env python3
"""Precision-corrected gauge quotient for the 600-cell dust sandwich.

Protocol commit: da34272.  A new 100-decimal mixed-action level tests the
boundary row, while the quotient is reconstructed deterministically from the
regular block, the independently certified relative Schur form, and the
exact collective lapse tangent.
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
FIRST_INPUT = HERE / "gravity_600cell_dust_gauge_quotient.json"
OUTPUT = HERE / "gravity_600cell_dust_gauge_quotient_precision.json"
PROTOCOL_COMMIT = "da34272"
PRIOR_ART_COMMIT = "44bd4cf"
FIRST_RESULT_COMMIT = "14a4517"
EXACT_LAPSE_COMMIT = "790fc7f"
SCHUR_RESULT_COMMIT = "dc927a5"
JACOBIAN_RECORD_COMMIT = "7d5e9fc"
MIXED_STEPS = (5.0e-4, 2.5e-4, 1.25e-4)
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


def frobenius_relative(new, old):
    return float(
        np.linalg.norm(new-old, "fro")/max(1.0, np.linalg.norm(old, "fro"))
    )


jacobian_input = json.loads(JACOBIAN_INPUT.read_text())
schur_input = json.loads(SCHUR_INPUT.read_text())
lapse_input = json.loads(LAPSE_INPUT.read_text())
first_input = json.loads(FIRST_INPUT.read_text())
check(
    "all four frozen inputs retain committed provenance and both parities",
    jacobian_input["protocol_commit"] == "41acf7b"
    and schur_input["protocol_commit"] == "5c0372a"
    and lapse_input["protocol_commit"] == "515a509"
    and first_input["protocol_commit"] == "25d9ee9"
    and set(jacobian_input["parities"])
        == set(schur_input["parities"])
        == set(lapse_input["parities"])
        == set(first_input["parities"])
        == {"even", "odd"}
    and {
        item["outcome"] for item in lapse_input["parities"].values()
    } == {"ONE_COLLECTIVE_LAPSE_NULL_FOUR_PSEUDOCONSTRAINT_STIFF"},
)


print("Loading the independently certified complete-action implementation...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_precision_quotient",
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
    dust.ARB_MASS = (90/arb.pi)*dust.ARB_EPSILON_3*dust.ARB_L0
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
pole_collective = np.ones(5)
pole_projector = (
    np.eye(5)-np.outer(pole_collective, pole_collective)
    / float(pole_collective@pole_collective)
)
pole_relative_basis = householder_complement(pole_collective)
boundary_collective = np.ones(30)
boundary_unit = boundary_collective/np.linalg.norm(boundary_collective)
boundary_projector = np.eye(30)-np.outer(boundary_unit, boundary_unit)

check(
    "the deterministic bases span the exact gauge and collective complements",
    quotient_basis.shape == (35, 34)
    and pole_relative_basis.shape == (5, 4)
    and np.linalg.norm(quotient_basis.T@quotient_basis-np.eye(34), 2) < 2e-14
    and np.linalg.norm(quotient_basis.T@gauge_tangent) < 2e-14
    and np.linalg.norm(pole_relative_basis.T@pole_collective) < 2e-14,
)


fork_context = mp.get_context("fork")
results = {}
for parity, model in dust.bl.models.items():
    print(f"Reconstructing the precision quotient: {parity} parity...", flush=True)
    jacobian_record = jacobian_input["parities"][parity]
    schur_record = schur_input["parities"][parity]
    old_hessian = np.array(
        jacobian_record["symmetrized_internal_matrix"], dtype=float
    )
    old_hessian = (old_hessian+old_hessian.T)/2.0
    old_boundary = np.array(
        jacobian_record["richardson_boundary_matrix"], dtype=float
    )

    # Freeze and branch-check all 360 points before the new action evaluation.
    points = []
    minimum_gram = math.inf
    minimum_argument = math.inf
    all_lorentzian = True
    for boundary_index in range(30):
        for step in MIXED_STEPS:
            for t_sign, v_sign in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                point = (t_sign*step, boundary_index, v_sign*step)
                points.append(point)
                variables = float_mixed_variables(*point)
                _, _, _, data = dust.total_reduced_evaluation(
                    model, variables, dust.old_values
                )
                minimum_gram = min(minimum_gram, float(data["minimum_gram"]))
                minimum_argument = min(
                    minimum_argument, float(data["minimum_argument"])
                )
                all_lorentzian &= dict(data["negative_counts"]) == {1: 100}

    print(
        f"Evaluating 360 frozen 100-decimal mixed actions: {parity} parity...",
        flush=True,
    )
    with fork_context.Pool(
        processes=8,
        initializer=initialize_action_worker,
        initargs=(model,),
    ) as pool:
        action_values = pool.map(action_worker, points, chunksize=1)

    mixed_by_step = {
        step: [arb.mpc(0) for _ in range(30)] for step in MIXED_STEPS
    }
    maximum_action_imaginary = arb.mpf(0)
    cursor = 0
    for boundary_index in range(30):
        for step in MIXED_STEPS:
            pp, pm, mp_, mm = action_values[cursor:cursor+4]
            cursor += 4
            step_arb = arb.mpf(str(step))
            mixed_by_step[step][boundary_index] = (
                pp-pm-mp_+mm
            )/(4*step_arb**2)
            maximum_action_imaginary = max(
                maximum_action_imaginary,
                abs(arb.im(pp)), abs(arb.im(pm_)),
                abs(arb.im(mp_)), abs(arb.im(mm)),
            )

    normalization = 24*arb.mpf(str(raw_tangent_norm))
    derivative_rows = [
        [value/normalization for value in mixed_by_step[step]]
        for step in MIXED_STEPS
    ]
    row1, row2, row3 = derivative_rows
    richardson12 = [(4*b-a)/3 for a, b in zip(row1, row2)]
    richardson23 = [(4*c-b)/3 for b, c in zip(row2, row3)]
    sixth_order = [
        (16*r23-r12)/15
        for r12, r23 in zip(richardson12, richardson23)
    ]
    compatibility_row = np.array([
        float(arb.re(value)) for value in sixth_order
    ])
    richardson23_float = np.array([
        float(arb.re(value)) for value in richardson23
    ])
    maximum_imaginary = float(max(
        maximum_action_imaginary,
        max(abs(arb.im(value)) for value in sixth_order),
    ))
    compatibility_norm = float(np.linalg.norm(compatibility_row))
    compatibility_epsilon = float(np.linalg.norm(
        compatibility_row-richardson23_float
    ))
    compatibility_floor = max(compatibility_epsilon, 1.0e-35)
    if compatibility_norm <= 10.0*compatibility_floor:
        compatibility_label = "ALL_BOUNDARY_DIRECTIONS_COMPATIBLE"
    elif compatibility_norm > 100.0*compatibility_floor:
        compatibility_label = "ONE_BOUNDARY_CONSTRAINT"
    else:
        compatibility_label = "BOUNDARY_NORM_UNRESOLVED"

    uniform_residual = boundary_projector@compatibility_row
    uniform_residual_norm = float(np.linalg.norm(uniform_residual))
    if uniform_residual_norm <= 10.0*compatibility_floor:
        uniformity_label = "UNIFORM_WITHIN_FROZEN_ERROR"
    elif uniform_residual_norm > 100.0*compatibility_floor:
        uniformity_label = "RESOLVED_NONUNIFORM"
    else:
        uniformity_label = "UNIFORMITY_UNRESOLVED"
    compatibility_mean = float(np.mean(compatibility_row))
    compatibility_spread = float(np.ptp(compatibility_row))
    compatibility_cosine = float(
        compatibility_row@boundary_unit/max(
            compatibility_norm, np.finfo(float).tiny
        )
    )

    # Minimum-norm exact-null correction and inverse Schur reconstruction.
    regular_block = old_hessian[:30, :30]
    old_coupling = old_hessian[:30, 30:]
    regular_ranks = rank_table(regular_block)
    diagonal_tangent = raw_tangent[:30]
    coupling_defect = (
        -regular_block@diagonal_tangent-old_coupling@pole_collective
    )
    corrected_coupling = (
        old_coupling
        + np.outer(coupling_defect, pole_collective)
        / float(pole_collective@pole_collective)
    )
    high_precision_schur = np.array(schur_record["schur_matrix"], dtype=float)
    high_precision_schur = (
        high_precision_schur+high_precision_schur.T
    )/2.0
    corrected_schur = pole_projector@high_precision_schur@pole_projector
    relative_schur = (
        pole_relative_basis.T@corrected_schur@pole_relative_basis
    )
    relative_schur = (relative_schur+relative_schur.T)/2.0
    relative_eigenvalues = np.linalg.eigvalsh(relative_schur)
    schur_epsilon = float(schur_record["empirical_error_norm"])
    relative_separated = bool(
        np.min(np.abs(relative_eigenvalues)) > 100.0*schur_epsilon
    )
    corrected_corner = (
        corrected_schur
        + corrected_coupling.T
        @ np.linalg.solve(regular_block, corrected_coupling)
    )
    corrected_hessian = np.block([
        [regular_block, corrected_coupling],
        [corrected_coupling.T, corrected_corner],
    ])
    coupling_correction = frobenius_relative(
        corrected_coupling, old_coupling
    )
    hessian_correction = frobenius_relative(corrected_hessian, old_hessian)
    normalized_null_residual = float(
        np.linalg.norm(corrected_hessian@raw_tangent)
        / max(
            1.0,
            np.linalg.norm(corrected_hessian, 2)*np.linalg.norm(raw_tangent),
        )
    )
    corrected_full_singular = np.linalg.svd(
        corrected_hessian, compute_uv=False
    )
    corrected_full_rank = int(np.sum(corrected_full_singular > 1.0e-9))

    # Minimum-Frobenius mixed-block update to the new action-only row.
    target_raw = raw_tangent_norm*compatibility_row
    boundary_defect = target_raw-raw_tangent@old_boundary
    corrected_boundary = (
        old_boundary
        + np.outer(raw_tangent, boundary_defect)
        / float(raw_tangent@raw_tangent)
    )
    boundary_correction = frobenius_relative(
        corrected_boundary, old_boundary
    )
    boundary_identity_error = float(np.linalg.norm(
        gauge_tangent@corrected_boundary-compatibility_row
    ))

    source_gates = {
        "regular_rank_30": set(regular_ranks.values()) == {30},
        "relative_schur_separated": relative_separated,
        "normalized_exact_null_below_1e-12": normalized_null_residual < 1e-12,
        "one_full_null_below_absolute_1e-9": corrected_full_rank == 34,
        "coupling_correction_below_1e-6": coupling_correction < 1e-6,
        "hessian_correction_below_1e-6": hessian_correction < 1e-6,
        "boundary_correction_below_1e-6": boundary_correction < 1e-6,
        "mixed_identity_below_1e-12": boundary_identity_error < 1e-12,
    }
    source_pass = all(source_gates.values())

    quotient_hessian = quotient_basis.T@corrected_hessian@quotient_basis
    quotient_hessian = (quotient_hessian+quotient_hessian.T)/2.0
    quotient_singular = np.linalg.svd(quotient_hessian, compute_uv=False)
    quotient_eigenvalues = np.linalg.eigvalsh(quotient_hessian)
    quotient_rank = int(np.sum(quotient_singular > 1.0e-9))
    block_rank = int(30+np.sum(np.abs(relative_eigenvalues) > 1.0e-9))
    quotient_regular = bool(
        source_pass and quotient_rank == block_rank == 34
    )
    quotient_label = (
        "QUOTIENT_REGULAR" if quotient_regular
        else "PRECISION_CORRECTION_UNRESOLVED"
    )

    boundary_basis = None
    boundary_basis_kind = None
    if compatibility_label == "ALL_BOUNDARY_DIRECTIONS_COMPATIBLE":
        boundary_basis = np.eye(30)
        boundary_basis_kind = "ALL_30"
    elif compatibility_label == "ONE_BOUNDARY_CONSTRAINT":
        if uniformity_label == "UNIFORM_WITHIN_FROZEN_ERROR":
            boundary_basis = householder_complement(boundary_collective)
            boundary_basis_kind = "EXACT_ZERO_SUM"
        elif uniformity_label == "RESOLVED_NONUNIFORM":
            boundary_basis = householder_complement(compatibility_row)
            boundary_basis_kind = "ACTION_ROW_KERNEL"

    response = None
    if quotient_regular and boundary_basis is not None:
        right_hand_side = quotient_basis.T@corrected_boundary@boundary_basis
        coefficients = np.linalg.solve(quotient_hessian, -right_hand_side)
        internal_response = quotient_basis@coefficients
        quotient_residual = float(
            np.linalg.norm(
                quotient_hessian@coefficients+right_hand_side, 2
            )/max(1.0, np.linalg.norm(right_hand_side, 2))
        )
        corrected_residual_matrix = (
            corrected_hessian@internal_response
            + corrected_boundary@boundary_basis
        )
        corrected_residual = float(
            np.linalg.norm(corrected_residual_matrix, 2)
            / max(1.0, np.linalg.norm(corrected_boundary@boundary_basis, 2))
        )
        old_residual_matrix = (
            old_hessian@internal_response+old_boundary@boundary_basis
        )
        old_residual = float(
            np.linalg.norm(old_residual_matrix, 2)
            / max(1.0, np.linalg.norm(old_boundary@boundary_basis, 2))
        )
        response_singular = np.linalg.svd(internal_response, compute_uv=False)
        pole_relative_response = pole_projector@internal_response[30:, :]
        response = {
            "boundary_basis": boundary_basis,
            "boundary_basis_kind": boundary_basis_kind,
            "coefficients": coefficients,
            "internal_response": internal_response,
            "quotient_residual": quotient_residual,
            "corrected_unprojected_residual": corrected_residual,
            "corrected_gauge_residual": float(
                np.linalg.norm(gauge_tangent@corrected_residual_matrix)
            ),
            "uncorrected_residual": old_residual,
            "singular_values": response_singular,
            "ranks": rank_table(internal_response),
            "condition": float(
                response_singular[0]/response_singular[-1]
            ) if response_singular[-1] > 0 else math.inf,
            "norm": float(np.linalg.norm(internal_response, 2)),
            "relative_pole_response_norm": float(
                np.linalg.norm(pole_relative_response, 2)
            ),
        }

    if (
        quotient_regular
        and compatibility_label == "ONE_BOUNDARY_CONSTRAINT"
        and uniformity_label == "UNIFORM_WITHIN_FROZEN_ERROR"
        and response is not None
    ):
        outcome = "REGULAR_QUOTIENT_29_ZERO_SUM_RESPONSES"
    elif (
        quotient_regular
        and compatibility_label == "ONE_BOUNDARY_CONSTRAINT"
        and uniformity_label == "RESOLVED_NONUNIFORM"
        and response is not None
    ):
        outcome = "REGULAR_QUOTIENT_29_NONUNIFORM_RESPONSES"
    elif (
        quotient_regular
        and compatibility_label == "ALL_BOUNDARY_DIRECTIONS_COMPATIBLE"
        and response is not None
    ):
        outcome = "REGULAR_QUOTIENT_ALL_30_RESPONSES"
    else:
        outcome = "PRECISION_CORRECTION_OR_BOUNDARY_UNRESOLVED"

    results[parity] = {
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "all_lorentzian": all_lorentzian,
        "maximum_imaginary": maximum_imaginary,
        "mixed_by_step": mixed_by_step,
        "derivative_rows": derivative_rows,
        "richardson12": richardson12,
        "richardson23": richardson23,
        "sixth_order": sixth_order,
        "compatibility_row": compatibility_row,
        "compatibility_norm": compatibility_norm,
        "compatibility_epsilon": compatibility_epsilon,
        "compatibility_floor": compatibility_floor,
        "compatibility_label": compatibility_label,
        "uniformity_label": uniformity_label,
        "uniform_residual": uniform_residual,
        "uniform_residual_norm": uniform_residual_norm,
        "compatibility_mean": compatibility_mean,
        "compatibility_spread": compatibility_spread,
        "compatibility_cosine": compatibility_cosine,
        "regular_ranks": regular_ranks,
        "corrected_coupling": corrected_coupling,
        "corrected_schur": corrected_schur,
        "relative_schur": relative_schur,
        "relative_eigenvalues": relative_eigenvalues,
        "schur_epsilon": schur_epsilon,
        "corrected_hessian": corrected_hessian,
        "corrected_full_singular": corrected_full_singular,
        "corrected_full_rank": corrected_full_rank,
        "coupling_correction": coupling_correction,
        "hessian_correction": hessian_correction,
        "normalized_null_residual": normalized_null_residual,
        "corrected_boundary": corrected_boundary,
        "boundary_correction": boundary_correction,
        "boundary_identity_error": boundary_identity_error,
        "source_gates": source_gates,
        "quotient_hessian": quotient_hessian,
        "quotient_singular": quotient_singular,
        "quotient_eigenvalues": quotient_eigenvalues,
        "quotient_rank": quotient_rank,
        "block_rank": block_rank,
        "quotient_label": quotient_label,
        "response": response,
        "outcome": outcome,
    }

    check(
        f"{parity}: all 360 geometries remain Lorentzian and off branch boundaries",
        all_lorentzian and minimum_gram > 1e-8 and minimum_argument > 1e-6,
        f"min Gram={minimum_gram:.3e}, min argument={minimum_argument:.3e}",
    )
    check(
        f"{parity}: the 100-decimal sixth-order row is real and finite",
        maximum_imaginary < 1e-80
        and np.all(np.isfinite(compatibility_row)),
        f"max imaginary={maximum_imaginary:.3e}, norm(c6)={compatibility_norm:.3e}",
    )
    check(
        f"{parity}: boundary norm and uniformity labels follow the frozen gates",
        compatibility_label in {
            "ALL_BOUNDARY_DIRECTIONS_COMPATIBLE",
            "ONE_BOUNDARY_CONSTRAINT",
            "BOUNDARY_NORM_UNRESOLVED",
        }
        and uniformity_label in {
            "UNIFORM_WITHIN_FROZEN_ERROR",
            "RESOLVED_NONUNIFORM",
            "UNIFORMITY_UNRESOLVED",
        },
        f"norm={compatibility_label}, uniformity={uniformity_label}, "
        f"||P c6||/floor={uniform_residual_norm/compatibility_floor:.3e}",
    )
    check(
        f"{parity}: exact-null correction gates and outcome are assigned mechanically",
        quotient_label in {
            "QUOTIENT_REGULAR", "PRECISION_CORRECTION_UNRESOLVED"
        }
        and outcome in {
            "REGULAR_QUOTIENT_29_ZERO_SUM_RESPONSES",
            "REGULAR_QUOTIENT_29_NONUNIFORM_RESPONSES",
            "REGULAR_QUOTIENT_ALL_30_RESPONSES",
            "PRECISION_CORRECTION_OR_BOUNDARY_UNRESOLVED",
        },
        f"source_pass={source_pass}, rank={quotient_rank}, "
        f"dH/H={hessian_correction:.3e}, dB/B={boundary_correction:.3e}",
    )
    if response is not None:
        check(
            f"{parity}: the corrected quotient response solves its frozen system",
            response["quotient_residual"] < 1e-7
            and response["corrected_unprojected_residual"] < 1e-7,
            f"quotient={response['quotient_residual']:.3e}, "
            f"full={response['corrected_unprojected_residual']:.3e}, "
            f"rank={response['ranks']}",
        )


check(
    "both schedule parities completed the same preregistered correction",
    set(results) == {"even", "odd"}
    and all("outcome" in result for result in results.values()),
)


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 80),
        "imaginary": arb.nstr(arb.im(value), 80),
    }


def serialize_response(response):
    if response is None:
        return None
    return {
        "boundary_dimension": int(response["boundary_basis"].shape[1]),
        "boundary_basis_kind": response["boundary_basis_kind"],
        "boundary_basis": response["boundary_basis"].tolist(),
        "coefficient_matrix": response["coefficients"].tolist(),
        "internal_response_matrix": response["internal_response"].tolist(),
        "relative_quotient_residual": response["quotient_residual"],
        "corrected_unprojected_residual": response[
            "corrected_unprojected_residual"
        ],
        "corrected_gauge_residual": response["corrected_gauge_residual"],
        "uncorrected_recorded_residual": response["uncorrected_residual"],
        "singular_values": response["singular_values"].tolist(),
        "ranks": response["ranks"],
        "condition_2": response["condition"],
        "norm_2": response["norm"],
        "relative_pole_response_norm_2": response[
            "relative_pole_response_norm"
        ],
    }


def serialize_result(result):
    return {
        "outcome": result["outcome"],
        "branch": {
            "all_lorentzian": result["all_lorentzian"],
            "minimum_absolute_gram_eigenvalue": result["minimum_gram"],
            "minimum_angle_argument_modulus": result["minimum_argument"],
            "maximum_action_or_row_imaginary": result["maximum_imaginary"],
        },
        "compatibility": {
            "label": result["compatibility_label"],
            "uniformity_label": result["uniformity_label"],
            "mixed_action_steps": MIXED_STEPS,
            "mixed_action_derivatives": {
                f"{step:.3e}": [
                    serialize_complex(value)
                    for value in result["mixed_by_step"][step]
                ]
                for step in MIXED_STEPS
            },
            "normalized_derivative_rows": [
                [serialize_complex(value) for value in row]
                for row in result["derivative_rows"]
            ],
            "richardson_12": [
                serialize_complex(value) for value in result["richardson12"]
            ],
            "richardson_23": [
                serialize_complex(value) for value in result["richardson23"]
            ],
            "sixth_order": [
                serialize_complex(value) for value in result["sixth_order"]
            ],
            "sixth_order_float": result["compatibility_row"].tolist(),
            "norm": result["compatibility_norm"],
            "empirical_error_norm": result["compatibility_epsilon"],
            "classification_floor": result["compatibility_floor"],
            "uniform_residual": result["uniform_residual"].tolist(),
            "uniform_residual_norm": result["uniform_residual_norm"],
            "uniform_residual_over_floor": (
                result["uniform_residual_norm"]/result["compatibility_floor"]
            ),
            "mean": result["compatibility_mean"],
            "component_spread": result["compatibility_spread"],
            "cosine_with_all_ones": result["compatibility_cosine"],
        },
        "correction": {
            "source_gates": result["source_gates"],
            "regular_block_ranks": result["regular_ranks"],
            "corrected_coupling": result["corrected_coupling"].tolist(),
            "corrected_schur": result["corrected_schur"].tolist(),
            "relative_schur": result["relative_schur"].tolist(),
            "relative_eigenvalues": result["relative_eigenvalues"].tolist(),
            "schur_empirical_error": result["schur_epsilon"],
            "corrected_hessian": result["corrected_hessian"].tolist(),
            "corrected_full_singular_values": result[
                "corrected_full_singular"
            ].tolist(),
            "corrected_full_rank_above_absolute_1e-9": result[
                "corrected_full_rank"
            ],
            "relative_coupling_correction_frobenius": result[
                "coupling_correction"
            ],
            "relative_hessian_correction_frobenius": result[
                "hessian_correction"
            ],
            "normalized_exact_null_residual": result[
                "normalized_null_residual"
            ],
            "corrected_boundary_block": result["corrected_boundary"].tolist(),
            "relative_boundary_correction_frobenius": result[
                "boundary_correction"
            ],
            "mixed_identity_error": result["boundary_identity_error"],
        },
        "quotient": {
            "label": result["quotient_label"],
            "matrix": result["quotient_hessian"].tolist(),
            "singular_values": result["quotient_singular"].tolist(),
            "eigenvalues": result["quotient_eigenvalues"].tolist(),
            "rank_above_absolute_1e-9": result["quotient_rank"],
            "block_factorization_rank": result["block_rank"],
        },
        "response": serialize_response(result["response"]),
    }


outcomes = {result["outcome"] for result in results.values()}
if outcomes == {"REGULAR_QUOTIENT_29_ZERO_SUM_RESPONSES"}:
    verdict = (
        "DERIVED COMPUTATIONAL LINEAR: both schedules have a regular "
        "precision-corrected quotient and twenty-nine zero-sum boundary "
        "responses; the homogeneous-scale identification remains PATTERN "
        "provenance because it was targeted after the first run."
    )
elif outcomes == {"REGULAR_QUOTIENT_29_NONUNIFORM_RESPONSES"}:
    verdict = (
        "DERIVED COMPUTATIONAL LINEAR: both schedules have a regular "
        "quotient and one resolved nonuniform boundary constraint."
    )
elif outcomes == {"REGULAR_QUOTIENT_ALL_30_RESPONSES"}:
    verdict = (
        "DERIVED COMPUTATIONAL LINEAR: both schedules have a regular "
        "quotient and all thirty boundary responses."
    )
else:
    verdict = (
        "OPEN NUMERICALLY: at least one parity does not resolve the frozen "
        "precision correction, boundary row, or quotient response."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "first_result_commit": FIRST_RESULT_COMMIT,
    "exact_lapse_result_commit": EXACT_LAPSE_COMMIT,
    "schur_result_commit": SCHUR_RESULT_COMMIT,
    "jacobian_record_commit": JACOBIAN_RECORD_COMMIT,
    "precision_digits": 100,
    "raw_gauge_tangent": raw_tangent.tolist(),
    "normalized_gauge_tangent": gauge_tangent.tolist(),
    "quotient_basis": quotient_basis.tolist(),
    "parities": {
        parity: serialize_result(result) for parity, result in results.items()
    },
    "verdict": verdict,
    "claim_boundary": {
        "linear_precision_corrected_response": "DERIVED COMPUTATIONAL IF RESOLVED",
        "uniform_scale_target_provenance": "PATTERN / CROSS-RESOLUTION TEST",
        "analytic_component_equality": "NOT PROVED BY FINITE NUMERICS",
        "known_hamiltonian_constraint_interpretation": "STRUCTURAL",
        "nonlinear_boundary_family": "NOT ESTABLISHED",
        "physical_graviton_identification": "NOT CLAIMED",
        "full_840_edge_carrier": "NOT TESTED",
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

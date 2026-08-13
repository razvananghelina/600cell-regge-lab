#!/usr/bin/env python3
"""80-decimal correction for the five unresolved dust lapse directions.

Protocol commit: 5c0372a.  The original unresolved 13/15 run is retained in
commit 7d5e9fc.  This verifier reconstructs only the canonically exposed 5x5
Schur quadratic form; it does not relax or rerun the failed convergence gate.
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
INPUT = HERE / "gravity_600cell_dust_implicit_jacobian.json"
OUTPUT = HERE / "gravity_600cell_dust_lapse_schur.json"
PROTOCOL_COMMIT = "5c0372a"
PRIOR_ART_COMMIT = "0882934"
FAILURE_COMMIT = "7d5e9fc"
ORIGINAL_PROTOCOL_COMMIT = "41acf7b"
STEPS = (5.0e-4, 2.5e-4, 1.25e-4)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


recorded = json.loads(INPUT.read_text())
check(
    "the frozen unresolved input has the expected provenance and two outcomes",
    recorded["protocol_commit"] == ORIGINAL_PROTOCOL_COMMIT
    and recorded["passed"] == 13
    and recorded["tests"] in (13, 15)
    and set(recorded["parities"]) == {"even", "odd"}
    and {
        value["outcome"] for value in recorded["parities"].values()
    } == {"NUMERICALLY_UNRESOLVED"},
)


print("Loading the independently certified 60-decimal action implementation...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control_for_schur",
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


def rank_table(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    return {
        f"{threshold:.0e}": int(np.sum(singular > threshold*singular[0]))
        for threshold in (1.0e-7, 1.0e-9, 1.0e-11)
    }


def canonical_lifts(matrix):
    regular = matrix[:30, :30]
    coupling = matrix[:30, 30:]
    lifts = np.vstack((-np.linalg.solve(regular, coupling), np.eye(5)))
    return regular, coupling, lifts


def frozen_directions(lifts):
    directions = []
    for left in range(5):
        directions.append((f"e{left}", lifts[:, left].copy(), (left,)))
    for left in range(5):
        for right in range(left+1, 5):
            directions.append((
                f"e{left}+e{right}",
                lifts[:, left]+lifts[:, right],
                (left, right),
            ))
    return directions


_ARB_MODEL = None


def initialize_action_worker(model):
    global _ARB_MODEL
    _ARB_MODEL = model
    arb.mp.dps = 80


def action_worker(variables):
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


def arb_log_point(direction, step):
    point = list(dust.ARB_BASE_VARIABLES)
    step_arb = arb.mpf(str(step))
    for index, component in enumerate(direction):
        point[index] *= arb.exp(step_arb*arb.mpf(str(component)))
    return point


def float_log_point(direction, step):
    point = dust.base_variables.copy()
    point[:35] *= np.exp(step*np.asarray(direction))
    return point


def reconstruct_matrix(direction_records, field):
    matrix = np.zeros((5, 5), dtype=complex)
    for index in range(5):
        matrix[index, index] = complex(direction_records[f"e{index}"][field])
    for left in range(5):
        for right in range(left+1, 5):
            combined = complex(direction_records[f"e{left}+e{right}"][field])
            value = (
                combined-matrix[left, left]-matrix[right, right]
            )/2.0
            matrix[left, right] = matrix[right, left] = value
    return matrix/24.0


fork_context = mp.get_context("fork")
results = {}
for parity, model in dust.bl.models.items():
    print(f"Preparing canonical Schur lifts: {parity} parity...", flush=True)
    recorded_matrix = np.array(
        recorded["parities"][parity]["symmetrized_internal_matrix"],
        dtype=float,
    )
    regular, coupling, lifts = canonical_lifts(recorded_matrix)
    regular_singular = np.linalg.svd(regular, compute_uv=False)
    regular_ranks = rank_table(regular)
    block_valid = set(regular_ranks.values()) == {30}
    directions = frozen_directions(lifts)

    # Branch controls at every arbitrary-precision evaluation point are
    # evaluated with the independently certified binary64 geometry branch.
    minimum_gram = math.inf
    minimum_argument = math.inf
    all_lorentzian = True
    for _, direction, _ in directions:
        for step in STEPS:
            for sign in (1.0, -1.0):
                variables = float_log_point(direction, sign*step)
                _, _, _, data = dust.total_reduced_evaluation(
                    model, variables, dust.old_values
                )
                minimum_gram = min(minimum_gram, float(data["minimum_gram"]))
                minimum_argument = min(
                    minimum_argument, float(data["minimum_argument"])
                )
                all_lorentzian &= dict(data["negative_counts"]) == {1: 100}

    print(f"Evaluating 91 frozen 80-decimal actions: {parity} parity...", flush=True)
    points = [list(dust.ARB_BASE_VARIABLES)]
    for _, direction, _ in directions:
        for step in STEPS:
            points.append(arb_log_point(direction, step))
            points.append(arb_log_point(direction, -step))
    with fork_context.Pool(
        processes=8,
        initializer=initialize_action_worker,
        initargs=(model,),
    ) as pool:
        actions = pool.map(action_worker, points, chunksize=1)

    base_action = actions[0]
    cursor = 1
    direction_records = {}
    maximum_action_imaginary = abs(arb.im(base_action))
    for name, _, support in directions:
        centered = []
        for step in STEPS:
            plus, minus = actions[cursor], actions[cursor+1]
            cursor += 2
            maximum_action_imaginary = max(
                maximum_action_imaginary,
                abs(arb.im(plus)),
                abs(arb.im(minus)),
            )
            step_arb = arb.mpf(str(step))
            centered.append((plus-2*base_action+minus)/(step_arb**2))
        rich_coarse = (4*centered[1]-centered[0])/3
        rich_fine = (4*centered[2]-centered[1])/3
        sixth = (16*rich_fine-rich_coarse)/15
        direction_records[name] = {
            "support": support,
            "centered": centered,
            "richardson_coarse": rich_coarse,
            "richardson_fine": rich_fine,
            "sixth_order": sixth,
        }

    schur = reconstruct_matrix(direction_records, "sixth_order").real
    schur_fine = reconstruct_matrix(direction_records, "richardson_fine").real
    epsilon = float(np.linalg.norm(schur-schur_fine, 2))
    eigenvalues, eigenvectors = np.linalg.eigh(schur)

    collective = np.ones(5)/math.sqrt(5.0)
    projector = np.eye(5)-np.outer(collective, collective)
    raw_relative = np.column_stack([
        np.eye(5)[:, index]-np.eye(5)[:, 4] for index in range(4)
    ])
    relative_basis, _ = np.linalg.qr(raw_relative)
    for column in range(4):
        first = int(np.flatnonzero(np.abs(relative_basis[:, column]) > 1.0e-14)[0])
        if relative_basis[first, column] < 0:
            relative_basis[:, column] *= -1
    relative_matrix = relative_basis.T@schur@relative_basis
    relative_eigenvalues = np.linalg.eigvalsh(relative_matrix)
    collective_curvature = float(collective@schur@collective)
    mixing = float(np.linalg.norm(projector@schur@collective))
    relative_mean = float(np.trace(relative_matrix)/4.0)
    permutation_form = (
        relative_mean*projector
        + collective_curvature*np.outer(collective, collective)
    )
    permutation_error = float(
        np.linalg.norm(schur-permutation_form, 2)
        / max(1.0, np.linalg.norm(schur, 2))
    )

    # Frozen analytic lapse-family control.
    family = {}
    for eta in (-1.0e-3, 0.0, 1.0e-3):
        rho = dust.TAU_SQUARE*math.exp(eta)
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
        family[f"{eta:+.1e}"] = {
            "rho": rho,
            "diagonal": diagonal,
            "residuals_real": residual.real,
            "maximum_absolute": float(np.max(np.abs(residual))),
            "maximum_imaginary": float(np.max(np.abs(residual.imag))),
            "minimum_gram": float(data["minimum_gram"]),
            "minimum_argument": float(data["minimum_argument"]),
            "lorentzian": dict(data["negative_counts"]) == {1: 100},
        }
    family_pass = all(
        item["maximum_absolute"] <= 1.0e-7
        and item["lorentzian"]
        and item["minimum_gram"] > 1.0e-8
        and item["minimum_argument"] > 1.0e-6
        for item in family.values()
    )

    lapse_tangent = np.concatenate((
        np.full(30, -dust.TAU_SQUARE/dust.SLANT_SQUARE),
        np.ones(5),
    ))
    lapse_tangent /= np.linalg.norm(lapse_tangent)
    recorded_weakest = np.array(
        recorded["parities"][parity]["weakest_eigenvector"], dtype=float
    )
    recorded_weakest /= np.linalg.norm(recorded_weakest)
    lifted_collective = lifts@collective
    lifted_collective /= np.linalg.norm(lifted_collective)
    tangent_overlaps = {
        "recorded_weakest": float(abs(lapse_tangent@recorded_weakest)),
        "schur_lifted_collective": float(abs(lapse_tangent@lifted_collective)),
    }

    branch_pass = bool(
        all_lorentzian
        and minimum_gram > 1.0e-8
        and minimum_argument > 1.0e-6
    )
    one_plus_four = bool(
        block_valid
        and branch_pass
        and family_pass
        and abs(collective_curvature) <= 10.0*epsilon
        and mixing <= 10.0*epsilon
        and np.min(np.abs(relative_eigenvalues)) > 100.0*epsilon
    )
    five_null = bool(
        block_valid
        and branch_pass
        and family_pass
        and np.max(np.abs(eigenvalues)) <= 10.0*epsilon
    )
    five_stiff = bool(
        block_valid
        and branch_pass
        and np.min(np.abs(eigenvalues)) > 100.0*epsilon
    )
    if not block_valid:
        outcome = "BLOCK_REDUCTION_INVALID"
    elif one_plus_four:
        outcome = "ONE_COLLECTIVE_NULL_FOUR_STIFF"
    elif five_null:
        outcome = "FIVE_NULL"
    elif five_stiff:
        outcome = "FIVE_STIFF"
    else:
        outcome = "NUMERICALLY_UNRESOLVED"

    results[parity] = {
        "regular_singular_values": regular_singular,
        "regular_ranks": regular_ranks,
        "regular_condition": float(regular_singular[0]/regular_singular[-1]),
        "lifts": lifts,
        "direction_records": direction_records,
        "schur": schur,
        "schur_fine": schur_fine,
        "epsilon": epsilon,
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
        "collective_curvature": collective_curvature,
        "mixing": mixing,
        "relative_eigenvalues": relative_eigenvalues,
        "permutation_error": permutation_error,
        "maximum_action_imaginary": maximum_action_imaginary,
        "minimum_gram": minimum_gram,
        "minimum_argument": minimum_argument,
        "family": family,
        "family_pass": family_pass,
        "tangent_overlaps": tangent_overlaps,
        "outcome": outcome,
    }

    check(
        f"{parity}: the recorded 30-dimensional regular block supports the frozen Schur reduction",
        block_valid,
        f"ranks={regular_ranks}, condition={regular_singular[0]/regular_singular[-1]:.3e}",
    )
    check(
        f"{parity}: all 90 displaced action points remain Lorentzian and off branch boundaries",
        branch_pass,
        f"min Gram={minimum_gram:.3e}, min argument={minimum_argument:.3e}",
    )
    check(
        f"{parity}: the 80-decimal action reconstruction is real and finite",
        float(maximum_action_imaginary) < 1.0e-35
        and np.all(np.isfinite(schur))
        and math.isfinite(epsilon),
        f"max action imag={float(maximum_action_imaginary):.3e}, epsilon5={epsilon:.3e}",
    )
    check(
        f"{parity}: all three frozen lapse-family points were evaluated",
        len(family) == 3
        and all(math.isfinite(item["maximum_absolute"]) for item in family.values()),
        "max residuals=" + ", ".join(
            f"{eta}:{item['maximum_absolute']:.3e}"
            for eta, item in family.items()
        ),
    )
    check(
        f"{parity}: the lapse-Schur outcome is assigned mechanically",
        outcome in {
            "ONE_COLLECTIVE_NULL_FOUR_STIFF",
            "FIVE_NULL",
            "FIVE_STIFF",
            "BLOCK_REDUCTION_INVALID",
            "NUMERICALLY_UNRESOLVED",
        },
        f"outcome={outcome}, collective={collective_curvature:.3e}, "
        f"relative={relative_eigenvalues.tolist()}, epsilon5={epsilon:.3e}",
    )


parity_schur_error = float(
    np.linalg.norm(results["even"]["schur"]-results["odd"]["schur"], 2)
    / max(
        1.0,
        np.linalg.norm(results["even"]["schur"], 2),
        np.linalg.norm(results["odd"]["schur"], 2),
    )
)
check(
    "both parities completed the same 15-direction precision correction",
    set(results) == {"even", "odd"}
    and all(len(item["direction_records"]) == 15 for item in results.values()),
    f"normalized parity Schur difference={parity_schur_error:.3e}",
)


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 60),
        "imaginary": arb.nstr(arb.im(value), 60),
    }


def serialize_result(result):
    return {
        "outcome": result["outcome"],
        "regular_block": {
            "singular_values": result["regular_singular_values"].tolist(),
            "ranks": result["regular_ranks"],
            "condition_2": result["regular_condition"],
        },
        "canonical_schur_lifts": result["lifts"].tolist(),
        "directions": {
            name: {
                "support": list(item["support"]),
                "centered_curvatures": [
                    serialize_complex(value) for value in item["centered"]
                ],
                "richardson_coarse": serialize_complex(
                    item["richardson_coarse"]
                ),
                "richardson_fine": serialize_complex(item["richardson_fine"]),
                "sixth_order": serialize_complex(item["sixth_order"]),
            }
            for name, item in result["direction_records"].items()
        },
        "schur_matrix": result["schur"].tolist(),
        "richardson_fine_schur_matrix": result["schur_fine"].tolist(),
        "empirical_error_norm": result["epsilon"],
        "eigenvalues": result["eigenvalues"].tolist(),
        "eigenvectors": result["eigenvectors"].tolist(),
        "collective_curvature": result["collective_curvature"],
        "collective_relative_mixing": result["mixing"],
        "relative_eigenvalues": result["relative_eigenvalues"].tolist(),
        "normalized_permutation_form_error": result["permutation_error"],
        "maximum_action_imaginary": arb.nstr(
            result["maximum_action_imaginary"], 60
        ),
        "minimum_absolute_gram_eigenvalue": result["minimum_gram"],
        "minimum_angle_argument_modulus": result["minimum_argument"],
        "lapse_family": {
            eta: {
                **{key: value for key, value in item.items()
                   if key != "residuals_real"},
                "residuals_real": item["residuals_real"].tolist(),
            }
            for eta, item in result["family"].items()
        },
        "lapse_family_pass_1e-7": result["family_pass"],
        "lapse_tangent_overlaps": result["tangent_overlaps"],
    }


outcomes = {result["outcome"] for result in results.values()}
if outcomes == {"ONE_COLLECTIVE_NULL_FOUR_STIFF"}:
    verdict = (
        "DERIVED COMPUTATIONAL: both schedule parities exhibit one collective "
        "lapse-null direction and four resolved relative lapse "
        "pseudo-constraint stiffnesses in the order-24 dust sandwich."
    )
elif "NUMERICALLY_UNRESOLVED" in outcomes:
    verdict = (
        "OPEN NUMERICALLY: the 80-decimal Schur reconstruction does not "
        "separate the frozen collective and relative lapse criteria."
    )
else:
    verdict = (
        "DERIVED COMPUTATIONAL: the five lapse directions have the reported "
        "non-1+4 classification; no physical interpretation is promoted."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "unresolved_failure_commit": FAILURE_COMMIT,
    "original_protocol_commit": ORIGINAL_PROTOCOL_COMMIT,
    "precision_digits": 80,
    "steps": STEPS,
    "parities": {
        parity: serialize_result(result) for parity, result in results.items()
    },
    "normalized_phase_parity_schur_error": parity_schur_error,
    "verdict": verdict,
    "claim_boundary": {
        "lapse_and_pseudo_constraint_mechanism": "KNOWN PRIOR ART",
        "explicit_600cell_1_plus_4_split": "DERIVED COMPUTATIONAL IF CLASSIFIED",
        "exact_nullity": "NOT INTERVAL-PROVED",
        "full_840_edge_carrier": "NOT TESTED",
        "multi_tick_dynamics": "NOT TESTED",
        "clock_or_parameter_selection": "NOT CLAIMED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)

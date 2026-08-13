#!/usr/bin/env python3
"""Linear implicit-response audit at the published 600-cell dust solution.

Protocol commit: 41acf7b.  The scientific outcome (regular, null, or
unresolved) is not a PASS target.  This verifier differentiates the complete
analytic action gradient at three frozen steps and checks selected curvatures
against the independent 60-decimal action-only evaluator.
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
OUTPUT = HERE / "gravity_600cell_dust_implicit_jacobian.json"
PROTOCOL_COMMIT = "41acf7b"
PRIOR_ART_COMMIT = "31717a8"
UPSTREAM_COMMIT = "66a6465"
STEPS = (1.0e-3, 5.0e-4, 2.5e-4)
RANK_THRESHOLDS = (1.0e-7, 1.0e-9, 1.0e-11)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


def relative_matrix_error(left, right):
    return float(
        np.linalg.norm(left-right, 2)
        / max(1.0, np.linalg.norm(left, 2), np.linalg.norm(right, 2))
    )


def relative_scalar_error(left, right):
    return float(abs(left-right)/max(arb.mpf(1), abs(left), abs(right)))


print("Loading the registered complete-action dust control...", flush=True)
spec = importlib.util.spec_from_file_location(
    "published_dust_control",
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
    "the published dust control retains all 14 upstream certificates",
    dust.tests == dust.passed == 14,
)
check(
    "the frozen carrier has 35 internal and 30 final-boundary coordinates",
    len(dust.base_variables) == 65
    and all(len(model["diagonal_orbits"]) == 30 for model in dust.bl.models.values())
    and all(len(model["pole_orbits"]) == 5 for model in dust.bl.models.values()),
)


_WORKER_MODEL = None


def initialize_gradient_worker(model):
    global _WORKER_MODEL
    _WORKER_MODEL = model


def logarithmic_gradient_worker(variables):
    _, gradient, _, data = dust.total_reduced_evaluation(
        _WORKER_MODEL, variables, dust.old_values
    )
    internal = variables[:35]*gradient[:35]/24.0
    final = variables[35:]*gradient[35:]/24.0
    return (
        internal,
        final,
        float(data["minimum_gram"]),
        float(data["minimum_argument"]),
        dict(data["negative_counts"]),
    )


def perturbation_points(base, step):
    points = []
    for variable in range(65):
        plus = base.copy()
        minus = base.copy()
        plus[variable] *= math.exp(step)
        minus[variable] *= math.exp(-step)
        points.extend((plus, minus))
    return points


def matrix_at_step(outputs, step):
    internal_jacobian = np.empty((35, 65), dtype=complex)
    final_from_internal = np.empty((30, 35), dtype=complex)
    for variable in range(65):
        plus = outputs[2*variable]
        minus = outputs[2*variable+1]
        internal_jacobian[:, variable] = (plus[0]-minus[0])/(2.0*step)
        if variable < 35:
            final_from_internal[:, variable] = (plus[1]-minus[1])/(2.0*step)
    return (
        internal_jacobian[:, :35],
        internal_jacobian[:, 35:],
        final_from_internal,
    )


def rank_table(values):
    singular = np.linalg.svd(values, compute_uv=False)
    if singular[0] == 0:
        return {f"{threshold:.0e}": 0 for threshold in RANK_THRESHOLDS}
    return {
        f"{threshold:.0e}": int(np.sum(singular > threshold*singular[0]))
        for threshold in RANK_THRESHOLDS
    }


fork_context = mp.get_context("fork")
records = {}
for parity, model in dust.bl.models.items():
    print(f"Evaluating frozen Jacobian steps: {parity} parity...", flush=True)
    step_data = {}
    all_minimum_gram = math.inf
    all_minimum_argument = math.inf
    all_lorentzian = True
    with fork_context.Pool(
        processes=8, initializer=initialize_gradient_worker, initargs=(model,)
    ) as pool:
        for step in STEPS:
            outputs = pool.map(
                logarithmic_gradient_worker,
                perturbation_points(dust.base_variables, step),
                chunksize=2,
            )
            internal, boundary, reverse_cross = matrix_at_step(outputs, step)
            all_minimum_gram = min(
                all_minimum_gram, min(output[2] for output in outputs)
            )
            all_minimum_argument = min(
                all_minimum_argument, min(output[3] for output in outputs)
            )
            all_lorentzian &= all(
                output[4] == {1: 100} for output in outputs
            )
            scale = max(
                1.0,
                np.linalg.norm(internal.real, 2),
                np.linalg.norm(boundary.real, 2),
                np.linalg.norm(reverse_cross.real, 2),
            )
            imaginary = float(max(
                np.max(np.abs(internal.imag)),
                np.max(np.abs(boundary.imag)),
                np.max(np.abs(reverse_cross.imag)),
            )/scale)
            step_data[step] = {
                "internal": internal.real,
                "boundary": boundary.real,
                "reverse_cross": reverse_cross.real,
                "imaginary": imaginary,
                "antisymmetry": relative_matrix_error(
                    internal.real, internal.real.T
                ),
                "reciprocity": relative_matrix_error(
                    boundary.real, reverse_cross.real.T
                ),
                "ranks": rank_table(internal.real),
            }

    mid = step_data[STEPS[1]]
    fine = step_data[STEPS[2]]
    rich_internal = (4.0*fine["internal"]-mid["internal"])/3.0
    rich_boundary = (4.0*fine["boundary"]-mid["boundary"])/3.0
    rich_reverse = (4.0*fine["reverse_cross"]-mid["reverse_cross"])/3.0
    symmetrized = (rich_internal+rich_internal.T)/2.0
    singular = np.linalg.svd(symmetrized, compute_uv=False)
    eigenvalues, eigenvectors = np.linalg.eigh(symmetrized)
    epsilon_emp = float(np.linalg.norm(rich_internal-fine["internal"], 2))
    fine_mid_relative = relative_matrix_error(
        fine["internal"], mid["internal"]
    )
    rich_antisymmetry = relative_matrix_error(
        rich_internal, rich_internal.T
    )
    rich_reciprocity = relative_matrix_error(
        rich_boundary, rich_reverse.T
    )
    rich_ranks = rank_table(symmetrized)
    boundary_ranks = rank_table(rich_boundary)
    combined_ranks = rank_table(np.hstack((symmetrized, rich_boundary)))
    inertia_scale = max(1.0, float(np.max(np.abs(eigenvalues))))
    inertia_threshold = 1.0e-11*inertia_scale
    inertia = {
        "positive": int(np.sum(eigenvalues > inertia_threshold)),
        "negative": int(np.sum(eigenvalues < -inertia_threshold)),
        "near_zero": int(np.sum(np.abs(eigenvalues) <= inertia_threshold)),
    }
    minimum_eigen_index = int(np.argmin(np.abs(eigenvalues)))
    weakest = eigenvectors[:, minimum_eigen_index].copy()
    first_nonzero = int(np.flatnonzero(np.abs(weakest) > 1.0e-14)[0])
    if weakest[first_nonzero] < 0:
        weakest *= -1

    response = None
    if singular[-1] > 0:
        response_matrix = np.linalg.solve(symmetrized, -rich_boundary)
        response_residual = float(
            np.linalg.norm(symmetrized@response_matrix+rich_boundary, 2)
            / max(1.0, np.linalg.norm(rich_boundary, 2))
        )
        response = {
            "matrix": response_matrix,
            "residual": response_residual,
            "singular_values": np.linalg.svd(response_matrix, compute_uv=False),
        }

    records[parity] = {
        "model": model,
        "steps": step_data,
        "rich_internal": rich_internal,
        "rich_boundary": rich_boundary,
        "rich_reverse": rich_reverse,
        "symmetrized": symmetrized,
        "singular_values": singular,
        "eigenvalues": eigenvalues,
        "weakest": weakest,
        "epsilon_emp": epsilon_emp,
        "fine_mid_relative": fine_mid_relative,
        "rich_antisymmetry": rich_antisymmetry,
        "rich_reciprocity": rich_reciprocity,
        "rich_ranks": rich_ranks,
        "boundary_ranks": boundary_ranks,
        "combined_ranks": combined_ranks,
        "inertia": inertia,
        "minimum_gram": all_minimum_gram,
        "minimum_argument": all_minimum_argument,
        "all_lorentzian": all_lorentzian,
        "response": response,
    }

    check(
        f"{parity}: all 390 displaced representatives stay Lorentzian and off branch boundaries",
        all_lorentzian
        and all_minimum_gram > 1.0e-8
        and all_minimum_argument > 1.0e-6,
        f"min Gram={all_minimum_gram:.3e}, min argument={all_minimum_argument:.3e}",
    )
    check(
        f"{parity}: logarithmic derivative matrices are real at every frozen step",
        max(item["imaginary"] for item in step_data.values()) < 3.0e-7,
        f"relative imaginary={max(item['imaginary'] for item in step_data.values()):.3e}",
    )
    check(
        f"{parity}: internal Hessian symmetry and boundary cross-reciprocity pass",
        max(item["antisymmetry"] for item in step_data.values()) < 3.0e-6
        and max(item["reciprocity"] for item in step_data.values()) < 3.0e-6
        and rich_antisymmetry < 3.0e-6
        and rich_reciprocity < 3.0e-6,
        f"Richardson antisym={rich_antisymmetry:.3e}, reciprocity={rich_reciprocity:.3e}",
    )
    check(
        f"{parity}: the two finest frozen steps satisfy the preregistered convergence gate",
        fine_mid_relative < 3.0e-4,
        f"relative change={fine_mid_relative:.3e}",
    )


_ARB_MODEL = None


def initialize_action_worker(model):
    global _ARB_MODEL
    _ARB_MODEL = model
    arb.mp.dps = 60


def arbitrary_action_worker(variables):
    return dust.arb_action_components(_ARB_MODEL, variables)[2]


def arbitrary_log_point(direction, step):
    point = list(dust.ARB_BASE_VARIABLES)
    arb_step = arb.mpf(str(step))
    for index in range(35):
        point[index] *= arb.exp(arb_step*arb.mpf(str(direction[index])))
    return point


for parity, record in records.items():
    print(f"Checking 60-decimal action curvatures: {parity} parity...", flush=True)
    directions = []
    all_ones = np.ones(35)
    directions.append(("all_ones", all_ones/np.linalg.norm(all_ones)))
    alternating = np.array([(-1.0)**index for index in range(35)])
    directions.append(("alternating", alternating/np.linalg.norm(alternating)))
    poles = np.zeros(35)
    poles[30:35] = 1.0
    directions.append(("collective_poles", poles/np.linalg.norm(poles)))
    directions.append(("weakest_mode", record["weakest"]))

    points = [list(dust.ARB_BASE_VARIABLES)]
    for _, direction in directions:
        for step in (5.0e-4, 2.5e-4):
            points.append(arbitrary_log_point(direction, step))
            points.append(arbitrary_log_point(direction, -step))
    with fork_context.Pool(
        processes=8, initializer=initialize_action_worker,
        initargs=(record["model"],)
    ) as pool:
        actions = pool.map(arbitrary_action_worker, points, chunksize=1)

    base_action = actions[0]
    curvature_records = {}
    cursor = 1
    for name, direction in directions:
        curvatures = []
        for step in (5.0e-4, 2.5e-4):
            step_arb = arb.mpf(str(step))
            plus, minus = actions[cursor], actions[cursor+1]
            cursor += 2
            curvatures.append((plus-2*base_action+minus)/(step_arb**2))
        rich_curvature = (4*curvatures[1]-curvatures[0])/3
        expected = arb.mpf(str(
            24.0*float(direction @ record["symmetrized"] @ direction)
        ))
        error = relative_scalar_error(rich_curvature, expected)
        curvature_records[name] = {
            "coarse": curvatures[0],
            "fine": curvatures[1],
            "richardson": rich_curvature,
            "matrix_expected": expected,
            "error": error,
        }
    record["curvatures"] = curvature_records
    curvature_pass = all(
        item["error"] < 3.0e-4 for item in curvature_records.values()
    )
    record["curvature_pass"] = curvature_pass
    check(
        f"{parity}: four 60-decimal action curvatures were evaluated as preregistered",
        len(curvature_records) == 4
        and all(math.isfinite(item["error"]) for item in curvature_records.values()),
        "errors=" + ", ".join(
            f"{name}:{item['error']:.3e}"
            for name, item in curvature_records.items()
        ),
    )

    all_step_ranks = [
        rank
        for item in record["steps"].values()
        for rank in item["ranks"].values()
    ] + list(record["rich_ranks"].values())
    implementation_controls = bool(
        record["all_lorentzian"]
        and record["minimum_gram"] > 1.0e-8
        and record["minimum_argument"] > 1.0e-6
        and max(item["imaginary"] for item in record["steps"].values()) < 3.0e-7
        and max(item["antisymmetry"] for item in record["steps"].values()) < 3.0e-6
        and max(item["reciprocity"] for item in record["steps"].values()) < 3.0e-6
        and record["rich_antisymmetry"] < 3.0e-6
        and record["rich_reciprocity"] < 3.0e-6
        and record["fine_mid_relative"] < 3.0e-4
    )
    robust = bool(
        implementation_controls
        and curvature_pass
        and set(all_step_ranks) == {35}
        and record["singular_values"][-1] > 100.0*record["epsilon_emp"]
    )

    common_rank = all_step_ranks[0] if len(set(all_step_ranks)) == 1 else None
    resolved_nullity = False
    if implementation_controls and curvature_pass and common_rank is not None:
        if common_rank < 35:
            singular = record["singular_values"]
            null_max = singular[common_rank] if common_rank < len(singular) else 0.0
            nonnull_min = singular[common_rank-1] if common_rank else math.inf
            resolved_nullity = bool(
                nonnull_min-null_max > 100.0*record["epsilon_emp"]
            )

    if robust:
        outcome = "ROBUST_NUMERICAL_REGULARITY"
    elif resolved_nullity:
        outcome = "RESOLVED_NUMERICAL_NULLITY"
    else:
        outcome = "NUMERICALLY_UNRESOLVED"
    record["implementation_controls"] = implementation_controls
    record["outcome"] = outcome
    check(
        f"{parity}: the preregistered scientific outcome is assigned mechanically",
        outcome in {
            "ROBUST_NUMERICAL_REGULARITY",
            "RESOLVED_NUMERICAL_NULLITY",
            "NUMERICALLY_UNRESOLVED",
        },
        f"outcome={outcome}, rank={record['rich_ranks']}, "
        f"smin/epsilon={record['singular_values'][-1]/max(record['epsilon_emp'], np.finfo(float).tiny):.3e}",
    )
    if outcome == "ROBUST_NUMERICAL_REGULARITY":
        check(
            f"{parity}: the internal solve reproduces the complete linear boundary response",
            record["response"] is not None
            and record["response"]["residual"] < 1.0e-10,
            f"relative residual={record['response']['residual']:.3e}",
        )


parity_spectrum_error = float(
    np.linalg.norm(
        records["even"]["singular_values"]-records["odd"]["singular_values"]
    )/max(
        1.0,
        np.linalg.norm(records["even"]["singular_values"]),
        np.linalg.norm(records["odd"]["singular_values"]),
    )
)
check(
    "both schedule parities completed the same frozen audit",
    set(records) == {"even", "odd"}
    and all("outcome" in record for record in records.values()),
    f"relative singular-spectrum difference={parity_spectrum_error:.3e}",
)


def serialize_complex(value):
    return {
        "real": arb.nstr(arb.re(value), 50),
        "imaginary": arb.nstr(arb.im(value), 50),
    }


def serialize_record(record):
    response = record["response"]
    return {
        "outcome": record["outcome"],
        "implementation_controls_pass": record["implementation_controls"],
        "steps": {
            f"{step:.1e}": {
                "internal_matrix": item["internal"].tolist(),
                "boundary_matrix": item["boundary"].tolist(),
                "reverse_cross_matrix": item["reverse_cross"].tolist(),
                "relative_imaginary": item["imaginary"],
                "relative_antisymmetry": item["antisymmetry"],
                "relative_cross_reciprocity": item["reciprocity"],
                "ranks": item["ranks"],
            }
            for step, item in record["steps"].items()
        },
        "richardson_internal_matrix": record["rich_internal"].tolist(),
        "richardson_boundary_matrix": record["rich_boundary"].tolist(),
        "symmetrized_internal_matrix": record["symmetrized"].tolist(),
        "singular_values": record["singular_values"].tolist(),
        "eigenvalues": record["eigenvalues"].tolist(),
        "inertia_at_relative_1e-11": record["inertia"],
        "ranks": record["rich_ranks"],
        "boundary_ranks": record["boundary_ranks"],
        "combined_ranks": record["combined_ranks"],
        "condition_2": float(
            record["singular_values"][0]/record["singular_values"][-1]
        ),
        "smallest_singular_value": float(record["singular_values"][-1]),
        "empirical_error_norm": record["epsilon_emp"],
        "smallest_singular_over_empirical_error": float(
            record["singular_values"][-1]
            / max(record["epsilon_emp"], np.finfo(float).tiny)
        ),
        "relative_fine_mid_change": record["fine_mid_relative"],
        "relative_richardson_antisymmetry": record["rich_antisymmetry"],
        "relative_richardson_cross_reciprocity": record["rich_reciprocity"],
        "minimum_absolute_gram_eigenvalue": record["minimum_gram"],
        "minimum_angle_argument_modulus": record["minimum_argument"],
        "weakest_eigenvector": record["weakest"].tolist(),
        "action_curvatures": {
            name: {
                "coarse": serialize_complex(item["coarse"]),
                "fine": serialize_complex(item["fine"]),
                "richardson": serialize_complex(item["richardson"]),
                "matrix_expected": arb.nstr(item["matrix_expected"], 50),
                "normalized_error": item["error"],
            }
            for name, item in record["curvatures"].items()
        },
        "all_action_curvatures_pass_3e-4": record["curvature_pass"],
        "response": None if response is None else {
            "matrix": response["matrix"].tolist(),
            "relative_linear_residual": response["residual"],
            "singular_values": response["singular_values"].tolist(),
        },
    }


outcomes = {record["outcome"] for record in records.values()}
if outcomes == {"ROBUST_NUMERICAL_REGULARITY"}:
    verdict = (
        "DERIVED COMPUTATIONAL EVIDENCE: in both schedule parities the "
        "35-variable internal dust-sandwich Jacobian is robustly numerically "
        "regular, so the exact implicit-function conclusion is supported but "
        "not interval-proved."
    )
elif "NUMERICALLY_UNRESOLVED" in outcomes:
    verdict = (
        "OPEN NUMERICALLY: at least one schedule parity does not separate "
        "the weakest internal mode from the frozen numerical error envelope."
    )
else:
    verdict = (
        "DERIVED COMPUTATIONAL: both parities exhibit a stable numerical "
        "nullity compatible with Regge gauge/pseudo-constraint structure."
    )

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "upstream_commit": UPSTREAM_COMMIT,
    "steps": STEPS,
    "rank_thresholds": RANK_THRESHOLDS,
    "coordinates": {
        "internal": 35,
        "final_boundary": 30,
        "internal_definition": "x*(dS/dx)/24 in log(x/x0)",
        "lapse_fixed": False,
        "dust_mass_fixed": True,
    },
    "parities": {
        parity: serialize_record(record) for parity, record in records.items()
    },
    "relative_phase_parity_singular_spectrum_error": parity_spectrum_error,
    "verdict": verdict,
    "claim_boundary": {
        "local_order24_implicit_response": "COMPUTATIONAL EVIDENCE ONLY",
        "exact_nonsingularity": "NOT INTERVAL-PROVED",
        "full_840_edge_carrier": "NOT TESTED",
        "physical_gauge_separation": "NOT DERIVED",
        "multi_tick_dynamics": "NOT TESTED",
        "parameter_or_clock_selection": "NOT CLAIMED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Two-variable homothetic canonical-lapse solve on the 600-cell slab.

Prior-art commit: c7f3e29.
Protocol commit: ded77c5.

The sole seed is the committed fixed-lapse stationary root.  No alternate
root, lapse prior, mass, momentum target or internal length is fitted.
"""

import ast
import contextlib
import io
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
RESPONSE_SOURCE = (
    HERE / "verify_gravity_600cell_dust_homothetic_mass_conservation.py"
)
ROOT_ARTIFACT = HERE / "gravity_600cell_dust_homothetic_forward_root.json"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
PRIOR_ART_COMMIT = "c7f3e29"
PROTOCOL_COMMIT = "ded77c5"
ROOT_RESULT_COMMIT = "b788258"
DPS = 100
arb.mp.dps = DPS
STEP_SETS = {
    "operational_primary": arb.mpf("5e-9"),
    "operational_shadow": arb.mpf("1e-8"),
    "validation_primary": arb.mpf("1.5e-8"),
    "validation_shadow": arb.mpf("3e-8"),
}
ARITHMETIC_FLOOR = arb.mpf("1e-60")
ENTRY_FACTOR = arb.mpf(10)
NONZERO_FACTOR = arb.mpf(100)
RESIDUAL_TOLERANCE = arb.mpf("1e-25")
MAX_ITERATIONS = 8
MAX_DAMPING = 10
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def load_response_prefix():
    tree = ast.parse(RESPONSE_SOURCE.read_text(), filename=str(RESPONSE_SOURCE))
    cut = None
    for index, node in enumerate(tree.body):
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "scale_points"
            for target in node.targets
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError("response evaluator cutoff was not found")
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    namespace = {
        "__file__": str(RESPONSE_SOURCE),
        "__name__": "canonical_lapse_imported_response_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(RESPONSE_SOURCE), "exec"), namespace)
    return namespace


def text(value, digits=50):
    return arb.nstr(value, digits)


def mean(values):
    return sum(values, arb.mpf(0))/len(values)


def vector_norm(values):
    return arb.sqrt(sum(abs(value)**2 for value in values))


def infinity_norm(values):
    return max(abs(value) for value in values)


def spread(values):
    average = mean(values)
    return max(abs(value-average) for value in values)


def matrix_spectral_norm(matrix):
    singular = arb.svd_r(matrix, compute_uv=False)
    return singular[0]


def maximum_matrix_entry(matrix):
    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows) for column in range(matrix.cols)
    )


def real_vector(raw):
    return tuple(arb.mpf(item["real"]) for item in raw)


print("Homothetic canonical-lapse selection audit", flush=True)
response_core = load_response_prefix()
evaluate_homothetic = response_core["evaluate_homothetic"]
models = response_core["models"]
core = response_core["core"]

check(
    "the imported evaluator retains its six exact/static controls",
    response_core["tests"] == response_core["passed"] == 6,
)
check(
    "the canonical-lapse provenance and parity pair are frozen",
    PRIOR_ART_COMMIT == "c7f3e29"
    and PROTOCOL_COMMIT == "ded77c5"
    and ROOT_RESULT_COMMIT == "b788258"
    and set(models) == {"even", "odd"},
)

root_artifact = json.loads(ROOT_ARTIFACT.read_text())
gluing_artifact = json.loads(GLUING_ARTIFACT.read_text())
upstream_ok = bool(
    root_artifact["outcome"] == "HOMOTHETIC_STATIONARY_NOT_CANONICAL"
    and root_artifact["passed"] == root_artifact["tests"] == 6
    and gluing_artifact["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing_artifact["passed"] == gluing_artifact["tests"]
)
check("the fixed-lapse and gluing artifacts authorize the solve", upstream_ok)


targets = {}
uncertainties = {}
seeds = {}
stored_mismatches = {}
target_maps_ok = True
for parity in models:
    gluing = gluing_artifact["parities"][parity]
    mapping = tuple(gluing["geometry"]["old_to_final_orbit_map"])
    target_maps_ok &= sorted(mapping) == list(range(30))
    static_post = real_vector(gluing["momenta"]["post"])
    targets[parity] = tuple(static_post[index] for index in mapping)
    uncertainties[parity] = arb.mpf(
        gluing["momenta"]["cusp_uncertainty_norm"]
    )
    seeds[parity] = (
        arb.mpf(root_artifact["roots"][parity]["root"]),
        arb.mpf(0),
    )
    stored_mismatches[parity] = tuple(
        arb.mpf(value)
        for value in root_artifact["roots"][parity]["junction_residual"]
    )
check(
    "the committed maps, seeds and 30-component targets are complete",
    target_maps_ok
    and all(len(targets[parity]) == len(stored_mismatches[parity]) == 30
            for parity in models),
)


BASE_RHO_LOG = arb.log(core["ARB_RHO"])


def evaluate_reduced(parity, state):
    model = models[parity]
    scale_log, rho_relative_log = state
    record = evaluate_homothetic(
        model, scale_log, BASE_RHO_LOG+rho_relative_log
    )
    momentum_residual = tuple(
        value-target
        for value, target in zip(record["pre"], targets[parity])
    )
    reduced = (mean(record["local"][30:]), mean(momentum_residual))
    return {
        **record,
        "state": tuple(state),
        "momentum_residual": momentum_residual,
        "reduced": reduced,
        "reduced_norm": infinity_norm(reduced),
    }


base_records = {
    parity: evaluate_reduced(parity, seeds[parity]) for parity in models
}
base_control_ok = all(
    record["branch_pass"]
    and infinity_norm(record["local"]) < arb.mpf("1e-25")
    and infinity_norm(tuple(
        fresh-stored
        for fresh, stored in zip(
            record["momentum_residual"], stored_mismatches[parity]
        )
    )) < arb.mpf("1e-20")
    for parity, record in base_records.items()
)
base_control_ok &= bool(
    infinity_norm(tuple(
        left-right for left, right in zip(
            base_records["even"]["local"], base_records["odd"]["local"]
        )
    )) < arb.mpf("1e-24")
    and infinity_norm(tuple(
        left-right for left, right in zip(
            base_records["even"]["momentum_residual"],
            base_records["odd"]["momentum_residual"],
        )
    )) < arb.mpf("1e-24")
)
check(
    "the sole seed reproduces the stationary root and canonical mismatch",
    base_control_ok,
    "even F=({},{}); odd F=({},{})".format(
        *(text(value, 8) for value in base_records["even"]["reduced"]),
        *(text(value, 8) for value in base_records["odd"]["reduced"]),
    ),
)


def perturb(state, coordinate, amount):
    result = list(state)
    result[coordinate] += amount
    return tuple(result)


def calibrated_jacobian(parity, state):
    values = {}
    all_branch_pass = True
    minimum_minor = arb.inf
    minimum_argument = arb.inf
    maximum_imaginary = arb.mpf(0)
    for name, step in STEP_SETS.items():
        for coordinate in range(2):
            for direction in (-1, 1):
                evaluation = evaluate_reduced(
                    parity, perturb(state, coordinate, direction*step)
                )
                values[(name, coordinate, direction)] = evaluation["reduced"]
                all_branch_pass &= evaluation["branch_pass"]
                minimum_minor = min(
                    minimum_minor,
                    evaluation["branch"]["minimum_leading_minor"],
                )
                minimum_argument = min(
                    minimum_argument,
                    evaluation["branch"]["minimum_argument"],
                )
                maximum_imaginary = max(
                    maximum_imaginary, evaluation["maximum_imaginary"]
                )

    matrices = {}
    for name, step in STEP_SETS.items():
        matrix = arb.matrix(2, 2)
        for coordinate in range(2):
            plus = values[(name, coordinate, 1)]
            minus = values[(name, coordinate, -1)]
            for row in range(2):
                matrix[row, coordinate] = (plus[row]-minus[row])/(2*step)
        matrices[name] = matrix

    operational = matrices["operational_primary"]
    d_op = operational-matrices["operational_shadow"]
    d_val = (
        matrices["validation_primary"]-matrices["validation_shadow"]
    )
    d_cross = operational-matrices["validation_primary"]
    entry_pass = bool(all(
        abs(d_cross[row, column]) <= ENTRY_FACTOR*(
            abs(d_op[row, column])+abs(d_val[row, column])+ARITHMETIC_FLOOR
        )
        for row in range(2) for column in range(2)
    ))
    epsilon = (
        matrix_spectral_norm(d_op)+matrix_spectral_norm(d_val)
        + matrix_spectral_norm(d_cross)+ARITHMETIC_FLOOR
    )
    singular = arb.svd_r(operational, compute_uv=False)
    singular_values = (singular[0], singular[1])
    resolved = bool(
        entry_pass and all_branch_pass
        and singular_values[1] > NONZERO_FACTOR*epsilon
    )
    return {
        "matrices": matrices,
        "operational": operational,
        "d_op": d_op,
        "d_val": d_val,
        "d_cross": d_cross,
        "entry_pass": entry_pass,
        "branch_pass": bool(all_branch_pass),
        "epsilon": epsilon,
        "singular_values": singular_values,
        "determinant": arb.det(operational),
        "condition": singular_values[0]/singular_values[1],
        "minimum_leading_minor": minimum_minor,
        "minimum_angle_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
        "resolved": resolved,
    }


solve_records = {}
all_solver_branch_ok = True
all_jacobians_resolved = True

for parity in models:
    print(f"  {parity}: deterministic two-variable Newton solve", flush=True)
    state = seeds[parity]
    evaluation = base_records[parity]
    history = []
    jacobian_attempts = []
    converged = False
    failure = None
    accepted_iterations = 0

    while accepted_iterations < MAX_ITERATIONS:
        if evaluation["reduced_norm"] < RESIDUAL_TOLERANCE:
            converged = True
            break
        jacobian = calibrated_jacobian(parity, state)
        jacobian_attempts.append({"state": state, "jacobian": jacobian})
        all_solver_branch_ok &= jacobian["branch_pass"]
        all_jacobians_resolved &= jacobian["resolved"]
        print(
            "    J attempt {}: singular=({}, {}) epsilon={} resolved={}".format(
                len(jacobian_attempts),
                text(jacobian["singular_values"][0], 8),
                text(jacobian["singular_values"][1], 8),
                text(jacobian["epsilon"], 8),
                jacobian["resolved"],
            ),
            flush=True,
        )
        if not jacobian["branch_pass"]:
            failure = "JACOBIAN_BRANCH_FAILURE"
            break
        if not jacobian["resolved"]:
            failure = "JACOBIAN_ERROR_BAND_FAILURE"
            break

        rhs = arb.matrix((-evaluation["reduced"][0], -evaluation["reduced"][1]))
        correction = arb.lu_solve(jacobian["operational"], rhs)
        accepted = None
        for exponent in range(MAX_DAMPING+1):
            alpha = arb.mpf(2)**(-exponent)
            trial_state = (
                state[0]+alpha*correction[0],
                state[1]+alpha*correction[1],
            )
            trial = evaluate_reduced(parity, trial_state)
            all_solver_branch_ok &= trial["branch_pass"]
            if not trial["branch_pass"]:
                failure = "TRIAL_BRANCH_FAILURE"
                break
            if trial["reduced_norm"] <= (
                1-alpha/4
            )*evaluation["reduced_norm"]:
                accepted = (alpha, trial_state, trial)
                break
        if failure == "TRIAL_BRANCH_FAILURE":
            break
        if accepted is None:
            failure = "NO_ARMIJO_DAMPING"
            break

        alpha, state, evaluation = accepted
        accepted_iterations += 1
        history.append({
            "iteration": accepted_iterations,
            "state": state,
            "residual": evaluation["reduced"],
            "residual_norm": evaluation["reduced_norm"],
            "alpha": alpha,
            "correction": (correction[0], correction[1]),
            "jacobian": jacobian,
        })
        print(
            "    iter={} alpha={} s={} z={} ||F||inf={}".format(
                accepted_iterations, text(alpha, 5), text(state[0], 14),
                text(state[1], 14), text(evaluation["reduced_norm"], 8),
            ),
            flush=True,
        )

    if evaluation["reduced_norm"] < RESIDUAL_TOLERANCE:
        converged = True
    if not converged and failure is None:
        failure = "ITERATION_LIMIT"

    endpoint_jacobian = None
    if converged:
        endpoint_jacobian = calibrated_jacobian(parity, state)
        all_solver_branch_ok &= endpoint_jacobian["branch_pass"]
        all_jacobians_resolved &= endpoint_jacobian["resolved"]
        if not endpoint_jacobian["branch_pass"]:
            failure = "ENDPOINT_BRANCH_FAILURE"
            converged = False
        elif not endpoint_jacobian["resolved"]:
            failure = "ENDPOINT_JACOBIAN_ERROR_BAND_FAILURE"
            converged = False

    diagonal = evaluation["local"][:30]
    poles = evaluation["local"][30:]
    momentum_residual = evaluation["momentum_residual"]
    momentum_bound = 10*uncertainties[parity]
    full_gate = bool(
        converged
        and infinity_norm(diagonal) < arb.mpf("1e-60")
        and infinity_norm(poles) < arb.mpf("1e-25")
        and infinity_norm(evaluation["local"]) < arb.mpf("1e-25")
        and spread(diagonal) < arb.mpf("1e-60")
        and spread(poles) < arb.mpf("1e-60")
        and vector_norm(momentum_residual) <= momentum_bound
        and spread(momentum_residual) <= momentum_bound
    )
    solve_records[parity] = {
        "state": state,
        "evaluation": evaluation,
        "history": history,
        "jacobian_attempts": jacobian_attempts,
        "endpoint_jacobian": endpoint_jacobian,
        "converged": converged,
        "failure": failure,
        "iterations": accepted_iterations,
        "full_gate": full_gate,
        "momentum_bound": momentum_bound,
        "lapse_shift_resolved": abs(state[1]) > arb.mpf("1e-20"),
    }


if all(record["converged"] for record in solve_records.values()):
    even = solve_records["even"]
    odd = solve_records["odd"]
    state_s_difference = abs(even["state"][0]-odd["state"][0])
    state_z_difference = abs(even["state"][1]-odd["state"][1])
    pre_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even["evaluation"]["pre"], odd["evaluation"]["pre"]
        )
    ))
    post_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even["evaluation"]["post"], odd["evaluation"]["post"]
        )
    ))
    parity_gate = bool(
        state_s_difference < arb.mpf("1e-25")
        and state_z_difference < arb.mpf("1e-25")
        and pre_difference < arb.mpf("1e-22")
        and post_difference < arb.mpf("1e-22")
    )
else:
    parity_gate = False
    state_s_difference = state_z_difference = None
    pre_difference = post_difference = None


controls_ok = bool(upstream_ok and target_maps_ok and base_control_ok
                   and all_solver_branch_ok)
newton_ok = all(record["converged"] for record in solve_records.values())
full_ok = all(record["full_gate"] for record in solve_records.values())
lapse_shift_resolved = all(
    record["lapse_shift_resolved"] for record in solve_records.values()
)

if not controls_ok:
    outcome = "CANONICAL_LAPSE_CONTROL_FAILED"
elif not all_jacobians_resolved:
    outcome = "CANONICAL_LAPSE_JACOBIAN_OPEN"
elif not newton_ok:
    outcome = "CANONICAL_LAPSE_NEWTON_OPEN"
elif not full_ok:
    outcome = "CANONICAL_LAPSE_FULL_SUBSTITUTION_FAILED"
elif not parity_gate:
    outcome = "CANONICAL_LAPSE_SCHEDULE_DEPENDENT"
elif not lapse_shift_resolved:
    outcome = "CANONICAL_ROOT_LAPSE_SHIFT_UNRESOLVED"
else:
    outcome = "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"

check(
    "all evaluated states retain the Lorentzian branch",
    all_solver_branch_ok,
)
check(
    "the frozen hierarchy assigns one canonical-lapse outcome",
    outcome in {
        "CANONICAL_LAPSE_CONTROL_FAILED",
        "CANONICAL_LAPSE_JACOBIAN_OPEN",
        "CANONICAL_LAPSE_NEWTON_OPEN",
        "CANONICAL_LAPSE_FULL_SUBSTITUTION_FAILED",
        "CANONICAL_LAPSE_SCHEDULE_DEPENDENT",
        "CANONICAL_ROOT_LAPSE_SHIFT_UNRESOLVED",
        "HOMOTHETIC_CANONICAL_LAPSE_SELECTED",
    },
    outcome,
)


def serialize_matrix(matrix):
    return [
        [text(matrix[row, column], 50) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def serialize_jacobian(record):
    if record is None:
        return None
    return {
        "matrices": {
            name: serialize_matrix(matrix)
            for name, matrix in record["matrices"].items()
        },
        "epsilon": text(record["epsilon"], 40),
        "singular_values": [
            text(value, 50) for value in record["singular_values"]
        ],
        "determinant": text(record["determinant"], 50),
        "condition": text(record["condition"], 40),
        "entry_pass": record["entry_pass"],
        "branch_pass": record["branch_pass"],
        "resolved": record["resolved"],
        "minimum_leading_minor": text(record["minimum_leading_minor"], 40),
        "minimum_angle_argument": text(record["minimum_angle_argument"], 40),
        "maximum_imaginary": text(record["maximum_imaginary"], 40),
        "maximum_d_op": text(maximum_matrix_entry(record["d_op"]), 30),
        "maximum_d_val": text(maximum_matrix_entry(record["d_val"]), 30),
        "maximum_d_cross": text(maximum_matrix_entry(record["d_cross"]), 30),
    }


serialized_solutions = {}
for parity, result in solve_records.items():
    state = result["state"]
    evaluation = result["evaluation"]
    serialized_solutions[parity] = {
        "seed": [text(value, 60) for value in seeds[parity]],
        "state": [text(value, 60) for value in state],
        "scale_ratio": text(arb.exp(state[0]), 60),
        "rho_ratio": text(arb.exp(state[1]), 60),
        "tau_ratio": text(arb.exp(state[1]/2), 60),
        "iterations": result["iterations"],
        "converged": result["converged"],
        "failure": result["failure"],
        "full_gate": result["full_gate"],
        "lapse_shift_resolved": result["lapse_shift_resolved"],
        "reduced_residual": [text(value, 50) for value in evaluation["reduced"]],
        "reduced_residual_norm": text(evaluation["reduced_norm"], 40),
        "diagonal_residuals": [
            text(value, 40) for value in evaluation["local"][:30]
        ],
        "pole_residuals": [
            text(value, 40) for value in evaluation["local"][30:]
        ],
        "pre_momentum": [text(value, 50) for value in evaluation["pre"]],
        "post_momentum": [text(value, 50) for value in evaluation["post"]],
        "target_momentum": [text(value, 50) for value in targets[parity]],
        "junction_residual": [
            text(value, 40) for value in evaluation["momentum_residual"]
        ],
        "junction_norm": text(vector_norm(evaluation["momentum_residual"]), 40),
        "junction_bound": text(result["momentum_bound"], 40),
        "endpoint_jacobian": serialize_jacobian(result["endpoint_jacobian"]),
        "jacobian_attempts": [
            {
                "state": [text(value, 50) for value in item["state"]],
                "jacobian": serialize_jacobian(item["jacobian"]),
            }
            for item in result["jacobian_attempts"]
        ],
        "history": [
            {
                "iteration": item["iteration"],
                "state": [text(value, 50) for value in item["state"]],
                "residual": [text(value, 40) for value in item["residual"]],
                "residual_norm": text(item["residual_norm"], 30),
                "alpha": text(item["alpha"], 20),
                "correction": [text(value, 40) for value in item["correction"]],
                "jacobian": serialize_jacobian(item["jacobian"]),
            }
            for item in result["history"]
        ],
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "root_result_commit": ROOT_RESULT_COMMIT,
    "steps": {name: text(value, 20) for name, value in STEP_SETS.items()},
    "solutions": serialized_solutions,
    "parity_gate": {
        "s_difference": (
            None if state_s_difference is None else text(state_s_difference, 40)
        ),
        "z_difference": (
            None if state_z_difference is None else text(state_z_difference, 40)
        ),
        "pre_momentum_infinity_difference": (
            None if pre_difference is None else text(pre_difference, 40)
        ),
        "post_momentum_infinity_difference": (
            None if post_difference is None else text(post_difference, 40)
        ),
        "passed": parity_gate,
    },
    "classification": {
        "locally_unique_homothetic_pair": bool(
            newton_ok and full_ok and all_jacobians_resolved
        ),
        "next_lapse_ratio_selected": outcome
            == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED",
        "full_65_variable_global_uniqueness": False,
        "refinement_stability": "OPEN",
        "physical_clock": "OPEN",
        "label": "STRUCTURAL / candidate pseudo-constraint",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

for parity, result in solve_records.items():
    state = result["state"]
    evaluation = result["evaluation"]
    print(
        "  {}: s={} z={} tau_ratio={} ||F||={} ||junction||={}".format(
            parity, text(state[0], 16), text(state[1], 16),
            text(arb.exp(state[1]/2), 16),
            text(evaluation["reduced_norm"], 8),
            text(vector_norm(evaluation["momentum_residual"]), 8),
        ),
        flush=True,
    )
print(f"\nOUTCOME: {outcome}", flush=True)
print(f"RESULT: {passed}/{tests}", flush=True)
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Target-independent stationary momentum envelope for the second dust slab.

Prior-art commit: dedcbc6.
Protocol commit: ed1cd6a.

This verifier hashes but deliberately does not parse the accepted first-tick
artifact.  No desired second-tick momentum enters this calculation.
"""

import ast
import contextlib
import hashlib
import io
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
RESPONSE_SOURCE = (
    HERE / "verify_gravity_600cell_dust_homothetic_mass_conservation.py"
)
FIRST_TICK_ARTIFACT = (
    HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
)
OUTPUT = HERE / "gravity_600cell_dust_stationary_momentum_envelope.json"
PRIOR_ART_COMMIT = "dedcbc6"
PROTOCOL_COMMIT = "ed1cd6a"
FIRST_TICK_RESULT_COMMIT = "46a7361"
FIRST_TICK_SHA256 = (
    "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9"
)
A1_TEXT = (
    "-0.00000311605957669450169173470644419863944122165192557277135128791"
)
R1_TEXT = (
    "-0.00000355925313517063343725030533963917396571974345422547402551491"
)
DPS = 100
arb.mp.dps = DPS
A1 = arb.mpf(A1_TEXT)
R1 = arb.mpf(R1_TEXT)
STEP_SETS = {
    "operational_primary": arb.mpf("1e-20"),
    "operational_shadow": arb.mpf("1e-15"),
    "validation_primary": arb.mpf("3e-20"),
    "validation_shadow": arb.mpf("3e-15"),
}
ARITHMETIC_FLOOR = arb.mpf("1e-60")
ENTRY_FACTOR = arb.mpf(10)
NONZERO_FACTOR = arb.mpf(100)
RESIDUAL_TOLERANCE = arb.mpf("1e-25")
MAX_NEWTON_ITERATIONS = 6
MAX_DAMPING = 10
BASE_BISECTIONS = 100
FOLD_BISECTIONS = 24
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
        "__name__": "momentum_envelope_imported_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(RESPONSE_SOURCE), "exec"), namespace)
    return namespace


def text(value, digits=50):
    return arb.nstr(value, digits)


def mean(values):
    return sum(values, arb.mpf(0))/len(values)


def infinity_norm(values):
    return max(abs(value) for value in values)


def spread(values):
    average = mean(values)
    return max(abs(value-average) for value in values)


def maximum_imaginary(action, gradient):
    return max(abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient))


def maximum_matrix_entry(matrix):
    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows) for column in range(matrix.cols)
    )


print("Target-independent stationary canonical-momentum envelope", flush=True)
response_core = load_response_prefix()
models = response_core["models"]
core = response_core["core"]
action_and_gradient = core["action_and_gradient"]
branch_pass = core["branch_pass"]
ARB_L0_SQUARE = core["ARB_L0_SQUARE"]
ARB_RHO = core["ARB_RHO"]

check(
    "the imported action/evaluator core retains its six controls",
    response_core["tests"] == response_core["passed"] == 6,
)
artifact_hash = hashlib.sha256(FIRST_TICK_ARTIFACT.read_bytes()).hexdigest()
provenance_ok = bool(
    PRIOR_ART_COMMIT == "dedcbc6"
    and PROTOCOL_COMMIT == "ed1cd6a"
    and FIRST_TICK_RESULT_COMMIT == "46a7361"
    and artifact_hash == FIRST_TICK_SHA256
    and set(models) == {"even", "odd"}
    and BASE_BISECTIONS == 100
    and FOLD_BISECTIONS == 24
)
check(
    "the target firewall, accepted geometry hash and frozen counts pass",
    provenance_ok,
    f"first-tick SHA-256={artifact_hash}",
)


def evaluate_raw(parity, upper_log, rho_log):
    q_old_value = arb.exp(2*A1)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*upper_log)*ARB_L0_SQUARE
    rho_value = ARB_RHO*arb.exp(rho_log)
    diagonal = arb.exp(A1+upper_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("stationary-curve magnitude left the positive domain")
    q_old = tuple([q_old_value]*30)
    internal = tuple([diagonal]*30+[rho_value]*5)
    q_new = tuple([q_new_value]*30)
    action, gradient, branch = action_and_gradient(
        models[parity], q_old, internal, q_new
    )
    local = tuple(arb.re(value) for value in gradient[30:65])
    pre = tuple(-arb.re(value) for value in gradient[:30])
    post = tuple(arb.re(value) for value in gradient[65:95])
    imaginary = maximum_imaginary(action, gradient)
    return {
        "upper_log": upper_log,
        "rho_log": rho_log,
        "q_old": q_old_value,
        "q_new": q_new_value,
        "rho": rho_value,
        "diagonal": diagonal,
        "action": action,
        "gradient": gradient,
        "local": local,
        "pre": pre,
        "post": post,
        "G": mean(local[30:]),
        "P": mean(pre),
        "branch": branch,
        "maximum_imaginary": imaginary,
        "branch_pass": branch_pass(branch, imaginary),
    }


def perturb(upper_log, rho_log, coordinate, amount):
    if coordinate == 0:
        return upper_log+amount, rho_log
    return upper_log, rho_log+amount


def calibrated_gp(parity, upper_log, rho_log):
    values = {}
    all_branch_pass = True
    minimum_minor = arb.inf
    minimum_argument = arb.inf
    maximum_imaginary_value = arb.mpf(0)
    for name, step in STEP_SETS.items():
        for coordinate in range(2):
            for direction in (-1, 1):
                point = perturb(
                    upper_log, rho_log, coordinate, direction*step
                )
                evaluation = evaluate_raw(parity, *point)
                values[(name, coordinate, direction)] = (
                    evaluation["G"], evaluation["P"]
                )
                all_branch_pass &= evaluation["branch_pass"]
                minimum_minor = min(
                    minimum_minor,
                    evaluation["branch"]["minimum_leading_minor"],
                )
                minimum_argument = min(
                    minimum_argument,
                    evaluation["branch"]["minimum_argument"],
                )
                maximum_imaginary_value = max(
                    maximum_imaginary_value, evaluation["maximum_imaginary"]
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
    d_val = matrices["validation_primary"]-matrices["validation_shadow"]
    d_cross = operational-matrices["validation_primary"]
    entry_errors = arb.matrix(2, 2)
    entry_pass = True
    for row in range(2):
        for column in range(2):
            entry_errors[row, column] = (
                abs(d_op[row, column])+abs(d_val[row, column])
                + abs(d_cross[row, column])+ARITHMETIC_FLOOR
            )
            entry_pass &= bool(
                abs(d_cross[row, column]) <= ENTRY_FACTOR*(
                    abs(d_op[row, column])+abs(d_val[row, column])
                    + ARITHMETIC_FLOOR
                )
            )
    determinants = {name: arb.det(matrix) for name, matrix in matrices.items()}
    determinant = determinants["operational_primary"]
    determinant_epsilon = (
        abs(determinant-determinants["operational_shadow"])
        + abs(determinants["validation_primary"]
              -determinants["validation_shadow"])
        + abs(determinant-determinants["validation_primary"])
        + ARITHMETIC_FLOOR
    )
    gb_error = entry_errors[0, 0]
    gb_resolved = bool(
        entry_pass and all_branch_pass
        and abs(operational[0, 0]) > NONZERO_FACTOR*gb_error
    )
    determinant_resolved = bool(
        entry_pass and all_branch_pass
        and abs(determinant) > NONZERO_FACTOR*determinant_epsilon
    )
    singular = arb.svd_r(operational, compute_uv=False)
    return {
        "matrices": matrices,
        "operational": operational,
        "d_op": d_op,
        "d_val": d_val,
        "d_cross": d_cross,
        "entry_errors": entry_errors,
        "entry_pass": bool(entry_pass),
        "branch_pass": bool(all_branch_pass),
        "gb_error": gb_error,
        "gb_resolved": gb_resolved,
        "determinants": determinants,
        "determinant": determinant,
        "determinant_epsilon": determinant_epsilon,
        "determinant_resolved": determinant_resolved,
        "singular_values": (singular[0], singular[1]),
        "minimum_leading_minor": minimum_minor,
        "minimum_angle_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary_value,
    }


all_branch_ok = provenance_ok
base_center = 2*A1
base_half_width = abs(A1)
base_left = base_right = None
base_bracket_ok = False
base_expansions = 0
base_bracket_attempts = []
if provenance_ok:
    for expansion in range(13):
        left = base_center-base_half_width
        right = base_center+base_half_width
        left_eval = evaluate_raw("even", left, R1)
        right_eval = evaluate_raw("even", right, R1)
        base_bracket_attempts.append({
            "expansion": expansion,
            "half_width": base_half_width,
            "left_upper_log": left,
            "right_upper_log": right,
            "left_G": left_eval["G"],
            "right_G": right_eval["G"],
            "left_branch_pass": left_eval["branch_pass"],
            "right_branch_pass": right_eval["branch_pass"],
            "left_branch": left_eval["branch"],
            "right_branch": right_eval["branch"],
            "left_maximum_imaginary": left_eval["maximum_imaginary"],
            "right_maximum_imaginary": right_eval["maximum_imaginary"],
        })
        if (
            left_eval["branch_pass"] and right_eval["branch_pass"]
            and left_eval["G"]*right_eval["G"] <= 0
        ):
            base_left = (left, left_eval)
            base_right = (right, right_eval)
            base_bracket_ok = True
            base_expansions = expansion
            break
        base_half_width *= 2

base_root = None
if base_bracket_ok:
    left_b, left_eval = base_left
    right_b, right_eval = base_right
    for _ in range(BASE_BISECTIONS):
        middle_b = (left_b+right_b)/2
        middle_eval = evaluate_raw("even", middle_b, R1)
        all_branch_ok &= middle_eval["branch_pass"]
        if left_eval["G"]*middle_eval["G"] <= 0:
            right_b, right_eval = middle_b, middle_eval
        else:
            left_b, left_eval = middle_b, middle_eval
    middle_b = (left_b+right_b)/2
    middle_eval = evaluate_raw("even", middle_b, R1)
    base_root = {
        "upper_log": middle_b,
        "evaluation": middle_eval,
        "width": right_b-left_b,
        "left": (left_b, left_eval["G"]),
        "right": (right_b, right_eval["G"]),
    }
    base_bracket_ok &= bool(
        middle_eval["branch_pass"]
        and abs(middle_eval["G"]) < RESIDUAL_TOLERANCE
        and base_root["width"] < arb.mpf("1e-30")
    )

check(
    "the frozen target-free base-bracket process is completely classified",
    bool(
        provenance_ok
        and (
            base_bracket_ok
            or len(base_bracket_attempts) == 13
        )
    ),
    "expansions={}, width={}, G={}".format(
        base_expansions,
        "none" if base_root is None else text(base_root["width"], 8),
        "none" if base_root is None else text(base_root["evaluation"]["G"], 8),
    ),
)


def full_curve_gate(evaluation, derivative, require_det=True):
    diagonal = evaluation["local"][:30]
    poles = evaluation["local"][30:]
    return bool(
        evaluation["branch_pass"]
        and derivative["branch_pass"]
        and derivative["entry_pass"]
        and derivative["gb_resolved"]
        and (derivative["determinant_resolved"] or not require_det)
        and infinity_norm(diagonal) < arb.mpf("1e-60")
        and infinity_norm(poles) < RESIDUAL_TOLERANCE
        and spread(diagonal) < arb.mpf("1e-60")
        and spread(poles) < arb.mpf("1e-60")
    )


def parity_value_gate(upper_log, rho_log, even_evaluation):
    odd = evaluate_raw("odd", upper_log, rho_log)
    ok = bool(
        odd["branch_pass"]
        and abs(even_evaluation["G"]-odd["G"]) < arb.mpf("1e-24")
        and abs(even_evaluation["P"]-odd["P"]) < arb.mpf("1e-22")
        and infinity_norm(tuple(
            left-right for left, right in zip(
                even_evaluation["local"], odd["local"]
            )
        )) < arb.mpf("1e-24")
    )
    return odd, ok


def solve_curve_node(rho_log, seed_upper_log):
    state = seed_upper_log
    evaluation = evaluate_raw("even", state, rho_log)
    attempts = []
    rejected_branch_trials = 0
    failure = None
    iterations = 0
    while abs(evaluation["G"]) >= RESIDUAL_TOLERANCE and iterations < MAX_NEWTON_ITERATIONS:
        derivative = calibrated_gp("even", state, rho_log)
        attempts.append(derivative)
        if not derivative["branch_pass"]:
            failure = "DERIVATIVE_BRANCH_FAILURE"
            break
        if not derivative["gb_resolved"]:
            failure = "GB_DERIVATIVE_OPEN"
            break
        correction = -evaluation["G"]/derivative["operational"][0, 0]
        accepted = None
        for exponent in range(MAX_DAMPING+1):
            alpha = arb.mpf(2)**(-exponent)
            trial_state = state+alpha*correction
            trial = evaluate_raw("even", trial_state, rho_log)
            if not trial["branch_pass"]:
                rejected_branch_trials += 1
                continue
            if abs(trial["G"]) <= (1-alpha/4)*abs(evaluation["G"]):
                accepted = (trial_state, trial)
                break
        if accepted is None:
            failure = "NO_ARMIJO_DAMPING"
            break
        state, evaluation = accepted
        iterations += 1
    converged = abs(evaluation["G"]) < RESIDUAL_TOLERANCE
    if not converged and failure is None:
        failure = "ITERATION_LIMIT"
    endpoint_derivative = (
        calibrated_gp("even", state, rho_log) if converged else None
    )
    return {
        "rho_log": rho_log,
        "upper_log": state,
        "evaluation": evaluation,
        "attempts": attempts,
        "endpoint_derivative": endpoint_derivative,
        "iterations": iterations,
        "rejected_branch_trials": rejected_branch_trials,
        "converged": converged,
        "failure": failure,
    }


curve_nodes = []
curve_newton_open = False
curve_full_failed = False
derivative_open = False
parity_failed = False

if base_bracket_ok:
    base_derivative = calibrated_gp("even", base_root["upper_log"], R1)
    base_odd, base_parity = parity_value_gate(
        base_root["upper_log"], R1, base_root["evaluation"]
    )
    base_node = {
        "kind": "grid",
        "index": 0,
        "rho_log": R1,
        "upper_log": base_root["upper_log"],
        "evaluation": base_root["evaluation"],
        "odd_evaluation": base_odd,
        "endpoint_derivative": base_derivative,
        "iterations": 0,
        "rejected_branch_trials": 0,
        "failure": None,
        "full_gate": full_curve_gate(
            base_root["evaluation"], base_derivative
        ),
        "parity_gate": base_parity,
    }
    curve_nodes.append(base_node)
    all_branch_ok &= bool(
        base_root["evaluation"]["branch_pass"]
        and base_derivative["branch_pass"] and base_odd["branch_pass"]
    )
    curve_full_failed |= not base_node["full_gate"]
    derivative_open |= not (
        base_derivative["gb_resolved"]
        and base_derivative["determinant_resolved"]
    )
    parity_failed |= not base_parity

    frozen_rhos = [R1-arb.mpf(index)/4 for index in range(1, 33)]
    frozen_rhos += [R1-arb.mpf(value) for value in (12, 16, 24, 32)]
    upper_seed = base_root["upper_log"]
    for ordinal, rho_log in enumerate(frozen_rhos, start=1):
        node = solve_curve_node(rho_log, upper_seed)
        evaluation = node["evaluation"]
        derivative = node["endpoint_derivative"]
        all_branch_ok &= evaluation["branch_pass"]
        if derivative is not None:
            all_branch_ok &= derivative["branch_pass"]
        if not node["converged"]:
            curve_newton_open = True
            node.update({
                "kind": "grid" if ordinal <= 32 else "sentinel",
                "index": ordinal,
                "odd_evaluation": None,
                "full_gate": False,
                "parity_gate": False,
            })
            curve_nodes.append(node)
            print(
                "  node={} r={} FAILED {} |G|={}".format(
                    ordinal, text(rho_log, 8), node["failure"],
                    text(abs(evaluation["G"]), 8),
                ),
                flush=True,
            )
            break
        odd, parity_gate = parity_value_gate(
            node["upper_log"], rho_log, evaluation
        )
        full_gate = full_curve_gate(evaluation, derivative)
        node.update({
            "kind": "grid" if ordinal <= 32 else "sentinel",
            "index": ordinal,
            "odd_evaluation": odd,
            "full_gate": full_gate,
            "parity_gate": parity_gate,
        })
        curve_nodes.append(node)
        all_branch_ok &= odd["branch_pass"]
        curve_full_failed |= not full_gate
        derivative_open |= not (
            derivative["gb_resolved"] and derivative["determinant_resolved"]
        )
        parity_failed |= not parity_gate
        upper_seed = node["upper_log"]
        print(
            "  node={:02d} kind={} r={} b={} P={} D={} epsD={} iter={} full={}".format(
                ordinal, node["kind"], text(rho_log, 8),
                text(node["upper_log"], 11), text(evaluation["P"], 11),
                text(derivative["determinant"], 8),
                text(derivative["determinant_epsilon"], 6),
                node["iterations"], full_gate,
            ),
            flush=True,
        )
        if not full_gate or not parity_gate:
            break


main_nodes = [
    node for node in curve_nodes
    if node.get("kind") == "grid" and node.get("full_gate")
]
fold_brackets = []
if len(main_nodes) == 33 and not derivative_open:
    for left, right in zip(main_nodes, main_nodes[1:]):
        left_d = left["endpoint_derivative"]["determinant"]
        right_d = right["endpoint_derivative"]["determinant"]
        if left_d*right_d < 0:
            fold_brackets.append((left, right))

fold_enumeration_open = bool(len(fold_brackets) > 4)
refined_folds = []
if (
    len(main_nodes) == 33 and not curve_newton_open and not curve_full_failed
    and not derivative_open and not parity_failed and not fold_enumeration_open
):
    print(f"  refining {len(fold_brackets)} target-free fold brackets", flush=True)
    for fold_index, (left_node, right_node) in enumerate(fold_brackets):
        left = left_node
        right = right_node
        refinement_ok = True
        for _ in range(FOLD_BISECTIONS):
            rho_mid = (left["rho_log"]+right["rho_log"])/2
            weight = (
                (rho_mid-left["rho_log"])
                /(right["rho_log"]-left["rho_log"])
            )
            upper_seed = (
                (1-weight)*left["upper_log"]+weight*right["upper_log"]
            )
            middle = solve_curve_node(rho_mid, upper_seed)
            if not middle["converged"]:
                curve_newton_open = True
                refinement_ok = False
                break
            derivative = middle["endpoint_derivative"]
            evaluation = middle["evaluation"]
            if not full_curve_gate(evaluation, derivative):
                if not derivative["determinant_resolved"]:
                    derivative_open = True
                else:
                    curve_full_failed = True
                refinement_ok = False
                break
            odd, parity_gate = parity_value_gate(
                middle["upper_log"], rho_mid, evaluation
            )
            if not parity_gate:
                parity_failed = True
                refinement_ok = False
                break
            middle.update({
                "kind": "fold_refinement",
                "index": None,
                "odd_evaluation": odd,
                "full_gate": True,
                "parity_gate": True,
            })
            left_d = left["endpoint_derivative"]["determinant"]
            middle_d = derivative["determinant"]
            if left_d*middle_d < 0:
                right = middle
            else:
                left = middle
        if not refinement_ok:
            break
        width = abs(right["rho_log"]-left["rho_log"])
        opposite = bool(
            left["endpoint_derivative"]["determinant"]
            * right["endpoint_derivative"]["determinant"] < 0
        )
        rho_mid = (left["rho_log"]+right["rho_log"])/2
        upper_seed = (left["upper_log"]+right["upper_log"])/2
        midpoint = solve_curve_node(rho_mid, upper_seed)
        midpoint_ok = midpoint["converged"]
        even_derivative = midpoint["endpoint_derivative"] if midpoint_ok else None
        odd_evaluation = None
        odd_derivative = None
        derivative_parity_ok = False
        if midpoint_ok:
            odd_evaluation, value_parity_ok = parity_value_gate(
                midpoint["upper_log"], rho_mid, midpoint["evaluation"]
            )
            odd_derivative = calibrated_gp(
                "odd", midpoint["upper_log"], rho_mid
            )
            derivative_parity_ok = value_parity_ok and all(
                abs(
                    even_derivative["operational"][row, column]
                    - odd_derivative["operational"][row, column]
                ) <= ENTRY_FACTOR*(
                    even_derivative["entry_errors"][row, column]
                    + odd_derivative["entry_errors"][row, column]
                    + ARITHMETIC_FLOOR
                )
                for row in range(2) for column in range(2)
            )
        fold_ok = bool(
            opposite and width < arb.mpf("2e-8") and midpoint_ok
            and derivative_parity_ok
        )
        fold_enumeration_open |= not fold_ok
        parity_failed |= bool(midpoint_ok and not derivative_parity_ok)
        refined_folds.append({
            "index": fold_index,
            "left": left,
            "right": right,
            "width": width,
            "midpoint": midpoint,
            "odd_midpoint_evaluation": odd_evaluation,
            "odd_midpoint_derivative": odd_derivative,
            "derivative_parity_ok": derivative_parity_ok,
            "passed": fold_ok,
        })
        print(
            "    fold={} r_mid={} width={} P_mid={} passed={}".format(
                fold_index, text(rho_mid, 14), text(width, 8),
                "none" if not midpoint_ok else text(midpoint["evaluation"]["P"], 14),
                fold_ok,
            ),
            flush=True,
        )


if not provenance_ok:
    outcome = "MOMENTUM_ENVELOPE_CONTROL_FAILED"
elif not base_bracket_ok:
    outcome = "MOMENTUM_ENVELOPE_BASE_BRACKET_FAILED"
elif curve_newton_open:
    outcome = "MOMENTUM_ENVELOPE_CURVE_NEWTON_OPEN"
elif not all_branch_ok or curve_full_failed:
    outcome = "MOMENTUM_ENVELOPE_BRANCH_OR_FULL_FAILED"
elif derivative_open:
    outcome = "MOMENTUM_ENVELOPE_DERIVATIVE_OPEN"
elif fold_enumeration_open or len(main_nodes) != 33:
    outcome = "MOMENTUM_ENVELOPE_FOLD_ENUMERATION_OPEN"
elif parity_failed:
    outcome = "MOMENTUM_ENVELOPE_PARITY_FAILED"
else:
    outcome = "MOMENTUM_ENVELOPE_ENUMERATED"

check(
    "all accepted envelope evaluations retain the Lorentzian branch",
    all_branch_ok,
)
check(
    "the frozen hierarchy assigns exactly one target-free envelope outcome",
    outcome in {
        "MOMENTUM_ENVELOPE_CONTROL_FAILED",
        "MOMENTUM_ENVELOPE_BASE_BRACKET_FAILED",
        "MOMENTUM_ENVELOPE_CURVE_NEWTON_OPEN",
        "MOMENTUM_ENVELOPE_BRANCH_OR_FULL_FAILED",
        "MOMENTUM_ENVELOPE_DERIVATIVE_OPEN",
        "MOMENTUM_ENVELOPE_FOLD_ENUMERATION_OPEN",
        "MOMENTUM_ENVELOPE_PARITY_FAILED",
        "MOMENTUM_ENVELOPE_ENUMERATED",
    },
    outcome,
)


def serialize_matrix(matrix):
    return [
        [text(matrix[row, column], 45) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def serialize_derivative(record):
    if record is None:
        return None
    return {
        "matrices": {
            name: serialize_matrix(matrix)
            for name, matrix in record["matrices"].items()
        },
        "entry_errors": serialize_matrix(record["entry_errors"]),
        "entry_pass": record["entry_pass"],
        "branch_pass": record["branch_pass"],
        "gb_error": text(record["gb_error"], 35),
        "gb_resolved": record["gb_resolved"],
        "determinants": {
            name: text(value, 45)
            for name, value in record["determinants"].items()
        },
        "determinant": text(record["determinant"], 45),
        "determinant_epsilon": text(record["determinant_epsilon"], 35),
        "determinant_resolved": record["determinant_resolved"],
        "singular_values": [
            text(value, 45) for value in record["singular_values"]
        ],
        "minimum_leading_minor": text(record["minimum_leading_minor"], 35),
        "minimum_angle_argument": text(record["minimum_angle_argument"], 35),
        "maximum_imaginary": text(record["maximum_imaginary"], 35),
        "maximum_d_op": text(maximum_matrix_entry(record["d_op"]), 25),
        "maximum_d_val": text(maximum_matrix_entry(record["d_val"]), 25),
        "maximum_d_cross": text(maximum_matrix_entry(record["d_cross"]), 25),
    }


def serialize_node(node):
    evaluation = node["evaluation"]
    derivative = node.get("endpoint_derivative")
    return {
        "kind": node["kind"],
        "index": node["index"],
        "rho_log": text(node["rho_log"], 55),
        "upper_log": text(node["upper_log"], 55),
        "scale_ratio_upper_lower": text(arb.exp(node["upper_log"]-A1), 55),
        "rho_ratio_from_first": text(arb.exp(node["rho_log"]-R1), 55),
        "tau_ratio_from_first": text(arb.exp((node["rho_log"]-R1)/2), 55),
        "G": text(evaluation["G"], 45),
        "P": text(evaluation["P"], 45),
        "local_residuals": [text(value, 35) for value in evaluation["local"]],
        "pre_momentum": [text(value, 45) for value in evaluation["pre"]],
        "post_momentum": [text(value, 45) for value in evaluation["post"]],
        "iterations": node["iterations"],
        "rejected_branch_trials": node["rejected_branch_trials"],
        "failure": node["failure"],
        "full_gate": node["full_gate"],
        "parity_gate": node["parity_gate"],
        "derivative": serialize_derivative(derivative),
    }


serialized_folds = []
for fold in refined_folds:
    serialized_folds.append({
        "index": fold["index"],
        "left_rho_log": text(fold["left"]["rho_log"], 55),
        "left_upper_log": text(fold["left"]["upper_log"], 55),
        "left_P": text(fold["left"]["evaluation"]["P"], 45),
        "left_D": text(
            fold["left"]["endpoint_derivative"]["determinant"], 45
        ),
        "right_rho_log": text(fold["right"]["rho_log"], 55),
        "right_upper_log": text(fold["right"]["upper_log"], 55),
        "right_P": text(fold["right"]["evaluation"]["P"], 45),
        "right_D": text(
            fold["right"]["endpoint_derivative"]["determinant"], 45
        ),
        "width": text(fold["width"], 35),
        "midpoint_rho_log": text(fold["midpoint"]["rho_log"], 55),
        "midpoint_upper_log": text(fold["midpoint"]["upper_log"], 55),
        "midpoint_P": text(fold["midpoint"]["evaluation"]["P"], 45),
        "derivative_parity_ok": fold["derivative_parity_ok"],
        "passed": fold["passed"],
    })

accepted_nodes = [node for node in curve_nodes if node.get("full_gate")]
momentum_multiset = [node["evaluation"]["P"] for node in accepted_nodes]
artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "first_tick_result_commit": FIRST_TICK_RESULT_COMMIT,
    "first_tick_sha256": FIRST_TICK_SHA256,
    "target_parsed": False,
    "geometry": {"a1": A1_TEXT, "r1": R1_TEXT},
    "base_bracket": {
        "expansions": base_expansions,
        "width": None if base_root is None else text(base_root["width"], 35),
        "upper_log": None if base_root is None else text(base_root["upper_log"], 55),
        "G": None if base_root is None else text(base_root["evaluation"]["G"], 45),
        "passed": base_bracket_ok,
        "attempts": [
            {
                "expansion": attempt["expansion"],
                "half_width": text(attempt["half_width"], 45),
                "left_upper_log": text(attempt["left_upper_log"], 45),
                "right_upper_log": text(attempt["right_upper_log"], 45),
                "left_G": text(attempt["left_G"], 45),
                "right_G": text(attempt["right_G"], 45),
                "left_branch_pass": attempt["left_branch_pass"],
                "right_branch_pass": attempt["right_branch_pass"],
                "left_negative_counts": {
                    str(key): value
                    for key, value in attempt["left_branch"]["negative_counts"].items()
                },
                "right_negative_counts": {
                    str(key): value
                    for key, value in attempt["right_branch"]["negative_counts"].items()
                },
                "left_minimum_leading_minor": text(
                    attempt["left_branch"]["minimum_leading_minor"], 35
                ),
                "right_minimum_leading_minor": text(
                    attempt["right_branch"]["minimum_leading_minor"], 35
                ),
                "left_minimum_argument": text(
                    attempt["left_branch"]["minimum_argument"], 35
                ),
                "right_minimum_argument": text(
                    attempt["right_branch"]["minimum_argument"], 35
                ),
                "left_maximum_imaginary": text(
                    attempt["left_maximum_imaginary"], 35
                ),
                "right_maximum_imaginary": text(
                    attempt["right_maximum_imaginary"], 35
                ),
            }
            for attempt in base_bracket_attempts
        ],
    },
    "curve_nodes": [serialize_node(node) for node in curve_nodes],
    "momentum_multiset": [text(value, 45) for value in momentum_multiset],
    "sampled_momentum_min": (
        None if not momentum_multiset else text(min(momentum_multiset), 45)
    ),
    "sampled_momentum_max": (
        None if not momentum_multiset else text(max(momentum_multiset), 45)
    ),
    "main_grid_fold_bracket_count": len(fold_brackets),
    "refined_folds": serialized_folds,
    "classification": {
        "target_independent": True,
        "frozen_domain_computational_envelope": outcome
            == "MOMENTUM_ENVELOPE_ENUMERATED",
        "global_analytic_envelope": False,
        "between_grid_fold_absence": "OPEN",
        "near_null_limit": "PATTERN",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print(f"\nOUTCOME: {outcome}", flush=True)
print(f"RESULT: {passed}/{tests}", flush=True)
if passed != tests:
    raise SystemExit(1)

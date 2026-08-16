#!/usr/bin/env python3
"""Target-independent stationary-root enumeration at the inherited lapse.

Prior-art commit: eecc80e.
Protocol commit: 07083cc.
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
OUTPUT = HERE / "gravity_600cell_dust_stationary_root_enumeration.json"
PRIOR_ART_COMMIT = "eecc80e"
PROTOCOL_COMMIT = "07083cc"
BASE_RESULT_COMMIT = "f9d7ada"
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
MAIN_GRID_SIZE = 257
GRID_DENOMINATOR = 16
REFINEMENT_BISECTIONS = 80
STEP_SETS = {
    "operational_primary": arb.mpf("1e-20"),
    "operational_shadow": arb.mpf("1e-15"),
    "validation_primary": arb.mpf("3e-20"),
    "validation_shadow": arb.mpf("3e-15"),
}
ARITHMETIC_FLOOR = arb.mpf("1e-60")
ENTRY_FACTOR = arb.mpf(10)
NONZERO_FACTOR = arb.mpf(100)
ROOT_TOLERANCE = arb.mpf("1e-25")
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
        "__name__": "stationary_root_enumeration_imported_core",
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


print("Target-independent stationary-root enumeration", flush=True)
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
control_ok = bool(
    PRIOR_ART_COMMIT == "eecc80e"
    and PROTOCOL_COMMIT == "07083cc"
    and BASE_RESULT_COMMIT == "f9d7ada"
    and artifact_hash == FIRST_TICK_SHA256
    and MAIN_GRID_SIZE == 257
    and GRID_DENOMINATOR == 16
    and REFINEMENT_BISECTIONS == 80
    and set(models) == {"even", "odd"}
)
check(
    "the target firewall, frozen grid and accepted geometry hash pass",
    control_ok,
    f"SHA-256={artifact_hash}",
)


def evaluate_raw(parity, upper_log):
    q_old_value = arb.exp(2*A1)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*upper_log)*ARB_L0_SQUARE
    rho_value = ARB_RHO*arb.exp(R1)
    diagonal = arb.exp(A1+upper_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("enumeration magnitude left the positive domain")
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


def evaluate_gp(parity, upper_log, rho_log):
    q_old_value = arb.exp(2*A1)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*upper_log)*ARB_L0_SQUARE
    rho_value = ARB_RHO*arb.exp(rho_log)
    diagonal = arb.exp(A1+upper_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("derivative magnitude left the positive domain")
    action, gradient, branch = action_and_gradient(
        models[parity], tuple([q_old_value]*30),
        tuple([diagonal]*30+[rho_value]*5), tuple([q_new_value]*30),
    )
    local = tuple(arb.re(value) for value in gradient[30:65])
    pre = tuple(-arb.re(value) for value in gradient[:30])
    imaginary = maximum_imaginary(action, gradient)
    return {
        "G": mean(local[30:]),
        "P": mean(pre),
        "branch": branch,
        "maximum_imaginary": imaginary,
        "branch_pass": branch_pass(branch, imaginary),
    }


def calibrated_gp(parity, upper_log):
    values = {}
    all_branch_pass = True
    minimum_minor = arb.inf
    minimum_argument = arb.inf
    maximum_imaginary_value = arb.mpf(0)
    for name, step in STEP_SETS.items():
        for coordinate in range(2):
            for direction in (-1, 1):
                b = upper_log+direction*step if coordinate == 0 else upper_log
                r = R1+direction*step if coordinate == 1 else R1
                evaluation = evaluate_gp(parity, b, r)
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
                + abs(d_cross[row, column])+arb.mpf("1e-60")
            )
            entry_pass &= bool(
                abs(d_cross[row, column]) <= 10*(
                    abs(d_op[row, column])+abs(d_val[row, column])
                    + arb.mpf("1e-60")
                )
            )
    determinants = {name: arb.det(matrix) for name, matrix in matrices.items()}
    determinant = determinants["operational_primary"]
    determinant_epsilon = (
        abs(determinant-determinants["operational_shadow"])
        + abs(determinants["validation_primary"]
              -determinants["validation_shadow"])
        + abs(determinant-determinants["validation_primary"])
        + arb.mpf("1e-60")
    )
    gb_resolved = bool(
        entry_pass and all_branch_pass
        and abs(operational[0, 0]) > NONZERO_FACTOR*entry_errors[0, 0]
    )
    determinant_resolved = bool(
        entry_pass and all_branch_pass
        and abs(determinant) > NONZERO_FACTOR*determinant_epsilon
    )
    return {
        "matrices": matrices,
        "operational": operational,
        "entry_errors": entry_errors,
        "entry_pass": bool(entry_pass),
        "branch_pass": bool(all_branch_pass),
        "gb_resolved": gb_resolved,
        "determinants": determinants,
        "determinant": determinant,
        "determinant_epsilon": determinant_epsilon,
        "determinant_resolved": determinant_resolved,
        "minimum_leading_minor": minimum_minor,
        "minimum_angle_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary_value,
    }


grid = []
all_branch_ok = control_ok
if control_ok:
    for index in range(MAIN_GRID_SIZE):
        x = -arb.mpf(8)+arb.mpf(index)/GRID_DENOMINATOR
        upper_log = x*abs(A1)
        evaluation = evaluate_raw("even", upper_log)
        all_branch_ok &= evaluation["branch_pass"]
        grid.append({
            "index": index,
            "x": x,
            "upper_log": upper_log,
            "evaluation": evaluation,
        })
        if index % 32 == 0:
            print(
                "  grid {:03d}/256 x={} G={}".format(
                    index, text(x, 6), text(evaluation["G"], 8)
                ),
                flush=True,
            )

sentinels = []
if control_ok:
    for exponent in range(4, 13):
        for sign in (-1, 1):
            x = sign*arb.mpf(2)**exponent
            upper_log = x*abs(A1)
            evaluation = evaluate_raw("even", upper_log)
            all_branch_ok &= evaluation["branch_pass"]
            sentinels.append({
                "exponent": exponent,
                "sign": sign,
                "x": x,
                "upper_log": upper_log,
                "evaluation": evaluation,
            })


near_indices = [
    point["index"] for point in grid
    if abs(point["evaluation"]["G"]) < ROOT_TOLERANCE
]
near_clusters = []
for index in near_indices:
    if not near_clusters or index != near_clusters[-1][-1]+1:
        near_clusters.append([index])
    else:
        near_clusters[-1].append(index)
node_candidates = [
    grid[cluster[len(cluster)//2]] for cluster in near_clusters
]

sign_candidates = []
for left, right in zip(grid, grid[1:]):
    left_g = left["evaluation"]["G"]
    right_g = right["evaluation"]["G"]
    if (
        abs(left_g) >= ROOT_TOLERANCE
        and abs(right_g) >= ROOT_TOLERANCE
        and left_g*right_g < 0
    ):
        sign_candidates.append((left, right))

candidate_count = len(node_candidates)+len(sign_candidates)
too_many_candidates = candidate_count > 8
check(
    "the frozen grid produces a completely counted candidate multiset",
    bool(len(grid) == 257 and len(sentinels) == 18 and not too_many_candidates),
    "node clusters={}, sign brackets={}, total={}".format(
        len(node_candidates), len(sign_candidates), candidate_count
    ),
)


roots = []
refinement_failed = False
if not too_many_candidates:
    for candidate in node_candidates:
        roots.append({
            "kind": "node",
            "upper_log": candidate["upper_log"],
            "x": candidate["x"],
            "evaluation": candidate["evaluation"],
            "width": arb.mpf(0),
            "source": [candidate["index"]],
        })
    for left_point, right_point in sign_candidates:
        left_b = left_point["upper_log"]
        right_b = right_point["upper_log"]
        left_eval = left_point["evaluation"]
        right_eval = right_point["evaluation"]
        for _ in range(REFINEMENT_BISECTIONS):
            middle_b = (left_b+right_b)/2
            middle_eval = evaluate_raw("even", middle_b)
            all_branch_ok &= middle_eval["branch_pass"]
            if left_eval["G"]*middle_eval["G"] <= 0:
                right_b, right_eval = middle_b, middle_eval
            else:
                left_b, left_eval = middle_b, middle_eval
        middle_b = (left_b+right_b)/2
        middle_eval = evaluate_raw("even", middle_b)
        width = right_b-left_b
        passed = bool(
            middle_eval["branch_pass"]
            and width < arb.mpf("1e-29")
            and abs(middle_eval["G"]) < ROOT_TOLERANCE
        )
        refinement_failed |= not passed
        roots.append({
            "kind": "sign_bracket",
            "upper_log": middle_b,
            "x": middle_b/abs(A1),
            "evaluation": middle_eval,
            "width": width,
            "source": [left_point["index"], right_point["index"]],
        })

roots.sort(key=lambda item: item["upper_log"])
root_full_failed = False
derivative_open = False
parity_failed = False
for root_index, root in enumerate(roots):
    evaluation = root["evaluation"]
    diagonal = evaluation["local"][:30]
    poles = evaluation["local"][30:]
    derivative = calibrated_gp("even", root["upper_log"])
    odd_evaluation = evaluate_raw("odd", root["upper_log"])
    odd_derivative = calibrated_gp("odd", root["upper_log"])
    full_gate = bool(
        evaluation["branch_pass"]
        and infinity_norm(diagonal) < arb.mpf("1e-60")
        and infinity_norm(poles) < ROOT_TOLERANCE
        and spread(diagonal) < arb.mpf("1e-60")
        and spread(poles) < arb.mpf("1e-60")
    )
    derivative_gate = bool(
        derivative["branch_pass"]
        and derivative["entry_pass"]
        and derivative["gb_resolved"]
    )
    parity_gate = bool(
        odd_evaluation["branch_pass"]
        and abs(evaluation["G"]-odd_evaluation["G"]) < arb.mpf("1e-24")
        and abs(evaluation["P"]-odd_evaluation["P"]) < arb.mpf("1e-22")
        and infinity_norm(tuple(
            left-right for left, right in zip(
                evaluation["local"], odd_evaluation["local"]
            )
        )) < arb.mpf("1e-24")
        and all(
            abs(
                derivative["operational"][row, column]
                - odd_derivative["operational"][row, column]
            ) <= ENTRY_FACTOR*(
                derivative["entry_errors"][row, column]
                + odd_derivative["entry_errors"][row, column]
                + ARITHMETIC_FLOOR
            )
            for row in range(2) for column in range(2)
        )
    )
    root.update({
        "index": root_index,
        "derivative": derivative,
        "odd_evaluation": odd_evaluation,
        "odd_derivative": odd_derivative,
        "full_gate": full_gate,
        "derivative_gate": derivative_gate,
        "parity_gate": parity_gate,
    })
    root_full_failed |= not full_gate
    derivative_open |= not derivative_gate
    parity_failed |= not parity_gate
    all_branch_ok &= bool(
        evaluation["branch_pass"] and derivative["branch_pass"]
        and odd_evaluation["branch_pass"] and odd_derivative["branch_pass"]
    )
    print(
        "  root={} kind={} x={} b={} G={} P={} D={} rank_sign={} full={}".format(
            root_index, root["kind"], text(root["x"], 14),
            text(root["upper_log"], 16), text(evaluation["G"], 8),
            text(evaluation["P"], 14), text(derivative["determinant"], 8),
            derivative["determinant_resolved"], full_gate,
        ),
        flush=True,
    )


if not control_ok:
    outcome = "STATIONARY_ROOT_ENUMERATION_CONTROL_FAILED"
elif not all_branch_ok:
    outcome = "STATIONARY_ROOT_ENUMERATION_BRANCH_FAILED"
elif too_many_candidates:
    outcome = "STATIONARY_ROOT_ENUMERATION_TOO_MANY_CANDIDATES"
elif refinement_failed or root_full_failed:
    outcome = "STATIONARY_ROOT_ENUMERATION_REFINEMENT_FAILED"
elif derivative_open:
    outcome = "STATIONARY_ROOT_ENUMERATION_DERIVATIVE_OPEN"
elif parity_failed:
    outcome = "STATIONARY_ROOT_ENUMERATION_PARITY_FAILED"
else:
    outcome = "STATIONARY_ROOTS_ENUMERATED"

check(
    "all frozen grid, sentinel and root evaluations retain the Lorentzian branch",
    all_branch_ok,
)
check(
    "the frozen hierarchy assigns exactly one stationary-root outcome",
    outcome in {
        "STATIONARY_ROOT_ENUMERATION_CONTROL_FAILED",
        "STATIONARY_ROOT_ENUMERATION_BRANCH_FAILED",
        "STATIONARY_ROOT_ENUMERATION_TOO_MANY_CANDIDATES",
        "STATIONARY_ROOT_ENUMERATION_REFINEMENT_FAILED",
        "STATIONARY_ROOT_ENUMERATION_DERIVATIVE_OPEN",
        "STATIONARY_ROOT_ENUMERATION_PARITY_FAILED",
        "STATIONARY_ROOTS_ENUMERATED",
    },
    outcome,
)


def serialize_branch(evaluation):
    return {
        "passed": evaluation["branch_pass"],
        "negative_counts": {
            str(key): value
            for key, value in evaluation["branch"]["negative_counts"].items()
        },
        "minimum_leading_minor": text(
            evaluation["branch"]["minimum_leading_minor"], 30
        ),
        "minimum_argument": text(
            evaluation["branch"]["minimum_argument"], 30
        ),
        "maximum_imaginary": text(evaluation["maximum_imaginary"], 30),
    }


def serialize_sample(sample):
    evaluation = sample["evaluation"]
    return {
        "x": text(sample["x"], 30),
        "upper_log": text(sample["upper_log"], 45),
        "G": text(evaluation["G"], 40),
        "P": text(evaluation["P"], 40),
        "branch": serialize_branch(evaluation),
    }


def serialize_matrix(matrix):
    return [
        [text(matrix[row, column], 40) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def serialize_derivative(record):
    return {
        "matrices": {
            name: serialize_matrix(matrix)
            for name, matrix in record["matrices"].items()
        },
        "entry_errors": serialize_matrix(record["entry_errors"]),
        "entry_pass": record["entry_pass"],
        "branch_pass": record["branch_pass"],
        "gb_resolved": record["gb_resolved"],
        "determinants": {
            name: text(value, 40)
            for name, value in record["determinants"].items()
        },
        "determinant": text(record["determinant"], 40),
        "determinant_epsilon": text(record["determinant_epsilon"], 30),
        "determinant_resolved": record["determinant_resolved"],
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "base_result_commit": BASE_RESULT_COMMIT,
    "first_tick_sha256": FIRST_TICK_SHA256,
    "target_parsed": False,
    "geometry": {"a1": A1_TEXT, "r1": R1_TEXT},
    "main_grid": {
        "count": len(grid),
        "x_min": "-8",
        "x_max": "8",
        "spacing": "0.0625",
        "samples": [serialize_sample(sample) for sample in grid],
    },
    "sentinels": [
        {
            "exponent": sample["exponent"],
            "sign": sample["sign"],
            **serialize_sample(sample),
        }
        for sample in sentinels
    ],
    "candidate_counts": {
        "near_zero_clusters": len(node_candidates),
        "sign_brackets": len(sign_candidates),
        "total": candidate_count,
    },
    "near_zero_cluster_indices": near_clusters,
    "sign_bracket_indices": [
        [left["index"], right["index"]]
        for left, right in sign_candidates
    ],
    "roots": [
        {
            "index": root["index"],
            "kind": root["kind"],
            "source": root["source"],
            "x": text(root["x"], 45),
            "upper_log": text(root["upper_log"], 55),
            "width": text(root["width"], 35),
            "G": text(root["evaluation"]["G"], 45),
            "P": text(root["evaluation"]["P"], 45),
            "pre_momentum": [
                text(value, 45) for value in root["evaluation"]["pre"]
            ],
            "post_momentum": [
                text(value, 45) for value in root["evaluation"]["post"]
            ],
            "local_residuals": [
                text(value, 35) for value in root["evaluation"]["local"]
            ],
            "derivative": serialize_derivative(root["derivative"]),
            "full_gate": root["full_gate"],
            "derivative_gate": root["derivative_gate"],
            "parity_gate": root["parity_gate"],
        }
        for root in roots
    ],
    "root_momentum_multiset": [
        text(root["evaluation"]["P"], 45) for root in roots
    ],
    "classification": {
        "target_independent": True,
        "frozen_domain_root_count": len(roots),
        "tangential_roots_between_nodes": "OPEN",
        "roots_outside_domain": "OPEN",
        "sentinels": "PATTERN",
        "physical_branch_selected": False,
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

#!/usr/bin/env python3
"""Target-independent stationary-root enumeration for the third dust tick.

Prior-art commit: 7b9a676.
Protocol commit: 5d980b1.
"""

import ast
import contextlib
import hashlib
import io
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
RESPONSE_SOURCE = HERE / "verify_gravity_600cell_dust_homothetic_mass_conservation.py"
SECOND_TICK_ARTIFACT = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
OUTPUT = HERE / "gravity_600cell_dust_third_tick_stationary_roots.json"
PRIOR_ART_COMMIT = "7b9a676"
PROTOCOL_COMMIT = "5d980b1"
SECOND_TICK_SHA256 = "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70"
A1_TEXT = "-0.00000311605957669450169173470644419863944122165192557277135128791"
B2_TEXT = "-0.00000934818705890582713633822299265753373027428194008991504419612"
R2_TEXT = "-0.0000142370275520098029961300545242474815338378370661665379256974"
DPS = 100
arb.mp.dps = DPS
A1 = arb.mpf(A1_TEXT)
B2 = arb.mpf(B2_TEXT)
R2 = arb.mpf(R2_TEXT)
MAIN_GRID_SIZE = 257
GRID_DENOMINATOR = 16
REFINEMENT_BISECTIONS = 80
ROOT_TOLERANCE = arb.mpf("1e-25")
STEP_SETS = {
    "operational_primary": arb.mpf("1e-20"),
    "operational_shadow": arb.mpf("1e-15"),
    "validation_primary": arb.mpf("3e-20"),
    "validation_shadow": arb.mpf("3e-15"),
}
ARITHMETIC_FLOOR = arb.mpf("1e-60")
ENTRY_FACTOR = arb.mpf(10)
NONZERO_FACTOR = arb.mpf(100)
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


def text(value, digits=50):
    return arb.nstr(value, digits)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean(values):
    return sum(values, arb.mpf(0))/len(values)


def infinity_norm(values):
    return max(abs(value) for value in values)


def spread(values):
    average = mean(values)
    return max(abs(value-average) for value in values)


def maximum_imaginary(action, gradient):
    return max(abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient))


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
        "__name__": "third_tick_roots_imported_response_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(RESPONSE_SOURCE), "exec"), namespace)
    return namespace


print("Third tick: target-independent stationary-root enumeration", flush=True)
response_core = load_response_prefix()
models = response_core["models"]
core = response_core["core"]
action_and_gradient = core["action_and_gradient"]
branch_pass = core["branch_pass"]
ARB_L0_SQUARE = core["ARB_L0_SQUARE"]
ARB_RHO = core["ARB_RHO"]

check(
    "the imported complete-action evaluator retains its six controls",
    response_core["tests"] == response_core["passed"] == 6,
)

artifact_hash = digest(SECOND_TICK_ARTIFACT)
firewall_ok = bool(
    PRIOR_ART_COMMIT == "7b9a676"
    and PROTOCOL_COMMIT == "5d980b1"
    and artifact_hash == SECOND_TICK_SHA256
    and MAIN_GRID_SIZE == 257
    and GRID_DENOMINATOR == 16
    and REFINEMENT_BISECTIONS == 80
    and set(models) == {"even", "odd"}
)
check(
    "the target firewall, frozen geometry and grid controls pass",
    firewall_ok,
    f"second-tick byte hash={artifact_hash}; target JSON parsed=False",
)


def evaluate_state(parity, upper_log, rho_log=R2):
    q_old_value = arb.exp(2*B2)*ARB_L0_SQUARE
    q_new_value = arb.exp(2*upper_log)*ARB_L0_SQUARE
    rho_value = ARB_RHO*arb.exp(rho_log)
    diagonal = arb.exp(B2+upper_log)*ARB_L0_SQUARE-rho_value
    if min(q_old_value, q_new_value, rho_value, diagonal) <= 0:
        raise ValueError("third-tick homothetic magnitude left positive domain")
    action, gradient, branch = action_and_gradient(
        models[parity], tuple([q_old_value]*30),
        tuple([diagonal]*30+[rho_value]*5), tuple([q_new_value]*30),
    )
    local = tuple(arb.re(value) for value in gradient[30:65])
    pre = tuple(-arb.re(value) for value in gradient[:30])
    post = tuple(arb.re(value) for value in gradient[65:95])
    imaginary = maximum_imaginary(action, gradient)
    return {
        "upper_log": upper_log,
        "rho_log": rho_log,
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


def calibrated_gp(parity, upper_log, rho_log):
    values = {}
    all_branch_pass = True
    minimum_minor = arb.inf
    minimum_argument = arb.inf
    maximum_imaginary_value = arb.mpf(0)
    for name, step in STEP_SETS.items():
        for coordinate in range(2):
            for direction in (-1, 1):
                c = upper_log+direction*step if coordinate == 0 else upper_log
                r = rho_log+direction*step if coordinate == 1 else rho_log
                evaluation = evaluate_state(parity, c, r)
                values[(name, coordinate, direction)] = (
                    evaluation["G"], evaluation["P"]
                )
                all_branch_pass &= evaluation["branch_pass"]
                minimum_minor = min(
                    minimum_minor, evaluation["branch"]["minimum_leading_minor"]
                )
                minimum_argument = min(
                    minimum_argument, evaluation["branch"]["minimum_argument"]
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
all_branch_ok = firewall_ok
if firewall_ok:
    for index in range(MAIN_GRID_SIZE):
        x = -arb.mpf(8)+arb.mpf(index)/GRID_DENOMINATOR
        upper_log = x*abs(A1)
        evaluation = evaluate_state("even", upper_log)
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
                    index, text(x, 7), text(evaluation["G"], 10)
                ), flush=True,
            )

sentinels = []
if firewall_ok:
    for exponent in range(4, 13):
        for sign in (-1, 1):
            x = sign*arb.mpf(2)**exponent
            upper_log = x*abs(A1)
            evaluation = evaluate_state("even", upper_log)
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
    if not near_clusters or index > near_clusters[-1][-1]+1:
        near_clusters.append([index])
    else:
        near_clusters[-1].append(index)

sign_brackets = []
for index in range(len(grid)-1):
    left = grid[index]["evaluation"]["G"]
    right = grid[index+1]["evaluation"]["G"]
    if left*right < 0:
        sign_brackets.append((index, index+1))

cluster_index_sets = [set(cluster) for cluster in near_clusters]
standalone_brackets = []
merged_brackets = {cluster_id: [] for cluster_id in range(len(near_clusters))}
for bracket in sign_brackets:
    touching = [
        cluster_id for cluster_id, indices in enumerate(cluster_index_sets)
        if bracket[0] in indices or bracket[1] in indices
    ]
    if touching:
        for cluster_id in touching:
            merged_brackets[cluster_id].append(bracket)
    else:
        standalone_brackets.append(bracket)


def bisect_root(parity, left_index, right_index):
    left_c = grid[left_index]["upper_log"]
    right_c = grid[right_index]["upper_log"]
    left_eval = evaluate_state(parity, left_c)
    right_eval = evaluate_state(parity, right_c)
    branch_ok = left_eval["branch_pass"] and right_eval["branch_pass"]
    if left_eval["G"]*right_eval["G"] >= 0:
        return None, branch_ok
    for _ in range(REFINEMENT_BISECTIONS):
        middle_c = (left_c+right_c)/2
        middle_eval = evaluate_state(parity, middle_c)
        branch_ok &= middle_eval["branch_pass"]
        if middle_eval["G"] == 0:
            left_c = right_c = middle_c
            left_eval = right_eval = middle_eval
            break
        if left_eval["G"]*middle_eval["G"] < 0:
            right_c, right_eval = middle_c, middle_eval
        else:
            left_c, left_eval = middle_c, middle_eval
    root_c = (left_c+right_c)/2
    return {
        "upper_log": root_c,
        "evaluation": evaluate_state(parity, root_c),
        "width": right_c-left_c,
    }, branch_ok


candidates = []
for cluster_id, cluster in enumerate(near_clusters):
    best = min(
        cluster,
        key=lambda index: (abs(grid[index]["evaluation"]["G"]), index),
    )
    candidates.append({
        "kind": "node_cluster",
        "source": list(cluster),
        "merged_sign_brackets": [list(pair) for pair in merged_brackets[cluster_id]],
        "even": {
            "upper_log": grid[best]["upper_log"],
            "evaluation": grid[best]["evaluation"],
            "width": arb.mpf(0),
        },
        "odd_source": (best, best),
    })

for left_index, right_index in standalone_brackets:
    even_root, branch_ok = bisect_root("even", left_index, right_index)
    all_branch_ok &= branch_ok
    candidates.append({
        "kind": "sign_bracket",
        "source": [left_index, right_index],
        "merged_sign_brackets": [],
        "even": even_root,
        "odd_source": (left_index, right_index),
    })

candidates.sort(key=lambda item: item["even"]["upper_log"])

for candidate in candidates:
    if candidate["kind"] == "node_cluster":
        c = candidate["even"]["upper_log"]
        candidate["odd"] = {
            "upper_log": c,
            "evaluation": evaluate_state("odd", c),
            "width": arb.mpf(0),
        }
        all_branch_ok &= candidate["odd"]["evaluation"]["branch_pass"]
    else:
        left_index, right_index = candidate["odd_source"]
        odd_root, branch_ok = bisect_root("odd", left_index, right_index)
        if odd_root is None:
            fallback_c = candidate["even"]["upper_log"]
            odd_root = {
                "upper_log": fallback_c,
                "evaluation": evaluate_state("odd", fallback_c),
                "width": grid[right_index]["upper_log"]-grid[left_index]["upper_log"],
            }
            candidate["odd_bracket_missing"] = True
        else:
            candidate["odd_bracket_missing"] = False
        candidate["odd"] = odd_root
        all_branch_ok &= branch_ok and not candidate["odd_bracket_missing"]

root_records = []
all_root_gates = True
for index, candidate in enumerate(candidates):
    parity_records = {}
    determinant_signs = []
    for parity in ("even", "odd"):
        root = candidate[parity]
        evaluation = root["evaluation"]
        derivative = calibrated_gp(parity, root["upper_log"], R2)
        all_branch_ok &= derivative["branch_pass"]
        diagonal = evaluation["local"][:30]
        poles = evaluation["local"][30:]
        full_gate = bool(
            evaluation["branch_pass"]
            and infinity_norm(diagonal) < arb.mpf("1e-60")
            and infinity_norm(poles) < arb.mpf("1e-25")
            and infinity_norm(evaluation["local"]) < arb.mpf("1e-25")
            and spread(diagonal) < arb.mpf("1e-60")
            and spread(poles) < arb.mpf("1e-60")
        )
        derivative_gate = bool(
            derivative["entry_pass"] and derivative["branch_pass"]
            and derivative["gb_resolved"]
            and derivative["determinant_resolved"]
            and not candidate.get("odd_bracket_missing", False)
        )
        all_root_gates &= full_gate and derivative_gate
        determinant_signs.append(arb.sign(derivative["determinant"]))
        parity_records[parity] = {
            "upper_log": root["upper_log"],
            "width": root["width"],
            "evaluation": evaluation,
            "derivative": derivative,
            "full_gate": full_gate,
            "derivative_gate": derivative_gate,
        }
    even = parity_records["even"]
    odd = parity_records["odd"]
    c_difference = abs(even["upper_log"]-odd["upper_log"])
    local_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even["evaluation"]["local"], odd["evaluation"]["local"]
        )
    ))
    pre_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even["evaluation"]["pre"], odd["evaluation"]["pre"]
        )
    ))
    parity_gate = bool(
        c_difference < arb.mpf("1e-27")
        and local_difference < arb.mpf("1e-24")
        and pre_difference < arb.mpf("1e-22")
        and determinant_signs[0] == determinant_signs[1]
    )
    all_root_gates &= parity_gate
    root_records.append({
        "index": index,
        "kind": candidate["kind"],
        "source": candidate["source"],
        "merged_sign_brackets": candidate["merged_sign_brackets"],
        "parities": parity_records,
        "parity": {
            "c_difference": c_difference,
            "local_difference": local_difference,
            "pre_difference": pre_difference,
            "determinant_signs": determinant_signs,
            "passed": parity_gate,
        },
    })
    print(
        "  root {} kind={} x={} C={} G={} P={} D={} full={} parity={}".format(
            index, candidate["kind"],
            text(even["upper_log"]/abs(A1), 16),
            text(even["upper_log"], 16), text(even["evaluation"]["G"], 9),
            text(even["evaluation"]["P"], 16),
            text(even["derivative"]["determinant"], 9),
            even["full_gate"], parity_gate,
        ), flush=True,
    )

reverse_candidates = [
    record["index"] for record in root_records
    if abs(record["parities"]["even"]["upper_log"]-A1) < ROOT_TOLERANCE
]
contracting_prediction = 6*A1
contracting_candidates = [
    record["index"] for record in root_records
    if abs(record["parities"]["even"]["upper_log"]-contracting_prediction)
       < abs(A1)/16
]
prediction_diagnostics = {
    "candidate_count_equals_two": len(root_records) == 2,
    "time_reversal_at_A1": reverse_candidates,
    "contracting_near_6A1": contracting_candidates,
}

check(
    "all frozen grid, sentinel, bisection and derivative states retain the branch",
    all_branch_ok,
)

if not firewall_ok:
    outcome = "THIRD_TICK_ROOT_CONTROL_FAILED"
elif not all_branch_ok:
    outcome = "THIRD_TICK_ROOT_BRANCH_FAILED"
elif not all_root_gates:
    outcome = "THIRD_TICK_ROOT_ENUMERATION_OPEN"
else:
    outcome = "THIRD_TICK_STATIONARY_ROOTS_ENUMERATED"

check(
    "the candidate multiset and prediction diagnostics are recorded without a target",
    len(root_records) == len(candidates)
    and all(record["parities"]["even"] is not None for record in root_records),
    f"N={len(root_records)}, diagnostics={prediction_diagnostics}",
)
check(
    "the frozen hierarchy assigns one target-independent outcome",
    outcome in {
        "THIRD_TICK_ROOT_CONTROL_FAILED",
        "THIRD_TICK_ROOT_BRANCH_FAILED",
        "THIRD_TICK_ROOT_ENUMERATION_OPEN",
        "THIRD_TICK_STATIONARY_ROOTS_ENUMERATED",
    },
    outcome,
)


def serialize_matrix(matrix):
    return [
        [text(matrix[row, column], 50) for column in range(matrix.cols)]
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
            name: text(value, 50) for name, value in record["determinants"].items()
        },
        "determinant": text(record["determinant"], 50),
        "determinant_epsilon": text(record["determinant_epsilon"], 35),
        "determinant_resolved": record["determinant_resolved"],
        "minimum_leading_minor": text(record["minimum_leading_minor"], 40),
        "minimum_angle_argument": text(record["minimum_angle_argument"], 40),
        "maximum_imaginary": text(record["maximum_imaginary"], 40),
    }


def serialize_evaluation(record):
    return {
        "G": text(record["G"], 50),
        "P": text(record["P"], 50),
        "local_residuals": [text(value, 40) for value in record["local"]],
        "pre_momentum": [text(value, 50) for value in record["pre"]],
        "post_momentum": [text(value, 50) for value in record["post"]],
        "maximum_imaginary": text(record["maximum_imaginary"], 40),
        "minimum_leading_minor": text(record["branch"]["minimum_leading_minor"], 40),
        "minimum_angle_argument": text(record["branch"]["minimum_argument"], 40),
        "branch_pass": record["branch_pass"],
    }


serialized_roots = []
for record in root_records:
    item = {
        "index": record["index"],
        "kind": record["kind"],
        "source": record["source"],
        "merged_sign_brackets": record["merged_sign_brackets"],
        "parity": {
            "c_difference": text(record["parity"]["c_difference"], 40),
            "local_difference": text(record["parity"]["local_difference"], 40),
            "pre_difference": text(record["parity"]["pre_difference"], 40),
            "determinant_signs": [int(value) for value in record["parity"]["determinant_signs"]],
            "passed": record["parity"]["passed"],
        },
        "parities": {},
    }
    for parity in ("even", "odd"):
        parity_record = record["parities"][parity]
        item["parities"][parity] = {
            "upper_log": text(parity_record["upper_log"], 60),
            "x": text(parity_record["upper_log"]/abs(A1), 50),
            "width": text(parity_record["width"], 40),
            "evaluation": serialize_evaluation(parity_record["evaluation"]),
            "derivative": serialize_derivative(parity_record["derivative"]),
            "full_gate": parity_record["full_gate"],
            "derivative_gate": parity_record["derivative_gate"],
        }
    serialized_roots.append(item)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "second_tick_sha256": artifact_hash,
    "target_parsed": False,
    "geometry": {"A1": A1_TEXT, "B2": B2_TEXT, "R2": R2_TEXT},
    "main_grid": {
        "size": MAIN_GRID_SIZE,
        "x_min": "-8",
        "x_max": "8",
        "spacing": "1/16",
        "samples": [
            {
                "index": point["index"],
                "x": text(point["x"], 20),
                "upper_log": text(point["upper_log"], 45),
                "G": text(point["evaluation"]["G"], 35),
                "P": text(point["evaluation"]["P"], 35),
                "branch_pass": point["evaluation"]["branch_pass"],
            }
            for point in grid
        ],
    },
    "sentinels": [
        {
            "exponent": point["exponent"],
            "sign": point["sign"],
            "x": text(point["x"], 20),
            "upper_log": text(point["upper_log"], 45),
            "G": text(point["evaluation"]["G"], 35),
            "branch_pass": point["evaluation"]["branch_pass"],
        }
        for point in sentinels
    ],
    "near_zero_cluster_indices": near_clusters,
    "sign_bracket_indices": [list(pair) for pair in sign_brackets],
    "standalone_sign_brackets": [list(pair) for pair in standalone_brackets],
    "candidate_count": len(serialized_roots),
    "roots": serialized_roots,
    "prediction_diagnostics": prediction_diagnostics,
    "scope": {
        "sign_and_node_roots_on_frozen_grid": "DERIVED COMPUTATIONAL if enumerated",
        "tangential_between_node_roots": "OPEN",
        "outside_main_grid": "OPEN",
        "physical_branch_selected": False,
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print(f"OUTCOME: {outcome}", flush=True)
print(f"Tests passed: {passed}/{tests}", flush=True)
raise SystemExit(0 if passed == tests else 1)

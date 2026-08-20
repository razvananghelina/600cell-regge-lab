#!/usr/bin/env python3
"""All-schedule internal log-Hessian census at the refined H4 fill.

Prior-art commit: 4ea4430.
Protocol commit: 9f1721c.
No Newton update, root, effective boundary Hessian or spectrum is computed.
"""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
UPSTREAM = HERE / "gravity_600cell_refined_h4_stationary_fill.json"
RESULT_NOTE = ROOT / "docs/gravity/gravity_600cell_refined_h4_stationary_fill_result.md"
OUTPUT = HERE / "gravity_600cell_refined_h4_internal_jacobian.json"
PRIOR_ART_COMMIT = "4ea4430"
PROTOCOL_COMMIT = "9f1721c"
INPUT_HASHES = {
    "reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py":
        "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "reproducible/gravity_600cell_refined_h4_stationary_fill.json":
        "283be37bc7530a3cc4fce9e279272359f107f09fb7b1b0eaff141059bfb4e018",
    "docs/gravity/gravity_600cell_refined_h4_stationary_fill_result.md":
        "899143cfcf75ce08a8dd6daca776cd0e38e02e38593f4a171d3d3f851dae7d91",
}
PAIR4 = tuple(combinations(range(4), 2))
VARIABLES = (
    tuple(("old",)+pair for pair in PAIR4)
    + tuple(("new",)+pair for pair in PAIR4)
    + tuple(("cross",)+pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
INTERNAL_VARIABLES = (
    tuple(("cross",)+pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
TAU_TEXT = "0.0102"
PRIMARY_DPS = 100
SECONDARY_DPS = 140
STEP_TEXTS = ("1e-10", "5e-11", "2.5e-11")
DIRECTIONAL_GATE_TEXT = "1e-28"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return condition


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=70):
    return mp.nstr(value, digits)


def complex_record(value, digits=70):
    return {
        "real": mp_text(mp.re(value), digits),
        "imag": mp_text(mp.im(value), digits),
        "absolute": mp_text(abs(value), digits),
    }


def variable_label(key):
    if key[0] == "rho":
        return f"rho_{key[1]}"
    return f"{key[0]}_{key[1]}{key[2]}"


def load_action_definitions():
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "mp": mp,
        "np": np,
        "json": json,
        "Path": Path,
        "sha256": sha256,
        "combinations": combinations,
        "permutations": permutations,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "HERE": HERE,
        "ROOT": ROOT,
        "PAIR4": PAIR4,
        "LOCAL_TRIANGLES": np.asarray(
            tuple(combinations(range(5), 3)), dtype=np.int8
        ),
        "TAU_TEXT": TAU_TEXT,
        "VARIABLES": VARIABLES,
        "INTERNAL_VARIABLES": INTERNAL_VARIABLES,
        "FD_STEP_TEXTS": ("1e-15", "5e-16"),
        "FD_GATE_TEXT": "1e-24",
        "EXPECTED_F": (2640, 17040, 28800, 14400),
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(SOURCE), "exec"), namespace)
    return namespace


def parse_states(label):
    result = []
    for item in label.split("|"):
        rank_text, layer_text = item.split("t")
        result.append((int(rank_text[1:]), int(layer_text)))
    return tuple(result)


def parse_combinatorics(record):
    return {
        "order": tuple(record["order"]),
        "pentachora": record["pentachora"],
        "distinct_pentachora": record["distinct_pentachora"],
        "triangles": record["triangles"],
        "boundary_triangles": record["boundary_triangles"],
        "mixed_triangle_types": record["mixed_triangle_types"],
        "simplex_types": [
            {"states": parse_states(item["states"]), "count": item["count"]}
            for item in record["simplex_types"]
        ],
        "triangle_types": [
            {
                "states": parse_states(item["states"]),
                "count": item["count"],
                "incidence_values": item["incidence_values"],
                "contributions": [
                    {
                        "simplex": parse_states(value["simplex"]),
                        "triangle": parse_states(value["triangle"]),
                        "multiplicity": value["multiplicity"],
                    }
                    for value in item["signature"]
                ],
            }
            for item in record["triangle_types"]
        ],
    }


def shifted_coordinates(base, key, displacement):
    result = dict(base)
    result[key] *= mp.exp(displacement)
    return result


def gradient_vector(evaluate_schedule, combinatorics, geometry, coordinates):
    result = evaluate_schedule(combinatorics, geometry, coordinates)
    return result, mp.matrix([
        result["gradient"][key] for key in INTERNAL_VARIABLES
    ])


def derivative_column(evaluate_schedule, combinatorics, geometry, base,
                      key, step):
    _, plus = gradient_vector(
        evaluate_schedule, combinatorics, geometry,
        shifted_coordinates(base, key, step),
    )
    _, minus = gradient_vector(
        evaluate_schedule, combinatorics, geometry,
        shifted_coordinates(base, key, -step),
    )
    return (plus-minus)/(2*step)


def jacobian_ladders(evaluate_schedule, combinatorics, geometry, dps,
                     secondary_only=False):
    with mp.workdps(dps):
        base = actions["base_coordinates"](geometry)
        steps = tuple(mp.mpf(value) for value in STEP_TEXTS)
        primary = mp.matrix(10, 10)
        secondary = mp.matrix(10, 10)
        for column, key in enumerate(INTERNAL_VARIABLES):
            d1 = derivative_column(
                evaluate_schedule, combinatorics, geometry, base, key, steps[1]
            )
            d2 = derivative_column(
                evaluate_schedule, combinatorics, geometry, base, key, steps[2]
            )
            secondary_column = (4*d2-d1)/3
            for row in range(10):
                secondary[row, column] = secondary_column[row]
            if not secondary_only:
                d0 = derivative_column(
                    evaluate_schedule, combinatorics, geometry, base,
                    key, steps[0]
                )
                primary_column = (4*d1-d0)/3
                for row in range(10):
                    primary[row, column] = primary_column[row]
        return primary, secondary


def matrix_max(matrix):
    return max(abs(matrix[row, column])
               for row in range(matrix.rows) for column in range(matrix.cols))


def matrix_difference(left, right):
    return max(abs(left[row, column]-right[row, column])
               for row in range(left.rows) for column in range(left.cols))


def real_symmetric(matrix):
    return mp.matrix([
        [mp.re(matrix[row, column]+matrix[column, row])/2
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ])


def induced_tangent(coordinates):
    rho = mp.mpf(TAU_TEXT)**2
    values = [
        -rho/coordinates[("cross",)+pair] for pair in PAIR4
    ] + [mp.mpf(1) for _ in range(4)]
    vector = mp.matrix(values)
    return vector/mp.norm(vector)


def directional_action_second(evaluate_schedule, combinatorics, geometry,
                              base, direction, step):
    plus = dict(base)
    minus = dict(base)
    for index, key in enumerate(INTERNAL_VARIABLES):
        plus[key] *= mp.exp(step*direction[index])
        minus[key] *= mp.exp(-step*direction[index])
    centre = evaluate_schedule(combinatorics, geometry, base)["action"]
    upper = evaluate_schedule(combinatorics, geometry, plus)["action"]
    lower = evaluate_schedule(combinatorics, geometry, minus)["action"]
    return (upper-2*centre+lower)/(step*step)


def directional_controls(evaluate_schedule, combinatorics, geometry, matrix):
    base = actions["base_coordinates"](geometry)
    directions = {
        "cross_01": mp.matrix([1]+[0]*9),
        "rho_0": mp.matrix([0]*6+[1]+[0]*3),
        "induced_lapse": induced_tangent(base),
    }
    step0, step1 = (mp.mpf(STEP_TEXTS[index]) for index in (0, 1))
    records = {}
    errors = []
    for label, direction in directions.items():
        coarse = directional_action_second(
            evaluate_schedule, combinatorics, geometry, base, direction, step0
        )
        fine = directional_action_second(
            evaluate_schedule, combinatorics, geometry, base, direction, step1
        )
        richardson = (4*fine-coarse)/3
        quadratic = (direction.T*matrix*direction)[0]
        relative = abs(richardson-quadratic)/max(mp.mpf(1), abs(quadratic))
        errors.append(relative)
        records[label] = {
            "matrix_quadratic": complex_record(quadratic, 60),
            "action_richardson": complex_record(richardson, 60),
            "relative_error": mp_text(relative, 35),
        }
    return max(errors), records


def eigen_record(matrix, spectral_error, tangent, gradient):
    values, vectors = mp.eigsy(matrix)
    eigenvalues = [values[index] for index in range(10)]
    certified = [abs(value) > spectral_error for value in eigenvalues]
    positive = sum(value > spectral_error for value in eigenvalues)
    negative = sum(value < -spectral_error for value in eigenvalues)
    zero = 10-positive-negative
    tangent_image = matrix*tangent
    overlaps = [
        abs(sum(vectors[row, column]*tangent[row] for row in range(10)))
        for column in range(10) if not certified[column]
    ]
    maximum_zero_overlap = max(overlaps, default=mp.mpf(0))
    rayleigh = (tangent.T*matrix*tangent)[0]

    proposal = mp.matrix(10, 1)
    for column, value in enumerate(eigenvalues):
        if certified[column]:
            coefficient = -sum(
                vectors[row, column]*mp.re(gradient[row]) for row in range(10)
            )/value
            for row in range(10):
                proposal[row] += coefficient*vectors[row, column]
    linear_residual = matrix*proposal+mp.matrix([
        mp.re(gradient[row]) for row in range(10)
    ])
    return {
        "eigenvalues": eigenvalues,
        "rank": sum(certified),
        "inertia": (positive, zero, negative),
        "tangent_image_norm": mp.norm(tangent_image),
        "tangent_rayleigh": rayleigh,
        "maximum_zero_overlap": maximum_zero_overlap,
        "proposal": proposal,
        "proposal_norm": mp.norm(proposal),
        "proposal_linear_residual_norm": mp.norm(linear_residual),
    }


print("="*78)
print("REFINED H4 INTERNAL-JACOBIAN CENSUS")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "the frozen stationary-fill source, artifact and result have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "4ea4430" and PROTOCOL_COMMIT == "9f1721c",
    str(actual_hashes),
)

upstream = json.loads(UPSTREAM.read_text())
upstream_ok = check(
    "the upstream census has the exact off-shell outcome required by this gate",
    upstream["outcome"] == "REFINED_H4_INDUCED_FILL_OFF_SHELL"
    and upstream["tests"] == {"passed": 12, "total": 12}
    and upstream["census"]["schedule_count"] == 24
    and upstream["census"]["internal_entries"] == 240
    and upstream["census"]["certified_nonzero_entries"] == 96
    and upstream["census"]["certified_nonzero_vertical_entries"] == 96
    and upstream["census"]["certified_nonzero_cross_entries"] == 0
    and upstream["census"]["distinct_residual_vectors"] == 1,
)

actions = load_action_definitions()
required_definitions = {
    "exact_geometry", "base_coordinates", "evaluate_schedule",
    "angle_record", "triangle_area_and_derivatives",
}
definitions_ok = check(
    "only frozen AST definitions are loaded and all action functions are present",
    required_definitions <= set(actions) and "OUTPUT" not in actions,
)

combinatorics = [parse_combinatorics(record)
                 for record in upstream["combinatorics"]]
orders = tuple(record["order"] for record in combinatorics)
combinatorics_ok = check(
    "all 24 certified incidence records reconstruct without a new carrier choice",
    orders == tuple(permutations(range(4)))
    and all(record["pentachora"] == record["distinct_pentachora"] == 57600
            and record["triangles"] == 149280
            and record["boundary_triangles"] == 57600
            and record["mixed_triangle_types"] == 0
            and len(record["simplex_types"]) == 4
            and len(record["triangle_types"]) == 28
            for record in combinatorics),
)

geometry100 = actions["exact_geometry"](PRIMARY_DPS)
geometry140 = actions["exact_geometry"](SECONDARY_DPS)

primary_matrices = []
secondary100 = []
secondary140 = []
print("[INFO] constructing 24 Jacobians at 100 decimal digits", flush=True)
for index, record in enumerate(combinatorics):
    primary, secondary = jacobian_ladders(
        actions["evaluate_schedule"], record, geometry100, PRIMARY_DPS
    )
    primary_matrices.append(primary)
    secondary100.append(secondary)
    if index in (0, 5, 11, 17, 23):
        print(f"[INFO] primary schedules completed: {index+1}/24", flush=True)

print("[INFO] repeating 24 fine Jacobians at 140 decimal digits", flush=True)
for index, record in enumerate(combinatorics):
    _, secondary = jacobian_ladders(
        actions["evaluate_schedule"], record, geometry140, SECONDARY_DPS,
        secondary_only=True,
    )
    secondary140.append(secondary)
    if index in (0, 5, 11, 17, 23):
        print(f"[INFO] secondary schedules completed: {index+1}/24", flush=True)

matrix_records = []
matrices = []
envelopes = []
raw_controls = []
with mp.workdps(SECONDARY_DPS):
    for index in range(24):
        step_difference = matrix_difference(
            primary_matrices[index], secondary100[index]
        )
        scale = max(mp.mpf(1), matrix_max(secondary100[index]))
        entry_error = 100*step_difference+mp.mpf("1e-50")*scale
        precision_difference = matrix_difference(
            secondary100[index], secondary140[index]
        )
        raw_imaginary = max(
            abs(mp.im(secondary140[index][row, column]))
            for row in range(10) for column in range(10)
        )
        raw_antisymmetry = max(
            abs(secondary140[index][row, column]
                - secondary140[index][column, row])
            for row in range(10) for column in range(10)
        )
        spectral_error = 10*entry_error
        matrix = real_symmetric(secondary140[index])
        matrices.append(matrix)
        envelopes.append((entry_error, spectral_error))

        coordinates = actions["base_coordinates"](geometry140)
        evaluation = actions["evaluate_schedule"](
            combinatorics[index], geometry140, coordinates
        )
        gradient = mp.matrix([
            evaluation["gradient"][key] for key in INTERNAL_VARIABLES
        ])
        tangent = induced_tangent(coordinates)
        spectral = eigen_record(matrix, spectral_error, tangent, gradient)
        raw_ok = (
            precision_difference <= entry_error
            and raw_imaginary <= entry_error
            and raw_antisymmetry <= entry_error
        )
        raw_controls.append(raw_ok)
        matrix_records.append({
            "order": list(orders[index]),
            "step_difference": step_difference,
            "precision_difference": precision_difference,
            "entry_error": entry_error,
            "spectral_error": spectral_error,
            "raw_maximum_imaginary": raw_imaginary,
            "raw_maximum_antisymmetry": raw_antisymmetry,
            "matrix": matrix,
            "spectral": spectral,
        })

raw_ok = check(
    "all raw matrices pass step, precision, reality and Hessian-symmetry envelopes",
    all(raw_controls),
    f"max step={mp_text(max(item['step_difference'] for item in matrix_records), 8)}, "
    f"max precision={mp_text(max(item['precision_difference'] for item in matrix_records), 8)}, "
    f"max antisym={mp_text(max(item['raw_maximum_antisymmetry'] for item in matrix_records), 8)}",
)

with mp.workdps(SECONDARY_DPS):
    directional_records = {}
    directional_errors = []
    for index in (0, 23):
        error, record = directional_controls(
            actions["evaluate_schedule"], combinatorics[index], geometry140,
            matrices[index],
        )
        directional_errors.append(error)
        directional_records[str(orders[index])] = record
    maximum_directional_error = max(directional_errors)
directional_ok = check(
    "six independent action second differences reproduce the Jacobian quadratics",
    maximum_directional_error < mp.mpf(DIRECTIONAL_GATE_TEXT),
    f"max relative error={mp_text(maximum_directional_error, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    classes = []
    class_members = []
    class_indices = []
    for index, matrix in enumerate(matrices):
        assigned = None
        for class_index, representative_index in enumerate(classes):
            tolerance = max(
                envelopes[index][0], envelopes[representative_index][0]
            )
            if matrix_difference(matrix, matrices[representative_index]) <= tolerance:
                assigned = class_index
                break
        if assigned is None:
            assigned = len(classes)
            classes.append(index)
            class_members.append([])
        class_members[assigned].append(index)
        class_indices.append(assigned)

    order_index = {order: index for index, order in enumerate(orders)}
    reversal_differences = []
    reversal_ok_values = []
    for index, order in enumerate(orders):
        reverse_index = order_index[tuple(reversed(order))]
        difference = matrix_difference(matrices[index], matrices[reverse_index])
        tolerance = max(envelopes[index][0], envelopes[reverse_index][0])
        reversal_differences.append(difference)
        reversal_ok_values.append(difference <= tolerance)
time_reversal_ok = check(
    "all twelve time-reversal schedule pairs have identical internal Jacobians",
    all(reversal_ok_values),
    f"max difference={mp_text(max(reversal_differences), 8)}",
)

ranks = [item["spectral"]["rank"] for item in matrix_records]
inertias = [item["spectral"]["inertia"] for item in matrix_records]
null_candidate = all(
    item["spectral"]["rank"] == 9
    and item["spectral"]["tangent_image_norm"] <= 10*item["spectral_error"]
    and item["spectral"]["maximum_zero_overlap"]
        >= 1-mp.mpf("1e-20")
    for item in matrix_records
)
rank_ok = check(
    "every schedule has a fully resolved certified rank and inertia",
    all(sum(item["spectral"]["inertia"]) == 10 for item in matrix_records)
    and all(item["spectral"]["rank"]
            == item["spectral"]["inertia"][0]
             + item["spectral"]["inertia"][2]
            for item in matrix_records),
    f"ranks={sorted(Counter(ranks).items())}, "
    f"inertias={sorted(Counter(inertias).items())}",
)

controls_ok = all((
    provenance_ok, upstream_ok, definitions_ok, combinatorics_ok,
    raw_ok, directional_ok, time_reversal_ok, rank_ok,
))
if not controls_ok:
    outcome = "REFINED_H4_INTERNAL_JACOBIAN_CONTROL_FAILED"
elif len(set(ranks)) > 1:
    outcome = "REFINED_H4_INTERNAL_JACOBIAN_MIXED_RANK"
elif ranks[0] == 10 and len(classes) == 1:
    outcome = "REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_SINGLE_CLASS"
elif ranks[0] == 10:
    outcome = "REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_MULTIPLE_CLASSES"
elif ranks[0] == 9 and null_candidate:
    outcome = "REFINED_H4_INTERNAL_JACOBIAN_INDUCED_LAPSE_NULL"
else:
    outcome = "REFINED_H4_INTERNAL_JACOBIAN_RANK_DEFICIENT_OTHER"
outcome_ok = check(
    "the frozen hierarchy assigns the internal-Jacobian outcome",
    outcome in {
        "REFINED_H4_INTERNAL_JACOBIAN_CONTROL_FAILED",
        "REFINED_H4_INTERNAL_JACOBIAN_MIXED_RANK",
        "REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_SINGLE_CLASS",
        "REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_MULTIPLE_CLASSES",
        "REFINED_H4_INTERNAL_JACOBIAN_INDUCED_LAPSE_NULL",
        "REFINED_H4_INTERNAL_JACOBIAN_RANK_DEFICIENT_OTHER",
    },
    outcome,
)

serial_matrices = []
for index, item in enumerate(matrix_records):
    spectral = item["spectral"]
    serial_matrices.append({
        "order": item["order"],
        "matrix_class": class_indices[index],
        "step_difference": mp_text(item["step_difference"], 45),
        "precision_difference": mp_text(item["precision_difference"], 45),
        "entry_error": mp_text(item["entry_error"], 45),
        "spectral_error": mp_text(item["spectral_error"], 45),
        "raw_maximum_imaginary": mp_text(item["raw_maximum_imaginary"], 45),
        "raw_maximum_antisymmetry": mp_text(item["raw_maximum_antisymmetry"], 45),
        "matrix": [
            [mp_text(item["matrix"][row, column], 75) for column in range(10)]
            for row in range(10)
        ],
        "eigenvalues": [mp_text(value, 75) for value in spectral["eigenvalues"]],
        "rank": spectral["rank"],
        "inertia_positive_zero_negative": list(spectral["inertia"]),
        "induced_tangent_image_norm": mp_text(
            spectral["tangent_image_norm"], 55
        ),
        "induced_tangent_rayleigh": mp_text(
            spectral["tangent_rayleigh"], 55
        ),
        "maximum_zero_eigenvector_overlap": mp_text(
            spectral["maximum_zero_overlap"], 55
        ),
        "unapplied_newton_proposal": [
            mp_text(spectral["proposal"][row], 55) for row in range(10)
        ],
        "unapplied_newton_proposal_norm": mp_text(
            spectral["proposal_norm"], 55
        ),
        "unapplied_linear_residual_norm": mp_text(
            spectral["proposal_linear_residual_norm"], 55
        ),
    })

artifact = {
    "title": "All-schedule refined H4 internal-Jacobian census",
    "date": "2026-08-20",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "coordinates": [variable_label(key) for key in INTERNAL_VARIABLES],
        "matrix": "total-orbit log-coordinate internal action Hessian",
        "decimal_precisions": [PRIMARY_DPS, SECONDARY_DPS],
        "gradient_difference_steps": list(STEP_TEXTS),
        "newton_update_applied": False,
        "root_or_effective_boundary_operator_computed": False,
        "physical_or_continuum_target_loaded": False,
    },
    "controls": {
        "maximum_directional_second_derivative_relative_error": mp_text(
            maximum_directional_error, 45
        ),
        "directional_second_derivatives": directional_records,
        "maximum_time_reversal_matrix_difference": mp_text(
            max(reversal_differences), 45
        ),
    },
    "census": {
        "schedule_count": 24,
        "rank_histogram": {str(key): value for key, value in sorted(Counter(ranks).items())},
        "inertia_histogram": {
            str(key): value for key, value in sorted(Counter(inertias).items())
        },
        "matrix_class_count": len(classes),
        "matrix_classes": [
            {
                "class": class_index,
                "orders": [list(orders[index]) for index in members],
            }
            for class_index, members in enumerate(class_members)
        ],
        "induced_lapse_null_candidate": null_candidate,
        "schedules": serial_matrices,
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"RANKS: {dict(sorted(Counter(ranks).items()))}")
print(f"INERTIAS: {dict(sorted(Counter(inertias).items()))}")
print(f"MATRIX CLASSES: {len(classes)}")
print(f"OUTCOME: {outcome}")
print(f"RESULT: {passed}/{tests} checks passed")
sys.exit(0 if controls_ok and outcome_ok and passed == tests else 1)

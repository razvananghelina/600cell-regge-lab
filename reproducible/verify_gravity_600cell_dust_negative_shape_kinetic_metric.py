#!/usr/bin/env python3
"""Canonical kinetic-metric Rouché audit of the negative-shape recurrence."""

import ast
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "verify_gravity_600cell_dust_negative_shape_root_count.py"
SOURCE_RESULT = HERE / "gravity_600cell_dust_negative_shape_root_count.json"
CARRIERS = HERE / "gravity_600cell_dust_negative_shape_root_count_carriers.json"
OUTPUT = HERE / "gravity_600cell_dust_negative_shape_kinetic_metric.json"
PRIOR_ART_COMMIT = "cdcaf8d"
PROTOCOL_COMMIT = "afe2c4e"
EXPORT_COMMIT = "6bd1a02"
EXPECTED_HASHES = {
    "source": "f7878c25ecd49291a6d3387ddca1b94944f9416b63a90ae2656bf7fce2fc9cca",
    "source_result": "7f71d680f4ba34da2f3a8af7c8e4b92668e1cf16e682acbf83f9d6bbec41b280",
    "carriers": "a0586e7e4a4b61d511703988f37f7245b29a424692e2777d5ebe0ea8b077d8c3",
}
MACHINE_EPSILON = np.finfo(float).eps
INITIAL_INTERVALS = 256
MAX_DEPTH = 32
MAX_INTERVALS = 2_000_000
EXPECTED_CELLS = 16
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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decode_matrix(record):
    matrix = np.asarray(record["real"], dtype=float) + 1j * np.asarray(
        record["imag"], dtype=float
    )
    if list(matrix.shape) != record["shape"]:
        raise RuntimeError("carrier matrix shape mismatch")
    return matrix


def load_frozen_contour_helpers():
    wanted = {"sf", "operator_norm", "quadratic_value", "contour_cover"}
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    body = [node for node in tree.body if isinstance(node, ast.FunctionDef)
            and node.name in wanted]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"frozen contour helper mismatch: {wanted-found}")
    exec(compile(ast.Module(body=body, type_ignores=[]), str(SOURCE), "exec"),
         globals())


def fixed_positive_sqrt(metric):
    values, vectors = la.eigh(metric)
    if values[0] <= 0:
        return values, None, None
    root = (vectors * np.sqrt(values)) @ vectors.conj().T
    inverse = (vectors * (1 / np.sqrt(values))) @ vectors.conj().T
    return values, root, inverse


def transform_coefficients(root, inverse, coefficients):
    return tuple(root @ matrix @ inverse for matrix in coefficients)


def unitary_reversal_control(coefficients):
    m = coefficients[0].shape[0]
    reversal = np.fliplr(np.eye(m, dtype=np.complex128))
    transformed = tuple(reversal @ matrix @ reversal for matrix in coefficients)
    maximum = 0.0
    maximum_floor = 0.0
    for j in range(INITIAL_INTERVALS):
        theta = 2 * math.pi * (j + 0.5) / INITIAL_INTERVALS
        z = complex(math.cos(theta), math.sin(theta))
        left = quadratic_value(*coefficients, z)
        right = quadratic_value(*transformed, z)
        left_values = la.svdvals(left)
        right_values = la.svdvals(right)
        difference = float(np.max(np.abs(left_values - right_values)))
        floor = float(
            1000 * MACHINE_EPSILON * m
            * max(1.0, operator_norm(left), operator_norm(right))
        )
        maximum = max(maximum, difference)
        maximum_floor = max(maximum_floor, floor)
    return maximum, maximum_floor, bool(maximum <= 100 * maximum_floor)


print("=" * 78)
print("600-CELL NEGATIVE SHAPE KINETIC-METRIC ROOT CERTIFICATE")
print("=" * 78)

hashes = {
    "source": sha256(SOURCE),
    "source_result": sha256(SOURCE_RESULT),
    "carriers": sha256(CARRIERS),
}
source_result = json.loads(SOURCE_RESULT.read_text())
carrier_source = json.loads(CARRIERS.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and source_result["outcome"] == "NEGATIVE_SHAPE_ROOT_COUNT_SAFETY_OPEN"
    and source_result["passed"] == source_result["tests"] == 13
    and source_result["counts"]["inherited_cells"] == EXPECTED_CELLS
    and source_result["counts"]["literal_cover"] == {"PASS": EXPECTED_CELLS}
    and source_result["counts"]["safety_100_cover"] == {"OPEN": EXPECTED_CELLS}
    and len(carrier_source["cells"]) == EXPECTED_CELLS
)
check("the target-disclosed carrier export has exact frozen provenance", provenance_ok)

prior_counts = {}
for parity, rows in source_result["parities"].items():
    for row in rows:
        for variant in row["variants"]:
            key = (parity, row["sector_index"], variant["variant"])
            prior_counts[key] = variant["transferred_counts"]
prior_count_ok = bool(
    len(prior_counts) == EXPECTED_CELLS
    and all(value == {"inside": 15, "on": 0, "outside": 15}
            for value in prior_counts.values())
)
check("all inherited similarity-invariant counts are exactly 15/0/15", prior_count_ok)

load_frozen_contour_helpers()
records = []
kinetic_counts = Counter()
literal_counts = Counter()
safety_counts = Counter()
sqrt_failures = scalar_failures = unitary_failures = transport_failures = 0
identifier_counts = Counter()

for cell in carrier_source["cells"]:
    key = (cell["parity"], cell["sector_index"], cell["variant"])
    identifier_counts[key] += 1
    a2, a1, a0 = (decode_matrix(cell[name]) for name in ("A2", "A1", "A0"))
    metric = decode_matrix(cell["B_negative"])
    epsilon_g = float(cell["epsilon_Gamma"])
    epsilon_o = float(cell["epsilon_Omega"])
    epsilon_b = float(cell["epsilon_B_negative"])
    m = metric.shape[0]
    values_b, root, inverse = fixed_positive_sqrt(metric)
    minimum_b = float(values_b[0])
    maximum_b = float(values_b[-1])
    positive_resolved = bool(minimum_b > 100 * epsilon_b)
    kinetic_counts["POSITIVE_RESOLVED" if positive_resolved else "OPEN"] += 1

    record = {
        "parity": key[0], "sector_index": key[1], "variant": key[2],
        "minimum_B_eigenvalue": sf(minimum_b),
        "maximum_B_eigenvalue": sf(maximum_b),
        "epsilon_B": sf(epsilon_b),
        "B_error_units": sf(minimum_b / epsilon_b if epsilon_b else math.inf),
        "kinetic_label": "POSITIVE_RESOLVED" if positive_resolved else "OPEN",
        "inherited_counts": prior_counts.get(key),
    }
    if not positive_resolved or root is None:
        record["status"] = "KINETIC_OPEN"
        records.append(record)
        continue

    identity = np.eye(m, dtype=np.complex128)
    norm_b = operator_norm(metric)
    norm_root = operator_norm(root)
    reconstruction_floor = float(
        1000 * MACHINE_EPSILON * m
        * max(1.0, norm_b, norm_root * norm_root)
    )
    root_residual = operator_norm(root @ root - metric)
    inverse_residual = operator_norm(root @ inverse - identity)
    sqrt_ok = bool(
        root_residual <= reconstruction_floor
        and inverse_residual <= reconstruction_floor
    )
    sqrt_failures += int(not sqrt_ok)
    condition_b = float(maximum_b / minimum_b)
    condition_s = float(math.sqrt(condition_b))
    transformed = transform_coefficients(root, inverse, (a2, a1, a0))
    epsilon_g_b = condition_s * epsilon_g
    epsilon_o_b = condition_s * epsilon_o
    transport_ok = bool(
        abs(epsilon_g_b - condition_s * epsilon_g) == 0
        and abs(epsilon_o_b - condition_s * epsilon_o) == 0
    )
    transport_failures += int(not transport_ok)

    literal = contour_cover(*transformed, epsilon_g_b, epsilon_o_b, 1)
    safety = contour_cover(*transformed, epsilon_g_b, epsilon_o_b, 100)
    literal_counts["PASS" if literal["passed"] else "OPEN"] += 1
    safety_counts["PASS" if safety["passed"] else "OPEN"] += 1

    scaled_values, scaled_root, scaled_inverse = fixed_positive_sqrt(7 * metric)
    scaled = transform_coefficients(scaled_root, scaled_inverse, (a2, a1, a0))
    scaled_condition = float(math.sqrt(scaled_values[-1] / scaled_values[0]))
    transform_floor = float(
        1000 * MACHINE_EPSILON * m
        * max(1.0, *(operator_norm(matrix) for matrix in transformed))
    )
    scalar_matrix_difference = max(
        operator_norm(left - right) for left, right in zip(transformed, scaled)
    )
    scalar_condition_difference = abs(scaled_condition - condition_s)
    condition_floor = float(
        1000 * MACHINE_EPSILON * m * max(1.0, condition_s)
    )
    scaled_literal = contour_cover(
        *scaled, scaled_condition * epsilon_g, scaled_condition * epsilon_o, 1
    )
    scaled_safety = contour_cover(
        *scaled, scaled_condition * epsilon_g, scaled_condition * epsilon_o, 100
    )
    scalar_ok = bool(
        scalar_matrix_difference <= 100 * transform_floor
        and scalar_condition_difference <= 100 * condition_floor
        and scaled_literal["passed"] == literal["passed"]
        and scaled_safety["passed"] == safety["passed"]
    )
    scalar_failures += int(not scalar_ok)

    unitary_difference, unitary_floor, unitary_ok = unitary_reversal_control(
        (a2, a1, a0)
    )
    unitary_failures += int(not unitary_ok)
    record.update({
        "status": "AUDITED",
        "root_residual": sf(root_residual),
        "inverse_residual": sf(inverse_residual),
        "reconstruction_floor": sf(reconstruction_floor),
        "square_root_control": sqrt_ok,
        "condition_B": sf(condition_b),
        "condition_S": sf(condition_s),
        "epsilon_Gamma_B": sf(epsilon_g_b),
        "epsilon_Omega_B": sf(epsilon_o_b),
        "transport_control": transport_ok,
        "literal_rouche_cover": literal,
        "safety_100_cover": safety,
        "scalar_rescaling_control": {
            "matrix_difference": sf(scalar_matrix_difference),
            "matrix_floor": sf(transform_floor),
            "condition_difference": sf(scalar_condition_difference),
            "condition_floor": sf(condition_floor),
            "literal_verdict_equal": scaled_literal["passed"] == literal["passed"],
            "safety_verdict_equal": scaled_safety["passed"] == safety["passed"],
            "passed": scalar_ok,
        },
        "unitary_reversal_control": {
            "maximum_singular_value_difference": sf(unitary_difference),
            "maximum_floor": sf(unitary_floor),
            "passed": unitary_ok,
        },
    })
    records.append(record)

identifiers_ok = bool(
    len(identifier_counts) == EXPECTED_CELLS
    and all(count == 1 for count in identifier_counts.values())
    and set(identifier_counts) == set(prior_counts)
)
check("all 16 carrier identifiers occur exactly once", identifiers_ok)
check(
    "all kinetic metrics are completely classified before transformation",
    sum(kinetic_counts.values()) == EXPECTED_CELLS,
    str(dict(kinetic_counts)),
)
check(
    "all constructed positive square roots and inverses satisfy their floors",
    sqrt_failures == 0,
    f"failures={sqrt_failures}",
)
check(
    "all Euclidean coefficient balls include the mandatory kinetic condition factor",
    transport_failures == 0,
    f"failures={transport_failures}",
)
check(
    "multiplying the kinetic metric by seven changes no transformed verdict",
    scalar_failures == 0,
    f"failures={scalar_failures}",
)
check(
    "the fixed unitary reversal preserves all sampled singular values",
    unitary_failures == 0,
    f"failures={unitary_failures}",
)
check(
    "all constructed literal and 100x covers are classified",
    sum(literal_counts.values()) == kinetic_counts["POSITIVE_RESOLVED"]
    and sum(safety_counts.values()) == kinetic_counts["POSITIVE_RESOLVED"],
    f"literal={dict(literal_counts)}, safety={dict(safety_counts)}",
)

controls_ok = bool(
    provenance_ok and prior_count_ok and identifiers_ok
    and sqrt_failures == scalar_failures == unitary_failures == transport_failures == 0
)
if not controls_ok:
    outcome = "NEGATIVE_SHAPE_KINETIC_METRIC_CONTROL_FAILED"
elif kinetic_counts["OPEN"]:
    outcome = "NEGATIVE_SHAPE_KINETIC_METRIC_OPEN"
elif literal_counts["OPEN"]:
    outcome = "NEGATIVE_SHAPE_KINETIC_LITERAL_OPEN"
elif len(safety_counts) > 1:
    outcome = "NEGATIVE_SHAPE_KINETIC_SCHEDULE_DEPENDENT"
elif safety_counts["OPEN"] == EXPECTED_CELLS:
    outcome = "NEGATIVE_SHAPE_KINETIC_SAFETY_OPEN"
elif safety_counts["PASS"] == EXPECTED_CELLS:
    outcome = "NEGATIVE_SHAPE_LOCAL_HYPERBOLIC_KINETIC_RESOLVED"
else:
    outcome = "NEGATIVE_SHAPE_KINETIC_MIXED_OPEN"

allowed = {
    "NEGATIVE_SHAPE_KINETIC_METRIC_CONTROL_FAILED",
    "NEGATIVE_SHAPE_KINETIC_METRIC_OPEN",
    "NEGATIVE_SHAPE_KINETIC_LITERAL_OPEN",
    "NEGATIVE_SHAPE_KINETIC_SCHEDULE_DEPENDENT",
    "NEGATIVE_SHAPE_KINETIC_SAFETY_OPEN",
    "NEGATIVE_SHAPE_LOCAL_HYPERBOLIC_KINETIC_RESOLVED",
    "NEGATIVE_SHAPE_KINETIC_MIXED_OPEN",
}
check("the preregistered kinetic-metric outcome tree is exhausted",
      outcome in allowed, outcome)

payload = {
    "title": "Canonical kinetic-metric Rouche certificate",
    "date": "2026-08-18",
    "classification": "DERIVED COMPUTATIONAL, TARGET-DISCLOSED",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "carrier_export_commit": EXPORT_COMMIT,
    "input_hashes": hashes,
    "protocol": {
        "metric": "unique positive square root of B_negative",
        "error_transport": "kappa_2(S) times Euclidean balls",
        "scalar_rescaling_control": 7,
        "unitary_control": "coordinate reversal",
        "interpretation_scope": "local frozen recurrence only",
    },
    "counts": {
        "kinetic": dict(kinetic_counts),
        "literal_cover": dict(literal_counts),
        "safety_100_cover": dict(safety_counts),
        "square_root_failures": sqrt_failures,
        "scalar_rescaling_failures": scalar_failures,
        "unitary_failures": unitary_failures,
        "transport_failures": transport_failures,
    },
    "cells": records,
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print()
print("kinetic labels:", dict(kinetic_counts))
print("literal covers:", dict(literal_counts))
print("100x covers:", dict(safety_counts))
print("outcome:", outcome)
print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)

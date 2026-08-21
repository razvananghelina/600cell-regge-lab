#!/usr/bin/env python3
"""Corrected adjudication of the frozen direct-action constrained H4 audit."""

from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FAILED_VERIFIER = (
    HERE / "verify_gravity_600cell_refined_h4_constrained_response_adversarial.py"
)
DIRECT = HERE / "gravity_600cell_refined_h4_constrained_response_adversarial.json"
FIRST_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_adversarial_first_result.md"
)
AUXILIARY_VERIFIER = (
    HERE / "verify_gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic.py"
)
AUXILIARY = (
    HERE / "gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic.json"
)
AUXILIARY_PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic_protocol.md"
)
AUXILIARY_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic_result.md"
)
PRIMARY = HERE / "gravity_600cell_refined_h4_constrained_response_corrected.json"
NULL_COUPLING = HERE / "gravity_600cell_refined_h4_null_coupling_adversarial.json"
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_corrected_adjudication_protocol.md"
)
RERUN_RECEIPT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_unchanged_rerun_receipt.md"
)
OUTPUT = (
    HERE / "gravity_600cell_refined_h4_constrained_response_corrected_adjudication.json"
)

PROTOCOL_COMMIT = "2a10336"
RERUN_COMMIT = "ab81207"
EXPECTED_HASHES = {
    "failed_verifier": "78f6b52f6f019a150a86ddadcb819b67c3757244c015687ab67f4649784ac53d",
    "direct": "a23ef4cc23d08ad8768f1df66789aa900cdb95a7f3529486df80697a53b1fe81",
    "first_result": "c4203c07b859ed323ee5049875d54d4894a1b815c3792b7c1d1de0e71677ad64",
    "auxiliary_verifier": "70beeffe19cf4b6e90a613d3936f9c30bd98021e0a7b6ae6b7e93d60c01c0bc4",
    "auxiliary": "f66177326afc3b3457a60b544745b739cbaa6b6d6e7f367b57d60f31eeeddeb7",
    "auxiliary_protocol": "2f6d9d72e04c4baf1dc385425ef7f26ba0a55f6249d6505db111aa21e0836405",
    "auxiliary_result": "afd3bb5bcec476bded8ea003c5749a83fa46c488b0d1c7da3d693da93fbe9423",
    "primary": "85adea23f6a19153f61f3ed066137a5e40ab77b8901d4cc81cfc4f864e0bc093",
    "null_coupling": "5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9",
    "protocol": "a4bf483b025abedcdcb17be7681dfbf5c19f590587a0d4a527f3ca745d5cf017",
    "rerun_receipt": "daa23f6b4cac09118a08768718996a92f77ad4a53d2cb86c81bce6fd2b8446f4",
}

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=70):
    return mp.nstr(value, digits, strip_zeros=False)


def matrix_from_text(rows):
    return mp.matrix([[mp.mpf(value) for value in row] for row in rows])


def matrix_max(matrix):
    if not matrix.rows or not matrix.cols:
        return mp.mpf(0)
    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows)
        for column in range(matrix.cols)
    )


def matrix_difference(left, right):
    return matrix_max(left - right)


def submatrix(matrix, rows, columns):
    return mp.matrix([
        [matrix[row, column] for column in columns]
        for row in rows
    ])


def kernel_basis(vector, pivot):
    columns = tuple(index for index in range(len(vector)) if index != pivot)
    if vector[pivot] == 0:
        raise ZeroDivisionError("kernel-basis pivot is zero")
    result = mp.matrix(len(vector), len(vector) - 1)
    for column, index in enumerate(columns):
        result[index, column] = 1
        result[pivot, column] = -vector[index] / vector[pivot]
    return result, columns


def reversal_matrix():
    result = mp.matrix(12, 12)
    for index in range(6):
        result[index, 6 + index] = 1
        result[6 + index, index] = 1
    return result


def transformed_envelope(matrix, envelope):
    maximum_column_sum = max(
        mp.fsum(abs(matrix[row, column]) for row in range(matrix.rows))
        for column in range(matrix.cols)
    )
    return maximum_column_sum ** 2 * envelope


def classify(matrices, envelopes):
    representatives = []
    memberships = []
    assignments = []
    for index, matrix in enumerate(matrices):
        assigned = None
        for class_index, representative in enumerate(representatives):
            if matrix_difference(matrix, matrices[representative]) <= (
                envelopes[index] + envelopes[representative]
            ):
                assigned = class_index
                break
        if assigned is None:
            assigned = len(representatives)
            representatives.append(index)
            memberships.append([])
        memberships[assigned].append(index)
        assignments.append(assigned)
    return representatives, memberships, assignments


print("=" * 78)
print("CORRECTED ADJUDICATION: DIRECT CONSTRAINED H4 RESPONSE")
print("=" * 78)

paths = {
    "failed_verifier": FAILED_VERIFIER,
    "direct": DIRECT,
    "first_result": FIRST_RESULT,
    "auxiliary_verifier": AUXILIARY_VERIFIER,
    "auxiliary": AUXILIARY,
    "auxiliary_protocol": AUXILIARY_PROTOCOL,
    "auxiliary_result": AUXILIARY_RESULT,
    "primary": PRIMARY,
    "null_coupling": NULL_COUPLING,
    "protocol": PROTOCOL,
    "rerun_receipt": RERUN_RECEIPT,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "the twice-reproduced direct data, diagnostic and corrected protocol have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PROTOCOL_COMMIT == "2a10336"
    and RERUN_COMMIT == "ab81207",
)

mp.mp.dps = 100
direct = json.loads(DIRECT.read_text())
auxiliary = json.loads(AUXILIARY.read_text())
null_coupling = json.loads(NULL_COUPLING.read_text())
historical_ok = check(
    "the historical audit remains the exact reproducible 15/17 control failure",
    direct["outcome"]
        == "ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED"
    and direct["tests"] == {"passed": 15, "total": 17}
    and direct["definitions"]["analytic_gradient_or_hessian_used"] is False
    and direct["definitions"]["primary_loaded_after_direct_census"] is True
    and direct["definitions"]["directions_per_complete_reconstruction"] == 210
    and direct["definitions"]["schedule_count"] == 24,
)

stationarity = auxiliary["stationarity"]
off_shell = auxiliary["off_shell_curvature"]
halving_ratios = tuple(mp.mpf(value) for value in off_shell["halving_ratios"])
diagnostic_ok = check(
    "the auxiliary diagnostic mechanically resolves both historical controls",
    auxiliary["outcome"]
        == "REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_FAILURES_RESOLVED"
    and auxiliary["tests"] == {"passed": 13, "total": 13}
    and stationarity["component_count_per_precision"] == 240
    and stationarity["all_zero_gates_pass"] is True
    and stationarity["all_precision_gates_pass"] is True
    and mp.mpf(stationarity["maximum_tenth_order_residual"]) < mp.mpf("1e-60")
    and off_shell["directions_per_schedule"] == 210
    and off_shell["resolved_parity_count"] == 5040
    and all(mp.mpf("1.99") <= value <= mp.mpf("2.01") for value in halving_ratios)
    and mp.mpf(off_shell["maximum_even_to_odd_ratio"]) < mp.mpf("1e-6")
    and mp.mpf(off_shell["maximum_relative_imaginary_complete_action"])
        < mp.mpf("1e-100")
    and mp.mpf(off_shell["maximum_angle_identity_residual"]) < mp.mpf("1e-100")
    and mp.mpf(off_shell["minimum_angle_argument"]) > mp.mpf("1e-2"),
)
failure_link_ok = check(
    "the diagnostic reproduces the failed derivative and off-shell curvature scales",
    abs(
        mp.mpf(stationarity["maximum_legacy_fourth_order_residual"])
        - mp.mpf(direct["controls"]["maximum_stationarity_residual"])
    ) < mp.mpf("1e-100")
    and abs(
        mp.mpf(off_shell["step_maximum_imaginary_curvatures"][0])
        - mp.mpf(direct["branch"]["maximum_imaginary_curvature"])
    ) < mp.mpf("1e-70"),
)

schedules = direct["census"]["schedules"]
orders = tuple(permutations(range(4)))
dimension_ok = check(
    "all frozen direct second variations, lifts, responses and spectra are complete",
    len(schedules) == 24
    and tuple(tuple(item["order"]) for item in schedules) == orders
    and all(
        len(item["restricted_second_variation"]) == 20
        and all(len(row) == 20 for row in item["restricted_second_variation"])
        and len(item["lift_coefficients"]) == 9
        and all(len(row) == 11 for row in item["lift_coefficients"])
        and len(item["response"]) == 11
        and all(len(row) == 11 for row in item["response"])
        and len(item["internal_eigenvalues"]) == 9
        for item in schedules
    ),
)

direct_numerics_ok = check(
    "all stored direct internal, solve, precision and polynomial controls pass",
    all(
        min(mp.mpf(value) for value in item["internal_eigenvalues"])
            > 100 * mp.mpf(item["second_variation_envelope"])
        and mp.mpf(item["solve_residual"]) < mp.mpf("1e-80")
        and mp.mpf(item["maximum_coordinate_displacement"]) < mp.mpf("2e-5")
        and mp.mpf(item["maximum_ladder_imaginary"])
            <= mp.mpf(item["second_variation_envelope"])
        for item in schedules
    )
    and all(
        mp.mpf(record["difference"]) <= mp.mpf(record["gate"])
        for record in direct["controls"]["precision_repeat"].values()
    )
    and mp.mpf(direct["controls"]["polynomial_error"])
        <= mp.mpf(direct["controls"]["polynomial_envelope"])
    and mp.mpf(direct["controls"]["wrong_polarization_error"])
        > mp.mpf("1e6") * mp.mpf(direct["controls"]["polynomial_envelope"]),
)

c = mp.matrix([
    mp.mpf(value) for value in null_coupling["compatibility"]["adversarial_row"]
])
p, p_columns = kernel_basis(c, 3)
reversal = reversal_matrix()
rp = reversal * p
t_reversal = submatrix(rp, p_columns, tuple(range(11)))
basis_error = max(matrix_max(c.T * p), matrix_difference(p * t_reversal, rp))
basis_ok = check(
    "the reversal action is reconstructed algebraically from the accepted compatibility row",
    basis_error < mp.mpf("1e-90"),
    f"max identity error={mp_text(basis_error, 8)}",
)

responses = [matrix_from_text(item["response"]) for item in schedules]
response_envelopes = [mp.mpf(item["response_envelope"]) for item in schedules]
order_index = {order: index for index, order in enumerate(orders)}
reversal_differences = []
reversal_gates = []
canonical_matrices = []
canonical_envelopes = []
for index, order in enumerate(orders):
    reverse_index = order_index[tuple(reversed(order))]
    transformed_reverse = t_reversal.T * responses[reverse_index] * t_reversal
    reversal_differences.append(matrix_difference(responses[index], transformed_reverse))
    reversal_gates.append(
        response_envelopes[index]
        + transformed_envelope(t_reversal, response_envelopes[reverse_index])
    )
    if order <= tuple(reversed(order)):
        canonical_matrices.append(responses[index])
        canonical_envelopes.append(response_envelopes[index])
    else:
        canonical_matrices.append(t_reversal.T * responses[index] * t_reversal)
        canonical_envelopes.append(
            transformed_envelope(t_reversal, response_envelopes[index])
        )

time_reversal_covariant = all(
    difference <= gate
    for difference, gate in zip(reversal_differences, reversal_gates)
)
representatives, memberships, assignments = classify(
    canonical_matrices, canonical_envelopes
)
print("[INFO] independently re-adjudicated direct classes:", flush=True)
for class_index, members in enumerate(memberships):
    print(
        f"[INFO] class {class_index}: "
        + ", ".join(str(orders[index]) for index in members),
        flush=True,
    )
class_census_ok = check(
    "time reversal and the target-free direct class census are recomputed without stored labels",
    len(assignments) == 24 and len(memberships) >= 1,
    f"reversal={time_reversal_covariant}, classes={len(memberships)}",
)

# Delayed target comparison: primary matrices are loaded only after the direct
# class census above has been rebuilt.
primary = json.loads(PRIMARY.read_text())
primary_ok = check(
    "the delayed primary target is the frozen clean one-class result",
    primary["outcome"] == "REFINED_H4_CONSTRAINED_RESPONSE_SINGLE_SCHEDULE_CLASS"
    and primary["tests"] == {"passed": 19, "total": 19}
    and primary["census"]["schedule_count"] == 24,
)
primary_differences = []
primary_gates = []
for index in range(24):
    primary_record = primary["census"]["schedules"][index]
    primary_matrix = matrix_from_text(primary_record["primary_response"])
    primary_envelope = mp.mpf(primary_record["response_envelopes"]["primary"])
    primary_differences.append(matrix_difference(responses[index], primary_matrix))
    primary_gates.append(10 * (response_envelopes[index] + primary_envelope))
primary_matches = [
    difference <= gate
    for difference, gate in zip(primary_differences, primary_gates)
]
primary_comparison_ok = check(
    "all 24 direct matrices independently re-match the analytic-gradient primary route",
    all(primary_matches),
    f"matches={sum(primary_matches)}/24, max fraction={mp_text(max(
        difference / gate for difference, gate in zip(primary_differences, primary_gates)
    ), 8)}",
)

corruption_size = mp.mpf("1e-6") * max(mp.mpf(1), matrix_max(canonical_matrices[0]))
corrupted = mp.matrix(canonical_matrices[0])
corrupted[0, 0] += corruption_size
corrupt_representatives, corrupt_memberships, _ = classify(
    [canonical_matrices[0], corrupted],
    [canonical_envelopes[0], canonical_envelopes[0]],
)
corrupt_primary_difference = matrix_difference(
    corrupted,
    matrix_from_text(primary["census"]["schedules"][0]["primary_response"]),
)
corruption_ok = check(
    "a freshly reconstructed one-entry corruption splits the class and fails primary comparison",
    len(corrupt_representatives) == len(corrupt_memberships) == 2
    and corrupt_primary_difference > primary_gates[0]
    and abs(corruption_size - mp.mpf(direct["controls"]["corruption_size"]))
        < mp.mpf("1e-55"),
    f"corruption/gate={mp_text(corrupt_primary_difference / primary_gates[0], 8)}",
)

controls_ok = all((
    provenance_ok, historical_ok, diagnostic_ok, failure_link_ok, dimension_ok,
    direct_numerics_ok, basis_ok, class_census_ok, primary_ok, corruption_ok,
))
if not controls_ok:
    outcome = "CORRECTED_ADJUDICATION_REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED"
elif (
    not time_reversal_covariant
    or len(memberships) != 1
    or not primary_comparison_ok
):
    outcome = "CORRECTED_ADJUDICATION_REFINED_H4_CONSTRAINED_RESPONSE_DISAGREEMENT"
else:
    outcome = "CORRECTED_ADJUDICATION_REFINED_H4_CONSTRAINED_RESPONSE_CORROBORATED"

expected_corroboration = (
    controls_ok
    and time_reversal_covariant
    and len(memberships) == 1
    and primary_comparison_ok
)
outcome_ok = check(
    "the corrected outcome follows mechanically from frozen controls and recomputed comparisons",
    (outcome.endswith("_CORROBORATED")) == expected_corroboration
    and (outcome.endswith("_DISAGREEMENT")) == (
        controls_ok and not expected_corroboration
    )
    and (outcome.endswith("_CONTROL_FAILED")) == (not controls_ok),
    outcome,
)

artifact = {
    "title": "Corrected adjudication of direct constrained H4 response",
    "date": "2026-08-21",
    "protocol_commit": PROTOCOL_COMMIT,
    "unchanged_rerun_commit": RERUN_COMMIT,
    "input_sha256": actual_hashes,
    "historical_audit": {
        "outcome_preserved": direct["outcome"],
        "tests_preserved": direct["tests"],
        "byte_identical_rerun": True,
        "direct_artifact_sha256": actual_hashes["direct"],
    },
    "auxiliary_resolution": {
        "outcome": auxiliary["outcome"],
        "tests": auxiliary["tests"],
        "maximum_tenth_order_stationarity_residual": (
            stationarity["maximum_tenth_order_residual"]
        ),
        "curvature_halving_ratios": off_shell["halving_ratios"],
        "resolved_curvature_parity_count": off_shell["resolved_parity_count"],
        "maximum_relative_imaginary_complete_action": (
            off_shell["maximum_relative_imaginary_complete_action"]
        ),
    },
    "recomputed_census": {
        "schedule_count": 24,
        "class_count": len(memberships),
        "classes": [
            {
                "indices": members,
                "orders": [orders[index] for index in members],
            }
            for members in memberships
        ],
        "assignments": assignments,
        "time_reversal_covariant": time_reversal_covariant,
        "maximum_reversal_difference": mp_text(max(reversal_differences)),
        "maximum_reversal_fraction": mp_text(max(
            difference / gate
            for difference, gate in zip(reversal_differences, reversal_gates)
        )),
        "primary_match_count": sum(primary_matches),
        "maximum_primary_difference": mp_text(max(primary_differences)),
        "maximum_primary_fraction": mp_text(max(
            difference / gate
            for difference, gate in zip(primary_differences, primary_gates)
        )),
        "minimum_internal_eigenvalue": mp_text(min(
            mp.mpf(value)
            for item in schedules for value in item["internal_eigenvalues"]
        )),
        "maximum_solve_residual": mp_text(max(
            mp.mpf(item["solve_residual"]) for item in schedules
        )),
    },
    "controls": {
        "basis_identity_error": mp_text(basis_error),
        "precision_repeat_count": len(direct["controls"]["precision_repeat"]),
        "corruption_class_count": len(corrupt_memberships),
        "corruption_primary_fraction": mp_text(
            corrupt_primary_difference / primary_gates[0]
        ),
    },
    "scope": {
        "new_action_or_response_matrix_evaluated": False,
        "stored_class_labels_trusted": False,
        "nonhomogeneous_operator_or_spectrum_computed": False,
        "full_suite_executed": False,
        "root_search_or_deferred_census_executed": False,
        "physical_constant_extracted": False,
    },
    "status_labels": {
        "homogeneous_H4_schedule_independence": (
            "DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED"
            if outcome.endswith("_CORROBORATED") else "OPEN"
        ),
        "nonhomogeneous_propagation": "OPEN",
        "tick_c_G_planck_particles": "NOT ESTABLISHED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Recomputed classes: {len(memberships)}")
print(f"Primary matches: {sum(primary_matches)}/24")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")

if passed != tests:
    raise SystemExit(1)


#!/usr/bin/env python3
"""Post-commit comparison of canonical-carrier acceleration with FLRW."""

from hashlib import sha256
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
INPUT = HERE / "gravity_600cell_projected_rank_edgewise_acceleration_blind.json"
OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_acceleration_comparison.json"
INPUT_SHA256 = "2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2"
BLIND_COMMIT = "9469e33"
TARGET = -0.5
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return condition


print("="*78)
print("CANONICAL PROJECTED-CARRIER ACCELERATION COMPARISON")
print("="*78)

actual_hash = sha256(INPUT.read_bytes()).hexdigest()
blind = json.loads(INPUT.read_text())
check(
    "the target-blind artifact was committed before comparison",
    actual_hash == INPUT_SHA256 and BLIND_COMMIT == "9469e33",
    f"input={actual_hash}, commit={BLIND_COMMIT}",
)
check(
    "Stage A loaded neither the continuum nor projected-red targets",
    blind["continuum_target_loaded"] is False
    and blind["projected_red_coefficients_loaded"] is False
    and blind["outcome"] ==
    "CANONICAL_CARRIER_ACCELERATION_COEFFICIENTS_DERIVED"
    and blind["passed"] == blind["tests"] == 10,
)

coefficients = blind["blind_coefficients"]
expected_names = ("projected_barycentric", "projected_rank_edgewise_2")
check(
    "the comparison has exactly the two preregistered coefficients",
    tuple(sorted(coefficients)) == tuple(sorted(expected_names)),
    str(coefficients),
)

errors = {name: abs(coefficients[name]-TARGET) for name in expected_names}
improvement = errors[expected_names[0]]/errors[expected_names[1]]
calibration_pass = errors[expected_names[1]] < errors[expected_names[0]]
check(
    "the canonical refinement strictly improves toward closed-dust FLRW",
    calibration_pass,
    f"errors={errors}, reduction={improvement:.12g}",
)

outcome = (
    "CANONICAL_CARRIER_HOMOGENEOUS_CALIBRATION_PASS"
    if calibration_pass else
    "CANONICAL_CARRIER_HOMOGENEOUS_CALIBRATION_FAILURE"
)
check(
    "the preregistered comparison hierarchy assigns one outcome",
    outcome in {
        "CANONICAL_CARRIER_HOMOGENEOUS_CALIBRATION_PASS",
        "CANONICAL_CARRIER_HOMOGENEOUS_CALIBRATION_FAILURE",
    },
    outcome,
)

artifact = {
    "blind_commit": BLIND_COMMIT,
    "blind_input_sha256": actual_hash,
    "continuum_closed_dust_half_step_coefficient": TARGET,
    "coefficients": coefficients,
    "absolute_errors": errors,
    "error_reduction_factor": improvement,
    "classification": "DERIVED_NUMERICAL_CONTROL",
    "asymptotic_order_status": "OPEN_ONE_REFINEMENT_ONLY",
    "inhomogeneous_dynamics_status": "OPEN",
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("="*78)
print(f"Errors: {errors}")
print(f"Reduction factor: {improvement:.12g}")
print(f"Outcome: {outcome}")
print(f"{passed}/{tests} checks passed")
sys.exit(0 if passed == tests and calibration_pass else 1)

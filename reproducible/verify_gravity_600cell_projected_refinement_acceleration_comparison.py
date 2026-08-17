#!/usr/bin/env python3
"""Disclosed continuum comparison of the frozen refined Regge coefficients."""

import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE / "gravity_600cell_projected_refinement_acceleration_blind.json"
OUTPUT = HERE / "gravity_600cell_projected_refinement_acceleration_comparison.json"
INPUT_SHA256 = (
    "640bc0dd3d6f1ae727f8113bf29514878874effffd14f539f5a43e3c3b18d069"
)
BLIND_RESULT_COMMIT = "98a1e1d"
TARGET = -0.5
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


artifact = json.loads(INPUT.read_text())
input_ok = bool(
    digest(INPUT) == INPUT_SHA256
    and artifact.get("outcome")
        == "PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_DERIVED"
    and artifact.get("tests") == artifact.get("passed") == 9
    and artifact.get("protocol")
        == "direct coefficient production; no refined continuum-target comparison performed"
    and BLIND_RESULT_COMMIT == "98a1e1d"
)
check(
    "the disclosed stage reads the exact committed blind artifact",
    input_ok,
    f"sha256={digest(INPUT)}",
)

variants = artifact["variants"]
coefficients = {
    variant: {
        level: float(artifact["records"][variant][f"level{level}"]
                     ["coefficient_audit"]["coefficient"])
        for level in range(3)
    }
    for variant in variants
}

comparisons = {}
finite_ok = True
for variant in variants:
    comparisons[variant] = {}
    for level in range(3):
        coefficient = coefficients[variant][level]
        ratio = coefficient/TARGET
        error = abs(ratio-1)
        finite_ok &= bool(np.isfinite(coefficient)
                          and np.isfinite(ratio) and np.isfinite(error))
        comparisons[variant][f"level{level}"] = {
            "coefficient": coefficient,
            "coefficient_over_friedmann": ratio,
            "relative_error": error,
        }
    error0 = comparisons[variant]["level0"]["relative_error"]
    error1 = comparisons[variant]["level1"]["relative_error"]
    error2 = comparisons[variant]["level2"]["relative_error"]
    comparisons[variant]["improvement_factors"] = {
        "level0_error_over_level1_error": error0/error1,
        "level1_error_over_level2_error": error1/error2,
    }
check("all disclosed ratios and errors are finite", finite_ok)

strict_improvement = {
    variant: bool(
        comparisons[variant]["level1"]["relative_error"]
            < comparisons[variant]["level0"]["relative_error"]
        and comparisons[variant]["level2"]["relative_error"]
            < comparisons[variant]["level1"]["relative_error"])
    for variant in variants
}
all_improve = all(strict_improvement.values())
check(
    "every registered regulator tower improves strictly at both levels",
    all_improve,
    str(strict_improvement),
)

spreads = {}
subdominant = {}
for level in (1, 2):
    values = np.asarray([coefficients[variant][level] for variant in variants])
    spread = float(np.max(values)-np.min(values))
    minimum_distance = float(np.min(abs(values-TARGET)))
    spreads[f"level{level}"] = {
        "coefficient_spread": spread,
        "minimum_distance_to_target": minimum_distance,
        "spread_over_minimum_distance": spread/minimum_distance,
    }
    subdominant[f"level{level}"] = spread < minimum_distance
regulator_robust = all(subdominant.values())
check(
    "the registered diagonal-choice spread is subdominant at both levels",
    regulator_robust,
    str({level: spreads[level]["spread_over_minimum_distance"]
         for level in spreads}),
)

if not input_ok or not finite_ok:
    outcome = "PROJECTED_REGGE_ACCELERATION_COMPARISON_OPEN"
elif not all_improve:
    outcome = "PROJECTED_REGGE_ACCELERATION_NOT_UNIVERSAL"
elif regulator_robust:
    outcome = "PROJECTED_REGGE_ACCELERATION_ROBUST_TWO_LEVEL_IMPROVEMENT"
else:
    outcome = "PROJECTED_REGGE_ACCELERATION_CHOICE_SENSITIVE_IMPROVEMENT"

expected_outcome = (
    "PROJECTED_REGGE_ACCELERATION_ROBUST_TWO_LEVEL_IMPROVEMENT"
    if input_ok and finite_ok and all_improve and regulator_robust
    else outcome
)
check(
    "the disclosed outcome follows the frozen hierarchy",
    outcome == expected_outcome and not outcome.endswith("OPEN"),
    f"outcome={outcome}",
)

payload = {
    "input": INPUT.name,
    "input_sha256": INPUT_SHA256,
    "blind_result_commit": BLIND_RESULT_COMMIT,
    "friedmann_coefficient": TARGET,
    "comparisons": comparisons,
    "strict_improvement": strict_improvement,
    "spreads": spreads,
    "subdominant_regulator_ambiguity": subdominant,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print("\nDisclosed continuum comparison:")
for variant in variants:
    row = comparisons[variant]
    print(
        f"  {variant}: "
        + ", ".join(
            f"level{level} ratio={row[f'level{level}']['coefficient_over_friedmann']:.12g} "
            f"error={row[f'level{level}']['relative_error']:.12g}"
            for level in range(3)
        )
    )
    factors = row["improvement_factors"]
    print(
        "    error-reduction factors: "
        f"{factors['level0_error_over_level1_error']:.9g}, "
        f"{factors['level1_error_over_level2_error']:.9g}"
    )
print(f"\nSummary: {passed}/{tests} checks passed")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")
raise SystemExit(0 if passed == tests else 1)

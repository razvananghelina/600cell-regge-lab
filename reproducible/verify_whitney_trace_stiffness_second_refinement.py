#!/usr/bin/env python3
"""Second barycentric refinement of exact Whitney trace stiffness.

Protocol commit 702fa5b froze the three complete control carriers, exact
iterated geometry, quotient solvers, dense calibration, residual gates, and
two-step comparison before the level-two spectra were evaluated.
"""

import json
from pathlib import Path

import sympy as sy

from whitney_trace_refinement_tools import (
    BLOCK_SIZE,
    LANCZOS_MAXITER,
    LANCZOS_TOLERANCE,
    LOBPCG_MAXITER,
    LOBPCG_TOLERANCE,
    audit_level,
    barycentric_refine,
    make_base_level,
)


OUTPUT = Path(__file__).with_name(
    "whitney_trace_stiffness_second_refinement.json"
)
CONTROL_CERTIFICATE = Path(__file__).with_name(
    "whitney_trace_stiffness.json"
)
PROTOCOL_COMMIT = "702fa5b"
EXPECTED_CONTROL_PROTOCOL = "b9a4104"
RANDOM_SEED = 60_020_260_811
RITZ_RESIDUAL_GATE = 1e-7
CALIBRATION_RELATIVE_GATE = 5e-7
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("SECOND REFINEMENT OF EXACT WHITNEY TRACE STIFFNESS")
print("=" * 78)

control = json.loads(CONTROL_CERTIFICATE.read_text())
check("the dense two-level control certificate has the frozen protocol",
      control["protocol_commit"] == EXPECTED_CONTROL_PROTOCOL)

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
level_zero = make_base_level(reference_vertices)
level_one = barycentric_refine(level_zero)
level_two = barycentric_refine(level_one)
levels = (level_zero, level_one, level_two)

expected_f_vectors = (
    [5, 10, 10, 5],
    [30, 150, 240, 120],
    [540, 3420, 5760, 2880],
)
check("all three exact f-vectors and duplicated dimensions are frozen values",
      all(list(map(len, level["cells"])) == expected
          for level, expected in zip(levels, expected_f_vectors))
      and [15 * len(level["top"]) for level in levels]
      == [75, 1800, 43200])

audits = []
for level in levels:
    print(f"\n-- level {level['level']} --")
    audits.append(audit_level(
        level, RANDOM_SEED + 100 * level["level"]
    ))

# Mandatory dense calibration of the identical quotient machinery.
calibration_errors = []
calibration_residuals = []
for level_index in range(2):
    targets = control["levels"][level_index]["trace_records"]
    for degree in range(3):
        observed = audits[level_index]["degree_records"][degree]
        target = targets[degree]
        calibration_errors.extend((
            abs(observed["positive_gap"] / target["positive_gap"] - 1.0),
            abs(
                observed["maximum_positive_eigenvalue"]
                / target["maximum_eigenvalue"] - 1.0
            ),
        ))
        calibration_residuals.append(
            observed["maximum_relative_ritz_residual"]
        )
maximum_calibration_error = max(calibration_errors)
maximum_calibration_residual = max(calibration_residuals)
check("levels zero and one reproduce every dense spectral edge",
      maximum_calibration_error < CALIBRATION_RELATIVE_GATE
      and maximum_calibration_residual < RITZ_RESIDUAL_GATE,
      f"max value error={maximum_calibration_error:.3e}, "
      f"max residual={maximum_calibration_residual:.3e}")

expected_rows = (
    [30, 30, 10],
    [720, 720, 240],
    [17280, 17280, 5760],
)
expected_ranks = (
    [15, 20, 10],
    [450, 570, 240],
    [10980, 13860, 5760],
)
check("all exact row and quotient-rank counts are reproduced",
      all(
          [record["constraint_rows"] for record in audit["degree_records"]]
          == rows
          and [record["constraint_rank"]
               for record in audit["degree_records"]] == ranks
          for audit, rows, ranks in zip(audits, expected_rows, expected_ranks)
      ))
check("every occurrence graph is connected and every shared-face metric agrees",
      all(
          record["all_occurrence_graphs_connected"]
          and record["face_metric_mismatches"] == 0
          for audit in audits for record in audit["degree_records"]
      ))
maximum_basis_residual = max(
    record["maximum_basis_orthonormality_residual"]
    for audit in audits for record in audit["degree_records"]
)
check("all exact row-image bases are orthonormal numerically",
      maximum_basis_residual < 1e-11,
      f"maximum residual={maximum_basis_residual:.3e}")
maximum_ritz_residual = max(
    record["maximum_relative_ritz_residual"]
    for audit in audits for record in audit["degree_records"]
)
check("all 18 extremal blocks meet the frozen Ritz residual gate",
      maximum_ritz_residual < RITZ_RESIDUAL_GATE,
      f"maximum residual={maximum_ritz_residual:.3e}")
check("all nine finite quotient gaps remain strictly positive",
      all(
          record["positive_gap"] > 1e-10
          for audit in audits for record in audit["degree_records"]
      ))

scales = [
    [record["a_over_gap"] for record in audit["degree_records"]]
    for audit in audits
]
ratios_01 = [scales[1][degree] / scales[0][degree]
             for degree in range(3)]
ratios_12 = [scales[2][degree] / scales[1][degree]
             for degree in range(3)]
flow_ratios = [ratios_12[degree] / ratios_01[degree]
               for degree in range(3)]
spread_01 = max(ratios_01) / min(ratios_01)
spread_12 = max(ratios_12) / min(ratios_12)
exact_repetition = all(
    abs(value - 1.0) < CALIBRATION_RELATIVE_GATE
    for value in flow_ratios
)
balance_improves = spread_12 < spread_01
check("the exact repetition label follows the frozen componentwise criterion",
      exact_repetition == all(
          abs(value - 1.0) < CALIBRATION_RELATIVE_GATE
          for value in flow_ratios
      ))
check("the balance-flow label follows the frozen spread comparison",
      balance_improves == (spread_12 < spread_01))

repetition_verdict = (
    "DERIVED NUMERICAL: the first-step factors repeat once on the control"
    if exact_repetition else
    "DERIVED NUMERICAL NEGATIVE: the first-step factors do not repeat exactly on the control"
)
balance_verdict = (
    "PATTERN: the second step flows toward a common degree scaling"
    if balance_improves else
    "PATTERN NEGATIVE: degree balance worsens at the second step"
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "control_protocol_commit": EXPECTED_CONTROL_PROTOCOL,
    "phenomenological_target_used": False,
    "candidate_count": 1,
    "solver": {
        "block_size": BLOCK_SIZE,
        "lobpcg_tolerance": LOBPCG_TOLERANCE,
        "lobpcg_maximum_iterations": LOBPCG_MAXITER,
        "degree_two_lanczos_tolerance": LANCZOS_TOLERANCE,
        "degree_two_lanczos_maximum_iterations": LANCZOS_MAXITER,
        "ritz_residual_gate": RITZ_RESIDUAL_GATE,
        "seed": RANDOM_SEED,
    },
    "calibration": {
        "maximum_value_relative_error": maximum_calibration_error,
        "maximum_ritz_relative_residual": maximum_calibration_residual,
    },
    "audits": audits,
    "two_step_comparison": {
        "scales_a_over_g": scales,
        "ratios_0_to_1": ratios_01,
        "ratios_1_to_2": ratios_12,
        "componentwise_flow_ratios": flow_ratios,
        "degree_spread_0_to_1": spread_01,
        "degree_spread_1_to_2": spread_12,
        "exact_factor_repetition": exact_repetition,
        "degree_balance_improves": balance_improves,
    },
    "verdicts": [
        repetition_verdict,
        balance_verdict,
        "OPEN: complete second-refined 600-cell and asymptotic law",
        "OPEN: absolute stiffness, chirality, and causal dynamics",
    ],
    "scope": (
        "Complete boundary-of-4-simplex control through two barycentric "
        "refinements. Not a complete second-refined 600-cell certificate."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured second-refinement certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("RATIOS_0_TO_1=" + str([float(value) for value in ratios_01]))
print("RATIOS_1_TO_2=" + str([float(value) for value in ratios_12]))
print("FLOW_RATIOS=" + str([float(value) for value in flow_ratios]))
print(f"SPREAD_0_TO_1={spread_01:.12g}")
print(f"SPREAD_1_TO_2={spread_12:.12g}")
print("REPETITION_VERDICT: " + repetition_verdict)
print("BALANCE_VERDICT: " + balance_verdict)
raise SystemExit(0 if passed == tests else 1)

#!/usr/bin/env python3
"""Whitney trace stiffness on the canonical rank-edgewise control tower.

Protocol commit c7e4335 froze the carrier, operator, solver gates, comparisons,
and labels before the new spectra were evaluated.
"""

import gc
import json
from pathlib import Path

import sympy as sy

from whitney_trace_refinement_tools import (
    audit_level,
    barycentric_refine,
    make_base_level,
    rank_edgewise_level,
)


OUTPUT = Path(__file__).with_name("whitney_rank_edgewise_stiffness.json")
CONTROL = Path(__file__).with_name(
    "whitney_trace_stiffness_second_refinement.json"
)
PROTOCOL_COMMIT = "c7e4335"
EXPECTED_CONTROL_PROTOCOL = "702fa5b"
RANDOM_SEED = 60_020_260_812
CALIBRATION_RELATIVE_GATE = 5e-7
RITZ_RESIDUAL_GATE = 1e-7
SHAPE_SCALING_GATE = 1e-10
REPETITION_GATE = 1e-6
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def checkpoint(audits, status):
    OUTPUT.write_text(json.dumps({
        "protocol_commit": PROTOCOL_COMMIT,
        "status": status,
        "completed_resolutions": [audit["edgewise_resolution"]
                                  for audit in audits],
        "audits": audits,
    }, indent=2, sort_keys=True) + "\n")


print("=" * 78)
print("WHITNEY TRACE STIFFNESS ON THE RANK-EDGEWISE TOWER")
print("=" * 78)

control = json.loads(CONTROL.read_text())
check("the frozen barycentric control certificate is available",
      control["protocol_commit"] == EXPECTED_CONTROL_PROTOCOL,
      f"control protocol={control['protocol_commit']}")

reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
base = make_base_level(reference_vertices)
ranked = barycentric_refine(base)

audits = []
expected_top = {1: 120, 2: 960, 4: 7680}
construction_ok = True
for resolution in (1, 2, 4):
    print(f"\n-- direct edgewise resolution k={resolution} --", flush=True)
    level = rank_edgewise_level(ranked, resolution)
    construction_ok &= len(level["top"]) == expected_top[resolution]
    print(
        f"f-vector={list(map(len, level['cells']))}, "
        f"top={len(level['top'])}",
        flush=True,
    )
    audit = audit_level(level, RANDOM_SEED + 100 * resolution)
    audit["edgewise_resolution"] = resolution
    audits.append(audit)
    checkpoint(audits, "IN_PROGRESS")
    del level
    gc.collect()

check("all preregistered top-cell counts are exact",
      construction_ok,
      f"observed={[audit['top_count'] for audit in audits]}")

calibration_errors = []
control_ranked = control["audits"][1]
for degree in range(3):
    observed = audits[0]["degree_records"][degree]
    target = control_ranked["degree_records"][degree]
    calibration_errors.extend((
        abs(observed["positive_gap"] / target["positive_gap"] - 1.0),
        abs(observed["maximum_positive_eigenvalue"]
            / target["maximum_positive_eigenvalue"] - 1.0),
    ))
maximum_calibration_error = max(calibration_errors)
check("k=1 reproduces all accepted barycentric spectral edges",
      maximum_calibration_error < CALIBRATION_RELATIVE_GATE,
      f"maximum relative error={maximum_calibration_error:.3e}")

geometry_ok = all(
    record["all_occurrence_graphs_connected"]
    and record["face_metric_mismatches"] == 0
    for audit in audits for record in audit["degree_records"]
)
check("every occurrence graph is connected and every face metric agrees",
      geometry_ok)

maximum_basis_residual = max(
    record["maximum_basis_orthonormality_residual"]
    for audit in audits for record in audit["degree_records"]
)
check("all row-image bases pass the orthonormality control",
      maximum_basis_residual < 1e-11,
      f"maximum residual={maximum_basis_residual:.3e}")

maximum_ritz_residual = max(
    record["maximum_relative_ritz_residual"]
    for audit in audits for record in audit["degree_records"]
)
check("all 18 extremal Ritz blocks pass the frozen residual gate",
      maximum_ritz_residual < RITZ_RESIDUAL_GATE,
      f"maximum residual={maximum_ritz_residual:.3e}")

shape_scaling_ratio = (
    audits[2]["maximum_local_dirac_norm"]
    / audits[1]["maximum_local_dirac_norm"]
)
check("the stable k=2 and k=4 shape set gives exact factor-two Dirac scaling",
      abs(shape_scaling_ratio - 2.0) < SHAPE_SCALING_GATE,
      f"a4/a2={shape_scaling_ratio:.12g}")

positive_gaps = all(
    record["positive_gap"] > 1e-10
    for audit in audits for record in audit["degree_records"]
)
check("all nine quotient gaps are strictly positive", positive_gaps)

scales = [[record["a_over_gap"] for record in audit["degree_records"]]
          for audit in audits]
ratios_12 = [scales[1][degree] / scales[0][degree]
             for degree in range(3)]
ratios_24 = [scales[2][degree] / scales[1][degree]
             for degree in range(3)]
flow_ratios = [ratios_24[degree] / ratios_12[degree]
               for degree in range(3)]
spread_12 = max(ratios_12) / min(ratios_12)
spread_24 = max(ratios_24) / min(ratios_24)
exact_repetition = all(abs(value - 1.0) < REPETITION_GATE
                       for value in flow_ratios)
balance_improves = spread_24 < spread_12

repetition_verdict = (
    "DERIVED NUMERICAL: exact refinement-factor repetition on the control"
    if exact_repetition else
    "DERIVED NUMERICAL NEGATIVE: refinement factors do not repeat on the control"
)
balance_verdict = (
    "PATTERN: degree balance improves on the shape-regular second step"
    if balance_improves else
    "PATTERN NEGATIVE: degree balance does not improve on the shape-regular second step"
)

barycentric_comparison = control["two_step_comparison"]
payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "status": "COMPLETE" if passed == tests else "FAILED_GATES",
    "phenomenological_target_used": False,
    "candidate_count": 1,
    "calibration_maximum_relative_error": maximum_calibration_error,
    "maximum_basis_orthonormality_residual": maximum_basis_residual,
    "maximum_ritz_relative_residual": maximum_ritz_residual,
    "shape_scaling_a4_over_a2": shape_scaling_ratio,
    "audits": audits,
    "two_step_comparison": {
        "scales_a_over_g": scales,
        "ratios_k1_to_k2": ratios_12,
        "ratios_k2_to_k4": ratios_24,
        "componentwise_flow_ratios": flow_ratios,
        "degree_spread_k1_to_k2": spread_12,
        "degree_spread_k2_to_k4": spread_24,
        "exact_factor_repetition": exact_repetition,
        "degree_balance_improves": balance_improves,
    },
    "barycentric_control_comparison": {
        "ratios_step_0_to_1": barycentric_comparison["ratios_0_to_1"],
        "ratios_step_1_to_2": barycentric_comparison["ratios_1_to_2"],
        "degree_spread_step_0_to_1": (
            barycentric_comparison["degree_spread_0_to_1"]
        ),
        "degree_spread_step_1_to_2": (
            barycentric_comparison["degree_spread_1_to_2"]
        ),
    },
    "verdicts": [
        repetition_verdict,
        balance_verdict,
        "OPEN: complete 600-cell rank-edgewise spectra",
        "OPEN: asymptotic continuum operator and absolute physical scales",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the complete structured certificate was written", OUTPUT.exists())

# Rewrite after the final certificate check so its count is included.
payload["status"] = "COMPLETE" if passed == tests else "FAILED_GATES"
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("RATIOS_K1_TO_K2=" + str([float(value) for value in ratios_12]))
print("RATIOS_K2_TO_K4=" + str([float(value) for value in ratios_24]))
print("FLOW_RATIOS=" + str([float(value) for value in flow_ratios]))
print(f"SPREAD_K1_TO_K2={spread_12:.12g}")
print(f"SPREAD_K2_TO_K4={spread_24:.12g}")
print("REPETITION_VERDICT: " + repetition_verdict)
print("BALANCE_VERDICT: " + balance_verdict)
raise SystemExit(0 if passed == tests else 1)


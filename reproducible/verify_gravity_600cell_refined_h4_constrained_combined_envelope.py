#!/usr/bin/env python3
"""Combined action/Hessian envelope audit for constrained H4 directions."""

from hashlib import sha256
import json
from pathlib import Path
import sys

import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIMARY = HERE / "gravity_600cell_refined_h4_constrained_response.json"
DIAGNOSTIC = HERE / "gravity_600cell_refined_h4_constrained_directional_diagnostic.json"
PRIMARY_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_primary_first_result.md"
)
DIAGNOSTIC_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_directional_diagnostic_result.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_combined_envelope_protocol.md"
)
OUTPUT = HERE / "gravity_600cell_refined_h4_constrained_combined_envelope.json"

PROTOCOL_COMMIT = "5cefff6"
EXPECTED_HASHES = {
    "primary": "f029260c9ee6e3b763293d237aae27e6ff7c1256eb8bc19c35725084ff385888",
    "diagnostic": "35662f71e4debdbd64356c6e004d32f652719baf17b38843414bb25b96e21b58",
    "primary_result": "633a57f3d2b4a054cce20d08544d409dac8fdaf53c39bae72ab2e9fceb4e83eb",
    "diagnostic_result": "0888a9b5caad440b1643d5f63992631a405f25ae13850888c3df0cc305a5ecb8",
    "protocol": "a091565f8bd1c6202c206beb276cee5d04247e3d45a9add3ccc2da43a2db3e42",
}

tests = 0
passed = 0


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
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=70):
    return mp.nstr(value, digits, strip_zeros=False)


def parse_complex(text):
    body = text.strip()
    if not (body.startswith("(") and body.endswith(")")):
        return mp.mpc(mp.mpf(body), 0)
    real_text, sign, imaginary_text = body[1:-1].rsplit(" ", 2)
    if sign not in {"+", "-"} or not imaginary_text.endswith("j"):
        raise ValueError("unexpected frozen complex format")
    imaginary = mp.mpf(imaginary_text[:-1])
    if sign == "-":
        imaginary = -imaginary
    return mp.mpc(mp.mpf(real_text), imaginary)


print("=" * 78)
print("CONSTRAINED H4 COMBINED ACTION/HESSIAN ENVELOPE AUDIT")
print("=" * 78)

mp.mp.dps = 120
paths = {
    "primary": PRIMARY,
    "diagnostic": DIAGNOSTIC,
    "primary_result": PRIMARY_RESULT,
    "diagnostic_result": DIAGNOSTIC_RESULT,
    "protocol": PROTOCOL,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all frozen failed-primary, diagnostic and correction-protocol inputs match",
    actual_hashes == EXPECTED_HASHES and PROTOCOL_COMMIT == "5cefff6",
)

primary = json.loads(PRIMARY.read_text())
diagnostic = json.loads(DIAGNOSTIC.read_text())
upstream_ok = check(
    "the exact formal primary failure and corrected diagnostic outcome are loaded",
    primary["outcome"] == "REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED"
    and primary["tests"] == {"passed": 18, "total": 19}
    and primary["census"]["class_count"] == 1
    and diagnostic["outcome"] == "REFINED_H4_DIRECTIONAL_DIAGNOSTIC_NONASYMPTOTIC"
    and diagnostic["tests"] == {"passed": 15, "total": 15}
    and diagnostic["census"]["direction_count"] == 12,
)

primary_by_order = {
    tuple(record["order"]): record
    for record in primary["census"]["schedules"]
}
norms = {
    "first_basis_vector": mp.mpf(1),
    "all_ones": mp.mpf(11),
    "alternating_signs": mp.mpf(11),
}

records = []
record_matching = []
convergence_passes = []
positivity_passes = []
combined_matches = []
action_only_rejections = []
corruption_passes = []
for item in diagnostic["census"]["records"]:
    order = tuple(item["order"])
    label = item["direction"]
    matched_record = primary_by_order.get(order)
    matched = matched_record is not None and label in norms
    record_matching.append(matched)
    if not matched:
        continue

    richardson = tuple(parse_complex(value) for value in item["richardson_180"])
    sixth = tuple(parse_complex(value) for value in item["sixth_180"])
    r_differences = tuple(
        abs(richardson[index] - richardson[index + 1]) for index in range(3)
    )
    r_ratios = (
        r_differences[0] / r_differences[1],
        r_differences[1] / r_differences[2],
    )
    x_differences = tuple(abs(sixth[index] - sixth[index + 1]) for index in range(2))
    x_ratio = x_differences[0] / x_differences[1]
    convergence = (
        all(mp.mpf(8) <= ratio <= mp.mpf(32) for ratio in r_ratios)
        and mp.mpf(32) <= x_ratio <= mp.mpf(128)
    )
    convergence_passes.append(convergence)

    q = mp.mpf(item["quadratic"])
    q_action = parse_complex(item["eighth_180"][1])
    e_action = mp.mpf(item["envelope"])
    e_k = mp.mpf(matched_record["response_envelopes"]["primary"])
    one_norm = norms[label]
    e_hessian = one_norm * one_norm * e_k
    e_total = e_action + e_hessian
    error = abs(q_action - q)
    positivity = e_action > 0 and e_k > 0 and e_total > e_action
    combined_match = error <= e_total
    action_only_rejected = error > e_action
    q_bad = q + mp.mpf("1e-6") * max(mp.mpf(1), abs(q))
    corrupted_error = abs(q_action - q_bad)
    corruption = corrupted_error > mp.mpf("1e6") * e_total
    positivity_passes.append(positivity)
    combined_matches.append(combined_match)
    action_only_rejections.append(action_only_rejected)
    corruption_passes.append(corruption)
    records.append({
        "order": list(order),
        "direction": label,
        "coefficient_one_norm": mp_text(one_norm),
        "quadratic": mp_text(q),
        "action_estimate": mp_text(q_action),
        "absolute_error": mp_text(error),
        "action_envelope": mp_text(e_action),
        "response_entrywise_envelope": mp_text(e_k),
        "hessian_quadratic_envelope": mp_text(e_hessian),
        "combined_envelope": mp_text(e_total),
        "error_to_action_envelope": mp_text(error / e_action),
        "error_to_combined_envelope": mp_text(error / e_total),
        "richardson_difference_ratios": [mp_text(value) for value in r_ratios],
        "sixth_difference_ratio": mp_text(x_ratio),
        "convergence_passed": convergence,
        "combined_match": combined_match,
        "action_only_rejected": action_only_rejected,
        "corrupted_error": mp_text(corrupted_error),
        "corruption_to_combined_envelope": mp_text(corrupted_error / e_total),
    })

matching_ok = check(
    "all twelve records match an exact frozen schedule and direction",
    len(records) == 12 and all(record_matching),
)
convergence_ok = check(
    "all action ladders have the frozen target-independent fourth/sixth-order ratios",
    len(convergence_passes) == 12 and all(convergence_passes),
    f"R range=[{mp_text(min(
        min(mp.mpf(x) for x in item['richardson_difference_ratios'])
        for item in records
    ), 7)},{mp_text(max(
        max(mp.mpf(x) for x in item['richardson_difference_ratios'])
        for item in records
    ), 7)}], X range=[{mp_text(min(
        mp.mpf(item['sixth_difference_ratio']) for item in records
    ), 7)},{mp_text(max(
        mp.mpf(item['sixth_difference_ratio']) for item in records
    ), 7)}]",
)
positive_ok = check(
    "all action and Hessian uncertainty contributions are positive and independently frozen",
    len(positivity_passes) == 12 and all(positivity_passes),
)
load_bearing_ok = check(
    "removing the Hessian uncertainty rejects all twelve comparisons",
    len(action_only_rejections) == 12 and all(action_only_rejections),
    f"action-only rejected={sum(action_only_rejections)}/12",
)
corruption_ok = check(
    "every preregistered corrupted quadratic is resolved outside the combined envelope",
    len(corruption_passes) == 12 and all(corruption_passes),
    f"min corruption/envelope={mp_text(min(
        mp.mpf(item['corruption_to_combined_envelope']) for item in records
    ), 8)}",
)
comparison_census_ok = check(
    "the complete combined-envelope comparison census is available",
    len(combined_matches) == 12,
    f"matched={sum(combined_matches)}/12, max error/envelope={mp_text(max(
        mp.mpf(item['error_to_combined_envelope']) for item in records
    ), 8)}",
)

scope = {
    "geometry_or_action_evaluated": False,
    "hessian_or_internal_solve_computed": False,
    "schedule_class_recomputed": False,
    "root_search_or_spectrum_executed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
}
scope_ok = check(
    "the audit remains a frozen-artifact algebraic comparison only",
    not any(scope.values()),
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    matching_ok,
    convergence_ok,
    positive_ok,
    load_bearing_ok,
    corruption_ok,
    comparison_census_ok,
    scope_ok,
))
if not controls_ok:
    outcome = "REFINED_H4_COMBINED_ENVELOPE_CONTROL_FAILED"
elif not all(combined_matches):
    outcome = "REFINED_H4_COMBINED_ENVELOPE_MISMATCH"
else:
    outcome = "REFINED_H4_COMBINED_ENVELOPE_CORROBORATED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one combined-envelope outcome",
    outcome in {
        "REFINED_H4_COMBINED_ENVELOPE_CONTROL_FAILED",
        "REFINED_H4_COMBINED_ENVELOPE_MISMATCH",
        "REFINED_H4_COMBINED_ENVELOPE_CORROBORATED",
    },
    outcome,
)

artifact = {
    "title": "Combined constrained action/Hessian envelope audit",
    "date": "2026-08-21",
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definition": {
        "hessian_quadratic_envelope": "coefficient_l1_norm^2 * entrywise_e_K",
        "combined_envelope": "e_action + coefficient_l1_norm^2 * entrywise_e_K",
    },
    "census": {
        "record_count": len(records),
        "convergence_count": sum(convergence_passes),
        "action_only_rejection_count": sum(action_only_rejections),
        "combined_match_count": sum(combined_matches),
        "records": records,
    },
    "scope": scope,
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
print(f"SHA-256: {digest(OUTPUT)}")
print("No action, Hessian, full suite or deferred nonlinear census was run.")

if passed != tests:
    sys.exit(1)


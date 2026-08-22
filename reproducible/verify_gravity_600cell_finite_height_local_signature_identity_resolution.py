#!/usr/bin/env python3
"""Exact resolution of the local-signature endpoint-identity width failure."""

import hashlib
import json
from pathlib import Path

import sympy as sp
from flint import arb, ctx


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = (
    ROOT
    / "docs"
    / "gravity"
    / "gravity_600cell_finite_height_local_signature_identity_resolution_protocol.md"
)
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_local_signature.json"
FIRST_OPEN_INPUT = (
    HERE / "gravity_600cell_finite_height_local_signature_adversarial.json"
)
MONOTONE_OPEN_INPUT = (
    HERE
    / "gravity_600cell_finite_height_local_signature_adversarial_resolution.json"
)
OUTPUT = (
    HERE / "gravity_600cell_finite_height_local_signature_identity_resolution.json"
)

PROTOCOL_COMMIT = "f2bf410"
PROTOCOL_SHA256 = (
    "f1b8f39a6f70ad1c803263b26f16159568dbebc098b26b7750976e2fbe56406a"
)
PRIMARY_SHA256 = (
    "9f524cc22df8cfb5083f372481b3efd19868252b85551d56378327eea7a6d613"
)
FIRST_OPEN_SHA256 = (
    "139dcee2e9ee021c131aae1090433fe16bd70c9f2b10ec52d32b0c5ebd7748a7"
)
MONOTONE_OPEN_SHA256 = (
    "70448c78be2156ef84fbaa986c543c6063bcca8ca4395ee77bdbf657ab2760d1"
)

ctx.dps = 240
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
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_sign(value):
    if value.lower() > 0:
        return 1
    if value.upper() < 0:
        return -1
    return 0


def arb_record(value):
    return {
        "pretty": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "contains_zero": bool(value.contains(0)),
        "sign": strict_sign(value),
    }


def series_asinh(value):
    return (value + (value * value + 1).sqrt()).log()


API = arb.pi()


def arb_epsilon(value):
    square = value * value
    return 2 * API - 5 * (
        (square + 2) / (2 * (square + 3))
    ).acos()


def arb_mu(value):
    return 180 * arb_epsilon(value) / (
        API * (value * value + 4).sqrt()
    )


primary = json.loads(PRIMARY_INPUT.read_text())
first_open = json.loads(FIRST_OPEN_INPUT.read_text())
monotone_open = json.loads(MONOTONE_OPEN_INPUT.read_text())
provenance_ok = bool(
    digest(PROTOCOL) == PROTOCOL_SHA256
    and digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and digest(FIRST_OPEN_INPUT) == FIRST_OPEN_SHA256
    and digest(MONOTONE_OPEN_INPUT) == MONOTONE_OPEN_SHA256
    and primary["outcome"] == "LOCAL_SIGNATURE_PRIMARY_CERTIFIED"
    and first_open["outcome"] == "LOCAL_SIGNATURE_ADVERSARIAL_OPEN"
    and monotone_open["outcome"]
    == "LOCAL_SIGNATURE_ADVERSARIAL_DISAGREEMENT_OPEN"
    and primary["passed"] == 10
    and primary["tests"] == 10
    and first_open["passed"] == 4
    and first_open["tests"] == 11
    and monotone_open["passed"] == 8
    and monotone_open["tests"] == 10
)
check(
    "the identity protocol and all three frozen result artifacts are exact",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}",
)


mu_symbol, p_symbol, m_symbol, pi_symbol, q_symbol = sp.symbols(
    "mu p m pi q", real=True
)
E_symbol = 4 * sp.pi * (mu_symbol - m_symbol) + q_symbol * (
    p_symbol - pi_symbol
)
h_symbol = (p_symbol - pi_symbol) / (2 * sp.pi * mu_symbol)
r_symbol = 2 * m_symbol / mu_symbol - 1
identity_remainder = sp.factor(
    r_symbol
    - (1 + h_symbol * q_symbol)
    + E_symbol / (2 * sp.pi * mu_symbol)
)
symbolic_ok = identity_remainder == 0
check(
    "the endpoint relation factors exactly through the root equation",
    symbolic_ok,
)


EXPECTED_COUNTS = {
    "root": (2, 1),
    "root/c0": (3, 2),
    "root/c0/c0": (2, 0),
    "root/c0/c1": (3, 1),
    "root/c0/c1/c0": (3, 1),
}
EXPECTED_TERMINALS = {
    "root/c0": None,
    "root/c0/c0": "DEAD",
    "root/c0/c1": None,
    "root/c0/c1/c0": None,
    "root/c0/c1/c0/c0": "ENTERED_D",
}

states = {}
physical_rows = []
tree_ok = True
certificates_ok = True
terminal_ok = True


def audit_state(state):
    global tree_ok, certificates_ok, terminal_ok
    path = state["path"]
    states[path] = state
    expected = EXPECTED_COUNTS.get(path)
    tree_ok &= bool(
        expected is not None
        and state["all_real_count"] == expected[0]
        and state["physical_count"] == expected[1]
    )
    for record in state["stationary_bisections"]:
        if record.get("exact_diagonal"):
            continue
        certificates_ok &= bool(
            record.get("certified")
            and record.get("monotonicity_certified")
        )
    for record in state["root_bisections"]:
        if record.get("exact_diagonal"):
            continue
        certificates_ok &= bool(
            record.get("certified")
            and record.get("monotonicity_certified")
        )
    for root in state["roots"]:
        if root["physical"]:
            physical_rows.append({"path": path, "root": root})
    for child in state["children"]:
        terminal_ok &= bool(
            child["path"] in EXPECTED_TERMINALS
            and child["terminal"] == EXPECTED_TERMINALS[child["path"]]
        )
        if child.get("next") is not None:
            audit_state(child["next"])


audit_state(monotone_open["tree"])
tree_audit_ok = bool(
    tree_ok
    and certificates_ok
    and terminal_ok
    and set(states) == set(EXPECTED_COUNTS)
    and len(physical_rows) == 5
)
check(
    "all five states retain certified roots, frozen counts and terminals",
    tree_audit_ok,
    f"states={len(states)}; physical_edges={len(physical_rows)}",
)


mu_rows = []
false_identity_rows = []
positivity_ok = True
negative_control_ok = True
for item in physical_rows:
    root = item["root"]
    q_value = arb(root["q"]["pretty"])
    h_value = arb(root["h"]["pretty"])
    r_value = arb(root["r"]["pretty"])
    mu_value = arb_mu(q_value)
    mu_positive = strict_sign(mu_value) > 0
    false_residual = r_value - (1 + 2 * h_value * q_value)
    false_rejected = strict_sign(false_residual) != 0
    positivity_ok &= mu_positive
    negative_control_ok &= false_rejected
    mu_rows.append(
        {
            "path": item["path"],
            "mu": arb_record(mu_value),
            "positive": mu_positive,
        }
    )
    false_identity_rows.append(
        {
            "path": item["path"],
            "residual": arb_record(false_residual),
            "rejected": false_rejected,
        }
    )
check(
    "mu is strictly positive on every certified physical root interval",
    positivity_ok and len(mu_rows) == 5,
)
check(
    "the false relation r=1+2*h*q is rejected on every physical edge",
    negative_control_ok and len(false_identity_rows) == 5,
)


stored_identity_ok = bool(
    len(monotone_open["transitions"]) == 5
    and all(
        row["identity"]["contains_zero"]
        for row in monotone_open["transitions"]
    )
)
comparison_ok = bool(
    len(monotone_open["post_construction_primary_comparison"]) == 5
    and all(
        row["passed"]
        for row in monotone_open["post_construction_primary_comparison"]
    )
    and monotone_open["preserved_open_endpoint_comparison"]
)
d_entry_ok = bool(
    len(monotone_open["D_entries"]) == 1
    and monotone_open["D_entries"][0]["passed"]
    and monotone_open["D_entries"][0]["hostile_126_rejected"]
)
artifact_audit_ok = stored_identity_ok and comparison_ok and d_entry_ok
check(
    "the stored zero-containing identities, comparisons and D entry are intact",
    artifact_audit_ok,
)


complete = bool(
    provenance_ok
    and symbolic_ok
    and tree_audit_ok
    and positivity_ok
    and negative_control_ok
    and artifact_audit_ok
)
outcome = (
    "LOCAL_SIGNATURE_ENDPOINT_IDENTITY_EXACTLY_RESOLVED"
    if complete
    else "LOCAL_SIGNATURE_ENDPOINT_IDENTITY_OPEN"
)
check(
    "the exact identity resolves the sole remaining adversarial failure",
    outcome == "LOCAL_SIGNATURE_ENDPOINT_IDENTITY_EXACTLY_RESOLVED",
    outcome,
)


artifact = {
    "provenance": {
        "protocol_commit": PROTOCOL_COMMIT,
        "protocol_sha256": PROTOCOL_SHA256,
        "primary_sha256": PRIMARY_SHA256,
        "first_open_sha256": FIRST_OPEN_SHA256,
        "monotone_open_sha256": MONOTONE_OPEN_SHA256,
    },
    "symbolic": {
        "identity": "r-(1+h*q)=-E/(2*pi*mu)",
        "remainder_zero": symbolic_ok,
    },
    "tree_audit": {
        "state_paths": sorted(states),
        "physical_edges": len(physical_rows),
        "certified": tree_audit_ok,
    },
    "mu_positive": mu_rows,
    "false_identity_controls": false_identity_rows,
    "stored_artifact_audit": {
        "all_five_numeric_identity_balls_contain_zero": stored_identity_ok,
        "primary_and_first_open_comparisons_pass": comparison_ok,
        "strict_D_entry_and_hostile_control_pass": d_entry_ok,
    },
    "claims": {
        "endpoint_identity": (
            "EXACTLY_RESOLVED" if complete else "OPEN"
        ),
        "local_signature": (
            "ADVERSARIALLY_CORROBORATED" if complete else "OPEN"
        ),
        "explicit_radius": "NOT_COMPUTED",
        "global_incoming_basin": "OPEN",
        "physical_selection_rule": "NOT_DERIVED",
        "nonhomogeneous_physics": "NOT_TESTED",
    },
    "passed": passed,
    "tests": tests,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print()
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
raise SystemExit(0 if passed == tests else 1)


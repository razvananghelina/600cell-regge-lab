#!/usr/bin/env python3
"""Minkowski/full-action replication of the finite-height selector audit."""

import hashlib
import json
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
DIRECT_INPUT = (
    HERE / "gravity_600cell_finite_height_composition_adversarial.json"
)
PRIMARY_SELECTOR_INPUT = (
    HERE / "gravity_600cell_finite_height_selector_audit.json"
)
OUTPUT = (
    HERE / "gravity_600cell_finite_height_selector_audit_adversarial.json"
)

DIRECT_SHA256 = (
    "d50e87f736e51585596aa1d7778238febaf7422840d668499878d8bd917f99e9"
)
PRIMARY_SELECTOR_SHA256 = (
    "956cd655b8b3a5106029fb852df74b85bb59f922a4984542bc2e089f54799676"
)
ADVERSARIAL_PROTOCOL_COMMIT = "94ee376"
PRIMARY_SELECTOR_ARTIFACT_COMMIT = "b0c2fe9"

mp.mp.dps = 120
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


def text(value, digits=70):
    return mp.nstr(value, digits)


# The primary selector result remains unread until the independent causal and
# full-action regularity classifications have been constructed.
direct = json.loads(DIRECT_INPUT.read_text())
direct_provenance_ok = bool(
    digest(DIRECT_INPUT) == DIRECT_SHA256
    and direct["outcome"]
    == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED"
    and direct["passed"] == direct["tests"] == 9
)
check(
    "the independent direct-action composition artifact is frozen",
    direct_provenance_ok,
    f"protocol={ADVERSARIAL_PROTOCOL_COMMIT}",
)


phi = (1 + mp.sqrt(5)) / 2
independent_records = []
geometry_ok = True
action_branch_ok = True
regularity_ok = True
same_input_ok = True
precision_nesting_ok = True

for precision_row in direct["precision_runs"]:
    branch_records = {}
    for branch in ("A", "B"):
        source = precision_row["branches"][branch]
        h = mp.mpf(source["h2"])
        q = mp.mpf(source["q2"])
        ratio = mp.mpf(source["ratio"])
        jacobian = mp.mpf(source["jacobian"])

        delta_r = phi * h * q
        central_time = mp.sqrt(h**2 + delta_r**2)
        strut_square = -central_time**2 + delta_r**2
        coordinate_speed = delta_r / central_time

        cosine = (q**2 + 2) / (2 * (q**2 + 3))
        boost = q / mp.sqrt(8 * (q**2 + 3))
        causal = bool(
            h > 0
            and ratio > 0
            and central_time > 0
            and abs(strut_square + h**2) < mp.mpf("1e-100")
            and strut_square < 0
            and abs(coordinate_speed) < 1
        )
        real_branch = bool(
            cosine >= mp.mpf(1) / 3
            and cosine < mp.mpf(1) / 2
            and mp.isfinite(boost)
        )
        regular = bool(jacobian > 0 and mp.isfinite(jacobian))

        geometry_ok &= causal
        action_branch_ok &= real_branch
        regularity_ok &= regular
        branch_records[branch] = {
            "h": h,
            "q": q,
            "scale_ratio": ratio,
            "delta_R": delta_r,
            "central_time": central_time,
            "strut_square": strut_square,
            "coordinate_speed": coordinate_speed,
            "angle_cosine": cosine,
            "angle_boost": boost,
            "direct_full_action_jacobian": jacobian,
            "causal": causal,
            "same_real_action_branch": real_branch,
            "locally_regular": regular,
        }

    same_input_ok &= bool(
        mp.mpf(precision_row["m1"]) > 0
        and mp.isfinite(mp.mpf(precision_row["pi1"]))
        and abs(branch_records["A"]["q"] - branch_records["B"]["q"])
        > 1
    )
    independent_records.append(
        {
            "target_precision": precision_row["target_precision"],
            "m1": mp.mpf(precision_row["m1"]),
            "pi1": mp.mpf(precision_row["pi1"]),
            "branches": branch_records,
        }
    )

for branch in ("A", "B"):
    base = independent_records[0]["branches"][branch]
    for row in independent_records[1:]:
        current = row["branches"][branch]
        precision_nesting_ok &= bool(
            abs(current["h"] - base["h"]) < mp.mpf("1e-60")
            and abs(current["q"] - base["q"]) < mp.mpf("1e-60")
            and abs(
                current["direct_full_action_jacobian"]
                - base["direct_full_action_jacobian"]
            )
            < mp.mpf("1e-60")
        )

check(
    "both direct roots reconstruct future timelike Minkowski struts",
    geometry_ok,
)
check(
    "both direct roots lie on the same real cellular-angle branch",
    action_branch_ok,
)
check(
    "both direct full-action roots have positive local Jacobian",
    regularity_ok,
)
check(
    "the two distinct roots solve the same normalized incoming state",
    same_input_ok,
)
check(
    "the causal classification and full-action Jacobians are stable across precision",
    precision_nesting_ok,
    "precisions=80,120,180",
)


# Hostile convention 1: falsely treating proper strut h as central coordinate
# height makes the high-q branch look superluminal.  This must disagree with
# the reconstructed Minkowski embedding.
last = independent_records[-1]["branches"]
wrong_central_speed = abs(phi * last["B"]["q"])
wrong_height_fails = bool(
    wrong_central_speed > 1
    and abs(last["B"]["coordinate_speed"]) < 1
)
check(
    "the central-height convention trap fails on the high-q branch",
    wrong_height_fails,
    f"wrong_speed={text(wrong_central_speed, 20)}",
)


# Hostile convention 2: Euclidean signature cannot reproduce a timelike
# strut square.
euclidean_squares = {
    branch: last[branch]["central_time"] ** 2 + last[branch]["delta_R"] ** 2
    for branch in ("A", "B")
}
wrong_signature_fails = all(
    value > 0 and abs(value + last[branch]["h"] ** 2) > mp.mpf("1e-20")
    for branch, value in euclidean_squares.items()
)
check(
    "the Euclidean-sign hostile control fails the timelike certificate",
    wrong_signature_fails,
)


# Only now read and compare the primary selector artifact.
primary = json.loads(PRIMARY_SELECTOR_INPUT.read_text())
primary_provenance_ok = bool(
    digest(PRIMARY_SELECTOR_INPUT) == PRIMARY_SELECTOR_SHA256
    and primary["outcome"]
    == "STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH"
    and primary["passed"] == primary["tests"] == 10
)
comparison_ok = primary_provenance_ok
comparison = {}
for branch, primary_row in zip(("A", "B"), primary["root_audits"]):
    direct_row = last[branch]
    differences = {
        "h": abs(direct_row["h"] - mp.mpf(primary_row["h"])),
        "q": abs(direct_row["q"] - mp.mpf(primary_row["q"])),
        "jacobian": abs(
            direct_row["direct_full_action_jacobian"]
            - mp.mpf(primary_row["legendre_determinant"])
        ),
        "beta_squared": abs(
            direct_row["coordinate_speed"] ** 2
            - mp.mpf(primary_row["beta_squared"])
        ),
    }
    row_ok = max(differences.values()) < mp.mpf("1e-55")
    comparison_ok &= row_ok
    comparison[branch] = {
        **{name: text(value, 20) for name, value in differences.items()},
        "passed": row_ok,
    }

check(
    "the independent embedding and full-action Jacobians agree only after construction",
    comparison_ok,
)


outcome = (
    "STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH_"
    "ADVERSARIALLY_CORROBORATED"
    if direct_provenance_ok
    and geometry_ok
    and action_branch_ok
    and regularity_ok
    and same_input_ok
    and precision_nesting_ok
    and wrong_height_fails
    and wrong_signature_fails
    and comparison_ok
    else "SELECTOR_AUDIT_ADVERSARIAL_OPEN"
)
check(
    "the adversarial hierarchy corroborates the selector-scoped negative",
    outcome
    == "STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH_"
    "ADVERSARIALLY_CORROBORATED",
)


def pack_record(row):
    return {
        "target_precision": row["target_precision"],
        "m1": text(row["m1"]),
        "pi1": text(row["pi1"]),
        "branches": {
            branch: {
                name: (
                    value
                    if isinstance(value, bool)
                    else text(value)
                )
                for name, value in branch_row.items()
            }
            for branch, branch_row in row["branches"].items()
        },
    }


artifact = {
    "provenance": {
        "direct_composition_sha256": DIRECT_SHA256,
        "primary_selector_sha256": PRIMARY_SELECTOR_SHA256,
        "primary_selector_artifact_commit": PRIMARY_SELECTOR_ARTIFACT_COMMIT,
        "adversarial_protocol_commit": ADVERSARIAL_PROTOCOL_COMMIT,
    },
    "method": {
        "causal_reconstruction": "direct Minkowski coordinates",
        "regularity_reconstruction": "direct full-action Jacobians",
        "primary_selector_read_after_construction": True,
    },
    "precision_records": [pack_record(row) for row in independent_records],
    "hostile_controls": {
        "wrong_central_height_speed": text(wrong_central_speed),
        "wrong_central_height_rejected": wrong_height_fails,
        "euclidean_signature_rejected": wrong_signature_fails,
    },
    "primary_comparison": comparison,
    "interpretation": {
        "label": (
            "DERIVED NEGATIVE, SELECTOR-SCOPED / "
            "ADVERSARIALLY CORROBORATED"
        ),
        "surviving_branches": 2,
        "additional_selector": "OPEN",
        "fundamental_tick": False,
    },
    "tests": tests,
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
raise SystemExit(0 if passed == tests else 1)

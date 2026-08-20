#!/usr/bin/env python3
"""Adversarial direct-minor replication of nonhomogeneous full rank."""

import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

from flint import acb, acb_mat, arb, ctx
import mpmath as mp
import numpy as np
from scipy import linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_canonical_precision_adversarial_protocol.md"
PRIMARY_RESULT = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_canonical_precision_primary_result.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_canonical_precision.py"
PRIMARY_INPUT = HERE / "gravity_600cell_full_scale_strut_canonical_precision.json"
OUTPUT = HERE / "gravity_600cell_full_scale_strut_canonical_precision_adversarial.json"

PROTOCOL_COMMIT = "ad03ede"
PRIMARY_COMMIT = "99efe9b"
EXPECTED_HASHES = {
    "protocol": "eab9fac2fd58a92fa1b3379849edde9b6983b77d5363ed9085bdcdb3d874fb62",
    "primary_result": "3a3d4e1b3d08a1707cd0873876576ec4a6307febb6cef3af74b8e3308d93db99",
    "primary_source": "b836cad394fe8a54644d514b6f31cb899ff5c3697b6c2a5f4edfc2b5f0ac5d62",
    "primary": "75351ae4dfde26dd75ed8faa927b0a49cd725d83c7629d4545268030b54e2706",
}
INPUTS = {
    "protocol": PROTOCOL,
    "primary_result": PRIMARY_RESULT,
    "primary_source": PRIMARY_SOURCE,
    "primary": PRIMARY_INPUT,
}

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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def midpoint_numpy(matrix):
    return np.array([
        [
            complex(float(mp.re(matrix[row, column])),
                    float(mp.im(matrix[row, column])))
            for column in range(matrix.cols)
        ]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


def selected_rows(matrix):
    midpoint = midpoint_numpy(matrix)
    columns = midpoint.shape[1]
    _, _, pivots = la.qr(midpoint.T, mode="economic", pivoting=True)
    return tuple(int(value) for value in pivots[:columns])


def component_ball(value, radius):
    radius_text = mp.nstr(radius, 80)
    real = arb(f"{mp.nstr(mp.re(value), 155)} +/- {radius_text}")
    imaginary = arb(f"{mp.nstr(mp.im(value), 155)} +/- {radius_text}")
    return acb(real, imaginary)


def direct_minor(matrix, rows, radius):
    columns = matrix.cols
    if len(rows) != columns:
        raise ValueError("a direct rank minor must be square")
    ball = acb_mat(
        columns,
        columns,
        [
            component_ball(matrix[rows[row], column], radius)
            for row in range(columns) for column in range(columns)
        ],
    )
    determinant = ball.det()
    return {
        "rows": list(rows),
        "contains_zero": bool(determinant.contains(0)),
        "determinant": str(determinant),
        "abs_lower": str(determinant.abs_lower()),
        "abs_upper": str(determinant.abs_upper()),
        "radius": str(determinant.rad()),
    }


def exact_direct_minor(matrix, rows):
    columns = matrix.cols
    ball = acb_mat(
        columns,
        columns,
        [
            acb(mp.nstr(mp.re(matrix[rows[row], column]), 80),
                mp.nstr(mp.im(matrix[rows[row], column]), 80))
            for row in range(columns) for column in range(columns)
        ],
    )
    determinant = ball.det()
    return bool(determinant.contains(0)), str(determinant)


print("=" * 78)
print("ADVERSARIAL DIRECT-MINOR NONHOMOGENEOUS RANK REPLICATION")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
primary_before = PRIMARY_INPUT.read_bytes()
primary = json.loads(primary_before)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and primary["outcome"] == "FULL_SCALE_STRUT_CANONICAL_HOMOGENEOUS_OPEN"
    and primary["passed"] == primary["tests"] == 17
    and primary["nonhomogeneous"]["all_full_rank"]
    and not primary["nonhomogeneous"]["binary_disagreement"]
    and set(primary["homogeneous_status"].values()) == {"OPEN"}
    and all(record["pass"] for record in primary["geometry_repair"].values())
)
check("all adversarial inputs retain frozen provenance", provenance_ok)

# Known pass/fail controls exercise the same QR selector and direct determinant.
identity = mp.matrix(7, 4)
for index in range(4):
    identity[index, index] = 1
identity_rows = selected_rows(identity)
identity_contains_zero, identity_determinant = exact_direct_minor(identity, identity_rows)
check(
    "the exact padded-identity direct minor excludes zero",
    not identity_contains_zero,
    f"rows={identity_rows}, determinant={identity_determinant}",
)

duplicate = mp.matrix(7, 5)
for row in range(7):
    for column in range(4):
        duplicate[row, column] = identity[row, column]
    duplicate[row, 4] = identity[row, 0]
duplicate_rows = selected_rows(duplicate)
duplicate_contains_zero, duplicate_determinant = exact_direct_minor(
    duplicate, duplicate_rows
)
check(
    "an exact duplicate column makes the direct minor contain zero",
    duplicate_contains_zero,
    f"rows={duplicate_rows}, determinant={duplicate_determinant}",
)

print("replaying the frozen primary construction under capture", flush=True)
captured = io.StringIO()
original_exit = sys.exit


def audited_exit(code=0):
    if code not in (None, 0):
        raise SystemExit(code)


try:
    sys.exit = audited_exit
    with contextlib.redirect_stdout(captured):
        frozen = runpy.run_path(str(PRIMARY_SOURCE))
finally:
    sys.exit = original_exit

primary_replay_ok = bool(
    frozen["tests"] == frozen["passed"] == 17
    and frozen["outcome"] == "FULL_SCALE_STRUT_CANONICAL_HOMOGENEOUS_OPEN"
    and PRIMARY_INPUT.read_bytes() == primary_before
    and digest(PRIMARY_INPUT) == EXPECTED_HASHES["primary"]
)
check("the frozen primary construction replays byte-identically", primary_replay_ok)

ctx.dps = 140
level_live = frozen["level_live"]
records = []
all_direct = True
nonhomogeneous_count = 0

for parity in ("even", "odd"):
    for sector_index, high in level_live["P160"][parity].items():
        if high["homogeneous"]:
            continue
        nonhomogeneous_count += 1
        low = level_live["P100"][parity][sector_index]
        for matrix_name, radius_name in (("D", "D_radius"), ("K", "K_radius")):
            low_rows = selected_rows(low[matrix_name])
            high_rows = selected_rows(high[matrix_name])
            low_to_high = direct_minor(
                high[matrix_name], low_rows, high[radius_name]
            )
            high_to_low = direct_minor(
                low[matrix_name], high_rows, low[radius_name]
            )
            local_ok = bool(
                not low_to_high["contains_zero"]
                and not high_to_low["contains_zero"]
            )
            all_direct &= local_ok
            records.append({
                "parity": parity,
                "sector_index": sector_index,
                "dimension": high["dimension"],
                "matrix": matrix_name,
                "P100_rows_P160_minor": low_to_high,
                "P160_rows_P100_minor": high_to_low,
                "pass": local_ok,
            })

expected_census = bool(
    nonhomogeneous_count == 12
    and len(records) == 24
    and sum(
        int(not record[direction]["contains_zero"])
        for record in records
        for direction in ("P100_rows_P160_minor", "P160_rows_P100_minor")
    ) == 48
)
check(
    "all 48 cross-precision nonhomogeneous direct minors exclude zero",
    all_direct and expected_census,
    f"sectors={nonhomogeneous_count}, matrices={len(records)}",
)

frozen_falsifiers_ok = bool(
    float(primary["nonhomogeneous"]["maximum_source_target_corruption_margin"]) > 1e100
    and all(
        record["D_P160_residual"] != "0.0"
        and record["K_P160_residual"] != "0.0"
        for record in primary["no_refit_validation"].values()
    )
    and all(
        item["passed"] and item["determinant_changed"]
        for item in primary["pole_deletion"].values()
    )
)
check("the frozen source/target and pole-deletion falsifiers remain active", frozen_falsifiers_ok)

outcome = (
    "NONHOMOGENEOUS_DIRECT_MINOR_REPLICATED"
    if passed == tests
    else "NONHOMOGENEOUS_DIRECT_MINOR_DISAGREEMENT"
)
check(
    "the preregistered adversarial hierarchy assigns one verdict",
    outcome in {
        "NONHOMOGENEOUS_DIRECT_MINOR_REPLICATED",
        "NONHOMOGENEOUS_DIRECT_MINOR_DISAGREEMENT",
    },
    outcome,
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "primary_commit": PRIMARY_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "controls": {
        "identity": {
            "rows": list(identity_rows),
            "contains_zero": identity_contains_zero,
            "determinant": identity_determinant,
        },
        "duplicate_column": {
            "rows": list(duplicate_rows),
            "contains_zero": duplicate_contains_zero,
            "determinant": duplicate_determinant,
        },
        "primary_replay_byte_identical": primary_replay_ok,
        "frozen_falsifiers_active": frozen_falsifiers_ok,
    },
    "nonhomogeneous_sector_count": nonhomogeneous_count,
    "direct_minor_certificate_count": sum(
        int(not record[direction]["contains_zero"])
        for record in records
        for direction in ("P100_rows_P160_minor", "P160_rows_P100_minor")
    ),
    "records": records,
    "classification": {
        "nonhomogeneous_canonical_intersection": (
            "ZERO; ADVERSARIALLY REPLICATED"
            if outcome == "NONHOMOGENEOUS_DIRECT_MINOR_REPLICATED"
            else "OPEN/DISAGREEMENT"
        ),
        "homogeneous_canonical_intersection": "OPEN; NOT RETESTED",
        "physical_interpretation": "NOT EVALUATED",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"TOTAL: {passed}/{tests} tests PASSED")
print(f"Artifact: {OUTPUT.name}")
if passed != tests:
    raise SystemExit(1)

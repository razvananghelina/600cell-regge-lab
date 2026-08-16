#!/usr/bin/env python3
"""Post-result audit of the committed complete-action continuation.

This verifier was written after the numerical result was known.  It checks
artifact integrity, recomputes the reported counts, and audits the frozen
35-dimensional carrier.  It does not upgrade the result's preregistered
provenance or turn a finite-grid pattern into a no-root theorem.
"""

from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESULT = HERE / "gravity_600cell_dust_action_solver_repair_continuation.json"
PRECISION = HERE / "gravity_600cell_dust_gauge_quotient_precision.json"
BOUNDARY = HERE / "gravity_global_boundary_legendre.json"

EXPECTED_SHA256 = (
    "14bac60c3c561db74290105d4049f15bd44ff056c597a72882aa16a4ee6d7719"
)
RESULT_COMMIT = "94675be"
PROTOCOL_COMMIT = "cf27934"

tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


payload = RESULT.read_bytes()
result = json.loads(payload)
precision = json.loads(PRECISION.read_text())
boundary = json.loads(BOUNDARY.read_text())

digest = hashlib.sha256(payload).hexdigest()
check(
    "the post-result auditor is bound to the committed machine artifact",
    digest == EXPECTED_SHA256,
    f"result commit={RESULT_COMMIT}, sha256={digest}",
)

check(
    "the continuation retains its preregistered provenance and precision",
    result["protocol_commit"] == PROTOCOL_COMMIT
    and result["prior_art_commit"] == "6b7f9e4"
    and result["frozen_result_commit"] == "64a13f6"
    and result["old_protocol_commit"] == "17f9560"
    and result["action_precision_digits"] == 100
    and result["tests"] == result["passed"] == 11,
)

check(
    "both target-blind derivative calibrations pass every frozen gate",
    all(
        result["calibration"][parity]["pass"]
        and all(result["calibration"][parity]["gates"].values())
        for parity in ("even", "odd")
    ),
)

cases = result["cases"]
states = [state for case in cases for state in case["grid_states"]]
state_keys = {
    (
        state["parity"], state["direction_index"], state["sign"], state["t"]
    )
    for state in states
}
expected_keys = {
    (parity, direction, sign, t)
    for parity in ("even", "odd")
    for direction in range(1, 5)
    for sign in (-1, 1)
    for t in (-0.10, -0.05, 0.0, 0.05, 0.10)
}
check(
    "the frozen look-elsewhere carrier contains exactly 16 cases and 80 states",
    len(cases) == 16 and len(states) == len(state_keys) == 80
    and state_keys == expected_keys,
)

outcomes = Counter(state["solver_outcome"] for state in states)
validated = [state for state in states if state["validation"] is not None]
unresolved = [state for state in states if state["validation"] is None]
check(
    "the three solver outcomes recompute to 63 validated plus 17 unresolved",
    outcomes == Counter({
        "SOLVER_TRANSVERSE_ZERO_CONSISTENT": 63,
        "NO_ROBUST_ACTION_DESCENT": 9,
        "SOLVER_ITERATION_LIMIT": 8,
    })
    and len(validated) == 63
    and len(unresolved) == 17,
    str(dict(outcomes)),
)

check(
    "validation occurs exactly for transverse-zero-consistent states",
    all(
        (state["validation"] is not None)
        == (state["solver_outcome"] == "SOLVER_TRANSVERSE_ZERO_CONSISTENT")
        for state in states
    ),
)

check(
    "all 63 independent validation rows pass every numerical and branch gate",
    all(
        all(state["validation"]["gates"].values())
        and state["validation"]["transverse_validated"]
        and state["validation"]["action_row"]["branch"]["pass"]
        and state["validation"]["action_row"]["imaginary_below_1e-70"]
        and state["validation"]["operational_validation_agreement"][
            "component_pass"
        ]
        and state["validation"]["operational_validation_agreement"][
            "transverse_pass"
        ]
        for state in validated
    ),
)

scalars = [
    float(state["validation"]["action_row"]["collective_scalar"])
    for state in validated
]
scalar_ratios = [state["validation"]["scalar_over_error"] for state in validated]
check(
    "every validated reduced scalar is resolved nonzero with positive sign",
    all(
        state["validation"]["scalar_label"]
        == "REDUCED_SCALAR_RESOLVED_NONZERO"
        and state["validation"]["scalar_sign"] == 1
        and state["scalar_label"] == "REDUCED_SCALAR_RESOLVED_NONZERO"
        for state in validated
    )
    and min(scalars) > 8.4e-17
    and max(scalars) < 1.04e-16
    and min(scalar_ratios) > 1.8e14,
    (
        f"g=[{min(scalars):.16e},{max(scalars):.16e}], "
        f"min |g|/error={min(scalar_ratios):.6e}"
    ),
)

check(
    "no unresolved state receives a scalar classification",
    all(state["scalar_label"] is None for state in unresolved),
)

case_outcomes = Counter(case["outcome"] for case in cases)
check(
    "six signed cases are grid-resolved and ten remain numerically open",
    case_outcomes == Counter({
        "NO_STATIONARY_POINT_IN_FROZEN_GRID_AND_BRACKETS": 6,
        "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED": 10,
    }),
    str(dict(case_outcomes)),
)

check(
    "no sign change forces a bisection and every preregistered hit count is zero",
    result["attempt_count"] == {
        "grid_states": 80,
        "mechanically_forced_bisection_midpoints": 0,
        "total": 80,
    }
    and result["hit_fractions"] == {
        "signed_cases": {"hits": 0, "total": 16},
        "direction_parity_pairs": {"hits": 0, "total": 8},
        "phase_contrasts": {"hits": 0, "total": 4},
    },
)

check(
    "the global label remains numerically unresolved rather than a no-root theorem",
    result["global_outcome"] == "ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED"
    and result["claim_boundary"] == {
        "finite_grid_result": "DERIVED COMPUTATIONAL IF VALIDATED",
        "positive_sign_provenance": "PATTERN-INFORMED",
        "continuous_interval_no_root": "NOT ESTABLISHED",
        "second_slab_and_full_carrier": "NOT TESTED",
        "physical_time_mass_speed_planck_units": "NOT DERIVED",
    },
)

check(
    "the upstream action distinguishes 35 internal constraints from 30 post-momenta",
    boundary["protocol_commit"] == "8c2482b"
    and boundary["tests"] == boundary["passed"] == 33
    and boundary["variables"]["internal"] == 35
    and boundary["variables"]["final_boundary"] == 30
    and boundary["labels"]["boundary_momenta"]
    == "DERIVED COMPUTATIONAL OFF SHELL",
)

quotient = np.asarray(precision["quotient_basis"], dtype=float)
zeta = (math.pi**2 * math.sqrt(2.0) / 50.0) ** (1.0 / 3.0)
l0_square = (zeta * 4.0 * 10.0 / (3.0 * math.pi)) ** 2
tau_square = 0.0102**2
minimum_singular = math.inf
full_ranks = []
for t in (-0.10, -0.05, 0.0, 0.05, 0.10):
    rho = tau_square * math.exp(t)
    raw = np.concatenate((
        np.full(30, -rho / (l0_square - rho)),
        np.ones(5),
    ))
    tangent = raw / np.linalg.norm(raw)
    full_basis = np.column_stack((quotient, tangent))
    singular_values = np.linalg.svd(full_basis, compute_uv=False)
    minimum_singular = min(minimum_singular, float(singular_values[-1]))
    full_ranks.append(int(np.linalg.matrix_rank(full_basis, tol=1e-12)))
check(
    "the 34 quotient directions plus the collective tangent span all 35 internals",
    quotient.shape == (35, 34)
    and full_ranks == [35] * 5
    and minimum_singular > 0.99999,
    f"minimum singular value={minimum_singular:.12f}",
)

print("-" * 78)
print(f"RESULT: {passed}/{tests} post-result audit checks passed")
print(
    "POST-RESULT AUDIT ONLY; 63 validated positive slopes, 17 unresolved "
    "states, and no stationary hit on the frozen grid."
)
raise SystemExit(0 if passed == tests else 1)

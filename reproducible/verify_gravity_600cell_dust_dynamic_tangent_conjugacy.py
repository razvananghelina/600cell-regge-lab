#!/usr/bin/env python3
"""Compare only the preregistered finite geometric conjugacy candidates.

Prior-art commit: 0129053.
Protocol commit: e006b8c.
Enumeration commit: 858512c.
"""

import hashlib
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
ENUMERATION = HERE / "gravity_600cell_dust_dynamic_tangent_conjugacy_enumeration.json"
TANGENT = HERE / "gravity_600cell_dust_dynamic_tangent.json"
OUTPUT = HERE / "gravity_600cell_dust_dynamic_tangent_conjugacy.json"

PRIOR_ART_COMMIT = "0129053"
PROTOCOL_COMMIT = "e006b8c"
ENUMERATION_COMMIT = "858512c"
ENUMERATION_SHA256 = "51b52457eba84ca1e41926b6e4fb1c51032f788b70bde916a3fb755d0323cb3e"
TANGENT_SHA256 = "1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5"
DPS = 100
STORE_ERROR = arb.mpf("3e-49")
CONSISTENT_FACTOR = arb.mpf(10)
FAIL_FACTOR = arb.mpf(100)

arb.mp.dps = DPS


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_from_rows(rows):
    return arb.matrix([[arb.mpf(value) for value in row] for row in rows])


def permutation_matrix(permutation):
    matrix = arb.matrix(30, 30)
    for source, target in enumerate(permutation):
        matrix[target, source] = 1
    return matrix


def phase_lift(permutation, momentum_sign):
    boundary = permutation_matrix(permutation)
    lift = arb.matrix(60, 60)
    for row in range(30):
        for column in range(30):
            lift[row, column] = boundary[row, column]
            lift[30+row, 30+column] = momentum_sign*boundary[row, column]
    return lift


def frobenius(matrix):
    return arb.sqrt(sum(abs(matrix[row, column])**2
                        for row in range(matrix.rows)
                        for column in range(matrix.cols)))


def classify(residual, uncertainty):
    if residual <= CONSISTENT_FACTOR*uncertainty:
        return "PASS"
    if residual > FAIL_FACTOR*uncertainty:
        return "FAIL"
    return "OPEN"


def number(value, digits=60):
    return arb.nstr(value, digits, strip_zeros=False)


enumeration = json.loads(ENUMERATION.read_text())
tangent = json.loads(TANGENT.read_text())

hashes = {
    "enumeration": digest(ENUMERATION),
    "tangent": digest(TANGENT),
}
provenance_ok = bool(
    hashes == {
        "enumeration": ENUMERATION_SHA256,
        "tangent": TANGENT_SHA256,
    }
    and enumeration.get("prior_art_commit") == PRIOR_ART_COMMIT
    and enumeration.get("protocol_commit") == PROTOCOL_COMMIT
    and enumeration.get("outcome") == "GEOMETRIC_CONJUGACY_CANDIDATES_ENUMERATED"
    and enumeration.get("passed") == enumeration.get("tests") == 7
    and enumeration.get("tangent_matrices_parsed") is False
    and enumeration.get("spectral_target_parsed") is False
)
tangent_ok = bool(
    tangent.get("prior_art_commit") == "25722d9"
    and tangent.get("protocol_commit") == "0bceb9b"
    and tangent.get("passed") == tangent.get("tests") == 12
    and tangent.get("number_of_maps") == 2
    and tangent.get("continuum_target_parsed") is False
    and tangent.get("speed_target_parsed") is False
    and tangent.get("full_720_edge_carrier") is False
)

matrices = {
    parity: matrix_from_rows(tangent["parities"][parity]["tangent_matrix"])
    for parity in ("even", "odd")
}
dimensions_ok = all(matrix.rows == matrix.cols == 60 for matrix in matrices.values())

deltas = {
    parity: arb.mpf(tangent["parities"][parity]["epsilon_t"])+STORE_ERROR
    for parity in ("even", "odd")
}
sigmas = {
    parity: (
        arb.mpf(tangent["parities"][parity]["full_spectrum"]["singular_values"][0])
        + arb.mpf(tangent["parities"][parity]["full_spectrum"]["epsilon_svd"])
    )
    for parity in ("even", "odd")
}
direct_uncertainty = arb.sqrt(60)*(deltas["even"]+deltas["odd"])
reversed_uncertainty = arb.sqrt(60)*(
    deltas["odd"]*sigmas["even"]
    + sigmas["odd"]*deltas["even"]
    + deltas["odd"]*deltas["even"]
)

omega = arb.matrix(60, 60)
for index in range(30):
    omega[index, 30+index] = 1
    omega[30+index, index] = -1

lift_controls = []
complete_attempts = []
boundary_attempts = []

for index, candidate in enumerate(enumeration["direct_slab_candidates"]):
    old_lift = phase_lift(candidate["old_to_old"], 1)
    final_lift = phase_lift(candidate["final_to_final"], 1)
    lift_controls.extend([
        frobenius(old_lift.T*omega*old_lift-omega),
        frobenius(final_lift.T*omega*final_lift-omega),
    ])
    residual = frobenius(
        matrices["odd"]*old_lift-final_lift*matrices["even"]
    )
    complete_attempts.append({
        "candidate_index": index,
        "kind": "DIRECT_COMPLETE_SLAB",
        "h4_action_count": candidate["h4_action_count"],
        "endpoint_permutations_equal": candidate["endpoint_permutations_equal"],
        "residual": residual,
        "uncertainty": direct_uncertainty,
        "ratio": residual/direct_uncertainty,
        "classification": classify(residual, direct_uncertainty),
    })

for index, candidate in enumerate(enumeration["reversed_slab_candidates"]):
    old_final = phase_lift(candidate["old_to_final"], -1)
    final_old = phase_lift(candidate["final_to_old"], -1)
    lift_controls.extend([
        frobenius(old_final.T*omega*old_final+omega),
        frobenius(final_old.T*omega*final_old+omega),
    ])
    residual = frobenius(
        matrices["odd"]*final_old*matrices["even"]-old_final
    )
    complete_attempts.append({
        "candidate_index": index,
        "kind": "REVERSED_COMPLETE_SLAB",
        "h4_action_count": candidate["h4_action_count"],
        "endpoint_permutations_equal": candidate["endpoint_permutations_equal"],
        "residual": residual,
        "uncertainty": reversed_uncertainty,
        "ratio": residual/reversed_uncertainty,
        "classification": classify(residual, reversed_uncertainty),
    })

for index, candidate in enumerate(enumeration["boundary_candidates"]):
    canonical = phase_lift(candidate["permutation"], 1)
    reversed_lift = phase_lift(candidate["permutation"], -1)
    lift_controls.extend([
        frobenius(canonical.T*omega*canonical-omega),
        frobenius(reversed_lift.T*omega*reversed_lift+omega),
    ])

    direct_residual = frobenius(
        matrices["odd"]*canonical-canonical*matrices["even"]
    )
    boundary_attempts.append({
        "candidate_index": index,
        "kind": "DIRECT_BOUNDARY",
        "sources": candidate["sources"],
        "h4_action_count": candidate["h4_action_count"],
        "residual": direct_residual,
        "uncertainty": direct_uncertainty,
        "ratio": direct_residual/direct_uncertainty,
        "classification": classify(direct_residual, direct_uncertainty),
    })

    reversed_residual = frobenius(
        matrices["odd"]*reversed_lift*matrices["even"]-reversed_lift
    )
    boundary_attempts.append({
        "candidate_index": index,
        "kind": "REVERSED_BOUNDARY",
        "sources": candidate["sources"],
        "h4_action_count": candidate["h4_action_count"],
        "residual": reversed_residual,
        "uncertainty": reversed_uncertainty,
        "ratio": reversed_residual/reversed_uncertainty,
        "classification": classify(reversed_residual, reversed_uncertainty),
    })

complete_counts = {
    label: sum(item["classification"] == label for item in complete_attempts)
    for label in ("PASS", "FAIL", "OPEN")
}
boundary_counts = {
    label: sum(item["classification"] == label for item in boundary_attempts)
    for label in ("PASS", "FAIL", "OPEN")
}

if not complete_attempts:
    complete_verdict = "NO_CROSS_PARITY_H4_SLAB_ISOMORPHISM"
elif complete_counts["PASS"]:
    complete_verdict = "COMPLETE_SLAB_COVARIANCE_DERIVED"
elif complete_counts["OPEN"]:
    complete_verdict = "COMPLETE_SLAB_COVARIANCE_OPEN"
else:
    complete_verdict = "COMPLETE_SLAB_COVARIANCE_REFUTED"

if boundary_counts["PASS"]:
    boundary_verdict = "BOUNDARY_INTERTWINER_STRUCTURAL"
elif boundary_counts["OPEN"]:
    boundary_verdict = "BOUNDARY_INTERTWINER_OPEN"
else:
    boundary_verdict = "BOUNDARY_INTERTWINER_REFUTED"

if complete_counts["PASS"]:
    outcome = "DYNAMIC_TANGENT_GEOMETRIC_COVARIANCE_DERIVED"
elif boundary_counts["PASS"]:
    outcome = "DYNAMIC_TANGENT_BOUNDARY_COVARIANCE_ONLY"
elif complete_counts["OPEN"] or boundary_counts["OPEN"]:
    outcome = "DYNAMIC_TANGENT_COVARIANCE_OPEN"
else:
    outcome = "DYNAMIC_TANGENT_ISOSPECTRALITY_UNEXPLAINED"

candidate_counts_ok = bool(
    enumeration["N_direct_slab"] == len(enumeration["direct_slab_candidates"])
    and enumeration["N_reversed_slab"] == len(enumeration["reversed_slab_candidates"])
    and enumeration["N_boundary"] == len(enumeration["boundary_candidates"]) == 60
)
lifts_ok = all(value < arb.mpf("1e-90") for value in lift_controls)
residuals_finite = all(
    arb.isfinite(item["residual"])
    and arb.isfinite(item["uncertainty"])
    and item["uncertainty"] > 0
    for item in complete_attempts+boundary_attempts
)
classifications_ok = all(
    item["classification"] in {"PASS", "FAIL", "OPEN"}
    for item in complete_attempts+boundary_attempts
)
denominator_ok = len(boundary_attempts) == 2*enumeration["N_boundary"] == 120

tests = [
    ("committed enumeration and tangent hashes", provenance_ok),
    ("blind tangent provenance and exclusions", tangent_ok),
    ("frozen candidate counts reproduced", candidate_counts_ok),
    ("both tangent matrices are 60 by 60", dimensions_ok),
    ("all permutation lifts are symplectic or anti-symplectic", lifts_ok),
    ("every frozen residual and uncertainty is finite", residuals_finite),
    ("every frozen candidate receives exactly one label", classifications_ok),
    ("boundary look-elsewhere denominator is 120", denominator_ok),
]
passed = sum(bool(ok) for _, ok in tests)
if passed != len(tests):
    outcome = "DYNAMIC_TANGENT_CONJUGACY_CONTROL_FAILED"


def serialize_attempt(item):
    return {
        **{key: value for key, value in item.items()
           if key not in {"residual", "uncertainty", "ratio"}},
        "residual": number(item["residual"]),
        "uncertainty": number(item["uncertainty"]),
        "ratio": number(item["ratio"]),
    }


payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "enumeration_commit": ENUMERATION_COMMIT,
    "input_sha256": hashes,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "full_720_edge_carrier": False,
    "N_complete": len(complete_attempts),
    "N_boundary_permutations": enumeration["N_boundary"],
    "N_boundary_attempts": len(boundary_attempts),
    "complete_counts": complete_counts,
    "boundary_counts": boundary_counts,
    "complete_verdict": complete_verdict,
    "boundary_verdict": boundary_verdict,
    "direct_uncertainty": number(direct_uncertainty),
    "reversed_uncertainty": number(reversed_uncertainty),
    "complete_attempts": [serialize_attempt(item) for item in complete_attempts],
    "boundary_attempts": [serialize_attempt(item) for item in boundary_attempts],
    "passed": passed,
    "tests": len(tests),
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"complete: {complete_verdict} {complete_counts}")
print(f"boundary: {boundary_verdict} {boundary_counts}")
if boundary_attempts:
    best = min(boundary_attempts, key=lambda item: item["ratio"])
    print(
        "best boundary diagnostic: index={} kind={} ratio={}".format(
            best["candidate_index"], best["kind"], number(best["ratio"], 10)
        )
    )
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) else 1)

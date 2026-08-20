#!/usr/bin/env python3
"""Classify the unique homogeneous weak line against the full pole equation."""

import hashlib
import json
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIOR = ROOT / "docs/gravity/gravity_600cell_homogeneous_pole_transversality_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_homogeneous_pole_transversality_protocol.md"
CONSOLIDATION = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial_result.md"
LAPSE_SOURCE = HERE / "verify_gravity_600cell_dust_homothetic_canonical_lapse.py"
LAPSE_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
PRIMARY_INPUT = HERE / "gravity_600cell_full_scale_strut_homogeneous_resolution.json"
ADVERSARIAL_INPUT = HERE / "gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial_p200g.json"
OUTPUT = HERE / "gravity_600cell_homogeneous_pole_transversality.json"

PROTOCOL_COMMIT = "72b0fad"
EXPECTED_HASHES = {
    "prior": "4b2edd2a520ab31dd581af5f7212967e0aac3dbb4a63c0beed118087174ce08b",
    "protocol": "5be9b1b3c940544f6165ffa4d2209ad1f742a804543901614abc2a2f9be55e06",
    "consolidation": "9a44770ee9abccdc47755e98c00bb355a833d4e0654fa304d9fe8b0cc203e895",
    "lapse_source": "8ae83004dcdeadfde27b91947a1c517915fa59af60807c9b128406a20c63508c",
    "lapse": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "primary": "70d7583756acdbee77893f98d57054ab074d9353a86247840cc1eb2c7b6be931",
    "adversarial": "fab74a26ae940cf0e65f26a4f6f167285cc269e282c40d7a630f37d65ba7ab07",
}
INPUTS = {
    "prior": PRIOR,
    "protocol": PROTOCOL,
    "consolidation": CONSOLIDATION,
    "lapse_source": LAPSE_SOURCE,
    "lapse": LAPSE_INPUT,
    "primary": PRIMARY_INPUT,
    "adversarial": ADVERSARIAL_INPUT,
}

mp.mp.dps = 120
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


def text(value, digits=70):
    return mp.nstr(value, digits)


def normalized(vector):
    norm = mp.sqrt(mp.fsum(abs(value) ** 2 for value in vector))
    return [value / norm for value in vector]


def projector_distance(left, right):
    left = normalized(left)
    right = normalized(right)
    return mp.sqrt(mp.fsum(
        abs(left[row] * mp.conj(left[column])
            - right[row] * mp.conj(right[column])) ** 2
        for row in range(2) for column in range(2)
    ))


print("=" * 78)
print("FULL POLE EQUATION ON THE UNIQUE HOMOGENEOUS WEAK LINE")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
lapse = json.loads(LAPSE_INPUT.read_text())
primary = json.loads(PRIMARY_INPUT.read_text())
adversarial = json.loads(ADVERSARIAL_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and lapse["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and lapse["passed"] == lapse["tests"] == 7
    and lapse["parity_gate"]["passed"]
    and primary["outcome"] == "HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE"
    and primary["passed"] == primary["tests"] == 10
    and adversarial["outcome"]
    == "HOMOGENEOUS_WEAK_POLE_LINE_REPLICATED_AFTER_CONTROL_REPAIR"
    and adversarial["passed"] == adversarial["tests"] == 7
)
check("all pole-transversality inputs retain frozen provenance", provenance_ok)

# A rank-one matrix has a pole row proportional to its momentum row, so the
# momentum-row kernel is also pole-null.  This is the negative control.
planted_a, planted_b, planted_c, planted_d = map(mp.mpf, (2, 4, 3, 6))
planted_v = (-planted_d, planted_c)
planted_momentum = planted_c * planted_v[0] + planted_d * planted_v[1]
planted_pole = planted_a * planted_v[0] + planted_b * planted_v[1]
check(
    "the planted rank-one matrix leaves its momentum kernel pole-null",
    planted_momentum == 0 and planted_pole == 0,
)

records = {}
projectors = {}
determinant_ok = True
convention_ok = True
for parity in ("even", "odd"):
    endpoint = lapse["solutions"][parity]["endpoint_jacobian"]
    matrix = endpoint["matrices"]["operational_primary"]
    a, b = (mp.mpf(value) for value in matrix[0])
    c, d = (mp.mpf(value) for value in matrix[1])
    epsilon = mp.mpf(endpoint["epsilon"])
    determinant = a * d - b * c
    stored_determinant = mp.mpf(endpoint["determinant"])
    determinant_bound = epsilon * (
        abs(a) + abs(b) + abs(c) + abs(d)
    ) + 2 * epsilon**2
    vector = (-d, c)
    momentum_derivative = c * vector[0] + d * vector[1]
    pole_derivative = a * vector[0] + b * vector[1]
    bridge = primary["bridges"][parity]
    primary_vector = (-mp.mpf(bridge["p_z"]), mp.mpf(bridge["p_s"]))
    line_distance = projector_distance(vector, primary_vector)
    projectors[parity] = vector

    parity_determinant_ok = bool(
        endpoint["resolved"] and endpoint["branch_pass"] and endpoint["entry_pass"]
        and abs(determinant) > 100 * determinant_bound
        and abs(determinant - stored_determinant) < mp.mpf("1e-45")
    )
    parity_convention_ok = bool(
        momentum_derivative == 0
        and pole_derivative == -determinant
        and line_distance < mp.mpf("1e-30")
    )
    determinant_ok &= parity_determinant_ok
    convention_ok &= parity_convention_ok
    records[parity] = {
        "jacobian": [[text(a), text(b)], [text(c), text(d)]],
        "epsilon": text(epsilon),
        "determinant": text(determinant),
        "determinant_error_bound": text(determinant_bound),
        "determinant_over_error_bound": text(abs(determinant) / determinant_bound),
        "weak_line_delta_s_delta_z": [text(vector[0]), text(vector[1])],
        "momentum_directional_derivative": text(momentum_derivative),
        "pole_directional_derivative": text(pole_derivative),
        "primary_line_projector_distance": text(line_distance),
        "determinant_certified": parity_determinant_ok,
        "convention_agrees": parity_convention_ok,
    }

parity_distance = projector_distance(projectors["even"], projectors["odd"])
parity_ok = bool(
    parity_distance < mp.mpf("1e-70")
    and mp.sign(mp.mpf(records["even"]["determinant"]))
    == mp.sign(mp.mpf(records["odd"]["determinant"]))
)
check("both endpoint determinants exclude zero inside their calibrated error", determinant_ok)
check("the exact weak-line convention is the momentum-row kernel", convention_ok)
check("both parities give the same transverse weak line", parity_ok,
      f"projector distance={text(parity_distance, 12)}")

if not provenance_ok or not parity_ok:
    outcome = "HOMOGENEOUS_POLE_TRANSVERSALITY_CONTROL_FAILED"
elif not convention_ok:
    outcome = "HOMOGENEOUS_POLE_TRANSVERSALITY_CONVENTION_DISAGREEMENT"
elif not determinant_ok:
    outcome = "HOMOGENEOUS_WEAK_LINE_FULLY_NULL_OPEN"
else:
    outcome = "HOMOGENEOUS_WEAK_LINE_TRANSVERSE_TO_POLE_EQUATION"

allowed = {
    "HOMOGENEOUS_POLE_TRANSVERSALITY_CONTROL_FAILED",
    "HOMOGENEOUS_POLE_TRANSVERSALITY_CONVENTION_DISAGREEMENT",
    "HOMOGENEOUS_WEAK_LINE_FULLY_NULL_OPEN",
    "HOMOGENEOUS_WEAK_LINE_TRANSVERSE_TO_POLE_EQUATION",
}
check("the preregistered pole-transversality hierarchy assigns one verdict",
      outcome in allowed, outcome)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "parities": records,
    "even_odd_projector_distance": text(parity_distance),
    "classification": {
        "weak_line_in_full_fixed_input_equations": (
            "TRANSVERSE/OFF-SHELL"
            if outcome == "HOMOGENEOUS_WEAK_LINE_TRANSVERSE_TO_POLE_EQUATION"
            else "OPEN"
        ),
        "accepted_nonstatic_endpoint": "PRESERVED",
        "local_fixed_input_solution_freedom": (
            "ZERO"
            if outcome == "HOMOGENEOUS_WEAK_LINE_TRANSVERSE_TO_POLE_EQUATION"
            else "OPEN"
        ),
        "discrete_evolution_map": "NOT REFUTED; NOT TESTED HERE",
        "free_tick_or_gauge": "NOT SUPPORTED BY THIS LINE",
        "tick_c_G_planck": "NOT DERIVED",
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


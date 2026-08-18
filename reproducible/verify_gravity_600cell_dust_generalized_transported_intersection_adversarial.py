#!/usr/bin/env python3
"""Independent float64 rank audit of the transported phase intersection.

Independence gate: 0d3a46c.
Preregistered protocol: 08ff16c.
"""

from collections import Counter
import contextlib
import hashlib
import io
import json
import math
from pathlib import Path
import runpy

import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ADVERSARIAL_SOURCE = (
    HERE / "verify_gravity_600cell_dust_generalized_phase_transport_adversarial.py"
)
ADVERSARIAL_ARTIFACT = (
    HERE / "gravity_600cell_dust_generalized_phase_transport_adversarial.json"
)
TANGENT_ARCHIVE = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"
EXACT_INTERSECTION = (
    HERE / "gravity_600cell_dust_generalized_transported_intersection.json"
)
OUTPUT = (
    HERE
    / "gravity_600cell_dust_generalized_transported_intersection_adversarial.json"
)

PRIOR_ART_COMMIT = "0d3a46c"
PROTOCOL_COMMIT = "08ff16c"
EXPECTED_HASHES = {
    "adversarial_source": (
        "f1cd1674af43573fd1c16b18bd37f7405093b580ce5cfa3ccad606ecb6a733cc"
    ),
    "adversarial_artifact": (
        "c33615ac6d0f3133e53077f46c5ee766b9c633d4d64c32124c24839c9c84c880"
    ),
    "tangent_archive": (
        "ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d"
    ),
    "exact_intersection": (
        "207cbe61bfaaf2b13d62cc3dbbb2ed5ea4931b7aab13cd47a8dd2802410c55c0"
    ),
}

PARITIES = ("even", "odd")
TARGET_SECTORS = (4, 5)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
MACHINE_EPSILON = np.finfo(float).eps
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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sf(value):
    return f"{float(value):.17e}"


def operator_norm(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    return float(singular[0]) if len(singular) else 0.0


def error_floor(matrix):
    return float(
        1000 * MACHINE_EPSILON * 60
        * max(1.0, operator_norm(matrix))
    )


def singular_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error):
        return "AUDIT_SINGULAR_OPEN"
    if value <= 10 * error:
        return "AUDIT_SINGULAR_ZERO_CONSISTENT"
    if value > 100 * error:
        return "AUDIT_SINGULAR_NONZERO_RESOLVED"
    return "AUDIT_SINGULAR_OPEN"


def leakage_spectrum(tangent, source, target):
    complement = la.null_space(target.conj().T)
    leakage = complement.conj().T @ tangent @ source
    singular = np.linalg.svd(leakage, compute_uv=False)
    error = error_floor(leakage)
    labels = [singular_label(value, error) for value in singular]
    return {
        "complement": complement,
        "leakage": leakage,
        "singular": singular,
        "error": error,
        "labels": labels,
        "counts": Counter(labels),
    }


hashes = {
    "adversarial_source": sha256(ADVERSARIAL_SOURCE),
    "adversarial_artifact": sha256(ADVERSARIAL_ARTIFACT),
    "tangent_archive": sha256(TANGENT_ARCHIVE),
    "exact_intersection": sha256(EXACT_INTERSECTION),
}
adversarial_artifact = json.loads(ADVERSARIAL_ARTIFACT.read_text())
exact_intersection = json.loads(EXACT_INTERSECTION.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and adversarial_artifact["passed"]
    == adversarial_artifact["tests"] == 8
    and adversarial_artifact["outcome"]
    == "ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED"
    and exact_intersection["passed"] == exact_intersection["tests"] == 7
    and exact_intersection["outcome"]
    == "TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL"
)
check("all frozen independent rank inputs retain exact provenance",
      provenance_ok, str(hashes))

print("[setup] replaying the independent Cholesky phase construction",
      flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    adversarial = runpy.run_path(str(ADVERSARIAL_SOURCE))
replay_ok = bool(
    adversarial["passed"] == adversarial["tests"] == 8
    and adversarial["outcome"]
    == "ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED"
    and sha256(ADVERSARIAL_ARTIFACT)
    == EXPECTED_HASHES["adversarial_artifact"]
    and len(adversarial["bases"]) == 32
)
check("the independent Cholesky phase construction replays byte-identically",
      replay_ok)

cell_records = []
singular_records = []
global_counts = Counter()
spectra_ok = True
gauge_ok = True
synthetic_ok = True
all_finite = True

with np.load(TANGENT_ARCHIVE) as tangent_archive:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for variant in VARIANTS:
                prefix = f"{parity}_sector{sector_index}_t2_{variant}"
                tangent = tangent_archive[f"{prefix}_midpoint"]
                old_basis = adversarial["bases"][
                    ("old", parity, sector_index, variant)
                ]
                shifted_basis = adversarial["bases"][
                    ("shifted", parity, sector_index, variant)
                ]
                source = adversarial["phase_basis"](
                    old_basis, "canonical_dual"
                )
                target = adversarial["phase_basis"](
                    shifted_basis, "canonical_dual"
                )
                actual = leakage_spectrum(tangent, source, target)
                counts = actual["counts"]
                global_counts.update(counts)

                positive_map = target @ source.conj().T
                positive = leakage_spectrum(positive_map, source, target)
                negative_map = actual["complement"] @ source.conj().T
                negative = leakage_spectrum(negative_map, source, target)

                phases0 = np.exp(1j * np.arange(30) / 7.0)
                phases1 = np.exp(-1j * np.arange(30) / 11.0)
                gauged_source = source[:, ::-1] * phases0
                gauged_target = target[:, ::-1] * phases1
                gauged = leakage_spectrum(
                    tangent, gauged_source, gauged_target
                )
                gauge_difference = float(np.max(np.abs(
                    actual["singular"] - gauged["singular"]
                )))

                cell_finite = bool(
                    np.all(np.isfinite(actual["singular"]))
                    and math.isfinite(actual["error"])
                    and math.isfinite(gauge_difference)
                )
                all_finite &= cell_finite
                spectra_ok &= bool(
                    actual["leakage"].shape == (30, 30)
                    and actual["complement"].shape == (60, 30)
                    and len(actual["singular"]) == 30
                    and all(
                        actual["singular"][index]
                        >= actual["singular"][index + 1]
                        for index in range(29)
                    )
                )
                gauge_ok &= gauge_difference <= 10 * actual["error"]
                synthetic_ok &= bool(
                    positive["counts"]["AUDIT_SINGULAR_ZERO_CONSISTENT"]
                    == 30
                    and negative["counts"][
                        "AUDIT_SINGULAR_NONZERO_RESOLVED"
                    ] == 30
                )
                resolved_rank = counts[
                    "AUDIT_SINGULAR_NONZERO_RESOLVED"
                ]
                cell_records.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "leakage_norm": sf(actual["singular"][0]),
                    "smallest_singular": sf(actual["singular"][-1]),
                    "error": sf(actual["error"]),
                    "smallest_error_units": sf(
                        actual["singular"][-1] / actual["error"]
                    ),
                    "numerical_rank": resolved_rank,
                    "zero_consistent": counts[
                        "AUDIT_SINGULAR_ZERO_CONSISTENT"
                    ],
                    "open": counts["AUDIT_SINGULAR_OPEN"],
                    "gauge_spectrum_difference": sf(gauge_difference),
                    "positive_control_rank": positive["counts"][
                        "AUDIT_SINGULAR_NONZERO_RESOLVED"
                    ],
                    "negative_control_rank": negative["counts"][
                        "AUDIT_SINGULAR_NONZERO_RESOLVED"
                    ],
                    "complete": bool(
                        cell_finite and resolved_rank == 30
                        and gauge_difference <= 10 * actual["error"]
                    ),
                })
                for index, (value, label) in enumerate(zip(
                    actual["singular"], actual["labels"]
                )):
                    singular_records.append({
                        "parity": parity,
                        "sector_index": sector_index,
                        "variant": variant,
                        "index": index,
                        "value": sf(value),
                        "error": sf(actual["error"]),
                        "error_units": sf(value / actual["error"]),
                        "label": label,
                    })

check("all 16 independent leakage spectra are finite and complete",
      spectra_ok and all_finite and len(cell_records) == 16)
check("all independent spectra are invariant under basis gauge",
      gauge_ok)
check("all positive-intersection and transverse controls discriminate",
      synthetic_ok)
census_ok = bool(
    len(singular_records) == 480
    and sum(global_counts.values()) == 480
)
check("all 480 independent singular values receive frozen labels",
      census_ok, str(dict(global_counts)))

controls_ok = bool(
    provenance_ok and replay_ok and spectra_ok and all_finite and gauge_ok
    and synthetic_ok and census_ok and len(cell_records) == 16
)
full_rank_cells = sum(item["numerical_rank"] == 30 for item in cell_records)
if not controls_ok:
    outcome = "ADVERSARIAL_INTERSECTION_CONTROL_FAILED"
elif full_rank_cells < 16:
    outcome = "ADVERSARIAL_INTERSECTION_DISAGREEMENT_OPEN"
else:
    outcome = "ADVERSARIAL_INTERSECTION_ZERO_CORROBORATED"

allowed = {
    "ADVERSARIAL_INTERSECTION_CONTROL_FAILED",
    "ADVERSARIAL_INTERSECTION_DISAGREEMENT_OPEN",
    "ADVERSARIAL_INTERSECTION_ZERO_CORROBORATED",
}
check("the preregistered adversarial rank hierarchy assigns one outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "exact_numeric_singular_values_read": False,
    "fitted_graph_used": False,
    "full_rank_cells": full_rank_cells,
    "global_label_counts": dict(global_counts),
    "cell_records": cell_records,
    "singular_records": singular_records,
    "classification": {
        "exact_zero_intersection": (
            "DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED"
            if outcome == "ADVERSARIAL_INTERSECTION_ZERO_CORROBORATED"
            else "OPEN PENDING DISAGREEMENT RESOLUTION"
        ),
        "this_replication": "STRUCTURAL INDEPENDENT CORROBORATION",
        "generalized_fiber_phase_route": (
            "CLOSED"
            if outcome == "ADVERSARIAL_INTERSECTION_ZERO_CORROBORATED"
            else "OPEN"
        ),
        "other_regge_perturbations": "NOT TESTED",
        "graph_propagator_dispersion_mass_inertia_or_speed": "NOT COMPUTED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"global independent labels: {dict(global_counts)}")
print(f"full-rank cells: {full_rank_cells}/16")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)


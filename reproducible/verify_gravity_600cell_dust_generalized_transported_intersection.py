#!/usr/bin/env python3
"""Certified rank census for the transported generalized phase intersection.

Prior-art/framing gate: 656b1d7.
Preregistered protocol: 1e6d8ce.
"""

from collections import Counter
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy

import mpmath as mp


HERE = Path(__file__).resolve().parent
PHASE_SOURCE = (
    HERE / "verify_gravity_600cell_dust_generalized_phase_transport.py"
)
PHASE_ARTIFACT = HERE / "gravity_600cell_dust_generalized_phase_transport.json"
ADVERSARIAL_ARTIFACT = (
    HERE / "gravity_600cell_dust_generalized_phase_transport_adversarial.json"
)
OUTPUT = (
    HERE / "gravity_600cell_dust_generalized_transported_intersection.json"
)

PRIOR_ART_COMMIT = "656b1d7"
PROTOCOL_COMMIT = "1e6d8ce"
EXPECTED_HASHES = {
    "phase_source": (
        "9c4c36b463a8faaa8d40b7db1b6b1852e3c04155c1b6ada4d02fbda747f6fcf3"
    ),
    "phase_artifact": (
        "45eb9a3e80ead758d9b3c2f8e1eccff44b06e2759251ab00c447aa53e6705743"
    ),
    "adversarial_artifact": (
        "c33615ac6d0f3133e53077f46c5ee766b9c633d4d64c32124c24839c9c84c880"
    ),
}

DPS = 100
PARITIES = ("even", "odd")
TARGET_SECTORS = (4, 5)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
mp.mp.dps = DPS
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


def smp(value, digits=30):
    value = mp.mpf(value)
    if mp.isinf(value):
        return "inf" if value > 0 else "-inf"
    if mp.isnan(value):
        return "nan"
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def singular_label(value, error):
    if not mp.isfinite(value) or not mp.isfinite(error):
        return "SINGULAR_OPEN"
    if value <= 10 * error:
        return "SINGULAR_ZERO_CONSISTENT"
    if value > 100 * error:
        return "SINGULAR_NONZERO_RESOLVED"
    return "SINGULAR_OPEN"


hashes = {
    "phase_source": sha256(PHASE_SOURCE),
    "phase_artifact": sha256(PHASE_ARTIFACT),
    "adversarial_artifact": sha256(ADVERSARIAL_ARTIFACT),
}
phase_artifact = json.loads(PHASE_ARTIFACT.read_text())
adversarial_artifact = json.loads(ADVERSARIAL_ARTIFACT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and phase_artifact["passed"] == phase_artifact["tests"] == 7
    and phase_artifact["outcome"] == "GENERALIZED_PHASE_TRANSPORT_REFUTED"
    and adversarial_artifact["passed"] == adversarial_artifact["tests"] == 8
    and adversarial_artifact["outcome"]
    == "ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED"
)
check("the accepted phase inputs retain exact provenance",
      provenance_ok, str(hashes))

print("[setup] replaying the accepted high-precision phase verifier", flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    phase = runpy.run_path(str(PHASE_SOURCE))
residual = phase["residual"]
phase_replay_ok = bool(
    phase["passed"] == phase["tests"] == 7
    and phase["outcome"] == "GENERALIZED_PHASE_TRANSPORT_REFUTED"
    and sha256(PHASE_ARTIFACT) == EXPECTED_HASHES["phase_artifact"]
    and len(residual["projectors"]) == 32
    and len(phase["tangent_cells"]) == 16
)
check("the accepted phase verifier replays byte-identically",
      phase_replay_ok)

committed_full = {
    (item["parity"], item["sector_index"], item["variant"]): item
    for item in phase_artifact["full_records"]
}
cell_records = []
singular_records = []
global_counts = Counter()
spectra_ok = True
norm_overlap_ok = True
structural_ok = True

for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for variant in VARIANTS:
            key = (parity, sector_index, variant)
            old = residual["projectors"][
                ("old", parity, sector_index, variant)
            ]
            shifted = residual["projectors"][
                ("shifted", parity, sector_index, variant)
            ]
            q0 = phase["block_diagonal_phase"](old["projector"])
            q1 = phase["block_diagonal_phase"](shifted["projector"])
            tangent = phase["acb_midpoint_to_mp"](
                phase["tangent_cells"][key]
            )
            residual_matrix = (mp.eye(60) - q1) * tangent * q0
            singular_matrix = mp.svd(residual_matrix, compute_uv=False)
            singulars = sorted(
                (mp.mpf(mp.re(singular_matrix[index])) for index in range(60)),
                reverse=True,
            )
            committed = committed_full[key]
            error = mp.mpf(committed["residual_error"])
            committed_norm = mp.mpf(committed["residual_norm"])
            norm_difference = abs(singulars[0] - committed_norm)
            cell_norm_overlap = norm_difference <= 10 * error
            norm_overlap_ok &= cell_norm_overlap
            finite = all(mp.isfinite(value) for value in singulars)
            ordered = all(
                singulars[index] >= singulars[index + 1]
                for index in range(59)
            )
            labels = [singular_label(value, error) for value in singulars]
            counts = Counter(labels)
            global_counts.update(counts)
            resolved_rank = counts["SINGULAR_NONZERO_RESOLVED"]
            lower_half_zero = all(
                label == "SINGULAR_ZERO_CONSISTENT"
                for label in labels[30:]
            )
            cell_structural = bool(
                resolved_rank <= 30 and lower_half_zero
            )
            structural_ok &= cell_structural
            spectra_ok &= bool(
                len(singulars) == 60 and finite and ordered
            )
            certified_dimension = 0 if resolved_rank == 30 else None
            cell_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "committed_residual_norm": smp(committed_norm),
                "recomputed_residual_norm": smp(singulars[0]),
                "residual_error": smp(error),
                "norm_difference": smp(norm_difference),
                "norm_overlap": bool(cell_norm_overlap),
                "resolved_rank_lower_bound": resolved_rank,
                "structural_rank_upper_bound": 30,
                "certified_intersection_dimension": certified_dimension,
                "zero_consistent": counts["SINGULAR_ZERO_CONSISTENT"],
                "open": counts["SINGULAR_OPEN"],
                "nonzero_resolved": resolved_rank,
                "lower_30_zero_consistent": bool(lower_half_zero),
                "complete": bool(
                    finite and ordered and cell_norm_overlap and cell_structural
                ),
            })
            for index, (value, label) in enumerate(zip(singulars, labels)):
                singular_records.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "index": index,
                    "value": smp(value),
                    "error": smp(error),
                    "error_units": smp(value / error),
                    "label": label,
                })

check("all 16 high-precision singular spectra are finite and ordered",
      spectra_ok and len(cell_records) == 16)
check("all recomputed spectral norms overlap their committed controls",
      norm_overlap_ok)
census_ok = bool(
    len(singular_records) == 960
    and sum(global_counts.values()) == 960
)
check("all 960 singular values receive frozen labels",
      census_ok, str(dict(global_counts)))
check("all spectra obey the structural rank-30 upper bound",
      structural_ok)

controls_ok = bool(
    provenance_ok and phase_replay_ok and spectra_ok and norm_overlap_ok
    and census_ok and structural_ok and len(cell_records) == 16
)
zero_cells = sum(
    item["certified_intersection_dimension"] == 0
    for item in cell_records
)
if not controls_ok:
    outcome = "TRANSPORTED_INTERSECTION_CONTROL_FAILED"
elif zero_cells == 16:
    outcome = "TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL"
else:
    outcome = "TRANSPORTED_INTERSECTION_DIMENSION_OPEN"

allowed = {
    "TRANSPORTED_INTERSECTION_CONTROL_FAILED",
    "TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL",
    "TRANSPORTED_INTERSECTION_DIMENSION_OPEN",
}
check("the preregistered intersection hierarchy assigns one outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "desired_rank_inspected": False,
    "fitted_graph_used": False,
    "cell_count": len(cell_records),
    "zero_intersection_cells": zero_cells,
    "global_label_counts": dict(global_counts),
    "cell_records": cell_records,
    "singular_records": singular_records,
    "classification": {
        "full_cotangent_phase_transport": "DERIVED COMPUTATIONAL REFUTATION",
        "transported_intersection": (
            "DERIVED COMPUTATIONAL ZERO"
            if outcome == "TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL"
            else "OPEN"
        ),
        "generalized_fiber_phase_route": (
            "CLOSED"
            if outcome == "TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL"
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
print(f"global singular labels: {dict(global_counts)}")
print(f"zero-intersection cells: {zero_cells}/16")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)


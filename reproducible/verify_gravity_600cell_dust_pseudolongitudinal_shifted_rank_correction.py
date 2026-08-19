#!/usr/bin/env python3
"""Target-disclosed logical correction of the shifted rank conjunction."""

import ast
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_pseudolongitudinal_shifted_rank_correction_protocol.md"
OPEN_RESULT = ROOT / "docs/gravity/gravity_600cell_dust_pseudolongitudinal_shifted_adversarial_result.md"
DIRECT_SOURCE = HERE / "verify_gravity_600cell_dust_pseudolongitudinal_shifted_adversarial.py"
DIRECT_JSON = HERE / "gravity_600cell_dust_pseudolongitudinal_shifted_adversarial.json"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_dust_pseudolongitudinal_shifted.py"
PRIMARY_JSON = HERE / "gravity_600cell_dust_pseudolongitudinal_shifted.json"
DIRECT_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_pseudolongitudinal_shifted_adversarial_protocol.md"
OUTPUT = HERE / "gravity_600cell_dust_pseudolongitudinal_shifted_rank_correction.json"

PROTOCOL_COMMIT = "5ba7830"
EXPECTED_HASHES = {
    "protocol": "c8cb3978eed5bff05b9ac215aee1080d278033d18d89a6696305f0cf3445f9ca",
    "open_result": "c97c72260df77bb0509a35cd815d6a68147afb02e01f9f711d35af125fc068b3",
    "direct_source": "72046f30f83e3af5192f60b108ea61e6b237dbecd11d615687b4c7c73417f521",
    "direct_json": "9e9f7253fd10422f3534914fae020857162862123fd4eae889e3570083552179",
    "primary_source": "e4c5bcc18007c1c0ba7fbd38e29dffcc33a526fd790dbfcba8defe2ae44b7ab2",
    "primary_json": "0480f5d49d24e0f5d8e4e95f0cf62b7d0d9242459ed2b8f6d8e835ecd6e103a7",
    "direct_protocol": "1dc9712a46b6ff6ac3c9b62e9d144f959f85622fe1ddbe2fc84de6ece3fa0982",
}
tests = passed = 0


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


paths = {
    "protocol": PROTOCOL,
    "open_result": OPEN_RESULT,
    "direct_source": DIRECT_SOURCE,
    "direct_json": DIRECT_JSON,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "direct_protocol": DIRECT_PROTOCOL,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all classifier-correction inputs have exact provenance",
      provenance_ok, str(hashes))

direct = json.loads(DIRECT_JSON.read_text())
primary = json.loads(PRIMARY_JSON.read_text())
direct_open_ok = bool(
    direct["outcome"] == "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_OPEN"
    and direct["passed"] == direct["tests"] == 18
    and direct["cells"] == direct["complete_cells"]
    == direct["target_rank_cells"] == 16
    and direct["changed_cells"] == 0
    and direct["centered_MV_archive_loaded"] is False
    and direct["primary_residual_loaded_only_post_census"] is True
)
check("the literal direct OPEN outcome and 16-cell census are preserved",
      direct_open_ok)

tree = ast.parse(DIRECT_SOURCE.read_text(), filename=str(DIRECT_SOURCE))
functions = {
    node.name: ast.get_source_segment(DIRECT_SOURCE.read_text(), node)
    for node in tree.body if isinstance(node, ast.FunctionDef)
}
relative_source = functions.get("direct_relative_cell", "")
source_shape_ok = bool(
    "left_bl, singular_bl, _ = la.svd(bl, full_matrices=False)"
    in relative_source
    and "(np.eye(25) - left_bl @ left_bl.conj().T) @ al"
    in relative_source
    and "rank_longitudinal" in relative_source
)
check("the frozen direct source constructs the orthogonal BL span residual",
      source_shape_ok)

records = direct["records"]
carrier_ok = bool(
    len(records) == 16
    and all(
        record["controls"]["exact_geometry_carrier"] is True
        and record["controls"]["kinetic_positive_definite_resolved"] is True
        and record["stiffness"]["sign_counts"]
        == {"NEGATIVE_RESOLVED": 15, "POSITIVE_RESOLVED": 10}
        and record["pseudolongitudinal"]["rank_longitudinal"] == 15
        and float(record["pseudolongitudinal"]["minimum_B_eigenvalue"]) > 0
        and min(
            float(record["pseudolongitudinal"]["AL_norm"]),
            float(record["pseudolongitudinal"]["GL_norm"]),
        ) > 1e-12
        and record["pseudolongitudinal"]["inequality_one"] is True
        and record["pseudolongitudinal"]["inequality_two"] is True
        for record in records
    )
)
check("all direct carriers, denominators, signs and norm inequalities pass",
      carrier_ok)

span_counts = Counter(
    record["pseudolongitudinal"]["span_label"] for record in records
)
comm_counts = Counter(
    record["pseudolongitudinal"]["commutator_label"] for record in records
)
direct_nonzero = bool(
    span_counts == {"NONZERO_RESOLVED": 16}
    and comm_counts == {"NONZERO_RESOLVED": 16}
)
direct_zero = bool(
    span_counts == {"ZERO_CONSISTENT": 16}
    and comm_counts == {"ZERO_CONSISTENT": 16}
)
check("the direct residual labels receive a complete non-majority census",
      direct_nonzero or direct_zero,
      f"span={dict(span_counts)}, commutator={dict(comm_counts)}")

auxiliary_open_preserved = all(
    record["pseudolongitudinal"]["augmented_rank"] == 15
    and 0 < float(record["pseudolongitudinal"]["sixteenth_augmented_singular"])
    < float(record["pseudolongitudinal"]["augmented_threshold"])
    for record in records
)
primary_agrees = bool(
    primary["outcome"] == "SHIFTED_PSEUDOLONGITUDINAL_DEFECT_PERSISTS"
    and primary["passed"] == primary["tests"] == 10
    and primary["label_counts"]["relative_span"]
    == {"NONZERO_RESOLVED": 16}
    and primary["label_counts"]["relative_commutator"]
    == {"NONZERO_RESOLVED": 16}
)
check("the auxiliary rank failure is preserved and the primary route agrees",
      auxiliary_open_preserved and primary_agrees)

controls_ok = bool(
    provenance_ok and direct_open_ok and source_shape_ok and carrier_ok
    and auxiliary_open_preserved and primary_agrees
)
if not controls_ok:
    outcome = "SHIFTED_PSEUDOLONGITUDINAL_RANK_CORRECTION_CONTROL_FAILED"
elif direct_nonzero:
    outcome = "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_RESIDUAL_CONFIRMATION"
elif direct_zero:
    outcome = "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_RESIDUAL_REFUTATION"
else:
    outcome = "SHIFTED_PSEUDOLONGITUDINAL_RANK_CORRECTION_OPEN"
allowed = {
    "SHIFTED_PSEUDOLONGITUDINAL_RANK_CORRECTION_CONTROL_FAILED",
    "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_RESIDUAL_CONFIRMATION",
    "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_RESIDUAL_REFUTATION",
    "SHIFTED_PSEUDOLONGITUDINAL_RANK_CORRECTION_OPEN",
}
check("the frozen correction hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "preserved_direct_outcome": direct["outcome"],
    "preserved_auxiliary_augmented_rank_counts":
        direct["pseudolongitudinal_label_counts"]["augmented_rank"],
    "direct_residual_label_counts": {
        "relative_span": dict(span_counts),
        "relative_commutator": dict(comm_counts),
    },
    "exact_rank_lemma":
        "rank([X,Y]) = rank(X) + rank((I-P_X)Y)",
    "classification": {
        "two_tick_temporal_noninvariance": (
            "CONFIRMED DERIVED COMPUTATIONAL / STRUCTURAL"
            if outcome.endswith("CONFIRMATION") else
            "REFUTED" if outcome.endswith("REFUTATION") else "OPEN"
        ),
        "original_direct_rank_conjunction": "PRESERVED OPEN",
        "curvature_refinement_continuum_or_physics": "NOT TESTED / OPEN",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
print("preserved direct outcome:", direct["outcome"])
print("direct span labels:", dict(span_counts))
print("direct commutator labels:", dict(comm_counts))
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)


#!/usr/bin/env python3
"""Compare committed third-tick stationary roots with the canonical target.

The scalar result direction was observed in a disclosed shell diagnostic after
root commit 3401137 and before specification commit 3f665fc.  This verifier is
therefore a reproducible exhaustive comparison, not a clean preregistration.
"""

import hashlib
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
ROOT_ARTIFACT = HERE / "gravity_600cell_dust_third_tick_stationary_roots.json"
SECOND_TICK_ARTIFACT = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_third_tick_stationary_target.json"
PRIOR_ART_COMMIT = "7b9a676"
ROOT_RESULT_COMMIT = "3401137"
SPECIFICATION_COMMIT = "3f665fc"
ROOT_SHA256 = "02d4589a7df0851c67a31fc0a41c5ef8851a82c758214c1c5e8729afddfe479f"
SECOND_TICK_SHA256 = "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70"
GLUING_SHA256 = "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77"
arb.mp.dps = 100
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
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(value):
    return arb.mpf(str(value))


def text(value, digits=55):
    return arb.nstr(value, digits)


def norm(values):
    return arb.sqrt(sum(value*value for value in values))


def spread(values):
    average = sum(values, arb.mpf(0))/len(values)
    return max(abs(value-average) for value in values)


roots = json.loads(ROOT_ARTIFACT.read_text())
second_tick = json.loads(SECOND_TICK_ARTIFACT.read_text())
gluing = json.loads(GLUING_ARTIFACT.read_text())

hashes_ok = bool(
    digest(ROOT_ARTIFACT) == ROOT_SHA256
    and digest(SECOND_TICK_ARTIFACT) == SECOND_TICK_SHA256
    and digest(GLUING_ARTIFACT) == GLUING_SHA256
)
provenance_ok = bool(
    PRIOR_ART_COMMIT == "7b9a676"
    and ROOT_RESULT_COMMIT == "3401137"
    and SPECIFICATION_COMMIT == "3f665fc"
    and roots["outcome"] == "THIRD_TICK_STATIONARY_ROOTS_ENUMERATED"
    and roots["passed"] == roots["tests"] == 5
    and roots["target_parsed"] is False
    and roots["candidate_count"] == len(roots["roots"]) == 2
    and second_tick["outcome"] == "SECOND_HOMOTHETIC_TICK_ACCEPTED"
    and second_tick["passed"] == second_tick["tests"] == 6
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing["passed"] == gluing["tests"] == 25
)
check(
    "the committed root-first provenance and exact input hashes pass",
    hashes_ok and provenance_ok,
    "root={} tick={} glue={}".format(
        digest(ROOT_ARTIFACT), digest(SECOND_TICK_ARTIFACT),
        digest(GLUING_ARTIFACT),
    ),
)

parities = ("even", "odd")
targets = {}
bounds = {}
maps = {}
map_ok = True
for parity in parities:
    mapping = tuple(gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"])
    post = tuple(
        number(value) for value in second_tick["solutions"][parity]["post_momentum"]
    )
    maps[parity] = mapping
    targets[parity] = tuple(post[index] for index in mapping)
    bounds[parity] = number(second_tick["solutions"][parity]["junction_bound"])
    map_ok &= bool(
        len(mapping) == 30 and sorted(mapping) == list(range(30))
        and len(post) == 30 and bounds[parity] > 0
    )
check(
    "both complete canonical targets use the independently stored maps",
    map_ok,
    f"bounds={text(bounds['even'])}/{text(bounds['odd'])}",
)

A1 = number(roots["geometry"]["A1"])
B2 = number(roots["geometry"]["B2"])
labels = []
label_ok = True
for root in roots["roots"]:
    c = number(root["parities"]["even"]["upper_log"])
    if c < B2:
        label = "CONTRACTING"
    elif (
        root["kind"] == "node_cluster"
        and abs(c-A1) < number("1e-25")
    ):
        label = "TIME_REVERSAL"
    else:
        label = "UNCLASSIFIED"
    labels.append(label)
    label_ok &= all(
        len(root["parities"][parity]["evaluation"]["pre_momentum"]) == 30
        for parity in parities
    )
label_ok &= labels == ["CONTRACTING", "TIME_REVERSAL"]
check(
    "structural labels are assigned without target proximity",
    label_ok,
    f"labels={labels}",
)

comparisons = []
hit_indices = []
comparison_ok = True
for root_index, root in enumerate(roots["roots"]):
    parity_records = {}
    root_hit = True
    for parity in parities:
        pre = tuple(
            number(value)
            for value in root["parities"][parity]["evaluation"]["pre_momentum"]
        )
        residual = tuple(left-right for left, right in zip(pre, targets[parity]))
        residual_norm = norm(residual)
        maximum = max(abs(value) for value in residual)
        residual_spread = spread(residual)
        hit = residual_norm <= bounds[parity]
        root_hit &= hit
        comparison_ok &= bool(
            len(residual) == 30 and arb.isfinite(residual_norm)
            and arb.isfinite(maximum) and arb.isfinite(residual_spread)
        )
        parity_records[parity] = {
            "target": [text(value) for value in targets[parity]],
            "residual": [text(value) for value in residual],
            "norm": text(residual_norm),
            "maximum_absolute_component": text(maximum),
            "component_spread": text(residual_spread),
            "junction_bound": text(bounds[parity]),
            "hit": bool(hit),
        }
    if root_hit:
        hit_indices.append(root_index)
    comparisons.append({
        "root_index": root_index,
        "structural_label": labels[root_index],
        "parities": parity_records,
        "hit_both_parities": bool(root_hit),
    })
check(
    "both roots are compared on all components and both parities",
    comparison_ok and len(comparisons) == 2,
)

controls_ok = bool(hashes_ok and provenance_ok and map_ok and label_ok)
if not controls_ok:
    outcome = "THIRD_TICK_TARGET_CONTROL_FAILED"
elif hit_indices:
    outcome = "STATIONARY_THIRD_TICK_HIT"
else:
    outcome = "STATIONARY_THIRD_TICK_NO_HIT"
classification_ok = bool(
    (not controls_ok and outcome == "THIRD_TICK_TARGET_CONTROL_FAILED")
    or (controls_ok and hit_indices and outcome == "STATIONARY_THIRD_TICK_HIT")
    or (controls_ok and not hit_indices and outcome == "STATIONARY_THIRD_TICK_NO_HIT")
)
check(
    "the disclosed deterministic hierarchy assigns the outcome",
    classification_ok,
    f"outcome={outcome}, hits={len(hit_indices)}/2",
)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "root_result_commit": ROOT_RESULT_COMMIT,
    "specification_commit": SPECIFICATION_COMMIT,
    "clean_comparison_preregistration": False,
    "process_disclosure": (
        "scalar residual direction observed after root commit and before "
        "comparison specification commit"
    ),
    "root_enumeration_target_firewall_intact": True,
    "input_sha256": {
        "stationary_roots": digest(ROOT_ARTIFACT),
        "second_tick": digest(SECOND_TICK_ARTIFACT),
        "gluing": digest(GLUING_ARTIFACT),
    },
    "candidate_count": 2,
    "hit_count": len(hit_indices),
    "hit_fraction": f"{len(hit_indices)}/2",
    "hit_indices": hit_indices,
    "comparisons": comparisons,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

for record in comparisons:
    even = record["parities"]["even"]
    print(
        "root {} {}: norm={} max={} bound={} hit={}".format(
            record["root_index"], record["structural_label"], even["norm"],
            even["maximum_absolute_component"], even["junction_bound"],
            record["hit_both_parities"],
        )
    )
print(f"OUTCOME: {outcome}")
print(f"Tests passed: {passed}/{tests}")
raise SystemExit(0 if passed == tests else 1)

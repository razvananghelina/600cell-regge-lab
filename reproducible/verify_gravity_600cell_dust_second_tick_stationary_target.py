#!/usr/bin/env python3
"""Compare the preregistered stationary roots with the second-tick target.

Prior-art commit: fcc4d7c.
Protocol commits: 9127b04, 28c4dd1.
"""

import hashlib
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
ROOT_ARTIFACT = HERE / "gravity_600cell_dust_stationary_root_enumeration.json"
TICK_ARTIFACT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_second_tick_stationary_target.json"
PRIOR_ART_COMMIT = "fcc4d7c"
PROTOCOL_COMMIT = "9127b04"
CONTROL_CORRECTION_COMMIT = "28c4dd1"
ROOT_SHA256 = "0ec5ba520ea25b39dd6cfd3c349d49fe480df2abee359854e1316b5af4d9fa2f"
TICK_SHA256 = "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9"
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
tick = json.loads(TICK_ARTIFACT.read_text())
gluing = json.loads(GLUING_ARTIFACT.read_text())

hashes_ok = bool(
    digest(ROOT_ARTIFACT) == ROOT_SHA256
    and digest(TICK_ARTIFACT) == TICK_SHA256
    and digest(GLUING_ARTIFACT) == GLUING_SHA256
)
provenance_ok = bool(
    PRIOR_ART_COMMIT == "fcc4d7c"
    and PROTOCOL_COMMIT == "9127b04"
    and CONTROL_CORRECTION_COMMIT == "28c4dd1"
    and roots.get("outcome") == "STATIONARY_ROOTS_ENUMERATED"
    and roots.get("passed") == roots.get("tests") == 5
    and roots.get("protocol_commit") == "07083cc"
    and roots.get("target_parsed") is False
    and len(roots.get("roots", [])) == 2
    and tick.get("outcome") == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and tick.get("passed") == tick.get("tests") == 7
    and tick.get("root_result_commit") == "b788258"
    and gluing.get("outcome") == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing.get("passed") == gluing.get("tests") == 25
)
check(
    "all three frozen artifacts and their provenance are exact",
    hashes_ok and provenance_ok,
    f"root={digest(ROOT_ARTIFACT)} tick={digest(TICK_ARTIFACT)} glue={digest(GLUING_ARTIFACT)}",
)

parities = ("even", "odd")
maps = {}
targets = {}
bounds = {}
map_ok = True
tick_ok = True
for parity in parities:
    mapping = tuple(gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"])
    first = tick["solutions"][parity]
    post = tuple(number(value) for value in first["post_momentum"])
    maps[parity] = mapping
    targets[parity] = tuple(post[index] for index in mapping)
    bounds[parity] = number(first["junction_bound"])
    map_ok &= len(mapping) == 30 and sorted(mapping) == list(range(30))
    tick_ok &= bool(
        first["converged"]
        and first["full_gate"]
        and len(post) == 30
        and bounds[parity] > 0
    )
check(
    "the independently derived maps construct two complete next-boundary targets",
    map_ok and tick_ok,
    f"bounds even/odd={text(bounds['even'])}/{text(bounds['odd'])}",
)

a1 = number(roots["geometry"]["a1"])
labels = []
label_ok = True
for index, root in enumerate(roots["roots"]):
    b = number(root["upper_log"])
    if b < a1:
        label = "CONTRACTING"
    elif root.get("kind") == "node" and abs(b) < arb.mpf("1e-25"):
        label = "TIME_REVERSAL"
    else:
        label = "UNCLASSIFIED"
    labels.append(label)
    label_ok &= len(root["pre_momentum"]) == 30
label_ok &= labels == ["CONTRACTING", "TIME_REVERSAL"]
check(
    "root labels follow geometry rather than target proximity",
    label_ok,
    f"labels={labels}, a1={text(a1)}",
)

comparisons = []
hit_indices = []
comparison_ok = True
for index, root in enumerate(roots["roots"]):
    pre = tuple(number(value) for value in root["pre_momentum"])
    parity_records = {}
    root_hit = True
    for parity in parities:
        residual = tuple(left-right for left, right in zip(pre, targets[parity]))
        residual_norm = norm(residual)
        maximum = max(abs(value) for value in residual)
        residual_spread = spread(residual)
        hit = residual_norm <= bounds[parity]
        root_hit &= hit
        comparison_ok &= bool(
            len(residual) == 30
            and arb.isfinite(residual_norm)
            and arb.isfinite(maximum)
            and arb.isfinite(residual_spread)
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
        hit_indices.append(index)
    comparisons.append({
        "root_index": index,
        "structural_label": labels[index],
        "upper_log": root["upper_log"],
        "pre_momentum": list(root["pre_momentum"]),
        "parities": parity_records,
        "hit_both_parities": bool(root_hit),
    })

check(
    "both committed roots are compared componentwise against both parity targets",
    comparison_ok and len(comparisons) == 2,
)

controls_ok = bool(hashes_ok and provenance_ok and map_ok and tick_ok and label_ok)
if not controls_ok:
    outcome = "SECOND_TICK_TARGET_CONTROL_FAILED"
elif hit_indices:
    outcome = "STATIONARY_SECOND_TICK_HIT"
else:
    outcome = "STATIONARY_SECOND_TICK_NO_HIT"

classification_ok = bool(
    (not controls_ok and outcome == "SECOND_TICK_TARGET_CONTROL_FAILED")
    or (controls_ok and hit_indices and outcome == "STATIONARY_SECOND_TICK_HIT")
    or (controls_ok and not hit_indices and outcome == "STATIONARY_SECOND_TICK_NO_HIT")
)
check(
    "the mechanical outcome follows the frozen hierarchy",
    classification_ok,
    f"outcome={outcome}, hits={len(hit_indices)}/2",
)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "control_correction_commit": CONTROL_CORRECTION_COMMIT,
    "input_sha256": {
        "stationary_roots": digest(ROOT_ARTIFACT),
        "first_tick": digest(TICK_ARTIFACT),
        "gluing": digest(GLUING_ARTIFACT),
    },
    "fixed_mass": True,
    "mass_recomputed_from_later_scale": False,
    "target_construction": "mapped accepted first-tick post momentum",
    "candidate_count": 2,
    "hit_count": len(hit_indices),
    "hit_fraction": f"{len(hit_indices)}/2",
    "hit_indices": hit_indices,
    "comparisons": comparisons,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

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

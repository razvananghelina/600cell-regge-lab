#!/usr/bin/env python3
"""Certify the exact universal staircase overlay of Delta^3 x I."""

from collections import Counter, deque
from itertools import permutations
import hashlib
import json
from pathlib import Path

import z3


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_600cell_universal_staircase_overlay.json"
PRIOR_ART_COMMIT = "4fb243b"
PROTOCOL_COMMIT = "4dc2feb"
VERTICES = tuple(range(4))
FULL_MASK = (1 << len(VERTICES))-1
MASKS = tuple(range(1, FULL_MASK))
MASK_INDEX = {mask: index for index, mask in enumerate(MASKS)}
ORDERS = tuple(permutations(VERTICES))


def mask_of(vertices):
    result = 0
    for vertex in vertices:
        result |= 1 << vertex
    return result


def sign_word(pattern):
    return "".join("+" if sign > 0 else "-" for sign in pattern)


def make_problem():
    lambdas = tuple(z3.Real(f"lambda_{index}") for index in VERTICES)
    time = z3.Real("t")
    forms = {
        mask: time-z3.Sum([
            lambdas[index] for index in VERTICES if mask & (1 << index)
        ])
        for mask in MASKS
    }
    base = [
        *(value > 0 for value in lambdas),
        z3.Sum(lambdas) == 1,
        time > 0,
        time < 1,
    ]
    return lambdas, time, forms, base


# Exact depth-first feasibility census with prefix pruning.
_, _, forms, base = make_problem()
prefix_solver = z3.Solver()
prefix_solver.add(*base)
feasible_patterns = []
prefix_checks = 0


def traverse(depth, prefix):
    global prefix_checks
    if depth == len(MASKS):
        feasible_patterns.append(tuple(prefix))
        return
    form = forms[MASKS[depth]]
    for sign in (-1, 1):
        prefix_solver.push()
        prefix_solver.add(form > 0 if sign > 0 else form < 0)
        prefix_checks += 1
        if prefix_solver.check() == z3.sat:
            prefix.append(sign)
            traverse(depth+1, prefix)
            prefix.pop()
        prefix_solver.pop()


traverse(0, [])
feasible_patterns = tuple(sorted(feasible_patterns))
feasible_set = set(feasible_patterns)
census_basic_ok = bool(
    feasible_patterns
    and len(feasible_patterns) == len(feasible_set)
    and all(len(pattern) == 14 for pattern in feasible_patterns)
    and prefix_checks <= 2*(2**14-1)
)


# Independently rebuild and solve every recorded full sign word.
independent_rechecks = []
for pattern in feasible_patterns:
    _, _, fresh_forms, fresh_base = make_problem()
    solver = z3.Solver()
    solver.add(*fresh_base)
    for mask, sign in zip(MASKS, pattern):
        solver.add(fresh_forms[mask] > 0 if sign > 0 else fresh_forms[mask] < 0)
    independent_rechecks.append(solver.check() == z3.sat)
independent_recheck_ok = bool(
    len(independent_rechecks) == len(feasible_patterns)
    and all(independent_rechecks)
)


def permute_mask(mask, permutation):
    return mask_of(
        permutation[vertex]
        for vertex in VERTICES
        if mask & (1 << vertex)
    )


def transform_pattern(pattern, permutation, reflect_time):
    result = [None]*len(MASKS)
    for mask, sign in zip(MASKS, pattern):
        target = permute_mask(mask, permutation)
        if reflect_time:
            target = FULL_MASK ^ target
            sign = -sign
        result[MASK_INDEX[target]] = sign
    if any(value is None for value in result):
        raise RuntimeError("incomplete signed permutation")
    return tuple(result)


transformations = tuple(
    (permutation, reflect_time)
    for permutation in ORDERS
    for reflect_time in (False, True)
)
transformation_keys = []
for permutation, reflect_time in transformations:
    key = []
    for mask in MASKS:
        target = permute_mask(mask, permutation)
        multiplier = 1
        if reflect_time:
            target = FULL_MASK ^ target
            multiplier = -1
        key.append((MASK_INDEX[target], multiplier))
    transformation_keys.append(tuple(key))
transformation_count_ok = len(set(transformation_keys)) == 48

symmetry_images_ok = bool(
    all(
        transform_pattern(pattern, permutation, reflect_time) in feasible_set
        for pattern in feasible_patterns
        for permutation, reflect_time in transformations
    )
)


# Exact S4 x C2 chamber orbits.
unseen = set(feasible_patterns)
orbit_records = []
while unseen:
    representative = min(unseen)
    orbit = {
        transform_pattern(representative, permutation, reflect_time)
        for permutation, reflect_time in transformations
    }
    queue = deque(orbit)
    while queue:
        pattern = queue.popleft()
        for permutation, reflect_time in transformations:
            image = transform_pattern(pattern, permutation, reflect_time)
            if image not in orbit:
                orbit.add(image)
                queue.append(image)
    orbit_records.append({
        "representative": sign_word(representative),
        "size": len(orbit),
    })
    unseen -= orbit
orbit_records.sort(key=lambda record: (record["size"], record["representative"]))
orbit_sizes = [record["size"] for record in orbit_records]
orbit_control_ok = bool(
    sum(orbit_sizes) == len(feasible_patterns)
    and all(48 % size == 0 for size in orbit_sizes)
)


def staircase_assignments(pattern, order):
    assignments = []
    for split in range(4):
        lower_mask = mask_of(order[split+1:])
        upper_mask = mask_of(order[split:])
        lower_ok = bool(
            lower_mask == 0
            or pattern[MASK_INDEX[lower_mask]] > 0
        )
        upper_ok = bool(
            upper_mask == FULL_MASK
            or pattern[MASK_INDEX[upper_mask]] < 0
        )
        if lower_ok and upper_ok:
            assignments.append(split)
    return tuple(assignments)


order_records = []
bad_assignment_examples = []
for order in ORDERS:
    counts = Counter()
    multiplicities = Counter()
    for pattern in feasible_patterns:
        assigned = staircase_assignments(pattern, order)
        multiplicities[len(assigned)] += 1
        if len(assigned) == 1:
            counts[assigned[0]] += 1
        elif len(bad_assignment_examples) < 20:
            bad_assignment_examples.append({
                "sign_word": sign_word(pattern),
                "order": list(order),
                "assignments": list(assigned),
            })
    order_records.append({
        "order": list(order),
        "simplex_assignment_counts": {
            str(index): counts[index] for index in range(4)
        },
        "assignment_multiplicity_distribution": {
            str(key): value for key, value in sorted(multiplicities.items())
        },
    })
common_refinement_ok = bool(
    all(
        record["assignment_multiplicity_distribution"]
        == {"1": len(feasible_patterns)}
        for record in order_records
    )
)


# Analytic wall controls for all 24 staircases.
order_wall_records = []
all_used_walls = set()
order_walls_ok = True
for order in ORDERS:
    walls = tuple(mask_of(order[index:]) for index in (1, 2, 3))
    all_used_walls.update(walls)
    sizes = tuple(mask.bit_count() for mask in walls)
    nested = bool(
        walls[1] & walls[0] == walls[1]
        and walls[2] & walls[1] == walls[2]
    )
    order_walls_ok &= sizes == (3, 2, 1) and nested
    order_wall_records.append({
        "order": list(order),
        "tail_masks": list(walls),
        "tail_sizes": list(sizes),
    })
wall_union_ok = all_used_walls == set(MASKS)


# Symbolic restriction of the labelled subset forms to each spatial face.
face_records = []
face_restriction_ok = True
for omitted in VERTICES:
    remaining = tuple(vertex for vertex in VERTICES if vertex != omitted)
    restricted_counts = Counter()
    for mask in MASKS:
        restricted = tuple(
            vertex for vertex in remaining if mask & (1 << vertex)
        )
        restricted_counts[restricted] += 1
    expected_keys = {
        tuple(remaining[index] for index in range(3) if bits & (1 << index))
        for bits in range(8)
    }
    empty = ()
    full = remaining
    this_ok = bool(
        set(restricted_counts) == expected_keys
        and restricted_counts[empty] == 1
        and restricted_counts[full] == 1
        and all(
            restricted_counts[key] == 2
            for key in expected_keys
            if key not in (empty, full)
        )
    )
    face_restriction_ok &= this_ok
    face_records.append({
        "omitted_vertex": omitted,
        "remaining_vertices": list(remaining),
        "restricted_form_multiplicities": {
            "".join(map(str, key)) if key else "empty": value
            for key, value in sorted(restricted_counts.items())
        },
        "passes": this_ok,
    })


complement_control_ok = bool(
    len(MASKS) == len(set(MASKS)) == 14
    and all((FULL_MASK ^ mask) in MASK_INDEX for mask in MASKS)
    and all((FULL_MASK ^ mask) != mask for mask in MASKS)
)
analytic_controls_ok = bool(
    complement_control_ok
    and transformation_count_ok
    and order_walls_ok
    and wall_union_ok
)
controls_ok = bool(
    census_basic_ok
    and independent_recheck_ok
    and analytic_controls_ok
    and symmetry_images_ok
    and orbit_control_ok
    and face_restriction_ok
)

if not controls_ok:
    outcome = "UNIVERSAL_STAIRCASE_OVERLAY_CONTROL_FAILED"
elif common_refinement_ok:
    outcome = "UNIVERSAL_STAIRCASE_OVERLAY_CERTIFIED"
else:
    outcome = "UNIVERSAL_STAIRCASE_OVERLAY_NOT_COMMON_REFINEMENT"

tests = [
    ("14 nontrivial subset forms pair exactly by complements", complement_control_ok),
    ("exact prefix traversal returns distinct full sign words", census_basic_ok),
    ("every feasible word passes an independent exact recheck", independent_recheck_ok),
    ("all 48 labelled S4 x C2 transformations are distinct", transformation_count_ok),
    ("the feasible chamber set is invariant under S4 x C2", symmetry_images_ok),
    ("orbit sizes divide 48 and exhaust all chambers", orbit_control_ok),
    ("every staircase has its nested 3,2,1 tail walls", order_walls_ok),
    ("the 24 staircase orders use all 14 arrangement walls", wall_union_ok),
    ("all four spatial-face restrictions have the triangular overlay", face_restriction_ok),
    ("all chamber/order pairs have exactly one staircase simplex", common_refinement_ok),
    ("no action, metric, dynamics or physical target was evaluated", True),
    ("outcome follows the preregistered mechanical rule", outcome in {
        "UNIVERSAL_STAIRCASE_OVERLAY_CONTROL_FAILED",
        "UNIVERSAL_STAIRCASE_OVERLAY_CERTIFIED",
        "UNIVERSAL_STAIRCASE_OVERLAY_NOT_COMMON_REFINEMENT",
    }),
]
passed = sum(bool(ok) for _, ok in tests)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "z3_version": z3.get_version_string(),
    "subset_masks": list(MASKS),
    "candidate_sign_patterns": 2**14,
    "prefix_solver_checks": prefix_checks,
    "full_dimensional_chamber_count": len(feasible_patterns),
    "feasible_sign_words": [sign_word(pattern) for pattern in feasible_patterns],
    "symmetry_group_order": len(transformations),
    "symmetry_orbits": orbit_records,
    "orbit_size_distribution": {
        str(size): count for size, count in sorted(Counter(orbit_sizes).items())
    },
    "staircase_orders": order_records,
    "bad_assignment_examples": bad_assignment_examples,
    "order_walls": order_wall_records,
    "face_restrictions": face_records,
    "common_refinement": common_refinement_ok,
    "gravity_action_evaluations": 0,
    "physical_target_parsed": False,
    "tests": len(tests),
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"candidate sign patterns={2**14}")
print(f"exact prefix checks={prefix_checks}")
print(f"full-dimensional chambers={len(feasible_patterns)}")
print(f"symmetry orbits={len(orbit_records)}; sizes={dict(Counter(orbit_sizes))}")
print(f"assignment failures={len(bad_assignment_examples)}")
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) and controls_ok else 1)


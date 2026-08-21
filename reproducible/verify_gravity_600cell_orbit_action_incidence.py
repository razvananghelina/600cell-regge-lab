#!/usr/bin/env python3
"""Exact flag-incidence audit of the reduced 600-cell slab action.

Prior-art commit: c14b5ac.
Protocol commit: c90833d.
Frozen direct/orbit disagreement commit: af862ab.

No existing action evaluator is modified.  Exact finite incidence data are
constructed before any stored action is compared.
"""

import ast
from collections import Counter
import contextlib
from fractions import Fraction
import hashlib
from itertools import combinations
import io
import json
from pathlib import Path
import re

import mpmath as mp


HERE = Path(__file__).resolve().parent
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_dust_canonical_continuation.py"
DIRECT_SOURCE = HERE / "verify_gravity_global_boundary_legendre.py"
PRIMARY_ARTIFACT = HERE / "gravity_600cell_tick_scale_covariance.json"
DIRECT_ARTIFACT = HERE / "gravity_600cell_tick_scale_covariance_precision.json"
REGULAR_ARTIFACT = HERE / "gravity_600cell_published_dust_control.json"
OUTPUT = HERE / "gravity_600cell_orbit_action_incidence.json"
PRIOR_ART_COMMIT = "c14b5ac"
PROTOCOL_COMMIT = "c90833d"
DISAGREEMENT_COMMIT = "af862ab"
DPS = 80
mp.mp.dps = DPS
ACTION_TOLERANCE = mp.mpf("1e-45")
IMAGINARY_TOLERANCE = mp.mpf("1e-60")
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


def text(value, digits=35):
    return mp.nstr(value, digits)


def relative_error(left, right):
    return abs(left - right) / max(mp.mpf(1), abs(left), abs(right))


def definition_prefix(path, cutoff_target):
    tree = ast.parse(path.read_text(), filename=str(path))
    cut = None
    for index, node in enumerate(tree.body):
        if cutoff_target == "first_print":
            if (
                isinstance(node, ast.Expr)
                and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name)
                and node.value.func.id == "print"
            ):
                cut = index
                break
        elif (
            isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == cutoff_target
                for target in node.targets
            )
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError(f"definition cutoff {cutoff_target} absent in {path}")
    namespace = {
        "__file__": str(path),
        "__name__": f"orbit_incidence_{path.stem}",
    }
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(path), "exec"), namespace)
    return namespace


print("Exact flag-incidence audit of the reduced slab action", flush=True)
primary_core = definition_prefix(PRIMARY_SOURCE, "first_print")
direct_core = definition_prefix(DIRECT_SOURCE, "CONTROLS")
primary_models = primary_core["models"]
direct_models = direct_core["models"]
provenance_ok = bool(
    PRIOR_ART_COMMIT == "c14b5ac"
    and PROTOCOL_COMMIT == "c90833d"
    and DISAGREEMENT_COMMIT == "af862ab"
    and primary_core["tests"] == primary_core["passed"] == 4
    and direct_core["tests"] == direct_core["passed"] == 7
    and set(primary_models) == set(direct_models) == {"even", "odd"}
)
check("both frozen carriers and the audit provenance are intact", provenance_ok)


def action_key(action):
    return tuple(int(value) for value in action)


def image_vertex(action, vertex):
    return int(action[vertex] if vertex < 120 else action[vertex - 120] + 120)


def image_cell(action, cell):
    return tuple(sorted(image_vertex(action, vertex) for vertex in cell))


def image_flag(action, flag):
    triangle, simplex = flag
    return image_cell(action, triangle), image_cell(action, simplex)


def partition(items, stabilizer, image):
    unseen = set(items)
    result = []
    while unseen:
        seed = min(unseen)
        orbit = frozenset(image(action, seed) for action in stabilizer)
        if not orbit <= set(items):
            raise RuntimeError("group image left the frozen finite set")
        result.append(orbit)
        unseen -= orbit
    return tuple(sorted(result, key=lambda orbit: (len(orbit), min(orbit))))


def normalized_stabilizer(model):
    return tuple(sorted(action_key(action) for action in model["stabilizer"]))


def triangle_set(model):
    return frozenset(
        triangle for orbit in model["triangle_orbits"] for triangle in orbit
    )


carrier_ok = True
carrier_details = {}
for parity in primary_models:
    left = primary_models[parity]
    right = direct_models[parity]
    final_left = dict(left["final_lookup"])
    final_right = {
        edge: variable - 35 for edge, variable in right["final_to_variable"].items()
    }
    one_ok = bool(
        left["slab"] == right["slab"]
        and normalized_stabilizer(left) == normalized_stabilizer(right)
        and left["old_edges"] == right["old_edges"]
        and left["internal_edges"] == right["internal_edges"]
        and left["new_edges"] == right["new_edges"]
        and triangle_set(left) == triangle_set(right)
        and left["boundary_triangles"] == right["boundary_triangles"]
        and left["old_lookup"] == right["old_to_orbit"]
        and left["edge_to_variable"] == right["edge_to_variable"]
        and left["edge_jacobian"] == right["edge_jacobian"]
        and final_left == final_right
    )
    carrier_ok &= one_ok
    carrier_details[parity] = {
        "slab_simplices": len(left["slab"]),
        "stabilizer_order": len(left["stabilizer"]),
        "triangles": len(triangle_set(left)),
        "maps_equal": one_ok,
    }
check(
    "the primary and direct carriers, coordinate maps and edge signs agree exactly",
    carrier_ok,
)


# Frozen states, reconstructed before incidence or artifact comparisons.
M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON_3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON_3 * L0
TAU = mp.mpf("0.0102")
RHO = TAU**2
SLANT_SQUARE = L0_SQUARE - RHO
I = mp.mpc(0, 1)


def perturb(values, modulus, center):
    return tuple(
        value * mp.exp(mp.mpf("1e-6") * ((index % modulus) - center))
        for index, value in enumerate(values)
    )


states = {
    "regular": {
        "old": tuple([L0_SQUARE] * 30),
        "variables": tuple([SLANT_SQUARE] * 30 + [RHO] * 5 + [L0_SQUARE] * 30),
    },
    "off_shell": {
        "old": perturb(tuple([L0_SQUARE] * 30), 7, 3),
        "variables": perturb(
            tuple([SLANT_SQUARE] * 30 + [RHO] * 5 + [L0_SQUARE] * 30),
            5,
            2,
        ),
    },
}


def edge_value(model, edge, state):
    edge = tuple(sorted(edge))
    if edge in model["old_to_orbit"]:
        return state["old"][model["old_to_orbit"][edge]]
    if edge in model["all_edge_to_variable"]:
        index = model["all_edge_to_variable"][edge]
        return mp.mpf(model["all_edge_jacobian"][edge]) * state["variables"][index]
    raise ValueError(f"edge missing from state expansion: {edge}")


invariance_ok = True
for parity, model in direct_models.items():
    all_edges = model["old_edges"] | model["internal_edges"] | model["new_edges"]
    for state in states.values():
        for action in model["stabilizer"]:
            for edge in all_edges:
                mapped = image_cell(action, edge)
                invariance_ok &= bool(
                    edge_value(model, edge, state) == edge_value(model, mapped, state)
                )
check(
    "the regular and frozen off-shell edge states are exactly stabilizer-invariant",
    invariance_ok,
)


enumerations = {}
enumeration_ok = True
shortcut_differs = False
for parity, model in direct_models.items():
    stabilizer = model["stabilizer"]
    triangles = triangle_set(model)
    simplices = frozenset(model["slab"])
    flags = frozenset(
        (tuple(triangle), tuple(simplex))
        for simplex in simplices
        for triangle in combinations(simplex, 3)
    )
    triangle_orbits = partition(triangles, stabilizer, image_cell)
    simplex_orbits = partition(simplices, stabilizer, image_cell)
    flag_orbits = partition(flags, stabilizer, image_flag)
    triangle_lookup = {
        item: index for index, orbit in enumerate(triangle_orbits) for item in orbit
    }
    simplex_lookup = {
        item: index for index, orbit in enumerate(simplex_orbits) for item in orbit
    }
    flag_lookup = {
        item: index for index, orbit in enumerate(flag_orbits) for item in orbit
    }
    records = []
    coefficients = []
    row_count_ok = True
    for flag_index, orbit in enumerate(flag_orbits):
        representative = min(orbit)
        triangle_index = triangle_lookup[representative[0]]
        simplex_index = simplex_lookup[representative[1]]
        triangle_orbit = triangle_orbits[triangle_index]
        divisible = len(orbit) % len(triangle_orbit) == 0
        coefficient = Fraction(len(orbit), len(triangle_orbit))
        coefficients.append(coefficient)
        incident_counts = []
        for triangle in triangle_orbit:
            incident_counts.append(sum(flag[0] == triangle for flag in orbit))
        row_count_ok &= bool(
            divisible
            and coefficient.denominator == 1
            and all(count == coefficient.numerator for count in incident_counts)
        )
        records.append({
            "flag_orbit": flag_index,
            "representative_triangle": list(representative[0]),
            "representative_simplex": list(representative[1]),
            "flag_orbit_size": len(orbit),
            "triangle_orbit": triangle_index,
            "triangle_orbit_size": len(triangle_orbit),
            "simplex_orbit": simplex_index,
            "simplex_orbit_size": len(simplex_orbits[simplex_index]),
            "coefficient": coefficient.numerator,
        })
    shortcut_counts = Counter()
    for simplex_orbit in model["simplex_orbits"]:
        simplex = min(simplex_orbit)
        for triangle in combinations(simplex, 3):
            shortcut_counts[flag_lookup[(tuple(triangle), tuple(simplex))]] += 1
    mismatches = []
    for index, coefficient in enumerate(coefficients):
        shortcut = shortcut_counts[index]
        if shortcut != coefficient:
            mismatches.append({
                "flag_orbit": index,
                "correct": coefficient.numerator,
                "shortcut": shortcut,
                "triangle_orbit": records[index]["triangle_orbit"],
                "simplex_orbit": records[index]["simplex_orbit"],
            })
    shortcut_differs |= bool(mismatches)
    stored_triangle_orbits = {
        frozenset(orbit) for orbit in model["triangle_orbits"]
    }
    stored_simplex_orbits = {
        frozenset(orbit) for orbit in model["simplex_orbits"]
    }
    counting_ok = bool(
        len(flags) == 24000
        and sum(map(len, flag_orbits)) == 24000
        and sum(
            coefficient.numerator * len(triangle_orbits[records[index]["triangle_orbit"]])
            for index, coefficient in enumerate(coefficients)
        ) == 24000
        and sum(shortcut_counts.values()) == 10 * len(simplex_orbits)
        and {frozenset(orbit) for orbit in triangle_orbits} == stored_triangle_orbits
        and {frozenset(orbit) for orbit in simplex_orbits} == stored_simplex_orbits
    )
    one_ok = bool(row_count_ok and counting_ok)
    enumeration_ok &= one_ok
    payload = {
        "parity": parity,
        "triangle_orbit_sizes": sorted(map(len, triangle_orbits)),
        "simplex_orbit_sizes": sorted(map(len, simplex_orbits)),
        "flag_orbit_sizes": sorted(map(len, flag_orbits)),
        "coefficient_multiset": {
            str(key): value for key, value in sorted(
                Counter(coefficient.numerator for coefficient in coefficients).items()
            )
        },
        "records": records,
        "mismatches": mismatches,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    enumerations[parity] = {
        "triangles": triangles,
        "simplices": simplices,
        "flags": flags,
        "triangle_orbits": triangle_orbits,
        "simplex_orbits": simplex_orbits,
        "flag_orbits": flag_orbits,
        "triangle_lookup": triangle_lookup,
        "flag_lookup": flag_lookup,
        "coefficients": tuple(coefficients),
        "shortcut_counts": shortcut_counts,
        "payload": payload,
        "sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "ok": one_ok,
    }
    check(
        f"{parity}: exact triangle, simplex and 24000-flag orbit counts close",
        one_ok,
        f"triangle orbits={len(triangle_orbits)}, simplex orbits={len(simplex_orbits)}, flag orbits={len(flag_orbits)}, mismatches={len(mismatches)}",
    )


def signed_volume_square(squared, vertices):
    vertices = list(vertices)
    dimension = len(vertices) - 1
    if dimension == 0:
        return mp.mpf(1)
    base = vertices[0]
    others = vertices[1:]
    gram = mp.matrix([
        [
            (squared[base][left] + squared[base][right] - squared[left][right]) / 2
            for right in others
        ]
        for left in others
    ])
    return mp.det(gram) / mp.factorial(dimension) ** 2


def simplex_squared(model, simplex, state):
    squared = [[mp.mpf(0) for _ in range(5)] for _ in range(5)]
    for left, right in combinations(range(5), 2):
        value = edge_value(model, (simplex[left], simplex[right]), state)
        squared[left][right] = squared[right][left] = value
    return squared


def log_minus(value):
    scale = max(mp.mpf(1), abs(value))
    if abs(mp.im(value)) < mp.mpf("1e-65") * scale:
        real = mp.re(value)
        if real < 0:
            return mp.log(-real) - I * mp.pi
        return mp.log(real)
    return mp.log(value)


def angle_record(squared):
    gram = mp.matrix([
        [
            (squared[0][left] + squared[0][right] - squared[left][right]) / 2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    inverse = gram**-1
    simplex_volume_square = signed_volume_square(squared, range(5))
    facet_volume_squares = {
        omitted: signed_volume_square(
            squared, [vertex for vertex in range(5) if vertex != omitted]
        )
        for omitted in range(5)
    }
    leading_minors = []
    for size in range(1, 5):
        principal = mp.matrix([
            [gram[left, right] for right in range(size)] for left in range(size)
        ])
        leading_minors.append(mp.det(principal))
    signs = [1] + [1 if value > 0 else -1 if value < 0 else 0 for value in leading_minors]
    negative_directions = None if 0 in signs else sum(
        left != right for left, right in zip(signs, signs[1:])
    )
    angles = {}
    minimum_argument = mp.inf
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge = tuple(
            vertex for vertex in range(5)
            if vertex not in (omitted_a, omitted_b)
        )
        hinge_volume_square = signed_volume_square(squared, hinge)
        derivative = mp.matrix(4, 4)
        opposite = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                derivative[left - 1, right - 1] = (
                    int({0, left} == opposite)
                    + int({0, right} == opposite)
                    - int(left != right and {left, right} == opposite)
                ) / 2
        product = inverse * derivative
        volume_derivative = simplex_volume_square * sum(
            product[index, index] for index in range(4)
        )
        denominator = (
            mp.sqrt(mp.mpc(facet_volume_squares[omitted_a]))
            * mp.sqrt(mp.mpc(facet_volume_squares[omitted_b]))
        )
        cosine = 16 * volume_derivative / denominator
        sine = -mp.mpf(4) / 3 * (
            mp.sqrt(mp.mpc(hinge_volume_square))
            * mp.sqrt(mp.mpc(simplex_volume_square))
        ) / denominator
        argument = cosine + I * sine
        minimum_argument = min(minimum_argument, abs(argument))
        angles[hinge] = -I * log_minus(argument)
    return angles, negative_directions, min(abs(value) for value in leading_minors), minimum_argument


def triangle_area_square(values):
    x, y, z = values
    return (2 * (x * y + x * z + y * z) - x**2 - y**2 - z**2) / 16


def reduced_actions(model, enumeration, state):
    triangle_orbits = enumeration["triangle_orbits"]
    flag_orbits = enumeration["flag_orbits"]
    triangle_lookup = enumeration["triangle_lookup"]
    exact = [
        mp.pi if min(orbit) in model["boundary_triangles"] else 2 * mp.pi
        for orbit in triangle_orbits
    ]
    shortcut = list(exact)
    cache = {}
    negative_counts = Counter()
    minimum_minor = mp.inf
    minimum_argument = mp.inf

    def angles_for(simplex):
        nonlocal minimum_minor, minimum_argument
        if simplex not in cache:
            angles, negative, minor, argument = angle_record(
                simplex_squared(model, simplex, state)
            )
            cache[simplex] = angles
            negative_counts[negative] += 1
            minimum_minor = min(minimum_minor, minor)
            minimum_argument = min(minimum_argument, argument)
        return cache[simplex]

    for flag_index, flag_orbit in enumerate(flag_orbits):
        triangle, simplex = min(flag_orbit)
        local_hinge = tuple(
            index for index, vertex in enumerate(simplex) if vertex in triangle
        )
        angle = angles_for(simplex)[local_hinge]
        triangle_index = triangle_lookup[triangle]
        exact[triangle_index] += (
            enumeration["coefficients"][flag_index].numerator * angle
        )
    for simplex_orbit in model["simplex_orbits"]:
        simplex = min(simplex_orbit)
        for local_hinge, angle in angles_for(simplex).items():
            triangle = tuple(sorted(simplex[index] for index in local_hinge))
            shortcut[triangle_lookup[triangle]] += angle

    actions = []
    for curvatures in (exact, shortcut):
        action_sum = mp.mpc(0)
        for index, orbit in enumerate(triangle_orbits):
            triangle = min(orbit)
            values = tuple(
                edge_value(model, edge, state)
                for edge in combinations(triangle, 2)
            )
            area = mp.sqrt(mp.mpc(triangle_area_square(values)))
            action_sum += len(orbit) * area * curvatures[index]
        gravitational = -I * action_sum
        dust = -(8 * mp.pi * MASS / 5) * sum(
            mp.sqrt(value) for value in state["variables"][30:35]
        )
        actions.append(gravitational + dust)
    branch_ok = bool(
        sum(negative_counts.values()) == len(cache)
        and negative_counts == Counter({1: len(cache)})
        and minimum_minor > mp.mpf("1e-20")
        and minimum_argument > mp.mpf("1e-6")
        and max(abs(mp.im(action)) for action in actions) < IMAGINARY_TOLERANCE
    )
    return {
        "flag": actions[0],
        "shortcut": actions[1],
        "evaluated_simplices": len(cache),
        "negative_counts": negative_counts,
        "minimum_minor": minimum_minor,
        "minimum_argument": minimum_argument,
        "branch_ok": branch_ok,
    }


# Finish all eight reduced values before loading any action artifact.
action_records = {}
action_branch_ok = True
for parity, model in direct_models.items():
    action_records[parity] = {}
    for name, state in states.items():
        print(f"  evaluating {parity}/{name}", flush=True)
        record = reduced_actions(model, enumerations[parity], state)
        action_records[parity][name] = record
        action_branch_ok &= record["branch_ok"]
check("all flag and shortcut representative evaluations retain the Lorentzian branch", action_branch_ok)


def parse_mpc(value):
    match = re.fullmatch(r"\((\S+)\s+([+-])\s+(\S+)j\)", value.strip())
    if match is None:
        raise ValueError(f"cannot parse mpc text: {value}")
    imaginary = mp.mpf(match.group(3))
    if match.group(2) == "-":
        imaginary = -imaginary
    return mp.mpc(mp.mpf(match.group(1)), imaginary)


primary = json.loads(PRIMARY_ARTIFACT.read_text())
direct = json.loads(DIRECT_ARTIFACT.read_text())
regular = json.loads(REGULAR_ARTIFACT.read_text())
artifact_ok = bool(
    primary.get("outcome") == "TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED"
    and primary.get("tests") == primary.get("passed") == 12
    and direct.get("outcome") == "TICK_SCALE_COVARIANCE_HIGH_PRECISION_DISAGREES"
    and direct.get("tests") == 10
    and direct.get("passed") == 7
    and regular.get("verdict")
        == (
            "DERIVED EXTERNAL CONTROL: the published time-symmetric dust "
            "sandwich solves all 35 complete one-slab orbit equations in both "
            "parities."
        )
    and regular.get("tests") == regular.get("passed") == 14
)

discriminator_ok = True
regular_ok = True
comparison_records = {}
for parity in direct_models:
    flag = action_records[parity]["off_shell"]["flag"]
    shortcut = action_records[parity]["off_shell"]["shortcut"]
    direct_target = parse_mpc(
        direct["parities"][parity]["states"]["1.0"]["action"]
    )
    shortcut_target = parse_mpc(primary["parities"][parity]["base"]["action"])
    flag_error = relative_error(flag, direct_target)
    shortcut_error = relative_error(shortcut, shortcut_target)
    disagreement_error = relative_error(
        shortcut - flag, shortcut_target - direct_target
    )
    one_discriminator_ok = bool(
        flag_error < ACTION_TOLERANCE
        and shortcut_error < ACTION_TOLERANCE
        and disagreement_error < ACTION_TOLERANCE
    )
    discriminator_ok &= one_discriminator_ok

    regular_flag = action_records[parity]["regular"]["flag"]
    regular_shortcut = action_records[parity]["regular"]["shortcut"]
    regular_target = mp.mpc(
        regular["parities"][parity]["arbitrary_precision_base_action_real"],
        regular["parities"][parity]["arbitrary_precision_base_action_imaginary"],
    )
    regular_equal_error = relative_error(regular_flag, regular_shortcut)
    regular_flag_error = relative_error(regular_flag, regular_target)
    regular_shortcut_error = relative_error(regular_shortcut, regular_target)
    one_regular_ok = bool(
        regular_equal_error < ACTION_TOLERANCE
        and regular_flag_error < ACTION_TOLERANCE
        and regular_shortcut_error < ACTION_TOLERANCE
    )
    regular_ok &= one_regular_ok
    comparison_records[parity] = {
        "off_shell": {
            "flag_action": text(flag, 65),
            "shortcut_action": text(shortcut, 65),
            "direct_target_error": text(flag_error, 15),
            "primary_target_error": text(shortcut_error, 15),
            "disagreement_reproduction_error": text(disagreement_error, 15),
            "pass": one_discriminator_ok,
        },
        "regular": {
            "flag_action": text(regular_flag, 65),
            "shortcut_action": text(regular_shortcut, 65),
            "flag_shortcut_error": text(regular_equal_error, 15),
            "flag_target_error": text(regular_flag_error, 15),
            "shortcut_target_error": text(regular_shortcut_error, 15),
            "pass": one_regular_ok,
        },
    }
    check(
        f"{parity}: exact flags reproduce the direct action and the shortcut reproduces the primary action",
        one_discriminator_ok,
        f"flag/direct={text(flag_error, 6)}, shortcut/primary={text(shortcut_error, 6)}",
    )
    check(
        f"{parity}: flag and shortcut reductions retain the regular published control",
        one_regular_ok,
        f"flag/shortcut={text(regular_equal_error, 6)}, flag/control={text(regular_flag_error, 6)}",
    )

check("all three frozen action artifacts have the required outcomes", artifact_ok)

controls_ok = bool(
    provenance_ok
    and carrier_ok
    and invariance_ok
    and enumeration_ok
    and action_branch_ok
    and artifact_ok
    and regular_ok
)
if not controls_ok:
    outcome = "ORBIT_ACTION_INCIDENCE_CONTROL_FAILED"
elif not shortcut_differs and not discriminator_ok:
    outcome = "ORBIT_ACTION_DIRECT_CONSTRUCTION_SUSPECT"
elif shortcut_differs and discriminator_ok:
    outcome = "ORBIT_ACTION_FLAG_MULTIPLICITY_BUG_DERIVED"
else:
    outcome = "ORBIT_ACTION_DISAGREEMENT_UNRESOLVED"

check(
    "the exact incidence audit identifies the source of the action disagreement",
    outcome == "ORBIT_ACTION_FLAG_MULTIPLICITY_BUG_DERIVED",
    outcome,
)

artifact = {
    "title": "Exact flag-incidence audit of the reduced 600-cell slab action",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "frozen_disagreement_commit": DISAGREEMENT_COMMIT,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "precision_digits": DPS,
    "carrier": carrier_details,
    "state_invariance": invariance_ok,
    "enumerations": {
        parity: {
            **data["payload"],
            "sha256": data["sha256"],
        }
        for parity, data in enumerations.items()
    },
    "actions": comparison_records,
    "interpretation": {
        "derived_if_outcome_three": (
            "the existing representative-simplex orbit reduction uses incorrect "
            "flag-incidence multiplicities off the regular symmetric locus"
        ),
        "scope_audit_required": True,
        "not_impugned": [
            "direct complete Regge action",
            "exact global scale homogeneity",
        ],
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

print(f"RESULT: {passed}/{tests}", flush=True)
print(f"OUTCOME: {outcome}", flush=True)
if passed != tests:
    raise SystemExit(1)

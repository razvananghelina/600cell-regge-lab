#!/usr/bin/env python3
"""Multi-level action diagnostic for the constrained H4 response control."""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
CURVATURE = HERE / "gravity_600cell_refined_local_curvature_mass.json"
PRIMARY = HERE / "gravity_600cell_refined_h4_constrained_response.json"
PRIMARY_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_primary_first_result.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_directional_diagnostic_protocol.md"
)
CELL600 = ROOT / "commons/cell600.py"
OUTPUT = HERE / "gravity_600cell_refined_h4_constrained_directional_diagnostic.json"

PROTOCOL_COMMIT = "0baeeee"
EXPECTED_HASHES = {
    "action_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "curvature": "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "primary": "f029260c9ee6e3b763293d237aae27e6ff7c1256eb8bc19c35725084ff385888",
    "primary_result": "633a57f3d2b4a054cce20d08544d409dac8fdaf53c39bae72ab2e9fceb4e83eb",
    "protocol": "6dcaabdf4aee95a32751a52b481a0c1d61e33e0f08a8a67e0755a9cdfd932606",
    "cell600": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}

PAIR4 = tuple(combinations(range(4), 2))
BOUNDARY_VARIABLES = (
    tuple(("old",) + pair for pair in PAIR4)
    + tuple(("new",) + pair for pair in PAIR4)
)
INTERNAL_VARIABLES = (
    tuple(("cross",) + pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
VARIABLES = BOUNDARY_VARIABLES + INTERNAL_VARIABLES
LOCAL_TRIANGLES = np.asarray(tuple(combinations(range(5), 3)), dtype=np.int8)
TAU_TEXT = "0.0102"
DPS_LEVELS = (140, 180)
STEP_TEXTS = ("1e-10", "5e-11", "2.5e-11", "1.25e-11", "6.25e-12")
DIRECTIONAL_INDICES = (0, 1, 22, 23)

tests = 0
passed = 0


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
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=70):
    return mp.nstr(value, digits, strip_zeros=False)


def parse_mpmath_complex(text):
    body = text.strip()
    if not (body.startswith("(") and body.endswith(")")):
        raise ValueError("stored complex value lacks parentheses")
    real_text, sign, imaginary_text = body[1:-1].rsplit(" ", 2)
    if sign not in {"+", "-"} or not imaginary_text.endswith("j"):
        raise ValueError("stored complex value has an unexpected form")
    imaginary = mp.mpf(imaginary_text[:-1])
    if sign == "-":
        imaginary = -imaginary
    return mp.mpc(mp.mpf(real_text), imaginary)


def load_action_definitions():
    tree = ast.parse(ACTION_SOURCE.read_text(), filename=str(ACTION_SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "mp": mp,
        "np": np,
        "json": json,
        "Path": Path,
        "sha256": sha256,
        "combinations": combinations,
        "permutations": permutations,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "HERE": HERE,
        "ROOT": ROOT,
        "PAIR4": PAIR4,
        "LOCAL_TRIANGLES": LOCAL_TRIANGLES,
        "TAU_TEXT": TAU_TEXT,
        "VARIABLES": VARIABLES,
        "INTERNAL_VARIABLES": INTERNAL_VARIABLES,
        "FD_STEP_TEXTS": ("1e-15", "5e-16"),
        "FD_GATE_TEXT": "1e-24",
        "EXPECTED_F": (2640, 17040, 28800, 14400),
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(ACTION_SOURCE), "exec"), namespace)
    return namespace


def matrix_from_text(rows):
    return mp.matrix([[mp.mpf(value) for value in row] for row in rows])


def matrix_max(matrix):
    return max(abs(value) for value in matrix)


def total_action(evaluate_schedule, combinatorics, geometry, coordinates, masses):
    gravitational = evaluate_schedule(
        combinatorics, geometry, coordinates
    )["action"]
    dust = -8 * mp.pi * mp.fsum(
        masses[rank] * mp.sqrt(coordinates["rho", rank])
        for rank in range(4)
    )
    return gravitational + dust


def centred_second(evaluate_schedule, combinatorics, geometry, base, masses,
                   direction, step, centre):
    plus = dict(base)
    minus = dict(base)
    for index, key in enumerate(VARIABLES):
        plus[key] *= mp.exp(step * direction[index])
        minus[key] *= mp.exp(-step * direction[index])
    upper = total_action(evaluate_schedule, combinatorics, geometry, plus, masses)
    lower = total_action(evaluate_schedule, combinatorics, geometry, minus, masses)
    return (upper - 2 * centre + lower) / (step * step)


def extrapolation(values):
    richardson = tuple(
        (4 * values[index + 1] - values[index]) / 3
        for index in range(4)
    )
    sixth = tuple(
        (16 * richardson[index + 1] - richardson[index]) / 15
        for index in range(3)
    )
    eighth = tuple(
        (64 * sixth[index + 1] - sixth[index]) / 63
        for index in range(2)
    )
    return richardson, sixth, eighth


def coefficient_directions():
    return {
        "first_basis_vector": mp.matrix([1] + [0] * 10),
        "all_ones": mp.matrix([1] * 11),
        "alternating_signs": mp.matrix([
            1 if index % 2 == 0 else -1 for index in range(11)
        ]),
    }


def polynomial_value(x):
    return (
        mp.mpf(7) * x * x / 2
        + 11 * x ** 4
        + 13 * x ** 6
        + 17 * x ** 8
        + 19 * x ** 10
    )


def polynomial_ladder(dps):
    with mp.workdps(dps):
        values = []
        for text_step in STEP_TEXTS:
            step = mp.mpf(text_step)
            values.append(
                (polynomial_value(step) - 2 * polynomial_value(0)
                 + polynomial_value(-step)) / (step * step)
            )
        return tuple(values), extrapolation(values)


print("=" * 78)
print("CONSTRAINED H4 DIRECTIONAL TRUNCATION DIAGNOSTIC")
print("=" * 78)

paths = {
    "action_source": ACTION_SOURCE,
    "curvature": CURVATURE,
    "primary": PRIMARY,
    "primary_result": PRIMARY_RESULT,
    "protocol": PROTOCOL,
    "cell600": CELL600,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "the failed primary artifact and diagnostic protocol have exact provenance",
    actual_hashes == EXPECTED_HASHES and PROTOCOL_COMMIT == "0baeeee",
)

curvature = json.loads(CURVATURE.read_text())
primary = json.loads(PRIMARY.read_text())
upstream_ok = check(
    "the diagnostic starts from exactly the frozen formal control failure",
    curvature["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and primary["outcome"] == "REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED"
    and primary["tests"] == {"passed": 18, "total": 19}
    and primary["census"]["class_count"] == 1
    and len(primary["controls"]["directional_records"]) == 12,
)

actions = load_action_definitions()
definitions_ok = check(
    "only frozen action definitions are loaded and no Hessian verifier is executed",
    {
        "tetrahedra_from_adjacency",
        "barycentric_chambers",
        "schedule_combinatorics",
        "exact_geometry",
        "base_coordinates",
        "evaluate_schedule",
    } <= set(actions)
    and "OUTPUT" not in actions,
)

parser_plus = parse_mpmath_complex("(1.25 + 2.5e-3j)")
parser_minus = parse_mpmath_complex("(-3.5 - 4.25e-7j)")
parser_ok = check(
    "the exact stored-complex parser handles both frozen signs",
    parser_plus == mp.mpc(mp.mpf("1.25"), mp.mpf("2.5e-3"))
    and parser_minus == mp.mpc(mp.mpf("-3.5"), mp.mpf("-4.25e-7")),
)

_, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
_, top, colours = actions["barycentric_chambers"](coarse_top)
orders = tuple(permutations(range(4)))
combinatorics = tuple(
    actions["schedule_combinatorics"](top, colours, order) for order in orders
)
topology_ok = check(
    "all four frozen schedules are reconstructed inside the complete 24-order census",
    len(combinatorics) == 24
    and all(
        combinatorics[index]["pentachora"] == 57600
        and combinatorics[index]["triangles"] == 149280
        for index in DIRECTIONAL_INDICES
    ),
)

# Frozen artifact reconstruction and every post-ladder comparison must retain
# at least the larger action precision; workdps blocks below temporarily round
# the action evaluations to their declared 140/180-digit levels.
mp.mp.dps = 180
direction_labels = coefficient_directions()
primary_direction_map = {
    (tuple(item["order"]), item["direction"]): item
    for item in primary["controls"]["directional_records"]
}

directions = {}
quadratics = {}
reconstruction_errors = []
maximum_coarse_displacement = mp.mpf(0)
maximum_fine_displacement = mp.mpf(0)
for index in DIRECTIONAL_INDICES:
    schedule_record = primary["census"]["schedules"][index]
    boundary_basis = matrix_from_text(primary["bases"]["primary_boundary_basis"])
    lift = matrix_from_text(schedule_record["primary_lift"])
    response = matrix_from_text(schedule_record["primary_response"])
    for label, coefficient in direction_labels.items():
        boundary = boundary_basis * coefficient
        internal = lift * coefficient
        full = mp.matrix(list(boundary) + list(internal))
        computed_q = (coefficient.T * response * coefficient)[0]
        frozen = primary_direction_map[(orders[index], label)]
        frozen_q = mp.mpf(frozen["response_quadratic"])
        reconstruction_errors.append(abs(computed_q - frozen_q))
        directions[index, label] = full
        quadratics[index, label] = frozen_q
        maximum_coarse_displacement = max(
            maximum_coarse_displacement, mp.mpf(STEP_TEXTS[0]) * matrix_max(full)
        )
        maximum_fine_displacement = max(
            maximum_fine_displacement, mp.mpf(STEP_TEXTS[-1]) * matrix_max(full)
        )
reconstruction_ok = check(
    "the twelve stored lifts and quadratic values are reconstructed without a Hessian",
    max(reconstruction_errors) < mp.mpf("1e-55"),
    f"max q error={mp_text(max(reconstruction_errors), 8)}",
)
displacement_ok = check(
    "all frozen coarse and finest coordinate displacements remain in the preregistered regime",
    maximum_coarse_displacement < mp.mpf("2e-5")
    and maximum_fine_displacement < mp.mpf("2e-6"),
    f"coarse={mp_text(maximum_coarse_displacement, 8)}, "
    f"fine={mp_text(maximum_fine_displacement, 8)}",
)

ladders = {}
for dps in DPS_LEVELS:
    with mp.workdps(dps):
        geometry = actions["exact_geometry"](dps)
        geometry["mass"] = mp.mpf(0)
        base = actions["base_coordinates"](geometry)
        masses = tuple(
            mp.mpf(value)
            for value in curvature["selected_rank_matter"]["total_masses"]
        )
        steps = tuple(mp.mpf(value) for value in STEP_TEXTS)
        for index in DIRECTIONAL_INDICES:
            centre = total_action(
                actions["evaluate_schedule"], combinatorics[index], geometry,
                base, masses,
            )
            for label in direction_labels:
                full = mp.matrix(directions[index, label])
                values = tuple(
                    centred_second(
                        actions["evaluate_schedule"], combinatorics[index],
                        geometry, base, masses, full, step, centre,
                    )
                    for step in steps
                )
                richardson, sixth, eighth = extrapolation(values)
                ladders[dps, index, label] = {
                    "centred": values,
                    "richardson": richardson,
                    "sixth": sixth,
                    "eighth": eighth,
                }
    print(f"[INFO] action ladders completed at {dps} digits", flush=True)

diagnostic_records = []
all_asymptotic = []
all_matches = []
all_imaginary = []
all_corruptions = []
original_reproduction_errors = []
for index in DIRECTIONAL_INDICES:
    for label in direction_labels:
        q = quadratics[index, label]
        low = ladders[140, index, label]
        high = ladders[180, index, label]
        final = high["eighth"][1]
        envelope = (
            100 * max(
                abs(high["eighth"][0] - high["eighth"][1]),
                abs(low["eighth"][1] - high["eighth"][1]),
            )
            + mp.mpf("1e-50") * max(mp.mpf(1), abs(final))
        )
        richardson_errors = tuple(abs(value - q) for value in high["richardson"])
        sixth_errors = tuple(abs(value - q) for value in high["sixth"])
        ratios = tuple(
            richardson_errors[j] / richardson_errors[j + 1]
            for j in range(3)
        )
        monotone_r = all(
            richardson_errors[j] > richardson_errors[j + 1]
            for j in range(3)
        )
        ratio_ok = all(mp.mpf(8) <= ratio <= mp.mpf(32) for ratio in ratios)
        monotone_x = all(
            sixth_errors[j] > sixth_errors[j + 1] for j in range(2)
        )
        asymptotic = monotone_r and ratio_ok and monotone_x
        final_error = abs(final - q)
        matched = final_error <= envelope
        imaginary = max(
            *(abs(mp.im(value)) for value in high["centred"]),
            *(abs(mp.im(value)) for value in high["richardson"]),
            *(abs(mp.im(value)) for value in high["sixth"]),
            *(abs(mp.im(value)) for value in high["eighth"]),
        )
        imaginary_ok = imaginary <= envelope
        q_bad = q + mp.mpf("1e-12") * max(mp.mpf(1), abs(q))
        corrupted_error = abs(final - q_bad)
        corruption_ok = corrupted_error > mp.mpf("1e6") * envelope
        frozen_action = parse_mpmath_complex(
            primary_direction_map[(orders[index], label)]["action_richardson"]
        )
        original_reproduction_errors.append(
            abs(low["richardson"][0] - frozen_action)
        )
        all_asymptotic.append(asymptotic)
        all_matches.append(matched)
        all_imaginary.append(imaginary_ok)
        all_corruptions.append(corruption_ok)
        diagnostic_records.append({
            "order": orders[index],
            "direction": label,
            "quadratic": q,
            "centred_180": high["centred"],
            "richardson_180": high["richardson"],
            "sixth_180": high["sixth"],
            "eighth_140": low["eighth"],
            "eighth_180": high["eighth"],
            "richardson_errors": richardson_errors,
            "sixth_errors": sixth_errors,
            "richardson_ratios": ratios,
            "asymptotic": asymptotic,
            "final_error": final_error,
            "envelope": envelope,
            "matched": matched,
            "maximum_imaginary": imaginary,
            "corrupted_error": corrupted_error,
        })

original_ok = check(
    "the first Richardson level reproduces the frozen failed action estimates",
    max(original_reproduction_errors) < mp.mpf("1e-45"),
    f"max error={mp_text(max(original_reproduction_errors), 8)}",
)

poly140 = polynomial_ladder(140)
poly180 = polynomial_ladder(180)
poly_final = poly180[1][2][1]
poly_envelope = (
    100 * max(
        abs(poly180[1][2][0] - poly180[1][2][1]),
        abs(poly140[1][2][1] - poly180[1][2][1]),
    )
    + mp.mpf("1e-50") * max(mp.mpf(1), abs(poly_final))
)
poly_error = abs(poly_final - 7)
polynomial_ok = check(
    "the exact even-polynomial control recovers its second derivative after the same ladder",
    poly_error <= poly_envelope,
    f"error={mp_text(poly_error, 8)}, envelope={mp_text(poly_envelope, 8)}",
)

precision_imaginary_ok = check(
    "all twelve final estimates are precision-resolved and real inside their envelopes",
    all(all_imaginary),
    f"max precision ratio={mp_text(max(
        abs(item['eighth_140'][1] - item['eighth_180'][1]) / item['envelope']
        for item in diagnostic_records
    ), 8)}",
)
corruption_ok = check(
    "every deliberately corrupted quadratic value is rejected",
    all(all_corruptions),
    f"min corruption/envelope={mp_text(min(
        item['corrupted_error'] / item['envelope'] for item in diagnostic_records
    ), 8)}",
)
asymptotic_census_ok = check(
    "the complete truncation census is available for all twelve directions",
    len(diagnostic_records) == 12 and len(all_asymptotic) == 12,
    f"asymptotic={sum(all_asymptotic)}/12, "
    f"ratio range=[{mp_text(min(min(item['richardson_ratios']) for item in diagnostic_records), 6)},"
    f"{mp_text(max(max(item['richardson_ratios']) for item in diagnostic_records), 6)}]",
)
match_census_ok = check(
    "the complete extrapolated Hessian/action comparison is available",
    len(all_matches) == 12,
    f"matched={sum(all_matches)}/12, max error/envelope={mp_text(max(
        item['final_error'] / item['envelope'] for item in diagnostic_records
    ), 8)}",
)

scope = {
    "hessian_or_internal_solve_computed": False,
    "schedule_class_recomputed": False,
    "root_search_or_nested_census_executed": False,
    "nonhomogeneous_operator_or_spectrum_computed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
}
scope_ok = check(
    "the diagnostic remains inside the frozen directional-control scope",
    not any(scope.values()),
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    definitions_ok,
    parser_ok,
    topology_ok,
    reconstruction_ok,
    displacement_ok,
    original_ok,
    polynomial_ok,
    precision_imaginary_ok,
    corruption_ok,
    asymptotic_census_ok,
    match_census_ok,
    scope_ok,
))

if not controls_ok:
    outcome = "REFINED_H4_DIRECTIONAL_DIAGNOSTIC_CONTROL_FAILED"
elif not all(all_asymptotic):
    outcome = "REFINED_H4_DIRECTIONAL_DIAGNOSTIC_NONASYMPTOTIC"
elif not all(all_matches):
    outcome = "REFINED_H4_DIRECTIONAL_HESSIAN_ACTION_MISMATCH"
else:
    outcome = "REFINED_H4_DIRECTIONAL_TRUNCATION_CONFIRMED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one directional-diagnostic outcome",
    outcome in {
        "REFINED_H4_DIRECTIONAL_DIAGNOSTIC_CONTROL_FAILED",
        "REFINED_H4_DIRECTIONAL_DIAGNOSTIC_NONASYMPTOTIC",
        "REFINED_H4_DIRECTIONAL_HESSIAN_ACTION_MISMATCH",
        "REFINED_H4_DIRECTIONAL_TRUNCATION_CONFIRMED",
    },
    outcome,
)

artifact = {
    "title": "Constrained H4 directional truncation diagnostic",
    "date": "2026-08-21",
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "decimal_precisions": list(DPS_LEVELS),
        "steps": list(STEP_TEXTS),
        "schedule_indices": list(DIRECTIONAL_INDICES),
        "direction_labels": list(direction_labels),
        "ladder": (
            "D centered second; R cancels h^2; X cancels h^4; "
            "Y cancels h^6"
        ),
    },
    "displacements": {
        "maximum_coarse": mp_text(maximum_coarse_displacement),
        "maximum_finest": mp_text(maximum_fine_displacement),
    },
    "polynomial_control": {
        "final": mp_text(poly_final),
        "error": mp_text(poly_error),
        "envelope": mp_text(poly_envelope),
    },
    "census": {
        "direction_count": 12,
        "asymptotic_count": sum(all_asymptotic),
        "matched_count": sum(all_matches),
        "records": [
            {
                "order": list(item["order"]),
                "direction": item["direction"],
                "quadratic": mp_text(item["quadratic"]),
                "centred_180": [mp_text(value) for value in item["centred_180"]],
                "richardson_180": [mp_text(value) for value in item["richardson_180"]],
                "sixth_180": [mp_text(value) for value in item["sixth_180"]],
                "eighth_140": [mp_text(value) for value in item["eighth_140"]],
                "eighth_180": [mp_text(value) for value in item["eighth_180"]],
                "richardson_errors": [
                    mp_text(value) for value in item["richardson_errors"]
                ],
                "sixth_errors": [mp_text(value) for value in item["sixth_errors"]],
                "richardson_ratios": [
                    mp_text(value) for value in item["richardson_ratios"]
                ],
                "asymptotic": item["asymptotic"],
                "final_error": mp_text(item["final_error"]),
                "envelope": mp_text(item["envelope"]),
                "matched": item["matched"],
                "maximum_imaginary": mp_text(item["maximum_imaginary"]),
                "corrupted_error": mp_text(item["corrupted_error"]),
            }
            for item in diagnostic_records
        ],
    },
    "scope": scope,
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
print(f"SHA-256: {digest(OUTPUT)}")
print("No Hessian, full suite or deferred nonlinear census was run.")

if passed != tests:
    sys.exit(1)

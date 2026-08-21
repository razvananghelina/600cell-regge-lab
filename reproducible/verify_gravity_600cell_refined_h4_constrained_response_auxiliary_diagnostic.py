#!/usr/bin/env python3
"""Diagnose the two auxiliary failures in the direct constrained H4 audit."""

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
FEASIBILITY = HERE / "gravity_600cell_refined_canonical_map_feasibility.json"
NULL_COUPLING = HERE / "gravity_600cell_refined_h4_null_coupling_adversarial.json"
FAILED_VERIFIER = (
    HERE / "verify_gravity_600cell_refined_h4_constrained_response_adversarial.py"
)
FAILED_ARTIFACT = (
    HERE / "gravity_600cell_refined_h4_constrained_response_adversarial.json"
)
FIRST_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_adversarial_first_result.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic_protocol.md"
)
CELL600 = ROOT / "commons/cell600.py"
OUTPUT = (
    HERE / "gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic.json"
)

PROTOCOL_COMMIT = "9e44dd7"
EXPECTED_HASHES = {
    "action_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "curvature": "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "feasibility": "ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e",
    "null_coupling": "5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9",
    "failed_verifier": "78f6b52f6f019a150a86ddadcb819b67c3757244c015687ab67f4649784ac53d",
    "failed_artifact": "a23ef4cc23d08ad8768f1df66789aa900cdb95a7f3529486df80697a53b1fe81",
    "first_result": "c4203c07b859ed323ee5049875d54d4894a1b815c3792b7c1d1de0e71677ad64",
    "protocol": "2f6d9d72e04c4baf1dc385425ef7f26ba0a55f6249d6505db111aa21e0836405",
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
STATIONARITY_DPS = (180, 220)
STATIONARITY_STEPS = (
    "4e-15", "2e-15", "1e-15", "5e-16", "2.5e-16", "1.25e-16"
)
BRANCH_DPS = 180
BRANCH_STEPS = ("1e-10", "5e-11", "2.5e-11")
BOUNDARY_PIVOT = 3
INTERNAL_PIVOT = 9

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


def matrix_max(matrix):
    return max(abs(value) for value in matrix) if matrix.rows * matrix.cols else mp.mpf(0)


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
        "FD_GATE_TEXT": "1e-60",
        "EXPECTED_F": (2640, 17040, 28800, 14400),
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(ACTION_SOURCE), "exec"), namespace)
    return namespace


def kernel_basis(vector, pivot):
    columns = tuple(index for index in range(len(vector)) if index != pivot)
    if vector[pivot] == 0:
        raise ZeroDivisionError("kernel-basis pivot is zero")
    result = mp.matrix(len(vector), len(vector) - 1)
    for column, index in enumerate(columns):
        result[index, column] = 1
        result[pivot, column] = -vector[index] / vector[pivot]
    return result


def analytic_tangent(base):
    tau = mp.mpf(TAU_TEXT)
    return mp.matrix(
        [-tau * tau / base[("cross",) + pair] for pair in PAIR4]
        + [mp.mpf(1)] * 4
    )


def direct_slice(boundary_basis, internal_basis):
    result = mp.matrix(22, 20)
    for row in range(12):
        for column in range(11):
            result[row, column] = boundary_basis[row, column]
    for row in range(10):
        for column in range(9):
            result[12 + row, 11 + column] = internal_basis[row, column]
    return result


def shifted_coordinates(base, direction, step):
    result = dict(base)
    for index, key in enumerate(VARIABLES):
        result[key] *= mp.exp(step * direction[index])
    return result


def new_tracker():
    return {
        "evaluations": 0,
        "maximum_angle_identity_residual": mp.mpf(0),
        "maximum_imaginary_curvature": mp.mpf(0),
        "maximum_relative_imaginary_action": mp.mpf(0),
        "minimum_angle_argument": mp.inf,
    }


def scalar_record(combinatorics, geometry, coordinates, masses, tracker):
    angle_lookup = {}
    maximum_identity = mp.mpf(0)
    minimum_argument = mp.inf
    for simplex_record in combinatorics["simplex_types"]:
        states = simplex_record["states"]
        squared = actions["simplex_squared"](states, coordinates)
        angles, identity, argument = actions["angle_record"](squared)
        maximum_identity = max(maximum_identity, identity)
        minimum_argument = min(minimum_argument, argument)
        for local_triangle, angle in angles.items():
            triangle_states = actions["canonical_states"](
                states[index] for index in local_triangle
            )
            angle_lookup[states, triangle_states] = angle

    gravitational_sum = mp.mpc(0)
    curvatures = {}
    for triangle_record in combinatorics["triangle_types"]:
        states = triangle_record["states"]
        boundary = len({layer for _, layer in states}) == 1
        curvature_value = mp.pi if boundary else 2 * mp.pi
        for contribution in triangle_record["contributions"]:
            curvature_value += contribution["multiplicity"] * angle_lookup[
                contribution["simplex"], contribution["triangle"]
            ]
        values = []
        for left, right in combinations(range(3), 2):
            key, jacobian = actions["edge_coordinate"](states[left], states[right])
            values.append(jacobian * coordinates[key])
        x, y, z = values
        area_square = (2 * (x * y + x * z + y * z) - x*x - y*y - z*z) / 16
        area = mp.sqrt(mp.mpc(area_square))
        gravitational_sum += triangle_record["count"] * area * curvature_value
        curvatures[states] = curvature_value

    action = -mp.j * gravitational_sum - 8 * mp.pi * mp.fsum(
        masses[rank] * mp.sqrt(coordinates["rho", rank])
        for rank in range(4)
    )
    maximum_imaginary_curvature = max(
        abs(mp.im(value)) for value in curvatures.values()
    )
    tracker["evaluations"] += 1
    tracker["maximum_angle_identity_residual"] = max(
        tracker["maximum_angle_identity_residual"], maximum_identity
    )
    tracker["maximum_imaginary_curvature"] = max(
        tracker["maximum_imaginary_curvature"], maximum_imaginary_curvature
    )
    tracker["maximum_relative_imaginary_action"] = max(
        tracker["maximum_relative_imaginary_action"],
        abs(mp.im(action)) / max(mp.mpf(1), abs(mp.re(action))),
    )
    tracker["minimum_angle_argument"] = min(
        tracker["minimum_angle_argument"], minimum_argument
    )
    return {
        "action": action,
        "curvatures": curvatures,
        "maximum_imaginary_curvature": maximum_imaginary_curvature,
    }


def derivative_extrapolation(values):
    richardson = tuple(
        (4 * values[index + 1] - values[index]) / 3
        for index in range(5)
    )
    sixth = tuple(
        (16 * richardson[index + 1] - richardson[index]) / 15
        for index in range(4)
    )
    eighth = tuple(
        (64 * sixth[index + 1] - sixth[index]) / 63
        for index in range(3)
    )
    tenth = tuple(
        (256 * eighth[index + 1] - eighth[index]) / 255
        for index in range(2)
    )
    return richardson, sixth, eighth, tenth


def scalar_first_ladder(schedule, geometry, base, masses, direction, tracker):
    values = []
    for step_text in STATIONARITY_STEPS:
        step = mp.mpf(step_text)
        plus = scalar_record(
            schedule, geometry, shifted_coordinates(base, direction, step),
            masses, tracker,
        )["action"]
        minus = scalar_record(
            schedule, geometry, shifted_coordinates(base, direction, -step),
            masses, tracker,
        )["action"]
        values.append((plus - minus) / (2 * step))
    richardson, sixth, eighth, tenth = derivative_extrapolation(tuple(values))
    envelope = (
        100 * abs(tenth[1] - tenth[0])
        + mp.mpf("1e-60") * max(mp.mpf(1), *(abs(value) for value in values))
    )
    legacy = (4 * values[3] - values[2]) / 3
    return {
        "centred": tuple(values),
        "richardson": richardson,
        "sixth": sixth,
        "eighth": eighth,
        "tenth": tenth,
        "final": tenth[1],
        "envelope": envelope,
        "legacy": legacy,
    }


def polynomial_first_ladder(linear):
    values = []
    for step_text in STATIONARITY_STEPS:
        step = mp.mpf(step_text)

        def value(x):
            return linear * x + 7*x*x + 11*x**4 + 13*x**6 + 17*x**8 + 19*x**10

        values.append((value(step) - value(-step)) / (2 * step))
    *_, tenth = derivative_extrapolation(tuple(values))
    envelope = (
        100 * abs(tenth[1] - tenth[0])
        + mp.mpf("1e-60") * max(mp.mpf(1), *(abs(value) for value in values))
    )
    return tenth[1], envelope


print("=" * 78)
print("CONSTRAINED H4 AUXILIARY FAILURE DIAGNOSTIC")
print("=" * 78)

paths = {
    "action_source": ACTION_SOURCE,
    "curvature": CURVATURE,
    "feasibility": FEASIBILITY,
    "null_coupling": NULL_COUPLING,
    "failed_verifier": FAILED_VERIFIER,
    "failed_artifact": FAILED_ARTIFACT,
    "first_result": FIRST_RESULT,
    "protocol": PROTOCOL,
    "cell600": CELL600,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "the failed audit, scalar action and diagnostic protocol have exact provenance",
    actual_hashes == EXPECTED_HASHES and PROTOCOL_COMMIT == "9e44dd7",
)

curvature = json.loads(CURVATURE.read_text())
null_coupling = json.loads(NULL_COUPLING.read_text())
failed = json.loads(FAILED_ARTIFACT.read_text())
upstream_ok = check(
    "the diagnostic starts from exactly the frozen two-control failure",
    curvature["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and null_coupling["outcome"]
        == "ADVERSARIAL_REFINED_H4_NULL_COUPLING_CORROBORATED"
    and failed["outcome"]
        == "ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED"
    and failed["tests"] == {"passed": 15, "total": 17}
    and failed["census"]["class_count"] == 1
    and failed["census"]["primary_match_count"] == 24,
)

actions = load_action_definitions()
definitions_ok = check(
    "only frozen geometric definitions are loaded and no response matrix is imported",
    {
        "tetrahedra_from_adjacency", "barycentric_chambers", "all_simplices",
        "schedule_combinatorics", "exact_geometry", "base_coordinates",
        "angle_record", "simplex_squared", "edge_coordinate",
    } <= set(actions)
    and "OUTPUT" not in actions,
)

_, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
_, top, colours = actions["barycentric_chambers"](coarse_top)
spatial_cells = actions["all_simplices"](tuple(map(tuple, top)))
orders = tuple(permutations(range(4)))
combinatorics = tuple(
    actions["schedule_combinatorics"](top, colours, order) for order in orders
)
topology_ok = check(
    "the same exact carrier and all 24 schedules are reconstructed",
    tuple(len(layer) for layer in spatial_cells) == (2640, 17040, 28800, 14400)
    and len(combinatorics) == 24
    and all(
        record["pentachora"] == 57600
        and record["triangles"] == 149280
        and len(record["triangle_types"]) == 28
        for record in combinatorics
    ),
)

mp.mp.dps = 230
geometries = {dps: actions["exact_geometry"](dps) for dps in STATIONARITY_DPS}
for geometry in geometries.values():
    geometry["mass"] = mp.mpf(0)
bases = {dps: actions["base_coordinates"](geometries[dps]) for dps in STATIONARITY_DPS}
masses = tuple(
    mp.mpf(value) for value in curvature["selected_rank_matter"]["total_masses"]
)
c = mp.matrix([
    mp.mpf(value) for value in null_coupling["compatibility"]["adversarial_row"]
])
n = analytic_tangent(bases[220])
p = kernel_basis(c, BOUNDARY_PIVOT)
q = kernel_basis(n, INTERNAL_PIVOT)
slice_matrix = direct_slice(p, q)
basis_error = max(matrix_max(c.T * p), matrix_max(n.T * q))
bases_ok = check(
    "the diagnostic slice satisfies the frozen compatibility identities",
    basis_error < mp.mpf("1e-180"),
    f"max identity error={mp_text(basis_error, 8)}",
)

with mp.workdps(220):
    even_final, even_envelope = polynomial_first_ladder(mp.mpf(0))
    nonzero_final, nonzero_envelope = polynomial_first_ladder(mp.mpf("1e-20"))
stationarity_controls_ok = check(
    "the derivative ladder accepts a known stationary action and rejects a nonzero gradient",
    abs(even_final) <= even_envelope
    and abs(nonzero_final) > nonzero_envelope
    and abs(nonzero_final - mp.mpf("1e-20")) < mp.mpf("1e-100"),
    f"zero={mp_text(even_final, 8)}, nonzero={mp_text(nonzero_final, 8)}",
)

stationarity = {}
stationarity_trackers = {}
for dps in STATIONARITY_DPS:
    tracker = new_tracker()
    stationarity_trackers[dps] = tracker
    with mp.workdps(dps):
        for schedule_index, schedule in enumerate(combinatorics):
            for internal_index in range(10):
                direction = mp.matrix(22, 1)
                direction[12 + internal_index] = 1
                stationarity[dps, schedule_index, internal_index] = scalar_first_ladder(
                    schedule, geometries[dps], bases[dps], masses, direction, tracker
                )
            if schedule_index in (5, 11, 17, 23):
                print(
                    f"[INFO] stationarity {dps} digits: {schedule_index + 1}/24",
                    flush=True,
                )

stationarity_passes = []
precision_passes = []
stationarity_fractions = []
precision_fractions = []
for schedule_index in range(24):
    for internal_index in range(10):
        left = stationarity[180, schedule_index, internal_index]
        right = stationarity[220, schedule_index, internal_index]
        for record in (left, right):
            stationarity_passes.append(abs(record["final"]) <= record["envelope"])
            stationarity_fractions.append(abs(record["final"]) / record["envelope"])
        difference = abs(left["final"] - right["final"])
        gate = left["envelope"] + right["envelope"]
        precision_passes.append(difference <= gate)
        precision_fractions.append(difference / gate)
maximum_legacy = max(
    abs(record["legacy"]) for record in stationarity.values()
)
maximum_final = max(abs(record["final"]) for record in stationarity.values())
stationarity_ok = check(
    "all 240 scalar-action internal derivatives converge to the frozen stationary scale",
    all(stationarity_passes),
    f"legacy max={mp_text(maximum_legacy, 8)}, final max={mp_text(maximum_final, 8)}, "
    f"max fraction={mp_text(max(stationarity_fractions), 8)}",
)
stationarity_precision_ok = check(
    "all stationarity ladders agree between 180 and 220 decimal digits",
    all(precision_passes),
    f"max fraction={mp_text(max(precision_fractions), 8)}",
)

branch_tracker = new_tracker()
branch_step_maxima = [mp.mpf(0)] * len(BRANCH_STEPS)
base_maximum = mp.mpf(0)
finest_parity_ratios = []
resolved_parity_count = 0
with mp.workdps(BRANCH_DPS):
    geometry = geometries[180]
    base = bases[180]
    directions = [slice_matrix[:, index] for index in range(20)]
    directions.extend(
        slice_matrix[:, left] + slice_matrix[:, right]
        for left, right in combinations(range(20), 2)
    )
    for schedule_index, schedule in enumerate(combinatorics):
        base_record = scalar_record(schedule, geometry, base, masses, branch_tracker)
        base_maximum = max(base_maximum, base_record["maximum_imaginary_curvature"])
        for direction in directions:
            for step_index, step_text in enumerate(BRANCH_STEPS):
                step = mp.mpf(step_text)
                plus = scalar_record(
                    schedule, geometry,
                    shifted_coordinates(base, direction, step),
                    masses, branch_tracker,
                )
                minus = scalar_record(
                    schedule, geometry,
                    shifted_coordinates(base, direction, -step),
                    masses, branch_tracker,
                )
                branch_step_maxima[step_index] = max(
                    branch_step_maxima[step_index],
                    plus["maximum_imaginary_curvature"],
                    minus["maximum_imaginary_curvature"],
                )
                if step_index == len(BRANCH_STEPS) - 1:
                    keys = set(base_record["curvatures"])
                    if keys != set(plus["curvatures"]) or keys != set(minus["curvatures"]):
                        raise RuntimeError("triangle curvature keys changed under perturbation")
                    odd_scale = max(
                        abs(
                            mp.im(plus["curvatures"][key])
                            - mp.im(minus["curvatures"][key])
                        )
                        for key in keys
                    )
                    even_error = max(
                        abs(
                            mp.im(plus["curvatures"][key])
                            + mp.im(minus["curvatures"][key])
                            - 2 * mp.im(base_record["curvatures"][key])
                        )
                        for key in keys
                    )
                    if odd_scale > mp.mpf("1e-100"):
                        resolved_parity_count += 1
                        finest_parity_ratios.append(even_error / odd_scale)
        if schedule_index in (5, 11, 17, 23):
            print(f"[INFO] curvature ladders: {schedule_index + 1}/24", flush=True)

branch_ratios = (
    branch_step_maxima[0] / branch_step_maxima[1],
    branch_step_maxima[1] / branch_step_maxima[2],
)

with mp.workdps(180):
    smooth_values = [abs(7*mp.mpf(step) + 11*mp.mpf(step)**3) for step in BRANCH_STEPS]
    smooth_ratios = (
        smooth_values[0] / smooth_values[1],
        smooth_values[1] / smooth_values[2],
    )
    smooth_parity = mp.mpf(0)
    discontinuous_values = [mp.mpf(1)] * 3
    discontinuous_ratios = (
        discontinuous_values[0] / discontinuous_values[1],
        discontinuous_values[1] / discontinuous_values[2],
    )
branch_controls_ok = check(
    "the smooth complex control passes and a discontinuous branch fails the halving gate",
    all(mp.mpf("1.99") <= ratio <= mp.mpf("2.01") for ratio in smooth_ratios)
    and smooth_parity < mp.mpf("1e-6")
    and not all(
        mp.mpf("1.99") <= ratio <= mp.mpf("2.01")
        for ratio in discontinuous_ratios
    ),
    f"smooth={tuple(mp_text(x, 8) for x in smooth_ratios)}, "
    f"bad={tuple(mp_text(x, 8) for x in discontinuous_ratios)}",
)

curvature_scaling_ok = check(
    "individual off-shell imaginary curvatures vanish linearly at the real background",
    base_maximum < mp.mpf("1e-100")
    and branch_step_maxima[2] > mp.mpf("1e-12")
    and all(mp.mpf("1.99") <= ratio <= mp.mpf("2.01") for ratio in branch_ratios),
    f"base={mp_text(base_maximum, 8)}, maxima={tuple(mp_text(x, 8) for x in branch_step_maxima)}, "
    f"ratios={tuple(mp_text(x, 8) for x in branch_ratios)}",
)
curvature_parity_ok = check(
    "every resolved finest-step curvature pair has the frozen smooth odd-leading behaviour",
    resolved_parity_count == 24 * 210
    and max(finest_parity_ratios) < mp.mpf("1e-6"),
    f"resolved={resolved_parity_count}/5040, max ratio={mp_text(max(finest_parity_ratios), 8)}",
)
complete_branch_ok = check(
    "the complete action stays real on a safe analytic branch throughout the diagnostic",
    branch_tracker["maximum_relative_imaginary_action"] < mp.mpf("1e-100")
    and branch_tracker["maximum_angle_identity_residual"] < mp.mpf("1e-100")
    and branch_tracker["minimum_angle_argument"] > mp.mpf("1e-2"),
    f"action imag={mp_text(branch_tracker['maximum_relative_imaginary_action'], 8)}, "
    f"identity={mp_text(branch_tracker['maximum_angle_identity_residual'], 8)}, "
    f"min arg={mp_text(branch_tracker['minimum_angle_argument'], 8)}",
)

known_controls_ok = all((
    provenance_ok, upstream_ok, definitions_ok, topology_ok, bases_ok,
    stationarity_controls_ok, branch_controls_ok,
))
physical_diagnostics_ok = all((
    stationarity_ok, stationarity_precision_ok, curvature_scaling_ok,
    curvature_parity_ok, complete_branch_ok,
))
if not known_controls_ok:
    outcome = "REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_DIAGNOSTIC_CONTROL_FAILED"
elif physical_diagnostics_ok:
    outcome = "REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_FAILURES_RESOLVED"
else:
    outcome = "REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_FAILURES_CONFIRMED"

outcome_ok = check(
    "the frozen auxiliary outcome follows mechanically from controls and diagnostics",
    (outcome.endswith("_RESOLVED")) == (known_controls_ok and physical_diagnostics_ok)
    and (outcome.endswith("_CONFIRMED")) == (known_controls_ok and not physical_diagnostics_ok)
    and (outcome.endswith("_CONTROL_FAILED")) == (not known_controls_ok),
    outcome,
)

artifact = {
    "title": "Constrained H4 response auxiliary failure diagnostic",
    "date": "2026-08-21",
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "stationarity": {
        "decimal_precisions": STATIONARITY_DPS,
        "steps": STATIONARITY_STEPS,
        "component_count_per_precision": 240,
        "maximum_legacy_fourth_order_residual": mp_text(maximum_legacy),
        "maximum_tenth_order_residual": mp_text(maximum_final),
        "maximum_zero_envelope_fraction": mp_text(max(stationarity_fractions)),
        "maximum_precision_fraction": mp_text(max(precision_fractions)),
        "all_zero_gates_pass": all(stationarity_passes),
        "all_precision_gates_pass": all(precision_passes),
    },
    "off_shell_curvature": {
        "decimal_precision": BRANCH_DPS,
        "steps": BRANCH_STEPS,
        "directions_per_schedule": 210,
        "base_maximum_imaginary_curvature": mp_text(base_maximum),
        "step_maximum_imaginary_curvatures": [
            mp_text(value) for value in branch_step_maxima
        ],
        "halving_ratios": [mp_text(value) for value in branch_ratios],
        "resolved_parity_count": resolved_parity_count,
        "maximum_even_to_odd_ratio": mp_text(max(finest_parity_ratios)),
        "maximum_relative_imaginary_complete_action": mp_text(
            branch_tracker["maximum_relative_imaginary_action"]
        ),
        "maximum_angle_identity_residual": mp_text(
            branch_tracker["maximum_angle_identity_residual"]
        ),
        "minimum_angle_argument": mp_text(branch_tracker["minimum_angle_argument"]),
        "evaluation_count": branch_tracker["evaluations"],
    },
    "controls": {
        "stationary_polynomial_final": mp_text(even_final),
        "nonstationary_polynomial_final": mp_text(nonzero_final),
        "nonstationary_polynomial_envelope": mp_text(nonzero_envelope),
        "smooth_branch_halving_ratios": [mp_text(value) for value in smooth_ratios],
        "discontinuous_branch_halving_ratios": [
            mp_text(value) for value in discontinuous_ratios
        ],
    },
    "interpretation": {
        "legacy_stationarity_failure": (
            "fourth-order finite-difference truncation"
            if stationarity_ok and stationarity_precision_ok else "unresolved"
        ),
        "individual_curvature_reality_failure": (
            "overstrong off-shell termwise reality gate"
            if curvature_scaling_ok and curvature_parity_ok and complete_branch_ok
            else "unresolved"
        ),
        "failed_audit_outcome_retroactively_changed": False,
        "corrected_adversarial_run_licensed": outcome.endswith("_RESOLVED"),
    },
    "scope": {
        "response_matrix_or_class_recomputed": False,
        "primary_response_loaded": False,
        "full_suite_executed": False,
        "root_search_or_deferred_census_executed": False,
        "physical_constant_extracted": False,
    },
    "status_labels": {
        "auxiliary_failure_diagnosis": (
            "DERIVED COMPUTATIONAL" if outcome.endswith("_RESOLVED") else "OPEN"
        ),
        "primary_single_class_acceptance": "OPEN",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")

if passed != tests:
    raise SystemExit(1)


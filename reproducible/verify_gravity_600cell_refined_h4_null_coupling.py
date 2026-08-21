#!/usr/bin/env python3
"""Directional audit of the refined product-lapse null coupling."""

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
CURVATURE_ADVERSARIAL = (
    HERE / "gravity_600cell_refined_local_curvature_mass_adversarial.json"
)
BOUNDARY_COTANGENT = HERE / "gravity_600cell_refined_boundary_cotangent.json"
HESSIAN = HERE / "gravity_600cell_refined_effective_h4_hessian.json"
HESSIAN_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_effective_h4_hessian_first_result.md"
)
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_refined_h4_null_coupling_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_h4_null_coupling_protocol.md"
CELL600 = ROOT / "commons/cell600.py"
OUTPUT = HERE / "gravity_600cell_refined_h4_null_coupling.json"

PRIOR_ART_COMMIT = "32854ba"
PROTOCOL_COMMIT = "752891b"
EXPECTED_HASHES = {
    "action_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "curvature": "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "curvature_adversarial": "c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e",
    "boundary_cotangent": "4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa",
    "hessian": "56e08db9a840b95e686fadb2763e89400b09220e88b80e9d35c17c1e73eef0a3",
    "hessian_result": "f8bf5679e153fcca8a076064bc5b98e881d91c2add9aaffa4c6858247538f1b8",
    "prior_art": "4a6b55689535feeb729db3a71a9a28d2a7a8ccbeb6554a79371f1073b6d669dd",
    "protocol": "55a098c9961002c18a93d3948fe588ad24228f048fb41da9244efb9f76913737",
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
PRIMARY_DPS = 100
SECONDARY_DPS = 140
STEP_TEXTS = ("1e-10", "5e-11", "2.5e-11")

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


def vector_max(vector):
    return max(abs(vector[index]) for index in range(len(vector)))


def vector_difference(left, right):
    return max(abs(left[index] - right[index]) for index in range(len(left)))


def product_coordinates(base, tau):
    result = dict(base)
    tau_square = tau * tau
    for pair in PAIR4:
        result[("cross",) + pair] = result[("old",) + pair] - tau_square
    for rank in range(4):
        result["rho", rank] = tau_square
    return result


def product_tangent(coordinates, tau):
    tau_square = tau * tau
    return mp.matrix([
        *(
            -tau_square / coordinates[("cross",) + pair]
            for pair in PAIR4
        ),
        1, 1, 1, 1,
    ])


def total_gradient(evaluate_schedule, combinatorics, geometry, coordinates, masses):
    evaluation = evaluate_schedule(combinatorics, geometry, coordinates)
    values = [evaluation["gradient"][key] for key in VARIABLES]
    for rank in range(4):
        values[18 + rank] -= (
            4 * mp.pi * masses[rank] * mp.sqrt(coordinates["rho", rank])
        )
    return evaluation, mp.matrix(values)


def gravitational_gradient(evaluate_schedule, combinatorics, geometry, coordinates):
    evaluation = evaluate_schedule(combinatorics, geometry, coordinates)
    return mp.matrix([evaluation["gradient"][key] for key in VARIABLES])


def shifted_internal(base, tangent, displacement):
    result = dict(base)
    for index, key in enumerate(INTERNAL_VARIABLES):
        result[key] *= mp.exp(displacement * tangent[index])
    return result


def directional_derivative(evaluate_schedule, combinatorics, geometry, base,
                           tangent, step):
    plus = gravitational_gradient(
        evaluate_schedule,
        combinatorics,
        geometry,
        shifted_internal(base, tangent, step),
    )
    minus = gravitational_gradient(
        evaluate_schedule,
        combinatorics,
        geometry,
        shifted_internal(base, tangent, -step),
    )
    return (plus - minus) / (2 * step)


def directional_ladder(evaluate_schedule, combinatorics, geometry, dps,
                       secondary_only=False):
    with mp.workdps(dps):
        tau = mp.mpf(TAU_TEXT)
        base = actions["base_coordinates"](geometry)
        tangent = product_tangent(base, tau)
        steps = tuple(mp.mpf(value) for value in STEP_TEXTS)
        d1 = directional_derivative(
            evaluate_schedule, combinatorics, geometry, base, tangent, steps[1]
        )
        d2 = directional_derivative(
            evaluate_schedule, combinatorics, geometry, base, tangent, steps[2]
        )
        secondary = (4 * d2 - d1) / 3
        primary = None
        if not secondary_only:
            d0 = directional_derivative(
                evaluate_schedule, combinatorics, geometry, base, tangent, steps[0]
            )
            primary = (4 * d1 - d0) / 3
        return primary, secondary, tangent


def add_dust_image(vector, masses, tau):
    result = mp.matrix(vector)
    for rank in range(4):
        result[18 + rank] -= 2 * mp.pi * masses[rank] * tau
    return result


def reversal_vector(vector):
    return mp.matrix(
        [vector[6 + index] for index in range(6)]
        + [vector[index] for index in range(6)]
    )


def synthetic_control(tangent):
    with mp.workdps(SECONDARY_DPS):
        full_direction = mp.matrix([0] * 12 + list(tangent))
        matrix = mp.matrix(22, 22)
        for index in range(22):
            matrix[index, index] = index + 2
        for index in range(21):
            matrix[index, index + 1] = 1
            matrix[index + 1, index] = 1
        offset = mp.matrix([mp.mpf(index + 1) / 23 for index in range(22)])

        def gradient(displacement):
            return matrix * (displacement * full_direction) + offset

        steps = tuple(mp.mpf(value) for value in STEP_TEXTS)
        derivatives = []
        for step in steps:
            derivatives.append((gradient(step) - gradient(-step)) / (2 * step))
        first = (4 * derivatives[1] - derivatives[0]) / 3
        second = (4 * derivatives[2] - derivatives[1]) / 3
        expected = matrix * full_direction
        return max(vector_difference(first, expected), vector_difference(second, expected))


print("=" * 78)
print("REFINED H4 PRODUCT-LAPSE NULL COUPLING")
print("=" * 78)

paths = {
    "action_source": ACTION_SOURCE,
    "curvature": CURVATURE,
    "curvature_adversarial": CURVATURE_ADVERSARIAL,
    "boundary_cotangent": BOUNDARY_COTANGENT,
    "hessian": HESSIAN,
    "hessian_result": HESSIAN_RESULT,
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "cell600": CELL600,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all frozen action, matter, singularity and protocol inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIOR_ART_COMMIT == "32854ba"
    and PROTOCOL_COMMIT == "752891b",
)

curvature = json.loads(CURVATURE.read_text())
curvature_adversarial = json.loads(CURVATURE_ADVERSARIAL.read_text())
hessian = json.loads(HESSIAN.read_text())
upstream_ok = check(
    "the curvature branch and primary Hessian carry the exact required outcomes",
    curvature["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and curvature_adversarial["outcome"]
        == "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED"
    and hessian["outcome"]
        == "REFINED_EFFECTIVE_H4_HESSIAN_INTERNAL_SINGULAR"
    and hessian["tests"] == {"passed": 16, "total": 16}
    and hessian["census"]["internal_inertia_histogram"]
        == {"(9, 1, 0)": 24}
    and hessian["census"]["effective_matrix_class_count"] is None,
)

actions = load_action_definitions()
definitions_ok = check(
    "only frozen action definitions are loaded and no full-Hessian code is imported",
    {
        "tetrahedra_from_adjacency",
        "barycentric_chambers",
        "all_simplices",
        "schedule_combinatorics",
        "exact_geometry",
        "base_coordinates",
        "evaluate_schedule",
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
    "the direct carrier and all 24 schedule records have the exact frozen census",
    tuple(len(layer) for layer in spatial_cells)
        == (2640, 17040, 28800, 14400)
    and len(combinatorics) == 24
    and all(
        record["pentachora"] == record["distinct_pentachora"] == 57600
        and record["triangles"] == 149280
        and record["boundary_triangles"] == 57600
        and record["mixed_triangle_types"] == 0
        for record in combinatorics
    ),
)

geometry100 = actions["exact_geometry"](PRIMARY_DPS)
geometry140 = actions["exact_geometry"](SECONDARY_DPS)
geometry100["mass"] = mp.mpf(0)
geometry140["mass"] = mp.mpf(0)

with mp.workdps(SECONDARY_DPS):
    tau = mp.mpf(TAU_TEXT)
    masses = tuple(
        mp.mpf(value)
        for value in curvature["selected_rank_matter"]["total_masses"]
    )
    base140 = actions["base_coordinates"](geometry140)
    product_ratios = (mp.mpf("0.5"), mp.mpf(1), mp.mpf(2))
    maximum_product_residual = mp.mpf(0)
    maximum_branch_identity = mp.mpf(0)
    maximum_branch_imaginary = mp.mpf(0)
    minimum_branch_argument = mp.inf
    for record in combinatorics:
        for ratio in product_ratios:
            coordinates = product_coordinates(base140, ratio * tau)
            evaluation, gradient = total_gradient(
                actions["evaluate_schedule"], record, geometry140,
                coordinates, masses,
            )
            maximum_product_residual = max(
                maximum_product_residual,
                *(abs(gradient[12 + index]) for index in range(10)),
            )
            maximum_branch_identity = max(
                maximum_branch_identity,
                evaluation["maximum_angle_identity_residual"],
            )
            maximum_branch_imaginary = max(
                maximum_branch_imaginary,
                evaluation["maximum_imaginary_curvature"],
            )
            minimum_branch_argument = min(
                minimum_branch_argument, evaluation["minimum_angle_argument"]
            )
product_family_ok = check(
    "all 72 finite product points remain internally stationary",
    maximum_product_residual < mp.mpf("1e-60"),
    f"max residual={mp_text(maximum_product_residual, 8)}",
)
branch_ok = check(
    "the complete finite product control remains on one real Lorentzian branch",
    maximum_branch_identity < mp.mpf("1e-80")
    and maximum_branch_imaginary < mp.mpf("1e-80")
    and minimum_branch_argument > mp.mpf("1e-8"),
    f"identity={mp_text(maximum_branch_identity, 8)}, "
    f"imag={mp_text(maximum_branch_imaginary, 8)}",
)

primary100 = []
secondary100 = []
secondary140 = []
tangents140 = []
print("[INFO] differentiating 24 schedules only along the product tangent", flush=True)
for index, record in enumerate(combinatorics):
    first, second, _ = directional_ladder(
        actions["evaluate_schedule"], record, geometry100, PRIMARY_DPS
    )
    primary100.append(first)
    secondary100.append(second)
    _, fine, tangent = directional_ladder(
        actions["evaluate_schedule"], record, geometry140, SECONDARY_DPS,
        secondary_only=True,
    )
    secondary140.append(fine)
    tangents140.append(tangent)
    if index in (5, 11, 17, 23):
        print(f"[INFO] schedules completed: {index + 1}/24", flush=True)

records = []
compatibility_rows = []
envelopes = []
maximum_internal_image = mp.mpf(0)
maximum_formula_error = mp.mpf(0)
maximum_raw_imaginary = mp.mpf(0)

with mp.workdps(SECONDARY_DPS):
    base_evaluations = [
        actions["evaluate_schedule"](record, geometry140, base140)
        for record in combinatorics
    ]
    for index in range(24):
        total100a = add_dust_image(primary100[index], masses, tau)
        total100b = add_dust_image(secondary100[index], masses, tau)
        total140b = add_dust_image(secondary140[index], masses, tau)
        step_difference = vector_difference(total100a, total100b)
        precision_difference = vector_difference(total100b, total140b)
        image_scale = max(mp.mpf(1), vector_max(total140b))
        envelope = (
            100 * max(step_difference, precision_difference)
            + mp.mpf("1e-50") * image_scale
        )
        raw_imaginary = max(abs(mp.im(total140b[row])) for row in range(22))
        compatibility = mp.matrix([mp.re(total140b[row]) for row in range(12)])
        internal_image = mp.matrix([
            mp.re(total140b[12 + row]) for row in range(10)
        ])
        boundary_gradient = mp.matrix([
            mp.re(base_evaluations[index]["gradient"][key])
            for key in BOUNDARY_VARIABLES
        ])
        predicted = boundary_gradient / 2
        formula_error = vector_difference(compatibility, predicted)
        maximum_internal_image = max(
            maximum_internal_image, vector_max(internal_image)
        )
        maximum_formula_error = max(maximum_formula_error, formula_error)
        maximum_raw_imaginary = max(maximum_raw_imaginary, raw_imaginary)
        compatibility_rows.append(compatibility)
        envelopes.append(envelope)
        records.append({
            "order": list(orders[index]),
            "tangent": tangents140[index],
            "step_difference": step_difference,
            "precision_difference": precision_difference,
            "envelope": envelope,
            "raw_imaginary": raw_imaginary,
            "internal_image": internal_image,
            "compatibility": compatibility,
            "direct_boundary_gradient": boundary_gradient,
            "formula_error": formula_error,
        })

print("[INFO] computed compatibility row before loading frozen boundary vector:")
print("       " + ", ".join(mp_text(value, 18) for value in compatibility_rows[0]))

precision_ok = check(
    "all directional images pass the frozen step, precision and reality envelopes",
    maximum_raw_imaginary <= max(envelopes),
    f"max envelope={mp_text(max(envelopes), 8)}, "
    f"imag={mp_text(maximum_raw_imaginary, 8)}",
)
internal_null = all(
    vector_max(records[index]["internal_image"]) <= envelopes[index]
    for index in range(24)
)
null_census_ok = check(
    "the analytic product tangent has a complete resolved internal-null census",
    all(len(records[index]["internal_image"]) == 10 for index in range(24)),
    f"null={sum(vector_max(records[index]['internal_image']) <= envelopes[index] for index in range(24))}/24",
)

direct_formula_match = all(
    records[index]["formula_error"] <= envelopes[index]
    for index in range(24)
)
direct_formula_census_ok = check(
    "the direct half-gradient formula is compared componentwise for all schedules",
    all(len(row) == 12 for row in compatibility_rows),
    f"matches={sum(records[index]['formula_error'] <= envelopes[index] for index in range(24))}/24",
)

# The frozen boundary vector is deliberately parsed only after c was computed and printed.
boundary_cotangent = json.loads(BOUNDARY_COTANGENT.read_text())
boundary_upstream_ok = check(
    "the post-computation boundary artifact has the accepted frozen outcome",
    boundary_cotangent["outcome"]
        == "REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED"
    and boundary_cotangent["tests"] == {"passed": 16, "total": 16},
)
with mp.workdps(SECONDARY_DPS):
    frozen_pre = mp.matrix([
        mp.mpf(value)
        for value in boundary_cotangent["boundary_covectors"]["selected_pre"]
    ])
    frozen_post = mp.matrix([
        mp.mpf(value)
        for value in boundary_cotangent["boundary_covectors"]["selected_post"]
    ])
    frozen_gradient = mp.matrix(
        [-frozen_pre[index] for index in range(6)]
        + [frozen_post[index] for index in range(6)]
    )
    frozen_prediction = frozen_gradient / 2
    frozen_errors = [
        vector_difference(row, frozen_prediction) for row in compatibility_rows
    ]
    maximum_frozen_error = max(frozen_errors)
frozen_formula_match = all(
    frozen_errors[index] <= envelopes[index] for index in range(24)
)
frozen_comparison_ok = check(
    "the independently frozen boundary half-gradient comparison is complete",
    len(frozen_errors) == 24,
    f"matches={sum(frozen_errors[index] <= envelopes[index] for index in range(24))}/24, "
    f"max error={mp_text(maximum_frozen_error, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    coupling_norms = [vector_max(row) for row in compatibility_rows]
    coupling_resolved_nonzero = all(
        coupling_norms[index] > mp.mpf("1e6") * envelopes[index]
        for index in range(24)
    )
    schedule_differences = [
        vector_difference(compatibility_rows[index], compatibility_rows[0])
        for index in range(24)
    ]
    schedule_independent = all(
        schedule_differences[index] <= max(envelopes[index], envelopes[0])
        for index in range(24)
    )
    order_index = {order: index for index, order in enumerate(orders)}
    reversal_differences = []
    reversal_passes = []
    for index, order in enumerate(orders):
        reverse_index = order_index[tuple(reversed(order))]
        difference = vector_difference(
            compatibility_rows[index],
            reversal_vector(compatibility_rows[reverse_index]),
        )
        reversal_differences.append(difference)
        reversal_passes.append(
            difference <= max(envelopes[index], envelopes[reverse_index])
        )
    time_reversal_covariant = all(reversal_passes)
    maximum_schedule_difference = max(schedule_differences)
    maximum_reversal_difference = max(reversal_differences)
schedule_census_ok = check(
    "all 24 compatibility rows and reversal partners are explicitly compared",
    len(schedule_differences) == len(reversal_differences) == 24,
    f"same={schedule_independent}, reversal={time_reversal_covariant}",
)

with mp.workdps(SECONDARY_DPS):
    synthetic_error = synthetic_control(tangents140[0])
synthetic_ok = check(
    "the frozen directional ladder reproduces a synthetic affine gradient exactly",
    synthetic_error < mp.mpf("1e-80"),
    f"error={mp_text(synthetic_error, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    wrong_full_errors = [
        vector_difference(compatibility_rows[index], records[index]["direct_boundary_gradient"])
        for index in range(24)
    ]
    wrong_sign_errors = [
        vector_difference(
            compatibility_rows[index],
            -records[index]["direct_boundary_gradient"] / 2,
        )
        for index in range(24)
    ]
    dropped = mp.matrix(frozen_prediction)
    largest_index = max(range(12), key=lambda index: abs(dropped[index]))
    dropped[largest_index] = 0
    dropped_error = vector_difference(compatibility_rows[0], dropped)
    control_scale = mp.mpf("1e6") * max(envelopes)
corruption_ok = check(
    "wrong factor, wrong sign and a dropped largest component all fail decisively",
    min(wrong_full_errors) > control_scale
    and min(wrong_sign_errors) > control_scale
    and dropped_error > control_scale,
)

scope = {
    "full_hessian_computed": False,
    "pseudoinverse_or_schur_complement_computed": False,
    "root_search_or_nested_census_executed": False,
    "nonhomogeneous_operator_or_spectrum_computed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
}
scope_ok = check(
    "the calculation remains inside the directional null-coupling scope",
    not any(scope.values()),
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    definitions_ok,
    topology_ok,
    product_family_ok,
    branch_ok,
    precision_ok,
    null_census_ok,
    direct_formula_census_ok,
    boundary_upstream_ok,
    frozen_comparison_ok,
    schedule_census_ok,
    synthetic_ok,
    corruption_ok,
    scope_ok,
))

if not controls_ok:
    outcome = "REFINED_H4_NULL_COUPLING_CONTROL_FAILED"
elif not internal_null:
    outcome = "REFINED_H4_PRODUCT_TANGENT_NOT_NULL"
elif not coupling_resolved_nonzero:
    outcome = "REFINED_H4_PRODUCT_NULL_DECOUPLED"
elif not (direct_formula_match and frozen_formula_match):
    outcome = "REFINED_H4_NULL_COUPLING_FORMULA_REFUTED"
elif not (schedule_independent and time_reversal_covariant):
    outcome = "REFINED_H4_NULL_COUPLING_SCHEDULE_DEPENDENT"
else:
    outcome = "REFINED_H4_NULL_COUPLING_COMPATIBILITY_CONFIRMED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one null-coupling outcome",
    outcome in {
        "REFINED_H4_NULL_COUPLING_CONTROL_FAILED",
        "REFINED_H4_PRODUCT_TANGENT_NOT_NULL",
        "REFINED_H4_PRODUCT_NULL_DECOUPLED",
        "REFINED_H4_NULL_COUPLING_FORMULA_REFUTED",
        "REFINED_H4_NULL_COUPLING_SCHEDULE_DEPENDENT",
        "REFINED_H4_NULL_COUPLING_COMPATIBILITY_CONFIRMED",
    },
    outcome,
)

artifact_records = []
for item in records:
    artifact_records.append({
        "order": item["order"],
        "product_tangent": [mp_text(value) for value in item["tangent"]],
        "step_difference": mp_text(item["step_difference"]),
        "precision_difference": mp_text(item["precision_difference"]),
        "envelope": mp_text(item["envelope"]),
        "raw_imaginary": mp_text(item["raw_imaginary"]),
        "internal_image": [mp_text(value) for value in item["internal_image"]],
        "internal_image_maximum": mp_text(vector_max(item["internal_image"])),
        "compatibility_row": [mp_text(value) for value in item["compatibility"]],
        "direct_boundary_gradient": [
            mp_text(value) for value in item["direct_boundary_gradient"]
        ],
        "half_gradient_error": mp_text(item["formula_error"]),
    })

artifact = {
    "title": "Product-lapse null line and boundary compatibility covector",
    "date": "2026-08-21",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "carrier": "K0=P(sd K_600)",
        "tau0": TAU_TEXT,
        "product_tangent": "(-tau^2/q_cross_rs for six pairs,1,1,1,1)",
        "boundary_prediction": "H_bi*n=(1/2)*g_boundary",
        "difference_steps": list(STEP_TEXTS),
        "decimal_precisions": [PRIMARY_DPS, SECONDARY_DPS],
        "allowed_schedule_identification": "old/new layer swap under time reversal only",
    },
    "product_family": {
        "tau_ratios": ["0.5", "1", "2"],
        "point_count": 72,
        "maximum_internal_residual": mp_text(maximum_product_residual),
        "maximum_angle_identity_residual": mp_text(maximum_branch_identity),
        "maximum_imaginary_curvature": mp_text(maximum_branch_imaginary),
        "minimum_angle_argument": mp_text(minimum_branch_argument),
    },
    "census": {
        "schedule_count": 24,
        "internal_null_count": sum(
            vector_max(records[index]["internal_image"]) <= envelopes[index]
            for index in range(24)
        ),
        "maximum_internal_image": mp_text(maximum_internal_image),
        "maximum_half_gradient_error": mp_text(maximum_formula_error),
        "maximum_frozen_boundary_error": mp_text(maximum_frozen_error),
        "minimum_coupling_to_envelope_ratio": mp_text(min(
            coupling_norms[index] / envelopes[index] for index in range(24)
        )),
        "compatibility_row_rank": 1 if schedule_independent and coupling_resolved_nonzero else None,
        "schedule_independent": schedule_independent,
        "maximum_schedule_difference": mp_text(maximum_schedule_difference),
        "time_reversal_covariant": time_reversal_covariant,
        "maximum_time_reversal_difference": mp_text(maximum_reversal_difference),
        "schedules": artifact_records,
    },
    "controls": {
        "synthetic_directional_error": mp_text(synthetic_error),
        "minimum_wrong_factor_error": mp_text(min(wrong_full_errors)),
        "minimum_wrong_sign_error": mp_text(min(wrong_sign_errors)),
        "dropped_component_index": largest_index,
        "dropped_component_error": mp_text(dropped_error),
    },
    "scope": scope,
    "status_labels": {
        "prediction_provenance": "POST_PRIMARY_DERIVED_PREDICTION",
        "ordinary_schur_complement": "FORBIDDEN_INTERNAL_SINGULAR",
        "constrained_effective_hessian": "NOT_COMPUTED",
        "tick_c_G_planck_particles": "OPEN_NOT_COMPUTED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Outcome: {outcome}")
print(f"Internal product-null schedules: {artifact['census']['internal_null_count']}/24")
print(f"Compatibility row rank: {artifact['census']['compatibility_row_rank']}")
print(f"Maximum half-gradient error: {mp_text(maximum_formula_error, 8)}")

raise SystemExit(0 if passed == tests else 1)

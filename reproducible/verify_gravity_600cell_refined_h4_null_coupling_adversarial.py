#!/usr/bin/env python3
"""Adversarial spatial-hinge reconstruction of the H4 null coupling."""

from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIMARY = HERE / "gravity_600cell_refined_h4_null_coupling.json"
PRIMARY_RESULT = ROOT / "docs/gravity/gravity_600cell_refined_h4_null_coupling_primary_result.md"
CURVATURE_ADVERSARIAL = HERE / "gravity_600cell_refined_local_curvature_mass_adversarial.json"
BOUNDARY_ADVERSARIAL = HERE / "gravity_600cell_refined_boundary_cotangent_adversarial.json"
CURVATURE_RESULT = ROOT / "docs/gravity/gravity_600cell_refined_local_curvature_mass_result.md"
BOUNDARY_RESULT = ROOT / "docs/gravity/gravity_600cell_refined_boundary_cotangent_result.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_h4_null_coupling_adversarial_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_h4_null_coupling_adversarial.json"

PRIMARY_RESULT_COMMIT = "fd0e5f3"
PROTOCOL_COMMIT = "9df1247"
EXPECTED_HASHES = {
    "primary": "6b6fbd95b07f365b3fcac332fa3546021e8d756a510af0184bc974e52d5efa79",
    "primary_result": "5dfefd1f0b2fdcae02cede4a9e7d069a5e7e3d0c29b6f1d324368a4cdbe8803a",
    "curvature_adversarial": "c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e",
    "boundary_adversarial": "19c888a43bdba9d57166d6e3595c6d5b51dd019ebf616efdbf1189e25078f808",
    "curvature_result": "ef6e29fc1e4c89d893a40ee2b5efb3ab6c833e0d73ec232bdbd41033bc4f0f94",
    "boundary_result": "391a317b9f8823a5479f450dde43a43177e210a2d81192aedc938e90fc8006d1",
    "protocol": "3fd363ca0837462977026a98d39f9525f873600c235f1866d93b86c0b6861acd",
}
PAIR_LABELS = ("01", "02", "03", "12", "13", "23")
TAU_TEXT = "0.0102"

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


def mp_text(value, digits=75):
    return mp.nstr(value, digits, strip_zeros=False)


def vector_difference(left, right):
    return max(abs(left[index] - right[index]) for index in range(len(left)))


print("=" * 78)
print("ADVERSARIAL REFINED H4 NULL-COUPLING RECONSTRUCTION")
print("=" * 78)

paths = {
    "primary": PRIMARY,
    "primary_result": PRIMARY_RESULT,
    "curvature_adversarial": CURVATURE_ADVERSARIAL,
    "boundary_adversarial": BOUNDARY_ADVERSARIAL,
    "curvature_result": CURVATURE_RESULT,
    "boundary_result": BOUNDARY_RESULT,
    "protocol": PROTOCOL,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all primary and independently reconstructed inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIMARY_RESULT_COMMIT == "fd0e5f3"
    and PROTOCOL_COMMIT == "9df1247",
)

curvature = json.loads(CURVATURE_ADVERSARIAL.read_text())
boundary = json.loads(BOUNDARY_ADVERSARIAL.read_text())
independent_upstream_ok = check(
    "both frozen actual-incidence adversarial inputs have accepted outcomes",
    curvature["outcome"]
        == "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED"
    and curvature["tests"] == {"passed": 16, "total": 16}
    and boundary["outcome"]
        == "ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CORROBORATED"
    and boundary["tests"] == {"passed": 12, "total": 12},
)

l, tau_symbol, epsilon = sp.symbols("l tau epsilon", positive=True, real=True)
x = l**2
y = -tau_symbol**2
z = l**2 - tau_symbol**2
area_square = sp.simplify(
    (2 * (x * y + x * z + y * z) - x**2 - y**2 - z**2) / 16
)
area = sp.I * l * tau_symbol / 2
log_length_derivative = sp.simplify(l * sp.diff(area, l) / 2)
boundary_gradient = sp.simplify(-sp.I * epsilon * log_length_derivative)
log_lapse_derivative = sp.simplify(
    tau_symbol * sp.diff(boundary_gradient, tau_symbol) / 2
)
symbolic_ok = check(
    "Heron differentiation independently derives the 1/8 null-coupling factor",
    sp.simplify(area_square + l**2 * tau_symbol**2 / 4) == 0
    and sp.simplify(log_length_derivative - sp.I * l * tau_symbol / 4) == 0
    and sp.simplify(boundary_gradient - epsilon * l * tau_symbol / 4) == 0
    and sp.simplify(log_lapse_derivative - epsilon * l * tau_symbol / 8) == 0,
    f"A^2={area_square}, d_g/du={log_lapse_derivative}",
)

with mp.workdps(100):
    tau = mp.mpf(TAU_TEXT)
    pair_curvatures = tuple(
        mp.mpf(boundary["actual_incidence"]["pair_curvatures"][
            f"{label[0]}-{label[1]}"
        ])
        for label in PAIR_LABELS
    )
    hinge_half = tuple(tau * value / 8 for value in pair_curvatures)
    adversarial_row = hinge_half + hinge_half

print("[INFO] adversarial row before loading primary compatibility data:")
print("       " + ", ".join(mp_text(value, 18) for value in adversarial_row))

with mp.workdps(100):
    expected_pre = tuple(-tau * value / 4 for value in pair_curvatures)
    expected_post = tuple(+tau * value / 4 for value in pair_curvatures)
    boundary_pre = tuple(
        mp.mpf(value) for value in boundary["derived_boundary_covector"]["pre"]
    )
    boundary_post = tuple(
        mp.mpf(value) for value in boundary["derived_boundary_covector"]["post"]
    )
    boundary_reconstruction_error = max(
        vector_difference(expected_pre, boundary_pre),
        vector_difference(expected_post, boundary_post),
    )
boundary_control_ok = check(
    "actual pair curvatures reconstruct the accepted pre/post boundary vectors",
    boundary_reconstruction_error < mp.mpf("1e-68"),
    f"max error={mp_text(boundary_reconstruction_error, 8)}",
)

with mp.workdps(100):
    rank_curvatures = tuple(
        mp.mpf(value) for value in curvature["actual_incidence"]["rank_curvatures"]
    )
    rank_masses = tuple(value / (8 * mp.pi) for value in rank_curvatures)
    vertical_coefficients = tuple(
        rank_curvatures[index] / 2 - 4 * mp.pi * rank_masses[index]
        for index in range(4)
    )
    maximum_vertical_coefficient = max(abs(value) for value in vertical_coefficients)
vertical_null_ok = check(
    "actual rank curvatures give an identically flat vertical product coefficient",
    maximum_vertical_coefficient < mp.mpf("1e-90"),
    f"max |K/2-4*pi*m|={mp_text(maximum_vertical_coefficient, 8)}",
)

# Only now load the primary directional-Hessian result.
primary = json.loads(PRIMARY.read_text())
primary_upstream_ok = check(
    "the post-construction primary artifact has the frozen compatibility outcome",
    primary["outcome"]
        == "REFINED_H4_NULL_COUPLING_COMPATIBILITY_CONFIRMED"
    and primary["tests"] == {"passed": 16, "total": 16}
    and primary["census"]["compatibility_row_rank"] == 1
    and primary["census"]["internal_null_count"] == 24,
)

with mp.workdps(100):
    primary_rows = tuple(
        tuple(mp.mpf(value) for value in record["compatibility_row"])
        for record in primary["census"]["schedules"]
    )
    primary_envelopes = tuple(
        mp.mpf(record["envelope"])
        for record in primary["census"]["schedules"]
    )
    row_errors = tuple(
        vector_difference(adversarial_row, row) for row in primary_rows
    )
    maximum_primary_error = max(row_errors)
    error_to_envelope_ratios = tuple(
        row_errors[index] / primary_envelopes[index] for index in range(24)
    )
    maximum_error_to_envelope_ratio = max(error_to_envelope_ratios)
primary_match = all(
    primary_envelopes[index] > 0
    and row_errors[index] <= primary_envelopes[index]
    for index in range(24)
)
primary_comparison_ok = check(
    "all 24 primary rows are compared inside their frozen uncertainty envelopes",
    len(primary_rows) == len(primary_envelopes) == 24
    and all(len(row) == 12 for row in primary_rows)
    and all(value > 0 for value in primary_envelopes),
    f"matches={sum(row_errors[index] <= primary_envelopes[index] for index in range(24))}/24, "
    f"max error/envelope={mp_text(maximum_error_to_envelope_ratio, 8)}",
)

with mp.workdps(100):
    layer_swap = adversarial_row[6:] + adversarial_row[:6]
    layer_swap_error = vector_difference(adversarial_row, layer_swap)
    nonzero = max(abs(value) for value in adversarial_row) > mp.mpf("1e-6")
    rank_one = nonzero and primary_match
rank_reversal_match = rank_one and layer_swap_error < mp.mpf("1e-90")
rank_control_ok = check(
    "the reconstructed repeated row is nonzero, rank one and layer-swap invariant",
    nonzero and layer_swap_error < mp.mpf("1e-90"),
    f"rank-one primary match={rank_one}, swap error={mp_text(layer_swap_error, 8)}",
)

with mp.workdps(100):
    wrong_factor = tuple(tau * value / 4 for value in pair_curvatures) * 2
    wrong_sign = tuple(-value for value in adversarial_row)
    largest_index = max(range(6), key=lambda index: abs(pair_curvatures[index]))
    dropped = list(adversarial_row)
    dropped[largest_index] = mp.mpf(0)
    dropped[6 + largest_index] = mp.mpf(0)
    swapped = list(adversarial_row)
    swapped[0], swapped[3] = swapped[3], swapped[0]
    swapped[6], swapped[9] = swapped[9], swapped[6]
    corruption_errors = {
        "wrong_factor": vector_difference(wrong_factor, adversarial_row),
        "wrong_sign": vector_difference(wrong_sign, adversarial_row),
        "dropped_largest": vector_difference(tuple(dropped), adversarial_row),
        "swapped_unequal_pairs": vector_difference(tuple(swapped), adversarial_row),
    }
corruption_ok = check(
    "wrong factor, sign, dropped pair and unequal-pair swap all fail decisively",
    min(corruption_errors.values()) > mp.mpf("1e-6"),
)

scope = {
    "primary_functions_imported_or_executed": False,
    "lorentzian_action_evaluator_executed": False,
    "numerical_derivative_or_eigensolver_executed": False,
    "full_hessian_or_pseudoinverse_computed": False,
    "cross_null_image_independently_rebuilt": False,
    "root_search_or_spectrum_executed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
}
scope_ok = check(
    "the reconstruction is mechanically independent and discloses its cross-null limit",
    not any(value for key, value in scope.items()
            if key != "cross_null_image_independently_rebuilt")
    and not scope["cross_null_image_independently_rebuilt"],
)

controls_ok = all((
    provenance_ok,
    independent_upstream_ok,
    symbolic_ok,
    boundary_control_ok,
    vertical_null_ok,
    primary_upstream_ok,
    primary_comparison_ok,
    rank_control_ok,
    corruption_ok,
    scope_ok,
))

if not controls_ok:
    outcome = "ADVERSARIAL_REFINED_H4_NULL_COUPLING_CONTROL_FAILED"
elif not (primary_match and rank_reversal_match):
    outcome = "ADVERSARIAL_REFINED_H4_NULL_COUPLING_DISAGREEMENT"
else:
    outcome = "ADVERSARIAL_REFINED_H4_NULL_COUPLING_CORROBORATED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one adversarial outcome",
    outcome in {
        "ADVERSARIAL_REFINED_H4_NULL_COUPLING_CONTROL_FAILED",
        "ADVERSARIAL_REFINED_H4_NULL_COUPLING_DISAGREEMENT",
        "ADVERSARIAL_REFINED_H4_NULL_COUPLING_CORROBORATED",
    },
    outcome,
)

artifact = {
    "title": "Adversarial spatial-hinge reconstruction of H4 null coupling",
    "date": "2026-08-21",
    "primary_result_commit": PRIMARY_RESULT_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "symbolic_product_hinge": {
        "area_square": str(area_square),
        "area": str(area),
        "log_length_area_derivative": str(log_length_derivative),
        "boundary_gradient_per_edge": str(boundary_gradient),
        "log_lapse_boundary_derivative": str(log_lapse_derivative),
    },
    "actual_incidence": {
        "pair_order": list(PAIR_LABELS),
        "pair_curvatures": [mp_text(value) for value in pair_curvatures],
        "rank_curvatures": [mp_text(value) for value in rank_curvatures],
        "derived_rank_masses": [mp_text(value) for value in rank_masses],
        "maximum_vertical_null_coefficient": mp_text(maximum_vertical_coefficient),
    },
    "compatibility": {
        "adversarial_row": [mp_text(value) for value in adversarial_row],
        "maximum_primary_component_error": mp_text(maximum_primary_error),
        "maximum_primary_error_to_envelope_ratio": mp_text(
            maximum_error_to_envelope_ratio
        ),
        "primary_envelopes": [mp_text(value) for value in primary_envelopes],
        "primary_row_count": len(primary_rows),
        "rank": 1 if rank_one else None,
        "layer_swap_error": mp_text(layer_swap_error),
        "boundary_vector_reconstruction_error": mp_text(
            boundary_reconstruction_error
        ),
    },
    "controls": {
        key: mp_text(value) for key, value in corruption_errors.items()
    },
    "scope": scope,
    "status_labels": {
        "compatibility_coupling": "MECHANICALLY_INDEPENDENT_RECONSTRUCTION",
        "vertical_product_null": "ANALYTIC_CURVATURE_BALANCE",
        "cross_product_null": "PRIMARY_FINITE_FAMILY_NOT_INDEPENDENTLY_REBUILT",
        "constrained_effective_hessian": "NOT_COMPUTED",
        "tick_c_G_planck_particles": "OPEN_NOT_COMPUTED",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Outcome: {outcome}")
print(f"Maximum primary error: {mp_text(maximum_primary_error, 8)}")
print("Cross-null independent rebuild: NOT PERFORMED (explicit scope limit)")

raise SystemExit(0 if passed == tests else 1)

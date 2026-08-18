#!/usr/bin/env python3
"""Target-disclosed serialization audit for shifted principal coefficients.

Prior-art commit: 52f337b.
Protocol commit: 86a9129.
No stiffness eigenvalue or sign field is accessed.
"""

from collections import Counter
import hashlib
import json
import math
from pathlib import Path

from flint import arb, acb, acb_mat, ctx
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
TWO_STEP_JSON = HERE / "gravity_600cell_dust_two_step_full_tangent.json"
TWO_STEP_NPZ = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"
TWO_STEP_SOURCE = HERE / "verify_gravity_600cell_dust_two_step_full_tangent.py"
SHIFTED_JACOBI_JSON = HERE / "gravity_600cell_dust_shifted_jacobi.json"
SHIFTED_JACOBI_NPZ = HERE / "gravity_600cell_dust_shifted_jacobi.npz"
SHIFTED_JACOBI_SOURCE = HERE / "verify_gravity_600cell_dust_shifted_jacobi.py"
SHIFTED_CENTERED_JSON = HERE / "gravity_600cell_dust_shifted_centered.json"
SHIFTED_CENTERED_NPZ = HERE / "gravity_600cell_dust_shifted_centered.npz"
SHIFTED_CENTERED_SOURCE = HERE / "verify_gravity_600cell_dust_shifted_centered.py"
SHIFTED_SHAPE_JSON = HERE / "gravity_600cell_dust_shifted_shape_stiffness.json"
SHIFTED_SHAPE_SOURCE = HERE / "verify_gravity_600cell_dust_shifted_shape_stiffness.py"
OUTPUT = HERE / "gravity_600cell_dust_shifted_precision_audit.json"

PRIOR_ART_COMMIT = "52f337b"
PROTOCOL_COMMIT = "86a9129"
EXPECTED_HASHES = {
    "two_step_json": "f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc",
    "two_step_npz": "ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d",
    "two_step_source": "c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717",
    "shifted_jacobi_json": "63b37b6000146d5d53dbbc01da5c9aba9a5e3373b8bc3830a404ef0f681ecf31",
    "shifted_jacobi_npz": "d2f507c4a2fa11c5d7a808849c199a986278516f422cf43654f6de153ab170d0",
    "shifted_jacobi_source": "fc070de1a5f89524f119fd09ae25611559cb0bb64defbd34486b7980f058470d",
    "shifted_centered_json": "265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47",
    "shifted_centered_npz": "c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8",
    "shifted_centered_source": "a3c45e3e636057d83a663d3248dd023f7d04ec6e544c698f9116307822be337a",
    "shifted_shape_json": "14fe5bc91e3ae4712c6ea19b8120785e2facd364e1ceb194009123fa353a4315",
    "shifted_shape_source": "031d0dd1cab45d0093015fcab7ce7b56e098a5742895eed71f5a531aee31c2a6",
}
PARITIES = ("even", "odd")
SECTORS = (4, 5)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
SIZE = 30
MACHINE_EPSILON = np.finfo(float).eps
ctx.dps = 80
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


def sf(value):
    return f"{float(value):.17e}"


def radius_ball(midpoint, radius):
    center = arb(format(float(midpoint), ".17g"))
    if radius <= 0:
        return center
    return center + arb(0, format(float(radius), ".17g"))


def reenclose_binary_matrix(midpoint, radii, include_half_ulp):
    rows, columns = midpoint.shape
    matrix = acb_mat(rows, columns)
    component_radii = np.empty((rows, columns), dtype=float)
    for row in range(rows):
        for column in range(columns):
            value = complex(midpoint[row, column])
            stored = float(radii[row, column])
            real_extra = 0.5 * abs(float(np.spacing(value.real))) if include_half_ulp else 0.0
            imag_extra = 0.5 * abs(float(np.spacing(value.imag))) if include_half_ulp else 0.0
            real_radius = stored + real_extra
            imag_radius = stored + imag_extra
            component_radii[row, column] = math.hypot(real_radius, imag_radius)
            matrix[row, column] = acb(
                radius_ball(value.real, real_radius),
                radius_ball(value.imag, imag_radius),
            )
    return matrix, component_radii


def submatrix(matrix, rows, columns):
    result = acb_mat(len(rows), len(columns))
    for out_row, row in enumerate(rows):
        for out_column, column in enumerate(columns):
            result[out_row, out_column] = matrix[row, column]
    return result


def identity(size):
    result = acb_mat(size, size)
    for index in range(size):
        result[index, index] = 1
    return result


def split_tangent(matrix):
    q = list(range(SIZE))
    p = list(range(SIZE, 2 * SIZE))
    return {
        "A": submatrix(matrix, q, q),
        "B": submatrix(matrix, q, p),
        "C": submatrix(matrix, p, q),
        "D": submatrix(matrix, p, p),
    }


def principal_s10(blocks):
    inverse_b = blocks["B"].solve(identity(SIZE))
    s11 = blocks["D"] * inverse_b
    return blocks["C"] - s11 * blocks["A"], blocks["B"].det()


def acb_midpoint_and_radii(matrix):
    midpoint = np.empty((matrix.nrows(), matrix.ncols()), dtype=np.complex128)
    radii = np.empty((matrix.nrows(), matrix.ncols()), dtype=float)
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            value = matrix[row, column]
            midpoint[row, column] = complex(
                float(value.real.mid()), float(value.imag.mid())
            )
            radii[row, column] = math.hypot(
                float(value.real.rad().upper()),
                float(value.imag.rad().upper()),
            )
    return midpoint, radii


def component_reenclosure_radii(midpoint, stored):
    real_half_ulp = 0.5 * np.abs(np.spacing(np.real(midpoint)))
    imaginary_half_ulp = 0.5 * np.abs(np.spacing(np.imag(midpoint)))
    return np.hypot(stored + real_half_ulp, stored + imaginary_half_ulp)


def operator_norm(matrix):
    values = la.svdvals(matrix)
    return float(values[0]) if len(values) else 0.0


def variant_record(payload, parity, sector_index, variant):
    sector = payload["parities"][parity][sector_index]
    matches = [item for item in sector["variants"] if item["variant"] == variant]
    if len(matches) != 1:
        raise RuntimeError("variant lookup is not unique")
    return matches[0]


def classification(serial_ratio, coefficient_ratio):
    if serial_ratio > 100 and coefficient_ratio > 10:
        return "SERIALIZATION_DOMINANT_RESOLVED"
    if serial_ratio < 10 or coefficient_ratio < 1:
        return "SERIALIZATION_NOT_DOMINANT"
    return "SERIALIZATION_MIXED"


print("=" * 78)
print("SHIFTED PRINCIPAL-FUNCTION BINARY SERIALIZATION AUDIT")
print("=" * 78)

paths = {
    "two_step_json": TWO_STEP_JSON,
    "two_step_npz": TWO_STEP_NPZ,
    "two_step_source": TWO_STEP_SOURCE,
    "shifted_jacobi_json": SHIFTED_JACOBI_JSON,
    "shifted_jacobi_npz": SHIFTED_JACOBI_NPZ,
    "shifted_jacobi_source": SHIFTED_JACOBI_SOURCE,
    "shifted_centered_json": SHIFTED_CENTERED_JSON,
    "shifted_centered_npz": SHIFTED_CENTERED_NPZ,
    "shifted_centered_source": SHIFTED_CENTERED_SOURCE,
    "shifted_shape_json": SHIFTED_SHAPE_JSON,
    "shifted_shape_source": SHIFTED_SHAPE_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
two_step_json = json.loads(TWO_STEP_JSON.read_text())
shifted_jacobi_json = json.loads(SHIFTED_JACOBI_JSON.read_text())
shifted_centered_json = json.loads(SHIFTED_CENTERED_JSON.read_text())
shifted_shape = json.loads(SHIFTED_SHAPE_JSON.read_text())
two_step = np.load(TWO_STEP_NPZ, allow_pickle=False)
shifted_jacobi = np.load(SHIFTED_JACOBI_NPZ, allow_pickle=False)
shifted_centered = np.load(SHIFTED_CENTERED_NPZ, allow_pickle=False)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and two_step_json["outcome"] == "TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED"
    and two_step_json["passed"] == two_step_json["tests"] == 16
    and len(two_step.files) == 448
    and shifted_jacobi_json["outcome"] == "SHIFTED_JACOBI_CERTIFIED"
    and shifted_jacobi_json["passed"] == shifted_jacobi_json["tests"] == 8
    and len(shifted_jacobi.files) == 560
    and shifted_centered_json["outcome"] == "SHIFTED_CENTERED_CERTIFIED"
    and shifted_centered_json["passed"] == shifted_centered_json["tests"] == 7
    and len(shifted_centered.files) == 560
    and shifted_shape["outcome"] == "SHIFTED_SHAPE_STIFFNESS_SIGN_OPEN"
    and shifted_shape["passed"] == shifted_shape["tests"] == 12
)
check("all frozen precision-audit inputs have exact provenance", provenance_ok, str(hashes))

own_source = Path(__file__).read_text()
forbidden_keys = (
    "A_" + "eigenvalues",
    "A_" + "sign_labels",
    "sign_" + "counts",
)
target_fields_absent = all(key not in own_source for key in forbidden_keys)
check("the audit source contains no stiffness-spectrum or sign-field access", target_fields_absent)

records = []
labels = Counter()
all_determinants = True
all_reproduced = True
all_downstream_sums = True

for parity in PARITIES:
    for sector_index in SECTORS:
        for variant in VARIANTS:
            tangent_prefix = f"{parity}_sector{sector_index}_t2_{variant}"
            tangent_midpoint = np.asarray(two_step[f"{tangent_prefix}_midpoint"])
            tangent_stored = np.asarray(two_step[f"{tangent_prefix}_radii"])
            certified_ball, certified_input_components = reenclose_binary_matrix(
                tangent_midpoint, tangent_stored, True
            )
            counterfactual_ball, counterfactual_input_components = reenclose_binary_matrix(
                tangent_midpoint, tangent_stored, False
            )
            certified_s10, certified_det = principal_s10(split_tangent(certified_ball))
            counterfactual_s10, counterfactual_det = principal_s10(
                split_tangent(counterfactual_ball)
            )
            all_determinants &= bool(
                not certified_det.contains(0)
                and not counterfactual_det.contains(0)
            )
            certified_midpoint, certified_radius = acb_midpoint_and_radii(certified_s10)
            _, counterfactual_radius = acb_midpoint_and_radii(counterfactual_s10)

            coefficient_prefix = f"{parity}_sector{sector_index}_{variant}"
            stored_kminus_midpoint = np.asarray(
                shifted_jacobi[f"{coefficient_prefix}_Kminus_midpoint"]
            )
            stored_kminus_radius = np.asarray(
                shifted_jacobi[f"{coefficient_prefix}_Kminus_radii"]
            )
            midpoint_reproduced = np.array_equal(
                certified_midpoint, stored_kminus_midpoint
            )
            radius_reproduced = np.allclose(
                certified_radius, stored_kminus_radius,
                rtol=1e-14, atol=0.0,
            )
            all_reproduced &= bool(midpoint_reproduced and radius_reproduced)

            certified_output_radius = float(la.norm(certified_radius, "fro"))
            counterfactual_output_radius = float(
                la.norm(counterfactual_radius, "fro")
            )
            serial_ratio = (
                certified_output_radius / counterfactual_output_radius
                if counterfactual_output_radius > 0 else math.inf
            )
            coefficient_radii = {}
            for name in ("Kminus", "Kzero", "Kplus"):
                coefficient_radii[name] = float(la.norm(
                    shifted_jacobi[f"{coefficient_prefix}_{name}_radii"], "fro"
                ))
            coefficient_ratio = coefficient_radii["Kminus"] / (
                coefficient_radii["Kzero"] + coefficient_radii["Kplus"]
            )
            label = classification(serial_ratio, coefficient_ratio)
            labels[label] += 1

            centered_v_midpoint = np.asarray(
                shifted_centered[f"{coefficient_prefix}_V_midpoint"]
            )
            centered_v_stored = np.asarray(
                shifted_centered[f"{coefficient_prefix}_V_radii"]
            )
            centered_v_radius = component_reenclosure_radii(
                centered_v_midpoint, centered_v_stored
            )
            h_v = (centered_v_midpoint + centered_v_midpoint.conj().T) / 2
            radius_hv = (centered_v_radius + centered_v_radius.T) / 2
            norm_hv = operator_norm(h_v)
            source_radius = float(la.norm(radius_hv, "fro"))
            arithmetic_hv = float(
                1000 * MACHINE_EPSILON * SIZE * max(1.0, norm_hv)
            )
            epsilon_hv = source_radius + arithmetic_hv
            stiffness_record = variant_record(
                shifted_shape, parity, sector_index, variant
            )
            eta_s = float(
                stiffness_record["shape_carrier"]["subspace_error_eta_S"]
            )
            carrier_lift = float(2 * eta_s * (norm_hv + epsilon_hv))
            arithmetic_restriction = arithmetic_hv
            reconstructed_restricted = float(
                epsilon_hv + carrier_lift + arithmetic_restriction
            )
            stored_restricted = float(
                stiffness_record["hermitian_pencil"]["restricted_A_error"]
            )
            downstream_tolerance = float(
                100 * MACHINE_EPSILON * max(1.0, stored_restricted)
            )
            downstream_agreement = bool(
                abs(reconstructed_restricted - stored_restricted)
                <= downstream_tolerance
            )
            all_downstream_sums &= downstream_agreement

            total = reconstructed_restricted
            records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "certified_B_determinant_excludes_zero": bool(
                    not certified_det.contains(0)
                ),
                "counterfactual_B_determinant_excludes_zero": bool(
                    not counterfactual_det.contains(0)
                ),
                "certified_input_component_radius_frobenius": sf(
                    la.norm(certified_input_components, "fro")
                ),
                "stored_only_input_component_radius_frobenius": sf(
                    la.norm(counterfactual_input_components, "fro")
                ),
                "certified_S10_radius_frobenius": sf(certified_output_radius),
                "stored_only_S10_radius_frobenius": sf(
                    counterfactual_output_radius
                ),
                "serialization_ratio": sf(serial_ratio),
                "committed_Kminus_reproduced": bool(
                    midpoint_reproduced and radius_reproduced
                ),
                "coefficient_radius_frobenius": {
                    name: sf(value) for name, value in coefficient_radii.items()
                },
                "Kminus_to_other_coefficient_radius_ratio": sf(
                    coefficient_ratio
                ),
                "classification": label,
                "downstream_restricted_error": {
                    "source_radius": sf(source_radius),
                    "first_arithmetic_floor": sf(arithmetic_hv),
                    "epsilon_HV": sf(epsilon_hv),
                    "carrier_lift": sf(carrier_lift),
                    "final_arithmetic_floor": sf(arithmetic_restriction),
                    "reconstructed_total": sf(reconstructed_restricted),
                    "stored_total": sf(stored_restricted),
                    "agreement_tolerance": sf(downstream_tolerance),
                    "agrees": downstream_agreement,
                    "fractions": {
                        "epsilon_HV": sf(epsilon_hv / total),
                        "carrier_lift": sf(carrier_lift / total),
                        "final_arithmetic": sf(arithmetic_restriction / total),
                    },
                },
            })

check("all 32 declared boundary-twist determinant balls exclude zero", all_determinants)
check("all 16 certified reconstructions reproduce committed Kminus balls", all_reproduced)
check(
    "all 16 downstream restricted-error sums reproduce independently",
    all_downstream_sums,
)
complete = len(records) == 16 and sum(labels.values()) == 16
check("all 16 serialization audits receive frozen classifications", complete, str(dict(labels)))

if not (
    provenance_ok
    and target_fields_absent
    and all_determinants
    and all_reproduced
    and all_downstream_sums
    and complete
):
    outcome = "SHIFTED_PRECISION_AUDIT_CONTROL_FAILED"
elif labels["SERIALIZATION_NOT_DOMINANT"]:
    outcome = "SHIFTED_PRECISION_SERIALIZATION_NOT_DOMINANT"
elif labels["SERIALIZATION_MIXED"]:
    outcome = "SHIFTED_PRECISION_ATTRIBUTION_MIXED"
else:
    outcome = "SHIFTED_PRECISION_BINARY_SERIALIZATION_DOMINANT"

allowed = {
    "SHIFTED_PRECISION_AUDIT_CONTROL_FAILED",
    "SHIFTED_PRECISION_SERIALIZATION_NOT_DOMINANT",
    "SHIFTED_PRECISION_ATTRIBUTION_MIXED",
    "SHIFTED_PRECISION_BINARY_SERIALIZATION_DOMINANT",
}
check("the preregistered hierarchy assigns the precision outcome", outcome in allowed, outcome)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "cells": len(records),
    "classification_counts": dict(labels),
    "records": records,
    "counterfactual_warning": (
        "stored-ball-only omits binary midpoint serialization uncertainty "
        "and is not a rigorous source enclosure"
    ),
    "spectrum_or_sign_field_loaded": False,
    "direct_high_precision_reconstruction_authorized": bool(
        outcome == "SHIFTED_PRECISION_BINARY_SERIALIZATION_DOMINANT"
    ),
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

serial_values = [float(item["serialization_ratio"]) for item in records]
coefficient_values = [
    float(item["Kminus_to_other_coefficient_radius_ratio"])
    for item in records
]
print("-" * 78)
print("SCIENTIFIC OUTCOME:", outcome)
print("classification counts:", dict(labels))
print("serialization ratio:", min(serial_values), "...", max(serial_values))
print("Kminus/other radius ratio:", min(coefficient_values), "...", max(coefficient_values))
print(f"{passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Blind reconstruction of the shifted three-slice 600-cell Jacobi operator.

Prior-art commit: 920ce5d.
Preregistered protocol commit: 9ee1f3f.
No spatial spectrum, desired degeneracy, speed or continuum target is loaded.
"""

from collections import Counter
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

from flint import arb, acb, acb_mat, ctx
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
SECOND_JSON = HERE / "gravity_600cell_dust_two_step_full_tangent.json"
SECOND_NPZ = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"
SECOND_SOURCE = HERE / "verify_gravity_600cell_dust_two_step_full_tangent.py"
THIRD_JSON = HERE / "gravity_600cell_dust_later_slab_tangent.json"
THIRD_NPZ = HERE / "gravity_600cell_dust_later_slab_tangent.npz"
THIRD_SOURCE = HERE / "verify_gravity_600cell_dust_later_slab_tangent.py"
OUTPUT = HERE / "gravity_600cell_dust_shifted_jacobi.json"
NUMERIC_OUTPUT = HERE / "gravity_600cell_dust_shifted_jacobi.npz"

PRIOR_ART_COMMIT = "920ce5d"
PROTOCOL_COMMIT = "9ee1f3f"
EXPECTED_HASHES = {
    "second_json": "f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc",
    "second_npz": "ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d",
    "second_source": "c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717",
    "third_json": "58a95d90d569b25a3aa396346f5198472c8aed706846b4057182b04c9f7480c4",
    "third_npz": "77b4dd54a5dcba9d1aa12870b361c9d7d7dde11ccaaa558361b9c1dc24768196",
    "third_source": "3740622494b88bafc55839e7a53d9c4a335d4783503d236b0f2e9637429d9eea",
}
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
DIMENSIONS = (3, 2, 2, 2, 1, 1, 1)
POSITION_DIMS = tuple(30 * value for value in DIMENSIONS)
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


def serialize_float(value):
    return f"{float(value):.17e}"


def radius_ball(midpoint, radius):
    center = arb(format(float(midpoint), ".17g"))
    if radius <= 0:
        return center
    return center + arb(0, format(float(radius), ".17g"))


def reenclose_binary_matrix(midpoint, radii):
    rows, columns = midpoint.shape
    matrix = acb_mat(rows, columns)
    for row in range(rows):
        for column in range(columns):
            value = complex(midpoint[row, column])
            stored = float(radii[row, column])
            real_radius = stored + 0.5 * abs(float(np.spacing(value.real)))
            imag_radius = stored + 0.5 * abs(float(np.spacing(value.imag)))
            matrix[row, column] = acb(
                radius_ball(value.real, real_radius),
                radius_ball(value.imag, imag_radius),
            )
    return matrix


def identity(size):
    result = acb_mat(size, size)
    for index in range(size):
        result[index, index] = 1
    return result


def submatrix(matrix, rows, columns):
    result = acb_mat(len(rows), len(columns))
    for out_row, row in enumerate(rows):
        for out_column, column in enumerate(columns):
            result[out_row, out_column] = matrix[row, column]
    return result


def negate(matrix):
    result = acb_mat(matrix.nrows(), matrix.ncols())
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            result[row, column] = -matrix[row, column]
    return result


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


def residual_stats(matrix):
    midpoint, radii = acb_midpoint_and_radii(matrix)
    contains_zero = all(
        matrix[row, column].contains(0)
        for row in range(matrix.nrows())
        for column in range(matrix.ncols())
    )
    return {
        "contains_zero_entrywise": contains_zero,
        "midpoint_frobenius": serialize_float(la.norm(midpoint, "fro")),
        "radius_frobenius": serialize_float(la.norm(radii, "fro")),
        "maximum_midpoint_modulus": serialize_float(np.max(np.abs(midpoint))),
        "maximum_radius": serialize_float(np.max(radii)),
    }


def determinant_record(value):
    return {
        "contains_zero": bool(value.contains(0)),
        "ball": str(value),
    }


def deterministic_npz(path, arrays):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            np.lib.format.write_array(
                buffer, np.asarray(arrays[name]), allow_pickle=False
            )
            info = zipfile.ZipInfo(
                f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0)
            )
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            archive.writestr(info, buffer.getvalue())
    temporary.replace(path)


def split_tangent(matrix, size):
    q = list(range(size))
    p = list(range(size, 2 * size))
    return {
        "A": submatrix(matrix, q, q),
        "B": submatrix(matrix, q, p),
        "C": submatrix(matrix, p, q),
        "D": submatrix(matrix, p, p),
    }


def reconstruct_principal(blocks):
    inv_b = blocks["B"].solve(identity(blocks["B"].nrows()))
    s00 = inv_b * blocks["A"]
    s01 = negate(inv_b)
    s11 = blocks["D"] * inv_b
    s10 = blocks["C"] - s11 * blocks["A"]
    residuals = {
        "S00_adjoint": s00 - s00.transpose().conjugate(),
        "S11_adjoint": s11 - s11.transpose().conjugate(),
        "S10_S01_adjoint": s10 - s01.transpose().conjugate(),
        "recover_A": blocks["A"] - blocks["B"] * s00,
        "recover_C": blocks["C"] - (s10 + s11 * blocks["A"]),
        "recover_D": blocks["D"] - s11 * blocks["B"],
    }
    return {"00": s00, "01": s01, "10": s10, "11": s11}, residuals


def concatenate_numpy(matrices):
    return np.hstack(matrices)


def matrix_family_uncertainty(midpoints, radii, names):
    primary = concatenate_numpy([midpoints[name][VARIANTS[0]] for name in names])
    variation = max(
        la.norm(
            concatenate_numpy([midpoints[name][variant] for name in names])
            - primary,
            2,
        )
        for variant in VARIANTS
    )
    radius = max(
        la.norm(
            concatenate_numpy([radii[name][variant] for name in names]), "fro"
        )
        for variant in VARIANTS
    )
    singular = la.svdvals(primary)
    floor = 10 * np.finfo(float).eps * max(1.0, float(singular[0]))
    return {
        "midpoint": primary,
        "singular": singular,
        "variation": float(variation),
        "radius": float(radius),
        "floor": float(floor),
        "epsilon": float(variation + radius + floor),
    }


def schedule_label(distance, epsilon):
    if not math.isfinite(epsilon):
        return "SCHEDULE_OPEN"
    if distance <= 10 * epsilon:
        return "SCHEDULE_ROBUST"
    if distance > 100 * epsilon:
        return "SCHEDULE_DEPENDENT"
    return "SCHEDULE_OPEN"


hashes = {
    "second_json": sha256(SECOND_JSON),
    "second_npz": sha256(SECOND_NPZ),
    "second_source": sha256(SECOND_SOURCE),
    "third_json": sha256(THIRD_JSON),
    "third_npz": sha256(THIRD_NPZ),
    "third_source": sha256(THIRD_SOURCE),
}
second_json = json.loads(SECOND_JSON.read_text())
second_npz = np.load(SECOND_NPZ)
third_json = json.loads(THIRD_JSON.read_text())
third_npz = np.load(THIRD_NPZ)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and second_json["outcome"] == "TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED"
    and second_json["passed"] == second_json["tests"] == 16
    and second_json["numeric_archive_arrays"] == len(second_npz.files) == 448
    and second_json["numeric_archive_sha256"] == EXPECTED_HASHES["second_npz"]
    and third_json["outcome"] == "LATER_SLAB_TANGENT_CERTIFIED"
    and third_json["passed"] == third_json["tests"] == 16
    and third_json["numeric_archive_arrays"] == len(third_npz.files) == 448
    and third_json["numeric_archive_sha256"] == EXPECTED_HASHES["third_npz"]
)
check("all target-free tangent inputs have exact frozen provenance",
      provenance_ok, str(hashes))

carrier_ok = all(
    tuple(
        item["dimension"] for item in second_json["parities"][parity]["sectors"]
    ) == DIMENSIONS
    and tuple(
        item["dimension"] for item in third_json["parities"][parity]["sectors"]
    ) == DIMENSIONS
    for parity in ("even", "odd")
)
check("both archives use the frozen seven-sector 720-position carrier",
      carrier_ok, f"positions={sum(d*d*30 for d in DIMENSIONS)}")

print("=" * 78)
print("BLIND SHIFTED THREE-SLICE REGGE JACOBI RECONSTRUCTION")
print("=" * 78)

records = {"even": [], "odd": []}
numeric_arrays = {}
all_twists = True
all_principal_identities = True
all_product_identities = True
twist_count = 0

for parity in ("even", "odd"):
    for sector_index, (dimension, size) in enumerate(zip(DIMENSIONS, POSITION_DIMS)):
        print(f"[{parity}] sector {sector_index + 1}/7 n={size}", flush=True)
        variant_records = {}
        matrix_midpoints = {
            name: {} for name in ("Kminus", "Kzero", "Kplus", "P", "Q")
        }
        matrix_radii = {
            name: {} for name in ("Kminus", "Kzero", "Kplus", "P", "Q")
        }
        b_midpoints = {"B2": {}, "B3": {}}
        b_radii = {"B2": {}, "B3": {}}

        for variant in VARIANTS:
            prefix2 = f"{parity}_sector{sector_index}"
            t2 = reenclose_binary_matrix(
                second_npz[f"{prefix2}_t2_{variant}_midpoint"],
                second_npz[f"{prefix2}_t2_{variant}_radii"],
            )
            t3 = reenclose_binary_matrix(
                third_npz[f"{prefix2}_t3_{variant}_midpoint"],
                third_npz[f"{prefix2}_t3_{variant}_radii"],
            )
            product = reenclose_binary_matrix(
                third_npz[f"{prefix2}_product_{variant}_midpoint"],
                third_npz[f"{prefix2}_product_{variant}_radii"],
            )
            two = split_tangent(t2, size)
            three = split_tangent(t3, size)
            composed = split_tangent(product, size)

            det2 = two["B"].det()
            det3 = three["B"].det()
            twist_count += 2
            twist_ok = not det2.contains(0) and not det3.contains(0)
            all_twists &= twist_ok

            principal2, residual2 = reconstruct_principal(two)
            principal3, residual3 = reconstruct_principal(three)
            kminus = principal2["10"]
            kzero = principal2["11"] + principal3["00"]
            kplus = principal3["01"]
            p_matrix = kplus.solve(negate(kzero))
            q_matrix = kplus.solve(negate(kminus))

            seam_residuals = {
                "implicit_q1": (
                    kminus + kzero * two["A"] + kplus * composed["A"]
                ),
                "implicit_p1": kzero * two["B"] + kplus * composed["B"],
                "solved_q1": p_matrix * two["A"] + q_matrix - composed["A"],
                "solved_p1": p_matrix * two["B"] - composed["B"],
            }
            principal_stats = {
                f"slab2_{name}": residual_stats(value)
                for name, value in residual2.items()
            }
            principal_stats.update({
                f"slab3_{name}": residual_stats(value)
                for name, value in residual3.items()
            })
            seam_stats = {
                name: residual_stats(value)
                for name, value in seam_residuals.items()
            }
            principal_ok = all(
                item["contains_zero_entrywise"]
                for item in principal_stats.values()
            )
            product_ok = all(
                item["contains_zero_entrywise"] for item in seam_stats.values()
            )
            all_principal_identities &= principal_ok
            all_product_identities &= product_ok

            matrices = {
                "Kminus": kminus,
                "Kzero": kzero,
                "Kplus": kplus,
                "P": p_matrix,
                "Q": q_matrix,
            }
            for name, matrix in matrices.items():
                midpoint, radii = acb_midpoint_and_radii(matrix)
                matrix_midpoints[name][variant] = midpoint
                matrix_radii[name][variant] = radii
                prefix = f"{parity}_sector{sector_index}_{variant}_{name}"
                numeric_arrays[f"{prefix}_midpoint"] = midpoint
                numeric_arrays[f"{prefix}_radii"] = radii

            for name, matrix in (("B2", two["B"]), ("B3", three["B"])):
                midpoint, radii = acb_midpoint_and_radii(matrix)
                b_midpoints[name][variant] = midpoint
                b_radii[name][variant] = radii

            variant_records[variant] = {
                "twist_determinants": {
                    "B2": determinant_record(det2),
                    "B3": determinant_record(det3),
                },
                "principal_residuals": principal_stats,
                "product_residuals": seam_stats,
            }

        family_k = matrix_family_uncertainty(
            matrix_midpoints, matrix_radii, ("Kminus", "Kzero", "Kplus")
        )
        family_pq = matrix_family_uncertainty(
            matrix_midpoints, matrix_radii, ("P", "Q")
        )
        b_diagnostics = {}
        for name in ("B2", "B3"):
            primary = b_midpoints[name][VARIANTS[0]]
            singular = la.svdvals(primary)
            variation = max(
                la.norm(b_midpoints[name][variant] - primary, 2)
                for variant in VARIANTS
            )
            radius = max(
                la.norm(b_radii[name][variant], "fro") for variant in VARIANTS
            )
            epsilon = variation + radius + (
                10 * np.finfo(float).eps * max(1.0, float(singular[0]))
            )
            b_diagnostics[name] = {
                "minimum_singular": serialize_float(singular[-1]),
                "maximum_singular": serialize_float(singular[0]),
                "condition": serialize_float(singular[0] / singular[-1]),
                "epsilon_singular": serialize_float(epsilon),
            }

        primary_matrices = {
            name: matrix_midpoints[name][VARIANTS[0]]
            for name in matrix_midpoints
        }
        denominator = max(
            la.norm(primary_matrices["Kplus"], "fro"),
            la.norm(primary_matrices["Kminus"], "fro"),
        )
        asymmetry = la.norm(
            primary_matrices["Kplus"]
            - primary_matrices["Kminus"].conj().T,
            "fro",
        ) / denominator
        records[parity].append({
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "position_dimension": size,
            "variants": variant_records,
            "twist_diagnostics": b_diagnostics,
            "matrix_frobenius_norms": {
                name: serialize_float(la.norm(matrix, "fro"))
                for name, matrix in primary_matrices.items()
            },
            "background_asymmetry": serialize_float(asymmetry),
            "K_family": {
                "singular_values": [serialize_float(value)
                                    for value in family_k["singular"]],
                "variation": serialize_float(family_k["variation"]),
                "radius": serialize_float(family_k["radius"]),
                "floor": serialize_float(family_k["floor"]),
                "epsilon": serialize_float(family_k["epsilon"]),
            },
            "PQ_family": {
                "singular_values": [serialize_float(value)
                                    for value in family_pq["singular"]],
                "variation": serialize_float(family_pq["variation"]),
                "radius": serialize_float(family_pq["radius"]),
                "floor": serialize_float(family_pq["floor"]),
                "epsilon": serialize_float(family_pq["epsilon"]),
            },
        })

check("all 112 boundary-twist determinant balls exclude zero",
      all_twists, f"determinants={twist_count}")
check("all reconstructed principal-function identities hold entrywise",
      all_principal_identities)
check("all implicit and solved product-equivalence identities hold entrywise",
      all_product_identities)

archive_ok = len(numeric_arrays) == 560
deterministic_npz(NUMERIC_OUTPUT, numeric_arrays)
check("the deterministic Jacobi archive contains exactly 560 arrays",
      archive_ok,
      f"arrays={len(numeric_arrays)}, sha={sha256(NUMERIC_OUTPUT)}")

schedule_records = []
for sector_index in range(7):
    left = records["even"][sector_index]
    right = records["odd"][sector_index]
    for family_name in ("K_family", "PQ_family"):
        left_singular = np.array([
            float(value) for value in left[family_name]["singular_values"]
        ])
        right_singular = np.array([
            float(value) for value in right[family_name]["singular_values"]
        ])
        distance = float(np.max(np.abs(left_singular - right_singular)))
        epsilon = (
            float(left[family_name]["epsilon"])
            + float(right[family_name]["epsilon"])
            + 10 * np.finfo(float).eps
            * max(1.0, float(left_singular[0]), float(right_singular[0]))
        )
        schedule_records.append({
            "sector_index": sector_index,
            "family": family_name,
            "distance": serialize_float(distance),
            "epsilon": serialize_float(epsilon),
            "label": schedule_label(distance, epsilon),
        })

schedule_counts = Counter(item["label"] for item in schedule_records)
schedule_complete = len(schedule_records) == 14
check("all fourteen target-free schedule comparisons receive labels",
      schedule_complete, str(dict(schedule_counts)))

controls_ok = provenance_ok and carrier_ok and archive_ok
if not controls_ok:
    outcome = "SHIFTED_JACOBI_CONTROL_FAILED"
elif not all_twists:
    outcome = "SHIFTED_JACOBI_TWIST_SINGULAR"
elif not all_principal_identities or not all_product_identities:
    outcome = "SHIFTED_JACOBI_VARIATIONAL_IDENTITY_FAILED"
elif schedule_counts["SCHEDULE_DEPENDENT"]:
    outcome = "SHIFTED_JACOBI_SCHEDULE_DEPENDENT"
elif schedule_counts["SCHEDULE_OPEN"]:
    outcome = "SHIFTED_JACOBI_SCHEDULE_OPEN"
else:
    outcome = "SHIFTED_JACOBI_CERTIFIED"

allowed = {
    "SHIFTED_JACOBI_CONTROL_FAILED",
    "SHIFTED_JACOBI_TWIST_SINGULAR",
    "SHIFTED_JACOBI_VARIATIONAL_IDENTITY_FAILED",
    "SHIFTED_JACOBI_SCHEDULE_DEPENDENT",
    "SHIFTED_JACOBI_SCHEDULE_OPEN",
    "SHIFTED_JACOBI_CERTIFIED",
}
check("the frozen hierarchy assigns the Jacobi outcome", outcome in allowed,
      outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "blind": True,
    "spatial_target_loaded": False,
    "outcome": outcome,
    "full_position_dimension": 720,
    "full_phase_dimension": 1440,
    "twist_determinants": twist_count,
    "numeric_archive": NUMERIC_OUTPUT.name,
    "numeric_archive_arrays": len(numeric_arrays),
    "numeric_archive_sha256": sha256(NUMERIC_OUTPUT),
    "schedule_counts": dict(schedule_counts),
    "schedule_comparisons": schedule_records,
    "parities": records,
    "classification": {
        "finite_shifted_three_slice_operator": "DERIVED COMPUTATIONAL",
        "general_variational_identity": "KNOWN",
        "wave_equation_interpretation": "OPEN",
        "dispersion_or_speed": "OPEN",
        "external_novelty": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
print(f"Numeric archive: {NUMERIC_OUTPUT}")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Blind centered mass--drift--stiffness census for the Regge Jacobi stencil.

Prior-art commit: 0b69556.
Preregistered protocol commit: 20cf463.
No spatial, continuum, dispersion, speed or particle target is loaded.
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
from scipy.optimize import linear_sum_assignment


HERE = Path(__file__).resolve().parent
SOURCE_JSON = HERE / "gravity_600cell_dust_three_slice_jacobi.json"
SOURCE_NPZ = HERE / "gravity_600cell_dust_three_slice_jacobi.npz"
SOURCE_VERIFIER = HERE / "verify_gravity_600cell_dust_three_slice_jacobi.py"
OUTPUT = HERE / "gravity_600cell_dust_centered_jacobi.json"
NUMERIC_OUTPUT = HERE / "gravity_600cell_dust_centered_jacobi.npz"

PRIOR_ART_COMMIT = "0b69556"
PROTOCOL_COMMIT = "20cf463"
EXPECTED_HASHES = {
    "json": "514e01937d621e82c240ea5cad621fb2bc699d09c4940b9be46fa1498152d90c",
    "npz": "63d95e79c11b25cada660f9a2422654eb92180263dad64e1cbf0ecc30b67d7f8",
    "verifier": "a875751eebb202dbb0b92780c3f48e7e275470360442f3ed1154e310cc36a884",
}
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
DIMENSIONS = (3, 2, 2, 2, 1, 1, 1)
POSITION_DIMS = tuple(30 * value for value in DIMENSIONS)
MATRIX_NAMES = ("M", "N", "V", "Gamma", "Omega")
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


def serialize_complex(value):
    return {
        "real": serialize_float(np.real(value)),
        "imaginary": serialize_float(np.imag(value)),
    }


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


def scale(matrix, factor):
    result = acb_mat(matrix.nrows(), matrix.ncols())
    scalar = acb(str(factor))
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            result[row, column] = scalar * matrix[row, column]
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
    return {
        "contains_zero_entrywise": all(
            matrix[row, column].contains(0)
            for row in range(matrix.nrows())
            for column in range(matrix.ncols())
        ),
        "midpoint_frobenius": serialize_float(la.norm(midpoint, "fro")),
        "radius_frobenius": serialize_float(la.norm(radii, "fro")),
        "maximum_midpoint_modulus": serialize_float(np.max(np.abs(midpoint))),
        "maximum_radius": serialize_float(np.max(radii)),
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


def family_uncertainty(midpoints, radii, hermitian=False):
    primary = midpoints[VARIANTS[0]]
    variation = max(
        la.norm(midpoints[variant] - primary, 2) for variant in VARIANTS
    )
    radius = max(la.norm(radii[variant], "fro") for variant in VARIANTS)
    norm = la.norm(primary, 2)
    floor = 10 * np.finfo(float).eps * max(1.0, float(norm))
    singular = la.svdvals(primary)
    result = {
        "primary": primary,
        "variation": float(variation),
        "radius": float(radius),
        "floor": float(floor),
        "epsilon": float(variation + radius + floor),
        "singular": singular,
    }
    if hermitian:
        result["eigenvalues"] = la.eigvalsh(primary)
    return result


def inertia_record(family, multiplicity):
    epsilon = family["epsilon"]
    values = family["eigenvalues"]
    positive = values > 100 * epsilon
    negative = values < -100 * epsilon
    zero = np.abs(values) < 10 * epsilon
    open_flags = ~(positive | negative | zero)
    if np.all(positive):
        label = "POSITIVE_DEFINITE"
    elif np.all(negative):
        label = "NEGATIVE_DEFINITE"
    elif np.any(positive) and np.any(negative):
        label = "INDEFINITE"
    else:
        label = "INERTIA_OPEN"
    counts = {
        "positive_minimal": int(np.sum(positive)),
        "negative_minimal": int(np.sum(negative)),
        "zero_minimal": int(np.sum(zero)),
        "open_minimal": int(np.sum(open_flags)),
    }
    counts.update({
        key.replace("minimal", "weighted"): int(multiplicity * value)
        for key, value in list(counts.items())
    })
    return {
        "label": label,
        "epsilon": serialize_float(epsilon),
        "minimum_eigenvalue": serialize_float(values[0]),
        "maximum_eigenvalue": serialize_float(values[-1]),
        "counts": counts,
    }


def optimal_spectral_distance(left, right):
    cost = np.abs(left[:, None] - right[None, :])
    rows, columns = linear_sum_assignment(cost)
    return float(np.max(cost[rows, columns]))


def operator_record(family, multiplicity, classify_reality=False):
    matrix = family["primary"]
    values, vectors = la.eig(matrix)
    try:
        eigen_condition = float(np.linalg.cond(vectors))
    except np.linalg.LinAlgError:
        eigen_condition = math.inf
    if math.isfinite(eigen_condition):
        eigen_epsilon = (
            eigen_condition * family["epsilon"]
            + 10 * np.finfo(float).eps
            * max(1.0, float(np.max(np.abs(values))))
        )
    else:
        eigen_epsilon = math.inf
    record = {
        "frobenius_norm": serialize_float(la.norm(matrix, "fro")),
        "minimum_singular": serialize_float(family["singular"][-1]),
        "maximum_singular": serialize_float(family["singular"][0]),
        "condition": serialize_float(
            family["singular"][0] / family["singular"][-1]
        ),
        "matrix_epsilon": serialize_float(family["epsilon"]),
        "eigenvector_condition": serialize_float(eigen_condition),
        "eigenvalue_epsilon": serialize_float(eigen_epsilon),
        "singular_values": [serialize_float(value) for value in family["singular"]],
        "eigenvalues": [serialize_complex(value) for value in values],
    }
    if classify_reality:
        if math.isfinite(eigen_epsilon):
            imaginary = np.abs(np.imag(values))
            real_consistent = imaginary < 10 * eigen_epsilon
            resolved_complex = imaginary > 100 * eigen_epsilon
            open_flags = ~(real_consistent | resolved_complex)
        else:
            real_consistent = np.zeros(len(values), dtype=bool)
            resolved_complex = np.zeros(len(values), dtype=bool)
            open_flags = np.ones(len(values), dtype=bool)
        counts = {
            "real_consistent_minimal": int(np.sum(real_consistent)),
            "resolved_complex_minimal": int(np.sum(resolved_complex)),
            "complex_open_minimal": int(np.sum(open_flags)),
        }
        counts.update({
            key.replace("minimal", "weighted"): int(multiplicity * value)
            for key, value in list(counts.items())
        })
        record["reality_counts"] = counts
    record["_eigen_array"] = values
    record["_eigen_epsilon_float"] = eigen_epsilon
    return record


def schedule_label(distance, epsilon):
    if not math.isfinite(epsilon):
        return "SCHEDULE_OPEN"
    if distance <= 10 * epsilon:
        return "SCHEDULE_ROBUST"
    if distance > 100 * epsilon:
        return "SCHEDULE_DEPENDENT"
    return "SCHEDULE_OPEN"


hashes = {
    "json": sha256(SOURCE_JSON),
    "npz": sha256(SOURCE_NPZ),
    "verifier": sha256(SOURCE_VERIFIER),
}
source_json = json.loads(SOURCE_JSON.read_text())
source_npz = np.load(SOURCE_NPZ)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and source_json["outcome"] == "THREE_SLICE_JACOBI_CERTIFIED"
    and source_json["numeric_archive_arrays"] == len(source_npz.files) == 560
    and source_json["numeric_archive_sha256"] == EXPECTED_HASHES["npz"]
    and source_json["spatial_target_loaded"] is False
)
check("the blind Jacobi source has exact frozen provenance",
      provenance_ok, str(hashes))

carrier_ok = all(
    tuple(item["irrep_dimension"] for item in source_json["parities"][parity])
    == DIMENSIONS
    and tuple(item["position_dimension"] for item in source_json["parities"][parity])
    == POSITION_DIMS
    for parity in ("even", "odd")
)
check("both schedules retain the complete seven-sector 720-position carrier",
      carrier_ok, f"positions={sum(30*d*d for d in DIMENSIONS)}")

print("=" * 78)
print("BLIND CENTERED JACOBI MASS--DRIFT--STIFFNESS CENSUS")
print("=" * 78)

records = {"even": [], "odd": []}
numeric_arrays = {}
all_mass_regular = True
all_centered_identities = True
mass_determinants = 0
aggregate_inertia = Counter()
aggregate_reality = Counter()

for parity in ("even", "odd"):
    for sector_index, (dimension, size) in enumerate(zip(DIMENSIONS, POSITION_DIMS)):
        print(f"[{parity}] sector {sector_index + 1}/7 n={size}", flush=True)
        mids = {name: {} for name in MATRIX_NAMES}
        rads = {name: {} for name in MATRIX_NAMES}
        h_mids = {}
        h_rads = {}
        variants_record = {}

        for variant in VARIANTS:
            base = f"{parity}_sector{sector_index}_{variant}"
            kminus = reenclose_binary_matrix(
                source_npz[f"{base}_Kminus_midpoint"],
                source_npz[f"{base}_Kminus_radii"],
            )
            kzero = reenclose_binary_matrix(
                source_npz[f"{base}_Kzero_midpoint"],
                source_npz[f"{base}_Kzero_radii"],
            )
            kplus = reenclose_binary_matrix(
                source_npz[f"{base}_Kplus_midpoint"],
                source_npz[f"{base}_Kplus_radii"],
            )
            m_matrix = scale(kminus + kplus, "0.5")
            n_matrix = scale(kplus - kminus, "0.5")
            v_matrix = kminus + kzero + kplus
            determinant = m_matrix.det()
            mass_determinants += 1
            mass_regular = not determinant.contains(0)
            all_mass_regular &= mass_regular
            gamma = m_matrix.solve(n_matrix)
            omega = m_matrix.solve(v_matrix)
            inverse_kminus = m_matrix.solve(kminus)
            inverse_kzero = m_matrix.solve(kzero)
            inverse_kplus = m_matrix.solve(kplus)
            one = identity(size)
            centered_residuals = {
                "Kminus": kminus - (m_matrix - n_matrix),
                "Kplus": kplus - (m_matrix + n_matrix),
                "Kzero": kzero - (v_matrix - scale(m_matrix, 2)),
                "normalized_Kminus": inverse_kminus - (one - gamma),
                "normalized_Kplus": inverse_kplus - (one + gamma),
                "normalized_Kzero": (
                    inverse_kzero - (omega - scale(one, 2))
                ),
            }
            residual_record = {
                name: residual_stats(value)
                for name, value in centered_residuals.items()
            }
            identity_ok = all(
                value["contains_zero_entrywise"]
                for value in residual_record.values()
            )
            all_centered_identities &= identity_ok
            matrices = {
                "M": m_matrix,
                "N": n_matrix,
                "V": v_matrix,
                "Gamma": gamma,
                "Omega": omega,
            }
            for name, matrix in matrices.items():
                midpoint, radii = acb_midpoint_and_radii(matrix)
                mids[name][variant] = midpoint
                rads[name][variant] = radii
                prefix = f"{parity}_sector{sector_index}_{variant}_{name}"
                numeric_arrays[f"{prefix}_midpoint"] = midpoint
                numeric_arrays[f"{prefix}_radii"] = radii
            h_matrix = scale(
                m_matrix + m_matrix.transpose().conjugate(), "0.5"
            )
            h_mids[variant], h_rads[variant] = acb_midpoint_and_radii(h_matrix)
            variants_record[variant] = {
                "mass_determinant_contains_zero": bool(determinant.contains(0)),
                "mass_determinant_ball": str(determinant),
                "centered_residuals": residual_record,
            }

        families = {
            name: family_uncertainty(mids[name], rads[name])
            for name in MATRIX_NAMES
        }
        h_family = family_uncertainty(h_mids, h_rads, hermitian=True)
        inertia = inertia_record(h_family, dimension)
        aggregate_inertia.update({
            key: value for key, value in inertia["counts"].items()
            if key.endswith("weighted")
        })
        gamma_record = operator_record(families["Gamma"], dimension)
        omega_record = operator_record(
            families["Omega"], dimension, classify_reality=True
        )
        aggregate_reality.update({
            key: value for key, value in omega_record["reality_counts"].items()
            if key.endswith("weighted")
        })
        primary = {name: mids[name][VARIANTS[0]] for name in ("M", "N", "V")}
        adjoint_ratios = {
            name: serialize_float(
                la.norm(matrix - matrix.conj().T, "fro")
                / max(1.0, la.norm(matrix, "fro"))
            )
            for name, matrix in primary.items()
        }
        records[parity].append({
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "position_dimension": size,
            "variants": variants_record,
            "mass_inertia": inertia,
            "adjoint_defect_ratios": adjoint_ratios,
            "matrix_frobenius_norms": {
                name: serialize_float(la.norm(mids[name][VARIANTS[0]], "fro"))
                for name in MATRIX_NAMES
            },
            "matrix_errors": {
                name: {
                    "variation": serialize_float(family["variation"]),
                    "radius": serialize_float(family["radius"]),
                    "floor": serialize_float(family["floor"]),
                    "epsilon": serialize_float(family["epsilon"]),
                }
                for name, family in families.items()
            },
            "Gamma": gamma_record,
            "Omega": omega_record,
        })

check("all 56 centered-mass determinant balls exclude zero",
      all_mass_regular, f"determinants={mass_determinants}")
check("all centered and normalized recurrence identities hold entrywise",
      all_centered_identities)

archive_ok = len(numeric_arrays) == 560
deterministic_npz(NUMERIC_OUTPUT, numeric_arrays)
check("the deterministic centered archive contains exactly 560 arrays",
      archive_ok,
      f"arrays={len(numeric_arrays)}, sha={sha256(NUMERIC_OUTPUT)}")

schedule_records = []
for sector_index in range(7):
    left = records["even"][sector_index]
    right = records["odd"][sector_index]
    for operator in ("Gamma", "Omega"):
        left_singular = np.array([
            float(value) for value in left[operator]["singular_values"]
        ])
        right_singular = np.array([
            float(value) for value in right[operator]["singular_values"]
        ])
        singular_distance = float(np.max(np.abs(left_singular - right_singular)))
        matrix_epsilon = (
            float(left[operator]["matrix_epsilon"])
            + float(right[operator]["matrix_epsilon"])
            + 10 * np.finfo(float).eps
            * max(1.0, float(left_singular[0]), float(right_singular[0]))
        )
        left_eigen = left[operator].pop("_eigen_array")
        right_eigen = right[operator].pop("_eigen_array")
        left_epsilon = left[operator].pop("_eigen_epsilon_float")
        right_epsilon = right[operator].pop("_eigen_epsilon_float")
        if operator == "Omega":
            eigen_distance = optimal_spectral_distance(left_eigen, right_eigen)
            eigen_epsilon = left_epsilon + right_epsilon
            eigen_label = schedule_label(eigen_distance, eigen_epsilon)
        else:
            eigen_distance = None
            eigen_epsilon = None
            eigen_label = "NOT_APPLICABLE"
        schedule_records.append({
            "sector_index": sector_index,
            "operator": operator,
            "singular_distance": serialize_float(singular_distance),
            "singular_epsilon": serialize_float(matrix_epsilon),
            "singular_label": schedule_label(singular_distance, matrix_epsilon),
            "eigen_distance": (
                serialize_float(eigen_distance) if eigen_distance is not None else None
            ),
            "eigen_epsilon": (
                serialize_float(eigen_epsilon) if eigen_epsilon is not None else None
            ),
            "eigen_label": eigen_label,
        })

primary_counts = Counter(item["singular_label"] for item in schedule_records)
secondary_counts = Counter(
    item["eigen_label"] for item in schedule_records
    if item["operator"] == "Omega"
)
schedule_complete = len(schedule_records) == 14
check("all fourteen blind centered schedule comparisons receive labels",
      schedule_complete,
      f"singular={dict(primary_counts)}, eigen={dict(secondary_counts)}")

controls_ok = provenance_ok and carrier_ok and archive_ok
resolved_dependent = any(
    item["singular_label"] == "SCHEDULE_DEPENDENT"
    or (item["operator"] == "Omega"
        and item["eigen_label"] == "SCHEDULE_DEPENDENT")
    for item in schedule_records
)
primary_open = any(
    item["singular_label"] == "SCHEDULE_OPEN" for item in schedule_records
)
if not controls_ok:
    outcome = "CENTERED_JACOBI_CONTROL_FAILED"
elif not all_mass_regular:
    outcome = "CENTERED_JACOBI_MASS_SINGULAR"
elif not all_centered_identities:
    outcome = "CENTERED_JACOBI_IDENTITY_FAILED"
elif resolved_dependent:
    outcome = "CENTERED_JACOBI_SCHEDULE_DEPENDENT"
elif primary_open:
    outcome = "CENTERED_JACOBI_SCHEDULE_OPEN"
else:
    outcome = "CENTERED_JACOBI_CERTIFIED"

allowed = {
    "CENTERED_JACOBI_CONTROL_FAILED",
    "CENTERED_JACOBI_MASS_SINGULAR",
    "CENTERED_JACOBI_IDENTITY_FAILED",
    "CENTERED_JACOBI_SCHEDULE_DEPENDENT",
    "CENTERED_JACOBI_SCHEDULE_OPEN",
    "CENTERED_JACOBI_CERTIFIED",
}
check("the frozen hierarchy assigns the centered-Jacobi outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "blind": True,
    "spatial_target_loaded": False,
    "outcome": outcome,
    "full_position_dimension": 720,
    "mass_determinants": mass_determinants,
    "numeric_archive": NUMERIC_OUTPUT.name,
    "numeric_archive_arrays": len(numeric_arrays),
    "numeric_archive_sha256": sha256(NUMERIC_OUTPUT),
    "aggregate_mass_inertia": dict(aggregate_inertia),
    "aggregate_omega_reality": dict(aggregate_reality),
    "schedule_label_counts": {
        "singular": dict(primary_counts),
        "eigenvalue": dict(secondary_counts),
    },
    "schedule_comparisons": schedule_records,
    "parities": records,
    "classification": {
        "centered_finite_operator": "DERIVED COMPUTATIONAL",
        "mass_hermitian_inertia": "STRUCTURAL DIAGNOSTIC",
        "generalized_stiffness_spectrum": "STRUCTURAL",
        "physical_wave_interpretation": "OPEN",
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
print(f"Mass inertia: {dict(aggregate_inertia)}")
print(f"Omega reality: {dict(aggregate_reality)}")
print(f"Artifact: {OUTPUT}")
print(f"Numeric archive: {NUMERIC_OUTPUT}")
if passed != tests:
    raise SystemExit(1)

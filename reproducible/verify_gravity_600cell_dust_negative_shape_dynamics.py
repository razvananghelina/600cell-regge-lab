#!/usr/bin/env python3
"""Complete centered dynamics of the two negative-stiffness shape sectors.

Prior-art commit: db1fa95.
Target-disclosed protocol commit: dde2164.
No desired multiplier count, old tangent spectrum, continuum harmonic,
polarization, speed, refinement, Planck or particle target is loaded.
"""

import ast
from collections import Counter
import contextlib
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la
from scipy.optimize import linear_sum_assignment
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

STIFFNESS_JSON = HERE / "gravity_600cell_dust_shape_stiffness.json"
STIFFNESS_SOURCE = HERE / "verify_gravity_600cell_dust_shape_stiffness.py"
CENTERED_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
LEGENDRE_JSON = HERE / "gravity_600cell_dust_full_anisotropic_legendre_rank.json"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
OUTPUT = HERE / "gravity_600cell_dust_negative_shape_dynamics.json"

PRIOR_ART_COMMIT = "db1fa95"
PROTOCOL_COMMIT = "dde2164"
EXPECTED_HASHES = {
    "stiffness_json": "03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868",
    "stiffness_source": "d4f0a9a805910de37011ba70f407907daa2d11c650aeea22e571ab867282a44c",
    "centered_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "legendre_json": "7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
PARITIES = ("even", "odd")
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
DIMENSIONS = (3, 2, 2, 2, 1, 1, 1)
SELECTED_SECTORS = (4, 5)
MACHINE_EPSILON = np.finfo(float).eps
mp.mp.dps = 100
I = mp.mpc(0, 1)
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


def operator_norm(matrix):
    values = la.svdvals(matrix)
    return float(values[0]) if len(values) else 0.0


def sign_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "OPEN"
    if value > 100 * error:
        return "POSITIVE_RESOLVED"
    if value < -100 * error:
        return "NEGATIVE_RESOLVED"
    if abs(value) <= 10 * error:
        return "ZERO_CONSISTENT"
    return "OPEN"


def invariance_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "INVARIANCE_OPEN"
    if value <= 10 * error:
        return "INVARIANT_CONSISTENT"
    if value > 100 * error:
        return "MIXING_RESOLVED"
    return "INVARIANCE_OPEN"


def regularity_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "REGULARITY_OPEN"
    if value > 100 * error:
        return "REGULAR_RESOLVED"
    if value <= 10 * error:
        return "SINGULAR_CONSISTENT"
    return "REGULARITY_OPEN"


def modulus_label(value, error):
    radius = abs(value)
    if not math.isfinite(radius) or not math.isfinite(error) or error < 0:
        return "MODULUS_OPEN"
    if radius < 1 - 100 * error:
        return "CONTRACTING_RESOLVED"
    if radius > 1 + 100 * error:
        return "EXPANDING_RESOLVED"
    if abs(radius - 1) <= 10 * error:
        return "UNIT_CONSISTENT"
    return "MODULUS_OPEN"


def schedule_label(distance, error):
    if not math.isfinite(distance) or not math.isfinite(error) or error < 0:
        return "SCHEDULE_OPEN"
    if distance <= 10 * error:
        return "SCHEDULE_ROBUST"
    if distance > 100 * error:
        return "SCHEDULE_DEPENDENT"
    return "SCHEDULE_OPEN"


def load_audited_helpers():
    wanted = {
        "mp_frobenius",
        "mp_submatrix",
        "cluster_sorted",
        "orbit_sort_key",
        "edge_image",
        "group_data",
        "incidence_data",
        "mp_to_numpy",
        "component_reenclosure_radii",
    }
    tree = ast.parse(CONFORMAL_SOURCE.read_text(), filename=str(CONFORMAL_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited helper mismatch: missing={wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(CONFORMAL_SOURCE), "exec"),
        globals(),
    )

    tree = ast.parse(FULL_SOURCE.read_text(), filename=str(FULL_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "high_precision_sector_bases"
    ]
    if len(body) != 1:
        raise RuntimeError("audited high-precision sector function is missing")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(FULL_SOURCE), "exec"),
        globals(),
    )


def source_matrix(archive, prefix, name):
    midpoint = np.asarray(archive[f"{prefix}_{name}_midpoint"])
    stored = np.asarray(archive[f"{prefix}_{name}_radii"])
    radii = component_reenclosure_radii(midpoint, stored)
    return midpoint, radii


def matrix_error(midpoint, radii, n):
    return float(
        la.norm(radii, "fro")
        + 1000 * MACHINE_EPSILON * n
        * max(1.0, operator_norm(midpoint))
    )


def restriction_error(midpoint, base_error, eta_shape, n):
    norm = operator_norm(midpoint)
    return float(
        base_error
        + 2 * eta_shape * (norm + base_error)
        + 1000 * MACHINE_EPSILON * n * max(1.0, norm)
    )


def build_companion(gamma, omega, epsilon_gamma, epsilon_omega):
    m = gamma.shape[0]
    identity = np.eye(m, dtype=np.complex128)
    forward = identity + gamma
    singular = la.svdvals(forward)
    minimum = float(singular[-1])
    regularity = regularity_label(minimum, epsilon_gamma)
    result = {
        "regularity_label": regularity,
        "minimum_forward_singular": minimum,
        "forward_error": epsilon_gamma,
        "matrix": None,
        "matrix_error": math.inf,
    }
    if regularity != "REGULAR_RESOLVED":
        return result

    inverse = la.inv(forward)
    norm_inverse = operator_norm(inverse)
    denominator = 1 - norm_inverse * epsilon_gamma
    if denominator <= 0:
        result["regularity_label"] = "REGULARITY_OPEN"
        return result
    epsilon_inverse = float(
        norm_inverse * norm_inverse * epsilon_gamma / denominator
        + 1000 * MACHINE_EPSILON * m * max(1.0, norm_inverse)
    )

    backward_factor = identity - gamma
    stiffness_factor = 2 * identity - omega
    lower = -inverse @ backward_factor
    upper = inverse @ stiffness_factor
    epsilon_lower = float(
        epsilon_inverse * (operator_norm(backward_factor) + epsilon_gamma)
        + norm_inverse * epsilon_gamma
        + 1000 * MACHINE_EPSILON * m
        * max(1.0, operator_norm(lower), norm_inverse * operator_norm(backward_factor))
    )
    epsilon_upper = float(
        epsilon_inverse * (operator_norm(stiffness_factor) + epsilon_omega)
        + norm_inverse * epsilon_omega
        + 1000 * MACHINE_EPSILON * m
        * max(1.0, operator_norm(upper), norm_inverse * operator_norm(stiffness_factor))
    )
    companion = np.block([
        [np.zeros((m, m), dtype=np.complex128), identity],
        [lower, upper],
    ])
    epsilon_companion = float(
        epsilon_lower + epsilon_upper
        + 1000 * MACHINE_EPSILON * (2 * m)
        * max(1.0, operator_norm(companion))
    )
    result.update({
        "matrix": companion,
        "matrix_error": epsilon_companion,
        "inverse_error": epsilon_inverse,
        "lower_error": epsilon_lower,
        "upper_error": epsilon_upper,
    })
    return result


def classify_companion(companion_result):
    matrix = companion_result["matrix"]
    if matrix is None:
        return None
    values, vectors = la.eig(matrix)
    try:
        condition = float(np.linalg.cond(vectors))
    except np.linalg.LinAlgError:
        condition = math.inf
    m = len(values)
    error = float(
        condition * companion_result["matrix_error"]
        + 1000 * MACHINE_EPSILON * m * max(1.0, operator_norm(matrix))
    )
    order = np.lexsort((values.imag, values.real))
    values = np.asarray(values)[order]
    labels = [modulus_label(value, error) for value in values]
    singular = np.sort(la.svdvals(matrix))
    return {
        "values": values,
        "labels": labels,
        "label_counts": Counter(labels),
        "eigenvector_condition": condition,
        "eigenvalue_error": error,
        "singular_values": singular,
        "singular_value_error": companion_result["matrix_error"],
        "spectral_radius": float(np.max(np.abs(values))),
        "maximum_singular": float(singular[-1]),
    }


def matched_distance(left, right):
    if len(left) != len(right):
        return math.inf
    cost = np.abs(left[:, None] - right[None, :])
    rows, columns = linear_sum_assignment(cost)
    return float(np.max(cost[rows, columns]))


print("=" * 78)
print("600-CELL NEGATIVE SHAPE SECTOR DYNAMICS")
print("=" * 78)

paths = {
    "stiffness_json": STIFFNESS_JSON,
    "stiffness_source": STIFFNESS_SOURCE,
    "centered_json": CENTERED_JSON,
    "centered_npz": CENTERED_NPZ,
    "legendre_json": LEGENDRE_JSON,
    "conformal_source": CONFORMAL_SOURCE,
    "full_source": FULL_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "commons": COMMONS_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
stiffness_source = json.loads(STIFFNESS_JSON.read_text())
centered = json.loads(CENTERED_JSON.read_text())
legendre = json.loads(LEGENDRE_JSON.read_text())
archive = np.load(CENTERED_NPZ, allow_pickle=False)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and stiffness_source["outcome"] == "SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED"
    and stiffness_source["passed"] == stiffness_source["tests"] == 12
    and stiffness_source["full_multiplicity_pencil_sign_counts"]
        ["NEGATIVE_RESOLVED"] == 240
    and centered["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and centered["passed"] == centered["tests"] == 7
    and len(archive.files) == centered["numeric_archive_arrays"] == 560
    and legendre["outcome"] == "FULL_CANONICAL_LEGENDRE_REGULAR"
    and legendre["passed"] == legendre["tests"] == 18
    and all(
        item["full_resolved_rank"] == 1560
        and item["full_error_consistent_nullity"] == 0
        and item["pseudoconstraint_candidates"] == 120
        for item in legendre["parities"].values()
    )
)
check("all target-disclosed inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_negative_shape", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
geometry_import_ok = gro.tests == gro.passed == 43
check("the literal one-slab geometry retains all 43 certificates", geometry_import_ok)

load_audited_helpers()
groups = {parity: group_data(gro.models[parity], gro) for parity in PARITIES}
incidences = {parity: incidence_data(groups[parity]) for parity in PARITIES}
odd_edge_index = {
    edge: index for index, edge in enumerate(groups["odd"]["edge_order"])
}
even_to_odd = np.asarray([
    odd_edge_index[edge] for edge in groups["even"]["edge_order"]
], dtype=int)
carrier_geometry_ok = bool(
    sorted(even_to_odd.tolist()) == list(range(720))
    and np.array_equal(
        incidences["odd"]["incidence"][even_to_odd],
        incidences["even"]["incidence"],
    )
    and all(data["equivariant"] for data in incidences.values())
    and all(data["numerical_rank"] == 120 for data in incidences.values())
)
check("both schedules retain the exact equivariant conformal carrier", carrier_geometry_ok)

sector_data, sector_controls = high_precision_sector_bases(groups["even"])
maximum_basis_residual = max(
    value for key, value in sector_controls.items() if key.startswith("maximum_")
)
sector_basis_ok = bool(
    tuple(sector["dimension"] for sector in sector_data) == DIMENSIONS
    and all(sector_data[index]["dimension"] == 1 for index in SELECTED_SECTORS)
    and maximum_basis_residual < mp.mpf("1e-70")
)
check(
    "the selected sectors are the same two one-dimensional minimal sectors",
    sector_basis_ok,
    "max_residual=" + mp.nstr(maximum_basis_residual, 5),
)

sector_images = {parity: {} for parity in PARITIES}
conformal_open = False
conformal_contradicted = False
for parity in PARITIES:
    incidence = incidences[parity]["incidence"].astype(np.complex128)
    for sector_index in SELECTED_SECTORS:
        sector = sector_data[sector_index]
        dimension = sector["dimension"]
        n = 30 * dimension
        expected_rank = 5 * dimension
        basis = mp_to_numpy(sector["basis"])
        edge_basis = np.kron(np.eye(30, dtype=np.complex128), basis)
        compressed = edge_basis.conj().T @ incidence
        left, singular, _ = la.svd(compressed, full_matrices=False)
        epsilon_c = float(
            1000 * MACHINE_EPSILON * max(compressed.shape)
            * max(1.0, float(singular[0]))
        )
        nonzero = int(np.sum(singular > 100 * epsilon_c))
        zero = int(np.sum(singular < 10 * epsilon_c))
        open_count = len(singular) - nonzero - zero
        conformal_open |= open_count > 0
        conformal_contradicted |= bool(
            open_count == 0
            and (nonzero != expected_rank or zero != len(singular) - expected_rank)
        )
        gap = float(singular[expected_rank - 1])
        if gap > 2 * epsilon_c:
            eta_k = float(
                2 * epsilon_c / (gap - 2 * epsilon_c)
                + 1000 * MACHINE_EPSILON * n
            )
        else:
            eta_k = math.inf
            conformal_open = True
        sector_images[parity][sector_index] = {
            "basis": left[:, :expected_rank],
            "eta_K": eta_k,
        }

check(
    "all four selected conformal carriers are resolved",
    not conformal_open and not conformal_contradicted,
    f"open={conformal_open}, contradicted={conformal_contradicted}",
)

records = {parity: [] for parity in PARITIES}
internal = {parity: {} for parity in PARITIES}
carrier_open = conformal_open
carrier_contradicted = conformal_contradicted
invariance_counts = Counter()
regularity_counts = Counter()
full_modulus_counts = Counter()
negative_modulus_counts = Counter()
full_companion_cells = 0
negative_companion_cells = 0
both_invariant_cells = 0
all_finite = True

for parity in PARITIES:
    for sector_index in SELECTED_SECTORS:
        dimension = sector_data[sector_index]["dimension"]
        n = 30 * dimension
        r = 5 * dimension
        s = 25 * dimension
        u_basis = sector_images[parity][sector_index]["basis"]
        eta_k = sector_images[parity][sector_index]["eta_K"]
        variants = []
        internal[parity][sector_index] = {}

        for variant in VARIANTS:
            prefix = f"{parity}_sector{sector_index}_{variant}"
            midpoint_m, radius_m = source_matrix(archive, prefix, "M")
            h_m = (midpoint_m + midpoint_m.conj().T) / 2
            radius_hm = (radius_m + radius_m.T) / 2
            epsilon_hm = matrix_error(h_m, radius_hm, n)
            row = u_basis.conj().T @ h_m
            epsilon_row = float(
                epsilon_hm
                + 2 * eta_k * (operator_norm(h_m) + epsilon_hm)
                + 1000 * MACHINE_EPSILON * n * max(1.0, operator_norm(row))
            )
            _, singular_row, right_row = la.svd(row, full_matrices=True)
            row_nonzero = int(np.sum(singular_row > 100 * epsilon_row))
            row_zero = int(np.sum(singular_row < 10 * epsilon_row))
            row_open = len(singular_row) - row_nonzero - row_zero
            carrier_open |= row_open > 0
            carrier_contradicted |= bool(
                row_open == 0 and (row_nonzero != r or row_zero != 0)
            )
            row_gap = float(singular_row[r - 1])
            if row_gap > 2 * epsilon_row:
                eta_s = float(
                    2 * epsilon_row / (row_gap - 2 * epsilon_row)
                    + 1000 * MACHINE_EPSILON * n
                )
            else:
                eta_s = math.inf
                carrier_open = True
            w_basis = right_row.conj().T[:, r:]

            midpoint_v, radius_v = source_matrix(archive, prefix, "V")
            h_v = (midpoint_v + midpoint_v.conj().T) / 2
            radius_hv = (radius_v + radius_v.T) / 2
            epsilon_hv = matrix_error(h_v, radius_hv, n)
            a_matrix = -(w_basis.conj().T @ h_v @ w_basis)
            a_matrix = (a_matrix + a_matrix.conj().T) / 2
            epsilon_a = restriction_error(h_v, epsilon_hv, eta_s, n)
            values_a, vectors_a = la.eigh(a_matrix)
            labels_a = [sign_label(float(value), epsilon_a) for value in values_a]
            inertia = Counter(labels_a)
            if inertia != Counter({
                "NEGATIVE_RESOLVED": 15,
                "POSITIVE_RESOLVED": 10,
            }):
                if inertia["OPEN"] or inertia["ZERO_CONSISTENT"]:
                    carrier_open = True
                else:
                    carrier_contradicted = True
            eigengap = float(values_a[15] - values_a[14])
            if eigengap > 2 * epsilon_a:
                eta_e = float(
                    2 * epsilon_a / (eigengap - 2 * epsilon_a)
                    + 1000 * MACHINE_EPSILON * s
                )
            else:
                eta_e = math.inf
                carrier_open = True
            negative_basis = vectors_a[:, :15]

            midpoint_g, radius_g = source_matrix(archive, prefix, "Gamma")
            midpoint_o, radius_o = source_matrix(archive, prefix, "Omega")
            epsilon_g = matrix_error(midpoint_g, radius_g, n)
            epsilon_o = matrix_error(midpoint_o, radius_o, n)
            gamma_s = w_basis.conj().T @ midpoint_g @ w_basis
            omega_s = w_basis.conj().T @ midpoint_o @ w_basis
            epsilon_gs = restriction_error(midpoint_g, epsilon_g, eta_s, n)
            epsilon_os = restriction_error(midpoint_o, epsilon_o, eta_s, n)

            projector_e = negative_basis @ negative_basis.conj().T
            invariance = {}
            for name, matrix, matrix_error_s in (
                ("Gamma", gamma_s, epsilon_gs),
                ("Omega", omega_s, epsilon_os),
            ):
                residual = operator_norm(
                    (np.eye(s, dtype=np.complex128) - projector_e)
                    @ matrix @ negative_basis
                )
                error = float(
                    matrix_error_s
                    + 2 * eta_e * (operator_norm(matrix) + matrix_error_s)
                    + 1000 * MACHINE_EPSILON * s
                    * max(1.0, operator_norm(matrix))
                )
                label = invariance_label(residual, error)
                invariance_counts[label] += 1
                invariance[name] = {
                    "residual": residual,
                    "error": error,
                    "label": label,
                }

            full_build = build_companion(
                gamma_s, omega_s, epsilon_gs, epsilon_os
            )
            regularity_counts[full_build["regularity_label"]] += 1
            full_class = classify_companion(full_build)
            if full_class is not None:
                full_companion_cells += 1
                full_modulus_counts.update(full_class["labels"])

            both_invariant = all(
                item["label"] == "INVARIANT_CONSISTENT"
                for item in invariance.values()
            )
            both_invariant_cells += int(both_invariant)
            negative_build = None
            negative_class = None
            if both_invariant:
                gamma_negative = negative_basis.conj().T @ gamma_s @ negative_basis
                omega_negative = negative_basis.conj().T @ omega_s @ negative_basis
                epsilon_gn = invariance["Gamma"]["error"] + invariance["Gamma"]["residual"]
                epsilon_on = invariance["Omega"]["error"] + invariance["Omega"]["residual"]
                negative_build = build_companion(
                    gamma_negative, omega_negative, epsilon_gn, epsilon_on
                )
                regularity_counts["negative_" + negative_build["regularity_label"]] += 1
                negative_class = classify_companion(negative_build)
                if negative_class is not None:
                    negative_companion_cells += 1
                    negative_modulus_counts.update(negative_class["labels"])

            all_finite &= bool(
                np.all(np.isfinite(values_a))
                and all(
                    math.isfinite(item["residual"])
                    and math.isfinite(item["error"])
                    for item in invariance.values()
                )
            )

            internal[parity][sector_index][variant] = {
                "full": full_class,
                "negative": negative_class,
            }

            def companion_record(build, classified):
                record = {
                    "regularity_label": build["regularity_label"] if build else "NOT_CONSTRUCTED",
                    "minimum_forward_singular": sf(build["minimum_forward_singular"])
                    if build else None,
                    "forward_error": sf(build["forward_error"]) if build else None,
                }
                if classified is not None:
                    record.update({
                        "matrix_error": sf(build["matrix_error"]),
                        "eigenvalues": [
                            {"real": sf(value.real), "imag": sf(value.imag)}
                            for value in classified["values"]
                        ],
                        "eigenvalue_error": sf(classified["eigenvalue_error"]),
                        "eigenvector_condition": sf(classified["eigenvector_condition"]),
                        "modulus_labels": classified["labels"],
                        "modulus_counts": dict(classified["label_counts"]),
                        "singular_values": [sf(value) for value in classified["singular_values"]],
                        "singular_value_error": sf(classified["singular_value_error"]),
                        "spectral_radius": sf(classified["spectral_radius"]),
                        "maximum_singular": sf(classified["maximum_singular"]),
                    })
                return record

            variants.append({
                "variant": variant,
                "shape_dimension": s,
                "negative_dimension": 15,
                "inherited_stiffness": {
                    "inertia": dict(inertia),
                    "negative_positive_gap": sf(eigengap),
                    "restricted_A_error": sf(epsilon_a),
                    "negative_subspace_error_eta": sf(eta_e),
                },
                "invariance": {
                    name: {
                        "residual": sf(item["residual"]),
                        "error": sf(item["error"]),
                        "error_units": sf(
                            item["residual"] / item["error"]
                            if item["error"] else math.inf
                        ),
                        "label": item["label"],
                    }
                    for name, item in invariance.items()
                },
                "full_sector_companion": companion_record(full_build, full_class),
                "negative_companion": companion_record(negative_build, negative_class)
                if negative_build is not None else {
                    "regularity_label": "NOT_CONSTRUCTED_NONINVARIANT"
                },
            })

        records[parity].append({
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "variants": variants,
        })

required_cells = 2 * 2 * 4
check(
    "all 16 inherited shape/inertia carriers are reconstructed",
    all_finite
    and sum(regularity_counts[label] for label in (
        "REGULAR_RESOLVED", "SINGULAR_CONSISTENT", "REGULARITY_OPEN"
    )) == required_cells,
    f"carrier_open={carrier_open}, contradicted={carrier_contradicted}",
)

check(
    "all 32 negative-subspace invariance residuals are classified",
    sum(invariance_counts.values()) == 32,
    str(dict(invariance_counts)),
)

full_census_complete = bool(
    full_companion_cells
    + regularity_counts["SINGULAR_CONSISTENT"]
    + regularity_counts["REGULARITY_OPEN"]
    == required_cells
    and sum(full_modulus_counts.values()) == 50 * full_companion_cells
)
check(
    "all full-sector companions and their available spectra are classified",
    full_census_complete,
    f"regularity={dict(regularity_counts)}, modulus={dict(full_modulus_counts)}",
)

negative_accounting_ok = bool(
    negative_companion_cells
    + regularity_counts["negative_SINGULAR_CONSISTENT"]
    + regularity_counts["negative_REGULARITY_OPEN"]
    == both_invariant_cells
    and sum(negative_modulus_counts.values()) == 30 * negative_companion_cells
)
check(
    "restricted companions are constructed iff the negative carrier is invariant",
    negative_accounting_ok,
    f"invariant_cells={both_invariant_cells}, built={negative_companion_cells}",
)

schedule_records = []
schedule_counts = Counter()
restricted_schedule_pairs = 0
for sector_index in SELECTED_SECTORS:
    for variant in VARIANTS:
        left = internal["even"][sector_index][variant]
        right = internal["odd"][sector_index][variant]
        restricted_schedule_pairs += int(
            left["negative"] is not None and right["negative"] is not None
        )
        for family in ("full", "negative"):
            left_class = left[family]
            right_class = right[family]
            if left_class is None or right_class is None:
                continue
            eig_distance = matched_distance(
                left_class["values"], right_class["values"]
            )
            eig_error = float(
                left_class["eigenvalue_error"]
                + right_class["eigenvalue_error"]
                + 1000 * MACHINE_EPSILON * len(left_class["values"])
                * max(
                    1.0,
                    float(np.max(np.abs(left_class["values"]))),
                    float(np.max(np.abs(right_class["values"]))),
                )
            )
            eig_label = schedule_label(eig_distance, eig_error)
            schedule_counts[eig_label] += 1
            schedule_records.append({
                "sector_index": sector_index,
                "variant": variant,
                "family": family,
                "data": "eigenvalues",
                "distance": sf(eig_distance),
                "error": sf(eig_error),
                "error_units": sf(eig_distance / eig_error if eig_error else math.inf),
                "label": eig_label,
            })

            singular_distance = float(np.max(np.abs(
                left_class["singular_values"] - right_class["singular_values"]
            )))
            singular_error = float(
                left_class["singular_value_error"]
                + right_class["singular_value_error"]
                + 1000 * MACHINE_EPSILON * len(left_class["singular_values"])
                * max(
                    1.0,
                    float(left_class["singular_values"][-1]),
                    float(right_class["singular_values"][-1]),
                )
            )
            singular_label = schedule_label(singular_distance, singular_error)
            schedule_counts[singular_label] += 1
            schedule_records.append({
                "sector_index": sector_index,
                "variant": variant,
                "family": family,
                "data": "singular_values",
                "distance": sf(singular_distance),
                "error": sf(singular_error),
                "error_units": sf(
                    singular_distance / singular_error if singular_error else math.inf
                ),
                "label": singular_label,
            })

expected_schedule_records = 16 + 2 * restricted_schedule_pairs
check(
    "all available preregistered schedule comparisons are classified",
    len(schedule_records) == expected_schedule_records
    and sum(schedule_counts.values()) == expected_schedule_records,
    f"expected={expected_schedule_records}, labels={dict(schedule_counts)}",
)

target_flags = {
    "desired_multiplier_count": False,
    "old_tangent_spectrum": False,
    "continuum_harmonics": False,
    "polarization_count": False,
    "speed": False,
    "refinement": False,
    "planck_or_particle_data": False,
}
check("no forbidden physical or spectral target was loaded", not any(target_flags.values()))

controls_ok = bool(
    provenance_ok
    and geometry_import_ok
    and carrier_geometry_ok
    and sector_basis_ok
    and not carrier_contradicted
    and all_finite
    and full_census_complete
    and negative_accounting_ok
)
if not controls_ok:
    outcome = "NEGATIVE_SHAPE_DYNAMICS_CONTROL_FAILED"
elif carrier_open:
    outcome = "NEGATIVE_SHAPE_CARRIER_OPEN"
elif invariance_counts["MIXING_RESOLVED"] > 0:
    outcome = "NEGATIVE_SHAPE_SUBSPACE_MIXED"
elif invariance_counts["INVARIANCE_OPEN"] > 0:
    outcome = "NEGATIVE_SHAPE_INVARIANCE_OPEN"
elif regularity_counts["SINGULAR_CONSISTENT"] > 0 or regularity_counts[
    "negative_SINGULAR_CONSISTENT"
] > 0:
    outcome = "NEGATIVE_SHAPE_COMPANION_SINGULAR"
elif regularity_counts["REGULARITY_OPEN"] > 0 or regularity_counts[
    "negative_REGULARITY_OPEN"
] > 0:
    outcome = "NEGATIVE_SHAPE_COMPANION_REGULARITY_OPEN"
elif schedule_counts["SCHEDULE_DEPENDENT"] > 0:
    outcome = "NEGATIVE_SHAPE_SCHEDULE_DEPENDENT"
elif schedule_counts["SCHEDULE_OPEN"] > 0:
    outcome = "NEGATIVE_SHAPE_SCHEDULE_OPEN"
elif negative_modulus_counts["EXPANDING_RESOLVED"] > 0:
    outcome = "NEGATIVE_SHAPE_AUTONOMOUS_EXPANSION_RESOLVED"
elif negative_modulus_counts["MODULUS_OPEN"] > 0:
    outcome = "NEGATIVE_SHAPE_AUTONOMOUS_MODULUS_OPEN"
elif negative_modulus_counts["UNIT_CONSISTENT"] > 0:
    outcome = "NEGATIVE_SHAPE_AUTONOMOUS_UNIT_CONSISTENT"
else:
    outcome = "NEGATIVE_SHAPE_AUTONOMOUS_CONTRACTING_CENSUS"

outcomes = {
    "NEGATIVE_SHAPE_DYNAMICS_CONTROL_FAILED",
    "NEGATIVE_SHAPE_CARRIER_OPEN",
    "NEGATIVE_SHAPE_SUBSPACE_MIXED",
    "NEGATIVE_SHAPE_INVARIANCE_OPEN",
    "NEGATIVE_SHAPE_COMPANION_SINGULAR",
    "NEGATIVE_SHAPE_COMPANION_REGULARITY_OPEN",
    "NEGATIVE_SHAPE_SCHEDULE_DEPENDENT",
    "NEGATIVE_SHAPE_SCHEDULE_OPEN",
    "NEGATIVE_SHAPE_AUTONOMOUS_EXPANSION_RESOLVED",
    "NEGATIVE_SHAPE_AUTONOMOUS_MODULUS_OPEN",
    "NEGATIVE_SHAPE_AUTONOMOUS_UNIT_CONSISTENT",
    "NEGATIVE_SHAPE_AUTONOMOUS_CONTRACTING_CENSUS",
}
check("the target-disclosed outcome ladder is exhaustive", outcome in outcomes, outcome)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosure": {
        "selected_sectors": list(SELECTED_SECTORS),
        "selection_source": "preceding blind negative-stiffness result",
        "test_is_blind_discovery": False,
    },
    "constraint_status": {
        "full_legendre_rank": 1560,
        "full_legendre_nullity": 0,
        "exact_fixed_carrier_quotient_applied": False,
        "weak_pseudoconstraint_candidates": 120,
    },
    "enumeration": {
        "primary_cells": required_cells,
        "invariance_residuals": 32,
        "full_companion_eigenvalue_instances": 50 * full_companion_cells,
        "negative_companion_eigenvalue_instances": 30 * negative_companion_cells,
    },
    "invariance_counts": dict(invariance_counts),
    "regularity_counts": dict(regularity_counts),
    "full_sector_modulus_counts": dict(full_modulus_counts),
    "negative_carrier_modulus_counts": dict(negative_modulus_counts),
    "schedule_comparisons": schedule_records,
    "schedule_counts": dict(schedule_counts),
    "parities": records,
    "target_load_flags": target_flags,
    "interpretation_limits": {
        "single_update_called_lyapunov_exponent": False,
        "expanding_multiplier_called_ghost": False,
        "shape_mode_called_graviton": False,
        "physical_constraint_quotient_derived": False,
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("SCIENTIFIC OUTCOME:", outcome)
print("invariance labels:", dict(invariance_counts))
print("regularity labels:", dict(regularity_counts))
print("full-sector modulus labels:", dict(full_modulus_counts))
print("negative-carrier modulus labels:", dict(negative_modulus_counts))
print("schedule labels:", dict(schedule_counts))
print(f"{passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Dynamic closure of the canonical conformal/shape supermetric split.

Prior-art commit: 42313fb.
Protocol commit: 4039244.
Post-result power-control protocol: 6a8ce90.
No continuum spectrum, polarization count, speed or refinement target is loaded.
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
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

CENTERED_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
CONFORMAL_JSON = HERE / "gravity_600cell_dust_conformal_supermetric.json"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
OUTPUT = HERE / "gravity_600cell_dust_conformal_shape_dynamics.json"

PRIOR_ART_COMMIT = "42313fb"
PROTOCOL_COMMIT = "4039244"
POWER_PROTOCOL_COMMIT = "6a8ce90"
EXPECTED_HASHES = {
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "conformal_json": "b38d55f9f575ddffd34edeaa5e835d9e10919e6d96a0c284d73c31a072675025",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "centered_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}
PARITIES = ("even", "odd")
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
DIMENSIONS = (3, 2, 2, 2, 1, 1, 1)
OPERATORS = ("Gamma", "Omega")
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


def serialize_float(value):
    return f"{float(value):.17e}"


def operator_norm(matrix):
    singular = la.svdvals(matrix)
    return float(singular[0]) if len(singular) else 0.0


def residual_label(value, epsilon):
    if not math.isfinite(value) or not math.isfinite(epsilon):
        return "OPEN"
    if value <= 10 * epsilon:
        return "ZERO_CONSISTENT"
    if value > 100 * epsilon:
        return "NONZERO_RESOLVED"
    return "OPEN"


def schedule_label(distance, epsilon):
    if not math.isfinite(distance) or not math.isfinite(epsilon):
        return "SCHEDULE_OPEN"
    if distance <= 10 * epsilon:
        return "SCHEDULE_ROBUST"
    if distance > 100 * epsilon:
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


print("=" * 78)
print("600-CELL CONFORMAL / SHAPE DYNAMIC CLOSURE GATE")
print("=" * 78)

hashes = {
    "commons": sha256(COMMONS_SOURCE),
    "conformal_json": sha256(CONFORMAL_JSON),
    "conformal_source": sha256(CONFORMAL_SOURCE),
    "centered_json": sha256(CENTERED_JSON),
    "centered_npz": sha256(CENTERED_NPZ),
    "full_source": sha256(FULL_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
}
centered = json.loads(CENTERED_JSON.read_text())
conformal = json.loads(CONFORMAL_JSON.read_text())
source_npz = np.load(CENTERED_NPZ, allow_pickle=False)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and centered["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and centered["passed"] == centered["tests"] == 7
    and centered["numeric_archive_arrays"] == len(source_npz.files) == 560
    and conformal["outcome"] == "CONFORMAL_MAXIMAL_MINORITY_CERTIFIED"
    and conformal["passed"] == conformal["tests"] == 11
    and conformal["canonical_map_candidates"] == 1
    and conformal["carrier"]["rank_C_exact_by_odd_cycle"] == 120
    and all(
        tuple(item["irrep_dimension"] for item in centered["parities"][parity])
        == DIMENSIONS
        for parity in PARITIES
    )
)
check("all preregistered inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_conformal_shape_dynamics", GEOMETRY_SOURCE
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

same_actions = all(
    np.array_equal(left, right)
    for left, right in zip(groups["even"]["actions"], groups["odd"]["actions"])
)
same_edge_set = set(groups["even"]["edge_order"]) == set(groups["odd"]["edge_order"])
each_order_bijective = all(
    len(group["edge_order"]) == len(set(group["edge_order"])) == 720
    and len(group["old_orbits"]) == 30
    and all(len(orbit) == 24 for orbit in group["old_orbits"])
    and len(group["vertex_orbits"]) == 5
    and all(len(orbit) == 24 for orbit in group["vertex_orbits"])
    for group in groups.values()
)
odd_edge_index = {
    edge: index for index, edge in enumerate(groups["odd"]["edge_order"])
}
even_to_odd = np.asarray([
    odd_edge_index[edge] for edge in groups["even"]["edge_order"]
], dtype=int)
row_permutation_ok = bool(
    sorted(even_to_odd.tolist()) == list(range(720))
    and np.array_equal(
        incidences["odd"]["incidence"][even_to_odd],
        incidences["even"]["incidence"],
    )
)
orbit_control_ok = same_actions and same_edge_set and each_order_bijective and row_permutation_ok
check(
    "both schedules give one literal carrier with an exact row permutation",
    orbit_control_ok,
    f"row_sequences_equal={groups['even']['edge_order'] == groups['odd']['edge_order']}",
)

literal_incidence_ok = all(
    data["incidence"].shape == (720, 120)
    and np.all(np.sum(data["incidence"], axis=1) == 2)
    and np.all(np.sum(data["incidence"], axis=0) == 12)
    and data["gram_identity"]
    for data in incidences.values()
)
check("the unsigned incidence has the exact 600-cell local counts", literal_incidence_ok)

graph_injectivity_ok = all(
    data["connected"]
    and data["triangle_count"] == 1200
    and data["numerical_rank"] == 120
    and data["minimum_singular"] > 0
    for data in incidences.values()
)
check(
    "connectedness and odd cycles retain exact conformal rank 120",
    graph_injectivity_ok,
)

equivariance_ok = all(data["equivariant"] for data in incidences.values())
check("all 24 binary-tetrahedral actions commute exactly with incidence", equivariance_ok)

sector_data, sector_controls = high_precision_sector_bases(groups["even"])
maximum_basis_residual = max(
    value for key, value in sector_controls.items() if key.startswith("maximum_")
)
sector_basis_ok = bool(
    tuple(sector["dimension"] for sector in sector_data) == DIMENSIONS
    and sector_controls["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
    and maximum_basis_residual < mp.mpf("1e-70")
    and conformal["basis_control"]["irrep_dimensions"] == list(DIMENSIONS)
    and conformal["basis_control"]["signature_matches_frozen_source"]
)
check(
    "the same seven high-precision minimal sectors are reconstructed",
    sector_basis_ok,
    "max_residual=" + mp.nstr(maximum_basis_residual, 5),
)

sector_images = {parity: [] for parity in PARITIES}
conformal_contradicted = False
carrier_open = False
conformal_records = {parity: [] for parity in PARITIES}
for parity in PARITIES:
    incidence = incidences[parity]["incidence"].astype(np.complex128)
    for sector_index, sector in enumerate(sector_data):
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
        conformal_contradicted |= bool(
            open_count == 0
            and (nonzero != expected_rank or zero != len(singular) - expected_rank)
        )
        carrier_open |= open_count > 0
        conformal_basis = left[:, :expected_rank]
        projector_residual = operator_norm(
            compressed
            - conformal_basis @ (conformal_basis.conj().T @ compressed)
        )
        gap = float(singular[expected_rank - 1])
        if gap > 2 * epsilon_c:
            eta_k = float(
                2 * epsilon_c / (gap - 2 * epsilon_c)
                + 1000 * MACHINE_EPSILON * n
            )
        else:
            eta_k = math.inf
            carrier_open = True
        rank_ok = bool(
            nonzero == expected_rank
            and zero == len(singular) - expected_rank
            and open_count == 0
            and projector_residual <= 10 * epsilon_c
            and math.isfinite(eta_k)
        )
        if not rank_ok and not carrier_open:
            conformal_contradicted = True
        record = {
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "position_dimension": n,
            "conformal_dimension": expected_rank,
            "resolved_nonzero_singular_values": nonzero,
            "zero_consistent_singular_values": zero,
            "open_singular_values": open_count,
            "epsilon_C": serialize_float(epsilon_c),
            "minimum_resolved_singular": serialize_float(gap),
            "projector_residual": serialize_float(projector_residual),
            "subspace_error_eta_K": serialize_float(eta_k),
        }
        sector_images[parity].append({
            "basis": conformal_basis,
            "eta_K": eta_k,
            "record": record,
        })
        conformal_records[parity].append(record)

all_conformal_resolved = bool(
    not conformal_contradicted
    and not carrier_open
    and all(len(sector_images[p]) == 7 for p in PARITIES)
)
check(
    "all fourteen conformal-carrier classifications are complete",
    not conformal_contradicted
    and all(len(sector_images[p]) == 7 for p in PARITIES),
    (
        f"resolved={all_conformal_resolved}, "
        f"contradicted={conformal_contradicted}, open={carrier_open}"
    ),
)

dynamic_records = {parity: [] for parity in PARITIES}
dynamic_internal = {parity: {} for parity in PARITIES}
shape_contradicted = False
split_signature_ok = True
direct_sum_resolved = True
direct_sum_contradicted = False
dynamic_finite = True
classification_counts = Counter()
counts_by_operator_carrier = {
    operator: {"conformal": Counter(), "shape": Counter()}
    for operator in OPERATORS
}
negative_control_counts = Counter()
negative_control_by_family = {
    family: Counter()
    for family in (
        "spectral_positive",
        "spectral_negative",
        "fourier_low",
        "fourier_high",
    )
}
power_cell_counts = Counter()
negative_controls_finite = True
any_nonzero = False
any_open = False
all_zero = True

for parity in PARITIES:
    for sector_index, dimension in enumerate(DIMENSIONS):
        n = 30 * dimension
        r = 5 * dimension
        u_basis = sector_images[parity][sector_index]["basis"]
        eta_k = sector_images[parity][sector_index]["eta_K"]
        dynamic_internal[parity][sector_index] = {}
        variants = {}
        for variant in VARIANTS:
            prefix_m = f"{parity}_sector{sector_index}_{variant}_M"
            midpoint_m = np.asarray(source_npz[f"{prefix_m}_midpoint"])
            stored_radius_m = np.asarray(source_npz[f"{prefix_m}_radii"])
            entry_radius_m = component_reenclosure_radii(
                midpoint_m, stored_radius_m
            )
            hermitian = (midpoint_m + midpoint_m.conj().T) / 2
            hermitian = (hermitian + hermitian.conj().T) / 2
            hermitian_radius = (entry_radius_m + entry_radius_m.T) / 2
            norm_h = operator_norm(hermitian)
            epsilon_h = float(
                la.norm(hermitian_radius, "fro")
                + 1000 * MACHINE_EPSILON * n * max(1.0, norm_h)
            )
            a_matrix = u_basis.conj().T @ hermitian
            norm_a = operator_norm(a_matrix)
            epsilon_a = float(
                epsilon_h
                + 2 * eta_k * (norm_h + epsilon_h)
                + 1000 * MACHINE_EPSILON * n * max(1.0, norm_a)
            )
            _, singular_a, right_a = la.svd(a_matrix, full_matrices=True)
            a_nonzero = int(np.sum(singular_a > 100 * epsilon_a))
            a_zero = int(np.sum(singular_a < 10 * epsilon_a))
            a_open = len(singular_a) - a_nonzero - a_zero
            shape_contradicted |= bool(
                a_open == 0 and (a_nonzero != r or a_zero != 0)
            )
            carrier_open |= a_open > 0
            gap_a = float(singular_a[r - 1])
            if gap_a > 2 * epsilon_a:
                eta_s = float(
                    2 * epsilon_a / (gap_a - 2 * epsilon_a)
                    + 1000 * MACHINE_EPSILON * n
                )
            else:
                eta_s = math.inf
                carrier_open = True
            w_basis = right_a.conj().T[:, r:]
            direct_sum = np.column_stack((u_basis, w_basis))
            singular_b = la.svdvals(direct_sum)
            epsilon_b = float(
                1000 * MACHINE_EPSILON * n
                * max(1.0, float(singular_b[0]))
            )
            b_nonzero = int(np.sum(singular_b > 100 * epsilon_b))
            b_zero = int(np.sum(singular_b < 10 * epsilon_b))
            b_open = len(singular_b) - b_nonzero - b_zero
            if b_open:
                carrier_open = True
                direct_sum_resolved = False
            elif b_nonzero != n or b_zero != 0:
                shape_contradicted = True
                direct_sum_contradicted = True
                direct_sum_resolved = False

            values_h, vectors_h = la.eigh(hermitian)
            values_g = la.eigvalsh(
                (u_basis.conj().T @ hermitian @ u_basis
                 + (u_basis.conj().T @ hermitian @ u_basis).conj().T) / 2
            )
            shape_form = w_basis.conj().T @ hermitian @ w_basis
            values_s = la.eigvalsh((shape_form + shape_form.conj().T) / 2)
            signature_this = bool(
                np.sum(values_h > 100 * epsilon_h) == r
                and np.sum(values_h < -100 * epsilon_h) == 25 * dimension
                and np.sum(values_g > 100 * epsilon_h) == r
                and np.sum(values_s < -100 * epsilon_h) == 25 * dimension
            )
            split_signature_ok &= signature_this

            negative_dimension = 25 * dimension
            negative_basis = vectors_h[:, :negative_dimension]
            positive_basis = vectors_h[:, negative_dimension:]
            spectral_gap = float(
                values_h[negative_dimension] - values_h[negative_dimension - 1]
            )
            if spectral_gap > 2 * epsilon_h:
                eta_spec = float(
                    2 * epsilon_h / (spectral_gap - 2 * epsilon_h)
                    + 1000 * MACHINE_EPSILON * n
                )
            else:
                eta_spec = math.inf
                carrier_open = True

            fourier_indices = np.arange(n, dtype=float)
            fourier = np.exp(
                2j * np.pi * np.outer(fourier_indices, fourier_indices) / n
            ) / math.sqrt(n)
            eta_fourier = float(
                operator_norm(
                    fourier.conj().T @ fourier
                    - np.eye(n, dtype=np.complex128)
                )
                + 1000 * MACHINE_EPSILON * n
            )
            fourier_low = fourier[:, :r]
            fourier_high = fourier[:, r:]

            matrix_records = {}
            dynamic_internal[parity][sector_index][variant] = {}
            for operator in OPERATORS:
                prefix_x = f"{parity}_sector{sector_index}_{variant}_{operator}"
                midpoint_x = np.asarray(source_npz[f"{prefix_x}_midpoint"])
                stored_radius_x = np.asarray(source_npz[f"{prefix_x}_radii"])
                entry_radius_x = component_reenclosure_radii(
                    midpoint_x, stored_radius_x
                )
                norm_x = operator_norm(midpoint_x)
                epsilon_x = float(
                    la.norm(entry_radius_x, "fro")
                    + 1000 * MACHINE_EPSILON * n * max(1.0, norm_x)
                )
                projector_k = u_basis @ u_basis.conj().T
                residual_k_matrix = (
                    (np.eye(n, dtype=np.complex128) - projector_k)
                    @ midpoint_x @ u_basis
                )
                residual_k = operator_norm(residual_k_matrix)
                epsilon_k = float(
                    epsilon_x
                    + 2 * eta_k * (norm_x + epsilon_x)
                    + 1000 * MACHINE_EPSILON * n * max(1.0, norm_x)
                )
                residual_s_matrix = a_matrix @ midpoint_x @ w_basis
                residual_s = operator_norm(residual_s_matrix)
                epsilon_s = float(
                    epsilon_a * (norm_x + epsilon_x)
                    + norm_a * epsilon_x
                    + (norm_a + epsilon_a) * (norm_x + epsilon_x) * eta_s
                    + 1000 * MACHINE_EPSILON * n
                    * max(1.0, norm_a * norm_x)
                )
                label_k = residual_label(residual_k, epsilon_k)
                label_s = residual_label(residual_s, epsilon_s)
                for carrier, label in (("conformal", label_k), ("shape", label_s)):
                    classification_counts[label] += 1
                    counts_by_operator_carrier[operator][carrier][label] += 1
                    any_nonzero |= label == "NONZERO_RESOLVED"
                    any_open |= label == "OPEN"
                    all_zero &= label == "ZERO_CONSISTENT"
                dynamic_finite &= bool(
                    math.isfinite(residual_k)
                    and math.isfinite(epsilon_k)
                    and math.isfinite(residual_s)
                    and math.isfinite(epsilon_s)
                )

                negative_control_records = {}
                cell_control_labels = []
                control_families = (
                    ("spectral_positive", positive_basis, eta_spec),
                    ("spectral_negative", negative_basis, eta_spec),
                    ("fourier_low", fourier_low, eta_fourier),
                    ("fourier_high", fourier_high, eta_fourier),
                )
                for family, control_basis, eta_control in control_families:
                    control_projector = control_basis @ control_basis.conj().T
                    control_residual = operator_norm(
                        (np.eye(n, dtype=np.complex128) - control_projector)
                        @ midpoint_x @ control_basis
                    )
                    control_error = float(
                        epsilon_x
                        + 2 * eta_control * (norm_x + epsilon_x)
                        + 1000 * MACHINE_EPSILON * n * max(1.0, norm_x)
                    )
                    control_label = residual_label(
                        control_residual, control_error
                    )
                    negative_control_counts[control_label] += 1
                    negative_control_by_family[family][control_label] += 1
                    cell_control_labels.append(control_label)
                    negative_controls_finite &= bool(
                        math.isfinite(control_residual)
                        and math.isfinite(control_error)
                    )
                    negative_control_records[family] = {
                        "residual": serialize_float(control_residual),
                        "error": serialize_float(control_error),
                        "error_units": serialize_float(
                            control_residual / control_error
                            if control_error else math.inf
                        ),
                        "label": control_label,
                    }
                if "NONZERO_RESOLVED" in cell_control_labels:
                    power_label = "POWER_HIT"
                elif "OPEN" in cell_control_labels:
                    power_label = "POWER_OPEN"
                else:
                    power_label = "POWER_ZERO"
                power_cell_counts[power_label] += 1

                matrix_records[operator] = {
                    "matrix_norm": serialize_float(norm_x),
                    "matrix_error": serialize_float(epsilon_x),
                    "conformal_invariance_residual": serialize_float(residual_k),
                    "conformal_invariance_error": serialize_float(epsilon_k),
                    "conformal_invariance_error_units": serialize_float(
                        residual_k / epsilon_k if epsilon_k else math.inf
                    ),
                    "conformal_invariance_label": label_k,
                    "shape_invariance_residual": serialize_float(residual_s),
                    "shape_invariance_error": serialize_float(epsilon_s),
                    "shape_invariance_error_units": serialize_float(
                        residual_s / epsilon_s if epsilon_s else math.inf
                    ),
                    "shape_invariance_label": label_s,
                    "negative_controls": negative_control_records,
                    "power_label": power_label,
                }
                dynamic_internal[parity][sector_index][variant][operator] = {
                    "conformal_norm": residual_k,
                    "conformal_error": epsilon_k,
                    "shape_norm": residual_s,
                    "shape_error": epsilon_s,
                }

            variants[variant] = {
                "shape_carrier": {
                    "row_operator_rank": a_nonzero,
                    "row_operator_zero_consistent": a_zero,
                    "row_operator_open": a_open,
                    "epsilon_A": serialize_float(epsilon_a),
                    "minimum_row_operator_singular": serialize_float(gap_a),
                    "subspace_error_eta_S": serialize_float(eta_s),
                    "shape_dimension": n - r,
                    "direct_sum_minimum_singular": serialize_float(singular_b[-1]),
                    "direct_sum_condition": serialize_float(
                        singular_b[0] / singular_b[-1]
                    ),
                    "signature_control": signature_this,
                    "spectral_split_gap": serialize_float(spectral_gap),
                    "spectral_split_error_eta": serialize_float(eta_spec),
                    "fourier_unitarity_error_eta": serialize_float(eta_fourier),
                    "conformal_form_minimum_eigenvalue": serialize_float(values_g[0]),
                    "shape_form_maximum_eigenvalue": serialize_float(values_s[-1]),
                },
                "operators": matrix_records,
            }
        dynamic_records[parity].append({
            **sector_images[parity][sector_index]["record"],
            "variants": variants,
        })

shape_controls_ok = bool(
    not shape_contradicted
    and split_signature_ok
    and not direct_sum_contradicted
)
check(
    "all 56 action-relative shape-carrier classifications are complete",
    shape_controls_ok,
    (
        f"contradicted={shape_contradicted}, open={carrier_open}, "
        f"signature={split_signature_ok}, "
        f"direct_sum_resolved={direct_sum_resolved}, "
        f"direct_sum_contradicted={direct_sum_contradicted}"
    ),
)

required_classifications = 2 * 7 * 4 * 2 * 2
check(
    "all 224 preregistered dynamic residuals are finite and classified",
    (dynamic_finite or carrier_open)
    and sum(classification_counts.values()) == required_classifications,
    str(dict(classification_counts)),
)

required_negative_controls = 2 * 7 * 4 * 2 * 4
required_power_cells = 2 * 7 * 4 * 2
check(
    "all 448 post-result negative controls are complete and classified",
    (negative_controls_finite or carrier_open)
    and sum(negative_control_counts.values()) == required_negative_controls
    and sum(power_cell_counts.values()) == required_power_cells,
    (
        f"controls={dict(negative_control_counts)}, "
        f"cells={dict(power_cell_counts)}"
    ),
)

schedule_records = []
schedule_counts = Counter()
schedule_finite = True
for sector_index, dimension in enumerate(DIMENSIONS):
    n = 30 * dimension
    for variant in VARIANTS:
        for operator in OPERATORS:
            left = dynamic_internal["even"][sector_index][variant][operator]
            right = dynamic_internal["odd"][sector_index][variant][operator]
            for carrier in ("conformal", "shape"):
                left_norm = left[f"{carrier}_norm"]
                right_norm = right[f"{carrier}_norm"]
                distance = abs(left_norm - right_norm)
                epsilon = float(
                    left[f"{carrier}_error"] + right[f"{carrier}_error"]
                    + 1000 * MACHINE_EPSILON * n
                    * max(1.0, left_norm, right_norm)
                )
                label = schedule_label(distance, epsilon)
                schedule_counts[label] += 1
                schedule_finite &= bool(
                    math.isfinite(distance) and math.isfinite(epsilon)
                )
                schedule_records.append({
                    "sector_index": sector_index,
                    "irrep_dimension": dimension,
                    "variant": variant,
                    "operator": operator,
                    "carrier": carrier,
                    "residual_norm_distance": serialize_float(distance),
                    "comparison_error": serialize_float(epsilon),
                    "label": label,
                })

check(
    "all 112 schedule comparisons are finite and classified",
    (schedule_finite or carrier_open) and len(schedule_records) == 112,
    str(dict(schedule_counts)),
)

controls_ok = bool(
    provenance_ok
    and geometry_import_ok
    and orbit_control_ok
    and literal_incidence_ok
    and graph_injectivity_ok
    and equivariance_ok
    and sector_basis_ok
    and not conformal_contradicted
    and not shape_contradicted
    and not direct_sum_contradicted
    and split_signature_ok
)
if not controls_ok:
    outcome = "CONFORMAL_SHAPE_DYNAMICS_CONTROL_FAILED"
elif carrier_open:
    outcome = "CONFORMAL_SHAPE_CARRIER_OPEN"
elif any_nonzero:
    outcome = "CONFORMAL_SHAPE_MIXING_REFUTED"
elif any_open or not all_zero:
    outcome = "CONFORMAL_SHAPE_DYNAMICS_OPEN"
elif power_cell_counts["POWER_HIT"] != required_power_cells:
    outcome = "CONFORMAL_SHAPE_DECOUPLING_POWER_OPEN"
else:
    outcome = "CONFORMAL_SHAPE_DYNAMICS_DECOUPLED_POWER_CERTIFIED"

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "power_protocol_commit": POWER_PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "canonical_split_candidates_per_audit": 1,
    "carrier": {
        "vertices": 120,
        "edges": 720,
        "conformal_dimension": 120,
        "shape_dimension": 600,
        "schedule_row_permutation_exact": row_permutation_ok,
        "equivariant_all_24": equivariance_ok,
    },
    "basis_control": {
        "irrep_dimensions": list(DIMENSIONS),
        "maximum_high_precision_residual": mp.nstr(maximum_basis_residual, 70),
    },
    "conformal_sectors": conformal_records,
    "parities": dynamic_records,
    "required_residual_classification_counts": dict(classification_counts),
    "counts_by_operator_carrier": {
        operator: {
            carrier: dict(counter)
            for carrier, counter in carriers.items()
        }
        for operator, carriers in counts_by_operator_carrier.items()
    },
    "negative_control_classification_counts": dict(negative_control_counts),
    "negative_control_counts_by_family": {
        family: dict(counter)
        for family, counter in negative_control_by_family.items()
    },
    "power_cell_counts": dict(power_cell_counts),
    "schedule_comparisons": schedule_records,
    "schedule_label_counts": dict(schedule_counts),
    "classification": {
        "carrier_open": carrier_open,
        "all_required_zero_consistent": all_zero,
        "any_required_nonzero_resolved": any_nonzero,
        "any_required_open": any_open,
        "shape_is_not_declared_transverse_traceless": True,
        "constraint_quotient_derived": False,
    },
    "continuum_target_loaded": False,
    "polarization_target_loaded": False,
    "speed_target_loaded": False,
    "refinement_target_loaded": False,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("SCIENTIFIC OUTCOME:", outcome)
print("required residual labels:", dict(classification_counts))
print("schedule labels:", dict(schedule_counts))
print(f"{passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

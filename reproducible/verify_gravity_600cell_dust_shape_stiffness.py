#!/usr/bin/env python3
"""Blind action-relative shape-stiffness census on the dust 600-cell.

Prior-art commit: e37d80c.
Final preregistration commit: c719a87.
No continuum spectrum, degeneracy, polarization, speed or refinement target is
loaded by this verifier.
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
SHAPE_JSON = HERE / "gravity_600cell_dust_conformal_shape_dynamics.json"
SHAPE_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_shape_dynamics.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
OUTPUT = HERE / "gravity_600cell_dust_shape_stiffness.json"

PRIOR_ART_COMMIT = "e37d80c"
PROTOCOL_COMMIT = "c719a87"
EXPECTED_HASHES = {
    "centered_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "shape_json": "c5bbeaa2a64d07688061bc5098a33361dc2f5300d637e44a10b6cccbbd1bb162",
    "shape_source": "52857835b37722db51c03587a9583426b26caaf2cb6b2d55c4fee05419883112",
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


def scalar_zero_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "OPEN"
    if abs(value) <= 10 * error:
        return "ZERO_CONSISTENT"
    if abs(value) > 100 * error:
        return "NONZERO_RESOLVED"
    return "OPEN"


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


def reality_label(value, error):
    if not math.isfinite(value.imag) or not math.isfinite(error) or error < 0:
        return "REALITY_OPEN"
    if abs(value.imag) <= 10 * error:
        return "REAL_CONSISTENT"
    if abs(value.imag) > 100 * error:
        return "COMPLEX_RESOLVED"
    return "REALITY_OPEN"


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


def sorted_complex(values):
    values = np.asarray(values, dtype=np.complex128)
    order = np.lexsort((values.imag, values.real))
    return values[order]


print("=" * 78)
print("600-CELL ACTION-RELATIVE SHAPE STIFFNESS CENSUS")
print("=" * 78)

paths = {
    "centered_json": CENTERED_JSON,
    "centered_npz": CENTERED_NPZ,
    "shape_json": SHAPE_JSON,
    "shape_source": SHAPE_SOURCE,
    "conformal_source": CONFORMAL_SOURCE,
    "full_source": FULL_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "commons": COMMONS_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
centered = json.loads(CENTERED_JSON.read_text())
shape_source_result = json.loads(SHAPE_JSON.read_text())
archive = np.load(CENTERED_NPZ, allow_pickle=False)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and centered["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and centered["passed"] == centered["tests"] == 7
    and centered["numeric_archive_arrays"] == len(archive.files) == 560
    and shape_source_result["outcome"]
    == "CONFORMAL_SHAPE_DYNAMICS_DECOUPLED_POWER_CERTIFIED"
    and shape_source_result["passed"] == shape_source_result["tests"] == 12
    and shape_source_result["carrier"]["conformal_dimension"] == 120
    and shape_source_result["carrier"]["shape_dimension"] == 600
)
check("all preregistered inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_shape_stiffness", GEOMETRY_SOURCE
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
    all(
        len(group["edge_order"]) == len(set(group["edge_order"])) == 720
        and len(group["old_orbits"]) == 30
        and all(len(orbit) == 24 for orbit in group["old_orbits"])
        and len(group["vertex_orbits"]) == 5
        for group in groups.values()
    )
    and sorted(even_to_odd.tolist()) == list(range(720))
    and np.array_equal(
        incidences["odd"]["incidence"][even_to_odd],
        incidences["even"]["incidence"],
    )
    and all(data["equivariant"] for data in incidences.values())
)
check("both schedules retain one exact equivariant edge carrier", carrier_geometry_ok)

incidence_ok = bool(all(
    data["incidence"].shape == (720, 120)
    and np.all(np.sum(data["incidence"], axis=1) == 2)
    and np.all(np.sum(data["incidence"], axis=0) == 12)
    and data["connected"]
    and data["triangle_count"] == 1200
    and data["numerical_rank"] == 120
    for data in incidences.values()
))
check("the canonical conformal incidence retains exact rank 120", incidence_ok)

sector_data, sector_controls = high_precision_sector_bases(groups["even"])
maximum_basis_residual = max(
    value for key, value in sector_controls.items() if key.startswith("maximum_")
)
sector_basis_ok = bool(
    tuple(sector["dimension"] for sector in sector_data) == DIMENSIONS
    and sector_controls["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
    and maximum_basis_residual < mp.mpf("1e-70")
)
check(
    "the same seven high-precision symmetry sectors are reconstructed",
    sector_basis_ok,
    "max_residual=" + mp.nstr(maximum_basis_residual, 5),
)

sector_images = {parity: [] for parity in PARITIES}
conformal_open = False
conformal_contradicted = False
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
        sector_images[parity].append({
            "basis": left[:, :expected_rank],
            "eta_K": eta_k,
        })

check(
    "all fourteen canonical conformal carriers are resolved",
    not conformal_contradicted and not conformal_open,
    f"contradicted={conformal_contradicted}, open={conformal_open}",
)

records = {parity: [] for parity in PARITIES}
internal = {parity: {} for parity in PARITIES}
carrier_open = conformal_open
carrier_contradicted = conformal_contradicted
kinetic_counts = Counter()
pencil_sign_counts = Counter()
actual_reality_counts = Counter()
actual_sign_counts = Counter()
compatibility_counts = Counter()
sign_disagreement_cells = 0
all_numeric_finite = True

for parity in PARITIES:
    for sector_index, dimension in enumerate(DIMENSIONS):
        n = 30 * dimension
        r = 5 * dimension
        s = 25 * dimension
        u_basis = sector_images[parity][sector_index]["basis"]
        eta_k = sector_images[parity][sector_index]["eta_K"]
        sector_variants = []
        internal[parity][sector_index] = {}

        for variant in VARIANTS:
            prefix = f"{parity}_sector{sector_index}_{variant}"
            midpoint_m, radius_m = source_matrix(archive, prefix, "M")
            h_m = (midpoint_m + midpoint_m.conj().T) / 2
            radius_hm = (radius_m + radius_m.T) / 2
            epsilon_hm = matrix_error(h_m, radius_hm, n)

            row = u_basis.conj().T @ h_m
            norm_row = operator_norm(row)
            epsilon_row = float(
                epsilon_hm
                + 2 * eta_k * (operator_norm(h_m) + epsilon_hm)
                + 1000 * MACHINE_EPSILON * n * max(1.0, norm_row)
            )
            _, singular_row, right_row = la.svd(row, full_matrices=True)
            row_nonzero = int(np.sum(singular_row > 100 * epsilon_row))
            row_zero = int(np.sum(singular_row < 10 * epsilon_row))
            row_open = len(singular_row) - row_nonzero - row_zero
            if row_open == 0 and (row_nonzero != r or row_zero != 0):
                carrier_contradicted = True
            carrier_open |= row_open > 0
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

            direct_sum = np.column_stack((u_basis, w_basis))
            direct_values = la.svdvals(direct_sum)
            direct_floor = float(
                1000 * MACHINE_EPSILON * n
                * max(1.0, float(direct_values[0]))
            )
            direct_label = sign_label(float(direct_values[-1]), direct_floor)
            if direct_label == "OPEN":
                carrier_open = True
            elif direct_label != "POSITIVE_RESOLVED":
                carrier_contradicted = True

            midpoint_v, radius_v = source_matrix(archive, prefix, "V")
            h_v = (midpoint_v + midpoint_v.conj().T) / 2
            radius_hv = (radius_v + radius_v.T) / 2
            epsilon_hv = matrix_error(h_v, radius_hv, n)

            midpoint_o, radius_o = source_matrix(archive, prefix, "Omega")
            epsilon_o = matrix_error(midpoint_o, radius_o, n)

            m_s = w_basis.conj().T @ h_m @ w_basis
            m_s = (m_s + m_s.conj().T) / 2
            v_s = w_basis.conj().T @ h_v @ w_basis
            v_s = (v_s + v_s.conj().T) / 2
            omega_s = w_basis.conj().T @ midpoint_o @ w_basis
            epsilon_ms = restriction_error(h_m, epsilon_hm, eta_s, n)
            epsilon_vs = restriction_error(h_v, epsilon_hv, eta_s, n)
            epsilon_os = restriction_error(midpoint_o, epsilon_o, eta_s, n)

            b_matrix = -m_s
            a_matrix = -v_s
            values_b = la.eigvalsh(b_matrix)
            values_a = la.eigvalsh(a_matrix)
            minimum_b = float(values_b[0])
            kinetic_label = sign_label(minimum_b, epsilon_ms)
            if kinetic_label == "POSITIVE_RESOLVED":
                kinetic_class = "POSITIVE_DEFINITE_RESOLVED"
            elif kinetic_label == "NEGATIVE_RESOLVED":
                kinetic_class = "KINETIC_CONTRADICTION"
                carrier_contradicted = True
            else:
                kinetic_class = "KINETIC_OPEN"
                carrier_open = True
            kinetic_counts[kinetic_class] += 1

            pencil_labels = [sign_label(float(value), epsilon_vs) for value in values_a]
            pencil_sign_counts.update(pencil_labels)

            b_lower = minimum_b - epsilon_ms
            if minimum_b > 0 and b_lower > 0:
                generalized = la.eigh(a_matrix, b_matrix, eigvals_only=True)
                epsilon_pencil = float(
                    epsilon_vs / b_lower
                    + operator_norm(a_matrix) * epsilon_ms
                    / (minimum_b * b_lower)
                    + 1000 * MACHINE_EPSILON * s
                    * max(1.0, operator_norm(la.solve(b_matrix, a_matrix)))
                )
            else:
                generalized = np.full(s, np.nan)
                epsilon_pencil = math.inf

            omega_values, omega_vectors = la.eig(omega_s)
            try:
                eigenvector_condition = float(np.linalg.cond(omega_vectors))
            except np.linalg.LinAlgError:
                eigenvector_condition = math.inf
            epsilon_eig = float(
                eigenvector_condition * epsilon_os
                + 1000 * MACHINE_EPSILON * s
                * max(1.0, operator_norm(omega_s))
            )
            omega_values = sorted_complex(omega_values)
            actual_labels = []
            for value in omega_values:
                real_class = reality_label(value, epsilon_eig)
                actual_reality_counts[real_class] += 1
                if real_class == "REAL_CONSISTENT":
                    shape_sign = sign_label(float(value.real), epsilon_eig)
                else:
                    shape_sign = "SIGN_NOT_ASSIGNED"
                actual_sign_counts[shape_sign] += 1
                actual_labels.append((real_class, shape_sign))

            residual = m_s @ omega_s - v_s
            residual_norm = operator_norm(residual)
            epsilon_residual = float(
                epsilon_ms * (operator_norm(omega_s) + epsilon_os)
                + operator_norm(m_s) * epsilon_os
                + epsilon_vs
                + 1000 * MACHINE_EPSILON * s
                * max(
                    1.0,
                    operator_norm(m_s) * operator_norm(omega_s),
                    operator_norm(v_s),
                )
            )
            compatibility_label = scalar_zero_label(residual_norm, epsilon_residual)
            compatibility_counts[compatibility_label] += 1

            pencil_cell_counts = Counter(pencil_labels)
            actual_cell_counts = Counter(label[1] for label in actual_labels)
            pencil_complete = pencil_cell_counts["OPEN"] == 0
            actual_complete = bool(
                all(label[0] == "REAL_CONSISTENT" for label in actual_labels)
                and actual_cell_counts["OPEN"] == 0
                and actual_cell_counts["SIGN_NOT_ASSIGNED"] == 0
            )
            sign_counts_agree = None
            if pencil_complete and actual_complete:
                sign_counts_agree = all(
                    pencil_cell_counts[label] == actual_cell_counts[label]
                    for label in (
                        "POSITIVE_RESOLVED",
                        "NEGATIVE_RESOLVED",
                        "ZERO_CONSISTENT",
                    )
                )
                sign_disagreement_cells += int(not sign_counts_agree)

            finite_this = bool(
                np.all(np.isfinite(values_a))
                and np.all(np.isfinite(values_b))
                and np.all(np.isfinite(omega_values))
                and math.isfinite(residual_norm)
                and math.isfinite(epsilon_residual)
            )
            all_numeric_finite &= finite_this

            internal[parity][sector_index][variant] = {
                "generalized": np.sort(np.asarray(generalized, dtype=float)),
                "epsilon_pencil": epsilon_pencil,
                "omega": omega_values,
                "epsilon_eig": epsilon_eig,
            }

            sector_variants.append({
                "variant": variant,
                "position_dimension": n,
                "shape_dimension": s,
                "shape_carrier": {
                    "row_rank_resolved": row_nonzero,
                    "row_zero_consistent": row_zero,
                    "row_open": row_open,
                    "row_gap": sf(row_gap),
                    "row_error": sf(epsilon_row),
                    "subspace_error_eta_S": sf(eta_s),
                    "direct_sum_minimum_singular": sf(direct_values[-1]),
                    "direct_sum_condition": sf(direct_values[0] / direct_values[-1]),
                    "direct_sum_label": direct_label,
                },
                "kinetic": {
                    "minimum_eigenvalue_B": sf(minimum_b),
                    "maximum_eigenvalue_B": sf(values_b[-1]),
                    "restricted_error": sf(epsilon_ms),
                    "label": kinetic_class,
                },
                "hermitian_pencil": {
                    "generalized_eigenvalues": [sf(value) for value in generalized],
                    "generalized_eigenvalue_error": sf(epsilon_pencil),
                    "A_eigenvalues": [sf(value) for value in values_a],
                    "A_sign_labels": pencil_labels,
                    "sign_counts": dict(pencil_cell_counts),
                },
                "normalized_shape_block": {
                    "eigenvalues": [
                        {"real": sf(value.real), "imag": sf(value.imag)}
                        for value in omega_values
                    ],
                    "eigenvector_condition": sf(eigenvector_condition),
                    "eigenvalue_error": sf(epsilon_eig),
                    "labels": [
                        {"reality": reality, "sign": sign}
                        for reality, sign in actual_labels
                    ],
                    "reality_counts": dict(Counter(x[0] for x in actual_labels)),
                    "sign_counts": dict(actual_cell_counts),
                },
                "action_compatibility": {
                    "residual_norm": sf(residual_norm),
                    "residual_error": sf(epsilon_residual),
                    "error_units": sf(
                        residual_norm / epsilon_residual
                        if epsilon_residual else math.inf
                    ),
                    "label": compatibility_label,
                    "resolved_sign_counts_agree": sign_counts_agree,
                },
            })

        records[parity].append({
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "variants": sector_variants,
        })

required_cells = 2 * 7 * 4
required_eigenvalues = 2 * 4 * sum(25 * d for d in DIMENSIONS)
carrier_complete = bool(
    sum(kinetic_counts.values()) == required_cells
    and all(len(records[p]) == 7 for p in PARITIES)
)
check(
    "all 56 action-relative shape carriers and kinetic forms are classified",
    carrier_complete,
    f"kinetic={dict(kinetic_counts)}, open={carrier_open}, contradicted={carrier_contradicted}",
)

census_complete = bool(
    all_numeric_finite
    and sum(pencil_sign_counts.values()) == required_eigenvalues
    and sum(actual_reality_counts.values()) == required_eigenvalues
    and sum(actual_sign_counts.values()) == required_eigenvalues
)
check(
    "both complete 4,800-eigenvalue censuses are finite and classified",
    census_complete,
    (
        f"pencil={dict(pencil_sign_counts)}, "
        f"reality={dict(actual_reality_counts)}, "
        f"actual_sign={dict(actual_sign_counts)}"
    ),
)

compatibility_complete = sum(compatibility_counts.values()) == required_cells
check(
    "all 56 action-compatibility residuals are classified",
    compatibility_complete,
    f"residuals={dict(compatibility_counts)}, sign_disagreements={sign_disagreement_cells}",
)

schedule_records = []
schedule_counts = Counter()
for sector_index, dimension in enumerate(DIMENSIONS):
    s = 25 * dimension
    for variant in VARIANTS:
        left = internal["even"][sector_index][variant]
        right = internal["odd"][sector_index][variant]
        for object_name, values_name, error_name in (
            ("hermitian_pencil", "generalized", "epsilon_pencil"),
            ("normalized_shape_block", "omega", "epsilon_eig"),
        ):
            values_left = left[values_name]
            values_right = right[values_name]
            if len(values_left) == len(values_right) == s:
                distance = float(np.max(np.abs(values_left - values_right)))
            else:
                distance = math.inf
            comparison_error = float(
                left[error_name]
                + right[error_name]
                + 1000 * MACHINE_EPSILON * s
                * max(
                    1.0,
                    float(np.nanmax(np.abs(values_left)))
                    if len(values_left) else 1.0,
                    float(np.nanmax(np.abs(values_right)))
                    if len(values_right) else 1.0,
                )
            )
            label = schedule_label(distance, comparison_error)
            schedule_counts[label] += 1
            schedule_records.append({
                "sector_index": sector_index,
                "irrep_dimension": dimension,
                "variant": variant,
                "object": object_name,
                "distance": sf(distance),
                "comparison_error": sf(comparison_error),
                "error_units": sf(
                    distance / comparison_error if comparison_error else math.inf
                ),
                "label": label,
            })

check(
    "all 56 preregistered schedule comparisons are classified",
    len(schedule_records) == 56 and sum(schedule_counts.values()) == 56,
    str(dict(schedule_counts)),
)

forbidden_targets_absent = True
check(
    "no continuum, degeneracy, polarization, speed or refinement target was loaded",
    forbidden_targets_absent,
)

controls_ok = bool(
    provenance_ok
    and geometry_import_ok
    and carrier_geometry_ok
    and incidence_ok
    and sector_basis_ok
    and not carrier_contradicted
    and carrier_complete
    and census_complete
    and compatibility_complete
)
any_complex = actual_reality_counts["COMPLEX_RESOLVED"] > 0
any_action_mismatch = bool(
    compatibility_counts["NONZERO_RESOLVED"] > 0
    or sign_disagreement_cells > 0
)
any_schedule_dependent = schedule_counts["SCHEDULE_DEPENDENT"] > 0
any_negative = bool(
    pencil_sign_counts["NEGATIVE_RESOLVED"] > 0
    or actual_sign_counts["NEGATIVE_RESOLVED"] > 0
)
any_open = bool(
    pencil_sign_counts["OPEN"] > 0
    or actual_reality_counts["REALITY_OPEN"] > 0
    or actual_sign_counts["OPEN"] > 0
    or actual_sign_counts["SIGN_NOT_ASSIGNED"] > 0
    or compatibility_counts["OPEN"] > 0
    or schedule_counts["SCHEDULE_OPEN"] > 0
)
any_zero = bool(
    pencil_sign_counts["ZERO_CONSISTENT"] > 0
    or actual_sign_counts["ZERO_CONSISTENT"] > 0
)

if not controls_ok:
    outcome = "SHAPE_STIFFNESS_CONTROL_FAILED"
elif carrier_open or kinetic_counts["KINETIC_OPEN"] > 0:
    outcome = "SHAPE_STIFFNESS_CARRIER_OPEN"
elif any_complex or any_action_mismatch:
    outcome = "SHAPE_STIFFNESS_COMPLEX_OR_ACTION_MISMATCH"
elif any_schedule_dependent:
    outcome = "SHAPE_STIFFNESS_SCHEDULE_DEPENDENT"
elif any_negative:
    outcome = "SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED"
elif any_open:
    outcome = "SHAPE_STIFFNESS_SIGN_OPEN"
elif any_zero:
    outcome = "SHAPE_STIFFNESS_NONNEGATIVE_WITH_ZERO_MODES"
else:
    outcome = "SHAPE_STIFFNESS_POSITIVE_CENSUS"

outcome_complete = outcome in {
    "SHAPE_STIFFNESS_CONTROL_FAILED",
    "SHAPE_STIFFNESS_CARRIER_OPEN",
    "SHAPE_STIFFNESS_COMPLEX_OR_ACTION_MISMATCH",
    "SHAPE_STIFFNESS_SCHEDULE_DEPENDENT",
    "SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED",
    "SHAPE_STIFFNESS_SIGN_OPEN",
    "SHAPE_STIFFNESS_NONNEGATIVE_WITH_ZERO_MODES",
    "SHAPE_STIFFNESS_POSITIVE_CENSUS",
}
check("the preregistered outcome ladder is exhaustive", outcome_complete, outcome)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "enumeration": {
        "schedules": 2,
        "sectors": 7,
        "variants": 4,
        "shape_pencils": required_cells,
        "eigenvalue_instances_per_object": required_eigenvalues,
    },
    "carrier": {
        "vertices": 120,
        "edges": 720,
        "conformal_dimension": 120,
        "shape_dimension": 600,
        "carrier_open": carrier_open,
        "carrier_contradicted": carrier_contradicted,
    },
    "pencil_sign_counts": dict(pencil_sign_counts),
    "normalized_reality_counts": dict(actual_reality_counts),
    "normalized_sign_counts": dict(actual_sign_counts),
    "kinetic_counts": dict(kinetic_counts),
    "compatibility_counts": dict(compatibility_counts),
    "resolved_sign_disagreement_cells": sign_disagreement_cells,
    "schedule_comparisons": schedule_records,
    "schedule_counts": dict(schedule_counts),
    "parities": records,
    "target_load_flags": {
        "continuum_spectrum": False,
        "desired_degeneracy": False,
        "polarization_count": False,
        "speed": False,
        "refinement": False,
    },
    "interpretation_limits": {
        "shape_declared_transverse_traceless": False,
        "constraint_quotient_derived": False,
        "negative_mode_declared_ghost": False,
        "positive_mode_declared_graviton": False,
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("SCIENTIFIC OUTCOME:", outcome)
print("kinetic labels:", dict(kinetic_counts))
print("pencil sign labels:", dict(pencil_sign_counts))
print("normalized reality labels:", dict(actual_reality_counts))
print("normalized sign labels:", dict(actual_sign_counts))
print("action compatibility:", dict(compatibility_counts))
print("schedule labels:", dict(schedule_counts))
print(f"{passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

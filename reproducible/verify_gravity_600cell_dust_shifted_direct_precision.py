#!/usr/bin/env python3
"""Target-disclosed direct high-precision shifted stiffness rank gate.

Prior-art commit: 52f337b.
Precision-attribution artifact commit: 67e99fb.
Preregistered protocol commit: 4db96b0.

The old rank-15 target and the shifted OPEN result are disclosed.  No old
stiffness eigenvalue or sign label is loaded.  The calculation removes only
the preregistered binary64 tangent-serialization interface.
"""

import ast
from collections import Counter, defaultdict
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
import json
import math
from pathlib import Path
import sys

from flint import acb, acb_mat, ctx
import mpmath as mp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIRST_TICK = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
SECOND_TICK = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
THIRD_TICK = HERE / "gravity_600cell_dust_third_tick_local_correction.json"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
OLD_CENTERED_JSON = HERE / "gravity_600cell_dust_shifted_centered.json"
OLD_CENTERED_NPZ = HERE / "gravity_600cell_dust_shifted_centered.npz"
PRECISION_AUDIT = HERE / "gravity_600cell_dust_shifted_precision_audit.json"
OUTPUT = HERE / "gravity_600cell_dust_shifted_direct_precision.json"

PRIOR_ART_COMMIT = "52f337b"
ATTRIBUTION_COMMIT = "67e99fb"
PROTOCOL_COMMIT = "4db96b0"
EXPECTED_HASHES = {
    "first_tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "second_tick": "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70",
    "third_tick": "ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "old_centered_json": "265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47",
    "old_centered_npz": "c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8",
    "precision_audit": "409e428ca05b4b6c6e380d7af6d84fc3834afc1fdb95065d6f0a05618e1d2cee",
}

DPS = 100
BALL_DPS = 80
mp.mp.dps = DPS
ctx.dps = BALL_DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-20"),
    "operational_shadow": mp.mpf("1e-15"),
    "validation_primary": mp.mpf("3e-20"),
    "validation_shadow": mp.mpf("3e-15"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
PARITIES = ("even", "odd")
TARGET_SECTORS = (4, 5)
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {hinge: index for index, hinge in enumerate(LOCAL_HINGES)}
I = mp.mpc(0, 1)
ARITHMETIC_FLOOR = mp.mpf("1e-70")
MACHINE_EPSILON = np.finfo(float).eps
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


def load_named_functions(source, wanted):
    tree = ast.parse(source.read_text(), filename=str(source))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch in {source.name}: {wanted-found}")
    exec(compile(ast.Module(body=body, type_ignores=[]), str(source), "exec"), globals())


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


def scale(matrix, factor):
    result = acb_mat(matrix.nrows(), matrix.ncols())
    scalar = acb(str(factor))
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            result[row, column] = scalar * matrix[row, column]
    return result


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


def all_entries_contain_zero(matrix):
    return all(
        matrix[row, column].contains(0)
        for row in range(matrix.nrows())
        for column in range(matrix.ncols())
    )


def operator_norm(matrix):
    values = la.svdvals(matrix)
    return float(values[0]) if len(values) else 0.0


def matrix_error(midpoint, radii, n):
    return float(
        la.norm(radii, "fro")
        + 1000 * MACHINE_EPSILON * n * max(1.0, operator_norm(midpoint))
    )


def restriction_error(midpoint, base_error, eta_shape, n):
    norm = operator_norm(midpoint)
    return float(
        base_error
        + 2 * eta_shape * (norm + base_error)
        + 1000 * MACHINE_EPSILON * n * max(1.0, norm)
    )


def sign_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "OPEN"
    if value > 100 * error:
        return "POSITIVE_RESOLVED"
    if value < -100 * error:
        return "NEGATIVE_RESOLVED"
    if abs(value) < 10 * error:
        return "ZERO_CONSISTENT"
    return "OPEN"


def scalar_zero_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "OPEN"
    if abs(value) < 10 * error:
        return "ZERO_CONSISTENT"
    if abs(value) > 100 * error:
        return "NONZERO_RESOLVED"
    return "OPEN"


def slab_index_data(model, a_old, a_new, r_value):
    data = dict(group_and_index_data(model, (mp.mpf(0), r_value)))
    rho = RHO0 * mp.exp(r_value)
    values = {
        "old": mp.exp(2 * a_old) * L0_SQUARE,
        "internal": mp.exp(a_old + a_new) * L0_SQUARE - rho,
        "pole": -rho,
        "new": mp.exp(2 * a_new) * L0_SQUARE,
    }
    if min(values["old"], values["internal"], rho, values["new"]) <= 0:
        raise RuntimeError("direct slab left the positive magnitude domain")
    data["rho"] = rho
    data["signed_base"] = tuple(values[kind] for kind in data["edge_kind"])
    return data, values


def boundary_mapping(index_data):
    mapping = []
    for old_type in range(30):
        shifted = tuple(
            tuple(vertex + 120 for vertex in edge)
            for edge in index_data["orbit_edges"][old_type]
        )
        matches = [
            final for final in range(30)
            if shifted == index_data["orbit_edges"][65 + final]
        ]
        if len(matches) != 1:
            raise RuntimeError(f"boundary mapping is not unique: {old_type}, {matches}")
        mapping.append(matches[0])
    return tuple(mapping)


def direct_and_old_matrix(matrix, old_archive, prefix, name):
    midpoint, stored = acb_midpoint_and_radii(matrix)
    radii = component_reenclosure_radii(midpoint, stored)
    old_midpoint = np.asarray(old_archive[f"{prefix}_{name}_midpoint"])
    old_stored = np.asarray(old_archive[f"{prefix}_{name}_radii"])
    old_radii = component_reenclosure_radii(old_midpoint, old_stored)
    tolerance = old_radii + radii
    difference = midpoint - old_midpoint
    overlap = bool(
        np.all(np.abs(np.real(difference)) <= tolerance)
        and np.all(np.abs(np.imag(difference)) <= tolerance)
    )
    new_norm = float(la.norm(radii, "fro"))
    old_norm = float(la.norm(old_radii, "fro"))
    ratio = old_norm / new_norm if new_norm > 0 else math.inf
    return midpoint, radii, overlap, ratio, old_norm, new_norm


paths = {
    "first_tick": FIRST_TICK,
    "second_tick": SECOND_TICK,
    "third_tick": THIRD_TICK,
    "geometry_source": GEOMETRY_SOURCE,
    "rank_source": RANK_SOURCE,
    "full_source": FULL_SOURCE,
    "conformal_source": CONFORMAL_SOURCE,
    "commons": COMMONS_SOURCE,
    "old_centered_json": OLD_CENTERED_JSON,
    "old_centered_npz": OLD_CENTERED_NPZ,
    "precision_audit": PRECISION_AUDIT,
}
hashes = {name: sha256(path) for name, path in paths.items()}
first_tick = json.loads(FIRST_TICK.read_text())
second_tick = json.loads(SECOND_TICK.read_text())
third_tick = json.loads(THIRD_TICK.read_text())
old_centered = json.loads(OLD_CENTERED_JSON.read_text())
old_archive = np.load(OLD_CENTERED_NPZ, allow_pickle=False)
precision_audit = json.loads(PRECISION_AUDIT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and first_tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and second_tick["outcome"] == "SECOND_HOMOTHETIC_TICK_ACCEPTED"
    and third_tick["outcome"] == "THIRD_HOMOTHETIC_TICK_ACCEPTED"
    and second_tick["fixed_mass"] is True
    and third_tick["fixed_mass"] is True
    and old_centered["outcome"] == "SHIFTED_CENTERED_CERTIFIED"
    and old_centered["numeric_archive_arrays"] == len(old_archive.files) == 560
    and precision_audit["outcome"]
    == "SHIFTED_PRECISION_BINARY_SERIALIZATION_DOMINANT"
    and precision_audit["direct_high_precision_reconstruction_authorized"] is True
)
check("all direct-precision inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_shifted_direct_precision", GEOMETRY_SOURCE
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

load_named_functions(RANK_SOURCE, {
    "orbit_sort_key", "augment_boundary_orbits", "log_minus",
    "signed_volume_square", "angle_record", "area_data",
    "extended_edge_image", "group_and_index_data", "prepare_geometry",
})
load_named_functions(FULL_SOURCE, {
    "norm2", "mp_frobenius", "mp_submatrix", "cluster_sorted",
    "high_precision_sector_bases", "high_precision_pattern_cache",
    "assemble_full_representative_kernels", "project_full_kernel",
    "mp_to_acb", "mp_matrix_to_acb", "acb_midpoint_and_radii",
    "expanded_types", "build_tangent_ball",
})
load_named_functions(CONFORMAL_SOURCE, {
    "edge_image", "group_data", "incidence_data", "mp_to_numpy",
    "component_reenclosure_radii",
})

M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON3 * L0
RHO0 = mp.mpf("0.0102") ** 2
models = {
    parity: augment_boundary_orbits(model)
    for parity, model in gro.models.items()
}

groups = {parity: group_data(gro.models[parity], gro) for parity in PARITIES}
incidences = {parity: incidence_data(groups[parity]) for parity in PARITIES}
incidence_ok = all(
    data["incidence"].shape == (720, 120)
    and data["equivariant"]
    and data["connected"]
    and data["numerical_rank"] == 120
    for data in incidences.values()
)
check("both exact equivariant incidence carriers retain rank 120", incidence_ok)

print("=" * 78)
print("DIRECT HIGH-PRECISION SHIFTED SHAPE-STIFFNESS RANK GATE")
print("=" * 78)

records = []
all_branch_controls = True
all_carrier_order = True
all_twists = True
all_principal_identities = True
all_overlap = True
all_precision_reduction = True
all_shape_carriers = True
all_kinetic = True
all_compatibility = True
aggregate_signs = Counter()

for parity in PARITIES:
    a1, _ = [mp.mpf(value) for value in first_tick["solutions"][parity]["state"]]
    a2, r2 = [
        mp.mpf(value)
        for value in second_tick["solutions"][parity]["state_absolute"]
    ]
    a3, r3 = [
        mp.mpf(value)
        for value in third_tick["solutions"][parity]["state_absolute"]
    ]
    slab_specs = (("slab2", a1, a2, r2), ("slab3", a2, a3, r3))
    slab_data = {}
    for slab_name, a_old, a_new, r_value in slab_specs:
        index_data, kind_values = slab_index_data(
            models[parity], a_old, a_new, r_value
        )
        geometry = prepare_geometry(models[parity], index_data)
        mapping = boundary_mapping(index_data)
        sectors, sector_control = high_precision_sector_bases(index_data)
        edge_order = tuple(
            edge for orbit in index_data["orbit_edges"][:30] for edge in orbit
        )
        carrier_order_ok = bool(
            edge_order == groups[parity]["edge_order"]
            and sorted(mapping) == list(range(30))
            and tuple(sector["dimension"] for sector in sectors)
            == (3, 2, 2, 2, 1, 1, 1)
            and sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        )
        all_carrier_order &= carrier_order_ok
        print(f"[{parity}] differentiating {slab_name} local patterns", flush=True)
        pattern_cache, branch = high_precision_pattern_cache(
            geometry["patterns"], kind_values
        )
        kernels, kernel_control = assemble_full_representative_kernels(
            index_data, geometry, pattern_cache
        )
        branch_ok = bool(
            branch["entry_pass"]
            and branch["base_negative_counts"] == Counter({1: 2400})
            and branch["minimum_leading_minor"] > 0
            and branch["minimum_argument"] > mp.mpf("1e-6")
            and kernel_control["maximum_imaginary"] < ARITHMETIC_FLOOR
            and len(set(kernel_control["nonzero_entries"].values())) == 1
        )
        all_branch_controls &= branch_ok
        slab_data[slab_name] = {
            "index_data": index_data,
            "mapping": mapping,
            "sectors": sectors,
            "kernels": kernels,
            "branch_ok": branch_ok,
            "carrier_order_ok": carrier_order_ok,
        }

    signature2 = tuple(
        (item["dimension"], mp.nstr(item["old_central_eigenvalue"], 70), item["splitter"])
        for item in slab_data["slab2"]["sectors"]
    )
    signature3 = tuple(
        (item["dimension"], mp.nstr(item["old_central_eigenvalue"], 70), item["splitter"])
        for item in slab_data["slab3"]["sectors"]
    )
    all_carrier_order &= signature2 == signature3

    for sector_index in TARGET_SECTORS:
        sector2 = slab_data["slab2"]["sectors"][sector_index]
        sector3 = slab_data["slab3"]["sectors"][sector_index]
        dimension = sector2["dimension"]
        if dimension != 1 or sector3["dimension"] != dimension:
            all_carrier_order = False
        n = 30 * dimension
        r = 5 * dimension
        s = 25 * dimension
        basis = mp_to_numpy(sector2["basis"])
        edge_basis = np.kron(np.eye(30, dtype=np.complex128), basis)
        incidence = incidences[parity]["incidence"].astype(np.complex128)
        compressed_incidence = edge_basis.conj().T @ incidence
        left_incidence, singular_incidence, _ = la.svd(
            compressed_incidence, full_matrices=False
        )
        epsilon_c = float(
            1000 * MACHINE_EPSILON * max(compressed_incidence.shape)
            * max(1.0, float(singular_incidence[0]))
        )
        incidence_rank = int(np.sum(singular_incidence > 100 * epsilon_c))
        incidence_open = int(np.sum(
            (singular_incidence >= 10 * epsilon_c)
            & (singular_incidence <= 100 * epsilon_c)
        ))
        gap_c = float(singular_incidence[r - 1])
        eta_k = (
            float(2 * epsilon_c / (gap_c - 2 * epsilon_c)
                  + 1000 * MACHINE_EPSILON * n)
            if gap_c > 2 * epsilon_c else math.inf
        )
        u_basis = left_incidence[:, :r]

        projected = {"slab2": {}, "slab3": {}}
        for slab_name, sector in (("slab2", sector2), ("slab3", sector3)):
            for variant in VARIANTS:
                block = project_full_kernel(
                    slab_data[slab_name]["kernels"][variant], sector
                )
                _, determinant, tangent, _ = build_tangent_ball(
                    block, dimension, slab_data[slab_name]["mapping"]
                )
                all_twists &= not determinant.contains(0)
                tangent_blocks = split_tangent(tangent, n)
                principal, residuals = reconstruct_principal(tangent_blocks)
                all_principal_identities &= all(
                    all_entries_contain_zero(value) for value in residuals.values()
                )
                projected[slab_name][variant] = principal

        for variant in VARIANTS:
            principal2 = projected["slab2"][variant]
            principal3 = projected["slab3"][variant]
            kminus = principal2["10"]
            kzero = principal2["11"] + principal3["00"]
            kplus = principal3["01"]
            m_ball = scale(kminus + kplus, "0.5")
            v_ball = kminus + kzero + kplus
            mass_regular = not m_ball.det().contains(0)
            omega_ball = m_ball.solve(v_ball)
            prefix = f"{parity}_sector{sector_index}_{variant}"
            m_mid, m_rad, m_overlap, m_ratio, old_m_rad, new_m_rad = (
                direct_and_old_matrix(m_ball, old_archive, prefix, "M")
            )
            v_mid, v_rad, v_overlap, v_ratio, old_v_rad, new_v_rad = (
                direct_and_old_matrix(v_ball, old_archive, prefix, "V")
            )
            omega_mid, omega_stored = acb_midpoint_and_radii(omega_ball)
            omega_rad = component_reenclosure_radii(omega_mid, omega_stored)
            overlap_ok = m_overlap and v_overlap
            reduction_ok = v_ratio > 100
            all_overlap &= overlap_ok
            all_precision_reduction &= reduction_ok

            h_m = (m_mid + m_mid.conj().T) / 2
            h_v = (v_mid + v_mid.conj().T) / 2
            radius_hm = (m_rad + m_rad.T) / 2
            radius_hv = (v_rad + v_rad.T) / 2
            epsilon_hm = matrix_error(h_m, radius_hm, n)
            epsilon_hv = matrix_error(h_v, radius_hv, n)
            epsilon_omega = matrix_error(omega_mid, omega_rad, n)

            row = u_basis.conj().T @ h_m
            _, singular_row, right_row = la.svd(row, full_matrices=True)
            norm_row = operator_norm(row)
            epsilon_row = float(
                epsilon_hm
                + 2 * eta_k * (operator_norm(h_m) + epsilon_hm)
                + 1000 * MACHINE_EPSILON * n * max(1.0, norm_row)
            )
            row_rank = int(np.sum(singular_row > 100 * epsilon_row))
            row_zero = int(np.sum(singular_row < 10 * epsilon_row))
            row_open = len(singular_row) - row_rank - row_zero
            row_gap = float(singular_row[r - 1])
            eta_s = (
                float(2 * epsilon_row / (row_gap - 2 * epsilon_row)
                      + 1000 * MACHINE_EPSILON * n)
                if row_gap > 2 * epsilon_row else math.inf
            )
            w_basis = right_row.conj().T[:, r:]
            direct_sum = np.column_stack((u_basis, w_basis))
            direct_singular = la.svdvals(direct_sum)
            direct_floor = float(
                1000 * MACHINE_EPSILON * n
                * max(1.0, float(direct_singular[0]))
            )
            direct_sum_ok = direct_singular[-1] > 100 * direct_floor
            shape_ok = bool(
                incidence_rank == r
                and incidence_open == 0
                and row_rank == r
                and row_zero == 0
                and row_open == 0
                and math.isfinite(eta_s)
                and direct_sum_ok
            )
            all_shape_carriers &= shape_ok

            m_s = w_basis.conj().T @ h_m @ w_basis
            m_s = (m_s + m_s.conj().T) / 2
            v_s = w_basis.conj().T @ h_v @ w_basis
            v_s = (v_s + v_s.conj().T) / 2
            omega_s = w_basis.conj().T @ omega_mid @ w_basis
            epsilon_ms = restriction_error(h_m, epsilon_hm, eta_s, n)
            epsilon_vs = restriction_error(h_v, epsilon_hv, eta_s, n)
            epsilon_os = restriction_error(omega_mid, epsilon_omega, eta_s, n)

            b_values = la.eigvalsh(-m_s)
            a_values = la.eigvalsh(-v_s)
            kinetic_ok = bool(
                mass_regular and float(b_values[0]) > 100 * epsilon_ms
            )
            all_kinetic &= kinetic_ok
            labels = [sign_label(float(value), epsilon_vs) for value in a_values]
            counts = Counter(labels)
            aggregate_signs.update(counts)

            residual = m_s @ omega_s - v_s
            residual_norm = operator_norm(residual)
            residual_error = float(
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
            compatibility = scalar_zero_label(residual_norm, residual_error)
            compatibility_ok = compatibility == "ZERO_CONSISTENT"
            all_compatibility &= compatibility_ok

            records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "position_dimension": n,
                "shape_dimension": s,
                "controls": {
                    "incidence_rank": incidence_rank,
                    "incidence_open": incidence_open,
                    "shape_row_rank": row_rank,
                    "shape_row_zero": row_zero,
                    "shape_row_open": row_open,
                    "shape_subspace_error": sf(eta_s),
                    "direct_sum_minimum_singular": sf(direct_singular[-1]),
                    "direct_sum_floor": sf(direct_floor),
                    "old_new_M_balls_overlap": m_overlap,
                    "old_new_V_balls_overlap": v_overlap,
                    "old_V_radius_frobenius": sf(old_v_rad),
                    "direct_V_radius_frobenius": sf(new_v_rad),
                    "V_precision_reduction_ratio": sf(v_ratio),
                    "M_precision_reduction_ratio": sf(m_ratio),
                    "mass_determinant_excludes_zero": mass_regular,
                    "kinetic_minimum_eigenvalue": sf(b_values[0]),
                    "kinetic_error": sf(epsilon_ms),
                    "kinetic_positive_definite_resolved": kinetic_ok,
                    "compatibility_residual": sf(residual_norm),
                    "compatibility_error": sf(residual_error),
                    "compatibility_label": compatibility,
                },
                "stiffness": {
                    "minimum_eigenvalue": sf(a_values[0]),
                    "maximum_eigenvalue": sf(a_values[-1]),
                    "restricted_error": sf(epsilon_vs),
                    "minimum_absolute_error_units": sf(
                        np.min(np.abs(a_values)) / epsilon_vs
                        if epsilon_vs else math.inf
                    ),
                    "sign_counts": dict(counts),
                    "labels": labels,
                    "eigenvalues": [sf(value) for value in a_values],
                },
            })

check("all four direct slab branch reconstructions pass", all_branch_controls)
check("both slabs retain the same ordered sector carrier", all_carrier_order)
check("all 32 direct boundary-twist determinant balls exclude zero", all_twists)
check("all 32 principal-function identity families contain zero entrywise",
      all_principal_identities)
check("all 16 direct M,V balls overlap the committed broad controls", all_overlap)
reduction_ratios = [float(item["controls"]["V_precision_reduction_ratio"])
                    for item in records]
check("all 16 final V radii improve by a factor greater than 100",
      all_precision_reduction,
      f"ratio={min(reduction_ratios):.6e}...{max(reduction_ratios):.6e}")
check("all 16 action-selected shape carriers are resolved", all_shape_carriers)
check("all 16 restricted kinetic forms are positive-definite resolved", all_kinetic)
check("all 16 action-compatibility residuals are zero-consistent", all_compatibility)

complete_cells = [
    record for record in records
    if record["stiffness"]["sign_counts"].get("OPEN", 0) == 0
    and record["stiffness"]["sign_counts"].get("ZERO_CONSISTENT", 0) == 0
]
changed_cells = [
    record for record in complete_cells
    if record["stiffness"]["sign_counts"].get("NEGATIVE_RESOLVED", 0) != 15
    or record["stiffness"]["sign_counts"].get("POSITIVE_RESOLVED", 0) != 10
]
target_cells = [
    record for record in complete_cells
    if record["stiffness"]["sign_counts"].get("NEGATIVE_RESOLVED", 0) == 15
    and record["stiffness"]["sign_counts"].get("POSITIVE_RESOLVED", 0) == 10
]
controls_ok = bool(
    provenance_ok and geometry_import_ok and incidence_ok
    and all_branch_controls and all_carrier_order and all_twists
    and all_principal_identities and all_overlap and all_precision_reduction
    and all_shape_carriers and all_kinetic and all_compatibility
    and len(records) == 16
)
if not controls_ok:
    outcome = "SHIFTED_DIRECT_PRECISION_CONTROL_FAILED"
elif changed_cells:
    outcome = "SHIFTED_DIRECT_NEGATIVE_RANK_CHANGED"
elif len(target_cells) < 16:
    outcome = "SHIFTED_DIRECT_NEGATIVE_RANK_OPEN"
else:
    outcome = "SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS"

allowed = {
    "SHIFTED_DIRECT_PRECISION_CONTROL_FAILED",
    "SHIFTED_DIRECT_NEGATIVE_RANK_CHANGED",
    "SHIFTED_DIRECT_NEGATIVE_RANK_OPEN",
    "SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS",
}
check("the preregistered direct-precision hierarchy assigns the outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "precision_attribution_commit": ATTRIBUTION_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "old_stiffness_eigenvalues_loaded": False,
    "old_stiffness_sign_labels_loaded": False,
    "target_sectors": list(TARGET_SECTORS),
    "cells": len(records),
    "outcome": outcome,
    "controls_ok": controls_ok,
    "complete_cells": len(complete_cells),
    "target_rank_cells": len(target_cells),
    "changed_cells": len(changed_cells),
    "aggregate_sign_counts": dict(aggregate_signs),
    "minimum_V_precision_reduction_ratio": sf(min(reduction_ratios)),
    "maximum_V_precision_reduction_ratio": sf(max(reduction_ratios)),
    "records": records,
    "classification": {
        "binary_serialization_attribution": "DERIVED COMPUTATIONAL",
        "shifted_negative_rank": (
            "DERIVED COMPUTATIONAL" if outcome.endswith("PERSISTS") else "OPEN"
        ),
        "common_projector_or_connection": "OPEN",
        "physical_inertia_or_instability": "OPEN",
        "wave_graviton_continuum_or_speed": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"complete cells: {len(complete_cells)}/16")
print(f"target 15+10 cells: {len(target_cells)}/16")
print(f"aggregate signs: {dict(aggregate_signs)}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

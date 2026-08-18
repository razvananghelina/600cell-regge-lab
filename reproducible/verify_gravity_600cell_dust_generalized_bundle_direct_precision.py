#!/usr/bin/env python3
"""Direct-precision generalized-mode bundle and closure audit.

Prior-art/framing commit: 6fdddc6.
Preregistered protocol commit: d177762.
"""

import ast
from collections import Counter, defaultdict
from itertools import combinations
import contextlib
import hashlib
import importlib.util
import io
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
OLD_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
OLD_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
NEW_JSON = HERE / "gravity_600cell_dust_shifted_centered.json"
NEW_NPZ = HERE / "gravity_600cell_dust_shifted_centered.npz"
DIRECT_RANK = HERE / "gravity_600cell_dust_shifted_direct_precision.json"
BROAD_GENERALIZED = HERE / "gravity_600cell_dust_generalized_mode_closure.json"
DIRECT_SOURCE = HERE / "verify_gravity_600cell_dust_shifted_direct_precision.py"
GENERALIZED_SOURCE = HERE / "verify_gravity_600cell_dust_generalized_mode_closure.py"
OUTPUT = HERE / "gravity_600cell_dust_generalized_bundle_direct_precision.json"

PRIOR_ART_COMMIT = "6fdddc6"
PROTOCOL_COMMIT = "d177762"
EXPECTED_HASHES = {
    "first_tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "second_tick": "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70",
    "third_tick": "ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "old_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "old_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "new_json": "265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47",
    "new_npz": "c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8",
    "direct_rank": "86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944",
    "broad_generalized": "53e046e2020a97fc992559546ce3d45479c0c0de7ce2e01322b09998ba85cf80",
    "direct_source": "1b54cd25899037fc66c2b58e01ef3bac267c6ebf2c6917d2a05ac4ac0feed1c5",
    "generalized_source": "0a84c8ec4fab1c9626d5e4c711f89c6f9638cf37c15ef6b0050d6b66dfdde6c1",
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
TIMES = ("old", "shifted")
TARGET_SECTORS = (4, 5)
MATRIX_NAMES = ("M", "N", "V", "Gamma", "Omega")
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


def operator_norm(matrix):
    values = la.svdvals(matrix)
    return float(values[0]) if len(values) else 0.0


def restriction_error(midpoint, base_error, eta_shape, n):
    norm = operator_norm(midpoint)
    return float(
        base_error
        + 2 * eta_shape * (norm + base_error)
        + 1000 * MACHINE_EPSILON * n * max(1.0, norm)
    )


def projector_identity_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error):
        return "GENERALIZED_FIBER_IDENTITY_OPEN"
    if value <= 10 * error:
        return "GENERALIZED_COMMON_FIBER_RESOLVED"
    if value > 100 * error:
        return "GENERALIZED_ROTATED_FIBER_RESOLVED"
    return "GENERALIZED_FIBER_IDENTITY_OPEN"


def leakage_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error):
        return "LEAKAGE_OPEN"
    if value <= 10 * error:
        return "LEAKAGE_ZERO_CONSISTENT"
    if value > 100 * error:
        return "LEAKAGE_NONZERO_RESOLVED"
    return "LEAKAGE_OPEN"


def source_ball(matrix):
    midpoint, stored = acb_midpoint_and_radii(matrix)
    radii = component_reenclosure_radii(midpoint, stored)
    return midpoint, radii


def broad_ball(archive, prefix, name):
    midpoint = np.asarray(archive[f"{prefix}_{name}_midpoint"])
    stored = np.asarray(archive[f"{prefix}_{name}_radii"])
    return midpoint, component_reenclosure_radii(midpoint, stored)


def balls_overlap(left_mid, left_rad, right_mid, right_rad):
    tolerance = left_rad + right_rad
    difference = left_mid - right_mid
    return bool(
        np.all(np.abs(np.real(difference)) <= tolerance)
        and np.all(np.abs(np.imag(difference)) <= tolerance)
    )


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


paths = {
    "first_tick": FIRST_TICK,
    "second_tick": SECOND_TICK,
    "third_tick": THIRD_TICK,
    "geometry_source": GEOMETRY_SOURCE,
    "rank_source": RANK_SOURCE,
    "full_source": FULL_SOURCE,
    "conformal_source": CONFORMAL_SOURCE,
    "commons": COMMONS_SOURCE,
    "old_json": OLD_JSON,
    "old_npz": OLD_NPZ,
    "new_json": NEW_JSON,
    "new_npz": NEW_NPZ,
    "direct_rank": DIRECT_RANK,
    "broad_generalized": BROAD_GENERALIZED,
    "direct_source": DIRECT_SOURCE,
    "generalized_source": GENERALIZED_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
first_tick = json.loads(FIRST_TICK.read_text())
second_tick = json.loads(SECOND_TICK.read_text())
third_tick = json.loads(THIRD_TICK.read_text())
old_json = json.loads(OLD_JSON.read_text())
new_json = json.loads(NEW_JSON.read_text())
direct_rank = json.loads(DIRECT_RANK.read_text())
broad_generalized = json.loads(BROAD_GENERALIZED.read_text())
old_archive = np.load(OLD_NPZ, allow_pickle=False)
new_archive = np.load(NEW_NPZ, allow_pickle=False)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and first_tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and second_tick["outcome"] == "SECOND_HOMOTHETIC_TICK_ACCEPTED"
    and third_tick["outcome"] == "THIRD_HOMOTHETIC_TICK_ACCEPTED"
    and second_tick["fixed_mass"] is True
    and third_tick["fixed_mass"] is True
    and old_json["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and new_json["outcome"] == "SHIFTED_CENTERED_CERTIFIED"
    and direct_rank["outcome"] == "SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS"
    and broad_generalized["outcome"]
    == "GENERALIZED_MODE_RECURRENCE_CLOSURE_CERTIFIED"
    and len(old_archive.files) == len(new_archive.files) == 560
)
check("all direct-bundle inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_direct_generalized_bundle", GEOMETRY_SOURCE
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
sector_data, sector_controls = high_precision_sector_bases(groups["even"])
sector_ok = bool(
    tuple(item["dimension"] for item in sector_data) == (3, 2, 2, 2, 1, 1, 1)
    and sector_controls["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
)
check("the exact incidence and deterministic sector carriers reconstruct",
      incidence_ok and sector_ok)

print("=" * 78)
print("DIRECT-PRECISION GENERALIZED-MODE BUNDLE AUDIT")
print("=" * 78)

matrix_cells = {}
matrix_records = []
hessian_records = []
principal_records = []
all_branch_controls = True
all_carrier_order = True
all_hessian_family = True
all_twists = True
all_principal = True
all_overlaps = True
all_mass_regular = True
v_reductions = []

for parity in PARITIES:
    a1, r1 = [mp.mpf(value) for value in first_tick["solutions"][parity]["state"]]
    a2, r2 = [
        mp.mpf(value)
        for value in second_tick["solutions"][parity]["state_absolute"]
    ]
    a3, r3 = [
        mp.mpf(value)
        for value in third_tick["solutions"][parity]["state_absolute"]
    ]
    slab_specs = (
        ("slab1", mp.mpf(0), a1, r1),
        ("slab2", a1, a2, r2),
        ("slab3", a2, a3, r3),
    )
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
        carrier_ok = bool(
            edge_order == groups[parity]["edge_order"]
            and sorted(mapping) == list(range(30))
            and tuple(item["dimension"] for item in sectors)
            == (3, 2, 2, 2, 1, 1, 1)
            and sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        )
        all_carrier_order &= carrier_ok
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
        }

    signatures = [
        tuple(
            (item["dimension"], mp.nstr(item["old_central_eigenvalue"], 70),
             item["splitter"])
            for item in slab_data[name]["sectors"]
        )
        for name in ("slab1", "slab2", "slab3")
    ]
    all_carrier_order &= signatures[0] == signatures[1] == signatures[2]

    for sector_index in TARGET_SECTORS:
        projected = {name: {} for name in ("slab1", "slab2", "slab3")}
        for slab_name in projected:
            sector = slab_data[slab_name]["sectors"][sector_index]
            dimension = sector["dimension"]
            all_carrier_order &= dimension == 1
            raw_blocks = {
                variant: project_full_kernel(
                    slab_data[slab_name]["kernels"][variant], sector
                )
                for variant in VARIANTS
            }
            primary_raw = raw_blocks[VARIANTS[0]]
            family_variation = max(
                mp_frobenius(raw_blocks[variant] - primary_raw)
                for variant in VARIANTS
            )
            for variant in VARIANTS:
                raw = raw_blocks[variant]
                defect = mp_frobenius((raw - raw.H) / 2)
                enveloped = bool(defect <= family_variation + ARITHMETIC_FLOOR)
                all_hessian_family &= enveloped
                block = (raw + raw.H) / 2
                _, determinant, tangent, _ = build_tangent_ball(
                    block, dimension, slab_data[slab_name]["mapping"]
                )
                twist_ok = not determinant.contains(0)
                all_twists &= twist_ok
                principal, residuals = reconstruct_principal(
                    split_tangent(tangent, 30)
                )
                failing = {
                    name: sum(
                        not residual[row, column].contains(0)
                        for row in range(residual.nrows())
                        for column in range(residual.ncols())
                    )
                    for name, residual in residuals.items()
                }
                principal_ok = not any(failing.values())
                all_principal &= principal_ok
                projected[slab_name][variant] = principal
                hessian_records.append({
                    "parity": parity,
                    "slab": slab_name,
                    "sector_index": sector_index,
                    "variant": variant,
                    "antihermitian_defect": mp.nstr(defect, 25),
                    "family_variation": mp.nstr(family_variation, 25),
                    "enveloped": enveloped,
                    "twist_regular": twist_ok,
                })
                principal_records.append({
                    "parity": parity,
                    "slab": slab_name,
                    "sector_index": sector_index,
                    "variant": variant,
                    "failing_entries": failing,
                })

        for time_name, lower, upper, archive in (
            ("old", "slab1", "slab2", old_archive),
            ("shifted", "slab2", "slab3", new_archive),
        ):
            for variant in VARIANTS:
                lower_principal = projected[lower][variant]
                upper_principal = projected[upper][variant]
                kminus = lower_principal["10"]
                kzero = lower_principal["11"] + upper_principal["00"]
                kplus = upper_principal["01"]
                balls = {
                    "M": scale(kminus + kplus, "0.5"),
                    "N": scale(kplus - kminus, "0.5"),
                    "V": kminus + kzero + kplus,
                }
                all_mass_regular &= not balls["M"].det().contains(0)
                balls["Gamma"] = balls["M"].solve(balls["N"])
                balls["Omega"] = balls["M"].solve(balls["V"])
                prefix = f"{parity}_sector{sector_index}_{variant}"
                matrices = {}
                for name in MATRIX_NAMES:
                    midpoint, radii = source_ball(balls[name])
                    old_midpoint, old_radii = broad_ball(
                        archive, prefix, name
                    )
                    overlap = balls_overlap(
                        midpoint, radii, old_midpoint, old_radii
                    )
                    direct_radius = float(la.norm(radii, "fro"))
                    broad_radius = float(la.norm(old_radii, "fro"))
                    reduction = (
                        broad_radius / direct_radius
                        if direct_radius > 0 else math.inf
                    )
                    all_overlaps &= overlap
                    if name == "V":
                        v_reductions.append(reduction)
                    matrices[name] = {"midpoint": midpoint, "radii": radii}
                    matrix_records.append({
                        "time": time_name,
                        "parity": parity,
                        "sector_index": sector_index,
                        "variant": variant,
                        "matrix": name,
                        "overlap": overlap,
                        "direct_radius_frobenius": sf(direct_radius),
                        "broad_radius_frobenius": sf(broad_radius),
                        "radius_reduction": sf(reduction),
                    })
                matrix_cells[(time_name, parity, sector_index, variant)] = matrices

check("all six direct slab branch reconstructions pass", all_branch_controls)
check("all slabs retain one common ordered sector carrier", all_carrier_order)
check("all 48 Hessian defects lie inside their derivative-family variation",
      all_hessian_family)
check("all 48 twists and principal identity families pass",
      all_twists and all_principal)
check("all 160 direct matrices overlap their broader serialized controls",
      all_overlaps, str(Counter(item["matrix"] for item in matrix_records
                                if not item["overlap"])))
v_precision_ok = len(v_reductions) == 32 and min(v_reductions) > 100
check("all 32 direct V radii improve by more than the inherited factor 100",
      v_precision_ok,
      f"ratio={min(v_reductions):.6e}...{max(v_reductions):.6e}")

# Complete four-variant finite-difference family envelopes.
for time_name in TIMES:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            for variant in VARIANTS:
                cell = matrix_cells[(time_name, parity, sector_index, variant)]
                for name in MATRIX_NAMES:
                    midpoint = cell[name]["midpoint"]
                    radii = cell[name]["radii"]
                    radius_norm = float(la.norm(radii, "fro"))
                    candidates = []
                    for other in VARIANTS:
                        other_cell = matrix_cells[
                            (time_name, parity, sector_index, other)
                        ][name]
                        candidates.append(
                            operator_norm(midpoint - other_cell["midpoint"])
                            + radius_norm
                            + float(la.norm(other_cell["radii"], "fro"))
                        )
                    arithmetic = float(
                        1000 * MACHINE_EPSILON * 30
                        * max(1.0, operator_norm(midpoint))
                    )
                    cell[name]["error"] = max(candidates) + arithmetic

# Deterministic conformal bases in each parity/sector.
conformal = {}
for parity in PARITIES:
    incidence = incidences[parity]["incidence"].astype(np.complex128)
    for sector_index in TARGET_SECTORS:
        basis = mp_to_numpy(sector_data[sector_index]["basis"])
        edge_basis = np.kron(np.eye(30, dtype=np.complex128), basis)
        compressed = edge_basis.conj().T @ incidence
        left, singular, _ = la.svd(compressed, full_matrices=False)
        epsilon_c = float(
            1000 * MACHINE_EPSILON * max(compressed.shape)
            * max(1.0, float(singular[0]))
        )
        rank_c = int(np.sum(singular > 100 * epsilon_c))
        gap_c = float(singular[4])
        eta_k = (
            float(2 * epsilon_c / (gap_c - 2 * epsilon_c)
                  + 1000 * MACHINE_EPSILON * 30)
            if gap_c > 2 * epsilon_c else math.inf
        )
        conformal[(parity, sector_index)] = {
            "basis": left[:, :5], "rank": rank_c, "eta": eta_k,
        }

projectors = {}
projector_records = []
all_projectors = True
for time_name in TIMES:
    for parity in PARITIES:
        for sector_index in TARGET_SECTORS:
            u_basis = conformal[(parity, sector_index)]["basis"]
            eta_k = conformal[(parity, sector_index)]["eta"]
            rank_c = conformal[(parity, sector_index)]["rank"]
            for variant in VARIANTS:
                cell = matrix_cells[(time_name, parity, sector_index, variant)]
                h_m = (cell["M"]["midpoint"]
                       + cell["M"]["midpoint"].conj().T) / 2
                h_v = (cell["V"]["midpoint"]
                       + cell["V"]["midpoint"].conj().T) / 2
                epsilon_hm = cell["M"]["error"]
                epsilon_hv = cell["V"]["error"]
                row = u_basis.conj().T @ h_m
                _, singular_row, right_row = la.svd(row, full_matrices=True)
                epsilon_row = float(
                    epsilon_hm
                    + 2 * eta_k * (operator_norm(h_m) + epsilon_hm)
                    + 1000 * MACHINE_EPSILON * 30
                    * max(1.0, operator_norm(row))
                )
                row_rank = int(np.sum(singular_row > 100 * epsilon_row))
                row_zero = int(np.sum(singular_row < 10 * epsilon_row))
                row_open = len(singular_row) - row_rank - row_zero
                row_gap = float(singular_row[4])
                eta_s = (
                    float(2 * epsilon_row / (row_gap - 2 * epsilon_row)
                          + 1000 * MACHINE_EPSILON * 30)
                    if row_gap > 2 * epsilon_row else math.inf
                )
                w_basis = right_row.conj().T[:, 5:]
                m_s = w_basis.conj().T @ h_m @ w_basis
                m_s = (m_s + m_s.conj().T) / 2
                v_s = w_basis.conj().T @ h_v @ w_basis
                v_s = (v_s + v_s.conj().T) / 2
                epsilon_ms = restriction_error(
                    h_m, epsilon_hm, eta_s, 30
                )
                epsilon_vs = restriction_error(
                    h_v, epsilon_hv, eta_s, 30
                )
                b_matrix = -m_s
                a_matrix = -v_s
                b_values = la.eigvalsh(b_matrix)
                minimum_b = float(b_values[0])
                maximum_b = float(b_values[-1])
                b_lower = minimum_b - epsilon_ms
                kinetic_ok = b_lower > 0 and minimum_b > 100 * epsilon_ms
                if kinetic_ok:
                    values, vectors = la.eigh(a_matrix, b_matrix)
                    epsilon_pencil = float(
                        epsilon_vs / b_lower
                        + operator_norm(a_matrix) * epsilon_ms
                        / (minimum_b * b_lower)
                        + 1000 * MACHINE_EPSILON * 25
                        * max(1.0, operator_norm(la.solve(b_matrix, a_matrix)))
                    )
                else:
                    values = np.full(25, np.nan)
                    vectors = np.full((25, 25), np.nan, dtype=np.complex128)
                    epsilon_pencil = math.inf
                negative = int(np.sum(values < 0))
                positive = int(np.sum(values > 0))
                gap = float(values[15] - values[14])
                eta_eig = (
                    float(2 * epsilon_pencil / (gap - 2 * epsilon_pencil)
                          + 1000 * MACHINE_EPSILON * 30)
                    if gap > 2 * epsilon_pencil else math.inf
                )
                eta_p = (
                    float(
                        2 * eta_s
                        + math.sqrt(maximum_b / b_lower) * eta_eig
                        + epsilon_ms / b_lower
                        + 1000 * MACHINE_EPSILON * 30
                    )
                    if kinetic_ok and math.isfinite(eta_eig) else math.inf
                )
                lifted = w_basis @ vectors[:, :15]
                q_basis, _ = la.qr(lifted, mode="economic")
                projector = q_basis @ q_basis.conj().T
                projector_residual = max(
                    operator_norm(projector - projector.conj().T),
                    operator_norm(projector @ projector - projector),
                )
                complete = bool(
                    rank_c == 5
                    and row_rank == 5 and row_zero == 0 and row_open == 0
                    and kinetic_ok and negative == 15 and positive == 10
                    and gap > 2 * epsilon_pencil and math.isfinite(eta_p)
                    and projector_residual
                    <= 1000 * MACHINE_EPSILON * 30
                )
                all_projectors &= complete
                key = (time_name, parity, sector_index, variant)
                projectors[key] = {"projector": projector, "eta": eta_p}
                projector_records.append({
                    "time": time_name,
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "negative": negative,
                    "positive": positive,
                    "kinetic_minimum": sf(minimum_b),
                    "kinetic_error": sf(epsilon_ms),
                    "generalized_gap": sf(gap),
                    "pencil_error": sf(epsilon_pencil),
                    "projector_error": sf(eta_p),
                    "complete": complete,
                })

check("all 32 direct Hermitian-definite projectors are resolved",
      all_projectors)

identity_records = []
identity_counts = Counter()
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for variant in VARIANTS:
            old = projectors[("old", parity, sector_index, variant)]
            new = projectors[("shifted", parity, sector_index, variant)]
            distance = operator_norm(new["projector"] - old["projector"])
            error = old["eta"] + new["eta"] + 1000 * MACHINE_EPSILON * 30
            label = projector_identity_label(distance, error)
            identity_counts[label] += 1
            identity_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "distance": sf(distance),
                "error": sf(error),
                "error_units": sf(distance / error if error else math.inf),
                "label": label,
            })

leakage_records = []
local_counts = Counter()
cross_counts = Counter()
all_finite = True
for projector_time in TIMES:
    for operator_time in TIMES:
        kind = "local" if projector_time == operator_time else "cross"
        for parity in PARITIES:
            for sector_index in TARGET_SECTORS:
                for variant in VARIANTS:
                    p_record = projectors[
                        (projector_time, parity, sector_index, variant)
                    ]
                    projector = p_record["projector"]
                    one = np.eye(30, dtype=np.complex128)
                    cell = matrix_cells[
                        (operator_time, parity, sector_index, variant)
                    ]
                    for operator_name in ("Gamma", "Omega"):
                        operator = cell[operator_name]["midpoint"]
                        epsilon_x = cell[operator_name]["error"]
                        residual = (one - projector) @ operator @ projector
                        residual_norm = operator_norm(residual)
                        residual_error = float(
                            epsilon_x
                            + (2 * p_record["eta"] + p_record["eta"]**2)
                            * (operator_norm(operator) + epsilon_x)
                            + 1000 * MACHINE_EPSILON * 30
                            * max(1.0, operator_norm(operator))
                        )
                        label = leakage_label(residual_norm, residual_error)
                        (local_counts if kind == "local" else cross_counts)[label] += 1
                        all_finite &= bool(
                            np.all(np.isfinite(operator))
                            and math.isfinite(residual_norm)
                            and math.isfinite(residual_error)
                        )
                        leakage_records.append({
                            "kind": kind,
                            "projector_time": projector_time,
                            "operator_time": operator_time,
                            "parity": parity,
                            "sector_index": sector_index,
                            "variant": variant,
                            "operator": operator_name,
                            "residual_norm": sf(residual_norm),
                            "residual_error": sf(residual_error),
                            "error_units": sf(
                                residual_norm / residual_error
                                if residual_error else math.inf
                            ),
                            "label": label,
                        })

identity_complete = len(identity_records) == 16
leakage_complete = (
    len(leakage_records) == 128
    and sum(local_counts.values()) == 64
    and sum(cross_counts.values()) == 64
)
check("all 16 direct old/shifted projector comparisons receive labels",
      identity_complete, str(dict(identity_counts)))
check("all 128 direct local/cross leakages receive labels",
      leakage_complete and all_finite,
      f"local={dict(local_counts)}, cross={dict(cross_counts)}")

controls_ok = bool(
    provenance_ok and geometry_import_ok and incidence_ok and sector_ok
    and all_branch_controls and all_carrier_order and all_hessian_family
    and all_twists and all_principal and all_overlaps and all_mass_regular
    and v_precision_ok and all_projectors and identity_complete
    and leakage_complete and all_finite
)
if not controls_ok:
    outcome = "DIRECT_GENERALIZED_BUNDLE_CONTROL_FAILED"
elif local_counts["LEAKAGE_NONZERO_RESOLVED"]:
    outcome = "DIRECT_GENERALIZED_LOCAL_CLOSURE_REFUTED"
elif local_counts["LEAKAGE_OPEN"]:
    outcome = "DIRECT_GENERALIZED_LOCAL_CLOSURE_OPEN"
elif len(identity_counts) != 1:
    outcome = "DIRECT_GENERALIZED_BUNDLE_NONUNIFORM"
elif identity_counts["GENERALIZED_ROTATED_FIBER_RESOLVED"] == 16:
    outcome = "DIRECT_GENERALIZED_BUNDLE_ROTATION_RESOLVED"
elif identity_counts["GENERALIZED_FIBER_IDENTITY_OPEN"] == 16:
    outcome = "DIRECT_GENERALIZED_BUNDLE_IDENTITY_OPEN"
elif cross_counts["LEAKAGE_NONZERO_RESOLVED"]:
    outcome = "DIRECT_GENERALIZED_COMMON_CROSS_CLOSURE_REFUTED"
elif cross_counts["LEAKAGE_OPEN"]:
    outcome = "DIRECT_GENERALIZED_COMMON_CROSS_CLOSURE_OPEN"
else:
    outcome = "DIRECT_GENERALIZED_COMMON_BUNDLE_RESOLVED"

allowed = {
    "DIRECT_GENERALIZED_BUNDLE_CONTROL_FAILED",
    "DIRECT_GENERALIZED_LOCAL_CLOSURE_REFUTED",
    "DIRECT_GENERALIZED_LOCAL_CLOSURE_OPEN",
    "DIRECT_GENERALIZED_BUNDLE_NONUNIFORM",
    "DIRECT_GENERALIZED_BUNDLE_ROTATION_RESOLVED",
    "DIRECT_GENERALIZED_BUNDLE_IDENTITY_OPEN",
    "DIRECT_GENERALIZED_COMMON_CROSS_CLOSURE_REFUTED",
    "DIRECT_GENERALIZED_COMMON_CROSS_CLOSURE_OPEN",
    "DIRECT_GENERALIZED_COMMON_BUNDLE_RESOLVED",
}
check("the preregistered direct-bundle hierarchy assigns the outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "matrix_records": matrix_records,
    "hessian_records": hessian_records,
    "principal_records": principal_records,
    "projector_records": projector_records,
    "identity_counts": dict(identity_counts),
    "identity_records": identity_records,
    "local_leakage_counts": dict(local_counts),
    "cross_leakage_counts": dict(cross_counts),
    "leakage_records": leakage_records,
    "classification": {
        "direct_generalized_bundle": "DERIVED COMPUTATIONAL",
        "bundle_transport": "OPEN",
        "reduced_propagator": "NOT COMPUTED",
        "dispersion_inertia_mass_or_speed": "NOT COMPUTED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"identity labels: {dict(identity_counts)}")
print(f"local leakage labels: {dict(local_counts)}")
print(f"cross leakage labels: {dict(cross_counts)}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

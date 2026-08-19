#!/usr/bin/env python3
"""Direct adversarial audit of shifted pseudo-longitudinal persistence."""

import ast
from collections import Counter, defaultdict
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations, permutations, product
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
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402

FIRST_TICK = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
SECOND_TICK = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
THIRD_TICK = HERE / "gravity_600cell_dust_third_tick_local_correction.json"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_pseudolongitudinal_shifted_adversarial_protocol.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_dust_pseudolongitudinal_shifted_prior_art.md"
PRIMARY_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_pseudolongitudinal_shifted_protocol.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_dust_pseudolongitudinal_shifted.py"
PRIMARY_JSON = HERE / "gravity_600cell_dust_pseudolongitudinal_shifted.json"
SHIFTED_DIRECT_SOURCE = HERE / "verify_gravity_600cell_dust_shifted_direct_precision.py"
SHIFTED_DIRECT_JSON = HERE / "gravity_600cell_dust_shifted_direct_precision.json"
CURRENT_DIRECT_SOURCE = HERE / "verify_gravity_600cell_dust_action_york_direct_precision.py"
OUTPUT = HERE / "gravity_600cell_dust_pseudolongitudinal_shifted_adversarial.json"

PROTOCOL_COMMIT = "1bfe9e9"
EXPECTED_HASHES = {
    "protocol": "1dc9712a46b6ff6ac3c9b62e9d144f959f85622fe1ddbe2fc84de6ece3fa0982",
    "prior_art": "740eefaee14ea3ff634f8cff237041cecd675c4ceaf3d5be6ccdb9a3778a57ef",
    "primary_protocol": "f09c7450f238774b48425674e5b0800d6ef8ea2fee9f45b82794a0edef0a2375",
    "primary_source": "e4c5bcc18007c1c0ba7fbd38e29dffcc33a526fd790dbfcba8defe2ae44b7ab2",
    "primary_json": "0480f5d49d24e0f5d8e4e95f0cf62b7d0d9242459ed2b8f6d8e835ecd6e103a7",
    "shifted_direct_source": "1b54cd25899037fc66c2b58e01ef3bac267c6ebf2c6917d2a05ac4ac0feed1c5",
    "shifted_direct_json": "86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944",
    "current_direct_source": "73d852d58b21a9a15306a565d5cf4fb998b159fadb82830739ab0996ac07270e",
    "first_tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "second_tick": "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70",
    "third_tick": "ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
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
    if abs(value) <= 10 * error:
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


def exact_vertices():
    phi = (1 + mp.sqrt(5)) / 2
    candidates = {}

    def add(values):
        norm = mp.sqrt(mp.fsum(value * value for value in values))
        exact = tuple(value / norm for value in values)
        key = tuple(round(float(value), 10) for value in exact)
        candidates[key] = exact

    for coordinate in range(4):
        for sign in (1, -1):
            values = [mp.mpf(0)] * 4
            values[coordinate] = mp.mpf(sign)
            add(values)
    for signs in product((1, -1), repeat=4):
        add([mp.mpf(sign) / 2 for sign in signs])
    base = (mp.mpf(0), mp.mpf(1) / 2, phi / 2, 1 / (2 * phi))
    even = []
    for order in permutations(range(4)):
        inversions = sum(
            order[left] > order[right]
            for left in range(4) for right in range(left + 1, 4)
        )
        if inversions % 2 == 0:
            even.append(order)
    for order in even:
        values = [base[order[index]] for index in range(4)]
        nonzero = [index for index, value in enumerate(values) if value]
        for signs in product((1, -1), repeat=3):
            signed = list(values)
            for index, sign in zip(nonzero, signs):
                signed[index] *= sign
            add(signed)
    keys = sorted(candidates)
    if len(keys) != 120:
        raise RuntimeError(f"exact 600-cell has {len(keys)} vertices")
    return [candidates[key] for key in keys], np.asarray(keys, dtype=float)


def compressed_geometry(index_data, sector, vertices):
    basis = sector["basis"]
    if sector["dimension"] != 1 or basis.rows != 24:
        raise RuntimeError("selected carrier is not a one-dimensional sector")
    phi = (1 + mp.sqrt(5)) / 2
    length_square = 2 - phi
    incidence = mp.matrix(30, 120)
    tangent = mp.matrix(30, 480)
    edge_order = []
    for orbit in range(30):
        orbit_edges = index_data["orbit_edges"][orbit]
        if len(orbit_edges) != 24:
            raise RuntimeError("old edge orbit is not free")
        for group, edge in enumerate(orbit_edges):
            left, right = edge
            edge_order.append(tuple(sorted((left, right))))
            weight = mp.conj(basis[group, 0])
            incidence[orbit, left] += weight
            incidence[orbit, right] += weight
            delta = [vertices[left][axis] - vertices[right][axis]
                     for axis in range(4)]
            coefficient = [2 * value / length_square for value in delta]
            for axis in range(4):
                tangent[orbit, 4 * left + axis] += weight * (
                    coefficient[axis] - vertices[left][axis]
                )
                tangent[orbit, 4 * right + axis] += weight * (
                    -coefficient[axis] - vertices[right][axis]
                )
    return incidence, tangent, tuple(edge_order)


def direct_relative_cell(m_mid, v_mid, m_radii, v_radii, c_matrix, d_matrix):
    h_m = (m_mid + m_mid.conj().T) / 2
    h_v = (v_mid + v_mid.conj().T) / 2
    left_c, singular_c, _ = la.svd(c_matrix, full_matrices=False)
    u_basis = left_c[:, :5]
    row = u_basis.conj().T @ h_m
    _, singular_row, right_row = la.svd(row, full_matrices=True)
    w_basis = right_row.conj().T[:, 5:]
    q_basis = np.column_stack((u_basis, w_basis))
    coefficients = la.solve(q_basis, d_matrix)
    longitudinal_raw = coefficients[5:, :]
    left_l, singular_l, _ = la.svd(longitudinal_raw, full_matrices=False)
    l_basis = left_l[:, :15]

    b_matrix = -(w_basis.conj().T @ h_m @ w_basis)
    a_matrix = -(w_basis.conj().T @ h_v @ w_basis)
    b_matrix = (b_matrix + b_matrix.conj().T) / 2
    a_matrix = (a_matrix + a_matrix.conj().T) / 2
    bl = b_matrix @ l_basis
    al = a_matrix @ l_basis
    left_bl, singular_bl, _ = la.svd(bl, full_matrices=False)
    generalized = la.solve(b_matrix, a_matrix)
    gl = generalized @ l_basis

    span_residual = operator_norm(
        (np.eye(25) - left_bl @ left_bl.conj().T) @ al
    )
    commutator_residual = operator_norm(
        (np.eye(25) - l_basis @ l_basis.conj().T) @ gl
    )
    al_norm = operator_norm(al)
    gl_norm = operator_norm(gl)
    relative_span = span_residual / al_norm if al_norm > 0 else math.inf
    relative_commutator = (
        commutator_residual / gl_norm if gl_norm > 0 else math.inf
    )

    b_values = la.eigvalsh(b_matrix)
    generalized_values = la.eigvalsh(a_matrix, b_matrix)
    norm_a = operator_norm(a_matrix)
    norm_g = operator_norm(generalized)
    kappa = max(
        float(np.linalg.cond(q_basis)),
        float(singular_l[0] / singular_l[14]),
        float(np.linalg.cond(b_matrix)),
        float(singular_bl[0] / singular_bl[-1]),
        float(norm_a / al_norm) if al_norm else math.inf,
        float(norm_g / gl_norm) if gl_norm else math.inf,
    )
    epsilon_m = matrix_error(h_m, (m_radii + m_radii.T) / 2, 30)
    epsilon_v = matrix_error(h_v, (v_radii + v_radii.T) / 2, 30)
    relative_floor = float(
        1000 * MACHINE_EPSILON * 30 * kappa
        + kappa * max(
            epsilon_m / max(1.0, operator_norm(h_m)),
            epsilon_v / max(1.0, operator_norm(h_v)),
        )
    )
    augmented = la.svdvals(np.column_stack((bl, al)))

    norm_b = operator_norm(b_matrix)
    inverse_norm_b = 1 / float(la.svdvals(b_matrix)[-1])
    inequality_floor = float(
        relative_floor * max(1.0, span_residual, commutator_residual)
    )
    inequality_one = bool(
        span_residual <= norm_b * commutator_residual + inequality_floor
    )
    inequality_two = bool(
        commutator_residual <= inverse_norm_b * span_residual + inequality_floor
    )
    return {
        "rank_C": int(np.sum(singular_c > 1e-6)),
        "rank_shape_row": int(np.sum(singular_row > 1e-6)),
        "rank_longitudinal": int(np.sum(singular_l > 1e-6)),
        "minimum_C_gap": float(singular_c[4]),
        "minimum_shape_gap": float(singular_row[4]),
        "minimum_longitudinal_gap": float(singular_l[14]),
        "minimum_B_eigenvalue": float(b_values[0]),
        "negative_stiffness": int(np.sum(generalized_values < 0)),
        "positive_stiffness": int(np.sum(generalized_values > 0)),
        "span_residual": span_residual,
        "commutator_residual": commutator_residual,
        "AL_norm": al_norm,
        "GL_norm": gl_norm,
        "relative_span": relative_span,
        "relative_commutator": relative_commutator,
        "kappa_relative": kappa,
        "relative_floor": relative_floor,
        "augmented_singular": augmented,
        "inequality_one": inequality_one,
        "inequality_two": inequality_two,
    }


def direct_matrix(matrix):
    midpoint, stored = acb_midpoint_and_radii(matrix)
    radii = component_reenclosure_radii(midpoint, stored)
    return midpoint, radii


paths = {
    "protocol": PROTOCOL,
    "prior_art": PRIOR_ART,
    "primary_protocol": PRIMARY_PROTOCOL,
    "primary_source": PRIMARY_SOURCE,
    "shifted_direct_source": SHIFTED_DIRECT_SOURCE,
    "shifted_direct_json": SHIFTED_DIRECT_JSON,
    "current_direct_source": CURRENT_DIRECT_SOURCE,
    "first_tick": FIRST_TICK,
    "second_tick": SECOND_TICK,
    "third_tick": THIRD_TICK,
    "geometry_source": GEOMETRY_SOURCE,
    "rank_source": RANK_SOURCE,
    "full_source": FULL_SOURCE,
    "conformal_source": CONFORMAL_SOURCE,
    "commons": COMMONS_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
first_tick = json.loads(FIRST_TICK.read_text())
second_tick = json.loads(SECOND_TICK.read_text())
third_tick = json.loads(THIRD_TICK.read_text())
shifted_direct = json.loads(SHIFTED_DIRECT_JSON.read_text())
provenance_ok = bool(
    hashes == {name: value for name, value in EXPECTED_HASHES.items()
               if name != "primary_json"}
    and first_tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and second_tick["outcome"] == "SECOND_HOMOTHETIC_TICK_ACCEPTED"
    and third_tick["outcome"] == "THIRD_HOMOTHETIC_TICK_ACCEPTED"
    and second_tick["fixed_mass"] is True
    and third_tick["fixed_mass"] is True
    and shifted_direct["outcome"] == "SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS"
    and shifted_direct["passed"] == shifted_direct["tests"] == 14
)
check("all direct adversarial inputs have exact frozen provenance",
      provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_shifted_pseudolongitudinal_direct", GEOMETRY_SOURCE
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

vertices_mp, vertex_keys = exact_vertices()
frozen_vertices, frozen_adjacency, _ = build_600cell()
exact_edge_set = {
    tuple(sorted((left, right)))
    for left in range(120) for right in range(left + 1, 120)
    if frozen_adjacency[left, right]
}
vertex_control = bool(
    np.array_equal(vertex_keys, frozen_vertices)
    and len(exact_edge_set) == 720
    and max(abs(mp.fsum(value * value for value in vertex) - 1)
            for vertex in vertices_mp) < mp.mpf("1e-90")
)
check("independent exact golden-ratio vertices recover the frozen 600-cell",
      vertex_control)

print("=" * 78)
print("DIRECT HIGH-PRECISION SHIFTED PSEUDO-LONGITUDINAL AUDIT")
print("=" * 78)

records = []
pseudo_records = []
all_branch_controls = True
all_carrier_order = True
all_hessian_symmetry_enveloped = True
all_twists = True
all_principal_identities = True
all_shape_carriers = True
all_kinetic = True
all_compatibility = True
all_exact_carriers = True
all_pseudo_carriers = True
all_pseudo_inequalities = True
aggregate_signs = Counter()
principal_diagnostics = []
hessian_symmetry_diagnostics = []

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
        c_mp, d_mp, exact_edge_order = compressed_geometry(
            slab_data["slab2"]["index_data"], sector2, vertices_mp
        )
        compressed_incidence = mp_to_numpy(c_mp)
        compressed_tangent = mp_to_numpy(d_mp)
        exact_carrier_ok = bool(
            set(exact_edge_order) == exact_edge_set
            and exact_edge_order == groups[parity]["edge_order"]
        )
        all_exact_carriers &= exact_carrier_ok
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
                raw_block = raw_blocks[variant]
                antihermitian = (raw_block - raw_block.H) / 2
                symmetry_defect = mp_frobenius(antihermitian)
                symmetry_enveloped = bool(
                    symmetry_defect <= family_variation + ARITHMETIC_FLOOR
                )
                all_hessian_symmetry_enveloped &= symmetry_enveloped
                hessian_symmetry_diagnostics.append({
                    "parity": parity,
                    "slab": slab_name,
                    "sector_index": sector_index,
                    "variant": variant,
                    "symmetry_defect_frobenius": mp.nstr(symmetry_defect, 25),
                    "family_variation_frobenius": mp.nstr(family_variation, 25),
                    "enveloped": symmetry_enveloped,
                })
                block = (raw_block + raw_block.H) / 2
                _, determinant, tangent, _ = build_tangent_ball(
                    block, dimension, slab_data[slab_name]["mapping"]
                )
                all_twists &= not determinant.contains(0)
                tangent_blocks = split_tangent(tangent, n)
                principal, residuals = reconstruct_principal(tangent_blocks)
                for residual_name, residual in residuals.items():
                    residual_midpoint, residual_radii = acb_midpoint_and_radii(residual)
                    failing_entries = sum(
                        not residual[row, column].contains(0)
                        for row in range(residual.nrows())
                        for column in range(residual.ncols())
                    )
                    all_principal_identities &= failing_entries == 0
                    principal_diagnostics.append({
                        "parity": parity,
                        "slab": slab_name,
                        "sector_index": sector_index,
                        "variant": variant,
                        "identity": residual_name,
                        "failing_entries": failing_entries,
                        "midpoint_frobenius": sf(
                            la.norm(residual_midpoint, "fro")
                        ),
                        "radius_frobenius": sf(
                            la.norm(residual_radii, "fro")
                        ),
                        "maximum_midpoint_modulus": sf(
                            np.max(np.abs(residual_midpoint))
                        ),
                        "maximum_radius": sf(np.max(residual_radii)),
                    })
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
            m_mid, m_rad = direct_matrix(m_ball)
            v_mid, v_rad = direct_matrix(v_ball)
            omega_mid, omega_stored = acb_midpoint_and_radii(omega_ball)
            omega_rad = component_reenclosure_radii(omega_mid, omega_stored)

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

            pseudo = direct_relative_cell(
                m_mid, v_mid, m_rad, v_rad,
                compressed_incidence, compressed_tangent,
            )
            pseudo_ok = bool(
                pseudo["rank_C"] == pseudo["rank_shape_row"] == 5
                and pseudo["rank_longitudinal"] == 15
                and min(
                    pseudo["minimum_C_gap"],
                    pseudo["minimum_shape_gap"],
                    pseudo["minimum_longitudinal_gap"],
                ) > 1e-6
                and pseudo["minimum_B_eigenvalue"] > 1e-3
                and pseudo["negative_stiffness"] == 15
                and pseudo["positive_stiffness"] == 10
                and min(pseudo["AL_norm"], pseudo["GL_norm"]) > 1e-12
                and math.isfinite(pseudo["relative_floor"])
            )
            all_pseudo_carriers &= pseudo_ok
            all_pseudo_inequalities &= bool(
                pseudo["inequality_one"] and pseudo["inequality_two"]
            )

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
                    "exact_geometry_carrier": exact_carrier_ok,
                    "direct_M_radius_frobenius": sf(la.norm(m_rad, "fro")),
                    "direct_V_radius_frobenius": sf(la.norm(v_rad, "fro")),
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
                    "base_restricted_error": sf(epsilon_vs),
                    "schedule_eigenvalue_variation": None,
                    "minimum_absolute_error_units": sf(
                        np.min(np.abs(a_values)) / epsilon_vs
                        if epsilon_vs else math.inf
                    ),
                    "sign_counts": dict(counts),
                    "labels": labels,
                    "eigenvalues": [sf(value) for value in a_values],
                },
                "_eigenvalues_float": [float(value) for value in a_values],
                "_base_error_float": epsilon_vs,
                "_pseudo": pseudo,
            })

for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        family = [
            record for record in records
            if record["parity"] == parity
            and record["sector_index"] == sector_index
        ]
        primary = next(
            record for record in family
            if record["variant"] == VARIANTS[0]
        )
        primary_values = np.asarray(primary["_eigenvalues_float"], dtype=float)
        schedule_variation = max(
            float(np.max(np.abs(
                np.asarray(record["_eigenvalues_float"], dtype=float)
                - primary_values
            )))
            for record in family
        )
        primary_pseudo = primary["_pseudo"]
        pseudo_errors = {}
        for field in ("relative_span", "relative_commutator"):
            pseudo_errors[field] = float(
                max(abs(record["_pseudo"][field] - primary_pseudo[field])
                    for record in family)
                + max(record["_pseudo"]["relative_floor"]
                      for record in family)
            )
        augmented_error = max(pseudo_errors.values())
        for record in family:
            effective_error = record["_base_error_float"] + schedule_variation
            values = np.asarray(record["_eigenvalues_float"], dtype=float)
            labels = [sign_label(float(value), effective_error) for value in values]
            counts = Counter(labels)
            aggregate_signs.update(counts)
            record["stiffness"]["restricted_error"] = sf(effective_error)
            record["stiffness"]["schedule_eigenvalue_variation"] = sf(
                schedule_variation
            )
            record["stiffness"]["minimum_absolute_error_units"] = sf(
                np.min(np.abs(values)) / effective_error
                if effective_error else math.inf
            )
            record["stiffness"]["sign_counts"] = dict(counts)
            record["stiffness"]["labels"] = labels
            pseudo = record["_pseudo"]
            pseudo["span_label"] = scalar_zero_label(
                pseudo["relative_span"], pseudo_errors["relative_span"]
            )
            pseudo["commutator_label"] = scalar_zero_label(
                pseudo["relative_commutator"],
                pseudo_errors["relative_commutator"],
            )
            threshold = float(
                100 * augmented_error * pseudo["augmented_singular"][0]
            )
            pseudo["augmented_rank"] = int(np.sum(
                pseudo["augmented_singular"] > threshold
            ))
            pseudo["augmented_threshold"] = threshold
            pseudo["sixteenth_augmented_singular"] = float(
                pseudo["augmented_singular"][15]
            )
            del pseudo["augmented_singular"]
            record["pseudolongitudinal_errors"] = {
                name: sf(value) for name, value in pseudo_errors.items()
            }
            record["pseudolongitudinal"] = {
                name: (sf(value) if isinstance(value, (float, np.floating))
                       else value)
                for name, value in pseudo.items()
            }
            pseudo_records.append(pseudo)
            del record["_eigenvalues_float"]
            del record["_base_error_float"]
            del record["_pseudo"]

check("all four direct slab branch reconstructions pass", all_branch_controls)
check("both slabs retain the same ordered sector carrier", all_carrier_order)
check("all 32 raw Hessian symmetry defects lie inside their schedule-family variation",
      all_hessian_symmetry_enveloped)
check("all 32 direct boundary-twist determinant balls exclude zero", all_twists)
check("all 32 principal-function identity families contain zero entrywise",
      all_principal_identities,
      str(Counter(
          item["identity"]
          for item in principal_diagnostics if item["failing_entries"]
      )))
check("all four exact sector carriers recover the complete ordered edge set",
      all_exact_carriers)
check("all 16 action-selected shape carriers are resolved", all_shape_carriers)
check("all 16 restricted kinetic forms are positive-definite resolved", all_kinetic)
check("all 16 action-compatibility residuals are zero-consistent", all_compatibility)
check("all 16 direct pseudo-longitudinal carriers and denominators are resolved",
      all_pseudo_carriers)
check("both exact residual inequalities hold in all 16 direct cells",
      all_pseudo_inequalities)

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

span_counts = Counter(record["span_label"] for record in pseudo_records)
commutator_counts = Counter(
    record["commutator_label"] for record in pseudo_records
)
rank_counts = Counter(record["augmented_rank"] for record in pseudo_records)
pseudo_complete = bool(
    len(pseudo_records) == 16
    and sum(span_counts.values()) == 16
    and sum(commutator_counts.values()) == 16
    and sum(rank_counts.values()) == 16
)
check("both direct residuals and augmented ranks classify all 16 cells",
      pseudo_complete,
      f"span={dict(span_counts)}, comm={dict(commutator_counts)}, "
      f"ranks={dict(rank_counts)}")

# No primary residual byte is read before the complete direct census above.
hashes["primary_json"] = sha256(PRIMARY_JSON)
primary = json.loads(PRIMARY_JSON.read_text())
primary_ok = bool(
    hashes == EXPECTED_HASHES
    and primary["outcome"] == "SHIFTED_PSEUDOLONGITUDINAL_DEFECT_PERSISTS"
    and primary["passed"] == primary["tests"] == 10
    and primary["label_counts"]["relative_span"]
    == {"NONZERO_RESOLVED": 16}
    and primary["label_counts"]["relative_commutator"]
    == {"NONZERO_RESOLVED": 16}
)
post_comparison = []
comparison_ok = primary_ok
if primary_ok:
    for record in records:
        cell = f"{record['parity']}_sector{record['sector_index']}"
        variant = record["variant"]
        direct = record["pseudolongitudinal"]
        archived = primary["cells"][cell]["variants"][variant]
        item = {
            "cell": cell,
            "variant": variant,
            "direct_over_primary_span": (
                float(direct["relative_span"])
                / float(archived["relative_span"])
            ),
            "direct_over_primary_commutator": (
                float(direct["relative_commutator"])
                / float(archived["relative_commutator"])
            ),
            "direct_minus_primary_span": (
                float(direct["relative_span"])
                - float(archived["relative_span"])
            ),
            "direct_minus_primary_commutator": (
                float(direct["relative_commutator"])
                - float(archived["relative_commutator"])
            ),
        }
        comparison_ok &= all(
            math.isfinite(value) for name, value in item.items()
            if name not in {"cell", "variant"}
        )
        post_comparison.append(item)
check("the primary artifact is opened only post-census and comparison is finite",
      comparison_ok)

controls_ok = bool(
    provenance_ok and geometry_import_ok and incidence_ok
    and vertex_control and all_branch_controls and all_carrier_order and all_twists
    and all_hessian_symmetry_enveloped and all_principal_identities
    and all_exact_carriers and all_shape_carriers and all_kinetic
    and all_compatibility and all_pseudo_carriers and all_pseudo_inequalities
    and len(records) == 16 and pseudo_complete and comparison_ok
)
persists = bool(
    controls_ok
    and len(target_cells) == 16
    and span_counts == {"NONZERO_RESOLVED": 16}
    and commutator_counts == {"NONZERO_RESOLVED": 16}
    and all(record["augmented_rank"] > 15 for record in pseudo_records)
)
refutes = bool(
    controls_ok
    and len(target_cells) == 16
    and span_counts == {"ZERO_CONSISTENT": 16}
    and commutator_counts == {"ZERO_CONSISTENT": 16}
    and all(record["augmented_rank"] == 15 for record in pseudo_records)
)
if not controls_ok:
    outcome = "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_CONTROL_FAILED"
elif persists:
    outcome = "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_PERSISTENCE_CONFIRMED"
elif refutes:
    outcome = "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_PERSISTENCE_REFUTED"
else:
    outcome = "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_OPEN"

allowed = {
    "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_CONTROL_FAILED",
    "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_PERSISTENCE_CONFIRMED",
    "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_PERSISTENCE_REFUTED",
    "SHIFTED_PSEUDOLONGITUDINAL_DIRECT_OPEN",
}
check("the preregistered direct adversarial hierarchy assigns the outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "centered_MV_archive_loaded": False,
    "primary_residual_loaded_only_post_census": True,
    "target_sectors": list(TARGET_SECTORS),
    "cells": len(records),
    "outcome": outcome,
    "controls_ok": controls_ok,
    "complete_cells": len(complete_cells),
    "target_rank_cells": len(target_cells),
    "changed_cells": len(changed_cells),
    "aggregate_sign_counts": dict(aggregate_signs),
    "pseudolongitudinal_label_counts": {
        "relative_span": dict(span_counts),
        "relative_commutator": dict(commutator_counts),
        "augmented_rank": {str(key): value for key, value in rank_counts.items()},
    },
    "post_census_comparison": [
        {
            name: (sf(value) if isinstance(value, (float, np.floating)) else value)
            for name, value in item.items()
        }
        for item in post_comparison
    ],
    "records": records,
    "hessian_symmetry_diagnostics": hessian_symmetry_diagnostics,
    "principal_identity_diagnostics": principal_diagnostics,
    "classification": {
        "two_tick_temporal_persistence": (
            "CONFIRMED DERIVED COMPUTATIONAL / STRUCTURAL"
            if persists else "REFUTED" if refutes else "OPEN"
        ),
        "curvature_or_refinement_limit": "NOT TESTED",
        "conservation_monotonicity_or_lifetime": "NOT TESTED",
        "physical_instability_or_gauge_recovery": "OPEN",
        "continuum_wave_graviton_or_speed": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"target 15+10 cells: {len(target_cells)}/16")
print(f"aggregate signs: {dict(aggregate_signs)}")
print(f"relative span labels: {dict(span_counts)}")
print(f"relative commutator labels: {dict(commutator_counts)}")
print(f"augmented ranks: {dict(rank_counts)}")
if post_comparison:
    span_ratios = [item["direct_over_primary_span"] for item in post_comparison]
    comm_ratios = [item["direct_over_primary_commutator"] for item in post_comparison]
    print("direct/primary span ratio range:", min(span_ratios), max(span_ratios))
    print("direct/primary commutator ratio range:", min(comm_ratios), max(comm_ratios))
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

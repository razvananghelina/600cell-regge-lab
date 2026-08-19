#!/usr/bin/env python3
"""Direct high-precision falsifier for the action-York longitudinal identity."""

import ast
from collections import Counter, defaultdict
import contextlib
from hashlib import sha256
import importlib.util
import io
from itertools import combinations, permutations, product
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_dust_action_york_direct_precision.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_dust_action_york_direct_precision_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_action_york_direct_precision_protocol.md"
PREVIOUS_RESULT = ROOT / "docs/gravity/gravity_600cell_dust_action_york_negative_result.md"
PREVIOUS_JSON = HERE / "gravity_600cell_dust_action_york_negative.json"
PREVIOUS_SOURCE = HERE / "verify_gravity_600cell_dust_action_york_negative.py"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
TWO_STEP_SOURCE = HERE / "verify_gravity_600cell_dust_two_step_full_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
FIRST_TICK = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
SECOND_TICK = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
COMMONS_SOURCE = ROOT / "commons/cell600.py"

PRIOR_ART_COMMIT = "1b076f2"
PROTOCOL_COMMIT = "167077b"
INITIAL_PROTOCOL_COMMIT = "852cd54"
EXPECTED_HASHES = {
    "prior_art": "6cff580ed92893dada0099ad3330afbadc4d7a4d67b34984747a699e20a3b057",
    "protocol": "e655d6025e790ff2beb653c5e9f4c2f38233606c3607f9219ae222bafdfed36e",
    "previous_result": "338b2fca275c8a7af0f866a680512b88b1be5698a028e0588580b47c9e463c87",
    "previous_json": "fd0763af779cb02d96f7e1d7a8856b117dd4bf2c9413f01de6246c597743df27",
    "previous_source": "370bae86c27e82f9dda4592e8db1774786a1d2c1919ed96e3fceb6e372e6be7b",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "two_step_source": "c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "first_tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "second_tick": "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70",
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}

DPS = 120
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-25"),
    "operational_shadow": mp.mpf("1e-18"),
    "validation_primary": mp.mpf("3e-25"),
    "validation_shadow": mp.mpf("3e-18"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
PARITIES = ("even", "odd")
SELECTED_SECTORS = (4, 5)
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {hinge: index for index, hinge in enumerate(LOCAL_HINGES)}
I = mp.mpc(0, 1)
ARITHMETIC_FLOOR = mp.mpf("1e-90")
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


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def sf(value):
    return f"{float(value):.17e}"


def operator_norm(matrix):
    singular = la.svdvals(matrix)
    return float(singular[0]) if len(singular) else 0.0


def zero_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "OPEN"
    if value <= 10 * error:
        return "ZERO_CONSISTENT"
    if value > 100 * error:
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


def load_named_functions(source, wanted):
    tree = ast.parse(source.read_text(), filename=str(source))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch in {source.name}: {wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(source), "exec"),
        globals(),
    )


load_named_functions(RANK_SOURCE, {
    "orbit_sort_key", "augment_boundary_orbits", "log_minus",
    "signed_volume_square", "angle_record", "area_data",
    "extended_edge_image", "group_and_index_data", "prepare_geometry",
})
load_named_functions(TANGENT_SOURCE, {
    "mp_frobenius", "mp_submatrix", "cluster_sorted",
    "high_precision_sector_bases", "high_precision_pattern_cache",
    "assemble_full_representative_kernels", "project_full_kernel",
    "expanded_types",
})


def mp_to_numpy(matrix):
    return np.asarray([
        [complex(float(mp.re(matrix[row, column])),
                 float(mp.im(matrix[row, column])))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


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
        raise RuntimeError("slab left the positive magnitude domain")
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


def build_tangent_mp(block, dimension, mapping):
    old = expanded_types(0, 30, dimension)
    internal = expanded_types(30, 65, dimension)
    new = expanded_types(65, 95, dimension)
    nd = 30 * dimension
    xd = 35 * dimension

    k_xx = mp_submatrix(block, internal, internal)
    k_xn = mp_submatrix(block, internal, new)
    k_ox = mp_submatrix(block, old, internal)
    k_on = mp_submatrix(block, old, new)
    j_matrix = mp.matrix(xd + nd, xd + nd)
    for row in range(xd):
        for column in range(xd):
            j_matrix[row, column] = k_xx[row, column]
        for column in range(nd):
            j_matrix[row, xd + column] = k_xn[row, column]
    for row in range(nd):
        for column in range(xd):
            j_matrix[xd + row, column] = -k_ox[row, column]
        for column in range(nd):
            j_matrix[xd + row, xd + column] = -k_on[row, column]

    rhs = mp.matrix(xd + nd, 2 * nd)
    k_xo = mp_submatrix(block, internal, old)
    k_oo = mp_submatrix(block, old, old)
    for row in range(xd):
        for column in range(nd):
            rhs[row, column] = -k_xo[row, column]
    for row in range(nd):
        for column in range(nd):
            rhs[xd + row, column] = k_oo[row, column]
        rhs[xd + row, nd + row] = 1

    solved = (j_matrix ** -1) * rhs
    y_x = mp.matrix(xd, 2 * nd)
    y_n = mp.matrix(nd, 2 * nd)
    for row in range(xd):
        for column in range(2 * nd):
            y_x[row, column] = solved[row, column]
    for row in range(nd):
        for column in range(2 * nd):
            y_n[row, column] = solved[xd + row, column]

    k_no = mp_submatrix(block, new, old)
    k_nx = mp_submatrix(block, new, internal)
    k_nn = mp_submatrix(block, new, new)
    direct = mp.matrix(nd, 2 * nd)
    for row in range(nd):
        for column in range(nd):
            direct[row, column] = k_no[row, column]
    p_post = direct + k_nx * y_x + k_nn * y_n
    raw = mp.matrix(2 * nd, 2 * nd)
    for row in range(nd):
        for column in range(2 * nd):
            raw[row, column] = y_n[row, column]
            raw[nd + row, column] = p_post[row, column]

    tangent = mp.matrix(2 * nd, 2 * nd)
    for next_type, final_type in enumerate(mapping):
        for component in range(dimension):
            target = next_type * dimension + component
            source = final_type * dimension + component
            for column in range(2 * nd):
                tangent[target, column] = raw[source, column]
                tangent[nd + target, column] = raw[nd + source, column]
    return tangent, j_matrix


def principal_blocks(tangent):
    nd = tangent.rows // 2
    q = list(range(nd))
    p = list(range(nd, 2 * nd))
    a = mp_submatrix(tangent, q, q)
    b = mp_submatrix(tangent, q, p)
    c = mp_submatrix(tangent, p, q)
    d = mp_submatrix(tangent, p, p)
    inverse_b = b ** -1
    return {
        "00": inverse_b * a,
        "01": -inverse_b,
        "10": c - d * inverse_b * a,
        "11": d * inverse_b,
    }, b


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
        raise RuntimeError("selected carrier is not a one-dimensional minimal sector")
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
            delta = [vertices[left][axis] - vertices[right][axis] for axis in range(4)]
            coefficient = [2 * value / length_square for value in delta]
            for axis in range(4):
                tangent[orbit, 4 * left + axis] += weight * (
                    coefficient[axis] - vertices[left][axis]
                )
                tangent[orbit, 4 * right + axis] += weight * (
                    -coefficient[axis] - vertices[right][axis]
                )
    return incidence, tangent, tuple(edge_order)


def analyze_cell(m_matrix, v_matrix, c_matrix, d_matrix):
    h_m = (m_matrix + m_matrix.conj().T) / 2
    h_v = (v_matrix + v_matrix.conj().T) / 2
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

    m_s = (w_basis.conj().T @ h_m @ w_basis)
    m_s = (m_s + m_s.conj().T) / 2
    v_s = (w_basis.conj().T @ h_v @ w_basis)
    v_s = (v_s + v_s.conj().T) / 2
    a_matrix = -v_s
    b_matrix = -m_s
    row_t = l_basis.conj().T @ b_matrix
    _, singular_t, right_t = la.svd(row_t, full_matrices=True)
    t_basis = right_t.conj().T[:, 15:]

    cross = operator_norm(l_basis.conj().T @ a_matrix @ t_basis)
    reduced_b = l_basis.conj().T @ b_matrix @ l_basis
    reduced_a = l_basis.conj().T @ a_matrix @ l_basis
    generator = la.solve(reduced_b, reduced_a)
    image_residual = operator_norm(
        a_matrix @ l_basis - b_matrix @ l_basis @ generator
    )

    generalized_values, generalized_vectors = la.eigh(a_matrix, b_matrix)
    negative_basis, _ = la.qr(generalized_vectors[:, :15], mode="economic")
    p_l = l_basis @ l_basis.conj().T
    p_negative = negative_basis @ negative_basis.conj().T
    projector_distance = operator_norm(p_l - p_negative)
    control_seed = np.column_stack((l_basis[:, 1:], t_basis[:, :1]))
    control_basis, _ = la.qr(control_seed, mode="economic")
    control_distance = operator_norm(
        control_basis @ control_basis.conj().T - p_negative
    )
    values_l = la.eigvalsh((reduced_a + reduced_a.conj().T) / 2)
    restricted_t = t_basis.conj().T @ a_matrix @ t_basis
    values_t = la.eigvalsh((restricted_t + restricted_t.conj().T) / 2)
    b_values = la.eigvalsh(b_matrix)

    kappa = max(
        float(np.linalg.cond(q_basis)),
        float(singular_l[0] / singular_l[14]),
        float(np.linalg.cond(reduced_b)),
    )
    scale = max(
        1.0, operator_norm(a_matrix), operator_norm(b_matrix),
        operator_norm(d_matrix),
    )
    arithmetic_floor = float(1000 * MACHINE_EPSILON * 30 * kappa * scale)
    return {
        "rank_C": int(np.sum(singular_c > 1e-6)),
        "shape_rank": int(np.sum(singular_row > 1e-6)),
        "longitudinal_rank": int(np.sum(singular_l > 1e-6)),
        "transverse_dimension": t_basis.shape[1],
        "minimum_C_gap": float(singular_c[4]),
        "minimum_shape_gap": float(singular_row[4]),
        "minimum_longitudinal_gap": float(singular_l[14]),
        "minimum_B": float(b_values[0]),
        "kappa_carrier": kappa,
        "arithmetic_floor": arithmetic_floor,
        "cross_residual": cross,
        "image_residual": image_residual,
        "projector_distance": projector_distance,
        "longitudinal_maximum": float(values_l[-1]),
        "transverse_minimum": float(values_t[0]),
        "rotated_control_distance": control_distance,
        "generalized_negative": int(np.sum(generalized_values < 0)),
        "generalized_positive": int(np.sum(generalized_values > 0)),
    }


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "previous_result": PREVIOUS_RESULT,
    "previous_json": PREVIOUS_JSON,
    "previous_source": PREVIOUS_SOURCE,
    "centered_npz": CENTERED_NPZ,
    "tangent_source": TANGENT_SOURCE,
    "two_step_source": TWO_STEP_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "rank_source": RANK_SOURCE,
    "first_tick": FIRST_TICK,
    "second_tick": SECOND_TICK,
    "commons": COMMONS_SOURCE,
}
hashes = {name: digest(path) for name, path in paths.items()}
previous = json.loads(PREVIOUS_JSON.read_text())
first_tick = json.loads(FIRST_TICK.read_text())
second_tick = json.loads(SECOND_TICK.read_text())
provenance_ok = hashes == EXPECTED_HASHES
check("all target-disclosed protocol inputs have exact provenance",
      provenance_ok, str(hashes))
previous_ok = bool(
    previous["outcome"] == "NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_OPEN"
    and previous["passed"] == 13 and previous["tests"] == 15
    and previous["projector_label_counts"] == {"EQUALITY_CONSISTENT": 16}
)
check("the preceding identity result remains the frozen 13/15 OPEN outcome",
      previous_ok)

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_action_york_direct", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check("the independent literal slab geometry retains all 43 controls",
      gro.tests == gro.passed == 43)

M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON3 * L0
RHO0 = mp.mpf("0.0102") ** 2
models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}

state_controls = {}
states_ok = True
for parity in PARITIES:
    first = first_tick["solutions"][parity]
    second = second_tick["solutions"][parity]
    first_residual = mp.mpf(first["reduced_residual_norm"])
    second_residual = mp.mpf(second["reduced_residual_norm"])
    first_junction = mp.mpf(first["junction_norm"])
    second_junction = mp.mpf(second["junction_norm"])
    first_bound = mp.mpf(first["junction_bound"])
    second_bound = mp.mpf(second["junction_bound"])
    ok = bool(
        first["converged"] and second["converged"]
        and first_residual < mp.mpf("1e-20")
        and second_residual < mp.mpf("1e-20")
        and first_junction <= first_bound and second_junction <= second_bound
    )
    states_ok &= ok
    state_controls[parity] = {
        "first_residual": first_residual,
        "second_residual": second_residual,
        "first_junction": first_junction,
        "second_junction": second_junction,
        "first_bound": first_bound,
        "second_bound": second_bound,
    }
check("both frozen backgrounds and canonical seams pass their residual bounds",
      states_ok)

vertices_mp, vertex_keys = exact_vertices()
frozen_vertices, frozen_adjacency, _ = build_600cell()
vertex_control = bool(
    np.array_equal(vertex_keys, frozen_vertices)
    and max(abs(mp.fsum(value * value for value in vertex) - 1)
            for vertex in vertices_mp) < mp.mpf("1e-110")
)
check("the exact golden-ratio vertices recover all 120 frozen labels",
      vertex_control)

print("=" * 78)
print("DIRECT 120-DIGIT ACTION-YORK RESIDUAL")
print("=" * 78)

direct_cells = {}
geometry_cells = {}
construction_controls = {}
all_geometry = True
all_branches = True
all_basis = True
all_legendre = True

for parity in PARITIES:
    print(f"[{parity}] reconstructing exact carrier and both local slabs", flush=True)
    first_solution = first_tick["solutions"][parity]
    second_solution = second_tick["solutions"][parity]
    a1, r1 = [mp.mpf(value) for value in first_solution["state"]]
    a2, r2 = [mp.mpf(value) for value in second_solution["state_absolute"]]
    index1, values1 = slab_index_data(models[parity], mp.mpf(0), a1, r1)
    index2, values2 = slab_index_data(models[parity], a1, a2, r2)
    geometry1 = prepare_geometry(models[parity], index1)
    geometry2 = prepare_geometry(models[parity], index2)
    mapping1 = boundary_mapping(index1)
    mapping2 = boundary_mapping(index2)
    sectors1, sector_control1 = high_precision_sector_bases(index1)
    sectors2, sector_control2 = high_precision_sector_bases(index2)

    exact_edges = {
        tuple(sorted((left, right)))
        for left in range(120) for right in range(left + 1, 120)
        if frozen_adjacency[left, right]
    }
    orbit_edges = {
        tuple(sorted(edge))
        for orbit in range(30) for edge in index1["orbit_edges"][orbit]
    }
    geometry_ok = bool(
        len(exact_edges) == len(orbit_edges) == 720
        and exact_edges == orbit_edges
        and sorted(mapping1) == sorted(mapping2) == list(range(30))
        and len(geometry1["patterns"]) == len(geometry2["patterns"]) == 20
    )
    all_geometry &= geometry_ok
    maxima1 = max(
        value for key, value in sector_control1.items() if key.startswith("maximum_")
    )
    maxima2 = max(
        value for key, value in sector_control2.items() if key.startswith("maximum_")
    )
    basis_ok = bool(
        tuple(sector["dimension"] for sector in sectors1) == (3, 2, 2, 2, 1, 1, 1)
        and tuple(sector["dimension"] for sector in sectors2) == (3, 2, 2, 2, 1, 1, 1)
        and maxima1 < mp.mpf("1e-90") and maxima2 < mp.mpf("1e-90")
    )
    all_basis &= basis_ok
    check(f"{parity}: exact edges and selected symmetry carriers are fixed",
          geometry_ok and basis_ok)

    slab_data = []
    for slab_name, index_data, geometry, values in (
        ("first", index1, geometry1, values1),
        ("second", index2, geometry2, values2),
    ):
        print(f"[{parity}] differentiating {slab_name}-slab patterns", flush=True)
        cache, branch = high_precision_pattern_cache(geometry["patterns"], values)
        kernels, kernel_control = assemble_full_representative_kernels(
            index_data, geometry, cache
        )
        branch_ok = bool(
            branch["entry_pass"]
            and branch["base_negative_counts"] == Counter({1: 2400})
            and branch["displaced_negative_counts"] == Counter({1: 1600})
            and branch["minimum_leading_minor"] > 0
            and branch["minimum_argument"] > mp.mpf("1e-6")
            and kernel_control["maximum_imaginary"] < mp.mpf("1e-90")
            and len(set(kernel_control["nonzero_entries"].values())) == 1
        )
        all_branches &= branch_ok
        check(f"{parity}: {slab_name} branch and four derivative scales converge",
              branch_ok)
        slab_data.append((kernels, branch, kernel_control))

    parity_legendre = True
    for sector_index in SELECTED_SECTORS:
        sector1 = sectors1[sector_index]
        sector2 = sectors2[sector_index]
        c_mp, d_mp, edge_order = compressed_geometry(
            index1, sector1, vertices_mp
        )
        geometry_cells[(parity, sector_index)] = {
            "C": mp_to_numpy(c_mp),
            "D": mp_to_numpy(d_mp),
            "edge_order": edge_order,
        }
        direct_cells[(parity, sector_index)] = {}
        print(f"[{parity}] sector {sector_index}: solving 8 direct Legendre maps",
              flush=True)
        for variant in VARIANTS:
            block1 = project_full_kernel(slab_data[0][0][variant], sector1)
            block2 = project_full_kernel(slab_data[1][0][variant], sector2)
            tangent1, j1 = build_tangent_mp(block1, 1, mapping1)
            tangent2, j2 = build_tangent_mp(block2, 1, mapping2)
            principal1, twist1 = principal_blocks(tangent1)
            principal2, twist2 = principal_blocks(tangent2)
            kminus = principal1["10"]
            kzero = principal1["11"] + principal2["00"]
            kplus = principal2["01"]
            m_matrix = (kminus + kplus) / 2
            v_matrix = kminus + kzero + kplus

            tangent1_np = mp_to_numpy(tangent1)
            tangent2_np = mp_to_numpy(tangent2)
            omega = np.block([
                [np.zeros((30, 30)), np.eye(30)],
                [-np.eye(30), np.zeros((30, 30))],
            ]).astype(np.complex128)
            symplectic1 = operator_norm(tangent1_np.conj().T @ omega @ tangent1_np - omega)
            symplectic2 = operator_norm(tangent2_np.conj().T @ omega @ tangent2_np - omega)
            minimum_j = min(
                la.svdvals(mp_to_numpy(j1))[-1],
                la.svdvals(mp_to_numpy(j2))[-1],
            )
            minimum_twist = min(
                la.svdvals(mp_to_numpy(twist1))[-1],
                la.svdvals(mp_to_numpy(twist2))[-1],
            )
            legendre_ok = bool(
                minimum_j > 1e-10 and minimum_twist > 1e-10
                and symplectic1 < 1e-7 and symplectic2 < 1e-7
            )
            parity_legendre &= legendre_ok
            direct_cells[(parity, sector_index)][variant] = {
                "M": mp_to_numpy(m_matrix),
                "V": mp_to_numpy(v_matrix),
                "minimum_J": float(minimum_j),
                "minimum_twist": float(minimum_twist),
                "symplectic": float(max(symplectic1, symplectic2)),
            }
    all_legendre &= parity_legendre
    check(f"{parity}: all selected direct Legendre maps remain regular and canonical",
          parity_legendre)
    construction_controls[parity] = {
        "geometry": geometry_ok,
        "basis": basis_ok,
        "branches": bool(all(
            item[1]["entry_pass"] for item in slab_data
        )),
        "legendre": parity_legendre,
    }

# The protocol permits opening this archive only after every direct M,V exists.
centered = np.load(CENTERED_NPZ)
archive_distances = []
for (parity, sector_index), variants in direct_cells.items():
    for variant, matrices in variants.items():
        prefix = f"{parity}_sector{sector_index}_{variant}"
        for name in ("M", "V"):
            archive_distances.append(operator_norm(
                matrices[name] - centered[f"{prefix}_{name}_midpoint"]
            ))
maximum_archive_distance = max(archive_distances)
archive_ok = maximum_archive_distance <= 1e-9
check("all direct M,V matrices reproduce the frozen centered midpoints",
      archive_ok, f"maximum distance={maximum_archive_distance:.3e}")

analysis = {}
carrier_ok = True
for key, variants in direct_cells.items():
    c_matrix = geometry_cells[key]["C"]
    d_matrix = geometry_cells[key]["D"]
    records = {
        variant: analyze_cell(
            matrices["M"], matrices["V"], c_matrix, d_matrix
        )
        for variant, matrices in variants.items()
    }
    primary = records["operational_primary"]
    fields = (
        "cross_residual", "image_residual", "projector_distance",
        "longitudinal_maximum", "transverse_minimum",
        "rotated_control_distance",
    )
    errors = {}
    for field in fields:
        variation = max(
            abs(record[field] - primary[field]) for record in records.values()
        )
        errors[field] = float(
            variation + max(record["arithmetic_floor"] for record in records.values())
        )
    for record in records.values():
        record["cross_label"] = zero_label(
            record["cross_residual"], errors["cross_residual"]
        )
        record["image_label"] = zero_label(
            record["image_residual"], errors["image_residual"]
        )
        record["projector_label"] = zero_label(
            record["projector_distance"], errors["projector_distance"]
        )
        record["longitudinal_sign"] = sign_label(
            record["longitudinal_maximum"], errors["longitudinal_maximum"]
        )
        record["transverse_sign"] = sign_label(
            record["transverse_minimum"], errors["transverse_minimum"]
        )
        record["control_label"] = zero_label(
            record["rotated_control_distance"], errors["rotated_control_distance"]
        )
        carrier_ok &= bool(
            record["rank_C"] == record["shape_rank"] == 5
            and record["longitudinal_rank"] == 15
            and record["transverse_dimension"] == 10
            and min(
                record["minimum_C_gap"], record["minimum_shape_gap"],
                record["minimum_longitudinal_gap"],
            ) > 1e-6
            and record["minimum_B"] > 1e-3
            and record["generalized_negative"] == 15
            and record["generalized_positive"] == 10
        )
    analysis[key] = {"records": records, "errors": errors}

check("all sixteen cells retain the geometry-selected 15+10 carrier and signs",
      carrier_ok)

all_records = [
    record
    for family in analysis.values() for record in family["records"].values()
]
cross_counts = Counter(record["cross_label"] for record in all_records)
image_counts = Counter(record["image_label"] for record in all_records)
projector_counts = Counter(record["projector_label"] for record in all_records)
longitudinal_counts = Counter(record["longitudinal_sign"] for record in all_records)
transverse_counts = Counter(record["transverse_sign"] for record in all_records)
control_counts = Counter(record["control_label"] for record in all_records)
complete_labels = bool(
    len(all_records) == 16
    and sum(cross_counts.values()) == sum(image_counts.values()) == 16
    and sum(projector_counts.values()) == 16
    and sum(longitudinal_counts.values()) == sum(transverse_counts.values()) == 16
    and sum(control_counts.values()) == 16
)
check("both direct invariance residuals receive all sixteen frozen labels",
      complete_labels,
      f"cross={dict(cross_counts)}, image={dict(image_counts)}")

controls_ok = bool(
    provenance_ok and previous_ok and gro.tests == gro.passed == 43
    and states_ok and vertex_control and all_geometry and all_basis
    and all_branches and all_legendre and archive_ok and carrier_ok
    and complete_labels
)
numerically_refuted = bool(
    controls_ok and all(
        record["cross_label"] == "NONZERO_RESOLVED"
        or record["image_label"] == "NONZERO_RESOLVED"
        for record in all_records
    )
)
numerically_resolved = bool(
    controls_ok
    and cross_counts == {"ZERO_CONSISTENT": 16}
    and image_counts == {"ZERO_CONSISTENT": 16}
    and projector_counts == {"ZERO_CONSISTENT": 16}
    and longitudinal_counts == {"NEGATIVE_RESOLVED": 16}
    and transverse_counts == {"POSITIVE_RESOLVED": 16}
    and control_counts == {"NONZERO_RESOLVED": 16}
)
if not controls_ok:
    outcome = "DIRECT_ACTION_YORK_CONTROL_FAILED"
elif numerically_refuted:
    outcome = "DIRECT_LONGITUDINAL_IDENTITY_NUMERICALLY_REFUTED"
elif numerically_resolved:
    outcome = "DIRECT_LONGITUDINAL_IDENTITY_NUMERICALLY_RESOLVED"
else:
    outcome = "DIRECT_LONGITUDINAL_IDENTITY_OPEN"
allowed = {
    "DIRECT_ACTION_YORK_CONTROL_FAILED",
    "DIRECT_LONGITUDINAL_IDENTITY_NUMERICALLY_REFUTED",
    "DIRECT_LONGITUDINAL_IDENTITY_NUMERICALLY_RESOLVED",
    "DIRECT_LONGITUDINAL_IDENTITY_OPEN",
}
check("the preregistered hierarchy assigns exactly one direct outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "initial_protocol_commit": INITIAL_PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "mpmath_decimal_precision": DPS,
    "derivative_steps": {name: mp.nstr(value, 20)
                         for name, value in DERIVATIVE_STEPS.items()},
    "state_controls": {
        parity: {name: mp.nstr(value, 50) for name, value in record.items()}
        for parity, record in state_controls.items()
    },
    "construction_controls": construction_controls,
    "maximum_archive_midpoint_distance": sf(maximum_archive_distance),
    "cells": {
        f"{parity}_sector{sector}": {
            "errors": {name: sf(value) for name, value in family["errors"].items()},
            "variants": {
                variant: {
                    name: (sf(value) if isinstance(value, (float, np.floating)) else value)
                    for name, value in record.items()
                }
                for variant, record in family["records"].items()
            },
        }
        for (parity, sector), family in analysis.items()
    },
    "label_counts": {
        "cross": dict(cross_counts),
        "image": dict(image_counts),
        "projector": dict(projector_counts),
        "longitudinal_sign": dict(longitudinal_counts),
        "transverse_sign": dict(transverse_counts),
        "rotated_control": dict(control_counts),
    },
    "classification": {
        "fixed_curved_slab_exact_longitudinal_identity": (
            "REFUTED DERIVED COMPUTATIONAL / STRUCTURAL"
            if numerically_refuted else
            "DERIVED COMPUTATIONAL / STRUCTURAL"
            if numerically_resolved else "OPEN"
        ),
        "formal_interval_or_symbolic_theorem": False,
        "pseudo_longitudinal_interpretation": (
            "DERIVED COMPUTATIONAL / STRUCTURAL"
            if numerically_refuted else "OPEN"
        ),
        "continuum_diffeomorphism_claim": "NOT TESTED",
        "polarization_speed_mass": "NOT COMPUTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
print("cross labels:", dict(cross_counts))
print("image labels:", dict(image_counts))
print("projector labels:", dict(projector_counts))
print("longitudinal signs:", dict(longitudinal_counts))
print("transverse signs:", dict(transverse_counts))
print("rotated controls:", dict(control_counts))
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

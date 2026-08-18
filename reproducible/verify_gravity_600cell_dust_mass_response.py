#!/usr/bin/env python3
"""Conserved inhomogeneous dust-mass response of the full 600-cell tick.

Prior-art commit: 6b587ea.
Protocol commit: a0b253d.
The 119+1 count and the earlier expanding/weak-lift proximity were disclosed.
No continuum, speed, Planck, particle, or Standard-Model target is loaded.
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


HERE = Path(__file__).resolve().parent
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
RANK_INPUT = HERE / "gravity_600cell_dust_full_anisotropic_legendre_rank.json"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
TANGENT_INPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
SCHUR_INPUT = HERE / "gravity_600cell_dust_full_lapse_schur.json"
CURVATURE_INPUT = HERE / "gravity_600cell_dust_internal_curvature_response.json"
CURVATURE_SOURCE = HERE / "verify_gravity_600cell_dust_internal_curvature_response.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
GLUING_INPUT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_mass_response.json"
NUMERIC_OUTPUT = HERE / "gravity_600cell_dust_mass_response.npz"

PRIOR_ART_COMMIT = "6b587ea"
PROTOCOL_COMMIT = "a0b253d"
EXPECTED_HASHES = {
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "rank": "7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "tangent": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "schur": "4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349",
    "curvature": "95b6edd8e21ad20a0db97a7c8e7027db7da6547b2b994ad1eb595cf2307f29dc",
    "curvature_source": "276982879fae5f8fa735f27a6fa30bfe965dc3e41c169d8a229a61c23511ae66",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "gluing": "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77",
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
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {hinge: index for index, hinge in enumerate(LOCAL_HINGES)}
I = mp.mpc(0, 1)
ARITHMETIC_FLOOR = mp.mpf("1e-70")
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
    value = float(value)
    if not math.isfinite(value):
        return str(value)
    return format(value, ".17g")


paths = {
    "tick": TICK_INPUT,
    "rank": RANK_INPUT,
    "rank_source": RANK_SOURCE,
    "tangent": TANGENT_INPUT,
    "tangent_numeric": TANGENT_NUMERIC,
    "tangent_source": TANGENT_SOURCE,
    "schur": SCHUR_INPUT,
    "curvature": CURVATURE_INPUT,
    "curvature_source": CURVATURE_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "gluing": GLUING_INPUT,
}
hashes = {name: sha256(path) for name, path in paths.items()}
tick = json.loads(TICK_INPUT.read_text())
rank_input = json.loads(RANK_INPUT.read_text())
tangent_input = json.loads(TANGENT_INPUT.read_text())
schur_input = json.loads(SCHUR_INPUT.read_text())
curvature_input = json.loads(CURVATURE_INPUT.read_text())
gluing = json.loads(GLUING_INPUT.read_text())
numeric = np.load(TANGENT_NUMERIC)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and rank_input["outcome"] == "FULL_CANONICAL_LEGENDRE_REGULAR"
    and tangent_input["outcome"] == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and tangent_input["numeric_archive_arrays"] == len(numeric.files) == 224
    and tangent_input["numeric_archive_sha256"] == hashes["tangent_numeric"]
    and schur_input["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
    and curvature_input["outcome"] == "STRONG_TANGENT_CURVATURE_INJECTIVE"
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
)
check("all preregistered inputs have exact frozen provenance", provenance_ok, str(hashes))


spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_mass_response", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check(
    "the direct one-slab geometry import retains all 43 certificates",
    gro.tests == gro.passed == 43,
)


M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON3 * L0
RHO0 = mp.mpf("0.0102") ** 2


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
    "orbit_sort_key",
    "augment_boundary_orbits",
    "log_minus",
    "signed_volume_square",
    "angle_record",
    "area_data",
    "extended_edge_image",
    "group_and_index_data",
    "prepare_geometry",
})
load_named_functions(TANGENT_SOURCE, {
    "mp_frobenius",
    "mp_submatrix",
    "cluster_sorted",
    "high_precision_sector_bases",
    "high_precision_pattern_cache",
    "assemble_full_representative_kernels",
    "project_full_kernel",
    "mp_to_numpy",
    "mp_to_acb",
    "mp_matrix_to_acb",
    "acb_midpoint_and_radii",
    "expanded_types",
})
load_named_functions(CURVATURE_SOURCE, {
    "triangle_area_square",
    "extended_triangle_image",
    "group_inverses",
    "triangle_response_data",
    "project_curvature_kernel",
})

models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}


def acb_zero_matrix(rows, columns):
    return acb_mat(rows, columns)


def acb_right_multiply_numpy(matrix, right):
    right_ball = acb_mat(
        right.shape[0], right.shape[1],
        [acb(sf(right[row, column]))
         for row in range(right.shape[0]) for column in range(right.shape[1])],
    )
    return matrix * right_ball


def canonical_j_matrix(block, dimension):
    old = expanded_types(0, 30, dimension)
    internal = expanded_types(30, 65, dimension)
    new = expanded_types(65, 95, dimension)
    nd = 30 * dimension
    xd = 35 * dimension
    j_matrix = mp.matrix(xd + nd, xd + nd)
    k_xx = mp_submatrix(block, internal, internal)
    k_xn = mp_submatrix(block, internal, new)
    k_ox = mp_submatrix(block, old, internal)
    k_on = mp_submatrix(block, old, new)
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
    return j_matrix


def forced_response_balls(
    block, curvature_block, dimension, weak_positions, mapping, coefficient
):
    old = expanded_types(0, 30, dimension)
    internal = expanded_types(30, 65, dimension)
    new = expanded_types(65, 95, dimension)
    nd = 30 * dimension
    xd = 35 * dimension
    source_dimension = 5 * dimension
    j_matrix = canonical_j_matrix(block, dimension)
    j_ball = mp_matrix_to_acb(j_matrix)

    source = acb_zero_matrix(xd + nd, source_dimension)
    for pole_number, position in enumerate(weak_positions):
        for component in range(dimension):
            source[position * dimension + component,
                   pole_number * dimension + component] = mp_to_acb(coefficient)
    solved = j_ball.solve(-source)

    y_internal = acb_mat(xd, source_dimension, [
        solved[row, column]
        for row in range(xd) for column in range(source_dimension)
    ])
    y_new = acb_mat(nd, source_dimension, [
        solved[xd + row, column]
        for row in range(nd) for column in range(source_dimension)
    ])

    k_nx = mp_matrix_to_acb(mp_submatrix(block, new, internal))
    k_nn = mp_matrix_to_acb(mp_submatrix(block, new, new))
    p_new = k_nx * y_internal + k_nn * y_new
    raw_phase = acb_mat(2 * nd, source_dimension, [
        (y_new[row, column] if row < nd else p_new[row - nd, column])
        for row in range(2 * nd) for column in range(source_dimension)
    ])
    phase = acb_zero_matrix(2 * nd, source_dimension)
    for next_type, final_type in enumerate(mapping):
        for component in range(dimension):
            target = next_type * dimension + component
            source_row = final_type * dimension + component
            for column in range(source_dimension):
                phase[target, column] = raw_phase[source_row, column]
                phase[nd + target, column] = raw_phase[nd + source_row, column]

    slab = acb_zero_matrix(95 * dimension, source_dimension)
    for row in range(xd):
        for column in range(source_dimension):
            slab[30 * dimension + row, column] = y_internal[row, column]
    for row in range(nd):
        for column in range(source_dimension):
            slab[65 * dimension + row, column] = y_new[row, column]
    curvature = mp_matrix_to_acb(curvature_block) * slab

    return {
        "determinant": j_ball.det(),
        "j_matrix": j_matrix,
        "source": source,
        "solved": solved,
        "phase": phase,
        "curvature": curvature,
    }


def zero_sum_bases(trivial, source_dimension):
    if not trivial:
        identity = np.eye(source_dimension)
        return identity, identity
    differences = np.zeros((5, 4))
    for column in range(4):
        differences[column, column] = 1
        differences[4, column] = -1
    q_qr, _ = la.qr(differences, mode="economic")
    q_null = la.null_space(np.ones((1, 5)), rcond=None)
    return q_qr, q_null


def ball_mid_radii_after_domain(matrix, domain):
    reduced = acb_right_multiply_numpy(matrix, domain)
    midpoint, radii = acb_midpoint_and_radii(reduced)
    return reduced, midpoint, radii


def calibrated_rank(variant_data, field):
    matrices = {name: item[f"{field}_midpoint"] for name, item in variant_data.items()}
    radii = {name: item[f"{field}_radii"] for name, item in variant_data.items()}
    operational = matrices["operational_primary"]
    epsilon_step = (
        la.norm(operational - matrices["operational_shadow"], 2)
        + la.norm(matrices["validation_primary"] - matrices["validation_shadow"], 2)
        + la.norm(operational - matrices["validation_primary"], 2)
    )
    epsilon_ball = max(la.norm(value, "fro") for value in radii.values())
    left, singular, right_h = la.svd(
        operational, full_matrices=False, lapack_driver="gesdd"
    )
    singular_gesvd = la.svd(
        operational, compute_uv=False, lapack_driver="gesvd"
    )
    residual = max((
        np.linalg.norm(operational @ right_h[index].conj() - value * left[:, index])
        + np.linalg.norm(operational.conj().T @ left[:, index]
                         - value * right_h[index].conj())
        for index, value in enumerate(singular)
    ), default=0.0)
    epsilon_svd = float(residual + np.max(np.abs(singular - singular_gesvd)))
    epsilon = float(epsilon_step + epsilon_ball + epsilon_svd + 1e-70)
    other = {
        name: la.svdvals(matrix) for name, matrix in matrices.items()
        if name != "operational_primary"
    }
    resolved = singular > 100 * epsilon
    zero = np.array([
        max([singular[index]] + [value[index] for value in other.values()])
        < 10 * epsilon
        for index in range(len(singular))
    ], dtype=bool)
    opened = ~(resolved | zero)
    return {
        "singular_values": singular,
        "epsilon_step": epsilon_step,
        "epsilon_ball": epsilon_ball,
        "epsilon_svd": epsilon_svd,
        "epsilon": epsilon,
        "resolved_count": int(np.sum(resolved)),
        "zero_count": int(np.sum(zero)),
        "open_count": int(np.sum(opened)),
    }


def orthonormal_columns(matrix):
    q, _ = la.qr(matrix, mode="economic")
    singular = la.svdvals(matrix)
    return q, singular


def subspace_distance(left, right):
    q_left, s_left = orthonormal_columns(left)
    q_right, s_right = orthonormal_columns(right)
    overlap = la.svdvals(q_left.conj().T @ q_right)
    minimum = min(1.0, max(0.0, float(np.min(overlap))))
    return math.sqrt(max(0.0, 1 - minimum**2)), minimum, s_left, s_right


def extreme_subspace(matrix, count, branch):
    eigenvalues, eigenvectors = la.eig(matrix)
    moduli = np.abs(eigenvalues)
    if branch == "plus":
        order = np.argsort(moduli)[::-1]
        selected = order[:count]
        boundary = math.sqrt(moduli[order[count - 1]] * moduli[order[count]])
        selector = lambda value: abs(value) > boundary
        gap = float(moduli[order[count - 1]] / moduli[order[count]])
    else:
        order = np.argsort(moduli)
        selected = order[:count]
        boundary = math.sqrt(moduli[order[count - 1]] * moduli[order[count]])
        selector = lambda value: abs(value) < boundary
        gap = float(moduli[order[count]] / moduli[order[count - 1]])
    _, vectors, selected_count = la.schur(matrix, output="complex", sort=selector)
    basis = vectors[:, :selected_count]
    direct, _ = orthonormal_columns(eigenvectors[:, selected])
    unselected = np.setdiff1d(np.arange(len(eigenvalues)), selected)
    separation = float(np.min(np.abs(
        eigenvalues[selected, None] - eigenvalues[None, unselected]
    )))
    direct_distance, _, _, _ = subspace_distance(basis, direct)
    return {
        "basis": basis,
        "direct_basis": direct,
        "selected_count": int(selected_count),
        "gap": gap,
        "boundary": boundary,
        "direct_distance": direct_distance,
        "eigenvector_condition": float(np.linalg.cond(eigenvectors)),
        "spectral_separation": separation,
    }


def comparison_record(variant_data, branch):
    distances = {}
    minimum_overlaps = {}
    bounds = []
    direct_discrepancies = []
    conditions = []
    for name, data in variant_data.items():
        phase = data["phase_midpoint"]
        extreme = data[f"extreme_{branch}"]
        distance, overlap, phase_singular, _ = subspace_distance(
            phase, extreme["basis"]
        )
        direct_distance, _, _, _ = subspace_distance(
            phase, extreme["direct_basis"]
        )
        distances[name] = distance
        minimum_overlaps[name] = overlap
        direct_discrepancies.append(max(
            abs(distance - direct_distance), extreme["direct_distance"]
        ))
        phase_radius = la.norm(data["phase_radii"], "fro")
        tangent_radius = la.norm(data["tangent_radii"], "fro")
        bounds.append(
            phase_radius / max(1e-300, float(phase_singular[-1]))
            + tangent_radius * extreme["eigenvector_condition"]
              / max(1e-300, extreme["spectral_separation"])
        )
        conditions.append(max(
            1.0,
            float(phase_singular[0] / phase_singular[-1]),
            extreme["eigenvector_condition"],
        ))
    op = distances["operational_primary"]
    epsilon_step = (
        abs(op - distances["operational_shadow"])
        + abs(distances["validation_primary"] - distances["validation_shadow"])
        + abs(op - distances["validation_primary"])
    )
    epsilon_binary = max(direct_discrepancies) + (
        10 * np.finfo(float).eps * max(conditions)
    )
    epsilon_ball = max(bounds)
    epsilon = float(epsilon_step + epsilon_binary + epsilon_ball + 1e-70)
    if op <= 10 * epsilon:
        label = "IDENTIFIED"
    elif op > 100 * epsilon:
        label = "SEPARATED"
    else:
        label = "NUMERICALLY_OPEN"
    return {
        "distance": op,
        "minimum_overlap": minimum_overlaps["operational_primary"],
        "angle_degrees": math.degrees(math.asin(min(1.0, max(0.0, op)))),
        "epsilon_step": epsilon_step,
        "epsilon_binary": epsilon_binary,
        "epsilon_ball": epsilon_ball,
        "epsilon": epsilon,
        "label": label,
        "variant_distances": distances,
    }


def stored_rank_sector(sector, stored, used):
    target = complex(
        float(mp.re(sector["old_central_eigenvalue"])),
        float(mp.im(sector["old_central_eigenvalue"])),
    )
    choices = []
    for index, item in enumerate(stored):
        if index in used:
            continue
        center = complex(
            float(item["central_eigenvalue"]["real"]),
            float(item["central_eigenvalue"]["imaginary"]),
        )
        choices.append((abs(center - target), index))
    _, index = min(choices)
    used.add(index)
    return index, stored[index]


def serialize_rank(record):
    return {
        "singular_values": [sf(value) for value in record["singular_values"]],
        "epsilon_step": sf(record["epsilon_step"]),
        "epsilon_ball": sf(record["epsilon_ball"]),
        "epsilon_svd": sf(record["epsilon_svd"]),
        "epsilon": sf(record["epsilon"]),
        "resolved_count": record["resolved_count"],
        "zero_count": record["zero_count"],
        "open_count": record["open_count"],
    }


def serialize_comparison(record):
    return {
        "distance": sf(record["distance"]),
        "minimum_overlap": sf(record["minimum_overlap"]),
        "angle_degrees": sf(record["angle_degrees"]),
        "epsilon_step": sf(record["epsilon_step"]),
        "epsilon_binary": sf(record["epsilon_binary"]),
        "epsilon_ball": sf(record["epsilon_ball"]),
        "epsilon": sf(record["epsilon"]),
        "label": record["label"],
        "variant_distances": {
            name: sf(value) for name, value in record["variant_distances"].items()
        },
    }


print("=" * 78)

records = {}
all_comparisons = []
numeric_arrays = {}
global_controls = provenance_ok and gro.tests == gro.passed == 43
finite_difference_errors = []
zero_response_exact = True
all_gap_controls = True
restored_phase = Counter()
restored_curvature = Counter()

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing action and curvature kernels", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    mapping = tuple(gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"])
    index_data = group_and_index_data(model, state)
    geometry = prepare_geometry(model, index_data)
    weak_positions = [
        position for position in range(35)
        if index_data["edge_kind"][24 * (30 + position)] == "pole"
    ]
    incidence_ok = bool(
        len(weak_positions) == 5
        and weak_positions == schur_input["parities"][parity]["weak_orbit_positions"]
        and 24 * len(weak_positions) == 120
    )
    check(
        f"{parity}: source carrier is exactly the 120 incidence-selected pole rows",
        incidence_ok,
        f"weak positions={weak_positions}",
    )

    sectors, basis_control = high_precision_sector_bases(index_data)
    basis_ok = bool(
        sorted(sector["dimension"] for sector in sectors) == [1, 1, 1, 2, 2, 2, 3]
        and sum(5 * sector["dimension"] ** 2 for sector in sectors) == 120
        and all(
            value < mp.mpf("1e-70")
            for key, value in basis_control.items() if key.startswith("maximum_")
        )
    )
    check(
        f"{parity}: high-precision 2T sectors exhaust all 120 mass sources",
        basis_ok,
        str([sector["dimension"] for sector in sectors]),
    )

    s = mp.mpf(state[0])
    kind_values = {
        "old": L0_SQUARE,
        "internal": mp.exp(s) * L0_SQUARE - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * s) * L0_SQUARE,
    }
    pattern_cache, branch_control = high_precision_pattern_cache(
        geometry["patterns"], kind_values
    )
    kernels, kernel_control = assemble_full_representative_kernels(
        index_data, geometry, pattern_cache
    )
    curvature_data = triangle_response_data(
        model, index_data, geometry, pattern_cache
    )
    action_ok = bool(
        branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-70")
        and curvature_data["orbit_labels_constant"]
        and len(curvature_data["internal_global_types"]) == 160
        and max(curvature_data["maximum_equivariance_residual"].values())
            < mp.mpf("1e-70")
    )
    check(
        f"{parity}: branch, action-Hessian and curvature-equivariance controls pass",
        action_ok,
        "kernel imag=" + mp.nstr(kernel_control["maximum_imaginary"], 5)
        + ", curvature residual="
        + mp.nstr(max(curvature_data["maximum_equivariance_residual"].values()), 5),
    )

    m0 = MASS / 120
    coefficient = -4 * mp.pi * m0 * mp.sqrt(index_data["rho"])
    for step in (mp.mpf("1e-20"), mp.mpf("3e-20")):
        plus = -4 * mp.pi * m0 * (1 + step) * mp.sqrt(index_data["rho"])
        minus = -4 * mp.pi * m0 * (1 - step) * mp.sqrt(index_data["rho"])
        estimate = (plus - minus) / (2 * step)
        finite_difference_errors.append(abs(estimate - coefficient))
    source_formula_ok = bool(
        coefficient != 0 and max(finite_difference_errors) < mp.mpf("1e-70")
    )
    check(
        f"{parity}: analytic relative-mass source passes two finite differences",
        source_formula_ok,
        "max error=" + mp.nstr(max(finite_difference_errors), 5),
    )

    used_rank = set()
    sector_records = []
    determinant_ok = True
    reproduction_errors = []
    source_ranks_ok = True
    zero_sum_control_distances = []
    uniform_control_ranks = []
    selection_ok = True

    for sector_index, sector in enumerate(sectors):
        dimension = sector["dimension"]
        trivial = bool(sector["constant_overlap"] > 1 - mp.mpf("1e-70"))
        expected_count = 4 if trivial else 5 * dimension
        domain, domain_alt = zero_sum_bases(trivial, 5 * dimension)
        print(
            f"[{parity}] sector {sector_index + 1}/7 d={dimension}: "
            f"Flint mass response, retained={expected_count}",
            flush=True,
        )

        blocks = {
            name: project_full_kernel(kernel, sector)
            for name, kernel in kernels.items()
        }
        curvature_blocks = {
            name: project_curvature_kernel(kernel, sector)
            for name, kernel in curvature_data["kernels"].items()
        }
        variant_data = {}
        for name in VARIANTS:
            response = forced_response_balls(
                blocks[name], curvature_blocks[name], dimension,
                weak_positions, mapping, coefficient,
            )
            determinant_ok &= not response["determinant"].contains(0)
            source_mid, _ = acb_midpoint_and_radii(response["source"])
            source_ranks_ok &= np.linalg.matrix_rank(source_mid) == 5 * dimension

            zero_source = acb_zero_matrix(65 * dimension, 5 * dimension)
            zero_solved = mp_matrix_to_acb(response["j_matrix"]).solve(zero_source)
            zero_mid, zero_rad = acb_midpoint_and_radii(zero_solved)
            zero_response_exact &= bool(
                np.max(np.abs(zero_mid)) == 0 and np.max(zero_rad) == 0
            )

            _, phase_mid, phase_rad = ball_mid_radii_after_domain(
                response["phase"], domain
            )
            _, curvature_mid, curvature_rad = ball_mid_radii_after_domain(
                response["curvature"], domain
            )
            _, solved_mid, solved_rad = ball_mid_radii_after_domain(
                response["solved"], domain
            )
            _, phase_alt, _ = ball_mid_radii_after_domain(
                response["phase"], domain_alt
            )
            alt_distance, _, _, _ = subspace_distance(phase_mid, phase_alt)
            zero_sum_control_distances.append(alt_distance)

            tangent = numeric[
                f"{parity}_sector{sector_index}_{name}_tangent_midpoint"
            ]
            tangent_radii = numeric[
                f"{parity}_sector{sector_index}_{name}_tangent_radii"
            ]
            extreme_plus = extreme_subspace(tangent, expected_count, "plus")
            extreme_minus = extreme_subspace(tangent, expected_count, "minus")
            selection_ok &= bool(
                extreme_plus["selected_count"] == expected_count
                and extreme_minus["selected_count"] == expected_count
            )
            all_gap_controls &= bool(
                extreme_plus["gap"] > 2 and extreme_minus["gap"] > 2
            )
            variant_data[name] = {
                "phase_midpoint": phase_mid,
                "phase_radii": phase_rad,
                "curvature_midpoint": curvature_mid,
                "curvature_radii": curvature_rad,
                "solved_midpoint": solved_mid,
                "solved_radii": solved_rad,
                "tangent_radii": tangent_radii,
                "extreme_plus": extreme_plus,
                "extreme_minus": extreme_minus,
            }

            if trivial:
                uniform = np.ones((5, 1)) / math.sqrt(5)
                _, uniform_mid, _ = ball_mid_radii_after_domain(
                    response["phase"], uniform
                )
                uniform_control_ranks.append(int(np.linalg.matrix_rank(uniform_mid)))

            prefix = f"{parity}_sector{sector_index}_{name}"
            numeric_arrays[f"{prefix}_phase_midpoint"] = phase_mid
            numeric_arrays[f"{prefix}_phase_radii"] = phase_rad
            numeric_arrays[f"{prefix}_curvature_midpoint"] = curvature_mid
            numeric_arrays[f"{prefix}_curvature_radii"] = curvature_rad

        phase_rank = calibrated_rank(variant_data, "phase")
        curvature_rank = calibrated_rank(variant_data, "curvature")
        solved_rank = calibrated_rank(variant_data, "solved")
        for label, count in (
            ("resolved", phase_rank["resolved_count"]),
            ("zero", phase_rank["zero_count"]),
            ("open", phase_rank["open_count"]),
        ):
            restored_phase[label] += count * dimension
        for label, count in (
            ("resolved", curvature_rank["resolved_count"]),
            ("zero", curvature_rank["zero_count"]),
            ("open", curvature_rank["open_count"]),
        ):
            restored_curvature[label] += count * dimension

        stored_index, stored_sector = stored_rank_sector(
            sector, rank_input["parities"][parity]["sectors"], used_rank
        )
        current_j_singular = la.svdvals(mp_to_numpy(
            canonical_j_matrix(blocks["operational_primary"], dimension)
        ))
        stored_j_singular = np.asarray([
            float(value) for value in stored_sector["singular_values"]
        ])
        reproduction_error = float(np.max(
            np.abs(current_j_singular - stored_j_singular)
            / np.maximum(1.0, np.abs(stored_j_singular))
        ))
        reproduction_errors.append(reproduction_error)

        comparisons = {}
        for branch in ("plus", "minus"):
            comparison = comparison_record(variant_data, branch)
            comparisons[branch] = comparison
            all_comparisons.append({
                "parity": parity,
                "sector_index": sector_index,
                "dimension": dimension,
                "trivial": trivial,
                "branch": branch,
                **comparison,
            })

        sector_records.append({
            "sector_index": sector_index,
            "stored_rank_sector_index": stored_index,
            "dimension": dimension,
            "trivial": trivial,
            "retained_mass_dimension": expected_count,
            "j_spectrum_reproduction_error": reproduction_error,
            "zero_sum_basis_distance_maximum": max(
                zero_sum_control_distances[-len(VARIANTS):]
            ),
            "minimum_plus_gap": min(
                data["extreme_plus"]["gap"] for data in variant_data.values()
            ),
            "minimum_minus_gap": min(
                data["extreme_minus"]["gap"] for data in variant_data.values()
            ),
            "solved_rank": serialize_rank(solved_rank),
            "phase_rank": serialize_rank(phase_rank),
            "curvature_rank": serialize_rank(curvature_rank),
            "comparisons": {
                name: serialize_comparison(value)
                for name, value in comparisons.items()
            },
        })

    reproduction_ok = bool(
        len(used_rank) == 7 and max(reproduction_errors) < 2e-10
    )
    basis_choice_ok = bool(max(zero_sum_control_distances) < 2e-7)
    uniform_ok = bool(uniform_control_ranks == [1] * len(VARIANTS))
    check(
        f"{parity}: every J determinant, source rank and frozen spectrum control passes",
        determinant_ok and source_ranks_ok and reproduction_ok,
        f"max spectrum error={max(reproduction_errors):.3e}",
    )
    check(
        f"{parity}: zero source, uniform source and alternate zero-sum basis controls pass",
        zero_response_exact and uniform_ok and basis_choice_ok,
        f"basis distance={max(zero_sum_control_distances):.3e}",
    )
    check(
        f"{parity}: all strong branch selections have fixed counts and gap above two",
        selection_ok and all_gap_controls,
        "minimum plus/minus="
        f"{min(item['minimum_plus_gap'] for item in sector_records):.3e}/"
        f"{min(item['minimum_minus_gap'] for item in sector_records):.3e}",
    )
    controls_ok = bool(
        incidence_ok and basis_ok and action_ok and source_formula_ok
        and determinant_ok and source_ranks_ok and reproduction_ok
        and zero_response_exact and uniform_ok and basis_choice_ok and selection_ok
    )
    global_controls &= controls_ok
    records[parity] = {
        "controls_ok": controls_ok,
        "weak_positions": weak_positions,
        "source_coefficient": coefficient,
        "maximum_j_spectrum_reproduction_error": max(reproduction_errors),
        "maximum_zero_sum_basis_distance": max(zero_sum_control_distances),
        "sectors": sector_records,
    }


label_counts = Counter(item["label"] for item in all_comparisons)
branch_counts = {
    branch: Counter(
        item["label"] for item in all_comparisons if item["branch"] == branch
    )
    for branch in ("plus", "minus")
}

phase_full = bool(
    restored_phase["resolved"] == 2 * 119
    and restored_phase["zero"] == restored_phase["open"] == 0
)
curvature_full = bool(
    restored_curvature["resolved"] == 2 * 119
    and restored_curvature["zero"] == restored_curvature["open"] == 0
)
all_plus_identified = branch_counts["plus"]["IDENTIFIED"] == 14
all_minus_identified = branch_counts["minus"]["IDENTIFIED"] == 14
all_plus_separated = branch_counts["plus"]["SEPARATED"] == 14
all_minus_separated = branch_counts["minus"]["SEPARATED"] == 14


def schedule_record(even, odd):
    difference = float(np.max(np.abs(
        np.asarray(even["singular_values"]) - np.asarray(odd["singular_values"])
    )))
    epsilon = float(even["epsilon"] + odd["epsilon"] + 1e-70)
    if difference <= 10 * epsilon:
        label = "SCHEDULE_ROBUST"
    elif difference > 100 * epsilon:
        label = "SCHEDULE_SEPARATED"
    else:
        label = "SCHEDULE_OPEN"
    return {"distance": difference, "epsilon": epsilon, "label": label}


schedule_records = []
used_odd = set()
for even_sector in records["even"]["sectors"]:
    choices = [
        (abs(even_sector["dimension"] - item["dimension"]), index)
        for index, item in enumerate(records["odd"]["sectors"])
        if index not in used_odd
        and item["stored_rank_sector_index"] == even_sector["stored_rank_sector_index"]
    ]
    _, odd_index = min(choices)
    used_odd.add(odd_index)
    odd_sector = records["odd"]["sectors"][odd_index]
    for field in ("phase_rank", "curvature_rank"):
        even_rank = {
            "singular_values": np.asarray([
                float(value) for value in even_sector[field]["singular_values"]
            ]),
            "epsilon": float(even_sector[field]["epsilon"]),
        }
        odd_rank = {
            "singular_values": np.asarray([
                float(value) for value in odd_sector[field]["singular_values"]
            ]),
            "epsilon": float(odd_sector[field]["epsilon"]),
        }
        schedule_records.append({
            "stored_rank_sector_index": even_sector["stored_rank_sector_index"],
            "dimension": even_sector["dimension"],
            "field": field,
            **schedule_record(even_rank, odd_rank),
        })

schedule_counts = Counter(item["label"] for item in schedule_records)

if not global_controls:
    outcome = "DUST_MASS_RESPONSE_CONTROL_FAILED"
elif not all_gap_controls or not phase_full or not curvature_full:
    outcome = "DUST_MASS_RESPONSE_MIXED_OR_OPEN"
elif all_plus_identified and all_minus_separated:
    outcome = "DUST_MASS_RESPONSE_EXPANDING_IDENTIFIED"
elif all_minus_identified and all_plus_separated:
    outcome = "DUST_MASS_RESPONSE_CONTRACTING_IDENTIFIED"
elif label_counts["SEPARATED"] == 28:
    outcome = "DUST_MASS_RESPONSE_BOTH_BRANCHES_SEPARATED"
else:
    outcome = "DUST_MASS_RESPONSE_MIXED_OR_OPEN"

check(
    "all 28 target-disclosed comparisons receive mechanical labels",
    len(all_comparisons) == 28 and sum(label_counts.values()) == 28,
    str(dict(label_counts)),
)
check(
    "the zero-total-mass phase and curvature rank ledgers restore 119 per schedule",
    phase_full and curvature_full,
    f"phase={dict(restored_phase)}, curvature={dict(restored_curvature)}",
)
check(
    "all 14 schedule comparisons receive calibrated labels",
    len(schedule_records) == 14 and sum(schedule_counts.values()) == 14,
    str(dict(schedule_counts)),
)
check(
    "the preregistered hierarchy assigns exactly one mass-response outcome",
    outcome in {
        "DUST_MASS_RESPONSE_CONTROL_FAILED",
        "DUST_MASS_RESPONSE_EXPANDING_IDENTIFIED",
        "DUST_MASS_RESPONSE_CONTRACTING_IDENTIFIED",
        "DUST_MASS_RESPONSE_BOTH_BRANCHES_SEPARATED",
        "DUST_MASS_RESPONSE_MIXED_OR_OPEN",
    },
    outcome,
)

np.savez_compressed(NUMERIC_OUTPUT, **numeric_arrays)
numeric_sha = sha256(NUMERIC_OUTPUT)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "classification": "DERIVED_COMPUTATIONAL_SOURCE_RESPONSE",
    "matter_status": "CONSERVED_SOURCE_PARAMETERS_NOT_CANONICAL_DUST_FIELDS",
    "target_disclosed": True,
    "known_119_count_disclosed": True,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "planck_target_parsed": False,
    "particle_target_parsed": False,
    "mass_source_dimension": 120,
    "zero_total_mass_dimension": 119,
    "comparison_attempts": 28,
    "comparison_label_counts": dict(sorted(label_counts.items())),
    "branch_label_counts": {
        branch: dict(sorted(counts.items())) for branch, counts in branch_counts.items()
    },
    "restored_phase_rank_ledger": dict(sorted(restored_phase.items())),
    "restored_curvature_rank_ledger": dict(sorted(restored_curvature.items())),
    "schedule_label_counts": dict(sorted(schedule_counts.items())),
    "schedule_comparisons": [{
        **item,
        "distance": sf(item["distance"]),
        "epsilon": sf(item["epsilon"]),
    } for item in schedule_records],
    "parities": {
        parity: {
            "controls_ok": record["controls_ok"],
            "weak_positions": record["weak_positions"],
            "source_coefficient": mp.nstr(record["source_coefficient"], 80),
            "maximum_j_spectrum_reproduction_error": sf(
                record["maximum_j_spectrum_reproduction_error"]
            ),
            "maximum_zero_sum_basis_distance": sf(
                record["maximum_zero_sum_basis_distance"]
            ),
            "sectors": record["sectors"],
        }
        for parity, record in records.items()
    },
    "numeric_archive": NUMERIC_OUTPUT.name,
    "numeric_archive_arrays": len(numeric_arrays),
    "numeric_archive_sha256": numeric_sha,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("=" * 78)
print(
    f"Comparisons: {dict(label_counts)}; "
    f"phase={dict(restored_phase)}; curvature={dict(restored_curvature)}",
    flush=True,
)
print(f"Schedules: {dict(schedule_counts)}", flush=True)
print(f"Artifact SHA-256: {sha256(OUTPUT)}", flush=True)
print(f"Numeric SHA-256: {numeric_sha}", flush=True)
print(f"Outcome: {outcome}", flush=True)
print(f"{passed}/{tests} checks passed", flush=True)
sys.exit(0 if passed == tests else 1)

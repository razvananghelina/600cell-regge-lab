#!/usr/bin/env python3
"""Confirmatory extreme-tangent versus lapse-subspace alignment audit.

Prior-art commit: 24d2ce6.
Protocol commit: e50a0ea.
The 119+1 count was disclosed before this comparison; no continuum target is loaded.
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
TANGENT_INPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
SCHUR_INPUT = HERE / "gravity_600cell_dust_full_lapse_schur.json"
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
OUTPUT = HERE / "gravity_600cell_dust_hyperbolic_lapse_alignment.json"

PRIOR_ART_COMMIT = "24d2ce6"
PROTOCOL_COMMIT = "e50a0ea"
EXPECTED_HASHES = {
    "tangent": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "schur": "4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349",
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


hashes = {
    "tangent": sha256(TANGENT_INPUT),
    "tangent_numeric": sha256(TANGENT_NUMERIC),
    "tangent_source": sha256(TANGENT_SOURCE),
    "schur": sha256(SCHUR_INPUT),
}
tangent_input = json.loads(TANGENT_INPUT.read_text())
schur_input = json.loads(SCHUR_INPUT.read_text())
tick = json.loads(TICK_INPUT.read_text())
numeric = np.load(TANGENT_NUMERIC)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and tangent_input["outcome"] == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and tangent_input["passed"] == tangent_input["tests"] == 19
    and tangent_input["numeric_archive_arrays"] == len(numeric.files) == 224
    and tangent_input["numeric_archive_sha256"] == hashes["tangent_numeric"]
    and schur_input["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
    and all(
        record["resolved_schur_rank"] == 120
        and record["schur_zero_count"] == 0
        and record["schur_open_count"] == 0
        for record in schur_input["parities"].values()
    )
)
check("the confirmatory comparison has exact frozen provenance", provenance_ok, str(hashes))


spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_hyperbolic_lapse", GEOMETRY_SOURCE
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
    "the directly imported slab geometry retains all 43 controls",
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

models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}


def response_and_lift_ball(block, dimension, weak_type_positions):
    old = expanded_types(0, 30, dimension)
    internal = expanded_types(30, 65, dimension)
    new = expanded_types(65, 95, dimension)
    nd = 30 * dimension
    xd = 35 * dimension
    size = xd + nd

    k_xx = mp_submatrix(block, internal, internal)
    k_xn = mp_submatrix(block, internal, new)
    k_ox = mp_submatrix(block, old, internal)
    k_on = mp_submatrix(block, old, new)
    j_matrix = mp.matrix(size, size)
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

    rhs = mp.matrix(size, 2 * nd)
    k_xo = mp_submatrix(block, internal, old)
    k_oo = mp_submatrix(block, old, old)
    for row in range(xd):
        for column in range(nd):
            rhs[row, column] = -k_xo[row, column]
    for row in range(nd):
        for column in range(nd):
            rhs[xd + row, column] = k_oo[row, column]
        rhs[xd + row, nd + row] = 1

    j_ball = mp_matrix_to_acb(j_matrix)
    det_j = j_ball.det()
    response_ball = j_ball.solve(mp_matrix_to_acb(rhs))
    response_midpoint, response_radii = acb_midpoint_and_radii(response_ball)

    weak = [
        position * dimension + component
        for position in weak_type_positions for component in range(dimension)
    ]
    strong = [index for index in range(size) if index not in set(weak)]
    a_mp = mp_submatrix(j_matrix, strong, strong)
    b_mp = mp_submatrix(j_matrix, strong, weak)
    solve_ball = mp_matrix_to_acb(a_mp).solve(mp_matrix_to_acb(b_mp))
    ordered_lift = acb_mat(size, len(weak))
    for row in range(len(strong)):
        for column in range(len(weak)):
            ordered_lift[row, column] = -solve_ball[row, column]
    for index in range(len(weak)):
        ordered_lift[len(strong) + index, index] = 1
    order = strong + weak
    lift_ball = acb_mat(size, len(weak))
    for ordered_row, original_row in enumerate(order):
        for column in range(len(weak)):
            lift_ball[original_row, column] = ordered_lift[ordered_row, column]
    lift_midpoint, lift_radii = acb_midpoint_and_radii(lift_ball)
    return {
        "det_j": det_j,
        "response_midpoint": response_midpoint,
        "response_radii": response_radii,
        "lift_midpoint": lift_midpoint,
        "lift_radii": lift_radii,
        "weak": weak,
        "strong": strong,
    }


def geometric_lapse_matrix(index_data, weak_positions, state):
    diagonal = mp.exp(mp.mpf(state[0])) * L0_SQUARE - index_data["rho"]
    coefficient = float(-index_data["rho"] / diagonal)
    weak_orbit_types = [30 + position for position in weak_positions]
    pole_edges = []
    for orbit_type in weak_orbit_types:
        pole_edges.extend(index_data["orbit_edges"][orbit_type])
    pole_to_column = {edge: index for index, edge in enumerate(pole_edges)}
    matrix = np.zeros((1560, 120), dtype=float)
    column_position = {
        orbit_type: position
        for position, orbit_type in enumerate(tuple(range(30, 65)) + tuple(range(65, 95)))
    }
    for pole, column in pole_to_column.items():
        pole_index = index_data["edge_to_index"][pole]
        pole_type, pole_group = divmod(pole_index, 24)
        matrix[24 * column_position[pole_type] + pole_group, column] = 1.0
        new_vertex = pole[1]
        for edge, global_index in index_data["edge_to_index"].items():
            if index_data["edge_kind"][global_index] == "internal" and edge[1] == new_vertex:
                orbit_type, group = divmod(global_index, 24)
                matrix[24 * column_position[orbit_type] + group, column] = coefficient
    return matrix, coefficient


def projected_geometric_lapse(geometry_matrix, sector):
    dimension = sector["dimension"]
    basis = mp_to_numpy(sector["basis"])
    full_basis = np.kron(np.eye(65), basis)
    pole_basis = np.kron(np.eye(5), basis)
    return full_basis.conj().T @ geometry_matrix @ pole_basis


def orthonormal_columns(matrix):
    q, r = la.qr(matrix, mode="economic")
    singular = la.svd(matrix, compute_uv=False, lapack_driver="gesvd")
    return q, singular


def subspace_distance(left, right):
    q_left, singular_left = orthonormal_columns(left)
    q_right, singular_right = orthonormal_columns(right)
    overlap = la.svd(q_left.conj().T @ q_right, compute_uv=False, lapack_driver="gesvd")
    minimum = min(1.0, max(0.0, float(np.min(overlap))))
    distance = math.sqrt(max(0.0, 1 - minimum**2))
    return distance, minimum, singular_left, singular_right


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
    _, schur_vectors, selected_count = la.schur(
        matrix, output="complex", sort=selector
    )
    schur_basis = schur_vectors[:, :selected_count]
    direct_basis, _ = orthonormal_columns(eigenvectors[:, selected])
    direct_distance, _, _, _ = subspace_distance(schur_basis, direct_basis)
    unselected = np.setdiff1d(np.arange(len(eigenvalues)), selected)
    spectral_separation = float(np.min(np.abs(
        eigenvalues[selected, None] - eigenvalues[None, unselected]
    )))
    return {
        "basis": schur_basis,
        "direct_basis": direct_basis,
        "selected_count": int(selected_count),
        "gap": gap,
        "boundary": boundary,
        "direct_distance": direct_distance,
        "eigenvector_condition": float(np.linalg.cond(eigenvectors)),
        "spectral_separation": spectral_separation,
    }


def comparison_record(variant_data, branch, candidate):
    distances = {}
    overlaps = {}
    binary_discrepancies = []
    ball_bounds = []
    condition_bounds = []
    for name, data in variant_data.items():
        response = data["response_midpoint"]
        extreme = data[f"extreme_{branch}"]
        transported = response @ extreme["basis"]
        transported_direct = response @ extreme["direct_basis"]
        target = data["lift_midpoint"] if candidate == "canonical" else data["geometric"]
        distance, overlap, singular_transport, singular_target = subspace_distance(
            transported, target
        )
        distances[name] = distance
        overlaps[name] = overlap
        binary_distance, _, _, _ = subspace_distance(
            transported, transported_direct
        )
        binary_discrepancies.append(max(binary_distance, extreme["direct_distance"]))
        response_radius = la.norm(data["response_radii"], "fro")
        tangent_radius = la.norm(data["tangent_radii"], "fro")
        target_radius = (
            la.norm(data["lift_radii"], "fro") if candidate == "canonical" else 0.0
        )
        ball_bounds.append(
            response_radius / max(1e-300, float(singular_transport[-1]))
            + target_radius / max(1e-300, float(singular_target[-1]))
            + tangent_radius * extreme["eigenvector_condition"]
              / max(1e-300, extreme["spectral_separation"])
        )
        condition_bounds.append(max(
            1.0,
            float(singular_transport[0] / singular_transport[-1]),
            float(singular_target[0] / singular_target[-1]),
            extreme["eigenvector_condition"],
        ))
    op = distances["operational_primary"]
    epsilon_step = (
        abs(op - distances["operational_shadow"])
        + abs(distances["validation_primary"] - distances["validation_shadow"])
        + abs(op - distances["validation_primary"])
    )
    epsilon_binary = max(binary_discrepancies) + (
        10 * np.finfo(float).eps * max(condition_bounds)
    )
    epsilon_ball = max(ball_bounds)
    epsilon = epsilon_step + epsilon_binary + epsilon_ball + 1e-70
    if op <= 10 * epsilon:
        label = "IDENTIFIED"
    elif op > 100 * epsilon:
        label = "SEPARATED"
    else:
        label = "NUMERICALLY_OPEN"
    return {
        "distance": op,
        "minimum_overlap": overlaps["operational_primary"],
        "angle_degrees": math.degrees(math.asin(min(1.0, max(0.0, op)))),
        "epsilon_step": epsilon_step,
        "epsilon_binary": epsilon_binary,
        "epsilon_ball": epsilon_ball,
        "epsilon_distance": epsilon,
        "label": label,
        "variant_distances": distances,
    }


def stored_sector_match(sector, stored, used):
    target = complex(
        float(mp.re(sector["old_central_eigenvalue"])),
        float(mp.im(sector["old_central_eigenvalue"])),
    )
    choices = []
    for index, item in enumerate(stored):
        if index in used:
            continue
        center = complex(
            float(item["old_central_eigenvalue"]["real"]),
            float(item["old_central_eigenvalue"]["imaginary"]),
        )
        choices.append((abs(center - target), index))
    _, index = min(choices)
    used.add(index)
    return index, stored[index]


def sf(value):
    return f"{float(value):.17e}"


print("=" * 78)
print("CONFIRMATORY HYPERBOLIC--LAPSE SUBSPACE ALIGNMENT")
print("=" * 78)

records = {}
global_controls = provenance_ok and gro.tests == gro.passed == 43
all_extreme_controls = True
all_comparisons = []

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing response kernels", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    index_data = group_and_index_data(model, state)
    geometry = prepare_geometry(model, index_data)
    weak_positions = [
        position for position in range(35)
        if index_data["edge_kind"][24 * (30 + position)] == "pole"
    ]
    stored_weak = schur_input["parities"][parity]["weak_orbit_positions"]
    weak_ok = bool(weak_positions == stored_weak and len(weak_positions) == 5)
    check(
        f"{parity}: five pole orbits are selected before spectral vectors",
        weak_ok,
        f"positions={weak_positions}",
    )

    sectors, sector_control = high_precision_sector_bases(index_data)
    basis_ok = bool(
        sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and all(
            value < mp.mpf("1e-70")
            for key, value in sector_control.items() if key.startswith("maximum_")
        )
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
    kernel_ok = bool(
        branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-70")
    )
    check(
        f"{parity}: response reconstruction retains basis, branch and reality controls",
        basis_ok and kernel_ok,
        f"kernel imag={mp.nstr(kernel_control['maximum_imaginary'], 5)}",
    )

    geometry_lapse, geometry_coefficient = geometric_lapse_matrix(
        index_data, weak_positions, state
    )
    used_stored = set()
    sector_records = []
    determinant_ok = True
    rank_ok = True
    reproduction_errors = []
    gaps_ok = True

    for sector_index, sector in enumerate(sectors):
        dimension = sector["dimension"]
        count = 5 * dimension
        print(
            f"[{parity}] sector {sector_index + 1}/7 d={dimension}: "
            "Flint responses and ordered Schur spaces",
            flush=True,
        )
        blocks = {
            name: project_full_kernel(kernel, sector)
            for name, kernel in kernels.items()
        }
        geometric = projected_geometric_lapse(geometry_lapse, sector)
        variant_data = {}
        for name, block in blocks.items():
            response = response_and_lift_ball(
                block, dimension, weak_positions
            )
            determinant_ok &= not response["det_j"].contains(0)
            tangent = numeric[
                f"{parity}_sector{sector_index}_{name}_tangent_midpoint"
            ]
            tangent_radii = numeric[
                f"{parity}_sector{sector_index}_{name}_tangent_radii"
            ]
            extreme_plus = extreme_subspace(tangent, count, "plus")
            extreme_minus = extreme_subspace(tangent, count, "minus")
            gaps_ok &= bool(
                extreme_plus["selected_count"] == count
                and extreme_minus["selected_count"] == count
                and extreme_plus["gap"] > 2
                and extreme_minus["gap"] > 2
            )
            _, response_singular = orthonormal_columns(response["response_midpoint"])
            _, lift_singular = orthonormal_columns(response["lift_midpoint"])
            _, geometric_singular = orthonormal_columns(geometric)
            rank_ok &= bool(
                response_singular[-1] > 1e-12
                and lift_singular[-1] > 1e-12
                and geometric_singular[-1] > 1e-12
            )
            variant_data[name] = {
                **response,
                "geometric": geometric,
                "tangent_radii": tangent_radii,
                "extreme_plus": extreme_plus,
                "extreme_minus": extreme_minus,
            }

        stored_index, stored_sector = stored_sector_match(
            sector, schur_input["parities"][parity]["sectors"], used_stored
        )
        current_distance, _, _, _ = subspace_distance(
            variant_data["operational_primary"]["lift_midpoint"], geometric
        )
        stored_distance = float(
            stored_sector["subspaces"]["canonical_vs_geometric"]["projector_distance"]
        )
        reproduction_error = abs(current_distance - stored_distance)
        reproduction_errors.append(reproduction_error)

        comparisons = {}
        for branch in ("plus", "minus"):
            for candidate in ("canonical", "geometric"):
                key = f"{branch}_vs_{candidate}"
                comparison = comparison_record(variant_data, branch, candidate)
                comparisons[key] = comparison
                all_comparisons.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "dimension": dimension,
                    "branch": branch,
                    "candidate": candidate,
                    **comparison,
                })

        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "stored_sector_index": stored_index,
            "constant_overlap": sector["constant_overlap"],
            "geometry_coefficient": geometry_coefficient,
            "canonical_geometric_distance": current_distance,
            "canonical_geometric_reproduction_error": reproduction_error,
            "plus_gap_minimum": min(
                data["extreme_plus"]["gap"] for data in variant_data.values()
            ),
            "minus_gap_minimum": min(
                data["extreme_minus"]["gap"] for data in variant_data.values()
            ),
            "plus_direct_schur_maximum": max(
                data["extreme_plus"]["direct_distance"] for data in variant_data.values()
            ),
            "minus_direct_schur_maximum": max(
                data["extreme_minus"]["direct_distance"] for data in variant_data.values()
            ),
            "comparisons": comparisons,
        })

    reproduction_ok = bool(
        len(used_stored) == 7 and max(reproduction_errors) < 2e-8
    )
    check(
        f"{parity}: all response/lapse ranks and Flint determinants pass",
        determinant_ok and rank_ok,
    )
    check(
        f"{parity}: prior canonical/geometric distances reproduce independently",
        reproduction_ok,
        f"max error={max(reproduction_errors):.3e}",
    )
    check(
        f"{parity}: all 56 extreme selections have fixed count and gap above two",
        gaps_ok,
        f"minimum gaps plus/minus="
        f"{min(item['plus_gap_minimum'] for item in sector_records):.3e}/"
        f"{min(item['minus_gap_minimum'] for item in sector_records):.3e}",
    )
    all_extreme_controls &= gaps_ok
    controls_ok = bool(
        weak_ok and basis_ok and kernel_ok and determinant_ok and rank_ok and reproduction_ok
    )
    global_controls &= controls_ok
    records[parity] = {
        "controls_ok": controls_ok,
        "weak_positions": weak_positions,
        "sectors": sector_records,
    }


label_counts = Counter(item["label"] for item in all_comparisons)
hit_fractions = {}
for branch in ("plus", "minus"):
    for candidate in ("canonical", "geometric"):
        selected = [
            item for item in all_comparisons
            if item["branch"] == branch and item["candidate"] == candidate
        ]
        hit_fractions[f"{branch}_vs_{candidate}"] = Counter(
            item["label"] for item in selected
        )


def all_identified(branch, candidate):
    selected = [
        item for item in all_comparisons
        if item["branch"] == branch and item["candidate"] == candidate
    ]
    return len(selected) == 14 and all(item["label"] == "IDENTIFIED" for item in selected)


geometric_branch = next((
    branch for branch in ("plus", "minus")
    if all_identified(branch, "canonical") and all_identified(branch, "geometric")
), None)
canonical_branch = next((
    branch for branch in ("plus", "minus")
    if all_identified(branch, "canonical")
), None)

if not global_controls:
    outcome = "HYPERBOLIC_LAPSE_ALIGNMENT_CONTROL_FAILED"
elif not all_extreme_controls:
    outcome = "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
elif geometric_branch is not None:
    outcome = "HYPERBOLIC_GEOMETRIC_LAPSE_IDENTIFIED"
elif canonical_branch is not None:
    outcome = "HYPERBOLIC_CANONICAL_WEAK_IDENTIFIED_ONLY"
elif label_counts["SEPARATED"] == 56:
    outcome = "HYPERBOLIC_LAPSE_ALIGNMENT_REFUTED"
else:
    outcome = "HYPERBOLIC_LAPSE_ALIGNMENT_MIXED_OR_OPEN"

check(
    "all 56 preregistered comparisons receive mechanical labels",
    len(all_comparisons) == 56 and sum(label_counts.values()) == 56,
    str(dict(label_counts)),
)
check(
    "the frozen hierarchy assigns the hyperbolic-lapse outcome",
    outcome in {
        "HYPERBOLIC_LAPSE_ALIGNMENT_CONTROL_FAILED",
        "HYPERBOLIC_EXTREME_SUBSPACE_OPEN",
        "HYPERBOLIC_GEOMETRIC_LAPSE_IDENTIFIED",
        "HYPERBOLIC_CANONICAL_WEAK_IDENTIFIED_ONLY",
        "HYPERBOLIC_LAPSE_ALIGNMENT_REFUTED",
        "HYPERBOLIC_LAPSE_ALIGNMENT_MIXED_OR_OPEN",
    },
    f"outcome={outcome}",
)


def public_comparison(item):
    return {
        "distance": sf(item["distance"]),
        "minimum_overlap": sf(item["minimum_overlap"]),
        "angle_degrees": sf(item["angle_degrees"]),
        "epsilon_step": sf(item["epsilon_step"]),
        "epsilon_binary": sf(item["epsilon_binary"]),
        "epsilon_ball": sf(item["epsilon_ball"]),
        "epsilon_distance": sf(item["epsilon_distance"]),
        "label": item["label"],
        "variant_distances": {
            name: sf(value) for name, value in item["variant_distances"].items()
        },
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "post_observed_count_match_disclosed": True,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "outcome": outcome,
    "identified_global_branch": geometric_branch or canonical_branch,
    "label_counts": dict(label_counts),
    "hit_fractions": {
        key: {"denominator": 14, **dict(counts)}
        for key, counts in hit_fractions.items()
    },
    "parities": {
        parity: {
            "controls_ok": record["controls_ok"],
            "weak_orbit_positions": record["weak_positions"],
            "sectors": [
                {
                    "sector_index": sector["sector_index"],
                    "irrep_dimension": sector["dimension"],
                    "stored_sector_index": sector["stored_sector_index"],
                    "constant_overlap": mp.nstr(sector["constant_overlap"], 70),
                    "geometry_coefficient": sf(sector["geometry_coefficient"]),
                    "canonical_geometric_distance": sf(
                        sector["canonical_geometric_distance"]
                    ),
                    "canonical_geometric_reproduction_error": sf(
                        sector["canonical_geometric_reproduction_error"]
                    ),
                    "minimum_expanding_gap": sf(sector["plus_gap_minimum"]),
                    "minimum_contracting_gap": sf(sector["minus_gap_minimum"]),
                    "maximum_expanding_direct_schur_distance": sf(
                        sector["plus_direct_schur_maximum"]
                    ),
                    "maximum_contracting_direct_schur_distance": sf(
                        sector["minus_direct_schur_maximum"]
                    ),
                    "comparisons": {
                        key: public_comparison(value)
                        for key, value in sector["comparisons"].items()
                    },
                }
                for sector in record["sectors"]
            ],
        }
        for parity, record in records.items()
    },
    "classification": {
        "exact_lapse_alignment": (
            "DERIVED COMPUTATIONAL"
            if outcome == "HYPERBOLIC_GEOMETRIC_LAPSE_IDENTIFIED"
            else "REFUTED_OR_OPEN"
        ),
        "small_angles": "PATTERN ONLY",
        "gauge_interpretation": "OPEN",
        "physical_instability": "OPEN",
        "curvature_observable": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Labels: {dict(label_counts)}")
for key, counts in hit_fractions.items():
    print(f"{key}: {dict(counts)} / 14")
for parity, record in records.items():
    distances = [
        item["distance"] for item in all_comparisons if item["parity"] == parity
    ]
    print(f"{parity}: distance range {min(distances):.6g} ... {max(distances):.6g}")
print(f"Results: {passed}/{tests} tests passed.")
sys.exit(0 if passed == tests else 1)

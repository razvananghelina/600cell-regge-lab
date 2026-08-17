#!/usr/bin/env python3
"""Identify the unique homogeneous internal-curvature kernel line.

Prior-art commit: 9177531.
Protocol commit: 53dc168.
Angular-resolution gate commit: c2cbcd3.
The near-minus target and the complete 20-comparison census are disclosed.
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
CURVATURE_INPUT = HERE / "gravity_600cell_dust_internal_curvature_response.json"
CURVATURE_SOURCE = HERE / "verify_gravity_600cell_dust_internal_curvature_response.py"
TANGENT_INPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
ALIGNMENT_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
OUTPUT = HERE / "gravity_600cell_dust_homogeneous_curvature_kernel.json"

PRIOR_ART_COMMIT = "9177531"
PROTOCOL_COMMIT = "53dc168"
RESOLUTION_CORRECTION_COMMIT = "c2cbcd3"
EXPECTED_HASHES = {
    "curvature": "95b6edd8e21ad20a0db97a7c8e7027db7da6547b2b994ad1eb595cf2307f29dc",
    "curvature_source": "276982879fae5f8fa735f27a6fa30bfe965dc3e41c169d8a229a61c23511ae66",
    "tangent": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "alignment_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
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


def norm2(matrix):
    singular = la.svdvals(np.asarray(matrix))
    return float(singular[0]) if len(singular) else 0.0


hashes = {
    "curvature": sha256(CURVATURE_INPUT),
    "curvature_source": sha256(CURVATURE_SOURCE),
    "tangent": sha256(TANGENT_INPUT),
    "tangent_numeric": sha256(TANGENT_NUMERIC),
    "tangent_source": sha256(TANGENT_SOURCE),
    "alignment_source": sha256(ALIGNMENT_SOURCE),
    "rank_source": sha256(RANK_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
    "tick": sha256(TICK_INPUT),
}
curvature_input = json.loads(CURVATURE_INPUT.read_text())
tangent_input = json.loads(TANGENT_INPUT.read_text())
tick = json.loads(TICK_INPUT.read_text())
numeric = np.load(TANGENT_NUMERIC)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and curvature_input["outcome"] == "STRONG_TANGENT_CURVATURE_INJECTIVE"
    and curvature_input["passed"] == curvature_input["tests"] == 14
    and curvature_input["full_label_counts"]
    == {"INJECTIVE": 12, "PARTIAL_OR_OPEN": 2}
    and tangent_input["outcome"] == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and tangent_input["numeric_archive_arrays"] == len(numeric.files) == 224
    and tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
)
check("all target-disclosed inputs have exact frozen provenance", provenance_ok, str(hashes))


spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_homogeneous_curvature_kernel", GEOMETRY_SOURCE
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
    "singular_record",
})
load_named_functions(ALIGNMENT_SOURCE, {
    "response_and_lift_ball",
    "geometric_lapse_matrix",
    "projected_geometric_lapse",
})

models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}


def orthonormal_columns(matrix):
    q, _ = la.qr(matrix, mode="economic")
    return q


def line_to_subspace_distance(line, subspace):
    vector = np.asarray(line, dtype=np.complex128).reshape(-1)
    vector = vector / np.linalg.norm(vector)
    basis = orthonormal_columns(np.asarray(subspace, dtype=np.complex128))
    residual = vector - basis @ (basis.conj().T @ vector)
    return float(np.linalg.norm(residual))


def line_distance(left, right):
    return line_to_subspace_distance(left, np.asarray(right).reshape(-1, 1))


def kernel_record(f_matrix):
    _, singular, vh = la.svd(f_matrix, full_matrices=False, lapack_driver="gesvd")
    line = vh.conj().T[:, -1]
    gap = float(singular[-2] - singular[-1])
    residual = float(np.linalg.norm(f_matrix @ line))
    binary_matrix_error = (
        10 * np.finfo(float).eps * max(f_matrix.shape) * max(1.0, norm2(f_matrix))
    )
    angular_error = (residual + binary_matrix_error) / max(1e-300, gap)
    return {
        "line": line / np.linalg.norm(line),
        "singular": singular,
        "gap": gap,
        "separation_ratio": float(singular[-2] / max(1e-300, singular[-1])),
        "residual": residual,
        "angular_error": angular_error,
    }


def near_minus_record(tangent, tangent_radius):
    eigenvalues, left, right = la.eig(tangent, left=True, right=True)
    distances = np.abs(eigenvalues + 1)
    order = np.argsort(distances)
    selected = order[:2]
    third = order[2]
    selected_values = eigenvalues[selected]
    contracting_index = selected[np.argmin(np.abs(selected_values))]
    expanding_index = selected[np.argmax(np.abs(selected_values))]
    contracting = right[:, contracting_index]
    expanding = right[:, expanding_index]
    contracting /= np.linalg.norm(contracting)
    expanding /= np.linalg.norm(expanding)

    boundary = math.sqrt(distances[order[1]] * distances[third])
    selector = lambda value: abs(value + 1) < boundary
    _, schur_vectors, selected_count = la.schur(
        tangent, output="complex", sort=selector
    )
    schur_plane = schur_vectors[:, :selected_count]
    direct_plane = orthonormal_columns(right[:, selected])
    plane_distance = max(
        line_to_subspace_distance(direct_plane[:, column], schur_plane)
        for column in range(direct_plane.shape[1])
    )

    conditions = {}
    errors = {}
    separations = {}
    for name, index in (
        ("contracting", contracting_index), ("expanding", expanding_index)
    ):
        overlap = abs(np.vdot(left[:, index], right[:, index]))
        condition = (
            np.linalg.norm(left[:, index]) * np.linalg.norm(right[:, index])
            / max(1e-300, overlap)
        )
        separation = float(np.min(np.abs(
            eigenvalues[index] - np.delete(eigenvalues, index)
        )))
        residual = float(np.linalg.norm(
            tangent @ right[:, index] - eigenvalues[index] * right[:, index]
        ))
        binary = 10 * np.finfo(float).eps * tangent.shape[0] * max(1.0, norm2(tangent))
        conditions[name] = float(condition)
        separations[name] = separation
        errors[name] = float(
            (residual + (binary + tangent_radius) * condition)
            / max(1e-300, separation)
        )

    plane_separation = float(
        distances[third] - distances[order[1]]
    )
    plane_error = float(
        plane_distance
        + (10 * np.finfo(float).eps * tangent.shape[0] * max(1.0, norm2(tangent))
           + tangent_radius) / max(1e-300, plane_separation)
    )
    reciprocal_defect = float(abs(
        eigenvalues[expanding_index]
        - 1 / np.conjugate(eigenvalues[contracting_index])
    ))
    return {
        "contracting_line": contracting,
        "expanding_line": expanding,
        "plane": schur_plane,
        "selected_count": int(selected_count),
        "selection_gap": float(distances[third] / distances[order[1]]),
        "contracting_eigenvalue": eigenvalues[contracting_index],
        "expanding_eigenvalue": eigenvalues[expanding_index],
        "reciprocal_defect": reciprocal_defect,
        "maximum_selected_imaginary": float(max(abs(np.imag(selected_values)))),
        "direct_schur_plane_distance": plane_distance,
        "line_errors": errors,
        "line_conditions": conditions,
        "line_separations": separations,
        "plane_error": plane_error,
    }


def comparison_record(variant_data, candidate):
    distances = {}
    source_errors = []
    target_errors = []
    binary_conditions = []
    for name, data in variant_data.items():
        if candidate == "near_minus_contracting":
            line = data["K"]
            target = data["near"]["contracting_line"].reshape(-1, 1)
            target_error = data["near"]["line_errors"]["contracting"]
        elif candidate == "near_minus_expanding":
            line = data["K"]
            target = data["near"]["expanding_line"].reshape(-1, 1)
            target_error = data["near"]["line_errors"]["expanding"]
        elif candidate == "near_minus_plane":
            line = data["K"]
            target = data["near"]["plane"]
            target_error = data["near"]["plane_error"]
        elif candidate in {
            "uniform_position_line", "uniform_momentum_line",
            "uniform_phase_plane", "pure_position", "pure_momentum",
        }:
            line = data["K"]
            target = data["exact_candidates"][candidate]
            target_error = 0.0
        elif candidate == "canonical_weak_lift":
            line = data["transported_K"]
            target = data["lift_midpoint"]
            target_error = data["lift_error"]
        elif candidate == "geometric_lapse":
            line = data["transported_K"]
            target = data["geometric"]
            target_error = 0.0
        else:
            raise ValueError(candidate)
        distances[name] = line_to_subspace_distance(line, target)
        source_errors.append(
            data["transported_error"]
            if candidate in {"canonical_weak_lift", "geometric_lapse"}
            else data["kernel_error"]
        )
        target_errors.append(target_error)
        binary_conditions.append(max(1.0, np.linalg.cond(target)))

    op = distances["operational_primary"]
    epsilon_step = (
        abs(op - distances["operational_shadow"])
        + abs(distances["validation_primary"] - distances["validation_shadow"])
        + abs(op - distances["validation_primary"])
    )
    epsilon_source = max(source_errors)
    epsilon_target = max(target_errors)
    epsilon_binary = (
        10 * np.finfo(float).eps * 160 * max(binary_conditions)
    )
    epsilon = (
        epsilon_step + epsilon_source + epsilon_target + epsilon_binary + 1e-70
    )
    if epsilon >= 1e-2:
        label = "NUMERICALLY_OPEN"
    elif op <= 10 * epsilon:
        label = "IDENTIFIED"
    elif op > 100 * epsilon:
        label = "SEPARATED"
    else:
        label = "NUMERICALLY_OPEN"
    return {
        "distance": op,
        "epsilon_distance": epsilon,
        "epsilon_step": epsilon_step,
        "epsilon_source": epsilon_source,
        "epsilon_target": epsilon_target,
        "epsilon_binary": epsilon_binary,
        "label": label,
        "variant_distances": distances,
    }


def invariance_record(variant_data):
    distances = {}
    errors = []
    multipliers = {}
    residuals = {}
    for name, data in variant_data.items():
        line = data["K"]
        image = data["tangent"] @ line
        distances[name] = line_distance(line, image)
        multiplier = np.vdot(line, image)
        multipliers[name] = multiplier
        residuals[name] = float(
            np.linalg.norm(image - multiplier * line) / max(1.0, norm2(data["tangent"]))
        )
        errors.append(
            data["kernel_error"] * (1 + norm2(data["tangent"]))
            / max(1e-300, np.linalg.norm(image))
            + 10 * np.finfo(float).eps * 60
        )
    op = distances["operational_primary"]
    epsilon_step = (
        abs(op - distances["operational_shadow"])
        + abs(distances["validation_primary"] - distances["validation_shadow"])
        + abs(op - distances["validation_primary"])
    )
    epsilon = epsilon_step + max(errors) + 1e-70
    if epsilon >= 1e-2:
        label = "NUMERICALLY_OPEN"
    elif op <= 10 * epsilon:
        label = "IDENTIFIED"
    elif op > 100 * epsilon:
        label = "SEPARATED"
    else:
        label = "NUMERICALLY_OPEN"
    return {
        "distance": op,
        "epsilon_distance": epsilon,
        "epsilon_step": epsilon_step,
        "label": label,
        "multiplier": multipliers["operational_primary"],
        "relative_residual": residuals["operational_primary"],
        "variant_distances": distances,
    }


def cross_schedule_record(left, right):
    distances = {
        name: line_distance(left[name]["K"], right[name]["K"])
        for name in VARIANTS
    }
    op = distances["operational_primary"]
    epsilon_step = (
        abs(op - distances["operational_shadow"])
        + abs(distances["validation_primary"] - distances["validation_shadow"])
        + abs(op - distances["validation_primary"])
    )
    epsilon = (
        epsilon_step
        + max(data["kernel_error"] for data in left.values())
        + max(data["kernel_error"] for data in right.values())
        + 10 * np.finfo(float).eps * 60 + 1e-70
    )
    if epsilon >= 1e-2:
        label = "NUMERICALLY_OPEN"
    elif op <= 10 * epsilon:
        label = "IDENTIFIED"
    elif op > 100 * epsilon:
        label = "SEPARATED"
    else:
        label = "NUMERICALLY_OPEN"
    return {
        "distance": op,
        "epsilon_distance": epsilon,
        "epsilon_step": epsilon_step,
        "label": label,
        "variant_distances": distances,
    }


def sf(value):
    return f"{float(value):.17e}"


def sc(value):
    return {"real": sf(np.real(value)), "imaginary": sf(np.imag(value))}


def public_comparison(record):
    return {
        "distance": sf(record["distance"]),
        "epsilon_distance": sf(record["epsilon_distance"]),
        "epsilon_step": sf(record["epsilon_step"]),
        "label": record["label"],
        "variant_distances": {
            name: sf(value) for name, value in record["variant_distances"].items()
        },
        **({
            "epsilon_source": sf(record["epsilon_source"]),
            "epsilon_target": sf(record["epsilon_target"]),
            "epsilon_binary": sf(record["epsilon_binary"]),
        } if "epsilon_source" in record else {}),
    }


print("=" * 78)
print("IDENTIFICATION OF THE HOMOGENEOUS INTERNAL-CURVATURE KERNEL LINE")
print("=" * 78)

candidate_names = (
    "near_minus_contracting",
    "near_minus_expanding",
    "near_minus_plane",
    "uniform_position_line",
    "uniform_momentum_line",
    "uniform_phase_plane",
    "pure_position",
    "pure_momentum",
    "canonical_weak_lift",
    "geometric_lapse",
)

records = {}
variant_records_by_parity = {}
global_controls = provenance_ok and gro.tests == gro.passed == 43
old_boundary_orderings = {}

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing the trivial 60-dimensional response", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    index_data = group_and_index_data(model, state)
    geometry = prepare_geometry(model, index_data)
    sectors, sector_control = high_precision_sector_bases(index_data)
    trivial_indices = [
        index for index, sector in enumerate(sectors)
        if sector["constant_overlap"] > mp.mpf("0.5")
    ]
    trivial_ok = len(trivial_indices) == 1
    sector_index = trivial_indices[0]
    sector = sectors[sector_index]
    old_boundary_orderings[parity] = index_data["orbit_edges"][:30]

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
    curvature_data = triangle_response_data(
        model, index_data, geometry, pattern_cache
    )
    hessian_kernels, hessian_control = assemble_full_representative_kernels(
        index_data, geometry, pattern_cache
    )
    weak_positions = [
        position for position in range(35)
        if index_data["edge_kind"][24 * (30 + position)] == "pole"
    ]
    geometric_full, geometry_coefficient = geometric_lapse_matrix(
        index_data, weak_positions, state
    )
    geometric = projected_geometric_lapse(geometric_full, sector)

    reconstruction_control = bool(
        trivial_ok and sector["dimension"] == 1
        and len(curvature_data["internal_global_types"]) == 160
        and branch_control["entry_pass"]
        and all(
            value < ARITHMETIC_FLOOR
            for value in curvature_data["maximum_derivative_imaginary"].values()
        )
        and all(
            value < ARITHMETIC_FLOOR
            for value in curvature_data["maximum_equivariance_residual"].values()
        )
        and hessian_control["maximum_imaginary"] < ARITHMETIC_FLOOR
        and len(weak_positions) == 5
        and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
    )
    check(
        f"{parity}: trivial sector, curvature kernel and five pole types reconstruct",
        reconstruction_control,
        f"sector={sector_index}, weak={weak_positions}",
    )

    hessian_blocks = {
        name: project_full_kernel(kernel, sector)
        for name, kernel in hessian_kernels.items()
    }
    curvature_blocks = {
        name: project_curvature_kernel(kernel, sector)
        for name, kernel in curvature_data["kernels"].items()
    }

    variant_data = {}
    determinants_ok = True
    near_controls_ok = True
    for name in VARIANTS:
        response = response_and_lift_ball(
            hessian_blocks[name], 1, weak_positions
        )
        determinants_ok &= not response["det_j"].contains(0)
        old_projection = np.zeros((30, 60), dtype=np.complex128)
        old_projection[:, :30] = np.eye(30)
        z_matrix = np.vstack((old_projection, response["response_midpoint"]))
        z_radii = np.vstack((
            np.zeros_like(old_projection, dtype=float), response["response_radii"]
        ))
        d_matrix = mp_to_numpy(curvature_blocks[name])
        f_matrix = d_matrix @ z_matrix
        kernel = kernel_record(f_matrix)
        tangent = numeric[f"{parity}_sector{sector_index}_{name}_tangent_midpoint"]
        tangent_radii = numeric[f"{parity}_sector{sector_index}_{name}_tangent_radii"]
        tangent_radius = float(la.norm(tangent_radii, "fro"))
        near = near_minus_record(tangent, tangent_radius)
        near_controls_ok &= bool(
            near["selected_count"] == 2
            and near["selection_gap"] > 10
            and abs(near["contracting_eigenvalue"]) < 1
            and abs(near["expanding_eigenvalue"]) > 1
            and near["maximum_selected_imaginary"] < 1e-8
            and near["reciprocal_defect"] < 1e-8
            and near["direct_schur_plane_distance"] < 1e-6
        )

        k_line = kernel["line"]
        transported = response["response_midpoint"] @ k_line
        transported_norm = max(1e-300, float(np.linalg.norm(transported)))
        y_norm = norm2(response["response_midpoint"])
        response_radius = float(la.norm(response["response_radii"], "fro"))
        transported_error = (
            y_norm * kernel["angular_error"] + response_radius
        ) / transported_norm
        lift_singular = la.svdvals(response["lift_midpoint"])
        lift_error = float(
            la.norm(response["lift_radii"], "fro")
            / max(1e-300, lift_singular[-1])
        )

        uniform_q = np.r_[np.ones(30), np.zeros(30)].reshape(-1, 1)
        uniform_p = np.r_[np.zeros(30), np.ones(30)].reshape(-1, 1)
        exact_candidates = {
            "uniform_position_line": uniform_q,
            "uniform_momentum_line": uniform_p,
            "uniform_phase_plane": np.hstack((uniform_q, uniform_p)),
            "pure_position": np.vstack((np.eye(30), np.zeros((30, 30)))),
            "pure_momentum": np.vstack((np.zeros((30, 30)), np.eye(30))),
        }
        variant_data[name] = {
            "F": f_matrix,
            "D_norm": norm2(d_matrix),
            "Z_radius": float(la.norm(z_radii, "fro")),
            "K": k_line,
            "kernel": kernel,
            "kernel_error": kernel["angular_error"],
            "tangent": tangent,
            "near": near,
            "exact_candidates": exact_candidates,
            "transported_K": transported,
            "transported_error": transported_error,
            "lift_midpoint": response["lift_midpoint"],
            "lift_error": lift_error,
            "geometric": geometric,
        }

    full_rank = singular_record(variant_data)
    rank_ok = bool(
        full_rank["resolved_rank"] == 59
        and full_rank["zero_count"] == 1
        and full_rank["open_count"] == 0
        and full_rank["columns"] == 60
        and all(
            data["kernel"]["separation_ratio"] > 1e4
            for data in variant_data.values()
        )
    )
    check(
        f"{parity}: the full response independently reproduces rank 59 and one line",
        determinants_ok and rank_ok,
        f"rank={full_rank['resolved_rank']}, zero={full_rank['zero_count']}",
    )
    check(
        f"{parity}: the disclosed near-minus pair and independent Schur plane pass",
        near_controls_ok,
        "gap=" + sf(min(
            data["near"]["selection_gap"] for data in variant_data.values()
        )),
    )

    comparisons = {
        candidate: comparison_record(variant_data, candidate)
        for candidate in candidate_names
    }
    invariance = invariance_record(variant_data)
    complete = len(comparisons) == 10 and all(
        record["label"] in {"IDENTIFIED", "SEPARATED", "NUMERICALLY_OPEN"}
        for record in comparisons.values()
    )
    check(
        f"{parity}: all ten candidate comparisons and tangent invariance are labelled",
        complete and invariance["label"] in {
            "IDENTIFIED", "SEPARATED", "NUMERICALLY_OPEN"
        },
        str(Counter(item["label"] for item in comparisons.values())),
    )

    controls_ok = bool(
        reconstruction_control and determinants_ok and rank_ok
        and near_controls_ok and complete
    )
    global_controls &= controls_ok
    records[parity] = {
        "controls_ok": controls_ok,
        "sector_index": sector_index,
        "weak_positions": weak_positions,
        "geometry_coefficient": geometry_coefficient,
        "full_rank": full_rank,
        "comparisons": comparisons,
        "invariance": invariance,
        "variant_data": variant_data,
    }
    variant_records_by_parity[parity] = variant_data


literal_ordering_ok = old_boundary_orderings["even"] == old_boundary_orderings["odd"]
cross_schedule = cross_schedule_record(
    variant_records_by_parity["even"], variant_records_by_parity["odd"]
)
check(
    "the two schedule kernels are compared in a literal common boundary ordering",
    literal_ordering_ok,
)
check(
    "the cross-schedule kernel comparison receives its calibrated label",
    cross_schedule["label"] in {"IDENTIFIED", "SEPARATED", "NUMERICALLY_OPEN"},
    f"label={cross_schedule['label']}, distance={cross_schedule['distance']:.3e}",
)


def both_identified(candidate):
    return all(
        records[parity]["comparisons"][candidate]["label"] == "IDENTIFIED"
        for parity in ("even", "odd")
    )


identified_eigenline = next((
    candidate for candidate in (
        "near_minus_contracting", "near_minus_expanding"
    )
    if both_identified(candidate)
    and all(records[p]["invariance"]["label"] == "IDENTIFIED"
            for p in ("even", "odd"))
), None)
identified_uniform_line = next((
    candidate for candidate in (
        "uniform_position_line", "uniform_momentum_line"
    ) if both_identified(candidate)
), None)
localized_subspaces = [
    candidate for candidate in candidate_names
    if candidate not in {
        "near_minus_contracting", "near_minus_expanding",
        "uniform_position_line", "uniform_momentum_line",
    } and both_identified(candidate)
]

if not global_controls:
    outcome = "HOMOGENEOUS_CURVATURE_KERNEL_CONTROL_FAILED"
elif cross_schedule["label"] == "SEPARATED":
    outcome = "HOMOGENEOUS_CURVATURE_KERNEL_SCHEDULE_DEPENDENT"
elif identified_eigenline is not None:
    outcome = "HOMOGENEOUS_CURVATURE_KERNEL_EIGENLINE_IDENTIFIED"
elif identified_uniform_line is not None:
    outcome = "HOMOGENEOUS_CURVATURE_KERNEL_LINE_IDENTIFIED"
elif localized_subspaces:
    outcome = "HOMOGENEOUS_CURVATURE_KERNEL_SUBSPACE_LOCALIZED"
else:
    outcome = "HOMOGENEOUS_CURVATURE_KERNEL_UNIDENTIFIED_OR_OPEN"

all_comparisons = [
    record for parity in ("even", "odd")
    for record in records[parity]["comparisons"].values()
]
label_counts = Counter(record["label"] for record in all_comparisons)
check(
    "the preregistered look-elsewhere ledger contains exactly 20 comparisons",
    len(all_comparisons) == 20 and sum(label_counts.values()) == 20,
    str(dict(label_counts)),
)
check(
    "the frozen hierarchy assigns the homogeneous-kernel outcome",
    outcome in {
        "HOMOGENEOUS_CURVATURE_KERNEL_CONTROL_FAILED",
        "HOMOGENEOUS_CURVATURE_KERNEL_SCHEDULE_DEPENDENT",
        "HOMOGENEOUS_CURVATURE_KERNEL_EIGENLINE_IDENTIFIED",
        "HOMOGENEOUS_CURVATURE_KERNEL_LINE_IDENTIFIED",
        "HOMOGENEOUS_CURVATURE_KERNEL_SUBSPACE_LOCALIZED",
        "HOMOGENEOUS_CURVATURE_KERNEL_UNIDENTIFIED_OR_OPEN",
    },
    f"outcome={outcome}",
)


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "resolution_correction_commit": RESOLUTION_CORRECTION_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "candidate_count_per_schedule": 10,
    "total_candidate_comparisons": 20,
    "outcome": outcome,
    "identified_eigenline": identified_eigenline,
    "identified_uniform_line": identified_uniform_line,
    "consistently_localized_subspaces": localized_subspaces,
    "label_counts": dict(label_counts),
    "cross_schedule": public_comparison(cross_schedule),
    "parities": {
        parity: {
            "controls_ok": record["controls_ok"],
            "trivial_sector_index": record["sector_index"],
            "weak_orbit_positions": record["weak_positions"],
            "geometry_coefficient": sf(record["geometry_coefficient"]),
            "full_response_rank": record["full_rank"]["resolved_rank"],
            "full_response_zero_count": record["full_rank"]["zero_count"],
            "full_response_open_count": record["full_rank"]["open_count"],
            "minimum_kernel_separation_ratio": sf(min(
                data["kernel"]["separation_ratio"]
                for data in record["variant_data"].values()
            )),
            "near_minus": {
                "contracting_eigenvalue": sc(
                    record["variant_data"]["operational_primary"]["near"][
                        "contracting_eigenvalue"
                    ]
                ),
                "expanding_eigenvalue": sc(
                    record["variant_data"]["operational_primary"]["near"][
                        "expanding_eigenvalue"
                    ]
                ),
                "selection_gap": sf(
                    record["variant_data"]["operational_primary"]["near"][
                        "selection_gap"
                    ]
                ),
                "reciprocal_defect": sf(
                    record["variant_data"]["operational_primary"]["near"][
                        "reciprocal_defect"
                    ]
                ),
                "direct_schur_plane_distance": sf(
                    record["variant_data"]["operational_primary"]["near"][
                        "direct_schur_plane_distance"
                    ]
                ),
            },
            "tangent_invariance": {
                **public_comparison(record["invariance"]),
                "multiplier": sc(record["invariance"]["multiplier"]),
                "relative_residual": sf(record["invariance"]["relative_residual"]),
            },
            "comparisons": {
                candidate: public_comparison(comparison)
                for candidate, comparison in record["comparisons"].items()
            },
        }
        for parity, record in records.items()
    },
    "classification": {
        "kernel_line_and_memberships": "DERIVED COMPUTATIONAL",
        "near_minus_comparison": "CONFIRMATORY TARGET-DISCLOSED",
        "time_or_gauge_interpretation": "OPEN",
        "nonlinear_integrability": "OPEN",
        "physical_mode_interpretation": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Candidate labels: {dict(label_counts)} / 20")
print(
    "Cross-schedule: "
    f"{cross_schedule['label']} at {cross_schedule['distance']:.6g} "
    f"+/- {cross_schedule['epsilon_distance']:.3g}"
)
for parity in ("even", "odd"):
    hits = [
        name for name, item in records[parity]["comparisons"].items()
        if item["label"] == "IDENTIFIED"
    ]
    print(
        f"{parity}: hits={hits}, T-invariance={records[parity]['invariance']['label']}"
    )
print(f"Results: {passed}/{tests} tests passed.")
sys.exit(0 if passed == tests else 1)

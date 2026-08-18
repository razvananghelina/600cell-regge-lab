#!/usr/bin/env python3
"""Certified unit-circle root count for the negative-shape recurrence.

Prior-art commit: 3467728.
Target-disclosed protocol commit: 925811e.
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

STIFFNESS_JSON = HERE / "gravity_600cell_dust_shape_stiffness.json"
STIFFNESS_SOURCE = HERE / "verify_gravity_600cell_dust_shape_stiffness.py"
CENTERED_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
LEGENDRE_JSON = HERE / "gravity_600cell_dust_full_anisotropic_legendre_rank.json"
DYNAMICS_JSON = HERE / "gravity_600cell_dust_negative_shape_dynamics.json"
DYNAMICS_SOURCE = HERE / "verify_gravity_600cell_dust_negative_shape_dynamics.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
OUTPUT = HERE / "gravity_600cell_dust_negative_shape_root_count.json"

PRIOR_ART_COMMIT = "3467728"
PROTOCOL_COMMIT = "925811e"
EXPECTED_HASHES = {
    "stiffness_json": "03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868",
    "stiffness_source": "d4f0a9a805910de37011ba70f407907daa2d11c650aeea22e571ab867282a44c",
    "centered_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "legendre_json": "7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226",
    "dynamics_json": "51ff46a529958491d7338c62d34612721352c502e926bb4c1df98aa477c9a854",
    "dynamics_source": "6e7659ca398037e806f9a35a9f3db3d6035f992a8655699b47a2519b0c37453e",
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
INITIAL_INTERVALS = 256
MAX_DEPTH = 32
MAX_INTERVALS = 2_000_000
WINDING_MESHES = (8192, 16384)
SAFETIES = (1, 100)
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


def load_audited_helpers():
    wanted = {
        "mp_frobenius", "mp_submatrix", "cluster_sorted", "orbit_sort_key",
        "edge_image", "group_data", "incidence_data", "mp_to_numpy",
        "component_reenclosure_radii",
    }
    tree = ast.parse(CONFORMAL_SOURCE.read_text(), filename=str(CONFORMAL_SOURCE))
    body = [node for node in tree.body if isinstance(node, ast.FunctionDef)
            and node.name in wanted]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited helper mismatch: missing={wanted-found}")
    exec(compile(ast.Module(body=body, type_ignores=[]),
                 str(CONFORMAL_SOURCE), "exec"), globals())

    tree = ast.parse(FULL_SOURCE.read_text(), filename=str(FULL_SOURCE))
    body = [node for node in tree.body if isinstance(node, ast.FunctionDef)
            and node.name == "high_precision_sector_bases"]
    if len(body) != 1:
        raise RuntimeError("audited high-precision sector function is missing")
    exec(compile(ast.Module(body=body, type_ignores=[]),
                 str(FULL_SOURCE), "exec"), globals())


def source_matrix(archive, prefix, name):
    midpoint = np.asarray(archive[f"{prefix}_{name}_midpoint"])
    stored = np.asarray(archive[f"{prefix}_{name}_radii"])
    return midpoint, component_reenclosure_radii(midpoint, stored)


def matrix_error(midpoint, radii, n):
    return float(la.norm(radii, "fro") + 1000 * MACHINE_EPSILON * n
                 * max(1.0, operator_norm(midpoint)))


def restriction_error(midpoint, base_error, eta_subspace, n):
    norm = operator_norm(midpoint)
    return float(base_error + 2 * eta_subspace * (norm + base_error)
                 + 1000 * MACHINE_EPSILON * n * max(1.0, norm))


def quadratic_value(a2, a1, a0, z):
    return a2 * (z * z) + a1 * z + a0


def contour_cover(a2, a1, a0, epsilon_g, epsilon_o, safety):
    """Cover the unit circle with the preregistered Lipschitz lower bound."""
    m = a0.shape[0]
    norm_a2, norm_a1, norm_a0 = map(operator_norm, (a2, a1, a0))
    evaluation_floor = float(1000 * MACHINE_EPSILON * m
                             * max(1.0, 2 * norm_a2 + norm_a1 + norm_a0))
    lipschitz_q = float(2 * norm_a2 + norm_a1)
    lipschitz_delta = float(2 * epsilon_g)
    width = 2 * math.pi / INITIAL_INTERVALS
    stack = [(j * width, (j + 1) * width, 0)
             for j in reversed(range(INITIAL_INTERVALS))]
    evaluated = certified = maximum_depth = 0
    weakest_lower = weakest_point_ratio = math.inf
    failure = None

    while stack:
        if evaluated >= MAX_INTERVALS:
            failure = {"kind": "RESOURCE_LIMIT", "evaluated": evaluated}
            break
        left, right, depth = stack.pop()
        theta = (left + right) / 2
        half_width = (right - left) / 2
        z = complex(math.cos(theta), math.sin(theta))
        minimum = float(la.svdvals(quadratic_value(a2, a1, a0, z))[-1])
        delta = float(epsilon_g * abs(z * z - 1) + epsilon_o)
        uncertainty = delta + evaluation_floor
        point_margin = float(minimum - safety * uncertainty)
        ratio = float(minimum / uncertainty) if uncertainty else math.inf
        lower = float(point_margin
                      - (lipschitz_q + safety * lipschitz_delta) * half_width)
        evaluated += 1
        maximum_depth = max(maximum_depth, depth)
        weakest_lower = min(weakest_lower, lower)
        weakest_point_ratio = min(weakest_point_ratio, ratio)
        diagnostic = {
            "theta": sf(theta), "sigma_min": sf(minimum),
            "delta": sf(delta), "evaluation_floor": sf(evaluation_floor),
            "point_margin": sf(point_margin), "ratio": sf(ratio),
            "depth": depth,
        }
        if point_margin <= 0:
            failure = {"kind": "POINTWISE_MARGIN_NONPOSITIVE", **diagnostic}
            break
        if lower > 0:
            certified += 1
            continue
        if depth >= MAX_DEPTH:
            failure = {"kind": "MAX_DEPTH", **diagnostic,
                       "interval_lower": sf(lower)}
            break
        middle = (left + right) / 2
        stack.append((middle, right, depth + 1))
        stack.append((left, middle, depth + 1))

    return {
        "passed": failure is None, "safety": safety,
        "evaluated_intervals": evaluated,
        "certified_leaf_intervals": certified,
        "maximum_depth": maximum_depth,
        "weakest_interval_lower": sf(weakest_lower),
        "weakest_sample_ratio": sf(weakest_point_ratio),
        "evaluation_floor": sf(evaluation_floor),
        "lipschitz_Q": sf(lipschitz_q),
        "lipschitz_delta": sf(lipschitz_delta),
        "failure": failure,
    }


def determinant_winding(a2, a1, a0, intervals, orientation=1):
    phases = []
    singular_samples = 0
    for j in range(intervals + 1):
        theta = orientation * 2 * math.pi * j / intervals
        z = complex(math.cos(theta), math.sin(theta))
        sign, _ = np.linalg.slogdet(quadratic_value(a2, a1, a0, z))
        if abs(sign) == 0:
            singular_samples += 1
            phases.append(float("nan"))
        else:
            phases.append(float(np.angle(sign)))
    phases = np.asarray(phases)
    if singular_samples:
        return {"intervals": intervals, "orientation": orientation,
                "singular_samples": singular_samples, "raw_winding": None,
                "corrected_winding": None,
                "maximum_principal_phase_increment": None}
    increments = np.angle(np.exp(1j * np.diff(phases)))
    raw = int(round(float(np.sum(increments) / (2 * math.pi))))
    return {
        "intervals": intervals, "orientation": orientation,
        "singular_samples": 0, "raw_winding": raw,
        "corrected_winding": raw if orientation == 1 else -raw,
        "maximum_principal_phase_increment": sf(np.max(np.abs(increments))),
    }


def midpoint_root_count(a2, a1, a0):
    m = a0.shape[0]
    identity = np.eye(m, dtype=np.complex128)
    companion = np.block([
        [np.zeros((m, m), dtype=np.complex128), identity],
        [la.solve(a2, -a0), la.solve(a2, -a1)],
    ])
    values, vectors = la.eig(companion)
    moduli = np.abs(values)
    inside, outside = int(np.sum(moduli < 1)), int(np.sum(moduli > 1))
    on = int(len(values) - inside - outside)
    norms = tuple(map(operator_norm, (a2, a1, a0)))
    maximum_residual = 0.0
    for index, value in enumerate(values):
        x = vectors[:m, index]
        xnorm = la.norm(x)
        denominator = (norms[0] * abs(value) ** 2
                       + norms[1] * abs(value) + norms[2]) * xnorm
        residual = la.norm(quadratic_value(a2, a1, a0, value) @ x)
        maximum_residual = max(
            maximum_residual,
            float(residual / denominator) if denominator else math.inf,
        )
    windings = [determinant_winding(a2, a1, a0, n)
                for n in WINDING_MESHES]
    fine = windings[-1]
    winding_ok = bool(
        all(item["singular_samples"] == 0 for item in windings)
        and len({item["corrected_winding"] for item in windings}) == 1
        and fine["corrected_winding"] == inside
        and float(fine["maximum_principal_phase_increment"]) < math.pi / 4
    )
    return {
        "companion_counts": {"inside": inside, "on": on, "outside": outside},
        "minimum_modulus_distance_from_one": sf(np.min(np.abs(moduli - 1))),
        "maximum_normalized_polynomial_residual": sf(maximum_residual),
        "residual_gate": maximum_residual < 1e-10,
        "windings": windings, "winding_gate": winding_ok,
    }


def synthetic_controls():
    identity = np.eye(2, dtype=np.complex128)
    zero = np.zeros((2, 2), dtype=np.complex128)
    free = contour_cover(identity, -2 * identity, identity, 0.0, 0.0, 1)
    unit = contour_cover(identity, zero, identity, 0.0, 0.0, 1)
    hyper = contour_cover(identity, -2.5 * identity, identity,
                          1e-12, 1e-12, 100)
    hyper_count = midpoint_root_count(identity, -2.5 * identity, identity)
    forward = determinant_winding(identity, -2.5 * identity, identity, 8192, 1)
    reverse = determinant_winding(identity, -2.5 * identity, identity, 8192, -1)
    ok = bool(
        not free["passed"] and not unit["passed"] and hyper["passed"]
        and hyper_count["companion_counts"]
        == {"inside": 2, "on": 0, "outside": 2}
        and hyper_count["residual_gate"] and hyper_count["winding_gate"]
        and forward["corrected_winding"] == reverse["corrected_winding"] == 2
    )
    return {
        "passed": ok, "free_literal_cover": free,
        "unit_literal_cover": unit, "hyperbolic_100x_cover": hyper,
        "hyperbolic_count": hyper_count,
        "hyperbolic_forward_winding": forward,
        "hyperbolic_reverse_winding": reverse,
    }


print("=" * 78)
print("600-CELL NEGATIVE SHAPE COEFFICIENT-BALL ROOT COUNT")
print("=" * 78)

paths = {
    "stiffness_json": STIFFNESS_JSON,
    "stiffness_source": STIFFNESS_SOURCE,
    "centered_json": CENTERED_JSON,
    "centered_npz": CENTERED_NPZ,
    "legendre_json": LEGENDRE_JSON,
    "dynamics_json": DYNAMICS_JSON,
    "dynamics_source": DYNAMICS_SOURCE,
    "conformal_source": CONFORMAL_SOURCE,
    "full_source": FULL_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "commons": COMMONS_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
stiffness_source = json.loads(STIFFNESS_JSON.read_text())
centered = json.loads(CENTERED_JSON.read_text())
legendre = json.loads(LEGENDRE_JSON.read_text())
dynamics = json.loads(DYNAMICS_JSON.read_text())
archive = np.load(CENTERED_NPZ, allow_pickle=False)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and stiffness_source["outcome"] == "SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED"
    and stiffness_source["passed"] == stiffness_source["tests"] == 12
    and centered["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and centered["passed"] == centered["tests"] == 7
    and len(archive.files) == centered["numeric_archive_arrays"] == 560
    and legendre["outcome"] == "FULL_CANONICAL_LEGENDRE_REGULAR"
    and legendre["passed"] == legendre["tests"] == 18
    and dynamics["outcome"] == "NEGATIVE_SHAPE_AUTONOMOUS_UNIT_CONSISTENT"
    and dynamics["passed"] == dynamics["tests"] == 12
)
check("all target-disclosed inputs have exact frozen provenance", provenance_ok)

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_negative_shape_root_count", GEOMETRY_SOURCE
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
)

controls = synthetic_controls()
check("the three polynomial controls and orientation reversal pass", controls["passed"])

records = {parity: [] for parity in PARITIES}
carrier_open = conformal_open
carrier_contradicted = conformal_contradicted
invariance_counts = Counter()
regularity_counts = Counter()
literal_counts = Counter()
safety_counts = Counter()
midpoint_failures = 0
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
            inertia = Counter(sign_label(float(value), epsilon_a) for value in values_a)
            if inertia != Counter({"NEGATIVE_RESOLVED": 15, "POSITIVE_RESOLVED": 10}):
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
                invariance[name] = {"residual": residual, "error": error,
                                    "label": label}

            both_invariant = all(
                item["label"] == "INVARIANT_CONSISTENT"
                for item in invariance.values()
            )
            if not both_invariant:
                carrier_open = True

            gamma_negative = negative_basis.conj().T @ gamma_s @ negative_basis
            omega_negative = negative_basis.conj().T @ omega_s @ negative_basis
            epsilon_gn = (invariance["Gamma"]["error"]
                          + invariance["Gamma"]["residual"])
            epsilon_on = (invariance["Omega"]["error"]
                          + invariance["Omega"]["residual"])
            identity = np.eye(15, dtype=np.complex128)
            a2 = identity + gamma_negative
            a1 = -2 * identity + omega_negative
            a0 = identity - gamma_negative

            minimum_a2 = float(la.svdvals(a2)[-1])
            leading_label = regularity_label(minimum_a2, epsilon_gn)
            regularity_counts[leading_label] += 1
            midpoint = midpoint_root_count(a2, a1, a0)
            midpoint_ok = bool(midpoint["residual_gate"] and midpoint["winding_gate"])
            midpoint_failures += int(not midpoint_ok)
            literal = contour_cover(a2, a1, a0, epsilon_gn, epsilon_on, 1)
            safety = contour_cover(a2, a1, a0, epsilon_gn, epsilon_on, 100)
            literal_counts["PASS" if literal["passed"] else "OPEN"] += 1
            safety_counts["PASS" if safety["passed"] else "OPEN"] += 1

            transferred = None
            if literal["passed"] and midpoint_ok and leading_label == "REGULAR_RESOLVED":
                transferred = dict(midpoint["companion_counts"])
                transferred["on"] = 0
            all_finite &= bool(math.isfinite(minimum_a2)
                               and math.isfinite(epsilon_gn)
                               and math.isfinite(epsilon_on))
            variants.append({
                "variant": variant, "negative_dimension": 15,
                "inherited_stiffness_inertia": dict(inertia),
                "negative_positive_gap": sf(eigengap),
                "negative_subspace_error_eta": sf(eta_e),
                "invariance": {
                    name: {"residual": sf(item["residual"]),
                           "error": sf(item["error"]), "label": item["label"]}
                    for name, item in invariance.items()
                },
                "coefficients": {
                    "norm_A2": sf(operator_norm(a2)),
                    "norm_A1": sf(operator_norm(a1)),
                    "norm_A0": sf(operator_norm(a0)),
                    "epsilon_Gamma": sf(epsilon_gn),
                    "epsilon_Omega": sf(epsilon_on),
                    "minimum_A2_singular": sf(minimum_a2),
                    "A2_regularity_label": leading_label,
                },
                "midpoint": midpoint,
                "literal_rouche_cover": literal,
                "safety_100_cover": safety,
                "transferred_counts": transferred,
            })

        records[parity].append({
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "variants": variants,
        })

required_cells = 16
carrier_ok = bool(not carrier_open and not carrier_contradicted and all_finite)
check(
    "all 16 inherited negative carriers are reconstructed without contradiction",
    carrier_ok,
    f"open={carrier_open}, contradicted={carrier_contradicted}",
)
check(
    "all 32 invariance tests remain invariant-consistent",
    invariance_counts == Counter({"INVARIANT_CONSISTENT": 32}),
    str(dict(invariance_counts)),
)
check(
    "all 16 leading coefficients remain regular-resolved",
    regularity_counts == Counter({"REGULAR_RESOLVED": required_cells}),
    str(dict(regularity_counts)),
)
check(
    "all midpoint companion and winding counts satisfy their preregistered gates",
    midpoint_failures == 0,
    f"failures={midpoint_failures}",
)
check(
    "all literal contour covers are completely classified",
    sum(literal_counts.values()) == required_cells,
    str(dict(literal_counts)),
)
check(
    "all 100x contour covers are completely classified",
    sum(safety_counts.values()) == required_cells,
    str(dict(safety_counts)),
)

all_records = [
    variant
    for parity_rows in records.values()
    for row in parity_rows
    for variant in row["variants"]
]
transferred = [item["transferred_counts"] for item in all_records]
schedule_robust = bool(
    all(item is not None for item in transferred)
    and len({(item["inside"], item["on"], item["outside"])
             for item in transferred}) == 1
)
controls_ok = bool(
    provenance_ok and geometry_import_ok and carrier_geometry_ok
    and sector_basis_ok and controls["passed"] and carrier_ok
    and invariance_counts == Counter({"INVARIANT_CONSISTENT": 32})
    and regularity_counts == Counter({"REGULAR_RESOLVED": required_cells})
)
if not controls_ok:
    outcome = "NEGATIVE_SHAPE_ROOT_COUNT_CONTROL_FAILED"
elif midpoint_failures:
    outcome = "NEGATIVE_SHAPE_ROOT_COUNT_MIDPOINT_INCONSISTENT"
elif literal_counts["OPEN"]:
    outcome = "NEGATIVE_SHAPE_ROOT_COUNT_COEFFICIENT_OPEN"
elif safety_counts["OPEN"]:
    outcome = "NEGATIVE_SHAPE_ROOT_COUNT_SAFETY_OPEN"
elif not schedule_robust:
    outcome = "NEGATIVE_SHAPE_ROOT_COUNT_SCHEDULE_DEPENDENT"
else:
    common = transferred[0]
    if common["inside"] > 0 and common["outside"] > 0:
        outcome = "NEGATIVE_SHAPE_LOCAL_HYPERBOLIC_RESOLVED"
    elif common["outside"] == 0:
        outcome = "NEGATIVE_SHAPE_LOCAL_NONEXPANDING_RESOLVED"
    else:
        outcome = "NEGATIVE_SHAPE_ROOT_COUNT_OTHER_RESOLVED"

allowed_outcomes = {
    "NEGATIVE_SHAPE_ROOT_COUNT_CONTROL_FAILED",
    "NEGATIVE_SHAPE_ROOT_COUNT_MIDPOINT_INCONSISTENT",
    "NEGATIVE_SHAPE_ROOT_COUNT_COEFFICIENT_OPEN",
    "NEGATIVE_SHAPE_ROOT_COUNT_SAFETY_OPEN",
    "NEGATIVE_SHAPE_ROOT_COUNT_SCHEDULE_DEPENDENT",
    "NEGATIVE_SHAPE_LOCAL_HYPERBOLIC_RESOLVED",
    "NEGATIVE_SHAPE_LOCAL_NONEXPANDING_RESOLVED",
    "NEGATIVE_SHAPE_ROOT_COUNT_OTHER_RESOLVED",
}
check("the preregistered outcome tree is exhausted exactly once",
      outcome in allowed_outcomes, outcome)

payload = {
    "title": "Coefficient-ball root count for the negative-shape recurrence",
    "date": "2026-08-18",
    "classification": "DERIVED COMPUTATIONAL, TARGET-DISCLOSED",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_hashes": hashes,
    "protocol": {
        "initial_intervals": INITIAL_INTERVALS,
        "maximum_depth": MAX_DEPTH,
        "maximum_intervals": MAX_INTERVALS,
        "winding_meshes": list(WINDING_MESHES),
        "safeties": list(SAFETIES),
        "contour": "unit circle, counterclockwise",
        "perturbation": "DeltaGamma*(z^2-1)+DeltaOmega*z",
        "interpretation_scope": "local frozen recurrence only",
    },
    "controls": controls,
    "counts": {
        "inherited_cells": required_cells,
        "literal_cover": dict(literal_counts),
        "safety_100_cover": dict(safety_counts),
        "midpoint_failures": midpoint_failures,
        "schedule_robust_transferred_count": schedule_robust,
    },
    "parities": records,
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print()
print("literal cover counts:", dict(literal_counts))
print("100x cover counts:", dict(safety_counts))
print("midpoint failures:", midpoint_failures)
print("outcome:", outcome)
print(f"RESULT: {passed}/{tests} tests passed")

if passed != tests:
    sys.exit(1)

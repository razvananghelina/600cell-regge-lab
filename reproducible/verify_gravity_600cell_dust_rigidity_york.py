#!/usr/bin/env python3
"""Dynamic closure test for the canonical 600-cell self-stress carrier.

Prior-art commit: a318d6e.
Protocol commit: f9b692e.
The 470/250 rigidity counts are theorem controls, not discovery targets.
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

from commons.cell600 import build_600cell  # noqa: E402


CENTERED_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
CONFORMAL_JSON = HERE / "gravity_600cell_dust_conformal_supermetric.json"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
OUTPUT = HERE / "gravity_600cell_dust_rigidity_york.json"

PRIOR_ART_COMMIT = "a318d6e"
PROTOCOL_COMMIT = "f9b692e"
EXPECTED_HASHES = {
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "conformal_json": "b38d55f9f575ddffd34edeaa5e835d9e10919e6d96a0c284d73c31a072675025",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "centered_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
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
MATRIX_NAMES = ("M", "N", "V", "Gamma", "Omega")
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

    full_tree = ast.parse(FULL_SOURCE.read_text(), filename=str(FULL_SOURCE))
    full_body = [
        node for node in full_tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "high_precision_sector_bases"
    ]
    if len(full_body) != 1:
        raise RuntimeError("audited high-precision sector function is missing")
    exec(
        compile(ast.Module(body=full_body, type_ignores=[]), str(FULL_SOURCE), "exec"),
        globals(),
    )


def rank_record(matrix):
    singular = la.svdvals(matrix)
    epsilon = (
        1000 * MACHINE_EPSILON * max(matrix.shape)
        * max(1.0, float(singular[0]) if len(singular) else 1.0)
    )
    nonzero = singular > 100 * epsilon
    zero = singular < 10 * epsilon
    open_flags = ~(nonzero | zero)
    rank = int(np.sum(nonzero))
    return {
        "singular": singular,
        "epsilon": float(epsilon),
        "rank": rank,
        "zero": int(np.sum(zero)),
        "open": int(np.sum(open_flags)),
    }


def public_rank(record):
    rank = record["rank"]
    singular = record["singular"]
    return {
        "resolved_rank": rank,
        "zero_consistent": record["zero"],
        "open": record["open"],
        "epsilon": serialize_float(record["epsilon"]),
        "minimum_resolved_singular": serialize_float(
            singular[rank - 1] if rank else 0.0
        ),
        "maximum_zero_singular": serialize_float(
            singular[rank] if rank < len(singular) else 0.0
        ),
    }


def build_rigidity(vertices, edges):
    squared_lengths = np.asarray([
        float(np.dot(vertices[left] - vertices[right],
                     vertices[left] - vertices[right]))
        for left, right in edges
    ])
    common_squared_length = float(np.mean(squared_lengths))
    rigidity = np.zeros((720, 480), dtype=float)
    radial = np.zeros((480, 120), dtype=float)
    tangent = np.zeros((480, 480), dtype=float)
    for vertex in range(120):
        radial[4 * vertex:4 * vertex + 4, vertex] = vertices[vertex]
        tangent[4 * vertex:4 * vertex + 4, 4 * vertex:4 * vertex + 4] = (
            np.eye(4) - np.outer(vertices[vertex], vertices[vertex])
        )
    for row, (left, right) in enumerate(edges):
        coefficient = 2 * (vertices[left] - vertices[right]) / common_squared_length
        rigidity[row, 4 * left:4 * left + 4] = coefficient
        rigidity[row, 4 * right:4 * right + 4] = -coefficient
    return rigidity, radial, tangent, squared_lengths, common_squared_length


def cross_label(value, epsilon):
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


print("=" * 78)
print("600-CELL RIGIDITY / SELF-STRESS DYNAMIC CLOSURE GATE")
print("=" * 78)

hashes = {
    "commons": sha256(COMMONS_SOURCE),
    "conformal_json": sha256(CONFORMAL_JSON),
    "conformal_source": sha256(CONFORMAL_SOURCE),
    "centered_json": sha256(CENTERED_JSON),
    "centered_npz": sha256(CENTERED_NPZ),
    "full_source": sha256(FULL_SOURCE),
    "rank_source": sha256(RANK_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
}
centered = json.loads(CENTERED_JSON.read_text())
conformal = json.loads(CONFORMAL_JSON.read_text())
source_npz = np.load(CENTERED_NPZ)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and centered["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and centered["passed"] == centered["tests"] == 7
    and centered["numeric_archive_arrays"] == len(source_npz.files) == 560
    and conformal["outcome"] == "CONFORMAL_MAXIMAL_MINORITY_CERTIFIED"
    and conformal["passed"] == conformal["tests"] == 11
    and all(
        tuple(item["irrep_dimension"] for item in centered["parities"][parity])
        == DIMENSIONS
        for parity in PARITIES
    )
)
check("all preregistered inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_rigidity_york", GEOMETRY_SOURCE
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
vertices, adjacency, _ = build_600cell()
vertices = vertices / np.linalg.norm(vertices, axis=1)[:, None]
coordinate_control_ok = bool(
    vertices.shape == (120, 4)
    and np.max(np.abs(vertices @ vertices.T - np.eye(120))[np.eye(120, dtype=bool)]) < 1e-12
    and np.array_equal(adjacency.astype(np.int8), incidences["even"]["adjacency"])
    and np.array_equal(adjacency.astype(np.int8), incidences["odd"]["adjacency"])
)
check("the embedded unit-quaternion vertices reproduce the literal edge graph", coordinate_control_ok)

rigidity = {}
for parity in PARITIES:
    result = build_rigidity(vertices, groups[parity]["edge_order"])
    rigidity[parity] = {
        "R": result[0],
        "J": result[1],
        "P": result[2],
        "lengths": result[3],
        "L2": result[4],
    }
    rigidity[parity]["D"] = rigidity[parity]["R"] @ rigidity[parity]["P"]

odd_edge_index = {
    edge: index for index, edge in enumerate(groups["odd"]["edge_order"])
}
even_to_odd = np.asarray([
    odd_edge_index[edge] for edge in groups["even"]["edge_order"]
], dtype=int)
row_rigidity_residual = float(la.norm(
    rigidity["odd"]["R"][even_to_odd] - rigidity["even"]["R"], 2
))
length_spread = max(float(np.ptp(rigidity[p]["lengths"])) for p in PARITIES)
row_floor = (
    1000 * MACHINE_EPSILON * 720
    * max(1.0, float(la.norm(rigidity["even"]["R"], 2)))
)
row_covariance_ok = bool(
    sorted(even_to_odd.tolist()) == list(range(720))
    and row_rigidity_residual <= 10 * row_floor
    and length_spread <= 10 * row_floor
)
check(
    "both schedule rigidity matrices are related by the exact edge permutation",
    row_covariance_ok,
    f"operator_residual={row_rigidity_residual:.3e}, length_spread={length_spread:.3e}",
)

global_rank_records = {}
radial_residuals = {}
global_controls_ok = True
for parity in PARITIES:
    r_matrix = rigidity[parity]["R"]
    d_matrix = rigidity[parity]["D"]
    c_matrix = incidences[parity]["incidence"].astype(float)
    rank_r = rank_record(r_matrix)
    rank_d = rank_record(d_matrix)
    rank_cd = rank_record(np.c_[c_matrix, d_matrix])
    radial_residual = float(la.norm(
        r_matrix @ rigidity[parity]["J"] - c_matrix, 2
    ))
    radial_epsilon = max(rank_r["epsilon"], rank_record(c_matrix)["epsilon"])
    intersection = 120 + rank_d["rank"] - rank_cd["rank"]
    stress_dimension = 720 - rank_r["rank"]
    control = bool(
        rank_r["rank"] == 470 and rank_r["open"] == 0
        and rank_d["rank"] == 354 and rank_d["open"] == 0
        and rank_cd["rank"] == 470 and rank_cd["open"] == 0
        and intersection == 4
        and stress_dimension == 250
        and radial_residual <= 10 * radial_epsilon
    )
    global_controls_ok &= control
    radial_residuals[parity] = radial_residual
    global_rank_records[parity] = {
        "R": public_rank(rank_r),
        "D_tangent": public_rank(rank_d),
        "C_plus_D": public_rank(rank_cd),
        "conformal_tangent_intersection_dimension": intersection,
        "self_stress_dimension": stress_dimension,
        "radial_identity_residual": serialize_float(radial_residual),
        "edge_squared_length_mean": serialize_float(rigidity[parity]["L2"]),
        "edge_squared_length_spread": serialize_float(
            np.ptp(rigidity[parity]["lengths"])
        ),
    }
check(
    "the literal matrices reproduce the theorem ranks 470/354/4/250",
    global_controls_ok,
    str({p: {k: v["resolved_rank"] for k, v in global_rank_records[p].items()
             if isinstance(v, dict) and "resolved_rank" in v}
         for p in PARITIES}),
)

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
    "the same seven high-precision minimal sectors are reconstructed",
    sector_basis_ok,
    "max_residual=" + mp.nstr(maximum_basis_residual, 5),
)

sector_carriers = {parity: [] for parity in PARITIES}
carrier_open = False
weighted_ranks = {parity: 0 for parity in PARITIES}
weighted_stresses = {parity: 0 for parity in PARITIES}
for parity in PARITIES:
    for sector_index, sector in enumerate(sector_data):
        dimension = sector["dimension"]
        n = 30 * dimension
        edge_basis = np.kron(
            np.eye(30, dtype=np.complex128), mp_to_numpy(sector["basis"])
        )
        compressed = edge_basis.conj().T @ rigidity[parity]["R"]
        left, singular, _ = la.svd(compressed, full_matrices=True)
        record = rank_record(compressed)
        rank = record["rank"]
        stress_dimension = n - rank
        range_basis = left[:, :rank]
        stress_basis = left[:, rank:]
        projector_residual = float(la.norm(
            compressed - range_basis @ (range_basis.conj().T @ compressed), 2
        ))
        if rank and singular[rank - 1] > 2 * record["epsilon"]:
            eta_s = float(
                2 * record["epsilon"]
                / (singular[rank - 1] - 2 * record["epsilon"])
                + 1000 * MACHINE_EPSILON * n
            )
        else:
            eta_s = math.inf
        sector_resolved = bool(
            record["open"] == 0
            and projector_residual <= 10 * record["epsilon"]
            and math.isfinite(eta_s)
        )
        carrier_open |= not sector_resolved
        weighted_ranks[parity] += dimension * rank
        weighted_stresses[parity] += dimension * stress_dimension
        sector_carriers[parity].append({
            "sector_index": sector_index,
            "dimension": dimension,
            "range_basis": range_basis,
            "stress_basis": stress_basis,
            "eta_S": eta_s,
            "record": {
                "sector_index": sector_index,
                "irrep_dimension": dimension,
                "position_dimension": n,
                "rigidity": public_rank(record),
                "self_stress_dimension_minimal": stress_dimension,
                "self_stress_dimension_weighted": dimension * stress_dimension,
                "projector_residual": serialize_float(projector_residual),
                "subspace_error_eta_S": serialize_float(eta_s),
            },
        })

restoration_ok = all(
    weighted_ranks[p] == 470 and weighted_stresses[p] == 250
    for p in PARITIES
)
check(
    "blind minimal ranks restore exactly 470 rigidity and 250 self-stress dimensions",
    restoration_ok,
    f"ranks={weighted_ranks}, stresses={weighted_stresses}, open={carrier_open}",
)

dynamic_records = {parity: [] for parity in PARITIES}
internal_restricted = {parity: {} for parity in PARITIES}
classification_counts = Counter()
any_required_nonzero = False
any_required_open = False
all_required_zero = True
all_dynamic_finite = True

for parity in PARITIES:
    for sector_index, dimension in enumerate(DIMENSIONS):
        n = 30 * dimension
        range_basis = sector_carriers[parity][sector_index]["range_basis"]
        stress_basis = sector_carriers[parity][sector_index]["stress_basis"]
        eta_s = sector_carriers[parity][sector_index]["eta_S"]
        variants = {}
        internal_restricted[parity][sector_index] = {}
        for variant in VARIANTS:
            matrices = {}
            internal_restricted[parity][sector_index][variant] = {}
            for name in MATRIX_NAMES:
                prefix = f"{parity}_sector{sector_index}_{variant}_{name}"
                midpoint = np.asarray(source_npz[f"{prefix}_midpoint"])
                stored_radius = np.asarray(source_npz[f"{prefix}_radii"])
                entry_radius = component_reenclosure_radii(midpoint, stored_radius)
                norm_x = float(la.norm(midpoint, 2))
                epsilon_x = float(
                    la.norm(entry_radius, "fro")
                    + 1000 * MACHINE_EPSILON * n * max(1.0, norm_x)
                )
                epsilon_cross = float(epsilon_x + 2 * eta_s * norm_x)
                left_cross = range_basis.conj().T @ midpoint @ stress_basis
                right_cross = stress_basis.conj().T @ midpoint @ range_basis
                left_norm = operator_norm(left_cross)
                right_norm = operator_norm(right_cross)
                left_label = cross_label(left_norm, epsilon_cross)
                right_label = cross_label(right_norm, epsilon_cross)
                restricted = stress_basis.conj().T @ midpoint @ stress_basis
                restricted_singular = la.svdvals(restricted)
                internal_restricted[parity][sector_index][variant][name] = {
                    "singular": restricted_singular,
                    "epsilon": epsilon_cross,
                }
                required_labels = []
                if name in ("M", "N", "V"):
                    required_labels.extend((left_label, right_label))
                else:
                    required_labels.append(left_label)
                for label in required_labels:
                    classification_counts[label] += 1
                    any_required_nonzero |= label == "NONZERO_RESOLVED"
                    any_required_open |= label == "OPEN"
                    all_required_zero &= label == "ZERO_CONSISTENT"
                all_dynamic_finite &= bool(
                    math.isfinite(left_norm)
                    and math.isfinite(right_norm)
                    and math.isfinite(epsilon_cross)
                    and np.all(np.isfinite(restricted_singular))
                )
                matrices[name] = {
                    "matrix_norm": serialize_float(norm_x),
                    "matrix_error": serialize_float(epsilon_x),
                    "cross_error": serialize_float(epsilon_cross),
                    "left_cross_norm": serialize_float(left_norm),
                    "left_cross_label": left_label,
                    "right_cross_norm": serialize_float(right_norm),
                    "right_cross_label": right_label,
                    "normalized_self_stress_leakage": (
                        serialize_float(left_norm)
                        if name in ("Gamma", "Omega") else None
                    ),
                    "self_stress_block_minimum_singular": serialize_float(
                        restricted_singular[-1] if len(restricted_singular) else 0.0
                    ),
                    "self_stress_block_maximum_singular": serialize_float(
                        restricted_singular[0] if len(restricted_singular) else 0.0
                    ),
                }
            variants[variant] = matrices
        dynamic_records[parity].append({
            **sector_carriers[parity][sector_index]["record"],
            "variants": variants,
        })

required_classifications = 2 * 7 * 4 * (2 * 3 + 2)
check(
    "all preregistered dynamic cross blocks are finite and classified",
    all_dynamic_finite and sum(classification_counts.values()) == required_classifications,
    f"counts={dict(classification_counts)}",
)

schedule_comparisons = []
schedule_counts = Counter()
schedule_finite = True
for sector_index, dimension in enumerate(DIMENSIONS):
    for variant in VARIANTS:
        for name in ("Gamma", "Omega"):
            left = internal_restricted["even"][sector_index][variant][name]
            right = internal_restricted["odd"][sector_index][variant][name]
            distance = float(np.max(np.abs(
                left["singular"] - right["singular"]
            ))) if len(left["singular"]) else 0.0
            binary_floor = (
                1000 * MACHINE_EPSILON * max(1, len(left["singular"]))
                * max(
                    1.0,
                    float(left["singular"][0]) if len(left["singular"]) else 0.0,
                    float(right["singular"][0]) if len(right["singular"]) else 0.0,
                )
            )
            epsilon = float(left["epsilon"] + right["epsilon"] + binary_floor)
            label = schedule_label(distance, epsilon)
            schedule_counts[label] += 1
            schedule_finite &= math.isfinite(distance) and math.isfinite(epsilon)
            schedule_comparisons.append({
                "sector_index": sector_index,
                "irrep_dimension": dimension,
                "variant": variant,
                "operator": name,
                "maximum_singular_spectrum_distance": serialize_float(distance),
                "comparison_error": serialize_float(epsilon),
                "label": label,
            })
check(
    "all 56 normalized self-stress schedule comparisons are finite",
    schedule_finite and len(schedule_comparisons) == 56,
    str(dict(schedule_counts)),
)

controls_ok = bool(
    provenance_ok
    and geometry_import_ok
    and coordinate_control_ok
    and row_covariance_ok
    and global_controls_ok
    and sector_basis_ok
    and restoration_ok
)
if not controls_ok:
    outcome = "RIGIDITY_YORK_CONTROL_FAILED"
elif carrier_open:
    outcome = "RIGIDITY_YORK_CARRIER_OPEN"
elif any_required_nonzero:
    outcome = "RIGIDITY_YORK_DECOUPLING_REFUTED"
elif any_required_open or not all_required_zero:
    outcome = "RIGIDITY_YORK_DECOUPLING_OPEN"
else:
    outcome = "RIGIDITY_YORK_DECOUPLING_CERTIFIED"

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "canonical_carrier_candidates": 1,
    "known_theorem_controls": global_rank_records,
    "weighted_sector_ranks": weighted_ranks,
    "weighted_self_stress_dimensions": weighted_stresses,
    "basis_control": {
        "irrep_dimensions": list(DIMENSIONS),
        "maximum_high_precision_residual": mp.nstr(maximum_basis_residual, 70),
    },
    "parities": dynamic_records,
    "required_cross_classification_counts": dict(classification_counts),
    "schedule_comparisons": schedule_comparisons,
    "schedule_label_counts": dict(schedule_counts),
    "classification": {
        "carrier_open": carrier_open,
        "all_required_zero_consistent": all_required_zero,
        "any_required_nonzero_resolved": any_required_nonzero,
        "any_required_open": any_required_open,
        "self_stress_is_not_declared_physical": True,
        "exact_constraint_quotient_derived": False,
    },
    "continuum_target_loaded": False,
    "polarization_target_loaded": False,
    "speed_target_loaded": False,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("SCIENTIFIC OUTCOME:", outcome)
print("required cross labels:", dict(classification_counts))
print("schedule labels:", dict(schedule_counts))
print(f"{passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

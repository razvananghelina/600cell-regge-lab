#!/usr/bin/env python3
"""Canonical vertex-conformal restriction of the centered Regge kinetic form.

Prior-art commit: 96dd1ff.
Protocol commit: 4d23b25.
Geometry-only orbit-order correction: 298035f.
No continuum spectrum, polarization count, speed or particle target is loaded.
"""

import ast
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
CENTERED_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
CENTERED_SOURCE = HERE / "verify_gravity_600cell_dust_centered_jacobi.py"
FULL_JSON = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
OUTPUT = HERE / "gravity_600cell_dust_conformal_supermetric.json"

PRIOR_ART_COMMIT = "96dd1ff"
PROTOCOL_COMMIT = "4d23b25"
PROTOCOL_CORRECTION_COMMIT = "298035f"
EXPECTED_HASHES = {
    "centered_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "centered_source": "359b8d7642746c2dc22e304353e3b83104874badd86755de4f8f9e6f25e56a20",
    "full_json": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
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


def serialize_complex(value):
    return {
        "real": serialize_float(np.real(value)),
        "imaginary": serialize_float(np.imag(value)),
    }


def mp_frobenius(matrix):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in matrix))


def mp_submatrix(matrix, rows, columns):
    return mp.matrix([[matrix[row, column] for column in columns] for row in rows])


def cluster_sorted(values, tolerance=mp.mpf("1e-70")):
    clusters = []
    for index, value in enumerate(values):
        if not clusters or abs(value - values[clusters[-1][0]]) > tolerance:
            clusters.append([index])
        else:
            clusters[-1].append(index)
    return clusters


def load_audited_sector_function():
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


def orbit_sort_key(orbit, phase):
    representative = min(orbit)
    logical = tuple(vertex % 120 for vertex in representative)
    phase_pair = tuple(sorted(phase[vertex] for vertex in logical))
    return phase_pair, tuple(sorted(orbit))


def edge_image(action, edge):
    return tuple(sorted((int(action[edge[0]]), int(action[edge[1]]))))


def group_data(model, gro):
    actions = tuple(sorted(
        tuple(int(value) for value in action) for action in model["stabilizer"]
    ))
    action_arrays = tuple(np.asarray(action, dtype=np.int16) for action in actions)
    action_index = {action: index for index, action in enumerate(actions)}
    identity = action_index[tuple(range(120))]
    table = np.empty((24, 24), dtype=np.int16)
    for left, action in enumerate(action_arrays):
        for right, other in enumerate(action_arrays):
            table[left, right] = action_index[tuple(action[other])]
    inverses = tuple(
        next(
            right for right in range(24)
            if table[left, right] == identity and table[right, left] == identity
        )
        for left in range(24)
    )
    orders = []
    for element in range(24):
        product = identity
        order = 0
        while True:
            order += 1
            product = int(table[product, element])
            if product == identity:
                break
        orders.append(order)
    unseen = set(range(24))
    classes = []
    while unseen:
        seed = min(unseen)
        conjugacy_class = tuple(sorted({
            int(table[table[group, seed], inverses[group]])
            for group in range(24)
        }))
        classes.append(conjugacy_class)
        unseen -= set(conjugacy_class)
    classes = tuple(sorted(
        classes, key=lambda item: (len(item), orders[item[0]], item)
    ))

    old_orbits = tuple(sorted(
        gro.orbit_partition(model["old_edges"], model["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, model["phase"]),
    ))
    orbit_edges = []
    for orbit in old_orbits:
        seed = min(orbit)
        ordered = tuple(edge_image(action, seed) for action in action_arrays)
        if len(set(ordered)) != 24 or set(ordered) != set(orbit):
            raise RuntimeError("old boundary edge orbit is not free regular")
        orbit_edges.append(ordered)
    edge_order = tuple(edge for orbit in orbit_edges for edge in orbit)

    vertex_unseen = set(range(120))
    vertex_orbits = []
    while vertex_unseen:
        seed = min(vertex_unseen)
        orbit = frozenset(int(action[seed]) for action in action_arrays)
        vertex_orbits.append(orbit)
        vertex_unseen -= orbit

    return {
        "actions": action_arrays,
        "table": table,
        "identity": identity,
        "orders": tuple(orders),
        "classes": classes,
        "old_orbits": old_orbits,
        "orbit_edges": tuple(orbit_edges),
        "edge_order": edge_order,
        "vertex_orbits": tuple(vertex_orbits),
    }


def incidence_data(group):
    edges = group["edge_order"]
    edge_index = {edge: index for index, edge in enumerate(edges)}
    incidence = np.zeros((720, 120), dtype=np.int8)
    adjacency = np.zeros((120, 120), dtype=np.int8)
    for row, (left, right) in enumerate(edges):
        incidence[row, left] = 1
        incidence[row, right] = 1
        adjacency[left, right] = adjacency[right, left] = 1

    sparse_incidence = sp.csr_matrix(incidence)
    equivariant = True
    for action in group["actions"]:
        vertex_rows = np.asarray(action, dtype=int)
        vertex_columns = np.arange(120)
        vertex_permutation = sp.csr_matrix(
            (np.ones(120, dtype=np.int8), (vertex_rows, vertex_columns)),
            shape=(120, 120),
        )
        edge_rows = np.asarray([
            edge_index[edge_image(action, edge)] for edge in edges
        ], dtype=int)
        edge_columns = np.arange(720)
        edge_permutation = sp.csr_matrix(
            (np.ones(720, dtype=np.int8), (edge_rows, edge_columns)),
            shape=(720, 720),
        )
        residual = sparse_incidence @ vertex_permutation - edge_permutation @ sparse_incidence
        if residual.nnz != 0:
            equivariant = False
            break

    reached = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for neighbour in np.flatnonzero(adjacency[vertex]):
            neighbour = int(neighbour)
            if neighbour not in reached:
                reached.add(neighbour)
                frontier.append(neighbour)
    triangle_count_sixfold = int(np.trace(
        adjacency.astype(np.int64) @ adjacency.astype(np.int64)
        @ adjacency.astype(np.int64)
    ))
    gram = incidence.T.astype(np.int64) @ incidence.astype(np.int64)
    numerical_singular = la.svdvals(incidence.astype(float))
    return {
        "incidence": incidence,
        "adjacency": adjacency,
        "equivariant": equivariant,
        "connected": len(reached) == 120,
        "triangle_count": triangle_count_sixfold // 6,
        "gram_identity": np.array_equal(
            gram, 12 * np.eye(120, dtype=np.int64) + adjacency
        ),
        "numerical_rank": int(np.linalg.matrix_rank(incidence.astype(float))),
        "minimum_singular": float(numerical_singular[-1]),
    }


def mp_to_numpy(matrix):
    return np.array([
        [complex(float(mp.re(matrix[row, column])), float(mp.im(matrix[row, column])))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


def component_reenclosure_radii(midpoint, stored):
    real_half_ulp = 0.5 * np.abs(np.spacing(np.real(midpoint)))
    imaginary_half_ulp = 0.5 * np.abs(np.spacing(np.imag(midpoint)))
    return np.hypot(stored + real_half_ulp, stored + imaginary_half_ulp)


def sign_counts(values, epsilon):
    positive = values > 100 * epsilon
    negative = values < -100 * epsilon
    zero = np.abs(values) < 10 * epsilon
    open_flags = ~(positive | negative | zero)
    return {
        "positive": int(np.sum(positive)),
        "negative": int(np.sum(negative)),
        "zero_consistent": int(np.sum(zero)),
        "open": int(np.sum(open_flags)),
    }


def threshold_label(value, epsilon, identified, separated):
    if not math.isfinite(value) or not math.isfinite(epsilon):
        return "OPEN"
    if value <= 10 * epsilon:
        return identified
    if value > 100 * epsilon:
        return separated
    return "OPEN"


print("=" * 78)
print("CANONICAL 600-CELL CONFORMAL SUPERMETRIC GATE")
print("=" * 78)

hashes = {
    "centered_json": sha256(CENTERED_JSON),
    "centered_npz": sha256(CENTERED_NPZ),
    "centered_source": sha256(CENTERED_SOURCE),
    "full_json": sha256(FULL_JSON),
    "full_source": sha256(FULL_SOURCE),
    "rank_source": sha256(RANK_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
}
centered = json.loads(CENTERED_JSON.read_text())
full = json.loads(FULL_JSON.read_text())
source_npz = np.load(CENTERED_NPZ)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and centered["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and centered["passed"] == centered["tests"] == 7
    and centered["numeric_archive_arrays"] == len(source_npz.files) == 560
    and centered["numeric_archive_sha256"] == EXPECTED_HASHES["centered_npz"]
    and full["outcome"] == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and full["passed"] == full["tests"] == 19
    and all(
        tuple(item["irrep_dimension"] for item in centered["parities"][parity])
        == DIMENSIONS
        for parity in PARITIES
    )
)
check("all preregistered inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_conformal_supermetric", GEOMETRY_SOURCE
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
check("the direct literal one-slab geometry retains all 43 certificates", geometry_import_ok)

load_audited_sector_function()
groups = {parity: group_data(gro.models[parity], gro) for parity in PARITIES}
incidences = {parity: incidence_data(groups[parity]) for parity in PARITIES}

same_actions = all(
    np.array_equal(left, right)
    for left, right in zip(groups["even"]["actions"], groups["odd"]["actions"])
)
same_edge_set = set(groups["even"]["edge_order"]) == set(groups["odd"]["edge_order"])
each_order_bijective = all(
    len(group["edge_order"]) == len(set(group["edge_order"])) == 720
    and len(group["old_orbits"]) == 30
    and all(len(orbit) == 24 for orbit in group["old_orbits"])
    and len(group["vertex_orbits"]) == 5
    and all(len(orbit) == 24 for orbit in group["vertex_orbits"])
    for group in groups.values()
)
odd_edge_index = {
    edge: index for index, edge in enumerate(groups["odd"]["edge_order"])
}
even_to_odd = np.asarray([
    odd_edge_index[edge] for edge in groups["even"]["edge_order"]
], dtype=int)
row_permutation_ok = bool(
    sorted(even_to_odd.tolist()) == list(range(720))
    and np.array_equal(
        incidences["odd"]["incidence"][even_to_odd],
        incidences["even"]["incidence"],
    )
)
orbit_control_ok = same_actions and same_edge_set and each_order_bijective and row_permutation_ok
check(
    "both schedules give the same literal carrier with an exact row permutation",
    orbit_control_ok,
    f"row_sequences_equal={groups['even']['edge_order'] == groups['odd']['edge_order']}",
)

literal_incidence_ok = all(
    data["incidence"].shape == (720, 120)
    and np.all(np.sum(data["incidence"], axis=1) == 2)
    and np.all(np.sum(data["incidence"], axis=0) == 12)
    and data["gram_identity"]
    for data in incidences.values()
)
check("the canonical unsigned incidence matrix has the exact 600-cell local counts", literal_incidence_ok)

graph_injectivity_ok = all(
    data["connected"]
    and data["triangle_count"] == 1200
    and data["numerical_rank"] == 120
    and data["minimum_singular"] > 0
    for data in incidences.values()
)
check(
    "connectedness plus odd cycles prove exact conformal-map injectivity",
    graph_injectivity_ok,
    "triangles=" + str({p: incidences[p]["triangle_count"] for p in PARITIES}),
)

equivariance_ok = all(data["equivariant"] for data in incidences.values())
check("all 24 binary-tetrahedral actions commute exactly with incidence", equivariance_ok)

sector_data, sector_controls = high_precision_sector_bases(groups["even"])
reconstructed_signature = tuple(
    (
        sector["dimension"],
        complex(float(mp.re(sector["old_central_eigenvalue"])),
                float(mp.im(sector["old_central_eigenvalue"]))),
        sector["splitter"],
    )
    for sector in sector_data
)
signature_matches = True
for parity in PARITIES:
    public = full["parities"][parity]["sectors"]
    if len(public) != 7:
        signature_matches = False
        continue
    for index, (dimension, center, splitter) in enumerate(reconstructed_signature):
        frozen_center = complex(
            float(public[index]["central_eigenvalue"]["real"]),
            float(public[index]["central_eigenvalue"]["imaginary"]),
        )
        signature_matches &= bool(
            public[index]["sector_index"] == index
            and public[index]["stored_sector_index"] == index
            and public[index]["dimension"] == dimension
            and public[index]["splitter_group_index"] == splitter
            and abs(frozen_center - center) < 1e-60
        )
maximum_basis_residual = max(
    value for key, value in sector_controls.items() if key.startswith("maximum_")
)
sector_basis_ok = bool(
    signature_matches
    and tuple(sector["dimension"] for sector in sector_data) == DIMENSIONS
    and sector_controls["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
    and maximum_basis_residual < mp.mpf("1e-70")
)
check(
    "seven high-precision minimal sectors match the frozen source by full signature",
    sector_basis_ok,
    "dims=" + str([sector["dimension"] for sector in sector_data])
    + ", max_residual=" + mp.nstr(maximum_basis_residual, 5),
)

sector_images = {parity: [] for parity in PARITIES}
all_sector_ranks_resolved = True
for parity in PARITIES:
    incidence = incidences[parity]["incidence"].astype(np.complex128)
    for sector_index, sector in enumerate(sector_data):
        dimension = sector["dimension"]
        basis = mp_to_numpy(sector["basis"])
        edge_basis = np.kron(np.eye(30, dtype=np.complex128), basis)
        compressed = edge_basis.conj().T @ incidence
        left, singular, _ = la.svd(compressed, full_matrices=False)
        epsilon_c = (
            100 * MACHINE_EPSILON * max(compressed.shape)
            * max(1.0, float(singular[0]))
        )
        expected_rank = 5 * dimension
        nonzero = int(np.sum(singular > 100 * epsilon_c))
        zero = int(np.sum(singular < 10 * epsilon_c))
        open_count = len(singular) - nonzero - zero
        conformal_basis = left[:, :expected_rank]
        projector_residual = la.norm(
            compressed - conformal_basis @ (conformal_basis.conj().T @ compressed),
            2,
        )
        rank_resolved = bool(
            nonzero == expected_rank
            and zero == len(singular) - expected_rank
            and open_count == 0
            and projector_residual <= 10 * epsilon_c
        )
        all_sector_ranks_resolved &= rank_resolved
        sector_images[parity].append({
            "sector_index": sector_index,
            "dimension": dimension,
            "compressed": compressed,
            "basis": conformal_basis,
            "rank_resolved": rank_resolved,
            "record": {
                "sector_index": sector_index,
                "irrep_dimension": dimension,
                "position_dimension": 30 * dimension,
                "expected_conformal_rank": expected_rank,
                "resolved_nonzero_singular_values": nonzero,
                "zero_consistent_singular_values": zero,
                "open_singular_values": open_count,
                "epsilon_C": serialize_float(epsilon_c),
                "minimum_resolved_singular": serialize_float(singular[expected_rank - 1]),
                "maximum_zero_singular": serialize_float(
                    singular[expected_rank] if expected_rank < len(singular) else 0.0
                ),
                "projector_residual": serialize_float(projector_residual),
            },
        })
check(
    "all fourteen minimal conformal images have the preregistered rank 5d",
    len(sector_images["even"]) == len(sector_images["odd"]) == 7
    and all(
        np.all(np.isfinite(item["compressed"]))
        for parity in PARITIES for item in sector_images[parity]
    ),
    "resolved=" + str(all_sector_ranks_resolved),
)

matrix_records = {parity: [] for parity in PARITIES}
internal = {parity: {} for parity in PARITIES}
all_full_exact = True
all_restrictions_positive = True
restriction_refuted = False
restriction_open = False
full_refuted = False
all_matrix_arithmetic_finite = True

for parity in PARITIES:
    for sector_index, dimension in enumerate(DIMENSIONS):
        n = 30 * dimension
        expected_conformal = 5 * dimension
        conformal_basis = sector_images[parity][sector_index]["basis"]
        variants = {}
        internal[parity][sector_index] = {}
        for variant in VARIANTS:
            prefix = f"{parity}_sector{sector_index}_{variant}_M"
            midpoint = np.asarray(source_npz[f"{prefix}_midpoint"])
            stored_radius = np.asarray(source_npz[f"{prefix}_radii"])
            entry_radius = component_reenclosure_radii(midpoint, stored_radius)
            hermitian = (midpoint + midpoint.conj().T) / 2
            hermitian = (hermitian + hermitian.conj().T) / 2
            hermitian_radius = (entry_radius + entry_radius.T) / 2
            arithmetic_floor = (
                1000 * MACHINE_EPSILON * n
                * max(1.0, float(la.norm(hermitian, 2)))
            )
            epsilon_h = float(la.norm(hermitian_radius, "fro") + arithmetic_floor)
            epsilon_g = float(epsilon_h + arithmetic_floor)
            restricted = conformal_basis.conj().T @ hermitian @ conformal_basis
            restricted = (restricted + restricted.conj().T) / 2
            full_values = la.eigvalsh(hermitian)
            restricted_values = la.eigvalsh(restricted)
            full_counts = sign_counts(full_values, epsilon_h)
            restricted_counts = sign_counts(restricted_values, epsilon_g)
            full_expected = bool(full_counts == {
                "positive": expected_conformal,
                "negative": 25 * dimension,
                "zero_consistent": 0,
                "open": 0,
            })
            restricted_positive = bool(restricted_counts == {
                "positive": expected_conformal,
                "negative": 0,
                "zero_consistent": 0,
                "open": 0,
            })
            all_full_exact &= full_expected
            all_restrictions_positive &= restricted_positive
            restriction_refuted |= restricted_counts["negative"] > 0
            restriction_open |= bool(
                restricted_counts["zero_consistent"] > 0
                or restricted_counts["open"] > 0
            )
            full_refuted |= bool(
                full_counts["positive"] > expected_conformal
                or full_counts["negative"] > 25 * dimension
            )
            if not full_expected and not full_refuted:
                restriction_open = True
            all_matrix_arithmetic_finite &= bool(
                np.all(np.isfinite(full_values))
                and np.all(np.isfinite(restricted_values))
                and math.isfinite(epsilon_h)
                and math.isfinite(epsilon_g)
            )
            internal[parity][sector_index][variant] = {
                "H": hermitian,
                "G": restricted,
                "full_values": full_values,
                "restricted_values": restricted_values,
                "epsilon_H": epsilon_h,
                "epsilon_G": epsilon_g,
            }
            variants[variant] = {
                "epsilon_H": serialize_float(epsilon_h),
                "epsilon_G": serialize_float(epsilon_g),
                "full_inertia": full_counts,
                "restricted_inertia": restricted_counts,
                "full_minimum_eigenvalue": serialize_float(full_values[0]),
                "full_maximum_eigenvalue": serialize_float(full_values[-1]),
                "restricted_minimum_eigenvalue": serialize_float(restricted_values[0]),
                "restricted_maximum_eigenvalue": serialize_float(restricted_values[-1]),
                "restricted_minimum_margin_error_units": serialize_float(
                    restricted_values[0] / epsilon_g
                ),
                "restricted_condition": serialize_float(
                    restricted_values[-1] / restricted_values[0]
                    if restricted_values[0] != 0 else math.inf
                ),
            }
        matrix_records[parity].append({
            **sector_images[parity][sector_index]["record"],
            "variants": variants,
        })

check(
    "all 56 full and restricted Hermitian eigencensuses are finite and classified",
    all_matrix_arithmetic_finite
    and sum(len(item["variants"]) for p in PARITIES for item in matrix_records[p]) == 56,
    f"full_exact={all_full_exact}, restricted_positive={all_restrictions_positive}",
)

secondary_records = {parity: [] for parity in PARITIES}
secondary_finite = True
for parity in PARITIES:
    for sector_index, dimension in enumerate(DIMENSIONS):
        item = internal[parity][sector_index]["operational_primary"]
        hermitian = item["H"]
        values, vectors = la.eigh(hermitian)
        positive_indices = np.flatnonzero(values > 100 * item["epsilon_H"])
        negative_indices = np.flatnonzero(values < -100 * item["epsilon_H"])
        conformal_basis = sector_images[parity][sector_index]["basis"]
        conformal_projector = conformal_basis @ conformal_basis.conj().T
        if len(positive_indices) == 5 * dimension and len(negative_indices) == 25 * dimension:
            positive_basis = vectors[:, positive_indices]
            positive_projector = positive_basis @ positive_basis.conj().T
            projector_distance = float(la.norm(
                conformal_projector - positive_projector, 2
            ))
            maximum_angle = math.asin(min(1.0, max(0.0, projector_distance)))
            gap = float(values[positive_indices[0]] - values[negative_indices[-1]])
            if gap > 2 * item["epsilon_H"]:
                epsilon_p = float(
                    2 * item["epsilon_H"] / (gap - 2 * item["epsilon_H"])
                    + 1000 * MACHINE_EPSILON * (30 * dimension)
                )
            else:
                epsilon_p = math.inf
        else:
            projector_distance = math.nan
            maximum_angle = math.nan
            gap = math.nan
            epsilon_p = math.inf
        leakage = float(la.norm(
            (np.eye(30 * dimension) - conformal_projector)
            @ hermitian @ conformal_basis,
            2,
        ))
        epsilon_l = float(
            item["epsilon_H"]
            + 1000 * MACHINE_EPSILON * (30 * dimension)
            * max(1.0, float(la.norm(hermitian, 2)))
        )
        projector_label = threshold_label(
            projector_distance, epsilon_p, "SPECTRAL_IDENTIFIED", "SPECTRAL_SEPARATED"
        )
        invariance_label = threshold_label(
            leakage, epsilon_l, "INVARIANT_IDENTIFIED", "INVARIANCE_SEPARATED"
        )
        secondary_finite &= bool(
            math.isfinite(leakage)
            and math.isfinite(epsilon_l)
            and (math.isfinite(projector_distance) or math.isinf(epsilon_p))
        )
        secondary_records[parity].append({
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "projector_distance": serialize_float(projector_distance),
            "maximum_principal_angle_radians": serialize_float(maximum_angle),
            "spectral_gap": serialize_float(gap),
            "projector_error": serialize_float(epsilon_p),
            "projector_label": projector_label,
            "invariance_leakage": serialize_float(leakage),
            "invariance_error": serialize_float(epsilon_l),
            "invariance_label": invariance_label,
        })
check(
    "fourteen operational spectral and invariance diagnostics are complete",
    secondary_finite and all(len(secondary_records[p]) == 7 for p in PARITIES),
)

schedule_comparisons = []
schedule_finite = True
for sector_index, dimension in enumerate(DIMENSIONS):
    for variant in VARIANTS:
        left = internal["even"][sector_index][variant]
        right = internal["odd"][sector_index][variant]
        distance = float(np.max(np.abs(
            left["restricted_values"] - right["restricted_values"]
        )))
        epsilon = float(left["epsilon_G"] + right["epsilon_G"])
        label = threshold_label(
            distance, epsilon, "SCHEDULE_ROBUST", "SCHEDULE_DEPENDENT"
        )
        schedule_finite &= math.isfinite(distance) and math.isfinite(epsilon)
        schedule_comparisons.append({
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "variant": variant,
            "maximum_ordered_restricted_eigenvalue_distance": serialize_float(distance),
            "comparison_error": serialize_float(epsilon),
            "label": label,
        })
check(
    "all 28 restricted schedule comparisons are finite and target-free",
    schedule_finite and len(schedule_comparisons) == 28,
)

controls_ok = bool(
    provenance_ok
    and geometry_import_ok
    and orbit_control_ok
    and literal_incidence_ok
    and graph_injectivity_ok
    and equivariance_ok
    and sector_basis_ok
)
if not controls_ok:
    outcome = "CONFORMAL_SUPERMETRIC_CONTROL_FAILED"
elif not all_sector_ranks_resolved:
    outcome = "CONFORMAL_SUPERMETRIC_SECTOR_OPEN"
elif restriction_refuted or full_refuted:
    outcome = "CONFORMAL_MAXIMAL_MINORITY_REFUTED"
elif restriction_open or not all_full_exact or not all_restrictions_positive:
    outcome = "CONFORMAL_SUPERMETRIC_RESTRICTION_OPEN"
else:
    outcome = "CONFORMAL_MAXIMAL_MINORITY_CERTIFIED"

projector_labels = {}
invariance_labels = {}
for parity in PARITIES:
    for item in secondary_records[parity]:
        projector_labels[item["projector_label"]] = (
            projector_labels.get(item["projector_label"], 0) + 1
        )
        invariance_labels[item["invariance_label"]] = (
            invariance_labels.get(item["invariance_label"], 0) + 1
        )

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "input_sha256": hashes,
    "canonical_map_candidates": 1,
    "carrier": {
        "vertices": 120,
        "edges": 720,
        "vertex_degree": 12,
        "triangles": incidences["even"]["triangle_count"],
        "rank_C_exact_by_odd_cycle": 120 if graph_injectivity_ok else None,
        "minimum_numeric_singular_C": serialize_float(
            incidences["even"]["minimum_singular"]
        ),
        "schedule_row_sequences_equal": bool(
            groups["even"]["edge_order"] == groups["odd"]["edge_order"]
        ),
        "schedule_row_permutation_exact": row_permutation_ok,
        "equivariant_all_24": equivariance_ok,
    },
    "basis_control": {
        "irrep_dimensions": list(DIMENSIONS),
        "maximum_high_precision_residual": mp.nstr(maximum_basis_residual, 70),
        "signature_matches_frozen_source": signature_matches,
    },
    "parities": matrix_records,
    "secondary_structural": secondary_records,
    "secondary_label_counts": {
        "projector": projector_labels,
        "invariance": invariance_labels,
    },
    "schedule_comparisons": schedule_comparisons,
    "classification": {
        "all_sector_ranks_resolved": all_sector_ranks_resolved,
        "all_full_inertias_exact": all_full_exact,
        "all_conformal_restrictions_positive": all_restrictions_positive,
        "restriction_refuted": restriction_refuted,
        "restriction_open": restriction_open,
        "full_refuted": full_refuted,
        "primary_statement": (
            "The canonical vertex-conformal image is a maximal subspace "
            "carrying the minority inertia sign on the fixed carrier."
            if outcome == "CONFORMAL_MAXIMAL_MINORITY_CERTIFIED" else None
        ),
        "secondary_is_coordinate_dependent": True,
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
print("secondary projector labels:", projector_labels)
print("secondary invariance labels:", invariance_labels)
print(f"{passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)


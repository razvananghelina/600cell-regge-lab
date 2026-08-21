#!/usr/bin/env python3
"""Geometric adversarial replication of the schedule-0 Coxeter certificate.

Protocol commit: b54f590.
The Coxeter action is reconstructed from ambient R4 reflections and the
blocks from explicit sparse Fourier isometries.  Primary block code is not
loaded.
"""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


SOURCE = HERE / "verify_gravity_600cell_refined_nonhomogeneous_internal_hessian.py"
ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
FROZEN = HERE / "gravity_600cell_refined_nonhomogeneous_internal_hessian.json"
CURVATURE = HERE / "gravity_600cell_refined_local_curvature_mass.json"
PRIMARY = HERE / "gravity_600cell_refined_nonhomogeneous_coxeter_blocks.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_block_adversarial_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_nonhomogeneous_coxeter_blocks_adversarial.json"

PROTOCOL_COMMIT = "b54f590"
EXPECTED_HASHES = {
    "reproducible/verify_gravity_600cell_refined_nonhomogeneous_internal_hessian.py":
        "2a96c8ce466d6d4e9be2cadc4ebf932b4e42eff16fc64fbfb08cd580d680879e",
    "reproducible/gravity_600cell_refined_nonhomogeneous_internal_hessian.json":
        "4a05968c68f8e6a35a1308ddf6114bb19b7106f214bfdcf798e7af2387bddec1",
    "reproducible/gravity_600cell_refined_local_curvature_mass.json":
        "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_blocks.json":
        "640f07a3d13ae3692761243cb62ace3ac2fd38f646d03b7750df05883d3f0267",
    "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_block_adversarial_protocol.md":
        "b9f33f1e95ec52ab607f158760ca90970b7157b128406a7756b5d549ba57da39",
}

PAIR4 = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR4)}
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_TRIANGLES = tuple(combinations(range(5), 3))
TRIANGLE_EDGES = tuple(combinations(range(3), 2))
TAU_TEXT = "0.0102"
DECIMAL_PRECISIONS = (100, 140)
DIFFERENCE_STEPS = ("1e-10", "5e-11", "2.5e-11")
EXPECTED_F = (2640, 17040, 28800, 14400)
EXPECTED_COUNTS = {
    "pentachora": 57600,
    "triangles": 149280,
    "boundary_edges": 34080,
    "internal_edges": 19680,
    "cross_edges": 17040,
    "vertical_edges": 2640,
    "all_edges": 53760,
}
WORD = (0, 1, 2, 3)
ORDER = 30
DIMENSION = 656
MATCH_GATE = 5e-8

tests = 0
passed = 0


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


def ftext(value):
    return format(float(value), ".17e")


def load_source_definitions():
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "ast": ast,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "sha256": sha256,
        "combinations": combinations,
        "permutations": permutations,
        "json": json,
        "Path": Path,
        "sys": sys,
        "mp": mp,
        "np": np,
        "sp": sp,
        "spla": spla,
        "build_600cell": build_600cell,
        "HERE": HERE,
        "ROOT": ROOT,
        "ACTION_SOURCE": ACTION_SOURCE,
        "PAIR4": PAIR4,
        "PAIR_INDEX": PAIR_INDEX,
        "LOCAL_EDGES": LOCAL_EDGES,
        "LOCAL_TRIANGLES": LOCAL_TRIANGLES,
        "TRIANGLE_EDGES": TRIANGLE_EDGES,
        "TAU_TEXT": TAU_TEXT,
        "DECIMAL_PRECISIONS": DECIMAL_PRECISIONS,
        "DIFFERENCE_STEPS": DIFFERENCE_STEPS,
        "EXPECTED_F": EXPECTED_F,
        "EXPECTED_COUNTS": EXPECTED_COUNTS,
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(SOURCE), "exec"), namespace)
    return namespace


def matrix_order(matrix, limit=60, gate=1e-8):
    power = np.eye(len(matrix))
    identity = power.copy()
    for exponent in range(1, limit + 1):
        power = matrix @ power
        if np.linalg.norm(power - identity, ord=np.inf) <= gate:
            return exponent
    return None


def permutation_order(mapping, limit=60):
    identity = np.arange(len(mapping), dtype=np.int32)
    power = identity.copy()
    for exponent in range(1, limit + 1):
        power = mapping[power]
        if np.array_equal(power, identity):
            return exponent
    return None


def cell_centres(vertices, cells):
    result = np.empty((len(cells), 4), dtype=np.float64)
    for index, cell in enumerate(cells):
        value = np.sum(vertices[np.asarray(cell, dtype=np.int32)], axis=0)
        result[index] = value / np.linalg.norm(value)
    return result


def chamber_reflections(centres, base_flag):
    base = centres[np.asarray(base_flag, dtype=np.int32)]
    reflections = []
    normals = []
    for rank in range(4):
        wall = np.delete(base, rank, axis=0)
        _, _, right = np.linalg.svd(wall, full_matrices=True)
        normal = right[-1]
        normal /= np.linalg.norm(normal)
        reflection = np.eye(4) - 2 * np.outer(normal, normal)
        normals.append(normal)
        reflections.append(reflection)
    return tuple(reflections), tuple(normals)


def sequential_product(reflections, word):
    result = np.eye(4)
    for rank in word:
        result = reflections[rank] @ result
    return result


def vertex_permutation(vertices, transformation):
    transformed = (transformation @ vertices.T).T
    dots = transformed @ vertices.T
    mapping = np.argmax(dots, axis=1).astype(np.int32)
    residuals = np.linalg.norm(transformed - vertices[mapping], axis=1)
    return mapping, float(residuals.max()), len(set(map(int, mapping))) == len(vertices)


def induced_cell_permutation(vertex_cells, vertex_mapping):
    index = {tuple(cell): position for position, cell in enumerate(vertex_cells)}
    result = np.empty(len(vertex_cells), dtype=np.int32)
    missing = 0
    for position, cell in enumerate(vertex_cells):
        image = tuple(sorted(int(vertex_mapping[vertex]) for vertex in cell))
        if image not in index:
            missing += 1
            result[position] = -1
        else:
            result[position] = index[image]
    bijective = missing == 0 and len(set(map(int, result))) == len(result)
    return result, bijective, missing


def induced_edge_permutation(source, schedule, cell_mapping, cell_count):
    mapped = schedule["internal_edges"].copy()
    for endpoint in range(2):
        layer = mapped[:, endpoint] >= cell_count
        base = mapped[:, endpoint] % cell_count
        mapped[:, endpoint] = cell_mapping[base] + layer * cell_count
    mapped.sort(axis=1)
    codes = source["encoded_edges"](schedule["internal_edges"], 2 * cell_count)
    mapped_codes = source["encoded_edges"](mapped, 2 * cell_count)
    positions = np.searchsorted(codes, mapped_codes)
    invariant = (
        np.all(positions < len(codes))
        and np.array_equal(codes[positions], mapped_codes)
        and len(set(map(int, positions))) == len(positions)
    )
    return positions.astype(np.int32), invariant


def cycle_coordinates(mapping):
    cycle_id = np.full(len(mapping), -1, dtype=np.int32)
    position = np.full(len(mapping), -1, dtype=np.int16)
    cycles = []
    for seed in range(len(mapping)):
        if cycle_id[seed] >= 0:
            continue
        cycle = []
        value = seed
        while cycle_id[value] < 0:
            cycle_id[value] = len(cycles)
            position[value] = len(cycle)
            cycle.append(value)
            value = int(mapping[value])
        if value != seed:
            raise RuntimeError("cycle traversal did not close at its seed")
        cycles.append(np.asarray(cycle, dtype=np.int32))
    return tuple(cycles), cycle_id, position


def covariance_defects(source, matrix, mapping):
    identity = np.arange(len(mapping), dtype=np.int32)
    power = identity.copy()
    defects = []
    for exponent in range(ORDER):
        if exponent:
            power = mapping[power]
        defects.append(source["sparse_row_norm"](
            matrix - matrix[power, :][:, power]
        ))
    return defects, float(np.mean(defects))


def explicit_fourier_isometry(cycle_id, position, character, sign=-1):
    rows = np.arange(len(cycle_id), dtype=np.int32)
    columns = cycle_id.copy()
    data = np.exp(
        sign * 2j * np.pi * character
        * position.astype(np.float64) / ORDER
    ) / np.sqrt(ORDER)
    return sp.coo_matrix(
        (data, (rows, columns)), shape=(len(rows), DIMENSION)
    ).tocsr()


def sparse_support_products(matrix, cycle_id):
    support_q = sp.coo_matrix(
        (
            np.ones(len(cycle_id)),
            (np.arange(len(cycle_id), dtype=np.int32), cycle_id),
        ),
        shape=(len(cycle_id), DIMENSION),
    ).tocsr()
    support_c = matrix.copy()
    support_c.data = np.ones_like(support_c.data)
    support_counts = support_q.T @ (support_c @ support_q)
    absolute_q = support_q / np.sqrt(ORDER)
    absolute_bound = absolute_q.T @ (abs(matrix) @ absolute_q)
    return support_counts, absolute_bound


def sparse_product_block(matrix, isometry, support_counts, absolute_bound):
    product = matrix.astype(np.complex128) @ isometry
    block_sparse = isometry.conj().T @ product
    block = block_sparse.toarray()
    hermitian_defect = float(np.max(np.sum(
        np.abs(block - block.conj().T), axis=1
    )))
    block = (block + block.conj().T) / 2

    counts = support_counts.toarray()
    absolute = absolute_bound.toarray()
    unit_roundoff = np.finfo(np.float64).eps / 2
    maximum_count = int(np.max(counts))
    gamma = maximum_count * unit_roundoff / (
        1 - maximum_count * unit_roundoff
    )
    roundoff_bound = float(
        (gamma + 20 * unit_roundoff)
        * np.max(np.sum(absolute, axis=1))
        + hermitian_defect
    )
    return block, product, {
        "maximum_contribution_count": maximum_count,
        "hermitian_defect_row_norm": ftext(hermitian_defect),
        "sparse_product_roundoff_bound": ftext(roundoff_bound),
    }


def lifted_residuals(matrix, isometry, product, eigenvalues, eigenvectors,
                     bordered=False, tangent=None):
    if bordered:
        top_vectors = eigenvectors[:-1]
        border_values = eigenvectors[-1]
        lifted_top = isometry @ top_vectors
        left_top = product @ top_vectors + tangent[:, None] * border_values
        left_border = tangent @ lifted_top
        residual_top = left_top - lifted_top * eigenvalues
        residual_border = left_border - border_values * eigenvalues
        return np.sqrt(
            np.sum(np.abs(residual_top) ** 2, axis=0)
            + np.abs(residual_border) ** 2
        )
    lifted = isometry @ eigenvectors
    left = product @ eigenvectors
    return np.linalg.norm(left - lifted * eigenvalues, axis=0)


def primary_spectrum(primary):
    values = []
    for record in primary["block_census"]["records"]:
        weight = int(record["spectral_weight"])
        for value in record["eigenvalues"]:
            values.extend([float(value)] * weight)
    return np.sort(np.asarray(values, dtype=np.float64))


print("=" * 78)
print("GEOMETRIC ADVERSARIAL COXETER REPLICATION FOR SCHEDULE 0")
print("=" * 78)

actual_hashes = {name: digest(ROOT / name) for name in EXPECTED_HASHES}
provenance_ok = check(
    "the adversarial protocol and all frozen inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES and PROTOCOL_COMMIT == "b54f590",
    str(actual_hashes),
)

source = load_source_definitions()
actions = source["load_action_definitions"]()
definitions_ok = check(
    "the primary block implementation is not loaded",
    "OUTPUT" not in source and "OUTPUT" not in actions,
)
frozen = json.loads(FROZEN.read_text())
curvature = json.loads(CURVATURE.read_text())

vertices, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
vertex_cells, top, colours = actions["barycentric_chambers"](coarse_top)
schedule = source["schedule_geometry"](actions, top, colours, WORD)
simplex_kinds, triangle_kinds = source["global_pattern_catalogue"]((schedule,))
geometries = {
    dps: actions["exact_geometry"](dps) for dps in DECIMAL_PRECISIONS
}
for geometry in geometries.values():
    geometry["mass"] = mp.mpf(0)
masses = tuple(
    mp.mpf(value)
    for value in curvature["selected_rank_matter"]["per_vertex_masses"]
)

print("Rebuilding the frozen schedule-0 matrix...", flush=True)
angle_cache, stencil_diagnostics = source["build_angle_cache"](
    actions, simplex_kinds, geometries
)
area_cache = source["build_area_cache"](triangle_kinds, geometries[140])
matrix, matrix_diagnostics = source["assemble_internal"](
    schedule, angle_cache, area_cache, masses, geometries[140]
)
matrix_digest = source["csr_digest"](matrix)
expected_digest = frozen["census"]["schedules"][0]["csr_sha256"]
matrix_ok = check(
    "the adversarial source matrix reproduces the frozen CSR digest",
    matrix_digest == expected_digest,
    f"actual={matrix_digest}, expected={expected_digest}",
)

centres = cell_centres(vertices, vertex_cells)
reflections, normals = chamber_reflections(centres, top[0])
pair_orders = tuple(tuple(
    1 if left == right else matrix_order(
        reflections[left] @ reflections[right], limit=10
    )
    for right in range(4)
) for left in range(4))
expected_orders = (
    (1, 3, 2, 2),
    (3, 1, 3, 2),
    (2, 3, 1, 5),
    (2, 2, 5, 1),
)
reflection_residual = max(
    float(np.linalg.norm(reflection @ reflection - np.eye(4), ord=np.inf))
    for reflection in reflections
)
single_reflection_order = matrix_order(reflections[0])

coxeter = sequential_product(reflections, WORD)
reverse_coxeter = reflections[0] @ reflections[1] @ reflections[2] @ reflections[3]
vertex_map, vertex_residual, vertex_bijective = vertex_permutation(
    vertices, coxeter
)
reverse_vertex_map, reverse_vertex_residual, reverse_vertex_bijective = (
    vertex_permutation(vertices, reverse_coxeter)
)
cell_map, cell_bijective, missing_cells = induced_cell_permutation(
    vertex_cells, vertex_map
)
reverse_cell_map, reverse_cell_bijective, reverse_missing_cells = (
    induced_cell_permutation(vertex_cells, reverse_vertex_map)
)
edge_map, edge_invariant = induced_edge_permutation(
    source, schedule, cell_map, len(vertex_cells)
)
reverse_edge_map, reverse_edge_invariant = induced_edge_permutation(
    source, schedule, reverse_cell_map, len(vertex_cells)
)

reflection_ok = check(
    "ambient R4 reflections independently reproduce H4 and the order-30 action",
    pair_orders == expected_orders
    and reflection_residual <= MATCH_GATE
    and single_reflection_order == 2
    and matrix_order(coxeter) == ORDER
    and matrix_order(reverse_coxeter) == ORDER
    and vertex_bijective
    and reverse_vertex_bijective
    and vertex_residual <= MATCH_GATE
    and reverse_vertex_residual <= MATCH_GATE
    and cell_bijective
    and reverse_cell_bijective
    and missing_cells == 0
    and reverse_missing_cells == 0
    and permutation_order(cell_map) == ORDER
    and permutation_order(reverse_cell_map) == ORDER
    and np.array_equal(colours[cell_map], colours)
    and np.array_equal(colours[reverse_cell_map], colours)
    and edge_invariant
    and reverse_edge_invariant
    and permutation_order(edge_map) == ORDER
    and permutation_order(reverse_edge_map) == ORDER,
    f"pair orders={pair_orders}, vertex residual={vertex_residual:.3e}",
)

inverse_relation = np.array_equal(
    reverse_edge_map[edge_map], np.arange(len(edge_map), dtype=np.int32)
)
convention_ok = check(
    "the reverse matrix word induces the inverse internal-edge action",
    inverse_relation,
)
negative_order_ok = check(
    "the single-reflection negative control has order two rather than thirty",
    single_reflection_order == 2 and single_reflection_order != ORDER,
)

cycles, cycle_id, cycle_position = cycle_coordinates(edge_map)
reverse_cycles, reverse_cycle_id, reverse_cycle_position = cycle_coordinates(
    reverse_edge_map
)
cycle_histogram = dict(sorted(Counter(map(len, cycles)).items()))
reverse_cycle_histogram = dict(sorted(Counter(map(len, reverse_cycles)).items()))
cycles_ok = check(
    "both geometric Coxeter conventions give exactly 656 free length-30 cycles",
    cycle_histogram == {30: DIMENSION}
    and reverse_cycle_histogram == {30: DIMENSION},
    f"forward={cycle_histogram}, reverse={reverse_cycle_histogram}",
)

operator_error = matrix_diagnostics["operator_error_row_bound"]
covariance_values, averaging_bound = covariance_defects(
    source, matrix, edge_map
)
reverse_covariance_values, reverse_averaging_bound = covariance_defects(
    source, matrix, reverse_edge_map
)
covariance_ok = check(
    "both geometric actions preserve the Hessian inside the local operator envelope",
    averaging_bound <= 100 * operator_error
    and reverse_averaging_bound <= 100 * operator_error,
    f"forward={averaging_bound:.3e}, reverse={reverse_averaging_bound:.3e}",
)

corrupted = matrix.copy().tolil()
corrupted[0, 0] += 1e-4
corrupted = corrupted.tocsr()
_, corrupted_bound = covariance_defects(source, corrupted, edge_map)
corruption_ok = check(
    "the preregistered diagonal corruption violates geometric covariance",
    corrupted_bound > 100 * operator_error,
    f"corrupted bound={corrupted_bound:.6e}",
)

# Character-independent support and absolute contribution matrices for the
# explicit sparse products.  These are not direct phase accumulations.
support_counts, absolute_bound = sparse_support_products(matrix, cycle_id)
reverse_support_counts, reverse_absolute_bound = sparse_support_products(
    matrix, reverse_cycle_id
)

tangent = source["product_tangent"](schedule, geometries[140])
block_records = []
blocks = {}
eigenvalue_sets = {}
all_diagonalized = True
all_isometric = True
all_separated = True
maximum_roundoff = 0.0
maximum_lifted_residual = 0.0
zero_candidates = []

for character in range(16):
    print(
        f"Explicit sparse-Q block and full evd spectrum k={character}/15...",
        flush=True,
    )
    isometry = explicit_fourier_isometry(
        cycle_id, cycle_position, character
    )
    gram = (isometry.conj().T @ isometry).toarray()
    isometry_defect = float(np.linalg.norm(
        gram - np.eye(DIMENSION), ord=np.inf
    ))
    all_isometric &= isometry_defect <= 1e-12
    block, product, diagnostics = sparse_product_block(
        matrix, isometry, support_counts, absolute_bound
    )
    bordered = character == 0
    if bordered:
        projection = np.asarray(isometry.conj().T @ tangent).ravel()
        complete_block = np.zeros((657, 657), dtype=np.complex128)
        complete_block[:656, :656] = block
        complete_block[:656, 656] = projection
        complete_block[656, :656] = projection.conj()
        block = complete_block

    try:
        values, vectors = la.eigh(
            block, driver="evd", check_finite=True, overwrite_a=False
        )
        residuals = lifted_residuals(
            matrix,
            isometry,
            product,
            values,
            vectors,
            bordered=bordered,
            tangent=tangent if bordered else None,
        )
    except Exception as error:
        all_diagonalized = False
        block_records.append({
            "character": character,
            "error": repr(error),
            "isometry_defect": ftext(isometry_defect),
        })
        continue

    roundoff = float(diagnostics["sparse_product_roundoff_bound"])
    maximum_roundoff = max(maximum_roundoff, roundoff)
    maximum_lifted_residual = max(
        maximum_lifted_residual, float(residuals.max())
    )
    gates = 100 * (
        operator_error + averaging_bound + roundoff + residuals
    )
    separated = np.abs(values) > gates
    all_separated &= bool(np.all(separated))
    nearest = int(np.argmin(np.abs(values)))
    blocks[character] = block
    eigenvalue_sets[character] = values
    block_records.append({
        "character": character,
        "dimension": len(values),
        "spectral_weight": 1 if character in (0, 15) else 2,
        "isometry_defect": ftext(isometry_defect),
        "minimum_absolute_eigenvalue": ftext(np.min(np.abs(values))),
        "nearest_eigenvalue": ftext(values[nearest]),
        "nearest_lifted_residual": ftext(residuals[nearest]),
        "nearest_gate": ftext(gates[nearest]),
        "all_eigenvalues_separated": bool(np.all(separated)),
        "maximum_lifted_residual": ftext(residuals.max()),
        "eigenvalues": [ftext(value) for value in values],
        "diagnostics": diagnostics,
    })
    for index in np.flatnonzero(~separated):
        zero_candidates.append({
            "character": character,
            "eigenvalue": ftext(values[index]),
            "lifted_residual": ftext(residuals[index]),
            "gate": ftext(gates[index]),
        })

blocks_ok = check(
    "all explicit Fourier matrices are isometries and all blocks are diagonalized",
    all_isometric
    and all_diagonalized
    and len(block_records) == 16
    and sum(
        record.get("dimension", 0) * record.get("spectral_weight", 0)
        for record in block_records
    ) == 19681,
)

# Reverse-action and wrong-phase convention attacks are evaluated after the
# decisive forward spectrum is complete.
reverse_q1 = explicit_fourier_isometry(
    reverse_cycle_id, reverse_cycle_position, 1
)
reverse_block1, _, reverse_diagnostics = sparse_product_block(
    matrix, reverse_q1, reverse_support_counts, reverse_absolute_bound
)
reverse_values1 = la.eigvalsh(reverse_block1, driver="evd")
reverse_sector_difference = float(np.max(np.abs(
    reverse_values1 - eigenvalue_sets.get(1, np.full(DIMENSION, np.inf))
)))
reverse_sector_gate = 100 * (
    operator_error + averaging_bound + reverse_averaging_bound
    + maximum_roundoff
    + float(reverse_diagnostics["sparse_product_roundoff_bound"])
)
reverse_ok = check(
    "the inverse geometric action reproduces the nontrivial sector spectrum",
    reverse_sector_difference <= reverse_sector_gate,
    f"difference={reverse_sector_difference:.3e}, gate={reverse_sector_gate:.3e}",
)

wrong_q1 = explicit_fourier_isometry(
    cycle_id, cycle_position, 1, sign=1
)
wrong_block1, _, wrong_diagnostics = sparse_product_block(
    matrix, wrong_q1, support_counts, absolute_bound
)
wrong_same_sector_difference = float(np.linalg.norm(
    wrong_block1 - blocks.get(1, np.zeros_like(wrong_block1)), ord=np.inf
))
wrong_phase_gate = 100 * (
    operator_error + averaging_bound + maximum_roundoff
    + float(wrong_diagnostics["sparse_product_roundoff_bound"])
)
wrong_phase_ok = check(
    "the unrelabeled opposite Fourier phase fails the same-sector block test",
    wrong_same_sector_difference > wrong_phase_gate,
    f"difference={wrong_same_sector_difference:.3e}, gate={wrong_phase_gate:.3e}",
)

trace_full = float(matrix.diagonal().sum())
trace_blocks = sum(
    record.get("spectral_weight", 0)
    * float(np.sum(eigenvalue_sets.get(record["character"], ())))
    for record in block_records
)
trace_difference = abs(trace_full - trace_blocks)
trace_gate = 100 * (
    19681 * maximum_lifted_residual
    + np.finfo(float).eps * max(abs(trace_full), 1.0)
)
trace_ok = check(
    "the adversarial weighted traces reproduce the full bordered trace",
    trace_difference <= trace_gate,
    f"difference={trace_difference:.3e}, gate={trace_gate:.3e}",
)

frobenius_full = np.sqrt(
    float(np.sum(matrix.data ** 2)) + 2 * float(np.dot(tangent, tangent))
)
frobenius_blocks = np.sqrt(sum(
    record.get("spectral_weight", 0)
    * float(np.linalg.norm(blocks[record["character"]], ord="fro") ** 2)
    for record in block_records if record["character"] in blocks
))
frobenius_difference = abs(frobenius_full - frobenius_blocks)
frobenius_gate = 100 * (
    np.sqrt(19681) * (
        averaging_bound + maximum_roundoff + maximum_lifted_residual
    )
    + np.finfo(float).eps * max(frobenius_full, 1.0)
)
parseval_ok = check(
    "the adversarial weighted Frobenius norm satisfies Parseval",
    frobenius_difference <= frobenius_gate,
    f"difference={frobenius_difference:.3e}, gate={frobenius_gate:.3e}",
)

adversarial_values = []
for record in block_records:
    for value in eigenvalue_sets.get(record["character"], ()):
        adversarial_values.extend(
            [float(value)] * record.get("spectral_weight", 0)
        )
adversarial_values = np.sort(np.asarray(adversarial_values, dtype=np.float64))

# The primary artifact is first parsed only after the adversarial full
# spectrum has been constructed.
primary = json.loads(PRIMARY.read_text())
primary_values = primary_spectrum(primary)
primary_common = float(
    primary["block_census"]["common_operator_and_averaging_error"]
)
primary_roundoff = float(primary["block_census"]["maximum_block_roundoff"])
primary_residual = float(primary["block_census"]["maximum_eigenpair_residual"])
primary_forward_bound = 100 * (
    primary_common + primary_roundoff + primary_residual
)
adversarial_forward_bound = 100 * (
    operator_error + averaging_bound + maximum_roundoff
    + maximum_lifted_residual
)
comparison_gate = primary_forward_bound + adversarial_forward_bound
if len(adversarial_values) == len(primary_values) == 19681:
    spectrum_difference = float(np.max(np.abs(
        adversarial_values - primary_values
    )))
else:
    spectrum_difference = np.inf
full_spectrum_ok = check(
    "the complete adversarial and primary spectra agree within both envelopes",
    spectrum_difference <= comparison_gate,
    f"difference={spectrum_difference:.3e}, gate={comparison_gate:.3e}",
)

construction_ok = all((
    provenance_ok,
    definitions_ok,
    matrix_ok,
    reflection_ok,
    convention_ok,
    negative_order_ok,
    cycles_ok,
    covariance_ok,
    corruption_ok,
    blocks_ok,
    reverse_ok,
    wrong_phase_ok,
    trace_ok,
    parseval_ok,
    full_spectrum_ok,
))
if not construction_ok:
    outcome = "ADVERSARIAL_CONSTRUCTION_INVALID"
elif all_separated:
    outcome = "SCHEDULE0_KERNEL_ADVERSARIALLY_CORROBORATED"
else:
    outcome = "PRIMARY_SCHEDULE0_KERNEL_NOT_CORROBORATED"

check(
    "the adversarial census follows the frozen verdict hierarchy",
    outcome in {
        "ADVERSARIAL_CONSTRUCTION_INVALID",
        "SCHEDULE0_KERNEL_ADVERSARIALLY_CORROBORATED",
        "PRIMARY_SCHEDULE0_KERNEL_NOT_CORROBORATED",
    }
    and (outcome == "ADVERSARIAL_CONSTRUCTION_INVALID") == (not construction_ok)
    and (
        outcome == "SCHEDULE0_KERNEL_ADVERSARIALLY_CORROBORATED"
    ) == (construction_ok and all_separated),
    outcome,
)

artifact = {
    "title": "Geometric adversarial replication of schedule-0 Coxeter blocks",
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "scope": {
        "schedule_indices": [0],
        "all_cyclic_sectors_included": True,
        "primary_block_code_loaded": False,
        "full_suite_run": False,
        "physical_target_loaded": False,
    },
    "matrix": {
        "csr_sha256": matrix_digest,
        "dimension": matrix.shape[0],
        "nnz": matrix.nnz,
        "operator_error_row_bound": ftext(operator_error),
        "stencil_diagnostics": stencil_diagnostics,
    },
    "geometric_action": {
        "reflection_normals": [
            [ftext(value) for value in normal] for normal in normals
        ],
        "pair_orders": [list(row) for row in pair_orders],
        "single_reflection_order": single_reflection_order,
        "coxeter_matrix_order": matrix_order(coxeter),
        "reverse_coxeter_matrix_order": matrix_order(reverse_coxeter),
        "vertex_matching_residual": ftext(vertex_residual),
        "reverse_vertex_matching_residual": ftext(reverse_vertex_residual),
        "cell_order": permutation_order(cell_map),
        "edge_order": permutation_order(edge_map),
        "cycle_histogram": {str(key): value for key, value in cycle_histogram.items()},
        "reverse_cycle_histogram": {
            str(key): value for key, value in reverse_cycle_histogram.items()
        },
        "reverse_is_inverse": bool(inverse_relation),
    },
    "covariance": {
        "forward_power_defects": [ftext(value) for value in covariance_values],
        "reverse_power_defects": [
            ftext(value) for value in reverse_covariance_values
        ],
        "forward_group_average_bound": ftext(averaging_bound),
        "reverse_group_average_bound": ftext(reverse_averaging_bound),
        "corrupted_group_average_bound": ftext(corrupted_bound),
    },
    "block_census": {
        "independent_block_count": len(block_records),
        "weighted_dimension": len(adversarial_values),
        "all_eigenvalues_separated": bool(all_separated),
        "maximum_sparse_product_roundoff": ftext(maximum_roundoff),
        "maximum_lifted_residual": ftext(maximum_lifted_residual),
        "records": block_records,
        "zero_candidates": zero_candidates,
        "reverse_sector_difference": ftext(reverse_sector_difference),
        "reverse_sector_gate": ftext(reverse_sector_gate),
        "wrong_phase_same_sector_difference": ftext(wrong_same_sector_difference),
        "wrong_phase_gate": ftext(wrong_phase_gate),
        "trace_full": ftext(trace_full),
        "trace_blocks": ftext(trace_blocks),
        "frobenius_full": ftext(frobenius_full),
        "frobenius_blocks": ftext(frobenius_blocks),
        "primary_full_spectrum_maximum_difference": ftext(spectrum_difference),
        "primary_full_spectrum_comparison_gate": ftext(comparison_gate),
        "smallest_32_by_absolute_value": [
            ftext(value) for value in adversarial_values[
                np.argsort(np.abs(adversarial_values))[:32]
            ]
        ] if len(adversarial_values) else [],
    },
    "status_labels": {
        "schedule0_kernel": (
            "DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED"
            if outcome == "SCHEDULE0_KERNEL_ADVERSARIALLY_CORROBORATED"
            else "OPEN"
        ),
        "other_schedules": "OPEN / NOT TESTED",
        "graviton_or_propagation": "OPEN / NOT TESTED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Geometric matching residual: {vertex_residual:.6e}")
print(f"All adversarial eigenvalues separated: {all_separated}")
print(f"Full-spectrum primary difference: {spectrum_difference:.6e}")
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests} passed")
print(f"Artifact: {OUTPUT.relative_to(ROOT)}")
sys.exit(0 if passed == tests else 1)

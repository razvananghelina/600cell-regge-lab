#!/usr/bin/env python3
"""Exhaustive C30 block certificate for refined internal schedule 0.

Prior-art commit: 295e90d.
Protocol commit: acfe795.
Every cyclic sector is diagonalized; no physical or continuum target is used.
"""

import ast
from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
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
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_block_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_block_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_nonhomogeneous_coxeter_blocks.json"

PRIOR_ART_COMMIT = "295e90d"
PROTOCOL_COMMIT = "acfe795"
EXPECTED_HASHES = {
    "reproducible/verify_gravity_600cell_refined_nonhomogeneous_internal_hessian.py":
        "2a96c8ce466d6d4e9be2cadc4ebf932b4e42eff16fc64fbfb08cd580d680879e",
    "reproducible/gravity_600cell_refined_nonhomogeneous_internal_hessian.json":
        "4a05968c68f8e6a35a1308ddf6114bb19b7106f214bfdcf798e7af2387bddec1",
    "reproducible/gravity_600cell_refined_local_curvature_mass.json":
        "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_block_prior_art.md":
        "e5db7f7fe28d778a11a40b91115b00d9d4beffe7207ce2d0d5a3bb6e3b2f2018",
    "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_block_protocol.md":
        "808145af3bd49d13f1140fff1ad177ccdb6bf029128e7a253ccf1a17757ddfc6",
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
COXETER_WORD = (0, 1, 2, 3)
COXETER_ORDER = 30
EXPECTED_CYCLES = 656
EXPECTED_SECTOR_DIMENSION = 656

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


def compose(left, right):
    """Permutation left after right, arrays mapping old index to new."""
    return left[right]


def permutation_order(mapping, limit=120):
    identity = np.arange(len(mapping), dtype=np.int32)
    power = identity.copy()
    for exponent in range(1, limit + 1):
        power = compose(mapping, power)
        if np.array_equal(power, identity):
            return exponent
    return None


def coloured_chamber_maps(top):
    maps = []
    for colour in range(4):
        residues = defaultdict(list)
        for chamber, row in enumerate(top):
            key = tuple(int(row[index]) for index in range(4) if index != colour)
            residues[key].append(chamber)
        if any(len(pair) != 2 for pair in residues.values()):
            raise RuntimeError("a coloured chamber residue does not have size two")
        mapping = np.empty(len(top), dtype=np.int32)
        for left, right in residues.values():
            mapping[left] = right
            mapping[right] = left
        maps.append(mapping)
    return tuple(maps)


def left_action_from_base(colour_maps, word):
    target = 0
    for colour in word:
        target = int(colour_maps[colour][target])
    image = np.full(len(colour_maps[0]), -1, dtype=np.int32)
    image[0] = target
    queue = deque((0,))
    consistent = True
    while queue:
        chamber = queue.popleft()
        image_chamber = int(image[chamber])
        for colour in range(4):
            neighbour = int(colour_maps[colour][chamber])
            image_neighbour = int(colour_maps[colour][image_chamber])
            if image[neighbour] < 0:
                image[neighbour] = image_neighbour
                queue.append(neighbour)
            elif image[neighbour] != image_neighbour:
                consistent = False
    return image, target, consistent


def induced_cell_map(top, chamber_map, cell_count):
    mapping = np.full(cell_count, -1, dtype=np.int32)
    contradictions = 0
    for chamber, image_chamber in enumerate(chamber_map):
        for rank in range(4):
            old = int(top[chamber, rank])
            new = int(top[image_chamber, rank])
            if mapping[old] < 0:
                mapping[old] = new
            elif mapping[old] != new:
                contradictions += 1
    well_defined = (
        contradictions == 0
        and np.all(mapping >= 0)
        and len(set(map(int, mapping))) == cell_count
    )
    return mapping, well_defined, contradictions


def induced_internal_edge_map(source, schedule, cell_map, cell_count):
    mapped = schedule["internal_edges"].copy()
    for endpoint in range(2):
        layer = mapped[:, endpoint] >= cell_count
        base = mapped[:, endpoint] % cell_count
        mapped[:, endpoint] = cell_map[base] + layer * cell_count
    mapped.sort(axis=1)
    edge_codes = source["encoded_edges"](
        schedule["internal_edges"], 2 * cell_count
    )
    mapped_codes = source["encoded_edges"](mapped, 2 * cell_count)
    positions = np.searchsorted(edge_codes, mapped_codes)
    invariant = (
        np.all(positions < len(edge_codes))
        and np.array_equal(edge_codes[positions], mapped_codes)
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
            raise RuntimeError("permutation traversal entered an earlier cycle")
        cycles.append(np.asarray(cycle, dtype=np.int32))
    return tuple(cycles), cycle_id, position


def covariance_defects(source, matrix, mapping):
    identity = np.arange(len(mapping), dtype=np.int32)
    power = identity.copy()
    defects = []
    for exponent in range(COXETER_ORDER):
        if exponent:
            power = compose(mapping, power)
        pulled = matrix[power, :][:, power]
        defects.append(source["sparse_row_norm"](matrix - pulled))
    return defects, float(np.mean(defects))


def fourier_block(matrix_coo, cycle_id, position, character, dimension):
    rows = matrix_coo.row
    columns = matrix_coo.col
    codes = (
        cycle_id[rows].astype(np.int64) * dimension
        + cycle_id[columns].astype(np.int64)
    )
    phase = np.exp(
        2j * np.pi * character
        * (position[rows].astype(np.float64)
           - position[columns].astype(np.float64)) / COXETER_ORDER
    )
    values = matrix_coo.data * phase / COXETER_ORDER
    block = np.zeros(dimension * dimension, dtype=np.complex128)
    absolute = np.zeros(dimension * dimension, dtype=np.float64)
    counts = np.zeros(dimension * dimension, dtype=np.int32)
    np.add.at(block, codes, values)
    np.add.at(absolute, codes, np.abs(values))
    np.add.at(counts, codes, 1)
    block = block.reshape(dimension, dimension)
    absolute = absolute.reshape(dimension, dimension)
    counts = counts.reshape(dimension, dimension)

    unit_roundoff = np.finfo(np.float64).eps / 2
    gamma = np.zeros_like(absolute)
    nonzero = counts > 0
    operations = counts[nonzero] + 8
    gamma[nonzero] = operations * unit_roundoff / (
        1 - operations * unit_roundoff
    )
    entry_error = (gamma + 4 * unit_roundoff) * absolute
    hermitian_defect = float(np.max(np.sum(np.abs(block - block.conj().T), axis=1)))
    block = (block + block.conj().T) / 2
    entry_error = (entry_error + entry_error.T) / 2
    roundoff_bound = float(np.max(np.sum(entry_error, axis=1)))
    return block, entry_error, {
        "maximum_term_count": int(counts.max()),
        "hermitian_defect_row_norm": ftext(hermitian_defect),
        "roundoff_row_bound": ftext(roundoff_bound),
    }


def tangent_projection(tangent, cycles, character):
    result = np.empty(len(cycles), dtype=np.complex128)
    root = np.sqrt(COXETER_ORDER)
    phases = np.exp(
        2j * np.pi * character * np.arange(COXETER_ORDER) / COXETER_ORDER
    )
    for index, cycle in enumerate(cycles):
        result[index] = np.dot(phases, tangent[cycle]) / root
    return result


def diagonalize_block(block, common_error, block_roundoff):
    eigenvalues, eigenvectors = np.linalg.eigh(block)
    residual_matrix = block @ eigenvectors - eigenvectors * eigenvalues
    residuals = np.linalg.norm(residual_matrix, axis=0)
    gates = 100 * (common_error + block_roundoff + residuals)
    separated = np.abs(eigenvalues) > gates
    nearest = int(np.argmin(np.abs(eigenvalues)))
    return eigenvalues, eigenvectors, residuals, separated, nearest


def reconstruct_sector_vector(coefficients, cycles, character, bordered=False):
    total = sum(map(len, cycles))
    result = np.zeros(total + int(bordered), dtype=np.complex128)
    phases = np.exp(
        -2j * np.pi * character * np.arange(COXETER_ORDER) / COXETER_ORDER
    ) / np.sqrt(COXETER_ORDER)
    for cycle_index, cycle in enumerate(cycles):
        result[cycle] = coefficients[cycle_index] * phases
    if bordered:
        result[-1] = coefficients[-1]
    return result


print("=" * 78)
print("EXHAUSTIVE COXETER BLOCKS FOR REFINED INTERNAL SCHEDULE 0")
print("=" * 78)

actual_hashes = {name: digest(ROOT / name) for name in EXPECTED_HASHES}
provenance_ok = check(
    "the prior-art gate, protocol and frozen inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIOR_ART_COMMIT == "295e90d"
    and PROTOCOL_COMMIT == "acfe795",
    str(actual_hashes),
)

source = load_source_definitions()
actions = source["load_action_definitions"]()
definitions_ok = check(
    "only frozen function definitions are loaded",
    "OUTPUT" not in source and "OUTPUT" not in actions,
)
frozen = json.loads(FROZEN.read_text())
curvature = json.loads(CURVATURE.read_text())

_, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
vertex_cells, top, colours = actions["barycentric_chambers"](coarse_top)
schedule = source["schedule_geometry"](
    actions, top, colours, COXETER_WORD
)
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

print("Building only the four schedule-0 angle stencils...", flush=True)
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
    "the complete schedule-0 CSR matrix reproduces the frozen digest",
    matrix_digest == expected_digest,
    f"actual={matrix_digest}, expected={expected_digest}",
)

colour_maps = coloured_chamber_maps(top)
coxeter_matrix = tuple(tuple(
    1 if left == right else permutation_order(
        compose(colour_maps[left], colour_maps[right]), 10
    )
    for right in range(4)
) for left in range(4))
expected_coxeter = (
    (1, 3, 2, 2),
    (3, 1, 3, 2),
    (2, 3, 1, 5),
    (2, 2, 5, 1),
)
left_action, base_target, left_consistent = left_action_from_base(
    colour_maps, COXETER_WORD
)
cell_map, cells_well_defined, cell_contradictions = induced_cell_map(
    top, left_action, len(vertex_cells)
)
edge_map, edges_invariant = induced_internal_edge_map(
    source, schedule, cell_map, len(vertex_cells)
)
action_ok = check(
    "the colour-preserving left Coxeter action is well defined and has order 30",
    coxeter_matrix == expected_coxeter
    and left_consistent
    and len(set(map(int, left_action))) == len(left_action)
    and permutation_order(left_action) == COXETER_ORDER
    and cells_well_defined
    and cell_contradictions == 0
    and permutation_order(cell_map) == COXETER_ORDER
    and np.array_equal(colours[cell_map], colours)
    and edges_invariant
    and permutation_order(edge_map) == COXETER_ORDER,
    f"base target={base_target}, relations={coxeter_matrix}",
)

# Convention-negative control: the raw right product is a chamber
# permutation but does not descend to fixed-rank cell residues.
right_action = np.arange(len(top), dtype=np.int32)
for colour in COXETER_WORD:
    right_action = compose(colour_maps[colour], right_action)
_, right_well_defined, right_contradictions = induced_cell_map(
    top, right_action, len(vertex_cells)
)
right_control_ok = check(
    "the forbidden right-product convention fails cell well-definedness",
    not right_well_defined and right_contradictions > 0,
    f"contradictions={right_contradictions}",
)

cycles, cycle_id, cycle_position = cycle_coordinates(edge_map)
cycle_histogram = dict(sorted(Counter(map(len, cycles)).items()))
sector_dimensions = [
    sum(1 for cycle in cycles if (character * len(cycle)) % COXETER_ORDER == 0)
    for character in range(COXETER_ORDER)
]
cycles_ok = check(
    "the edge action has exactly 656 free C30 cycles and 30 complete sectors",
    cycle_histogram == {30: EXPECTED_CYCLES}
    and sector_dimensions == [EXPECTED_SECTOR_DIMENSION] * COXETER_ORDER
    and sum(sector_dimensions) == matrix.shape[0],
    f"cycles={cycle_histogram}, dimensions={sector_dimensions}",
)

operator_error = matrix_diagnostics["operator_error_row_bound"]
covariance_values, averaging_bound = covariance_defects(
    source, matrix, edge_map
)
covariance_ok = check(
    "the group-averaging distance is inside the declared local operator envelope",
    averaging_bound <= 100 * operator_error,
    f"average={averaging_bound:.6e}, gate={100*operator_error:.6e}",
)

corrupted = matrix.copy().tolil()
corrupted[0, 0] += 1e-4
corrupted = corrupted.tocsr()
_, corrupted_bound = covariance_defects(source, corrupted, edge_map)
corruption_ok = check(
    "a single preregistered diagonal corruption violates Coxeter covariance",
    corrupted_bound > 100 * operator_error,
    f"corrupted average={corrupted_bound:.6e}",
)

tangent = source["product_tangent"](schedule, geometries[140])
tangent_projections = [
    tangent_projection(tangent, cycles, character)
    for character in range(COXETER_ORDER)
]
unit_roundoff = np.finfo(np.float64).eps / 2
gamma30 = 30 * unit_roundoff / (1 - 30 * unit_roundoff)
tangent_gate = 100 * gamma30 * max(
    float(np.max([np.sum(np.abs(tangent[cycle])) for cycle in cycles])),
    np.finfo(float).tiny,
)
tangent_noninvariant = max(
    float(np.max(np.abs(tangent_projections[character])))
    for character in range(1, COXETER_ORDER)
)
tangent_ok = check(
    "the analytic duration tangent lies only in the invariant Coxeter sector",
    tangent_noninvariant <= tangent_gate,
    f"maximum noninvariant={tangent_noninvariant:.6e}, gate={tangent_gate:.6e}",
)

matrix_coo = matrix.tocoo()
matrix_coo.sum_duplicates()
common_error = operator_error + averaging_bound
block_records = []
blocks = {}
eigenvalue_sets = {}
eigenvector_sets = {}
all_diagonalized = True
all_separated = True
maximum_block_roundoff = 0.0
maximum_eigen_residual = 0.0
zero_candidates = []

for character in range(16):
    print(f"Building and diagonalizing exhaustive sector k={character}/15...", flush=True)
    block, entry_error, diagnostics = fourier_block(
        matrix_coo, cycle_id, cycle_position,
        character, EXPECTED_SECTOR_DIMENSION,
    )
    if character == 0:
        projection = tangent_projections[0]
        bordered = np.zeros((657, 657), dtype=np.complex128)
        bordered[:656, :656] = block
        bordered[:656, 656] = projection
        bordered[656, :656] = projection.conj()
        block = bordered

        tangent_error = (
            gamma30 + 4 * unit_roundoff
        ) * np.asarray([
            np.sum(np.abs(tangent[cycle])) / np.sqrt(COXETER_ORDER)
            for cycle in cycles
        ])
        bordered_error = np.zeros((657, 657), dtype=np.float64)
        bordered_error[:656, :656] = entry_error
        bordered_error[:656, 656] = tangent_error
        bordered_error[656, :656] = tangent_error
        entry_error = bordered_error
        diagnostics["roundoff_row_bound"] = ftext(
            float(np.max(np.sum(entry_error, axis=1)))
        )

    block_roundoff = float(diagnostics["roundoff_row_bound"])
    maximum_block_roundoff = max(maximum_block_roundoff, block_roundoff)
    try:
        values, vectors, residuals, separated, nearest = diagonalize_block(
            block, common_error, block_roundoff
        )
    except Exception as error:
        all_diagonalized = False
        block_records.append({
            "character": character,
            "dimension": len(block),
            "error": repr(error),
        })
        continue
    maximum_eigen_residual = max(maximum_eigen_residual, float(residuals.max()))
    all_separated &= bool(np.all(separated))
    blocks[character] = block
    eigenvalue_sets[character] = values
    eigenvector_sets[character] = vectors
    record = {
        "character": character,
        "dimension": len(block),
        "spectral_weight": 1 if character in (0, 15) else 2,
        "minimum_absolute_eigenvalue": ftext(np.min(np.abs(values))),
        "nearest_index": nearest,
        "nearest_eigenvalue": ftext(values[nearest]),
        "nearest_residual": ftext(residuals[nearest]),
        "nearest_gate": ftext(100 * (
            common_error + block_roundoff + residuals[nearest]
        )),
        "all_eigenvalues_separated": bool(np.all(separated)),
        "maximum_eigenpair_residual": ftext(residuals.max()),
        "eigenvalues": [ftext(value) for value in values],
        "diagnostics": diagnostics,
    }
    block_records.append(record)
    if not np.all(separated):
        for index in np.flatnonzero(~separated):
            full = reconstruct_sector_vector(
                vectors[:, index], cycles, character, bordered=(character == 0)
            )
            if character == 0:
                column = sp.csr_matrix(tangent.reshape(-1, 1))
                full_matrix = sp.bmat(
                    [[matrix, column], [column.T, sp.csr_matrix((1, 1))]],
                    format="csr",
                )
            else:
                full_matrix = matrix
            residual = np.linalg.norm(
                full_matrix @ full - values[index] * full
            )
            zero_candidates.append({
                "character": character,
                "eigenvalue": ftext(values[index]),
                "block_residual": ftext(residuals[index]),
                "full_residual_against_unaveraged_matrix": ftext(residual),
            })

diagonalization_ok = check(
    "all 16 independent blocks are diagonalized exhaustively",
    all_diagonalized
    and len(block_records) == 16
    and sum(
        record["dimension"] * record.get("spectral_weight", 0)
        for record in block_records
    ) == 19681,
)

# Explicitly build the conjugate partner k=29 rather than infer it.
block29, error29, diagnostics29 = fourier_block(
    matrix_coo, cycle_id, cycle_position, 29, EXPECTED_SECTOR_DIMENSION
)
if 1 in blocks:
    record1 = next(
        record for record in block_records if record.get("character") == 1
    )
    conjugate_difference = float(np.max(np.sum(
        np.abs(block29 - blocks[1].conj()), axis=1
    )))
    conjugate_gate = 100 * (
        float(diagnostics29["roundoff_row_bound"])
        + float(record1["diagnostics"]["roundoff_row_bound"])
        + averaging_bound
    )
else:
    conjugate_difference = np.inf
    conjugate_gate = 0.0
conjugate_ok = check(
    "the explicitly constructed k=29 block is conjugate to k=1",
    conjugate_difference <= conjugate_gate,
    f"difference={conjugate_difference:.6e}, gate={conjugate_gate:.6e}",
)

trace_full = float(matrix.diagonal().sum())
trace_blocks = sum(
    record["spectral_weight"] * float(np.sum(eigenvalue_sets[record["character"]]))
    for record in block_records if record["character"] in eigenvalue_sets
)
trace_error = abs(trace_full - trace_blocks)
trace_gate = 100 * (
    19681 * np.finfo(float).eps * max(1.0, abs(trace_full))
    + 19681 * maximum_eigen_residual
)
trace_ok = check(
    "the weighted block traces reproduce the full bordered trace",
    trace_error <= trace_gate,
    f"difference={trace_error:.6e}, gate={trace_gate:.6e}",
)

frobenius_full = np.sqrt(
    float(np.sum(matrix.data ** 2)) + 2 * float(np.dot(tangent, tangent))
)
frobenius_blocks = np.sqrt(sum(
    record["spectral_weight"] * float(np.linalg.norm(
        blocks[record["character"]], ord="fro"
    ) ** 2)
    for record in block_records if record["character"] in blocks
))
frobenius_difference = abs(frobenius_full - frobenius_blocks)
frobenius_gate = 100 * (
    np.sqrt(19681) * averaging_bound
    + np.sqrt(19681) * maximum_block_roundoff
    + np.finfo(float).eps * max(frobenius_full, 1.0)
)
parseval_ok = check(
    "the weighted block Frobenius norm satisfies Parseval in the covariance envelope",
    frobenius_difference <= frobenius_gate,
    f"difference={frobenius_difference:.6e}, gate={frobenius_gate:.6e}",
)

complete_values = []
for character in range(16):
    weight = 1 if character in (0, 15) else 2
    for value in eigenvalue_sets.get(character, ()):
        complete_values.extend([float(value)] * weight)
complete_values = np.asarray(complete_values)
complete_values = complete_values[np.argsort(np.abs(complete_values))]
old_values = np.asarray([
    float(value) for value in
    frozen["census"]["time_reversal_pairs"][0]["bordered_spectrum"]
    ["eigen_runs"][0]["eigenvalues_nearest_zero"]
])
old_ritz = float(
    frozen["census"]["time_reversal_pairs"][0]["bordered_spectrum"]
    ["maximum_ritz_residual"]
)
old_gate = 100 * (
    common_error + maximum_block_roundoff + maximum_eigen_residual + old_ritz
)
old_reproduction_difference = float(np.max(
    np.abs(complete_values[:8] - old_values)
)) if len(complete_values) == 19681 else np.inf
old_control_ok = check(
    "the exhaustive blocks reproduce the eight old Ritz values post hoc",
    old_reproduction_difference <= old_gate,
    f"difference={old_reproduction_difference:.6e}, gate={old_gate:.6e}",
)

construction_ok = all((
    provenance_ok,
    definitions_ok,
    matrix_ok,
    action_ok,
    right_control_ok,
    cycles_ok,
    covariance_ok,
    corruption_ok,
    tangent_ok,
    diagonalization_ok,
    conjugate_ok,
    trace_ok,
    parseval_ok,
    old_control_ok,
))
if not construction_ok:
    outcome = "COXETER_BLOCK_CONSTRUCTION_INVALID"
elif all_separated:
    outcome = "SCHEDULE0_BORDERED_NONSINGULAR_NUMERICALLY_CERTIFIED"
else:
    outcome = "SCHEDULE0_ADDITIONAL_ZERO_COMPATIBLE_SECTOR"

check(
    "the complete block census follows the preregistered verdict hierarchy",
    outcome in {
        "COXETER_BLOCK_CONSTRUCTION_INVALID",
        "SCHEDULE0_BORDERED_NONSINGULAR_NUMERICALLY_CERTIFIED",
        "SCHEDULE0_ADDITIONAL_ZERO_COMPATIBLE_SECTOR",
    }
    and (outcome == "COXETER_BLOCK_CONSTRUCTION_INVALID") == (not construction_ok)
    and (
        outcome == "SCHEDULE0_BORDERED_NONSINGULAR_NUMERICALLY_CERTIFIED"
    ) == (construction_ok and all_separated),
    outcome,
)

artifact = {
    "title": "Exhaustive C30 blocks of the refined internal schedule-0 Hessian",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "scope": {
        "schedule_indices": [0],
        "all_cyclic_sectors_included": True,
        "full_suite_run": False,
        "old_sparse_census_rerun": False,
        "physical_target_loaded": False,
    },
    "matrix": {
        "csr_sha256": matrix_digest,
        "dimension": matrix.shape[0],
        "nnz": matrix.nnz,
        "operator_error_row_bound": ftext(operator_error),
        "local_simplex_pattern_count": len(simplex_kinds),
        "local_triangle_pattern_count": len(triangle_kinds),
        "stencil_diagnostics": stencil_diagnostics,
    },
    "coxeter_action": {
        "word": list(COXETER_WORD),
        "base_target": base_target,
        "coxeter_relations": [list(row) for row in coxeter_matrix],
        "chamber_order": permutation_order(left_action),
        "cell_order": permutation_order(cell_map),
        "internal_edge_order": permutation_order(edge_map),
        "cycle_histogram": {str(key): value for key, value in cycle_histogram.items()},
        "sector_dimensions": sector_dimensions,
        "right_convention_well_defined": bool(right_well_defined),
        "right_convention_contradictions": right_contradictions,
    },
    "covariance": {
        "power_defect_row_norms": [ftext(value) for value in covariance_values],
        "group_averaging_distance_bound": ftext(averaging_bound),
        "acceptance_gate": ftext(100 * operator_error),
        "corrupted_group_averaging_distance_bound": ftext(corrupted_bound),
    },
    "tangent": {
        "maximum_noninvariant_projection": ftext(tangent_noninvariant),
        "gate": ftext(tangent_gate),
        "invariant_projection_norm": ftext(np.linalg.norm(tangent_projections[0])),
    },
    "block_census": {
        "independent_block_count": len(block_records),
        "weighted_dimension": len(complete_values),
        "common_operator_and_averaging_error": ftext(common_error),
        "maximum_block_roundoff": ftext(maximum_block_roundoff),
        "maximum_eigenpair_residual": ftext(maximum_eigen_residual),
        "all_observed_eigenvalues_separated": bool(all_separated),
        "records": block_records,
        "zero_candidates": zero_candidates,
        "conjugate_k1_k29_difference": ftext(conjugate_difference),
        "conjugate_gate": ftext(conjugate_gate),
        "trace_full": ftext(trace_full),
        "trace_blocks": ftext(trace_blocks),
        "frobenius_full": ftext(frobenius_full),
        "frobenius_blocks": ftext(frobenius_blocks),
        "old_ritz_reproduction_difference": ftext(old_reproduction_difference),
        "old_ritz_reproduction_gate": ftext(old_gate),
        "smallest_32_complete_eigenvalues": [
            ftext(value) for value in complete_values[:32]
        ],
    },
    "status_labels": {
        "method": "KNOWN finite-group Fourier block diagonalization",
        "schedule0_kernel": (
            "DERIVED COMPUTATIONAL / NUMERICALLY CERTIFIED"
            if outcome == "SCHEDULE0_BORDERED_NONSINGULAR_NUMERICALLY_CERTIFIED"
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
print(f"Group-averaging bound: {averaging_bound:.6e}")
print(f"All exhaustive block eigenvalues separated: {all_separated}")
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests} passed")
print(f"Artifact: {OUTPUT.relative_to(ROOT)}")
sys.exit(0 if passed == tests else 1)

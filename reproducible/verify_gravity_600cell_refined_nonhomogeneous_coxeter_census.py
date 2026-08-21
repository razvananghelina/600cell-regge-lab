#!/usr/bin/env python3
"""Exhaustive primary C30 block census for all schedule representatives.

Protocol commit: 2eef0e1.  No continuum or physical target is loaded.
The normal mode recomputes from scratch.  ``--resume`` explicitly continues
an exact-provenance partial checkpoint written by this fixed implementation.
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
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_refined_nonhomogeneous_coxeter_blocks.py"
ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
CURVATURE = HERE / "gravity_600cell_refined_local_curvature_mass.json"
MANIFEST = HERE / "gravity_600cell_refined_nonhomogeneous_csr_manifest.json"
PRIMARY_RESULT = HERE / "gravity_600cell_refined_nonhomogeneous_coxeter_blocks.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_block_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_census_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_nonhomogeneous_coxeter_census.json"

PROTOCOL_COMMIT = "2eef0e1"
SCHEDULE0_ACCEPTANCE_COMMIT = "66e47a7"
EXPECTED_HASHES = {
    "reproducible/verify_gravity_600cell_refined_nonhomogeneous_internal_hessian.py":
        "2a96c8ce466d6d4e9be2cadc4ebf932b4e42eff16fc64fbfb08cd580d680879e",
    "reproducible/verify_gravity_600cell_refined_nonhomogeneous_coxeter_blocks.py":
        "82d29fd5a04eb29cf3fd3c04f95c7fa8489fc0cd6102f2e87f0a3fbb9847a1d4",
    "reproducible/gravity_600cell_refined_local_curvature_mass.json":
        "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "reproducible/gravity_600cell_refined_nonhomogeneous_csr_manifest.json":
        "a6d884a30563bffe1074c964ae27bd2877876fd14a6b4aa169b11ee5eb8a8f1f",
    "reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_blocks.json":
        "640f07a3d13ae3692761243cb62ace3ac2fd38f646d03b7750df05883d3f0267",
    "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_block_prior_art.md":
        "e5db7f7fe28d778a11a40b91115b00d9d4beffe7207ce2d0d5a3bb6e3b2f2018",
    "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_census_protocol.md":
        "0f4699475998a533cff07b6b9806dcb60af7a4ec4fb081fa8a8353dafd103405",
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


def load_primary_functions():
    """Load definitions only; the frozen primary top level is not executed."""
    tree = ast.parse(PRIMARY_SOURCE.read_text(), filename=str(PRIMARY_SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "ast": ast,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "deque": deque,
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
        "SOURCE": SOURCE,
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
        "COXETER_WORD": COXETER_WORD,
        "COXETER_ORDER": COXETER_ORDER,
        "EXPECTED_CYCLES": EXPECTED_CYCLES,
        "EXPECTED_SECTOR_DIMENSION": EXPECTED_SECTOR_DIMENSION,
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(PRIMARY_SOURCE), "exec"), namespace)
    return namespace


def matrix_local_record(source, matrix, tangent, diagnostics, expected_digest):
    operator_error = float(diagnostics["operator_error_row_bound"])
    gradient_gate = 100 * max(
        float(diagnostics["gradient_forward_error"]), np.finfo(float).tiny
    )
    reality_gate = 100 * max(operator_error, np.finfo(float).tiny)
    null_residual = float(np.linalg.norm(matrix @ tangent, ord=np.inf))
    multiplication_error = (
        np.finfo(np.float64).eps
        * source["sparse_row_norm"](matrix)
        * np.linalg.norm(tangent, ord=1)
    )
    null_gate = 100 * (operator_error + multiplication_error)
    actual_digest = source["csr_digest"](matrix)
    stationary = (
        float(diagnostics["gradient_maximum_absolute"]) <= gradient_gate
        and float(diagnostics["gradient_maximum_imaginary"]) <= gradient_gate
    )
    reciprocal = (
        float(diagnostics["raw_hessian_imaginary_row_norm"]) <= reality_gate
        and float(diagnostics["raw_hessian_antisymmetric_row_norm"]) <= reality_gate
    )
    return {
        "csr_sha256": actual_digest,
        "expected_csr_sha256": expected_digest,
        "digest_matches": actual_digest == expected_digest,
        "stationary": bool(stationary),
        "gradient_gate": ftext(gradient_gate),
        "real_symmetric": bool(reciprocal),
        "reality_gate": ftext(reality_gate),
        "product_tangent_residual": ftext(null_residual),
        "product_tangent_gate": ftext(null_gate),
        "product_tangent_is_null": bool(null_residual <= null_gate),
        "operator_error_row_bound": ftext(operator_error),
        "nnz": int(matrix.nnz),
    }


def spectrum_from_primary_artifact(artifact):
    values = []
    for record in artifact["block_census"]["records"]:
        weight = int(record["spectral_weight"])
        for value in record["eigenvalues"]:
            values.extend([float(value)] * weight)
    return np.sort(np.asarray(values, dtype=np.float64))


def write_checkpoint(input_hashes, records, complete=False, outcome="IN_PROGRESS"):
    artifact = {
        "title": "Exhaustive Coxeter census for all refined schedule representatives",
        "protocol_commit": PROTOCOL_COMMIT,
        "schedule0_acceptance_commit": SCHEDULE0_ACCEPTANCE_COMMIT,
        "implementation_sha256": digest(Path(__file__).resolve()),
        "input_sha256": input_hashes,
        "run_complete": bool(complete),
        "scope": {
            "representative_indices": [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14],
            "all_cyclic_sectors_included": True,
            "old_sparse_census_rerun": False,
            "old_sparse_spectral_fields_loaded": False,
            "full_suite_run": False,
            "physical_target_loaded": False,
        },
        "representatives": records,
        "outcome": outcome,
        "tests": {"passed": passed, "total": tests},
    }
    temporary = OUTPUT.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    temporary.replace(OUTPUT)


print("=" * 78)
print("EXHAUSTIVE COXETER CENSUS FOR ALL SCHEDULE REPRESENTATIVES")
print("=" * 78)

actual_hashes = {name: digest(ROOT / name) for name in EXPECTED_HASHES}
provenance_ok = check(
    "the protocol and every frozen input have exact provenance",
    actual_hashes == EXPECTED_HASHES and PROTOCOL_COMMIT == "2eef0e1",
    str(actual_hashes),
)

manifest = json.loads(MANIFEST.read_text())
representative_indices = tuple(manifest["representative_indices"])
pairs = tuple(tuple(pair) for pair in manifest["time_reversal_pairs"])
manifest_schedules = {int(record["index"]): record for record in manifest["schedules"]}
manifest_ok = check(
    "the structural manifest contains 24 digests, 12 reversal pairs and no spectrum",
    manifest["spectral_fields_included"] is False
    and manifest["source_artifact_outcome"] == "LOCAL_EXTENSION_INVALID"
    and len(manifest_schedules) == 24
    and len(set(record["csr_sha256"] for record in manifest_schedules.values())) == 24
    and len(pairs) == 12
    and tuple(pair[0] for pair in pairs) == representative_indices,
)

primary = load_primary_functions()
source = primary["load_source_definitions"]()
actions = source["load_action_definitions"]()
definitions_ok = check(
    "only frozen function definitions are loaded",
    "OUTPUT" not in primary and "OUTPUT" not in source and "OUTPUT" not in actions,
)

_, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
vertex_cells, top, colours = actions["barycentric_chambers"](coarse_top)
spatial_cells = actions["all_simplices"](tuple(map(tuple, top)))
orders = tuple(permutations(range(4)))
schedules = tuple(
    source["schedule_geometry"](actions, top, colours, order) for order in orders
)
simplex_kinds, triangle_kinds = source["global_pattern_catalogue"](schedules)
topology_ok = check(
    "all 24 schedules and their time reversals have the frozen complete carrier",
    tuple(len(layer) for layer in spatial_cells) == EXPECTED_F
    and len(schedules) == 24
    and all(schedule["counts"] == EXPECTED_COUNTS for schedule in schedules)
    and all(tuple(schedules[i]["order"]) == tuple(manifest_schedules[i]["order"])
            for i in range(24))
    and all(tuple(reversed(orders[left])) == orders[right] for left, right in pairs),
    f"simplex patterns={len(simplex_kinds)}, triangle patterns={len(triangle_kinds)}",
)

geometries = {dps: actions["exact_geometry"](dps) for dps in DECIMAL_PRECISIONS}
for geometry in geometries.values():
    geometry["mass"] = mp.mpf(0)
curvature = json.loads(CURVATURE.read_text())
masses = tuple(
    mp.mpf(value) for value in curvature["selected_rank_matter"]["per_vertex_masses"]
)

print("Building the shared target-free local stencil catalogue...", flush=True)
angle_cache, stencil_diagnostics = source["build_angle_cache"](
    actions, simplex_kinds, geometries
)
area_cache = source["build_area_cache"](triangle_kinds, geometries[140])
stencil_ok = check(
    "all shared local stencils remain on the frozen Lorentzian branch",
    stencil_diagnostics["all_displaced_lorentzian"]
    and mp.mpf(stencil_diagnostics["minimum_logarithm_argument"]) > 0
    and mp.mpf(stencil_diagnostics["minimum_gram_leading_minor"]) > 0,
    str(stencil_diagnostics),
)

colour_maps = primary["coloured_chamber_maps"](top)
coxeter_relations = tuple(tuple(
    1 if left == right else primary["permutation_order"](
        primary["compose"](colour_maps[left], colour_maps[right]), 10
    )
    for right in range(4)
) for left in range(4))
expected_relations = (
    (1, 3, 2, 2),
    (3, 1, 3, 2),
    (2, 3, 1, 5),
    (2, 2, 5, 1),
)
left_action, base_target, left_consistent = primary["left_action_from_base"](
    colour_maps, COXETER_WORD
)
cell_map, cells_well_defined, cell_contradictions = primary["induced_cell_map"](
    top, left_action, len(vertex_cells)
)
action_ok = check(
    "the fixed colour-preserving left action realizes the H4 Coxeter element",
    coxeter_relations == expected_relations
    and left_consistent
    and len(set(map(int, left_action))) == len(left_action)
    and primary["permutation_order"](left_action) == COXETER_ORDER
    and cells_well_defined
    and cell_contradictions == 0
    and primary["permutation_order"](cell_map) == COXETER_ORDER
    and np.array_equal(colours[cell_map], colours),
    f"base target={base_target}, relations={coxeter_relations}",
)

right_action = np.arange(len(top), dtype=np.int32)
for colour in COXETER_WORD:
    right_action = primary["compose"](colour_maps[colour], right_action)
_, right_well_defined, right_contradictions = primary["induced_cell_map"](
    top, right_action, len(vertex_cells)
)
right_control_ok = check(
    "the forbidden right-product convention fails cell well-definedness",
    not right_well_defined and right_contradictions > 0,
    f"contradictions={right_contradictions}",
)

resume = "--resume" in sys.argv[1:]
if any(argument != "--resume" for argument in sys.argv[1:]):
    raise SystemExit("usage: verifier [--resume]")
records = []
if resume and OUTPUT.exists():
    checkpoint = json.loads(OUTPUT.read_text())
    checkpoint_indices = [record["representative_schedule"] for record in checkpoint.get("representatives", [])]
    if (
        checkpoint.get("protocol_commit") != PROTOCOL_COMMIT
        or checkpoint.get("implementation_sha256") != digest(Path(__file__).resolve())
        or checkpoint.get("input_sha256") != actual_hashes
        or checkpoint_indices != list(representative_indices[:len(checkpoint_indices)])
    ):
        raise RuntimeError("the requested checkpoint does not have exact frozen provenance")
    records = checkpoint["representatives"]
    print(f"Resuming after {len(records)} exact-provenance representatives.", flush=True)

primary_schedule0 = None

for pair_number, (forward_index, reverse_index) in enumerate(pairs):
    if pair_number < len(records):
        record = records[pair_number]
        check(
            f"pair {pair_number}: cached source/action controls passed",
            record["construction_controls_passed"],
        )
        check(
            f"pair {pair_number}: cached exhaustive block construction completed",
            record["block_construction_passed"],
            f"all separated={record['all_eigenvalues_separated']}",
        )
        continue

    print("-" * 78, flush=True)
    print(
        f"Pair {pair_number + 1}/12: schedules {forward_index}/{reverse_index}, "
        f"order={orders[forward_index]}",
        flush=True,
    )
    pair_data = []
    local_records = []
    for schedule_index in (forward_index, reverse_index):
        schedule = schedules[schedule_index]
        print(f"Assembling complete schedule {schedule_index}...", flush=True)
        matrix, diagnostics = source["assemble_internal"](
            schedule, angle_cache, area_cache, masses, geometries[140]
        )
        tangent = source["product_tangent"](schedule, geometries[140])
        local_record = matrix_local_record(
            source, matrix, tangent, diagnostics,
            manifest_schedules[schedule_index]["csr_sha256"],
        )
        pair_data.append((schedule, matrix, tangent, diagnostics))
        local_records.append(local_record)

    forward, reverse = pair_data
    reversal = source["reversal_permutation"](
        forward[0], reverse[0], len(colours)
    )
    reverse_pulled = reverse[1][reversal, :][:, reversal]
    reversal_difference = source["sparse_row_norm"](forward[1] - reverse_pulled)
    reversal_gate = 100 * (
        float(forward[3]["operator_error_row_bound"])
        + float(reverse[3]["operator_error_row_bound"])
    )
    reversal_ok = reversal_difference <= reversal_gate

    schedule, matrix, tangent, diagnostics = forward
    edge_map, edges_invariant = primary["induced_internal_edge_map"](
        source, schedule, cell_map, len(vertex_cells)
    )
    cycles, cycle_id, cycle_position = primary["cycle_coordinates"](edge_map)
    cycle_histogram = dict(sorted(Counter(map(len, cycles)).items()))
    sector_dimensions = [
        sum(1 for cycle in cycles if (character * len(cycle)) % COXETER_ORDER == 0)
        for character in range(COXETER_ORDER)
    ]
    cycles_ok = (
        edges_invariant
        and primary["permutation_order"](edge_map) == COXETER_ORDER
        and cycle_histogram == {30: EXPECTED_CYCLES}
        and sector_dimensions == [EXPECTED_SECTOR_DIMENSION] * COXETER_ORDER
        and sum(sector_dimensions) == len(tangent)
    )

    operator_error = float(diagnostics["operator_error_row_bound"])
    covariance_values, averaging_bound = primary["covariance_defects"](
        source, matrix, edge_map
    )
    covariance_ok = averaging_bound <= 100 * operator_error
    corrupted = matrix.copy().tolil()
    corrupted[0, 0] += 1e-4
    corrupted = corrupted.tocsr()
    _, corrupted_bound = primary["covariance_defects"](source, corrupted, edge_map)
    corruption_ok = corrupted_bound > 100 * operator_error

    tangent_projections = [
        primary["tangent_projection"](tangent, cycles, character)
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
    tangent_sector_ok = tangent_noninvariant <= tangent_gate

    construction_controls_passed = all(
        record["digest_matches"]
        and record["stationary"]
        and record["real_symmetric"]
        and record["product_tangent_is_null"]
        for record in local_records
    ) and reversal_ok and cycles_ok and covariance_ok and corruption_ok and tangent_sector_ok
    check(
        f"pair {pair_number}: both sources, reversal, action and covariance controls pass",
        construction_controls_passed,
        f"reversal={reversal_difference:.3e}/{reversal_gate:.3e}, "
        f"covariance={averaging_bound:.3e}, corruption={corrupted_bound:.3e}",
    )

    block_records = []
    blocks = {}
    eigenvalue_sets = {}
    all_diagonalized = True
    all_separated = bool(construction_controls_passed)
    maximum_block_roundoff = 0.0
    maximum_eigen_residual = 0.0
    zero_candidates = []
    matrix_coo = matrix.tocoo()
    matrix_coo.sum_duplicates()
    common_error = operator_error + averaging_bound

    if construction_controls_passed:
        for character in range(16):
            print(
                f"Pair {pair_number + 1}/12: exhaustive sector k={character}/15...",
                flush=True,
            )
            block, entry_error, block_diagnostics = primary["fourier_block"](
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
                tangent_error = (gamma30 + 4 * unit_roundoff) * np.asarray([
                    np.sum(np.abs(tangent[cycle])) / np.sqrt(COXETER_ORDER)
                    for cycle in cycles
                ])
                bordered_error = np.zeros((657, 657), dtype=np.float64)
                bordered_error[:656, :656] = entry_error
                bordered_error[:656, 656] = tangent_error
                bordered_error[656, :656] = tangent_error
                entry_error = bordered_error
                block_diagnostics["roundoff_row_bound"] = ftext(
                    float(np.max(np.sum(entry_error, axis=1)))
                )

            block_roundoff = float(block_diagnostics["roundoff_row_bound"])
            maximum_block_roundoff = max(maximum_block_roundoff, block_roundoff)
            try:
                values, vectors, residuals, separated, nearest = primary["diagonalize_block"](
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
            gates = 100 * (common_error + block_roundoff + residuals)
            ratios = np.abs(values) / gates
            margin_index = int(np.argmin(ratios))
            maximum_eigen_residual = max(maximum_eigen_residual, float(residuals.max()))
            all_separated &= bool(np.all(separated))
            blocks[character] = block
            eigenvalue_sets[character] = values
            block_records.append({
                "character": character,
                "dimension": len(block),
                "spectral_weight": 1 if character in (0, 15) else 2,
                "minimum_absolute_eigenvalue": ftext(np.min(np.abs(values))),
                "minimum_margin_ratio": ftext(ratios[margin_index]),
                "minimum_margin_eigenvalue": ftext(values[margin_index]),
                "minimum_margin_gate": ftext(gates[margin_index]),
                "nearest_index": nearest,
                "nearest_eigenvalue": ftext(values[nearest]),
                "nearest_residual": ftext(residuals[nearest]),
                "nearest_gate": ftext(gates[nearest]),
                "all_eigenvalues_separated": bool(np.all(separated)),
                "maximum_eigenpair_residual": ftext(residuals.max()),
                "eigenvalues": [ftext(value) for value in values],
                "diagnostics": block_diagnostics,
            })
            if not np.all(separated):
                for index in np.flatnonzero(~separated):
                    full = primary["reconstruct_sector_vector"](
                        vectors[:, index], cycles, character,
                        bordered=(character == 0),
                    )
                    if character == 0:
                        column = sp.csr_matrix(tangent.reshape(-1, 1))
                        full_matrix = sp.bmat(
                            [[matrix, column], [column.T, sp.csr_matrix((1, 1))]],
                            format="csr",
                        )
                    else:
                        full_matrix = matrix
                    residual = np.linalg.norm(full_matrix @ full - values[index] * full)
                    zero_candidates.append({
                        "character": character,
                        "eigenvalue": ftext(values[index]),
                        "gate": ftext(gates[index]),
                        "block_residual": ftext(residuals[index]),
                        "full_residual_against_unaveraged_matrix": ftext(residual),
                    })

    dimension_ok = (
        all_diagonalized
        and len(block_records) == 16
        and sum(
            record["dimension"] * record.get("spectral_weight", 0)
            for record in block_records
        ) == 19681
    )
    wrong_weight_dimension = sum(
        record["dimension"]
        * (0 if record.get("character") == 15 else record.get("spectral_weight", 0))
        for record in block_records
    )
    omission_control_ok = dimension_ok and wrong_weight_dimension != 19681
    if dimension_ok:
        block29, _, diagnostics29 = primary["fourier_block"](
            matrix_coo, cycle_id, cycle_position, 29, EXPECTED_SECTOR_DIMENSION
        )
        record1 = next(record for record in block_records if record["character"] == 1)
        conjugate_difference = float(np.max(np.sum(
            np.abs(block29 - blocks[1].conj()), axis=1
        )))
        conjugate_gate = 100 * (
            float(diagnostics29["roundoff_row_bound"])
            + float(record1["diagnostics"]["roundoff_row_bound"])
            + averaging_bound
        )
        conjugate_ok = conjugate_difference <= conjugate_gate

        trace_full = float(matrix.diagonal().sum())
        trace_blocks = sum(
            record["spectral_weight"] * float(np.sum(eigenvalue_sets[record["character"]]))
            for record in block_records
        )
        trace_difference = abs(trace_full - trace_blocks)
        trace_gate = 100 * (
            19681 * np.finfo(float).eps * max(1.0, abs(trace_full))
            + 19681 * maximum_eigen_residual
        )
        trace_ok = trace_difference <= trace_gate

        frobenius_full = np.sqrt(
            float(np.sum(matrix.data ** 2)) + 2 * float(np.dot(tangent, tangent))
        )
        frobenius_blocks = np.sqrt(sum(
            record["spectral_weight"]
            * float(np.linalg.norm(blocks[record["character"]], ord="fro") ** 2)
            for record in block_records
        ))
        frobenius_difference = abs(frobenius_full - frobenius_blocks)
        frobenius_gate = 100 * (
            np.sqrt(19681) * averaging_bound
            + np.sqrt(19681) * maximum_block_roundoff
            + np.finfo(float).eps * max(frobenius_full, 1.0)
        )
        parseval_ok = frobenius_difference <= frobenius_gate
    else:
        conjugate_difference = np.inf
        conjugate_gate = 0.0
        conjugate_ok = False
        trace_full = trace_blocks = trace_difference = np.nan
        trace_gate = 0.0
        trace_ok = False
        frobenius_full = frobenius_blocks = frobenius_difference = np.nan
        frobenius_gate = 0.0
        parseval_ok = False

    full_values = []
    for character in range(16):
        weight = 1 if character in (0, 15) else 2
        for value in eigenvalue_sets.get(character, ()):
            full_values.extend([float(value)] * weight)
    full_values = np.sort(np.asarray(full_values, dtype=np.float64))

    schedule0_comparison = None
    schedule0_ok = True
    if forward_index == 0 and len(full_values) == 19681:
        # The frozen primary result is deliberately parsed only now, after the
        # new complete schedule-0 spectrum has been constructed.
        primary_schedule0 = json.loads(PRIMARY_RESULT.read_text())
        primary_values = spectrum_from_primary_artifact(primary_schedule0)
        primary_common = float(
            primary_schedule0["block_census"]["common_operator_and_averaging_error"]
        )
        primary_roundoff = float(
            primary_schedule0["block_census"]["maximum_block_roundoff"]
        )
        primary_residual = float(
            primary_schedule0["block_census"]["maximum_eigenpair_residual"]
        )
        comparison_gate = 100 * (
            common_error + maximum_block_roundoff + maximum_eigen_residual
            + primary_common + primary_roundoff + primary_residual
        )
        comparison_difference = float(np.max(np.abs(full_values - primary_values)))
        schedule0_ok = comparison_difference <= comparison_gate
        schedule0_comparison = {
            "maximum_difference": ftext(comparison_difference),
            "gate": ftext(comparison_gate),
            "passes": bool(schedule0_ok),
        }

    block_construction_passed = all((
        construction_controls_passed,
        dimension_ok,
        omission_control_ok,
        conjugate_ok,
        trace_ok,
        parseval_ok,
        schedule0_ok,
    ))
    check(
        f"pair {pair_number}: all exhaustive blocks and global invariants complete",
        block_construction_passed,
        f"all separated={all_separated}, zero candidates={len(zero_candidates)}",
    )

    record = {
        "pair_number": pair_number,
        "representative_schedule": forward_index,
        "reverse_schedule": reverse_index,
        "representative_order": list(orders[forward_index]),
        "reverse_order": list(orders[reverse_index]),
        "source_matrices": {
            "representative": local_records[0],
            "reverse": local_records[1],
            "reversal_maximum_row_difference": ftext(reversal_difference),
            "reversal_gate": ftext(reversal_gate),
            "reversal_congruent": bool(reversal_ok),
        },
        "coxeter_action": {
            "order": primary["permutation_order"](edge_map),
            "cycle_histogram": {str(key): value for key, value in cycle_histogram.items()},
            "sector_dimensions": sector_dimensions,
            "invariant": bool(edges_invariant),
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
            "weighted_dimension": len(full_values),
            "k15_omitted_weighted_dimension": wrong_weight_dimension,
            "sector_omission_detected": bool(omission_control_ok),
            "common_operator_and_averaging_error": ftext(common_error),
            "maximum_block_roundoff": ftext(maximum_block_roundoff),
            "maximum_eigenpair_residual": ftext(maximum_eigen_residual),
            "records": block_records,
            "zero_candidates": zero_candidates,
            "conjugate_k1_k29_difference": ftext(conjugate_difference),
            "conjugate_gate": ftext(conjugate_gate),
            "trace_full": ftext(trace_full),
            "trace_blocks": ftext(trace_blocks),
            "trace_difference": ftext(trace_difference),
            "trace_gate": ftext(trace_gate),
            "frobenius_full": ftext(frobenius_full),
            "frobenius_blocks": ftext(frobenius_blocks),
            "frobenius_difference": ftext(frobenius_difference),
            "frobenius_gate": ftext(frobenius_gate),
        },
        "schedule0_primary_comparison": schedule0_comparison,
        "construction_controls_passed": bool(construction_controls_passed),
        "block_construction_passed": bool(block_construction_passed),
        "all_eigenvalues_separated": bool(all_separated),
    }
    records.append(record)
    write_checkpoint(actual_hashes, records)

all_construction = (
    provenance_ok and manifest_ok and definitions_ok and topology_ok and stencil_ok
    and action_ok and right_control_ok
    and len(records) == 12
    and all(record["construction_controls_passed"] for record in records)
    and all(record["block_construction_passed"] for record in records)
)
all_separated = len(records) == 12 and all(
    record["all_eigenvalues_separated"] for record in records
)
if not all_construction:
    outcome = "COXETER_CENSUS_CONSTRUCTION_INVALID"
elif all_separated:
    outcome = "ALL_12_REPRESENTATIVE_KERNELS_ARE_PRODUCT_DURATION_LINES_PRIMARY"
else:
    outcome = "ADDITIONAL_INTERNAL_KERNEL_CANDIDATE_FOUND"

verdict_ok = check(
    "the complete census follows the preregistered verdict hierarchy",
    outcome in {
        "COXETER_CENSUS_CONSTRUCTION_INVALID",
        "ALL_12_REPRESENTATIVE_KERNELS_ARE_PRODUCT_DURATION_LINES_PRIMARY",
        "ADDITIONAL_INTERNAL_KERNEL_CANDIDATE_FOUND",
    }
    and (outcome == "COXETER_CENSUS_CONSTRUCTION_INVALID") == (not all_construction)
    and (
        outcome == "ALL_12_REPRESENTATIVE_KERNELS_ARE_PRODUCT_DURATION_LINES_PRIMARY"
    ) == (all_construction and all_separated),
    outcome,
)

write_checkpoint(actual_hashes, records, complete=True, outcome=outcome)
artifact = json.loads(OUTPUT.read_text())
artifact["global_construction"] = {
    "spatial_f_vector": list(EXPECTED_F),
    "schedule_count": len(schedules),
    "representative_count": len(records),
    "simplex_pattern_count": len(simplex_kinds),
    "triangle_pattern_count": len(triangle_kinds),
    "stencil_diagnostics": stencil_diagnostics,
    "coxeter_relations": [list(row) for row in coxeter_relations],
    "base_target": base_target,
    "right_convention_well_defined": bool(right_well_defined),
    "right_convention_contradictions": right_contradictions,
}
artifact["status_labels"] = {
    "method": "KNOWN finite-group Fourier block diagonalization",
    "schedule0_kernel": "DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED",
    "other_schedule_kernels": (
        "DERIVED COMPUTATIONAL / PRIMARY ONLY"
        if outcome == "ALL_12_REPRESENTATIVE_KERNELS_ARE_PRODUCT_DURATION_LINES_PRIMARY"
        else "OPEN"
    ),
    "physical_tick": "OPEN / NOT SELECTED",
    "propagation_c_G_planck": "OPEN / NOT TESTED",
    "external_novelty": "OPEN",
}
artifact["tests"] = {"passed": passed, "total": tests}
temporary = OUTPUT.with_suffix(".json.tmp")
temporary.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
temporary.replace(OUTPUT)

print("-" * 78)
print(f"Representatives completed: {len(records)}/12")
print(f"All construction controls passed: {all_construction}")
print(f"All bordered eigenvalues separated: {all_separated}")
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests} passed")
print(f"Artifact: {OUTPUT.relative_to(ROOT)}")
sys.exit(0 if passed == tests and verdict_ok else 1)

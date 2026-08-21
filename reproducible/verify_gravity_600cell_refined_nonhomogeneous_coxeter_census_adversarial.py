#!/usr/bin/env python3
"""Geometric adversarial replication of the complete Coxeter census.

Protocol commit: 3e8e797.  Primary spectra are not parsed until every one of
the eleven new adversarial spectra has been constructed.
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
ADVERSARIAL_SOURCE = HERE / "verify_gravity_600cell_refined_nonhomogeneous_coxeter_blocks_adversarial.py"
ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
CURVATURE = HERE / "gravity_600cell_refined_local_curvature_mass.json"
MANIFEST = HERE / "gravity_600cell_refined_nonhomogeneous_csr_manifest.json"
PRIMARY = HERE / "gravity_600cell_refined_nonhomogeneous_coxeter_census.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_census_adversarial_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_nonhomogeneous_coxeter_census_adversarial.json"

PROTOCOL_COMMIT = "3e8e797"
PRIMARY_RESULT_COMMIT = "d0adbac"
EXPECTED_HASHES = {
    "reproducible/verify_gravity_600cell_refined_nonhomogeneous_internal_hessian.py":
        "2a96c8ce466d6d4e9be2cadc4ebf932b4e42eff16fc64fbfb08cd580d680879e",
    "reproducible/verify_gravity_600cell_refined_nonhomogeneous_coxeter_blocks_adversarial.py":
        "40d7bba8f339b466ba3b1157dfe42ff7925b5133d140842312819d758e9bb48d",
    "reproducible/gravity_600cell_refined_local_curvature_mass.json":
        "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "reproducible/gravity_600cell_refined_nonhomogeneous_csr_manifest.json":
        "a6d884a30563bffe1074c964ae27bd2877876fd14a6b4aa169b11ee5eb8a8f1f",
    "reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_census.json":
        "5a50b8179fc272d75c0811dcf34fb8c6a464e564729956bdd6baaa82a6f058b6",
    "docs/gravity/gravity_600cell_refined_nonhomogeneous_coxeter_census_adversarial_protocol.md":
        "e34fee430cf141b2cbadeb7629b87bd88d06e8e071a6a5e734e89fe61ef13dd9",
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
REPRESENTATIVES = (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14)

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


def load_adversarial_functions():
    """Load definitions only; never execute the frozen schedule-0 top level."""
    tree = ast.parse(
        ADVERSARIAL_SOURCE.read_text(), filename=str(ADVERSARIAL_SOURCE)
    )
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
        "la": la,
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
        "WORD": WORD,
        "ORDER": ORDER,
        "DIMENSION": DIMENSION,
        "MATCH_GATE": MATCH_GATE,
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(ADVERSARIAL_SOURCE), "exec"), namespace)
    return namespace


def weighted_spectrum(record):
    values = []
    for block in record["block_census"]["records"]:
        weight = int(block["spectral_weight"])
        for value in block["eigenvalues"]:
            values.extend([float(value)] * weight)
    return np.sort(np.asarray(values, dtype=np.float64))


def write_checkpoint(input_hashes, geometry_record, records,
                     complete=False, outcome="IN_PROGRESS"):
    artifact = {
        "title": "Geometric adversarial replication of the complete Coxeter census",
        "protocol_commit": PROTOCOL_COMMIT,
        "primary_result_commit": PRIMARY_RESULT_COMMIT,
        "implementation_sha256": digest(Path(__file__).resolve()),
        "input_sha256": input_hashes,
        "run_complete": bool(complete),
        "scope": {
            "schedule_indices": list(REPRESENTATIVES),
            "all_cyclic_sectors_included": True,
            "primary_block_code_loaded": False,
            "primary_spectra_parsed_after_blind_phase": bool(complete),
            "old_sparse_census_rerun": False,
            "old_sparse_spectral_fields_loaded": False,
            "full_suite_run": False,
            "physical_target_loaded": False,
        },
        "geometric_action": geometry_record,
        "representatives": records,
        "outcome": outcome,
        "tests": {"passed": passed, "total": tests},
    }
    temporary = OUTPUT.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    temporary.replace(OUTPUT)


print("=" * 78)
print("GEOMETRIC ADVERSARIAL REPLICATION OF THE COMPLETE COXETER CENSUS")
print("=" * 78)

actual_hashes = {name: digest(ROOT / name) for name in EXPECTED_HASHES}
provenance_ok = check(
    "the adversarial protocol and every frozen input have exact provenance",
    actual_hashes == EXPECTED_HASHES and PROTOCOL_COMMIT == "3e8e797",
    str(actual_hashes),
)

manifest = json.loads(MANIFEST.read_text())
manifest_schedules = {int(record["index"]): record for record in manifest["schedules"]}
manifest_ok = check(
    "the structural manifest fixes all eleven new source matrices without spectra",
    manifest["spectral_fields_included"] is False
    and tuple(i for i in manifest["representative_indices"] if i != 0) == REPRESENTATIVES
    and all(index in manifest_schedules for index in REPRESENTATIVES),
)

adversary = load_adversarial_functions()
source = adversary["load_source_definitions"]()
actions = source["load_action_definitions"]()
definitions_ok = check(
    "no primary action, block or diagonalization implementation is loaded",
    "OUTPUT" not in adversary and "OUTPUT" not in source and "OUTPUT" not in actions,
)

vertices, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
vertex_cells, top, colours = actions["barycentric_chambers"](coarse_top)
spatial_cells = actions["all_simplices"](tuple(map(tuple, top)))
orders = tuple(permutations(range(4)))
schedules = tuple(
    source["schedule_geometry"](actions, top, colours, order) for order in orders
)
simplex_kinds, triangle_kinds = source["global_pattern_catalogue"](schedules)
topology_ok = check(
    "the complete carrier and all schedule geometries have the frozen census",
    tuple(len(layer) for layer in spatial_cells) == EXPECTED_F
    and len(schedules) == 24
    and all(schedule["counts"] == EXPECTED_COUNTS for schedule in schedules)
    and all(tuple(schedules[i]["order"]) == tuple(manifest_schedules[i]["order"])
            for i in range(24)),
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
    "the shared local stencils remain on the frozen Lorentzian branch",
    stencil_diagnostics["all_displaced_lorentzian"]
    and mp.mpf(stencil_diagnostics["minimum_logarithm_argument"]) > 0
    and mp.mpf(stencil_diagnostics["minimum_gram_leading_minor"]) > 0,
    str(stencil_diagnostics),
)

centres = adversary["cell_centres"](vertices, vertex_cells)
reflections, normals = adversary["chamber_reflections"](centres, top[0])
pair_orders = tuple(tuple(
    1 if left == right else adversary["matrix_order"](
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
single_reflection_order = adversary["matrix_order"](reflections[0])
coxeter = adversary["sequential_product"](reflections, WORD)
reverse_coxeter = reflections[0] @ reflections[1] @ reflections[2] @ reflections[3]
vertex_map, vertex_residual, vertex_bijective = adversary["vertex_permutation"](
    vertices, coxeter
)
reverse_vertex_map, reverse_vertex_residual, reverse_vertex_bijective = (
    adversary["vertex_permutation"](vertices, reverse_coxeter)
)
cell_map, cell_bijective, missing_cells = adversary["induced_cell_permutation"](
    vertex_cells, vertex_map
)
reverse_cell_map, reverse_cell_bijective, reverse_missing_cells = (
    adversary["induced_cell_permutation"](vertex_cells, reverse_vertex_map)
)
geometry_ok = check(
    "ambient R4 reflections independently reproduce H4 and both order-30 actions",
    pair_orders == expected_orders
    and reflection_residual <= MATCH_GATE
    and single_reflection_order == 2
    and adversary["matrix_order"](coxeter) == ORDER
    and adversary["matrix_order"](reverse_coxeter) == ORDER
    and vertex_bijective and reverse_vertex_bijective
    and vertex_residual <= MATCH_GATE and reverse_vertex_residual <= MATCH_GATE
    and cell_bijective and reverse_cell_bijective
    and missing_cells == 0 and reverse_missing_cells == 0
    and adversary["permutation_order"](cell_map) == ORDER
    and adversary["permutation_order"](reverse_cell_map) == ORDER
    and np.array_equal(colours[cell_map], colours)
    and np.array_equal(colours[reverse_cell_map], colours),
    f"pair orders={pair_orders}, matching residual={vertex_residual:.3e}",
)
negative_order_ok = check(
    "the single-reflection negative control has order two rather than thirty",
    single_reflection_order == 2 and single_reflection_order != ORDER,
)

geometry_record = {
    "reflection_normals": [[ftext(value) for value in normal] for normal in normals],
    "pair_orders": [list(row) for row in pair_orders],
    "single_reflection_order": single_reflection_order,
    "coxeter_matrix_order": adversary["matrix_order"](coxeter),
    "reverse_coxeter_matrix_order": adversary["matrix_order"](reverse_coxeter),
    "vertex_matching_residual": ftext(vertex_residual),
    "reverse_vertex_matching_residual": ftext(reverse_vertex_residual),
    "cell_order": adversary["permutation_order"](cell_map),
    "reverse_cell_order": adversary["permutation_order"](reverse_cell_map),
}

resume = "--resume" in sys.argv[1:]
if any(argument != "--resume" for argument in sys.argv[1:]):
    raise SystemExit("usage: verifier [--resume]")
records = []
if resume and OUTPUT.exists():
    checkpoint = json.loads(OUTPUT.read_text())
    checkpoint_indices = [record["schedule_index"] for record in checkpoint.get("representatives", [])]
    if (
        checkpoint.get("protocol_commit") != PROTOCOL_COMMIT
        or checkpoint.get("implementation_sha256") != digest(Path(__file__).resolve())
        or checkpoint.get("input_sha256") != actual_hashes
        or checkpoint_indices != list(REPRESENTATIVES[:len(checkpoint_indices)])
        or checkpoint.get("scope", {}).get("primary_spectra_parsed_after_blind_phase")
    ):
        raise RuntimeError("the requested checkpoint is not an exact blind-phase prefix")
    records = checkpoint["representatives"]
    print(f"Resuming after {len(records)} exact-provenance representatives.", flush=True)

for ordinal, schedule_index in enumerate(REPRESENTATIVES):
    if ordinal < len(records):
        record = records[ordinal]
        check(
            f"schedule {schedule_index}: cached geometry and matrix controls pass",
            record["construction_controls_passed"],
        )
        check(
            f"schedule {schedule_index}: cached explicit-Q census is complete",
            record["block_construction_passed"],
            f"all separated={record['all_eigenvalues_separated']}",
        )
        continue

    schedule = schedules[schedule_index]
    print("-" * 78, flush=True)
    print(
        f"Adversarial representative {ordinal + 1}/11: schedule {schedule_index}, "
        f"order={orders[schedule_index]}",
        flush=True,
    )
    matrix, matrix_diagnostics = source["assemble_internal"](
        schedule, angle_cache, area_cache, masses, geometries[140]
    )
    matrix_digest = source["csr_digest"](matrix)
    matrix_ok = matrix_digest == manifest_schedules[schedule_index]["csr_sha256"]
    operator_error = float(matrix_diagnostics["operator_error_row_bound"])
    tangent = source["product_tangent"](schedule, geometries[140])
    null_residual = float(np.linalg.norm(matrix @ tangent, ord=np.inf))
    null_multiplication_error = (
        np.finfo(np.float64).eps
        * source["sparse_row_norm"](matrix)
        * np.linalg.norm(tangent, ord=1)
    )
    null_gate = 100 * (operator_error + null_multiplication_error)

    edge_map, edge_invariant = adversary["induced_edge_permutation"](
        source, schedule, cell_map, len(vertex_cells)
    )
    reverse_edge_map, reverse_edge_invariant = adversary["induced_edge_permutation"](
        source, schedule, reverse_cell_map, len(vertex_cells)
    )
    inverse_relation = np.array_equal(
        reverse_edge_map[edge_map], np.arange(len(edge_map), dtype=np.int32)
    )
    cycles, cycle_id, cycle_position = adversary["cycle_coordinates"](edge_map)
    reverse_cycles, reverse_cycle_id, reverse_cycle_position = adversary["cycle_coordinates"](
        reverse_edge_map
    )
    cycle_histogram = dict(sorted(Counter(map(len, cycles)).items()))
    reverse_cycle_histogram = dict(sorted(Counter(map(len, reverse_cycles)).items()))
    action_ok = (
        edge_invariant and reverse_edge_invariant and inverse_relation
        and adversary["permutation_order"](edge_map) == ORDER
        and adversary["permutation_order"](reverse_edge_map) == ORDER
        and cycle_histogram == {30: DIMENSION}
        and reverse_cycle_histogram == {30: DIMENSION}
    )

    covariance_values, averaging_bound = adversary["covariance_defects"](
        source, matrix, edge_map
    )
    reverse_covariance_values, reverse_averaging_bound = adversary["covariance_defects"](
        source, matrix, reverse_edge_map
    )
    covariance_ok = (
        averaging_bound <= 100 * operator_error
        and reverse_averaging_bound <= 100 * operator_error
    )
    corrupted = matrix.copy().tolil()
    corrupted[0, 0] += 1e-4
    corrupted = corrupted.tocsr()
    _, corrupted_bound = adversary["covariance_defects"](
        source, corrupted, edge_map
    )
    corruption_ok = corrupted_bound > 100 * operator_error
    construction_controls_passed = all((
        matrix_ok,
        null_residual <= null_gate,
        action_ok,
        covariance_ok,
        corruption_ok,
    ))
    check(
        f"schedule {schedule_index}: source, geometric action and covariance controls pass",
        construction_controls_passed,
        f"covariance={averaging_bound:.3e}, corruption={corrupted_bound:.3e}",
    )

    support_counts, absolute_bound = adversary["sparse_support_products"](
        matrix, cycle_id
    )
    reverse_support_counts, reverse_absolute_bound = adversary["sparse_support_products"](
        matrix, reverse_cycle_id
    )
    block_records = []
    blocks = {}
    eigenvalue_sets = {}
    all_diagonalized = True
    all_isometric = True
    maximum_isometry_defect = 0.0
    all_separated = bool(construction_controls_passed)
    maximum_roundoff = 0.0
    maximum_lifted_residual = 0.0
    maximum_noninvariant_tangent = 0.0
    zero_candidates = []

    if construction_controls_passed:
        for character in range(16):
            print(
                f"Schedule {schedule_index}: explicit-Q evd sector k={character}/15...",
                flush=True,
            )
            isometry = adversary["explicit_fourier_isometry"](
                cycle_id, cycle_position, character
            )
            gram = (isometry.conj().T @ isometry).toarray()
            isometry_defect = float(np.linalg.norm(
                gram - np.eye(DIMENSION), ord=np.inf
            ))
            maximum_isometry_defect = max(
                maximum_isometry_defect, isometry_defect
            )
            all_isometric &= isometry_defect <= 1e-12
            projection = np.asarray(isometry.conj().T @ tangent).ravel()
            if character:
                maximum_noninvariant_tangent = max(
                    maximum_noninvariant_tangent, float(np.max(np.abs(projection)))
                )
            block, product, diagnostics = adversary["sparse_product_block"](
                matrix, isometry, support_counts, absolute_bound
            )
            bordered = character == 0
            if bordered:
                complete_block = np.zeros((657, 657), dtype=np.complex128)
                complete_block[:656, :656] = block
                complete_block[:656, 656] = projection
                complete_block[656, :656] = projection.conj()
                block = complete_block
            try:
                values, vectors = la.eigh(
                    block, driver="evd", check_finite=True, overwrite_a=False
                )
                residuals = adversary["lifted_residuals"](
                    matrix, isometry, product, values, vectors,
                    bordered=bordered, tangent=tangent if bordered else None,
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
            ratios = np.abs(values) / gates
            margin_index = int(np.argmin(ratios))
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
                "minimum_margin_ratio": ftext(ratios[margin_index]),
                "minimum_margin_eigenvalue": ftext(values[margin_index]),
                "minimum_margin_gate": ftext(gates[margin_index]),
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

        # The remaining sectors are conjugate for the real source matrix, but
        # the protocol requires every Q_k itself to be constructed and tested.
        for character in range(16, ORDER):
            isometry = adversary["explicit_fourier_isometry"](
                cycle_id, cycle_position, character
            )
            gram = (isometry.conj().T @ isometry).toarray()
            isometry_defect = float(np.linalg.norm(
                gram - np.eye(DIMENSION), ord=np.inf
            ))
            maximum_isometry_defect = max(
                maximum_isometry_defect, isometry_defect
            )
            all_isometric &= isometry_defect <= 1e-12
            projection = np.asarray(isometry.conj().T @ tangent).ravel()
            maximum_noninvariant_tangent = max(
                maximum_noninvariant_tangent,
                float(np.max(np.abs(projection))),
            )

    dimension_ok = (
        all_diagonalized and all_isometric and len(block_records) == 16
        and sum(
            record.get("dimension", 0) * record.get("spectral_weight", 0)
            for record in block_records
        ) == 19681
    )
    tangent_unit_roundoff = np.finfo(np.float64).eps / 2
    tangent_gamma = 30 * tangent_unit_roundoff / (1 - 30 * tangent_unit_roundoff)
    tangent_gate = 100 * tangent_gamma * max(
        float(np.max([np.sum(np.abs(tangent[cycle])) for cycle in cycles])),
        np.finfo(float).tiny,
    )
    tangent_sector_ok = maximum_noninvariant_tangent <= tangent_gate

    if dimension_ok:
        reverse_q1 = adversary["explicit_fourier_isometry"](
            reverse_cycle_id, reverse_cycle_position, 1
        )
        reverse_block1, _, reverse_diagnostics = adversary["sparse_product_block"](
            matrix, reverse_q1, reverse_support_counts, reverse_absolute_bound
        )
        reverse_values1 = la.eigvalsh(reverse_block1, driver="evd")
        reverse_sector_difference = float(np.max(np.abs(
            reverse_values1 - eigenvalue_sets[1]
        )))
        reverse_sector_gate = 100 * (
            operator_error + averaging_bound + reverse_averaging_bound
            + maximum_roundoff
            + float(reverse_diagnostics["sparse_product_roundoff_bound"])
        )
        reverse_ok = reverse_sector_difference <= reverse_sector_gate

        trace_full = float(matrix.diagonal().sum())
        trace_blocks = sum(
            record["spectral_weight"]
            * float(np.sum(eigenvalue_sets[record["character"]]))
            for record in block_records
        )
        trace_difference = abs(trace_full - trace_blocks)
        trace_gate = 100 * (
            19681 * maximum_lifted_residual
            + np.finfo(float).eps * max(abs(trace_full), 1.0)
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
            np.sqrt(19681) * (
                averaging_bound + maximum_roundoff + maximum_lifted_residual
            )
            + np.finfo(float).eps * max(frobenius_full, 1.0)
        )
        parseval_ok = frobenius_difference <= frobenius_gate
    else:
        reverse_sector_difference = np.inf
        reverse_sector_gate = 0.0
        reverse_ok = False
        trace_full = trace_blocks = trace_difference = np.nan
        trace_gate = 0.0
        trace_ok = False
        frobenius_full = frobenius_blocks = frobenius_difference = np.nan
        frobenius_gate = 0.0
        parseval_ok = False

    wrong_phase_record = None
    wrong_phase_ok = True
    if schedule_index == REPRESENTATIVES[0] and dimension_ok:
        wrong_q1 = adversary["explicit_fourier_isometry"](
            cycle_id, cycle_position, 1, sign=1
        )
        wrong_block1, _, wrong_diagnostics = adversary["sparse_product_block"](
            matrix, wrong_q1, support_counts, absolute_bound
        )
        wrong_difference = float(np.linalg.norm(
            wrong_block1 - blocks[1], ord=np.inf
        ))
        wrong_gate = 100 * (
            operator_error + averaging_bound + maximum_roundoff
            + float(wrong_diagnostics["sparse_product_roundoff_bound"])
        )
        wrong_phase_ok = wrong_difference > wrong_gate
        wrong_phase_record = {
            "same_sector_difference": ftext(wrong_difference),
            "gate": ftext(wrong_gate),
            "detected": bool(wrong_phase_ok),
        }

    block_construction_passed = all((
        construction_controls_passed,
        dimension_ok,
        tangent_sector_ok,
        reverse_ok,
        trace_ok,
        parseval_ok,
        wrong_phase_ok,
    ))
    check(
        f"schedule {schedule_index}: explicit-Q blocks and invariants are complete",
        block_construction_passed,
        f"all separated={all_separated}, zero candidates={len(zero_candidates)}",
    )

    record = {
        "schedule_index": schedule_index,
        "order": list(orders[schedule_index]),
        "matrix": {
            "csr_sha256": matrix_digest,
            "expected_csr_sha256": manifest_schedules[schedule_index]["csr_sha256"],
            "digest_matches": bool(matrix_ok),
            "dimension": matrix.shape[0],
            "nnz": matrix.nnz,
            "operator_error_row_bound": ftext(operator_error),
            "product_tangent_residual": ftext(null_residual),
            "product_tangent_gate": ftext(null_gate),
        },
        "geometric_edge_action": {
            "forward_order": adversary["permutation_order"](edge_map),
            "reverse_order": adversary["permutation_order"](reverse_edge_map),
            "forward_cycle_histogram": {str(key): value for key, value in cycle_histogram.items()},
            "reverse_cycle_histogram": {str(key): value for key, value in reverse_cycle_histogram.items()},
            "reverse_is_inverse": bool(inverse_relation),
        },
        "covariance": {
            "forward_power_defects": [ftext(value) for value in covariance_values],
            "reverse_power_defects": [ftext(value) for value in reverse_covariance_values],
            "forward_group_average_bound": ftext(averaging_bound),
            "reverse_group_average_bound": ftext(reverse_averaging_bound),
            "acceptance_gate": ftext(100 * operator_error),
            "corrupted_group_average_bound": ftext(corrupted_bound),
        },
        "tangent": {
            "maximum_noninvariant_projection": ftext(maximum_noninvariant_tangent),
            "gate": ftext(tangent_gate),
        },
        "block_census": {
            "independent_block_count": len(block_records),
            "weighted_dimension": sum(
                record.get("dimension", 0) * record.get("spectral_weight", 0)
                for record in block_records
            ),
            "all_eigenvalues_separated": bool(all_separated),
            "maximum_sparse_product_roundoff": ftext(maximum_roundoff),
            "maximum_lifted_residual": ftext(maximum_lifted_residual),
            "maximum_fourier_isometry_defect_all_30": ftext(maximum_isometry_defect),
            "records": block_records,
            "zero_candidates": zero_candidates,
            "reverse_sector_difference": ftext(reverse_sector_difference),
            "reverse_sector_gate": ftext(reverse_sector_gate),
            "wrong_phase_control": wrong_phase_record,
            "trace_full": ftext(trace_full),
            "trace_blocks": ftext(trace_blocks),
            "trace_difference": ftext(trace_difference),
            "trace_gate": ftext(trace_gate),
            "frobenius_full": ftext(frobenius_full),
            "frobenius_blocks": ftext(frobenius_blocks),
            "frobenius_difference": ftext(frobenius_difference),
            "frobenius_gate": ftext(frobenius_gate),
        },
        "primary_comparison": None,
        "construction_controls_passed": bool(construction_controls_passed),
        "block_construction_passed": bool(block_construction_passed),
        "all_eigenvalues_separated": bool(all_separated),
    }
    records.append(record)
    write_checkpoint(actual_hashes, geometry_record, records)

# The complete primary artifact is first parsed here, after all eleven
# adversarial spectra have been constructed or recovered from a blind-phase
# exact-provenance checkpoint.
primary = json.loads(PRIMARY.read_text())
primary_records = {
    int(record["representative_schedule"]): record
    for record in primary["representatives"]
}
all_comparisons = len(records) == len(REPRESENTATIVES)
for record in records:
    schedule_index = record["schedule_index"]
    adversarial_values = weighted_spectrum(record)
    primary_record = primary_records[schedule_index]
    primary_values = weighted_spectrum(primary_record)
    primary_blocks = primary_record["block_census"]
    adversarial_blocks = record["block_census"]
    primary_forward_bound = 100 * (
        float(primary_blocks["common_operator_and_averaging_error"])
        + float(primary_blocks["maximum_block_roundoff"])
        + float(primary_blocks["maximum_eigenpair_residual"])
    )
    adversarial_forward_bound = 100 * (
        float(record["matrix"]["operator_error_row_bound"])
        + float(record["covariance"]["forward_group_average_bound"])
        + float(adversarial_blocks["maximum_sparse_product_roundoff"])
        + float(adversarial_blocks["maximum_lifted_residual"])
    )
    comparison_gate = primary_forward_bound + adversarial_forward_bound
    if len(adversarial_values) == len(primary_values) == 19681:
        comparison_difference = float(np.max(np.abs(
            adversarial_values - primary_values
        )))
    else:
        comparison_difference = np.inf
    comparison_ok = comparison_difference <= comparison_gate
    record["primary_comparison"] = {
        "maximum_difference": ftext(comparison_difference),
        "gate": ftext(comparison_gate),
        "passes": bool(comparison_ok),
    }
    all_comparisons &= comparison_ok
    check(
        f"schedule {schedule_index}: complete adversarial spectrum reproduces primary",
        comparison_ok,
        f"difference={comparison_difference:.3e}, gate={comparison_gate:.3e}",
    )

all_construction = all((
    provenance_ok,
    manifest_ok,
    definitions_ok,
    topology_ok,
    stencil_ok,
    geometry_ok,
    negative_order_ok,
    len(records) == len(REPRESENTATIVES),
    all(record["construction_controls_passed"] for record in records),
    all(record["block_construction_passed"] for record in records),
    all_comparisons,
))
all_separated = len(records) == len(REPRESENTATIVES) and all(
    record["all_eigenvalues_separated"] for record in records
)
if not all_construction:
    outcome = "ADVERSARIAL_COXETER_CENSUS_CONSTRUCTION_INVALID"
elif all_separated:
    outcome = "ALL_24_SCHEDULE_KERNELS_ADVERSARIALLY_CORROBORATED"
else:
    outcome = "PRIMARY_COMPLETE_KERNEL_CLAIM_NOT_CORROBORATED"

verdict_ok = check(
    "the adversarial census follows the preregistered verdict hierarchy",
    outcome in {
        "ADVERSARIAL_COXETER_CENSUS_CONSTRUCTION_INVALID",
        "ALL_24_SCHEDULE_KERNELS_ADVERSARIALLY_CORROBORATED",
        "PRIMARY_COMPLETE_KERNEL_CLAIM_NOT_CORROBORATED",
    }
    and (outcome == "ADVERSARIAL_COXETER_CENSUS_CONSTRUCTION_INVALID") == (not all_construction)
    and (
        outcome == "ALL_24_SCHEDULE_KERNELS_ADVERSARIALLY_CORROBORATED"
    ) == (all_construction and all_separated),
    outcome,
)

write_checkpoint(actual_hashes, geometry_record, records, complete=True, outcome=outcome)
artifact = json.loads(OUTPUT.read_text())
artifact["global_construction"] = {
    "spatial_f_vector": list(EXPECTED_F),
    "schedule_count": len(schedules),
    "adversarial_representative_count": len(records),
    "simplex_pattern_count": len(simplex_kinds),
    "triangle_pattern_count": len(triangle_kinds),
    "stencil_diagnostics": stencil_diagnostics,
}
artifact["status_labels"] = {
    "matrix_definition": "SHARED FROZEN INPUT / NOT INDEPENDENTLY DERIVED",
    "all_schedule_kernels": (
        "DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED"
        if outcome == "ALL_24_SCHEDULE_KERNELS_ADVERSARIALLY_CORROBORATED"
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
print(f"Adversarial representatives completed: {len(records)}/11")
print(f"All construction and comparison controls passed: {all_construction}")
print(f"All adversarial eigenvalues separated: {all_separated}")
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests} passed")
print(f"Artifact: {OUTPUT.relative_to(ROOT)}")
sys.exit(0 if passed == tests and verdict_ok else 1)

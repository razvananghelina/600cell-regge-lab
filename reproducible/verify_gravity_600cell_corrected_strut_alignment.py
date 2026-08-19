#!/usr/bin/env python3
"""Compare the preregistered corrected strut carrier with frozen dynamics."""

from collections import Counter
import contextlib
import hashlib
import io
import json
import math
from pathlib import Path
import runpy

import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_corrected_strut_alignment.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_corrected_strut_alignment_protocol.md"
CORRECTED_SOURCE = HERE / "verify_gravity_600cell_corrected_strut_carrier.py"
CORRECTED_INPUT = HERE / "gravity_600cell_corrected_strut_carrier.json"
OLD_ALIGNMENT_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
OLD_ALIGNMENT_INPUT = HERE / "gravity_600cell_dust_hyperbolic_lapse_alignment.json"
TANGENT_INPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
SCHUR_INPUT = HERE / "gravity_600cell_dust_full_lapse_schur.json"
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"

PROTOCOL_COMMIT = "d019ba6"
EXPECTED_HASHES = {
    "protocol": "920c7206018fffb7fc45e180e1b264090c45b320f582926b13a9cdc48ba6270a",
    "corrected_source": "80f0a17960adee496fe7d51678ea99849280ecd3fca6254efc8acd3753aad348",
    "corrected_input": "e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7",
    "old_alignment_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "old_alignment_input": "a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff",
    "tangent": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "schur": "4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349",
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
}
INPUT_PATHS = {
    "protocol": PROTOCOL,
    "corrected_source": CORRECTED_SOURCE,
    "corrected_input": CORRECTED_INPUT,
    "old_alignment_source": OLD_ALIGNMENT_SOURCE,
    "old_alignment_input": OLD_ALIGNMENT_INPUT,
    "tangent": TANGENT_INPUT,
    "tangent_numeric": TANGENT_NUMERIC,
    "tangent_source": TANGENT_SOURCE,
    "schur": SCHUR_INPUT,
    "tick": TICK_INPUT,
}

mp.mp.dps = 100
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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sf(value):
    return f"{float(value):.17e}"


def serialize_comparison(record):
    return {
        "distance": sf(record["distance"]),
        "minimum_overlap": sf(record["minimum_overlap"]),
        "angle_degrees": sf(record["angle_degrees"]),
        "epsilon_step": sf(record["epsilon_step"]),
        "epsilon_binary": sf(record["epsilon_binary"]),
        "epsilon_ball": sf(record["epsilon_ball"]),
        "epsilon_carrier": sf(record.get("epsilon_carrier", 0.0)),
        "epsilon_distance": sf(record["epsilon_distance"]),
        "label": record["label"],
        "variant_distances": {
            key: sf(value) for key, value in record["variant_distances"].items()
        },
    }


print("=" * 78)
print("CONFIRMATORY CORRECTED-STRUT DYNAMIC ALIGNMENT")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUT_PATHS.items()}
corrected_input = json.loads(CORRECTED_INPUT.read_text())
old_input_before = json.loads(OLD_ALIGNMENT_INPUT.read_text())
tangent_input = json.loads(TANGENT_INPUT.read_text())
schur_input = json.loads(SCHUR_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and corrected_input["outcome"] == "CORRECTED_STRUT_CARRIER_FROZEN"
    and corrected_input["passed"] == corrected_input["tests"] == 13
    and corrected_input["candidate_count_per_parity"] == 1
    and corrected_input["target_artifacts_loaded"] is False
    and old_input_before["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
    and old_input_before["passed"] == old_input_before["tests"] == 14
    and tangent_input["outcome"] == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and schur_input["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
)
check("all corrected and dynamic inputs retain exact frozen provenance", provenance_ok, str(hashes))

print("reconstructing the frozen old response audit as an ordering control", flush=True)
old_stdout = io.StringIO()
with contextlib.redirect_stdout(old_stdout):
    old = runpy.run_path(str(OLD_ALIGNMENT_SOURCE))
old_reproduction_ok = bool(
    old["tests"] == old["passed"] == 14
    and old["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
    and digest(OLD_ALIGNMENT_INPUT) == EXPECTED_HASHES["old_alignment_input"]
)
check("the frozen 14-control old alignment audit reproduces byte-identically", old_reproduction_ok)

models = old["models"]
tick = old["tick"]
numeric = old["numeric"]
VARIANTS = old["VARIANTS"]
weak_artifact_errors = max(
    float(corrected_input["parities"][parity]["spectrum"][key])
    for parity in ("even", "odd")
    for key in (
        "maximum_full_relative_discrepancy",
        "maximum_restricted_relative_discrepancy",
    )
)


def committed_row_map(parity):
    result = {}
    for row in corrected_input["parities"][parity]["rows"]:
        edge = tuple(row["edge"])
        result[edge] = {
            int(item["column"]): (item["role"], item["value"])
            for item in row["coefficients"]
        }
    return result


def corrected_domain_matrix(parity, index_data, weak_positions, state):
    row_map = committed_row_map(parity)
    weak_orbit_types = [30 + position for position in weak_positions]
    pole_edges = []
    for orbit_type in weak_orbit_types:
        pole_edges.extend(index_data["orbit_edges"][orbit_type])
    logical_to_column = {
        int(edge[0]): column for column, edge in enumerate(pole_edges)
    }
    if set(logical_to_column) != set(range(120)):
        raise RuntimeError("pole ordering does not cover logical vertices")
    column_position = {
        orbit_type: position
        for position, orbit_type in enumerate(tuple(range(30, 65)) + tuple(range(65, 95)))
    }

    lam = mp.exp(mp.mpf(state[0]))
    rho = index_data["rho"]
    q_diag = lam * old["L0_SQUARE"] - rho
    kappa = rho / ((lam - 1) * q_diag)
    direct_strings = {
        "source_kappa": mp.nstr(kappa, 70),
        "target_minus_lambda_kappa": mp.nstr(-lam * kappa, 70),
        "one": "1",
    }

    matrix = np.zeros((1560, 120), dtype=float)
    direct_matrix = mp.matrix(1560, 120)
    coverage = set()
    serialized_exact = True
    pole_identity = True
    for edge, global_index in index_data["edge_to_index"].items():
        kind = index_data["edge_kind"][global_index]
        if kind not in {"internal", "pole", "new"}:
            continue
        orbit_type, group = divmod(global_index, 24)
        if orbit_type not in column_position:
            continue
        row = 24 * column_position[orbit_type] + group
        edge = tuple(map(int, edge))
        coverage.add(edge)
        stored = row_map.get(edge)
        if stored is None:
            serialized_exact = False
            continue
        for logical, (role, value) in stored.items():
            serialized_exact &= value == direct_strings[role]
            column = logical_to_column[logical]
            coefficient = {
                "one": mp.mpf(1),
                "source_kappa": kappa,
                "target_minus_lambda_kappa": -lam * kappa,
            }[role]
            direct_matrix[row, column] = coefficient
            matrix[row, column] = float(coefficient)
        if kind == "pole":
            logical = edge[0]
            pole_identity &= bool(
                len(stored) == 1
                and logical in stored
                and stored[logical][0] == "one"
                and matrix[row, logical_to_column[logical]] == 1.0
            )

    coverage_ok = bool(
        len(coverage) == 1560 and coverage == set(row_map)
        and len(row_map) == 1560
    )
    old_geometric, old_coefficient = old["geometric_lapse_matrix"](
        index_data, weak_positions, state
    )
    collective_error = float(np.max(np.abs(
        matrix @ np.ones(120) - old_geometric @ np.ones(120)
    )))
    return matrix, direct_matrix, old_geometric, {
        "coverage_ok": coverage_ok,
        "serialized_formula_exact": serialized_exact,
        "pole_identity": pole_identity,
        "collective_error": collective_error,
        "kappa": kappa,
        "old_collective_coefficient": old_coefficient,
        "logical_to_column": logical_to_column,
    }


def projected_carrier(matrix, sector):
    dimension = sector["dimension"]
    basis = old["mp_to_numpy"](sector["basis"])
    full_basis = np.kron(np.eye(65), basis)
    pole_basis = np.kron(np.eye(5), basis)
    return full_basis.conj().T @ matrix @ pole_basis


def orthonormal_columns(matrix):
    q, _ = la.qr(matrix, mode="economic")
    singular = la.svd(matrix, compute_uv=False, lapack_driver="gesvd")
    return q, singular


def subspace_distance(left, right):
    q_left, singular_left = orthonormal_columns(left)
    q_right, singular_right = orthonormal_columns(right)
    overlap = la.svd(
        q_left.conj().T @ q_right,
        compute_uv=False,
        lapack_driver="gesvd",
    )
    minimum = min(1.0, max(0.0, float(np.min(overlap))))
    distance = math.sqrt(max(0.0, 1 - minimum**2))
    return distance, minimum, singular_left, singular_right


def classify_distance(distances, binary_terms, ball_terms, condition_terms):
    op = distances["operational_primary"]
    epsilon_step = (
        abs(op - distances["operational_shadow"])
        + abs(distances["validation_primary"] - distances["validation_shadow"])
        + abs(op - distances["validation_primary"])
    )
    epsilon_carrier = weak_artifact_errors
    epsilon_binary = max(binary_terms, default=0.0) + (
        10 * np.finfo(float).eps * max(condition_terms, default=1.0)
    )
    epsilon_ball = max(ball_terms, default=0.0)
    epsilon = epsilon_step + epsilon_binary + epsilon_ball + epsilon_carrier + 1e-70
    if op <= 10 * epsilon:
        label = "IDENTIFIED"
    elif op > 100 * epsilon:
        label = "SEPARATED"
    else:
        label = "NUMERICALLY_OPEN"
    return {
        "distance": op,
        "minimum_overlap": None,
        "angle_degrees": math.degrees(math.asin(min(1.0, max(0.0, op)))),
        "epsilon_step": epsilon_step,
        "epsilon_binary": epsilon_binary,
        "epsilon_ball": epsilon_ball,
        "epsilon_carrier": epsilon_carrier,
        "epsilon_distance": epsilon,
        "label": label,
        "variant_distances": distances,
    }


def canonical_comparison(variant_data, corrected):
    distances = {}
    overlaps = {}
    ball_terms = []
    conditions = []
    for name, data in variant_data.items():
        distance, overlap, singular_lift, singular_corrected = subspace_distance(
            data["lift_midpoint"], corrected
        )
        distances[name] = distance
        overlaps[name] = overlap
        ball_terms.append(
            la.norm(data["lift_radii"], "fro")
            / max(1e-300, float(singular_lift[-1]))
        )
        conditions.append(max(
            1.0,
            float(singular_lift[0] / singular_lift[-1]),
            float(singular_corrected[0] / singular_corrected[-1]),
        ))
    result = classify_distance(distances, [], ball_terms, conditions)
    result["minimum_overlap"] = overlaps["operational_primary"]
    return result


def extreme_comparison(variant_data, branch, corrected):
    distances = {}
    overlaps = {}
    binary_terms = []
    ball_terms = []
    conditions = []
    for name, data in variant_data.items():
        response = data["response_midpoint"]
        extreme = data[f"extreme_{branch}"]
        transported = response @ extreme["basis"]
        transported_direct = response @ extreme["direct_basis"]
        distance, overlap, singular_transport, singular_corrected = subspace_distance(
            transported, corrected
        )
        distances[name] = distance
        overlaps[name] = overlap
        binary_distance, _, _, _ = subspace_distance(
            transported, transported_direct
        )
        binary_terms.append(max(binary_distance, extreme["direct_distance"]))
        ball_terms.append(
            la.norm(data["response_radii"], "fro")
            / max(1e-300, float(singular_transport[-1]))
            + la.norm(data["tangent_radii"], "fro")
              * extreme["eigenvector_condition"]
              / max(1e-300, extreme["spectral_separation"])
        )
        conditions.append(max(
            1.0,
            float(singular_transport[0] / singular_transport[-1]),
            float(singular_corrected[0] / singular_corrected[-1]),
            extreme["eigenvector_condition"],
        ))
    result = classify_distance(distances, binary_terms, ball_terms, conditions)
    result["minimum_overlap"] = overlaps["operational_primary"]
    return result


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


records = {}
all_comparisons = []
all_controls = provenance_ok and old_reproduction_ok
all_gap_gates = True
all_edge_controls = True
all_response_controls = True
all_reproduction_controls = True
nonuniform_distinct = True
corruption_sector_changes = []
corruption_comparison_changes = []

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing corrected carrier and dynamic sectors", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    index_data = old["group_and_index_data"](model, state)
    geometry = old["prepare_geometry"](model, index_data)
    weak_positions = [
        position for position in range(35)
        if index_data["edge_kind"][24 * (30 + position)] == "pole"
    ]
    weak_ok = weak_positions == schur_input["parities"][parity]["weak_orbit_positions"]
    corrected, corrected_mp, old_geometric, carrier_control = corrected_domain_matrix(
        parity, index_data, weak_positions, state
    )
    edge_ok = bool(
        weak_ok and carrier_control["coverage_ok"]
        and carrier_control["serialized_formula_exact"]
        and carrier_control["pole_identity"]
        and carrier_control["collective_error"] < 1e-14
    )
    all_edge_controls &= edge_ok
    check(
        f"{parity}: edge labels reconstruct the unique corrected matrix and collective lapse",
        edge_ok,
        f"collective error={carrier_control['collective_error']:.3e}",
    )

    sectors, sector_control = old["high_precision_sector_bases"](index_data)
    basis_ok = bool(
        sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and all(
            value < mp.mpf("1e-70")
            for key, value in sector_control.items() if key.startswith("maximum_")
        )
    )
    s_value = mp.mpf(state[0])
    kind_values = {
        "old": old["L0_SQUARE"],
        "internal": mp.exp(s_value) * old["L0_SQUARE"] - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * s_value) * old["L0_SQUARE"],
    }
    pattern_cache, branch_control = old["high_precision_pattern_cache"](
        geometry["patterns"], kind_values
    )
    kernels, kernel_control = old["assemble_full_representative_kernels"](
        index_data, geometry, pattern_cache
    )
    kernel_ok = bool(
        branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-70")
    )
    check(
        f"{parity}: target response basis, branch and reality controls pass",
        basis_ok and kernel_ok,
    )

    # Frozen corruption: reverse the source/target values on the first diagonal.
    corrupt = corrected.copy()
    row_map = committed_row_map(parity)
    first_diagonal = min(
        edge for edge, roles in row_map.items()
        if len(roles) == 2 and edge[0] < 120 <= edge[1]
    )
    global_index = index_data["edge_to_index"][first_diagonal]
    orbit_type, group = divmod(global_index, 24)
    column_position = {
        orbit: position
        for position, orbit in enumerate(tuple(range(30, 65)) + tuple(range(65, 95)))
    }
    corrupt_row = 24 * column_position[orbit_type] + group
    source_logical = first_diagonal[0]
    target_logical = first_diagonal[1] - 120
    source_column = carrier_control["logical_to_column"][source_logical]
    target_column = carrier_control["logical_to_column"][target_logical]
    corrupt[corrupt_row, source_column], corrupt[corrupt_row, target_column] = (
        corrupt[corrupt_row, target_column], corrupt[corrupt_row, source_column]
    )

    used_stored = set()
    sector_records = []
    determinant_ok = True
    ranks_ok = True
    reproductions = []
    selections_ok = True
    parity_gaps_ok = True

    for sector_index, sector in enumerate(sectors):
        dimension = sector["dimension"]
        count = 5 * dimension
        print(f"[{parity}] sector {sector_index + 1}/7 d={dimension}", flush=True)
        corrected_sector = projected_carrier(corrected, sector)
        corrupt_sector = projected_carrier(corrupt, sector)
        old_sector = old["projected_geometric_lapse"](old_geometric, sector)
        _, corrected_singular = orthonormal_columns(corrected_sector)
        ranks_ok &= corrected_singular[-1] > 1e-12

        old_new_distance, _, _, _ = subspace_distance(old_sector, corrected_sector)
        if sector["constant_overlap"] < 1 - mp.mpf("1e-70"):
            nonuniform_distinct &= old_new_distance > 1e-10
        corrupt_distance, _, _, _ = subspace_distance(
            corrected_sector, corrupt_sector
        )
        corruption_sector_changes.append(corrupt_distance)

        blocks = {
            name: old["project_full_kernel"](kernel, sector)
            for name, kernel in kernels.items()
        }
        variant_data = {}
        for name, block in blocks.items():
            response = old["response_and_lift_ball"](
                block, dimension, weak_positions
            )
            determinant_ok &= not response["det_j"].contains(0)
            tangent = numeric[
                f"{parity}_sector{sector_index}_{name}_tangent_midpoint"
            ]
            tangent_radii = numeric[
                f"{parity}_sector{sector_index}_{name}_tangent_radii"
            ]
            extreme_plus = old["extreme_subspace"](tangent, count, "plus")
            extreme_minus = old["extreme_subspace"](tangent, count, "minus")
            selections_ok &= bool(
                extreme_plus["selected_count"] == count
                and extreme_minus["selected_count"] == count
            )
            parity_gaps_ok &= bool(
                extreme_plus["gap"] > 2 and extreme_minus["gap"] > 2
            )
            variant_data[name] = {
                **response,
                "tangent_radii": tangent_radii,
                "extreme_plus": extreme_plus,
                "extreme_minus": extreme_minus,
            }

        stored_index, stored_sector = stored_sector_match(
            sector, schur_input["parities"][parity]["sectors"], used_stored
        )
        old_distance, _, _, _ = subspace_distance(
            variant_data["operational_primary"]["lift_midpoint"], old_sector
        )
        stored_distance = float(
            stored_sector["subspaces"]["canonical_vs_geometric"]["projector_distance"]
        )
        reproductions.append(abs(old_distance - stored_distance))

        comparisons = {
            "corrected_vs_canonical": canonical_comparison(
                variant_data, corrected_sector
            ),
            "corrected_vs_plus": extreme_comparison(
                variant_data, "plus", corrected_sector
            ),
            "corrected_vs_minus": extreme_comparison(
                variant_data, "minus", corrected_sector
            ),
        }
        corrupt_canonical = canonical_comparison(variant_data, corrupt_sector)
        corruption_comparison_changes.append(abs(
            comparisons["corrected_vs_canonical"]["distance"]
            - corrupt_canonical["distance"]
        ))
        for name, comparison in comparisons.items():
            all_comparisons.append({
                "parity": parity,
                "sector_index": sector_index,
                "dimension": dimension,
                "comparison": name,
                **comparison,
            })
        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "stored_sector_index": stored_index,
            "constant_overlap": sf(sector["constant_overlap"]),
            "corrected_rank": int(np.linalg.matrix_rank(corrected_sector)),
            "old_corrected_projector_distance": sf(old_new_distance),
            "corrupted_projector_distance": sf(corrupt_distance),
            "minimum_plus_gap": sf(min(
                data["extreme_plus"]["gap"] for data in variant_data.values()
            )),
            "minimum_minus_gap": sf(min(
                data["extreme_minus"]["gap"] for data in variant_data.values()
            )),
            "comparisons": {
                name: serialize_comparison(record)
                for name, record in comparisons.items()
            },
        })

    reproduction_ok = len(used_stored) == 7 and max(reproductions) < 2e-8
    response_ok = bool(
        basis_ok and kernel_ok and determinant_ok and ranks_ok
        and selections_ok
    )
    all_response_controls &= response_ok
    all_reproduction_controls &= reproduction_ok
    all_gap_gates &= parity_gaps_ok
    check(
        f"{parity}: all corrected ranks, Flint determinants and response selections pass",
        response_ok,
    )
    check(
        f"{parity}: the old target ordering reproduces every committed control distance",
        reproduction_ok,
        f"maximum error={max(reproductions):.3e}",
    )
    check(
        f"{parity}: every fixed-count extreme selection receives its frozen gap verdict",
        selections_ok,
        f"all gaps >2={parity_gaps_ok}",
    )
    all_controls &= bool(edge_ok and response_ok and reproduction_ok)
    records[parity] = {
        "edge_controls": carrier_control,
        "all_gap_gates_pass": parity_gaps_ok,
        "weak_positions": weak_positions,
        "corrupted_edge": list(first_diagonal),
        "sectors": sector_records,
    }

label_counts = Counter(item["label"] for item in all_comparisons)
comparison_counts = {}
for name in (
    "corrected_vs_canonical",
    "corrected_vs_plus",
    "corrected_vs_minus",
):
    comparison_counts[name] = dict(Counter(
        item["label"] for item in all_comparisons if item["comparison"] == name
    ))

ledger_ok = len(all_comparisons) == 42 and sum(label_counts.values()) == 42
check("the preregistered look-elsewhere ledger contains exactly 42 comparisons", ledger_ok, str(comparison_counts))

distinct_ok = nonuniform_distinct
check("the corrected and old geometric carriers are distinct in every non-uniform sector", distinct_ok)

corruption_ok = bool(
    max(corruption_sector_changes) > 1e-10
    and max(corruption_comparison_changes) > 1e-10
)
check(
    "the frozen source-target corruption changes a sector and a target comparison",
    corruption_ok,
    f"max sector/comparison changes={max(corruption_sector_changes):.3e}/"
    f"{max(corruption_comparison_changes):.3e}",
)
all_controls &= bool(ledger_ok and distinct_ok and corruption_ok)


def all_label(comparison, label):
    selected = [
        item for item in all_comparisons if item["comparison"] == comparison
    ]
    return len(selected) == 14 and all(item["label"] == label for item in selected)


if not all_controls:
    outcome = "CORRECTED_STRUT_ALIGNMENT_CONTROL_FAILED"
elif not all_gap_gates:
    outcome = "CORRECTED_STRUT_EXTREME_SELECTION_OPEN"
elif all_label("corrected_vs_canonical", "IDENTIFIED") and all_label("corrected_vs_plus", "IDENTIFIED"):
    outcome = "CORRECTED_STRUT_CANONICAL_AND_PLUS_IDENTIFIED"
elif all_label("corrected_vs_canonical", "IDENTIFIED") and all_label("corrected_vs_minus", "IDENTIFIED"):
    outcome = "CORRECTED_STRUT_CANONICAL_AND_MINUS_IDENTIFIED"
elif all_label("corrected_vs_canonical", "IDENTIFIED"):
    outcome = "CORRECTED_STRUT_CANONICAL_IDENTIFIED_ONLY"
elif label_counts == Counter({"SEPARATED": 42}):
    outcome = "CORRECTED_STRUT_ALIGNMENT_REFUTED"
else:
    outcome = "CORRECTED_STRUT_ALIGNMENT_MIXED_OR_OPEN"

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "corrected_carrier_commit": "dab941b",
    "input_sha256": hashes,
    "candidate_count_per_parity": 1,
    "comparison_count": len(all_comparisons),
    "comparison_counts": comparison_counts,
    "label_counts": dict(label_counts),
    "all_extreme_gap_gates_pass": all_gap_gates,
    "parities": records,
    "corruption": {
        "maximum_sector_projector_change": sf(max(corruption_sector_changes)),
        "maximum_canonical_comparison_change": sf(max(corruption_comparison_changes)),
        "rejected": corruption_ok,
    },
    "classification": {
        "carrier_provenance": "TARGET-BLIND FORMULA, TARGET-DRIVEN HYPOTHESIS",
        "subspace_result": "CONFIRMATORY",
        "gauge_interpretation": "OPEN",
        "curvature_response": "OPEN",
        "physical_instability": "OPEN",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")

print("-" * 78)
print(outcome)
print("labels", dict(label_counts))
print(f"TOTAL: {passed}/{tests} tests PASSED")
if outcome not in {
    "CORRECTED_STRUT_EXTREME_SELECTION_OPEN",
    "CORRECTED_STRUT_CANONICAL_AND_PLUS_IDENTIFIED",
    "CORRECTED_STRUT_CANONICAL_AND_MINUS_IDENTIFIED",
    "CORRECTED_STRUT_CANONICAL_IDENTIFIED_ONLY",
    "CORRECTED_STRUT_ALIGNMENT_REFUTED",
    "CORRECTED_STRUT_ALIGNMENT_MIXED_OR_OPEN",
} or passed != tests:
    raise SystemExit(1)


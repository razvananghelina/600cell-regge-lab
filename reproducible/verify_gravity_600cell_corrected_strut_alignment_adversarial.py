#!/usr/bin/env python3
"""Mechanically different audit of corrected-strut/dynamic separation.

The primary comparison used QR, overlap singular values and ordered Schur
vectors.  This audit uses Hermitian Gram polar factors, eigenvalues of explicit
projector differences and direct tangent eigenvectors.
"""

from collections import Counter
import contextlib
import hashlib
import io
import json
import math
from pathlib import Path
import runpy
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_corrected_strut_alignment_adversarial_protocol.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_corrected_strut_alignment.py"
PRIMARY_INPUT = HERE / "gravity_600cell_corrected_strut_alignment.json"
CORRECTED_INPUT = HERE / "gravity_600cell_corrected_strut_carrier.json"
OLD_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
OLD_INPUT = HERE / "gravity_600cell_dust_hyperbolic_lapse_alignment.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
OUTPUT = HERE / "gravity_600cell_corrected_strut_alignment_adversarial.json"

PROTOCOL_COMMIT = "19ea7d3"
PRIMARY_ARTIFACT_COMMIT = "7ef7a7b"
EXPECTED_HASHES = {
    "protocol": "8496e24125d300e9fccf3741e625ded94812ebd109f8010d7b459d5941dbc882",
    "primary_source": "d79f39380e4480aa2599d6ad0d6f56dc599268f510fe2ded1f59b9b585fb2b70",
    "primary_input": "5652b1371563ff11919be130af15f5b48850e2cc65a50ec35e5de85fdb587f90",
    "corrected_input": "e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7",
    "old_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "old_input": "a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
}
INPUTS = {
    "protocol": PROTOCOL,
    "primary_source": PRIMARY_SOURCE,
    "primary_input": PRIMARY_INPUT,
    "corrected_input": CORRECTED_INPUT,
    "old_source": OLD_SOURCE,
    "old_input": OLD_INPUT,
    "tangent_numeric": TANGENT_NUMERIC,
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


def polar_basis(matrix):
    """Orthonormalize through the Hermitian Gram eigendecomposition only."""
    gram = matrix.conj().T @ matrix
    gram = (gram + gram.conj().T) / 2
    values, vectors = la.eigh(gram, driver="evr")
    floor = float(values[-1]) * 1e-13
    if values[0] <= floor:
        raise RuntimeError(
            f"polar Gram matrix is not safely positive: {values[0]} <= {floor}"
        )
    inverse_root = (vectors / np.sqrt(values)[None, :]) @ vectors.conj().T
    basis = matrix @ inverse_root
    residual = float(la.norm(
        basis.conj().T @ basis - np.eye(matrix.shape[1]), ord=2
    ))
    condition = math.sqrt(float(values[-1] / values[0]))
    return basis, residual, condition, values


def projector_distance(left, right):
    q_left, residual_left, condition_left, _ = polar_basis(left)
    q_right, residual_right, condition_right, _ = polar_basis(right)
    difference = q_left @ q_left.conj().T - q_right @ q_right.conj().T
    difference = (difference + difference.conj().T) / 2
    eigenvalues = la.eigvalsh(difference, driver="evr")
    distance = float(np.max(np.abs(eigenvalues)))
    return {
        "distance": distance,
        "orthonormality_residual": max(residual_left, residual_right),
        "condition": max(condition_left, condition_right),
        "q_left": q_left,
        "q_right": q_right,
    }


def direct_extreme(matrix, count, branch):
    values, vectors = la.eig(matrix)
    moduli = np.abs(values)
    if branch == "plus":
        order = np.argsort(moduli)[::-1]
        selected = order[:count]
        gap = float(moduli[order[count - 1]] / moduli[order[count]])
    else:
        order = np.argsort(moduli)
        selected = order[:count]
        gap = float(moduli[order[count]] / moduli[order[count - 1]])
    basis = vectors[:, selected]
    residual = float(la.norm(
        matrix @ basis - basis * values[selected][None, :], ord="fro"
    ) / max(1.0, la.norm(matrix, ord="fro") * la.norm(basis, ord="fro")))
    return basis, gap, residual


def frozen_corrected_matrix(parity, index_data, weak_positions, committed):
    """Read exact committed rows without calling the primary verifier."""
    orbit_order = tuple(range(30, 65)) + tuple(range(65, 95))
    orbit_position = {orbit: position for position, orbit in enumerate(orbit_order)}
    weak_types = [30 + position for position in weak_positions]
    pole_edges = [
        tuple(map(int, edge))
        for orbit in weak_types for edge in index_data["orbit_edges"][orbit]
    ]
    logical_to_column = {
        edge[0]: column for column, edge in enumerate(pole_edges)
    }
    if set(logical_to_column) != set(range(120)):
        raise RuntimeError("frozen pole ordering does not cover 120 vertices")

    matrix = np.zeros((1560, 120), dtype=float)
    row_edges = set()
    for record in committed["parities"][parity]["rows"]:
        edge = tuple(map(int, record["edge"]))
        row_edges.add(edge)
        global_index = index_data["edge_to_index"][edge]
        orbit, group = divmod(global_index, 24)
        row = 24 * orbit_position[orbit] + group
        for item in record["coefficients"]:
            logical = int(item["column"])
            matrix[row, logical_to_column[logical]] = float(mp.mpf(item["value"]))

    covered = {
        tuple(map(int, edge))
        for edge, global_index in index_data["edge_to_index"].items()
        if index_data["edge_kind"][global_index] in {"internal", "pole", "new"}
    }
    return matrix, logical_to_column, bool(row_edges == covered and len(covered) == 1560)


def projected(matrix, sector, old):
    dimension = sector["dimension"]
    basis = old["mp_to_numpy"](sector["basis"])
    return (
        np.kron(np.eye(65), basis).conj().T
        @ matrix
        @ np.kron(np.eye(5), basis)
    )


def positive_and_negative_controls(candidate):
    count = candidate.shape[1]
    triangular = np.eye(count, dtype=complex)
    triangular[np.arange(count - 1), np.arange(1, count)] = 0.01
    positive = projector_distance(candidate, candidate @ triangular)["distance"]

    q, _, _, _ = polar_basis(candidate)
    projector = q @ q.conj().T
    complement = None
    complement_norm = -1.0
    for index in range(candidate.shape[0]):
        vector = np.zeros(candidate.shape[0], dtype=complex)
        vector[index] = 1
        vector -= projector @ vector
        norm = float(la.norm(vector))
        if norm > complement_norm:
            complement = vector
            complement_norm = norm
    complement /= complement_norm
    corrupted = candidate.copy().astype(complex)
    corrupted[:, 0] += 0.05 * complement
    negative = projector_distance(candidate, corrupted)["distance"]
    return positive, negative


print("=" * 78)
print("ADVERSARIAL CORRECTED-STRUT SEPARATION AUDIT")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
primary = json.loads(PRIMARY_INPUT.read_text())
committed = json.loads(CORRECTED_INPUT.read_text())
old_before = json.loads(OLD_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and primary["outcome"] == "CORRECTED_STRUT_EXTREME_SELECTION_OPEN"
    and primary["label_counts"] == {"SEPARATED": 42}
    and primary["passed"] == primary["tests"] == 15
    and committed["outcome"] == "CORRECTED_STRUT_CARRIER_FROZEN"
    and old_before["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
)
check("all post-result adversarial inputs have exact frozen provenance", provenance_ok)

print("rebuilding the frozen response operators", flush=True)
captured = io.StringIO()
original_exit = sys.exit


def audited_exit(code=0):
    if code not in (None, 0):
        raise SystemExit(code)


try:
    sys.exit = audited_exit
    with contextlib.redirect_stdout(captured):
        old = runpy.run_path(str(OLD_SOURCE))
finally:
    sys.exit = original_exit
old_reproduction_ok = bool(
    old["tests"] == old["passed"] == 14
    and old["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
    and digest(OLD_INPUT) == EXPECTED_HASHES["old_input"]
)
check("the frozen response audit still reproduces byte-identically", old_reproduction_ok)

records = {}
all_distances = []
all_primary_errors = []
all_conjugation_errors = []
all_orthogonality = []
all_conditions = []
all_eigen_residuals = []
all_positive_controls = []
all_negative_controls = []
all_corrupt_sector_changes = []
all_corrupt_target_changes = []
all_nonlast_gaps = []
all_last_gaps = []
all_geometry_ok = True

models = old["models"]
tick = old["tick"]
numeric = old["numeric"]

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing independent polar comparisons", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    index_data = old["group_and_index_data"](model, state)
    geometry = old["prepare_geometry"](model, index_data)
    weak_positions = [
        position for position in range(35)
        if index_data["edge_kind"][24 * (30 + position)] == "pole"
    ]
    corrected, logical_to_column, coverage_ok = frozen_corrected_matrix(
        parity, index_data, weak_positions, committed
    )
    sectors, sector_control = old["high_precision_sector_bases"](index_data)
    sector_bases_ok = bool(
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
    geometry_ok = bool(
        coverage_ok and sector_bases_ok and branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-70")
    )
    all_geometry_ok &= geometry_ok
    check(f"{parity}: independent carrier and response geometry reconstruct", geometry_ok)

    # Frozen role-reversal corruption on the first lower-to-upper diagonal.
    row_records = {
        tuple(map(int, record["edge"])): record
        for record in committed["parities"][parity]["rows"]
    }
    first_diagonal = min(
        edge for edge, record in row_records.items()
        if len(record["coefficients"]) == 2 and edge[0] < 120 <= edge[1]
    )
    corrupt = corrected.copy()
    orbit_order = tuple(range(30, 65)) + tuple(range(65, 95))
    orbit_position = {orbit: position for position, orbit in enumerate(orbit_order)}
    global_index = index_data["edge_to_index"][first_diagonal]
    orbit, group = divmod(global_index, 24)
    row = 24 * orbit_position[orbit] + group
    source_column = logical_to_column[first_diagonal[0]]
    target_column = logical_to_column[first_diagonal[1] - 120]
    corrupt[row, source_column], corrupt[row, target_column] = (
        corrupt[row, target_column], corrupt[row, source_column]
    )

    parity_records = []
    for sector_index, sector in enumerate(sectors):
        dimension = sector["dimension"]
        count = 5 * dimension
        print(f"[{parity}] sector {sector_index + 1}/7 d={dimension}", flush=True)
        corrected_sector = projected(corrected, sector, old)
        corrupt_sector = projected(corrupt, sector, old)
        positive, negative = positive_and_negative_controls(corrected_sector)
        all_positive_controls.append(positive)
        all_negative_controls.append(negative)
        corrupt_sector_distance = projector_distance(
            corrected_sector, corrupt_sector
        )["distance"]
        all_corrupt_sector_changes.append(corrupt_sector_distance)

        block = old["project_full_kernel"](
            kernels["operational_primary"], sector
        )
        response = old["response_and_lift_ball"](
            block, dimension, weak_positions
        )
        tangent = numeric[
            f"{parity}_sector{sector_index}_operational_primary_tangent_midpoint"
        ]
        plus_basis, plus_gap, plus_residual = direct_extreme(
            tangent, count, "plus"
        )
        minus_basis, minus_gap, minus_residual = direct_extreme(
            tangent, count, "minus"
        )
        all_eigen_residuals.extend((plus_residual, minus_residual))
        if sector_index == 6:
            all_last_gaps.extend((plus_gap, minus_gap))
        else:
            all_nonlast_gaps.extend((plus_gap, minus_gap))

        targets = {
            "corrected_vs_canonical": response["lift_midpoint"],
            "corrected_vs_plus": response["response_midpoint"] @ plus_basis,
            "corrected_vs_minus": response["response_midpoint"] @ minus_basis,
        }
        comparison_records = {}
        for name, target in targets.items():
            audit = projector_distance(corrected_sector, target)
            conjugate = projector_distance(
                corrected_sector.conj(), target.conj()
            )["distance"]
            primary_distance = float(
                primary["parities"][parity]["sectors"][sector_index]
                ["comparisons"][name]["distance"]
            )
            primary_error = abs(audit["distance"] - primary_distance)
            conjugation_error = abs(audit["distance"] - conjugate)
            all_distances.append(audit["distance"])
            all_primary_errors.append(primary_error)
            all_conjugation_errors.append(conjugation_error)
            all_orthogonality.append(audit["orthonormality_residual"])
            all_conditions.append(audit["condition"])
            comparison_records[name] = {
                "polar_projector_distance": sf(audit["distance"]),
                "primary_distance": sf(primary_distance),
                "primary_difference": sf(primary_error),
                "conjugation_difference": sf(conjugation_error),
                "maximum_orthonormality_residual": sf(
                    audit["orthonormality_residual"]
                ),
                "maximum_condition": sf(audit["condition"]),
            }

        corrupt_target_change = abs(
            projector_distance(corrected_sector, targets["corrected_vs_canonical"])["distance"]
            - projector_distance(corrupt_sector, targets["corrected_vs_canonical"])["distance"]
        )
        all_corrupt_target_changes.append(corrupt_target_change)
        parity_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "plus_gap": sf(plus_gap),
            "minus_gap": sf(minus_gap),
            "plus_eigen_residual": sf(plus_residual),
            "minus_eigen_residual": sf(minus_residual),
            "positive_control_distance": sf(positive),
            "negative_control_distance": sf(negative),
            "corrupted_sector_distance": sf(corrupt_sector_distance),
            "corrupted_target_distance_change": sf(corrupt_target_change),
            "comparisons": comparison_records,
        })
    records[parity] = parity_records

    check(
        f"{parity}: polar positive and negative controls discriminate spans",
        max(all_positive_controls[-7:]) < 2e-11
        and min(all_negative_controls[-7:]) > 1e-3,
        f"positive max={max(all_positive_controls[-7:]):.3e}, "
        f"negative min={min(all_negative_controls[-7:]):.3e}",
    )
    check(
        f"{parity}: all 21 direct/polar comparisons remain above 0.99",
        min(all_distances[-21:]) > 0.99,
        f"distance range={min(all_distances[-21:]):.9f}..{max(all_distances[-21:]):.9f}",
    )
    check(
        f"{parity}: all 21 distances reproduce the primary result",
        max(all_primary_errors[-21:]) < 2e-6,
        f"maximum difference={max(all_primary_errors[-21:]):.3e}",
    )
    check(
        f"{parity}: gap and direct-eigenvector controls retain the frozen verdict",
        min(all_nonlast_gaps[-24:]) > 2
        and max(all_last_gaps[-2:]) < 2
        and max(all_eigen_residuals[-14:]) < 2e-12,
        f"nonlast min={min(all_nonlast_gaps[-24:]):.6f}, "
        f"last max={max(all_last_gaps[-2:]):.6f}, "
        f"eigen residual={max(all_eigen_residuals[-14:]):.3e}",
    )

even_odd_errors = []
for sector_index in range(7):
    for name in (
        "corrected_vs_canonical", "corrected_vs_plus", "corrected_vs_minus"
    ):
        even = float(records["even"][sector_index]["comparisons"][name]["polar_projector_distance"])
        odd = float(records["odd"][sector_index]["comparisons"][name]["polar_projector_distance"])
        even_odd_errors.append(abs(even - odd))

conventions_ok = bool(
    max(all_conjugation_errors) < 2e-7
    and max(even_odd_errors) < 2e-7
    and max(all_orthogonality) < 2e-11
)
check(
    "complex-conjugation and staircase conventions preserve all distances",
    conventions_ok,
    f"conjugation={max(all_conjugation_errors):.3e}, "
    f"even/odd={max(even_odd_errors):.3e}, "
    f"orthonormality={max(all_orthogonality):.3e}",
)

corruption_ok = bool(
    max(all_corrupt_sector_changes) > 1e-10
    and max(all_corrupt_target_changes) > 1e-10
)
check(
    "the frozen source-target corruption changes a sector and a target distance",
    corruption_ok,
    f"sector={max(all_corrupt_sector_changes):.3e}, "
    f"target={max(all_corrupt_target_changes):.3e}",
)

all_controls = bool(
    provenance_ok and old_reproduction_ok and all_geometry_ok
    and max(all_positive_controls) < 2e-11
    and min(all_negative_controls) > 1e-3
    and min(all_distances) > 0.99
    and max(all_primary_errors) < 2e-6
    and min(all_nonlast_gaps) > 2 and max(all_last_gaps) < 2
    and max(all_eigen_residuals) < 2e-12
    and conventions_ok and corruption_ok
)
if all_controls:
    outcome = "CORRECTED_STRUT_SEPARATION_ADVERSARIALLY_CORROBORATED"
elif provenance_ok and old_reproduction_ok and all_geometry_ok:
    outcome = "CORRECTED_STRUT_SEPARATION_ADVERSARIAL_DISAGREEMENT"
else:
    outcome = "CORRECTED_STRUT_SEPARATION_ADVERSARIALLY_OPEN"

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "primary_artifact_commit": PRIMARY_ARTIFACT_COMMIT,
    "input_sha256": hashes,
    "method": "HERMITIAN_GRAM_POLAR_PLUS_PROJECTOR_EIGENVALUES_AND_DIRECT_TANGENT_EIGENVECTORS",
    "comparison_count": len(all_distances),
    "distance_range": [sf(min(all_distances)), sf(max(all_distances))],
    "maximum_primary_difference": sf(max(all_primary_errors)),
    "maximum_conjugation_difference": sf(max(all_conjugation_errors)),
    "maximum_even_odd_difference": sf(max(even_odd_errors)),
    "maximum_orthonormality_residual": sf(max(all_orthogonality)),
    "maximum_direct_eigenvector_residual": sf(max(all_eigen_residuals)),
    "condition_range": [sf(min(all_conditions)), sf(max(all_conditions))],
    "positive_control_maximum": sf(max(all_positive_controls)),
    "negative_control_minimum": sf(min(all_negative_controls)),
    "corruption": {
        "maximum_sector_change": sf(max(all_corrupt_sector_changes)),
        "maximum_target_distance_change": sf(max(all_corrupt_target_changes)),
    },
    "nonlast_gap_minimum": sf(min(all_nonlast_gaps)),
    "last_gap_maximum": sf(max(all_last_gaps)),
    "primary_outcome_preserved": "CORRECTED_STRUT_EXTREME_SELECTION_OPEN",
    "classification": {
        "subspace_separation": "DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED" if all_controls else "OPEN",
        "complete_extreme_branch_selection": "OPEN",
        "gauge_interpretation": "OPEN",
        "physical_instability": "OPEN",
    },
    "parities": records,
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"distance range {min(all_distances):.9f}..{max(all_distances):.9f}")
print(f"maximum primary difference {max(all_primary_errors):.3e}")
print(f"TOTAL: {passed}/{tests} tests PASSED")
if outcome != "CORRECTED_STRUT_SEPARATION_ADVERSARIALLY_CORROBORATED" or passed != tests:
    raise SystemExit(1)

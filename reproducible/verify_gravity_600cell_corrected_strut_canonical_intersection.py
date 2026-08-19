#!/usr/bin/env python3
"""Calibrated census of corrected-strut/canonical-graph intersections."""

from collections import Counter
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIOR = ROOT / "docs/gravity/gravity_600cell_corrected_strut_canonical_intersection_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_corrected_strut_canonical_intersection_protocol.md"
CORRECTED = HERE / "gravity_600cell_corrected_strut_carrier.json"
PRIMARY_ALIGNMENT = HERE / "gravity_600cell_corrected_strut_alignment.json"
ADVERSARIAL_ALIGNMENT = HERE / "gravity_600cell_corrected_strut_alignment_adversarial.json"
OLD_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
OLD_INPUT = HERE / "gravity_600cell_dust_hyperbolic_lapse_alignment.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
OUTPUT = HERE / "gravity_600cell_corrected_strut_canonical_intersection.json"

PRIOR_COMMIT = "27b8e85"
PROTOCOL_COMMIT = "786cae7"
EXPECTED_HASHES = {
    "prior": "a56c5042be8596876ca1b2b45085069bcb1da7d1d859ac12d98f5e7e8420fc4f",
    "protocol": "f6cbe0f14f78cf04d1ebd96374389a554fc53b5a3ec79a6daca3e5875c8fc456",
    "corrected": "e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7",
    "primary_alignment": "5652b1371563ff11919be130af15f5b48850e2cc65a50ec35e5de85fdb587f90",
    "adversarial_alignment": "3b0fd6da76195279f1beac540c326c61eff5e3172a63bb89baf69502254c5b1f",
    "old_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "old_input": "a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
}
INPUTS = {
    "prior": PRIOR,
    "protocol": PROTOCOL,
    "corrected": CORRECTED,
    "primary_alignment": PRIMARY_ALIGNMENT,
    "adversarial_alignment": ADVERSARIAL_ALIGNMENT,
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


def frozen_corrected_matrix(parity, index_data, weak_positions, committed):
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
        raise RuntimeError("pole ordering does not cover the 120 vertices")

    matrix = np.zeros((1560, 120), dtype=float)
    stored_edges = set()
    for record in committed["parities"][parity]["rows"]:
        edge = tuple(map(int, record["edge"]))
        stored_edges.add(edge)
        global_index = index_data["edge_to_index"][edge]
        orbit, group = divmod(global_index, 24)
        row = 24 * orbit_position[orbit] + group
        for item in record["coefficients"]:
            logical = int(item["column"])
            matrix[row, logical_to_column[logical]] = float(mp.mpf(item["value"]))

    expected_edges = {
        tuple(map(int, edge))
        for edge, global_index in index_data["edge_to_index"].items()
        if index_data["edge_kind"][global_index] in {"internal", "pole", "new"}
    }
    coverage_ok = stored_edges == expected_edges and len(expected_edges) == 1560
    return matrix, logical_to_column, coverage_ok


def projected(matrix, sector, old):
    dimension = sector["dimension"]
    basis = old["mp_to_numpy"](sector["basis"])
    return (
        np.kron(np.eye(65), basis).conj().T
        @ matrix
        @ np.kron(np.eye(5), basis)
    )


def singular_values(matrix):
    return la.svd(
        matrix, compute_uv=False, lapack_driver="gesdd", check_finite=True
    )


def classify_singulars(values, epsilon):
    labels = []
    for value in values:
        if value <= 10 * epsilon:
            labels.append("ZERO")
        elif value > 100 * epsilon:
            labels.append("NONZERO")
        else:
            labels.append("OPEN")
    return labels


def kernel_projector(matrix, nullity):
    if nullity == 0:
        return np.zeros((matrix.shape[1], matrix.shape[1]), dtype=complex)
    _, _, vh = la.svd(matrix, full_matrices=False, lapack_driver="gesdd")
    basis = vh[-nullity:, :].conj().T
    return basis @ basis.conj().T


def serialize_complex_matrix(matrix):
    return [
        [
            {"real": sf(value.real), "imaginary": sf(value.imag)}
            for value in row
        ]
        for row in matrix
    ]


print("=" * 78)
print("CORRECTED-STRUT / CANONICAL-GRAPH INTERSECTION CENSUS")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
committed = json.loads(CORRECTED.read_text())
primary = json.loads(PRIMARY_ALIGNMENT.read_text())
adversarial = json.loads(ADVERSARIAL_ALIGNMENT.read_text())
old_before = json.loads(OLD_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and committed["outcome"] == "CORRECTED_STRUT_CARRIER_FROZEN"
    and primary["outcome"] == "CORRECTED_STRUT_EXTREME_SELECTION_OPEN"
    and primary["label_counts"] == {"SEPARATED": 42}
    and adversarial["outcome"]
    == "CORRECTED_STRUT_SEPARATION_ADVERSARIALLY_CORROBORATED"
    and old_before["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
)
check("all intersection-census inputs have exact frozen provenance", provenance_ok)

print("reconstructing frozen response graphs", flush=True)
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
old_ok = bool(
    old["tests"] == old["passed"] == 14
    and old["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
    and digest(OLD_INPUT) == EXPECTED_HASHES["old_input"]
)
check("the frozen response audit reproduces byte-identically", old_ok)

models = old["models"]
tick = old["tick"]
variants = old["VARIANTS"]
carrier_relative_error = max(
    float(committed["parities"][parity]["spectrum"][key])
    for parity in ("even", "odd")
    for key in (
        "maximum_full_relative_discrepancy",
        "maximum_restricted_relative_discrepancy",
    )
)

records = {}
all_resolved = True
all_graph_controls = True
all_synthetic_controls = True
all_image_controls = True
all_basis_controls = True
all_conjugation_errors = []
all_corruption_matrix_changes = []
all_corruption_singular_changes = []
all_nullities = []

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing 7 graph differences", flush=True)
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
        coverage_ok
        and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and all(
            value < mp.mpf("1e-70")
            for key, value in sector_control.items() if key.startswith("maximum_")
        )
        and branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-70")
    )
    check(f"{parity}: carrier, sectors and response geometry reconstruct", geometry_ok)
    all_graph_controls &= geometry_ok

    # Frozen source/target-role corruption.
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

    sector_records = []
    parity_resolved = True
    parity_graph_ok = True
    parity_synthetic_ok = True
    parity_image_ok = True
    parity_basis_ok = True
    for sector_index, sector in enumerate(sectors):
        dimension = sector["dimension"]
        count = 5 * dimension
        print(f"[{parity}] sector {sector_index + 1}/7 d={dimension}", flush=True)
        g_matrix = projected(corrected, sector, old)
        corrupt_matrix = projected(corrupt, sector, old)
        weak_rows = [
            position * dimension + component
            for position in weak_positions for component in range(dimension)
        ]
        pole_g_error = float(la.norm(
            g_matrix[weak_rows, :] - np.eye(count), ord=np.inf
        ))

        variant_data = {}
        for name in variants:
            block = old["project_full_kernel"](kernels[name], sector)
            response = old["response_and_lift_ball"](
                block, dimension, weak_positions
            )
            variant_data[name] = response
        canonical = {
            name: data["lift_midpoint"] for name, data in variant_data.items()
        }
        pole_c_errors = [
            float(la.norm(matrix[weak_rows, :] - np.eye(count), ord=np.inf))
            for matrix in canonical.values()
        ]
        ranks_ok = bool(
            np.linalg.matrix_rank(g_matrix) == count
            and all(np.linalg.matrix_rank(matrix) == count for matrix in canonical.values())
        )
        graph_ok = bool(
            ranks_ok and pole_g_error < 1e-13 and max(pole_c_errors) < 1e-13
        )
        parity_graph_ok &= graph_ok

        differences = {
            name: g_matrix - matrix for name, matrix in canonical.items()
        }
        operational = differences["operational_primary"]
        lift_radius = max(
            float(la.norm(data["lift_radii"], ord="fro"))
            for data in variant_data.values()
        )
        step_error = max(
            float(la.norm(matrix - operational, ord=2))
            for matrix in differences.values()
        )
        carrier_error = carrier_relative_error * float(la.norm(g_matrix, ord=2))
        condition = max(
            float(np.linalg.cond(g_matrix)),
            *(float(np.linalg.cond(matrix)) for matrix in canonical.values()),
        )
        matrix_norm = max(
            float(la.norm(g_matrix, ord=2)),
            *(float(la.norm(matrix, ord=2)) for matrix in canonical.values()),
        )
        binary_error = 50 * np.finfo(float).eps * matrix_norm * condition
        epsilon = lift_radius + step_error + carrier_error + binary_error + 1e-70

        variant_singulars = {
            name: singular_values(matrix) for name, matrix in differences.items()
        }
        variant_labels = {
            name: classify_singulars(values, epsilon)
            for name, values in variant_singulars.items()
        }
        variant_nullities = {
            name: labels.count("ZERO") for name, labels in variant_labels.items()
        }
        actual_labels = variant_labels["operational_primary"]
        resolved = bool(
            all("OPEN" not in labels for labels in variant_labels.values())
            and len(set(variant_nullities.values())) == 1
        )
        nullity = variant_nullities["operational_primary"] if resolved else None
        parity_resolved &= resolved
        if nullity is not None:
            all_nullities.append(nullity)

        # Synthetic controls under the exact same absolute threshold.
        positive_labels = classify_singulars(
            singular_values(np.zeros_like(operational)), epsilon
        )
        injection = np.zeros_like(operational)
        nonpole_rows = [index for index in range(operational.shape[0]) if index not in weak_rows]
        injection[nonpole_rows[:count], :] = np.eye(count)
        negative_labels = classify_singulars(singular_values(injection), epsilon)
        synthetic_ok = bool(
            positive_labels.count("ZERO") == count
            and negative_labels.count("NONZERO") == count
        )
        parity_synthetic_ok &= synthetic_ok

        # Same coefficient change on both graphs: D -> D T.
        transform = np.eye(count, dtype=complex)
        transform[np.arange(count - 1), np.arange(1, count)] = 0.01
        transformed_labels = classify_singulars(
            singular_values(operational @ transform), epsilon
        )
        transformed_nullity = transformed_labels.count("ZERO")
        basis_ok = bool(
            "OPEN" not in transformed_labels
            and nullity is not None and transformed_nullity == nullity
        )
        parity_basis_ok &= basis_ok

        # Independent image-intersection formula.
        joined = np.hstack((g_matrix, canonical["operational_primary"]))
        joined_labels = classify_singulars(singular_values(joined), epsilon)
        joined_rank = joined_labels.count("NONZERO")
        image_intersection = 2 * count - joined_rank
        image_ok = bool(
            "OPEN" not in joined_labels
            and nullity is not None and image_intersection == nullity
        )
        parity_image_ok &= image_ok

        conjugate_error = float(np.max(np.abs(
            singular_values(operational) - singular_values(operational.conj())
        )))
        all_conjugation_errors.append(conjugate_error)
        corruption_matrix_change = float(la.norm(
            (corrupt_matrix - canonical["operational_primary"]) - operational,
            ord=2,
        ))
        corruption_singular_change = float(np.max(np.abs(
            singular_values(corrupt_matrix - canonical["operational_primary"])
            - variant_singulars["operational_primary"]
        )))
        all_corruption_matrix_changes.append(corruption_matrix_change)
        all_corruption_singular_changes.append(corruption_singular_change)

        projector = (
            kernel_projector(operational, nullity)
            if nullity is not None else None
        )
        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "coefficient_dimension": count,
            "pole_identity_errors": {
                "corrected": sf(pole_g_error),
                "canonical_maximum": sf(max(pole_c_errors)),
            },
            "matrix_uncertainty": {
                "lift_ball": sf(lift_radius),
                "step": sf(step_error),
                "carrier": sf(carrier_error),
                "binary_condition": sf(binary_error),
                "total": sf(epsilon),
            },
            "operational_singular_values": [
                sf(value) for value in variant_singulars["operational_primary"]
            ],
            "operational_labels": actual_labels,
            "variant_nullities": variant_nullities,
            "resolved": resolved,
            "nullity": nullity,
            "smallest_singular_over_uncertainty": sf(
                variant_singulars["operational_primary"][-1] / epsilon
            ),
            "kernel_projector": (
                serialize_complex_matrix(projector) if projector is not None else None
            ),
            "image_formula_intersection": image_intersection,
            "basis_transform_nullity": transformed_nullity,
            "positive_control_nullity": positive_labels.count("ZERO"),
            "negative_control_nullity": negative_labels.count("ZERO"),
            "conjugation_singular_error": sf(conjugate_error),
            "corruption_matrix_change": sf(corruption_matrix_change),
            "corruption_singular_change": sf(corruption_singular_change),
        })

    all_resolved &= parity_resolved
    all_graph_controls &= parity_graph_ok
    all_synthetic_controls &= parity_synthetic_ok
    all_image_controls &= parity_image_ok
    all_basis_controls &= parity_basis_ok
    check(
        f"{parity}: all literal graph and full-rank controls pass",
        parity_graph_ok,
    )
    check(
        f"{parity}: all seven nullities receive calibrated verdicts",
        parity_resolved,
        str([record["nullity"] for record in sector_records]),
    )
    check(
        f"{parity}: synthetic and coefficient-basis controls pass",
        parity_synthetic_ok and parity_basis_ok,
    )
    check(
        f"{parity}: stacked-image formula reproduces every nullity",
        parity_image_ok,
    )
    records[parity] = sector_records

conjugation_ok = max(all_conjugation_errors) < 1e-13
check(
    "complex conjugation preserves all singular spectra",
    conjugation_ok,
    f"maximum error={max(all_conjugation_errors):.3e}",
)
corruption_ok = bool(
    max(all_corruption_matrix_changes) > carrier_relative_error
    and max(all_corruption_singular_changes) > carrier_relative_error
)
check(
    "the frozen source-target corruption changes a matrix and singular value",
    corruption_ok,
    f"matrix={max(all_corruption_matrix_changes):.3e}, "
    f"singular={max(all_corruption_singular_changes):.3e}",
)

hard_controls = bool(
    provenance_ok and old_ok and all_graph_controls
    and all_synthetic_controls and all_image_controls and all_basis_controls
    and conjugation_ok and corruption_ok
)
if not hard_controls:
    outcome = "CORRECTED_STRUT_CANONICAL_INTERSECTION_CONTROL_FAILED"
elif not all_resolved:
    outcome = "CORRECTED_STRUT_CANONICAL_INTERSECTION_NUMERICALLY_OPEN"
else:
    outcome = "CORRECTED_STRUT_CANONICAL_INTERSECTION_RESOLVED"

payload = {
    "prior_commit": PRIOR_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "sector_count": 14,
    "nullities": all_nullities,
    "nullity_counts": dict(Counter(all_nullities)),
    "all_resolved": all_resolved,
    "maximum_conjugation_singular_error": sf(max(all_conjugation_errors)),
    "corruption": {
        "maximum_matrix_change": sf(max(all_corruption_matrix_changes)),
        "maximum_singular_change": sf(max(all_corruption_singular_changes)),
    },
    "classification": {
        "intersection_dimensions": "DERIVED COMPUTATIONAL" if all_resolved else "OPEN",
        "gauge_interpretation": "OPEN",
        "physical_interpretation": "OPEN",
        "full_scale_plus_strut_carrier": "OPEN",
    },
    "parities": records,
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print("nullities", all_nullities)
print(f"TOTAL: {passed}/{tests} tests PASSED")
if outcome not in {
    "CORRECTED_STRUT_CANONICAL_INTERSECTION_RESOLVED",
    "CORRECTED_STRUT_CANONICAL_INTERSECTION_NUMERICALLY_OPEN",
} or passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""QR/Frobenius adversarial audit of zero corrected-strut intersection."""

import ast
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
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_corrected_strut_canonical_intersection_adversarial_protocol.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_corrected_strut_canonical_intersection.py"
PRIMARY_INPUT = HERE / "gravity_600cell_corrected_strut_canonical_intersection.json"
CORRECTED = HERE / "gravity_600cell_corrected_strut_carrier.json"
OLD_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
OLD_INPUT = HERE / "gravity_600cell_dust_hyperbolic_lapse_alignment.json"
OUTPUT = HERE / "gravity_600cell_corrected_strut_canonical_intersection_adversarial.json"

PROTOCOL_COMMIT = "2afb0c7"
PRIMARY_ARTIFACT_COMMIT = "b64fd83"
EXPECTED_HASHES = {
    "protocol": "7356d0bfa01a42c0c373a07f092d79cba80e918ec6c2880296b40a1522e98ce2",
    "primary_source": "efedf6c360766168b019ecd94c24fe4159bd23ffb07c71574e0bb34bfd8d7af7",
    "primary_input": "422d8d8cb0fc0d72d842e3bf79609d4d985da6237c58e7c699b5f9cc21b65cec",
    "corrected": "e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7",
    "old_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "old_input": "a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff",
}
INPUTS = {
    "protocol": PROTOCOL,
    "primary_source": PRIMARY_SOURCE,
    "primary_input": PRIMARY_INPUT,
    "corrected": CORRECTED,
    "old_source": OLD_SOURCE,
    "old_input": OLD_INPUT,
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


def load_nondeterminative_helpers():
    wanted = {"frozen_corrected_matrix", "projected"}
    tree = ast.parse(PRIMARY_SOURCE.read_text(), filename=str(PRIMARY_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    if {node.name for node in body} != wanted:
        raise RuntimeError("primary geometry helper set changed")
    namespace = {"np": np, "mp": mp}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(PRIMARY_SOURCE), "exec"), namespace)
    return namespace["frozen_corrected_matrix"], namespace["projected"]


def qr_lower_bound(matrix):
    """No SVD/eigenvalue call: lower bound from pivoted QR and ||R^-1||_F."""
    rows, columns = matrix.shape
    if columns == 0:
        return {
            "lower": float("inf"), "residual": 0.0,
            "unitarity": 0.0, "inverse_frobenius": 0.0,
        }
    q, r, pivots = la.qr(
        matrix, mode="economic", pivoting=True, check_finite=True
    )
    scale = max(1.0, float(la.norm(matrix, ord="fro")))
    residual = float(la.norm(matrix[:, pivots] - q @ r, ord="fro") / scale)
    unitarity = float(la.norm(
        q.conj().T @ q - np.eye(columns), ord="fro"
    ))
    diagonal = np.abs(np.diag(r))
    if len(diagonal) != columns or np.min(diagonal) == 0:
        lower = 0.0
        inverse_frobenius = float("inf")
    else:
        try:
            inverse = la.solve_triangular(
                r, np.eye(columns, dtype=complex), lower=False,
                check_finite=True,
            )
            inverse_frobenius = float(la.norm(inverse, ord="fro"))
            lower = 1.0 / inverse_frobenius
        except la.LinAlgError:
            lower = 0.0
            inverse_frobenius = float("inf")
    return {
        "lower": lower,
        "residual": residual,
        "unitarity": unitarity,
        "inverse_frobenius": inverse_frobenius,
        "minimum_r_diagonal": float(np.min(diagonal)) if len(diagonal) else 0.0,
        "rows": rows,
        "columns": columns,
    }


def uncertainty(g_matrix, canonical, variant_data, carrier_relative_error):
    differences = {
        name: g_matrix - matrix for name, matrix in canonical.items()
    }
    operational = differences["operational_primary"]
    lift_radius = max(
        float(la.norm(data["lift_radii"], ord="fro"))
        for data in variant_data.values()
    )
    step = max(
        float(la.norm(matrix - operational, ord=2))
        for matrix in differences.values()
    )
    carrier = carrier_relative_error * float(la.norm(g_matrix, ord=2))
    # Condition surrogates come from QR rather than singular values.
    graph_qr = [qr_lower_bound(g_matrix)] + [
        qr_lower_bound(matrix) for matrix in canonical.values()
    ]
    graph_norms = [float(la.norm(g_matrix, ord="fro"))] + [
        float(la.norm(matrix, ord="fro")) for matrix in canonical.values()
    ]
    condition_surrogate = max(
        norm * record["inverse_frobenius"]
        for norm, record in zip(graph_norms, graph_qr)
    )
    binary = float(
        50 * np.finfo(float).eps * max(graph_norms) * condition_surrogate
    )
    return {
        "lift_ball": lift_radius,
        "step": step,
        "carrier": carrier,
        "binary_condition": binary,
        "total": float(lift_radius + step + carrier + binary + 1e-70),
        "condition_surrogate": condition_surrogate,
        "graph_qr": graph_qr,
        "differences": differences,
    }


print("=" * 78)
print("ADVERSARIAL ZERO CORRECTED-STRUT INTERSECTION AUDIT")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
primary = json.loads(PRIMARY_INPUT.read_text())
committed = json.loads(CORRECTED.read_text())
old_before = json.loads(OLD_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and primary["outcome"] == "CORRECTED_STRUT_CANONICAL_INTERSECTION_RESOLVED"
    and primary["nullity_counts"] == {"0": 14}
    and primary["passed"] == primary["tests"] == 14
    and committed["outcome"] == "CORRECTED_STRUT_CARRIER_FROZEN"
    and old_before["outcome"] == "HYPERBOLIC_EXTREME_SUBSPACE_OPEN"
)
check("all disclosed adversarial inputs have exact frozen provenance", provenance_ok)

frozen_corrected_matrix, projected = load_nondeterminative_helpers()

print("reconstructing response graphs for QR audit", flush=True)
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

carrier_relative_error = max(
    float(committed["parities"][parity]["spectrum"][key])
    for parity in ("even", "odd")
    for key in (
        "maximum_full_relative_discrepancy",
        "maximum_restricted_relative_discrepancy",
    )
)
variants = old["VARIANTS"]
records = {}
all_certified = True
all_geometry = True
all_qr_residuals = []
all_unitarity = []
all_robust_ratios = []
all_primary_consistency = []
all_positive_failures = []
all_negative_passes = []
all_transform_passes = []
all_row_errors = []
all_conjugation_errors = []
all_corruption_changes = []

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing 7 QR certificates", flush=True)
    model = old["models"][parity]
    state = old["tick"]["solutions"][parity]["state"]
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
    all_geometry &= geometry_ok
    check(f"{parity}: carrier, sectors and response geometry reconstruct", geometry_ok)

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
    parity_certified = True
    parity_controls = True
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

        variant_data = {}
        for name in variants:
            block = old["project_full_kernel"](kernels[name], sector)
            variant_data[name] = old["response_and_lift_ball"](
                block, dimension, weak_positions
            )
        canonical = {
            name: data["lift_midpoint"] for name, data in variant_data.items()
        }
        pole_error = max(
            float(la.norm(g_matrix[weak_rows, :] - np.eye(count), ord=np.inf)),
            *(float(la.norm(matrix[weak_rows, :] - np.eye(count), ord=np.inf))
              for matrix in canonical.values()),
        )
        error = uncertainty(
            g_matrix, canonical, variant_data, carrier_relative_error
        )
        epsilon = error["total"]

        graph_certificates = error["graph_qr"]
        graphs_full = all(
            record["lower"] - epsilon > 100 * epsilon
            for record in graph_certificates
        )
        variant_certificates = {
            name: qr_lower_bound(matrix)
            for name, matrix in error["differences"].items()
        }
        robust_lowers = {
            name: record["lower"] - epsilon
            for name, record in variant_certificates.items()
        }
        certified = bool(
            pole_error < 1e-13 and graphs_full
            and all(value > 100 * epsilon for value in robust_lowers.values())
        )
        parity_certified &= certified
        all_robust_ratios.extend(
            value / epsilon for value in robust_lowers.values()
        )
        all_qr_residuals.extend(
            record["residual"] for record in variant_certificates.values()
        )
        all_unitarity.extend(
            record["unitarity"] for record in variant_certificates.values()
        )

        # One-sided consistency only; the stored SVD is not used to certify.
        primary_sigma = float(
            primary["parities"][parity][sector_index]
            ["operational_singular_values"][-1]
        )
        primary_consistency = bool(
            variant_certificates["operational_primary"]["lower"]
            <= primary_sigma * (1 + 2e-7) + epsilon
        )
        all_primary_consistency.append(primary_consistency)

        zero_certificate = qr_lower_bound(np.zeros_like(g_matrix))
        positive_fails = zero_certificate["lower"] == 0.0
        injection = np.zeros_like(g_matrix)
        nonpole_rows = [index for index in range(g_matrix.shape[0]) if index not in weak_rows]
        injection[nonpole_rows[:count], :] = np.eye(count)
        negative_certificate = qr_lower_bound(injection)
        negative_passes = bool(
            negative_certificate["lower"] - epsilon > 100 * epsilon
        )

        transform = np.eye(count, dtype=complex)
        transform[np.arange(count - 1), np.arange(1, count)] = 0.01
        transformed = qr_lower_bound(
            error["differences"]["operational_primary"] @ transform
        )
        transform_passes = bool(
            transformed["lower"] - epsilon > 100 * epsilon
        )

        operational = error["differences"]["operational_primary"]
        phases = np.where(np.arange(operational.shape[0]) % 2 == 0, 1.0, -1.0)
        row_changed = operational[::-1, :] * phases[:, None]
        row_lower = qr_lower_bound(row_changed)["lower"]
        raw_lower = variant_certificates["operational_primary"]["lower"]
        row_error = abs(row_lower - raw_lower) / max(raw_lower, 1e-300)
        conjugate_lower = qr_lower_bound(operational.conj())["lower"]
        conjugation_error = abs(conjugate_lower - raw_lower) / max(raw_lower, 1e-300)

        corrupted_lower = qr_lower_bound(
            corrupt_matrix - canonical["operational_primary"]
        )["lower"]
        corruption_change = abs(corrupted_lower - raw_lower)

        controls = bool(
            primary_consistency and positive_fails and negative_passes
            and transform_passes and row_error < 2e-7
            and conjugation_error < 2e-7
            and max(record["residual"] for record in variant_certificates.values()) < 2e-12
            and max(record["unitarity"] for record in variant_certificates.values()) < 2e-12
        )
        parity_controls &= controls
        all_positive_failures.append(positive_fails)
        all_negative_passes.append(negative_passes)
        all_transform_passes.append(transform_passes)
        all_row_errors.append(row_error)
        all_conjugation_errors.append(conjugation_error)
        all_corruption_changes.append(corruption_change)

        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "coefficient_dimension": count,
            "pole_identity_error": sf(pole_error),
            "matrix_uncertainty": {
                key: sf(value) for key, value in error.items()
                if key in {"lift_ball", "step", "carrier", "binary_condition", "total", "condition_surrogate"}
            },
            "variant_raw_lower_bounds": {
                name: sf(record["lower"])
                for name, record in variant_certificates.items()
            },
            "variant_robust_lower_bounds": {
                name: sf(value) for name, value in robust_lowers.items()
            },
            "minimum_robust_lower_over_uncertainty": sf(
                min(robust_lowers.values()) / epsilon
            ),
            "maximum_qr_residual": sf(
                max(record["residual"] for record in variant_certificates.values())
            ),
            "maximum_q_unitarity_residual": sf(
                max(record["unitarity"] for record in variant_certificates.values())
            ),
            "primary_smallest_singular": sf(primary_sigma),
            "primary_one_sided_consistency": primary_consistency,
            "positive_intersection_control_rejected": positive_fails,
            "negative_intersection_control_accepted": negative_passes,
            "basis_transform_accepted": transform_passes,
            "row_convention_relative_error": sf(row_error),
            "conjugation_relative_error": sf(conjugation_error),
            "corruption_lower_bound_change": sf(corruption_change),
            "certified_full_column_rank": certified,
        })

    all_certified &= parity_certified
    check(
        f"{parity}: all seven QR/Frobenius lower bounds certify zero intersection",
        parity_certified,
        f"minimum ratio={min(float(record['minimum_robust_lower_over_uncertainty']) for record in sector_records):.3e}",
    )
    check(
        f"{parity}: synthetic, basis and convention controls pass",
        parity_controls,
    )
    records[parity] = sector_records

corruption_ok = max(all_corruption_changes) > carrier_relative_error
check(
    "the frozen role corruption changes an adversarial QR lower bound",
    corruption_ok,
    f"maximum change={max(all_corruption_changes):.3e}",
)

all_controls = bool(
    provenance_ok and old_ok and all_geometry
    and all(all_positive_failures) and all(all_negative_passes)
    and all(all_transform_passes) and all(all_primary_consistency)
    and max(all_row_errors) < 2e-7
    and max(all_conjugation_errors) < 2e-7
    and max(all_qr_residuals) < 2e-12
    and max(all_unitarity) < 2e-12
    and corruption_ok
)
if not all_controls:
    outcome = "CORRECTED_STRUT_ZERO_INTERSECTION_ADVERSARIALLY_OPEN"
elif all_certified:
    outcome = "CORRECTED_STRUT_ZERO_INTERSECTION_ADVERSARIALLY_CORROBORATED"
else:
    outcome = "CORRECTED_STRUT_ZERO_INTERSECTION_ADVERSARIAL_DISAGREEMENT"

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "primary_artifact_commit": PRIMARY_ARTIFACT_COMMIT,
    "input_sha256": hashes,
    "method": "COLUMN_PIVOTED_QR_AND_INVERSE_R_FROBENIUS_LOWER_BOUND",
    "sector_count": 14,
    "minimum_robust_lower_over_uncertainty": sf(min(all_robust_ratios)),
    "maximum_qr_residual": sf(max(all_qr_residuals)),
    "maximum_q_unitarity_residual": sf(max(all_unitarity)),
    "maximum_row_convention_relative_error": sf(max(all_row_errors)),
    "maximum_conjugation_relative_error": sf(max(all_conjugation_errors)),
    "maximum_corruption_lower_bound_change": sf(max(all_corruption_changes)),
    "classification": {
        "pure_strut_canonical_intersection": (
            "DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED"
            if all_certified and all_controls else "OPEN"
        ),
        "scale_plus_strut_intersection": "OPEN",
        "gauge_interpretation": "OPEN",
        "physical_interpretation": "OPEN",
    },
    "parities": records,
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"minimum robust lower/uncertainty {min(all_robust_ratios):.3e}")
print(f"TOTAL: {passed}/{tests} tests PASSED")
if outcome != "CORRECTED_STRUT_ZERO_INTERSECTION_ADVERSARIALLY_CORROBORATED" or passed != tests:
    raise SystemExit(1)

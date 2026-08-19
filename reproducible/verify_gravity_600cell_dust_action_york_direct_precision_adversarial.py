#!/usr/bin/env python3
"""Mechanically independent audit of the pseudo-longitudinal residual."""

import ast
from collections import Counter
import contextlib
from hashlib import sha256
import importlib.util
import io
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_dust_action_york_direct_precision_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_action_york_direct_precision_adversarial_protocol.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_dust_action_york_direct_precision.py"
PRIMARY_JSON = HERE / "gravity_600cell_dust_action_york_direct_precision.json"
PRIMARY_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_action_york_direct_precision_protocol.md"
CENTERED_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
RIGIDITY_SOURCE = HERE / "verify_gravity_600cell_dust_rigidity_york.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
COMMONS_SOURCE = ROOT / "commons/cell600.py"

PROTOCOL_COMMIT = "c65ed2a"
EXPECTED_HASHES = {
    "protocol": "a11ff07846ff62335f9b029af1b680007e68b3d165e8aedb114d71dd907ce0a9",
    "primary_source": "73d852d58b21a9a15306a565d5cf4fb998b159fadb82830739ab0996ac07270e",
    "primary_json": "d57351e852ab40eb7809397c84e5f57ff58e5ae0bd31f9dcaf87efdc84be76b5",
    "primary_protocol": "e655d6025e790ff2beb653c5e9f4c2f38233606c3607f9219ae222bafdfed36e",
    "centered_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "rigidity_source": "deba8d9f9bca4a5848134943ec77544e5487d44a59c44234f632b6f2aeb51382",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
PARITIES = ("even", "odd")
VARIANTS = (
    "operational_primary", "operational_shadow",
    "validation_primary", "validation_shadow",
)
SELECTED_SECTORS = (4, 5)
MACHINE_EPSILON = np.finfo(float).eps
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
    return sha256(path.read_bytes()).hexdigest()


def sf(value):
    return f"{float(value):.17e}"


def operator_norm(matrix):
    singular = la.svdvals(matrix)
    return float(singular[0]) if len(singular) else 0.0


def zero_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "OPEN"
    if value <= 10 * error:
        return "ZERO_CONSISTENT"
    if value > 100 * error:
        return "NONZERO_RESOLVED"
    return "OPEN"


def load_functions(source, wanted):
    tree = ast.parse(source.read_text(), filename=str(source))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"missing audited helpers in {source.name}: {wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(source), "exec"),
        globals(),
    )


load_functions(CONFORMAL_SOURCE, {
    "mp_frobenius", "mp_submatrix", "cluster_sorted", "orbit_sort_key",
    "edge_image", "group_data", "incidence_data", "mp_to_numpy",
})
load_functions(FULL_SOURCE, {"high_precision_sector_bases"})
load_functions(RIGIDITY_SOURCE, {"build_rigidity"})


def calibrated_rank(matrix):
    singular = la.svdvals(matrix)
    error = float(
        1000 * MACHINE_EPSILON * max(matrix.shape)
        * max(1.0, float(singular[0]) if len(singular) else 1.0)
    )
    rank = int(np.sum(singular > 100 * error))
    open_count = int(np.sum((singular >= 10 * error) & (singular <= 100 * error)))
    return rank, open_count, error


def adversarial_cell(m_matrix, v_matrix, c_matrix, d_matrix):
    h_m = (m_matrix + m_matrix.conj().T) / 2
    h_v = (v_matrix + v_matrix.conj().T) / 2

    q_c, r_c, piv_c = la.qr(c_matrix, mode="economic", pivoting=True)
    diagonal_c = np.abs(np.diag(r_c))
    u_basis = q_c[:, :5]
    row = u_basis.conj().T @ h_m
    q_row, r_row = la.qr(row.conj().T, mode="full")
    w_basis = q_row[:, 5:]
    q_sum = np.column_stack((u_basis, w_basis))
    longitudinal_raw = la.solve(q_sum, d_matrix)[5:, :]
    q_l, r_l, piv_l = la.qr(
        longitudinal_raw, mode="economic", pivoting=True
    )
    diagonal_l = np.abs(np.diag(r_l))
    l_basis = q_l[:, :15]

    m_s = w_basis.conj().T @ h_m @ w_basis
    v_s = w_basis.conj().T @ h_v @ w_basis
    b_matrix = -(m_s + m_s.conj().T) / 2
    a_matrix = -(v_s + v_s.conj().T) / 2
    bl = b_matrix @ l_basis
    al = a_matrix @ l_basis
    q_bl, r_bl, piv_bl = la.qr(bl, mode="economic", pivoting=True)
    span_residual = operator_norm(
        (np.eye(25) - q_bl @ q_bl.conj().T) @ al
    )
    generalized = la.solve(b_matrix, a_matrix)
    commutator_residual = operator_norm(
        (np.eye(25) - l_basis @ l_basis.conj().T)
        @ generalized @ l_basis
    )
    augmented = la.svdvals(np.column_stack((bl, al)))

    singular_bl = la.svdvals(bl)
    kappa = max(
        float(diagonal_l[0] / diagonal_l[14]),
        float(np.linalg.cond(b_matrix)),
        float(singular_bl[0] / singular_bl[-1]),
    )
    scale = max(
        1.0, operator_norm(a_matrix), operator_norm(b_matrix),
        operator_norm(d_matrix),
    )
    arithmetic_floor = float(1000 * MACHINE_EPSILON * 30 * kappa * scale)
    augmented_rank = int(np.sum(augmented > 100 * arithmetic_floor))
    return {
        "rank_C": int(np.sum(diagonal_c > 1e-6)),
        "rank_shape_row": int(np.sum(np.abs(np.diag(r_row)) > 1e-6)),
        "rank_longitudinal": int(np.sum(diagonal_l > 1e-6)),
        "minimum_longitudinal_QR_diagonal": float(diagonal_l[14]),
        "minimum_B_eigenvalue": float(la.eigvalsh(b_matrix)[0]),
        "kappa": kappa,
        "arithmetic_floor": arithmetic_floor,
        "span_residual": span_residual,
        "commutator_residual": commutator_residual,
        "augmented_rank": augmented_rank,
        "sixteenth_augmented_singular": float(augmented[15]),
    }


paths = {
    "protocol": PROTOCOL,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "primary_protocol": PRIMARY_PROTOCOL,
    "centered_json": CENTERED_JSON,
    "centered_npz": CENTERED_NPZ,
    "rigidity_source": RIGIDITY_SOURCE,
    "conformal_source": CONFORMAL_SOURCE,
    "full_source": FULL_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "commons": COMMONS_SOURCE,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all adversarial inputs have exact frozen provenance",
      provenance_ok, str(hashes))
primary = json.loads(PRIMARY_JSON.read_text())
centered_json = json.loads(CENTERED_JSON.read_text())
primary_ok = bool(
    primary["outcome"] == "DIRECT_LONGITUDINAL_IDENTITY_NUMERICALLY_REFUTED"
    and primary["passed"] == primary["tests"] == 17
    and primary["label_counts"]["cross"] == {"NONZERO_RESOLVED": 16}
    and primary["label_counts"]["image"] == {"NONZERO_RESOLVED": 16}
    and centered_json["outcome"] == "CENTERED_JACOBI_CERTIFIED"
)
check("the primary refutation and independent centered source retain their outcomes",
      primary_ok)

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_action_york_adversarial", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check("the independently imported slab geometry retains 43/43 controls",
      gro.tests == gro.passed == 43)

groups = {parity: group_data(gro.models[parity], gro) for parity in PARITIES}
incidences = {parity: incidence_data(groups[parity]) for parity in PARITIES}
vertices, adjacency, _ = build_600cell()
vertices = vertices / np.linalg.norm(vertices, axis=1)[:, None]
rigidity = {}
global_controls = {}
global_ok = True
for parity in PARITIES:
    r_matrix, radial, tangent, lengths, length_square = build_rigidity(
        vertices, groups[parity]["edge_order"]
    )
    d_matrix = r_matrix @ tangent
    c_matrix = incidences[parity]["incidence"].astype(float)
    rank_r, open_r, _ = calibrated_rank(r_matrix)
    rank_d, open_d, _ = calibrated_rank(d_matrix)
    rank_cd, open_cd, _ = calibrated_rank(np.column_stack((c_matrix, d_matrix)))
    intersection = 120 + rank_d - rank_cd
    ok = bool(
        rank_r == 470 and rank_d == 354 and rank_cd == 470
        and intersection == 4 and open_r == open_d == open_cd == 0
        and np.array_equal(adjacency.astype(np.int8),
                           incidences[parity]["adjacency"])
    )
    global_ok &= ok
    rigidity[parity] = {"C": c_matrix, "D": d_matrix}
    global_controls[parity] = {
        "rank_R": rank_r, "rank_D": rank_d, "rank_C_plus_D": rank_cd,
        "intersection_C_D": intersection,
    }
check("the binary rigidity route independently recovers ranks 470/354/4",
      global_ok, str(global_controls))

sector_data, sector_controls = high_precision_sector_bases(groups["even"])
basis_ok = bool(
    tuple(sector["dimension"] for sector in sector_data) == (3, 2, 2, 2, 1, 1, 1)
    and max(value for key, value in sector_controls.items()
            if key.startswith("maximum_")) < mp.mpf("1e-70")
)
check("the independent carrier rebuild retains all seven symmetry sectors",
      basis_ok)

geometry = {}
carrier_geometry_ok = True
for parity in PARITIES:
    for sector_index in SELECTED_SECTORS:
        basis = mp_to_numpy(sector_data[sector_index]["basis"])
        edge_basis = np.kron(np.eye(30, dtype=np.complex128), basis)
        c_matrix = edge_basis.conj().T @ rigidity[parity]["C"]
        d_matrix = edge_basis.conj().T @ rigidity[parity]["D"]
        geometry[(parity, sector_index)] = {"C": c_matrix, "D": d_matrix}
        carrier_geometry_ok &= c_matrix.shape == (30, 120) and d_matrix.shape == (30, 480)
check("all four adversarial minimal carriers have the expected literal sizes",
      carrier_geometry_ok)

archive = np.load(CENTERED_NPZ)
families = {}
carrier_ok = True
for key, carrier in geometry.items():
    parity, sector_index = key
    records = {}
    for variant in VARIANTS:
        prefix = f"{parity}_sector{sector_index}_{variant}"
        records[variant] = adversarial_cell(
            archive[f"{prefix}_M_midpoint"],
            archive[f"{prefix}_V_midpoint"],
            carrier["C"], carrier["D"],
        )
    primary_record = records["operational_primary"]
    errors = {}
    for field in ("span_residual", "commutator_residual"):
        errors[field] = float(
            max(abs(record[field] - primary_record[field])
                for record in records.values())
            + max(record["arithmetic_floor"] for record in records.values())
        )
    for record in records.values():
        record["span_label"] = zero_label(
            record["span_residual"], errors["span_residual"]
        )
        record["commutator_label"] = zero_label(
            record["commutator_residual"], errors["commutator_residual"]
        )
        carrier_ok &= bool(
            record["rank_C"] == record["rank_shape_row"] == 5
            and record["rank_longitudinal"] == 15
            and record["minimum_longitudinal_QR_diagonal"] > 1e-6
            and record["minimum_B_eigenvalue"] > 1e-3
            and math.isfinite(record["kappa"])
        )
    families[key] = {"records": records, "errors": errors}

check("QR reconstruction independently retains every 5+25 and 15+10 carrier",
      carrier_ok)
all_records = [
    record for family in families.values() for record in family["records"].values()
]
span_counts = Counter(record["span_label"] for record in all_records)
commutator_counts = Counter(record["commutator_label"] for record in all_records)
rank_counts = Counter(record["augmented_rank"] for record in all_records)
labels_complete = bool(
    len(all_records) == 16
    and sum(span_counts.values()) == sum(commutator_counts.values()) == 16
)
check("both independent adversarial residuals classify all sixteen cells",
      labels_complete,
      f"span={dict(span_counts)}, commutator={dict(commutator_counts)}, ranks={dict(rank_counts)}")

ratios = {}
for (parity, sector), family in families.items():
    adversarial = family["records"]["operational_primary"]
    primary_record = primary["cells"][f"{parity}_sector{sector}"]["variants"]["operational_primary"]
    ratios[f"{parity}_sector{sector}"] = {
        "span_to_primary_image": (
            adversarial["span_residual"] / float(primary_record["image_residual"])
        ),
        "commutator_to_primary_image": (
            adversarial["commutator_residual"] / float(primary_record["image_residual"])
        ),
    }
ratio_finite = all(
    math.isfinite(value)
    for record in ratios.values() for value in record.values()
)
check("the new residuals admit finite comparisons with the primary artifact",
      ratio_finite)

controls_ok = bool(
    provenance_ok and primary_ok and gro.tests == gro.passed == 43
    and global_ok and basis_ok and carrier_geometry_ok and carrier_ok
    and labels_complete and ratio_finite
)
confirmed = bool(
    controls_ok
    and span_counts == {"NONZERO_RESOLVED": 16}
    and commutator_counts == {"NONZERO_RESOLVED": 16}
    and all(record["augmented_rank"] > 15 for record in all_records)
)
refuted = bool(
    controls_ok and any(
        record["span_label"] == "ZERO_CONSISTENT"
        or record["commutator_label"] == "ZERO_CONSISTENT"
        for record in all_records
    )
)
if not controls_ok:
    outcome = "ADVERSARIAL_PSEUDOLONGITUDINAL_CONTROL_FAILED"
elif confirmed:
    outcome = "ADVERSARIAL_DIRECT_REFUTATION_CONFIRMED"
elif refuted:
    outcome = "ADVERSARIAL_DIRECT_REFUTATION_REFUTED"
else:
    outcome = "ADVERSARIAL_DIRECT_REFUTATION_OPEN"
allowed = {
    "ADVERSARIAL_PSEUDOLONGITUDINAL_CONTROL_FAILED",
    "ADVERSARIAL_DIRECT_REFUTATION_CONFIRMED",
    "ADVERSARIAL_DIRECT_REFUTATION_REFUTED",
    "ADVERSARIAL_DIRECT_REFUTATION_OPEN",
}
check("the frozen adversarial hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "primary_outcome": primary["outcome"],
    "global_controls": global_controls,
    "cells": {
        f"{parity}_sector{sector}": {
            "errors": {name: sf(value) for name, value in family["errors"].items()},
            "variants": {
                variant: {
                    name: (sf(value) if isinstance(value, (float, np.floating)) else value)
                    for name, value in record.items()
                }
                for variant, record in family["records"].items()
            },
            "primary_comparison_ratios": {
                name: sf(value)
                for name, value in ratios[f"{parity}_sector{sector}"].items()
            },
        }
        for (parity, sector), family in families.items()
    },
    "label_counts": {
        "span": dict(span_counts),
        "commutator": dict(commutator_counts),
        "augmented_rank": {str(key): value for key, value in rank_counts.items()},
    },
    "classification": {
        "primary_numerical_refutation": (
            "INDEPENDENTLY CONFIRMED DERIVED COMPUTATIONAL / STRUCTURAL"
            if confirmed else "REFUTED" if refuted else "OPEN"
        ),
        "formal_symbolic_or_interval_theorem": False,
        "continuum_diffeomorphism_claim": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
print("span labels:", dict(span_counts))
print("commutator labels:", dict(commutator_counts))
print("augmented ranks:", dict(rank_counts))
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

#!/usr/bin/env python3
"""Scale-invariant correction of the independent action-York audit."""

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
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_dust_action_york_scale_invariant_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_action_york_scale_invariant_adversarial_protocol.md"
OLD_SOURCE = HERE / "verify_gravity_600cell_dust_action_york_direct_precision_adversarial.py"
OLD_JSON = HERE / "gravity_600cell_dust_action_york_direct_precision_adversarial.json"
OLD_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_action_york_direct_precision_adversarial_protocol.md"
OLD_RESULT = ROOT / "docs/gravity/gravity_600cell_dust_action_york_direct_precision_adversarial_result.md"
CENTERED_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
RIGIDITY_SOURCE = HERE / "verify_gravity_600cell_dust_rigidity_york.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
COMMONS_SOURCE = ROOT / "commons/cell600.py"

PROTOCOL_COMMIT = "199696f"
EXPECTED_HASHES = {
    "protocol": "36f9373d07744712efdce875a6b0a1a4d5be7aa236cc29cc7f5255e6f77c00e9",
    "old_source": "719d46bc519e152c06fd0bc064962532ae759303fcfd69759c89dd5cd7bc8352",
    "old_json": "e39203741513f128a208f22896abef53daa12db089ee7e43abf9c90643fc579b",
    "old_protocol": "a11ff07846ff62335f9b029af1b680007e68b3d165e8aedb114d71dd907ce0a9",
    "old_result": "b9dd475651a1a840d7bfc07e32a8dab72c32f591b9db5e2614b06e667384624e",
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
    return rank, open_count


def relative_cell(m_matrix, v_matrix, c_matrix, d_matrix):
    h_m = (m_matrix + m_matrix.conj().T) / 2
    h_v = (v_matrix + v_matrix.conj().T) / 2
    q_c, r_c, _ = la.qr(c_matrix, mode="economic", pivoting=True)
    diagonal_c = np.abs(np.diag(r_c))
    u_basis = q_c[:, :5]
    row = u_basis.conj().T @ h_m
    q_row, r_row = la.qr(row.conj().T, mode="full")
    w_basis = q_row[:, 5:]
    q_sum = np.column_stack((u_basis, w_basis))
    longitudinal_raw = la.solve(q_sum, d_matrix)[5:, :]
    q_l, r_l, _ = la.qr(longitudinal_raw, mode="economic", pivoting=True)
    diagonal_l = np.abs(np.diag(r_l))
    l_basis = q_l[:, :15]

    b_matrix = -(w_basis.conj().T @ h_m @ w_basis)
    a_matrix = -(w_basis.conj().T @ h_v @ w_basis)
    b_matrix = (b_matrix + b_matrix.conj().T) / 2
    a_matrix = (a_matrix + a_matrix.conj().T) / 2
    bl = b_matrix @ l_basis
    al = a_matrix @ l_basis
    q_bl, r_bl, _ = la.qr(bl, mode="economic", pivoting=True)
    generalized = la.solve(b_matrix, a_matrix)
    gl = generalized @ l_basis

    span_residual = operator_norm(
        (np.eye(25) - q_bl @ q_bl.conj().T) @ al
    )
    commutator_residual = operator_norm(
        (np.eye(25) - l_basis @ l_basis.conj().T) @ gl
    )
    al_norm = operator_norm(al)
    gl_norm = operator_norm(gl)
    if min(al_norm, gl_norm) <= 1e-12:
        relative_span = relative_commutator = math.inf
    else:
        relative_span = span_residual / al_norm
        relative_commutator = commutator_residual / gl_norm

    singular_bl = la.svdvals(bl)
    norm_a = operator_norm(a_matrix)
    norm_g = operator_norm(generalized)
    kappa_relative = max(
        float(diagonal_l[0] / diagonal_l[14]),
        float(np.linalg.cond(b_matrix)),
        float(singular_bl[0] / singular_bl[-1]),
        float(norm_a / al_norm),
        float(norm_g / gl_norm),
    )
    relative_floor = float(1000 * MACHINE_EPSILON * 30 * kappa_relative)
    augmented = la.svdvals(np.column_stack((bl, al)))

    norm_b = operator_norm(b_matrix)
    inverse_norm_b = 1 / float(la.svdvals(b_matrix)[-1])
    inequality_floor = float(
        1000 * MACHINE_EPSILON * 30 * kappa_relative
        * max(1.0, span_residual, commutator_residual)
    )
    inequality_one = bool(
        span_residual <= norm_b * commutator_residual + inequality_floor
    )
    inequality_two = bool(
        commutator_residual <= inverse_norm_b * span_residual + inequality_floor
    )
    return {
        "rank_C": int(np.sum(diagonal_c > 1e-6)),
        "rank_shape_row": int(np.sum(np.abs(np.diag(r_row)) > 1e-6)),
        "rank_longitudinal": int(np.sum(diagonal_l > 1e-6)),
        "minimum_longitudinal_QR_diagonal": float(diagonal_l[14]),
        "minimum_B_eigenvalue": float(la.eigvalsh(b_matrix)[0]),
        "span_residual": span_residual,
        "commutator_residual": commutator_residual,
        "AL_norm": al_norm,
        "GL_norm": gl_norm,
        "relative_span": relative_span,
        "relative_commutator": relative_commutator,
        "kappa_relative": kappa_relative,
        "relative_floor": relative_floor,
        "augmented_singular": augmented,
        "inequality_one": inequality_one,
        "inequality_two": inequality_two,
    }


paths = {
    "protocol": PROTOCOL,
    "old_source": OLD_SOURCE,
    "old_json": OLD_JSON,
    "old_protocol": OLD_PROTOCOL,
    "old_result": OLD_RESULT,
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
check("all scale-invariant audit inputs have exact provenance",
      provenance_ok, str(hashes))
old = json.loads(OLD_JSON.read_text())
old_ok = bool(
    old["outcome"] == "ADVERSARIAL_DIRECT_REFUTATION_REFUTED"
    and old["passed"] == old["tests"] == 10
    and old["label_counts"]["span"] == {"NONZERO_RESOLVED": 16}
    and old["label_counts"]["commutator"] == {"ZERO_CONSISTENT": 16}
    and old["label_counts"]["augmented_rank"] == {"24": 16}
)
check("the contradictory predecessor is preserved literally",
      old_ok)

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_york_scale_invariant", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check("the literal geometry import retains all 43 controls",
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
    rank_r, open_r = calibrated_rank(r_matrix)
    rank_d, open_d = calibrated_rank(d_matrix)
    rank_cd, open_cd = calibrated_rank(np.column_stack((c_matrix, d_matrix)))
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
        "rank_R": rank_r, "rank_D": rank_d,
        "rank_C_plus_D": rank_cd, "intersection_C_D": intersection,
    }
check("the corrected audit independently retains ranks 470/354/4",
      global_ok, str(global_controls))

sector_data, sector_controls = high_precision_sector_bases(groups["even"])
basis_ok = bool(
    tuple(sector["dimension"] for sector in sector_data) == (3, 2, 2, 2, 1, 1, 1)
    and max(value for key, value in sector_controls.items()
            if key.startswith("maximum_")) < mp.mpf("1e-70")
)
check("all seven symmetry sectors are reconstructed independently",
      basis_ok)

geometry = {}
for parity in PARITIES:
    for sector_index in SELECTED_SECTORS:
        basis = mp_to_numpy(sector_data[sector_index]["basis"])
        edge_basis = np.kron(np.eye(30, dtype=np.complex128), basis)
        geometry[(parity, sector_index)] = {
            "C": edge_basis.conj().T @ rigidity[parity]["C"],
            "D": edge_basis.conj().T @ rigidity[parity]["D"],
        }

archive = np.load(CENTERED_NPZ)
families = {}
carrier_ok = denominator_ok = inequalities_ok = True
for key, carrier in geometry.items():
    parity, sector_index = key
    records = {}
    for variant in VARIANTS:
        prefix = f"{parity}_sector{sector_index}_{variant}"
        record = relative_cell(
            archive[f"{prefix}_M_midpoint"],
            archive[f"{prefix}_V_midpoint"],
            carrier["C"], carrier["D"],
        )
        records[variant] = record
        carrier_ok &= bool(
            record["rank_C"] == record["rank_shape_row"] == 5
            and record["rank_longitudinal"] == 15
            and record["minimum_longitudinal_QR_diagonal"] > 1e-6
            and record["minimum_B_eigenvalue"] > 1e-3
        )
        denominator_ok &= min(record["AL_norm"], record["GL_norm"]) > 1e-12
        inequalities_ok &= record["inequality_one"] and record["inequality_two"]

    primary = records["operational_primary"]
    errors = {}
    for field in ("relative_span", "relative_commutator"):
        errors[field] = float(
            max(abs(record[field] - primary[field]) for record in records.values())
            + max(record["relative_floor"] for record in records.values())
        )
    augmented_error = max(errors.values())
    for record in records.values():
        record["span_label"] = zero_label(
            record["relative_span"], errors["relative_span"]
        )
        record["commutator_label"] = zero_label(
            record["relative_commutator"], errors["relative_commutator"]
        )
        threshold = float(
            100 * augmented_error * record["augmented_singular"][0]
        )
        record["augmented_rank"] = int(np.sum(
            record["augmented_singular"] > threshold
        ))
        record["augmented_threshold"] = threshold
        record["sixteenth_augmented_singular"] = float(
            record["augmented_singular"][15]
        )
        del record["augmented_singular"]
    families[key] = {"records": records, "errors": errors}

check("every corrected carrier and both normalization denominators are resolved",
      carrier_ok and denominator_ok)
check("both exact conditioning inequalities hold in all sixteen cells",
      inequalities_ok)

all_records = [
    record for family in families.values() for record in family["records"].values()
]
span_counts = Counter(record["span_label"] for record in all_records)
commutator_counts = Counter(record["commutator_label"] for record in all_records)
rank_counts = Counter(record["augmented_rank"] for record in all_records)
complete = len(all_records) == 16
check("both relative residuals and augmented ranks classify all sixteen cells",
      complete,
      f"span={dict(span_counts)}, commutator={dict(commutator_counts)}, ranks={dict(rank_counts)}")

controls_ok = bool(
    provenance_ok and old_ok and gro.tests == gro.passed == 43
    and global_ok and basis_ok and carrier_ok and denominator_ok
    and inequalities_ok and complete
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
        and record["commutator_label"] == "ZERO_CONSISTENT"
        and record["augmented_rank"] == 15
        for record in all_records
    )
)
if not controls_ok:
    outcome = "SCALE_INVARIANT_ADVERSARIAL_CONTROL_FAILED"
elif confirmed:
    outcome = "SCALE_INVARIANT_DIRECT_REFUTATION_CONFIRMED"
elif refuted:
    outcome = "SCALE_INVARIANT_DIRECT_REFUTATION_REFUTED"
else:
    outcome = "SCALE_INVARIANT_DIRECT_REFUTATION_OPEN"
allowed = {
    "SCALE_INVARIANT_ADVERSARIAL_CONTROL_FAILED",
    "SCALE_INVARIANT_DIRECT_REFUTATION_CONFIRMED",
    "SCALE_INVARIANT_DIRECT_REFUTATION_REFUTED",
    "SCALE_INVARIANT_DIRECT_REFUTATION_OPEN",
}
check("the corrected frozen hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "preserved_predecessor_outcome": old["outcome"],
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
        }
        for (parity, sector), family in families.items()
    },
    "label_counts": {
        "relative_span": dict(span_counts),
        "relative_commutator": dict(commutator_counts),
        "augmented_rank": {str(key): value for key, value in rank_counts.items()},
    },
    "classification": {
        "independent_primary_refutation": (
            "CONFIRMED DERIVED COMPUTATIONAL / STRUCTURAL"
            if confirmed else "REFUTED" if refuted else "OPEN"
        ),
        "predecessor_absolute_classifier": "REFUTED BY SCALE DEPENDENCE",
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
print("relative span labels:", dict(span_counts))
print("relative commutator labels:", dict(commutator_counts))
print("augmented ranks:", dict(rank_counts))
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

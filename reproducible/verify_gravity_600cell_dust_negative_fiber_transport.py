#!/usr/bin/env python3
"""Target-disclosed old/shifted negative-fiber transport gate.

Prior-art/framing commit: 42dc0e2.
Preregistered protocol commit: 07de221.
No fitted basis, Procrustes, polar or permutation alignment is used.
"""

import ast
from collections import Counter
import contextlib
import hashlib
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
OLD_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
OLD_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
OLD_SHAPE = HERE / "gravity_600cell_dust_shape_stiffness.json"
NEW_JSON = HERE / "gravity_600cell_dust_shifted_centered.json"
NEW_NPZ = HERE / "gravity_600cell_dust_shifted_centered.npz"
NEW_SHAPE = HERE / "gravity_600cell_dust_shifted_shape_stiffness.json"
DIRECT = HERE / "gravity_600cell_dust_shifted_direct_precision.json"
TANGENT_JSON = HERE / "gravity_600cell_dust_two_step_full_tangent.json"
TANGENT_NPZ = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_two_step_full_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
OUTPUT = HERE / "gravity_600cell_dust_negative_fiber_transport.json"

PRIOR_ART_COMMIT = "42dc0e2"
PROTOCOL_COMMIT = "07de221"
EXPECTED_HASHES = {
    "old_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "old_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "old_shape": "03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868",
    "new_json": "265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47",
    "new_npz": "c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8",
    "new_shape": "14fe5bc91e3ae4712c6ea19b8120785e2facd364e1ceb194009123fa353a4315",
    "direct": "86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944",
    "tangent_json": "f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc",
    "tangent_npz": "ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d",
    "tangent_source": "c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
PARITIES = ("even", "odd")
TARGET_SECTORS = (4, 5)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sf(value):
    return f"{float(value):.17e}"


def operator_norm(matrix):
    values = la.svdvals(matrix)
    return float(values[0]) if len(values) else 0.0


def matrix_error(midpoint, radii, n):
    return float(
        la.norm(radii, "fro")
        + 1000 * MACHINE_EPSILON * n * max(1.0, operator_norm(midpoint))
    )


def restriction_error(midpoint, base_error, eta_shape, n):
    norm = operator_norm(midpoint)
    return float(
        base_error
        + 2 * eta_shape * (norm + base_error)
        + 1000 * MACHINE_EPSILON * n * max(1.0, norm)
    )


def load_audited_helpers():
    wanted = {
        "mp_frobenius", "mp_submatrix", "cluster_sorted",
        "orbit_sort_key", "edge_image", "group_data", "incidence_data",
        "mp_to_numpy", "component_reenclosure_radii",
    }
    tree = ast.parse(CONFORMAL_SOURCE.read_text(), filename=str(CONFORMAL_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited conformal helper mismatch: {wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(CONFORMAL_SOURCE), "exec"),
        globals(),
    )
    tree = ast.parse(FULL_SOURCE.read_text(), filename=str(FULL_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "high_precision_sector_bases"
    ]
    if len(body) != 1:
        raise RuntimeError("audited high-precision sector function is missing")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(FULL_SOURCE), "exec"),
        globals(),
    )


def source_matrix(archive, prefix, name):
    midpoint = np.asarray(archive[f"{prefix}_{name}_midpoint"])
    stored = np.asarray(archive[f"{prefix}_{name}_radii"])
    radii = component_reenclosure_radii(midpoint, stored)
    return midpoint, radii


def projector_class(distance, error):
    if not math.isfinite(distance) or not math.isfinite(error) or error < 0:
        return "FIBER_IDENTITY_OPEN"
    if distance <= 10 * error:
        return "COMMON_FIBER_RESOLVED"
    if distance > 100 * error:
        return "ROTATED_FIBER_RESOLVED"
    return "FIBER_IDENTITY_OPEN"


def leakage_class(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "LEAKAGE_OPEN"
    if value <= 10 * error:
        return "LEAKAGE_ZERO_CONSISTENT"
    if value > 100 * error:
        return "LEAKAGE_NONZERO_RESOLVED"
    return "LEAKAGE_OPEN"


paths = {
    "old_json": OLD_JSON,
    "old_npz": OLD_NPZ,
    "old_shape": OLD_SHAPE,
    "new_json": NEW_JSON,
    "new_npz": NEW_NPZ,
    "new_shape": NEW_SHAPE,
    "direct": DIRECT,
    "tangent_json": TANGENT_JSON,
    "tangent_npz": TANGENT_NPZ,
    "tangent_source": TANGENT_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "full_source": FULL_SOURCE,
    "conformal_source": CONFORMAL_SOURCE,
    "commons": COMMONS_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
old_json = json.loads(OLD_JSON.read_text())
old_shape = json.loads(OLD_SHAPE.read_text())
new_json = json.loads(NEW_JSON.read_text())
new_shape = json.loads(NEW_SHAPE.read_text())
direct = json.loads(DIRECT.read_text())
tangent_json = json.loads(TANGENT_JSON.read_text())
old_archive = np.load(OLD_NPZ, allow_pickle=False)
new_archive = np.load(NEW_NPZ, allow_pickle=False)
tangent_archive = np.load(TANGENT_NPZ, allow_pickle=False)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and old_json["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and old_json["numeric_archive_arrays"] == len(old_archive.files) == 560
    and old_shape["outcome"] == "SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED"
    and new_json["outcome"] == "SHIFTED_CENTERED_CERTIFIED"
    and new_json["numeric_archive_arrays"] == len(new_archive.files) == 560
    and new_shape["outcome"] == "SHIFTED_SHAPE_STIFFNESS_SIGN_OPEN"
    and direct["outcome"] == "SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS"
    and direct["target_rank_cells"] == direct["cells"] == 16
    and tangent_json["outcome"] == "TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED"
    and tangent_json["numeric_archive_arrays"] == len(tangent_archive.files) == 448
)
check("all fiber-transport inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_negative_fiber_transport", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
geometry_ok = gro.tests == gro.passed == 43
check("the literal one-slab geometry retains all 43 certificates", geometry_ok)

load_audited_helpers()
groups = {parity: group_data(gro.models[parity], gro) for parity in PARITIES}
incidences = {parity: incidence_data(groups[parity]) for parity in PARITIES}
carrier_geometry_ok = all(
    len(groups[parity]["edge_order"]) == 720
    and incidences[parity]["incidence"].shape == (720, 120)
    and incidences[parity]["equivariant"]
    and incidences[parity]["connected"]
    and incidences[parity]["numerical_rank"] == 120
    and tangent_json["parities"][parity]["boundary_mapping"] == list(range(30))
    for parity in PARITIES
)
check("both schedules retain the exact rank-120 carrier and identity seam map",
      carrier_geometry_ok)

sector_data, sector_controls = high_precision_sector_bases(groups["even"])
sector_basis_ok = bool(
    tuple(item["dimension"] for item in sector_data) == (3, 2, 2, 2, 1, 1, 1)
    and sector_controls["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
    and max(
        value for key, value in sector_controls.items()
        if key.startswith("maximum_")
    ) < mp.mpf("1e-70")
)
check("the same deterministic seven-sector basis is reconstructed", sector_basis_ok)

print("=" * 78)
print("NEGATIVE-FIBER IDENTITY AND CANONICAL-TANGENT TRANSPORT")
print("=" * 78)

projectors = {"old": {}, "shifted": {}}
projector_records = []
projector_controls = True
for parity in PARITIES:
    incidence = incidences[parity]["incidence"].astype(np.complex128)
    for sector_index in TARGET_SECTORS:
        sector = sector_data[sector_index]
        dimension = sector["dimension"]
        n = 30 * dimension
        r = 5 * dimension
        s = 25 * dimension
        basis = mp_to_numpy(sector["basis"])
        edge_basis = np.kron(np.eye(30, dtype=np.complex128), basis)
        compressed = edge_basis.conj().T @ incidence
        left, singular, _ = la.svd(compressed, full_matrices=False)
        epsilon_c = float(
            1000 * MACHINE_EPSILON * max(compressed.shape)
            * max(1.0, float(singular[0]))
        )
        rank_c = int(np.sum(singular > 100 * epsilon_c))
        gap_c = float(singular[r - 1])
        eta_k = (
            float(2 * epsilon_c / (gap_c - 2 * epsilon_c)
                  + 1000 * MACHINE_EPSILON * n)
            if gap_c > 2 * epsilon_c else math.inf
        )
        u_basis = left[:, :r]
        for variant in VARIANTS:
            key = (parity, sector_index, variant)
            for time_name, archive in (("old", old_archive), ("shifted", new_archive)):
                prefix = f"{parity}_sector{sector_index}_{variant}"
                m_mid, m_rad = source_matrix(archive, prefix, "M")
                v_mid, v_rad = source_matrix(archive, prefix, "V")
                h_m = (m_mid + m_mid.conj().T) / 2
                h_v = (v_mid + v_mid.conj().T) / 2
                epsilon_hm = matrix_error(h_m, (m_rad + m_rad.T) / 2, n)
                epsilon_hv = matrix_error(h_v, (v_rad + v_rad.T) / 2, n)

                row = u_basis.conj().T @ h_m
                _, singular_row, right_row = la.svd(row, full_matrices=True)
                epsilon_row = float(
                    epsilon_hm
                    + 2 * eta_k * (operator_norm(h_m) + epsilon_hm)
                    + 1000 * MACHINE_EPSILON * n
                    * max(1.0, operator_norm(row))
                )
                row_rank = int(np.sum(singular_row > 100 * epsilon_row))
                row_zero = int(np.sum(singular_row < 10 * epsilon_row))
                row_open = len(singular_row) - row_rank - row_zero
                row_gap = float(singular_row[r - 1])
                eta_s = (
                    float(2 * epsilon_row / (row_gap - 2 * epsilon_row)
                          + 1000 * MACHINE_EPSILON * n)
                    if row_gap > 2 * epsilon_row else math.inf
                )
                w_basis = right_row.conj().T[:, r:]
                v_s = w_basis.conj().T @ h_v @ w_basis
                v_s = (v_s + v_s.conj().T) / 2
                epsilon_vs = restriction_error(h_v, epsilon_hv, eta_s, n)
                values, vectors = la.eigh(-v_s)
                midpoint_negative = int(np.sum(values < 0))
                midpoint_positive = int(np.sum(values > 0))
                cluster_gap = float(values[15] - values[14])
                arithmetic_floor = float(
                    1000 * MACHINE_EPSILON * n
                    * max(1.0, operator_norm(v_s))
                )
                eta_eig = (
                    float(2 * epsilon_vs / (cluster_gap - 2 * epsilon_vs)
                          + arithmetic_floor)
                    if cluster_gap > 2 * epsilon_vs else math.inf
                )
                eta_p = float(2 * eta_s + eta_eig + arithmetic_floor)
                negative_basis = w_basis @ vectors[:, :15]
                projector = negative_basis @ negative_basis.conj().T
                hermitian_residual = operator_norm(projector - projector.conj().T)
                idempotent_residual = operator_norm(projector @ projector - projector)
                complete = bool(
                    rank_c == r
                    and row_rank == r and row_zero == 0 and row_open == 0
                    and midpoint_negative == 15 and midpoint_positive == 10
                    and cluster_gap > 2 * epsilon_vs
                    and math.isfinite(eta_p)
                    and hermitian_residual <= arithmetic_floor
                    and idempotent_residual <= arithmetic_floor
                )
                projector_controls &= complete
                projectors[time_name][key] = {
                    "projector": projector,
                    "eta": eta_p,
                }
                projector_records.append({
                    "time": time_name,
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "shape_rank": row_rank,
                    "midpoint_negative": midpoint_negative,
                    "midpoint_positive": midpoint_positive,
                    "cluster_gap": sf(cluster_gap),
                    "restricted_error": sf(epsilon_vs),
                    "shape_subspace_error": sf(eta_s),
                    "spectral_projector_error": sf(eta_eig),
                    "full_projector_error": sf(eta_p),
                    "hermitian_residual": sf(hermitian_residual),
                    "idempotent_residual": sf(idempotent_residual),
                    "complete": complete,
                })

check("all 32 old/shifted rank-15 projectors have separated certified gaps",
      projector_controls)

identity_records = []
identity_counts = Counter()
leakage_records = []
leakage_counts = Counter()
all_finite = True
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        n = 30
        one = np.eye(n, dtype=np.complex128)
        for variant in VARIANTS:
            key = (parity, sector_index, variant)
            old = projectors["old"][key]
            new = projectors["shifted"][key]
            p0 = old["projector"]
            p1 = new["projector"]
            floor_p = float(1000 * MACHINE_EPSILON * n)
            identity_distance = operator_norm(p1 - p0)
            identity_error = old["eta"] + new["eta"] + floor_p
            identity_label = projector_class(identity_distance, identity_error)
            identity_counts[identity_label] += 1
            identity_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "projector_distance": sf(identity_distance),
                "projector_error": sf(identity_error),
                "error_units": sf(
                    identity_distance / identity_error
                    if identity_error else math.inf
                ),
                "label": identity_label,
            })

            tangent_prefix = f"{parity}_sector{sector_index}_t2_{variant}"
            tangent_mid = np.asarray(
                tangent_archive[f"{tangent_prefix}_midpoint"]
            )
            tangent_rad = component_reenclosure_radii(
                tangent_mid,
                np.asarray(tangent_archive[f"{tangent_prefix}_radii"]),
            )
            blocks = {
                "A": tangent_mid[:n, :n],
                "B": tangent_mid[:n, n:],
                "C": tangent_mid[n:, :n],
                "D": tangent_mid[n:, n:],
            }
            block_radii = {
                "A": tangent_rad[:n, :n],
                "B": tangent_rad[:n, n:],
                "C": tangent_rad[n:, :n],
                "D": tangent_rad[n:, n:],
            }
            residuals = {
                "A": (one - p1) @ blocks["A"] @ p0,
                "B": (one - p1) @ blocks["B"] @ np.conj(p0),
                "C": (one - np.conj(p1)) @ blocks["C"] @ p0,
                "D": (one - np.conj(p1)) @ blocks["D"] @ np.conj(p0),
            }
            for block_name in ("A", "B", "C", "D"):
                block = blocks[block_name]
                epsilon_x = matrix_error(block, block_radii[block_name], n)
                residual_norm = operator_norm(residuals[block_name])
                residual_error = float(
                    epsilon_x
                    + (old["eta"] + new["eta"] + old["eta"] * new["eta"])
                    * (operator_norm(block) + epsilon_x)
                    + 1000 * MACHINE_EPSILON * n
                    * max(1.0, operator_norm(block))
                )
                label = leakage_class(residual_norm, residual_error)
                leakage_counts[label] += 1
                all_finite &= bool(
                    np.all(np.isfinite(block))
                    and math.isfinite(residual_norm)
                    and math.isfinite(residual_error)
                )
                leakage_records.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "block": block_name,
                    "residual_norm": sf(residual_norm),
                    "residual_error": sf(residual_error),
                    "error_units": sf(
                        residual_norm / residual_error
                        if residual_error else math.inf
                    ),
                    "label": label,
                })

identity_complete = len(identity_records) == 16 and sum(identity_counts.values()) == 16
leakage_complete = len(leakage_records) == 64 and sum(leakage_counts.values()) == 64
check("all 16 identity-aligned projector comparisons receive frozen labels",
      identity_complete, str(dict(identity_counts)))
check("all 64 canonical-tangent leakage blocks receive frozen labels",
      leakage_complete and all_finite, str(dict(leakage_counts)))

controls_ok = bool(
    provenance_ok and geometry_ok and carrier_geometry_ok and sector_basis_ok
    and projector_controls and identity_complete and leakage_complete and all_finite
)
if not controls_ok:
    outcome = "NEGATIVE_FIBER_TRANSPORT_CONTROL_FAILED"
elif leakage_counts["LEAKAGE_NONZERO_RESOLVED"]:
    outcome = "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED"
elif leakage_counts["LEAKAGE_OPEN"]:
    outcome = "NEGATIVE_FIBER_TANGENT_CLOSURE_OPEN"
else:
    outcome = "NEGATIVE_FIBER_TANGENT_CLOSURE_CERTIFIED"

allowed = {
    "NEGATIVE_FIBER_TRANSPORT_CONTROL_FAILED",
    "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED",
    "NEGATIVE_FIBER_TANGENT_CLOSURE_OPEN",
    "NEGATIVE_FIBER_TANGENT_CLOSURE_CERTIFIED",
}
check("the preregistered transport hierarchy assigns the outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "fitted_alignment_used": False,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "projector_records": projector_records,
    "identity_classification_counts": dict(identity_counts),
    "identity_records": identity_records,
    "leakage_classification_counts": dict(leakage_counts),
    "leakage_records": leakage_records,
    "classification": {
        "persistent_negative_rank": "DERIVED COMPUTATIONAL",
        "identity_aligned_fiber": "STRUCTURAL",
        "canonical_tangent_phase_fiber_closure": (
            "DERIVED COMPUTATIONAL REFUTATION"
            if outcome.endswith("REFUTED") else "OPEN"
        ),
        "alternative_constraint_derived_phase_lift": "OPEN",
        "physical_inertia_instability_graviton_or_wave": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"identity labels: {dict(identity_counts)}")
print(f"tangent leakage labels: {dict(leakage_counts)}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

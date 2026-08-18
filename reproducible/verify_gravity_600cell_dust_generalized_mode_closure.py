#!/usr/bin/env python3
"""Target-disclosed Hermitian-definite generalized-mode closure gate.

Prior-art/framing commit: 52ec90a.
Preregistered protocol commit: dfd833e.
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
TRANSPORT = HERE / "gravity_600cell_dust_negative_fiber_transport.json"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
OUTPUT = HERE / "gravity_600cell_dust_generalized_mode_closure.json"

PRIOR_ART_COMMIT = "52ec90a"
PROTOCOL_COMMIT = "dfd833e"
EXPECTED_HASHES = {
    "old_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "old_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "old_shape": "03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868",
    "new_json": "265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47",
    "new_npz": "c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8",
    "new_shape": "14fe5bc91e3ae4712c6ea19b8120785e2facd364e1ceb194009123fa353a4315",
    "direct": "86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944",
    "transport": "d630bf07066f88c35eee5a62a80ec1f43399a95ea882a43528289220c67f4599",
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
    return midpoint, component_reenclosure_radii(midpoint, stored)


def identity_class(distance, error):
    if not math.isfinite(distance) or not math.isfinite(error):
        return "GENERALIZED_FIBER_IDENTITY_OPEN"
    if distance <= 10 * error:
        return "GENERALIZED_COMMON_FIBER_RESOLVED"
    if distance > 100 * error:
        return "GENERALIZED_ROTATED_FIBER_RESOLVED"
    return "GENERALIZED_FIBER_IDENTITY_OPEN"


def leakage_class(value, error):
    if not math.isfinite(value) or not math.isfinite(error):
        return "LEAKAGE_OPEN"
    if value <= 10 * error:
        return "LEAKAGE_ZERO_CONSISTENT"
    if value > 100 * error:
        return "LEAKAGE_NONZERO_RESOLVED"
    return "LEAKAGE_OPEN"


paths = {
    "old_json": OLD_JSON, "old_npz": OLD_NPZ, "old_shape": OLD_SHAPE,
    "new_json": NEW_JSON, "new_npz": NEW_NPZ, "new_shape": NEW_SHAPE,
    "direct": DIRECT, "transport": TRANSPORT,
    "geometry_source": GEOMETRY_SOURCE, "full_source": FULL_SOURCE,
    "conformal_source": CONFORMAL_SOURCE, "commons": COMMONS_SOURCE,
}
hashes = {name: sha256(path) for name, path in paths.items()}
old_json = json.loads(OLD_JSON.read_text())
old_shape = json.loads(OLD_SHAPE.read_text())
new_json = json.loads(NEW_JSON.read_text())
new_shape = json.loads(NEW_SHAPE.read_text())
direct = json.loads(DIRECT.read_text())
transport = json.loads(TRANSPORT.read_text())
old_archive = np.load(OLD_NPZ, allow_pickle=False)
new_archive = np.load(NEW_NPZ, allow_pickle=False)
transport_pattern = Counter(
    (item["block"], item["label"]) for item in transport["leakage_records"]
)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and old_json["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and old_json["numeric_archive_arrays"] == len(old_archive.files) == 560
    and old_shape["outcome"] == "SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED"
    and new_json["outcome"] == "SHIFTED_CENTERED_CERTIFIED"
    and new_json["numeric_archive_arrays"] == len(new_archive.files) == 560
    and new_shape["outcome"] == "SHIFTED_SHAPE_STIFFNESS_SIGN_OPEN"
    and direct["outcome"] == "SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS"
    and transport["outcome"] == "NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED"
    and transport_pattern == Counter({
        ("A", "LEAKAGE_ZERO_CONSISTENT"): 16,
        ("B", "LEAKAGE_NONZERO_RESOLVED"): 16,
        ("C", "LEAKAGE_ZERO_CONSISTENT"): 16,
        ("D", "LEAKAGE_NONZERO_RESOLVED"): 16,
    })
)
check("all generalized-mode inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_generalized_mode", GEOMETRY_SOURCE
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
carrier_ok = all(
    incidences[parity]["equivariant"]
    and incidences[parity]["connected"]
    and incidences[parity]["numerical_rank"] == 120
    for parity in PARITIES
)
sector_data, sector_controls = high_precision_sector_bases(groups["even"])
sector_ok = bool(
    tuple(item["dimension"] for item in sector_data) == (3, 2, 2, 2, 1, 1, 1)
    and sector_controls["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
)
check("the exact carrier and deterministic sector basis are reconstructed",
      carrier_ok and sector_ok)

print("=" * 78)
print("KINETIC--STIFFNESS GENERALIZED-MODE RECURRENCE CLOSURE")
print("=" * 78)

projectors = {"old": {}, "shifted": {}}
projector_records = []
leakage_records = []
leakage_counts = Counter()
controls = True
all_finite = True

for parity in PARITIES:
    incidence = incidences[parity]["incidence"].astype(np.complex128)
    for sector_index in TARGET_SECTORS:
        dimension = sector_data[sector_index]["dimension"]
        n, r, s = 30 * dimension, 5 * dimension, 25 * dimension
        basis = mp_to_numpy(sector_data[sector_index]["basis"])
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
                matrices = {}
                errors = {}
                for name in ("M", "V", "Gamma", "Omega"):
                    midpoint, radii = source_matrix(archive, prefix, name)
                    matrices[name] = midpoint
                    errors[name] = matrix_error(midpoint, radii, n)

                h_m = (matrices["M"] + matrices["M"].conj().T) / 2
                h_v = (matrices["V"] + matrices["V"].conj().T) / 2
                epsilon_hm = errors["M"]
                epsilon_hv = errors["V"]
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
                m_s = w_basis.conj().T @ h_m @ w_basis
                m_s = (m_s + m_s.conj().T) / 2
                v_s = w_basis.conj().T @ h_v @ w_basis
                v_s = (v_s + v_s.conj().T) / 2
                epsilon_ms = restriction_error(h_m, epsilon_hm, eta_s, n)
                epsilon_vs = restriction_error(h_v, epsilon_hv, eta_s, n)
                b_matrix = -m_s
                a_matrix = -v_s
                b_values = la.eigvalsh(b_matrix)
                minimum_b = float(b_values[0])
                maximum_b = float(b_values[-1])
                b_lower = minimum_b - epsilon_ms
                kinetic_ok = b_lower > 0 and minimum_b > 100 * epsilon_ms
                if kinetic_ok:
                    values, vectors = la.eigh(a_matrix, b_matrix)
                    epsilon_pencil = float(
                        epsilon_vs / b_lower
                        + operator_norm(a_matrix) * epsilon_ms
                        / (minimum_b * b_lower)
                        + 1000 * MACHINE_EPSILON * s
                        * max(1.0, operator_norm(la.solve(b_matrix, a_matrix)))
                    )
                else:
                    values = np.full(s, np.nan)
                    vectors = np.full((s, s), np.nan, dtype=np.complex128)
                    epsilon_pencil = math.inf
                negative = int(np.sum(values < 0))
                positive = int(np.sum(values > 0))
                gap = float(values[15] - values[14])
                floor = float(1000 * MACHINE_EPSILON * n)
                eta_eig = (
                    float(2 * epsilon_pencil / (gap - 2 * epsilon_pencil) + floor)
                    if gap > 2 * epsilon_pencil else math.inf
                )
                eta_p = float(
                    2 * eta_s
                    + math.sqrt(maximum_b / b_lower) * eta_eig
                    + epsilon_ms / b_lower
                    + floor
                ) if kinetic_ok and math.isfinite(eta_eig) else math.inf
                lifted = w_basis @ vectors[:, :15]
                q_basis, _ = la.qr(lifted, mode="economic")
                projector = q_basis @ q_basis.conj().T
                projector_residual = max(
                    operator_norm(projector - projector.conj().T),
                    operator_norm(projector @ projector - projector),
                )
                complete = bool(
                    rank_c == r and row_rank == r and row_zero == 0 and row_open == 0
                    and kinetic_ok and negative == 15 and positive == 10
                    and gap > 2 * epsilon_pencil and math.isfinite(eta_p)
                    and projector_residual <= floor
                )
                controls &= complete
                projectors[time_name][key] = {"projector": projector, "eta": eta_p}
                projector_records.append({
                    "time": time_name,
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "negative": negative,
                    "positive": positive,
                    "kinetic_minimum": sf(minimum_b),
                    "kinetic_error": sf(epsilon_ms),
                    "generalized_gap": sf(gap),
                    "pencil_error": sf(epsilon_pencil),
                    "projector_error": sf(eta_p),
                    "complete": complete,
                })

                one = np.eye(n, dtype=np.complex128)
                for operator_name in ("Gamma", "Omega"):
                    operator = matrices[operator_name]
                    residual = (one - projector) @ operator @ projector
                    residual_norm = operator_norm(residual)
                    epsilon_x = errors[operator_name]
                    residual_error = float(
                        epsilon_x
                        + (2 * eta_p + eta_p**2)
                        * (operator_norm(operator) + epsilon_x)
                        + 1000 * MACHINE_EPSILON * n
                        * max(1.0, operator_norm(operator))
                    )
                    label = leakage_class(residual_norm, residual_error)
                    leakage_counts[(time_name, operator_name, label)] += 1
                    all_finite &= bool(
                        np.all(np.isfinite(operator))
                        and math.isfinite(residual_norm)
                        and math.isfinite(residual_error)
                    )
                    leakage_records.append({
                        "time": time_name,
                        "parity": parity,
                        "sector_index": sector_index,
                        "variant": variant,
                        "operator": operator_name,
                        "residual_norm": sf(residual_norm),
                        "residual_error": sf(residual_error),
                        "error_units": sf(
                            residual_norm / residual_error
                            if residual_error else math.inf
                        ),
                        "label": label,
                    })

check("all 32 Hermitian-definite generalized projectors are resolved", controls)

identity_records = []
identity_counts = Counter()
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for variant in VARIANTS:
            key = (parity, sector_index, variant)
            old = projectors["old"][key]
            new = projectors["shifted"][key]
            distance = operator_norm(new["projector"] - old["projector"])
            error = old["eta"] + new["eta"] + 1000 * MACHINE_EPSILON * 30
            label = identity_class(distance, error)
            identity_counts[label] += 1
            identity_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "distance": sf(distance),
                "error": sf(error),
                "error_units": sf(distance / error if error else math.inf),
                "label": label,
            })

identity_complete = len(identity_records) == 16
leakage_complete = len(leakage_records) == 64 and sum(leakage_counts.values()) == 64
check("all 16 generalized old/shifted fiber comparisons receive labels",
      identity_complete, str(dict(identity_counts)))
check("all 64 Gamma/Omega generalized-fiber leakages receive labels",
      leakage_complete and all_finite, str(dict(leakage_counts)))

flat_leakage_counts = Counter(item["label"] for item in leakage_records)
controls_ok = bool(
    provenance_ok and geometry_ok and carrier_ok and sector_ok and controls
    and identity_complete and leakage_complete and all_finite
)
if not controls_ok:
    outcome = "GENERALIZED_MODE_CONTROL_FAILED"
elif flat_leakage_counts["LEAKAGE_NONZERO_RESOLVED"]:
    outcome = "GENERALIZED_MODE_RECURRENCE_CLOSURE_REFUTED"
elif flat_leakage_counts["LEAKAGE_OPEN"]:
    outcome = "GENERALIZED_MODE_RECURRENCE_CLOSURE_OPEN"
else:
    outcome = "GENERALIZED_MODE_RECURRENCE_CLOSURE_CERTIFIED"

allowed = {
    "GENERALIZED_MODE_CONTROL_FAILED",
    "GENERALIZED_MODE_RECURRENCE_CLOSURE_REFUTED",
    "GENERALIZED_MODE_RECURRENCE_CLOSURE_OPEN",
    "GENERALIZED_MODE_RECURRENCE_CLOSURE_CERTIFIED",
}
check("the preregistered generalized-mode hierarchy assigns the outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "projector_records": projector_records,
    "identity_counts": dict(identity_counts),
    "identity_records": identity_records,
    "leakage_counts": {
        "|".join(key): value for key, value in sorted(leakage_counts.items())
    },
    "leakage_records": leakage_records,
    "classification": {
        "generalized_pencil_fiber": "DERIVED COMPUTATIONAL",
        "recurrence_closure": (
            "DERIVED COMPUTATIONAL REFUTATION"
            if outcome.endswith("REFUTED") else "OPEN"
        ),
        "riccati_or_constraint_phase_lift": "OPEN",
        "particle_inertia_mass_or_dispersion": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"identity labels: {dict(identity_counts)}")
print(f"leakage labels: {dict(flat_leakage_counts)}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)

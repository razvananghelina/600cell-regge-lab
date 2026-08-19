#!/usr/bin/env python3
"""Target-disclosed action-weighted longitudinal audit of negative shape modes."""

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


OUTPUT = HERE / "gravity_600cell_dust_action_york_negative.json"
CENTERED_JSON = HERE / "gravity_600cell_dust_centered_jacobi.json"
CENTERED_NPZ = HERE / "gravity_600cell_dust_centered_jacobi.npz"
STIFFNESS_JSON = HERE / "gravity_600cell_dust_shape_stiffness.json"
RIGIDITY_JSON = HERE / "gravity_600cell_dust_rigidity_york.json"
SHAPE_SOURCE = HERE / "verify_gravity_600cell_dust_shape_stiffness.py"
RIGIDITY_SOURCE = HERE / "verify_gravity_600cell_dust_rigidity_york.py"
NEGATIVE_SOURCE = HERE / "verify_gravity_600cell_dust_negative_shape_dynamics.py"
CONFORMAL_SOURCE = HERE / "verify_gravity_600cell_dust_conformal_supermetric.py"
FULL_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
COMMONS_SOURCE = ROOT / "commons" / "cell600.py"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_dust_action_york_negative_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_dust_action_york_negative_protocol.md"

PRIOR_ART_COMMIT = "9bd990a"
PROTOCOL_COMMIT = "7165752"
EXPECTED_HASHES = {
    "prior_art": "e5728865b8498c5750cdbf45d9d93938c530ba40d95c441452f34c02bf00cd1d",
    "protocol": "d4532d626c12c4f2d8af39fded1f59a663329da73fa3d5243c12e4135828cc8f",
    "centered_json": "fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56",
    "centered_npz": "1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef",
    "stiffness_json": "03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868",
    "rigidity_json": "251851c08f81ba2f0d41c2d0da428ab11f1ba918b9cb59e0a1e347143c883981",
    "shape_source": "d4f0a9a805910de37011ba70f407907daa2d11c650aeea22e571ab867282a44c",
    "rigidity_source": "deba8d9f9bca4a5848134943ec77544e5487d44a59c44234f632b6f2aeb51382",
    "negative_source": "6e7659ca398037e806f9a35a9f3db3d6035f992a8655699b47a2519b0c37453e",
    "conformal_source": "d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4",
    "full_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "commons": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
PARITIES = ("even", "odd")
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)
DIMENSIONS = (3, 2, 2, 2, 1, 1, 1)
SELECTED_SECTORS = (4, 5)
MACHINE_EPSILON = np.finfo(float).eps
COORDINATE_ENVELOPE = 100 * 5.7e-11
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
    values = la.svdvals(matrix)
    return float(values[0]) if len(values) else 0.0


def zero_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "OPEN"
    if value <= 10 * error:
        return "ZERO_CONSISTENT"
    if value > 100 * error:
        return "NONZERO_RESOLVED"
    return "OPEN"


def equality_label(value, error):
    label = zero_label(value, error)
    return {
        "ZERO_CONSISTENT": "EQUALITY_CONSISTENT",
        "NONZERO_RESOLVED": "SEPARATED_RESOLVED",
        "OPEN": "EQUALITY_OPEN",
    }[label]


def sign_label(value, error):
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return "OPEN"
    if value > 100 * error:
        return "POSITIVE_RESOLVED"
    if value < -100 * error:
        return "NEGATIVE_RESOLVED"
    if abs(value) <= 10 * error:
        return "ZERO_CONSISTENT"
    return "OPEN"


def rank_record(matrix, error):
    singular = la.svdvals(matrix)
    nonzero = singular > 100 * error
    zero = singular < 10 * error
    return {
        "rank": int(np.sum(nonzero)),
        "zero": int(np.sum(zero)),
        "open": int(len(singular)-np.sum(nonzero)-np.sum(zero)),
        "singular": singular,
        "error": float(error),
    }


def load_audited_helpers():
    wanted = {
        "mp_frobenius",
        "mp_submatrix",
        "cluster_sorted",
        "orbit_sort_key",
        "edge_image",
        "group_data",
        "incidence_data",
        "mp_to_numpy",
        "component_reenclosure_radii",
    }
    tree = ast.parse(CONFORMAL_SOURCE.read_text(), filename=str(CONFORMAL_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    if {node.name for node in body} != wanted:
        raise RuntimeError("audited conformal helper set changed")
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
        raise RuntimeError("audited sector helper changed")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(FULL_SOURCE), "exec"),
        globals(),
    )

    tree = ast.parse(RIGIDITY_SOURCE.read_text(), filename=str(RIGIDITY_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "build_rigidity"
    ]
    if len(body) != 1:
        raise RuntimeError("audited rigidity helper changed")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(RIGIDITY_SOURCE), "exec"),
        globals(),
    )


def source_matrix(archive, prefix, name):
    midpoint = np.asarray(archive[f"{prefix}_{name}_midpoint"])
    stored = np.asarray(archive[f"{prefix}_{name}_radii"])
    radii = component_reenclosure_radii(midpoint, stored)
    return midpoint, radii


def matrix_error(midpoint, radii, n):
    return float(
        la.norm(radii, "fro")
        + 1000 * MACHINE_EPSILON * n
        * max(1.0, operator_norm(midpoint))
    )


def restriction_error(midpoint, base_error, eta_shape, n):
    norm = operator_norm(midpoint)
    return float(
        base_error
        + 2 * eta_shape * (norm + base_error)
        + 1000 * MACHINE_EPSILON * n * max(1.0, norm)
    )


print("="*78)
print("ACTION-WEIGHTED LONGITUDINAL IDENTITY OF NEGATIVE SHAPE MODES")
print("="*78)

paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "centered_json": CENTERED_JSON,
    "centered_npz": CENTERED_NPZ,
    "stiffness_json": STIFFNESS_JSON,
    "rigidity_json": RIGIDITY_JSON,
    "shape_source": SHAPE_SOURCE,
    "rigidity_source": RIGIDITY_SOURCE,
    "negative_source": NEGATIVE_SOURCE,
    "conformal_source": CONFORMAL_SOURCE,
    "full_source": FULL_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "commons": COMMONS_SOURCE,
}
hashes = {name: digest(path) for name, path in paths.items()}
centered = json.loads(CENTERED_JSON.read_text())
stiffness = json.loads(STIFFNESS_JSON.read_text())
rigidity_input = json.loads(RIGIDITY_JSON.read_text())
archive = np.load(CENTERED_NPZ, allow_pickle=False)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and PRIOR_ART_COMMIT == "9bd990a"
    and PROTOCOL_COMMIT == "7165752"
    and centered["outcome"] == "CENTERED_JACOBI_CERTIFIED"
    and centered["passed"] == centered["tests"] == 7
    and len(archive.files) == centered["numeric_archive_arrays"] == 560
    and stiffness["outcome"] == "SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED"
    and stiffness["passed"] == stiffness["tests"] == 12
    and rigidity_input["outcome"] == "RIGIDITY_YORK_DECOUPLING_REFUTED"
    and rigidity_input["passed"] == rigidity_input["tests"] == 9
)
check("all frozen action, stiffness and geometry inputs have exact provenance",
      provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_action_york", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check("the independently imported complete slab geometry retains 43/43 controls",
      gro.tests == gro.passed == 43)

load_audited_helpers()
groups = {parity: group_data(gro.models[parity], gro) for parity in PARITIES}
incidences = {parity: incidence_data(groups[parity]) for parity in PARITIES}
vertices, _, _ = build_600cell()

rigidity = {}
global_controls = {}
for parity in PARITIES:
    r_matrix, radial, tangent, lengths, length_square = build_rigidity(
        vertices, groups[parity]["edge_order"]
    )
    d_matrix = r_matrix @ tangent
    c_matrix = incidences[parity]["incidence"].astype(float)
    rank_c = int(np.linalg.matrix_rank(c_matrix))
    rank_r = int(np.linalg.matrix_rank(r_matrix))
    rank_d = int(np.linalg.matrix_rank(d_matrix))
    rank_cd = int(np.linalg.matrix_rank(np.column_stack((c_matrix, d_matrix))))
    intersection = rank_c+rank_d-rank_cd
    rigidity[parity] = {
        "R": r_matrix,
        "D": d_matrix,
        "C": c_matrix,
    }
    global_controls[parity] = {
        "rank_C": rank_c,
        "rank_R": rank_r,
        "rank_D": rank_d,
        "intersection_C_D": intersection,
        "length_spread": float(np.ptp(lengths)),
        "radial_conformal_residual": operator_norm(r_matrix @ radial-c_matrix),
    }
check(
    "the full conformal/rigidity/tangent carrier has ranks 120/470/354 and overlap 4",
    all(
        record["rank_C"] == 120 and record["rank_R"] == 470
        and record["rank_D"] == 354 and record["intersection_C_D"] == 4
        and record["length_spread"] < 1e-9
        and record["radial_conformal_residual"] < 1e-8
        for record in global_controls.values()
    ),
    str(global_controls),
)

sector_data, sector_controls = high_precision_sector_bases(groups["even"])
maximum_basis_residual = max(
    value for key, value in sector_controls.items() if key.startswith("maximum_")
)
sector_ok = bool(
    tuple(sector["dimension"] for sector in sector_data) == DIMENSIONS
    and maximum_basis_residual < mp.mpf("1e-70")
)
check("all seven minimal binary-tetrahedral sectors are reconstructed",
      sector_ok, "maximum residual="+mp.nstr(maximum_basis_residual, 5))

sector_geometry = {parity: {} for parity in PARITIES}
geometry_open = False
for parity in PARITIES:
    incidence = rigidity[parity]["C"].astype(np.complex128)
    d_full = rigidity[parity]["D"].astype(np.complex128)
    for sector_index, sector in enumerate(sector_data):
        dimension = sector["dimension"]
        n = 30*dimension
        r = 5*dimension
        basis = mp_to_numpy(sector["basis"])
        edge_basis = np.kron(np.eye(30, dtype=np.complex128), basis)
        compressed_c = edge_basis.conj().T @ incidence
        left, singular_c, _ = la.svd(compressed_c, full_matrices=False)
        epsilon_c = float(
            1000*MACHINE_EPSILON*max(compressed_c.shape)
            * max(1.0, float(singular_c[0]))
            + COORDINATE_ENVELOPE
        )
        rank_c = int(np.sum(singular_c > 100*epsilon_c))
        zero_c = int(np.sum(singular_c < 10*epsilon_c))
        open_c = len(singular_c)-rank_c-zero_c
        gap_c = float(singular_c[r-1])
        eta_k = (
            float(2*epsilon_c/(gap_c-2*epsilon_c)
                  + 1000*MACHINE_EPSILON*n)
            if gap_c > 2*epsilon_c else math.inf
        )
        geometry_open |= bool(rank_c != r or open_c or not math.isfinite(eta_k))
        sector_geometry[parity][sector_index] = {
            "U": left[:, :r],
            "D": edge_basis.conj().T @ d_full,
            "eta_K": eta_k,
            "rank_C": rank_c,
        }
check("every sector has the exact 5d conformal image before stiffness is loaded",
      not geometry_open)

records = {parity: [] for parity in PARITIES}
internal = {}
dimension_census = {}
projector_counts = Counter()
cross_counts = Counter()
sign_counts = Counter()
dynamic_counts = Counter()
control_counts = Counter()
all_projection_identities = True
all_kinetic = True
all_finite = True

for parity in PARITIES:
    for sector_index, dimension in enumerate(DIMENSIONS):
        n = 30*dimension
        r = 5*dimension
        s = 25*dimension
        u_basis = sector_geometry[parity][sector_index]["U"]
        d_matrix = sector_geometry[parity][sector_index]["D"]
        eta_k = sector_geometry[parity][sector_index]["eta_K"]
        sector_variants = []

        for variant in VARIANTS:
            prefix = f"{parity}_sector{sector_index}_{variant}"
            midpoint_m, radius_m = source_matrix(archive, prefix, "M")
            midpoint_v, radius_v = source_matrix(archive, prefix, "V")
            midpoint_g, radius_g = source_matrix(archive, prefix, "Gamma")
            midpoint_o, radius_o = source_matrix(archive, prefix, "Omega")
            h_m = (midpoint_m+midpoint_m.conj().T)/2
            h_v = (midpoint_v+midpoint_v.conj().T)/2
            radius_hm = (radius_m+radius_m.T)/2
            radius_hv = (radius_v+radius_v.T)/2
            epsilon_hm = matrix_error(h_m, radius_hm, n)
            epsilon_hv = matrix_error(h_v, radius_hv, n)
            epsilon_g = matrix_error(midpoint_g, radius_g, n)
            epsilon_o = matrix_error(midpoint_o, radius_o, n)

            row = u_basis.conj().T @ h_m
            _, singular_row, right_row = la.svd(row, full_matrices=True)
            epsilon_row = float(
                epsilon_hm
                + 2*eta_k*(operator_norm(h_m)+epsilon_hm)
                + 1000*MACHINE_EPSILON*n*max(1.0, operator_norm(row))
            )
            row_rank = int(np.sum(singular_row > 100*epsilon_row))
            row_zero = int(np.sum(singular_row < 10*epsilon_row))
            row_open = len(singular_row)-row_rank-row_zero
            row_gap = float(singular_row[r-1])
            eta_s = (
                float(2*epsilon_row/(row_gap-2*epsilon_row)
                      + 1000*MACHINE_EPSILON*n)
                if row_gap > 2*epsilon_row else math.inf
            )
            w_basis = right_row.conj().T[:, r:]
            q_basis = np.column_stack((u_basis, w_basis))
            q_condition = float(np.linalg.cond(q_basis))
            coefficients = la.solve(q_basis, d_matrix)
            longitudinal_raw = coefficients[r:, :]
            longitudinal_input_error = float(
                (2*q_condition+1)*eta_s
                * max(1.0, operator_norm(longitudinal_raw))
                + COORDINATE_ENVELOPE*q_condition
                * max(1.0, operator_norm(d_matrix))
                + 1000*MACHINE_EPSILON*max(q_basis.shape)
                * max(1.0, operator_norm(longitudinal_raw))
            )
            longitudinal_rank = rank_record(
                longitudinal_raw, longitudinal_input_error
            )
            l_dimension = longitudinal_rank["rank"]
            left_l, singular_l, _ = la.svd(longitudinal_raw, full_matrices=False)
            l_basis = left_l[:, :l_dimension]
            l_gap = float(singular_l[l_dimension-1]) if l_dimension else 0.0
            eta_l = (
                float(2*longitudinal_input_error
                      /(l_gap-2*longitudinal_input_error)
                      + 1000*MACHINE_EPSILON*s)
                if l_gap > 2*longitudinal_input_error else math.inf
            )

            inverse_q = la.inv(q_basis)
            selector = np.zeros((s, n), dtype=np.complex128)
            selector[:, r:] = np.eye(s)
            p_shape = w_basis @ selector @ inverse_q
            identity_residuals = {
                "idempotent": operator_norm(p_shape@p_shape-p_shape),
                "kills_conformal": operator_norm(p_shape@u_basis),
                "lands_in_shape": operator_norm(u_basis.conj().T@h_m@p_shape),
            }
            projection_error = float(
                2*eta_s*q_condition*(1+operator_norm(h_m))
                + epsilon_hm*q_condition
                + COORDINATE_ENVELOPE
                + 1000*MACHINE_EPSILON*n
                * max(1.0, operator_norm(p_shape), operator_norm(h_m))
            )
            projection_labels = {
                name: zero_label(value, projection_error)
                for name, value in identity_residuals.items()
            }
            all_projection_identities &= all(
                label == "ZERO_CONSISTENT"
                for label in projection_labels.values()
            )

            m_s = (w_basis.conj().T@h_m@w_basis)
            m_s = (m_s+m_s.conj().T)/2
            v_s = (w_basis.conj().T@h_v@w_basis)
            v_s = (v_s+v_s.conj().T)/2
            a_matrix = -v_s
            b_matrix = -m_s
            epsilon_ms = restriction_error(h_m, epsilon_hm, eta_s, n)
            epsilon_vs = restriction_error(h_v, epsilon_hv, eta_s, n)
            b_values = la.eigvalsh(b_matrix)
            minimum_b = float(b_values[0])
            maximum_b = float(b_values[-1])
            b_lower = minimum_b-epsilon_ms
            kinetic_ok = minimum_b > 100*epsilon_ms and b_lower > 0
            all_kinetic &= kinetic_ok

            # B-orthogonal complement of the geometry-selected longitudinal image.
            row_t = l_basis.conj().T @ b_matrix
            _, singular_t, right_t = la.svd(row_t, full_matrices=True)
            t_dimension = s-l_dimension
            t_basis = right_t.conj().T[:, l_dimension:]
            t_residual = operator_norm(l_basis.conj().T@b_matrix@t_basis)
            eta_t = float(
                2*eta_l + epsilon_ms/max(b_lower, 1e-300)
                + 1000*MACHINE_EPSILON*s
            ) if math.isfinite(eta_l) and b_lower > 0 else math.inf

            if (parity, sector_index, variant) not in internal:
                internal[(parity, sector_index, variant)] = {}
            if sector_index not in dimension_census:
                dimension_census[sector_index] = set()
            dimension_census[sector_index].add((l_dimension, t_dimension))

            selected_record = None
            if sector_index in SELECTED_SECTORS and kinetic_ok:
                generalized_values, generalized_vectors = la.eigh(a_matrix, b_matrix)
                epsilon_pencil = float(
                    epsilon_vs/b_lower
                    + operator_norm(a_matrix)*epsilon_ms/(minimum_b*b_lower)
                    + 1000*MACHINE_EPSILON*s
                    * max(1.0, operator_norm(la.solve(b_matrix, a_matrix)))
                )
                generalized_labels = [
                    sign_label(float(value), epsilon_pencil)
                    for value in generalized_values
                ]
                generalized_inertia = Counter(generalized_labels)
                generalized_gap = float(
                    generalized_values[15]-generalized_values[14]
                )
                eta_eig = (
                    float(2*epsilon_pencil/(generalized_gap-2*epsilon_pencil)
                          + 1000*MACHINE_EPSILON*s)
                    if generalized_gap > 2*epsilon_pencil else math.inf
                )
                eta_p = (
                    float(
                        2*eta_s + math.sqrt(maximum_b/b_lower)*eta_eig
                        + epsilon_ms/b_lower
                        + 1000*MACHINE_EPSILON*n
                    )
                    if math.isfinite(eta_eig) else math.inf
                )
                e_lifted = generalized_vectors[:, :15]
                e_basis, _ = la.qr(e_lifted, mode="economic")
                p_l = l_basis@l_basis.conj().T
                p_e = e_basis@e_basis.conj().T
                projector_distance = operator_norm(p_l-p_e)
                comparison_error = float(
                    eta_l+eta_p+1000*MACHINE_EPSILON*s
                )
                projector_label = equality_label(
                    projector_distance, comparison_error
                )
                projector_counts[projector_label] += 1

                cross_a = operator_norm(l_basis.conj().T@a_matrix@t_basis)
                cross_b = t_residual
                cross_error = float(
                    epsilon_vs
                    + 2*(eta_l+eta_t)*(operator_norm(a_matrix)+epsilon_vs)
                    + 1000*MACHINE_EPSILON*s
                    * max(1.0, operator_norm(a_matrix))
                )
                cross_b_error = float(
                    epsilon_ms
                    + 2*(eta_l+eta_t)*(operator_norm(b_matrix)+epsilon_ms)
                    + 1000*MACHINE_EPSILON*s
                    * max(1.0, operator_norm(b_matrix))
                )
                cross_a_label = zero_label(cross_a, cross_error)
                cross_b_label = zero_label(cross_b, cross_b_error)
                cross_counts["A_"+cross_a_label] += 1
                cross_counts["B_"+cross_b_label] += 1

                a_l = (l_basis.conj().T@a_matrix@l_basis)
                a_l = (a_l+a_l.conj().T)/2
                a_t = (t_basis.conj().T@a_matrix@t_basis)
                a_t = (a_t+a_t.conj().T)/2
                values_l = la.eigvalsh(a_l)
                values_t = la.eigvalsh(a_t)
                restriction_sign_error = float(
                    epsilon_vs
                    + 2*max(eta_l, eta_t)*(operator_norm(a_matrix)+epsilon_vs)
                    + 1000*MACHINE_EPSILON*s
                    * max(1.0, operator_norm(a_matrix))
                )
                l_sign = sign_label(float(values_l[-1]), restriction_sign_error)
                t_sign = sign_label(float(values_t[0]), restriction_sign_error)
                sign_counts["L_"+l_sign] += 1
                sign_counts["T_"+t_sign] += 1

                gamma_s = w_basis.conj().T@midpoint_g@w_basis
                omega_s = w_basis.conj().T@midpoint_o@w_basis
                epsilon_gs = restriction_error(midpoint_g, epsilon_g, eta_s, n)
                epsilon_os = restriction_error(midpoint_o, epsilon_o, eta_s, n)
                dynamic = {}
                for name, matrix, base_error in (
                    ("Gamma", gamma_s, epsilon_gs),
                    ("Omega", omega_s, epsilon_os),
                ):
                    residual = operator_norm(
                        (np.eye(s)-p_l)@matrix@l_basis
                    )
                    error = float(
                        base_error
                        + (2*eta_l+eta_l**2)
                        * (operator_norm(matrix)+base_error)
                        + 1000*MACHINE_EPSILON*s
                        * max(1.0, operator_norm(matrix))
                    )
                    label = zero_label(residual, error)
                    dynamic_counts[name+"_"+label] += 1
                    dynamic[name] = {
                        "residual": sf(residual),
                        "error": sf(error),
                        "label": label,
                    }

                # Same-dimensional control: replace one longitudinal direction
                # by the first independently selected B-orthogonal direction.
                control_seed = np.column_stack((l_basis[:, 1:], t_basis[:, :1]))
                control_basis, _ = la.qr(control_seed, mode="economic")
                p_control = control_basis@control_basis.conj().T
                control_distance = operator_norm(p_control-p_e)
                control_label = equality_label(
                    control_distance, comparison_error
                )
                control_counts[control_label] += 1

                selected_record = {
                    "generalized_inertia": dict(generalized_inertia),
                    "generalized_gap": sf(generalized_gap),
                    "pencil_error": sf(epsilon_pencil),
                    "projector_distance": sf(projector_distance),
                    "projector_error": sf(comparison_error),
                    "projector_label": projector_label,
                    "A_cross": sf(cross_a),
                    "A_cross_error": sf(cross_error),
                    "A_cross_label": cross_a_label,
                    "B_cross": sf(cross_b),
                    "B_cross_error": sf(cross_b_error),
                    "B_cross_label": cross_b_label,
                    "longitudinal_maximum_eigenvalue": sf(values_l[-1]),
                    "transverse_minimum_eigenvalue": sf(values_t[0]),
                    "restriction_sign_error": sf(restriction_sign_error),
                    "longitudinal_sign_label": l_sign,
                    "transverse_sign_label": t_sign,
                    "dynamic": dynamic,
                    "rotated_control_distance": sf(control_distance),
                    "rotated_control_label": control_label,
                }

            all_finite &= bool(
                math.isfinite(q_condition)
                and math.isfinite(longitudinal_input_error)
                and np.all(np.isfinite(singular_l))
                and math.isfinite(projection_error)
            )
            sector_variants.append({
                "variant": variant,
                "shape_dimension": s,
                "longitudinal_dimension": l_dimension,
                "transverse_dimension": t_dimension,
                "longitudinal_rank_open": longitudinal_rank["open"],
                "longitudinal_rank_error": sf(longitudinal_input_error),
                "longitudinal_minimum_singular": sf(l_gap),
                "longitudinal_subspace_error": sf(eta_l),
                "shape_subspace_error": sf(eta_s),
                "direct_sum_condition": sf(q_condition),
                "projection_identities": {
                    name: {
                        "residual": sf(value),
                        "error": sf(projection_error),
                        "label": projection_labels[name],
                    }
                    for name, value in identity_residuals.items()
                },
                "kinetic_minimum": sf(minimum_b),
                "kinetic_error": sf(epsilon_ms),
                "selected_comparison": selected_record,
            })

        records[parity].append({
            "sector_index": sector_index,
            "irrep_dimension": dimension,
            "variants": sector_variants,
        })

census_uniform = all(len(values) == 1 for values in dimension_census.values())
census_pairs = {
    sector: next(iter(values)) if len(values) == 1 else sorted(values)
    for sector, values in dimension_census.items()
}
weighted_longitudinal = sum(
    DIMENSIONS[sector]*census_pairs[sector][0] for sector in range(7)
) if census_uniform else -1
weighted_transverse = sum(
    DIMENSIONS[sector]*census_pairs[sector][1] for sector in range(7)
) if census_uniform else -1
check(
    "the all-sector longitudinal census is uniform across schedules and variants",
    census_uniform and all_finite,
    str(census_pairs),
)
check(
    "multiplicity restoration gives exactly 350 longitudinal plus 250 transverse shapes",
    weighted_longitudinal == 350 and weighted_transverse == 250,
    f"longitudinal={weighted_longitudinal}, transverse={weighted_transverse}",
)
check(
    "all oblique action-shape projector identities are zero-consistent",
    all_projection_identities,
)
check(
    "all 56 shape kinetic restrictions remain positive resolved",
    all_kinetic,
)

required_selected = 2*2*4
selected_dimensions_ok = all(
    census_pairs[sector] == (15, 10) for sector in SELECTED_SECTORS
)
check(
    "the disclosed sectors have geometry-selected dimensions 15 plus 10",
    selected_dimensions_ok,
    str({sector: census_pairs[sector] for sector in SELECTED_SECTORS}),
)
check(
    "all selected generalized-negative projectors equal the longitudinal carrier",
    projector_counts["EQUALITY_CONSISTENT"] == required_selected
    and sum(projector_counts.values()) == required_selected,
    str(dict(projector_counts)),
)
check(
    "all selected A and B longitudinal/transverse cross blocks vanish",
    cross_counts["A_ZERO_CONSISTENT"] == required_selected
    and cross_counts["B_ZERO_CONSISTENT"] == required_selected
    and sum(cross_counts.values()) == 2*required_selected,
    str(dict(cross_counts)),
)
check(
    "longitudinal stiffness is negative and its complement positive in every selected cell",
    sign_counts["L_NEGATIVE_RESOLVED"] == required_selected
    and sign_counts["T_POSITIVE_RESOLVED"] == required_selected
    and sum(sign_counts.values()) == 2*required_selected,
    str(dict(sign_counts)),
)
check(
    "Gamma and Omega preserve every selected longitudinal carrier",
    dynamic_counts["Gamma_ZERO_CONSISTENT"] == required_selected
    and dynamic_counts["Omega_ZERO_CONSISTENT"] == required_selected
    and sum(dynamic_counts.values()) == 2*required_selected,
    str(dict(dynamic_counts)),
)
check(
    "the same-dimensional rotated control is separated in every selected cell",
    control_counts["SEPARATED_RESOLVED"] == required_selected
    and sum(control_counts.values()) == required_selected,
    str(dict(control_counts)),
)

controls_ok = bool(
    provenance_ok and gro.tests == gro.passed == 43 and sector_ok
    and census_uniform and all_projection_identities and all_kinetic and all_finite
    and weighted_longitudinal == 350 and weighted_transverse == 250
    and selected_dimensions_ok
)
positive_identity = bool(
    controls_ok
    and projector_counts["EQUALITY_CONSISTENT"] == required_selected
    and cross_counts["A_ZERO_CONSISTENT"] == required_selected
    and cross_counts["B_ZERO_CONSISTENT"] == required_selected
    and sign_counts["L_NEGATIVE_RESOLVED"] == required_selected
    and sign_counts["T_POSITIVE_RESOLVED"] == required_selected
    and dynamic_counts["Gamma_ZERO_CONSISTENT"] == required_selected
    and dynamic_counts["Omega_ZERO_CONSISTENT"] == required_selected
    and control_counts["SEPARATED_RESOLVED"] == required_selected
)
resolved_refutation = bool(
    projector_counts["SEPARATED_RESOLVED"]
    or cross_counts["A_NONZERO_RESOLVED"]
    or cross_counts["B_NONZERO_RESOLVED"]
    or dynamic_counts["Gamma_NONZERO_RESOLVED"]
    or dynamic_counts["Omega_NONZERO_RESOLVED"]
)
if positive_identity:
    outcome = "NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_RESOLVED"
elif controls_ok and resolved_refutation:
    outcome = "NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_REFUTED"
else:
    outcome = "NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_OPEN"

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "continuum_target_loaded": False,
    "polarization_target_loaded": False,
    "speed_or_mass_target_loaded": False,
    "global_carrier_controls": global_controls,
    "all_sector_dimension_census": {
        str(sector): {
            "irrep_dimension": DIMENSIONS[sector],
            "longitudinal": census_pairs[sector][0],
            "transverse": census_pairs[sector][1],
        }
        for sector in range(7)
    } if census_uniform else {},
    "weighted_longitudinal_dimension": weighted_longitudinal,
    "weighted_transverse_dimension": weighted_transverse,
    "parities": records,
    "projector_label_counts": dict(projector_counts),
    "cross_label_counts": dict(cross_counts),
    "sign_label_counts": dict(sign_counts),
    "dynamic_label_counts": dict(dynamic_counts),
    "rotated_control_label_counts": dict(control_counts),
    "classification": {
        "action_weighted_longitudinal_identity": (
            "DERIVED COMPUTATIONAL / STRUCTURAL" if positive_identity
            else "REFUTED" if resolved_refutation else "OPEN"
        ),
        "exact_gauge_interpretation": "OPEN",
        "physical_tensor_quotient": "OPEN",
        "continuum_refinement": "OPEN",
        "dispersion_speed_mass": "NOT COMPUTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
print(f"all-sector dimensions: {census_pairs}")
print(f"projector labels: {dict(projector_counts)}")
print(f"cross labels: {dict(cross_counts)}")
print(f"sign labels: {dict(sign_counts)}")
print(f"dynamic labels: {dict(dynamic_counts)}")
print(f"control labels: {dict(control_counts)}")
if passed != tests:
    raise SystemExit(1)

